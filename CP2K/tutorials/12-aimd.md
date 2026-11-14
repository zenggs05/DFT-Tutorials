# 第12章 从头算分子动力学（AIMD）

本章讲解如何用 CP2K Quickstep 做 Born–Oppenheimer 从头算分子动力学（ab initio molecular dynamics, AIMD）：`RUN_TYPE MD` 与 `&MOTION/&MD` 的写法，NVE/NVT/NPT_I/NPT_F/LANGEVIN 系综的选择，CSVR 与 Nose–Hoover 链恒温器（`REGION GLOBAL/MOLECULE/MASSIVE`），NPT 必需的 `STRESS_TENSOR ANALYTICAL` 与提高 CUTOFF，ASPC 波函数外推（`EXTRAPOLATION ASPC`、`EXTRAPOLATION_ORDER 3`）与 OT + `OUTER_SCF`，轨迹、速度、力、晶胞与重启文件的输出控制，`-1.ener` 各列（Kin、Temp、Pot、Cons Qty）的含义，`&EXT_RESTART` 续算。之后讲 AIMD 的质量检查（守恒量漂移、温度分布）和后处理（RDF、MSD/扩散系数、VACF/振动态密度、偶极自相关红外谱），最后介绍用 REFTRAJ 重算轨迹，包括 CP2K 2026.2 起的 `EVAL` 写法和从 EXTXYZ 读取晶胞。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出完整的 CP2K AIMD 输入，说清 `ENSEMBLE`、`STEPS`、`TIMESTEP`、`TEMPERATURE` 以及恒温器、压浴各关键字的作用；
- 按体系选择步长（含 H 体系取 0.5 fs）、恒温器类型和时间常数，知道什么时候必须用 NVE、什么时候不能用 `REGION MASSIVE`；
- 正确设置 NPT：`STRESS_TENSOR ANALYTICAL`、更高的 CUTOFF、`CELL_REF`，并会读晶胞与压力输出；
- 用 ASPC 外推和 OT + `OUTER_SCF` 让每个 MD 步只需少量 SCF 迭代；
- 读懂 `-1.ener`、`-pos-1.xyz`、`-vel-1.xyz`、`-frc-1.xyz`、`-1.cell`、`-1.restart`、`-RESTART.wfn` 等输出；
- 用 `&EXT_RESTART` 在作业超时或换系综后无缝续算；
- 判断一条 AIMD 轨迹是否可信（能量守恒、温度分布、SCF 行为），并从轨迹算出 RDF、扩散系数、振动态密度和红外谱；
- 用 REFTRAJ 对已有轨迹换方法重算能量和力（例如为机器学习势打标签），并处理 2026.2 的语法变化。

| 任务 | 推荐做法 | 本章位置 |
|---|---|---|
| 液体、溶液的结构与热力学 | NVT（CSVR）平衡 → NVT 生产 | §4、§6 |
| 检验步长与 SCF 精度 | 短 NVE，看守恒量漂移 | §7 |
| 扩散、振动谱、红外谱等动力学性质 | NVT 平衡后切换到 NVE 采样 | §7、§11 |
| 密度、热膨胀、相变 | NPT_I（液体）或 NPT_F（晶体、界面） | §9 |
| 快速、稳健的升温或随机热浴 | LANGEVIN 系综 | §8 |
| 换泛函/基组重算轨迹、为 MLP 准备标签 | REFTRAJ + `EVAL ENERGY_FORCES` | §12 |
| 稀有事件（反应能垒、相变自由能） | 普通 AIMD 采样不到，用元动力学或约束 MD | [第15章](15-metadynamics-plumed.md) |
| 0 K 下的反应路径与过渡态 | NEB/Dimer | [第13章](13-neb-dimer.md) |

本章所有输入都在 [../examples/12-aimd-water/README.md](../examples/12-aimd-water/README.md)，体系是 32 个水分子的立方周期盒子（96 个原子），工作流如下：

```text
gen_water_box.py -> water32.xyz
        |
        v
water32_nvt_csvr.inp (NVT equilibration, CSVR 50 fs, 2 ps)
        |  water32_nvt-1.restart + water32_nvt-RESTART.wfn
        +--> md_restart.inp        (NVT production, CSVR 200 fs, 5 ps) --> reftraj_eval.inp (TZV2P re-evaluation)
        +--> water32_nvt_nose.inp  (alternative production: Nose-Hoover chain)
        +--> water32_nve.inp       (NVE, energy conservation, velocities + dipole for VDOS/IR)
        +--> water32_npt_i.inp     (NPT_I, 1 atm, density)                 --> reftraj_npt.inp (variable cell)
```

> [!NOTE]
> 32 个水分子是教学用的最小合理体系。AIMD 每一步都要做一次完整的 SCF，成本远高于经典 MD：这个例子 2 ps 的平衡段在几十个核上也要数小时（取决于硬件）。先把 `STEPS` 改成几十步试跑，读出 `MD| CPU time per MD step` 后再估算总机时。正式研究中盒子大小、模拟时长都要按待测性质的关联长度和关联时间来定。

---

## 2. 原理速览

### 2.1 Born–Oppenheimer MD 与速度 Verlet 积分

BOMD 假设电子始终处于当前核构型下的基态。每个时间步先对固定的核坐标 $\mathbf R(t)$ 做 SCF，得到能量 $E[\mathbf R(t)]$ 和 Hellmann–Feynman 力（GPW 下还包括基组相关的 Pulay 项）$\mathbf F_I = -\partial E/\partial \mathbf R_I$，再用速度 Verlet 推进核运动：

$$
\mathbf R_I(t+\Delta t) = \mathbf R_I(t) + \mathbf v_I(t)\,\Delta t + \frac{\mathbf F_I(t)}{2M_I}\Delta t^2
$$

$$
\mathbf v_I(t+\Delta t) = \mathbf v_I(t) + \frac{\mathbf F_I(t) + \mathbf F_I(t+\Delta t)}{2M_I}\Delta t
$$

速度 Verlet 是时间可逆的辛积分器。只要力是精确的、步长足够小，总能量 $E_\mathrm{tot} = E_\mathrm{pot} + \sum_I \tfrac12 M_I v_I^2$ 只会在一个常数附近小幅振荡，不会系统漂移。实际计算中 SCF 不可能完全收敛，力带有噪声，这正是 AIMD 能量漂移的主要来源（§2.5、§10）。

### 2.2 步长怎么定

步长 $\Delta t$ 要能分辨体系中最快的核运动，经验上取最快振动周期的 1/10–1/20。水中 O–H 伸缩约 3700 cm⁻¹，对应周期

$$
T = \frac{1}{\tilde\nu c} = \frac{1}{3700\ \mathrm{cm^{-1}}\times 3.0\times10^{10}\ \mathrm{cm\,s^{-1}}} \approx 9\ \mathrm{fs}
$$

所以**含 H 的体系取 0.5 fs**（CP2K 的默认值也是 0.5 fs）。只含重元素的固体可用 1–2 fs。若只关心构型采样、不关心 H 的真实动力学，可以把 H 的质量改成氘（`&KIND H` 中写 `MASS 2.014`），振动频率约降为原来的 $1/\sqrt2$，步长可放宽到约 1 fs。高温、反应性体系或初始结构较差时要用更小的步长。

### 2.3 恒温器：CSVR、Nose–Hoover 链与 Langevin

**CSVR**（canonical sampling through velocity rescaling，Bussi–Donadio–Parrinello 2007）每一步按如下随机微分方程重新标定总动能 $K$：

$$
dK = (\bar K - K)\frac{dt}{\tau} + 2\sqrt{\frac{K\bar K}{N_f\,\tau}}\,dW
$$

其中 $\bar K = \tfrac12 N_f k_\mathrm B T_0$ 是目标动能，$N_f$ 是自由度数，$\tau$ 就是 `&CSVR/TIMECON`，$dW$ 是 Wiener 噪声。第一项让动能向目标值弛豫，第二项给出正确的正则涨落。CSVR 严格产生正则分布，$\tau$ 小时耦合强、升降温快，$\tau$ 大时对动力学扰动小。官方 MD 文档把它列为日常 NVT 的首选。

**Nose–Hoover 链**在哈密顿量中引入一串热浴自由度 $\{\eta_k, p_{\eta_k}\}$，第一个热浴与粒子动能耦合，后面的热浴依次耦合前一个，以克服单个 Nose–Hoover 热浴在小体系或刚性体系中不遍历的问题。它是确定性的，对动力学的扰动较平滑，但离平衡较远时会出现大幅温度振荡。CP2K 中 `LENGTH` 是链长（默认 3），`TIMECON` 是热浴特征时间（默认 1000 fs）。

**Langevin 动力学**给每个原子加摩擦力和随机力：

$$
M_I\ddot{\mathbf R}_I = \mathbf F_I - \gamma M_I\dot{\mathbf R}_I + \boldsymbol\xi_I(t),\qquad
\langle \xi_{I\alpha}(t)\,\xi_{J\beta}(t')\rangle = 2\gamma M_I k_\mathrm B T_0\,\delta_{IJ}\delta_{\alpha\beta}\delta(t-t')
$$

`&LANGEVIN/GAMMA`（单位 fs⁻¹）就是 $\gamma$。Langevin 对每个自由度都“局部”控温，稳健但会改变扩散等动力学性质。

`&THERMOSTAT/REGION` 决定热浴接在哪里：`GLOBAL` 整个体系一个热浴（默认，对动力学扰动最小）；`MOLECULE` 每个分子种类一个；`MASSIVE` 每个自由度一个（平衡最快，但会严重破坏动力学性质）；`DEFINED` 按 `&DEFINE_REGION` 自定义分区。

### 2.4 压浴：NPT_I 与 NPT_F

NPT 系综在恒温器之外再加一个压浴（CP2K 实现的是 Martyna–Tobias–Klein 类型的扩展系统方法），把晶胞也当作动力学变量。压浴的驱动力是瞬时压力（由应力张量给出）与目标压力之差。`NPT_I` 只允许体积各向同性变化，适合液体；`NPT_F` 允许晶胞形状和体积都变化，适合晶体、界面和各向异性材料。`&BAROSTAT/TIMECON` 控制晶胞运动的特征时间。

AIMD 中的应力来自 DFT，平面波（GPW 的辅助网格）在晶胞变化时会带来 **Pulay 应力**：截断能不够高时，压力系统偏低、晶胞被“压缩”。因此 NPT 必须：

1. 在 `&FORCE_EVAL` 设 `STRESS_TENSOR ANALYTICAL`，否则 CP2K 根本不计算维里，NPT 直接中止（§常见问题 3）；
2. 把 CUTOFF 提高到明显高于 NVT 的值（例如 400 → 600 Ry），截断能的收敛测试方法见 [第4章](04-cutoff-convergence.md)；
3. 用 `&CELL_REF` 指定一个略大的参考晶胞，让 FFT 网格点数按参考晶胞确定并在模拟中保持不变，避免盒子膨胀时有效截断能下降。

### 2.5 波函数外推：ASPC

相邻 MD 步的结构非常接近，上一步的波函数或密度矩阵就是很好的初猜。`&QS/EXTRAPOLATION` 选择外推方式，默认值就是 `ASPC`（always stable predictor–corrector，Kolafa 2004；Quickstep 中的实现见 VandeVondele 2005）。ASPC 用前 $K+2$ 步的结果按固定系数外推出预测值，再与本步一次 SCF 的结果按系数 $\omega = (K+2)/(2K+3)$ 混合；$K$ 就是 `EXTRAPOLATION_ORDER`（默认 3，官方建议 2–4）。它牺牲一点初猜精度，换来长时间积分的稳定性。效果好时，每个 MD 步只需 3–8 次 SCF 迭代（经验范围，以你的输出为准）。

外推只在连续 MD 中有意义。帧间隔很大的 REFTRAJ 重算（§12）应改用 `USE_PREV_WF`，或干脆用 `SCF_GUESS ATOMIC`。

### 2.6 守恒量

`-1.ener` 第 6 列 `Cons Qty` 是当前系综下应当守恒的量：

- NVE：$E_\mathrm{pot} + K$；
- Nose–Hoover NVT/NPT：再加上热浴（和压浴）自由度的动能与势能；
- CSVR：再加上恒温器累计从体系取走的能量（Bussi 2007 定义的“有效能量”）；
- Langevin：摩擦和随机力使守恒量不再严格守恒，主要看温度是否正常。

无论哪个系综，`Cons Qty` 的**系统漂移**都是检验步长和 SCF 精度最灵敏的指标。

---

## 3. 示例体系与准备

### 3.1 初始结构

[../examples/12-aimd-water/README.md](../examples/12-aimd-water/README.md) 中的 `gen_water_box.py` 只用 Python 标准库，把 32 个水分子的 O 原子放在 fcc 格点上、分子随机取向，盒长按 0.997 g/cm³ 算得 $L = 9.8654$ Å：

```bash
cd CP2K/examples/12-aimd-water
python gen_water_box.py --seed 2026 -o water32.xyz
head -3 water32.xyz
```

这只是“没有原子重叠、密度正确”的起点，**没有氢键网络**，前 1–2 ps 的数据只能作平衡用。正式研究更好的做法是先用经典力场（如 GROMACS、LAMMPS）或 Packmol 生成并平衡初始结构，再交给 AIMD。对于晶体或表面，先做几何优化（[第7章](07-geo-cell-opt.md)）能显著降低初始受力，避免开局温度飙升。

### 3.2 电子结构设置

全部示例使用同一套设置：PBE-D3(BJ)；`DZVP-MOLOPT-PBE-GTH-q1/q6` 基组配 `GTH-PBE-q1/q6` 赝势（`BASIS_MOLOPT_UZH`/`POTENTIAL_UZH`，见 [第3章](03-gpw-basis-pseudo.md)）；`CUTOFF 400 / REL_CUTOFF 60` Ry（NPT 用 600 Ry，收敛测试见 [第4章](04-cutoff-convergence.md)）；OT（DIIS + `FULL_SINGLE_INVERSE`）+ `OUTER_SCF`，`EPS_SCF 1E-6`（SCF 策略见 [第5章](05-scf-strategy.md)）；ASPC 外推，阶数 3；步长 0.5 fs。纯 PBE 水结构偏强、扩散偏慢是已知问题（§11.5），本章只把它当作教学基线。

> [!IMPORTANT]
> AIMD 的力精度由 `EPS_SCF` 和 `EPS_DEFAULT` 共同决定。`EPS_SCF 1E-6` 配合 ASPC 是官方和大量文献使用的起点，但 `EPS_SCF` 设得越松，力噪声越大，能量漂移越明显。正式生产前一定要做 §7 的 NVE 测试。

---

## 4. 实战一：NVT 平衡（CSVR 恒温器）

第一步从 `water32.xyz` 出发，用强耦合的 CSVR（`TIMECON 50 fs`）把体系快速带到 300 K 并初步建立氢键网络。完整输入 `water32_nvt_csvr.inp`：

```fortran
# Step 1: NVT equilibration of 32 H2O with PBE-D3(BJ), CSVR thermostat
#   starting structure: water32.xyz (fcc oxygens + random orientations, NOT equilibrated)
#   2 ps with a short thermostat time constant (50 fs) to reach 300 K quickly.
# Continue with md_restart.inp (production NVT, TIMECON 200 fs).

&GLOBAL
  PROJECT     water32_nvt
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE    NVT
    STEPS       4000              ! 4000 x 0.5 fs = 2 ps
    TIMESTEP    0.5               ! fs, light H atoms
    TEMPERATURE 300.0             ! K, target and initial velocities
    &THERMOSTAT
      TYPE   CSVR
      REGION GLOBAL
      &CSVR
        TIMECON 50.0              ! fs, strong coupling for equilibration
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &VELOCITIES ON
      &EACH
        MD 10
      &END EACH
    &END VELOCITIES
    &FORCES ON
      &EACH
        MD 10
      &END EACH
    &END FORCES
    &RESTART
      BACKUP_COPIES 2
      &EACH
        MD 20
      &END EACH
    &END RESTART
    &RESTART_HISTORY
      &EACH
        MD 1000
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     400
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT         1.0E-12
      EXTRAPOLATION       ASPC    ! default in recent versions, written explicitly
      EXTRAPOLATION_ORDER 3
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   30
      &OT
        MINIMIZER      DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
      &PRINT
        &RESTART
          BACKUP_COPIES 0
          &EACH
            MD 10
          &END EACH
        &END RESTART
        &RESTART_HISTORY OFF
        &END RESTART_HISTORY
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 9.8654 9.8654 9.8654
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   water32.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### 4.1 `&MD` 与恒温器关键字

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `RUN_TYPE` | `MD` | 在 `&GLOBAL` 中，开启分子动力学 |
| `ENSEMBLE` | `NVT` | 默认 `NVE`。可选 `NVE`、`NVT`、`NPT_I`、`NPT_F`、`LANGEVIN`、`REFTRAJ`，以及较专门的 `NPE_I`、`NPE_F`、`ISOKIN`、`MSST`、`MSST_DAMPED`、`HYDROSTATICSHOCK`、`NVT_ADIABATIC`、`NPT_IA` |
| `STEPS` | 视需要 | 本次运行要走的步数，从 `STEP_START_VAL`（续算时即重启文件中的步号）开始计。默认只有 3 步，**一定要写** |
| `MAX_STEPS` | 可选 | 从第 1 步算起的总步数上限（默认 10⁹），适合续算时给整条轨迹设硬上限；到达上限那一步的输出默认不写，见 §6 |
| `TIMESTEP` | 0.5 fs（含 H） | 单位 fs，默认 0.5 |
| `TEMPERATURE` | 目标温度 | 单位 K，默认 300。既是 NVT/NPT 的目标温度，也用来生成初始随机速度（续算读入速度时不再用于初始化） |
| `&THERMOSTAT/TYPE` | `CSVR` | 默认是 `NOSE`。可选 `NOSE`、`CSVR`、`GLE`、`AD_LANGEVIN` |
| `&THERMOSTAT/REGION` | `GLOBAL` | 默认 `GLOBAL`。另有 `MOLECULE`、`MASSIVE`、`DEFINED`、`THERMAL`、`NONE`，见 §8 |
| `&CSVR/TIMECON` | 平衡 50 fs，通用 100 fs，生产 ≥200 fs | 单位 fs，**默认 1000 fs**。越小耦合越强 |
| `TEMP_TOL` | 不用 | 温度偏离超过该值时粗暴地重标速度。官方文档已注明该关键字过时，建议用短 `TIMECON` 的 CSVR 代替 |
| `COMVEL_TOL`、`ANGVEL_ZERO` | 孤立体系才用 | 去除质心平动和整体转动；周期体系一般不需要 |

### 4.2 电子结构部分的关键字

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `&QS/EXTRAPOLATION` | `ASPC` | 默认值。其他选项有 `PS`、`USE_PREV_WF`、`USE_PREV_P`、`LINEAR_P`、`GEXT_PROJ`、`GEXT_PROJ_QTR` 等。官方建议常规 BOMD 先用 ASPC |
| `&QS/EXTRAPOLATION_ORDER` | 3 | ASPC/PS 一般取 2–4；GEXT 系列可取 4–10 |
| `&QS/EPS_DEFAULT` | 1.0E-12 | 积分筛选阈值；太松会引入力噪声 |
| `&SCF/EPS_SCF` | 1.0E-6 | AIMD 常用值；NVE 漂移大时收紧到 1E-7 |
| `&SCF/MAX_SCF` | 20–30 | OT 内层最大迭代次数 |
| `&OUTER_SCF` | `MAX_SCF` 8–10，`EPS_SCF` 与内层相同 | 内层没收敛时重新构造预条件子再来一轮，显著提高 OT 的稳健性 |
| `&OT/PRECONDITIONER` | `FULL_SINGLE_INVERSE` 或 `FULL_ALL` | 绝缘体常用；金属或小带隙体系改用对角化 + 展宽（[第5章](05-scf-strategy.md)） |
| `&SCF/&PRINT/&RESTART/&EACH/MD` | 10 | 每 10 个 MD 步写一次 `-RESTART.wfn`，减少 I/O；`BACKUP_COPIES 0` 不保留旧 wfn |
| `SCF_GUESS` | 首次 `ATOMIC`，续算 `RESTART` | 续算时配合 `WFN_RESTART_FILE_NAME` |

> [!WARNING]
> 不要在 MD 中设 `IGNORE_CONVERGENCE_FAILURE`。自 2024.1 起，SCF 不收敛默认会中止计算，这是对 AIMD 的保护：未收敛步的力是错的，会把错误悄悄带进后面所有步。遇到偶发不收敛，应当按 §常见问题 1 调整 SCF 设置，再从最近的重启文件续算。

### 4.3 运行

```bash
cd CP2K/examples/12-aimd-water
export CP2K_DATA_DIR=/path/to/cp2k/data
cp2k.psmp -c water32_nvt_csvr.inp                       # syntax check only
mpirun -np 32 cp2k.psmp -i water32_nvt_csvr.inp -o water32_nvt_csvr.out
# or on a Slurm cluster:
sbatch job_slurm.sh water32_nvt_csvr.inp
```

示例目录里的 `job_slurm.sh` 每个作业只跑一个输入，配合 `sbatch --dependency=afterok:<jobid>` 把平衡、生产、NVE、NPT 串起来，命令见示例 README。通用的 Slurm/PBS 脚本是本仓库的 `../scripts/job_cp2k_slurm.sh` 和 `../scripts/job_cp2k_pbs.sh`，调度系统的基础见 [../../Common/tutorials/02-job-scheduler.md](../../Common/tutorials/02-job-scheduler.md)。

---

## 5. 读懂 MD 输出

### 5.1 主输出中的 `MD|` 块

每个 MD 步结束时，主输出会打印一个 `MD|` 块（`PRINT_LEVEL LOW` 下也会打印）：

```text
 MD| ***************************************************************************
 MD| Step number                                                            1200
 MD| Time [fs]                                                        600.000000
 MD| Conserved quantity [hartree]                           -0.550812345678E+03
 MD| ---------------------------------------------------------------------------
 MD|                                          Instantaneous             Averages
 MD| CPU time per MD step [s]                         4.112345             4.098765
 MD| Energy drift per atom [K]          0.123456789012E+01   0.876543210987E+00
 MD| Potential energy [hartree]        -0.550912345678E+03  -0.550910123456E+03
 MD| Kinetic energy [hartree]           0.135678901234E+00   0.136789012345E+00
 MD| Temperature [K]                             298.765432           301.234567
 MD| ***************************************************************************
```

（示意输出，数值仅供参考）

- `Conserved quantity`：即 §2.6 的守恒量；
- `Energy drift per atom [K]`：守恒量相对第一步的偏离，除以原子数并换算成开尔文（1 meV ≈ 11.6 K）。“Averages”列是运行平均。它是最快捷的漂移指标；
- `CPU time per MD step [s]`：用来估算总机时；
- NPT 时还会出现 `Pressure [bar]`、`Barostat temperature [K]`、`Cell volume [ang^3]`、`Cell lengths [ang]` 等行。

```bash
grep "MD| Temperature"                water32_nvt_csvr.out | tail -n 5
grep "MD| Energy drift per atom"      water32_nvt_csvr.out | tail -n 5
grep "MD| CPU time per MD step"       water32_nvt_csvr.out | tail -n 1
grep -c "SCF run converged"           water32_nvt_csvr.out    # one per MD step
# histogram of SCF iterations per MD step ("*** SCF run converged in N steps ***")
grep "SCF run converged" water32_nvt_csvr.out | awk '{print $6}' | sort -n | uniq -c
```

如果最后一条命令在 `PRINT_LEVEL LOW` 下统计不到结果，可以在 `&SCF/&PRINT` 中加 `&PROGRAM_RUN_INFO ON`，或者改用 `PRINT_LEVEL MEDIUM` 做短测试。

### 5.2 输出文件一览

文件名中的 `water32_nvt` 来自 `PROJECT`，`-1` 是迭代层级编号：

| 文件 | 由谁控制 | 内容与单位 |
|---|---|---|
| `water32_nvt-1.ener` | `&MD/&PRINT/&ENERGY`（默认打开，每步） | 每步一行，各列见 §5.3 |
| `water32_nvt-pos-1.xyz` | `&MOTION/&PRINT/&TRAJECTORY` | 坐标（Å），**不折回盒子**（便于算 MSD）；标题行含 `i = 步号, time = fs, E = 势能(Ha)` |
| `water32_nvt-vel-1.xyz` | `&VELOCITIES`（默认级别 HIGH，需写 `ON`） | 速度，默认单位 bohr/au_time（1 au_time ≈ 0.0242 fs） |
| `water32_nvt-frc-1.xyz` | `&FORCES`（默认级别 HIGH，需写 `ON`） | 力，Hartree/bohr |
| `water32_npt-1.cell` | `&CELL`（默认 HIGH，需写 `ON`） | 步号、时间、三个晶胞矢量（Å）和体积（Å³） |
| `water32_npt-1.stress` | `&STRESS`（默认 HIGH，需写 `ON`） | 应力张量（bar） |
| `water32_nvt-1.restart` | `&RESTART`（默认每 20 步） | 运动重启文件，ASCII，本身就是完整的 CP2K 输入 |
| `water32_nvt-1.restart.bak-1` … | `BACKUP_COPIES` | 旧版本的重启文件 |
| `water32_nvt-1_1000.restart` | `&RESTART_HISTORY`（默认每 500 步） | 按步号保留的历史重启点，用于回退 |
| `water32_nvt-RESTART.wfn` | `&SCF/&PRINT/&RESTART` | 波函数，二进制 |

`&MOTION/&PRINT` 中各打印键都可以用 `&EACH/MD n` 设输出间隔。`&TRAJECTORY` 还有 `FORMAT` 选项：默认 `XMOL`（即普通 xyz），另有 `DCD`、`PDB`、`ATOMIC`，以及 2026.1 起的 `EXTXYZ`（标题行含 `Lattice="..."` 晶胞信息，ASE 和 OVITO 可直接读）。

### 5.3 `-1.ener` 各列含义

```text
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
```

| 列 | 名称 | 含义 |
|---|---|---|
| 1 | `Step Nr.` | MD 步号（续算时接着上次的步号，除非 `RESTART_COUNTERS F`） |
| 2 | `Time[fs]` | 模拟时间 |
| 3 | `Kin.[a.u.]` | 核动能（Hartree） |
| 4 | `Temp[K]` | 瞬时温度 $T = 2K/(N_f k_\mathrm B)$ |
| 5 | `Pot.[a.u.]` | 势能，即本步 DFT 总能量（Hartree） |
| 6 | `Cons Qty[a.u.]` | 守恒量（Hartree），定义随系综变化，见 §2.6 |
| 7 | `UsedTime[s]` | 本步墙钟时间 |

快速统计平均温度和涨落（跳过前 1000 步）：

```bash
awk '!/^#/ && $1>=1000 {n++; s+=$4; ss+=$4*$4} END {m=s/n; printf "T = %.1f +- %.1f K (%d steps)\n", m, sqrt(ss/n-m*m), n}' water32_nvt-1.ener
```

画温度、势能和守恒量随时间的曲线，可以用本仓库的 `../scripts/md_analyze.py`（`.ener` 统计、漂移、RDF、MSD、VACF、红外谱，用 `python ../scripts/md_analyze.py -h` 查看用法），也可以直接用 gnuplot：

```bash
gnuplot -p -e "set xlabel 'time / fs'; plot 'water32_nvt-1.ener' u 2:4 w l t 'T / K'"
```

> [!TIP]
> 平衡阶段，温度曲线应在 1–2 ps 内从过冲回落并围绕 300 K 波动，势能应停止单调下降、转为围绕平均值涨落。只有这两个量都“平了”，才算进入可以采样的阶段。

---

## 6. 实战二：续算与生产段（`&EXT_RESTART`）

作业超时、想延长轨迹、想换恒温器参数，都通过 `&EXT_RESTART` 续算。它从上一次的 `.restart` 文件中读出坐标、速度、恒温器状态和步号计数器，再配合 `.wfn` 作为 SCF 初猜。`md_restart.inp` 从平衡段的终点出发，换成弱耦合（`TIMECON 200 fs`）做 5 ps 生产：

```fortran
&GLOBAL
  PROJECT     water32_prod
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&EXT_RESTART
  RESTART_FILE_NAME water32_nvt-1.restart
  ! defaults restart everything that is meaningful: positions, velocities,
  ! thermostat, counters. Written out here for clarity.
  RESTART_POS        T
  RESTART_VEL        T
  RESTART_THERMOSTAT T
  RESTART_COUNTERS   T
&END EXT_RESTART

&MOTION
  &MD
    ENSEMBLE    NVT
    STEPS       10000             ! 10000 MORE steps (5 ps) after the restart point
    ! MAX_STEPS (cap on the TOTAL step counter) is deliberately not set: the step
    ! that reaches MAX_STEPS is treated as a stop, and without FORCE_LAST its
    ! output and restart file are not written (see chapter 12, section 6).
    TIMESTEP    0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE   CSVR
      REGION GLOBAL
      &CSVR
        TIMECON 200.0             ! fs, gentle coupling for production
      &END CSVR
    &END THERMOSTAT
  &END MD
  ! &PRINT: same as step 1 (TRAJECTORY/VELOCITIES/FORCES every 10 steps,
  !         RESTART every 20 steps, RESTART_HISTORY every 2000 steps)
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    WFN_RESTART_FILE_NAME water32_nvt-RESTART.wfn
    ! ... MGRID, QS, XC identical to step 1 ...
    &SCF
      SCF_GUESS RESTART
      ! ... OT, OUTER_SCF, PRINT identical to step 1 ...
    &END SCF
  &END DFT
  ! &SUBSYS identical to step 1; positions are taken from the restart file
&END FORCE_EVAL
```

上面省略的部分与 §4 完全相同，完整文件见示例目录中的 `md_restart.inp`。

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `RESTART_FILE_NAME` | `<project>-1.restart` | 也可以填 `RESTART_HISTORY` 留下的 `<project>-1_<step>.restart` 回退到某个历史点 |
| `RESTART_DEFAULT` | T（默认） | 所有可重启的量默认都重启；设为 F 后只重启显式写成 T 的量 |
| `RESTART_POS` / `RESTART_VEL` | T | 坐标、速度 |
| `RESTART_THERMOSTAT` | 同类恒温器续算时为 T；换恒温器类型时为 F | Nose 链/CSVR 的内部状态 |
| `RESTART_BAROSTAT` | NPT 续算时为 T（默认） | 压浴的质量和速度；重启文件中没有 `&BAROSTAT`（例如从 NVT 起跑 NPT）时自动跳过 |
| `RESTART_CELL` | NPT 续 NPT 时为 T（默认）；**从定容运行的重启文件起跑 NPT 并新加 `CELL_REF` 时必须设为 F** | 用重启文件中的 `&CELL` **整段替换**输入中的 `&CELL`（含 `CELL_REF` 子段）。NVT 的重启文件里没有 `CELL_REF`，不设 F 时你在输入中新写的 `CELL_REF` 会被静默丢掉 |
| `RESTART_COUNTERS` | 续同一条轨迹为 T；开始新阶段为 F | 步号和时间；为 F 时从 0 开始计 |
| `RESTART_RANDOMG` | 默认 T | 随机数发生器状态（CSVR、Langevin 用到），保证续算的随机序列连续 |
| `WFN_RESTART_FILE_NAME` | `<project>-RESTART.wfn` | 在 `&DFT` 中，配合 `SCF_GUESS RESTART` |

几条要点：

1. **新阶段用新的 `PROJECT`**（这里是 `water32_prod`），各阶段的输出互不覆盖，平衡段丢弃起来也方便。如果只是因为超时要把同一条轨迹接着跑，可以沿用原 `PROJECT`，此时 CP2K 会接着往 `-pos-1.xyz`、`-1.ener` 等文件末尾追加。
2. 续算点在最后一个写重启文件的步上（默认每 20 步）。作业被杀时，`.ener` 和轨迹中可能已经多写了几步，接着跑就会出现**重复帧**。分析前要按步号去重；`cp2k-output-tools` 提供的命令行工具 `xyz_restart_cleaner`（用法 `xyz_restart_cleaner old.xyz new.xyz`）就是做这件事的。
3. `STEPS` 计的是**本次还要走的步数**，`MAX_STEPS` 计的是**总步号上限**。上例从 4000 步续算，`STEPS 10000` 走到第 14000 步。`MAX_STEPS` 适合“同一个输入反复提交、直到总步数够了为止”的用法（`STEPS` 设大，`MAX_STEPS` 设为目标总步数）。但要注意：按 2026.2 的 `md_run.F`，计数器达到 `MAX_STEPS` 的那一步被当作“停止”处理，和 `WALLTIME`/`EXIT` 一样，默认不写该步的 `.ener`、轨迹和重启文件，除非加 `&MD/&PRINT/FORCE_LAST T`；起始步号已经 ≥ `MAX_STEPS` 时则直接报错 `maximum step number smaller than initial step value`。所以示例中只用 `STEPS`，没有写 `MAX_STEPS`。
4. `.restart` 文件本身就是完整输入，也可以直接 `cp2k.psmp -i water32_nvt-1.restart` 运行，但它会继承原来的 `PROJECT` 和所有设置，不便修改。更推荐像上面这样写一个新输入并通过 `&EXT_RESTART` 引用，修改重启文件的辅助脚本见 `../scripts/restart_update.py`（用 `python ../scripts/restart_update.py -h` 查看用法）。
5. 换赝势或基组的 q 值、改变原子数或原子顺序后，`.wfn` 不再适用；这时把 `SCF_GUESS` 改回 `ATOMIC`。

> [!TIP]
> 防止作业在写重启文件前被调度器杀掉：在 `&GLOBAL` 中设 `WALLTIME`（秒数或 `HH:MM:SS`），比队列时限少几十分钟，CP2K 到时会在当前 MD 步结束后正常退出（输出末尾有 `PROGRAM ENDED AT`）。也可以在运行目录中新建一个名为 `EXIT`（或 `<PROJECT>.EXIT`）的空文件，CP2K 检测到后会在当前步结束时停止。注意：这两种方式默认**不会**为最后一步补写输出和重启文件，重启点仍是最近一次按 `&RESTART/&EACH` 写出的步；如果希望停止时补写最后一步，在 `&MOTION/&MD/&PRINT` 中加 `FORCE_LAST T`（`&RESTART` 打印键默认已带 `ADD_LAST NUMERIC`）。停下后检查 `.restart` 中的 `STEP_START_VAL`，再按本节方法续算。

---

## 7. 实战三：NVE 与能量守恒检查

NVE 没有恒温器，是检验步长和 SCF 精度最直接的办法，也是计算扩散、振动谱、红外谱等动力学性质时首选的采样系综。`water32_nve.inp` 从平衡态出发跑 2 ps，每 2 步（1 fs）输出一次速度和体系偶极矩：

```fortran
&GLOBAL
  PROJECT     water32_nve
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&EXT_RESTART
  RESTART_FILE_NAME  water32_nvt-1.restart
  RESTART_POS        T
  RESTART_VEL        T
  RESTART_THERMOSTAT F
  RESTART_COUNTERS   F
&END EXT_RESTART

&MOTION
  &MD
    ENSEMBLE    NVE
    STEPS       4000              ! 2 ps
    TIMESTEP    0.5
    TEMPERATURE 300.0             ! not used when velocities are restarted
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 2
      &END EACH
    &END TRAJECTORY
    &VELOCITIES ON
      &EACH
        MD 2
      &END EACH
    &END VELOCITIES
    &RESTART
      BACKUP_COPIES 2
      &EACH
        MD 20
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    ! ... BASIS/POTENTIAL files, MGRID, QS, SCF (SCF_GUESS RESTART), XC as in md_restart.inp ...
    WFN_RESTART_FILE_NAME water32_nvt-RESTART.wfn
    &PRINT
      ! total dipole of the periodic cell (Berry phase), appended every 2 MD steps
      &MOMENTS ON
        PERIODIC T
        FILENAME =water32_nve_dipole.dat
        &EACH
          MD 2
        &END EACH
      &END MOMENTS
    &END PRINT
  &END DFT
  ! ... &SUBSYS as in step 1 ...
&END FORCE_EVAL
```

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `ENSEMBLE NVE` | — | 无恒温器，`&THERMOSTAT` 段不起作用 |
| `RESTART_THERMOSTAT F` | — | NVE 不需要热浴状态 |
| `&VELOCITIES` `&EACH MD 2` | 1–2 步 | VACF/VDOS 需要密集的速度输出：间隔 $\Delta\tau$ 决定可分辨的最高频率 $\tilde\nu_\mathrm{max} = 1/(2c\Delta\tau)$，1 fs 对应约 16700 cm⁻¹，足够 |
| `&MOMENTS` `PERIODIC T` | — | 周期体系用 Berry 相公式计算总偶极（默认级别 HIGH，需写 `ON`）；`FILENAME =name` 中的 `=` 表示文件名原样使用 |
| `&MOMENTS/&EACH/MD` | 与速度相同 | 偶极自相关求红外谱用 |

检查能量守恒：

```bash
# drift of the conserved quantity per atom and per ps (96 atoms)
awk '!/^#/ {if (!f) {c0=$6; t0=$2; f=1}; c=$6; t=$2} END {printf "drift = %.2e Ha/atom/ps\n", (c-c0)/96/((t-t0)/1000)}' water32_nve-1.ener
grep "MD| Energy drift per atom" water32_nve.out | tail -n 3
```

判据（经验值，取决于研究目标）：

- `Cons Qty` 只有快速的小幅振荡，没有单调趋势；
- 长时间漂移小于约 1 meV/atom/ps（≈ 3.7×10⁻⁵ Ha/atom/ps，≈ 11.6 K/atom/ps）；
- 温度围绕平均值波动且没有系统上升，NVE 的平均温度由初始状态决定，不一定正好是 300 K。

漂移偏大时依次检查：`EPS_SCF` 收紧到 1E-7、`EPS_DEFAULT` 收紧到 1E-14、CUTOFF 是否收敛、步长减到 0.25 fs 看漂移是否按 $\Delta t^2$ 下降。如果只有收紧 `EPS_SCF` 才有效，说明问题在力噪声而不在积分器。

---

## 8. 恒温器的选择：CSVR、Nose–Hoover 链与 Langevin

### 8.1 Nose–Hoover 链

示例目录中的 `water32_nvt_nose.inp` 是生产段的替代写法，只有 `&MOTION/&MD` 和 `&EXT_RESTART` 与 `md_restart.inp` 不同：

```fortran
&EXT_RESTART
  RESTART_FILE_NAME  water32_nvt-1.restart
  RESTART_POS        T
  RESTART_VEL        T
  RESTART_THERMOSTAT F            ! CSVR state cannot be reused by a Nose chain
  RESTART_COUNTERS   F
&END EXT_RESTART

&MOTION
  &MD
    ENSEMBLE    NVT
    STEPS       10000             ! 5 ps
    TIMESTEP    0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE   NOSE
      REGION GLOBAL
      &NOSE
        LENGTH  3                 ! chain length (default 3)
        YOSHIDA 3                 ! Suzuki-Yoshida integrator order (default 3)
        MTS     2                 ! multiple time steps for the chain (default 2)
        TIMECON 200.0             ! fs (default 1000 fs); smaller = stronger coupling
      &END NOSE
    &END THERMOSTAT
  &END MD
  ! &PRINT as in md_restart.inp
&END MOTION
```

| 关键字 | 默认值 | 说明 |
|---|---|---|
| `&NOSE/LENGTH` | 3 | 链长。3–4 足够，更长没有明显好处 |
| `&NOSE/TIMECON` | 1000 fs | 热浴特征时间。太短会和体系的振动模式共振，太长则控温很慢；常见范围 100–1000 fs |
| `&NOSE/YOSHIDA` | 3 | 热浴运动方程的 Suzuki–Yoshida 积分阶数 |
| `&NOSE/MTS` | 2 | 热浴积分的多时间步数 |

### 8.2 `REGION` 的影响

| `REGION` | 热浴个数 | 适用 | 注意 |
|---|---|---|---|
| `GLOBAL` | 1 | 默认；需要动力学性质时的唯一合理选择 | 远离平衡时温度在各部分间分配不均，平衡较慢 |
| `MOLECULE` | 每个分子种类一个 | 多组分体系（溶质 + 溶剂）各自控温 | 需要 CP2K 能识别分子（由 `&TOPOLOGY` 连接性生成） |
| `MASSIVE` | 每个自由度一个 | 只求结构和热力学平均、快速平衡、强振动耦合差的体系 | **破坏扩散系数、VACF、红外谱等动力学量** |
| `DEFINED` | 由 `&DEFINE_REGION` 给定 | 界面、表面吸附体系中给不同区域单独控温 | 每个区域原子数不能太少，否则温度涨落极大 |

检查能量均分：在 `&MD` 中写 `TEMP_KIND T`，并在 `&MD/&PRINT` 中打开 `&TEMP_KIND ON`，可以按元素分别输出温度。用 `GLOBAL` 恒温器时，若 H 和 O 的平均温度长期不一致，说明体系还没平衡，或者恒温器耦合太弱。

### 8.3 Langevin 系综

Langevin 不是 `&THERMOSTAT` 的一种，而是独立的系综：

```fortran
&MOTION
  &MD
    ENSEMBLE    LANGEVIN
    STEPS       4000
    TIMESTEP    0.5
    TEMPERATURE 300.0
    &LANGEVIN
      GAMMA 0.001                 ! fs^-1, i.e. friction time 1 ps
    &END LANGEVIN
  &END MD
&END MOTION
```

`GAMMA` 的单位是 fs⁻¹，默认 0（即无摩擦）。值越小越接近 NVE，越大控温越强。`NOISY_GAMMA` 是 Kühne 等人的第二代 Car–Parrinello MD 中用来补偿噪声力的附加摩擦项，普通 BOMD 不需要设置。

---

## 9. 实战四：NPT 模拟（NPT_I 与 NPT_F）

### 9.1 NPT_I：液态水的密度

`water32_npt_i.inp` 从 NVT 平衡态出发，在 1 atm、300 K 下跑 3 ps 的各向同性 NPT。与 NVT 输入相比，变化的部分如下：

```fortran
&EXT_RESTART
  RESTART_FILE_NAME  water32_nvt-1.restart
  RESTART_POS        T
  RESTART_VEL        T
  RESTART_THERMOSTAT F            ! new thermostat settings below
  RESTART_COUNTERS   F
  RESTART_CELL       F            ! keep the &CELL of THIS input: the NVT restart
                                  ! has no CELL_REF and would replace the whole &CELL
&END EXT_RESTART

&MOTION
  &MD
    ENSEMBLE    NPT_I
    STEPS       6000              ! 3 ps
    TIMESTEP    0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE   CSVR
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
    &BAROSTAT
      PRESSURE [bar] 1.01325
      TIMECON  [fs]  300.0
    &END BAROSTAT
  &END MD
  &PRINT
    &CELL ON                      ! water32_npt-1.cell: step, time, cell vectors, volume
      &EACH
        MD 10                     ! same stride as TRAJECTORY, needed by reftraj_npt.inp
      &END EACH
    &END CELL
    &STRESS ON
      &EACH
        MD 10
      &END EACH
    &END STRESS
    ! ... TRAJECTORY every 10 steps, RESTART every 20 steps ...
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL        ! mandatory for NPT
  &DFT
    &MGRID
      CUTOFF     600              ! raised from 400 Ry to reduce the Pulay stress
      REL_CUTOFF 60
    &END MGRID
    ! ... rest as in step 1 ...
  &END DFT
  &SUBSYS
    &CELL
      ABC 9.8654 9.8654 9.8654
      PERIODIC XYZ
      &CELL_REF
        ABC 10.4 10.4 10.4        ! fixes the FFT grid for the whole run
      &END CELL_REF
    &END CELL
    ! ... TOPOLOGY and KIND as in step 1 ...
  &END SUBSYS
&END FORCE_EVAL
```

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `ENSEMBLE` | `NPT_I` / `NPT_F` | 各向同性 / 全柔性晶胞 |
| `STRESS_TENSOR` | `ANALYTICAL` | 在 `&FORCE_EVAL` 中，默认 `NONE`。NPT 必须设置，Quickstep GPW 支持解析应力 |
| `&BAROSTAT/PRESSURE` | 目标压力 | 默认单位 bar，默认 0。可以用方括号写单位，如 `[bar]`、`[GPa]` |
| `&BAROSTAT/TIMECON` | 300–1000 fs | 默认 1000 fs。太短晶胞剧烈振荡，太长体积弛豫极慢 |
| `&BAROSTAT/VIRIAL` | `XYZ`（默认） | 仅 `NPT_F`：只让指定的笛卡尔方向弛豫，例如界面体系只放开 `Z` |
| `&BAROSTAT/TEMPERATURE` | 不写 | 压浴初始温度，默认等于系综温度 |
| `CUTOFF` | 比 NVT 高 30–50% | 降低 Pulay 应力；用 [第7章](07-geo-cell-opt.md) 中 CELL_OPT 的做法检查压力是否收敛 |
| `&CELL_REF` | 比初始晶胞大 3–5% | 让 FFT 网格固定在参考晶胞上，盒子胀大时有效截断能不会降低 |
| `&EXT_RESTART/RESTART_CELL` | `F`（从 NVT 重启文件起跑时） | 否则重启文件中不含 `CELL_REF` 的 `&CELL` 会整段替换本输入的 `&CELL`，`CELL_REF` 被静默丢弃。NPT 续 NPT 时保持默认 T，NPT 的重启文件会写出 `CELL_REF` |

检查密度：

```bash
# last column of the .cell file = volume in Angstrom^3; skip the first 1000 steps
awk '!/^#/ && $1>=1000 {n++; v+=$NF} END {V=v/n; printf "<V> = %.1f A^3, rho = %.4f g/cm^3\n", V, 32*18.015/0.602214076/V}' water32_npt-1.cell
grep "MD| Pressure \[bar\]" water32_npt_i.out | tail -n 3
```

密度换算 $\rho = N M/(N_A V)$：$V$ 以 Å³ 计、$\rho$ 以 g/cm³ 计时，$N_A\times10^{-24} = 0.6022$，即上面命令中的常数。

> [!WARNING]
> 96 个原子的盒子压力涨落可达数千 bar，瞬时压力本身没有意义，要看长时间平均。小盒子 NPT 的体积会剧烈振荡，至少要数 ps 才能得到有意义的平均密度。官方文档的建议是：小晶胞先单独做 CELL_OPT，再做定容 NVT；或者换更大的超胞后再做 NPT。PBE-D3 水的密度与实验 0.997 g/cm³ 存在偏差，结果以你的计算为准。

### 9.2 NPT_F：晶体与界面

把 `ENSEMBLE` 改成 `NPT_F`，晶胞的 6 个自由度（3 个长度、3 个夹角）都可以变化。适用于有限温度下的晶格常数、热膨胀、应力诱导相变。要点：

- 晶体先在 0 K 下 CELL_OPT（[第7章](07-geo-cell-opt.md)），再从优化后的晶胞开始 NPT_F；
- 晶胞形状变化时可能出现整体转动，通常不影响结果，分析时注意 `.cell` 中的矢量方向；
- 界面或 slab + 液体体系常用 `VIRIAL Z`，只让法向弛豫；真空 slab 不能做 NPT（真空方向没有压力定义）；
- 输出中的 `MD| Cell angles [deg]` 和 `.cell` 文件给出晶胞随时间的变化。

---

## 10. AIMD 质量检查清单

| 检查项 | 方法 | 合格标准（经验值） |
|---|---|---|
| 守恒量漂移 | NVE 中 `Cons Qty` 线性拟合斜率；或 `MD| Energy drift per atom` | ≲ 1 meV/atom/ps；没有阶跃 |
| 温度分布 | NVT 生产段温度直方图 | 平均值 ≈ 目标温度；标准差 ≈ $T_0\sqrt{2/N_f}$ |
| 能量均分 | `TEMP_KIND`，H 与 O 的平均温度 | 两者一致（差别在统计误差内） |
| SCF 行为 | 每步 SCF 迭代数、是否触发外层 SCF | 稳定在少数几次；不应逐渐变多 |
| 平衡 | 势能、温度、体积（NPT）的滑动平均 | 不再有趋势；丢弃平衡段后再统计 |
| 结构合理性 | 可视化最后一帧；最短原子间距 | 没有分子解离、原子逃逸或异常靠近 |
| 统计误差 | 分块平均（block averaging） | 各块平均值之差远小于你关心的效应 |

正则系综中动能服从自由度为 $N_f$ 的 $\chi^2$ 分布，温度的相对涨落为

$$
\frac{\sigma_T}{\langle T\rangle} = \sqrt{\frac{2}{N_f}},\qquad N_f = 3N - 3\ (\text{周期体系，扣除质心平动})
$$

本例 $N = 96$，$N_f = 285$，300 K 下 $\sigma_T \approx 25$ K。如果温度直方图明显比这个窄（例如用了 `TEMP_TOL` 硬重标）或明显更宽，说明恒温器设置有问题。

```bash
# histogram of the temperature (10 K bins) after equilibration
awk '!/^#/ && $1>=1000 {b=int($4/10)*10; h[b]++} END {for (k in h) print k, h[k]}' water32_prod-1.ener | sort -n
```

---

## 11. 轨迹后处理：RDF、MSD、VACF 与红外谱

本仓库的 `../scripts/md_analyze.py` 实现了下面的全部分析（读 `-1.ener`、`-pos-1.xyz`、`-vel-1.xyz` 和偶极文件），用 `python ../scripts/md_analyze.py -h` 查看用法。对本章示例的典型调用如下（在示例目录中运行；xyz 轨迹本身不含晶胞，定容轨迹用 `--abc` 给出，NPT 轨迹用 `--cell` 读 `.cell` 文件）：

```bash
S=../../scripts/md_analyze.py
python $S ener water32_nve-1.ener --natoms 96 -o nve_ener.png          # T, E_pot, Cons Qty drift
python $S rdf  water32_prod-pos-1.xyz --abc 9.8654 9.8654 9.8654 --pairs O-O O-H
python $S rdf  water32_npt-pos-1.xyz --cell water32_npt-1.cell --pairs O-O
python $S msd  water32_nve-pos-1.xyz --abc 9.8654 9.8654 9.8654 --elements O
python $S vacf water32_nve-vel-1.xyz --elements H O --max-lag-fs 1000  # VACF + VDOS
python $S ir   water32_nve_dipole.dat --dt-fs 1.0                      # dipole every 2 x 0.5 fs
```

下面给出原理和判断标准，便于你核对任何工具的结果。其他常用工具：TRAVIS（<http://www.travis-analyzer.de/>，对 AIMD 的 RDF、光谱分析非常全面）、MDAnalysis、ASE、VMD。

### 11.1 径向分布函数（RDF）

$$
g_{AB}(r) = \frac{V}{N_A N_B}\left\langle \sum_{i\in A}\sum_{j\in B,\,j\neq i} \frac{\delta(r - r_{ij})}{4\pi r^2}\right\rangle
$$

周期体系要用最小镜像距离，且 $r$ 最大只能取到 $L/2$。一个只依赖 numpy 的 O–O RDF 示例（立方盒子、定容轨迹）：

```python
import numpy as np

L = 9.8654                       # cubic box length (Angstrom)
frames = []
with open("water32_prod-pos-1.xyz") as f:
    while True:
        line = f.readline()
        if not line.strip():
            break
        n = int(line)
        f.readline()             # title line: i = ..., time = ..., E = ...
        rows = [f.readline().split() for _ in range(n)]
        sym = np.array([r[0] for r in rows])
        pos = np.array([[float(x) for x in r[1:4]] for r in rows])
        frames.append(pos[sym == "O"])

nbin, rmax = 200, L / 2
hist = np.zeros(nbin)
for o in frames:
    d = o[:, None, :] - o[None, :, :]
    d -= L * np.round(d / L)     # minimum image convention
    r = np.linalg.norm(d, axis=-1)[np.triu_indices(len(o), k=1)]
    hist += np.histogram(r, bins=nbin, range=(0.0, rmax))[0]

nO = len(frames[0])
edges = np.linspace(0.0, rmax, nbin + 1)
shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
ideal = len(frames) * nO / 2 * (nO - 1) / L ** 3 * shell   # ideal-gas pair count
r_mid = 0.5 * (edges[1:] + edges[:-1])
np.savetxt("rdf_OO.dat", np.c_[r_mid, hist / ideal], header="r[A] g_OO(r)")
```

参考：X 射线衍射给出的室温水 O–O 第一峰位置约 2.80 Å（Skinner 等，2013）。PBE 类 GGA 在 300 K 下普遍给出过高、过尖的第一峰（过度结构化），这是泛函本身的已知问题，并非输入错误。32 个水分子的盒子太小，RDF 只能看到约 4.9 Å 以内。

### 11.2 均方位移与扩散系数

$$
\mathrm{MSD}(t) = \left\langle |\mathbf r_i(t_0+t) - \mathbf r_i(t_0)|^2 \right\rangle_{i,t_0},\qquad
D = \lim_{t\to\infty}\frac{\mathrm{MSD}(t)}{6t}
$$

- CP2K 的 `-pos-1.xyz` 默认**不折回盒子**，可以直接计算 MSD；如果坐标被其他工具折回过，必须先展开；
- 只在 MSD 随时间线性增长的区间做拟合，跳过最初的弹道区（几十 fs）；
- 多个时间原点 $t_0$ 取平均能显著降低噪声；
- 周期盒子会使扩散系数偏小，Yeh–Hummer 有限尺寸校正为 $D_\infty = D_\mathrm{PBC} + \dfrac{\xi k_\mathrm B T}{6\pi\eta L}$，$\xi \approx 2.837$，$\eta$ 为剪切黏度；
- 实验值：298 K 水的自扩散系数约 2.3×10⁻⁹ m²/s（Mills 1973）。PBE 水的扩散明显偏慢。

REFTRAJ 自带 `&MSD` 子段（`MSD_PER_KIND`、`MSD_PER_MOLKIND` 等），可以在不重算能量（`EVAL NONE`）的情况下用 CP2K 直接算 MSD，详见输入参考。

### 11.3 速度自相关函数与振动态密度

$$
C_v(t) = \frac{\langle \mathbf v_i(0)\cdot\mathbf v_i(t)\rangle}{\langle \mathbf v_i(0)\cdot\mathbf v_i(0)\rangle},\qquad
\mathrm{VDOS}(\omega) \propto \int_{-\infty}^{\infty} C_v(t)\,e^{-i\omega t}\,dt
$$

VACF 的积分还给出扩散系数（Green–Kubo）：$D = \tfrac13\int_0^\infty \langle\mathbf v_i(0)\cdot\mathbf v_i(t)\rangle dt$。VDOS 中的峰对应体系的振动模式（水约 3400 cm⁻¹ 的 O–H 伸缩、约 1600 cm⁻¹ 的弯曲、低于 1000 cm⁻¹ 的平动/摆动）。按元素分别计算（只取 H 或只取 O）有助于指认峰。

- 速度至少每 1–2 步输出一次（§7）；
- 用 NVE 轨迹，或者 `REGION GLOBAL`、大 `TIMECON` 的 CSVR 轨迹；`MASSIVE` 恒温器会把谱线严重展宽；
- 频率分辨率 $\Delta\tilde\nu \approx 1/(c\,t_\mathrm{corr})$，相关时间窗取 0.5–1 ps 时约 33–67 cm⁻¹，再加窗函数（如 Hann）平滑。

### 11.4 偶极自相关与红外谱

红外吸收与总偶极 $\mathbf M(t)$ 的自相关函数相关（谐振近似下的量子修正因子已包含在系数中）：

$$
\alpha(\omega)\,n(\omega) \propto \omega^2 \int_{-\infty}^{\infty} \langle \mathbf M(0)\cdot\mathbf M(t)\rangle\, e^{-i\omega t}\,dt
= \int_{-\infty}^{\infty} \langle \dot{\mathbf M}(0)\cdot\dot{\mathbf M}(t)\rangle\, e^{-i\omega t}\,dt
$$

实际计算中常用右边的形式：先对偶极序列做有限差分得到 $\dot{\mathbf M}$，再做自相关和傅里叶变换。

- `water32_nve.inp` 通过 `&MOMENTS` `PERIODIC T` 每 1 fs 把 Berry 相总偶极写入 `water32_nve_dipole.dat`；
- 周期体系的 Berry 相偶极只确定到“极化量子” $e\mathbf L$ 的整数倍，相邻帧之间可能出现跳变。求导前要把跳变展开（相邻帧差值超过半个量子时加减一个量子）；
- 另一种做法是用最大局域 Wannier 函数中心（`&DFT/&LOCALIZE` 与 `&PRINT/&WANNIER_CENTERS`）构造分子偶极，可以得到分子分辨的谱，但成本更高；
- 谱线位置受泛函影响很大，PBE 的 O–H 伸缩峰一般红移，只做定性比较。

### 11.5 结果与文献比较时的注意事项

32 个水分子、几 ps 的 PBE-D3 轨迹，统计误差和有限尺寸效应都很大。与实验比较前要考虑：泛函误差（GGA 水过度结构化）、核量子效应（经典核忽略零点运动，常见做法是略微升高模拟温度或用路径积分 MD）、盒子大小和模拟时长。所有数值以你的计算为准，不要把本章给出的参考值当作你算出的结果。

---

## 12. 实战五：REFTRAJ 重算轨迹

`ENSEMBLE REFTRAJ` 不积分运动方程，而是逐帧读入一条已有轨迹，对每帧做指定的计算。典型用途：

- 用更高精度的方法（更大基组、杂化泛函）对 GGA 轨迹的若干帧重算能量和力，评估方法误差或为机器学习势（DeePMD、MACE 等）准备训练标签；
- 对轨迹计算额外的性质（偶极、电荷、cube 等，在 `&FORCE_EVAL/&DFT/&PRINT` 中打开）；
- 不做电子结构计算，只用 CP2K 计算 MSD 等结构量（`EVAL NONE`）。

### 12.1 定晶胞轨迹：`reftraj_eval.inp`

```fortran
# Syntax for CP2K >= 2026.2:   EVAL ENERGY_FORCES
#   (EVAL exists since 2024.1; the old logical keywords
#    EVAL_ENERGY_FORCES / EVAL_FORCES were removed in 2026.2.)
# Frames: 1, 101, 201, ..., 901 of water32_prod-pos-1.xyz (10 frames).
# LAST_SNAPSHOT must equal FIRST_SNAPSHOT + STRIDE * (number of strides).

&GLOBAL
  PROJECT     water32_reftraj
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE REFTRAJ
    STEPS    10                   ! ignored: reset to LAST-FIRST+1 = 901, every STRIDE-th -> 10 frames
    &REFTRAJ
      TRAJ_FILE_NAME water32_prod-pos-1.xyz
      FIRST_SNAPSHOT 1
      LAST_SNAPSHOT  901
      STRIDE         100
      EVAL           ENERGY_FORCES
    &END REFTRAJ
  &END MD
  &PRINT
    &FORCES ON
      &EACH
        MD 1
      &END EACH
    &END FORCES
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    ! ... basis/potential files and PBE-D3(BJ) &XC as before ...
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT   1.0E-12
      EXTRAPOLATION USE_PREV_WF   ! frames are 500 fs apart: only reuse the previous WFN
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      ! ... OT + OUTER_SCF as before, &PRINT &RESTART OFF ...
    &END SCF
  &END DFT
  &SUBSYS
    ! &CELL (ABC 9.8654, PERIODIC XYZ) and &TOPOLOGY (water32.xyz: defines atoms
    ! and kinds; positions come from the trajectory) as before
    &KIND H
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 关键字 | 推荐值 | 说明 |
|---|---|---|
| `TRAJ_FILE_NAME` | 轨迹文件 | 默认 `reftraj.xyz`。原子数和元素顺序必须与 `&SUBSYS` 一致 |
| `FIRST_SNAPSHOT` / `LAST_SNAPSHOT` / `STRIDE` | — | 帧编号从 1 开始。`LAST_SNAPSHOT` 默认 0，**必须设置**，且必须恰好等于 `FIRST_SNAPSHOT + STRIDE × 间隔数`，否则会读到文件末尾而报错 |
| `STEPS` | — | 设置了 `LAST_SNAPSHOT` 时被内部改写为 `LAST_SNAPSHOT − FIRST_SNAPSHOT + 1`，再每隔 `STRIDE` 取一帧，实际计算的帧数为 $(\text{LAST}-\text{FIRST})/\text{STRIDE}+1$；写不写都行 |
| `EVAL` | `NONE` / `ENERGY` / `ENERGY_FORCES` | 默认 `NONE`（只读结构，用于 MSD 等分析）。为 MLP 打标签用 `ENERGY_FORCES` |
| `WRAP` | `NONE`（默认） | `POSITIVE` / `CENTRAL` 把原子折回盒子；计算 MSD 时会被忽略 |
| `VARIABLE_VOLUME` / `CELL_FILE_NAME` | 变晶胞轨迹用 | 见 §12.2 |
| `EXTRAPOLATION` | `USE_PREV_WF` | 帧间隔很大，ASPC 的历史外推没有意义 |

运行与检查：

```bash
mpirun -np 32 cp2k.psmp -i reftraj_eval.inp -o reftraj_eval.out
cat water32_reftraj-1.ener                      # Pot. column = new energies, one line per frame
grep -c "i =" water32_reftraj-frc-1.xyz         # number of force frames, should be 10
```

REFTRAJ 中 `.ener` 每帧一行：第 1、2 列是从轨迹标题行（`i = ..., time = ...`）读到的原始步号和时间，即该帧在生产段中的 MD 步号，`Pot.` 列是新算出的能量，`Kin.` 和 `Temp` 为 0（REFTRAJ 不积分运动方程，没有核动能）。新的力在 `water32_reftraj-frc-1.xyz` 中。主输出中的 `MD| Energy drift per atom` 只在 `EVAL NONE` 时有意义（此时它是轨迹标题行中原能量相对第一帧的差，每原子，K）；`EVAL ENERGY`/`ENERGY_FORCES` 时这一行恒为 0，不要拿它判断什么。

**版本差异**：

| CP2K 版本 | 写法 |
|---|---|
| ≤ 2023.2 | `EVAL_ENERGY_FORCES T`（只算能量）；再加 `EVAL_FORCES T` 才算力 |
| 2024.1 – 2026.1 | 新写法 `EVAL NONE/ENERGY/ENERGY_FORCES` 与旧的逻辑关键字并存 |
| ≥ 2026.2 | 只认 `EVAL`；旧输入会报 `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ` |

旧关键字与新写法的对应关系：只有 `EVAL_ENERGY_FORCES T` → `EVAL ENERGY`；有 `EVAL_FORCES T` → `EVAL ENERGY_FORCES`。

### 12.2 变晶胞（NPT）轨迹：`reftraj_npt.inp`

NPT 轨迹每帧的晶胞不同，必须让 REFTRAJ 逐帧更新晶胞。有两种方法：

**方法 A（各版本通用）**：`VARIABLE_VOLUME T` + `CELL_FILE_NAME`，从 `.cell` 文件逐帧读入晶胞：

```fortran
&MOTION
  &MD
    ENSEMBLE REFTRAJ
    STEPS    11                   ! ignored: reset to LAST-FIRST+1 = 501, every STRIDE-th -> 11 frames
    &REFTRAJ
      TRAJ_FILE_NAME  water32_npt-pos-1.xyz
      CELL_FILE_NAME  water32_npt-1.cell
      VARIABLE_VOLUME T
      FIRST_SNAPSHOT  101
      LAST_SNAPSHOT   601
      STRIDE          50
      EVAL            ENERGY_FORCES
    &END REFTRAJ
  &END MD
&END MOTION
! &FORCE_EVAL: same electronic structure as water32_npt_i.inp (600 Ry, CELL_REF),
!   without STRESS_TENSOR, SCF_GUESS ATOMIC, EXTRAPOLATION USE_PREV_WF
```

前提是 NPT 运行中 `&TRAJECTORY` 和 `&CELL` 的输出间隔相同（示例中都是 `MD 10`），因为 CP2K 会逐帧核对两个文件中的步号。

**方法 B（2026.2 起）**：NPT 运行时把轨迹写成 `&TRAJECTORY FORMAT EXTXYZ`（2026.1 起可输出），标题行带 `Lattice="..."`。2026.2 的 REFTRAJ 能从 EXTXYZ 标题行直接读出每帧晶胞，此时不需要 `.cell` 文件；输入参考中也写明，只有在轨迹文件中找不到晶胞信息时才会去读 `CELL_FILE_NAME`。这对 ASE 等工具生成的 EXTXYZ 轨迹同样有用。

> [!NOTE]
> 为机器学习势准备数据时，常见做法是先用便宜的设置跑 AIMD 采样构型，再用 REFTRAJ 以高精度设置重算能量和力。`cp2kdata`（附带 dpdata 插件）可以把 CP2K 的 AIMD/REFTRAJ 输出转成 DeePMD 训练集；新版本输出格式变化后，先用少量数据验证解析结果是否正确。

---

## 13. 性能与资源

- 每个 MD 步的成本 ≈ 一次单点 SCF 的成本 × (外推后的 SCF 迭代数 / 从头 SCF 的迭代数)。先跑 20–50 步，读 `MD| CPU time per MD step` 估算总机时；
- 输出频率要克制：每步写速度、力和 wfn，I/O 可能占用明显的时间；
- MPI/OpenMP 混合并行、核数选择见 [第1章](01-install.md) 与作业脚本；AIMD 是强扩展问题，核数加倍不一定带来速度加倍，先做 2–3 个核数的小测试；
- 杂化泛函 AIMD 成本高得多，必须用 ADMM（[第9章](09-hybrid-admm.md)），并先用 GGA 平衡再切换；
- 金属或小带隙体系用 OT 往往难以收敛（2026.2 中 k 点计算也不支持 OT），一般改用对角化 + `&SMEAR` + 足够的 `ADDED_MOS`（[第5章](05-scf-strategy.md)）。Γ 点计算中 ASPC 外推仍可照常使用；k 点下的波函数外推是 2026.2 才加入的功能（见 2026.2 release notes），旧版本做 k 点 AIMD 时要留意每步的 SCF 迭代数。

---

## 参数选择与收敛建议

| 项目 | 建议 | 理由 |
|---|---|---|
| `TIMESTEP` | 含 H：0.5 fs；氘代或只含重元素：1–2 fs | 最快振动周期的 1/10–1/20 |
| `EPS_SCF` | 1E-6（生产），NVE 漂移大时 1E-7 | 力噪声决定能量漂移 |
| `EPS_DEFAULT` | 1E-12 | 积分筛选误差也是力噪声来源 |
| 外推 | `ASPC` + `EXTRAPOLATION_ORDER 3` | 默认且稳定；不稳定时降到 2 |
| SCF 算法 | 绝缘体：OT + `OUTER_SCF`；金属：对角化 + 展宽 | 每步 SCF 次数少且稳健 |
| `CUTOFF` | NVT/NVE 用单点收敛值；NPT 提高 30–50% 并加 `CELL_REF` | Pulay 应力 |
| 恒温器 | 平衡：CSVR 50–100 fs；生产：CSVR ≥ 200 fs，`REGION GLOBAL` | 官方文档推荐 CSVR 作日常首选 |
| Nose 链 | `LENGTH 3`，`TIMECON` 100–1000 fs | 仅在体系已基本平衡时使用 |
| 高温 / 退火 | 分几段逐步改变 `TEMPERATURE`，每段先短平衡；降温时 `TIMECON` 不要太短 | 温度突变容易让 SCF 失败或结构崩溃 |
| 动力学性质 | NVT 平衡后切 NVE，或用 `REGION GLOBAL`、大 `TIMECON` 的 CSVR | `MASSIVE`、Langevin 会改变动力学 |
| 压浴 | `NPT_I`（液体）/ `NPT_F`（晶体），`TIMECON` 300–1000 fs | 太小振荡，太大弛豫慢 |
| 输出频率 | 轨迹 5–10 步；动力学性质的速度/偶极 1–2 步；重启 20 步；`RESTART_HISTORY` 1000–2000 步 | 平衡分析需要与 I/O 成本 |
| 平衡段 | ≥ 1–2 ps（简单液体），更复杂体系更长 | 丢弃后再统计 |
| 生产段 | 结构性质 ≥ 5–10 ps；扩散、谱 ≥ 10–20 ps 或多条独立轨迹 | 统计误差 |
| 体系大小 | 大于待测性质的关联长度；RDF 只能看到 L/2 | 有限尺寸效应 |

---

## 常见问题与报错

**1. `SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`**
原因：某一步 SCF 在 `MAX_SCF × OUTER_SCF/MAX_SCF` 次迭代内没收敛，常见于开局结构差、温度过高、`EPS_SCF` 太严而 `MAX_SCF` 太小、金属体系用了 OT。解决：加大 `OUTER_SCF/MAX_SCF`；OT 预条件子换 `FULL_ALL`；金属改对角化 + 展宽；从 `RESTART_HISTORY` 中较早的重启点续算。不要设 `IGNORE_CONVERGENCE_FAILURE` 继续 MD。

**2. 温度一开始就飙升到数千 K，或者分子解离、原子飞出盒子。**
原因：初始结构中有原子重叠或极短的键（受力巨大），或者步长太大。解决：先做几何优化；检查 xyz 中是否有重复原子；含 H 体系用 0.5 fs；平衡段用短 `TIMECON` 的 CSVR 或 Langevin；必要时先在低温下短时间平衡再升温。

**3. `Virial evaluation not requested for this run in the input file! You may consider to switch on the virial evaluation with the keyword: STRESS_TENSOR. Be sure the method you are using can compute the virial!`**
原因：`NPT_I`/`NPT_F`/`NPE_*` 系综需要应力张量，但 `&FORCE_EVAL` 中没有 `STRESS_TENSOR`。解决：加 `STRESS_TENSOR ANALYTICAL`，同时提高 CUTOFF。

**4. NPT 中体积持续收缩，或密度明显偏高。**
原因：CUTOFF 太低引起的 Pulay 应力；或者盒子太小、压浴 `TIMECON` 太短导致剧烈振荡。解决：CUTOFF 提高到压力收敛的值（例如 600 Ry 以上）并加 `CELL_REF`；`TIMECON` 加大到 500–1000 fs；换更大的盒子。

**5. NVE 中守恒量持续漂移。**
原因：力噪声（`EPS_SCF` 太松、`EPS_DEFAULT` 太松、截断能不足）或步长太大。解决：按 §7 末尾的顺序逐项排查。

**6. 续算后轨迹或 `.ener` 中出现重复的步号。**
原因：作业被杀时最后几步已写入输出，但重启文件停在更早的步上。解决：分析前按步号去重（保留后一次运行的数据），或用 `cp2k-output-tools` 的 `xyz_restart_cleaner`；用 `WALLTIME` 让 CP2K 自己正常退出（必要时配合 `FORCE_LAST T`，见 §6）。

**7. 续算后 SCF 第一步迭代很多，或能量跳变。**
原因：`.wfn` 与当前体系不匹配（改了基组、赝势 q 值、原子顺序），或用的是另一次运行的 `.wfn`。解决：确认 `.restart` 和 `.wfn` 来自同一次运行；体系有变化时改回 `SCF_GUESS ATOMIC`。

**8. 2026.2 中运行旧 REFTRAJ 输入报 `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ`（或 `EVAL_FORCES`）。**
原因：2026.2 删除了这两个逻辑关键字。解决：改成 `EVAL ENERGY`（原来只有 `EVAL_ENERGY_FORCES T`）或 `EVAL ENERGY_FORCES`（原来有 `EVAL_FORCES T`）。`EVAL` 自 2024.1 起可用，新写法可以在 2024.1 及以后的版本中通用。

**9. `Reached the end of the Trajectory  frames in the TRAJECTORY file. Number of missing frames (n).`**
原因：`LAST_SNAPSHOT` 超出了轨迹中的帧数，或者不等于 `FIRST_SNAPSHOT + STRIDE × 间隔数`。解决：用 `grep -c "i =" traj.xyz` 统计帧数，按公式重设 `LAST_SNAPSHOT`。

**10. `Number of atoms for step ... in the trajectory file does not match the reference configuration` 或 `Atomic configuration from trajectory file does not match the reference configuration`。**
原因：REFTRAJ 读入的轨迹与 `&SUBSYS` 中的原子数或元素顺序不一致（例如可视化软件导出时重排了原子）。解决：`&TOPOLOGY` 中使用生成该轨迹时的同一个坐标文件；不要手工编辑轨迹中的元素符号。

**11. REFTRAJ 变晶胞重算时因步号不一致而中止。**
原因：`.cell` 与轨迹的输出间隔不同，两文件的帧对不上。CP2K 用断言（`CPASSERT`）核对两者的步号，所以报错信息只是一个指向 `integrator.F` 的 assertion failure，不会说明原因。解决：重新运行时让 `&CELL` 与 `&TRAJECTORY` 使用相同的 `&EACH/MD`；或者在 2026.2 中改用 EXTXYZ 轨迹（§12.2 方法 B）。如果 `FIRST_SNAPSHOT` 超过了轨迹（或 `.cell`）中的帧数，则会报 `Reached the end of the trajectory file for REFTRAJ. Number of steps skipped equal to the number of steps present in the file.`（`.cell` 文件对应 `... end of the cell file ...`）。

**12. 用 `REGION MASSIVE` 算出的扩散系数比 `GLOBAL` 小得多、VDOS 峰很宽。**
原因：每个自由度单独控温，等于加了强随机扰动，破坏了动力学关联。解决：动力学性质用 NVE 或 `GLOBAL` 的 CSVR（大 `TIMECON`）轨迹计算。

**13. 从 NVT 重启文件起跑 NPT，输入里写的 `&CELL_REF` 没有生效（NPT 的 `.restart` 中找不到 `&CELL_REF`）。**
原因：`&EXT_RESTART` 默认 `RESTART_CELL T`，会用重启文件中的 `&CELL` 整段替换输入中的 `&CELL`。NVT 运行没有 `CELL_REF`，它的重启文件里也就没有，于是你新写的 `CELL_REF` 被静默丢弃，没有任何报错。解决：在 `&EXT_RESTART` 中加 `RESTART_CELL F`（示例 `water32_npt_i.inp` 就是这样写的）。使用了参考晶胞时，CP2K 写出的 `.restart` 会包含 `&CELL_REF` 子段，可以用 `grep -c CELL_REF water32_npt-1.restart` 确认。

**14. 续算一开始就报 `maximum step number smaller than initial step value`。**
原因：`RESTART_COUNTERS T` 读入的起始步号已经 ≥ `MAX_STEPS`（例如重复提交了一个写有 `MAX_STEPS` 的续算输入，而上次已经跑完）。解决：调大 `MAX_STEPS` 或删掉它，只用 `STEPS` 控制本次步数。

---

## 练习

1. 把 `water32_nve.inp` 的 `TIMESTEP` 分别改为 0.25、0.5、1.0 fs，各跑 1 ps，比较守恒量漂移和振荡幅度，验证漂移大致按 $\Delta t^2$ 变化，并找出 1.0 fs 下漂移变大的原因。
2. 固定步长 0.5 fs，把 `EPS_SCF` 分别设为 1E-5、1E-6、1E-7，比较每步 SCF 迭代数、每步耗时和 NVE 漂移，找出性价比最高的设置。
3. 把 `EXTRAPOLATION` 改成 `USE_PREV_WF` 和 `PS`，与 `ASPC` 比较每步平均 SCF 迭代数。
4. 分别用 CSVR（`TIMECON` 50、200、1000 fs）和 Nose–Hoover 链（`water32_nvt_nose.inp`）跑生产段，比较温度分布的标准差与 §10 中的理论值 25 K，以及 O–O RDF 是否一致。
5. 用 NVE 轨迹分别通过 MSD 和 VACF 积分计算 H₂O 的自扩散系数，比较两者，并加上 Yeh–Hummer 校正后与实验值 2.3×10⁻⁹ m²/s 比较（水的剪切黏度可取实验值约 0.89 mPa·s）。
6. 用 `water32_nve_dipole.dat` 计算红外谱，指认 O–H 伸缩、H–O–H 弯曲和低频的摆动带，并与 VDOS 对比哪些峰在红外中更强。
7. 把 `water32_npt_i.inp` 的 CUTOFF 分别改为 400 和 600 Ry，比较 3 ps 内的平均密度，体会 Pulay 应力的影响。
8. 用 `reftraj_eval.inp` 得到 10 帧的 DZVP 与 TZV2P 能量差和力的均方根偏差，判断 DZVP 对这个体系是否够用。

---

## 参考资料

**官方文档**

- Molecular Dynamics（2026.2，含系综、恒温器、外推、平衡和验证清单）：<https://manual.cp2k.org/cp2k-2026_2-branch/methods/sampling/molecular_dynamics.html>
- 输入参考 `&MOTION/&MD`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD.html>
- `&THERMOSTAT`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD/THERMOSTAT.html>（子段 `CSVR`、`NOSE`、`GLE`、`AD_LANGEVIN`、`DEFINE_REGION` 在同目录下）
- `&BAROSTAT`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD/BAROSTAT.html>
- `&LANGEVIN`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD/LANGEVIN.html>
- `&REFTRAJ`（含 `EVAL`、`WRAP`、`VARIABLE_VOLUME`、`&MSD`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD/REFTRAJ.html>
- `&MOTION/&PRINT/&TRAJECTORY`（`FORMAT` 含 `EXTXYZ`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/PRINT/TRAJECTORY.html>
- `&EXT_RESTART`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/EXT_RESTART.html>
- `&QS`（`EXTRAPOLATION`、`EXTRAPOLATION_ORDER`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS.html>
- `&DFT/&PRINT/&MOMENTS`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/MOMENTS.html>
- Restarting：<https://manual.cp2k.org/trunk/methods/restarting.html>
- How to make a SCF run converge：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- CP2K 2026.2 release notes（REFTRAJ `EVAL` 替代旧关键字、EXTXYZ 读晶胞）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- CP2K 2026.1 release notes（EXTXYZ 轨迹输出）：<https://github.com/cp2k/cp2k/releases/tag/v2026.1>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/` 下大量 MD、REFTRAJ 回归测试输入（如 `tests/Fist/regtest-1-2/deca_ala_reftraj*.inp`），是查关键字真实写法的最佳来源
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，有 CP2K AIMD、REFTRAJ、constrained MD 与 DeePMD 联用的实践经验
- [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata)：CP2K 输出解析，带 dpdata 插件，可把 AIMD 数据转为 DeePMD 训练集
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：`cp2kparse` 解析输出；`xyz_restart_cleaner` 去除续算产生的重复帧
- [molmod/psiflow](https://github.com/molmod/psiflow)：以 CP2K 为参考计算的大规模 MLP 训练工作流
- [chenggroup/ai2-kit](https://github.com/chenggroup/ai2-kit)：以 CP2K 为标注引擎的势函数训练与作业管理工具链
- [MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis)：通用的轨迹分析库（RDF、MSD 等）

**本仓库脚本**：`../scripts/md_analyze.py`（`.ener` 统计与漂移、RDF、MSD、VACF、红外谱）、`../scripts/restart_update.py`（修改重启文件）、`../scripts/cp2k_outparse.py`（解析主输出）、`../scripts/job_cp2k_slurm.sh` 与 `../scripts/job_cp2k_pbs.sh`（作业脚本）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- G. Bussi, D. Donadio, M. Parrinello, *J. Chem. Phys.* **126**, 014101 (2007)，CSVR 恒温器：<https://doi.org/10.1063/1.2408420>
- G. J. Martyna, M. L. Klein, M. Tuckerman, *J. Chem. Phys.* **97**, 2635 (1992)，Nose–Hoover 链：<https://doi.org/10.1063/1.463940>
- G. J. Martyna, D. J. Tobias, M. L. Klein, *J. Chem. Phys.* **101**, 4177 (1994)，恒压 MD 算法：<https://doi.org/10.1063/1.467468>
- J. Kolafa, *J. Comput. Chem.* **25**, 335 (2004)，ASPC 外推：<https://doi.org/10.1002/jcc.10385>
- J. VandeVondele *et al.*, *Comput. Phys. Commun.* **167**, 103 (2005)，Quickstep：<https://doi.org/10.1016/j.cpc.2004.12.014>
- T. D. Kühne, M. Krack, F. R. Mohamed, M. Parrinello, *Phys. Rev. Lett.* **98**, 066401 (2007)，第二代 Car–Parrinello MD：<https://doi.org/10.1103/PhysRevLett.98.066401>
- R. Resta, *Phys. Rev. Lett.* **80**, 1800 (1998)，周期体系的位置算符与 Berry 相极化：<https://doi.org/10.1103/PhysRevLett.80.1800>
- I.-C. Yeh, G. Hummer, *J. Phys. Chem. B* **108**, 15873 (2004)，扩散系数的有限尺寸校正：<https://doi.org/10.1021/jp0477147>
- L. B. Skinner *et al.*, *J. Chem. Phys.* **138**, 074506 (2013)，室温水 O–O RDF 的 X 射线衍射测定：<https://doi.org/10.1063/1.4790861>
- R. Mills, *J. Phys. Chem.* **77**, 685 (1973)，水的自扩散系数：<https://doi.org/10.1021/j100624a025>
- M. J. Gillan, D. Alfè, A. Michaelides, *J. Chem. Phys.* **144**, 130901 (2016)，DFT 描述水的综述：<https://doi.org/10.1063/1.4944633>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](11-cube-charge-workfunction.md) | [返回目录](../README.md) | [下一章 →](13-neb-dimer.md)
