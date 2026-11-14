# 第15章 从头算分子动力学与机器学习力场（Ab Initio MD & Machine-Learned Force Fields in VASP）

本章讲解用 VASP 做从头算分子动力学（ab initio molecular dynamics, AIMD）：`IBRION = 0`、时间步长 `POTIM`、NVE/NVT/NpT 系综与恒温器（Nosé–Hoover、Langevin、CSVR 等，`MDALGO`）、Parrinello–Rahman 恒压，以及 XDATCAR 轨迹的 RDF、MSD/扩散系数、振动态密度分析。后半部分系统介绍 VASP 内置机器学习力场（machine-learned force field, MLFF，VASP ≥ 6.3）：`ML_MODE = train / refit / run / select / delta` 的完整流程、`ML_CTIFOR`、`ML_MB`、`ML_RCUT1/2`、`ML_ESTBLOCK` 等关键参数、ML_LOGFILE 的解读与验证方法，并简要比较 GRACE、MACE、NEP 等其他机器学习势。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据体系选择时间步长、系综与恒温器，写出 NVT 与 NpT 的完整 INCAR；
- 读懂 OSZICAR 中 MD 行的各列，判断体系是否已经平衡、能量是否守恒；
- 从 XDATCAR（定胞/变胞）计算 g(r)、配位数、MSD、自扩散系数和振动态密度；
- 走通 MLFF 的"训练 → 续训 → refit → 生产运行"流程，正确处理 `ML_ABN → ML_AB`、`ML_FFN → ML_FF` 的文件复制；
- 从 ML_LOGFILE 判断训练是否健康，在生产运行中用溢出因子（spilling factor）监控外推风险；
- 了解版本差异（`ML_ISTART` → `ML_MODE`，`ML_IERR` → `ML_ESTBLOCK`，6.6 的 delta 模式、GRACE 力场和溢出因子默认开启）。

| 需求 | 推荐方案 | 本章 |
|---|---|---|
| 几 ps、几十到一两百个原子，需要 DFT 精度 | AIMD（`IBRION = 0`） | §3–§5 |
| 几十到几百 ps、数百原子以上，或需要大量统计 | VASP MLFF：在线训练 + 快速预测模式 | §6–§8 |
| 已有大规模通用预训练模型，只做预测 | GRACE（VASP ≥ 6.6，实验性）或外部 MLIP | §9 |
| 想用便宜的泛函跑 MD、又想接近高精度泛函 | `ML_MODE = delta`（VASP ≥ 6.6.0） | §9.2 |
| 百万原子级、GPU 加速 | 用 VASP 产生训练数据，训练 NEP / MACE 等后在 GPUMD / LAMMPS 中运行 | §9.3 |

> [!NOTE]
> 本章默认你熟悉 [第3章](03-convergence.md) 的收敛测试和 [第4章](04-relax-eos.md) 的 Pulay 应力。MD 是"很多次单点计算"，每一步的 ENCUT、k 点、SCF 收敛都会直接影响轨迹质量。

---

## 2. 原理速览

### 2.1 Born–Oppenheimer 分子动力学

AIMD 的每一步都先做一次（或部分）电子自洽，得到 Hellmann–Feynman 力，再用 Verlet 类算法推进原子位置与速度。计算量 ≈（MD 步数）×（每步 SCF 时间），所以：

- 减少每步电子步数很重要：`NELMIN = 4` 保证最少的自洽步，`MAXMIX = 40` 复用上一步的混合历史（仅适用于普通 AIMD，**不要**用于 MLFF 训练）；
- 用 `vasp_gam` 做 Γ 点计算，`LREAL = Auto`，`PREC = Normal`，`ISYM = 0`；
- 定容 MD 中 ENCUT 取 ENMAX 附近即可；变胞（NpT）MD 必须至少提高约 30%，以抑制 Pulay 应力。

### 2.2 时间步长

时间步长应为体系中最快振动周期的 1/10–1/20。经验值：

| 体系 | `POTIM`（fs） | 说明 |
|---|---|---|
| 含 H（O–H、C–H、N–H 伸缩约 9–11 fs 周期） | 0.5（官方建议不超过约 0.7） | 或把 H 的质量改为氘（`POMASS`），可放宽到约 1 fs，但模拟的是重水/氘代体系的动力学 |
| 含 O、C、N 等轻元素，无 H | 1.0–1.5 | — |
| Si 及更重的元素 | 2–3 | 官方提到 Si 用 3 fs 也可以工作 |

检查方法：用 NVE（或很弱的恒温器）跑几百步，看 OSZICAR 中的守恒量 `E=` 是否漂移；漂移明显说明 `POTIM` 太大或 `EDIFF` 太松。

### 2.3 系综与恒温器（`MDALGO`）

| `MDALGO` | 恒温器 | 系综 | 必须设置的标签 | 备注 |
|---|---|---|---|---|
| 0 | Nosé–Hoover（旧实现） | NVT / NVE | `SMASS`（−3 为 NVE） | 已弃用：输出坐标总被折回盒子，无法算 MSD |
| 1 | Andersen | NVT / NVE | `ANDERSEN_PROB`（0.0 为 NVE） | 随机碰撞 |
| 2 | Nosé–Hoover | NVT | `SMASS` | 确定性；是 `MDALGO = 4` 链长为 1 的特例 |
| 3 | Langevin | NVT / NpT / NpH | `LANGEVIN_GAMMA`（每元素一个值，ps⁻¹）；NpT 另需 `LANGEVIN_GAMMA_L`、`PMASS` | 随机性强、遍历性好，训练 MLFF 的首选 |
| 4 | Nosé–Hoover 链 | NVT | `NHC_NCHAINS`（1–5）、`NHC_PERIOD`（以 MD 步数计，默认 40） | `NHC_NCHAINS = 0` 为 NVE |
| 5 | CSVR（Bussi 速度重标度） | NVT | `CSVR_PERIOD` | VASP ≥ 6.4.3；`CSVR_PERIOD = 0` 为 NVE |
| 13 | 多个 Andersen 恒温器 | NVT | `NSUBSYS`、`TSUBSYS`、`PSUBSYS` | 不同子体系不同温度 |

所有 `MDALGO` 选项都需要 VASP 以 `-Dtbdyn` 编译（5.4.4 之后的官方 makefile.include 模板默认包含）。

其他常用 MD 标签：

| 标签 | 含义 |
|---|---|
| `IBRION = 0` | 分子动力学 |
| `NSW` | MD 步数 |
| `POTIM` | 时间步长（fs）；`IBRION = 0` 时**没有默认值，必须设置** |
| `TEBEG` / `TEEND` | 起始 / 终止温度（K），线性升温；未提供初速度时按 TEBEG 的 Maxwell–Boltzmann 分布随机生成 |
| `ISIF` | 0/2 为定容（2 会计算应力，可监控压力）；3 为变胞（NpT，仅 `MDALGO = 3`）；MD 中不支持 4 和 7 |
| `NBLOCK` / `KBLOCK` | 每 NBLOCK 步写一次 XDATCAR 并累计对关联函数；每 NBLOCK×KBLOCK 步把平均后的 PCDAT 写出 |
| `RANDOM_SEED` | 三个整数，固定随机数种子以复现轨迹（初速度、Langevin 随机力） |
| `PSTRESS` | NpT 的目标压力（kB），1 bar = 0.001 kB |

> [!NOTE]
> 温度的定义与自由度有关：Nosé–Hoover 与 Andersen 守恒质心动量，自由度为 3N − 3（减去约束）；Langevin 不守恒质心，自由度多 3 个。OUTCAR 中打印的温度已按此处理。

### 2.4 NpT：Parrinello–Rahman + Langevin

`MDALGO = 3` + `ISIF = 3` 使用 Parrinello–Rahman 方法让晶格自由度也参与动力学，并对原子和晶格分别施加 Langevin 恒温：

- `LANGEVIN_GAMMA`：原子摩擦系数（ps⁻¹），每个 POSCAR 元素一个值；
- `LANGEVIN_GAMMA_L`：晶格摩擦系数（ps⁻¹）；
- `PMASS`：晶格的虚拟质量（amu，默认 1000）。太大则晶格响应很慢、采样低效；太小则每步晶格变化过大可能出现数值问题。官方给出立方盒子的晶格振动频率 ω_L = √(3BL/W)（B 为体模量、L 为盒长、W 为 PMASS），盒长加倍时 PMASS 也应加倍以保持 ω_L 不变；
- `PSTRESS`：目标外压（kB）。

> [!WARNING]
> 液体等缺乏长程有序的体系在完全自由的 NpT 中可能发生不可逆的盒子剪切变形（盒子被"压扁"成一层原子）。必须用 `ICONST` 文件约束晶格形状，或对正交盒子使用 `LATTICE_CONSTRAINTS = .FALSE. .FALSE. .TRUE.`（只让一个方向像活塞一样伸缩）。

### 2.5 MD 的输出文件

| 文件 | 内容 |
|---|---|
| OSZICAR | 每个 MD 步一行：`T=` 温度、`E=` 守恒总能（含离子动能和恒温器能量）、`F=` 势能面上的自由能（即 TOTEN）、`E0=`、`EK=` 离子动能、`SP=` / `SK=` 恒温器势能/动能 |
| XDATCAR | 每 NBLOCK 步的分数坐标；变胞（ISIF = 3）时每帧前重复写出晶格 |
| OUTCAR | 每步的受力、应力、压力（MD 中含离子动能贡献）、温度 |
| REPORT | `MDALGO ≥ 1` 时的详细 MD 信息（恒温器、约束、随机数种子等） |
| PCDAT | 对关联函数（NBLOCK×KBLOCK 步平均） |
| CONTCAR | 最后的结构**和速度**，可直接作为下一段 MD 的 POSCAR |
| vaspout.h5 | HDF5 输出（需编译开启），py4vasp 读取能量、结构、对关联函数 |

---

## 3. 实战一：液态 Si 的 NVT AIMD

完整输入在 [../examples/15-aimd-si-liquid-nvt/README.md](../examples/15-aimd-si-liquid-nvt/README.md)。

### 3.1 结构准备

- 超胞：64 个 Si 原子（2×2×2 金刚石惯用胞的原子排列），只用 Γ 点 → `vasp_gam`；
- 密度：盒长 10.511 Å，对应液态 Si 在熔点附近的密度约 2.57 g/cm³（近似值）。NVT 的结果依赖于所选密度，若密度不确定，应先做 NpT 平衡（§5）；
- 策略：先在 3000 K 熔化（避免从完美晶体过热不足而残留晶序），再降到 1800 K 做生产运行。

超胞也可以用 ASE 或 VASPKIT 生成：

```python
from ase.build import bulk
from ase.io import write
atoms = bulk("Si", "diamond", a=5.2557, cubic=True).repeat((2, 2, 2))   # 64 atoms
write("POSCAR", atoms, format="vasp", direct=True, sort=True)
```

### 3.2 第 1 段：熔化（INCAR.melt）

```bash
# Stage 1: melt the 64-atom crystal at 3000 K (NVT, Nose-Hoover)
SYSTEM  = Si64 AIMD stage 1 melt 3000 K

# electronic (fixed volume: ENCUT ~ ENMAX is acceptable)
PREC    = Normal
ENCUT   = 250
EDIFF   = 1E-5
ALGO    = Fast
LREAL   = Auto
ISMEAR  = 0
SIGMA   = 0.10         # liquid Si is metallic
NELMIN  = 4            # at least 4 SCF steps per MD step
MAXMIX  = 40           # reuse mixer history between MD steps (plain AIMD only)
ISYM    = 0            # no symmetry in MD

# molecular dynamics
IBRION  = 0
MDALGO  = 2            # Nose-Hoover thermostat
SMASS   = 1.0          # Nose mass parameter (see VASP wiki SMASS / NVT ensemble)
ISIF    = 2            # compute stress (pressure) but keep the cell fixed
TEBEG   = 3000
TEEND   = 3000
POTIM   = 2.0          # fs; Si is heavy (no H)
NSW     = 1500         # 3 ps
NBLOCK  = 1            # write every step to XDATCAR
RANDOM_SEED = 248489752 0 0

# output
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 取值 | 为什么 |
|---|---|---|
| `IBRION` | 0 | 分子动力学 |
| `MDALGO` / `SMASS` | 2 / 1.0 | Nosé–Hoover NVT；SMASS 控制恒温器"质量"（温度振荡频率），应与体系的特征振动频率相近 |
| `POTIM` | 2.0 | Si 较重，2 fs 安全；含 H 体系改为 0.5 |
| `TEBEG`/`TEEND` | 3000 | 熔化阶段用较高温度 |
| `ISIF` | 2 | 定容；计算应力以便监控压力 |
| `ENCUT` | 250 | 定容 MD 中约 ENMAX 即可 |
| `EDIFF` | 1E-5 | MD 中常用；能量漂移明显时收紧到 1E-6 |
| `NELMIN` / `MAXMIX` | 4 / 40 | 保证每步最少自洽次数、复用混合历史以减少电子步 |
| `ISMEAR`/`SIGMA` | 0 / 0.1 | 液态 Si 为金属性；也可用 Fermi 展宽（`ISMEAR = -1`，SIGMA 取 k_BT） |
| `ISYM` | 0 | MD 中关闭对称性 |
| `NBLOCK` | 1 | 每步都写 XDATCAR（计算 VDOS 需要）；长轨迹可增大 |
| `RANDOM_SEED` | 固定 | 初速度可复现 |
| `LWAVE`/`LCHARG` | `.FALSE.` | 避免每步写大文件 |

### 3.3 第 2 段：生产运行（INCAR.prod）

只修改温度和步数，其余保持不变：

```bash
TEBEG   = 1800
TEEND   = 1800
NSW     = 5000         # 10 ps; discard the first 1-2 ps as equilibration
```

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR      # 仓库不提供 POTCAR
cp INCAR.melt INCAR && mpirun -np 32 vasp_gam > melt.out
cp CONTCAR POSCAR                                     # CONTCAR 带速度，接着跑
cp INCAR.prod INCAR && mpirun -np 32 vasp_gam > prod.out
```

### 3.4 运行监控与平衡判断

```bash
grep "T=" OSZICAR | tail -5                              # 最近几步的温度和能量
awk '/T=/{print $1, $3, $5, $7}' OSZICAR > thermo.txt    # 步数 温度 守恒能E 势能F
grep -E "external pressure|total pressure" OUTCAR | tail -3
grep -cE "^(DAV|RMM):" OSZICAR                           # 电子步总数，除以 MD 步数得到平均值
```

OSZICAR 中一行 MD 输出形如（示意输出，数值仅供参考）：

```text
   500 T=  1795. E= -.33012345E+03 F= -.34190123E+03 E0= -.34189987E+03  EK= 0.14612E+02 SP= 0.27E+01 SK= 0.23E+00
```

判断平衡与质量：

- 温度在 TEBEG 附近涨落，没有系统性漂移；64 原子的相对温度涨落约 √(2/(3N)) ≈ 10%，属于正常；
- 势能 F 的滑动平均不再变化（熔化或相变后会出现台阶）；
- 守恒量 E 在 NVT 中应近似守恒（Nosé–Hoover 的扩展哈密顿量），长时间漂移说明 POTIM 太大或 SCF 不够准；
- 用块平均（block averaging）估计统计误差：把轨迹分成若干块分别求平均，看块平均的方差随块长是否收敛。

---

## 4. 轨迹分析：RDF、MSD、扩散系数与振动态密度

### 4.1 XDATCAR 与坐标展开

XDATCAR 中是分数坐标。`MDALGO ≥ 1` 时坐标通常不被折回盒子，但为保险起见，分析脚本应按最小像约定对相邻两帧的位移做"展开"（unwrap）：Δs = s(t+1) − s(t) − round(Δs)，再累加得到连续轨迹。这要求相邻两帧之间原子位移小于半个盒长（NBLOCK 不要太大）。`MDALGO = 0` 总是把坐标折回盒子，官方明确指出这会使 MSD 无法计算。

### 4.2 径向分布函数 g(r) 与配位数

$$
g_{AB}(r) = \frac{V}{N_A N_B}\left\langle \sum_{i\in A}\sum_{j\in B} \frac{\delta(r - r_{ij})}{4\pi r^2} \right\rangle,\qquad
n_{AB}(r_c) = 4\pi\rho_B\int_0^{r_c} g_{AB}(r)\,r^2\,dr
$$

配位数取到 g(r) 第一个极小值处。截断半径不应超过最短晶胞高度的一半。

### 4.3 均方位移与扩散系数

$$
\mathrm{MSD}(t) = \Big\langle \big|\mathbf{r}_i(t_0+t)-\mathbf{r}_i(t_0)\big|^2 \Big\rangle_{i,t_0},\qquad
D = \lim_{t\to\infty}\frac{\mathrm{MSD}(t)}{6t}
$$

- 对多个时间原点 t0 取平均可显著降低噪声；
- 只拟合 MSD 的线性区：跳过开头的弹道区（MSD ∝ t²），也不要用接近总时长的长时间尾部（统计点太少）；
- 单位：1 Å²/fs = 0.1 cm²/s；
- 周期性边界下 D 有有限尺寸效应，液体可用 Yeh–Hummer 修正 D_∞ = D_PBC + 2.837297 k_BT/(6πηL)（η 为黏度，L 为盒长）；
- 固体中 MSD 很快饱和到一个平台（振动幅度），只有发生扩散（如离子导体、液体）时才线性增长。

### 4.4 振动态密度（VDOS）

由速度自相关函数（VACF）的傅里叶变换得到：

$$
C_{vv}(t) = \frac{\langle \mathbf{v}(t_0+t)\cdot\mathbf{v}(t_0)\rangle}{\langle \mathbf{v}(t_0)\cdot\mathbf{v}(t_0)\rangle},\qquad
g(\nu) \propto \int C_{vv}(t)\cos(2\pi\nu t)\,dt
$$

XDATCAR 不含速度，可以用展开后坐标的中心差分得到速度（需要 `NBLOCK = 1`）。VDOS 包含全部非谐效应，可与 [第12章](12-phonon.md) 的谐振声子 DOS 对比。

### 4.5 分析工具

示例脚本 `analyze_md.py`（只依赖 numpy、matplotlib；支持定胞和变胞 XDATCAR）：

```bash
python3 analyze_md.py --potim 2.0 --skip 500 --pair Si-Si           # 温度/能量、g(r)、MSD、D
python3 analyze_md.py --potim 2.0 --skip 500 --vdos                 # 再加振动态密度
python3 analyze_md.py --potim 0.5 --skip 2000 --pair O-O --species O  # 水：O-O g(r)，O 的扩散
```

输出格式（示意输出，数值仅供参考）：

```text
5000 frames, 64 atoms (Si64)
analysed frames: 4500  (dt between frames = 2 fs, total 9.00 ps); <V> = 1161.37 A^3
<T> = 1801.3 +- 180.2 K (after skip)
RDF Si-Si: first peak at 2.4x A, first minimum ~3.0x A, coordination n(r_min) = 6.x
MSD Si: D = x.xxe-04 cm^2/s  (fit 1.80-5.40 ps)
```

液态 Si 的参考量级（文献值，近似）：g(r) 第一峰约 2.4–2.5 Å，配位数约 6，熔点附近 D 为 10⁻⁴ cm²/s 量级。

其他工具：

| 工具 | 用法要点 |
|---|---|
| 仓库脚本 [../scripts/md_analysis.py](../scripts/md_analysis.py) | 温度/能量/体积曲线、RDF、MSD、扩散系数、ML_LOGFILE 误差曲线；用 `python ../scripts/md_analysis.py -h` 查看用法 |
| py4vasp | `calc.energy[:].plot("temperature")` 画温度随步数变化；`calc.structure[:]` 取整条轨迹；`calc.pair_correlation.plot()` 画 VASP 计算的对关联函数 |
| ASE | `ase.io.read("XDATCAR", index=":", format="vasp-xdatcar")` 读入整条轨迹，可再写成 extxyz 交给 OVITO、MDAnalysis 等 |
| VASPKIT | 721 计算 MSD（需要参考结构 POSCAR.ref），722 计算 RDF |
| OVITO | 可视化轨迹、配位分析、共近邻分析 |

```python
from py4vasp import Calculation
calc = Calculation.from_path(".")
calc.energy[:].plot("temperature")          # 选择关键字可用 calc.energy.selections() 查看
print(calc.energy.selections())
traj = calc.structure[:]                     # 整条轨迹
calc.pair_correlation.plot()
```

---

## 5. 实战二：液态水的 NpT 分子动力学

完整输入在 [../examples/15-aimd-water-npt/README.md](../examples/15-aimd-water-npt/README.md)，其中 `build_water_box.py` 可生成任意大小的水盒子（fcc 格点 + 随机取向）。

### 5.1 INCAR

```bash
SYSTEM  = (H2O)32 NpT Langevin 300 K PBE-D3(BJ)

# electronic
PREC    = Normal
ENCUT   = 550          # >= 1.3 x ENMAX(O)=400 eV: the basis is fixed while the cell changes
EDIFF   = 1E-5
ALGO    = Fast
LREAL   = Auto
ISMEAR  = 0
SIGMA   = 0.05
NELMIN  = 4
MAXMIX  = 40
ISYM    = 0
IVDW    = 12           # DFT-D3 with Becke-Johnson damping

# molecular dynamics, NpT
IBRION  = 0
MDALGO  = 3            # Langevin thermostat
ISIF    = 3            # variable cell -> NpT (Parrinello-Rahman)
TEBEG   = 300
TEEND   = 300
POTIM   = 0.5          # fs; O-H stretching period is ~9 fs
NSW     = 6000         # 3 ps; extend for production
LANGEVIN_GAMMA   = 10.0 10.0   # ps^-1, one value per species (O H), POSCAR order
LANGEVIN_GAMMA_L = 10.0        # ps^-1, lattice degrees of freedom
PMASS   = 1000         # fictitious lattice mass (amu); test 100-1000
PSTRESS = 0.0          # target pressure in kB (1 bar = 0.001 kB)
NBLOCK  = 1
RANDOM_SEED = 486686595 0 0

# output
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 说明 |
|---|---|
| `MDALGO = 3` + `ISIF = 3` | Langevin + Parrinello–Rahman，VASP 中实现 NpT 的方式 |
| `LANGEVIN_GAMMA` | 个数必须等于 POSCAR（POTCAR）中的元素种类数；设为 0 的元素不受恒温器作用 |
| `LANGEVIN_GAMMA_L`、`PMASS` | 晶格的摩擦系数与虚拟质量；不设置时默认值为 0（摩擦），会得到不同的系综 |
| `ENCUT = 550` | 变胞 MD 的基组在运行中固定，体积变化会引入 Pulay 应力；官方建议比定容计算至少高 30%，并经常重启 |
| `IVDW = 12` | 纯 GGA 水明显过度结构化、密度偏低，色散修正是常见改进（结果仍依赖泛函） |
| `POTIM = 0.5` | 含 H |

### 5.2 用 ICONST 保持立方盒子

```text
LA 1 2 0
LA 1 3 0
LA 2 3 0
LR 1 0
LR 2 0
LR 3 0
S  1  0  0  0  0  0  0
S  0  1  0  0  0  0  0
S  0  0  1  0  0  0  0
S  0  0  0  1 -1  0  0
S  0  0  0  1  0 -1  0
S  0  0  0  0  1 -1  0
```

前 6 行定义原始坐标（3 个晶格夹角 `LA`、3 个晶格长度 `LR`），最后一列 0 表示"约束"。一旦定义了复合坐标 `S`，原始坐标本身只作为基，不再单独约束；这里 6 个 `S` 行分别约束三个夹角不变、三个边长两两相等（a1 − a2 = 0 等），于是盒子保持立方、只允许体积变化。`ICONST` 只能用于 `IBRION = 0`。

### 5.3 检查与注意事项

```bash
grep "volume of cell" OUTCAR | awk 'NR>1{print NR-1, 957.27/$5}' > density.dat   # 32 H2O: M = 576.48 amu
tail density.dat
python3 ../15-aimd-si-liquid-nvt/analyze_md.py --potim 0.5 --skip 2000 --pair O-O --species O
```

- 密度：平衡后在某个值附近涨落，32 分子体系涨落约 ±0.03 g/cm³；体积若持续单调漂移，先检查 ENCUT（Pulay 应力）与 PMASS。
- 结构：O–O g(r) 第一峰约 2.7–2.8 Å（实验约 2.8 Å，近似值）。
- 动力学：`LANGEVIN_GAMMA = 10` ps⁻¹ 有利于快速平衡，但会人为减慢扩散；计算扩散系数时，平衡后换成较弱的摩擦或 Nosé–Hoover / CSVR 的 NVT。
- 计算量：96 原子、550 eV、0.5 fs 步长的 AIMD 每 ps 需要数十个节点小时，这正是下面引入 MLFF 的动机。

---

## 6. 机器学习力场（MLFF）原理速览

### 6.1 方法要点

VASP 的 MLFF（Jinnouchi, Karsai, Kresse 等）是基于核函数的势：

- **描述符**：每个原子周围的原子密度展开为径向（两体，截断 `ML_RCUT1`，默认 8 Å）与角向（三体，截断 `ML_RCUT2`，默认 5 Å）描述符；
- **能量模型**：原子能量是当前局域环境与一组"局域参考构型"（local reference configurations）之间核函数的线性组合；
- **贝叶斯线性回归**：拟合的同时给出预测误差估计（Bayesian error estimate, BEE）；
- **在线学习（on-the-fly）**：MD 的每一步先用力场预测能量、力及其误差；若最大力误差超过阈值 `ML_CTIFOR`，就做一次 DFT，把这一结构加入训练集，并从中挑选新的局域参考构型，必要时重新拟合力场；否则直接用力场推进 MD。随着训练进行，DFT 调用越来越少；
- **稀疏化**：用类似 CUR 的算法去掉近似线性相关的局域参考构型（`ML_EPS_LOW` 控制程度）；
- **refit**：训练结束后用 SVD 正则化重新拟合（`ML_MODE = refit` 自动设置 `ML_IALGO_LINREG = 4` 与 `ML_LFAST = .TRUE.`），得到支持"快速预测模式"的力场，官方给出的速度提升约 20–100 倍，代价是不能再给出贝叶斯误差估计；
- **溢出因子（spilling factor）**：衡量当前原子环境与局域参考构型的相似度，取值 0（被很好描述）到 1（外推）；它不依赖贝叶斯误差，可以在快速模式的生产运行中监控外推风险。

### 6.2 文件

| 文件 | 输入/输出 | 内容 |
|---|---|---|
| `ML_AB` | 输入 | 训练数据库（结构、DFT 能量、力、应力，及局域参考构型列表）——续训或 refit 时由上一轮的 `ML_ABN` 复制而来 |
| `ML_ABN` | 输出 | 本次运行结束时的完整训练数据库 |
| `ML_FF` | 输入 | 力场（`ML_MODE = run`/`delta` 时读取）——由 `ML_FFN` 复制而来 |
| `ML_FFN` | 输出 | 本次运行生成的力场（二进制，带 ASCII 文件头，可查看是否支持快速模式） |
| `ML_LOGFILE` | 输出 | 设置摘要、内存估计、每步的 STATUS / ERR / BEEF / LCONF 等 |
| `ML_REG` | 输出 | 训练集上 DFT 与 MLFF 预测值的逐项对比 |

### 6.3 `ML_MODE` 与版本

| `ML_MODE`（≥ 6.4.0） | 旧版 `ML_ISTART`（6.3） | 作用 |
|---|---|---|
| `train` | 0（无 ML_AB）/ 1（有 ML_AB） | 在线训练：没有 ML_AB 时从零开始，有 ML_AB 时在其基础上续训 |
| `select` | 3 | 从 ML_AB 中的数据重新挑选局域参考构型（可用于外部或预先计算的 DFT 数据集） |
| `refit` | 4 | 基于 ML_AB 重新拟合，生成快速预测模式的力场（**强烈推荐**在生产运行前做） |
| `refitbayesian` | 4 | 用贝叶斯回归 refit，可保留误差估计但更慢更不准（已弃用） |
| `run` | 2 | 只用 ML_FF 预测，不做 DFT、不学习 |
| `delta` | — | VASP ≥ 6.6.0：MD 由"DFT 结果 + ML_FF 预测"之和驱动，力场学习两个理论层次之间的差值 |
| `none` | — | 忽略该标签 |

相关版本变化：

- ML 力场自 VASP 6.3.0 起可用（`ML_LMLFF`）；`ML_MODE` 自 6.4.0 起取代 `ML_ISTART`；
- 快速预测模式（`ML_LFAST`）、`ML_OUTBLOCK`、`ML_OUTPUT_MODE` 自 6.4.0 起；`ML_OUTBLOCK` 的实现到 6.4.3 才完整，生产运行建议 ≥ 6.4.3；
- 溢出因子自 6.5.0 引入；6.5.1 起 `ML_IERR` 更名为 `ML_ESTBLOCK`（旧名仍可用）；
- 6.6.0：溢出因子始终开启，超过 0.9 时计算终止；`ML_OUTBLOCK` 默认 10；新增 `delta` 模式、GRACE 力场接口、ML_LOGFILE 中的分元素误差与 NORME（归一化误差）行。

---

## 7. 实战三：Si 的 MLFF 全流程（train → refit → run）

完整输入在 [../examples/15-mlff-si-train-run/README.md](../examples/15-mlff-si-train-run/README.md)。目标：得到一个能在 800 K 附近做晶体 Si NpT 生产运行的力场。

### 7.1 第 1 步：在线训练（INCAR_train）

```bash
# Stage 1: on-the-fly training (VASP >= 6.4 for ML_MODE; 6.3 used ML_ISTART)
SYSTEM  = Si64 MLFF on-the-fly training NpT 200-1000 K

# ab initio part: identical in every training/continuation run!
PREC    = Normal
ENCUT   = 350          # NpT: >= 1.3 x the fixed-volume cutoff (ENMAX Si = 245 eV)
EDIFF   = 1E-6         # forces must be accurate: they are the training labels
ALGO    = Fast
LREAL   = Auto
ISMEAR  = 0
SIGMA   = 0.05
NELM    = 100
ISYM    = 0
# do NOT set MAXMIX > 0 in MLFF training
LWAVE   = .FALSE.
LCHARG  = .FALSE.

# molecular dynamics
IBRION  = 0
MDALGO  = 3            # Langevin: stochastic, good phase-space sampling
ISIF    = 3            # NpT: cell fluctuations make the force field more robust
TEBEG   = 200
TEEND   = 1000
POTIM   = 2.0
NSW     = 10000
LANGEVIN_GAMMA   = 10.0
LANGEVIN_GAMMA_L = 10.0
PMASS   = 1000
RANDOM_SEED = 486686595 0 0

# machine-learned force field
ML_LMLFF = .TRUE.
ML_MODE  = train       # starts from scratch without ML_AB, continues if ML_AB exists
```

训练设置背后的"最佳实践"（官方 Best practices 页面的要点）：

| 建议 | 原因 |
|---|---|
| 升温训练：从较低温度升到比目标温度高约 30% | 探索更大的相空间，力场更稳健；但不要远高于目标温度（误差随温度增大） |
| 优先用 NpT（`ISIF = 3`），液体加形状约束 | 晶胞涨落提高鲁棒性；液体需防止盒子塌缩 |
| NVT 训练用 Langevin，避免 NVE | Langevin 的随机性有利于遍历；NVE 采样太窄 |
| `POTIM`：含 H ≤ 0.7 fs，含 O ≤ 1.5 fs，Si 等重元素可到 3 fs | 与 AIMD 相同的稳定性要求 |
| 不要设 `MAXMIX > 0` | DFT 之间可能隔了成百上千步，复用混合历史常导致 SCF 不收敛 |
| 关闭对称性 `ISYM = 0` | 与普通 MD 相同 |
| NpT 时 ENCUT 至少比定容高 30%，并经常重启续训 | 重启会重新初始化 PAW 基组，减小 Pulay 应力 |
| 训练与续训之间**不改** DFT 参数和 POTCAR | 力场只能拟合一个势能面，混入不同设置的数据会破坏一致性 |
| 表面、孤立分子：`ML_WTSIF` 设为很小的值（如 1E-10） | 真空方向没有应力，不应训练应力 |
| 复杂体系分步训练：先体相，再表面，再分子+表面 | 大量 DFT 计算发生在便宜的子体系上 |
| SCF 未在 NELM 步内收敛的结构不会被加入训练集 | 这是安全措施；频繁出现说明 SCF 设置需要调整 |

运行：

```bash
mkdir 1_train && cp INCAR_train 1_train/INCAR && cp POSCAR KPOINTS POTCAR 1_train/
cd 1_train && mpirun -np 32 vasp_std > vasp.out
```

### 7.2 监控 ML_LOGFILE

```bash
grep "Total memory consumption" ML_LOGFILE                     # 启动后先确认内存估计
grep STATUS ML_LOGFILE | grep -v "#" | awk '{print $3}' | sort | uniq -c   # 各状态的步数
grep STATUS ML_LOGFILE | grep -E "learning|critical" | grep -v "#" | tail
grep ERR ML_LOGFILE | grep -v "#" | tail -3                    # 训练集 RMSE：能量(eV/atom) 力(eV/A) 应力(kB)
grep BEEF ML_LOGFILE | grep -v "#" | awk '{print $2, $4, $6}' | tail   # 步数、最大力 BEE、当前 ML_CTIFOR
grep LCONF ML_LOGFILE | grep -v "#" | tail -1                  # 各元素的局域参考构型数
```

ML_LOGFILE 的数据行都以关键字开头，列名写在文件开头的注释头里：

| 关键字 | 列（按官方文件头） | 怎么看 |
|---|---|---|
| `STATUS` | nstep, state, is, doabin, iff, nsample, ngenff | state：`accurate`（力场足够准）、`threshold`（超阈值，采集 DFT 数据）、`learning`（用新数据重训）、`critical`（误差过大，强制 DFT + 学习）、`predict`（只预测） |
| `ERR` | nstep, rmse_energy, rmse_force, rmse_stress | 训练集上的真实 RMSE |
| `BEEF` | nstep, bee_energy, bee_max_force, bee_ave_force, threshold, bee_max_stress, bee_ave_stress | 第 4 列（最大力 BEE）与第 6 列（阈值）用于判断何时触发 DFT |
| `LCONF` | nstep，各元素的局域参考构型数 | 数目随学习增长 |
| `THRUPD` / `THRHIST` | 阈值更新与历史 | — |
| `NORME` | 归一化 RMSE（%，VASP ≥ 6.6） | 误差除以训练数据的标准差，便于跨体系比较 |

健康的训练过程：

- 开始时几乎每步都是 `critical`/`learning`，随后 DFT 调用迅速变稀疏，最终只偶尔学习；
- 如果 50 步之后仍然每步都在学习，或很早就完全停止学习，都说明设置有问题（官方经验）；
- 升温训练中 ERR 缓慢上升是正常的；突然跳升通常意味着晶胞异常形变或发生了相变，要检查 XDATCAR/CONTCAR；
- 贝叶斯误差（BEEF）一般比真实误差（ERR）小，它是"相对"指标。

示例脚本 `parse_ml_logfile.py` 自动读取文件头中的列名，统计各状态步数和 DFT 次数，导出 `ERR.dat`、`BEEF.dat` 等并作图：

```bash
python3 parse_ml_logfile.py -f 1_train/ML_LOGFILE --plot
```

精度参考（官方经验值）：300–1000 K 下能量误差小于约 1 meV/atom、力误差 30–100 meV/Å 通常可接受；液体比简单固体需要更多局域参考构型（约 2000–4000 对比 500–1000）才能达到相近精度。

### 7.3 续训：扩展力场的适用范围

需要覆盖更高温度、其他相、缺陷或表面时，在原有数据基础上续训：

```bash
mkdir 1b_continue && cd 1b_continue
cp ../1_train/ML_ABN ML_AB              # 关键：把上一轮的输出数据库改名为 ML_AB
cp ../INCAR_train INCAR                 # DFT 部分保持不变，只改 MD 设置（温度、步数等）
cp ../new_structure/POSCAR .            # 新结构可以有不同的原子数甚至不同的元素
cp ../KPOINTS ../POTCAR .
mpirun -np 32 vasp_std > vasp.out
```

> [!WARNING]
> 最常见的错误是忘记 `cp ML_ABN ML_AB`：直接复用旧的 ML_AB 会丢掉上一轮新采集的所有数据，续训后的力场反而更差。另外，`ML_MODE = train` 在目录中存在 ML_AB 时**总是**续训；想从零开始就删掉 ML_AB。

`ML_CTIFOR` 在续训时可以设为上一轮最后的阈值，跳过开头不必要的采样；但从液体转到晶体时不要沿用液体的阈值（液体误差大，阈值过高会导致晶体部分完全不学习）。

### 7.4 第 2 步：refit（INCAR_refit）

```bash
# Stage 2: refit the collected data (ML_AB) for the fast prediction mode.
SYSTEM  = Si64 MLFF refit for fast mode
PREC    = Normal
ENCUT   = 350
ISMEAR  = 0
SIGMA   = 0.05
ISYM    = 0
LWAVE   = .FALSE.
LCHARG  = .FALSE.

ML_LMLFF = .TRUE.
ML_MODE  = refit       # sets ML_LFAST = .TRUE., NSW = 1 and an SVD-based fit
```

```bash
mkdir 2_refit && cp INCAR_refit 2_refit/INCAR && cp POSCAR KPOINTS POTCAR 2_refit/
cp 1_train/ML_ABN 2_refit/ML_AB
cd 2_refit && mpirun -np 32 vasp_std > vasp.out
grep ERR ML_LOGFILE | grep -v "#" | tail -1      # refit 后的训练集误差
```

refit 阶段也是调超参数的地方（每次都从默认值出发）：`ML_RCUT1`/`ML_RCUT2`、径向/角向基函数数 `ML_MRB1`/`ML_MRB2`、`ML_LMAX2`、SVD 正则化 `ML_SIGW0`、稀疏化 `ML_EPS_LOW`。官方提示：对液体，`ML_LMAX2 = 2` 与 `ML_RCUT2 = 4` 往往拟合更好；只含 H、B、C、N、O、F 的有机体系，`ML_RCUT2` 取 4 Å 左右可提高学习效率。

### 7.5 第 3 步：生产运行（INCAR_run）

```bash
# Stage 3: production MD with the refitted force field only (no DFT).
IBRION  = 0
MDALGO  = 3
ISIF    = 3
TEBEG   = 800
TEEND   = 800
POTIM   = 2.0
NSW     = 100000       # 200 ps
LANGEVIN_GAMMA   = 1.0 # weaker friction for production
LANGEVIN_GAMMA_L = 10.0
PMASS   = 1000

ML_LMLFF    = .TRUE.
ML_MODE     = run
ML_OUTBLOCK = 10       # write output every 10 steps (default 10 since VASP 6.6)
ML_ESTBLOCK = 50       # spilling factor every 50 steps (VASP <= 6.5.0: ML_IERR)
```

```bash
mkdir 3_run && cp INCAR_run 3_run/INCAR && cp KPOINTS POTCAR 3_run/
cp 1_train/CONTCAR 3_run/POSCAR
cp 2_refit/ML_FFN 3_run/ML_FF                    # 关键：refit 得到的 ML_FFN 改名为 ML_FF
cd 3_run && mpirun -np 32 vasp_std > vasp.out
```

| 关键词 | 说明 |
|---|---|
| `ML_MODE = run` | 只用 ML_FF 预测；仍会初始化 KS 轨道（不使用），部分 DFT 参数（如 ENCUT）会被自动调小 |
| `ML_OUTBLOCK` | 每隔多少步写一次输出（OSZICAR、OUTCAR、XDATCAR、vaspout.h5 等）；`NBLOCK` 会被自动设为 max(NBLOCK, ML_OUTBLOCK)。快速模式每步计算时间与写文件时间相当，必须减小输出频率 |
| `ML_ESTBLOCK` | 每隔多少步计算溢出因子（ML_LOGFILE 中的 `SFF` 行）；只在 `run` 模式中可自由设置，其他模式设为非 1 的值会报错。长 MD 建议 20–100 |
| `ML_OUTPUT_MODE` | `run` 模式默认 0：不计算对关联函数、不写 PCDAT |
| 并行 | MLFF 部分只有 MPI 并行（`run`/`select` 不能用 OpenMP 或 GPU 加速）；需要 scaLAPACK 编译 |

溢出因子接近 0 表示当前环境被训练集很好地覆盖；持续接近 1 表示力场在外推，结果需要谨慎对待。VASP 6.6 起溢出因子始终计算，超过 0.9 时计算会终止——这时应把相应结构所在的条件（温度、压力、相）加入续训。

> [!TIP]
> 力场可以直接用于更大的超胞（计算量与原子数线性相关），例如把 64 原子的 POSCAR 扩成 3×3×3 后做热导、扩散等需要大体系的统计。

### 7.6 验证力场

训练误差小不等于力场可靠，至少还要做：

1. **测试集误差**：从 MLFF 生产轨迹中等间隔抽取若干帧，分别做 DFT 单点（参数与训练**完全相同**）和 MLFF 单点，比较能量、力、应力；py4vasp 提供 `MLFFErrorAnalysis` 辅助比较（见官方 Best practices 页面的示例）。
2. **物理性质对比**：与 DFT 或实验比较晶格常数、声子谱（[第12章](12-phonon.md)）、弹性常数（[第13章](13-elastic.md)）、RDF、相对相能量、缺陷形成能等。
3. **统计收敛**：用不同 `RANDOM_SEED` 跑多条轨迹，时间平均应一致；不一致说明采样不足。

---

## 8. MLFF 关键参数速查

| 标签 | 默认/典型值 | 作用与调节建议 |
|---|---|---|
| `ML_LMLFF` | `.FALSE.` | 总开关；为 `.FALSE.` 时所有 `ML_` 标签被忽略 |
| `ML_MODE` | `train`（ML_LMLFF 为真且未指定时） | 见 §6.3；拼写错误会直接报错退出 |
| `ML_CTIFOR` | 自动更新（初值随 `ML_CALGO` 而定） | 力误差阈值（eV/Å）。采样太少可调低；训练充分后 300–500 K 约 0.02、1000–2000 K 约 0.06 是常见的可达值 |
| `ML_ICRITERIA` | 1（train） | 阈值更新方式：1 用学习后误差的平均（配合 `ML_CX = 0`），2 用全部误差的滑动平均（配合 `ML_CX ≈ 0.1–0.3`，易过采样），0 固定阈值 |
| `ML_CX` | 0 | 阈值 = 历史 BEE 平均 ×(1 + ML_CX)；负值增加采样，定温训练时可设 −0.1 |
| `ML_MB` | 自动（与 NSW、原子数有关） | 每元素局域参考构型上限；液体、无定形、多相体系可能需要增大 |
| `ML_LBASIS_DISCARD` | `.TRUE.` | 达到 ML_MB 时丢弃部分参考构型继续运行（而非停止）；需检查精度 |
| `ML_MCONF` | 自动 | 训练结构数上限；达到时 VASP 停止并提示增大 |
| `ML_MCONF_NEW` | — | 每次学习前缓存的候选结构数；`select` 太慢时可增大到 10–16 |
| `ML_RCUT1` / `ML_RCUT2` | 8.0 / 5.0 Å | 径向 / 角向描述符截断；增大更准但更慢 |
| `ML_EPS_LOW` | — | 稀疏化阈值；增大会删去更多参考构型，不建议超过 1E-7 |
| `ML_WTSIF` | — | 应力权重；表面/分子设为极小值（如 1E-10） |
| `ML_DESC_TYPE` | 0 | 1 为元素约化的三体描述符，多元素体系可显著提速（≥ 6.4.3；所有模式可用需 ≥ 6.5.0） |
| `ML_OUTBLOCK` | 10（6.6） | 生产运行的输出间隔 |
| `ML_ESTBLOCK`（旧名 `ML_IERR`） | 等于 ML_OUTBLOCK | 溢出因子/误差估计的输出间隔 |
| `ML_OUTPUT_MODE` | 0（run） | 关闭对关联函数计算以提速 |

> [!NOTE]
> 同一元素在不同化学环境中（如表面与体相、不同氧化态）可以拆成不同"元素"训练：在 POSCAR 中分组并起不同名字（如 `O1`、`O2`，最多两个字符），POTCAR 中重复对应条目。这能提高精度，但计算量随元素种类数（近似二次）增加。

---

## 9. 外部力场、Δ 学习与其他机器学习势

### 9.1 GRACE 预训练力场（VASP ≥ 6.6，实验性）

VASP 6.6 可以通过 VASPml 库在 `ML_MODE = run` 下使用预训练的 GRACE 力场驱动 MD、结构优化等，需要编译时加入 `-Dlibvaspml` 与 `-DVASPML_ENABLE_GRACE`（以及 TensorFlow、cppflow），模型可用 `tensorpotential` 包的 `grace_models list` / `grace_models download 模型名` 下载：

```bash
# INCAR fragment for a GRACE foundation model (VASP >= 6.6, experimental)
ML_LMLFF       = .TRUE.
ML_LIB         = .TRUE.
ML_MODE        = run
ML_TYPE        = grace
ML_GRACE_MODEL = GRACE-2L-OMAT     # name of a downloaded model
```

注意：

- GRACE 预测会自动使用可用的 GPU，但目前不支持 MPI 并行，必须 `mpirun -np 1` 启动；
- VASP 6.6.0 与 `tensorpotential ≥ 0.5.8` 生成的模型不兼容（known issue 95），6.6.1 已修复；
- 模型目录需要 `metadata.json`（较新的 tensorpotential 会自动提供）；
- 通用预训练模型的精度因体系而异，用于研究前同样要按 §7.6 验证，并注意其训练数据所用的泛函（往往不是你的 PBE 设置）。

### 9.2 Δ 学习：`ML_MODE = delta`（VASP ≥ 6.6.0）

`delta` 模式下，MD 的能量与力 = 当前 INCAR 的从头算结果 + ML_FF 的预测。若 ML_FF 学的是"高精度方法 − 低精度方法"的差值（例如 HSE06 − PBE，或密 k 点 − 稀 k 点），就能以接近低精度方法的代价得到接近高精度方法的轨迹。差值力场的训练数据需要两种理论层次在相同结构上的计算结果，具体组织方式见官方 `ML_MODE` 页面。

### 9.3 与其他机器学习势的关系

| 方案 | 特点 | 与 VASP 的关系 |
|---|---|---|
| VASP MLFF（本章） | 在线主动学习、贝叶斯误差估计、与 DFT 同一程序，上手最快 | 原生 |
| GRACE | 图原子团簇展开，有通用预训练模型 | VASP 6.6 可直接调用（实验性） |
| NEP（GPUMD） | 神经演化势，GPU 上速度极快，适合百万原子、热输运 | 用 VASP 产生训练集（ML_AB / OUTCAR → extxyz） |
| MACE、NequIP/Allegro、DeePMD 等 | 等变/深度神经网络势，精度高，生态丰富 | 同上，在 LAMMPS / ASE 中运行 |

把 VASP 的训练数据迁移到其他框架时，可以用 [utf/pymlff](https://github.com/utf/pymlff) 读写和合并 ML_AB，或用 ASE / pymatgen 从 OUTCAR、vasprun.xml 导出 extxyz。本仓库的 [../../Scripts/README.md](../../Scripts/README.md) 收录了 NEP / GPUMD 训练与数据集处理脚本；机器学习势的整体学习路线见 [Common 第7章](../../Common/tutorials/07-learning-path.md)，高通量工作流见 [Common 第6章](../../Common/tutorials/06-workflow-automation.md)。

---

## 参数选择与收敛建议

| 参数 | AIMD 推荐 | MLFF 训练推荐 | 说明 |
|---|---|---|---|
| `POTIM` | 含 H 0.5；轻元素 1–1.5；重元素 2–3 fs | 同左 | 用 NVE 测试能量守恒 |
| `ENCUT` | 定容 ≈ ENMAX；NpT ≥ 1.3×定容值 | 与生产所需精度一致；NpT 同样提高 | Pulay 应力 |
| `EDIFF` | 1E-5（漂移时 1E-6） | 1E-6 左右 | 训练数据是"标签"，要更准 |
| `NELMIN` / `MAXMIX` | 4 / 40 | 4 / **不设** | MLFF 中 DFT 步之间相隔很远 |
| k 点 | 大超胞 Γ 点（`vasp_gam`） | 足以给出准确的力（半导体小超胞用 2×2×2） | — |
| 恒温器 | NVT：Nosé–Hoover / NHC / CSVR；NpT：Langevin | Langevin（NVT 或 NpT） | 动力学性质用弱耦合恒温器 |
| `LANGEVIN_GAMMA` | 平衡 5–10 ps⁻¹，生产 ≤ 1 ps⁻¹ | 10 ps⁻¹ | 摩擦会减慢扩散 |
| `PMASS` | 默认 1000 起测试 | 同左 | 盒子越大 PMASS 越大 |
| `NBLOCK` / `ML_OUTBLOCK` | 1（需要 VDOS）或 5–10 | 生产运行 10–100 | 控制文件大小 |
| 模拟时长 | 平衡 1–2 ps + 生产 ≥ 10 ps | 训练 1–2 万步以上；生产按需要 | 用块平均评估误差 |
| 训练温度 | — | 目标温度 ×1.3 左右的升温 | 不要远超目标温度 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决办法 |
|---|---|---|
| MD 一开始就崩溃 | `IBRION = 0` 时 `POTIM` 没有默认值 | 在 INCAR 中显式设置 `POTIM` |
| 温度失控上升、守恒能量明显漂移 | 步长太大、EDIFF 太松、每步 SCF 步数不足 | 减小 POTIM，`EDIFF = 1E-6`，`NELMIN = 4–6` |
| `MDALGO` 不起作用或报不认识的标签 | VASP 编译时没有 `-Dtbdyn` | 用官方模板重新编译（5.4.4 后默认包含） |
| MSD 很快饱和、像固体一样，但期望是液体 | 温度不够高或时间太短，结构仍为晶体 | 先高温熔化再降温；检查 g(r) |
| 用 `MDALGO = 0` 的轨迹算出的 MSD 不对 | 旧 Nosé–Hoover 实现把坐标折回盒子 | 改用 `MDALGO = 2/4`，或按 §4.1 展开坐标 |
| NpT 中液体盒子被压扁、严重剪切 | 液体缺乏抗剪切能力 | `ICONST` 约束形状，或 `LATTICE_CONSTRAINTS` |
| NpT 体积持续漂移 | Pulay 应力（ENCUT 不够）、PMASS 不当 | ENCUT +30%，分段重启；测试 PMASS |
| PCDAT / XDATCAR 输出频率和预期不符 | `KBLOCK × NBLOCK > NSW` 时 VASP 会自动重置；MLFF run 中 NBLOCK 至少为 ML_OUTBLOCK | 调整 KBLOCK、NBLOCK |
| 极长 MD 中出现 NaN 或收敛错误 | NSW×NBLOCK 等乘积超过 2³¹−1（known issue 23，6.6.0 修复） | 升级，或减小相关乘积 |
| 6.6.0 中自旋极化 r2SCAN/SCAN 的 MD 受力错误 | known issue 94 | 升级到 6.6.1，或给 6.6.0 打官方 pot.F 补丁 |
| MLFF 续训后力场变差 | 没有把 ML_ABN 复制为 ML_AB，或改动了 DFT 参数/POTCAR | `cp ML_ABN ML_AB`，保持 INCAR 的 DFT 部分完全一致 |
| `ML_MODE = run` 报找不到 ML_FF 或类型不匹配 | 忘记 `cp ML_FFN ML_FF`，或文件损坏 | 复制 refit 得到的 ML_FFN；必要时从最新 ML_AB 重新 refit |
| 训练中 SCF 频繁不收敛、电子步异常多 | 设置了 `MAXMIX > 0` | 删除 MAXMIX |
| 训练很早就停止学习（全是 `accurate`） | `ML_CTIFOR` 太大（如沿用了液体的阈值） | 降低 ML_CTIFOR 或回到默认 |
| 每一步都在学习、几乎没有力场步 | 体系变形/相变、温度过高、DFT 设置不一致 | 检查结构与温度；分阶段训练 |
| 报错要求增大 `ML_MCONF`，或局域参考构型达到 `ML_MB` | 训练结构 / 参考构型数达到上限 | 增大 ML_MCONF / ML_MB；从已有 ML_AB 续训不会丢数据 |
| 设置 `ML_ESTBLOCK`（≠ 1）后在 train 模式报错 | ML_ESTBLOCK 只能在 run 模式自由设置 | 训练时不要设置 |
| 6.6 生产运行因溢出因子 > 0.9 中止 | 力场在严重外推 | 把对应条件加入续训，再 refit |
| 纯 MLFF 计算提示 FFT 网格不足的警告 | 6.5.1 中的误报（known issue 41） | 只做 MLFF 步时可忽略 |
| 6.5.0 多节点 + `-Duse_shmem` 训练结果错误 | known issue 28，6.5.1 修复 | 升级或单节点运行 |
| GRACE 报 `Node 'StatefulPartitionedCall' ... does not have output 1` | 6.6.0 与新版 tensorpotential 模型不兼容（known issue 95） | 升级到 6.6.1 或使用 tensorpotential 0.5.7 的模型 |

---

## 练习

1. 对液态 Si 分别用 `POTIM` = 1、2、4 fs 做 1 ps 的 NVE（`MDALGO = 2` 改为 `MDALGO = 4` + `NHC_NCHAINS = 0`），比较守恒能量的漂移。
2. 用不同 `RANDOM_SEED` 跑三条 1800 K 的液态 Si 轨迹，估计 D 的统计误差。
3. 用 `build_water_box.py --cells 3` 生成 108 个水分子，比较 32 与 108 分子体系的 O–O g(r) 与扩散系数（考虑 Yeh–Hummer 修正）。
4. 按 §7 训练 Si 的 MLFF，再在 1000–1800 K 的液态 Si 上续训；比较只用晶体数据与加入液体数据后，力场在液态 Si 上的测试误差。
5. 用 refit 后的力场计算 Si 的晶格常数、体模量与弹性常数（第13章方法），与 DFT 结果比较。

---

## 参考资料

**官方文档**

- [VASP Wiki: IBRION](https://www.vasp.at/wiki/index.php/IBRION)、[POTIM](https://www.vasp.at/wiki/index.php/POTIM)、[MDALGO](https://www.vasp.at/wiki/index.php/MDALGO)、[SMASS](https://www.vasp.at/wiki/index.php/SMASS)、[TEBEG](https://www.vasp.at/wiki/index.php/TEBEG)、[NBLOCK](https://www.vasp.at/wiki/index.php/NBLOCK)、[RANDOM_SEED](https://www.vasp.at/wiki/index.php/RANDOM_SEED)
- [VASP Wiki: NVT ensemble](https://www.vasp.at/wiki/index.php/NVT_ensemble)、[NpT ensemble](https://www.vasp.at/wiki/index.php/NpT_ensemble)、[Langevin thermostat](https://www.vasp.at/wiki/index.php/Langevin_thermostat)、[LANGEVIN_GAMMA](https://www.vasp.at/wiki/index.php/LANGEVIN_GAMMA)、[LANGEVIN_GAMMA_L](https://www.vasp.at/wiki/index.php/LANGEVIN_GAMMA_L)、[PMASS](https://www.vasp.at/wiki/index.php/PMASS)、[PSTRESS](https://www.vasp.at/wiki/index.php/PSTRESS)
- [VASP Wiki: ICONST](https://www.vasp.at/wiki/index.php/ICONST)、[LATTICE_CONSTRAINTS](https://www.vasp.at/wiki/index.php/LATTICE_CONSTRAINTS)、[MAXMIX](https://www.vasp.at/wiki/index.php/MAXMIX)
- [VASP Wiki: ML_MODE](https://www.vasp.at/wiki/index.php/ML_MODE)、[ML_LMLFF](https://www.vasp.at/wiki/index.php/ML_LMLFF)、[ML_CTIFOR](https://www.vasp.at/wiki/index.php/ML_CTIFOR)、[ML_ICRITERIA](https://www.vasp.at/wiki/index.php/ML_ICRITERIA)、[ML_MB](https://www.vasp.at/wiki/index.php/ML_MB)、[ML_MCONF](https://www.vasp.at/wiki/index.php/ML_MCONF)、[ML_RCUT1](https://www.vasp.at/wiki/index.php/ML_RCUT1)、[ML_RCUT2](https://www.vasp.at/wiki/index.php/ML_RCUT2)、[ML_ESTBLOCK](https://www.vasp.at/wiki/index.php/ML_ESTBLOCK)、[ML_OUTBLOCK](https://www.vasp.at/wiki/index.php/ML_OUTBLOCK)、[ML_LOGFILE](https://www.vasp.at/wiki/index.php/ML_LOGFILE)
- [VASP Wiki: Machine learning force field calculations: Basics](https://www.vasp.at/wiki/index.php/Machine_learning_force_field_calculations:_Basics) — 五步流程（准备 → 训练 → 续训 → refit → 生产）
- [VASP Wiki: Best practices for machine-learned force fields](https://www.vasp.at/wiki/index.php/Best_practices_for_machine-learned_force_fields) — 训练设置、精度经验值、监控与测试
- [VASP Wiki: Machine learning force field: Theory](https://www.vasp.at/wiki/index.php/Machine_learning_force_field:_Theory) — 描述符、贝叶斯回归、溢出因子、稀疏化
- [VASP Wiki: Running GRACE force fields in VASP](https://www.vasp.at/wiki/index.php/Running_GRACE_force_fields_in_VASP)
- [VASP Changelog](https://www.vasp.at/wiki/index.php/Changelog)、[Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — 版本变化与 issue 23、28、41、94、95
- [VASP 官方教程：分子动力学 / 机器学习力场](https://www.vasp.at/tutorials/latest/)
- [py4vasp 文档](https://vasp.at/py4vasp/latest/)

**相关 GitHub 仓库（仅链接）**

- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) — 读取 vaspout.h5 中的能量、结构、对关联函数与 MLFF 误差
- [utf/pymlff](https://github.com/utf/pymlff) — 读写、合并 VASP MLFF 训练集 ML_AB
- [gabkrenzer/vaspMD](https://github.com/gabkrenzer/vaspMD) — AIMD 与 VASP MLFF 分子动力学流程（基准测试、热容）
- [Trebonius91/utils4VASP](https://github.com/Trebonius91/utils4VASP) — 体相/界面建模与分析脚本，涵盖 ML-FF
- [PytLab/VASPy](https://github.com/PytLab/VASPy) — 轻量级 VASP 文件（含 XDATCAR）读写库
- [Ionizing/rsgrad](https://github.com/Ionizing/rsgrad) — 快速查看 OUTCAR 收敛过程与各步受力
- [ishikawa-group/vasp-tutorial](https://github.com/ishikawa-group/vasp-tutorial) — 含 MD 一节的循序渐进教程（日文）
- [Bit-Part-Young/scripts](https://github.com/Bit-Part-Young/scripts) — 含 MLIP、MD 后处理的中文脚本合集

**关键文献**

- R. Jinnouchi, J. Lahnsteiner, F. Karsai, G. Kresse, M. Bokdam, *Phys. Rev. Lett.* 122, 225701 (2019)（贝叶斯在线学习力场）. DOI: [10.1103/PhysRevLett.122.225701](https://doi.org/10.1103/PhysRevLett.122.225701)
- R. Jinnouchi, F. Karsai, G. Kresse, *Phys. Rev. B* 100, 014105 (2019)（在线机器学习力场与熔点）. DOI: [10.1103/PhysRevB.100.014105](https://doi.org/10.1103/PhysRevB.100.014105)
- S. Nosé, *J. Chem. Phys.* 81, 511 (1984). DOI: [10.1063/1.447334](https://doi.org/10.1063/1.447334)；W. G. Hoover, *Phys. Rev. A* 31, 1695 (1985). DOI: [10.1103/PhysRevA.31.1695](https://doi.org/10.1103/PhysRevA.31.1695)
- M. Parrinello, A. Rahman, *Phys. Rev. Lett.* 45, 1196 (1980). DOI: [10.1103/PhysRevLett.45.1196](https://doi.org/10.1103/PhysRevLett.45.1196)；*J. Appl. Phys.* 52, 7182 (1981). DOI: [10.1063/1.328693](https://doi.org/10.1063/1.328693)
- G. Bussi, D. Donadio, M. Parrinello, *J. Chem. Phys.* 126, 014101 (2007)（CSVR 恒温器）. DOI: [10.1063/1.2408420](https://doi.org/10.1063/1.2408420)
- I.-C. Yeh, G. Hummer, *J. Phys. Chem. B* 108, 15873 (2004)（扩散系数的有限尺寸修正）. DOI: [10.1021/jp0477147](https://doi.org/10.1021/jp0477147)
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 132, 154104 (2010)（DFT-D3）. DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)；S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* 32, 1456 (2011)（BJ 阻尼）. DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- A. Bochkarev, Y. Lysogorskiy, R. Drautz, *Phys. Rev. X* 14, 021036 (2024)（GRACE）. DOI: [10.1103/PhysRevX.14.021036](https://doi.org/10.1103/PhysRevX.14.021036)
- Z. Fan et al., *Phys. Rev. B* 104, 104309 (2021)（NEP 神经演化势）. DOI: [10.1103/PhysRevB.104.104309](https://doi.org/10.1103/PhysRevB.104.104309)
- D. Frenkel, B. Smit, *Understanding Molecular Simulation*, Academic Press（块平均、MSD、VACF 等分析方法）

---

[← 上一章](14-optics-dielectric.md) | [返回目录](../README.md) | [下一章 →](16-charge-bonding.md)
