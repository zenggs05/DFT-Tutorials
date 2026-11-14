# 第15章 增强采样：元动力学、PLUMED 与约束 MD（Enhanced Sampling with CP2K: Metadynamics, PLUMED & Constrained MD）

本章讲解怎样在 CP2K 中计算稀有事件的自由能：在 `&FORCE_EVAL/&SUBSYS/&COLVAR` 中定义集体变量（`DISTANCE`、`COORDINATION`、`ANGLE`、`TORSION`、`COMBINE_COLVAR`），用 `&MOTION/&FREE_ENERGY` 的 `METHOD METADYN` 做标准与 well-tempered 元动力学（`DO_HILLS`、`NT_HILLS`、`WW`、`&METAVAR`/`SCALE`、`&WALL`、`WELL_TEMPERED` + `DELTA_T`/`WTGAMMA`），由 `HILLS`/`COLVAR` 输出重构自由能面，通过 `USE_PLUMED` + `PLUMED_INPUT_FILE` 接入 PLUMED，以及用 `&CONSTRAINT/&COLLECTIVE`（`COLVAR`、`TARGET`、`INTERMOLECULAR`）做 blue moon 约束 MD 并对拉格朗日乘子积分。全章以 NaCl 在水中的解离为例，最后给出“GFN1-xTB 预采样 → DFT 精修”的完整流程。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清楚自由能、集体变量（collective variable, CV）和偏置势之间的关系，判断一个问题该用元动力学、约束 MD 还是普通 AIMD；
- 在 CP2K 中定义常用 CV，知道它们的内部单位（bohr、弧度、无量纲）以及 `SCALE`、`TARGET`、墙位置为什么都要用内部单位；
- 写出完整的 well-tempered 元动力学输入，看懂 `-COLVAR.metadynLog` 和 `-HILLS.metadynLog`，重构 1D/2D 自由能面并检查收敛；
- 在编译了 PLUMED 的 CP2K 中用 `plumed.dat` 完成同样的计算，并知道原生实现与 PLUMED 各自的长处；
- 用 blue moon 约束 MD 逐窗口计算平均力，积分得到自由能曲线，并理解 $2k_BT\ln r$ 熵项；
- 用便宜的 xTB 预先摸清 CV 范围、偏置参数和窗口起点，再把计算量集中到 DFT 上。

| 任务 | 推荐方法 | 本章示例 |
|---|---|---|
| 不知道反应路径，想先“逼”体系越过能垒 | well-tempered 元动力学（1–2 个 CV） | `nacl_dissociation_wtmetad.inp` |
| 需要用 PLUMED 的 CV 库（路径 CV、`RMSD`、自定义函数）或 OPES 等新方法 | CP2K + PLUMED | `plumed_metad.inp` + `plumed.dat` |
| 反应坐标明确且是一维的（键长、配位数） | blue moon 约束 MD / 伞形采样 | `constrained_md_distance.inp` |
| 隐藏的慢变量（溶剂化壳层） | 2D 元动力学 | `nacl_wtmetad_2d.inp` |
| DFT 太贵，参数没把握 | xTB 预采样 → DFT 精修 | 第8节 |
| 只有能垒（0 K）、反应物和产物已知 | NEB/Dimer | [第13章](13-neb-dimer.md) |

示例输入集：[../examples/15-metadynamics-nacl/README.md](../examples/15-metadynamics-nacl/README.md)。本章默认你已经掌握 [第12章](12-aimd.md) 的 AIMD 设置（恒温器、ASPC 外推、续算），xTB 的细节见 [第16章](16-qmmm-xtb.md)。

> [!IMPORTANT]
> 增强采样的结果对 CV 的选择极其敏感。CV 选得不好，偏置再强也只会把体系推到不相关的区域。每个新体系都应先跑一段无偏 MD 看 CV 的行为，再决定偏置参数（第3.7节）。

---

## 2. 原理速览

### 2.1 自由能与集体变量

设 $s(\mathbf R)$ 是原子坐标的某个函数（CV），在温度 $T$ 的正则系综中，它的自由能定义为

$$
F(s) = -k_BT\ln P(s) + C,\qquad P(s) = \langle \delta\big(s - s(\mathbf R)\big)\rangle
$$

两个状态之间的自由能差和能垒都可以从 $F(s)$ 上读出。困难在于：能垒比 $k_BT$ 高很多时，普通 MD 在可行的模拟时间（DFT 下通常是几十 ps）内几乎看不到跨越事件，$P(s)$ 在能垒附近得不到采样。增强采样方法都是在 CV 空间里加偏置或加约束，再在统计上把偏置的影响扣除。

### 2.2 标准元动力学

元动力学（metadynamics）每隔一段时间在当前 CV 位置放一个高斯函数，累积成随时间增长的偏置势：

$$
V(\mathbf s,t) = \sum_{j:\,t_j<t} W\prod_{k=1}^{n_\mathrm{CV}} \exp\!\left[-\frac{\big(s_k - s_k(t_j)\big)^2}{2\sigma_k^2}\right]
$$

这正是 CP2K 手册对 `&METAVAR/SCALE` 给出的表达式：`WW` 是高度 $W$（hartree），`SCALE` 是宽度 $\sigma_k$（CV 的内部单位），`NT_HILLS` 是放置间隔（MD 步数）。偏置把已访问过的区域“填平”，体系被迫去探索新的区域。充分填平后

$$
F(\mathbf s) \approx -V(\mathbf s, t\to\infty) + C
$$

缺点是 $V$ 不收敛，而是围绕真实值振荡，误差与 $W$ 和放置频率有关；一直跑下去还会把体系推到物理上无关的高能区域。

### 2.3 Well-tempered 元动力学

Well-tempered（WT）元动力学让高斯高度随已累积的偏置衰减：

$$
W_j = W_0\exp\!\left[-\frac{V(\mathbf s(t_j),t_j)}{k_B\Delta T}\right]
$$

长时间极限下偏置收敛到

$$
V(\mathbf s,t\to\infty) = -\frac{\Delta T}{T+\Delta T}F(\mathbf s) + C,\qquad
F(\mathbf s) = -\frac{T+\Delta T}{\Delta T}\,V(\mathbf s,t\to\infty) + C
$$

CV 的采样分布变成 $P_\mathrm{bias}(\mathbf s)\propto\exp[-F(\mathbf s)/(\gamma k_BT)]$，其中偏置因子

$$
\gamma = \frac{T+\Delta T}{T}
$$

也就是说，能垒在偏置下被“降低”为原来的 $1/\gamma$。CP2K 中 `DELTA_T` 直接给 $\Delta T$，`WTGAMMA` 给 $\gamma$，两者只需写一个；都写时必须满足 $\Delta T = (\gamma-1)T$，否则程序中止（见“常见问题”）。这里的 $T$ 取自 `&MD/TEMPERATURE`；`&METADYN/TEMPERATURE` 是扩展拉格朗日方案用的 CV 温度，与 WT 无关。

> [!NOTE]
> CP2K 在 WT 模式下写入 `HILLS` 文件（以及重启文件）的高度已经乘过 $(T+\Delta T)/\Delta T$（源码 `src/metadynamics.F`），所以无论是否 well-tempered，自由能都等于 HILLS 中所有高斯求和后取负号。PLUMED 的 `HILLS` 文件对 WT 也做了同样的缩放，`plumed sum_hills` 可直接使用。

### 2.4 约束 MD（blue moon）

把 CV 用 SHAKE/RATTLE 严格固定在 $\xi$ 处做 MD，约束力的拉格朗日乘子 $\lambda$ 的系综平均给出自由能梯度。对一般的 CV，需要 Sprik–Ciccotti 的质量度规修正：

$$
\frac{\mathrm dA}{\mathrm d\xi} = \frac{\big\langle Z^{-1/2}\,[-\lambda + k_BT\,G]\big\rangle_{\xi}}{\big\langle Z^{-1/2}\big\rangle_{\xi}},\qquad
Z = \sum_i \frac{1}{m_i}\left|\frac{\partial\xi}{\partial\mathbf r_i}\right|^2
$$

其中 $G$ 是含 $\xi$ 二阶导数的修正项，$\lambda$ 采用 CP2K 输出的符号约定（官方 NaCl 练习中写作 $A(d) = -\int\langle\lambda\rangle\,\mathrm dx$）。**对单个原子间距离**，$Z = 1/m_1 + 1/m_2$ 是常数，$G = 0$，公式简化为

$$
\frac{\mathrm dA}{\mathrm dr} = -\langle\lambda\rangle_r,\qquad A(r) = -\int_{r_0}^{r}\langle\lambda\rangle_{r'}\,\mathrm dr'
$$

这样得到的 $A(r) = -k_BT\ln P(r)$ 包含了球壳体积带来的 Jacobian 项。若要与经典模拟中常用的平均力势 $w(r) = -k_BT\ln g(r)$ 对比，需要加上

$$
w(r) = A(r) + 2k_BT\ln r + C
$$

示例脚本 `bluemoon_integrate.py` 的 `--pmf` 选项就是做这一步。

### 2.5 方法对比

| | 标准 MetaD | WT-MetaD | Blue moon | 伞形采样 |
|---|---|---|---|---|
| 需要预先知道路径 | 否 | 否 | 是（一维坐标） | 是 |
| 收敛性 | 振荡，不收敛 | 平滑收敛 | 每个窗口独立收敛 | 每个窗口独立收敛 |
| CV 维数 | 1–3 | 1–3 | 通常 1 | 通常 1 |
| 并行方式 | 单条轨迹（可多 walker） | 同左 | 窗口完全独立 | 窗口完全独立 |
| 后处理 | 高斯求和 | 高斯求和 / 重加权 | 平均力积分 | WHAM/MBAR |
| CP2K 关键字 | `&METADYN` | `WELL_TEMPERED` | `&COLLECTIVE` | `&COLLECTIVE/&RESTRAINT` |

---

## 3. CP2K 中的集体变量 `&COLVAR`

### 3.1 位置、编号和单位

CV 定义在 `&FORCE_EVAL/&SUBSYS` 中，每个 `&COLVAR` 段定义一个 CV。其他地方（`&METAVAR/COLVAR`、`&COLLECTIVE/COLVAR`）按**输入文件中出现的顺序**用 1, 2, 3… 引用。原子序号从 1 开始，与坐标文件中的顺序一致。

| CV | 内部单位 | 典型 `SCALE` |
|---|---|---|
| `DISTANCE` | bohr | 0.1–0.4 bohr |
| `ANGLE`、`TORSION` | rad | 0.1–0.3 rad |
| `COORDINATION` | 无量纲 | 0.05–0.3 |
| `COMBINE_COLVAR` | 由函数决定 | 视情况而定 |

> [!WARNING]
> `&METAVAR/SCALE`、`&WALL/POSITION` 和 `&WALL/&QUADRATIC/K` 在手册中的单位标为 `internal_cp2k` 或 `hartree`，即距离按 bohr、角度按弧度计算。把 0.2 Å 写成 `SCALE 0.2`，实际宽度只有 0.106 Å。`&COLLECTIVE/TARGET` 可以写单位前缀，如 `TARGET [angstrom] 3.20`；不写时同样按内部单位。

### 3.2 `DISTANCE`

```fortran
&COLVAR
  &DISTANCE
    ATOMS 1 2          ! Na and Cl
  &END DISTANCE
&END COLVAR
```

在周期性体系中取最小镜像距离。`AXIS Z` 只计算 z 分量（例如吸附物离表面的高度），配合 `SIGN T` 保留正负号。`ATOMS` 也可以引用 `&POINT` 子段定义的几何中心。

### 3.3 `COORDINATION`

配位数用一个光滑的开关函数统计“成键”原子对：

$$
\mathrm{CN} = \frac{1}{N_\mathrm{from}}\sum_{i\in\mathrm{from}}\sum_{j\in\mathrm{to}}\frac{1-(r_{ij}/R_0)^{NN}}{1-(r_{ij}/R_0)^{ND}}
$$

```fortran
&COLVAR
  &COORDINATION
    ATOMS_FROM 1           ! Na
    KINDS_TO   O           ! all water oxygens
    R0 [angstrom] 3.0
    NN 12
    ND 24
  &END COORDINATION
&END COLVAR
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ATOMS_FROM` / `KINDS_FROM` | 中心原子 | 两者二选一，`KINDS_*` 按元素（kind）选 |
| `ATOMS_TO` / `KINDS_TO` | 配体原子 | 同上 |
| `R0` | 第一配位壳层与第二壳层之间（RDF 第一极小附近） | 默认 3.0 **bohr**，建议写 `[angstrom]` |
| `NN` | 6–12 | 默认 6，越大开关越陡 |
| `ND` | 2×`NN` | 默认 12；`ND` > `NN` 才能保证远处衰减到 0 |

开关函数太“软”（`NN 6`，`ND 12`）会把第二壳层也算进来一部分；太“硬”则导数集中在 $R_0$ 附近一个很窄的区间里，偏置力会变得很尖锐。拿到无偏轨迹后，可以对 Na–O 距离直接算这个函数，看第一壳层原子的贡献是否接近 1、第二壳层是否接近 0。`ATOMS_TO_B`/`KINDS_TO_B` 与 `R0_B`、`NN_B`、`ND_B` 用于两个配位数相乘的 CV。

### 3.4 `ANGLE` 与 `TORSION`

```fortran
&COLVAR
  &ANGLE
    ATOMS 3 1 4            ! angle 3-1-4, atom 1 at the apex
  &END ANGLE
&END COLVAR
&COLVAR
  &TORSION
    ATOMS 5 7 9 15         ! dihedral 5-7-9-15
  &END TORSION
&END COLVAR
```

两者都以弧度计算。`TORSION` 是周期性 CV，CP2K 在计算高斯与当前 CV 的差值时会自动折回到 $[-\pi,\pi]$；`fes_from_hills.py` 需要用 `--periodic 6.2831853` 告诉它周期。（这段原子序号只是示意，NaCl 示例中没有用到。）

### 3.5 `COMBINE_COLVAR`

把已有的 CV 用解析函数组合成一个新 CV，例如质子转移坐标 $\delta = d(\mathrm{O_1\!-\!H}) - d(\mathrm{O_2\!-\!H})$：

```fortran
&COLVAR
  &COMBINE_COLVAR
    &COLVAR
      &DISTANCE
        ATOMS 1 3          ! O1-H
      &END DISTANCE
    &END COLVAR
    &COLVAR
      &DISTANCE
        ATOMS 2 3          ! O2-H
      &END DISTANCE
    &END COLVAR
    FUNCTION D1-D2
    VARIABLES D1 D2
  &END COMBINE_COLVAR
&END COLVAR
```

| 关键词 | 说明 |
|---|---|
| 内嵌 `&COLVAR` | 参与组合的 CV，按出现顺序对应 `VARIABLES` |
| `FUNCTION` | 支持 `+ - * / ^`（或 `**`）和 Fortran 内置函数（`exp`、`log`、`sqrt`、`tanh` 等） |
| `VARIABLES` | 变量名，与内嵌 `&COLVAR` 的顺序一一对应 |
| `PARAMETERS` / `VALUES` | 函数中的常数参数名及数值 |
| `DX`、`ERROR_LIMIT` | 用 Ridders 方法数值求导的步长（默认 0.1 bohr）和误差上限；误差超限时打印警告 |

> [!TIP]
> 组合 CV 的导数靠数值差分。函数里如果出现 `sqrt`、除法之类的奇点，或量纲混杂（距离 + 配位数），偏置力会出问题。能用一个内置 CV 表达的，就不要用 `COMBINE_COLVAR`。

### 3.6 其他内置 CV

`&COLVAR` 还提供 `ANGLE_PLANE_PLANE`、`DISTANCE_POINT_PLANE`、`GYRATION_RADIUS`、`RMSD`、`REACTION_PATH`/`DISTANCE_FROM_PATH`（路径 CV）、`HYDRONIUM_SHELL`/`HYDRONIUM_DISTANCE`/`ACID_HYDRONIUM_*`（质子转移）、`POPULATION`、`QPARM`、`RING_PUCKERING`、`HBP`、`WC`、`U`（势能本身）等，完整列表与关键字见 [`&COLVAR` 输入参考](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html)。PLUMED 的 CV 库更大，需要时见第6节。

### 3.7 先跑无偏 MD，看 CV 怎么涨落

选 `SCALE` 的依据是 CV 在初始盆地中自发的涨落幅度。CP2K 的做法是照常打开 `&METADYN`，但令 `DO_HILLS F`：不放任何高斯，只在每步计算并输出 CV。示例 `cv_monitor.inp` 的核心部分：

```fortran
&FREE_ENERGY
  METHOD METADYN
  &METADYN
    DO_HILLS F
    &METAVAR
      COLVAR 1
      SCALE 0.35
    &END METAVAR
    &METAVAR
      COLVAR 2
      SCALE 0.20
    &END METAVAR
    &PRINT
      &COLVAR
        COMMON_ITERATION_LEVELS 3
        &EACH
          MD 1
        &END EACH
      &END COLVAR
    &END PRINT
  &END METADYN
&END FREE_ENERGY
```

```bash
cd CP2K/examples/15-metadynamics-nacl
mpirun -np 16 cp2k.psmp -i cv_monitor.inp -o cv_monitor_xtb.out
# columns 2 and 3: CV1 (bohr) and CV2; mean and standard deviation after 1 ps (2000 steps)
awk 'NR>2000{n++; s+=$2; ss+=$2*$2; c+=$3; cc+=$3*$3}
     END{printf "d: %.3f +- %.3f bohr   CN: %.3f +- %.3f\n", s/n, sqrt(ss/n-(s/n)^2), c/n, sqrt(cc/n-(c/n)^2)}' \
    nacl_cvmon_xtb-COLVAR.metadynLog
```

经验法则：`SCALE` 取无偏涨落标准差的 1/3–1/2。宽度太大，自由能面上的细节（窄极小、肩峰）会被抹平；太小则需要更多高斯才能填平同样的体积。

如果已经有一条 AIMD 轨迹，也可以用 `ENSEMBLE REFTRAJ` 重新读入轨迹，只计算 CV 而不算能量（官方元动力学教程就是这样做的）。`&REFTRAJ/EVAL` 自 2024.1 起可用，默认值就是 `NONE`；2023.x 及更早的版本只有逻辑关键字 `EVAL_ENERGY_FORCES`（写 `.FALSE.`）。2024.1–2025.2 两种写法都接受，2026.2 删除了旧关键字，含旧关键字的输入在解析阶段就会中止：

```fortran
&MD
  ENSEMBLE REFTRAJ
  STEPS 1000                       ! number of frames to process
  &REFTRAJ
    TRAJ_FILE_NAME nacl_equil_xtb-pos-1.xyz
    FIRST_SNAPSHOT 1
    LAST_SNAPSHOT 1000
    STRIDE 1
    EVAL NONE                      ! CP2K >= 2024.1 (only form accepted by 2026.2)
    ! EVAL_ENERGY_FORCES .FALSE.   ! CP2K <= 2023.x (removed in 2026.2)
  &END REFTRAJ
&END MD
```

REFTRAJ 的更多用法见 [第12章](12-aimd.md)。

---

## 4. 示例体系与准备：NaCl 水溶液

示例体系是 1 个 Na⁺、1 个 Cl⁻ 和 62 个水分子，放在 12.42 Å 的立方盒子中（188 个原子，约 1 g/cm³）。离子对在水中解离是增强采样的经典教学问题：CV 简单（Na–Cl 距离），但溶剂化壳层的重排是隐藏的慢自由度。CP2K 官方练习也用这个体系做约束 MD（见参考资料）。

所有输入都把 `&FORCE_EVAL` 放在 include 文件里，用预处理变量切换理论层次：

```fortran
@SET LEVEL ${LEVEL-dft}
@SET COORD ${COORD-nacl_eq.xyz}
...
@INCLUDE 'nacl_${LEVEL}_force_eval.inc'
```

命令行 `-E LEVEL=xtb` 选 `nacl_xtb_force_eval.inc`（GFN1-xTB），`-E LEVEL=dft` 选 `nacl_dft_force_eval.inc`（PBE-D3(BJ)）。两个文件中的 `&COLVAR` 完全一致：

```fortran
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &QS
      METHOD XTB
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
      &XTB
        GFN_TYPE 1
        DO_EWALD T
      &END XTB
    &END QS
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.42 12.42 12.42
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    ! COLVAR 1: Na-Cl distance (internal unit: bohr)
    &COLVAR
      &DISTANCE
        ATOMS 1 2
      &END DISTANCE
    &END COLVAR
    ! COLVAR 2: coordination number of Na by water oxygens
    &COLVAR
      &COORDINATION
        ATOMS_FROM 1
        KINDS_TO O
        R0 [angstrom] 3.0
        NN 12
        ND 24
      &END COORDINATION
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL
```

（上面省略了 SCF 重启文件的打印设置，完整文件见示例目录。）

| 设置 | xTB 层次 | DFT 层次 |
|---|---|---|
| 方法 | GFN1-xTB（CP2K 原生，`GFN_TYPE 1`） | PBE + D3(BJ) |
| 基组/赝势 | — | DZVP-MOLOPT-PBE-GTH（`BASIS_MOLOPT_UZH`）+ GTH-PBE（`POTENTIAL_UZH`），Na q9、Cl q7、O q6、H q1 |
| 网格 | — | `CUTOFF 400`、`REL_CUTOFF 50` |
| SCF | OT + `OUTER_SCF`，`EPS_SCF 1E-6` | 同左，`EPS_DEFAULT 1E-12` |
| 外推 | ASPC，阶数 3 | ASPC，阶数 3 |
| 用途 | 平衡、CV 涨落、预采样、挑选窗口起点 | 生产计算 |

> [!NOTE]
> `DO_EWALD T` 在 2026.1 及更早版本中是周期性 xTB 必需的；2026.2 手册把它标为已弃用（周期性改由 `&CELL/PERIODIC` 推断），但仍然接受。CP2K 原生 xTB 的 `CHECK_ATOMIC_CHARGES` 默认为 `T`：偏置把体系推到极端构型、原子电荷超出化学上合理的范围时，计算会主动中止，这时应先检查结构，而不是直接关掉检查。GFN1-xTB 对液态水只是定性描述，第8节会讨论这对预采样意味着什么。

### 4.1 平衡：带谐振约束的 NVT

格点生成的起始结构离平衡很远，先用 xTB 跑 5 ps NVT。为了防止离子对在平衡阶段就自发分开，用一个谐振**约束势**（restraint，不是严格约束）把 $d$ 拉在 2.80 Å 附近：

```fortran
&CONSTRAINT
  &COLLECTIVE
    COLVAR 1
    INTERMOLECULAR T
    TARGET [angstrom] 2.80
    &RESTRAINT
      K 0.05
    &END RESTRAINT
  &END COLLECTIVE
&END CONSTRAINT
```

`&RESTRAINT` 段一出现，`&COLLECTIVE` 就从 SHAKE 约束变成谐振势 $K(\xi-\xi_0)^2$（注意没有 1/2）。`K` 用内部单位，0.05 Ha/bohr² ≈ 469 kJ mol⁻¹ Å⁻²。

```bash
mpirun -np 16 cp2k.psmp -i nacl_equil.inp -o nacl_equil_xtb.out
grep "PROGRAM ENDED AT" nacl_equil_xtb.out
bash last_frame.sh nacl_equil_xtb-pos-1.xyz nacl_eq.xyz
```

检查 `nacl_equil_xtb-1.ener` 中温度在 300 K 附近波动、势能不再单调下降（各列含义见 [第12章](12-aimd.md)）。

---

## 5. 实战一：well-tempered 元动力学（CP2K 原生）

### 5.1 完整输入

`nacl_dissociation_wtmetad.inp`：

```fortran
@SET LEVEL ${LEVEL-dft}
@SET COORD ${COORD-nacl_eq.xyz}

&GLOBAL
  PROJECT nacl_wtmetad_${LEVEL}
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL

! &EXT_RESTART
!   RESTART_FILE_NAME nacl_wtmetad_dft-1.restart
! &END EXT_RESTART

@INCLUDE 'nacl_${LEVEL}_force_eval.inc'

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS ${STEPS-40000}
    TIMESTEP 0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE CSVR
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
  &FREE_ENERGY
    METHOD METADYN
    &METADYN
      DO_HILLS T
      NT_HILLS 40
      WW 4.0E-4
      WELL_TEMPERED T
      DELTA_T 1500.0
      &METAVAR
        COLVAR 1
        SCALE 0.35
        &WALL
          TYPE QUADRATIC
          POSITION 10.6
          &QUADRATIC
            DIRECTION WALL_PLUS
            K 0.05
          &END QUADRATIC
        &END WALL
      &END METAVAR
      &PRINT
        &COLVAR
          COMMON_ITERATION_LEVELS 3
          &EACH
            MD 1
          &END EACH
        &END COLVAR
        &HILLS
          COMMON_ITERATION_LEVELS 3
          &EACH
            MD 1
          &END EACH
        &END HILLS
      &END PRINT
    &END METADYN
  &END FREE_ENERGY
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
      &EACH
        MD 100
      &END EACH
    &END RESTART
    &RESTART_HISTORY OFF
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

### 5.2 关键词解释

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `&FREE_ENERGY/METHOD` | `METADYN` | 默认值就是 `METADYN`；另有 `UI`（伞形积分）和 `AC`（炼金术变换） |
| `DO_HILLS` | `T` | 默认 `F`，即只监测 CV、不加偏置 |
| `NT_HILLS` | 20–100 步 | 放置间隔（默认 30）。应长于 CV 的弛豫时间，否则偏置“追着”体系跑 |
| `WW` | 0.5–2 kJ/mol（2E-4 – 8E-4 Ha） | 高斯高度，单位 hartree。**默认 0.1 Ha（约 263 kJ/mol），必须改** |
| `WELL_TEMPERED` | `T` | 打开 WT，此时必须给 `DELTA_T` 或 `WTGAMMA` |
| `DELTA_T` | $(\gamma-1)T$ | 单位 K。示例 1500 K，即 $\gamma = 6$ |
| `WTGAMMA` | 5–15 | 与 `DELTA_T` 二选一；`nacl_wtmetad_2d.inp` 用 `WTGAMMA 8` |
| `&METAVAR/COLVAR` | CV 编号 | 对应 `&SUBSYS` 中第几个 `&COLVAR` |
| `&METAVAR/SCALE` | 涨落的 1/3–1/2 | 高斯宽度，内部单位（别名 `WIDTH`）。0.35 bohr = 0.185 Å |
| `&WALL/TYPE` | `QUADRATIC` | 另有 `REFLECTIVE`（反转 CV 速度）、`QUARTIC`、`GAUSSIAN` |
| `&WALL/POSITION` | 盒长一半以内 | 内部单位。10.6 bohr = 5.61 Å < 6.21 Å |
| `&QUADRATIC/DIRECTION` | `WALL_PLUS` | `WALL_PLUS` 作用于 CV > POSITION 一侧，`WALL_MINUS` 作用于 CV < POSITION 一侧 |
| `&QUADRATIC/K` | 0.01–0.1 | 势为 $K(\mathrm{CV}-\mathrm{POS})^2$ |
| `&PRINT/&COLVAR`、`&HILLS` | 每步 | 写出 `<PROJECT>-COLVAR.metadynLog` 和 `<PROJECT>-HILLS.metadynLog`；`COMMON_ITERATION_LEVELS 3` 让内容写进同一个文件 |
| `&THERMOSTAT` | CSVR，`TIMECON` 100 fs | 偏置不断向体系注入能量，必须用 NVT |

为什么需要墙：$d$ 超过盒长一半后，最小镜像会让离子“看到”另一个周期像，CV 不再有物理意义；而且 WT-MetaD 会不停地把体系推向从未去过的大 $d$ 区域。墙势不计入 HILLS，所以重构出的 $F(s)$ 在墙外不可信，分析时只看墙内区间。

参数之间的关系可以这样估算：每 20 fs 放一个 1.05 kJ/mol 的高斯，CV 方向宽 0.185 Å，要填平的区间约 3 Å。WT 下高斯会逐渐变矮，最终偏置约为 $(\gamma-1)/\gamma \approx 83\%$ 的 $-F$。

### 5.3 运行

```bash
# xTB 预采样（100 ps），几个小时级别
mpirun -np 16 cp2k.psmp -E LEVEL=xtb -E STEPS=200000 \
       -i nacl_dissociation_wtmetad.inp -o wtmetad_xtb.out
# DFT（默认 20 ps），需要较多计算资源
mpirun -np 64 cp2k.psmp -i nacl_dissociation_wtmetad.inp -o wtmetad_dft.out
```

Slurm 下用 `sbatch --export=ALL,STEP=wtmetad,LEVEL=xtb,STEPS=200000 job_slurm.sh`（`STEPS` 可省略，省略时用输入文件中的默认步数）。作业脚本的写法见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)，CP2K 部分的通用模板是 [../scripts/job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh)（用法写在文件头的注释里）。

### 5.4 输出文件

`<PROJECT>-COLVAR.metadynLog`：每步一行。没有扩展拉格朗日时，各列依次为：

| 列 | 内容 |
|---|---|
| 1 | 时间（fs） |
| 2 … $n+1$ | 各 CV 的瞬时值（内部单位） |
| 接下来 $n$ 列 | 高斯偏置作用在各 CV 上的“力” $-\partial V_\mathrm{hills}/\partial s_k$（官方文档称为梯度；源码中的 `ff_hills`） |
| 接下来 $n$ 列 | 墙势作用在各 CV 上的“力”（`ff_walls`） |
| 倒数第 2 列 | 高斯偏置势的值（hartree） |
| 最后一列 | 墙势的值（hartree） |

`<PROJECT>-HILLS.metadynLog`：每个高斯一行：时间、各 CV 的中心、各 CV 的宽度、高度（hartree，WT 下已乘 $(T+\Delta T)/\Delta T$）。

> [!WARNING]
> CP2K 用 Fortran 格式 `f13.5` 写 HILLS（源码 `src/metadynamics.F`），高度以 hartree 为单位只保留 5 位小数：$4.8\times10^{-4}$ Ha 写成 `0.00048`，WT 后期约 $3\times10^{-5}$ Ha 的高斯相对误差可达 15%。HILLS 适合监控和快速作图；定量重构请改用重启文件 `<PROJECT>-1.restart` 中的 `&SPAWNED_HILLS_POS`/`_SCALE`/`_HEIGHT` 段，那里保存的是全精度数值（`graph` 工具和 `fes_from_hills.py` 都能直接读它，见第5.5节）。

```bash
P=nacl_wtmetad_xtb
tail -3 $P-COLVAR.metadynLog
wc -l $P-HILLS.metadynLog                     # number of hills ~ STEPS / NT_HILLS
awk '{printf "%10.1f %8.3f %12.4e\n", $1, $2*0.529177, $NF}' $P-HILLS.metadynLog | sed -n '1p;100p;1000p;$p'
awk '{print $1/1000, $2*0.529177}' $P-COLVAR.metadynLog > d_vs_t.dat   # ps vs Angstrom
```

示意输出（数值仅供参考）：

```text
      20.0    2.813   4.8000e-04
    2000.0    3.402   3.2000e-04
   20000.0    4.876   1.0000e-04
  100000.0    3.127   4.0000e-05
```

（第3列是 HILLS 中缩放后的高度，由于上面说的 5 位小数，只有 1–2 位有效数字。第一个高斯放在偏置为零处，高度为 $W_0\cdot\gamma/(\gamma-1) = 4.0\times10^{-4}\times 1.2 = 4.8\times10^{-4}$ Ha；其余数字是虚构的，只用来说明高度应随时间下降。）

### 5.5 重构自由能面

示例脚本 `fes_from_hills.py` 对 HILLS 中的高斯求和：

```bash
python3 fes_from_hills.py nacl_wtmetad_xtb-HILLS.metadynLog --bohr2ang --blocks 4 -o fes_xtb.dat --png fes_xtb.png
```

- `--bohr2ang`：把 CV 和宽度从 bohr 换成 Å；
- `--blocks 4`：分别用前 1/4、2/4、3/4 和全部高斯计算 $F(s)$，写成 `fes_xtb.dat` 的不同列；
- 屏幕上会打印高斯数、第一个和最后一个高斯的高度，以及曲线上的局部极小和极大。

把第一个参数换成重启文件，脚本会自动识别并读取全精度的 `&SPAWNED_HILLS_*` 段（此时没有时间信息，`--blocks` 按高斯序号分段）：

```bash
python3 fes_from_hills.py nacl_wtmetad_xtb-1.restart --bohr2ang --blocks 4 -o fes_xtb_restart.dat
```

CP2K 自带的 `graph` 程序（源码 `src/metadyn_tools/graph.F`，可执行文件名为 `graph.psmp` 或 `graph.ssmp`，取决于构建）也能从重启文件重构自由能面，对 2–3 个 CV 尤其方便：

```bash
graph.psmp -cp2k -file nacl_wtmetad_xtb-1.restart -ndim 1 -out fes_graph.dat
```

`-cp2k` 表示输入是 CP2K 重启文件（高斯保存在其中的 `SPAWNED_HILLS_*` 段），`-ndim` 是 CV 数，`-stride N` 每 N 个高斯输出一次自由能，可以用来看收敛过程。默认输出文件是 `fes.dat`。只输出部分 CV 时用 `-ndw`（例如 `-ndim 2 -ndw 1`）；此时 `graph` 默认取另一维上的**最小值**作为投影，要做 Boltzmann 积分需加 `-integrated_fes`（不带参数运行 `graph.psmp` 会打印全部选项）。

**收敛判据**（缺一不可）：

1. CV 在各个盆地之间来回了多次（`d_vs_t.dat` 画出来像“方波”），而不是只跨越了一次；
2. WT 下高斯高度已经降到初始值的 10–20% 以下；
3. `--blocks` 给出的后两条曲线在关心的区间内（墙内）差别在 $k_BT$ 量级以内；
4. 两个极小之间的自由能差 $\Delta F$ 随时间不再漂移。可以把 $\Delta F(t)$ 画出来，它应在某个值附近趋于平稳。

> [!CAUTION]
> 在 DFT 层次上，20–50 ps 的 1D WT-MetaD 常常只能得到“定性正确”的曲线。如果 CV 只跨越过一两次，**不要**报告能垒的具体数值。把收敛性检查的图放进论文的补充材料，是增强采样研究的基本要求（见 PLUMED NEST 的倡议，参考资料）。

### 5.6 续算

元动力学的高斯历史保存在 `.restart` 文件的 `&METADYN` 段中，续算时：

1. 取消输入中 `&EXT_RESTART` 的注释，`RESTART_FILE_NAME` 指向最新的 `<PROJECT>-1.restart`；
2. `STEPS` 是**本次还要走的**步数（源码 `md_run.F` 中的循环是 `DO istep = 1, simpar%nsteps`），步号计数器则从重启文件中恢复。例如第一段走了 40000 步，续算时写 `STEPS 40000` 就会走到第 80000 步；需要给整条轨迹设硬上限时，在 `&MD` 中加 `MAX_STEPS 80000`（与 [第12章](12-aimd.md) 的续算约定相同）；
3. CP2K 默认以追加方式打开输出文件，`HILLS` 和 `COLVAR` 会接在原文件后面。运行后检查时间列是否连续，必要时去掉重复行。

也可以在续算时修改 `WW`、`NT_HILLS`、`SCALE`。官方教程指出，这在某个事件发生后 CV 的动力学发生变化时很有用。续算流程的通用细节见 [第12章](12-aimd.md)；辅助脚本 [../scripts/restart_update.py](../scripts/restart_update.py) 的用法用 `python ../scripts/restart_update.py -h` 查看。

### 5.7 二维元动力学：距离 + 配位数

只偏置 $d$ 时，Na⁺ 的水合壳层重排可能跟不上离子对的分离，自由能曲线会出现滞后（正向、反向跨越时的形状不一样）。把 Na–O 配位数（COLVAR 2）也加入偏置，就是 `nacl_wtmetad_2d.inp`：

```fortran
&METADYN
  DO_HILLS T
  NT_HILLS 40
  WW 4.0E-4
  WELL_TEMPERED T
  WTGAMMA 8.0
  &METAVAR
    COLVAR 1
    SCALE 0.35
    &WALL
      TYPE QUADRATIC
      POSITION 10.6
      &QUADRATIC
        DIRECTION WALL_PLUS
        K 0.05
      &END QUADRATIC
    &END WALL
  &END METAVAR
  &METAVAR
    COLVAR 2
    SCALE 0.20
    &WALL
      TYPE QUADRATIC
      POSITION 3.5
      &QUADRATIC
        DIRECTION WALL_MINUS
        K 0.01
      &END QUADRATIC
    &END WALL
  &END METAVAR
  ...
&END METADYN
```

2D 需要填平的“体积”大得多，高斯数通常要多一个数量级，所以示例默认在 xTB 上跑 200 ps。重构时：

```bash
python3 fes_from_hills.py nacl_wt2d_xtb-HILLS.metadynLog --bohr2ang 1 -T 300 -o fes_d.dat --out2d fes2d.dat
graph.psmp -cp2k -file nacl_wt2d_xtb-1.restart -ndim 2 -out fes2d_graph.dat
```

`--bohr2ang 1` 只转换第 1 个 CV（距离），配位数保持无量纲；`-T 300` 用于把 $F(d,\mathrm{CN})$ 投影到一维：$F(d) = -k_BT\ln\int e^{-F(d,\mathrm{CN})/k_BT}\,\mathrm d\mathrm{CN}$。`fes2d.dat` 是 gnuplot `splot` 格式。

> [!TIP]
> 一维投影**不能**简单地把 2D 高斯只按第一个 CV 求和，必须按上式做 Boltzmann 积分。另外，`SCALE 0` 在 CP2K 中表示该方向上无限宽，这可以用来在多维元动力学中组合几个一维偏置。

### 5.8 其他变体（简介）

- **标准元动力学**：去掉 `WELL_TEMPERED`/`DELTA_T` 即可。只建议用于快速探索，不用于定量计算；
- **扩展拉格朗日元动力学**（`LAGRANGE T`，`&METAVAR` 中的 `MASS`、`LAMBDA`，Iannuzzi–Laio–Parrinello 2003）：引入与 CV 耦合的辅助粒子。它**不能**与 WT 同时使用；
- **自适应放置**（`MIN_DISP` > 0，可配合 `MIN_NT_HILLS`）：设了 `MIN_DISP` 就不能再用 WT（程序中止，见“常见问题”第3条）；
- **多 walker**（`&MULTIPLE_WALKERS`）：多条轨迹共享同一组高斯；
- **高斯尾部截断**（`HILL_TAIL_CUTOFF`、`P_EXPONENT`、`Q_EXPONENT`）：让高斯远端更快衰减到零。

这些选项的完整说明见 [`&METADYN` 输入参考](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/FREE_ENERGY/METADYN.html)。

---

## 6. 实战二：用 PLUMED 驱动元动力学

### 6.1 什么时候用 PLUMED

| | CP2K 原生 `&METADYN` | PLUMED |
|---|---|---|
| 是否需要额外编译 | 否 | 是（`-DCP2K_USE_PLUMED=ON` 或 toolchain `--with-plumed`） |
| CV 种类 | `&COLVAR` 中的几十种 | 数百种，支持自定义函数、路径 CV、机器学习 CV（另需模块） |
| 偏置方法 | MetaD、WT-MetaD、多 walker、约束/约束势 | MetaD、WT、OPES、伞形采样、移动约束势、VES 等（OPES、VES 需在编译 PLUMED 时启用相应模块） |
| 输入可移植性 | 只能在 CP2K 中使用 | 同一 `plumed.dat` 可用于 GROMACS、LAMMPS 等 |
| 后处理 | 自己写脚本 / `graph` | `plumed sum_hills`、`plumed driver`、重加权工具 |

当 CP2K 原生 CV 不够用，或者希望用经典 MD 程序先在力场上调好 `plumed.dat` 再拿到 AIMD 上用时，就选 PLUMED。

### 6.2 检查 CP2K 是否带 PLUMED

```bash
cp2k.psmp -v | grep -o plumed2       # prints "plumed2" if PLUMED support is compiled in
```

- 源码构建：CMake 选项 `-DCP2K_USE_PLUMED=ON`（它依赖 `CP2K_USE_MPI`，只有 MPI 构建才能打开）；toolchain 用 `--with-plumed=install`（或指向系统已有的 PLUMED）。toolchain 同样只在 MPI 构建中启用 PLUMED，不用 MPI 时会给出警告 "Not using MPI, so PLUMED is disabled."；
- 发行版（conda-forge、Docker 镜像等）是否带 PLUMED 因版本而异，以 `-v` 的输出为准；
- 运行时还需要能找到 PLUMED 动态库（`LD_LIBRARY_PATH` 或 `source` PLUMED 的环境脚本）。

安装与构建细节见 [第1章](01-install.md)。

### 6.3 CP2K 输入

`plumed_metad.inp` 的 `&FREE_ENERGY` 段只有两行：

```fortran
&FREE_ENERGY
  METHOD METADYN
  &METADYN
    USE_PLUMED T
    PLUMED_INPUT_FILE plumed.dat
  &END METADYN
&END FREE_ENERGY
```

`USE_PLUMED T` 时，CP2K 只读取 `PLUMED_INPUT_FILE`（默认 `./plumed.dat`），`&METADYN` 中的其他关键字都会被忽略。include 文件里的 `&COLVAR` 也不再起作用（保留它们没有坏处）。其余部分（`&MD`、恒温器、输出）与原生版本完全相同。

### 6.4 `plumed.dat`

```text
UNITS LENGTH=A ENERGY=kj/mol TIME=fs

d:  DISTANCE ATOMS=1,2
cn: COORDINATION GROUPA=1 GROUPB=3-186:3 R_0=3.0 NN=12 MM=24

metad: METAD ARG=d SIGMA=0.185 HEIGHT=1.05 PACE=40 BIASFACTOR=6 TEMP=300 FILE=HILLS GRID_MIN=1.0 GRID_MAX=8.0 GRID_BIN=700

uwall: UPPER_WALLS ARG=d AT=5.6 KAPPA=500.0

PRINT ARG=d,cn,metad.bias,uwall.bias STRIDE=20 FILE=COLVAR
```

| 行 | 说明 |
|---|---|
| `UNITS` | 只影响 `plumed.dat` 中的数值和 PLUMED 输出文件的单位。CP2K 在初始化时会把自己的单位换算告诉 PLUMED（源码中调用 `setMDLengthUnits`、`setMDEnergyUnits`、`setMDTimeUnits`），因此这里可以放心使用 Å 和 kJ/mol |
| `DISTANCE` | 原子序号从 1 开始，与 CP2K 一致；默认考虑周期性边界 |
| `COORDINATION` | `GROUPB=3-186:3` 表示 3, 6, …, 186（全部水 O）；开关函数与 CP2K 的 `NN`/`ND` 相同（PLUMED 中称 `NN`/`MM`） |
| `METAD` | `SIGMA`/`HEIGHT`/`PACE`/`BIASFACTOR` 分别对应 `SCALE`/`WW`/`NT_HILLS`/$\gamma$；数值与原生输入等价（0.35 bohr = 0.185 Å，4.0E-4 Ha = 1.05 kJ/mol，γ = 6）。`GRID_*` 让偏置在格点上计算，长时间模拟时更快 |
| `UPPER_WALLS` | 谐振上限墙 $\kappa(d-a)^2$，500 kJ mol⁻¹ Å⁻² ≈ 0.05 Ha/bohr² |
| `PRINT` | 每 20 步把 CV 和偏置写到 `COLVAR` |

> [!IMPORTANT]
> CP2K 初始化 PLUMED 时只传单位、时间步长、原子数、MPI 通信器和输入文件名（`src/metadynamics.F` 中的 `metadyn_initialise_plumed`），**不传温度**（没有调用 `setKbT`）。PLUMED 手册规定用 `BIASFACTOR` 做 well-tempered 时必须同时给 `TEMP`，所以 `plumed.dat` 中的 `TEMP=300` 不能省略，并且要与 `&MD/TEMPERATURE` 一致。PLUMED 的日志固定写在 `PLUMED.OUT` 中，出错时先看它。

### 6.5 运行与后处理

```bash
mpirun -np 64 cp2k.psmp -i plumed_metad.inp -o plumed_metad.out
grep -i "error" PLUMED.OUT
plumed sum_hills --hills HILLS --mintozero --outfile fes_plumed.dat
plumed sum_hills --hills HILLS --stride 100 --mintozero      # fes_0.dat, fes_1.dat, ... for convergence
```

`sum_hills` 输出的 CV 单位与 `UNITS` 一致（这里是 Å），能量单位是 kJ/mol。PLUMED 还提供与时间无关的重加权估计（Tiwary–Parrinello，`METAD ... CALC_RCT`）和 `REWEIGHT_BIAS`，可以把元动力学轨迹重加权到任意其他 CV 上，详见 PLUMED 手册。

用同样的步数分别跑原生版本和 PLUMED 版本，两条 $F(d)$ 应在统计误差内一致。这是检查单位换算是否正确的好办法。

续算时要注意：PLUMED 的偏置历史保存在它自己的 `HILLS` 文件里，不在 CP2K 的 `.restart` 中。用 `&EXT_RESTART` 续算 CP2K 的坐标和速度时，还要在 `plumed.dat` 的 `UNITS` 行之后加一行 `RESTART`，让 PLUMED 读回已有的 `HILLS` 并以追加方式写输出；否则它会把旧文件备份成 `bck.*` 再从零开始放高斯。

---

## 7. 实战三：blue moon 约束 MD

### 7.1 一个窗口的完整输入

`constrained_md_distance.inp`：

```fortran
@SET LEVEL ${LEVEL-dft}
@SET COORD ${COORD-nacl_eq.xyz}
@SET TARGET ${TARGET-2.80}

&GLOBAL
  PROJECT nacl_bm_${LEVEL}
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL

@INCLUDE 'nacl_${LEVEL}_force_eval.inc'

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS ${STEPS-10000}
    TIMESTEP 0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE CSVR
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
  &CONSTRAINT
    CONSTRAINT_INIT T
    SHAKE_TOLERANCE 1.0E-8
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR T
      TARGET [angstrom] ${TARGET}
    &END COLLECTIVE
    &LAGRANGE_MULTIPLIERS ON
      COMMON_ITERATION_LEVELS 1
      FILENAME =bluemoon.LagrangeMultLog
    &END LAGRANGE_MULTIPLIERS
  &END CONSTRAINT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
      &EACH
        MD 100
      &END EACH
    &END RESTART
    &RESTART_HISTORY OFF
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

### 7.2 关键词解释

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `&COLLECTIVE/COLVAR` | 1 | 被约束的 CV 编号 |
| `TARGET` | 窗口位置 | 带单位前缀 `[angstrom]`；不写单位时按 bohr |
| `INTERMOLECULAR` | `T` | 约束涉及的原子属于不同“分子”时必须打开，此时 COLVAR 中的序号是全局原子序号 |
| `MOLECULE` / `MOLNAME` | — | 分子内约束时，指定施加约束的分子类型（二者互斥），约束作用于该类型的分子 |
| `CONSTRAINT_INIT` | `T` | 在第一步之前就把初始坐标和速度投影到约束面上。默认只在第一步之后施加；起点离 `TARGET` 较远时必须打开 |
| `SHAKE_TOLERANCE` | 1E-8 | 默认 1E-6。拉格朗日乘子是要被平均的物理量，收紧容差可以减少噪声 |
| `&LAGRANGE_MULTIPLIERS` | `ON` | 打印键，默认级别 `HIGH`，在 `PRINT_LEVEL LOW` 下需要显式 `ON` |
| `FILENAME` | `=bluemoon.LagrangeMultLog` | 以 `=` 开头表示文件名完全按所写的来，不加项目名前缀 |
| `TARGET_GROWTH`、`TARGET_LIMIT` | — | 让 `TARGET` 随时间线性变化（slow growth / 拉伸），本章不用 |

> [!WARNING]
> CP2K 根据成键距离自动生成“分子”（`&TOPOLOGY` 的默认连接性）。Na⁺ 与第一壳层水分子的距离很短，可能被归到同一个分子里，而且归属会随模拟变化。用 `INTERMOLECULAR T` 并给出全局原子序号，可以完全避开这个问题。

`bluemoon.LagrangeMultLog` 中每一步有两行，格式如下（示意输出，数值仅供参考）：

```text
Shake  Lagrangian Multipliers:           -0.004213887
Rattle Lagrangian Multipliers:           -0.004198521
```

单位是原子单位（hartree/bohr）。对单个距离约束，取所有 `Shake` 行的平均即可，见第2.4节。

### 7.3 运行所有窗口

```bash
# chained: window k starts from the last frame of window k-1
bash run_bluemoon.sh
# cheap test at xTB level
LEVEL=xtb STEPS=20000 bash run_bluemoon.sh
# just create the directories and look
DRYRUN=1 bash run_bluemoon.sh
```

默认窗口为 2.60–5.20 Å，间隔 0.2 Å，共 14 个，每个 5 ps。每个窗口目录 `window_<t>/` 中都有 `bm.out`、`bluemoon.LagrangeMultLog` 和轨迹。串联方式（`CHAIN=1`）的好处是相邻窗口结构接近，`CONSTRAINT_INIT` 只需要很小的位移；缺点是必须依次运行。窗口相互独立的并行方式见第8节。

检查：

```bash
grep -L "PROGRAM ENDED AT" window_*/bm.out                       # windows that did not finish
grep -l "NOT converged" window_*/bm.out                          # SHAKE/RATTLE or SCF problems
for w in window_*/; do echo "$w $(grep -c '^Shake' $w/bluemoon.LagrangeMultLog)"; done
```

### 7.4 积分得到自由能曲线

```bash
python3 bluemoon_integrate.py window_*/ --skip 2000 --pmf -T 300 -o bluemoon_fes.dat --png bluemoon_fes.png
```

脚本对每个窗口去掉前 2000 条记录（1 ps）作为平衡段，计算自由能梯度 $\mathrm dA/\mathrm dr = -\langle\lambda\rangle$（下文沿用习惯叫“平均力”，注意它与物理意义上的平均力 $-\mathrm dA/\mathrm dr$ 差一个符号）和 5 块的块平均误差，按距离排序后用梯形法积分。`--pmf` 额外输出 $w(r) = A(r) + 2k_BT\ln r$。输出列为：`r(A)`、`dA/dr`（kJ mol⁻¹ Å⁻¹）、误差、`A(r)`、`w(r)`、样本数。

判读要点：

- 平均力在接触离子对极小处从负变正（$A$ 先降后升）；在能垒处平均力为零；
- 相邻窗口的平均力变化很陡时（通常在能垒附近），要在中间加窗口；
- 误差主要来自溶剂重排的慢弛豫，误差大的窗口应延长，而不是减小间隔；
- $A(r)$ 和 $w(r)$ 的差别在大 $r$ 处可以达到数 kJ/mol，与文献对比时要弄清对方报告的是哪一个。

### 7.5 伞形采样：用约束势代替约束

把 `&COLLECTIVE` 中加上 `&RESTRAINT`/`K`（第4.1节）就得到一个伞形采样窗口：CV 在 `TARGET` 附近涨落而不是被固定。各窗口的 CV 时间序列可以用 WHAM 或 MBAR 合并（例如 Grossfield 的 WHAM 程序或 pymbar），CP2K 自身也提供 `&FREE_ENERGY` 的 `METHOD UI`（umbrella integration）。伞形采样不需要质量度规修正，适合复杂 CV；代价是每个窗口要有足够的 CV 分布重叠。

### 7.6 一般 CV 的 blue moon

约束角度、配位数、组合 CV，或者同时有多个约束、固定原子时，$Z$ 不再是常数，必须做第2.4节的度规修正。修正需要每一步的 $\partial\xi/\partial\mathbf r_i$ 和二阶导数，CP2K 只输出 $\lambda$，不计算这些量，需要自己从轨迹计算。CP2K 的 trunk 文档新增了 [Constrained molecular dynamics](https://manual.cp2k.org/trunk/methods/sampling/constrained_dynamics.html) 一节，给出了完整公式，并介绍了主仓库 `tools/blue_moon/` 中的独立后处理脚本（处理单个约束的距离、角度、二面角、点到平面距离及其线性组合）。这个脚本不在 2026.2 的发行版中，需要从 GitHub 主分支获取，使用前请先读它的 README 了解适用范围。不想做这些修正时，更简单的做法是改用伞形采样或 WT-MetaD。示例中的 `bluemoon_integrate.py` 只适用于**单个原子间距离**。

---

## 8. xTB 预采样 + DFT 精修

DFT 元动力学最大的浪费是“试错”：`SCALE` 太小、墙放错、CV 不对，几十 ps 的 DFT 时间就白费了。推荐的流程如下：

| 步骤 | 层次 | 输入 | 得到什么 |
|---|---|---|---|
| 1. 平衡 | xTB | `nacl_equil.inp` | 合理的起始结构 `nacl_eq.xyz` |
| 2. 无偏 MD | xTB（可再用 DFT 跑 2–4 ps 核对） | `cv_monitor.inp` | CV 涨落 → `SCALE` |
| 3. WT-MetaD 预采样 | xTB，100 ps 量级 | `nacl_dissociation_wtmetad.inp -E LEVEL=xtb` | CV 的实际范围、墙的位置、是否需要第二个 CV、全 CV 区间的结构 |
| 4a. DFT WT-MetaD | DFT | 同一输入，`-E LEVEL=dft` | 定量的 $F(s)$，参数已经调好 |
| 4b. 或 DFT blue moon | DFT | `pick_windows.py` → `CHAIN=0 bash run_bluemoon.sh` | 每个窗口独立，可并行 |

### 8.1 从 xTB 轨迹为 DFT 窗口挑起点

```bash
python3 pick_windows.py nacl_wtmetad_xtb-pos-1.xyz --box 12.42 \
    --targets 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20
CHAIN=0 bash run_bluemoon.sh                                   # all windows here, or
for t in 2.60 2.80 3.00; do
  sbatch --export=ALL,STEP=bluemoon,LEVEL=dft,CHAIN=0,TARGETS=$t job_slurm.sh
done
```

脚本按最小镜像计算每一帧的 $d$(Na–Cl)，为每个目标距离挑最接近的一帧，写入 `window_<t>/start.xyz`，偏差超过 0.1 Å 时给出提示。窗口目录名直接使用你给的字符串，所以要和 `TARGETS` 写法一致（`2.60` 而不是 `2.6`）。DFT 窗口开始时 `CONSTRAINT_INIT T` 会把离子精确移到目标距离；但 xTB 的溶剂结构与 DFT 不同，每个窗口的平衡段（`--skip`）要留足，建议至少 1 ps。

### 8.2 需要注意什么

- **xTB 只用来“探路”**。GFN1-xTB 对液态水的结构和氢键强度只是定性描述，离子对的自由能曲线与 DFT 可能有明显差别。参数和起点可以从 xTB 继承，自由能数值不可以；
- **CV 在两个层次上要相同**：两个 include 文件中的 `&COLVAR` 保持一致，才能直接对比；
- **单位换算**：xTB 和 DFT 的 HILLS 都是 hartree/bohr，同一个脚本可以处理；
- **用 DFT 能量重加权 xTB 轨迹**：从原理上可以用 REFTRAJ（2026.2 用 `EVAL ENERGY`）算出 xTB 帧的 DFT 能量，再做自由能微扰 $\Delta A = -k_BT\ln\langle e^{-\beta(E_\mathrm{DFT}-E_\mathrm{xTB})}\rangle$。但在 188 个原子的凝聚相体系中，两种方法能量差的涨落远大于 $k_BT$，指数平均几乎不收敛。只在小体系或者能量差涨落很小时才考虑，否则直接用 DFT 采样；
- **更进一步**：用 DFT 标注 xTB/元动力学采到的结构，训练机器学习势再做长时间增强采样，是目前常见的做法（相关工具见 [第17章](17-tools-ecosystem.md)）。

---

## 9. 结果解读与文献对照

NaCl 在水中的平均力势在经典力场和 AIMD 研究中都有大量报道。曲线的定性形状是：

- 接触离子对（CIP）极小，大约在 2.7–2.9 Å；
- 3.5–4 Å 附近有一个去溶剂化能垒；
- 溶剂分隔离子对（SSIP）极小，大约在 4.5–5.2 Å；
- 更远处趋于平坦（在小盒子里受周期像影响）。

Timko、Bucher 和 Kuyucak 用 AIMD 约束 MD 研究了这个问题（J. Chem. Phys. 132, 114510 (2010)），是与本章示例最接近的参考。对比时注意：(1) 对方报告的是 $A(r)$ 还是 $w(r)$；(2) 泛函和色散校正；(3) 盒子大小和水分子数；(4) 每个窗口的模拟时间。本章示例使用 62 个水分子和较短的采样时间，主要用于教学，结果不应直接用来与文献比较能垒的具体数值。

---

## 参数选择与收敛建议

| 参数 | 建议 | 理由 |
|---|---|---|
| CV | 先用无偏 MD 和化学直觉确认它能区分反应物、产物和过渡区 | 偏置无法弥补选错的 CV |
| CV 个数 | 1–2 个；3 个以上用 PLUMED 的 OPES 或路径 CV | 需要填平的体积随维数指数增长 |
| `SCALE` | 无偏涨落标准差的 1/3–1/2 | 兼顾分辨率和填充速度 |
| `WW` | 0.2–0.5 $k_BT$（300 K 下约 0.5–1.2 kJ/mol） | 太高会损坏动力学，太低填充太慢 |
| `NT_HILLS` | 10–50 fs 对应的步数 | 应长于 CV 的自相关时间 |
| $\gamma$（`WTGAMMA`） | 使 $\gamma k_BT$ 与要跨越的能垒相当，常用 5–15 | 能垒在偏置下约降为 $1/\gamma$ |
| 墙 | 放在感兴趣区间之外、盒长一半之内；分析时舍去墙外区域 | 墙势不在 HILLS 中 |
| 恒温器 | CSVR，`TIMECON` 50–200 fs | 偏置持续注入能量 |
| 时间步长 | 0.5 fs（含 H）；氘化后可用 1 fs | 与 [第12章](12-aimd.md) 相同 |
| Blue moon 窗口 | 0.1–0.2 Å 间隔；能垒附近加密；每个窗口 ≥ 5 ps（DFT） | 平均力变化剧烈处积分误差大 |
| SCF | `EPS_SCF` ≤ 1E-6，ASPC 外推，OT + `OUTER_SCF` | 力的噪声会被偏置放大 |
| 统计 | 至少两段独立数据（或块平均）给出误差 | 只有一条曲线时无法判断收敛 |

---

## 常见问题与报错

**1. `When using Well-Tempered metadynamics, DELTA_T (or WTGAMMA) should be explicitly specified.`**
原因：写了 `WELL_TEMPERED T`，但没有给 $\Delta T$ 或 $\gamma$。解决：加 `DELTA_T 1500.0` 或 `WTGAMMA 6.0`。

**2. `Inconsistency between DELTA_T and WTGAMMA (both specified): please, verify that DELTA_T=(WTGAMMA-1)*TEMPERATURE`**
原因：两个都写了但互相矛盾（这里的温度是 `&MD/TEMPERATURE`）。解决：只保留一个。

**3. `Well-Tempered metadynamics not possible with extended-lagrangian formulation.`（或 `... with Adaptive hills.`）**
原因：WT 与 `LAGRANGE T` 或 `MIN_DISP` 同时使用。解决：二者取其一。

**4. `An error occurred in the specification of COLVAR for METAVAR. Specified COLVAR #(3) is larger than the maximum number of COLVARS defined in the SUBSYS (2) !`**
原因：`&METAVAR/COLVAR` 的编号超过了 `&SUBSYS` 中 `&COLVAR` 段的个数。常见于切换 include 文件后 CV 个数不同。解决：核对编号，让各层次的 include 文件定义相同的 CV。

**5. `Requested to use plumed for metadynamics, but cp2k was not compiled with plumed support.`**
原因：`USE_PLUMED T`，但 CP2K 没有编译 PLUMED 接口。解决：用 `cp2k.psmp -v` 确认 cp2kflags 中有 `plumed2`，没有就重新编译（第6.2节）。

**6. `NO PLUMED IN YOUR LD_LIBRARY_PATH?`**
原因：CP2K 编译时带了 PLUMED，但运行时找不到 PLUMED 动态库。解决：`source` PLUMED 或 toolchain 生成的环境脚本，或把 `libplumedKernel.so` 所在目录加入 `LD_LIBRARY_PATH`。

**7. `Shake NOT converged in ... iterations in the intermolecular constraint. CP2K continues but results could be meaningless.`**
原因：约束无法满足，常见于起点离 `TARGET` 太远而没有 `CONSTRAINT_INIT`，或时间步长太大、SCF 力有噪声。解决：打开 `CONSTRAINT_INIT T`，从相邻窗口的最后一帧起步，必要时减小 `TIMESTEP`。出现这条警告的窗口，其拉格朗日乘子不能用。

**8. 自由能曲线在大 $d$ 处一直下降，没有平台。**
原因：墙的位置超过了盒长一半，或者根本没有墙，最小镜像让 CV 失去意义。解决：墙放在 $L/2$ 以内，并舍去墙外区域；需要更大距离就用更大的盒子。

**9. 高斯宽度看起来比预期小很多（自由能曲线锯齿很多）。**
原因：`SCALE` 按 Å 写了，实际单位是 bohr。解决：所有 `SCALE`、`POSITION` 都换算成内部单位（1 Å = 1.8897 bohr）。

**10. WT-MetaD 跑了很久，高斯高度几乎没有下降。**
原因：CV 一直在探索新区域，可能是墙太远、CV 没有边界（例如 $d$ 在溶液中可以无限增大），或者 $\gamma$ 太大。解决：加墙或缩小区间，适当减小 $\gamma$。

**11. 正向和反向跨越得到的自由能曲线形状不同（滞后）。**
原因：有隐藏的慢变量没有被偏置（这里是水合壳层）。解决：加第二个 CV（`nacl_wtmetad_2d.inp`），或降低沉积速率（增大 `NT_HILLS`、减小 `WW`）。

**12. 2026.2 中重新使用旧的 REFTRAJ 输入，解析阶段就报错，指向 `EVAL_ENERGY_FORCES`（或 `EVAL_FORCES`）。**
原因：2026.2 删除了这两个逻辑关键字（PR #5401）。解决：改为 `&REFTRAJ/EVAL NONE|ENERGY|ENERGY_FORCES`（第3.7节）。`EVAL` 自 2024.1 起可用，新写法在 2024.1 及以后的版本中通用。

**13. SCF 不收敛导致元动力学中止，报 `SCF run NOT converged`。**
原因：偏置把体系推到了不常见的构型。解决：见 [第5章](05-scf-strategy.md)，适当增加 `MAX_SCF`/`OUTER_SCF`，检查 `CHARGE`/自旋；不要用 `IGNORE_CONVERGENCE_FAILURE` 掩盖问题，否则偏置会作用在错误的力上。更多报错见 [第18章](18-troubleshooting.md)。

**14. 在同一目录重跑后，自由能曲线明显变深，HILLS 的时间列从头又出现了一遍。**
原因：CP2K 以追加（APPEND）方式打开打印文件，旧的 `*-HILLS.metadynLog`、`*-COLVAR.metadynLog`、`*.LagrangeMultLog` 不会被覆盖，两次运行的数据叠在一起。解决：重跑前删除或改名旧文件（`run_bluemoon.sh` 会自动删除各窗口的旧 `bluemoon.LagrangeMultLog`）。续算时追加是正常的，只需去掉重复的行。

**15. PLUMED 版本一启动就中止，`PLUMED.OUT` 提示 well-tempered 需要温度。**
原因：`plumed.dat` 的 `METAD` 用了 `BIASFACTOR` 却没写 `TEMP`，而 CP2K 不会把温度传给 PLUMED（第6.4节）。解决：在 `METAD` 行加 `TEMP=300`（与 `&MD/TEMPERATURE` 一致）。

---

## 练习

1. 用 `cv_monitor.inp` 测量 $d$(Na–Cl) 和 CN 的涨落，按 1/3 和 1/2 标准差各设一套 `SCALE`，在 xTB 上各跑 50 ps WT-MetaD，比较两条 $F(d)$ 的分辨率和收敛速度。
2. 在 xTB 上分别用 $\gamma$ = 3、6、12 运行 1D WT-MetaD，画出 $d(t)$ 和高斯高度随时间的变化，讨论 $\gamma$ 对跨越频率和收敛的影响。
3. 把 `nacl_dissociation_wtmetad.inp` 的 `WELL_TEMPERED` 去掉（标准元动力学），用 `fes_from_hills.py --blocks 8` 观察自由能曲线怎样振荡而不收敛。
4. 用同样的参数分别跑原生和 PLUMED 版本（xTB 层次即可），比较两条 $F(d)$，验证第6.4节中的单位换算。
5. 用 `pick_windows.py` 从 xTB 轨迹挑起点，在 xTB 层次上做 blue moon 并积分；与同一层次 WT-MetaD 的曲线对比（记得区分 $A(r)$ 与 $w(r)$）。
6. 运行 2D WT-MetaD，比较投影得到的 $F(d)$ 与 1D 结果。哪一个的滞后更小？
7. 为下面的问题设计 CV 并写出 `&COLVAR` 段：(a) 水中 HCl 的质子转移（`COMBINE_COLVAR`）；(b) 甲醇在 Pt(111) 表面的脱附（`DISTANCE AXIS Z`）；(c) 丁烷的 C–C–C–C 二面角旋转（`TORSION`）。

---

## 参考资料

**官方文档**

- Metadynamics 方法页（甲酸/水在金红石表面的例子，`graph` 用法）：<https://manual.cp2k.org/trunk/methods/sampling/metadynamics.html>
- 输入参考 `&METADYN`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/FREE_ENERGY/METADYN.html>；`&METAVAR`、`&WALL`、`&PRINT/&HILLS` 在同一目录下
- 输入参考 `&COLVAR`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR.html>；`COORDINATION`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR/COORDINATION.html>；`COMBINE_COLVAR`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/COLVAR/COMBINE_COLVAR.html>
- 输入参考 `&CONSTRAINT`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/CONSTRAINT.html>；`&COLLECTIVE`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/CONSTRAINT/COLLECTIVE.html>；`&LAGRANGE_MULTIPLIERS`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/CONSTRAINT/LAGRANGE_MULTIPLIERS.html>
- Constrained molecular dynamics 方法页（trunk；乘子输出格式、blue moon 修正、`TARGET_GROWTH`、`tools/blue_moon`）：<https://manual.cp2k.org/trunk/methods/sampling/constrained_dynamics.html>
- 输入参考 `&REFTRAJ`（2026.2 的 `EVAL`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/MD/REFTRAJ.html>
- Molecular Dynamics：<https://manual.cp2k.org/trunk/methods/sampling/molecular_dynamics.html>；xTB：<https://manual.cp2k.org/trunk/methods/semiempiricals/xtb.html>
- CP2K 官方练习“Free Energy Profile of NaCl Dissociation”（约束 MD 与拉格朗日乘子积分）：<https://www.cp2k.org/exercises:2014_ethz_mmm:nacl_free_energy>
- CP2K 2026.2 release notes：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- PLUMED 文档：<https://www.plumed.org/doc>（`METAD`、`UNITS`、`COORDINATION`、`UPPER_WALLS`、`RESTART` 和 `sum_hills` 各页）；PLUMED-NEST（可复现增强采样输入库）：<https://www.plumed-nest.org/>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/Fist/regtest-15/metadyn*.inp`（各种元动力学选项的回归测试）、`tests/Fist/regtest-plumed2/`（PLUMED 接口测试）、`src/metadynamics.F`（WT 高度缩放）、`src/metadyn_tools/graph.F`（`graph` 工具）
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库，含元动力学专题（无 LICENSE，只能链接）
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含 CP2K constrained MD 与 PMF 计算
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含元动力学和 xTB 章节（2025.2 验证）
- [ishikawa-group/cp2k-tutorial](https://github.com/ishikawa-group/cp2k-tutorial)：日文入门教程，含元动力学
- [svarner9/CP2K-Workshop](https://github.com/svarner9/CP2K-Workshop)：用 ASE 和 PySAGES 驱动 CP2K 做增强采样
- [JunboLu/CP2K_kit](https://github.com/JunboLu/CP2K_kit)：含热力学积分等功能的 CP2K 工具箱

**本仓库脚本**：示例目录中的 [fes_from_hills.py](../examples/15-metadynamics-nacl/fes_from_hills.py)、[bluemoon_integrate.py](../examples/15-metadynamics-nacl/bluemoon_integrate.py)、[pick_windows.py](../examples/15-metadynamics-nacl/pick_windows.py)、[run_bluemoon.sh](../examples/15-metadynamics-nacl/run_bluemoon.sh)；CP2K 通用脚本 [../scripts/md_analyze.py](../scripts/md_analyze.py)（MD 轨迹分析）、[../scripts/restart_update.py](../scripts/restart_update.py)（续算）、[../scripts/job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh)（作业模板）。Python 脚本用 `python <脚本> -h` 查看用法，shell 脚本的用法写在文件头注释里。

**关键文献**

- A. Laio, M. Parrinello, *PNAS* **99**, 12562 (2002)，元动力学：<https://doi.org/10.1073/pnas.202427399>
- A. Barducci, G. Bussi, M. Parrinello, *Phys. Rev. Lett.* **100**, 020603 (2008)，well-tempered 元动力学：<https://doi.org/10.1103/PhysRevLett.100.020603>
- M. Iannuzzi, A. Laio, M. Parrinello, *Phys. Rev. Lett.* **90**, 238302 (2003)，扩展拉格朗日元动力学：<https://doi.org/10.1103/PhysRevLett.90.238302>
- O. Valsson, P. Tiwary, M. Parrinello, *Annu. Rev. Phys. Chem.* **67**, 159 (2016)：<https://doi.org/10.1146/annurev-physchem-040215-112229>
- G. Bussi, A. Laio, *Nat. Rev. Phys.* **2**, 200 (2020)，元动力学综述：<https://doi.org/10.1038/s42254-020-0153-0>
- P. Tiwary, M. Parrinello, *J. Phys. Chem. B* **119**, 736 (2015)，与时间无关的重加权估计：<https://doi.org/10.1021/jp504920s>
- E. A. Carter, G. Ciccotti, J. T. Hynes, R. Kapral, *Chem. Phys. Lett.* **156**, 472 (1989)，blue moon 约束动力学：<https://doi.org/10.1016/S0009-2614(89)87314-2>
- M. Sprik, G. Ciccotti, *J. Chem. Phys.* **109**, 7737 (1998)，约束 MD 的自由能与度规修正：<https://doi.org/10.1063/1.477419>
- J. Timko, D. Bucher, S. Kuyucak, *J. Chem. Phys.* **132**, 114510 (2010)，NaCl 在水中解离的 AIMD：<https://doi.org/10.1063/1.3360310>
- G. A. Tribello *et al.*, *Comput. Phys. Commun.* **185**, 604 (2014)，PLUMED 2：<https://doi.org/10.1016/j.cpc.2013.09.018>
- The PLUMED consortium, *Nat. Methods* **16**, 670 (2019)，PLUMED-NEST：<https://doi.org/10.1038/s41592-019-0506-8>
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.* **13**, 1989 (2017)，GFN1-xTB：<https://doi.org/10.1021/acs.jctc.7b00118>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](14-vibration-thermo.md) | [返回目录](../README.md) | [下一章 →](16-qmmm-xtb.md)
