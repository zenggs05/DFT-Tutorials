# 第9章 杂化泛函与 ADMM（Hybrid Functionals & ADMM in CP2K）

本章讲在 CP2K 中做杂化泛函（PBE0、HSE06、B3LYP）计算的完整流程：`&XC/&HF` 的 `FRACTION`、`&SCREENING`、`&MEMORY`、`&INTERACTION_POTENTIAL` 各自控制什么；周期体系为什么必须用截断库仑（`TRUNCATED`）或短程（`SHORTRANGE`）算符，`CUTOFF_RADIUS` 与晶胞尺寸的关系；如何用辅助密度矩阵方法 ADMM（`&AUXILIARY_DENSITY_MATRIX_METHOD`、`BASIS_SET AUX_FIT`、`BASIS_ADMM_UZH`）把 Hartree–Fock 交换的开销降低一个量级。示例覆盖孤立水分子 PBE0、32 水分子盒子 PBE0-TC-LRC（PBE → PBE0 重启）和 64 原子 Si 的 HSE06，并介绍 k 点杂化 RI-HFXk 以及 CP2K 2026.2 新增的 ACE 加速和 HSE 短程截断半径关键字。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出分子和周期体系的 PBE0 / HSE06 / B3LYP 输入，并说清每个 `&HF` 子段的作用；
- 根据晶胞尺寸选择 `CUTOFF_RADIUS`，看懂 “larger than half the minimal cell dimension” 警告；
- 正确配置 ADMM：选择辅助基组、`ADMM_TYPE`、`EXCH_CORRECTION_FUNC`，并评估 ADMM 误差；
- 选用与泛函匹配的基组和赝势（UZH 的 `-PBE0-`/`-HYB-` 系列）；
- 用“先 PBE、再杂化”的重启链和 `SCREEN_ON_INITIAL_P` 加速收敛；
- 估算并设置 `MAX_MEMORY`，合理安排 MPI × OpenMP，读懂 `HFX_MEM_INFO` 输出；
- 了解 k 点杂化（RI-HFXk，≥ 2024.1）、ACE 加速（2026.2）和 libGint GPU 后端（2026.2）的适用条件。

什么时候值得用杂化泛函：

| 目标 | GGA（PBE） | 杂化（PBE0/HSE06） |
|---|---|---|
| 带隙、能级对齐、缺陷能级 | 系统性低估 | 明显改善，常用 |
| 极化子、局域电荷（TiO2、CeO2 中的电子/空穴） | 常过度离域 | 更合理（或用 DFT+U，见[第8章](08-dispersion-spin-dftu.md)） |
| 反应能垒 | 常偏低 | 通常更好 |
| 液态水结构（AIMD） | 过度结构化 | 配合色散校正有改善 |
| 金属、大体系结构优化 | 首选 | 代价高、收益小 |

计算量参考：同一体系的杂化泛函单点通常是 GGA 的数倍到数十倍（取决于 ADMM、截断半径、基组弥散程度），第一步 SCF 由于要计算并存储双电子积分（ERI）尤其慢。示例输入集：[../examples/09-hybrid-admm/README.md](../examples/09-hybrid-admm/README.md)。

---

## 2. 原理速览

### 2.1 杂化泛函

PBE0 把 25% 的 Hartree–Fock（HF）精确交换混入 PBE：

$$
E_{xc}^\mathrm{PBE0} = \tfrac14\,E_x^\mathrm{HF} + \tfrac34\,E_x^\mathrm{PBE} + E_c^\mathrm{PBE}
$$

HSE06 只在短程部分混入精确交换，用误差函数分割库仑算符 $1/r = \mathrm{erfc}(\omega r)/r + \mathrm{erf}(\omega r)/r$：

$$
E_{xc}^\mathrm{HSE06} = \tfrac14\,E_x^\mathrm{HF,SR}(\omega) + \tfrac34\,E_x^\mathrm{PBE,SR}(\omega) + E_x^\mathrm{PBE,LR}(\omega) + E_c^\mathrm{PBE},\qquad \omega = 0.11\ \text{bohr}^{-1}
$$

CP2K 中 HF 交换能为

$$
E_x^\mathrm{HF} = -\frac12\sum_{\mu\nu\lambda\sigma} P_{\mu\sigma}P_{\nu\lambda}\,(\mu\nu|g(r_{12})|\lambda\sigma)
$$

其中 $g(r)$ 是交换算符（`POTENTIAL_TYPE`），$(\mu\nu|g|\lambda\sigma)$ 是四中心双电子积分（ERI）。

### 2.2 周期体系的难点：截断库仑

周期体系中 $1/r$ 是长程的，Γ 点 HF 交换的晶格求和收敛很慢且依赖周期像。CP2K 的标准做法（Guidon 等，2009）是用截断库仑算符：

$$
g_\mathrm{TC}(r) = \begin{cases} 1/r, & r \le R_c \\ 0, & r > R_c \end{cases}
$$

在最小像（minimum image）约定下，$R_c$ 必须不超过晶胞最小宽度的一半。CP2K 用晶面间距（而不是晶胞边长）判断，对斜晶胞更严格；k 点计算则按 Born–von Kármán 超胞（k 网格 × 晶胞）判断。被截掉的长程部分可以用 PBE 交换空穴模型补回（`PBE_HOLE_T_C_LR`），这就是 PBE0-TC-LRC。HSE06 的短程算符本身会衰减，但在 ω = 0.11 bohr⁻¹ 时衰减较慢，同样需要考虑截断。

### 2.3 积分筛选与存储

Schwarz 不等式 $|(\mu\nu|\lambda\sigma)| \le \sqrt{(\mu\nu|\mu\nu)}\sqrt{(\lambda\sigma|\lambda\sigma)}$ 用来跳过可忽略的积分（`EPS_SCHWARZ`）；再乘上密度矩阵元的大小（`SCREEN_ON_INITIAL_P`）还能筛掉更多。第一步计算出的 ERI 压缩后存在内存（`MAX_MEMORY`）或磁盘（`MAX_DISK_SPACE`）中，后续 SCF 步直接复用；内存不够的部分每步重新计算（“on the fly”），会显著变慢。

### 2.4 ADMM

ADMM 把主基组的密度矩阵 $\mathbf P$ 投影到一个小而紧凑的辅助基组上得到 $\hat{\mathbf P}$，在辅助基组中算 HF 交换，再用一个 GGA 交换泛函修正两者的差：

$$
E_x^\mathrm{HF}[\mathbf P] \approx E_x^\mathrm{HF}[\hat{\mathbf P}] + \Big(E_x^\mathrm{DFT}[\mathbf P] - E_x^\mathrm{DFT}[\hat{\mathbf P}]\Big)
$$

ERI 的数量随基函数数目的四次方增长，辅助基组小很多，所以 ADMM 往往能把 HFX 部分加速一个数量级，而能量误差一般在可接受范围内（需要测试，见第 6.4 节）。

---

## 3. 杂化泛函的输入骨架

### 3.1 `&XC` 的两种写法

**写法一：快捷方式。** `&XC_FUNCTIONAL PBE0` 和 `&XC_FUNCTIONAL B3LYP` 会自动设置 GGA 部分的缩放，并把 `&HF/FRACTION` 设为 0.25（PBE0）或 0.20（B3LYP）（2026.2 源码 `input_cp2k_check.F` 中的展开逻辑）。但截断半径、算符类型、筛选阈值和内存都**不会**自动设置，仍需要自己写 `&HF` 段：

```text
&XC
  &XC_FUNCTIONAL PBE0
  &END XC_FUNCTIONAL
  &HF
    &SCREENING
      EPS_SCHWARZ 1.0E-6
    &END SCREENING
    &MEMORY
      MAX_MEMORY 2000
    &END MEMORY
  &END HF
&END XC
```

**写法二：手动组合（推荐，显式且便于修改）。**

```text
&XC
  &XC_FUNCTIONAL
    &PBE
      SCALE_X 0.75
      SCALE_C 1.0
    &END PBE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 0.25
    ...
  &END HF
&END XC
```

> [!WARNING]
> 使用 LibXC 的杂化泛函（如 `&HYB_GGA_XC_HSE06`、`&HYB_GGA_XC_PBEH`）时，CP2K **不会**自动生成 `&HF` 段（输入手册原文：“CP2K does NOT set it up automatically not even for LibXC functionals”）。必须自己写 `FRACTION` 和对应的 `&INTERACTION_POTENTIAL`，否则只算了 GGA 部分。

### 3.2 `&HF` 关键字表

| 位置 | 关键字 | 默认值 | 建议 | 说明 |
|---|---|---|---|---|
| `&HF` | `FRACTION` | 1.0 | PBE0/HSE06: 0.25 | 精确交换比例。默认 1.0 对应纯 HF，杂化泛函务必改写 |
| `&HF/&INTERACTION_POTENTIAL` | `POTENTIAL_TYPE` | `COULOMB` | 分子 `COULOMB`；周期 `TRUNCATED`；HSE `SHORTRANGE` | 还有 `LONGRANGE`、`MIX_CL`、`MIX_CL_TRUNC`、`GAUSSIAN` 等 |
| 同上 | `CUTOFF_RADIUS` | 无（Å） | ≤ L/2 | `TRUNCATED` 必需；2026.2 起对 `SHORTRANGE` 也生效 |
| 同上 | `T_C_G_DATA` | `t_c_g.dat` | 默认 | 截断库仑所需的数据表，在 data 目录中 |
| 同上 | `OMEGA` | 0.0 | HSE06: 0.11 | 范围分离参数，单位 bohr⁻¹ |
| `&HF/&SCREENING` | `EPS_SCHWARZ` | 1.0E-10 | 1E-6 ~ 1E-8 | 积分筛选阈值，越小越准越贵 |
| 同上 | `EPS_SCHWARZ_FORCES` | 1.0E-6 | 约等于期望的力精度 | 力计算中的筛选阈值 |
| 同上 | `SCREEN_ON_INITIAL_P` | F | 有重启波函数时 T | 用初始密度矩阵筛选积分，**需要从重启文件起步** |
| 同上 | `SCREEN_P_FORCES` | T | 默认 | 力计算中用密度矩阵筛选 |
| `&HF/&MEMORY` | `MAX_MEMORY` | 512 | 按节点内存计算 | **每个 MPI 进程**可用于 HFX 的内存（MiB），OpenMP 线程共享 |
| 同上 | `EPS_STORAGE_SCALING` | 1.0 | 0.1 | 压缩存储阈值为 `EPS_SCHWARZ × EPS_STORAGE_SCALING` |
| 同上 | `MAX_DISK_SPACE` | 0 | 一般不用 | 允许把 ERI 存到磁盘（MiB），配合 `STORAGE_LOCATION` |
| `&HF/&HF_INFO` | — | MEDIUM 级别打印 | — | 打印 HFX 相关信息 |
| `&HF/&ACE` | `ACTIVE`、`REBUILD_FREQUENCY` | F、20 | 见第 9 节 | 2026.2 新增，仅对角化 SCF |
| `&HF` | `HFX_LIBRARY` | `LIBINT` | — | 2026.2 新增 `LIBGINT`（CUDA GPU），需编译支持 |

### 3.3 基组与赝势要与泛函匹配

GTH 赝势是针对特定泛函生成的。UZH 系列为杂化泛函提供了专门的赝势和基组：

| 用途 | 基组（`BASIS_MOLOPT_UZH`） | 赝势（`POTENTIAL_UZH`） | ADMM 辅助基组（`BASIS_ADMM_UZH`） |
|---|---|---|---|
| PBE | `DZVP-MOLOPT-PBE-GTH-q6`（别名 `-GGA-`） | `GTH-PBE-q6`（别名 `GTH-GGA-q6`） | — |
| PBE0 / 其它杂化 | `DZVP-MOLOPT-PBE0-GTH-q6`（别名 `DZVP-MOLOPT-HYB-GTH-q6`） | `GTH-PBE0-q6`（别名 `GTH-HYB-q6`） | `admm-dz-q6`、`admm-dzp-q6`、`admm-tzp-q6`、`admm-tz2p-q6` |

（以 O 为例，其它元素把 `q6` 换成对应价电子数；同一 KIND 中基组和赝势的 `-qN` 必须一致。）

HSE06 也是以 PBE 为基础的杂化，通常同样使用 `-PBE0-`/`-HYB-` 赝势；这是近似做法，正式研究中可以与 GTH-PBE 赝势的结果对比。老的 `BASIS_MOLOPT` 族可以搭配 `BASIS_ADMM` 中的 `cFIT3`、`cpFIT3`、`FIT3` 等辅助基组（H–Ar 常见元素），`BASIS_ADMM_MOLOPT` 补充了 Li、Al、K 等元素的 `cFITn`/`FITn`。

### 3.4 色散校正与杂化泛函

杂化泛函同样需要色散校正。`&VDW_POTENTIAL` 的写法与第8章相同，只是 `REFERENCE_FUNCTIONAL` 要换成杂化泛函名（如 `PBE0`、`B3LYP`）。手写组合泛函时 CP2K 无法自动识别，务必显式写出。

---

## 4. 例 1：孤立水分子 PBE0 + ADMM

### 4.1 输入文件

[../examples/09-hybrid-admm/h2o_pbe0_admm.inp](../examples/09-hybrid-admm/h2o_pbe0_admm.inp)：

```text
@SET PROJ h2o_pbe0_admm

&GLOBAL
  PROJECT ${PROJ}
  RUN_TYPE GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    BASIS_SET_FILE_NAME BASIS_ADMM_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &AUXILIARY_DENSITY_MATRIX_METHOD
      METHOD BASIS_PROJECTION
      ADMM_PURIFICATION_METHOD MO_DIAG
      EXCH_CORRECTION_FUNC DEFAULT
    &END AUXILIARY_DENSITY_MATRIX_METHOD
    &MGRID
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
      EXTRAPOLATION ASPC
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 30
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL
        &PBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE COULOMB
        &END INTERACTION_POTENTIAL
        &SCREENING
          EPS_SCHWARZ 1.0E-8
          SCREEN_ON_INITIAL_P F
        &END SCREENING
        &MEMORY
          MAX_MEMORY 1000
        &END MEMORY
      &END HF
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &COORD
      O  0.000000  0.000000  0.119000
      H  0.000000  0.763000 -0.477000
      H  0.000000 -0.763000 -0.477000
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET ORB TZVP-MOLOPT-PBE0-GTH-q6
      BASIS_SET AUX_FIT admm-dzp-q6
      POTENTIAL GTH-PBE0-q6
    &END KIND
    &KIND H
      BASIS_SET ORB TZVP-MOLOPT-PBE0-GTH-q1
      BASIS_SET AUX_FIT admm-dzp-q1
      POTENTIAL GTH-PBE0-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 100
  &END GEO_OPT
&END MOTION
```

### 4.2 要点

| 设置 | 说明 |
|---|---|
| 两个 `BASIS_SET_FILE_NAME` | 主基组来自 `BASIS_MOLOPT_UZH`，ADMM 辅助基组来自 `BASIS_ADMM_UZH`，可以写多行 |
| `BASIS_SET ORB ...` / `BASIS_SET AUX_FIT ...` | 同一 KIND 中用类型标签区分主基组（`ORB`，默认）和 ADMM 辅助基组（`AUX_FIT`） |
| `POTENTIAL_TYPE COULOMB` | 孤立分子（`PERIODIC NONE`）不需要截断 |
| `EPS_SCHWARZ 1.0E-8` | 小分子积分少，可以设得更紧 |
| `METHOD BASIS_PROJECTION` + `ADMM_PURIFICATION_METHOD MO_DIAG` | 经典 ADMM1 组合，手册建议 MD 中使用 `MO_DIAG` 和 `BASIS_PROJECTION` |
| `EXCH_CORRECTION_FUNC DEFAULT` | 按 `INTERACTION_POTENTIAL` 自动选择匹配的 PBE 交换修正（见第 6.3 节） |

运行与检查：

```bash
mpirun -np 4 cp2k.psmp -i h2o_pbe0_admm.inp -o h2o_pbe0_admm.out
grep "Hartree-Fock Exchange energy:" h2o_pbe0_admm.out | tail -1
grep "ADMM|" h2o_pbe0_admm.out | head
grep "HFX_MEM_INFO" h2o_pbe0_admm.out | head -12
grep "GEOMETRY OPTIMIZATION COMPLETED" h2o_pbe0_admm.out
```

能量分解中会出现 `Hartree-Fock Exchange energy:` 一行；`ADMM|` 开头的行汇报 ADMM 设置。PBE0 下水分子的 O–H 键长通常比 PBE 略短（更接近实验值 0.957 Å），可以与第7章的 PBE 结果对比。

---

## 5. 例 2：周期水盒子 PBE0-TC-LRC（PBE → PBE0 两步）

### 5.1 为什么分两步

杂化泛函的 SCF 每一步都很贵。先用 PBE 收敛得到波函数，再以它为初猜开始 PBE0，可以：

1. 大幅减少杂化 SCF 的步数；
2. 打开 `SCREEN_ON_INITIAL_P T`：用初始（PBE）密度矩阵筛选积分，减少 ERI 数量和内存。没有重启文件时不能打开这个选项。

两步必须使用**相同的主基组和赝势**（这里都用 `-PBE0-` 系列），否则 PBE 的波函数不是合法的重启文件。

### 5.2 第一步：PBE

[../examples/09-hybrid-admm/water32_pbe_step1.inp](../examples/09-hybrid-admm/water32_pbe_step1.inp) 是普通的 PBE 单点（OT，Γ 点），`PROJECT water32_pbe`，结束后得到 `water32_pbe-RESTART.wfn`。

```bash
mpirun -np 32 cp2k.psmp -i water32_pbe_step1.inp -o water32_pbe_step1.out
ls -l water32_pbe-RESTART.wfn
```

### 5.3 第二步：PBE0-TC-LRC + ADMM

[../examples/09-hybrid-admm/water32_pbe0_admm_truncated.inp](../examples/09-hybrid-admm/water32_pbe0_admm_truncated.inp) 的核心部分：

```text
@SET RC 4.5
...
&DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
  BASIS_SET_FILE_NAME BASIS_ADMM_UZH
  POTENTIAL_FILE_NAME POTENTIAL_UZH
  WFN_RESTART_FILE_NAME water32_pbe-RESTART.wfn
  &AUXILIARY_DENSITY_MATRIX_METHOD
    METHOD BASIS_PROJECTION
    ADMM_PURIFICATION_METHOD MO_DIAG
    EXCH_CORRECTION_FUNC DEFAULT
  &END AUXILIARY_DENSITY_MATRIX_METHOD
  &SCF
    SCF_GUESS RESTART
    EPS_SCF 1.0E-6
    MAX_SCF 20
    &OT
      MINIMIZER DIIS
      PRECONDITIONER FULL_SINGLE_INVERSE
    &END OT
    &OUTER_SCF
      EPS_SCF 1.0E-6
      MAX_SCF 10
    &END OUTER_SCF
  &END SCF
  &XC
    &XC_FUNCTIONAL
      &PBE
        SCALE_X 0.75
        SCALE_C 1.0
      &END PBE
      &PBE_HOLE_T_C_LR
        CUTOFF_RADIUS ${RC}
        SCALE_X 0.25
      &END PBE_HOLE_T_C_LR
    &END XC_FUNCTIONAL
    &HF
      FRACTION 0.25
      &INTERACTION_POTENTIAL
        POTENTIAL_TYPE TRUNCATED
        CUTOFF_RADIUS ${RC}
        T_C_G_DATA t_c_g.dat
      &END INTERACTION_POTENTIAL
      &SCREENING
        EPS_SCHWARZ 1.0E-6
        SCREEN_ON_INITIAL_P T
      &END SCREENING
      &MEMORY
        MAX_MEMORY 2500
      &END MEMORY
    &END HF
  &END XC
&END DFT
```

| 设置 | 说明 |
|---|---|
| `CUTOFF_RADIUS 4.5` | 盒子边长 9.8528 Å，L/2 = 4.93 Å，取略小的 4.5 Å |
| `&PBE_HOLE_T_C_LR` | 用 PBE 交换空穴模型补上 $r > R_c$ 被截掉的 25% 交换，`CUTOFF_RADIUS` 必须与 HF 段相同 |
| `&PBE SCALE_X 0.75` | 短程部分仍是 PBE0 的 75% PBE 交换 |
| `SCREEN_ON_INITIAL_P T` | 依赖 `SCF_GUESS RESTART` 读入的 PBE 密度矩阵 |
| `MAX_MEMORY 2500` | 每个 MPI 进程 2.5 GB；16 个进程共约 40 GB，按节点内存调整 |
| `RUN_TYPE ENERGY_FORCE` | 同时算力，可检查力的筛选设置 |

> [!NOTE]
> 不加 `PBE_HOLE_T_C_LR` 也可以运行（即“截断 PBE0”），在 $R_c$ 足够大时差别不大，CP2K 自带的部分基准输入就是这样写的。加上长程修正后，结果对 $R_c$ 的依赖更弱，是 Guidon 等推荐的 PBE0-TC-LRC 形式。

运行与检查：

```bash
export OMP_NUM_THREADS=2
mpirun -np 16 cp2k.psmp -i water32_pbe0_admm_truncated.inp -o water32_pbe0.out
grep -i "WARNING" water32_pbe0.out | sort | uniq -c
grep "HFX_MEM_INFO| Number of sph. ERI's calculated on the fly" water32_pbe0.out | head -1
grep "HFX_MEM_INFO| Total memory consumption ERI's RAM" water32_pbe0.out | head -1
grep "Hartree-Fock Exchange energy:" water32_pbe0.out | tail -1
grep "ENERGY| Total FORCE_EVAL" water32_pbe0.out
```

判断标准：

- 不出现 `cutoff radius ... is larger than half the minimal cell dimension` 警告；
- `Number of sph. ERI's calculated on the fly` 为 0（或很小）说明 `MAX_MEMORY` 足够，所有积分都存在内存里；数值很大就意味着每步都在重算积分，应增大 `MAX_MEMORY`、减少每节点 MPI 进程数（改用更多 OpenMP 线程）或增加节点；
- SCF 从 PBE 重启后通常十几步内收敛；
- 第一步 SCF 远慢于之后各步（ERI 计算），这是正常的。

> [!TIP]
> 示例的 `water32.xyz` 是在 fcc 格点上随机取向生成的初始构型，没有平衡，只用于演示输入。真实研究应先用经典 MD 或 PBE-D3 AIMD 平衡结构（见第12章）。

---

## 6. ADMM 详解

### 6.1 `ADMM_TYPE` 快捷方式

`&AUXILIARY_DENSITY_MATRIX_METHOD/ADMM_TYPE` 会一次性设好 `METHOD`、`ADMM_PURIFICATION_METHOD` 和 `EXCH_SCALING_MODEL`。根据 2026.2 源码（`cp_control_utils.F`），对应关系为：

| `ADMM_TYPE` | `METHOD` | `ADMM_PURIFICATION_METHOD` | `EXCH_SCALING_MODEL` | 文献 |
|---|---|---|---|---|
| `NONE`（默认） | 用户自己设置 | 用户自己设置 | 用户自己设置 | — |
| `ADMM1` | `BASIS_PROJECTION` | 取输入中的值 | `NONE` | Guidon 2010 |
| `ADMM2` | `BASIS_PROJECTION` | `NONE` | `NONE` | Guidon 2010 |
| `ADMMS` | `CHARGE_CONSTRAINED_PROJECTION` | `NONE` | `MERLOT` | Merlot 2014 |
| `ADMMP` | `BASIS_PROJECTION` | `NONE` | `MERLOT` | Merlot 2014 |
| `ADMMQ` | `CHARGE_CONSTRAINED_PROJECTION` | `NONE` | `NONE` | Merlot 2014 |

同时写了 `ADMM_TYPE` 和单独的 `METHOD` 等关键字时，后者会被覆盖，CP2K 给出警告。官方方法文档的示例使用 `ADMM_TYPE ADMMS`，并指出 `ADMM2` 在 ground-state 之外的方法（如部分激发态、后 HF）中支持最广。本章示例分别使用了显式的 ADMM1 组合（`BASIS_PROJECTION` + `MO_DIAG`）和 `ADMM_TYPE ADMM2`。

### 6.2 辅助基组怎么选

| 辅助基组文件 | 名称示例 | 搭配 |
|---|---|---|
| `BASIS_ADMM_UZH` | `admm-dz-qN`、`admm-dzp-qN`、`admm-tzp-qN`、`admm-tz2p-qN` | UZH 主基组（新项目推荐） |
| `BASIS_ADMM` | `cFIT3`、`cpFIT3`、`FIT3`、`pFIT3`、`aug-cFIT3` … | 老的 MOLOPT 主基组 |
| `BASIS_ADMM_MOLOPT` | `cFIT4`、`FIT5` …（Li、Al、K 等元素） | MOLOPT 主基组 |

辅助基组本身就是近似的一部分：太小则交换修正项变大、精度下降；太大则加速有限。经验做法是从 `admm-dzp` 起步，在一个小的代表性体系上与 `admm-tzp` 和不用 ADMM 的结果比较。

### 6.3 `EXCH_CORRECTION_FUNC`

| 取值 | 含义 |
|---|---|
| `DEFAULT`（默认） | 按 HF 的 `POTENTIAL_TYPE` 自动选择 PBE 类修正：`COULOMB` → PBE 交换；`SHORTRANGE` → 短程 PBE 交换（`XWPBE`，带相同 ω）；`TRUNCATED` → `PBE_HOLE_T_C_LR`（带相同 $R_c$）与 `XWPBE` 的组合（2026.2 源码 `hfx_admm_utils.F`） |
| `PBEX` | 总是用完整的 PBE 交换，不考虑算符 |
| `OPTX`、`BECKE88X` | OPTX 或 B88 交换（与 B3LYP 等 B88 型交换配合） |
| `*_LIBXC` | 对应泛函的 LibXC 实现 |
| `NONE` | 不做修正 |

官方文档的原则是：修正泛函应与主 XC 设置中的交换部分一致。PBE 类杂化 + 截断/短程算符时，用 `DEFAULT` 最省心。

### 6.4 检查 ADMM 误差

1. 取一个小的代表性体系（如几个水分子或小超胞）；
2. 分别用 ADMM（`admm-dzp`、`admm-tzp`）和不用 ADMM 计算目标性质（反应能、带隙、力）；
3. ADMM 误差远小于你关心的物理效应时才在大体系上使用 ADMM。

ADMM 只加速交换积分，主基组、网格（`CUTOFF`）和 SCF 阈值的收敛测试一个都不能省。

---

## 7. 例 3：HSE06 与短程截断半径（Si64）

### 7.1 输入要点

[../examples/09-hybrid-admm/si_hse06_admm.inp](../examples/09-hybrid-admm/si_hse06_admm.inp) 把 PBE 预计算和 HSE06 放在同一个文件里，用 `@SET STAGE PBE` / `HSE06` 切换，HSE06 阶段从 PBE 波函数重启。HSE06 部分：

```text
&XC
  &XC_FUNCTIONAL
    &PBE
      SCALE_X 0.0
      SCALE_C 1.0
    &END PBE
    &XWPBE
      SCALE_X -0.25
      SCALE_X0 1.0
      OMEGA 0.11
    &END XWPBE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 0.25
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE SHORTRANGE
      OMEGA 0.11
      CUTOFF_RADIUS ${RC}
    &END INTERACTION_POTENTIAL
    &SCREENING
      EPS_SCHWARZ 1.0E-6
      SCREEN_ON_INITIAL_P T
    &END SCREENING
    &MEMORY
      MAX_MEMORY 3000
    &END MEMORY
  &END HF
&END XC
```

`&XWPBE` 是 ωPBE 交换：`SCALE_X0 1.0` 给出完整的 PBE 交换，`SCALE_X -0.25` 减去 25% 的短程 PBE 交换；`&PBE` 只保留关联（`SCALE_X 0.0`）。合起来正好是 2.1 节的 HSE06 表达式。也可以用 LibXC 的 `&HYB_GGA_XC_HSE06` 代替 `&PBE` + `&XWPBE`，但 `&HF` 段（`FRACTION 0.25`、`SHORTRANGE`、`OMEGA 0.11`）仍要自己写。

### 7.2 `CUTOFF_RADIUS` 与版本差异

| 版本 | `SHORTRANGE` 下的截断半径 |
|---|---|
| ≤ 2026.1 | 忽略输入的 `CUTOFF_RADIUS`，由 $\mathrm{erfc}(\omega r)/r = \epsilon_\mathrm{Schwarz}$ 自动求解 |
| 2026.2（PR #4945） | 显式给出 `CUTOFF_RADIUS` 时使用输入值；省略时仍按上式自动求解 |

按 PR #4945 的说明，ω = 0.11 时 `EPS_SCHWARZ 1.0E-6` 自动得到约 14.28 Å，`1.0E-4` 得到约 10.52 Å，都超过很多晶胞的 L/2，于是出现：

```text
Periodic Hartree Fock calculation requested with the use of a truncated or shortrange
potential. The cutoff radius ( x A) is larger than half the minimal cell dimension ( y A).
This may lead to unphysical total energies. Reduce the cutoff radius in order to avoid
possible problems.
```

2026.2 中显式设置比自动值小的 `CUTOFF_RADIUS` 时，会出现另一条警告，提醒该值小于满足 $\mathrm{erfc}(\omega r)/r = \epsilon_\mathrm{Schwarz}$ 所需的半径（PR #4945 给出的示例输出为 `the requested cutoff radius 10.000 A is smaller than what is necessary to satisfy erfc(omega*r)/r = EPS_SCHWARZ at r = 14.279 A`）。这条警告是在提醒你：截断半径成了一个需要做收敛测试的参数。

本例 64 原子 Si 的 L = 10.862 Å，L/2 = 5.43 Å，取 `RC 5.0`。截断在 5 Å 时 erfc 尚未完全衰减，因此：

- 正式计算应比较 64 原子与 216 原子超胞（或不同 $R_c$）下的目标性质；
- 小原胞、需要致密 k 网格时，改用第 8 节的 RI-HFXk；
- 2026.1 及以前的版本中，本输入的 `CUTOFF_RADIUS` 不起作用，需要用更大的超胞或 RI-HFXk。

### 7.3 运行

```bash
sed 's/^@SET STAGE .*/@SET STAGE PBE/' si_hse06_admm.inp > si64_pbe.inp
mpirun -np 32 cp2k.psmp -i si64_pbe.inp -o si64_pbe.out            # writes si64_PBE-RESTART.wfn
export OMP_NUM_THREADS=2
mpirun -np 16 cp2k.psmp -i si_hse06_admm.inp -o si64_hse06.out
grep -i -A3 "cutoff radius" si64_hse06.out | head -12
grep "Hartree-Fock Exchange energy:" si64_hse06.out | tail -1
```

> [!NOTE]
> Γ 点超胞的 HOMO–LUMO 差并不等于 Si 的间接带隙（导带底在 Γ–X 线上约 0.85 X 处，不一定正好折叠到超胞的 Γ 点）。带隙要通过 DOS 或能带计算获得，见[第10章](10-dos-pdos-band.md)。作为参照，Si 的实验带隙约 1.17 eV（低温），HSE06 通常给出接近实验的值，PBE 则明显偏低。

---

## 8. k 点杂化：RI-HFXk（CP2K ≥ 2024.1）

四中心 ERI 方法只能做 Γ 点。小原胞配致密 k 网格时，CP2K 提供基于局域原子 RI 基组的 RI-HFXk（Bussy & Hutter, 2024），输入手册中 `&HF/&RI` 的 `KP_*` 关键字从 2024.1 开始出现。官方文档说明：其主要开销是构建实空间交换矩阵，因此成本基本与 k 点数无关，k 点很多（约 1000）时才因每个 k 点的对角化而线性增长；ADMM 可以直接配合使用。

以下片段展示 Si 原胞（2 原子）PBE0 + RI-HFXk 的关键设置（思路示意，未在本仓库实测；官方文档的例子使用全电子 pob 基组和 GAPW，赝势 + GPW 的组合请先在小体系上验证，并用 `cp2k.psmp -c` 检查语法）：

```text
&DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
  BASIS_SET_FILE_NAME BASIS_ADMM_UZH
  POTENTIAL_FILE_NAME POTENTIAL_UZH
  SORT_BASIS EXP                 ! sort basis functions by exponent (data locality)
  AUTO_BASIS RI_HFX MEDIUM       ! generate the RI basis automatically
  WFN_RESTART_FILE_NAME si_pbe-RESTART.kp
  &AUXILIARY_DENSITY_MATRIX_METHOD
    ADMM_TYPE ADMMS
  &END AUXILIARY_DENSITY_MATRIX_METHOD
  &QS
    EPS_PGF_ORB 1.0E-5           ! keep equal to HF/RI/EPS_PGF_ORB
  &END QS
  &KPOINTS
    SCHEME MONKHORST-PACK 6 6 6
  &END KPOINTS
  &SCF
    SCF_GUESS RESTART
    EPS_DIIS 0.05
    &DIAGONALIZATION
      ALGORITHM STANDARD
    &END DIAGONALIZATION
  &END SCF
  &XC
    &XC_FUNCTIONAL
      &PBE
        SCALE_X 0.75
      &END PBE
    &END XC_FUNCTIONAL
    &HF
      FRACTION 0.25
      &RI
        KP_NGROUPS 8             ! MPI sub-groups; total ranks must be divisible by it
        EPS_PGF_ORB 1.0E-5
        &PRINT
          KP_RI_MEMORY_ESTIMATE
        &END PRINT
      &END RI
      &INTERACTION_POTENTIAL
        POTENTIAL_TYPE TRUNCATED
        CUTOFF_RADIUS 5.0
      &END INTERACTION_POTENTIAL
    &END HF
  &END XC
&END DFT
```

要点（均来自官方 RI-HFXk 文档）：

- k 点计算必须用对角化；
- 截断半径的约束换成“半径为 $R_c$ 的球要放得进 Born–von Kármán 超胞”：本例 FCC 原胞晶面间距约 3.14 Å，6×6×6 网格的 BvK 超胞半宽约 9.4 Å，$R_c$ = 5 Å 满足要求。文档指出多数体系 $R_c$ 从约 6 Å 起就已收敛，但应自行测试；
- `KP_NGROUPS` 是最重要的性能参数（有效到约 16 组，32 组仍合理），代价是内存复制，总 MPI 进程数必须能被它整除；
- RI 度量（`RI_METRIC`）保持默认（与 HFX 算符相同）即可；
- 从 PBE 的 `.kp` 重启可减少 SCF 步数。由于 HFX 需要比 PBE 多得多的周期像，文档建议 PBE 预计算时也用较紧的 `EPS_PGF_ORB`；
- 2026.2 及以前，k 点计算不能用 OT。

---

## 9. CP2K 2026.2 新功能：ACE 与 libGint

### 9.1 ACE（Adaptively Compressed Exchange）

ACE（Lin, 2016）在一次完整的 HFX 构建之后，用低秩投影算符 $K_\mathrm{ACE} = -\mathbf W\mathbf W^{\mathsf T}$ 代替完整交换矩阵，之后若干 SCF 步只需要这个投影，从而减少 HFX 构建次数。2026.2 把它接入了 HFX/ADMM 的基态计算路径（PR #5238，含力计算）。

```text
&SCF
  &DIAGONALIZATION              ! ACE requires diagonalization-based SCF (not OT)
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
    ALPHA 0.3
  &END MIXING
&END SCF
...
&HF
  FRACTION 0.25
  &ACE
    ACTIVE .TRUE.
    REBUILD_FREQUENCY 20        ! rebuild projectors every N SCF steps (default 20)
  &END ACE
  ...
&END HF
```

| 关键字 | 默认值 | 说明 |
|---|---|---|
| `ACTIVE` | F | 打开 ACE |
| `REBUILD_FREQUENCY` | 20 | 每 N 个 SCF 步重建一次投影算符；N = 1 等价于完整 HFX（仅用于测试）；手册建议生产计算 10–20。结构变化时总会重建 |

> [!IMPORTANT]
> 手册明确说明 ACE 需要基于对角化的 SCF，与 OT **不兼容**（OT 不构造 MO 系数）。本章用 OT 的示例不能直接打开 ACE，需要先把 SCF 改成 `&DIAGONALIZATION` + `&MIXING`。ACE 是 2026.2 的新功能，建议先在小体系上与不开 ACE 的结果对比。

### 9.2 libGint（CUDA HFX）

2026.2 新增 `&HF/HFX_LIBRARY`，可选 `LIBINT`（默认，CPU）、`LIBGINT`（CUDA GPU 加速 HFX）和 `BOTH`（临时调试选项，两者都算并比较 Fock 矩阵）。使用 `LIBGINT` 需要编译时接入 libGint（见[第1章](01-install.md)的 GPU 构建说明）。

---

## 10. 内存、并行与性能

### 10.1 估算 `MAX_MEMORY`

1. 先跑一次（或用小体系）看输出中的 `HFX_MEM_INFO` 块：`Est. max. program size before HFX [MiB]` 表示 HFX 之前程序已占用的内存，`Total memory consumption ERI's RAM [MiB]` 是 ERI 实际用量，`Number of sph. ERI's calculated on the fly` 是放不下而需要每步重算的积分数；
2. 每个 MPI 进程可用内存 ≈（节点内存 × 0.8 − 系统占用）/ 每节点 MPI 进程数，再减去 `Est. max. program size before HFX`，即可作为 `MAX_MEMORY` 的上限；
3. 内存不够时，优先**减少每节点 MPI 进程、增加 OpenMP 线程**：`MAX_MEMORY` 是按 MPI 进程计的，同一进程内的线程共享这块内存；其次增加节点数（ERI 分布到更多进程上）；最后才考虑放宽 `EPS_SCHWARZ`。

示例（每节点 128 核、512 GB）：

```bash
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
export OMP_NUM_THREADS=4
# per MPI rank: about 512*0.8/32 = 12.8 GB -> MAX_MEMORY ~ 10000 (MiB), minus the pre-HFX size
```

### 10.2 降低开销的手段（按推荐顺序）

| 手段 | 效果 | 代价 / 注意 |
|---|---|---|
| ADMM | HFX 加速约一个量级 | 需检查辅助基组误差 |
| 先 PBE 再重启 + `SCREEN_ON_INITIAL_P T` | SCF 步数和积分数都减少 | 需要重启文件 |
| 合理的 `CUTOFF_RADIUS` | 截断越小积分越少 | 必须 ≤ L/2，并做收敛测试 |
| `EPS_SCHWARZ` 1E-6 | 比默认 1E-10 少很多积分 | 不要再放宽；力精度看 `EPS_SCHWARZ_FORCES` |
| `EPS_STORAGE_SCALING 0.1` | 压缩存储，节省内存 | 精度影响小 |
| MPI × OpenMP 混合 | 每进程可用内存增大 | 线程数一般 2–8 |
| ACE（2026.2） | 减少 HFX 构建次数 | 需对角化，新功能 |
| RI-HFXk（k 点） | 小原胞 + 致密 k 网格时远比超胞高效 | 需对角化，参数较多 |

### 10.3 在 MD 和优化中使用杂化泛函

- `SCREEN_ON_INITIAL_P T` 在 MD 第一步需要重启文件提供密度矩阵；
- 结构优化中每一步都要重新计算 ERI（原子移动后积分变了），开销与单点相近；
- AIMD 常用 ADMM + 截断 PBE0 + D3，并把杂化泛函 AIMD 的步数控制在必要范围内，先用 GGA 平衡，见第12章。

---

## 参数选择与收敛建议

| 参数 | 推荐值 | 收敛/检查方法 |
|---|---|---|
| `FRACTION` | PBE0/HSE06 0.25，B3LYP 0.20 | 不要依赖默认值 1.0 |
| `POTENTIAL_TYPE` | 分子 `COULOMB`；周期 `TRUNCATED`；HSE `SHORTRANGE` | — |
| `CUTOFF_RADIUS`（Å） | ≤ L/2（Γ 点）或 ≤ BvK 超胞半宽（k 点），常见 4–6 Å | 改变 $R_c$ ±1 Å 比较目标性质；配 `PBE_HOLE_T_C_LR` 减弱依赖 |
| `EPS_SCHWARZ` | 1E-6（周期）、1E-8（小分子） | 与 1E-7/1E-8 对比能量 |
| `EPS_SCHWARZ_FORCES` | 约等于期望的力精度（默认 1E-6） | 优化/MD 中检查力 |
| `MAX_MEMORY`（MiB/进程） | 尽量让 on-the-fly ERI 为 0 | 看 `HFX_MEM_INFO` |
| ADMM 辅助基组 | `admm-dzp-qN` 起步 | 与 `admm-tzp`、无 ADMM 对比 |
| `ADMM_TYPE` | `ADMM2` / `ADMMS`，或显式 `BASIS_PROJECTION` + `MO_DIAG` | 同一项目内保持一致 |
| `EXCH_CORRECTION_FUNC` | `DEFAULT` | 截断/短程算符时尤其重要 |
| 基组 / 赝势 | UZH `-PBE0-`（`-HYB-`）系列 + `GTH-PBE0-qN` | `-qN` 一致 |
| 初猜 | PBE 收敛后 `SCF_GUESS RESTART` | 主基组、赝势、原子顺序保持一致 |

---

## 常见问题与报错

**1. `Periodic Hartree Fock calculation requested with the use of a truncated or shortrange potential. The cutoff radius (x A) is larger than half the minimal cell dimension (y A). This may lead to unphysical total energies.`**
原因：$R_c$ 超过了晶胞最小晶面间距的一半（k 点计算为 BvK 超胞，对应的警告以 `K-point Hartree-Fock calculation requested ...` 开头）。HSE06 在 2026.1 及以前会自动取一个很大的 $R_c$，因此特别常见。
解决：减小 `CUTOFF_RADIUS`（截断库仑；2026.2 起也可用于 `SHORTRANGE`），或扩大超胞 / 增加 k 点。该信息是警告而非中止，但不处理会导致能量不可靠。

**2. 杂化计算内存溢出（`Out of memory`、作业被 OOM killer 杀掉）或第一步极慢**
原因：ERI 在第一步计算并按 `MAX_MEMORY`（每 MPI 进程）存储，进程太多 × `MAX_MEMORY` 太大超过了节点内存；或 `MAX_MEMORY` 太小导致大量积分 on-the-fly 重算。
解决：按 10.1 节重新估算；减少每节点 MPI 进程、增加 OpenMP 线程；用 ADMM；从 PBE 重启并打开 `SCREEN_ON_INITIAL_P T`；不要把 `EPS_SCHWARZ` 设得过松来“省内存”。

**3. 打开 `SCREEN_ON_INITIAL_P T` 后能量异常或 SCF 不收敛**
原因：没有可用的重启文件（`SCF_GUESS ATOMIC`），初始密度矩阵质量差，筛掉了不该筛的积分。
解决：只在 `SCF_GUESS RESTART` 且重启文件来自同一体系、同一基组时打开。

**4. 结果看起来像纯 GGA，`Hartree-Fock Exchange energy` 为 0 或不出现**
原因：用 LibXC 杂化泛函（如 `&HYB_GGA_XC_HSE06`）但没写 `&HF` 段。
解决：补上 `&HF`（`FRACTION`、`&INTERACTION_POTENTIAL`）。

**5. 结果像纯 HF（带隙极大）**
原因：写了 `&HF` 但忘了 `FRACTION`，默认值是 1.0。
解决：显式写 `FRACTION 0.25`，并检查 GGA 交换的缩放（PBE0 中 `SCALE_X 0.75`）。

**6. 2026.2 的输入拿到 2026.1 及以前的版本运行，报 unknown keyword `HFX_LIBRARY` 或 unknown subsection `ACE`**
原因：这两项是 2026.2 新增的。
解决：旧版本中删除；跨版本输入用 `@IF` 包起来，运行前 `cp2k.psmp -c` 检查。

**7. `The requested basis set <admm-dzp-q6> for element <O> was not found ...`**
原因：没有把 `BASIS_ADMM_UZH` 加到 `BASIS_SET_FILE_NAME`，或辅助基组的 `-qN` 与赝势不一致。
解决：再加一行 `BASIS_SET_FILE_NAME BASIS_ADMM_UZH`；核对 `-qN`。

**8. PBE → PBE0 重启后 SCF 反而更难收敛**
原因：两步使用了不同的主基组或赝势（如 PBE 用 `-PBE-`、PBE0 用 `-PBE0-`），读入的波函数与当前基组不匹配；或原子顺序变了。
解决：两步使用相同的主基组、赝势和坐标（本章示例都这样做）。

**9. ACE 打开后报错或无效**
原因：SCF 使用了 OT。
解决：改用 `&DIAGONALIZATION` + `&MIXING`。

---

## 练习

1. 对 `h2o_pbe0_admm.inp` 分别用 `admm-dz`、`admm-dzp`、`admm-tzp` 和不用 ADMM（删除 `&AUXILIARY_DENSITY_MATRIX_METHOD` 和 `AUX_FIT` 行）计算单点能量和耗时，整理 ADMM 误差与加速比。
2. 在 32 水盒子上，把 `RC` 分别设为 3.0、3.5、4.0、4.5 Å，比较带和不带 `PBE_HOLE_T_C_LR` 时总能量对 $R_c$ 的依赖。
3. 把 water32 的 PBE0 计算分别用 `SCF_GUESS ATOMIC`（`SCREEN_ON_INITIAL_P F`）和 PBE 重启（`SCREEN_ON_INITIAL_P T`）运行，比较 SCF 步数、`HFX_MEM_INFO` 中的 ERI 数量和总耗时。
4. 固定总核数，分别用 64×1、32×2、16×4 的 MPI×OpenMP 组合运行 water32 PBE0，记录 `Number of sph. ERI's calculated on the fly` 和总时间，找出最优组合。
5. （2026.2）对 Si64 分别取 `RC` = 4.0、4.5、5.0、5.4 Å 计算 HSE06 总能量，并与 3×3×3（216 原子）超胞、`RC` = 7 Å 的结果比较。
6. （2026.2）把 h2o 示例的 SCF 改成对角化，打开 `&ACE`，比较 `REBUILD_FREQUENCY` 为 1、10、20 时的能量和耗时。

---

## 参考资料

**官方文档**

- Hartree–Fock Exchange（总览）：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/index.html>
- HFX with ADMM：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/admm.html>
- HFX-RI for Γ-Point (non-periodic)：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/ri_gamma.html>
- HFX-RI with k-Points（RI-HFXk）：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/ri_kpoints.html>
- `&XC/&HF` 输入参考（2026.2，含 `ACE`、`HFX_LIBRARY`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC/HF.html>
- `&HF/&INTERACTION_POTENTIAL`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC/HF/INTERACTION_POTENTIAL.html>
- `&AUXILIARY_DENSITY_MATRIX_METHOD`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/AUXILIARY_DENSITY_MATRIX_METHOD.html>
- Pseudopotentials（UZH `-PBE0-`/`-HYB-` 系列）：<https://manual.cp2k.org/trunk/methods/dft/pseudopotentials.html>
- CP2K v2026.2 release notes（ACE、HFX 短程截断半径、libGint）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- PR #4945（`SHORTRANGE` 下的 `CUTOFF_RADIUS`）：<https://github.com/cp2k/cp2k/pull/4945>
- Sobereva：CP2K 杂化泛函计算要点（中文）：<http://sobereva.com/690>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-hfx-*`、`regtest-admm-*`、`regtest-hybrid-ace-gs` 和 `benchmarks/QS_reference` 中有杂化泛函与 ADMM 的标准输入
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库，含 RI-HFX 专题
- [compchem-cybertraining/Tutorials_CP2K](https://github.com/compchem-cybertraining/Tutorials_CP2K)：含钙钛矿的杂化泛函与大体系性能示例（基于 v6.1/7.1，部分关键字已过时）
- [csc-training/sscc-cp2k](https://github.com/csc-training/sscc-cp2k)：PBE 下做 NEB，再用 PBE0 重新收敛 MEP 的完整案例
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含 PBE0 章节
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含杂化泛函计算经验
- [juerghutter/BASIS](https://github.com/juerghutter/BASIS)：Jürg Hutter 维护的基组集合，可了解 MOLOPT/UZH 与 ADMM 基组来源
- [DCM-Uni-Paderborn/UZH-protocol](https://github.com/DCM-Uni-Paderborn/UZH-protocol)：UZH 基组与赝势的验证数据

**本仓库脚本**：`../scripts/cp2k_outparse.py`（提取能量、警告汇总、耗时）、`../scripts/job_cp2k_slurm.sh`（MPI × OpenMP 作业模板）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- C. Adamo, V. Barone, *J. Chem. Phys.* **110**, 6158 (1999)，PBE0：<https://doi.org/10.1063/1.478522>
- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* **118**, 8207 (2003)，HSE：<https://doi.org/10.1063/1.1564060>
- A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, G. E. Scuseria, *J. Chem. Phys.* **125**, 224106 (2006)，HSE06 参数：<https://doi.org/10.1063/1.2404663>
- M. Guidon, F. Schiffmann, J. Hutter, J. VandeVondele, *J. Chem. Phys.* **128**, 214104 (2008)，杂化泛函 AIMD：<https://doi.org/10.1063/1.2931945>
- M. Guidon, J. Hutter, J. VandeVondele, *J. Chem. Theory Comput.* **5**, 3010 (2009)，周期 HFX 与截断库仑：<https://doi.org/10.1021/ct900494g>
- M. Guidon, J. Hutter, J. VandeVondele, *J. Chem. Theory Comput.* **6**, 2348 (2010)，ADMM：<https://doi.org/10.1021/ct1002225>
- P. Merlot *et al.*, *J. Chem. Phys.* **141**, 094104 (2014)，电荷约束 ADMM（ADMMS/P/Q）：<https://doi.org/10.1063/1.4894267>
- L. Lin, *J. Chem. Theory Comput.* **12**, 2242 (2016)，ACE：<https://doi.org/10.1021/acs.jctc.6b00092>
- A. Bussy, J. Hutter, *J. Chem. Phys.* **160**, 064116 (2024)，RI-HFXk：<https://doi.org/10.1063/5.0189659>
- M. Iannuzzi *et al.*, The CP2K Program Package Made Simple, *J. Phys. Chem. B* **130**, 1237 (2026)：<https://doi.org/10.1021/acs.jpcb.5c05851>

---

[← 上一章](08-dispersion-spin-dftu.md) | [返回目录](../README.md) | [下一章 →](10-dos-pdos-band.md)
