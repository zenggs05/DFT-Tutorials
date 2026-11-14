# 第13章 过渡态搜索：NEB 与 Dimer（Transition-State Search in CP2K: CI-NEB & Dimer）

本章讲解用 CP2K 搜索过渡态（transition state, TS）、计算反应能垒的完整流程：先严格优化初末态，再用 `RUN_TYPE BAND` 做 NEB/CI-NEB（climbing-image nudged elastic band）得到最小能量路径（MEP），然后用 `GEO_OPT TYPE TRANSITION_STATE` 的 Dimer 方法精修鞍点，最后用频率计算确认结构只有一个虚频。全章以 H 原子在 Cu(111) 表面 fcc → hcp 空位之间的扩散为例，同时说明 `NUMBER_OF_REPLICA × NPROC_REP` 的并行设置、各副本的 `-BAND*.out` 日志，以及 CP2K 2026.2 新增的 NEB 能量剖面图、`FINAL` 结构文件和从 molden 文件初始化 Dimer 的功能。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 NEB、IT-NEB、CI-NEB 与 Dimer 的原理和各自的适用条件；
- 准备满足 NEB 要求的初态/末态：同一个 `&FORCE_EVAL`、同样的固定原子、相同的原子顺序，并且都已严格优化；
- 写出 `RUN_TYPE BAND` 的完整输入，理解 `BAND_TYPE`、`NUMBER_OF_REPLICA`、`K_SPRING`、`NPROC_REP`、`ALIGN_FRAMES`/`ROTATE_FRAMES`、`&OPTIMIZE_BAND`、`&CONVERGENCE_CONTROL` 等关键字；
- 按 `NUMBER_OF_REPLICA × NPROC_REP` 分配 MPI 进程，知道去哪里找各副本的日志和报错；
- 读懂 BAND 输出中的距离与能量表，判断 NEB 是否收敛，提取能垒和过渡态结构；
- 用 Dimer 方法从一个近似 TS 出发收敛到鞍点（2026.2 起可以用振动模式初始化 dimer 方向）；
- 用频率计算验证过渡态（只有一个虚频），并给出零点能校正后的能垒和谐振子过渡态理论速率的估算方法。

| 方法 | 需要什么 | 得到什么 | CP2K 写法 |
|---|---|---|---|
| NEB / IT-NEB | 初态 + 末态（可加中间猜测） | 离散的最小能量路径，能垒近似值 | `RUN_TYPE BAND`，`BAND_TYPE IT-NEB` |
| CI-NEB | 同上 | 最高点被推到鞍点，能垒更准确 | `BAND_TYPE CI-NEB` + `&CI_NEB` |
| Dimer | 一个近似 TS（来自 NEB、化学直觉或扫描） | 精确鞍点，不需要末态 | `RUN_TYPE GEO_OPT` + `TYPE TRANSITION_STATE` + `METHOD DIMER` |
| 频率验证 | 收敛的 TS 结构 | 虚频个数、ZPE、指前因子 | `RUN_TYPE VIBRATIONAL_ANALYSIS`，见 [第14章](14-vibration-thermo.md) |
| 自由能能垒（有限温度、溶剂） | 合适的集体变量 | 自由能面 | 元动力学 / 约束 MD，见 [第15章](15-metadynamics-plumed.md) |
| 酶或溶液中的大体系反应 | QM/MM 拓扑 | QM/MM 下的 MEP 与 TS | 见 [第16章](16-qmmm-xtb.md) |

本章的完整输入在 [../examples/13-neb-dimer-h-cu111/README.md](../examples/13-neb-dimer-h-cu111/README.md)。流程如下：

```text
initial.xyz ─┐                                            ┌─> ts_guess_vib.inp (只动 H 的部分 Hessian)
             ├─ geo_opt_endpoint.inp ─> h_diffusion_cineb.inp ─┤        │ molden 文件
final.xyz  ──┘   (端点严格优化)          (CI-NEB, 7 个副本)     └─> dimer_ts.inp (Dimer 精修)
                                                                         │
                                                第14章 ts_freq_check.inp <┘ (确认只有 1 个虚频)
```

> [!NOTE]
> 示例体系是 3×3×3 的 Cu(111) slab（27 个 Cu + 1 个 H），只算 Γ 点，晶格常数取约等于 PBE 值的 3.63 Å。这是为教学压缩过的设置：正式研究请用你自己 `CELL_OPT` 得到的晶格常数（[第7章](07-geo-cell-opt.md)），把表面 k 点收敛到 3×3×1 或更密（[第6章](06-periodic-surface-kpoints.md)），并检查 slab 层数和真空厚度。

---

## 2. 原理速览

### 2.1 鞍点与最小能量路径

在势能面 $E(\mathbf R)$ 上，反应物和产物是两个局部极小点。连接它们、在每一点上与梯度方向平行的路径称为**最小能量路径**（MEP）。MEP 上能量最高的点是一阶鞍点，也就是过渡态，它满足：

1. 梯度为零：$\nabla E(\mathbf R^\ddagger) = 0$；
2. Hessian $\mathbf H = \partial^2 E/\partial \mathbf R^2$（质量加权后）恰有**一个**负本征值，对应的本征矢就是反应坐标方向。

振动分析中，负本征值表现为一个虚频（CP2K 把它打印成负数），其余 $3N-7$（非线性分子）或去掉固定原子后的全部模式都应为实频。

能垒定义为 $E_a = E(\mathbf R^\ddagger) - E(\mathbf R_\mathrm{IS})$。在谐振子过渡态理论（hTST）中，速率常数为

$$
k = \nu^{*}\exp\!\left(-\frac{E_a}{k_\mathrm B T}\right),\qquad
\nu^{*} = \frac{\prod_{i=1}^{3N}\nu_i^\mathrm{IS}}{\prod_{j=1}^{3N-1}\nu_j^{\ddagger}}
$$

其中 $\nu^{*}$ 是 Vineyard 指前因子：分子是初态全部实频的乘积，分母是 TS 除虚频外所有实频的乘积。表面扩散的 $\nu^{*}$ 通常在 $10^{12}$–$10^{13}\ \mathrm{s^{-1}}$ 量级（经验范围，以你的计算为准）。

### 2.2 NEB：弹簧串起的一组副本

NEB 把路径离散成 $M$ 个副本（replica / image）$\mathbf R_1,\dots,\mathbf R_M$，其中 $\mathbf R_1$、$\mathbf R_M$ 是固定的端点。对中间副本 $i$，先估计路径切向 $\hat\tau_i$，再构造 NEB 力：

$$
\mathbf F_i^\mathrm{NEB} = -\nabla E(\mathbf R_i)\big|_\perp + \mathbf F_i^{s}\big|_\parallel,\qquad
\mathbf F_i^{s}\big|_\parallel = k\left(|\mathbf R_{i+1}-\mathbf R_i| - |\mathbf R_i-\mathbf R_{i-1}|\right)\hat\tau_i
$$

“nudged”的意思是：真实力只保留垂直于路径的分量（让副本落到 MEP 上），弹簧力只保留沿路径的分量（让副本沿路径均匀分布）。这样弹簧不会把路径“拉直”而抄近道，真实力也不会让副本滑向两端的极小点。弹簧常数 $k$ 就是 CP2K 的 `K_SPRING`。

切向的取法决定了 NEB 的变体。CP2K 的 `BAND_TYPE` 支持：

| `BAND_TYPE` | 含义 | 说明 |
|---|---|---|
| `IT-NEB`（默认） | improved tangent NEB | Henkelman & Jónsson (2000) 的改进切向，避免路径出现扭结，是日常首选 |
| `CI-NEB` | climbing image NEB | 在 IT-NEB 基础上让能量最高的副本“爬坡”到鞍点 |
| `B-NEB` | bisection NEB | 原始（二分）切向 |
| `D-NEB` | doubly nudged elastic band | 保留部分垂直弹簧力，常用于 L-BFGS 类优化器 |
| `EB` | elastic band（哈密顿形式） | 不做 nudging 的弹性带 |
| `SM` | string method | 字符串方法，靠重参数化而不是弹簧保持等距 |

### 2.3 CI-NEB：让最高点爬到鞍点

普通 NEB 中，能量最高的副本只是“离鞍点最近的那个离散点”，能垒会被低估，低估的程度取决于副本密度。CI-NEB（Henkelman, Uberuaga & Jónsson, 2000）对能量最高的副本 $i_\mathrm{max}$ 去掉弹簧力，并把真实力沿切向的分量反号：

$$
\mathbf F_{i_\mathrm{max}} = -\nabla E(\mathbf R_{i_\mathrm{max}}) + 2\,\big(\nabla E(\mathbf R_{i_\mathrm{max}})\cdot\hat\tau\big)\,\hat\tau
$$

于是它在垂直方向上往下走，在切向上往上爬，收敛后正好停在鞍点上。CP2K 的 `&CI_NEB/NSTEPS_IT`（默认 5）规定先做多少步普通 IT-NEB，再打开爬坡。路径还很粗糙时过早爬坡，容易选错最高点，所以初始路径较差时可以把这个值加大。

### 2.4 Dimer：只用一阶导数找鞍点

Dimer 方法（Henkelman & Jónsson, 1999）用两个相距很近的点 $\mathbf R_{1,2} = \mathbf R_0 \pm \Delta R\,\hat N$ 组成一个“哑铃”。沿 $\hat N$ 方向的曲率可以只用两端的力估算：

$$
C = \frac{(\mathbf F_2-\mathbf F_1)\cdot\hat N}{2\,\Delta R}
$$

每个迭代分两步：

1. **转动**：在固定中心 $\mathbf R_0$ 的前提下旋转 $\hat N$，使曲率 $C$ 最小，$\hat N$ 就逼近最低曲率模式（在 TS 附近就是虚频模式）；
2. **平移**：用修正力 $\mathbf F^\dagger = \mathbf F_0 - 2(\mathbf F_0\cdot\hat N)\hat N$ 移动中心。这相当于沿 $\hat N$ 往上爬、在其他方向往下走。

Dimer 不需要末态，也不需要 Hessian，每步只要两三次力计算。它的成败取决于**初始结构是否已在鞍点附近**，以及**初始 dimer 方向是否合理**。CP2K 默认随机生成初始方向，2026.2 起可以从振动分析的 molden 文件中读入某个（或某几个）振动模式作为初始方向（`INITIALIZATION_METHOD MOLDEN`），这正是本章 §8 的做法。CP2K 还实现了约束平移的 K-Dimer 变体（`K-DIMER T`，参数 `BETA`，文献见输入参考）。

---

## 3. 准备端点：同一套设置下的严格优化

NEB 对端点有三条硬性要求：

1. **同一个 `&FORCE_EVAL`**：基组、赝势、CUTOFF、泛函、色散、SCF 设置、晶胞都必须相同，否则端点能量和中间副本能量无法比较；
2. **同样的约束**：几何优化中固定了哪些原子，NEB 里也要固定同样的原子；
3. **相同的原子顺序和一致的周期映像**：第 $k$ 个原子在初态和末态必须是“同一个原子”。CP2K 按坐标逐原子线性插值，原子顺序错了，插出来的中间结构就是原子互相穿越的荒谬构型。

本例用一个包含文件 `cu111_h_force_eval.inc` 把 `&FORCE_EVAL` 共享给所有步骤（`@INCLUDE` 的用法见 [第2章](02-input-structure.md)），从机制上保证第 1 条：

```fortran
! cu111_h_force_eval.inc - shared FORCE_EVAL for H on a Cu(111) 3x3x3 slab
! The including file must define COORD; EPS_SCF is optional (default 1.0E-7).
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    SURFACE_DIPOLE_CORRECTION T
    SURF_DIP_DIR Z
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF ${EPS_SCF-1.0E-7}
      MAX_SCF 200
      ADDED_MOS 60
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD BROYDEN_MIXING
        ALPHA 0.2
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
      &PRINT
        &RESTART
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 7.70039 7.70039 22.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Cu
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

几点说明：

- Cu 是金属，所以用对角化 + Broyden 混合 + Fermi–Dirac 展宽，`ADDED_MOS 60` 覆盖费米面以上的分数占据（298 个价电子，149 个占据轨道）。这些设置的原理见 [第5章](05-scf-strategy.md)。OT 在 2026.2 及以前不支持展宽，金属 NEB 不要用 OT。
- `EPS_SCF` 默认 1.0E-7，比单点计算更紧：NEB 和 Dimer 的收敛判据是力，SCF 噪声会直接变成力噪声。写成 `${EPS_SCF-1.0E-7}` 是为了让[第14章](14-vibration-thermo.md)的频率验证输入在包含同一文件时把它收紧到 1.0E-8。
- 只有一侧吸附 H，slab 两侧不对称，所以打开偶极校正（[第6章](06-periodic-surface-kpoints.md)）。
- 坐标文件名用变量 `${COORD}`，由主输入在 `@INCLUDE` 之前 `@SET`。

端点优化的主输入 `geo_opt_endpoint.inp`：

```fortran
@SET COORD ${COORD-initial.xyz}
@SET PROJ ${PROJ-opt_initial}

&GLOBAL
  PROJECT ${PROJ}
  RUN_TYPE GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

@INCLUDE 'cu111_h_force_eval.inc'

&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 300
    MAX_FORCE 2.0E-4
    RMS_FORCE 1.0E-4
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..9
    &END FIXED_ATOMS
  &END CONSTRAINT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
    &END RESTART
  &END PRINT
&END MOTION
```

`${COORD-initial.xyz}` 表示“变量 COORD 未定义时取默认值 initial.xyz”，因此同一个输入可以用命令行 `-E` 切换到末态：

```bash
cd CP2K/examples/13-neb-dimer-h-cu111
mpirun -np 32 cp2k.psmp -i geo_opt_endpoint.inp -o opt_initial.out
mpirun -np 32 cp2k.psmp -E COORD=final.xyz -E PROJ=opt_final -i geo_opt_endpoint.inp -o opt_final.out

grep "GEOMETRY OPTIMIZATION COMPLETED" opt_initial.out opt_final.out
grep "ENERGY| Total FORCE_EVAL" opt_initial.out | tail -1
grep "ENERGY| Total FORCE_EVAL" opt_final.out   | tail -1

# 取最后一帧作为 NEB 端点
bash last_frame.sh opt_initial-pos-1.xyz initial_opt.xyz
bash last_frame.sh opt_final-pos-1.xyz   final_opt.xyz
```

`MAX_FORCE 2.0E-4` Ha/bohr（约 0.01 eV/Å）比默认的 `4.5E-4` 更紧。端点优化得不充分，NEB 两端会残留力，能垒就有系统误差。2026.2 起，优化结束时还会写出 `<PROJECT>-FINAL-1_<步数>.xyz`（扩展 XYZ）和同名 `.cif`，也可以直接拿来用（[第7章](07-geo-cell-opt.md)）。

**检查原子映射。** 下面这段 Python 逐个原子比较两个端点的位移。正常情况下只有 H 和少数表面原子明显移动；如果大量原子位移都在 1 Å 以上，说明原子顺序或周期映像不一致：

```python
import numpy as np

def read_xyz(fn):
    lines = open(fn).read().split("\n")
    n = int(lines[0])
    sym = [l.split()[0] for l in lines[2:2 + n]]
    xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines[2:2 + n]])
    return sym, xyz

s1, a = read_xyz("initial_opt.xyz")
s2, b = read_xyz("final_opt.xyz")
assert s1 == s2, "element order differs!"
d = np.linalg.norm(b - a, axis=1)
for i in np.argsort(d)[::-1][:5]:
    print(f"atom {i + 1:3d} {s1[i]:2s} moved {d[i]:.3f} A")
```

> [!WARNING]
> CP2K 在 NEB 中直接用笛卡尔坐标插值，不会自动替你选最近的周期映像。如果某个原子在一个端点中被“卷回”到晶胞另一侧（坐标差接近一个晶格矢量），插值路径会穿过整个晶胞。先在可视化软件里叠加两个端点检查，必要时手动把坐标平移一个晶格矢量。

---

## 4. CI-NEB 输入逐段讲解

完整输入 `h_diffusion_cineb.inp`（`&FORCE_EVAL` 来自上一节的包含文件）：

```fortran
@SET COORD ${INIT-initial_opt.xyz}

&GLOBAL
  PROJECT h_cu111_neb
  RUN_TYPE BAND
  PRINT_LEVEL LOW
&END GLOBAL

@INCLUDE 'cu111_h_force_eval.inc'

&MOTION
  &BAND
    BAND_TYPE CI-NEB
    NUMBER_OF_REPLICA 7
    NPROC_REP ${NPROC_REP-8}
    K_SPRING 0.05
    ALIGN_FRAMES F
    ROTATE_FRAMES F
    &CI_NEB
      NSTEPS_IT 10
    &END CI_NEB
    &OPTIMIZE_BAND
      OPT_TYPE DIIS
      &DIIS
        MAX_STEPS 300
        N_DIIS 5
      &END DIIS
    &END OPTIMIZE_BAND
    &CONVERGENCE_CONTROL
      MAX_FORCE 1.0E-3
      RMS_FORCE 5.0E-4
      MAX_DR 5.0E-3
      RMS_DR 2.5E-3
    &END CONVERGENCE_CONTROL
    &REPLICA
      COORD_FILE_NAME ${INIT-initial_opt.xyz}
    &END REPLICA
    &REPLICA
      COORD_FILE_NAME ${FINAL-final_opt.xyz}
    &END REPLICA
    &PROGRAM_RUN_INFO ON
      INITIAL_CONFIGURATION_INFO
      ! CP2K >= 2026.2 only:
      ! PLOT_REL_ENERGY T
    &END PROGRAM_RUN_INFO
    &CONVERGENCE_INFO ON
    &END CONVERGENCE_INFO
  &END BAND
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..9
    &END FIXED_ATOMS
  &END CONSTRAINT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
    &END RESTART
  &END PRINT
&END MOTION
```

### 4.1 `&BAND` 主段

| 关键词 | 本例取值 | 默认值 | 说明 |
|---|---|---|---|
| `RUN_TYPE BAND`（`&GLOBAL`） | — | — | 打开 BAND 类任务（NEB 族和 string method） |
| `BAND_TYPE` | `CI-NEB` | `IT-NEB` | 见 §2.2 的表 |
| `NUMBER_OF_REPLICA` | 7 | 10 | **包括两个端点**的副本总数。可以大于 `&REPLICA` 段的个数，缺的副本会自动生成：每次在相邻副本距离最大的地方插入一个中点（线性插值），直到凑够数目 |
| `NPROC_REP` | 8 | 1 | 每个副本使用的 MPI 进程数，见 §5 |
| `K_SPRING`（别名 `K`） | 0.05 | 0.02 | 弹簧常数，原子单位。太小，副本会向两端极小点滑；太大，路径变“硬”，垂直方向收敛慢。0.02–0.1 都常见 |
| `ALIGN_FRAMES` | `F` | `T` | 计算开始时把各副本整体平移/转动对齐 |
| `ROTATE_FRAMES` | `F` | `T` | 每一步计算 RMSD 并旋转副本使其最小 |
| `OPTIMIZE_END_POINTS`（`&OPTIMIZE_BAND` 内） | 不写 | `F` | 端点是否一起优化。端点已用同一 `FORCE_EVAL` 优化过时可设 `T`，此时端点受力每步被清零 |
| `PROC_DIST_TYPE` | 不写 | `BLOCKED` | 进程到副本的映射方式（`BLOCKED`/`INTERLEAVED`） |
| `USE_COLVARS`、`POT_TYPE` | 不写 | `F`、`FULL` | 在集体变量子空间中做投影 NEB/自由能 NEB，属于高级用法 |

> [!IMPORTANT]
> 对 slab、固定了原子的体系或任何周期性晶胞，**一定要把 `ALIGN_FRAMES` 和 `ROTATE_FRAMES` 设为 `F`**。这两个选项是为孤立分子设计的：转动副本可以去掉整体转动对路径长度的贡献，但在周期晶胞中整体转动没有物理意义，还会让固定原子离开原位。孤立分子的 NEB 可以保留默认的 `T`。

### 4.2 `&CI_NEB`

只有一个关键词 `NSTEPS_IT`（默认 5）：先做多少步 IT-NEB，再打开爬坡。本例取 10，让路径先大致成形。若初始路径来自粗糙的线性插值，或者路径上可能有多个能量相近的峰，可以加大到 20–30，或者干脆先跑一个 `BAND_TYPE IT-NEB` 的任务，再用重启文件切换到 CI-NEB（§7）。

### 4.3 `&OPTIMIZE_BAND`：优化器

| 关键词 | 取值 | 说明 |
|---|---|---|
| `OPT_TYPE` | `DIIS`（默认）/ `MD` | DIIS：先做最速下降（带线搜索），再切到 DIIS 加速；MD：用阻尼分子动力学（带退火）优化 |
| `&DIIS/MAX_STEPS` | 300 | 最大优化步数，默认只有 100，DFT 级 NEB 常常不够 |
| `&DIIS/N_DIIS` | 5 | DIIS 历史向量数（默认 7）。路径振荡时减小 |
| `&DIIS/STEPSIZE` | 默认 1.0 | 线搜索初始步长，DIIS 不稳定时可减小 |
| `&DIIS/MAX_STEPSIZE` | 默认 2.0 | 线搜索最大步长，初始结构很差时减小以稳定 |
| `&DIIS/NO_LS` | 默认 `F` | 最速下降阶段不做线搜索，适合配合较小的 `STEPSIZE` 处理极不平衡的起始结构 |
| `&DIIS/MAX_SD_STEPS` | 默认 1 | 切换到 DIIS 前最多做多少步最速下降 |
| `&DIIS/CHECK_DIIS` | 默认 `T` | 对 DIIS 步做一系列检查后才接受；大量 DIIS 步被拒绝时可设 `F` |

`OPT_TYPE MD` 的写法见官方回归测试 `tests/NEB/regtest-1/2gly_CI-NEB.inp`：在 `&MD` 中设 `MAX_STEPS`、`TIMESTEP`、`TEMPERATURE`，再用 `&VEL_CONTROL` 的 `ANNEALING`（如 0.95）逐步降温。它比 DIIS 更稳，但通常需要更多步。

### 4.4 `&CONVERGENCE_CONTROL`：收敛判据

| 关键词 | 本例 | 默认值 | 单位 |
|---|---|---|---|
| `MAX_FORCE` | 1.0E-3 | 4.5E-4 | Ha/bohr（1.0E-3 Ha/bohr ≈ 0.051 eV/Å） |
| `RMS_FORCE` | 5.0E-4 | 3.0E-4 | Ha/bohr |
| `MAX_DR` | 5.0E-3 | 2.0E-4 | bohr |
| `RMS_DR` | 2.5E-3 | 1.0E-4 | bohr |

这里的力是 NEB 力（垂直真实力 + 平行弹簧力；爬坡副本是修正后的力），并对所有可动副本统计。默认的位移判据（2.0E-4/1.0E-4 bohr）对 DFT 级 NEB 偏严，噪声大时很难满足，所以本例把位移判据放宽，让力判据起主导作用。约 0.05 eV/Å 的 `MAX_FORCE` 是表面 NEB 常用的收敛标准；如果之后还要用 Dimer 精修，NEB 收敛到 0.1 eV/Å 左右再切换也可以。

### 4.5 `&REPLICA`：给出路径

每个 `&REPLICA` 段给一个副本的坐标，可以用 `COORD_FILE_NAME`（xyz 文件，单位 Å）或内嵌的 `&COORD` 段。段的顺序就是路径顺序：第一个是初态，最后一个是末态，中间可以插入任意多个猜测结构。注意：

- 每个 xyz 文件必须有两行头（原子数 + 标题行），**标题行不能为空**，否则报 `Number of lines in XYZ format not equal to the number of atoms ... Very probably the line with title is missing or is empty`；
- `&FORCE_EVAL/&SUBSYS` 里仍然要给一套坐标（这里用 `@SET COORD` 指向初态），它决定原子种类和数目，具体位置会被各副本坐标覆盖；
- 所有副本的原子数必须与 `&SUBSYS` 一致，否则程序在读入时直接断言失败。

**自动插值够不够？** CP2K 的补点是逐次二分的笛卡尔线性插值。本例的 H 在两个空位之间几乎沿直线移动，线性插值已经足够。插值完成后，如果最大与最小副本间距之比仍大于 1.5，CP2K 会在输出中提示：

```text
 NEB| After interpolating replica, the maximum distance between adjacent replica
 NEB| is still larger than 1.5 times the minimum distance. ...
```

对于有转动、键断裂重组的反应（如分子在表面解离），线性插值常常导致原子过近、SCF 难收敛。这时可以用 ASE 的 IDPP（image dependent pair potential）插值生成中间副本，再逐个写成 `&REPLICA`：

```python
from ase.io import read, write
try:
    from ase.mep import NEB          # ASE >= 3.23
except ImportError:
    from ase.neb import NEB          # older ASE

cell = [7.70039, 7.70039, 22.0, 90, 90, 60]
ini, fin = read("initial_opt.xyz"), read("final_opt.xyz")
for a in (ini, fin):
    a.set_cell(cell)
    a.pbc = True
images = [ini] + [ini.copy() for _ in range(5)] + [fin]
neb = NEB(images)
neb.interpolate("idpp", mic=True)
for i, im in enumerate(images[1:-1], start=2):
    write(f"replica_{i}.xyz", im)        # replica_2.xyz ... replica_6.xyz
```

然后在 `&BAND` 中按顺序写 7 个 `&REPLICA` 段（初态、replica_2 … replica_6、末态），`NUMBER_OF_REPLICA` 等于 7。

### 4.6 `&PROGRAM_RUN_INFO` 与 `&CONVERGENCE_INFO`

这两个打印段控制 BAND 每一步的汇总表和收敛判据表。2026.2 把它们的默认打印级别从 `MEDIUM` 降到 `LOW`，并把收敛判据表从 `-BAND*.out` 挪进主输出。为了在 2023.x–2026.2 上都能在主输出中看到汇总表，本例显式写了 `ON`。`INITIAL_CONFIGURATION_INFO` 会额外打印初始路径的构建信息。`PLOT_REL_ENERGY T`（**仅 2026.2+**）会在主输出里画一个横向的相对能量示意图，并用 `X`/`x` 标出局部极大/极小。

---

## 5. 并行：`NUMBER_OF_REPLICA × NPROC_REP`

BAND 任务把全部 MPI 进程分成若干组，每组 `NPROC_REP` 个进程，负责一个副本的能量和力。组数等于 `总进程数 / NPROC_REP`（取整）。理想情况是

$$
N_\mathrm{MPI} = \texttt{NUMBER\_OF\_REPLICA} \times \texttt{NPROC\_REP}
$$

这时所有副本（**包括两个端点**）同时计算，每一步 BAND 迭代只需一次“单点能 + 力”的墙钟时间。

| 总进程数 | `NPROC_REP` | 组数 | 每步的计算轮数 | 评价 |
|---|---|---|---|---|
| 56 | 8 | 7 | 1 | 最佳：一组一个副本 |
| 28 | 4 | 7 | 1 | 也很好，单个副本更慢但节点更少 |
| 32 | 8 | 4 | 2（4 + 3） | 第二轮有一组空闲，效率约 7/8 |
| 30 | 8 | 3 | 3（3 + 3 + 1），另有 6 个进程闲置 | 差：主输出警告 `[6] processors will be wasted!` |
| 8 | 8（或默认 1） | 1（或 8） | 7 | 能跑，但每步要串行算 7 个副本（`NPROC_REP 1` 时每个副本只有 1 个进程） |

进程数不能整除时，CP2K 会给出警告：

```text
 *** WARNING ... Number of totally requested processors (30) is not compatible with the
 number of processors requested per replica (8) and the number of replicas (3) .
 [6] processors will be wasted! ***
```

（这里的 “number of replicas” 指副本**组**数，不是 `NUMBER_OF_REPLICA`。）

实践建议：

- 先用单点计算测出“一个副本用多少进程最划算”（并行效率开始明显下降之前的进程数），把它作为 `NPROC_REP`，再乘以副本数申请资源；
- 混合并行时，`NPROC_REP` 只数 MPI 进程，OpenMP 线程另算（`OMP_NUM_THREADS`）；
- 端点在 `NUMBER_OF_REPLICA` 中占两个名额并分到进程，即使 `OPTIMIZE_END_POINTS F` 也一样（输入参考的说明），资源规划时要算上它们；
- 本例的 [job_slurm.sh](../examples/13-neb-dimer-h-cu111/job_slurm.sh) 根据 `SLURM_NTASKS` 自动设置 `-E NPROC_REP=$((NTASKS/7))`；通用的 MPI × OpenMP 模板见 [job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh) 和 [公共篇第2章](../../Common/tutorials/02-job-scheduler.md)。

---

## 6. 运行、监控与结果提取

### 6.1 运行与输出文件

```bash
mpirun -np 56 cp2k.psmp -i h_diffusion_cineb.inp -o h_diffusion_cineb.out
```

| 文件 | 内容 |
|---|---|
| `h_diffusion_cineb.out` | 主输出：BAND 设置、每步汇总表（距离、能量）、收敛判据（2026.2 起）、`BAND TASK COMPLETED`（2026.2 起） |
| `h_cu111_neb-BAND1.out` … `-BAND7.out` | 每个副本自己的 Quickstep 输出（SCF 迭代、能量、警告、报错）。副本数 ≥ 10 时编号补零，如 `BAND01` |
| `h_cu111_neb-pos-Replica_nr_1-1.xyz` … | 每个副本的结构随 BAND 步的变化（轨迹） |
| `h_cu111_neb-1.restart` | 重启文件，包含所有副本的当前坐标 |
| `h_cu111_neb-*.ener` | 由 `&BAND/&ENERGY` 控制，每行：步数、各副本能量、相邻副本距离 |
| `h_cu111_neb-FINAL-1.xyz` | **2026.2+**：结束时（无论是否收敛）所有副本的最终结构，扩展 XYZ 格式，标题行带 `Replica=`、`Energy=`、`Ener_loc_max=`、`Ener_loc_min=` |

> [!TIP]
> 主输出里只有汇总信息。某个副本 SCF 不收敛、基组找不到这类报错，往往只出现在对应的 `-BAND*.out` 里。排错时先 `grep -l "ABORT\|WARNING" h_cu111_neb-BAND*.out`。

### 6.2 读懂每一步的汇总表

每个 BAND 步在主输出中打印一个汇总块（示意输出，数值仅供参考）：

```text
 *******************************************************************************
 BAND TYPE                     =                                          CI-NEB
 BAND TYPE OPTIMIZATION        =                                            DIIS
 STEP NUMBER                   =                                              37
 NUMBER OF NEB REPLICA         =                                               7
 DISTANCES REP =        0.468712        0.471035        0.470188        0.470563
                        0.471322        0.468894
 ENERGIES [au] =    -1328.471580    -1328.469871    -1328.465904    -1328.463722
                    -1328.465911    -1328.469920    -1328.471663
 BAND TOTAL ENERGY [au]        =                           -9299.27858100000000
 *******************************************************************************
```

- `DISTANCES REP`：相邻副本之间的距离（原子单位 bohr，是全部原子坐标组成的 3N 维向量的距离），应大致相等；
- `ENERGIES [au]`：各副本能量（Hartree），第一个和最后一个是端点；
- `BAND TOTAL ENERGY`：各副本能量之和加弹簧能，只用于监控，没有物理意义。

收敛判据块（2026.2 在主输出中；2026.1 及以前在 `-BAND*.out` 中，且数字位数较少）：

```text
 *******************************************************************************
 RMS DISPLACEMENT              =      0.0006112035     [    0.0025000000 ] (YES)
 MAX DISPLACEMENT              =      0.0019876101     [    0.0050000000 ] (YES)
 RMS FORCE                     =      0.0002010442     [    0.0005000000 ] (YES)
 MAX FORCE                     =      0.0007733208     [    0.0010000000 ] (YES)
 *******************************************************************************
```

（示意输出，数值仅供参考）四项都为 `YES` 时 BAND 收敛，2026.2 会再打印一个 `BAND TASK COMPLETED` 横幅。

常用的 grep：

```bash
grep -c "STEP NUMBER" h_diffusion_cineb.out                 # 已完成的 BAND 步数
grep -A3 "ENERGIES \[au\]" h_diffusion_cineb.out | tail -4  # 最新的能量
grep "MAX FORCE" h_diffusion_cineb.out | tail -3           # 2026.2：收敛趋势
grep "MAX FORCE" h_cu111_neb-BAND*.out | tail -3           # 2026.1 及以前
grep "BAND TASK COMPLETED" h_diffusion_cineb.out           # 2026.2：收敛标志
```

### 6.3 能垒与能量曲线

示例目录中的 `neb_profile.py`（只依赖标准库，作图需要 matplotlib）读取主输出中最后一个汇总块，打印每个副本的累计反应坐标、相对能量和正/逆向能垒：

```bash
python3 neb_profile.py h_diffusion_cineb.out --png mep.png
```

本仓库还有功能更全的 [neb_plot.py](../scripts/neb_plot.py)（收集各副本能量、画 MEP、导出过渡态 xyz），用 `python ../scripts/neb_plot.py -h` 查看用法。

判断结果是否可信：

1. 收敛判据全部满足，且最后几十步能量不再系统性变化；
2. 能量最高的副本位于路径内部（不是端点或紧挨端点的副本），两端端点能量与单独 GEO_OPT 的结果一致（差值应在 SCF 精度以内）；
3. 曲线只有一个极大值。如果出现两个峰、中间夹着一个谷（2026.2 的 `PLOT_REL_ENERGY` 会标出 `x`），说明路径经过了一个**中间体**：把这个谷底副本拿出来单独优化，再把反应拆成两段分别做 NEB；
4. 副本间距大致均匀，鞍点附近没有明显的“空档”。

### 6.4 取出过渡态结构

**2026.2+**：直接从 `FINAL` 文件中取 `Ener_loc_max=T` 的那一帧（本例每帧 28 个原子，加两行头）：

```bash
grep -o 'Replica=[0-9]* Energy=[^ ]*' h_cu111_neb-FINAL-1.xyz    # 各副本编号与能量
grep -B1 -A28 "Ener_loc_max=T" h_cu111_neb-FINAL-1.xyz > ts_from_neb.xyz
```

如果有多个局部极大，上面的命令会取出多帧，请手动挑选能量最高的那个。

**2026.1 及以前**：先从汇总表中找到能量最高的副本编号（例如 4），再取它轨迹的最后一帧：

```bash
bash last_frame.sh h_cu111_neb-pos-Replica_nr_4-1.xyz ts_from_neb.xyz
```

CI-NEB 收敛后，爬坡副本本身就是过渡态的良好近似；后续可以直接做频率验证（§9），或者先用 Dimer 精修（§8）。

---

## 7. 重启与常见调整

**从中断处继续。** BAND 的重启文件 `<PROJECT>-1.restart` 记录了所有副本的坐标。在原输入中加入

```fortran
&EXT_RESTART
  RESTART_FILE_NAME h_cu111_neb-1.restart
&END EXT_RESTART
```

即可从最后保存的路径继续（官方回归测试 `tests/NEB/regtest-3/2gly_IT-NEB-res.inp` 就是这种写法）。各副本的波函数文件按副本的项目名保存（形如 `h_cu111_neb-BAND4-RESTART.wfn`），保留它们可以加快续算的第一步 SCF。

**先 IT-NEB，再 CI-NEB。** 第一个任务用 `BAND_TYPE IT-NEB` 收敛到 0.1–0.2 eV/Å，第二个任务改成 `CI-NEB` 并用 `&EXT_RESTART` 读入路径，这比一开始就爬坡更稳。

**加副本。** 想在已有路径上加密，可以把重启文件中的副本导出成 xyz，再在原来的 `&REPLICA` 列表中插入新结构，并相应增大 `NUMBER_OF_REPLICA`。

**优化不稳定。** 能量剧烈跳动或 DIIS 步大量被拒绝时，依次尝试：减小 `N_DIIS`（如 3）；减小 `STEPSIZE`/`MAX_STEPSIZE`；设 `NO_LS`；改用 `OPT_TYPE MD` 加退火；检查是否有副本的 SCF 实际上没收敛（看对应的 `-BAND*.out`）。

---

## 8. Dimer 方法精修过渡态

### 8.1 什么时候用 Dimer

- 已经有了靠近鞍点的结构（CI-NEB 的爬坡副本、化学直觉搭的结构、约束扫描的最高点），想要更精确的鞍点；
- 不知道或不关心末态（例如探索一个吸附物可能的分解通道）；
- 体系很大，NEB 的多副本太贵：Dimer 每步只需 2–3 次力计算，并行时可以把全部进程给一个结构。

Dimer 不能保证收敛到“你想要的”鞍点，它会找离起点最近、沿初始方向最容易到达的那个一阶鞍点。所以起点和初始方向越好，结果越可靠。

### 8.2 第一步（2026.2）：用部分 Hessian 得到虚频模式

`ts_guess_vib.inp` 对 TS 猜测结构做一次很便宜的振动分析：固定全部 27 个 Cu，只位移 H 的 3 个坐标（中心差分共 6 次力计算）。它的唯一目的，是把 H 的“反应坐标”模式写进 molden 文件：

```fortran
@SET COORD ${COORD-ts_guess.xyz}

&GLOBAL
  PROJECT ts_guess_vib
  RUN_TYPE VIBRATIONAL_ANALYSIS
  PRINT_LEVEL LOW
&END GLOBAL

@INCLUDE 'cu111_h_force_eval.inc'

&VIBRATIONAL_ANALYSIS
  DX 0.01
  NPROC_REP ${NPROC_REP-8}
  &PRINT
    &MOLDEN_VIB ON
    &END MOLDEN_VIB
  &END PRINT
&END VIBRATIONAL_ANALYSIS

&MOTION
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..27
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

```bash
mpirun -np 48 cp2k.psmp -E COORD=ts_from_neb.xyz -i ts_guess_vib.inp -o ts_guess_vib.out
grep "VIB|Frequency" ts_guess_vib.out        # 应有 1 个负值（虚频）
ls ts_guess_vib-VIBRATIONS-1.mol
```

示例目录中的 `ts_guess.xyz` 是按几何关系放在桥位的猜测结构，可以在没有 NEB 结果时直接使用；有 CI-NEB 结果时，请用 `-E COORD=ts_from_neb.xyz` 换成爬坡副本。

CP2K 写 molden 时，固定原子的位移分量全部写成 0，所以 molden 中仍然列出全部 28 个原子，可以直接用于全体系的 Dimer。各模式按频率从低到高编号，虚频（负值）排在最前面，因此“模式 1”就是虚频模式。振动分析的细节见 [第14章](14-vibration-thermo.md)。

### 8.3 第二步：Dimer 输入

```fortran
@SET COORD ${COORD-ts_guess.xyz}

&GLOBAL
  PROJECT h_cu111_dimer
  RUN_TYPE GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

@INCLUDE 'cu111_h_force_eval.inc'

&MOTION
  &GEO_OPT
    TYPE TRANSITION_STATE
    OPTIMIZER CG
    MAX_ITER 300
    MAX_FORCE 2.0E-4
    RMS_FORCE 1.0E-4
    &CG
      &LINE_SEARCH
        TYPE 2PNT
      &END LINE_SEARCH
    &END CG
    &TRANSITION_STATE
      METHOD DIMER
      &DIMER
        DR [angstrom] 0.01
        ANGLE_TOLERANCE [deg] 2.0
        INTERPOLATE_GRADIENT T
        INITIALIZATION_METHOD MOLDEN
        VIB_MOLDEN_NAME ts_guess_vib-VIBRATIONS-1.mol
        VIB_INDEX 1
        VIB_WEIGHT 1.0
        &ROT_OPT
          OPTIMIZER CG
          &CG
            &LINE_SEARCH
              TYPE 2PNT
            &END LINE_SEARCH
          &END CG
        &END ROT_OPT
      &END DIMER
    &END TRANSITION_STATE
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..9
    &END FIXED_ATOMS
  &END CONSTRAINT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
    &END RESTART
  &END PRINT
&END MOTION
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `&GEO_OPT/TYPE` | `TRANSITION_STATE` | 默认 `MINIMIZATION`；改成 TS 搜索后需要 `&TRANSITION_STATE` 段 |
| `&GEO_OPT/OPTIMIZER` | `CG` | Dimer 的平移步用哪个优化器。官方回归测试（`tests/QS/regtest-gpw-2-3/hcn_ts*.inp` 等）一律用 CG + `2PNT` 线搜索，建议照此设置 |
| `MAX_FORCE`/`RMS_FORCE` | 2.0E-4/1.0E-4 | 收敛判据与普通 GEO_OPT 相同，单位 Ha/bohr。TS 后面要做频率，宜收紧 |
| `&TRANSITION_STATE/METHOD` | `DIMER` | 目前唯一的选项 |
| `DR` | 0.01 Å | dimer 半长度 $\Delta R$（默认 0.01 Å）。太大曲率估计不准，太小受 SCF 噪声影响 |
| `ANGLE_TOLERANCE` | 2° | 转动收敛判据：转动角小于该值即停止转动（默认约 5°，即 0.0873 rad） |
| `INTERPOLATE_GRADIENT` | `T`（默认） | 尽量插值梯度，每次转动省一次力计算 |
| `INITIALIZATION_METHOD` | `MOLDEN`（**2026.2+**） | 默认 `RANDOM`。`MOLDEN` 时用 `VIB_INDEX` 指定的振动模式按 `VIB_WEIGHT` 线性组合成初始方向 |
| `VIB_MOLDEN_NAME` | molden 文件 | 由 `VIBRATIONAL_ANALYSIS/PRINT/MOLDEN_VIB` 生成 |
| `VIB_INDEX`、`VIB_WEIGHT` | `1`、`1.0` | 可写多个，如 `VIB_INDEX 1 2` + `VIB_WEIGHT 0.4 0.6` |
| `&DIMER_VECTOR` | 不写 | 显式给出初始方向（每行一个原子的 3 个分量）；存在时覆盖 `INITIALIZATION_METHOD`，重启文件中会自动写入 |
| `K-DIMER`、`BETA` | 不写 | 打开约束平移的 K-Dimer 及其开关函数指数 |
| `&ROT_OPT` | CG + 2PNT | 转动子优化的设置，结构与 `&GEO_OPT` 相同 |

> [!NOTE]
> **CP2K ≤ 2026.1** 没有 `INITIALIZATION_METHOD`、`VIB_*` 这几个关键词（会报未知关键字），请删掉这四行。这时初始方向是随机的；更好的办法是自己生成 `&DIMER_VECTOR`，例如用 CI-NEB 中爬坡副本两侧的副本之差（即路径切向）：
>
> ```python
> import numpy as np
> def xyz(fn):
>     l = open(fn).read().split("\n"); n = int(l[0])
>     return np.array([[float(v) for v in s.split()[1:4]] for s in l[2:2 + n]])
> tau = xyz("replica_5.xyz") - xyz("replica_3.xyz")   # neighbours of the climbing image
> tau /= np.linalg.norm(tau)
> print("        &DIMER_VECTOR")
> for v in tau:
>     print(f"          {v[0]:12.8f} {v[1]:12.8f} {v[2]:12.8f}")
> print("        &END DIMER_VECTOR")
> ```
>
> 把输出粘贴进 `&DIMER` 段即可（行数必须等于原子数，否则报 `Size of input DIMER_VECTOR inconsistent with natom * 3`）。固定原子对应的分量会被程序自动置零，向量也会被自动归一化。

### 8.4 运行与判断

```bash
mpirun -np 32 cp2k.psmp -i dimer_ts.inp -o dimer_ts.out
grep "DIMER|" dimer_ts.out | head                  # 初始化方式、固定原子数
grep "OPT| Local curvature" dimer_ts.out | tail -5 # 应为负值
grep "OPT| Maximum gradient " dimer_ts.out | tail -3
grep "GEOMETRY OPTIMIZATION COMPLETED" dimer_ts.out
bash last_frame.sh h_cu111_dimer-pos-1.xyz ts_dimer.xyz
```

输出由两类块交替组成：`OPTIMIZING DIMER ROTATION` 块（转动步，打印 `OPT| Local curvature`、`OPT| Total rotational force` 和角度收敛信息），以及普通的平移步块（`OPT| Step number`、`OPT| Maximum gradient` 等，与 GEO_OPT 相同）。判断要点：

- **曲率必须为负**并大致稳定。如果曲率始终为正，说明起点不在鞍点附近（处在势阱里），Dimer 只会沿最平缓的方向爬，未必到达你要的鞍点。这时应换更好的起点（CI-NEB 爬坡副本）；
- 能量应收敛到接近 CI-NEB 最高点的值（Dimer 结果通常略高一点，因为 CI-NEB 未完全收敛）；
- 最终结构中 H 应位于桥位附近，且两端 Cu–H 距离近似相等（本体系的对称性）；
- 2026.2 也会写出 `h_cu111_dimer-FINAL-1_<步数>.xyz`。

---

## 9. 验证过渡态：只有一个虚频

无论结构来自 CI-NEB 还是 Dimer，最后都要做振动分析。本仓库第14章的示例 [ts_freq_check.inp](../examples/14-vibration-thermo/ts_freq_check.inp) 就是针对本例 TS 的频率计算，原理和设置见 [第14章 §7](14-vibration-thermo.md)。要点：

1. 输入的力场设置与 NEB/Dimer 相同，`EPS_SCF` 收紧到 1E-8 左右；
2. 至少放开 H 和与它成键的表面原子；只放开 H 的部分 Hessian 便宜，但忽略了表面原子的耦合，作为快速检查足够，正式发表请做更大的可动区域并检验收敛；
3. 输出中 `VIB|Frequency (cm^-1)` 行应只有**一个负值**，其模式（用 molden 可视化）应把 H 从 fcc 推向 hcp；
4. 若出现额外的小虚频（几十 cm⁻¹）：先怀疑数值噪声（SCF、`EPS_DEFAULT`、`DX`、CUTOFF），再怀疑结构没收敛干净。若出现两个明显的虚频，说明找到的是二阶鞍点，应沿第二个虚频模式扰动后重新做 Dimer。

**模拟 IRC（“伪 IRC”）。** CP2K 的 `RUN_TYPE` 中没有 IRC（内禀反应坐标）任务。实践中可以沿虚频模式把 TS 结构分别向正、负方向推 0.05–0.1 Å，然后各做一次普通 GEO_OPT，检查两边是否分别落到初态和末态。如果两边都回到同一个极小点，说明这个鞍点连接的不是你要的两个态。

**零点能校正后的能垒。** 有了 IS 和 TS 的频率（TS 去掉虚频），

$$
E_a^\mathrm{ZPE} = \left[E^\ddagger + \sum_{j\neq\mathrm{imag}}\tfrac12 h\nu_j^\ddagger\right] - \left[E^\mathrm{IS} + \sum_{i}\tfrac12 h\nu_i^\mathrm{IS}\right]
$$

对含 H 的反应，ZPE 校正通常在几十 meV 到 0.1 eV 量级，不可忽略。Vineyard 指前因子用 §2.1 的公式，两边要用**同一组可动原子**。计算方法和脚本见 [第14章](14-vibration-thermo.md)。

---

## 参数选择与收敛建议

| 项目 | 建议 | 理由 |
|---|---|---|
| 端点优化 | `MAX_FORCE` ≤ 2E-4 Ha/bohr，与 NEB 同一 `FORCE_EVAL`、同样固定原子 | 端点残余力会直接进入能垒误差 |
| `EPS_SCF` | 1E-7（NEB）、1E-7–1E-8（Dimer/频率） | 力噪声决定 NEB/Dimer 能收敛到多紧 |
| `NUMBER_OF_REPLICA` | 简单扩散 5–9；复杂反应 9–16；路径长度/副本数 ≲ 0.5–1 Å（3N 维距离） | 副本太少，鞍点附近分辨不够；太多则成本线性增加 |
| `K_SPRING` | 0.02–0.1（默认 0.02） | 副本分布不均或向端点滑动时增大 |
| `BAND_TYPE` | 先 `IT-NEB` 粗收敛，再 `CI-NEB` | 避免过早爬坡选错最高点 |
| `NSTEPS_IT` | 5–30 | 初始路径越差，取得越大 |
| 收敛判据 | `MAX_FORCE` 5E-4–1E-3 Ha/bohr，位移判据适当放宽 | 0.025–0.05 eV/Å 是常用标准 |
| 周期体系 | `ALIGN_FRAMES F`、`ROTATE_FRAMES F` | 周期晶胞中整体转动没有意义 |
| 并行 | 进程数 = 副本数 × `NPROC_REP` | 否则有进程空转或多轮串行 |
| Dimer | CG + `2PNT`；`DR` 0.005–0.01 Å；`ANGLE_TOLERANCE` 1–5° | 与官方回归测试一致；转动收敛太松会导致平移方向偏差 |
| Dimer 初始方向 | 2026.2：`INITIALIZATION_METHOD MOLDEN`；旧版：用 NEB 切向写 `&DIMER_VECTOR` | 随机方向在大体系中会严重拖慢收敛甚至扭曲结构 |
| 验证 | 频率计算恰有一个虚频，伪 IRC 两侧回到 IS/FS | 鞍点的定义 |
| 金属表面 | 对角化 + 展宽，k 点收敛，固定底层 | OT 不支持展宽（≤ 2026.2） |

---

## 常见问题与报错

**1. `Number of totally requested processors (30) is not compatible with the number of processors requested per replica (8) and the number of replicas (3) . [6] processors will be wasted!`**
原因：总进程数不是 `NPROC_REP` 的整数倍。解决：让 MPI 进程数 = `NUMBER_OF_REPLICA × NPROC_REP`（或至少是 `NPROC_REP` 的整数倍，且组数整除副本数）。

**2. `Number of lines in XYZ format not equal to the number of atoms. Error in XYZ format for REPLICA coordinates. Very probably the line with title is missing or is empty.`**
原因：`&REPLICA/COORD_FILE_NAME` 指向的 xyz 文件缺少标题行、标题行为空，或原子数写错。解决：确保第二行有任意非空文字（`last_frame.sh` 生成的文件保留了 CP2K 的标题行）。

**3. 主输出只有寥寥几行，或者看不到汇总表/收敛信息。**
原因：2026.1 及以前，`&PROGRAM_RUN_INFO` 和 `&CONVERGENCE_INFO` 默认在 `MEDIUM` 级别才打印，而 BAND 常用 `PRINT_LEVEL LOW`；收敛信息在 `-BAND*.out` 中。解决：像本例一样显式写 `&PROGRAM_RUN_INFO ON`、`&CONVERGENCE_INFO ON`，或在 `-BAND*.out` 中查找。

**4. 某个副本报 `SCF run NOT converged`，整个 BAND 中止。**
原因：线性插值产生了原子过近的构型，或者金属副本的 SCF 设置不够稳。报错只出现在对应的 `-BAND*.out` 中。解决：检查该副本结构；用 IDPP 插值或手动给出中间副本；按 [第5章](05-scf-strategy.md) 调整混合与展宽。不要用 `IGNORE_CONVERGENCE_FAILURE` 掩盖问题。

**5. 路径上出现两个能量峰，或者最高副本紧挨着端点。**
原因：路径经过中间体；或端点本身不是极小点（没优化好、优化到了别的构型）。解决：把谷底副本拿出来单独优化，拆成两段 NEB；重新检查端点优化。

**6. slab 的副本被整体转动、固定原子“动了”，或者能量出现不连续跳变。**
原因：周期体系保留了默认的 `ALIGN_FRAMES T`/`ROTATE_FRAMES T`。解决：两者都设为 `F`。

**7. 插值路径中某个原子穿过了整个晶胞。**
原因：两个端点中该原子处于不同的周期映像。解决：按 §3 的脚本逐原子检查位移，把出问题的坐标平移一个晶格矢量，使两端坐标连续。

**8. Dimer 报 `INITIALIZATION_METHOD MOLDEN requires VIB_MOLDEN_NAME` 或 `Could not found <vibration 1> from molden file`。**
原因：没给 molden 文件名，或文件不是 `VIBRATIONAL_ANALYSIS` 写出的格式。解决：用 `&MOLDEN_VIB` 生成的 `*-VIBRATIONS-1.mol`，文件名与路径要正确。

**9. Dimer 报 `The norm of the dimer vector is 0! Calculation cannot proceed further.`**
原因：初始 dimer 向量的所有非零分量都落在被固定的原子上（例如振动分析放开的原子在 Dimer 中被固定了），置零后向量为零。解决：保证振动分析的可动原子在 Dimer 中也可动。

**10. Dimer 曲率一直为正，结构滑回极小点或跑向别的构型。**
原因：起点离鞍点太远，或初始方向不对。解决：用 CI-NEB 的爬坡副本作起点，用虚频模式初始化方向；减小步长（CG 线搜索已较保守），必要时先用约束扫描找到更好的起点。

**11. 频率验证出现 2 个以上虚频。**
原因：二阶鞍点，或数值噪声。解决：小于 50 cm⁻¹ 的虚频先通过收紧 `EPS_SCF`、`EPS_DEFAULT`、增大 CUTOFF、调整 `DX` 排查；大的第二虚频则沿该模式扰动后重新做 Dimer。

---

## 练习

1. 把 `NUMBER_OF_REPLICA` 改成 5 和 9，比较能垒和收敛步数；再把 `BAND_TYPE` 改为 `IT-NEB`，比较最高副本能量与 CI-NEB 的差别，体会“普通 NEB 低估能垒”。
2. 分别用 `K_SPRING` 0.01、0.05、0.2 运行，观察 `DISTANCES REP` 的均匀程度和收敛速度。
3. 用 28 个 MPI 进程分别设置 `NPROC_REP` 为 4、7、8，记录每个 BAND 步的墙钟时间，验证 §5 的表。
4. 不做 NEB，直接从 `ts_guess.xyz` 出发，分别用随机初始方向（删掉 MOLDEN 相关四行）和 MOLDEN 初始化跑 Dimer，比较收敛所需的力计算次数（需要 CP2K 2026.2）。
5. 对收敛的 TS 做“伪 IRC”：沿虚频模式正负各推 0.1 Å 再优化，确认分别回到 fcc 和 hcp 空位。
6. 加上 `&KPOINTS`（MONKHORST-PACK 3 3 1），重新计算端点和 TS 的单点能，评估 Γ 点近似对能垒的影响。

---

## 参考资料

**官方文档**

- 输入参考 `&MOTION/&BAND`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/BAND.html>（子段 `OPTIMIZE_BAND`、`CI_NEB`、`CONVERGENCE_CONTROL`、`REPLICA`、`PROGRAM_RUN_INFO`、`FINAL_BAND` 在同目录下）
- 输入参考 `&GEO_OPT/&TRANSITION_STATE/&DIMER`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/GEO_OPT/TRANSITION_STATE/DIMER.html>
- Geometry and Cell Optimization（含 `TYPE TRANSITION_STATE` 与输出文件说明）：<https://manual.cp2k.org/trunk/methods/optimization/geometry_and_cell_opt.html>
- Nudged Elastic Band（官方页面目前只给出练习链接）：<https://manual.cp2k.org/trunk/methods/optimization/nudged_elastic_band.html>
- 旧版练习：NEB（SN2 反应，PM6）：<https://www.cp2k.org/exercises:common:neb>
- CP2K 2026.2 release notes（NEB 新输出、Dimer 的 molden 初始化）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- PR #5382 “NEB: revamp output”：<https://github.com/cp2k/cp2k/pull/5382>；PR #5312 “Dimer: add control of initialization method and implement read from molden”：<https://github.com/cp2k/cp2k/pull/5312>
- Restarting：<https://manual.cp2k.org/trunk/methods/restarting.html>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/NEB/`（IT/CI/B/D-NEB、DIIS/MD 优化器、重启）、`tests/QS/regtest-gpw-2-3/hcn_ts*.inp`（Quickstep Dimer、K-Dimer、固定原子）、`tests/Fist/regtest-5-vib/argon_dimer_molden_read.inp`（molden 初始化 Dimer）是查关键字真实写法的最佳来源
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，有 CP2K NEB 的实践经验
- [csc-training/sscc-cp2k](https://github.com/csc-training/sscc-cp2k)：一条完整的反应研究链：PBE 下做 NEB → TS 振动分析 → PBE0 重新收敛 MEP → 溶液中 AIMD
- [arvpinto/enzyme_neb_cp2k](https://github.com/arvpinto/enzyme_neb_cp2k)：酶反应 QM/MM 的 CI-NEB、Dimer 与频率验证教程（CC0），附能量与频率分析脚本
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：NEB 初始化与分析的中文脚本思路（无 license，仅参考）
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：`cp2kparse` 解析输出为 JSON/YAML

**本仓库脚本**：[neb_plot.py](../scripts/neb_plot.py)（MEP 作图、导出 TS 结构）、[vib_thermo.py](../scripts/vib_thermo.py)（频率与 ZPE）、[cp2k_outparse.py](../scripts/cp2k_outparse.py)（解析 OPT| 块）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- G. Henkelman, H. Jónsson, *J. Chem. Phys.* **113**, 9978 (2000)，改进切向 NEB（IT-NEB）：<https://doi.org/10.1063/1.1323224>
- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* **113**, 9901 (2000)，CI-NEB：<https://doi.org/10.1063/1.1329672>
- G. Henkelman, H. Jónsson, *J. Chem. Phys.* **111**, 7010 (1999)，Dimer 方法：<https://doi.org/10.1063/1.480097>
- S. A. Trygubenko, D. J. Wales, *J. Chem. Phys.* **120**, 2082 (2004)，doubly nudged elastic band：<https://doi.org/10.1063/1.1636455>
- W. E, W. Ren, E. Vanden-Eijnden, *Phys. Rev. B* **66**, 052301 (2002)，string method：<https://doi.org/10.1103/PhysRevB.66.052301>
- S. Smidstrup, A. Pedersen, K. Stokbro, H. Jónsson, *J. Chem. Phys.* **140**, 214106 (2014)，IDPP 插值：<https://doi.org/10.1063/1.4878664>
- G. H. Vineyard, *J. Phys. Chem. Solids* **3**, 121 (1957)，谐振子过渡态理论：<https://doi.org/10.1016/0022-3697(57)90059-8>
- K-Dimer：*J. Chem. Phys.* **141**, 164111 (2014)（CP2K 输入参考中 `K-DIMER` 关键字所引文献）
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](12-aimd.md) | [返回目录](../README.md) | [下一章 →](14-vibration-thermo.md)
