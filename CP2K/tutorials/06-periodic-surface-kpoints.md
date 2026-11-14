# 第6章 周期性体系、表面模型与 k 点（Periodic Systems, Slab Models & k-Points in CP2K）

本章讲 CP2K 中周期性体系的三项基本设置：用 `&CELL` 定义晶胞（`ABC`/`ALPHA_BETA_GAMMA`、`A/B/C` 向量、`CELL_FILE_NAME`、`MULTIPLE_UNIT_CELL`），让 `PERIODIC` 与 `&POISSON` 求解器（PERIODIC、MT、WAVELET、ANALYTIC）保持一致，以及用 `&KPOINTS` 做 Monkhorst–Pack 布里渊区采样。我们用三个可运行示例演示孤立分子（H2O + MT）、带偶极校正的表面 slab（Pt(111)+O，`SURFACE_DIPOLE_CORRECTION`）和体相 Si 原胞的 k 点收敛测试，并说明 CP2K 2026.2 新增的 k 点对称约化、k 点波函数外推等功能，以及“什么时候 Γ 点就够、什么时候必须收敛 k 网格”。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 用三种方式写 `&CELL`，知道它们的优先级，避开 `SYMMETRY HEXAGONAL` 和 `MULTIPLE_UNIT_CELL` 的常见陷阱；
- 为体相、分子、slab、纳米线选择匹配的 `PERIODIC` 和 `POISSON_SOLVER`；
- 搭建一个带真空层、偶极校正和固定底层原子的 slab 输入，并检查真空是否足够；
- 用 `&KPOINTS` 设置 Monkhorst–Pack 网格，读懂 `BRILLOUIN|` 输出，完成 k 点收敛测试；
- 了解 2026.2 中 k 点相关功能的状态（哪些支持、哪些是实验性、哪些仍不支持）。

| 体系类型 | 本章示例 | 关键设置 |
|---|---|---|
| 体相晶体（小原胞） | `si_primitive_kpoints.inp`、`kpoint_convergence_template.inp` | `PERIODIC XYZ` + `&KPOINTS` + 对角化 |
| 孤立分子 | `h2o_isolated_mt.inp` | `PERIODIC NONE` + `POISSON_SOLVER MT` + `&CENTER_COORDINATES` |
| 表面/吸附 | `pt111_slab_dipole.inp` | 真空 ≥ 15 Å + `SURFACE_DIPOLE_CORRECTION` + k 点 n×n×1 |

示例输入集：[../examples/06-periodic-surface-kpoints/README.md](../examples/06-periodic-surface-kpoints/README.md)。

---

## 2. 原理速览

### 2.1 周期边界与布里渊区积分

周期体系的单电子态满足 Bloch 定理，用晶体动量 $\mathbf k$ 标记。电子密度、总能量等都要在第一布里渊区上积分，数值上用加权求和近似：

$$
\frac{1}{\Omega_\mathrm{BZ}}\int_\mathrm{BZ} f(\mathbf k)\,\mathrm d\mathbf k \;\approx\; \sum_{\mathbf k} w_{\mathbf k}\, f(\mathbf k),\qquad \sum_{\mathbf k} w_{\mathbf k} = 1
$$

倒格矢长度与实空间晶胞长度成反比：$|\mathbf b_i| = 2\pi/a_i$（正交晶胞）。晶胞越大，布里渊区越小，需要的 k 点越少。超胞足够大时只取 $\mathbf k = 0$（Γ 点）即可，这也是 CP2K 最常用、最高效的模式，所有功能都支持。

**超胞与 k 网格的等价性**：在 $N_1\times N_2\times N_3$ 超胞上只算 Γ 点，采样到的 k 点与在原胞上用以 Γ 为中心的 $N_1\times N_2\times N_3$ 网格完全相同（能带折叠）。所以“扩胞”和“加 k 点”是同一件事的两种做法：前者能用 OT 和全部功能，后者对小胞更省。

### 2.2 静电边界条件

Hartree 势由 Poisson 方程 $\nabla^2 V_\mathrm H = -4\pi\rho$ 决定。三维周期方向用 FFT 求解最自然；如果某个方向在物理上不是周期的（分子的三个方向、slab 的法向），就需要专门的 Green 函数或“去耦合”方法，否则体系会和它的周期像发生静电相互作用：

- **偶极–偶极作用**：极性分子或非对称 slab 的周期像之间存在虚假的电场，误差随盒子尺寸以 $1/L^3$（分子）衰减；对 slab 则表现为真空区存在一个人为的线性电势；
- **净电荷**：3D 周期下带电体系要引入均匀补偿背景，能量随晶胞大小缓慢收敛。

CP2K 提供两类处理：直接换成非周期的 Poisson 求解器（MT、WAVELET、ANALYTIC、MULTIPOLE、IMPLICIT），或者保持 3D 周期，对 slab 加偶极校正（Neugebauer–Scheffler / Bengtsson 方案，即 `SURFACE_DIPOLE_CORRECTION`）。

---

## 3. `&CELL`：晶胞的几种写法

### 3.1 三种写法与优先级

```text
&CELL                                   ! (1) lengths + angles
  ABC              5.614428 5.614428 26.0     ! Angstrom (default unit)
  ALPHA_BETA_GAMMA 90.0 90.0 60.0             ! alpha = angle(B,C), beta = angle(A,C), gamma = angle(A,B)
  PERIODIC XYZ
&END CELL

&CELL                                   ! (2) explicit Cartesian vectors (columns of h)
  A 5.614428 0.000000  0.0
  B 2.807214 4.862237  0.0
  C 0.000000 0.000000 26.0
&END CELL

&CELL                                   ! (3) read from a file
  CELL_FILE_NAME   structure.cif
  CELL_FILE_FORMAT CIF                  ! CP2K | CIF | XSC; 2026.2 adds EXTXYZ | XYZ | PDB
&END CELL
```

| 关键字 | 默认 | 说明 |
|---|---|---|
| `ABC` | — | 三个晶胞向量的长度（Å）。约定 A 沿 x 轴、B 在 xy 平面内 |
| `ALPHA_BETA_GAMMA`（别名 `ANGLES`） | 90 90 90 [deg] | α 是 B 与 C 的夹角，β 是 A 与 C，γ 是 A 与 B |
| `A`、`B`、`C` | — | 显式笛卡尔向量（构成 h 矩阵的列）。写了 `ABC` 或 `CELL_FILE_NAME` 时会被忽略 |
| `CELL_FILE_NAME` / `CELL_FILE_FORMAT` | — / `CP2K` | 从文件读取晶胞；≤2026.1 支持 CP2K/CIF/XSC，2026.2 新增 EXTXYZ（读取 `Lattice=`）、XYZ（EXTXYZ 别名）、PDB（读 CRYST1） |
| `PERIODIC` | `XYZ` | 原子体系的周期方向（近邻列表、坐标折回）。静电部分的周期性另由 `&POISSON/PERIODIC` 指定 |
| `SYMMETRY` | `NONE` | 强制晶胞对称性（通过长度和角度重建 h 矩阵），不改变原子坐标 |
| `MULTIPLE_UNIT_CELL` | 1 1 1 | 把晶胞当作单胞沿三个晶胞向量方向复制；只影响晶胞，坐标的复制要另设 `&TOPOLOGY/MULTIPLE_UNIT_CELL` |
| `CANONICALIZE` | `AUTO` | **2026.2 新增**：是否把一般取向的输入晶胞转换到“A 沿 x、B 在 xy 面”的内部约定 |
| `&CELL_REF` | — | 参考晶胞，用于 CELL_OPT/NPT 时固定网格（见[第4章](04-cutoff-convergence.md)） |

2026.2 的输入参考给出的解析优先级是：`CELL_FILE_NAME` > `ABC`/`ALPHA_BETA_GAMMA` > `A`/`B`/`C`。2026.2 起，如果三者都没写，还会尝试从 `&TOPOLOGY/COORD_FILE_NAME` 指定的文件（如 CIF、EXTXYZ、PDB）中读取晶胞。

> [!WARNING]
> `SYMMETRY HEXAGONAL` 是 `HEXAGONAL_GAMMA_60` 的别名，即 **γ = 60°**；需要 γ = 120° 时必须写 `HEXAGONAL_GAMMA_120`。两种设定描述的是同一个六方格子，但晶胞向量的方向不同，坐标必须与之对应。从 Materials Project 下载的六方结构通常是 γ = 120°。不确定时干脆不写 `SYMMETRY`，直接给出 `ABC` 和 `ALPHA_BETA_GAMMA`。

### 3.2 分数坐标

```text
&COORD
  SCALED                     ! coordinates below are fractional (in units of A, B, C)
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
&END COORD
```

`SCALED` 的分数坐标与晶胞取向无关，最不容易出错。笛卡尔坐标（默认单位 Å）必须与 `ABC`/`ALPHA_BETA_GAMMA` 所隐含的取向一致，即 A 沿 x、B 在 xy 平面。从其它软件导出的笛卡尔坐标若来自不同取向，要么改用分数坐标，要么用 `A`/`B`/`C` 原样写出对方的晶胞向量。

### 3.3 MULTIPLE_UNIT_CELL：两处必须同时设置

用原胞直接生成超胞：

```text
&SUBSYS
  &CELL
    ABC              3.840297 3.840297 3.840297
    ALPHA_BETA_GAMMA 60 60 60
    MULTIPLE_UNIT_CELL 2 2 2           ! enlarges the CELL ...
  &END CELL
  &TOPOLOGY
    MULTIPLE_UNIT_CELL 2 2 2           ! ... and replicates the COORDINATES; both are required
  &END TOPOLOGY
  &COORD                               ! Cartesian coordinates (Angstrom) of the 2-atom unit cell
    Si 0.000000 0.000000 0.000000
    Si 1.920149 1.108598 0.783897
  &END COORD
  ...
&END SUBSYS
```

两个关键字分工不同：`&CELL/MULTIPLE_UNIT_CELL` 只放大晶胞，`&TOPOLOGY/MULTIPLE_UNIT_CELL` 只复制坐标（按单胞的三个晶胞向量平移），两者的数值必须相同：

- 只写 `&TOPOLOGY` 中的那个（或两处数值不同），程序会中止并报 `The input parameters for SUBSYS%TOPOLOGY%MULTIPLE_UNIT_CELL and SUBSYS%CELL%MULTIPLE_UNIT_CELL have to agree.`；
- 只写 `&CELL` 中的那个，程序**不会报错**，但得到的是一个放大的盒子，里面只有原胞那几个原子。这种错误最隐蔽，务必查看输出中的原子数。

上例的坐标是按 `ABC` + `ALPHA_BETA_GAMMA` 的取向（A 沿 x、B 在 xy 平面）换算出的笛卡尔坐标，复制时直接加上单胞向量的整数倍，含义最明确。

### 3.4 检查晶胞

```bash
cp2k.psmp -e input.inp | grep -A8 "&CELL"    # 展开后的晶胞设置
grep -i "CELL|" output.out | head -20         # 程序实际使用的晶胞向量、长度、角度和体积
```

建议用 VESTA、OVITO 或 ASE GUI 连同周期像一起看一遍结构：真空方向不应有原子“贴”到盒子边界，周期方向上相邻像之间既不能重叠，也不能出现人为的空隙。建模和格式转换见 [Common 第3章 结构准备](../../Common/tutorials/03-structure-modeling.md) 和 [Common 第4章 ASE 与 pymatgen](../../Common/tutorials/04-ase-pymatgen.md)。

---

## 4. 周期性与 Poisson 求解器的搭配

CP2K 中有两个独立的 `PERIODIC`：`&SUBSYS/&CELL/PERIODIC` 管原子体系（近邻列表、坐标折回），`&DFT/&POISSON/PERIODIC` 管静电 Green 函数。**两者通常必须一致**。官方文档明确提醒：不要靠让两者不一致来“压掉”真空方向的相互作用。

| 体系 | `&CELL PERIODIC` | `&POISSON PERIODIC` | `POISSON_SOLVER` | 要点 |
|---|---|---|---|---|
| 体相晶体、液体 | `XYZ` | `XYZ` | `PERIODIC`（默认） | 只能用于 3D 周期 |
| 分子、团簇 | `NONE` | `NONE` | `MT`（常用）、`WAVELET`、`MULTIPOLE`、`ANALYTIC` | MT 要求盒子 ≥ 电荷分布尺寸的 2 倍 |
| Slab（方式 A，最常用） | `XYZ` | `XYZ` | `PERIODIC` + `SURFACE_DIPOLE_CORRECTION` | 需要足够的真空，法向必须沿某个笛卡尔轴 |
| Slab（方式 B，真 2D） | `XY` | `XY` | `MT`、`ANALYTIC`、`WAVELET` | WAVELET 的 2D 组合：≤2026.1 只允许 `XZ`，2026.2 允许 XY/XZ/YZ |
| 纳米线 | 如 `Z` | 如 `Z` | `ANALYTIC`（1D）或 `IMPLICIT` | IMPLICIT 的边界条件要在其子段中另设 |

各求解器的适用范围（据 2026.2 输入参考整理）：

| `POISSON_SOLVER` | 维度 | 说明 |
|---|---|---|
| `PERIODIC` | 3D | 标准倒空间解法 |
| `ANALYTIC` | 0D、1D、2D | 倒空间解析 Green 函数，收敛较慢 |
| `MT`（Martyna–Tuckerman） | 0D、2D | 盒子必须至少是（含尾部的）电子密度尺寸的 2 倍，否则结果**完全错误且不报错**；`&MT/ALPHA`（默认 7.0）、`&MT/REL_CUTOFF`（默认 2.0）一般不用改 |
| `MULTIPOLE` | 0D | 每个原子用一个高斯拟合总电荷 |
| `WAVELET` | 0D、2D、3D | 只要求非周期方向的盒面上密度为零，不需要特别大的盒子；输入参考注明需要 `&GLOBAL PREFERRED_FFT_LIBRARY FFTSG` |
| `IMPLICIT` | 0D–3D | 广义边界条件与介电模型；其 `BOUNDARY_CONDITIONS` 独立设置，默认是周期的 |

> [!IMPORTANT]
> 带电体系的能量与边界条件绑定：3D `PERIODIC` 会自动加均匀补偿背景，能量随晶胞大小以代数方式缓慢变化；0D 求解器（MT/WAVELET）可以给出孤立带电团簇的有限能量；带净电荷的 slab 或纳米线在没有对电极或补偿模型时，能量甚至不存在与真空厚度无关的极限。不同 Poisson 求解器算出的带电体系总能量**不在同一参考上**，不能直接相减。

---

## 5. 孤立分子：H2O + MT 求解器（h2o_isolated_mt.inp）

### 5.1 完整输入

```text
# Isolated H2O molecule: non-periodic cell + Martyna-Tuckerman Poisson solver
# Change BOX to test the box-size convergence (10, 12, 14, 16 Angstrom).
@SET BOX 12.0

&GLOBAL
  PROJECT     h2o_mt
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC       NONE          ! must match &CELL PERIODIC
      POISSON_SOLVER MT            ! box must be >= 2x the extent of the electron density
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   30
      &OT ON
        MINIMIZER      DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC ${BOX} ${BOX} ${BOX}
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   h2o.xyz
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES          ! put the molecule at the box centre (default: cell/2)
      &END CENTER_COORDINATES
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

关键点：

| 设置 | 作用 |
|---|---|
| `&CELL PERIODIC NONE` 与 `&POISSON PERIODIC NONE` | 两处一致，原子和静电都按孤立体系处理 |
| `POISSON_SOLVER MT` | Martyna–Tuckerman 去耦合，计算开销与周期计算接近 |
| `ABC ${BOX}` | 盒子边长。H2O 的电子密度（含尾部）约在 5 Å 量级的范围内，12 Å 是比较稳妥的起点 |
| `&CENTER_COORDINATES` | 把分子平移到盒子中心，避免密度“漏”到盒子边界；`CENTER_POINT` 可指定其它中心 |

### 5.2 盒子尺寸收敛

MT 盒子太小时结果会出错，而且**不会报错**。所以对新分子要做一次盒子尺寸测试：

```bash
for L in 10 12 14 16; do
  sed -E "s/^@SET[[:space:]]+BOX[[:space:]].*/@SET BOX ${L}.0/" h2o_isolated_mt.inp > box_${L}.inp
  cp2k.psmp -i box_${L}.inp -o box_${L}.out
  printf "%3s A  %s\n" "$L" "$(grep 'ENERGY| Total FORCE_EVAL' box_${L}.out | awk '{print $NF}')"
done
```

能量应当进入平台（例如相邻两点相差 < 10⁻⁵ Ha）。注意盒子变化时网格点数也会变，测试应在已收敛的 CUTOFF 下进行，以免把网格误差和盒子误差混在一起。

### 5.3 换成 WAVELET 求解器

WAVELET 对盒子的要求较低，只需要非周期方向的盒面上密度为零：

```text
&GLOBAL
  ...
  PREFERRED_FFT_LIBRARY FFTSG      ! required by the WAVELET solver according to the input reference
&END GLOBAL
...
&POISSON
  PERIODIC       NONE
  POISSON_SOLVER WAVELET
  &WAVELET
    SCF_TYPE 40                    ! default; order of the interpolating scaling functions
  &END WAVELET
&END POISSON
```

> [!NOTE]
> 把分子放进一个 `PERIODIC XYZ` 的大盒子里，用默认的周期求解器计算，是很多旧教程的做法。对非极性分子误差往往可以接受，但对水这类极性分子，周期像之间的偶极相互作用要到盒子很大时才消失。做孤立分子的结合能、电离能、电子亲和能等，请用本节的非周期设置。

---

## 6. 表面 slab：Pt(111) + O 吸附与偶极校正（pt111_slab_dipole.inp）

### 6.1 建模要点

- **层数**：金属表面一般 4–6 层，底部 1–2 层固定在体相位置，模拟半无限体相；
- **真空**：≥ 15 Å，以真空中电子密度降到可忽略为准；功函数等性质需要更厚的真空并做收敛；
- **表面超胞**：吸附物覆盖度决定横向大小（本例 p(2×2)，0.25 ML）；
- **k 点**：表面内两个方向按晶格长度取 n×n，法向（真空方向）取 1；
- **非对称 slab**：只有一侧有吸附物，或上下表面不等价时，存在净偶极，需要偶极校正。

本例的结构用 ASE 生成（晶格常数取 3.97 Å，接近常见的 PBE 值；实验值约 3.92 Å）：

```python
from ase.build import fcc111, add_adsorbate

slab = fcc111("Pt", size=(2, 2, 4), a=3.97, vacuum=None, orthogonal=False)
add_adsorbate(slab, "O", height=1.2, position="fcc")      # O in the fcc hollow
cell = slab.get_cell()
cell[2] = [0.0, 0.0, 26.0]                                # c along z, ~18 A of vacuum
slab.set_cell(cell)
slab.translate([0.0, 0.0, 1.0])                           # keep atoms away from z = 0
print(slab.cell.cellpar())                                # -> 5.614 5.614 26.0 90 90 60
slab.write("pt111_o_p2x2.xyz", format="xyz")
```

得到的晶胞是 γ = 60° 的六方表面胞，c 轴沿 z，正好满足偶极校正“表面法向平行于某个笛卡尔轴”的要求。示例目录中的 `pt111_o_p2x2.xyz` 就是这样生成的（初始结构，未优化）。

### 6.2 完整输入

```text
# Pt(111) p(2x2) 4-layer slab with O in the fcc hollow site
# 3D-periodic cell + surface dipole correction along z; k-points 4x4x1; metal -> smearing
# Set RUNTYPE to GEO_OPT to relax the top two layers and O (bottom two layers are fixed).
@SET RUNTYPE ENERGY_FORCE

&GLOBAL
  PROJECT     pt111_o
  RUN_TYPE    ${RUNTYPE}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    SURFACE_DIPOLE_CORRECTION T    ! asymmetric slab -> compensate the dipole
    SURF_DIP_DIR              Z    ! surface normal (must be a Cartesian axis)
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 1  ! one k-point along the vacuum direction
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   200
      ADDED_MOS 60                 ! 294 electrons -> 147 occupied MOs (+~40 %)
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD  BROYDEN_MIXING
        ALPHA   0.2
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC              5.614428 5.614428 26.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   pt111_o_p2x2.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Pt
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
      POTENTIAL GTH-PBE-q18
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL

&MOTION                            ! only used when RUN_TYPE = GEO_OPT
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER  200
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1 2 3 4 5 6 7 8         ! bottom two Pt layers
      COMPONENTS_TO_FIX XYZ
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

### 6.3 关键字说明

| 关键字 | 位置 | 默认 | 本例 | 说明 |
|---|---|---|---|---|
| `SURFACE_DIPOLE_CORRECTION`（别名 `SURF_DIP`） | `&DFT` | F | T | 对非对称 slab 在真空区引入补偿偶极层，消除周期像之间的虚假电场 |
| `SURF_DIP_DIR` | `&DFT` | Z | Z | 表面法向，只能是 X、Y、Z 之一；C 轴必须严格沿该方向 |
| `SURF_DIP_POS` | `&DFT` | -1 | 默认 | 偶极层位置（Å）；-1 表示自动取平面平均密度最小处 |
| `SURF_DIP_SWITCH` | `&DFT` | F | 不用 | 实验性：SCF 难收敛时自动关闭偶极校正 |
| `CORE_CORR_DIP` | `&DFT` | F | 不用 | 与 `&CORE_CORRECTION` 等带电处理配合时使用 |
| `SCHEME MONKHORST-PACK 4 4 1` | `&KPOINTS` | — | — | 表面内 4×4，真空方向 1 |
| `ADDED_MOS` | `&SCF` | 0 | 60 | 金属 slab 需要足够的空轨道；基函数总数为 16 × 26 + 13 = 429 |
| `&FIXED_ATOMS LIST` | `&MOTION/&CONSTRAINT` | — | 1–8 | 底部两层；只在 `GEO_OPT`/`MD` 时生效 |

### 6.4 运行与检查

```bash
mpirun -np 32 cp2k.psmp -i pt111_slab_dipole.inp -o pt111_slab_dipole.out

grep "BRILLOUIN| List of Kpoints" pt111_slab_dipole.out
grep -A1 "Total dipole moment perpendicular" pt111_slab_dipole.out | tail -2
grep -i "smearing includes the last MO"     pt111_slab_dipole.out   # must be empty
grep "SCF run converged"                    pt111_slab_dipole.out
grep "ENERGY| Total FORCE_EVAL"             pt111_slab_dipole.out
```

打开偶极校正后，每个 SCF 步都会打印 `Total dipole moment perpendicular to the slab [electrons-Angstroem]:`，可以借此确认校正确实生效，并看到 slab 偶极的大小。若真空太薄，程序会中止并报：

```text
 Dipole correction needs more vacuum space above the surface
```

（2026.2 源码中的判据是：在偶极层所在的平面上，总电荷密度的平面平均绝对值不小于 10⁻⁵（原子单位）。）若 C 轴不沿 `SURF_DIP_DIR`，或 A、B 在 z 方向有分量，则报 `Dipole correction only for surface perpendicular to one Cartesian axis`。

**真空与偶极校正的验证**：分别用 c = 22、26、30 Å 计算，比较能量和受力；再输出 V_HARTREE cube，沿 z 做平面平均。真空区的静电势应当是两段平台，平台之差对应 slab 两侧功函数之差。做法见 [第11章](11-cube-charge-workfunction.md) 和本仓库的 [cube_tools.py](../scripts/cube_tools.py)（`python ../scripts/cube_tools.py -h`）。

> [!TIP]
> 2D 周期加真正 2D 的 Poisson 求解器（方式 B：`PERIODIC XY` + `POISSON_SOLVER MT` 或 `ANALYTIC`）是偶极校正之外的另一种选择，物理上没有周期像问题。2026.2 修复了 2D 体系在 ANALYTIC/MT 下的解析应力。不过方式 B 用 MT 时同样要求法向盒子 ≥ 2 倍 slab 电荷分布，计算量通常更大。日常的吸附、反应计算用方式 A 即可。

---

## 7. `&KPOINTS`：k 点采样

### 7.1 关键字一览（以 2026.2 为准）

| 关键字 | 默认 | 说明 | 版本 |
|---|---|---|---|
| `SCHEME` | — | `NONE`、`GAMMA`、`MONKHORST-PACK n1 n2 n3`、`MACDONALD n1 n2 n3 s1 s2 s3`、`GENERAL` | 全部 |
| `KPOINT` | — | `GENERAL` 方案下逐行给出 `x y z w`（可重复），权重会自动归一化 | 全部 |
| `UNITS` | `B_VECTOR` | `KPOINT` 坐标单位：倒格矢分数，或 `CART_ANGSTROM`/`CART_BOHR`（2π/长度） | 全部 |
| `SYMMETRY` | F | 用原子（空间群）对称性约化 k 点；2026.2 文档标注为**实验性**功能 | 关键字各版本都有；“k 点对称约化”列在 2026.2 的新功能中 |
| `FULL_GRID` | F | 用完整网格，不做约化 | 全部 |
| `WAVEFUNCTIONS` | `COMPLEX` | `REAL` 只适用于 Γ 和波函数可取实数的特殊点 | 全部 |
| `PARALLEL_GROUP_SIZE` | -1 | 每个 k 点组的 MPI 进程数：-1 自动取满足约束的最小值，0 表示全部进程，n 表示恰好 n 个 | 全部 |
| `EPS_GEO` → `EPS_SYMMETRY` | 1.0E-6 | 对称性判定精度；2026.2 改名为 `EPS_SYMMETRY`，`EPS_GEO` 仍可作别名 | 2026.2 改名 |
| `GAMMA_CENTERED` | F | 生成以 Γ 为中心的 MP/MacDonald 网格（只对偶数分割有意义） | **2026.2 新增** |
| `INVERSION_SYMMETRY_ONLY` | F | 只用 k → −k（时间反演）约化 | **2026.2 新增** |
| `SYMMETRY_BACKEND` / `SYMMETRY_REDUCTION_METHOD` | `K290` | 可选 `SPGLIB`（需编译时启用 spglib） | **2026.2 新增** |
| `VERBOSE` | F | 输出更详细的 k 点信息 | 全部 |

几点容易混淆的地方：

- **完全不写 `&KPOINTS`** 就是传统的 Γ 点计算（`SCHEME NONE`），可以使用 OT 和全部功能，重启文件为 `.wfn`；
- **`SCHEME GAMMA`** 虽然也只有一个 Γ 点，但走的是 k 点代码路径，重启文件为 `.kp`。大超胞 Γ 点计算请直接省略 `&KPOINTS`；
- **偶数 MP 网格不包含 Γ 点**。CP2K 的 Monkhorst–Pack 网格对偶数分割是相对 Γ 平移的；需要包含 Γ 时，2026.2 可用 `GAMMA_CENTERED T`，旧版本改用奇数网格（奇数分割的 MP 网格自然包含 Γ）；
- **默认的约化方式**：`SYMMETRY F` 加 `FULL_GRID F` 时，规则网格上会把 k 与 −k 合并（时间反演），偶数 n×n×n 网格约化为 n³/2 个点。要得到严格的完整网格参考，需同时设 `SYMMETRY F` 和 `FULL_GRID T`。

### 7.2 体相 Si 原胞（si_primitive_kpoints.inp）

```text
# Bulk Si, 2-atom primitive fcc cell, Monkhorst-Pack k-point sampling
# k-point runs require diagonalization in CP2K <= 2026.2 (no OT).
&GLOBAL
  PROJECT     si2_kp
  RUN_TYPE    ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     350
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 8 8 8
      SYMMETRY      F              ! default; atomic symmetry reduction is experimental
      FULL_GRID     F              ! default; k and -k are combined
      WAVEFUNCTIONS COMPLEX        ! default for k-point runs
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   80
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD  BROYDEN_MIXING
        ALPHA   0.4
        NBUFFER 8
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC              3.840297 3.840297 3.840297   ! a/sqrt(2) with a = 5.431 Angstrom
      ALPHA_BETA_GAMMA 60.0 60.0 60.0
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Si  0.00  0.00  0.00
      Si  0.25  0.25  0.25
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

运行与检查：

```bash
mpirun -np 8 cp2k.psmp -i si_primitive_kpoints.inp -o si_primitive_kpoints.out
grep "BRILLOUIN|" si_primitive_kpoints.out | head -12
```

`BRILLOUIN|` 块的前几行大致如下（示意输出，具体行随版本和设置而变）：

```text
 BRILLOUIN| K-point scheme                                          Monkhorst-Pack
 BRILLOUIN| K-Point grid                                               8    8    8
 BRILLOUIN| K-point symmetry accuracy                              0.10000000E-05
 BRILLOUIN| K-Point point group symmetrization                                OFF
 BRILLOUIN| Wavefunction type                                             COMPLEX
 BRILLOUIN| List of Kpoints [2 Pi/Bohr]                                       256
```

`List of Kpoints` 后的数字就是实际计算的 k 点数。8×8×8 偶数网格只做时间反演配对时为 8³/2 = 256，计算量大致与它成正比。随后是每个 k 点的权重和坐标列表。

### 7.3 k 点并行

k 点计算的 MPI 进程会被分成若干组，每组负责一部分 k 点。`PARALLEL_GROUP_SIZE` 的约束是：组大小能整除总进程数，组数能整除 k 点数。默认值 -1 会自动选择满足约束的最小组大小，k 点多时效率最好。若 k 点数是质数或与进程数不匹配，可以调整进程数，或者显式给出 `PARALLEL_GROUP_SIZE`。MPI 与 OpenMP 的配比见 [第18章](18-troubleshooting.md)。

---

## 8. k 点收敛测试（kpoint_convergence_template.inp + run_kpoint_scan.sh）

### 8.1 模板与脚本

模板与 7.2 节的输入基本相同，区别是网格由 `@SET NK` 控制、打开解析应力以便同时检查压力（所以 `RUN_TYPE` 必须是 `ENERGY_FORCE`），并关闭重启文件写出：

```text
# k-point convergence template: bulk Si primitive cell, MONKHORST-PACK NK x NK x NK
# run_kpoint_scan.sh rewrites the "@SET NK" line.
@SET NK 4

&GLOBAL
  PROJECT     si2_kconv
  RUN_TYPE    ENERGY_FORCE         ! ENERGY alone computes neither forces nor stress
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL         ! pressure converges more slowly than the energy
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     350
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK ${NK} ${NK} ${NK}
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-8
      MAX_SCF   100
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD  BROYDEN_MIXING
        ALPHA   0.4
        NBUFFER 8
      &END MIXING
      &PRINT
        &RESTART OFF
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
      ABC              3.840297 3.840297 3.840297
      ALPHA_BETA_GAMMA 60.0 60.0 60.0
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Si  0.00  0.00  0.00
      Si  0.25  0.25  0.25
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
  &PRINT
    &STRESS_TENSOR ON
    &END STRESS_TENSOR
  &END PRINT
&END FORCE_EVAL
```

```bash
bash run_kpoint_scan.sh                          # NK = 2 3 4 5 6 8 10 12
NKS="2 4 6 8 10" CP2K="mpirun -np 8 cp2k.psmp" bash run_kpoint_scan.sh
```

脚本为每个 NK 建立 `kp_<NK>/` 目录，运行后汇总到 `kpoint_scan.csv`。各列为：网格、实际 k 点数（`BRILLOUIN| List of Kpoints`）、总能量、相对最密网格的能量差（meV/原子）、压力（`STRESS| 1/3 Trace`，单位 bar）。本仓库的通用脚本 [conv_kpoints.py](../scripts/conv_kpoints.py) 能按晶格长度和 k 点密度自动生成网格序列（`python ../scripts/conv_kpoints.py -h`）。

### 8.2 判读

下面是一组**示意输出，数值仅供参考**：

```text
NK   nkpts  energy_Ha        dE_meV_per_atom  pressure_bar
2    4      -7.8730154215    109.6052         -21410.2
3    14     -7.8792338672    25.0010          -2104.7
4    32     -7.8806231904    6.0974           -545.8
5    63     -7.8809741118    1.3229           -150.3
6    108    -7.8810458392    0.3471           -58.9
8    256    -7.8810685021    0.0388           -14.2
10   500    -7.8810711297    0.0031           -9.8
12   864    -7.8810713573    0.0000           -9.6
```

（示意输出，数值仅供参考；k 点数按“偶数网格 n³/2、奇数网格 (n³+1)/2”的时间反演配对规则推算。）

判读要点：

1. 能量差在 6×6×6 时已低于 1 meV/原子；若要求 0.1 meV/原子，需要 8×8×8；
2. 压力收敛得更慢，到 10×10×10 才稳定在 1 kbar 以内。所以 CELL_OPT 要按压力选网格；
3. 偶数和奇数网格的收敛行为可能交替（是否包含 Γ 点），比较时最好只在同一类网格之间对比；
4. 半导体 Si 收敛较快；**金属**需要更密的网格，而且 k 点与展宽宽度要一起收敛（展宽越小，需要的 k 点越多）。

### 8.3 网格的经验起点

可以用“k 点长度” $R_k$ 估计每个方向的分割数：$n_i \approx \lceil R_k / a_i \rceil$，其中 $a_i$ 是该方向的晶格长度（Å）。下表只用来确定扫描起点，不能代替测试：

| 体系 | $R_k$（Å） | 例子 |
|---|---|---|
| 宽带隙绝缘体、分子晶体 | 15–25 | MgO 原胞（a ≈ 3.0 Å）：6×6×6 起 |
| 半导体 | 25–35 | Si 原胞（a ≈ 3.84 Å）：8×8×8 |
| 金属 | 40–60 | Al 原胞（a ≈ 2.86 Å）：14×14×14 起 |
| Slab（面内） | 同体相；法向取 1 | Pt(111) p(2×2)（a ≈ 5.6 Å）：金属 → 6×6×1 至 8×8×1 |

---

## 9. CP2K 2026.2 中的 k 点功能

### 9.1 新增功能（2026.2 release notes）

- **k 点对称约化**（`SYMMETRY T`，实验性）：在时间反演之外，再用原子结构的空间群操作合并等价 k 点。要求晶胞符合 CP2K 的下三角约定（用 `ABC` + `ALPHA_BETA_GAMMA` 或合适的 CIF 定义最稳妥），否则会警告并回退到只用时间反演；需要复数波函数；
- **k 点波函数外推**：GEO_OPT 和 MD 中，k 点计算现在也能用 `&QS EXTRAPOLATION` 加速每一步的 SCF（以前只有 Γ 点可以）；
- **k 点 Lowdin 布居分析**和 **k 点下的 Mulliken DFT+U**（`PLUS_U_METHOD MULLIKEN`），后者见 [第8章](08-dispersion-spin-dftu.md)；
- **带展宽和 k 点投影的 DOS/PDOS**（新的 `&DOS` 输入段），见 [第10章](10-dos-pdos-band.md)；
- k 点下的 Harris 泛函能量修正。

启用对称约化时，官方要求先与严格的完整网格做对比验证：

```text
&KPOINTS                           ! production candidate (experimental in 2026.2)
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY T
  WAVEFUNCTIONS COMPLEX
&END KPOINTS

&KPOINTS                           ! strict full-mesh reference for validation
  SCHEME MONKHORST-PACK 8 8 8
  SYMMETRY  F
  FULL_GRID T
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
```

能量、力、应力（以及你关心的其它性质）都一致后，才把约化版本用于生产计算。几何优化或 MD 中原子移动会改变对称性，程序会按当前构型重新判定，不可约 k 点集合可能随之变化；需要保持空间群时，配合 `KEEP_SPACE_GROUP` 使用（见 [第7章](07-geo-cell-opt.md)）。

### 9.2 k 点功能兼容性（2026.2）

| 功能 | 状态 |
|---|---|
| 标准对角化 | 支持 |
| OT | **不支持**（trunk 已支持，见 [第5章](05-scf-strategy.md)） |
| 原子对称约化 | 实验性，需对照验证 |
| 波函数外推 | 支持（2026.2 新增） |
| 其它对角化方法（Davidson、Krylov、OT 对角化） | 不支持 |
| 杂化泛函 | 通过 RI-HFXk 支持（见 [第9章](09-hybrid-admm.md)） |
| DFT+U | 有限：只支持 Mulliken |
| SCCS 隐式溶剂 | 未经验证，谨慎使用 |
| CDFT、XAS/RIXS、线性响应/DFPT、周期性电场 | 不支持 |
| TDDFPT | 有限：只支持 `KERNEL NONE` |
| GW | 走单独的工作流，不依赖 `&DFT/&KPOINTS` |
| Molden 输出 | k 点计算不可用（可改用 `.mokp` 输出） |

---

## 10. Γ 点还是 k 点？

| 情况 | 建议 |
|---|---|
| 分子、团簇（孤立） | 只用 Γ（不写 `&KPOINTS`），非周期 Poisson 求解器 |
| 液体、溶液、无定形（≥ 100 原子，盒长 ≥ 12–15 Å） | Γ 点 + OT，AIMD 的标准设置 |
| 大超胞绝缘体/半导体（各方向 ≥ 15–20 Å） | Γ 点通常足够；用一次 2×2×2 网格验证 |
| 小原胞半导体/绝缘体 | 必须收敛 k 网格（第 8 节） |
| 小原胞金属 | 必须收敛 k 网格，并与展宽一起收敛；对角化 + 展宽 |
| 金属 slab | 面内 k 点必须收敛，法向取 1；或者用足够大的表面超胞加少量 k 点 |
| 需要 OT 或 k 点不支持的功能（CDFT、XAS 等） | 扩大超胞，只用 Γ 点 |

> [!TIP]
> 判断“大超胞只用 Γ 点是否够”，最直接的方法是：同一个超胞分别用 Γ 和 MP 2×2×2 各算一次单点，比较每原子能量和最大力的差别。差别在你的精度要求之内，就可以放心只用 Γ。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 判据 |
|---|---|---|
| 分子盒子（MT） | ≥ 2 倍电子密度尺寸，常用 12–20 Å；`&CENTER_COORDINATES` | 盒长变化时能量 < 10⁻⁵ Ha |
| Slab 真空 | ≥ 15 Å；功函数计算 ≥ 20 Å | 真空区平面平均势出现平台；能量、受力随真空变化可忽略 |
| Slab 层数 | 金属 4–6 层，底部 1–2 层固定 | 表面能或吸附能随层数变化 < 目标精度 |
| 偶极校正 | 非对称 slab 必开：`SURFACE_DIPOLE_CORRECTION T`，`SURF_DIP_DIR Z` | 输出中的偶极矩稳定；不报真空不足 |
| k 网格（半导体） | $R_k$ ≈ 25–35 Å 起扫 | ΔE < 1 meV/原子；CELL_OPT 用 ΔP < 1 kbar |
| k 网格（金属） | $R_k$ ≈ 40–60 Å 起扫，与展宽一起收敛 | 同上，并检查费米能稳定 |
| 对称约化 | 2026.2 中默认关闭；开启前与 `FULL_GRID T` 对照 | 能量、力、应力一致 |
| 晶胞写法 | `ABC` + `ALPHA_BETA_GAMMA` + 分数坐标最不易出错 | `cp2k.psmp -e` 回显检查 |

---

## 常见问题与报错

**1. `An invalid set of cell vectors was specified. The cell volume is too small`**（或 `The determinant det(h) is too small`）
原因：`&CELL` 漏写、写错单位，`A`/`B`/`C` 中有向量为 0 或共线，或者 `CELL_FILE_NAME` 读到的晶胞为空。解决：检查 `ABC` 和 `ALPHA_BETA_GAMMA`（单位默认 Å 和度），用 `cp2k.psmp -e` 查看展开后的晶胞。

**2. `GEOMETRY wrong or EMAX_SPLINE too small`**
原因：原子重叠或间距过近，常见于周期像之间：xyz 不含晶胞信息、晶胞取得比结构小、扩胞或切面后出现重复原子、CIF 中的无序位点没有处理。解决：用可视化软件连同周期像检查结构和晶胞。

**2b. `The input parameters for SUBSYS%TOPOLOGY%MULTIPLE_UNIT_CELL and SUBSYS%CELL%MULTIPLE_UNIT_CELL have to agree.`**
原因：只在 `&TOPOLOGY` 中设了 `MULTIPLE_UNIT_CELL`，或两处数值不同。解决：两处写成相同的值（见 3.3 节）。反过来只在 `&CELL` 中设置时不会报错，但原子没有被复制，要检查输出中的原子数。

**3. `OT not possible with kpoint calculations`**
原因：2026.2 及以前，OT 不能与 `&KPOINTS` 共存（包括 `SCHEME GAMMA`）。解决：k 点计算改用对角化；大超胞只用 Γ 点时直接删掉 `&KPOINTS`。

**4. `Dipole correction needs more vacuum space above the surface`**
原因：偶极层位置的平面平均密度太大，真空不够或 slab 跨过了晶胞边界。解决：增大 c，把 slab 平移到晶胞中下部（不要跨越 z = 0 和 z = c 边界），必要时用 `SURF_DIP_POS` 指定真空中的位置。

**5. `Dipole correction only for surface perpendicular to one Cartesian axis`**
原因：C 轴不沿 `SURF_DIP_DIR`，或 A、B 在该方向有分量（例如从其它软件拷贝的倾斜晶胞）。解决：把表面晶胞转到“A、B 在 xy 平面、C 沿 z”，最简单的做法是用 `ABC` + `ALPHA_BETA_GAMMA 90 90 γ` 加分数坐标。

**6. MT 计算的能量随盒子大小明显变化，但没有任何报错。**
原因：盒子小于 2 倍电子密度尺寸，MT 结果不可靠。解决：按 5.2 节做盒子测试；分子较大时改用 WAVELET。

**7. 六方体系的结构“被扭曲”，或者对称性不对。**
原因：`SYMMETRY HEXAGONAL` 的含义是 γ = 60°，而坐标是按 γ = 120° 写的。解决：改成 `HEXAGONAL_GAMMA_120`，或删掉 `SYMMETRY`，并使用分数坐标。

**8. 改了 `&KPOINTS` 设置后重启失败、回退到 ATOMIC 初猜。**
原因：Γ 点（`.wfn`）和 k 点（`.kp`）的重启文件不能互换，加上或去掉 `&KPOINTS` 都会改变格式。`.kp` 重启对对称性设置也有要求：用 `FULL_GRID T` 写出的 `.kp` 可以用于 `FULL_GRID F` 的计算，反之亦然，但其它设置要一致。解决：格式不匹配时重新从 ATOMIC 起步，或者先用相同的 k 点设置跑一个便宜的计算。

**9. k 点计算并行效率很差。**
原因：k 点组的划分不合理（组数不能整除 k 点数，于是退化为较少的组）。解决：调整 MPI 进程数，使它与 k 点数有较多公因子，或者显式设置 `PARALLEL_GROUP_SIZE`。

**10. slab 计算的功函数、吸附能随真空厚度漂移。**
原因：非对称 slab 没有开偶极校正，或者真空太薄。解决：`SURFACE_DIPOLE_CORRECTION T`、`SURF_DIP_DIR Z`，真空 ≥ 15 Å，法向 k 点取 1，并做真空收敛。

---

## 练习

1. 用 `MULTIPLE_UNIT_CELL 2 2 2`（两处都设）把 Si 原胞扩成 16 原子超胞，只用 Γ 点计算，与原胞 2×2×2 Γ 中心网格（2026.2 可用 `GAMMA_CENTERED T`）的每原子能量比较，验证 2.1 节的等价关系。
2. 对 `h2o_isolated_mt.inp` 做 10–16 Å 的盒子测试；再改成 `PERIODIC XYZ` 加默认周期求解器重复一遍，比较两组曲线，体会周期像偶极相互作用的大小。
3. 关掉 `pt111_slab_dipole.inp` 中的 `SURFACE_DIPOLE_CORRECTION`，比较总能量和 O 原子受力；再把 O 原子删掉（对称的清洁 slab），重复比较。
4. 用 `run_kpoint_scan.sh` 完成 Si 原胞的 k 点收敛，分别找出能量和压力收敛到 1 meV/原子、1 kbar 所需的网格；如果你使用 2026.2，再对 8×8×8 打开 `SYMMETRY T`，与 `FULL_GRID T` 的结果对比，记录 k 点数和耗时。

---

## 参考资料

**官方文档**

- K-Points（2026.2 版）：<https://manual.cp2k.org/cp2k-2026_2-branch/methods/dft/k-points.html>
- Electrostatics and Poisson Solvers：<https://manual.cp2k.org/trunk/methods/dft/electrostatics/index.html>
- 输入参考 `&KPOINTS`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/KPOINTS.html>
- 输入参考 `&POISSON`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html>
- 输入参考 `&CELL`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/CELL.html>
- 输入参考 `&DFT`（`SURFACE_DIPOLE_CORRECTION` 等）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT.html>
- Restarting（`.wfn` 与 `.kp`）：<https://manual.cp2k.org/trunk/methods/restarting.html>
- CP2K 2026.2 release notes：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-kp-*` 等 regtest 是 k 点、对称约化和重启写法的权威参考；k 点支持状态跟踪见 issue #4854
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：中文 Wiki，含 slab 建模与偶极校正的实践经验
- [ishikawa-group/cp2k-tutorial](https://github.com/ishikawa-group/cp2k-tutorial)：日文入门教程，Pt 表面加苯分子的表面体系示例
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含周期 Si 章节
- [RobinLi-95/CP2K_tools](https://github.com/RobinLi-95/CP2K_tools)：V_HARTREE 平面平均、功函数、ghost layer slab 等脚本
- [funkymunkycool/Cube-Toolz](https://github.com/funkymunkycool/Cube-Toolz)：cube 文件平面平均与差分，可用于检查真空区电势
- [nanotech-empa/cp2k-spm-tools](https://github.com/nanotech-empa/cp2k-spm-tools)：表面体系的 STM/STS 模拟（进阶）
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：k 点收敛脚本思路（无 license，仅参考）

**关键文献**

- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* **13**, 5188 (1976)：<https://doi.org/10.1103/PhysRevB.13.5188>
- G. J. Martyna, M. E. Tuckerman, *J. Chem. Phys.* **110**, 2810 (1999)，MT 求解器：<https://doi.org/10.1063/1.477923>
- L. Genovese *et al.*, *J. Chem. Phys.* **125**, 074105 (2006)，自由边界小波 Poisson 求解器：<https://doi.org/10.1063/1.2335442>
- L. Genovese *et al.*, *J. Chem. Phys.* **127**, 054704 (2007)，表面边界条件的小波 Poisson 求解器：<https://doi.org/10.1063/1.2754685>
- J. Neugebauer, M. Scheffler, *Phys. Rev. B* **46**, 16067 (1992)，偶极校正：<https://doi.org/10.1103/PhysRevB.46.16067>
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999)，slab 偶极校正：<https://doi.org/10.1103/PhysRevB.59.12301>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](05-scf-strategy.md) | [返回目录](../README.md) | [下一章 →](07-geo-cell-opt.md)
