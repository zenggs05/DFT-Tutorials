# 第11章 实空间性质：cube 文件、电荷分析与功函数（Cube Files, Charge Analysis & Work Function in CP2K）

本章讲如何从 CP2K Quickstep 计算中导出实空间数据并做定量分析，内容包括：电子密度、Hartree 静电势、分子轨道和 ELF 的 cube 文件（`&E_DENSITY_CUBE`、`&V_HARTREE_CUBE`、`&MO_CUBES`、`&ELF_CUBE`）；Mulliken、Löwdin、Hirshfeld 原子电荷，以及借助外部程序做 Bader 分析；输出 molden 文件交给 Multiwfn 做波函数分析；由 `V_HARTREE` 平面平均计算表面功函数（work function，含偶极校正）；差分电荷密度（charge density difference）。最后简要介绍 CP2K 2026.2 新增的 3D 周期 slab 恒电势（巨正则）SCF。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 按需要打开各类 cube 输出，理解 `STRIDE`、`NHOMO`/`NLUMO`、`&EACH` 对文件大小和输出时机的影响；
- 读懂 cube 文件头，知道其中的单位（bohr、e/bohr³、Hartree）和数据排列顺序；
- 根据问题选择电荷分析方法，并清楚各方法的局限；
- 为 Multiwfn 生成可用的 molden 文件，知道 2026.2 新增的 `WRITE_CELL`、`WRITE_PSEUDO`、`NLUMO` 解决了哪些老问题；
- 用偶极校正后的 Hartree 势计算 slab 的功函数，并判断真空层是否够厚；
- 正确地做差分电荷密度：同一网格、同一坐标、同一参数。

| 想做的事 | 用什么 | 本章小节 |
|---|---|---|
| 画电子密度、自旋密度等值面 | `&E_DENSITY_CUBE` | 3 |
| 画 HOMO/LUMO 等轨道 | `&MO_CUBES`（Γ 点） | 3 |
| 原子电荷、电荷转移 | Mulliken / Löwdin / Hirshfeld（输出中直接给出）；Bader（外部程序读密度 cube） | 4 |
| ADCH、CM5、AIM、IGMH、键级等波函数分析 | `&MO_MOLDEN` + Multiwfn | 5 |
| 功函数、表面偶极、能级对齐 | `&V_HARTREE_CUBE` + `SURFACE_DIPOLE_CORRECTION` | 6 |
| 吸附、成键引起的电荷重新分布 | 三个 `E_DENSITY` cube 相减 | 7 |
| 给定电极电势下的表面计算 | 2026.2 的 `&SCF/&GCE` | 8 |

> [!NOTE]
> 示例输入集在 [../examples/11-cube-charge-workfunction/README.md](../examples/11-cube-charge-workfunction/README.md)。周期 slab 的建模、真空层和 k 点设置见 [第6章](06-periodic-surface-kpoints.md)；DOS/PDOS 见 [第10章](10-dos-pdos-band.md)；通用可视化（VESTA、Multiwfn、出版级作图）见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md)。

---

## 2. 原理速览

### 2.1 cube 文件：CP2K 写出的是什么

cube 是 Gaussian 定义的体数据格式，CP2K 沿用了它：

```text
-Quickstep-                                        <- comment line 1
 ELECTRON DENSITY                                  <- comment line 2 (what is stored)
    3    0.000000    0.000000    0.000000          <- number of atoms, origin (bohr)
   80    0.283459    0.000000    0.000000          <- n1 and voxel vector along a (bohr)
   80    0.000000    0.283459    0.000000          <- n2 and voxel vector along b
   80    0.000000    0.000000    0.283459          <- n3 and voxel vector along c
    8    6.000000   11.338    11.338    11.555     <- Z, valence charge, x, y, z (bohr)
    1    1.000000   ...
  1.23456E-05  2.34567E-05 ...                     <- data, 6 per line, z fastest
```

（示意，数值仅供参考。）要点：

- 长度单位是 **bohr**（1 bohr = 0.529177 Å），网格可以是斜交的（voxel 矢量平行于晶胞矢量）；
- 数据顺序是 x 最慢、z 最快：对每个 (ix, iy) 连续写出所有 iz，每 6 个数换行。Python 中读成一维数组后 `reshape(n1, n2, n3)` 即可；
- 原子行的第 2 列是 CP2K 写入的**有效核电荷（价电子数）**，不是原子序数；
- 密度只含价电子（GTH 赝势）。在 GPW 中原子核附近的密度是平滑的，这对可视化和 Bader 分析都可以接受，但不能拿它谈“芯区”性质。全电子密度需要 GAPW（见 [第3章](03-gpw-basis-pseudo.md)）。

### 2.2 网格密度与文件大小

cube 的网格就是 CP2K 最细的实空间网格（由 `&MGRID/CUTOFF` 和晶胞决定），再按 `STRIDE` 抽稀。网格间距约为 $h\approx\pi/\sqrt{E_\mathrm{cut}[\mathrm{Ry}]}$ bohr：`CUTOFF 400` 时约 0.157 bohr。12 Å（22.7 bohr）立方盒子每个方向约 145 个点（CP2K 会再取一个便于 FFT 的整数），`STRIDE 1` 时约 300 万个数据点；每个数占 13 个字符，一个 cube 约 40 MB。`STRIDE 2 2 2` 把点数减到 1/8（约 5 MB）。以上是粗略估算。

| 用途 | 建议 `STRIDE` |
|---|---|
| 只看等值面、做图 | 2 2 2（默认）或更大 |
| Bader 分析、电荷积分、平面平均求功函数、差分密度 | 1 1 1 |
| 大体系、MD 中每隔若干步输出 | 3 3 3 以上，或用 `MAX_FILE_SIZE_MB` 限制（见 3.2 节） |

### 2.3 原子电荷：没有“唯一正确”的定义

原子电荷不是量子力学可观测量，不同划分方法的数值可以差很多：

| 方法 | 划分对象 | 优点 | 局限 |
|---|---|---|---|
| Mulliken | 基函数空间，重叠布居平分 | 计算零成本，输出默认就有 | 强烈依赖基组，弥散函数多时数值失真甚至无意义 |
| Löwdin | 对称正交化后的基函数空间 | 比 Mulliken 稳定一些 | 仍依赖基组 |
| Hirshfeld | 实空间，按孤立原子密度比例划分 | 对基组不敏感，数值偏小 | 电荷转移被系统性低估；Hirshfeld-I（自洽迭代）可部分改善 |
| Bader（AIM） | 实空间，按密度梯度零通量面划分 | 物理图像清晰，适合比较电荷转移趋势 | 需要细网格；只有价电子密度时，零通量面在原子核附近可能不够准；数值偏大 |
| ADCH、CM5 等 | 基于 Hirshfeld 的修正 | 能较好重现偶极矩 | 需要 Multiwfn 读 molden 文件 |

实用原则：**比较趋势时只用一种方法，并在同一基组、同一泛函下比较**；论文中写清楚用的是哪种电荷。

### 2.4 功函数与偶极校正

功函数定义为把电子从费米能级移到表面外真空中所需的最小能量：

$$
\Phi=V_\mathrm{vac}-E_F
$$

其中 $V_\mathrm{vac}$ 是真空区中电子静电势能的平台值，$E_F$ 是同一计算中的费米能级。周期 slab 计算中，静电势的零点由周期性边界条件决定（平均静电势为零），**$E_F$ 本身的绝对值没有意义，只有 $V_\mathrm{vac}-E_F$ 有意义**。

不对称 slab（单面吸附、极性表面）在法向有净偶极，周期性重复会在真空中产生人为电场，两侧真空能级不再平直。偶极校正（Bengtsson 1999）在真空中加一个锯齿形的补偿势，抵消人为电场，使两侧各自出现平台；两侧平台之差就是表面偶极引起的势阶。CP2K 中用 `SURFACE_DIPOLE_CORRECTION T` + `SURF_DIP_DIR Z` 打开，补偿势会加到 Hartree 势上，因此写出的 `V_HARTREE` cube 已经包含了校正。

### 2.5 差分电荷密度

$$
\Delta\rho(\mathbf r)=\rho_{AB}(\mathbf r)-\rho_{A}(\mathbf r)-\rho_{B}(\mathbf r)
$$

$A$、$B$ 两部分的原子位置必须与 $AB$ 中**完全一致**（直接从 AB 结构中删掉另一部分，不做优化），三次计算的晶胞、`CUTOFF`、`REL_CUTOFF` 和 `STRIDE` 也必须完全一致，这样网格点才能一一对应。沿法向的平面平均 $\Delta\bar\rho(z)$ 及其积分

$$
\Delta Q(z)=\int_{0}^{z}\Delta\bar\rho(z')\,dz'
$$

给出电荷沿 z 的重新分布（charge displacement curve）。

---

## 3. 实战一：H2O 的密度、Hartree 势、轨道和 ELF cube

### 3.1 输入文件

`h2o_density_mo_cubes.inp` 的 `&PRINT` 部分（其余与[第10章](10-dos-pdos-band.md)的 H2O 输入相同，OT + MT 求解器）：

```text
    &PRINT
      ! valence electron density (e/bohr^3); STRIDE 1 = full grid, needed for Bader
      &E_DENSITY_CUBE
        STRIDE 1 1 1
      &END E_DENSITY_CUBE
      ! Hartree potential of electrons + ions (hartree), electron-energy sign convention
      &V_HARTREE_CUBE
        STRIDE 2 2 2
      &END V_HARTREE_CUBE
      ! HOMO-1, HOMO, LUMO, LUMO+1 (with OT the LUMOs are computed after the SCF)
      &MO_CUBES
        NHOMO  2
        NLUMO  2
        STRIDE 2 2 2
      &END MO_CUBES
      ! electron localization function, values between 0 and 1
      &ELF_CUBE
        STRIDE 2 2 2
      &END ELF_CUBE
      &MULLIKEN ON
      &END MULLIKEN
      &LOWDIN ON
      &END LOWDIN
      &HIRSHFELD ON
      &END HIRSHFELD
      ! dipole moment; PERIODIC F because the cell is not periodic
      &MOMENTS ON
        PERIODIC F
      &END MOMENTS
    &END PRINT
```

这些打印键的默认级别多数是 `HIGH` 或更高（`E_DENSITY_CUBE`、`V_HARTREE_CUBE`、`MO_CUBES`、`ELF_CUBE`、`TOT_DENSITY_CUBE`、`LOWDIN`、`MOMENTS` 都是 HIGH），在默认的 `PRINT_LEVEL MEDIUM` 下不会输出，**写出该段即开启**。`MULLIKEN` 和 `HIRSHFELD` 的默认级别是 MEDIUM，所以 `PRINT_LEVEL MEDIUM` 时本来就会打印；在 `PRINT_LEVEL LOW` 的 MD 或大体系计算中，要像上面这样显式写 `ON`。

### 3.2 各 cube 段的关键字

| 段 / 关键字 | 默认值 | 说明 |
|---|---|---|
| `&E_DENSITY_CUBE/STRIDE` | 2 2 2 | 可给 3 个数（x、y、z）或 1 个数 |
| `&E_DENSITY_CUBE/DENSITY_INCLUDE` | `TOTAL_HARD_APPROX` | 只对 GAPW 有意义：`TOTAL_DENSITY` 输出含硬密度的总密度，`SOFT_DENSITY` 只输出软部分 |
| `&V_HARTREE_CUBE` | — | 由电子 + 离子总电荷产生的静电势，单位 Hartree。**符号约定与物理静电势相反**，即它是电子在该处的静电势能（与本征值同一约定）。只适用于 GPW |
| `&MO_CUBES/NHOMO`、`NLUMO` | 1、0 | 输出最高的 NHOMO 个占据轨道、最低的 NLUMO 个空轨道（−1 = 全部）。OT 计算中 LUMO 在 SCF 后单独求解 |
| `&MO_CUBES/HOMO_LIST` | — | 按编号指定要输出的占据轨道，优先于 `NHOMO` |
| `&MO_CUBES/WRITE_CUBE` | T | 设为 F 时不写 cube，只计算并打印本征值和 `HOMO - LUMO gap [eV]`（OT 计算中拿 LUMO 能量的便捷办法） |
| `MAX_FILE_SIZE_MB` | 0（不限） | 自动选 stride 使文件不超过给定大小。`MO_CUBES` 自 2026.1 起可用；2026.2 中 `V_HARTREE_CUBE`、`TOT_DENSITY_CUBE` 也有；`E_DENSITY_CUBE` 和 `ELF_CUBE` 没有 |
| `&ELF_CUBE/DENSITY_CUTOFF` | 1.0E-10 | 密度低于此值的点不计算 ELF（真空区） |
| `&TOT_DENSITY_CUBE` | — | 电子 + 原子芯的总电荷密度（电子为正、芯为负） |
| `APPEND` | F | 每次输出时追加到已有文件（MD 中做“cube 轨迹”） |
| `&EACH` 子段 | — | 控制输出频率，见 3.4 节 |

### 3.3 运行与输出

```bash
cd CP2K/examples/11-cube-charge-workfunction
mpirun -np 4 cp2k.psmp -i h2o_density_mo_cubes.inp -o h2o_density_mo_cubes.out
ls -1 h2o_cube-*.cube
grep "HOMO - LUMO gap" h2o_density_mo_cubes.out
```

文件命名规则为 `<PROJECT>-<中间名>-<迭代编号>.cube`，单点计算的迭代编号通常是 `1_0`：

| 文件 | 内容 |
|---|---|
| `h2o_cube-ELECTRON_DENSITY-1_0.cube` | 价电子密度（e/bohr³）；UKS 计算还会另写 `...-SPIN_DENSITY-1_0.cube`（α − β） |
| `h2o_cube-v_hartree-1_0.cube` | Hartree 势（Hartree） |
| `h2o_cube-WFN_00004_1-1_0.cube` | 第 4 个 MO（H2O 的 HOMO），`_1` 表示自旋通道 1 |
| `h2o_cube-ELF_S1-1_0.cube` | ELF，`S1` 表示自旋通道 1 |

检查要点：

- 对密度 cube 积分（数据之和 × voxel 体积）应得到价电子数 8；用本仓库的 `python ../scripts/cube_tools.py -h` 查看积分功能，或用第 7 节的 `chgdiff.py` 读取方式自己算；
- VESTA 打开 `WFN_00004_1` 取 ±0.05（单位 bohr$^{-3/2}$）左右的等值面，应看到垂直于分子平面（yz 平面）的 O 2p 型轨道（1b1）；轨道 cube 存的是波函数值（有正负），不是密度；
- ELF 在 O–H 键中点和 O 的孤对电子区接近 1。

> [!TIP]
> 想要在 VMD 中显示 cube：`vmd h2o_cube-ELECTRON_DENSITY-1_0.cube`，在 Graphics → Representations 中选 Isosurface。Multiwfn 也能直接读 cube 做格点数据的平面图、积分和差值运算。

### 3.4 在优化或 MD 中控制输出时机

GEO_OPT 中，cube 默认会在每一步优化后输出，既浪费磁盘又容易覆盖。只要最终结构的 cube 时：

```text
      &E_DENSITY_CUBE
        STRIDE 2 2 2
        ADD_LAST NUMERIC         ! also print at the last step, file name gets the step number
        &EACH
          GEO_OPT 0              ! never during the optimization ...
        &END EACH
      &END E_DENSITY_CUBE
```

`ADD_LAST NUMERIC` 让最后一步一定输出，`&EACH GEO_OPT 0` 关闭中间步。MD 中可以用 `&EACH MD 50` 每 50 步输出一次，配合 `APPEND` 写成多帧文件（VMD 可读多帧 cube 的轨迹，但文件会非常大）。

---

## 4. 电荷分析：Mulliken、Löwdin、Hirshfeld 与 Bader

### 4.1 Mulliken 与 Löwdin

`PRINT_LEVEL MEDIUM` 时 Mulliken 布居默认打印。表头（RKS）为：

```text
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
```

UKS 时变为 `Atomic population (alpha,beta) Net charge  Spin moment`，最后一行是 `# Total charge and spin`。`Spin moment` 一列是判断磁性体系（如反铁磁 NiO 的两套 Ni 子晶格）自旋是否收敛到预期状态的最快方法。

| 关键字 | 说明 |
|---|---|
| `&MULLIKEN/PRINT_GOP` | 额外打印每个基函数的轨道布居（gross orbital population） |
| `&MULLIKEN/PRINT_ALL` | 打印完整的净 AO 和重叠布居矩阵（大体系慎用） |
| `&LOWDIN`（`PRINT_GOP`、`PRINT_ALL`） | Löwdin 布居，默认级别 HIGH，需显式打开；2026.2 起支持 k 点 |

### 4.2 Hirshfeld 与 Hirshfeld-I

```text
      &HIRSHFELD ON
        SHAPE_FUNCTION   GAUSSIAN   ! or DENSITY (atomic densities as shape functions)
        REFERENCE_CHARGE ATOMIC     ! or MULLIKEN
        SELF_CONSISTENT  F          ! T = iterative Hirshfeld-I
      &END HIRSHFELD
```

| 关键字 | 默认 | 说明 |
|---|---|---|
| `SHAPE_FUNCTION` | `GAUSSIAN` | 原子划分函数的类型：球形高斯，或 `DENSITY`（孤立原子密度） |
| `REFERENCE_CHARGE` | `ATOMIC` | 划分函数对应的原子电荷：中性原子，或用 Mulliken 电荷 |
| `SELF_CONSISTENT` | F | T 时做 Hirshfeld-I 自洽迭代（只缩放整个形状函数） |
| `USER_RADIUS` + `ATOMIC_RADII` | F | 自定义高斯的原子半径（Å） |

输出表头为 `Hirshfeld Charges`，列为 `#Atom  Element  Kind  Ref Charge  Population  Net charge`（UKS 时多一列 `Spin moment`），最后一行是 `Total Charge`。

```bash
grep -A6 "Mulliken Population Analysis" h2o_density_mo_cubes.out
grep -A6 "Hirshfeld Charges"           h2o_density_mo_cubes.out
```

### 4.3 Bader 分析：用外部程序读密度 cube

CP2K 本身不做 Bader 划分，常用 Henkelman 组的 `bader` 程序（<https://theory.cm.utexas.edu/henkelman/code/bader/>）读取 `E_DENSITY_CUBE`：

```bash
# 1) the density cube must be written with STRIDE 1 1 1 (full grid)
bader h2o_cube-ELECTRON_DENSITY-1_0.cube
# 2) per-atom electron populations are in ACF.dat (column CHARGE)
cat ACF.dat
```

注意：

- `ACF.dat` 中的 `CHARGE` 是划分到每个原子盆的**电子数**，原子净电荷 = 该原子的价电子数（赝势的 q 值，如 O 为 6、H 为 1、Pt 为 18）− `CHARGE`；cube 原子行第 2 列就是这个价电子数；
- Bader 结果对网格很敏感：要用 `STRIDE 1`，并检查提高 `CUTOFF` 后电荷是否稳定（变化 < 0.01 e）；
- 价电子密度在原子核附近没有尖峰，个别情况下会出现“非核吸引子”或盆划分异常；可以先用 Multiwfn 的 AIM 模块（读 molden，见第 5 节）交叉检验；
- 周期体系、斜交晶胞的 cube 也能直接用 `bader` 读取。

### 4.4 其他实空间方法：Voronoi 与偶极矩

- `&PRINT/&VORONOI`：基于 LibVori 的 Voronoi 积分（需编译时链接 libvori），可输出原子电荷和分子电多极矩，适合 MD 轨迹（`MOLECULAR_PROPERTIES`、`OUTPUT_EMP`）；
- `&PRINT/&MOMENTS`：体系偶极矩。孤立体系用 `PERIODIC F`；周期体系默认 `PERIODIC T`，按 Berry 相位公式计算，结果只确定到一个晶格矢量量子的整数倍，MD 中可能出现跳变（官方说明），做 IR 谱时要做去跳处理（见 [第12章](12-aimd.md)）。

```bash
grep -A2 "Dipole moment \[Debye\]" h2o_density_mo_cubes.out
```

---

## 5. 输出 molden 文件给 Multiwfn

Multiwfn 读 CP2K 的 molden 文件后，可以做 ADCH/CM5/RESP 电荷、AIM 拓扑、IGMH/IRI 弱相互作用、键级、轨道成分、DOS 等大量分析。sobereva 在《详谈使用 CP2K 产生给 Multiwfn 用的 molden 格式的波函数文件》（<http://sobereva.com/651>，2026-11 更新）中详细说明了各版本的做法，下面按 2026.2 总结。

### 5.1 输入

`h2o_molden_for_multiwfn.inp` 的 `&PRINT` 部分：

```text
    &PRINT
      &MO_MOLDEN
        NDIGITS      9            ! keep coefficients down to 1E-9 (default 3 is for pictures only)
        GTO_KIND     SPHERICAL
        NLUMO        10           ! >= 2026.2: 10 virtual orbitals (OT computes them after SCF)
        WRITE_PSEUDO T            ! >= 2026.2: [Pseudo] block with valence charges
        ! WRITE_CELL T            ! >= 2026.2: add [Cell] block, use for PERIODIC systems
      &END MO_MOLDEN
    &END PRINT
```

| 关键字 | 默认 | 版本 | 说明 |
|---|---|---|---|
| `NDIGITS` | 3 | 全部 | 保留的有效位数；绝对值小于 $10^{-\mathrm{NDIGITS}}$ 的系数不写。默认 3 只够画图，**做波函数分析用 9** |
| `GTO_KIND` | `SPHERICAL` | 全部 | 球谐或笛卡尔 GTO |
| `NLUMO` | 0 | ≥ 2026.2 | 写入的空轨道数。OT 计算在 SCF 后补算；对角化计算中不足时自动加大 `ADDED_MOS` |
| `WRITE_CELL` | F | ≥ 2026.2 | 写 `[Cell]` 段。周期体系必须提供晶胞信息，否则 Multiwfn 会把它当孤立体系分析 |
| `WRITE_PSEUDO` | F | ≥ 2026.2 | 写 `[Pseudo]` 段，记录每个原子的有效核电荷；Multiwfn 据此做原子电荷分析，并自动补上芯电子密度（EDFlib） |
| `UNIT` | `BOHR` | ≥ 2026.2 | 坐标与晶胞的单位（`BOHR`/`ANGSTROM`） |
| `MARK_GHOST` | F | ≥ 2026.2 | 把 ghost 原子的原子序数写成 0 |

输出文件为 `h2o_molden-MOS-1_0.molden`（`FILENAME` 默认 `MOS`）。

### 5.2 旧版本（≤ 2026.1）的变通办法

在 2026.1 及以前，`&MO_MOLDEN` 只有 `NDIGITS`、`GTO_KIND` 等少数关键字，需要手工处理（详见 sobereva.com/651）：

- **空轨道**：用对角化 SCF 并设 `ADDED_MOS`；OT 计算写出的 molden 只有占据轨道，且没有轨道能量。常用做法是先用 OT 收敛，再以 `SCF_GUESS RESTART` 读入 `RESTART.wfn` 做一次对角化单点并写 molden；
- **晶胞**：在 molden 文件开头手工插入 `[Cell]` 段（三行晶胞矢量，Å），或在同目录放 `[Cell].txt` 让 Multiwfn 读取；
- **价电子数**：手工修改 `[Atoms]` 段第 3 列为有效核电荷，或插入 `[Nval]` 段。

### 5.3 注意事项

- **k 点计算无法输出 molden**（官方 K-Points 文档：Molden output is not available for k-point calculations），周期体系要用 Γ 点超胞；
- 半经验方法（xTB、DFTB）不能输出 molden；
- MOLOPT 基组是广义收缩基，原函数很多，在 Multiwfn 中做实空间函数（AIM、IGMH、ELF）分析耗时较高；Mulliken 型分析对弥散函数敏感，sobereva 建议这类分析用 MOLOPT-SR 等弥散较弱的基组；
- 用 Multiwfn 生成 CP2K 输入（主菜单输入 `cp2k`，见 <http://sobereva.com/587>）时，它会根据你的选择自动加上 `&MO_MOLDEN` 等设置；生成后先用 `cp2k.psmp -c x.inp` 按当前版本检查，详见 [第17章](17-tools-ecosystem.md)。使用 Multiwfn 的结果发表时请引用 Multiwfn 的论文（见参考资料）。

---

## 6. 实战二：Pt(111) 表面的功函数

### 6.1 输入文件

`pt111_workfunction.inp`：

```text
&GLOBAL
  PROJECT     pt111_wf
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    SURFACE_DIPOLE_CORRECTION T
    SURF_DIP_DIR Z
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 6 6 1
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   300
      ADDED_MOS 50                ! 144 occupied MOs (16 Pt x 18 e / 2)
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.1               ! small mixing against charge sloshing in slabs
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
    &PRINT
      ! full grid along z gives the smoothest planar average
      &V_HARTREE_CUBE
        STRIDE 1 1 1
      &END V_HARTREE_CUBE
      &E_DENSITY_CUBE
        STRIDE 2 2 2
      &END E_DENSITY_CUBE
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 5.5494 5.5494 30.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   pt111_2x2_4L.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Pt
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
      POTENTIAL GTH-PBE-q18
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 关键字 | 本例取值 | 说明 |
|---|---|---|
| `SURFACE_DIPOLE_CORRECTION` | T | 开启表面偶极校正（别名 `SURF_DIP`）。只实现了法向平行于某个笛卡尔轴的 slab |
| `SURF_DIP_DIR` | Z | 表面法向 |
| `SURF_DIP_POS` | 默认 −1 | 补偿势的位置（Å，沿法向）；默认放在网格上电子密度最小处，一般不用改 |
| 晶胞 | 5.549 × 5.549 × 30 Å，γ = 60° | p(2×2) 表面胞，4 层 slab 厚约 6.8 Å，真空约 23 Å。c 轴平行于 z，平面平均只需对 xy 求平均 |
| `&KPOINTS 6 6 1` | — | 真空方向 1 个 k 点；k 点计算在 2026.2 中必须用对角化 |
| `&SMEAR` + `ADDED_MOS 50` | 500 K | 金属表面需要展宽，`ADDED_MOS` 约为占据态数的 1/3，覆盖费米能以上足够多的态 |
| `ALPHA 0.1` | — | slab 的长真空方向容易出现电荷“晃荡”（charge sloshing），用较小的混合系数 |
| `&V_HARTREE_CUBE STRIDE 1 1 1` | — | 平面平均对 z 方向的分辨率最敏感，用完整网格 |

> [!IMPORTANT]
> 结构中用的是 Pt 的实验晶格常数 3.924 Å。正式计算应使用同一泛函、同一基组优化得到的体相晶格常数来切表面，否则 slab 内存在应变，会使功函数产生系统偏差。表层原子的弛豫也会改变功函数（通常是零点几 eV 的量级以内），需要先做 GEO_OPT 并固定底部一到两层（见 [第7章](07-geo-cell-opt.md)）。

### 6.2 运行与求功函数

```bash
cd CP2K/examples/11-cube-charge-workfunction
sbatch job_slurm.sh pt111_workfunction.inp           # or: mpirun -np 16 cp2k.psmp -i ... -o ...
grep "Fermi energy:" pt111_workfunction.out | tail -n 1
python workfunction.py pt111_wf-v_hartree-1_0.cube pt111_workfunction.out
```

`workfunction.py`（示例目录中，只需 numpy）的做法是：读 cube，对每个 z 层的 xy 网格求平均并乘以 27.2114 换算成 eV；从输出最后一个 `Fermi energy:` 行读取费米能（Hartree）；在 slab 上下表面外 4–6 Å 处各取一个 2 Å 的窗口求平均作为两侧真空能级，并给出窗口内的波动（spread）。核心只有几行：

```python
v = values.reshape(n1, n2, n3)               # cube data, hartree
vz = v.mean(axis=(0, 1)) * 27.211386         # planar average along z, eV
phi = vz[vacuum_window].mean() - e_fermi_ev  # work function, eV
```

本仓库还有两个可用的工具：`CP2K/scripts/cube_tools.py`（平面平均、差分、积分；`python ../scripts/cube_tools.py -h` 查看用法），以及仓库作者收录的 `Scripts/analysis/extract_vacuum_level_cube.py`（读 `<项目名>-v_hartree-1_0.cube`，取平面平均势的最大值作为真空能级，说明见 [Scripts/README.md](../../Scripts/README.md) 第 9.1 节）。

### 6.3 如何判断结果可信

1. **平台是否平直**：把 `vh_planar.dat` 的第 1、3 列画出来（$V(z)-E_F$），真空区应有明显的平台，窗口内的 spread 应在 0.01 eV 量级。若平台倾斜，说明真空太薄，或不对称 slab 没开偶极校正；
2. **两侧是否一致**：对称的干净 slab 两侧平台应相同；有吸附物的一侧会因表面偶极而不同，两侧各自对应一个功函数；
3. **收敛性**：分别测试真空厚度（15、20、25 Å）、层数（4、5、6 层）和面内 k 点（6×6、9×9、12×12），功函数变化应 < 0.05 eV；
4. **符号**：CP2K 的 `V_HARTREE` 是电子的静电势能，所以真空平台的数值**高于** $E_F$，$\Phi$ 为正。若算出负值，检查是否误把势取了负号，或者单位没有统一（cube 是 Hartree，`Fermi energy:` 也是 Hartree）；
5. **与文献比较**：Pt(111) 的实验功函数约 5.9–6.1 eV（不同测量差异较大，见 Derry 等 2015 年的推荐值综述）；PBE 一般低估零点几 eV。以上只是数量级参考，以你的收敛结果为准。

> [!WARNING]
> 若报错 `Dipole correction needs more vacuum space above the surface`，说明 CP2K 找到的“最低密度层”上电子密度仍不小于 1E-5，即真空层太薄，或者 slab 跨越了晶胞边界、两侧分子挤在一起。加厚真空，或用 `SURF_DIP_POS` 指定补偿势放置的位置。

### 6.4 另一种做法：二维周期 Poisson 求解器

把 `&CELL PERIODIC XY` 与 `&POISSON PERIODIC XY` 配合使用支持二维边界的求解器（如 `ANALYTIC`、`MT`、`WAVELET`），可以从根本上去掉 z 方向的周期像，不再需要偶极校正。代价是对盒子尺寸和求解器有额外要求，与 k 点、应力等功能的兼容性也要逐项确认（见 [第6章](06-periodic-surface-kpoints.md)）。大多数表面研究仍采用 3D 周期 + 偶极校正。

---

## 7. 实战三：CO/Cu(111) 的差分电荷密度

### 7.1 工作流

示例目录 `co_on_cu_charge_diff/`（说明见 [../examples/11-cube-charge-workfunction/co_on_cu_charge_diff/README.md](../examples/11-cube-charge-workfunction/co_on_cu_charge_diff/README.md)）包含三个输入，它们通过 `@INCLUDE dft_common.inc` 共用同一套 `&DFT` 设置：

```text
# full.inp (system AB); slab.inp and co.inp differ only in PROJECT,
# the coordinate file, the &KIND list and the value of NADD
@SET NADD 40

&GLOBAL
  PROJECT     full
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
@INCLUDE dft_common.inc
  &END DFT
  &SUBSYS
    &CELL
      ABC 5.1124 5.1124 28.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   cu111_co_full.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Cu
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
    ...
```

`dft_common.inc` 中固定了 `CUTOFF 500`、`REL_CUTOFF 60`、4×4×1 k 点、偶极校正、Fermi–Dirac 展宽和 `&E_DENSITY_CUBE STRIDE 1 1 1`，并用 `ADDED_MOS ${NADD}` 接收各输入的 `@SET` 变量（CO 单独计算时只有 26 个基函数、5 个占据态，取 10；含 Cu 的体系取 40）。

```bash
cd CP2K/examples/11-cube-charge-workfunction/co_on_cu_charge_diff
cp2k.psmp -c full.inp       # syntax check, expands @INCLUDE
bash run_all.sh             # three single points + chgdiff.py
```

`chgdiff.py` 读入三个 cube，检查三者的格点数和 voxel 矢量是否一致（不一致直接退出），写出 `diff.cube` 和平面平均文件 `diff_planar.dat`（列：z [Å]、$\Delta\bar\rho(z)$ [e/Å]、$\Delta Q(z)$ [e]），并打印 $\int\Delta\rho\,d\mathbf r$ 作为自检（应接近 0）。

### 7.2 结果解读

- **三维图**：在 VESTA 中打开 `diff.cube`，正、负等值面（如 ±0.003 e/bohr³）分别表示电子聚集和电子流失。对 CO/金属体系，常见的图像是 C 端 5σ 区域电子减少（σ 给予）、CO 的 2π* 区域和 Cu–C 之间电子增加（π 反馈），具体以计算为准；
- **一维曲线**：$\Delta\bar\rho(z)$ 显示电荷沿法向的移动；$\Delta Q(z)$ 在某个 z 处的值等于该平面以下净增加的电子数，常取吸附物与表面之间 $\Delta Q$ 的极值作为“电荷转移量”的一种估计；
- **与原子电荷对比**：同时比较 AB 与 A、B 中 CO 的 Bader 或 Hirshfeld 电荷差（`dft_common.inc` 已开 `&HIRSHFELD`），两种方法得到的转移量不会相同，但符号和趋势应一致。

> [!CAUTION]
> 三个最常见的错误：(1) 对 A、B 分别做了几何优化，坐标不再与 AB 重合；(2) 在某个输入中用了 `&TOPOLOGY/&CENTER_COORDINATES`，把子体系重新居中了；(3) 三个计算的 `CUTOFF`、晶胞或 `STRIDE` 不同，网格对不上。后两种情况，`chgdiff.py` 会检测到网格不一致而退出；第一种不会报错，但结果毫无意义。

---

## 8. 2026.2 新功能：3D 周期 slab 的恒电势（巨正则）SCF

电化学界面上，电极电势（即功函数）由外部电源固定，电极上的净电荷随吸附状态改变。CP2K 2026.2 新增了面向 slab 模型的巨正则 SCF（PR #5467），在 SCF 迭代中调节体系的电子数，使功函数趋向用户给定的目标值；同时新增了平面反电荷模型（planar counter charge，PR #5064）来屏蔽 slab 的净电荷，以及配套的新 Pulay 混合（PR #5407）。

根据 2026.2 的输入说明，目前的实现有以下限制：

- 只用于**三维周期**计算，且必须配合 `&DFT/&PLANAR_COUNTER_CHARGE`；
- 需要 Fermi–Dirac 展宽和标准对角化；
- 目标功函数是针对**对称 slab** 模型定义的。

一个最小的设置骨架如下（关键字取自 2026.2 输入参考，数值仅为示意，需要按体系调整；可对照官方 regtest `tests/QS/regtest-gcscf/`）：

```text
  &DFT
    ...
    &PLANAR_COUNTER_CHARGE ON
      PARALLEL_PLANE XY        ! surface parallel to the xy plane
      DIST_EDGE      [angstrom] 1.5   ! distance of the counter-charge Gaussian from the cell boundary
      GAU_C          [angstrom] 1.5   ! spread of the Gaussian
    &END PLANAR_COUNTER_CHARGE
    &SCF
      ADDED_MOS 40
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &GCE ON
        TARGET_WORKFUNCTION [eV] 4.5   ! default 0.16 hartree
        MIXING_COEF 0.3                ! fraction of the target mixed in each step
      &END GCE
    &END SCF
  &END DFT
```

| 关键字 | 默认 | 说明 |
|---|---|---|
| `&SCF/&GCE` | 关闭 | 打开巨正则 SCF |
| `TARGET_WORKFUNCTION` | 0.16 Hartree | 目标功函数（可写 `[eV]`） |
| `MIXING_COEF` | 0.3 | 每步把目标功函数按此比例混入上一步的功函数 |
| `&PLANAR_COUNTER_CHARGE/PARALLEL_PLANE` | `XY` | 表面所在平面 |
| `DIST_EDGE`、`GAU_C` | 无默认，需给出 | 反电荷高斯的位置（距晶胞边界）与宽度，单位 Å |

输出中会报告每步的功函数、与目标值之差、电荷和 GCE 自由能。按照输入说明：GCE 电荷为正表示缺电子；`&DFT/CHARGE` 定义参考电荷态；GCE 自由能等于 DFT 总能量加上目标功函数与电子数项的乘积。这是 2026.2 刚加入的功能，使用前请先在小体系上复现官方 regtest，并关注后续版本的文档更新。

---

## 参数选择与收敛建议

| 项目 | 推荐起点 | 判据 / 说明 |
|---|---|---|
| cube `STRIDE` | 画图 2；定量分析 1 | Bader、积分、平面平均、差分一律用 1 |
| `&MO_MOLDEN/NDIGITS` | 9 | 默认 3 只够画图 |
| `CUTOFF`（实空间分析） | 比能量收敛值高约 50–100 Ry | Bader 电荷变化 < 0.01 e |
| slab 真空层 | ≥ 15 Å，功函数建议 ≥ 20 Å | 平面平均势在真空中出现平台 |
| slab 层数 | 金属 ≥ 4 层，氧化物按化学计量单元 | 功函数变化 < 0.05 eV |
| 面内 k 点 | 金属 p(2×2) 胞 ≥ 6×6 | 功函数、费米能稳定 |
| 偶极校正 | 不对称 slab 必开 | 两侧平台都平直 |
| 差分密度 | 三次计算完全相同的晶胞、`CUTOFF`、`STRIDE`、k 点 | $\int\Delta\rho\approx0$ |
| 电荷方法 | 比较趋势时固定一种 | 论文写明方法与基组 |

---

## 常见问题与报错

**1. 运行结束但没有任何 cube 文件。**
原因：cube 打印键默认级别为 HIGH，只写了 `PRINT_LEVEL MEDIUM` 而没写对应的 `&..._CUBE` 段；或者 GEO_OPT/MD 中 `&EACH` 设成了 0。解决：显式写出 `&E_DENSITY_CUBE ... &END`，需要时加 `ADD_LAST NUMERIC`。

**2. k 点计算中没有轨道 cube，输出有关于 MO cube 与 k 点的警告。**
原因：`MO_CUBES` 不支持 k 点（molden 同样不支持）。解决：用 Γ 点超胞计算轨道；密度和 Hartree 势 cube 在 k 点下可以正常输出。

**3. `Dipole correction needs more vacuum space above the surface`。**
原因：真空层太薄或结构跨越边界，偶极校正找不到电子密度足够小的层。解决：加厚真空（≥ 15 Å），检查 slab 是否被周期边界切开，必要时设 `SURF_DIP_POS`。

**4. 功函数为负或大得离谱。**
原因：误把 `V_HARTREE` 取负号（CP2K 已是电子势能约定）；单位混用（cube 与 `Fermi energy:` 都是 Hartree）；平台取在了 slab 内部或偶极层附近。解决：画出 $V(z)-E_F$ 曲线检查平台位置。

**5. 不对称 slab 两侧真空能级相差很大，功函数“取哪一侧”？**
两侧各对应一个表面的功函数，二者之差等于表面偶极引起的势阶。报告时注明是哪一侧（例如“吸附物一侧”）。若没开偶极校正，两侧都不可信。

**6. 差分电荷密度出现大面积、与成键无关的正负区域，或 $\int\Delta\rho$ 明显不为 0。**
原因：A、B 坐标与 AB 不一致，或三个计算参数不同、电子数不守恒（例如 CO 单独计算时设了 `CHARGE`）。解决：从 AB 结构直接删原子生成 A、B；核对三份输出中的电子数（`grep "Number of electrons"`）。

**7. Bader 电荷之和与价电子总数不符，或出现大量非核吸引子。**
原因：网格太粗（`STRIDE` > 1）或 `CUTOFF` 偏低。解决：`STRIDE 1 1 1`，提高 `CUTOFF` 复查；用 Multiwfn 的 AIM 分析交叉检验。

**8. Multiwfn 载入 CP2K molden 后，原子电荷明显不对，或周期体系的分析结果在盒子边界附近异常。**
原因：molden 中缺少有效核电荷或晶胞信息。解决：2026.2 加 `WRITE_PSEUDO T` 和 `WRITE_CELL T`；旧版本按 sobereva.com/651 手工插入 `[Nval]`/`[Cell]`。

**9. molden 中没有空轨道，或 OT 计算的轨道能量全为 0。**
原因：旧版本 OT 不写空轨道和轨道能量。解决：2026.2 用 `&MO_MOLDEN/NLUMO`；旧版本改用对角化单点（从 OT 的 `RESTART.wfn` 起步）并设 `ADDED_MOS`。

**10. cube 文件太大，磁盘被写满。**
原因：大盒子 + `STRIDE 1` + MD 中频繁输出。解决：增大 `STRIDE`，用 `&EACH` 降低频率；`MO_CUBES`（≥ 2026.1）、`V_HARTREE_CUBE`（2026.2）可用 `MAX_FILE_SIZE_MB` 限制大小。

---

## 练习

1. 对 H2O 分别用 `STRIDE 1`、`2`、`3` 输出密度，比较对 cube 积分得到的电子数和 Bader 电荷，体会网格对定量分析的影响。
2. 把 H2O 换成 O2 三重态（`UKS`、`MULTIPLICITY 3`，见[第8章](08-dispersion-spin-dftu.md)），输出 `SPIN_DENSITY` cube，并对照 Mulliken 输出中的 `Spin moment` 列。
3. 对 Pt(111) 示例分别取 3、4、5、6 层和 15、20、25 Å 真空，画出功函数随层数和真空厚度的变化。
4. 在 Pt(111) 的一侧放一个 CO（顶位），分别在开、关偶极校正时计算两侧的真空能级，解释二者的差别。
5. 对 CO/Cu(111) 示例比较 Hirshfeld、Bader（`bader` 程序）和 $\Delta Q(z)$ 给出的电荷转移量。

---

## 参考资料

**官方文档**

- 输入参考 `&DFT/&PRINT`（各 cube、电荷分析、molden 打印键，2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/PRINT.html>
- 输入参考 `&V_HARTREE_CUBE`：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE.html>
- 输入参考 `&SCF/&GCE`（巨正则 SCF，2026.2+）：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/SCF/GCE.html>
- K-Points（molden 不支持 k 点等说明）：<https://manual.cp2k.org/trunk/methods/dft/k-points.html>
- CP2K v2026.2 发布说明（恒电势 SCF #5467、平面反电荷 #5064、新 Pulay 混合 #5407）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- CP2K v2026.1 发布说明（`MO_CUBES` 的 `MAX_FILE_SIZE_MB`）：<https://github.com/cp2k/cp2k/releases/tag/v2026.1>
- Bader 程序（Henkelman 组）：<https://theory.cm.utexas.edu/henkelman/code/bader/>
- sobereva：详谈使用 CP2K 产生给 Multiwfn 用的 molden 格式的波函数文件：<http://sobereva.com/651>
- sobereva：使用 Multiwfn 非常方便地创建 CP2K 程序的输入文件：<http://sobereva.com/587>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-gcscf/`（恒电势 SCF）、`tests/QS/regtest-slab/`（偶极校正）等 regtest 输入
- [funkymunkycool/Cube-Toolz](https://github.com/funkymunkycool/Cube-Toolz)：cube 差分、平面平均、积分、平移，适用于 CP2K 的 `ELECTRON_DENSITY` 与 `v_hartree` cube
- [RobinLi-95/CP2K_tools](https://github.com/RobinLi-95/CP2K_tools)：V_HARTREE 平面平均（vtotav）、功函数计算、ghost layer slab 等脚本
- [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata)：读 cube、pdos、输出文件的后处理库
- [CommonCE/cp2kkit](https://github.com/CommonCE/cp2kkit)：类似 VASPKIT 的工具箱，含 cube 差分与平面平均（预发布版，无 license）
- [nanotech-empa/cp2k-spm-tools](https://github.com/nanotech-empa/cp2k-spm-tools)：读 `.wfn` 做 STM/STS 模拟和轨道杂化分析
- [ljpyztd-cyber/CP2K_Forge](https://github.com/ljpyztd-cyber/CP2K_Forge)：调用 Multiwfn 生成 CP2K 输入的 Windows GUI（2026.1/2026.2）
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含 slab 与电化学界面相关的 CP2K 专题

**关键文献**

- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999)，表面偶极校正：<https://doi.org/10.1103/PhysRevB.59.12301>
- F. L. Hirshfeld, *Theor. Chim. Acta* **44**, 129 (1977)，Hirshfeld 划分：<https://doi.org/10.1007/BF00549096>
- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* **36**, 354 (2006)，快速 Bader 划分算法：<https://doi.org/10.1016/j.commatsci.2005.04.010>
- A. D. Becke, K. E. Edgecombe, *J. Chem. Phys.* **92**, 5397 (1990)，ELF：<https://doi.org/10.1063/1.458517>
- T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012)，Multiwfn：<https://doi.org/10.1002/jcc.22885>
- T. Lu, *J. Chem. Phys.* **161**, 082503 (2024)，Multiwfn 综述：<https://doi.org/10.1063/5.0216272>
- G. N. Derry, M. E. Kern, E. H. Worth, *J. Vac. Sci. Technol. A* **33**, 060801 (2015)，金属表面功函数推荐值：<https://doi.org/10.1116/1.4934685>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](10-dos-pdos-band.md) | [返回目录](../README.md) | [下一章 →](12-aimd.md)
