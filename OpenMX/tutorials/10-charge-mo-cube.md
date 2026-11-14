# 第10章 电荷分析、分子轨道与 cube 可视化（Charge Analysis, Molecular Orbitals & Cube Visualization with OpenMX）

本章讲解如何用 OpenMX 输出并分析实空间数据：总电荷密度、自旋密度、形变密度与 Kohn–Sham 势的 Gaussian cube 文件（tden/sden/dden/v0/vhart），Mulliken、Voronoi、ESP 拟合与自然布居分析（NPA）四类原子电荷，分子轨道与晶体轨道（HOMO/LUMO cube），吸附体系的电荷密度差（diff_gcube），以及 VESTA / XCrySDen / OpenMX Viewer 可视化和 `HS.fileout` 生成的 scfout 文件（OpenMX charge analysis, Mulliken / Voronoi / ESP charges, molecular orbitals, charge density difference, cube visualization tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 用 `level.of.fileout` 控制 OpenMX 写出哪些 cube 文件，读懂 cube 文件头，用 Python 对 cube 积分核对电子数；
- 从 `System.Name.out` 中读出 Mulliken 电荷及其轨道分解，打开 `Voronoi.charge`，用后处理程序 `esp` 计算 ESP 拟合电荷，用 `NBO.switch` 做 NPA，并理解各方法为什么给出不同的数；
- 用 `MO.fileout`、`num.HOMOs`、`num.LUMOs` 输出分子的 HOMO/LUMO，用 `MO.Nkpoint`、`<MO.kpoint` 输出周期体系指定 k 点的晶体轨道（实部 `_r` / 虚部 `_i`）；
- 按 AB − A − B 的流程（`scf.fixed.grid` + `add_gcube`/`diff_gcube`）计算吸附体系的电荷密度差，并做平面平均与电荷转移积分；
- 用 `OutData.bin.flag` + `bin2txt` 处理大文件，用 `cube2xsf` 转换周期 XSF 格式；
- 理解 `HS.fileout on` 生成的 `System.Name.scfout` 是 jx、polB、kSpin、calB、Z2FH、TB2J、DeepH 等后处理的共同入口。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| 看成键、孤对电子、形变密度 | `level.of.fileout 1`，看 `dden.cube` | §3、§4 |
| 比较原子电荷、电荷转移趋势 | Mulliken + Voronoi + ESP/NPA 对照，看**差值**而不是绝对值 | §4 |
| 分子前线轨道、反应活性位点 | `MO.fileout on`（Cluster） | §5 |
| 晶体能带边的轨道特征 | `MO.fileout on` + `<MO.kpoint`（Band） | §6 |
| 吸附/界面电荷重新分布 | AB − A − B 电荷密度差 + 平面平均 | §7 |
| 上千原子、cube 写盘太慢 | `OutData.bin.flag on` 或 `level.of.fileout 0` | §8 |
| 交换耦合、极化、拓扑、机器学习哈密顿量 | `HS.fileout on` → scfout | §10 |

> [!NOTE]
> 本章假设你已经能独立完成一次 SCF 计算（[第3章](03-input-first-calc.md)），并理解基组与收敛参数（[第4章](04-pao-vps.md)、[第5章](05-scf-kpoints-convergence.md)）。态密度/PDOS 见 [第9章](09-dos-pdos-cohp-optics.md)，VASP 下的同类分析（Bader、CDD、ELF）可对照 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)。

---

## 2. 原理速览

### 2.1 电子密度、形变密度与势

在 LCPAO（线性组合赝原子轨道）框架下，Kohn–Sham 轨道写成 $\psi_{n\mathbf k}(\mathbf r)=\sum_{i\alpha}c_{n\mathbf k,i\alpha}\,\phi_{i\alpha}(\mathbf r-\boldsymbol\tau_i)$（周期体系再乘 Bloch 相因子并对晶格平移求和）。价电子密度为

$$
\rho(\mathbf r)=\sum_{\sigma}\sum_{n\mathbf k} w_{\mathbf k}\, f_{n\mathbf k\sigma}\,|\psi_{n\mathbf k\sigma}(\mathbf r)|^2
=\sum_{\sigma}\sum_{i\alpha,j\beta} \rho^{\sigma}_{i\alpha,j\beta}\,\phi_{i\alpha}(\mathbf r)\phi_{j\beta}(\mathbf r),
$$

其中 $\rho^{\sigma}_{i\alpha,j\beta}$ 是密度矩阵。OpenMX 写出的几个实空间量：

- **tden**：总价电子密度 $\rho=\rho_\uparrow+\rho_\downarrow$（赝势框架，不含芯电子）；
- **sden**：自旋密度 $\rho_\uparrow-\rho_\downarrow$（共线自旋极化时）；
- **dden**：形变密度（deformation density）$\rho-\sum_i\rho_i^{\rm atom}$，即相对于"孤立中性原子密度简单叠加"的变化，最适合看成键电荷积累和孤对电子；
- **v0 / v1**：上/下自旋的 Kohn–Sham 势（不含非局域赝势部分；若 `scf.ProExpn.VNA on`，也不含 VNA 势）；
- **vhart**：由差分密度 $\delta\rho=\rho-\rho^{\rm atom}$ 产生的 Hartree 势 $\delta V_{\rm H}$；
- **vxc0 / vxc1**：上/下自旋交换关联势（`level.of.fileout 2`）。

> [!IMPORTANT]
> `dden`（相对孤立原子）与本章 §7 的"电荷密度差"（相对两个孤立片段 A、B）是两回事：前者反映"原子 → 分子/晶体"的全部成键变化，后者只反映 A 与 B 之间相互作用引起的变化。论文里画吸附电荷转移图时用的是后者。

### 2.2 原子电荷没有唯一定义

把连续的 $\rho(\mathbf r)$ 分给各个原子本质上带有人为约定，不同方法的数值可以差很多，**只能在同一方法内部比较趋势**。OpenMX 4.0 提供以下几种（手册第 35、36 章）：

| 方法 | 原理 | 优点 | 主要缺点 | OpenMX 中如何得到 |
|---|---|---|---|---|
| Mulliken | 把重叠布居 $D_{\mu\nu}S_{\nu\mu}$ 平分给两个原子：$N_A=\sum_{\mu\in A}\sum_\nu D_{\mu\nu}S_{\nu\mu}$ | 默认输出，可分解到每个轨道、每个自旋 | 强烈依赖基组，大基组/弥散函数时可能失真 | 默认写入 `.out` |
| Voronoi | 在模糊 Voronoi 多面体（fuzzy cell）内积分密度 | 只依赖结构，致密固体中较稳定，还给出 Voronoi 体积 | 不考虑原子半径，分子中常与化学直觉相反 | `Voronoi.charge on` |
| ESP 拟合 | 拟合点电荷使其再现分子外部的静电势 | 小分子中与化学直觉吻合，适合给力场用 | 体相/大分子内部原子缺少采样点，不适用 | `esp` 后处理程序 |
| NPA | 由自然原子轨道（NAO）的占据数给出自然布居 | 对基组较不敏感，化学解释清晰 | 需额外计算；大体系要用 O(N) 模式 | `NBO.switch on1/on2` |
| CWF | 用最接近 Wannier 函数（4.0 新增）的占据数 | 与局域轨道图像对应 | 需先构造 CWF | 见 [第14章](14-wannier.md) |
| Bader | 按密度零通量面划分（外部程序） | 常用于固体电荷转移 | OpenMX 只给价电子密度，需谨慎 | 外部 `bader` 读 `tden.cube` |

OpenMX 4.0 手册（第 35 章）用甲烷给出了一个说明问题的对照（数值摘自手册）：Mulliken 给 C 约 5.02 个价电子（带负电），ESP 给 C 约 −1.07 的有效电荷（也带负电），而 Voronoi 给 C 只有约 2.18 个电子、H 约 1.46 个电子——C 反而带正电，与化学直觉相反，原因是 Voronoi 多面体完全不考虑 C、H 原子尺寸差异。这正说明：**电荷的绝对值没有意义，变化趋势才有意义**。

### 2.3 分子轨道与晶体轨道

分子（Cluster）计算中轨道是实函数，HOMO/LUMO 可以直接画等值面（正负相位用两种颜色）。周期（Band）计算中 Bloch 态一般是复函数，OpenMX 把实部和虚部分别写成 `_r.cube` 和 `_i.cube`；它们各自依赖于任意的整体相位，物理上有意义的是 $|\psi|^2=\psi_r^2+\psi_i^2$。只有在 Γ 等时间反演不变点（且无 SOC）时，才能通过选择相位得到实波函数。

### 2.4 电荷密度差

吸附/界面体系常用

$$
\Delta\rho(\mathbf r)=\rho_{AB}(\mathbf r)-\rho_{A}(\mathbf r)-\rho_{B}(\mathbf r),
$$

其中 A、B 两个片段必须保持在复合体系 AB 中的几何结构，并使用**同一个晶胞、同一套实空间网格（包括网格原点）和相同的计算参数**。沿表面法向做平面平均 $\Delta\lambda(z)=\iint\Delta\rho\,dx\,dy$，再积分 $\Delta Q(z)=\int_{-\infty}^{z}\Delta\lambda(z')dz'$，就得到"跨过某个平面的电子转移量"。

---

## 3. cube 文件：输出控制与格式

### 3.1 level.of.fileout 决定写出哪些文件

`level.of.fileout` 默认是 1。各级别输出（据 OpenMX 4.0 手册第 9 章）：

| 文件 | 条件 | 内容 | 单位 |
|---|---|---|---|
| `System.Name.out`、`.xyz`、`.md`、`.md2`、`.cif`、`.ene`、`_rst/` | 任何级别（含 0） | 能量、Mulliken、偶极、力、结构、重启文件 | — |
| `System.Name.tden.cube` | `level.of.fileout` ≥ 1 | 总价电子密度 | e/Bohr³ |
| `System.Name.sden.cube` | ≥ 1 且自旋极化（共线） | 自旋密度 ρ↑−ρ↓ | e/Bohr³ |
| `System.Name.dden.cube` | ≥ 1 | 形变密度 ρ − Σρ_atom | e/Bohr³ |
| `System.Name.v0.cube` / `.v1.cube` | ≥ 1（v1 仅自旋极化） | KS 势（不含非局域部分） | Hartree |
| `System.Name.vhart.cube` | ≥ 1 | 差分密度的 Hartree 势 | Hartree |
| `System.Name.vxc0.cube` / `.vxc1.cube` | = 2 | 交换关联势 | Hartree |
| `System.Name.grid` | = 2 | 实空间积分网格 | — |
| `System.Name.nc.xsf`、`.ncsden.xsf`、`.nco.xsf` | ≥ 1 且 `scf.SpinPolarization NC` | 非共线自旋/轨道磁矩矢量（XCrySDen 用） | — |
| `System.Name.homo*/lumo*.cube` | `MO.fileout on` | 分子/晶体轨道 | Bohr^(-3/2) |
| `System.Name.scfout` | `HS.fileout on` | 哈密顿量、重叠、密度矩阵（二进制） | — |

> [!TIP]
> 结构优化、MD、高通量筛选时把 `level.of.fileout` 设为 0，可省下大量磁盘与 I/O 时间；最后在优化好的结构上用 `scf.restart on` 或 `pp`（见 [第5章](05-scf-kpoints-convergence.md)）补一次 `level.of.fileout 1` 的计算即可。

### 3.2 cube 格式速读

OpenMX 写出的 cube 是标准 Gaussian cube 文本格式：

```text
 SYS1
 SYS1
    3   -7.937866   -7.937866   -7.937866
   90    0.176397    0.000000    0.000000
   90    0.000000    0.176397    0.000000
   90    0.000000    0.000000    0.176397
    8    0.6xxxxx    0.000000    0.000000    0.000000
    1    0.3xxxxx    1.430178    1.107485    0.000000
    1    0.3xxxxx   -1.430178    1.107485    0.000000
  1.234E-03  1.456E-03 ...
```

（示意输出，数值仅供参考）

| 行 | 含义 |
|---|---|
| 1–2 | 注释行。密度/势文件为 `SYS1`；轨道文件会写出该本征值（绝对值与相对化学势的值，Hartree） |
| 3 | 原子数、网格原点 (x, y, z)，单位 Bohr。原点就是 `.std` 中打印的 `Grid_Origin` |
| 4–6 | 沿 a、b、c 方向的网格点数 N1、N2、N3 及每一格的位移矢量（Bohr） |
| 原子行 | 原子序数、一个实数（OpenMX 源码中写的是该原子的价电荷减去其 Mulliken 布居，即净电荷）、坐标（Bohr） |
| 数据 | 按 n1（最慢）、n2、n3（最快）顺序排列，每行 6 个数 |

几个实用检查：

```bash
head -n 6 h2o.tden.cube          # grid numbers and voxel vectors
grep -c "" h2o.tden.cube         # rough size check (lines)
ls -lh *.cube                    # cube files are large: N1*N2*N3 numbers each
```

### 3.3 用 Python 积分 cube，核对电子数

对 `tden.cube` 积分应得到体系的价电子数（例如 H2O：O_PBE19 有 6 个价电子，两个 H 各 1 个，共 8 个）。下面这段独立小脚本只依赖 numpy：

```python
import sys
import numpy as np

def read_cube(fname):
    with open(fname) as f:
        f.readline(); f.readline()
        tok = f.readline().split()
        natom = abs(int(tok[0])); origin = np.array(tok[1:4], float)
        n, vec = [], []
        for _ in range(3):
            tok = f.readline().split()
            n.append(int(tok[0])); vec.append([float(x) for x in tok[1:4]])
        atoms = [f.readline().split() for _ in range(natom)]
        data = np.array(f.read().split(), float).reshape(n)
    return origin, np.array(vec), atoms, data

origin, vec, atoms, rho = read_cube(sys.argv[1])
dV = abs(np.linalg.det(vec))            # voxel volume in Bohr^3
print("grid =", rho.shape, " voxel volume (Bohr^3) =", dV)
print("integral =", rho.sum() * dV)     # number of electrons for tden.cube
```

```bash
python cube_integrate.py h2o.tden.cube
```

积分值与理论电子数相差 1e-2 量级以内是正常的（网格离散化误差）；若相差很大，说明网格太粗（提高 `scf.energycutoff`）或分子跨过了晶胞边界又被截断（Cluster 计算请把盒子取大）。对 `dden.cube` 积分应接近 0；对 `sden.cube` 积分应接近总自旋磁矩（μB）。更完整的 cube 工具（平面平均、密度差、积分）见仓库脚本 [`omx_cube_tools.py`](../scripts/omx_cube_tools.py)，用 `python ../scripts/omx_cube_tools.py -h` 查看用法。

---

## 4. 实战一：水分子的四种原子电荷（H2O_cube.dat）

完整文件见 [../examples/10-h2o-cube-charge/README.md](../examples/10-h2o-cube-charge/README.md)。

### 4.1 输入文件

```text
#
# H2O molecule: cube output and four kinds of atomic charges
# OpenMX 4.0 (patch 4.0.1); DFT_DATA19
#
System.CurrrentDirectory    ./
System.Name                 h2o
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            2          # 2: also vxc cubes, .grid and LCAO coefficients

Species.Number              2
<Definition.of.Atomic.Species
  O   O6.0-s2p2d1   O_PBE19
  H   H6.0-s2p1     H_PBE19
Definition.of.Atomic.Species>

Atoms.Number                3
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1  O    0.000000   0.000000   0.000000   3.0  3.0
  2  H    0.756810   0.586060   0.000000   0.5  0.5
  3  H   -0.756810   0.586060   0.000000   0.5  0.5
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  12.0   0.0   0.0
   0.0  12.0   0.0
   0.0   0.0  12.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.EigenvalueSolver        Cluster
scf.energycutoff            250.0
scf.maxIter                 100
scf.criterion               1.0e-8
scf.ElectronicTemperature   300.0
scf.Mixing.Type             RMM-DIIS
scf.Init.Mixing.Weight      0.30
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.40
scf.Mixing.History          7
scf.Mixing.StartPulay       5

MD.Type                     Nomd

Voronoi.charge              on         # Voronoi charges -> h2o.out
NBO.switch                  on1        # natural population analysis -> standard output
HS.fileout                  off
```

### 4.2 关键词说明

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `level.of.fileout` | 2 | 1 已能输出 tden/dden/v0/vhart；2 额外输出 vxc0/vxc1、`.grid`，并在 Cluster 计算中把 LCAO 系数写进 `.out` |
| `Atoms.UnitVectors` | 12 Å 立方盒 | Cluster 计算若不给晶胞，OpenMX 会自动定一个"刚好让基函数不重叠"的盒子；画 cube 等值面时盒子大一点更好看，ESP 拟合也需要分子外有足够空间 |
| `scf.energycutoff` | 250 Ry | 实空间网格越密，cube 越平滑、积分越准；Hartree 势与 ESP 对网格较敏感，建议不低于 200 Ry |
| `scf.EigenvalueSolver` | Cluster | 孤立分子必须用 Cluster；`MO.fileout` 与 NBO `on1` 都要求 Cluster |
| `Voronoi.charge` | on | 默认 off；结果追加到 `.out` 的 "Voronoi charges" 段 |
| `NBO.switch` | on1 | Cluster 计算用全密度矩阵做 NPA，结果打印到标准输出（`.std`）；Krylov O(N) 大体系用 `on2` |
| 坐标行第 6、7 列 | O 3.0 3.0，H 0.5 0.5 | 上/下自旋初始电子数，和必须等于 VPS 价电子数（O_PBE19 = 6，H_PBE19 = 1） |

### 4.3 运行与输出

```bash
mpirun -np 4 openmx H2O_cube.dat -nt 1 > h2o.std 2>&1
ls h2o.*
# h2o.out h2o.std h2o.xyz h2o.tden.cube h2o.dden.cube h2o.v0.cube h2o.vhart.cube
# h2o.vxc0.cube h2o.grid ... h2o_rst/
grep "Utot\." h2o.out
grep -A 12 "Dipole moment (Debye)" h2o.out
```

SCF 是否收敛：看 `h2o.out` 中 "SCF history" 最后一步的 `NormRD` 是否足够小（远小于 1e-6），以及 `.std` 中迭代次数是否小于 `scf.maxIter`。偶极矩：实验气相值约 1.85 D，GGA-PBE + Standard 基组一般给出相近数值；偏差很大时先检查几何与盒子大小。

### 4.4 Mulliken 电荷与轨道分解

```bash
grep -A 8 "Mulliken populations" h2o.out | head -12
grep -A 40 "Decomposed Mulliken populations" h2o.out | less
```

"Mulliken populations" 表中每个原子给出 Up spin、Down spin、Sum、Diff 四列。**Sum 是电子布居，不是电荷**：原子净电荷 = VPS 价电子数 − Sum。例如若 O 的 Sum 为 6.6，则 O 的 Mulliken 电荷约 −0.6 e。表末的

```text
 Sum of MulP: up   =     4.00000 down          =     4.00000
              total=     8.00000 ideal(neutral)=     8.00000
```

（示意输出）用来核对总电子数。"Decomposed Mulliken populations" 把布居分解到每个基函数：

| 行 | 含义 |
|---|---|
| `s 0`、`s 1` | 第 1、第 2 组 s 径向函数（对应基组记号 s2 中的两组） |
| `px py pz` | p 轨道的三个实球谐分量 |
| `d3z^2-r^2 dx^2-y^2 dxy dxz dyz` | d 轨道，**OpenMX 的顺序**（与 VASP 的 dxy、dyz、dz2、dxz、dx2-y2 不同） |
| `f5z^2-3r^2 ... f3yx^2-y^3` | f 轨道 7 个分量 |
| `sum over m` | 对磁量子数求和（例如某组 p 的 px+py+pz） |
| `sum over m+mul` | 再对多组径向函数（multiplicity）求和，即该角动量的总布居 |

> [!WARNING]
> Mulliken 布居强烈依赖基组：把 O6.0-s2p2d1 换成 O6.0-s3p2d2，O 的电荷可能明显变化。比较不同体系时必须用相同的 PAO 与 VPS；论文中报告 Mulliken 电荷时应注明基组。

### 4.5 Voronoi 电荷

```bash
grep -A 14 "Voronoi charges" h2o.out
```

输出包括各自旋的 Voronoi 电荷总和、由 Voronoi 电荷给出的总自旋磁矩（"Total spin magnetic moment (muB) by Voronoi charges"）、每个原子的 Up/Down/Sum/Diff，以及 **Voronoi 体积（Å³）**。OpenMX 4.0 用模糊边界（fuzzy-cell partitioning）构造 Voronoi 多面体。在分子中 H 原子的 Voronoi 多面体会延伸到真空中，"抢到"较多密度，所以 Voronoi 电荷常常与化学直觉相反；在致密固体（合金、金属间化合物）中，Voronoi 电荷与体积则是有用的局域结构指标。自旋极化磁性体系的局域磁矩，手册也推荐用 Voronoi 划分（多组基函数时 Mulliken 自旋矩偏大）。

### 4.6 ESP 拟合电荷（esp 程序）

ESP 拟合需要 `level.of.fileout` ≥ 1 时生成的两个文件：`System.Name.out` 和 `System.Name.vhart.cube`。先编译后处理程序（4.0 中 `esp` 是 `source/` 下 makefile 的一个目标，编译后自动复制到 `work/`）：

```bash
cd ~/openmx4.0/source
make esp                     # executable is copied to ../work/esp
```

然后在算例目录运行（参数是 System.Name，不带扩展名）：

```bash
~/openmx4.0/work/esp h2o -c 0 -s 1.4 2.0 | tee h2o.esp
grep "Fitting Effective Charge" h2o.esp
grep "RMS" h2o.esp
```

| 选项 | 含义 | 建议 |
|---|---|---|
| `-c 0` | 只约束总电荷守恒 | **推荐** |
| `-c 1` | 同时约束总电荷与偶极矩 | 大分子中可能出现不合理的点电荷 |
| `-s 1.4 2.0` | 采样壳层：各原子 vdW 半径乘以 1.4 与 2.0 所围成的两层壳之间的网格点 | 与 Merz–Kollman 方案的常用设置一致 |

判断：ESP 电荷之和应等于体系总电荷（中性分子为 0）；输出的偶极矩应与 `.out` 中 DFT 偶极接近；最后一行给出拟合 RMS。ESP 只适用于小/中等分子，不适用于体相（内部原子在分子外表面附近没有采样点）。

### 4.7 自然布居分析（NPA）

本例已设 `NBO.switch on1`，结果在标准输出中：

```bash
grep -B 2 -A 6 "Total  :" h2o.std     # NP of every atom and their sum
```

首先打印每个原子的自然布居（NP，电子数）及总和，随后对每个原子列出 NAO 的占据数、能量（Hartree）与对应的 LCPAO 系数。相关参数：

| 关键词 | 默认 | 说明 |
|---|---|---|
| `NBO.switch` | off | `on1`：Cluster 计算、全密度矩阵；`on2`：O(N) Krylov 计算（`scf.EigenvalueSolver Krylov`），只分析指定原子 |
| `NBO.Num.CenterAtoms` + `<NBO.CenterAtoms` | — | 仅 `on2` 有效：要分析的原子数与原子序号列表 |
| `NAO.threshold` | 0.85 | 区分"占据 NAO"与 Rydberg NAO 的阈值（每个自旋轨道的电子数）；得到不合理的占据 NAO 数时再调 |

官方 work 目录中的 `EC_NAO.dat`（碳酸乙烯酯分子）和 `SiO2_NAO.dat`（648 原子 SiO2，`on2`）是两个可直接对照的例子。

### 4.8 结果怎么读：方法对照

运行完后把四种电荷整理成一张表（建议用 [`omx_parse_out.py`](../scripts/omx_parse_out.py) 自动解析 `.out`，用 `-h` 查看用法）：

| 原子 | Mulliken | Voronoi | ESP（-c 0） | NPA |
|---|---|---|---|---|
| O | 负 | 可能接近 0 或为负 | 负 | 负 |
| H | 正 | 可能为负（多面体伸入真空） | 正 | 正 |

（上表只给出符号趋势，具体数值以你的计算为准。）

合理的结论应该是"O 从 H 处得到电子"——Mulliken、ESP、NPA 一致，而 Voronoi 在分子中不可靠。若要研究化学反应前后某原子的电荷**变化**，选定一种方法后在同一基组、同一参数下比较差值。

### 4.9 Bader 与 CWF

- **Bader**：Henkelman 组的 `bader` 程序可直接读 Gaussian cube：`bader h2o.tden.cube`，输出 `ACF.dat`。注意 OpenMX 的 `tden` 只含价电子（赝密度），核附近没有芯电子峰，Bader 分区在重原子周围可能出现非核吸引子；结果中的电子数应与 VPS 价电子数对照理解，必要时用更密网格（更高 `scf.energycutoff`）。
- **CWF**：4.0 新增的最接近 Wannier 函数可以给出"原子有效电荷与自旋"，用法见 [第14章](14-wannier.md)。

---

## 5. 实战二：苯分子的 HOMO / LUMO（Benzene_MO.dat）

完整文件见 [../examples/10-benzene-mo/README.md](../examples/10-benzene-mo/README.md)。

### 5.1 输入文件

```text
#
# Benzene molecule: HOMO/LUMO cube files (cluster calculation)
#
System.CurrrentDirectory    ./
System.Name                 benzene
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            1

Species.Number              2
<Definition.of.Atomic.Species
  C   C6.0-s2p2d1   C_PBE19
  H   H6.0-s2p1     H_PBE19
Definition.of.Atomic.Species>

Atoms.Number                12
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1  C    1.39700   0.00000   0.00000   2.0  2.0
  2  C    0.69850   1.20984   0.00000   2.0  2.0
  3  C   -0.69850   1.20984   0.00000   2.0  2.0
  4  C   -1.39700   0.00000   0.00000   2.0  2.0
  5  C   -0.69850  -1.20984   0.00000   2.0  2.0
  6  C    0.69850  -1.20984   0.00000   2.0  2.0
  7  H    2.48100   0.00000   0.00000   0.5  0.5
  8  H    1.24050   2.14861   0.00000   0.5  0.5
  9  H   -1.24050   2.14861   0.00000   0.5  0.5
 10  H   -2.48100   0.00000   0.00000   0.5  0.5
 11  H   -1.24050  -2.14861   0.00000   0.5  0.5
 12  H    1.24050  -2.14861   0.00000   0.5  0.5
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  16.0   0.0   0.0
   0.0  16.0   0.0
   0.0   0.0  12.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.EigenvalueSolver        Cluster
scf.energycutoff            220.0
scf.maxIter                 100
scf.criterion               1.0e-8
scf.Mixing.Type             RMM-DIIS
scf.Init.Mixing.Weight      0.30
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.40
scf.Mixing.History          7
scf.Mixing.StartPulay       5

MD.Type                     Nomd

MO.fileout                  on         # write HOMO/LUMO cube files
num.HOMOs                   2          # HOMO and HOMO-1
num.LUMOs                   2          # LUMO and LUMO+1
MO.Nkpoint                  1          # must be 1 for Cluster
<MO.kpoint
  0.0  0.0  0.0
MO.kpoint>
```

几何取近似实验键长（C–C 1.397 Å、C–H 1.084 Å），未做优化；要发表结果请先按 [第6章](06-optimization-eos-neb.md) 优化。

### 5.2 关键词说明

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `MO.fileout` | on | 默认 off。只在 `scf.EigenvalueSolver` 为 Cluster、Band（以及少数对角化类求解器）时有效；用 DC、Krylov 等 O(N) 求解器时 OpenMX 会打印 "MO.fileout=ON is not supported ..." 并自动改为 OFF |
| `num.HOMOs` / `num.LUMOs` | 2 / 2 | 输出的占据/空轨道个数。苯的 HOMO（e1g）与 LUMO（e2u）都是二重简并，所以至少各取 2 才能看到完整的一对 |
| `MO.Nkpoint` | 1 | Cluster 计算必须为 1，否则报 "MO.Nkpoint should be 1 for the cluster calculation." |
| `<MO.kpoint` | 0 0 0 | Cluster 时写 Γ 即可 |

### 5.3 输出文件命名

| 文件 | 含义 |
|---|---|
| `benzene.homo0_0.cube` | 自旋 0（上自旋；非自旋极化时只有 0）的 HOMO |
| `benzene.homo0_1.cube` | HOMO−1 |
| `benzene.lumo0_0.cube` | LUMO |
| `benzene.lumo0_1.cube` | LUMO+1 |

命名规则：`homo{自旋}_{序号}`，自旋 0 = up、1 = down；序号 0、1、2 分别对应 HOMO、HOMO−1、HOMO−2（LUMO 类似）。每个轨道 cube 的第 1 行注释写着该本征值（绝对值与相对化学势，Hartree）：

```bash
head -n 2 benzene.homo0_0.cube
grep -A 25 "Eigenvalues (Hartree)" benzene.out | head -30
grep "HOMO =" benzene.out
```

苯共有 30 个价电子（6×4 + 6×1），占据 15 个空间轨道，所以 `.out` 中应看到 `HOMO = 15`，第 14、15 个本征值近似简并，第 16、17 个也近似简并。HOMO–LUMO 能隙 = (ε16 − ε15) × 27.2114 eV；GGA 的 KS 能隙系统性低于实验光学/电学能隙，不要直接与实验比较绝对值。

### 5.4 LCAO 系数

`level.of.fileout 2` 时，Cluster 计算会把全部 LCAO 系数写进 `.out`（"LCAO coefficients for up (U) and down (D) spins"）。对苯来说，HOMO/LUMO 应基本只由 C 的 pz（`0 pz`、`1 pz`）组成——这是 π 轨道最直接的数值证据。周期体系要输出 LCAO 系数，则需打开 `MO.fileout` 并在 `<MO.kpoint` 指定 k 点（不是靠 `level.of.fileout`）。

### 5.5 在 VESTA 中画轨道

1. File → Open → `benzene.homo0_0.cube`；
2. Properties → Isosurfaces：Isosurface level 取 0.02–0.05（Bohr^-3/2），勾选显示正、负两种颜色（Positive / Negative）；
3. 再 Edit → Edit Data → Volumetric data → Import 叠加 `benzene.homo0_1.cube`，对比简并的两个轨道；
4. Utilities → 2D Data Display 可以画任意截面（例如分子平面上方 0.5 Å 的平面）。

> [!TIP]
> 二重简并轨道的任意线性组合仍是本征态，所以 `homo0_0` 与 `homo0_1` 的形状可能与教科书上的标准图"转了一个角度"，这不是错误。把两者的平方和 $|\psi_{14}|^2+|\psi_{15}|^2$ 画出来就具有完整的六重对称性。

---

## 6. 周期体系的晶体轨道（MO.kpoint）

周期体系（`scf.EigenvalueSolver Band`）可以在任意 k 点输出晶体轨道。以金刚石结构 Si 的 VBM/CBM 为例，在 SCF 输入中加入：

```text
MO.fileout          on
num.HOMOs           3          # triply degenerate VBM at Gamma (without SOC)
num.LUMOs           1
MO.Nkpoint          2          # k point 0 = Gamma, k point 1 = near the CBM on Gamma-X
<MO.kpoint
  0.000  0.000  0.000
  0.425  0.425  0.000
MO.kpoint>
```

k 点坐标以 `Atoms.UnitVectors` 对应的倒格矢为基（与 `Band.kpath` 相同的约定，不受 `<Band.KPath.UnitCell` 影响）。对于以 a/2(0,1,1)、a/2(1,0,1)、a/2(1,1,0) 为基矢的 fcc 原胞，Γ–X 方向的点写成 (t, t, 0) 形式即可；确切的 CBM 位置以你的能带计算为准（见 [第8章](08-band.md)）。

输出文件：

| 文件 | 含义 |
|---|---|
| `si.homo0_0_0_r.cube` / `si.homo0_0_0_i.cube` | 第 0 个 k 点、自旋 0、HOMO 的实部/虚部 |
| `si.homo1_0_0_r.cube` | 第 1 个 k 点的 HOMO 实部 |
| `si.lumo1_0_0_r.cube` | 第 1 个 k 点的 LUMO 实部 |

命名规则 `homo{k序号}_{自旋}_{轨道序号}_{r|i}`。画概率密度时把两者平方相加：

```python
# |psi|^2 from real and imaginary parts (reuse read_cube() from Section 3.3)
o, v, at, re = read_cube("si.lumo1_0_0_r.cube")
o, v, at, im = read_cube("si.lumo1_0_0_i.cube")
dens = re**2 + im**2
print("norm in cell =", dens.sum() * abs(np.linalg.det(v)))
```

写回 cube 时保持文件头不变即可（或者直接用 [`omx_cube_tools.py`](../scripts/omx_cube_tools.py)）。同一 k 点的 LCAO 系数会写入 `.out`，可据此判断能带边的原子/轨道成分（例如 Si 的 VBM 以 p 为主）。

---

## 7. 实战三：CO/Cu(111) 吸附的电荷密度差

完整文件与一键脚本见 [../examples/10-co-cu111-cdd/README.md](../examples/10-co-cu111-cdd/README.md)。

### 7.1 流程

```text
(1) AB : CO + Cu(111) slab           -> co_cu_AB.tden.cube   + Grid_Origin in .std
(2) A  : CO only, same cell, same geometry, scf.fixed.grid = Grid_Origin of (1)
(3) B  : slab only, same cell, same geometry, scf.fixed.grid = Grid_Origin of (1)
(4) add_gcube A.cube B.cube A_B.cube ;  diff_gcube AB.cube A_B.cube dAB.cube
(5) planar average of dAB along z, cumulative charge transfer
```

为什么需要 `scf.fixed.grid`：OpenMX 会根据原子位置自动选择实空间网格的原点（`Grid_Origin`）。A、B 体系的原子不同，自动原点一般与 AB 不同，两个 cube 的网格点就不重合，逐点相减毫无意义。手册第 68 章给出的办法是：从 AB 计算的标准输出中读出 `Grid_Origin xxx yyy zzz`（单位 Bohr），在 A、B 输入中写

```text
scf.fixed.grid   xxx yyy zzz
```

晶胞与 `scf.energycutoff` 相同，则网格点数 N1×N2×N3 也相同。

### 7.2 复合体系 AB 的输入

模型：Cu 晶格常数取实验值约 3.61 Å，(111) 面 2×2 表面晶胞、3 层，CO 以 C 端竖直吸附在顶位（Cu–C 1.90 Å、C–O 1.15 Å，均为未优化的近似值，仅用于演示流程）。

```text
#
# CO on Cu(111) 2x2, 3 layers, top site: composite system AB
#
System.CurrrentDirectory    ./
System.Name                 co_cu_AB
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            1

Species.Number              3
<Definition.of.Atomic.Species
  Cu  Cu6.0S-s2p1d1   Cu_PBE19S
  C   C6.0-s2p2d1     C_PBE19
  O   O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>

Atoms.Number                14
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1  Cu   0.00000   0.00000    8.00000   5.5  5.5
  2  Cu   1.27633   2.21066    8.00000   5.5  5.5
  3  Cu   2.55266   0.00000    8.00000   5.5  5.5
  4  Cu   3.82898   2.21066    8.00000   5.5  5.5
  5  Cu   1.27633   0.73689   10.08423   5.5  5.5
  6  Cu   2.55266   2.94755   10.08423   5.5  5.5
  7  Cu   3.82898   0.73689   10.08423   5.5  5.5
  8  Cu   5.10531   2.94755   10.08423   5.5  5.5
  9  Cu   2.55266   1.47378   12.16847   5.5  5.5
 10  Cu   3.82898   3.68444   12.16847   5.5  5.5
 11  Cu   5.10531   1.47378   12.16847   5.5  5.5
 12  Cu   6.38164   3.68444   12.16847   5.5  5.5
 13  C    2.55266   1.47378   14.06847   2.0  2.0
 14  O    2.55266   1.47378   15.21847   3.0  3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  5.10531   0.00000   0.00000
  2.55266   4.42133   0.00000
  0.00000   0.00000  25.00000
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.EigenvalueSolver        Band
scf.Kgrid                   6 6 1
scf.energycutoff            220.0
scf.ElectronicTemperature   500.0
scf.maxIter                 200
scf.criterion               1.0e-7
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       10

MD.Type                     Nomd
```

A（`co_cu_A`）只保留第 13、14 号原子（编号改为 1、2），B（`co_cu_B`）只保留 1–12 号 Cu 原子；**晶胞、k 点、截断能、电子温度、混合参数全部不变**，并各加一行 `scf.fixed.grid`。

| 关键参数 | 取值 | 理由 |
|---|---|---|
| Cu 基组 / 赝势 | Cu6.0S-s2p1d1 / Cu_PBE19S（11 价电子） | DFT_DATA19 推荐表中的 Quick 组合，演示用足够；定量研究改用 Standard（Cu6.0S-s3p2d1） |
| `scf.Kgrid` | 6 6 1 | 2×2 表面晶胞（约 5.1 Å）；沿真空方向 1 个 k 点 |
| `scf.ElectronicTemperature` | 500 K | 金属板，适当展宽帮助收敛；三步必须一致 |
| `scf.Mixing.Type` | RMM-DIISK | 金属板模型最稳的混合方案 |
| 真空 | 约 17.8 Å | 由 c = 25 Å 减去板厚与 CO 高度；足够大，PAO 截断半径 6 Bohr ≈ 3.2 Å |

### 7.3 运行三个计算

```bash
mpirun -np 8 openmx CO_on_Cu111_AB.dat -nt 1 > co_cu_AB.std 2>&1
grep "Grid_Origin" co_cu_AB.std          # e.g. Grid_Origin  x  y  z  (Bohr)

ORIG=$(grep "Grid_Origin" co_cu_AB.std | tail -1 | awk '{print $2, $3, $4}')
for s in A B; do
  sed "s|^#scf.fixed.grid.*|scf.fixed.grid   ${ORIG}|" CO_on_Cu111_${s}.dat > ${s}_run.dat
  mpirun -np 8 openmx ${s}_run.dat -nt 1 > co_cu_${s}.std 2>&1
done

head -n 6 co_cu_AB.tden.cube | tail -n 4   # these 4 lines must be identical ...
head -n 6 co_cu_A.tden.cube  | tail -n 4   # ... in all three cube files
head -n 6 co_cu_B.tden.cube  | tail -n 4
```

示例目录中的 `run_cdd.sh` 把上述步骤（含下一步的差分）写成了一个 Slurm 作业。

### 7.4 add_gcube / diff_gcube

两个工具的源码在 `openmx4.0/source/`（手册说明在 work 目录中提供；若找不到可 `find ~/openmx4.0 -name diff_gcube.c`），用 gcc 单独编译：

```bash
cd ~/openmx4.0/source
gcc diff_gcube.c -lm -o diff_gcube
gcc add_gcube.c  -lm -o add_gcube
cp diff_gcube add_gcube ~/bin/            # or anywhere in PATH
```

```bash
add_gcube  co_cu_A.tden.cube  co_cu_B.tden.cube  co_cu_A_B.tden.cube
diff_gcube co_cu_AB.tden.cube co_cu_A_B.tden.cube co_cu_dAB.cube
```

- `add_gcube` 与 `diff_gcube` 会逐项比较两个文件头。A 与 B 原子数不同，会打印 `Found a difference in the number of atoms`——这是预期的、无害的提示（输出文件沿用第一个文件的原子列表，所以 `diff_gcube` 要把 AB 放在第一个）。
- 若打印 `Found a difference in the number of grid on a-axis` 或 `... x-coordinate of the origin`，说明网格没有对齐（忘了 `scf.fixed.grid`、晶胞或截断能不一致），得到的差分图是错误的，必须重算。
- 旧版本中这两个程序在发现差异后可能询问 "Are you sure you want to continue?, yes(1) or no(2)"；仅在确认只是原子数不同时回答 1。

同样的命令也可用于自旋密度（`sden.cube`）之差，或者同一体系在两个电场/两种自旋方向下的密度之差（手册第 66 章）。

### 7.5 平面平均与电荷转移量

示例目录提供的 `cdd_planar.py`（numpy；matplotlib 可选）沿 c 轴对差分密度做平面积分并累计：

```bash
python cdd_planar.py co_cu_dAB.cube --png cdd_z.png > cdd_z.dat
```

输出三列：z（Å）、$\Delta\lambda(z)$（e/Å，平面积分后的线密度）、$\Delta Q(z)$（e，从晶胞底部累计到 z 的电子数变化）。读图方法：

- $\Delta\lambda(z)>0$ 表示电子在该高度积累，<0 表示亏损；
- 在 Cu 顶层与 C 原子之间某个平面（常取二者中点）读 $\Delta Q$，其绝对值就是"跨过该平面"的净电子转移估计；整个晶胞积分应回到约 0（电荷守恒检查）；
- 对 CO/金属体系，经典的 Blyholder 图像是 5σ 向金属给电子、金属 d 向 CO 2π* 反馈；Δρ 等值面上通常能看到 CO 轴向（σ）区域亏损、垂直于轴的 π 区域积累，金属表面原子附近重新分布。具体以你的计算结果为准。

### 7.6 等值面

在 VESTA 中打开 `co_cu_dAB.cube`，Isosurface level 取 ±0.001 ～ ±0.003 e/Bohr³，正值（积累）与负值（亏损）用两种颜色；Boundary 设置为 x、y 方向各 0–2 可看到周期延拓。论文作图请在图注中写明等值面数值与单位（cube 中为 e/Bohr³，换算 e/Å³ 需乘以 6.748）。

---

## 8. 大文件：OutData.bin.flag、bin2txt 与 cube2xsf

上千原子、大盒子的计算中，每个 cube 文件可能有几百 MB，文本写盘本身就很慢。可以改为二进制输出：

```text
OutData.bin.flag    on        # default off
```

此时 cube 文件扩展名变为 `.cube.bin`（非共线的 `nccube` 同理）。转换回文本：

```bash
cd ~/openmx4.0/source
gcc bin2txt.c -lm -o bin2txt
gcc cube2xsf.c -o cube2xsf
cd /path/to/calc
bin2txt *.bin                     # -> text cube files
cube2xsf big.tden.cube            # cube -> periodic XSF (big.tden.xsf)
cube2xsf big.tden.cube.bin        # cube2xsf also reads .cube.bin directly
```

| 工具 | 输入 | 输出 | 用途 |
|---|---|---|---|
| `bin2txt` | `*.bin` | 文本 cube | 在 I/O 慢的机器上先写二进制、事后转换 |
| `cube2xsf` | `.cube` 或 `.cube.bin`（可多个） | `.xsf` | XCrySDen 标准格式，显式保留周期性；把周期体系的密度/轨道交给其它可视化软件时更稳妥 |

> [!TIP]
> 如果只需要能量和电荷，干脆用 `level.of.fileout 0` 不写 cube；需要密度时再用 `scf.restart pp` 读重启文件、只做一次对角化输出（[第5章](05-scf-kpoints-convergence.md)）。

---

## 9. 可视化：VESTA、XCrySDen、OpenMX Viewer

| 工具 | 适合 | 读入 | 本章用法 |
|---|---|---|---|
| [VESTA](https://jp-minerals.org/vesta/en/) | 等值面、截面、出版级作图 | `.cube`、`.xsf` | 密度、轨道、差分密度；Utilities → 2D Data Display 画截面 |
| [XCrySDen](http://www.xcrysden.org/) | 周期 XSF、矢量场 | `.xsf`、`.cube` | 非共线自旋矢量 `nc.xsf`/`ncsden.xsf`：Display → Forces 显示箭头 |
| [OpenMX Viewer](https://www.openmx-square.org/viewer/) | 浏览器端、免安装（另有离线 Windows 版） | `.cube`、`.xsf`、`.md`、`.Band`、`.bxsf` 等 | 快速查看等值面、用第二个 cube 给等值面着色、查看 MD/优化轨迹中的电荷与自旋 |

OpenMX Viewer（2026-05-26 大更新）支持在一个 cube 的等值面上用另一个兼容 cube 的数值着色——例如把 `v0.cube`（静电势）映射到 `tden.cube` 的等值面上得到"分子静电势图"，这对判断亲电/亲核位点很直观。拖放到浏览器的文件在本地处理，不上传服务器。

常用等值面数值（经验范围，按体系调整）：

| 数据 | 单位 | 建议等值面 |
|---|---|---|
| `tden.cube` | e/Bohr³ | 0.02–0.1（看分子外形）；0.15–0.3（看成键区域） |
| `dden.cube` | e/Bohr³ | ±0.005–±0.02 |
| 电荷密度差 dAB | e/Bohr³ | ±0.001–±0.003 |
| `sden.cube` | e/Bohr³ | ±0.002–±0.01 |
| 轨道 `homo/lumo` | Bohr^-3/2 | ±0.02–±0.05 |

更系统的出版级作图与配色技巧见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md)。

---

## 10. HS.fileout 与 scfout：后处理生态的入口

`HS.fileout on` 会把 SCF 收敛后的 Kohn–Sham 哈密顿量、重叠矩阵、密度矩阵（以及位置、动量算符矩阵元等）以二进制写入 `System.Name.scfout`：

```text
HS.fileout     on        # default off
```

4.0 的后处理程序读 scfout 时会打印 "The file format of the SCFOUT file:  3" 以及它支持的功能（jx、polB、kSpin、Z2FH、calB）。官方提供 `analysis_example` 演示如何用 `read_scfout.h` / `read_scfout()` 读取：

```bash
cd ~/openmx4.0/source && make analysis_example    # copied to ../work
~/openmx4.0/work/analysis_example h2o.scfout > HS.out
```

| 下游程序 | 需要的计算类型 | 用途 | 本教程章节 |
|---|---|---|---|
| `jx` | 共线自旋极化 | 交换耦合 J_ij | [第11章](11-spin-noncollinear.md) |
| [TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX) | 自旋极化（含 SOC） | 交换、DMI、各向异性交换、磁振子 | [第11章](11-spin-noncollinear.md) |
| `kSpin`、`calB`、`Z2FH` | 非共线（+SOC） | 自旋织构、Berry 曲率/Chern 数、Z2 | [第12章](12-soc-topology.md) |
| [OMXsigmaxy](https://github.com/hikaruri/OMXsigmaxy)、[HWC](https://github.com/Ncmexp2717/HWC) | 视功能而定 | 反常霍尔电导、混合 Wannier 中心 | [第12章](12-soc-topology.md) |
| `polB` | 任意 | Berry 相宏观极化 | [第15章](15-berry-efield.md) |
| [DeepH-pack](https://github.com/mzjb/DeepH-pack)、[HamGNN](https://github.com/QuantumLab-ZY/HamGNN) | 任意 | 机器学习哈密顿量数据集 | [第18章](18-parallel-on-errors.md) |
| [HopTB.jl](https://github.com/HopTB/HopTB.jl) | 任意 | 紧束缚模型、Berry 曲率、光学/非线性响应 | — |

> [!WARNING]
> scfout 是二进制文件，与 OpenMX 版本、编译器字节序相关。第三方工具（尤其是基于 3.9 源码开发的解析器）读取 4.0 生成的 scfout 前，请先用小体系验证；跨机器拷贝时注意字节序（[OMXTool](https://github.com/Ncmexp2717/OMXTool) 提供了 scfout 字节序转换）。scfout 文件大小与 k 点数无关，但随原子数、每个原子的基函数数以及轨道相互作用范围增长，大体系可达数 GB。

---

## 参数选择与收敛建议

| 量 | 推荐设置 | 收敛/检查方法 |
|---|---|---|
| `scf.energycutoff`（cube 质量） | 200–300 Ry；ESP、Bader、差分密度用上限 | 积分 tden 与价电子数差 < 0.01 e；差分密度图无网格噪声 |
| Cluster 盒子 | 分子边缘到盒子边 ≥ PAO 截断半径 + 3 Å | 增大盒子后 ESP 电荷与偶极基本不变 |
| 基组（Mulliken） | 同一研究中固定 PAO/VPS | Quick→Standard→Precise 看电荷**差值**是否稳定 |
| `scf.criterion` | 电荷/轨道分析 1e-8；差分密度 1e-8 以下 | 三个计算（AB/A/B）用同一标准 |
| 差分密度 | 同晶胞、同 `scf.energycutoff`、同 k 点、同电子温度、`scf.fixed.grid` | 三个 cube 的第 3–6 行完全一致；ΔQ 全晶胞积分 ≈ 0 |
| ESP 壳层 | `-c 0 -s 1.4 2.0` | 改为 `-s 1.2 2.0` 等看电荷是否稳健；RMS 足够小 |
| NPA | `NAO.threshold` 默认 0.85 | 占据 NAO 数与路易斯结构预期一致 |
| 轨道等值面 | ±0.02–0.05 | 相邻简并轨道一起画或画平方和 |
| 磁盘 I/O | 大体系 `OutData.bin.flag on` 或 `level.of.fileout 0` | 用 `bin2txt`/`cube2xsf` 事后转换 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 没有任何 `.cube` 文件 | `level.of.fileout 0` | 改为 1（或 2）；或用 `scf.restart pp` 补算 |
| `MO.fileout=ON is not supported in case of scf.EigenvalueSolver=DC` 后跟 `MO.fileout is changed to OFF` | O(N)/DC 类求解器不输出 MO | 小体系改用 Cluster/Band；大体系考虑 NBO `on2` 或 CWF 做局域分析 |
| `MO.Nkpoint should be 1 for the cluster calculation.` | Cluster 计算写了多个 MO k 点 | `MO.Nkpoint 1`，`<MO.kpoint` 只留 `0 0 0` |
| 只有 `homo0_*`，没有 `homo1_*` | 非自旋极化计算只有一个自旋通道 | 正常；自旋极化计算才有 `homo1_*`（down 自旋） |
| 周期体系轨道图"看起来不对称" | 只画了 `_r` 或 `_i`，受任意整体相位影响 | 画 `_r² + _i²` |
| `diff_gcube` 打印 `Found a difference in x-coordinate of the origin` | A/B 未使用 AB 的 `Grid_Origin` | 在 A、B 输入中加 `scf.fixed.grid` 并重算；数值原样复制（Bohr） |
| `diff_gcube` 打印 `Found a difference in the number of grid on c-axis` | 晶胞或 `scf.energycutoff` 不一致 | 三个输入用完全相同的晶胞与截断能（必要时都写同一个 `scf.Ngrid`） |
| 差分密度在原子核附近出现强烈正负交替 | A/B 的几何与 AB 不一致，或收敛不足 | 坐标逐位复制；提高 `scf.criterion` |
| `esp` 报 `Invalid argument` | 参数格式错误（例如把扩展名写进了名字） | `esp h2o -c 0 -s 1.4 2.0`，确保 `h2o.out` 与 `h2o.vhart.cube` 在当前目录 |
| ESP 电荷离谱、RMS 很大 | 盒子太小、分子外采样点不足，或体系是周期固体 | 增大盒子；ESP 只用于分子 |
| `esp` 找不到 `vhart.cube` | `level.of.fileout 0` | 用 1 或 2 重算 |
| Voronoi 电荷 H 带负电 | 分子中 Voronoi 多面体不考虑原子尺寸 | 属方法特性；分子改用 Mulliken/ESP/NPA |
| NPA 大体系内存不足 | `on1` 使用全密度矩阵 | Krylov 求解器 + `NBO.switch on2` + `<NBO.CenterAtoms` 只分析局部 |
| VESTA 打开大 cube 卡死 | 文件过大 | 先在集群上用 Python 截取子区域或降采样；或改看平面平均 |
| 积分 tden 远小于价电子数 | 分子跨越 Cluster 盒子边界被截断，或网格太粗 | 把分子放在盒子中央、增大盒子与截断能 |
| 非共线 + 仅 Γ 点（Cluster）时 HOMO 序号不对 | 4.0 中 `Mulliken_Charge.c` 的已知问题 | 打 patch 4.0.1（见 [第2章](02-install.md)） |

---

## 练习

1. **基组与电荷**：把 H2O 的基组从 O6.0-s2p2d1/H6.0-s2p1 换成 Precise（O6.0-s3p2d2/H7.0-s2p2d1），比较 Mulliken、Voronoi、ESP、NPA 四种电荷的变化幅度，哪一种对基组最敏感？
2. **积分检查**：对 `h2o.tden.cube`、`h2o.dden.cube` 分别积分，验证前者约等于 8、后者约等于 0；把 `scf.energycutoff` 从 250 Ry 降到 120 Ry 再试一次。
3. **苯的 π 轨道**：画出 HOMO−1、HOMO、LUMO、LUMO+1，并计算 $|\psi_{14}|^2+|\psi_{15}|^2$，验证其六重对称性；在 `.out` 中找出它们的 LCAO 系数主要分量。
4. **吸附位点对比**：把 CO 移到 fcc 空位（three-fold hollow），重复电荷密度差流程，比较 Cu 顶层与 C 之间平面上的 ΔQ。
5. **自旋密度**：对 O2 分子（三重态，总自旋磁矩 2 μB：每个 O 的初始占据写 3.5 2.5）做自旋极化计算，画 `sden.cube`，并比较 Mulliken 与 Voronoi 的自旋磁矩。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 9 章 Output files：<https://www.openmx-square.org/openmx_man4.0/s9_output-files.html>
- 第 21 章 Visualization：<https://www.openmx-square.org/openmx_man4.0/s21_visualization.html>
- 第 33 章 LCAO coefficients：<https://www.openmx-square.org/openmx_man4.0/s33_lcao-coeff.html>
- 第 34 章 Molecular orbitals：<https://www.openmx-square.org/openmx_man4.0/s34_molecular-orb.html>
- 第 35 章 Charge analysis（Mulliken / Voronoi / ESP / CWF）：<https://www.openmx-square.org/openmx_man4.0/s35_charge-analysis.html>
- 第 36 章 Natural population analysis：<https://www.openmx-square.org/openmx_man4.0/s36_npa.html>
- 第 66 章 diff_gcube：<https://www.openmx-square.org/openmx_man4.0/s66_cube-diff.html>
- 第 68 章 相互作用引起的电荷密度差：<https://www.openmx-square.org/openmx_man4.0/s68_charge-diff.html>
- 第 70 章 scfout 开发者接口：<https://www.openmx-square.org/openmx_man4.0/s70_dev-interface.html>
- 第 76 章 二进制输出、第 77 章 cube→XSF：<https://www.openmx-square.org/openmx_man4.0/s76_binary-output.html>，<https://www.openmx-square.org/openmx_man4.0/s77_cube-to-xsf.html>
- OpenMX Viewer：<https://www.openmx-square.org/viewer/>
- Patch 4.0.1 说明：<https://www.openmx-square.org/bugfixed/26May08/README.txt>
- Henkelman 组 Bader 程序：<https://theory.cm.utexas.edu/henkelman/code/bader/>

**相关 GitHub 仓库（仅链接，未复制代码）**

- [mzjb/DeepH-pack](https://github.com/mzjb/DeepH-pack)：用 OpenMX 哈密顿量（HS 输出）训练深度学习哈密顿量的代表性下游工具。
- [QuantumLab-ZY/HamGNN](https://github.com/QuantumLab-ZY/HamGNN)：支持 OpenMX 数据的等变图神经网络哈密顿量，附 OpenMX 后处理流程。
- [pulkin/openmx-hks](https://github.com/pulkin/openmx-hks)：读取 `.hks` 并导出哈密顿量/重叠矩阵与结构，了解二进制格式可参考。
- [nnardoto/LowFoldIn](https://github.com/nnardoto/LowFoldIn)：原生读取 3.9/4.0 的 scfout 做能带与能带展开。
- [Ncmexp2717/OMXTool](https://github.com/Ncmexp2717/OMXTool)：Tcl/Tk 图形工具，含 scfout 字节序转换。
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：基于 OpenMX 哈密顿量的紧束缚后处理（Berry 曲率、光学响应）。
- [samlipton/pyfesm](https://github.com/samlipton/pyfesm)：轻量 Python 解析器（SCF、Mulliken 等），可借鉴输出解析思路。
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，含轨道投影与磁矩可视化示例。

**关键文献**

- R. S. Mulliken, J. Chem. Phys. 23, 1833 (1955). DOI: [10.1063/1.1740588](https://doi.org/10.1063/1.1740588)
- A. D. Becke, J. Chem. Phys. 88, 2547 (1988)（模糊胞划分）. DOI: [10.1063/1.454033](https://doi.org/10.1063/1.454033)
- U. C. Singh, P. A. Kollman, J. Comput. Chem. 5, 129 (1984)（ESP 拟合）. DOI: [10.1002/jcc.540050204](https://doi.org/10.1002/jcc.540050204)
- A. E. Reed, R. B. Weinstock, F. Weinhold, J. Chem. Phys. 83, 735 (1985)（NPA）. DOI: [10.1063/1.449486](https://doi.org/10.1063/1.449486)
- T. Ohwaki, M. Otani, T. Ozaki, J. Chem. Phys. 140, 244105 (2014)（OpenMX 中大体系 NBO/NAO 方法）. DOI: [10.1063/1.4884119](https://doi.org/10.1063/1.4884119)
- G. Henkelman, A. Arnaldsson, H. Jónsson, Comput. Mater. Sci. 36, 354 (2006)（Bader 算法）. DOI: [10.1016/j.commatsci.2005.04.010](https://doi.org/10.1016/j.commatsci.2005.04.010)
- K. Momma, F. Izumi, J. Appl. Crystallogr. 44, 1272 (2011)（VESTA 3）. DOI: [10.1107/S0021889811038970](https://doi.org/10.1107/S0021889811038970)
- A. Kokalj, Comput. Mater. Sci. 28, 155 (2003)（XCrySDen）. DOI: [10.1016/S0927-0256(03)00104-6](https://doi.org/10.1016/S0927-0256(03)00104-6)
- T. Ozaki, Phys. Rev. B 67, 155108 (2003)（OpenMX 的 PAO 变分优化）. DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)

[← 上一章](09-dos-pdos-cohp-optics.md) | [返回目录](../README.md) | [下一章 →](11-spin-noncollinear.md)
