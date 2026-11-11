# 第12章 自旋轨道耦合与拓扑性质（Spin-Orbit Coupling & Topological Properties with OpenMX）

本章讲解如何在 OpenMX 中进行自旋轨道耦合（SOC）计算并分析能带拓扑：非共线框架下的 `scf.SpinOrbit.Coupling`、DFT_DATA19 中 j 依赖赝势的全相对论与标量相对论处理、`scf.SO.factor` 调节 SOC 强度；以 GaAs（Γ/L 点劈裂）、单层 MoS2（K 谷自旋劈裂）和 Bi2Se3（能带反转）为例做定量验证；再用 `HS.fileout` 生成的 scfout 配合 kSpin（自旋织构）、calB（Berry 曲率与 Chern 数，FHS 方法）、Z2FH（Z2 不变量）以及 OMXsigmaxy、HWC 等工具研究拓扑性质，并介绍 Rashba 体系与费米面（bxsf/frmsf）的输出（OpenMX spin-orbit coupling, band inversion, Z2 invariant, Chern number, Berry curvature, spin texture tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 正确设置 SOC 计算（`scf.SpinPolarization nc` + `scf.SpinOrbit.Coupling on`），理解为什么共线计算不能加 SOC；
- 用同一套 j 依赖赝势只切换 `scf.SpinOrbit.Coupling on/off`，做全相对论与标量相对论的严格对比；用 `scf.SO.factor` 按元素、按角动量缩放 SOC；
- 从 `.Band` 文件定量读出 SOC 劈裂（GaAs 的 Δ0、MoS2 的 K 点价带劈裂），并与文献/实验对照；
- 通过 SOC 开/关和 SOC 强度扫描识别 Bi2Se3 的能带反转，用 Z2FH 得到 Z2 = (1; 000)；
- 用 kSpin 计算动量空间自旋织构，用 calB 计算 Berry 曲率与 Chern 数；
- 知道反常霍尔电导、Wannier 电荷中心、Rashba 参数、费米面等分析该用哪些工具和输出文件。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| 重元素（Bi、Pb、Pt、W、I…）的能带 | nc + SOC，必要时与 SOC off 对比 | §3 |
| 半导体价带顶劈裂、自旋分裂 | `.Band` + `soc_splitting.py` | §3、§4 |
| 无反演中心体系的自旋劈裂、能谷物理 | SOC 能带 + kSpin | §4、§7 |
| 拓扑绝缘体判定 | SOC 能带 + Z2FH（非磁） | §5、§6 |
| 磁性拓扑、量子反常霍尔 | calB（Chern 数）、OMXsigmaxy | §8、§9 |
| Rashba 表面、费米面自旋织构 | kSpin FermiLoop、FermiSurfer 输出 | §7、§10 |
| 磁各向异性、轨道磁矩 | 见 [第11章](11-spin-noncollinear.md) §8、§9 | — |
| 表面态、Wannier 紧束缚模型 | 见 第14章 | — |

> [!NOTE]
> 本章以 [第8章](08-band.md) 的能带计算、[第11章](11-spin-noncollinear.md) 的非共线输入为基础；scfout 文件的通用说明见 [第10章](10-charge-mo-cube.md) §10。VASP 中的 SOC 设置可对照 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)。

---

## 2. 原理速览

### 2.1 SOC 与 j 依赖赝势

原子中的自旋轨道耦合 $H_{\rm SO}=\xi(r)\,\mathbf L\cdot\mathbf S$ 把同一 l 的态劈裂为 $j=l\pm 1/2$ 两组。OpenMX 使用的赝势由 ADPACK 求解原子 Dirac 方程生成，对每个 l 分别给出 $j=l+1/2$ 与 $j=l-1/2$ 两个通道（j 依赖赝势），质量–速度项、Darwin 项与 SOC 都包含在这一框架内：

$$
\hat V_{\rm NL}=\sum_{l}\sum_{j=l\pm1/2}\sum_{m_j}|\,l j m_j\rangle\, V_{lj}\,\langle l j m_j\,| .
$$

SOC 不是作为微扰后加的，而是通过把二分量旋量的 α、β 分量直接耦合起来自洽地包含在 Kohn–Sham 方程中。若关闭 SOC（`scf.SpinOrbit.Coupling off`），OpenMX 读入 j 依赖赝势时会按简并度自动平均

$$
V_l^{\rm SR}=\frac{(l+1)\,V_{l,\,l+1/2}+l\,V_{l,\,l-1/2}}{2l+1},
$$

得到标量相对论赝势。因此**同一套 DFT_DATA19 赝势既可做全相对论也可做标量相对论计算**，二者的差别只来自 SOC 本身，这是做对比的理想条件。

> [!IMPORTANT]
> 共线计算（`scf.SpinPolarization off/on`）只能做标量相对论处理；SOC 只能在非共线框架（`nc`）下打开。OpenMX 直接对价电子使用 Schrödinger 型方程，价电子的直接相对论效应被忽略，但主要的间接效应（经由芯电子）已包含在赝势中。

### 2.2 对称性、Kramers 简并与自旋劈裂

- **时间反演对称（T）**：$E_\uparrow(\mathbf k)=E_\downarrow(-\mathbf k)$。在时间反演不变动量点（TRIM，如 Γ、L、M）每个能级至少二重简并（Kramers 简并）。
- **空间反演（P）+ T**：整个布里渊区所有能带保持二重简并，SOC 只会劈裂轨道简并（如 GaAs 若有反演中心时的 Γ8/Γ7）。
- **缺少反演中心**（闪锌矿、单层 MoS2、表面、界面）：SOC 使一般 k 点上的自旋简并解除——体相的 Dresselhaus 型、表面/界面的 Rashba 型。单层 MoS2 在 K、K' 的价带劈裂约 0.15 eV，两谷自旋相反（自旋–能谷锁定）。

### 2.3 Berry 相、Chern 数与 Z2

Bloch 态 $|u_{n\mathbf k}\rangle$ 的 Berry 联络与 Berry 曲率为

$$
\mathbf A_n(\mathbf k)=-i\langle u_{n\mathbf k}|\nabla_{\mathbf k}|u_{n\mathbf k}\rangle,\qquad
\mathbf F_n=\nabla_{\mathbf k}\times\mathbf A_n .
$$

二维（或三维中的一个 k 平面）所有占据带的 Chern 数 $C=\frac{1}{2\pi}\sum_n^{\rm occ}\int F_{n,z}\,d^2k$ 是整数，霍尔电导 $\sigma_{xy}=-\frac{e^2}{h}C$。时间反演对称体系 C 必为 0，此时用 Z2 不变量区分平庸（0）与拓扑（1）绝缘体；三维体系由 6 个时间反演不变平面定义 4 个独立的 Z2 指标 $(\nu_0;\nu_1\nu_2\nu_3)$，$\nu_0=1$ 为强拓扑绝缘体。

OpenMX 的 calB 与 Z2FH 采用 Fukui–Hatsugai–Suzuki（FHS）方法：在离散 k 网格上用相邻 k 点的重叠矩阵行列式 $U_{\Delta\mathbf k}(\mathbf k)=\det\langle u_n(\mathbf k)|u_m(\mathbf k+\Delta\mathbf k)\rangle$ 计算每个小方格（plaquette）上的 Berry 曲率 $F={\rm Im}\log U_{12}U_{23}U_{34}U_{41}$，再求和得到 Chern 数；Z2FH 在半个布里渊区上做积分，并在边界上按平移、时间反演、Kramers 简并三种方式修正规范。这类方法即使网格较粗也能给出严格的整数。

---

## 3. SOC 计算的基本设置与 GaAs 验证（GaAs_SOC_band.dat）

完整文件见 [../examples/12-gaas-soc-band/README.md](../examples/12-gaas-soc-band/README.md)。

### 3.1 必要关键词

| 关键词 | 取值 | 说明 |
|---|---|---|
| `scf.SpinPolarization` | `nc` | SOC 的前提。非磁体系坐标行写 7 列（上下占据相等）即可，Euler 角列可省略 |
| `scf.SpinOrbit.Coupling` | `on` | 默认 off；off 时 j 依赖赝势自动平均为标量相对论 |
| `scf.XcType` | GGA-PBE / LSDA-CA / LSDA-PW | 非共线计算不能用非自旋的 `LDA` |
| VPS | DFT_DATA19 的 `*_PBE19`、`*_CA19` | 均为 j 依赖赝势，可直接用于 SOC；PAO 与 VPS 必须同属 DFT_DATA19 |
| `<scf.SO.factor` | 可选 | 按 species 与 l 缩放 SOC 强度（§3.4） |

相关报错：共线 + SOC 报 `Spin-orbit coupling is not supported for collinear DFT calculations.`；`scf.SpinPolarization off` + SOC 报 `scf.SpinOrbit.Coupling must be OFF when scf.SpinPolarization=OFF`。

### 3.2 GaAs 输入文件

```text
System.CurrrentDirectory    ./
System.Name                 gaas_soc
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            0

Species.Number              2
<Definition.of.Atomic.Species
  Ga   Ga7.0-s3p2d2   Ga_PBE19
  As   As7.0-s3p2d2   As_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Ga   0.00  0.00  0.00   6.5  6.5
  2  As   0.25  0.25  0.25   7.5  7.5
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  2.8265   2.8265   0.0000
  2.8265   0.0000   2.8265
  0.0000   2.8265   2.8265
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on
scf.EigenvalueSolver        Band
scf.Kgrid                   8 8 8
scf.energycutoff            200.0
scf.ElectronicTemperature   300.0
scf.maxIter                 200
scf.criterion               1.0e-9
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.40
scf.Mixing.History          30
scf.Mixing.StartPulay       6
MD.Type                     Nomd

Band.dispersion             on
Band.Nkpath                 4
<Band.kpath
  30  0.500 0.500 0.500   0.000 0.000 0.000   L  G
  35  0.000 0.000 0.000   0.500 0.500 0.000   G  X
  12  0.500 0.500 0.000   0.625 0.625 0.250   X  U
  37  0.750 0.375 0.375   0.000 0.000 0.000   K  G
Band.kpath>
```

要点：

- a = 5.653 Å（实验值）；原胞基矢顺序为 a/2(1,1,0)、a/2(1,0,1)、a/2(0,1,1)，此时 L = (½,½,½)、X = (½,½,0)、U = (⅝,⅝,¼)、K = (¾,⅜,⅜)（以原胞倒格矢为基）。换了基矢顺序就要重新推导，或改用 `<Band.KPath.UnitCell`（[第8章](08-band.md)）；
- Ga_PBE19（13 个价电子，含 3d）与 As_PBE19（15 个价电子）都含 d 半芯态，占据写 6.5/6.5 与 7.5/7.5；
- 绝缘体用 300 K 电子温度即可；`scf.criterion` 取 1e-9 让本征值足够精确；
- 无 SOC 对照文件 `GaAs_noSOC_band.dat` 只把 `scf.SpinOrbit.Coupling` 改为 `off`（仍是 `nc`），保证两者的基组、哈密顿量维度、k 点完全一致。

### 3.3 运行与定量分析

```bash
mpirun -np 8 openmx GaAs_SOC_band.dat   -nt 1 > gaas_soc.std   2>&1
mpirun -np 8 openmx GaAs_noSOC_band.dat -nt 1 > gaas_nosoc.std 2>&1
python soc_splitting.py gaas_soc.Band   --at G L
python soc_splitting.py gaas_nosoc.Band --at G L
```

`soc_splitting.py`（示例目录提供）直接读 `.Band`（Hartree，未平移），换算到以化学势为零点的 eV，并在指定标签的 k 点按简并度分组。输出形如：

```text
  k = G [0. 0. 0.]   (energy eV, degeneracy)
        0.5xxx  x2
     ---------- E_F ----------
       -0.0xxx  x4
       -0.3xxx  x2
     splitting of the two highest valence groups : 0.3xxx eV
```

（示意输出，数值仅供参考）

判断标准：

| k 点 | 无 SOC | 有 SOC | 参考值 |
|---|---|---|---|
| Γ 价带顶 | 6 重简并（Γ15 × 自旋） | 4 重 Γ8 + 下方 2 重 Γ7 | Δ0：实验约 0.34 eV；手册表 6（LDA）0.344 eV |
| L 价带顶 | 4 重简并 | 两组 2 重 | 手册表 6：L3v 劈裂 0.213 eV |
| 带隙 | — | 比无 SOC 小约 Δ0/3（Γ8 上移） | GGA 严重低估绝对带隙（实验 1.52 eV） |

用 `bandgnu13 gaas_soc.Band` 生成 `gaas_soc.GNUBAND`（能量已换算为 eV、以化学势为零点），或用 `../scripts/omx_band_plot.py`（`-h` 查看用法）把 SOC 与无 SOC 能带画在同一张图上。

> [!TIP]
> 若想复现手册表 6 的数值，可参照手册设置（LDA：`scf.XcType LSDA-CA` 配 `*_CA19` 赝势，Ga7.0-s2p2d2 / As7.0-s2p2d2，12³ k 点，140 Ry）。patch 4.0.1 修正了官方 work 目录里的 `GaAs.dat`，正是为了让手册 38.1 节的结果可以复现。

### 3.4 调节 SOC 强度：scf.SO.factor

OpenMX（3.9 起）可以直接用 DFT_DATA19 的赝势按角动量缩放 SOC，无需生成特殊赝势：

```text
<scf.SO.factor
  Ga  s 1.0 p 2.0 d 1.0 f 1.0
  As  s 1.0 p 2.0 d 1.0 f 1.0
scf.SO.factor>
```

- 第 1 列是 `Definition.of.Atomic.Species` 中的 species 名，其后依次是角动量符号与缩放因子；1.0 对应真实原子的 SOC 强度；
- **必须为体系中所有 species 写满 s、p、d、f 四项**，即使该元素没有 f 基函数；
- 典型用途：(1) 把某元素的 p 因子从 0 连续调到 1，观察拓扑相变中能隙的闭合与重开（§5.3）；(2) 分辨哪个元素/哪个轨道的 SOC 主导了某个劈裂（例如只关掉 Se 的 SOC 看 Bi2Se3 的变化）；(3) 构造"放大 SOC"的模型体系做教学演示。

### 3.5 成本与收敛

- 非共线计算的哈密顿量维度是共线的 2 倍，稠密对角化成本约 8 倍，内存约 4 倍；k 点通常不能利用自旋简并；
- SOC 引起的能量与劈裂通常在 meV–百 meV 量级，`scf.criterion` 建议 1e-8～1e-10；
- 对重元素，基组中的 p、d（乃至 f）极化函数影响 SOC 劈裂，建议在 Standard 与 Precise 之间做一次对比；
- 先做无 SOC 的非共线（或共线）收敛测试确定截断能与 k 点，再打开 SOC，可以节省大量机时。

---

## 4. 单层 MoS2：无反演中心的自旋劈裂（MoS2_SOC_band.dat）

完整文件见 [../examples/12-mos2-soc-band/README.md](../examples/12-mos2-soc-band/README.md)。

### 4.1 结构与路径

六角晶胞取 a1 = (a, 0, 0)、a2 = (−a/2, √3a/2, 0)（夹角 120°），c = 20 Å 提供真空；Mo 在 (0, 0, ½)，两个 S 在 (⅓, ⅔, ½ ± Δz)。此约定下 M = (½, 0, 0)、K = (⅓, ⅓, 0)、K' = −K。输入的关键部分：

```text
Species.Number              2
<Definition.of.Atomic.Species
  Mo   Mo7.0-s3p2d2     Mo_PBE19
  S    S7.0-s2p2d1f1    S_PBE19
Definition.of.Atomic.Species>

Atoms.Number                3
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Mo   0.0000000000  0.0000000000  0.5000000   7.0  7.0
  2  S    0.3333333333  0.6666666667  0.5782500   3.0  3.0
  3  S    0.3333333333  0.6666666667  0.4217500   3.0  3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
   3.160000   0.000000   0.000000
  -1.580000   2.736640   0.000000
   0.000000   0.000000  20.000000
Atoms.UnitVectors>

scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on
scf.Kgrid                   12 12 1
HS.fileout                  on

Band.dispersion             on
Band.Nkpath                 4
<Band.kpath
  40  0.0000000000  0.0000000000  0.0   0.5000000000  0.0000000000  0.0   G   M
  25  0.5000000000  0.0000000000  0.0   0.3333333333  0.3333333333  0.0   M   K
  45  0.3333333333  0.3333333333  0.0   0.0000000000  0.0000000000  0.0   K   G
  45  0.0000000000  0.0000000000  0.0  -0.3333333333 -0.3333333333  0.0   G  Kp
Band.kpath>
```

（几何为近似实验值，未优化；Mo_PBE19 有 14 个价电子、S_PBE19 有 6 个。）

### 4.2 结果解读

```bash
python ../12-gaas-soc-band/soc_splitting.py mos2_soc.Band --at G K Kp
```

- K 点价带顶分裂为两个非简并能级，劈裂约 0.15 eV（文献 GGA 值约 148 meV，Zhu et al. 2011）；K' 处劈裂大小相同；
- Γ 是 TRIM 点，Kramers 简并使每个能级仍为二重；
- 导带底在 K 点也有很小（几 meV 量级）的劈裂，对计算参数敏感；
- K 与 K' 处同一能量的态自旋相反，可用 kSpin 验证（§7）。这一"自旋–能谷锁定"是能谷电子学的基础（Xiao et al. 2012）。

> [!NOTE]
> 若 K 点劈裂明显偏小，先检查是否真的打开了 SOC（`.out` 开头回显的输入）、基组是否包含足够的 Mo d 函数；若明显偏大，检查晶格常数与 Mo–S 距离是否合理。

---

## 5. Bi2Se3：能带反转（Bi2Se3_SOC_band.dat）

完整文件见 [../examples/12-bi2se3-soc-z2/README.md](../examples/12-bi2se3-soc-z2/README.md)。

### 5.1 结构

Bi2Se3 属 R-3m，菱方原胞 5 个原子：Bi 在 ±(u,u,u)，Se1 在 (0,0,0)，Se2 在 ±(v,v,v)。示例取文献常用的近似实验值 a = 4.138 Å、c = 28.64 Å、u = 0.4008、v = 0.2117，菱方基矢 (a/2, −a/(2√3), c/3)、(0, a/√3, c/3)、(−a/2, −a/(2√3), c/3)：

```text
Species.Number              2
<Definition.of.Atomic.Species
  Bi   Bi8.0-s3p2d2f1   Bi_PBE19
  Se   Se7.0-s3p2d2     Se_PBE19
Definition.of.Atomic.Species>

Atoms.Number                5
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Bi   0.4008  0.4008  0.4008   7.5  7.5
  2  Bi   0.5992  0.5992  0.5992   7.5  7.5
  3  Se   0.0000  0.0000  0.0000   3.0  3.0
  4  Se   0.2117  0.2117  0.2117   3.0  3.0
  5  Se   0.7883  0.7883  0.7883   3.0  3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
   2.069000  -1.194538   9.546667
   0.000000   2.389075   9.546667
  -2.069000  -1.194538   9.546667
Atoms.UnitVectors>

scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on
scf.Kgrid                   8 8 8

Band.dispersion             on
Band.Nkpath                 3
<Band.kpath
  30  0.5 0.5 0.5   0.5 0.5 0.0   Z  F
  40  0.5 0.5 0.0   0.0 0.0 0.0   F  G
  40  0.0 0.0 0.0   0.5 0.0 0.0   G  L
Band.kpath>
```

菱方倒格矢基下 Z = (½,½,½)、F = (½,½,0)、L = (½,0,0)，与 Zhang et al.（Nat. Phys. 2009）使用的标记一致。总价电子数 2×15 + 3×6 = 48，所以非共线计算中第 48 条带是价带顶、第 49 条是导带底。

### 5.2 SOC 开/关对比

```bash
mpirun -np 16 openmx Bi2Se3_SOC_band.dat -nt 1 > bi2se3_soc.std 2>&1
sed -e 's/^System.Name .*/System.Name   bi2se3_nosoc/' \
    -e 's/^scf.SpinOrbit.Coupling .*/scf.SpinOrbit.Coupling   off/' \
    Bi2Se3_SOC_band.dat > Bi2Se3_noSOC_band.dat
mpirun -np 16 openmx Bi2Se3_noSOC_band.dat -nt 1 > bi2se3_nosoc.std 2>&1
python ../12-gaas-soc-band/soc_splitting.py bi2se3_soc.Band   --at G Z F L
python ../12-gaas-soc-band/soc_splitting.py bi2se3_nosoc.Band --at G Z F L
```

看什么：

- 有 SOC：带隙约 0.3 eV 量级（实验约 0.3 eV），带边在 Γ 附近；Γ 附近最高价带呈"驼峰/M 形"——这是 Bi-p 与 Se-p 特征在 Γ 处互换（反转）后，与远离 Γ 处正常次序的能带杂化的结果；
- 无 SOC：Γ 处为正常次序（导带底以 Bi-p 为主、价带顶以 Se-p 为主），能隙与能带形状明显不同；
- 能带形状只是"线索"，严格判断要靠 Z2 不变量（§6）或宇称/Wannier 电荷中心分析（§9）。

### 5.3 SOC 强度扫描：能隙闭合与重开

在 SOC 计算的输入中加入（每个 species 写满 s–f）：

```text
<scf.SO.factor
  Bi  s 1.0 p 0.5 d 1.0 f 1.0
  Se  s 1.0 p 1.0 d 1.0 f 1.0
scf.SO.factor>
```

把 Bi 的 p 因子依次取 0.0、0.25、0.5、0.75、1.0（每个值一个 System.Name），用 `soc_splitting.py --at G` 记录 Γ 点能隙（第 48、49 条带之差）。拓扑相变的特征是：Γ 点能隙随 SOC 增大先减小、在某个临界值附近闭合、然后重新打开；在两侧分别计算 Z2，应从 0 变为 1。这个"数值实验"比单看最终能带更有说服力。

### 5.4 Γ 点的轨道成分（可选）

要直接看 Γ 点价带顶/导带底的原子成分，可用 kSpin 的 MulPOnly 模式（§7）：在 Γ 点只算第 47–50 条带，输出按原子与 s/p/d 分解的 Mulliken 布居（`AtomMulP`、`MulP_p` 等），再用 MulPCalc 提取。Γ 点坐标为 (0, 0, 0)，k 点文件不涉及单位问题：

```text
4
0.0 0.0 0.0 47
0.0 0.0 0.0 48
0.0 0.0 0.0 49
0.0 0.0 0.0 50
```

对比 SOC 开/关两种情况下第 48 与第 49 条带上 Bi 与 Se 的 p 成分，即可直接看到"反转"。

---

## 6. Z2 不变量：Z2FH（Bi2Se3_Z2.dat）

### 6.1 流程

1. 非共线 + SOC 的 SCF，打开 `HS.fileout on` 生成 scfout（示例中按官方 `Bi2Se3-Z2.dat` 的做法写了 `scf.Generation.Kpoint Regular`）；
2. 编译 Z2FH：`cd ~/openmx4.0/source && make Z2FH`（复制到 work）；
3. 运行：

```bash
mpirun -np 16 openmx Bi2Se3_Z2.dat -nt 1 > bi2se3_z2.std 2>&1
mpirun -np 4 ~/openmx4.0/work/Z2FH bi2se3_z2.scfout < Z2FH.in > Z2FH.out
cat Z2.dat
```

不带输入重定向运行时，Z2FH 会交互地依次询问以下三项；`Z2FH.in` 就是这三行的答案：

| 行 | 本例 | 含义 |
|---|---|---|
| 1 | `5` | 半个积分区域（即布里渊区的四分之一）上每个方向的网格数；5 表示 5×5 的方格 |
| 2 | `1` | 1 = 计算全部 6 个时间反演不变面 (x0, xπ, y0, yπ, z0, zπ)；0 = 只算 4 个独立面 |
| 3 | `1` | 从第几个面开始（1–6 对应 k1=0、k1=G1/2、k2=0、k2=G2/2、k3=0、k3=G3/2），用于中断后续算 |

### 6.2 输出与判断

| 文件 | 内容 |
|---|---|
| `Z2.dat` | 第 1 行：6 个面的 Z2；第 2 行：(nu0, nu1, nu2, nu3) |
| `LCNum*.dat` / `LCNum*.pl` | 各面上整数场 n(k)（取 0、±1）及 gnuplot 脚本（`gnuplot LCNum1.pl`） |
| `temporal_12345.input` | 从 scfout 中恢复的 SCF 输入副本 |

Bi2Se3 应得到 (ν0; ν1ν2ν3) = (1; 000)，即强拓扑绝缘体。手册给出的格式示例为

```text
Z2 invariant:(x0,xpi,y0,ypi,z0,zpi)=(1.000000,-0.000000,1.000000,-0.000000,1.000000,-0.000000)
Z2 invariant:(nu0,nu1,nu2,nu3):(1,0,0,0)
```

（摘自 OpenMX 4.0 手册第 59.2 节）

注意事项：

- 6 个面的结果必须满足 $x_0+x_\pi=y_0+y_\pi=z_0+z_\pi\ ({\rm mod}\ 2)$，否则说明网格太粗或带隙未打开；
- n(k) 本身依赖规范，不同机器/编译选项下可能不同，但 Z2 应可重复；
- 把网格数从 5 加到 8、10，确认 Z2 不变；
- Z2FH 只适用于**非磁**的非共线计算（时间反演对称）；磁性体系没有 Z2，用 Chern 数（§8）；
- 体系必须是绝缘体（至少在每个 k 点上第 N 与第 N+1 条带之间有直接能隙）；
- 官方 work 目录中 `Bi2Se3_trivial-Z2.dat`、`Bi2Se3_weak-Z2.dat` 只改变了内坐标，就分别得到平庸相和弱拓扑相；`Bi111-Z2.dat` 是二维拓扑绝缘体 Bi(111) 薄膜的例子——拓扑相对结构非常敏感，务必使用可靠结构。

---

## 7. 自旋织构：kSpin

kSpin（Kotaka、Yamaguchi、Ishii 开发）从非共线 + SOC 的 scfout 计算每个本征态、每个 k 点的 2×2 动量空间自旋密度矩阵，得到自旋织构 ⟨σx⟩、⟨σy⟩、⟨σz⟩，并可分解到原子与轨道。只能用于非共线计算。

### 7.1 编译与运行

```bash
cd ~/openmx4.0/source && make kSpin MulPCalc     # copy to ../work if not copied automatically
# SCF with: scf.SpinPolarization NC, scf.SpinOrbit.Coupling ON, HS.fileout ON
./kSpin input_with_kspin_keywords.dat
mpirun -np 4 ./kSpin input_with_kspin_keywords.dat   # MPI over k points
```

kSpin 的关键词写在 .dat 文件中（可以与 OpenMX 输入放在同一个文件，openmx 会忽略它们）：

| 关键词 | 说明 |
|---|---|
| `Filename.scfout` | 读取的 scfout 文件名 |
| `Filename.outdata` | 输出文件前缀 |
| `Calc.Type` | `FermiLoop` / `GridCalc` / `BandDispersion` / `MulPOnly`（默认） |
| `Energy.Range` | 能量范围（eV，相对化学势）。FermiLoop 中两值相同，表示等能面；其它模式为搜索能带的窗口 |
| `Search.kCentral` | 搜索中心 k 点（写法同 `Band.kpath`） |
| `Calc.Type.3mesh` | 计算平面：1 = ka–kb，2 = kb–kc，3 = kc–ka |
| `kRange.3mesh` | 平面上的 k 范围（两个半宽） |
| `k-plane.1stStep` / `k-plane.2ndStep` | 两级搜索的网格划分 |
| `Calc.Bandbyband`、`Calc.Band.Min`、`Calc.Band.Max` | 直接指定能带序号范围 |
| `MulP.Vec.Scale` | 作图时自旋矢量的缩放 |
| `Filename.kpointdata` | MulPOnly 模式的 k 点文件（第 1 行为 k 点数，其后每行 kx ky kz 能带序号） |

| 模式 | 用途 | 主要输出 |
|---|---|---|
| FermiLoop | 等能面（如费米面）上的自旋织构，典型用于 Rashba | `*.FermiSurf_YY`（k 与能量）、`*.Pxyz_YY`（⟨σ⟩）、`*.plotexample`（gnuplot） |
| GridCalc | 二维 k 网格上指定能量窗内各能带的自旋 | 同上类文件 |
| BandDispersion | 沿 `Band.kpath` 路径的原子分辨自旋密度矩阵 | `*.BandXX_s`、`*.AMulPBand*` |
| MulPOnly | 用户指定 k 点与能带 | `*.AtomMulP`、`*.MulP_s/p/d...` |

MulPCalc 从 AtomMulP / MulP_xx / AMulPBand 文件中提取指定原子、指定轨道的数据（例如只看表面一侧原子的贡献）。

### 7.2 示例

- **MoS2（本仓库）**：`MoS2_SOC_band.dat` 末尾已写好 BandDispersion 所需关键词，SCF 后运行 `./kSpin MoS2_SOC_band.dat`，沿 Γ–M–K–Γ–K' 得到各能带的自旋期望值；K 与 K' 附近价带顶 ⟨σz⟩ 符号相反。
- **Au(111) 表面 Rashba（官方 work 目录）**：`Au111Surface_FL.dat`（FermiLoop）、`Au111Surface_GC.dat`（GridCalc）、`Au111Surface_BD.dat`（BandDispersion）、`Au111Surface_MO.dat`（MulPOnly）四个输入覆盖全部模式；FermiLoop 完成后 `gnuplot Au111Surface_FL.plotexample` 即得手册图 78 那样的 Rashba 自旋织构。

`Pxyz_YY` 文件的前三列是 kx、ky、kz（Bohr⁻¹），后三列是 ⟨σx⟩、⟨σy⟩、⟨σz⟩；`FermiSurf_YY` 的第 4、5 列是能带序号与能量（eV）。

---

## 8. Chern 数与 Berry 曲率：calB

### 8.1 流程

```bash
cd ~/openmx4.0/source && make calB
# SCF: scf.SpinPolarization NC (SOC on/off as needed), HS.fileout on
./calB System.Name.scfout                         # interactive
./calB System.Name.scfout < calB.in > calB.out
mpirun -np 4 ./calB System.Name.scfout < calB.in > calB.out
```

`calB.in`（四行）：

```text
1
0
0 0 1
100 100
```

| 行 | 含义 |
|---|---|
| 1 | 1 = 对所有能带求和计算 Berry 曲率；0 = 对孤立能带组计算（只在该组能带与其它能带完全分离时有效） |
| 2 | 计入的能带数（从最低能带数起）；0 = 全部价带 |
| 3 | 计算所在的 k 平面，例如 `0 0 1` 表示 k3 = 0 平面 |
| 4 | 该平面的网格划分，如 100×100 |

输出 `BerryC*.dat`（* 为计入的能带数）：文件头含网格、能带数与 `#ChernNumber`，其后三列为 k1、k2 与 Berry 曲率 F（Å²）；`splot "BerryC8.dat" w l` 可画三维图。官方示例 `Graphene-Chern.dat`（石墨烯，`HS.fileout on`）展示了 K、K' 附近的大 Berry 曲率。

### 8.2 物理检查

- 时间反演对称体系（如无磁性的石墨烯、Bi2Se3）所有占据带的总 Chern 数必为 0——K 与 K' 的 Berry 曲率等大反号。拿到非零值时首先怀疑网格太粗或带隙闭合；
- 非零 Chern 数需要破坏时间反演（磁性 + SOC），例如量子反常霍尔体系；此时用自旋极化的非共线 + SOC 计算；
- 只有当所选能带与其它能带全程有能隙时，Chern 数才是良定义的整数；金属中 Berry 曲率的积分对应的是反常霍尔电导（§9），而不是整数。

---

## 9. 反常霍尔电导与 Wannier 电荷中心

OpenMX 本身之外，有几类基于 scfout 的社区工具（均为外部程序，按各自说明编译使用）：

| 工具 | 输入 | 用途 | 说明 |
|---|---|---|---|
| [OMXsigmaxy](https://github.com/hikaruri/OMXsigmaxy) | scfout（`HS.fileout on`） | 基于扩展 FHS 方法计算反常霍尔电导 σxy | 与 OpenMX 使用相同的 CC/LIB 编译；适合铁磁金属（如 bcc Fe）的本征 AHC |
| [HWC](https://github.com/Ncmexp2717/HWC) | scfout | 混合 Wannier 电荷中心（WCC 演化）、层极化 | 复制进 OpenMX source 后 `make HWC`；WCC 随 k 的"缠绕"可判断 Z2 与 Chern 数 |
| [pymx](https://github.com/knhr0727/pymx) | scfout | Python 构建 H(k)，示例含 Bi2Se3 能带、宇称本征值与 WCC（Wilson loop） | 适合理解宇称判据（Fu–Kane）的实现思路 |
| [HopTB.jl](https://github.com/HopTB/HopTB.jl) | OpenMX 哈密顿量 | 紧束缚模型上的 Berry 曲率、AHC、光学与非线性响应 | Julia |
| Wannier90 / WannierTools | OpenMX 的 Wannier90 格式输出 | 表面态、Wilson loop、Weyl 点 | 见 第14章 的 MLWF/CWF 接口 |

> [!TIP]
> 对有反演中心的体系（如 Bi2Se3），Fu–Kane 宇称判据给出 $(-1)^{\nu_0}=\prod_{i=1}^{8}\prod_{m}\xi_{2m}(\Gamma_i)$（8 个 TRIM 点上占据 Kramers 对的宇称乘积），是 Z2FH 结果的独立检验。pymx 的示例实现了这一计算，可作为思路参考（不要复制其代码，自己实现时只需 scfout 中的哈密顿量与重叠矩阵）。

---

## 10. Rashba 体系与费米面

### 10.1 Rashba 参数

表面、界面或极性体相（如 BiTeI）中，Rashba 型 SOC 使能带变为 $E_\pm(k)=\frac{\hbar^2k^2}{2m^*}\pm\alpha_R|k|$。从 SOC 能带读出两个量即可：动量偏移 $k_0$（劈裂能带极值相对 Γ 的偏移）与 Rashba 能 $E_R$（极值相对交叉点的深度），则

$$
\alpha_R=\frac{2E_R}{k_0}.
$$

在 `.Band` 中沿 Γ 附近加密 k 点（例如在 ±0.1 倒格矢范围内取 100 个点），用 `soc_splitting.py` 或自写脚本读出极值位置。板模型中上下两个表面的 Rashba 态会简并（官方 `Au111Surface23_FL.dat` 的例子中出现两组 Rashba 带），可用 kSpin + MulPCalc 只取一侧表面原子的贡献，或用 ESM/偶极修正区分（[第16章](16-surface-esm-workfunction.md)）。

### 10.2 费米面文件（bxsf / frmsf）

在 DOS 计算中加入 `FermiSurfer.fileout`：

```text
Dos.fileout                 on
Dos.Erange                  -20.0  20.0
Dos.Kgrid                   61 61 61      # smooth Fermi surfaces need many k points
FermiSurfer.fileout         on            # default off
```

| 计算类型 | XCrySDen 用 | FermiSurfer 用 |
|---|---|---|
| 共线自旋极化 | `System.Name.FermiSurf0.bxsf`（上）、`FermiSurf1.bxsf`（下） | `System.Name.FermiSurf_s0_aA.frmsf`、`_s1_aA.frmsf`（A 为原子序号，用于按原子着色） |
| 非共线（含 SOC） | `System.Name.FermiSurf.bxsf` | `System.Name.FermiSurf_aA.frmsf` |

> [!WARNING]
> **patch 4.0.1 改变了 bxsf/frmsf 的单位**：之前版本为 Hartree 与 Bohr，打补丁后按 XCrySDen 的约定改为 eV 与 Å。若你的后处理脚本按旧单位换算，费米面能量与 k 尺度会差一个常数因子；混用不同版本的结果时务必核对。使用 `DosGauss.fileout`（而非 `Dos.fileout`）时不会生成 FermiSurfer 文件。

OpenMX Viewer（2026-05 更新）可以直接读 BXSF/FRMSF，显示第一布里渊区、特殊 k 点、费米能平移、截面，还能计算闭合轨道面积与对应的 de Haas–van Alphen 频率。官方示例 `Cafcc_FS.dat`（fcc Ca）可用来练习。

---

## 参数选择与收敛建议

| 参数 | 推荐 | 说明 |
|---|---|---|
| `scf.SpinPolarization` | `nc` | SOC 的前提 |
| 赝势 / 基组 | DFT_DATA19（j 依赖）+ Standard；重元素与 Precise 对比 | 极化函数影响 SOC 劈裂 |
| `scf.criterion` | 1e-9（劈裂、拓扑）；1e-8 一般能带 | 能量尺度小，标准要严 |
| `scf.Kgrid` | 先用无 SOC 计算收敛，再沿用到 SOC | SOC 很少改变 k 点收敛性 |
| 混合 | 非磁 SOC：RMM-DIISK，初始权重 0.1；磁性 + SOC：更小权重（0.01–0.05） | 磁性 + SOC 收敛难，参见 [第11章](11-spin-noncollinear.md) |
| 对照计算 | 同一 VPS、同一 nc 框架，只切换 `scf.SpinOrbit.Coupling` | 唯一变量是 SOC |
| `scf.SO.factor` | 每个 species 写满 s–f | 用于拓扑相变扫描、贡献分解 |
| Z2FH 网格 | 第 1 行 5 起步，加到 8–10 验证 | 6 个面满足一致性条件 |
| calB 网格 | 100×100 起步 | Chern 数为整数且随网格不变 |
| 费米面 `Dos.Kgrid` | 50³–80³ | 光滑费米面需要很多 k 点 |
| kSpin MPI | 进程数 ≤ 对应模式的 k 点数 | FermiLoop 受 `k-plane.2ndStep` 限制 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `Spin-orbit coupling is not supported for collinear DFT calculations.` | 共线（`on`）+ SOC | `scf.SpinPolarization nc` |
| `scf.SpinOrbit.Coupling must be OFF when scf.SpinPolarization=OFF` | 非自旋 + SOC | 同上 |
| SOC 打开了但能带看不到劈裂 | 未生效（拼写错误被当作未知关键词忽略）、或该处本应简并（Kramers、反演+时间反演） | 在 `.out` 的输入回显中确认；对照 §2.2 的对称性判断 |
| 劈裂数值与文献差很多 | 基组太小、晶格常数不同、泛函不同 | 用 Standard/Precise 对比；核对结构；比较时注意 LDA/GGA 差异 |
| `.Band` 中价带顶简并度"不对" | 简并分组容差不合适，或 k 点不在高对称点上 | 调整 `soc_splitting.py --tol`；检查 `Band.kpath` 端点坐标 |
| 能带标签错位（fcc/六方/菱方） | `Band.kpath` 坐标基与原胞基矢顺序不匹配 | 重新推导，或用 `<Band.KPath.UnitCell`、OpenMX Viewer 自动生成路径 |
| Z2.dat 中 x0+xπ ≠ y0+yπ | 网格太粗，或体系不是绝缘体 | 增大 `Z2FH.in` 第 1 行；检查带隙 |
| Z2FH 用于磁性体系结果无意义 | 时间反演被破坏 | 改用 calB 计算 Chern 数 |
| 拓扑结论与文献相反 | 结构（内坐标、晶格常数）或泛函不同 | 用实验或可靠优化结构；对 SOC 强度做扫描确认相变位置 |
| calB 得到非整数 Chern 数 | 所选能带与其它能带有交叠（金属） | 只对孤立能带组计算；金属请计算 AHC |
| kSpin 报读不到 scfout | `HS.fileout` 未开或 `Filename.scfout` 写错 | 确认 SCF 时 `HS.fileout on`，文件名含 `.scfout` |
| kSpin 输出为空 | `Energy.Range` 内没有能带 | 能量相对化学势（eV），适当放宽窗口 |
| nc + SOC + CDDF 光学计算段错误，或 `Parameter 10 was incorrect on entry to ZHEEVX` | 4.0 已知问题 | 参考 [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix) 的分析，或关注官方补丁 |
| 费米面能量/尺度不对 | patch 4.0.1 前后 bxsf/frmsf 单位不同 | 打补丁后为 eV/Å，按版本处理 |
| 自旋螺旋 + SOC 报错 | 广义 Bloch 定理与 SOC 不兼容 | 二者不能同时使用（见 [第11章](11-spin-noncollinear.md) §10） |
| SOC 计算内存不足 | 维度加倍 | 减少每节点 MPI 进程、改用 MPI×OpenMP、或减小基组（见 第18章） |

---

## 练习

1. **GaAs 基组收敛**：分别用 Quick（Ga7.0-s2p2d1 / As7.0-s3p2d1）、Standard、Precise 基组计算 Γ 点 Δ0，与实验 0.34 eV 比较。
2. **放大的 SOC**：在 GaAs 中用 `scf.SO.factor` 把 As 的 p 因子设为 2.0（其余 1.0），Δ0 大约变为多少？再只放大 Ga 的 p，比较两者对 Δ0 的贡献。
3. **MoS2 → WS2**：把 Mo 换成 W（W7.0-s3p2d2f1 / W_PBE19，12 个价电子，注意初始占据与晶格常数），比较 K 点价带劈裂的变化。
4. **Bi2Se3 拓扑相变**：按 §5.3 扫描 Bi 的 p 因子，画出 Γ 点能隙随因子的变化，并在相变两侧各算一次 Z2。
5. **石墨烯 Berry 曲率**：用官方 `Graphene-Chern.dat` 运行 calB，确认总 Chern 数为 0，并画出 K、K' 附近反号的 Berry 曲率。
6. **费米面**：对 bcc Fe（[第11章](11-spin-noncollinear.md) 的输入，改为 nc + SOC）输出 `FermiSurf.bxsf`，在 OpenMX Viewer 或 XCrySDen 中查看。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 37 章 Non-collinear DFT：<https://www.openmx-square.org/openmx_man4.0/s37_ncl-dft.html>
- 第 38 章 Relativistic effects（38.1 全相对论、38.2 SOC 强度控制、38.3 标量相对论）：<https://www.openmx-square.org/openmx_man4.0/s38_rel-effects.html>
- 第 56 章 动量空间自旋织构（kSpin）：<https://www.openmx-square.org/openmx_man4.0/s56_spin-texture-k.html>
- 第 58 章 Chern 数与 Berry 曲率（calB，FHS 方法）：<https://www.openmx-square.org/openmx_man4.0/s58_chern-berry-fhs.html>
- 第 59 章 Z2 不变量（Z2FH）：<https://www.openmx-square.org/openmx_man4.0/s59_z2-fh.html>
- 第 65 章 Fermi surfaces：<https://www.openmx-square.org/openmx_man4.0/s65_fermi-surface.html>
- 第 55.5 节 非共线/SOC 下的能带展开：<https://www.openmx-square.org/openmx_man4.0/s55_5_ncl.html>
- 第 22 章 能带结构：<https://www.openmx-square.org/openmx_man4.0/s22_band.html>
- Patch 4.0.1 说明（bxsf/frmsf 单位、Mulliken HOMO 修复、GaAs.dat）：<https://www.openmx-square.org/bugfixed/26May08/README.txt>
- DFT_DATA19 数据库：<https://www.openmx-square.org/vps_pao2019/>
- kSpin 技术说明：<http://www.openmx-square.org/tech_notes/note_kSpin-1_0.pdf>
- OpenMX Viewer：<https://www.openmx-square.org/viewer/>
- FermiSurfer：<https://mitsuaki1987.github.io/fermisurfer/>

**相关 GitHub 仓库（仅链接）**

- [hikaruri/OMXsigmaxy](https://github.com/hikaruri/OMXsigmaxy)：基于扩展 FHS 方法从 scfout 计算反常霍尔电导。
- [hikaruri/Z2-support-OpenMX](https://github.com/hikaruri/Z2-support-OpenMX)：OpenMX 自带 Z2FH 的使用支持页。
- [Ncmexp2717/HWC](https://github.com/Ncmexp2717/HWC)：读 scfout 计算混合 Wannier 中心与层极化。
- [knhr0727/pymx](https://github.com/knhr0727/pymx)：Python 读 scfout 构建 H(k)，示例含 Bi2Se3 宇称本征值与 Wilson loop。
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：基于 OpenMX 哈密顿量的 Berry 曲率、光学与非线性响应。
- [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix)：记录 4.0 中 nc+SOC+CDDF 的崩溃问题与补丁。
- [Ncmexp2717/OMXTool](https://github.com/Ncmexp2717/OMXTool)：预览 .Band、给能带着色（手册 kSpin 一节提到用于标出 Rashba 带）。
- [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS)：浏览器端能带 + PDOS 查看器。
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，含 Bi2Se3 SOC 示例。
- [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-)：印尼语教程，含 SOC、Rashba、Berry 曲率、反常霍尔、边缘态、Z2 等分支。

**关键文献**

- G. Theurich, N. A. Hill, Phys. Rev. B 64, 073106 (2001)（非共线框架下的 j 依赖赝势与 SOC）. DOI: [10.1103/PhysRevB.64.073106](https://doi.org/10.1103/PhysRevB.64.073106)
- T. Fukui, Y. Hatsugai, H. Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005)（FHS 方法）. DOI: [10.1143/JPSJ.74.1674](https://doi.org/10.1143/JPSJ.74.1674)
- T. Fukui, Y. Hatsugai, J. Phys. Soc. Jpn. 76, 053702 (2007)（Z2 的格点计算）. DOI: [10.1143/JPSJ.76.053702](https://doi.org/10.1143/JPSJ.76.053702)
- H. Sawahata, N. Yamaguchi, H. Kotaka, F. Ishii, Jpn. J. Appl. Phys. 57, 030309 (2018)（Z2FH）. DOI: [10.7567/JJAP.57.030309](https://doi.org/10.7567/JJAP.57.030309)
- H. Kotaka, F. Ishii, M. Saito, Jpn. J. Appl. Phys. 52, 035204 (2013)（kSpin）. DOI: [10.7567/JJAP.52.035204](https://doi.org/10.7567/JJAP.52.035204)
- L. Fu, C. L. Kane, Phys. Rev. B 76, 045302 (2007)（宇称判据）. DOI: [10.1103/PhysRevB.76.045302](https://doi.org/10.1103/PhysRevB.76.045302)
- H. Zhang et al., Nat. Phys. 5, 438 (2009)（Bi2Se3 家族的拓扑绝缘体预言）. DOI: [10.1038/nphys1270](https://doi.org/10.1038/nphys1270)
- Y. Xia et al., Nat. Phys. 5, 398 (2009)（Bi2Se3 的 ARPES 观测）. DOI: [10.1038/nphys1274](https://doi.org/10.1038/nphys1274)
- Z. Y. Zhu, Y. C. Cheng, U. Schwingenschlögl, Phys. Rev. B 84, 153402 (2011)（单层 TMD 的 SOC 自旋劈裂）. DOI: [10.1103/PhysRevB.84.153402](https://doi.org/10.1103/PhysRevB.84.153402)
- D. Xiao, G.-B. Liu, W. Feng, X. Xu, W. Yao, Phys. Rev. Lett. 108, 196802 (2012)（自旋–能谷耦合）. DOI: [10.1103/PhysRevLett.108.196802](https://doi.org/10.1103/PhysRevLett.108.196802)
- I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001)（III–V 半导体能带参数，含 GaAs 的 Δso）. DOI: [10.1063/1.1368156](https://doi.org/10.1063/1.1368156)
- D. Xiao, M.-C. Chang, Q. Niu, Rev. Mod. Phys. 82, 1959 (2010)（Berry 相综述）. DOI: [10.1103/RevModPhys.82.1959](https://doi.org/10.1103/RevModPhys.82.1959)
- M. Z. Hasan, C. L. Kane, Rev. Mod. Phys. 82, 3045 (2010)（拓扑绝缘体综述）. DOI: [10.1103/RevModPhys.82.3045](https://doi.org/10.1103/RevModPhys.82.3045)
- M. Kawamura, Comput. Phys. Commun. 239, 197 (2019)（FermiSurfer）. DOI: [10.1016/j.cpc.2019.01.017](https://doi.org/10.1016/j.cpc.2019.01.017)

[← 上一章](11-spin-noncollinear.md) | [返回目录](../README.md) | [下一章 →](13-dft-u.md)
