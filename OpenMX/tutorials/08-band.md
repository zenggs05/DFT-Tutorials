# 第8章 能带结构（含自动 k 路径与能带展开）（Band Structure, Automatic k-Path & Band Unfolding with OpenMX）

本章讲解 OpenMX 4.0 能带结构计算的完整流程：`Band.dispersion`、`Band.Nkpath`、`<Band.kpath` 的写法，fcc/bcc/六方晶格为什么要用 `<Band.KPath.UnitCell`，如何用 OpenMX Viewer、cif2omx 或 seekpath 自动生成高对称 k 路径，用 `bandgnu13` 与 Python 解析 `.Band` 文件、提取带隙与 VBM/CBM，处理自旋极化与非共线能带，以及用 4.0 的 `scf.restart pp` 只做后处理。最后介绍超胞能带展开（band unfolding：`Unfolding.Electronic.Band`、`<Unfolding.ReferenceVectors`、`<Unfolding.Map`、`<Unfolding.kpoint`），把缺陷/掺杂超胞的能带"展开"回原胞布里渊区，便于与 ARPES 对比。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出正确的 `<Band.kpath`，理解"分点数 + 起点 + 终点 + 两个标签"每一列的含义；
- 知道 k 坐标以哪组倒格矢为单位，在 fcc、bcc、六方、菱方体系中正确使用 `<Band.KPath.UnitCell`；
- 用 OpenMX Viewer / cif2omx 自动生成 Setyawan–Curtarolo 路径，并能检查它是否对应你的晶胞；
- 用 `bandgnu13` 生成 gnuplot 数据，或用 Python 直接读 `.Band`（Hartree）得到以化学势为零点的 eV 能带、带隙、VBM/CBM 位置；
- 区分非自旋、共线自旋（`BANDDAT1/2`）与非共线/SOC（单文件、本征值加倍）三种输出；
- 用 `scf.restart pp`（OpenMX ≥ 4.0）在不重复 SCF、不覆盖重启文件的前提下补算能带；
- 对超胞做能带展开，读懂 `unfold_totup`、`unfold_orbup` 并作图。

| 任务 | 本章位置 | 示例 |
|---|---|---|
| 半导体能带、带隙 | §4、§7 | [算例 8-1 金刚石与 Si](../examples/08-band-cdia-si/README.md) |
| fcc/bcc 高对称点标签 | §5 | 同上 + Fe |
| 自动 k 路径 | §6 | — |
| 二维材料、自旋极化金属 | §8 | [算例 8-2 石墨烯与 bcc Fe](../examples/08-band-graphene-fe/README.md) |
| 缺陷/掺杂超胞的能带展开 | §10 | [算例 8-3 Si 超胞展开](../examples/08-unfold-si-supercell/README.md) |

> [!NOTE]
> 能带计算的前提是 SCF 已经收敛且参数合理：基组与赝势见[第4章](04-pao-vps.md)，SCF 与 k 点收敛见[第5章](05-scf-kpoints-convergence.md)，晶格常数最好先按[第6章](06-optimization-eos-neb.md)优化。态密度在[第9章](09-dos-pdos-cohp-optics.md)，SOC 与拓扑能带在[第12章](12-soc-topology.md)。

---

## 2. 原理速览

### 2.1 LCAO 下的能带

OpenMX 用数值赝原子轨道（PAO）$\phi_{i\alpha}$ 作基组。对每个 $\mathbf{k}$，Bloch 和
$\chi_{i\alpha}^{\mathbf{k}}=\frac{1}{\sqrt{N}}\sum_{\mathbf{R}}e^{i\mathbf{k}\cdot\mathbf{R}}\phi_{i\alpha}(\mathbf{r}-\boldsymbol{\tau}_i-\mathbf{R})$
把 Kohn–Sham 方程化为广义本征值问题

$$
H(\mathbf{k})\,c_\mu(\mathbf{k})=\varepsilon_\mu(\mathbf{k})\,S(\mathbf{k})\,c_\mu(\mathbf{k})
$$

SCF 收敛后，沿高对称路径上的一系列 $\mathbf{k}$ 点对角化 $H(\mathbf{k})$ 即得能带。两个直接后果：

1. **能带数等于基函数总数**。例如 `C6.0-s2p2d1` 每个原子 2×1+2×3+1×5 = 13 个基函数，金刚石原胞 2 个原子就是 26 条能带。
2. **价带和低能导带可靠，高能导带依赖基组**。远高于费米能的导带受基组完备性限制；关心导带时应比较 Standard 与 Precise 基组，必要时在空隙处放置空原子（官方 SiC 展开示例就在六角形中心放了带长尾基函数的空原子 `E` 来改善导带）。

### 2.2 倒格矢与 k 路径

倒格矢 $\mathbf{b}_i\cdot\mathbf{a}_j=2\pi\delta_{ij}$。`<Band.kpath` 中的 $(k_1,k_2,k_3)$ 表示 $\mathbf{k}=k_1\mathbf{b}_1+k_2\mathbf{b}_2+k_3\mathbf{b}_3$。**同一个物理 k 点，在原胞倒格矢下和惯用胞倒格矢下的坐标不同**——这是能带标签错位最常见的原因（§5）。

### 2.3 带隙

- 直接带隙：VBM 与 CBM 在同一 k 点；间接带隙：不在同一 k 点。
- 能带图只沿路径采样，严格的 VBM/CBM 可能不在路径上（如低对称性体系），需要配合态密度或全布里渊区网格确认。
- GGA（PBE）系统性低估带隙；OpenMX 4.0 的 `scf.XcType` 只有 `LDA`、`LSDA-CA`、`LSDA-PW`、`GGA-PBE`，不支持杂化泛函，需要更准确的带隙时可考虑 DFT+U（[第13章](13-dft-u.md)）或其他程序。

---

## 3. 能带计算的流程与关键词

OpenMX 的能带在**同一次运行**中、SCF 收敛之后自动计算，然后用后处理程序转换格式：

```text
.dat (Band.dispersion on) --openmx--> System.Name.Band --bandgnu13--> .GNUBAND + .BANDDAT1(/2) --gnuplot/Python--> 图
```

| 关键词 | 取值 / 默认 | 说明 |
|---|---|---|
| `Band.dispersion` | `on` / `off`（默认 `off`） | 打开能带计算；只在 `scf.EigenvalueSolver Band` 时生效 |
| `Band.Nkpath` | 整数（默认 0） | 路径段数，必须等于 `<Band.kpath` 的行数 |
| `<Band.kpath` | 每行：`分点数 k1 k2 k3 k1' k2' k3' 标签1 标签2` | 分点数包含两个端点；标签是不含空格的字符串，Γ 习惯写 `g` 或 `G` |
| `<Band.KPath.UnitCell` | 3 行晶格矢量，单位同 `Atoms.UnitVectors.Unit` | 若给出，k 路径坐标按这组晶格矢量的倒格矢解释；不给则按 `Atoms.UnitVectors` |
| `scf.eigen.lib` | `elpa1` / `elpa2`（4.0 中 Band 计算默认 `elpa2`） | 对角化库，一般不用改 |

> [!TIP]
> 分点数决定能带曲线的平滑程度。建议让各段的点数大致正比于段长，例如每 0.1 × (2π/a) 取 3 个点；需要求有效质量时在极值点附近单独加密一段。相邻两段共享的端点会被计算两次，这对结果没有影响。

---

## 4. 算例一：金刚石的能带

完整文件见 [../examples/08-band-cdia-si/README.md](../examples/08-band-cdia-si/README.md)。

### 4.1 完整输入 `Cdia_band.dat`

```text
# Diamond (C, fcc primitive cell): SCF + band structure along the
# Setyawan-Curtarolo FCC path G-X-W-K-G-L-U-W-L-K|U-X. OpenMX 4.0.
System.CurrrentDirectory    ./
System.Name                 cdia
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: dir containing PAO/ and VPS/
level.of.stdout             1
level.of.fileout            1

Species.Number              1
<Definition.of.Atomic.Species
  C   C6.0-s2p2d1   C_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  C   0.00  0.00  0.00   2.0 2.0
  2  C   0.25  0.25  0.25   2.0 2.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  1.7835  1.7835  0.0000
  1.7835  0.0000  1.7835
  0.0000  1.7835  1.7835
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            220.0
scf.maxIter                 100
scf.EigenvalueSolver        Band
scf.Kgrid                   10 10 10
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.criterion               1.0e-8

MD.Type                     Nomd

Band.dispersion             on
<Band.KPath.UnitCell
  3.567  0.000  0.000
  0.000  3.567  0.000
  0.000  0.000  3.567
Band.KPath.UnitCell>
Band.Nkpath                 10
<Band.kpath
  30   0.00 0.00 0.00    1.00 0.00 0.00    g  X
  15   1.00 0.00 0.00    1.00 0.50 0.00    X  W
  11   1.00 0.50 0.00    0.75 0.75 0.00    W  K
  32   0.75 0.75 0.00    0.00 0.00 0.00    K  g
  26   0.00 0.00 0.00    0.50 0.50 0.50    g  L
  18   0.50 0.50 0.50    1.00 0.25 0.25    L  U
  11   1.00 0.25 0.25    1.00 0.50 0.00    U  W
  21   1.00 0.50 0.00    0.50 0.50 0.50    W  L
  18   0.50 0.50 0.50    0.75 0.75 0.00    L  K
  11   1.00 0.25 0.25    1.00 0.00 0.00    U  X
Band.kpath>
```

### 4.2 关键行解释

| 内容 | 说明 |
|---|---|
| `Atoms.UnitVectors` | fcc 原胞，a = 3.567 Å（实验值），三个矢量为 (a/2)(1,1,0) 等 |
| `<Band.KPath.UnitCell` | 惯用立方胞（a = 3.567 Å）。于是 k 坐标以 $2\pi/a$ 为单位：X = (1,0,0)，L = (½,½,½)，W = (1,½,0)，K = (¾,¾,0)，U = (1,¼,¼) |
| 第 1 段 `30 0 0 0 1 0 0 g X` | Γ→X 取 30 个点（含两端），标签 `g`（Γ）和 `X` |
| 第 9、10 段 | L→K 之后从 U 重新开始到 X，这是路径中的"断点" K\|U。`bandgnu13` 把两段在横轴上首尾相接（零长度跳跃），作图时两个标签画在同一位置 |
| `scf.Kgrid 10 10 10` | SCF 的布里渊区网格，与能带路径无关；先按第 5 章做收敛测试 |

### 4.3 运行与 `bandgnu13`

```bash
cd OpenMX/examples/08-band-cdia-si
mpirun -np 8 openmx Cdia_band.dat -nt 1 | tee Cdia_band.std
ls cdia.Band                                   # written after the SCF converged
```

`bandgnu13` 不随 `make install` 编译，需要单独编译一次：

```bash
cd /path/to/openmx4.0/source
gcc bandgnu13.c -lm -o bandgnu13
cp bandgnu13 ~/bin/                 # or copy it into the working directory
cd -  && bandgnu13 cdia.Band        # -> cdia.GNUBAND, cdia.BANDDAT1
gnuplot -p cdia.GNUBAND
```

生成的文件：

| 文件 | 内容 |
|---|---|
| `cdia.GNUBAND` | gnuplot 脚本：坐标轴、高对称点刻度、调用数据文件 |
| `cdia.BANDDAT1` | （上自旋或非自旋）能带数据：两列"路径长度（Bohr⁻¹）  E − μ（eV）" |
| `cdia.BANDDAT2` | 仅共线自旋极化计算：下自旋 |

> [!IMPORTANT]
> 单位约定：`System.Name.Band` 中的本征值和化学势都是 **Hartree、未平移**；`bandgnu13` 输出的 `BANDDAT` 已换算为 **eV，并以化学势为零点**（换算因子 27.2113845）。横轴是沿路径累计的 k 空间距离，单位 Bohr⁻¹（含 2π）。

`BANDDAT1` 的结构（读源码 `bandgnu13.c` 可知）：文件开头先是若干个"竖线"数据块（每个路径分段处一条、横坐标相同的竖线，供 gnuplot 画分隔线），之后才是能带——每条能带按路径段输出，块与块之间用两个空行分隔。用 Python 读 `BANDDAT` 时要先跳过横坐标不变的竖线块（[算例 8-3 的 `plot_unfold.py`](../examples/08-unfold-si-supercell/plot_unfold.py) 里有现成实现）；更推荐直接读 `.Band`（§7）。

修改 gnuplot 外观时，编辑 `.GNUBAND` 中的 `set yrange`、`set xtics` 等即可；也可以用 [Ncmexp2717/OMXTool](https://github.com/Ncmexp2717/OMXTool) 这类 GUI 工具预览 `.Band`。

### 4.4 结果判读

（示意输出，数值仅供参考）

```text
VBM =  -x.xxxx eV at G k=(0.000 0.000 0.000) (spin 0, band 4)
CBM =   x.xxxx eV at k=(0.xxx 0.xxx 0.xxx) (spin 0, band 5)
band gap along the path = 4.1xxx eV (indirect)
```

金刚石是间接带隙：VBM 在 Γ（三重简并），CBM 在 Γ–X 线上靠近 X 处。PBE 带隙文献值约 4.1–4.2 eV，实验值 5.48 eV。8 个价电子占据 4 条能带，所以 VBM 是第 4 条带、CBM 是第 5 条带。

---

## 5. `<Band.KPath.UnitCell` 与高对称点

### 5.1 为什么 fcc/bcc 一定要注意

OpenMX 手册（第 22 章）明确建议：对体心、面心、底心和菱方晶格，按能带计算的惯例用 `<Band.KPath.UnitCell` 定义倒格矢。做个对比实验：把 `Cdia_band.dat` 中的 `<Band.KPath.UnitCell` 块删掉而不改坐标，`1.00 0.00 0.00` 就变成了"原胞倒格矢 $\mathbf{b}_1$"，它不是 X 点（实际上是另一个倒格点方向上的点），图上标着 X 的位置画的是别的地方的能带。

| 情况 | k 坐标的单位 |
|---|---|
| 写了 `<Band.KPath.UnitCell` | 该块中晶格矢量的倒格矢 |
| 没写 | `Atoms.UnitVectors`（SCF 晶胞）的倒格矢 |

`.Band` 文件中写出的路径坐标已经换算为 SCF 晶胞的倒格矢单位（源码中 `kpath_changeunit` 完成换算），所以在 `.Band` 里看到的数字与输入不同是正常的。

### 5.2 常用高对称点（供手写路径时查对）

**fcc**（`<Band.KPath.UnitCell` 为惯用立方胞，单位 $2\pi/a$）：

| 点 | 坐标 |
|---|---|
| Γ | (0, 0, 0) |
| X | (1, 0, 0) |
| L | (½, ½, ½) |
| W | (1, ½, 0) |
| K | (¾, ¾, 0) |
| U | (1, ¼, ¼) |

Setyawan–Curtarolo 推荐路径：Γ–X–W–K–Γ–L–U–W–L–K\|U–X。

**bcc**（惯用立方胞，单位 $2\pi/a$）：

| 点 | 坐标 |
|---|---|
| Γ | (0, 0, 0) |
| H | (1, 0, 0) |
| N | (½, ½, 0) |
| P | (½, ½, ½) |

推荐路径：Γ–H–N–Γ–P–H\|P–N。

**六方 / 二维六方**（原胞，$\mathbf{a}_2=(-a/2,\sqrt3a/2,0)$，即 120° 约定）：

| 点 | 坐标 |
|---|---|
| Γ | (0, 0, 0) |
| M | (½, 0, 0) |
| K | (⅓, ⅓, 0) |
| A | (0, 0, ½) |
| L | (½, 0, ½) |
| H | (⅓, ⅓, ½) |

> [!WARNING]
> 若六方晶格矢量取 60° 约定 $\mathbf{a}_2=(a/2,\sqrt3a/2,0)$，K 点坐标变为 (⅔, ⅓, 0)（或等价点），(⅓, ⅓, 0) 就不再是 K。每次都要根据自己的晶格矢量核对，或用自动工具生成（§6）。

**菱方**（如 Bi₂Se₃）以及单斜、三斜等低对称体系，高对称点依赖晶格参数，强烈建议用 §6 的工具生成。

### 5.3 算例：Si 的能带

`Si_band.dat` 与金刚石结构相同，只把 a 换成 5.43 Å，路径取 L–Γ–X–W–K–Γ（与 §10 的能带展开算例相同，便于叠加对比）：

```text
Band.dispersion             on
<Band.KPath.UnitCell
  5.430  0.000  0.000
  0.000  5.430  0.000
  0.000  0.000  5.430
Band.KPath.UnitCell>
Band.Nkpath                 5
<Band.kpath
  26   0.50 0.50 0.50    0.00 0.00 0.00    L  g
  30   0.00 0.00 0.00    1.00 0.00 0.00    g  X
  15   1.00 0.00 0.00    1.00 0.50 0.00    X  W
  11   1.00 0.50 0.00    0.75 0.75 0.00    W  K
  32   0.75 0.75 0.00    0.00 0.00 0.00    K  g
Band.kpath>
```

这里原胞取 $\mathbf{a}_1=(0,a/2,a/2)$、$\mathbf{a}_2=(a/2,0,a/2)$、$\mathbf{a}_3=(a/2,a/2,0)$，原子在分数坐标 (0,0,0) 与 (¼,¼,¼)。PBE 下 Si 是约 0.5–0.6 eV 的间接带隙（实验约 1.17 eV），CBM 在 Γ–X 线上约 85% 处。

---

## 6. 自动生成 k 路径

### 6.1 OpenMX Viewer（推荐）

OpenMX Viewer 是官方的浏览器端工具（2026-05 大幅更新，也有离线 Windows 版）：https://www.openmx-square.org/viewer/

1. 拖入 CIF / POSCAR / OpenMX `.dat` 等结构文件（文件在本地浏览器中处理，不上传服务器），也可以从 Materials Project / COD 直接获取结构；
2. 打开 "Brillouin zone" 显示，查看第一布里渊区和红色的推荐路径（基于 Setyawan–Curtarolo 方案）；
3. 导出 OpenMX 输入文件时，`Band.Nkpath` 和 `<Band.kpath` 会按当前显示的推荐路径自动写入。

手册（第 22.1 节、10.2 节）提醒：推荐路径与所选晶胞有关，最终仍要检查 `<Band.KPath.UnitCell` 和 `<Band.kpath` 是否对应你想要的晶胞；需要与文献完全一致的路径、或要在某个交叉点/极值点附近加密时，在自动路径基础上手动修改。

### 6.2 cif2omx（命令行，适合批量）

cif2omx 是 4.0 起正式文档化的 CIF→`.dat` 工具，源码在 `openmx4.0/source`：

```bash
cd /path/to/openmx4.0/source
gcc -std=c99 -O2 -Wall -Wextra -pedantic -DCIF2OPENMX_MAIN cif2omx.c -lm -o cif2omx
./cif2omx Si_mp149.cif -o Si.dat -b standard -k 0.25 -e 220
./cif2omx --help
```

其中 `-b quick|standard|precise` 选择基组档次，`-k` 按倒空间间距（1/Å，含 2π）设置 `scf.Kgrid`（默认 0.33，与 OpenMX Viewer 一致），`-e` 写入 `scf.energycutoff`（不给时写 220.0），`-c xyz|fractional` 选择坐标格式。手册说明 cif2omx 生成 `.dat` 时也可以自动生成 k 路径，可作为批量能带计算的起点；生成后用 OpenMX Viewer 打开检查结构与路径。

### 6.3 seekpath 与仓库脚本

也可以用 Python 的 seekpath（HPKOT 约定）或 pymatgen 的 `HighSymmKpath` 生成路径。仓库脚本 [omx_make_input.py](../scripts/omx_make_input.py) 可从 POSCAR/CIF 生成 `.dat` 并可选写入 `<Band.kpath`，用 `python ../scripts/omx_make_input.py -h` 查看用法。

> [!CAUTION]
> seekpath 会先把晶胞**标准化**（可能旋转、换轴或换成另一种原胞），它给出的 k 坐标对应标准化后的晶胞。要么用 seekpath 返回的晶胞做 SCF，要么把 k 坐标变换回你自己的晶胞；不同工具（Setyawan–Curtarolo、HPKOT）对同一晶系的标签和路径也可能不同，与文献对比时要说明所用约定。

---

## 7. `.Band` 文件格式与 Python 解析

### 7.1 文件格式

`System.Name.Band` 是纯文本（由 `Band_DFT_kpath.c` 写出、`bandgnu13.c` 读取）：

```text
n_eig  spin_flag  mu(Hartree)                  <- line 1; spin_flag 0: non-spin or nc, 1: collinear
b11 b12 b13 b21 b22 b23 b31 b32 b33            <- reciprocal vectors of the SCF cell (Bohr^-1)
Nkpath                                          <- number of segments
n  k1 k2 k3  k1' k2' k3'  lab1 lab2             <- Nkpath lines (k already in SCF-cell units)
n  k1 k2 k3                                     <- for each segment, each k point, each spin:
e1 e2 ... en                                    <-   n eigenvalues in Hartree
...
```

共线自旋极化时，每个 k 点先写上自旋块、再写下自旋块；非共线（`nc`，含 SOC）时 `spin_flag` 为 0，但本征值个数是 2 × 基函数数。

### 7.2 解析与带隙

示例目录中的 [band_parse.py](../examples/08-band-cdia-si/band_parse.py) 完整实现了解析、带隙判定和作图，核心部分如下：

```python
import numpy as np
HA2EV = 27.2113845

tok = open("si.Band").read().split()
neig, spin_flag, mu = int(tok[0]), int(tok[1]), float(tok[2])
rtv = np.array(tok[3:12], float).reshape(3, 3)          # rows: b1, b2, b3 (Bohr^-1)
nseg = int(tok[12]); p = 13
segs = []
for _ in range(nseg):
    segs.append((int(tok[p]), tok[p + 7], tok[p + 8])); p += 9
nspin = 2 if spin_flag == 1 else 1
kf, E = [], [[] for _ in range(nspin)]
for nk, _, _ in segs:
    for _ in range(nk):
        for s in range(nspin):
            n = int(tok[p]); k = np.array(tok[p + 1:p + 4], float); p += 4
            E[s].append(np.array(tok[p:p + n], float)); p += n
            if s == 0:
                kf.append(k)
E = (np.array(E) - mu) * HA2EV                            # [spin, k, band], E_F = 0
```

运行：

```bash
python3 band_parse.py si.Band --emin -13 --emax 6 --plot si_band.png
```

（示意输出，数值仅供参考）

```text
5 segments, 114 k points, 1 spin channel(s), 26 bands
ticks: L=0.0000  G=0.xxxx  X=x.xxxx  W=x.xxxx  K=x.xxxx  G=x.xxxx
E_F (chemical potential) = -x.xxxx eV (set to 0)
VBM =  -0.xxxx eV at G k=(0.000 0.000 0.000) (spin 0, band 4)
CBM =   0.xxxx eV at k=(0.xxx 0.xxx 0.xxx) (spin 0, band 5)
band gap along the path = 0.5xxx eV (indirect)
smallest direct gap along the path = x.xxxx eV
note: VBM/CBM are searched only on the k points of the path
```

判定逻辑：以化学势为零点，所有 ≤ 0 的本征值视为占据。若某条能带在路径上同时有 ≤ 0 和 > 0 的值，就是穿过费米能的金属/半金属能带，此时不报告带隙。对半导体，OpenMX 在有限电子温度下确定的化学势位于带隙之中，所以"VBM < 0 < CBM"。

仓库还提供功能更全的 [omx_band_plot.py](../scripts/omx_band_plot.py)（自动从 `.dat` 读路径标签、画自旋分辨能带），用 `python ../scripts/omx_band_plot.py -h` 查看用法。浏览器端的 [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS) 可以同时显示能带与 DOS/PDOS、标记带隙与 VBM/CBM，并能读取展开能带。

---

## 8. 自旋极化、非共线与二维体系

### 8.1 铁磁 bcc Fe（共线自旋）

[算例 8-2](../examples/08-band-graphene-fe/README.md) 的 `Fe_bcc_spin_band.dat` 要点：

```text
<Definition.of.Atomic.Species
  Fe   Fe5.5H-s3p2d1   Fe_PBE19H        # 16 valence electrons (3s3p semicore)
Definition.of.Atomic.Species>
<Atoms.SpeciesAndCoordinates
  1  Fe  0.0  0.0  0.0   10.0  6.0      # up/down initial electrons: 4 muB initial moment
Atoms.SpeciesAndCoordinates>
<Atoms.UnitVectors
  -1.435   1.435   1.435
   1.435  -1.435   1.435
   1.435   1.435  -1.435
Atoms.UnitVectors>
scf.SpinPolarization        on
scf.ElectronicTemperature   500.0
scf.Kgrid                   14 14 14
<Band.KPath.UnitCell
  2.870  0.000  0.000
  0.000  2.870  0.000
  0.000  0.000  2.870
Band.KPath.UnitCell>
Band.Nkpath                 6
<Band.kpath
  30   0.00 0.00 0.00    1.00 0.00 0.00    g  H
  21   1.00 0.00 0.00    0.50 0.50 0.00    H  N
  21   0.50 0.50 0.00    0.00 0.00 0.00    N  g
  26   0.00 0.00 0.00    0.50 0.50 0.50    g  P
  26   0.50 0.50 0.50    1.00 0.00 0.00    P  H
  15   0.50 0.50 0.50    0.50 0.50 0.00    P  N
Band.kpath>
```

- 自旋极化计算（`scf.SpinPolarization on`）后，`bandgnu13 fe.Band` 生成 `fe.BANDDAT1`（上自旋）和 `fe.BANDDAT2`（下自旋），`fe.GNUBAND` 会把两者画在一起；
- `band_parse.py` 用红/蓝两色分别画上、下自旋，并提示有能带穿过费米能（金属）；
- 交换劈裂可在 Γ 或 H 点比较同一 d 带上下自旋的能量差来估计；磁矩看 `.out` 中的 Mulliken 分析（第 11 章详述，注意 OpenMX 自旋磁矩的符号约定）。

### 8.2 非共线与 SOC

`scf.SpinPolarization nc`（可加 `scf.SpinOrbit.Coupling on`）时，`.Band` 只有一个通道、本征值个数加倍，`bandgnu13` 只生成 `BANDDAT1`。能带的自旋劈裂、自旋织构分析见[第12章](12-soc-topology.md)。

### 8.3 二维体系：石墨烯

`Graphene_band.dat` 采用 120° 约定的六方晶格（a = 2.46 Å，c 方向 20 Å 真空），两个 C 在分数坐标 (0,0,½) 与 (⅔,⅓,½)：

```text
scf.Kgrid                   18 18 1
Band.Nkpath                 3
<Band.kpath
  30   0.000000 0.000000 0.0    0.500000 0.000000 0.0    g  M
  18   0.500000 0.000000 0.0    0.333333 0.333333 0.0    M  K
  35   0.333333 0.333333 0.0    0.000000 0.000000 0.0    K  g
Band.kpath>
```

要点：真空方向 k 网格取 1、路径第三分量为 0；`scf.Kgrid` 取 3 的倍数时，4.0 默认的 Γ 中心网格（`scf.Generation.Kpoint GCenter`）恰好包含 K 点，有利于描述狄拉克锥附近的占据。π/π* 带在 K 点相接于费米能附近。

> [!NOTE]
> OpenMX 4.0 默认 `scf.Generation.Kpoint GCenter`（总包含 Γ 与 TRIM 点），3.9 使用 `Regular` 网格。用 4.0 复现 3.9 的结果时请显式写 `scf.Generation.Kpoint Regular`（见[第5章](05-scf-kpoints-convergence.md)）。这只影响 SCF 网格，不影响能带路径本身。

---

## 9. 只做后处理：`scf.restart pp`（OpenMX ≥ 4.0）

已有收敛的 SCF（`System.Name_rst/` 存在），想换一条路径或加密分点重新画能带时，不必重跑 SCF：

```bash
cp Si_band.dat Si_band_pp.dat
# keep the SAME System.Name (si) so that si_rst/ is read; then edit Si_band_pp.dat:
#   add            scf.restart   pp
#   change         Band.Nkpath / <Band.kpath  (new path or more points)
mpirun -np 8 openmx Si_band_pp.dat -nt 1 | tee Si_band_pp.std
bandgnu13 si.Band
```

`scf.restart pp`（post processing）读入重启用的哈密顿量和电荷密度，自动设置 `scf.maxIter = 1`、`MD.Type = nomd`，只做一次对角化后进行能带/DOS 等后处理，并且**不更新** `_rst` 中的重启文件（手册 17.1 节）。注意新的 `.Band` 会覆盖旧文件，需要保留时先改名。

在 3.9 中没有 `pp`，常用做法是 `scf.restart on`：从重启文件出发，SCF 一两步即收敛，但会覆盖 `_rst` 文件。

---

## 10. 能带展开（Band Unfolding）

### 10.1 为什么需要展开

用超胞计算缺陷、掺杂、界面或无序时，超胞的布里渊区更小，原胞的能带被"折叠"进去，能带图变成一团乱麻，无法与 ARPES 直接比较。OpenMX 实现了 Lee、Yamada-Takamura 与 Ozaki 提出的 LCAO 能带展开方法：把超胞本征态投影到用户指定的参考胞（通常是原胞）的布里渊区上，给出每个 $(\mathbf{k},E)$ 的谱权重。完美晶体展开后权重为 0 或 1、精确重现原胞能带；引入缺陷后权重变小、能带展宽，反映有效能带结构的"模糊程度"。该功能支持共线、非共线（含 SOC）以及 +U 计算（无需额外关键词），并且可输出按轨道分解的谱权重。4.0 改进了展开代码的内存占用与速度。

### 10.2 关键词

| 关键词 | 默认 | 说明 |
|---|---|---|
| `Unfolding.Electronic.Band` | `off` | 打开能带展开；只支持 `Band` 与 `Cluster` 求解器 |
| `Unfolding.LowerBound` / `Unfolding.UpperBound` | −10 / 10 | 分析的能量窗口（eV，相对费米能） |
| `<Unfolding.ReferenceVectors` | 默认用 `Atoms.UnitVectors` | 参考胞的三个晶格矢量，单位同 `Atoms.UnitVectors.Unit`；不写时不展开，只做轨道组成分析 |
| `<Unfolding.Map` | 每个原子自成一组 | 每行"原子序号 组号"，第一列必须按 1, 2, 3… 列出全部原子；同组原子对应参考胞中的同一个原子位置；组号为正整数，不必从 1 开始、不必升序 |
| `Unfolding.Nkpoint` | — | `<Unfolding.kpoint` 的行数 |
| `<Unfolding.kpoint` | — | 每行"标签 k1 k2 k3"，以**参考胞**的倒格矢为单位；按顺序连接成路径 |
| `Unfolding.desired_totalnkpt` | — | 路径上希望的总 k 点数，程序按等间距划分并保证包含所有列出的点 |
| `<Unfolding.ReferenceOrigin` | 自动 | 参考胞原点（单位同 `Atoms.UnitVectors.Unit`）；默认按"同组原子尽量一一对应到参考胞"和"含原子的参考胞数最少"两条规则自动确定，表面或强畸变体系可能需要手动指定 |

### 10.3 算例：8 原子 Si 立方胞展开到 fcc 原胞

[算例 8-3](../examples/08-unfold-si-supercell/README.md) 的 `Si_supercell_unfold.dat`（关键部分）：

```text
Atoms.Number                8
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Si  0.00  0.00  0.00   2.0 2.0
  2  Si  0.00  0.50  0.50   2.0 2.0
  3  Si  0.50  0.00  0.50   2.0 2.0
  4  Si  0.50  0.50  0.00   2.0 2.0
  5  Si  0.25  0.25  0.25   2.0 2.0
  6  Si  0.25  0.75  0.75   2.0 2.0
  7  Si  0.75  0.25  0.75   2.0 2.0
  8  Si  0.75  0.75  0.25   2.0 2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  5.430  0.000  0.000
  0.000  5.430  0.000
  0.000  0.000  5.430
Atoms.UnitVectors>
scf.EigenvalueSolver        Band
scf.Kgrid                   6 6 6

Unfolding.Electronic.Band   on
Unfolding.LowerBound        -13.0
Unfolding.UpperBound          6.0
<Unfolding.ReferenceVectors
  0.000  2.715  2.715
  2.715  0.000  2.715
  2.715  2.715  0.000
Unfolding.ReferenceVectors>
<Unfolding.Map
  1 1
  2 1
  3 1
  4 1
  5 2
  6 2
  7 2
  8 2
Unfolding.Map>
Unfolding.Nkpoint           6
<Unfolding.kpoint
  L   0.500  0.500  0.500
  G   0.000  0.000  0.000
  X   0.000  0.500  0.500
  W   0.250  0.500  0.750
  K   0.375  0.375  0.750
  G   0.000  0.000  0.000
Unfolding.kpoint>
Unfolding.desired_totalnkpt 120
```

设计思路：

- 参考胞取 fcc 原胞 $\mathbf{a}_1=(0,a/2,a/2)$ 等，与 `Si_band.dat` 完全相同；
- 金刚石结构由两套 fcc 子格组成：原子 1–4 在第一套（组 1），5–8 在第二套（组 2），每组在参考胞中对应一个原子；
- k 点以参考胞（原胞）倒格矢为单位。由惯用胞坐标换算：X = (1,0,0)·2π/a → (0, ½, ½)，W = (1,½,0) → (¼, ½, ¾)，K = (¾,¾,0) → (⅜, ⅜, ¾)，L 仍为 (½, ½, ½)；
- 这样展开路径 L–Γ–X–W–K–Γ 与 `Si_band.dat` 的路径在物理上完全相同，横轴（Bohr⁻¹）也一致，可直接叠加。

掺杂变体 `Si7Ge1_unfold.dat` 把原子 1 换成 Ge，并在 `<Unfolding.Map` 中给它新的组号 3（手册对杂质的推荐做法）。

### 10.4 输出与作图

| 文件 | 内容 |
|---|---|
| `System.Name.unfold_totup`（`_totdn`） | 共线计算的总谱权重；列：路径距离（Bohr⁻¹）、E − E_F（eV）、权重。非自旋计算也写 `_totup` |
| `System.Name.unfold_tot` | 非共线计算的总谱权重 |
| `System.Name.unfold_orbup`（`_orbdn`、`_orb`） | 按轨道分解的谱权重：前两列同上，之后每列一个轨道；列号 = `.out` 中轨道序号 + 2 |
| `System.Name.unfold_plotexample` | 自动生成的 gnuplot 脚本（含高对称点刻度） |

```bash
gnuplot -p -e "load 'si8_unfold.unfold_plotexample'"
python3 plot_unfold.py si8_unfold --band ../08-band-cdia-si/si.BANDDAT1 --out si8_unfold.png
grep -i -n "orbital weights" -A 30 si7ge1_unfold.out | head -40     # orbital order
python3 plot_unfold.py si7ge1_unfold --orb-cols 3,4,5,6 --out ge_orbitals.png
```

强度图（intensity map）：手册提供 `intensity_map.c`，对谱权重做洛伦兹展宽后输出网格数据：

```bash
gcc intensity_map.c -lm -o intensity_map            # in openmx4.0/source
./intensity_map si7ge1_unfold.unfold_totup -c 3 -k 0.1 -e 0.1 -l -13 -u 6 > si7ge1_intmap.txt
# -c column of the weight, -k delta_k (Bohr^-1), -e delta_E (eV), -l/-u energy window (eV)
```

然后在 gnuplot 中 `set pm3d map; splot 'si7ge1_intmap.txt'`。手册强调强度图的绝对值没有物理意义，只看相对分布。

**判据**：完美 Si 的展开结果中，权重 ≈ 1 的点应精确落在原胞能带上，权重 ≈ 0 的点对应折叠进来的"假"能带。若对不上，按顺序检查：参考矢量是否正确、映射是否把同一子格的原子归为同组、k 点是否按参考胞倒格矢书写、两条路径是否一致。

---

## 11. 延伸：费米面与其他工具

- **费米面**：在 DOS 计算（`Dos.fileout on`，较密的 `Dos.Kgrid`，如 61×61×61）中加 `FermiSurfer.fileout on`，会输出 `System.Name.FermiSurf0.bxsf`（XCrySDen）和按原子分解的 `*.frmsf`（FermiSurfer）；自旋极化时有 `FermiSurf0/1`，非共线时为 `FermiSurf.bxsf`。**patch 4.0.1 起 bxsf/frmsf 的单位从 Hartree/Bohr 改为 eV/Å**，旧脚本若仍按原单位处理会得到错误数值。DOS 相关设置见[第9章](09-dos-pdos-cohp-optics.md)。
- **有效质量**：在 VBM/CBM 附近加一段短而密的路径（如 ±0.05 × 2π/a、21 个点），用 `band_parse.py` 读出后做抛物线拟合 $E=\hbar^2k^2/2m^*$。
- **Wannier 插值能带**：见[第14章](14-wannier.md)。
- **能带的轨道成分**：不做展开时（不写 `<Unfolding.ReferenceVectors`），`unfold_orbup` 就是原胞能带的轨道分解（"fat band"），手册 55.1 节用它区分 SiC 的 σ 与 π 带。

---

## 参数选择与收敛建议

| 参数 | 推荐 | 说明 |
|---|---|---|
| `scf.Kgrid` | 按第 5 章收敛（半导体 1 meV/atom 判据，金属更密） | 能带形状对 SCF 网格不太敏感，但费米能位置和金属能带需要充分收敛 |
| 路径分点数 | 每 0.1 × (2π/a) 约 3 点；极值附近加密 | 点太少会漏掉 CBM 位置、能带交叉看起来像反交叉 |
| 基组 | Standard；关心导带时比较 Precise | 能带数 = 基函数数；高能导带依赖基组完备性 |
| `scf.energycutoff` | 200–300 Ry | 与 SCF 一致即可 |
| `scf.ElectronicTemperature` | 半导体 300 K；金属 300–1000 K | 只影响占据与化学势，不改变本征值的计算方式 |
| `Unfolding.desired_totalnkpt` | 100–300 | 强度图需要更多（手册示例用 300） |
| `Unfolding.LowerBound/UpperBound` | 只覆盖关心的窗口 | 窗口越大，输出文件和计算量越大 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 没有生成 `System.Name.Band` | `Band.dispersion` 未开、`Band.Nkpath` 为 0，或求解器不是 `Band`（如 `Cluster`） | 设 `Band.dispersion on`、`scf.EigenvalueSolver Band` |
| `Format error near Band.kpath>` | `Band.Nkpath` 与行数不一致，或某行少了标签/坐标 | 每行 9 项：分点数、6 个坐标、2 个标签 |
| `<Band.kpath is necessary.` | `Band.Nkpath` > 0 但没写 `<Band.kpath` 块（或块名拼错） | 补上块，检查 `<` 和 `>` 的位置 |
| X、L、W 处能带与文献对不上 | k 坐标的单位与预想不同 | fcc/bcc/菱方用 `<Band.KPath.UnitCell`；或用 OpenMX Viewer 生成路径 |
| 六方体系看不到狄拉克锥 | 60°/120° 约定混用，(⅓,⅓,0) 不是 K | 按 §5.2 核对或用自动工具 |
| `bandgnu13: command not found` / 编译失败 | 没有单独编译，或用 `make` 编译 | `gcc bandgnu13.c -lm -o bandgnu13` |
| 自己画的能带比 `bandgnu13` 的图偏移几 eV / 数量级不对 | 直接用了 `.Band` 中的 Hartree 值，没减化学势 | (E − μ) × 27.2113845 eV；或直接用 `BANDDAT` |
| Python 读 `BANDDAT` 后出现竖线 | `BANDDAT` 开头有画分隔线用的数据块 | 跳过横坐标不变的块，或改读 `.Band` |
| 自旋极化只看到一套能带 | 只画了 `BANDDAT1` | 同时画 `BANDDAT2` |
| `scf.restart pp` 后能带全错或报找不到重启文件 | `System.Name` 与 SCF 时不同，或 `_rst` 被删 | 保持相同 `System.Name`，保留 `System.Name_rst/` |
| 高能导带"奇怪"、能带数太少 | LCAO 能带数等于基函数数，高能区依赖基组 | 增大基组（Precise）或加空原子；只解读费米能附近几 eV |
| 致密金属（fcc/bcc）大基组出现伪本征值 | 基组过完备导致数值不稳定 | 换较小的优化基组，或 `scf.ProExpn.VNA off` 并提高截断能（手册 Known issues） |
| `Unfolding.fileout=ON is not supported in case of scf.EigenvalueSolver=...` | 展开只支持 `Band`/`Cluster` | 改求解器 |
| `Format error in Unfolding.Map! (Values cannot be negative.)` | 第一列没有按 1, 2, 3… 顺序列出全部原子，或组号为负 | 每个原子一行、按序号排列 |
| 完美超胞展开后权重不是 0/1 | 参考矢量/映射/k 点单位错误，或参考原点不当 | 逐项检查（§10.4 判据）；表面体系试 `<Unfolding.ReferenceOrigin` |

---

## 练习

1. 删除 `Cdia_band.dat` 中的 `<Band.KPath.UnitCell` 块，其余不变，重新计算并与原图比较，找出"X"实际对应的是哪个点。
2. 用 OpenMX Viewer 打开 `Si_band.dat`，导出带推荐路径的新输入，与本章手写路径比较。
3. 对 Si 在 CBM 附近（Γ–X 上约 0.85 X 处）加一段密集路径，拟合纵向有效质量，并与文献值（约 0.9 mₑ 量级）比较。
4. 用 `scf.restart pp` 在不重跑 SCF 的情况下，把 Si 能带的路径改为 Γ–X–U\|K–Γ–L。
5. 把 `Si7Ge1_unfold.dat` 中的 Ge 换成一个空位（删除原子 1，`Atoms.Number 7`，并相应修改映射），观察带隙中的缺陷态和谱权重变化；比较不加自旋极化与加自旋极化的结果。

---

## 参考资料

**官方文档**

- [OpenMX 4.0 手册 第 22 章 Band structure](https://www.openmx-square.org/openmx_man4.0/s22_band.html) 与 [22.1 自动 k 路径](https://www.openmx-square.org/openmx_man4.0/s22_1_auto-k-path.html)
- [手册 10.2 用 OpenMX Viewer 准备 k 路径](https://www.openmx-square.org/openmx_man4.0/s10_2_kpath.html)、[8.3 cif2omx](https://www.openmx-square.org/openmx_man4.0/s8_3_gen-dat-files.html)
- [手册 8.2 关键词（Band.dispersion、Band.KPath.UnitCell、Band.kpath）](https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html)
- [手册 第 55 章 Band unfolding method](https://www.openmx-square.org/openmx_man4.0/s55_band-unfolding.html)：[55.1](https://www.openmx-square.org/openmx_man4.0/s55_1_analysis.html)、[55.2](https://www.openmx-square.org/openmx_man4.0/s55_2_unfolding.html)、[55.3](https://www.openmx-square.org/openmx_man4.0/s55_3_ref-cell.html)、[55.4 强度图](https://www.openmx-square.org/openmx_man4.0/s55_4_intensity-map.html)
- [手册 17.1 scf.restart（on/pp/c2n）](https://www.openmx-square.org/openmx_man4.0/s17_1_restarting-general.html)、[第 15 章 k 点采样（GCenter/Regular）](https://www.openmx-square.org/openmx_man4.0/s15_k-sampling.html)、[第 65 章 费米面](https://www.openmx-square.org/openmx_man4.0/s65_fermi-surface.html)
- [手册 第 79 章 Known issues（基组过完备）](https://www.openmx-square.org/openmx_man4.0/s79_known-issues.html)
- [Patch 4.0.1 README（bxsf/frmsf 单位变更）](https://www.openmx-square.org/bugfixed/26May08/README.txt)
- [OpenMX Viewer](https://www.openmx-square.org/viewer/)

**相关 GitHub 仓库（仅链接，未复制其代码）**

- [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS) — 单 HTML 浏览器工具：能带 + DOS/PDOS、展开能带、带隙/VBM/CBM 标记
- [Ncmexp2717/OMXTool](https://github.com/Ncmexp2717/OMXTool) — Tcl/Tk 图形工具，预览 `.Band` 并导出 gnuplot 图
- [nnardoto/LowFoldIn](https://github.com/nnardoto/LowFoldIn) — 直接读取 3.9/4.0 的 `.scfout` 做能带与能带展开
- [mitsuaki1987/cif2input](https://github.com/mitsuaki1987/cif2input) — CIF→OpenMX 输入（含原胞与高对称路径），可与 cif2omx 对比
- [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx) — POSCAR→`.dat` 转换，可参考其基组与 k 网格设定逻辑
- [KazumeNishidate/tube](https://github.com/KazumeNishidate/tube) — 纳米管结构生成器，附一维体系能带计算示例
- [hikaruri/DFT_startup](https://github.com/hikaruri/DFT_startup) — 含用 gnuplot 画三维能带的入门资料
- [pranabdas/openmx](https://github.com/pranabdas/openmx) — 英文 OpenMX 教程站，含 Si/金刚石能带示例
- [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath) — HPKOT 约定的高对称路径生成
- 本仓库：[omx_band_plot.py](../scripts/omx_band_plot.py)、[omx_make_input.py](../scripts/omx_make_input.py)、公共篇 [第5章 后处理与可视化](../../Common/tutorials/05-postprocess-visualization.md)、[第3章 结构准备](../../Common/tutorials/03-structure-modeling.md)

**关键文献**

- T. Ozaki, *Phys. Rev. B* 67, 155108 (2003)（变分优化的数值原子轨道）. DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)；T. Ozaki, H. Kino, *Phys. Rev. B* 69, 195113 (2004). DOI: [10.1103/PhysRevB.69.195113](https://doi.org/10.1103/PhysRevB.69.195113)
- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* 49, 299 (2010)（高通量能带路径）. DOI: [10.1016/j.commatsci.2010.05.010](https://doi.org/10.1016/j.commatsci.2010.05.010)
- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* 128, 140 (2017)（seekpath/HPKOT）. DOI: [10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015)
- C.-C. Lee, Y. Yamada-Takamura, T. Ozaki, *J. Phys.: Condens. Matter* 25, 345501 (2013)（LCAO 能带展开，OpenMX 的实现）. DOI: [10.1088/0953-8984/25/34/345501](https://doi.org/10.1088/0953-8984/25/34/345501)
- W. Ku, T. Berlijn, C.-C. Lee, *Phys. Rev. Lett.* 104, 216401 (2010). DOI: [10.1103/PhysRevLett.104.216401](https://doi.org/10.1103/PhysRevLett.104.216401)；V. Popescu, A. Zunger, *Phys. Rev. B* 85, 085201 (2012). DOI: [10.1103/PhysRevB.85.085201](https://doi.org/10.1103/PhysRevB.85.085201)
- M. Kawamura, *Comput. Phys. Commun.* 239, 197 (2019)（FermiSurfer）. DOI: [10.1016/j.cpc.2019.01.017](https://doi.org/10.1016/j.cpc.2019.01.017)

---

[← 上一章](07-md.md) | [返回目录](../README.md) | [下一章 →](09-dos-pdos-cohp-optics.md)
