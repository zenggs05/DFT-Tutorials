# 第9章 态密度 DOS/PDOS、COHP 与光学性质（DOS, PDOS, COHP/COOP & Optical Properties with OpenMX）

本章讲解 OpenMX 4.0 中与"能量分辨的电子结构"有关的后处理：用 `Dos.fileout`、`Dos.Erange`、`Dos.Kgrid` 输出本征值与本征矢，用 `DosMain` 以四面体法或高斯展宽计算总态密度（DOS）和投影态密度（PDOS），读懂 `*.PDOS.Tetrahedron.atomN.d1` 等文件的命名与列含义，并用 Python 按元素和角动量求和；用 4.0 新增的 `cohp` 程序做 COHP/COOP 成键分析；用 CDDF（Kubo–Greenwood）计算光学电导率、介电函数、吸收与反射谱；最后介绍 BoltzTraP/BoltzTraP2 输运接口。算例包括金刚石、单层 MoS₂、反铁磁 NiO 和硅。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 DOS 计算的关键词，理解 `Dos.Erange` 以化学势为零点、`Dos.Kgrid` 应比 `scf.Kgrid` 更密；
- 交互式或用 here-doc 非交互地运行 `DosMain`，选择四面体法/高斯展宽、总 DOS/PDOS 和原子；
- 读懂 PDOS 文件名中的 `s1`、`p1–p3`、`d1–d5`，知道 d 轨道在 OpenMX 中的顺序，正确处理自旋极化文件的列（下自旋带负号）；
- 理解"积分 DOS 从 `Dos.Erange` 下限开始积分"带来的陷阱；
- 用 `COHP.calc.flag`、`<COHP.pairs` 和 `cohp` 后处理得到 −COHP、ICOHP，判断成键/反键；
- 用 `CDDF.start` 计算 ε(ω)、σ(ω)，知道其 k 点收敛要求和"仅支持共线"的限制；
- 把 OpenMX 结果交给 BoltzTraP2 计算 Seebeck 系数与电导率。

| 任务 | 关键词 / 程序 | 本章 | 示例 |
|---|---|---|---|
| 总 DOS、PDOS | `Dos.fileout` + `DosMain` | §3–§6 | [算例 9-1 金刚石](../examples/09-dos-pdos-cdia/README.md) |
| 二维材料、磁性氧化物的轨道分析 | 同上 + Python 求和 | §7–§8 | [算例 9-2 MoS₂ 与 NiO](../examples/09-pdos-mos2-nio/README.md) |
| 大体系 / 很多 k 点 | `DosGauss.fileout` | §9 | — |
| 化学键分析 | `COHP.calc.flag` + `cohp`（≥ 4.0） | §10 | [算例 9-3 Si](../examples/09-cohp-cddf-si/README.md) |
| 光学性质 | `CDDF.start` | §11 | 同上 |
| 热电输运 | BoltzTraP / BoltzTraP2 接口 | §12 | — |

> [!NOTE]
> 本章默认你已能完成 SCF（[第5章](05-scf-kpoints-convergence.md)）与能带计算（[第8章](08-band.md)）。DOS 与能带最好在同一组参数下计算并画在一起；费米面输出（`FermiSurfer.fileout`）在第 8 章 §11 介绍，Mulliken/Voronoi 电荷在[第10章](10-charge-mo-cube.md)。

---

## 2. 原理速览

### 2.1 态密度与展宽方法

$$
g(E)=\sum_{n}\int_\mathrm{BZ}\frac{d^3k}{V_\mathrm{BZ}}\,\delta\big(E-\varepsilon_{n\mathbf{k}}\big)
$$

离散 k 网格上有两种近似 δ 函数的办法：

- **四面体法（tetrahedron）**：把布里渊区剖分为四面体，在四面体内线性插值本征值。没有人为展宽参数，带边、范霍夫奇点更锐利，是周期体系的首选；
- **高斯展宽**：OpenMX 的高斯函数定义为 $\exp[-(E/a)^2]$，$a$（eV）由用户给出。分子（只有 Γ 点）和 O(N) 计算只能用这种方法。

### 2.2 OpenMX 的 PDOS 是 Mulliken 型投影

在非正交的 PAO 基组下，第 $i$ 个原子第 $\alpha$ 个轨道的 PDOS 定义为

$$
g_{i\alpha}(E)=\sum_{n\mathbf{k}}w_\mathbf{k}\,\mathrm{Re}\Big[c^{*}_{n\mathbf{k},i\alpha}\big(S(\mathbf{k})c_{n\mathbf{k}}\big)_{i\alpha}\Big]\,\delta(E-\varepsilon_{n\mathbf{k}})
$$

对所有轨道求和恰好等于总 DOS；积分到费米能即 Mulliken 布居。因为是 Mulliken 划分，个别能量处可能出现很小的负值，s/p/d 的积分值也不是整数——手册（第 23.1 节）特别提醒，金刚石 s 与 p 轨道积分到费米能并不严格等于 1 和 3，因为电荷在 s、p、d 之间重新分配。

### 2.3 COHP 与 COOP

把能量分辨的密度矩阵 $\rho_{i\alpha j\beta}(E)$ 分别乘以哈密顿矩阵元与重叠矩阵元：

$$
\mathrm{COHP}_{i\alpha,j\beta}(E)=\rho_{i\alpha j\beta}(E)\,h_{i\alpha j\beta},\qquad
\mathrm{COOP}_{i\alpha,j\beta}(E)=\rho_{i\alpha j\beta}(E)\,s_{i\alpha j\beta}
$$

对所有原子对、轨道、晶格平移和自旋求和并积分到化学势，COHP 给出能带能量，COOP 给出电子数。**符号约定（手册第 24.1 节）：COHP 为负表示成键、为正表示反键；COOP 为正表示成键、为负表示反键**。作图时习惯画 −COHP，使成键态位于右侧。

### 2.4 光学电导率与介电函数

OpenMX 在独立粒子近似下用 Kubo–Greenwood 公式计算光学电导率张量 $\sigma_{\alpha\beta}(\omega)$（由动量矩阵元与带间跃迁得到），再给出介电函数

$$
\varepsilon_{\alpha\beta}(\omega)=\delta_{\alpha\beta}+\frac{i\,\sigma_{\alpha\beta}(\omega)}{\varepsilon_0\,\omega}
$$

以及由它导出的折射率 $n+i\kappa=\sqrt{\varepsilon}$、消光、吸收、透射和反射。不包含局域场效应与激子效应。

---

## 3. DOS 计算的流程与关键词

```text
.dat (Dos.fileout on) --openmx--> System.Name.Dos.val (text) + System.Name.Dos.vec (binary)
                      --DosMain--> System.Name.DOS.Tetrahedron, System.Name.PDOS.Tetrahedron.atomN[.lm] ...
```

| 关键词 | 默认 | 说明 |
|---|---|---|
| `Dos.fileout` | `off` | 输出 DOS 所需的本征值（`.Dos.val`，文本）和本征矢（`.Dos.vec`，二进制） |
| `Dos.Erange` | `-20 20` | 能量范围（eV），**以化学势为零点**；下限决定积分 DOS 的起点（§6） |
| `Dos.Kgrid` | 等于 `scf.Kgrid` | DOS 用的 k 网格，应比 SCF 更密；SCF 收敛后只做一次对角化 |
| `DosGauss.fileout` | `off` | 在 OpenMX 内部即时做高斯展宽，不写 `.Dos.vec`，省内存（§9）；与 `Dos.fileout` 互斥 |
| `DosGauss.Num.Mesh` | 200 | `DosGauss` 模式下 `Dos.Erange` 的划分数 |
| `DosGauss.Width` | 0.2 | `DosGauss` 模式下的展宽 $a$（eV） |

DOS 部分可以与 SCF 在同一次运行中完成；若 SCF 已收敛，也可以保持 `System.Name` 不变、加 `scf.restart pp`（OpenMX ≥ 4.0）只做一次对角化来输出 DOS 数据（见[第8章 §9](08-band.md)）。

> [!TIP]
> `.Dos.vec` 的大小约为（k 点数 × 能量窗口内的本征态数 × 基函数数）个数，密网格 + 大体系时可达数十 GB，`DosMain` 读入时也要同样多的内存。只需要总 DOS 或高斯展宽 PDOS 时，用 `DosGauss.fileout`；只关心费米能附近时，缩小 `Dos.Erange`。

---

## 4. 算例一：金刚石的 DOS 与 PDOS

完整文件见 [../examples/09-dos-pdos-cdia/README.md](../examples/09-dos-pdos-cdia/README.md)。

### 4.1 完整输入 `Cdia_dos.dat`

```text
# Diamond (C): SCF + data for total DOS and PDOS (Dos.fileout). OpenMX 4.0.
System.CurrrentDirectory    ./
System.Name                 cdia_dos
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

Dos.fileout                 on
Dos.Erange                  -30.0  20.0    # eV, relative to the chemical potential; lower bound below all valence bands
Dos.Kgrid                   16 16 16       # denser than scf.Kgrid
```

| 关键词 | 取值 | 理由 |
|---|---|---|
| `Dos.Erange` | −30 ~ 20 eV | 金刚石价带宽度约 21–22 eV，化学势位于带隙中，下限取 −30 eV 才能完整包含价带底，积分 DOS 在费米能处才等于 8 |
| `Dos.Kgrid` | 16×16×16 | 四面体法对网格要求比 SCF 高；收敛判据见本章末表 |
| `scf.Kgrid` | 10×10×10 | 只用于 SCF |

### 4.2 运行 OpenMX

```bash
cd OpenMX/examples/09-dos-pdos-cdia
mpirun -np 8 openmx Cdia_dos.dat -nt 1 | tee cdia_dos.std
ls -lh cdia_dos.Dos.val cdia_dos.Dos.vec
```

### 4.3 编译 `DosMain`

```bash
cd /path/to/openmx4.0/source
make DosMain            # or "make all" to build all post-processing tools
cp DosMain ~/bin/       # or copy it to the working directory
```

用官方 deb 包安装时 `DosMain` 已在 `/opt/openmx4.0` 下，但直接运行可能报 `error while loading shared libraries: libmpifort.so.12`，先执行：

```bash
export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
```

然后直接调用 `DosMain`（不要写成 `./DosMain`）。

### 4.4 交互式运行 `DosMain`

```bash
DosMain cdia_dos.Dos.val cdia_dos.Dos.vec
```

程序依次提问（下面的提示文字摘自 4.0 `DosMain.c` 中的打印语句，"→" 后是要输入的内容）：

```text
Which method do you use?, Tetrahedron(1), Gaussian Broadening(2)   -> 1
Do you want Dos(1) or PDos(2)?                                       -> 2
Which atoms for PDOS : (1,...,2), e.g. 1 2                           -> 1 2
```

- 选高斯展宽（2）时还会问 `Please input a value of gaussian (double) (eV)`，输入 $a$，如 `0.1`；
- 若 `.Dos.val` 中的 k 网格是 1×1×1（分子、Cluster 计算），或数据来自 DC 等 O(N) 方法，程序自动采用高斯展宽，**不再询问方法**，但仍询问展宽值；
- 体系只有 1 个原子时不询问原子编号；
- 一次运行只能选 DOS 或 PDOS，二者都要就运行两次。

### 4.5 非交互运行（作业脚本中必需）

`DosMain` 从标准输入读取回答，可以用 here-doc 或答案文件：

```bash
# total DOS, tetrahedron
DosMain cdia_dos.Dos.val cdia_dos.Dos.vec << EOF
1
1
EOF

# PDOS of atoms 1 and 2, tetrahedron
DosMain cdia_dos.Dos.val cdia_dos.Dos.vec < dosmain_answers.txt   # file: 1 / 2 / "1 2"
```

示例目录中的 `run_dosmain.sh` 会先读取 `.Dos.val` 头部的 `mode`、`Kgrid`、`atomnum`，自动决定要回答哪些问题，然后依次生成总 DOS 与 PDOS：

```bash
bash run_dosmain.sh cdia_dos tetra            # tetrahedron, all atoms
bash run_dosmain.sh cdia_dos gauss 0.1 1      # Gaussian a = 0.1 eV, PDOS of atom 1 only
```

仓库脚本 ../scripts/omx_dos_auto.py 把"驱动 DosMain + 按元素/角动量求和 + 作图"合在一起，用 `python ../scripts/omx_dos_auto.py -h` 查看用法。

### 4.6 输出文件与列

| 文件 | 内容 |
|---|---|
| `cdia_dos.DOS.Tetrahedron` | 总 DOS |
| `cdia_dos.PDOS.Tetrahedron.atom1` | 原子 1 所有轨道之和 |
| `cdia_dos.PDOS.Tetrahedron.atom1.s1` | 原子 1 的 s 轨道（所有 s 径向函数之和） |
| `cdia_dos.PDOS.Tetrahedron.atom1.p1` / `p2` / `p3` | px / py / pz |
| `cdia_dos.PDOS.Tetrahedron.atom1.d1` … `d5` | 5 个 d 轨道（顺序见 §5.1） |
| `*.DOS.Gaussian`、`*.PDOS.Gaussian.*` | 选择高斯展宽时的对应文件 |

列的含义（手册 23.1 节；自旋列的符号据 4.0 `DosMain.c` 源码）：

| 计算类型 | 第 1 列 | 第 2 列 | 第 3 列 | 第 4 列 | 第 5 列 |
|---|---|---|---|---|---|
| 非自旋 / 非共线 | E − μ（eV） | DOS（states/eV，已含自旋简并因子 2） | 积分 DOS | — | — |
| 共线自旋极化 | E − μ（eV） | 上自旋 DOS | **下自旋 DOS 取负号** | 上自旋积分 | 下自旋积分 |

> [!WARNING]
> 手册只写"第 2、3 列分别为上、下自旋"，而 4.0 `DosMain.c` 写出的下自旋 DOS 带负号（便于直接画成上下对称的图）。自己处理数据时先 `head` 看一眼符号，求和、积分时对下自旋取绝对值。

### 4.7 检查结果

```bash
# integrated DOS at E = 0 (column 3), should be ~8 electrons for diamond
awk '$1<=0 {n=$3} END{print "N(E_F) =", n}' cdia_dos.DOS.Tetrahedron
python3 sum_pdos.py cdia_dos --dat Cdia_dos.dat --plot --emin -28 --emax 15
```

（示意输出，数值仅供参考）

```text
method: Tetrahedron, 18 files, spin-polarized: False
C d                          -> cdia_dos.pdos_C_d.dat              N(E<E_F) = 0.0xx
C p                          -> cdia_dos.pdos_C_p.dat              N(E<E_F) = x.xxx
C s                          -> cdia_dos.pdos_C_s.dat              N(E<E_F) = x.xxx
max |sum(PDOS) - DOS| = 0.00xx states/eV (small only if PDOS was done for ALL atoms)
```

判据：

- 总 DOS 在约 4 eV 宽的范围内为零（PBE 带隙，见第 8 章）；价带宽约 21–22 eV，化学势大致位于带隙中间，所以价带底出现在约 −23 eV 附近；
- s、p、d 积分之和约为 8（两个 C 原子），与 `.out` 中的 Mulliken 布居一致；d 轨道（极化函数）有少量布居是正常的；
- 对**所有原子**做了 PDOS 时，PDOS 之和应等于总 DOS（脚本打印的最大偏差应很小）。

---

## 5. PDOS 的轨道命名与求和

### 5.1 轨道顺序

文件名中的数字是该角动量下的磁量子数序号（从 1 开始），所有径向函数（多重 zeta）已经求和。OpenMX 使用实球谐函数，顺序与 `.out` 中 "Decomposed Mulliken populations" 表一致（核对自 4.0 源码 `Mulliken_Charge.c`）：

| 文件后缀 | 轨道 |
|---|---|
| `s1` | s |
| `p1`, `p2`, `p3` | px, py, pz |
| `d1` … `d5` | d(3z²−r²), d(x²−y²), dxy, dxz, dyz |
| `f1` … `f7` | f(5z²−3r²)z 类型等，按 `.out` 表中的顺序 |

轨道的 x、y、z 指笛卡尔坐标轴，即 `Atoms.UnitVectors` 所在的坐标系。分析晶体场劈裂时，要让坐标轴与局域配位多面体对齐：例如 NiO 的晶格矢量写在立方坐标系中，$e_g$ = d1 + d2，$t_{2g}$ = d3 + d4 + d5。

### 5.2 用 Python 按元素和角动量求和

示例目录的 `sum_pdos.py` 读取 `.dat` 得到"原子 → 物种 → 元素"的对应关系，然后把同一元素（或 `--by species` 同一物种、`--by atom` 单个原子）同一角动量的文件相加，输出 `前缀.pdos_元素_l.dat`，并打印费米能以下的积分电子数。核心逻辑：

```python
import glob, re, numpy as np
pat = re.compile(r"\.atom(\d+)\.([spdfg])(\d+)$")
elem = {1: "C", 2: "C"}                          # atom index -> element (read from the .dat)
sums = {}
for f in glob.glob("cdia_dos.PDOS.Tetrahedron.atom*.*"):
    m = pat.search(f)
    if not m:
        continue                                 # skip the per-atom total files
    atom, l = int(m.group(1)), m.group(2)
    d = np.loadtxt(f)
    dos = d[:, 1:3] if d.shape[1] >= 5 else d[:, 1:2]   # spin: up, -down
    key = (elem[atom], l)
    sums[key] = sums.get(key, 0.0) + dos
energy = d[:, 0]
```

常用命令：

```bash
python3 sum_pdos.py cdia_dos --dat Cdia_dos.dat --plot
python3 sum_pdos.py nio_afm  --dat NiO_AFM_pdos.dat --by atom --resolve-m --plot   # m-resolved d orbitals
```

---

## 6. 积分 DOS 的陷阱

`DosMain` 输出的积分 DOS 是从 `Dos.Erange` 的**下限**开始累加的，因此：

1. 若下限高于某些占据态（典型是半芯态：Mo 4s/4p、Ni 3p、Fe_PBE19H 的 3s/3p 等，它们位于费米能以下数十 eV），这些电子不会被计入，0 eV 处的积分值会小于价电子总数；
2. 高斯展宽的尾巴会让积分值在能量窗口边缘不准确；
3. 自旋极化时上、下自旋分别积分。

正确做法：

- 需要"总电子数自洽检查"时，把 `Dos.Erange` 下限设到所有占据态以下（或者只比较两个能量之间积分值的**差**）；
- 求某个轨道/原子的占据数，用 `.out` 中的 Mulliken 布居（或第 10 章的其他电荷分析），而不是积分 PDOS；
- 例如单层 MoS₂（Mo_PBE19 有 14 个价电子、S_PBE19 有 6 个）在 `Dos.Erange -20 10` 下，0 eV 处的积分值会明显小于 26，差值大致对应被排除在窗口外的 Mo 4s/4p 半芯电子——这不是错误。

---

## 7. 算例二：单层 MoS₂ 的 PDOS

完整文件见 [../examples/09-pdos-mos2-nio/README.md](../examples/09-pdos-mos2-nio/README.md)。关键部分：

```text
<Definition.of.Atomic.Species
  Mo   Mo7.0-s3p2d2     Mo_PBE19
  S    S7.0-s2p2d1f1    S_PBE19
Definition.of.Atomic.Species>
<Atoms.SpeciesAndCoordinates
  1  Mo   0.00000   1.82443   10.00000   7.0  7.0
  2  S    1.58000   0.91221   11.56000   3.0  3.0
  3  S    1.58000   0.91221    8.44000   3.0  3.0
Atoms.SpeciesAndCoordinates>
<Atoms.UnitVectors
   3.16000   0.00000   0.00000
  -1.58000   2.73664   0.00000
   0.00000   0.00000  20.00000
Atoms.UnitVectors>
scf.Kgrid                   12 12 1
Dos.fileout                 on
Dos.Erange                  -20.0  10.0
Dos.Kgrid                   30 30 1
```

```bash
cd OpenMX/examples/09-pdos-mos2-nio
mpirun -np 8 openmx MoS2_monolayer_pdos.dat -nt 1 | tee mos2.std
bash ../09-dos-pdos-cdia/run_dosmain.sh mos2 tetra 0.1
python3 ../09-dos-pdos-cdia/sum_pdos.py mos2 --dat MoS2_monolayer_pdos.dat --resolve-m --plot --emin -8 --emax 5
```

判读：

- 三棱柱配位下 Mo 4d 分裂为 $a_1'$（d3z²−r²，即 d1）、$e'$（dx²−y²、dxy，即 d2、d3）和 $e''$（dxz、dyz，即 d4、d5）。价带顶与导带底附近 Mo d 占主导，并与 S p 杂化；
- PBE 下单层 MoS₂ 是 K 点的直接带隙，文献中 PBE 值约 1.6–1.7 eV（依赖晶格常数），实验光致发光约 1.9 eV；本例晶格常数取近似实验值、未弛豫，数值只作参考；
- 二维体系的 DOS 单位是"每个超胞"，与真空厚度无关，但 k 网格沿 c 只能取 1。

---

## 8. 算例三：反铁磁 NiO 的自旋分辨 PDOS

NiO 采用沿 [111] 交替的 AFM-II 序，用含 2 个 Ni、2 个 O 的菱方晶胞：

```text
<Definition.of.Atomic.Species
  Ni   Ni6.0S-s3p2d1   Ni_PBE19S        # 16 valence electrons
  O    O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Ni  0.00  0.00  0.00   9.0  7.0
  2  Ni  0.50  0.50  0.50   7.0  9.0
  3  O   0.25  0.25  0.25   3.0  3.0
  4  O   0.75  0.75  0.75   3.0  3.0
Atoms.SpeciesAndCoordinates>
<Atoms.UnitVectors
  4.170  2.085  2.085
  2.085  4.170  2.085
  2.085  2.085  4.170
Atoms.UnitVectors>
scf.SpinPolarization        on
Dos.Erange                  -20.0  10.0
Dos.Kgrid                   14 14 14
```

- 反铁磁初态由坐标行最后两列给出：Ni1 上/下 = 9/7（+2 μB），Ni2 = 7/9（−2 μB），两者之和都等于 `Ni_PBE19S` 的 16 个价电子；O 为 3/3；
- 晶格矢量 $a(1,\tfrac12,\tfrac12)$ 等都是面心立方格矢，晶胞体积为立方胞的一半，含 2 个 NiO 式量；Ni1 与 Ni2 分别位于相邻的 (111) 面上；
- 晶格矢量写在立方坐标系中，因此 d1 + d2 是 $e_g$、d3 + d4 + d5 是 $t_{2g}$。

```bash
mpirun -np 16 openmx NiO_AFM_pdos.dat -nt 1 | tee nio_afm.std
grep -i "spin moment" nio_afm.out | tail -n 3     # then inspect the Mulliken section of nio_afm.out
bash ../09-dos-pdos-cdia/run_dosmain.sh nio_afm tetra 0.1
python3 ../09-dos-pdos-cdia/sum_pdos.py nio_afm --dat NiO_AFM_pdos.dat --by atom --resolve-m --plot --emin -9 --emax 5
```

判读（量级参考）：

- 总自旋磁矩应接近 0（反铁磁），两个 Ni 的 Mulliken 磁矩大小相等、符号相反（注意 OpenMX 自旋磁矩的符号约定，见[第11章](11-spin-noncollinear.md)）；
- Ni1 的 $e_g$ 上自旋占满、下自旋空，$t_{2g}$ 两个自旋都占据；Ni2 则上下互换；价带顶有明显的 O 2p 与 Ni 3d 杂化；
- GGA（无 +U）严重低估 NiO 带隙（文献 PBE 值约 0.5–1 eV 量级，实验约 4 eV）和磁矩（PBE 约 1.2–1.4 μB，实验约 1.6–1.9 μB）。加 Hubbard U 的做法见[第13章](13-dft-u.md)。

---

## 9. 大体系、很多 k 点与分子

**`DosGauss.fileout`**（手册 23.2 节）：OpenMX 在计算过程中直接做高斯展宽，不写巨大的 `.Dos.vec`：

```text
Dos.fileout          off        # must be off: the two are mutually exclusive
DosGauss.fileout     on
DosGauss.Num.Mesh    400        # divisions of Dos.Erange (default 200)
DosGauss.Width       0.1        # eV, a in exp(-(E/a)^2) (default 0.2)
Dos.Erange           -20.0 10.0
Dos.Kgrid            12 12 12
```

之后仍用 `DosMain` 做 DOS/PDOS（它识别出数据来自高斯展宽模式，不再询问方法），但**不能**用四面体法，也不会生成 FermiSurfer 文件。该模式只支持 Band 计算。

**分子（Cluster）**：`.Dos.val` 中 k 网格为 1×1×1，`DosMain` 自动采用高斯展宽，只询问展宽值。分子"DOS"本质上是展宽后的分子轨道能级谱，展宽一般取 0.1–0.3 eV；HOMO/LUMO 的能量与形状见第 10 章的 `MO.fileout`。

**O(N) 方法（DC 等）**：同样只能用高斯展宽（手册说明 O(N) 计算也支持 DOS 输出）。

---

## 10. COHP 与 COOP（OpenMX ≥ 4.0）

### 10.1 关键词

COHP/COOP 与 DOS 同时计算，所以必须打开 `Dos.fileout`：

| 关键词 | 默认 | 说明 |
|---|---|---|
| `Dos.fileout` | `off` | 必须为 `on`，否则报 `In case of the COHP calculation, Dos.fileout should be switched on.` |
| `Dos.Erange` / `Dos.Kgrid` | — | 也决定 COHP 的能量范围和 k 网格 |
| `COHP.calc.flag` | `off` | 打开 COHP/COOP（COOP 自动同时计算，没有单独的关键词） |
| `COHP.num.pairs` | 1 | 原子对数目 |
| `<COHP.pairs` | — | 每行 `序号 原子1 原子2 l1 l2 l3`，原子 2 位于平移 $l_1\mathbf{a}_1+l_2\mathbf{a}_2+l_3\mathbf{a}_3$ 的晶胞中 |

### 10.2 算例：Si 的最近邻与次近邻

`Si_COHP.dat`（Si 原胞，$\mathbf{a}_1=(0,a/2,a/2)$ 等，原子 1 在原点、原子 2 在 (¼,¼,¼)）中：

```text
Dos.fileout                 on
Dos.Erange                  -15.0  10.0
Dos.Kgrid                   16 16 16

COHP.calc.flag              on
COHP.num.pairs              3
<COHP.pairs
  1   1  2   0  0  0
  2   1  2  -1  0  0
  3   1  1   1  0  0
COHP.pairs>
```

- 对 1：原子 1 与同一晶胞中的原子 2，距离 $\sqrt3a/4\approx2.35$ Å（最近邻）；
- 对 2：原子 1 与 $-\mathbf{a}_1$ 晶胞中的原子 2，位置 $(a/4,-a/4,-a/4)$，也是最近邻，结果应与对 1 相同（对称性检查）；
- 对 3：原子 1 与 $+\mathbf{a}_1$ 晶胞中的原子 1，距离 $a/\sqrt2\approx3.84$ Å（次近邻）。

### 10.3 运行与 `cohp` 后处理

```bash
cd OpenMX/examples/09-cohp-cddf-si
mpirun -np 8 openmx Si_COHP.dat -nt 1 | tee si_cohp.std
ls si_cohp.cohp                              # binary file written after the run

cd /path/to/openmx4.0/source && make cohp && cp cohp ~/bin/ && cd -
cohp si_cohp.cohp
```

交互提问（提示文字摘自 4.0 `cohp.c`）：

```text
Which method do you use?, Tetrahedron(1), Gaussian Broadening(2)   -> 1
Please input the number of energy meshes                             -> 3000
Which pair for COHP : (1 to 3)                                       -> 1
```

选高斯展宽时，在能量网格数之前还会询问展宽值（eV）。每次运行处理一个原子对，批处理：

```bash
for PAIR in 1 2 3; do printf '1\n3000\n%d\n' $PAIR | cohp si_cohp.cohp > cohp_pair$PAIR.log; done
```

### 10.4 输出文件

以 `System.Name = si_cohp`、原子对 1、自旋 0、四面体法为例：

| 文件 | 含义 |
|---|---|
| `si_cohp_i1_s0_tetra.cohp` | 两个原子所有轨道求和：$\mathrm{COHP}_{ij}(E)$ |
| `si_cohp_i1_o0_s0_tetra.cohp` … | 原子 1 的第 0、1、2…个基函数对原子 2 所有基函数求和 |
| `si_cohp_i1_o0_o0_s0_tetra.cohp` … | 完全轨道分辨 |
| 同名 `.coop` | COOP，格式相同 |
| 高斯展宽时 | `_tetra` 换成 `_gauss` |

列：能量（eV，化学势为零点）、COHP（无量纲）、积分到该能量的 COHP。轨道编号 `o0, o1, …` 从 0 开始，按基函数顺序（s 的各径向函数、p 的各径向函数 × px/py/pz……），可对照 `.out` 中的轨道分解表。

### 10.5 作图与判读

```bash
python3 plot_cohp_cddf.py cohp si_cohp_i1_s0_tetra.cohp si_cohp_i2_s0_tetra.cohp si_cohp_i3_s0_tetra.cohp --out si_cohp.png
```

（示意输出，数值仅供参考）

```text
si_cohp_i1_s0_tetra.cohp                      ICOHP(E_F) =  -x.xxxx   (bonding if negative)
si_cohp_i2_s0_tetra.cohp                      ICOHP(E_F) =  -x.xxxx   (bonding if negative)
si_cohp_i3_s0_tetra.cohp                      ICOHP(E_F) =  -0.xxxx   (bonding if negative)
```

- Si 最近邻：化学势以下几乎全是成键态（−COHP > 0），导带中出现反键态——这正是共价半导体的特征（手册图 19 给出了同样的结论）；
- 对 1 与对 2 的 ICOHP 应相同；次近邻（对 3）的 |ICOHP| 小得多；
- 手册还给出 ScN、TiN、VN 的例子（`work/ScN_COHP.dat` 等）：随 d 电子增加，费米能附近开始占据反键态。

> [!NOTE]
> COHP 的绝对数值依赖于基组（非正交 PAO 下矩阵元的划分方式），更适合比较同一套基组下不同键、不同体系的**趋势**。与平面波代码 + LOBSTER 的 ICOHP 数值不要直接对比。

---

## 11. 光学性质：CDDF（光学电导率与介电函数）

### 11.1 关键词

| 关键词 | 默认 | 说明 |
|---|---|---|
| `CDDF.start` | `off` | 打开光学电导率与介电函数计算 |
| `CDDF.FWHM` | 0.2 | 洛伦兹展宽的半高全宽（eV） |
| `CDDF.maximum_energy` | 10.0 | 输出的最大光子能量（eV），范围从 0 开始 |
| `CDDF.additional_maximum_energy` | 0.0 | 额外纳入计算的能量范围（eV），考虑窗口外态因展宽带来的贡献 |
| `CDDF.frequency.grid.total_number` | 10000 | 频率网格点数；步长 = 最大能量 / 点数 |
| `CDDF.Kgrid` | 等于 `scf.Kgrid` | 光学计算的 k 网格，可以（也应该）比 SCF 密得多 |
| `CDDF.shift_ChemP` | 0.0 | **4.0 新增**：把化学势平移 μ → μ + shift（eV），可用于粗略考察载流子掺杂对跃迁的影响；它不是"剪刀算符" |

> [!IMPORTANT]
> - 4.0 删除了 `CDDF.material_type`（3.9 中用于区分绝缘体/金属）。旧输入中的该行在 4.0 中不再被读取，建议删掉，以免误以为它仍起作用。
> - 手册 62.1 节明确：**目前只支持共线计算**（非自旋或 `scf.SpinPolarization on`），非共线扩展"将在近期实现"。社区报告 4.0 中 nc + SOC + CDDF 会段错误或报 `Parameter 10 was incorrect on entry to ZHEEVX` / `Parameter 9 was incorrect on entry to DSTEVX`，分析与补丁见 [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix)，或等待官方补丁。
> - SCF 用 Γ 点（`scf.Kgrid 1 1 1`，Cluster）而 `CDDF.Kgrid` 不是 1×1×1 时，程序报 `You want to 'Cluster' calculation for SCF, but 'Band' calculation for CDDF calculation.` 并退出：SCF 与光学计算必须用同一种求解器。

### 11.2 算例：硅

`Si_CDDF.dat` 与 `Si_COHP.dat` 结构相同，光学部分：

```text
CDDF.start                         on
CDDF.FWHM                          0.2
CDDF.maximum_energy               10.0
CDDF.additional_maximum_energy     1.0
CDDF.frequency.grid.total_number  2000
CDDF.Kgrid                        30 30 30
CDDF.shift_ChemP                   0.0
```

```bash
mpirun -np 32 openmx Si_CDDF.dat -nt 1 | tee si_cddf.std
grep -n "Optical calculations" si_cddf.std        # start / end (with elapsed time)
ls si_cddf.cd_re si_cddf.cd_im si_cddf.df_re si_cddf.df_im si_cddf.absorption \
   si_cddf.extinction si_cddf.transmission si_cddf.reflection si_cddf.refractive_index
head -n 4 si_cddf.df_im                           # header gives quantity, unit and column index
```

SCF 收敛后标准输出会打印 `<Optical calculations start>`，结束时打印 `<Optical calculations end, time=... (s)>`。输出文件（手册 62.2 节）：

| 文件后缀 | 内容 |
|---|---|
| `.cd_re` / `.cd_im` | 光学电导率张量实部 / 虚部（S/m） |
| `.df_re` / `.df_im` | 介电函数张量实部 / 虚部 |
| `.absorption`、`.extinction`、`.transmission`、`.reflection`、`.refractive_index` | 由电导率导出的各量（单位见文件头） |

每个文件的列：光子能量（eV）、xx、xy、xz、yx、yy、yz、zx、zy、zz、(xx+yy+zz)/3。

```bash
python3 plot_cohp_cddf.py optics si_cddf --out si_optics.png
```

（示意输出，数值仅供参考）

```text
energy grid: 0.0000 - 10.0000 eV, 2001 points
eps1 at the lowest energy (0.000 eV) = 1x.xxx  (electronic, no local-field effects)
eps2 peaks near: 3.xx eV, 4.xx eV
```

### 11.3 判读与收敛

- 实验上 Si 的 ε₂ 有两个主峰 E₁ ≈ 3.4 eV、E₂ ≈ 4.3 eV（Aspnes & Studna）。PBE 低估带隙，计算谱整体红移约 0.5 eV；ε₁(0)（电子部分的静态介电常数）实验约 11.7–11.9，独立粒子近似下的 PBE 结果通常偏大一些；
- **k 点收敛最关键**：手册 62.4 节对 Si 从 10×10×10 测到 100×100×100，指出需要很密的网格才能使介电函数收敛；本例用 30×30×30 只为控制计算量，正式工作至少测试 50×50×50；
- `CDDF.FWHM` 越大谱越平滑但峰越宽；比较不同 k 网格时保持 FWHM 不变；
- 金属的带内（Drude）贡献需要特别注意，解读低能区前请阅读手册 62 章的说明与示例（`work/cddf_example/Febcc-Col_k30x30x30.dat` 等）；
- 安装后可用自动测试检查光学模块：`mpirun -np 56 ./openmx -runtestCDDF -nt 2`，结果在 `runtestCDDF.result`，介电函数差值 ≤ 0.1 即视为正常。

---

## 12. 输运：BoltzTraP / BoltzTraP2 接口

OpenMX 提供到 BoltzTraP（`work/MX_TRAP.sh`）和 4.0 新增的 BoltzTraP2（`work/MX_TRAP2.sh`）的接口，在恒定弛豫时间近似下计算 Seebeck 系数、电导率 σ/τ、电子热导率与霍尔系数。BoltzTraP2 接口支持非自旋、共线自旋与非共线计算。以下流程按手册 63.2 节：

```bash
# 1) SCF with a dense k grid; follow the settings of work/Si_BoltzTraP.dat
#    (the converter reads the eigenvalues from the .out file)
mpirun -np 28 openmx Si_BoltzTraP.dat | tee Si_BoltzTraP.std

# 2) BoltzTraP2 + the OpenMX loader (once)
conda create --name boltztrap2 -c conda-forge
conda activate boltztrap2
conda install -c conda-forge boltztrap2
#    copy work/dft_forOpenMXloader.py -> <env>/lib/python3.x/site-packages/BoltzTraP2/dft.py
#    copy work/io_forOpenMXloader.py  -> <env>/lib/python3.x/site-packages/BoltzTraP2/io.py
#    (back up the original dft.py / io.py first; with NumPy 2.x replace np.infty by np.inf in io.py)

# 3) convert: answer the prompt with Si_BoltzTraP.out
bash MX_TRAP2.sh          # -> Si_BoltzTraP/{Si_BoltzTraP.MXenergy, .MXstructure, .out}

# 4) interpolate and integrate (commands as given in the OpenMX manual)
btp2 interpolate -m 25 Si_BoltzTraP
cp interpolation.bt2 Si_BoltzTraP.bt2
btp2 integrate -t Si_BoltzTraP.bt2 300
#    -> Si_BoltzTraP.trace / .condtens / .halltens / .btj
```

要点：

- 加载器使用 `.MXenergy`、`.MXstructure` 扩展名，以免与 BoltzTraP2 自带的加载器冲突；
- 插值倍数（`-m`）、温度和近似方式都会影响结果，要做收敛检查；不同 BoltzTraP2 版本的命令行参数可能不同，请用 `btp2 interpolate -h`、`btp2 integrate -h` 核对；
- 旧版 BoltzTraP 1.2.5 路线：`sh MX_TRAP.sh` 生成 `.energy`（自旋极化为 `.energyup/.energydn`、非共线为 `.energyso`）、`.struct`、`.intrans`，再运行 `x_trans BoltzTraP`（相应加 `-up`、`-dn`、`-so`）；
- 透射谱型（NEGF）输运与热电见[第17章](17-negf-transport.md)。

---

## 参数选择与收敛建议

| 参数 | 推荐 | 检查方法 |
|---|---|---|
| `Dos.Kgrid` | 半导体约为 `scf.Kgrid` 的 1.5–2 倍；金属更密 | 加密后 DOS 峰形、带边不再变化 |
| `Dos.Erange` | 下限低于所有关心的占据态 | 0 eV 处积分值等于窗口内电子数 |
| 四面体 / 高斯 | 周期体系用四面体；分子、O(N)、大体系用高斯 | 高斯宽度 0.05–0.2 eV，比较两种方法的峰位 |
| 高斯宽度 $a$ | 0.05–0.1 eV（固体），0.1–0.3 eV（分子） | 太大抹平带边，太小出现锯齿 |
| `cohp` 能量网格数 | 2000–5000 | 与 DOS 分辨率一致即可 |
| `CDDF.Kgrid` | Si 类半导体 ≥ 50×50×50；金属更密 | 不同网格下 ε₂ 峰位与峰高收敛 |
| `CDDF.FWHM` | 0.1–0.2 eV | 与实验分辨率或展宽对应 |
| `CDDF.maximum_energy` + `additional` | 覆盖关心范围 + 1 eV | 上端 1 eV 内的数据受截断影响，舍弃 |
| BoltzTraP2 `-m` | 10–25 | 增大到结果不变 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `error while loading shared libraries: libmpifort.so.12` | 官方 deb 包的运行库不在搜索路径 | `export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH`，直接调用 `DosMain` |
| 没有 `.Dos.val`/`.Dos.vec` | `Dos.fileout` 未打开，或 SCF 中途失败 | 检查输入与 `.std` 末尾 |
| `Dos.fileout and DosGauss.fileout are mutually exclusive.` | 两者同时为 `on` | 只开其一 |
| `DosGauss.fileout is supported for only band calculation.` | Cluster 等求解器下使用了 `DosGauss.fileout` | 改用 `Dos.fileout` |
| `DosMain` 因内存不足崩溃 | `.Dos.vec` 太大 | 减小 `Dos.Erange`/`Dos.Kgrid`，或改用 `DosGauss.fileout` |
| 0 eV 处积分 DOS 不等于价电子数 | 从 `Dos.Erange` 下限开始积分，半芯态在窗口外；高斯尾巴 | 下移下限；用 Mulliken 布居求占据数（§6） |
| 自旋极化 PDOS 的下自旋是负数 | 4.0 `DosMain` 写出 −DOS↓ | 作图时直接使用；积分、求和时取绝对值 |
| 找不到 `d1`… 文件 | 该原子的基组没有 d 轨道（如 `H6.0-s2p1`），或该原子未被选入 PDOS | 检查 PAO 定义与 `DosMain` 中输入的原子编号 |
| 非交互运行 `DosMain` 卡住或结果错乱 | 答案顺序与提问不符（Cluster/DC 时不问方法；单原子不问原子编号） | 用 `run_dosmain.sh`，或先交互运行一次确认提问顺序 |
| `In case of the COHP calculation, Dos.fileout should be switched on.` | 开了 `COHP.calc.flag` 没开 `Dos.fileout` | 打开 `Dos.fileout` |
| `Format error for COHP.pairs` | 行数不等于 `COHP.num.pairs`，或每行不是 6 个整数 | 核对格式 `序号 原子1 原子2 l1 l2 l3` |
| 等价键的 ICOHP 不相等 | 晶胞平移 `l1 l2 l3` 写错，选到了别的原子 | 用距离核对（可在 VESTA/OpenMX Viewer 中量） |
| `You want to 'Cluster' calculation for SCF, but 'Band' calculation for CDDF calculation.` | `scf.Kgrid 1 1 1` 而 `CDDF.Kgrid` 更密 | SCF 也用 Band 与多个 k 点，或让二者一致 |
| nc/SOC + CDDF 段错误、`Parameter 10 was incorrect on entry to ZHEEVX` | 4.0 的 CDDF 仅支持共线 | 用共线计算；或参考 klsk007/openmx-soc-cddf-fix |
| 光谱锯齿、峰位随网格跳动 | `CDDF.Kgrid` 不够密 | 加密到收敛（手册测试到 100×100×100） |
| BoltzTraP2 报 `np.infty` 相关错误 | NumPy 2.x 已删除 `np.infty` | 在加载器 `io.py` 中替换为 `np.inf` |
| 费米面文件数值不对 | patch 4.0.1 起 bxsf/frmsf 单位改为 eV/Å | 更新后处理脚本的单位换算 |

---

## 练习

1. 对金刚石分别用 `Dos.Kgrid` = 8、16、24 和四面体法计算 DOS，观察带边与价带中范霍夫奇点的收敛；再用高斯展宽 a = 0.05、0.2 eV 对比。
2. 把 `Cdia_dos.dat` 的 `Dos.Erange` 下限改为 −15 eV，比较 0 eV 处的积分 DOS，解释差别。
3. 对 MoS₂ 画出 Mo 的 d3z²−r²、(dx²−y² + dxy)、(dxz + dyz) 三组 PDOS，判断价带顶和导带底的主要轨道成分，并与第 8 章能带的 K、Γ 点对应。
4. 在 NiO 上用 `sum_pdos.py --by atom --resolve-m` 分别画 Ni1 与 Ni2 的 $e_g$/$t_{2g}$ 自旋分辨 PDOS，验证两者镜像对称；再用第 13 章的 +U 输入比较带隙变化。
5. 给 Si 的 COHP 再加一对第三近邻，比较 ICOHP 随距离的衰减；用 `CDDF.Kgrid` = 20、30、50 计算 ε₂，记录峰位与峰高的变化。

---

## 参考资料

**官方文档**

- [OpenMX 4.0 手册 第 23 章 Density of states](https://www.openmx-square.org/openmx_man4.0/s23_dos.html)：[23.1 标准方法与 DosMain](https://www.openmx-square.org/openmx_man4.0/s23_1_standard-method.html)、[23.2 DosGauss](https://www.openmx-square.org/openmx_man4.0/s23_2_many-k.html)
- [手册 第 24 章 COHP and COOP](https://www.openmx-square.org/openmx_man4.0/s24_cohp-coop.html)：[24.1 定义与符号](https://www.openmx-square.org/openmx_man4.0/s24_1_overview.html)、[24.2 计算流程](https://www.openmx-square.org/openmx_man4.0/s24_2_procedure.html)
- [手册 第 62 章 Optical conductivity and dielectric function](https://www.openmx-square.org/openmx_man4.0/s62_optical-diel.html)：[62.2 Si 示例](https://www.openmx-square.org/openmx_man4.0/s62_2_si.html)、[62.3 关键词](https://www.openmx-square.org/openmx_man4.0/s62_3_keywords.html)、[62.4 基准测试](https://www.openmx-square.org/openmx_man4.0/s62_4_benchmark.html)、[62.7 自动测试](https://www.openmx-square.org/openmx_man4.0/s62_7_auto-test.html)
- [手册 第 63 章 BoltzTraP/BoltzTraP2 接口](https://www.openmx-square.org/openmx_man4.0/s63_boltztrap.html)：[63.1](https://www.openmx-square.org/openmx_man4.0/s63_1_boltztrap.html)、[63.2](https://www.openmx-square.org/openmx_man4.0/s63_2_boltztrap2.html)
- [手册 8.2 关键词（Dos.fileout、Dos.Erange、Dos.Kgrid）](https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html)、[第 35 章 电荷分析](https://www.openmx-square.org/openmx_man4.0/s35_charge-analysis.html)
- [OpenMX What's new（4.0：COHP/COOP、BoltzTraP2、CDDF.shift_ChemP）](https://www.openmx-square.org/whatsnew/whatsnew.html)
- [Installation via a Debian package（deb 包运行库路径）](https://www.openmx-square.org/openmx_man4.0/s4_11_deb-install.html)
- [BoltzTraP2 官方仓库与文档（GitLab）](https://gitlab.com/sousaw/BoltzTraP2)

**相关 GitHub 仓库（仅链接，未复制其代码）**

- [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS) — 浏览器中同时显示能带与 DOS/PDOS，含 PDOS 积分与出版级样式
- [pranabdas/openmx](https://github.com/pranabdas/openmx) — 英文教程站：Si/金刚石 DOS、轨道投影、NiO 示例
- [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix) — 记录 4.0 中 nc + SOC + CDDF 报错/段错误的分析与补丁
- [samlipton/pyfesm](https://github.com/samlipton/pyfesm) — 轻量 Python 解析器，可参考其输出文件解析思路
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl) — 由 OpenMX 哈密顿量构建紧束缚模型，计算光学与非线性响应
- [hikaruri/OMXsigmaxy](https://github.com/hikaruri/OMXsigmaxy) — 由 `.scfout` 计算反常霍尔电导
- [TnakaTnakaTnaka/TEP_tools](https://github.com/TnakaTnakaTnaka/TEP_tools) — 由透射谱计算 Seebeck、ZT 的思路参考
- 本仓库：../scripts/omx_dos_auto.py、../scripts/omx_parse_out.py、第5章 后处理与可视化、[相关仓库总表](../../Resources/github-repos.md)

**关键文献**

- G. Lehmann, M. Taut, *Phys. Status Solidi B* 54, 469 (1972)（线性四面体法）. DOI: [10.1002/pssb.2220540211](https://doi.org/10.1002/pssb.2220540211)；P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* 49, 16223 (1994). DOI: [10.1103/PhysRevB.49.16223](https://doi.org/10.1103/PhysRevB.49.16223)
- R. S. Mulliken, *J. Chem. Phys.* 23, 1833 (1955). DOI: [10.1063/1.1740588](https://doi.org/10.1063/1.1740588)
- R. Dronskowski, P. E. Blöchl, *J. Phys. Chem.* 97, 8617 (1993)（COHP）. DOI: [10.1021/j100135a014](https://doi.org/10.1021/j100135a014)；T. Hughbanks, R. Hoffmann, *J. Am. Chem. Soc.* 105, 3528 (1983)（COOP）. DOI: [10.1021/ja00349a027](https://doi.org/10.1021/ja00349a027)
- C.-C. Lee, Y.-T. Lee, M. Fukuda, T. Ozaki, *Phys. Rev. B* 98, 115115 (2018)；Y.-T. Lee, C.-C. Lee, M. Fukuda, T. Ozaki, *Phys. Rev. B* 102, 075143 (2020)（OpenMX 中光学电导率的实现，源码注释引用）. DOI: [10.1103/PhysRevB.98.115115](https://doi.org/10.1103/PhysRevB.98.115115)、[10.1103/PhysRevB.102.075143](https://doi.org/10.1103/PhysRevB.102.075143)
- D. E. Aspnes, A. A. Studna, *Phys. Rev. B* 27, 985 (1983)（Si 等半导体的实验介电函数）. DOI: [10.1103/PhysRevB.27.985](https://doi.org/10.1103/PhysRevB.27.985)
- K. F. Mak, C. Lee, J. Hone, J. Shan, T. F. Heinz, *Phys. Rev. Lett.* 105, 136805 (2010)（单层 MoS₂ 直接带隙）. DOI: [10.1103/PhysRevLett.105.136805](https://doi.org/10.1103/PhysRevLett.105.136805)
- G. A. Sawatzky, J. W. Allen, *Phys. Rev. Lett.* 53, 2339 (1984)（NiO 带隙）. DOI: [10.1103/PhysRevLett.53.2339](https://doi.org/10.1103/PhysRevLett.53.2339)
- G. K. H. Madsen, D. J. Singh, *Comput. Phys. Commun.* 175, 67 (2006)（BoltzTraP）. DOI: [10.1016/j.cpc.2006.03.007](https://doi.org/10.1016/j.cpc.2006.03.007)；G. K. H. Madsen, J. Carrete, M. J. Verstraete, *Comput. Phys. Commun.* 231, 140 (2018)（BoltzTraP2）. DOI: [10.1016/j.cpc.2018.05.010](https://doi.org/10.1016/j.cpc.2018.05.010)

---

[← 上一章](08-band.md) | [返回目录](../README.md) | [下一章 →](10-charge-mo-cube.md)
