# 第5章 能带结构与态密度（含 KPOINTS_OPT）（Band Structure & DOS with VASP, incl. KPOINTS_OPT）

本章以金刚石 Si 为例，完整演示 VASP 能带结构（band structure）与态密度（DOS/PDOS）的计算：传统的"自洽 + ICHARG=11 非自洽"两步法、VASP 6.3 起的 KPOINTS_OPT 一步法，高对称 k 路径的生成（seekpath、pymatgen、VASPKIT、ASE、sumo），LORBIT=11 投影态密度与 fat band，BANDGAP 标签直接报告带隙，以及用 py4vasp、sumo、VASPKIT 和自写 matplotlib 脚本作图。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 用两步法（SCF → `ICHARG = 11`）算出正确的 PBE 能带和 DOS，并知道两步之间哪些参数**必须一致**；
- 用 `KPOINTS_OPT`（VASP ≥ 6.3，推荐 ≥ 6.6.1）在一次运行里同时完成自洽和能带路径计算；
- 为任意晶体生成高对称路径，并理解"标准原胞 vs 输入晶胞"的约定问题；
- 设置 `LORBIT = 11`、`NEDOS`、`ISMEAR = -5` 得到平滑的 TDOS/PDOS，读懂 DOSCAR；
- 用 `BANDGAP`（VASP ≥ 6.4.3）、py4vasp、VASPKIT 读取带隙、VBM/CBM 位置；
- 用 py4vasp、sumo、VASPKIT 或自己的 Python 脚本画能带、fat band 和 PDOS。

| 你要算的 | 泛函 | 推荐方法 | 小节 |
|---|---|---|---|
| PBE/LDA/PBE+U 能带 | 半局域 | 两步法 或 KPOINTS_OPT | §5、§6 |
| HSE06 / PBE0 能带 | 杂化 | 零权重 k 点 或 KPOINTS_OPT（**不能用 ICHARG=11**） | [第6章](06-hybrid-metagga.md) |
| SCAN / r2SCAN / mBJ 能带 | meta-GGA | 同杂化泛函（需要规则网格计算动能密度） | [第6章](06-hybrid-metagga.md) |
| DOS / PDOS / d 带中心 | 任意 | 密网格 + `ISMEAR = -5` + `LORBIT = 11` | §7 |
| SOC 能带、自旋纹理 | — | `vasp_ncl`，见 [第7章](07-magnetism-u-soc.md) | — |
| 超胞能带反折叠、有效质量、费米面 | — | 见 [第17章](17-advanced-electronic.md) | — |

---

## 2. 原理速览

### 2.1 能带与态密度

周期体系的 Kohn–Sham 本征值 $\varepsilon_{n\mathbf k}$ 依赖于布里渊区中的波矢 $\mathbf k$ 与带指标 $n$。**能带图**是沿一条连接高对称点的路径画出的 $\varepsilon_{n\mathbf k}$；**态密度**则是对整个布里渊区积分：

$$
g(E) = \sum_n \int_{\mathrm{BZ}} \frac{d^3k}{\Omega_{\mathrm{BZ}}}\, \delta(E - \varepsilon_{n\mathbf k}),
$$

投影态密度（PDOS）在 δ 函数前乘上波函数投影到原子 $I$ 的 $lm$ 轨道上的权重 $|\langle Y^{I}_{lm}|\psi_{n\mathbf k}\rangle|^2$。实际计算中 δ 函数用 Gaussian/Methfessel–Paxton 展宽，或者用四面体方法（`ISMEAR = -5`，Blöchl 修正）做积分；后者对 DOS 来说最平滑、无需选择 SIGMA。

### 2.2 为什么能带要分两步

对 LDA/GGA 而言，哈密顿量完全由电荷密度 $\rho(\mathbf r)$ 决定。电荷密度需要在**均匀网格**上对整个布里渊区积分才能得到；而能带路径上的 k 点既不均匀、也没有正确的权重，**不能**用来做自洽。所以标准流程是：

1. 在均匀 k 网格上自洽，得到并保存 `CHGCAR`；
2. 读入 `CHGCAR` 并**固定**电荷密度（`ICHARG = 11`），在路径 k 点上只对角化哈密顿量。

`KPOINTS_OPT` 本质上把这两步合到同一次运行中：自洽完成后，VASP 在 `KPOINTS_OPT` 指定的 k 点上再做一次非自洽对角化。

> [!IMPORTANT]
> 杂化泛函和 meta-GGA 的哈密顿量不只依赖电荷密度（还依赖规则网格上的轨道或动能密度），因此**绝对不能**对它们用 `ICHARG = 11`。它们的能带见 [第6章](06-hybrid-metagga.md)。

### 2.3 费米能级与带隙

- 金属：费米能级 $E_F$ 由占据数决定，必须取自**均匀网格**的计算（SCF 或密网格 DOS），能带路径计算给出的 E-fermi 不可靠。
- 半导体：$E_F$ 在带隙中的位置并不唯一，VASP 默认放在带隙中偏下的某个位置；从 VASP 6.4 起可以设 `EFERMI = MIDGAP` 让它落在带隙中间。作图时半导体通常把 **VBM 设为 0**。
- 间接带隙 $E_g = \min_{\mathbf k}\varepsilon_{c}(\mathbf k) - \max_{\mathbf k}\varepsilon_{v}(\mathbf k)$；直接带隙是同一 k 点上导带与价带的最小间距。Si 的 CBM 位于 Γ–X 线上约 85% 处，**不在**常用的均匀网格点上，所以只看 SCF 网格得到的带隙会偏大。

---

## 3. 工作流总览

```text
            POSCAR (relaxed, Chapter 4)
                     |
          1-scf: uniform Gamma mesh, LCHARG=.TRUE.  --> CHGCAR, E-fermi, BANDGAP
             /                     \
 2-band: ICHARG=11 +            3-dos: ICHARG=11 + dense mesh
   line-mode KPOINTS              ISMEAR=-5, LORBIT=11, NEDOS
   --> EIGENVAL/PROCAR/            --> DOSCAR, vasprun.xml
       vasprun.xml/vaspout.h5

 or, in ONE run (VASP >= 6.3):
 4-kpoints-opt: KPOINTS (uniform) + KPOINTS_OPT (line mode)
   --> SCF results + band path in vaspout.h5 (kpoints_opt), vasprun.xml, PROCAR_OPT
```

示例输入全部在 [../examples/05-band-dos-si/README.md](../examples/05-band-dos-si/README.md)，按目录 `1-scf`、`2-band-icharg11`、`3-dos`、`4-band-kpoints-opt` 组织，`job.slurm` 会按顺序全部跑完。

---

## 4. 第一步：自洽计算

### 4.1 输入

POSCAR 使用第4章优化得到的晶格（示例中取 a = 5.47 Å，请替换为你自己的 CONTCAR）：

```text
Si diamond primitive cell (use your relaxed lattice constant)
5.47
  0.0  0.5  0.5
  0.5  0.0  0.5
  0.5  0.5  0.0
Si
2
Direct
  0.00  0.00  0.00
  0.25  0.25  0.25
```

KPOINTS：

```text
Gamma-centered 8x8x8 for SCF
0
Gamma
8 8 8
0 0 0
```

INCAR：

```text
SYSTEM  = Si SCF (charge density for band and DOS)

PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.

NSW     = 0
IBRION  = -1

LCHARG  = .TRUE.     # CHGCAR is needed by the ICHARG=11 runs
LWAVE   = .FALSE.
```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `ENCUT`、`PREC` | 与后续步骤**完全相同** | 两步法的硬性要求 |
| `ISMEAR`/`SIGMA` | 0 / 0.01 | 半导体；金属用 `ISMEAR = 1`、`SIGMA = 0.1~0.2`，或 `ISMEAR = -5` 得到更准的 E-fermi |
| `LCHARG` | .TRUE. | 必须写出 CHGCAR |
| `EDIFF` | 1E-7 | 电荷密度越精确，后续非自洽能带越可靠 |

可选标签：

- `BANDGAP = COMPACT | WEIGHT | KPOINT`（VASP ≥ 6.4.3）：控制 OUTCAR 中带隙报告的详细程度，默认 `COMPACT`；
- `EFERMI = MIDGAP`（VASP ≥ 6.4）：半导体的费米能放在带隙中间，官方推荐。

### 4.2 运行与检查

```bash
cd VASP/examples/05-band-dos-si/1-scf
cp ../POTCAR .            # POTCAR built once in the parent directory
mpirun -np 8 vasp_std > vasp.out

grep "E-fermi" OUTCAR | tail -1      # Fermi level, needed later for plotting
ls -l CHGCAR                         # must exist and be non-empty
grep -A12 "Band structure" OUTCAR    # band-gap report written by VASP >= 6.4.3
```

BANDGAP 报告的结构如下（示意输出，数值仅供参考）：

```text
 Band structure
 --------------
                        spin independent
 val. band max:               5.62
 cond. band min:              6.30
 fundamental gap:             0.68
 VBM @ kpoint:       0.0000   0.0000   0.0000
 CBM @ kpoint:       0.3750   0.0000   0.3750
 ...
 direct gap:                  2.55
```

> [!NOTE]
> SCF 网格上报告的"fundamental gap"只是在网格点上找到的最大/最小值。Si 的 CBM 不在 8×8×8 网格点上，所以这里的带隙会略大于真实值；以能带路径或更密网格的结果为准。

---

## 5. 高对称 k 路径

### 5.1 约定问题：先选好"晶胞"

高对称点的分数坐标是**相对于某个特定原胞的倒格矢**写的。不同工具使用不同约定：

| 工具 | 约定 | 会不会改变你的晶胞 |
|---|---|---|
| Setyawan–Curtarolo（pymatgen `HighSymmKpath` 默认、AFLOW） | 标准原胞 | 要求输入就是它的标准原胞；否则用 `kpath.prim` 作为 POSCAR |
| Hinuma 等（SeeK-path、VASPKIT 303 的体相路径） | 晶体学标准原胞 | 会输出新的原胞（VASPKIT 写出 `PRIMCELL.vasp`），可能旋转/重排晶格矢量 |
| Latimer–Munro（pymatgen `path_type="latimer_munro"`） | 适用于任意输入晶胞 | 不需要换晶胞 |
| ASE `cell.bandpath()` | 自动把特殊点变换到**你的**晶胞 | 不需要换晶胞 |
| sumo-kgen | 默认 Bradley–Cracknell，可选 `--seekpath`/`--pymatgen` | 可能写出标准化的 POSCAR，按其提示使用 |

> [!WARNING]
> 最常见的错误：用工具 A 的高对称点坐标，却用工具 B（或自己手搭）的晶胞做计算。结果是能带"看起来怪怪的"——带隙位置不对、本该简并的地方劈裂。**规则：生成路径的工具给出了新原胞，就用那个原胞重新做 SCF。** 至少检查一次：把路径 k 点换算成笛卡尔坐标，看高对称点是否落在布里渊区表面的正确位置。

### 5.2 Si（fcc）的路径

本章 POSCAR 的晶格矢量 $(0,\tfrac12,\tfrac12)a$、$(\tfrac12,0,\tfrac12)a$、$(\tfrac12,\tfrac12,0)a$ 恰好是 Setyawan–Curtarolo 的 fcc 原胞约定，高对称点（倒格矢分数坐标）为：

| 点 | 分数坐标 | 笛卡尔坐标（单位 $2\pi/a$） |
|---|---|---|
| Γ | (0, 0, 0) | (0, 0, 0) |
| X | (0.5, 0, 0.5) | (0, 1, 0) |
| L | (0.5, 0.5, 0.5) | (0.5, 0.5, 0.5) |
| W | (0.5, 0.25, 0.75) | (0.5, 1, 0) |
| K | (0.375, 0.375, 0.75) | (0.75, 0.75, 0) |
| U | (0.625, 0.25, 0.625) | (0.25, 1, 0.25) |

金刚石结构常用路径 L–Γ–X–U|K–Γ，Line-mode KPOINTS 如下（标签必须是 ASCII，用 `GAMMA` 表示 Γ）：

```text
Si band path L-G-X-U|K-G
40              ! points per segment (end points included)
Line-mode
Reciprocal
  0.500  0.500  0.500   L
  0.000  0.000  0.000   GAMMA

  0.000  0.000  0.000   GAMMA
  0.500  0.000  0.500   X

  0.500  0.000  0.500   X
  0.625  0.250  0.625   U

  0.375  0.375  0.750   K
  0.000  0.000  0.000   GAMMA
```

格式要点：第 2 行是每段的点数（包含两个端点）；第 3 行首字母 `L` 表示 line mode；第 4 行 `Reciprocal`（或 Fractional，即非 `C/K` 开头）表示分数坐标；之后每两行定义一段，空行只是为了可读。U 与 K 是等价点，"U|K" 处的跳跃在作图时不增加横坐标距离。

### 5.3 自动生成路径的几种方法

**SeeK-path**（会给出它的标准原胞，务必一起使用）：

```python
import seekpath
from ase.io import read

atoms = read("POSCAR")
res = seekpath.get_path((atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers), with_time_reversal=True)
print("Bravais lattice:", res["bravais_lattice_extended"])
print("path:", res["path"])                       # list of (start, end) label pairs
for k, v in res["point_coords"].items():
    print("%-8s %8.4f %8.4f %8.4f" % (k, *v))
print("primitive lattice (use this cell!):\n", res["primitive_lattice"])
```

**pymatgen**：

```python
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp import Kpoints

s = Structure.from_file("POSCAR")
kpath = HighSymmKpath(s, path_type="setyawan_curtarolo")   # or "hinuma", "latimer_munro"
kpath.prim.to(filename="POSCAR_prim")                      # standard primitive cell for this path
Kpoints.automatic_linemode(40, kpath).write_file("KPOINTS_line")
```

**VASPKIT**（1.5.x，命令行模式）：

```bash
vaspkit -task 303      # 3D bulk: writes KPATH.in and PRIMCELL.vasp (and HIGH_SYMMETRY_POINTS)
cp PRIMCELL.vasp POSCAR   # if the cell changed, redo the SCF with this POSCAR
cp KPATH.in KPOINTS       # for the ICHARG=11 run (or KPOINTS_OPT)
# 2D materials: vaspkit -task 302 ; 1D: 301
```

**ASE**（直接对你的晶胞给出路径）：

```python
from ase.io import read
atoms = read("POSCAR")
bp = atoms.cell.bandpath("LGXU,KG", npoints=160)    # ',' marks a discontinuity
print(bp.special_points)                            # coordinates in the basis of YOUR reciprocal cell
```

**sumo**：在含 POSCAR 的目录运行 `sumo-kgen`（默认 Bradley–Cracknell 路径，`--seekpath`/`--pymatgen` 切换约定，`--density` 控制点密度，默认 60）。

**本仓库脚本**：`../scripts/vasp_inputgen.py` 可以为 band 任务生成 Line-mode KPOINTS 和 KPOINTS_OPT，参数用 `python ../scripts/vasp_inputgen.py -h` 查看。

### 5.4 二维材料的路径

二维六角晶格（如 MoS2、石墨烯）常用 Γ–M–K–Γ，z 分量为 0。K 点的分数坐标取决于两个面内晶格矢量的夹角：

| 面内夹角 | M | K |
|---|---|---|
| 120°（如第4章 MoS2 的 POSCAR） | (0.5, 0, 0) | (1/3, 1/3, 0) |
| 60° | (0.5, 0, 0) | (1/3, 2/3, 0) 或 (2/3, 1/3, 0) |

VASPKIT 302 会按二维约定生成路径；也可以用 ASE `bandpath` 让程序自己判断。

---

## 6. 方法 A：传统两步法（ICHARG = 11）

### 6.1 INCAR

把 SCF 的 INCAR 复制过来，只改动下面标注的部分：

```text
SYSTEM  = Si band structure, non-self-consistent (ICHARG=11)

PREC    = Accurate
ENCUT   = 400        # MUST be identical to the SCF run
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.

ICHARG  = 11         # read CHGCAR and keep the density fixed
NBANDS  = 16         # 4 occupied bands for Si; add empty bands to see the conduction band
LORBIT  = 11         # lm-decomposed PROCAR (fat bands)

NSW     = 0
IBRION  = -1
LCHARG  = .FALSE.    # do not overwrite the SCF CHGCAR
LWAVE   = .FALSE.
```

KPOINTS 使用 §5.2 的 Line-mode 文件。

### 6.2 两步之间必须一致的参数

| 必须一致 | 原因 |
|---|---|
| POSCAR（晶格、原子顺序）、POTCAR | CHGCAR 只对同一结构、同一赝势有意义 |
| `ENCUT`、`PREC`（FFT 网格 NGX/NGY/NGZ） | 网格不同，读 CHGCAR 会插值或报错 |
| `GGA`/`METAGGA`、`LASPH`、vdW 设置 | 改变哈密顿量 |
| `ISPIN`、`MAGMOM`（磁性体系）、`LNONCOLLINEAR`/`LSORBIT` | 自旋密度也在 CHGCAR 中 |
| `LDAU*` 全部标签、**`LMAXMIX`** | DFT+U 的 on-site 占据矩阵只写到 `LMAXMIX` 阶；d 元素必须 `LMAXMIX = 4`、f 元素 `6`，**SCF 时就要设好** |

可以（也应该）改动的：`ICHARG`、`KPOINTS`、`NBANDS`、`LORBIT`、`ISMEAR`/`SIGMA`（对能带本征值无影响）、输出控制标签。

### 6.3 运行与检查

```bash
cd ../2-band-icharg11
cp ../1-scf/{POSCAR,POTCAR,CHGCAR} .
mpirun -np 8 vasp_std > vasp.out

grep "NKPTS" OUTCAR | head -1      # 160 k-points = 4 segments x 40
ls EIGENVAL PROCAR vasprun.xml
```

一个很实用的一致性检查：路径上的 Γ 点本征值应与 SCF 中 Γ 点的本征值**完全相同**（差别 < 1 meV）。

```bash
# SCF, first k-point (Gamma), bands 1-6
grep -A7 "k-point     1 :" ../1-scf/OUTCAR | tail -6
# band run: the 40th k-point is GAMMA (end of segment L-GAMMA)
grep -A7 "k-point    40 :" OUTCAR | tail -6
```

如果两者差几十 meV，说明两步之间某个参数不一致（最常见：`ENCUT`、`LMAXMIX`、`MAGMOM`），或 CHGCAR 不是这个结构的。

> [!WARNING]
> `ICHARG = 11` 这一步的 OUTCAR 里也有 E-fermi，但它是由路径 k 点算出来的，**不可靠**（金属尤其明显）。作图时请使用 SCF（或密网格 DOS 计算）的 E-fermi。VASP 也允许在 `ICHARG = 11` 计算中用 `EFERMI = <数值>` 直接固定费米能。

---

## 7. 方法 B：KPOINTS_OPT 一步法（VASP ≥ 6.3）

### 7.1 原理与文件

`KPOINTS_OPT` 是一个可选输入文件，格式与 KPOINTS 完全相同。VASP 先用 `KPOINTS` 中的 k 点完成自洽，然后**在同一次运行中**对 `KPOINTS_OPT` 中的 k 点做一次非自洽计算。

- `KPOINTS` 必须是**均匀网格**（自洽用）；`KPOINTS_OPT` 通常写 Line-mode 路径，也可以写更密的均匀网格（用于 DOS）。
- 目录中存在 `KPOINTS_OPT` 就会自动读取；不想用时设 `LKPOINTS_OPT = .FALSE.`，不必删文件。
- 结果写入：`vaspout.h5` 中带 `kpoints_opt` 标记的数据组；`vasprun.xml` 中的 `<eigenvalues_kpoints_opt>`、`<projected_kpoints_opt>`（以及 `comment="kpoints_opt"` 的 DOS 段）；`LORBIT ≥ 10` 时写 `PROCAR_OPT`（格式同 PROCAR）。请从这些文件读取路径能带，不要从 EIGENVAL 读。
- 内存不足时用 `KPOINTS_OPT_NKBATCH` 减小每批同时处理的 k 点数（默认等于自洽的不可约 k 点数）。
- VASP ≥ 6.6.0 可以不写 `KPOINTS_OPT` 文件，而用 `KSPACING_OPT`（Å⁻¹，默认 0.5）自动生成一个更密的均匀网格，用于得到比自洽网格更细的 DOS。

### 7.2 版本说明

| 版本 | 与 KPOINTS_OPT 相关的变化 |
|---|---|
| 6.3.0 | 引入 KPOINTS_OPT |
| 6.4.3 | 用 `vasp_gam` 运行时给出明确的报错（`vasp_gam` 不支持 KPOINTS_OPT，请用 `vasp_std`） |
| 6.5.0 | 修复：`NCORE ≠ 1` 时结果不正确、`LORBIT = 10/11/12` 时崩溃、非共线计算崩溃、vasprun.xml 投影错误；新增 `KPOINTS_OPT_MODE` |
| 6.5.0/6.5.1 | `PROCAR_OPT` 损坏（known issue 45，6.6.0 修复） |
| 6.6.0 | 新增 `KSPACING_OPT`；`ISPIN = 2` + `LORBIT ≥ 11` 时 vaspout.h5 缺第二自旋通道的 PDOS（known issue 69） |
| 6.6.1 | 修复 issue 69 |

> [!TIP]
> 结论：非自旋极化体系 6.5.0 起即可放心使用；需要 `PROCAR_OPT` 时用 ≥ 6.6.0；自旋极化 + 投影时用 **6.6.1**。

### 7.3 输入

INCAR（与 SCF 相同，加上投影和空带）：

```text
SYSTEM  = Si SCF + band path in one run (KPOINTS_OPT, VASP >= 6.3)

PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.

NBANDS  = 16
LORBIT  = 11         # PROCAR_OPT + projections in vaspout.h5

NSW     = 0
IBRION  = -1
LCHARG  = .TRUE.
LWAVE   = .FALSE.
```

KPOINTS：与 §4 相同的 Γ 中心 8×8×8。KPOINTS_OPT：与 §5.2 的 Line-mode 文件内容相同。

```bash
cd ../4-band-kpoints-opt
cp ../POTCAR ../1-scf/POSCAR .
mpirun -np 8 vasp_std > vasp.out       # vasp_std, NOT vasp_gam
ls PROCAR_OPT vaspout.h5
```

### 7.4 两种方法怎么选

| | ICHARG = 11 两步法 | KPOINTS_OPT |
|---|---|---|
| 版本 | 所有版本 | ≥ 6.3（建议 ≥ 6.5.0 / 6.6.1） |
| 作业数 | 2（SCF + band） | 1 |
| 误操作风险 | 两步参数不一致、CHGCAR 拷错 | 很小 |
| 输出 | EIGENVAL、PROCAR、vasprun.xml | vaspout.h5、vasprun.xml、PROCAR_OPT |
| 杂化 / meta-GGA | 不可用 | 可用（第6章） |
| 后处理工具兼容性 | 几乎所有工具 | py4vasp、pymatgen（自动识别）；老工具可能只认 EIGENVAL |

---

## 8. 态密度 DOS / PDOS

### 8.1 INCAR 与 KPOINTS

DOS 需要对整个布里渊区积分，k 网格要比 SCF **更密**。推荐在 SCF 的 CHGCAR 上做一次 `ICHARG = 11` 的密网格计算：

```text
SYSTEM  = Si DOS on a dense mesh (ICHARG=11)

PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
LREAL   = .FALSE.

ICHARG  = 11
ISMEAR  = -5         # tetrahedron method with Bloechl corrections (Gamma-centered mesh, >= 4 k-points)
LORBIT  = 11         # site- and lm-projected DOS in DOSCAR
NEDOS   = 3001       # energy grid points between EMIN and EMAX
EMIN    = -15
EMAX    = 10
NBANDS  = 16

NSW     = 0
IBRION  = -1
LCHARG  = .FALSE.
LWAVE   = .FALSE.
```

```text
Gamma-centered 16x16x16 for DOS
0
Gamma
16 16 16
0 0 0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ISMEAR` | −5 | 最平滑的 DOS；要求 Γ 中心网格且不可约 k 点 ≥ 4。分子、slab 的 Γ 点计算不能用 |
| `LORBIT` | 11 | PAW 投影，不需要 `RWIGS`；DOSCAR 中给出每个原子的 s、p_y、p_z、p_x、d_xy… 分量 |
| `NEDOS` | 2001–3001 | 默认只有 301，峰会很粗糙；能量间隔 = (EMAX−EMIN)/NEDOS |
| `EMIN`/`EMAX` | 覆盖感兴趣的能区 | 绝对能量（不是相对 E-fermi），先看 SCF 的 E-fermi 再定 |
| `NBANDS` | 足够覆盖 EMAX | 空带太少，导带部分的 DOS 会被截断 |

> [!TIP]
> 金属体系的 DOS 同样推荐 `ISMEAR = -5`（只算 DOS 时没问题），但**金属的结构优化不要用 −5**（受力不准）。如果担心四面体方法在某些体系（大超胞、k 点很少）出问题，可退回 `ISMEAR = 0`、`SIGMA = 0.05` 并适当加密网格。

### 8.2 DOSCAR 格式

```text
line 1      : NIONS NIONS 1(partial DOS present) NCDIJ
line 2-5    : volume, lattice lengths, POTIM / TEBEG / 'CAR' / SYSTEM
line 6      : EMAX  EMIN  NEDOS  E-fermi  1.0
next NEDOS  : energy  DOS  integrated-DOS                     (ISPIN=1)
              energy  DOS(up) DOS(dn)  intDOS(up) intDOS(dn)  (ISPIN=2)
then for each atom: one header line like line 6, followed by NEDOS lines
              LORBIT=11, ISPIN=1 : energy s py pz px dxy dyz dz2 dxz dx2-y2 (+ f columns if present)
              LORBIT=11, ISPIN=2 : every orbital has an (up, down) pair of columns
```

注意：DOSCAR 中的能量**没有**减去费米能；单位是 states/eV（每个晶胞）；PDOS 是 states/atom/eV。PDOS 之和一般小于 TDOS，因为 PAW 球之间的间隙区电子没有被投影到任何原子上。

### 8.3 检查

```bash
grep "E-fermi" OUTCAR | tail -1
sed -n 6p DOSCAR                      # EMAX EMIN NEDOS E-fermi
# integrated DOS at E-fermi should equal NELECT (8 for Si)
awk 'NR==6{ef=$4} NR>6 && NR<=6+3001 && $1<=ef {n=$3} END{print "integrated DOS at E_F =", n}' DOSCAR
grep NELECT OUTCAR
```

### 8.4 用 Python 画 TDOS 与 Si 的 s/p PDOS

```python
#!/usr/bin/env python3
"""Plot total DOS and s/p PDOS of all atoms from DOSCAR (ISPIN=1, LORBIT=11)."""
import numpy as np
import matplotlib.pyplot as plt

with open("DOSCAR") as fh:
    lines = fh.readlines()
nions = int(lines[0].split()[0])
emax, emin, nedos, efermi = [float(x) for x in lines[5].split()[:4]]
nedos = int(nedos)

tot = np.array([l.split() for l in lines[6:6 + nedos]], float)
energy, tdos = tot[:, 0] - efermi, tot[:, 1]

s = np.zeros(nedos)
p = np.zeros(nedos)
start = 6 + nedos
for ion in range(nions):
    blk = np.array([l.split() for l in lines[start + 1:start + 1 + nedos]], float)
    s += blk[:, 1]
    p += blk[:, 2:5].sum(axis=1)          # py + pz + px
    start += nedos + 1

fig, ax = plt.subplots(figsize=(5, 3.2))
ax.fill_between(energy, tdos, color="0.85", label="total")
ax.plot(energy, s, lw=1.2, label="s (all atoms)")
ax.plot(energy, p, lw=1.2, label="p (all atoms)")
ax.axvline(0, ls="--", lw=0.8, color="k")
ax.set_xlim(-13, 6)
ax.set_ylim(0, None)
ax.set_xlabel(r"$E - E_F$ (eV)")
ax.set_ylabel("DOS (states/eV/cell)")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("dos.png", dpi=200)
```

功能更完整的版本（按元素/原子/轨道选择、t2g/eg、自旋上下对称作图、d 带中心）见本仓库 `../scripts/plot_dos.py`，用 `python ../scripts/plot_dos.py -h` 查看用法。

---

## 9. 读取带隙、VBM 与 CBM

| 方法 | 命令 | 备注 |
|---|---|---|
| VASP 自带（≥ 6.4.3） | `grep -A12 "Band structure" OUTCAR` | 由 `BANDGAP` 标签控制详细程度；也写入 vaspout.h5 |
| py4vasp | `calc.bandgap.fundamental()`、`calc.bandgap.direct()` | 读取 vaspout.h5 |
| VASPKIT | `vaspkit -task 911` | 输出带隙类型、VBM/CBM 及其 k 点 |
| pymatgen | `Vasprun("vasprun.xml").get_band_structure(line_mode=True).get_band_gap()` | 返回 energy、direct、transition |
| 本仓库脚本 | `python ../scripts/plot_band.py -h` | 作图同时报告带隙 |

> [!NOTE]
> 用于报告带隙时，请使用**能带路径 + 均匀网格**中的极值（极值可能在路径之外，比如低对称体系）。参考值：Si 实验间接带隙约 1.17 eV（低温），PBE 一般只给出约 0.5–0.6 eV——这是半局域泛函的系统性低估，不是计算错误。要得到更合理的带隙，请看 [第6章](06-hybrid-metagga.md)。

---

## 10. 投影能带（fat band）与 PROCAR

`LORBIT = 11` 时，PROCAR（或 KPOINTS_OPT 的 PROCAR_OPT）给出每个 k 点、每条带在每个原子 s/p/d 轨道上的投影权重：

```text
k-point     1 :    0.50000000 0.50000000 0.50000000     weight = 0.00625000

band     1 # energy   -9.3   # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
    1  0.xxx  0.xxx  ...
    2  0.xxx  0.xxx  ...
tot    0.xxx  0.xxx  ...
```

（示意格式，数值省略）

把权重作为线宽或颜色叠加到能带上就是 fat band。常用工具：

- py4vasp：`calc.band.plot("Si(s), Si(p)")`、`calc.band.plot("1(p)")`（第 1 个原子的 p 轨道），`width` 参数控制最大线宽；
- pyprocar：功能最全的 PROCAR 绘图工具（fat band、自旋纹理、能带反折叠、费米面），见 [pyprocar 文档](https://pyprocar.readthedocs.io/)；
- sumo-bandplot：支持投影能带，选项见 `sumo-bandplot -h`；
- VASPKIT：`212`（选定原子）、`213`（按元素）、`214`（选定原子求和）；
- 本仓库 `../scripts/plot_band.py` 也支持 PROCAR 的 fat band。

---

## 11. 作图：py4vasp、sumo、VASPKIT、pymatgen 与自写脚本

### 11.1 py4vasp（官方 Python 后处理，读取 vaspout.h5）

需要 VASP 编译时打开 HDF5（`-DVASP_HDF5`，见 [第1章](01-setup-compile.md)）。安装：`pip install py4vasp`。在 Jupyter 中：

```python
import py4vasp

# two-step method: band run + Fermi energy from the SCF run
scf  = py4vasp.Calculation.from_path("1-scf")
band = py4vasp.Calculation.from_path("2-band-icharg11")
ef = scf.dos.read()["fermi_energy"]
g = band.band.plot(fermi_energy=ef)      # energies shifted so that E_F = 0
g.yrange = (-13, 6)
g                                         # displays in a notebook; g.show() elsewhere

band.band.plot("Si(s), Si(p)", fermi_energy=ef)   # fat bands (needs LORBIT=11)

# one-step method
opt = py4vasp.Calculation.from_path("4-band-kpoints-opt")
opt.band.plot("kpoints_opt")              # band path stored with the kpoints_opt label
print(opt.bandgap.fundamental(), opt.bandgap.direct())

# DOS
dos = py4vasp.Calculation.from_path("3-dos")
dos.dos.plot("s, p")
dos.dos.to_csv(filename="dos.csv")        # export for Origin/Excel
```

示例目录中的 `plot_with_py4vasp.ipynb` 包含上述全部单元。

### 11.2 sumo（出版级命令行绘图）

```bash
pip install sumo
cd 2-band-icharg11 && sumo-bandplot          # reads vasprun.xml (+ KPOINTS) in the current directory
cd ../3-dos && sumo-dosplot                  # total + projected DOS
sumo-bandplot -h ; sumo-dosplot -h           # projection, energy range, colours, output format ...
```

### 11.3 VASPKIT（交互或命令行）

```bash
cd 2-band-icharg11 && vaspkit -task 211      # BAND.dat, BAND_REFORMATTED.dat, KLABELS, BAND_GAP (+ band.png if matplotlib is configured)
cd ../3-dos && vaspkit -task 111             # TDOS.dat
vaspkit -task 113                            # PDOS for each element: PDOS_Si.dat
```

> [!NOTE]
> VASPKIT 读取 E-fermi 的方式可能与你的预期不同（例如在 band 目录中读该目录的 DOSCAR/OUTCAR）。对两步法，请确认它使用的是 SCF 的 E-fermi，必要时手动平移。

### 11.4 pymatgen

```python
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter

efermi = Vasprun("1-scf/vasprun.xml", parse_potcar_file=False).efermi
vr = Vasprun("2-band-icharg11/vasprun.xml", parse_potcar_file=False)
bs = vr.get_band_structure(kpoints_filename="2-band-icharg11/KPOINTS", line_mode=True, efermi=efermi)
print(bs.get_band_gap())
ax = BSPlotter(bs).get_plot(ylim=(-6, 6))
ax.figure.savefig("band_pmg.png", dpi=200)

# KPOINTS_OPT runs are detected automatically:
bs_opt = Vasprun("4-band-kpoints-opt/vasprun.xml", parse_potcar_file=False).get_band_structure(line_mode=True)
```

### 11.5 不依赖任何专用库：从 EIGENVAL 画能带

```python
#!/usr/bin/env python3
"""Minimal band plot from EIGENVAL + line-mode KPOINTS + POSCAR (ISPIN=1 or 2).
Usage: python band_eigenval.py  E_fermi_from_SCF"""
import sys
import numpy as np
import matplotlib.pyplot as plt

efermi = float(sys.argv[1])

# reciprocal lattice (with 2*pi) from POSCAR
pos = open("POSCAR").read().split("\n")
scale = float(pos[1])
A = np.array([[float(x) for x in pos[i].split()[:3]] for i in (2, 3, 4)]) * scale
B = 2 * np.pi * np.linalg.inv(A).T

# EIGENVAL
ev = open("EIGENVAL").read().split("\n")
ispin = int(ev[0].split()[3])
nelect, nk, nb = [int(float(x)) for x in ev[5].split()[:3]]
kfrac, bands = [], []
i = 7
for _ in range(nk):
    kfrac.append([float(x) for x in ev[i].split()[:3]])
    rows = np.array([ev[i + 1 + j].split() for j in range(nb)], float)
    bands.append(rows[:, 1:1 + ispin])
    i += nb + 2
kfrac, bands = np.array(kfrac), np.array(bands)           # bands: (nk, nb, ispin)

# labels and segments from the line-mode KPOINTS
kp = [l.split() for l in open("KPOINTS").read().split("\n")[4:] if l.strip()]
labels = [l[3] if len(l) > 3 else "" for l in kp]
npts = nk // (len(kp) // 2)

# cumulative path length; no jump between segments (e.g. U|K)
kcart = kfrac @ B
dist = np.zeros(nk)
for k in range(1, nk):
    dist[k] = dist[k - 1] + (0.0 if k % npts == 0 else np.linalg.norm(kcart[k] - kcart[k - 1]))
ticks, ticklabels = [dist[0]], [labels[0]]
for s in range(len(kp) // 2):
    end = dist[(s + 1) * npts - 1]
    lab = labels[2 * s + 1]
    if s + 1 < len(kp) // 2 and labels[2 * s + 2] != lab:
        lab = lab + "|" + labels[2 * s + 2]
    ticks.append(end)
    ticklabels.append(lab)
ticklabels = ["$\\Gamma$" if t.upper() in ("G", "GAMMA") else t for t in ticklabels]

fig, ax = plt.subplots(figsize=(4.5, 3.6))
for sp, color in zip(range(ispin), ("C0", "C3")):          # spin up: blue, spin down: red
    for s in range(len(kp) // 2):                           # draw segment by segment (no line across U|K)
        seg = slice(s * npts, (s + 1) * npts)
        ax.plot(dist[seg], bands[seg, :, sp] - efermi, color=color, lw=1.0)
for t in ticks:
    ax.axvline(t, color="0.6", lw=0.6)
ax.axhline(0, ls="--", lw=0.8, color="k")
ax.set_xticks(ticks)
ax.set_xticklabels(ticklabels)
ax.set_xlim(dist[0], dist[-1])
ax.set_ylim(-13, 6)
ax.set_ylabel(r"$E - E_F$ (eV)")
fig.tight_layout()
fig.savefig("band.png", dpi=200)
```

用法：`python band_eigenval.py $(grep E-fermi ../1-scf/OUTCAR | tail -1 | awk '{print $3}')`。本仓库的 `../scripts/plot_band.py` 是功能完整的版本（vasprun.xml/KPOINTS_OPT、自旋、零权重 k 点剔除、fat band、带隙报告、导出 Origin 可读的 dat），用 `python ../scripts/plot_band.py -h` 查看用法。更多作图规范见 Common 第5章。

---

## 12. 自旋极化、金属与其他体系的注意事项

- **自旋极化**（`ISPIN = 2`）：EIGENVAL、DOSCAR、PROCAR 中都有上下两套数据；SCF 与能带计算的 `MAGMOM` 必须一致，否则可能读入 CHGCAR 后又收敛到别的磁态。自旋极化 + KPOINTS_OPT + PDOS 请用 6.6.1。详见 [第7章](07-magnetism-u-soc.md)。
- **DFT+U**：SCF 时就设 `LMAXMIX = 4`（d）/`6`（f），否则 `ICHARG = 11` 得到的能带与自洽结果不一致。
- **金属**：E-fermi 取自均匀网格（SCF 或密网格 DOS）；能带图上费米能处的能带交叉要用足够密的路径点。
- **二维材料 / slab**：k 网格 z 方向取 1；能带路径 z 分量为 0；真空方向的 DOS 不需要特殊处理，但功函数需要真空能级参考，见 [第9章](09-surface-workfunction.md)。
- **超胞**（缺陷、合金）：能带折叠后难以解读，用能带反折叠（easyunfold、VaspBandUnfolding），见 [第17章](17-advanced-electronic.md)。
- **GPU 版本**：KPOINTS_OPT 在 GPU（OpenACC）上同样可用；DOS 的 `ISMEAR = -5` 与并行无关。

---

## 参数选择与收敛建议

| 参数 | SCF | 能带（ICHARG=11 / KPOINTS_OPT） | DOS |
|---|---|---|---|
| k 点 | 收敛测试值（Si：Γ 8×8×8） | Line-mode，每段 30–60 点（按段长度可变） | SCF 的 1.5–2 倍密度（Si：16×16×16） |
| `ISMEAR`/`SIGMA` | 半导体 0/0.01–0.05；金属 1/0.1–0.2 或 −5 | 0/0.01（不影响本征值） | −5 |
| `NBANDS` | 默认 | 占据带数 + 若干空带（Si：16） | 覆盖 EMAX |
| `LORBIT` | 不需要 | 11（fat band） | 11 |
| `NEDOS` | 默认 | — | 2001–3001 |
| `EDIFF` | 1E-6 ~ 1E-7 | 同 SCF | 同 SCF |
| `LMAXMIX` | d：4；f：6（DFT+U/磁性） | 与 SCF 相同 | 与 SCF 相同 |
| `ICHARG` | 默认（2） | 11（两步法） | 11 |
| 可执行程序 | vasp_std | vasp_std（KPOINTS_OPT 不支持 vasp_gam） | vasp_std |

收敛检查：把 SCF 的 k 网格从 8×8×8 加到 12×12×12，路径能带变化应 < 10 meV；把 DOS 网格加密一级，主要峰位变化应 < 0.05 eV。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `Tetrahedron method fails (number of k-points < 4)` | `ISMEAR = -5` 但不可约 k 点太少（分子、slab Γ 点、很小的网格） | 加密 Γ 中心网格，或改用 `ISMEAR = 0`、`SIGMA = 0.05` |
| `Routine TETIRR needs special values` | 四面体方法与非 Γ 中心/不规则网格不兼容 | 改用 Γ 中心均匀网格；Line-mode 路径计算不要用 `ISMEAR = -5` |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角/三方体系用了 MP 偶数网格等与晶格对称性不相容的网格 | 改用 Γ 中心网格 |
| 能带与 SCF 在 Γ 点的本征值不一致 | 两步之间 ENCUT/PREC/POTCAR/LMAXMIX/MAGMOM 不一致，或 CHGCAR 不对应 | 复制 SCF 的 INCAR，只改 ICHARG、KPOINTS、NBANDS、LORBIT；重新拷贝 CHGCAR |
| `WARNING: The number of bands has been changed from X to Y` | 为了并行效率，VASP 自动增加了 NBANDS | 对能带作图无害；需要固定 NBANDS 时（LOBSTER、光学）设为并行组数的整数倍 |
| 能带"形状奇怪"、简并错误、带隙位置不在预期的高对称点 | 高对称点坐标与晶胞约定不一致 | 用生成路径工具输出的原胞重新 SCF；或用 ASE `bandpath` / Latimer–Munro 路径 |
| 金属的费米面附近能带整体平移 | 用了 ICHARG=11 路径计算的 E-fermi | 使用 SCF/DOS 计算的 E-fermi，或 `EFERMI = <值>` 固定 |
| 导带看不到 / DOS 在高能处突然截断 | NBANDS 太少 | 增加 NBANDS |
| DOS 锯齿严重、峰很粗 | `NEDOS` 太小或 k 点太稀 | NEDOS 加到 2001–3001，加密网格，用 `ISMEAR = -5` |
| PDOS 之和小于 TDOS | 间隙区电子未投影到 PAW 球内 | 正常现象，定性分析即可 |
| 用 `vasp_gam` 跑 KPOINTS_OPT 报错 | Γ-only 版本不支持 KPOINTS_OPT | 改用 `vasp_std` |
| HSE/SCAN 计算中用 `ICHARG = 11` | 杂化、meta-GGA 不能固定电荷密度做非自洽 | 改用零权重 k 点或 KPOINTS_OPT，见 [第6章](06-hybrid-metagga.md) |
| 自旋极化 KPOINTS_OPT 的 PDOS 缺少自旋向下通道（6.6.0） | known issue 69 | 升级到 6.6.1 |

---

## 练习

1. 用两步法与 KPOINTS_OPT 分别计算 Si 能带，把两者叠在一张图上，最大偏差应 < 1 meV。
2. 从能带上读出 Si 的 CBM 位置（Γ–X 线上的相对位置）与间接带隙，并与 SCF 网格上 BANDGAP 报告的值比较，解释差异。
3. 把 DOS 的 `NEDOS` 分别设为 301 和 3001，`ISMEAR` 分别设为 0（SIGMA = 0.05）和 −5，对比 DOS 曲线。
4. 对第4章优化好的单层 MoS2，用 VASPKIT 302 或 ASE 生成 Γ–M–K–Γ 路径，计算 PBE 能带，判断是直接还是间接带隙。
5. 用 py4vasp 画 Si 的 fat band，说明价带顶和导带底分别以哪种轨道（s 或 p）为主。

---

## 参考资料

**官方文档**

- [VASP Wiki: KPOINTS](https://www.vasp.at/wiki/index.php/KPOINTS) — 规则网格、显式列表、Line-mode 格式
- [VASP Wiki: KPOINTS_OPT](https://www.vasp.at/wiki/index.php/KPOINTS_OPT)、[KSPACING_OPT](https://www.vasp.at/wiki/index.php/KSPACING_OPT)、[PROCAR_OPT](https://www.vasp.at/wiki/index.php/PROCAR_OPT)、[KPOINTS_OPT_NKBATCH](https://www.vasp.at/wiki/index.php/KPOINTS_OPT_NKBATCH)
- [VASP Wiki: ICHARG](https://www.vasp.at/wiki/index.php/ICHARG)、[LMAXMIX](https://www.vasp.at/wiki/index.php/LMAXMIX)
- [VASP Wiki: LORBIT](https://www.vasp.at/wiki/index.php/LORBIT)、[DOSCAR](https://www.vasp.at/wiki/index.php/DOSCAR)、[NEDOS](https://www.vasp.at/wiki/index.php/NEDOS)、[ISMEAR](https://www.vasp.at/wiki/index.php/ISMEAR)
- [VASP Wiki: BANDGAP](https://www.vasp.at/wiki/index.php/BANDGAP)、[EFERMI](https://www.vasp.at/wiki/index.php/EFERMI)
- [VASP Changelog](https://www.vasp.at/wiki/index.php/Changelog) 与 [Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — KPOINTS_OPT 相关修复（issue 45、69）
- [py4vasp 文档](https://vasp.at/py4vasp/latest/) — `band`、`dos`、`bandgap`
- [VASPKIT 教程（1.5.x）](https://vaspkit.com/tutorials.html) — 任务 303/302、211–214、111–115、911
- [sumo 文档](https://smtg-bham.github.io/sumo/)、[SeeK-path 文档](https://seekpath.readthedocs.io/)、[pymatgen.io.vasp](https://pymatgen.org/pymatgen.io.vasp.html)、[pyprocar 文档](https://pyprocar.readthedocs.io/)

**相关 GitHub 仓库（仅链接）**

- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) — VASP 官方后处理库，读取 vaspout.h5
- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) — sumo-kgen / sumo-bandplot / sumo-dosplot
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — `HighSymmKpath`、`Vasprun.get_band_structure`（自动识别 KPOINTS_OPT）、`BSPlotter`
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar) — fat band、自旋纹理、费米面
- [vaspkit/VASPKIT_NatureProtocols](https://github.com/vaspkit/VASPKIT_NatureProtocols) — VASPKIT 论文配套示例
- [zjwang11/IRVSP](https://github.com/zjwang11/IRVSP) — 能带不可约表示，用于对称性标注
- [cpfpengfei/VASP-Examples](https://github.com/cpfpengfei/VASP-Examples) — 初学者 Si 能带/DOS 笔记

**关键文献**

- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* 49, 299 (2010). DOI: [10.1016/j.commatsci.2010.05.010](https://doi.org/10.1016/j.commatsci.2010.05.010)
- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* 128, 140 (2017). DOI: [10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015)
- J. M. Munro, K. Latimer, M. K. Horton, S. Dwaraknath, K. A. Persson, *npj Comput. Mater.* 6, 112 (2020)（Latimer–Munro 路径）. DOI: [10.1038/s41524-020-00383-7](https://doi.org/10.1038/s41524-020-00383-7)
- P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* 49, 16223 (1994)（四面体方法）. DOI: [10.1103/PhysRevB.49.16223](https://doi.org/10.1103/PhysRevB.49.16223)
- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* 13, 5188 (1976). DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)
- V. Wang et al., *Comput. Phys. Commun.* 267, 108033 (2021)（VASPKIT）. DOI: [10.1016/j.cpc.2021.108033](https://doi.org/10.1016/j.cpc.2021.108033)
- A. M. Ganose, A. J. Jackson, D. O. Scanlon, *J. Open Source Softw.* 3, 717 (2018)（sumo）. DOI: [10.21105/joss.00717](https://doi.org/10.21105/joss.00717)
- U. Herath et al., *Comput. Phys. Commun.* 251, 107080 (2020)（PyProcar）. DOI: [10.1016/j.cpc.2019.107080](https://doi.org/10.1016/j.cpc.2019.107080)

---

[← 上一章](04-relax-eos.md) | [返回目录](../README.md) | [下一章 →](06-hybrid-metagga.md)
