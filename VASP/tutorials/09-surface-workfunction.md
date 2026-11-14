# 第9章 表面模型、表面能、功函数与偶极校正（Surface Slabs, Surface Energy, Work Function & Dipole Correction with VASP）

本章讲解用 VASP 做表面计算的完整流程：用 ASE / pymatgen SlabGenerator / VASPKIT 切表面，选择层数、真空和 k 点，用 Selective dynamics 固定底层；计算表面能（对称 slab、单侧弛豫、Fiorentini–Methfessel 线性拟合、非化学计量表面）；用 LDIPOL/IDIPOL/DIPOL 做偶极校正，用 LVHAR 输出静电势并沿法向平面平均求功函数，以及 VASP 6.4.3+ 的 LVACPOTAV、WRT_POTENTIAL 与 py4vasp 新方法（VASP surface slab / surface energy / work function / dipole correction tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 用 ASE、pymatgen 或 VASPKIT 从体相结构生成任意 (hkl) 表面的 slab，判断对称/不对称、极性/非极性终止面；
- 对层数、真空、k 点做收敛测试，正确设置 Selective dynamics；
- 计算并正确报告表面能（J/m²），知道体相参考能量为何必须"完全一致"；
- 理解偶极校正的物理图像，会设置 `LDIPOL`、`IDIPOL`、`DIPOL`，处理偶极校正带来的 SCF 收敛问题；
- 用 LOCPOT 平面平均或 `LVACPOTAV` 求功函数，并能判断真空能级是否可信；
- 对表面重构、表面态和半导体表面钝化有基本认识。

| 计算目标 | 关键设置 | 本章小节 |
|---|---|---|
| 搭建 slab | 终止面、层数、真空、居中、固定底层 | §3–§5 |
| 表面能 | 对称 slab + 一致的体相参考；线性拟合 | §6 |
| 不对称 slab（单侧吸附、底层固定、极性面） | `LDIPOL = .TRUE.`、`IDIPOL = 3`、`DIPOL` | §7 |
| 功函数 | `LVHAR` + 平面平均；≥ 6.4.3 可用 `LVACPOTAV` | §8 |
| 重构、表面态、钝化 | 合适的超胞、投影能带、赝氢 | §9 |

> [!NOTE]
> 表面计算是吸附与催化（[第10章](10-adsorption-catalysis.md)）、过渡态（[第11章](11-neb-dimer.md)）的基础。本章使用的晶格常数必须来自你自己用**相同 ENCUT、POTCAR、泛函**优化得到的体相结构（[第4章](04-relax-eos.md)）。

---

## 2. 原理速览

### 2.1 slab 模型

平面波程序要求三维周期性，因此表面用"slab + 真空"的超胞表示：沿表面法向是若干原子层，层与其周期镜像之间隔着真空。一个 slab 有上下两个表面，模型的三个关键参数是：

- **表面超胞**：如 p(1×1)、p(2×2)、(√3×√3)R30°，决定吸附覆盖度和可描述的重构；
- **层数**：要足够厚，使中间层接近体相、上下表面互不影响；
- **真空厚度**：要足够大，使镜像之间的相互作用（波函数重叠、静电）可以忽略。

### 2.2 表面能

对上下两个表面完全相同（对称）的化学计量 slab：

$$
\gamma=\frac{E_{\rm slab}(N)-N\,\varepsilon_{\rm bulk}}{2A}
$$

N 是 slab 原子数（或化学式单元数），$\varepsilon_{\rm bulk}$ 是每原子（每化学式）体相能量，A 是表面超胞面积，分母的 2 对应两个表面。底部固定、只弛豫顶部的 slab，两个表面不同，需要拆分：

$$
\gamma_{\rm unrel}=\frac{E_{\rm slab}^{\rm unrel}-N\varepsilon_{\rm bulk}}{2A},\qquad
\gamma_{\rm rel}=\frac{E_{\rm slab}^{\rm rel}-N\varepsilon_{\rm bulk}}{A}-\gamma_{\rm unrel}
$$

非化学计量（或极性）slab 的表面能依赖化学势，例如 TiO₂ slab：

$$
\gamma(\mu_{\rm O})=\frac{E_{\rm slab}-N_{\rm Ti}E_{\rm TiO_2}^{\rm bulk}-(N_{\rm O}-2N_{\rm Ti})\,\mu_{\rm O}}{2A}
$$

其中 $\mu_{\rm O}$ 的范围由 $\tfrac12E_{\rm O_2}+\tfrac12\Delta H_f({\rm TiO_2})\le\mu_{\rm O}\le\tfrac12E_{\rm O_2}$ 限定（贫氧到富氧）。单位换算：1 eV/Å² = 16.02 J/m²。

### 2.3 功函数

功函数是把电子从费米能级移到"表面附近、但足够远"的真空所需的最小能量：

$$
\Phi=V_{\rm vac}-E_F
$$

$V_{\rm vac}$ 是真空中平坦区域的静电势（以 eV 计，即 $-e\phi$ 已乘上），$E_F$ 取自 OUTCAR。功函数是**表面**性质：同一晶体的不同晶面、不同终止、有无吸附，功函数都不同。

### 2.4 偶极校正

若 slab 两侧不对称（例如只在一侧吸附），它沿法向带有净偶极矩 $p_z$（每个表面超胞）。在三维周期边界下，这些偶极层周期排列，会在真空中产生一个人为的均匀电场，使真空区的电势出现斜率，既没有确定的真空能级，又会扭曲受力。偶极校正在真空中加一个与之抵消的锯齿形势（偶极层），结果是：两侧真空各自平坦，中间出现一个台阶，台阶高度就是两侧功函数之差：

$$
\Delta\Phi=\frac{e\,p_z}{\varepsilon_0A}\approx 180.95\times\frac{p_z\,[e\,\text{Å}]}{A\,[\text{Å}^2]}\ \text{eV}
$$

（Neugebauer & Scheffler, PRB 46, 16067 (1992)；Bengtsson, PRB 59, 12301 (1999)。）

---

## 3. 切表面：建模工具

通用的建模知识（数据库下载、超胞、表面、缺陷、异质结）见 [Common 第3章](../../Common/tutorials/03-structure-modeling.md) 与 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。这里给出与 VASP 表面计算直接相关的写法。

### 3.1 ASE

ASE 为常见晶面提供了现成函数（`fcc100/fcc110/fcc111/bcc100/bcc110/hcp0001` 等），任意 (hkl) 用 `surface()`：

```python
import numpy as np
from ase.build import bulk, surface, fcc111
from ase.constraints import FixAtoms
from ase.io import write

# 1) ready-made builder for fcc(111): p(2x2), 4 layers, 10 A vacuum on EACH side
slab = fcc111("Pt", size=(2, 2, 4), a=3.97, vacuum=10.0)

# 2) general (hkl) from any bulk; 'layers' counts repeat units along the normal
cu = bulk("Cu", "fcc", a=3.63, cubic=True)
s211 = surface(cu, (2, 1, 1), layers=6, vacuum=10.0)

# fix the bottom two layers -> written as Selective dynamics F F F
z = slab.positions[:, 2]
layers = np.unique(np.round(z, 2))
slab.set_constraint(FixAtoms(mask=z < layers[1] + 0.1))
write("POSCAR", slab, format="vasp", direct=True, sort=True)
```

> [!NOTE]
> `fcc111(..., vacuum=10.0)` 是在 slab **上下各加** 10 Å（总真空 20 Å），并且自动把 slab 居中。ASE 的 `fcc111` 默认给出 60° 夹角的六角表面胞，这对 VASP 没有问题（k 网格用 Γ 中心即可）。

### 3.2 pymatgen SlabGenerator

pymatgen 能自动找出某 (hkl) 的所有不同终止面，并判断 slab 是否对称、是否极性：

```python
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar

# conventional fcc Pt; better: bulk = Structure.from_file("CONTCAR") of YOUR relaxed bulk
bulk = Structure(Lattice.cubic(3.97), ["Pt"] * 4,
                 [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

gen = SlabGenerator(bulk, miller_index=(1, 1, 1),   # Miller index w.r.t. the INPUT cell
                    min_slab_size=4,                 # in units of (111) planes ...
                    min_vacuum_size=7,               # ... also in planes: 7 x 2.29 A ~ 16 A
                    in_unit_planes=True,
                    center_slab=True)
slabs = gen.get_slabs()                              # one Slab per distinct termination
print("terminations:", len(slabs))
slab = slabs[0].get_orthogonal_c_slab()             # make c perpendicular to the surface
slab.make_supercell([2, 2, 1])                       # p(2x2) surface cell
print(slab.composition, "| symmetric:", slab.is_symmetric(), "| polar:", slab.is_polar())
print("a, b, c = %.3f %.3f %.3f A, area = %.2f A^2" % (*slab.lattice.abc, slab.surface_area))

# fix the bottom half of the slab
z = slab.cart_coords[:, 2]
layers = np.unique(np.round(z, 2))
print("layers at z =", layers)
sd = [[False] * 3 if zi < layers[len(layers) // 2] - 0.01 else [True] * 3 for zi in z]
Poscar(slab, selective_dynamics=sd).write_file("POSCAR_Pt111_pmg")
```

运行输出（pymatgen 2026.10.24 实测）：

```text
terminations: 1
Pt16 | symmetric: True | polar: False
a, b, c = 5.614 5.614 25.213 A, area = 27.30 A^2
layers at z = [ 9.17 11.46 13.75 16.04]
```

要点：

- `in_unit_planes=True` 时，`min_slab_size` 和 `min_vacuum_size` **都**以"晶面层数"为单位；设为 `False` 时两者都以 Å 为单位。
- Miller 指数是相对于**输入晶胞**的。输入原胞和惯用胞时，同一个 (111) 指的不是同一族晶面——建议输入惯用胞。
- 多元化合物通常有多个终止面（`len(slabs) > 1`），要逐个检查化学计量、对称性和极性。`get_slabs(symmetrize=True)` 可尝试生成上下对称的 slab。
- 自动化的层数/真空收敛测试可参考 [SMTG-Bham/surfaxe](https://github.com/SMTG-Bham/surfaxe)（基于 pymatgen）。

### 3.3 VASPKIT

VASPKIT 1.5.x 的相关功能（编号以 [vaspkit.com](https://vaspkit.com/) 文档为准，不同版本可能变动）：

| 功能号 | 作用 |
|---|---|
| 803 | 按 Miller 指数切表面（Build Surface by Specified Miller Indices） |
| 801 | 在指定方向加真空（Build Vacuum Slab in Specified Direction） |
| 804 | 由两个 slab 建异质结 |
| 402 | 按层固定原子（生成 Selective dynamics 的 POSCAR_FIX） |
| 426 / 427 | LOCPOT 的平面平均 / 宏观平均势（功函数、能带对齐） |

```bash
vaspkit -task 803      # interactive prompts follow; or run 'vaspkit' and choose from the menu
```

### 3.4 终止面、对称性与极性

- **对称 slab**（上下表面等价）：没有净偶极，可以直接用 §2.2 的公式除以 2A；做表面能首选。
- **不对称 slab**：底部固定、单侧吸附、或两侧终止不同。需要偶极校正（§7），表面能需拆分公式。
- **极性表面**：按 Tasker 分类（J. Phys. C 12, 4977 (1979)），由交替带电原子面组成、每个重复单元有净偶极的是 Type III 极性面（如岩盐 (111)、纤锌矿 (0001)），理想截断的总能会随厚度发散，必须通过重构、空位、吸附或转移半层原子来补偿偶极。pymatgen 的 `is_polar()` 可做初步判断。

---

## 4. 层数、真空与 k 点

| 参数 | 起点建议 | 收敛判据 |
|---|---|---|
| 层数 | 金属 4–6 层；表面能 5–11 层；氧化物按"重复单元"计 | γ 变化 < 0.01–0.02 J/m²；功函数变化 < 0.02–0.05 eV；吸附能变化 < 0.02–0.05 eV |
| 真空 | ≥ 15 Å（功函数 15–20 Å） | 总能 < 1 meV/atom；真空区电势平台足够宽 |
| 面内 k 点 | 与体相相同的 k 点密度（倒空间长度 × 网格数相近） | 能量/功函数收敛 |
| 法向 k 点 | 1 | — |
| k 网格类型 | Γ 中心（六角/斜方晶胞必须） | — |
| smearing | 金属 `ISMEAR = 1/2`，`SIGMA ≈ 0.1–0.2`；半导体 `ISMEAR = 0`，`SIGMA ≈ 0.05` | T*S < 1 meV/atom；比较能量用 `energy(sigma->0)` |
| `LREAL` | 原子多时可用 `Auto`；但与体相参考保持一致 | — |

> [!WARNING]
> **不要对 slab 用 `ISIF = 3`**（真空区没有应力支撑，晶胞会塌缩）；slab 一律 `ISIF = 2`，面内晶格常数来自体相优化。也不要用 `ISMEAR = -5` 做 slab 弛豫：它给出的受力不准，且 k 点少时直接报 `Tetrahedron method fails`。

另外两条容易被忽略的规则：

1. **原子不要跨越 z 方向的周期边界**：所有原子的 z 分数坐标应在 (0, 1) 内且远离 0 和 1（slab 居中），否则平面平均、偶极校正、可视化都会出问题。
2. **slab 中心附近要像体相**：检查中心层的层间距、电荷、磁矩是否接近体相值，这是层数是否足够的直接证据。

---

## 5. 实战一：Pt(111) slab 的弛豫（固定底部两层）

完整输入见 [../examples/09-pt111-slab/README.md](../examples/09-pt111-slab/README.md)。

### 5.1 输入文件

POSCAR（ASE `fcc111` 生成，p(2×2)，4 层，a = 3.97 Å，上下各约 10 Å 真空；底部两层 `F F F`）：

```text
Pt(111) p(2x2), 4 layers, a = 3.97 A (PBE-like), 20 A vacuum, bottom 2 layers fixed
1.0
    5.6144278426   0.0000000000   0.0000000000
    2.8072139213   4.8622371394   0.0000000000
    0.0000000000   0.0000000000  26.8762417060
  Pt
  16
Selective dynamics
Direct
    0.0000000000   0.0000000000   0.3720758322   F F F
    0.5000000000   0.0000000000   0.3720758322   F F F
    0.0000000000   0.5000000000   0.3720758322   F F F
    0.5000000000   0.5000000000   0.3720758322   F F F
    0.1666666667   0.1666666667   0.4573586107   F F F
    0.6666666667   0.1666666667   0.4573586107   F F F
    0.1666666667   0.6666666667   0.4573586107   F F F
    0.6666666667   0.6666666667   0.4573586107   F F F
    0.8333333333   0.3333333333   0.5426413893   T T T
    0.3333333333   0.3333333333   0.5426413893   T T T
    0.8333333333   0.8333333333   0.5426413893   T T T
    0.3333333333   0.8333333333   0.5426413893   T T T
    0.0000000000   0.0000000000   0.6279241678   T T T
    0.5000000000   0.0000000000   0.6279241678   T T T
    0.0000000000   0.5000000000   0.6279241678   T T T
    0.5000000000   0.5000000000   0.6279241678   T T T
```

`Selective dynamics` 行必须在 `Direct` 之前；每个原子后面三个 `T/F` 分别对应 x、y、z 是否允许移动。

KPOINTS：

```text
Gamma-centred 8x8x1 for the p(2x2) Pt(111) slab (1 k-point along the surface normal)
0
Gamma
8 8 1
0 0 0
```

INCAR：

```text
SYSTEM  = Pt(111) p(2x2) 4-layer slab, relaxation of the top 2 layers

# ---- electronic ----
PREC    = Accurate
ENCUT   = 400          # use the SAME value for the bulk reference (surface energy)
EDIFF   = 1E-6
NELM    = 200
LREAL   = .FALSE.      # 16 atoms: reciprocal-space projection is affordable
ISMEAR  = 1            # metal: Methfessel-Paxton
SIGMA   = 0.10
ISPIN   = 1            # Pt is non-magnetic

# ---- ionic: only the 'T T T' atoms of POSCAR move ----
IBRION  = 2
ISIF    = 2            # never relax the cell of a slab
NSW     = 100
EDIFFG  = -0.02        # eV/A on the free atoms

# ---- dipole correction along z (bottom fixed / top relaxed -> slightly asymmetric) ----
LDIPOL  = .TRUE.
IDIPOL  = 3
DIPOL   = 0.5 0.5 0.5  # centre of mass of the slab in direct coordinates (z matters)

LWAVE   = .TRUE.       # reused by the work-function run
LCHARG  = .TRUE.
```

| 关键词 | 值 | 为什么 |
|---|---|---|
| `ISIF` | 2 | 只动原子；晶胞（面内晶格 + 真空）固定 |
| `EDIFFG` | −0.02 | 表面弛豫常用 0.01–0.03 eV/Å；固定原子的受力不参与判据 |
| `LDIPOL/IDIPOL` | `.TRUE.` / 3 | 底层固定、顶层弛豫，两侧略有不同；加上偶极校正更稳妥 |
| `DIPOL` | 0.5 0.5 0.5 | slab 已居中；只有 z 分量有意义，必须位于 slab 内部 |
| `LREAL` | `.FALSE.` | 小体系用倒空间投影更准；大 slab 可用 `Auto`，但体相参考要一致 |

```bash
cd ../examples/09-pt111-slab
cat $VASP_PP_PATH/potpaw_PBE/Pt/POTCAR > POTCAR
sbatch job.slurm
```

### 5.2 结果检查

```bash
grep "reached required accuracy" OUTCAR
grep "F=" OSZICAR | tail -3
grep -i "dipol" OUTCAR | tail -3
```

弛豫后的层间距（相对体相 (111) 面间距 $a/\sqrt3$）：

```python
import numpy as np
from ase.io import read

a = read("CONTCAR")
z = np.sort(a.positions[:, 2])
groups = np.split(z, np.where(np.diff(z) > 0.5)[0] + 1)    # a gap > 0.5 A starts a new layer
zl = np.array([g.mean() for g in groups])
d0 = 3.97 / np.sqrt(3)                                     # bulk (111) spacing
for i, di in enumerate(np.diff(zl)):
    print(f"d{i+1}{i+2} = {di:.3f} A  ({100 * (di - d0) / d0:+.2f} %)")
```

对未弛豫的 POSCAR 运行，所有层间距都等于 2.292 Å（+0.00%）；弛豫后只有顶部两层之间（d34）以及 d23 会变化，固定的 d12 保持不变。金属表面的层间弛豫通常只有百分之一量级。

---

## 6. 表面能

### 6.1 体相参考：误差的主要来源

表面能是两个大数之差除以面积，任何 slab 与体相之间不一致的设置都会被放大：

- ENCUT、PREC、POTCAR、`LREAL`、`ADDGRID`、泛函、vdW、smearing 必须完全相同；
- **k 点**：slab 的面内 k 网格与体相的 k 点密度不同，会使 $E_{\rm slab}(N)-N\varepsilon_{\rm bulk}$ 中残留一个正比于 N 的误差，表现为 γ 随层数**线性漂移**、永不收敛（Boettger 的经典问题）。

两种解决办法：

1. **取向体相（oriented bulk）**：用与 slab 面内晶胞相同、法向为若干层的体相晶胞计算 $\varepsilon_{\rm bulk}$，面内 k 网格与 slab 完全相同，法向 k 点按密度选取。
2. **Fiorentini–Methfessel 线性拟合**（J. Phys.: Condens. Matter 8, 6525 (1996)）：对一系列层数做线性拟合 $E_{\rm slab}(N)=2A\gamma+N\varepsilon$，斜率给出与 slab 设置自洽的 $\varepsilon$，截距给出 γ。

### 6.2 实战二：Cu(100) 表面能与层数收敛

完整输入见 [../examples/09-cu100-surface-energy/README.md](../examples/09-cu100-surface-energy/README.md)。流程：

```bash
cd ../examples/09-cu100-surface-energy
cat $VASP_PP_PATH/potpaw_PBE/Cu_pv/POTCAR > POTCAR
bash run.sh bulk       # fcc Cu, ISIF = 3, ENCUT = 500
bash run.sh slabs      # obulk/ + slab_N/ (relaxed) + slab_N_unrel/ for N = 3,5,7,9,11
bash run.sh collect    # surface_energy.py
```

- `make_slabs.py` 读 `bulk/CONTCAR` 求晶格常数，用 ASE `fcc100` 生成 1×1 slab（总真空 15 Å），**固定中心层**，上下两侧对称弛豫，可以直接用 $\gamma=(E-N\varepsilon)/2A$；
- `obulk/` 是 1×1 面内、2 层的取向体相，面内 16×16、法向 12 个 k 点；slab 用 16×16×1；
- 所有 INCAR 的 ENCUT（500 eV）、`ISMEAR = 1`、`SIGMA = 0.10` 相同，能量取 `energy(sigma->0)`。

`surface_energy.py` 的输出格式（示意，数值仅供参考）：

```text
oriented bulk: 2 atoms, eps_bulk = -x.xxxxxx eV/atom

== relaxed slabs ==
   N    E_slab (eV)  gamma (eV/A^2)  gamma (J/m^2)
   3     -xx.xxxxxx         0.09xxx         1.4xxx
   ...
linear fit over N = 7, 9, 11: eps_fit = ... (bulk ref ..., diff +0.0x meV)
Fiorentini-Methfessel gamma = 0.09xxx eV/A^2 = 1.4xxx J/m^2
```

判断标准：

- 取向体相参考下 γ(N) 应在 N ≥ 7 后基本不变；
- 线性拟合的 $\varepsilon_{\rm fit}$ 与 obulk 的 $\varepsilon_{\rm bulk}$ 相差应 < 1 meV/atom；
- 弛豫 γ 略低于未弛豫 γ；
- 参考量级：PBE 的 Cu(100) 表面能约 1.4–1.5 J/m²（Tran et al., Sci. Data 3, 160080 (2016)，可在 Materials Project 的 Crystalium 数据库查询）；实验"平均表面能"约 1.8 J/m²——PBE 系统性低估表面能是已知现象。

### 6.3 单侧弛豫的 slab

吸附研究中常见"底部固定、顶部弛豫"的 slab。此时先算同一 slab 的未弛豫能量 $E^{\rm unrel}$，再用 §2.2 的拆分公式求弛豫表面的 γ。快速计算：

```python
EV_A2_TO_J_M2 = 16.021766
def gamma_one_side(e_rel, e_unrel, n, eps_bulk, area):
    """Surface energy (J/m^2) of the relaxed side of a slab whose other side is bulk-terminated."""
    g_unrel = (e_unrel - n * eps_bulk) / (2.0 * area)
    g_rel = (e_rel - n * eps_bulk) / area - g_unrel
    return g_unrel * EV_A2_TO_J_M2, g_rel * EV_A2_TO_J_M2
```

### 6.4 Wulff 形貌（简介）

得到多个晶面的 γ 后，可用 pymatgen 的 `WulffShape`（`pymatgen.analysis.wulff`）构建平衡晶体形貌并计算各晶面所占面积比例，常用于纳米颗粒模型和催化活性晶面分析。

---

## 7. 偶极校正：LDIPOL / IDIPOL / DIPOL

### 7.1 标签

| 标签 | 推荐值 | 说明 |
|---|---|---|
| `LDIPOL` | `.TRUE.` | 在局部势中加入线性（带电体系还有二次）修正，**同时修正势、受力和能量**；需要配合 `IDIPOL` |
| `IDIPOL` | 3（表面法向沿 a3） | 1–3：只沿第 1/2/3 个晶格矢量计算偶极；4：三个方向都算（孤立分子用） |
| `DIPOL` | slab 质心（分数坐标） | 计算偶极矩的参考中心；表面只用法向分量。不设时 VASP 按平面平均电荷密度最小处自动确定（对正交晶胞较可靠） |
| `EPSILON` | 1（默认） | 体系所在介质的介电常数 |
| `EFIELD` | eV/Å | 沿 `IDIPOL`（1–3）方向施加外电场，官方建议同时打开 `LDIPOL`；注意 VASP 的电场方向定义与通常约定相反（电子沿 EFIELD 方向移动） |
| `LMONO` | 不要设置 | 会抑制偶极修正的输出 |

> [!IMPORTANT]
> 只写 `IDIPOL` 而不写 `LDIPOL`，VASP 只对**总能**做偶极/四极修正，势和受力不变——这对功函数和受力都不够。表面计算请写 `LDIPOL = .TRUE.` + `IDIPOL = 3`。

### 7.2 DIPOL 放在哪里

DIPOL 定义了偶极矩积分的中心；偶极层（电势台阶）出现在与 DIPOL 相隔半个晶胞的位置。因此：

- DIPOL 的 z 分量应位于 slab 的质心附近，这样台阶落在真空正中，离两侧表面都远；
- slab 居中时取 `DIPOL = 0.5 0.5 0.5`；slab 若放在晶胞底部（z ≈ 0.05–0.35），DIPOL 的 z 分量要改成约 0.2，否则台阶会落进原子层里；
- 结构弛豫中质心略有移动无妨，但不要让它移出 slab。

计算质心的分数坐标：

```python
from ase.io import read
atoms = read("POSCAR")
print("DIPOL = %.4f %.4f %.4f" % tuple(atoms.get_center_of_mass(scaled=True)))
```

### 7.3 收敛问题与对策

偶极校正常使 SCF 变慢：电荷需要在 slab 两侧之间转移，迭代中来回振荡（charge sloshing）。VASP Wiki 的建议：

```text
AMIN   = 0.01     # damp long-wavelength charge sloshing
NELM   = 200      # allow more electronic steps
EDIFF  = 1E-6     # an unconverged dipole correction gives wrong forces
```

仍不收敛时：先不加 `LDIPOL` 把波函数收敛，再读入 WAVECAR（`ISTART = 1`）打开偶极校正；同时检查真空是否足够（真空太小时电荷会"穿过"真空）。

### 7.4 带电 slab 与外电场

- 对带电体系（设置 `NELECT`），VASP 的势修正目前**只对立方超胞实现**，非立方晶胞 + `LDIPOL` 会直接停止。带电 slab 可考虑 VASP ≥ 6.5.0 的库仑核截断（`KERNEL_TRUNCATION`，见 [第8章](08-vdw-2d.md) §8.5）或隐式溶剂/补偿电荷方案（[第10章](10-adsorption-catalysis.md)）。
- `EFIELD` 施加外电场时，VASP < 6.6.0 不会自动去掉被电场破坏的反演对称性（known issue 50），结果可能错误；旧版本请设置 `ISYM = 0`。

### 7.5 输出

```bash
grep -i "dipol" OUTCAR | tail -5      # dipole moment along the selected direction
```

由输出的偶极矩 $p_z$（e·Å）和表面面积 A，可以用 §2.4 的公式估算两侧功函数之差，与 §8 平面平均得到的台阶高度互相验证。

---

## 8. 功函数

### 8.1 标准做法：LVHAR + 平面平均

1. 静态计算（通常接在弛豫之后），不对称 slab 打开偶极校正；
2. `LVHAR = .TRUE.`：LOCPOT 只写离子势 + Hartree 势（不含交换关联势）。交换关联势在真空中衰减很慢，只用静电部分可以在较小的真空下得到平坦的真空能级。注意 `LVHAR = .TRUE.` 会把 `LVTOT` 设为 `.FALSE.`；
3. 沿法向对 LOCPOT 做平面平均，取真空平坦区的平均值为 $V_{\rm vac}$；
4. $E_F$ 取 OUTCAR 中最后一个 `E-fermi`；$\Phi=V_{\rm vac}-E_F$。

INCAR（[../examples/09-pt111-workfunction-dipole/README.md](../examples/09-pt111-workfunction-dipole/README.md)）：

```text
SYSTEM    = Pt(111) slab: work function with dipole correction (static)

# ---- electronic ----
PREC      = Accurate
ENCUT     = 400
EDIFF     = 1E-6          # an unconverged dipole correction gives wrong potentials
NELM      = 200
LREAL     = .FALSE.
ISMEAR    = 1
SIGMA     = 0.10
ISPIN     = 1
AMIN      = 0.01          # damps charge sloshing between the two slab sides

# ---- static ----
IBRION    = -1
NSW       = 0

# ---- dipole correction along a3 ----
LDIPOL    = .TRUE.        # correct potential AND forces (IDIPOL alone: energy only)
IDIPOL    = 3
DIPOL     = 0.5 0.5 0.5   # centre of mass of the slab (direct coordinates)

# ---- potential output ----
LVHAR     = .TRUE.        # LOCPOT = ionic + Hartree (no xc): flattest vacuum level
LVACPOTAV = .TRUE.        # VASP >= 6.4.3: vacuum levels of both sides in OUTCAR
#WRT_POTENTIAL = hartree ionic   # VASP >= 6.4.3: also store potentials in vaspout.h5

LWAVE     = .FALSE.
LCHARG    = .FALSE.
```

平面平均的核心只有几行（LOCPOT 中 x 下标变化最快，因此用 Fortran 顺序重排）：

```python
import numpy as np
# pot: flat list of NGX*NGY*NGZ values read after the grid line of LOCPOT
pot = np.array(values).reshape((ngx, ngy, ngz), order="F")
v_z = pot.mean(axis=(0, 1))           # planar average along a3
z = np.arange(ngz) * c_length / ngz
```

算例中的 [plot_locpot_planar.py](../examples/09-pt111-workfunction-dipole/plot_locpot_planar.py) 实现了完整流程：读 LOCPOT 与 E-fermi、自动识别最大的无原子区间和偶极台阶、分别平均两侧真空、报告真空区残余斜率并作图：

```bash
python plot_locpot_planar.py LOCPOT --outcar OUTCAR            # --axis 3 --margin 3.0 by default
```

输出格式（示意，数值仅供参考）：

```text
dipole step     : yes, -0.xxx eV
E_F             : -x.xxxx eV
V_vac upper side: x.xxxx eV (std 0.00xx, slope +0.000x eV/A)
V_vac lower side: x.xxxx eV (std 0.00xx, slope -0.000x eV/A)
Phi upper       : x.xxxx eV
Phi lower       : x.xxxx eV
Delta Phi       : +0.xxxx eV (upper - lower)
```

仓库的通用脚本 [../scripts/locpot_workfunction.py](../scripts/locpot_workfunction.py) 提供平面平均、宏观平均和功函数计算，用 `python ../scripts/locpot_workfunction.py -h` 查看用法。

> [!WARNING]
> 用 ASE 的 `VaspChargeDensity` 读 LOCPOT 要小心：它按 CHGCAR 的约定把数据除以晶胞体积，读势文件时需要再乘回体积，否则数值完全不对。直接按文件格式读取（如上面的脚本）最稳妥。

### 8.2 VASP ≥ 6.4.3：LVACPOTAV

`LVACPOTAV = .TRUE.` 让 VASP 在沿 `IDIPOL` 方向的无场区（平面平均电场小于 `VACPOTFLAT`，默认 0.1 eV/Å，且平均电荷密度近乎为零）中求平均，直接在 OUTCAR 给出两侧真空能级：

```bash
grep upper OUTCAR
grep "E-fermi" OUTCAR | tail -1
```

输出格式（数值仅为示意）：

```text
 vacuum level on the upper side and lower side of the slab         x.xxx         x.xxx
```

功函数 = 对应一侧的真空能级 − E_F。它是后处理步骤，也可以对已收敛的计算用 `ALGO = None` 续算得到。官方建议与 `LVHAR` 一起使用。注意事项：

- 只支持 `IDIPOL = 1–3`；
- 偶极体系必须同时打开偶极校正，否则不存在无场区；
- 若真空太小，OUTCAR 会出现警告：

```text
 |     Did not find any points to average over, which means that no vacuum     |
 |     field-free region was found. Please increase the size of  your cell     |
```

或"最小电荷密度乘以体积大于 1E-1"的警告——解决办法是增大真空、收紧 `EDIFF`。增大 `VACPOTFLAT` 可以快速粗估，减小它可以得到更严格的值。

### 8.3 WRT_POTENTIAL 与 py4vasp

VASP ≥ 6.4.3 还可以把势写进 `vaspout.h5`（需要 HDF5 编译）：

```text
WRT_POTENTIAL = hartree ionic     # or: total / xc ; several can be listed
```

用 py4vasp 读取（0.11.x）：

```python
import py4vasp
calc = py4vasp.Calculation.from_path(".")
pot = calc.potential.read()          # dict with 'hartree', 'ionic', ... on the FFT grid
print(calc.workfunction)             # needs IDIPOL + LVACPOTAV: vacuum levels, E_F, VBM/CBM
calc.workfunction.plot()             # averaged potential along the IDIPOL direction
```

`WRT_POTENTIAL` 写出的 Hartree 势与离子势各自的 G = 0 分量被设为零，单独看任一项没有意义，应看两者之和。

### 8.4 其他工具

- **VASPKIT 426 / 427**：LOCPOT 平面平均与宏观平均（宏观平均在界面/异质结能带对齐中常用）；
- **[MacroDensity](https://github.com/WMD-group/MacroDensity)**：平面平均、宏观平均、真空能级与能带对齐的 Python 工具；
- **py4vasp**：见 §8.3。

### 8.5 实战三：干净与 O 覆盖 Pt(111) 的功函数

```bash
cd ../examples/09-pt111-workfunction-dipole
cat $VASP_PP_PATH/potpaw_PBE/Pt/POTCAR > POTCAR_Pt
cat $VASP_PP_PATH/potpaw_PBE/{Pt,O}/POTCAR > POTCAR_PtO
bash run.sh prepare        # clean/ (relaxed slab), O_nodip/, O_dip/
bash run.sh collect
```

预期现象（定性）：

| 计算 | 真空区电势 | 结论 |
|---|---|---|
| `clean/` | 两侧平坦，ΔΦ ≈ 0 | 对称良好；PBE 的 Pt(111) 功函数文献值多在 5.6–5.9 eV（近似；实验约 5.9–6.1 eV；系统性数据见 Tran et al., Surf. Sci. 687, 48 (2019)） |
| `O_nodip/` | 有斜率，脚本报 `WARNING: potential is not flat` | 无偶极校正，功函数不确定 |
| `O_dip/` | 两侧平坦，中间有台阶 | O 一侧功函数更高，ΔΦ 与 OUTCAR 偶极矩按 §2.4 估算一致 |

### 8.6 功函数的收敛与常见误区

- 对层数、真空、面内 k 点分别收敛；对金属，E_F 对 k 点和 smearing 比较敏感；
- 半导体的 E_F 落在带隙中某处，数值随 smearing 和缺陷而变，"功函数"意义不大——改报电离能（V_vac − VBM）和电子亲和能（V_vac − CBM）；
- 比较不同吸附覆盖度的功函数变化 ΔΦ 时，各计算的 slab 厚度、真空和 k 点要一致；
- 不要用含交换关联势的 `LVTOT` LOCPOT 直接取"真空值"，除非真空非常大、已确认平坦。

---

## 9. 表面重构与表面态（简介）

### 9.1 重构

许多表面在理想截断后会重构以降低能量，例如 Si(001)-(2×1) 二聚体、Si(111)-(7×7)、Au(111) 的 22×√3 "人字形"重构、Pt(100)/Au(100) 的 hex 重构、GaAs(001) 的多种重构。计算要点：

- 超胞必须与重构周期相容（如 Si(001) 至少 2×1，研究二聚体倾斜需 p(2×2) 或 c(4×2)）；
- 比较不同重构的能量要换算到相同面积（如每 1×1 单元），化学计量不同时还要引入化学势（§2.2）；
- 初始结构要"打破对称"（如二聚体略微倾斜），否则优化可能停在高对称的鞍点。

### 9.2 半导体表面的钝化

只想研究一侧表面时，可以用赝氢（分数价电子的 H，如 `H.75`、`H1.25`、`H.5`、`H1.5`，VASP 赝势库提供 0.25–1.75 价的一系列）饱和底部的悬挂键，消除底面的表面态。例如对 III–V 半导体，阳离子悬挂键用 1.25 价的 H，阴离子悬挂键用 0.75 价的 H（对应 3/4 与 5/4 个电子的悬挂键）。注意：POTCAR 中价电子数的位数有限（例如三个 `H.33` 加起来是 0.99 而不是 1），需要时用 `NELECT` 修正总电子数。

### 9.3 表面态

- 用较厚的 slab 计算能带（[第5章](05-band-dos.md)，二维 k 路径如 Γ̄–M̄–K̄–Γ̄），把能带投影到表层原子（`LORBIT = 11`，PROCAR）与中心层原子，表层权重高且落在体相投影能带"空隙"中的就是表面态；
- 体相投影能带：对体相沿法向 k_z 扫描，把所有 k_z 的本征值叠加到表面布里渊区上；
- slab 太薄时上下两个表面的表面态会相互耦合而劈裂，增加层数应使其简并；
- 经典例子：Cu(111)、Au(111) 在 Γ̄ 附近的 Shockley 表面态；加入 SOC 后 Au(111) 表面态出现 Rashba 劈裂（[第7章](07-magnetism-u-soc.md) §8）。投影能带可用 [pyprocar](https://github.com/romerogroup/pyprocar) 或 [sumo](https://github.com/SMTG-Bham/sumo) 绘制。

---

## 参数选择与收敛建议

| 项目 | 建议 | 备注 |
|---|---|---|
| 体相晶格常数 | 用与 slab 相同的 ENCUT/POTCAR/泛函/vdW 优化 | 面内晶格从这里来 |
| 层数 | 收敛 γ、Φ 或吸附能；金属 ≥ 4 层起步 | 中心层应接近体相 |
| 固定层 | 底部 1/2 左右（如 4 层固定 2 层）；对称 slab 固定中心层 | 固定层越多，越需要更厚的 slab |
| 真空 | ≥ 15 Å；功函数 15–20 Å | `LVHAR` 可减小所需真空 |
| k 点 | Γ 中心，面内与体相等密度，法向 1 | 表面能用取向体相参考 |
| smearing | 金属 MP（`ISMEAR = 1`、`SIGMA` 0.1–0.2）；半导体 `ISMEAR = 0` | 能量用 `energy(sigma->0)` |
| 表面能参考 | 取向体相或线性拟合 | γ 随 N 线性漂移说明参考不一致 |
| 偶极校正 | 不对称 slab 一律 `LDIPOL = .TRUE.`、`IDIPOL = 3`、`DIPOL` = 质心 | 配合 `AMIN = 0.01`、`EDIFF = 1E-6` |
| 功函数 | `LVHAR = .TRUE.`；≥ 6.4.3 加 `LVACPOTAV` | 检查真空区斜率 < 0.01 eV/Å |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 功函数计算中真空区电势不平坦（有斜率） | 不对称 slab 没有偶极校正；偶极校正未收敛；真空太小 | `LDIPOL = .TRUE.`、`IDIPOL = 3`、`DIPOL` = 质心；收紧 `EDIFF`；增大真空 |
| 打开 `LDIPOL` 后 SCF 振荡、迟迟不收敛 | 电荷在两侧之间来回振荡 | `AMIN = 0.01`、增大 `NELM`；先不加 LDIPOL 收敛再续算 |
| 电势台阶出现在原子层中间，功函数离谱 | `DIPOL` 不在 slab 质心（如 slab 在晶胞底部却用 0.5） | 用质心的分数坐标设置 DIPOL；或把 slab 居中 |
| `Did not find any points to average over, which means that no vacuum field-free region was found ...` | `LVACPOTAV` 找不到无场真空区 | 增大真空；确认偶极校正打开；收紧 EDIFF |
| `The minimum charge density times volume of the cell along the axis of the dipole correction is larger 1E-1 ...` | 真空中仍有明显电荷密度 | 增大真空、提高精度 |
| 带电 slab 打开 `LDIPOL` 后 VASP 停止 | 带电体系的势修正只支持立方超胞 | 用立方超胞；或改用库仑核截断等方案 |
| 表面能随层数线性增大/减小、不收敛 | 体相参考与 slab 设置不一致（k 点、smearing、ENCUT、LREAL） | 用取向体相参考或 Fiorentini–Methfessel 线性拟合 |
| slab 用 `ISIF = 3` 后真空塌缩 | 真空方向没有应力支撑 | slab 一律 `ISIF = 2` |
| `Tetrahedron method fails (number of k-points < 4)` | slab 法向只有 1 个 k 点，`ISMEAR = -5` 不适用 | 改用 `ISMEAR = 0/1` |
| `VERY BAD NEWS! internal error in subroutine IBZKPT ...` | 六角表面胞用了 Monkhorst–Pack 偶数网格 | 改为 Gamma 中心网格 |
| 原子"穿过"z 边界，可视化中 slab 被劈成两半 | slab 未居中、坐标接近 0 或 1 | 居中 slab（`atoms.center(axis=2)`）后重新计算 |
| 极性表面能量随厚度发散 | Tasker III 型极性面的宏观偶极 | 重构、空位或吸附补偿偶极；检查 pymatgen `is_polar()` |
| `EFIELD` 计算结果对称得不合理（< 6.6.0） | 外电场方向的反演对称没有被去掉（known issue 50） | `ISYM = 0` 或升级 |

---

## 练习

1. 用 ASE 或 pymatgen 生成 Pt(111) 的 3、4、5、6、7 层 slab，计算功函数随层数的收敛曲线。
2. 在 Cu(100) 算例中，用"原胞体相能量（16×16×16 k 点）"代替取向体相作为参考，观察 γ(N) 的线性漂移，再用线性拟合修正。
3. 对 `O_dip/` 的 OUTCAR 读出偶极矩 $p_z$，用 $\Delta\Phi \approx 180.95\,p_z/A$ 估算功函数差，与平面平均得到的台阶比较。
4. 把 Pt(111) slab 平移到晶胞底部（z ≈ 0.05–0.40），分别用 `DIPOL = 0.5 0.5 0.5` 和正确的质心计算，比较平面平均势的形状。
5. 用 pymatgen 生成 TiO₂ 金红石 (110) 的所有终止面，判断哪些是化学计量且对称的，写出 γ(μ_O) 的表达式。

---

## 参考资料

**VASP 官方文档**

- LDIPOL：<https://www.vasp.at/wiki/index.php/LDIPOL>；IDIPOL：<https://www.vasp.at/wiki/index.php/IDIPOL>；DIPOL：<https://www.vasp.at/wiki/index.php/DIPOL>
- 功函数计算指南：<https://www.vasp.at/wiki/index.php/Computing_the_work_function>
- LVHAR：<https://www.vasp.at/wiki/index.php/LVHAR>；LVACPOTAV：<https://www.vasp.at/wiki/index.php/LVACPOTAV>；VACPOTFLAT：<https://www.vasp.at/wiki/index.php/VACPOTFLAT>；WRT_POTENTIAL：<https://www.vasp.at/wiki/index.php/WRT_POTENTIAL>
- 静电修正：<https://www.vasp.at/wiki/index.php/Electrostatic_corrections>
- 赝势选择（含分数价赝氢）：<https://www.vasp.at/wiki/index.php/Choosing_pseudopotentials>
- py4vasp：<https://vasp.at/py4vasp/latest/>
- VASP 表面教程：<https://www.vasp.at/tutorials/latest/>

**相关 GitHub 仓库（仅链接）**

- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：`SlabGenerator`、`WulffShape`、表面分析
- [SMTG-Bham/surfaxe](https://github.com/SMTG-Bham/surfaxe)：自动生成 slab、层数/真空收敛、表面能与静电势分析
- [WMD-group/MacroDensity](https://github.com/WMD-group/MacroDensity)：LOCPOT 平面/宏观平均、真空能级与能带对齐
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：`calc.potential`、`calc.workfunction`
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book)：VASP + ASE 表面章节（表面能、功函数、偶极校正）
- [skelton-group/VASP-Workshop-Chester-2019](https://github.com/skelton-group/VASP-Workshop-Chester-2019)：NaCl(100) 表面与吸附案例（CC-BY）
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar)：按原子层投影的能带，识别表面态

**关键文献**

- J. Neugebauer, M. Scheffler, PRB 46, 16067 (1992)（偶极校正）. DOI: [10.1103/PhysRevB.46.16067](https://doi.org/10.1103/PhysRevB.46.16067)
- L. Bengtsson, PRB 59, 12301 (1999)（偶极校正）. DOI: [10.1103/PhysRevB.59.12301](https://doi.org/10.1103/PhysRevB.59.12301)
- V. Fiorentini, M. Methfessel, J. Phys.: Condens. Matter 8, 6525 (1996)（表面能线性拟合）. DOI: [10.1088/0953-8984/8/36/005](https://doi.org/10.1088/0953-8984/8/36/005)
- P. W. Tasker, J. Phys. C 12, 4977 (1979)（离子晶体表面稳定性/极性面分类）. DOI: [10.1088/0022-3719/12/22/036](https://doi.org/10.1088/0022-3719/12/22/036)
- W. Sun, G. Ceder, Surf. Sci. 617, 53 (2013)（slab 的高效生成与收敛）. DOI: [10.1016/j.susc.2013.05.016](https://doi.org/10.1016/j.susc.2013.05.016)
- R. Tran et al., Sci. Data 3, 160080 (2016)（元素晶体表面能数据库）. DOI: [10.1038/sdata.2016.80](https://doi.org/10.1038/sdata.2016.80)
- R. Tran et al., Surf. Sci. 687, 48 (2019)（元素晶体各向异性功函数）. DOI: [10.1016/j.susc.2019.05.002](https://doi.org/10.1016/j.susc.2019.05.002)
- S. P. Ong et al., Comput. Mater. Sci. 68, 314 (2013)（pymatgen）. DOI: [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)

[← 上一章](08-vdw-2d.md) | [返回目录](../README.md) | [下一章 →](10-adsorption-catalysis.md)
