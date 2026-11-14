# 第3章 结构准备：数据库下载、建模（超胞/表面/缺陷/掺杂/异质结）与格式转换（Structure Preparation: Databases, Supercells, Slabs, Defects, Doping, Heterostructures & Format Conversion）

任何 DFT 计算都从一个“对的结构”开始。本章讲清楚如何用 Materials Project 新版 API（mp-api / pymatgen `MPRester`）、COD、PubChem 等数据库获取结构，如何用 spglib 的 `symprec` 判断对称性、得到原胞与标准惯用胞（`get_primitive_standard_structure`），如何用 ASE 与 pymatgen 搭建超胞（`make_supercell`）、表面 slab（`ase.build.surface`、`SlabGenerator`，含真空层与固定底层）、吸附结构（`AdsorbateSiteFinder`）、空位/替位/间隙缺陷与掺杂、二维异质结，以及如何用 RDKit/Avogadro 建分子，最后用 `ase.io` 把 CIF 转成 VASP `POSCAR`、xyz、Gaussian `.gjf`、CP2K `&COORD` 和 OpenMX `Atoms.SpeciesAndCoordinates`，并用 VESTA 检查键长和对称性（CIF/POSCAR/xyz/gjf conversion, supercell, slab, adsorption, defects, heterostructure, molecules）。

> [!NOTE]
> 本章所有脚本都在 `Common/examples/03-*` 下，按主题分为 5 个目录（见 1.2 节）。脚本只负责“生成结构”，不调用任何 DFT 程序，在笔记本上就能运行。ASE、pymatgen 更系统的用法（计算器、输入集、输出解析）见 [第4章](04-ase-pymatgen.md)；作业提交见 [第2章](02-job-scheduler.md)。

---

## 1. 本章目标与适用场景

### 1.1 学完本章你应当能够

| 目标 | 具体内容 | 对应小节 |
|---|---|---|
| 从数据库拿结构 | MP 新 API key、按 mp-id/化学式检索、按 E_hull 选稳定相；COD/ICSD 的实验 CIF；PubChem 分子 | 3 |
| 判断对称性 | spglib `symprec` 扫描；原胞、标准惯用胞；Wyckoff 位置与不等价位点 | 4 |
| 格式转换 | CIF → POSCAR / xyz / extxyz / gjf / CP2K / OpenMX，知道每种格式丢掉了什么信息 | 5 |
| 超胞 | 对角扩胞、一般变换矩阵、近立方超胞、√3×√3 R30° | 6 |
| 表面 slab | Miller 指数、终止面、极性、层数与真空、固定底层（Selective dynamics） | 7 |
| 吸附 | top/bridge/hollow 位点、初始高度、吸附能所需的三个结构 | 8 |
| 缺陷与掺杂 | 空位、替位、间隙、双掺杂对称去重、随机掺杂、带电缺陷的注意事项 | 9 |
| 二维异质结 | 晶格匹配、应变与扭转角、层间距与真空 | 10 |
| 分子 | SMILES → 3D、电荷与多重度检查、分子放进周期盒子 | 11 |
| 检查 | VESTA/ASE 检查键长、原子重叠、真空、元素顺序 | 12 |

适用场景：刚开始一个新体系、要从文献或数据库复现一个结构、要搭表面/缺陷/界面模型，或者要把同一个结构交给不同软件计算。

### 1.2 配套示例一览

| 示例目录 | 内容 | 主要工具 |
|---|---|---|
| [03-fetch-convert-si](../examples/03-fetch-convert-si/README.md) | MP 下载 Si、对称性分析、一键转四款软件输入 | mp-api、pymatgen、spglib、ASE |
| [03-slab-adsorption](../examples/03-slab-adsorption/README.md) | Cu(111) slab、TiO2(110) 终止面筛选、吸附位点 | ASE `fcc111`、`SlabGenerator`、`AdsorbateSiteFinder` |
| [03-defect-doping-mgo](../examples/03-defect-doping-mgo/README.md) | MgO 空位、替位、间隙、双掺杂去重 | ASE、spglib、`StructureMatcher` |
| [03-heterostructure-2d](../examples/03-heterostructure-2d/README.md) | 石墨烯/MoS2 等六方二维异质结的晶格匹配 | ASE `graphene`、`mx2`、`make_supercell` |
| [03-molecules](../examples/03-molecules/README.md) | SMILES → 3D，分子 → gjf/CP2K/POSCAR/OpenMX | RDKit、ASE |

此外，Common 板块提供两个通用辅助脚本：[convert_any.py](../scripts/convert_any.py)（任意格式互转，可选原胞/标准胞、超胞、加真空、按层固定）和 [kpath_gen.py](../scripts/kpath_gen.py)（标准高对称路径），用 `python ../scripts/convert_any.py -h`、`python ../scripts/kpath_gen.py -h` 查看用法，参数说明见 [Common/scripts/README.md](../scripts/README.md)。本章的核心转换功能在 [cif_to_all.py](../examples/03-fetch-convert-si/cif_to_all.py) 和 [xyz_to_inputs.py](../examples/03-molecules/xyz_to_inputs.py) 中已经完整实现，不依赖这两个脚本。仓库作者自己的建模脚本（随机掺杂替换、球面碳片段）在 [Scripts/modeling](../../Scripts/README.md)，第 9、10 节会用到。

---

## 2. 原理速览

### 2.1 晶胞、分数坐标与周期性

晶体由三个晶格矢量 $\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3$（按行组成矩阵 $\mathbf{h}$）和晶胞内的原子位置描述。分数坐标（POSCAR 的 `Direct`、OpenMX 的 `FRAC`）与笛卡尔坐标的关系为

$$
\mathbf{r} = x\,\mathbf{a}_1 + y\,\mathbf{a}_2 + z\,\mathbf{a}_3,\qquad V=\left|\mathbf{a}_1\cdot(\mathbf{a}_2\times\mathbf{a}_3)\right| .
$$

周期性边界条件下，平面波程序（VASP、CP2K 的 GPW、OpenMX）**总是**周期的：分子、表面也必须放进一个晶胞，用真空隔开周期像。Gaussian 默认是孤立分子，只有写了平移矢量 `TV` 才做周期计算。

### 2.2 超胞变换矩阵

任意超胞都可以写成整数矩阵 $\mathbf{P}$ 作用在原晶胞上：

$$
\mathbf{h}_{\text{super}} = \mathbf{P}\,\mathbf{h},\qquad N_{\text{super}} = |\det \mathbf{P}|\;N_{\text{prim}} .
$$

对角矩阵 `diag(2,2,2)` 就是常见的 2×2×2 扩胞。六方面内的 $\sqrt3\times\sqrt3\,R30^\circ$ 超胞（$\det\mathbf{P}=3$）的矩阵**取决于面内晶格夹角的约定**：

- $\gamma=120^\circ$ 约定（ASE 的 `graphene()`、`mx2()`，多数二维材料数据库）：$\mathbf{P}=\begin{pmatrix}2&1&0\\-1&1&0\\0&0&1\end{pmatrix}$；
- $\gamma=60^\circ$ 约定（ASE 的 `fcc111(..., orthogonal=False)`、`hcp0001`）：$\mathbf{P}=\begin{pmatrix}1&1&0\\-1&2&0\\0&0&1\end{pmatrix}$。

用错约定得到的就不是 $\sqrt3$ 超胞（例如把第一个矩阵用在 60° 晶胞上，第一条边长变成 $\sqrt7\,a$）。动手前先用 `atoms.cell.angles()` 看一下 $\gamma$。

### 2.3 Miller 指数、层间距与 slab

(hkl) 面的面间距 $d_{hkl}=2\pi/|\mathbf{G}_{hkl}|$，立方晶系为 $d_{hkl}=a/\sqrt{h^2+k^2+l^2}$，例如 fcc(111) 的相邻原子层间距为 $a/\sqrt3$。slab 模型由“若干原子层 + 真空层”组成；表面能

$$
\gamma = \frac{E_{\text{slab}} - N\,E_{\text{bulk}}}{2A}
$$

要求上下两个表面相同（对称 slab）、化学计量与体相一致，否则 $\gamma$ 没有明确意义。

### 2.4 对称性容差 symprec

数值结构（尤其是弛豫后的 CONTCAR）里原子位置有 $10^{-4}$～$10^{-2}$ Å 量级的噪声。spglib 判断两个原子“对称等价”的依据是：对称操作后位置差小于 `symprec`（单位 Å）。`symprec` 太小，弛豫后结构会降到低对称甚至 P1；太大，会把真实的畸变“抹平”成高对称结构。

### 2.5 缺陷周期像与晶格失配

超胞中的点缺陷会与其周期像相互作用，误差随超胞边长 $L$ 衰减（带电缺陷的静电项约按 $1/L$ 衰减），所以缺陷超胞通常要求 $L\gtrsim 10$ Å 并做尺寸收敛。两种材料搭界面时，上层的应变

$$
\varepsilon = \frac{|\mathbf{v}_{\text{bottom}}|}{|\mathbf{v}_{\text{top}}|} - 1
$$

其中 $\mathbf{v}$ 是两层各自的公共超胞矢量（Zur–McGill 匹配思想），一般控制在 1～2% 以内。

### 2.6 各格式保存了什么

| 格式 | 晶胞 | 周期性 | 元素 | 固定信息 | 磁矩/电荷 | 典型用途 |
|---|---|---|---|---|---|---|
| CIF | 有 | 隐含 3D | 有 | 无 | 无（可有占位度） | 数据库、文献交换 |
| POSCAR/CONTCAR | 有 | 3D | 第 6 行（VASP 5+） | `Selective dynamics` | 无（在 INCAR） | VASP |
| xyz | **无** | 无 | 有 | 无 | 无 | 分子、可视化 |
| extxyz | 有（`Lattice=`） | `pbc=` | 有 | 可存 | 可存 | ASE、MLIP 数据集 |
| Gaussian .gjf | 仅 `TV` 行 | `TV` 行数 | 有 | 冻结标志 | 电荷/多重度 | Gaussian |
| CP2K `&CELL`/`&COORD` | 有 | `PERIODIC` | 有 | `&FIXED_ATOMS` | `&KIND` | CP2K |
| OpenMX `.dat` | `Atoms.UnitVectors` | 3D | 物种名 | `MD.Fixed.XYZ` | up/down 电子数 | OpenMX |

记住这张表就能避开最常见的转换错误：**xyz 不存晶胞**，CIF 不存固定信息，POSCAR 不存元素到 POTCAR 的对应（顺序要自己保证）。

---

## 3. 从数据库获取结构

### 3.1 常用数据库

| 数据库 | 内容 | 获取方式 | 备注 |
|---|---|---|---|
| [Materials Project](https://next-gen.materialsproject.org/) | 十余万种无机晶体的 DFT 弛豫结构、能量、能带 | 网页、mp-api、pymatgen | 需要免费 API key；结构是 PBE(+U)/r2SCAN 弛豫的 |
| [COD](https://www.crystallography.net/cod/) | 开放的实验晶体结构（含大量有机/金属有机） | 网页下载 CIF、REST | 免费、可再分发 |
| [ICSD](https://icsd.products.fiz-karlsruhe.de/) | 最全的无机实验结构 | 机构订阅 | 有版权，不要把 CIF 放进公开仓库 |
| [OQMD](https://oqmd.org/)、[AFLOW](https://aflow.org/)、[JARVIS](https://jarvis.nist.gov/) | 高通量 DFT 数据库 | 网页、API、OPTIMADE | 与 MP 互补，可交叉验证 |
| [C2DB](https://c2db.fysik.dtu.dk/) | 二维材料单层结构与性质 | 网页、ASE db | 二维/异质结建模首选 |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | 分子（SMILES、3D SDF） | 网页、PUG REST | Gaussian 分子起点 |
| [CCDC/CSD](https://www.ccdc.cam.ac.uk/) | 有机与金属有机晶体 | 订阅（部分免费查询） | 分子晶体、MOF |
| [OPTIMADE](https://www.optimade.org/) | 统一查询上述多个数据库的标准 API | REST、`optimade-python-tools` | 一次查询多个库 |

### 3.2 Materials Project 新版 API

> [!IMPORTANT]
> MP 旧版 REST API（`materialsproject.org/rest/v2`）已于 **2025 年 9 月下线**，旧的 16 位 key 失效。现在必须使用 dashboard 上显示的 **32 位新 key**，客户端用 `mp-api`，或新版 pymatgen 自带的轻量 `MPRester`（不依赖 mp-api）。参见 [MP API 入门文档](https://docs.materialsproject.org/downloading-data/using-the-api/getting-started)。

第一步，登录 <https://next-gen.materialsproject.org/dashboard> 复制 key，写进环境变量（`~/.bashrc`，见 [第1章](01-linux-hpc.md) 第 8.4 节）：

```bash
pip install -U mp-api pymatgen
export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    # used by mp-api
export PMG_MAPI_KEY=$MP_API_KEY                       # used by pymatgen.ext.matproj
unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT               # remove old endpoints pointing to rest/v2 (mp-api / pymatgen)
```

最小的 Python 用法如下（mp-api）：

```python
# fetch_min.py - minimal mp-api example (new API, 2025+)
from mp_api.client import MPRester

with MPRester() as mpr:                      # key read from $MP_API_KEY
    docs = mpr.materials.summary.search(
        formula="TiO2",
        fields=["material_id", "formula_pretty", "symmetry",
                "energy_above_hull", "band_gap"])
    for d in sorted(docs, key=lambda d: d.energy_above_hull):
        print(d.material_id, d.symmetry.symbol, round(d.energy_above_hull, 3))
    st = mpr.get_structure_by_material_id("mp-2657")          # primitive cell
    conv = mpr.get_structure_by_material_id(
        "mp-2657", conventional_unit_cell=True)                # conventional cell
st.to(filename="POSCAR_mp-2657", fmt="poscar")
```

| 关键词/参数 | 推荐值 | 说明 |
|---|---|---|
| `formula` | `"TiO2"` | 按约化化学式检索；也可用 `chemsys="Ti-O"`、`material_ids=[...]` |
| `fields` | 只取需要的字段 | 减少下载量；pymatgen 的 `MPRester` 对应参数名是 `_fields` |
| `energy_above_hull` | 先选 0 或最小 | 0 表示在凸包上（热力学稳定）；同一化学式往往有多个多晶型 |
| `conventional_unit_cell` | 表面/缺陷建模用 `True` | 默认返回原胞；惯用胞更直观，切表面更方便 |

更完整的命令行脚本在 [fetch_mp_structure.py](../examples/03-fetch-convert-si/fetch_mp_structure.py)，同时支持两种客户端：

```bash
cd Common/examples/03-fetch-convert-si
python fetch_mp_structure.py --formula Si                      # list candidates, E_hull ascending
python fetch_mp_structure.py mp-149                            # -> mp-149.cif, POSCAR_mp-149
python fetch_mp_structure.py mp-149 --conventional --prefix Si_conv
python fetch_mp_structure.py mp-2657 --backend pymatgen        # use pymatgen.ext.matproj instead
```

检查与成功标准：

```bash
python fetch_mp_structure.py mp-149 | grep -E "Space group|Sites"     # expect Fd-3m (227), 2 sites
head -7 POSCAR_mp-149                                                  # element line "Si", count "2"
```

> [!WARNING]
> MP 的晶格常数来自 PBE（或 r2SCAN）弛豫：PBE 对 Si 等多数共价/离子晶体会把晶格常数高估约 1%。**正式计算前要用你自己的泛函、ENCUT、POTCAR 重新优化体相**（见 [VASP 第4章](../../VASP/tutorials/04-relax-eos.md)），后续 slab、缺陷、异质结都基于这个优化值。

### 3.3 实验 CIF 的坑：无序、占位与氢原子

从 ICSD/COD 下载的实验 CIF 经常带有：

- **部分占位**（`_atom_site_occupancy < 1`）：合金、固溶体、无序离子。DFT 需要有序结构，pymatgen 写 POSCAR 时会报 `ValueError: Disordered structure with partial occupancies cannot be converted into POSCAR!`。解决：用 `OrderDisorderedStructureTransformation`（小体系枚举）或 SQS（icet、ATAT `mcsqs`）构造有序近似；
- **缺失氢原子**：X 射线结构常不含 H，或 H 位置不准，需要补 H 再优化；
- **对称展开**：CIF 只列不等价原子，读取器根据空间群展开。不同程序展开时的容差不同，偶尔会出现重叠原子，读入后务必检查原子数（`len(atoms)`）是否等于 $Z\times$ 化学式原子数。

```python
# read a CIF with pymatgen (new API: parse_structures replaces get_structures)
from pymatgen.io.cif import CifParser
st = CifParser("TiO2_exp.cif").parse_structures(primitive=False)[0]
print(st.composition, len(st), st.is_ordered)    # is_ordered=False -> partial occupancy present
```

---

## 4. 对称性分析：symprec、原胞与标准惯用胞

### 4.1 为什么要关心对称性

- k 点高对称路径（能带）按**标准原胞**定义，用错晶胞会画出错误的能带（板块辅助脚本 [kpath_gen.py](../scripts/kpath_gen.py) 基于 SeeK-path/spglib 生成标准路径，用 `python ../scripts/kpath_gen.py -h` 查看用法；也可直接用 `seekpath.get_path()` 或 VASPKIT；能带计算见 [VASP 第5章](../../VASP/tutorials/05-band-dos.md)）；
- 原胞原子数最少、计算最便宜；惯用胞便于切表面、扩缺陷超胞；
- 不等价位点（Wyckoff 位置）决定了缺陷、掺杂、吸附**需要算几种构型**；
- VASP 的 `SYMPREC`（默认 1e-5）、`ISYM` 也依赖同样的容差概念，弛豫后对称性“丢失”会导致 k 点数变化。

### 4.2 用 spglib 与 pymatgen 分析

```python
# symmetry_min.py - spglib and pymatgen on the same structure
import spglib
from ase.io import read

atoms = read("Si.cif")
cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)
for symprec in (1e-5, 1e-3, 1e-1):
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)  # spglib >= 2.5: attribute access
    print(symprec, ds.international, ds.number)
prim = spglib.standardize_cell(cell, to_primitive=True, symprec=1e-3)   # (lattice, positions, numbers)
conv = spglib.standardize_cell(cell, to_primitive=False, symprec=1e-3)

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
sga = SpacegroupAnalyzer(Structure.from_file("Si.cif"), symprec=0.01)
ps = sga.get_primitive_standard_structure()       # standard primitive (Setyawan-Curtarolo convention)
cs = sga.get_conventional_standard_structure()
print(sga.get_space_group_symbol(), len(ps), len(cs))
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `symprec`（spglib） | 理想结构 1e-5～1e-3；弛豫后结构 1e-3～1e-2；实验 CIF 1e-2～1e-1 | 单位 Å；spglib 函数默认 1e-5 |
| `symprec`（pymatgen `SpacegroupAnalyzer`） | 默认 0.01 | 另有 `angle_tolerance`（默认 5°） |
| `to_primitive` | 能带、声子：`True`；切面、缺陷：`False` | `standardize_cell` 同时会“理想化”晶格（`no_idealize=False`） |
| `get_primitive_standard_structure()` | 用于 pymatgen 的 `HighSymmKpath` | 与 spglib/SeeK-path 的原胞**取向约定可能不同**，k 路径要和晶胞配套 |

> [!TIP]
> spglib 2.5 起 `get_symmetry_dataset` 返回数据类对象，用 `ds.international`、`ds.wyckoffs` 访问；老代码里的 `ds["international"]` 会给出弃用警告。示例脚本按新写法编写。

完整脚本 [check_symmetry.py](../examples/03-fetch-convert-si/check_symmetry.py) 扫描 5 个 `symprec`，打印原胞/惯用胞和 Wyckoff 位置：

```bash
python check_symmetry.py Si.cif
python check_symmetry.py Si.cif --rattle 0.01     # emulate a noisy CONTCAR
python check_symmetry.py CONTCAR                  # check your relaxed structure
```

成功标准：`Si.cif` 在所有 `symprec` 下都是 `Fd-3m (227)`，原胞 2 原子、|a| = 5.431/√2 ≈ 3.840 Å，惯用胞 8 原子；加 `--rattle 0.01` 后小 `symprec` 降为 `P1`。若你的 CONTCAR 在 `symprec = 1e-3` 下已经降低对称，先确认是不是真实的物理畸变（如 Jahn–Teller、铁电位移），再决定是否对称化。

---

## 5. 格式转换：一个结构，四款软件

### 5.1 ase.io 的常用格式名

| 文件 | `format=` | 读 | 写 | 说明 |
|---|---|---|---|---|
| CIF | `cif` | ✓ | ✓ | 写出时为 P1 全原子列表 |
| POSCAR/CONTCAR | `vasp` | ✓ | ✓ | `direct=True` 分数坐标，`sort=True` 按元素分组 |
| OUTCAR / vasprun.xml | `vasp-out` / `vasp-xml` | ✓ | — | 取最后一帧或全部帧 |
| xyz / extxyz | `xyz` / `extxyz` | ✓ | ✓ | xyz 无晶胞；extxyz 存 `Lattice`、`pbc` 和任意数组 |
| Gaussian 输入 / 输出 | `gaussian-in` / `gaussian-out` | ✓ | ✓/— | 周期体系写 `TV` 行 |
| CP2K | `cp2k-restart`、`cp2k-dcd` | ✓ | — | ASE **不写** CP2K 输入，需自己写 `&CELL/&COORD` 或用 `COORD_FILE_NAME` |
| OpenMX | — | — | — | `ase.io` 不提供 OpenMX 输入格式，块需自己写（或用 ASE 的 OpenMX 计算器） |
| PDB、mol、sdf | `proteindatabank`、`mol`、`sdf` | ✓ | ✓ | 分子（键信息、电荷等不一定保留） |

命令行一行转换：`ase convert Si.cif POSCAR -o vasp`（`-i/-o` 指定输入/输出格式，`-f` 覆盖已有文件）。

### 5.2 一键转换脚本

[cif_to_all.py](../examples/03-fetch-convert-si/cif_to_all.py) 把一个结构写成全部六种文件，并保证 ASCII + LF 换行：

```bash
python cif_to_all.py Si.cif                       # keep the 8-atom cell -> converted/
python cif_to_all.py Si.cif --primitive -o prim   # spglib standard primitive cell
python cif_to_all.py TiO2.cif --valence Ti=12 O=6 # OpenMX valence electrons for new elements
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `--primitive` | 能带/声子用 | 调用 `spglib.standardize_cell(to_primitive=True)` |
| `--symprec` | 1e-3 | 原胞约化的容差 |
| `--valence El=N` | 按 DFT_DATA19 VPS 填 | OpenMX 每个原子的 up/down 初始电子数之和 |
| `-o` | `converted` | 输出目录 |

板块通用脚本 [convert_any.py](../scripts/convert_any.py) 另外支持超胞（`--supercell`）、加真空（`--vacuum`）和按层固定原子（`--fix-layers`、`--fix-below`），用 `python ../scripts/convert_any.py -h` 查看用法；不想引入额外脚本时，超胞、真空和固定层也可按第 6、7 节的几行代码完成。

### 5.3 各软件的结构写法要点

**VASP POSCAR**（格式详见 [VASP 第2章](../../VASP/tutorials/02-input-output-files.md)）：

```text
Si8 converted from Si.cif
 1.0000000000000000
     5.4310000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.4310000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    5.4310000000000000
 Si
  8
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.0000000000000000  0.5000000000000000  0.5000000000000000
  (6 more lines)
```

第 6 行元素顺序必须与 POTCAR 拼接顺序一致；ASE 写出时用 `sort=True` 让同种元素连续，避免出现 `O Ti O` 这样的重复元素行。

**Gaussian**：分子直接写笛卡尔坐标；周期体系在坐标后加 1～3 行 `TV`（平移矢量，单位 Å），Gaussian 据此做 PBC 计算。文件必须以空行结尾（[Gaussian 第2章](../../Gaussian/tutorials/02-input-anatomy.md)）。

**CP2K**：坐标写在 `&FORCE_EVAL/&SUBSYS` 中，可直接内嵌，也可引用外部 xyz（输入结构见 [CP2K 第2章](../../CP2K/tutorials/02-input-structure.md)）：

```text
&SUBSYS
  &CELL
    A [angstrom] 5.431 0.000 0.000
    B [angstrom] 0.000 5.431 0.000
    C [angstrom] 0.000 0.000 5.431
    PERIODIC XYZ
  &END CELL
  &TOPOLOGY
    COORD_FILE_NAME Si.xyz       ! plain xyz, Angstrom
    COORD_FILE_FORMAT XYZ
  &END TOPOLOGY
  &KIND Si
    BASIS_SET DZVP-MOLOPT-SR-GTH
    POTENTIAL GTH-PBE-q4
  &END KIND
&END SUBSYS
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `A/B/C` 或 `ABC` | 晶格矢量 | `A/B/C` 直接给三个矢量（任意晶胞，最不易出错）；`ABC` 给三个边长，非正交晶胞还要配 `ALPHA_BETA_GAMMA`（默认 90 90 90） |
| `PERIODIC` | `XYZ`（体相）、`XY`（表面）、`NONE`（分子） | 非 `XYZ` 时 `&DFT/&POISSON` 的 `PERIODIC`/`PSOLVER` 要配套 |
| `COORD_FILE_FORMAT` | `XYZ` | xyz 不含晶胞，晶胞仍由 `&CELL` 给出 |
| `&KIND` | 每种元素一个 | 基组与 GTH 赝势选择见 [CP2K 第3章](../../CP2K/tutorials/03-gpw-basis-pseudo.md) |

**OpenMX**（[OpenMX 第3章](../../OpenMX/tutorials/03-input-first-calc.md)）：

```text
Atoms.Number  2
Atoms.SpeciesAndCoordinates.Unit  FRAC
<Atoms.SpeciesAndCoordinates
   1  Si   0.000000000  0.000000000  0.000000000   2.00  2.00
   2  Si   0.250000000  0.250000000  0.250000000   2.00  2.00
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit  Ang
<Atoms.UnitVectors
   0.000000000  2.715500000  2.715500000
   2.715500000  0.000000000  2.715500000
   2.715500000  2.715500000  0.000000000
Atoms.UnitVectors>
```

每行依次是序号、物种名、三个坐标、初始 up 与 down 电子数；两者之和必须等于该物种 VPS 的价电子数（Si_PBE19 为 4），物种名必须在 `<Definition.of.Atomic.Species` 中定义。

### 5.4 分子格式：Open Babel 与 Avogadro

分子之间的格式转换用 Open Babel 最快：

```bash
obabel h2o.xyz -O h2o.gjf                 # xyz -> Gaussian input (edit route line afterwards)
obabel -:"c1ccccc1" --gen3d -O benzene.xyz  # SMILES -> 3D
obabel mol.sdf -O mol.xyz
```

转换后用 `file h2o.gjf` 确认是 ASCII，用 `tail -c 2 h2o.gjf | od -c` 确认以空行结尾。

---

## 6. 超胞：make_supercell

### 6.1 对角扩胞与一般变换

```python
# supercell_min.py
import numpy as np
from ase.build import bulk, make_supercell, find_optimal_cell_shape

prim = bulk("Si", "diamond", a=5.431)               # 2-atom fcc primitive cell
sc1 = prim.repeat((3, 3, 3))                        # 54 atoms, still a rhombohedral cell
P = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])  # primitive fcc -> conventional cube
cube = make_supercell(prim, P)                      # 8 atoms, a = 5.431 A
sc2 = make_supercell(prim, 2 * P)                   # 64-atom cube, L = 10.86 A
P_opt = find_optimal_cell_shape(prim.cell, target_size=32, target_shape="sc")
sc3 = make_supercell(prim, P_opt)                   # most cubic 64-atom cell
print(len(sc1), len(cube), len(sc2), len(sc3), sc2.cell.lengths())
sc3.write("POSCAR_sc", format="vasp", direct=True)     # checked in 6.2
```

pymatgen 的等价写法：`st.make_supercell([[2,0,0],[0,2,0],[0,0,2]])`（原地修改）或 `st * (2, 2, 2)`（返回新结构）。

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `P` | 整数矩阵，$\det P>0$ | $\det P<0$ 会得到左手系晶胞，后续程序可能出错；交换矩阵的两行即可变回右手系 |
| `target_size` / `target_shape` | 原胞个数（不是原子数）/ `"sc"`（近立方）或 `"fcc"` | 缺陷、声子超胞希望各方向周期像距离接近；矩阵元搜索范围由 `lower_limit=-2`、`upper_limit=2` 控制 |
| `wrap` | 默认 `True` | 把原子折回晶胞内 |
| 超胞尺寸 | 缺陷/声子 ≥ 10 Å；AIMD ≥ 10～12 Å | 以“最短周期像距离”衡量，不是看原子数 |

### 6.2 如何检查

```bash
python - <<'EOF'
from ase.io import read
from ase.geometry import get_distances
a = read("POSCAR_sc")
d = get_distances(a.positions, cell=a.cell, pbc=True)[1]
d[d == 0] = 99
print(len(a), a.cell.lengths().round(3), "min distance %.3f A" % d.min())
EOF
```

成功标准：原子数 = $|\det P|\times$ 原胞原子数；最短原子间距与体相一致（Si 约 2.35 Å），没有重叠原子；最短晶格矢量长度满足 6.1 节的尺寸要求。

---

## 7. 表面 slab：终止面、真空层与固定底层

### 7.1 用 ASE 搭简单金属表面

[build_cu111_slab.py](../examples/03-slab-adsorption/build_cu111_slab.py) 的核心只有几行：

```python
from ase.build import fcc111
from ase.constraints import FixAtoms

slab = fcc111("Cu", size=(3, 3, 4), a=3.63, vacuum=15.0 / 2, orthogonal=False)
tags = slab.get_tags()                      # 1 = top layer ... 4 = bottom layer
slab.set_constraint(FixAtoms(mask=tags > 2))  # fix the 2 bottom layers
slab.write("POSCAR_Cu111_slab", format="vasp", direct=True)   # writes Selective dynamics
```

```bash
cd Common/examples/03-slab-adsorption
python build_cu111_slab.py                               # defaults: 3x3, 4 layers, 15 A, fix 2
python build_cu111_slab.py --layers 6 --fix 3 --vacuum 20 -o POSCAR_Cu111_6L
grep -c "F   F   F" POSCAR_Cu111_slab                    # 18 fixed atoms
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `a` | 自己优化的体相值 | 示例用 3.63 Å（PBE 量级），实验约 3.615 Å |
| `size` | 面内 ≥ 2×2；吸附 3×3 或 4×4 | 吸附物之间距离 ≥ 7～8 Å，避免横向相互作用 |
| 层数 | 金属 4～6 层，氧化物 3～5 个三明治层 | 用表面能或吸附能对层数收敛 |
| `vacuum` | **每侧**的真空（ASE 约定） | 总真空 = 2 × vacuum；表面计算总真空 ≥ 15 Å |
| 固定层数 | 底部 1/2～1/3 | 只弛豫上层，下层保持体相；对称 slab 则两侧都放开、中间固定 |
| `orthogonal` | `False`（菱形最小胞） | `True` 得到矩形胞，需偶数行 |

> [!WARNING]
> ASE 中 `fcc111(..., vacuum=x)`、`surface(..., vacuum=x)`、`atoms.center(vacuum=x, axis=2)` 的 `x` 都是加在 slab **两侧**的，总真空是 `2x`。很多“真空只有 7.5 Å”的错误就来自这里。

### 7.2 任意晶面：ase.build.surface

```python
from ase.build import bulk, surface
from ase.io import write
conv = bulk("Si", "diamond", a=5.431, cubic=True)      # MUST be the conventional cell
s = surface(conv, (1, 1, 0), layers=6, vacuum=7.5)    # 6 repeat units, 15 A total vacuum
print(len(s), len(set(s.positions[:, 2].round(3))))   # 48 atoms in 12 atomic (110) planes
write("POSCAR_Si110", s, format="vasp", direct=True, sort=True)
```

`surface()` 要求输入是惯用胞，不能是原胞；它只给出一种终止面。注意 `layers` 数的是沿法向的**重复单元**（ASE 文档称 equivalent layers），不一定等于原子层数：上例 Si(110) 每个重复单元含 2 个原子面，`layers=6` 实际是 12 个原子面、48 个原子（面内晶胞 √2a × a ≈ 7.68 × 5.43 Å）。写完后务必用 `len(set(z))` 或 VESTA 数一下原子层。对多元化合物，用 `from ase.build.surfaces_with_termination import surfaces_with_termination`（按指定的顶层元素筛选终止面）或 pymatgen 的 `SlabGenerator` 枚举不同终止面。

### 7.3 化合物表面：pymatgen SlabGenerator

多元化合物（氧化物、钙钛矿）的同一 Miller 面有多种“切法”（shift），得到不同终止面，有的不符合化学计量、有的有垂直偶极（极性面）。[build_slab_TiO2_110.py](../examples/03-slab-adsorption/build_slab_TiO2_110.py) 演示完整筛选：

```python
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator

bulk = Structure.from_spacegroup("P4_2/mnm", Lattice.tetragonal(4.594, 2.959),
                                 ["Ti", "O"], [[0, 0, 0], [0.3048, 0.3048, 0]])
gen = SlabGenerator(bulk, (1, 1, 0), min_slab_size=12.0, min_vacuum_size=15.0,
                    center_slab=True, lll_reduce=True)
slabs = gen.get_slabs()            # one Slab per termination (shift)
good = [s for s in slabs if s.composition.reduced_formula == "TiO2"
        and s.is_symmetric() and not s.is_polar()]
slab = good[0].get_orthogonal_c_slab()   # c perpendicular to the surface
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `min_slab_size` | 10～15 Å（或 `in_unit_planes=True` 时按层数） | 最小厚度，实际 slab 取整到完整层 |
| `min_vacuum_size` | 15 Å | 这是**总**真空（与 ASE 约定不同） |
| `center_slab` | `True` | slab 居中，原子不跨越 c 方向边界 |
| `lll_reduce` | `True` | 面内晶胞约化为更“方”的形状 |
| `primitive` | 默认 `True` | 面内取最小表面胞 |
| `get_slabs(bonds=..., max_broken_bonds=...)` | 如 `{("Ti","O"): 2.2}` | 排除切断过多 Ti–O 键的终止面 |
| `is_symmetric()` / `is_polar()` | `True` / `False` | 上下对称且无净偶极，表面能才有定义 |
| `get_orthogonal_c_slab()` | 用于 VASP 偶极修正 | 让 c 轴垂直于表面 |

```bash
python build_slab_TiO2_110.py                          # 1x1, ~12 A, fix 2 trilayers
python build_slab_TiO2_110.py --thick 16 --repeat 2 1 --fix-trilayers 2
```

成功标准：屏幕列出的候选中至少有一个 `TiO2 symmetric=True polar=False`；输出打印的 `vacuum gap` ≥ 15 Å；底部固定原子数等于所选三明治层的原子数。批量枚举所有低指数面可用 `generate_all_slabs` 与 `get_symmetrically_distinct_miller_indices`。

### 7.4 表面计算的配套设置

| 问题 | VASP | CP2K | OpenMX |
|---|---|---|---|
| 不对称 slab / 吸附导致的偶极 | `LDIPOL = .TRUE.`，`IDIPOL = 3` | `&DFT` 中 `SURFACE_DIPOLE_CORRECTION .TRUE.` | ESM 方法 |
| 法向方向的 k 点 | KPOINTS 第三个数取 1（法向沿 c） | 同左（或只用 Γ 点） | 与法向对应的那个数取 1；用 ESM 时法向必须沿 **a 轴**，即 `scf.Kgrid 1 n n` |
| 固定原子 | POSCAR `Selective dynamics` | `&MOTION/&CONSTRAINT/&FIXED_ATOMS` | `<MD.Fixed.XYZ` |
| 法向约定 | 任意，习惯沿 c | 偶极修正 `SURF_DIP_DIR` 默认 Z | ESM 要求法向沿 a 轴 |

> [!WARNING]
> ASE、pymatgen 生成的 slab 法向都沿 c（z）。交给 OpenMX 做 ESM 计算前，要把晶格矢量和坐标轮换成“法向沿 a（x）”，例如 ASE 中 `slab.set_cell(slab.cell[[2, 0, 1]][:, [2, 0, 1]], scale_atoms=False); slab.positions = slab.positions[:, [2, 0, 1]]`（同时轮换矢量顺序和笛卡尔分量，保持右手系），然后再写 `Atoms.UnitVectors`。前提是原来的 c 严格垂直于表面、沿 z（pymatgen slab 先调用 `get_orthogonal_c_slab()`），轮换后 a 轴才与 x 平行（`ESM.direction` 默认 x）。详见 [OpenMX 第16章](../../OpenMX/tutorials/16-surface-esm-workfunction.md)。

具体计算流程见 [VASP 第9章](../../VASP/tutorials/09-surface-workfunction.md)、[CP2K 第6章](../../CP2K/tutorials/06-periodic-surface-kpoints.md)、[OpenMX 第16章](../../OpenMX/tutorials/16-surface-esm-workfunction.md)。

---

## 8. 吸附位点：AdsorbateSiteFinder 与 add_adsorbate

### 8.1 pymatgen 自动找位点

`AdsorbateSiteFinder` 先按高度识别表面原子，再在其上方构造 ontop（单原子）、bridge（两原子中点）、hollow（三原子中心）位点，并用 slab 的对称性去掉等价位点：

```python
from pymatgen.core import Molecule, Structure
from pymatgen.analysis.adsorption import AdsorbateSiteFinder   # also pymatgen.core.adsorption in new pymatgen-core

slab = Structure.from_file("POSCAR_Cu111_slab")
asf = AdsorbateSiteFinder(slab, selective_dynamics=False)       # keep F/T flags read from POSCAR
sites = asf.find_adsorption_sites(distance=1.85, positions=("ontop", "bridge", "hollow"))
co = Molecule(["C", "O"], [[0, 0, 0], [0, 0, 1.15]])            # C binds to the surface
st = asf.add_adsorbate(co, sites["hollow"][0], translate=False, reorient=True)
st.to(filename="POSCAR_CO_hollow", fmt="poscar")
# or all at once: asf.generate_adsorption_structures(co, repeat=[1, 1, 1], translate=False)
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `distance` | 1.2～2.0 Å | 吸附原子与表面原子平面的初始高度；只是初猜 |
| `selective_dynamics` | 已有 F/T 时用 `False` | `True` 会按“表面/次表面”重新设定固定，覆盖你原来的设置 |
| `symm_reduce`、`near_reduce` | 默认 1e-2 | 对称去重和近邻去重阈值 |
| `translate` | `False` 时以分子坐标原点（这里是 C）落在位点上 | `True` 时以最低原子质心为原点 |
| `reorient` | `True` | 让分子 z 轴沿表面法向 |
| `repeat` | 需要更大面内超胞时设置 | 降低覆盖度 |

### 8.2 ASE 命名位点

ASE 的 `fcc111/fcc100/hcp0001` 等函数生成的 slab 自带命名位点（`ontop`、`bridge`、`fcc`、`hcp` 等），用 `add_adsorbate(slab, "O", 1.2, "fcc")` 即可放置。[adsorption_sites.py](../examples/03-slab-adsorption/adsorption_sites.py) 同时演示两种方式：

```bash
python adsorption_sites.py                          # pymatgen part + ASE part
python adsorption_sites.py --skip-pymatgen          # ASE only
ls ads_*.vasp | wc -l
```

### 8.3 吸附能所需的结构

$$
E_{\text{ads}} = E_{\text{slab+ads}} - E_{\text{slab}} - E_{\text{ads,gas}}
$$

三个结构必须用**相同的晶胞、ENCUT、k 点、泛函**；气相分子放在同样大小（或足够大）的盒子里。成功标准：吸附物在表面上方而非穿过 c 方向边界出现在 slab 底部；吸附物 `T T T`，slab 底层仍为 `F F F`；fcc 与 hcp 两种 hollow 都算到。吸附自由能与 CHE 模型见 [VASP 第10章](../../VASP/tutorials/10-adsorption-catalysis.md)。

---

## 9. 缺陷与掺杂：空位、替位、间隙

### 9.1 基本操作

| 缺陷 | ASE 写法 | pymatgen 写法 |
|---|---|---|
| 空位 V_X | `del atoms[i]` | `st.remove_sites([i])` |
| 替位 A_B | `atoms[i].symbol = "Al"` | `st.replace(i, "Al")` 或 `st[i] = "Al"` |
| 间隙 X_i | `atoms.append("Li"); atoms.positions[-1] = r` | `st.append("Li", frac, coords_are_cartesian=False)` |
| 不等价位点 | `spglib.get_symmetry_dataset(...).equivalent_atoms` | `SpacegroupAnalyzer(st).get_symmetrized_structure().equivalent_indices` |
| 构型去重 | — | `StructureMatcher().group_structures(list)` |

### 9.2 MgO 示例

[make_defects.py](../examples/03-defect-doping-mgo/make_defects.py) 生成完美超胞、V_O、Al_Mg、Li_i 和对称不等价的 Al–Al 双掺杂构型：

```bash
cd Common/examples/03-defect-doping-mgo
python make_defects.py                   # 2x2x2 -> 64 atoms, L = 8.42 A (practice only)
python make_defects.py --n 3 --a 4.25    # 216 atoms, L = 12.75 A (use YOUR optimized a)
for f in MgO_*.vasp; do echo "$f: $(sed -n 6p $f) | $(sed -n 7p $f)"; done
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `--a` | 自己优化的晶格常数 | 默认 4.212 Å 接近室温实验值 |
| `--n` | 使 L ≥ 10 Å，并做 n 的收敛 | 2 → 8.4 Å 只适合练习 |
| `--max-pairs` | 4～6 | 列出最近的若干 Al–Al 距离壳层 |

成功标准：完美超胞空间群 `Fm-3m (225)`，只有 Mg、O 两组不等价位点（所以 V_O 只需算一种）；`MgO_VO.vasp` 为 Mg32O31；Li_i 与最近邻距离 √3/4·a ≈ 1.82 Å；双掺杂构型按 Al–Al 距离从 a/√2 ≈ 2.98 Å 起递增。

### 9.3 缺陷形成能与带电缺陷

$$
E_f[X^q] = E_{\text{tot}}[X^q] - E_{\text{tot}}[\text{bulk}] - \sum_i n_i\mu_i + q\,(E_{\text{VBM}} + E_F) + E_{\text{corr}}
$$

$n_i$ 为加入（+）或移走（−）的第 $i$ 种原子数，$\mu_i$ 为化学势，$E_{\text{corr}}$ 为带电超胞的有限尺寸修正（Freysoldt、Kumagai 等方案）。建模层面要注意：

- 带电缺陷在 VASP 中通过改 `NELECT` 实现，超胞仍是中性原子排布；
- 缺陷附近常有对称性破缺（极化子、Jahn–Teller），从完美对称位置出发可能卡在亚稳态。对初始结构做小幅随机位移（rattle）或用 [ShakeNBreak](https://github.com/SMTG-Bham/ShakeNBreak) 系统搜索；
- 系统性的缺陷计算（生成、电荷态、修正、转变能级）推荐 [doped](https://github.com/SMTG-Bham/doped)，pymatgen 的缺陷模块已拆分为独立包 `pymatgen-analysis-defects`。

### 9.4 随机掺杂与合金

掺杂浓度较高或要模拟无序合金时，单个手工构型不具有代表性：

- 仓库作者的 [replace_dopant_symbols.py](../../Scripts/modeling/replace_dopant_symbols.py)：读入 `POSCAR`，扩胞后按数量随机替换目标元素，可一次生成多个随机构型（输出 extxyz），也可按“每个原胞拷贝替换数相同”的方式分布掺杂。参数在文件头修改，详见 [Scripts/README.md](../../Scripts/README.md) 第 5.1 节；
- SQS（special quasirandom structures）：用 icet 或 ATAT `mcsqs` 生成在短程关联上最接近随机合金的超胞；
- 少量取代的全部不等价构型：pymatgen `OrderDisorderedStructureTransformation` 或 enumlib 枚举，再用 `StructureMatcher` 去重。

---

## 10. 二维材料与范德华异质结

### 10.1 晶格匹配

两层六方材料的公共超胞：底层取超胞矢量 $\mathbf{v}_1=i_1\mathbf{a}_1+j_1\mathbf{a}_2$，顶层取 $\mathbf{v}_2$，当 $|\mathbf{v}_1|\approx|\mathbf{v}_2|$ 时匹配。在 $\gamma=120^\circ$ 约定下 $|\mathbf{v}|=a\sqrt{i^2-ij+j^2}$，超胞矩阵 $\begin{pmatrix}i&j\\-j&i-j\end{pmatrix}$ 含 $n=i^2-ij+j^2$ 个原胞（60° 约定下对应 $i^2+ij+j^2$，见 2.2 节）。[build_heterostructure.py](../examples/03-heterostructure-2d/build_heterostructure.py) 枚举所有组合，按原子数排序：

```bash
cd Common/examples/03-heterostructure-2d
python build_heterostructure.py --pair graphene-MoS2 --list           # candidates only
python build_heterostructure.py --pair graphene-MoS2 --pick 0         # build the smallest one
python build_heterostructure.py --pair MoS2-WS2                       # 1x1 on 1x1, +0.32 % strain
sed -n 6,7p POSCAR_graphene_MoS2
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `--tol` | 0.01～0.03 | 上层最大应变；越小超胞越大 |
| `--natoms-max` | 按算力，100～300 | 超过就不列出 |
| `--dist` | 3.3～3.5 Å | 初始层间距，之后必须用 vdW 修正弛豫 |
| `--vacuum` | ≥ 15～20 Å | c 方向总真空 |

成功标准：候选表中应变和扭转角符合预期；生成的结构 γ = 120°，两层完整铺满晶胞，层间距约为 `--dist`；正式计算前比较不同堆垛（0°/60° 扭转、AA/AB）的能量。

### 10.2 其他工具与注意事项

- 非六方晶格、体相/体相界面：pymatgen `ZSLGenerator` 与 `CoherentInterfaceBuilder`（`pymatgen.analysis.interfaces`），算法来自 Zur–McGill；
- 单层结构可从 [C2DB](https://c2db.fysik.dtu.dk/) 下载，或用 ASE 的 `graphene()`、`mx2()` 直接生成；
- 应变会改变带隙和能带对齐，论文中要注明加在哪一层、多大；
- vdW 修正（DFT-D3、optB86b-vdW 等）和二维材料的 k 点、真空、偶极设置见 [VASP 第8章](../../VASP/tutorials/08-vdw-2d.md)；
- 需要弯曲的碳片段（研究曲率效应）时，可用仓库作者的 [fabricate_spherical_carbon.py](../../Scripts/modeling/fabricate_spherical_carbon.py)：在平面六角点阵上截取圆盘，再投影到给定半径的球面，输出 XYZ 和 CIF（[Scripts/README.md](../../Scripts/README.md) 第 5.2 节）。

---

## 11. 分子建模：Avogadro、RDKit 与电荷/多重度

### 11.1 三种起点

| 方式 | 适合 | 说明 |
|---|---|---|
| [Avogadro 2](https://two.avogadro.cc/) 手绘 | 小分子、配合物、需要直观调整的结构 | 免费开源，可用内置力场预优化，并生成 Gaussian、CP2K 等程序的输入模板 |
| RDKit SMILES → 3D | 有机分子、批量处理 | ETKDG 构象生成 + MMFF94 预优化，可脚本化 |
| PubChem/CSD 下载 | 已知化合物 | 3D SDF 用 Open Babel 或 ASE 转 xyz |

### 11.2 RDKit 脚本

[smiles_to_3d.py](../examples/03-molecules/smiles_to_3d.py) 生成多个构象、力场优化后保留最低能量者，并把电荷和多重度写进 xyz 注释行：

```bash
cd Common/examples/03-molecules
python smiles_to_3d.py "c1ccccc1" -o benzene_rdkit.xyz
python smiles_to_3d.py "CC(=O)Nc1ccc(O)cc1" -n 30 -o paracetamol.xyz
python smiles_to_3d.py "[NH4+]" -o nh4.xyz            # charge +1 read from SMILES
head -2 paracetamol.xyz                                 # SMILES="..." charge=0 mult=1 MMFF94_E=...
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `-n` | 刚性分子 10，柔性分子 50～300 | 构象数；更系统的搜索用 CREST 或 aqme |
| `--seed` | 固定值 | 保证可复现 |
| 力场 | MMFF94，缺参数时自动换 UFF | 只是 DFT 优化的起点 |

### 11.3 一个分子 → 四款软件

[xyz_to_inputs.py](../examples/03-molecules/xyz_to_inputs.py) 在写文件之前检查电荷/多重度与电子数奇偶是否一致，再写出 `.gjf`、CP2K 孤立分子片段（`PERIODIC NONE`）、VASP 立方盒子 POSCAR 和 OpenMX 块：

```bash
python xyz_to_inputs.py h2o.xyz                     # 10 electrons, singlet: OK
python xyz_to_inputs.py h2o.xyz --mult 2            # stops: impossible multiplicity
python xyz_to_inputs.py paracetamol.xyz --task "opt freq" --nproc 16
grep -A2 "&CELL" h2o_cp2k_subsys.inc
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `--charge` / `--mult` | 自己确认 | 电子数为偶数时多重度只能是奇数，反之亦然 |
| `--padding` | 7.5～10 Å | 盒子每侧真空；VASP 盒子边长 = 分子尺寸 + 2 × padding |
| `--method/--basis/--extra` | 如 `B3LYP`、`def2SVP`、`EmpiricalDispersion=GD3BJ` | 只是模板，按 [Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md) 选择 |

孤立分子在周期性程序中的设置：VASP 只用 Γ 点（`vasp_gam`），极性分子加 `LDIPOL = .TRUE.`、`IDIPOL = 4`；CP2K 用 `PERIODIC NONE` 并在 `&POISSON` 中选 `PSOLVER MT`（盒子约为分子尺寸的 2 倍）或 `WAVELET`；OpenMX 用 `scf.Kgrid 1 1 1`。仓库自带的 [h2o.xyz](../examples/03-molecules/h2o.xyz)（实验气相几何）和 [benzene.gjf](../examples/03-molecules/benzene.gjf)（D6h 起始结构，opt + freq）可直接作为练习输入。

---

## 12. 建模后检查清单：VESTA、ASE 与几行脚本

提交计算前花两分钟检查，能省下几天排错时间：

| 检查项 | 方法 | 合格标准 |
|---|---|---|
| 原子数与组成 | `sed -n 6,7p POSCAR`；`atoms.get_chemical_formula()` | 与设计一致，元素顺序与 POTCAR 一致 |
| 最短原子间距 | 6.2 节脚本；VESTA 距离工具 | 无 < 0.8 Å 的重叠；键长与常见值相差 < 10% |
| 对称性 | `check_symmetry.py` | 体相恢复期望空间群；slab/缺陷对称性降低在预期之内 |
| 真空 | `c - (zmax - zmin)` | 表面、二维 ≥ 15 Å；分子每侧 ≥ 7 Å |
| 固定原子 | `grep -c "F   F   F" POSCAR` | 数目、位置（底层）正确 |
| 原子跨边界 | VESTA 中看 slab 是否被 c 方向边界切开 | slab 连续、居中 |
| 换行与编码 | `file POSCAR` | `ASCII text`，没有 `with CRLF line terminators` |

**VESTA**（稳定版 3.5.8，开发版 3.90.x 提供 Apple Silicon、Windows on Arm 原生构建）：直接打开 CIF、POSCAR、xyz；用工具栏的距离/角度工具测键长键角，在 Bonds 对话框里设定成键截断，异常长或短的键会一目了然；用 `File → Export Data` 可另存为 VASP 等格式。**ASE GUI**（`ase gui POSCAR`）适合快速查看，还能显示约束原子和轨迹。

---

## 参数选择与收敛建议

| 项目 | 起步值 | 收敛/检查方法 |
|---|---|---|
| 体相晶格常数 | 用自己的泛函优化（不用 MP/实验值直接建 slab） | E–V 拟合或 ISIF=3 弛豫，见 [VASP 第4章](../../VASP/tutorials/04-relax-eos.md) |
| `symprec` | 理想 1e-5～1e-3；CONTCAR 1e-3～1e-2；实验 CIF 1e-2～1e-1 | 扫描几个值，空间群稳定的区间即为合理取值 |
| 缺陷超胞 | 最短周期像距离 ≥ 10 Å | 至少两个尺寸比较形成能，变化 < 0.05～0.1 eV |
| slab 层数 | 金属 4～6 层；氧化物 3～5 个重复单元 | 表面能变化 < 0.01 J/m²（或吸附能 < 0.02～0.05 eV） |
| 真空 | ≥ 15 Å（总） | 功函数/静电势在真空区平台平坦 |
| 固定层 | 底部 1/2～1/3 | 放开一层后吸附能变化很小 |
| 吸附超胞 | 吸附物间距 ≥ 7～8 Å | 不同覆盖度对比 |
| 异质结应变 | < 1～2% | 比较两个候选超胞的能带对齐/结合能 |
| 层间距 | 3.3～3.5 Å 起步 | vdW 修正下弛豫，比较不同堆垛 |
| 分子盒子 | 每侧 ≥ 7～8 Å | 盒子增大 2 Å，能量变化 < 1 meV |
| 构象 | 力场多构象筛选 | 低能构象再用 DFT 优化比较 |

---

## 常见问题与报错

| 现象 / 报错原文 | 原因 | 解决 |
|---|---|---|
| `ValueError: Invalid or old API key. Please obtain an updated API key at https://materialsproject.org/dashboard.`（pymatgen `MPRester`） | 使用了旧的 16 位 key，或没有设置 `PMG_MAPI_KEY` | 到 dashboard 复制 32 位新 key，`export PMG_MAPI_KEY=...` 或 `MP_API_KEY=...` |
| mp-api 请求返回 401/403、连接错误 | 旧 key、旧地址 `rest/v2`（2025-09 下线） | `pip install -U mp-api`，删除 `MP_API_ENDPOINT` 环境变量，换新 key |
| `ImportError: cannot import name 'Structure' from 'pymatgen'` | 旧代码写法 | 改为 `from pymatgen.core import Structure` |
| `ValueError: Disordered structure with partial occupancies cannot be converted into POSCAR!` | CIF 含 `_atom_site_occupancy < 1` | 构造有序结构（枚举、SQS），见 3.3 节 |
| OpenMX ESM 计算几何检查报错或结果异常，slab 由 ASE/pymatgen 生成 | slab 法向沿 c，而 ESM 要求法向沿 a 轴 | 按 7.4 节轮换晶格矢量与坐标，`scf.Kgrid` 改为 `1 n n` |
| 读入 CIF 后原子数比预期多，有两个原子几乎重合 | 坐标位数太少（如 0.333 代替 1/3），对称展开时等价位点没有被合并 | 换用 pymatgen `CifParser`（`site_tolerance`、`frac_tolerance` 控制合并，默认 1e-4，可适当放宽），或在 VESTA 中打开后导出 P1 结构；检查 `len(atoms)` 与最短距离 |
| VASP 输出 `The distance between some ions is very small. Please check the nearest-neighbor list in the OUTCAR file.` | 原子重叠：重复原子、吸附物穿过边界、超胞变换出错 | 用 6.2 节脚本找出最短距离原子，修正结构 |
| VASP 结果离谱但不报错，OUTCAR 中 `POTCAR:` 顺序与 POSCAR 第 6 行不一致 | POSCAR 元素顺序与 POTCAR 拼接顺序不同 | 写 POSCAR 时 `sort=True`，按第 6 行顺序拼 POTCAR |
| slab 计算功函数异常、静电势在真空区不平 | 真空太薄（ASE `vacuum` 只给了一侧值）或未做偶极修正 | 总真空 ≥ 15 Å；`LDIPOL/IDIPOL = 3` |
| 表面能随层数振荡、两侧表面不同 | slab 不对称或非化学计量、极性面 | `SlabGenerator` 中选 `is_symmetric()` 且 `not is_polar()` 的终止面 |
| 弛豫后 spglib 给出 `P1` | `symprec` 太小，或真实畸变 | 放宽到 1e-3～1e-2 再判断；确认是否物理畸变 |
| Gaussian：`The combination of multiplicity 2 and 10 electrons is impossible.` | 电荷/多重度与电子数奇偶矛盾 | 用 `xyz_to_inputs.py` 预检；修改电荷或多重度 |
| Gaussian：`Small interatomic distances encountered` 后 `Error termination` | 坐标中有重叠原子，或单位搞错（Bohr 当 Å） | 检查结构；`.xyz` 默认单位是 Å |
| `sbatch`/VASP 读不对文件，`file POSCAR` 显示 `with CRLF line terminators` | Windows 编辑或保存 | `dos2unix POSCAR`；示例脚本都以 LF 写文件 |
| `smiles_to_3d.py` 报 `Embedding failed; try more conformers or another seed` | RDKit `EmbedMultipleConfs` 返回空列表，构象嵌入失败（大环、笼状、金属配合物） | 增大 `-n`、换 `--seed`；配合物改用 Avogadro 手工搭建 |
| OpenMX 计算电子数不对或磁矩异常 | `Atoms.SpeciesAndCoordinates` 中 up+down ≠ VPS 价电子数 | 按 DFT_DATA19 中 VPS 的价电子数填写（脚本用 `--valence` 指定） |

---

## 练习

1. 用 `fetch_mp_structure.py --formula TiO2` 列出 MP 中的全部 TiO2 条目，找出金红石和锐钛矿的 mp-id 与 E_hull，分别下载惯用胞，用 `check_symmetry.py` 确认空间群（P4₂/mnm 与 I4₁/amd）。
2. 对 `Si.cif` 加 `--rattle 0.005/0.02/0.05`，记录每种情况下能恢复 `Fd-3m` 的最小 `symprec`，总结位移幅度与 `symprec` 的关系。
3. 用 6.1 节的方法，从 Si 原胞出发分别生成 54 原子的 3×3×3 菱面体超胞和 64 原子的立方超胞，计算两者的最短周期像距离，说明缺陷计算应该选哪个。
4. 用 `build_cu111_slab.py` 生成 4、5、6 层 Cu(111)，再用 `adsorption_sites.py` 在每个 slab 上放 O@fcc，为“吸附能对层数收敛”的测试准备全部输入（不必计算）。
5. 修改 `build_slab_TiO2_110.py`，打印 `get_slabs()` 找到的所有终止面中哪些是极性的，并用 VESTA 观察极性与非极性终止面的差别。
6. 用 `make_defects.py --n 3` 生成 MgO 双 Al 掺杂构型，统计 `StructureMatcher` 去重前后的构型数。
7. 用 `build_heterostructure.py --pair graphene-hBN --list` 比较 `--tol 0.01` 与 `--tol 0.03` 时最小候选的原子数，体会“应变—原子数”的权衡。
8. 用 `smiles_to_3d.py` 生成乙醇和甘氨酸（中性与两性离子）的结构，再用 `xyz_to_inputs.py` 写出 Gaussian 输入，检查电荷和多重度。

---

## 参考资料

官方文档：

- Materials Project API 入门（新 key、mp-api）：<https://docs.materialsproject.org/downloading-data/using-the-api/getting-started>
- ASE 文档（2026-06 起迁到 docs.ase-lib.org）：<https://docs.ase-lib.org/>；文件读写 `ase.io`：<https://docs.ase-lib.org/ase/io/io.html>；表面构建：<https://docs.ase-lib.org/ase/build/surface.html>；版本说明：<https://docs.ase-lib.org/releasenotes.html>
- pymatgen 文档：<https://pymatgen.org/>；`pymatgen.core`（Structure、surface）：<https://pymatgen.org/pymatgen.core.html>；`pymatgen.analysis`（adsorption、interfaces、structure_matcher）：<https://pymatgen.org/pymatgen.analysis.html>；更新日志：<https://pymatgen.org/CHANGES.html>
- spglib 文档：<https://spglib.readthedocs.io/>；SeeK-path：<https://seekpath.readthedocs.io/>
- RDKit 文档：<https://www.rdkit.org/docs/>；Open Babel 文档：<https://openbabel.org/docs/>；Avogadro 2：<https://two.avogadro.cc/>
- VESTA 下载（含 3.90.x 开发版）：<https://jp-minerals.org/vesta/en/download.html>
- VASP Wiki POSCAR：<https://www.vasp.at/wiki/index.php/POSCAR>；LDIPOL：<https://www.vasp.at/wiki/index.php/LDIPOL>
- Gaussian 周期边界条件（PBC）：<https://gaussian.com/pbc/>
- CP2K 输入手册 `&SUBSYS/&TOPOLOGY`：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/TOPOLOGY.html>
- OpenMX 官网与手册：<https://www.openmx-square.org/>
- 数据库：COD <https://www.crystallography.net/cod/>、C2DB <https://c2db.fysik.dtu.dk/>、OPTIMADE <https://www.optimade.org/>、PubChem <https://pubchem.ncbi.nlm.nih.gov/>

相关 GitHub 仓库（只链接，不复制代码）：

- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — Structure、SlabGenerator、AdsorbateSiteFinder、界面与 MP 客户端；2026 年核心拆分为 pymatgen-core。
- [materialsproject/api](https://github.com/materialsproject/api) — mp-api 客户端，旧 API 已于 2025-09 下线。
- [ase/ase（GitLab）](https://gitlab.com/ase/ase) — `ase.io`、`ase.build`（surface、make_supercell、graphene、mx2）。
- [spglib/spglib](https://github.com/spglib/spglib) — 空间群识别、原胞/标准胞。
- [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath) — 标准原胞与高对称 k 路径。
- [SMTG-Bham/doped](https://github.com/SMTG-Bham/doped) — 点缺陷超胞、电荷态、有限尺寸修正一体化工具。
- [OpenChemistry/avogadroapp](https://github.com/OpenChemistry/avogadroapp) — 开源分子编辑器，GaussView 的免费替代。
- [jvalegre/aqme](https://github.com/jvalegre/aqme) — 从 SMILES 出发的构象搜索与 QM 输入批量生成。
- [hsiaoyi0504/awesome-cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics) — RDKit、Open Babel 等化学信息学工具列表。
- [computron/pymatgen_tutorials](https://github.com/computron/pymatgen_tutorials) — pymatgen 作者的结构、对称性、MP API 视频与 notebook。
- [materialsproject/workshop](https://github.com/materialsproject/workshop) — MP 官方 workshop（pymatgen、MP API）。
- [janosh/matterviz](https://github.com/janosh/matterviz) — 浏览器中交互查看晶体结构，VESTA 的轻量替代之一。

关键文献：

- A. H. Larsen et al., The atomic simulation environment—a Python library for working with atoms, *J. Phys.: Condens. Matter* **29**, 273002 (2017). DOI: [10.1088/1361-648X/aa680e](https://doi.org/10.1088/1361-648X/aa680e)
- S. P. Ong et al., Python Materials Genomics (pymatgen), *Comput. Mater. Sci.* **68**, 314 (2013). DOI: [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)
- A. Jain et al., The Materials Project, *APL Mater.* **1**, 011002 (2013). DOI: [10.1063/1.4812323](https://doi.org/10.1063/1.4812323)
- A. Togo, K. Shinohara, I. Tanaka, Spglib: a software library for crystal symmetry search, *Sci. Technol. Adv. Mater.: Methods* **4**, 2384822 (2024). DOI: [10.1080/27660400.2024.2384822](https://doi.org/10.1080/27660400.2024.2384822)
- Y. Hinuma et al., Band structure diagram paths based on crystallography, *Comput. Mater. Sci.* **128**, 140 (2017). DOI: [10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015)
- R. Tran et al., Surface energies of elemental crystals, *Sci. Data* **3**, 160080 (2016). DOI: [10.1038/sdata.2016.80](https://doi.org/10.1038/sdata.2016.80)
- J. H. Montoya, K. A. Persson, A high-throughput framework for determining adsorption energies on solid surfaces, *npj Comput. Mater.* **3**, 14 (2017). DOI: [10.1038/s41524-017-0017-z](https://doi.org/10.1038/s41524-017-0017-z)
- A. Zur, T. C. McGill, Lattice match: An application to heteroepitaxy, *J. Appl. Phys.* **55**, 378 (1984). DOI: [10.1063/1.333084](https://doi.org/10.1063/1.333084)
- C. Freysoldt et al., First-principles calculations for point defects in solids, *Rev. Mod. Phys.* **86**, 253 (2014). DOI: [10.1103/RevModPhys.86.253](https://doi.org/10.1103/RevModPhys.86.253)
- A. Zunger, S.-H. Wei, L. G. Ferreira, J. E. Bernard, Special quasirandom structures, *Phys. Rev. Lett.* **65**, 353 (1990). DOI: [10.1103/PhysRevLett.65.353](https://doi.org/10.1103/PhysRevLett.65.353)
- S. R. Kavanagh et al., doped: Python toolkit for robust and repeatable charged defect supercell calculations, *J. Open Source Softw.* **9**, 6433 (2024). DOI: [10.21105/joss.06433](https://doi.org/10.21105/joss.06433)
- S. Riniker, G. A. Landrum, Better informed distance geometry: using what we know to improve conformation generation, *J. Chem. Inf. Model.* **55**, 2562 (2015). DOI: [10.1021/acs.jcim.5b00654](https://doi.org/10.1021/acs.jcim.5b00654)
- M. D. Hanwell et al., Avogadro: an advanced semantic chemical editor, visualization, and analysis platform, *J. Cheminform.* **4**, 17 (2012). DOI: [10.1186/1758-2946-4-17](https://doi.org/10.1186/1758-2946-4-17)
- N. M. O'Boyle et al., Open Babel: An open chemical toolbox, *J. Cheminform.* **3**, 33 (2011). DOI: [10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)
- K. Momma, F. Izumi, VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data, *J. Appl. Cryst.* **44**, 1272 (2011). DOI: [10.1107/S0021889811038970](https://doi.org/10.1107/S0021889811038970)

[← 上一章](02-job-scheduler.md) | [返回目录](../README.md) | [下一章 →](04-ase-pymatgen.md)
