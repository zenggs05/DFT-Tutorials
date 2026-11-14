# 第4章 ASE 与 pymatgen 实战：统一驱动四款软件、输入集与输出解析（ASE & pymatgen in Practice: Driving VASP, Gaussian, CP2K and OpenMX, Input Sets and Output Parsing）

本章用 Python 把四款软件串起来。ASE 负责用同一个 `Atoms` 对象、同一套优化器（BFGS/FIRE）、EOS 拟合和 NEB（IDPP 插值）驱动 VASP、Gaussian、CP2K、OpenMX；pymatgen 负责用 `MPRelaxSet`、`MPStaticSet`、`MPNonSCFSet`、`MPHSEBSSet`、`MP24RelaxSet`、`MatPESStaticSet` 等输入集生成 VASP 输入，用 `Vasprun`、`Outcar`、`Oszicar`、`Chgcar`、`Locpot` 解析输出，再结合 Materials Project 新版 API 做相图和 $E_\text{hull}$。最后用 ASE、cclib、dpdata 跨程序读取能量、力和结构，批量检查作业状态，并把结果归档到 ase db、CSV/HDF5，上传 NOMAD 或通过 OPTIMADE 检索（Atoms, Calculator, optimizers, NEB, ase db, pymatgen input sets, Vasprun parsing, phase diagram, cclib, dpdata, job status, data archiving）。

> [!NOTE]
> 本章示例分四个目录，每个目录的 README 都写明了文件说明、运行步骤和检查方法：
> [04-ase-emt-basics](../examples/04-ase-emt-basics/README.md)（不装 DFT 软件也能跑的 ASE 练习）、
> [04-ase-calculators](../examples/04-ase-calculators/README.md)（ASE 驱动四款软件）、
> [04-pymatgen-io](../examples/04-pymatgen-io/README.md)（pymatgen 输入集、解析、相图、Gaussian/CP2K I/O、体数据运算、OPTIMADE）、
> [04-output-parsing](../examples/04-output-parsing/README.md)（跨程序解析、状态扫描、归档、MLIP 数据）。
> 本章内容基于 ASE 3.29.0、pymatgen 2026.10.x（pymatgen-core 2026.10.x）、cclib 1.8.1、dpdata 1.1.0，目标程序为 VASP 6.4–6.6、Gaussian 16 C.01/C.02、CP2K 2026.2、OpenMX 4.0。Python 环境的搭建见 [第1章](01-linux-hpc.md)。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 结构操作 | 用 ASE `Atoms` 建模、加 `FixAtoms` 约束，用 `ase.io.read/write` 在 POSCAR、CIF、xyz、extxyz、gjf、traj 之间转换 |
| 统一驱动 | 用 `Vasp`、`Gaussian`、`CP2K`、`OpenMX` 四个计算器对同一结构算能量和力，并说清每个计算器的环境变量和单位约定 |
| 优化与拟合 | 用 BFGS/FIRE 做约束优化，理解 `fmax`；用 `EquationOfState` 拟合晶格常数和体模量 |
| 过渡态 | 用 `ase.mep.NEB` + IDPP 插值 + 爬坡图像做 NEB，显式指定 `method`（ASE 3.28 起默认改为 `improvedtangent`） |
| VASP 输入集 | 用 pymatgen 输入集生成 relax→static→band/DOS→HSE 一整套输入，并用 `user_incar_settings` 覆盖默认值 |
| 输出解析 | 用 `Vasprun`/`Outcar`/`Oszicar`/`Chgcar`/`Locpot` 提取能量、带隙、磁矩、功函数，画能带和 DOS |
| 热力学稳定性 | 用 `PhaseDiagram` 计算形成能、$E_\text{hull}$ 和分解产物，用 `MPRester` 查询 summary |
| 跨程序 | 用 pymatgen 写 Gaussian/CP2K 输入；用 ASE、cclib、dpdata 统一读取四款软件的输出 |
| 数据管理 | 批量判断作业状态；归档到 ase db/CSV/HDF5；了解 NOMAD 上传与 OPTIMADE 查询；知道哪些文件不能公开 |

适用场景与工具选择：

| 场景 | 推荐 | 原因 |
|---|---|---|
| 单个体系、需要精细控制每个关键词 | 手写输入（见各软件章节） | 透明，便于排错 |
| 同一流程换不同程序/不同体系重复做 | ASE 计算器 | 一套 Python 代码，改一行换程序 |
| VASP 高通量、要与 Materials Project 数据对比 | pymatgen 输入集 | 参数与 MP 一致，可直接进相图 |
| 分子体系的 Gaussian 输出分析 | cclib、pymatgen `GaussianOutput` | 频率、轨道能、几何一次读出 |
| 机器学习势训练集 | dpdata、ASE extxyz | 能量、力、维里统一格式 |
| 数百个作业的状态检查与归档 | 状态扫描脚本 + ase db | 标准库即可运行，结果可查询 |

> [!TIP]
> ASE 与 pymatgen 不冲突：`pymatgen.io.ase.AseAtomsAdaptor` 可以在 `Atoms` 与 `Structure` 之间互转。常见搭配是用 pymatgen 建模和分析对称性，用 ASE 驱动计算和优化。

---

## 2. 原理速览

### 2.1 ASE 的三层结构

```text
 Atoms（结构：symbols, positions[Å], cell[Å], pbc, constraints, info）
   |  atoms.calc = 某个 Calculator
   v
 Calculator（Vasp / Gaussian / CP2K / OpenMX / EMT / MACE ...）
   |  写输入 -> 调用外部程序 -> 读输出，返回 energy[eV], forces[eV/Å], stress[eV/Å^3]
   v
 Dynamics / Optimizer / NEB / EOS / ase.db（只依赖能量和力，与具体程序无关）
```

ASE 内部单位统一为 eV、Å、eV/Å，计算器负责与程序的原生单位互相换算（Gaussian、CP2K 内部用 Hartree/Bohr，OpenMX 输入用 Ry/Hartree）。常用换算：

$$1\ \text{Ha} = 27.211386\ \text{eV},\quad 1\ \text{Ry} = 13.605693\ \text{eV},\quad 1\ \text{eV/Å}^3 = 160.2177\ \text{GPa}$$

### 2.2 优化收敛判据

ASE 优化器在最大原子力小于阈值时停止：

$$\max_i \lVert \mathbf{F}_i \rVert < f_\text{max}$$

被 `FixAtoms` 固定的原子不参与判断。BFGS 用近似 Hessian 做拟牛顿步，近平衡结构收敛快；FIRE 是带“惯性”的分子动力学式最速下降，初始结构很差或势能面很平时更稳。VASP 中对应的判据是 `EDIFFG < 0`（同样按力判断），Gaussian 用最大力、均方根力和位移四项判据。

### 2.3 状态方程（EOS）

对一组 $(V, E)$ 数据拟合 Birch–Murnaghan 方程：

$$E(V) = E_0 + \frac{9V_0B_0}{16}\left\{\left[\left(\tfrac{V_0}{V}\right)^{2/3}-1\right]^3 B_0' + \left[\left(\tfrac{V_0}{V}\right)^{2/3}-1\right]^2\left[6-4\left(\tfrac{V_0}{V}\right)^{2/3}\right]\right\}$$

得到平衡体积 $V_0$、体模量 $B_0$ 及其压力导数 $B_0'$。ASE 的 `EquationOfState.fit()` 返回的 $B$ 单位是 eV/Å³，乘以 160.2177（或除以 `ase.units.kJ` 再乘 $10^{24}$）得到 GPa。

### 2.4 NEB 与爬坡图像

NEB 在初末态之间放 $N$ 个图像，每个图像受到“垂直于路径的真实力 + 沿路径的弹簧力”：

$$\mathbf{F}_i^\text{NEB} = -\nabla E(\mathbf{R}_i)\big|_\perp + \mathbf{F}_i^\text{s}\big|_\parallel ,\qquad \mathbf{F}_i^\text{s}\big|_\parallel = k\left(\lvert\mathbf{R}_{i+1}-\mathbf{R}_i\rvert - \lvert\mathbf{R}_i-\mathbf{R}_{i-1}\rvert\right)\hat{\boldsymbol\tau}_i$$

切线 $\hat{\boldsymbol\tau}_i$ 的取法就是 ASE `method` 参数的区别：`improvedtangent` 是 Henkelman–Jónsson 的改进切线，`aseneb` 是 ASE 早期的自有实现。爬坡图像（CI）去掉能量最高图像的弹簧力并把沿切线的真实力反向：

$$\mathbf{F}_\text{CI} = -\nabla E + 2\left(\nabla E\cdot\hat{\boldsymbol\tau}\right)\hat{\boldsymbol\tau}$$

IDPP（image dependent pair potential）插值让中间图像的原子间距在初末态之间平滑过渡，比线性插值更少出现原子“穿过”彼此的初始路径。

### 2.5 形成能与凸包能

对化合物 $\mathrm{A}_x\mathrm{B}_y$，每原子形成能：

$$\Delta E_f = \frac{E(\mathrm{A}_x\mathrm{B}_y) - x\,\mu_\mathrm{A} - y\,\mu_\mathrm{B}}{x+y}$$

所有相的 $\Delta E_f$ 在组成空间的下凸包构成相图。某相的凸包能 $E_\text{hull} = \Delta E_f - \Delta E_f^\text{hull}$（同组成处凸包上的值），稳定相为 0，亚稳相为正值。只有**计算设置一致**（同一泛函、U 值、POTCAR、能量修正方案）的能量才能放进同一张相图。

---

## 3. 安装与版本要点

```bash
conda activate dft                       # 第1章建立的 Python 3.12 环境
pip install "ase>=3.28" "pymatgen>=2026.3" mp-api "cclib==1.8.1" dpdata pandas tables
python -c "import ase; print('ase', ase.__version__)"
python -c "import pymatgen.core as pc; print('pymatgen-core', pc.__version__)"
python -c "import cclib, dpdata; print('cclib', cclib.__version__, 'dpdata', dpdata.__version__)"
ase info --calculators                    # 列出 ASE 能找到的计算器及其配置
```

| 版本变化 | 影响 | 处理 |
|---|---|---|
| ASE 3.28.0（2026-03）起 NEB 默认 `method` 由 `aseneb` 改为 `improvedtangent` | 不写 `method` 会给出 UserWarning，旧脚本结果可能变化 | 始终显式写 `NEB(images, method="improvedtangent")`；复现旧结果写 `method="aseneb"` |
| ASE 3.26 起有实验性的 `ase._4` 命名空间 | 下划线开头的包随时可能变 | 生产脚本不要 import `ase._4` |
| ASE 文档 2026-06 迁到 <https://docs.ase-lib.org/> | 旧链接（ase-lib.org、wiki.fysik.dtu.dk）可能失效 | 收藏新地址；源码在 GitLab（gitlab.com/ase/ase） |
| pymatgen v2026.3.23 起核心代码拆为 `pymatgen-core` 单独发布 | `pip install pymatgen` 不受影响，导入路径不变 | 统一写 `from pymatgen.core import Structure` |
| `from pymatgen import Structure` 早已移除 | `ImportError: cannot import name 'Structure' from 'pymatgen'` | 改为 `from pymatgen.core import Structure` |
| pymatgen-core v2026.7.31 新增 `MatPESRelaxSet`，并从 `MP24RelaxSet`/`MPScanRelaxSet` 等去掉 `ENAUG` | 新旧版本生成的 INCAR 不同 | 论文中写明 pymatgen 版本；需要时用 `user_incar_settings` 显式设定 |
| `Oszicar` 遇到 VASP 的 `******` 溢出值返回 `nan` 并警告，不再抛异常 | 能量可能是 nan | 读完后检查 `math.isnan` |
| `BaderAnalysis` 的 `parse_atomic_densities` 参数已弃用，计划 2026-12-31 后删除 | 旧脚本将来报错 | 不要再传该参数 |
| Materials Project 旧版 REST API 于 2025-09 下线 | 旧 key、旧 URL 报 403 或连接错误 | 到 dashboard 取 32 位新 key，`export MP_API_KEY=...` |

> [!IMPORTANT]
> cclib 2.0a1 是 alpha 版，解析 Gaussian 16 输出有已知 bug（`ValueError: invalid literal for int() with base 10: 'electrons'`），生产环境固定 `cclib==1.8.1`。

本节用到的辅助脚本：[env_check.sh](../scripts/env_check.sh) 检查环境变量和 Python 包是否齐全，用 `bash ../scripts/env_check.sh -h` 查看用法。

---

## 4. ASE 基础：Atoms、约束与文件读写

### 4.1 建一个带约束的表面模型

```python
from ase.build import fcc111, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.io import read, write

slab = fcc111("Pt", size=(3, 3, 4), vacuum=7.5)      # 4 layers, 15 A vacuum in total
add_adsorbate(slab, molecule("CO"), height=1.9, position="fcc", mol_index=1)  # atom 1 of CO is C
# tag: 1 = top layer ... 4 = bottom layer, adsorbate = 0
slab.set_constraint(FixAtoms(mask=[a.tag >= 3 for a in slab]))   # fix bottom 2 layers
print(slab)                                          # Atoms(symbols='Pt36OC', pbc=[True, True, False], ...)
print(slab.get_chemical_formula(), len(slab), slab.cell.cellpar().round(3))
write("POSCAR", slab, format="vasp", direct=True, sort=True)     # Selective dynamics F F F
write("slab.cif", slab)
write("slab.traj", slab)                             # keeps constraints and tags
```

| 对象/方法 | 说明 |
|---|---|
| `atoms.positions`、`atoms.cell`、`atoms.pbc` | 笛卡尔坐标（Å）、晶胞矩阵（Å）、各方向是否周期 |
| `atoms.get_scaled_positions()` | 分数坐标 |
| `atoms.set_cell(cell, scale_atoms=True)` | 改晶胞同时按比例移动原子（EOS 扫描用） |
| `atoms.center(vacuum=5.0)` | 分子居中并在每侧留 5 Å 真空 |
| `atoms.repeat((2, 2, 1))` / `ase.build.make_supercell` | 超胞 |
| `FixAtoms(indices=[...])` 或 `FixAtoms(mask=[...])` | 固定原子；写 POSCAR 时自动加 `Selective dynamics` |
| `atoms.info`、`atoms.arrays` | 自定义元数据与逐原子数组（extxyz 会一并保存） |

### 4.2 常用格式

`ase.io.read(file, index=..., format=...)` 中 `index=-1` 取最后一帧，`index=":"` 取全部帧。

| format | 读/写 | 典型文件 | 备注 |
|---|---|---|---|
| `vasp` | 读写 | POSCAR、CONTCAR | `direct=True` 写分数坐标；`sort=True` 按元素排序 |
| `vasp-out` | 读 | OUTCAR | 每个离子步的能量、力、应力；`energy` 是 `energy(sigma->0)`，`free_energy` 是 `TOTEN` |
| `vasp-xml` | 读 | vasprun.xml | 同上，另含 k 点和本征值 |
| `vasp-xdatcar` | 读写 | XDATCAR | AIMD 轨迹 |
| `gaussian-in` / `gaussian-out` | 写 / 读 | .com/.gjf、.log | 读 log 得到每步几何、能量、力（需要 `force` 关键词） |
| `cp2k-restart`、`cp2k-dcd` | 读 | `*-1.restart`、`*.dcd` | 主输出中的能量用正则或 pymatgen 读取 |
| `xyz` / `extxyz` | 读写 | `*-pos-1.xyz`、训练集 | extxyz 保存晶胞、能量、力，是 MLIP 的通用格式 |
| `cif` | 读写 | 数据库下载的结构 | 部分占位（无序）CIF 需先在 pymatgen 中处理 |
| `traj` / `db` | 读写 | ASE 轨迹、ase.db | 保存约束、计算结果，`ase gui` 可直接打开 |

```bash
ase info --formats | head -40              # list all formats known to ASE
ase convert CONTCAR final.cif              # one-line format conversion
ase gui OUTCAR@:                           # view all ionic steps of a VASP run
ase gui g16/h2o_opt.log                    # Gaussian optimisation trajectory
```

OpenMX 的输出没有注册为 `ase.io` 格式，用计算器模块里的读取函数：`from ase.calculators.openmx.reader import read_openmx; atoms = read_openmx("omx/si")`（自动读 `si.out`、`si.dat`、`si.log` 等）。结构建模的更多技巧（数据库下载、表面、缺陷、异质结）见 [第3章](03-structure-modeling.md)；批量格式互转可用辅助脚本 [convert_any.py](../scripts/convert_any.py)，用 `python ../scripts/convert_any.py -h` 查看用法。

---

## 5. 先用 EMT 跑通流程：EOS、约束优化、NEB 与 ase db

示例目录 [04-ase-emt-basics](../examples/04-ase-emt-basics/README.md) 使用 ASE 自带的 EMT 经验势，在笔记本上几秒到一分钟就能跑完。流程跑通后，把 `EMT()` 换成第 6 节的 DFT 计算器即可，其余代码不变。

```bash
cd Common/examples/04-ase-emt-basics
export MPLBACKEND=Agg                       # no display on clusters
python ase_eos_cu.py
python ase_relax_constraints.py
python ase_neb_idpp.py --nimages 5 --method improvedtangent
python ase_db_demo.py
```

> [!WARNING]
> EMT 只适用于 Al、Cu、Ag、Au、Ni、Pd、Pt 等少数金属，结果**不是** DFT 结果，只用来练习流程，不能与实验比较。

### 5.1 EOS：晶格常数与体模量

[ase_eos_cu.py](../examples/04-ase-emt-basics/ase_eos_cu.py) 把 Cu 原胞缩放 0.97–1.03 倍取 7 个点，核心几行是：

```python
from ase.eos import EquationOfState
from ase.units import kJ
eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
v0, e0, B = eos.fit()                # B in eV/A^3
print(B / kJ * 1.0e24, "GPa")
eos.plot("cu_eos.png")
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| 缩放范围 | 晶格常数 ±3%（体积约 ±9%） | 太窄拟合噪声大，太宽偏离谐和区 |
| 点数 | 7–11 | 至少 5 点；最低点应在区间内部 |
| `eos` | `birchmurnaghan`（固体常用），也可 `vinet`、`murnaghan`、`sj` | 不同形式的 $V_0$ 差别应很小，可用来检查拟合是否可靠 |
| DFT 计算器 | 每个体积一个 `directory` | `ase_eos_cu.py --calc vasp` 中每点写到 `eos_00/` … `eos_06/`，互不覆盖 |

成功标准：`cu_eos.png` 中曲线穿过所有点、最低点在内部；用 VASP 计算时还要做 ENCUT 收敛（体积变化会引入 Pulay 应力，见 [VASP 第4章](../../VASP/tutorials/04-relax-eos.md)）。

### 5.2 约束优化：BFGS 与 FIRE

[ase_relax_constraints.py](../examples/04-ase-emt-basics/ase_relax_constraints.py) 在 Cu(100) 上放一个吸附原子，固定下两层，分别用 BFGS 和 FIRE 优化：

```python
opt = BFGS(atoms, trajectory="relax_bfgs.traj", logfile="relax_bfgs.log")
converged = opt.run(fmax=0.02, steps=500)     # True if fmax was reached within steps
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `fmax` | 普通结构 0.02–0.05 eV/Å；声子、弹性前 ≤ 0.01 eV/Å | 与 VASP `EDIFFG = -0.02` 等价的判据 |
| `steps` | 200–500 | 防止不收敛时无限运行；`run()` 返回值判断是否真的收敛 |
| `trajectory` | 必写 | 中断后 `read("relax.traj", -1)` 取最后结构续算 |
| 优化器 | 先 BFGS；初始结构很差或出现振荡时用 FIRE | `BFGSLineSearch`、`LBFGS` 适合大体系 |
| 晶胞优化 | `FrechetCellFilter(atoms)` 包一层再交给优化器 | ASE ≥ 3.23 推荐，替代旧的 `ExpCellFilter`/`UnitCellFilter` |

检查：`relax_bfgs.log` 最后一行的 `fmax` 小于阈值；`CONTCAR_bfgs` 中固定原子带 `F F F`。

### 5.3 NEB：IDPP 插值 + 爬坡图像

[ase_neb_idpp.py](../examples/04-ase-emt-basics/ase_neb_idpp.py) 计算吸附原子在相邻 hollow 位之间的扩散势垒，关键步骤：

```python
from ase.mep import NEB, NEBTools
images = [initial] + [initial.copy() for _ in range(5)] + [final]
for img in images[1:-1]:
    img.calc = EMT()                         # every image needs its own calculator
neb = NEB(images, k=0.1, climb=False, method="improvedtangent")
neb.interpolate(method="idpp")               # IDPP instead of linear
FIRE(neb, trajectory="neb.traj").run(fmax=0.10)    # stage 1: plain NEB
neb.climb = True
FIRE(neb, trajectory="neb.traj", append_trajectory=True).run(fmax=0.03)  # stage 2: CI-NEB
barrier, de = NEBTools(images).get_barrier()
```

| 参数 | 推荐值 | 说明 |
|---|---|---|
| `method` | `"improvedtangent"`（显式写出） | 3.28 起的默认值；`"aseneb"` 仅用于复现旧结果 |
| `k` | 0.1 eV/Å² 左右 | 弹簧常数；图像分布不均时调大 |
| 图像数 | 3–7 | 路径长或有中间体时增加 |
| `interpolate(method="idpp")` | 推荐 | 线性插值常导致原子过近、首步力极大 |
| 两段优化 | 先不爬坡到 0.1 eV/Å，再开 `climb=True` 到 0.03–0.05 eV/Å | 直接爬坡容易让 CI 爬到错误的鞍点 |
| `NEB(..., parallel=True)` | 每个图像一组 MPI 进程 | DFT 计算时可并行；或直接用 VTST/VASP 内置 NEB |

检查：`ase nebplot neb.traj` 或 `neb_barrier.png` 中最高点在中间图像；初末态等价时反应能应为 0。DFT 级别的 NEB 一般直接用程序自身的实现（[VASP 第11章](../../VASP/tutorials/11-neb-dimer.md)、[CP2K 第13章](../../CP2K/tutorials/13-neb-dimer.md)、[OpenMX 第6章](../../OpenMX/tutorials/06-optimization-eos-neb.md)），ASE 的 NEB 适合 MLIP 预探索路径或生成 IDPP 初始图像。

### 5.4 ase db 与 ase gui

[ase_db_demo.py](../examples/04-ase-emt-basics/ase_db_demo.py) 把 6 种金属的晶格常数和体模量写进 SQLite 数据库：

```python
from ase.db import connect
db = connect("metals.db")
db.write(atoms, metal="Cu", a0=3.59, B_GPa=134.3, method="EMT",
         data={"volumes": vols, "energies": ens})   # data: arrays, not searchable
for row in db.select("a0>3.9", sort="a0"):
    print(row.formula, row.a0, row.energy)
atoms = db.get(metal="Cu").toatoms()
```

```bash
ase db metals.db                                     # list all rows
ase db metals.db 'a0>3.9' -c id,metal,a0,B_GPa       # filter + choose columns
ase db metals.db metal=Cu -L 0 --csv > cu.csv        # export as CSV
ase gui metals.db@metal=Cu                           # open selected rows in the GUI
```

| 要点 | 说明 |
|---|---|
| 键值对（key_value_pairs） | 可查询的标量；`energy`、`forces`、`magmom`、`fmax`、`formula`、`natoms`、`id`、`calculator` 以及元素符号（如 `Cu`）等是保留名（见 `ase.db.core.reserved_keys`），用作键会报 `ValueError: Bad key: ...` |
| `data=` | 任意数组/字典，不可查询，用于保存 EOS 点、DOS 等 |
| 后端 | `.db` 为 SQLite（单用户）；多人同时写入时用 PostgreSQL |
| 幂等 | 写入前先 `db.count(...)` 判断是否已存在，脚本可以重复运行 |

示例 README 中记录的 EMT 结果供核对流程（EMT 势，不是 DFT 结果）：Cu 的 $a_0 \approx 3.59$ Å、$B \approx 134$ GPa；BFGS 与 FIRE 最终能量相差 < 1 meV，BFGS 步数更少；NEB 的反应能为 0。

---

## 6. 用 ASE 驱动四款软件

示例目录：[04-ase-calculators](../examples/04-ase-calculators/README.md)。四个脚本都先支持 `--dry-run`（只写输入、不运行程序），**先检查生成的输入，再提交计算**。

### 6.1 告诉 ASE 怎么运行程序

ASE 按“计算器参数 `command` → 配置文件 → 环境变量”的顺序寻找运行命令。最简单的是在 `~/.bashrc` 或作业脚本里设置环境变量：

```bash
# VASP: VASP_PP_PATH points to the directory that CONTAINS potpaw_PBE/
export VASP_PP_PATH=$HOME/pp                       # $HOME/pp/potpaw_PBE/Si/POTCAR
export ASE_VASP_COMMAND="mpirun -np 16 vasp_std"
# Gaussian: PREFIX is replaced by the label
export ASE_GAUSSIAN_COMMAND="g16 < PREFIX.com > PREFIX.log"
export GAUSS_SCRDIR=/tmp/$USER
# CP2K: shell mode (-s) keeps one cp2k process alive for many calls
export CP2K_DATA_DIR=/path/to/cp2k/data
export ASE_CP2K_COMMAND="mpirun -np 4 cp2k.psmp -s"
# OpenMX
export OPENMX_DFT_DATA_PATH=/path/to/openmx4.0/DFT_DATA19
export ASE_OPENMX_COMMAND=openmx
```

也可以写在 `~/.config/ase/config.ini` 的 `[environment]` 段（路径可由 `ASE_CONFIG_PATH` 改变）：

```ini
[environment]
VASP_PP_PATH = /home/user/pp
ASE_VASP_COMMAND = mpirun -np 16 vasp_std
ASE_GAUSSIAN_COMMAND = g16 < PREFIX.com > PREFIX.log
```

> [!NOTE]
> config.ini 使用 Python `configparser` 的扩展插值（`ExtendedInterpolation`），值里不要出现 `$` 变量；需要 `$SLURM_NTASKS` 这类变量时，在作业脚本里 `export` 环境变量。ASE 3.29 仍然读取环境变量，但会发出 `ASEEnvDeprecationWarning`（`Loaded ... from environment. Please use configfile.`，属于 DeprecationWarning，默认不显示），官方推荐的方向是配置文件。`ase info --calculators` 可查看哪些计算器已配置。

### 6.2 VASP：两种优化方式

[ase_vasp_relax.py](../examples/04-ase-calculators/ase_vasp_relax.py) 以体相 Si 为例：

```python
from ase.build import bulk
from ase.calculators.vasp import Vasp
si = bulk("Si", "diamond", a=5.40)
si.calc = Vasp(directory="si_relax_vasp", xc="PBE", setups="recommended",
               encut=520, prec="Accurate", kpts=(8, 8, 8), gamma=True,
               ismear=0, sigma=0.05, ediff=1e-7, lreal=False,
               ibrion=2, isif=3, nsw=60, ediffg=-0.005,
               lwave=False, lcharg=False, ncore=4)
e = si.get_potential_energy()       # runs VASP once; atoms updated from CONTCAR
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `xc` | `"PBE"` | 同时决定 POTCAR 目录（`potpaw_PBE`）和 `GGA` 标签 |
| `setups` | `"recommended"` | 采用 VASP 推荐的 POTCAR 变体（如 `Ti_sv`、`Fe_pv`）；也可传字典逐元素指定 |
| `encut` | ≥ 1.3×ENMAX；变胞优化用更高值 | 变胞时 Pulay 应力需要更高截断能 |
| `kpts`、`gamma` | `(8,8,8)`、`True` | Γ 中心网格；先做收敛测试（[VASP 第3章](../../VASP/tutorials/03-convergence.md)） |
| `ibrion/isif/nsw/ediffg` | `2/3/60/-0.005` | `--mode vasp`：VASP 自己优化，ASE 只运行一次 |
| `ibrion=-1, nsw=0` | — | `--mode ase`：VASP 只算单点，ASE 的 `BFGS + FrechetCellFilter` 负责优化 |
| `directory` | 每个任务一个目录 | 避免输出互相覆盖 |

```bash
python ase_vasp_relax.py --dry-run           # writes si_relax_vasp/{INCAR,KPOINTS,POSCAR,POTCAR}
head -1 si_relax_vasp/POTCAR                 # check which POTCAR was picked
sbatch job_ase_vasp.slurm                    # runs --mode vasp inside Slurm
grep "reached required accuracy" si_relax_vasp/OUTCAR
grep "free  energy   TOTEN" si_relax_vasp/OUTCAR | tail -1
```

成功标准：OUTCAR 中出现 `reached required accuracy`；PBE 下 Si 的立方晶格常数约 5.47 Å（文献 PBE 值 5.47 Å，实验 5.431 Å，见 Haas *et al.*, DOI: 10.1103/PhysRevB.79.085104）。

> [!TIP]
> 两种模式怎么选：结构优化、MD 优先 `--mode vasp`（VASP 内部的优化器能重用波函数，速度快）；需要 ASE 特有功能（约束、自定义势能面、与 MLIP 混合）时才用 `--mode ase`。作业脚本模板见 [job_ase_vasp.slurm](../examples/04-ase-calculators/job_ase_vasp.slurm)：Python 进程本身是串行的，ASE 每次调用 VASP 时才用 `ASE_VASP_COMMAND` 中的 `srun`/`mpirun` 启动并行程序。

### 6.3 Gaussian：用 Gaussian 自己的优化器

[ase_gaussian_h2o_opt.py](../examples/04-ase-calculators/ase_gaussian_h2o_opt.py) 先优化 H₂O，再在同一理论级别做频率：

```python
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
calc = Gaussian(label="g16/h2o_opt", chk="h2o.chk", mem="4GB", nprocshared=8,
                method="B3LYP", basis="def2SVP", EmpiricalDispersion="GD3BJ",
                scf="tight", charge=0, mult=1)
GaussianOptimizer(h2o, calc).run(fmax="tight", steps=100)   # -> opt=(tight,maxcycle=100)
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `label` | `"g16/h2o_opt"` | 生成 `g16/h2o_opt.com` 与 `.log` |
| `mem`、`nprocshared`、`chk` | 与作业申请资源一致 | 写成 `%mem`、`%nprocshared`、`%chk` |
| `method`、`basis` | 如 `B3LYP`、`def2SVP` | 路由行 `# B3LYP/def2SVP` |
| 其它关键词 | `EmpiricalDispersion="GD3BJ"`、`scf="tight"` | 任意 `keyword=value` 都会写入路由行 |
| `GaussianOptimizer` | `fmax="tight"`、`steps=100` | 调用 Gaussian 内置优化（冗余内坐标），远比 ASE 优化器逐步调用 Gaussian 高效 |

```bash
python ase_gaussian_h2o_opt.py --dry-run     # print g16/h2o_opt.com only
python ase_gaussian_h2o_opt.py
tail -1 g16/h2o_opt.log                      # Normal termination of Gaussian 16 ...
grep -c "imaginary frequencies" g16/h2o_freq.log   # 0 for a minimum
```

ASE 不解析频率；频率、轨道能、热化学用 cclib 或 pymatgen 读取（第 10、11 节）。Gaussian 输入的各部分含义见 [Gaussian 第2章](../../Gaussian/tutorials/02-input-anatomy.md)。

### 6.4 CP2K：shell 模式与 `inp` 模板

[ase_cp2k_water_sp.py](../examples/04-ase-calculators/ase_cp2k_water_sp.py) 计算孤立水分子的能量和力：

```python
from ase.calculators.cp2k import CP2K
calc = CP2K(label="h2o", xc="PBE",
            basis_set="DZVP-MOLOPT-SR-GTH", basis_set_file="BASIS_MOLOPT",
            pseudo_potential="GTH-PBE", potential_file="GTH_POTENTIALS",
            cutoff=400 * 13.605693,          # ASE expects eV: 400 Ry
            max_scf=30, stress_tensor=False, inp=INP)
with calc:                                   # closes the cp2k shell at the end
    h2o.calc = calc
    print(h2o.get_potential_energy(), h2o.get_forces())
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `cutoff` | 400 Ry × 13.605693（**以 eV 输入**） | 对应 `&MGRID CUTOFF`；先做收敛测试（[CP2K 第4章](../../CP2K/tutorials/04-cutoff-convergence.md)） |
| `basis_set` / `pseudo_potential` | `DZVP-MOLOPT-SR-GTH` / `GTH-PBE` | 文件名用 `basis_set_file`、`potential_file` 指定，路径由 `CP2K_DATA_DIR` 提供 |
| `stress_tensor` | 孤立分子设 `False` | 默认 `True`，会计算应力 |
| `inp` | CP2K 输入片段字符串 | 简单参数覆盖不到的设置（OT、OUTER_SCF、EPS_DEFAULT、REL_CUTOFF 等）合并进生成的输入 |
| `pbc=False` | 孤立分子 | ASE 自动加 `POISSON PERIODIC NONE` 与 `PSOLVER MT` |
| `debug=True` | 排错时 | 打印 ASE 与 cp2k shell 之间的通信 |

```bash
python ase_cp2k_water_sp.py                  # single point
python ase_cp2k_water_sp.py --relax          # + ASE BFGS relaxation
grep "ENERGY| Total FORCE_EVAL" h2o.out | tail -1
grep -c "SCF run NOT converged" h2o.out      # must be 0
```

> [!WARNING]
> ASE 生成的输入只含 `&GLOBAL`、`&FORCE_EVAL` 等基本部分，几何优化、MD 由 ASE 负责，因此不涉及 CP2K 2026.2 删除的 `&MOTION/&CELL_OPT/TYPE`。但如果你在 `inp` 模板里加了 `&DFT/&PRINT/&PDOS` 这类 2026.2 已改结构的段落，需按新语法 `&DFT/&PRINT/&DOS`（内嵌 `&PDOS`）书写，见 [CP2K 第10章](../../CP2K/tutorials/10-dos-pdos-band.md)。

### 6.5 OpenMX：注意单位换算

[ase_openmx_si.py](../examples/04-ase-calculators/ase_openmx_si.py)：

```python
from ase.calculators.openmx import OpenMX
from ase.units import Ha, Ry
calc = OpenMX(label="omx/si", xc="PBE",
              energy_cutoff=220 * Ry,          # eV in, written as scf.energycutoff 220 (Ry)
              convergence=1.0e-7 * Ha,         # eV in, written as scf.criterion 1e-7 (Hartree)
              kpts=(8, 8, 8), maxiter=100, eigensolver="Band",
              dft_data_year="19",
              definition_of_atomic_species=[["Si", "Si7.0-s2p2d1", "Si_PBE19"]],
              mpi={"processes": 8, "threads": 1},
              system_currentdirectory="./")
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `energy_cutoff` | `220 * Ry` | **以 eV 输入**，写出时换算回 Ry |
| `convergence` | `1e-7 * Ha` | **以 eV 输入**，写出为 Hartree 的 `scf.criterion` |
| `dft_data_year` | `"19"` | 使用 DFT_DATA19（OpenMX 3.9/4.0） |
| `definition_of_atomic_species` | `[元素, PAO, VPS]` | PAO/VPS 选择见 [OpenMX 第4章](../../OpenMX/tutorials/04-pao-vps.md) |
| `mpi` | `{"processes": N, "threads": 1}` | ASE 拼成 `mpirun -np N openmx si.dat -nt 1` |
| `system_currentdirectory` | `"./"` | ASE 默认写入绝对路径，显式设为 `./` 让输出与 `.dat` 在同一目录 |

```bash
python ase_openmx_si.py --dry-run            # print omx/si.dat
python ase_openmx_si.py --np 8
grep "Utot\." omx/si.out
```

ASE 的 OpenMX 接口只覆盖常见关键词；DFT+U、非共线、NEGF、Wannier 等高级功能请直接手写 `.dat`（见 OpenMX 各章），ASE 只用来读结果。

### 6.6 MLIP 作为计算器（前沿）

MACE、DeePMD、NEP 等机器学习势都提供 ASE 计算器。常见用法是：用通用势（如 MACE-MP）对初始结构做预优化、为 NEB 生成合理的初始路径，然后再交给 DFT 精修。代码层面只需把 `atoms.calc = EMT()` 换成对应的计算器，参考 [ACEsuit/mace](https://github.com/ACEsuit/mace) 的文档。预优化结果不能直接当作 DFT 结果报告。

---

## 7. pymatgen：结构与 VASP 输入集

### 7.1 Structure 与对称性

```python
from pymatgen.core import Lattice, Structure        # NOT: from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

s = Structure.from_file("POSCAR")                     # POSCAR / CIF / ... detected by name
sga = SpacegroupAnalyzer(s, symprec=0.01)
print(sga.get_space_group_symbol(), sga.get_space_group_number())
prim = sga.get_primitive_standard_structure()
prim.to(filename="POSCAR_prim")                       # format chosen by file name
atoms = AseAtomsAdaptor.get_atoms(prim)               # pymatgen -> ASE
s2 = AseAtomsAdaptor.get_structure(atoms)             # ASE -> pymatgen
```

`symprec` 默认 0.01 Å；结构来自实验或弛豫后略有畸变时可放宽到 0.05–0.1。原胞/标准胞、超胞、表面的构建见 [第3章](03-structure-modeling.md)，能带路径可用 [kpath_gen.py](../scripts/kpath_gen.py) 生成（`python ../scripts/kpath_gen.py -h` 查看用法）。

### 7.2 输入集一览

| 输入集 | 用途 | 主要默认值（以 pymatgen 2026.10 的 YAML 为准） |
|---|---|---|
| `MPRelaxSet` | MP 标准 PBE(+U) 结构优化 | `ENCUT=520`、`ISIF=3`、`ISPIN=2`、`ISMEAR=-5`、`EDIFF_PER_ATOM=5e-5`、`LREAL=Auto`、`reciprocal_density=64`、POTCAR 为 `PBE`，含过渡金属氧化物/氟化物的 U 值 |
| `MPStaticSet` | 静态 SCF（写 CHGCAR，给后续非自洽计算用） | 继承 MPRelaxSet，`NSW=0`、`LCHARG=True`、`LORBIT=11` |
| `MPNonSCFSet` | 非自洽能带/DOS | `mode="line"`（高对称路径）、`"uniform"`（DOS 密网格）、`"boltztrap"`；`ICHARG=11` |
| `MPHSEBSSet` | HSE06 能带（自洽，零权重 k 点） | `mode` 为 `"gap"`（默认）、`"line"`、`"uniform"`、`"uniform_dense"` |
| `MPScanRelaxSet` | r2SCAN 结构优化 | `METAGGA=R2SCAN`；`auto_kspacing=True`，从上一步计算（`from_prev_calc` 或传入 `bandgap`）的带隙确定 `KSPACING` |
| `MP24RelaxSet` / `MP24StaticSet` | MP 2024 新标准（r2SCAN） | `ENCUT=680`、`LMAXMIX=6`、POTCAR 为 `PBE_64`；没有带隙信息时 `KSPACING=0.22`，有上一步计算时按带隙自动设定 |
| `MatPESStaticSet` / `MatPESRelaxSet` | MatPES 数据集标准（默认 PBE），适合 MLIP 训练数据 | `ENCUT=680`、`KSPACING=0.22`、`PBE_64`；`MatPESRelaxSet` 需 pymatgen-core ≥ 2026.7.31 |

通用参数：

| 参数 | 说明 |
|---|---|
| `user_incar_settings={...}` | 覆盖或新增 INCAR 标签，优先级最高；值为 `None` 表示删除该标签，如 `{"LREAL": None}` |
| `user_kpoints_settings={"reciprocal_density": 200}` | 改 k 点密度；也可传 `Kpoints` 对象 |
| `user_potcar_settings={"Fe": "Fe"}` | 改某元素的 POTCAR 变体 |
| `user_potcar_functional` | 如 `"PBE_54"`、`"PBE_64"` |
| `XxxSet.from_prev_calc(dir, ...)` | 从上一步目录继承结构、`ISPIN`、`MAGMOM`、`NBANDS` 等，并可复制 CHGCAR |
| `write_input(dir, potcar_spec=True)` | 只写 `POTCAR.spec`（元素符号列表），不需要 POTCAR 库，适合在笔记本上检查 |

### 7.3 一套完整的 Si 流程

[pmg_write_inputset.py](../examples/04-pymatgen-io/pmg_write_inputset.py) 生成 `01_relax` … `07_matpes_static` 七个目录：

```python
relax = MPRelaxSet(s, user_incar_settings={"ENCUT": 520, "EDIFFG": -0.01,
                                            "ISMEAR": 0, "SIGMA": 0.05,
                                            "NCORE": 4, "LWAVE": False},
                   user_kpoints_settings={"reciprocal_density": 200})
relax.write_input("01_relax", potcar_spec=True)
MPStaticSet(s, user_incar_settings={"ENCUT": 520, "NCORE": 4}).write_input("02_static", potcar_spec=True)
MPNonSCFSet(s, mode="line", kpoints_line_density=40).write_input("03_band", potcar_spec=True)
MPNonSCFSet(s, mode="uniform", reciprocal_density=400, nedos=3001).write_input("04_dos", potcar_spec=True)
MPHSEBSSet(s, mode="line").write_input("05_hse_band", potcar_spec=True)
```

真实流程中第 2 步之后应写成 `MPStaticSet.from_prev_calc("01_relax")`、`MPNonSCFSet.from_prev_calc("02_static", mode="line")`，让结构、自旋、能带数从上一步继承。

```bash
cd Common/examples/04-pymatgen-io
python pmg_write_inputset.py
cat 01_relax/INCAR 03_band/INCAR | grep -E "ENCUT|ISPIN|ICHARG|ISMEAR"
cat 01_relax/POTCAR.spec
```

检查要点：`03_band/INCAR` 中 `ICHARG = 11`；`05_hse_band` 中 `LHFCALC = True`；`06_mp24_relax` 中 `METAGGA = R2scan`、`KSPACING = 0.22`。

> [!WARNING]
> `MPRelaxSet` 默认 `ISPIN = 2`，并按元素给 `MAGMOM`（配置中未列出的元素用 0.6 μB）。Si 这类非磁体系可加 `"ISPIN": 1` 省一半时间；磁性体系请按 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md) 自己设定初始磁矩。

### 7.4 配置 POTCAR 库

要让 `write_input` 写出真正的 POTCAR，先用 `pmg config` 把 VASP 发布的赝势目录整理成 pymatgen 的格式：

```bash
pmg config -p /path/to/vasp_pp ~/pmg_psp       # vasp_pp contains potpaw_PBE/ potpaw_PBE.54/ ...
pmg config --add PMG_VASP_PSP_DIR ~/pmg_psp
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54   # only if your default should be PBE_54
python pmg_write_inputset.py --potcar
```

配置写在 `~/.config/.pmgrc.yaml`（旧版本为 `~/.pmgrc.yaml`；可用环境变量 `PMG_CONFIG_FILE` 指定其他位置），以 `PMG_` 开头的同名环境变量会覆盖文件中的设置。

> [!CAUTION]
> POTCAR 受 VASP 许可证保护，不要提交到 Git 或公开上传。同理，`vaspout.h5`（VASP 6.5 起加强的 HDF5 输出）中也包含 POTCAR 数据，分享前要确认已去除。

---

## 8. pymatgen 解析 VASP 输出

### 8.1 Vasprun、Outcar、Oszicar

[pmg_parse_vasprun.py](../examples/04-pymatgen-io/pmg_parse_vasprun.py) 读一个计算目录并汇总：

```python
from pymatgen.io.vasp.outputs import Vasprun, Outcar, Oszicar
vr = Vasprun("vasprun.xml", parse_potcar_file=False, exception_on_bad_xml=False)
print(vr.final_energy, vr.converged_electronic, vr.converged_ionic, vr.efermi)
gap, cbm, vbm, is_direct = vr.eigenvalue_band_properties
oc = Outcar("OUTCAR");   print(oc.total_mag)          # total magnetisation (muB)
osz = Oszicar("OSZICAR"); print(osz.final_energy, len(osz.ionic_steps))
```

| 参数/属性 | 说明 |
|---|---|
| `exception_on_bad_xml=False` | 作业中断、`vasprun.xml` 不完整时，保留已完整写出的离子步，只给警告；默认 `True` 会报 XML ParseError |
| `parse_potcar_file=False` | 不去找 POTCAR 校验（换了机器时避免报错） |
| `parse_projected_eigen=True` | 读取投影本征值（画 fat band 需要），占内存大，默认 `False` |
| `vr.final_energy` | `energy(sigma->0)`，单位 eV |
| `vr.converged_electronic / converged_ionic` | 最后一步 SCF 是否收敛 / 离子步是否在 NSW 内收敛 |
| `vr.eigenvalue_band_properties` | 按占据数估计带隙；金属可能给出虚假的小带隙，判断金属性用能带对象 |
| `Oszicar.final_energy` | 遇到 `******` 溢出值时为 `nan`（2026 新版行为） |

```bash
python pmg_parse_vasprun.py path/to/run              # summary
python pmg_parse_vasprun.py path/to/band_run --band  # band.png (line-mode KPOINTS required)
python pmg_parse_vasprun.py path/to/dos_run --dos    # dos.png (LORBIT = 11)
python pmg_parse_vasprun.py path/to/slab --locpot    # planar average + work function
```

### 8.2 能带、DOS 与带隙

```python
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
bs = vr.get_band_structure(kpoints_filename="KPOINTS", line_mode=True, efermi="smart")
print(bs.is_metal(), bs.get_band_gap())          # {'energy':..., 'direct':..., 'transition': 'X-\\Gamma'}
BSPlotter(bs).get_plot(ylim=(-6, 6)).figure.savefig("band.png", dpi=300)
dp = DosPlotter(sigma=0.05)
dp.add_dos("Total", vr.complete_dos)
dp.add_dos_dict(vr.complete_dos.get_element_dos())
dp.get_plot(xlim=(-8, 6)).figure.savefig("dos.png", dpi=300)
```

- 能带必须读**非自洽 line 模式**的 `vasprun.xml`，`KPOINTS` 要带高对称点标签；`efermi="smart"` 会把费米能放到带隙中合理位置，避免半导体被误判为金属。
- 用 `KPOINTS_OPT`（VASP ≥ 6.3）计算的能带，目录中要有 `KPOINTS_OPT` 文件，细节见 [VASP 第5章](../../VASP/tutorials/05-band-dos.md)。
- 新版 pymatgen 的 `get_plot()` 返回 matplotlib `Axes`，用 `.figure.savefig()` 保存。出版级作图和 sumo/pyprocar 的用法见 [第5章](05-postprocess-visualization.md)。

### 8.3 Chgcar 与 Locpot：体数据运算

`Chgcar`、`Locpot` 都继承自 `VolumetricData`，支持 `+`、`-` 运算和沿某一轴的平面平均。[pmg_chgcar_ops.py](../examples/04-pymatgen-io/pmg_chgcar_ops.py) 实现了两个常见需求：

```bash
# 1) Bader reference density: AECCAR0 + AECCAR2 (needs LAECHG = .TRUE.)
python pmg_chgcar_ops.py sum AECCAR0 AECCAR2          # -> CHGCAR_sum
bader CHGCAR -ref CHGCAR_sum
# 2) charge density difference rho(AB) - rho(A) - rho(B)
python pmg_chgcar_ops.py diff AB/CHGCAR A/CHGCAR B/CHGCAR   # -> CHGCAR_diff, diff_z.dat
```

CHGCAR 中存的是 $\rho \cdot V_\text{cell}$，脚本输出的平面积分差分电荷为

$$\Delta\lambda(z) = \int\!\!\int \Delta\rho(x,y,z)\,\mathrm{d}x\,\mathrm{d}y = \frac{\langle \Delta\rho V\rangle_{xy}}{h},\qquad h = V/A_{xy}$$

| 检查项 | 成功标准 |
|---|---|
| 网格 | 三次计算的 FFT 网格一致；否则脚本报 `FFT grids differ`，需显式设置相同的 `NGXF/NGYF/NGZF`、ENCUT 和晶胞 |
| 电荷守恒 | `net charge of d_rho` 接近 0 |
| Bader | `ACF.dat` 中 CHARGE 是价电子数，净电荷 = ZVAL − CHARGE（详见 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)） |

功函数：`Locpot.from_file("LOCPOT").get_average_along_axis(2)` 给出沿 z 的平面平均静电势，真空区平台值减去费米能即功函数（`LVHAR = .TRUE.`，偶极修正见 [VASP 第9章](../../VASP/tutorials/09-surface-workfunction.md)）。

---

## 9. 相图、凸包能与 Materials Project

### 9.1 新版 API 与 key

```bash
export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    # 32-character key from the MP dashboard
unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT                # remove any old .../rest/v2 endpoint (mp-api / pymatgen)
```

pymatgen 自带的 `pymatgen.ext.matproj.MPRester` 不依赖 mp-api，key 长度不是 32 位时直接报 `ValueError: Invalid or old API key`，服务器地址取自 `PMG_MAPI_ENDPOINT`（默认 `https://api.materialsproject.org/`）；也可以安装 `mp-api` 使用功能更全的 `mp_api.client.MPRester`，它读取的是 `MP_API_KEY` 和 `MP_API_ENDPOINT`。

```python
from pymatgen.ext.matproj import MPRester
with MPRester() as mpr:                       # reads PMG_MAPI_KEY from .pmgrc.yaml or env
    docs = mpr.summary_search(formula="LiFePO4",
                              _fields=["material_id", "energy_above_hull",
                                       "formation_energy_per_atom", "band_gap"])
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])
```

> [!NOTE]
> pymatgen 的 `MPRester()` 不带参数时从 pymatgen 设置中的 `PMG_MAPI_KEY` 取 key（`pmg config --add PMG_MAPI_KEY <key>` 或同名环境变量）；示例脚本中显式读取 `MP_API_KEY` 再传进去，两种方式都可以。

### 9.2 构建相图

[pmg_phase_diagram_LiFePO4.py](../examples/04-pymatgen-io/pmg_phase_diagram_LiFePO4.py)：

```python
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
pd = PhaseDiagram(entries)
for e in candidates:                                   # all LiFePO4 polymorphs
    print(e.entry_id, pd.get_e_above_hull(e), pd.get_form_energy_per_atom(e))
decomp, ehull = pd.get_decomp_and_e_above_hull(best)  # decomposition products
```

```bash
python pmg_phase_diagram_LiFePO4.py               # first run: download and cache entries
python pmg_phase_diagram_LiFePO4.py --offline     # reuse LiFePO_entries.json
python pmg_phase_diagram_LiFePO4.py --offline --my-vasprun myrun/vasprun.xml   # your own MP-compatible run
python pmg_phase_diagram_LiFePO4.py --offline --my-energy <E_total_eV> --my-formula Li4Fe4P4O16  # uncorrected, rough test only
```

| 要点 | 说明 |
|---|---|
| 能量修正 | `get_entries_in_chemsys` 返回的条目已按 MP2020 兼容方案（`MaterialsProject2020Compatibility`）修正 GGA/GGA+U 混合能量 |
| 加入自己的计算 | 必须用与 MP 一致的设置（`MPRelaxSet`、相同 U 值和 POTCAR），并对自己的 `ComputedEntry` 做同样的兼容性处理：脚本的 `--my-vasprun` 用 `Vasprun.get_computed_entry()` 生成条目，再交给 `MaterialsProject2020Compatibility().process_entry()`，设置不兼容时返回 `None` 并退出；`--my-energy` 不做修正，只能粗略参考。否则 $E_\text{hull}$ 没有意义 |
| 判据 | $E_\text{hull} = 0$ 为稳定；0–50 meV/atom 常被视为“可能合成”的亚稳相，这只是经验范围 |
| 作图 | `PDPlotter` 支持 2–4 元；四元系可画三元子系统，示例输出 `LiFeO_pd.png` |

---

## 10. pymatgen 处理 Gaussian 与 CP2K

### 10.1 Gaussian

[pmg_gaussian_input.py](../examples/04-pymatgen-io/pmg_gaussian_input.py)：

```python
from pymatgen.io.gaussian import GaussianInput, GaussianOutput
gi = GaussianInput(mol, charge=0, spin_multiplicity=1,
                   functional="B3LYP", basis_set="def2SVP",
                   route_parameters={"opt": "", "freq": "",
                                     "EmpiricalDispersion": "GD3BJ", "scf": "tight"},
                   link0_parameters={"%mem": "8GB", "%nprocshared": "16", "%chk": "benzene.chk"},
                   dieze_tag="#P")
go = GaussianOutput("benzene_opt_freq.log")
print(go.properly_terminated, go.final_energy, go.stationary_type)   # final_energy in Hartree
```

| 参数 | 说明 |
|---|---|
| `route_parameters` | 值为 `""` 写成裸关键词，字符串写成 `keyword=value` |
| `link0_parameters` | 键带 `%` |
| `dieze_tag` | `#P` 输出更详细，便于解析 |
| 空行结尾 | Gaussian 要求输入以空行结束，脚本写文件时确保末尾有 `\n\n` |
| `GaussianOutput.frequencies` | 每次频率计算一组 `[{"frequency": ..., ...}]`，负值为虚频 |

```bash
python pmg_gaussian_input.py                 # writes benzene_opt_freq.gjf
python pmg_gaussian_input.py --parse benzene_opt_freq.log
```

### 10.2 CP2K

[pmg_cp2k_input.py](../examples/04-pymatgen-io/pmg_cp2k_input.py) 用 `pymatgen.io.cp2k.sets.StaticSet` 生成水分子输入，再用 `Cp2kInput.update()` 关掉 UKS、改为非周期、删掉打印段：

```python
ci = Cp2kInput.from_str(cs.get_str())
ci.update({"FORCE_EVAL": {"STRESS_TENSOR": "NONE",
                          "DFT": {"UKS": False, "MULTIPLICITY": 1,
                                  "POISSON": {"PERIODIC": "NONE", "PSOLVER": "MT"}},
                          "SUBSYS": {"CELL": {"PERIODIC": "NONE"}}}})
ci.write_file("h2o_static.inp")
```

```bash
python pmg_cp2k_input.py
cp2k.psmp -c h2o_static.inp                  # syntax check only, no calculation
python pmg_cp2k_input.py --parse h2o.out     # Cp2kOutput: version, run type, energy (eV)
```

> [!WARNING]
> pymatgen 的 CP2K 输入集面向高通量，默认打开 UKS；非 k 点计算还会写入旧式的 `&DFT/&PRINT/&PDOS`（`DftSet` 的 `print_pdos=True` 默认值）。CP2K 2026.2 已把 PDOS 移到 `&DFT/&PRINT/&DOS` 下（内嵌 `&PDOS`），旧写法在解析阶段就会中止（`unknown subsection PDOS of section PRINT`），所以示例脚本删掉了整个 `&DFT/&PRINT` 段；需要 PDOS 时按新语法手动添加，或构造输入集时传 `print_pdos=False`。pymatgen 2026.10 的 `CellOptSet` 不写 `&MOTION/&CELL_OPT/TYPE`；但如果你沿用旧模板或手动加了 `TYPE DIRECT_CELL_OPT`，2026.2 会报 `The specified keyword 'TYPE' is not available anymore`。生成后一定先用 `cp2k.psmp -c` 检查，变化细节见 [CP2K 第7章](../../CP2K/tutorials/07-geo-cell-opt.md) 和 [CP2K 第10章](../../CP2K/tutorials/10-dos-pdos-band.md)。

---

## 11. 跨程序输出解析：ASE、cclib、dpdata

### 11.1 一个入口读四款软件

[parse_all_outputs.py](../examples/04-output-parsing/parse_all_outputs.py) 自动识别文件类型，统一返回能量（eV）、最大力和末结构（ASE `Atoms`），可汇总写成 extxyz：

```bash
cd Common/examples/04-output-parsing
python parse_all_outputs.py relax/OUTCAR g16/h2o_freq.log cp2k/h2o.out omx/si.out --extxyz final.xyz
```

| 程序 | 读取方式 | 得到 |
|---|---|---|
| VASP | `read("OUTCAR", format="vasp-out", index=-1)` 或 `vasp-xml` | 能量（sigma→0）、力、应力、结构 |
| Gaussian | `read(log, format="gaussian-out", index=-1)`；cclib 补充频率、HOMO/LUMO | 能量、力（需 `force`）、结构 |
| CP2K | 主输出用正则匹配 `ENERGY| Total FORCE_EVAL`（2025.1 起单位写 `[hartree]`，此前为 `[a.u.]:`）；结构读 `PROJECT-pos-1.xyz` | 能量、末结构 |
| OpenMX | `ase.calculators.openmx.reader.read_openmx(label)` | 能量（`Utot.`）、力、结构 |

### 11.2 cclib：Gaussian 输出的通用解析

```python
import cclib
data = cclib.io.ccread("h2o_freq.log")
print(data.scfenergies[-1])            # eV in cclib 1.8.1
print(data.vibfreqs)                   # cm^-1; negative = imaginary
print(data.moenergies[0][data.homos[0]])   # HOMO energy (eV)
```

```bash
pip install "cclib==1.8.1"
ccget scfenergies vibfreqs h2o_freq.log     # command-line quick look
ccwrite xyz h2o_opt.log                     # last geometry -> xyz
```

### 11.3 dpdata：为机器学习势准备数据

[dpdata_convert.py](../examples/04-output-parsing/dpdata_convert.py) 用 `dpdata.LabeledSystem` 读带标签（能量、力、维里）的帧，抽样并划分训练/测试集：

```bash
python dpdata_convert.py aimd_*/OUTCAR --fmt vasp/outcar --stride 5 --test-frac 0.1 --out dataset
```

| `--fmt` | 来源 |
|---|---|
| `vasp/outcar`、`vasp/xml` | VASP 优化或 AIMD |
| `cp2k/output`、`cp2k/aimd_output` | CP2K 单点/优化主输出；AIMD 目录 |
| `gaussian/log` | 带 `force` 关键词的 Gaussian log |
| `openmx/md` | OpenMX MD |

输出 `deepmd/npy`（DeePMD-kit）和 `extxyz`（MACE、NequIP、GPUMD 等）。检查：脚本打印每个输入的帧数、原子数和能量范围；能量范围异常大通常说明混入了 SCF 未收敛的帧。

> [!IMPORTANT]
> 同一训练集里的所有帧必须来自**相同的 DFT 设置**（泛函、截断能、赝势、smearing）；混用不同设置的能量会让势函数学到虚假的能量偏移。

---

## 12. 作业状态检查与数据归档

### 12.1 判断作业是否真正完成

| 程序 | 正常结束标志 | 收敛/错误标志 | 位置 |
|---|---|---|---|
| VASP | `General timing and accounting` | 优化收敛：`reached required accuracy`；SCF 步数等于 NELM 说明电子步未收敛；`ZBRENT: fatal error`、`Error EDDDAV`、`VERY BAD NEWS` 等报错 | OUTCAR；报错多写在标准输出（`vasp.out`、`slurm-*.out`） |
| Gaussian | `Normal termination of Gaussian` | `Error termination via Lnk1e in .../l502.exe`（链接号指示出错模块）；`-- Stationary point found.` | .log 末尾 |
| CP2K | `PROGRAM ENDED AT` | `SCF run NOT converged`；`ABORT`；`GEOMETRY OPTIMIZATION COMPLETED` | 主输出 |
| OpenMX | `.out` 中出现 `Elapsed.Time.`；标准输出末尾 `The calculation was normally finished.` | 最后一个 `SCF history at MD=` 块的 `SCF=` 行数达到 `scf.maxIter`（默认 40）即 SCF 未收敛；MD/优化时 `.md` 中每步的 `scf_conv=0/1` | `.out`、`.md` 与标准输出（`.std`；ASE 运行时为 `.log`） |

```bash
grep -L "reached required accuracy" */OUTCAR          # relaxations NOT converged
grep -l "Error termination" */*.log                   # failed Gaussian jobs
grep -l "SCF run NOT converged" */*.out               # CP2K SCF failures
```

[job_status_scanner.py](../examples/04-output-parsing/job_status_scanner.py) 只用标准库，递归扫描目录，把上述规则汇总成表格，可导出 CSV。目录中的 `mock_runs/` 是**人为构造的标记文件**，只用来测试扫描器：

```bash
python job_status_scanner.py mock_runs
python job_status_scanner.py ~/projects/TiO2 --csv status.csv
```

（示意输出，数值来自模拟文件，没有物理意义）

```text
path                   code     status          energy(eV) steps  note
cp2k_scf_fail/h2o.out  CP2K     ERROR                                SCF run NOT converged;
g16_error/bad.log      Gaussian ERROR                                failed in l502
vasp_zbrent            VASP     ERROR           -50.10000000         ZBRENT (ionic step bracketing failed);
```

日常使用中，与调度器结合的作业监控可用辅助脚本 [job_status.py](../scripts/job_status.py)（`python ../scripts/job_status.py -h` 查看用法），作业脚本可用 [mkjob.py](../scripts/mkjob.py) 生成，调度器细节见 [第2章](02-job-scheduler.md)。自动纠错重跑（custodian）和工作流见 [第6章](06-workflow-automation.md)。

### 12.2 归档：ase db、CSV、HDF5

[ase_db_collect.py](../examples/04-output-parsing/ase_db_collect.py) 遍历项目目录，把**已正常结束**的 VASP/Gaussian 计算写入 ase db，并导出 CSV 和 extxyz：

```bash
python ase_db_collect.py ~/projects/TiO2 --db tio2.db --csv tio2.csv --extxyz tio2.xyz
ase db tio2.db code=vasp -c id,formula,energy,max_force,dir
```

需要 HDF5 时，可以在 pandas 中把汇总表写成 HDF5（需要 `tables` 包）：

```python
import pandas as pd
df = pd.read_csv("tio2.csv")
df.to_hdf("tio2.h5", key="runs", mode="w")        # requires PyTables (pip install tables)
print(pd.read_hdf("tio2.h5", "runs").head())
```

| 方式 | 适合 | 注意 |
|---|---|---|
| ase db（SQLite） | 结构 + 能量 + 自定义键，可按条件查询 | 保留键名不能自定义；并发写入用 PostgreSQL |
| CSV | 作图、表格、Excel | 不保存结构 |
| HDF5 | 大量数组（DOS、轨迹） | 需 `tables`/`h5py`；与他人交换时说明字段含义 |
| extxyz | 结构 + 能量/力，MLIP 通用 | 单位 eV、Å |

### 12.3 公开数据：NOMAD 与 OPTIMADE

- **NOMAD**（<https://nomad-lab.eu>）：注册账户后，把计算目录打包成 zip 在网页上传，NOMAD 会自动识别 VASP、Gaussian、CP2K、OpenMX 等程序的输出并生成可检索条目，发表时可获得 DOI。也可用官方的 `nomad-utility-workflows` 包在 Python 中上传（如 `upload_files_to_nomad(filename=..., url="test")`），官方教程建议先在 test 部署上练习，生产部署上的发布不可撤销。
- **OPTIMADE**：Materials Project、OQMD、Alexandria、COD、NOMAD 等数据库实现的统一查询接口。[pmg_optimade_query.py](../examples/04-pymatgen-io/pmg_optimade_query.py) 用 pymatgen 的 `OptimadeRester` 同时查询多个数据库并保存 CIF：

```bash
python pmg_optimade_query.py --providers mp oqmd --formula TiO2 --max 5
python pmg_optimade_query.py --filter 'elements HAS ALL "Ga","N" AND nsites<=8'
```

OPTIMADE 中的 `chemical_formula_hill` 用 Hill 顺序（TiO2 写成 `O2Ti`），脚本会自动转换。不同数据库的结构来源和泛函不同，下载的结构可以用作初始结构，能量不能直接混用。

> [!CAUTION]
> 打包上传或公开分享前删除 POTCAR，并确认 `vaspout.h5` 中不含赝势数据；WAVECAR、CHGCAR 体积大，一般不需要上传。Gaussian、VASP 的可执行文件和源码同样不能分发。

---

## 参数选择与收敛建议

| 项目 | 建议 |
|---|---|
| ASE 优化 `fmax` | 普通结构 0.02–0.05 eV/Å；声子、弹性、NEB 端点 ≤ 0.01 eV/Å；同时检查能量变化 |
| 用 ASE 还是程序自身优化器 | VASP/Gaussian/CP2K 的内置优化器效率更高，优先用；ASE 优化器用于特殊约束、MLIP 混合流程 |
| EOS | 7–11 点，±3% 晶格常数；换 2 种 EOS 形式比较 $V_0$；DFT 下用足够高的 ENCUT |
| NEB | 显式 `method`；IDPP 初始路径；先普通 NEB 后 CI；图像 3–7 个 |
| pymatgen 输入集 | 用 `potcar_spec=True` 先检查；所有覆盖写进 `user_incar_settings`，并在论文里写明 pymatgen 版本 |
| 与 MP 数据比较 | 完全使用 `MPRelaxSet`（或 MP24）设置与相同的能量修正方案 |
| 解析中断的作业 | `Vasprun(..., exception_on_bad_xml=False)` 读能量，但仍要用 OUTCAR 判断是否收敛 |
| 状态判断 | 不要只看“有没有能量”：同时检查正常结束标志、SCF 收敛、离子步收敛 |
| 版本固定 | 用 `environment.yml` 锁定 `ase`、`pymatgen`、`cclib==1.8.1` 版本（[第1章](01-linux-hpc.md)） |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `Please set either command in calculator or one of the following environment variables (prioritized as follows): ASE_VASP_COMMAND, VASP_COMMAND, VASP_SCRIPT` | ASE 不知道如何运行 VASP | `export ASE_VASP_COMMAND="mpirun -np 16 vasp_std"`，或写进 config.ini 的 `[environment]` |
| `RuntimeError: VASP_PP_PATH does not exist: ...` 或找不到 POTCAR | 路径写错，或指向了 `potpaw_PBE` 本身 | `VASP_PP_PATH` 应指向**包含** `potpaw_PBE/` 的上一级目录 |
| `EnvironmentError: No configuration of openmx.  Missing section [openmx] in configuration` | ASE ≥ 3.23 的文件型计算器找不到命令 | 设置 `ASE_OPENMX_COMMAND`（或其它程序的 `ASE_<NAME>_COMMAND`），或在 config.ini 中配置 |
| NEB 运行时出现 `UserWarning: The default method has changed from 'aseneb' to 'improvedtangent'...` | ASE ≥ 3.28 未显式指定 `method` | 写 `NEB(images, method="improvedtangent")`；复现旧结果用 `"aseneb"` |
| `ValueError: Bad key: energy` | ase db 中用了保留键名（`energy`、`fmax`、`id`、元素符号等） | 改名，如 `max_force`、`E_ads`、`metal` |
| `ImportError: cannot import name 'Structure' from 'pymatgen'` | 旧写法 | `from pymatgen.core import Structure` |
| `No POTCAR for Fe_pv with functional PBE found. Please set the PMG_VASP_PSP_DIR ...`，新版为 `FileNotFoundError: You do not have the right POTCAR with functional='PBE' and symbol='Fe_pv' ...` | 未配置 pymatgen 的 POTCAR 目录或缺该版本 | `pmg config -p` 整理后 `pmg config --add PMG_VASP_PSP_DIR`；需要 PBE_54/PBE_64 时设 `PMG_DEFAULT_FUNCTIONAL` 或 `user_potcar_functional` |
| `ValueError: Invalid or old API key. Please obtain an updated API key at https://materialsproject.org/dashboard.` / HTTP 403 | 旧版 MP key 或旧地址（旧 API 2025-09 下线） | 用 dashboard 上的 32 位新 key；`unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT` |
| `xml.etree.ElementTree.ParseError` 解析 vasprun.xml | 作业被杀，XML 不完整 | `Vasprun(..., exception_on_bad_xml=False)`，或用 `Outcar`/`Oszicar`；然后检查作业是否真正收敛 |
| `Oszicar.final_energy` 为 `nan` 并有警告 | OSZICAR 中出现 `******` 溢出 | 体系能量发散或格式溢出，查看 OUTCAR 中对应离子步 |
| 能带图中半导体被画成金属、带隙为 0 | 用了自洽计算的费米能或 vasprun 不是 line 模式 | `get_band_structure(..., line_mode=True, efermi="smart")`，读非自洽能带的 vasprun.xml |
| `ValueError: invalid literal for int() with base 10: 'electrons'`（cclib） | cclib 2.0a1 解析 Gaussian 16 的 bug | `pip install "cclib==1.8.1"` |
| `ValueError: Data have different keys! Maybe one is spin-polarized and the other is not?` | Chgcar 相减时一个是自旋极化、一个不是 | 子体系与整体用相同 `ISPIN` 重算 |
| CP2K 2026.2 运行 pymatgen 生成的输入时报 `unknown subsection PDOS of section PRINT` | pymatgen 输入集默认写旧式 `&DFT/&PRINT/&PDOS` | 删掉该段或改写成 `&DFT/&PRINT/&DOS` 内嵌 `&PDOS`，见 [CP2K 第10章](../../CP2K/tutorials/10-dos-pdos-band.md) |
| CP2K 2026.2 报 `The specified keyword 'TYPE' is not available anymore` | 旧模板中的 `&MOTION/&CELL_OPT/TYPE`（pymatgen 2026.10 的 `CellOptSet` 本身不写这一行） | 删除该关键词即可（晶胞优化固定为直接法），见 [CP2K 第7章](../../CP2K/tutorials/07-geo-cell-opt.md) |
| 状态扫描显示 VASP `RUNNING/KILLED` 但作业已不在队列中 | 超时或节点故障，OUTCAR 没有 `General timing` | 检查 `slurm-*.out`；从 CONTCAR 续算 |
| `/bin/bash^M: bad interpreter` 或 `sbatch: error: Batch script contains DOS line breaks (\r\n)` | 在 Windows 中编辑了脚本 | `dos2unix job.sh *.py`，见[第1章](01-linux-hpc.md) |

---

## 练习

1. 运行 `04-ase-emt-basics` 中的四个脚本。把 `ase_eos_cu.py` 的 `--eos` 分别换成 `vinet`、`murnaghan`，比较 $a_0$ 和 $B$ 的差别。
2. 在 `ase_neb_idpp.py` 中分别使用 `--method improvedtangent` 和 `--method aseneb`，并把 IDPP 改为线性插值（`neb.interpolate()`），比较收敛步数和势垒。
3. 用 `ase_vasp_relax.py --dry-run` 和 `pmg_write_inputset.py` 分别为 Si 生成输入，对比两者的 INCAR、KPOINTS 和 POTCAR 选择，说明差异来自哪里。
4. 有 VASP 许可证的同学：用 `--mode vasp` 和 `--mode ase` 各优化一次 Si，比较晶格常数、离子步数和总耗时。
5. 用 `pmg_phase_diagram_LiFePO4.py` 找出 LiFePO4 各多形体的 $E_\text{hull}$，再查 MP summary 中的 `energy_above_hull`，两者是否一致？如果不一致，可能的原因是什么？（提示：`get_entries_in_chemsys` 只保留经 MP2020 方案修正的 GGA/GGA+U 条目，而网站 summary 中的 `energy_above_hull` 可能来自混合了 r2SCAN 数据的凸包，请查阅 MP 文档中 thermo type 的说明。）
6. 在 `04-output-parsing/mock_runs` 中新建一个目录，按 Gaussian “Normal termination” 和 CP2K “PROGRAM ENDED AT” 的格式写入标记行，运行 `job_status_scanner.py`，确认能被正确识别。
7. 用 `pmg_optimade_query.py` 从两个数据库下载 GaN 结构，用 `SpacegroupAnalyzer` 比较它们的空间群和晶格常数。

---

## 参考资料

官方文档：

- ASE 文档（2026-06 起迁到 docs.ase-lib.org）：<https://docs.ase-lib.org/>；版本说明：<https://docs.ase-lib.org/releasenotes.html>
- ASE 计算器与配置文件：<https://docs.ase-lib.org/ase/calculators/calculators.html>
- pymatgen 文档：<https://pymatgen.org/>；更新日志：<https://pymatgen.org/CHANGES.html>；`pymatgen.io.vasp` API：<https://pymatgen.org/pymatgen.io.vasp.html>
- Materials Project API 入门（新 key）：<https://docs.materialsproject.org/downloading-data/using-the-api/getting-started>
- cclib 文档：<https://cclib.github.io/>
- dpdata 仓库与文档：<https://github.com/deepmodeling/dpdata>
- VASP Wiki：vaspout.h5：<https://vasp.at/wiki/Vaspout.h5>
- NOMAD API 上传教程：<https://nomad-lab.eu/prod/v1/docs/tutorial/upload_publish_api.html>
- OPTIMADE：<https://www.optimade.org/>
- spglib 文档：<https://spglib.readthedocs.io/>

相关 GitHub 仓库（只链接，不复制代码）：

- [ase/ase（GitLab）](https://gitlab.com/ase/ase) — ASE 官方仓库，统一的 Atoms 与计算器接口。
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — 输入集、解析、相图、MP 客户端；2026 年核心拆为 pymatgen-core。
- [materialsproject/api](https://github.com/materialsproject/api) — mp-api 客户端。
- [computron/pymatgen_tutorials](https://github.com/computron/pymatgen_tutorials) — pymatgen 作者的视频 + notebook 教程。
- [materialsproject/workshop](https://github.com/materialsproject/workshop) — MP 官方 workshop 课程。
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book) — 经典的 ASE + VASP 实战书。
- [cclib/cclib](https://github.com/cclib/cclib) — 量化程序输出解析。
- [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata) — 多程序数据格式互转。
- [ACEsuit/mace](https://github.com/ACEsuit/mace) — 可作为 ASE 计算器的通用机器学习势。
- [nglviewer/nglview](https://github.com/nglviewer/nglview) — 在 Jupyter 中显示 ASE/pymatgen 结构与轨迹。
- [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools) — CP2K 输入校验与 JSON/YAML 互转。
- [Quantum-Accelerators/quacc](https://github.com/Quantum-Accelerators/quacc) — 基于 ASE 的多程序工作流。
- [PytLab/VASPy](https://github.com/PytLab/VASPy) — 轻量的 VASP 文件读写库（中文开发者）。

关键文献：

- A. H. Larsen *et al.*, The atomic simulation environment—a Python library for working with atoms, *J. Phys.: Condens. Matter* **29**, 273002 (2017). DOI: 10.1088/1361-648X/aa680e
- S. P. Ong *et al.*, Python Materials Genomics (pymatgen), *Comput. Mater. Sci.* **68**, 314 (2013). DOI: 10.1016/j.commatsci.2012.10.028
- A. Jain *et al.*, The Materials Project, *APL Mater.* **1**, 011002 (2013). DOI: 10.1063/1.4812323
- A. Wang *et al.*, A framework for quantifying uncertainty in DFT energy corrections, *Sci. Rep.* **11**, 15496 (2021). DOI: 10.1038/s41598-021-94550-5
- G. Henkelman, B. P. Uberuaga, H. Jónsson, A climbing image nudged elastic band method, *J. Chem. Phys.* **113**, 9901 (2000). DOI: 10.1063/1.1329672
- G. Henkelman, H. Jónsson, Improved tangent estimate in the nudged elastic band method, *J. Chem. Phys.* **113**, 9978 (2000). DOI: 10.1063/1.1323224
- S. Smidstrup *et al.*, Improved initial guess for minimum energy path calculations, *J. Chem. Phys.* **140**, 214106 (2014). DOI: 10.1063/1.4878664
- E. Bitzek *et al.*, Structural relaxation made simple, *Phys. Rev. Lett.* **97**, 170201 (2006). DOI: 10.1103/PhysRevLett.97.170201
- P. Haas, F. Tran, P. Blaha, Calculation of the lattice constant of solids with semilocal functionals, *Phys. Rev. B* **79**, 085104 (2009). DOI: 10.1103/PhysRevB.79.085104
- N. M. O'Boyle *et al.*, cclib: a library for package-independent computational chemistry algorithms, *J. Comput. Chem.* **29**, 839 (2008). DOI: 10.1002/jcc.20823
- C. W. Andersen *et al.*, OPTIMADE, an API for exchanging materials data, *Sci. Data* **8**, 217 (2021). DOI: 10.1038/s41597-021-00974-z
- C. Draxl, M. Scheffler, The NOMAD laboratory: from data sharing to artificial intelligence, *J. Phys. Mater.* **2**, 036001 (2019). DOI: 10.1088/2515-7639/ab13bb

[← 上一章](03-structure-modeling.md) | [返回目录](../README.md) | [下一章 →](05-postprocess-visualization.md)
