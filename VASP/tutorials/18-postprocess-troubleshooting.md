# 第18章 后处理生态、自动化与报错排查（VASP Post-processing, Automation & Troubleshooting）

本章把散落在前面各章的后处理与运维知识系统化：VASPKIT 1.5.1、pymatgen、ASE、py4vasp 四大后处理工具各自擅长什么、怎样用命令行和 Python 脚本批量提取数据；VESTA、OVITO、sumo、pyprocar 等可视化工具的分工；Shell 批处理、Slurm 作业数组、custodian 自动纠错和 atomate2 工作流的自动化方案；以及一套"错误信息 → 原因 → 修复"的 VASP 报错排查方法和对照表。最后说明公开仓库中处理 POTCAR、大文件和 AI 工具的合规要求（VASP 2026 年 AI 工具使用规范）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据任务选择合适的后处理工具，知道它们读取哪些文件、对 VASP 版本有什么要求；
- 用 VASPKIT 的命令行/管道模式批量处理几十个计算目录；
- 用 pymatgen 解析 vasprun.xml/OUTCAR、生成标准输入集，用 ASE 驱动 VASP 和做 NEB 插值，用 py4vasp 读取 vaspout.h5；
- 用作业数组、依赖链、custodian、atomate2 把重复劳动自动化；
- 遇到报错时按固定流程定位原因，并能从 OUTCAR 的 WARNING/ADVICE 中读出 VASP 的建议；
- 正确配置 `.gitignore`，遵守 VASP 许可证与 AI 工具使用规范。

| 场景 | 推荐工具 | 本章小节 |
|---|---|---|
| 不想写代码，快速出 DOS/能带/功能性数据 | VASPKIT | §3 |
| 写分析脚本、批量入库、生成输入集 | pymatgen | §4 |
| Python 中直接调用 VASP、结构建模、NEB 插值 | ASE | §5 |
| VASP 6.x 新输出（vaspout.h5）、KPOINTS_OPT、STM、MLFF 误差 | py4vasp | §6 |
| 出版级能带/DOS 图 | sumo、pyprocar、自写 matplotlib | §7 |
| 几十到上千个计算 | Shell + 作业数组、custodian、atomate2 | §8 |
| 计算报错或结果异常 | 排查流程 + 对照表 | §9、常见问题与报错 |
| 公开代码与数据 | `.gitignore`、合规清单 | §10 |

配套算例：[examples/18-postprocess-toolkit](../examples/18-postprocess-toolkit/README.md)。

---

## 2. 后处理生态地图

| 工具 | 当前版本（2026-11） | 主要读取 | 擅长 | 安装 | 许可 |
|---|---|---|---|---|---|
| VASPKIT | 1.5.1 | POSCAR、DOSCAR、PROCAR、EIGENVAL、CHGCAR、LOCPOT、OUTCAR、XDATCAR | 一键式菜单：输入生成、k 路径、DOS/能带、电荷/势、光学、MD、弹性 | SourceForge 下载二进制 | 学术免费 |
| pymatgen | 2026.10.24（自动依赖 pymatgen-core） | vasprun.xml、OUTCAR、CHGCAR、LOCPOT、PROCAR | 结构与分析 API、输入集（MPRelaxSet 等）、相图、与 Materials Project 对接 | `pip install pymatgen` | MIT 类（见仓库） |
| ASE | 3.29.0 | 几乎所有结构格式，vasprun.xml/OUTCAR | 计算器接口、建模、NEB、优化器、EOS | `pip install ase` | LGPL |
| py4vasp | 0.11.3 | vaspout.h5（需 `-DVASP_HDF5`） | VASP 官方：能带/DOS/声子/介电/STM/MLFF/电声 | `pip install py4vasp` | Apache-2.0 |
| sumo | 3.0.0 | vasprun.xml | 出版级能带、DOS、光学、声子图；k 路径生成 | `pip install sumo` | MIT |
| pyprocar | 6.5.0 | PROCAR、vasprun.xml 等 | 投影能带、自旋纹理、费米面、反折叠 | `pip install pyprocar` | GPL-3.0 |
| LobsterPy | 0.6.1 | LOBSTER 输出 | 自动成键分析（第16章） | `pip install lobsterpy` | BSD-3 |
| VESTA | — | CHGCAR、PARCHG、ELFCAR、LOCPOT、POSCAR | 晶体结构与体积数据可视化 | 官网下载 | 免费 |
| OVITO | — | XDATCAR、POSCAR、extxyz | MD 轨迹可视化与分析 | 官网下载 / `pip install ovito` | 基础版免费 |
| p4vasp | （停更） | vasprun.xml | 老牌 GUI，Python 2 时代 | — | GPL-2.0 |
| rsgrad | — | OUTCAR、PROCAR 等 | Rust 编写的快速 OUTCAR 监控（收敛、受力、磁矩） | GitHub 发布页 | MIT |

> [!TIP]
> 选择原则：**交互式快速出结果用 VASPKIT；需要可重复、可批量的分析用 Python（pymatgen/ASE/py4vasp）；VASP 6.4 以后新增的输出（KPOINTS_OPT、部分电荷 h5、声子 h5、MLFF）优先用 py4vasp。** 同一个量最好用两种工具交叉验证一次（例如 VASPKIT 911 与 pymatgen 的带隙）。

本仓库自带的脚本（由各章作者编写，用 `-h` 查看用法）：[plot_band.py](../scripts/plot_band.py)、[plot_dos.py](../scripts/plot_dos.py)、[chgdiff.py](../scripts/chgdiff.py)、[locpot_workfunction.py](../scripts/locpot_workfunction.py)、optics_post.py、[elastic_post.py](../scripts/elastic_post.py)、[md_analysis.py](../scripts/md_analysis.py)、[check_vasp.py](../scripts/check_vasp.py)、[vasp_inputgen.py](../scripts/vasp_inputgen.py)、[conv_test.sh](../scripts/conv_test.sh)、job_vasp_slurm.sh 等。通用的 ASE/pymatgen 用法见 Common 第4章，跨软件的作图见 Common 第5章。

---

## 3. VASPKIT 1.5.1

### 3.1 安装与配置

```bash
tar -zxvf vaspkit.1.5.1.linux.x64.tar.gz      # 文件名以下载页为准
cd vaspkit.1.5.1
cp -f how_to_set_environment_variable ~/.vaspkit
vi ~/.vaspkit        # 至少修改 PBE_PATH（赝势库路径）；作图需设置 PYTHON_BIN 并打开 PLOT_MATPLOTLIB
echo 'export PATH=$HOME/vaspkit.1.5.1/bin:$PATH' >> ~/.bashrc
```

`~/.vaspkit` 中 `PBE_PATH` 指向 PBE 赝势库，`RECOMMENDED_POTCAR = .TRUE.` 时 103 会按推荐后缀拼接 POTCAR。

### 3.2 四种运行方式

| 方式 | 示例 | 适用 |
|---|---|---|
| 交互菜单 | `vaspkit` → 输入 `11` → `113` | 第一次使用、探索功能 |
| 命令行参数 | `vaspkit -task 102 -kpr 0.04` | 部分功能支持（`vaspkit -help` 查看） |
| 管道 | `echo -e "303\n" \| vaspkit` | 所有功能，按菜单顺序依次输入 |
| 输入文件 | `vaspkit < cmd.in` | 输入较多时更清晰 |

### 3.3 常用功能

| 一级编号 | 内容 | 典型子功能 |
|---|---|---|
| 01 | 输入文件生成 | 102 KPOINTS（按 2π×Å⁻¹ 分辨率）、103 POTCAR、105 cif→POSCAR |
| 03 | k 路径 | 302 二维、303 体相（写 KPATH.in 与 PRIMCELL.vasp）、304 Wannier90、305 phonopy |
| 11 | 态密度 | 111 总 DOS、113 分元素、115 指定原子与轨道 |
| 21 | 能带 | 211 能带、213 投影能带；25x 为杂化泛函能带 |
| 31 | 电荷 | 312 自旋密度、314 差分电荷、316 平面平均、325 STM |
| 42 | 势 | 426 平面平均势、427 宏观平均、428/429 格式转换 |
| 71 | 光学 | 711 体相线性光学、710 二维 |
| 72 | 分子动力学 | 722 MSD、723 扩散系数、726 RDF、728 振动 DOS |
| 91 | 半导体 | 911 带隙、912 有效质量 |
| 92 | 二维材料 | 921 居中、922 改真空层、927 相对真空能级的带边 |
| 95 | 声子 | 1.5.0 起新增从 OUTCAR 读取声子相关数据 |

完整速查表见 [vaspkit_cli_cheatsheet.md](../examples/18-postprocess-toolkit/vaspkit_cli_cheatsheet.md)。

### 3.4 批处理

```bash
# 每个子目录提取带隙
for d in */; do
  gap=$(cd "$d" && echo -e "911\n" | vaspkit 2>/dev/null | grep -i "band gap" | head -1)
  echo "$d  $gap"
done

# 每个子目录导出分元素 PDOS（生成 PDOS_*.dat）
for d in dos_*/; do (cd "$d" && echo -e "11\n113\n" | vaspkit > /dev/null); done
```

> [!WARNING]
> VASPKIT 的屏幕输出格式会随版本变化，用 `sed -n '40p'` 这类按行号截取的写法在升级后容易失效。脚本化提取时优先读取 VASPKIT 写出的数据文件（如 `BAND_GAP`、`PDOS_*.dat`），或改用 Python 工具。

---

## 4. pymatgen

> [!NOTE]
> 2026 年起 pymatgen 拆分为 `pymatgen-core`（核心数据结构与各软件的文件读写）和 `pymatgen`（分析与应用）。`pip install pymatgen` 会自动安装两者，`from pymatgen.io.vasp import Vasprun` 等导入路径不变。只需要读写文件时可以只装 `pymatgen-core`。

### 4.1 解析输出

```python
from pymatgen.io.vasp import Vasprun, Outcar

vr = Vasprun("vasprun.xml", parse_potcar_file=False)   # 不去找 POTCAR
print(vr.converged_electronic, vr.converged_ionic)     # 两者都为 True 才能用
print(vr.final_energy, vr.efermi)
gap, cbm, vbm, is_direct = vr.eigenvalue_band_properties
st = vr.final_structure                                 # pymatgen Structure

oc = Outcar("OUTCAR")
print(oc.total_mag, oc.run_stats["Elapsed time (sec)"])
print([m["tot"] for m in oc.magnetization])             # 每个原子的磁矩（需要 LORBIT）
```

`Outcar` 还提供 `final_energy`（energy(sigma→0)）、`final_fr_energy`（TOTEN）、`elastic_tensor`（kBar）、`drift` 等属性；读取某些特定块需要调用 `read_*` 方法，具体见 API 文档。

### 4.2 能带与 DOS 作图

```python
import matplotlib
matplotlib.use("Agg")
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

bs = Vasprun("band/vasprun.xml", parse_potcar_file=False).get_band_structure(
        kpoints_filename="band/KPOINTS", line_mode=True)
print(bs.get_band_gap())                  # {'energy':..., 'direct':..., 'transition':...}
BSPlotter(bs).get_plot(zero_to_efermi=True, ylim=(-6, 6)).figure.savefig("band.png", dpi=200)

cdos = Vasprun("dos/vasprun.xml", parse_potcar_file=False).complete_dos
p = DosPlotter(sigma=0.05)
p.add_dos("Total", cdos)
p.add_dos_dict(cdos.get_element_dos())
p.get_plot(xlim=(-8, 6)).figure.savefig("dos.png", dpi=200)
print(cdos.get_band_center())             # d-band center（默认 d 轨道，相对 E_F）
```

完整脚本：[pymatgen_parse_vasprun.py](../examples/18-postprocess-toolkit/pymatgen_parse_vasprun.py)。

### 4.3 输入集

```python
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPNonSCFSet

s = Structure.from_file("POSCAR")
MPRelaxSet(s, user_incar_settings={"ENCUT": 520, "EDIFFG": -0.01}).write_input("relax")
static = MPStaticSet.from_prev_calc("relax")          # 读上一步结果继承结构与磁矩
static.write_input("static")
```

- 写 POTCAR 需要先配置赝势库：`pmg config -p <解压的赝势目录> <整理后的目录>`，再 `pmg config --add PMG_VASP_PSP_DIR <整理后的目录>`；
- Materials Project 的输入集是为高通量一致性设计的，参数（ENCUT、k 点密度、U 值）不一定是你体系的最优选择，正式研究仍需第3章的收敛测试；
- 输入集类在不同版本间有默认参数调整，论文中应注明 pymatgen 版本。

---

## 5. ASE

### 5.1 读写与轨迹

```python
from ase.io import read, write

atoms = read("CONTCAR")                       # 自动识别 VASP 格式
traj = read("OUTCAR", index=":")              # 每个离子步的结构、能量、受力
write("relax.extxyz", traj)                   # 带能量/受力的 extxyz，可直接给 OVITO 或 MLIP 训练
md = read("XDATCAR", index="::10")            # 每 10 帧取一帧
write("POSCAR_new", atoms, format="vasp", direct=True, sort=True)
```

### 5.2 Vasp 计算器

```bash
export ASE_VASP_COMMAND="mpirun -np 16 vasp_std"
export VASP_PP_PATH=/path/to/potentials      # 其下有 potpaw_PBE/ 等目录
```

```python
from ase.build import bulk
from ase.calculators.vasp import Vasp

atoms = bulk("Si", "diamond", a=5.431)
atoms.calc = Vasp(directory="si_scf", xc="PBE", encut=400, kpts=(8, 8, 8), gamma=True,
                  ismear=0, sigma=0.05, prec="Accurate")
print(atoms.get_potential_energy())          # 写输入、运行 VASP、读结果
```

- 关键词与 INCAR 标签同名（小写）；`setups={"Ga": "_d"}` 选择 POTCAR 后缀；
- ASE 会按元素重新排序原子并写 `ase-sort.dat`，在 ASE 外部处理该目录时要注意原子顺序；
- `calc.write_input(atoms)` 只写输入不运行（同样需要 `VASP_PP_PATH`）。

### 5.3 NEB 插值

```python
from ase.io import read, write
from ase.mep import NEB          # ASE >= 3.23；更早版本为 ase.neb

ini, fin = read("POSCAR_IS"), read("POSCAR_FS")
images = [ini] + [ini.copy() for _ in range(5)] + [fin]
NEB(images, method="improvedtangent").interpolate("idpp")
for i, im in enumerate(images):
    write(f"{i:02d}/POSCAR", im, format="vasp", direct=True)
```

IDPP 插值比线性插值更不容易让原子靠得过近，生成的 `00…06` 目录可直接用于 VTST NEB（[第11章](11-neb-dimer.md)）。完整示例见 [ase_vasp_calculator_example.py](../examples/18-postprocess-toolkit/ase_vasp_calculator_example.py)。

---

## 6. py4vasp

### 6.1 前提

- VASP 编译时 `CPP_OPTIONS += -DVASP_HDF5` 并链接 HDF5（[第1章](01-setup-compile.md)），运行后得到 `vaspout.h5`；
- py4vasp 要求 VASP > 6.2，官方建议使用与 VASP 发布时间相近的 py4vasp 版本；
- 只当库使用、不需要可视化依赖时可以安装 `py4vasp-core`。

### 6.2 基本用法

```python
import py4vasp

calc = py4vasp.Calculation.from_path(".")     # 或 Calculation.from_file("backup.h5")
calc.structure.to_POSCAR()                     # 结构字符串；to_ase() 转 ASE
calc.energy.read()                             # {'free energy TOTEN':..., 'energy without entropy':..., 'energy(sigma->0)':...}
calc.energy[:].plot()                          # 优化/MD 全部离子步
calc.band.plot()                               # 能带；"kpoints_opt" 选择 KPOINTS_OPT 路径
calc.band.plot("Ga(d), As(p)")                 # 投影能带（LORBIT = 11）
calc.dos.plot("s, p, d")                       # 投影 DOS
calc.bandgap.print()
calc.dos.to_csv(filename="dos.csv")            # 导出数据
calc.dos.to_image(filename="dos.png")          # 导出图片（需要 kaleido）
```

选择语法：原子用元素名或编号（`Si`、`1`、`1:4`），轨道用 `s`、`px`、`d` 等，自旋用 `up`/`down`/`total`，可以嵌套和逗号分隔（如 `Sr(s, p)`），还可以用 `kpoints_opt`/`kpoints_wan` 选择不同的 k 点集合。不确定时调用 `selections()`。

### 6.3 常用量

| 属性 | 内容 | 相关章节 |
|---|---|---|
| `structure`、`energy`、`force`、`stress` | 结构与能量、受力、应力（可按离子步切片） | 第4章 |
| `band`、`dos`、`bandgap`、`projector` | 能带、DOS、带隙、投影信息 | 第5章 |
| `density`、`potential`、`workfunction` | 电荷密度、局域势、功函数 | 第9、16章 |
| `partial_density` | 部分电荷与 STM（6.5+，`LPARDH5`） | 第16章 |
| `phonon.band`、`phonon.dos`、`force_constant`、`born_effective_charge` | 声子与介电相关 | 第12章 |
| `dielectric_function`、`dielectric_tensor` | 光学与介电 | 第14章 |
| `elastic_modulus`、`piezoelectric_tensor` | 弹性与压电 | 第13章 |
| `electron_phonon.*` | 电声耦合（6.5+） | — |

完整笔记本：[py4vasp_quickstart.ipynb](../examples/18-postprocess-toolkit/py4vasp_quickstart.ipynb)。

---

## 7. 可视化

| 需求 | 工具 | 做法 |
|---|---|---|
| 晶体结构、键长键角、多面体 | VESTA | 直接打开 POSCAR/CONTCAR；Utilities 中有键长、键角统计 |
| 电荷密度、差分电荷、ELF、部分电荷 | VESTA | 打开 CHGCAR/PARCHG/ELFCAR，设置等值面；2D Data Display 切片 |
| MD 轨迹、扩散路径、配位分析 | OVITO | 打开 XDATCAR 或 ASE 导出的 extxyz；修改器中有 RDF、配位数、位移等 |
| NEB 路径动画 | OVITO/VESTA | 把 00…NN 的结构合并成 extxyz（ASE）后播放 |
| 出版级能带/DOS | sumo | `sumo-bandplot`、`sumo-dosplot`、`sumo-bandstats`（带边与有效质量）、`sumo-optplot`、`sumo-phonon-bandplot`；k 路径用 `sumo-kgen`，杂化泛函用 `sumo-kgen --hybrid`；参数用 `-h` 查看 |
| 投影能带、自旋纹理 | pyprocar | `pyprocar.bandsplot(code="vasp", dirname=".", mode="parametric", ...)` |
| 快速检查收敛 | rsgrad、[check_vasp.py](../scripts/check_vasp.py) | 读取 OUTCAR 中每一步的能量、受力、磁矩 |

> [!TIP]
> 论文作图时把数据导出为纯文本（py4vasp `to_csv`、VASPKIT 的 `*.dat`、pymatgen 的数组），再用统一的 matplotlib 样式或 Origin 作图，比直接截取各工具的默认图更容易保持风格一致。

---

## 8. 自动化

### 8.1 Shell 批处理：参数扫描

```bash
#!/bin/bash
# lattice-constant scan: one directory per value, the same INCAR/KPOINTS/POTCAR everywhere
for a in 5.38 5.40 5.42 5.44 5.46 5.48; do
  d=a_$a; mkdir -p $d
  cp INCAR KPOINTS POTCAR $d/
  sed "2s/.*/$a/" POSCAR > $d/POSCAR          # line 2 of POSCAR = scaling factor
done
ls -d a_* > dirs.txt
```

收集结果：

```bash
for d in $(cat dirs.txt); do
  e=$(grep "energy(sigma->0)" $d/OUTCAR | tail -1 | awk '{print $NF}')
  v=$(grep "volume of cell" $d/OUTCAR | tail -1 | awk '{print $NF}')
  echo "$d $v $e"
done | sort -k2 -n > ev.dat
```

### 8.2 Slurm 作业数组与依赖链

```bash
N=$(wc -l < dirs.txt)
sbatch --array=1-${N}%10 array_job.slurm        # 最多同时运行 10 个
```

作业数组模板 [array_job.slurm](../examples/18-postprocess-toolkit/array_job.slurm) 用 `SLURM_ARRAY_TASK_ID` 从 `dirs.txt` 取目录，并自动跳过已正常结束的计算，重复提交是安全的。

前后依赖的任务（自洽 → 非自洽能带 → 后处理）用依赖链：

```bash
jid1=$(sbatch --parsable scf.slurm)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 band.slurm)
sbatch --dependency=afterok:$jid2 post.slurm
```

调度器的更多用法见 Common 第2章；仓库的 job_vasp_slurm.sh / job_vasp_pbs.sh 是通用模板。

### 8.3 作业健康检查

```bash
python ../scripts/check_vasp.py -h          # 仓库脚本：步数、是否收敛、最大受力、磁矩、耗时与常见报错的中文提示
grep -c "General timing" */OUTCAR           # 快速判断哪些作业正常结束
grep -l "reached required accuracy" */OUTCAR
```

### 8.4 custodian：自动检测并纠正错误

custodian 包装 VASP 的运行：每次运行后用错误处理器扫描输出，发现已知错误时按预设策略修改 INCAR/KPOINTS（例如 ZBRENT → 用 CONTCAR 续算并收紧 EDIFF；BRMIX → 依次尝试续算、Kerker 混合、改网格、关闭对称性；EDDDAV → `ALGO = All`），备份旧文件后重新运行，所有操作记录在 `custodian.json` 中。

```python
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob

handlers = [VaspErrorHandler(), UnconvergedErrorHandler()]
jobs = [VaspJob(["mpirun", "-np", "32", "vasp_std"])]       # stdout -> vasp.out
Custodian(handlers, jobs, max_errors=5).run()
```

完整脚本（含双重弛豫模式）见 [custodian_run.py](../examples/18-postprocess-toolkit/custodian_run.py)。custodian 的 `VaspErrorHandler` 源码本身就是一份很好的"VASP 报错字符串大全"。

> [!WARNING]
> 自动纠错可能悄悄改变物理设置（例如把 `ISMEAR = -5` 改成 0、关闭对称性、改变 k 网格类型）。批量计算结束后务必检查 `custodian.json`，确认同一批数据的设置仍然一致。

### 8.5 atomate2：工作流

```python
from jobflow import run_locally
from pymatgen.core import Structure
from atomate2.vasp.flows.core import RelaxBandStructureMaker

flow = RelaxBandStructureMaker().make(Structure.from_file("POSCAR"))
run_locally(flow, create_folders=True)       # 双重弛豫 -> 均匀网格与线模式能带
```

- 运行前在 `~/.atomate2.yaml` 中设置 `VASP_CMD`（或环境变量 `ATOMATE2_VASP_CMD`），并按 §4.3 配置 pymatgen 的 POTCAR 路径；
- atomate2 内置弛豫、能带、弹性、声子、缺陷、LOBSTER、MLFF 等工作流，底层用 custodian 运行 VASP；
- `run_locally` 适合在一个作业中顺序执行；成百上千个任务应使用 jobflow-remote 或 FireWorks 管理队列与数据库，见 Common 第6章。

示例脚本：[atomate2_relax_band.py](../examples/18-postprocess-toolkit/atomate2_relax_band.py)。

---

## 9. 报错排查方法

### 9.1 五步定位

```bash
# 1) 作业是否正常结束？
grep -c "General timing" OUTCAR              # 0 表示中途退出或被杀
tail -20 slurm-*.out vasp.out                # 调度器与 VASP 的最后输出

# 2) 标准输出里有没有明确的错误？
grep -n -i -E "error|very bad news|bug" vasp.out | head

# 3) OUTCAR 中的警告与建议框
grep -n "WARNING" OUTCAR | head -20
grep -n -A12 "ADVICE" OUTCAR | head -60      # 含 ADVICE 字样的提示框（例如大超胞时建议 LREAL = Auto）

# 4) 收敛过程：电子步是否到达 NELM、能量是否振荡
grep -E "^ *(DAV|RMM|CG|SDA) *:" OSZICAR | awk '{print $2}' | sort -n | tail -1   # 最大电子步数
grep " F= " OSZICAR | tail -5

# 5) 结构是否合理：最近邻距离、是否原子重叠
grep -A3 "nearest neighbor table" OUTCAR | head
```

几点经验：

- VASP 的警告框（大字号的 WARNING）不一定导致崩溃，但几乎都值得读一遍，例如 `wrap around errors must be expected`（建议 `PREC = Accurate`）、关于 NBANDS 和并行的提示；
- 同一个报错可能有多个原因，修改时**一次只改一个参数**，并在新目录中测试；
- 修改前先排除最常见的低级错误：POTCAR 顺序、INCAR 中的 Tab/CRLF、MAGMOM 个数、KPOINTS 网格类型；
- 用 `cp CONTCAR POSCAR` 续算前确认 CONTCAR 非空且格式完整（作业被杀时 CONTCAR 可能只写了一半）。

### 9.2 可以忽略的输出

| 输出 | 含义 |
|---|---|
| `WARNING: random wavefunctions but no delay for mixing, default for NELMDL` | 从随机波函数开始且未设置 NELMDL，VASP 自动处理 |
| `The number of bands has been changed from X to Y` | 为适配并行增加了 NBANDS，一般计算无影响（LOBSTER、LOPTICS、续训 MLFF 等需要固定 NBANDS 时除外） |
| 第一个离子步的电子步数很多 | 从随机波函数开始的正常现象 |

### 9.3 版本相关的已知问题

升级或排查奇怪结果时先查 VASP [Known issues](https://www.vasp.at/wiki/index.php/Known_issues)。与本教程相关的例子：6.6.0 中自旋极化 meta-GGA 受力错误（issue 94，6.6.1 修复或打 `pot.F.patch`）；`ISPIN = 2` + `KPOINTS_OPT` 时第二自旋通道的 PDOS 不写入 vaspout.h5（issue 69，6.6.1 修复）；`KPOINTS_WAN` 在高对称晶格上的插值错误（issue 105，见第17章）。

---

## 10. 仓库规范与合规

### 10.1 不能公开的内容

| 内容 | 原因 | 正确做法 |
|---|---|---|
| POTCAR（任何版本的赝势文件） | VASP 许可证保护的保密材料 | 只写"如何生成"（`cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR` 或 vaspkit 103） |
| VASP 源码、补丁后的源码 | 同上 | 只描述修改步骤；第三方补丁（VTST、VASPsol）链接到其官方仓库 |
| 官方 `makefile.include` 模板原文 | 随 VASP 分发 | 只列出需要修改的行 |
| WAVECAR、CHGCAR、vaspout.h5 等大文件 | 体积大，且可能包含可重构的赝势信息 | 写入 `.gitignore`；需要共享数据时用 NOMAD、Materials Cloud 等数据平台 |

本仓库根目录的 `.gitignore` 已排除 POTCAR、WAVECAR、CHGCAR、CHG；更完整的模板见 [.gitignore_vasp](../examples/18-postprocess-toolkit/.gitignore_vasp)。如果不小心提交了 POTCAR，仅仅在下一次提交中删除是不够的——它仍保留在 Git 历史中，需要改写历史（如 `git filter-repo`）并在远程强制更新，必要时联系平台删除缓存。

### 10.2 VASP 的 AI 工具使用规范（2026-06-16）

VASP 官方于 2026 年 6 月 16 日发布了 AI 工具使用指南，核心原则是：**AI 可以帮助你使用 VASP，但不能让 AI 摄入 VASP 软件本身**（源码、PP 数据库即 POTCAR 及相关材料）。要点：

- **允许**：用 AI 帮助编写提交脚本和后处理脚本；分析、解读 OUTCAR、h5、日志等计算输出；一般性地讨论方法和工作流；在自托管或明确配置了"不用于训练"的企业级服务中做编译、调试、代码协助；
- **禁止**：把 VASP 源文件或你的修改粘贴到公共、免费或消费级模型中"解释"或"重构"；把 PP 数据库上传到没有不训练承诺的工具；用 VASP 源码训练、微调或构建任何 AI 模型或功能等价物；把仓库访问权限或登录凭据交给 AI 工具；
- 计算结果、日志和输出数据（只要不泄露源码或 POTCAR 内容）不受限制，可以正常分析和发表；
- 实用检验：发送任何与 VASP 软件相关的内容之前，先确认该服务承诺不保留、不用于训练你的输入，并且这些设置已默认生效或已手动开启。

> [!CAUTION]
> 上面是对官方指南的概括，不构成法律意见。具体义务以你所在机构的 VASP 许可协议和官方原文为准：<https://vasp.at/info/post/>。

### 10.3 计算数据的整理建议

- 每个计算目录保留：INCAR、KPOINTS、POSCAR、CONTCAR、OSZICAR、压缩后的 OUTCAR 或 vasprun.xml（`xz -9`），以及一份说明 POTCAR 选择（`grep TITEL POTCAR` 的输出）的文本；
- 在项目根目录用表格记录每个目录的目的、VASP 版本、关键参数和结论；
- 发表时把输入文件和关键输出上传到数据平台并在论文中给出 DOI。

---

## 工具与参数选择建议

| 问题 | 建议 |
|---|---|
| VASP ≥ 6.4 且编译了 HDF5 | 优先 py4vasp（KPOINTS_OPT、部分电荷、声子、MLFF 输出最完整），pymatgen/VASPKIT 作交叉验证 |
| 没有 HDF5 或 VASP 5.4.4 | pymatgen（vasprun.xml）+ VASPKIT（DOSCAR/PROCAR/EIGENVAL） |
| 只需带隙/DOS 数字 | VASPKIT 911/11x，或 pymatgen `eigenvalue_band_properties` |
| 需要大规模入库 | pymatgen + JSON/数据库，或 atomate2 自带的数据存储 |
| 计算数 < 50 | Shell 循环 + 作业数组 + check 脚本 |
| 计算数 50–1000 | 作业数组 + custodian；或 atomate2 + jobflow-remote |
| 计算数 > 1000、多步依赖 | atomate2/FireWorks/AiiDA 等工作流引擎（Common 第6章） |
| 结果需要论文复现 | 记录所有软件版本（VASP、POTCAR 版本、pymatgen、phonopy 等），脚本放进仓库，数据上传平台 |

### 提交前检查清单

- [ ] POSCAR 元素行与 POTCAR 的 `TITEL` 顺序一致；
- [ ] INCAR 无 Tab、无 CRLF（`file INCAR`、`dos2unix INCAR`）；MAGMOM 个数正确；
- [ ] ENCUT、k 点来自收敛测试；六角/二维体系使用 Γ 中心网格；
- [ ] `ISMEAR`/`SIGMA` 与体系类型匹配，`ISMEAR = -5` 时 k 点 ≥ 4 且为 Γ 中心；
- [ ] 续算时 CHGCAR/WAVECAR 与当前设置兼容；
- [ ] 并行参数（NCORE、KPAR）是进程数的因子；GPU 版每块 GPU 一个 rank；
- [ ] 作业脚本中有 `ulimit -s unlimited`，输出重定向到 `vasp.out`；
- [ ] 大文件不在 Git 跟踪范围内。

---

## 常见问题与报错

下表列出最常见的错误；更完整的对照表（含资源与环境类问题）见 [error_table.md](../examples/18-postprocess-toolkit/error_table.md)。

| 错误信息（关键片段） | 原因 | 修复 |
|---|---|---|
| `Error reading item 'MAGMOM' from file INCAR` | MAGMOM 个数与原子数不符（非共线时每原子 3 个分量）；INCAR 有 Tab/CRLF | 用 `N*value` 写满；`vasp_ncl` 写 3×NIONS 个数；`dos2unix INCAR` |
| `Tetrahedron method fails (number of k-points < 4)` / `Routine TETIRR needs special values` | `ISMEAR = -5` 需要至少 4 个 k 点且 Γ 中心 | 分子、slab、Γ 点计算改 `ISMEAR = 0`、`SIGMA = 0.05`；金属优化不用 −5 |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角/三方体系用了偶数 Monkhorst–Pack 网格 | 改 Γ 中心网格，或调整网格与对称性相容 |
| `ZBRENT: fatal error in bracketing` | 线搜索失败：离极小点远或受力噪声大 | `cp CONTCAR POSCAR`；`EDIFF = 1E-6`；减小 POTIM；`IBRION = 1` 或 VTST FIRE（`IBRION = 3`、`POTIM = 0`、`IOPT = 7`） |
| `WARNING: Sub-Space-Matrix is not hermitian in DAV` / `Error EDDDAV: Call to ZHEGV failed` | Davidson 失败：原子重叠、步长大、`LREAL` 误差 | 检查结构；`ALGO = All`（6.5.1+ 可加 `ISEARCH = 1`）；`LREAL = .FALSE.`；从头开始 |
| `BRMIX: very serious problems the old and the new charge density differ` | 读入了不匹配的 CHGCAR，对称性判定不一致 | 删除 CHGCAR/WAVECAR 重算；`ISYM = 0` 或调整 `SYMPREC` |
| `internal error in subroutine PRICEL` / `Found some non-integer element in rotation matrix` / `SGRCON` | 结构对称性有噪声 | spglib/vaspkit 整理结构；调整 `SYMPREC`；必要时 `ISYM = 0` |
| `WARNING: The number of bands has been changed from X to Y` | VASP 为并行自动加大 NBANDS | 需要固定 NBANDS 时设为能带组数的整数倍 |
| `ERROR RSPHER` / `REAL_OPTLAY: internal error` | 实空间投影失败 | `LREAL = .FALSE.` |
| `LAPACK: Routine ZPOTRF failed` | 原子过近或步长过大；小体系 NCORE 过大 | 检查结构；减小 POTIM；`ISYM = 0`；减小 NCORE |
| `ELF: KPAR>1 not implemented` | ELF 不支持 k 点并行 | `KPAR = 1` |
| `PEAD routines do not work for NCORE` | LEPSILON/LCALCEPS/DFPT 不支持能带并行 | 删除 NCORE/NPAR |
| `One of the lattice vectors is very long (>50 A), but AMIN is rather large` | 长真空方向 | `AMIN = 0.01` |
| `Spin polarized Harris functional dynamics is a good joke` | `ICHARG ≥ 10` + `ISPIN = 2` + `NSW > 0` | 非自洽计算 `NSW = 0` |
| `forrtl: severe (174): SIGSEGV` | 栈空间或内存不足 | `ulimit -s unlimited`；减少每节点进程数 |
| entropy T*S 每原子超过约 1 meV | SIGMA 过大 | 减小 SIGMA；比较能量用 `energy(sigma->0)` |
| ISIF=3 优化后晶格偏小/偏大 | Pulay 应力 | `ENCUT ≥ 1.3×ENMAX`，CONTCAR 再优化直至体积不变 |
| 6.6.0 自旋极化 SCAN/r2SCAN 受力错误 | Known issue 94 | 升级到 6.6.1 或打 `pot.F.patch` |
| py4vasp 报找不到数据 / 没有 vaspout.h5 | VASP 未开启 HDF5，或该量未写入（如部分电荷需 `LPARDH5`） | 重新编译加 `-DVASP_HDF5`；按需设置对应标签；或改用 pymatgen |
| pymatgen 读 vasprun.xml 报解析错误 | 作业中断，XML 不完整 | 确认作业正常结束；必要时用 OUTCAR/OSZICAR 提取数据 |
| ASE 报 `No pseudopotential for ...` | 未设置 `VASP_PP_PATH` 或目录结构不对 | `VASP_PP_PATH` 指向包含 `potpaw_PBE/` 的目录 |

---

## 练习

1. **三种工具对带隙**：对第5章的 Si 能带/DOS 计算，分别用 VASPKIT 911、pymatgen `eigenvalue_band_properties` 和 py4vasp `calc.bandgap` 读带隙，解释可能的差别（k 点集合、费米能定义）。
2. **作业数组**：为 GaAs 生成 7 个晶格常数目录，用 `array_job.slurm` 提交，用 §8.1 的命令收集 E–V 数据，再用 ASE 的 `EquationOfState` 拟合。
3. **custodian 实验**：故意在 INCAR 中设 `ISMEAR = -5` 并只用 Γ 点 KPOINTS，用 `custodian_run.py` 运行，查看 `custodian.json` 记录了什么修改。
4. **报错复现**：对六角 MoS2 使用 4×4×1 Monkhorst–Pack 网格，观察 VASP 的输出；再改为 Γ 中心网格对比。
5. **仓库体检**：把 `.gitignore_vasp` 加到自己的项目，用 `git status --ignored` 检查哪些大文件被正确忽略。

---

## 参考资料

**官方文档**

- VASP Wiki 总入口：<https://www.vasp.at/wiki/index.php/The_VASP_Manual>；[Known issues](https://www.vasp.at/wiki/index.php/Known_issues)；[Changelog](https://www.vasp.at/wiki/index.php/Changelog)
- VASP News（含 2026-06-16 AI 工具使用指南）：<https://vasp.at/info/post/>
- py4vasp 文档：<https://vasp.at/py4vasp/latest/>
- VASPKIT 文档（1.5.x）：<https://vaspkit.com/>；下载：<https://sourceforge.net/projects/vaspkit/>
- pymatgen `pymatgen.io.vasp`：<https://pymatgen.org/pymatgen.io.vasp.html>
- ASE VASP 计算器：<https://docs.ase-lib.org/ase/calculators/vasp.html>
- sumo 文档：<https://smtg-bham.github.io/sumo/>
- custodian 文档：<https://materialsproject.github.io/custodian/>
- atomate2 文档：<https://materialsproject.github.io/atomate2/>

**相关 GitHub 仓库**（只链接，不转载代码）

- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：VASP 官方 Python 后处理库
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：pymatgen（核心 I/O 位于 pymatgen-core）
- [materialsproject/custodian](https://github.com/materialsproject/custodian)：VaspErrorHandler 收录大量真实报错字符串与修复策略
- [materialsproject/atomate2](https://github.com/materialsproject/atomate2)：高通量 VASP 工作流
- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo)：出版级能带/DOS 作图与 k 路径
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar)：投影能带、自旋纹理、费米面
- [orest-d/p4vasp](https://github.com/orest-d/p4vasp)：老牌 vasprun.xml 可视化 GUI（已基本停更）
- [Ionizing/rsgrad](https://github.com/Ionizing/rsgrad)：快速查看 OUTCAR 收敛过程、受力与磁矩
- [PytLab/VASPy](https://github.com/PytLab/VASPy)：轻量级 VASP 文件读写库
- [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script)：中文社区常用脚本集合（含二进制与官方 PDF，只能链接不能转载）
- [vaspkit/VASPKIT_NatureProtocols](https://github.com/vaspkit/VASPKIT_NatureProtocols)：VASPKIT Nature Protocols 论文配套示例
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book)：经典的 VASP + ASE 教程书
- [Bit-Part-Young/scripts](https://github.com/Bit-Part-Young/scripts)：中文注释的 ASE/pymatgen/phonopy/vaspkit 脚本集
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)；本仓库作者的 NEP/GPUMD/数据集脚本见 [Scripts/README.md](../../Scripts/README.md)

**关键文献**

- V. Wang et al., *Comput. Phys. Commun.* 267, 108033 (2021)（VASPKIT）. DOI: [10.1016/j.cpc.2021.108033](https://doi.org/10.1016/j.cpc.2021.108033)
- S. P. Ong et al., *Comput. Mater. Sci.* 68, 314 (2013)（pymatgen，也是 custodian 的推荐引用）. DOI: [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)
- A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017)（ASE）. DOI: [10.1088/1361-648X/aa680e](https://doi.org/10.1088/1361-648X/aa680e)
- A. M. Ganose, A. J. Jackson, D. O. Scanlon, *J. Open Source Softw.* 3, 717 (2018)（sumo）. DOI: [10.21105/joss.00717](https://doi.org/10.21105/joss.00717)
- U. Herath et al., *Comput. Phys. Commun.* 251, 107080 (2020)（PyProcar）. DOI: [10.1016/j.cpc.2019.107080](https://doi.org/10.1016/j.cpc.2019.107080)
- A. M. Ganose et al., *Digital Discovery* (2025)（atomate2）. DOI: [10.1039/D5DD00019J](https://doi.org/10.1039/D5DD00019J)
- K. Momma, F. Izumi, *J. Appl. Cryst.* 44, 1272 (2011)（VESTA）. DOI: [10.1107/S0021889811038970](https://doi.org/10.1107/S0021889811038970)
- A. Stukowski, *Modelling Simul. Mater. Sci. Eng.* 18, 015012 (2010)（OVITO）. DOI: [10.1088/0965-0393/18/1/015012](https://doi.org/10.1088/0965-0393/18/1/015012)

---

[← 上一章](17-advanced-electronic.md) | [返回目录](../README.md)
