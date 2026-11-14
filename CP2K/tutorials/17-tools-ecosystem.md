# 第17章 自动化与工具生态：ASE、cp2k-input-tools、Multiwfn、cp2kdata（Automation & Tool Ecosystem: ASE, cp2k-input-tools, Multiwfn, cp2kdata）

本章介绍 CP2K 周边常用的自动化工具，以及它们各自适用的场合。内容包括：ASE CP2K 计算器（shell 模式、`ASE_CP2K_COMMAND`、`inp` 模板，以及旧版 ASE 对 `cp2k_shell` 的断言）；cp2k-input-tools（`cp2klint`、`cp2kgen` 参数扫描、`fromcp2k`/`tocp2k` 的 YAML 输入）；用 Multiwfn 一键生成 CP2K 输入，以及生成后按 2026.2 语法要做的检查；用 cp2kdata、cp2k-output-tools 解析输出，用 dpdata 转成 DeePMD 训练集；aiida-cp2k 2.2.1 的高通量工作流；Slurm/PBS 作业模板。全章都围绕一个原则：工具负责生成和搬运，**`cp2k.psmp -c` 和收敛测试负责判断对错**。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 CP2K shell 模式的工作方式，按 ASE 版本正确设置 `ASE_CP2K_COMMAND`，必要时用 wrapper 绕过旧版 ASE 的断言；
- 用 ASE 的 `CP2K` 计算器做分子优化和晶体 E–V 曲线，把 k 点、OT、UZH 基组等复杂设置放进 `inp` 模板；
- 用 `cp2klint`、`cp2kgen`、`fromcp2k`/`tocp2k`、`cp2kget` 检查、批量生成、转换和读取输入，并知道它们内置的输入结构版本落后于 2026.2；
- 用 Multiwfn 在几步之内生成 CP2K 输入，并按当前版本检查和修改生成的输入；
- 用 cp2kdata/cp2k-output-tools 解析输出，用 dpdata 生成 DeePMD 训练集；
- 了解 AiiDA（aiida-cp2k）适合什么场合，会写 Slurm/PBS 作业脚本，包括运行 ASE 驱动脚本的作业。

| 任务 | 推荐工具 | 本章位置 | 示例文件 |
|---|---|---|---|
| 检查语法（以所用版本为准） | `cp2k.psmp -c` / `-e` | §3 | — |
| Python 里循环、优化、EOS、与 MLIP 混用 | ASE `CP2K` 计算器 | §4 | `ase_h2o_opt.py`、`ase_si_eos.py` |
| 批量改一个或几个参数（CUTOFF、k 点、晶格常数） | `cp2kgen`，或预处理器 `-E` | §3、§5 | `cp2kgen_cutoff_scan.sh` |
| 用 YAML/JSON 管理输入 | `tocp2k` / `fromcp2k` | §5 | `inp_from_yaml.yaml` |
| 从 cif/POSCAR/gjf 快速得到“能跑”的输入 | Multiwfn（`cp2k` 功能）、本仓库 [`cp2k_inpgen.py`](../scripts/cp2k_inpgen.py) | §6 | `Si.POSCAR`、`inp_migration_check.py` |
| 解析能量、力、优化步、cube、pdos | cp2kdata、cp2k-output-tools、本仓库 [`cp2k_outparse.py`](../scripts/cp2k_outparse.py) | §7 | `cp2kdata_to_dpdata.py` |
| 生成机器学习势训练集 | cp2kdata + dpdata | §7.3 | `cp2kdata_to_dpdata.py` |
| 高通量、自动重启、可追溯 | AiiDA + aiida-cp2k | §8 | `aiida_h2o_energy.py` |
| 提交作业 | Slurm / PBS 模板 | §9 | `job_slurm.sh`、`job_pbs.sh` |

示例目录：[../examples/17-automation-tools/README.md](../examples/17-automation-tools/README.md)。ASE 的通用用法（Atoms、优化器、EOS、NEB）见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)，本章只讲与 CP2K 相关的部分。

**本章涉及的版本**（2026-11 查询 PyPI/GitHub 的结果）：

| 工具 | 版本 | 发布时间 | Python | 备注 |
|---|---|---|---|---|
| CP2K | 2026.2 | 2026-07-15 | — | `-s/--shell`、`-c`、`-e`、`-E`、`--xml` 均可用；2026.2 删除了 Cython Python 绑定 |
| ASE | 3.29.0 | 2026-06-21 | ≥3.10 | 默认命令改为 `cp2k.psmp -s`；3.28 及以前默认是 `cp2k_shell`，并断言命令中含该字符串 |
| cp2k-input-tools | 0.9.3 | 2026-10-26 | ≥3.10 | 内置的输入结构为 **CP2K 2024.1** |
| cp2k-output-tools | 0.5.0（PyPI） | 2021-07 | ≥3.6 | GitHub 仓库仍在维护，但 PyPI 上的版本较旧 |
| cp2kdata | 0.7.3（PyPI） | 2025-05 | ≥3.8 | 带 dpdata 插件（`cp2kdata/e_f`、`cp2kdata/md`）；PyPI 版读不到 CP2K ≥ 2025.1 的能量行，需用 GitHub 主分支，见 §7.1 |
| dpdata | 1.1.0 | 2026-08-19 | ≥3.10 | `deepmd/npy` 等格式 |
| aiida-cp2k | 2.2.1 | 2026-11-03 | ≥3.10 | 需要 aiida-core ≥2.7、<3 |
| Multiwfn | 官网最新版 | 持续更新 | — | 生成 CP2K 输入的菜单会随版本调整 |

> [!IMPORTANT]
> 每个工具都自带一份“它所理解的 CP2K”：ASE 只知道自己会写哪些关键字，cp2k-input-tools 内置的是 2024.1 的输入结构，Multiwfn 跟随作者的更新节奏。只有你所用的 `cp2k.psmp` 才知道 2026.2 的真实语法。所以凡是工具生成的输入，都要先过一遍 `cp2k.psmp -c`。

---

## 2. 原理速览

### 2.1 CP2K 输入是一棵树

所有工具都把 CP2K 输入看成一棵树：**段（section）**可以嵌套，可以带一个**段参数**（如 `&KIND O` 中的 `O`、`&OT ON` 中的 `ON`），有的段可以重复（`&KIND`、`&FORCE_EVAL`）；**关键字（keyword）**是叶子。树中的位置可以用路径表示：

```text
force_eval/dft/mgrid/cutoff          -> &FORCE_EVAL / &DFT / &MGRID / CUTOFF
force_eval/subsys/kind/0/basis_set   -> 第 1 个 &KIND 里的 BASIS_SET（重复段用下标）
```

cp2k-input-tools 的 `cp2kgen`/`cp2kget`、aiida-cp2k 的参数字典、ASE 内部的 `InputSection`，用的都是这种表示。区别只在于它们用哪一版的“合法关键字表”来校验：CP2K 可以用 `cp2k.psmp --xml` 导出当前版本的完整输入结构（`cp2k_input.xml`）。

### 2.2 shell 模式：ASE 如何与 CP2K 对话

`cp2k.psmp -s`（或 `--shell`）启动一个常驻的 CP2K 进程，从 stdin 读命令、向 stdout 写结果。ASE 的 `CP2K` 计算器依次发送下面这些命令（命令名取自 ASE 3.29 源码）：

| 阶段 | shell 命令 | 作用 |
|---|---|---|
| 启动 | 等待 `* READY`，再发 `VERSION`、`HARSH` | 确认 shell 版本；HARSH 表示遇到任何错误立即停止 |
| 建立力场环境 | `WRITE_FILE`、`LOAD <inp> <out>` | 写出 ASE 生成的输入并加载，得到 force_env 编号 |
| 每个几何结构 | `UNITS_EV_A`、`SET_CELL`、`SET_POS`（或 `SET_POS_FILE`） | 以 eV/Å 为单位传入晶胞和坐标 |
| 计算 | `EVAL_EF`、`GET_E`、`GET_F`、`GET_STRESS` | 计算能量和力，并取回结果 |
| 结束 | `DESTROY`、`EXIT` | 释放 force_env，退出进程 |

由此可以得出几点：

1. CP2K 进程常驻，波函数可以在相邻几何结构之间外推，因此 ASE 优化每一步的 SCF 比每次重新启动快得多；
2. 通信走 stdin/stdout，**stdout 上最好不要有别的输出**（例如 wrapper 里的 `echo`、`module load` 打印的信息）。ASE ≤ 3.28 启动后要求第一行就是 `* READY`，多出任何一行都会触发 `AssertionError`；ASE 3.29 会跳过 `* READY` 之前最多 1000 行的启动信息（源码注释举的例子是 SIRIUS 构建打印的 banner），但 `* READY` 之后混入的输出仍会打乱对话；
3. MPI 运行时只有 0 号进程读 stdin，所以启动器必须把 stdin 交给 0 号进程（`mpirun`/`mpiexec` 默认如此）。

### 2.3 单位

ASE 内部使用 eV、Å、eV/Å，CP2K 内部使用 Hartree、bohr。ASE 的 `cutoff` 参数以 **eV** 为单位，会写成 `CUTOFF [eV] ...`。自己解析 CP2K 输出时（cp2kdata 的 `Cp2kOutput` 返回原生单位）要换算：

$$
1\,E_\mathrm h = 27.211386\ \mathrm{eV},\qquad 1\,\mathrm{Ry} = 13.605693\ \mathrm{eV},\qquad 1\,E_\mathrm h/a_0 = 51.422067\ \mathrm{eV/\AA}
$$

### 2.4 状态方程拟合

§4.6 用 ASE 的 `EquationOfState` 拟合三阶 Birch–Murnaghan 方程：

$$
E(V) = E_0 + \frac{9V_0B_0}{16}\left\{\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^3 B_0' + \left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^2\left[6-4\left(\frac{V_0}{V}\right)^{2/3}\right]\right\}
$$

ASE 给出的 $B_0$ 单位是 eV/Å³，乘以 160.2177 换算为 GPa。

---

## 3. CP2K 自带的命令行工具

先介绍 CP2K 自带的命令行选项，它们不需要额外安装（以下均为 2026.2 中 `cp2k.psmp -h` 列出的选项）：

| 命令 | 作用 | 典型用途 |
|---|---|---|
| `cp2k.psmp -c x.inp` | 语法检查，不运行 | **任何工具生成的输入都先跑一次** |
| `cp2k.psmp -e x.inp > full.inp` | 输出预处理后的输入，并把默认值全部显式写出 | 查默认值，确认 `-E` 变量是否生效 |
| `cp2k.psmp -E CUT=500 -i x.inp -o x.out` | 设置预处理器变量的初值（`--set`） | 不借助外部工具做参数扫描 |
| `cp2k.psmp --xml` | 在当前目录写出 `cp2k_input.xml` | 给 cp2k-input-tools 等工具提供当前版本的输入结构 |
| `cp2k.psmp -v` | 版本号和 cp2kflags（编译了哪些库） | 写在作业日志里，便于复现 |
| `cp2k.psmp -s` | 交互式 shell 模式 | ASE 计算器、i-PI 风格的外部驱动 |

最简单的参数扫描不需要任何 Python 包。模板里写 `${CUT-400}`（带默认值，见[第2章 §6.2](02-input-structure.md)），然后循环调用：

```bash
for c in 300 400 500 600; do
  mkdir -p c$c && cp template.inp c$c/ && cd c$c
  cp2k.psmp -E CUT=$c -c template.inp 2>&1 | grep -q "SUCCESS, the input could be parsed correctly" \
    || { echo "syntax error at $c"; exit 1; }
  cp2k.psmp -E CUT=$c -i template.inp -o run.out
  cd ..
done
grep "ENERGY| Total FORCE_EVAL" c*/run.out
```

`-c` 通过时会打印 `SUCCESS, the input could be parsed correctly.`，失败时打印 `ERROR,   the input could *NOT* be parsed correctly.`（2026.2 源码 `src/start/cp2k.F`）。脚本里用 grep 检查这一行，比依赖退出码更稳妥，本章的作业模板也都这样写。

> [!TIP]
> CP2K 也会把 shell 模式绑定到程序名：如果可执行文件（或指向它的符号链接）的名字以 `cp2k_shell` 开头，CP2K 启动后会自动进入 shell 模式（见 2026.2 源码 `src/start/cp2k.F`）。§4.2 的第二种绕过办法就利用了这一点。

---

## 4. ASE CP2K 计算器

### 4.1 安装与启动命令

```bash
python -m pip install "ase>=3.23"          # Python >= 3.10
python -c "import ase; print(ase.__version__)"
export CP2K_DATA_DIR=/path/to/cp2k/data    # BASIS_MOLOPT_UZH / POTENTIAL_UZH 所在目录
```

ASE 按下面的顺序决定用什么命令启动 CP2K：构造函数参数 `command=` → 类变量 `CP2K.command` → ASE 配置文件 `[cp2k]` 段的 `cp2k_shell` 项 → 环境变量 `ASE_CP2K_COMMAND` → 内置默认值。两个大版本的差别如下：

| ASE 版本 | 内置默认命令 | 对命令的要求 | 在 CMake 构建的 CP2K（只有 `cp2k.psmp`/`cp2k.ssmp`）上 |
|---|---|---|---|
| ≤ 3.28 | `cp2k_shell` | 源码里有 `assert 'cp2k_shell' in command` | `/bin/sh: cp2k_shell: not found`；直接写 `cp2k.psmp -s` 会触发 `AssertionError` |
| ≥ 3.29 | `cp2k.psmp -s` | 无断言 | `export ASE_CP2K_COMMAND="cp2k.psmp -s"` 或 `"mpirun -np 4 cp2k.psmp -s"` |

（对照 ASE 3.23/3.27/3.28/3.29 的 `ase/calculators/cp2k.py` 源码核实。）

### 4.2 旧版 ASE：wrapper 或符号链接

旧集群上常见 ASE ≤ 3.28 和 CMake 构建的 CP2K 2026.x 共存，有两种绕过办法：

**办法一：文件名含 `cp2k_shell` 的 wrapper 脚本**（[ase_cp2k_shell_wrapper.sh](../examples/17-automation-tools/ase_cp2k_shell_wrapper.sh)）：

```bash
#!/bin/bash
# no echo before exec: ASE reads stdout
NP=${CP2K_SHELL_NP:-1}
if [ "$NP" -gt 1 ]; then
  exec mpirun -np "$NP" cp2k.psmp --shell
else
  exec cp2k.psmp --shell
fi
```

```bash
chmod +x ase_cp2k_shell_wrapper.sh
export ASE_CP2K_COMMAND="$PWD/ase_cp2k_shell_wrapper.sh"   # 路径里含 cp2k_shell，满足断言
export CP2K_SHELL_NP=4 OMP_NUM_THREADS=1
```

**办法二：符号链接**。名字以 `cp2k_shell` 开头的链接会让 CP2K 自动进入 shell 模式：

```bash
mkdir -p ~/bin && ln -s "$(command -v cp2k.psmp)" ~/bin/cp2k_shell.psmp
export ASE_CP2K_COMMAND="mpirun -np 4 $HOME/bin/cp2k_shell.psmp"
```

此外还可以直接升级到 ASE ≥ 3.29，或者使用不依赖 shell 模式的第三方实现 [ase-cp2k-noshell](https://github.com/LouieSlocombe/ase-cp2k-noshell)。

### 4.3 计算器参数

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `label` | `"h2o_ase"` | 写成 `PROJECT`，输入/输出文件名为 `<label>.inp/.out` |
| `command` | `None` | 为 `None` 时按 §4.1 的顺序查找 |
| `inp` | 输入模板字符串 | **计算器参数表达不了的设置都放在这里**（OT、k 点、`REL_CUTOFF`、UZH 的 `&KIND`、D3、+U……） |
| `xc` | `"PBE"`；混合泛函用 `None` | 默认 `LDA`。普通泛函名会写成 `&XC_FUNCTIONAL / &PBE`，以 `XC_` 开头的名字按 libxc 处理；PBE0/HSE06 等设为 `None`，在模板里写完整的 `&XC` |
| `cutoff` | `500 * Rydberg`（`from ase.units import Rydberg`） | **单位是 eV**，默认 400 Ry；写成 `CUTOFF [eV]` |
| `basis_set` / `basis_set_file` | `None` / `"BASIS_MOLOPT_UZH"` | 默认值 `DZVP-MOLOPT-SR-GTH` / `BASIS_MOLOPT`。UZH 基组名里含 -qN，各元素不同，因此设为 `None`，改在模板的 `&KIND` 中指定 |
| `pseudo_potential` / `potential_file` | `None` / `"POTENTIAL_UZH"` | 默认 `auto`：LDA、PADE、BP、BLYP、PBE 会自动推断为 `GTH-<xc>`，其他泛函退回 `GTH-PBE` 并发出警告 |
| `max_scf` | 30（OT 内层）～100（对角化） | 同时写入 `&SCF` 和 `&LS_SCF` 的 `MAX_SCF`；外层循环在模板里写 `&OUTER_SCF` |
| `uks` | 开壳层设为 `True` | 写成 `UNRESTRICTED_KOHN_SHAM ON` |
| `charge` | 整数 | 非零时写成 `CHARGE` |
| `multiplicity` | 例如 O2 取 3 | 写成 `MULTIPLICITY`；默认 `None`，即不写 |
| `stress_tensor` | 孤立分子设 `False` | 默认 `True`，会写 `STRESS_TENSOR ANALYTICAL` 并打印应力 |
| `poisson_solver` | `"auto"` | 当 `atoms.pbc` 全为 False 时，自动加 `&POISSON PERIODIC NONE` 和 `PSOLVER MT` |
| `print_level` | `"MEDIUM"` | 默认 `LOW`。需要用 cp2kdata 解析时必须设为 MEDIUM |
| `set_pos_file` | 原子数 > 100 的 MPI 运行设为 `True` | 通过文件传坐标，避免某些 MPI 库在 stdin 传大量数据时卡住；需要 CP2K ≥ 2024.2 |
| `debug` | 排错时设为 `True` | 打印与 shell 之间的全部通信 |

ASE 会**自动写入**以下内容：`GLOBAL/PROJECT`、`GLOBAL/PRINT_LEVEL`、`FORCE_EVAL/METHOD`、`STRESS_TENSOR`、`BASIS_SET_FILE_NAME`、`POTENTIAL_FILE_NAME`、`MGRID/CUTOFF`、`SCF/MAX_SCF`、`XC`、`UKS`、`MULTIPLICITY`、`CHARGE`、`POISSON`（孤立体系）、`SUBSYS/COORD`、`SUBSYS/CELL`（A/B/C 与 PERIODIC），以及每个元素的 `&KIND`。模板里如果再写同名关键字，`add_keyword` 会抛出 `Exception: Keyword CUTOFF already present in section MGRID` 这样的错误。

> [!WARNING]
> 模板中不要使用 `@SET`/`@INCLUDE` 等预处理指令。ASE 自己解析模板再重新写出，并不认识预处理器。模板里也不要放 `&MOTION`：几何优化和 MD 由 ASE 的优化器或积分器负责，CP2K 只做单点的能量和力。

### 4.4 实战一：H2O 结构优化（ase_h2o_opt.py）

[ase_h2o_opt.py](../examples/17-automation-tools/ase_h2o_opt.py) 的核心部分如下（完整脚本带参数解析和 ASE 版本检查）：

```python
from ase.build import molecule
from ase.calculators.cp2k import CP2K
from ase.optimize import BFGS
from ase.units import Rydberg

INP = """
&FORCE_EVAL
  &DFT
    &MGRID
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      EPS_SCF 1.0E-7
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
  &END DFT
  &SUBSYS
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""

atoms = molecule("H2O")
atoms.center(vacuum=5.0)
atoms.pbc = False                      # -> PERIODIC NONE + MT Poisson solver
calc = CP2K(label="h2o_ase", inp=INP, xc="PBE", cutoff=500 * Rydberg, max_scf=30,
            basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH",
            basis_set=None, pseudo_potential=None, stress_tensor=False,
            print_level="MEDIUM")
with calc:                             # closes the CP2K shell at the end
    atoms.calc = calc
    BFGS(atoms, trajectory="h2o_opt.traj", logfile="h2o_opt.log").run(fmax=0.01)
```

```bash
cd CP2K/examples/17-automation-tools
export ASE_CP2K_COMMAND="cp2k.psmp -s"
python ase_h2o_opt.py --fmax 0.01 --cutoff 500
cp2k.psmp -c h2o_ase.inp                       # 复核 ASE 生成的输入
grep -c "SCF run NOT converged" h2o_ase.out    # 应为 0
grep "ENERGY| Total FORCE_EVAL" h2o_ase.out | tail -3
tail -3 h2o_opt.log                            # BFGS: step, energy, fmax
```

成功的标准：`h2o_opt.log` 最后一行的 fmax < 0.01 eV/Å，输出中没有 `SCF run NOT converged`，打印的键长和键角与 PBE 的典型值（O–H ≈ 0.97 Å，∠HOH ≈ 104°）在同一量级。注意 `fmax=0.01 eV/Å` 约为 $1.9\times10^{-4}\,E_\mathrm h/a_0$，比 CP2K `GEO_OPT` 的默认 `MAX_FORCE`（$4.5\times10^{-4}$）更严格。

> [!NOTE]
> 与在 CP2K 里直接做 `RUN_TYPE GEO_OPT`（[第7章](07-geo-cell-opt.md)）相比，ASE 驱动的好处是可以使用 ASE 的优化器（BFGS、FIRE、Sella 鞍点搜索）和约束（`FixAtoms`、`FixBondLength`），也可以与 MLIP 混合使用；代价是每一步都要经过 Python 通信，并且失去 CP2K 自带的重启文件（`-1.restart`）。纯 DFT 的常规优化仍推荐用 CP2K 自身的 `GEO_OPT`。

### 4.5 开壳层、带电体系和杂化泛函

```python
# O2 triplet
calc = CP2K(label="o2", inp=INP_O, xc="PBE", uks=True, multiplicity=3, charge=0,
            basis_set=None, pseudo_potential=None, stress_tensor=False,
            basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH")

# hybrid functional: xc=None, full &XC (and ADMM, &HF) in the template
calc = CP2K(label="h2o_pbe0", inp=INP_PBE0, xc=None, ...)
```

杂化泛函的 `&XC`、`&HF`、`&AUXILIARY_DENSITY_MATRIX_METHOD` 和截断库仑的写法见[第9章](09-hybrid-admm.md)。模板里写了 `&XC` 时，`xc` 必须设为 `None`，否则 ASE 会在 `&XC_FUNCTIONAL` 中再加一个 `&PBE` 子段，结果变成两个泛函相加。

### 4.6 实战二：体相 Si 的 E–V 曲线（ase_si_eos.py）

周期体系的要点：k 点、对角化、混合都写在模板里；`atoms.pbc=True` 时 ASE 不会加 MT 求解器。2026.2 及以前的版本中，k 点只能和对角化一起用（OT 不支持 k 点）。模板如下（`{k}` 由脚本替换）：

```text
&FORCE_EVAL
  &DFT
    &MGRID
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK {k} {k} {k}
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-8
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
    &END SCF
  &END DFT
  &SUBSYS
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `SCHEME MONKHORST-PACK` | 6 6 6 起，需做收敛测试 | 2 原子原胞；k 点收敛方法见[第6章](06-periodic-surface-kpoints.md)，也可用 [`conv_kpoints.py`](../scripts/conv_kpoints.py) |
| `&DIAGONALIZATION` | `ON`，`ALGORITHM STANDARD` | 2026.2 及以前 k 点必须用对角化 |
| `&MIXING` | Broyden，`ALPHA` 0.4 | 对角化时的密度混合 |
| `EPS_SCF` | 1e-8 | E–V 曲线上相邻点的能量差只有 meV 量级，SCF 必须收得更紧 |
| `stress_tensor`（计算器参数） | `False` | 拟合 E–V 只需要能量 |

```bash
python ase_si_eos.py --kpts 6 --cutoff 500 --npoints 7 --range 0.03
cat si_eos.dat
grep -l "SCF run NOT converged" si_eos_*.out     # 应无输出
grep "BRILLOUIN|" si_eos_00.out | head -4        # 确认 k 点生效：K-point scheme Monkhorst-Pack, K-Point grid 6 6 6
```

脚本对每个体积**新建一个计算器**。原因是：如果只改晶胞，ASE 会向正在运行的 shell 发送 `SET_CELL`，沿用已有的 force_env。为每个体积单独建立 force_env 才是通常意义上“固定 CUTOFF”的 E–V 曲线，每个点也彼此独立。输出示例：

```text
a = 5.2681 A  V/atom = 18.2771 A^3  E/atom = -107.xxxxxx eV
...
V0 = 20.4xx A^3/atom, a0 = 5.4xx A, B = 8x.x GPa
```

（示意输出，数值仅供参考）

PBE 下 Si 的晶格常数文献值约为 5.47 Å，略高于实验值 5.43 Å（Haas, Tran, Blaha, *Phys. Rev. B* **79**, 085104 (2009)）。若与此偏差超过 0.02 Å，先检查 k 点、CUTOFF 和基组的收敛性。CP2K 原生的晶胞优化见[第7章](07-geo-cell-opt.md)。

### 4.7 并行、作业与读写 CP2K 文件

- **并行**：把 MPI 启动器写进命令，例如 `ASE_CP2K_COMMAND="mpirun -np 32 cp2k.psmp -s"`。Python 进程本身是串行的，只占一个核。在 Slurm 下运行，见 §9 的 `job_slurm.sh --ase`。
- **大体系**：原子数超过 100 且使用 psmp 时，ASE 会警告 stdin 传坐标可能卡住，此时设置 `set_pos_file=True`（CP2K ≥ 2024.2）。
- **读写文件**：ASE 能直接读取 CP2K 的重启文件和 DCD 轨迹。2026.1 起 CP2K 还能输出 extended XYZ 轨迹，可以用 `ase.io.read(..., format="extxyz")` 读取。

```python
from ase.io import read
final = read("geo-1.restart", format="cp2k-restart")           # from a CP2K GEO_OPT run
traj = read("md-pos-1.dcd", index=":", format="cp2k-dcd")      # needs ref_atoms for symbols
```

---

## 5. cp2k-input-tools：检查、扫描与 YAML 输入

### 5.1 安装与命令一览

```bash
python -m pip install "cp2k-input-tools[yaml]"     # [yaml] 提供 ruamel.yaml；[lsp] 提供语言服务器
```

| 命令 | 作用 | 输入结构来源 |
|---|---|---|
| `cp2klint x.inp` | 语法检查，给出行号和出错位置 | 内置（2024.1） |
| `cp2kgen x.inp "path=[v1,v2]" ...` | 按表达式批量生成输入；`--zip` 表示多个表达式逐项配对，默认做笛卡尔积 | 内置 |
| `fromcp2k x.inp` | 转成 JSON（默认）、`-f yaml`、`-f aiida-cp2k-calc`（生成 aiida-cp2k 运行脚本模板） | **可用 `--xml` 或环境变量 `FROMCP2K_XML_DEFINITION` 指定** |
| `tocp2k [--yaml] x.json` | JSON/YAML 转回 CP2K 输入 | 内置 |
| `cp2kget x.inp "force_eval/dft/mgrid/cutoff"` | 按路径取值，常用于从 `-1.restart` 读出晶胞 | 内置 |
| `cp2k-datafile-lint basis FILE` | 检查和美化基组/赝势数据文件（格式参数为 `basis` 或 `pseudo`） | — |
| `cp2k-language-server` | 编辑器（VS Code、Vim 等）的补全和诊断后端 | 内置 |

> [!CAUTION]
> 0.9.3 内置的 `cp2k_input.xml` 是 **CP2K 2024.1** 的输入结构（文件头写着 `CP2K version 2024.1`）。因此 2024.1 之后新增的写法（例如 2025.1 起 `&XTB` 的 `GFN_TYPE`、2026.2 在 `&DFT/&PRINT/&DOS` 里内嵌的 `&PDOS`/`&CURVE`/`&BROADEN`）会被 `cp2klint`、`cp2kgen`、`tocp2k` 报错；而 2026.2 已经删除的 `&CELL_OPT TYPE`、`&PRINT &PDOS`、REFTRAJ 的 `EVAL_ENERGY_FORCES` 反而能通过。REFTRAJ 的新关键字 `EVAL` 在 2024.1 中已经存在，所以两种写法它都接受。`cp2klint` 通过了只说明它符合 2024.1 的语法，**不能代替** `cp2k.psmp -c`。

`fromcp2k` 支持用 `--xml` 指定外部的输入结构，可以借此按你所用的 CP2K 版本检查输入（是否完全兼容 2026.2 导出的 XML，请自行验证）：

```bash
cp2k.psmp --xml                                   # writes cp2k_input.xml (current version)
fromcp2k --xml cp2k_input.xml my.inp > /dev/null && echo "parsed with the 2026.2 schema"
export FROMCP2K_XML_DEFINITION=$PWD/cp2k_input.xml   # 之后 fromcp2k 默认使用这份 XML
```

### 5.2 cp2klint

```bash
cp2klint h2o_cutoff_template.inp
# All done! Happy calculating!                    <- 通过
```

出错时打印 `Syntax error: <原因>, in <文件>:`，再给出行号和 `^` 标记的位置；如果是关键字重复，还会显示“previous definition”所在的行。它的优点是不需要安装 CP2K，在笔记本上就能用（例如在 VS Code 里配合 `cp2k-language-server` 实时提示）。

### 5.3 YAML 驱动输入：tocp2k / fromcp2k

“简化格式”的规则（与 `fromcp2k` 默认输出一致）：

| 规则 | 例子 |
|---|---|
| 段和关键字用小写的嵌套映射 | `force_eval: {dft: {mgrid: {cutoff: 500}}}` |
| 段参数写在键 `_` 下 | `xc_functional: {_: PBE}`、`ot: {_: true}` |
| `&COORD` 这类“无名行”写在键 `"*"` 下 | `coord: {"*": ["O 5.0 5.0 5.12", ...]}` |
| 重复段可以写成列表，参数唯一时也可以写成以参数为键的映射 | `kind: [{_: H, ...}, {_: O, ...}]` 或 `kind: {H: {...}, O: {...}}` |
| 多值关键字写成列表 | `abc: [10.0, 10.0, 10.0]` |
| 段名与关键字同名时，段的键加 `+` 前缀 | `+print`（极少用到） |

完整例子见 [inp_from_yaml.yaml](../examples/17-automation-tools/inp_from_yaml.yaml)，其中的片段：

```yaml
force_eval:
  dft:
    scf:
      ot:
        _: true
        minimizer: diis
        preconditioner: full_all
    xc:
      xc_functional:
        _: PBE
  print:
    forces:
      _: "ON"           # print keys take ON/OFF/LOW/MEDIUM..., not true/false
```

```bash
tocp2k --yaml inp_from_yaml.yaml > h2o_from_yaml.inp
cp2k.psmp -c h2o_from_yaml.inp
fromcp2k -f yaml h2o_from_yaml.inp > roundtrip.yaml       # 反向转换，可用于比较
fromcp2k -f aiida-cp2k-calc h2o_from_yaml.inp > run_aiida.py   # aiida-cp2k 脚本模板（§8）
```

YAML 的好处是可以用 Python 程序化地修改，也方便做版本管理和比较（`diff`）。缺点是转换时注释和预处理指令都会丢失，输出的段落顺序也与原来不同。

### 5.4 cp2kgen：参数扫描

[cp2kgen_cutoff_scan.sh](../examples/17-automation-tools/cp2kgen_cutoff_scan.sh) 先检查模板，再用 `--zip` 把 CUTOFF 和 PROJECT 名一一配对，然后把每个输入放进单独的目录：

```bash
cp2kgen --zip h2o_cutoff_template.inp \
  "force_eval/dft/mgrid/cutoff=[300,400,500,600]" \
  "global/project_name=[c300,c400,c500,c600]"
# Writing 'h2o_cutoff_template-cutoff_300-project_name_c300.inp'...
```

| 选项/规则 | 说明 |
|---|---|
| `"path=[a,b,c]"` | 方括号列表，产生多个值；不带方括号就只替换为一个值 |
| 默认（`--no-zip`） | 多个表达式之间做笛卡尔积，例如 4 个 CUTOFF × 3 个 REL_CUTOFF = 12 个输入 |
| `--zip` | 逐项配对，各个列表的长度必须相同 |
| `-E VAR=val` | 在解析模板前给预处理器变量赋值 |
| 输出文件名 | `<模板名>-<关键字>_<值>-...inp`，写在**当前目录** |
| 副作用 | `@SET`/`@INCLUDE` 已展开、注释已丢失；文件头带有生成命令 |

```bash
DRYRUN=1 bash cp2kgen_cutoff_scan.sh            # 只生成，先检查 scan_cutoff/c*/run.inp
bash cp2kgen_cutoff_scan.sh                     # 生成 + cp2k -c + 运行 + 汇总
cat scan_cutoff/summary.dat                     # CUTOFF, E_total[Ha], dE[meV]
```

判断收敛的标准（与[第4章](04-cutoff-convergence.md)相同）：相对最大 CUTOFF 的能量差 $|E(c)-E(c_\mathrm{max})|/N_\mathrm{atom} < 1$ meV/atom，同时每级网格分到的高斯函数数目（`count for grid`）分布合理。本仓库的 [conv_cutoff.py](../scripts/conv_cutoff.py) 把 CUTOFF 和 REL_CUTOFF 的两步扫描、网格计数提取和作图都做成了一条命令（在示例目录中用 `python ../../scripts/conv_cutoff.py -h` 查看用法）。

### 5.5 Python API

```python
from cp2k_input_tools.parser import CP2KInputParserSimplified
from cp2k_input_tools.generator import CP2KInputGenerator

with open("h2o_cutoff_template.inp") as fh:
    tree = CP2KInputParserSimplified().parse(fh)          # nested dict (validated, 2024.1 schema)
tree["force_eval"]["dft"]["mgrid"]["cutoff"] = 600
with open("h2o_600.inp", "w") as fh:
    for line in CP2KInputGenerator().line_iter(tree):
        fh.write(line + "\n")
```

`CP2KInputParser` 和 `CP2KInputParserSimplified` 都接受 `xmlspec=` 参数，可以指向 `cp2k.psmp --xml` 导出的文件。

---

## 6. Multiwfn 生成 CP2K 输入

Multiwfn（<http://sobereva.com/multiwfn/>）主功能 100 的子功能 2 可以生成多种程序的输入。其中的 CP2K 功能也可以在主菜单直接输入 `cp2k` 进入。它会根据载入的结构和你在菜单中的选择，写出一个带大量注释、基本可以直接运行的输入。官方介绍见 sobereva 的文章（<http://sobereva.com/587>，2021 年首发，2025-03 更新），本节只概括流程和需要注意的地方。

### 6.1 流程与支持的结构格式

1. 启动 Multiwfn，载入含结构的文件；
2. 在主菜单输入 `cp2k`；
3. 输入新输入文件的路径，直接回车则使用默认名（与结构文件同名，后缀为 `.inp`）；
4. 在菜单中设置理论方法、基组、任务类型、k 点、色散等；
5. 选 `0` 生成文件。

能**同时提供坐标和晶胞**的格式（缺少晶胞信息时，Multiwfn 把体系当作孤立分子，在分子外加一定距离作为盒子）：

| 格式 | 要求 |
|---|---|
| cif | — |
| pdb / pqr | 含 `CRYST1` 行 |
| CP2K `.inp` / `.restart` | 后缀必须是这两种 |
| GROMACS `.gro` | — |
| mol2 | 含 `@<TRIPOS>CRYSIN` |
| Gaussian `.gjf`（含 `Tv`）、PBC 任务的 `.fch` | GaussView 可以直接建周期模型 |
| xyz | 第二行写 `Tv_1: ... Tv_2: ... Tv_3: ...`，或 extended xyz 的 `Lattice="..."` |
| VASP POSCAR/CHGCAR/LOCPOT 等 | **文件名必须含相应字样**，如 `Si.POSCAR` |
| molden（含 `[Cell]`）、mwfn、QE 输入、Turbomole | 详见 sobereva.com/587 |

### 6.2 常用菜单项

下表按 sobereva.com/587（2025-03 版）整理。**编号可能随 Multiwfn 版本变化，请以屏幕显示为准**：

| 选项 | 作用 | 常用的下一步 |
|---|---|---|
| `1` | 理论方法 | 例如 `-3` 为 PBEsol；也有 B97M-rV、GFN1-xTB 等 |
| `2` | 基组 | 例如 TZVP-MOLOPT-GTH；选高斯基组（如 6-31G*）时使用 GAPW |
| `3` | 色散校正 | `2` 为 DFT-D3(BJ) |
| `4` | 开启 OT | — |
| `6` | 占据数展宽（smearing） | 金属必选 |
| `8` | k 点 | 输入如 `3,3,3` |
| `-1` | 任务类型 | `3` 定胞优化，`4` 原子+晶胞优化，`5` 振动分析，`6` MD，`7` 过渡态（Dimer） |
| `-2` | 输出 molden 文件 | 供 Multiwfn 做波函数分析（[第11章](11-cube-charge-workfunction.md)） |
| `-4` | 原子电荷 | 例如 REPEAT |
| `-7` | 周期性 | 如 `XY`（二维材料、表面） |
| `-9` | 其他设置 | 例如晶胞重复次数（生成 `MULTIPLE_UNIT_CELL`） |
| `-11` | 几何操作界面 | 例如构建超胞 |
| `9`（选了优化任务后） | 冻结原子 | 例如 `1,5,9-12` |
| `0` | 生成输入文件 | — |

例：用本章的 `Si.POSCAR` 生成 6×6×6 k 点的单点输入，以及 Cu 晶胞优化（PBEsol + smearing）的输入：

```text
Multiwfn Si.POSCAR        Multiwfn Cu.cif
cp2k                      cp2k
Si_sp.inp                 Cu_cellopt.inp
8                         1
6,6,6                     -3
0                         8
                          3,3,3
                          -1
                          4
                          6
                          0
```

批量生成时，可以把按键序列写进文本文件，用 `Multiwfn x.cif < keys.txt` 喂给 Multiwfn。序列请先交互式走一遍确认。

需要 cube 文件（电子密度、静电势、轨道）时，先在菜单里找相应的输出选项（sobereva.com/587 没有给出固定编号，以屏幕显示为准）。也可以生成输入后，按[第11章](11-cube-charge-workfunction.md)手动加入 `&E_DENSITY_CUBE`、`&V_HARTREE_CUBE` 或 `&MO_CUBES`。要做波函数分析，用 `-2` 输出 molden 文件更直接：CP2K 算完后会生成 `<项目名>-MOS-1_0.molden`，可直接载入 Multiwfn。

### 6.3 生成后必须做的检查

Multiwfn 生成的输入“完成度很高”，但仍然要针对你的体系和 CP2K 版本做检查：

```bash
cd CP2K/examples/17-automation-tools
python inp_migration_check.py Si_sp.inp       # 版本相关写法的静态检查（默认目标 2026.2）
cp2k.psmp -c Si_sp.inp                         # 最终的语法检查
grep -n "CUTOFF\|REL_CUTOFF\|EPS_SCF\|BASIS_SET\|POTENTIAL\|KPOINTS" Si_sp.inp
```

| 检查项 | 原因 | 处理 |
|---|---|---|
| `&MOTION/&CELL_OPT/TYPE` | 2026.2 删除了该关键字，运行时报 `The specified keyword 'TYPE' is not available anymore` | 删掉这一行（2026.2 固定使用 DIRECT_CELL_OPT），见[第7章](07-geo-cell-opt.md) |
| `&DFT/&PRINT/&PDOS` | 2026.2 改为 `&PRINT/&DOS/&PDOS`，旧写法报 `unknown subsection PDOS of section PRINT` | 按[第10章](10-dos-pdos-band.md)改写；反过来，新写法放到 2026.1 及以前的版本也会报错 |
| `CUTOFF`/`REL_CUTOFF`/`EPS_SCF` | Multiwfn 给的是经验值 | 按[第4章](04-cutoff-convergence.md)做收敛测试 |
| 基组与赝势 | 默认是 `DZVP-MOLOPT-SR-GTH` 一类的老基组 | 精度要求高时换 TZV2P 或 UZH 基组，并核对 -qN 一致 |
| 杂化泛函 | 周期体系会自动加截断库仑势和 `SCREEN_ON_INITIAL_P` | 核对 `CUTOFF_RADIUS` ≤ 最小晶胞宽度的一半（k 点计算按 Born–von Kármán 超胞），先用 PBE 收敛的波函数重启；必读 sobereva.com/690 和[第9章](09-hybrid-admm.md) |
| k 点 + OT | 2026.2 及以前不支持 | 改用对角化 + 混合（+ smearing） |
| `PERIODIC` | 表面、二维材料要设 `XY` 等 | `&CELL` 与 `&POISSON` 的 PERIODIC 要一致 |

[inp_migration_check.py](../examples/17-automation-tools/inp_migration_check.py) 只用 Python 标准库，能静态检查上表中与版本相关的项目（`--target 2025.2` 可检查“新语法放到旧版本”的问题）、PERIODIC 一致性、-qN 一致性和 HFX 截断半径。它是 `cp2k.psmp -c` 的补充，不能代替后者。

> [!IMPORTANT]
> 使用 Multiwfn 生成输入并发表论文时，作者要求在文中注明输入借助 Multiwfn 生成，并引用 Multiwfn 的论文：T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012)，<https://doi.org/10.1002/jcc.22885>；T. Lu, *J. Chem. Phys.* **161**, 082503 (2024)，<https://doi.org/10.1063/5.0216272>。

### 6.4 其他输入生成器

| 工具 | 说明 |
|---|---|
| 本仓库 [cp2k_inpgen.py](../scripts/cp2k_inpgen.py) | 从 xyz/cif/POSCAR 生成完整输入，自动写 CELL/PERIODIC，核对 UZH 基组与赝势的 -qN；`--version`（默认 2026.2）决定输出新的还是旧的 DOS/PDOS 语法（在示例目录中用 `python ../../scripts/cp2k_inpgen.py -h` 查看用法） |
| [CP2K_Forge](https://github.com/ljpyztd-cyber/CP2K_Forge) | Windows 图形界面，底层调用 Multiwfn，面向 2026.1/2026.2；同样需要引用 Multiwfn |
| [cp2kkit](https://github.com/CommonCE/cp2kkit) | 类似 VASPKIT 的工具箱（读输出、cube 处理、生成输入），目前是预发布版，没有许可证 |
| [xcp2k](https://github.com/superstar54/xcp2k) | ASE 风格的接口，可以直接操作全部输入段 |
| pymatgen `pymatgen.io.cp2k` | `Cp2kInput`/`Cp2kOutput` 等类，可从 `Structure` 生成输入；它内置的关键字同样可能落后于 2026.2，生成后照样要 `cp2k.psmp -c`。pymatgen 的通用用法见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md) |

---

## 7. 输出解析：cp2kdata、cp2k-output-tools 与 dpdata

### 7.1 cp2kdata

```bash
python -m pip install cp2kdata          # PyPI 0.7.3，同时安装 dpdata 和插件
# CP2K >= 2025.1 的输出：改装 GitHub 主分支（原因见下面的 WARNING）
python -m pip install --upgrade "git+https://github.com/robinzyb/cp2kdata.git"
```

> [!WARNING]
> CP2K 2025.1 新增了 `&FORCE_EVAL/&PRINT/&PROGRAM_RUN_INFO/ENERGY_UNIT`，总能量行随之从 `ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:` 变成 `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]`（2026.2 源码 `force_env_methods.F`，单位名后面没有冒号）。cp2kdata 0.7.3（PyPI 上的最新版）用正则 `energy \S{6}:` 匹配这一行，因此读 2025.1 及以后的输出时 `get_energies_list()` 返回 `None`，dpdata 的 `cp2kdata/e_f` 插件随后会因 `None * 单位换算` 而报 `TypeError`。GitHub 主分支已经放宽了这个正则，也增加了对 `FORCES| Atomic forces` 新格式的解析。无论装哪个版本，都请先用一个已知结果的算例核对能量和力。

```python
from cp2kdata import Cp2kOutput
out = Cp2kOutput("run.out")             # MD 等任务还会读同目录的 -1.ener / -pos-1.xyz / -frc-1.xyz
print(out)                               # summary: version, run type, atoms, frames, force/stress present
print(out.get_version_string(), out.get_run_type())
e = out.get_energies_list()             # Hartree
f = out.get_atomic_forces_list()        # Hartree/bohr, shape (nframes, natoms, 3)
info = out.get_geo_opt_info()           # GEO_OPT: per-step convergence data
freqs = out.get_vib_freq_list()         # VIBRATIONAL_ANALYSIS: cm^-1
```

| 类 / 命令 | 用途 |
|---|---|
| `Cp2kOutput` | 主输出：能量、力、应力、Mulliken 电荷与自旋、优化信息、振动频率 |
| `Cp2kCube` | cube 文件：`get_pav`（平面平均）、`get_mav`（宏观平均）、`get_integration`、`+`/`-` 做差分、`write_cube` |
| `Cp2kPdos` | `.pdos` 文件的展宽和求和（文档基于旧的 `&PRINT &PDOS` 输出格式） |
| `cp2kdata cube view --cube_file x.cube` | 在终端里直接看平面平均曲线（适合没有图形界面的集群） |

> [!WARNING]
> cp2kdata 0.7.x 拒绝解析 `PRINT_LEVEL LOW` 的输出，会报 `please provide cp2k output file with MEDIUM print level`。单点计算还要打印原子力（`&FORCE_EVAL/&PRINT/&FORCES ON`）。CP2K 每次改输出格式，解析器都可能失效（上面的能量行就是一例），所以换 CP2K 版本后先用一个已知结果的算例验证解析是否正确。`Cp2kPdos` 是按旧的 `&PRINT &PDOS` 文件写的，能否读 2026.2 `&DOS` 输出的文件尚未验证；本仓库的 [pdos_plot.py](../scripts/pdos_plot.py) 同时支持新旧两种 PDOS 格式。

### 7.2 cp2k-output-tools

```bash
python -m pip install cp2k-output-tools
cp2kparse run.out                               # JSON: version, energies, Mulliken, timings ...
cp2kparse -f yaml run.out
cp2kparse -k "energies/total force_eval" run.out   # 只取一个值
cp2kparse -f highlight run.out | less -R        # 高亮显示哪些行被解析器匹配到了
cp2k_bs2csv Si.bs                               # band structure -> Si.bs.set-N.csv
xyz_restart_cleaner md-pos-1.xyz md-pos-clean.xyz   # 去掉 MD 续算产生的重复帧
cp2k_pdos --help                                # PDOS 高斯展宽 -> CSV
```

PyPI 上的最新版本是 0.5.0（2021 年），对新版 CP2K 输出中的某些段落可能不认识。`-f highlight` 能直观显示哪些行解析成功了，适合用来确认。

### 7.3 dpdata：生成 DeePMD 训练集

cp2kdata 安装后会注册两个 dpdata 格式：

| fmt | 输入 | 要求 |
|---|---|---|
| `cp2kdata/e_f` | 单个 `ENERGY_FORCE` 输出文件 | `PRINT_LEVEL MEDIUM`，`&FORCES ON`；需要维里时加 `STRESS_TENSOR ANALYTICAL` 和 `&STRESS_TENSOR ON` |
| `cp2kdata/md` | MD 目录 + `cp2k_output_name=` | `-pos-1.xyz`、`-frc-1.xyz`、`-1.ener`（以及 NPT 时的 `-1.cell`）的输出频率必须相同（`&EACH MD n`） |

```python
import dpdata
s = dpdata.LabeledSystem("c500/run.out", fmt="cp2kdata/e_f")
md = dpdata.LabeledSystem("aimd_dir", cp2k_output_name="md.out", fmt="cp2kdata/md")
s.to_deepmd_npy("dp_data/sp")                  # or: s.to("deepmd/npy", "dp_data/sp")
md.to_deepmd_npy("dp_data/md")                 # systems with identical atoms can be joined with +
```

示例脚本 [cp2kdata_to_dpdata.py](../examples/17-automation-tools/cp2kdata_to_dpdata.py) 把这些步骤封装成命令行工具：

```bash
python cp2kdata_to_dpdata.py ef scan_cutoff/c*/run.out --deepmd dp_h2o
python cp2kdata_to_dpdata.py md aimd_dir --output md.out --deepmd dp_md
ls dp_h2o/*/set.000/        # box.npy coord.npy energy.npy force.npy
```

成功的标准：每个文件报告的帧数大于 0（SCF 未收敛的单点会被插件当作空数据跳过），`type_map.raw` 中的元素顺序符合预期，能量单位已换算为 eV。

> [!CAUTION]
> 训练集里的所有数据必须使用**同一套** DFT 参数（泛函、色散、基组、CUTOFF、EPS_SCF、k 点）。从 AIMD 中抽帧时，要隔开足够长的相关时间。不同 CUTOFF 的数据混在一起会让模型学到噪声。AIMD 的设置和轨迹后处理见[第12章](12-aimd.md)（本仓库的 [`md_analyze.py`](../scripts/md_analyze.py) 可以画温度和守恒量漂移）。

### 7.4 选哪个解析器

| 需求 | 推荐 |
|---|---|
| 零依赖、在集群上快速查看 SCF 是否收敛、优化每一步、WARNING 汇总 | 本仓库 [cp2k_outparse.py](../scripts/cp2k_outparse.py)（支持 `--json`；在示例目录中用 `python ../../scripts/cp2k_outparse.py -h` 查看用法） |
| Python 里拿到数组（能量、力、应力、频率），以及 cube/pdos 处理 | cp2kdata |
| 输出转成 JSON/YAML 存档，能带转 CSV | cp2k-output-tools |
| 机器学习势训练集 | cp2kdata + dpdata |
| 几何优化收敛图 | 本仓库 [plot_opt.py](../scripts/plot_opt.py) |

---

## 8. AiiDA 与 aiida-cp2k 简介

[AiiDA](https://aiida.readthedocs.io/projects/aiida-core/en/stable/) 把每一次计算的输入、输出、代码版本和计算机都记录在数据库里（provenance graph），由守护进程（daemon）负责提交、监控、取回和失败重试。[aiida-cp2k](https://github.com/aiidateam/aiida-cp2k) 2.2.1（2026-11-03 发布，需要 aiida-core ≥2.7）提供：

| 入口 | 名称 | 说明 |
|---|---|---|
| CalcJob | `cp2k`（`Cp2kCalculation`） | 输入：`parameters`（Dict）、`structure`、`kpoints`、`file`（附加文件）、`parent_calc_folder`（重启）、`settings`；输出：`output_parameters`、`output_structure`、`output_trajectory`、`output_bands` |
| WorkChain | `cp2k.base`（`Cp2kBaseWorkChain`） | 包装 CalcJob，对部分已知错误（如计算未完成）自动重启 |
| 解析器 | `cp2k_base_parser`、`cp2k_advanced_parser` | 能量、优化步、能带等 |

最小流程（[aiida_h2o_energy.py](../examples/17-automation-tools/aiida_h2o_energy.py)）：

```bash
python -m pip install aiida-cp2k
verdi presto                                   # 快速建立 profile（正式使用时配 PostgreSQL/RabbitMQ）
verdi computer setup                           # 交互式：集群地址、调度器（core.slurm / core.pbspro）
verdi computer configure core.ssh <computer>
verdi code create core.code.installed          # 交互式：可执行文件 cp2k.psmp，默认插件 cp2k
python aiida_h2o_energy.py cp2k@mycluster --nprocs 4
verdi process list -a
verdi calcjob outputcat <PK> | tail            # 查看 aiida.out
```

参数字典的格式就是 §5.3 的树：段名为键，段参数写在 `"_"` 下，重复段写成列表：

```python
"SCF": {"EPS_SCF": 1.0e-7, "OT": {"_": "ON", "PRECONDITIONER": "FULL_ALL"}},
"SUBSYS": {"CELL": {"PERIODIC": "NONE"},
           "KIND": [{"_": "O", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6", "POTENTIAL": "GTH-PBE-q6"}]},
```

已有的 `.inp` 可以用 `fromcp2k -f aiida-cp2k-calc x.inp` 直接生成这样的脚本模板。注意以下几点（来自 aiida-cp2k 2.2.1 源码）：

- 插件强制使用 `PROJECT aiida`，因此输出文件是 `aiida.out`、`aiida-1.restart`、`aiida-pos-1.xyz` 等；
- 插件从 `StructureData` 写入晶胞向量 A/B/C 和坐标文件（`&TOPOLOGY COORD_FILE_NAME`），但**不会把 pbc 转换成 PERIODIC**，孤立分子要自己在参数里写 `&CELL PERIODIC NONE` 和 `&POISSON`；
- 基组和赝势默认取远程 CP2K 安装目录里的 data；如果希望把它们也纳入 provenance，就用 `file` 命名空间上传 `SinglefileData`。

**什么时候值得用 AiiDA**：上百上千个结构的高通量筛选、需要自动重试和可追溯性的项目、团队共享数据库。如果只是十几个计算，§3/§5 的脚本加 §9 的作业模板就够了。

---

## 9. 作业脚本：Slurm 与 PBS

示例目录提供了两个通用模板，本仓库的 [job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh)、[job_cp2k_pbs.sh](../scripts/job_cp2k_pbs.sh) 是功能更全的版本（作业调度器的基础见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)）。

| 功能 | `job_slurm.sh` | `job_pbs.sh` |
|---|---|---|
| 提交单个输入 | `sbatch job_slurm.sh x.inp` | `qsub -v INP=x.inp job_pbs.sh` |
| 运行 ASE 驱动脚本 | `sbatch job_slurm.sh --ase ase_si_eos.py` | `qsub -v ASE_SCRIPT=ase_si_eos.py job_pbs.sh` |
| MPI × OpenMP | `--ntasks-per-node` × `--cpus-per-task` | `select=...:mpiprocs=..:ompthreads=..`，并设 `NTHREADS` |
| 提交前检查 | `cp2k.psmp -c` 的输出里没有 `SUCCESS` 就退出 | 同左 |
| 结束后检查 | 有无 `PROGRAM ENDED AT`，打印 MPI 进程数和最终能量，没完成时列出重启文件 | 同左 |

`job_slurm.sh` 的关键部分：

```bash
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}   # srun >= 22.05 no longer inherits --cpus-per-task
chk=$(cp2k.psmp -c "$inp" 2>&1 || true)               # capture first: avoids SIGPIPE under pipefail
grep -q "SUCCESS, the input could be parsed" <<< "$chk" || exit 1
rm -f "$out"                                         # -o appends to an existing file
srun --cpu-bind=cores cp2k.psmp -i "$inp" -o "$out"
grep -q "PROGRAM ENDED AT" "$out"
grep -m1 "DBCSR| MPI: Number of processes" "$out"   # must equal the requested ranks
```

扫描目录用数组作业（每个子目录一个任务）：

```bash
#!/bin/bash
#SBATCH --array=0-7
#SBATCH --ntasks=8
dirs=(scan_cutoff/c*/)
cd "${dirs[$SLURM_ARRAY_TASK_ID]}"
srun cp2k.psmp -i run.inp -o run.out
```

> [!TIP]
> 在作业中运行 ASE 脚本时，MPI 启动器要写在 `ASE_CP2K_COMMAND` 里（模板默认用 `mpirun -np $SLURM_NTASKS cp2k.psmp -s`）。如果改用 `srun`，要确认 stdin 能送到 0 号进程，具体看所用 Slurm 版本 `man srun` 中 `--input` 的说明。不要在计算节点上用 `module load` 之类会向 stdout 打印信息的命令包住 CP2K。

---

## 10. 串起来：一个典型的自动化流程

| 步骤 | 工具 | 产物 | 检查点 |
|---|---|---|---|
| 1. 建模 | Multiwfn / [`cp2k_inpgen.py`](../scripts/cp2k_inpgen.py) / ASE | `x.inp` | `inp_migration_check.py`、`cp2k.psmp -c` |
| 2. 收敛测试 | `cp2kgen` 或 [`conv_cutoff.py`](../scripts/conv_cutoff.py)、[`conv_kpoints.py`](../scripts/conv_kpoints.py) | `summary.dat` | < 1 meV/atom |
| 3. 生产计算 | `job_slurm.sh` / `job_pbs.sh` / AiiDA | `.out`、`-1.restart` | `PROGRAM ENDED AT`、无 `SCF run NOT converged` |
| 4. 结构与性质 | ASE（EOS、约束优化）或 CP2K 原生 `GEO_OPT`/`CELL_OPT` | `si_eos.dat`、优化后的结构 | 与文献/实验对比 |
| 5. 解析与归档 | [`cp2k_outparse.py`](../scripts/cp2k_outparse.py)、cp2kdata、cp2k-output-tools | JSON/CSV | 单位、帧数 |
| 6. 训练集（可选） | dpdata | `deepmd/npy` | 参数一致、帧数、元素顺序 |

---

## 参数选择与收敛建议

- **工具的版本和 CP2K 的版本要分别记录**：在作业日志里写下 `cp2k.psmp -v`、`ase.__version__`、`pip show cp2k-input-tools cp2kdata` 的输出。工具升级可能悄悄改变它生成的输入（例如 ASE 的默认命令、Multiwfn 的默认参数）。
- **ASE 的 `cutoff` 以 eV 为单位**：500 Ry 要写成 `500 * Rydberg`。直接写 `cutoff=500` 只有 500 eV ≈ 37 Ry，计算会严重不准，却不会报错。
- **ASE 计算的 SCF 精度**：ASE 的优化器依赖力的精度。`fmax=0.01 eV/Å` 时 `EPS_SCF` 建议 ≤ 1e-7；EOS、声子等差分类任务用 1e-8。
- **cp2kgen 扫描时只改一个变量**，其余参数（`EPS_DEFAULT`、`REL_CUTOFF`、k 点）固定不变；每个输入放在单独的目录里，避免重启文件互相覆盖。
- **训练集数据**：所有标签计算的 DFT 参数完全一致；输出用 `PRINT_LEVEL MEDIUM`；必须打印力（需要时还有应力）；SCF 必须收敛，不要使用 `IGNORE_CONVERGENCE_FAILURE`。
- **Multiwfn 生成的输入**：把它当作“初稿”，CUTOFF、REL_CUTOFF、k 点、基组都要按第4、6章重新收敛。
- **ASE vs CP2K 原生**：常规优化、MD、NEB 优先用 CP2K 原生功能（有重启文件，效率高）；ASE 适合跨程序比较、自定义约束、与 MLIP 混合的工作流。

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `/bin/sh: cp2k_shell: not found` | ASE ≤ 3.28 的默认命令是 `cp2k_shell`，而 CMake 构建的 CP2K 没有这个可执行文件 | 按 §4.2 用 wrapper 或符号链接，或升级 ASE ≥ 3.29 并设 `ASE_CP2K_COMMAND="cp2k.psmp -s"` |
| `AssertionError`，位于 `assert 'cp2k_shell' in command` | ASE ≤ 3.28 要求命令字符串含 `cp2k_shell` | 同上 |
| `RuntimeError: Did not receive "* READY" after starting CP2K shell.`（ASE 3.29） | CP2K 启动失败（环境没加载、找不到动态库、MPI 不匹配），或 stdin 没有送到 0 号进程 | 在同一终端手动运行命令，看报错；`which cp2k.psmp`、`ldd`；把 `srun` 换成 `mpirun` |
| ASE ≤ 3.28 在 `Cp2kShell.__init__`/`expect`/`recv` 里抛出没有说明的 `AssertionError` | 同上；或者 `* READY` 之前混入了其他输出（3.28 不跳过启动信息） | 先按上一行排查；去掉 wrapper 和环境脚本中打印到 stdout 的内容（重定向到 stderr），或升级到 ASE 3.29 |
| `Cannot determine version of CP2K shell.` | `* READY` 之后、版本号之前混入了其他输出，或 CP2K 太旧（错误信息要求 CP2K 3.0 以上） | 去掉所有打印到 stdout 的内容（需要的话重定向到 stderr）；`debug=True` 查看收到了什么 |
| `Exception: Keyword CUTOFF already present in section MGRID` | 模板与计算器参数重复 | 从模板删掉这个关键字，或把对应的计算器参数设为 `None` |
| ASE 卡住不动（大体系、psmp） | 通过 stdin 传大量坐标时，部分 MPI 实现会阻塞 | `set_pos_file=True`（CP2K ≥ 2024.2） |
| ASE 计算结果比 CP2K 直接算的差很多 | `cutoff` 按 Ry 填了，或者 `pseudo_potential='auto'` 选到了不同的赝势 | `cutoff=N*Rydberg`；用模板 `&KIND` 明确指定基组和赝势，并对比 `<label>.inp` |
| `cp2klint` 对 2024.1 之后的新写法（`&DOS` 内嵌的 `&PDOS`/`&CURVE`、`&XTB` 的 `GFN_TYPE`）报错 | 内置的是 2024.1 的输入结构 | 以 `cp2k.psmp -c` 为准；或者 `cp2k.psmp --xml` 后用 `fromcp2k --xml` 检查 |
| `cp2klint` 通过，但 CP2K 2026.2 报 `The specified keyword 'TYPE' is not available anymore` | 同上，旧语法在 2024.1 中合法 | 删除 `&CELL_OPT` 里的 `TYPE`；先运行 `inp_migration_check.py` |
| `unknown subsection PDOS of section PRINT` | 2026.2 的 DOS 输入已重构 | 改成 `&PRINT/&DOS/&PDOS`（[第10章](10-dos-pdos-band.md)） |
| REFTRAJ 报 `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ`（或 `EVAL_FORCES`） | 2026.2 删除了这两个关键字 | 改用 `EVAL ENERGY` 或 `EVAL ENERGY_FORCES`（2024.1 起就有，新旧版本通用；见[第12章](12-aimd.md)） |
| cp2kgen/tocp2k 的单位报错（如 `[fs]` 被识别成 femtosiemens） | cp2k-input-tools 0.9.3 + pint 0.26 的单位解析问题（本仓库 2026-11-10 测试时发现，与 Python 哈希种子有关） | 去掉方括号单位，改用该关键字的默认单位；或者固定 `PYTHONHASHSEED`；以 `cp2k.psmp -c` 为准 |
| cp2kgen 生成的输入找不到（文件名不对） | 输出写在**当前目录**，文件名是 `<模板名>-<关键字>_<值>.inp` | 在模板所在目录运行，或者按 `cp2kgen_cutoff_scan.sh` 的方式移动文件 |
| cp2kdata 报 `please provide cp2k output file with MEDIUM print level` | 输出是 `PRINT_LEVEL LOW` | 改成 MEDIUM 重新计算 |
| cp2kdata 的 `get_energies_list()` 返回 `None`；dpdata 报 `TypeError: unsupported operand type(s) for *: 'NoneType' and 'float'` | PyPI 上的 cp2kdata 0.7.3 不认识 CP2K ≥ 2025.1 的能量行 `energy [hartree]` | 改装 GitHub 主分支（§7.1），并用已知算例核对 |
| dpdata 得到 0 帧 | SCF 未收敛、没有打印力，或者 MD 各文件的输出频率不一致 | 检查 `SCF run NOT converged`；加 `&FORCES ON`；统一 `&EACH MD n` |
| Multiwfn 生成的输入在 2026.2 中报错 | 输入是按旧语法生成的，或者 Multiwfn 版本太旧 | 升级 Multiwfn；运行 `inp_migration_check.py` 和 `cp2k.psmp -c` |
| Multiwfn 生成的是孤立分子的盒子 | 载入的文件没有晶胞信息（例如普通 xyz、不含 `CRYST1` 的 pdb） | 改用 cif/POSCAR，或带 `Tv` 的 gjf |
| AiiDA 计算失败，输出中找不到基组 | 远程 CP2K 的 data 目录里没有该文件，或者文件名写错 | 核对 `BASIS_SET_FILE_NAME`；或者用 `file` 命名空间上传（`--with-data-files`） |
| AiiDA 算孤立分子时出现周期性镜像相互作用 | 插件不会根据 pbc 设置 PERIODIC | 在参数里写 `&CELL PERIODIC NONE` 和 `&POISSON PERIODIC NONE` + `POISSON_SOLVER MT` |
| 输出中出现多次 `DBCSR| MPI: Number of processes 1` | MPI 启动器与 CP2K 编译时用的 MPI 不一致，实际上运行了多个串行进程 | [第18章](18-troubleshooting.md)；`which mpirun` 与 `ldd cp2k.psmp` 要对应同一套 MPI |

---

## 练习

1. 在你的集群上分别用 ASE 3.28 和 3.29（各建一个 venv）运行 `ase_h2o_opt.py`：前者用 wrapper，后者直接用 `cp2k.psmp -s`，比较两者生成的 `h2o_ase.inp` 是否相同。
2. 把 `ase_h2o_opt.py` 的 `cutoff` 改成 `500`（忘记乘 `Rydberg`），观察能量和键长的变化。从 `h2o_ase.out` 里找出实际使用的 CUTOFF。
3. 用 `cp2kgen` 对 `h2o_cutoff_template.inp` 同时扫描 CUTOFF（300–800）和 REL_CUTOFF（40、60、80），比较 `--zip` 与默认（笛卡尔积）各生成多少个输入，并画出二维收敛图。
4. 用 `cp2k.psmp --xml` 导出你所用版本的 `cp2k_input.xml`，写一个含 2026.2 `&DOS/&PDOS` 的输入，分别用 `cp2klint`、`fromcp2k --xml cp2k_input.xml` 和 `cp2k.psmp -c` 检查，记录三者的结论。
5. 用 Multiwfn 由 `Si.POSCAR` 生成 `CELL_OPT` 输入（任务 `4`），然后运行 `inp_migration_check.py`，根据提示修改，直到 `cp2k.psmp -c` 通过。
6. 运行 `ase_si_eos.py`，k 点分别取 4、6、8，比较 $a_0$ 和 $B$ 的变化。再用 CP2K 原生的 `CELL_OPT`（[第7章](07-geo-cell-opt.md)）优化同一体系，比较两种方法得到的 $a_0$。
7. 从[第12章](12-aimd.md)的 AIMD 轨迹中每隔 20 步抽一帧，对这些帧做 `ENERGY_FORCE` 单点，然后用 `cp2kdata_to_dpdata.py` 生成 `deepmd/npy` 数据集，检查帧数和 `type_map.raw`。
8. 用 `fromcp2k -f aiida-cp2k-calc` 把[第2章](02-input-structure.md)的 H2O 输入转成 AiiDA 脚本，补上 `code` 和资源设置后提交，用 `verdi process show` 查看 provenance。

---

## 参考资料

**官方文档与源码**

- CP2K 命令行选项（`-c/-e/-E/-s/--xml`）见 2026.2 源码 `src/start/cp2k.F`：<https://github.com/cp2k/cp2k/blob/v2026.2/src/start/cp2k.F>
- CP2K 2026.2 发布说明（删除 Cython 绑定、`CELL_OPT TYPE`，DOS 重构）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- CP2K 输入参考（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT.html>；REFTRAJ：<https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/MD/REFTRAJ.html>
- CP2K 新的 DOS 接口：<https://manual.cp2k.org/trunk/methods/electronic_structure/dos.html>
- ASE CP2K 计算器（文档字符串即官方文档）：<https://gitlab.com/ase/ase/-/blob/master/ase/calculators/cp2k.py>
- cp2k-input-tools：<https://github.com/cp2k/cp2k-input-tools>
- cp2k-output-tools：<https://github.com/cp2k/cp2k-output-tools>
- cp2kdata：<https://github.com/robinzyb/cp2kdata>（dpdata 插件说明：<https://github.com/robinzyb/cp2kdata/blob/master/docs/dpdata_plugin.md>）；PyPI：<https://pypi.org/project/cp2kdata/>
- dpdata：<https://github.com/deepmodeling/dpdata>；文档：<https://docs.deepmodeling.com/projects/dpdata/en/master/>
- aiida-cp2k：<https://github.com/aiidateam/aiida-cp2k>；文档：<https://aiida-cp2k.readthedocs.io/en/latest/>
- AiiDA：<https://aiida.readthedocs.io/projects/aiida-core/en/stable/>
- Multiwfn 主页：<http://sobereva.com/multiwfn/>；用 Multiwfn 生成 CP2K 输入：<http://sobereva.com/587>；CP2K 杂化泛函计算要点：<http://sobereva.com/690>；CP2K 生成 molden 给 Multiwfn：<http://sobereva.com/651>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools)：cp2klint、cp2kgen、fromcp2k/tocp2k、cp2kget、语言服务器
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：cp2kparse、cp2k_bs2csv、cp2k_pdos、xyz_restart_cleaner
- [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata)：输出、cube、pdos 解析，以及 dpdata 插件
- [aiidateam/aiida-cp2k](https://github.com/aiidateam/aiida-cp2k)：`examples/single_calculations/` 下有 DFT、GEO_OPT、能带、重启等示例
- [KMNitesh05/cp2k_2_deepmdkit](https://github.com/KMNitesh05/cp2k_2_deepmdkit)：把 AIMD 的 pos/frc/ener/cell 转成 DeepMD npy，可与 dpdata 对照使用
- [superstar54/xcp2k](https://github.com/superstar54/xcp2k)：ASE 风格、覆盖全部输入段的接口
- [LouieSlocombe/ase-cp2k-noshell](https://github.com/LouieSlocombe/ase-cp2k-noshell)：不依赖 shell 模式的 ASE 计算器
- [svarner9/CP2K-Workshop](https://github.com/svarner9/CP2K-Workshop)：用 ASE + PySAGES 驱动 CP2K 做增强采样
- [ljpyztd-cyber/CP2K_Forge](https://github.com/ljpyztd-cyber/CP2K_Forge)：基于 Multiwfn 的 Windows 图形界面
- [CommonCE/cp2kkit](https://github.com/CommonCE/cp2kkit)：类似 VASPKIT 的 CP2K 工具箱（预发布版）

**本仓库脚本**：[cp2k_inpgen.py](../scripts/cp2k_inpgen.py)（生成输入）、[conv_cutoff.py](../scripts/conv_cutoff.py)/[conv_kpoints.py](../scripts/conv_kpoints.py)（收敛测试）、[cp2k_outparse.py](../scripts/cp2k_outparse.py)（零依赖解析）、[plot_opt.py](../scripts/plot_opt.py)（优化收敛图）、[restart_update.py](../scripts/restart_update.py)（从重启文件续算）、[job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh)/[job_cp2k_pbs.sh](../scripts/job_cp2k_pbs.sh)（作业模板）、[pdos_plot.py](../scripts/pdos_plot.py)、[band_plot.py](../scripts/band_plot.py)、[md_analyze.py](../scripts/md_analyze.py)、[cube_tools.py](../scripts/cube_tools.py)、[neb_plot.py](../scripts/neb_plot.py)、[vib_thermo.py](../scripts/vib_thermo.py)。它们都在 `CP2K/scripts/` 下，在示例目录中用 `python ../../scripts/<name> -h` 查看用法。

**关键文献**

- A. H. Larsen *et al.*, *J. Phys.: Condens. Matter* **29**, 273002 (2017)，ASE：<https://doi.org/10.1088/1361-648X/aa680e>
- S. P. Huber *et al.*, *Sci. Data* **7**, 300 (2020)，AiiDA 1.0：<https://doi.org/10.1038/s41597-020-00638-4>
- T. Lu, F. Chen, *J. Comput. Chem.* **33**, 580 (2012)，Multiwfn：<https://doi.org/10.1002/jcc.22885>
- T. Lu, *J. Chem. Phys.* **161**, 082503 (2024)，Multiwfn 综述：<https://doi.org/10.1063/5.0216272>
- H. Wang, L. Zhang, J. Han, W. E, *Comput. Phys. Commun.* **228**, 178 (2018)，DeePMD-kit：<https://doi.org/10.1016/j.cpc.2018.03.016>
- P. Haas, F. Tran, P. Blaha, *Phys. Rev. B* **79**, 085104 (2009)，半局域泛函的晶格常数（含 PBE 的 Si）：<https://doi.org/10.1103/PhysRevB.79.085104>
- F. Birch, *Phys. Rev.* **71**, 809 (1947)，Birch–Murnaghan 状态方程：<https://doi.org/10.1103/PhysRev.71.809>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](16-qmmm-xtb.md) | [返回目录](../README.md) | [下一章 →](18-troubleshooting.md)
