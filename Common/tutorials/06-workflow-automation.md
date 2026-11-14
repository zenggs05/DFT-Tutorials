# 第6章 工作流自动化与高通量：custodian、atomate2、jobflow、AiiDA（Workflow Automation & High-Throughput: custodian, atomate2, jobflow, AiiDA）

计算从几个结构增加到几十、几百个以后，靠手动提交、`grep` 报错、改 INCAR、重提就管不过来了。本章按“自动纠错 → 工作流编写 → 集群调度与数据库 → 完整溯源”的顺序，介绍 custodian（VaspErrorHandler、UnconvergedErrorHandler 自动处理 ZBRENT/EDDDAV 等报错）、jobflow（`@job`、`Flow`、`run_locally`）、atomate2（RelaxBandStructureMaker、PhononMaker、VaspLobsterMaker 与 `~/.atomate2.yaml` 中的 `VASP_CMD`）、jobflow-remote 1.x（Slurm 提交 + MongoDB 管理）、AiiDA（`verdi presto`、aiida-vasp/aiida-cp2k/aiida-gaussian 的 WorkChain），以及更轻量的 quacc 配方和 dpdispatcher 批量提交，最后按“个人课题 / 课题组高通量 / 需要完整溯源”三种场景给出选型建议。

> [!NOTE]
> 本章默认你已经会写 Slurm 作业脚本（[第2章](02-job-scheduler.md)），也会用 ASE/pymatgen 读写结构和输出（[第4章](04-ase-pymatgen.md)）。四个配套示例目录：[06-custodian-vasp-si](../examples/06-custodian-vasp-si/README.md)、[06-jobflow-atomate2](../examples/06-jobflow-atomate2/README.md)、[06-aiida-cp2k-vasp](../examples/06-aiida-cp2k-vasp/README.md)、[06-quacc-dpdispatcher](../examples/06-quacc-dpdispatcher/README.md)。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 自动纠错 | 用 custodian 包装 VASP，读懂 `custodian.json`，知道每个 handler 会改哪些 INCAR 参数 |
| 写工作流 | 用 jobflow 的 `@job`/`Flow` 把多个计算串成依赖图，并用 `run_locally` 在本地跑通 |
| 用现成工作流 | 用 atomate2 的 Maker 跑弛豫 + 能带、声子、LOBSTER，并通过 powerup 批量修改 INCAR |
| 上集群 | 配置 jobflow-remote 项目文件，让 runner 自动上传、提交、下载、入库 |
| 完整溯源 | 用 `verdi presto` 建库，用 aiida-cp2k/aiida-vasp/aiida-gaussian 运行并查看溯源图 |
| 轻量方案 | 用 quacc 配方或 dpdispatcher 完成“小批量、无数据库”的自动化 |
| 选型 | 根据课题规模、团队和数据要求选择工具组合 |

**适用场景**：收敛测试、掺杂/吸附位点扫描、成分筛选、声子和弹性常数这类“结构多、步骤固定”的计算；需要给论文或数据库提供可复现记录的项目。**不适用**：只算两三个结构的探索性工作。这时一个 Slurm 依赖链（[第2章](02-job-scheduler.md)第8节）加上 [job_status.py](../scripts/job_status.py)（`python ../scripts/job_status.py -h` 查看用法）就够了，引入数据库和守护进程反而更麻烦。

---

## 2. 原理速览

### 2.1 为什么需要自动纠错

设一批 $N$ 个独立计算，每个作业因为数值问题（ZBRENT、EDDDAV、SCF 不收敛……）失败的概率为 $p$，则全部一次成功的概率为

$$
P_{\text{all}} = (1-p)^N .
$$

取 $N=200$、$p=0.1$，$P_{\text{all}}\approx 7\times10^{-10}$，平均要手动处理 $Np = 20$ 个失败作业。如果有一个程序能识别报错并自动修正，每次修正的成功率为 $r$，最多尝试 $m$ 次，则单个作业最终失败的概率降为

$$
p_m = p\,(1-r)^m .
$$

取 $r=0.7$、$m=3$，$p_m \approx 2.7\times10^{-3}$，200 个作业中平均只剩约 0.5 个需要人工处理（这是一个示意模型，$p$、$r$ 因体系而异）。custodian 做的就是这件事。

### 2.2 工作流 = 有向无环图

工作流把计算步骤表示为有向无环图（DAG）$G=(V,E)$：节点 $v\in V$ 是一个作业（job），边 $(u\to v)\in E$ 表示 $v$ 要用 $u$ 的输出（结构、电荷密度、目录）。一个作业只有在所有父节点都完成后才能运行。工作流引擎负责三件事：

1. **依赖解析**：把“引用”（如“上一步弛豫后的结构”）在运行时替换为真实数据；
2. **状态管理**：记录每个作业处于 READY/RUNNING/COMPLETED/FAILED 中的哪个状态，失败后可重跑；
3. **数据存储**：把结果写入数据库，方便查询和复现。

### 2.3 工具分层

```text
 层次            作用                         代表工具
 ---------------------------------------------------------------------------
 4 溯源/数据库   每个输入输出都是数据库节点     AiiDA (PostgreSQL/SQLite + broker)
 3 调度/队列     批量上传、提交、监控、下载     jobflow-remote, dpdispatcher, FireWorks
 2 工作流定义    把步骤写成 DAG                jobflow + atomate2, quacc, AiiDA WorkChain
 1 单步纠错      运行程序、识别报错、改输入重跑  custodian (VASP/CP2K/Gaussian/QChem/LOBSTER)
 0 程序本体      VASP, CP2K, Gaussian, OpenMX   + Slurm/PBS/LSF
```

atomate2 的每个 VASP 步骤内部调用 custodian（第 1 层），工作流用 jobflow 描述（第 2 层），在集群上推荐交给 jobflow-remote（第 3 层）。AiiDA 自己覆盖第 2–4 层，纠错由各插件 WorkChain 中的 process handler 完成。quacc 只负责第 2 层，第 3 层可以接 Parsl、Dask、Prefect、jobflow 等引擎。

---

## 3. 安装与版本

截至 2026-11-14 的 PyPI 版本：

| 包 | 版本 | Python 要求 | 说明 |
|---|---|---|---|
| custodian | 2026.10.21 | ≥ 3.11 | 含 `custodian.vasp`、`custodian.cp2k`、`custodian.gaussian`、`custodian.qchem`、`custodian.lobster` |
| jobflow | 0.4.0 | ≥ 3.10 | atomate2 的工作流层 |
| atomate2 | 0.1.7 | **≥ 3.12** | 可选依赖：`[phonons]`（phonopy、seekpath）、`[lobster]`（lobsterpy、ijson） |
| jobflow-remote | 1.0.1 | ≥ 3.10（排除 3.14.1） | 需要 MongoDB |
| aiida-core | 2.9.4 | ≥ 3.10 | 2.6 起提供 `verdi presto` |
| aiida-vasp | 5.1.0 | ≥ 3.10 | WorkChain：`vasp.vasp`、`vasp.relax`、`vasp.bands`、`vasp.converge`、`vasp.neb` |
| aiida-cp2k | 2.2.1 | ≥ 3.10 | 需要 aiida-core ≥ 2.7 |
| aiida-gaussian | 2.2.0 | ≥ 3.7（声明值） | 用 cclib（`>=1.8,<=2.0`）解析输出；2024-09 发布，之后未更新 |
| quacc | 1.8.2 | ≥ 3.11 | 基于 ASE 的配方库 |
| dpdispatcher | 1.1.0 | ≥ 3.10 | 只做提交与文件传输 |

建议建两个独立的 conda 环境，避免依赖冲突：

```bash
# Materials Project stack: custodian / jobflow / atomate2 / jobflow-remote / quacc
conda create -n wf python=3.12 -y
conda activate wf
pip install custodian "atomate2[phonons,lobster]" jobflow-remote quacc dpdispatcher

# AiiDA stack
conda create -n aiida python=3.12 -y
conda activate aiida
pip install aiida-core aiida-vasp aiida-cp2k aiida-gaussian     # aiida-gaussian pulls in cclib (1.8.x)
```

pymatgen 生成 POTCAR 需要先配置赝势目录（atomate2、custodian 的部分修正、quacc 都会用到）：

```bash
pmg config -p /path/to/extracted_potpaw /path/to/pymatgen_potcars
pmg config --add PMG_VASP_PSP_DIR /path/to/pymatgen_potcars
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54      # only if you use the PBE.54 set
```

可以用本板块的 [env_check.sh](../scripts/env_check.sh) 检查 Python 版本、`PMG_VASP_PSP_DIR`、`ulimit -s` 等设置（`bash ../scripts/env_check.sh -h` 查看用法）。

> [!WARNING]
> 工作流工具更新很快，教程里的参数名可能在一两个版本后变化。生产项目开始时请用 `pip freeze > requirements-lock.txt` 记录版本（[第1章](01-linux-hpc.md)第11节），并在升级前先用小体系重跑一遍。

---

## 4. custodian：给 VASP 加上自动纠错

### 4.1 工作机制

custodian 的核心对象只有三种：

- **Job**：怎么运行程序。`VaspJob` 负责备份输入（`*.orig`）、启动 VASP（标准输出写入 `vasp.out`）、结束后给输出加后缀；
- **ErrorHandler**：怎么发现并修正错误。`is_monitor = True` 的 handler（如 `VaspErrorHandler`）在 VASP **运行期间**定期检查；VASP **结束后**，custodian 再用全部 handler（包括 `UnconvergedErrorHandler` 这类非 monitor handler）检查一遍；
- **Validator**：结束后检查结果是否可用（如 `vasprun.xml` 能否解析），不通过就报错退出。

运行流程：

```text
Custodian.run()
  for job in jobs:
      job.setup()                 # backup INCAR -> INCAR.orig, apply settings_override
      job.run()                   # Popen(vasp_cmd) > vasp.out
      every polling_time_step s:  # 10 s by default
          every monitor_freq polls: run monitor handlers (VaspErrorHandler, ...)
          if error: kill VASP -> tar.gz backup -> handler.correct() -> restart job
      after VASP exits: run all handlers -> correct & restart if needed
      validators -> job.postprocess()   # copy/move outputs to *<suffix>
  write custodian.json            # every correction is logged
```

每次修正前，handler 会把关键文件（INCAR、KPOINTS、POSCAR、OUTCAR、CONTCAR、OSZICAR、vasprun.xml、vasp.out、std_err.txt）打包成 `error.1.tar.gz`、`error.2.tar.gz`……，方便事后查看出错现场。

### 4.2 VASP handler 与修正策略

以下内容整理自 custodian 2026.10.21 源码 `custodian/vasp/handlers.py`。“报错原文”是 handler 在 `vasp.out` 中匹配的字符串。

| handler / 错误键 | 报错原文（匹配片段） | 自动修正（按触发次数递进） |
|---|---|---|
| `VaspErrorHandler` / `zbrent` | `ZBRENT: fatal error in bracketing`、`ZBRENT: fatal internal in` | CONTCAR→POSCAR；EDIFF 收紧（大于 1e-6 时设为 1e-6，否则除以 10，到 1e-8 为止）；NELMIN = 8。第二次触发时改 `IBRION = 1`。`vtst_fixes=True` 时改用 VTST FIRE（`IOPT = 7, IBRION = 3, POTIM = 0`） |
| `VaspErrorHandler` / `edddav` | `Error EDDDAV: Call to ZHEGV failed` | CONTCAR→POSCAR；删除 CHGCAR；`ALGO = All, ISEARCH = 1`；电子数少于 CPU 核数时 NCORE 加倍（并去掉 NPAR） |
| `VaspErrorHandler` / `eddrmm` | `WARNING in EDDRMM: call to ZHEGV failed` | CONTCAR→POSCAR；`ALGO` 为 Fast/VeryFast 时改 Normal，否则 POTIM 减半；删除 CHGCAR/WAVECAR |
| `VaspErrorHandler` / `brmix` | `BRMIX: very serious problems` | 依次尝试 `ISTART = 1` 重跑、`IMIX = 1`（Kerker 混合）、调整 KPOINTS 网格类型……（INCAR 中设置了 NELECT 的带电体系会跳过） |
| `VaspErrorHandler` / `tet` | `Tetrahedron method fails` 等 | `ISMEAR = 0, SIGMA = 0.05` |
| `VaspErrorHandler` / `subspacematrix` | `WARNING: Sub-Space-Matrix is not hermitian in DAV` | 先 `LREAL = .FALSE.`，再 `PREC = Accurate` |
| `VaspErrorHandler` / `too_few_bands` | `TOO FEW BANDS` | NBANDS 增大约 10% |
| `VaspErrorHandler` / `pssyevx`、`pdsyevx` | `ERROR in subspace rotation PSSYEVX` | `ALGO = Normal` |
| `UnconvergedErrorHandler` | VASP 结束后 `vasprun.xml` 显示未收敛 | 电子步：`ALGO` 按 VeryFast→Fast→Normal→All 逐级更换（`ISMEAR < 0` 时不换 All）；meta-GGA（`ISMEAR ≥ 0`）直接改 All；杂化泛函按 All→Damped（`TIME = 0.5`）；阶梯走完后改混合参数（`ISTART = 1, ALGO = Normal, NELMDL = -6, BMIX = 0.001, AMIX_MAG = 0.8, BMIX_MAG = 0.001`）。离子步：`IBRION = 1` 并 CONTCAR→POSCAR |
| `NonConvergingErrorHandler` | 最近 10 个离子步的电子步都达到 NELM | 与上面类似的 ALGO 阶梯，最后改线性混合参数 |
| `PotimErrorHandler` | 离子步能量突变（默认每原子 > 1 eV） | POTIM 减半；POTIM < 0.2 时改阻尼 MD（`IBRION = 3, SMASS = 0.75`） |
| `PositiveEnergyErrorHandler` | OSZICAR 中总能量为正 | `ALGO` 不是 Normal 时改 Normal；已是 Normal 且不是静态计算时 POTIM 减半 |
| `FrozenJobErrorHandler` | `vasp.out` 超过 `timeout` 秒未更新（默认 21600 s） | `ALGO = Fast` 时改 Normal；否则在离子弛豫（NSW > 0）中 CONTCAR→POSCAR 后重启；静态计算无法修复 |
| `MeshSymmetryErrorHandler` | VASP 结束且未收敛时，`vasp.out` 中有 `Reciprocal lattice and k-lattice belong to different class of lattices.` | Monkhorst-Pack 网格改为 Γ 中心；已经是 Γ 中心时设 `ISYM = 0` |
| `WalltimeHandler` | 剩余时间少于 max(3 × 最长离子步耗时, `buffer_time`) | 写 STOPCAR（`LSTOP = .TRUE.`），让 VASP 在当前离子步结束后正常退出 |
| `StdErrHandler` | `std_err.txt` 中的 `Allocation would exceed memory limit`、`internal error in GENERATE_KPOINTS_TRANS` 等 | 内存不足时 KPAR 减半；`GENERATE_KPOINTS_TRANS` 错误时设 `ISYM = 0` |

> [!IMPORTANT]
> custodian 的修正是“经验规则”，不保证物理上合理。例如 `ISMEAR = 0` 对金属的 DOS 可能不是最佳选择，`IBRION = 1` 对初始结构很差的体系可能更慢。论文数据要检查 `custodian.json`，确认最终用的参数，并在方法部分如实写出。手动排错思路见 [VASP 第18章](../../VASP/tutorials/18-postprocess-troubleshooting.md)。

### 4.3 Python 写法：弛豫 → 弛豫 → 静态

示例输入（[06-custodian-vasp-si](../examples/06-custodian-vasp-si/README.md)）是 Si 原胞全弛豫：

```text
# Si full relaxation (ions + cell) run under custodian
SYSTEM  = Si-custodian-demo
PREC    = Accurate
ENCUT   = 400        # > 1.3 x ENMAX(Si, about 245 eV): required for ISIF = 3 (Pulay stress)
EDIFF   = 1E-6
NELM    = 100
ALGO    = Fast       # custodian falls back to Normal / All if the SCF does not converge
ISMEAR  = 0
SIGMA   = 0.05
IBRION  = 2
ISIF    = 3
NSW     = 40
EDIFFG  = -0.01
LREAL   = .FALSE.
LWAVE   = .FALSE.
LCHARG  = .TRUE.
```

[custodian_vasp_run.py](../examples/06-custodian-vasp-si/custodian_vasp_run.py) 的核心部分：

```python
import os, shlex
from custodian.custodian import Custodian
from custodian.vasp.handlers import (
    VaspErrorHandler, MeshSymmetryErrorHandler, UnconvergedErrorHandler,
    NonConvergingErrorHandler, PotimErrorHandler, PositiveEnergyErrorHandler,
    FrozenJobErrorHandler, StdErrHandler, WalltimeHandler)
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator

vasp_cmd = shlex.split(os.environ.get("VASP_CMD", "srun vasp_std"))   # must be a list

handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(), UnconvergedErrorHandler(),
            NonConvergingErrorHandler(), PotimErrorHandler(), PositiveEnergyErrorHandler(),
            FrozenJobErrorHandler(timeout=7200), StdErrHandler()]
if "CUSTODIAN_WALLTIME" in os.environ:
    handlers.append(WalltimeHandler(wall_time=int(os.environ["CUSTODIAN_WALLTIME"]),
                                    buffer_time=600))

copy_contcar = {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}}
jobs = [
    VaspJob(vasp_cmd, final=False, suffix=".relax1"),
    VaspJob(vasp_cmd, final=False, suffix=".relax2", backup=False,
            settings_override=[copy_contcar]),
    VaspJob(vasp_cmd, final=True, suffix=".static", backup=False,
            settings_override=[copy_contcar,
                {"dict": "INCAR",
                 "action": {"_set": {"NSW": 0, "IBRION": -1, "ISIF": 2, "ALGO": "Normal"}}}]),
]
Custodian(handlers, jobs, validators=[VasprunXMLValidator(), VaspFilesValidator()],
          max_errors=10, max_errors_per_job=5,
          polling_time_step=10, monitor_freq=30).run()
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `vasp_cmd` | `["srun", "vasp_std"]` | **参数列表**，不是字符串；用 `shlex.split` 转换 |
| `final` | 最后一步 True | False：输出**复制**为 `OUTCAR.relax1` 等；True：输出**移动**为带后缀的文件 |
| `suffix` | `.relax1` / `.static` | 输出文件后缀 |
| `backup` | 第一步 True | 把原始输入保存为 `INCAR.orig` 等 |
| `settings_override` | 动作列表 | 运行前修改文件（`_file_copy`）或 INCAR 字典（`_set`/`_unset`） |
| `max_errors` | 5–10 | 整个序列允许的修正总次数（custodian 默认只有 1），超过后抛出 `MaxCorrectionsError` |
| `max_errors_per_job` | 3–5 | 单个 VaspJob 的修正上限（默认等于 `max_errors`） |
| `polling_time_step` | 10 s（默认） | 检查 VASP 是否仍在运行的间隔 |
| `monitor_freq` | 30（默认） | 每 30 次轮询（约 5 分钟）运行一次 monitor handler |
| `FrozenJobErrorHandler(timeout=)` | 7200–21600 s | 输出多久不更新视为卡死（默认 21600 s）；大体系杂化泛函要设大一些 |
| `WalltimeHandler(wall_time=, buffer_time=)` | 与 `--time` 一致 / 600 s | 到时前写 STOPCAR，避免被 Slurm 强杀后留下不完整的 vasprun.xml。`buffer_time` 默认 300 s；计时从 handler 创建时开始，所以 `wall_time` 应覆盖整个序列 |

> [!TIP]
> 为什么要弛豫两次？`ISIF = 3` 改变晶胞时平面波基组固定在初始晶胞上，第一次弛豫结束时的有效截断能与名义 ENCUT 不同（Pulay 应力）。从 CONTCAR 重新开始第二次弛豫可以消除这一偏差，最后的静态计算在固定晶胞上给出准确能量。原理见 [VASP 第4章](../../VASP/tutorials/04-relax-eos.md)。

### 4.4 YAML 写法与 `cstdn`

不想写 Python 时，可以用 [custodian_spec.yaml](../examples/06-custodian-vasp-si/custodian_spec.yaml) 描述同样的序列：

```yaml
jobs:
- jb: custodian.vasp.jobs.VaspJob
  params:
    final: False
    suffix: .relax
- jb: custodian.vasp.jobs.VaspJob
  params:
    final: True
    suffix: .static
    backup: False
    settings_override:
    - {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}}
    - {"dict": "INCAR", "action": {"_set": {"NSW": 0, "IBRION": -1}}}

jobs_common_params:
  $vasp_cmd: ["mpirun", "-np", "$SLURM_NTASKS", "vasp_std"]

handlers:
- hdlr: custodian.vasp.handlers.VaspErrorHandler
- hdlr: custodian.vasp.handlers.MeshSymmetryErrorHandler
- hdlr: custodian.vasp.handlers.UnconvergedErrorHandler
- hdlr: custodian.vasp.handlers.NonConvergingErrorHandler
- hdlr: custodian.vasp.handlers.PotimErrorHandler

validators:
- vldr: custodian.vasp.validators.VasprunXMLValidator

custodian_params:
  max_errors: 10
```

以 `$` 开头的键（`$vasp_cmd`）中的值会按环境变量展开，所以 `$SLURM_NTASKS` 在作业里会被替换为实际核数。运行命令为 `cstdn run custodian_spec.yaml`；`cstdn example` 会打印官方示例。

### 4.5 作业脚本与运行

[job_custodian_slurm.sh](../examples/06-custodian-vasp-si/job_custodian_slurm.sh) 的要点：custodian 运行在分配到的第一个节点上，由它启动 VASP。

```bash
#!/bin/bash
#SBATCH --job-name=si-custodian
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00
#SBATCH --partition=cpu

set -eo pipefail                             # no "-u": module/conda scripts use unset variables
ulimit -s unlimited
export OMP_NUM_THREADS=1
module purge
module load intel-oneapi vasp/6.5.1          # example module names
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wf
cd "$SLURM_SUBMIT_DIR"
if [ ! -f POTCAR ]; then
    : "${VASP_PP_PATH:?POTCAR missing and VASP_PP_PATH is not set}"
    cat "$VASP_PP_PATH/potpaw_PBE/Si/POTCAR" > POTCAR
fi
command -v dos2unix >/dev/null && dos2unix -q INCAR KPOINTS POSCAR

export VASP_CMD="srun vasp_std"
export CUSTODIAN_WALLTIME=7200               # seconds, same as --time
python custodian_vasp_run.py
```

```bash
bash -n job_custodian_slurm.sh      # syntax check
sbatch job_custodian_slurm.sh
```

作业脚本也可以用 [mkjob.py](../scripts/mkjob.py) 生成后再加上 custodian 这一行（`python ../scripts/mkjob.py -h` 查看用法）。

### 4.6 检查结果

```bash
ls *.relax1 *.relax2 *.static custodian.json error.*.tar.gz 2>/dev/null
grep -c "reached required accuracy" OUTCAR.relax1 OUTCAR.relax2
grep "free  energy   TOTEN" OUTCAR.static | tail -1
python show_custodian_log.py custodian.json          # readable list of corrections
```

[show_custodian_log.py](../examples/06-custodian-vasp-si/show_custodian_log.py) 只用标准库，可以在没装 custodian 的电脑上查看日志。**成功标准**：两个弛豫的 OUTCAR 各含一次 `reached required accuracy`；`custodian.json` 中最后一个 job 没有 `max_errors` 或 `validator` 失败标记；若出现过修正，逐条确认修正后的参数可以接受。

### 4.7 其他程序

custodian 还提供 `custodian.cp2k`（`UnconvergedScfErrorHandler`、`DivergingScfErrorHandler`、`FrozenJobErrorHandler`、`WalltimeHandler` 等）和 `custodian.gaussian`（`GaussianErrorHandler` 识别 `Convergence failure`、`Optimization stopped`、`Out-of-memory error in routine` 等，`GaussianJob` 负责运行）。atomate2 的 CP2K 工作流就是通过 `custodian.cp2k` 运行的。Gaussian 和 CP2K 的手动排错见 [Gaussian 第18章](../../Gaussian/tutorials/18-errors-capstone.md) 和 [CP2K 第18章](../../CP2K/tutorials/18-troubleshooting.md)。

---

## 5. jobflow：把计算写成依赖图

### 5.1 三个核心概念

| 概念 | 写法 | 含义 |
|---|---|---|
| Job | 用 `@job` 装饰函数，调用时**不执行**，只返回 `Job` 对象 | 一个计算步骤 |
| OutputReference | `job.output`、`job.output["energy"]` | 对未来输出的“引用”，运行时才解析成真实值 |
| Flow | `Flow([job1, job2], output=...)` | 一组 Job 及其依赖关系，可以嵌套 |

依赖关系不需要手工指定：只要 job B 的参数里出现了 job A 的 `output`，jobflow 就会自动加上 A → B 的边。

### 5.2 不装 DFT 也能跑的例子：EMT 状态方程

[jobflow_emt_eos.py](../examples/06-jobflow-atomate2/jobflow_emt_eos.py) 用 ASE 的 EMT 势做 7 个单点计算，再“扇入”到一个拟合作业，结构与 ENCUT 扫描、EOS 拟合完全相同：

```python
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.eos import EquationOfState
from ase.units import kJ
from jobflow import Flow, job, run_locally


@job
def energy_at(a: float, element: str = "Cu") -> dict:
    atoms = bulk(element, "fcc", a=a)
    atoms.calc = EMT()
    return {"a": a, "volume": atoms.get_volume(), "energy": atoms.get_potential_energy()}


@job
def fit_eos(points: list) -> dict:
    vols = np.array([p["volume"] for p in points])
    ens = np.array([p["energy"] for p in points])
    v0, e0, b = EquationOfState(vols, ens, eos="birchmurnaghan").fit()
    return {"a0": (4.0 * v0) ** (1 / 3), "v0": v0, "e0": e0, "B_GPa": b / kJ * 1.0e24}


scan = [energy_at(a) for a in np.linspace(0.95, 1.05, 7) * 3.6]
fit = fit_eos([j.output for j in scan])           # fan-in: 7 parents -> 1 child
flow = Flow([*scan, fit], output=fit.output, name="Cu EMT EOS")

responses = run_locally(flow, create_folders=True, ensure_success=True)
print(responses[fit.uuid][1].output)              # {uuid: {index: Response}}
```

```bash
pip install jobflow ase
python jobflow_emt_eos.py
python jobflow_emt_eos.py --folders     # one job_<timestamp>-<5 random digits> folder per job
```

| 参数（`run_locally`） | 推荐值 | 说明 |
|---|---|---|
| `create_folders` | True（真实计算） | 每个 job 在独立目录 `job_<UTC 时间戳>-<5 位随机数>` 中运行，避免文件互相覆盖（默认 False） |
| `root_dir` | 默认当前目录 | job 目录的父目录 |
| `ensure_success` | True | 有 job 失败时抛出异常，而不是静默返回 |
| `store` | 默认读 `~/.jobflow.yaml` | 没配置时用内存数据库（MemoryStore），程序结束结果即丢失 |

返回值是 `{uuid: {index: Response}}` 字典；`index` 通常为 1，job 被替换（replace）后会递增。**成功标准**：打印出 `a0` 与 `B`，没有异常。EMT 只用于验证流程，数值没有物理意义。

### 5.3 把结果存进 MongoDB

需要长期保存、按条件查询时，在 `~/.jobflow.yaml` 中配置 `JOB_STORE`（模板 [jobflow.yaml.example](../examples/06-jobflow-atomate2/jobflow.yaml.example)）：

```yaml
JOB_STORE:
  docs_store:
    type: MongoStore
    database: dft_tutorial
    collection_name: outputs
    host: localhost
    port: 27017
  additional_stores:
    data:
      type: GridFSStore
      database: dft_tutorial
      collection_name: outputs_blobs
      host: localhost
      port: 27017
```

（数据库开启认证时再加 `username`、`password`，完整模板见示例文件。）

`docs_store` 存小文档（能量、结构、参数），`additional_stores.data` 用 GridFS 存能带、DOS、电荷密度等大对象。也可以用环境变量 `JOBFLOW_CONFIG_FILE` 指向其他位置的配置文件。

---

## 6. atomate2：现成的 VASP 工作流

atomate2 把常见计算写成 **Maker**（带默认参数的数据类），`maker.make(structure)` 返回一个 jobflow `Flow`。每个 VASP 步骤在内部：用 pymatgen 输入集写 INCAR/KPOINTS/POSCAR/POTCAR → 调用 custodian 运行（默认 handler 包括 `VaspErrorHandler`、`UnconvergedErrorHandler`、`PotimErrorHandler`、`FrozenJobErrorHandler`、`LargeSigmaHandler`、`IncorrectSmearingHandler`、`KspacingMetalHandler` 等）→ 解析成 emmet 的 `TaskDoc`。

### 6.1 `~/.atomate2.yaml`

[atomate2.yaml.example](../examples/06-jobflow-atomate2/atomate2.yaml.example)：

```yaml
VASP_CMD: srun vasp_std
VASP_GAMMA_CMD: srun vasp_gam
VASP_NCL_CMD: srun vasp_ncl
VASP_CUSTODIAN_MAX_ERRORS: 5
VASP_HANDLE_UNSUCCESSFUL: error
VASP_INCAR_UPDATES:
  NCORE: 4
VASP_RUN_BADER: false
LOBSTER_CMD: lobster
CP2K_CMD: srun cp2k.psmp
```

| 键 | 默认值 | 说明 |
|---|---|---|
| `VASP_CMD` | `vasp_std` | 完整启动命令。先展开环境变量，再按 shell 规则拆分后交给 custodian，所以可以写 `srun -n $SLURM_NTASKS vasp_std` |
| `VASP_GAMMA_CMD` / `VASP_NCL_CMD` | `vasp_gam` / `vasp_ncl` | 只有 Γ 点时、非共线计算时使用 |
| `VASP_CUSTODIAN_MAX_ERRORS` | 5 | 每个 job 的 custodian 修正上限 |
| `VASP_HANDLE_UNSUCCESSFUL` | `error` | VASP 正常退出但未收敛时：`error` 报错；`true` 标记完成但停止后续 job；`false` 继续执行 |
| `VASP_INCAR_UPDATES` | `{}` | 加到**所有** INCAR 上的参数，适合放与集群相关的并行设置 |
| `CUSTODIAN_SCRATCH_DIR` | None | 在节点本地盘上运行 |
| `VASP_RUN_BADER` | false | 解析时调用 `bader`（需在 PATH 中） |

也可以用环境变量 `ATOMATE2_CONFIG_FILE` 指定配置文件，或用 `ATOMATE2_` 前缀的环境变量覆盖单个键。

### 6.2 弛豫 + 能带：`RelaxBandStructureMaker`

[atomate2_bandstructure_local.py](../examples/06-jobflow-atomate2/atomate2_bandstructure_local.py)：

```python
from jobflow import run_locally
from pymatgen.core import Lattice, Structure
from atomate2.vasp.flows.core import RelaxBandStructureMaker
from atomate2.vasp.powerups import add_metadata_to_flow, update_user_incar_settings

a = 5.43 / 2**0.5
si = Structure(Lattice.from_parameters(a, a, a, 60, 60, 60),
               ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])

flow = RelaxBandStructureMaker().make(si)
flow = update_user_incar_settings(flow, {"NCORE": 4, "LWAVE": False})
flow = add_metadata_to_flow(flow, {"project": "dft-tutorial", "system": "Si"})

for job, parents in flow.iterflow():          # inspect the DAG before running
    print(job.name, len(parents))

responses = run_locally(flow, create_folders=True, ensure_success=True)
```

工作流结构：`DoubleRelaxMaker`（两次弛豫，作业名带 ` 1`/` 2` 后缀）→ `static` → `non-scf uniform`（态密度）+ `non-scf line`（高对称线能带）。后面的步骤通过 `prev_dir` 读取前一步的 CHGCAR。

```bash
python atomate2_bandstructure_local.py --dry-run     # print the job graph only
python atomate2_bandstructure_local.py               # inside a Slurm job
```

常用 powerup（都返回修改后的新 Flow）：

| powerup | 作用 |
|---|---|
| `update_user_incar_settings(flow, {...})` | 修改所有 VASP job 的 INCAR |
| `update_user_kpoints_settings(flow, {...})` | 修改 k 点设置 |
| `update_user_potcar_settings` / `update_user_potcar_functional` | 更换赝势或泛函版本（如 `PBE_54`） |
| `add_metadata_to_flow(flow, {...})` | 给数据库文档加标签，方便查询 |
| `update_vasp_custodian_handlers(flow, handlers)` | 替换默认的 custodian handler（传空元组即关闭纠错） |

**检查**：每个 `job_*` 目录中有 gzip 压缩的输出（`OUTCAR.gz`、`vasprun.xml.gz`……）和 `custodian.json.gz`（默认 `VASP_ZIP_FILES: atomate`，custodian 日志也会被压缩）；`TaskDoc.output.bandgap`、`TaskDoc.output.energy` 给出带隙和能量；若 `VASP_HANDLE_UNSUCCESSFUL: error` 下有作业未收敛，`run_locally` 会因 `ensure_success=True` 抛出异常。能带作图见 [VASP 第5章](../../VASP/tutorials/05-band-dos.md)。

### 6.3 声子与 LOBSTER

[atomate2_more_flows.py](../examples/06-jobflow-atomate2/atomate2_more_flows.py) 默认只列出作业，加 `--run` 才真正运行：

```python
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.lobster import VaspLobsterMaker

phonon_flow = PhononMaker(
    min_length=15.0,        # supercell edges >= 15 A (default 20 A)
    displacement=0.01,      # finite displacement in A
    born_maker=None,        # Si is non-polar: skip the dielectric/Born-charge job
).make(si)

lobster_flow = VaspLobsterMaker().make(gaas)   # relax -> static (ISYM = 0, more bands) -> LOBSTER
```

| Maker 参数 | 默认 | 说明 |
|---|---|---|
| `PhononMaker.min_length` | 20.0 Å | 超胞最短边；越大越准、越贵 |
| `PhononMaker.displacement` | 0.01 Å | 有限位移大小 |
| `PhononMaker.born_maker` | `DielectricMaker` | 极性材料保留，用于非解析项修正；非极性材料设为 None 可省一步 |
| `PhononMaker.store_force_constants` | False（VASP 版） | 是否把力常数写进数据库 |
| `VaspLobsterMaker.delete_wavecars` | True | LOBSTER 完成后删除 WAVECAR 省空间 |

> [!NOTE]
> 声子工作流里的位移计算是**动态生成**的：`generate_phonon_displacements` 在运行时才根据弛豫后的结构生成位移超胞，`run_phonon_displacements` 再用 jobflow 的 `Response(replace=...)` 把自己替换成 N 个位移计算。所以运行前用 `iterflow()` 列出的只是顶层作业（弛豫、静态、超胞大小、生成位移、运行位移、后处理等），看不到具体有多少个位移计算。

```bash
pip install "atomate2[phonons,lobster]"
python atomate2_more_flows.py phonon            # list jobs
python atomate2_more_flows.py lobster --run     # run inside a Slurm job
```

声子和 COHP 的原理与手动做法见 [VASP 第12章](../../VASP/tutorials/12-phonon.md) 和 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)。

> [!WARNING]
> `run_locally` 在**当前进程**里依次运行所有 job：如果一个声子工作流有 20 个位移计算，它们会在同一个 Slurm 作业中串行执行。这样做适合学习和小体系；真实高通量请交给 jobflow-remote（第7节），让每个 job 成为独立的 Slurm 作业并行排队。

---

## 7. jobflow-remote 1.x：提交到集群并用 MongoDB 管理

### 7.1 架构

jobflow-remote 涉及三类机器（可以是同一台）：

```text
 USER    : 创建 Flow 并写入数据库 (submit_flow)、查看状态 (jf job list)
 RUNNER  : 守护进程 (jf runner start)，负责上传输入、提交、轮询、下载输出、入库
 WORKER  : 计算集群，作业在这里真正运行（只需装 jobflow-remote，不需要连数据库）
 MongoDB : 队列库（作业状态）+ 结果库（JobStore）；USER 和 RUNNER 必须能连上
```

最常见的两种部署方式：（1）**全在一处**：数据库和 runner 都在集群登录节点上可用（需管理员允许常驻进程）；（2）**工作站方式**：数据库和 runner 在课题组工作站上，通过免密 SSH 连接集群。

### 7.2 项目文件

放在 `~/.jfremote/<项目名>.yaml`。[jfremote_project.yaml](../examples/06-jobflow-atomate2/jfremote_project.yaml) 的主要部分：

```yaml
name: dft_tutorial

workers:
  hpc_cpu:
    type: remote
    scheduler_type: slurm
    host: login.cluster.example.edu
    user: myname
    work_dir: /scratch/myname/jfr_runs     # absolute path, no $VARIABLES
    pre_run: |
      source ~/miniforge3/etc/profile.d/conda.sh
      conda activate wf
      module load intel-oneapi vasp/6.5.1
      ulimit -s unlimited
      export OMP_NUM_THREADS=1
    timeout_execute: 60
    max_jobs: 20
    resources:
      partition: cpu
      nodes: 1
      ntasks_per_node: 32
      cpus_per_task: 1
      time: "08:00:00"
  local_shell:
    type: local
    scheduler_type: shell
    work_dir: /home/myname/jfr_local
    timeout_execute: 60

queue:
  store:
    type: MongoStore
    host: localhost
    port: 27017
    database: jfr_dft_tutorial
    username: dft_user
    password: change_me
    collection_name: jobs

jobstore:
  docs_store:
    type: MongoStore
    host: localhost
    port: 27017
    database: jfr_dft_tutorial
    username: dft_user
    password: change_me
    collection_name: outputs

exec_config:
  vasp6:
    modules: [intel-oneapi, vasp/6.5.1]
    export: {OMP_NUM_THREADS: 1}
    pre_run: ulimit -s unlimited
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `workers.<名>.type` | `remote` / `local` | remote 通过 SSH（Fabric）连接，需免密登录（[第1章](01-linux-hpc.md)第4节） |
| `scheduler_type` | `slurm` / `pbs` / `shell` | 由 qtoolkit 生成作业脚本；`shell` 不排队直接运行 |
| `work_dir` | 集群上的绝对路径 | 每个 job 在其下建一个子目录 |
| `resources` | 分区、节点、核数、时限 | 键名对应 qtoolkit 的 Slurm 模板（`partition`、`nodes`、`ntasks_per_node`、`cpus_per_task`、`time`、`account`、`qos`、`mem` 等），提交时可按 job 覆盖 |
| `max_jobs` | 10–50 | 队列中同时存在的作业上限，避免超过集群配额 |
| `pre_run` | 激活环境、加载模块 | 写进每个作业脚本 |
| `queue.store` | `MongoStore` | **必须**是 MongoStore 及其子类；不同项目不要共用同一个队列库 |
| `jobstore` | 与 jobflow 的 `JOB_STORE` 同格式 | jobflow-remote 会**忽略** `~/.jobflow.yaml` |
| `exec_config` | 按程序命名 | 提交时用 `exec_config="vasp6"` 选择 |

`jf project generate --full <名字>` 可以生成带全部选项和注释的模板。

### 7.3 初始化、提交与监控

```bash
pip install jobflow-remote
mkdir -p ~/.jfremote && cp jfremote_project.yaml ~/.jfremote/dft_tutorial.yaml
chmod 600 ~/.jfremote/dft_tutorial.yaml        # contains database passwords
jf project check --errors                      # tests worker SSH, queue store, jobstore
jf admin reset                                 # first time only: wipes the queue database!
jf runner start
python submit_jfremote.py POSCAR               # adds the flow to the queue
jf job list                                    # states of all jobs
jf flow list
jf job info 3                                  # details / error message of job 3
jf job output 3                                # output document of a completed job
```

提交脚本 [submit_jfremote.py](../examples/06-jobflow-atomate2/submit_jfremote.py) 的关键一行：

```python
from jobflow_remote import submit_flow
ids = submit_flow(flow, worker="hpc_cpu", exec_config="vasp6")
```

`submit_flow` 只把作业写进数据库，**不运行**，真正的执行由 runner 完成。作业状态的正常顺序为：

```text
WAITING -> READY -> CHECKED_OUT -> UPLOADED -> SUBMITTED -> RUNNING
        -> RUN_FINISHED -> DOWNLOADED -> COMPLETED
```

异常状态：`FAILED`（job 函数本身出错，如 custodian 达到修正上限）、`REMOTE_ERROR`（上传、提交、下载等环节失败，如网络中断）、`PAUSED`、`STOPPED`、`USER_STOPPED`。

| 情况 | 命令 |
|---|---|
| `REMOTE_ERROR`（网络、SSH 临时故障） | `jf job retry <id>`；批量：`jf job retry -s REMOTE_ERROR` |
| `FAILED`（计算出错，已修改参数） | `jf job rerun <id>`（作业回到 READY，子作业回到 WAITING）；`retry` 只用于 `REMOTE_ERROR` 或重置远程操作的尝试次数 |
| 修改单个作业的资源/执行配置 | `jf job set resources ...`、`jf job set exec-config ...`、`jf job set worker ...` |
| 作业被锁住（runner 异常退出） | `jf job list -l`（`--locked`）查看，确认 runner 已停止后 `jf admin unlock -jid <id>` |
| 修改了项目文件 | `jf runner stop` 后 `jf runner start`，或 `jf runner restart`（runner 只在启动时读取配置） |
| 查看 runner 状态 | `jf runner status`、`jf runner info` |

> [!CAUTION]
> jobflow-remote 1.x 的 `jf job rerun` 默认会**删除** WORKER 上该作业（及状态被重置的子作业）的运行目录，再从 READY 重新开始。如果加了 `--no-delete`（`-nd`），或者你手动在原目录里重跑，旧的 CONTCAR、WAVECAR、`custodian.json` 会留在**同一目录**里，可能被误读。这种情况下重跑前先清理或移走旧的 VASP 输出。另外要确认 WORKER 上 `~/.atomate2.yaml` 中的 `VASP_CMD` 正确。

---

## 8. AiiDA：带完整溯源的工作流

### 8.1 核心思想

AiiDA 把**每一个**输入（结构、参数、赝势、代码）、计算过程和输出都存成数据库节点，并用边记录“谁由谁生成”，形成溯源图（provenance graph）。即使几年后，也能从一个能量值追溯到它的完整计算链。代价是需要维护数据库（SQLite 或 PostgreSQL）和消息中间件（broker），学习曲线也更陡。

| 概念 | 说明 |
|---|---|
| Profile | 一套独立的数据库 + 文件仓库 + broker 配置 |
| Computer | 运行计算的机器（localhost 或通过 SSH 连接的集群）+ 调度器类型 |
| Code | 某台 Computer 上的一个可执行程序，绑定一个 CalcJob 插件 |
| CalcJob | 一次程序运行（如 `Cp2kCalculation`），自动生成输入、提交、取回、解析 |
| WorkChain | 由多个 CalcJob/WorkChain 组成的工作流，可带自动纠错的 process handler |
| daemon | 后台进程，执行 `submit()` 提交的工作；`run()` 则在当前进程中阻塞执行 |

### 8.2 `verdi presto` 快速建库

```bash
conda activate aiida
verdi presto                       # SQLite profile + localhost computer, no questions asked
verdi status                       # check profile, storage, broker, daemon
verdi presto --use-postgres        # production: PostgreSQL storage (server must be running)
```

`verdi presto` 从 aiida-core 2.6 开始提供。在 2.9.4 中，它会先检测本机是否运行 RabbitMQ，检测不到时改用 ZeroMQ broker（随 daemon 自动启动，`--use-zeromq` 可跳过检测）；`--no-broker` 则不配置 broker，此时只能用 `verdi run` 阻塞运行，不能用 daemon 提交。更早的版本在没有 RabbitMQ 时会建立无 broker 的 profile。

### 8.3 注册 Computer 和 Code

示例 YAML 位于 [06-aiida-cp2k-vasp](../examples/06-aiida-cp2k-vasp/README.md)：

```bash
# codes on localhost (created by verdi presto)
verdi code create core.code.installed --config cp2k_code.yaml
verdi code create core.code.installed --config vasp_code.yaml
verdi code create core.code.installed --config gaussian_code.yaml

# a remote Slurm cluster
verdi computer setup --config computer_slurm.yaml
verdi computer configure core.ssh hpc --config computer_ssh.yaml
verdi computer test hpc
```

`cp2k_code.yaml`：

```yaml
label: cp2k
description: CP2K 2026.2 psmp (example)
default_calc_job_plugin: cp2k
computer: localhost
filepath_executable: /opt/cp2k/exe/local/cp2k.psmp
with_mpi: true
prepend_text: |
  export OMP_NUM_THREADS=1
  ulimit -s unlimited
append_text: ''
```

| 关键词 | 说明 |
|---|---|
| `default_calc_job_plugin` | `cp2k`、`vasp.vasp`、`gaussian` 等插件入口名 |
| `filepath_executable` | 程序在**该 Computer 上**的绝对路径 |
| `with_mpi` | CP2K/VASP 为 true；Gaussian 为共享内存并行，设为 false |
| `prepend_text` | 写进作业脚本、在程序启动前执行（module、环境变量、ulimit） |
| `mpirun_command`（Computer） | 如 `srun` 或 `mpirun -np {tot_num_mpiprocs}` |
| `work_dir`（Computer） | 远程运行目录，`{username}` 会被替换 |

### 8.4 aiida-cp2k：参数写成嵌套字典

[aiida_cp2k_sp.py](../examples/06-aiida-cp2k-vasp/aiida_cp2k_sp.py) 计算 H2O 的 PBE 单点能。CP2K 输入的 section 对应字典的键，section 参数（如 `&XC_FUNCTIONAL PBE` 中的 `PBE`）用键 `"_"` 表示：

```python
parameters = orm.Dict({
    "GLOBAL": {"RUN_TYPE": "ENERGY", "PRINT_LEVEL": "MEDIUM"},
    "FORCE_EVAL": {
        "METHOD": "Quickstep",
        "DFT": {
            "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
            "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
            "QS": {"EPS_DEFAULT": 1.0e-12},
            "MGRID": {"CUTOFF": 400, "REL_CUTOFF": 50},
            "POISSON": {"PERIODIC": "NONE", "PSOLVER": "MT"},
            "SCF": {"SCF_GUESS": "ATOMIC", "EPS_SCF": 1.0e-6, "MAX_SCF": 50,
                    "OT": {"MINIMIZER": "DIIS", "PRECONDITIONER": "FULL_ALL"}},
            "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
        },
        "SUBSYS": {
            "CELL": {"PERIODIC": "NONE"},
            "KIND": [
                {"_": "O", "BASIS_SET": "DZVP-MOLOPT-SR-GTH", "POTENTIAL": "GTH-PBE-q6"},
                {"_": "H", "BASIS_SET": "DZVP-MOLOPT-SR-GTH", "POTENTIAL": "GTH-PBE-q1"},
            ],
        },
    },
})
builder = CalculationFactory("cp2k").get_builder()
builder.code = orm.load_code("cp2k@localhost")
builder.structure = orm.StructureData(ase=atoms)    # coordinates + cell are written by the plugin
builder.parameters = parameters
builder.metadata.options.resources = {"num_machines": 1, "num_mpiprocs_per_machine": 4}
results, node = run_get_node(builder)
print(results["output_parameters"]["energy"])      # Hartree
```

插件会生成 `aiida.inp`、坐标文件 `aiida.coords.xyz`，并根据 `StructureData` 写出 `CELL` 的 A/B/C，所以字典中不要再写这几项。CP2K 输入各 section 的含义见 [CP2K 第2章](../../CP2K/tutorials/02-input-structure.md)。

> [!WARNING]
> CP2K 2026.2 改动了部分输入语法：删除了 `&MOTION/&CELL_OPT/TYPE`，`&DFT/&PRINT/&PDOS` 改为 `&DFT/&PRINT/&DOS`（内含 `&PDOS`/`&CURVE`/`&BROADEN`），`&REFTRAJ` 的 `EVAL_ENERGY_FORCES`/`EVAL_FORCES` 改为 `EVAL`。aiida-cp2k 的输入字典由你自己写，旧教程里的字典要按新语法改；atomate2 通过 pymatgen 生成 CP2K 输入，而 pymatgen-core 2026.10.31 的 `activate_dos()`/`activate_pdos()` 仍在 `FORCE_EVAL/DFT/PRINT/PDOS` 下写 PDOS。在 2026.2 上批量运行前，先生成一个输入，确认其中没有已删除的写法。旧输入保留 `&CELL_OPT/TYPE` 时，CP2K 会报 `The specified keyword 'TYPE' is not available anymore`。

### 8.5 aiida-vasp：协议驱动的 WorkChain

aiida-vasp 5.x 提供 `vasp.vasp`（单次计算 + 纠错重启）、`vasp.relax`（反复弛豫直到体积收敛，最后做静态计算）、`vasp.bands`、`vasp.converge`、`vasp.neb` 等 WorkChain。先上传一次赝势族：

```bash
aiida-vasp potcar uploadfamily --path /path/to/potpaw_PBE.54 --name PBE.54 --description "PBE.54"
aiida-vasp potcar listfamilies
```

[aiida_vasp_relax.py](../examples/06-aiida-cp2k-vasp/aiida_vasp_relax.py)：

```python
VaspRelaxWorkChain = WorkflowFactory("vasp.relax")
builder = VaspRelaxWorkChain.get_builder_from_protocol(
    code=orm.load_code("vasp_std@localhost"),
    structure=orm.StructureData(ase=bulk("Si", "diamond", a=5.43)),
    protocol="balanced",
    overrides={"vasp": {"parameters": {"incar": {"encut": 400, "ispin": 1, "ncore": 4,
                                               "gga": "pe"}}}},   # protocol default is PBEsol
    relax_settings={"force_cutoff": 0.01, "steps": 40},
    options={"resources": {"num_machines": 1, "num_mpiprocs_per_machine": 32},
             "max_wallclock_seconds": 3600, "queue_name": "cpu"},
)
results, node = run_get_node(builder)
relaxed = results["relax"]["structure"]
```

| 参数 | 说明 |
|---|---|
| `protocol` | aiida-vasp 5.1.0 的 `vasp.yaml` 只定义了 `balanced`（默认）：`PBE.54` 赝势族、`PREC = Accurate`、`ENCUT = 550`、`kpoints_spacing = 0.05`、`ISPIN = 2`、`ISMEAR = 0`/`SIGMA = 0.05`，并且 **`GGA = PS`（PBEsol）**。要用 PBE，必须在 `overrides` 中写 `"gga": "pe"`。`relax.yaml` 另有一个只放宽弛豫判据的 `fast`（`force_cutoff = 0.05`、`convergence_volume = 0.05`），但内部的 VaspWorkChain 也会去 `vasp.yaml` 里找同名协议，找不到会报 `ValueError`。需要宽松判据时，保持 `balanced`，直接用 `relax_settings` 设置 `force_cutoff`、`convergence_volume` |
| `overrides` | 按输入命名空间嵌套，合并到协议上；INCAR 键用小写 |
| `relax_settings` | `force_cutoff`（eV/Å，默认 0.03）、`steps`（默认 60）、`convergence_volume`（体积相对变化阈值，默认 0.01）、`convergence_max_iterations`（默认 5）、`shape`/`volume`/`positions` 开关等 |
| `options` | 递归设置到所有内部 CalcJob 的 `metadata.options` |

协议默认值会随版本变化。实际使用时用 `verdi process show <PK>` 查看最终写出的 INCAR 参数，不要凭记忆。

### 8.6 aiida-gaussian：SCF 失败自动换 QC

[aiida_gaussian_base.py](../examples/06-aiida-cp2k-vasp/aiida_gaussian_base.py) 使用 `GaussianBaseWorkChain`（入口名 `gaussian.base`）：

```python
parameters = orm.Dict({
    "link0_parameters": {"%chk": "aiida.chk", "%mem": "4000MB", "%nprocshared": 8},
    "functional": "B3LYP",
    "basis_set": "6-31G(d)",
    "charge": 0,
    "multiplicity": 1,
    "route_parameters": {"opt": None, "scf": {"conver": 8, "maxcycle": 128}},
})
builder = WorkflowFactory("gaussian.base").get_builder()
builder.gaussian.code = orm.load_code("g16@localhost")
builder.gaussian.structure = orm.StructureData(ase=ase.build.molecule("H2O"))
builder.gaussian.parameters = parameters
builder.gaussian.metadata.options.resources = {"num_machines": 1, "tot_num_mpiprocs": 8}
builder.gaussian.metadata.options.max_memory_kb = int(1.25 * 4000) * 1024
```

内置 process handler：SCF 不收敛时先改 `SCF=YQC` 重算，仍失败再改 `SCF=XQC`，两者都失败则放弃（重写 `scf` 选项时只保留 `conver`，`maxcycle` 等会被丢掉）；默认关闭的 `handle_asytop_error`（用 `builder.handler_overrides = {"handle_asytop_error": True}` 打开）在出现 ASyTop 错误时把坐标四舍五入到小数点后 4 位再重算；正常结束标志缺失（如超时）的作业不会自动重启。手动处理 SCF 不收敛的方法见 [Gaussian 第4章](../../Gaussian/tutorials/04-sp-scf.md)。输出由 cclib 解析，`output_parameters["scfenergies"]` 的单位是 eV。

### 8.7 运行与检查

```bash
verdi run aiida_cp2k_sp.py                       # blocking
verdi daemon start 2                             # 2 daemon workers
python aiida_vasp_relax.py --submit              # through the daemon
verdi process list -a -p 1                       # all processes of the last day
verdi process report <PK>                        # WorkChain log: restarts and handler actions
verdi process show <PK>
verdi calcjob outputcat <PK>                     # main output file of a CalcJob
verdi node graph generate <PK>                   # provenance graph (PDF)
```

**成功标准**：`verdi process list` 中显示 `Finished [0]`。非零退出码表示进程结束但计算未成功，用 `verdi process report` 查看是哪一步、哪个 handler 放弃了。

> [!CAUTION]
> AiiDA 会把 POTCAR 的内容存进数据库仓库。用 `verdi archive create` 导出并公开分享（如上传到 Materials Cloud）之前，必须确认导出内容中不包含受许可证保护的赝势文件。

---

## 9. 轻量方案：quacc 配方与 dpdispatcher

### 9.1 quacc：ASE 计算器 + 配方

quacc 的“配方”就是带 `@job`/`@flow` 装饰器的 Python 函数，内部调用 ASE 计算器，返回普通字典（`atoms`、`results`、`dir_name` 等）。不设置工作流引擎（默认 `WORKFLOW_ENGINE = None`）时，调用函数就立即执行，不需要数据库；设置为 `jobflow`、`parsl`、`dask`、`prefect` 等之后，同一段代码会交给对应的引擎调度。

[quacc_emt_demo.py](../examples/06-quacc-dpdispatcher/quacc_emt_demo.py)：

```python
from ase.build import bulk, fcc111
from quacc.recipes.emt.core import relax_job, static_job

bulk_res = relax_job(bulk("Cu", "fcc", a=3.6), relax_cell=True)
e_bulk = bulk_res["results"]["energy"] / len(bulk_res["atoms"])

a_relaxed = bulk_res["atoms"].cell.cellpar()[0] * 2**0.5      # primitive -> cubic a
slab_res = relax_job(fcc111("Cu", size=(2, 2, 4), vacuum=8.0, a=a_relaxed),
                     opt_params={"fmax": 0.02})
single = static_job(slab_res["atoms"])
gamma = (single["results"]["energy"] - len(single["atoms"]) * e_bulk) \
        / (2 * slab_res["atoms"].cell.area(2))
```

把 import 换成 `quacc.recipes.vasp.core`（`relax_job(atoms, preset="DefaultSetGGA", relax_cell=...)`）或 `quacc.recipes.gaussian.core`（`relax_job(atoms, charge, spin_multiplicity, xc=..., basis=...)`），同样的调用方式就会运行 VASP 或 Gaussian。设置写在 `~/.quacc.yaml` 中（[quacc.yaml.example](../examples/06-quacc-dpdispatcher/quacc.yaml.example)）：

| 键 | 说明 |
|---|---|
| `WORKFLOW_ENGINE` | `null`（立即执行）/ `jobflow` / `parsl` / `dask` / `prefect` / `ray` / `redun`；只能通过配置文件或 `QUACC_WORKFLOW_ENGINE` 设置 |
| `RESULTS_DIR` / `SCRATCH_DIR` | 结果目录 / 临时运行目录 |
| `VASP_PARALLEL_CMD` + `VASP_CMD` | 拼接成 `"<VASP_PARALLEL_CMD> <VASP_CMD>"`，如 `srun` + `vasp_std` |
| `VASP_USE_CUSTODIAN` | 默认 true，VASP 经由 custodian 运行 |
| `VASP_PP_PATH` | 含 `potpaw_PBE` 的上级目录（与 ASE 约定相同） |
| `GAUSSIAN_CMD` | 默认 `g16` |

quacc 的优点是**一套 ASE 接口覆盖多个程序**，并且引擎可以替换。CP2K 和 OpenMX 目前没有 quacc 配方，可以用 ASE 计算器（[第4章第6节](04-ase-pymatgen.md)）自己写一个 `@job` 函数。

### 9.2 dpdispatcher：只管提交

dpdispatcher 不定义工作流，只负责把一批目录变成作业、提交、轮询、取回文件，适合收敛测试、机器学习势数据标注等“同一命令跑很多目录”的任务。基本用法已在 [第2章第11节](02-job-scheduler.md) 介绍；本章的 [dpdisp_encut_scan.py](../examples/06-quacc-dpdispatcher/dpdisp_encut_scan.py) 从自己的电脑通过 SSH 向集群提交 ENCUT 扫描：

```json
{
    "batch_type": "Slurm",
    "context_type": "SSHContext",
    "local_root": "./",
    "remote_root": "/scratch/myname/dpdisp_work",
    "remote_profile": {
        "hostname": "login.cluster.example.edu",
        "username": "myname",
        "port": 22,
        "key_filename": "/home/myname/.ssh/id_ed25519",
        "timeout": 30
    },
    "retry_count": 2
}
```

```bash
python dpdisp_encut_scan.py --template ../06-custodian-vasp-si --dry-run   # upload + scripts only
python dpdisp_encut_scan.py --template ../06-custodian-vasp-si             # submit, wait, collect E0
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `context_type` | `SSHContext` / `LocalContext` | 从本地电脑提交 / 在集群上直接提交 |
| `remote_root` | 集群工作盘绝对路径 | 任务目录上传到这里 |
| `remote_profile.key_filename` | 私钥的**绝对路径** | dpdispatcher 1.1.0 只对它做 `os.path.abspath`，不展开 `~`，写成 `~/.ssh/...` 会报 `... not found!` |
| `retry_count` | 1–3（默认 3） | 作业失败后自动重新提交的次数（**原样重提**，不修改输入），用完后报错 |
| `group_size`（resources） | 1 或更大 | 每个 Slurm 脚本打包几个任务 |
| `backward_files`（Task） | OUTCAR、OSZICAR、vasp.out | 只取回小文件 |

> [!TIP]
> dpdispatcher 原样重提不会修正输入。需要纠错时，把 `Task(command=...)` 换成调用 custodian 的脚本（第4节），两者可以组合使用：dpdispatcher 负责批量提交，custodian 负责每个作业内部的纠错。

---

## 10. 选型建议

| 场景 | 推荐组合 | 理由 | 不推荐 |
|---|---|---|---|
| **个人课题**：几十个结构、固定的 2–3 步 | Slurm 依赖链/作业数组（[第2章](02-job-scheduler.md)）+ **custodian**；或 **quacc**（无引擎模式）；批量用 **dpdispatcher** | 不需要数据库，学习成本最低；custodian 处理大部分数值报错 | 一开始就装 MongoDB/PostgreSQL |
| **课题组高通量**：几百到几万个计算，多人共用 | **atomate2 + jobflow-remote + MongoDB**（VASP 为主）；或 **quacc + jobflow/Parsl** | 现成 Maker 覆盖弛豫、能带、声子、弹性、LOBSTER；数据库可按成分、性质查询；runner 自动管理排队上限 | 每人各写一套 shell 脚本 |
| **需要完整溯源**：要发布数据库、要求审计、跨程序复杂工作流 | **AiiDA**（PostgreSQL + RabbitMQ）+ aiida-vasp / aiida-cp2k / aiida-gaussian | 每个数字都能追溯到输入和代码版本；可导出到 Materials Cloud | 用 SQLite profile 跑生产 |

补充原则：

1. **以程序为主选择**：VASP 为主选 atomate2；CP2K 为主两者都可以（atomate2 的 CP2K 支持与 aiida-cp2k）；Gaussian 为主可以用 aiida-gaussian、quacc 的 Gaussian 配方或 `custodian.gaussian`；OpenMX 没有专门的工作流插件，用 dpdispatcher + 自己的脚本，或者用 jobflow `@job` 包装 ASE 的 OpenMX 计算器。
2. **以团队为主选择**：有人维护数据库服务器就选 jobflow-remote 或 AiiDA；没有就选 custodian + dpdispatcher。
3. **先小后大**：任何框架都先用 1 个结构跑通、检查输入文件，再扩到 10 个，最后才是全部。
4. **数据是长期资产**：无论选什么，都要保存原始输入、程序版本、`custodian.json`/溯源图，并用 [job_status.py](../scripts/job_status.py) 或数据库查询定期统计失败率。

---

## 参数选择与收敛建议

| 项目 | 建议 | 原因 |
|---|---|---|
| custodian `max_errors` | 5–10（每个 job 3–5） | 太小容易放弃可修复的作业；太大会在根本错误（结构不合理、赝势错误）上浪费机时 |
| `FrozenJobErrorHandler.timeout` | 至少为单个离子步耗时的 3 倍 | 太短会误杀正常但慢的杂化泛函计算 |
| `WalltimeHandler` | `buffer_time` ≥ 1 个离子步耗时 | 留出时间写完 CONTCAR、WAVECAR |
| 修正后的参数 | 弛豫完成后用统一参数重做静态计算 | 不同结构可能被改成不同的 ALGO/ISMEAR，能量比较前要统一 |
| atomate2 输入集 | 一个项目只用一个输入集，并用 powerup 统一修改 | 混用 MP 旧集、`MPScanRelaxSet`、`MatPESRelaxSet` 得到的能量不可直接比较 |
| `VASP_INCAR_UPDATES` | 只放并行参数（NCORE、KPAR） | 物理参数应写在 Maker 或 powerup 里，便于从文档中追溯 |
| jobflow-remote `max_jobs` | 不超过集群单用户排队上限的 70–80% | 留余量给手动作业，避免被管理员限制 |
| 收敛测试 | 先用本章工具（或[第2章](02-job-scheduler.md)第7节的作业数组）自动化 [VASP 第3章](../../VASP/tutorials/03-convergence.md)的 ENCUT/k 点扫描 | 高通量的参数必须对最难收敛的体系也足够 |
| 数据库备份 | MongoDB 用 `mongodump`，AiiDA 用 `verdi storage backup` | 数据库损坏后无法恢复溯源信息 |
| 版本 | 项目开始时锁定版本，升级前重跑基准 | 输入集默认值和 handler 行为会随版本变化 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `MaxCorrectionsError: Max errors reached: 10.` | custodian 修正总次数已达上限 | 解开最后一个 `error.N.tar.gz` 查看现场；多为结构不合理、磁矩初值或赝势问题，需要人工处理 |
| `MaxCorrectionsPerJobError: Max errors per job reached: 5.` | 同一个 VaspJob 反复出错 | 同上；查看 `custodian.json` 是否在两种修正之间来回切换 |
| `NonRecoverableError: Unrecoverable error for handler: ...` | handler 识别了错误但没有可用的修正 | 按报错类型手动修改（参考 [VASP 第18章](../../VASP/tutorials/18-postprocess-troubleshooting.md)） |
| `ValidationError: Validation failed: VasprunXMLValidator` | VASP 结束但 `vasprun.xml` 不完整或无法解析 | 检查是否被 Slurm 强杀、磁盘是否已满；设置 `WalltimeHandler` |
| `FileNotFoundError: [Errno 2] No such file or directory: 's'`，或 VASP 立即退出、`vasp.out` 为空 | `vasp_cmd` 写成了字符串（VaspJob 会把它 `tuple()` 成单个字符），或 `srun`/模块在作业环境中不可用 | 用 `shlex.split()` 转成列表；在作业脚本里 `module load`，并确认 `which vasp_std` |
| `VaspErrorHandler` 从不触发 | 用户自己把 VASP 输出重定向到别的文件 | 让 VaspJob 写默认的 `vasp.out`，或用 `output_filename=` 告诉 handler |
| `FileNotFoundError: You do not have the right POTCAR with functional='PBE' and symbol='Fe_pv' ...`（旧版 pymatgen 为 `No POTCAR for Fe_pv with functional PBE found ...`） | pymatgen 未配置赝势目录或缺少该版本 | `pmg config -p ...` + `pmg config --add PMG_VASP_PSP_DIR ...`；需要 PBE_54 时设置 `PMG_DEFAULT_FUNCTIONAL` |
| `ImportError: cannot import name 'Structure' from 'pymatgen'` | 旧代码写法 | `from pymatgen.core import Structure` |
| 在 Python 3.11 环境中 `pip install atomate2` 装上了旧版本，导入 Maker 或参数时报错 | atomate2 0.1.x 要求 Python ≥ 3.12，pip 会跳过不兼容的版本，退回到旧版 | `pip show atomate2` 查看版本；新建 Python 3.12 环境 |
| atomate2 作业全在一个 Slurm 作业里串行跑、很慢 | `run_locally` 的设计如此 | 改用 jobflow-remote，每个 job 独立排队 |
| `run_locally(..., ensure_success=True)` 抛出 `RuntimeError: Flow did not finish running successfully` | 某个 job 失败（常见：VASP 未收敛且 `VASP_HANDLE_UNSUCCESSFUL: error`） | 在对应 `job_*` 目录查看 `custodian.json` 和 OUTCAR |
| `jf project check --errors` 中 worker 一项失败 | SSH 不能免密登录、`host`/`user` 写错、需要跳板机 | 先确认 `ssh <host>` 不需要密码；跳板机写进 `~/.ssh/config` 或用 `gateway` |
| 作业停在 `REMOTE_ERROR` | 上传/提交/下载时网络或文件系统出错 | `jf job info <id>` 查看原因，修复后 `jf job retry <id>` |
| 修改项目文件后不生效 | runner 只在启动时读取配置 | `jf runner stop && jf runner start` |
| 重跑后读到旧结果 | 用 `jf job rerun --no-delete` 或手动在原目录重跑，旧 CONTCAR/WAVECAR 仍在 | 重跑前清理旧 VASP 文件，或使用默认的 `jf job rerun`（会删除旧运行目录，见 7.3 节） |
| `verdi status` 显示 broker 无法连接 RabbitMQ | 没装或没启动 RabbitMQ | aiida ≥ 2.6 先用 `verdi presto` 跑通（2.9 会自动改用 ZeroMQ）；生产环境再安装并启动 RabbitMQ |
| `verdi process list` 显示 `Finished [非零]` | WorkChain 的纠错用尽或遇到无法处理的错误 | `verdi process report <PK>` 看日志，`verdi calcjob outputcat` 看原始输出 |
| `verdi computer test hpc` 失败 | SSH 配置、`work_dir` 无写权限、调度器命令不可用 | 逐项看测试输出；确认 `work_dir` 已存在且可写 |
| `paramiko.ssh_exception.AuthenticationException: Authentication failed.`（dpdispatcher/AiiDA） | 私钥路径或用户名错误、密钥权限过宽 | `chmod 600 ~/.ssh/id_ed25519`，确认 `key_filename` 和 `username` |
| dpdispatcher 报 `OSError: /当前目录/~/.ssh/id_ed25519 not found!` | `key_filename` 写成了 `~/...`，dpdispatcher 不展开 `~` | 改成绝对路径，如 `/home/myname/.ssh/id_ed25519` |
| `sbatch: error: Batch script contains DOS line breaks (\r\n)` | 作业脚本或 YAML 在 Windows 下编辑过 | `dos2unix` 后再提交（[第1章](01-linux-hpc.md)第10节） |
| CP2K 2026.2 报 `The specified keyword 'TYPE' is not available anymore` | aiida-cp2k 的输入字典或自己的模板中仍有已删除的 `&MOTION/&CELL_OPT/TYPE` | 删除该行即可（2026.2 的晶胞优化总是直接法，同时更新坐标和晶胞），见 [CP2K 第7章](../../CP2K/tutorials/07-geo-cell-opt.md) |

---

## 练习

1. 在 [06-custodian-vasp-si](../examples/06-custodian-vasp-si/README.md) 中把 INCAR 改为 `ALGO = VeryFast`、`NELM = 6`，提交后用 `show_custodian_log.py` 观察 `UnconvergedErrorHandler`（和 `NonConvergingErrorHandler`）是如何逐级修改 `ALGO` 和混合参数的。对照 4.2 节的表格，写出每一步修正的理由。由于 handler 不会改 NELM，这个作业最终会以 `MaxCorrectionsPerJobError` 或 `NonRecoverableError` 结束，请说明为什么“自动纠错”救不了这种输入错误。
2. 修改 `jobflow_emt_eos.py`：增加一个 `@job`，把 EOS 拟合得到的 `a0` 传给一个新的单点计算，验证能量是否为扫描中的最低值。画出这个 Flow 的 DAG（父子关系）。
3. 用 `atomate2_more_flows.py phonon` 列出 Si 声子工作流的顶层作业，指出哪一个会在运行时展开成位移计算（6.3 节的 NOTE）。再用 phonopy 的 Python API（`Phonopy(unitcell, supercell_matrix=...)`、`generate_displacements(distance=0.01)`、`supercells_with_displacements`）对 Si 原胞分别取 3×3×3 和 5×5×5 超胞（边长约 11.5 Å 和 19.2 Å），统计超胞原子数和位移数，解释 `min_length` 增大时计算量主要从哪里增加（提示：Si 的位移数由对称性决定）。
4. 用 `verdi presto` 建一个测试 profile，运行 `aiida_cp2k_sp.py`，然后用 `verdi node graph generate` 生成溯源图，指出图中哪些节点是输入、哪些是输出、哪个是 CalcJob。
5. 用 dpdispatcher 把 ENCUT 扫描改成 k 点扫描（KPOINTS 从 4×4×4 到 12×12×12），并把 `command` 换成调用 custodian 的脚本。
6. 结合你所在课题组的情况（人数、主要程序、是否有服务器、是否需要公开数据），按第10节写一份 200 字以内的工具选型说明。

---

## 参考资料

官方文档：

- custodian：<https://materialsproject.github.io/custodian/>
- atomate2：<https://materialsproject.github.io/atomate2/>；VASP 使用说明：<https://materialsproject.github.io/atomate2/user/codes/vasp.html>
- jobflow：<https://materialsproject.github.io/jobflow/>
- jobflow-remote：<https://matgenix.github.io/jobflow-remote/>
- AiiDA：<https://aiida.readthedocs.io/>；插件注册表：<https://aiida.net/plugin-registry/>
- AiiDA-VASP：<https://aiida-vasp.readthedocs.io/>
- quacc：<https://quantum-accelerators.github.io/quacc/>
- dpdispatcher：<https://docs.deepmodeling.com/projects/dpdispatcher/>
- pymatgen：<https://pymatgen.org/>；更新日志：<https://pymatgen.org/CHANGES.html>
- ASE 计算器：<https://docs.ase-lib.org/ase/calculators/calculators.html>
- Slurm sbatch：<https://slurm.schedmd.com/sbatch.html>

相关 GitHub 仓库（只链接，不复制代码）：

- [materialsproject/custodian](https://github.com/materialsproject/custodian) — VASP/CP2K/Gaussian/QChem 作业的自动纠错框架（MIT）。
- [materialsproject/atomate2](https://github.com/materialsproject/atomate2) — 基于 jobflow 的模块化工作流（弛豫、能带、声子、弹性、LOBSTER、缺陷）。
- [materialsproject/jobflow](https://github.com/materialsproject/jobflow) — `@job`/`Flow` 工作流编写库。
- [Matgenix/jobflow-remote](https://github.com/Matgenix/jobflow-remote) — 2026 年发布 1.0，把 jobflow/atomate2 工作流提交到 Slurm/PBS。
- [aiidateam/aiida-core](https://github.com/aiidateam/aiida-core) — 带完整溯源的工作流引擎。
- [aiida-vasp/aiida-vasp](https://github.com/aiida-vasp/aiida-vasp) — AiiDA 的 VASP 插件（relax、bands、NEB、收敛测试 WorkChain）。
- [aiidateam/aiida-cp2k](https://github.com/aiidateam/aiida-cp2k) — AiiDA 的 CP2K 插件。
- [nanotech-empa/aiida-gaussian](https://github.com/nanotech-empa/aiida-gaussian) — AiiDA 的 Gaussian 插件。
- [aiidateam/aiida-tutorials](https://github.com/aiidateam/aiida-tutorials) — AiiDA 官方教程源码。
- [Quantum-Accelerators/quacc](https://github.com/Quantum-Accelerators/quacc) — 基于 ASE 的配方式工作流平台（BSD-3-Clause）。
- [deepmodeling/dpdispatcher](https://github.com/deepmodeling/dpdispatcher) — Python 批量生成、提交、监控 Slurm/PBS/LSF 作业。
- [materialsproject/workshop](https://github.com/materialsproject/workshop) — Materials Project 官方 workshop，含 atomate/jobflow 课程。
- [tilde-lab/awesome-materials-informatics](https://github.com/tilde-lab/awesome-materials-informatics) — 材料信息学工具列表（AiiDA、atomate、数据库）。

关键文献：

- S. P. Ong et al., Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis（含 custodian 的介绍）, *Comput. Mater. Sci.* **68**, 314 (2013). DOI: 10.1016/j.commatsci.2012.10.028
- K. Mathew et al., Atomate: A high-level interface to generate, execute, and analyze computational materials science workflows, *Comput. Mater. Sci.* **139**, 140 (2017). DOI: 10.1016/j.commatsci.2017.07.030
- A. M. Ganose et al., Atomate2: modular workflows for materials science, *Digital Discovery* **4**, 1944 (2025). DOI: 10.1039/D5DD00019J
- A. S. Rosen et al., Jobflow: Computational Workflows Made Simple, *J. Open Source Softw.* **9**, 5995 (2024). DOI: 10.21105/joss.05995
- S. P. Huber et al., AiiDA 1.0, a scalable computational infrastructure for automated reproducible workflows and data provenance, *Sci. Data* **7**, 300 (2020). DOI: 10.1038/s41597-020-00638-4
- M. Uhrin et al., Workflows in AiiDA: Engineering a high-throughput, event-based engine for robust and modular computational workflows, *Comput. Mater. Sci.* **187**, 110086 (2021). DOI: 10.1016/j.commatsci.2020.110086
- A. Jain et al., FireWorks: a dynamic workflow system designed for high-throughput applications, *Concurr. Comput. Pract. Exp.* **27**, 5037 (2015). DOI: 10.1002/cpe.3505

[← 上一章](05-postprocess-visualization.md) | [返回目录](../README.md) | [下一章 →](07-learning-path.md)
