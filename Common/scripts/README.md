# Common 辅助脚本（Common helper scripts）

本目录收录四款软件（VASP、Gaussian 16、CP2K、OpenMX）都能用的通用脚本，涵盖作业脚本生成、高对称 k 路径、结构格式转换、批量作业状态汇总和计算环境自检。Python 脚本都可以用 `python <脚本> -h` 查看完整参数，Bash 脚本用 `bash env_check.sh -h` 查看。下文列出的参数全部取自脚本里的 argparse 定义和 usage 文本。

> [!IMPORTANT]
> 本 README 里出现的终端输出都是**示意输出（数值仅供参考）**，不是真实运行结果。脚本按 VASP 6.4–6.6、Gaussian 16 C.01/C.02、CP2K 2026.2（兼顾 2023.x–2025.x）、OpenMX 4.0/3.9 的输出格式编写。第一次在自己的集群上使用时，建议先拿一个小体系核对。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [mkjob.py](mkjob.py) | 按“软件 + 调度器 + 核数”生成 Slurm/PBS/LSF 作业脚本，包含 module、MPI × OpenMP、ulimit、CRLF 检查、scratch 目录和运行后检查 | Python ≥ 3.9，仅标准库 | [第2章](../tutorials/02-job-scheduler.md)、[第6章](../tutorials/06-workflow-automation.md) |
| [kpath_gen.py](kpath_gen.py) | 用 SeeK-path 或 ASE 生成标准高对称路径，一次输出 VASP KPOINTS/KPOINTS_OPT、杂化泛函零权重 KPOINTS、CP2K `&KPOINT_SET` 和 OpenMX `Band.kpath` | numpy、ase；seekpath + spglib（默认方法需要） | [第3章](../tutorials/03-structure-modeling.md)、[第4章](../tutorials/04-ase-pymatgen.md)、[第5章](../tutorials/05-postprocess-visualization.md) |
| [convert_any.py](convert_any.py) | 万能结构转换：CIF/POSCAR/xyz/extxyz/gjf/CP2K/OpenMX 互转，可选原胞/标准胞、超胞、加真空、按层固定原子 | numpy、ase；spglib（`--primitive/--conventional` 需要） | [第3章](../tutorials/03-structure-modeling.md)、[第4章](../tutorials/04-ase-pymatgen.md)、[第5章](../tutorials/05-postprocess-visualization.md)、[第7章](../tutorials/07-learning-path.md) |
| [job_status.py](job_status.py) | 批量扫描目录，判断作业是否正常结束、收敛、报错、运行中或中断，并汇总能量、步数和耗时 | Python ≥ 3.9，仅标准库；`--xlsx` 需要 pandas + openpyxl | [第2章](../tutorials/02-job-scheduler.md)、[第4章](../tutorials/04-ase-pymatgen.md)、[第6章](../tutorials/06-workflow-automation.md)、[第7章](../tutorials/07-learning-path.md) |
| [env_check.sh](env_check.sh) | 计算环境自检：检查 Python 和关键包版本、VASP_PP_PATH/PMG_VASP_PSP_DIR/MP_API_KEY、~/.vaspkit、ulimit -s、OMP_* 以及常用程序，每个问题都附修复建议 | bash ≥ 4、coreutils；python（用于报告版本） | [第1章](../tutorials/01-linux-hpc.md)、[第4章](../tutorials/04-ase-pymatgen.md)、[第5章](../tutorials/05-postprocess-visualization.md)、[第6章](../tutorials/06-workflow-automation.md)、[第7章](../tutorials/07-learning-path.md) |

## 安装依赖

```bash
# 推荐：使用仓库提供的 conda 环境（见第1章）
mamba env create -f ../examples/01-env-setup/environment.yml
conda activate dft

# 或者只装这几个脚本需要的包
pip install numpy ase spglib seekpath     # kpath_gen.py, convert_any.py
pip install pandas openpyxl               # 仅 job_status.py --xlsx 需要
```

`mkjob.py`、`job_status.py`（不加 `--xlsx`）和 `env_check.sh` 不需要任何第三方包，可以直接在集群登录节点上使用。上面的 `environment.yml` 见 [01-env-setup 示例](../examples/01-env-setup/README.md)。

> [!TIP]
> 脚本可以直接复制到集群上使用。如果是在 Windows 上编辑后上传的，先运行 `dos2unix *.py *.sh`。带 CRLF 换行的 bash 脚本会报 `$'\r': command not found`。

---

## 1. mkjob.py：生成作业脚本

### 用途

为 VASP / Gaussian 16 / CP2K / OpenMX 生成 Slurm、PBS（PBS Pro/OpenPBS，可选 Torque 语法）或 LSF 作业脚本。每次调用只写一个 LF 换行的脚本文件，不改动其他文件。生成的脚本包含以下内容：

- **调度器头部**：
  - Slurm：`#SBATCH` 的 `--nodes/--ntasks-per-node/--cpus-per-task/--mem/--time`，GPU 用 `--gres=gpu:N`。
  - PBS Pro：`select=N:ncpus=..:mpiprocs=..:ompthreads=..[:mem][:ngpus]`；加 `--torque` 时改成 `nodes=N:ppn=P`。
  - LSF：`-n`、`span[ptile=..]`、`-W [hour:]minute`。`rusage[mem]` 和 `affinity` 与站点配置有关，写成注释行。
- **USER SETTINGS 区块**：集中列出 `EXE`、`INPUT`、`OUTPUT`、`NP`、`NTHREADS`、`SCRATCH_ROOT`，方便手工修改。
- **运行环境**：
  - `module purge/load` 行（未指定时写占位行）和可选的 `source` 行。
  - `ulimit -s unlimited`，失败时退回到硬上限。
  - `OMP_NUM_THREADS`（Slurm 下取自 `SLURM_CPUS_PER_TASK`）；线程数大于 1 时额外设置 `OMP_PLACES/OMP_PROC_BIND`，并设置 `OMP_STACKSIZE`。
- **运行前检查**：
  - 输入文件存在，且不是 CRLF 换行（提示 `dos2unix`）。
  - VASP：INCAR/POSCAR/POTCAR 存在；INCAR 里没有 `KSPACING` 时要求 KPOINTS 存在；比较 POTCAR 数据集数与 POSCAR 元素种类数。
  - Gaussian：`%NProcShared/%CPU` 与申请的核数是否一致，有没有写 `%Mem`，文件末尾有没有空行。
  - CP2K：提示 `CP2K_DATA_DIR`。
  - OpenMX：`DATA.PATH` 目录是否存在。
- **启动命令**：
  - MPI 程序用 `srun --cpu-bind=cores` 或 `mpirun`（PBS 下加 `-machinefile $PBS_NODEFILE`）。
  - Gaussian 用 `g16 < in > log`。
  - CP2K 用 `cp2k.psmp -i ... -o ...`。
  - OpenMX 用 `openmx in.dat -nt $OMP_NUM_THREADS`。
- **可选 scratch 目录**：Gaussian 设置 `GAUSS_SCRDIR`，并用 trap 在退出时清理。其他软件把提交目录复制到 scratch 中运行，退出时（EXIT）再复制回来；多节点作业会给出警告。
- **运行后检查**：
  - VASP：查找 timing 块；`IBRION=1/2/3` 且 `NSW>0` 时还检查 `reached required accuracy`，并 grep 已知报错。
  - Gaussian：`Normal termination` 还是 `Error termination`，以及虚频。
  - CP2K：`PROGRAM ENDED AT`、`ABORT`、`SCF run NOT converged`、优化相关信息，以及进程/线程数行。
  - OpenMX：System.Name.out 中的 `Elapsed.Time.` 和 .md 中的 `scf_conv=0`。

### 依赖

Python ≥ 3.9，只用标准库（argparse、`string.Template`）。模板以 `@@` 作分隔符，所以 bash 里的 `$` 不需要转义。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `PROGRAM`（位置参数） | — | `vasp`、`gaussian`（别名 `g16`、`gau`）、`cp2k`（别名 `cp2k.psmp`）、`openmx`（别名 `omx`） |
| `-s, --scheduler {slurm,pbs,lsf}` | `slurm` | 调度系统 |
| `-i, --input FILE` | VASP 不需要；其他：`input.gjf/.inp/.dat` | 主输入文件。除 VASP 外必须给出，文件名不能含空格或引号 |
| `--output FILE` | `vasp.out`；其他为 `<输入名>.log/.out/.std` | 程序输出/日志文件 |
| `-N, --nodes N` | 1 | 节点数。Gaussian 固定为 1 |
| `-n, --ntasks-per-node N` | 32（Gaussian：1） | 每节点 MPI 进程数 |
| `-c, --threads N` | 1（Gaussian：核数，默认 16） | 每进程 OpenMP 线程数，即 cpus-per-task |
| `--mem 180G` | 不写 | 每节点内存。Gaussian 会提示把 `%Mem` 设为约 80% |
| `-t, --time` | `24:00:00` | 墙钟时间，格式 `D-HH:MM:SS`、`HH:MM:SS`、`HH:MM` 或分钟数 |
| `-q, --queue, --partition` | — | 分区（Slurm）或队列（PBS、LSF） |
| `-A, --account` | — | 计费账户或项目 |
| `--qos` | — | Slurm QOS |
| `-J, --name` | 程序名或输入文件主干 | 作业名 |
| `--gpus N` | 0 | 每节点 GPU 数（VASP GPU 版每张卡 1 个进程） |
| `--exclusive` | 关 | 独占节点 |
| `--torque` | 关 | PBS 使用 Torque 语法 `nodes=N:ppn=P` |
| `--module MOD` | — | 要加载的 module，可重复 |
| `--purge` | 关 | 在 module load 前加 `module purge` |
| `--source FILE` | — | 需要 source 的文件（如 CP2K toolchain 的 `install/setup`），可重复 |
| `--exe EXE` | `vasp_std`、`g16`、`cp2k.psmp`、`openmx` | 可执行文件 |
| `--launcher {srun,mpirun}` | Slurm 用 `srun`，其他用 `mpirun` | MPI 启动器。`srun` 只能和 Slurm 一起用，搭配 PBS/LSF 时脚本直接报错 |
| `--scratch DIR` | 不用 | scratch 根目录，如 `'${TMPDIR:-/tmp}'` 或 `/scratch` |
| `--omp-stacksize` | `512m` | 写入 `OMP_STACKSIZE` 的值 |
| `-o, --output-script FILE` | `job_<program>_<scheduler>.sh` | 写出的脚本文件名；`-` 表示打印到标准输出 |
| `-f, --force` | 关 | 覆盖已存在的脚本文件 |

### 用法示例

```bash
# VASP：2 节点 x 64 进程，Slurm 分区 cpu，加载 vasp/6.5.1
python mkjob.py vasp -N 2 -n 64 -t 24:00:00 -q cpu --module "vasp/6.5.1"

# 同样的 VASP 作业改用 PBS Pro，每节点 200G 内存
python mkjob.py vasp -s pbs -N 1 -n 64 --mem 200G -o job.pbs

# Gaussian：单节点 16 核，40G 内存，使用节点本地 scratch
python mkjob.py gaussian -i benzene.gjf -c 16 --mem 40G --scratch '${TMPDIR:-/tmp}'

# CP2K：MPI x OpenMP 混合并行（每节点 32 进程 x 2 线程），source toolchain 环境
python mkjob.py cp2k -i h2o_md.inp -N 2 -n 32 -c 2 --source ~/cp2k/tools/toolchain/install/setup

# OpenMX：LSF 队列 normal
python mkjob.py openmx -i Si.dat -n 32 -c 2 -s lsf -q normal

# 只打印到终端，先检查再保存
python mkjob.py cp2k -i opt.inp --launcher mpirun -o -
```

生成后，根据调度器用 `sbatch job.sh`、`qsub job.sh` 或 `bsub < job.sh` 提交。脚本结束时会打印提交命令。

### 输入与输出

- 输入：只有命令行参数。生成脚本时不读取 INCAR、.gjf 等文件，这些检查写进作业脚本，在作业运行时执行。
- 输出：一个带可执行权限的作业脚本。文件已存在时不覆盖，除非加 `-f`。

### 局限与注意事项

- 分区/队列名、账户、module 名、GPU 写法（`--gres`、`ngpus`、`-gpu`）和 MPI 启动方式都和具体集群有关，一定要对照本集群的用户手册修改。
- Gaussian 16（非 Linda 版）只能单节点运行：脚本强制 1 节点、1 个进程，`-c` 就是 `%NProcShared` 的值。
- VASP 用 `-c > 1` 时需要 OpenMP 编译版本，此时 VASP 会自己把 NCORE 设为 1。
- 多节点作业使用 `--scratch` 时，该目录必须在所有节点都能访问的文件系统上（节点本地 `/tmp` 不行）。
- 更多调度器知识和完整的手写模板见 [第2章](../tutorials/02-job-scheduler.md) 和 [02-job-templates 示例](../examples/02-job-templates/README.md)；各软件的专用脚本见 [VASP/scripts](../../VASP/scripts/README.md)、[Gaussian/scripts](../../Gaussian/scripts/README.md)、[CP2K/scripts](../../CP2K/scripts/README.md)、[OpenMX/scripts/job_omx_slurm.sh](../../OpenMX/scripts/job_omx_slurm.sh)。

---

## 2. kpath_gen.py：高对称 k 路径

### 用途

用 SeeK-path（HPKOT 约定，默认）或 ASE（Setyawan–Curtarolo 约定）生成能带计算用的标准高对称路径，一次写出三款周期性软件的输入片段。

> [!WARNING]
> k 点坐标只有和它所对应的晶胞一起使用才有意义。默认方法（seekpath）给出的坐标对应 SeeK-path 的**标准化原胞**，这个原胞可能相对你的输入做了旋转或轴重排，它会写入 `POSCAR_prim`。SCF 和能带计算都要用这个晶胞（可以用 [convert_any.py](convert_any.py) 转换成 CP2K/OpenMX 格式）。如果想沿用自己的晶胞，加 `--orig-cell`（需要 seekpath ≥ 2.1）。

### 依赖

```bash
pip install numpy ase seekpath spglib    # --method ase 时可以不装 seekpath/spglib
```

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `STRUCTURE`（位置参数） | — | 任意 ASE 可读的结构文件（POSCAR、CIF、extxyz 等），必须带三维晶胞 |
| `--format FMT` | 自动识别 | 识别失败时指定 ASE 格式名，如 `vasp`、`cif` |
| `--method {seekpath,ase}` | `seekpath` | 路径约定 |
| `--orig-cell` | 关 | 仅 seekpath：把路径表示在输入晶胞的倒格子基下（需要 seekpath ≥ 2.1，用 `hasattr` 检测） |
| `--symprec` | `1e-3` | spglib 对称性容差（Å） |
| `--no-time-reversal` | 关 | 仅 seekpath：不假定时间反演对称（体系没有反演中心时会加入 −k 段） |
| `-n, --npoints N` | 40 | 每段点数，用于 line-mode KPOINTS、CP2K `NPOINTS` 和 OpenMX |
| `--density PTS` | 不用 | 每 1/Å 的点数，用于 OpenMX 和显式零权重列表（按段长分配） |
| `--hse-npoints N` | 20 | 零权重列表中每段的点数 |
| `--ibzkpt FILE` | — | SCF 计算得到的 IBZKPT，指定后写出 `KPOINTS_hse` |
| `-o, --outdir DIR` | `kpath` | 输出目录 |
| `--print` | 关 | 同时把各个文本块打印到终端 |

### 用法示例

```bash
python kpath_gen.py POSCAR                       # seekpath，每段 40 点
python kpath_gen.py Si.cif -n 60 -o si_kpath
python kpath_gen.py CONTCAR --orig-cell          # 路径用输入晶胞的倒格子基表示
python kpath_gen.py POSCAR --ibzkpt scf/IBZKPT --hse-npoints 15
python kpath_gen.py POSCAR --density 30          # 按 1/Å 密度分配点数（CP2K/OpenMX/HSE 列表）
python kpath_gen.py POSCAR --method ase          # 不装 seekpath，用 ASE 的 bandpath
python kpath_gen.py POSCAR --print               # 同时在终端打印所有文本块
```

### 输出文件（写在 `--outdir` 下）

| 文件 | 内容 | 用法 |
|---|---|---|
| `POSCAR_prim` | seekpath 标准化原胞，即 k 点坐标对应的晶胞（`--orig-cell` 或 `--method ase` 时不写） | SCF 与能带计算都用这个晶胞 |
| `KPOINTS_band` | VASP line-mode KPOINTS | 非自洽能带计算（`ICHARG = 11`） |
| `KPOINTS_OPT` | 与 `KPOINTS_band` 内容相同（KPOINTS_OPT 采用 KPOINTS 格式，需要 VASP ≥ 6.3） | 在 SCF 计算中直接得到能带，例如 HSE06 |
| `KPOINTS_hse` | 只在给出 `--ibzkpt` 时写出：IBZKPT 的权重 k 点，加上带标签的零权重路径点 | 经典的杂化泛函能带算法 |
| `cp2k_kpoint_set.inc` | 每段连续子路径一个 `&KPOINT_SET`（`UNITS B_VECTOR`、`NPOINTS`、`SPECIAL_POINT 名称 x y z`） | 放进 `&DFT/&PRINT/&BAND_STRUCTURE` |
| `openmx_band_kpath.txt` | `Band.dispersion`、对应晶胞的 `<Band.KPath.UnitCell`、`Band.Nkpath`、`<Band.kpath`（每行 `n k1 k2 k3 k1' k2' k3' L1 L2`，GAMMA 写作 G） | 粘贴到 OpenMX 的 .dat |
| `kpath_info.txt` | 空间群、Bravais 格子、高对称点（分数坐标和 \|k\|）、各段长度和警告 | 核对路径用 |

配套的能带作图脚本：VASP 用 [plot_band.py](../../VASP/scripts/plot_band.py)，OpenMX 用 [omx_band_plot.py](../../OpenMX/scripts/omx_band_plot.py)。

### 局限与注意事项

- `--method ase` 只利用晶格对称性，会忽略原子降低的对称性，坐标对应输入晶胞。需要严格的晶体对称性时用 seekpath。
- 如果输入是超胞，脚本会在 `kpath_info.txt` 里提示：这时标签并不是超胞布里渊区的真正高对称点。
- `KPOINTS_hse` 只有在 IBZKPT 来自**同一个晶胞**（以及同样的 k 网格和 ISYM）时才有效。`--ibzkpt` 只接受倒空间坐标（Reciprocal）的显式列表。
- 推荐阅读 [第3章](../tutorials/03-structure-modeling.md)（原胞与标准化）和 [第5章](../tutorials/05-postprocess-visualization.md)（能带作图）。参考文献：Hinuma et al., Comput. Mater. Sci. 128, 140 (2017)，<https://doi.org/10.1016/j.commatsci.2016.10.015>；Setyawan & Curtarolo, Comput. Mater. Sci. 49, 299 (2010)，<https://doi.org/10.1016/j.commatsci.2010.05.010>。SeeK-path 项目：<https://github.com/giovannipizzi/seekpath>。

---

## 3. convert_any.py：万能结构转换

### 用途

在 CIF、POSCAR、xyz/extxyz、Gaussian .gjf、CP2K `&CELL/&COORD` 和 OpenMX .dat 之间互相转换，转换时可以同时取原胞或标准惯用胞、建超胞、加真空、按层固定原子。固定原子信息会写成对应格式的约束：Selective dynamics、`&FIXED_ATOMS`、`<MD.Fixed.XYZ` 或 Gaussian 冻结码。

处理顺序固定为：读取 → `--primitive/--conventional` → `--supercell` → `--vacuum/--center` → `--wrap` → `--sort` → 固定原子 → 写出。

**可读格式**

- 所有 ASE 能读的格式（CIF、POSCAR/CONTCAR、vasprun.xml、OUTCAR、xyz、extxyz、.gjf/.com、Gaussian .log、PDB 等）。
- CP2K 输入/restart 文件（`*.inp`、`*.restart`、`*.inc`）：读取 `&CELL`（A/B/C 或 ABC + ALPHA_BETA_GAMMA，支持 `[angstrom]/[bohr]` 单位和 PERIODIC）和 `&COORD`（笛卡尔或 SCALED、UNIT）。Fe1、O_w 这类 kind 名通过 `&KIND/ELEMENT` 或开头字母映射到元素。
- OpenMX 输入（`*.dat`、`*.dat#`）：`Atoms.UnitVectors(.Unit)`、`Atoms.SpeciesAndCoordinates(.Unit Ang|AU|FRAC)`，元素从 VPS 名称推断（如 `Fe_PBE19` → Fe）。

**可写格式**（由文件名或 `-t` 决定）

| 格式 | 文件名 | 说明 |
|---|---|---|
| `vasp` | `POSCAR`、`CONTCAR`、`*.vasp`、`*POSCAR*` | VASP 5 格式，Direct 坐标，有固定原子时写 `Selective dynamics` |
| `cif` | `*.cif` | 通过 ASE 写出 |
| `xyz` | `*.xyz` | 普通 xyz，不含晶胞 |
| `extxyz` | `*.extxyz` | extended xyz，含晶胞和 pbc |
| `gjf` | `*.gjf`、`*.com` | Gaussian 输入，周期性体系写 `Tv` 行，固定原子冻结码为 −1，末尾带空行 |
| `cp2k` | `*.inc`、`*.inp` | 供 `&SUBSYS` 中 `@INCLUDE` 使用的 `&CELL` + `&COORD`；有固定原子时另写 `<name>_fixed.inc`（`&FIXED_ATOMS`，放进 `&MOTION/&CONSTRAINT`） |
| `openmx` | `*.dat` | `Atoms.*` 块（有固定原子时加 `<MD.Fixed.XYZ`），species 定义以注释模板形式给出 |
| `ase:<格式>` | 任意 | 其他 ASE 支持的格式 |

### 依赖

```bash
pip install numpy ase spglib     # spglib 只有 --primitive/--conventional 需要
```

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `INPUT`（位置参数） | — | 输入结构文件 |
| `-o, --output FILE [FILE ...]` | — | 一个或多个输出文件，格式由文件名判断 |
| `-t, --to FMT [FMT ...]` | — | 输出格式，文件名取 `<输入文件主干>` 加对应扩展名（如 `Si.vasp`、`Si.inc`）；只给一个 `-o` 时，`-t` 指定该文件的格式 |
| `--in-format` | 自动 | 输入格式：`cp2k`、`openmx` 或任意 ASE 格式名 |
| `--index` | `-1` | 轨迹或 OUTCAR 中取哪一帧（ASE 索引写法） |
| `--primitive` / `--conventional` | 关 | spglib 标准化原胞或标准惯用胞（二者互斥） |
| `--symprec` | `1e-3` | spglib 容差（Å） |
| `--supercell N ...` | — | 3 个整数（`2 2 1` 或 `2x2x1`），或 9 个整数按行给出 3×3 矩阵 |
| `--vacuum ANG` | — | 沿 `--vacuum-axis` 的**总**真空厚度，两侧各一半 |
| `--vacuum-axis` | `c` | `--vacuum/--center` 的方向（a、b、c） |
| `--center` | 关 | 沿 `--vacuum-axis` 居中，不改变晶胞 |
| `--box ANG` | 10 | 没有晶胞的分子写成 POSCAR/CP2K/OpenMX/CIF 时，在分子尺寸外加的盒子余量 |
| `--wrap` | 关 | 把原子移回晶胞内 |
| `--sort` | 关 | 按元素排序（写 POSCAR 时推荐） |
| `--fix-layers N` | — | 固定沿 `--fix-axis` 最低的 N 层 |
| `--fix-below ANG` | — | 固定高度低于 ANG 的原子（高度从晶胞原点算起） |
| `--fix-indices LIST` | — | 按 1 起始的序号固定，如 `1-8,12` |
| `--fix-axis` | `c` | 表面法向（a、b、c） |
| `--layer-tol` | 0.5 | 沿法向间距小于该值（Å）的原子归为同一层 |
| `--unfix` | 关 | 丢弃输入中读到的所有固定标记 |
| `--route` | `B3LYP/6-31G(d) Opt`；周期性：`PBEPBE/6-31G(d) SCF=Tight` | Gaussian route 行 |
| `--charge` / `--mult` | 0 / 1 | 电荷与自旋多重度 |
| `--nproc` / `--mem` | 8 / `8GB` | `%NProcShared` 与 `%Mem` |
| `--no-tv` | 关 | 不写 `Tv` 平移矢量行 |
| `--valence El=N ...` | 部分轻元素内置 | OpenMX 赝势价电子数，初始自旋上/下各取 N/2，如 `Ti=12 O=6` |
| `--omx-unit {Ang,FRAC}` | `Ang` | `Atoms.SpeciesAndCoordinates.Unit` |
| `-f, --force` | 关 | 覆盖已存在的输出文件 |

### 用法示例

```bash
python convert_any.py Si.cif -o POSCAR --primitive
python convert_any.py Si.cif -t vasp cp2k openmx gjf          # Si.vasp, Si.inc, Si.dat, Si.gjf
python convert_any.py POSCAR -o slab.vasp --supercell 3 3 1 --vacuum 15 --fix-layers 2
python convert_any.py slab.vasp -o slab.inc --fix-below 8.0   # CP2K + slab_fixed.inc
python convert_any.py h2o.xyz -o h2o.gjf --route "B3LYP/6-31G(d) Opt Freq" --nproc 8 --mem 16GB
python convert_any.py h2o.xyz -o POSCAR --box 12              # 分子放进 12 A 余量的盒子
python convert_any.py cp2k_opt-1.restart -o CONTCAR.vasp --sort
python convert_any.py Si.dat# -o relaxed.cif                  # OpenMX 结束后的 .dat# -> CIF
python convert_any.py OUTCAR --index -1 -o final.xyz
```

终端会打印读入的化学式、原子数、pbc、晶胞长度、固定原子数，以及写出后的体积和密度，用来快速核对。

### 局限与注意事项

- `--primitive/--conventional` 和一般矩阵超胞（9 个整数）会丢弃固定标记，脚本会给出提示。
- xyz、extxyz、CIF 以及 `ase:` 格式不保存固定标记。
- 写出的 CP2K `&FIXED_ATOMS` 的 `COMPONENTS_TO_FIX` 对应笛卡尔 X/Y/Z，只有正交晶胞时才等同于晶格方向；它只对 GEO_OPT、CELL_OPT 和 MD 起作用。`&KIND` 的基组和赝势只以注释模板给出，必须自己确认。
- OpenMX 的 species 定义（PAO/VPS）以注释模板写出。没有内置价电子数的元素，自旋列会写成 `?.?`，要用 `--valence` 补上（参考 DFT_DATA19 网页上的 VPS 表，很多金属包含半芯态）。
- Gaussian 的冻结码 −1 只在 `Opt` 中生效。
- 详细的建模流程见 [第3章](../tutorials/03-structure-modeling.md) 和 [03-fetch-convert-si 示例](../examples/03-fetch-convert-si/README.md)；ASE/pymatgen 的写法见 [第4章](../tutorials/04-ase-pymatgen.md)。

---

## 4. job_status.py：批量作业状态汇总

### 用途

递归扫描一个或多个目录，识别其中的 VASP / Gaussian 16 / CP2K / OpenMX 作业，用一张表汇总状态、能量、离子步（优化步）数和耗时。文件按行流式读取，几 GB 的 OUTCAR 也能处理。

| 软件 | 识别方式 | 判断依据 |
|---|---|---|
| VASP 6.x | 含 OUTCAR 的目录 | 结束：`General timing and accounting`。耗时：`Elapsed time (sec):`；未结束时累加 `LOOP+` 的 real time。能量：最后一个 `energy(sigma->0)`。步数：`free  energy   TOTEN` 行数。收敛：`IBRION=1/2/3` 且 `NSW>0` 时要求 `reached required accuracy`，并报告电子步跑满 NELM 的离子步。报错：ZBRENT、EDDDAV/EDDRMM（ZHEGV）、Sub-Space-Matrix、BRMIX、VERY BAD NEWS、Tetrahedron method fails、ZPOTRF、SIGSEGV、Slurm TIME LIMIT/OOM 等，在 OUTCAR 和 vasp.out、slurm-*.out 等标准输出文件中查找，VASP 6 跨行的方框报错会先合并再匹配 |
| Gaussian 16 | `*.log/*.out` 含 Gaussian 标识 | 每个 Link1 步的 `Normal termination` / `Error termination`（以及 link 编号）、最后一个 `SCF Done`、`Step number N out of a maximum of M`、`Optimization completed` / `Number of steps exceeded`、虚频、`Elapsed time:`（各 Link1 步累加） |
| CP2K | `*.out/*.log` 含 `PROGRAM STARTED AT` 或 `CP2K\|` | `PROGRAM ENDED AT`、`[ABORT]`、`SCF run NOT converged`、`ENERGY\| Total FORCE_EVAL ( QS ) energy [hartree]`（2025.1+）或 `[a.u.]`（旧版）、`OPT\| Step number`（2023.1+）或 `Informations at step`（旧版）、`MD\| Step number`、`GEOMETRY OPTIMIZATION COMPLETED` / `MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED`、timing 表中的 CP2K 行 |
| OpenMX 4.0/3.9 | `System.Name.out` | `Elapsed.Time.` 表示正常结束；最后一个 `Utot.`（Hartree）；System.Name.md 中每步的 `scf_conv=0/1`；步数与输入 .dat 中的 `MD.maxIter` 比较，优化停在 `MD.maxIter` 时会标记 |

状态取值：

| 状态 | 含义 |
|---|---|
| `DONE` | 正常结束并已收敛（或者不需要收敛判断） |
| `NOT_CONVERGED` | 正常结束，但 SCF 或几何优化没有收敛 |
| `ERROR` | 找到报错信息，或 Error termination |
| `RUNNING` | 未结束，且输出在 `--running-min` 分钟内有更新 |
| `INCOMPLETE` | 未结束，且输出已经很久没有更新（被杀、超时、节点故障等） |

### 依赖

Python ≥ 3.9 标准库。只有 `--xlsx` 需要 `pip install pandas openpyxl`。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `ROOTS ...`（位置参数） | `.` | 要扫描的目录，也可以是单个输出文件 |
| `--depth N` | 不限 | 每个根目录下最大扫描深度 |
| `--codes` | `vasp,gaussian,cp2k,openmx` | 只识别其中一部分软件，逗号分隔 |
| `--unit {eV,Ha}` | `eV` | 终端表格中的能量单位 |
| `--sort {path,status,energy,time}` | `path` | 排序方式 |
| `--only-bad` | 关 | 只显示 ERROR / NOT_CONVERGED / INCOMPLETE |
| `--running-min` | 10 | 输出在多少分钟内有更新算作 RUNNING |
| `--csv FILE` | — | 把表格写成 CSV |
| `--json FILE` | — | 把全部字段写成 JSON |
| `--xlsx FILE` | — | 写 Excel 文件（需要 pandas + openpyxl） |
| `--no-notes` | 关 | 终端不显示 notes 列 |

### 用法示例

```bash
python job_status.py                            # 扫描当前目录树
python job_status.py relax/ scf/ band/ --depth 2
python job_status.py ~/projects/TiO2 --only-bad
python job_status.py . --csv status.csv --json status.json
python job_status.py . --unit Ha --sort energy
python job_status.py . --xlsx status.xlsx       # 需要 pandas + openpyxl
```

（示意输出，数值仅供参考）

```text
path            code     status        task       E(eV)              steps        time  notes
relax/OUTCAR    vasp     DONE          relax              -43.381200     12    00:41:05
scf/h2o.log     gaussian ERROR         RB3LYP           -2078.957000      7    00:03:10  Error termination in l9999

2 job(s): ERROR=1, DONE=1
```

### 输入与输出

- 终端：每行一个作业，最后一行是各状态的计数。
- CSV/JSON：字段为 `path, code, status, task, energy_eV, steps, time_s, version, notes`。注意文件中的能量**始终是 eV**（`energy_eV`），`--unit` 只影响终端显示。
- 退出码：0 表示全部为 DONE/RUNNING；1 表示至少有一个 ERROR/NOT_CONVERGED/INCOMPLETE；2 表示没有找到作业（或 `--codes` 写错）。可以直接用在 shell 判断或工作流脚本里。

### 局限与注意事项

- 隐藏目录（以 `.` 开头）和以 `_rst` 结尾的目录会被跳过。名为 `vasp.out` 和以 `slurm-` 开头的文件不会被当成 Gaussian/CP2K/OpenMX 输出，只在查找 VASP 报错时使用。
- VASP 按目录识别：一个目录里只看一个 OUTCAR，VASP 的 NEB 子目录（00、01……）会分别列出。
- OpenMX 的 `.out` 文件如果头部没有 OpenMX 字样，要靠同目录下 `System.Name` 与文件名主干一致的 .dat 才能识别。
- RUNNING/INCOMPLETE 只看文件修改时间，不查询调度器；要确认作业是否还在队列里，用 `squeue -u $USER`、`qstat` 或 `bjobs`（见 [第2章](../tutorials/02-job-scheduler.md)）。
- 单个作业的深入检查可以用各软件板块的脚本，例如 [VASP/scripts](../../VASP/scripts/README.md) 中的 `check_vasp.py`、[Gaussian/scripts](../../Gaussian/scripts/README.md) 中的 `g16_log_summary.py`。批量工作流中的用法见 [第6章](../tutorials/06-workflow-automation.md)。

---

## 5. env_check.sh：计算环境自检

### 用途

在登录节点、WSL 或工作站上做一次环境体检，每一项都给出 `[ OK ]`、`WARN` 或 `FAIL`，有问题时附上修复命令。脚本只读取文件和环境变量，**不修改任何东西**。

| 检查组 | 检查内容 |
|---|---|
| System | 可见 CPU 核数、bash 版本（需要 ≥ 4） |
| Python | 解释器是否可用，版本是否 ≥ 3.9，当前 conda 环境（在 `base` 中工作会给出警告） |
| Python packages | ase、pymatgen、phonopy、numpy、matplotlib、spglib（缺少时 WARN）；sumo、pyprocar、py4vasp、scipy、seekpath、mp-api（可选，缺少时只提示）。版本通过 `importlib.metadata` 读取，不 import 包 |
| VASP / pymatgen / MP | `VASP_PP_PATH`（应该是**包含** `potpaw_PBE/` 的目录）；`PMG_VASP_PSP_DIR`（环境变量或 `~/.config/.pmgrc.yaml` / `~/.pmgrc.yaml`）及其中的 `POT_*` 子目录；.pmgrc.yaml 是否对其他用户可读；`MP_API_KEY` / `PMG_MAPI_KEY` 长度（只显示长度，不显示密钥）；是否残留旧的 `MP_API_ENDPOINT` |
| VASPKIT | `~/.vaspkit` 是否存在、其中的 `PBE_PATH` 是否存在、是否为 CRLF 换行，vaspkit 是否在 PATH 中 |
| Stack size and OpenMP | `ulimit -s`（软/硬上限）、`OMP_NUM_THREADS` 是否设置以及是否超过核数、`OMP_STACKSIZE`（Multiwfn 处理大体系约需 200M） |
| Programs and schedulers | vasp_std/vasp_gam/vasp_ncl、g16/formchk（以及 `g16root`、`GAUSS_SCRDIR`）、cp2k.psmp（以及 `CP2K_DATA_DIR`）、openmx、Multiwfn（以及 `Multiwfnpath`）、mpirun/srun、sbatch/qsub/bsub、dos2unix |

### 依赖

bash ≥ 4 和 coreutils；系统里有 python 时才能报告 Python 和包的版本。

### 参数

| 参数 | 说明 |
|---|---|
| `-p, --python PY` | 要检查的 Python 解释器（默认先找 `python3`，再找 `python`） |
| `-q, --quiet` | 只打印 WARN、FAIL 行和最后的汇总 |
| `--no-color` | 纯文本输出（标准输出不是终端时自动启用） |
| `-h, --help` | 显示帮助 |

### 用法示例

```bash
bash env_check.sh                                   # 完整报告
bash env_check.sh -q                                # 只看问题
bash env_check.sh -p ~/miniforge3/envs/dft/bin/python
bash env_check.sh --no-color > env_report.txt       # 保存报告，提问时附上
```

（示意输出，数值仅供参考）

```text
== Stack size and OpenMP ==
[WARN] ulimit -s = 8192 kB (hard limit unlimited); Fortran codes (VASP, CP2K, Multiwfn) may crash with SIGSEGV
       -> fix: put 'ulimit -s unlimited' in every job script (and in ~/.bashrc for interactive tools)
```

### 输入与输出

- 输入：当前 shell 的环境变量，以及 `~/.vaspkit`、`~/.pmgrc.yaml` 等配置文件。
- 输出：打印到终端的报告。退出码 0 表示没有 FAIL（允许有 WARN），1 表示至少有一个 FAIL；参数错误时退出码为 2。

### 局限与注意事项

- 在登录节点上运行时，看到的是登录节点的环境。计算节点的 module、`ulimit` 和核数可能不同，可以在作业脚本里也运行一次（例如 `bash env_check.sh -q --no-color`）。
- 只检查程序是否在 PATH 中、环境变量是否设置，不验证程序能否正确运行。
- 修复建议中提到的配置模板（`environment.yml`、`vaspkit_rc.example`、`pmgrc.yaml.example`）在 [01-env-setup 示例](../examples/01-env-setup/README.md) 中；环境搭建的完整步骤见 [第1章](../tutorials/01-linux-hpc.md)。

---

## 许可与反馈

这些脚本属于 [DFT-Tutorials](https://github.com/zenggs05/DFT-Tutorials)，采用 MIT 协议。发现解析错误（尤其是新版本软件的输出格式变化）时，欢迎在 GitHub Issues 中附上出错的输出片段。

[返回目录](../../README.md)
