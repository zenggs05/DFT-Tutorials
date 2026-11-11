# 第17章 作业脚本、并行与批量计算（Running Gaussian 16 on HPC: Job Scripts, Parallelism & Batch Jobs）

本章讲解如何在超算集群上高效、稳妥地运行 Gaussian 16：Slurm 与 PBS 作业脚本模板，`%nprocshared`/`%cpu`、`%mem` 与调度系统申请资源的对应关系，`GAUSS_CDEF`、`GAUSS_MDEF` 等环境变量与 `Default.Route` 的用法，每个作业独立的 scratch 目录，Linda 跨节点并行和 GPU 计算，批量提交、作业数组与依赖作业，作业监控，以及 `Opt=Restart`、`IRC=Restart`、`Restart` 关键词等断点续算方法。最后用实测方法说明为什么“核数越多越快”并不成立。关键词：Gaussian Slurm script、PBS、nprocshared、%mem、GAUSS_SCRDIR、GAUSS_CDEF、Linda、Gaussian GPU、job array、Opt=Restart。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 写出并看懂 Gaussian 16 的 Slurm 与 PBS 作业脚本，知道每一行 `#SBATCH`/`#PBS` 的作用；
2. 让 Gaussian 实际使用的核数、内存与调度系统分配的资源一致，避免超额占用（oversubscription）和内存不足；
3. 理解 Link 0 行、命令行选项、环境变量、`Default.Route` 四种设置方式的优先级，选择适合集群的写法；
4. 为每个作业建立独立的 scratch 目录并在结束时清理，必要时用 `%RWF` 拆分大文件；
5. 了解 Linda 多节点并行和 GPU 计算的适用范围与写法；
6. 批量提交几十上百个输入，使用作业数组和依赖作业串联多步计算；
7. 监控作业进度，作业中断后用 `Opt=Restart`、`IRC=Restart`、`Restart` 继续计算；
8. 用缩放测试为自己的体系选择合适的核数。

**前置知识**：Linux 命令行、集群登录与文件传输（[Common 第1章](../../Common/tutorials/01-linux-hpc.md)），调度系统的一般概念与模板（Common 第2章），Gaussian 安装与环境变量（[第1章](01-setup-first-job.md)），Link 0 命令（[第2章](02-input-anatomy.md)）。

本章配套示例：[../examples/17-hpc-job-scripts/README.md](../examples/17-hpc-job-scripts/README.md)（Slurm/PBS/GPU/Linda 模板、批量提交、作业数组、依赖作业、缩放测试、`Default.Route` 示例和断点续算输入）。

> [!NOTE]
> 集群之间差异很大：模块名、分区名、scratch 路径、是否允许 `--exclusive`、单用户最多排队多少作业，都要以你所在集群的使用手册为准。本章脚本中所有“site specific”的行都需要按实际情况修改。

---

## 2. 原理速览：Gaussian 怎样并行

### 2.1 三种并行方式

| 方式 | Link 0 | 作用范围 | 适用 |
|---|---|---|---|
| 共享内存并行（SMP） | `%CPU=0-15` 或 `%NProcShared=16` | 一个节点内，一个进程的多个线程 | 几乎所有作业类型，**首选** |
| Linda 集群并行 | `%LindaWorkers=node1,node2` | 多个节点，每个节点一个（或几个）worker | 大体系 HF/DFT 能量、梯度、频率；MP2 能量与梯度；需单独的 Linda 授权 |
| GPU | `%GPUCPU=0-3=0,1,16,17` | 一个节点内的 NVIDIA GPU | 较大分子的 HF/DFT 能量、梯度、频率（含激发态），仅 Linux |

官方建议：节点内用 SMP，节点之间才用 Linda；不要在一个节点上开多个 Linda worker 代替 SMP。

### 2.2 内存是所有线程共享的

`%mem` 是整个作业（一个节点）的动态内存总量，Gaussian 把它平均分给各个线程。Gaussian 16 默认只有 800 MB（100 MW），对稍大的计算远远不够。官方的建议是：

- 50 个原子以上或 500 个基函数以上的计算，每个核 4 GB 或更多；
- 内存应随核数线性增加：1 核用 4 GB 合适的作业，8 核时用 32 GB 才能保持效率；只增加核数而不增加内存，大规模并行时性能不会好；
- 大的频率计算和 CCSD/EOM-CCSD 需要把一部分物理内存留给操作系统做磁盘缓存，官方建议此时只把整机内存的 50–70% 交给 Gaussian；
- **请求超过物理可用内存会导致性能极差或作业被杀掉**。

### 2.3 并行效率与 Amdahl 定律

如果一个程序中比例为 $p$ 的部分可以并行，用 $N$ 个核的理论加速比为

$$
S(N)=\frac{1}{(1-p)+p/N}
$$

即使 $p=0.95$，$N\to\infty$ 时最大加速比也只有 20。Gaussian 中积分计算、DFT 格点积分、CPHF 等部分并行性很好，但还有串行部分、线程同步和内存带宽的限制。小分子（几十个原子、def2-SVP）在 8–16 核之后加速比明显变平；大分子在内存充足、线程绑定到具体核的前提下，官方称并行效率可以保持到 64 核或更多。所以**核数要按体系大小测出来，而不是越多越好**（第 13 节）。

### 2.4 磁盘：scratch 文件

Gaussian 运行时会在 `GAUSS_SCRDIR` 指定的目录中生成以进程号命名的临时文件：

| 文件 | 内容 | Link 0 |
|---|---|---|
| `Gau-NNNN.rwf` | read-write 文件，通常最大 | `%RWF` |
| `Gau-NNNN.int` | 双电子积分（默认为空） | `%Int` |
| `Gau-NNNN.d2e` | 积分导数（默认为空） | `%D2E` |
| `Gau-NNNN.skr` | 临时文件 | — |
| `Gau-NNNN.chk` | 未命名的检查点 | `%Chk` |

作业正常结束时，未命名的临时文件会被删除；**作业被杀掉或异常结束时不会删除**，久而久之会塞满磁盘。`GAUSS_SCRDIR` 没有设置时，scratch 文件放在当前目录。

---

## 3. 四种设置方式与优先级

Gaussian 16 的大多数运行参数都可以用四种方式指定，优先级从高到低为：

1. 输入文件中的 Link 0 行（`%` 开头）；
2. `g16` 命令行选项；
3. 环境变量；
4. `Default.Route` 文件（先在当前目录查找，再在 `GAUSS_EXEDIR` 的目录中查找）。

| 含义 | Link 0 | 命令行 | 环境变量 | Default.Route |
|---|---|---|---|---|
| 内存 | `%Mem` | `-m` | `GAUSS_MDEF` | `-M-` |
| 核列表（推荐） | `%CPU` | `-c` | `GAUSS_CDEF` | `-C-` |
| 核数（已弃用，仍可用） | `%NProcShared` | `-p` | `GAUSS_PDEF` | `-P-` |
| GPU 与控制核 | `%GPUCPU` | `-g` | `GAUSS_GDEF` | `-G-` |
| Linda worker 列表 | `%LindaWorkers` | `-w` | `GAUSS_WDEF` | `-W-` |
| 启动 worker 的命令 | `%SSH=ssh` / `%UseSSH` | `-s` | `GAUSS_SDEF` | `-S-` |
| route 默认关键词 | — | `-r` | `GAUSS_RDEF` | `-R-` |
| 检查点文件 | `%Chk` | `-y` | `GAUSS_YDEF` | — |
| read-write 文件 | `%RWF` | `-z` | `GAUSS_ZDEF` | — |
| formchk 默认选项 | — | — | `GAUSS_FDEF` | `-F-` |

命令行和环境变量这两种方式是 Rev. B.01 引入的。例如下面三种写法等价：

```bash
# 1) Link 0 lines inside job.gjf:   %CPU=0-15   %Mem=48GB
# 2) command line
g16 -c="0-15" -m=48GB job.gjf
# 3) environment variables (typical in job scripts)
export GAUSS_CDEF="0-15"
export GAUSS_MDEF="48GB"
g16 < job.gjf > job.log
```

`%CPU` 与 `%NProcShared` 的区别：`%CPU=0-15` 把 16 个线程**绑定**到 0–15 号核上，减少线程在核之间迁移造成的缓存失效，官方推荐用它取代 `%NProcShared`；`%NProcShared=16` 只指定线程数，由操作系统决定放在哪些核上。`%CPU` 还支持步长写法，例如 `%CPU=0-47/2` 表示使用 0, 2, 4, … 号核。

> [!IMPORTANT]
> **推荐策略**：输入文件里不写 `%nprocshared`/`%cpu`/`%mem`，由作业脚本根据调度系统实际分配的资源设置 `GAUSS_CDEF`（或 `GAUSS_PDEF`）和 `GAUSS_MDEF`。这样同一个 .gjf 拿到任何机器、任何队列都能直接用，资源只需在 `#SBATCH` 行改一处。反过来，如果输入里保留了 Link 0 资源行，它们会**覆盖**环境变量，这是“明明申请了 32 核却只跑 8 核”或“申请 16 核却跑了 64 线程”的常见原因。

### 3.1 Default.Route

`Default.Route` 是一个文本文件，为所有作业设置默认值。放在作业目录中只影响该目录的作业，放在 `$g16root/g16` 中影响所有用户（通常由管理员维护）。示例（`Default.Route.example`）：

```text
-M- 4GB
-R- MaxDisk=100GB
```

- `-M- 4GB`：把默认内存从 800 MB 提高到 4 GB，防止忘写 `%mem` 的作业因内存过小而极慢或失败；
- `-R- MaxDisk=100GB`：route 默认关键词，这里声明 scratch 可用 100 GB（MP2 等方法会参考它安排算法）。

旧版本文档中常见 `-#-` 写法表示 route 默认值，G16 官方等价表中写作 `-R-`。不建议在共享节点的 `Default.Route` 中写 `-C-`，因为调度系统每次分配的核编号可能不同。

---

## 4. Slurm 作业脚本

### 4.1 最小可用脚本

```bash
#!/bin/bash
#SBATCH --job-name=h2o
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out

source $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=/tmp/$USER/g16_$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR
cd $SLURM_SUBMIT_DIR
g16 < h2o.gjf > h2o.log
```

配套的输入里写 `%nprocshared=8`、`%mem=12GB`（16 GB 的 75%）。

| 行 | 含义 | 为什么这样写 |
|---|---|---|
| `--nodes=1` | 一个节点 | SMP 并行不能跨节点 |
| `--ntasks=1` | 一个任务（进程） | Gaussian 是一个主进程加多个线程，不是 MPI 程序 |
| `--cpus-per-task=8` | 这个任务用 8 个核 | 必须与 `%nprocshared`/`%CPU` 一致 |
| `--mem=16G` | 节点上的内存上限 | Gaussian 的 `%mem` 取它的 75–85% |
| `--time` | 最长运行时间 | 超时会被杀掉；chk 中保留了已完成的优化步，可续算 |
| `--output=%x-%j.out` | 调度系统的标准输出，`%x` 为作业名，`%j` 为作业号 | 与 Gaussian 的 .log 分开 |

> [!WARNING]
> 不要写成 `--ntasks=8`。那样 Slurm 认为你要 8 个单核进程，在启用了 cgroup 的集群上，Gaussian 的 8 个线程可能全挤在 1 个核上，速度慢几倍。

### 4.2 推荐模板 `g16_slurm.sh`

示例目录中的 [`g16_slurm.sh`](../examples/17-hpc-job-scripts/g16_slurm.sh) 在最小脚本的基础上加入了以下功能：

1. 从 `/proc/self/status` 读取 Slurm 实际分配给作业的核列表（`Cpus_allowed_list`），设置 `GAUSS_CDEF`，把线程精确绑定到这些核上；如果列表中的核数与 `--cpus-per-task` 不一致（例如 cgroup 中包含超线程），退回到 `GAUSS_PDEF`；
2. 按 `SLURM_MEM_PER_NODE`（或 `SLURM_MEM_PER_CPU` × 核数）的 80% 设置 `GAUSS_MDEF`；
3. 检查输入：CRLF 换行自动转换并备份，文件末尾缺空行时补上，发现 Link 0 资源行时打印警告；
4. 为每个作业建立独立的 scratch 目录，退出时自动删除（可设 `KEEP_SCR_ON_FAIL=1` 让失败作业保留 scratch）；
5. Gaussian 在后台运行并 `wait`，这样收到 `scancel` 或超时的 SIGTERM 时，清理代码能立即执行；
6. 结束后检查 `Normal termination`，自动对输入中所有 `%chk` 文件做 formchk；失败时打印出错的 Link 附近的内容，并以非零状态退出（依赖作业 `afterok` 需要这一点）。

核心部分如下（完整脚本见示例目录）：

```bash
# cores: pin threads to exactly the CPUs Slurm allocated
NCPU=${SLURM_CPUS_PER_TASK:-1}
CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
NLIST=$(expand_cpus "$CPULIST" | wc -l)          # "0-3,8" -> 5
if [ -n "$CPULIST" ] && [ "$NLIST" -eq "$NCPU" ]; then
    export GAUSS_CDEF=$CPULIST
else
    unset GAUSS_CDEF; export GAUSS_PDEF=$NCPU
fi

# memory: ~80% of what Slurm granted (SLURM_MEM_PER_NODE is in MB)
MEM_MB=${SLURM_MEM_PER_NODE:-4096}
export GAUSS_MDEF="$(( MEM_MB * 80 / 100 ))MB"

# private scratch, removed at exit
export GAUSS_SCRDIR="$SCR_BASE/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR"
trap cleanup EXIT

g16 < "$INPUT" > "$LOG" &
G16PID=$!
wait "$G16PID"
```

使用：

```bash
sbatch g16_slurm.sh benzene_opt_freq.gjf
sbatch --cpus-per-task=32 --mem=120G --time=48:00:00 g16_slurm.sh big.gjf   # override resources
SCR_BASE=/scratch/$USER sbatch g16_slurm.sh job.gjf                          # another scratch disk
```

命令行上的 `sbatch` 选项优先于脚本中的 `#SBATCH` 行，所以一个模板可以适配不同大小的作业。

> [!TIP]
> 为什么要绑定到 `Cpus_allowed_list` 而不是直接写 `%CPU=0-15`？在多个作业共享一个节点时，Slurm 可能把 16–31 号核分给你的作业。如果输入中写死 `%CPU=0-15`，Gaussian 会试图把线程绑到不属于本作业的核上。读取 `Cpus_allowed_list` 可以保证列表与分配一致。如果整节点独占（`--exclusive`），`0-(N-1)` 也是安全的。

### 4.3 提交与管理

```bash
sbatch g16_slurm.sh job.gjf          # submit; prints "Submitted batch job 123456"
squeue -u $USER                      # my jobs: PD = pending, R = running
scontrol show job 123456             # details: nodes, resources, reason for pending
scancel 123456                       # cancel
sacct -j 123456 --format=JobID,JobName,Elapsed,MaxRSS,State,ExitCode   # after the job
seff 123456                          # CPU and memory efficiency (if installed)
```

`seff` 会给出 CPU 效率和内存效率。CPU 效率远低于 100%，说明核申请多了或线程没有用起来（检查 Link 0 是否覆盖了环境变量）；内存效率很低，说明 `--mem` 可以减少，让作业更容易排上队。

### 4.4 超线程

Gaussian 官方建议关闭超线程；不能关闭时，每个物理核只用一个线程。Slurm 中可以加 `#SBATCH --hint=nomultithread`（是否生效取决于集群配置）。在作业中运行 `lscpu | grep -i "thread"` 可以查看每个核的线程数。

---

## 5. PBS/Torque、PBS Pro 与 LSF

### 5.1 PBS 模板

示例目录中的 [`g16_pbs.sh`](../examples/17-hpc-job-scripts/g16_pbs.sh) 核心部分：

```bash
#PBS -N g16
#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
source $g16root/g16/bsd/g16.profile
export GAUSS_PDEF=${PBS_NUM_PPN:-${NCPUS:-16}}
export GAUSS_MDEF=52428MB                       # ~80% of 64 GB
export GAUSS_SCRDIR=/tmp/$USER/g16_${PBS_JOBID%%.*}
mkdir -p $GAUSS_SCRDIR
g16 < $INPUT > ${INPUT%.*}.log
```

```bash
qsub -v INPUT=job.gjf g16_pbs.sh     # pass the input name as a variable
qstat -u $USER
qdel <jobid>
```

| 项目 | Torque | PBS Pro |
|---|---|---|
| 申请 1 节点 16 核 | `#PBS -l nodes=1:ppn=16` | `#PBS -l select=1:ncpus=16:mem=64gb` |
| 作业中的核数变量 | `PBS_NUM_PPN` | `NCPUS` |
| 提交目录 | `PBS_O_WORKDIR` | `PBS_O_WORKDIR` |
| 依赖作业 | `qsub -W depend=afterok:<jobid>` | 同左 |

### 5.2 LSF 简要对照

```bash
#BSUB -J g16
#BSUB -n 16
#BSUB -R "span[hosts=1]"          # all 16 slots on one host (required for SMP)
#BSUB -W 24:00
#BSUB -o %J.out
export GAUSS_PDEF=$LSB_DJOB_NUMPROC
```

LSF 的内存申请写法（`rusage[mem=...]` 的单位、按核还是按作业计）因站点配置而异，请查阅本集群文档。更多调度系统对照见 Common 第2章。

---

## 6. 核数与内存怎么配

### 6.1 基本规则

| 规则 | 说明 |
|---|---|
| `%nprocshared`（或 `%CPU` 列表的核数）= `--cpus-per-task` | 多了会超额占用，少了浪费 |
| `%mem` ≈ `--mem` 的 75–85% | 给 Gaussian 可执行文件本身、I/O 缓冲和系统留出余量；`%mem` 只是 Gaussian 的工作数组，进程实际占用的内存会比它大 |
| 内存随核数增加 | 约每核 2–4 GB 起步，大体系每核 4 GB 以上 |
| 大频率、CCSD 留出磁盘缓存 | 独占整节点时只给 Gaussian 整机内存的 50–70% |
| 多个小作业优于一个大作业 | 吞吐量优先时，16 个 8 核作业通常比 4 个 32 核作业完成得快 |

### 6.2 经验配置（起点，需按体系调整）

| 体系与任务 | 核数 | `--mem` | `%mem` |
|---|---|---|---|
| ≤ 30 原子，def2-SVP 优化 + 频率 | 8 | 16 GB | 12–13 GB |
| 30–80 原子，def2-SVP 优化 + 频率 | 16 | 32–64 GB | 26–52 GB |
| 30–80 原子，def2-TZVP 单点 | 16–32 | 64 GB | 52 GB |
| 100+ 原子，频率或 TD-DFT 优化 | 32–64 | 128–256 GB | 100–200 GB |
| 双杂化（DSD-PBEP86）单点，大基组 | 16–32 | 128 GB 以上 | 100 GB 以上，并确保 scratch 充足 |

表中的数字是常见的起始配置，不是测试结果。真正的最优值要用第 13 节的缩放测试确定。频率计算的内存需求可以用 Gaussian 自带的 `freqmem` 工具估计。

### 6.3 内存不当的典型表现

| 现象 | 原因 | 处理 |
|---|---|---|
| 作业一开始就报 `galloc: could not allocate memory.` | `%mem` 超过了节点当前可用的内存 | 降低 `%mem`，或申请更多内存 |
| 作业被调度系统杀掉，`sacct` 显示 `OUT_OF_MEMORY` 或类似状态 | 进程实际内存超出 `--mem`（`%mem` 设得太接近 `--mem`） | `%mem` 改为 `--mem` 的 75–80% |
| 程序提示需要更多内存的 words，或自动降低并行线程数 | `%mem` 太小（忘写时只有 800 MB） | 增加 `%mem`；确认没有被 Link 0 覆盖 |
| 作业不报错但极慢，`top` 中可见大量交换（swap） | 请求的内存超过物理内存 | 降低 `%mem` |

---

## 7. Scratch 目录与磁盘

### 7.1 每个作业一个 scratch 目录

```bash
export GAUSS_SCRDIR="/scratch/$USER/g16_$SLURM_JOB_ID"   # or a node-local disk such as /tmp
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT
```

要点：

- 优先使用**节点本地磁盘**（`/tmp`、`/local`、`/scratch` 等，以集群文档为准），速度比共享文件系统快，也不会拖慢其他用户；
- 目录名带作业号，同一节点上的多个作业互不干扰；
- `GAUSS_SCRDIR` 指向的目录**必须存在**，否则作业会立即失败；
- `trap ... EXIT` 保证正常结束、出错退出、被 `scancel` 时都会清理；`${GAUSS_SCRDIR:?}` 防止变量为空时误删其他目录；
- 作业开始前用 `df -h "$GAUSS_SCRDIR"` 检查剩余空间。大体系频率计算的 rwf 可以有几十到上百 GB。

### 7.2 chk 放在哪里

`%chk=job.chk` 写相对路径时，chk 位于作业的工作目录（通常是共享文件系统上的提交目录）。这是推荐做法：节点故障或作业被杀时，chk 仍然保留，可以续算。不要把唯一的 chk 放在会被清理的本地 scratch 中。

### 7.3 用 %RWF 拆分大文件

```text
%RWF=/local1/g16/,60GB,/local2/g16/,-1
%NoSave
%Chk=bigfreq.chk
```

- 每个位置后跟该段文件的最大大小，`-1` 表示可用剩余全部空间（只应用于最后一段）；只写目录时必须以 `/` 结尾，Gaussian 会自动生成文件名；
- `%NoSave` 之前命名的文件在正常结束时会被删除，之后命名的文件（这里是 chk）会保留；异常结束时所有文件都保留，便于续算；
- `MaxDisk=` 关键词告诉 Gaussian 可用的 scratch 大小，MP2 等方法据此选择算法；SCF 能量、梯度和频率的磁盘用量是固定的，不受 `MaxDisk` 控制。

---

## 8. Linda 跨节点并行

### 8.1 什么时候值得用

- 单节点的核数和内存已经不够，而计算是 HF/DFT 的能量、梯度、频率，或 MP2 能量与梯度，或 TD-DFT 能量与梯度；
- 你的 Gaussian 是 Linda 版本（需要单独授权；C.01 起需要 Linda 9.2）；
- CCSD、EOM-CCSD 的能量和优化只支持 SMP，不能跨节点；数值导数类作业（如 DFT 非谐频率 `Freq=Anharm`）会把不同方向的导数分给不同节点。

### 8.2 前提条件

1. 所有节点上 `$g16root/g16` 路径相同、都能访问；
2. 节点之间可以免密码 ssh（Linda 默认用 rsh，用 `%UseSSH` 或 `GAUSS_SDEF=ssh` 改为 ssh）；
3. 每个节点都有本地 scratch（`GAUSS_SCRDIR` 在每个节点上都要存在）；
4. 已执行 `g16.profile`，`LD_LIBRARY_PATH` 能找到 Linda 的库。

### 8.3 写法

```text
%Mem=96GB
%CPU=0-31
%LindaWorkers=node01,node02
#p B3LYP/def2TZVP Freq
```

- `%CPU` 和 `%Mem` **对每个节点分别生效**：上面表示两个节点各开 32 个线程、各用 96 GB；
- 主节点（启动 g16 的节点）总会启动一个 worker，无论它是否出现在列表中；
- 数值导数类作业（`Freq=Anharm`、CCSD 频率等）需要在主节点上多开一个收集结果的 worker，例如 `%LindaWorkers=node01:2,node02`；
- 不要使用已废弃的 `%NProcLinda`。

在 Slurm 中，节点名在运行时才知道，所以用环境变量传入（示例 [`g16_linda_slurm.sh`](../examples/17-hpc-job-scripts/g16_linda_slurm.sh)）：

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive

export GAUSS_WDEF=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export GAUSS_SDEF=ssh
export GAUSS_CDEF="0-31"
export GAUSS_MDEF="96GB"
srun --ntasks-per-node=1 mkdir -p "$GAUSS_SCRDIR"     # scratch on every node
g16 < job.gjf > job.log
```

worker 启动不了时，`export GAUSS_LFLAGS="-vv"`（或 Link 0 中加 `%DebugLinda`）打印详细的 Linda 信息。运行中在 worker 节点上用 `top` 应能看到 `lNNN.exel` 形式的进程。

---

## 9. GPU 计算

### 9.1 支持范围

| 项目 | 说明 |
|---|---|
| 平台 | 仅 Linux |
| GPU 型号 | NVIDIA K40、K80、P100（B.01 起）、V100（C.01 起）、A100（**C.02 起**） |
| 方法 | HF 与 DFT 的能量、梯度、频率（基态和激发态） |
| 效果 | 较大的分子才明显；小作业、MP2、CCSD 等后 SCF 方法基本无收益 |

### 9.2 %CPU 与 %GPUCPU

```text
%CPU=0-31
%GPUCPU=0-3=0,1,16,17
%Mem=384GB
```

- `%GPUCPU=GPU列表=控制核列表`：这里用 GPU 0–3，分别由 0、1、16、17 号核控制；也可写成 `%GPUCPU=0,1,2,3=0,1,16,17`；
- 控制核**包含在** `%CPU` 中，但在 GPU 并行的部分只负责驱动 GPU，不参与计算；
- 每个 GPU 必须有自己的控制核，不能共用；控制核应位于 GPU 所连接的那颗 CPU 上，用 `nvidia-smi topo -m` 的 “CPU Affinity” 列确认；
- 内存：Gaussian 把 `%mem` 平均分给每个线程，而控制 GPU 的线程需要与 GPU 计算量相当的主机内存。官方给出的经验是 16 GB 显存的 GPU 给 Gaussian 11–12 GB，因此总内存应为“线程数 × 每个 GPU 需要的内存”，例如 4 个 CPU 线程、2 个 16 GB GPU 时用 4 × 12 GB；
- 估算收益：控制核不参与计算。官方的例子是 32 核 + 8 GPU、每个 GPU 相当于 5 个核时，理论加速约 2 倍（不是 8 × 5 倍）。

### 9.3 Slurm 中的 GPU 脚本

示例 [`g16_gpu_slurm.sh`](../examples/17-hpc-job-scripts/g16_gpu_slurm.sh) 申请 `--gres=gpu:2`，读取分配到的核列表作为 `GAUSS_CDEF`，用前两个核作为控制核组成 `GAUSS_GDEF`，并在日志中打印 `nvidia-smi -L` 和 `nvidia-smi topo -m`，方便核对控制核是否合理。Slurm 一般通过 `CUDA_VISIBLE_DEVICES` 只暴露分配给作业的 GPU，作业内看到的编号通常从 0 开始；请以作业中 `nvidia-smi` 的实际输出为准。运行中可以在节点上执行 `nvidia-smi` 查看 GPU 利用率，确认 Gaussian 确实在使用 GPU。

---

## 10. 批量提交与依赖作业

### 10.1 循环提交整个目录

[`batch_submit_all_gjf.sh`](../examples/17-hpc-job-scripts/batch_submit_all_gjf.sh) 对目录中每个 .gjf 调用 `g16_slurm.sh`，并且：

- 跳过已经 `Normal termination` 的作业（可以反复运行，只补交失败或未完成的）；
- 跳过已在队列中的同名作业；
- 队列中作业数达到上限时等待（许多集群限制单用户排队作业数）。

```bash
bash batch_submit_all_gjf.sh conformers/ 40
SBATCH_OPTS="--cpus-per-task=8 --mem=16G" bash batch_submit_all_gjf.sh conformers/
```

### 10.2 作业数组

上百个同类小作业用作业数组更省事，调度器只看到一个作业：

```bash
ls conf_*.gjf > joblist.txt
sbatch --array=1-$(wc -l < joblist.txt)%10 g16_array_slurm.sh   # at most 10 running at once
```

[`g16_array_slurm.sh`](../examples/17-hpc-job-scripts/g16_array_slurm.sh) 取出 `joblist.txt` 的第 `$SLURM_ARRAY_TASK_ID` 行交给 `g16_slurm.sh`。数组中每个任务有自己的 `SLURM_JOB_ID`，scratch 目录不会冲突。

### 10.3 依赖作业

优化完成后自动开始高精度单点：

```bash
jid1=$(sbatch --parsable g16_slurm.sh benzene_opt_freq.gjf)
sbatch --dependency=afterok:$jid1 g16_slurm.sh benzene_sp_tz.gjf
```

`afterok` 表示前一个作业**以退出码 0 结束**才开始。所以作业脚本必须在 Gaussian 失败时返回非零退出码，`g16_slurm.sh` 已经这样做了；如果你的脚本最后一行是 `echo done` 之类的命令，即使 Gaussian 报错，作业也会被认为成功。完整示例见 [`chain_opt_sp_dependency.sh`](../examples/17-hpc-job-scripts/chain_opt_sp_dependency.sh)。PBS 中写作 `qsub -W depend=afterok:<jobid> ...`。

### 10.4 一个节点上同时跑多个小作业

有的集群按整节点计费。这时可以在一个作业中并行运行几个小 Gaussian 作业，每个绑定不同的核：

```bash
GAUSS_CDEF=0-7   GAUSS_MDEF=14GB g16 < a.gjf > a.log &
GAUSS_CDEF=8-15  GAUSS_MDEF=14GB g16 < b.gjf > b.log &
GAUSS_CDEF=16-23 GAUSS_MDEF=14GB g16 < c.gjf > c.log &
GAUSS_CDEF=24-31 GAUSS_MDEF=14GB g16 < d.gjf > d.log &
wait
```

各作业的内存之和要小于节点内存；它们共用同一个 `GAUSS_SCRDIR` 没有问题（文件名带进程号）。

### 10.5 用脚本生成输入

批量构象、批量底物的输入不要手工编辑。仓库脚本 [`gjf_builder.py`](../scripts/gjf_builder.py) 可以从 xyz/log/gjf 按任务模板生成输入（`python ../scripts/gjf_builder.py -h` 查看用法）；`sub_g16.sh` 是带自动内存换算、scratch 管理和批量提交功能的提交封装（`bash ../scripts/sub_g16.sh -h`）。用 ASE/pymatgen 驱动 Gaussian 见 Common 第4章，高通量工作流（AiiDA 等）见 Common 第6章。

---

## 11. 作业监控

```bash
tail -f job.log                                   # follow the output (Ctrl-C to stop)
grep "SCF Done" job.log | tail -n 5               # last SCF energies
grep -c "Step number" job.log                     # optimization steps so far
grep -A 5 "Converged?" job.log | tail -n 6        # latest convergence table
grep -c "Normal termination" job.log              # finished job steps
grep "Error termination" job.log                  # failure, and in which link
grep "Elapsed time" job.log                       # wall time of each finished step
```

优化过程中 `Converged?` 表列出 4 个判据（Maximum/RMS Force、Maximum/RMS Displacement）及阈值，全部为 `YES` 即收敛（[第5章](05-opt.md)）。用仓库脚本 [`opt_monitor.py`](../scripts/opt_monitor.py) 可以把每一步的能量和 4 个判据画成曲线，判断是振荡还是缓慢下降；[`g16_log_summary.py`](../scripts/g16_log_summary.py) 可以把一批 log 的终止状态、能量、虚频、热化学汇总成表格（都用 `-h` 查看用法）。

确认并行是否生效：

```bash
ssh <node>            # if your site allows logging in to nodes running your jobs
top -u $USER          # the lNNN.exe process should show ~N x 100% CPU for N threads
```

如果 16 核作业的 CPU 占用只有 100% 左右，通常是输入中的 `%nprocshared=1`（或遗漏）覆盖了环境变量，或者当前正在运行的 Link 本身是串行的（某些 Link 短时间内只用一个核是正常的）。

---

## 12. 断点续算

### 12.1 几何优化：Opt=Restart

优化因超时、节点故障或达到最大步数而中断时，chk 中保存了最后的结构和 Hessian 信息：

```text
%chk=benzene_opt_freq.chk
#p Opt=(Restart,MaxCycles=200)

```

- route 中**只写 Opt 关键词及其选项**（加上 `Restart`），方法、基组和结构都从 chk 读取，不需要标题和坐标；
- 原作业的选项（如 `TS`）必须再写一遍；可以修改 `MaxCycles` 等参数；
- 想在续算的第一步重新计算 Hessian，必须写 `CalcFC`（不论原作业是否写过）；
- 原作业是 `Opt Freq` 时，续算的只是优化部分，频率要另写一步（示例 `benzene_opt_restart.gjf` 用 `--Link1--` 接了一个 `Freq Geom=AllCheck Guess=Read`）。

另一种常用做法是从最后的结构重新开始一个新优化：`Geom=Check Guess=Read`（读 chk 中最后的结构和轨道），或用 [`log2xyz.py`](../scripts/log2xyz.py) 从 log 中取出最后一步或能量最低的一步。如果优化是因 L9999（步数用完）而停止，先用 [`opt_monitor.py`](../scripts/opt_monitor.py) 判断原因，再决定是续算还是改设置（[第5章](05-opt.md)、[第18章](18-errors-capstone.md)）。

### 12.2 IRC：IRC=Restart

```text
%chk=irc.chk
#p IRC=(Restart,MaxPoints=60)

```

`IRC=Restart` 既可以续算没有完成的 IRC，也可以为已经完成的 IRC 增加路径点（例如加大 `MaxPoints`）。

### 12.3 频率：Restart 关键词与 %RWF

解析频率的中间数据很大，不存放在 chk 中，所以续算依赖**命名并保留的 read-write 文件**。原作业这样写（`benzene_freq_rwf.gjf`）：

```text
%rwf=benzene_freq.rwf
%nosave
%chk=benzene_freq.chk
#p B3LYP EmpiricalDispersion=GD3BJ def2TZVP Freq
...
```

作业中断后，续算输入（`benzene_freq_restart.gjf`）只需要：

```text
%rwf=benzene_freq.rwf
%nosave
%chk=benzene_freq.chk
#p Restart

```

- `%nosave` 写在 `%rwf` 之后、`%chk` 之前：正常结束时删除 rwf、保留 chk；异常结束时全部保留；
- rwf 可能非常大，应放在空间充足、作业结束后**不会被自动清理**的位置。用 `g16_slurm.sh` 时注意：它会在退出时清理 `GAUSS_SCRDIR`，所以不要把要续算的 rwf 放在那里；
- `Restart` 关键词适用于解析频率（包括 ONIOM 的 ROA、VCD）、CCSD/EOM-CCSD、NMR、`Polar=OptRot` 等长时间作业；
- 数值频率（`Freq=Numer`）则直接从 chk 续算：重复原 route，加 `Freq=(Numer,Restart)`；
- 官方特别指出：可续算的是“执行中被打断”的作业；因收敛失败或步数用完而报错退出的作业需要人工处理，不能简单续算。

### 12.4 续算清单

| 作业类型 | 续算方法 | 需要的文件 |
|---|---|---|
| 几何优化、柔性扫描、TS 优化 | `Opt=(Restart,...)`，保留原选项 | chk |
| IRC | `IRC=(Restart,...)` | chk |
| 解析频率、NMR、CCSD 等 | `%RWF` + `%NoSave` + `#p Restart` | rwf + chk |
| 数值频率 | `Freq=(Numer,Restart)` | chk |
| 单点能 | 直接重算，`Guess=Read` 读入已有轨道可加快 SCF | chk |

---

## 13. 并行效率测试：选择合适的核数

不同体系、不同方法的最佳核数差别很大，最可靠的办法是自己测。示例中的 [`scaling_test.sh`](../examples/17-hpc-job-scripts/scaling_test.sh) 在**同一次分配**中依次用 1、2、4、8、16、32 个核运行同一个输入，记录墙钟时间，计算加速比和并行效率：

$$
S(N)=\frac{T(1)}{T(N)},\qquad E(N)=\frac{S(N)}{N}
$$

```bash
sbatch scaling_test.sh benzene_scaling_freq.gjf "1 2 4 8 16 32"
cat scaling_benzene_scaling_freq.dat
```

输出文件的格式（示意输出，数值仅为占位）：

```text
 cores     wall_s  speedup efficiency
     1      ...       1.00       1.00
     2      ...        ...        ...
   ...
```

判断标准：

- 并行效率 E(N) 降到 0.6–0.7 以下时，再加核基本是浪费，排队时间还会变长；
- 用同一套测试换一个更大的体系（例如 def2-TZVP、上百个原子），通常会看到效率随体系增大而提高；
- 测试时内存保持不变（脚本中 `GAUSS_MDEF` 固定），否则无法区分“核多了”和“内存多了”的影响；
- 苯这样的小分子正是用来说明问题的：它在少量核之后就几乎不再加速。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 说明 |
|---|---|---|
| 调度申请 | `--nodes=1 --ntasks=1 --cpus-per-task=N` | SMP 作业的标准写法 |
| 核数 | 小分子 8–16；大分子按缩放测试 | 吞吐量优先时宁可多个小作业 |
| 内存 | `%mem` = `--mem` 的 75–85%；每核 2–4 GB 起 | 大频率/CCSD 独占节点时给 50–70% |
| 资源写法 | 输入不写 `%nprocshared`/`%mem`，脚本设 `GAUSS_CDEF`/`GAUSS_MDEF` | 或者两边保持严格一致 |
| 绑核 | `GAUSS_CDEF` 取 `Cpus_allowed_list` | 共享节点上不要写死 `0-(N-1)` |
| scratch | 节点本地盘，每作业独立目录，`trap` 清理 | 先 `df -h` 检查空间 |
| chk | 放在提交目录（共享文件系统） | 节点故障后仍能续算 |
| 依赖作业 | `--dependency=afterok`，脚本失败时非零退出 | 否则失败作业也会触发后续作业 |
| 多节点 | 节点内 SMP + 节点间 Linda | 需要 Linda 授权和免密 ssh |
| GPU | 大分子 HF/DFT；`%CPU` 包含控制核；内存按线程数 × 单 GPU 需求 | A100 需要 C.02 |

---

## 常见问题与报错

**1. 申请了 32 核，作业只用 1 个核（CPU 占用约 100%）**

- 原因：输入中的 `%nprocshared=1`、`%nprocshared=4` 等覆盖了 `GAUSS_CDEF`/`GAUSS_PDEF`；或者既没写 Link 0 也没设环境变量（默认单核）。
- 处理：删掉输入中的资源行，让脚本统一设置；或使两者一致。`g16_slurm.sh` 会对发现的资源行打印警告。

**2. 作业比在工作站上还慢，节点负载远超核数**

- 原因：`%nprocshared` 大于 `--cpus-per-task`（例如 64 对 16），多个线程挤在少数核上；或者用了 `--ntasks=N` 而不是 `--cpus-per-task=N`。
- 处理：两者保持一致；按 4.1 节的写法申请资源。

**3. `galloc: could not allocate memory.`**

- 原因：`%mem` 超过节点可用内存。
- 处理：`%mem` 设为 `--mem` 的 75–85%；检查同一节点上自己的其他作业是否也占用了大量内存。

**4. 作业被调度系统杀掉，log 突然中断，没有 `Error termination`**

- 原因：超出 `--mem`（内存溢出被杀）或超出 `--time`。
- 处理：`sacct -j <jobid> --format=JobID,State,MaxRSS,Elapsed` 查看状态；内存问题降低 `%mem` 或增加 `--mem`；超时则用第 12 节的方法续算，并申请更长时间。

**5. `Erroneous write. Write -1 instead of ...`**

- 原因：scratch 磁盘写满或超出配额。
- 处理：`df -h $GAUSS_SCRDIR` 检查；换到更大的磁盘；清理以前异常结束遗留的 `Gau-*` 文件；用 `%RWF` 拆分到多个磁盘。

**6. 作业一启动就失败，提示无法打开 scratch 文件**

- 原因：`GAUSS_SCRDIR` 指向的目录不存在或没有写权限（例如 Linda 作业只在主节点上建了目录）。
- 处理：作业开始时 `mkdir -p "$GAUSS_SCRDIR"`；多节点作业用 `srun --ntasks-per-node=1 mkdir -p ...` 在每个节点上建立。

**7. `g16: command not found`，或提示找不到 Link 1（`l1.exe`）**

- 原因：作业脚本中没有 `source $g16root/g16/bsd/g16.profile`（计算节点不会继承登录节点的环境），或 `g16root`/`GAUSS_EXEDIR` 设置错误。
- 处理：在脚本中显式 `module load` 或设置 `g16root` 并 source profile。

**8. `sbatch: error: Batch script contains DOS line breaks (\r\n)` 或 `/bin/bash^M: bad interpreter`**

- 原因：脚本在 Windows 下编辑，带有 CRLF 换行。
- 处理：`dos2unix g16_slurm.sh` 或 `sed -i 's/\r$//' g16_slurm.sh`。.gjf 同理（`g16_slurm.sh` 会自动转换输入文件）。

**9. 依赖作业在前一个作业失败后仍然运行**

- 原因：前一个作业脚本以退出码 0 结束（最后一条命令成功）。
- 处理：在脚本中检查 `Normal termination`，失败时 `exit 1`。

**10. Linda 作业只在主节点上计算**

- 原因：不是 Linda 版本；ssh 需要密码；`GAUSS_WDEF`/`%LindaWorkers` 为空或节点名不对。
- 处理：确认安装版本；测试 `ssh node02 hostname` 无需密码；设置 `GAUSS_LFLAGS="-vv"` 查看 worker 启动信息。

**11. GPU 作业没有加速，甚至更慢**

- 原因：体系太小；控制核选得离 GPU 太远；`%mem` 不足；或方法本身（MP2、CCSD）不使用 GPU。
- 处理：只对较大的 HF/DFT 作业使用 GPU；按 `nvidia-smi topo -m` 选择控制核；按 9.2 节估算内存。

---

## 练习

1. 用 `g16_slurm.sh` 提交 `benzene_opt_freq.gjf`，在 `.out` 文件中找到 `GAUSS_CDEF` 和 `GAUSS_MDEF` 的实际取值，并用 `seff` 查看 CPU 效率。
2. 在 `benzene_opt_freq.gjf` 中故意加入 `%nprocshared=2`，再用 16 核提交，比较墙钟时间和 `seff` 的 CPU 效率。
3. 用 `scaling_test.sh` 分别测试苯（def2-TZVP 频率）和一个 50 原子以上分子的缩放性，画出加速比曲线并比较。
4. 写一个 `joblist.txt`，用作业数组提交 10 个构象的优化，限制同时运行 3 个。
5. 提交 `benzene_opt_freq.gjf` 时设置很短的 `--time`（如 2 分钟）让它超时，然后用 `benzene_opt_restart.gjf` 续算完成。
6. 用 `chain_opt_sp_dependency.sh` 串联优化和单点；再把优化输入故意改错（例如电荷与多重度不匹配），观察单点作业是否会被取消。

---

## 参考资料

**Gaussian 官方文档**

- [Running Gaussian 16](https://gaussian.com/running/)：环境变量、scratch 文件、`%RWF` 拆分、内存、SMP 与 Linda 并行
- [Link 0 Commands](https://gaussian.com/link0/)：`%Mem`、`%CPU`、`%NProcShared`、`%GPUCPU`、`%LindaWorkers`、`%RWF`、`%NoSave`，以及四种设置方式的等价表
- [g16 Command Line Options](https://gaussian.com/options/)：`-m`、`-c`、`-g`、`-w`、`-r`、`-fchk` 等
- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)：GPU（V100 自 C.01、A100 自 C.02）、Parallel Performance 说明、内存建议
- [Opt 关键词（Restart）](https://gaussian.com/opt/)、[IRC 关键词（Restart）](https://gaussian.com/irc/)、[Freq 关键词](https://gaussian.com/freq/)
- [Restart 关键词](https://gaussian.com/restart/)、[MaxDisk 关键词](https://gaussian.com/maxdisk/)

**调度系统文档**

- [Slurm sbatch](https://slurm.schedmd.com/sbatch.html)、[Slurm job arrays](https://slurm.schedmd.com/job_array.html)
- [OpenPBS 文档](https://openpbs.org/)

**相关 GitHub 仓库**（仅链接，未复制代码）

- [smleechem/gaussjob_slurm](https://github.com/smleechem/gaussjob_slurm)：每个作业独立的 `GAUSS_SCRDIR`，捕获信号后自动清理 scratch
- [MauricioCafiero/Gaussian16-Scripts-and-Functions](https://github.com/MauricioCafiero/Gaussian16-Scripts-and-Functions)：自动同步输入中的 `%mem`/`%nprocshared` 与作业请求，并提供批量提交（CC0）
- [TonksGroupUMN/Tonks-Gaussian-Scripts](https://github.com/TonksGroupUMN/Tonks-Gaussian-Scripts)：Slurm 与 PBS 两种提交脚本
- [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash)：成熟的 g16 bash 工具集（排队提交、testrt 检查 route 等）
- [njzjz/gaussianrunner](https://github.com/njzjz/gaussianrunner)：从 SMILES 批量生成并并行运行 Gaussian
- [lenhanpham/Gaussian-Extractor](https://github.com/lenhanpham/Gaussian-Extractor)：并行批量提取 log 数据、检查作业状态
- [nanotech-empa/aiida-gaussian](https://github.com/nanotech-empa/aiida-gaussian)：在 AiiDA 中运行 Gaussian 并记录数据溯源
- [thisisntnathan/dftCourse](https://github.com/thisisntnathan/dftCourse)：Cornell 的 Gaussian + SLURM 入门课程

**延伸阅读**

- G. M. Amdahl, Validity of the single processor approach to achieving large scale computing capabilities, *AFIPS Conf. Proc.* **1967**, 30, 483. DOI: [10.1145/1465482.1465560](https://doi.org/10.1145/1465482.1465560)

---

[← 上一章](16-postprocess-multiwfn.md) | [返回目录](../README.md) | [下一章 →](18-errors-capstone.md)
