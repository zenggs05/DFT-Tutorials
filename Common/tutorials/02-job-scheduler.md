# 第2章 HPC 作业调度：Slurm/PBS/LSF 模板与并行设置（HPC Job Scheduling: Slurm/PBS/LSF Templates and Parallel Settings）

本章讲解如何把 VASP、Gaussian、CP2K、OpenMX 的计算正确地交给集群调度器：Slurm（`sbatch`）、PBS（`qsub`）、LSF（`bsub`）三者的命令与资源写法对照，`--nodes`/`--ntasks-per-node`/`--cpus-per-task`/`--mem`/`--time` 的含义，`mpirun` 与 `srun` 的区别，`OMP_NUM_THREADS` 与 `ulimit -s unlimited` 为什么必须写。随后分别给出四款软件的并行设置（VASP 的 `vasp_std/vasp_gam/vasp_ncl` 与 NCORE/KPAR，Gaussian 的 `%NProcShared`/`%Mem`/`GAUSS_SCRDIR`，CP2K 的 `cp2k.psmp` MPI×OpenMP 混合并行，OpenMX 的 `-nt`），以及作业数组参数扫描、`--dependency=afterok` 依赖链、本地 scratch、断点续算、`squeue`/`sacct`/`scancel` 监控和 dpdispatcher 批量提交（Slurm job script, PBS, LSF, MPI/OpenMP hybrid, job array, job dependency, restart）。

> [!NOTE]
> 登录集群、上传文件、`module` 系统与 CRLF 换行问题见 [第1章](01-linux-hpc.md)。本章的模板与具体集群无关，分区名、账户、本地盘路径、`module load` 的内容都必须按你所在集群的用户手册修改。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 读懂资源请求 | 说清 `--nodes`、`--ntasks-per-node`、`--cpus-per-task`、`--mem`、`--time`、`--partition` 各管什么，算出作业实际占用的核数 |
| 跨调度器迁移 | 把一个 Slurm 脚本改写成 PBS Pro/Torque 或 LSF 脚本（反之亦然） |
| 正确启动并行程序 | 选择 `srun` 或 `mpirun`；设置 `OMP_NUM_THREADS`、`ulimit -s unlimited`、线程绑定 |
| 四款软件的并行设置 | VASP 的 NCORE/KPAR 与核数匹配；Gaussian 的核数、内存与 scratch；CP2K 与 OpenMX 的 MPI×OpenMP |
| 批量与流程 | 用作业数组做参数扫描，用依赖链串联弛豫→静态→能带，用 dpdispatcher 批量提交 |
| 稳健运行 | 使用每个作业独立的 scratch 目录，在时限前让程序体面退出并续算 |
| 监控与排错 | 用 `squeue`、`sacct`、`scontrol`、`scancel` 查看和管理作业，读懂常见的调度器报错 |

配套示例：

- [examples/02-job-templates](../examples/02-job-templates/README.md)：通用模板（纯 MPI、MPI×OpenMP、共享内存、数组、PBS、LSF）、VASP 专用模板、自动续算模板、`sacct` 汇总脚本；
- [examples/02-array-chain-si](../examples/02-array-chain-si/README.md)：Si 的 ENCUT 作业数组、弛豫→静态→能带依赖链和 dpdispatcher 示例。

各软件更完整的专用脚本：[VASP/scripts/job_vasp_slurm.sh](../../VASP/scripts/job_vasp_slurm.sh)、[Gaussian/scripts/sub_g16.sh](../../Gaussian/scripts/sub_g16.sh)、[CP2K/scripts/job_cp2k_slurm.sh](../../CP2K/scripts/job_cp2k_slurm.sh)、[OpenMX/scripts/job_omx_slurm.sh](../../OpenMX/scripts/job_omx_slurm.sh)。本板块的 `mkjob.py` 可以按“调度器 + 软件 + 核数”生成作业脚本，用 `python ../scripts/mkjob.py -h` 查看用法（脚本位于 [Common/scripts/mkjob.py](../scripts/mkjob.py)）。

---

## 2. 原理速览：调度器、资源与并行层次

### 2.1 调度器做什么

集群上的计算节点由调度器统一分配。你提交的是一个**批处理脚本**：脚本开头的特殊注释行（`#SBATCH`、`#PBS`、`#BSUB`）描述资源需求，后面是普通的 bash 命令。调度器按优先级（fair-share、排队时间、作业规模）和空闲资源决定何时在哪些节点上启动它。

```text
 sbatch job.sh ──> 排队 (PD, Pending) ──> 运行 (R) ──> 结束 (CD 完成 / F 失败 / TO 超时 / OOM / CA 取消)
                     │                       │
                     └ REASON: Priority,      └ 脚本在第一个分配到的节点上执行；
                       Resources, Dependency     srun/mpirun 再把并行进程铺到所有分配的节点上
```

### 2.2 三个层次：节点、任务、线程

Slurm 的术语最清楚，PBS 和 LSF 只是换了名字：

- **节点（node）**：一台物理服务器，有若干 CPU 插槽（socket），每个插槽又分成一个或多个 NUMA 域；
- **任务（task）**：一个进程，通常就是一个 MPI rank；
- **每任务 CPU 数（cpus-per-task）**：一个进程可以使用的核数，用于 OpenMP 线程。

作业实际占用的核数：

$$
N_\text{cores} = N_\text{nodes}\times N_\text{tasks-per-node}\times N_\text{cpus-per-task}
$$

| 并行方式 | 典型程序 | `--ntasks-per-node` | `--cpus-per-task` | `OMP_NUM_THREADS` |
|---|---|---|---|---|
| 纯 MPI | VASP（CPU 版）、OpenMX、CP2K（`popt` 方式） | 每节点物理核数，如 64 | 1 | 1 |
| MPI×OpenMP 混合 | `cp2k.psmp`、`openmx -nt N`、VASP GPU 版 | 如 16 | 如 4 | = cpus-per-task |
| 共享内存（单进程多线程） | Gaussian 16、Multiwfn | 1（`--ntasks=1`） | 如 32 | Gaussian 用 `%NProcShared`/`%CPU` 指定 |

### 2.3 并行效率与“核数不是越多越好”

设程序中可并行部分的比例为 $p$，Amdahl 定律给出 $N$ 核上的理论加速比：

$$
S(N)=\frac{1}{(1-p)+p/N}
$$

并行效率 $E(N)=T_1/(N\,T_N)$。实际中通信开销会让 $E$ 随 $N$ 下降得更快：两原子的 Si 用 256 核往往比 16 核还慢。经验规则是先做一次**扩展性测试**（16/32/64/128 核各跑一次相同的短计算，比较每个电子步或每个离子步的耗时），选效率仍在 70% 以上的最大核数（见“参数选择与收敛建议”）。

### 2.4 内存是按节点给的

一个节点上的所有进程共享该节点的内存。纯 MPI 作业每个进程都持有一份数据副本，进程越多总内存越大；混合并行和 VASP 的 KPAR 也会改变内存需求。Slurm 中 `--mem` 表示**每节点**内存，`--mem-per-cpu` 表示每个分配到的 CPU 的内存，两者只能选一个。被系统因内存不足杀掉的作业状态为 `OUT_OF_MEMORY`（OOM）。

---

## 3. Slurm、PBS、LSF 对照表

### 3.1 命令

| 操作 | Slurm | PBS Pro / OpenPBS（Torque 差异） | LSF |
|---|---|---|---|
| 提交 | `sbatch job.sh` | `qsub job.sh` | `bsub < job.sh`（注意 `<`） |
| 交互式作业 | `salloc -N1 -n4 -t 1:00:00` 或 `srun -n4 --pty bash` | `qsub -I -l select=1:ncpus=4` | `bsub -Is -n 4 bash` |
| 查看自己的作业 | `squeue -u $USER` | `qstat -u $USER` | `bjobs` |
| 作业详情 | `scontrol show job <id>` | `qstat -f <id>` | `bjobs -l <id>` |
| 已结束作业 | `sacct -j <id>` | `qstat -x <id>`（PBS Pro）；Torque 用 `tracejob` | `bhist -l <id>`、`bjobs -a` |
| 取消 | `scancel <id>` | `qdel <id>` | `bkill <id>` |
| 暂缓/释放排队中的作业 | `scontrol hold/release <id>` | `qhold`/`qrls <id>` | `bstop`/`bresume <id>` |
| 队列/分区 | `sinfo -s` | `qstat -Q` | `bqueues` |
| 节点状态 | `sinfo -N -l` | `pbsnodes -a` | `bhosts` |
| 查看运行中作业的输出 | 直接 `tail -f` 输出文件 | 视配置，常在作业结束后才拷回 | `bpeek <id>` |

### 3.2 资源写法

| 含义 | Slurm（`#SBATCH`） | PBS Pro（`#PBS`） | Torque（`#PBS`） | LSF（`#BSUB`） |
|---|---|---|---|---|
| 作业名 | `--job-name=si` | `-N si` | `-N si` | `-J si` |
| 队列/分区 | `--partition=cpu` | `-q workq` | `-q batch` | `-q normal` |
| 节点与进程 | `--nodes=2 --ntasks-per-node=64` | `-l select=2:ncpus=64:mpiprocs=64` | `-l nodes=2:ppn=64` | `-n 128 -R "span[ptile=64]"` |
| 每进程线程 | `--cpus-per-task=4` | `select=...:ompthreads=4` | — | `-R "affinity[core(4)]"`（视站点配置） |
| 墙钟时限 | `--time=1-00:00:00` | `-l walltime=24:00:00` | `-l walltime=24:00:00` | `-W 24:00`（时:分） |
| 内存 | `--mem=180G`（每节点） | `select=...:mem=180gb` | `-l mem=180gb` | `-M`、`-R "rusage[mem=...]"`（单位由站点决定） |
| 标准输出 | `--output=%x-%j.out` | `-o file`，`-j oe` 合并 | 同左 | `-o %J.out` |
| 计费账户 | `--account=proj` | `-A proj` | `-A proj` | `-P proj` |
| 作业数组 | `--array=0-6%4` | `-J 0-6` | `-t 0-6%4` | `-J "name[1-7]%4"` |
| 依赖 | `--dependency=afterok:<id>` | `-W depend=afterok:<id>` | 同左 | `-w "done(<id>)"` |

### 3.3 作业内可用的环境变量

| 含义 | Slurm | PBS | LSF |
|---|---|---|---|
| 提交目录 | `SLURM_SUBMIT_DIR`（作业默认就在此目录启动） | `PBS_O_WORKDIR`（作业在 `$HOME` 启动，**必须** `cd`） | `LS_SUBCWD` |
| 作业号 | `SLURM_JOB_ID` | `PBS_JOBID` | `LSB_JOBID` |
| 节点列表 | `SLURM_JOB_NODELIST` | `PBS_NODEFILE`（文件，每个槽一行） | `LSB_HOSTS`、`LSB_DJOB_HOSTFILE` |
| 总进程数 | `SLURM_NTASKS` | `wc -l < $PBS_NODEFILE` | `LSB_DJOB_NUMPROC` |
| 每任务 CPU | `SLURM_CPUS_PER_TASK`（仅在指定 `--cpus-per-task` 时才有） | `OMP_NUM_THREADS`（由 `ompthreads` 设置） | — |
| 数组下标 | `SLURM_ARRAY_TASK_ID` | `PBS_ARRAY_INDEX`（Torque: `PBS_ARRAYID`） | `LSB_JOBINDEX` |

完整的三种调度器模板见 [slurm_mpi.sh](../examples/02-job-templates/slurm_mpi.sh)、[pbs_mpi.sh](../examples/02-job-templates/pbs_mpi.sh)、[lsf_mpi.sh](../examples/02-job-templates/lsf_mpi.sh)。

---

## 4. 第一个 Slurm 作业：逐行读懂模板

### 4.1 最小可用的纯 MPI 脚本

下面是一个精简后的纯 MPI 作业，适用于 VASP、OpenMX 或纯 MPI 方式运行的 CP2K。带 CRLF 检查和信号处理的完整版本见 [slurm_mpi.sh](../examples/02-job-templates/slurm_mpi.sh)。

```bash
#!/bin/bash
#SBATCH --job-name=si-scf
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

module purge
module load intel-oneapi-mpi intel-oneapi-mkl vasp/6.5.1   # site specific

ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}" || exit 1
echo "Job ${SLURM_JOB_ID}: ${SLURM_NTASKS} ranks on ${SLURM_JOB_NODELIST}"

srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 vasp_std > vasp.out 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `#!/bin/bash` | 第一行 | 必须是第一行；`#SBATCH` 行必须位于第一条可执行命令之前，之后的 `#SBATCH` 会被当成普通注释 |
| `--job-name` | 有意义的短名 | `squeue` 中显示；`%x` 在输出文件名中代表它 |
| `--partition` | 按集群手册 | 分区（队列）；`sinfo -s` 查看；不同分区的节点核数、内存、时限可能不同 |
| `--nodes` | 1 起步 | 节点数；先确认单节点的扩展性再加节点 |
| `--ntasks-per-node` | 每节点物理核数 | 纯 MPI：每个核一个 rank；`lscpu` 中 `Core(s) per socket × Socket(s)` |
| `--cpus-per-task` | 纯 MPI 为 1 | 混合并行时为每个 rank 的线程数 |
| `--time` | 预计时间 × 1.5 | 格式 `HH:MM:SS` 或 `D-HH:MM:SS`；时间短的作业更容易被 backfill 提前调度 |
| `--mem` | `0` 或具体值 | 每节点内存；`--mem=0` 在多数站点表示“该节点全部内存”，但有的站点禁止，以手册为准 |
| `--output`/`--error` | `%x-%j.out` | `%j` 作业号、`%x` 作业名、`%A`/`%a` 数组主作业号/下标 |
| `--account`、`--qos` | 按集群要求 | 许多集群必须指定计费账户，否则提交报错 |
| `--exclusive` | 大作业建议 | 独占节点，避免与他人的作业抢内存带宽 |

提交与检查：

```bash
sbatch --test-only job.sh      # 只检查资源请求是否合法并估计开始时间，不真正提交
sbatch job.sh                  # 输出: Submitted batch job 123456
squeue -u $USER                # ST 列: PD 排队, R 运行, CG 正在结束
squeue -u $USER --start        # 调度器估计的开始时间（仅供参考）
tail -f si-scf-123456.out
```

命令行参数优先于脚本中的 `#SBATCH`：`sbatch --nodes=2 --time=4:00:00 job.sh` 可以临时修改而不必改脚本。

### 4.2 `mpirun` 与 `srun` 的区别

| | `srun` | `mpirun` / `mpiexec` |
|---|---|---|
| 属于 | Slurm | MPI 库（Intel MPI、Open MPI、MPICH） |
| 进程数与节点 | 自动继承作业分配，可用 `--ntasks`、`--cpus-per-task` 显式给出 | 大多数 MPI 库也能从 Slurm/PBS/LSF 读到分配；PBS 下有时需 `-machinefile $PBS_NODEFILE` |
| 进程管理 | 由 slurmstepd 启动、计费和清理，`sacct` 能看到每个 step 的内存与耗时 | 由 MPI 自己的启动器（hydra、PRRTE）启动 |
| 前提 | MPI 库需带 Slurm 的 PMI-2 或 PMIx 支持；`srun --mpi=list` 查看可用类型 | 只要 `module load` 的 MPI 与编译程序时用的一致即可 |
| 绑定 | `--cpu-bind=cores`，`--hint=nomultithread` 不用超线程 | Open MPI：`--bind-to core --map-by ...`；Intel MPI：`I_MPI_PIN_DOMAIN` |

建议：集群手册推荐什么就用什么。如果 `srun` 启动后输出被重复打印 N 次、每个进程都显示自己是 rank 0（或 VASP 报告 `running on 1 total cores` 却出现 64 份输出），说明 MPI 与 Slurm 的 PMI 不匹配，改用 `mpirun` 或询问管理员应使用的 `--mpi=pmi2`/`--mpi=pmix`。Intel MPI 不自带 Slurm 的 PMI 实现，Slurm 官方 MPI 指南给出的做法是在 `srun` 前 `export I_MPI_PMI_LIBRARY=/path/to/slurm/lib/libpmi2.so`（路径由站点决定，问管理员或 `find / -name "libpmi2.so*" 2>/dev/null`）。

> [!TIP]
> 在 `srun` 中显式写 `--cpus-per-task="${SLURM_CPUS_PER_TASK}"` 对所有 Slurm 版本都安全。Slurm 22.05–23.02 的 `srun` 不再从 `sbatch`/`salloc` 继承 `--cpus-per-task`（官方当时给出的替代办法是 `export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}`），混合并行时一个 rank 的全部线程会被挤在同一个核上；23.11 起取消了 `SRUN_CPUS_PER_TASK`，`srun` 重新读取 `SLURM_CPUS_PER_TASK`，恢复为继承。用 `sinfo --version` 查看集群的 Slurm 版本。

### 4.3 每个作业都要写的三行环境设置

```bash
ulimit -s unlimited          # stack size: avoids SIGSEGV in Fortran codes
export OMP_NUM_THREADS=1     # pure MPI: one thread per rank
export OMP_STACKSIZE=512m    # only for OpenMP runs: private stack per thread
```

- **`ulimit -s unlimited`**：VASP、CP2K、OpenMX 的 Fortran/C 代码在栈上放大数组，默认 8 MB 的栈会导致 `forrtl: severe (174): SIGSEGV, segmentation fault occurred`。有的站点不允许设为 unlimited，可以改为较大的数值（`ulimit -s 2000000`，单位 KB）。
- **`OMP_NUM_THREADS=1`**：纯 MPI 作业中，MKL、OpenBLAS 等库在不设该变量时可能在每个 rank 内再开满全部核的线程，64 个 rank × 64 个线程会让节点严重超载，计算反而慢十几倍。混合并行时设为 `${SLURM_CPUS_PER_TASK}`。
- **`OMP_STACKSIZE`**：只影响 OpenMP 线程的私有栈（`ulimit -s` 管不到它），CP2K/OpenMX 混合并行或 Multiwfn 处理大体系时建议设置。

另外两条习惯：`module purge` 后再 `module load` 与编译时完全一致的编译器/MPI/数学库（见 [第1章第9节](01-linux-hpc.md)）；`cd "${SLURM_SUBMIT_DIR}" || exit 1`，切换失败就立即退出，避免在错误的目录里覆盖文件。

### 4.4 MPI×OpenMP 混合并行与线程绑定

```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16     # MPI ranks per node
#SBATCH --cpus-per-task=4        # threads per rank: 16 x 4 = 64 cores

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
srun --ntasks="${SLURM_NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores cp2k.psmp -i in.inp > out.log 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `--cpus-per-task` | 2、4、8 | 让一个 rank 的全部线程落在同一个 NUMA 域内（`lscpu | grep NUMA` 查看每个域的核数） |
| `OMP_PLACES` | `cores` | 每个线程占一个物理核 |
| `OMP_PROC_BIND` | `close` | 同一 rank 的线程放在相邻的核上 |
| `--cpu-bind` | `cores` | Slurm 把每个 rank 绑定到它的那组核 |
| `--hint` | `nomultithread` | 节点开启超线程时，避免两个线程挤在同一个物理核的两个逻辑核上 |

完整模板（含 Open MPI 与 Intel MPI 的等价写法）见 [slurm_hybrid_mpi_omp.sh](../examples/02-job-templates/slurm_hybrid_mpi_omp.sh)。

---

## 5. 四款软件的并行设置

### 5.1 VASP：选对可执行文件，让 NCORE/KPAR 与核数匹配

VASP 编译后有三个可执行文件：

| 可执行文件 | 用途 | 说明 |
|---|---|---|
| `vasp_std` | 一般计算 | 默认选择 |
| `vasp_gam` | 只用 Γ 点（大超胞、分子、AIMD） | 利用 Γ 点波函数为实数的性质，内存和时间明显减少；KPOINTS 必须只含 Γ 点 |
| `vasp_ncl` | 非共线磁性、自旋轨道耦合（`LNONCOLLINEAR`、`LSOCORB`） | 计算量约为 `vasp_std` 的数倍 |

VASP 的 MPI 进程按层次分组：

$$
N_\text{ranks} = \text{IMAGES}\times\text{KPAR}\times N_\text{band groups}\times \text{NCORE}
$$

（不做 NEB 时 IMAGES = 1。）

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `NCORE` | 默认 1；快速网络的集群取 2 到每节点核数之间，VASP Wiki 的经验值：约 100 原子取 4，400 原子以上 12–16 | 每个能带由 NCORE 个进程合作做 FFT；应是每节点核数的因子，使 FFT 组的通信留在一个节点内（最好一个 NUMA 域内）；千兆以太网或 ≤ 8 核的小作业用 1 |
| `KPAR` | 能整除总进程数，并尽量整除不可约 k 点数 | k 点并行：每组 = 总进程数 / KPAR 个进程负责一个 k 点；内存随 KPAR 增加；Γ 点计算只能是 1 |
| `NPAR` | 不设 | 旧标签，与 NCORE 二选一；同时出现时 NPAR 优先 |
| `NSIM` | 默认 4（CPU） | 一次同时优化的能带数；GPU 版要设大很多 |

> [!IMPORTANT]
> 当 `OMP_NUM_THREADS > 1` 或使用 GPU 版本时，VASP 会自动把 NCORE 重置为 1，FFT 的并行改由线程/GPU 完成。VASP 的 OpenMP 与 GPU 运行方式和并行参数的详细讨论见 [VASP 第1章第13节](../../VASP/tutorials/01-setup-compile.md)。

一个 2 节点 × 64 核、每节点 8 个 NUMA 域、100 原子表面、6 个不可约 k 点的例子：

```text
# INCAR (parallel part)
NCORE = 8      ! one band per NUMA domain (8 cores)
KPAR  = 2      ! two k-point groups, one per node: 128 / 2 = 64 ranks per group
```

[slurm_vasp.sh](../examples/02-job-templates/slurm_vasp.sh) 在运行前自动检查 NCORE 是否整除每节点进程数、KPAR 是否整除总进程数，选择可执行文件用环境变量：

```bash
sbatch slurm_vasp.sh                              # vasp_std
sbatch --export=ALL,VARIANT=gam slurm_vasp.sh     # vasp_gam
sbatch --export=ALL,VARIANT=ncl slurm_vasp.sh     # vasp_ncl
```

运行后检查 VASP 实际采用的并行划分：

```bash
grep -E "running|distrk|distr:" OUTCAR | head
```

```text
 running   128 mpi-ranks, with    1 threads/rank, on    2 nodes
 distrk:  each k-point on   64 cores,    2 groups
 distr:  one band on NCORE=   8 cores,    8 groups
```

（示意输出，数值仅供参考；不同 6.x 版本的措辞略有差异。）

成功标准：`running` 行的 rank 数等于 `--nodes × --ntasks-per-node`、threads/rank 为 1（纯 MPI）；`distrk` 的组数等于 KPAR，`distr` 中的 NCORE 等于 INCAR 中的设置。NEB 计算还要求总进程数是 `IMAGES` 的整数倍，否则报 `M_divide: can not subdivide N nodes by M`。

### 5.2 Gaussian 16：单节点共享内存，核数、内存、scratch 三者一致

Gaussian 16 的普通版本只在一个节点内用线程并行（跨节点需要 Linda 版本）。作业请求与输入文件必须对应：

```text
%NProcShared=16
%Mem=32GB
%Chk=benzene.chk
# B3LYP/6-31G(d) Opt Freq

benzene opt+freq

0 1
...
```

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G

export GAUSS_SCRDIR="${TMPDIR:-/tmp}/${USER}/g16_${SLURM_JOB_ID}"
mkdir -p "${GAUSS_SCRDIR}"
trap 'rm -rf "${GAUSS_SCRDIR}"' EXIT
g16 < benzene.gjf > benzene.log 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `--nodes` / `--ntasks` | 1 / 1 | 共享内存程序只能用一个节点、一个进程 |
| `--cpus-per-task` | 8–32 | 与 `%NProcShared` 相等。Gaussian 官方推荐用 `%CPU` 列出核号并绑定线程，但核号必须是调度器**实际分配给本作业**的核：与他人共享节点时 Slurm 可能给你 16–31 号核，写死 `%CPU=0-15` 会绑到别人的核上。共享节点用 `%NProcShared`，或从 `grep Cpus_allowed_list /proc/self/status` 读取分配的核号再生成 `%CPU`；独占节点（`--exclusive`）时 `%CPU=0-15` 才安全 |
| `%Mem` | `--mem` 的 75–85% | `%Mem` 只是 Gaussian 自己的工作内存，程序本身和 I/O 缓冲还要额外占用；`%Mem` 等于 `--mem` 时容易被 OOM 杀掉 |
| `GAUSS_SCRDIR` | 计算节点**本地盘**上每个作业独立的目录 | 存放 `.rwf`、`.int`、`.d2e` 等大型临时文件；放在共享文件系统上会拖慢计算并影响他人 |
| `%Chk` | 放在提交目录 | `.chk` 需要保留用于续算和后处理，不能放在会被删除的 scratch 中 |

> [!WARNING]
> 不要在 Gaussian 作业中设置 `OMP_NUM_THREADS=1` 来“限制线程”，Gaussian 的核数只由 `%NProcShared`/`%CPU`（或对应的 `GAUSS_PDEF`/`GAUSS_CDEF` 环境变量）控制。多个 Gaussian 作业共享节点时，各自的核数加起来不能超过分配给你的核数。

完整模板见 [slurm_shared_memory.sh](../examples/02-job-templates/slurm_shared_memory.sh)（自动比对 `%NProcShared` 与 `--cpus-per-task`）；更多 Gaussian 专用内容（`%CPU`、`Default.Route`、Linda、GPU、作业数组、`Opt=Restart`）见 [Gaussian 第17章](../../Gaussian/tutorials/17-hpc-jobs.md) 和 [sub_g16.sh](../../Gaussian/scripts/sub_g16.sh)。

### 5.3 CP2K：`cp2k.psmp` 的 MPI×OpenMP

CP2K 2026.x 的生产计算用 `cp2k.psmp`（MPI + OpenMP）。`cp2k.popt` 在新版本中只是指向 psmp 的链接，以这个名字启动时线程数固定为 1，相当于纯 MPI（见 [CP2K 第1章](../../CP2K/tutorials/01-install.md)）。

```bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores OMP_PROC_BIND=close OMP_STACKSIZE=512m
srun --ntasks="${SLURM_NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
     cp2k.psmp -i h2o_md.inp -o h2o_md.out
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| ranks × threads | 每节点 = 物理核数 | 从“每 rank 1–2 线程”开始测试，再试 4、8 |
| `OMP_NUM_THREADS` | = `--cpus-per-task` | 不设置时 OpenMP 可能开满全部核，造成超载 |
| `-i` / `-o` | 输入 / 输出文件 | 也可以 `cp2k.psmp in.inp > out.log` |
| `&GLOBAL/WALLTIME` | 比作业时限短 10–20 分钟 | 到时 CP2K 主动结束并写出重启文件（见第 9 节） |

运行后确认实际的进程与线程数：

```bash
grep -E "Total number of message passing processes|Number of threads for this process" h2o_md.out
```

两个数应分别等于 `SLURM_NTASKS` 和 `--cpus-per-task`。大体系、杂化泛函（ADMM）和 k 点计算对内存敏感，增加每 rank 线程数、减少 rank 数可以降低内存。CP2K 专用脚本见 [job_cp2k_slurm.sh](../../CP2K/scripts/job_cp2k_slurm.sh)。

### 5.4 OpenMX：`mpirun openmx in.dat -nt N`

OpenMX 的 OpenMP 线程数用命令行参数 `-nt` 指定；手册说明不写 `-nt` 时线程数为 1，即普通的纯 MPI 并行。建议让 `OMP_NUM_THREADS` 与 `-nt` 相同，避免数学库自行开线程：

```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=512m
mpirun -np "${SLURM_NTASKS}" openmx Si.dat -nt "${OMP_NUM_THREADS}" > Si.std 2>&1
# or: srun --ntasks=... --cpus-per-task=... openmx Si.dat -nt ${OMP_NUM_THREADS}
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `-np`（MPI 进程数） | 周期体系可取到 自旋数 × 不可约 k 点数 或更多 | 4.0 的能带（周期）计算依次在自旋、k 点、本征态三层上做 MPI 并行；团簇计算在自旋和本征态上并行；O(N) 方法按原子分配。4.0 取消了 3.x 中“进程数不能超过原子数”的限制，并改为三维区域分解，不必再为并行效率挑选晶胞取向 |
| `-nt`（线程数） | 1、2、4 | 进程数 × 线程数 ≤ 物理核数；手册的初步测试表明混合并行主要节省内存，速度与纯 MPI 相近 |
| `DATA.PATH` | 指向 DFT_DATA19 | 输入中的赝势与基组目录，计算节点上必须可访问 |

> [!NOTE]
> 进程数大于“自旋数 × k 点数”时，OpenMX 4.0 每个 k 点只需对角化一次（否则为节省内存要对角化两次），因此跨过这个阈值时加速比可能超过理想值；这是算法切换造成的，不代表并行效率真的超过 100%。见 OpenMX 4.0 手册第 28.1 节。

先用 `-runtest` 确认安装正确（见 [OpenMX 第2章](../../OpenMX/tutorials/02-install.md)），再用自己的输入、把 `scf.maxIter` 设小，在目标节点上比较不同 `-np`/`-nt` 组合的耗时。OpenMX 专用脚本见 [job_omx_slurm.sh](../../OpenMX/scripts/job_omx_slurm.sh)。

---

## 6. PBS 与 LSF 模板

### 6.1 PBS Pro / OpenPBS

```bash
#!/bin/bash
#PBS -N si-scf
#PBS -q workq
#PBS -l select=2:ncpus=64:mpiprocs=64:ompthreads=1:mem=200gb
#PBS -l walltime=24:00:00
#PBS -j oe

module purge
module load vasp/6.5.1          # site specific
ulimit -s unlimited
export OMP_NUM_THREADS=1
cd "${PBS_O_WORKDIR}" || exit 1          # PBS starts in $HOME!
NP=$(wc -l < "${PBS_NODEFILE}")
mpirun -np "${NP}" -machinefile "${PBS_NODEFILE}" vasp_std > vasp.out 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `select=N` | 节点数 | 申请 N 个“chunk”，通常一个 chunk 对应一个节点 |
| `ncpus` | 每节点核数 | 每个 chunk 的核数 |
| `mpiprocs` | 纯 MPI 时等于 `ncpus` | 决定 `$PBS_NODEFILE` 中每个节点出现几行，即 MPI 进程数 |
| `ompthreads` | 每进程线程数 | PBS Pro 会据此设置 `OMP_NUM_THREADS` |
| `-l place=scatter:excl` | 大作业可选 | 每个 chunk 放在不同节点并独占 |
| `-j oe` | 推荐 | 合并标准输出与标准错误 |

Torque 的写法是 `#PBS -l nodes=2:ppn=64`。完整模板见 [pbs_mpi.sh](../examples/02-job-templates/pbs_mpi.sh)（默认 `EXE=vasp_std`，即 PBS 版的 VASP 脚本）。

### 6.2 LSF

```bash
#!/bin/bash
#BSUB -J si-scf
#BSUB -q normal
#BSUB -n 128
#BSUB -R "span[ptile=64]"
#BSUB -W 24:00
#BSUB -o %J.out
#BSUB -e %J.err

module purge
module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
cd "${LS_SUBCWD:-$PWD}" || exit 1
mpirun -np "${LSB_DJOB_NUMPROC}" vasp_std > vasp.out 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `-n` | 总槽数 | LSF 按总数申请，不是每节点数 |
| `span[ptile=64]` | 每节点核数 | 每个主机放 64 个槽；共享内存程序用 `span[hosts=1]` |
| `-W` | `HH:MM` | 墙钟时限，注意没有秒 |
| `bsub < job.sh` | 必须用 `<` | `bsub job.sh` 会把脚本当作命令执行，`#BSUB` 行不会被解析 |

完整模板见 [lsf_mpi.sh](../examples/02-job-templates/lsf_mpi.sh)。

---

## 7. 作业数组：参数扫描

### 7.1 为什么用数组

ENCUT、k 点、晶格常数、CP2K 的 CUTOFF、Gaussian 的一组构象……这类“同一脚本、参数不同”的计算适合做成作业数组：一次 `sbatch` 提交 N 个子任务，每个子任务通过 `SLURM_ARRAY_TASK_ID` 取自己的参数；`%4` 限制同时运行的数量，不会一下子占满配额。

### 7.2 Si 的 ENCUT 扫描

[slurm_array_encut.sh](../examples/02-array-chain-si/slurm_array_encut.sh) 的核心部分：

```bash
#SBATCH --job-name=si-encut
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --array=0-6%4               # 7 values, at most 4 running at once
#SBATCH --output=%x-%A_%a.out

VALUES=(250 300 350 400 450 500 550)
E=${VALUES[$SLURM_ARRAY_TASK_ID]}
D=encut_${E}
mkdir -p "${D}"
cp POSCAR KPOINTS POTCAR "${D}/"
sed "s/@ENCUT@/${E}/g" INCAR.encut.template > "${D}/INCAR"
cd "${D}" && srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 vasp_std > vasp.out 2>&1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `--array=0-6%4` | 下标范围 + 并发上限 | 也可写 `1,3,5`、`0-20:2`（步长 2） |
| `%A` / `%a` | 输出文件名中 | 数组主作业号 / 子任务下标，避免各子任务覆盖同一个输出文件 |
| `SLURM_ARRAY_TASK_ID` | 作业内变量 | 子任务下标，用作数组 `VALUES` 的索引 |
| `@ENCUT@` | 模板占位符 | `sed` 替换后写入每个子目录，其余参数保持完全一致 |

```bash
sbatch slurm_array_encut.sh
squeue -u $USER -r                       # -r: one line per array task
scancel 123456_3                         # cancel only task 3
scontrol update JobId=123456 ArrayTaskThrottle=2   # change the %N limit later
sbatch --array=2,5 slurm_array_encut.sh  # rerun failed values only
bash collect_encut.sh                    # summary: E/atom and dE vs. 550 eV
```

同一扫描的 PBS 与 LSF 写法见 [pbs_array_encut.sh](../examples/02-array-chain-si/pbs_array_encut.sh)（`#PBS -J 0-6`、`$PBS_ARRAY_INDEX`）和 [lsf_array_encut.sh](../examples/02-array-chain-si/lsf_array_encut.sh)（`#BSUB -J "si-encut[1-7]%4"`、`$LSB_JOBINDEX`，下标从 1 开始）。通用的“模板目录 + 占位符”数组脚本见 [slurm_array.sh](../examples/02-job-templates/slurm_array.sh)，收敛判据与作图见 [VASP 第3章](../../VASP/tutorials/03-convergence.md)。

> [!TIP]
> 很多集群限制数组的最大长度（`scontrol show config | grep MaxArraySize`）和每个用户的排队作业数。几百个只需几分钟的小计算，更适合在**一个**作业里用循环依次运行，或交给 dpdispatcher 的 `group_size` 打包（第 11 节）。

---

## 8. 依赖链：弛豫 → 静态 → 能带

能带计算依赖静态计算的 CHGCAR，静态计算依赖弛豫后的结构。用 `--dependency` 一次性提交三步，前一步成功结束（退出码为 0）后一步才开始：

```bash
j1=$(sbatch --parsable --export=ALL,STEP=relax  job_step.sh); j1=${j1%%;*}   # drop ";cluster"
j2=$(sbatch --parsable --export=ALL,STEP=static --dependency=afterok:${j1} --kill-on-invalid-dep=yes job_step.sh); j2=${j2%%;*}
j3=$(sbatch --parsable --export=ALL,STEP=band   --dependency=afterok:${j2} --kill-on-invalid-dep=yes job_step.sh); j3=${j3%%;*}
echo "relax=${j1} static=${j2} band=${j3}"
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `--parsable` | 必用 | `sbatch` 只输出作业号（多集群时为 `jobid;cluster`），便于存入变量 |
| `afterok:<id>` | 串联计算 | 前一作业**退出码为 0** 才开始；其他类型：`afterany`（无论成败）、`afternotok`（失败后运行，可用于清理或报警）、`aftercorr`（数组之间按下标一一对应）、`singleton`（同名作业依次运行） |
| `--kill-on-invalid-dep=yes` | 推荐 | 依赖不可能满足时自动取消；否则作业一直挂在队列里，`REASON` 显示 `DependencyNeverSatisfied` |
| `--export=ALL,STEP=...` | 传参 | 把变量 `STEP` 传进作业，同一脚本完成三个步骤 |

关键在于**作业脚本的退出码必须反映计算是否真的成功**。VASP 未收敛时通常也以 0 退出，所以 [job_step.sh](../examples/02-array-chain-si/job_step.sh) 在每一步后主动检查：

```bash
finished_ok() { grep -q "General timing and accounting" OUTCAR 2>/dev/null; }
finished_ok || exit 2
grep -q "reached required accuracy" OUTCAR || { echo "relaxation NOT converged" >&2; exit 3; }
```

一键提交脚本见 [chain_relax_static_band.sh](../examples/02-array-chain-si/chain_relax_static_band.sh)，在登录节点运行 `bash chain_relax_static_band.sh` 即可。PBS 的写法是 `qsub -W depend=afterok:${j1} job.sh`，LSF 是 `bsub -w "done(${j1})" < job.sh`。

检查：`squeue -u $USER` 中后两步的 `REASON` 为 `(Dependency)`；`sacct -j ${j1},${j2},${j3} --format=JobID,JobName,State,ExitCode` 在全部成功时显示三个 `COMPLETED  0:0`。

---

## 9. scratch 目录与断点续算

### 9.1 三类存储

| 存储 | 典型路径 | 特点 | 放什么 |
|---|---|---|---|
| 家目录 | `/home/$USER` | 小、有备份 | 脚本、源码、配置 |
| 共享工作盘（并行文件系统） | `/work`、`/data`、`/scratch/$USER` | 大、快、通常不备份，可能定期清理 | 输入和需要保留的输出 |
| 计算节点本地盘 | `$TMPDIR`、`/tmp`、`$SLURM_TMPDIR`、`/local` | 最快、作业结束后可能被清空 | Gaussian 的 `.rwf`、CP2K/VASP 的临时大文件 |

做法：在作业里建一个以作业号命名的目录，结束时（包括出错和被取消时）清理：

```bash
SCR="${TMPDIR:-/tmp}/${USER}/${SLURM_JOB_ID}"
mkdir -p "${SCR}" || exit 1
trap 'rm -rf "${SCR}"' EXIT
```

> [!CAUTION]
> 如果把整个计算放在本地盘上运行，必须在脚本结束前把结果拷回共享盘（`trap` 中先 `cp` 再 `rm`），而且时限被强制终止时 `trap` 不一定有机会执行完。对 VASP、CP2K、OpenMX，通常直接在共享工作盘上运行更安全；只有 Gaussian 这类产生大量临时 I/O 的程序才必须用本地 scratch。

### 9.2 在时限前体面退出

`#SBATCH --signal=B:USR1@900` 让 Slurm 在时限前 900 秒向批处理脚本发送 `USR1`，脚本中用 `trap` 接住它，再通知程序自行停止：

| 程序 | 体面停止的方法 | 续算方法 |
|---|---|---|
| VASP | 写 `STOPCAR`：`LSTOP = .TRUE.`（完成当前离子步后停止）或 `LABORT = .TRUE.`（完成当前电子步后停止） | `cp CONTCAR POSCAR`；保留 WAVECAR 可加快第一步自洽 |
| CP2K | `&GLOBAL` 中 `WALLTIME`；或在工作目录 `touch EXIT`（只结束当前的 MD 或几何优化循环可用 `EXIT_MD`、`EXIT_GEO`，也可加项目名前缀，如 `h2o_md.EXIT`） | `&EXT_RESTART` + `RESTART_FILE_NAME <PROJECT>-1.restart`，波函数用 `.wfn` |
| Gaussian | 不支持中途体面停止；靠 `.chk` 续算 | `%Chk` 放在提交目录；`Opt=Restart`、`IRC=Restart`、频率计算用 `Restart` 关键词（需要 `.rwf`） |
| OpenMX | 设置合理的 `MD.maxIter`、`scf.maxIter` 分段运行 | 电子结构 `scf.restart on`（读 `<System.Name>_rst/`）；结构可用输出的最终坐标重建输入 |

VASP 的完整实现见 [slurm_resubmit_vasp.sh](../examples/02-job-templates/slurm_resubmit_vasp.sh)：收到 `USR1` 后写 STOPCAR，VASP 完成当前离子步后正常退出；若尚未收敛，脚本把 CONTCAR 复制为 POSCAR、旧的 OUTCAR 改名保存，并重新提交自己（最多 `MAX_ROUND` 轮）。

```bash
#SBATCH --signal=B:USR1@900
stop_vasp() { echo "LSTOP = .TRUE." > STOPCAR; }
trap stop_vasp USR1
srun ... vasp_std > vasp.out 2>&1 &
PID=$!
while true; do wait "${PID}"; status=$?; kill -0 "${PID}" 2>/dev/null || break; done
```

注意 `srun ... &` 加 `wait` 的写法：bash 只有在等待后台进程时才能立即响应信号；如果程序在前台运行，`trap` 要等程序结束后才会执行，那就没有意义了。`B:` 表示信号只发给批处理脚本（不发给 VASP 进程本身）；sbatch 手册说明信号最多可能提前 60 秒发出。`LSTOP` 要等当前离子步结束才生效，所以提前量必须大于一个离子步的耗时（大体系的杂化泛函一步可能超过 15 分钟，相应加大 `@900`）。

CP2K 的对应写法（2026.2 与旧版语法相同）：

```text
&GLOBAL
  PROJECT  h2o_md
  RUN_TYPE MD
  ! job limit 24:00:00 minus 20 min; the manual requires double quotes
  ! around a WALLTIME value that is followed by an inline comment
  WALLTIME 23:40:00
&END GLOBAL
&EXT_RESTART
  RESTART_FILE_NAME h2o_md-1.restart
&END EXT_RESTART
```

`WALLTIME` 可以写秒数或 `HH:MM:SS`；CP2K 手册注明，值后面若紧跟行内注释，必须写成 `WALLTIME "23:40:00"  ! ...`，所以上面把注释单独放一行。第一次运行时不写 `&EXT_RESTART`，续算时再加上（细节见 [CP2K 第2章](../../CP2K/tutorials/02-input-structure.md) 和 [CP2K 第7章](../../CP2K/tutorials/07-geo-cell-opt.md)）。Gaussian 续算的完整示例见 [Gaussian 第17章](../../Gaussian/tutorials/17-hpc-jobs.md)，OpenMX 的 `scf.restart` 见 [OpenMX 第3章](../../OpenMX/tutorials/03-input-first-calc.md)。

---

## 10. 监控与管理作业

### 10.1 常用命令

```bash
squeue -u $USER                                   # my jobs
squeue -u $USER -o '%.10i %.12j %.2t %.10M %.6D %R' # id, name, state, time, nodes, reason/nodelist
squeue -u $USER --start                           # estimated start times of pending jobs
scontrol show job 123456                          # everything about one job (WorkDir, NumCPUs, Reason ...)
sacct -j 123456 --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,MaxRSS
sacct -u $USER -S 2026-11-01 -X                   # all main jobs since a date
sstat -a -j 123456 --format=JobID,MaxRSS,AveCPU   # running job, all steps (srun-launched steps)
seff 123456                                       # CPU / memory efficiency (contrib tool, if installed)
scancel 123456                                    # cancel one job
scancel -u $USER --state=PENDING                  # cancel all my pending jobs
scancel -n si-encut                               # cancel by job name
scontrol hold 123456 ; scontrol release 123456    # keep a pending job from starting / let it go
scontrol update JobId=123456 TimeLimit=12:00:00   # users can usually only DECREASE the limit
sinfo -s ; sshare -U ; sprio -u $USER             # partitions, fair-share, priority
```

[sacct_report.sh](../examples/02-job-templates/sacct_report.sh) 把最近几天的作业汇总成表，并单独列出 `FAILED`、`TIMEOUT`、`OUT_OF_MEMORY`、`NODE_FAIL` 的作业。要批量判断几十个计算目录是否“真正算完”（VASP 的 `reached required accuracy`、Gaussian 的 `Normal termination`、CP2K 的 SCF 是否收敛），用本板块的 [job_status.py](../scripts/job_status.py)，`python ../scripts/job_status.py -h` 查看用法。

### 10.2 读懂排队原因（`REASON` 列）

| REASON | 含义 | 处理 |
|---|---|---|
| `Priority` | 有优先级更高的作业 | 等待；缩短 `--time` 可能被 backfill 提前 |
| `Resources` | 你是下一个，等待资源释放 | 等待 |
| `Dependency` | 依赖的作业未结束 | 正常 |
| `DependencyNeverSatisfied` | 依赖的作业已失败 | `scancel`；下次加 `--kill-on-invalid-dep=yes` |
| `QOSMaxJobsPerUserLimit`、`AssocGrpCPULimit` 等 | 触及用户/组的配额 | 等前面的作业结束，或减少并发（数组 `%N`） |
| `ReqNodeNotAvail` | 请求的节点在维护或保留中 | 换分区或等待维护结束 |

### 10.3 作业结束后的检查清单

1. `sacct -j <id> --format=JobID,State,ExitCode,Elapsed,MaxRSS`：`State` 为 `COMPLETED`、`ExitCode` 为 `0:0`；
2. 程序自己的完成标志：VASP OUTCAR 末尾的 `General timing and accounting`、Gaussian 的 `Normal termination`、CP2K 输出末尾的 `PROGRAM ENDED AT`、OpenMX `.out` 中的计算耗时汇总；
3. 收敛标志：VASP 弛豫的 `reached required accuracy`，CP2K 中没有 `SCF run NOT converged`；
4. `MaxRSS` 与 `--mem` 对比，下次申请合理的内存；`seff` 显示的 CPU 效率明显低于 80% 时，检查是否进程/线程数设置错误。

---

## 11. dpdispatcher：用 Python 批量提交

[dpdispatcher](https://github.com/deepmodeling/dpdispatcher)（DeepModeling 社区，LGPL-3.0）根据 `Machine`（调度器类型、本地或 SSH 远程）和 `Resources`（节点、核数、队列、模块）自动生成 Slurm/PBS/LSF 脚本，提交、轮询状态，并把指定文件拷回。它适合几十到几千个结构的高通量计算或机器学习势的数据集标注，也被 DP-GEN 等工具使用。

[dpdispatcher_encut.py](../examples/02-array-chain-si/dpdispatcher_encut.py) 的核心：

```python
machine = Machine.load_from_dict({
    "batch_type": "Slurm",            # also: PBS, Torque, LSF, Shell
    "context_type": "LocalContext",   # SSHContext for a remote cluster
    "local_root": "./",
    "remote_root": "./dp_work",
    "remote_profile": {"symlink": False},
})
resources = Resources.load_from_dict({
    "number_node": 1, "cpu_per_node": 8, "gpu_per_node": 0,
    "queue_name": "cpu", "group_size": 1,
    "custom_flags": ["#SBATCH --time=00:30:00"],
    "module_list": [], "envs": {"OMP_NUM_THREADS": 1},
})
task = Task(command="mpirun -np 8 vasp_std", task_work_path="dp_encut_400/",
            forward_files=["INCAR", "POSCAR", "KPOINTS", "POTCAR"],
            backward_files=["OUTCAR", "OSZICAR"])
Submission(work_base="./", machine=machine, resources=resources,
           task_list=[task]).run_submission()
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `batch_type` | `Slurm` / `PBS` / `Torque` / `LSF` / `Shell` | 生成哪种调度器的脚本；`Shell` 直接在本机后台运行 |
| `context_type` | `LocalContext` / `SSHContext` | 在集群上运行用前者；从自己电脑经 SSH 提交用后者（需在 `remote_profile` 中填写主机、用户名、密钥） |
| `group_size` | 1 或更大 | 每个作业脚本打包几个任务；大量小任务时增大，可减少排队作业数 |
| `custom_flags` | 时限、账户等 | 直接写进生成脚本的额外 `#SBATCH` 行 |
| `forward_files` / `backward_files` | 输入 / 需要的输出 | 只拷回需要的文件，不要拷回 WAVECAR、CHGCAR |

```bash
python dpdispatcher_encut.py --dry-run    # stage directories only, nothing is submitted
python dpdispatcher_encut.py              # submit and wait (run it inside tmux)
```

`run_submission()` 会一直等待到所有任务结束，因此要在 tmux/screen 中运行（见 [第1章第6节](01-linux-hpc.md)）。中断后再次运行同一脚本，dpdispatcher 会根据记录恢复对已提交作业的跟踪，而不是重复提交。更复杂的自动纠错与工作流（custodian、atomate2、jobflow-remote、AiiDA）见[第6章](06-workflow-automation.md)。

---

## 参数选择与收敛建议

### 并行规模的选择流程

1. 先用**一个节点**、短时间（`NELM`/`NSW` 设小，或 CP2K 的 `MAX_SCF` 设小）做基准；
2. 依次测试 16、32、64 核（一个节点内），再测试 2、4 个节点；混合并行再测试每 rank 1、2、4 线程；
3. 比较每个电子步/离子步的平均耗时（VASP：`grep LOOP OUTCAR`；CP2K：输出末尾的 timing 表；Gaussian：`grep "Elapsed time" *.log`），计算并行效率 $E=T_\text{ref}N_\text{ref}/(T_N N)$；
4. 选效率 ≥ 70% 的最大规模；排队时间也算在内，“少核、快排上”往往比“多核、等一天”更快拿到结果。

### 经验速查表

| 情形 | 建议 |
|---|---|
| VASP 小原胞（< 20 原子）、k 点多 | 1 节点，`NCORE = 1–4`，KPAR 取能整除进程数和 k 点数的较大值 |
| VASP 100–300 原子、少量 k 点 | 1–4 节点，`NCORE` = NUMA 域核数，`KPAR` = 节点数（内存允许时） |
| VASP 只有 Γ 点 | `vasp_gam`，`KPAR = 1` |
| VASP NEB | 总进程数 = IMAGES × 每个 image 的进程数 |
| Gaussian | 1 节点，8–32 核；`%Mem` ≈ `--mem` 的 80%；本地 scratch |
| CP2K | 每节点 ranks × threads = 物理核数，先试每 rank 2 线程；`WALLTIME` 比时限短 10–20 分钟 |
| OpenMX | 周期体系 `-np` 可取到 自旋数 × k 点数 或更多，`-nt` 1–4（主要用来省内存）；先用小 `scf.maxIter` 测速 |
| 参数扫描 | 作业数组 + `%N` 并发上限；很多小任务打包到同一作业 |
| 长时间计算 | `--signal` + 体面停止 + 自动续投；不要申请远超需要的时限 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `sbatch: error: Batch script contains DOS line breaks (\r\n)` | 脚本在 Windows 下编辑 | `dos2unix job.sh`；编辑器设置 LF 换行（[第1章](01-linux-hpc.md)第10节） |
| `/bin/bash^M: bad interpreter: No such file or directory` | 同上，直接 `bash` 运行时出现 | 同上 |
| `sbatch: error: Batch job submission failed: Requested node configuration is not available` | 每节点核数、内存或 GPU 超过该分区节点的实际配置 | `sinfo -N -o "%N %c %m %P"` 查看真实配置，改 `--ntasks-per-node`、`--mem` |
| `sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified` | 没写或写错 `--account` | 用 `sacctmgr show assoc user=$USER format=account,partition` 查看可用账户 |
| `sbatch: error: Batch job submission failed: Job violates accounting/QOS policy (job submit limit, user's size and/or time limits)` | 超过分区时限、核数上限或排队作业数上限 | 缩短 `--time`、减少节点，或用数组 `%N` 控制数量 |
| `slurmstepd: error: *** JOB 123456 ON cn001 CANCELLED AT ... DUE TO TIME LIMIT ***` | 超过 `--time` | 增加时限；使用第 9 节的体面停止与续算 |
| `slurmstepd: error: Detected 1 oom_kill event ... Some of the step tasks have been OOM Killed.`（不同版本措辞略有差异），`sacct` 状态 `OUT_OF_MEMORY` | 内存不足 | 增加 `--mem` 或节点数；VASP 降低 KPAR、用 `vasp_gam`；Gaussian 降低 `%Mem` 使其留出余量 |
| `forrtl: severe (174): SIGSEGV, segmentation fault occurred` | 栈太小；MPI 与 OpenMP 混用不当；运行时模块与编译时不一致 | `ulimit -s unlimited`、`OMP_NUM_THREADS=1`（纯 MPI）、`OMP_STACKSIZE`；核对 `module list` |
| 输出重复 N 次，每个进程都认为自己是 rank 0 | `srun` 与 MPI 库的 PMI 不匹配，或 `mpirun` 属于另一个 MPI | 改用编译时的 MPI 的 `mpirun`；或 `srun --mpi=pmix`/`pmi2`（`srun --mpi=list` 查看），Intel MPI 设置 `I_MPI_PMI_LIBRARY` |
| `There are not enough slots available in the system to satisfy the N slots that were requested by the application`（Open MPI） | `mpirun -np` 大于作业分配的槽数 | 用 `-np ${SLURM_NTASKS}`；PBS 下检查 `mpiprocs` |
| 计算比预期慢十几倍，`top` 中负载远超核数 | 没设 `OMP_NUM_THREADS`，数学库在每个 rank 内开满线程 | 纯 MPI 设 `OMP_NUM_THREADS=1`；混合并行设为 `--cpus-per-task` |
| `M_divide: can not subdivide 64 nodes by 3`（VASP） | 进程数不能被 KPAR 或 IMAGES 整除 | 调整核数或 KPAR/IMAGES，使其整除 |
| Gaussian：`Erroneous write. Write -1 instead of ...` | scratch 所在磁盘写满或不可写 | `df -h $GAUSS_SCRDIR`；换到更大的本地盘；`MaxDisk` 限制 |
| Gaussian：`galloc: could not allocate memory.` | `%Mem` 大于作业实际可用内存 | `%Mem` 设为 `--mem` 的 75–85% |
| PBS 作业找不到输入文件 | 没有 `cd $PBS_O_WORKDIR`，作业在 `$HOME` 运行 | 脚本中加 `cd "${PBS_O_WORKDIR}"` |
| LSF 作业忽略了 `#BSUB` 设置 | 用了 `bsub job.sh` | 用 `bsub < job.sh` |
| 依赖作业一直排队，REASON 为 `DependencyNeverSatisfied` | 前一作业失败或被取消 | `scancel`；提交时加 `--kill-on-invalid-dep=yes` |
| 前一步 VASP 未收敛，后续作业却照常运行 | VASP 未收敛时退出码仍为 0 | 在脚本中检查 `reached required accuracy` 并 `exit` 非零值（第 8 节） |
| 作业在登录节点上运行（`sbatch` 后立刻有输出、没有作业号） | 用 `bash job.sh` 代替了 `sbatch job.sh` | 立即 `Ctrl-C`；用 `sbatch` 提交 |

---

## 练习

1. 在你的集群上运行 `sinfo -s`、`sinfo -N -o "%N %c %m %P" | head` 和 `lscpu`，写出：分区名、每节点物理核数、内存、NUMA 域数、是否开启超线程。据此修改 [slurm_mpi.sh](../examples/02-job-templates/slurm_mpi.sh) 头部，并用 `sbatch --test-only` 验证。
2. 把 [slurm_vasp.sh](../examples/02-job-templates/slurm_vasp.sh) 改写成 PBS Pro 版本和 LSF 版本（参考 [pbs_mpi.sh](../examples/02-job-templates/pbs_mpi.sh)、[lsf_mpi.sh](../examples/02-job-templates/lsf_mpi.sh)），保留其中的 NCORE/KPAR 检查。
3. 对一个约 50 原子的 VASP 体系，分别用 16、32、64 核（1 节点）和 2 节点运行 `NSW = 0`、`NELM = 10` 的单点计算，记录每个电子步的平均耗时，计算并行效率，并比较 `NCORE = 1, 4, 8, 16` 的差异。
4. 用 [02-array-chain-si](../examples/02-array-chain-si/README.md) 完成 ENCUT 扫描和弛豫→静态→能带依赖链；故意把 `INCAR.relax` 的 `NSW` 改为 1（只做一个离子步；实验晶格常数下 PBE 的应力不为零，`ISIF = 3` 时达不到收敛判据），观察依赖链如何被自动取消，并用 `sacct` 查看三个作业的状态与退出码。
5. 写一个 Gaussian 作业：`--cpus-per-task=16`、`--mem=40G`，`%NProcShared` 和 `%Mem` 按本章规则设置，`GAUSS_SCRDIR` 指向本地盘的作业目录；运行中登录计算节点（`srun --jobid=<id> --overlap --pty bash` 或站点允许的方式）用 `du -sh $GAUSS_SCRDIR` 观察 scratch 的增长。
6. 对 CP2K 的同一个输入，比较“64 ranks × 1 thread”“32 × 2”“16 × 4”三种组合的总耗时，并确认输出中 `Number of threads for this process` 与设置一致。

---

## 参考资料

官方文档：

- Slurm `sbatch` 手册：<https://slurm.schedmd.com/sbatch.html>
- Slurm `srun` 手册：<https://slurm.schedmd.com/srun.html>
- Slurm `sacct` 手册：<https://slurm.schedmd.com/sacct.html>
- Slurm 作业数组：<https://slurm.schedmd.com/job_array.html>
- Slurm MPI 用户指南（PMI-2/PMIx、`srun` 与各 MPI 库）：<https://slurm.schedmd.com/mpi_guide.html>
- Slurm CPU 绑定与多核支持：<https://slurm.schedmd.com/mc_support.html>
- Slurm 与 PBS/LSF 等的命令对照（Rosetta Stone）：<https://slurm.schedmd.com/rosetta.pdf>
- OpenPBS 文档：<https://openpbs.org/>（源码与发行说明：<https://github.com/openpbs/openpbs>）
- IBM Spectrum LSF 10.1 文档（`bsub`、`bjobs` 等命令参考）：<https://www.ibm.com/docs/en/spectrum-lsf/10.1.0>
- VASP Wiki：Parallelization、NCORE、KPAR、STOPCAR：<https://www.vasp.at/wiki/index.php/Category:Parallelization>，<https://www.vasp.at/wiki/index.php/NCORE>，<https://www.vasp.at/wiki/index.php/KPAR>，<https://www.vasp.at/wiki/index.php/STOPCAR>
- Gaussian 16 Link 0 命令（`%NProcShared`、`%CPU`、`%Mem`）：<https://gaussian.com/link0/>；运行环境与 `GAUSS_SCRDIR`：<https://gaussian.com/running/>
- CP2K 手册 `&GLOBAL`（`WALLTIME`）：<https://manual.cp2k.org/trunk/CP2K_INPUT/GLOBAL.html>；重启说明：<https://manual.cp2k.org/trunk/methods/restarting.html>
- OpenMX 4.0 手册第 28 节（并行与 `-nt`）：<https://www.openmx-square.org/openmx_man4.0/s28_parallel.html>，其中 28.1 MPI 并行：<https://www.openmx-square.org/openmx_man4.0/s28_1_mpi-parallel.html>，28.2 MPI/OpenMP 混合并行：<https://www.openmx-square.org/openmx_man4.0/s28_2_mpi-openmp-hybrid.html>
- Gaussian 16 `Default.Route` 与 `GAUSS_CDEF`/`GAUSS_PDEF`/`GAUSS_MDEF`：<https://gaussian.com/defroute/>
- dpdispatcher 文档：<https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/>

相关 GitHub 仓库（只链接，不复制代码）：

- [deepmodeling/dpdispatcher](https://github.com/deepmodeling/dpdispatcher) — 用 Python 生成并批量提交 Slurm/PBS/LSF 作业。
- [accre/SLURM](https://github.com/accre/SLURM) — Vanderbilt ACCRE 的 Slurm 示例脚本集（MPI、OpenMP、作业数组、GPU），可参考其写法。
- [dstdev/awesome-hpc](https://github.com/dstdev/awesome-hpc) — HPC 调度器、MPI、性能分析工具的 awesome 列表。
- [northword/wiki](https://github.com/northword/wiki) — 中文学习笔记，含 PBS 作业管理与 VASP 入门。
- [Matgenix/jobflow-remote](https://github.com/Matgenix/jobflow-remote) — 把 jobflow/atomate2 工作流提交到 Slurm/PBS 集群（[第6章](06-workflow-automation.md)）。

文献：

- G. M. Amdahl, "Validity of the single processor approach to achieving large scale computing capabilities", *AFIPS Conference Proceedings* **30**, 483 (1967). DOI: [10.1145/1465482.1465560](https://doi.org/10.1145/1465482.1465560)
- A. B. Yoo, M. A. Jette, M. Grondona, "SLURM: Simple Linux Utility for Resource Management", *Lecture Notes in Computer Science* **2862**, 44 (2003). DOI: [10.1007/10968987_3](https://doi.org/10.1007/10968987_3)
- T. D. Kühne et al., "CP2K: An electronic structure and molecular dynamics software package – Quickstep", *J. Chem. Phys.* **152**, 194103 (2020). DOI: [10.1063/5.0007045](https://doi.org/10.1063/5.0007045)

[← 上一章](01-linux-hpc.md) | [返回目录](../README.md) | [下一章 →](03-structure-modeling.md)
