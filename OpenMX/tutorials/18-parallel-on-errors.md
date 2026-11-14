# 第18章 并行、大体系 O(N)、作业脚本与常见错误排查（Parallelization, Large-Scale O(N) Methods, Job Scripts & Troubleshooting in OpenMX）

本章讲怎样把 OpenMX 4.0 高效地跑在集群上，以及出了问题怎样快速定位。内容包括：MPI × OpenMP 混合并行（`mpirun -np N openmx x.dat -nt M`、`OMP_NUM_THREADS`、`ulimit -s unlimited`），Slurm/PBS 作业模板与官方 deb 包的运行环境，`scf.eigen.lib elpa1|elpa2`，IO 与内存优化（`level.of.fileout 0`、`OutData.bin.flag on`、`memory.usage.fileout`）；大体系的 O(N) 线性标度方法 DC、DC-LNO、Krylov 与低阶标度方法，`orderN.HoppingRanges`、`orderN.KrylovH.order` 等关键词，以及 `-runtestL`/`-runtestL2` 基准测试；dc1394 维护的 GPU 版 openmx4.0_gpu、4.0 新增的 OpenMX API（C/Fortran/Python），`.scfout` 生态（TB2J、HopTB、DeepH/HamGNN、AiiDA、cif2x）。最后给出按编译、输入、SCF、优化、后处理、并行分类的常见错误速查表（OpenMX parallel performance, O(N) divide-and-conquer, troubleshooting）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出正确的混合并行命令，说清 `-np`、`-nt`、`OMP_NUM_THREADS`、核数之间的关系，避免超订（oversubscription）；
- 根据求解器（Cluster/Band/O(N)）和 k 点数选择 MPI 进程数，用基准测试找到本机最快的配置；
- 使用并修改 Slurm/PBS 作业模板，在作业结束后自动判断是否成功；
- 用 `level.of.fileout`、`OutData.bin.flag`、`memory.usage.fileout` 控制大体系的 IO 与内存；
- 判断什么时候该用 O(N) 方法，为 DC、DC-LNO、Krylov 选择 `orderN.HoppingRanges` 并验证精度，用“O(N) 自洽 + 一次对角化”得到大体系的本征值和 DOS；
- 了解 GPU 版、OpenMX API 与 `.scfout` 下游工具的用途和限制；
- 遇到报错时，先按速查表归类，再到对应章节找详细解释。

| 场景 | 本章对应小节 |
|---|---|
| 第一次在集群上提交 OpenMX | §3、§4 |
| 加了核反而变慢、节点负载异常 | §3.2–§3.4 |
| 计算正常但磁盘被 cube/重启文件写满 | §5 |
| 体系超过 1000 原子，常规对角化太慢 | §6 |
| 几百原子的分子/团簇，想用好 ELPA | §7 |
| 有 GPU 节点 | §8 |
| 想在自己的程序里调用 OpenMX、批量跑小计算 | §9 |
| 需要哈密顿量做 Jx、紧束缚、机器学习 | §10 |
| 碰到报错，不知道属于哪一类 | “常见问题与报错” |

本章示例：

- [18-job-scripts](../examples/18-job-scripts/README.md)：Slurm/PBS 模板、运行检查 `check_run.sh`、并行配置基准 `bench_np_nt.sh`；
- [18-on-si-dc-krylov](../examples/18-on-si-dc-krylov/README.md)：Si64 上 DC/DC-LNO/Krylov 的精度扫描，Si512 生产计算，O(N)+一次对角化；
- [18-c60-cluster-elpa](../examples/18-c60-cluster-elpa/README.md)：C60 的 Cluster 对角化（ELPA1/ELPA2）与 DC 对比；
- [18-api-python](../examples/18-api-python/README.md)：用 Python API 在一个 MPI 作业中同时跑多个输入。

> [!NOTE]
> 编译与安装（含 runtest 验证）见 [第2章](02-install.md)；SCF 参数与 `scf.restart` 见 [第5章](05-scf-kpoints-convergence.md)；集群与作业调度的通用知识见 [公共篇 作业调度](../../Common/tutorials/02-job-scheduler.md)。

---

## 2. 原理速览

### 2.1 计算量从哪里来

LCPAO 方法的核心是在每个 k 点求解广义本征值问题 $H\mathbf c=\varepsilon S\mathbf c$。矩阵维数 $N_b$ 等于全部原子的基函数个数之和（例如 `Si7.0-s2p1d1` 每个 Si 有 10 个基函数，512 个 Si 原子 $N_b=5120$）。常规对角化的耗时和内存分别为

$$
T_{\text{diag}}\propto N_k\,N_s\,N_b^{3},\qquad M\propto N_b^{2},
$$

$N_k$ 是不可约 k 点数，$N_s$ 是自旋重数（非自旋和非共线为 1，共线自旋极化为 2；非共线时矩阵维数翻倍）。一个复数稠密矩阵占 $16N_b^2$ 字节，$N_b=13000$ 时约 2.7 GB（由 ScaLAPACK 分布在各进程上）。哈密顿量矩阵元、实空间网格（Poisson、交换关联）等部分则大致随原子数线性增长。因此几百原子以内对角化与网格部分相当，上千原子后对角化成为瓶颈。

### 2.2 OpenMX 的并行层次

手册第 28 章给出的实现要点：

- **Cluster 计算**：MPI 先分配到自旋，再分配到本征态；本征值问题由 ELPA 求解（`scf.eigen.lib elpa1|elpa2`，Cluster 默认 `elpa1`）。
- **Band 计算**：MPI 依次分配到自旋、k 点、本征态；4.0 中 Band 默认 `elpa2`。当进程数 ≤ 自旋重数 × k 点数时，每个 k 点要对角化两次（第一次找化学势、第二次直接求密度矩阵，以节省存波函数的内存）；进程数更多时只需对角化一次，因此并行加速比可能“超线性”。
- **O(N) 计算**：天然按原子并行，DC-LNO 还有原子、自旋、对角化三层并行。
- **数据分布**：4.0 使用三维区域分解，并行效率与晶胞矢量的取法基本无关；进程数也不再受原子数限制。
- **OpenMP**：`-nt M` 指定每个 MPI 进程的线程数，不写时为 1。手册指出混合并行在内存上更省，速度与纯 MPI 接近，最优选择依赖机器，建议用基准测试确认。

> [!NOTE]
> 与 3.9 的差别：Band 计算默认 `elpa2` 是 4.0 起的改动（手册第 28.1.3 节“Since OpenMX Ver. 4.0, the default choice has been ELPA2”）；默认 k 网格 `scf.Generation.Kpoint GCenter`、`scf.restart pp`、OpenMX API 也都是 4.0 新增，3.9 中没有。用 3.9 复现文献结果或做基准对比时，请显式写出 `scf.eigen.lib` 和 k 网格类型。

共线计算的不可约 k 点数（手册第 28.1.3 节；非共线时 $N_k=N_1N_2N_3$）：

$$
N_k^{\text{GCenter}}=\frac{N_1N_2N_3+2^{m}}{2},\qquad
N_k^{\text{Regular}}=\begin{cases}\dfrac{N_1N_2N_3-1}{2}+1 & N_1,N_2,N_3\ \text{全为奇数}\\[2mm]\dfrac{N_1N_2N_3}{2} & \text{至少一个为偶数}\end{cases}
$$

$m$ 是 $(N_1,N_2,N_3)$ 中偶数的个数。例如共线、`scf.Kgrid 4 4 4`、默认 GCenter 时 $N_k=(64+8)/2=36$。

并行加速受串行部分（磁盘 IO、负载不均）限制。设可并行部分占比为 $f$，Amdahl 定律给出 $p$ 个进程的加速比上限

$$
S(p)=\frac{1}{(1-f)+f/p}.
$$

手册中 512 原子金刚石 DC 计算用 128 进程加速约 84 倍，原因正是每个进程分到的原子变少后，IO 等非并行部分占比上升。

### 2.3 O(N) 方法：把大问题拆成许多小团簇

**分治法（DC）**：对每个原子 $i$，取以它为中心、半径为 `orderN.HoppingRanges` 的球内原子构成截断团簇，求解团簇的本征值问题 $H^{(i)}\mathbf c^{(i)}_\mu=\varepsilon^{(i)}_\mu S^{(i)}\mathbf c^{(i)}_\mu$，只保留投影到中心原子上的部分：

$$
D_i(E)=\sum_{\mu}\sum_{\alpha\in i}\sum_{j\beta}c^{(i)*}_{\mu,i\alpha}\,S_{i\alpha,j\beta}\,c^{(i)}_{\mu,j\beta}\,\delta\!\left(E-\varepsilon^{(i)}_\mu\right),
\qquad
N_e=\sum_i\int D_i(E)\,f\!\left(\frac{E-\mu}{k_BT}\right)dE .
$$

第二式确定共同的化学势 $\mu$。每个团簇的大小与体系总原子数无关，团簇个数等于原子数 $N$，所以总耗时 $\propto N$。误差来自截断：有能隙的固体和分子随团簇增大近似指数收敛，金属收敛较慢（手册图 23）。

**DC-LNO**：在团簇外层用局域自然轨道（LNO）代替 PAO，降低团簇矩阵维数。LNO 由矩阵 $\Lambda_{0i}=\sum_{\mathbf Rj}\rho_{0i,\mathbf Rj}S_{\mathbf Rj,0i}$ 的本征向量给出，本征值可理解为 LNO 的占据数，保留占据数 ≥ `orderN.LNO.Occ.Cutoff` 的 LNO。手册在三种 O(N) 方法中**推荐 DC-LNO**：稳健、适用于包括金属在内的多种材料，并行效率最高。

**Krylov 子空间法**：把团簇的向量空间映射到维数为 $q$（`orderN.KrylovH.order`）的 Krylov 子空间

$$
\mathcal K_q=\operatorname{span}\left\{|K_0\rangle,\ (S^{-1}H)|K_0\rangle,\ \dots,\ (S^{-1}H)^{q-1}|K_0\rangle\right\},
$$

再在子空间里对角化。对金属比 DC 更高效，对共价/离子体系耗时相当（手册图 28）。

**数值精确的低阶标度方法**：不做截断近似，对一维、二维、三维体系的耗时分别为 $O(N(\log N)^2)$、$O(N^2)$、$O(N^{7/3})$，只支持共线的 Cluster 或 Γ 点周期计算，前因子大，只有在用很多核时才有优势（手册第 48 章）。

> [!IMPORTANT]
> O(N) 方法给出的是自洽密度、总能和力，**不提供波函数**。能带、分子轨道 cube、精细的 DOS 需要在 O(N) 自洽密度上再做一次常规对角化（§6.5）。

---

## 3. MPI × OpenMP 混合并行

### 3.1 运行命令与环境变量

```bash
export OMP_NUM_THREADS=2          # must equal -nt
export OMP_STACKSIZE=512M         # private stack of each OpenMP thread
ulimit -s unlimited               # stack of the main thread
mpirun -np 32 openmx si512.dat -nt 2 > si512.std 2>&1
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| `-np N` | 见 §3.2 | MPI 进程数 |
| `-nt M` | 1、2 或 4 | 每个 MPI 进程的 OpenMP 线程数，不写时为 1 |
| `OMP_NUM_THREADS` | 与 `-nt` 相同 | MKL/OpenBLAS 等库按它开线程；不设置时可能每个进程开满全部核 |
| `OMP_STACKSIZE` | `512M` | OpenMP 线程私有栈；大体系、OpenMP 运行出现段错误时尤其重要 |
| `ulimit -s` | `unlimited` | 主线程栈大小；写在作业脚本里，在 `mpirun` 之前执行 |
| `OMP_PLACES` / `OMP_PROC_BIND` | `cores` / `close` | 线程绑定到相邻物理核 |
| `I_MPI_PIN_DOMAIN`（Intel MPI） | `omp` | 每个进程的绑定域大小等于 `OMP_NUM_THREADS` |
| Open MPI 绑定 | `--map-by ppr:P:node:PE=M --bind-to core` | 每节点 P 个进程，每个进程 M 个核 |

**硬性规则**：每个节点上 `进程数 × 线程数 ≤ 物理核数`。超线程（SMT）一般不计入。

### 3.2 进程数怎么选

| 求解器 | 建议 | 理由 |
|---|---|---|
| `Band`（周期体系） | 进程数取 `自旋重数 × N_k` 的约数或倍数；进程多于 `自旋重数 × N_k` 时每个 k 点只对角化一次 | k 点并行几乎无通信，效率最高 |
| `Cluster`（分子） | 中等进程数 + 多线程；非自旋极化小分子不宜用太多进程 | 先分自旋、再分本征态，后者通信较多 |
| O(N)（DC/DC-LNO/Krylov） | 进程数可以接近原子数；DC-LNO 可以更多（原子 × 自旋 × 对角化） | 按原子天然并行；负载不均时效率下降 |
| `NEGF` | 见 [第17章](17-negf-transport.md) 与手册第 46.8 节 | 能量点与 k 点并行 |
| NEB | 见 [第6章 §7.5](06-optimization-eos-neb.md) | 按映像分组 |

内存不够时，首先减少每个节点的 MPI 进程数、增大 `-nt`（总核数不变），因为很多数据在每个进程中各存一份；其次增加节点数。

### 3.3 用基准测试找到最快配置

手册明确建议用基准测试确定纯 MPI 与混合并行的取舍。做法：取实际要算的体系，把 `scf.maxIter` 改成 5–10，在同一节点上扫描不同组合：

```bash
cd OpenMX/examples/18-c60-cluster-elpa
sed 's/^scf.maxIter .*/scf.maxIter                 10/' C60_cluster_elpa.dat > c60_short.dat
bash ../18-job-scripts/bench_np_nt.sh c60_short.dat "16x1 8x2 4x4 2x8"
column -s, -t bench_c60_elpa1.csv
```

`bench_np_nt.sh` 为每个组合改写 `System.Name` 后运行，从 `.out` 末尾的耗时表中取 `Total Computational Time`（以及 `Diagonalization` 分项，行名以你的版本输出为准），写入 CSV。判断标准：

1. 所有正常结束的组合 `Utot` 必须一致（相差 1e-8 Hartree 量级）；
2. 同一核数下取总时间最短者；
3. 再在 2 倍、4 倍核数下重复，确认加速比是否值得多占用资源。

官方的 `-runtest` 也可以当作简易基准：它在 `runtest.result` 末尾给出总耗时，见 [第2章 §7.2](02-install.md)。

### 3.4 检查实际运行状态

```bash
ssh <compute-node>
top -H -u $USER                     # threads: each openmx rank should show NT busy threads
ps -o pid,psr,pcpu,comm -u $USER | grep openmx    # psr = core the process runs on
grep -A30 "Computational Time" si512.out          # time spent in each part of the code
```

- 每个 `openmx` 进程的 CPU 占用约为 `NT × 100%`；若远超这个值，多半是 BLAS 自行开了线程（`OMP_NUM_THREADS` 没设）；
- 若所有进程挤在少数几个核上（`psr` 相同），是绑定选项错误，检查 `I_MPI_PIN_DOMAIN` 或 `--bind-to`；
- `.out` 末尾的耗时表列出各部分（对角化、哈密顿量构造、Poisson、混合、力等）的最小/最大耗时，最大值远大于最小值说明负载不均。

---

## 4. 作业脚本：Slurm 与 PBS

### 4.1 Slurm 模板

示例目录 [18-job-scripts](../examples/18-job-scripts/README.md) 中的 `job_slurm.sh`（完整内容）：

```bash
#!/bin/bash
#SBATCH --job-name=openmx
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32          # MPI processes per node
#SBATCH --cpus-per-task=2             # OpenMP threads per MPI process (= -nt)
#SBATCH --time=24:00:00
#SBATCH --partition=cpu               # site specific
#SBATCH --output=%x-%j.log
#
# OpenMX 4.0 hybrid MPI x OpenMP job (Slurm).
#   sbatch job_slurm.sh input.dat
# Total cores = nodes x ntasks-per-node x cpus-per-task (here 2 x 32 x 2 = 128).
#
# Environment switches (export before sbatch, or edit here):
#   OPENMX=/path/to/openmx       executable (default: openmx in PATH)
#   USE_DEB=1                    official Debian package in /opt/openmx4.0
#   LAUNCH=srun                  use srun instead of mpirun (MPI built with Slurm PMI support)
#   KEEP_SCRATCH=0               remove *_rst/ after a successful run
set -u
INPUT=${1:?usage: sbatch job_slurm.sh input.dat}
[ -f "${INPUT}" ] || { echo "input ${INPUT} not found"; exit 1; }

# ---------- environment ----------
# module purge
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl      # Intel toolchain
# module load gcc openmpi openblas scalapack fftw                           # GNU toolchain
if [ "${USE_DEB:-0}" = "1" ]; then
    # runtime libraries of the official .deb (fixes 'libmpifort.so.12: cannot open shared object file')
    export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi
OPENMX=${OPENMX:-openmx}

NP=${SLURM_NTASKS:-1}
NT=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${NT}          # keep identical to -nt (MKL/OpenBLAS threads follow OpenMP)
export OMP_STACKSIZE=512M             # large private stacks for OpenMP threads
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited                   # avoid segmentation faults from deep stacks
# Intel MPI: one pinning domain per MPI process, sized by OMP_NUM_THREADS
export I_MPI_PIN_DOMAIN=omp

NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
STD=${NAME}.std
echo "host $(hostname)  NP=${NP}  NT=${NT}  input=${INPUT}  name=${NAME}  start $(date)"

# ---------- run ----------
if [ "${LAUNCH:-mpirun}" = "srun" ]; then
    # some Slurm versions (e.g. 22.05-23.02) do not pass --cpus-per-task from sbatch to srun:
    # give it explicitly, which is safe for every version
    export SRUN_CPUS_PER_TASK=${NT}
    srun --cpus-per-task="${NT}" --cpu-bind=cores "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1
else
    # Open MPI binding (ignored by Intel MPI, which uses I_MPI_PIN_DOMAIN):
    #   mpirun --map-by ppr:${SLURM_NTASKS_PER_NODE}:node:PE=${NT} --bind-to core ...
    mpirun -np "${NP}" "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1
fi

# ---------- check (do not rely on the exit code: grep the final message) ----------
if grep -q "The calculation was normally finished." "${STD}"; then
    echo "finished normally $(date)"
    grep -h "Utot\." "${NAME}.out" | tail -n 1
    grep -h "Total Computational Time" "${NAME}.out" | tail -n 1
    if [ "${KEEP_SCRATCH:-1}" = "0" ]; then rm -rf "${NAME}_rst"; fi
else
    echo "OpenMX did NOT finish normally; last lines of ${STD}:" >&2
    tail -n 30 "${STD}" >&2
    exit 1
fi
```

| 项目 | 推荐值 | 说明 |
|---|---|---|
| `--ntasks-per-node × --cpus-per-task` | = 每节点物理核数 | 例如 64 核节点：32×2、16×4、64×1 |
| `-nt` | `$SLURM_CPUS_PER_TASK` | 由调度器的分配推出，避免手写不一致 |
| `System.Name` 提取 | `awk` 不区分大小写 | OpenMX 关键词不区分大小写，`.std` 用 System.Name 命名便于 `check_run.sh` 对应 `.out` |
| 成功判据 | `.std` 中有 `The calculation was normally finished.` | 不依赖程序退出码；正常结束也不等于 SCF 收敛，还要看 `.out`（[第3章](03-input-first-calc.md)） |
| `KEEP_SCRATCH` | 默认保留 `_rst/` | 需要续算、`scf.restart pp` 补算能带/DOS 时必须保留 |

提交与检查：

```bash
sbatch job_slurm.sh Si512_DC-LNO.dat
squeue -u $USER
bash check_run.sh si512_dclno.std
```

`check_run.sh` 对每个 `.std` 报告：是否正常结束、匹配到的已知报错（如 `Found an undefined species`、`ERROR key=`、`Segmentation fault`、`error while loading shared libraries`、`Invalid communicator`）、最后一个 SCF 块的步数与 `NormRD`、`|dUele|`、`Utot`、总耗时，有问题时返回 1。

### 4.2 PBS 模板

```bash
#PBS -l select=2:ncpus=64:mpiprocs=32:ompthreads=2
```

PBS Pro/OpenPBS 中 `mpiprocs` 是每个 chunk 的 MPI 进程数，`ompthreads` 会被写入 `OMP_NUM_THREADS`；`$PBS_NODEFILE` 每行对应一个 MPI 进程，`job_pbs.sh` 据此计算 `-np`，并以 `OMP_NUM_THREADS` 作为 `-nt`。提交：`qsub -v INPUT=input.dat job_pbs.sh`。Torque 的写法（`-l nodes=2:ppn=64`）见脚本头部注释，此时线程数要自己设置，并保证 `$PBS_NODEFILE` 行数 × 线程数不超过分配的核数。Slurm 与 PBS 的资源写法对照见 [公共篇 作业调度](../../Common/tutorials/02-job-scheduler.md)。仓库还提供可放进 `PATH` 长期使用的通用模板 [job_omx_slurm.sh](../scripts/job_omx_slurm.sh)、[job_omx_pbs.sh](../scripts/job_omx_pbs.sh)，用 `bash ../scripts/job_omx_slurm.sh -h` 查看用法。

### 4.3 官方 deb 包的运行环境

deb 包安装在 `/opt/openmx4.0`，自带运行 OpenMX 所需的 Intel oneAPI 运行库（[第2章 §8](02-install.md)、手册第 4.11 节）。它面向 Ubuntu/Debian/WSL 的个人电脑与工作站（安装需要 `sudo apt`）；共享集群上一般没有 root 权限，仍应按第2章从源码编译。在装了 deb 包的机器上写脚本时需要：

```bash
export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
command -v openmx DosMain               # call the installed commands by name, not ./openmx
DosMain si512_dclno.Dos.val si512_dclno.Dos.vec
```

并在输入中写 `DATA.PATH /opt/openmx4.0/DFT_DATA19`（以实际安装位置为准，目录下应同时有 `PAO/` 和 `VPS/`）。缺少 `LD_LIBRARY_PATH` 时，`DosMain` 等后处理工具会报 `error while loading shared libraries: libmpifort.so.12: cannot open shared object file`。手册第 4.11 节说明：deb 包装好后即可单进程运行；**在另外安装了 Intel oneAPI 的环境中**才能做 MPI 并行（`mpirun -np 2 openmx ./Methane.dat`），建议先用单进程确认安装，再上并行。`job_slurm.sh` 中用 `USE_DEB=1` 打开这一设置。

### 4.4 长作业：续算、在线调参与墙钟时间

- **SCF 续算**：保留 `System.Name_rst/`，下次运行设 `scf.restart on`（[第5章 §8](05-scf-kpoints-convergence.md)）；
- **在线调参/停止**：运行目录下建立 `<System.Name>_SCF_keywords` 文件即可修改混合参数或停止 SCF（[第5章 §7](05-scf-kpoints-convergence.md)），适合在墙钟时间到之前体面地结束；
- **优化/MD 续算**：每步更新的 `输入文件名#` 含最新坐标，改名后接着算（[第6章 §3.6](06-optimization-eos-neb.md)、[第7章](07-md.md)）；
- **墙钟时间**：大体系先用 `scf.maxIter 3` 的短作业估算单步耗时，再乘以预计步数并留 30% 余量。

---

## 5. IO 与内存优化

大体系最容易被忽视的瓶颈是磁盘：每个 cube 文件的大小与实空间网格点数成正比，几千原子时单个 cube 可达数 GB。

```text
level.of.stdout          1
level.of.fileout         0        # 0: minimal output, no cube/grid files; 1: standard (default); 2: more
OutData.bin.flag         on       # large files in binary form (default off)
memory.usage.fileout     on       # write System.Name.memory0, .memory1, ... (default off)
HS.fileout               off      # .scfout only when a post-processing code needs it
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `level.of.fileout` | 大体系、MD、优化用 `0` | 0 不写 cube 和网格文件；1 为默认标准输出；2 输出更多 |
| `OutData.bin.flag` | IO 慢的机器上 `on` | 大文件以二进制写出；之后用 `work/bin2txt.c` 转换：`gcc bin2txt.c -lm -o bin2txt`，`./bin2txt *.bin` |
| `memory.usage.fileout` | 估计内存时 `on` | 每个 MPI 进程写一个 `*.memory#` 文件，末行 `Memory: total` 为该进程统计到的总量；手册说明清单并不完整 |
| `HS.fileout` | 按需 | `.scfout` 含 H、S、密度矩阵等，大体系文件很大 |
| `System.Name_rst/` | 生产计算后可删除 | 需要 `scf.restart` 时保留；`KEEP_SCRATCH=0` 自动删除 |

内存估计的实用流程：

```bash
# short run with memory.usage.fileout on and scf.maxIter 1
grep "Memory: total" si512_dclno.memory*        # per-process totals (incomplete list)
sacct -j <jobid> --format=JobID,MaxRSS,Elapsed  # Slurm: real peak memory per task
```

把单进程峰值乘以每节点进程数，与节点内存比较；超过 80% 时减少每节点进程数（增大 `-nt`）或增加节点。

> [!TIP]
> 大体系先用 `level.of.fileout 0` 完成自洽和优化，需要电荷密度或势的 cube 时，再在最终结构上用 `scf.restart pp`（[第5章 §8](05-scf-kpoints-convergence.md)）设 `level.of.fileout 1` 补算一次。

---

## 6. 大体系 O(N)：DC / DC-LNO / Krylov

### 6.1 什么时候用 O(N)

| 体系规模与条件 | 推荐求解器 | 依据 |
|---|---|---|
| ≤ 几百原子 | `Band`（周期）/`Cluster`（分子） | 精确、可得波函数 |
| ~1000 原子，有数百核 | 仍可用常规对角化 | 手册 `-runtestL2`：1000 原子级体系每个 SCF 步约 7–17 s（192 进程），结构优化可行 |
| > 1000 原子，或需要长时间 MD | `DC-LNO`（首选） | 手册推荐；金属、绝缘体都适用，需团簇约 200–300 原子 |
| 共价/有能隙体系，团簇可以较小 | `DC` | 稳定，误差随团簇增大近似指数下降 |
| 简单金属 | `Krylov` 或 `DC-LNO` | Krylov 对金属比 DC 高效 |
| 低维（一维/二维）大体系，共线，Γ 点，核很多 | 低阶标度方法 | 数值精确，只有在核数多时才有优势 |

手册给出的经验：DC 在约 500 原子以上才可能比常规对角化快（依赖体系与核数）。所以 O(N) 不是“越早用越好”，而是要先用小超胞验证精度，再在大体系上用。

### 6.2 第一步：小超胞上的精度测试（Si64）

示例目录 [18-on-si-dc-krylov](../examples/18-on-si-dc-krylov/README.md) 用 64 原子 Si 超胞，把三种 O(N) 方法与 4×4×4 k 点的常规计算比较。DC 输入 `Si64_DC.dat`（坐标块共 64 行，此处只列出前 3 行和最后 1 行）：

```text
# Si64 (2x2x2 conventional diamond supercell, a = 5.43 A), Quick basis Si7.0-s2p1d1.
# O(N) divide-conquer (DC). Scan orderN.HoppingRanges with scan_hopping.sh. OpenMX 4.0.

System.CurrrentDirectory    ./
System.Name                 si64_dc
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: dir containing PAO/ and VPS/
level.of.stdout             1
level.of.fileout            0          # no cube files for large systems

Species.Number              1
<Definition.of.Atomic.Species
  Si   Si7.0-s2p1d1   Si_PBE19
Definition.of.Atomic.Species>

Atoms.Number                64
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
     1  Si  0.000000  0.000000  0.000000   2.0 2.0
     2  Si  0.000000  0.250000  0.250000   2.0 2.0
     3  Si  0.250000  0.000000  0.250000   2.0 2.0
    64  Si  0.875000  0.875000  0.625000   2.0 2.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  10.8600   0.0000   0.0000
   0.0000  10.8600   0.0000
   0.0000   0.0000  10.8600
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            200.0
scf.maxIter                 100
scf.EigenvalueSolver        DC
orderN.HoppingRanges        7.0        # Ang
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          15
scf.Mixing.StartPulay       8
scf.criterion               1.0e-8

MD.Type                     Nomd
```

另外三个输入只有求解器部分不同：

```text
# Si64_Band_ref.dat (reference)
scf.EigenvalueSolver        Band
scf.Kgrid                   4 4 4
scf.eigen.lib               elpa2

# Si64_DC-LNO.dat
scf.EigenvalueSolver        DC-LNO
orderN.HoppingRanges        8.0
orderN.LNO.Buffer           0.2
orderN.LNO.Occ.Cutoff       0.1
scf.dclno.threading         on

# Si64_Krylov.dat
scf.EigenvalueSolver        Krylov
orderN.HoppingRanges        7.0
orderN.KrylovH.order        400
orderN.Exact.Inverse.S      on
orderN.Recalc.Buffer        on
orderN.Expand.Core          on
```

运行扫描（每种方法 R = 7, 8, 9, 10, 11 Å）：

```bash
cd OpenMX/examples/18-on-si-dc-krylov
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat
NP=16 NT=1 bash scan_hopping.sh
column -s, -t hopping_scan.csv
```

`scan_hopping.sh` 先算参考值，再为每个 R 生成 `si64_dc_R7.dat` 等输入并运行，把 $\Delta E=(U_{\text{tot}}^{O(N)}-U_{\text{tot}}^{\text{Band}})/64$ 换算为 meV/atom。输出格式如下（示意输出，数值仅供参考）：

```text
solver   R_ang  Utot_Hartree  dE_meV_per_atom  wall_s  FNAN_SNAN_line
DC       7      ...           ...              ...     ...
DC       9      ...           ...              ...     ...
```

成功判据：

- 所有计算正常结束且 SCF 收敛（`grep "normally finished" *.std`，并看 `.out` 中 SCF 历史）；
- |ΔE| 随 R 增大而减小，对 Si 这类有能隙体系近似指数下降（手册图 23）；
- 取 |ΔE| 低于你的容差（如 1–5 meV/atom，取决于你要比较的能量差）的最小 R；
- 若要算力或做 MD，再比较 `.out` 中的原子力，力的误差同样要满足要求。

> [!NOTE]
> 为什么 O(N) 计算不写 `scf.Kgrid`？截断团簇在实空间中按周期性构建，已包含周期像原子，相当于在实空间完成了布里渊区积分；对应的常规计算必须用收敛的 k 点作为参考，不能拿 Γ 点结果比较。

### 6.3 关键词详解

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `scf.EigenvalueSolver` | `DC` / `DC-LNO` / `Krylov` | 手册第 8.2 节的选项：`DC`、`DC-LNO`、`Krylov`、`ON2`（低阶标度）、`Cluster`、`Band`（另有 NEGF 专用的 `NEGF`） |
| `orderN.HoppingRanges` | DC：6.0–7.0 Å（分子，手册建议）；DC-LNO：使团簇含 200–300 原子 | 以每个原子为中心的截断球半径（Å）。三种 O(N) 方法共用 |
| `orderN.LNO.Buffer` | 0.1–0.2（默认 0.2） | 外层分配 LNO 的范围。0.0 表示所有 SNAN 原子都用 LNO；1.0 表示都用 PAO，等同 DC |
| `orderN.LNO.Occ.Cutoff` | 0.1（默认） | 保留占据数 ≥ 此值的 LNO |
| `<LNOs.Num` | 例如 `Si 4`、`H 1` | 直接指定每种元素的 LNO 个数；**优化和 MD 推荐使用**，避免 LNO 个数在过程中变化 |
| `scf.dclno.threading` | 用 `-nt > 1` 时 `on` | DC-LNO 最底层（团簇对角化）用 OpenMP 并行 |
| `orderN.KrylovH.order` | 例如 400（手册示例） | 哈密顿量 Krylov 子空间维数 |
| `orderN.Exact.Inverse.S` | 共价体系 `on`（默认） | `on`：每个团簇严格求 $S^{-1}$；`off`：用 Krylov 近似，维数由 `orderN.KrylovS.order` 给出（默认为 `KrylovH.order` 的 4 倍） |
| `orderN.Recalc.Buffer` | `on`（默认） | 每个 SCF 步重算缓冲矩阵；`off` 只在第一步计算 |
| `orderN.Expand.Core` | 共价体系 `on`（默认） | `on`：核心区为半径 $1.2r_{\min}$ 内的原子（$r_{\min}$ 为最近邻距离）；简单金属设 `off` 可能收敛更好 |
| `<orderN.FNAN+SNAN` | 变晶胞计算时使用 | 每个原子直接给定团簇原子数，见 §6.6 |

DC-LNO 的 `<LNOs.Num` 块写法（物种名与 `Definition.of.Atomic.Species` 一致；块内不要写注释）：

```text
<LNOs.Num
  Si   4
LNOs.Num>
```

> [!WARNING]
> 团簇中的原子分两类：基函数与中心原子有重叠的近邻（FNAN）和由 `orderN.HoppingRanges` 额外纳入的原子（SNAN），两者的个数会写在标准输出中（`grep -i SNAN x.std`）。`orderN.HoppingRanges` 小于基函数重叠距离时，团簇实际上只由 FNAN 决定，再减小 R 不会改变结果。`Si7.0` 的截断半径为 7.0 bohr（3.70 Å），两个 Si 原子在 7.4 Å 以内就有重叠，所以本例的扫描从 7 Å 开始。

### 6.4 第二步：Si512 生产计算

`Si512_DC.dat`、`Si512_DC-LNO.dat` 与 Si64 输入的写法完全相同，只是超胞为 4×4×4（21.72 Å，512 行坐标），`orderN.HoppingRanges` 按扫描结果填写（文件中预设 8.0 Å，请修改）。更大的超胞用示例中的 `make_si_supercell.py` 生成（`python3 make_si_supercell.py -h` 查看用法），例如 8×8×8 共 4096 原子。

```bash
sbatch ../18-job-scripts/job_slurm.sh Si512_DC-LNO.dat
bash ../18-job-scripts/check_run.sh si512_dclno.std
grep "Utot\." si512_dclno.out
awk '/Total Computational Time/{print $NF}' si512_dclno.out
```

检查要点：

- 理想晶体中每个原子等价，O(N) 方法的截断团簇也相同，所以 `Utot/512` 应与 Si64 在相同 R 下的 `Utot/64` 非常接近（剩余差异来自 SCF 收敛精度和两个晶胞自动选取的实空间网格不完全对应；需要严格比较时用 `scf.Ngrid` 让两者网格间距相同）——这是检查超胞输入是否写对的简单办法；
- 用 `bench_np_nt.sh` 或不同节点数的短作业，记录单步 SCF 耗时随核数和原子数的变化，验证近似线性的标度。

### 6.5 第三步：O(N) 自洽 + 一次常规对角化

手册第 29.2 节的做法：先用 O(N) 得到自洽密度，再在这个密度上做一次常规对角化，得到波函数和本征值。`Si512_onestep_band.dat` 与 `Si512_DC-LNO.dat` 的区别只有：

```text
System.Name                 si512_dclno          # must be the same as the O(N) run
scf.maxIter                 1          # one diagonalization only
scf.EigenvalueSolver        Band
scf.Kgrid                   1 1 1
scf.restart                 on         # read si512_dclno_rst/

Dos.fileout                 on
Dos.Erange                  -20.0  10.0
Dos.Kgrid                   1 1 1
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `System.Name` | 与 O(N) 计算相同 | `scf.restart` 读取 `System.Name_rst/` |
| `scf.maxIter` | 1 | 只对角化一次 |
| `scf.EigenvalueSolver` | `Band`（周期）/`Cluster`（分子） | 常规对角化 |
| `scf.Kgrid` | 大超胞 `1 1 1` | 512 原子超胞的 Γ 点相当于 8 原子胞的 4×4×4 网格 |
| `Dos.fileout` / `MO.fileout` | 按需 | 手册示例用 `MO.fileout on`、`num.HOMOs 2`、`num.LUMOs 2`、`MO.Nkpoint 1` 输出 Γ 点 HOMO/LUMO |

```bash
cp si512_dclno.out si512_dclno_ON.out          # the second run overwrites the .out file
mpirun -np 64 openmx Si512_onestep_band.dat -nt 1 > si512_onestep.std 2>&1
ls si512_dclno.Dos.val si512_dclno.Dos.vec
```

之后用 `DosMain` 处理（[第9章](09-dos-pdos-cohp-optics.md)）。只有 Γ 一个 k 点时四面体法没有意义，`DosMain` 中选 Gaussian Broadening（2）并给一个 0.1–0.2 eV 的展宽。4.0 的 `scf.restart pp`（自动 `scf.maxIter=1`、`MD.Type=nomd`，且不覆盖重启文件）是同类功能（[第5章 §8](05-scf-kpoints-convergence.md)）；手册在第 29.2 节演示的是 `scf.restart on` + `scf.maxIter 1` 的写法，本例沿用它。

### 6.6 变晶胞计算：固定团簇大小 `<orderN.FNAN+SNAN`

用 `orderN.HoppingRanges` 定义团簇时，改变晶格常数会让球内原子数突变，E–V 曲线出现台阶（与 [第5章 §4.3](05-scf-kpoints-convergence.md) 的网格跳变类似）。此时改为直接给出每个原子的团簇原子数，行数必须等于原子数：

```text
<orderN.FNAN+SNAN
  1  80
  2  80
  3  80
  4  80
orderN.FNAN+SNAN>
```

第一列是 `Atoms.SpeciesAndCoordinates` 中的序号，第二列是 FNAN+SNAN；程序会为每个原子调整球半径使团簇原子数等于给定值。上面只是 4 原子体系的写法示意，实际数值先从平衡结构计算的标准输出中读出。同时用 `scf.Ngrid` 固定实空间网格。

### 6.7 低阶标度方法与 `-runtestL` / `-runtestL2` 基准

**低阶标度方法**：手册第 8.2 节的关键词列表把它写作 `scf.EigenvalueSolver ON2`，而第 48 章正文与示例使用 `cluster2`；极点数由 `scf.Npoles.ON2` 控制（电子温度高于 300 K 时 100 个极点已足够）。使用前先查看自己版本附带的示例：

```bash
grep -i "EigenvalueSolver\|Npoles" ~/openmx4.0/work/C60_LO.dat
```

以该文件中的写法为准。手册估计在 100 核以上时，与常规方法的耗时交叉点约在 300 原子（依赖维度）；它只支持共线的 Cluster 或 Γ 点计算，不支持 DOS 和波函数输出。

**大体系基准**：

```bash
cd ~/openmx4.0/work
mpirun -np 112 openmx -runtestL -nt 2      # 16 inputs in work/large_example -> runtestL.result
mpirun -np 128 openmx -runtestL2 -nt 4     # 7 inputs (~1000 atoms) in work/large2_example
```

`runtestL` 的 16 个输入包括 DC-LNO（`CG15c-DC-LNO.dat`）、DC 金刚石（`DIA512-1.dat`）、Mn12、MCCN 等；`runtestL2` 包括 C1000、Fe1000、GRA1024、Ih-Ice1200、Pt500、R-TiO2-1050、Si1000。`diff Utot`、`diff Force` 在小数点后第 7 位以内即为正确；`Elapsed time` 用来比较机器、编译选项和并行配置。手册列出的 `runtestL2` 参考结果（192 MPI × 1 线程，AMD 集群）总耗时约 11300 s，单个 SCF 步约 7–17 s，可以作为你自己集群的对照。核数少时这些测试可能因内存不足而段错误，只在集群上运行。

---

## 7. 常规对角化的大体系：Cluster + ELPA（C60 示例）

示例 [18-c60-cluster-elpa](../examples/18-c60-cluster-elpa/README.md) 用理想截角二十面体 C60（`C6.0-s2p2d1`，780 个基函数）比较 ELPA1、ELPA2 和 DC。`C60_cluster_elpa.dat` 的关键部分：

```text
System.Name                 c60_elpa1
level.of.fileout            0          # no cube files: this is a timing/accuracy test

Species.Number              1
<Definition.of.Atomic.Species
  C    C6.0-s2p2d1   C_PBE19
Definition.of.Atomic.Species>

Atoms.Number                60
Atoms.SpeciesAndCoordinates.Unit   Ang
# 60 coordinate lines (C-C = 1.44 A, centred in an 18 A box), see the example file

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            200.0
scf.maxIter                 100
scf.EigenvalueSolver        Cluster
scf.eigen.lib               elpa1      # elpa1|elpa2; cluster default = elpa1
scf.Mixing.Type             RMM-DIISK
scf.criterion               1.0e-8
```

```bash
cd OpenMX/examples/18-c60-cluster-elpa
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat
NP=8 NT=1 bash run_compare.sh
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `scf.eigen.lib` | Cluster 默认 `elpa1`；Band 默认 `elpa2` | ELPA 随源码内置，无需外部库。手册称 Cluster 中两者速度几乎相同，Band 中 ELPA2 略快 |
| `scf.EigenvalueSolver` | `Cluster` | 分子、团簇（放在足够大的盒子里）；周期体系只用 Γ 点时请用 `Band` + `scf.Kgrid 1 1 1` |
| `level.of.fileout` | 基准测试用 0 | 避免 IO 干扰计时 |

预期：ELPA1 与 ELPA2 的 `Utot` 一致到 SCF 精度；DC(7.0 Å) 与常规对角化相差 mHartree 量级。作为参照，手册表 5 用最小基组得到 C60 常规对角化 −343.90129 Hartree、DC(7.0) −343.90004 Hartree；本例基组不同，只比较差值的量级。对这么小的体系 DC 并不更快，本例的目的是看误差。

---

## 8. GPU 版本：openmx4.0_gpu

官方下载页链接了由 Hiroyuki Kawai（新潟大学，GitHub：dc1394）维护的 GPU 版 [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)（GPL-3.0-or-later）。下载页列出的是 v1.1.1（2026-08-25），GitHub README 已描述 v2.0 的分布式多 GPU 功能，请以你下载版本的 README 为准。根据其 README：

- 加速范围：Band 与 Cluster 的稠密本征值求解（共线和非共线）、O(N) 求解器（DC、DC-LNO、Krylov），以及哈密顿量/重叠矩阵构造、电荷密度与轨道网格、力、电荷混合等；
- 后端：NVIDIA CUDA、AMD HIP，以及不需要 GPU 的 CPU 构建，均用 CMake 选择；
- 该分支还扩展了 `DFT_DATA19` 中镧系元素的 VPS/PAO。

构建（README 中的示例，GPU 架构按你的卡修改：A100 为 `sm_80`、H100 为 `sm_90`）：

```bash
git clone https://github.com/dc1394/openmx4.0_gpu.git
cd openmx4.0_gpu
cmake -S . -B build/cuda -DOPENMX_GPU_BACKEND=CUDA -DOPENMX_GPU_ARCH=sm_80
cmake --build build/cuda
cmake --install build/cuda            # default: puts openmx into work/
cd work && mpirun -np 4 ./openmx -runtest -nt 1
```

运行与输入：

| 设置 | 取值 | 说明 |
|---|---|---|
| `scf.eigen.lib` | `gpusolver`（该分支默认） | cuSOLVER；放不下时在单节点上自动用多 GPU MAGMA（需含 MAGMA 的构建） |
| `scf.eigen.lib` | `gpusolver2` | ELPA GPU 内核 + COSMA 的分布式求解，**只支持 `scf.EigenvalueSolver cluster`** |
| `scf.eigen.lib` | `gpusolver3` | 强制使用 MAGMA 路径 |
| `scf.eigen.lib` | `elpa1` / `elpa2` | 回到 CPU 路径 |
| `scf.Gpu.Num` | 如 4 | 单节点可用 GPU 数上限 |
| `OPENMX_GPU=0`（环境变量） | — | 在有 GPU 的机器上强制走 CPU，便于用同一输入做对比 |
| NVIDIA MPS | 多个 MPI 进程共享一块 GPU 时启用 | `nvidia-cuda-mps-control -d` 启动，`echo quit \| nvidia-cuda-mps-control` 停止 |

> [!CAUTION]
> GPU 版是第三方维护的衍生版本。其 README 说明在 FP64 性能弱的消费级显卡上，大矩阵乘法默认用 INT8 张量核模拟 FP64（GEMMul8），线程归约顺序也会让结果不能逐位复现。正式使用前务必跑 `-runtest`，并用你自己的体系与官方 CPU 版比较总能和力；论文中写明所用的版本和求解器。编译依赖（NVIDIA HPC SDK、ROCm 等）较复杂，README 还提供了 Docker 镜像。

---

## 9. OpenMX API 与开发者接口

### 9.1 OpenMX API（4.0 新增）

4.0 允许外部程序把 OpenMX 当作库调用（手册第 71 章）。核心调用只有一个：`openmx_api(comm, input_file, nthreads)`，在给定的 MPI 通信器上运行一个输入文件。

| 语言 | 生成 | 调用形式 |
|---|---|---|
| C | `make api` → `libopenmx_api.a`；`make example_openmx_api` | `openmx_api(subcomm, "Methane.dat", Nth)`，头文件 `openmx_api.h` |
| Fortran | `make api`；`make example_openmx_api_f08` | `use openmx_api_f08`，`openmx_run(comm, input_file, nth)` |
| Python | `make clean` 后 `make pyapi` → `pyopenmx_api*.so` | `pyopenmx_api.run(subcomm.py2f(), input_file, nth)`，需要 mpi4py |

官方示例把 `MPI_COMM_WORLD` 分成两个子通信器，分别运行 `Methane.dat` 和 `Cdia.dat`，所以至少需要 2 个进程：

```bash
cd ~/openmx4.0/source
make api && make example_openmx_api
cd ../work && cp ../source/example_openmx_api .
mpirun -np 4 ./example_openmx_api
```

示例 [18-api-python](../examples/18-api-python/README.md) 中的 `api_batch.py` 把这一思路推广到任意多个输入：`mpirun -np 8 python3 api_batch.py --nt 1 a.dat b.dat c.dat d.dat` 把 8 个进程分成 4 组，各跑一个输入，适合大量互不相关的小计算一次排队完成。

要点（手册第 71.4–71.5 节）：

- 线程数由 API 参数传入，但实际上限仍受 `OMP_NUM_THREADS` 和调度器分配的核数约束；
- 链接时使用与编译 OpenMX 相同的 MPI、编译器和数学库；Python 版中 mpi4py 必须链接到同一个 MPI，否则可能启动后立即中止；
- `make pyapi` 需要 `-fPIC` 的目标文件，第一次构建前必须 `make clean`；
- 常见构建错误：`python3-config: command not found`、找不到 `Python.h`（缺 Python 开发包），`import pyopenmx_api` 时的共享库错误（Python 版本或运行库路径不匹配）。

### 9.2 `.scfout` 与 `analysis_example`、`elpa_col`/`elpa_noncol`

`HS.fileout on` 时，SCF 结束会写出 `System.Name.scfout`，其中包含哈密顿量、重叠矩阵、密度矩阵、位置算符和动量算符的矩阵元（手册第 70 章）。

```bash
cd ~/openmx4.0/work
make analysis_example
# si64_band.scfout exists only if HS.fileout on was added to Si64_Band_ref.dat
./analysis_example si64_band.scfout > HS.out          # dump H, S, DM ...
make elpa_col                                          # collinear; elpa_noncol for non-collinear
./elpa_col si64_band.scfout 0.0 0.0 0.0                # eigenvalues at k = (0,0,0) in reciprocal-lattice units
```

`analysis_example` 演示了用 `read_scfout.h` 中的 `read_scfout()` 读取数据，是开发后处理程序的起点（polB、jx、kSpin 都基于它）；`elpa_col`/`elpa_noncol` 在读入 `.scfout` 后用 ELPA1/ELPA2 并行对角化，命令行后三个数是 k 点的倒格矢分量。

---

## 10. 生态：`.scfout` 解析与下游工具

| 工具 | 用途 | 需要的 OpenMX 输出 | 备注 |
|---|---|---|---|
| [TB2J](https://github.com/mailhexu/TB2J)（[OpenMX 文档](https://tb2j.readthedocs.io/en/latest/src/openmx.html)）、[TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX) | 磁力定理计算海森堡交换 J、DMI、各向异性交换、磁振子 | `.scfout`（`HS.fileout on`），磁性计算 | 与 [第11章](11-spin-noncollinear.md) 的 `jx` 互补 |
| [HopTB.jl](https://github.com/HopTB/HopTB.jl) | 由 OpenMX 哈密顿量建紧束缚模型，算 Berry 曲率、光学、非线性响应 | `.scfout` | Julia |
| [DeepH-pack](https://github.com/mzjb/DeepH-pack)、[DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3) | 机器学习 DFT 哈密顿量 | OpenMX 计算的哈密顿量与重叠矩阵数据集 | 文档基于 OpenMX 3.9 及其后处理工具 |
| [HamGNN](https://github.com/QuantumLab-ZY/HamGNN) | 图神经网络预测哈密顿量 | OpenMX 数据（附 `openmx_postprocess`/`read_openmx` 编译流程） | 同上 |
| [sisl](https://sisl.readthedocs.io/en/latest/api/io/openmx.html) | 读 OpenMX 输出、紧束缚与输运分析 | 视功能而定 | Python |
| [aiida-plugin-openmx](https://github.com/zeleznyj-htp/aiida-plugin-openmx)、[aiida-openmx](https://github.com/azadoks/aiida-openmx) | AiiDA 工作流：自动提交、重启、数据溯源 | 输入/输出 | 前者较新（含 jx 工作流），后者较旧 |
| [cif2x](https://github.com/issp-center-dev/cif2x) | CIF + 模板批量生成 OpenMX 输入 | — | ISSP，基于 pymatgen |
| `cif2omx`（4.0 自带） | 由 CIF 生成 `.dat` | — | 手册第 8.3 节 |

> [!WARNING]
> 这些工具大多针对 OpenMX 3.9 开发和测试。用 4.0 的 `.scfout` 或输出前，先查看工具的版本说明，并在一个小体系上与 OpenMX 自身结果（例如 `elpa_col` 给出的本征值、`jx` 给出的交换常数）对照，确认读取正确。另外注意 4.0 默认 `scf.Generation.Kpoint GCenter`，与 3.9 的 k 网格不同（[第5章 §5.3](05-scf-kpoints-convergence.md)）。

仓库的 [`omx_parse_out.py`](../scripts/omx_parse_out.py) 可把 `.out`/`.std` 中的 Utot、SCF 历史、Mulliken 电荷、力、优化历史导出为 JSON/CSV，适合接入自己的工作流；用 `python ../scripts/omx_parse_out.py -h` 查看用法。

---

## 参数选择与收敛建议

1. **先小后大**：任何大体系计算之前，先在小超胞上确定基组、截断能（[第4章](04-pao-vps.md)、[第5章](05-scf-kpoints-convergence.md)），再确定 O(N) 的 `orderN.HoppingRanges`（§6.2）。
2. **求解器**：几百原子以内用 Band/Cluster；约 1000 原子且核数充足时常规对角化仍可行；更大或长 MD 用 DC-LNO（团簇 200–300 原子为精度与效率的平衡点），共价体系也可用 DC，金属可试 Krylov。
3. **O(N) 容差**：`orderN.HoppingRanges` 取 ΔE ≤ 1–5 meV/atom 的最小值；比较的能量差越小，容差越严。优化/MD 用 `<LNOs.Num` 固定 LNO 数；变晶胞用 `<orderN.FNAN+SNAN` 固定团簇大小并用 `scf.Ngrid` 固定网格。
4. **并行**：`进程 × 线程 ≤ 物理核数`；`OMP_NUM_THREADS = -nt`；Band 计算的进程数配合 `自旋重数 × N_k`；内存紧张时增大 `-nt`、减少每节点进程。每换一台机器或一类体系，用 `bench_np_nt.sh` 做一次 5–10 步的短基准。
5. **ELPA**：一般保持默认（Cluster `elpa1`、Band `elpa2`）；只有基准显示明显差别时才修改。
6. **IO**：大体系用 `level.of.fileout 0`，IO 慢时 `OutData.bin.flag on`；`HS.fileout` 只在需要 `.scfout` 时打开。
7. **稳定性**：作业脚本固定写 `ulimit -s unlimited` 和 `OMP_STACKSIZE=512M`；按 `.std` 末尾判断成功，按 `.out` 的 SCF 历史判断收敛。
8. **可复现**：记录 OpenMX 版本（含 patch）、编译器与库、`-np`/`-nt`、`scf.eigen.lib`、GPU 版本号；混合并行与 GPU 计算的结果可能有末位差异，比较能量差时用同一配置。

---

## 常见问题与报错

本节先给出排查流程，再按类别列出速查表。每一行只给出要点，详细解释在“详见”一列的章节中。

**排查流程**：

1. `tail -n 40 x.std`：程序是否正常结束？最后一条信息是什么？
2. 在作业日志和 `.std` 中搜索关键报错：`bash check_run.sh x.std`；
3. 正常结束但结果可疑：看 `.out` 中的 SCF 历史（`NormRD`、步数是否等于 `scf.maxIter`）、`Utot`、Mulliken 电荷之和；
4. 崩溃但没有报错信息：多半是内存或栈（段错误、`Killed`、`oom-kill`），查 `sacct`/`dmesg` 和 §5 的内存估计；
5. 只在并行时出错：减少到 `-np 1 -nt 1` 运行同一输入，判断是输入问题还是并行/环境问题；
6. 仍无头绪：用官方 `work/` 中最接近的示例改起，逐项替换为自己的设置，找出引起问题的那一项。

### 编译与安装

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| `multiple definition of 'SpeBasisName'; ... first defined here` | GCC ≥ 10 默认 `-fno-common` | `CC` 加 `-fcommon`，`make clean` 后重编 | [第2章](02-install.md) |
| `Error: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8))` | gfortran ≥ 10 的实参类型检查（`elpa1.f90` 等） | `FC` 加 `-fallow-argument-mismatch` | [第2章](02-install.md) |
| `icc: command not found` / `ifort: command not found`（`mpiicc`、`mpiifort` 找不到底层编译器） | 经典编译器已移除：`icc` 自 oneAPI 2024.0、`ifort` 自 2025.0 | 改用 `mpiicx`/`mpiifx`，加 `-fcommon -Wno-implicit-function-declaration` | [第2章](02-install.md) |
| `/usr/bin/ld: cannot find -lifcore`；去掉后 `undefined reference to 'mpi_abort_'` | GNU 工具链混用了 Intel 库 | 去掉 `-lifcore -liomp5 -lmkl_blacs_intelmpi_lp64`，加 `-lmpi_mpifh -lgfortran`；或用 MKL Link Line Advisor | [第2章](02-install.md) |
| `bandgnu13` 编译失败 | 它需要单独编译 | `gcc bandgnu13.c -lm -o bandgnu13` | [第8章](08-band.md) |
| ARM（aarch64）上 Krylov 相关源文件编译出错 | 默认使用 x86 SSE 指令 | `CC` 中加 `-Dnosse` | [第2章](02-install.md) |
| `runtest.result` 中 `diff Utot` > 1e-7 | 激进优化、库不匹配、混用 MPI | 去掉 fast-math 类选项重编；`ldd openmx` 检查 | [第2章](02-install.md) |

### 运行环境、并行与作业

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| `error while loading shared libraries: libmpifort.so.12: cannot open shared object file` | deb 包运行库不在搜索路径 | `export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH`，直接用命令名调用 | §4.3 |
| `error while loading shared libraries: libmkl_*.so` / `libiomp5.so` | 计算节点上没有加载 oneAPI 环境 | 作业脚本中 `module load` 或 `source setvars.sh` | [第2章](02-install.md) |
| `Segmentation fault` / `signal 11`，尤其大体系或 `-nt > 1` | 栈不足、内存不足 | `ulimit -s unlimited`、`OMP_STACKSIZE=512M`；减少 `-nt` 或每节点进程数；增加节点 | §3.1、§5 |
| 作业被杀，日志有 `Killed`、`oom-kill`、`out of memory` | 超出节点内存 | `memory.usage.fileout on` 估计；减少每节点进程；O(N) 减小团簇 | §5 |
| 多个进程各自从头计算、输出重复 | 用了与编译时不同的 MPI 的 `mpirun` | `which mpirun`，加载与编译一致的 MPI | [第2章](02-install.md) |
| 节点负载远超核数、比少用核还慢 | `OMP_NUM_THREADS` 未设，BLAS 自开线程；或 `-nt` 与申请不符 | `OMP_NUM_THREADS=-nt`；`--cpus-per-task` 与 `-nt` 一致 | §3.1、§4.1 |
| `srun` 启动后每个进程只用 1 个核 | 部分 Slurm 版本（如 22.05–23.02）的 `srun` 不从作业继承 `--cpus-per-task` | `srun --cpus-per-task=$SLURM_CPUS_PER_TASK`（对所有版本都安全），或 `export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK` | §4.1、[公共篇 作业调度](../../Common/tutorials/02-job-scheduler.md) |
| 加核后几乎不再提速 | IO 等串行部分、负载不均；k 点数太少 | `level.of.fileout 0`；Band 计算按 `自旋重数 × N_k` 选进程数 | §2.2、§3.2 |
| 自旋螺旋计算 `Fatal error in internal_Comm_size: Invalid communicator` | 4.0 默认 GCenter 与广义 Bloch 冲突 | `scf.Generation.Kpoint Regular` | [第11章](11-spin-noncollinear.md) |
| `mpirun has detected an attempt to run as root.` | Docker 中以 root 运行 Open MPI | 用普通用户；必要时 `--allow-run-as-root` | [第2章](02-install.md) |
| Python API：`import pyopenmx_api` 失败或启动即中止 | `.so` 与当前 Python 不匹配，或 mpi4py 与 OpenMX 的 MPI 不同 | 用同一套 MPI/编译器重建；检查 `mpi4py.get_config()` | §9.1 |

### 输入文件

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| 找不到 `.pao`/`.vps`，或 `can not open .../PAO/xxx.pao` | `DATA.PATH` 默认 `../DFT_DATA19` 是相对路径 | 写绝对路径，指向含 `PAO/` 和 `VPS/` 的目录（不要写到 `PAO/`） | [第3章](03-input-first-calc.md) |
| `Found an undefined species name ...` 后跟 `Please check your input file` | 物种名不一致；或文件为 UTF-16、带 BOM、全角空格、CRLF | 核对名称；`iconv -f UTF-16 -t ASCII`、`dos2unix` | [第3章](03-input-first-calc.md) |
| 输出没写到预期目录 | 把 `System.CurrrentDirectory`（官方拼写三个 r）写成两个 r | 按官方拼写 | [第1章](01-overview.md) |
| `ERROR key=...` 后跟 `errors in the inputfile` | 选项值非法（如 `scf.XcType PBE`） | 按手册写合法值 | [第3章](03-input-first-calc.md) |
| `Format error of the sequential number ...` / 块内读错 | 坐标块序号不连续、块内写了注释或空行 | 块内只写数据行 | [第3章](03-input-first-calc.md)、[第11章](11-spin-noncollinear.md) |
| `Species.Number may be wrong.` / `Atoms.Number may be wrong.` | 缺关键词或值与块中行数不符 | 补上并核对行数 | [第3章](03-input-first-calc.md) |
| 修改了关键词却没有生效 | 拼写错误被当作未知关键词忽略；或文件中有重复关键词 | `grep -in "关键词" x.dat`；在 `.out` 的输入回显中确认 | [第3章](03-input-first-calc.md) |
| 初始电子数之和与价电子数不符、磁矩异常 | 坐标行第 6、7 列之和不等于 VPS 价电子数 | 按 VPS 核对（如 `Fe_PBE19H`=16、`Fe_PBE19S`=14） | [第11章](11-spin-noncollinear.md) |
| 半芯态异常、能量离谱 | PAO 与 VPS 来自不同数据库版本 | 统一使用 DFT_DATA19 | [第4章](04-pao-vps.md) |
| `MO.fileout=ON is not supported in case of scf.EigenvalueSolver=DC` | O(N) 求解器不输出波函数 | 用 §6.5 的一次对角化 | §6.5、[第10章](10-charge-mo-cube.md) |

### SCF

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| SCF 步数等于 `scf.maxIter`，`NormRD` 振荡 | 金属、磁性、非共线、大体系 | `RMM-DIISK`（+U 用 `RMM-DIISH`），`scf.Mixing.History 30`，`StartPulay 10–30`，减小 `scf.Max.Mixing.Weight`，适当提高电子温度；运行中用 `<System.Name>_SCF_keywords` 调整，再 `scf.restart on` | [第5章](05-scf-kpoints-convergence.md) |
| 4.0 与 3.9 的结果对不上 | 4.0 默认 GCenter k 点 | 显式 `scf.Generation.Kpoint Regular` | [第5章](05-scf-kpoints-convergence.md) |
| 致密金属用大基组时出现伪本征值、能量异常 | 基组过完备 | 换较小的优化基组，或 `scf.ProExpn.VNA off` 并提高 `scf.energycutoff` | [第4章](04-pao-vps.md) |
| E–V 曲线有锯齿 | 实空间网格数随晶格常数跳变；O(N) 时团簇原子数跳变 | `scf.Ngrid` 固定网格；O(N) 加 `<orderN.FNAN+SNAN` | [第5章](05-scf-kpoints-convergence.md)、§6.6 |
| O(N) 计算 SCF 难收敛或能量明显偏离 | 团簇太小（金属尤甚） | 增大 `orderN.HoppingRanges`；金属改 DC-LNO/Krylov；Krylov 增大 `orderN.KrylovH.order`，简单金属试 `orderN.Exact.Inverse.S off`、`orderN.Expand.Core off` | §6.2、§6.3 |
| DC-LNO 的 MD/优化中能量不连续 | LNO 个数随结构变化 | 用 `<LNOs.Num` 固定 | §6.3 |
| DFT+U 不生效 | `scf.SpinPolarization off`，或 `Hubbard.U.values` 项数与 PAO 不符 | 用 `on`/`nc`；按 PAO 列出全部轨道 | [第13章](13-dft-u.md) |
| SOC 无效或报 `Spin-orbit coupling is not supported for collinear DFT calculations.` | 共线计算加了 SOC | `scf.SpinPolarization nc` + `scf.SpinOrbit.Coupling on` | [第12章](12-soc-topology.md) |

### 结构优化与 MD

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| 优化长时间不收敛或停在鞍点 | 弱相互作用体系、DIIS 过早 | 增大 `MD.Opt.StartDIIS`，先用最速下降；截断能提高到 300–400 Ry；必要时 `MD.Opt.criterion` 放宽到 5e-4 | [第6章](06-optimization-eos-neb.md) |
| 续算后立即结束 | `.dat#` 中 `MD.Current.Iter` 已到 `MD.maxIter` | 增大 `MD.maxIter` | [第7章](07-md.md) |
| 晶胞优化报错或不支持 | 晶胞优化只支持共线计算，且不能与 DFT-D 同时用 | 改为共线；vdW 体系手动扫描晶格常数 | [第6章](06-optimization-eos-neb.md) |
| `Non-collinear calculation is not compatible with stress calculation.` | NPT/晶胞优化需要应力 | 只用共线计算 | [第7章](07-md.md) |
| NVE 守恒量漂移 | 步长太大、SCF 不够收敛、网格太粗 | 减小步长、收紧 `scf.criterion`、提高截断能 | [第7章](07-md.md) |
| `.ene`/`.md` 中出现重复步号 | 同名重跑时以追加方式写文件 | 删除旧文件或换 `System.Name` | [第7章](07-md.md) |

### 后处理

| 报错 / 现象 | 原因 | 解决 | 详见 |
|---|---|---|---|
| 没有 cube 文件 | `level.of.fileout 0` | 用 `scf.restart pp` + `level.of.fileout 1` 补算 | §5、[第10章](10-charge-mo-cube.md) |
| 能带标签对不上（fcc/bcc/六方） | `Band.kpath` 以倒格矢为单位，用了原胞 | 加 `<Band.KPath.UnitCell`，或用 OpenMX Viewer/cif2omx 生成路径 | [第8章](08-band.md) |
| 自己画的能带偏移/单位不对 | `.Band` 是 Hartree 且未减化学势 | (E − μ) × 27.2114 eV，或读 `BANDDAT` | [第8章](08-band.md) |
| 0 eV 处积分 DOS 不等于电子数 | `DosMain` 从 `Dos.Erange` 下限积分、高斯展宽 | 下限放到所有价带以下，或比较两能量处积分之差 | [第9章](09-dos-pdos-cohp-optics.md) |
| nc+SOC+CDDF 段错误或 `Parameter 10 was incorrect on entry to ZHEEVX` | 4.0 已知问题 | 参考 [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix) 或关注官方补丁 | [第9章](09-dos-pdos-cohp-optics.md) |
| 费米面文件数值差一个固定倍数 | patch 4.0.1 起 bxsf/frmsf 单位改为 eV/Å | 更新后处理脚本的单位换算 | [第1章](01-overview.md) |
| 非共线 + 仅 Γ 点时 HOMO 序号错误 | 4.0 原版 `Mulliken_Charge.c` 的问题 | 打 patch 4.0.1 | [第2章](02-install.md) |
| ESM 报 `<ESM:Warning> The a-axis ...` | 表面法向不是第一条晶格矢量 a，或 b、c 与 a 不正交 | 把法向矢量放到 a（晶格矢量循环置换为 c、a、b），b、c 与 a 正交；只改 `ESM.direction` 而不重排矢量仍会报此警告 | [第16章](16-surface-esm-workfunction.md) |
| NEGF 理想体系透射不是整数 | 电极 k 点不足、几何不满足重叠条件 | 加密 Step 1 沿 a 的 k 点；检查 `<TRAN_Calc_GridBound>` | [第17章](17-negf-transport.md) |
| 下游工具读 `.scfout` 报错或结果离谱 | 没开 `HS.fileout`；工具只支持 3.9 的格式 | `HS.fileout on`；先在小体系上对照验证 | §10 |

---

## 练习

1. 在你的集群上，用 `bench_np_nt.sh` 对 C60（`scf.maxIter 10`）扫描单节点上至少 4 种“进程 × 线程”组合，画出总时间柱状图，并确认所有组合的 `Utot` 一致。
2. 运行 `scan_hopping.sh`，画出三种 O(N) 方法的 |ΔE|（meV/atom，对数纵轴）随 R 的变化，再画出单步耗时随 R 的变化，给出你认为最合适的 R。
3. 把 Si64 的晶格常数改为 5.30、5.43、5.56 Å，分别用 DC（固定 `orderN.HoppingRanges`）和 DC + `<orderN.FNAN+SNAN`（数值取 5.43 Å 结构的输出）计算，比较两条 E(a) 曲线的平滑程度。
4. 用 `make_si_supercell.py` 生成 Si216、Si512、Si1000 附近（如 5×5×5 = 1000）的 DC-LNO 输入，在相同核数下记录单步 SCF 耗时，验证是否近似线性标度；再用常规 Band（Γ 点）算 Si216 和 Si512，比较交叉点。
5. 对 Si512 完成“DC-LNO 自洽 + 一次 Band 对角化”，用 `DosMain` 得到 DOS，与 Si64 4×4×4 常规计算的 DOS 比较带隙附近的形状。
6. 打开 `memory.usage.fileout on`，比较 `-np 32 -nt 1` 与 `-np 8 -nt 4` 时单进程和整个节点的内存占用。

---

## 参考资料

**官方文档（OpenMX 4.0 手册）**

- 第 6、7 章 自动测试与大体系测试（`-runtest`、`-runtestL`）：<https://www.openmx-square.org/openmx_man4.0/s6_auto-test.html>，<https://www.openmx-square.org/openmx_man4.0/s7_auto-test-large.html>
- 第 8.2 节 关键词：<https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html>
- 第 27 章 O(N) 方法：<https://www.openmx-square.org/openmx_man4.0/s27_on-methods.html>
  - 27.1 DC：<https://www.openmx-square.org/openmx_man4.0/s27_1_dc.html>
  - 27.2 DC-LNO：<https://www.openmx-square.org/openmx_man4.0/s27_2_dc-lno.html>
  - 27.3 Krylov：<https://www.openmx-square.org/openmx_man4.0/s27_3_krylov.html>
  - 27.4 FNAN+SNAN：<https://www.openmx-square.org/openmx_man4.0/s27_4_fnan-snan.html>
- 第 28 章 并行：<https://www.openmx-square.org/openmx_man4.0/s28_parallel.html>（28.1 MPI：<https://www.openmx-square.org/openmx_man4.0/s28_1_mpi-parallel.html>；28.2 混合并行：<https://www.openmx-square.org/openmx_man4.0/s28_2_mpi-openmp-hybrid.html>）
- 第 29 章 大体系计算（`-runtestL2`、O(N)+对角化）：<https://www.openmx-square.org/openmx_man4.0/s29_large-scale.html>
- 第 48 章 数值精确的低阶标度方法：<https://www.openmx-square.org/openmx_man4.0/s48_low-order-diag.html>
- 第 70 章 开发者接口（`.scfout`、`elpa_col`）：<https://www.openmx-square.org/openmx_man4.0/s70_dev-interface.html>
- 第 71 章 OpenMX API：<https://www.openmx-square.org/openmx_man4.0/s71_openmx-api.html>
- 第 75 章 内存使用分析：<https://www.openmx-square.org/openmx_man4.0/s75_memory-usage.html>
- 第 76 章 二进制输出：<https://www.openmx-square.org/openmx_man4.0/s76_binary-output.html>
- 第 79 章 已知问题：<https://www.openmx-square.org/openmx_man4.0/s79_known-issues.html>
- 第 4.11 节 deb 包安装：<https://www.openmx-square.org/openmx_man4.0/s4_11_deb-install.html>
- 下载页（4.0、patch 4.0.1、GPU 版、deb、Docker）：<https://www.openmx-square.org/download.html>
- What's new：<https://www.openmx-square.org/whatsnew/whatsnew.html>
- OpenMX 论坛：<https://www.openmx-square.org/forum/patio.cgi>
- Ozaki 讲座《Practical guide to OpenMX calculations》（2023）：<https://t-ozaki.issp.u-tokyo.ac.jp/lectures-2023/OpenMX_Practical_Guide.pdf>
- Intel oneMKL Link Line Advisor：<https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html>
- Slurm `srun` 文档（`--cpus-per-task`、`SRUN_CPUS_PER_TASK`）：<https://slurm.schedmd.com/srun.html>

**相关 GitHub 仓库（只链接，不复制）**

- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：4.0 GPU 版（CUDA/HIP/CPU，CMake）；3.9.9 的 GPU 版见 [dc1394/openmx3.9.9_gpu](https://github.com/dc1394/openmx3.9.9_gpu)
- [mfukudaQED/openmx_installer](https://github.com/mfukudaQED/openmx_installer)：下载、打补丁、改 makefile、编译的一键脚本
- [atsushi-m-ito/openmx-patch-oneapi](https://github.com/atsushi-m-ito/openmx-patch-oneapi)：用 icx/ifx 编译 OpenMX 的补丁说明
- [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix)：4.0 nc+SOC+CDDF 问题的分析与补丁
- [mailhexu/TB2J](https://github.com/mailhexu/TB2J)、[mailhexu/TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX)：交换参数与磁振子
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：紧束缚后处理
- [mzjb/DeepH-pack](https://github.com/mzjb/DeepH-pack)、[Xiaoxun-Gong/DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3)、[QuantumLab-ZY/HamGNN](https://github.com/QuantumLab-ZY/HamGNN)：机器学习哈密顿量
- [issp-center-dev/cif2x](https://github.com/issp-center-dev/cif2x)、[zeleznyj-htp/aiida-plugin-openmx](https://github.com/zeleznyj-htp/aiida-plugin-openmx)、[azadoks/aiida-openmx](https://github.com/azadoks/aiida-openmx)：输入生成与工作流
- [hirtatsu/Calculations](https://github.com/hirtatsu/Calculations)：日文实操笔记，含超算上的编译与 PBS 作业脚本
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站

**关键文献**

- T. Ozaki, *Phys. Rev. B* **74**, 245101 (2006)，O(N) Krylov 子空间方法：<https://doi.org/10.1103/PhysRevB.74.245101>
- T. Ozaki, *Phys. Rev. B* **82**, 075131 (2010)，数值精确的低阶标度对角化方法：<https://doi.org/10.1103/PhysRevB.82.075131>
- T. Ozaki, M. Fukuda, G. Jiang, *Phys. Rev. B* **98**, 245137 (2018)，DC-LNO 方法：<https://doi.org/10.1103/PhysRevB.98.245137>
- W. Yang, *Phys. Rev. Lett.* **66**, 1438 (1991)，分治法（divide and conquer）：<https://doi.org/10.1103/PhysRevLett.66.1438>
- T. V. T. Duy, T. Ozaki, *Comput. Phys. Commun.* **185**, 777 (2014)，OpenMX 的三维区域分解并行：<https://doi.org/10.1016/j.cpc.2013.11.008>
- A. Marek *et al.*, *J. Phys.: Condens. Matter* **26**, 213201 (2014)，ELPA 本征值求解库：<https://doi.org/10.1088/0953-8984/26/21/213201>
- H. Li *et al.*, *Nat. Comput. Sci.* **2**, 367 (2022)，DeepH：<https://doi.org/10.1038/s43588-022-00265-6>

---

[← 上一章](17-negf-transport.md) | [返回目录](../README.md)
