# 第1章 环境准备、VASP 编译安装（CPU/GPU）与作业提交（Installing & Compiling VASP 6 on CPU/GPU and Submitting Jobs）

本章从零开始讲解 VASP 6.x 的编译安装：如何选择 `arch/makefile.include.*` 模板、如何开启 HDF5 / Wannier90 / DFT-D4 / libbeef 等可选功能、如何用 NVHPC 编译 OpenACC GPU 版、如何打 VTST 与 VASPsol 补丁，以及如何在 Slurm/PBS 集群上提交 vasp_std / vasp_gam / vasp_ncl 作业并设置 KPAR、NCORE、NSIM 等并行参数。读完后你应能独立完成“VASP compile + install + job submission”的全流程，并用官方测试集和 Si 小算例验证安装。

> [!IMPORTANT]
> VASP 是商业授权软件。**源码、POTCAR 赝势库、随源码分发的 `arch/makefile.include.*` 模板都受许可证保护**，不能放进公开仓库，也不能上传到消费级 AI 工具（VASP 官方 2026-06-16 发布了 AI 工具使用规范）。本章只描述“需要改哪些变量、为什么改”，不复制官方模板原文；示例目录里也只有我们自己写的“改动清单”。另外注意：`vaspout.h5` 内含 POTCAR 全文，分享计算结果时同样不能随意公开。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| 理解软件栈 | 编译器、MPI、BLAS/LAPACK/ScaLAPACK、FFTW、HDF5 各自的作用 |
| 选对模板 | oneAPI（ifx+MKL）、GNU（FOSS）、AOCC/AOCL、NVHPC（OpenACC GPU）、Cray/oneAPI（OpenMP offload GPU） |
| 开启可选功能 | `-DVASP_HDF5`、`-DVASP2WANNIER90`、`-Dlibbeef`、`-DDFTD4`、`-DSDFTD3`、`-DLIBMBD`、`-DPLUGINS`、`-Dlibvaspml` |
| 编译三个程序 | `make DEPS=1 -j N all` 得到 vasp_std / vasp_gam / vasp_ncl |
| 打补丁 | VTST（CI-NEB/Dimer/IOPT 优化器）、VASPsol（隐式溶剂） |
| 验证与运行 | `make test`、Si 小算例、Slurm/PBS 作业脚本（CPU 与 GPU） |
| 并行调优 | KPAR、NCORE、NSIM 的选择规则与基准测试方法 |

适用读者：已经拿到 VASP 许可证（或所在课题组有许可证）、需要在个人工作站、WSL2 或 HPC 集群上自行编译的用户；只想“用现成模块”的读者可以直接跳到第 12 节（运行与作业提交）和第 13 节（并行参数）。

前置知识：Linux 命令行、环境变量、`module` 系统和作业调度的基础，请先阅读公共章节 [第1章 计算环境搭建](../../Common/tutorials/01-linux-hpc.md) 与 [第2章 HPC 作业调度](../../Common/tutorials/02-job-scheduler.md)。

---

## 2. 原理速览：VASP 是怎样“跑起来”的

### 2.1 软件栈分层

VASP 是基于平面波基组和 PAW 方法的 Fortran 程序。一次计算的大部分时间花在三类运算上：

1. **三维 FFT**：在实空间和倒空间之间变换波函数与电荷密度，需要 FFTW 或 MKL 提供的 FFTW 接口；
2. **稠密线性代数**：子空间对角化、正交化，需要 BLAS / LAPACK，以及分布式的 ScaLAPACK；
3. **进程间通信**：能带、k 点、FFT 网格的并行分配，需要 MPI；GPU 版还需要 NCCL（NVIDIA）或 RCCL（AMD）。

因此一个可用的 VASP 构建 = **Fortran/C/C++ 编译器 + MPI + BLAS/LAPACK/ScaLAPACK + FFTW**（必需）+ **HDF5**（强烈推荐，py4vasp 依赖它）+ 各种可选库。并行总进程数满足（官方 Optimizing the parallelization 页面给出的关系）：

$$
N_{\text{ranks}} = N_{\text{band groups}} \times \text{NCORE} \times \text{KPAR} \times \text{IMAGES} \times \cdots
$$

### 2.2 三个可执行程序

同一份源码通过不同的预处理选项编译出三个程序（旧版本靠 `-DNGZhalf`、`-DwNGZhalf` 等手动控制，6.x 中这些选项已废弃，改为 `make std/gam/ncl` 目标自动处理）：

| 程序 | 波函数 | 适用场景 | 不能用于 |
|---|---|---|---|
| `vasp_std` | 复数，共线自旋 | 绝大多数体相、表面、磁性（ISPIN=2）计算 | SOC、非共线磁性 |
| `vasp_gam` | 实数（只含 Γ 点），利用 $\psi(-\mathbf{r})$ 对称性 | 大超胞、分子、AIMD 中只用 Γ 点的计算；内存约减半、速度明显提升 | 多 k 点、KPOINTS_OPT、LOBSTER |
| `vasp_ncl` | 二分量旋量 | 自旋轨道耦合（LSORBIT）、非共线磁性（LNONCOLLINEAR） | —（可以算共线体系，但更慢） |

> [!TIP]
> 判断用哪个程序最简单的办法：KPOINTS 只有 Γ 点（`Gamma 1 1 1`）→ `vasp_gam`；INCAR 里有 `LSORBIT = .TRUE.` 或 `LNONCOLLINEAR = .TRUE.` → `vasp_ncl`；其余 → `vasp_std`。用错程序时 VASP 通常会直接报错退出（例如用 vasp_std 跑 SOC）。

### 2.3 版本选择（截至 2026-11）

| 版本 | 发布日期 | 与本教程相关的主要变化 | 建议 |
|---|---|---|---|
| **6.6.1** | 2026-07-20 | 修复 6.6.0 中自旋极化（ISPIN=2）meta-GGA 受力错误（known issue 94）、OpenACC 下 LCALCEPS 的 Born 有效电荷错误（issue 93）、KPOINTS_OPT + ISPIN=2 的 PDOS 错误（issue 69）、2D 库仑截断噪声（issue 72）、GRACE 兼容问题 | **首选** |
| 6.6.0 | 2026-03-10 | Intel/AMD GPU 的 OpenMP offload（beta）、simple-DFT-D3（IVDW=15）、适配 DFT-D4 4.x 新 API、CHECKPOINT_FD、KSPACING_OPT、ML_MODE=delta、GRACE（实验性）、新增 CMake 构建 | 若只能用 6.6.0，必须打官方 `pot.F.patch`（`patch -p1 < pot.F.patch`）后重编 |
| 6.5.1 | 2025-03-11 | ML_IERR 更名 ML_ESTBLOCK；ISEARCH=1 改进 ALGO=A 线搜索；LNONCOLLINEAR 优先于 ISPIN=2 | 可用 |
| 6.5.0 | 2024-12-17 | 电声耦合、Python 插件（-DPLUGINS）、部分电荷写入 vaspout.h5；OpenACC 预编译选项由 `-D_OPENACC` 改为 `-DACC_OFFLOAD`（另需 `-DNVCUDA`） | 可用 |
| 6.4.3 | 2024-03-19 | BANDGAP、WRT_POTENTIAL、libMBD、ALGO=Dav；**首次提供 ifx（oneAPI LLVM 编译器）的模板** | 老集群常见 |
| 6.3.x | 2022 | KPOINTS_OPT、MLFF；旧 CUDA-C GPU 版被彻底移除 | 建议升级 |

> [!WARNING]
> 6.6.0 的自旋极化 SCAN/r2SCAN 等 meta-GGA **受力是错的**（能量正确），会直接影响结构优化和 MD。请升级到 6.6.1；许可证不覆盖 6.6.1 的用户，按官方说明给 6.6.0 打 `pot.F.patch`。

---

## 3. 准备工作

### 3.1 了解你的硬件

编译选项（如 `-xHOST`）和并行参数（NCORE、KPAR）都和硬件有关，先查清楚：

```bash
lscpu | grep -E "Model name|Socket|Core|Thread|NUMA"   # CPU 型号、插槽、核数、NUMA 域
numactl --hardware                                      # 每个 NUMA 域包含哪些核
free -g                                                 # 内存
nvidia-smi                                              # NVIDIA GPU 型号、显存、驱动/CUDA 版本
module avail 2>&1 | grep -i -E "intel|oneapi|gcc|openmpi|nvhpc|hdf5|fftw"
```

记录下：每节点物理核数、NUMA 域大小（例如每 8 核一个 NUMA 域）、是否超线程、GPU 型号与数量。

### 3.2 选择工具链

| 场景 | 推荐工具链 | 对应模板（arch/） | 备注 |
|---|---|---|---|
| Intel CPU 集群 | oneAPI：ifx + icx + Intel MPI + MKL | `makefile.include.oneapi`（纯 MPI）/ `oneapi_omp`（MPI+OpenMP） | 6.4.3 起提供 ifx 模板；oneAPI 2025 已移除 ifort |
| AMD EPYC | GCC + OpenMPI + AOCL，或 AOCC + AOCL | `gnu_ompi_aocl(_omp)`、`aocc_ompi_aocl(_omp)` | 也可以用 oneAPI，但要改 `-xHOST` |
| 个人电脑 / WSL2 | GCC + OpenMPI + OpenBLAS + ScaLAPACK + FFTW（系统包管理器安装） | `gnu` / `gnu_omp` | 官方有 Debian/Ubuntu/Fedora/Rocky/macOS 的逐步说明 |
| NVIDIA GPU | NVHPC（≥ 21.2）+ 自带 OpenMPI/NCCL + FFTW 或 MKL | `nvhpc_acc`、`nvhpc_omp_acc`、`nvhpc_ompi_mkl_omp_acc` | 避开 NVHPC 22.1/22.2（与 OpenMP 同用有严重 bug） |
| AMD GPU（6.6 beta） | HPE Cray CCE ≥ 19.0 + ROCm | `cray_omp_off` | 只在 HPE 集群上可用 |
| Intel GPU（6.6 beta） | ifx ≥ 2025.3.0 + oneMKL + Intel MPI | `oneapi_omp_off` | 旧版 ifx 能编译但会慢或出错，2025.3.0 是硬性要求 |

> [!NOTE]
> VASP 不能在 Windows 上原生运行。Windows 用户请先安装 WSL2 + Ubuntu，再按 Linux 流程编译（官方 Installing VASP.6.X.X 页面也是这样建议的）。

### 3.3 获取源码与目录结构

许可证持有者从 VASP Portal 下载 `vasp.6.6.1.tgz`，并用随附的校验和文件验证：

```bash
mkdir -p ~/software && cd ~/software
sha256sum vasp.6.6.1.tgz              # 与 Portal 提供的校验值比对
tar -xzf vasp.6.6.1.tgz
cd vasp.6.6.1
ls                                     # arch  bin  build  src  testsuite  tools  makefile ...
ls arch | head -40                     # 所有 makefile.include.* 模板
```

| 目录 | 内容 |
|---|---|
| `arch/` | 各工具链的 `makefile.include.*` 模板 |
| `src/` | Fortran 源码、低层 makefile、隐藏文件 `.objects`（对象列表，VTST 补丁要改它） |
| `src/lib`、`src/parser`、`src/fftlib`、`src/vaspml` | 内部库；`vaspml` 为机器学习力场 C++ 库（≥ 6.5.0） |
| `build/` | std/gam/ncl 三个版本分别在子目录中编译 |
| `bin/` | 编译生成的可执行文件 |
| `testsuite/` | 官方正确性测试集 |
| `tools/` | HDF5 相关的 Python 小工具 |

### 3.4 安装依赖库

**方案 A：Ubuntu 24.04 / WSL2（最省事，性能非最优）**

```bash
sudo apt-get update
sudo apt install rsync make build-essential g++ gfortran \
     libopenblas-dev libopenmpi-dev libscalapack-openmpi-dev \
     libfftw3-dev libhdf5-openmpi-dev
```

这套组合与官方 Personal computer installation 页面一致，对应模板 `makefile.include.gnu_omp`。系统 HDF5 的库和头文件位于 `/usr/lib/x86_64-linux-gnu/hdf5/openmpi/` 与 `/usr/include/hdf5/openmpi/`，稍后在 makefile.include 中直接写这两个路径即可。

**方案 B：Intel oneAPI（集群常用）**

```bash
# 集群上通常是 module load；自装 oneAPI 则 source setvars.sh
source /opt/intel/oneapi/setvars.sh
ifx --version            # Fortran（LLVM 后端）
icx --version            # C
mpiifx -v 2>/dev/null || mpiifort -fc=ifx -v   # Intel MPI 的 Fortran 包装器
echo $MKLROOT            # MKL 根目录，模板会用到
```

**HDF5 源码编译（任意工具链通用）**：HDF5 的 Fortran 模块文件（`hdf5.mod`）与编译器绑定，**必须用编译 VASP 的同一个 Fortran 编译器编译 HDF5**，否则会出现 “Cannot read module file” 一类错误。

```bash
# 以 oneAPI 为例；GNU 工具链把 CC/FC 换成 gcc/gfortran，NVHPC 换成 nvc/nvfortran
cd ~/software && tar -xzf hdf5-1.14.5.tar.gz && cd hdf5-1.14.5
CC=icx FC=ifx ./configure --prefix=$HOME/opt/hdf5-1.14.5-oneapi --enable-fortran
make -j 16 && make install
ls $HOME/opt/hdf5-1.14.5-oneapi/lib | grep hdf5_fortran   # 应看到 libhdf5_fortran.*
```

> [!TIP]
> VASP 只需要串行 HDF5 的 Fortran 接口（官方 Ubuntu 示例用的是 OpenMPI 版的系统包，自己编译时串行版即可）。如果链接的是共享库，运行时要把 `$HDF5_ROOT/lib` 加入 `LD_LIBRARY_PATH`，或者在 `LLIBS` 里加 `-Wl,-rpath,$(HDF5_ROOT)/lib` 把路径写进可执行文件。

---

## 4. 选择 makefile.include 模板

官方模板按编译器分组（以 6.6.x 的 arch/ 为准）：

| 模板 | 并行方式 | 数学库 | 说明 |
|---|---|---|---|
| `makefile.include.oneapi` | MPI | MKL | **ifx**，Intel CPU 首选 |
| `makefile.include.oneapi_omp` | MPI + OpenMP | MKL | 每节点核数很多（≥ 64）时可考虑 |
| `makefile.include.oneapi_omp_off` | MPI + OpenMP offload | oneMKL | Intel 数据中心 GPU（6.6 beta） |
| `makefile.include.intel`、`intel_omp`、`intel_ompi_mkl_omp` | MPI(+OpenMP) | MKL | 旧的 ifort（classic），新 oneAPI 已无 ifort |
| `makefile.include.gnu`、`gnu_omp` | MPI(+OpenMP) | OpenBLAS + ScaLAPACK + FFTW | 自由开源工具链 |
| `makefile.include.gnu_ompi_mkl_omp` | OpenMPI + OpenMP | MKL | GCC 编译器 + MKL |
| `makefile.include.gnu_ompi_aocl(_omp)` | OpenMPI(+OpenMP) | AOCL | AMD CPU |
| `makefile.include.aocc_ompi_aocl(_omp)` | OpenMPI(+OpenMP) | AOCL | AMD 编译器 |
| `makefile.include.nvhpc`、`nvhpc_omp`、`nvhpc_ompi_mkl_omp` | CPU 版 | NVHPC 自带 / MKL | NVHPC 编 CPU 版 |
| `makefile.include.nvhpc_acc`、`nvhpc_omp_acc`、`nvhpc_ompi_mkl_omp_acc` | OpenACC GPU | NVHPC 自带 / MKL | NVIDIA GPU |
| `makefile.include.cray`、`cray_omp`、`cray_omp_off` | Cray CPU / AMD GPU offload | Cray LibSci / ROCm | HPE Cray 机器 |
| `nec_aurora`、`fujitsu_a64fx(_omp)`、`amdflang(_omp)` | — | — | 特殊平台 |

> [!IMPORTANT]
> **永远使用与源码同一版本发布的模板**。例如 6.5.0 把 OpenACC 的预编译选项从 `-D_OPENACC` 改成了 `-DACC_OFFLOAD`（另需 `-DNVCUDA`），拿 6.4 的 makefile.include 去编 6.6.1 会出现各种奇怪的编译或运行错误。

选择流程可以概括为：

1. 有 NVIDIA GPU 并打算用 GPU 跑 → `nvhpc_*acc`（宿主机是 Intel CPU 时优先 `nvhpc_ompi_mkl_omp_acc`）；
2. 纯 CPU、集群有 oneAPI → `oneapi`（或 `oneapi_omp`）；
3. AMD CPU → `gnu_ompi_aocl` 或 `aocc_ompi_aocl`；
4. 个人电脑 → `gnu_omp`（配合系统包）。

> [!NOTE]
> “纯 MPI”还是“MPI+OpenMP”？官方说明：只用 1 个 OpenMP 线程时，MPI+OpenMP 版本与纯 MPI 版本几乎没有性能差别，所以编译 `_omp` 版本并没有坏处；但用 OpenMP 时如果不用 MKL 的 FFT，**应同时启用 fftlib**（模板中有现成注释块），否则多线程 FFT 规划开销会很大。

---

## 5. CPU 版编译实战（oneAPI：ifx + MKL + Intel MPI）

### 5.1 复制模板

```bash
cd ~/software/vasp.6.6.1
cp arch/makefile.include.oneapi makefile.include
```

### 5.2 需要改动的地方

打开 `makefile.include`，文件后半部分有一行 `## Customize as of this point!`，大部分修改都在这一行之后。按下表逐项检查（完整改动清单见示例文件 [makefile.include.oneapi.changes](../examples/01-setup-compile/makefile.include.oneapi.changes)）：

| 位置 / 变量 | 默认含义 | 需要怎么改 | 原因 |
|---|---|---|---|
| `FC`、`FCL`（文件前部） | 用 Intel MPI 的 Fortran 包装器调用 ifx | 确认调用的是 ifx：`mpiifort -fc=ifx` 或 `mpiifx`（Intel MPI ≥ 2021.11 提供） | oneAPI 2025 已移除 ifort，仍调用 ifort 会报 `ifort: command not found` |
| `CC_LIB`、`CXX_PARS` | icx / icpx | 一般不改 | lib 与 parser 两个内部库用 C/C++ 编译 |
| `VASP_TARGET_CPU` | `-xHOST`（针对编译机 CPU 优化） | 在计算节点上编译且为 Intel CPU：保持；登录节点与计算节点代际不同或为 AMD CPU：改为显式目标，如 `-march=core-avx2`（AVX2，AMD Zen 也可用）或 `-xCORE-AVX512` | `-xHOST` 生成的指令集可能在别的 CPU 上不可用或很慢 |
| `MKLROOT` | `?=` 默认值 | 通常不用改，`setvars.sh` / `module load` 已经导出 | 用 `echo $MKLROOT` 检查 |
| HDF5 区块 | 默认注释掉 | 取消注释（或在文件末尾追加），设置 `HDF5_ROOT` | vaspout.h5、py4vasp、电声耦合等功能需要 |
| 其他可选库 | 默认关闭 | 按需追加，见第 6 节 | — |

一个只做“最小必要改动 + HDF5”的追加片段如下（**这是我们自己写的追加内容，不是官方模板**；路径请换成你的）：

```makefile
# ---- appended to the end of makefile.include (after the Customize marker) ----
# Target CPU: keep -xHOST only if you compile on the Intel node that runs VASP
VASP_TARGET_CPU = -xHOST

# HDF5 (strongly recommended; needed by py4vasp and several features)
CPP_OPTIONS += -DVASP_HDF5
HDF5_ROOT   ?= $(HOME)/opt/hdf5-1.14.5-oneapi
LLIBS       += -L$(HDF5_ROOT)/lib -lhdf5_fortran -Wl,-rpath,$(HDF5_ROOT)/lib
INCS        += -I$(HDF5_ROOT)/include
```

> [!TIP]
> 模板中的路径变量大多用 `?=` 赋值，意思是“如果环境中还没有这个变量才用默认值”。所以也可以不改文件，直接 `export HDF5_ROOT=...` 后再 `make`；不过把路径写进 makefile.include 更利于复现，建议把最终的 makefile.include 保存在**私有**位置（不要放进公开仓库）。

### 5.3 编译

```bash
cd ~/software/vasp.6.6.1
make veryclean                      # 换了 makefile.include 后先清理旧的 build/
make DEPS=1 -j 16 all 2>&1 | tee make.log
# 只编部分版本：make DEPS=1 -j 16 std gam ncl
```

- `DEPS=1` 让 make 自动分析 Fortran 模块依赖，配合 `-j N` 并行编译；
- `all` 等价于 `std gam ncl` 三个目标；
- 在 16 核机器上完整编译三个版本一般需要十几分钟到半小时（取决于编译器和可选库）。

### 5.4 检查编译产物

```bash
ls -lh bin/                                   # vasp_std  vasp_gam  vasp_ncl
grep -i -E " error |error:" make.log | head   # 应为空
ldd bin/vasp_std | grep -E "mkl|mpi|hdf5"     # 确认链接到了预期的库（静态链接的库不会出现）
ldd bin/vasp_std | grep "not found"           # 应为空，否则运行时会找不到共享库
echo 'export PATH=$HOME/software/vasp.6.6.1/bin:$PATH' >> ~/.bashrc
```

### 5.5 GNU / AOCL 工具链的不同之处

以 Ubuntu 系统包为例（`cp arch/makefile.include.gnu_omp makefile.include`），在 `## Customize as of this point!` 之后：

| 变量 | 改法 | 说明 |
|---|---|---|
| `OPENBLAS_ROOT` / `BLASPACK` | 注释掉 `OPENBLAS_ROOT`，令 `BLASPACK = -lopenblas` | 系统包已在默认搜索路径 |
| `SCALAPACK_ROOT` / `SCALAPACK` | 注释掉 `SCALAPACK_ROOT`，令 `SCALAPACK = -lscalapack-openmpi` | Ubuntu 包的库名带 `-openmpi` |
| `FFTW_ROOT` | 注释掉，`LLIBS += -lfftw3 -lfftw3_omp`，`INCS += -I/usr/include` | OpenMP 版需要 `fftw3_omp` |
| HDF5 | `CPP_OPTIONS += -DVASP_HDF5`，`LLIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lhdf5_fortran`，`INCS += -I/usr/include/hdf5/openmpi/` | 系统 HDF5 路径 |
| `VASP_TARGET_CPU` | GNU 模板默认 `-march=native` | 跨机器运行时改为具体指令集，如 `-march=znver3` |

GNU 工具链常见的两个坑：

- GCC ≥ 10 对参数类型检查更严格，若出现 `Type mismatch between actual argument` 一类错误，在 `FFLAGS` 中加入 `-fallow-argument-mismatch`（多数 6.x 模板已经带上）；
- 官方 Toolchains 页面说明：GCC > 14.2 删除了部分遗留格式支持，截至 6.5.1 需要在编译选项中加 `-std=legacy` 才能编译。

用 OpenMP 版运行前务必设置 `export OMP_NUM_THREADS=1`（纯 MPI 运行时），否则每个 MPI 进程会按核数开线程，造成严重的超额订阅（oversubscription）。

---

## 6. 可选功能：预编译选项与外部库

所有可选功能都通过“`CPP_OPTIONS += -D...` + 指定库路径 + 链接”三步开启。下表汇总了本教程后续章节会用到的选项（选项名来自官方 Precompiler options 与 makefile.include 页面）：

| 预编译选项 | 功能 / 对应 INCAR 标签 | 外部依赖 | 版本与注意事项 | 用到的章节 |
|---|---|---|---|---|
| `-DVASP_HDF5` | 写 vaspout.h5 / vaspwave.h5；py4vasp | HDF5（Fortran 接口） | ≥ 6.2.0；强烈推荐，部分功能（如电声耦合）强制要求 | 全书 |
| `-DVASP2WANNIER90` | `LWANNIER90`、`WANNIER90_WIN` | Wannier90 的 `libwannier.a` | 6.2 起只支持 Wannier90 2.x/3.x；3.x 要编成**串行库**（make.inc 中不要启用 `COMMS=mpi`，`make lib`） | 第17章 |
| `-Dlibbeef` | `GGA = BF`（BEEF-vdW） | libbeef | GPU 版不支持 | 第8、10章 |
| `-DDFTD4` | `IVDW = 13`（DFT-D4） | dftd4 + multicharge + mctc-lib | DFT-D4 4.x 新 API 需 VASP ≥ 6.6.0；DFT-D4 ≤ 3.7.0 配 6.6 时改用 `-DDFTD4_API_V3` | 第8章 |
| `-DSDFTD3` | `IVDW = 15`（simple-DFT-D3） | s-dftd3 + mctc-lib | VASP ≥ 6.6.0 | 第8章 |
| `-DLIBMBD` | `IVDW = 14`、`LIBMBD_METHOD` | libMBD | 6.4.3 起；libMBD **不要**开 ScaLAPACK/MPI 编译（`-DENABLE_SCALAPACK_MPI=OFF`），否则 NEB 结果错误（known issue 74） | 第8章 |
| `-DUSELIBXC` | Libxc 泛函 | Libxc ≥ 5.2.0，配置时加 `--disable-fhc` | GPU 版不支持 | 第6章 |
| `-DPLUGINS` | Python 插件（`PLUGINS/...` 嵌套标签） | Python 3 + pybind11 | VASP ≥ 6.5；官方建议在 conda 环境（python 3.10）中编译与运行 | — |
| `-Dlibvaspml` | VASPml C++ 机器学习库 | C++17 MPI 编译器 | 实验性；`-DVASPML_ENABLE_GRACE`（≥ 6.6.0）启用 GRACE 力场 | 第15章 |
| `-DELPA` | ELPA 本征值求解器 | ELPA | 可选 | — |
| `-D_OPENMP` | MPI + OpenMP 混合并行 | 编译器 OpenMP 选项（如 `-qopenmp`） | `_omp` 模板已包含 | 第13节 |
| `-Dsol_compat` | VASPsol 兼容开关 | VASPsol 的 `solvation.F` | 见第 10 节 | 第10章 |

下面是一个**同时开启多个可选库**的追加片段示例（每个库对应的 `*_ROOT` 变量名与链接参数取自官方 makefile.include 页面；库目录是 `lib` 还是 `lib64` 取决于你的安装方式，请以实际目录为准）：

```makefile
# ---- optional features, appended after the Customize marker ----
# Wannier90 3.x (serial libwannier.a built with "make lib")
CPP_OPTIONS    += -DVASP2WANNIER90
WANNIER90_ROOT ?= $(HOME)/opt/wannier90-3.1.0
LLIBS          += -L$(WANNIER90_ROOT)/lib -lwannier

# libbeef (BEEF-vdW, GGA = BF); not supported by the GPU ports
CPP_OPTIONS  += -Dlibbeef
LIBBEEF_ROOT ?= $(HOME)/opt/libbeef
LLIBS        += -L$(LIBBEEF_ROOT)/lib -lbeef

# DFT-D4 (IVDW = 13, dftd4 >= 4.0 with VASP >= 6.6.0) and simple-DFT-D3 (IVDW = 15)
CPP_OPTIONS += -DDFTD4 -DSDFTD3
DFTD4_ROOT  ?= $(HOME)/opt/dftd4
SDFTD3_ROOT ?= $(HOME)/opt/s-dftd3
MCTC_ROOT   ?= $(HOME)/opt/mctc-lib
LLIBS       += -L$(DFTD4_ROOT)/lib64 -ldftd4 -lmulticharge
LLIBS       += -L$(SDFTD3_ROOT)/lib64 -ls-dftd3
LLIBS       += -L$(MCTC_ROOT)/lib64 -lmctc-lib
INCS        += -I$(DFTD4_ROOT)/include -I$(SDFTD3_ROOT)/include -I$(MCTC_ROOT)/include

# libMBD (IVDW = 14); build libMBD with -DENABLE_SCALAPACK_MPI=OFF
CPP_OPTIONS += -DLIBMBD
LIBMBD_ROOT ?= $(HOME)/opt/libmbd
LLIBS       += -L$(LIBMBD_ROOT)/build/src -lmbd
INCS        += -I$(LIBMBD_ROOT)/build/src/modules/
```

各外部库的获取与编译要点：

| 库 | 获取 | 编译要点 |
|---|---|---|
| Wannier90 | https://github.com/wannier-developers/wannier90 | 从 `config/` 复制与编译器匹配的 `make.inc`，去掉 `COMMS=mpi`，`make lib` 生成 `libwannier.a`；`-L` 指向该文件所在目录 |
| libbeef | https://github.com/vossjo/libbeef | 用与 VASP 相同的 C 编译器 configure/make/make install |
| mctc-lib、dftd4、simple-dftd3 | https://github.com/grimme-lab/mctc-lib 、https://github.com/dftd4/dftd4 、https://github.com/dftd3/simple-dftd3 | **先装 mctc-lib**，再装 dftd4 与 s-dftd3（meson 或 CMake），都用同一个 Fortran 编译器 |
| libMBD | https://github.com/libmbd/libmbd | CMake，加 `-DENABLE_SCALAPACK_MPI=OFF`；若生成共享库，运行时把 `build/src` 加入 `LD_LIBRARY_PATH`，或用 `-DBUILD_SHARED_LIBS=OFF` 生成静态库 |
| Libxc | https://gitlab.com/libxc/libxc | `./configure --prefix=... --disable-fhc --disable-shared` |

> [!NOTE]
> 开启 `-DPLUGINS` 的步骤（≥ 6.5）：`conda create -n vasp_plugin python=3.10` → `conda activate vasp_plugin` → 在 `src/plugins` 中 `pip install .` → 在 makefile.include 中加入 `-DPLUGINS` 以及由 `python3-config --ldflags --embed` 与 `python3 -m pybind11 --includes` 生成的链接/头文件参数 → 在**同一个 conda 环境中**编译；运行时还要把 `$(python3-config --prefix)/lib` 加入 `LD_LIBRARY_PATH`。

---

## 7. GPU 版编译

### 7.1 NVIDIA GPU：OpenACC + NVHPC

**软件要求**（官方 GPU ports of VASP 页面）：

- NVIDIA HPC SDK ≥ 21.2（避开 22.1 和 22.2：与 OpenMP 同用时有严重 bug，22.3 已修复）；6.6.1 的官方验证工具链之一是 nvhpc-25.1 + CUDA 12.6 + MKL 2025.0.1；
- BLAS/LAPACK/ScaLAPACK、CUDA 组件、CUDA-aware OpenMPI、NCCL（≥ 2.7.8，强烈推荐）均随 NVHPC 提供；**FFTW 需要自己装**（或使用 `nvhpc_ompi_mkl_omp_acc` 模板，由 MKL 提供 CPU 端 FFT/BLAS）；
- 驱动至少支持 CUDA 10.0；用 `nvidia-smi`、`nvcc --version` 查看。

**硬件**：官方测试过 P100、V100、A30/A100、H100 等数据中心卡。消费级显卡（GeForce）技术上可运行，但双精度（FP64）性能弱、显存小且无 ECC，官方不建议用于生产。

**模板与需要修改的项**（完整清单见 [makefile.include.nvhpc_acc.changes](../examples/01-setup-compile/makefile.include.nvhpc_acc.changes)）：

```bash
module load nvhpc/25.1          # 名称因集群而异
cp arch/makefile.include.nvhpc_ompi_mkl_omp_acc makefile.include   # Intel CPU 宿主机推荐
# 或 cp arch/makefile.include.nvhpc_acc makefile.include           # 纯 NVHPC + FFTW
```

| 修改项 | 怎么改 | 说明 |
|---|---|---|
| GPU 架构 `-gpu=ccXX` | 改成你的 GPU 计算能力：V100 → `cc70`，A100/A30 → `cc80`，L40/L4 → `cc89`，H100/H200 → `cc90`；在计算节点上编译可用 `ccnative` | 列多个（如 `cc80,cc90`）可移植性更好但编译更慢 |
| CUDA 版本 `cudaXX.Y` | 与 NVHPC 自带的 CUDA 版本一致 | 不一致会链接失败或运行时报驱动不兼容 |
| `NVROOT` | 模板一般用 `which nvfortran` 自动推导；推导失败时手动设为 `.../hpc_sdk/Linux_x86_64/<版本>` | QD 库（四精度模拟）路径依赖它 |
| `FFTW_ROOT` 或 `MKLROOT` | 指向 CPU 端 FFT 库 | GPU 版仍有部分代码在 CPU 上跑 |
| HDF5 | 同 CPU 版，但 **HDF5 必须用 nvfortran 编译** | Fortran 模块文件与编译器绑定 |
| 预编译选项 | 保留模板默认的 `-DACC_OFFLOAD -DNVCUDA -DUSENCCL -Dqd_emulate` 等 | 6.5.0 前的旧写法是 `-D_OPENACC` |
| 不要添加 | `-Dlibbeef`、`-DUSELIBXC` | GPU 版不支持 |

编译与运行前的环境：

```bash
make veryclean && make DEPS=1 -j 16 all 2>&1 | tee make_gpu.log
# 运行时需要 QD 库路径（NVHPC >= 21.2 要手动加入）
export LD_LIBRARY_PATH=$NVROOT/compilers/extras/qd/lib:$LD_LIBRARY_PATH
```

**GPU 版的功能边界**：除 RPA 相关（GW、ACFDT）外的大部分功能都已移植到 OpenACC；未移植的功能仍可在 GPU 版中使用，只是在 CPU 上执行。GPU 版会**自动把 NCORE 重置为 1**。截至 2026-11 仍开放的 GPU 相关已知问题举例：6.6.0 起 NVIDIA GPU 上杂化泛函配合 ALGO=Normal/DAV 偶发崩溃（issue 99，可改用 ALGO=Damped）；非共线线性响应（LEPSILON、IBRION=7/8）在 GPU 上可能出错（issue 114，改用 CPU 并设 `OMP_NUM_THREADS=1`）。使用前请查看 [Known issues](https://www.vasp.at/wiki/index.php/Known_issues)。

### 7.2 Intel / AMD GPU：OpenMP offload（VASP ≥ 6.6.0，beta）

| 项目 | Intel GPU（Max 1100/1550 等） | AMD GPU（MI210/MI250X/MI300A） |
|---|---|---|
| 模板 | `makefile.include.oneapi_omp_off` | `makefile.include.cray_omp_off` |
| 编译器 | ifx ≥ 2025.3.0（硬性要求） | HPE Cray CCE ≥ 19.0 |
| 数学库 | oneMKL | ROCm（hipFFT、hipBLAS、hipSolver；官方测过 6.1/6.3/6.4） |
| MPI | Intel MPI 或 MPICH（GPU-aware） | cray-mpich ≥ 8 |
| 关键环境变量 | `I_MPI_OFFLOAD=1`、`OMP_STACKSIZE=2048m`；`MKL_NUM_THREADS` 最好不设，让它跟随 `OMP_NUM_THREADS` | `MPICH_GPU_SUPPORT_ENABLED=1`、`OMP_STACKSIZE=2048m`、`OMP_NUM_THREADS=8`、`OMP_WAIT_POLICY=PASSIVE`、MI300A 上 `HSA_XNACK=0` |

6.6.0 在 AMD GPU 上多卡运行需要每个 rank 用 `ROCR_VISIBLE_DEVICES` 只看到自己的 GPU（issue 70，6.6.1 已修复）；6.6.x 在 AMD GPU 上 rocFFT 内核缓存未初始化导致 FFT 规划很慢（issue 115，结果不受影响，官方提供补丁）。线性响应和 BSE 等功能在 OpenMP offload 版中尚未移植。

---

## 8. CMake 构建（VASP ≥ 6.6.0，可选）

6.6.0 起官方提供了 CMake 构建系统，作为传统 makefile.include 的替代（需要 CMake ≥ 3.24），它会自动探测编译器和库：

```bash
cd ~/software/vasp.6.6.1
git clone -b 6.6.x https://github.com/vasp-dev/cmake.git
bash cmake/setup.sh                 # 在源码目录中建立 CMakeLists.txt 的符号链接
mkdir cmake_build && cd cmake_build
cmake .. -DVASP_OPENMP=ON -DVASP_HDF5=ON
make -j 16                          # 默认构建全部；也可 make vasp_std
make test                           # 测试集会被复制到构建目录中运行
```

> [!NOTE]
> 可用的 CMake 选项以 `vasp-dev/cmake` 仓库 README 为准。VTST、VASPsol 这类“改源码”的补丁与 CMake 构建能否顺利配合，取决于补丁是否同时修改了 CMake 的源文件列表；本教程的补丁步骤按传统 makefile 流程描述，更稳妥。

---

## 9. 打 VTST 补丁（CI-NEB、Dimer、IOPT 优化器）

[VTST Tools](https://theory.cm.utexas.edu/vtsttools/) 是 Henkelman 组开发的过渡态工具，包含 CI-NEB、Dimer、Lanczos、动力学矩阵，以及 IOPT 系列优化器（LBFGS、FIRE 等）。VASP 自带的 NEB 没有 climbing image（`LCLIMB` 是 VTST 的标签），所以做过渡态的用户基本都会编译带 VTST 的版本。Fortran 补丁源码在 [henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode)（Apache-2.0），按 VASP 版本分目录（截至 2026-11 有 `vtstcode6.4.3`、`vtstcode6.5.1`、`vtstcode6.6.0`、`vtstcode6.6.1`）。

### 9.1 步骤概览

```bash
cd ~/software
git clone https://github.com/henkelmangroup/vtstcode.git
cd ~/software/vasp.6.6.1
cp -r src src.orig                                  # 1) 备份原始源码
cp -r ../vtstcode/vtstcode6.6.1/* src/              # 2) 复制补丁文件（chain.F 会被替换）
grep -n "chain.o" src/.objects                      # 3) 找到 chain.o 在对象列表中的位置
grep -n "CHAIN_FORCE\|chain_init" src/main.F        # 4) 找到 main.F 中需要改的两处
grep -n "^LIB\|^dependencies" src/makefile          # 5) 找到 src/makefile 中的 LIB 与 dependencies
```

### 9.2 逐项修改（依据 VTST 官方安装说明）

1. **`src/.objects`**：在 `SOURCE` 列表中、`chain.o` **之前**加入以下对象（vtstcode 6.3 及以后的目录）：

   ```text
   bfgs.o dynmat.o instanton.o lbfgs.o sd.o cg.o dimer.o bbm.o \
   fire.o lanczos.o neb.o qm.o \
   pyamff_fortran/*.o ml_pyamff.o \
   opt.o \
   ```

   顺序要求：`dynmat.o neb.o dimer.o lanczos.o instanton.o` 必须在 `chain.o` 之前；各优化器对象（sd、cg、qm、lbfgs、bfgs、fire）以及 `ml_pyamff.o` 和 `pyamff_fortran/*.o` 必须在优化器驱动 `opt.o` 之前。

2. **`src/main.F` 第一处**（VTST v2.04 起需要，用于 solid-state NEB）：在 `CALL CHAIN_FORCE(...)` 的参数列表中，在 `TIFOR,` 之后、`LATT_CUR%A,` 之前插入 `TSIF,`。改完后这一调用的参数顺序应为 `T_INFO%NIONS, DYN%POSION, TOTEN, TIFOR, TSIF, LATT_CUR%A, LATT_CUR%B, IO%IU6`。

3. **`src/main.F` 第二处**（VASP ≥ 6.2）：把 `IF (LCHAIN) CALL chain_init( T_INFO, IO)` 改为无条件调用 `CALL chain_init( T_INFO, IO)`。

4. **`src/makefile`**（vtstcode 6.3 起，因为多了机器学习优化器 PyAMFF 的 Fortran 库）：在定义库目录的 `LIB` 变量末尾追加 `pyamff_fortran`（例如 `LIB= lib parser` → `LIB= lib parser pyamff_fortran`；你的版本中该行若还有其他条目，保留它们，只在末尾追加）。如果并行编译（`-j`），还要把依赖规则改为 `dependencies: sources libs`，保证库先于主程序编译。

5. **重新编译**：

   ```bash
   make veryclean
   make DEPS=1 -j 16 all 2>&1 | tee make_vtst.log
   ```

示例目录中的 [vtst_patch_steps.md](../examples/01-setup-compile/vtst_patch_steps.md) 给出了更详细的检查点，[check_vtst_patch.sh](../examples/01-setup-compile/check_vtst_patch.sh) 可以在编译前自动检查上述 4 处修改是否到位（只读检查，不改文件）。

### 9.3 验证 VTST 是否编进去了

```bash
# 可执行文件中应能找到 VTST 模块的符号（不同编译器的符号修饰不同，忽略大小写搜索即可）
nm bin/vasp_std | grep -i -c "dimer"        # > 0 说明 dimer 模块已链接
nm bin/vasp_std | grep -i -c "fire"         # > 0
```

真正的功能验证放在 [第11章 过渡态搜索](11-neb-dimer.md)：用 `IBRION = 3`、`POTIM = 0`、`IOPT = 1` 跑一个小 NEB，确认能正常进行。VTST 的 Perl/Python 后处理脚本（`nebmake.pl`、`nebresults.pl`、`chgsum.pl` 等）在 [henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts)，克隆后加入 `PATH` 即可：

```bash
git clone https://github.com/henkelmangroup/vtstscripts.git ~/software/vtstscripts
echo 'export PATH=$HOME/software/vtstscripts:$PATH' >> ~/.bashrc
```

> [!CAUTION]
> VTST 只在 INCAR 中设置了 `IMAGES` 或 `ICHAIN` 时才会被调用，不影响其他计算；但它替换了 `chain.F`，**必须使用与 VASP 版本对应的 vtstcode 目录**。用旧版 vtstcode 配新版 VASP 是编译失败最常见的原因。

---

## 10. 打 VASPsol 补丁（隐式溶剂）

[VASPsol](https://github.com/henniggroup/VASPsol)（Hennig 组，Apache-2.0）实现了线性极化连续介质模型（`LSOL`、`EB_K`、`TAU` 等标签），用于电催化与溶液界面（见 [第10章](10-adsorption-catalysis.md)）。VASP 6 自带一个空的 `solvation.F` 骨架，打补丁就是用 VASPsol 的实现替换它。按 VASPsol README 中 “VASP6 versions ≥ 6.1.0” 一节：

```bash
cd ~/software
git clone https://github.com/henniggroup/VASPsol.git
cp VASPsol/src/solvation.F ~/software/vasp.6.6.1/src/     # 覆盖 VASP 自带的骨架文件
grep -n "solvation.o\|pot.o" ~/software/vasp.6.6.1/src/.objects   # 确认 solvation.o 在 pot.o 之前
# 在 makefile.include 中追加：CPP_OPTIONS += -Dsol_compat
# 需要 FERMI_SHIFT 输出时，再按 README 打 VASPsol/patches/ 中对应 6.x 的补丁
cd ~/software/vasp.6.6.1 && make veryclean && make DEPS=1 -j 16 all
```

> [!NOTE]
> 原版 VASPsol 仓库已基本停止开发，新功能（非线性/非局域电解质模型、恒电势）转移到 [VASPsol++](https://github.com/VASPsol/VASPsol)。VASPsol++ 的源码补丁是**按 VASP 版本提供的**（README 中列出的是 5.4.4 与 6.3.2，分“带 VTST”和“不带 VTST”两种），需要邮件向作者索取；用新版 VASP 前请确认是否有对应补丁。VTST 与原版 VASPsol 修改的是不同文件，可以同时打。

---

## 11. 验证安装：官方测试集与 Si 小算例

### 11.1 官方测试集（testsuite）

```bash
cd ~/software/vasp.6.6.1
# 默认用 "mpirun -np 4 <root>/bin/vasp_*" 运行；如需修改：
export VASP_TESTSUITE_EXE_STD="mpirun -np 4 $PWD/bin/vasp_std"
export VASP_TESTSUITE_EXE_GAM="mpirun -np 4 $PWD/bin/vasp_gam"
export VASP_TESTSUITE_EXE_NCL="mpirun -np 4 $PWD/bin/vasp_ncl"
export OMP_NUM_THREADS=1
make test            # FAST 类测试，官方估计 4 核约 1.5 小时
# make test_all      # 全部测试，耗时很长
grep -c ERROR testsuite/testsuite.log     # 0 表示全部通过
tail -30 testsuite/testsuite.log          # 失败的测试会列在文件末尾
```

常用的筛选变量（均以 `VASP_TESTSUITE_` 开头）：`VASP_TESTSUITE_TESTS="bulk_GaAs_ACFDT"` 只跑指定测试；`VASP_TESTSUITE_SKIP_WAN90=Y` 跳过需要 Wannier90 的测试；`VASP_TESTSUITE_SKIP_GW=Y`、`VASP_TESTSUITE_SKIP_HYB=Y` 等跳过某类测试。在集群上建议把 `make test` 写进一个作业脚本提交，而不是在登录节点运行。

### 11.2 Si 小算例（2 分钟内完成）

示例目录 [test_Si/](../examples/01-setup-compile/README.md) 提供了 2 原子金刚石 Si 的 INCAR/POSCAR/KPOINTS（POTCAR 需自备）。它的用途是**快速确认程序能跑、并行结果一致**：

```bash
cd VASP/examples/01-setup-compile/test_Si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR      # 路径见第2章
mpirun -np 1 vasp_std > out.np1 2>&1 && grep "energy(sigma->0)" OUTCAR | tail -1
mpirun -np 4 vasp_std > out.np4 2>&1 && grep "energy(sigma->0)" OUTCAR | tail -1
grep "General timing and accounting" OUTCAR          # 出现即为正常结束
```

判断标准：

- `OUTCAR` 末尾出现 `General timing and accounting informations for this job`；
- 1 核与 4 核的 `energy(sigma->0)` 差异在 $10^{-5}$ eV 量级以内（与 EDIFF 同量级）；
- 2 原子 Si（PBE，a = 5.43 Å）的总能量在 −10.8 eV 左右（经验量级，仅供参考，以你自己的结果为准）。

示例 README 中还给出了 vasp_gam（Γ 点）与 vasp_ncl 的一致性测试。

---

## 12. 运行 VASP 与作业提交

### 12.1 本机或交互式节点运行

```bash
ulimit -s unlimited                 # 解除栈大小限制，避免段错误
export OMP_NUM_THREADS=1            # 纯 MPI 运行
export OMP_STACKSIZE=512m           # 若用 OpenMP 线程
mpirun -np 4 vasp_std > vasp.out 2>&1 &
tail -f OSZICAR                     # 实时查看电子步/离子步
```

> [!WARNING]
> 不要在集群登录节点上直接跑 VASP。登录节点是共享资源，管理员通常会杀掉长时间占用 CPU 的进程。

### 12.2 Slurm：CPU 作业

下面是完整的单节点 CPU 作业脚本（示例文件 [job_cpu.slurm](../examples/01-setup-compile/job_cpu.slurm)）。分区名、模块名、可执行文件路径需按你的集群修改；更通用的模板见仓库脚本 [`job_vasp_slurm.sh`](../scripts/job_vasp_slurm.sh)（用法见脚本内说明或 `bash ../scripts/job_vasp_slurm.sh -h`）。

```bash
#!/bin/bash
#SBATCH --job-name=vasp_cpu
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ---------------- environment (site specific) ----------------
module purge
module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
export PATH=$HOME/software/vasp.6.6.1/bin:$PATH
VASP_EXE=${VASP_EXE:-vasp_std}          # vasp_std / vasp_gam / vasp_ncl

ulimit -s unlimited
export OMP_NUM_THREADS=1

# ---------------- sanity checks ----------------
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in INCAR POSCAR POTCAR; do
    [ -s "$f" ] || { echo "ERROR: $f is missing or empty"; exit 1; }
done
if [ ! -s KPOINTS ] && ! grep -qi "^ *KSPACING" INCAR; then
    echo "WARNING: no KPOINTS file and no KSPACING in INCAR (VASP will use KSPACING = 0.5)"
fi
if grep -q $'\t' INCAR || grep -U -q $'\r' INCAR; then
    echo "WARNING: INCAR contains TAB or CR characters; run dos2unix INCAR"
fi

echo "Job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS MPI ranks, start $(date)"
mpirun -np "$SLURM_NTASKS" "$VASP_EXE" > vasp.out 2>&1
echo "End $(date)"

# ---------------- quick post-run check ----------------
if grep -q "General timing and accounting" OUTCAR 2>/dev/null; then
    echo "VASP finished normally"
else
    echo "VASP did NOT finish normally, check vasp.out and OUTCAR"
fi
grep "reached required accuracy" OUTCAR | tail -1
grep "energy(sigma->0)" OUTCAR | tail -1
```

| 行 | 作用 |
|---|---|
| `--ntasks-per-node=64` | 每节点 MPI 进程数 = 物理核数（不要算超线程） |
| `--cpus-per-task=1` | 每个 MPI 进程 1 个核（纯 MPI） |
| `ulimit -s unlimited` | 官方建议：避免大数组栈分配导致段错误 |
| `mpirun -np $SLURM_NTASKS` | Intel MPI / OpenMPI 都能识别 Slurm 分配；若集群要求用 `srun` 启动，请按集群文档配置 PMI |
| 提交前检查 | 缺 POTCAR 或 KPOINTS 是新手最常见的“秒退” |

提交与监控：`sbatch job_cpu.slurm`、`squeue -u $USER`、`scancel <jobid>`、`sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS`。

### 12.3 Slurm：NVIDIA GPU 作业

GPU 版的运行规则（官方 GPU ports 页面）：**每块 GPU 只用 1 个 MPI rank**，剩余 CPU 核用 OpenMP 线程；INCAR 中 `KPAR` 设为 GPU 数（前提是不可约 k 点数大致能被整除），`NSIM` 比 CPU 版大得多（默认 4 太小，可试 16–32），NCORE 会被自动重置为 1。完整脚本见 [job_gpu.slurm](../examples/01-setup-compile/job_gpu.slurm)：

```bash
#!/bin/bash
#SBATCH --job-name=vasp_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4          # = number of GPUs per node (1 rank per GPU)
#SBATCH --cpus-per-task=8            # OpenMP threads per rank
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out

module purge
module load nvhpc/25.1
export PATH=$HOME/software/vasp.6.6.1-gpu/bin:$PATH
# Root of the HPC SDK (the variable name differs between module systems)
NVROOT=${NVHPC_ROOT:-$(dirname "$(dirname "$(dirname "$(which nvfortran)")")")}
export LD_LIBRARY_PATH=$NVROOT/compilers/extras/qd/lib:$LD_LIBRARY_PATH

ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export MKL_THREADING_LAYER=INTEL     # only for builds linked against MKL

cd "$SLURM_SUBMIT_DIR"
grep -E "^ *(KPAR|NSIM|NCORE)" INCAR || echo "Hint: set KPAR = number of GPUs and NSIM = 16-32"
nvidia-smi --query-gpu=name,memory.total --format=csv

mpirun -np "$SLURM_NTASKS" --map-by ppr:"$SLURM_NTASKS_PER_NODE":node:PE="$SLURM_CPUS_PER_TASK" \
       --bind-to core --report-bindings vasp_std > vasp.out 2>&1
```

> [!TIP]
> `NVHPC_ROOT` 变量名因 module 而异（有的叫 `NVHPC`、`NVHPC_HOME`），所以脚本在变量不存在时从 `nvfortran` 的位置倒推 SDK 根目录（`.../Linux_x86_64/<版本>`）。`--report-bindings` 会在输出中打印每个 rank 绑定到哪些核，第一次运行时务必检查一下。

### 12.4 PBS / Torque 作业

示例文件 [job_pbs.sh](../examples/01-setup-compile/job_pbs.sh)（仓库通用模板为 [`job_vasp_pbs.sh`](../scripts/job_vasp_pbs.sh)）：

```bash
#!/bin/bash
#PBS -N vasp_cpu
#PBS -q workq
#PBS -l nodes=1:ppn=32
## PBS Pro syntax instead of the line above:
## #PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -j oe

module purge
module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
export PATH=$HOME/software/vasp.6.6.1/bin:$PATH

cd "$PBS_O_WORKDIR"
ulimit -s unlimited
export OMP_NUM_THREADS=1
NP=$(wc -l < "$PBS_NODEFILE")
mpirun -np "$NP" -machinefile "$PBS_NODEFILE" vasp_std > vasp.out 2>&1
```

PBS 常用命令：`qsub job_pbs.sh`、`qstat -u $USER`、`qdel <jobid>`。

### 12.5 运行环境变量清单

| 变量 | 推荐值 | 说明 |
|---|---|---|
| `ulimit -s` | `unlimited` | 主线程栈大小 |
| `OMP_NUM_THREADS` | 纯 MPI：`1`；混合/GPU：每 rank 线程数 | 务必显式设置 |
| `OMP_STACKSIZE` | CPU：`512m`；OpenMP offload GPU：`2048m` | 每个 OpenMP 线程的私有栈 |
| `OMP_PLACES` / `OMP_PROC_BIND` | `cores` / `close` | 线程绑定在父 rank 附近 |
| `I_MPI_PIN_DOMAIN` | `omp`（Intel MPI + OpenMP） | 为每个 rank 预留 OMP_NUM_THREADS 个核 |
| `MKL_THREADING_LAYER` | `INTEL`（链接 MKL 的 GPU/混合版本） | 保证 MKL 使用 Intel OpenMP 运行时 |
| `LD_LIBRARY_PATH` | 包含 HDF5、QD、libMBD 等共享库目录 | 运行时找不到 `.so` 会立即退出 |

---

## 13. 并行参数：KPAR、NCORE、NSIM

### 13.1 三个标签的含义

| 标签 | 默认值 | 含义 | 经验规则（官方建议 + 实践） |
|---|---|---|---|
| `KPAR` | 1 | 把 MPI 进程分成 KPAR 组，每组负责一部分 k 点 | 内存允许时增大；应能整除总进程数，并尽量整除不可约 k 点数；最好让每个节点/NUMA 域各自负责一组 k 点 |
| `NCORE` | 1 | 每个能带由多少个进程合作完成 FFT | 取每节点核数的因子，约 $\sqrt{N_{\text{ranks}}/\text{KPAR}}$；设为一个 NUMA 域的核数往往很好；小体系/慢网络用 1；约 100 原子用 4，400 原子以上可到 12–16 |
| `NSIM` | 4 | RMM-DIIS / Davidson 中同时处理的能带数 | CPU 上一般不改；GPU 上增大到 16–32 |
| `NPAR` | 可用进程数 | 旧标签，NPAR × NCORE = 可用进程数 | **不要与 NCORE 同时设置**（同时出现时 NPAR 优先） |

> [!WARNING]
> 使用 OpenMP 线程（OMP_NUM_THREADS > 1）或任何 GPU 版本时，VASP 会**自动把 NCORE 重置为 1**，此时 FFT 并行由线程或 GPU 完成；在 INCAR 里写 NCORE 不会生效。

### 13.2 先做“空跑”，再选参数

用 `ALGO = None`（不做电子优化，只完成 FFT 规划和 k 点约化）快速拿到 NBANDS、NKPTS 和平面波数：

```bash
cp INCAR INCAR.prod
echo "ALGO = None" >> INCAR             # 临时加一行
mpirun -np 64 vasp_std > dry.out 2>&1
grep -E "NKPTS|NBANDS" OUTCAR | head -3   # 不可约 k 点数与能带数
grep "dimension x,y,z NGX" OUTCAR         # FFT 网格
mv INCAR.prod INCAR
```

然后挑 3–4 组候选参数，各跑几个电子步（如 `NELM = 6`、`NSW = 0`），比较 `LOOP:` 计时：

```bash
grep "LOOP:" OUTCAR | awk '{s+=$NF; n++} END {print "mean real time per SCF step:", s/n, "s"}'
```

### 13.3 参考配置

| 体系 | 硬件 | 建议起点 | 理由 |
|---|---|---|---|
| 2 原子 Si，8×8×8 Γ 网格（29 个不可约 k 点） | 1 节点 32 核 | 小原胞时 `NCORE = 1`，KPAR 取能整除 32 的较大值（如 8 或 16） | 能带少、k 点多，k 点并行最有效；29 是素数，可按官方技巧把 IBZKPT 复制为 KPOINTS 并补零权重点凑成整除 |
| 100 原子表面，Γ + 少量 k 点 | 2 节点 × 64 核，每 8 核一个 NUMA 域 | `KPAR = 2`（每节点一组）、`NCORE = 8` | k 点组不跨节点，FFT 组不跨 NUMA 域 |
| 300 原子体系，只有 Γ 点 | 4 节点 × 64 核 | `vasp_gam`、`KPAR = 1`、`NCORE = 8`–`16` | 只能并行能带与 FFT |
| 50 原子体系 | 1 节点 4 × A100 | 4 个 rank、每 rank 8 线程、`KPAR = 4`、`NSIM = 16`–`32` | 每 GPU 一个 rank，k 点并行 |

> [!NOTE]
> 进程数改变时 VASP 会自动增大 NBANDS 使其能被并行组数整除，并在 OUTCAR/标准输出中给出 `The number of bands has been changed` 一类提示。普通计算无妨，但做 LOBSTER、LOPTICS、MLFF 等对 NBANDS 敏感的计算时，应固定并行设置并手动设置 NBANDS（见相关章节）。

---

## 参数选择与收敛建议（编译与运行速查表）

| 项目 | 推荐 | 说明 |
|---|---|---|
| VASP 版本 | 6.6.1 | 6.6.0 需打 pot.F.patch |
| CPU 模板 | `oneapi`（Intel）、`gnu_ompi_aocl`（AMD）、`gnu_omp`（个人电脑） | 与源码同版本 |
| GPU 模板 | `nvhpc_ompi_mkl_omp_acc`（Intel 宿主）/ `nvhpc_acc` | NVHPC ≥ 21.2，避开 22.1/22.2 |
| 必开选项 | `-DVASP_HDF5` | py4vasp、vaspout.h5 |
| 常用可选 | `-DVASP2WANNIER90`、`-DDFTD4`/`-DSDFTD3`、`-Dlibbeef`、VTST、VASPsol | 按研究方向 |
| 编译命令 | `make veryclean; make DEPS=1 -j N all` | 换配置先清理 |
| 安装验证 | `make test` + Si 小算例 | 并行结果一致 |
| CPU 运行 | 纯 MPI、`OMP_NUM_THREADS=1`、`ulimit -s unlimited` | — |
| GPU 运行 | 1 rank/GPU、OpenMP 线程、`KPAR` = GPU 数、`NSIM` = 16–32 | NCORE 自动为 1 |
| NCORE（CPU） | 每节点核数的因子，≈ NUMA 域大小 | 小体系用 1 |
| KPAR | 整除进程数与（尽量）整除不可约 k 点数 | 受内存限制 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `ifort: command not found`（编译时） | oneAPI 2025 起不再提供 ifort，旧模板或 `mpiifort` 默认调用 ifort | 用 6.4.3+ 的 `makefile.include.oneapi`（ifx）；把 FC/FCL 改为 `mpiifort -fc=ifx` 或 `mpiifx` |
| 编译 `fftw3d.F` 时找不到 `fftw3.f` | `INCS` 没有包含 FFTW 头文件目录 | MKL：`INCS += -I$(MKLROOT)/include/fftw`；FFTW：`INCS += -I$(FFTW_ROOT)/include` |
| `Type mismatch between actual argument ...`（gfortran） | GCC ≥ 10 参数检查更严格 | `FFLAGS += -fallow-argument-mismatch` |
| 编译 HDF5 相关文件时提示无法读取 `hdf5.mod`（或模块“由不同版本编译器生成”） | HDF5 由另一种/另一版本编译器编译 | 用编译 VASP 的同一 Fortran 编译器重编 HDF5 |
| 链接时 `undefined reference to ...h5...` | `LLIBS` 缺 `-lhdf5_fortran`（或还需 `-lhdf5`）、库路径不对 | 检查 `HDF5_ROOT/lib` 下的库名，补全 `-L` 与 `-l` |
| 链接时找不到 `wannier_setup` 等符号 | 设置了 `-DVASP2WANNIER90` 但没有链接 `libwannier.a`，或 W90 3.x 编成了 MPI 版本 | 检查 `WANNIER90_ROOT`，按串行方式 `make lib` |
| 运行时 `error while loading shared libraries: libxxx.so: cannot open shared object file` | 共享库不在 `LD_LIBRARY_PATH` | `export LD_LIBRARY_PATH=...` 或链接时加 `-Wl,-rpath,...`；用 `ldd bin/vasp_std | grep "not found"` 排查 |
| `forrtl: severe (174): SIGSEGV, segmentation fault occurred` | 栈空间不足、内存不足、或输入有严重问题 | `ulimit -s unlimited`、`OMP_STACKSIZE=512m`；增加节点/减少每节点进程数；检查 POSCAR |
| `mpirun noticed that process rank N ... exited on signal 9 (Killed)` 或 Slurm 日志出现 `oom-kill` | 内存不足被系统杀掉 | 减少每节点进程数、降低 KPAR、用 `vasp_gam`、增加节点 |
| 编译带 VTST 的版本时报 `CHAIN_FORCE` 参数个数/类型不匹配 | `main.F` 没有按要求加入 `TSIF` | 按第 9.2 节修改 main.F，`make veryclean` 后重编 |
| VTST 编译时找不到 PyAMFF 相关模块或对象 | `.objects` 中顺序不对，或 `src/makefile` 的 `LIB` 未加 `pyamff_fortran` | 按第 9.2 节第 1、4 步修改；并行编译时设 `dependencies: sources libs` |
| GPU 版速度比 CPU 版还慢 | 每块 GPU 跑了多个 rank、体系太小、NSIM 太小 | 1 rank/GPU；KPAR = GPU 数；NSIM 调大；小体系标准 DFT 本来就不适合 GPU |
| 6.6.0 做自旋极化 r2SCAN 优化时结构“越优化越怪” | known issue 94：ISPIN=2 meta-GGA 受力错误 | 升级 6.6.1 或打 pot.F.patch |
| 纯 MPI 运行却比预期慢很多、负载异常高 | OpenMP 版程序未设 `OMP_NUM_THREADS=1`，线程超额订阅 | 作业脚本中显式 `export OMP_NUM_THREADS=1` |

---

## 练习

1. 在你的集群上查出：每节点物理核数、NUMA 域数量、可用的 oneAPI / GCC / NVHPC 版本，并据此为 VASP 6.6.1 选一个模板，写出你需要修改的变量清单。
2. 编译带 HDF5 的 vasp_std，运行示例 Si 算例，分别用 1、2、4、8 个进程运行并比较 `energy(sigma->0)` 与 `LOOP:` 时间，画出加速比曲线。
3. 对同一个 30–50 原子的体系，比较 `NCORE = 1, 2, 4, 8` 和 `KPAR = 1, 2, 4` 的组合，找出最快的设置。
4. （有 GPU 的读者）比较 1 块 GPU（1 rank × 8 线程）与 1 个 CPU 节点在同一算例上的速度，并测试 `NSIM = 4, 16, 32` 的影响。
5. 打 VTST 补丁并运行第 9.3 节的检查命令，确认 dimer/fire 相关符号已链接。

---

## 参考资料

**官方文档**

- VASP Manual 总入口：https://www.vasp.at/wiki/index.php/The_VASP_Manual
- Installing VASP.6.X.X（makefile 与 CMake 两种方式）：https://www.vasp.at/wiki/index.php/Installing_VASP.6.X.X
- makefile.include（模板列表与 HDF5/Wannier90/Libxc/libbeef/DFT-D4/libMBD/Plugins/VASPml 配置）：https://www.vasp.at/wiki/index.php/Makefile.include
- Precompiler options：https://www.vasp.at/wiki/index.php/Precompiler_options
- Compiler options / Linking to libraries：https://www.vasp.at/wiki/index.php/Compiler_options ，https://www.vasp.at/wiki/index.php/Linking_to_libraries
- Toolchains（官方验证过的编译器与库版本）：https://www.vasp.at/wiki/index.php/Toolchains
- Personal computer installation（Debian/Ubuntu/Fedora/Rocky/macOS）：https://www.vasp.at/wiki/index.php/Personal_computer_installation
- GPU ports of VASP：https://www.vasp.at/wiki/index.php/GPU_ports_of_VASP
- Validation tests（testsuite 与环境变量）：https://www.vasp.at/wiki/index.php/Validation_tests
- Combining MPI and OpenMP：https://www.vasp.at/wiki/index.php/Combining_MPI_and_OpenMP
- Optimizing the parallelization / NCORE / KPAR / NSIM：https://www.vasp.at/wiki/index.php/Optimizing_the_parallelization
- Changelog 与 Known issues：https://www.vasp.at/wiki/index.php/Changelog ，https://www.vasp.at/wiki/index.php/Known_issues
- VASP 新闻（版本发布、AI 工具使用规范）：https://vasp.at/info/post/
- VTST 安装说明：https://theory.cm.utexas.edu/vtsttools/installation.html
- py4vasp 文档：https://vasp.at/py4vasp/latest/

**相关 GitHub 仓库**（仅链接，不转载代码）

- [henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode)：VTST Fortran 补丁源码，按 VASP 版本分目录（Apache-2.0）。
- [henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts)：nebmake.pl、nebresults.pl、chgsum.pl 等 VTST 脚本。
- [henniggroup/VASPsol](https://github.com/henniggroup/VASPsol) 与 [VASPsol/VASPsol](https://github.com/VASPsol/VASPsol)：隐式溶剂模型及其新版 VASPsol++。
- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90)：`-DVASP2WANNIER90` 需要的 libwannier.a。
- [dftd4/dftd4](https://github.com/dftd4/dftd4)、[dftd3/simple-dftd3](https://github.com/dftd3/simple-dftd3)：IVDW=13/15 的外部库。
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：官方 Python 后处理库，依赖 HDF5 版 VASP。
- [northword/wiki](https://github.com/northword/wiki)：中文学习笔记站，含 VTST 版 VASP 编译记录，可对照参考。
- [ishikawa-group/vasp-tutorial](https://github.com/ishikawa-group/vasp-tutorial)：附编译笔记的循序渐进教程（日文）。

**关键文献**

- G. Kresse, J. Furthmüller, Phys. Rev. B 54, 11169 (1996). DOI: 10.1103/PhysRevB.54.11169
- G. Kresse, J. Furthmüller, Comput. Mater. Sci. 6, 15 (1996). DOI: 10.1016/0927-0256(96)00008-0
- P. E. Blöchl, Phys. Rev. B 50, 17953 (1994). DOI: 10.1103/PhysRevB.50.17953
- G. Kresse, D. Joubert, Phys. Rev. B 59, 1758 (1999). DOI: 10.1103/PhysRevB.59.1758
- G. Henkelman, B. P. Uberuaga, H. Jónsson, J. Chem. Phys. 113, 9901 (2000). DOI: 10.1063/1.1329672
- K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, R. G. Hennig, J. Chem. Phys. 140, 084106 (2014). DOI: 10.1063/1.4865107
- F. Wende, M. Marsman, J. Kim, F. Vasilev, Z. Zhao, T. Steinke, Int. J. Quantum Chem. 118, e25851 (2018). DOI: 10.1002/qua.25851

[返回目录](../README.md) | [下一章 →](02-input-output-files.md)
