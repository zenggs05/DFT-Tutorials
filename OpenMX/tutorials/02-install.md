# 第2章 编译与安装（Intel / GCC / ScaLAPACK / ELPA / WSL / deb / Docker）（Compiling and Installing OpenMX 4.0）

本章给出 OpenMX 4.0 的完整安装方案：从源码编译（Intel oneAPI 的 `mpiicx`/`mpiifx` + MKL 路线，以及 GCC + OpenMPI + ScaLAPACK + OpenBLAS + FFTW 路线），打 patch 4.0.1，理解内置 ELPA 本征值求解器，用 `-runtest` 验证安装；也介绍免编译的官方 Debian 包（WSL/Ubuntu）、Docker 镜像、MateriApps LIVE! 与 macOS（Apple Silicon）的安装要点，以及在 HPC 集群上部署 OpenMX（OpenMX installation, build from source, makefile CC/FC/LIB）时的检查方法和常见编译错误。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 OpenMX 4.0 的依赖：MPI（必需）、ScaLAPACK + BLACS 或 oneAPI MKL（必需）、FFTW3（或 MKL 自带的 FFTW3 接口），以及**源码内置**的 ELPA；
- 正确修改 `source/makefile` 中的 `CC`、`FC`、`LIB`，用 Intel oneAPI 或 GCC 工具链编译出 `openmx` 和全部后处理程序；
- 打 patch 4.0.1，并用 `Methane.dat` 和 `-runtest` 判断安装是否正确；
- 在没有编译条件时，用官方 deb 包或 Docker 镜像快速开始计算；
- 看懂并解决最常见的编译、链接和运行时库错误。

| 你的情况 | 推荐方式 | 本章小节 |
|---|---|---|
| Windows 笔记本，只想尽快开始 | WSL2 + Ubuntu 24.04 + **官方 deb 包** | §8.1 |
| Windows/WSL，以后要改 makefile、编后处理工具或用 API | WSL2 + Intel oneAPI 源码编译（4.0 默认 makefile 就是为此配置的） | §4 |
| Linux 工作站，不想装 Intel 工具链 | GCC + OpenMPI + ScaLAPACK + OpenBLAS + FFTW | §5 |
| macOS（Apple Silicon） | MacPorts 编译（论坛帖）或 Docker 镜像 | §5.5、§8.2 |
| 非 Ubuntu 的 Linux、想隔离环境 | Docker 镜像（x86_64/ARM） | §8.2 |
| 超算集群（生产计算） | 用集群提供的 Intel 或 GCC 模块自行编译 | §9 |
| 需要复现 3.9 的旧结果 | 另装一份 3.9 + patch 3.9.9 | §10 |

> [!NOTE]
> Linux 基础操作、conda、SSH 与文件传输见公共篇 [第1章 计算环境搭建](../../Common/tutorials/01-linux-hpc.md)；Slurm/PBS 作业脚本见公共篇 [第2章 HPC 作业调度](../../Common/tutorials/02-job-scheduler.md)。

---

## 2. 原理速览：OpenMX 的软件栈与构建方式

### 2.1 依赖关系

| 组件 | 是否必需 | Intel 路线 | GCC 路线 | 说明 |
|---|---|---|---|---|
| C 编译器 | 必需 | `icx`（经 `mpiicx` 包装） | `gcc`（经 `mpicc` 包装） | 主程序是 C |
| Fortran 编译器 | 必需 | `ifx`（经 `mpiifx` 包装） | `gfortran`（经 `mpif90`/`mpifort` 包装） | 内置 ELPA 是 Fortran |
| MPI | **必需** | Intel MPI | Open MPI（或 MPICH） | 4.0 **不支持**串行版 |
| ScaLAPACK + BLACS | **必需** | MKL（`-lmkl_scalapack_lp64`、`-lmkl_blacs_*`） | `libscalapack-openmpi` | 不支持只用 LAPACK/BLAS 编译 |
| BLAS/LAPACK | 必需 | MKL | OpenBLAS | — |
| FFTW3 | 必需 | MKL 自带 FFTW3 接口（`-I${MKLROOT}/include/fftw`） | `libfftw3` | 4.0 只支持 FFTW3 |
| ELPA | 必需 | **源码内置**（`source/elpa-2018.05.001`、`elpa1.f90` 等） | 同左 | 不需要外部 libelpa |
| OpenMP | 可选 | `-qopenmp` | `-fopenmp` | 没有 OpenMP 时在 `CC` 中加 `-Dnoomp` |
| 系统 `gcc` | `make all` 需要 | — | — | makefile 用 `gcc` 直接编译 `bandgnu13`、`cif2omx`、`gcube2oned`、`diff_gcube` 等小工具 |

### 2.2 内置 ELPA 与运行时选择

OpenMX 4.0 用 ELPA（Eigenvalue SoLvers for Petascale Applications）做高度并行的本征值求解，ELPA 的源码随 OpenMX 一起编译，**不需要也不应该**链接系统里的 libelpa。运行时用关键词选择 ELPA1 或 ELPA2：

```text
scf.eigen.lib      elpa1    # elpa1|elpa2; default for Cluster calculations = elpa1
scf.eigen.lib      elpa2    # default for Band calculations since Ver. 4.0
```

手册的基准测试表明，分子（Cluster）计算中两者速度相近；能带（Band）计算中 ELPA2 略快，因此 4.0 起能带计算默认 `elpa2`。一般不需要改动。

### 2.3 makefile 的结构

OpenMX 没有 configure/CMake，只有一个 `source/makefile`。你需要关心的只有三个（Intel 路线还有 `MKLROOT`）变量：

```text
CC  = C compiler (MPI wrapper) + optimization + OpenMP + include paths
FC  = Fortran compiler (MPI wrapper) + optimization + OpenMP
LIB = ScaLAPACK/BLACS + BLAS/LAPACK + FFTW3 + Fortran runtime + MPI Fortran bindings + -lm -ldl ...
```

makefile 的开头部分既有一组**生效的**定义，也有大量被 `#` 注释掉的各平台示例（ISSP 超算、WSL 笔记本、AMD EPYC 等）；后面的目标文件列表之后还有 `CC += -I$(LIBELPADIR)`、`FC += -I$(LIBELPADIR)` 两行，用来把内置 ELPA 的目录加入头文件搜索路径。修改时只动生效的那一组 `MKLROOT`/`CC`/`FC`/`LIB`，不要动 `+=` 这两行。

> [!WARNING]
> 不要用 `make install CC="..." FC="..." LIB="..."` 这种命令行覆盖的方式。GNU make 中命令行变量会让 makefile 里对同名变量的所有赋值（**包括 `+=` 追加**）失效，`-I$(LIBELPADIR)` 就加不上，编译 ELPA 时会出错。请直接编辑 makefile（或用本章示例脚本自动编辑）。

### 2.4 make 目标

| 命令 | 作用 |
|---|---|
| `make install` | 编译 `openmx` 并复制到 `../work/` |
| `make all` | 编译 `openmx` 和全部后处理程序（`DosMain`、`cohp`、`jx`、`polB`、`calB`、`Z2FH`、`kSpin`、`bandgnu13`、`cif2omx`、`gcube2oned`、`diff_gcube`、`md2axsf`、`cube2xsf`、`bin2txt`、`esp`、`analysis_example` 等，以 makefile 中的 `UTIL` 为准），并全部复制到 `../work/` |
| `make DosMain`、`make cohp`、`make polB` ... | 只编译某一个后处理程序 |
| `make clean` | 删除目标文件、`.mod` 与可执行文件；**改动编译选项后一定要先 `make clean`** |

---

## 3. 下载与打补丁

### 3.1 下载

```bash
mkdir -p ~/software && cd ~/software
wget https://www.openmx-square.org/openmx4.0.tar.gz                       # ~182 MB
wget https://www.openmx-square.org/bugfixed/26May08/patch4.0.1.tar.gz
tar zxvf openmx4.0.tar.gz
ls openmx4.0                                                              # source  work  DFT_DATA19
```

> [!TIP]
> 在 WSL 中，请把源码和计算目录放在 Linux 文件系统（如 `~/software`、`~/work`）下，不要放在 `/mnt/c/...`。跨文件系统访问 Windows 盘的 I/O 很慢，会明显拖慢编译和计算。

### 3.2 打 patch 4.0.1

按补丁 README 的步骤（补丁包含 `Band_DFT_Dosout.c`、`Mulliken_Charge.c`、`GaAs.dat` 三个文件）：

```bash
cp ./patch4.0.1.tar.gz openmx4.0/source
cd openmx4.0/source
tar zxvf patch4.0.1.tar.gz        # overwrites Band_DFT_Dosout.c and Mulliken_Charge.c
mv GaAs.dat ../work/              # corrected input for manual Sec. 38.1
```

之后正常编译即可；如果此前已经编译过，需要 `make clean` 后重新 `make install`、`make all`。补丁修正的内容见 [第1章 §3.4](01-overview.md)。

> [!IMPORTANT]
> patch 4.0.1 把 bxsf/frmsf 费米面文件的单位从 Hartree/Bohr 改成了 eV/Å。打补丁前后得到的费米面文件不能混用同一套后处理脚本。

---

## 4. 路线 A：Intel oneAPI（WSL/Ubuntu 或集群）

4.0 的默认 makefile 就是为“WSL + Intel oneAPI”配置的，使用 LLVM 架构的新编译器 `icx`/`ifx` 及其 MPI 包装 `mpiicx`/`mpiifx`。

### 4.1 安装 oneAPI（个人电脑/WSL）

需要 C/C++ 编译器、Fortran 编译器、MKL 和 Intel MPI。用 Intel 官方 APT 源安装（包名以 Intel 官方安装文档为准）：

```bash
# add Intel's APT repository
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update

# components needed by OpenMX (or install the whole Base + HPC toolkits)
sudo apt install -y intel-oneapi-compiler-dpcpp-cpp intel-oneapi-compiler-fortran \
                    intel-oneapi-mkl-devel intel-oneapi-mpi-devel
sudo apt install -y make gcc            # system gcc is needed by 'make all' for small tools

# load the environment (add this line to ~/.bashrc if you use it often)
source /opt/intel/oneapi/setvars.sh
which mpiicx mpiifx icx ifx
echo $MKLROOT
```

> [!NOTE]
> Intel 在 oneAPI 2024 中移除了经典 C 编译器 `icc`，在 2025 中移除了经典 Fortran 编译器 `ifort`。网上大量 3.9 时代的教程使用 `mpiicc`/`mpiifort` + `-ip -no-prec-div`，在新版 oneAPI 上已不可用，请改用 `mpiicx`/`mpiifx`。

### 4.2 makefile 设置

打开 `openmx4.0/source/makefile`，找到**未被注释**的 `MKLROOT`、`CC`、`FC`、`LIB` 定义（可能用 `\` 续行跨多行），整块注释掉，然后写入下面的内容（与本章示例中的 [`makefile.oneapi.example`](../examples/02-install-build/makefile.oneapi.example) 相同）：

```make
MKLROOT = /opt/intel/oneapi/mkl/latest
CC  = mpiicx -O3 -xHOST -fcommon -Wno-implicit-function-declaration -qopenmp \
      -I${MKLROOT}/include -I${MKLROOT}/include/fftw
FC  = mpiifx -O3 -xHOST -qopenmp
LIB = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread \
      -lmkl_core -lifcore -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl
```

逐项解释：

| 选项 | 作用 | 何时修改 |
|---|---|---|
| `MKLROOT` | MKL 安装目录。makefile 中的赋值会覆盖环境变量 `$MKLROOT` | 若 oneAPI 不在 `/opt/intel/oneapi`，改成 `echo $MKLROOT` 的输出 |
| `mpiicx` / `mpiifx` | Intel MPI 包装的 icx/ifx | 集群若用 Open MPI + Intel 编译器，改为 `mpicc`/`mpifort` 并确认其底层是 icx/ifx |
| `-O3` | 优化等级 | — |
| `-xHOST` | 针对**编译所在机器**的 CPU 指令集优化 | AMD CPU 改为 `-march=x86-64-v3`（手册建议）；登录节点与计算节点 CPU 不同时改用通用选项 |
| `-fcommon` | 允许多个源文件中出现同名全局变量的“暂定定义” | 必需，否则链接报 `multiple definition` |
| `-Wno-implicit-function-declaration` | 关闭“隐式函数声明”诊断（新编译器默认把它当错误） | 必需；仍报错时可再加 `-Wno-error=implicit-int` |
| `-qopenmp` | 启用 OpenMP | 不要 OpenMP 时删去并在 CC 中加 `-Dnoomp` |
| `-I${MKLROOT}/include/fftw` | 使用 MKL 自带的 FFTW3 接口头文件 | 用外部 FFTW3 时改为其 include 目录，并在 LIB 中加 `-lfftw3` |
| `-lmkl_scalapack_lp64 ... -lmkl_core` | MKL 的 ScaLAPACK、LP64 接口、OpenMP 线程层与核心库 | 一般不改 |
| `-lmkl_blacs_intelmpi_lp64` | 与 **Intel MPI** 配套的 BLACS | 用 Open MPI 时改为 `-lmkl_blacs_openmpi_lp64` |
| `-lifcore` | ifx 的 Fortran 运行时库 | 仅 Intel 编译器路线需要 |
| `-liomp5` | Intel OpenMP 运行时 | — |

4.0 手册给出的默认 `CC`/`FC` 还带有 `-ipo`（过程间优化，链接时间明显变长）和 `-DLEAK_DETECT`（开发者用的内存泄漏检测，运行结束时打印每个进程的分配/泄漏统计）。生产使用时两者都可以去掉；若想与官方默认完全一致，也可以保留。

> [!TIP]
> MKL 的链接行与 MPI 实现、整数接口（LP64/ILP64）、线程层有关。拿不准时用 [Intel oneMKL Link Line Advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html) 生成：选 Linux、你的编译器、LP64、OpenMP 线程（Intel 编译器选 Intel OpenMP）、勾选 ScaLAPACK 和 BLACS 并选择 MPI 实现。

### 4.3 编译

```bash
cd ~/software/openmx4.0/source
source /opt/intel/oneapi/setvars.sh
make clean
make install 2>&1 | tee make_install.log      # builds ../work/openmx
make all     2>&1 | tee make_all.log          # builds DosMain, cohp, jx, polB, bandgnu13, ... and copies them to ../work
ls -l ../work/openmx ../work/DosMain ../work/bandgnu13
```

编译 `openmx` 在普通笔记本上大约需要十几分钟到半小时（带 `-ipo` 时更久）。保守起见不要加 `-j`：有用户报告并行 make 时会因 ELPA 的 Fortran 模块依赖顺序而失败；若你想尝试并行编译，失败后 `make clean` 再串行编译即可。

编译成功的判断：

```bash
grep -iE " error|error:" make_install.log | head        # should print nothing
ldd ../work/openmx | grep "not found"                   # should print nothing
```

### 4.4 一键脚本

本章示例 [02-install-build](../examples/02-install-build/README.md) 提供 `build_openmx_oneapi.sh`：自动下载、打补丁、备份并改写 makefile（AMD CPU 自动把 `-xHOST` 换成 `-march=x86-64-v3`）、编译并可选地运行 runtest。仓库脚本 [`omx_build.sh`](../scripts/omx_build.sh) 提供同样功能的通用版本，用 `bash ../scripts/omx_build.sh -h` 查看用法。

---

## 5. 路线 B：GCC + OpenMPI + ScaLAPACK + OpenBLAS + FFTW

这条路线全部使用发行版软件包，不需要 Intel 工具链，适合 Ubuntu/Debian 工作站和 WSL。

### 5.1 安装依赖（Ubuntu 24.04 / Debian 12–13）

```bash
sudo apt update
sudo apt install -y build-essential gfortran \
    openmpi-bin libopenmpi-dev \
    libscalapack-openmpi-dev libopenblas-dev libfftw3-dev \
    wget
mpicc --version; mpif90 --version; mpirun --version
```

> [!IMPORTANT]
> Ubuntu/Debian 的 ScaLAPACK 库文件名是 **`libscalapack-openmpi.so`**（来自 `libscalapack-openmpi-dev`），并没有 `libscalapack.so`。所以链接时要写 `-lscalapack-openmpi`；写成 `-lscalapack` 会报 `cannot find -lscalapack`。其他发行版（如 Fedora/RHEL 的 `scalapack-openmpi-devel`、Spack、自己编译的 ScaLAPACK）通常是 `-lscalapack`，并需要用 `-L` 指定目录。

### 5.2 makefile 设置

同样先注释掉 makefile 中生效的 `MKLROOT`/`CC`/`FC`/`LIB`，再写入（即示例中的 [`makefile.gcc.example`](../examples/02-install-build/makefile.gcc.example)）：

```make
CC  = mpicc  -O3 -fopenmp -fcommon -Wno-implicit-function-declaration
FC  = mpif90 -O3 -fopenmp -fallow-argument-mismatch
LIB = -lscalapack-openmpi -lopenblas -lfftw3 -lmpi_mpifh -lgfortran -lpthread -lm -ldl
```

| 选项 | 作用 |
|---|---|
| `-fcommon` | GCC ≥ 10 默认 `-fno-common`，不加会在链接时报 `multiple definition of 'SpeBasisName'` 等 |
| `-Wno-implicit-function-declaration` | OpenMX 部分源文件调用函数前没有原型声明；GCC ≥ 14 默认把它当错误 |
| `-fallow-argument-mismatch` | gfortran ≥ 10 对 MPI/BLAS 调用的实参类型检查更严，编译内置 ELPA（`elpa1.f90` 等）会报 `Type mismatch`，此选项把错误降为警告 |
| `-lscalapack-openmpi` | Debian/Ubuntu 的 ScaLAPACK（含 BLACS） |
| `-lopenblas` | BLAS + LAPACK（Debian 的 OpenBLAS 包含优化过的 LAPACK） |
| `-lfftw3` | FFTW3（头文件在 `/usr/include`，无需 `-I`） |
| `-lmpi_mpifh` | Open MPI 的 Fortran 绑定；主程序用 `mpicc` 链接 Fortran 写的 ELPA 时需要 |
| `-lgfortran` | gfortran 运行时库 |

可选优化：在**编译机与计算机 CPU 相同**时，可在 `CC`、`FC` 中加 `-march=native`；ARM（aarch64）平台需要在 `CC` 中加 `-Dnosse`（O(N) Krylov 例程默认使用 x86 的 SSE 指令）。

> [!NOTE]
> 不同 GCC 版本可能还需要额外的兼容选项。GCC ≥ 14 把隐式 int、不兼容指针类型等也变成了默认错误，如果遇到相应报错，可追加 `-Wno-error=implicit-int -Wno-error=incompatible-pointer-types`；GCC ≥ 15 默认使用 C23 标准，如遇到与旧式函数声明有关的新错误，可追加 `-std=gnu17`。示例脚本 `build_openmx_gcc.sh` 会根据 `mpicc -dumpversion` 自动添加这些选项。

### 5.3 编译与检查

```bash
cd ~/software/openmx4.0/source
make clean
make install 2>&1 | tee make_install.log
make all     2>&1 | tee make_all.log
ldd ../work/openmx | grep -E "scalapack|openblas|fftw|mpi|not found"
```

`ldd` 输出中应能看到 `libscalapack-openmpi.so`、`libopenblas.so`、`libfftw3.so`、`libmpi.so`，且没有 `not found`。

> [!TIP]
> OpenMX 自己用 OpenMP 做线程并行。用 OpenBLAS 编译时，建议在运行脚本中设置 `export OPENBLAS_NUM_THREADS=1`，避免 BLAS 内部再开线程造成超额订阅（oversubscription）。

### 5.4 其他 Linux 发行版

| 发行版 | 依赖包（参考，具体包名以发行版仓库为准） | 注意 |
|---|---|---|
| Fedora / RHEL 系（需 EPEL） | `gcc gcc-gfortran openmpi-devel scalapack-openmpi-devel openblas-devel fftw-devel` | 先 `module load mpi/openmpi-x86_64`；库通常在 `/usr/lib64/openmpi/lib`，LIB 写 `-L/usr/lib64/openmpi/lib -lscalapack` |
| Arch / Manjaro | `gcc-fortran openmpi scalapack openblas fftw` | `-lscalapack` |
| 无 root 的服务器 | Spack：`spack install netlib-scalapack fftw openblas` 后 `spack load` | 用 `spack location -i <pkg>` 找到 `-L`/`-I` 路径 |

### 5.5 macOS（Apple Silicon）

官方论坛帖 [How to Install OpenMX Ver. 4.0 on macOS (Apple Silicon)](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3524) 给出了基于 MacPorts 的完整步骤，要点是：

1. 用 MacPorts 安装带 gcc13 + OpenBLAS + Open MPI 变体的 `scalapack`、`openmpi-gcc13`，以及带 OpenMP/MPI 变体的 `fftw-3`；
2. `CC` 使用 MacPorts 的 `mpicc-openmpi-gcc13`，并加 `-Dnosse -fcommon -fopenmp` 和 `-I/opt/local/include`；
3. `FC` 使用 `mpif90-openmpi-gcc13` 并加 `-fallow-argument-mismatch -fopenmp`；
4. `LIB` 用 `-L/opt/local/lib` 链接 `scalapack`、`lapack`、`blas`、`fftw3`、`gfortran`、`mpi_mpifh`；
5. `make all` 后在 `work/` 下用 `mpirun-openmpi-gcc13 -np 2 openmx -runtest` 验证。

请以论坛原帖的具体命令为准。不想编译的 Mac 用户可以直接用 §8.2 的 Docker 镜像（支持 ARM）。

---

## 6. 编译产物与环境设置

### 6.1 产物一览

编译成功后，`openmx4.0/work/` 下会出现：

| 文件 | 来源 | 用途 |
|---|---|---|
| `openmx` | `make install` | 主程序 |
| `DosMain`、`cohp` | `make all` | DOS/PDOS、COHP/COOP |
| `bandgnu13` | `make all`（makefile 中直接调用 `gcc bandgnu13.c -lm`） | 能带数据转 gnuplot 格式 |
| `polB`、`jx`、`calB`、`Z2FH`、`kSpin` | `make all` | 极化、交换耦合、Berry 曲率/Chern 数、Z2、自旋织构 |
| `cif2omx` | `make all` | CIF → dat |
| `gcube2oned`、`diff_gcube`、`add_gcube`、`cube2xsf`、`md2axsf`、`bin2txt` | `make all` | cube/轨迹/二进制格式处理 |
| `analysis_example`、`esp` 等 | `make all` | `.scfout` 读取示例、ESP 电荷 |

手册中也给出了单独编译小工具的方法，例如 `gcc bandgnu13.c -lm -o bandgnu13`，然后把可执行文件复制到需要的目录。

### 6.2 把程序加入 PATH

`work/` 里有大量示例输入，直接把它加入 `PATH` 也可以，但更整洁的做法是建一个 `bin/`：

```bash
OMX=~/software/openmx4.0
mkdir -p $OMX/bin
cd $OMX/work
for f in openmx DosMain cohp bandgnu13 polB jx calB Z2FH kSpin cif2omx \
         gcube2oned diff_gcube add_gcube cube2xsf md2axsf bin2txt esp analysis_example; do
  [ -x "$f" ] && cp "$f" $OMX/bin/
done
echo 'export OPENMX_ROOT=$HOME/software/openmx4.0' >> ~/.bashrc
echo 'export PATH=$OPENMX_ROOT/bin:$PATH'           >> ~/.bashrc
```

示例目录中的 `env_openmx.sh` 把这些环境设置（含 Intel 运行库、栈大小、线程数）写在一起，需要时 `source` 即可。

> [!IMPORTANT]
> OpenMX **不读取**任何“数据目录”环境变量。在自己的工作目录中计算时，必须在每个 `.dat` 里写 `DATA.PATH`（绝对路径），例如 `DATA.PATH  /home/user/software/openmx4.0/DFT_DATA19`。

---

## 7. 验证安装

### 7.1 第一步：Methane.dat

手册第 5 节规定的标准测试。`work/Methane.dat` 使用默认的 `DATA.PATH ../DFT_DATA19`，所以要在 `work/` 下运行：

```bash
cd ~/software/openmx4.0/work
export OMP_NUM_THREADS=1
mpirun -np 1 ./openmx Methane.dat -nt 1 | tee met.std
grep "Utot\." met.out
tail -n 3 met.std
```

判断标准：

- `met.std` 最后一行是 `The calculation was normally finished.`；
- 生成 `met.out`、`met.xyz`、`met.ene`、`met.md`、`met.md2`、`met.cif`、`met.tden.cube`、`met.v0.cube`、`met.vhart.cube`、`met.dden.cube` 和目录 `met_rst/`；
- 手册第 5 节给出的参考总能量为 `Utot. -8.033540953061` Hartree（16 步 SCF 收敛，单核 3.6 GHz Xeon 约 5 秒）。你的结果应与之在 $10^{-8}$ Hartree 量级以内一致。

### 7.2 第二步：自动测试 `-runtest`

```bash
cd ~/software/openmx4.0/work
mpirun -np 4 ./openmx -runtest -nt 1          # MPI only
# mpirun -np 4 ./openmx -runtest -nt 2        # MPI/OpenMP hybrid (needs 8 cores)
cat runtest.result
```

`-runtest` 依次计算 `work/input_example/` 中的 14 个输入（Benzene、C60、CO、Cr2、Crys-MnO、GaAs、Glycine、Graphite4、H2O-EF、H2O、HMn、Methane、Mol_MnO、Ndia2），并把总能量与力相对参考数据的绝对差写入 `runtest.result`。每个算例有一组 `Elapsed time(s)`、`diff Utot`、`diff Force`。

判断标准（手册第 6 节）：**`diff Utot` 和 `diff Force` 在小数点后第 7 位以内（即 < 1×10⁻⁷）就说明安装正确**。由于浮点运算顺序依赖于机器和进程数，看到 1×10⁻¹¹ 量级的非零差值是正常的。

```text
（示意输出，数值仅供参考）
1
input_example/Benzene.dat
Elapsed time(s)= 4.93
diff Utot= 0.000000000000
diff Force= 0.000000000000
...
Total elapsed time (s) 99.01
```

可以用示例目录中的 `check_runtest.sh` 自动判断：

```bash
bash /path/to/OpenMX/examples/02-install-build/check_runtest.sh runtest.result      # threshold 1e-7 by default
```

`runtest.result` 中的 `Elapsed time` 也可以用来比较不同编译选项、不同机器的速度；`work/input_example/` 中还附带了几个平台上生成的 `runtest.result` 供参考。

### 7.3 可选：大体系基准 `-runtestL`

```bash
mpirun -np 112 ./openmx -runtestL -nt 2
```

`-runtestL` 计算 `work/large_example/` 中的 16 个大体系（如 DIA512、Pt63、Ni63-O64 等）并写入 `runtestL.result`。手册中用 128 MPI × 2 线程也需要约 20 分钟，且核数少时可能因内存不足而段错误，**只适合在集群上做性能基准**，不是安装验证的必需步骤。

> [!NOTE]
> 源码中还提供了针对特定功能的测试选项（如 `-runtestNEGF`、`-runtestWF`、`-runtestCWF`、`-runtestCDDF`），对应 `work/` 下的 `negf_example`、`wf_example`、`cwf_example`、`cddf_example` 等目录。在正式使用 NEGF、Wannier、光学等功能前，可以先跑相应测试确认结果与参考一致（这些选项在手册中介绍较少，以源码和对应目录中的 `runtest*.result` 为参考）。

### 7.4 并行与线程检查

- `-nt` 指定每个 MPI 进程的 OpenMP 线程数，不写时为 1。建议同时设置 `export OMP_NUM_THREADS=<与 -nt 相同>`，避免库函数自行开满线程。
- 进程数 × 线程数不要超过物理核数；在笔记本上跑 runtest 时，`-np 4 -nt 1` 是稳妥的选择。
- Open MPI 报 `There are not enough slots available in the system ...` 说明请求的进程数超过了它识别到的核数，减少 `-np`，或在确认核数足够时加 `--oversubscribe`。

---

## 8. 免编译安装

### 8.1 官方 Debian 包（Ubuntu/Debian/WSL）

官方 deb 包面向 x86_64（Intel 与 AMD CPU 均可），并**自带运行 OpenMX 所需的 oneAPI 运行库**，不需要事先安装 Intel oneAPI。手册第 4.11 节给出三个版本：

```bash
# (a) Ubuntu 24.04/26.04 and Debian 12/13, with Open MPI
wget -O ./openmx_4.0.1-1_amd64.deb https://www.openmx-square.org/deb/openmx_4.0.1-1_amd64.deb
sudo apt install -y ./openmx_4.0.1-1_amd64.deb

# (b) Ubuntu 24.04, with Intel oneAPI
wget -O ./openmx4.0.1-wsl_amd64.deb https://www.openmx-square.org/deb/openmx4.0.1-wsl_amd64-Ubuntu24.04.deb
sudo apt install -y ./openmx4.0.1-wsl_amd64.deb

# (c) Ubuntu 26.04, with Intel oneAPI
wget -O ./openmx4.0.1-wsl_amd64.deb https://www.openmx-square.org/deb/openmx4.0.1-wsl_amd64-Ubuntu26.04.deb
sudo apt install -y ./openmx4.0.1-wsl_amd64.deb
```

安装后程序和工具位于 `/opt/openmx4.0`，`openmx`、`DosMain` 等可直接用命令名调用。官方建议的首次测试：

```bash
mkdir -p ~/openmx_work && cd ~/openmx_work
cp /opt/openmx4.0/work/Methane.dat ./
ls /opt/openmx4.0                            # locate DFT_DATA19 (usually /opt/openmx4.0/DFT_DATA19)
sed -i 's#^ *DATA.PATH.*##' Methane.dat      # remove an old DATA.PATH line if present
echo "DATA.PATH   /opt/openmx4.0/DFT_DATA19" >> Methane.dat
openmx ./Methane.dat | tee met.std           # serial run first
grep "Utot\." met.out
```

确认与手册第 5 节结果一致后，再尝试 MPI 并行（如 `mpirun -np 2 openmx ./Methane.dat | tee met.std`）。手册说明：对于 oneAPI 版 deb 包，MPI 并行需要另外安装 Intel oneAPI 环境；先用非 MPI 方式确认安装，再按需使用并行。

> [!WARNING]
> 运行 deb 包里的 `DosMain` 等后处理程序时，如果出现
> `error while loading shared libraries: libmpifort.so.12: cannot open shared object file`，
> 先把包内的运行库目录加入搜索路径，再**直接用命令名**调用（不要写成 `./DosMain`）：
>
> ```bash
> export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
> DosMain met.Dos.val met.Dos.vec
> ```

想在 deb 安装上跑 `-runtest`：`/opt` 下的 `work/` 一般没有写权限，可以把 `work/` 复制到家目录，并在旁边建一个指向数据库的符号链接，使默认的 `../DFT_DATA19` 依然有效：

```bash
mkdir -p ~/omx_runtest && cd ~/omx_runtest
cp -r /opt/openmx4.0/work .
ln -s /opt/openmx4.0/DFT_DATA19 DFT_DATA19     # adjust if your DFT_DATA19 is elsewhere
cd work && mpirun -np 2 openmx -runtest
```

查看或卸载已安装的包：`dpkg -l | grep -i openmx`，然后 `sudo apt remove <包名>`。

### 8.2 Docker 镜像（Linux / macOS，x86_64 与 ARM）

官方论坛帖 [Created a Docker image for OpenMX Ver. 4.0](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3522) 发布了基于 Ubuntu 22.04、同时支持 x86_64 和 ARM 的镜像：

```bash
mkdir -p ~/host_openmx_work
docker run --shm-size=4gb --rm -it -v ~/host_openmx_work:/work dc1394/openmx4.0-ubuntu22.04:0.1

# inside the container
cd openmx_work
mpirun -np 2 openmx Methane.dat | tee met.std
cp met.std met.out /work/          # /work is the bind-mounted host directory
exit
```

要点：

- `--shm-size=4gb`：MPI 进程间通信使用共享内存，Docker 默认的 `/dev/shm` 很小，不加会导致 MPI 出错；
- `--rm` 表示退出后删除容器，**结果必须保存到挂载目录 `/work`**；
- 在容器内处理自己的体系时，把 `.dat` 放在 `/work` 下计算，`DATA.PATH` 指向容器内的 DFT_DATA19（可用 `find / -maxdepth 4 -type d -name DFT_DATA19 2>/dev/null` 查找）；
- 发布者实测同一台机器上 runtest 在容器内的吞吐约为宿主机的 70%，适合入门和中小体系，大规模生产计算仍建议原生编译。

在没有 Docker、但有 Apptainer/Singularity 的集群上，可以尝试 `apptainer pull openmx4.sif docker://dc1394/openmx4.0-ubuntu22.04:0.1` 后用 `apptainer exec openmx4.sif ...` 在单节点内运行（这是通用的容器用法，跨节点 MPI 需要额外配置，未经官方文档验证）。

### 8.3 MateriApps LIVE!

[MateriApps LIVE!](https://cmsi.github.io/MateriAppsLive/) 是预装了 OpenMX 等计算材料软件的 Linux 环境，提供 VirtualBox 版、Docker 版等，适合课堂演示或零基础体验。其中的 OpenMX 版本可能不是最新的，用于正式计算前请检查版本。

### 8.4 GPU 版

- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：官方下载页链接的 4.0 GPU 版（v1.1.1，2026-08-25），用 CMake 构建，可选 CUDA、HIP 或 CPU 后端，详细的构建选项与基准见仓库 README；
- [dc1394/openmx3.9.9_gpu](https://github.com/dc1394/openmx3.9.9_gpu)：3.9.9 的 GPU 移植版。

GPU 版不是官方 CPU 版的一部分，使用前请用 runtest 或自己的体系与 CPU 版结果对比。

---

## 9. 在 HPC 集群上编译与部署

### 9.1 选择与加载模块

```bash
module avail 2>&1 | grep -iE "intel|oneapi|mkl|openmpi|gcc|fftw|scalapack"
# Intel route (names differ between clusters)
module load intel-oneapi-compilers intel-oneapi-mkl intel-oneapi-mpi
# or GCC route
module load gcc openmpi netlib-scalapack openblas fftw
module list
```

### 9.2 检查 MPI 包装器

编译前确认 `mpicc` 背后是哪个编译器、链接了哪个 MPI：

```bash
# Open MPI
mpicc --showme:command      # underlying C compiler
mpicc --showme:compile      # added compile flags
mpicc --showme:link         # added link flags
# Intel MPI / MPICH
mpiicx -show
mpicc -show
```

常见的组合与 `LIB` 中 BLACS 的对应关系：

| 编译器 | MPI | MKL BLACS 库 |
|---|---|---|
| icx/ifx | Intel MPI（`mpiicx`/`mpiifx`） | `-lmkl_blacs_intelmpi_lp64` |
| icx/ifx | Open MPI（`mpicc`/`mpifort` 指向 icx/ifx） | `-lmkl_blacs_openmpi_lp64`，可能还需 `-lmpi_mpifh` 等 Open MPI Fortran 库 |
| gcc/gfortran | Open MPI + MKL | `-lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp`（建议用 Link Line Advisor 生成） |
| gcc/gfortran | Open MPI + 开源库 | 不用 MKL，按 §5.2 |

`source/makefile` 开头的注释中保留了开发者在 ISSP 超算等机器上使用过的完整设置，是很好的参考。

### 9.3 部署建议

- 把 OpenMX 装在共享文件系统（如 `/home` 或 `/apps`）上，计算节点才能访问；`DATA.PATH` 写绝对路径；
- 登录节点与计算节点 CPU 不同时，不要用 `-xHOST`/`-march=native`，改用目标机器支持的指令集（如 `-march=core-avx2`、`-march=x86-64-v3`）；
- 运行时加载与编译时**相同**的模块，并用编译时同一套 MPI 的 `mpirun`/`srun` 启动；
- 作业脚本中设置 `ulimit -s unlimited`、`export OMP_NUM_THREADS=<-nt 的值>`、`export OMP_STACKSIZE=512M`（大体系或 OpenMP 运行出现段错误时尤其重要）；
- 在计算节点上（交互作业）跑一次 `-runtest`，确认与登录节点结果一致。

作业脚本模板见仓库 [`job_omx_slurm.sh`](../scripts/job_omx_slurm.sh)、[`job_omx_pbs.sh`](../scripts/job_omx_pbs.sh) 和 [第18章](18-parallel-on-errors.md)。

---

## 10. 需要 3.9 时：并行安装旧版本

为了复现用 3.9 发表的结果，可以把 3.9 装在另一个目录，与 4.0 互不干扰：

```bash
cd ~/software
wget https://www.openmx-square.org/openmx3.9.tar.gz
wget https://www.openmx-square.org/bugfixed/21Oct17/patch3.9.9.tar.gz
tar zxvf openmx3.9.tar.gz
cp patch3.9.9.tar.gz openmx3.9/source && cd openmx3.9/source
tar zxvf patch3.9.9.tar.gz
# follow the README of patch 3.9.9 for files that must be moved to ../work
```

3.9 同样需要 MPI + ScaLAPACK/BLACS；编译选项与 4.0 基本相同（GCC 路线同样需要 `-fcommon` 与 `-fallow-argument-mismatch`；Intel 新编译器同样需要 `-fcommon -Wno-implicit-function-declaration`，也可参考 [atsushi-m-ito/openmx-patch-oneapi](https://github.com/atsushi-m-ito/openmx-patch-oneapi) 给出的 icx/ifx 补丁）。两个版本的 `DFT_DATA19` 内容基本相同，但请让每个版本使用各自目录下的数据库。

---

## 参数选择与收敛建议（编译选项速查表）

| 场景 | CC | FC | LIB 要点 |
|---|---|---|---|
| WSL / Intel CPU / oneAPI ≥ 2024 | `mpiicx -O3 -xHOST -fcommon -Wno-implicit-function-declaration -qopenmp -I${MKLROOT}/include -I${MKLROOT}/include/fftw` | `mpiifx -O3 -xHOST -qopenmp` | MKL ScaLAPACK + `-lmkl_blacs_intelmpi_lp64 -lifcore -liomp5` |
| AMD CPU / oneAPI | 同上，`-xHOST` → `-march=x86-64-v3` | 同上 | 同上 |
| Ubuntu 24.04 / GCC 13 | `mpicc -O3 -fopenmp -fcommon -Wno-implicit-function-declaration` | `mpif90 -O3 -fopenmp -fallow-argument-mismatch` | `-lscalapack-openmpi -lopenblas -lfftw3 -lmpi_mpifh -lgfortran` |
| GCC ≥ 14 | 再加 `-Wno-error=implicit-int -Wno-error=incompatible-pointer-types`（遇到相应错误时） | 同上 | 同上 |
| GCC ≥ 15 | 再加 `-std=gnu17`（遇到 C23 相关错误时） | 同上 | 同上 |
| ARM（aarch64，Linux/macOS） | 加 `-Dnosse` | 同上 | 按平台的 ScaLAPACK/BLAS |
| 无 OpenMP | 去掉 `-qopenmp`/`-fopenmp`，加 `-Dnoomp` | 去掉 OpenMP 选项 | 去掉 `-liomp5` 等 |
| 开发调试 | 加 `-DLEAK_DETECT`（内存泄漏统计） | — | — |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `multiple definition of 'SpeBasisName'; ... first defined here` | GCC ≥ 10 默认 `-fno-common`，`openmx_common.h` 中的全局变量在多个目标文件中重复定义 | `CC` 加 `-fcommon`，`make clean` 后重编 |
| `Error: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8))`（编译 `elpa1.f90` 等） | gfortran ≥ 10 的实参类型检查 | `FC` 加 `-fallow-argument-mismatch` |
| `error: implicit declaration of function '...' [-Wimplicit-function-declaration]` | icx 或 GCC ≥ 14 把隐式声明当错误 | `CC` 加 `-Wno-implicit-function-declaration`（必要时再加 `-Wno-error=implicit-int`）；icx/ifx 也可参考 atsushi-m-ito 的补丁 |
| `mpiicc: command not found`、`icc: not found`、`ifort: command not found` | oneAPI 2024/2025 起已移除 icc/ifort | 改用 `mpiicx`/`mpiifx`，去掉 `-ip -no-prec-div` 等旧选项 |
| `/usr/bin/ld: cannot find -lifcore` | 用 GCC 工具链却保留了 Intel 的库 | GCC 路线去掉 `-lifcore`、`-liomp5`、`-lmkl_blacs_intelmpi_lp64` |
| 去掉 `-lifcore` 后出现 `undefined reference to 'mpi_abort_'`、`'mpi_sendrecv_'` | 缺少 Open MPI 的 Fortran 绑定 | 加 `-lmpi_mpifh -lgfortran`；MKL 配 Open MPI 时 BLACS 用 `-lmkl_blacs_openmpi_lp64` |
| `/usr/bin/ld: cannot find -lscalapack` | Debian/Ubuntu 的库名是 `libscalapack-openmpi.so` | 改为 `-lscalapack-openmpi`（或自建符号链接，不推荐） |
| `fatal error: fftw3.h: No such file or directory` | 没装 FFTW 开发包，或 MKL 路线缺 `-I${MKLROOT}/include/fftw` | `apt install libfftw3-dev`，或补上 include 路径 |
| `makefile:NN: *** missing separator.  Stop.` | 编辑 makefile 时破坏了续行（`\` 后的行没有一起注释），或 recipe 行的 Tab 被换成空格 | 用备份的 makefile 重新修改；续行要整体注释；不要用会把 Tab 转空格的编辑器 |
| 编译 ELPA 文件时找不到 `.h`/`.X90` 等包含文件 | 用 `make CC=...` 命令行覆盖导致 `CC += -I$(LIBELPADIR)` 失效 | 改为编辑 makefile 中的定义 |
| `make all` 在 `bandgnu13` 处失败：`gcc: command not found` | 只装了 Intel 编译器，系统没有 gcc | `sudo apt install gcc`，或手动用 icx 编译该工具 |
| 运行时 `error while loading shared libraries: libmkl_*.so` / `libiomp5.so` | 运行环境没有加载 oneAPI | `source /opt/intel/oneapi/setvars.sh`，或在作业脚本中 `module load`；也可在 LIB 中加 `-Wl,-rpath,${MKLROOT}/lib/intel64` |
| deb 包：`libmpifort.so.12: cannot open shared object file`（运行 DosMain 等） | 工具找不到包内运行库 | `export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH`，直接用命令名 `DosMain` 调用 |
| `mpirun has detected an attempt to run as root.` | 在 Docker 等环境中以 root 运行 Open MPI | 用普通用户运行；确有需要时加 `--allow-run-as-root` |
| `runtest.result` 中 `diff Utot` 大于 1e-7 | 编译器过度优化（如激进的 fast-math）、库不匹配、混用了不同版本的 MPI | 先用 `-O2`/`-O3` 且不加 fast-math 类选项重编；检查 `ldd`；确认 `mpirun` 与编译所用 MPI 一致 |
| runtest/大体系运行中段错误（`Segmentation fault`） | 栈空间或内存不足 | `ulimit -s unlimited`、`export OMP_STACKSIZE=512M`，减少 `-nt`，增加节点 |
| 运行时多个进程各自从头计算、输出重复 | 用了与编译时不同的 MPI 的 `mpirun`（每个进程都以为自己是 rank 0） | 用 `which mpirun` 检查，加载与编译一致的 MPI 模块 |

---

## 练习

1. 在 WSL 或 Linux 上按 §5 用 GCC 路线编译 OpenMX 4.0 + patch 4.0.1，记录 `make install` 用时，并用 `ldd` 列出 `openmx` 链接的数学库。
2. 运行 `Methane.dat`，把你的 `Utot` 与手册参考值 −8.033540953061 Hartree 比较，计算差值。
3. 分别用 `-np 1 -nt 1`、`-np 2 -nt 1`、`-np 4 -nt 1`、`-np 2 -nt 2` 运行 `-runtest`，从 `runtest.result` 中比较 `Total elapsed time`，并确认所有 `diff` 都小于 1e-7。
4. （有 Intel 工具链的读者）再用 oneAPI 路线编译一份，比较两份可执行文件在 runtest 上的速度。
5. 安装官方 deb 包，用 `DosMain` 验证 `LD_LIBRARY_PATH` 的设置是否必要。
6. 阅读 `source/makefile` 开头注释中的各平台示例，找出 ISSP 超算上使用 Open MPI + MKL 时选用的 BLACS 库名。

---

## 参考资料

**官方文档**

- 安装（第 4 节目录）：https://www.openmx-square.org/openmx_man4.0/s4_installation.html
- 4.1 依赖库：https://www.openmx-square.org/openmx_man4.0/s4_1_required-libraries.html
- 4.3 MPI 版与默认 makefile：https://www.openmx-square.org/openmx_man4.0/s4_3_mpi-version.html
- 4.4 MPI/OpenMP 混合版：https://www.openmx-square.org/openmx_man4.0/s4_4_mpi-openmp-hybrid.html
- 4.8 安装提示（`mpicc --showme`、MKL Link Line Advisor）：https://www.openmx-square.org/openmx_man4.0/s4_8_hints-installation.html
- 4.10 WSL 源码安装：https://www.openmx-square.org/openmx_man4.0/s4_10_wsl-install.html
- 4.11 Debian 包安装：https://www.openmx-square.org/openmx_man4.0/s4_11_deb-install.html
- 第 5 节测试计算（Methane.dat 参考输出）：https://www.openmx-square.org/openmx_man4.0/s5_test-calc.html
- 第 6、7 节自动测试 runtest / runtestL：https://www.openmx-square.org/openmx_man4.0/s6_auto-test.html ，https://www.openmx-square.org/openmx_man4.0/s7_auto-test-large.html
- 第 28 节并行（`-nt`、`scf.eigen.lib`）：https://www.openmx-square.org/openmx_man4.0/s28_parallel.html
- 下载页：https://www.openmx-square.org/download.html ；patch 4.0.1 README：https://www.openmx-square.org/bugfixed/26May08/README.txt
- WSL 安装说明（3.9 时代 PDF，仍可参考 WSL/oneAPI 部分）：https://www.openmx-square.org/Installation_OpenMX_WSL_EN.pdf
- 论坛：macOS（Apple Silicon）安装 https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3524 ；Docker 镜像 https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3522
- Intel oneMKL Link Line Advisor：https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html

**相关 GitHub 仓库**（只链接，不转载其内容）

- [mfukudaQED/openmx_installer](https://github.com/mfukudaQED/openmx_installer)：下载 + 打补丁 + 改 makefile + 编译的一键脚本，给出 oneAPI 2025 之后 `mpiicx`/`mpiifx` 需要的 `-fcommon -Wno-error=implicit-function-declaration -Wno-error=implicit-int` 等关键选项。
- [atsushi-m-ito/openmx-patch-oneapi](https://github.com/atsushi-m-ito/openmx-patch-oneapi)：icc 被移除后用 icx/ifx 编译 OpenMX 的补丁与说明（日文），4.0 手册第 81 节也有收录。
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：`scripts/` 中有 Ubuntu 24.04 上 GCC 与 Intel 两条路线的 4.0 + patch 4.0.1 编译脚本（仅供参考阅读）。
- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：4.0 GPU 版（CMake；CUDA/HIP/CPU 后端）。
- [cometscome/OpenMX_Docker](https://github.com/cometscome/OpenMX_Docker)：最小 Dockerfile 示例（需要 `--shm-size`），可作为自建镜像的参考。
- [Ncmexp2717/OpenMX-deb](https://github.com/Ncmexp2717/OpenMX-deb)：早期 3.9 的非官方 deb 包（数据路径改为 `/usr/share/openmx/DFT_DATA19`），4.0 已有官方 deb，可作历史对照。
- [cmsi/MateriAppsLive](https://github.com/cmsi/MateriAppsLive)：预装 OpenMX 的 Linux 环境。
- [hikaruri/DFT_startup](https://github.com/hikaruri/DFT_startup)：macOS 上用 MacPorts 编译 OpenMX 3.9 + ScaLAPACK 的笔记。

**关键论文**

- T. Ozaki, *Phys. Rev. B* 67, 155108 (2003). DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)
- A. Marek *et al.*, The ELPA library: scalable parallel eigenvalue solutions for electronic structure theory and computational science, *J. Phys.: Condens. Matter* 26, 213201 (2014). DOI: [10.1088/0953-8984/26/21/213201](https://doi.org/10.1088/0953-8984/26/21/213201)

---

[← 上一章](01-overview.md) | [返回目录](../README.md) | [下一章 →](03-input-first-calc.md)
