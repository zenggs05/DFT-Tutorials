# 第1章 CP2K 简介、版本与安装（Introduction, Versions & Installation of CP2K）

本章介绍 CP2K 的定位（Quickstep 模块、GPW/GAPW 方法，擅长大体系 DFT、从头算分子动力学 AIMD 与 QM/MM），梳理 2025.1 → 2026.2 的版本线和可执行文件命名（ssmp/psmp），并给出五种安装方式的完整命令：官方静态二进制、conda-forge、Docker/Podman/Apptainer 容器、toolchain + CMake 源码编译（CP2K installation, build from source）、Spack。最后讲解 `CP2K_DATA_DIR` 数据目录、`cp2k.psmp` 的命令行选项与 MPI×OpenMP 运行方式，并用一个 H2O 冒烟测试验证安装是否可用。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 CP2K 与 VASP/QE（平面波）、Gaussian（分子高斯基）的区别，知道什么问题适合用 CP2K；
- 读懂版本号（`2026.2`）、可执行文件后缀（`ssmp`/`psmp`/`popt`）和 `cp2k.psmp -v` 打印的 `cp2kflags`；
- 根据自己的机器（笔记本、工作站、无 root 权限的集群、Mac、Windows）选择合适的安装方式，并完成安装；
- 正确设置 `CP2K_DATA_DIR`，理解 CP2K 查找基组/赝势文件的顺序；
- 用 `cp2k.psmp -i x.inp -o x.out`、`mpirun`/`srun` 加 `OMP_NUM_THREADS` 运行 CP2K，并用冒烟测试判断安装是否真正可用（包括 MPI 是否配对正确）。

| 你的情况 | 推荐安装方式 | 本章小节 |
|---|---|---|
| 只想 10 分钟内跑通第一个算例（Linux x86_64） | 官方静态二进制 `cp2k-2026.2-Linux-gnu-x86_64.ssmp` | §4 |
| 个人 Linux 工作站 / WSL2，单节点 | conda-forge 的 `cp2k`（自带 MPI） | §5 |
| 集群有 Apptainer，或 Mac（Apple Silicon） | 官方容器 `cp2k/cp2k` | §6 |
| 集群正式生产计算，追求性能、多节点 | toolchain + CMake 源码编译，或用集群管理员提供的 module | §7 |
| 需要 GPU、DLA-Future、SIRIUS 等复杂依赖 | Spack（`make_cp2k.sh` 或 `spack install cp2k`） | §8 |
| Windows 笔记本教学演示 | WSL2 + conda/静态二进制；社区一键安装包仅限教学 | §9 |

> [!NOTE]
> 本教程以 **CP2K 2026.2**（2026-07-15 发布，截至 2026-11 的最新稳定版）为准。很多集群仍在用 2023.x–2025.x，凡是输入语法或构建方式有差异的地方，本系列都会同时给出新旧写法并标明版本。Linux 基础、conda、集群登录和文件传输见 [第1章（公共篇）计算环境搭建](../../Common/tutorials/01-linux-hpc.md)，作业调度系统见 [第2章（公共篇）HPC 作业调度](../../Common/tutorials/02-job-scheduler.md)。

---

## 2. CP2K 是什么：定位与核心方法

### 2.1 一句话定位

CP2K 是一个开源（GPL-2.0-or-later）的原子尺度模拟程序包，核心电子结构模块叫 **Quickstep**。Quickstep 用 **高斯基组表示 Kohn–Sham 轨道、用平面波（规则实空间网格 + FFT）表示电子密度和静电势**，这就是 GPW（Gaussian and Plane Waves）方法；它的全电子扩展是 GAPW（Gaussian and Augmented Plane Waves）。两种方法的原理在 [第3章](03-gpw-basis-pseudo.md) 详细展开。

在 GPW 中，Kohn–Sham 轨道写成原子中心收缩高斯函数的线性组合：

$$
\psi_i(\mathbf r)=\sum_{\mu} c_{\mu i}\,\varphi_\mu(\mathbf r),\qquad
n(\mathbf r)=\sum_{\mu\nu} P_{\mu\nu}\,\varphi_\mu(\mathbf r)\varphi_\nu(\mathbf r)
$$

密度 $n(\mathbf r)$ 被“投影”（collocate）到规则网格上，经 FFT 在倒空间求解 Poisson 方程得到 Hartree 势，再积分回高斯基。这样既保留了高斯基“紧凑、稀疏、适合大体系”的优点，又借用了平面波“求静电势快”的优点。

### 2.2 CP2K 擅长什么、不擅长什么

| 方面 | CP2K 的特点 | 对比 |
|---|---|---|
| 体系规模 | 几百到几千原子的凝聚相体系（液体、界面、溶液、MOF、蛋白活性中心）是强项；矩阵稀疏 + OT 方法，内存和时间随原子数增长平缓 | 平面波程序（VASP/QE）对同样大小的体系通常更慢、更吃内存 |
| 从头算分子动力学（AIMD） | 波函数外推（ASPC）、OT 快速收敛，适合几十 ps 的 Born–Oppenheimer MD | 见 [第12章](12-aimd.md) |
| 多尺度 | 原生 QM/MM、xTB/DFTB 半经验、经典力场（FIST）、机器学习势接口 | 见 [第16章](16-qmmm-xtb.md) |
| 增强采样 | 内置元动力学、约束 MD，可接 PLUMED | 见 [第15章](15-metadynamics-plumed.md) |
| 杂化泛函 | 实空间 HFX + 截断库仑 + ADMM，凝聚相杂化泛函 AIMD 可行 | 见 [第9章](09-hybrid-admm.md) |
| 小原胞金属、精细能带 | 可以做（k 点 + 对角化），但不是传统强项；2026.2 及以前 OT 不支持 k 点 | 小原胞金属、能带/费米面研究，平面波程序往往更方便 |
| 基组 | 高斯基组需要同时收敛“基组质量”和“网格截断”两个维度 | 平面波只需收敛 ENCUT |

> [!TIP]
> 判断是否该用 CP2K 的简单标准：**体系大（≥ 100 原子）、以 Gamma 点为主、要做 MD 或需要 QM/MM**，就优先考虑 CP2K；**原胞小、要密 k 点、要能带/费米面细节**，优先考虑平面波程序（本仓库的 [VASP 篇](../../VASP/README.md)）。

### 2.3 程序结构：FORCE_EVAL 的几种 METHOD

CP2K 的输入以“力的计算引擎”（`&FORCE_EVAL`）和“原子怎么动”（`&MOTION`）为两大支柱，详见 [第2章](02-input-structure.md)。`&FORCE_EVAL` 中 `METHOD` 的常用取值（以 2026.2 输入参考为准）：

| METHOD | 含义 | 典型用途 |
|---|---|---|
| `QS`（= `QUICKSTEP`） | Quickstep 电子结构（DFT、HF、xTB、DFTB、半经验） | 本教程绝大多数章节 |
| `FIST` | 经典分子力学 | QM/MM 的 MM 部分、力场 MD |
| `QMMM` | QM/MM 耦合 | 酶反应、溶液化学 |
| `SIRIUS` | 通过 SIRIUS 库做平面波 DFT | 与平面波结果对照 |
| `NNP` | 神经网络势 | 机器学习势 MD |
| `MIXED` | 多个 FORCE_EVAL 的组合 | 特殊耦合、能量差 |

---

## 3. 版本线与命名规则

### 3.1 发布节奏与近期版本

CP2K 每年发布两个版本，版本号为 `年份.序号`：`YYYY.1` 通常在 1 月，`YYYY.2` 通常在 7 月。近期版本的关键变化：

| 版本 | 发布日期 | 与教程相关的关键变化 |
|---|---|---|
| 2025.1 | 2025-01-02 | 输出总能量行改为 `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]`（2024.x 及以前为 `... energy [a.u.]:`） |
| 2025.2 | 2025-07-23 | **最后一个带 GNU Makefile/arch 文件的版本**；新增 GFN-xTB/Efield-xTB、TDDFPT 激子描述符等 |
| 2026.1 | 2026-01-06 | **移除 Makefile 和 arch 文件，只能用 CMake 构建**；DBCSR 变为必需的外部依赖；移除 QUIP；新增 `BASIS_AUG_MOLOPT`、`BASIS_RI_AUG_MOLOPT` 等 |
| **2026.2** | **2026-07-15** | **当前稳定版**。新增 `build_cp2k.sh`（toolchain 构建）和 `make_cp2k.sh`（Spack 构建）两个脚本；删除 `&MOTION/&CELL_OPT/TYPE`；DOS/PDOS 输入重构为 `&PRINT/&DOS`；REFTRAJ 改用 `EVAL`；删除 Cython Python 绑定；不再支持 GCC 8；数据目录新增带版本号的 UZH 基组/赝势快照 |
| trunk（master） | 持续更新 | 已有 OT + k 点、OT + 展宽等 2026.2 之后的新功能，预计进入下一版本 |

2026.2 的破坏性变更会让旧输入直接报错，涉及的章节如下（迁移清单见 [第18章](18-troubleshooting.md)）：

| 变更 | 旧写法的报错 | 处理 | 详见 |
|---|---|---|---|
| 删除 `&CELL_OPT` 中的 `TYPE` | `found an unknown keyword TYPE in section CELL_OPT` | 删掉 `TYPE DIRECT_CELL_OPT` 一行 | [第7章](07-geo-cell-opt.md) |
| `&DFT/&PRINT/&PDOS` 改为 `&DFT/&PRINT/&DOS/&PDOS` | `unknown subsection PDOS of section PRINT` | 改用新的 `&DOS` 段 | [第10章](10-dos-pdos-band.md) |
| REFTRAJ 删除 `EVAL_ENERGY_FORCES`/`EVAL_FORCES` | 未知关键字 | 改用 `&REFTRAJ` 中的 `EVAL` | [第12章](12-aimd.md) |

> [!IMPORTANT]
> 写论文时请在方法部分注明 CP2K 的**确切版本号**（如 “CP2K 2026.2”）以及所用基组/赝势文件名。不同版本的默认值、基组文件内容可能不同（[第3章](03-gpw-basis-pseudo.md) 会看到 2026.2 的 `BASIS_MOLOPT_UZH` 已更新了部分元素）。各版本对应的手册入口见 <https://manual.cp2k.org/trunk/versions.html>，2026.2 手册为 <https://manual.cp2k.org/cp2k-2026_2-branch/index.html>。

### 3.2 可执行文件后缀

CP2K 可执行文件名形如 `cp2k.<后缀>`，后缀由“串行/并行”和“优化/调试”两部分组成：

| 可执行文件 | 含义 | 用法 |
|---|---|---|
| `cp2k.ssmp` | serial + OpenMP：单进程多线程，无 MPI | 笔记本、单节点小算例、静态二进制 |
| `cp2k.psmp` | parallel + OpenMP：MPI + OpenMP 混合 | **集群生产计算的标准选择** |
| `cp2k.sdbg` / `cp2k.pdbg` | 对应的调试版（Debug 构建） | 开发和查 bug，速度慢 |
| `cp2k.sopt` / `cp2k.popt` | 指向 `ssmp`/`psmp` 的符号链接 | 以这个名字启动时，程序**强制每进程 1 个 OpenMP 线程**（相当于 `OMP_NUM_THREADS=1`） |
| `cp2k` | 部分发行版（conda-forge、容器）提供的链接，指向 `cp2k.psmp` | 与 `cp2k.psmp` 相同 |

> [!NOTE]
> 旧教程和旧作业脚本里常见 `cp2k.popt`（纯 MPI）。在新版本里它只是 `cp2k.psmp` 的别名加上“线程数固定为 1”，所以 `mpirun -np 28 cp2k.popt` 与 `OMP_NUM_THREADS=1 mpirun -np 28 cp2k.psmp` 等价。

### 3.3 `cp2k.psmp -v`：查看版本与编译特性

```bash
cp2k.psmp -v
```

官方构建文档给出的一个示例输出（省略了编译选项部分）：

```text
 CP2K version 2026.2
 Source code revision git:c92cc08
 cp2kflags: omp libint fftw3 libxc elpa parallel scalapack mpi_f08 cosma libxs spglib openblas libdftd4 s_dftd3 mctc-lib tblite libvori libbqb
 compiler: GCC version 14.3.1 20251022 (Red Hat 14.3.1-4)
 compiler target: cpuid   1002 (x86_avx2)
```

`cp2kflags` 告诉你这个二进制“能做什么”。常见标志：

| 标志 | 含义 | 缺少时的后果 |
|---|---|---|
| `omp` | 启用 OpenMP | 只能单线程 |
| `parallel`、`mpi_f08`、`scalapack` | MPI 并行版本（psmp） | 只能单进程 |
| `libint` | libint2 四中心积分 | **不能做杂化泛函/HFX、MP2/RPA**（[第9章](09-hybrid-admm.md)） |
| `libxc` | libxc 泛函库 | 只能用 CP2K 内置的少数泛函（PBE、BLYP 等内置实现仍可用） |
| `fftw3` | FFTW3 快速傅里叶变换 | 性能差；2026.2 说明未来可能把 FFTW3 接口变成硬依赖 |
| `elpa` | ELPA 本征值求解器 | 大体系对角化变慢（可用 ScaLAPACK） |
| `cosma` | COSMA 矩阵乘法 | 部分大规模计算变慢 |
| `spglib` | 空间群对称性库 | 不能用 `KEEP_SPACE_GROUP` 等依赖对称性的功能（[第7章](07-geo-cell-opt.md)） |
| `libdftd4`、`s_dftd3`、`mctc-lib` | DFT-D4 / simple-DFT-D3 外部库 | 无法用外部库版本的 D4（[第8章](08-dispersion-spin-dftu.md)） |
| `tblite` | tblite 库（GFN2-xTB） | 只能用原生 GFN0/GFN1（[第16章](16-qmmm-xtb.md)） |
| `libvori`、`libbqb` | Voronoi 积分、BQB 轨迹压缩 | 相应分析不可用 |

PLUMED、DeePMD、GPU（CUDA/HIP）等可选依赖若已编译进去，也会以各自的标志出现在这一行。**拿到一个新的 CP2K（尤其是集群 module）时，先看 `-v` 输出，再决定能否跑你的任务。**

---

## 4. 方式一：官方静态二进制（最快上手）

每个 CP2K Release 页面都附带一个静态链接的 `ssmp` 二进制，2026.2 的文件名是 `cp2k-2026.2-Linux-gnu-x86_64.ssmp`。它不依赖系统上的任何库，下载、赋予执行权限即可运行。

### 4.1 下载与运行

```bash
mkdir -p ~/opt/cp2k-2026.2/bin && cd ~/opt/cp2k-2026.2

# 1) 下载静态二进制（约几十到上百 MB）
wget https://github.com/cp2k/cp2k/releases/download/v2026.2/cp2k-2026.2-Linux-gnu-x86_64.ssmp \
     -O bin/cp2k-2026.2-Linux-gnu-x86_64.ssmp
chmod +x bin/cp2k-2026.2-Linux-gnu-x86_64.ssmp

# 2) 建一个短名字的链接，方便后续脚本统一调用 cp2k.ssmp
ln -s cp2k-2026.2-Linux-gnu-x86_64.ssmp bin/cp2k.ssmp
export PATH=$HOME/opt/cp2k-2026.2/bin:$PATH

# 3) 检查版本和编译特性
cp2k.ssmp -v
```

### 4.2 准备数据目录（基组、赝势）

静态二进制只有一个可执行文件，**不附带 `data/` 目录**（基组、赝势、DFT-D3 参数等）。最稳妥的做法是从同版本源码包中只解压 `data/`，并用 `CP2K_DATA_DIR` 指向它：

```bash
cd ~/opt/cp2k-2026.2
wget https://github.com/cp2k/cp2k/releases/download/v2026.2/cp2k-2026.2.tar.bz2
sha256sum cp2k-2026.2.tar.bz2          # 与 Release 页面公布的 SHA256 对比，一致才使用
tar -xjf cp2k-2026.2.tar.bz2 cp2k-2026.2/data   # 只解压 data/ 子目录
mv cp2k-2026.2/data ./data && rmdir cp2k-2026.2

export CP2K_DATA_DIR=$HOME/opt/cp2k-2026.2/data
ls $CP2K_DATA_DIR | head
```

把两条 `export` 写进 `~/.bashrc` 或本章示例中的 [`env_cp2k.sh`](../examples/01-install-check/env_cp2k.sh)。

> [!WARNING]
> 静态二进制是 **ssmp（无 MPI）**，只能在单节点上用 OpenMP 线程并行，且为了兼容性采用保守的编译选项，性能不如针对本机优化的源码编译版本。它适合学习、调试输入和小体系计算；多节点生产计算请用 §7/§8 或集群 module。

### 4.3 适用范围小结

| 优点 | 局限 |
|---|---|
| 无需编译、无需 root，几分钟可用 | 仅 Linux x86_64；只有 OpenMP 并行 |
| 与官方测试环境一致，适合复现问题 | 不能跨节点；性能非最优 |
| 可作为“参考答案”核对集群上自编译版本的结果 | 可选库（如 PLUMED）是否包含要看 `-v` |

---

## 5. 方式二：conda-forge

conda-forge 提供带 MPI 的 CP2K 包，适合个人工作站和 WSL2。

### 5.1 安装

```bash
# 列出可用版本和构建（build string 中包含 mpi_mpich / mpi_openmpi）
conda search conda-forge::cp2k

# 推荐：装进独立环境，避免与其他软件的依赖冲突
conda create -n cp2k_env conda-forge::cp2k
conda activate cp2k_env
cp2k -h -v
```

conda-forge 的 CP2K 2026.2 同时有 MPICH 和 OpenMPI 两种构建（求解器默认优先选择 MPICH 构建）。需要指定某一种时，用构建字符串通配：

```bash
conda create -n cp2k_ompi "conda-forge::cp2k=2026.2=mpi_openmpi_*"
```

安装后环境中有 `cp2k.psmp`、`cp2k.popt` 和指向 `cp2k.psmp` 的 `cp2k`，数据文件位于 `$CONDA_PREFIX/share/cp2k/data`。建议把它显式写进环境变量（`conda env config vars set` 会在每次 `conda activate` 时自动设置）：

```bash
conda activate cp2k_env
conda env config vars set CP2K_DATA_DIR=$CONDA_PREFIX/share/cp2k/data
conda deactivate && conda activate cp2k_env
echo $CP2K_DATA_DIR
```

### 5.2 运行

conda 环境自带 MPI，**必须用同一环境里的 `mpirun`/`mpiexec`**：

```bash
conda activate cp2k_env
which mpirun cp2k.psmp          # 两者都应在 $CONDA_PREFIX/bin 下
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out   # OpenMPI 构建
mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out # MPICH 构建
```

> [!WARNING]
> - conda-forge 的 2026.2 只有 **linux-64 和 linux-aarch64**；**osx-64 停留在 7.1.0**，Mac 用户用 conda 会装到一个非常老的版本。Mac 请用 §6 的容器或 Homebrew（`brew install cp2k`）。
> - 在集群上用 conda 版 MPI 跨节点运行，通常无法利用 InfiniBand 等高速网络的最佳配置，性能可能明显下降。conda 版适合单节点；多节点生产计算优先用集群 module 或源码编译。
> - 不要在已经 `module load` 了其他 MPI 的 shell 里混用 conda 的 CP2K，否则极易出现 §11 中的 “MPI 错配”。

---

## 6. 方式三：容器（Docker / Podman / Apptainer）

官方在 Docker Hub 发布 `cp2k/cp2k` 镜像。容器内已包含 CP2K、MPI 和完整的 `data/` 目录，**不需要设置 `CP2K_DATA_DIR`**。镜像 tag 会标明版本、MPI 实现（`mpich`/`openmpi`）、目标 CPU、是否 CUDA 以及 `psmp`，例如 `2026.2_mpich_aarch64_psmp`；不确定时用 `latest`（最新发布版，MPICH，通用 x86_64）。完整 tag 列表见 <https://hub.docker.com/r/cp2k/cp2k/tags>。

### 6.1 Podman / Docker（工作站、Mac）

官方推荐 Podman（无守护进程，默认以普通用户身份运行，输出文件归属当前用户）：

```bash
podman pull docker.io/cp2k/cp2k:latest
podman run --rm docker.io/cp2k/cp2k:latest cp2k -h -v

# 在当前目录运行 h2o.inp：当前目录挂载到容器内 /mnt（容器的工作目录）
cd ~/work/h2o
podman run --rm --shm-size=1g -v "$PWD:/mnt" docker.io/cp2k/cp2k:latest \
    mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k -i h2o.inp -o h2o.out
```

用 Docker 时写法相同，建议加 `-u $(id -u):$(id -g)`，否则输出文件归 root 所有：

```bash
docker pull cp2k/cp2k:latest
docker run -it --rm --shm-size=1g -v $PWD:/mnt -u $(id -u):$(id -g) cp2k/cp2k:latest \
    mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k -i h2o.inp -o h2o.out
```

| 参数 | 作用 |
|---|---|
| `--rm` | 运行结束后删除容器 |
| `--shm-size=1g` | 增大共享内存，MPI 进程间通信需要；太小会报 `/dev/shm` 相关错误 |
| `-v "$PWD:/mnt"` | 把当前目录挂到容器的 `/mnt`；输入、输出都在宿主机当前目录 |
| `-genv OMP_NUM_THREADS=1` | MPICH 传环境变量的写法；OpenMPI 镜像用 `-x OMP_NUM_THREADS=1` |
| `-it` | 允许 `Ctrl-C` 中断计算 |

**Mac（Apple Silicon）**：先 `brew install podman`，再 `podman machine init --cpus <N> --memory <MiB>` 和 `podman machine start` 创建 Linux 虚拟机，然后拉取 aarch64 镜像：

```bash
podman pull docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp
cd ~/cp2k-work/h2o        # 必须位于 /Users 之下，Podman machine 才能看到
podman run --rm --shm-size=1g -v "$PWD:/mnt" docker.io/cp2k/cp2k:2026.2_mpich_aarch64_psmp \
    mpiexec -n 4 -genv OMP_NUM_THREADS=2 cp2k -i h2o.inp -o h2o.out
```

aarch64 镜像在 Apple Silicon 上原生运行，无需指令集模拟。完整步骤见官方 macOS 指南 <https://manual.cp2k.org/trunk/getting-started/macos.html>。

### 6.2 Apptainer（集群首选的容器方式）

大多数集群不允许 Docker，但提供 Apptainer（原 Singularity）：

```bash
apptainer pull docker://cp2k/cp2k:latest          # 生成 cp2k_latest.sif
apptainer run cp2k_latest.sif cp2k -h -v

# 单节点：用容器内的 MPI
apptainer run -B $PWD cp2k_latest.sif \
    mpiexec -n 4 -genv OMP_NUM_THREADS=2 cp2k -i h2o.inp -o h2o.out

# 多节点（Slurm）：用宿主机的 srun 启动，每个进程运行容器里的 cp2k
srun -n 4 apptainer run -B $PWD cp2k_latest.sif cp2k -i h2o.inp -o h2o.out
```

跨节点时“宿主机 MPI + 容器 MPI”必须兼容（如都是 MPICH ABI），否则可能无法启动；官方说明这种方式性能最好，但不保证在所有集群上可用。

### 6.3 容器自检：run_tests

镜像内置完整回归测试：

```bash
# Apptainer：默认 8 个任务，可用 --maxtasks 调整（最好是 4 的倍数）
apptainer run cp2k_latest.sif run_tests --maxtasks 8 --workbasedir .
# Podman / Docker
podman run -it --rm --shm-size=1g docker.io/cp2k/cp2k:latest run_tests
```

完整测试需要较长时间；平时只需跑 §11 的 H2O 冒烟测试。

> [!NOTE]
> 从 2026.2 起，构建 CP2K 容器的 Dockerfile 已并入主仓库的 `tools/docker/`（如 `Dockerfile.test_spack_psmp-gcc15`）；原 `cp2k-containers` 仓库主要保留 Apptainer/Docker 的使用说明：<https://github.com/cp2k/cp2k-containers>。

---

## 7. 方式四：源码编译（toolchain + CMake）

源码编译能针对本机 CPU 优化、链接集群的 MPI 与数学库，是集群生产计算的首选。**2026.1 起只支持 CMake**；2026.2 新增 `tools/toolchain/build_cp2k.sh`，可以根据 toolchain 的配置自动生成 CMake 选项并完成编译安装。

### 7.1 依赖与前提

| 项目 | 要求（2026.2 文档） |
|---|---|
| 编译器 | 支持 Fortran 2008（`-std=f2008`）与 C11 的 Fortran/C 编译器；**2026.2 不再支持 GCC 8** |
| 必需库 | DBCSR、BLAS、LAPACK、FFTW、libwignernj ≥ 0.8 |
| MPI 版本额外需要 | MPI、ScaLAPACK |
| 常用可选库 | libint（杂化泛函）、libxc、ELPA、COSMA、spglib、libvori、DFT-D4、tblite、PLUMED、SIRIUS 等 |
| 系统工具 | `python3`、`git`、`make`、`cmake`、`wget`、`tar`、`bzip2`、`xz` 等 |

> [!TIP]
> toolchain 默认会自行下载并编译这些依赖，**不需要 root 权限**。离线集群可以在有网的机器上先下载好压缩包，放进 `tools/toolchain/build/` 目录；toolchain 发现文件名和 sha256 匹配的包时不会再去下载。

### 7.2 获取源码

**优先使用 Release 页面上的版本压缩包 `cp2k-2026.2.tar.bz2`**，不要用 GitHub 自动生成的 “Source code” 压缩包（后者对 ≤ 2025.2 的版本缺少捆绑的子模块）：

```bash
cd ~/src
wget https://github.com/cp2k/cp2k/releases/download/v2026.2/cp2k-2026.2.tar.bz2
sha256sum cp2k-2026.2.tar.bz2          # 与 Release 页面的 SHA256 一致才继续
tar -xjf cp2k-2026.2.tar.bz2
cd cp2k-2026.2                           # 以下称 CP2K_ROOT
```

> [!WARNING]
> 官方文档提示：在带 tag 的 Git 克隆目录里构建 toolchain 后再用 CMake 配置，DBCSR 的版本检测可能误读外层 CP2K 仓库信息，导致 CMake 找 DBCSR 时报错（GitHub issue #5184）。**发布版请用版本压缩包而不是 git tag。** 开发版才用 `git clone https://github.com/cp2k/cp2k.git`。

### 7.3 第一步：install_cp2k_toolchain.sh 安装依赖

toolchain 脚本**必须在 `tools/toolchain/` 目录内执行**（在其他目录调用会报 `Incorrect execution location`）：

```bash
cd tools/toolchain
./install_cp2k_toolchain.sh --help        # 完整选项
./install_cp2k_toolchain.sh --help-hpc    # 面向集群用户的额外提示
```

按场景选择一种配置（`-j` 为并行编译进程数）：

```bash
# A) 纯净环境：依赖全部由 toolchain 下载编译（包括 GCC 和 MPI），最省心但最慢
./install_cp2k_toolchain.sh -j 16 --install-all

# B) 集群：复用系统/模块提供的 GCC 和 OpenMPI，其余依赖由 toolchain 编译
module load gcc/13 openmpi/4.1           # 模块名以你的集群为准
./install_cp2k_toolchain.sh -j 16 --with-gcc=system --with-openmpi=system --mpi-mode=openmpi

# C) 登录节点与计算节点 CPU 不同：放弃 native 优化，换取可移植性
./install_cp2k_toolchain.sh -j 16 --with-gcc=system --with-openmpi=system --target-cpu=generic

# D) NVIDIA GPU（以 A100 为例）
./install_cp2k_toolchain.sh -j 16 --enable-cuda --gpu-ver=A100

# 只想看最终生效的配置、不真正下载编译：
./install_cp2k_toolchain.sh --with-gcc=system --with-openmpi=system --dry-run
```

| 选项 | 含义 |
|---|---|
| `--with-PKG=install` | 由 toolchain 下载并编译该包（单写 `--with-PKG` 等价于 `=install`） |
| `--with-PKG=system` | 从系统路径（含已 `module load` 的模块）寻找该包 |
| `--with-PKG=<path>` | 使用指定路径下已安装的包 |
| `--with-PKG=no` | 不使用该包 |
| `--install-all` | 把所有 `--with-PKG` 设为 `install`（其后的选项可以再单独覆盖） |
| `--mpi-mode=mpich/openmpi/intelmpi/no` | MPI 类型；`no` 表示不要 MPI（只编译 ssmp） |
| `--math-mode=mkl/openblas/...` | 核心数学库；存在 `$MKLROOT` 时默认用 MKL，否则 OpenBLAS |
| `--target-cpu=native/generic/haswell/...` | 编译优化目标，默认 `native` |
| `--enable-cuda` + `--gpu-ver=A100/H100/...` | 开启 CUDA 并指定 GPU 架构 |
| `--libint-lmax=4/5/6/7` | libint 支持的最高角动量，默认 5；越高编译越久 |
| `--no-check-certificate` | wget 证书校验失败时的权宜之计（每个包仍会做 sha256 校验） |

依赖默认编译到 `tools/toolchain/build/`，安装到 `tools/toolchain/install/`（可用 `--install-dir=/abs/path` 改到源码树外）。

> [!CAUTION]
> `--target-cpu=native` 会针对**执行 toolchain 的那台机器**优化。如果在新 CPU 的登录节点上编译、在老 CPU 的计算节点上运行，程序会以 `Illegal instruction` 崩溃。解决：在计算节点上（交互作业）编译，或使用 `--target-cpu=generic`/与计算节点匹配的架构名。

### 7.4 第二步：build_cp2k.sh 编译安装（2026.2 新增）

toolchain 结束时会提示运行 `build_cp2k.sh`。在 **2026.2** 中：

```bash
# 仍在 tools/toolchain/ 目录下
./build_cp2k.sh --help
./build_cp2k.sh -j 16                                  # 安装到 <CP2K_ROOT>/install
./build_cp2k.sh -j 16 --prefix /opt/cp2k/2026.2        # 2026.2：--prefix 后跟绝对路径（空格分隔）
./build_cp2k.sh --dry-run                               # 只打印生成的 CMake 选项
```

编译完成后，`<CP2K_ROOT>/build/` 是 CMake 构建目录（可删除以节省空间），`<CP2K_ROOT>/install/` 是安装目录，里面有 `bin/cp2k.psmp`、库文件、`share/` 下的数据文件，以及**运行前必须 source 的环境文件 `cp2k_env`**：

```bash
source ~/src/cp2k-2026.2/install/cp2k_env     # 每个新 shell、每个作业脚本都要执行
which cp2k.psmp
cp2k.psmp -v
```

> [!IMPORTANT]
> **版本差异（`build_cp2k.sh` 的选项）**：
> - **2026.2**：选项为 `-j N`、`--prefix /abs/path`、`--dry-run`、`--debug`、`--build-static`、`--rebuild-only`；**没有 `--preset`**，CPU 优化目标由 toolchain 的 `--target-cpu` 决定。
> - **trunk/master（2026.2 之后）**：改为 `--prefix=/abs/path` 写法，并新增 `--preset=<name>`（如 `--preset=native-gnu-x86_64`，可用 `cmake -S <CP2K_ROOT> --list-presets` 查看）。网上教程里的 `build_cp2k.sh --preset native-gnu-x86_64` 指的是开发版。
> - **2026.1**：还没有 `build_cp2k.sh`，需要按 toolchain 最后打印的提示手动执行 `cmake`（见 §7.6）。
>
> 以你手上版本的 `./build_cp2k.sh --help` 为准。

### 7.5 第三步：回归测试

```bash
source ~/src/cp2k-2026.2/install/cp2k_env
cd ~/src/cp2k-2026.2
./tests/do_regtest.py --help
# 冒烟测试（只跑少量代表性算例，几分钟）
./tests/do_regtest.py install/bin psmp --maxtasks 8 --smoketest
# 完整回归测试（数千个输入，耗时数小时）
./tests/do_regtest.py install/bin psmp --maxtasks 32
```

测试结束会汇总 `correct / wrong / failed` 数目。少量 `wrong`（数值超出容差）在非官方测试平台上并不罕见，但大量 `failed` 通常说明环境或依赖有问题。

### 7.6 手动 CMake（自管依赖，或 2026.1）

依赖已经由 module/Spack/自己准备好时，可以直接调用 CMake。下面是一个启用 MPI 和常用库的示例（选项名以 2026.2 文档为准，按需增减）：

```bash
cd ~/src/cp2k-2026.2
cmake -S . -B build -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k/2026.2 \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_COSMA=ON
cmake --build build --target install -j 32
```

| 选项 | 说明 |
|---|---|
| `-S . -B build` | 源码树根目录与构建目录（禁止在源码目录内原地构建） |
| `-GNinja` | 用 Ninja 生成器；GPU 构建强烈建议 |
| `--preset native-gnu-x86_64` | CMake preset：对 GNU 编译器加 `-march=native`；`cmake -S . --list-presets` 列出全部，`none` 表示不加额外优化 |
| `-DCMAKE_BUILD_TYPE=Release/Debug` | Debug 会生成 `psmp`/`ssmp` 的调试版 `pdbg`/`sdbg` |
| `-DCP2K_DATA_DIR=/abs/path` | 编译进程序的数据目录，默认是安装前缀下的 `share/cp2k/data` |
| `-DCP2K_USE_ACCEL=CUDA` + `-DCP2K_WITH_GPU=A100` | CUDA 加速及 GPU 型号 |
| `-DCP2K_USE_EVERYTHING=ON` | 打开全部可选依赖（需要它们都已就位） |

> [!NOTE]
> **≤ 2025.2 的旧构建方式**（Makefile + arch 文件）：toolchain 结束后执行 `cp install/arch/local.* ../../arch/`、`source install/setup`，回到源码根目录 `make -j 32 ARCH=local VERSION="ssmp psmp"`，可执行文件在 `exe/local/`。2026.1 起该方式已被删除，旧集群上的 arch 文件不能用于新版本。

---

## 8. 方式五：Spack（make_cp2k.sh 与 spack install）

Spack 能自动解决复杂的依赖链（GPU、DLA-Future、PEXSI、SIRIUS 等 toolchain 覆盖不全的组合），官方文档也提示 toolchain 未来“可能被精简或退役”，转而以 Spack 为主。

### 8.1 make_cp2k.sh（2026.2 新增，在 CP2K_ROOT 执行）

`make_cp2k.sh` 在源码根目录下建立一个本地 Spack（`spack/` 子目录），装好依赖后再用 CMake 编译 CP2K：

```bash
cd ~/src/cp2k-2026.2
./make_cp2k.sh --help
./make_cp2k.sh -j 32                                   # 默认：几乎全部功能，psmp，MPICH
./make_cp2k.sh -j 32 -mpi openmpi                      # 改用 OpenMPI
./make_cp2k.sh -j 32 -gm A100                          # NVIDIA A100 GPU
./make_cp2k.sh -j 32 -bd -df all -ef libint -ef spglib -ef elpa   # 最小构建再按需加功能
./make_cp2k.sh -j 32 -t "--maxtasks 16 --smoketest"    # 构建后立即跑回归测试
```

| 选项（2026.2，脚本版本 2.0） | 含义 |
|---|---|
| `-bd`/`--build_deps` | 删除 `spack/` 目录，从头重建全部依赖 |
| `-cv psmp/ssmp/...` | 要构建的 CP2K 版本，默认 `psmp` |
| `-df`/`-ef <包或功能>` | 关闭/开启某个包（如 `libint`、`elpa`、`plumed`、`tblite`）或功能 |
| `-gm <型号>` | GPU 型号（P100、V100、A100、H100、H200、GH200 等） |
| `-mpi mpich/openmpi/no` | MPI 类型，默认 `mpich` |
| `-ip <路径>` | 安装路径，默认 `./install` |
| `-t "TESTOPTS"` | 构建成功后运行 `tests/do_regtest.py` 并传入这些参数 |
| `-ue` | 使用宿主系统已有的软件包（更快，但可能拉进过旧的 gcc/python） |

> [!NOTE]
> 开发版（trunk）的 `make_cp2k.sh` 还增加了 `-preset`、`-ase`、`-gromacs` 等选项，包名也有变化（如 trunk 写 `libint2`，2026.2 写 `libint`）。**一律以 `./make_cp2k.sh --help` 的输出为准。**

用 Spack 构建的 CP2K 通过 `install/bin/launch` 启动（它会设置 `PATH`、`LD_LIBRARY_PATH`，并在未设置时默认 `OMP_NUM_THREADS=2`），回归测试脚本是 `install/bin/run_tests`：

```bash
~/src/cp2k-2026.2/install/bin/launch cp2k.psmp -v
~/src/cp2k-2026.2/install/bin/launch mpirun -np 4 cp2k.psmp -i h2o.inp -o h2o.out
~/src/cp2k-2026.2/install/bin/run_tests
```

### 8.2 直接使用 spack install

已有 Spack 环境（很多超算中心维护了 Spack）时：

```bash
spack info cp2k                          # 查看可用版本和变体（variants）
spack install cp2k +libint +libxc
spack install cp2k +libint +libxc +cuda cuda_arch=80      # A100 = compute capability 8.0
spack install cp2k ^intel-oneapi-mkl +cluster             # 用 oneAPI MKL 作为数学库
spack load cp2k                          # 装了多个变体时要写具体，如 spack load cp2k +cuda
cp2k.psmp -v
```

`+变体` 开启、`~变体` 关闭，`^依赖` 约束依赖包，`@版本` 指定版本（如 `cp2k@2026.2`，前提是你的 Spack 仓库已收录该版本，用 `spack info cp2k` 确认）。

---

## 9. Windows 与 macOS

| 平台 | 推荐做法 | 说明 |
|---|---|---|
| Windows 10/11 | **WSL2 + Ubuntu**，然后在 WSL 里用 §4 静态二进制或 §5 conda | `wsl --install -d Ubuntu`（管理员 PowerShell）；计算文件放在 WSL 的 Linux 文件系统（如 `~/work`）中，比放在 `/mnt/c` 下快得多 |
| Windows | Docker Desktop / Podman Desktop + §6 容器 | 适合已有容器环境的用户 |
| Windows | 社区一键安装包（如 <https://github.com/Shxnek/cp2k-windows-installer>） | **仅限教学演示**；正式计算用 Linux、WSL 或容器 |
| macOS（Apple Silicon） | Podman + `cp2k/cp2k:2026.2_mpich_aarch64_psmp` 镜像，或 `brew install cp2k` | 不要用 conda（osx-64 停在 7.1.0） |

> [!WARNING]
> 在 Windows 记事本里编辑过的作业脚本带有 CRLF 换行，提交到 Linux 后会报 `/bin/bash^M: bad interpreter`。用 `dos2unix job.sh` 或在 VS Code 右下角把换行改成 LF。输入文件也建议统一为 LF、纯 ASCII。

---

## 10. CP2K_DATA_DIR 与 data/ 目录

### 10.1 CP2K 怎样找到基组和赝势文件

输入里 `BASIS_SET_FILE_NAME BASIS_MOLOPT` 只写了文件名。对于需要读取的已有文件（基组、赝势、`dftd3.dat`、`@INCLUDE` 的文件等），CP2K 按以下顺序查找（2026.2 源码 `cp_files.F`/`cp_data_dir.c`）：

1. **按所写名字直接打开**：绝对路径，或相对于**当前工作目录**的路径；
2. 若不存在，再到**数据目录**下找：环境变量 `CP2K_DATA_DIR` 一旦设置就覆盖编译时写入的默认路径（CMake 选项 `-DCP2K_DATA_DIR`，默认为安装前缀下的 `share/cp2k/data`）。

两处都找不到时报错：

```text
The specified OLD file <BASIS_MOLOPT_UZH> cannot be opened. It does not exist. Data directory path: /path/that/cp2k/used
```

报错末尾的 `Data directory path` 就是 CP2K 实际使用的数据目录。正常运行时，输出文件头部也会打印这一行：

```bash
grep "CP2K| Data directory path" h2o.out
```

> [!TIP]
> 由第 1 步可知：把自定义的基组文件（例如从文献补充材料拿到的 `MY_BASIS`）**直接放在计算目录**，输入写 `BASIS_SET_FILE_NAME MY_BASIS` 即可，无需改数据目录。也可以写绝对路径。

### 10.2 data/ 目录中的关键文件

| 文件 | 内容 | 用在哪 |
|---|---|---|
| `BASIS_MOLOPT` | 经典 MOLOPT 基组：`SZV/DZVP/TZVP/TZV2P(X)-MOLOPT-GTH`、`SZV/DZVP-MOLOPT-SR-GTH` | 与 `GTH_POTENTIALS` 配对，大量文献使用 |
| `BASIS_MOLOPT_UZH` | UZH 协议基组：`DZVP-MOLOPT-PBE-GTH-q6`（别名 `DZVP-MOLOPT-GGA-GTH-q6`）等，另含全电子 `SVP/TZVPP/QZVPP-MOLOPT-PBE-ae` | **新项目推荐**，与 `POTENTIAL_UZH` 成对使用（[第3章](03-gpw-basis-pseudo.md)） |
| `BASIS_MOLOPT_UZH_2026.1`、`BASIS_MOLOPT_UZH_2026.2` | 2026.2 新增的带版本号快照；不带版本号的文件指向最新快照 | 保证结果可复现 |
| `POTENTIAL_UZH`（及 `_2026.1`/`_2026.2` 快照） | UZH 协议 GTH 赝势：`GTH-PBE-qN`、`GTH-SCAN-qN`、`GTH-PBE0-qN`，以及全电子 `ALL` | 与 UZH 基组配对 |
| `GTH_POTENTIALS` | 经典 GTH 赝势（`GTH-PBE-q6` 等） | 与 `BASIS_MOLOPT` 配对 |
| `POTENTIAL` | 更大的赝势合集（含多种泛函的 GTH 与 `ALL`） | `POTENTIAL_FILE_NAME` 的默认值 |
| `BASIS_SET` | 老的 GTH 基组合集 | `BASIS_SET_FILE_NAME` 的默认值 |
| `BASIS_ADMM_UZH`、`BASIS_ADMM` | ADMM 辅助基组（`admm-dzp-q6`、`cFIT3` 等） | 杂化泛函（[第9章](09-hybrid-admm.md)） |
| `EMSL_BASIS_SETS`、`ALL_BASIS_SETS` | Pople（`6-31Gx` 即 6-31G*）、Dunning 等全电子基组 | GAPW 全电子计算 |
| `dftd3.dat` | DFT-D3 参数 | `&VDW_POTENTIAL`（[第8章](08-dispersion-spin-dftu.md)） |
| `t_c_g.dat` | 截断库仑势表 | 周期性杂化泛函（[第9章](09-hybrid-admm.md)） |
| `rVV10_kernel_table.dat`、`vdW_kernel_table.dat` | 非局域 vdW 泛函核函数表 | rVV10 / vdW-DF |
| `xTB0_parameters`、`xTB1_parameters` | 原生 GFN0/GFN1-xTB 参数 | [第16章](16-qmmm-xtb.md) |
| `xc_section/` | 预定义 `&XC` 片段（如 `PBE.sec`），配合预处理器 `@XCTYPE` 使用 | [第2章](02-input-structure.md) |

> [!NOTE]
> 本仓库 `CP2K/data/` 下存有一份数据文件副本（其中 UZH 文件与 2025.2/2026.1 发布版一致，说明见 [CP2K/data/README.md](../data/README.md)），方便离线查阅基组名称。**实际计算请使用与所运行 CP2K 版本配套的 data 目录**，不要把不同版本的数据文件混用。

查某个元素有哪些基组/赝势，最直接的方法是 `grep`：

```bash
grep -E "^ *Ti " $CP2K_DATA_DIR/BASIS_MOLOPT_UZH      # Ti 的 UZH 基组
grep -E "^Cu " $CP2K_DATA_DIR/POTENTIAL_UZH            # Cu 的赝势有 q1/q11/q19 三种
```

---

## 11. 运行 CP2K：命令行与并行

### 11.1 命令行选项

以下选项取自 2026.2 源码中的帮助信息（`cp2k.psmp -h`）：

| 选项 | 作用 |
|---|---|
| `-i <file>` | 输入文件；若输入文件是最后一个参数，`-i` 可省略 |
| `-o <file>` | 输出文件（默认输出到屏幕）。**注意：若文件已存在，新的输出会追加到末尾，不会覆盖** |
| `-c`, `--check` | 只做输入语法检查，不运行。成功时打印 `SUCCESS, the input could be parsed correctly.` |
| `-e`, `--echo` | 回显经过预处理、**所有默认值显式写出**的完整输入，不运行（查默认值的利器） |
| `-E name=value`, `--set name=value` | 设置预处理器变量的初值（[第2章](02-input-structure.md)） |
| `-s`, `--shell` | 交互式 shell 模式（ASE 等外部程序驱动 CP2K 时使用，见 [第17章](17-tools-ecosystem.md)） |
| `--xml` | 把整个输入结构导出为 `cp2k_input.xml`（供 cp2k-input-tools、pycp2k 等工具使用） |
| `-v`, `--version` | 打印版本、源码版本号、`cp2kflags` 与编译器信息 |
| `-h`, `--help` | 打印帮助 |
| `-r`, `--run` | 与 `-c` 等同时使用时仍强制运行 |
| `--license` | 打印许可证 |
| `--mpi-mapping <method>` | 对 MPI 进程做重排映射（高级） |

> [!WARNING]
> `-o` 是**追加写入**。同一个目录里重复运行 `cp2k.psmp -i x.inp -o x.out`，旧输出会留在文件前半部分，`grep` 时容易拿到上一轮的结果。重跑前先 `rm -f x.out` 或改名；作业脚本里也建议这样做（本章示例脚本已处理）。

想同时在屏幕和文件里看到输出，用 `tee`（这时不要再加 `-o`）：

```bash
cp2k.ssmp -i h2o.inp | tee h2o.out
```

### 11.2 并行方式：MPI × OpenMP

`cp2k.psmp` 同时支持 MPI 进程和 OpenMP 线程：**MPI 进程数 × 每进程线程数 = 分配到的 CPU 核数**。

```bash
# 单机，OpenMP 版（ssmp）：只用线程
export OMP_NUM_THREADS=8
cp2k.ssmp -i h2o.inp -o h2o.out

# 单机，纯 MPI：4 进程 × 1 线程（OpenMPI 写法）
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out

# MPICH / Intel MPI 写法
mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k.psmp -i h2o.inp -o h2o.out

# 混合并行：8 进程 × 4 线程 = 32 核
export OMP_NUM_THREADS=4
mpirun -np 8 -x OMP_NUM_THREADS cp2k.psmp -i big.inp -o big.out
```

经验法则：

- 内存够用时，**以 MPI 为主**（每进程 1–2 线程）通常最快；官方也指出并非所有功能都针对 OpenMP 优化过。
- 内存紧张（杂化泛函、大体系对角化）时，减少每节点 MPI 进程数、增加线程数，可以降低总内存。
- 运行前设置 `ulimit -s unlimited`，避免大体系下 Fortran 栈溢出导致的段错误。
- 不同 MPI×OpenMP 组合的性能差别可能很大，生产计算前用一个代表性输入跑 2–3 种组合比较 timing 表（[第18章](18-troubleshooting.md)）。

### 11.3 集群作业脚本（Slurm 最小示例）

```bash
#!/bin/bash
#SBATCH --job-name=h2o
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00

source ~/src/cp2k-2026.2/install/cp2k_env        # 或 module load / conda activate
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}   # Slurm >= 22.05 的 srun 不再继承 --cpus-per-task
ulimit -s unlimited

rm -f h2o.out
srun --cpu-bind=cores cp2k.psmp -i h2o.inp -o h2o.out
grep -q "PROGRAM ENDED AT" h2o.out && echo "CP2K finished normally"
```

**计算节点上的环境不会自动继承登录节点的设置**，`source cp2k_env`/`module load`/`conda activate` 必须写进作业脚本。更完整的模板（MPI×OpenMP 参数化、结束检查、重启提示）见本节 [`job_cp2k_slurm.sh`](../scripts/job_cp2k_slurm.sh) 与 [`job_cp2k_pbs.sh`](../scripts/job_cp2k_pbs.sh)（用法见脚本头部注释）；Slurm/PBS/LSF 的通用知识见 [公共篇第2章](../../Common/tutorials/02-job-scheduler.md)。

---

## 12. 安装自检：环境脚本与 H2O 冒烟测试

本章示例目录 [`examples/01-install-check/`](../examples/01-install-check/README.md) 提供三个文件：

| 文件 | 作用 |
|---|---|
| [`env_cp2k.sh`](../examples/01-install-check/env_cp2k.sh) | 按安装方式（static / conda / toolchain / spack / module）设置 `PATH`、`CP2K_DATA_DIR`、`CP2K_EXE`，用 `source` 加载 |
| [`check_install.sh`](../examples/01-install-check/check_install.sh) | 自动检查：可执行文件、版本与 `cp2kflags`、数据目录关键文件、输入语法检查、运行冒烟测试、检测 MPI 错配 |
| [`h2o_smoke_test.inp`](../examples/01-install-check/h2o_smoke_test.inp) | 孤立 H2O 分子 PBE 单点能，数秒内完成 |

### 12.1 冒烟测试输入

下面是去掉注释后的内容，完整带注释的文件见 [`h2o_smoke_test.inp`](../examples/01-install-check/h2o_smoke_test.inp)：

```text
&GLOBAL
  PROJECT h2o_smoke
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QS
  &DFT
    ! Two basis files are listed on purpose: this also checks that the data
    ! directory contains both the classic and the UZH basis libraries.
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 300
      REL_CUTOFF 40
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 50
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 8.0 8.0 8.0
      PERIODIC NONE
    &END CELL
    &COORD
      O   0.000000   0.000000   0.117300
      H   0.000000   0.756990  -0.468540
      H   0.000000  -0.756990  -0.468540
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

各部分的含义在 [第2章](02-input-structure.md) 逐行讲解；这里只需知道它是一个“孤立分子 + 非周期 Poisson 求解器 + PBE + DZVP 基组”的最小计算。

### 12.2 运行自检脚本

```bash
cd CP2K/examples/01-install-check
source env_cp2k.sh                 # 先按你的安装方式编辑文件开头的 CP2K_METHOD
bash check_install.sh              # 默认 1 个进程
NP=4 bash check_install.sh         # 用 4 个 MPI 进程跑冒烟测试，同时检查 MPI 是否配对
```

（示意输出，数值仅供参考）

```text
[ OK ] executable : /home/user/src/cp2k-2026.2/install/bin/cp2k.psmp
[ OK ] version    : CP2K version 2026.2
[INFO] cp2kflags  : omp libint fftw3 libxc elpa parallel scalapack mpi_f08 ...
[ OK ] data dir   : /home/user/src/cp2k-2026.2/install/share/cp2k/data
[ OK ] syntax check (-c) passed
[ OK ] SCF converged      : *** SCF run converged in    11 steps ***
[ OK ] total energy       : ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -17.16...
[ OK ] MPI ranks reported by DBCSR: 4 (expected 4)
[ OK ] program ended normally
```

### 12.3 怎样判断安装“真正可用”

| 检查项 | 看什么 | 异常说明 |
|---|---|---|
| 程序正常结束 | 输出末尾有 `PROGRAM ENDED AT` | 没有：查看 stderr 和输出中的 `ABORT` 框 |
| SCF 收敛 | `*** SCF run converged in N steps ***` | 冒烟测试几乎不可能不收敛；若不收敛多半是坏的数学库 |
| 总能量 | `grep "ENERGY| Total FORCE_EVAL" h2o_smoke.out` | 与静态二进制或容器版本的结果在 $10^{-6}$ Ha 量级一致为宜 |
| 数据目录 | `CP2K| Data directory path` | 指向了不存在的旧路径时，基组文件会找不到 |
| MPI 进程数 | `DBCSR| MPI: Number of processes` 应**只出现一次**且等于 `-np` | 出现多次且每次为 1：MPI 错配，见下节 |
| 警告数 | `The number of warnings for this run is : N` | 冒烟测试应为 0 或很少；逐条看 `*** WARNING` |

也可以用本节的解析脚本汇总输出（SCF 是否收敛、能量、警告、耗时）：`python ../scripts/cp2k_outparse.py -h` 查看用法。

---

## 安装方式与运行参数选择建议

| 项目 | 推荐值/做法 | 说明 |
|---|---|---|
| 学习、调试输入 | 静态二进制或 conda | 与集群版本保持同一个大版本号（如都用 2026.2） |
| 集群生产计算 | 管理员 module 或 toolchain/Spack 源码编译的 `cp2k.psmp` | 用集群自己的 MPI 才能发挥高速网络性能 |
| 并行布局 | 先试 MPI 为主（`OMP_NUM_THREADS=1` 或 2） | 杂化泛函、大内存任务改为每节点少进程多线程 |
| 数据目录 | 显式 `export CP2K_DATA_DIR=...` | 并在输出头部核对 `CP2K| Data directory path` |
| 栈大小 | `ulimit -s unlimited` | 写进作业脚本 |
| 输出文件 | 每次运行前删除旧 `.out` | `-o` 是追加写入 |
| 记录版本 | 保存一份 `cp2k.psmp -v` 的输出 | 写论文、对比结果时用 |
| 编译优化 | 登录/计算节点 CPU 不同时用 `--target-cpu=generic` 或在计算节点编译 | 避免 `Illegal instruction` |

---

## 常见问题与报错

| 报错或现象 | 原因 | 解决 |
|---|---|---|
| `error while loading shared libraries: libxxx.so: cannot open shared object file` | 运行环境没有加载编译时的依赖库（作业脚本里漏了 `source install/cp2k_env` 或 `module load`） | 在作业脚本中显式加载环境；`ldd $(which cp2k.psmp) \| grep "not found"` 查缺哪个库 |
| 输出中 `DBCSR| MPI: Number of processes  1` 出现多次，`PROGRAM STARTED AT` 也重复出现，输出交错 | **MPI 错配**：启动用的 `mpirun` 与 CP2K 链接的 MPI 不是同一套，结果是 N 个互相独立的串行进程 | `which mpirun` 与 `ldd cp2k.psmp \| grep -i mpi` 对照；toolchain 版先 `source cp2k_env`，conda 版先 `conda activate`；不要混用 module MPI 和 conda MPI |
| `Illegal instruction (core dumped)` | 用 `native` 优化编译的程序跑在指令集更老的 CPU 上 | 在计算节点重新编译，或 toolchain 加 `--target-cpu=generic` |
| `cannot execute binary file: Exec format error` | 架构不符（如在 aarch64 机器上运行 x86_64 静态二进制） | 换对应架构的容器镜像或源码编译 |
| `The specified OLD file <BASIS_MOLOPT_UZH> cannot be opened. It does not exist. Data directory path: ...` | 数据目录未设置或指向错误/过旧的位置 | 设置正确的 `CP2K_DATA_DIR`；或把文件放到计算目录 |
| `ERROR: ... Incorrect execution location.`（toolchain） | 没有在 `tools/toolchain/` 目录内执行脚本 | `cd tools/toolchain` 后再运行 |
| toolchain 下载失败、证书错误 | 集群无外网或 wget 太旧 | 在外网机器下载所需压缩包放入 `tools/toolchain/build/`；或升级 wget；最后手段 `--no-check-certificate` |
| CMake 报找不到 DBCSR 或 DBCSR 版本信息错误 | 在带 tag 的 git 克隆目录中构建（issue #5184） | 发布版改用 `cp2k-2026.2.tar.bz2` 压缩包 |
| 编译时报 GCC 版本不受支持 | 2026.2 不再支持 GCC 8 | `module load` 更新的 GCC，或让 toolchain 安装 GCC |
| Mac 上 `conda install cp2k` 只装到 7.1.0 | conda-forge 的 osx-64 停在 7.1.0 | 改用 Podman + aarch64 镜像或 `brew install cp2k` |
| 容器里报 `Input file not found` | 输入不在挂载目录内（Mac 上必须在 `/Users` 下；Linux 上忘了 `-v $PWD:/mnt` 或 `-B $PWD`） | 在输入所在目录执行，并正确挂载 |
| 容器中 CP2K 被杀死、退出码 137 | 容器/虚拟机内存不足（OOM） | 增大 Podman machine 内存或减少 MPI 进程数 |
| MPI 报 `/dev/shm` 共享内存错误 | 容器未加 `--shm-size=1g` | 加上该参数 |
| `/bin/bash^M: bad interpreter` | 作业脚本是 Windows CRLF 换行 | `dos2unix job.sh` |
| 输出里同时有新旧两次运行的内容 | `-o` 追加写入 | 运行前删除旧输出 |

---

## 练习

1. 下载 2026.2 静态二进制并解压同版本 `data/`，运行 `cp2k.ssmp -v`，把 `cp2kflags` 中每个标志的含义对照 §3.3 的表格写下来；判断它能否做 PBE0 杂化泛函计算。
2. 在同一台机器上分别用静态二进制和 conda（或容器）运行 `h2o_smoke_test.inp`，比较两者的总能量与耗时（`grep "ENERGY| Total FORCE_EVAL"`、timing 表最后的 `CP2K` 行）。
3. 故意不加载 MPI 环境（或在 conda 环境外调用系统 `mpirun`）运行 `NP=4 bash check_install.sh`，观察 MPI 错配的症状，然后修复。
4. 把 `CP2K_DATA_DIR` 设为一个不存在的目录再运行冒烟测试，阅读报错信息中的 `Data directory path`；再把 `GTH_POTENTIALS` 复制到计算目录，确认 CP2K 优先读取当前目录的文件。
5. 用 `cp2k.psmp -e h2o_smoke_test.inp > echoed.inp` 导出“所有默认值都显式写出”的输入，找出 `&MGRID` 中 `NGRIDS`、`PROGRESSION_FACTOR` 的默认值，以及 `&SCF` 中 `SCF_GUESS` 的默认值。

---

## 参考资料

**官方文档**

- CP2K 手册（trunk，含方法文档与输入参考）：<https://manual.cp2k.org/trunk/>
- CP2K 2026.2 版本手册：<https://manual.cp2k.org/cp2k-2026_2-branch/index.html>；各版本手册索引：<https://manual.cp2k.org/trunk/versions.html>
- 发行版安装（conda-forge、Docker/Podman、apt、brew、EasyBuild、Spack……）：<https://manual.cp2k.org/trunk/getting-started/distributions.html>
- 源码构建（CMake、toolchain + build_cp2k.sh、make_cp2k.sh）：<https://manual.cp2k.org/trunk/getting-started/build-from-source.html>
- 用 Spack 构建：<https://manual.cp2k.org/trunk/getting-started/build-with-spack.html>
- macOS 安装（Podman + aarch64 镜像）：<https://manual.cp2k.org/trunk/getting-started/macos.html>
- 第一个计算：<https://manual.cp2k.org/trunk/getting-started/first-calculation.html>
- 官方排错指南：<https://manual.cp2k.org/trunk/getting-started/troubleshooting.html>
- v2026.2 发布说明：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>；v2026.1 发布说明：<https://github.com/cp2k/cp2k/releases/tag/v2026.1>
- Docker Hub 镜像：<https://hub.docker.com/r/cp2k/cp2k>；conda-forge 包：<https://anaconda.org/conda-forge/cp2k>
- CP2K 用户论坛（Google Group）：<https://groups.google.com/g/cp2k>

**相关 GitHub 仓库**（只链接，不复制内容；更多见 [Resources/github-repos.md](../../Resources/github-repos.md)）

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：官方源码；`docs/` 是手册源文件，`tests/` 下数千个回归测试输入是查关键字真实写法最权威的来源，`data/` 是基组与赝势。
- [cp2k/cp2k-containers](https://github.com/cp2k/cp2k-containers)：Apptainer/Docker 使用说明（含 `run_tests` 自检与 CUDA 镜像）；2026.2 起 Dockerfile 迁入主仓库 `tools/docker/`。
- [conda-forge/cp2k-feedstock](https://github.com/conda-forge/cp2k-feedstock)：conda-forge 配方，可查看各平台提供的版本与 MPI 变体。
- [hfp/xconfigure](https://github.com/hfp/xconfigure)：Intel 编译器下编译 CP2K 的历史配方（基于旧 Makefile，仅作参考）。
- [Shxnek/cp2k-windows-installer](https://github.com/Shxnek/cp2k-windows-installer)：中文 Windows 一键安装包，仅适合教学演示。
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含 CP2K 入门与集群使用经验。
- [EPCCed/intro_to_cp2k](https://github.com/EPCCed/intro_to_cp2k)：EPCC/ARCHER2 入门课，含并行性能练习。

**关键论文**

- T. D. Kühne *et al.*, CP2K: An electronic structure and molecular dynamics software package – Quickstep: Efficient and accurate electronic structure calculations, *J. Chem. Phys.* 152, 194103 (2020). DOI: [10.1063/5.0007045](https://doi.org/10.1063/5.0007045)
- M. Iannuzzi *et al.*, The CP2K Program Package Made Simple, *J. Phys. Chem. B* 130, 1237–1310 (2026). DOI: [10.1021/acs.jpcb.5c05851](https://doi.org/10.1021/acs.jpcb.5c05851)
- J. VandeVondele *et al.*, Quickstep: Fast and accurate density functional calculations using a mixed Gaussian and plane waves approach, *Comput. Phys. Commun.* 167, 103–128 (2005). DOI: [10.1016/j.cpc.2004.12.014](https://doi.org/10.1016/j.cpc.2004.12.014)
- G. Lippert, J. Hutter, M. Parrinello, A hybrid Gaussian and plane wave density functional scheme, *Mol. Phys.* 92, 477–487 (1997). DOI: [10.1080/002689797170220](https://doi.org/10.1080/002689797170220)

---

[返回目录](../README.md) | [下一章 →](02-input-structure.md)
