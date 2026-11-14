# 示例 02：编译安装 OpenMX 4.0 + patch 4.0.1 并验证（Build & Verify OpenMX）

配套教程：[第2章 编译与安装（Intel / GCC / ScaLAPACK / ELPA / WSL / deb / Docker）](../../tutorials/02-install.md)

## 目的

1. 提供 GCC 与 Intel oneAPI 两条路线的 makefile 片段（只含 `MKLROOT`/`CC`/`FC`/`LIB`，不含官方 makefile 原文）；
2. 提供一键脚本：下载源码与补丁 → 打 patch 4.0.1 → 备份并改写 makefile → `make install` + `make all` → 运行 `-runtest`；
3. 提供 `runtest.result` 的自动判读脚本和运行环境设置脚本。

## 文件说明

| 文件 | 说明 |
|---|---|
| `makefile.gcc.example` | GCC + Open MPI + ScaLAPACK + OpenBLAS + FFTW3 的 `CC`/`FC`/`LIB`（Ubuntu 24.04、Debian 12–13 的库名），附各选项用途注释 |
| `makefile.oneapi.example` | Intel oneAPI（`mpiicx`/`mpiifx`）+ MKL + Intel MPI 的 `MKLROOT`/`CC`/`FC`/`LIB` |
| `build_openmx_gcc.sh` | GCC 路线一键脚本；可选 `INSTALL_DEPS=1` 用 apt 安装依赖；按 GCC 版本自动追加兼容选项，aarch64 自动加 `-Dnosse` |
| `build_openmx_oneapi.sh` | oneAPI 路线一键脚本；自动 `source setvars.sh`，AMD CPU 自动把 `-xHOST` 换成 `-march=x86-64-v3` |
| `check_runtest.sh` | 读取 `runtest.result`（或 `runtestL.result`），逐个检查 `diff Utot`、`diff Force` 是否小于 1e-7，输出 PASS/FAIL |
| `env_openmx.sh` | 运行环境（`source` 使用）：PATH、可选的 oneAPI 与 deb 运行库、线程与栈设置，以及把输入文件 `DATA.PATH` 改成本机路径的函数 `omx_setdata` |

两个构建脚本修改 makefile 的方式相同：先把原文件备份为 `makefile.orig`，把其中**生效的** `MKLROOT`/`CC`/`FC`/`LIB` 定义（含 `\` 续行）注释掉，再把新的定义插到文件开头；`CC += -I$(LIBELPADIR)`、`FC += -I$(LIBELPADIR)` 两行保持不变（内置 ELPA 需要它们）。

## 运行步骤

### GCC 路线（Ubuntu 24.04 / Debian / WSL）

```bash
cd /path/to/DFT-Tutorials/OpenMX/examples/02-install-build
INSTALL_DEPS=1 bash build_openmx_gcc.sh          # first time: also installs apt packages (sudo)
# later rebuilds:
FORCE=1 bash build_openmx_gcc.sh
# custom location, skip the test:
PREFIX=/apps/openmx RUN_TEST=0 bash build_openmx_gcc.sh
```

### Intel oneAPI 路线（WSL 或集群）

```bash
source /opt/intel/oneapi/setvars.sh               # or: module load intel-oneapi-...
bash build_openmx_oneapi.sh
# AMD CPU or heterogeneous cluster: choose the CPU flag yourself
ARCH_FLAG="-march=core-avx2" bash build_openmx_oneapi.sh
# Intel compilers + Open MPI + MKL (cluster):
BLACS_LIB=-lmkl_blacs_openmpi_lp64 bash build_openmx_oneapi.sh   # also make sure mpiicx/mpiifx wrap your MPI
```

### 手动编译（不用脚本）

```bash
cd ~/software/openmx4.0/source
cp makefile makefile.orig
# edit makefile: comment out the active MKLROOT/CC/FC/LIB, paste makefile.gcc.example (or .oneapi.example)
make clean && make install && make all
```

### 验证

```bash
cd ~/software/openmx4.0/work
mpirun -np 1 ./openmx Methane.dat -nt 1 | tee met.std     # Utot should be -8.033540953061 Ha (manual Sec. 5)
mpirun -np 4 ./openmx -runtest -nt 1
bash /path/to/DFT-Tutorials/OpenMX/examples/02-install-build/check_runtest.sh runtest.result
```

### 日常使用

```bash
OPENMX_ROOT=$HOME/software/openmx4.0 source env_openmx.sh     # Intel build: add OMX_USE_ONEAPI=1
omx_setdata *.dat                                             # rewrite DATA.PATH in your inputs
mpirun -np 4 openmx Methane.dat -nt 1 | tee met.std
```

deb 包用户：`OPENMX_ROOT=/opt/openmx4.0 OMX_DEB_RUNTIME=1 source env_openmx.sh`（数据库目录若不是 `/opt/openmx4.0/DFT_DATA19`，另设 `OMX_DATA`）。

## 预期结果 / 如何检查

- 构建脚本最后应打印 `built: work/DosMain`、`built: work/bandgnu13` 等行，且 `ldd work/openmx` 中没有 `not found`；
- `$PREFIX/openmx4.0/source/make_install.log`、`make_all.log` 保存了完整编译日志，`$PREFIX/openmx4.0/INSTALL_NOTES.txt` 记录了编译器版本与所用选项；
- `check_runtest.sh` 对 14 个测试逐行给出 `ok`/`FAIL`，最后一行为 `RESULT: PASS` 即安装正确（判据：`diff Utot` 与 `diff Force` 均 < 1e-7，见手册第 6 节）。输出格式如下（示意输出，数值仅供参考）：

```text
  #  input                                       time(s)            diff Utot           diff Force  status
  1  input_example/Benzene.dat                      4.93       0.000000000000       0.000000000000  ok
...
tests parsed: 14   failed (>= 1e-7): 0   total elapsed time: 99.01 s
RESULT: PASS
```

## 注意事项

- 脚本默认安装到 `$HOME/software/openmx4.0`；目录已存在时会停止，需 `FORCE=1` 才会删除重建；
- **不要**用 `make CC=... FC=... LIB=...` 的命令行方式覆盖变量，否则 makefile 中的 `CC += -I$(LIBELPADIR)` 会失效；
- Ubuntu/Debian 的 ScaLAPACK 库名是 `-lscalapack-openmpi`；其他系统一般是 `-lscalapack`，可用 `SCALAPACK_LIB="-L/path/lib -lscalapack"` 覆盖；
- 编译 `openmx` 需要十几分钟到半小时；脚本不使用并行 make（有用户报告并行编译时 ELPA 模块依赖出错）；
- 在 Docker 容器中以 root 运行时，Open MPI 会拒绝启动 `mpirun`（`mpirun has detected an attempt to run as root`），请改用普通用户；
- `-xHOST` 与 `-march=native` 只适合编译机与计算机 CPU 相同的情况；
- 脚本面向 Linux（含 WSL）编写，编写时验证了 bash 语法、makefile 改写逻辑与 `check_runtest.sh` 的解析（用手册中的 runtest 输出格式），但无法覆盖所有发行版与编译器版本的实际编译；遇到编译错误请对照[教程第2章](../../tutorials/02-install.md)“常见问题与报错”。仓库中的通用版本 [`omx_build.sh`](../../scripts/omx_build.sh) 提供相同功能（用法见 `bash ../../scripts/omx_build.sh -h`）。
