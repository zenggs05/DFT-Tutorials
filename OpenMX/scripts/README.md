# OpenMX 辅助脚本（OpenMX helper scripts）

本目录收录 OpenMX 4.0 计算从编译到后处理的辅助脚本：一键编译并跑 runtest（`omx_build.sh`），从结构文件生成输入（`omx_make_input.py`），写初始磁构型（`omx_magconf.py`），截断能/k 点/基组收敛测试（`omx_conv_test.py`），状态方程拟合（`omx_eos_fit.py`），输出解析（`omx_parse_out.py`），能带、DOS/PDOS、cube 文件、轨迹和 NEGF 输运的后处理（`omx_band_plot.py`、`omx_dos_auto.py`、`omx_cube_tools.py`、`omx_traj.py`、`omx_negf_plot.py`），以及 Slurm/PBS 作业模板。输出格式均按 OpenMX 4.0（+ patch 4.0.1）源码和手册核对，同时尽量兼容 3.9 的写法（差异见各节说明）。

> [!NOTE]
> 所有 Python 脚本都支持 `-h` 查看完整参数（带子命令的脚本用 `python <脚本> <子命令> -h`），文件开头的 docstring 也写了格式说明和用法示例。三个 shell 脚本可用 `bash <脚本> -h` 查看用法。本 README 中的参数均直接取自各脚本的 argparse 定义或 `usage()` 函数，命令行示例只演示用法，没有附带运行结果。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [omx_build.sh](omx_build.sh) | 下载 openmx4.0 + patch4.0.1，按 GCC 或 oneAPI 工具链改写 makefile，编译全部工具并跑 `openmx -runtest`，汇总 diff Utot / diff Force | bash、wget/curl、make；GCC+MPI+ScaLAPACK+OpenBLAS+FFTW3 或 Intel oneAPI | [第2章](../tutorials/02-install.md) |
| [omx_make_input.py](omx_make_input.py) | 从 POSCAR/CIF/xyz 生成 `.dat`：内置 DFT_DATA19 推荐基组表与价电子数、初始自旋占据（FM/AFM/MAGMOM）、按 k 间距设 `scf.Kgrid`、可选 seekpath 能带路径 | numpy、ase（或 pymatgen）；可选 spglib、seekpath | [第3章](../tutorials/03-input-first-calc.md)、[第4章](../tutorials/04-pao-vps.md)、[第8章](../tutorials/08-band.md) |
| [omx_magconf.py](omx_magconf.py) | 改写坐标行第 6、7 列的上/下自旋初始电子数（FM/AFM/逐原子），非共线时写第 8–13 列，并检查每个原子的电子数守恒 | 仅标准库 | [第11章](../tutorials/11-spin-noncollinear.md)、[第13章](../tutorials/13-dft-u.md) |
| [omx_conv_test.py](omx_conv_test.py) | 截断能 / k 网格 / k 间距 / PAO 基组收敛测试：由模板批量生成输入、依次运行、解析 Utot、输出 CSV 与 ΔE 曲线并给出推荐值 | 标准库 + numpy；作图需 matplotlib | [第4章](../tutorials/04-pao-vps.md)、[第5章](../tutorials/05-scf-kpoints-convergence.md) |
| [omx_eos_fit.py](omx_eos_fit.py) | 按晶格缩放生成一组输入（固定 `scf.Ngrid`）、运行并拟合 Birch–Murnaghan / Murnaghan / Vinet 状态方程；也可拟合 `MD.Type EvsLC` 的 `.EvsLC` 文件 | numpy；作图需 matplotlib；Murnaghan/Vinet 需 scipy | [第6章](../tutorials/06-optimization-eos-neb.md) |
| [omx_parse_out.py](omx_parse_out.py) | 解析 `System.Name.out` / `.std` / `.md`：Utot 与能量分解、化学势、SCF 历史与收敛判断、Mulliken 电荷与磁矩、原子力、优化历史、网格与耗时，输出摘要/JSON/CSV；也可作模块导入 | 仅标准库 | [第3章](../tutorials/03-input-first-calc.md)、[第5章](../tutorials/05-scf-kpoints-convergence.md)、[第10章](../tutorials/10-charge-mo-cube.md)、[第18章](../tutorials/18-parallel-on-errors.md) |
| [omx_band_plot.py](omx_band_plot.py) | 读 `.Band`（Hartree）或 bandgnu13 的 `BANDDAT1/2`（eV），画自旋分辨能带，报告 VBM/CBM/带隙，可叠加对比（如有/无 SOC） | numpy、matplotlib | [第8章](../tutorials/08-band.md)、[第12章](../tutorials/12-soc-topology.md)、[第14章](../tutorials/14-wannier.md) |
| [omx_dos_auto.py](omx_dos_auto.py) | 非交互驱动 DosMain（四面体/高斯展宽，DOS/PDOS，原子列表），再把 `*.PDOS.*` 按元素/物种/原子与 s/p/d/f 求和，处理自旋列并作图 | numpy、matplotlib；DosMain | [第9章](../tutorials/09-dos-pdos-cohp-optics.md)、[第13章](../tutorials/13-dft-u.md) |
| [omx_cube_tools.py](omx_cube_tools.py) | cube 工具箱：信息、积分、平面平均、功函数（替代 gcube2oned）、差分电荷（替代 diff_gcube/add_gcube）、线性组合、\|ψ\|²、Bader 前处理 | numpy；`--plot` 需 matplotlib | [第10章](../tutorials/10-charge-mo-cube.md)、[第15章](../tutorials/15-berry-efield.md)、[第16章](../tutorials/16-surface-esm-workfunction.md) |
| [omx_traj.py](omx_traj.py) | 读 `.md` 轨迹：能量、最大力、温度、晶胞体积随步数变化，导出 extxyz / ASE `.traj` / dpdata（DeePMD 等） | numpy、matplotlib；可选 ase、dpdata | [第6章](../tutorials/06-optimization-eos-neb.md)、[第7章](../tutorials/07-md.md) |
| [omx_negf_plot.py](omx_negf_plot.py) | NEGF Step 3 后处理：k 平均透射谱、Landauer 电导、Seebeck 系数、电子热导、k 分辨透射/电导图、I–V 曲线 | numpy、matplotlib | [第17章](../tutorials/17-negf-transport.md) |
| [job_omx_slurm.sh](job_omx_slurm.sh) | Slurm 作业模板（MPI × OpenMP，`-nt` 与 `OMP_NUM_THREADS` 一致），运行后检查是否正常结束和 SCF 是否达到 `scf.maxIter` | bash、awk、Slurm | [第18章](../tutorials/18-parallel-on-errors.md)、[作业调度通用章节](../../Common/tutorials/02-job-scheduler.md) |
| [job_omx_pbs.sh](job_omx_pbs.sh) | PBS Pro / OpenPBS / Torque 作业模板，功能同上 | bash、awk、PBS/Torque | [第18章](../tutorials/18-parallel-on-errors.md)、[作业调度通用章节](../../Common/tutorials/02-job-scheduler.md) |

## 安装依赖与通用约定

```bash
# Python >= 3.9
pip install numpy matplotlib          # band_plot / dos_auto / cube_tools / traj / negf_plot / conv_test / eos_fit
pip install ase                       # omx_make_input.py (or: pip install pymatgen); omx_traj.py --traj
pip install spglib seekpath           # optional: omx_make_input.py --primitive / --kpath
pip install scipy                     # optional: omx_eos_fit.py --eos murnaghan|vinet
pip install dpdata                    # optional: omx_traj.py --dpdata
# omx_parse_out.py and omx_magconf.py use the standard library only
```

通用约定：

- **脚本之间的依赖**：请把整个 `scripts/` 目录放在一起使用。`omx_make_input.py` 提供内置 DFT_DATA19 表和 `valence_of()`、`recommended_basis()`、`default_vps()`、`read_vps_valence()` 等函数，`omx_magconf.py`、`omx_parse_out.py`、`omx_conv_test.py` 会从同目录导入它们。`omx_conv_test.py` 和 `omx_eos_fit.py` 用同目录的 `omx_parse_out.py` 读取输出（`omx_conv_test.py` 找不到时退回到简单的正则解析）。
- **价电子数**：凡是需要价电子数的脚本，都优先从 `DATA.PATH/VPS/<vps>.vps` 中的 `valence.electron` 行读取（`DATA.PATH` 取自输入文件或 `--data-path`），找不到 VPS 文件时才用内置表。因此在集群上运行这些脚本时，建议给出真实的 DFT_DATA19 路径。
- **作图**：matplotlib 只在作图时导入，默认用 Agg 后端直接保存图片（无图形界面的集群也能用），加 `--show` 才弹出交互窗口。多数脚本有 `--no-plot`，此时只打印结果、写数据文件。
- **单位**：OpenMX 的 `.out`、`.Band`、`.md` 中能量为 Hartree，力为 Hartree/Bohr，cube 文件长度为 Bohr。脚本会在输出中写明单位；`omx_traj.py` 和 `omx_cube_tools.py` 按 CODATA 2018（1 Ha = 27.211386245988 eV，1 Bohr = 0.529177210903 Å）换算，`omx_band_plot.py` 与 bandgnu13 一致，用 27.2113845 eV/Ha。
- **运行程序的脚本**：`omx_conv_test.py`、`omx_eos_fit.py` 用 `<mpirun> <openmx> <input> -nt <nt>` 依次运行 OpenMX，启动器默认取环境变量 `MPIRUN`（否则 `mpirun -np 4`），可执行文件取 `OPENMX`（否则 `openmx`）。它们适合在交互节点或一个作业脚本里整体执行，不会自己提交作业。
- **版本**：第一次在新版本 OpenMX 上使用时，建议先用 `omx_parse_out.py` 看一个小算例的解析结果是否完整，再用于批量任务。

> [!IMPORTANT]
> 脚本生成或改写的输入文件（`omx_make_input.py`、`omx_magconf.py`、`omx_conv_test.py`、`omx_eos_fit.py`）都应先人工检查一遍：截断能、k 网格、基组等级、自旋设置和混合参数必须按 [第4章](../tutorials/04-pao-vps.md)、[第5章](../tutorials/05-scf-kpoints-convergence.md) 做收敛测试，不能直接把默认值当作收敛值。

---

## 1. omx_build.sh — 一键编译 OpenMX 4.0

**用途**：[examples/02-install-build](../examples/02-install-build/README.md) 中 `build_openmx_gcc.sh` / `build_openmx_oneapi.sh` 的通用版本，把 [第2章](../tutorials/02-install.md) 的编译流程合成一条命令：

1. 下载 `openmx4.0.tar.gz` 和 `patch4.0.1.tar.gz`（优先 wget，没有时用 curl；文件已存在或 `--offline` 时跳过下载），并检查压缩包是否完整；
2. 解压到 `PREFIX/openmx4.0`，在 `source/` 中解开补丁，把补丁里修正过的 `GaAs.dat` 移到 `work/`；
3. 备份 `source/makefile` 为 `makefile.orig`，注释掉其中有效的 `MKLROOT`/`CC`/`FC`/`LIB` 定义（含反斜杠续行），在文件开头插入所选工具链的定义；`CC += ...`、`FC += ...` 这类行（内置 ELPA 用）保持不变；
4. `make clean`、`make install`、`make all`，日志写到 `source/make_install.log` 和 `source/make_all.log`；若 `make all` 没有生成 bandgnu13，再单独用 `gcc bandgnu13.c -lm` 编译；用 `ldd` 检查 openmx 是否有找不到的共享库，并逐个报告 `openmx DosMain bandgnu13 cohp polB jx cif2omx` 是否已生成；编译参数记录到 `PREFIX/openmx4.0/INSTALL_NOTES.txt`；
5. 在 `work/` 中运行 `$MPIRUN -np N ./openmx -runtest -nt 1`，用内嵌的 awk 解析 `runtest.result`，逐个算例列出 diff Utot 和 diff Force，两者都小于 1e-7 记为 ok，最后给出 `RUNTEST: PASS` 或 `RUNTEST: FAIL`（手册第 6 节的判据是差值在 1e-7 左右以内）。

**工具链细节**：

| 工具链 | 编译器 | 自动添加的选项 | 链接库 |
|---|---|---|---|
| `gcc`（默认） | `mpicc` / `mpif90` | 总是加 `-fopenmp -Wno-implicit-function-declaration`；GCC ≥ 10 加 `-fcommon`；gfortran ≥ 10 加 `-fallow-argument-mismatch`；GCC ≥ 14 加 `-Wno-error=implicit-int -Wno-error=incompatible-pointer-types`；GCC ≥ 15 加 `-std=gnu17`；aarch64 加 `-Dnosse` | ScaLAPACK（`ldconfig` 能找到 `libscalapack-openmpi` 时用 `-lscalapack-openmpi`，否则 `-lscalapack`）+ `-lopenblas -lfftw3`；Fortran MPI 库按 MPI 种类选择：Open MPI（`mpif90 --showme:link` 可用）用 `-lmpi_mpifh`，MPICH 类用 `-lmpifort` |
| `oneapi` | `mpiicx` / `mpiifx` | `-fcommon -Wno-implicit-function-declaration -qopenmp`；CPU 选项在 AMD 上默认 `-march=x86-64-v3`，其他 CPU 默认 `-xHOST` | MKL（ScaLAPACK、BLACS、FFTW3 接口）+ Intel MPI；找不到 `mpiicx` 时自动 `source <oneapi-root>/setvars.sh`，`MKLROOT` 未设置时报错退出 |

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `--toolchain gcc\|oneapi` | 编译器/库组合（`gcc`） |
| `--prefix DIR` | 安装根目录，源码放在 `DIR/openmx4.0`（`$HOME/software`） |
| `--src-dir DIR` | 压缩包的存放/下载目录（`PREFIX/src`） |
| `--offline` | 不下载，压缩包必须已在 `--src-dir` 中 |
| `--force` | 先删除已有的 `PREFIX/openmx4.0`（不加时目录已存在会报错退出） |
| `--opt "FLAGS"` | 优化选项（`-O3`；不要用 `-ffast-math`） |
| `--arch-flag FLAG` | oneapi：CPU 选项（见上表）；gcc：可选，如 `-march=native`（仅当编译机与计算节点 CPU 相同时使用） |
| `--scalapack-lib "LIBS"` | gcc：ScaLAPACK 链接选项（自动判断，见上表） |
| `--blas-lib "LIBS"` | gcc：BLAS/LAPACK 链接选项（`-lopenblas`） |
| `--fftw-lib "LIBS"` | gcc：FFTW3 链接选项（`-lfftw3`） |
| `--extra-lib "LIBS"` | 额外链接选项，如 `-L/opt/fftw/lib` |
| `--blacs-lib LIB` | oneapi：MKL BLACS 库（`-lmkl_blacs_intelmpi_lp64`；配 Open MPI 时改为 `-lmkl_blacs_openmpi_lp64`） |
| `--oneapi-root DIR` | oneapi：含 `setvars.sh` 的安装目录（环境变量 `ONEAPI_ROOT`，否则 `/opt/intel/oneapi`） |
| `--np N` | runtest 使用的 MPI 进程数（4） |
| `--no-test` | 跳过 `openmx -runtest` |
| `-h`/`--help` | 显示帮助 |
| 环境变量 `MPIRUN` | runtest 用的启动器（`mpirun`，可包含额外选项） |

**依赖**：bash、wget 或 curl、tar、make、系统 gcc（`make all` 编译小工具时要用）。gcc 工具链还需要 MPI（Open MPI 或 MPICH）、ScaLAPACK、OpenBLAS（或其他 BLAS/LAPACK）、FFTW3，Ubuntu 24.04 上为 `build-essential gfortran openmpi-bin libopenmpi-dev libscalapack-openmpi-dev libopenblas-dev libfftw3-dev`；oneapi 工具链需要 Intel oneAPI 2024 及以上（icx/ifx、Intel MPI、MKL）。

**用法示例**：

```bash
bash omx_build.sh -h
bash omx_build.sh --toolchain gcc
bash omx_build.sh --toolchain oneapi --prefix /apps/openmx --np 8
bash omx_build.sh --toolchain gcc --scalapack-lib "-L/opt/scalapack/lib -lscalapack" --no-test
bash omx_build.sh --toolchain gcc --offline --src-dir ~/tarballs --force
MPIRUN="mpirun --oversubscribe" bash omx_build.sh --toolchain gcc --np 4
```

**输出**：可执行文件在 `PREFIX/openmx4.0/work/`（openmx、DosMain、bandgnu13 等）。脚本结束时会提示把 `PREFIX/openmx4.0/work` 加入 `PATH`，并在输入文件中使用 `DATA.PATH  PREFIX/openmx4.0/DFT_DATA19`；oneapi 工具链运行前还需 `source setvars.sh`。

**局限**：只面向 Linux（含 WSL）；无法覆盖所有发行版与编译器版本组合，编译失败时先看 `make_install.log`，再对照 [第2章](../tutorials/02-install.md) 的“常见问题与报错”。不编译 GPU 版，也不处理外部 ELPA。runtest 判 FAIL 时，脚本会提示把 `work/runtest.result` 与 `work/input_example/` 中的参考结果对比。

---

## 2. omx_make_input.py — 从结构文件生成 OpenMX 输入

**用途**：读入 ASE（或 pymatgen）能读取的结构文件（POSCAR、CIF、xyz……），写出一个完整的 OpenMX 4.0 `.dat` 文件，包括：

- `System.CurrrentDirectory`（官方拼写就是三个 r）、`System.Name`、`DATA.PATH`、`level.of.stdout/fileout`；
- `Species.Number` 与 `<Definition.of.Atomic.Species`：每个元素一个物种，PAO 和 VPS 取自内置表（摘自官方 4.0 手册 Sec. 12.3 的 Table 1–2，即 Database Ver. 2019，共 82 个 VPS 条目，每条含价电子数和 Quick/Standard/Precise 三档推荐基组，已包含 4.0 中 Sn 的基组变更）；
- `Atoms.SpeciesAndCoordinates`：周期体系用 `FRAC`，分子用 `Ang`；第 6、7 列为上/下自旋初始电子数，两者之和等于 VPS 价电子数；
- `Atoms.UnitVectors`（Å）；分子取原子范围加两侧各 `--vacuum` 的长方盒，并用 `scf.EigenvalueSolver Cluster`；
- `scf.*`：交换关联、自旋、电子温度、截断能、`scf.maxIter`、求解器（周期体系 `Band`）、按 k 间距算出的 `scf.Kgrid`、`scf.Generation.Kpoint GCenter`（只在 `--omx-version 4.0` 时写）、RMM-DIISK 混合参数、`scf.criterion`；
- `MD.Type`（不是 `nomd` 时另写 `MD.maxIter 200` 和 `MD.Opt.criterion 3.0e-4`），可选 `Band.dispersion` + `Band.kpath`（`--kpath`）和 `Dos.fileout`（`--dos`）。

**依赖**：`pip install numpy ase`（或用 pymatgen 代替 ase 读结构）；`--primitive` 需要 `spglib`，`--kpath` 需要 `seekpath`。

**用法示例**：

```bash
python omx_make_input.py POSCAR -o si.dat --data-path /opt/openmx4.0/DFT_DATA19
python omx_make_input.py Fe.cif -o fe.dat --spin on --moment Fe=3.0 --etemp 1000
python omx_make_input.py NiO_afm.vasp --spin on --moment Ni=1.8 --afm Ni -o nio_afm.dat
python omx_make_input.py POSCAR --magmom "2*3.0 2*-3.0 4*0" --spin on -o mag.dat
python omx_make_input.py POSCAR --kpath --kpath-density 30 -o band.dat
python omx_make_input.py h2o.xyz --vacuum 6 -o h2o.dat               # molecule -> Cluster
python omx_make_input.py POSCAR --level precise --basis Si=Si7.0-s3p3d2 -o si_precise.dat
python omx_make_input.py POSCAR --basis Fe=Fe6.0S-s3p2d2f1:Fe_PBE19S:14 -o fe_soft.dat   # PAO:VPS:valence
python omx_make_input.py POSCAR --md-type RF --kgrid 8 8 8 --omx-version 3.9 -o opt39.dat
```

**参数**（括号内为默认值）：

| 类别 | 参数 | 说明 |
|---|---|---|
| 通用 | `structure` | 结构文件（位置参数） |
| | `--format FMT` | 强制指定 ASE 文件格式（如 `vasp`、`cif`、`xyz`） |
| | `-o/--output FILE` | 输出 `.dat`（不给时写到标准输出，摘要信息写到标准错误） |
| | `--name NAME` | `System.Name`（输出文件名去掉扩展名，否则 `omx`） |
| | `--data-path DIR` | `DATA.PATH`（`../DFT_DATA19`）；该目录在本机存在时，会检查 PAO/VPS 文件是否存在，并从 VPS 文件读价电子数 |
| | `--omx-version 4.0\|3.9` | 4.0（默认）写 `scf.Generation.Kpoint GCenter`，3.9 不写该关键词 |
| 基组/赝势 | `--level quick\|standard\|precise` | 推荐基组等级（`standard`） |
| | `--soft` | Fe、Co、Ni、Cu、Zn 改用不含半芯态的 `*_PBE19S`（默认 `*_PBE19H`） |
| | `--basis EL=PAO[:VPS[:VALENCE]]` | 覆盖某元素的基组，可重复；表中没有的元素（如 Eu、Gd、Tb、Yb）必须三项都给，价电子数从 VPS 文件中查 |
| | `--xc` | `GGA-PBE`（默认）/ `LDA` / `LSDA-CA` / `LSDA-PW`；LDA 类使用 `*_CA19` 赝势 |
| 自旋 | `--spin off\|on\|nc` | `scf.SpinPolarization`（给了非零磁矩时默认 `on`，否则 `off`） |
| | `--moment EL=MU` | 某元素的初始磁矩（μB，= 上 − 下电子数），可重复 |
| | `--afm EL` | 该元素按文件顺序交替正负号，可重复 |
| | `--magmom 'N*val ...'` | VASP MAGMOM 风格的逐原子磁矩，优先于 `--moment` |
| SCF | `--ecut` | `scf.energycutoff`（220 Ry） |
| | `--kspacing` | 由 k 间距定 `scf.Kgrid`（0.20 Å⁻¹，含 2π） |
| | `--kgrid N1 N2 N3` | 显式 `scf.Kgrid`，优先于 `--kspacing` |
| | `--etemp` | `scf.ElectronicTemperature`（300 K；金属常用 500–1000 K） |
| | `--criterion` | `scf.criterion`（1e-8 Hartree） |
| | `--maxiter` | `scf.maxIter`（200） |
| | `--charge` | `scf.system.charge`（正值 = 空穴掺杂；非零时才写） |
| 其他 | `--md-type` | `MD.Type`，如 `nomd`（默认）、`RF`、`EF`、`BFGS`、`RFC5` |
| | `--primitive` | 用 spglib 约化为原胞 |
| | `--symprec` | spglib/seekpath 的对称性容差（1e-3） |
| | `--kpath` | 用 seekpath 生成 `Band.kpath`，并**把结构替换为 seekpath 的标准化原胞**，保证路径标签正确（Γ 写成 `g`） |
| | `--kpath-density` | 能带路径上每 Å⁻¹ 的 k 点数（25） |
| | `--dos` | 加 `Dos.fileout on`、`Dos.Erange -20 20`、`Dos.Kgrid` = 2 × `scf.Kgrid` |
| | `--vacuum` | 仅分子：盒子每侧的真空层（6 Å） |
| | `--fileout 0\|1\|2` | `level.of.fileout`（1；0 = 不输出 cube 文件） |

**行为要点**：

- **价电子数**：`--data-path` 下存在对应 VPS 文件时，从其中的 `valence.electron` 行读取（推荐做法），否则用内置表。LDA（`*_CA19`）赝势的价电子数不在表中，没有 VPS 文件时会借用 PBE19 表的值并打印警告，提示用 `grep valence.electron` 核对。
- **k 网格**：`n_i = max(1, ceil(|b_i| / kspacing))`（`|b_i|` 含 2π），只有结构的 pbc 标记为非周期的方向才直接取 1。POSCAR/CIF 读入后三个方向都是周期的，slab 的真空方向按间距也可能得到 2，此时请用 `--kgrid N1 N2 1` 显式指定。默认 0.20 Å⁻¹ 只是起点，收敛性用 `omx_conv_test.py kspacing` 检验。
- **磁矩**：每个原子的上/下电子数 = (价电子数 ± 磁矩)/2，脚本会检查磁矩不超过价电子数。给了磁矩却用 `--spin off` 时，磁矩被忽略（上 = 下）并给出警告。
- **混合参数**：固定写入 `RMM-DIISK`，初始/最小/最大混合权重 0.10/0.001/0.30，`scf.Mixing.History 30`，`StartPulay 10`，`EveryPulay 1`。难收敛体系请按 [第5章](../tutorials/05-scf-kpoints-convergence.md) 调整。
- 输出文件为纯 ASCII、LF 换行，可直接在 Linux 上使用。

**局限**：每个元素只建一个物种，不能区分同一元素的不同物种（如 AFM 时需要不同 U 值的情况，请手工拆分物种）；不写 DFT+U、SOC、vdW、ESM 等关键词；`--spin nc` 时坐标行只写到第 7 列，非共线的 Euler 角等列请再用 `omx_magconf.py --nc` 写入；分子盒只支持正交盒。与官方 `cif2omx` 的对比见 [第3章](../tutorials/03-input-first-calc.md) 和 [第4章](../tutorials/04-pao-vps.md)。

---

## 3. omx_magconf.py — 初始自旋构型

**用途**：对已有的 `.dat` 改写 `Atoms.SpeciesAndCoordinates` 每行的第 6、7 列（上/下自旋初始电子数），实现铁磁、反铁磁或逐原子自定义的初始磁构型，并保证每个原子的“上 + 下 = 价电子数”。非共线计算（`--nc`）额外写第 8–13 列：第 8、9 列为自旋磁矩的 Euler 角 θ、φ（度），第 10、11 列为轨道磁矩的角度，第 12 列为约束/Zeeman 项开关（1/0），第 13 列为 DFT+U 轨道极化开关（on/off）。规则依据手册 Sec. 8.2 与 Sec. 37；与官方 `Cr2_CNC.dat` 示例一致，允许“下 > 上”的负磁矩与角度同时使用。

**依赖**：仅标准库。价电子数依次取自 `--valence`、VPS 文件（输入中的 `DATA.PATH` 或 `--data-path`）、同目录 `omx_make_input.py` 的内置表。

**用法示例**：

```bash
# ferromagnetic bcc Fe, 3 muB per Fe (species labels or element symbols are accepted)
python omx_magconf.py fe.dat --mode fm --moment Fe=3.0 -o fe_fm.dat
# AFM: alternate signs over the Ni atoms in file order, O stays non-magnetic
python omx_magconf.py nio.dat --mode afm --moment Ni=1.8 -o nio_afm.dat
# AFM with an explicit sign pattern for atoms 1-4, overwrite the file (keeps nio.dat.bak)
python omx_magconf.py cell.dat --mode afm --moment Cr=3 --atoms 1-4 --signs "+ - - +" --inplace
# per-atom list (VASP MAGMOM style)
python omx_magconf.py mno.dat --mode custom --magmom "2*4.5 2*-4.5 4*0" -o mno_afm.dat
# non-collinear: spins in the xz plane, constraint switch on for magnetic atoms
python omx_magconf.py cr2.dat --mode afm --moment Cr=2 --nc --theta 90 --phi 0 --constrain -o cr2_nc.dat
# non-collinear with per-atom angles
python omx_magconf.py mno.dat --mode fm --moment Mn=4.5 --nc --angles 2:180:0 --angles 4:180:0 -o mno_nc.dat
# only check the electron count of an existing file
python omx_magconf.py nio.dat --check
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `input` | OpenMX 输入文件（位置参数） |
| `-o/--output FILE` / `--inplace` / `--block-only` | 三选一：写到新文件 / 覆盖原文件（保留 `.bak`）/ 只打印坐标块；都不给时完整文件写到标准输出 |
| `--mode fm\|afm\|custom\|none` | `fm`：同号；`afm`：交替或按 `--signs`；`custom`：按 `--magmom`；`none`：非磁（上 = 下）（`fm`） |
| `--moment KEY=MU` | 物种标签或元素符号的初始磁矩（μB），可重复 |
| `--atoms '1-4,9'` | 只对这些原子设置 fm/afm（默认：所有在 `--moment` 中出现的物种/元素的原子） |
| `--signs '+ - - +'` | afm 的符号序列（也可写 `'1 -1 -1 1'`），按所选原子的文件顺序（默认每个元素内交替） |
| `--magmom 'N*val ...'` | custom 模式的逐原子磁矩 |
| `--valence KEY=N` | 覆盖某物种/元素的价电子数，可重复 |
| `--data-path DIR` | DFT_DATA19 目录，用于读 VPS 中的 `valence.electron`（默认用输入中的 `DATA.PATH`，目录存在时） |
| `--nc` | 非共线：写第 8–13 列并把 `scf.SpinPolarization` 设为 `nc` |
| `--theta` / `--phi` | nc：自旋 Euler 角（0 / 0 度） |
| `--angles I:THETA:PHI` | nc：逐原子的自旋角，可重复 |
| `--constrain` | nc：磁性原子第 12 列写 1（约束/Zeeman），其余写 0 |
| `--orbpol` | nc：磁性原子第 13 列写 on（DFT+U 轨道极化） |
| `--no-set-spin` | 不改 `scf.SpinPolarization`（默认有非零磁矩时设为 on/nc） |
| `--check` | 只报告当前占据和电子数检查结果 |

**输出与检查**：标准错误中打印每个原子的价电子数、上/下电子数和磁矩，以及总价电子数、总初始磁矩和价电子数来源。afm 模式下某元素被选中的原子数为奇数时会警告（无法完全抵消）。

> [!NOTE]
> OpenMX 输出的磁矩符号约定与电子的物理磁矩相反（见 [第11章](../tutorials/11-spin-noncollinear.md)），FM/AFM 只看相对符号。带电体系用 `scf.system.charge` 处理，不要改第 6、7 列的电子数之和。非共线计算的耗时约为共线计算的 4 倍。

**局限**：不改物种定义（AFM 子格需要不同 U 值或不同 PAO 时请手工拆物种）；不写 `<Hubbard.U.values` 等 DFT+U 块（见 [第13章](../tutorials/13-dft-u.md)）；初始磁矩只是 SCF 的起点，最终磁构型要看输出的 Mulliken 磁矩。

---

## 4. omx_conv_test.py — 收敛测试驱动

**用途**：实现 [第5章](../tutorials/05-scf-kpoints-convergence.md) 的截断能与 k 点收敛测试，以及 [第4章](../tutorials/04-pao-vps.md) 的基组测试。由一个模板输入批量生成 `<workdir>/<prefix>_<tag>.dat`，依次运行 OpenMX，解析 Utot，计算每原子能量差 ΔE（meV/atom），写 CSV、画图并给出推荐值。

**依赖**：标准库 + numpy（`kspacing` 模式用）；作图需 matplotlib；运行需 OpenMX 可执行文件。同目录的 `omx_parse_out.py`（解析）和 `omx_make_input.py`（`basis` 模式的等级名）会被自动导入。

**模式**（第一个位置参数）：

| 模式 | `--values` 的含义 | 示例 |
|---|---|---|
| `ecut` | `scf.energycutoff`（Ry）；模板中有效的 `scf.Ngrid` 行会被注释掉（它会覆盖截断能） | `150 180 220 260 300 350` |
| `kgrid` | k 网格：`8`（= 8 8 8）或 `8x8x1` | `4 6 8 10 12 14` |
| `kspacing` | k 间距（Å⁻¹，含 2π），`n_i = ceil(|b_i|/spacing)`，得到相同网格的间距自动去重 | `0.40 0.30 0.25 0.20 0.15` |
| `basis` | PAO 字符串（配合 `--species`；模板只有一个物种时可省略），或 `quick`/`standard`/`precise`（查 DFT_DATA19 表，作用于全部物种或只作用于 `--species`） | `Si7.0-s2p1d1 Si7.0-s2p2d1 Si7.0-s3p3d2` |

**用法示例**：

```bash
python omx_conv_test.py ecut Si_scf.dat --values 120 150 180 220 260 300 350 400
python omx_conv_test.py kgrid Si_scf.dat --values 4 6 8 10 12 14 16 --mpirun "mpirun -np 8" --nt 2
python omx_conv_test.py kgrid slab.dat --values 4x4x1 6x6x1 8x8x1 10x10x1 12x12x1
python omx_conv_test.py kspacing Si_scf.dat --values 0.4 0.3 0.25 0.2 0.15 0.12
python omx_conv_test.py basis Si_scf.dat --species Si --values Si7.0-s2p1d1 Si7.0-s2p2d1 Si7.0-s3p3d2
python omx_conv_test.py basis Si_scf.dat --values quick standard precise
python omx_conv_test.py ecut Si_scf.dat --values 150 220 300 --dry-run          # inputs only
python omx_conv_test.py ecut Si_scf.dat --values 150 220 300 --collect          # parse only
python omx_conv_test.py kgrid Si_scf.dat --values 6 8 10 12 --collect --ref 12x12x12 --thresh 0.5
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `mode` / `template` | 模式与模板输入（位置参数）；模板中可放占位符 `@VALUE@` |
| `--values ...` | 要测试的值（必需，至少两个不同的值） |
| `--species LABEL` | basis 模式：要替换 PAO 的物种标签 |
| `--workdir DIR` | 输入/输出目录（`conv_<mode>`） |
| `--prefix NAME` | `System.Name` 前缀（模板的 `System.Name`） |
| `--mpirun "CMD"` | MPI 启动器及其选项（`$MPIRUN` 或 `mpirun -np 4`；`''` 表示串行） |
| `--openmx EXE` | OpenMX 可执行文件（`$OPENMX` 或 `openmx`） |
| `--nt N` | 每个 MPI 进程的 OpenMP 线程数，作为 `-nt` 传入，同时设置 `OMP_NUM_THREADS`（1） |
| `--data-path DIR` | 覆盖所有生成输入的 `DATA.PATH` |
| `--keep-fileout` | 保留模板的 `level.of.fileout`（默认改为 0，不写 cube 文件） |
| `--dry-run` | 只写输入 |
| `--collect` | 只解析已有输出，不运行 |
| `--rerun` | 已正常结束的点也重新运行 |
| `--thresh` | 收敛阈值（1.0 meV/atom） |
| `--ref LABEL` | 参考值的标签（默认为排序后的最后一个值，见下） |
| `--natoms N` | 每个晶胞的原子数（默认读 `Atoms.Number`） |
| `--csv` / `--plot` | CSV 与 PNG 文件名（`<workdir>/<mode>_conv.csv`、`<workdir>/<mode>_conv.png`） |
| `--no-plot` / `--show` | 不画图 / 弹出交互窗口 |

**工作方式**：

1. **改写输入**：按关键词名（不区分大小写）修改 `scf.energycutoff`、`scf.Kgrid` 或物种行的 PAO；如果模板含 `@VALUE@`（ecut/kgrid/kspacing 模式），则改为替换这个占位符。`System.Name` 设为 `<prefix>_<tag>`；相对路径的 `DATA.PATH` 会换成绝对路径，因为计算在 `<workdir>` 中进行；`level.of.fileout` 改为 0。
2. **运行**：每个点执行 `<mpirun> <openmx> <name>.dat -nt <nt> > <name>.std`；`.std` 中已有 “The calculation was normally finished.” 的点会被跳过（除非 `--rerun`）。
3. **解析**：从 `<name>.out` 读 Utot、SCF 是否收敛、实空间网格和耗时（通过 `omx_parse_out.py`，导入失败时退回正则解析）。
4. **参考值与推荐值**：ΔE = (Utot − Utot_ref)/Natoms。参考值是列表中的最后一个点：`ecut` 按截断能升序排列，`kspacing` 按间距降序排列（最密的网格在最后），`kgrid` 和 `basis` **保持命令行给出的顺序**，所以请按从粗到精的顺序写 `--values`，或用 `--ref` 指定标签（如 `300`、`12x12x12`、`0.15`、`Si7.0-s3p3d2`）。比较始终按标签进行，不按 N1×N2×N3 比较。推荐值是满足“从它开始、所有更精确的点都在 `--thresh` 以内且 SCF 收敛”的最便宜的值。

**输出**：屏幕表格；CSV 列为 `value, name, Utot_Ha, dE_meV_per_atom, scf_converged, finished, grid, elapsed_s, reference`；PNG 为上下两栏：上栏是带 ±阈值色带的 ΔE，下栏是 |ΔE| 的对数坐标（ecut/kspacing 横轴为数值，kgrid/basis 横轴为标签）。若找不到满足阈值的值，会提示扩大扫描范围；有 SCF 未收敛的点时会警告。

> [!TIP]
> kgrid/kspacing 模式会检查 `scf.EigenvalueSolver`：只有 `Band` 求解器才使用 `scf.Kgrid`，模板用 `Cluster` 时会给出警告。截断能测试中实空间网格数是离散变化的，ΔE 曲线可能出现台阶，CSV 的 `grid` 列可以帮助判断。

**局限**：只比较总能量，不检查力或应力；各点串行运行，不提交作业（在作业脚本中调用即可）；对 4.0 的 `GCenter` 与 3.9 的 `Regular` 网格不做区分，跨版本比较前请统一 `scf.Generation.Kpoint`。

---

## 5. omx_eos_fit.py — 状态方程拟合

**用途**：实现 [第6章](../tutorials/06-optimization-eos-neb.md) 的 E–V 曲线方法。把收敛过的单点输入按晶格均匀缩放，生成 `<prefix>_s<scale>.dat`，依次运行，然后拟合状态方程，给出 E0、V0、B0（GPa）、B0′ 和晶格常数 a0。

**关键点**：所有体积必须使用**相同的实空间网格**（手册 Sec. 14.2），否则网格数随体积跳变会使能量曲线出现台阶。网格来源按优先级为：

1. `--ngrid N1 N2 N3`；
2. `--ngrid-from FILE.out`（读取 “Num. of grids of a-, b-, and c-axes = ...”）；
3. 模板中已有的 `scf.Ngrid`；
4. 自动探测：对**最大的**晶胞以 `scf.maxIter 1`、`MD.Type nomd` 跑一次（`<prefix>_gridprobe`），取其网格，保证每个体积的截断能都不低于要求值。

缩放时所有晶格矢量乘以 s；`Ang`/`AU` 坐标同时缩放，`FRAC` 坐标不变。三阶 Birch–Murnaghan 用 V^(−2/3) 的三次多项式精确拟合（只需 numpy）；Murnaghan 与 Vinet 用 `scipy.optimize.curve_fit`。B0 按 1 eV/Å³ = 160.21766 GPa 换算。

**依赖**：numpy；作图需 matplotlib；`--eos murnaghan|vinet` 需 scipy；运行需 OpenMX。

**用法示例**：

```bash
python omx_eos_fit.py Si_scf.dat --range -3 3 7 --a-ref 5.43 --mpirun "mpirun -np 8"
python omx_eos_fit.py Fe_bcc.dat --scales 0.96 0.98 1.00 1.02 1.04 --ngrid 32 32 32
python omx_eos_fit.py Si_scf.dat --range -3 3 7 --dry-run --ngrid 24 24 24      # inputs only
python omx_eos_fit.py Si_scf.dat --range -3 3 7 --collect --a-ref 5.43          # fit existing runs
python omx_eos_fit.py Si_scf.dat --range -4 4 9 --collect --eos vinet --a-ref 5.43
python omx_eos_fit.py --evslc cdia_evslc.EvsLC --natoms 2 --a-ref-factor 4      # fit an EvsLC file
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `template` | 收敛过的单点输入（位置参数；用 `--evslc` 时可省略） |
| `--scales S ...` / `--range MIN% MAX% N` | 二选一：缩放因子列表，或线应变百分比范围与点数（至少 4 个点） |
| `--evslc FILE` | 直接拟合已有的 `System.Name.EvsLC` 文件（不运行） |
| `--ngrid N1 N2 N3` / `--ngrid-from OUT` | 固定 `scf.Ngrid` / 从 `.out` 中读取网格 |
| `--workdir DIR` | 输入/输出目录（`eos`） |
| `--prefix NAME` | `System.Name` 前缀（模板的 `System.Name`） |
| `--mpirun` / `--openmx` / `--nt` | 同 `omx_conv_test.py`（`$MPIRUN` 或 `mpirun -np 4` / `$OPENMX` 或 `openmx` / 1） |
| `--data-path DIR` | 覆盖生成输入的 `DATA.PATH`（相对路径会改为绝对路径） |
| `--keep-fileout` | 保留 `level.of.fileout`（默认设为 0） |
| `--dry-run` | 只写输入（此时网格必须已知：`--ngrid`、`--ngrid-from`、模板的 `scf.Ngrid`，或已有的 gridprobe 输出） |
| `--collect` / `--rerun` | 只拟合已有输出 / 已完成的点也重新运行 |
| `--eos bm3\|murnaghan\|vinet` | 状态方程形式（`bm3`，三阶 Birch–Murnaghan） |
| `--a-ref A` | 缩放因子 1.0 对应的晶格常数（Å），报告 a0 = A × s0，s0 = (V0/V_ref)^(1/3) |
| `--a-ref-factor F` | EvsLC 模式：a0 = (F × V0)^(1/3)（fcc 原胞 F = 4，bcc 原胞 F = 2，惯用立方胞 F = 1） |
| `--natoms N` | 每个晶胞的原子数（默认读 `Atoms.Number`），用于输出每原子体积 |
| `--csv` / `--plot` | 默认 `<workdir>/eos.csv`、`<workdir>/eos.png` |
| `--no-plot` / `--show` | 不画图 / 弹出交互窗口 |

**输出**：屏幕上打印 E0（eV/cell）、V0（Å³/cell 及 Å³/atom）、B0（GPa）、B0′、s0、a0 和拟合残差 RMS（meV/cell）；V0 落在采样体积之外或靠近边缘时给出警告，提示扩展或平移缩放范围。CSV 每行为 `label, volume_Ang3, energy_eV, fit_eV, residual_meV`，末尾一行为拟合参数。

**局限**：模板的 `MD.Type` 若是晶胞优化或扫描类（`OptC*`、`RFC*`、`EvsLC` 等）会报错退出；若是原子弛豫（如 `RF`），每个体积都会弛豫内坐标，并取最后一个 Utot。只做各向同性缩放，不适用于需要同时优化 c/a 的体系（请参考第6章的晶胞优化）。

---

## 6. omx_parse_out.py — 输出文件解析

**用途**：解析 OpenMX 的主输出 `System.Name.out`，必要时结合 `.std` 和 `.md`，给出一份文字摘要，或导出 JSON/CSV。其他脚本（`omx_conv_test.py`、`omx_eos_fit.py`）把它当作模块导入。识别的内容（4.0 格式，3.9 同样适用）：

- 回显的输入关键词、截断能与网格（“Required cutoff energy”“Num. of grids of a-, b-, and c-axes”）；
- 每个 MD/优化步的 SCF 历史（`SCF= n NormRD= ... Uele= ...`）；
- `Total energy (Hartree)` 块中的 Utot 及各能量分量、化学势；Cluster 计算的 HOMO 与本征值；
- Mulliken 布居（Up/Down/Sum/Diff，非共线时还有 θ/φ）、总自旋磁矩、偶极矩；
- 原子坐标与力、结构/晶胞优化历史表、耗时（“Elapsed.Time.”“Total Computational Time”）；
- `.std` 中的正常结束标志 “The calculation was normally finished.”，以及 `.out` 不完整时（程序中途崩溃）`<DFT>` 行里的 Uele/dUele/NormRD；`.md` 中每步的 `scf_conv=0/1`。

**SCF 收敛判断**：OpenMX 达到收敛判据或用完 `scf.maxIter` 步都会退出 SCF 循环。某个 SCF 块的步数达到 `scf.maxIter`（从回显输入中读取，缺省按 40）或 `.md` 中对应步为 `scf_conv=0` 时，标记为“未收敛”。

**电荷与磁矩**：Mulliken 电荷 = VPS 价电子数 − Mulliken 布居（Sum 列）；价电子数来自 `--data-path` 或回显输入中的 `DATA.PATH` 下的 VPS 文件，找不到时用 `omx_make_input.py` 的内置表。自旋磁矩为 OpenMX 约定的 Diff 列（符号与电子的物理磁矩相反，见 [第11章](../tutorials/11-spin-noncollinear.md)）。

**依赖**：仅标准库。

**用法示例**：

```bash
python omx_parse_out.py met.out
python omx_parse_out.py si_ecut_*.out --summary-csv ecut_summary.csv
python omx_parse_out.py fe.out --json fe.json --csv-dir fe_tables
python omx_parse_out.py crashed.std                       # .std only (crashed runs)
python omx_parse_out.py nio.out --data-path /opt/openmx4.0/DFT_DATA19 --strict
python omx_parse_out.py *.out -q --strict || echo "some runs failed or did not converge"
```

```python
# as a module (run from this folder or add it to PYTHONPATH)
from omx_parse_out import parse_out
res = parse_out("met.out")
print(res["utot_Ha"], res["scf_converged"])
```

**参数**：

| 参数 | 说明 |
|---|---|
| `files` | 一个或多个 `.out`（或 `.std`）文件（位置参数） |
| `--data-path DIR` | DFT_DATA19 目录，用于读价电子数（默认取回显输入中的 `DATA.PATH`） |
| `--json FILE` | 把全部解析结果写成 JSON（`-` 表示标准输出） |
| `--csv-dir DIR` | 每个文件写 `<name>_scf.csv`、`_energies.csv`、`_mulliken.csv`、`_forces.csv`、`_opt.csv` |
| `--summary-csv FILE` | 每个输入文件一行的汇总表，列为 `file, system_name, finished, scf_converged, scf_nsteps_last, scf_maxIter, utot_Ha, utot_eV, chemical_potential_eV, total_spin_moment_muB, gap_eV, max_force_eV_Ang, ngrid, input_energycutoff, input_kgrid, elapsed_s` |
| `-q/--quiet` | 不打印文字摘要 |
| `--strict` | 只要有一个计算没有正常结束或最后一次 SCF 未收敛，退出码就为 2（适合在作业脚本里判断） |

**局限**：Mulliken 表在文字摘要中最多显示 50 个原子（完整数据见 `--json`/`--csv-dir`）；对 NEGF、Wannier 等专门输出不做解析；不同打印级别下的个别字段可能缺失，此时对应键不出现在结果中。

---

## 7. omx_band_plot.py — 能带图与带隙

**用途**：读取 [第8章](../tutorials/08-band.md) 中 `Band.dispersion on` 产生的 `System.Name.Band`，或 bandgnu13 输出的 `System.Name.BANDDAT1`（自旋极化时还有 `BANDDAT2`），画能带并报告 VBM、CBM 和带隙。

- `.Band` 文件中能量为 Hartree、未平移；脚本用 27.2113845 eV/Ha（与 bandgnu13 相同）换算，并以化学势 μ 为零点。横轴为累计路径长度（Bohr⁻¹，含 2π），段与段之间不插入间距，与 bandgnu13 的算法一致。
- `.Band` 头部的 `spin_flag = 1` 表示共线自旋（每个 k 点先上自旋后下自旋）；`spin_flag = 0` 既可能是非自旋，也可能是非共线/SOC 计算（此时本征值个数为基组数的 2 倍，每条能带只容纳 1 个电子），脚本会打印提示。
- `BANDDAT1/2` 已是 eV 且以化学势为零点；开头 Nkpath−1 个块是分隔竖线，会被跳过。

**标签**：`.Band` 自带 `<Band.kpath` 的标签；`BANDDAT1` 优先读同名 `.GNUBAND` 中的 `set xtics (...)`，否则用 `--dat` 指定的输入文件中的 `Band.Nkpath` / `<Band.kpath`；`--labels` 覆盖一切（连续路径给 Nkpath+1 个，否则给 2×Nkpath 个起止标签）。`g`、`G`、`Gamma`、`GM` 都画成希腊字母 Γ。

**带隙分析**：E − μ ≤ 0 的态算作占据（或用 `--nocc N` 指定每个自旋通道的占据带数，非共线计算中 N = 电子数）。若某条带跨过 E_F，则报告为沿路径金属性。VBM/CBM **只在路径上的 k 点中寻找**，真实带隙还需结合密 k 网格的 DOS 确认。

**依赖**：`pip install numpy matplotlib`。

**用法示例**：

```bash
python omx_band_plot.py si.Band
python omx_band_plot.py si.Band --emin -13 --emax 6 -o si_band.png --csv si_band.csv
python omx_band_plot.py fe.Band --spin both --emin -8 --emax 4
python omx_band_plot.py gaas_soc.Band --compare gaas.Band --compare-label "no SOC" --label "SOC"
python omx_band_plot.py cdia.BANDDAT1                     # labels from cdia.GNUBAND
python omx_band_plot.py cdia.BANDDAT1 --dat Cdia_band.dat
python omx_band_plot.py si.Band --labels "L G X W K G" --ref vbm
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `bandfile` | `System.Name.Band` 或 `System.Name.BANDDAT1`（位置参数） |
| `--dat FILE` | 从该输入文件的 `<Band.kpath` 读标签 |
| `--labels "L G X ..."` | 覆盖标签 |
| `--emin` / `--emax` | 能量范围（−6 / 6 eV） |
| `--ref mu\|vbm` | 能量零点：化学势（默认）或 VBM（仅限有带隙的体系） |
| `--shift` | 所有能带额外平移（0 eV） |
| `--nocc N` | 每个自旋通道的占据带数（代替 E ≤ μ 判据） |
| `--spin both\|up\|down` | 共线自旋时画哪个通道（`both`） |
| `--compare FILE` / `--compare-label` | 第二个 `.Band`/`BANDDAT1`，以虚线叠加（同一路径；路径长度不同时自动缩放横轴）/ 其图例（`compare`） |
| `--label` / `--title` / `--color` | 主数据图例 / 图标题 / 非自旋能带颜色 |
| `--figsize W H` / `--dpi` | 图尺寸（5.0 4.5 英寸）/ 分辨率（300） |
| `-o/--out` | 图片名（`<stem>_band.png`） |
| `--csv FILE` | 导出横坐标（`x_Bohr^-1`）和各能带能量（eV，自旋极化时列名带 `up_`/`dn_` 前缀） |
| `--no-plot` / `--show` | 只打印分析结果 / 弹出交互窗口 |

**局限**：不做能带展开（unfolding，见 [第8章](../tutorials/08-band.md) 的示例）和轨道投影（fat band）；`--compare` 要求两份数据的路径相同。SOC 与拓扑分析见 [第12章](../tutorials/12-soc-topology.md)，Wannier 插值能带对比见 [第14章](../tutorials/14-wannier.md)。

---

## 8. omx_dos_auto.py — 非交互 DosMain 与 PDOS 求和

**用途**：DosMain 运行时会交互式地提问（展宽方法、DOS 还是 PDOS、投影哪些原子），不便于批量处理。本脚本先读 `System.Name.Dos.val` 的头部（mode、Kgrid、Nspin、NonCol、atomnum、`<WhatSpecies`），按 4.0 `DosMain.c` 的提问顺序准备好答案并喂给 DosMain，然后把生成的 PDOS 文件按元素（或物种、原子）与 s/p/d/f 求和、作图。对应 [第9章](../tutorials/09-dos-pdos-cohp-optics.md)。

DosMain 的提问方式取决于 `.Dos.val` 的 mode：mode 1 且 Kgrid ≠ 1 1 1 时会问“Tetrahedron(1), Gaussian Broadening(2)”（选高斯时再问展宽）；mode 1 且 Kgrid = 1 1 1，或 mode 5（DC/Krylov）时自动用高斯展宽，只问展宽；mode 6（NEGF）、7（DosGauss.fileout）不问方法。之后问 “Dos(1) or PDos(2)”，atomnum > 1 时再问原子列表。

**依赖**：`pip install numpy matplotlib`；DosMain 可执行文件（`make DosMain` 或 `omx_build.sh` 编译，可用 `--dosmain` 或环境变量 `DOSMAIN` 指定路径）。

**用法示例**：

```bash
python omx_dos_auto.py si                                   # tetrahedron, total DOS + PDOS of all atoms
python omx_dos_auto.py nio_afm --method gauss --width 0.1 --dat NiO_AFM_pdos.dat
python omx_dos_auto.py mos2 --atoms "1-3" --by atom --resolve-m --emin -8 --emax 5
python omx_dos_auto.py si --no-run --dat Si.dat             # only sum existing PDOS files
python omx_dos_auto.py si --dry-run                         # print the DosMain answers only
python omx_dos_auto.py si --csv si_pdos.csv --no-plot
DOSMAIN=/opt/openmx4.0/work/DosMain python omx_dos_auto.py si
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `prefix` | `System.Name`（可带路径），需要 `<prefix>.Dos.val` 和 `<prefix>.Dos.vec`（位置参数） |
| `--dosmain EXE` | DosMain 可执行文件（`$DOSMAIN` 或 PATH 中的 `DosMain`） |
| `--method tetra\|gauss` | DosMain 询问时选择的方法（`tetra`） |
| `--width` | 高斯展宽参数 a（exp(−(E/a)²)，0.1 eV） |
| `--what both\|dos\|pdos` | 运行 DosMain 计算总 DOS、PDOS 或两者（`both`） |
| `--atoms "1 2 5-8"` | PDOS 投影的原子（默认全部） |
| `--no-run` | 不调用 DosMain，只处理已有文件 |
| `--dry-run` | 只打印将要输入给 DosMain 的答案 |
| `--dat FILE` | 用于把原子对应到元素的输入文件（默认依次找 `<base>.dat`、`<base>.out`） |
| `--by element\|species\|atom` | PDOS 分组方式（`element`） |
| `--resolve-m` | 另外按 m 分量求和（px、dxy……） |
| `--emin` / `--emax` | 作图范围（−10 / 5 eV） |
| `--csv FILE` | 所有求和曲线写到一个 CSV，代替多个 `.dat` 文件 |
| `-o/--out` | 图片名（`<base>_pdos.png`） |
| `--no-plot` / `--show` | 不画图 / 弹出交互窗口 |

**输入/输出**：

- DosMain 在 `.Dos.val` 所在目录中运行，日志写到 `<base>.dosmain_dos.log` / `<base>.dosmain_pdos.log`。它用文件名**第一个点之前**的部分命名输出（如 `si.Dos.val` → `si.DOS.Tetrahedron`、`si.PDOS.Tetrahedron.atom1.p2`），`System.Name` 含点时脚本会提示。
- PDOS 文件 `<base>.PDOS.<Method>.atomN.<l><m>`（DosMain 已对径向多重度求和）。3 列为非自旋（E、DOS、积分）；5 列为自旋（E、DOS_up、−DOS_down、两个积分），下自旋 DOS 以**负值**写出（非共线计算按上/下旋量投影时同样是 5 列）。能量以化学势为零点。目录里有多种方法的 PDOS 文件时，取最新的一组。
- 输出 `<base>_pdos_<group>_<l>.dat`（或 `--csv` 指定的单个 CSV）和图片；屏幕上打印各组在 E_F 以下的积分电子数，以及 ΣPDOS 与总 DOS 的最大偏差（只有投影了全部原子时才应接近 0）。

**局限**：积分从 `Dos.Erange` 的下限开始，与价电子数比较时要考虑能量窗口；元素名取自 VPS 名（如 `Ni_PBE19` → Ni，空原子 `E` → E），找不到输入文件时退回 `.Dos.val` 中的物种编号。COHP、光学性质请用 OpenMX 自带的 cohp 等工具，见 [第9章](../tutorials/09-dos-pdos-cohp-optics.md)；DFT+U 前后 d/f 态的对比见 [第13章](../tutorials/13-dft-u.md)。

---

## 9. omx_cube_tools.py — cube 文件工具箱

**用途**：处理 OpenMX 写出的 Gaussian cube 文件（`tden`/`sden`/`dden`/`v0`/`v1`/`vhart` 以及分子轨道 cube），可代替 OpenMX 自带的 gcube2oned、diff_gcube、add_gcube 等小工具，并补充功函数和 Bader 前处理。对应 [第10章](../tutorials/10-charge-mo-cube.md)、[第16章](../tutorials/16-surface-esm-workfunction.md)。

**OpenMX cube 约定**（按 4.0 `OutData.c` 核对）：长度单位 Bohr；密度单位 e/Bohr³；势（v0、v1、vhart、vxc）单位 Hartree；数据按 a（外层）→ b → c（内层）循环，每行 6 个数，格式 `%13.3E`，即只有 4 位有效数字，所以 v0 的平面平均精度约为 1e-4（相对），足够把功函数算到 0.01 eV 量级，但不要期望更高精度。原子行写的是原子序数（空原子为 0）和净 Mulliken 电荷，而不是核电荷。

> [!WARNING]
> 计算功函数请用 `v0.cube`（或 `v1.cube`，即不含非局域部分的 Kohn-Sham 势）。`vhart.cube` 是**差分电荷密度**产生的 Hartree 势，不是总静电势（见 OpenMX 手册第 50 章），不能直接用来求真空能级。

**子命令与参数**：

| 子命令 | 作用 | 参数（括号内为默认值） |
|---|---|---|
| `info` | 头部、网格、晶胞、最大/最小值与积分 | `cube` |
| `integrate` | 积分（tden 得电子数，sden 得 μB），可只积某轴方向的一段分数坐标 | `cube`、`--axis`（c）、`--range LO HI`（分数坐标 [LO, HI)）、`--abs`（积分绝对值） |
| `avg` | 沿 a/b/c 的平面平均，可加宏观平均 | `cube`、`--axis`（c，也可写 1/2/3）、`--unit none\|eV`（none；eV 时把 Hartree 势乘以 27.2114）、`--macro`（宏观平均窗口，Å；0 = 关闭）、`-o/--out`（`<cube>_avg_<axis>.dat`）、`--plot`、`--show` |
| `workfunc` | 由 v0.cube 和化学势求功函数 Φ = V_vac − μ | `cube`、`--axis`（表面法向，c）、`--mu`（Hartree）与 `--outfile`（从 `System.Name.out` 读 “Chemical potential (Hartree)”）二选一、`--left LO HI` / `--right LO HI`（两侧真空区窗口，Å，从 cube 原点算起；默认取晶胞长度的 5–15% 和 85–95%）、`--dat`（写 z、V、V−μ，eV）、`--plot`、`--show` |
| `diff` | ρ(AB) − ρ(A) − ρ(B) − …，写新 cube，可加平面积分的 Δρ(z) 与电荷位移曲线 ΔQ(z) | `cubes`（先 AB，再依次是要减去的）、`-o/--out`、`--axis`、`--dat`、`--plot`、`--show` |
| `combine` | 任意线性组合 c1·cube1 + c2·cube2 + … | `cubes`、`--coef`（每个 cube 一个系数，必需）、`-o/--out`（必需） |
| `psi2` | 周期体系分子轨道的 \|ψ\|² = re² + im² | `real`（`..._r.cube`）、`imag`（`..._i.cube`，可省略）、`-o/--out`（必需） |
| `bader-prep` | 为 Henkelman 的 bader 程序准备 cube：可把多个 cube 相加，从原子列表中删除空原子（Z = 0），可把小的负密度置零 | `cubes`、`-o/--out`（必需）、`--keep-empty`、`--clip-negative` |

平面平均的坐标是从 cube 原点到“另外两个网格矢量所在平面”的垂直距离，正交晶胞时就是沿该轴的长度。

**依赖**：`pip install numpy`；`--plot` 需要 matplotlib。

**用法示例**：

```bash
python omx_cube_tools.py info h2o.tden.cube
python omx_cube_tools.py integrate h2o.tden.cube
python omx_cube_tools.py integrate fe.sden.cube --axis c --range 0.0 0.5
python omx_cube_tools.py avg al111.v0.cube --axis a --unit eV -o al111_v0_avg.dat --plot avg.png
python omx_cube_tools.py workfunc al111_esm_on1_2e.v0.cube --axis a --outfile al111_esm_on1_2e.out --plot wf.png
python omx_cube_tools.py workfunc slab.v0.cube --axis c --mu -0.1523 --left 2 6 --right 40 44
python omx_cube_tools.py diff co_cu_AB.tden.cube co_cu_A.tden.cube co_cu_B.tden.cube -o cdd.cube --axis c --plot cdd.png
python omx_cube_tools.py combine x.up.cube x.dn.cube --coef 1 1 -o x.tot.cube
python omx_cube_tools.py psi2 si.lumo1_0_0_r.cube si.lumo1_0_0_i.cube -o lumo_k1.cube
python omx_cube_tools.py bader-prep h2o.tden.cube -o h2o_bader.cube
```

**局限**：`diff`/`combine` 要求所有 cube 的网格和晶胞完全一致（AB、A、B 必须用相同晶胞和 `scf.Ngrid` 计算，见 [第10章](../tutorials/10-charge-mo-cube.md) 的 [CO/Cu(111) 示例](../examples/10-co-cu111-cdd/README.md)）；功函数要求真空区足够厚、势已经变平，ESM 或偶极校正的设置见 [第16章](../tutorials/16-surface-esm-workfunction.md)；Berry 相/外电场下的势分析见 [第15章](../tutorials/15-berry-efield.md)。

---

## 10. omx_traj.py — 优化/MD 轨迹分析与导出

**用途**：读取结构优化或分子动力学写出的 `System.Name.md`（OpenMX 以追加方式每步写一帧），打印摘要、作能量/最大力/温度随步数的变化图，并导出 extxyz、ASE `.traj` 或 dpdata 格式（DeePMD 等机器学习势的训练数据）。对应 [第6章](../tutorials/06-optimization-eos-neb.md) 和 [第7章](../tutorials/07-md.md)。

**输入格式**（4.0 `iterout.c`）：每帧第一行为原子数；第二行含 `time=`、`Energy=`（Hartree）、MD 时的 `Temperature=`、`Cell_Vectors=`（9 个数，Å）和 `scf_conv=`；之后每个原子一行：元素、x y z（Å）、力（Hartree/Bohr）、速度（m/s）、净电荷、磁矩、θ、φ。速度只在 NVE / NVT_VS / NVT_VS2 / NVT_VS4 / NVT_NH / NVT_Langevin 时写出，优化和 NPT 积分器中为 0；NVT_VS4 按原子组打印多个 `Temperature[i]=`。列数较少的旧文件尽量读取坐标。

由于 `.md` 是追加写入的，同一个 `System.Name` 重复运行会把第二条轨迹接在第一条后面：脚本发现时间倒退时会警告，`--last-run` 只保留最后一次运行的帧。若存在 `<prefix>.ene`（仅 MD），会把其中的 Utot+Ukc（守恒量）一起作图，并报告每原子的漂移。

**依赖**：`pip install numpy matplotlib`；`--traj` 需要 ase；`--dpdata` 需要 dpdata（它用 `fmt="openmx/md"` 读取 `<prefix>.dat` + `<prefix>.md`）。

**用法示例**：

```bash
python omx_traj.py h2o_rf.md                              # summary + h2o_rf_traj.png
python omx_traj.py h2o_rf --fmax-line 0.0003 --force-unit au
python omx_traj.py H2O_NVE --extxyz H2O_NVE.extxyz --skip 100 --stride 2 --velocities
python omx_traj.py Si64_NVT_NH --traj si64.traj           # needs ASE
python omx_traj.py Si64_NVT_NH --dpdata deepmd_si64 --skip 200 --stride 5
python omx_traj.py Si64_NVT_NH --csv si64_steps.csv --no-plot
python omx_traj.py slab_opt --exclude "1-8" --last-run    # ignore fixed atoms for Fmax
```

**参数**（括号内为默认值）：

| 参数 | 说明 |
|---|---|
| `md` | `System.Name.md` 或 `System.Name`（位置参数） |
| `--skip N` | 丢弃前 N 帧（平衡段）（0） |
| `--stride N` | 导出时每 N 帧取一帧（1） |
| `--last-run` | 只保留最后一次时间重置之后的帧（追加运行时） |
| `--exclude "1-8"` | 计算最大力时忽略的原子（如固定原子） |
| `--force-unit ev\|au` | 最大力的显示单位：eV/Å（默认）或 Hartree/Bohr（与 `MD.Opt.criterion` 直接可比） |
| `--fmax-line` | 在力的子图上画一条判据线（单位同 `--force-unit`） |
| `--extxyz FILE` | 写扩展 XYZ（能量 eV，力 eV/Å） |
| `--velocities` | 在 extxyz 中加入速度（m/s） |
| `--traj FILE` | 写 ASE `.traj` 文件 |
| `--nopbc` | 把结构标记为非周期（Cluster 计算） |
| `--dpdata DIR` | 用 dpdata 导出到该目录 |
| `--dpdata-fmt` | dpdata 输出格式（`deepmd/npy`） |
| `--csv FILE` | 每步一行：step、time、E、dE、Fmax、T |
| `--ene FILE` | 要画的 `<prefix>.ene`（默认自动查找） |
| `-o/--out` | 图片名（`<prefix>_traj.png`） |
| `--no-plot` / `--show` | 不画图 / 弹出交互窗口 |

**输出**：摘要包括帧数、原子数、首末晶胞体积、能量（首/末/最小及所在帧，末减首折合 meV/atom）、首末最大力、温度均值与标准差（MD），以及 `scf_conv=0`（SCF 未收敛）的帧数。导出的单位为能量 eV、力 eV/Å（1 Ha = 27.211386245988 eV，1 Bohr = 0.529177210903 Å）。

**局限**：`--last-run` 不作用于 dpdata 导出（dpdata 自己读取整个文件，帧数不一致时脚本会提示）；不计算 RDF、MSD 等结构/动力学性质（可把 extxyz 交给 ASE、MDAnalysis 等工具，见 [Common/tutorials/05-postprocess-visualization.md](../../Common/tutorials/05-postprocess-visualization.md)）；NEB 的能量曲线请看 [第6章](../tutorials/06-optimization-eos-neb.md)。

---

## 11. omx_negf_plot.py — NEGF 输运后处理

**用途**：处理 [第17章](../tutorials/17-negf-transport.md) NEGF Step 3 的输出，文件格式按 4.0 `TRAN_Main_Analysis.c` / `TRAN_Main_Analysis_NC.c` 核对：

- `System.Name.tran<i2>_<i3>`：`NEGF.tran.Kgrid` 中每个 (k2, k3) 点一个文件，能量以**左电极化学势**为零点。共线/非自旋为 9 列（iw、ReE/ImE（Ha）、ReE/ImE（eV）、ReT_up、ImT_up、ReT_down、ImT_down），非共线为 7 列。非自旋计算中 “down” 列是 “up” 的拷贝（每个自旋通道的 T），因此共线情形总透射为 T_up + T_down，非共线为 T。文件覆盖完整 k 网格（k 与 −k），直接对所有文件平均就是正确的 k 平均。本征通道输出（`.traneval*`、`.tranevec*`、`.tranec*`）会被自动跳过。
- `System.Name.current`：偏压（V）与平均电流（A，共线时分上/下自旋）及各 k 点的电流。
- `System.Name.conductance`：平均电导（以每自旋通道的 G0 = e²/h 为单位，非共线时一个值）及各 k 点的电导。

**trans 子命令用到的公式**（E 以 eV 计）：

$$L_n=\int T_{\mathrm{tot}}(E)\,(E-\mu)^n\left(-\frac{\partial f}{\partial E}\right)dE,\qquad G=\frac{e^2}{h}L_0,\qquad S=-\frac{L_1}{eTL_0},\qquad \kappa_{\mathrm{el}}=\frac{1}{hT}\left(L_2-\frac{L_1^2}{L_0}\right)$$

$$I=\frac{e}{h}\int T_{\mathrm{tot}}(E)\,[f_L(E)-f_R(E)]\,dE$$

`--temp 0` 时电导直接取 (e²/h)·T_tot(μ)。由透射谱算电流只有在 T(E) 来自**同一偏压**的计算时才有意义。

**依赖**：`pip install numpy matplotlib`。

**用法示例**：

```bash
python omx_negf_plot.py trans negf-chain                       # all negf-chain.tran*_* files
python omx_negf_plot.py trans negf-graphene.tran*_* --temp 300 --mu -0.2 0 0.2 -o T.png --csv T.csv
python omx_negf_plot.py trans negf-chain --total --log --emin -3 --emax 3
python omx_negf_plot.py kmap negf-graphene --energy 0.0 -o Tk.png
python omx_negf_plot.py cond negf-graphene.conductance -o Gk.png
python omx_negf_plot.py iv negf-chain-V*.current -o iv.png --csv iv.csv
```

**子命令与参数**（括号内为默认值；每个子命令都有 `-o/--out`、`--no-plot`、`--show`）：

| 子命令 | 作用 | 参数 |
|---|---|---|
| `trans` | k 平均透射谱，Landauer 电导、Seebeck 系数、电子热导及由 T(E) 积分的电流 | `files`（`System.Name` 前缀或 `.tran*_*` 文件）、`--temp`（300 K；0 = 零温电导）、`--mu`（相对 μ_L 的化学势，可给多个，0 eV）、`--total`（画自旋求和后的 T）、`--emin`/`--emax`、`--log`（透射取对数坐标）、`--csv`（E、T_up、T_down、T_total）；默认图 `transmission.png` |
| `kmap` | 某一能量下的 k 分辨透射（k 网格为一维时画 T(k, E)） | `files`、`--energy`（E − μ_L，0 eV）、`--csv`（i2、i3、k2、k3、T_total）；默认图 `transmission_k.png` |
| `cond` | 读 `.conductance`，画 k 分辨电导 | `file`；默认图 `conductance_k.png` |
| `iv` | 由一组偏压的 `.current` 文件画 I–V 曲线和微分电导 dI/dV | `files`（可用通配符）、`--csv`（V、I_up、I_down、I_total、dI/dV）；默认图 `iv.png` |

**输出说明**：`trans` 对每个 μ 打印 T_tot(μ)、G（e²/h 与 μS）、S（μV/K）、κ_el（pW/K）和 Lorenz 比 L/L0；μ ± 10 k_BT 超出 T(E) 的能量范围时会警告。`iv` 中 I_tot = I_up + I_down（非自旋时 OpenMX 在两列写相同的值），偏压符号按文件头中的 μ_L − μ_R。k 点文件数与索引范围不符时会提示文件可能缺失。

**局限**：不计算功率因子或 ZT（需要声子热导）；不做本征通道可视化；热电量的精度受 Step 3 能量网格（`NEGF.tran.energyrange`、`NEGF.tran.energydiv`）密度的限制。

---

## 12. job_omx_slurm.sh / job_omx_pbs.sh — 作业模板

**用途**：MPI × OpenMP 混合并行运行 OpenMX 的作业模板，按手册第 28 章的建议设置：总核数 = MPI 进程数 × 线程数，`OMP_NUM_THREADS` 与 `-nt` 保持一致，`ulimit -s unlimited` 加 `OMP_STACKSIZE`（默认 512M）避免大体系段错误，同时设置 `OMP_PLACES=cores`、`OMP_PROC_BIND=close`、`MKL_NUM_THREADS`、`OPENBLAS_NUM_THREADS=1`、`I_MPI_PIN_DOMAIN=omp`。运行前检查：输入文件和可执行文件是否存在、输入是否为 Windows（CRLF）换行、`DATA.PATH` 下是否有 `PAO/` 和 `VPS/`、能否读到 `System.Name`。输出用 `tee` 写到 `System.Name.std`。

**成功判据**：`System.Name.std` 中出现 “The calculation was normally finished.”（mpirun 的退出码并不可靠），并比较 `.out` 中最后一个 SCF 块的步数与 `scf.maxIter`。退出码：0 = 正常结束且最后一次 SCF 收敛，1 = 运行失败（打印 `.std` 最后 30 行），2 = SCF 用满 `scf.maxIter`（多半未收敛）。正常结束时打印 Utot、最后 SCF 块步数和耗时。重启目录 `System.Name_rst/` 默认保留。

**依赖**：bash、awk、OpenMX 可执行文件及其 MPI 运行库；分别需要 Slurm 或 PBS Pro / OpenPBS / Torque。使用前先编辑头部的 `#SBATCH` / `#PBS` 资源行，以及“user-editable variables”中的 `module load` / `source setvars.sh` 环境行。

| 变量 | Slurm 模板 | PBS 模板 |
|---|---|---|
| `INPUT` | 第一个命令行参数，或环境变量 `INPUT`（`input.dat`） | `qsub -v INPUT=...`（`input.dat`） |
| `OPENMX` | `openmx`（不在 PATH 中时写完整路径） | 同左 |
| `OMX_DATA` | 设置后用 `sed -i` 把输入中已有的 `DATA.PATH` 行**原地改写**为该目录 | 同左 |
| 启动方式 | `LAUNCH=mpirun`（默认）或 `srun`（需 MPI 支持 Slurm PMI；会设置 `SRUN_CPUS_PER_TASK`，因为 Slurm ≥ 22.05 的 srun 不再继承 `--cpus-per-task`）；`MPI_EXTRA` 附加选项 | `mpirun` + 由 `$PBS_NODEFILE` 生成的 machine file；`MPI_EXTRA` 附加选项 |
| MPI 进程数 | `SLURM_NTASKS`；Slurm 之外用 `NP`（4） | 默认为 machine file 行数，`NP` 可覆盖；PBS 之外为 4 |
| 线程数 | `SLURM_CPUS_PER_TASK`；Slurm 之外用 `NT`（1） | `NT`（默认取 `OMP_NUM_THREADS`，否则 1） |
| 调度器差异 | — | `PBS_FLAVOR=pro`（默认，`select=...:mpiprocs=...:ompthreads=...`，nodefile 每个 MPI 进程一行）或 `torque`（`nodes=N:ppn=P`，nodefile 每核一行，脚本按每节点 核数/NT 生成 MPI 槽位） |
| 其他 | `CLEAN_RST=1`：成功后删除 `System.Name_rst/`；`USE_DEB=1`：把官方 deb 包的 `/opt/openmx4.0/runtime/lib` 加入 `LD_LIBRARY_PATH`（解决 “libmpifort.so.12: cannot open shared object file”） | 同左 |

```bash
sbatch job_omx_slurm.sh input.dat
sbatch --nodes=2 --ntasks-per-node=32 --cpus-per-task=2 job_omx_slurm.sh si64.dat
INPUT=fe.dat OMX_DATA=/opt/openmx4.0/DFT_DATA19 sbatch job_omx_slurm.sh
NP=4 NT=1 bash job_omx_slurm.sh Methane.dat           # also works without Slurm (workstation)
qsub -v INPUT=input.dat job_omx_pbs.sh
qsub -v INPUT=si64.dat,OMX_DATA=/opt/openmx4.0/DFT_DATA19 job_omx_pbs.sh
qsub -l select=2:ncpus=64:mpiprocs=32:ompthreads=2 -v INPUT=si64.dat job_omx_pbs.sh
bash job_omx_slurm.sh -h                              # print the usage text
```

> [!TIP]
> 进程数与线程数的搭配（大体系少进程多线程、O(N) 方法的并行方式）和常见运行错误见 [第18章](../tutorials/18-parallel-on-errors.md)；带运行检查 `check_run.sh` 和并行基准 `bench_np_nt.sh` 的更多模板见 [examples/18-job-scripts](../examples/18-job-scripts/README.md)；调度系统基础见 [Common/tutorials/02-job-scheduler.md](../../Common/tutorials/02-job-scheduler.md)。

**局限**：一个作业只运行一个输入；NEGF 的三步、能带/DOS 的后处理等多步流程需要自己在脚本末尾追加命令。

---

## 典型工作流

```bash
# 0) build OpenMX once
bash omx_build.sh --toolchain gcc --prefix $HOME/software --np 4
export PATH=$HOME/software/openmx4.0/work:$PATH
DATA=$HOME/software/openmx4.0/DFT_DATA19
# 1) draft an input from a structure file and review it
python omx_make_input.py Si.cif -o si.dat --data-path $DATA --kspacing 0.25
# 2) converge cutoff and k-mesh
python omx_conv_test.py ecut si.dat --values 150 180 220 260 300 350
python omx_conv_test.py kgrid si.dat --values 6 8 10 12 14
# 3) lattice constant from the equation of state
python omx_eos_fit.py si.dat --range -3 3 7 --a-ref 5.43
# 4) production run + health check
sbatch job_omx_slurm.sh si.dat
python omx_parse_out.py si.out --strict
# 5) band structure and DOS/PDOS (inputs with Band.kpath / Dos.fileout on)
python omx_make_input.py Si.cif -o si_band.dat --data-path $DATA --kpath --dos
python omx_band_plot.py si_band.Band --emin -13 --emax 6
python omx_dos_auto.py si_band --by element
```

## 相关资源

- OpenMX 官方手册（4.0）：https://www.openmx-square.org/openmx_man4.0/index.html
- 关键词列表：https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html
- Database (2019) 推荐基组表：https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html
- 下载页（4.0、patch 4.0.1、deb、Docker）：https://www.openmx-square.org/download.html
- 其他社区工具（只链接，未复制代码）：
  - [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文教程站，含 GCC/Intel 编译脚本与 Python 后处理示例；
  - [mfukudaQED/openmx_installer](https://github.com/mfukudaQED/openmx_installer)：下载、打补丁、改 makefile 的一键编译脚本，可与 `omx_build.sh` 对照；
  - [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx)、[pohao82/vasp-omx-moment-setter](https://github.com/pohao82/vasp-omx-moment-setter)：POSCAR 转 `.dat` 与初始磁矩设置的另一种实现；
  - [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS)：浏览器中查看能带 + DOS/PDOS 的单文件工具；
  - [TnakaTnakaTnaka/TEP_tools](https://github.com/TnakaTnakaTnaka/TEP_tools)：NEGF 透射谱的热电后处理（功率因子、ZT）；
  - [mailhexu/TB2J](https://github.com/mailhexu/TB2J)：由 OpenMX 结果计算磁交换参数。
- 结构建模、ASE/pymatgen 的通用用法见 [Common/tutorials/04-ase-pymatgen.md](../../Common/tutorials/04-ase-pymatgen.md)。

> [!CAUTION]
> 这些脚本按 OpenMX 4.0 的源码和手册中的输出格式编写，但无法覆盖所有打印级别、版本和体系组合，发布前也没有在真实计算上逐一测试。用在新体系上时，请先拿一个小算例与手工 `grep` 的结果对照；发现解析错误时，欢迎在 GitHub Issues 中附上输出片段反馈。

[返回目录](../README.md)
