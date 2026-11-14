# VASP 辅助脚本（VASP helper scripts）

本目录收录 VASP 计算的前处理、作业提交、收敛测试和后处理脚本，覆盖作业模板（Slurm/PBS）、输入生成（INCAR/KPOINTS/KPOINTS_OPT/POTCAR）、作业健康检查、能带/DOS 作图、功函数、差分电荷、NEB、phonopy 声子、弹性常数、光学性质和 AIMD/MLFF 分析。所有 Python 脚本都支持 `-h` 查看完整参数，Bash 脚本用 `bash <脚本> -h` 查看帮助。下面的参数都取自脚本源码中的 argparse 定义或 usage 文本。

> [!NOTE]
> 所有脚本都按 VASP 6.4–6.6 的真实输出格式编写，但没有在每个版本、每种编译方式上逐一验证。第一次用某个脚本时，建议先在小体系上和手工 `grep` 的结果核对一遍。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [job_vasp_slurm.sh](job_vasp_slurm.sh) | Slurm 作业模板：自动选 std/gam/ncl，MPI×OpenMP，运行前后检查 | bash, Slurm | [第1章](../tutorials/01-setup-compile.md)、[第4章](../tutorials/04-relax-eos.md) |
| [job_vasp_pbs.sh](job_vasp_pbs.sh) | PBS Pro / Torque 作业模板，功能同上 | bash, PBS | [第1章](../tutorials/01-setup-compile.md) |
| [vasp_inputgen.py](vasp_inputgen.py) | 按任务类型由 POSCAR 生成带注释的 INCAR/KPOINTS/KPOINTS_OPT/POTCAR | numpy；可选 ase、seekpath、spglib | [第2章](../tutorials/02-input-output-files.md)、[第4章](../tutorials/04-relax-eos.md)、[第5章](../tutorials/05-band-dos.md) |
| [vasp_incar_writer.py](vasp_incar_writer.py) | 为非磁、磁性、DFT+U 静态（DOS）计算写 INCAR.static | 仅标准库 | [第7章](../tutorials/07-magnetism-u-soc.md) |
| [check_vasp.py](check_vasp.py) | 作业健康检查：收敛、受力、磁矩、耗时、报错诊断 | 仅标准库（有 numpy 时会用） | [第2章](../tutorials/02-input-output-files.md)、[第18章](../tutorials/18-postprocess-troubleshooting.md) |
| [conv_test.sh](conv_test.sh) | ENCUT/k 点/KSPACING/SIGMA/任意 INCAR 参数的收敛测试驱动 | bash ≥ 4, awk | [第3章](../tutorials/03-convergence.md) |
| [conv_analyze.py](conv_analyze.py) | 分析收敛测试：dE/atom、受力、T*S 表格与曲线 | numpy, matplotlib | [第3章](../tutorials/03-convergence.md) |
| [plot_band.py](plot_band.py) | 能带图（含 KPOINTS_OPT、HSE 零权重、fat band），带隙/VBM/CBM | numpy, matplotlib | [第5章](../tutorials/05-band-dos.md)、[第6章](../tutorials/06-hybrid-metagga.md)、[第17章](../tutorials/17-advanced-electronic.md) |
| [plot_dos.py](plot_dos.py) | 总 DOS/PDOS 作图，d 带中心与带宽 | numpy, matplotlib | [第5章](../tutorials/05-band-dos.md)、[第10章](../tutorials/10-adsorption-catalysis.md) |
| [locpot_workfunction.py](locpot_workfunction.py) | LOCPOT 平面/宏观平均，真空能级与功函数 | numpy, matplotlib | [第9章](../tutorials/09-surface-workfunction.md) |
| [chgdiff.py](chgdiff.py) | 差分电荷密度、平面平均电荷转移、AECCAR 求和（Bader） | numpy, matplotlib | [第16章](../tutorials/16-charge-bonding.md)、[第10章](../tutorials/10-adsorption-catalysis.md) |
| [neb_tools.py](neb_tools.py) | NEB 插值（线性/IDPP）、原子顺序检查、能垒收集与样条拟合、XYZ 动画 | numpy, matplotlib, scipy；可选 ase | [第11章](../tutorials/11-neb-dimer.md) |
| [phonopy_workflow.sh](phonopy_workflow.sh) | phonopy 有限位移/DFPT 全流程：建目录、提交、检查、收集、作图、虚频检测 | bash ≥ 4, phonopy, seekpath | [第12章](../tutorials/12-phonon.md) |
| [elastic_post.py](elastic_post.py) | 弹性张量：kBar→GPa、Voigt 重排、Born 稳定性、VRH 模量、德拜温度、2D | numpy；可选 ase、spglib | [第13章](../tutorials/13-elastic.md) |
| [optics_post.py](optics_post.py) | 介电函数 → n、k、吸收系数、反射率、能量损失、Tauc 图；静态介电张量与 Born 电荷 | numpy, matplotlib | [第14章](../tutorials/14-optics-dielectric.md) |
| [md_analysis.py](md_analysis.py) | AIMD/MLFF 分析：温度/能量曲线、RDF、MSD 与扩散系数、ML_LOGFILE 误差 | numpy, matplotlib；可选 ase | [第15章](../tutorials/15-md-mlff.md) |

## 安装依赖

```bash
# 所有 Python 脚本都要求 Python >= 3.9
pip install numpy matplotlib scipy        # 大多数作图/分析脚本
pip install ase spglib seekpath           # 可选：IDPP 插值、k 路径、晶系识别、原子质量
pip install phonopy                       # phonopy_workflow.sh（phonopy 自带 phonopy-bandplot 等命令）
```

Bash 脚本需要 bash ≥ 4 和 awk。`conv_test.sh`、`phonopy_workflow.sh` 的提交功能需要 `sbatch` 或 `qsub`。

> [!TIP]
> 可以把本目录加入 `PATH`（`export PATH=$HOME/DFT-Tutorials/VASP/scripts:$PATH`），之后在任何计算目录下直接调用 `check_vasp.py`、`plot_band.py` 等。Python 脚本的首行是 `#!/usr/bin/env python3`。

一个典型的工作流如下：

```text
vasp_inputgen.py  ->  conv_test.sh + conv_analyze.py  ->  job_vasp_slurm.sh / job_vasp_pbs.sh
      |                                                           |
      v                                                           v
 （按任务生成输入）                                    check_vasp.py（检查收敛与报错）
                                                                  |
              plot_band.py / plot_dos.py / locpot_workfunction.py / chgdiff.py /
              neb_tools.py / phonopy_workflow.sh / elastic_post.py / optics_post.py / md_analysis.py
```

---

## 1. job_vasp_slurm.sh — Slurm 作业模板

**用途**：VASP 6.x 的 Slurm 提交脚本。它会自动选择 `vasp_std`/`vasp_gam`/`vasp_ncl`，按 MPI×OpenMP 设置环境变量，运行前检查输入文件，运行后检查是否正常结束、是否收敛以及有没有常见报错。`conv_test.sh` 和 `phonopy_workflow.sh` 默认把它复制到每个子目录，文件名改为 `job.sh`。

**依赖**：bash、Slurm、VASP 可执行文件。可选 `check_vasp.py`（运行后生成详细报告）。

**用法**：

```bash
sbatch job_vasp_slurm.sh                   # VARIANT=auto，自动选择可执行文件
sbatch job_vasp_slurm.sh ncl               # 强制 vasp_ncl（可选 std | gam | ncl | auto）
sbatch --export=ALL,VARIANT=gam,LAUNCHER=srun job_vasp_slurm.sh
sbatch -N 2 --ntasks-per-node=32 --cpus-per-task=2 job_vasp_slurm.sh   # MPI x OpenMP
sbatch --export=ALL,RESTART=1 job_vasp_slurm.sh   # 续算结构优化：CONTCAR -> POSCAR
bash job_vasp_slurm.sh -h                  # 打印帮助
```

**可调变量**（直接改文件，或用 `sbatch --export=ALL,VAR=...` 传入）：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `VARIANT` | `auto` | `auto`/`std`/`gam`/`ncl`；也可以作为第一个命令行参数传入 |
| `VASP_BIN_DIR` | 空（使用 `PATH`） | VASP 可执行文件所在目录，例如 `$HOME/software/vasp.6.5.1/bin` |
| `LAUNCHER` | `mpirun` | `mpirun`（Intel MPI / Open MPI）或 `srun` |
| `MPI_EXTRA` | 空 | 附加给启动器的选项，例如 `--report-bindings` |
| `STDOUT_FILE` | `vasp.out` | VASP 标准输出文件 |
| `RESTART` | `0` | `1`：先把原 POSCAR 备份成 `POSCAR.N`，再把非空的 CONTCAR 复制为 POSCAR |
| `CHECK_VASP` | 空 | `check_vasp.py` 的路径；设置后运行结束会自动调用它 |
| `PYTHON` | `python3` | 调用 `check_vasp.py` 用的解释器 |
| `OMP_STACKSIZE` | `512m` | OpenMP 线程栈大小 |

计算资源由文件头部的 `#SBATCH` 行决定：`--nodes`、`--ntasks-per-node`（MPI 进程数）、`--cpus-per-task`（每个进程的 OpenMP 线程数）和 `--time`。`--partition`/`--account`/`--mem`/`--gres` 默认注释掉（`##SBATCH`），需要时删掉一个 `#` 再按集群修改。环境加载部分（`module load ...`）也要按你的安装方式改，只保留一个方案。

**主要行为**：

- **自动选可执行文件**：INCAR 中 `LSORBIT` 或 `LNONCOLLINEAR = .TRUE.` 时用 `vasp_ncl`；没有 `KSPACING`，且 KPOINTS 是 Gamma 中心的 1×1×1 自动网格时用 `vasp_gam`；其余情况用 `vasp_std`。
- **MPI×OpenMP**：`OMP_NUM_THREADS` 取自 `SLURM_CPUS_PER_TASK`，同时设置 `OMP_STACKSIZE=512m`、`OMP_PLACES=cores`、`OMP_PROC_BIND=close`（依据 VASP wiki 的 [Combining MPI and OpenMP](https://www.vasp.at/wiki/index.php/Combining_MPI_and_OpenMP) 页面），并设置 `SRUN_CPUS_PER_TASK`。线程绑定方面，Open MPI 用 `--map-by slot:PE=n --bind-to core`，Intel MPI 用 `I_MPI_PIN_DOMAIN=omp`，`srun` 用 `--cpu-bind=cores`。每个进程多于 1 个线程时，脚本会提示 VASP 内部使用 NCORE = 1（同一 wiki 页面）。脚本还会执行 `ulimit -s unlimited`。
- **运行前检查**：INCAR/POSCAR/POTCAR 是否缺失或为空；是否既没有 KPOINTS 也没有 `KSPACING`；INCAR 中是否有 CR（Windows 换行）或 TAB；POSCAR 第 6 行的元素顺序是否与 POTCAR 中 `VRHFIN` 的顺序一致；`NCORE`、`KPAR` 是否整除 MPI 进程数。
- **运行后检查**：OUTCAR 是否有 `General timing and accounting`；`NSW > 0` 且 `IBRION = 1–3` 时是否出现 `reached required accuracy`；OSZICAR 中有多少离子步的 SCF 达到了 NELM；`vasp.out` 和 `slurm-<jobid>.out` 中是否有已知报错（ZBRENT、BRMIX、EDDDAV、Sub-Space-Matrix、BZINTS、`DUE TO TIME LIMIT`、oom-kill 等）。

**输出与退出码**：VASP 的标准输出写入 `vasp.out`，检查结果写入 `slurm-<jobid>.out`。退出码 0 表示正常结束，1 表示有问题（崩溃、超时、结构优化未收敛、VASP 返回非零值）。

**局限**：

- `vasp_gam` 的自动判断只识别"第 2 行为 0、第 3 行以 G 开头、网格为 1 1 1"的自动网格，显式列出的 Γ 点不会被识别为 Gamma-only。
- 元素顺序检查需要 VASP 5 格式的 POSCAR（第 6 行是元素名）。
- 线程绑定按 `mpirun --version` 的输出区分 Open MPI 和 Intel MPI，其他 MPI 实现要自己通过 `MPI_EXTRA` 设置。
- GPU（OpenACC）版本只在注释里给了 `--gres` 示例，绑定方式请按集群文档调整。

---

## 2. job_vasp_pbs.sh — PBS Pro / Torque 作业模板

**用途**：PBS 版的 VASP 作业模板，自动选择可执行文件、运行前检查、运行后检查都和 Slurm 版相同。默认启用 PBS Pro 的 `#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1`；Torque 用户把它换成注释里给出的 `#PBS -l nodes=1:ppn=64`。

**依赖**：bash、PBS Pro 或 Torque、VASP 可执行文件。

**用法**：

```bash
qsub job_vasp_pbs.sh                                   # 自动选择可执行文件
qsub -v VARIANT=ncl job_vasp_pbs.sh                    # std | gam | ncl | auto
qsub -v NTHREADS=2 -l select=2:ncpus=64:mpiprocs=32:ompthreads=2 job_vasp_pbs.sh
qsub -v RESTART=1 job_vasp_pbs.sh                      # 续算：CONTCAR -> POSCAR
qsub -v VARIANT=ncl,NTHREADS=2,RESTART=1 job_vasp_pbs.sh
bash job_vasp_pbs.sh -h
```

> [!IMPORTANT]
> `qsub` 不会把命令行参数传给脚本，所以 `qsub job_vasp_pbs.sh ncl` 不起作用，变量一律用 `-v VAR=value` 传入。

**可调变量**：`VARIANT`、`NTHREADS`（默认 1）、`VASP_BIN_DIR`、`LAUNCHER`（默认 `mpirun`）、`USE_HOSTFILE`（默认 1，即向 mpirun 传 `-machinefile $PBS_NODEFILE`；MPI 没有编译 PBS tm 支持时需要它）、`MPI_EXTRA`、`STDOUT_FILE`、`RESTART`、`CHECK_VASP`、`PYTHON`、`OMP_STACKSIZE`，含义同 Slurm 版。

**与 Slurm 版的区别**：

- MPI 进程数取 `$PBS_NODEFILE` 的行数（PBS Pro 由每个 chunk 的 `mpiprocs` 决定，Torque 由 `ppn` 决定）。
- PBS Pro 会按 `ompthreads` 设置 `OMP_NUM_THREADS`，而 `ompthreads` 默认等于 `ncpus`。所以脚本根据 `NTHREADS` 显式设置 `OMP_NUM_THREADS`，请保证 `mpiprocs × NTHREADS = ncpus`。
- 使用 Open MPI 时，OpenMP 相关环境变量通过 `-x` 传给各个进程。

**局限**：同 Slurm 版。另外，不同站点的 PBS 配置（队列名、`place=` 设置）差异较大，提交前请对照本集群的使用手册。

---

## 3. vasp_inputgen.py — 按任务生成输入文件

**用途**：读取 POSCAR，为指定任务生成带逐行注释的 INCAR、按 k 点间距生成 KPOINTS、生成能带用的线模式 KPOINTS 或 KPOINTS_OPT，并从 `$VASP_PP_PATH` 拼接 POTCAR。仓库本身不包含任何 POTCAR。

支持的任务：`relax`、`scf`、`band`、`dos`、`hse`、`dftu`、`soc`、`neb`、`phonon`、`md`、`mlff`、`optics`、`elastic`。

**依赖**：`pip install numpy`。可选：`ase`（k 路径与 Bravais 格子识别）、`seekpath`（`--path-backend seekpath`）、`spglib`（晶系识别）。

**用法**：

```bash
python vasp_inputgen.py relax                          # 读取当前目录的 POSCAR
python vasp_inputgen.py scf -p POSCAR -o scf --kspacing 0.2
python vasp_inputgen.py band -o band                   # ICHARG=11 + 线模式 KPOINTS
python vasp_inputgen.py band --kpoints-opt -o band     # SCF + KPOINTS_OPT（VASP >= 6.3）
python vasp_inputgen.py dos --ispin 2 --magmom "2*3.0 4*0"
python vasp_inputgen.py hse --band --ibzkpt ../pbe/IBZKPT   # HSE 能带：零权重 k 点
python vasp_inputgen.py dftu --hubbard "Ni=6.2" --magmom "2 -2 2*0"
python vasp_inputgen.py soc --saxis "0 0 1" --magmom "0 0 2.2"
python vasp_inputgen.py neb --images 5 -p 00/POSCAR
python vasp_inputgen.py md --temp 1000 --nsw 5000 --ensemble nvt
python vasp_inputgen.py mlff --ml-mode train --temp 600
python vasp_inputgen.py optics --optics-mode loptics
python vasp_inputgen.py relax --relax-2d-cell --vdw d3bj   # 二维材料只优化面内晶格
python vasp_inputgen.py relax --stdout                 # 只打印到屏幕，不写文件
```

**主要参数**：

| 类别 | 参数 |
|---|---|
| 文件 | `-p/--poscar`（默认 POSCAR）、`-o/--outdir`（默认 `.`）、`--elements`（VASP 4 格式的 POSCAR 用）、`-f/--force`、`--stdout`、`--lang {en,zh,none}`（INCAR 注释语言，默认 `en`，纯 ASCII）、`--system` |
| 电子结构 | `--encut`（默认取 `--encut-factor` × 最大 ENMAX，系数默认 1.3）、`--ediff`、`--functional {PBE,PBESOL,RPBE,REVPBE,PW91,AM05,SCAN,R2SCAN}`、`--metal`、`--ismear`、`--sigma`、`--lreal`、`--vdw {d2,d3,d3bj,d4,sd3,ts,mbd}`、`--ispin`、`--magmom`、`--nupdown`、`--hubbard "Ni=6.2,Fe=5.3"`、`--saxis`、`--soc-from-chgcar`、`--nkred`、`--bandgap {COMPACT,WEIGHT,KPOINT}`、`--nbands`、`--write-wavecar` |
| 离子/MD/NEB | `--isif`、`--relax-2d-cell`、`--nsw`、`--ediffg`、`--images`（默认 5）、`--iopt {1,2,3,4,7}`（默认 3）、`--climb/--no-climb`、`--no-vtst`、`--temp`、`--temp-end`、`--ensemble {nvt,npt,nve}`、`--potim`、`--pstress`、`--ml-mode {train,select,refit,run}`、`--dfpt`、`--born`、`--optics-mode {loptics,lepsilon,lcalceps}`、`--ionic` |
| k 点 | `--kspacing`（单位 1/Å，含 2π，与 VASP 的 KSPACING 约定相同）、`--kmesh "n1 n2 n3"`、`--kmode {auto,gamma,mp}`、`--kspacing-in-incar`、`--vacuum`（默认 8 Å）、`--band`、`--kpoints-opt`、`--ibzkpt`、`--nkpath`（默认 40）、`--path-backend {ase,seekpath}`、`--symprec` |
| POTCAR | `--potcar {auto,skip}`、`--pp-path`（默认 `$VASP_PP_PATH`）、`--pp-dir`、`--pp-flavor {recommended,legacy,gw}`、`--potcar-map "Fe=Fe,O=O_h"` |

**输出**：输出目录下的 INCAR、KPOINTS（或 KPOINTS_OPT），以及 POTCAR（`--potcar skip` 时不生成）。六方/三方晶胞、二维 slab（真空方向只取 1 个 k 点）和分子会自动使用 Gamma 中心网格。

**局限**：

- 生成的是起点参数，ENCUT、k 点等仍需按[第3章](../tutorials/03-convergence.md)做收敛测试。
- POTCAR 按脚本内置的推荐赝势表选择（`--pp-flavor legacy` 用旧表），请对照 VASP wiki 的 [Recommended pseudopotentials](https://www.vasp.at/wiki/index.php/Recommended_pseudopotentials) 和课题组习惯核实。
- `--vdw d4`/`sd3` 需要 VASP 编译时链接了相应的外部库。
- seekpath 给出的路径对应 seekpath 的标准原胞，可能和输入的 POSCAR 取向不同。

---

## 4. vasp_incar_writer.py — 静态计算 INCAR.static

**用途**：为常见的过渡金属氧化物等体系写静态（DOS）计算用的 `INCAR.static`，分为非磁、磁性（给 `-m`）和 DFT+U（给 `-hu`）三种情况。脚本只写电子结构相关的标签。默认值沿用原始脚本：PBEsol（`GGA = PS`）、`ENCUT = 600`、`EDIFF = 1e-7`、`ISMEAR = -5`、`LORBIT = 11`、`NEDOS = 3001`，以及 `ISTART = 1`/`ICHARG = 1`（接着结构优化继续算）。

**依赖**：仅 Python 标准库。

**用法**：

```bash
python vasp_incar_writer.py                          # 非磁，PBEsol
python vasp_incar_writer.py -m 2*0_2_-2_6*0          # 自旋极化，MAGMOM 用 '_' 分隔
python vasp_incar_writer.py -m 2_-2_2*0 -hu          # + DFT+U，使用内置 U 值表
python vasp_incar_writer.py -m 2_-2_2*0 -hu Ni=6.2,Mn=3.9
python vasp_incar_writer.py -xc PBE --encut 520 --ncore 8 -o INCAR
python vasp_incar_writer.py --scratch -f             # ISTART=0, ICHARG=2，覆盖已有文件
```

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-name` | POSCAR 第一行 | SYSTEM 名称 |
| `-m` | 无（ISPIN = 1） | MAGMOM，用 `_` 或空格分隔，支持 `N*value`；脚本会核对长度与原子数是否一致 |
| `-hu` | `False` | 只写 `-hu` 或 `-hu True` 时用内置 Dudarev U 值（Co 3.32、Cr 3.7、Fe 5.3、Mn 3.9、Mo 4.38、Ni 6.2、V 3.25、W 6.2 eV）；`-hu Ni=6.2,Mn=3.9` 显式指定 |
| `-p/--poscar` | `POSCAR` | 结构文件；LDAUL/LDAUU/LDAUJ 按其中的元素行逐元素生成 |
| `-o/--output` | `INCAR.static` | 输出文件 |
| `-xc/--functional` | `PBEsol` | `PBESOL`、`PBE`、`RPBE`、`REVPBE`、`PW91`、`AM05`、`R2SCAN`、`SCAN`（meta-GGA 会同时写 `LMIXTAU = .TRUE.`） |
| `--encut` / `--prec` / `--ediff` / `--nelm` | 600 / Accurate / 1e-7 / 40 | |
| `--nedos` / `--emin` / `--emax` | 3001 / -12 / 12 | `--no-window` 不写 EMIN/EMAX |
| `--ncore` / `--kpar` | 4 / 不写 | `--ncore 0` 表示不写 NCORE |
| `--scratch` | 关 | 从头算（ISTART = 0，ICHARG = 2） |
| `-f/--force` | 关 | 覆盖已有输出文件 |

**局限**：内置 U 值只是常用的文献参考值，需要按体系和所研究的性质确定（见[第7章](../tutorials/07-magnetism-u-soc.md)关于线性响应 U 的部分）。`ISMEAR = -5` 要求 k 点足够多（至少 4 个不可约 k 点），否则 VASP 会报四面体方法相关的错误。

---

## 5. check_vasp.py — 作业健康检查

**用途**：检查正在运行或已经结束的 VASP 作业，报告以下内容：

- VASP 版本、任务类型（IBRION/ISIF/NSW）、离子步数，以及每个离子步的电子步数（达到 NELM 的会标为未收敛）。
- 是否出现 `reached required accuracy`、NSW 是否用完、作业是已结束还是仍在运行。
- 最后的能量（TOTEN、`energy(sigma->0)`）、每原子的 T*S（用来检查 SIGMA 是否合适）、E-fermi，以及 ISIF ≥ 3 时的外压（含 Pulay stress）。
- 最大受力：考虑 POSCAR/CONTCAR 的 Selective dynamics，把力投影到允许移动的方向上（对非正交晶胞同样准确），并与 EDIFFG 比较。
- 总磁矩和最大的位点磁矩（共线、非共线均支持）。
- 耗时、平均每个离子步的耗时、内存。
- VASP 6 的 WARNING/ADVICE 框，以及 OUTCAR 和标准输出（`vasp.out`、`slurm-*.out`、`*.o<jobid>` 等）中的已知报错，每条都附有中文修复建议（ZBRENT、BRMIX、EDDDAV、Sub-Space-Matrix、IBZKPT、Tetrahedron/BZINTS、PRICEL、ZPOTRF、超时、内存不足等）。

**依赖**：Python ≥ 3.9，仅标准库（装了 numpy 时会使用）。

**用法**：

```bash
python check_vasp.py                        # 当前目录
python check_vasp.py relax/ scf/ band/      # 多个作业
python check_vasp.py */ --brief             # 每个目录输出一行
python check_vasp.py neb/                   # 自动识别 00..NN 图像目录
python check_vasp.py -o vasp.out --json report.json
python check_vasp.py --strict && echo converged   # 没有收敛时退出码为 1
python check_vasp.py -v                     # 列出每个离子步的 SCF 步数
```

**参数**：`dirs`（作业目录，默认 `.`）、`-o/--stdout FILE`（指定要扫描的标准输出文件，可重复，默认自动识别）、`--brief`、`--json FILE`、`--strict`、`--top N`（显示的最大位点磁矩数，默认 5）、`--max-boxes N`（最多显示的 WARNING/ADVICE 框数，默认 8）、`--running-age SEC`（OUTCAR 在多少秒内被修改过就判为"运行中"，默认 900）、`-v/--verbose`、`--version`。

**退出码**：0 = 正常；1 = 发现错误（或在 `--strict` 下有作业未收敛）；2 = 没有找到 VASP 输出。

**局限**：报错识别靠字符串匹配。VASP 新版本改写报错文字后，可能只能给出通用提示。修复建议是经验性的，处理前请先看 OUTCAR 中完整的报错上下文（见[第18章](../tutorials/18-postprocess-troubleshooting.md)）。

---

## 6. conv_test.sh — 收敛测试驱动

**用途**：为每个测试值建一个子目录（如 `conv_encut/encut_400/`），复制 INCAR/POSCAR/KPOINTS/POTCAR（以及可选的作业脚本），然后修改被测参数。具体做法：先删除 INCAR 中原有的同名标签（包括 `A = 1; B = 2` 这种写在同一行的），去掉 Windows 的 CR 字符；k 点测试时写 Gamma 或 MP 网格的 KPOINTS，KSPACING 测试时写 `KSPACING`/`KGAMMA`。默认改成单点计算（`NSW = 0`，`IBRION = -1`，`LWAVE`/`LCHARG = .FALSE.`）。最后写出 `conv_info.txt`（供 `conv_analyze.py` 读取）和 `run_all.sh`，可选择逐个提交或在本机依次运行。

**依赖**：bash ≥ 4、awk、coreutils；`-s` 需要 `sbatch` 或 `qsub`。

**用法**：

```bash
bash conv_test.sh -t encut  -v "300 350 400 450 500 550" -j job_vasp_slurm.sh -s
bash conv_test.sh -t kpoints -v "4 6 8 10 12" --mp         # n -> n n n
bash conv_test.sh -t kpoints -v "6x6x1 9x9x1 12x12x1"      # 显式网格（二维）
bash conv_test.sh -t kpoints -v "6 9 12" --2d              # n -> n n 1
bash conv_test.sh -t kspacing -v "0.5 0.4 0.3 0.25 0.2 0.15"
bash conv_test.sh -t sigma -v "0.2 0.1 0.05 0.02 0.01" --ismear 1
bash conv_test.sh -t LREAL -v "Auto .FALSE."               # 任意 INCAR 标签
bash conv_test.sh -t encut -v "300 400 500" -r "mpirun -np 8 vasp_std"
bash conv_test.sh -t encut -v "400 500" -c "sbatch -p debug" -j job_vasp_slurm.sh
bash conv_test.sh -t encut -v "400 500" -n                 # 只打印，不执行
```

**参数**：

| 参数 | 说明 |
|---|---|
| `-t, --type TYPE` | 必填：`encut`、`kpoints`、`kspacing`、`sigma`，或任意 INCAR 标签名 |
| `-v, --values "LIST"` | 必填：空格分隔的测试值 |
| `-i/-p/-k/-P` | 模板 INCAR / POSCAR / KPOINTS / POTCAR（默认同名文件；缺少 POTCAR 时只警告） |
| `-j, --job FILE` | 作业脚本，复制到每个目录并命名为 `job.sh` |
| `-d, --outdir DIR` | 父目录（默认 `conv_<type>`） |
| `-s, --submit` | 逐个提交（自动识别 sbatch/qsub） |
| `-c, --submit-cmd CMD` | 指定提交命令，同时隐含 `-s` |
| `-r, --run CMD` | 在本机依次运行，例如 `"mpirun -np 8 vasp_std"` |
| `--gamma` / `--mp` | k 点测试的网格类型（默认 Gamma） |
| `--2d` | 标量 k 值 n 写成 n n 1（slab、二维材料） |
| `--ismear N` | SIGMA 测试时一并设置 ISMEAR |
| `--keep-ionic` / `--keep-io` | 保留模板中的 NSW/IBRION / LWAVE/LCHARG |
| `-f, --force` / `-n, --dry-run` / `-h, --help` | 覆盖已有目录 / 只打印 / 帮助 |

**输出**：`<outdir>/<type>_<value>/`、`<outdir>/conv_info.txt`、`<outdir>/run_all.sh`。全部算完后运行 `python conv_analyze.py <outdir>` 分析。

**局限**：标签替换基于文本处理，INCAR 中复杂的续行写法需要人工检查。`kspacing` 测试会删除 KPOINTS 文件。弹性常数等对 ENCUT 更敏感的性质，应以目标性质而不是总能量作为收敛判据（见[第13章](../tutorials/13-elastic.md)）。

---

## 7. conv_analyze.py — 收敛测试分析

**用途**：读取 `conv_test.sh` 生成的目录（或任意一组 VASP 计算目录），列表输出每个计算的以下内容：`energy(sigma->0)`（`--energy toten` 时用 TOTEN）、每原子能量、相对参考计算的偏差（meV/atom）、最大受力及其偏差、外压、每原子 T*S、不可约 k 点数和耗时。从某个值起 |dE| 一直低于 `--tol`，就把这个值报告为收敛值。SIGMA 测试还要求 |T*S| 低于 `--tstol`。

被测参数的来源依次为：`conv_info.txt`、目录名（`encut_400`、`kpoints_8x8x8`、`sigma_0.05` 等）、`--param` 加 OUTCAR。参考计算默认取最精确的那个（最大 ENCUT、最密 k 网格、最小 KSPACING、最小 SIGMA）。

**依赖**：`pip install numpy matplotlib`

**用法**：

```bash
python conv_analyze.py conv_encut
python conv_analyze.py conv_kpoints --tol 0.5 --ftol 0.01
python conv_analyze.py encut_* --param ENCUT -o encut_test
python conv_analyze.py conv_sigma --energy toten --tstol 0.5
python conv_analyze.py conv_kpoints --ref 8x8x8     # 与 --ref 8 或 --ref 512 等价
python conv_analyze.py conv_encut --show            # 弹出交互窗口
```

**参数**：`paths`（测试父目录或各计算目录）、`--param TAG`、`--energy e0|toten`（默认 `e0`）、`--tol`（meV/atom，默认 1）、`--ftol`（eV/Å，可选）、`--ref VALUE`、`--tstol`（meV/atom，默认 1）、`--no-constraints`（计算受力时忽略 Selective dynamics）、`-o/--output PREFIX`（默认 `conv_<param>`）、`--no-plot`、`--show`、`--dpi`（默认 200）。

**k 点测试的 `--ref`**：k 网格的横坐标是 N1×N2×N3。`--ref` 可以写网格（`8x8x8`）、乘积（`512`）或单个 N（`8`，对应 8×8×8，或 slab 的 8×8×1，即所有大于 1 的分量都等于 N）。单个数字先按 N 匹配，匹配不到再按乘积匹配。

**输出**：屏幕上的表格、`<prefix>.dat`（空白分隔的表格）和 `<prefix>.png`。

**局限**：判据是"相对最精确计算的偏差"，所以测试范围本身必须足够宽。脚本只读 OUTCAR/OSZICAR，不检查各次计算的其他参数是否一致。

---

## 8. plot_band.py — 能带图与带隙

**用途**：从 `vasprun.xml` 或 `EIGENVAL`（加 POSCAR/CONTCAR 提供倒格子）读取本征值并画能带图，报告带隙、VBM、CBM，分自旋通道，并判断直接/间接带隙。支持以下情况：

- VASP ≥ 6.3 的 KPOINTS_OPT（`vasprun.xml` 中的 `<eigenvalues_kpoints_opt>` 块）。
- 自旋极化（ISPIN = 2）和非共线计算。
- HSE 零权重 k 点能带（自动去掉带权重的 SCF 网格点）。
- 高对称点标签可以来自线模式 KPOINTS/KPOINTS_OPT、显式 k 点列表中的标签、vaspkit 的 `KPATH.in` 或 `--labels`。
- 用 PROCAR/PROCAR_OPT（LORBIT = 11/12）画 fat band。
- 导出 Origin/gnuplot 可用的列数据。

**依赖**：`pip install numpy matplotlib`

**用法**：

```bash
python plot_band.py                                   # ./vasprun.xml（标签取自 ./KPOINTS）
python plot_band.py -f band/vasprun.xml --efermi scf/OUTCAR -o Si_band.png
python plot_band.py -f EIGENVAL --kpoints KPOINTS --poscar POSCAR --efermi 5.76
python plot_band.py -f vasprun.xml --opt              # KPOINTS_OPT 能带（VASP >= 6.3）
python plot_band.py -f vasprun.xml --kpath KPATH.in   # HSE 零权重（标签来自 vaspkit）
python plot_band.py --procar PROCAR --proj "Mo:d" --proj "S:p" --ylim -3 3
python plot_band.py --labels "L G X U|K G" --zero vbm --dat band.dat
```

**参数**：

| 参数 | 说明 |
|---|---|
| `-f, --file` | `vasprun.xml` 或 `EIGENVAL`（默认先找 `./vasprun.xml`，再找 `./EIGENVAL`） |
| `--kpoints` / `--kpath` / `--labels` | 标签来源；`--kpoints none` 禁用；`--labels` 中用竖线分隔不连续的路径段 |
| `--poscar` | EIGENVAL 模式下提供晶格和元素名 |
| `--efermi` | 数值，或从文件（OUTCAR、vasprun.xml、DOSCAR）中读取；ICHARG = 11 的计算请传 SCF 的 OUTCAR |
| `--zero {efermi,vbm,none}` | 能量零点（默认 efermi） |
| `--opt` / `--no-opt` | 强制使用 / 忽略 KPOINTS_OPT 块 |
| `--keep-weighted` | 存在零权重点时也保留带权重的点 |
| `--gap-method {auto,occ,efermi}`、`--gap-tol`（默认 0.05 eV）、`--nocc N` | 判断占据态的方式 |
| `--procar`、`--proj "ATOMS:ORBITALS"`（可重复）、`--fat-scale`（默认 60） | fat band；ATOMS 为元素、序号范围（`1-4,7`）或 `all`，ORBITALS 为 `s,p,d,f,t2g,eg,px,dz2,...` 或 `all` |
| `--ylim`（默认 -6 6）、`--figsize`、`--dpi`、`--title`、`-o/--output`（默认 `band.png`） | 作图 |
| `--dat FILE`、`--no-plot`、`--show` | 导出数据 / 只输出带隙 / 交互窗口 |

**局限**：非自洽（ICHARG = 11）计算中的 E-fermi 不可靠，请用 `--efermi` 传入 SCF 的 OUTCAR。金属和半金属的"带隙"判断取决于 `--gap-method` 和 `--gap-tol`。有效质量、能带反折叠等更复杂的分析见[第17章](../tutorials/17-advanced-electronic.md)。

---

## 9. plot_dos.py — 态密度与 d 带中心

**用途**：读取 DOSCAR（LORBIT = 10 或 11/12），画总 DOS 和投影 DOS，支持 ISPIN = 1/2 和非共线计算（非共线时只用 PDOS 的总分量）。投影可以按元素、原子序号和轨道组合选择，自旋向下的 DOS 画成负值。还可以计算任意选择的 d 带中心（一阶矩）、RMS 带宽（二阶矩）和填充度，并可额外做高斯展宽、导出列数据。

**依赖**：`pip install numpy matplotlib`

**用法**：

```bash
python plot_dos.py                                   # ./DOSCAR，按元素画 PDOS
python plot_dos.py -f DOSCAR -p "Ni:d" -p "O:p" --xlim -8 4 -o NiO_pdos.png
python plot_dos.py -p "Ni:t2g" -p "Ni:eg" --smear 0.05
python plot_dos.py -p "1-4,9:eg" -p "Pt:dz2"
python plot_dos.py --by-element --by-orbital spd --dat dos.dat
python plot_dos.py --dband "Pt" --dband-window -10 2 --no-plot
```

**参数**：`-f/--file`（默认 `DOSCAR`）、`--poscar`、`-p/--pdos "ATOMS:ORBITALS"`（可重复；轨道可以是 `s,p,d,f,t2g,eg,all` 或 `px`、`dz2`、`x2-y2` 等 lm 名称）、`--by-element`、`--by-orbital ORBS`、`--efermi`、`--ncl`、`--smear`（eV，默认 0）、`--dband SEL`（可重复）、`--dband-window EMIN EMAX`（相对 E-fermi）、`--xlim`（默认 -8 4）、`--ylim`、`--no-tdos`、`--tdos-scale`、`--figsize`、`--dpi`、`--title`、`-o/--output`（默认 `dos.png`）、`--dat`、`--no-plot`、`--show`。

**输入/输出**：E-fermi 取自 DOSCAR 第 6 行，元素名取自同目录的 CONTCAR/POSCAR。是否为非共线计算由同目录的 OUTCAR/INCAR 判断，或用 `--ncl` 强制指定。输出图片和可选的 `.dat` 列数据。

**局限**：t2g/eg 的划分以 VASP 的笛卡尔坐标轴为准，八面体没有对齐坐标轴时没有物理意义。d 带中心对能量窗口和 NEDOS 都敏感，比较不同体系时要用相同的设置（见[第10章](../tutorials/10-adsorption-catalysis.md)）。

---

## 10. locpot_workfunction.py — 功函数

**用途**：对 LOCPOT（或任何 CHGCAR 格式的体数据文件）沿指定轴做平面平均，可选再做宏观平均（一次或两次滑动平均，用于界面能带对齐）。脚本根据原子位置定位真空区，在其中寻找平台，计算功函数 Φ = V_vac − E_F。有偶极修正（`LDIPOL = .TRUE.`、`IDIPOL = 3`）时，会分别报告两侧的真空能级和两个功函数。E-fermi 取自 OUTCAR 中最后一个 `E-fermi`。对 VASP ≥ 6.4.3 中 `LVACPOTAV = .TRUE.` 的计算，脚本会把 OUTCAR 里 VASP 自己给出的上下两侧真空能级一起列出来对照。真空区电势不平时会给出警告。

**依赖**：`pip install numpy matplotlib`

**用法**：

```bash
python locpot_workfunction.py                          # ./LOCPOT、./OUTCAR，沿 z 轴
python locpot_workfunction.py -f LOCPOT --outcar OUTCAR -o wf.png --dat planar.dat
python locpot_workfunction.py --efermi -1.234 --axis z --macro 2.27
python locpot_workfunction.py --macro 3.2 2.9          # 双重宏观平均（界面）
```

**参数**：`-f/--file`（默认 `LOCPOT`）、`--outcar`、`--efermi`、`--axis`（默认 z）、`--macro L [L2]`（滑动平均的周期，Å）、`--grad-tol`（平台判定的梯度阈值，默认 0.01）、`--min-len`（平台最小长度，Å，默认 1.5）、`--figsize`、`--dpi`、`--title`、`-o/--output`（默认 `workfunction.png`）、`--dat`、`--no-plot`、`--show`。

**局限**：LOCPOT 必须用 `LVHAR = .TRUE.` 写出（只包含静电势，即 Hartree 项加离子项），脚本只用第一个数据块。真空层太薄，或不对称 slab 没有加偶极修正时，结果不可靠（见[第9章](../tutorials/09-surface-workfunction.md)）。

---

## 11. chgdiff.py — 差分电荷密度与 Bader 辅助

**用途**：有四个子命令。

- `diff`：检查各个 CHGCAR 的 FFT 网格和晶胞是否一致，片段中的原子相对整体结构移动了位置时给出警告，然后计算 Δρ = ρ(AB) − ρ(A) − ρ(B) − …，写出 VESTA、vaspkit、bader 都能读的 CHGCAR 格式文件，并画出沿某一轴的平面平均 Δλ(z)（e/Å）和累积电荷转移 ΔQ(z)。
- `sum`：把多个 CHGCAR 格式文件相加，典型用法是 AECCAR0 + AECCAR2 → `CHGCAR_sum`，作为 `bader CHGCAR -ref CHGCAR_sum` 的参考密度。
- `planar`：对任意 CHGCAR 格式文件做平面平均。
- `check`：只比较几个文件的网格、晶胞和原子位置。

**依赖**：`pip install numpy matplotlib`（matplotlib 只在作图时需要）

**用法**：

```bash
python chgdiff.py diff AB/CHGCAR slab/CHGCAR mol/CHGCAR -o CHGCAR_diff
python chgdiff.py diff AB/CHGCAR slab/CHGCAR mol/CHGCAR --axis z --dat cdd_z.dat
python chgdiff.py sum AECCAR0 AECCAR2 -o CHGCAR_sum      # 然后：bader CHGCAR -ref CHGCAR_sum
python chgdiff.py planar CHGCAR_diff --axis z --plot planar.png
python chgdiff.py check AB/CHGCAR slab/CHGCAR mol/CHGCAR
python chgdiff.py diff AB/CHGCAR A/CHGCAR B/CHGCAR --block mag   # 磁化密度差
```

**参数**：

- 各子命令共用：`files`、`--block total|mag`（自旋极化时选总密度块或磁化密度块，默认 total）、`--tol-cell`（Å，默认 1e-4）、`--tol-pos`（默认 0.02）。
- `diff`/`sum`：`-o/--output`。
- `diff`/`planar` 的作图选项：`--axis`、`--dat`、`--plot`（默认分别为 `chgdiff_planar.png`、`planar.png`）、`--no-plot`、`--figsize`、`--dpi`、`--show`。

**局限**：输出文件不包含 augmentation occupancies（VESTA 和 bader 不需要）。三个计算必须使用完全相同的晶胞、FFT 网格（相同的 ENCUT、PREC，或显式给出 NGXF 等）和原子位置，详见[第16章](../tutorials/16-charge-bonding.md)。

---

## 12. neb_tools.py — NEB 全流程辅助

**用途**：VASP/VTST NEB 计算的辅助工具，有以下子命令。

- `interp`：用线性或 IDPP（ASE）插值生成 `00 .. N+1` 目录（与 VTST 的 `nebmake.pl` 布局相同）。保留初态的 Selective dynamics，默认使用最小镜像约定。
- `check IS FS`：检查初末态的元素、原子顺序和晶胞是否一致，给出每个原子的位移，列出位移可疑地大的原子（往往是原子顺序错了），并用匈牙利算法求最优重排（`--write-reordered`）。
- `check DIR`：给出相邻图像的间距（类似 VTST 的 `dist.pl`）、每个图像中最短的原子间距和路径总长。
- `collect`：收集每个图像的能量和 VTST 切向力，从几何结构计算反应坐标。有 VTST 力时用基于力的三次 Hermite 样条，没有时用自然三次样条。输出正反向能垒和反应能，写出 `neb.dat`、`spline.dat` 和 PNG 图。
- `xyz`：导出整条带的扩展 XYZ 动画，坐标经过展开，原子不会跨晶胞边界跳动。

**依赖**：`pip install numpy matplotlib scipy`；`--method idpp` 还需要 `pip install ase`。

**用法**：

```bash
python neb_tools.py interp POSCAR_IS POSCAR_FS -n 5                  # 线性插值
python neb_tools.py interp POSCAR_IS POSCAR_FS -n 5 --method idpp \
       --is-outcar is/OUTCAR --fs-outcar fs/OUTCAR
python neb_tools.py check POSCAR_IS POSCAR_FS --write-reordered POSCAR_FS_new
python neb_tools.py check .                  # 00..NN 的图像间距
python neb_tools.py collect . -o neb         # neb.dat, spline.dat, neb.png
python neb_tools.py xyz . -o movie.xyz
```

**参数**：

| 子命令 | 参数 |
|---|---|
| `interp` | `initial final`、`-n/--nimages N`（必填，= IMAGES）、`--method {linear,idpp}`、`-o/--outdir`、`--no-mic`、`--keep-unwrapped`、`--allow-cell-change`（SS-NEB）、`--is-outcar`、`--fs-outcar`、`--idpp-steps`（默认 200）、`--idpp-fmax`（默认 0.05）、`--dmin`（默认 0.7 Å）、`-f/--force` |
| `check` | `paths`（IS FS 两个文件，或一个 DIR）、`--move`（默认 0.3 Å）、`--max-disp`（默认 2 Å）、`--top`（默认 15）、`--write-reordered FILE`、`--poscar`、`--dmin` |
| `collect` | `dir`（默认 `.`）、`--energy {e0,toten}`、`--poscar`、`-o/--output`（默认 `neb`）、`--points`（默认 300）、`--title`、`--no-plot`、`--show`、`--dpi` |
| `xyz` | `dir`、`-o/--output`（默认 `neb_movie.xyz`）、`--poscar` |

**局限**：切向力取自 VTST 输出的 `NEB: projections on to tangent (spring, REAL)` 行，最大 NEB 力取自 `FORCES: max atom, RMS` 行。不用 VTST 时这些行不存在，只能对能量做普通样条。端点 `00` 和 `NN` 中要放入初末态优化得到的 OUTCAR（和 `nebresults.pl` 的要求相同）。MPI 进程数必须能被 IMAGES 整除。建议先跑普通 NEB，再开 LCLIMB（见[第11章](../tutorials/11-neb-dimer.md)）。

---

## 13. phonopy_workflow.sh — phonopy 声子工作流

**用途**：把 VASP + phonopy 的有限位移（fd）或 DFPT 声子计算拆成若干子命令。

| 子命令 | 作用 |
|---|---|
| `prepare` | 运行 `phonopy-init -d`（phonopy ≥ 4）或 `phonopy -d`（2.x/3.x），生成 SPOSCAR 和 POSCAR-NNN。fd 方式下为每个位移建一个目录（`disp-001` …，内含 POSCAR、INCAR、KPOINTS、POTCAR、job.sh）；dfpt 方式下建 `dfpt/`，以 SPOSCAR 为 POSCAR（INCAR 中必须设 IBRION = 7 或 8）。给了 `--born-incar` 时再建 `born/`（原胞，LEPSILON 或 LCALCEPS），用于非解析项修正。设置保存到 `phonopy_workflow.cfg`，后续步骤会读取 |
| `submit` | 自动识别 sbatch/qsub，提交所有还没有完整 vasprun.xml 的目录 |
| `check` | 检查 vasprun.xml 是否以 `</modeling>` 结尾；fd 方式下检查是否恰好 1 个离子步、SCF 是否低于 NELM；DFPT 检查是否有 hessian；`born/` 检查是否有 Born 有效电荷 |
| `collect` | 生成 FORCE_SETS（fd）或 FORCE_CONSTANTS（dfpt）；有 `born/` 时用 `phonopy-vasp-born` 生成 BORN |
| `conf` | 写出 `band.conf`、`mesh.conf`（DOS/PDOS）和 `tprop.conf`（热力学性质） |
| `post` | 算能带、DOS/PDOS 和热力学性质并作图，用 `phonopy-bandplot --gnuplot` 导出 `band.dat`，再检测虚频 |
| `imag` | 在 `band.yaml`/`mesh.yaml` 中查找虚频（负频率） |
| `all` | `prepare` + `submit` |

**依赖**：bash ≥ 4、awk、phonopy（2.x/3.x 或 ≥ 4.0）、seekpath（Python 模块，`BAND = AUTO` 要用）；`submit` 需要 sbatch 或 qsub。

```bash
pip install phonopy seekpath
```

**用法**：

```bash
bash phonopy_workflow.sh prepare -c POSCAR --dim "2 2 2" -j job_vasp_slurm.sh
bash phonopy_workflow.sh submit
bash phonopy_workflow.sh check
bash phonopy_workflow.sh collect
bash phonopy_workflow.sh post --mesh "24 24 24" --pdos "1, 2"
bash phonopy_workflow.sh imag --imag-tol 0.2
# DFPT + Born 有效电荷（面心原胞的 PRIMITIVE_AXES）
bash phonopy_workflow.sh prepare -m dfpt -c POSCAR-unitcell -i INCAR_dfpt -k KPOINTS_dfpt \
     --born-incar INCAR_born --born-kpoints KPOINTS_born --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"
# 手动指定能带路径
bash phonopy_workflow.sh post --band "0 0 0 1/2 0 1/2 1/2 1/2 1/2" --labels '$\Gamma$ X L'
```

**参数**：

| 适用子命令 | 参数 |
|---|---|
| prepare / all | `-m/--method {fd,dfpt}`、`-c/--cell FILE`（默认 POSCAR）、`-d/--dim "N N N"`（3 个或 9 个整数，默认 `"2 2 2"`）、`--pa AXES`（默认 auto；也可以是 P、F、I…或 9 个数）、`-i/--incar`、`-k/--kpoints`、`-P/--potcar`、`-j/--job FILE`（默认按识别到的调度系统，选本目录下的 `job_vasp_slurm.sh` 或 `job_vasp_pbs.sh`）、`--amplitude A`（Å，phonopy 默认 0.01）、`--pm`、`--prefix NAME`（默认 `disp`）、`--born-incar FILE`、`--born-kpoints FILE`、`-f/--force` |
| submit | `-s/--submit-cmd CMD`，例如 `"sbatch -p debug"` |
| conf / post / imag | `--band "PATH"`（默认 `BAND = AUTO`）、`--labels "..."`、`--band-points N`（默认 101）、`--mesh "N N N"`（默认 `"24 24 24"`）、`--pdos "1 2, 3"`（原胞中的原子序号分组）、`--tmax T`（默认 1000 K）、`--tstep T`（默认 10 K）、`--imag-tol F`（THz，默认 0.1） |
| 其他 | `--phonopy-old`（强制使用 2.x/3.x 语法）、`-n/--dry-run`、`-h/--help` |

**版本识别**：脚本通过 `phonopy-init` 命令是否存在来判断是否为 phonopy ≥ 4。≥ 4 时把 `-d`、`-f`、`--fc` 交给 `phonopy-init`，这一拆分按 phonopy 4.x 的迁移说明核对过；否则使用 2.x/3.x 的 `phonopy -d/-f/--fc` 写法。

**`prepare` 的输入检查**：fd 方式下，若 INCAR 中 IBRION ≠ −1 或 NSW ≠ 0 会警告。DFPT 方式下，设置了 NCORE ≠ 1、NPAR 或 IVDW 会警告（DFPT 不支持 vdW 修正，应改用有限位移）。LREAL 不是 `.FALSE.` 时会警告（实空间投影会给力带来噪声）。INCAR 中有 MAGMOM 时会提醒它对应的是原胞，需要改成超胞的。

**输出**：`phonopy_disp.yaml`、`SPOSCAR`、`disp-NNN/` 或 `dfpt/`、`born/`、`FORCE_SETS`/`FORCE_CONSTANTS`、`BORN`、`*.conf`、`band.yaml`、`mesh.yaml`、`total_dos.dat`、`projected_dos.dat`、`thermal_properties.yaml`、PDF 图和 `band.dat`。

**退出码**：0 = 正常；1 = 出错、`check` 发现未完成的计算，或 `imag` 发现虚频。

**局限**：超胞的 KPOINTS 要按超胞大小缩小（保持与原胞相同的 k 点密度），脚本不会自动换算。`BAND = AUTO` 的路径对应 seekpath 的标准原胞。`phonopy-vasp-born` 只接受数值形式的 `--pa` 矩阵：`--pa` 不是数值时，born/POSCAR 会被当作原胞，如果它其实是带心或惯用胞，需要用 `--pa "9 个数"` 重新运行 `collect`。原理和收敛建议见[第12章](../tutorials/12-phonon.md)。

---

## 14. elastic_post.py — 弹性常数后处理

**用途**：解析 OUTCAR 中 `IBRION = 6`、`ISIF >= 3` 计算的 `ELASTIC MODULI ... (kBar)` 块（TOTAL、clamped-ion/symmetrized、离子弛豫贡献），完成以下处理：

- 把单位从 kBar 换成 GPa（除以 10），并把 VASP 的 XX YY ZZ XY YZ ZX 顺序重排为 Voigt 顺序 11 22 33 23 13 12。
- Born 力学稳定性：总是做本征值判据，另外按晶系做 Mouhat & Coudert 的显式判据（[PRB 90, 224104 (2014)](https://doi.org/10.1103/PhysRevB.90.224104)）。
- Voigt/Reuss/Hill 体模量和剪切模量、杨氏模量、泊松比、Pugh 比 B/G、通用各向异性指数 A^U、Zener 比和 Cauchy 压（立方晶系）、Chen 经验维氏硬度。
- 给出结构文件时，计算密度、纵波/横波/平均声速和德拜温度（Anderson 方法）。
- 二维材料（`--2d`）：面内刚度（N/m，即 C × c）、二维杨氏模量和泊松比，以及二维稳定性条件。

**依赖**：`pip install numpy`；可选 `pip install ase spglib`（原子质量、晶系识别）。

**用法**：

```bash
python elastic_post.py                                # ./OUTCAR（+ ./CONTCAR 或 ./POSCAR）
python elastic_post.py -f OUTCAR --poscar POSCAR --system cubic
python elastic_post.py --block symmetrized            # clamped-ion（刚性离子）模量
python elastic_post.py --matrix Cij_GPa.txt --system hexagonal --density 1.74
python elastic_post.py --2d --poscar CONTCAR          # 单层，单位 N/m
python elastic_post.py --json elastic.json --dat Cij_voigt.dat
```

**参数**：`-f/--file`（默认 `OUTCAR`）、`--block total|symmetrized|clamped|ionic`（默认 total，即包含离子弛豫）、`--matrix FILE`（从文本文件读 6×6 矩阵）、`--matrix-units GPa|kBar`、`--matrix-vasp-order`、`--poscar`、`--system auto|cubic|hexagonal|tetragonal|tetragonal2|trigonal|trigonal2|orthorhombic|monoclinic|triclinic`（tetragonal2 = 点群 4、−4、4/m；trigonal2 = 点群 3、−3）、`--density`（g/cm³）、`--2d`、`--c-length`（Å）、`--json`、`--dat`。

**局限**：晶系的确定顺序是：`--system`、spglib（需要结构文件）、最后才根据矩阵的形式猜测；最后一种不可靠。弹性常数对 ENCUT 和 k 点很敏感，应先做收敛测试（见[第13章](../tutorials/13-elastic.md)）。

---

## 15. optics_post.py — 光学与介电后处理

**用途**：有三个子命令。

- `optics`：读取 LOPTICS（以及 GW/BSE）计算的频率相关介电函数。数据来自 `vasprun.xml` 的 `<dielectricfunction>` 块：VASP 6 中标为 `density-density`/`current-current`，VASP 5 中没有标签的块按 density、current 的顺序读取。`vasprun.xml` 不完整时，可以改从 OUTCAR 读取。由 ε = ε1 + iε2（默认取 xx、yy、zz 的各向同性平均）计算折射率 n、消光系数 k、吸收系数 α（cm⁻¹）、正入射反射率 R、能量损失函数 −Im(1/ε)，以及给出 `--slab-length` 时的薄层吸光度 A。还可以画 Tauc 图 (αE)^(1/r)，并用 `--fit` 把线性段外推到能量轴。
- `static`：读取 LEPSILON/LCALCEPS 计算的静态介电张量（OUTCAR 中所有 `... DIELECTRIC TENSOR ...` 块，包括 IBRION = 5–8 的离子贡献）和 Born 有效电荷，输出 ε0 = ε∞ + ε_ion、平均值和主值、每个原子的 Z*，以及电荷中性求和规则的检查结果。也可以读 `vasprun.xml`。
- `list`：列出文件中找到的介电函数块和静态张量。

**依赖**：`pip install numpy matplotlib`（matplotlib 只在作图时需要）

**用法**：

```bash
python optics_post.py optics                              # ./vasprun.xml，density-density
python optics_post.py optics run/vasprun.xml --kind current
python optics_post.py optics --component zz --emax 8
python optics_post.py optics --tauc direct --fit 3.0 3.6
python optics_post.py optics --component inplane --slab-length 20   # 二维：吸光度
python optics_post.py optics OUTCAR                       # vasprun.xml 不完整时的备用方案
python optics_post.py static                              # ./OUTCAR（LEPSILON / LCALCEPS）
python optics_post.py static ionic/OUTCAR --json eps.json
python optics_post.py list vasprun.xml
```

**参数**：

- `optics`：`file`（默认 `vasprun.xml`）、`--kind density|current`（默认 density）、`--block SUBSTR`、`--component avg|inplane|xx|yy|zz|xy|yz|zx`、`--emin`（默认 0 eV）、`--emax`（默认 10 eV）、`--tauc direct|indirect`、`--fit EMIN EMAX`、`--slab-length LZ`、`--onset`（吸收边阈值，默认 1e4 cm⁻¹）、`-o/--prefix`（默认 `optics`）、`--no-plot`、`--show`、`--dpi`。
- `static`：`file`（默认 `OUTCAR`）、`--asr-tol`（默认 0.1）、`--json`。
- `list`：`file`。

**输出**（`optics`）：`<prefix>.dat`（E、波长、ε1、ε2、n、k、α、R、L[、A]）、`<prefix>_tensor.dat`（全部 6 个分量）、`<prefix>.png`；加 `--tauc` 时另有 `<prefix>_tauc.dat/png`。

**局限**：α、n、k、R 都是体相量，slab 或二维材料的这些值依赖真空层厚度，应改用 `--slab-length` 给出的吸光度。独立粒子近似的谱既不含激子效应，也不含声子辅助的间接跃迁。Tauc 外推得到的带隙依赖所选的拟合区间（见[第14章](../tutorials/14-optics-dielectric.md)）。

---

## 16. md_analysis.py — 分子动力学与 MLFF 分析

**用途**：有五个子命令。

- `thermo`：温度、能量（OSZICAR）和晶胞体积（XDATCAR）随时间的变化，平衡后的平均值，以及守恒能量 E 的漂移。
- `rdf`：元素对的径向分布函数 g(r) 和累积配位数，支持定容和变胞 XDATCAR。
- `msd`：按元素计算均方位移（展开坐标、去除质心漂移、通过 FFT 使用所有时间原点），由 Einstein 关系 D = slope / (2d) 求扩散系数。
- `mlff`：解析 ML_LOGFILE（VASP ≥ 6.3 的机器学习力场），画 ERR（训练 RMSE）和 BEEF（贝叶斯误差及阈值）曲线，并统计 STATUS 各状态出现的次数。
- `all`：对目录中输入文件齐全的分析全部运行一遍。

时间步长依次取自 `--potim`、OUTCAR、INCAR（都没有时默认 1 fs）。XDATCAR 每 NBLOCK 步写一帧，帧间隔由 `Direct configuration=` 后面的编号和 NBLOCK 推出，也可以用 `--dt` 覆盖。

**依赖**：`pip install numpy matplotlib`；可选 `pip install ase`（用原子质量做质心修正）。

**用法**：

```bash
python md_analysis.py thermo                          # ./OSZICAR（+ ./XDATCAR）
python md_analysis.py rdf --pair Si-Si --skip-ps 2 --rmax 6
python md_analysis.py rdf --pair O-H --pair O-O -d water_npt/
python md_analysis.py msd --skip-ps 2 --fit-range 1 5 --species Li
python md_analysis.py msd --species Li --dim 2 --dat msd.dat   # 层状材料的二维扩散
python md_analysis.py mlff -f ML_LOGFILE
python md_analysis.py all -d run1/ --skip-ps 1
```

**参数**：

- 所有子命令共用：`-d/--dir`（默认 `.`）、`--xdatcar`、`--oszicar`、`--potim`（fs）、`--dt`、`--skip N`（丢弃的帧数/步数）与 `--skip-ps T`（丢弃的平衡时间，ps）二选一、`--stride`、`--prefix`、`--figsize`、`--dpi`、`--no-plot`、`--show`。
- `rdf`：`--pair`（可重复，例如 `Si-Si`、`O-H` 或 `all`）、`--rmax`（默认取晶胞宽度的一半，且不超过 10 Å）、`--dr`（默认 0.02 Å）、`--dat`。
- `msd`：`--species`、`--fit-range T0 T1`（ps）、`--dim 1|2|3`、`--keep-com`、`--dat`。
- `mlff`：`-f/--file`（默认 `<dir>/ML_LOGFILE`）。
- `all`：可以同时使用 rdf、msd、mlff 的参数。

**输出**：`<prefix>thermo.png`、`<prefix>rdf.png`、`<prefix>msd.png`、`<prefix>mlff.png`，以及可选的 `.dat` 文件。

**局限**：扩散系数需要足够长的轨迹，MSD 在拟合区间内要呈线性。短轨迹得到的 D 误差很大，应报告多段或多条轨迹的统计误差。NVT 恒温器会影响动力学性质，详见[第15章](../tutorials/15-md-mlff.md)。

---

## 常见问题

| 现象 | 原因 | 解决 |
|---|---|---|
| `/usr/bin/env: 'python3\r': No such file or directory` | 脚本在 Windows 下保存成了 CRLF 换行 | `dos2unix *.py *.sh`，或用 `git config core.autocrlf input` 重新检出 |
| `ModuleNotFoundError: No module named 'matplotlib'` | 缺少依赖 | `pip install numpy matplotlib scipy`，集群上可以用 `pip install --user` 或 conda 环境 |
| 在集群登录节点作图时报 display 相关错误 | 没有图形界面 | 不要加 `--show`；脚本默认把图存成文件 |
| `qsub job_vasp_pbs.sh ncl` 没有生效 | qsub 不把命令行参数传给脚本 | 改用 `qsub -v VARIANT=ncl job_vasp_pbs.sh` |
| `phonopy_workflow.sh post` 报 `BAND = AUTO` 失败 | 没有安装 seekpath | `pip install seekpath`，或用 `--band` 手动给出路径 |

## 贡献与反馈

脚本有 bug，或者在某个 VASP 版本上解析出错，欢迎到 [GitHub Issues](https://github.com/zenggs05/DFT-Tutorials/issues) 反馈。请附上 VASP 版本、出错命令和相关输出片段（不要上传 POTCAR）。所有脚本采用 MIT 协议。

[返回目录](../README.md)
