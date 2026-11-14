# 第18章 常见错误排查、性能调优与版本迁移（Troubleshooting, Performance Tuning & Version Migration）

本章是 CP2K 篇的“排障手册”：按“输入解析 → 数据文件 → 结构与晶胞 → SCF → 运动层 → 运行环境”的顺序整理常见报错原文（unknown keyword/subsection、basis set/potential not found、`GEOMETRY wrong or EMAX_SPLINE too small`、cell volume too small、`Use the LSD option for an odd number of electrons`、`SCF run NOT converged`、HFX `CUTOFF_RADIUS` 警告、MPI 错配），并用一个故意写错 6 处的 fcc Cu 输入演示逐步排错。随后讲性能调优（MPI×OpenMP 配比、`OMP_NUM_THREADS`、`PREFERRED_DIAG_LIBRARY` ELPA/ScaLAPACK、GPU 构建、HFX `MAX_MEMORY` 与 `NPROC_REP`）、重启策略（`-1.restart`、`-RESTART.wfn`/`.kp`、`RESTART_HISTORY`、`&EXT_RESTART` 各开关），以及 2025.x → 2026.2 的迁移清单（删除 `CELL_OPT TYPE`、`&PDOS` 改为 `&DOS/&PDOS`、REFTRAJ `EVAL`、只能 CMake 构建、删除 Cython 绑定）。最后给出一份结果可信度自检表。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 拿到一个失败的 CP2K 任务时，按固定流程定位问题所在的阶段，而不是盲目改参数；
- 认出最常见的十几种报错原文，知道每种报错的原因与修法；
- 识别“程序跑完了但结果不可信”的情形（电子数不对、展宽空轨道不足、MPI 错配导致串行运行等）；
- 为一个具体体系选出合理的 MPI×OpenMP 布局，读懂 timing 报告，处理杂化泛函的内存问题；
- 正确续算中断的优化/MD/NEB 任务，知道 `.wfn` 与 `.kp`、`-1.restart` 与 `&EXT_RESTART` 各自的作用和限制；
- 把 2025.x/2026.1 的输入、作业脚本和编译流程迁移到 CP2K 2026.2；
- 在论文或报告发表前，用一张自检表确认计算参数已经收敛、结果与文献或其他程序可比。

| 你遇到的情况 | 先看本章哪一节 |
|---|---|
| 程序一启动就中止，输出很短 | §3 排查流程、§4.1 输入解析错误 |
| 报找不到基组/赝势/数据文件 | §4.2 |
| 报结构或晶胞错误 | §4.3 |
| SCF 不收敛、奇数电子、HFX 警告 | §4.4 及[第5章](05-scf-strategy.md)、[第9章](09-hybrid-admm.md) |
| 输出重复、找不到动态库、`Illegal instruction`、被系统杀掉 | §4.5、§6 |
| 算得太慢、内存不够 | §7 |
| 任务被时限打断，需要续算 | §8 |
| 旧输入在 2026.2 上报错 | §9 和迁移清单 |
| 结果与文献差得多 | §10 |

本章示例目录：[../examples/18-troubleshooting/README.md](../examples/18-troubleshooting/README.md)。其中有故意出错的输入 `broken_scf.inp` 与修正版 `broken_scf_fixed.inp`、一对新旧语法的 Si 晶胞优化输入、晶胞优化续算输入、输入迁移扫描脚本 `check_migration.sh`、输出快速诊断脚本 `quick_diagnose.sh`、MPI×OpenMP 基准脚本 `job_mpi_omp_hybrid.sh`，以及迁移清单 [migrate_2025_to_2026_2.md](../examples/18-troubleshooting/migrate_2025_to_2026_2.md)。

> [!NOTE]
> 本章所有“示意输出”只演示格式和要找的关键行，数值不是真实计算结果。报错原文以 CP2K 2026.2 源码和官方文档为准（见各处说明）；不同版本的措辞可能略有差别，排查时用关键片段 `grep` 即可。

---

## 2. 原理速览

### 2.1 CP2K 的一次运行分几个阶段

理解报错出现在哪个阶段，是排错的第一步：

| 阶段 | 做什么 | 这一阶段的典型报错 |
|---|---|---|
| ① 预处理与解析 | 展开 `@INCLUDE`/`@SET`/`@IF`，按输入结构（`cp2k.psmp --xml` 能导出）逐段核对关键字 | unknown keyword/subsection、`not available anymore`、类型不符 |
| ② 读取数据文件 | 打开 `BASIS_SET_FILE_NAME`、`POTENTIAL_FILE_NAME`、`dftd3.dat`、`t_c_g.dat` 等 | 文件打不开、基组/赝势名找不到 |
| ③ 体系初始化 | 建立晶胞、邻居表、多重网格、分布矩阵，统计电子数 | 晶胞体积过小、`GEOMETRY wrong or EMAX_SPLINE too small`、`Use the LSD option for an odd number of electrons` |
| ④ SCF | 构造 KS 矩阵、OT 或对角化 + 混合 | `OT not possible with kpoint calculations`、`SCF run NOT converged` |
| ⑤ 运动层 | GEO_OPT/CELL_OPT/MD/BAND/振动分析的外层循环 | 达到 `MAX_ITER`、副本日志里的错误 |
| ⑥ 收尾 | 写重启文件、timing 报告、`PROGRAM ENDED AT` | 磁盘满、被作业系统杀掉时没有这一段 |

CP2K 内部有两类消息：`cp_warn` 打印 `*** WARNING in 文件名:行号 :: ... ***` 后继续运行；`cp_abort` 打印一个带 `[ABORT]` 的方框后终止全部 MPI 进程。方框右下角是源码位置（文件名:行号），向开发者或论坛求助时要一并提供。

（示意输出，方框格式与 2026.2 一致，源码行号仅供参考）

```text
 *******************************************************************************
 *   ___                                                                       *
 *  /   \                                                                      *
 * [ABORT]                                                                     *
 *  \___/    SCF run NOT converged. To continue the calculation regardless,    *
 *    |                please set the keyword IGNORE_CONVERGENCE_FAILURE.      *
 *  O/|                                                                        *
 * /| |                                                                        *
 * / \                                                         qs_scf.F:NNN    *
 *******************************************************************************
```

> [!IMPORTANT]
> 用 `-o x.out` 运行时，**正常输出写进 x.out，但部分错误信息写到标准错误（stderr）**。在 Slurm 下 stderr 进 `slurm-<jobid>.out` 或 `#SBATCH --error` 指定的文件。只看 `x.out` 末尾找不到原因时，一定要看作业日志。NEB、振动分析等多副本任务的报错往往只出现在 `<PROJECT>-BAND*.out`、`<PROJECT>-r-*.out` 之类的副本日志中（官方 Troubleshooting 页对此有专门说明）。

### 2.2 并行模型：MPI × OpenMP

`cp2k.psmp` 是 MPI + OpenMP 混合并行程序。设每个节点 $n_\mathrm{rank}$ 个 MPI 进程、每进程 $n_\mathrm{thr}$ 个线程，则应满足

$$
n_\mathrm{rank}\times n_\mathrm{thr} = N_\mathrm{core,node}
$$

时间可以粗略地拆成串行部分、可并行部分和通信：

$$
T(N) \approx T_\mathrm{s} + \frac{T_\mathrm{p}}{N} + T_\mathrm{comm}(N),\qquad
E(N) = \frac{T(N_0)\,N_0}{T(N)\,N}
$$

$E(N)$ 是相对于参考核数 $N_0$ 的并行效率。MPI 进程多，稀疏矩阵乘法（DBCSR）和 FFT 的分布更细，但通信量和**每个进程复制的数据**也增多；线程多，内存占用降低，但并非所有子程序都有良好的 OpenMP 实现。所以最优布局取决于体系和方法，只能测出来（§7）。

内存方面，杂化泛函的 ERI 存储按 **MPI 进程**计：

$$
M_\mathrm{node} \approx n_\mathrm{rank}\,\left(M_\mathrm{base} + M_\mathrm{HFX}\right),\qquad M_\mathrm{HFX}\le \texttt{MAX\_MEMORY}
$$

同一进程的线程共享这块内存，因此“少进程多线程”是解决 HFX 内存不足的首选（[第9章](09-hybrid-admm.md) §10）。

### 2.3 重启文件的两条线

CP2K 有两类相互独立的重启文件（官方 Restarting 页）：

- **运动重启**（motion restart）：`<PROJECT>-1.restart`，格式化文本，本身就是一个完整的 CP2K 输入，包含坐标、速度、晶胞、恒温器/恒压器状态、NEB 副本、步数计数器等；由 `&MOTION/&PRINT/&RESTART` 和 `&RESTART_HISTORY` 写出；只有带 `&MOTION` 的任务（GEO_OPT、CELL_OPT、MD、BAND 等）才写。
- **波函数重启**：Γ 点为 `<PROJECT>-RESTART.wfn`，显式 `&KPOINTS` 为 `<PROJECT>-RESTART.kp`；二进制，不同版本、不同机器之间不保证兼容；由 `&SCF/&PRINT/&RESTART` 写出。

续算一个长任务通常两者都要用：运动重启保证轨迹或优化连续，波函数重启让第一步 SCF 快速收敛（§8）。

---

## 3. 排查的基本流程

### 3.1 五步法

```bash
# 0) Is the job still running, or was it killed? (scheduler view)
squeue -u $USER ; sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed

# 1) Did CP2K finish normally?
grep -c "PROGRAM ENDED AT" x.out          # 0 = abnormal termination
tail -40 x.out                            # the [ABORT] box, if any
tail -40 slurm-<jobid>.out                # stderr: MPI/library/OOM messages

# 2) Syntax check with the SAME binary you run in production
cp2k.psmp -c x.inp                        # parse only, then stop
cp2k.psmp -e x.inp > x.expanded.inp       # echo the fully expanded input (@SET/@INCLUDE/@IF resolved)

# 3) Environment: version, features, data dir, MPI
cp2k.psmp -v | head -20                   # version string + cp2kflags
grep -E "CP2K\| Data directory path|DBCSR\| MPI: Number of processes|DBCSR\| OMP: Current number of threads" x.out

# 4) Warnings and physical sanity
grep "The number of warnings for this run is" x.out
grep "\*\*\* WARNING" x.out | sort | uniq -c | sort -rn | head
grep -E "Number of electrons|Number of occupied orbitals" x.out | head -4
```

| 步骤 | 检查什么 | 判据 |
|---|---|---|
| 0 | 作业状态与退出码 | `OUT_OF_MEMORY`、`TIMEOUT`、`NODE_FAIL` 都不是 CP2K 输入问题 |
| 1 | `PROGRAM ENDED AT` | 没有这一行就说明不是正常结束；先读 `[ABORT]` 方框，再读 stderr |
| 2 | `cp2k.psmp -c` | 能通过才说明输入语法对**这个版本**有效；`-e` 展开后检查变量是否替换正确 |
| 3 | `cp2k.psmp -v` 与输出头 | `cp2kflags` 中有没有所需功能（如 `libint`、`elpa`、`spglib`、`libdftd4`）；数据目录是否正确；MPI 进程数只出现一次且等于你申请的数目 |
| 4 | 警告和电子数 | 警告逐条读；电子数与化学式一致 |

示例目录中的 `quick_diagnose.sh` 把第 1、3、4 步常用的 `grep` 打包在一起，一次输出“是否正常结束、ABORT 内容、MPI/OMP、SCF 收敛次数、警告汇总、电子数、内存峰值、最耗时的子程序”。它只读文本，不调用 CP2K：

```bash
bash ../examples/18-troubleshooting/quick_diagnose.sh x.out
```

更完整的输出解析（逐步能量、力、Mulliken 电荷、优化每步的收敛指标，支持 `--json`）用本篇的 [`cp2k_outparse.py`](../scripts/cp2k_outparse.py)：`python ../scripts/cp2k_outparse.py -h` 查看用法。

### 3.2 一次只改一处

排错时最常见的反模式是“一口气改五个参数，跑通了也不知道是哪个起作用”。建议：

1. 把失败的输入复制一份作为基线，不在原文件上改；
2. 每次只改一处，跑最小的代表性任务（单点、`MAX_SCF` 适当小、`RUN_TYPE ENERGY`）；
3. 把每次改动和结果记在一张表里（改了什么 → 报错/SCF 步数/能量）；
4. 跑通后回头确认：被你改过的参数（如 `EPS_SCF`、展宽温度、`ADDED_MOS`）对目标性质是否已收敛。

§5 的练习就按这个流程逐个修掉 6 处错误。

---

## 4. 报错目录（按阶段）

### 4.1 输入解析阶段

CP2K 2026.2 解析器（`input_parsing.F`）产生的几类报错：

| 报错原文 | 原因 | 修正 |
|---|---|---|
| `found an unknown keyword XXX in section YYY` | 关键字拼错、放错了段，或该版本没有这个关键字 | 报错前输出中通常有 `Possible matches for unknown keyword` 列表，给出拼写相近的候选及得分；对照**同版本**输入参考修改 |
| `unknown subsection XXX of section YYY` | 子段拼错或位置不对，例如 2026.2 下写旧的 `&DFT/&PRINT/&PDOS` | 同上，输出中有 `Possible matches for unknown subsection` |
| `The specified keyword 'XXX' is not available anymore: <说明>` | 该关键字在当前版本**已被删除**，解析器保留了它的名字以便给出迁移说明 | 按说明删除或替换（2026.2 中的 `CELL_OPT/TYPE` 即此类，见 §9） |
| `The specified keyword 'XXX' is deprecated and may be removed in a future version: <说明>`（警告） | 关键字仍可用但已标为过时，例如 2026.2 中 `&XTB`/`&DFTB` 的 `DO_EWALD` | 计划迁移；需要兼容旧版本时暂时保留 |
| `A (integer\|floating point\|string) type object was expected, found ...` 并给出 File/Line/Column | 关键字后的值类型不对，或坐标文件格式不对（例如 `COORD_FILE_FORMAT XYZ` 却给了 CIF） | 按报错的行列号修改；XYZ 第一行必须是原子数，第二行是注释 |

> [!TIP]
> 2026.2 中被删除的 `CELL_OPT/TYPE` 给出的完整报错是：
>
> `The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE has been removed because cell optimizations now always use DIRECT_CELL_OPT.`
>
> 而 2026.2 中被删除的 REFTRAJ 关键字 `EVAL_ENERGY_FORCES`/`EVAL_FORCES` 没有保留迁移说明，报的是普通的 `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ`。

用预处理器时还要注意：`@IF`/`@ENDIF` 必须成对；`${VAR}` 未定义会直接报错，可用 `${VAR-default}` 给默认值（示例 `broken_scf_fixed.inp` 用的就是 `${NEWDOS-1}`）；`-E VAR=val`（或 `--set`）在命令行覆盖变量。`cp2k.psmp -e x.inp` 打印展开后的输入，是检查预处理结果最直接的办法。

### 4.2 数据文件：基组、赝势与参数文件

| 报错原文 | 原因 | 修正 |
|---|---|---|
| `The specified OLD file <BASIS_MOLOPT_UZH> cannot be opened. It does not exist.`（后面附 `Data directory path`） | 数据目录未设置或指向旧版本；文件名拼错 | `export CP2K_DATA_DIR=/path/to/cp2k/data`；在输出头部核对 `CP2K| Data directory path`；或写文件的绝对路径 |
| `The requested basis set <DZVP-MOLOPT-SR-GTH> for element <Cu> was not found in the basis set files BASIS_MOLOPT_UZH` | 该名字只存在于另一个文件（这里是旧的 `BASIS_MOLOPT`） | 改用文件中真实存在的名字，如 `DZVP-MOLOPT-PBE-GTH-q11`；或多写一行 `BASIS_SET_FILE_NAME BASIS_MOLOPT` |
| `The requested potential <GTH-PBE-q11> for element <Cu> was not found in the potential files ...` | `POTENTIAL_FILE_NAME` 不对（`GTH_POTENTIALS`、`POTENTIAL`、`POTENTIAL_UZH` 收录的不完全相同），或该元素没有这个价电子数 | `grep -nE "^ *Cu " $CP2K_DATA_DIR/POTENTIAL_UZH` 查看可用的 qN |
| 无报错，但电子数不对 | 基组的 `-qN` 与赝势的 `-qN` 不一致：**CP2K 不按名字核对二者**，电子数由赝势决定 | 每次新体系都 `grep "Number of electrons"`，与 $\sum_i Z_\mathrm{val,i} - \mathrm{CHARGE}$ 对照 |
| D3 参数文件打不开，或泛函无法识别 D3 参数 | `PARAMETER_FILE_NAME dftd3.dat` 不在当前目录或数据目录；组合/自定义泛函 | 确认数据目录；显式写 `REFERENCE_FUNCTIONAL PBE`（[第8章](08-dispersion-spin-dftu.md)） |

查基组和赝势名字最可靠的方法是直接在数据文件里搜：

```bash
grep -nE "^ *Cu " $CP2K_DATA_DIR/BASIS_MOLOPT_UZH | head   # element line: "Cu DZVP-MOLOPT-PBE-GTH-q11 ..."
grep -nE "^ *Cu " $CP2K_DATA_DIR/POTENTIAL_UZH    | head   # "Cu GTH-PBE-q1 / q11 / q19 ..."
grep -nE "^ *Cu " $CP2K_DATA_DIR/BASIS_MOLOPT      | head   # old file: lines start with a blank (" Cu DZVP-MOLOPT-SR-GTH ...")
```

`^ *` 是必要的：旧的 `BASIS_MOLOPT` 等文件里元素名前有一个空格，只写 `"^Cu "` 会什么也找不到。

成功标准：输出中 `Number of electrons` 等于预期值（例如 4 个 Cu 用 q11 赝势：44），各 KIND 的 `Basis set` 和 `Potential` 行与你写的一致。

### 4.3 结构与晶胞

| 报错原文 | 原因 | 修正 |
|---|---|---|
| `GEOMETRY wrong or EMAX_SPLINE too small` | 建立邻居表时发现原子重叠或距离过近。名字里虽有 MM 关键字 `EMAX_SPLINE`，纯 QM 计算也会报 | 用可视化软件检查结构。官方列出的三类常见原因：xyz 等格式不含晶胞信息；重新定义晶胞、切面或扩胞时出现重复原子；CIF 中的无序（部分占据）位点没有处理 |
| `An invalid set of cell vectors was specified. The cell volume is too small` 或 `The determinant det(h) is too small` | `&CELL` 写错：漏写某个矢量、单位不对（默认 Å）、三个矢量共面 | 检查 `ABC`（可加 `ALPHA_BETA_GAMMA`）或 `A`/`B`/`C` 三个矢量的数值；`MULTIPLE_UNIT_CELL` 要与 `&TOPOLOGY` 中的同名关键字一致 |
| 孤立分子能量随盒子大小漂移，或 MT 求解器报盒子太小 | `&CELL PERIODIC NONE` 与 `&POISSON PERIODIC`/`POISSON_SOLVER` 不一致；MT 要求盒子至少是电荷分布尺寸的 2 倍 | 两处 `PERIODIC` 一致，`POISSON_SOLVER MT` 或 `WAVELET`；`&TOPOLOGY/&CENTER_COORDINATES` 把分子放在盒子中央（[第6章](06-periodic-surface-kpoints.md)） |
| slab 的功函数、吸附能随真空厚度漂移 | 不对称 slab 的偶极在周期方向上相互作用 | `SURFACE_DIPOLE_CORRECTION T`、`SURF_DIP_DIR Z`，真空 ≥ 15 Å，非周期方向 k 点取 1（[第11章](11-cube-charge-workfunction.md)） |

> [!WARNING]
> 结构问题经常“不报错但结果错”：例如 xyz 坐标的单位是 bohr 却当 Å 读入、`SCALED` 坐标漏写 `SCALED`。每个新体系第一次运行时，都用 `&MOTION/&PRINT/&TRAJECTORY` 或 `FINAL_STRUCTURE`（2026.2）输出的结构在 VESTA/OVITO/VMD 中看一眼，并核对输出中 `ATOMIC COORDINATES` 段的最近邻距离。

### 4.4 电子结构与 SCF

| 报错原文 | 原因 | 修正 |
|---|---|---|
| `Use the LSD option for an odd number of electrons` | 电子数为奇数却没开自旋极化（LSD 是 UKS 的别名） | `&DFT` 中加 `UKS` 并设正确 `MULTIPLICITY`；本应闭壳层的体系出现此错，多半是少了/多了 H 原子或 `CHARGE` 写错 |
| `LSD: try to use a different multiplicity` 一类消息 | `MULTIPLICITY` 与电子数不相容 | 按化学重设 `CHARGE`/`MULTIPLICITY`（[第8章](08-dispersion-spin-dftu.md) §6） |
| `OT not possible with kpoint calculations` | 2026.2 及以前，`&OT` 不能与 `&KPOINTS` 同时使用 | 改用 `&DIAGONALIZATION` + `&MIXING`（金属再加 `&SMEAR` 与 `ADDED_MOS`）；trunk（2026.2 之后）已开始支持 OT + k 点 |
| `SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.` | 自 2024.1 起 SCF 不收敛默认中止 | 按 [第5章](05-scf-strategy.md) §10 的流程处理；`IGNORE_CONVERGENCE_FAILURE` 只作最后手段，**不要在 MD 中使用** |
| 无报错，但金属的能量和力异常、SCF 振荡 | `&SMEAR` 时 `ADDED_MOS` 太小，费米分布尾部被截断 | `ADDED_MOS` 取占据态数的 10–30%；输出中最高 MO 的占据数应接近 0 |
| `Periodic Hartree Fock calculation requested with the use of a truncated or shortrange potential. The cutoff radius ( x A) is larger than half the minimal cell dimension ( y A). This may lead to unphysical total energies. ...` | 截断库仑/短程势的 `CUTOFF_RADIUS` 超过最小晶胞尺寸的一半 | `CUTOFF_RADIUS` ≤ L/2，否则扩大超胞（[第9章](09-hybrid-admm.md) §7） |
| 无报错，但 2026.2 与旧版本的展宽结果（熵项、能量）不同 | `&SMEAR/METHOD` 的默认值：2025.2/2026.1 的输入参考为 `ENERGY_WINDOW`，2026.2 为 `GAUSSIAN` | 始终显式写 `METHOD FERMI_DIRAC`（或你需要的方法）和对应的温度/宽度 |

SCF 不收敛的处理顺序（详见第5章）：

1. 先确认不是输入问题：结构、`CHARGE`/`MULTIPLICITY`/`UKS`、电子数、基组与赝势的 qN；
2. 有能隙体系（Γ 点）：OT + `&OUTER_SCF`（内层 `MAX_SCF` 20–30，外层 8–16），`PRECONDITIONER FULL_ALL` 或 `FULL_SINGLE_INVERSE`，`MINIMIZER DIIS`（困难时换 `CG`）；
3. 金属、小带隙或 k 点：对角化 + `BROYDEN_MIXING`（`ALPHA` 0.1–0.2，`NBUFFER` 8 左右）+ `&SMEAR` + 足够的 `ADDED_MOS`；
4. 从便宜计算（小基组、低截断、PBE）的 `RESTART.wfn` 起步；
5. 依然不行，回到第 1 步怀疑结构：极端的键长、重叠原子、错误的磁性初猜，都会让 SCF 无法收敛。

### 4.5 运行环境与资源

| 报错或现象 | 原因 | 修正 |
|---|---|---|
| 输出中 `DBCSR| MPI: Number of processes 1` 出现多次，`PROGRAM STARTED AT` 也重复出现、行与行交错 | **MPI 错配**：`mpirun` 与 CP2K 链接的 MPI 不是同一套，N 个串行进程各自运行并写同一个文件 | §6 |
| `error while loading shared libraries: libxxx.so: cannot open shared object file` | 计算节点上没有加载运行环境 | 作业脚本中 `source install/cp2k_env`、`module load` 或 `conda activate` |
| `Illegal instruction` | 用 `-march=native` 编译的程序跑在指令集更老的 CPU 上 | 在计算节点上编译，或 toolchain 用 `--target-cpu=generic`（[第1章](01-install.md) §7.3） |
| 无报错，作业被杀，`sacct` 显示 `OUT_OF_MEMORY` 或 stderr 中有 `Killed`/`oom-kill` | 内存不足（常见于 HFX、大体系对角化、MO cube 输出） | §7.5：减少每节点 MPI 进程、增加线程或节点；设 `MAX_MEMORY` |
| 大体系启动不久后段错误（`Segmentation fault`、`SIGSEGV`） | 栈空间不足最常见 | 作业脚本中 `ulimit -s unlimited`；必要时 `export OMP_STACKSIZE=64m`（按需调整） |
| 程序“卡住”，输出很久不更新 | 可能在做第一步 HFX ERI 计算（正常但很慢）、内存交换、或 MPI 死锁 | `top`/`htop` 看 CPU 占用与内存；HFX 首步本来就慢，不要过早杀掉；ELPA 相关的问题可试 `PREFERRED_DIAG_LIBRARY SCALAPACK`（官方 Troubleshooting 页提到可规避 issue #4484 一类的问题） |
| 输出中出现 `************`、`NaN`、`Inf` | 星号是 Fortran 定宽格式放不下数值；NaN/Inf 是数值发散 | 星号提示数值异常大/小，要回头检查；NaN/Inf 说明计算已失控，从结构和 SCF 设置查起 |
| 报 `CPASSERT failed`，只给出源码位置 | 内部一致性检查失败 | 先确认输入完全符合手册；仍出现则带最小复现输入到 GitHub issue 或 Google Group 报告 |

### 4.6 运动层（优化、MD、NEB、振动）

| 报错或现象 | 原因 | 修正 |
|---|---|---|
| `*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***` | `MAX_ITER` 用完仍未收敛 | 看 `OPT|` 块哪个判据卡住，从 `-1.restart` 续算而不是从头开始（[第7章](07-geo-cell-opt.md)） |
| 优化中能量和力在几个值之间振荡 | 力的数值噪声大于收敛阈值 | `EPS_SCF` 1E-7、`EPS_DEFAULT` 1E-12、提高 `CUTOFF` |
| CELL_OPT/NPT 的压力不准、晶胞不断缩小 | 截断不够导致 Pulay 应力 | `CUTOFF` 明显高于单点计算（如 ≥ 600 Ry），`STRESS_TENSOR ANALYTICAL` |
| AIMD 温度失控、守恒量漂移 | 步长太大、SCF 太松、外推设置不当 | 含 H 时 `TIMESTEP` 0.5 fs，`EPS_SCF` ≤ 1E-6，`EXTRAPOLATION ASPC`（[第12章](12-aimd.md)） |
| NEB 主输出里看不到错误，任务却停了 | 报错写在各副本的日志里 | 查看 `<PROJECT>-BAND*.out`；总进程数应为 `NUMBER_OF_REPLICA × NPROC_REP`（[第13章](13-neb-dimer.md)） |
| 振动分析出现虚频 | 优化未严格收敛、SCF 不够紧 | `MAX_FORCE` ≤ 1E-4 Ha/bohr，收紧 `EPS_SCF`，`DX` 0.01 bohr 左右（[第14章](14-vibration-thermo.md)） |

---

## 5. 实战一：修复一个有 6 处错误的输入（broken_scf.inp）

### 5.1 故意出错的输入

[`broken_scf.inp`](../examples/18-troubleshooting/broken_scf.inp) 是 fcc Cu 4 原子惯用胞（a = 3.615 Å，8×8×8 k 点）的 PBE 单点计算，里面混入了实际工作中最常见的 6 类错误，用 (A)–(F) 标注：

```text
# DELIBERATELY BROKEN INPUT - debugging exercise of chapter 18. Do not use for production.
&GLOBAL
  PROJECT  cu_broken
  RUN_TYPE ENERGY
&END GLOBAL

# leftover from an old cell optimisation input
&MOTION
  &CELL_OPT
    # (A) removed in CP2K 2026.2
    TYPE DIRECT_CELL_OPT
    KEEP_ANGLES T
  &END CELL_OPT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &KPOINTS
      SCHEME MONKHORST-PACK 8 8 8
    &END KPOINTS
    &SCF
      EPS_SCF 1.0E-7
      MAX_SCF 30
      # (D) OT cannot be combined with k-points (CP2K <= 2026.2)
      &OT ON
        MINIMIZER DIIS
      &END OT
      # (F) a metal without smearing, without unoccupied states, without a good mixer
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      # (B) pre-2026.2 PDOS syntax
      &PDOS
        NLUMO -1
      &END PDOS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 3.615 3.615 3.615
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Cu  0.0  0.0  0.0
      Cu  0.0  0.5  0.5
      Cu  0.5  0.0  0.5
      Cu  0.5  0.5  0.0
    &END COORD
    &KIND Cu
      # (C) this basis name exists in BASIS_MOLOPT, not in BASIS_MOLOPT_UZH
      BASIS_SET DZVP-MOLOPT-SR-GTH
      # (E) q1 pseudopotential, but every Cu GTH basis in BASIS_MOLOPT_UZH is a q11 basis
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

（上面为便于阅读删去了文件头部的几行注释，完整文件见示例目录。）

### 5.2 按报错顺序逐个修正

```bash
cd CP2K/examples/18-troubleshooting
export CP2K_DATA_DIR=/path/to/cp2k/data
cp broken_scf.inp work.inp                    # never edit the original
cp2k.psmp -c work.inp                         # round 1: parser errors (A, B on 2026.2)
bash check_migration.sh work.inp              # text scan for removed/changed syntax
# then: fix one problem, run again
mpirun -np 4 cp2k.psmp -i work.inp -o work.out; tail -25 work.out
```

在 CP2K 2026.2 上，错误按下面的顺序暴露（解析器遇到第一个错误就中止，所以每修好一个，下一个才会出现）：

| 顺序 | 标记 | 阶段 | 典型报错（节选） | 修正 |
|---|---|---|---|---|
| 1 | A | 解析输入 | `The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE has been removed because cell optimizations now always use DIRECT_CELL_OPT.` | 删除 `TYPE`；本例是单点计算，整个 `&MOTION` 都可以删 |
| 2 | B | 解析输入 | `unknown subsection PDOS of section PRINT` | 改用 `&PRINT/&DOS/&PDOS`（2026.2），或删除 |
| 3 | C | 读取基组 | `The requested basis set <DZVP-MOLOPT-SR-GTH> for element <Cu> was not found in the basis set files BASIS_MOLOPT_UZH` | 改为 `DZVP-MOLOPT-PBE-GTH-q11` |
| 4 | D | SCF 初始化 | `OT not possible with kpoint calculations` | 删除 `&OT`，改用对角化 |
| 5 | E | 无报错 | ——（CP2K 不按名字核对 qN） | `GTH-PBE-q1` → `GTH-PBE-q11`，并核对电子数 |
| 6 | F | SCF 迭代 | 很可能 `SCF run NOT converged`，或收敛到错误的整数占据态 | `&SMEAR`（显式 `METHOD`）、`ADDED_MOS`、Broyden 混合、`MAX_SCF 100` |

在 2026.1 及更早版本上，A 和 B 不是错误，练习从 C 开始。

> [!CAUTION]
> 第 5 处（E）最危险：如果只修了 C 而赝势仍是 q1（`POTENTIAL_UZH` 中确实有 `GTH-PBE-q1`），CP2K **不会因为基组与赝势的 qN 不一致而报错**，修好 D、F 后计算可以照常结束，但每个 Cu 只有 1 个价电子（4 个原子共 4 个电子），能带、磁性、表面性质全错。这是“正常结束 ≠ 结果正确”的典型例子。养成习惯：新体系第一次计算后 `grep "Number of electrons"`。

### 5.3 修正版输入

[`broken_scf_fixed.inp`](../examples/18-troubleshooting/broken_scf_fixed.inp) 的完整内容如下。DOS 输出用预处理变量 `NEWDOS` 在新旧语法之间切换，同一个文件在 2025.x/2026.1 和 2026.2 上都能用：

```text
&GLOBAL
  PROJECT     cu_fixed
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 8 8 8
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   100
      ADDED_MOS 12                 ! 22 occupied orbitals (44 electrons); keep the top MO empty
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.2
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD                 FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
      &PRINT
        # k-point runs write cu_fixed-RESTART.kp (not .wfn)
        &RESTART ON
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
@IF ${NEWDOS-1} == 1
      # CP2K >= 2026.2: PDOS lives inside &DOS, broadened curves via &CURVE
      &DOS
        NLUMO -1
        &PDOS
        &END PDOS
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO FERMI
          &BROADEN
            WIDTH [eV] 0.2
          &END BROADEN
        &END CURVE
      &END DOS
@ENDIF
@IF ${NEWDOS-1} == 0
      # CP2K <= 2026.1: &PDOS directly under &DFT/&PRINT
      &PDOS
        NLUMO -1
      &END PDOS
@ENDIF
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 3.615 3.615 3.615
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Cu  0.0  0.0  0.0
      Cu  0.0  0.5  0.5
      Cu  0.5  0.0  0.5
      Cu  0.5  0.5  0.0
    &END COORD
    &KIND Cu
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `BASIS_SET` | `DZVP-MOLOPT-PBE-GTH-q11` | `BASIS_MOLOPT_UZH` 中 Cu 的 DZVP 基组；名字中的 q11 必须与赝势一致 |
| `POTENTIAL` | `GTH-PBE-q11` | 11 个价电子（3d¹⁰4s¹）；4 个原子共 44 个电子 |
| `&DIAGONALIZATION ALGORITHM` | `STANDARD` | k 点计算在 2026.2 及以前只能用对角化 |
| `&MIXING METHOD` | `BROYDEN_MIXING` | 金属的首选；`ALPHA` 0.1–0.2 比默认 0.4 稳定 |
| `NBUFFER` | 8 | Broyden 使用的历史步数 |
| `&SMEAR METHOD` | `FERMI_DIRAC` | 显式写出，避免 2026.2 默认值变化（`GAUSSIAN`）带来的差异 |
| `ELECTRONIC_TEMPERATURE` | 300–1000 K | 越高越好收敛，但电子熵项越大；要报告外推到 0 K 的能量或检查温度收敛 |
| `ADDED_MOS` | 12 | 22 个占据轨道 × 约 50%，保证最高轨道占据数接近 0；大体系一般取 10–30% |
| `MAX_SCF` | 100 | 对角化 + 混合通常需要更多步 |
| `EPS_DEFAULT` | 1.0E-12 | 积分与网格精度，保证能量噪声远小于 `EPS_SCF` 对应的误差 |
| `NEWDOS`（预处理变量） | 1（默认） | `-E NEWDOS=0` 时使用 2026.1 及以前的 `&PDOS` 语法 |

运行与检查：

```bash
mpirun -np 8 cp2k.psmp -i broken_scf_fixed.inp -o fixed.out                  # CP2K >= 2026.2
mpirun -np 8 cp2k.psmp -E NEWDOS=0 -i broken_scf_fixed.inp -o fixed.out      # CP2K <= 2026.1
grep "Number of electrons" fixed.out | head -2           # expect 44
grep -E "SCF run converged|Fermi energy|Electronic entropic energy" fixed.out | tail -3
ls cu_fixed-RESTART.kp
ls cu_fixed*curve*.pdos cu_fixed*curve*.dos              # 2026.2 only: broadened k-point PDOS/DOS curves
```

成功标准：电子数 44；SCF 在 `MAX_SCF` 内收敛；打印费米能与电子熵；生成 `.kp` 重启文件；在 2026.2 上生成展宽后的 DOS/PDOS 曲线文件（文件名细节见[第10章](10-dos-pdos-band.md) §4.3、§7.2，用通配符匹配）。

> [!NOTE]
> **k 点计算不会打印“smearing includes the last MO”警告。** 这条警告（`Fermi-Dirac smearing includes the last MO => Add more MOs for proper smearing.`）来自 Γ 点的占据数例程；显式 `&KPOINTS` 时占据数由另一套 k 点例程计算，2026.2 源码中没有对应的检查。所以本例只能用收敛测试判断 `ADDED_MOS` 是否足够：把 `ADDED_MOS` 加大（如 12 → 20），总能量和 `Fermi energy` 不再变化（差异远小于目标精度）即可（练习 3）。
>
> 另外，**2026.1 及以前的 `&PDOS` 不支持 k 点**：`-E NEWDOS=0` 时程序只打印一条 PDOS 未实现的警告，不写 `.pdos` 文件（[第10章](10-dos-pdos-band.md) §3）。k 点 PDOS 只能在 2026.2 的 `&DOS/&PDOS` + `&CURVE` 下得到，而且只有展宽曲线、没有逐态文件；此时输出中的 `State-resolved k-point PDOS output is not implemented yet` 警告属于正常现象。

---

## 6. MPI 环境错配与运行环境检查

### 6.1 症状

官方 Troubleshooting 页给出的对比：用 `mpirun -np 4` 正常启动时，输出中

```text
 DBCSR| MPI: Number of processes                                               4
 DBCSR| OMP: Current number of threads                                         1
```

只出现一次。如果同样的两行出现 4 次且进程数都是 1，同时 `PROGRAM STARTED AT`、`PROGRAM PROCESS ID` 等行也成倍出现并相互交错，说明实际运行的是 4 个互相独立的串行 CP2K，它们在以不确定的顺序写同一个输出文件。此时计算仍可能“正常结束”，但速度没有提升、输出混乱，对于需要写文件的任务甚至会互相覆盖。

```bash
grep -c "DBCSR| MPI: Number of processes" x.out     # must be 1 for a single FORCE_EVAL run
grep -c "PROGRAM STARTED AT" x.out                  # must be 1
```

> [!NOTE]
> NEB、振动分析等多副本任务会在**副本日志**里各自打印并行信息，主输出中的 `DBCSR| MPI` 行数仍应为 1。统计时只看主输出。

### 6.2 原因与诊断

根本原因是：启动用的 `mpirun`/`mpiexec`/`srun` 所属的 MPI 库，与 `cp2k.psmp` 链接的 MPI 库不是同一套（例如 conda 装的 CP2K 链接 conda 的 MPICH，却用集群 module 的 OpenMPI `mpirun` 启动）。

```bash
which mpirun mpiexec srun
mpirun --version 2>&1 | head -2                # Open MPI / MPICH / Intel MPI ?
ldd $(which cp2k.psmp) | grep -iE "mpi|not found"
cp2k.psmp -v | grep -i cp2kflags               # should contain "parallel" and "mpi..." flags
```

`ldd` 显示的 `libmpi.so`（OpenMPI）或 `libmpich.so`/`libmpi.so.12`（MPICH 系）必须与 `mpirun --version` 的类型一致，而且来自同一安装目录。

### 6.3 修正

- toolchain/CMake 构建：作业脚本中先 `source <CP2K_ROOT>/install/cp2k_env`（toolchain 自带的 MPI 也由它加载），再用同一套 `mpirun`；
- conda 安装：先 `conda activate`，用环境里的 `mpiexec`；不要与 module 的 MPI 混用；
- 集群 module：只加载管理员给出的 CP2K module 及其依赖，不要在 `~/.bashrc` 里额外 `export PATH`/`LD_LIBRARY_PATH` 指向其他 MPI；
- Slurm 下优先用 `srun`（需要 MPI 编译时支持 PMI/PMIx，具体以集群文档为准）。

作业脚本的通用写法见 [第1章](01-install.md) §11.3；本篇的 Slurm/PBS 模板 [`job_cp2k_slurm.sh`](../scripts/job_cp2k_slurm.sh)、[`job_cp2k_pbs.sh`](../scripts/job_cp2k_pbs.sh) 已包含环境加载、MPI×OpenMP 参数化、结束检查与重启提示（用法见脚本头部注释）。Linux 与集群基础见 [公共篇第1章](../../Common/tutorials/01-linux-hpc.md)。

---

## 7. 性能调优

### 7.1 先读 timing 报告

每次正常结束的输出末尾都有一张 timing 表，按子程序的总耗时排序。名为 `CP2K` 的行是总时间：

（示意输出，数值仅供参考）

```text
 -------------------------------------------------------------------------------
 -                                T I M I N G                                  -
 -------------------------------------------------------------------------------
 SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    0.020    0.031  612.410  612.415
 qs_energies                          1  2.0    0.001    0.002  598.205  598.210
 scf_env_do_scf                       1  3.0    0.000    0.001  590.117  590.120
 ...
```

```bash
grep -A30 "T I M I N G" x.out | head -40
grep "MEMORY| Estimated peak process memory" x.out
```

阅读要点：

- `TOTAL TIME` 是包含子调用的时间，`SELF TIME` 是不含子调用的时间；找瓶颈时看 `SELF TIME` 大的行。`&GLOBAL/&TIMINGS` 中 `SORT_BY_SELF_TIME T` 可按自身时间排序，`THRESHOLD`（默认 0.02，即 2%）控制显示的最小占比；
- `AVERAGE` 与 `MAXIMUM` 相差很大，说明进程之间负载不均衡；`&TIMINGS REPORT_MAXLOC T` 可报告最慢的进程；
- 常见热点的大致归类：稀疏矩阵乘法（名字中含 `multiply`，DBCSR/DBM）、FFT（含 `fft`）、密度映射与势积分（含 `rho`/`integrate`）、对角化（含 `diag`/`syevd`/`elpa`）、HFX（含 `hfx`）。具体子程序名随版本变化，以你的输出为准。

| 热点 | 可能的优化 |
|---|---|
| 对角化 | 大体系优先 OT（Γ 点、有能隙）；用 ELPA；k 点用 `PARALLEL_GROUP_SIZE` 调整 k 点并行 |
| FFT、网格映射 | 检查 `CUTOFF` 是否远超收敛所需；`REL_CUTOFF` 是否过高；FFTW3 是否被使用（`GLOBAL| Preferred FFT library`） |
| 稀疏矩阵乘法 | MPI 进程数取能形成接近方形进程网格的数；更换 MPI×OpenMP 布局 |
| HFX | ADMM、`SCREEN_ON_INITIAL_P`、`MAX_MEMORY` 足够（积分不必每步重算）、2026.2 的 ACE 选项、GPU（libGint） |

### 7.2 MPI×OpenMP 配比

经验规则（需要用自己的体系验证）：

| 场景 | 建议起点 | 理由 |
|---|---|---|
| GGA、Γ 点、OT，几十到几百原子 | 每核 1 个 MPI 进程，或每进程 2 线程 | DBCSR 与 FFT 的 MPI 并行效率高 |
| GGA + 大量 k 点 | 以 MPI 为主；用 `&KPOINTS/PARALLEL_GROUP_SIZE` 控制每个 k 点使用的进程数 | k 点之间几乎独立，并行效率高；该值必须能整除总进程数 |
| 杂化泛函（HFX） | 每节点较少 MPI 进程 × 4–8 线程 | `MAX_MEMORY` 按进程计，线程共享内存 |
| 内存紧张的大体系对角化 | 减少每节点 MPI 进程，增加线程 | 降低复制数据的总量 |
| NEB、振动分析 | 总进程数 = 副本数 × `NPROC_REP`；每个副本内部再按上面规则 | 副本之间几乎独立 |

运行时的环境变量：

```bash
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}   # threads per MPI rank
export OMP_PLACES=cores                            # one thread per physical core
export OMP_PROC_BIND=close                         # keep the threads of a rank together
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}  # Slurm >= 22.05: srun does not inherit --cpus-per-task
ulimit -s unlimited
srun --cpu-bind=cores cp2k.psmp -i x.inp -o x.out
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| `OMP_NUM_THREADS` | 与 `--cpus-per-task` 相同 | 不设时，有些环境会让每个进程都开满全部核的线程，导致严重超订 |
| `OMP_PLACES` / `OMP_PROC_BIND` | `cores` / `close` | 线程绑定到物理核，避免迁移 |
| `--cpu-bind=cores` | 开 | 进程绑定到核；与 MPI 库自己的绑定选项二选一 |
| 超线程 | 一般不用 | 每个物理核只放一个线程 |
| `cp2k.popt` | —— | 指向 `psmp` 的别名，但强制每进程 1 线程（[第1章](01-install.md) §3.2） |

> [!WARNING]
> 总线程数超过物理核数（超订）是最常见的性能杀手：例如 `--ntasks-per-node=64` 却忘了设 `OMP_NUM_THREADS`。输出中 `DBCSR| OMP: Current number of threads` 必须等于你想要的值。

### 7.3 MPI×OpenMP 基准脚本

[`job_mpi_omp_hybrid.sh`](../examples/18-troubleshooting/job_mpi_omp_hybrid.sh) 在一个整节点上依次运行若干组合（默认 `N×1`、`N/2×2`、`N/4×4`、`N/8×8`），每个组合放在单独的 `bench_<布局>/` 目录中，最后汇总到 `bench_summary.dat`。它先用 `cp2k.psmp -c` 检查输入，并核对输出中报告的 MPI 进程数，错配时在 `status` 列标出：

```bash
sbatch job_mpi_omp_hybrid.sh                                       # default input: broken_scf_fixed.inp
sbatch job_mpi_omp_hybrid.sh si_cellopt_2026_2.inp "32x1 16x2 8x4"  # your input and layouts
cat bench_summary.dat
```

（示意输出，数值仅供参考）

```text
layout      ranks  threads     walltime_s status
64x1           64        1          301.2 ok
32x2           32        2          288.7 ok
16x4           16        4          322.5 ok
8x8             8        8          401.9 ok
```

做基准的原则：用与生产计算**同体系、同参数**的短任务（几个 SCF 步、几步 MD/优化），而不是玩具分子；每个组合至少跑两次看波动；比较的是总墙钟时间（`CP2K` 行的最后一列）。多节点强扩展测试可在此基础上把 `--nodes` 改为 1、2、4……，计算 §2.2 的并行效率 $E(N)$，一般在 $E$ 降到 0.6–0.7 以下时就不值得继续加节点了。

### 7.4 对角化库：PREFERRED_DIAG_LIBRARY

`&GLOBAL` 中的相关关键字（2026.2 输入参考）：

| 关键词 | 默认值 | 可选值 / 说明 |
|---|---|---|
| `PREFERRED_DIAG_LIBRARY` | `ELPA` | `ELPA`、`SCALAPACK`（别名 `SL`）、`CUSOLVER`（cuSOLVERMp）、`DLAF`（DLA-Future）；所选库不可用时退回 ScaLAPACK |
| `ELPA_NEIGVEC_MIN` | 64 | 本征矢数小于此值时用 ScaLAPACK（小矩阵 ELPA 不划算） |
| `ELPA_KERNEL` | `AUTO` | 一般不改；有 GPU 时有 `NVIDIA_GPU`、`AMD_GPU` 等选项 |
| `ELPA_ONE_STAGE` | `F` | 使用一阶段求解器 |
| `DLAF_NEIGVEC_MIN` | 1024 | DLA-Future 的最小规模 |
| `&FM_DIAG_SETTINGS` `PARAMETER_A`/`PARAMETER_X` | 4 / 60 | 决定参与对角化的进程数，一般不改 |

```text
&GLOBAL
  PROJECT  big_diag
  RUN_TYPE ENERGY
  PREFERRED_DIAG_LIBRARY ELPA     ! default; switch to SCALAPACK to work around ELPA problems
&END GLOBAL
```

```bash
grep -E "GLOBAL\| Preferred (diagonalization lib|FFT library)" x.out
cp2k.psmp -v | grep -i cp2kflags | tr ' ' '\n' | grep -iE "elpa|dlaf|cusolver|scalapack"
```

建议：程序带 ELPA 时保持默认；遇到对角化阶段的崩溃、结果异常或明显变慢，先改成 `SCALAPACK` 对比一次，若问题消失，就在作业中固定使用 ScaLAPACK，并向开发者报告。

### 7.5 内存：HFX、NPROC_REP 与其他

| 内存大户 | 控制方法 |
|---|---|
| HFX 的 ERI 存储 | `&HF/&MEMORY/MAX_MEMORY`（每个 MPI 进程，MiB，默认 512）；`MAX_DISK_SPACE` 可把一部分存到磁盘（`STORAGE_LOCATION` 指定目录）；减少每节点进程数、增加线程 |
| 多副本任务（NEB、振动分析） | `&BAND/NPROC_REP`、`&VIBRATIONAL_ANALYSIS/NPROC_REP`：每个副本的进程数。每个副本相当于一个只用 `NPROC_REP` 个进程的单点计算，`NPROC_REP` 越小，每个进程分到的数据越多、内存越高；内存不足时减少同时运行的副本数、增大 `NPROC_REP` |
| cube 输出 | 大体系的 `MO_CUBES`、`E_DENSITY_CUBE` 可能很大；用 `STRIDE` 降低分辨率，只输出需要的轨道 |
| 输入结构本身 | `&GLOBAL/SAVE_MEM T` 释放暂时不用的输入段，节省少量内存 |

杂化泛函的内存估算步骤（详见[第9章](09-hybrid-admm.md) §10）：

```bash
grep -A12 "HFX_MEM_INFO" x.out | head -20          # memory report of the HFX module
grep "calculated on the fly" x.out                   # >0: MAX_MEMORY too small, integrals recomputed each step
grep "MEMORY| Estimated peak process memory" x.out
```

```text
&XC
  &HF
    &MEMORY
      MAX_MEMORY 6000            ! MiB per MPI rank: (node RAM * 0.8) / ranks_per_node - other usage
    &END MEMORY
    &SCREENING
      EPS_SCHWARZ          1.0E-6
      SCREEN_ON_INITIAL_P  T     ! needs a converged GGA wavefunction as SCF_GUESS RESTART
    &END SCREENING
  &END HF
&END XC
```

成功标准：`calculated on the fly` 为 0 或很小；`Estimated peak process memory` × 每节点进程数小于节点物理内存的约 80%；第一步之后每个 SCF 步的耗时明显低于第一步。

### 7.6 GPU 构建与运行

CP2K 的 GPU 加速覆盖稀疏矩阵乘法（DBCSR/DBM）、网格映射（`&GLOBAL/&GRID BACKEND`，默认 `AUTO`）、部分 FFT、对角化（cuSOLVERMp、DLA-Future）以及 2026.2 新增的 libGint（CUDA 上的 HFX，`&HF/HFX_LIBRARY LIBGINT`）。GPU 版本必须在编译时开启：

| 构建方式 | GPU 相关设置（2026.2） | 参考 |
|---|---|---|
| CMake | `-DCP2K_USE_ACCEL=CUDA` 加 `-DCP2K_WITH_GPU=A100`（按显卡型号） | [第1章](01-install.md) §7.6 |
| Spack 脚本 | `./make_cp2k.sh -gm A100` | [第1章](01-install.md) §8 |
| 容器 | 官方 CUDA 镜像（见 cp2k-containers） | [第1章](01-install.md) |

运行建议：

- 一般每块 GPU 配 1–2 个 MPI 进程，其余 CPU 核用 OpenMP 线程；
- 用 `CUDA_VISIBLE_DEVICES` 或调度系统的 `--gpus-per-task` 把进程映射到 GPU；
- `cp2k.psmp -v` 的 `cp2kflags` 中应出现 CUDA/HIP 相关的 offload 标志；没有就说明这是 CPU 版本；
- 小体系或纯 GGA 中等体系在 GPU 上未必更快，同样要做 §7.3 的基准测试。

### 7.7 不改布局也能提速的设置

| 设置 | 做法 | 注意 |
|---|---|---|
| `CUTOFF`/`REL_CUTOFF` | 用收敛测试得到的最小值（[第4章](04-cutoff-convergence.md)） | 晶胞优化、应力需要更高截断 |
| 基组 | 测试阶段用 DZVP，确认结论后再用 TZV2P 检验 | 旧 `BASIS_MOLOPT` 中的 `-SR` 短程基组在凝聚相中更稳定、更快；`BASIS_MOLOPT_UZH` 不用这个后缀 |
| SCF 初猜 | MD/优化用 `EXTRAPOLATION ASPC`；新计算从低精度 `.wfn` 起步 | 2026.2 起 k 点也支持波函数外推 |
| `EPS_DEFAULT` | 1E-10（默认）到 1E-12 | 过紧会变慢，过松会让力有噪声 |
| OT vs 对角化 | 大的有能隙 Γ 点体系优先 OT | k 点或金属只能对角化（2026.2） |
| 输出 | 降低 `PRINT_LEVEL`、减少每步 cube/轨迹输出频率（`&EACH`） | 大体系的 I/O 不能忽视 |

---

## 8. 重启策略

### 8.1 重启相关文件一览

| 文件 | 由谁写出 | 格式 | 用途 |
|---|---|---|---|
| `<P>-1.restart` | `&MOTION/&PRINT/&RESTART`（默认开） | 文本，完整输入 | 续算优化/MD/NEB：直接作为输入，或被 `&EXT_RESTART` 读取 |
| `<P>-1.restart.bak-N` | 同上，`BACKUP_COPIES` 控制份数 | 文本 | 最近几步的备份 |
| `<P>-1_<步数>.restart` | `&MOTION/&PRINT/&RESTART_HISTORY`（按 `&EACH` 的间隔） | 文本 | 保留历史快照，例如每 200 步一份 |
| `<P>-RESTART.wfn` | `&SCF/&PRINT/&RESTART`（Γ 点，默认开） | 二进制 | `SCF_GUESS RESTART` |
| `<P>-RESTART.wfn.bak-N` | 同上（`BACKUP_COPIES` 默认 1） | 二进制 | `SCF_GUESS HISTORY_RESTART` 会用到 |
| `<P>-RESTART.kp` | 同上（显式 `&KPOINTS`） | 二进制 | k 点计算的重启（存的是各自旋、各晶胞像的密度矩阵） |

用 `file` 命令可以区分：运动重启文件显示 `ASCII text`，波函数文件不是文本。

### 8.2 运动重启：直接运行或 &EXT_RESTART

**方法一**：把 `-1.restart` 当输入直接运行。它包含原输入的全部设置（以及显式写出的默认值）、最新坐标/速度/晶胞和 `STEP_START_VAL`：

```bash
cp -p si8_cellopt-1.restart si8_cellopt-1.restart.keep   # keep a copy, it will be overwritten
cp2k.psmp -i si8_cellopt-1.restart -o si8_cellopt-r1.out
```

注意项目名不变，新的轨迹和能量文件会**追加或覆盖**到原文件上；重要数据先备份到别的目录。

**方法二**：保留（或修改）原输入，加 `&EXT_RESTART` 只读取需要的部分。这种方式适合“续算时顺便改参数”，例如收紧收敛判据、换项目名、从优化结果起跑 MD。示例 [`si_cellopt_restart.inp`](../examples/18-troubleshooting/si_cellopt_restart.inp) 续算 `si_cellopt_2026_2.inp` 的晶胞优化，同时读取 `.kp` 波函数：

```text
&GLOBAL
  PROJECT     si8_cellopt_r2
  RUN_TYPE    CELL_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&EXT_RESTART
  RESTART_FILE_NAME si8_cellopt-1.restart
  RESTART_DEFAULT   F
  RESTART_POS       T
  RESTART_CELL      T
  RESTART_COUNTERS  T
&END EXT_RESTART

&MOTION
  &CELL_OPT
    OPTIMIZER          BFGS
    KEEP_ANGLES        T
    EXTERNAL_PRESSURE  [bar] 1.0
    PRESSURE_TOLERANCE [bar] 50.0
    MAX_ITER           100
  &END CELL_OPT
  &PRINT
    &RESTART ON
      BACKUP_COPIES 2
    &END RESTART
    &RESTART_HISTORY ON
      &EACH
        CELL_OPT 10
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME   BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME   POTENTIAL_UZH
    WFN_RESTART_FILE_NAME si8_cellopt-RESTART.kp
    &MGRID
      CUTOFF     600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 4
    &END KPOINTS
    &SCF
      SCF_GUESS RESTART
      EPS_SCF   1.0E-7
      MAX_SCF   80
      ADDED_MOS 8
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA  0.3
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 5.40 5.40 5.40
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Si  0.00  0.00  0.00
      Si  0.00  0.50  0.50
      Si  0.50  0.00  0.50
      Si  0.50  0.50  0.00
      Si  0.25  0.25  0.25
      Si  0.25  0.75  0.75
      Si  0.75  0.25  0.75
      Si  0.75  0.75  0.25
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
  &PRINT
    &STRESS_TENSOR ON
    &END STRESS_TENSOR
  &END PRINT
&END FORCE_EVAL
```

`&CELL` 和 `&COORD` 仍要写（解析器需要它们），但 `RESTART_CELL T`、`RESTART_POS T` 会用重启文件中的值覆盖。`&EXT_RESTART` 的常用开关（2026.2 输入参考）：

| 关键词 | 默认 | 说明 |
|---|---|---|
| `RESTART_FILE_NAME` | —— | 要读取的重启文件（任何 CP2K 输入格式的文件都可以） |
| `RESTART_DEFAULT` | `T` | 所有未显式设置的开关的默认值；设为 `F` 后只恢复你显式打开的项 |
| `RESTART_POS` | 随 `RESTART_DEFAULT` | 原子坐标 |
| `RESTART_VEL` | 同上 | 速度（MD） |
| `RESTART_CELL` | 同上 | 晶胞（以及参考晶胞） |
| `RESTART_COUNTERS` | 同上 | 优化/MD 的步数计数器（`STEP_START_VAL`）。换项目重新编号时设 `F` |
| `RESTART_THERMOSTAT` / `RESTART_BAROSTAT` / `RESTART_BAROSTAT_THERMOSTAT` | 同上 | 恒温器、恒压器状态；更换系综或恒温器类型时设 `F` |
| `RESTART_RANDOMG` | 同上 | 随机数发生器状态（Langevin、CSVR 等随机恒温器的可重复性） |
| `RESTART_AVERAGES` | 同上 | MD 的累计平均量 |
| `RESTART_BAND` | 同上 | NEB 各副本的坐标和速度 |
| `RESTART_DIMER` | 同上 | Dimer 方法的状态 |
| `RESTART_METADYNAMICS` / `RESTART_WALKERS` | 同上 | 已沉积的 hills / 多 walker 信息（元动力学续算必须恢复） |
| `RESTART_CONSTRAINT` | 同上 | 约束/限制信息 |
| `CUSTOM_PATH` | —— | 把重启文件中指定路径的整个输入段复制到当前输入 |
| `BINARY_RESTART_FILE_NAME` | —— | 配合 `&RESTART/SPLIT_RESTART_FILE` 写出的二进制部分，用于特别大的体系 |

常见场景：

| 场景 | `&EXT_RESTART` 设置 |
|---|---|
| 同条件续跑 MD/优化 | 只写 `RESTART_FILE_NAME`（全部恢复） |
| 从 GEO_OPT/CELL_OPT 结果开始 MD | `RESTART_DEFAULT F`，`RESTART_POS T`，`RESTART_CELL T` |
| NVT 平衡后转 NPT 生产 | 恢复 `POS`、`VEL`、`CELL`、`THERMOSTAT`；`RESTART_BAROSTAT F`；生产段重新计数时 `RESTART_COUNTERS F` |
| 元动力学续算 | 默认全部恢复；确保 `RESTART_METADYNAMICS` 没有被关掉 |

> [!TIP]
> [`restart_update.py`](../scripts/restart_update.py) 可以把 `-1.restart` 或 `-pos-1.xyz` 最后一帧中的坐标/晶胞写回原输入或 xyz，用于“改了很多设置、只想要最新结构”的情况：`python ../scripts/restart_update.py -h` 查看用法。

### 8.3 RESTART_HISTORY：保留历史快照

`-1.restart` 只有一份（加 `BACKUP_COPIES` 份备份），会不断被覆盖。长 MD 或可能需要回退的优化，应打开 `RESTART_HISTORY`：

```text
&MOTION
  &PRINT
    &RESTART ON
      BACKUP_COPIES 1
      &EACH
        MD 1
      &END EACH
    &END RESTART
    &RESTART_HISTORY ON
      &EACH
        MD 500
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

这样每 500 步额外写一份 `<P>-1_<步数>.restart`。`&EACH/MD 1` 让 `-1.restart` 每步更新，任务被意外杀掉时损失最小；对极大的体系，每步写文本重启文件本身会成为 I/O 负担，可以适当放大间隔。

### 8.4 波函数重启：.wfn 与 .kp

```text
&DFT
  WFN_RESTART_FILE_NAME project-RESTART.wfn    ! or project-RESTART.kp for explicit &KPOINTS
  &SCF
    SCF_GUESS RESTART
  &END SCF
&END DFT
```

- 不写 `WFN_RESTART_FILE_NAME` 时，按 `PROJECT` 和 `&SCF/&PRINT/&RESTART/FILENAME` 推断文件名；
- 文件不存在时**不会中止**，而是打印警告 `User requested to restart the wavefunction from the file named: ... This file does not exist. ... Calculation continues using ATOMIC GUESS.` 后从原子初猜开始。续算后第一步 SCF 步数异常多时，先 `grep -i "does not exist"`；
- `SCF_GUESS HISTORY_RESTART` 读取 `.wfn` 及其 `.bak-N` 序列，配合 `&QS/EXTRAPOLATION` 使用。

兼容性规则（官方 Restarting 页）：

| 不一致的量 | Γ 点（`.wfn`） | 显式 k 点（`.kp`） |
|---|---|---|
| 原子数 | 拒绝，退回其他初猜 | 拒绝 |
| 自旋数 | 仅当重启文件自旋数少于输入、且输入多重度为 1 时允许 | 拒绝 |
| 基函数数目（换基组） | 可以 | 拒绝 |
| MO 数目 | 截断或补零；但重启文件的占据数多于输入可用 MO 时拒绝 | —— |
| 晶胞像数目 | —— | 允许 |

`.wfn` 与 `.kp` 由两套独立的代码读写，**不能互换**：Γ 点计算的 `.wfn` 不能给 k 点计算用，反之亦然。k 点计算中 `FULL_GRID ON` 写出的 `.kp` 包含全部 k 点，重启灵活性最大；官方文档说明不同 `FULL_GRID` 设置之间可以互相重启。

> [!CAUTION]
> 波函数文件是未格式化的二进制，**不保证跨版本、跨机器兼容**。升级 CP2K 后，旧的 `.wfn`/`.kp` 若读入失败或 SCF 行为异常，直接改回 `SCF_GUESS ATOMIC` 重新开始。运动重启文件是文本，跨版本一般可用，但要先用新版本 `cp2k.psmp -c` 检查（例如旧的 `-1.restart` 中含有 `TYPE DIRECT_CELL_OPT`，在 2026.2 上同样会报错）。

### 8.5 防止“白跑”：WALLTIME 与作业脚本

- `&GLOBAL/WALLTIME`（秒或 `HH:MM:SS`）设为比调度系统时限略短（例如短 10–20 分钟），CP2K 在运行中的检查点（如每个 MD/优化步）发现超时后停止，并正常收尾写出重启文件和 timing 表；
- 作业脚本结束时检查 `PROGRAM ENDED AT`，没有就打印续算提示（本篇 [`job_cp2k_slurm.sh`](../scripts/job_cp2k_slurm.sh) 已经这样做）；
- 每次续算都换一个输出文件名（`-o x-r1.out`、`x-r2.out`），不要覆盖前一段的输出；
- MD 续算后检查温度、势能、守恒量在接口处是否连续（无跳变），这是重启是否正确的最直接证据。

---

## 9. 2025.x → 2026.2 迁移

### 9.1 迁移清单（摘要）

完整清单（17 项，含构建、Python 工具）见 [migrate_2025_to_2026_2.md](../examples/18-troubleshooting/migrate_2025_to_2026_2.md)，依据是 v2026.1、v2026.2 的发布说明和对应版本的输入参考。最关键的几项：

| # | 变化 | 版本 | 旧输入在新版本上的现象 | 迁移方法 |
|---|---|---|---|---|
| 1 | 删除 `&MOTION/&CELL_OPT/TYPE`（PR #5257） | 2026.2 | `The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE has been removed because cell optimizations now always use DIRECT_CELL_OPT.` | 删掉 `TYPE` 行。旧默认值本来就是 `DIRECT_CELL_OPT`，删掉后在旧版本上行为不变（[第7章](07-geo-cell-opt.md) §8） |
| 2 | DOS/PDOS 输入重构（PR #5326） | 2026.2 | `unknown subsection PDOS of section PRINT` | `&DFT/&PRINT/&DOS` 内嵌 `&PDOS`（`COMPONENTS`）、`&LDOS`、`&R_LDOS`；`NLUMO` 移到 `&DOS`；展宽曲线用 `&CURVE` + `&BROADEN`（[第10章](10-dos-pdos-band.md)） |
| 3 | REFTRAJ 删除 `EVAL_ENERGY_FORCES`/`EVAL_FORCES`（PR #5401） | 2026.2 | `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ` | `&MOTION/&MD/&REFTRAJ` 中写 `EVAL ENERGY` 或 `EVAL ENERGY_FORCES`（可选值 `NONE`/`ENERGY`/`ENERGY_FORCES`） |
| 4 | `&SMEAR/METHOD` 默认值 | 2026.2 | 不报错，结果可能变化 | 显式写 `METHOD` |
| 5 | 删除 Makefile 与 arch 文件，只能用 CMake；DBCSR 成为必需的外部依赖 | 2026.1 | 旧的 `make ARCH=... VERSION=psmp` 不再可用 | toolchain + `build_cp2k.sh`（2026.2 新增），或 Spack 的 `make_cp2k.sh`（2026.2 新增），或手动 CMake |
| 6 | 删除 Cython Python 绑定（PR #5541） | 2026.2 | 依赖旧绑定的脚本无法导入 | 改用 ASE 的 CP2K 计算器（shell 模式）或 AiiDA 等 |
| 7 | 不再支持 GCC 8；FFTW3 接口将来可能成为硬依赖 | 2026.2 | 旧编译器构建失败 | 用较新的 GCC，带 FFTW3（或 MKL/AOCL）编译 |

反过来，**2026.2 的新写法在 2026.1 及以前同样会失败**（例如 `&DOS` 内嵌 `&PDOS` 在旧版本上报 `unknown subsection PDOS of section DOS`）。教程、脚本和组内模板要注明适用版本。

### 9.2 一对可以 diff 的输入

示例目录中的 [`legacy_si_cellopt_2025.inp`](../examples/18-troubleshooting/legacy_si_cellopt_2025.inp)（≤ 2026.1）与 [`si_cellopt_2026_2.inp`](../examples/18-troubleshooting/si_cellopt_2026_2.inp)（2026.2）是同一个 Si 晶胞优化，除文件头注释外的全部差异如下（`diff legacy_si_cellopt_2025.inp si_cellopt_2026_2.inp` 的输出节选）：

```diff
14d15
<     TYPE               DIRECT_CELL_OPT
57,60c58,70
<       &PDOS
<         NLUMO      -1
<         COMPONENTS
<       &END PDOS
---
>       &DOS
>         NLUMO -1
>         &PDOS
>           COMPONENTS
>         &END PDOS
>         &CURVE
>           ENERGY_UNIT EV
>           ENERGY_ZERO HOCO
>           &BROADEN
>             WIDTH [eV] 0.1
>           &END BROADEN
>         &END CURVE
>       &END DOS
```

| 新增关键词（2026.2） | 本例取值 | 说明 |
|---|---|---|
| `&DOS/NLUMO` | -1 | 包含全部空轨道（原来写在 `&PDOS` 中） |
| `&PDOS/COMPONENTS` | 开 | 按角动量分量输出 |
| `&CURVE/ENERGY_UNIT` | `EV` | 展宽曲线的能量单位（`HARTREE`/`EV`） |
| `&CURVE/ENERGY_ZERO` | `HOCO` | 能量零点：`AUTO`/`ABSOLUTE`/`FERMI`/`HOCO`；半导体常用价带顶 |
| `&BROADEN/WIDTH` | 0.1 eV | 展宽的半高全宽（FWHM），默认 0.1 eV |

> [!NOTE]
> 这个 Si 例子用了 k 点。旧写法在 2026.1 及以前虽然能通过解析、计算也能完成，但旧 `&PDOS` 不支持 k 点，只打印一条警告、不写 `.pdos` 文件；新写法在 2026.2 上输出的是展宽后的 k 点 PDOS 曲线（见[第10章](10-dos-pdos-band.md) §7.2）。迁移不只是“换个段名”，有时输出内容本身也变了。

### 9.3 批量扫描旧输入：check_migration.sh

[`check_migration.sh`](../examples/18-troubleshooting/check_migration.sh) 是一个基于文本的启发式扫描（awk），按段的嵌套关系查找 2026.2 中删除或改动的写法：`CELL_OPT/TYPE`、`&DFT/&PRINT/&PDOS`、REFTRAJ 的 `EVAL_*`（ERROR），`&SMEAR` 缺少 `METHOD`、`&OT` 与 `&KPOINTS` 共存、`&XTB`/`&DFTB` 的 `DO_EWALD`（WARN）。它不展开 `@INCLUDE`，`@IF` 块内的发现只作 INFO。

```bash
bash check_migration.sh legacy_si_cellopt_2025.inp si_cellopt_2026_2.inp
bash check_migration.sh $(find ~/projects -name "*.inp")      # scan a whole project tree
```

（示意输出，格式与脚本一致）

```text
legacy_si_cellopt_2025.inp:14: [ERROR] &MOTION/&CELL_OPT/TYPE was removed in 2026.2 (PR #5257); delete this line (DIRECT_CELL_OPT is now the only mode)
legacy_si_cellopt_2025.inp:57: [ERROR] &DFT/&PRINT/&PDOS is the pre-2026.2 syntax; 2026.2 expects &PRINT/&DOS/&PDOS (PR #5326)
legacy_si_cellopt_2025.inp: 2 error(s), 0 warning(s)
si_cellopt_2026_2.inp: 0 error(s), 0 warning(s)
```

扫描只是辅助，**最终以目标版本的 `cp2k.psmp -c` 为准**。Multiwfn、旧版 pycp2k、cp2k-input-tools 等生成器产生的输入，也要在 2026.2 下先 `-c` 再提交。

### 9.4 同一套输入兼容新旧版本

需要在新旧两个版本上运行同一个输入时（例如组里的集群还没升级），用预处理器分支：

```text
&PRINT
@IF ${NEWDOS-1} == 1
  &DOS
    NLUMO -1
    &PDOS
    &END PDOS
  &END DOS
@ENDIF
@IF ${NEWDOS-1} == 0
  &PDOS
    NLUMO -1
  &END PDOS
@ENDIF
&END PRINT
```

```bash
cp2k.psmp -i x.inp -o x.out                 # 2026.2 (NEWDOS defaults to 1)
cp2k.psmp -E NEWDOS=0 -i x.inp -o x.out     # 2025.x / 2026.1
```

`CELL_OPT/TYPE` 不需要分支：直接删除即可兼容所有版本。

### 9.5 构建与作业脚本的迁移

| 旧流程（≤ 2025.2） | 新流程（2026.2） |
|---|---|
| toolchain → `cp install/arch/local.* ../../arch/` → `make -j ARCH=local VERSION=psmp` → `exe/local/cp2k.psmp` | toolchain → `./build_cp2k.sh -j N`（可加 `--prefix /abs/path`；2026.2 的脚本**没有** `--preset` 选项，CPU 优化目标由 toolchain 的 `--target-cpu` 决定）→ `install/bin/cp2k.psmp` |
| `source tools/toolchain/install/setup` | `source <CP2K_ROOT>/install/cp2k_env` |
| 自写 arch 文件 | 手动 CMake：`cmake -S . -B build ... && cmake --build build --target install` |
| —— | Spack：`./make_cp2k.sh -j N`（可加 `-mpi openmpi`、`-gm A100` 等） |

> [!NOTE]
> 网上教程中的 `build_cp2k.sh --preset native-gnu-x86_64` 写法属于 2026.2 之后的开发版（trunk）。具体选项一律以你手上版本的 `./build_cp2k.sh --help` 为准（[第1章](01-install.md) §7.4）。

作业脚本迁移：可执行文件路径从 `exe/<ARCH>/` 变为 `install/bin/`；环境加载从 `setup` 变为 `cp2k_env`；在脚本里保存一份 `cp2k.psmp -v` 的输出，记录版本与 `cp2kflags`。

### 9.6 升级后的结果复现

1. 选 2–3 个代表性计算（一个分子、一个周期体系、一个你的主力体系），在旧版本与 2026.2 上用**同一份输入**（必要时用 §9.4 的分支）各跑一遍；
2. 比较总能量、力、能隙、SCF 步数；差异应远小于你的目标精度，若不是，先查发布说明中与该功能相关的改动（默认值变化、算法修复，例如 2026.1 修复了 GEO_OPT CG 的 3PNT 线搜索和振动分析中的力常数）；
3. 数据文件也要固定：据 2026.x 的数据目录说明，`data/` 中提供带版本号的 `BASIS_MOLOPT_UZH_2026.2`、`POTENTIAL_UZH_2026.2` 等快照（先 `ls $CP2K_DATA_DIR | grep UZH_` 确认你的安装里有没有）。长期项目建议在输入中写带版本号的文件名，或把所用数据文件复制到项目目录并记录来源，避免数据文件更新导致结果漂移；
4. 在项目记录中写下 CP2K 版本、`cp2kflags`、数据文件名与版本。

---

## 10. 结果可信度自检

“程序正常结束”只说明没有崩溃。投稿或得出结论前，用下面的清单逐项确认。

### 10.1 数值参数收敛表

把每个参数的测试结果整理成表格，放进论文的 SI。建议的测试范围与判据：

| 参数 | 测试方法 | 常用判据 | 工具 |
|---|---|---|---|
| `CUTOFF` | 固定 `REL_CUTOFF`，从 300 Ry 起每次加 50–100 Ry | 能量变化 < 1 meV/atom；力 < 1 meV/Å 量级；晶胞优化看应力 | [`conv_cutoff.py`](../scripts/conv_cutoff.py)（[第4章](04-cutoff-convergence.md)） |
| `REL_CUTOFF` | 固定 `CUTOFF`，40–80 Ry | 同上；各级网格的 `count for grid` 分布合理 | `conv_cutoff.py` |
| k 点 | 逐步加密 | 能量 < 1 meV/atom；金属还要看费米能、DOS | [`conv_kpoints.py`](../scripts/conv_kpoints.py)（[第6章](06-periodic-surface-kpoints.md)） |
| 基组 | DZVP → TZV2P | 目标性质（能垒、吸附能）变化小于需要的精度 | 手动 |
| 展宽 | 降低 `ELECTRONIC_TEMPERATURE` | 0 K 外推能量、力的变化 | 手动 |
| 真空层 / 盒子 | 15 → 20 → 25 Å | 功函数、吸附能 | [`cube_tools.py`](../scripts/cube_tools.py)（[第11章](11-cube-charge-workfunction.md)） |
| 超胞大小 | 2×2 → 3×3 等 | 缺陷/吸附物之间的相互作用 | 手动 |
| `EPS_SCF` / `EPS_DEFAULT` | 收紧一个量级 | 力和能量差不变 | 手动 |
| HFX `CUTOFF_RADIUS` | 不同 $R_c$ 或不同超胞 | 能隙、反应能 | [第9章](09-hybrid-admm.md) |

用 `python ../scripts/conv_cutoff.py -h`、`python ../scripts/conv_kpoints.py -h` 查看两个收敛测试脚本的用法；批量汇总多个输出时用 `cp2k_outparse.py --json` 输出再整理成表。

### 10.2 物理合理性

| 检查项 | 做法 |
|---|---|
| 电子数、自旋 | `Number of electrons`、`MULTIPLICITY`；UKS 计算查看 `<S^2>` 和 Mulliken 自旋布居是否符合预期的磁性构型 |
| 结构 | 优化后的键长与实验/文献的典型值比较；检查有没有原子“飞出”表面或穿过周期边界 |
| 优化与频率 | 极小点无虚频；过渡态只有一个虚频且振动方向对应反应坐标（[第14章](14-vibration-thermo.md)） |
| MD | 守恒量漂移、温度分布、平衡时间是否足够 |
| 能量分解 | 吸附能、反应能用**同一套**参数（同截断、同基组、同 k 点密度、同泛函与色散）计算的能量相减 |
| BSSE | 高斯基组存在基组重叠误差；弱相互作用的吸附能/结合能用 counterpoise 校正或更大基组检验 |

### 10.3 与文献和其他程序对比

- **与实验比**：先弄清楚泛函本身的系统误差。例如 PBE 通常略为高估晶格常数（Si 的 PBE 晶格常数在 5.47 Å 左右，实验约 5.43 Å，参见 Haas 等人 2009 年对多种固体的系统比较），低估带隙；与实验差 1% 的晶格常数不一定是“算错了”；
- **与同泛函的文献比**：对比时确认泛函、色散校正、赝势/全电子、自旋与 DFT+U 设置一致；
- **与平面波程序（VASP、Quantum ESPRESSO）比**：CP2K 用高斯基组 + GTH 赝势，与 PAW/超软赝势的平面波结果之间存在基组和赝势带来的差异。比较**相对量**（形成能、吸附能、能垒、晶格常数）比比较总能量有意义（不同程序的总能量零点本来就不同）；差异的大小取决于基组和赝势，差异明显时先在 CP2K 一侧换更大的基组（DZVP → TZV2P）检查结果是否向平面波结果靠拢，再检查赝势的价电子数与对方是否相当；
- **公开的验证基准**：Lejaeghere 等人的 Δ 测试和 Bosoni 等人 2024 年的验证工作，系统比较了多个 DFT 程序在元素晶体上的状态方程差异，可作为“不同程序之间应当一致到什么程度”的参照。UZH 基组与赝势的验证数据见 DCM-Uni-Paderborn/UZH-protocol 仓库。

> [!IMPORTANT]
> 报告结果时写清楚：CP2K 版本、泛函与色散、基组与赝势（含 qN 和数据文件版本）、`CUTOFF`/`REL_CUTOFF`、k 点或 Γ 点、展宽方法与温度、`EPS_SCF`、优化收敛判据。这些是别人复现你结果的最少信息。

---

## 参数选择与收敛建议

| 项目 | 建议 | 说明 |
|---|---|---|
| 新输入第一次运行 | `cp2k.psmp -c` → 小规模单点 → 正式任务 | 语法错误在登录节点几秒就能发现 |
| 输出检查 | `PROGRAM ENDED AT`、警告数、电子数、`DBCSR| MPI` 只出现一次 | `quick_diagnose.sh` 一次看完 |
| SCF 不收敛 | 先查结构与电子数，再按 OT/对角化分支调整 | `IGNORE_CONVERGENCE_FAILURE` 只作最后手段 |
| 金属 | 对角化 + Broyden + `&SMEAR`（显式 `METHOD`）+ `ADDED_MOS` 10–30% | 检查最高 MO 占据数接近 0 |
| 并行布局 | 以 MPI 为主起步；HFX、内存紧张时少进程多线程 | 用 `job_mpi_omp_hybrid.sh` 实测 |
| 对角化库 | 默认 ELPA；异常时对比 `SCALAPACK` | `GLOBAL| Preferred diagonalization lib.` 核对 |
| HFX 内存 | `MAX_MEMORY` 按（节点内存 × 0.8）/ 每节点进程数估算 | `calculated on the fly` 应为 0 |
| 多副本任务 | 总进程 = 副本数 × `NPROC_REP` | 报错看副本日志 |
| 重启 | `-1.restart` 每步更新 + `RESTART_HISTORY` 定期快照；`.wfn`/`.kp` 不互换 | 续算后检查连续性 |
| 版本 | 生产计算固定一个版本；升级后先复现 2–3 个代表性计算 | 输入、作业脚本、数据文件一起迁移 |

---

## 常见问题与报错

**1. `The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE has been removed because cell optimizations now always use DIRECT_CELL_OPT.`**
原因：2026.2 删除了 `&CELL_OPT/TYPE`，旧输入、旧 `-1.restart` 或旧生成器写出的输入中还有这一行。
解决：删除 `TYPE` 行。若原来用的是 `TYPE GEO_OPT` 或 `TYPE MD`，行为会改变，需要重新评估工作流（[第7章](07-geo-cell-opt.md) §8）。

**2. `unknown subsection PDOS of section PRINT`（2026.2）或 `unknown subsection PDOS of section DOS`（≤ 2026.1）**
原因：DOS/PDOS 输入在 2026.2 重构，新旧写法互不兼容。
解决：按运行的版本选择写法，或用 §9.4 的 `@IF` 分支。

**3. `found an unknown keyword EVAL_ENERGY_FORCES in section REFTRAJ`**
原因：2026.2 删除了 `EVAL_ENERGY_FORCES`/`EVAL_FORCES`。
解决：改为 `EVAL ENERGY_FORCES`（或 `EVAL ENERGY`）。`EVAL` 自 2024.1 起就已存在，2024.1–2026.1 中旧关键字仍可用但会打印 deprecated 警告，因此可以先在旧版本上改好。

**4. `The requested basis set <...> for element <...> was not found in the basis set files ...`**
原因：基组名不在所列文件中（MOLOPT 与 MOLOPT_UZH 的命名不同）。
解决：在数据文件中 `grep "^元素 "` 查真实名字，或增加一行 `BASIS_SET_FILE_NAME`。

**5. 程序正常结束，但电子数不对**
原因：基组与赝势的 qN 不一致，或 `CHARGE` 写错；CP2K 不会因为名字中的 qN 不同而报错。
解决：统一 qN，核对 `Number of electrons`。

**6. `Use the LSD option for an odd number of electrons`**
原因：奇数电子未开自旋极化。
解决：加 `UKS` 与正确的 `MULTIPLICITY`；闭壳层体系出现此错时检查结构是否缺/多原子。

**7. `GEOMETRY wrong or EMAX_SPLINE too small`**
原因：原子重叠或距离过近，多源于结构转换时丢失晶胞、扩胞/切面产生重复原子、CIF 无序位点未处理。
解决：可视化检查并去重；核对 `&CELL` 与坐标是否匹配。

**8. `An invalid set of cell vectors was specified. The cell volume is too small`**
原因：`&CELL` 漏写、单位错误或矢量共面。
解决：检查 `ABC`/`A`/`B`/`C`，默认单位 Å。

**9. `OT not possible with kpoint calculations`**
原因：2026.2 及以前 OT 不支持 k 点。
解决：改用对角化 + 混合（+ 展宽）。

**10. `SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`**
原因：SCF 在 `MAX_SCF`（及外层循环）内没有达到 `EPS_SCF`。
解决：见 §4.4 与[第5章](05-scf-strategy.md) §10；不要靠 `IGNORE_CONVERGENCE_FAILURE` 跑 MD 或优化。

**11. `Periodic Hartree Fock calculation requested with the use of a truncated or shortrange potential. The cutoff radius ... is larger than half the minimal cell dimension ...`**
原因：HFX 截断半径超过最小晶胞尺寸的一半。
解决：减小 `CUTOFF_RADIUS` 或扩大超胞，并把 $R_c$ 作为收敛参数测试（[第9章](09-hybrid-admm.md)）。

**12. 输出中多次出现 `DBCSR| MPI: Number of processes 1`**
原因：MPI 错配，实际是多个串行进程。
解决：§6，统一 `mpirun` 与 CP2K 链接的 MPI；作业脚本中加载正确环境。

**13. `error while loading shared libraries: ...`**
原因：计算节点没有加载运行环境。
解决：作业脚本中 `source install/cp2k_env`/`module load`/`conda activate`。

**14. 杂化泛函任务被系统杀掉（OOM）或第一步极慢**
原因：ERI 存储超出内存，或 `MAX_MEMORY` 太小导致每步重算。
解决：§7.5；减少每节点 MPI 进程、增加线程，使用 ADMM，先用 PBE 收敛并 `SCREEN_ON_INITIAL_P T`。

**15. 续算后第一步 SCF 很慢，日志中有 `This file does not exist ... Calculation continues using ATOMIC GUESS`**
原因：波函数重启文件名不对，或 Γ 点/k 点文件类型不匹配（`.wfn` vs `.kp`）。
解决：显式写 `WFN_RESTART_FILE_NAME`，确认文件存在且类型与当前计算一致。

**16. 续算的 MD 能量或温度在接口处跳变**
原因：没有恢复速度或恒温器状态、换了参数、波函数初猜不同导致 SCF 收敛到别的态。
解决：检查 `&EXT_RESTART` 的开关；保持 `EPS_SCF` 等参数一致；必要时同时提供波函数重启文件。

**17. 程序在对角化阶段崩溃或结果异常（带 ELPA 的构建）**
原因：ELPA 与特定硬件/编译器组合的问题。
解决：`&GLOBAL/PREFERRED_DIAG_LIBRARY SCALAPACK` 对比；问题消失就固定使用 ScaLAPACK 并报告。

**18. `Illegal instruction`**
原因：`-march=native` 编译的程序跑在较老的 CPU 上。
解决：在计算节点编译，或 `--target-cpu=generic`。

---

## 练习

1. 按 §5 的流程修复 `broken_scf.inp`。记录每一轮的报错原文、你的修改、下一轮的结果，填写成一张表；最后用 `diff` 把你的版本与 `broken_scf_fixed.inp` 对比，解释所有差异。
2. 故意只修 A–D 而保留 `GTH-PBE-q1`（如果你的数据文件中有 Cu 的 q1 赝势），运行后用 `grep "Number of electrons"` 和费米能比较 q1 与 q11 的结果，体会“正常结束但结果错误”。
3. 在 `broken_scf_fixed.inp` 中把 `ADDED_MOS` 依次改为 2、6、12、20，比较 SCF 步数、总能量和费米能，找出结果不再随 `ADDED_MOS` 变化的最小值（k 点计算不会打印“smearing includes the last MO”警告，只能靠这种测试判断）。再去掉 `&KPOINTS`、改用 Γ 点 2×2×2 超胞重复一次，观察 `ADDED_MOS` 太小时是否出现这条警告。
4. 用 `job_mpi_omp_hybrid.sh` 在你的集群上对一个生产体系做 MPI×OpenMP 基准，至少 4 种布局；再选最快的布局做 1、2、4 节点的强扩展测试，计算并行效率 $E(N)$。
5. 运行 `si_cellopt_2026_2.inp` 时在第 3–4 步手动取消任务，然后用 `si_cellopt_restart.inp` 续算。比较续算前后 `OPT|` 中的步数编号、体积和压力是否连续；再把 `RESTART_COUNTERS` 改为 `F` 重做一次，观察差别。
6. 找一个你自己以前在 2025.x 上写的输入（或 Multiwfn 生成的输入），先用 `check_migration.sh` 扫描，再用 2026.2 的 `cp2k.psmp -c` 检查，把它迁移到 2026.2，并用 §9.4 的方法让它同时兼容旧版本。
7. 为你当前课题的一个主要结论整理一张 §10.1 风格的收敛表，并写出与至少一篇文献或另一个程序的对比。

---

## 参考资料

**官方文档**

- Troubleshooting（官方报错目录：输入/输出去向、程序卡住、MPI 错配、星号/NaN、CPASSERT、`GEOMETRY wrong`、LSD 消息）：<https://manual.cp2k.org/trunk/getting-started/troubleshooting.html>
- Restarting CP2K Calculations（运动重启、`RESTART_HISTORY`、`&EXT_RESTART`、`.wfn`/`.kp` 兼容规则）：<https://manual.cp2k.org/trunk/methods/restarting.html>
- How to make a SCF run converge：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- How to Converge the CUTOFF and REL_CUTOFF：<https://manual.cp2k.org/trunk/methods/dft/cutoff.html>
- Hartree-Fock exchange / ADMM：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/index.html>
- Density of States（2026.2 新 `&DOS` 接口）：<https://manual.cp2k.org/trunk/methods/electronic_structure/dos.html>
- Build from Source（CMake、toolchain + build_cp2k.sh、Spack）：<https://manual.cp2k.org/trunk/getting-started/build-from-source.html>
- Build with Spack：<https://manual.cp2k.org/trunk/getting-started/build-with-spack.html>
- 输入参考（2026.2）：`&GLOBAL`（`PREFERRED_DIAG_LIBRARY`、`WALLTIME`、`&TIMINGS`）<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/GLOBAL.html>；`&EXT_RESTART` <https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/EXT_RESTART.html>
- CP2K 2026.2 Manual：<https://manual.cp2k.org/cp2k-2026_2-branch/index.html>；各版本手册：<https://manual.cp2k.org/trunk/versions.html>
- 发布说明：v2026.2 <https://github.com/cp2k/cp2k/releases/tag/v2026.2>；v2026.1 <https://github.com/cp2k/cp2k/releases/tag/v2026.1>
- CP2K Google Group（用户论坛，提问前先读官方 SUPPORT 说明）：<https://groups.google.com/g/cp2k>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：源码与 `tests/` 回归测试输入，查关键字真实写法的最权威来源；`docs/` 是手册源文件
- [cp2k/cp2k-containers](https://github.com/cp2k/cp2k-containers)：官方容器（含 CUDA 镜像），可快速得到一个“已知正确”的环境来排除本地编译问题
- [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools)：`cp2klint` 输入检查（内置结构较旧，最终以 `cp2k.psmp -c` 为准）
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：`cp2kparse` 把输出转成 JSON/YAML；`xyz_restart_parser` 去除 MD 续算产生的重复帧
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例，含 `converging_grid` 收敛测试范例
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：中文 Wiki，CP2K 实践经验
- [szt2006/cp2k-aimd-skill](https://github.com/szt2006/cp2k-aimd-skill)：中文“症状到处方”对照表（license 不明，仅链接）
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：restart 回写、OT/DIAG 切换、自动识别 PBS/Slurm 的提交脚本思路（无 license，仅参考）
- [EPCCed/intro_to_cp2k](https://github.com/EPCCed/intro_to_cp2k)：ARCHER2 入门课，含并行性能练习
- [hfp/xconfigure](https://github.com/hfp/xconfigure)：Intel 编译器下编译 CP2K 的历史配方（Makefile 时代，仅作参考）
- [DCM-Uni-Paderborn/UZH-protocol](https://github.com/DCM-Uni-Paderborn/UZH-protocol)：UZH 基组与赝势的验证数据

**关键文献**

- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>
- J. VandeVondele *et al.*, *Comput. Phys. Commun.* **167**, 103 (2005)，Quickstep：<https://doi.org/10.1016/j.cpc.2004.12.014>
- The CP2K Program Package Made Simple, *J. Phys. Chem. B* (2025)：<https://doi.org/10.1021/acs.jpcb.5c05851>
- U. Borštnik *et al.*, *Parallel Computing* **40**, 47 (2014)，DBCSR 稀疏矩阵库：<https://doi.org/10.1016/j.parco.2014.03.012>
- A. Marek *et al.*, *J. Phys.: Condens. Matter* **26**, 213201 (2014)，ELPA 本征值求解器：<https://doi.org/10.1088/0953-8984/26/21/213201>
- M. Guidon, J. Hutter, J. VandeVondele, *J. Chem. Theory Comput.* **5**, 3010 (2009)，周期体系 HFX 与截断库仑：<https://doi.org/10.1021/ct900494g>
- M. Guidon, J. Hutter, J. VandeVondele, *J. Chem. Theory Comput.* **6**, 2348 (2010)，ADMM：<https://doi.org/10.1021/ct1002225>
- K. Lejaeghere *et al.*, *Science* **351**, aad3000 (2016)，DFT 程序间的可重复性（Δ 测试）：<https://doi.org/10.1126/science.aad3000>
- E. Bosoni *et al.*, *Nat. Rev. Phys.* **6**, 45 (2024)，DFT 实现精度的验证工作流：<https://doi.org/10.1038/s42254-023-00655-3>
- P. Haas, F. Tran, P. Blaha, *Phys. Rev. B* **79**, 085104 (2009)，多种泛函下固体晶格常数的系统比较：<https://doi.org/10.1103/PhysRevB.79.085104>

---

[← 上一章](17-tools-ecosystem.md) | [返回目录](../README.md)
