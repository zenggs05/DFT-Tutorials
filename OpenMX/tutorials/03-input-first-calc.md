# 第3章 输入文件结构与第一个计算（甲烷分子）（OpenMX Input File and First Calculation: Methane）

本章以甲烷分子（CH4）为例，逐行讲解 OpenMX 输入文件 `.dat` 的语法与核心关键词：`System.CurrrentDirectory`、`System.Name`、`DATA.PATH`、`Definition.of.Atomic.Species`（PAO 基组 + VPS 赝势）、`Atoms.SpeciesAndCoordinates`（含上/下自旋初始电子数）、`Atoms.UnitVectors` 与分子计算用的 `scf.EigenvalueSolver Cluster`。随后运行第一个 OpenMX SCF 计算（OpenMX first calculation, methane），读懂标准输出中每一步的 NormRD、dUele，读懂 `.out` 中的总能量分解、本征值、Mulliken 电荷、偶极矩与原子力，并把同样的方法推广到 H2O 和 CO 分子。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出一个语法正确的 OpenMX 输入文件，知道关键词、块关键词（`<Name ... Name>`）、注释和大小写规则；
- 正确设置 `DATA.PATH`，让 OpenMX 在任意目录下都能找到 PAO 与 VPS 文件；
- 根据 VPS 的价电子数填写坐标行的上/下自旋初始电子数；
- 为分子选择足够大的盒子，并用 `Cluster` 求解器计算；
- 运行 `mpirun -np 4 openmx Methane.dat -nt 1 | tee met.std`，判断计算是否正常结束、SCF 是否收敛；
- 从 `.out` 中提取总能量、HOMO/LUMO、Mulliken 电荷、偶极矩、原子力，并换算成常用单位。

本章示例：[03-first-calc-methane](../examples/03-first-calc-methane/README.md)，包含 `Methane.dat`、`Methane_manual_check.dat`、`H2O.dat`、`CO.dat`、Slurm 作业脚本与结果检查脚本。

> [!NOTE]
> 本章默认你已按 [第2章](02-install.md) 安装好 OpenMX 4.0（源码编译、deb 包或 Docker 均可），并且 `openmx` 在 `PATH` 中或你知道它的完整路径。

---

## 2. 原理速览：一次 OpenMX SCF 计算做了什么

### 2.1 Kohn–Sham 方程与 SCF 循环

OpenMX 用 PAO 基组 $\{\chi_\mu\}$ 展开 Kohn–Sham 轨道，求解广义本征值问题 $H\mathbf c=\varepsilon S\mathbf c$。哈密顿量依赖于电子密度，因此需要自洽迭代：

```text
initial density (superposition of atomic densities, spin set by columns 6-7 of the coordinate lines)
   └─► build H[rho] on the real-space grid (Hartree by FFT, XC on grid points)
        └─► solve H c = e S c  (Cluster: Gamma only; Band: k points)
             └─► occupy states with Fermi-Dirac at scf.ElectronicTemperature -> new density matrix
                  └─► mix old/new (scf.Mixing.Type) -> converged?  no -> back to "build H"
                                                         yes -> total energy, forces, Mulliken, output files
```

每一步 OpenMX 都会打印两个最重要的量：

- **Uele**：能带能量（占据态本征值之和，单位 Hartree）；**dUele** 是相邻两步 Uele 之差的绝对值。SCF 的收敛条件是 `dUele < scf.criterion`；
- **NormRD**：混合方案中残差（输入与输出之差）的范数，越接近 0 说明越自洽。手册只写了 dUele 判据；从源码（`DFT.c`）看，判定收敛时还会同时检查混合残差是否已经足够小，所以只是 dUele 偶然很小、NormRD 仍很大时不会被误判为收敛。

### 2.2 分子计算为什么也要“盒子”

OpenMX 用 FFT 解 Poisson 方程，所以即使是孤立分子，也是放在周期性盒子里计算的。对分子：

- 用 `scf.EigenvalueSolver Cluster`，只考虑 Γ 点，不需要 k 点；
- 盒子要足够大，使中心盒子里的基函数与相邻镜像中的基函数**不重叠**。粗略判据：盒长 $L \gtrsim$ 分子在该方向上的尺寸 $+\,2r_c$，其中 $r_c$ 是所用 PAO 的最大截断半径（6.0 Bohr ≈ 3.18 Å）；
- 满足上述条件后，镜像之间仍有经典静电相互作用。对非极性分子影响很小；对极性分子可以适当加大盒子；对带电分子，可用 `scf.coulomb.cutoff on`（精确库仑截断法，手册第 61 节，截断半径为最短晶格矢量长度的一半，所以盒子要足够大）。

### 2.3 总能量的组成

`.out` 中会给出总能量 `Utot` 及其分解（手册第 5 节）：

$$
U_{\rm tot}=U_{\rm kin}+U_{\rm H0}+U_{\rm H1}+U_{\rm na}+U_{\rm nl}+U_{\rm xc0}+U_{\rm xc1}+U_{\rm core}+U_{\rm hub}+U_{\rm cs}+U_{\rm zs}+U_{\rm zo}+U_{\rm ef}+U_{\rm vdW}
$$

各项含义见 §7.2。实际研究中通常只用 `Utot`（以及优化时的力），能量分解主要用于检查与理解。

---

## 3. 输入文件 `.dat` 的语法规则

OpenMX 的输入是一种“自由格式”文本：

| 规则 | 说明 | 例子 |
|---|---|---|
| 普通关键词 | 一行：`关键词  值`（值可以有多个，用空格分隔） | `scf.energycutoff   200.0` |
| 块关键词 | 以 `<Name` 开始、以 `Name>` 结束，中间每行一条记录 | `<Atoms.UnitVectors` ... `Atoms.UnitVectors>` |
| 顺序 | 关键词顺序任意，可有空行 | — |
| 大小写 | 关键词名和选项值**不区分大小写** | `scf.eigenvaluesolver cluster` 与 `scf.EigenvalueSolver Cluster` 等价 |
| 注释 | `#` 之后的内容是注释；值后面也可以跟注释 | `scf.maxIter  100   # default=40` |
| 未知关键词 | **被静默忽略**，不报错 | 拼错的关键词会悄悄使用默认值 |
| 重复关键词 | 4.0 的解析器按行查找关键词，遇到第一处即返回，所以**以第一次出现的为准** | 修改参数时不要只在文件末尾追加一行 |
| 非法的值 | 打印 `ERROR key=...`，读完输入后提示 `errors in the inputfile` 并退出 | `scf.XcType PBE`（应为 `GGA-PBE`） |

> [!WARNING]
> **`System.CurrrentDirectory` 的官方拼写是三个 r**，这是 OpenMX 的历史写法。写成 `System.CurrentDirectory` 时，官方版本会把它当作未知关键词忽略，输出仍写到默认的 `./`。（个别第三方分支额外接受两个 r 的写法，不要依赖。）

> [!IMPORTANT]
> 输入文件必须是**纯 ASCII（或不带 BOM 的 UTF-8）+ LF 换行**。在 Windows 记事本、Word 中编辑后常见的问题有：UTF-16 编码、BOM、全角空格、CRLF 换行，它们会导致物种名无法识别（`Found an undefined species name ...`）等莫名其妙的错误。检查与修复：
>
> ```bash
> file Methane.dat                              # should say: ASCII text
> grep -nP '[^\x00-\x7F]' Methane.dat           # list lines containing non-ASCII characters
> iconv -f UTF-16 -t UTF-8 bad.dat > good.dat   # if the file is UTF-16
> dos2unix good.dat                             # CRLF -> LF
> ```

---

## 4. 逐段讲解：甲烷输入文件

下面是示例目录中的 `Methane.dat`（GGA-PBE、Standard 基组、200 Ry），我们按段落讲解。完整文件：

```text
#
# Methane (CH4) single-point SCF with OpenMX 4.0
# DFT-Tutorials, OpenMX chapter 3
#

#
# File names and output control
#
System.CurrrentDirectory      ./          # official spelling: three 'r'
System.Name                   met
DATA.PATH                     /path/to/openmx4.0/DFT_DATA19   # EDIT: absolute path that contains PAO/ and VPS/
level.of.stdout               1           # -1..2, default=1
level.of.fileout              1           # 0..2,  default=1 (1 writes cube files)

#
# Definition of atomic species
#
Species.Number                2
<Definition.of.Atomic.Species
  H    H6.0-s2p1      H_PBE19
  C    C6.0-s2p2d1    C_PBE19
Definition.of.Atomic.Species>

#
# Atoms (tetrahedral CH4, C-H = 1.089 Ang, centred in a 10 Ang cubic box)
#
Atoms.Number                  5
Atoms.SpeciesAndCoordinates.Unit   Ang    # Ang|AU|FRAC
<Atoms.SpeciesAndCoordinates
  1   C    5.000000   5.000000   5.000000    2.0  2.0
  2   H    5.628734   5.628734   5.628734    0.5  0.5
  3   H    4.371266   4.371266   5.628734    0.5  0.5
  4   H    4.371266   5.628734   4.371266    0.5  0.5
  5   H    5.628734   4.371266   4.371266    0.5  0.5
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit        Ang         # Ang|AU
<Atoms.UnitVectors
  10.0   0.0   0.0
   0.0  10.0   0.0
   0.0   0.0  10.0
Atoms.UnitVectors>

#
# SCF or electronic system
#
scf.XcType                    GGA-PBE     # LDA|LSDA-CA|LSDA-PW|GGA-PBE
scf.SpinPolarization          off         # On|Off|NC
scf.ElectronicTemperature     300.0       # K, default=300
scf.energycutoff              200.0       # Ry; always set it explicitly
scf.maxIter                   100         # default=40
scf.EigenvalueSolver          Cluster     # Cluster for isolated molecules
scf.Kgrid                     1 1 1       # not used by the Cluster solver
scf.Mixing.Type               RMM-DIIS    # Simple|RMM-DIIS|GR-Pulay|Kerker|RMM-DIISK|RMM-DIISV|RMM-DIISH
scf.Init.Mixing.Weight        0.30        # default=0.30
scf.Min.Mixing.Weight         0.001       # default=0.001
scf.Max.Mixing.Weight         0.40        # default=0.40
scf.Mixing.History            7           # default=5
scf.Mixing.StartPulay         5           # default=6
scf.criterion                 1.0e-8      # Hartree, default=1.0e-6

#
# MD or geometry optimization
#
MD.Type                       Nomd        # single point
MD.maxIter                    1
```

### 4.1 文件名与输出控制

| 关键词 | 本例取值 | 默认值 | 说明 |
|---|---|---|---|
| `System.CurrrentDirectory` | `./` | `./` | 输出文件写到哪个目录（三个 r！）。改成子目录时建议事先 `mkdir` |
| `System.Name` | `met` | — | 所有输出文件的前缀：`met.out`、`met.xyz`、`met_rst/` ... 同一目录中不同计算要用不同的名字，否则会互相覆盖 |
| `DATA.PATH` | 绝对路径 | `../DFT_DATA19` | 指向**同时包含 `PAO/` 和 `VPS/`** 的目录。OpenMX 会在 `DATA.PATH/PAO/` 中找 `H6.0.pao`、在 `DATA.PATH/VPS/` 中找 `H_PBE19.vps` |
| `level.of.stdout` | 1 | 1 | 标准输出详细程度：−1 不输出，0 最少（每个 SCF 步一行），1 标准，2 面向开发者 |
| `level.of.fileout` | 1 | 1 | 文件输出量：0 最少（不写 cube 与网格文件），1 标准（写电荷密度、势的 cube），2 额外写交换关联势 cube 与网格文件 |

> [!TIP]
> 运行示例前，先把 `DATA.PATH` 改成你自己的路径，例如：
>
> ```bash
> sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/software/openmx4.0/DFT_DATA19#" *.dat
> grep -H "^DATA.PATH" *.dat
> ```
>
> deb 包用户通常是 `/opt/openmx4.0/DFT_DATA19`（用 `ls /opt/openmx4.0` 确认）。

### 4.2 原子种类：基组与赝势

```text
Species.Number                2
<Definition.of.Atomic.Species
  H    H6.0-s2p1      H_PBE19
  C    C6.0-s2p2d1    C_PBE19
Definition.of.Atomic.Species>
```

| 列 | 含义 | 本例 |
|---|---|---|
| 第 1 列 | 物种名（自己取，坐标块中引用它） | `H`、`C` |
| 第 2 列 | PAO 文件名（不带 `.pao`）+ `-` + 每个角动量的径向函数个数 | `H6.0-s2p1`：`H6.0.pao`（$r_c$ = 6.0 Bohr），2 个 s、1 个 p |
| 第 3 列 | VPS 文件名（不带 `.vps`） | `H_PBE19`：氢的 GGA-PBE 赝势（DFT_DATA19） |

`Species.Number` 必须等于块中的行数。基函数个数：

| 物种 | 基组 | 基函数数 |
|---|---|---|
| H | `H6.0-s2p1` | $2\times1+1\times3=5$ |
| C | `C6.0-s2p2d1` | $2\times1+2\times3+1\times5=13$ |
| CH4 合计 | — | $13+4\times5=33$ |

这里选的是 DFT_DATA19 推荐表中的 **Standard** 档（H 的 Quick/Standard/Precise 分别为 `H5.0-s2`、`H6.0-s2p1`、`H7.0-s2p2d1`；C 为 `C6.0-s2p2`、`C6.0-s2p2d1`、`C6.0-s3p2d2`）。怎样选基组、什么是收缩写法 `s2>1`、同一元素定义多个物种等，见 [第4章](04-pao-vps.md)。

> [!IMPORTANT]
> 第 3 列 VPS 的泛函必须与 `scf.XcType` 一致：GGA-PBE 计算用 `*_PBE19`，LDA 计算用 `*_CA19`。PAO 与 VPS 都应来自同一数据库版本（DFT_DATA19），不要与旧版 DFT_DATA13 的文件混用。

### 4.3 原子坐标与初始电子数

```text
Atoms.Number                  5
Atoms.SpeciesAndCoordinates.Unit   Ang    # Ang|AU|FRAC
<Atoms.SpeciesAndCoordinates
  1   C    5.000000   5.000000   5.000000    2.0  2.0
  2   H    5.628734   5.628734   5.628734    0.5  0.5
  ...
Atoms.SpeciesAndCoordinates>
```

| 列 | 含义 |
|---|---|
| 1 | 原子序号，必须从 1 开始连续编号（否则报 `Format error of the sequential number ...`） |
| 2 | 物种名，必须在 `Definition.of.Atomic.Species` 中定义过 |
| 3–5 | x、y、z 坐标，单位由 `Atoms.SpeciesAndCoordinates.Unit` 决定：`Ang`（Å）、`AU`（Bohr）或 `FRAC`（以 a、b、c 为基的分数坐标，0–1 之外的值读入后会自动折回） |
| 6 | 该原子的**初始上自旋电子数** |
| 7 | 该原子的**初始下自旋电子数** |

第 6、7 列之和应等于该原子所用 VPS 的价电子数，可以用 `grep valence.electron` 查到：

```bash
grep valence.electron $HOME/software/openmx4.0/DFT_DATA19/VPS/{H,C,O}_PBE19.vps
```

| VPS | 价电子数 | 非磁性初值（上 下） |
|---|---|---|
| `H_PBE19` | 1 | `0.5 0.5` |
| `C_PBE19` | 4 | `2.0 2.0` |
| `O_PBE19` | 6 | `3.0 3.0` |

体系的总电子数由各原子 VPS 的价电子数之和（再减去 `scf.system.charge`）决定：CH4 为 $4+4\times1=8$，占据 4 个轨道，所以 `.out` 中会出现 `HOMO = 4`。第 6、7 列只决定**初始**的原子密度和自旋分布，填错不一定会报错，但会让初始密度偏离中性、SCF 起步变差，自旋极化计算还可能收敛到错误的磁态。磁性体系如何用这两列设置初始磁矩（以及非共线计算的第 8 列以后），见 [第11章](11-spin-noncollinear.md)。

> [!TIP]
> 本例把 C 放在盒子中心 (5, 5, 5) Å，只是为了让 cube 文件在可视化时分子位于盒子中间；放在原点附近也完全可以（官方 `work/Methane.dat` 就把 C 放在原点）。周期性盒子中，坐标平移不改变能量。

### 4.4 晶胞（盒子）

```text
Atoms.UnitVectors.Unit        Ang
<Atoms.UnitVectors
  10.0   0.0   0.0
   0.0  10.0   0.0
   0.0   0.0  10.0
Atoms.UnitVectors>
```

三行分别是晶格矢量 **a**、**b**、**c**。检查盒子是否够大：CH4 中 H 原子在每个方向上的跨度约 1.26 Å，加上两倍的 $r_c$（2 × 3.18 Å），约 7.6 Å < 10 Å，满足基函数不重叠的条件。

对 Cluster 计算，如果**省略** `<Atoms.UnitVectors` 块，OpenMX 会根据原子位置和 PAO 截断半径自动确定一个最小盒子，并在标准输出中以 `<Set_Cluster_UnitCell>` 开头的行打印出来（手册第 69 节），你可以参考它来设置自己的盒子。但为了结果可重复、便于比较能量，建议总是显式给出盒子。

### 4.5 SCF 参数

| 关键词 | 本例 | 默认 | 说明与选择理由 |
|---|---|---|---|
| `scf.XcType` | `GGA-PBE` | — | 4.0 可选 `LDA`、`LSDA-CA`、`LSDA-PW`、`GGA-PBE` |
| `scf.SpinPolarization` | `off` | `off` | 闭壳层分子；自旋极化用 `on`，非共线用 `nc`。`LDA` 只能配 `off` |
| `scf.ElectronicTemperature` | 300 K | 300 K | 费米–狄拉克展宽的电子温度。对有能隙的分子影响可忽略 |
| `scf.energycutoff` | 200 Ry | 4.0 手册：220 Ry | 实空间积分网格的截断能，不是基组截断能；显式设置，并做收敛测试（[第5章](05-scf-kpoints-convergence.md)） |
| `scf.maxIter` | 100 | 40 | SCF 最大步数，到达后即使未收敛也会停止 |
| `scf.EigenvalueSolver` | `Cluster` | — | 分子用 `Cluster`；周期体系用 `Band`；大体系 O(N) 用 `DC`、`DC-LNO`、`Krylov` 等 |
| `scf.Kgrid` | `1 1 1` | — | Cluster 计算不使用；保留它便于以后改成 Band |
| `scf.Mixing.Type` | `RMM-DIIS` | — | 小分子用 RMM-DIIS 很快；周期体系、金属、磁性体系通常用 `RMM-DIISK` |
| `scf.Init/Min/Max.Mixing.Weight` | 0.30 / 0.001 / 0.40 | 同左 | 初始、最小、最大混合权重 |
| `scf.Mixing.History` | 7 | 5 | Pulay 型混合参考的历史步数；难收敛体系手册建议约 30 |
| `scf.Mixing.StartPulay` | 5 | 6 | 从第几步开始 Pulay 型混合；难收敛体系建议 10–30 |
| `scf.criterion` | 1.0e-8 | 1.0e-6 | dUele 的收敛阈值（Hartree）。单点能与性质计算 1e-8 足够，结构优化可更严 |

### 4.6 MD / 结构优化部分

`MD.Type Nomd` 表示只做一次 SCF（单点能），`MD.maxIter 1`。结构优化（`MD.Type RF`、`EF`、`BFGS` 等）和分子动力学见 [第6章](06-optimization-eos-neb.md)、[第7章](07-md.md)。

### 4.7 官方 `work/Methane.dat` 有何不同

官方示例为了让安装测试只需几秒，用的是最小基组 `H5.0-s1`、`C5.0-s1p1`（CH4 共 8 个基函数）和 120 Ry 截断能，C 位于原点。示例目录中的 `Methane_manual_check.dat` 按相同的物理设置编写（基组、赝势、坐标、盒子、泛函、截断能、电子温度一致），运行后得到的 `Utot` 可以直接与手册第 5 节的参考值 −8.033540953061 Hartree 比较，用来确认你的安装和数据库没有问题。最小基组只适合做测试，不要用于研究。

---

## 5. 运行计算

### 5.1 准备工作目录

```bash
mkdir -p ~/omx_runs/ch03 && cd ~/omx_runs/ch03
cp /path/to/DFT-Tutorials/OpenMX/examples/03-first-calc-methane/*.dat .
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/software/openmx4.0/DFT_DATA19#" *.dat
```

可选：用第1章示例中的检查工具确认物种定义、基函数数和初始电子数是否自洽（用法见 `python list_dft_data.py -h`）：

```bash
python /path/to/DFT-Tutorials/OpenMX/examples/01-overview-tour/list_dft_data.py check-dat Methane.dat
```

### 5.2 运行命令

```bash
export OMP_NUM_THREADS=1
mpirun -np 4 openmx Methane.dat -nt 1 | tee met.std
```

| 部分 | 含义 |
|---|---|
| `mpirun -np 4` | 启动 4 个 MPI 进程 |
| `openmx Methane.dat` | 主程序 + 输入文件（若 `openmx` 不在 `PATH` 中，写完整路径） |
| `-nt 1` | 每个 MPI 进程的 OpenMP 线程数；不写时为 1 |
| `tee met.std`（经管道） | 标准输出同时显示在屏幕上并保存到 `met.std` |

几点说明：

- 对 CH4 这样的小分子，Cluster 计算的 MPI 并行主要分配在本征态上，用 1–4 个进程即可，进程再多反而更慢；
- 后台运行：`nohup mpirun -np 4 openmx Methane.dat -nt 1 > met.std 2>&1 &`，用 `tail -f met.std` 查看进度；
- 在集群上请通过作业系统提交，示例目录中的 `job_omx_slurm.sh` 会依次运行三个分子；通用模板见 [`job_omx_slurm.sh`](../scripts/job_omx_slurm.sh) 与 [第18章](18-parallel-on-errors.md)。

### 5.3 生成的文件

计算正常结束后（`level.of.fileout 1`），当前目录中会出现：

| 文件 | 内容 |
|---|---|
| `met.std` | 标准输出（`tee` 保存）：SCF 过程 |
| `met.out` | 主输出：输入文件副本、SCF 历史、总能量及分解、本征值、Mulliken 电荷、偶极矩、坐标与力、计算时间 |
| `met.xyz` | 最终结构（xyz 格式，可用 OpenMX Viewer、XCrySDen 打开） |
| `met.md` | 每个 MD/优化步的结构（多帧 xyz），单点计算只有一帧 |
| `met.md2` | 最后一步的结构（带你定义的物种名） |
| `met.ene` | 每个 MD 步的各种数值 |
| `met.cif` | 初始结构的 CIF |
| `met.tden.cube` | 总电子密度 |
| `met.dden.cube` | 差分电子密度（相对孤立原子密度叠加） |
| `met.v0.cube` | Kohn–Sham 势（上自旋，不含非局域部分） |
| `met.vhart.cube` | 差分电荷产生的 Hartree 势 |
| `met_rst/` | 重启文件目录（`scf.restart on`/`pp` 时读取） |

自旋极化计算还会有 `met.sden.cube`（自旋密度）和 `met.v1.cube`（下自旋势）；`level.of.fileout 2` 还会写 `met.vxc0.cube`、`met.vxc1.cube` 和网格文件；`level.of.fileout 0` 则不写任何 cube。

---

## 6. 读懂标准输出 `.std`

### 6.1 读入阶段

标准输出开头是欢迎横幅（含版本号），随后是读入输入文件、读取 PAO/VPS、计算重叠矩阵等信息。读取数据库时会出现类似下面的行（示意输出，具体内容以实际为准）：

```text
<SetPara_DFT>  PAOs of species H were normally found.
<SetPara_DFT>  PAOs of species C were normally found.
<SetPara_DFT>  VPSs of species H were normally found.
               H_PBE19.vps is j-dependent.
               In case of scf.SpinOrbit.Coupling=off,
               j-dependent pseudo potentials are averaged by j-degeneracy,
               which corresponds to a scalar relativistic treatment.
```

最后几行说明：DFT_DATA19 的 VPS 是全相对论（$j$ 依赖）的；不开自旋轨道耦合时，程序会按 $j$ 简并度平均，相当于**标量相对论**处理。这是正常信息，不是错误。

### 6.2 每个 SCF 步

`level.of.stdout 1` 时，每个 SCF 步都以 `MD= 1  SCF= n` 的分隔行开头，以能量和残差结束（示意输出，数值仅供参考）：

```text
******************* MD= 1  SCF= 7 *******************
  ...
<DFT>  Total Spin Moment (muB) =  0.000000000000
<DFT>  Mixing_weight= 0.300000000000
<DFT>  Uele   =     -3.889312346886  dUele     =      0.002941147152
<DFT>  NormRD =      0.000838640097  Criterion =      0.000000010000
```

| 量 | 含义 | 看什么 |
|---|---|---|
| `Total Spin Moment` | 总自旋磁矩（μB） | 闭壳层分子应为 0 |
| `Mixing_weight` | 当前混合权重 | 自适应变化，过小说明在“慢慢挪” |
| `Uele` | 能带能量（Hartree） | 逐步稳定 |
| `dUele` | 与上一步 Uele 之差 | 收敛时小于 `scf.criterion` |
| `NormRD` | 混合残差范数 | 应在震荡后快速下降几个数量级 |
| `Criterion` | 即 `scf.criterion` | — |

`level.of.stdout 0` 时，每步只打印一行 `<DFT> MD=   1 SCF=   7   NormRD = ...  Criterion = ...`，适合大体系长时间计算。

常用检查命令：

```bash
grep "dUele" met.std | tail -n 5                 # last few SCF steps
grep "NormRD" met.std | awk '{print $4}' | tail  # NormRD history
grep -c "SCF=" met.std                           # rough number of SCF steps printed
tail -n 2 met.std                                # should end with: The calculation was normally finished.
```

### 6.3 如何判断 SCF 收敛了

1. `met.std` 最后出现 `The calculation was normally finished.`——这只说明程序正常结束，**不等于 SCF 收敛**；
2. 查看 `met.out` 的 SCF 历史：最后一步的序号应**小于** `scf.maxIter`，最后几步的 `NormRD` 应比第 1 步（1.0）下降许多个数量级，`Uele` 在 `scf.criterion` 的精度上不再变化；
3. 如果 SCF 步数正好等于 `scf.maxIter`，说明是被步数上限截断的，结果不可靠，需要增大 `scf.maxIter` 或调整混合参数（[第5章](05-scf-kpoints-convergence.md)）。

```bash
grep -E "^ +SCF= +[0-9]+ +NormRD" met.out | tail -n 3
grep -E "^ +SCF= +[0-9]+ +NormRD" met.out | wc -l     # number of SCF steps in the MD=1 history
grep -i "scf.maxIter" Methane.dat
```

---

## 7. 读懂主输出 `.out`

`met.out` 依次包含：输入文件副本 → SCF 历史 → 总能量及分解 → 化学势 → 本征值 → Mulliken 电荷 → 偶极矩 → 坐标与力 → 分数坐标 → 计算时间。

### 7.1 SCF 历史

```text
***********************************************************
                  SCF history at MD= 1
***********************************************************
   SCF=   1  NormRD=  1.000000000000  Uele= -3.523169099737
   SCF=   2  NormRD=  0.181517404414  Uele= -3.686855123752
   ...
```

（以上数值摘自手册第 5 节官方 `Methane.dat` 的参考输出，该算例 16 步收敛。）

### 7.2 总能量及其分解

```bash
sed -n '/Total energy (Hartree) at MD = 1/,/Chemical potential/p' met.out
grep "Utot\." met.out
```

手册第 5 节给出的官方 `Methane.dat`（最小基组、120 Ry）参考值节选如下：

```text
  Uele.         -3.889316014119
  Ukin.          5.533759381311
  ...
  Utot.         -8.033540953061
```

| 项 | 含义 |
|---|---|
| `Uele` | 能带能量（占据本征值之和），不是总能量 |
| `Ukin` | 动能 |
| `UH0`、`UH1` | 屏蔽库仑能的静电部分、差分电子–电子库仑能 |
| `Una` | 中性原子势能量 |
| `Unl` | 非局域赝势能量 |
| `Uxc0`、`Uxc1` | 上、下自旋的交换关联能 |
| `Ucore` | 芯–芯库仑能 |
| `Uhub` | DFT+U 能量 |
| `Ucs` | 自旋方向约束能 |
| `Uzs`、`Uzo` | 自旋、轨道磁矩的 Zeeman 能 |
| `Uef` | 外电场能量 |
| `UvdW` | 半经验色散（DFT-D2/D3）能量 |
| `Utot` | 总能量 = 上述各项（除 Uele）之和，**做能量比较时用它** |
| `UpV`、`Enpy` | 压强项与焓（变胞计算时有意义；单点计算 `Enpy` = `Utot`） |

换算成 eV（1 Hartree = 27.211386 eV）：

```bash
grep "Utot\." met.out | awk '{printf "Utot = %.10f Ha = %.6f eV\n", $2, $2*27.211386}'
```

> [!NOTE]
> 赝势计算的总能量绝对值没有直接物理意义，只有在**相同的** VPS、基组、截断能、盒子/k 点设置下的能量**差**才有意义（反应能、结合能、相对稳定性等）。不同程序（OpenMX 与 VASP）之间的总能量绝对值不能比较。

### 7.3 化学势与本征值

```bash
grep -A 16 "Eigenvalues (Hartree) for SCF KS-eq." met.out
```

```text
（手册第 5 节官方参考输出节选）
   Chemical Potential (Hartree) =   0.00000000000000
   Number of States             =   8.00000000000000
   HOMO =  4
   Eigenvalues
                Up-spin            Down-spin
          1  -0.69897506409419  -0.69897506409419
          2  -0.41523055770122  -0.41523055770122
          3  -0.41523055765756  -0.41523055765756
          4  -0.41522182760629  -0.41522182760629
          5   0.21221759601986   0.21221759601986
          ...
```

- `HOMO = 4`：8 个价电子占据 4 个轨道；第 5 个是 LUMO；
- 本征值单位是 Hartree；第 2–4 个态是三重简并的 $t_2$ 轨道（数值在 1e-5 Ha 内相同，微小差别来自实空间网格破坏了严格的四面体对称性）；
- 对有能隙的分子，`Chemical Potential` 只是落在 HOMO 与 LUMO 之间的某个位置（由电子数守恒确定，不一定在能隙正中），判断能隙时请直接用 HOMO 与 LUMO 本征值；
- 官方最小基组算例中 HOMO–LUMO 差约 0.627 Ha（约 17 eV），这么大的数值主要是最小基组造成的，不具物理意义——这也说明了**测试用的基组不能用于研究**。

### 7.4 Mulliken 电荷

```bash
sed -n '/Mulliken populations/,/Decomposed Mulliken/p' met.out
```

```text
（手册第 5 节官方参考输出节选）
                    Up spin      Down spin     Sum           Diff
      1    C      2.509748760  2.509748760   5.019497520   0.000000000
      2    H      0.372562810  0.372562810   0.745125620   0.000000000
      ...
 Sum of MulP: up   =     4.00000 down          =     4.00000
              total=     8.00000 ideal(neutral)=     8.00000
```

- 每个原子的上/下自旋 Mulliken 布居、总和（Sum）与差值（Diff，即局域自旋磁矩，单位 μB）；
- `Sum of MulP: total` 必须等于 `ideal(neutral)`（中性体系的价电子总数，带电体系会相应改变），这是检查电子数设置的好方法；
- 原子电荷 = VPS 价电子数 − Mulliken 布居。例如 C：4 − 5.019 = −1.019 e，H：1 − 0.745 = +0.255 e；
- 其后的 `Decomposed Mulliken populations` 给出按轨道（s、px、py、pz、d…，以及多个径向函数 `multiple` 0、1…）分解的布居，在分析 PDOS 时用于确认轨道顺序（[第9章](09-dos-pdos-cohp-optics.md)）；
- Mulliken 电荷强烈依赖基组，只适合做同一基组下的定性比较；更多电荷分析方法见 [第10章](10-charge-mo-cube.md)。

### 7.5 偶极矩

```bash
grep -A 8 "Dipole moment (Debye)" met.out
```

输出 `Absolute D`（偶极矩大小）及 Total、Core、Electron、Back ground 的 x/y/z 分量，单位 Debye。CH4 有四面体对称性，偶极矩应接近 0；H2O 和 CO 则不为 0（见 §9）。

### 7.6 坐标与原子力

```bash
awk '/<coordinates.forces/,/coordinates.forces>/' met.out
```

块中每行依次是：序号、物种、x、y、z（Å）、Fx、Fy、Fz（Hartree/Bohr）。用途：

- 判断结构是否接近平衡：力的最大值与 `MD.Opt.criterion`（默认 3e-4 Ha/Bohr ≈ 0.015 eV/Å）比较；
- 检查对称性：CH4 四个 H 上的力大小应相同、方向沿 C–H 键。

> [!TIP]
> 官方最小基组算例中 H 原子上的力分量达到 0.06 Ha/Bohr 量级（手册参考输出），远大于结构优化的判据 3e-4 Ha/Bohr，说明在那个基组下 1.09 Å 的 C–H 键长远不是平衡键长。“能算完”不等于“算对了”：要比较结构或能量，先做基组和截断能的收敛测试，再做结构优化（[第4章](04-pao-vps.md)、[第6章](06-optimization-eos-neb.md)）。

### 7.7 计算时间

`.out` 最后的 `Computational Time (second)` 段列出总时间以及 `readfile`、`DFT`、`Set_Hamiltonian`、`Poisson`、`Diagonalization` 等各部分耗时（各进程的最小/最大值），用于判断瓶颈和并行效率。

### 7.8 用 Python 提取关键结果

下面是一个只用标准库的最小例子（在示例目录中运行），提取 `Utot`、HOMO/LUMO 与偶极矩并换算成 eV。更完整的解析（SCF 历史、Mulliken、力、优化历史，输出 JSON/CSV）请用仓库脚本 [`omx_parse_out.py`](../scripts/omx_parse_out.py)，用法见 `python ../scripts/omx_parse_out.py -h`。

```python
import re
import sys

HA2EV = 27.211386

def parse_out(fname):
    text = open(fname, encoding="ascii", errors="replace").read()
    res = {}
    m = re.search(r"^\s*Utot\.\s+(-?\d+\.\d+)", text, re.M)
    res["Utot_Ha"] = float(m.group(1)) if m else None
    blk = text.split("Eigenvalues (Hartree) for SCF KS-eq.")[-1]
    homo = int(re.search(r"HOMO\s*=\s*(\d+)", blk).group(1))
    eig = {int(a): float(b) for a, b, c in re.findall(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", blk, re.M)}
    res["HOMO_eV"] = eig[homo] * HA2EV
    res["LUMO_eV"] = eig.get(homo + 1, float("nan")) * HA2EV
    m = re.search(r"Absolute D\s+(\d+\.\d+)", text)
    res["dipole_D"] = float(m.group(1)) if m else None
    return res

for f in sys.argv[1:]:
    r = parse_out(f)
    print(f"{f}: Utot = {r['Utot_Ha']:.8f} Ha, HOMO = {r['HOMO_eV']:.3f} eV, "
          f"LUMO = {r['LUMO_eV']:.3f} eV, gap = {r['LUMO_eV'] - r['HOMO_eV']:.3f} eV, dipole = {r['dipole_D']} D")
```

保存为 `quick_parse.py` 后运行 `python quick_parse.py met.out h2o.out co.out`。

---

## 8. 可视化 cube 文件

`met.tden.cube`（总电子密度）和 `met.dden.cube`（差分密度）是标准的 Gaussian cube 格式，可以用：

- **OpenMX Viewer**（https://www.openmx-square.org/viewer/）：把 `met.xyz` 或 cube 文件拖进浏览器窗口，即可显示结构和等值面；
- **VESTA**：File → Open 打开 `met.dden.cube`，在 Properties → Isosurfaces 中设置正、负等值面（如 ±0.005 e/Bohr³），可以看到 C–H 键区域的电荷积累；
- **XCrySDen**：先用 OpenMX 自带的 `cube2xsf` 转成 XSF（`cube2xsf met.tden.cube` 生成 `.xsf` 文件，手册第 77 节），再 `xcrysden --xsf <生成的 .xsf 文件>`，在 Tools → Data Grid 中设置等值面。

差分密度 $\Delta\rho=\rho_{\rm SCF}-\sum_i\rho_i^{\rm atom}$ 直观地显示成键引起的电荷重新分布，是理解化学键的好工具。更多 cube 处理（密度差、平面平均、MO 可视化）见 [第10章](10-charge-mo-cube.md)，公共可视化方法见[公共篇第5章](../../Common/tutorials/05-postprocess-visualization.md)。

---

## 9. 举一反三：H2O 与 CO

示例目录中另外提供了两个分子，结构均取实验几何，设置与 `Methane.dat` 相同（Standard 基组、200 Ry、`scf.criterion 1.0e-8`），盒子取 12 Å 以减小极性分子镜像之间的静电相互作用。

### 9.1 H2O

```text
<Definition.of.Atomic.Species
  H    H6.0-s2p1      H_PBE19
  O    O6.0-s2p2d1    O_PBE19
Definition.of.Atomic.Species>
...
<Atoms.SpeciesAndCoordinates
  1   O    6.000000   6.000000   5.850000    3.0  3.0
  2   H    6.756950   6.000000   6.435870    0.5  0.5
  3   H    5.243050   6.000000   6.435870    0.5  0.5
Atoms.SpeciesAndCoordinates>
```

- 几何：O–H = 0.9572 Å、∠HOH = 104.52°（气相水分子的实验平衡结构）；
- `O_PBE19` 的价电子数为 6，所以 O 的初值是 `3.0 3.0`；总电子数 $6+2=8$，`HOMO = 4`；
- 偶极矩：实验值约 1.85 D。手册第 14 节在官方 `H2O.dat` 设置下（优化后的结构）得到约 1.83 D；本例使用实验几何，结果应在 1.8–1.9 D 附近（随基组和截断能略有变化）。偶极方向沿 z 轴（分子的 C2 轴）。

### 9.2 CO 与分子轨道输出

```text
<Atoms.SpeciesAndCoordinates
  1   C    6.000000   6.000000   5.436000    2.0  2.0
  2   O    6.000000   6.000000   6.564000    3.0  3.0
Atoms.SpeciesAndCoordinates>

MO.fileout                    on          # write molecular orbitals as cube files
num.HOMOs                     2
num.LUMOs                     2
```

- 键长 1.128 Å（实验值）；总电子数 $4+6=10$，`HOMO = 5`（5σ 轨道）；
- `MO.fileout on` 时，Cluster 计算会额外写出 `co.homo0_0.cube`（HOMO）、`co.homo0_1.cube`（HOMO−1）、`co.lumo0_0.cube`（LUMO）、`co.lumo0_1.cube`（LUMO+1）。文件名中第一个数字是自旋（0 = 上自旋），第二个是从 HOMO 往下/从 LUMO 往上数的序号。CO 的 HOMO−1 与 LUMO 都是二重简并的 π 型轨道，可以在 VESTA 中观察它们的节面；
- CO 的偶极矩很小（实验约 0.11 D，负端在 C 上），对基组和泛函非常敏感，是检验计算细节的经典例子；
- Kohn–Sham HOMO 本征值的负值并不等于电离能：GGA 的 $-\varepsilon_{\rm HOMO}$ 通常明显低于实验电离能（CO 约 14.0 eV、H2O 约 12.6 eV）。要计算电离能，应计算中性与阳离子的总能量差（ΔSCF，可配合 `scf.coulomb.cutoff on`，见手册第 61 节）。

### 9.3 批量运行与检查

```bash
for inp in Methane.dat H2O.dat CO.dat; do
  name=$(awk 'tolower($1)=="system.name"{print $2; exit}' $inp)
  mpirun -np 4 openmx $inp -nt 1 > ${name}.std 2>&1
done
bash check_results.sh            # summary: normal end, SCF steps, Utot, HOMO/LUMO, dipole
```

---

## 10. 运行中调整与重启（预告）

- **运行中修改 SCF 参数**：在计算目录中创建名为 `System.Name_SCF_keywords` 的文件（本例为 `met_SCF_keywords`），写入 `scf.maxIter`、`scf.Min.Mixing.Weight`、`scf.Max.Mixing.Weight`、`scf.Kerker.factor`、`scf.Mixing.StartPulay` 等新值。OpenMX 每个 SCF 步都会检查该文件，读到后在标准输出中打印 `The keywords for SCF iteration are renewed by ./met_SCF_keywords.`；把 `scf.maxIter` 设为负值可以让程序优雅退出。
- **重启**：`scf.restart on` 读取 `met_rst/` 中的重启文件作为初始哈密顿量与密度（`System.Name` 必须与上一次相同）；4.0 新增的 `scf.restart pp` 只做一次对角化用于 DOS/能带等后处理，且不覆盖重启文件。

详细用法见 [第5章](05-scf-kpoints-convergence.md)。

---

## 参数选择与收敛建议

| 参数 | 分子计算建议 | 说明 |
|---|---|---|
| 盒子大小 | 分子尺寸 + 2 $r_c$ 以上；极性分子再加大（12–15 Å 起） | 带电分子考虑 `scf.coulomb.cutoff on` |
| `scf.EigenvalueSolver` | `Cluster` | 分子不需要 k 点 |
| 基组 | 测试用 Quick，结果用 Standard/Precise | 见 [第4章](04-pao-vps.md) |
| `scf.energycutoff` | 150–300 Ry，做收敛测试 | 手册中 H2O 的结构与偶极矩在约 90 Ry 已收敛，但不同体系不同 |
| `scf.criterion` | 1e-8（单点、性质）；1e-9 或更严（优化、数值微分） | 默认 1e-6 偏松 |
| `scf.maxIter` | 100 | 闭壳层小分子一般 20–40 步内收敛 |
| `scf.Mixing.Type` | `RMM-DIIS`（小分子）；不收敛时试 `RMM-DIISK` | 见 [第5章](05-scf-kpoints-convergence.md) |
| `scf.SpinPolarization` | 闭壳层 `off`；自由基、O2 等开壳层 `on` 并给不对称初值 | 见 [第11章](11-spin-noncollinear.md) |
| `level.of.fileout` | 需要 cube 时 1；批量/大体系 0 | 节省磁盘 |
| 进程数 | 小分子 1–4 个 MPI 进程 | 进程过多反而变慢 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `input_open: can not open Methane.dat` | 输入文件名或路径错误 | 检查当前目录与文件名（区分大小写） |
| 读入 PAO/VPS 时报错退出，或 `Set_Cluster_UnitCell: can not open .../PAO/xxx.pao`（省略盒子时） | `DATA.PATH` 不对（默认的 `../DFT_DATA19` 在你的目录下不存在），或 PAO/VPS 名称拼写错误 | `DATA.PATH` 写绝对路径并指向含 `PAO/`、`VPS/` 的目录；`ls $DATA/PAO/H6.0.pao $DATA/VPS/H_PBE19.vps` 确认文件存在 |
| `Found an undefined species name X` / `Please check your input file` | 坐标块中的物种名未在 `Definition.of.Atomic.Species` 中定义；或文件含全角字符、BOM、UTF-16 | 核对物种名；`file x.dat`、`grep -nP '[^\x00-\x7F]' x.dat` 检查编码 |
| `Species.Number may be wrong.` | 没写 `Species.Number` 或值 ≤ 0 | 补上，且与块中行数一致 |
| `Atoms.Number may be wrong.` | 没写 `Atoms.Number` 或值 ≤ 0 | 补上，且与坐标行数一致 |
| `Format error of the sequential number N in <Atoms.SpeciesAndCoordinates` | 原子序号不连续、重复或行数与 `Atoms.Number` 不符 | 从 1 开始连续编号；核对行数 |
| `ERROR key=scf.XcType` 等，随后 `errors in the inputfile` | 选项值非法（如写成 `PBE`） | 按手册写合法值（`GGA-PBE`） |
| 输出文件出现在意料之外的目录或被覆盖 | `System.CurrrentDirectory` 拼成两个 r 被忽略；不同计算用了相同 `System.Name` | 用官方拼写；每个计算用不同的 `System.Name` |
| 修改了某个关键词但似乎没有生效 | 文件中有两处同名关键词，第一处生效；或关键词名拼错被忽略 | `grep -in "关键词" x.dat` 检查重复；对照手册拼写 |
| `.std` 末尾没有 `The calculation was normally finished.` | 程序中途退出（输入错误、内存不足、作业被杀） | 看 `.std` 的最后几十行和作业系统日志 |
| SCF 步数等于 `scf.maxIter` | 未收敛就停止了 | 增大 `scf.maxIter`；减小混合权重；增大 `scf.Mixing.History`/`StartPulay`；见第5章 |
| `Sum of MulP total` 与 `ideal(neutral)` 不一致 | 设置了 `scf.system.charge`，或计算未收敛 | 中性体系应相等；检查收敛 |
| 总能量与手册参考值不同 | 基组、截断能、坐标、盒子、泛函任一不同；或 PAO/VPS 版本不同 | 用 `Methane_manual_check.dat` 复现参考值；逐项核对设置 |
| 极性分子的偶极矩随盒子大小明显变化 | 镜像之间的静电相互作用 | 增大盒子，检查收敛 |
| 运行时 `Segmentation fault` | 栈空间不足、内存不足 | `ulimit -s unlimited`、`export OMP_STACKSIZE=512M`；减少进程数 |

---

## 练习

1. 运行 `Methane_manual_check.dat`，把 `Utot` 与手册参考值 −8.033540953061 Hartree 比较。
2. 运行 `Methane.dat`，比较它与最小基组算例的 `Utot`、HOMO 本征值和 H 原子上的力，说明基组对结果的影响。
3. 在 `Methane.dat` 中把 `System.CurrrentDirectory` 故意写成 `System.CurrentDirectory  ./out/`（事先 `mkdir out`），运行后观察输出文件写在哪里，体会“未知关键词被忽略”。
4. 把 `Methane.dat` 中 C 原子的初始电子数改为 `2.5 1.5`、把 `scf.SpinPolarization` 改为 `on`，观察 SCF 过程中 `Total Spin Moment` 的变化以及最终是否回到 0。
5. 把盒子从 10 Å 依次改为 8、12、14 Å，记录 `Utot`；再对 H2O 做同样测试并记录偶极矩，比较两者对盒子大小的敏感程度。
6. 把 `scf.energycutoff` 从 100 Ry 逐步增加到 400 Ry（步长 50 Ry），画出 CH4 的 `Utot` 随截断能的变化，观察它是否单调（提示：参见 [第5章](05-scf-kpoints-convergence.md) 与手册第 14 节）。
7. 用 VESTA 打开 `co.homo0_0.cube` 与 `co.lumo0_0.cube`，判断哪个是 σ 型、哪个是 π 型轨道。
8. （进阶）计算 H 原子和 C 原子（需自旋极化，并给出合适的初始磁矩）的能量，估算 CH4 的原子化能，并与实验值比较，思考 BSSE 与基组的影响。

---

## 参考资料

**官方文档**

- 第 5 节 测试计算（`Methane.dat` 参考输出）：https://www.openmx-square.org/openmx_man4.0/s5_test-calc.html
- 第 8 节 输入文件、8.1 甲烷示例：https://www.openmx-square.org/openmx_man4.0/s8_input-file.html ，https://www.openmx-square.org/openmx_man4.0/s8_1_example-methane-molecule.html
- 8.2 关键词说明：https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html
- 第 9 节 输出文件：https://www.openmx-square.org/openmx_man4.0/s9_output-files.html
- 第 11 节 泛函：https://www.openmx-square.org/openmx_man4.0/s11_functionals.html
- 12.1 基组概述、12.3 推荐基组表：https://www.openmx-square.org/openmx_man4.0/s12_1_basis-func-overview.html ，https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html
- 13.1 赝势与价电子数：https://www.openmx-square.org/openmx_man4.0/s13_1_valence-pseudo.html
- 14.1 截断能收敛（CH4、H2O）：https://www.openmx-square.org/openmx_man4.0/s14_1_convergence.html
- 16.3 运行中调整 SCF 参数：https://www.openmx-square.org/openmx_man4.0/s16_3_otf-scf.html
- 17.1 重启（`scf.restart on|off|pp|c2n`）：https://www.openmx-square.org/openmx_man4.0/s17_1_restarting-general.html
- 第 61 节 电离能与电子亲和能（`scf.coulomb.cutoff`）：https://www.openmx-square.org/openmx_man4.0/s61_ip-ea.html
- 第 69 节 盒子大小自动确定：https://www.openmx-square.org/openmx_man4.0/s69_auto-cell-size.html
- OpenMX Viewer：https://www.openmx-square.org/viewer/

**相关 GitHub 仓库**（只链接，不转载其内容）

- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，含输入文件入门与 Python 后处理示例。
- [hirtatsu/Calculations](https://github.com/hirtatsu/Calculations)：日文实操笔记，`OpenMX/` 目录中逐段讲解输入文件并演示用 tmux/nohup 运行与检查日志。
- [samlipton/pyfesm](https://github.com/samlipton/pyfesm)：轻量 Python 解析器，读取 SCF/体系信息、能带、Mulliken 等，可参考其输出解析思路。
- [wxia529/acoord](https://github.com/wxia529/acoord)：在 VS Code 中编辑与可视化 `.dat` 结构，保留参数块。
- [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx)：POSCAR → `.dat` 转换（自动设置 k 网格、基组档次、初始磁矩）。

**关键论文**

- T. Ozaki, *Phys. Rev. B* 67, 155108 (2003). DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)
- T. Ozaki, H. Kino, *Phys. Rev. B* 69, 195113 (2004). DOI: [10.1103/PhysRevB.69.195113](https://doi.org/10.1103/PhysRevB.69.195113)
- T. Ozaki, H. Kino, *Phys. Rev. B* 72, 045121 (2005)（总能量分解与投影展开）. DOI: [10.1103/PhysRevB.72.045121](https://doi.org/10.1103/PhysRevB.72.045121)
- I. Morrison, D. M. Bylander, L. Kleinman, *Phys. Rev. B* 47, 6728 (1993). DOI: [10.1103/PhysRevB.47.6728](https://doi.org/10.1103/PhysRevB.47.6728)
- J. P. Perdew, K. Burke, M. Ernzerhof, Generalized gradient approximation made simple, *Phys. Rev. Lett.* 77, 3865 (1996). DOI: [10.1103/PhysRevLett.77.3865](https://doi.org/10.1103/PhysRevLett.77.3865)

---

[← 上一章](02-install.md) | [返回目录](../README.md) | [下一章 →](04-pao-vps.md)
