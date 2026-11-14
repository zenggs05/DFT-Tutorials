# 第1章 环境准备与第一个 Gaussian 计算（Getting Started with Gaussian 16 and GaussView 6）

本章从零开始讲 Gaussian 16 量子化学计算：Linux 服务器与 Windows（G16W）上的安装与环境变量（g16root、GAUSS_SCRDIR），GaussView 6 建模和生成 gjf 输入文件，再完成第一个水分子 HF 单点能和 B3LYP 几何优化。读完本章，你能判断一个 Gaussian 作业是否正常结束（Normal termination / Error termination），读懂 log 输出里最重要的几个段落，也知道 .gjf、.log、.chk、.fchk、.cub 各是什么文件。

---

## 1. 本章目标与适用场景

| 你将学会 | 对应小节 |
|---|---|
| 区分 Gaussian 16 的各个修订版（A.03/B.01/C.01/C.02）以及 Windows 版与 Linux 版 | §2 |
| 在 Linux 上配置 `g16root`、`g16.profile`、`GAUSS_SCRDIR`，并检查安装是否可用 | §4 |
| 用命令行运行 Gaussian（`g16 < job.gjf > job.log`），或在 Windows 上用 G16W/GaussView 提交 | §5 |
| 用 GaussView 6 画分子、Clean、Symmetrize、设置计算并导出 gjf | §6 |
| 运行第一个单点能和第一个几何优化 | §7、§10 |
| 判断作业是否成功，失败时定位出错的 Link | §8 |
| 读懂 log 中的 Input/Standard orientation、SCF Done、Mulliken 电荷、偶极矩、archive 段 | §9 |
| 认识 Gaussian 的各种文件类型，避开 CRLF 换行符和中文路径这类常见坑 | §11、§12 |

适用读者：第一次接触 Gaussian 的本科生和研究生，或者从 VASP/CP2K 等周期性程序转过来做分子计算的读者。Linux 命令行基础（cd、ls、vi/nano、ssh、scp）请先看 [第1章 计算环境搭建](../../Common/tutorials/01-linux-hpc.md)，集群作业调度见 [第2章 HPC 作业调度](../../Common/tutorials/02-job-scheduler.md)。

> [!NOTE]
> Gaussian 是商业软件，需要所在单位购买许可证（站点许可或课题组许可）。本教程不提供任何安装包，只讲已获得合法授权后怎么配置和使用。

---

## 2. Gaussian 16 与 GaussView 6 概览

### 2.1 版本沿革

到 2026 年 11 月，Gaussian 的最新版本仍是 **Gaussian 16 Rev. C.02**，官方没有发布更新的大版本。各修订版要点如下：

| 修订版 | 主要变化（摘自官方 Release Notes） | 对日常使用的影响 |
|---|---|---|
| A.03 | G16 首个广泛发布版本；相对 G09，默认积分格点改为 UltraFine，积分精度改为 $10^{-12}$，SCRF 默认改为对称形式 IEFPCM，物理常数改用 2010 版，默认 `%Mem` 改为 800 MB | 与 G09 结果比较时会有小差别 |
| B.01 | 新泛函 M08HX、MN15、MN15L、PW6B95、PW6B95D3；双杂化 DSDPBEP86、PBE0DH、PBEQIDH；PM7；`Pop=DCT`；`Geom=GIC`；TD 态 `Freq=Raman`；`Opt=Recalc`；Link 0 命令有了等价的命令行选项和环境变量；支持 P100 GPU | 第3章、第10章会用到这些新功能 |
| C.01（2019） | 支持 NBO 7（`Pop=NBO7` 等）；`Pop=(MK,Resp=N)`；`Pop=SaveCM5/SaveHirshfeld`；工具程序支持 `-m` 内存参数（如 `formchk -m=1gb`）；`Opt=NGoUp/NGoDown/MaxEStep`；chk 文件带版本号；Linda 9.2；V100 GPU | 目前最常见的版本 |
| C.02 | 功能与 C.01 相同，增加 NVIDIA A100 GPU 支持，并改用 NVIDIA SDK 21.3 编译器构建 | 有 A100 才需要特意升级 |

查看自己用的是哪个修订版：每个 log 文件开头都会打印版本信息，例如：

```text
 ******************************************
 Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019
                10-Nov-2026
 ******************************************
```

（示意输出，日期与平台字段仅供参考）`ES64L` 表示 x86-64 Linux 版本，`G16RevC.01` 就是修订版号。发表文章时，引用格式中的 Revision 要与这里一致。

```bash
# Linux: show the revision used by an existing log file
grep -m1 "Gaussian 16:" job.log
```

### 2.2 Windows 版（G16W）与 Linux 版的区别

| 项目 | Gaussian 16W（Windows） | Gaussian 16（Linux/Unix） |
|---|---|---|
| 典型使用场景 | 个人电脑、教学、小分子 | 服务器、超算集群、生产计算 |
| 提交方式 | G16W 图形界面（File → Open 打开 gjf 后运行），或在 GaussView 中直接 Submit | 命令行 `g16`，配合 Slurm/PBS 等调度系统 |
| 输出扩展名 | 通常为 `.out` | 通常为 `.log` |
| 并行 | 单机共享内存并行 | 共享内存并行；另购 Linda 可跨节点；Linux 下支持 NVIDIA GPU |
| 常见坑 | 路径含中文/空格、杀毒软件拦截、`%mem` 超过可用物理内存 | 没有 source `g16.profile`、scratch 目录不存在或空间不足 |
| 输入/输出格式 | 与 Linux 版完全相同，gjf 可以互相拷贝 | 同左 |

> [!TIP]
> 一个常见的高效组合是：Windows 上用 GaussView 建模、检查结构和看结果，计算放到 Linux 服务器上跑。两边的 .gjf 和 .log 通用，二进制 .chk 不通用，跨平台要先转成 .fchk（见 §11）。

### 2.3 GaussView 6

GaussView 是 Gaussian 官方的图形界面，目前版本为 GaussView 6（最新的小版本是 6.1.1）。它的三类用途：

1. **建模**：画分子、调键长键角、Clean 预优化、对称化、定义片段（Fragment）、设置 ONIOM 分层和冻结原子。
2. **生成输入并提交**：Calculate → Gaussian Calculation Setup 填写各个标签页后，保存 gjf 或直接提交。
3. **结果可视化**：打开 .log/.out/.fchk，查看优化轨迹、能量曲线、振动动画、IR/Raman/UV-Vis/NMR 谱、分子轨道、静电势面、扫描与 IRC 曲线。

GaussView 也要单独授权。不想用 GaussView 时，可以用免费工具替代：Avogadro 2（建模、导出 Gaussian 输入）、Multiwfn + VMD（波函数分析和作图，见 [第16章](16-postprocess-multiwfn.md)），以及 IQmol、Chemcraft 等。

---

## 3. 原理速览：Gaussian 是怎么"跑"一个作业的

### 3.1 一次计算求的是什么

对一个固定的分子构型（核坐标 $\{\mathbf{R}_A\}$），Gaussian 求解电子的 Hartree–Fock 方程或 Kohn–Sham 方程，得到总能量 $E(\{\mathbf{R}_A\})$，单位是 Hartree（$E_h$）：

$$
E_{\text{tot}} = E_{\text{elec}} + V_{NN},\qquad V_{NN}=\sum_{A<B}\frac{Z_A Z_B}{R_{AB}}
$$

分子轨道 $\psi_i$ 展开为原子中心高斯函数（基组）$\chi_\mu$ 的线性组合：$\psi_i=\sum_\mu C_{\mu i}\chi_\mu$。迭代求解系数 $C$，直到能量和密度自洽，这一步叫 SCF。几何优化在此基础上求能量对坐标的梯度，并移动原子去找能量极小点；频率计算再求二阶导数（Hessian）。方法和基组怎么选见 [第3章](03-method-basis.md)。

常用能量换算（2010 CODATA，即 G16 使用的常数）：

| 1 Hartree 等于 | 数值 |
|---|---|
| kcal/mol | 627.5095 |
| kJ/mol | 2625.500 |
| eV | 27.21138 |

### 3.2 Link（链接程序）结构

Gaussian 由许多独立的可执行程序组成，叫做 **Link**，文件名是 `lNNN.exe`。route 行决定按什么顺序调用哪些 Link。出错时日志末尾会写 `Error termination via Lnk1e in .../l502.exe`，告诉你是哪个 Link 出的错。熟悉几个常见 Link 的分工，看报错会快很多：

| Link | 主要职责 | 常见报错举例 |
|---|---|---|
| L0（`g16` 本身） | 读取 `%` 开头的 Link 0 命令，启动后续 Link | 找不到 l1.exe（环境变量没设好） |
| L1 | 解析 route 行，生成内部"非标准 route" | `QPErr --- A syntax error was detected in the input line.` |
| L101 | 读取 title、电荷/多重度、分子坐标 | `End of file in ZSymb.` |
| L103 | Berny 几何优化器 | `FormBX had a problem`、`Linear angle in Bend` |
| L202 | 坐标重定向、对称性识别、距离矩阵 | `Small interatomic distances encountered` |
| L301 | 构建基组、ECP，检查电荷与多重度，设置色散参数 | 多重度与电子数不匹配、`Atomic number out of range for ... basis set` |
| L401 | 产生初猜轨道 | 读取 chk 失败 |
| L502 | SCF 迭代 | `Convergence failure -- run terminated.` |
| L601 | 布居分析（Mulliken、偶极矩等） | — |
| L9999 | 写 archive 段、结束作业 | 优化步数超限 `Number of steps exceeded` 在这里报出 |

完整报错大全见 [第18章](18-errors-capstone.md)。

---

## 4. Linux 安装与环境配置

### 4.1 安装（管理员或个人目录安装）

官方发布的二进制包是 `*.tbJ` 文件（xz 压缩的 tar 包）。`g16root` 指向**包含 `g16` 目录的那一层**，不是 `g16` 目录本身。例如程序装在 `/opt/g16/`，则 `g16root=/opt`。

```bash
# 1) choose an install prefix; g16root is the PARENT of the g16 directory
export g16root=$HOME/soft          # e.g. installs into $HOME/soft/g16
mkdir -p $g16root && cd $g16root

# 2) unpack the distribution (file name differs by platform/revision)
tar xvJf /path/to/G16-C02-xxxx.tbJ

# 3) restrict permissions (Gaussian license requires that only licensed
#    users can run it); use your own group name instead of 'g16users'
chgrp -R g16users g16      # skip if you install only for yourself
chmod -R o-rwx g16

# 4) run the install script shipped with the binary
cd g16
./bsd/install
```

> [!WARNING]
> 不要先在 Windows 上解压 .tbJ 再传到 Linux。Windows 解压会丢掉可执行权限和符号链接，安装后会出现各种 "Permission denied" 或找不到文件的错误。直接把 .tbJ 原样传到 Linux 再解压。

### 4.2 环境变量

把下面几行加到 `~/.bashrc`（csh/tcsh 用户改用 `setenv` 并 `source $g16root/g16/bsd/g16.login`）：

```bash
# ---- Gaussian 16 ----
export g16root=$HOME/soft                  # directory that CONTAINS g16/
export GAUSS_SCRDIR=/scratch/$USER/g16scr  # large, fast, local disk
mkdir -p "$GAUSS_SCRDIR"
source $g16root/g16/bsd/g16.profile        # sets GAUSS_EXEDIR, PATH, etc.
```

| 变量 | 由谁设置 | 作用 |
|---|---|---|
| `g16root` | 用户 | Gaussian 安装位置（`g16` 目录的上一级） |
| `GAUSS_SCRDIR` | 用户 | 临时文件（`Gau-12345.rwf` 等）存放目录；不设时用当前工作目录 |
| `GAUSS_EXEDIR` | `g16.profile` | 各个 `lNNN.exe` 所在目录，默认包含 `$g16root/g16` |
| `G16BASIS` | `g16.profile` | 内置基组文件（.gbs）目录，配合 `@` 包含文件使用（第3章） |
| `GAUSS_MDEF`、`GAUSS_CDEF` | 用户（可选） | 默认内存、默认 CPU 列表，相当于每个作业都写了 `%Mem`、`%CPU` |

`GAUSS_SCRDIR` 的要求：

- **目录必须事先存在**，否则作业一启动就失败。
- **空间要足够**。MP2、双杂化、大体系频率计算的 .rwf 可能有几十 GB。
- **尽量用计算节点本地盘**（`/scratch`、`/tmp`、`/local`），不要用网络共享的 home 目录，否则 I/O 会拖慢计算。
- 作业被 kill 或异常退出时，scratch 里的 `Gau-*` 文件**不会自动删除**，要定期清理：

```bash
ls -lh $GAUSS_SCRDIR
du -sh $GAUSS_SCRDIR
# delete leftover scratch files older than 3 days (make sure no job is running!)
find $GAUSS_SCRDIR -name 'Gau-*' -mtime +3 -delete
```

集群上的推荐做法是每个作业建一个独立的 scratch 子目录，作业结束时自动清理，模板见 [第17章](17-hpc-jobs.md)。

### 4.3 验证安装

```bash
source ~/.bashrc
which g16                 # should print $g16root/g16/g16
echo $GAUSS_EXEDIR        # should contain $g16root/g16
ls $g16root/g16/l1.exe    # the first link must exist
df -h $GAUSS_SCRDIR       # free space of the scratch disk
testrt "#p B3LYP/6-31G(d) Opt Freq"   # route syntax checker; prints the internal route
```

`testrt` 是 Gaussian 自带的 route 语法检查工具。route 合法时它打印对应的"非标准 route"（`1/18=20,...;` 这类数字行），有语法错误时打印 `QPErr` 和出错位置。在作业提交前用它检查 route 行是个好习惯。

> [!NOTE]
> 命令行选项与环境变量是 B.01 引入的：`g16 -m=8GB -c="0-7" job.gjf` 等价于在输入里写 `%Mem=8GB` 和 `%CPU=0-7`，环境变量 `GAUSS_MDEF`、`GAUSS_CDEF` 的作用也一样。优先级从高到低：输入文件中的 Link 0 行 > 命令行选项 > 环境变量 > `Default.Route` 文件。`%NProcShared` 对应的 `-p`/`GAUSS_PDEF` 已被官方标为过时（deprecated），推荐改用 `-c`/`GAUSS_CDEF`。细节见 [第2章](02-input-anatomy.md) 和 [第17章](17-hpc-jobs.md)。

### 4.4 Windows（G16W）安装注意事项

- **安装路径和工作目录都不要包含中文或空格**，例如用 `D:\G16W` 和 `D:\calc\h2o`，不要用 `C:\Program Files\` 或 `桌面\新建文件夹`。路径问题是 Windows 下出现 "Link 1 not found" 或 `Severe Error Message #2070` 等无法启动错误的头号原因。
- 杀毒软件（包括系统自带的 Defender）可能隔离 `l*.exe`，必要时把 G16W 目录加入白名单。
- 在 G16W 的设置（Preferences）里设置 scratch 目录，同样不要有中文和空格。
- `%mem` 不要超过物理内存的一半到 2/3，Windows 下内存申请过大同样会导致启动失败。

### 4.5 WSL 的情况

在 Windows 10/11 的 WSL2 Ubuntu 中，可以按 §4.1–4.3 安装 **Linux 版** Gaussian，前提是你的许可证允许在该机器上使用 Linux 版。WSL 的配置见 [第1章 计算环境搭建](../../Common/tutorials/01-linux-hpc.md)。建议把计算目录放在 WSL 自己的文件系统（如 `~/calc`）里，不要放在 `/mnt/c/...` 下，后者 I/O 很慢。

---

## 5. 运行 Gaussian 的几种方式

### 5.1 Linux 命令行

```bash
# style 1: redirect stdin/stdout (works with any file name/extension)
g16 < h2o_hf_sp.gjf > h2o_hf_sp.log

# style 2: give the job name; reads h2o_hf_sp.gjf, writes h2o_hf_sp.log
g16 h2o_hf_sp.gjf

# run in background and keep it alive after you log out
nohup g16 < h2o_b3lyp_opt.gjf > h2o_b3lyp_opt.log 2>&1 &

# override memory/cores from the command line (B.01+)
g16 -m=4GB -c="0-3" < h2o_b3lyp_opt.gjf > h2o_b3lyp_opt.log

# follow the output while it runs
tail -f h2o_b3lyp_opt.log
```

两种写法的区别：`g16 job.gjf` 会自动把输出命名为 `job.log`；重定向写法可以随意指定输出名。两种都要在**输入文件所在目录**下执行，因为 `%chk=xxx.chk` 这类相对路径以当前目录为基准。

> [!IMPORTANT]
> 共享登录节点（login node）上只能跑几秒钟的小测试。正式计算必须通过 Slurm/PBS 提交到计算节点，否则会影响其他用户，也可能被管理员终止。本章示例目录里提供了最小的 Slurm 脚本，完整模板见 [第17章](17-hpc-jobs.md)。

### 5.2 Windows

- **G16W 界面**：File → Open 选择 .gjf，确认输出文件名后点 Run；进度会显示在窗口里。
- **GaussView**：Calculate → Gaussian Calculation Setup 设置好后点 **Submit**，GaussView 会先让你保存 .gjf，再调用 G16W 运行。运行结束后可以直接打开结果。

### 5.3 本教程配套脚本

- 示例目录 [01-first-job-h2o](../examples/01-first-job-h2o/README.md) 里的 `run_g16_local.sh`：在本机依次运行若干 gjf，并检查每个作业是否 Normal termination。
- 仓库脚本 [`g16_log_summary.py`](../scripts/g16_log_summary.py)：批量汇总 log 的结束状态、最终能量、虚频数等，用 `python ../scripts/g16_log_summary.py -h` 查看用法。
- 仓库脚本 [`sub_g16.sh`](../scripts/sub_g16.sh)：集群提交封装（[第17章](17-hpc-jobs.md)详细介绍）。

---

## 6. GaussView 6 建模实战

下面以乙醇为例（水分子太简单，体现不出 Clean 和对称化的作用）。

### 6.1 画分子

1. 打开 GaussView，主窗口叫 *Builder*，空白的分子窗口会自动打开。
2. 点工具栏的 **Element Fragment**（元素周期表图标），选 C，再在下方选 sp3 碳的片段（带 4 个 H 的甲烷样式）。
3. 在分子窗口空白处单击，放下第一个 CH4。再在其中一个 H 上单击，就会用一个新的 CH3 替换这个 H，得到乙烷。
4. 切换元素为 O（选带 2 个价位的 O 片段），单击乙烷的一个 H，替换成 OH，得到乙醇。
5. 也可以用 **Ring Fragment**（苯环等）和 **R-Group Fragment**（官能团）快速搭建复杂分子。
6. 用 **Modify Bond / Angle / Dihedral** 手动调整几何参数，用 **Inquire** 查看距离和角度。

> [!TIP]
> 环较多或较复杂的分子，可以先在 ChemDraw 或在线工具里画好，导出 .mol/.sdf 再用 GaussView 打开；也可以从 PubChem 下载 3D SDF。数据库下载与格式转换见 第3章 结构准备。

### 6.2 Clean 与对称化

- **Edit → Clean**：用简单的力场规则快速整理几何（键长、键角归一化）。手画的结构一定要 Clean 一次，否则初始结构太差，优化会走很多步，甚至失败。Clean 不是量子化学优化，只是让初始结构"像样"。
- **Edit → Symmetrize**（或 **Edit → Point Group** 打开点群对话框）：识别并强制分子的点群对称性。乙醇在 Clean 后通常接近 Cs 对称。

对称化要注意两点：

1. Gaussian 的几何优化**保持初始结构的点群**。如果把一个实际上不对称的分子强行对称化，优化可能停在高对称的鞍点上。所以 Opt 后一定要做频率计算，确认没有虚频（[第6章](06-freq-thermo.md)）。
2. 对称性能节省计算量。小分子有明显对称性（H2O 是 C2v，苯是 D6h）时对称化是好事；对称性只是"近似"时，不要用过大的容差硬凑。

### 6.3 Calculate → Gaussian Calculation Setup 各标签页

| 标签页 | 作用 | 第一次使用时的建议 |
|---|---|---|
| **Job Type** | 任务类型：Energy（单点）、Optimization、Frequency、Opt+Freq、IRC、Scan、NMR、Stability 等，以及相应子选项（如优化到 Minimum 还是 TS） | 选 Optimization 或 Opt+Freq |
| **Method** | Ground State / Excited State；方法（HF、DFT 及泛函、MP2 等）；基组下拉组合；电荷与自旋多重度；Restricted/Unrestricted；底部的 **Additional Keywords** 可以填任意额外关键词 | 先选 DFT / B3LYP / 6-31G(d)，在 Additional Keywords 里加 `EmpiricalDispersion=GD3BJ` |
| **Title** | 标题行（不参与计算，只出现在输出里） | 写清体系与目的 |
| **Link 0** | 内存（%mem）、共享处理器数（%nprocshared）、检查点文件名（%chk）；GaussView 6 新增 **OldChk** 字段 | 务必勾选保存 chk，并按机器实际情况填内存和核数 |
| **General** | 常用通用选项：是否写连接表、Counterpoise、写片段信息、忽略对称性等 | 一般保持默认 |
| **Guess** | 初猜方法（Default、Huckel、Read 等），**Use fragments** 对应 `Guess=Fragment` | 默认即可 |
| **NBO** | 选择 NBO 分析类型 | 第12章再用 |
| **PBC** | 周期性边界条件（Gaussian 的 PBC 功能） | 分子计算不用 |
| **Solvation** | 隐式溶剂模型（PCM/SMD）与溶剂选择 | 第9章再用 |
| **Add. Inp.** | 坐标之后的附加输入段（如 Gen 基组、ModRedundant 等） | 第3章、第5章会用到 |
| **Preview** | 预览将要写出的完整输入文件 | **保存前一定要看一眼** |

GaussView 6 还可以基于已有 chk 设置后续作业，例如从 chk 读结构、读初猜、读力常数。对应关键词 `Geom=Check`、`Guess=Read`、`Opt=ReadFC`，第2章会讲它们的文本写法。

### 6.4 保存 gjf

在 Calculation Setup 中设置好以后，用 File → Save 保存，文件类型选 *Gaussian Input Files (*.gjf *.com)*。也可以直接点 Submit，GaussView 会先提示保存。保存时注意：

- 勾选保存对话框中的 **Write Cartesians**，坐标会写成笛卡尔坐标，更直观，也便于脚本处理。
- GaussView 默认在 route 中加 `Geom=Connectivity`，并在坐标后写一段**连接表**（每行是原子编号加成键原子与键级）。这段对普通 DFT 计算没有用，但删的时候必须连同 route 里的 `Geom=Connectivity` 一起删，否则会读错（详见 [第2章](02-input-anatomy.md)）。
- 文件名和路径只用英文、数字、下划线和短横线。

GaussView 导出的乙醇输入大致如下（坐标仅示意）：

```text
%nprocshared=4
%mem=4GB
%chk=ethanol.chk
# opt freq b3lyp/6-31g(d) empiricaldispersion=gd3bj geom=connectivity

Title Card Required

0 1
 C                 -0.97710000    0.44017000    0.00000000
 C                  0.54056000    0.39978000    0.00000000
 ...

 1 2 1.0 4 1.0 5 1.0 6 1.0
 2 3 1.0 7 1.0 8 1.0
 ...

```

> [!TIP]
> GaussView 的默认标题是 "Title Card Required"，请改成有意义的描述（例如 `ethanol opt freq, gas phase`），以后翻几百个 log 时会感谢自己。

### 6.5 用 GaussView 看结果

打开 .log 或 .fchk 后，**Results** 菜单里有：

| 菜单项 | 用途 |
|---|---|
| Summary | 方法、基组、能量、偶极矩、点群等汇总 |
| Charge Distribution | 原子电荷着色显示 |
| Optimization | 优化过程中能量和梯度随步数的曲线，可以播放结构变化 |
| Vibrations | 频率列表、振动动画、IR/Raman 谱 |
| Surfaces/Contours | 分子轨道、电子密度、静电势面（需要 chk/fchk） |
| Scan、IRC/Path | 扫描和 IRC 的能量曲线 |
| UV-VIS、NMR | 激发态和 NMR 结果 |

---

## 7. 第一个计算：水分子 HF/6-31G(d) 单点能

### 7.1 输入文件

文件 [`h2o_hf_sp.gjf`](../examples/01-first-job-h2o/h2o_hf_sp.gjf)：

```text
%chk=h2o_hf_sp.chk
%mem=1GB
%nprocshared=2
#p HF/6-31G(d) SP

water single point, my first Gaussian job

0 1
O     0.000000    0.000000    0.119262
H     0.000000    0.763239   -0.477047
H     0.000000   -0.763239   -0.477047

```

逐行解释：

| 行 | 内容 | 说明 |
|---|---|---|
| 1 | `%chk=h2o_hf_sp.chk` | Link 0 命令：保存二进制检查点文件，用于后续读取结构和波函数、作图 |
| 2 | `%mem=1GB` | 最多使用 1 GB 内存（G16 默认 800 MB）；单位直接跟在数字后面，不加空格 |
| 3 | `%nprocshared=2` | 使用 2 个 CPU 核做共享内存并行 |
| 4 | `#p HF/6-31G(d) SP` | route 行：`#p` 表示输出更详细，方法/基组为 HF/6-31G(d)，任务类型 SP（单点，也是默认任务，可省略） |
| 5 | （空行） | route 段以空行结束 |
| 6 | `water single point, ...` | title 段，最多 5 行，只做注释 |
| 7 | （空行） | title 段以空行结束 |
| 8 | `0 1` | 电荷 0，自旋多重度 1（单重态） |
| 9–11 | 元素 + x y z | 笛卡尔坐标，单位 Å（初始结构来自 G2 数据集的几何） |
| 12 | （空行） | 分子坐标段以空行结束；**文件最后必须有空行** |

> [!WARNING]
> 漏掉最后的空行是新手最常见的错误，结果是 L101 报 `End of file in ZSymb.`。养成习惯：每个 gjf 的末尾至少留一个空行。

### 7.2 运行

```bash
cd Gaussian/examples/01-first-job-h2o
g16 < h2o_hf_sp.gjf > h2o_hf_sp.log
# or: bash run_g16_local.sh h2o_hf_sp.gjf
```

这个作业在普通笔记本上只需一两秒。

### 7.3 检查结果

```bash
tail -n 3 h2o_hf_sp.log                 # last lines: Normal termination?
grep "SCF Done" h2o_hf_sp.log           # final SCF energy (Hartree)
grep -A4 "Mulliken charges:" h2o_hf_sp.log
grep -A1 "Dipole moment" h2o_hf_sp.log
```

（示意输出，数值仅供参考）

```text
 SCF Done:  E(RHF) =  -76.0107465     A.U. after    9 cycles
 Mulliken charges:
               1
     1  O   -0.866
     2  H    0.433
     3  H    0.433
 Dipole moment (field-independent basis, Debye):
    X=              0.0000    Y=              0.0000    Z=             -2.1990  Tot=              2.1990
```

能量约 −76.01 Hartree 量级、偶极矩约 2 Debye 量级就说明计算正常。实验偶极矩是 1.85 D（NIST CCCBDB），HF 一般会高估。

---

## 8. 判断作业是否成功

### 8.1 Normal termination

正常结束的作业，最后一行是：

```text
 Normal termination of Gaussian 16 at Mon Nov 10 10:00:00 2026.
```

注意：

- `Opt Freq` 这类组合任务在内部会拆成两个作业步（优化完成后自动做频率），所以**正常结束时会有两行** Normal termination。用 `--Link1--` 写的多步作业也一样，每步一行。
- 判断"整个作业成功"要看**最后一行**，同时确认没有 `Error termination`：

```bash
grep -c "Normal termination" job.log     # number of finished job steps
grep "Error termination" job.log         # should print nothing
tail -n 1 job.log
```

- 没有 Normal termination，也没有 Error termination：作业可能还在运行，也可能被调度系统 kill（超时、超内存）或节点宕机。去看调度系统的输出文件（如 `slurm-12345.out`）。

### 8.2 Error termination：定位出错的 Link

失败时，日志末尾形如（示意输出）：

```text
 Error termination via Lnk1e in /opt/g16/l502.exe at Mon Nov 10 10:05:12 2026.
 Job cpu time:       0 days  0 hours  3 minutes 12.0 seconds.
 Elapsed time:       0 days  0 hours  0 minutes 48.5 seconds.
```

排查步骤：

1. 从 `lNNN.exe` 看是哪个 Link 出错（对照 §3.2 的表）。
2. **往上翻 10–30 行**，真正的错误信息在 Error termination 上面，例如 `Convergence failure -- run terminated.`：

```bash
tail -n 40 job.log
grep -n "Error termination" job.log
grep -B 15 "Error termination" job.log | head -40
```

3. 按报错信息查原因。最常见的几种：

| 出错 Link | 典型报错 | 首先检查 |
|---|---|---|
| l1.exe | `QPErr --- A syntax error was detected in the input line.` | route 拼写、括号是否配对；用 `testrt` 检查 |
| l101.exe | `End of file in ZSymb.` | 空行是否齐全、电荷多重度行是否存在、是否为 CRLF 换行 |
| l202.exe | `Small interatomic distances encountered` | 原子重叠、坐标单位错（Bohr/Å）、重复原子 |
| l301.exe | `The combination of multiplicity 1 and 17 electrons is impossible.` | 电荷与多重度 |
| l502.exe | `Convergence failure -- run terminated.` | 结构、电荷、多重度；SCF 收敛技巧见 [第4章](04-sp-scf.md) |
| l9999.exe | `Number of steps exceeded` | 优化不收敛，见 [第5章](05-opt.md) |

仓库脚本 [`g16_error_doctor.py`](../scripts/g16_error_doctor.py) 能自动识别出错的 Link 和常见报错，并给出中文建议（用 `python ../scripts/g16_error_doctor.py -h` 查看用法）。报错大全见 [第18章](18-errors-capstone.md)，也可以查社区整理的 [Solutions for Gaussian Error Messages](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/)。

---

## 9. 输出文件导览

Gaussian 的 log 是纯文本，按程序执行顺序写出。下面按出现顺序介绍最重要的几段（摘录均为示意输出，数值仅供参考）。

### 9.1 头部：版本、Link 0 与 route 回显

```text
 Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019
 %chk=h2o_hf_sp.chk
 %mem=1GB
 %nprocshared=2
 Will use up to    2 processors via shared memory.
 ----------------------------
 #p HF/6-31G(d) SP
 ----------------------------
 1/38=1,57=2/1;
 2/12=2,17=6,18=5,40=1/2;
 3/6=3,11=9,25=1,30=1/1,2,3;
 ...
```

route 下面的数字行就是 L1 生成的"非标准 route"（`overlay/选项=值/要执行的 links;`）。一般不用看懂它，但它说明 route 已经被正确解析。

### 9.2 Input orientation 与 Standard orientation

```text
                          Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.119262
      2          1           0        0.000000    0.763239   -0.477047
      3          1           0        0.000000   -0.763239   -0.477047
 ---------------------------------------------------------------------
 ...
 Full point group                 C2V     NOp   4
 ...
                         Standard orientation:
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.119262
 ...
```

- **Input orientation**：你输入的坐标（优化过程中每一步都会打印一次）。
- **Standard orientation**：Gaussian 为利用对称性，把分子平移到电荷中心并旋转到标准取向后的坐标。**分子轨道、偶极矩分量、多极矩等都在标准取向下给出。**
- 加 `NoSymm` 关键词后，Gaussian 不再旋转分子，只打印 Input orientation，所有量都在输入取向下给出。对比结构或使用外部坐标系时有用（第2章）。
- `Full point group` 行给出识别出的点群。

> [!TIP]
> 提取优化后的结构时，取**最后一个** `Standard orientation`（或 `Input orientation`，用 NoSymm 时）块。可以用仓库脚本 [`log2xyz.py`](../scripts/log2xyz.py)（`python ../scripts/log2xyz.py -h`），或直接用 GaussView 打开 log，它默认显示最后一个结构。

### 9.3 基组信息与 SCF

```text
 Standard basis: 6-31G(d) (6D, 7F)
 ...
    19 basis functions,    36 primitive gaussians,    19 cartesian basis functions
     5 alpha electrons        5 beta electrons
       nuclear repulsion energy         9.1681932561 Hartrees.
 ...
 SCF Done:  E(RHF) =  -76.0107465184     A.U. after    9 cycles
```

- `basis functions` 是基函数数目 $N$，决定计算量：HF/DFT 大致按 $N^3$–$N^4$ 增长，MP2 按 $N^5$。
- `5 alpha electrons 5 beta electrons`：可以用来核对电荷和多重度是否符合预期。
- `SCF Done`：**SCF 收敛后的总能量**，单位 Hartree（A.U.）。括号里的 `RHF`、`RB3LYP`、`UB3LYP` 表明方法和限制/非限制类型。HF/DFT 单点的最终能量就是这一行；MP2、双杂化等后 SCF 方法的最终能量在其后另行打印（第3章）。
- 用 `#p` 时，SCF 的每一轮迭代都会打印出来，便于诊断不收敛（第4章）。

### 9.4 轨道能级与布居分析

```text
 Orbital symmetries:
       Occupied  (A1) (A1) (B2) (A1) (B1)
 ...
 The electronic state is 1-A1.
 Alpha  occ. eigenvalues --  -20.55848  -1.34616  -0.71421  -0.57040  -0.49790
 Alpha virt. eigenvalues --    0.21197   0.30472   1.01672   1.09280   1.13456
 ...
 Mulliken charges:
               1
     1  O   -0.866
     2  H    0.433
     3  H    0.433
 Sum of Mulliken charges =   0.00000
```

- 轨道能量单位为 Hartree。最后一个 `occ.` 值是 HOMO，第一个 `virt.` 值是 LUMO。
- Mulliken 电荷强烈依赖基组，只适合定性比较。更可靠的电荷（NPA、Hirshfeld/CM5、ESP 电荷）见 [第12章](12-population-charges.md)。

### 9.5 偶极矩与多极矩

```text
 Dipole moment (field-independent basis, Debye):
    X=              0.0000    Y=              0.0000    Z=             -2.1990  Tot=              2.1990
```

偶极矩的分量在标准取向下给出，`Tot` 是大小（Debye）。

### 9.6 archive 段（@ 段）

作业最后会打印一段以 `1\1\` 开头、以 `\\@` 结束的压缩摘要，称为 archive entry：

```text
 1\1\GINC-NODE01\SP\RHF\6-31G(d)\H2O1\USER\10-Nov-2026\0\\#p HF/6-31G(d)
  SP\\water single point, my first Gaussian job\\0,1\O,0.,0.,0.119262\H
 ,0.,0.763239,-0.477047\H,0.,-0.763239,-0.477047\\Version=ES64L-G16RevC.0
 1\State=1-A1\HF=-76.0107465\RMSD=...\Dipole=0.,0.,-0.8651...\PG=C02V [C2
 (O1),SGV(H2)]\\@
```

它浓缩了任务类型、方法/基组、分子式、route、title、几何、能量（`HF=` 字段在 DFT 计算中同样表示 SCF 能量）、点群等信息。很多脚本和数据库就是解析这一段的。注意它每行固定宽度、会把数字截断换行，**用 grep 抓能量时不要抓 archive 段，抓 `SCF Done` 行更可靠**。

### 9.7 时间与结束

```text
 Job cpu time:       0 days  0 hours  0 minutes  1.6 seconds.
 Elapsed time:       0 days  0 hours  0 minutes  0.9 seconds.
 File lengths (MBytes):  RWF=      5 Int=      0 D2E=      0 Chk=      1 Scr=      1
 Normal termination of Gaussian 16 at Mon Nov 10 10:00:00 2026.
```

`Job cpu time` 是所有核的 CPU 时间之和，`Elapsed time` 是实际墙钟时间，两者之比大致反映并行效率。

### 9.8 常用 grep 速查

| 想看什么 | 命令 |
|---|---|
| 是否正常结束 | `tail -n 1 job.log`、`grep -c "Normal termination" job.log` |
| 所有 SCF 能量（优化每步一行） | `grep "SCF Done" job.log` |
| 最终 SCF 能量 | `awk '/SCF Done/{e=$5} END{print e}' job.log` |
| 基函数数目 | `grep -m1 "basis functions" job.log` |
| 电子数（核对电荷多重度） | `grep -m1 "alpha electrons" job.log` |
| 点群 | `grep "Full point group" job.log` |
| 优化收敛表 | `grep -A4 "Item               Value" job.log` |
| 是否找到驻点 | `grep "Stationary point found" job.log` |
| 频率（第6章） | `grep "Frequencies --" job.log` |
| 用时 | `grep -e "Job cpu time" -e "Elapsed time" job.log` |

批量处理多个 log 时，用 [`g16_log_summary.py`](../scripts/g16_log_summary.py) 会更方便（`python ../scripts/g16_log_summary.py -h`）。

---

## 10. 第一个几何优化：水分子 B3LYP/6-31G(d)

### 10.1 输入

文件 [`h2o_b3lyp_opt.gjf`](../examples/01-first-job-h2o/h2o_b3lyp_opt.gjf)：

```text
%chk=h2o_b3lyp_opt.chk
%mem=1GB
%nprocshared=2
#p B3LYP/6-31G(d) Opt

water geometry optimization, B3LYP

0 1
O     0.000000    0.000000    0.110000
H     0.000000    0.800000   -0.470000
H     0.000000   -0.800000   -0.470000

```

与单点相比，route 只换了方法并加了 `Opt`。初始结构故意给得粗糙一些（O–H 约 0.99 Å），这样能看到优化过程。

| 关键词 | 含义 | 为什么这样选 |
|---|---|---|
| `B3LYP` | 最常用的杂化泛函（20% HF 交换） | 入门首选，文献数据最多；正式研究建议加色散校正（第3章） |
| `6-31G(d)` | Pople 劈裂价层基组 + 重原子 d 极化函数 | 小而够用的入门基组；Gaussian 中 6-31G(d) 默认使用 6 个笛卡尔 d 函数（6D） |
| `Opt` | 几何优化到能量极小点 | 默认用冗余内坐标和 Berny 算法，最多步数按体系大小自动设置 |

### 10.2 运行与监控

```bash
g16 < h2o_b3lyp_opt.gjf > h2o_b3lyp_opt.log
grep "SCF Done" h2o_b3lyp_opt.log          # one line per optimization step
grep -c "SCF Done" h2o_b3lyp_opt.log       # number of steps so far
```

每个优化步打印一张收敛判据表（示意输出，数值仅供参考）：

```text
         Item               Value     Threshold  Converged?
 Maximum Force            0.000021     0.000450     YES
 RMS     Force            0.000015     0.000300     YES
 Maximum Displacement     0.000040     0.001800     YES
 RMS     Displacement     0.000032     0.001200     YES
 Predicted change in Energy=-1.068893D-09
 Optimization completed.
    -- Stationary point found.
```

四项都是 `YES`，或者力已远小于阈值（详见 [第5章](05-opt.md)），就会出现 `Optimization completed.` 和 `-- Stationary point found.`，随后打印优化后的内坐标（示意输出）：

```text
                           !   Optimized Parameters   !
                           ! (Angstroms and Degrees)  !
 ! Name  Definition              Value          Derivative Info.                !
 ! R1    R(1,2)                  0.9686         -DE/DX =    0.0                 !
 ! R2    R(1,3)                  0.9686         -DE/DX =    0.0                 !
 ! A1    A(2,1,3)              103.67           -DE/DX =    0.0                 !
```

### 10.3 结果合理吗

| 量 | 本计算（B3LYP/6-31G(d)，示意） | 实验值 |
|---|---|---|
| r(O–H) / Å | 约 0.97 | 0.9572 |
| ∠H–O–H / ° | 约 103.7 | 104.52 |

实验值取自 NIST CCCBDB（https://cccbdb.nist.gov/）。差异来自方法和基组本身，以及实验平衡结构与计算模型的差别；CCCBDB 也列出了各计算级别的参考值，可以用来核对自己的结果。

> [!IMPORTANT]
> `Stationary point found` 只说明找到了**驻点**（梯度为零），不能证明它是能量极小点。要确认是极小点，必须在**同一方法/基组**下做频率计算，确认没有虚频，这就是为什么实际工作中几乎总写 `Opt Freq`（[第6章](06-freq-thermo.md)）。

### 10.4 取出优化后的结构

- GaussView 打开 log，默认显示最后一帧，File → Save 可以另存为新的 gjf。
- 命令行可以取最后一个 Standard orientation 块：

```bash
# print the last "Standard orientation" block: atomic number, x, y, z (Angstrom)
awk '/Standard orientation:/{s=NR; n=0; next}
     s && NR>s+4 {if ($0 ~ /^ -+$/) {s=0; next}; n++;
                  line[n]=sprintf("%3d %12.6f %12.6f %12.6f", $2, $4, $5, $6)}
     END{for (i=1; i<=n; i++) print line[i]}' h2o_b3lyp_opt.log

# last SCF energy only (5th field of the "SCF Done" line)
awk '/SCF Done/{e=$5} END{print e}' h2o_b3lyp_opt.log
```

Standard orientation 块的结构是：标题行、虚线、两行表头、虚线、原子行、虚线，所以上面的 awk 从标题后第 5 行开始读，遇到虚线停止，每遇到新块就清零，最后剩下的就是最后一个结构。它输出的是原子序数而不是元素符号；需要标准 xyz 文件时，推荐用 [`log2xyz.py`](../scripts/log2xyz.py)（`python ../scripts/log2xyz.py -h`）。
- 下一步计算可以直接从 chk 读结构（`Geom=Check`），完全不必手抄坐标，见 [第2章](02-input-anatomy.md)。

---

## 11. Gaussian 文件类型一览

| 扩展名 | 类型 | 说明 |
|---|---|---|
| `.gjf` / `.com` | 输入（文本） | 两者完全等价；GaussView 默认用 .gjf，Linux 用户也常用 .com |
| `.log` / `.out` | 输出（文本） | Linux 默认 .log，G16W 常用 .out；内容格式相同 |
| `.chk` | 二进制检查点 | 保存几何、基组、轨道、力常数等，用于 `Geom=Check`、`Guess=Read`、`Opt=ReadFC`、作图；**依赖平台和版本**，不能跨机器/跨操作系统直接用 |
| `.fchk` / `.fch` | 格式化检查点（文本） | `formchk job.chk job.fchk` 生成；可跨平台，是 GaussView、Multiwfn、VMD 等的通用输入 |
| `.cub` / `.cube` | 格点数据（文本） | `cubegen` 生成的轨道、密度、静电势格点数据，用于作等值面图 |
| `.wfn` / `.wfx` | 波函数文件 | `Output=WFX` 等输出，供 Multiwfn、AIMAll 做 QTAIM 等分析 |
| `Gau-NNNNN.rwf/.int/.d2e/.skr/.inp` | scratch 临时文件 | 位于 `GAUSS_SCRDIR`，正常结束自动删除；异常退出会残留 |
| `.gbs` | 基组文件 | Gen 基组的常用扩展名（第3章） |

chk 与 fchk 的相互转换：

```bash
formchk job.chk job.fchk          # binary -> formatted (portable)
formchk -m=1gb big.chk big.fchk   # C.01+: give the utility more memory
unfchk job.fchk job.chk           # formatted -> binary (on the new machine)
chkchk job.chk                    # show route/title stored in a chk (and job status, B.01+)
```

chk 文件带版本号是 C.01 引入的。旧版本生成的 chk 可以用 `c8616` 升级，或用 `unfchk -fixver` 处理。更多工具（`cubegen`、`freqchk`、`newzmat`）见 [第16章](16-postprocess-multiwfn.md)。

---

## 12. 编辑器与文本格式：避开"看不见"的错误

### 12.1 推荐编辑器

- **VSCode** 加 Gaussian 输入文件语法高亮插件，例如 [mizu-bai/VSCode-Gaussian-Input-File-Support](https://github.com/mizu-bai/VSCode-Gaussian-Input-File-Support)。关键词有颜色，拼写错误（route 中的 QPErr）更容易发现。配合 Remote-SSH 插件可以直接编辑服务器上的文件。
- 服务器上用 `vim`/`nano` 即可。Windows 上不要用系统记事本保存 gjf（可能加 BOM 和 CRLF）。

### 12.2 换行符（CRLF vs LF）

Windows 文本的行尾是 `\r\n`（CRLF），Linux 是 `\n`（LF）。在 Windows 上编辑、再上传到 Linux 的 gjf，行尾多出的 `\r` 会让"空行"不再是空行，可能引起 `End of file in ZSymb`、关键词无法识别等奇怪错误。

```bash
file job.gjf                     # "... with CRLF line terminators" means trouble
dos2unix job.gjf                 # convert in place (if dos2unix is installed)
sed -i 's/\r$//' job.gjf         # same thing with sed
grep -c $'\r' job.gjf            # count lines that still contain CR
```

VSCode 右下角可以把 CRLF 切换为 LF。也可以在 VSCode 设置里把 `files.eol` 设为 `\n`。

### 12.3 编码与中文

- gjf 文件**只用 ASCII 字符**。title 行、注释（`!` 开头）也不要写中文，全角括号、全角空格、中文引号都可能导致解析错误。
- 文件名和目录名避免中文、空格和特殊符号（`( ) & #` 等）。Windows 下路径含中文或空格是启动失败（Link 1 not found / Error #2070）的最常见原因。
- title 段官方建议避免 `@ # !` 等特殊字符（`@` 在输入中表示包含外部文件，`!` 表示注释）。

### 12.4 推荐的项目目录组织

```text
project_x/
  00_struct/          # initial structures (xyz, mol, cif)
  01_opt_freq/        # B3LYP-D3(BJ)/def2-SVP opt+freq: *.gjf *.log *.chk
  02_sp_highlevel/    # high-level single points
  03_analysis/        # fchk, cube, Multiwfn outputs, figures
  notes.md            # what was run, which level, problems met
```

文件名中写明体系和任务，例如 `ts1_b3lypd3_def2svp_optfreq.gjf`。一个 gjf 对应一个 log 和一个 chk，三者同名，日后追溯最方便。

---

## 参数选择与收敛建议

入门阶段的默认参数建议：

| 参数 | 入门推荐值 | 说明 |
|---|---|---|
| `%mem` | 每核 0.5–2 GB，小分子 2–4 GB | G16 默认 800 MB；设得比物理内存大会严重变慢或报 `galloc` 错误；集群上取申请内存的 75–85% |
| `%nprocshared` | 小分子 2–8 核 | 不超过机器或作业申请的核数；小分子核数越多不一定越快（第17章） |
| `%chk` | 每个作业都写 | 便于续算、读结构、作图 |
| route 前缀 | `#p` | 输出更详细，便于诊断 |
| 方法/基组（练习用） | HF/6-31G(d) → B3LYP/6-31G(d) | 正式研究用带色散校正的组合，见第3章 |
| 任务 | 先 `Opt Freq` 再看结果 | 优化完必须确认无虚频 |
| `GAUSS_SCRDIR` | 本地盘，剩余空间 ≥ 20 GB | 大作业需要更多 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `g16: command not found` | 没有 source `g16.profile`，或 `g16root` 设错 | 检查 `~/.bashrc`；`g16root` 应为 `g16` 目录的**上一级** |
| `Permission denied` | 当前用户不在 Gaussian 所属组，或在 Windows 上解压后再上传导致权限丢失 | 联系管理员加组；重新在 Linux 上解压安装 |
| 作业一开始就失败，日志几乎为空，或提示无法打开 scratch 文件 | `GAUSS_SCRDIR` 指向的目录不存在或没有写权限 | `mkdir -p $GAUSS_SCRDIR`，检查权限 |
| `Erroneous write. Write -1 instead of ...` | scratch 磁盘或配额已满 | 清理 `Gau-*` 文件，换更大的 scratch，或用 `%rwf` 拆分（第17章） |
| `galloc: could not allocate memory.` | `%mem` 大于实际可用内存 | 减小 `%mem`，或向调度系统申请更多内存 |
| `End of file in ZSymb.`（l101.exe） | 缺空行（title 后或文件末尾）、缺电荷多重度行、CRLF 换行 | 补空行，`dos2unix` |
| `QPErr --- A syntax error was detected in the input line.`（l1.exe） | route 拼写错误、括号不配对、关键词不存在 | 看报错行下方的标记位置；用 `testrt "#p ..."` 检查 |
| `The combination of multiplicity 1 and 17 electrons is impossible.`（l301.exe） | 电子数奇偶与多重度不符 | 按第2章重新计算电荷和多重度 |
| Windows：`Link 1 not found`、`Severe Error Message #2070` | 路径含中文或空格、杀毒软件拦截、环境未配好、`%mem` 过大 | 改用纯英文路径，加白名单，减小 `%mem` |
| 跑完没有 Normal termination，也没有 Error termination | 被调度系统 kill（超时、超内存）或节点故障 | 查看调度系统的输出/错误文件；用 chk 续算（第5章、第17章） |
| GaussView 打开 log 只看到一个结构 | 作业是单点，或者只打开了 .chk | 优化作业用 log 打开，才能在 Results → Optimization 看轨迹 |

---

## 练习

1. 运行 `h2o_hf_sp.gjf`，记录 SCF 能量、HOMO/LUMO 能量和偶极矩；把基组改为 `6-311+G(d,p)` 再算一次，比较能量和基函数数目的变化。
2. 运行 `h2o_b3lyp_opt.gjf`，数一数优化走了多少步（`grep -c "SCF Done"`），把优化后的键长、键角与 CCCBDB 上的实验值和 B3LYP/6-31G* 计算值比较。
3. 故意删掉 `h2o_hf_sp.gjf` 最后的空行，或者把多重度改成 2，运行并观察报错出现在哪个 Link，以及 Error termination 上方的错误信息。
4. 在 GaussView 中搭建乙醇，Clean 后保存为 gjf，用文本编辑器打开，找出 `Geom=Connectivity` 和文件末尾的连接表。
5. 用 `formchk` 把 `h2o_b3lyp_opt.chk` 转成 fchk，在 GaussView 中打开并显示 HOMO（Results → Surfaces/Contours）。

---

## 参考资料

**官方文档**

- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)：修订版新功能、G09 到 G16 默认值变化、命令行与环境变量等价关系
- [Running Gaussian](https://gaussian.com/running/)：环境变量、scratch 文件、内存与并行
- [g16 Command Line Options](https://gaussian.com/options/)
- [Link 0 Commands](https://gaussian.com/link0/)
- [About Gaussian 16 Input](https://gaussian.com/input/)
- [Utility Programs](https://gaussian.com/utils/)、[testrt](https://gaussian.com/testrt/)、[formchk](https://gaussian.com/formchk/)
- [GaussView 6](https://gaussian.com/gaussview6/)、[GaussView 6.1.1 Release Notes](https://gaussian.com/gv611rn/)、[GaussView 6 Job Setup additions](https://gaussian.com/gv6jobset/)
- UNIX 二进制安装说明：https://gaussian.com/g16/g16bin_install.pdf
- [NIST CCCBDB](https://cccbdb.nist.gov/)：小分子实验与计算参考数据

**相关 GitHub 仓库与资源**（仅外链，更多见 [相关仓库汇总](../../Resources/github-repos.md)）

- [mizu-bai/VSCode-Gaussian-Input-File-Support](https://github.com/mizu-bai/VSCode-Gaussian-Input-File-Support)：VSCode 的 gjf 语法高亮插件
- [thisisntnathan/dftCourse](https://github.com/thisisntnathan/dftCourse)（[网站](https://thisisntnathan.github.io/dftCourse/)）：Cornell 的 Gaussian + Slurm 入门课程
- [ishikawa-group/gaussian-tutorial](https://github.com/ishikawa-group/gaussian-tutorial)：结构清晰的 Gaussian 讲义
- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian log 的标准库
- [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash)：g16 bash 工具集（getenergy、testroute 等）
- [liyuanhe211/Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)：Gaussian 报错大全（CC BY-NC-SA，仅外链参考）
- [计算化学公社 Gaussian 报错原帖](http://bbs.keinsci.com/thread-4829-1-1.html)（中文）

**书籍与论文**

- M. J. Frisch et al., *Gaussian 16, Revision C.01*, Gaussian, Inc., Wallingford CT, 2016（引用时把 Revision 换成你实际使用的版本）。
- J. B. Foresman, Æ. Frisch, *Exploring Chemistry with Electronic Structure Methods*, 3rd ed., Gaussian, Inc., 2015：官方入门教材。

---

[返回目录](../README.md) | [下一章 →](02-input-anatomy.md)
