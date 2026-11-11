# 第2章 输入文件结构详解（Gaussian Input File Anatomy: Link 0, Route, Molecule Specification and Link1）

本章逐段拆解 Gaussian 16 输入文件（.gjf/.com）：Link 0 命令（%chk、%oldchk、%mem、%nprocshared、%cpu）、route 关键词行（#p、Opt=(TS,CalcFC)、IOp）、标题、电荷与自旋多重度、笛卡尔坐标与 Z-matrix、冻结标记与片段标签（Fragment）、Geom=Connectivity 连接表、Geom=Check/AllCheck 与 Guess=Read 读取 chk、--Link1-- 多步作业，以及 NoSymm 和坐标重定向。读完本章，你能手写或修改任何常规 Gaussian 输入，并看懂 GaussView 生成的输入里每一行在做什么。

---

## 1. 本章目标与适用场景

| 你将学会 | 对应小节 |
|---|---|
| 输入文件各段的顺序与空行规则 | §2 |
| Link 0 命令的写法、单位和优先级 | §3 |
| route 行的写法：`#`/`#p`/`#t`、关键词选项、IOp | §4 |
| 正确设置电荷与自旋多重度，并做电子数奇偶校验 | §6 |
| 笛卡尔坐标、Z-matrix、混合格式；原子参数（Fragment、Iso）和冻结标记 | §7 |
| 坐标之后的附加输入段（连接表、ModRedundant、Gen 基组等）的顺序 | §8 |
| Geom=Connectivity 连接表的含义，以及删除时的注意事项 | §9 |
| 从 chk 读结构和波函数（Geom=Check、Geom=AllCheck、Guess=Read、%oldchk） | §10 |
| 用 --Link1-- 写"优化+频率 → 高精度单点"的多步作业 | §11 |
| 输入单位、对称性、标准取向与 NoSymm | §12 |

本章示例输入都在 [examples/02-input-anatomy](../examples/02-input-anatomy/README.md)。

---

## 2. 原理速览：输入文件的整体结构

### 2.1 各段顺序

Gaussian 按"段"读取输入，**段与段之间用空行分隔**。完整顺序如下：

```text
%chk=job.chk                 <- Link 0 commands (lines starting with %), NO blank line after
%mem=8GB
%nprocshared=16
#p B3LYP/6-31G(d) Opt Freq   <- route section (one or more lines starting with #)
                             <- blank line
title line(s), at most 5     <- title section
                             <- blank line
0 1                          <- charge and spin multiplicity
C   0.000  0.000  0.000      <- molecule specification
H   0.000  0.000  1.090
...
                             <- blank line
additional input sections    <- e.g. connectivity, ModRedundant, Gen basis, ECP ...
                             <- each one ends with a blank line
                             <- the FILE must end with a blank line
```

| 段 | 是否必需 | 结束标志 | 说明 |
|---|---|---|---|
| Link 0（`%` 行） | 否 | **不需要空行**，直接接 route | 文件名、内存、核数等资源设置 |
| route（`#` 行） | 是 | 空行 | 方法、基组、任务类型和各种选项 |
| title | 是（`Geom=AllCheck` 时省略） | 空行 | 最多 5 行，只做注释，不参与计算 |
| 电荷/多重度 + 分子坐标 | 是（`Geom=AllCheck` 时省略；`Geom=Check` 时只写电荷多重度行） | 空行 | 见 §6、§7 |
| 附加输入段 | 视关键词而定 | 通常是空行 | 顺序固定，见 §8 |

> [!IMPORTANT]
> 两条最重要的规则：① Link 0 行与 route 行之间**不能**有空行；② 文件**最后必须**有空行。违反规则时，L101 常报 `End of file in ZSymb.`，或者把下一段内容当成本段来读，出现莫名其妙的错误。

### 2.2 通用语法规则

根据官方 *About Gaussian 16 Input* 文档：

- 输入**不区分大小写**，格式自由。`B3LYP`、`b3lyp`、`B3lyp` 都一样。
- 同一行内的各项可以用空格、制表符、逗号或斜杠分隔，多个空格视为一个分隔符。所以 `0 1` 和 `0,1` 都可以，`O 0.0 0.0 0.0` 和 `O,0.0,0.0,0.0` 也都可以。
- 关键词的选项有四种写法，等价：`Opt=TS`、`Opt(TS)`、`Opt=(TS,CalcFC)`、`Opt(TS,CalcFC)`。多个选项必须放在括号里。
- 关键词和选项可以缩写到**在整个 Gaussian 系统中唯一**的最短形式（例如 `SCF=Conventional` 可以写成 `SCF=Conven`，但不能写成 `Conv`，因为还有 `Convergence`）。为了可读性和可移植性，**不建议缩写**。
- `!` 开始注释，可以出现在行内任何位置，也可以单独成行。
- `@filename` 把外部文件的内容插入到当前位置，写成 `@filename/N` 时不在输出开头回显该文件内容。常用于引用基组文件（第3章）。

---

## 3. Link 0 命令

### 3.1 常用命令一览

| 命令 | 示例 | 作用与说明 |
|---|---|---|
| `%Chk` | `%chk=ts1.chk` | 检查点文件的名称和位置，可带路径。**几乎每个作业都应写**，便于续算、读结构、作图 |
| `%OldChk` | `%oldchk=opt.chk` | 作业步开始时，把这个 chk 的内容**复制**到 `%chk` 指定的新 chk 中，旧文件保持不变。推荐用它从上一个作业读取数据 |
| `%Mem` | `%mem=8GB` | 动态内存上限。单位可以是 KB、MB、GB、TB 或 KW、MW、GW、TW（W = 8 字节的 word），数字和单位之间不留空格。不写单位时按 word 计。默认 800 MB（100 MW） |
| `%NProcShared` | `%nprocshared=16` | 共享内存并行使用的核数，由操作系统选择具体是哪些核 |
| `%CPU` | `%cpu=0-15` | 指定具体的核编号列表，可用逗号和范围，例如 `0-7,16-23`；官方推荐它取代 `%NProcShared` |
| `%GPUCPU` | `%gpucpu=0-1=0,1` | GPU 列表 `=` 控制核列表，仅 Linux + NVIDIA GPU；控制核也要包含在 `%CPU` 中（第17章） |
| `%LindaWorkers` | `%lindaworkers=node1,node2` | Linda 跨节点并行（另需 Linda 许可） |
| `%RWF` | `%rwf=/scratch/a/,50GB,/scratch/b/,-1` | 指定（或拆分到多个磁盘）读写文件 .rwf；只写目录时必须以 `/` 结尾，`-1` 表示不限大小 |
| `%NoSave` | `%nosave` | 正常结束后删除**在它之前**命名的 scratch 文件；在它之后命名的文件会保留 |
| `%Save` | `%save` | 结束时保留 scratch 文件（调试用） |
| `%KJob` | `%kjob l502 2` | 第 2 次运行完 Link 502 后停止（调试用） |

`%NoSave` 的典型写法：保留 chk、删除 rwf。

```text
%rwf=/scratch/user/water
%nosave
%chk=water
```

`%rwf` 在 `%nosave` 之前，作业成功后会被删除；`%chk` 在之后，会被保留。作业**异常结束时所有文件都会保留**，方便续算。

### 3.2 %nprocshared 与 %cpu 怎么选

| 情形 | 推荐 |
|---|---|
| 独占整个节点或自己的工作站 | `%cpu=0-N` 或 `%nprocshared=N`，都可以 |
| 集群共享节点（多个作业共用一个节点） | `%nprocshared=N`，N 等于申请的核数。`%cpu=0-15` 会把线程绑到 0–15 号核，而调度系统分给你的可能是 16–31 号核，结果反而更慢甚至出错 |
| 用脚本统一管理 | 输入中不写，在作业脚本中设置环境变量 `GAUSS_CDEF`（或命令行 `-c`）；见 [第17章](17-hpc-jobs.md) |

> [!WARNING]
> `%nprocshared` 不要大于作业实际申请到的核数，否则会严重超订（oversubscription），速度可能慢几倍。例如 Slurm 申请 `--cpus-per-task=16`，输入里就写 `%nprocshared=16`，不要写 64。`%mem` 取申请内存的 75–85%，给程序本身留出余量。

### 3.3 Link 0、命令行、环境变量与 Default.Route

同一个设置有四种来源，优先级从高到低（官方 *Link 0 Commands → Equivalencies*）：

| 优先级 | 方式 | 示例 |
|---|---|---|
| 1（最高） | 输入文件中的 Link 0 行 | `%mem=8GB`、`%cpu=0-7` |
| 2 | g16 命令行选项（B.01 起） | `g16 -m=8GB -c="0-7" job.gjf` |
| 3 | 环境变量 | `export GAUSS_MDEF=8GB`、`export GAUSS_CDEF="0-7"` |
| 4 | `Default.Route` 文件（当前目录或 `$GAUSS_EXEDIR`） | `-M- 8GB`、`-C- 0-7` |

对应关系：`%Mem` ↔ `-m` ↔ `GAUSS_MDEF` ↔ `-M-`；`%CPU` ↔ `-c` ↔ `GAUSS_CDEF` ↔ `-C-`；`%GPUCPU` ↔ `-g` ↔ `GAUSS_GDEF` ↔ `-G-`；`%LindaWorkers` ↔ `-w` ↔ `GAUSS_WDEF` ↔ `-W-`。`%NProcShared` 对应的 `-p`/`GAUSS_PDEF`/`-P-` 已被官方标为 deprecated。还有 `-r`/`GAUSS_RDEF` 可以给所有作业追加默认 route 关键词。

---

## 4. Route 段

### 4.1 `#`、`#p`、`#t`

route 段以 `#` 开头，可以写多行，以空行结束。

| 写法 | 输出详细程度 | 建议 |
|---|---|---|
| `#` 或 `#N` | 正常输出 | 可以 |
| `#P` | 额外输出：每轮 SCF 迭代、各 Link 的进出与计时等 | **推荐**。诊断 SCF 不收敛等问题时非常有用，文件也大不了多少 |
| `#T` | 精简输出（terse） | 不推荐，排查问题时信息不够 |

### 4.2 关键词、选项与"模型化学"

route 的核心是"方法/基组"（Gaussian 称为 model chemistry）加任务类型：

```text
#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt=(MaxCycles=200) Freq SCRF=(SMD,Solvent=Water)
   |____________|  |_______________________| |__________________| |__| |______________________|
   method/basis      dispersion correction    job type + options  job   implicit solvent
```

- **方法/基组**用斜杠连接，例如 `B3LYP/6-31G(d)`、`M062X/def2TZVP`、`HF/STO-3G`。也可以分开写成 `B3LYP 6-31G(d)`，但斜杠写法最清楚。不写基组时默认用 STO-3G。
- **任务类型**：`SP`（单点，默认）、`Opt`、`Freq`、`IRC`、`Scan`、`Stable`、`NMR`、`Polar` 等，可以组合，例如 `Opt Freq`。Gaussian 会把 `Opt Freq` 自动拆成两个作业步（优化完成后自动做频率），所以 log 中会有两行 Normal termination。
- 方法名前可以加 `R`（限制性）、`U`（非限制性）、`RO`（限制性开壳层），例如 `UB3LYP`。不加时，单重态默认用 R，其他多重度默认用 U（第13章）。

常见 route 示例：

| route | 含义 |
|---|---|
| `#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt Freq` | 有机分子优化 + 频率的常用写法 |
| `#p M062X/def2TZVP EmpiricalDispersion=GD3 SCRF=(SMD,Solvent=Water)` | 溶液中的高精度单点 |
| `#p wB97XD/def2SVP Opt=(TS,CalcFC,NoEigenTest) Freq` | 过渡态优化 + 频率（第7章） |
| `#p B3LYP/GenECP EmpiricalDispersion=GD3BJ Opt` | 含重元素的混合基组（第3章） |
| `#p TD=(NStates=20) CAM-B3LYP/def2TZVP` | 前 20 个激发态（第10章） |
| `#p B3LYP/def2SVP Opt=ModRedundant` | 带约束的优化，坐标后要写 ModRedundant 段（第5章、第8章） |
| `#p wB97XD/def2TZVP Geom=Check Guess=Read` | 从 chk 读结构和初猜做单点（§10） |

> [!TIP]
> 提交前用 `testrt "#p B3LYP/6-31G(d) Opt=(TS,CalcFC) Freq"` 检查语法。有语法错误时它会打印 `QPErr` 和出错位置（只报第一个错误）。注意 `testrt` 只检查语法，不检查关键词组合在化学上是否合理。

### 4.3 IOp：内部选项

`IOp(overlay/option=value)` 直接设置 Gaussian 内部的"非标准 route"选项，多个用逗号分隔：

```text
#p TD=(NStates=10) CAM-B3LYP/def2SVP IOp(9/40=4)
#p HF/6-31G(d) Pop=MK IOp(6/33=2,6/41=10,6/42=17)
```

第一行让 TD 输出打印更多的跃迁组态系数（第10章），第二行输出 RESP 拟合所需的静电势格点（第12章）。官方 DFT 文档还用 `IOp(3/76=...)`、`IOp(3/77=...)`、`IOp(3/78=...)` 自定义杂化泛函中各部分的比例。

官方 IOp 文档的几条规则：

- `IOp(ov/op=val)` 把 overlay `ov` 的第 `op` 号选项设为 `val`，并且作用于该 overlay 在这个作业步中的**每一次**出现；
- 显式写在 route 中的 IOp **不会传递**到自动生成的后续作业步。例如 `Opt Freq` 中的频率步收不到你写的 IOp；想让优化和频率都用某个非默认格点，必须用 `Int=Grid=...` 这样的正式关键词，而不能用 IOp；
- 使用 IOp 后，log 末尾**不再生成 archive 段**。

> [!CAUTION]
> IOp 的风险：① 选项含义只在 IOps Reference 和少数关键词页面中说明，可能随版本变化；② 它不在关键词合法性检查范围内，写错往往**不报错**，只是悄悄得到错误结果；③ 在 `Opt Freq` 等组合任务中只对第一步生效，容易造成前后两步设置不一致。只使用官方文档或可靠文献中**明确给出**的 IOp，在论文 SI 中写明完整 route；能用正式关键词实现的功能就不要用 IOp。

---

## 5. Title 段

- 至少一行，最多 5 行，以空行结束。
- Gaussian 不解释其内容，只原样写入输出和 archive 段。
- 写清"体系 + 任务 + 条件"，例如 `TS1 for SN2 Cl- + CH3Br, SMD water`。
- 官方建议避免使用 `@`、`#`、`!` 和反斜杠等特殊字符（`@` 表示包含文件，`!` 表示注释）。只用 ASCII 字符，不要写中文。
- 用 `Geom=AllCheck` 时 title 从 chk 中读取，输入中**不能再写** title（§10）。

---

## 6. 电荷与自旋多重度

### 6.1 定义

分子坐标段的第一行是**净电荷**和**自旋多重度**：

$$
2S+1,\qquad S=\frac{N_\alpha-N_\beta}{2}
$$

电子总数为

$$
N_e=\sum_A Z_A-q
$$

其中 $Z_A$ 是原子序数，$q$ 是净电荷（阳离子为正）。**奇偶校验**：$N_e$ 为偶数时，多重度只能是奇数（1、3、5…）；$N_e$ 为奇数时，多重度只能是偶数（2、4、6…）。不满足时 L301 报错：

```text
 The combination of multiplicity 1 and 17 electrons is impossible.
```

使用赝势（ECP）不影响奇偶性，因为 ECP 替换掉的芯电子数总是偶数。

### 6.2 示例

| 体系 | $\sum Z$ | 电荷 | $N_e$ | 多重度 | 电荷/多重度行 |
|---|---|---|---|---|---|
| H2O | 10 | 0 | 10 | 1（单重态） | `0 1` |
| OH· 自由基 | 9 | 0 | 9 | 2（二重态） | `0 2` |
| OH⁻ | 9 | −1 | 10 | 1 | `-1 1` |
| NH4⁺ | 11 | +1 | 10 | 1 | `1 1` |
| O2 基态 | 16 | 0 | 16 | 3（三重态，$^3\Sigma_g^-$） | `0 3` |
| O2 单重态（激发态） | 16 | 0 | 16 | 1 | `0 1` |
| CH2 卡宾 | 8 | 0 | 8 | 3（基态）或 1 | `0 3` / `0 1` |
| 苯自由基阴离子 | 42 | −1 | 43 | 2 | `-1 2` |
| [Fe(H2O)6]³⁺ 高自旋（d⁵） | 86 | +3 | 83 | 6 | `3 6` |
| [Cu(NH3)4]²⁺（d⁹） | 69 | +2 | 67 | 2 | `2 2` |

> [!TIP]
> 过渡金属配合物的多重度往往不止一种可能（高自旋/低自旋）。应计算几种可能的多重度并比较能量，同时检查 $\langle S^2\rangle$ 和波函数稳定性，详见 [第13章](13-openshell-stability.md)。

### 6.3 片段的电荷与多重度

Counterpoise（BSSE 校正，[第14章](14-weak-bsse.md)）和 `Guess=Fragment`（碎片初猜/broken-symmetry，[第13章](13-openshell-stability.md)）需要为每个片段分别指定电荷与多重度，格式为：

```text
total_charge total_mult  frag1_charge frag1_mult  frag2_charge frag2_mult ...
```

例如两个中性单重态分子组成的复合物写 `0 1 0 1 0 1`。原子用 `(Fragment=n)` 标注归属，片段编号必须从 1 开始连续编号。`Guess=Fragment` 中负的多重度表示该片段的未成对电子取 β 自旋。ONIOM 的电荷多重度行有自己的格式（按层给出多对），见 [第15章](15-oniom.md)。

### 6.4 用脚本做奇偶校验

示例目录中的 [`check_charge_mult.py`](../examples/02-input-anatomy/check_charge_mult.py) 读取 gjf 的电荷/多重度行和原子列表，计算电子数并检查奇偶性（支持笛卡尔坐标、Z-matrix、`C1`、`O(Fragment=1)`、`C-CT-0.1` 这类标签，以及 `Bq` 鬼原子）。核心逻辑只有几行：

```python
mol = sections(path)[2]                        # molecule specification block
nums = [int(x) for x in re.split(r"[\s,]+", mol[0].strip())]
charge, mult = nums[0], nums[1]                # total charge / multiplicity
nelec = -charge
for line in mol[1:]:
    f = line.split()
    if f[0].lower().rstrip(":") in ("variables", "constants"):
        break                                  # Z-matrix variable section
    if len(f) == 2 or ("=" in f[0] and "(" not in f[0]):
        continue                               # "R1 0.97" or "R1=0.97"
    nelec += atomic_number(f[0].split(",")[0])
ok = nelec % 2 != mult % 2                     # even N_e <-> odd multiplicity
```

```bash
python check_charge_mult.py *.gjf
```

仓库的输入生成器 [`gjf_builder.py`](../scripts/gjf_builder.py) 生成输入时也会做同样的检查（`python ../scripts/gjf_builder.py -h`）。

---

## 7. 分子坐标

### 7.1 笛卡尔坐标

最常用、最直观的格式：`元素 x y z`，默认单位 Å。

```text
0 1
C     1.168181   -0.400382    0.000000
C     0.000000    0.559462    0.000000
O    -1.190083   -0.227669    0.000000
H    -1.946623    0.381525    0.000000
...
```

完整示例见 [`ethanol_cartesian.gjf`](../examples/02-input-anatomy/ethanol_cartesian.gjf)。元素可以写符号，也可以写原子序数（`6` 等于 `C`）；符号后可以跟编号作为标签（`C1`、`C2`，标签总长不超过 4 个字符）。

### 7.2 Z-matrix（内坐标）

每行用键长、键角、二面角相对已定义的原子来确定新原子的位置：

```text
Element  atom1 bond-length  atom2 bond-angle  atom3 dihedral
```

过氧化氢示例（[`h2o2_zmatrix.gjf`](../examples/02-input-anatomy/h2o2_zmatrix.gjf)）：

```text
0 1
H
O 1 R1
O 2 R2 1 A1
H 3 R1 2 A1 1 D1
  Variables:
R1 0.97
R2 1.45
A1 100.0
D1 115.0
```

| 行 | 含义 |
|---|---|
| `H` | 第 1 个原子，放在原点 |
| `O 1 R1` | 第 2 个原子，与原子 1 相距 R1 |
| `O 2 R2 1 A1` | 第 3 个原子，与原子 2 相距 R2，键角 ∠(3-2-1) = A1 |
| `H 3 R1 2 A1 1 D1` | 第 4 个原子，与原子 3 相距 R1，键角 ∠(4-3-2) = A1，二面角 (4-3-2-1) = D1 |
| `Variables:` 段 | 给变量赋值；两个 O–H 共用变量 R1，表示二者相等 |

要点：

- 键角必须在 0° 到 180° 之间（不能正好等于 180°），二面角可以取任意值。直线型片段（如 C≡C–H）要借助哑原子 `X` 处理。
- 原子引用既可以用序号（从 1 开始，电荷多重度行算第 0 行），也可以用标签（如 `O1`）。
- 变量段之后还可以写 `Constants:` 段，其中的变量在优化中保持固定。
- 默认的 `Opt` 会把 Z-matrix 转换成冗余内坐标来优化；`Opt=Z-matrix` 才严格按 Z-matrix 变量优化。现在做约束优化更推荐 `Opt=ModRedundant`（[第5章](05-opt.md)）。
- 刚性扫描（`Scan` 关键词）正是在 Z-matrix 变量上定义的（[第8章](08-scan.md)）。

哑原子示例（固定 NH3 的三重轴）：

```text
0 1
N
X 1 1.0
H 1 NH 2 HNX
H 1 NH 2 HNX 3  120.0
H 1 NH 2 HNX 3 -120.0
  Variables:
NH 1.01
HNX 68.0
```

`X` 不计入原子编号（ModRedundant 等处的原子序号不包括哑原子），也不带电子。

### 7.3 混合格式

笛卡尔行和 Z-matrix 行可以出现在同一个分子说明中，例如先用笛卡尔坐标给出骨架，再用内坐标挂上 H：

```text
0 1
C     0.000000    0.000000    0.000000
O     0.000000    0.000000    1.205000
H 1 1.10 2 121.8
H 1 1.10 2 121.8 3 180.0
```

### 7.4 原子附加参数

在元素符号后用括号添加参数，多个参数用逗号分隔：

| 参数 | 示例 | 用途 |
|---|---|---|
| `Fragment=n` | `O(Fragment=1)` | 片段归属，用于 Counterpoise、`Guess=Fragment`（§6.3） |
| `Iso=n` | `H(Iso=2)` | 同位素质量（整数时取该同位素的精确质量），影响频率和热化学 |
| `Spin=n` | `C(Iso=13,Spin=3)` | 核自旋（单位 1/2），用于 NMR/ESR 相关性质 |
| `ZNuc=n`、`QMom=n` 等 | — | 修改核电荷、核四极矩等，少用 |

示例：`C(Iso=13,Fragment=2) 0.00 0.00 0.00`。

分子力学和 ONIOM 计算还会在元素后加原子类型和部分电荷，例如 `C-CT-0.32`；`O-Bq` 表示 O 的**鬼原子**（有基函数和格点，但没有核电荷和电子）。详见 [第14章](14-weak-bsse.md)、[第15章](15-oniom.md)。

### 7.5 冻结标记

笛卡尔坐标的元素后可以加一个整数**冻结标记**：

```text
0 1
C    0   0.000000   0.000000   0.000000    ! 0 (or omitted): optimized
C   -1   1.540000   0.000000   0.000000    ! -1: frozen in Opt
...
```

| 标记 | 含义 |
|---|---|
| 省略或 `0` | 正常优化 |
| `-1` | 该原子在几何优化中固定不动 |

冻结标记只对 `Opt` 生效，单点计算会忽略它。表面簇模型固定边缘原子、蛋白活性位点截取模型固定骨架原子时常用。更灵活的约束（固定某个键长、键角、二面角）用 `Opt=ModRedundant`，见 [第5章](05-opt.md)。GaussView 中可以用 Edit → Atom List 批量设置冻结。

### 7.6 片段标签完整示例

[`fragment_labels_example.gjf`](../examples/02-input-anatomy/fragment_labels_example.gjf) 是水–氨氢键复合物的 Counterpoise 单点：

```text
%chk=h2o_nh3_cp.chk
%mem=2GB
%nprocshared=4
#p B3LYP/6-31+G(d,p) EmpiricalDispersion=GD3BJ Counterpoise=2

water ammonia complex, counterpoise corrected single point

0 1 0 1 0 1
O(Fragment=1)     0.000000    0.000000   -2.980000
H(Fragment=1)     0.000000    0.000000   -2.010000
H(Fragment=1)     0.929422    0.000000   -3.220365
N(Fragment=2)     0.000000    0.000000    0.000000
H(Fragment=2)     0.940000    0.000000    0.380000
H(Fragment=2)    -0.470000    0.814064    0.380000
H(Fragment=2)    -0.470000   -0.814064    0.380000

```

`Counterpoise=2` 表示 2 个片段；电荷多重度行依次是整体、片段 1、片段 2。每个原子都必须标注片段，漏标会导致报错或得到错误的 BSSE。输出的读法和结合能的计算见 [第14章](14-weak-bsse.md)。

---

## 8. 附加输入段及其顺序

很多关键词需要在分子坐标之后再读一段输入。多个附加段同时存在时，**顺序是固定的**（摘自官方 *Section Ordering* 表，只列常用部分）：

| 顺序 | 附加段 | 触发关键词 | 末尾空行 |
|---|---|---|---|
| 1 | 连接表 | `Geom=Connectivity` / `Geom=ModConnectivity` | 是 |
| 2 | 冻结原子修改 | `Geom=ReadOpt` | 是 |
| 3 | 内坐标修改（约束、扫描） | `Opt=ModRedundant` | 是 |
| 4 | 第 2/第 3 个分子（title + 坐标） | `Opt=QST2` / `Opt=QST3` | 是 |
| 5 | 同位素质量 | `Freq=ReadIsotopes` | 是 |
| 6 | 基组 | `Gen`、`GenECP`、`ExtraBasis` | 是 |
| 7 | ECP | `Pseudo=Read`、`GenECP` | 是 |
| 8 | PCM 额外参数 | `SCRF=Read` | 是 |
| 9 | NBO 输入 | `Pop=NBORead` 等 | — |

例如同时使用 `Geom=Connectivity`、`Opt=ModRedundant`、`GenECP` 时：

```text
0 1
... coordinates ...
                        <- blank
 1 2 1.0 3 1.0          <- connectivity table
 2
 3
                        <- blank
B 1 2 F                 <- ModRedundant input
                        <- blank
C H O 0                 <- Gen basis section
6-31G(d)
****
Pd 0
LANL2DZ
****
                        <- blank
Pd 0                    <- ECP section
LANL2DZ
                        <- blank (end of file)
```

---

## 9. Geom=Connectivity 与连接表

### 9.1 格式

GaussView 保存 gjf 时默认在 route 中加 `Geom=Connectivity`，并在坐标后写连接表：

```text
 1 2 2.0 3 1.0 5 1.0
 2
 3 4 1.0
 4
 5 6 1.0 7 1.0 8 1.0
 6
 7
 8
```

每行格式为 `M N1 Order1 N2 Order2 ...`：原子 M 与 N1 以键级 Order1 成键，依此类推。每个原子占一行（即使没有新增的键也要写出编号），每个键只在编号较小的原子那一行出现一次。键级 1.0、1.5（芳香）、2.0、3.0。上例是乙酸（C1=O2，C1–O3，O3–H4，C1–C5，C5–H6/H7/H8），完整输入见 [`geom_connectivity_example.gjf`](../examples/02-input-anatomy/geom_connectivity_example.gjf)。

### 9.2 什么时候需要、什么时候可以删

| 计算 | 连接表的作用 |
|---|---|
| 普通 DFT/HF 的单点、优化、频率 | 基本不起作用（冗余内坐标会根据几何自动生成），可以删 |
| 分子力学（UFF、Amber）和 ONIOM 的 MM 部分 | **必需**，决定成键项和原子类型 |
| 难以自动识别成键的体系（金属配位、弱配位） | 可以用它或 `Geom=ModConnectivity` 显式补充键，影响优化所用的内坐标 |

> [!WARNING]
> 删除连接表时必须**同时**删掉 route 中的 `Geom=Connectivity`，反之亦然。只删其中之一，Gaussian 会把下一段（如 Gen 基组或 ModRedundant）当成连接表来读，或者等待一段不存在的输入，从而报错。增删原子后，连接表的行数和编号也必须同步修改。

`Geom=ModConnectivity` 只写需要修改的行，键级写 `-1` 表示删除这根键：`8 4 1.0 5 2.0 9 -1`。

---

## 10. 从检查点文件读取：Geom=Check、Geom=AllCheck、Guess=Read、%OldChk

多步计算时，不要从 log 中手抄坐标。直接从 chk 读取不会损失精度，也不容易出错。

| 写法 | 从 chk 读什么 | 输入中还要写什么 |
|---|---|---|
| `Geom=Check` | 几何（包括 Z-matrix 变量） | title、电荷/多重度行（**不写坐标**） |
| `Geom=AllCheck` | 几何、电荷/多重度、title | **只写 route**（以及关键词需要的附加段），不写 title 和电荷多重度 |
| `Guess=Read` | 分子轨道（作为 SCF 初猜） | 与 Geom 选项独立；换了基组或方法时，Gaussian 会把旧轨道投影到新基组上 |
| `Opt=ReadFC` / `IRC=ReadFC` | 力常数（Hessian） | 第5章、第7章 |
| `Geom=(Check,Step=N)` | 失败或中断的优化中第 N 步的结构 | 从中间步重启优化 |

`%chk` 与 `%oldchk` 的配合：

```text
%oldchk=step1_opt.chk        <- read from here (left unchanged)
%chk=step2_sp_tz.chk         <- copy of it is made here, results are written here
#p wB97XD/def2TZVP Geom=AllCheck Guess=Read

```

- 只写 `%chk=step1_opt.chk` 也能读，但新计算的结果会**覆盖**这个 chk。如果新计算失败，旧的优化结果和波函数可能就没了。
- 用 `%oldchk` 读、`%chk` 写，旧文件保持不变，推荐这样做。
- chk 与 Gaussian 版本和平台有关。跨机器使用要先在原机器上 `formchk` 生成 fchk，再在新机器上 `unfchk` 转回 chk。C.01 起 chk 带版本号，旧 chk 可用 `c8616` 升级或 `unfchk -fixver` 处理（[第16章](16-postprocess-multiwfn.md)）。

> [!NOTE]
> `Geom=AllCheck` 默认沿用 chk 中的输入取向和连接信息。C.01 起，可以用 `Geom=(AllCheck,UseStandardOrientation)` 改用 chk 中的标准取向，用 `Geom=(AllCheck,GenConnectivity)` 重新生成连接关系。

---

## 11. --Link1-- 多步作业

### 11.1 规则

一个输入文件中可以写多个作业步，用一行 `--Link1--` 分隔：

- `--Link1--` 的**前面必须有一个空行**（即上一步的输入已经完整结束）。
- 每一步都是完整的输入：Link 0 → route → title → 分子说明 → 附加段。
- Link 0 设置按作业步生效。稳妥的做法是**每一步都写全** `%chk`、`%mem`、`%nprocshared`，而且各步的 `%chk` 名称要一致（或用 `%oldchk` 显式指定读取哪个 chk），否则 `Geom=Check` 找不到上一步的结果。
- 前一步失败时，后面的步骤也不会执行。

### 11.2 示例：优化 + 频率，再做高精度单点

[`link1_opt_then_sp.gjf`](../examples/02-input-anatomy/link1_opt_then_sp.gjf)：

```text
%chk=h2co_link1.chk
%mem=2GB
%nprocshared=4
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

formaldehyde opt and freq, step 1

0 1
O     0.000000    0.000000    0.683501
C     0.000000    0.000000   -0.536614
H     0.000000    0.934390   -1.124164
H     0.000000   -0.934390   -1.124164

--Link1--
%chk=h2co_link1.chk
%mem=2GB
%nprocshared=4
#p wB97XD/def2TZVP Geom=Check Guess=Read SP

formaldehyde single point at the optimized geometry, step 2

0 1

```

逐段解释：

| 部分 | 说明 |
|---|---|
| 第 1 步 route | B3LYP-D3(BJ)/def2-SVP 优化 + 频率。`Opt Freq` 本身就会在内部拆成两步，所以这个文件实际运行 3 个作业步，log 中会有 3 行 Normal termination |
| 坐标后的空行 | 结束第 1 步的分子说明；同时满足 "`--Link1--` 前必须有空行" |
| 第 2 步 `%chk` | 与第 1 步同名，才能读到优化后的结构 |
| `Geom=Check` | 从 chk 读优化后的几何；输入中仍需 title 和 `0 1` 行 |
| `Guess=Read` | 读第 1 步的轨道作初猜（会自动从 def2-SVP 投影到 def2-TZVP） |
| `wB97XD/def2TZVP` | 在更高水平上做单点，写作 ωB97X-D/def2-TZVP//B3LYP-D3(BJ)/def2-SVP（第3章） |
| 末尾空行 | `0 1` 行后的空行结束分子说明段，同时也是文件末尾空行 |

检查：

```bash
grep -c "Normal termination" link1_opt_then_sp.log   # expect 3
grep "SCF Done" link1_opt_then_sp.log | tail -n 1    # the wB97XD energy
```

用 `Geom=AllCheck` 写第 2 步可以更短（不要 title 和电荷多重度行）：

```text
--Link1--
%chk=h2co_link1.chk
%mem=2GB
%nprocshared=4
#p wB97XD/def2TZVP Geom=AllCheck Guess=Read SP

```

> [!TIP]
> `Opt Freq` 合写与分两步（先 Opt，再用 `Geom=AllCheck Guess=Read Freq` 读 chk）是等价的，只要两步的方法、基组、溶剂模型、积分格点完全一致。频率计算必须与优化在同一水平下进行，否则驻点判断和热化学都没有意义（[第6章](06-freq-thermo.md)）。

---

## 12. 单位、对称性与坐标重定向

### 12.1 输入单位

默认长度单位是 Å，角度单位是度。`Units` 关键词可以修改：

| 关键词 | 含义 |
|---|---|
| `Units=Ang`（默认） | 长度 Å |
| `Units=AU` | 长度 Bohr（原子单位，1 Bohr = 0.529177 Å） |
| `Units=Rad` | 角度弧度 |

```text
#p HF/6-31G(d) Units=AU

water in bohr

0 1
O     0.000000    0.000000    0.225372
H     0.000000    1.442312   -0.901488
H     0.000000   -1.442312   -0.901488

```

把 Bohr 坐标当成 Å 读入（或反之）是 L202 `Small interatomic distances encountered` 等错误的常见原因：例如 Å 坐标被当作 Bohr 读入，所有键长都会缩短约一半。从其他程序（如 VASP 中以分数坐标表示的结构、某些以 Bohr 输出的 xyz）拿结构时要特别注意单位。

### 12.2 标准取向与 NoSymm

Gaussian 默认会识别分子点群，并把分子平移、旋转到**标准取向**（Standard orientation），以便利用对称性节省计算量。因此：

- 输出中的轨道、偶极矩分量、多极矩、cube 文件都在标准取向坐标系中；
- 优化每一步打印的 `Input orientation` 和 `Standard orientation` 可能不同。

`NoSymm`（`NoSymmetry` 的常用简写）的作用是**不重定向分子**，所有计算都在输入取向下进行（程序仍会识别点群，但不利用对称性）。`Symmetry=None` 则完全跳过点群识别与对称性处理。

| 场景 | 是否加 NoSymm |
|---|---|
| 常规优化、频率 | 一般不加，对称性能加速 |
| 需要与外部坐标系对应（加外电场、比较不同结构的偶极矢量、拼接分子） | 加 |
| 柔性扫描、IRC、反应路径（避免结构在步骤之间被旋转、编号对应混乱） | 常加 |
| ONIOM、外部程序接口 | 常加 |
| 对称性导致波函数被约束在错误的态上（第4章、第13章） | 加，或用 `Guess=Mix`、`Stable=Opt` |
| `NoSymm` 本身在识别点群时报错 | 改用 `Symmetry=None` |

> [!NOTE]
> 优化会保持初始结构的点群。如果初始结构对称性偏高，可能停在鞍点上，必须用频率计算确认。想打破对称性，可以手动把结构轻微扰动（例如移动某个原子 0.02 Å），或用 `NoSymm`。

---

## 13. 本章示例一览

| 文件 | 展示的知识点 |
|---|---|
| [`ethanol_cartesian.gjf`](../examples/02-input-anatomy/ethanol_cartesian.gjf) | 完整的 Link 0 + route + title + 笛卡尔坐标；B3LYP-D3(BJ) Opt Freq |
| [`h2o2_zmatrix.gjf`](../examples/02-input-anatomy/h2o2_zmatrix.gjf) | Z-matrix + `Variables:` 段，共用变量 |
| [`link1_opt_then_sp.gjf`](../examples/02-input-anatomy/link1_opt_then_sp.gjf) | `--Link1--`、`Geom=Check`、`Guess=Read`，high//low 两步 |
| [`fragment_labels_example.gjf`](../examples/02-input-anatomy/fragment_labels_example.gjf) | `(Fragment=n)` 标签和多片段电荷多重度行，`Counterpoise=2` |
| [`geom_connectivity_example.gjf`](../examples/02-input-anatomy/geom_connectivity_example.gjf) | GaussView 风格的 `Geom=Connectivity` 连接表 |
| [`check_charge_mult.py`](../examples/02-input-anatomy/check_charge_mult.py) | 电子数奇偶校验小工具 |

运行方法与预期结果见 [示例 README](../examples/02-input-anatomy/README.md)。

---

## 参数选择与收敛建议

| 项目 | 建议 | 理由 |
|---|---|---|
| route 前缀 | `#p` | 输出更全，诊断方便 |
| `%chk` | 每个作业都写，用有意义的名字 | 续算、读结构、作图都离不开它 |
| 读取旧结果 | `%oldchk` 读 + `%chk` 写 + `Geom=AllCheck Guess=Read` | 不覆盖旧 chk，也不必手抄坐标 |
| `%mem` | 申请内存的 75–85% | 给程序开销留余量，避免被调度系统 kill |
| `%nprocshared` | 等于申请的核数 | 避免超订 |
| 坐标格式 | 笛卡尔坐标；只在刚性扫描或需要变量约束时用 Z-matrix | 直观、便于脚本处理 |
| 连接表 | 普通 QM 计算可以不要；要删就连同 `Geom=Connectivity` 一起删 | 避免读错 |
| 关键词 | 写全称，不缩写，避免 IOp | 可读、可移植、可复现 |
| 多步作业 | 每步写全 Link 0；`--Link1--` 前留空行 | 避免读不到 chk |
| 对称性 | 默认开启；扫描、IRC、需要固定坐标系时 `NoSymm` | 兼顾效率和可控性 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `End of file in ZSymb.`（l101） | title 后或文件末尾缺空行；缺电荷多重度行；CRLF 换行 | 补空行；`dos2unix` |
| `QPErr --- A syntax error was detected in the input line.`（l1） | route 中关键词拼错、括号不配对、选项不属于该关键词 | 看输出中标记的位置；`testrt` |
| `The combination of multiplicity 1 and 17 electrons is impossible.`（l301） | 电荷或多重度与电子数奇偶性不符 | 按 §6 重新计算；用 `check_charge_mult.py` 检查 |
| `Small interatomic distances encountered`（l202） | 原子重叠、重复原子、Bohr/Å 单位混淆 | GaussView 中检查结构；核对 `Units` |
| 连接表相关的读取错误，或 Gen/ModRedundant 段被当成别的段读取 | 删了连接表但保留了 `Geom=Connectivity`（或反之），原子数改了但连接表没改 | 两者同时删除或同时保留并更新 |
| `Geom=Check` 失败，提示 chk 中找不到几何或文件不存在 | `%chk` 名与上一步不同；上一步失败没生成 chk；跨平台的 chk | 统一 `%chk` 名或用 `%oldchk`；跨机器先 `formchk` 再 `unfchk` |
| 第 2 步没有执行 | `--Link1--` 前缺空行，或前一步出错 | 补空行；先解决前一步的错误 |
| 结构在优化后"转了方向"，与输入坐标对不上 | 标准取向 | 用 `NoSymm`，或取输出中的 `Input orientation` |
| Counterpoise 或 `Guess=Fragment` 报错 | 片段电荷多重度行写错、原子漏标 Fragment、片段编号不连续 | 按 §6.3、§7.6 检查 |
| 冻结标记不起作用 | 单点计算中冻结标记会被忽略；或写在了 Z-matrix 行上 | 只在 Opt 中有效；Z-matrix 中用 `Constants:` 段 |

---

## 练习

1. 把 [`ethanol_cartesian.gjf`](../examples/02-input-anatomy/ethanol_cartesian.gjf) 改写成 Z-matrix 形式（可以用 GaussView 另存时不勾选 Write Cartesians，或用 `newzmat` 工具），比较两种输入优化得到的能量是否一致。
2. 在 `ethanol_cartesian.gjf` 中把两个碳原子的冻结标记设为 `-1`，重新优化，观察 C–C 键长是否保持不变。
3. 运行 [`link1_opt_then_sp.gjf`](../examples/02-input-anatomy/link1_opt_then_sp.gjf)，数一数 Normal termination 的行数，并解释为什么是 3。然后把第 2 步改写成 `%oldchk` + `Geom=AllCheck` 的形式。
4. 故意把 `fragment_labels_example.gjf` 的电荷多重度行改成 `0 1`（只写总体），运行并观察报错信息。
5. 计算下列体系的电子数和合理的多重度，并用 `check_charge_mult.py` 验证：NO、NO⁺、CN⁻、[Co(NH3)6]³⁺（低自旋）、Mn²⁺（高自旋）。

---

## 参考资料

**官方文档**

- [About Gaussian 16 Input](https://gaussian.com/input/)：语法规则、各段顺序（Section Ordering）、多步作业
- [Link 0 Commands](https://gaussian.com/link0/)：全部 `%` 命令与等价关系表
- [g16 Command Line Options](https://gaussian.com/options/)
- [Molecule Specifications](https://gaussian.com/molspec/)：电荷多重度、原子参数、片段、MM 类型、鬼原子
- [Constructing Z-Matrices](https://gaussian.com/zmat/)：Z-matrix、混合格式、哑原子
- [Geom keyword](https://gaussian.com/geom/)：Check、AllCheck、Step、Connectivity、ModConnectivity、ModRedundant 语法
- [Guess keyword](https://gaussian.com/guess/)、[Symmetry keyword](https://gaussian.com/symmetry/)、[Units keyword](https://gaussian.com/units/)
- [testrt](https://gaussian.com/testrt/)：route 语法检查
- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)：C.01 的 Geom=GenConnectivity、Geom=UseStandardOrientation，chk 版本号

**相关 GitHub 仓库**（仅外链，更多见 [相关仓库汇总](../../Resources/github-repos.md)）

- [mizu-bai/VSCode-Gaussian-Input-File-Support](https://github.com/mizu-bai/VSCode-Gaussian-Input-File-Support)：gjf 语法高亮，减少 route 拼写错误
- [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash)：含调用 testrt 的 route 检查脚本和从 chk 派生新输入的工具
- [Chempower-ZZH/auto_gjf](https://github.com/Chempower-ZZH/auto_gjf)：xyz/log 批量转 gjf 的轻量脚本
- [chunxiangzheng/gaussian_log_file_converter](https://github.com/chunxiangzheng/gaussian_log_file_converter)：从 log 提取结构转 xyz/gjf
- [maserasgroup-repo/pyssian](https://github.com/maserasgroup-repo/pyssian)：面向对象的 Gaussian 输入/输出解析库
- [theochem/iodata](https://github.com/theochem/iodata)：读写 fchk/log/wfn 等格式，并能生成 Gaussian 输入
- [thisisntnathan/dftCourseCodeRepo](https://github.com/thisisntnathan/dftCourseCodeRepo)：入门课程配套的 Gaussian 输入示例（CC0）

**书籍**

- J. B. Foresman, Æ. Frisch, *Exploring Chemistry with Electronic Structure Methods*, 3rd ed., Gaussian, Inc., 2015：第 1–2 章讲输入文件与 Z-matrix。

---

[← 上一章](01-setup-first-job.md) | [返回目录](../README.md) | [下一章 →](03-method-basis.md)
