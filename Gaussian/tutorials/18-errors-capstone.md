# 第18章 常见报错排查大全与完整研究流程案例（Gaussian Error Messages & a Complete Research Workflow）

本章分两部分。第一部分是 Gaussian 16 报错排查手册：先讲如何从 `Error termination via Lnk1e in .../lXXX.exe` 定位出错的 Link，再按“输入语法、几何、收敛、资源与环境”四类，逐条解释 L1 QPErr、L101 End of file in ZSymb、L103 FormBX、L202 distance matrix、L301 多重度/基组/色散/GenECP、L502 Convergence failure 与 CalDSu、L9999 Number of steps exceeded、L123 corrector、galloc、Erroneous write、Link 1 not found 等常见错误的原因和修正方法。第二部分用丁二烯 + 乙烯 Diels–Alder 反应演示一个完整的研究流程：CREST 构象搜索 → B3LYP-D3(BJ)/def2-SVP(SMD) 优化与频率 → 过渡态与 IRC → ωB97X-D/def2-TZVP(SMD) 单点 → GoodVibes qRRHO 与 1 M 标准态校正 → Boltzmann 加权 → 自由能剖面图，最后给出结果自查清单。关键词：Gaussian error termination、L502 convergence failure、L9999、L301、L101、Gaussian troubleshooting、reaction mechanism workflow、free energy profile、GoodVibes、CREST。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 看到任何 Gaussian 报错时，先确定出错的 Link，再找到真正的错误信息，而不是只盯着最后一行；
2. 知道常见错误属于哪一类（输入、几何、收敛、资源/环境），按对应的思路处理；
3. 对 SCF 不收敛、优化步数用完、IRC corrector 失败等“收敛类”问题，先诊断原因再选择方法，不盲目加大 `MaxCycle`；
4. 独立完成一个反应机理研究的完整计算流程，并得到可以写进论文的自由能剖面；
5. 用自查清单检查计算结果，避免常见的“算完了但结果不可信”的情况。

| 你遇到的情况 | 本章对应小节 |
|---|---|
| 作业秒退，日志很短 | 3 输入语法类 |
| 结构一开始就报错，或优化中途报坐标问题 | 4 几何与坐标类 |
| SCF、优化、IRC 跑了很久后失败 | 5 收敛类 |
| 作业被杀、没有报错信息、内存/磁盘相关报错 | 6 资源与环境类 |
| 只想快速查表 | 7 报错速查表 |
| 想知道一个完整课题该怎么算 | 8 完整研究流程 |
| 投稿前检查结果 | 9 自查清单 |

配套示例：

| 示例目录 | 内容 |
|---|---|
| [../examples/18-error-cases/README.md](../examples/18-error-cases/README.md) | 每个文件一个典型错误（L301 多重度、L101 缺空行、L1 QPErr、L301 R6DS8、GenECP 缺 ECP 段）及修正版，难收敛 SCF 的处理流程 |
| [../examples/18-capstone-diels-alder/README.md](../examples/18-capstone-diels-alder/README.md) | Diels–Alder 反应完整流程：CREST、优化频率、TS、IRC、高精度单点、GoodVibes、剖面图 |

本章大量引用前面各章的细节：SCF 收敛见 [第4章](04-sp-scf.md)，优化见 [第5章](05-opt.md)，频率与热化学见 [第6章](06-freq-thermo.md)，TS 与 IRC 见 [第7章](07-ts-irc.md)，扫描见 [第8章](08-scan.md)，溶剂见 [第9章](09-solvation.md)，开壳层与稳定性见 [第13章](13-openshell-stability.md)，作业资源见 [第17章](17-hpc-jobs.md)。

---

## 2. 定位错误的通用方法

### 2.1 先看最后几行

正常结束的作业，最后一行是：

```text
 Normal termination of Gaussian 16 at Tue Nov 10 10:00:00 2026.
```

多步作业（`Opt Freq`、`--Link1--`）每一步结束都会打印一次，所以 `grep -c "Normal termination" job.log` 应该等于作业步数。

出错的作业，最后会有类似这样的行（示意输出）：

```text
 Error termination via Lnk1e in /opt/g16/l502.exe at Tue Nov 10 10:00:00 2026.
 Job cpu time:       0 days  0 hours 12 minutes 34.5 seconds.
```

`l502.exe` 中的数字就是出错的 Link（程序模块）编号。**真正的错误原因在这一行上方几行到几十行**，不在这一行本身。最快的查看方法：

```bash
tail -n 30 job.log
grep -B 10 "Error termination" job.log | head -n 20
```

### 2.2 Link 编号速查

| Link | 负责的工作 | 在这里报错通常意味着 |
|---|---|---|
| L1 | 读取并解析 route，安排后续 Link | route 语法错误、关键词拼写错误或互相冲突 |
| L101 | 读取标题、电荷与多重度、分子坐标 | 空行缺失、电荷多重度行缺失、坐标格式错误 |
| L103 | Berny 优化器（构建内坐标、判断收敛） | 内坐标在线性角附近失效、TS 优化 Hessian 曲率不对 |
| L123 | IRC 积分 | IRC 的 corrector 步失败 |
| L202 | 坐标重定向、对称性、原子间距离矩阵 | 原子重叠、单位错误 |
| L301 | 建立基组、ECP、色散参数 | 多重度与电子数矛盾、基组不支持该元素、Gen/ECP 格式错误、色散组合不支持 |
| L401 | 初始猜测 | `Guess=Read` 读不到 chk、chk 与当前体系不匹配 |
| L502 | SCF 迭代 | SCF 不收敛、数值积分不准 |
| L508 | 二次收敛 SCF（`SCF=QC`） | QC 也不收敛 |
| L801 | 积分变换（冻芯设置等） | 后 HF、双杂化、TD 计算中的冻芯问题 |
| L914 | CIS/TD-DFT 激发态 | 激发态迭代问题 |
| L1002 | CPHF（频率、NMR、极化率） | CPHF 不收敛或数值问题 |
| L9999 | 结束作业、写 archive | 优化或扫描没有在给定步数内收敛 |

### 2.3 四类错误

| 类别 | 特征 | 典型 Link | 处理思路 |
|---|---|---|---|
| 输入语法 | 作业几秒内结束，日志很短 | L1、L101、L301 | 对照输入格式逐行检查；用 `testrt` 检查 route |
| 几何与坐标 | 一开始或优化中途报坐标/距离问题 | L202、L103 | 在 GaussView 中看结构；换坐标系 |
| 数值收敛 | 运行了一段时间后失败 | L502、L9999、L123、L1002 | 先诊断（看能量变化趋势），再对症调整 |
| 资源与环境 | 内存、磁盘报错，或作业被杀、没有任何报错 | 任意 Link，或没有 Error termination | 核对 `%mem`、scratch、作业时间与调度系统日志 |

### 2.4 工具与外部资源

- 仓库脚本 [`g16_error_doctor.py`](../scripts/g16_error_doctor.py)：自动识别出错的 Link 和常见报错，输出中文解释与建议，并可起草续算输入供人工检查（`python ../scripts/g16_error_doctor.py -h` 查看用法）；
- [`opt_monitor.py`](../scripts/opt_monitor.py)：画出每一步的能量和 4 个收敛判据，用于诊断 L9999；
- Gaussian 自带的 `testrt`：只检查 route 语法而不运行计算。把整条 route 作为一个参数（用引号括起来）传给它，例如 `testrt "#p B3LYP/6-31G(d) Opt Freq"`；有语法错误时会报 QPErr；
- [Solutions for Gaussian Error Messages](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/)：目前最全的报错整理（按 Link 编号），中英文；
- [计算化学公社 Gaussian 报错原帖](http://bbs.keinsci.com/thread-4829-1-1.html)：中文讨论与案例。

> [!TIP]
> 求助时请贴出：完整的输入文件（至少 route、电荷多重度和结构）、`Error termination` 之上 20–30 行、Gaussian 版本（日志开头或 archive 中的 `ES64L-G16RevC.01` 之类字样）。只贴最后一行，别人无法判断原因。

---

## 3. 输入语法类错误

### 3.1 L1：`QPErr --- A syntax error was detected in the input line.`

```text
 #p B3LYP/6-31G(d) Opt=(MaxCycles=100 Freq
 ...
 QPErr --- A syntax error was detected in the input line.
```

（示意输出）Gaussian 会在出错位置下方打印一个标记，指出解析失败的位置。

| 常见原因 | 例子 | 修正 |
|---|---|---|
| 括号不配对 | `Opt=(TS,CalcFC` | `Opt=(TS,CalcFC)` |
| 关键词或选项拼写错误 | `Freqq`、`EmpricalDispersion` | 对照 [关键词列表](https://gaussian.com/keywords/) |
| 泛函或基组名称写法不对 | 文献写法 def2-SVP、ωB97X-D、M06-2X 对应的 Gaussian 关键词是 `def2SVP`、`wB97XD`、`M062X` | 按 Gaussian 关键词列表的写法（不带连字符）书写，避免名称无法识别（[第3章](03-method-basis.md)） |
| 选项之间缺逗号或多了空格 | `SCF=(XQC MaxCycle=200)` | `SCF=(XQC,MaxCycle=200)` |
| route 前没有 `#` | 直接写 `B3LYP/6-31G(d)` | route 以 `#`、`#p` 或 `#t` 开头 |
| 中文全角字符 | 全角括号 `（`、全角等号 | 只用半角 ASCII |

同属 L1 的还有 `Illegal IType or MSType generated by parse.`，一般是关键词组合不合法（例如对不支持的方法请求某种计算），处理思路相同：简化 route，逐个加回关键词定位问题。示例：[`err_l1_qperr_route.gjf`](../examples/18-error-cases/err_l1_qperr_route.gjf)。

### 3.2 L101：`End of file in ZSymb.`

Gaussian 的每个输入段都以空行结束。程序读分子坐标时一直读到空行为止；如果文件在最后一个原子之后直接结束，就会读到文件结尾（EOF）而报这个错。

| 原因 | 修正 |
|---|---|
| 坐标后没有空行（最常见） | 文件末尾至少留一个空行 |
| 标题行与电荷多重度行之间缺空行，或缺少电荷多重度行 | 按“route / 空行 / 标题 / 空行 / 电荷多重度 / 坐标 / 空行”的顺序检查 |
| `--Link1--` 之后的步骤缺少标题或电荷多重度（用 `Geom=Check` 时仍需要） | 补上；或改用 `Geom=AllCheck`（此时不能写标题和电荷多重度） |
| Windows 编辑器保存的文件末尾问题 | 用 `cat -A job.gjf \| tail` 查看末尾是否有 `$` 空行 |

示例：[`err_l101_missing_blank_line.gjf.txt`](../examples/18-error-cases/err_l101_missing_blank_line.gjf.txt) 与修正版 [`err_l101_fixed.gjf`](../examples/18-error-cases/err_l101_fixed.gjf)。

### 3.3 L301：`The combination of multiplicity 1 and 9 electrons is impossible.`

电子数 N = 各原子序数之和 − 总电荷，多重度 M = 2S + 1。N 为奇数时 M 只能是 2、4、6……；N 为偶数时 M 只能是 1、3、5……。

| 体系 | 电子数 | 正确的电荷 多重度 |
|---|---|---|
| CH₃·（甲基自由基） | 6 + 3 = 9 | `0 2` |
| CH₃⁺ | 9 − 1 = 8 | `1 1` |
| CH₃⁻ | 9 + 1 = 10 | `-1 1` |
| O₂ 基态 | 16 | `0 3` |
| NO | 15 | `0 2` |
| [FeCl₄]⁻ 高自旋 Fe(III) | 26 + 68 + 1 = 95 | `-1 6` |

- 用了赝势时仍按全部电子计数（Gaussian 内部会扣除 ECP 替代的内层电子）；
- 开壳层体系用 `UB3LYP` 等非限制性方法，并检查 `<S**2>`（[第13章](13-openshell-stability.md)）；
- 仓库脚本 [`gjf_builder.py`](../scripts/gjf_builder.py) 生成输入时会检查电荷、多重度与电子数是否一致。

示例：[`err_l301_wrong_multiplicity.gjf`](../examples/18-error-cases/err_l301_wrong_multiplicity.gjf)。

### 3.4 L301：`Atomic number out of range for 6-31G basis set.`

Pople 基组（6-31G、6-311G 系列）只覆盖到较轻的元素，不包含第二、三过渡系金属和 I 等重元素。

| 处理 | 写法 |
|---|---|
| 全体系换成 def2 系列（Rb 之后的元素自动使用 def2-ECP） | `def2SVP`、`def2TZVP` |
| 重元素用 LANL2DZ 或 SDD，其余元素用 Pople 基组 | `GenECP`，见 3.6 |
| 全体系用 LANL2DZ 或 SDD | 直接写 `LANL2DZ`（轻元素部分质量较差，只适合粗略计算） |

### 3.5 L301：`R6DS8: Unable to choose the S8 parameter`

`EmpiricalDispersion=GD3BJ` 只对 Gaussian 内置了 BJ 阻尼参数的泛函有效（如 B3LYP、PBE0 等）。对没有参数的泛函会在 L301 报这个错。

| 泛函 | 正确的色散写法 |
|---|---|
| B3LYP、PBE1PBE（PBE0）等 | `EmpiricalDispersion=GD3BJ` |
| M06-2X 等 Minnesota 泛函 | `EmpiricalDispersion=GD3`（零阻尼 D3） |
| ωB97X-D（`wB97XD`）、B97-D3（`B97D3`）、B2PLYP-D3（`B2PLYPD3`）、PW6B95-D3（`PW6B95D3`）、APFD | 已经包含色散，**不要**再加 `EmpiricalDispersion` |

完整的参数表见 Gaussian 的 [DFT 方法文档](https://gaussian.com/dft/)。没有内置参数的泛函需要色散校正时，可以用 [pyDFTD3](https://github.com/patonlab/pyDFTD3) 等程序单独计算色散能。示例：[`err_l301_m062x_gd3bj.gjf`](../examples/18-error-cases/err_l301_m062x_gd3bj.gjf)。

### 3.6 L301：Gen/GenECP 格式错误

常见报错：

- `End of file reading basis center.`：基组段没有正常结束（缺 `****`、缺结尾空行）；
- `EOF while reading ECP pointer card.`：程序要读 ECP 段时遇到了文件结尾，常见于 ECP 段缺失或缺少结尾空行。

正确格式（[PdCl₄]²⁻，Cl 用 6-31+G(d)，Pd 用 LANL2DZ 基组和 ECP）：

```text
#p B3LYP/GenECP SP

[PdCl4]2- with GenECP

-2 1
Pd    0.000000    0.000000    0.000000
Cl    2.330000    0.000000    0.000000
Cl   -2.330000    0.000000    0.000000
Cl    0.000000    2.330000    0.000000
Cl    0.000000   -2.330000    0.000000

Cl 0
6-31+G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ

```

检查清单：

1. 只要有任何原子用 ECP，就写 `GenECP`（或 `Gen Pseudo=Read`），不能只写 `Gen`；
2. 基组段中每个元素块以 `****` 结束；分子中出现的每种元素都要列出，且只列一次；
3. 基组段与 ECP 段之间空一行；ECP 段只列需要 ECP 的元素；
4. 整个文件以空行结束（最好多留一行）；
5. 从 [Basis Set Exchange](https://www.basissetexchange.org/) 下载 Gaussian 格式的基组时，下载的文本本身可能已经包含 `****`，粘贴时注意不要重复或遗漏；
6. 元素块的首行可以写元素符号或原子编号，后面的 `0` 是结束标记。

示例：[`err_genecp_bad_format.gjf`](../examples/18-error-cases/err_genecp_bad_format.gjf)（缺 ECP 段）与 [`err_genecp_fixed.gjf`](../examples/18-error-cases/err_genecp_fixed.gjf)。Gen/GenECP 的完整说明见 [第3章](03-method-basis.md)。

### 3.7 GaussView 生成的输入与连接表

GaussView 默认在 route 中写 `Geom=Connectivity`，并在坐标后附上连接表（每行“原子编号 + 相连原子及键级”）。常见问题：

- 删除了连接表但保留了 `Geom=Connectivity`，或在文本编辑器中增删原子后连接表没有同步更新，导致读入错误；
- 手工修改原子顺序后，`ModRedundant` 约束、扫描坐标和 `Phase=` 中的原子编号全部错位。

处理：手工编辑结构后，要么同时删除 `Geom=Connectivity` 和连接表，要么回到 GaussView 中修改并重新保存。

---

## 4. 几何与坐标类错误

### 4.1 L202：`Small interatomic distances encountered` / `Problem with the distance matrix.`

原子之间距离过近（通常远小于 0.5 Å）或完全重叠。常见原因：

- 复制粘贴坐标时重复了某几行；
- 单位错误：从其他程序得到的 Bohr 坐标当成 Å 使用（或反之），需要 `Units=Bohr`；
- Z-matrix 中的变量写错；
- 从晶体结构、MD 快照截取团簇时，周期性镜像原子重复。

处理：在 GaussView 中打开输入，看是否有原子叠在一起；用 Tools → Bond 或测量功能检查最短距离。`Small interatomic distances encountered:` 之后列出的原子编号就是需要检查的原子对。

### 4.2 L103：`FormBX had a problem` / `Linear angle in Bend` / `Linear angle in Tors`

优化默认使用冗余内坐标。当某个键角接近 180°（炔烃、腈、M–C≡O、某些 TS 中的 X···H···Y）时，键角和二面角的定义发生奇异，Gaussian 打印 `Bend failed for angle  1 - 11 - 3`、`Tors failed for dihedral ...` 等信息后报错。

| 方法 | 写法 | 说明 |
|---|---|---|
| 从最后的结构重启 | `Geom=Check Guess=Read Opt` | Gaussian 会按新结构重新生成内坐标，很多时候直接解决 |
| 改用笛卡尔坐标优化 | `Opt=Cartesian` | 收敛可能慢一些，但没有线性角问题 |
| 加哑原子 | 在线性片段旁放一个 `X` 哑原子，用 Z-matrix 定义 | 传统方法，适合手工构建的线性分子 |
| 扫描时避免经过 180° | 例如把二面角扫描范围调整为不跨越线性构型 | 见 [第8章](08-scan.md) |

### 4.3 结构合理但“不对”的情况

这些不一定报错，但结果没有意义：

- 原子没连成你以为的分子（例如两个碎片距离太远，优化成两个独立分子）；
- 优化后跳到另一个构象或发生了意外的反应（如质子转移）；
- 对称性限制使结构停在高对称的鞍点上（频率出现虚频）。

优化结束后**一定**要在 GaussView 中看最终结构，频率中有虚频时播放其振动动画。

---

## 5. 收敛类错误

### 5.1 L502：`Convergence failure -- run terminated.`

SCF 在最大迭代次数内没有收敛。用 `#p` 时日志会打印每一轮迭代的能量变化，先看趋势：

```bash
grep -E "^ (Cycle|E=)" job.log | tail -n 40
```

| 能量变化趋势 | 可能原因 | 推荐处理（按顺序） |
|---|---|---|
| 剧烈振荡 | 初猜差、HOMO–LUMO 能隙小 | `SCF=XQC`；`Guess=Read` 从小基组或相关闭壳层离子读入；`SCF=VShift=500`；`SCF=Fermi` |
| 一开始就很怪（能量异常高或乱跳） | 电荷/多重度错误、结构不合理、基组/ECP 写错 | 先检查输入，修正后重算 |
| 单调缓慢下降 | 体系大、收敛慢（较少见） | 只有这种情况才加大 `SCF=(MaxCycle=512)` |
| 很快收敛到某个值但随后的计算有问题 | 收敛到激发态或不稳定的解 | `Stable=Opt`（[第13章](13-openshell-stability.md)） |

要点：

- **不要一上来就加大 MaxCycle**。振荡的 SCF 迭代再多也不会收敛，只是浪费机时；
- `SCF=XQC` 在常规 DIIS 没有收敛时自动切换到二次收敛算法，是最常用的第一步；如果 QC 本身也失败，会在 L508 报 `Convergence failure`，此时去掉 QC 类选项，改用其他方法；
- 好的初猜往往最有效：先用小基组（def2-SVP）或闭壳层的阳/阴离子收敛，再 `Guess=Read`；
- 开壳层过渡金属体系还可以试 `Guess=Mix`、`Guess=Fragment`、`SCF=NoIncFock`；
- 示例 [`err_l502_hard_scf_fixed.gjf`](../examples/18-error-cases/err_l502_hard_scf_fixed.gjf) 演示了“小基组 + `SCF=XQC` + `Stable=Opt` → 大基组 `Guess=Read`”的标准流程。更多方法见 [第4章](04-sp-scf.md)。

### 5.2 L502（或 L1002）：`Inaccurate quadrature in CalDSu.`

DFT 数值积分得到的电子数与实际电子数偏差过大。

1. 先检查低级错误：赝势基组忘了写 ECP、元素用错基组、结构中有原子重叠或距离异常；
2. 加密积分格点：G16 默认已是 `Int=UltraFine`，可以改为 `Int=SuperFine`（同一项目中所有计算要统一格点）；
3. `SCF=NoVarAcc` 让 SCF 全程使用完整的积分精度；
4. 极弥散的基函数与拥挤的结构同时出现时容易出问题，考虑去掉不必要的弥散函数。

### 5.3 L9999：`Number of steps exceeded`

```text
 Optimization stopped.
    -- Number of steps exceeded,  NStep= 100
 ...
 Error termination via Lnk1e in /opt/g16/l9999.exe at ...
```

（示意输出）优化在 `MaxCycles` 步内没有收敛。**先画图诊断**：GaussView 打开 log → Results → Optimization，或用 [`opt_monitor.py`](../scripts/opt_monitor.py) 画出能量和 4 个判据随步数的变化。

| 诊断结果 | 处理 |
|---|---|
| 能量单调下降、力在减小，只是还没收敛 | `Opt=(Restart,MaxCycles=200)` 续算，或 `Geom=Check Guess=Read Opt=(MaxCycles=200)` 从最后结构重新开始 |
| 能量来回振荡 | `Opt=CalcFC`（计算准确的初始 Hessian）；减小步长 `Opt=(MaxStep=10)`；`Int=SuperFine` 降低数值噪声 |
| 力已经很小但位移不收敛（平坦势能面，常见于弱相互作用复合物、甲基转动） | 检查 `Converged?` 表；Gaussian 在力远小于阈值时会接受收敛；必要时 `Opt=CalcFC`、`Int=SuperFine` |
| 结构越走越离谱（键断开、原子飞走） | 初始结构或电荷多重度有问题，或连接性错误；重新构建 |
| TS 优化在两个区域间跳来跳去 | 初猜质量差，回到扫描或 QST2/QST3（[第7章](07-ts-irc.md)） |

### 5.4 L123：`Maximum number of corrector steps exceded.`（IRC）

（原文中的 “exceded” 拼写如此）IRC 的 corrector 步没有收敛。

1. 确认 TS 本身是收敛的、恰好 1 个虚频，且 IRC 用的是同一水平的力常数（`RCFC`）；
2. 加大 `IRC=(MaxCycle=100)`；
3. 改用 LQA 积分器：`IRC=(RCFC,LQA)`；或用 `IRC=(RCFC,GS2,StepSize=5)`（GS2 是 G03 的默认算法，配合更小的步长）；
4. 让 IRC 定期重算 Hessian：`IRC=(CalcFC,ReCalc=5)`（需要方法有解析二阶导数）；
5. IRC 走到很平的区域时 corrector 失败，往往已经接近端点：此时可以直接取最后一个点做 `Opt` 确认连接关系。

### 5.5 TS 优化：`Wrong number of Negative eigenvalues: Desired= 1 Actual= 0`

TS 优化要求 Hessian 有且只有一个负本征值。初始 Hessian（默认是估计的）曲率不对时会报这个错。

- 用 `Opt=(TS,CalcFC,NoEigenTest)`：`CalcFC` 计算准确的初始 Hessian，`NoEigenTest` 让优化在负本征值个数不对时继续；
- 初猜更接近 TS：用柔性扫描的能量最高点（[第8章](08-scan.md)），或 QST2/QST3；
- 无论哪种方法得到的 TS，最后都必须用频率确认**恰好 1 个虚频**，且虚频振动对应反应坐标，并用 IRC 确认两端（[第7章](07-ts-irc.md)）。

### 5.6 “收敛了”但有小虚频

优化正常结束，频率中却出现一个几十 cm⁻¹ 的小虚频：

- 频率必须与优化使用**完全相同**的方法、基组、溶剂模型和积分格点，否则得到的不是驻点；
- 收紧收敛：`Opt=Tight` 或 `VeryTight`，配合 `Int=SuperFine`；
- 沿虚频方向位移结构后重新优化（GaussView 中按虚频振动方向手动调整，或用 [pyQRC](https://github.com/patonlab/pyQRC)）；
- 低频的甲基转动等模式残留很小的虚频（如 −10 cm⁻¹ 量级）通常是数值噪声，但应在论文中说明处理方式。

### 5.7 其他收敛与数值问题

| 报错 | Link | 原因与处理 |
|---|---|---|
| `Convergence failure.`（使用 `SCF=QC`/`XQC` 时） | L508 | 二次收敛 SCF 也失败；去掉 QC 类选项，改用 `Guess=Read`、`VShift` 等方法 |
| `Excessive mixing of frozen core and valence orbitals.` | L801 | 后 HF、双杂化或 TD 计算中的冻芯设置不合理；改用全电子相关（如 `MP2(Full)`），或调整冻芯设置 |
| CPHF 相关的不收敛 | L1002 | 频率、NMR、极化率中求解 CPHF 失败；检查 SCF 是否稳定、格点是否足够 |

---

## 6. 资源与环境类错误

### 6.1 内存

| 报错或现象 | 原因 | 处理 |
|---|---|---|
| `galloc:  could not allocate memory.` | `%mem` 超过节点当前可用内存 | 降低 `%mem`；`%mem` 取调度系统申请内存的 75–85% |
| 作业被调度系统杀掉，日志中断，没有 `Error termination` | 实际内存超出申请值 | `sacct` 查看状态；`%mem` 不要接近 `--mem` |
| 提示需要更多内存的 words，或自动减少并行线程数 | `%mem` 太小（G16 默认只有 800 MB） | 增加 `%mem`；检查是否被 Link 0 或 `Default.Route` 覆盖 |

详见 [第17章](17-hpc-jobs.md) 第 6 节。

### 6.2 磁盘：`Erroneous write. Write -1 instead of 4096.`

（示意输出，数字因作业而异）写 scratch 文件失败，通常是磁盘满或超出配额。

1. `df -h $GAUSS_SCRDIR` 检查剩余空间，`quota -s` 检查配额；
2. 清理以前异常结束遗留的 `Gau-*.rwf` 等文件；
3. 把 `GAUSS_SCRDIR` 换到更大的本地盘；用 `%RWF=/disk1/,50GB,/disk2/,-1` 拆分；
4. `GAUSS_SCRDIR` 指向的目录不存在时，作业会立即失败，作业脚本中要先 `mkdir -p`。

### 6.3 `segmentation violation`（Linux）与 `Severe Error Message #2070`（Windows）

这两种信息本身**不说明原因**，只表示程序异常退出。要往上看日志中最后一个具体的错误或警告。常见的根本原因包括：内存设置超过物理内存、磁盘写满、输入中的特殊错误、个别版本的程序缺陷。排查顺序：先看上方的具体信息；降低 `%mem` 和核数重试；换一台机器或换 Gaussian 修订版本验证是否可重复。

### 6.4 `Link 1 not found` / 找不到 `l1.exe`

Gaussian 无法启动它的第一个程序模块：

| 平台 | 常见原因 | 处理 |
|---|---|---|
| Linux | 没有 `source $g16root/g16/bsd/g16.profile`；`g16root` 或 `GAUSS_EXEDIR` 设置错误；作业脚本中没有加载环境 | 在作业脚本中显式加载环境（[第1章](01-setup-first-job.md)） |
| Windows | 输入文件或工作目录的路径中含有中文或空格；安装目录不对；杀毒软件隔离了可执行文件；`%mem` 设得过大 | 把文件放到纯英文、无空格的路径（如 `D:\g16work\`）；检查 G16W 首选项中的工作目录与 scratch 设置；为 Gaussian 目录添加杀毒软件白名单 |

### 6.5 chk 相关：`NtrErr Called from FileIO.`

这一行之前往往有一大段数字，真正的报错在更上方。常见原因：`Geom=Check`/`Guess=Read` 指定的 chk 不存在或名字写错（多步作业中各步 `%chk` 不一致）；chk 来自其他平台或其他 Gaussian 版本；chk 已损坏（上一个作业写 chk 时被杀）。处理：核对文件名，用 `%oldchk` 读、`%chk` 写；跨平台或跨版本时先 `formchk` 再 `unfchk`，C.01 之前的 chk 用 `c8616` 升级（[第16章](16-postprocess-multiwfn.md)）。

---

## 7. 报错速查表

| 报错信息（或关键片段） | Link | 一句话原因 | 首选处理 |
|---|---|---|---|
| `QPErr --- A syntax error was detected in the input line.` | L1 | route 语法错误 | 按标记位置检查括号、拼写；`testrt` |
| `Illegal IType or MSType generated by parse.` | L1 | 关键词组合不合法 | 简化 route，逐个加回 |
| `End of file in ZSymb.` | L101 | 坐标后缺空行或缺电荷多重度行 | 补空行，检查分段 |
| `FormBX had a problem.` / `Linear angle in Bend.` / `Linear angle in Tors.` | L103 | 内坐标遇到线性角 | 从最后结构重启；`Opt=Cartesian` |
| `Wrong number of Negative eigenvalues` | L103 | TS 优化的 Hessian 曲率不对 | `Opt=(TS,CalcFC,NoEigenTest)`；改进初猜 |
| `Maximum number of corrector steps exceded.` | L123 | IRC corrector 失败 | 加大 `MaxCycle`；`LQA`；`GS2` + 小步长 |
| `Small interatomic distances encountered` / `Problem with the distance matrix.` | L202 | 原子重叠、单位错误 | 在 GaussView 中检查结构 |
| `The combination of multiplicity M and N electrons is impossible.` | L301 | 多重度与电子数奇偶不符 | 修正电荷与多重度 |
| `Atomic number out of range for ... basis set.` | L301 | 基组不支持该元素 | def2 系列、LANL2DZ、SDD 或 GenECP |
| `R6DS8: Unable to choose the S8 parameter` | L301 | 该泛函没有所选色散参数 | Minnesota 泛函用 `GD3`；自带色散的泛函不加 |
| `End of file reading basis center.` | L301 | Gen 基组段格式错误或缺结尾空行 | 检查 `****` 和空行 |
| `EOF while reading ECP pointer card.` | L301 | ECP 段缺失或缺结尾空行 | 补 ECP 段和空行 |
| `Convergence failure -- run terminated.` | L502 | SCF 不收敛 | 检查输入 → `SCF=XQC` → `Guess=Read` → `VShift` |
| `Inaccurate quadrature in CalDSu.` | L502/L1002 | DFT 数值积分不准 | 检查基组/ECP/结构；`Int=SuperFine` |
| `Convergence failure.`（QC） | L508 | QC SCF 也失败 | 去掉 QC，换其他方法 |
| `Excessive mixing of frozen core and valence orbitals.` | L801 | 冻芯不合理 | 全电子相关或调整冻芯 |
| `Number of steps exceeded, NStep= ...` | L9999 | 优化步数用完 | 先诊断，再续算或改设置 |
| `galloc:  could not allocate memory.` | 任意 | `%mem` 超过可用内存 | 降低 `%mem` |
| `Erroneous write. Write -1 instead of ...` | 任意 | 磁盘满 | 清理或更换 scratch |
| `NtrErr Called from FileIO.` | 任意 | chk/rwf 读写问题 | 检查 chk 名称、版本、平台 |
| `segmentation violation` / `Severe Error Message #2070` | 任意 | 通用的异常退出 | 往上找具体原因 |
| `Link 1 not found` | L1 之前 | 环境变量、路径问题 | source g16.profile；纯英文路径 |
| 无 `Error termination`，日志突然中断 | — | 被调度系统杀掉（时间、内存） | `sacct`/`qstat -f` 查看原因；续算 |

---

## 8. 完整研究流程案例：丁二烯 + 乙烯 Diels–Alder 反应

### 8.1 研究问题与流程总览

研究问题：在甲苯溶液中，1,3-丁二烯与乙烯发生 [4+2] 环加成生成环己烯的活化自由能 ΔG‡ 和反应自由能 ΔG<sub>r</sub> 是多少？

这是一个只有 16 个原子的经典周环反应，计算量小，但包含了一个完整课题的所有环节：

| 步骤 | 做什么 | 工具与水平 | 输入文件（示例目录） |
|---|---|---|---|
| 0 | 确定方法、基组、溶剂和整体方案 | — | — |
| 1 | 构象搜索 | CREST，GFN2-xTB/ALPB(toluene) | `crest_conformers.sh` |
| 2 | 反应物、产物优化 + 频率 | B3LYP-D3(BJ)/def2-SVP，SMD(toluene) | `reactant_opt.gjf`、`reactant2_opt.gjf`、`product_opt.gjf` |
| 3 | TS 优化 + 频率 | 同上，`Opt=(TS,CalcFC,NoEigenTest)` | `ts_opt.gjf` |
| 4 | IRC 验证 | 同上，`IRC=(RCFC,...)` | `irc.gjf` |
| 5 | 高精度单点 | ωB97X-D/def2-TZVP，SMD(toluene) | `sp_highlevel.gjf`、`make_sp_inputs.sh` |
| 6 | 热化学：qRRHO + 1 M 标准态 | GoodVibes（或 Shermo、`thermo_qrrho.py`） | `goodvibes_run.sh` |
| 7 | 多构象时 Boltzmann 加权 | GoodVibes `--boltz` 等 | — |
| 8 | 自由能剖面图 | GoodVibes `--pes`，或 `plot_profile.py` | `da_profile.yaml`、`plot_profile.py` |

整个流程的 Gaussian 部分可以用 [`run_workflow_slurm.sh`](../examples/18-capstone-diels-alder/run_workflow_slurm.sh) 一次提交，每一步检查失败即停止。

### 8.2 第 0 步：选择计算方案

| 决策 | 本例选择 | 理由 |
|---|---|---|
| 优化与频率 | B3LYP-D3(BJ)/def2-SVP | 有机反应中结构和频率的常用组合，成本低；D3(BJ) 改善非共价作用（[第3章](03-method-basis.md)） |
| 溶剂 | SMD(toluene)，优化和频率都在溶剂中做 | 结构与溶剂一致；甲苯是 Diels–Alder 反应的常见溶剂（[第9章](09-solvation.md)） |
| 单点 | ωB97X-D/def2-TZVP，SMD(toluene) | 范围分离杂化 + 内含色散，是有机反应能垒单点的常用选择；也可以用双杂化 DSD-PBEP86/def2-TZVP 做对照。正式研究前应在小模型上做泛函测试（8.11 节） |
| 热校正 | 低水平频率 + qRRHO（Grimme，截断 100 cm⁻¹） | TS 和弱复合物中的低频模式会使谐振子熵偏大 |
| 标准态 | 1 mol/L | 溶液反应；双分子反应的能垒受影响 |
| 温度 | 298.15 K | 与实验条件一致时再改 |

论文中写作 ωB97X-D/def2-TZVP(SMD)//B3LYP-D3(BJ)/def2-SVP(SMD)，“//”左边是单点水平，右边是结构优化水平。

> [!IMPORTANT]
> 同一个反应剖面中的所有物种必须使用**完全相同**的方案：相同的泛函、基组、溶剂模型、积分格点、温度和标准态。只要有一个物种不同（例如某个中间体单独加了 `Int=SuperFine`，或忘记加 SMD），得到的相对能量就没有意义。

### 8.3 第 1 步：CREST 构象搜索

丁二烯有两个构象：s-trans（能量最低）和 gauche-s-cis（发生 Diels–Alder 反应的构象，能量更高）。如果直接用 s-cis 作为反应物参考，会低估能垒。对于更复杂的底物（取代基、柔性链），构象数可以有几十上百个，必须系统搜索。

```bash
NP=8 bash crest_conformers.sh          # crest butadiene_strans.xyz --gfn2 --alpb toluene -T 8 --imtdgc
```

- CREST 3.x 没有默认运行类型，需要用 `--imtdgc` 等参数显式指定（`--quick`、`--squick` 等会隐式选择 iMTD-GC）；
- `--alpb toluene` 在构象搜索阶段就考虑溶剂；默认能量窗口为 6 kcal/mol（`--ewin` 修改）；
- 结果 `crest_conformers.xyz` 按能量排序，脚本把每个构象转成 `conf_NN_opt.gjf`；
- GFN2-xTB 的能量排序只是初筛：**所有在窗口内的构象都要用 DFT 重新优化**，再按 DFT 自由能排序；
- TS 也有构象（取代反应物的 endo/exo、侧链朝向），复杂体系可以在固定反应中心的情况下对 TS 做约束构象搜索（CREST 支持约束，见其文档）。

### 8.4 第 2 步：反应物与产物的优化和频率

```bash
g16 < reactant_opt.gjf  > reactant_opt.log
g16 < reactant2_opt.gjf > reactant2_opt.log
g16 < product_opt.gjf   > product_opt.log
grep "imaginary frequencies" reactant_opt.log reactant2_opt.log product_opt.log   # expect no output
```

route（三个文件相同）：

```text
#p B3LYP EmpiricalDispersion=GD3BJ def2SVP SCRF=(SMD,Solvent=Toluene) Opt Freq
```

检查：两步都 `Normal termination`；0 个虚频；结构合理（丁二烯为平面 s-trans，环己烯为半椅式）。

### 8.5 第 3 步：过渡态

`ts_opt.gjf` 中的初猜是手工搭建的 C<sub>s</sub> 结构：s-cis 丁二烯平面在下，乙烯在上，两根正在形成的 C–C 键都设为 2.27 Å，末端 CH₂ 的 H 略微向外弯曲。

```text
#p B3LYP EmpiricalDispersion=GD3BJ def2SVP SCRF=(SMD,Solvent=Toluene) Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq
```

检查：

```bash
grep "imaginary frequencies" ts_opt.log      # exactly: 1 imaginary frequencies
grep -m 1 "Frequencies --" ts_opt.log        # the first value is negative
```

在 GaussView 中播放虚频：两根 C···C 键应同步伸缩，对应协同的 [4+2] 成键。如果得到 0 个虚频（优化到了极小点）或 2 个虚频，回到 [第7章](07-ts-irc.md) 的方法：用 `Geom=GIC` 定义两根键长之和做柔性扫描（[第8章](08-scan.md)），以能量最高点为初猜；或用 QST2。

### 8.6 第 4 步：IRC

```text
%oldchk=ts_opt.chk
%chk=irc.chk
#p B3LYP EmpiricalDispersion=GD3BJ def2SVP SCRF=(SMD,Solvent=Toluene) IRC=(RCFC,Phase=(1,6),MaxPoints=50,StepSize=10,MaxCycle=50) Geom=Check Guess=Read
```

- `RCFC` 读取 TS 频率计算保存在 chk 中的笛卡尔力常数；
- `Phase=(1,6)` 规定 C1–C6 距离增大为正向（走向反应物），使正反方向的含义明确；
- 默认双向积分；`MaxPoints` 为每个方向的点数。

检查：`Summary of reaction path following` 表中 TS 位于中间，两侧能量单调下降；一端两根键都断开（s-cis 丁二烯 + 乙烯的松散复合物），另一端成环（环己烯）。如需严格确认，把两端结构取出（GaussView 中另存，或用 [`log2xyz.py`](../scripts/log2xyz.py)）再做 `Opt Freq`，与 `product_opt` 和反应物比较。

> [!NOTE]
> IRC 反应物端得到的是 s-cis 丁二烯与乙烯的弱复合物，它不是能垒的参考态。能垒参考的是分离的、各自处于最低能构象的反应物（s-trans 丁二烯 + 乙烯），并在 1 M 标准态下计算。

### 8.7 第 5 步：高精度单点

```bash
bash make_sp_inputs.sh          # reactant_opt_TZ.gjf, reactant2_opt_TZ.gjf, ts_opt_TZ.gjf, product_opt_TZ.gjf
for r in reactant_opt reactant2_opt ts_opt product_opt; do g16 < ${r}_TZ.gjf > ${r}_TZ.log; done
```

模板 `sp_highlevel.gjf`：

```text
%oldchk=ts_opt.chk
%chk=ts_opt_TZ.chk
#p wB97XD/def2TZVP SCRF=(SMD,Solvent=Toluene) Geom=AllCheck Guess=Read

```

`Geom=AllCheck` 从 chk 中读取结构、标题、电荷和多重度，所以输入中不写这些内容。文件名采用 `X.log` 与 `X_TZ.log` 的配对方式，GoodVibes 用 `--spc TZ` 自动配对。

### 8.8 第 6 步：热化学——qRRHO 与 1 M 标准态

溶液中每个物种的标准自由能：

$$
G^\circ(1\,\mathrm{M}) = E_{\mathrm{SP}}^{\text{high}} + G_{\mathrm{corr}}^{\mathrm{qRRHO}}(\text{low}) + RT\ln(24.46)
$$

- $E_{\mathrm{SP}}^{\text{high}}$：ωB97X-D/def2-TZVP(SMD) 单点能；
- $G_{\mathrm{corr}}^{\mathrm{qRRHO}}$：B3LYP-D3(BJ)/def2-SVP(SMD) 频率给出的热校正，低频模式（< 100 cm⁻¹）的振动熵用 Grimme 的准刚性转子插值代替谐振子；
- $RT\ln(24.46)$：从 1 atm 理想气体标准态（Gaussian 的默认）换算到 1 mol/L，298.15 K 时为 +1.89 kcal/mol（24.46 L 是 298.15 K、1 atm 下 1 mol 理想气体的体积）。

```bash
bash goodvibes_run.sh
# inside: goodvibes reactant_opt.log reactant2_opt.log ts_opt.log product_opt.log \
#           --spc TZ --temp 298.15 --conc 1.0 -f 100 --qs grimme --pes da_profile.yaml
```

GoodVibes 输出的表格中，`qh-G(T)_SPC` 一列即上式的 G°（Hartree）。脚本在运行前会检查所有作业正常结束、极小点 0 个虚频、TS 恰好 1 个虚频。GoodVibes 会根据识别到的方法和基组自动选择频率缩放因子（`-v 1.0` 关闭缩放），论文中要写明。

反应能垒和反应自由能：

$$
\Delta G^\ddagger = G^\circ(\mathrm{TS}) - G^\circ(\text{diene}) - G^\circ(\text{ene}),\qquad
\Delta G_r = G^\circ(\mathrm{P}) - G^\circ(\text{diene}) - G^\circ(\text{ene})
$$

单位换算：1 E<sub>h</sub> = 627.5095 kcal/mol = 2625.50 kJ/mol = 27.2114 eV；298.15 K 时 RT = 0.5925 kcal/mol。

> [!WARNING]
> 双分子反应（A + B → TS）中，反应物一侧有两个物种、TS 只有一个，所以标准态换算不会抵消：ΔG‡(1 M) = ΔG‡(1 atm) − 1.89 kcal/mol。单分子步骤（如中间体 → TS）中两边各一个物种，换算会抵消。忘记或重复这一校正是常见错误。

其他工具：[Shermo](http://sobereva.com/soft/shermo/)（Sobereva 开发，可在 settings.ini 中设置低频处理方法与浓度，详见其手册）、[OpenThermo](https://github.com/lenhanpham/OpenThermo)，以及仓库的教学脚本 [`thermo_qrrho.py`](../scripts/thermo_qrrho.py)（`-h` 查看用法）。不同程序的 qRRHO 实现细节（截断频率、插值函数、转动惯量的取法）略有差别，同一项目中只用一种。

### 8.9 第 7 步：Boltzmann 加权

当某个物种有多个低能构象时（本例中丁二烯的 s-cis 构象能量较高，贡献很小；更复杂的底物则不然），有两种处理方式：

1. **只用最低能构象**：最常用，适用于最低构象明显占优的情况；
2. **构象系综自由能**：把各构象看作同一物种的不同状态，

$$
p_i=\frac{e^{-G_i/RT}}{\sum_j e^{-G_j/RT}},\qquad
G_{\text{ens}} = -RT\ln\sum_i e^{-G_i/RT}
$$

$G_{\text{ens}}$ 总是不高于最低构象的 $G$，差值反映构象熵的贡献。光谱（ECD、NMR）等性质用 $p_i$ 加权平均（[第10章](10-tddft.md)、[第11章](11-nmr-props.md)）。GoodVibes 的 `--boltz` 输出各构象的 Boltzmann 权重，reaction-profile 文档中也可以按构象汇总（见其文档中的 gconf 选项）。

### 8.10 第 8 步：自由能剖面图

GoodVibes ≥ 4.4 按 `da_profile.yaml` 中的定义输出各点相对于反应物的 ΔG 与 Δqh-G（kcal/mol）；安装了绘图依赖时，加 `--pes-plot da_profile.svg` 可直接出图。也可以把数值填入仓库示例中的简单脚本：

```bash
python plot_profile.py R=0 TS=<dG_TS> P=<dG_P> -o da_profile.png
```

作图规范：横轴为反应坐标（不需要数值），纵轴标明 ΔG（kcal/mol）、温度、标准态和计算水平；TS 标上 ‡；数值保留一位小数。

### 8.11 与实验比较

该反应在气相中的实验活化能约为 27.5 kcal/mol（Rowley & Steiner, 1951，见参考资料）。比较时注意：

- 实验 E<sub>a</sub> 与 ΔG‡ 是不同的量。由过渡态理论，气相双分子反应有 E<sub>a</sub> = ΔH‡ + 2RT；ΔG‡ 还包含 −TΔS‡，双分子反应的活化熵是很大的负值，所以 ΔG‡ 明显高于 ΔH‡；
- 本例在 SMD(甲苯) 中计算，而实验数据来自气相，溶剂对这一非极性反应的影响不大但不为零；
- 由 ΔG‡ 可以用 Eyring 方程估计速率常数：

$$
k=\frac{k_\mathrm{B}T}{h}\,e^{-\Delta G^\ddagger/RT}
$$

298.15 K 时 $k_\mathrm{B}T/h \approx 6.21\times10^{12}\ \mathrm{s^{-1}}$，ΔG‡ 每变化 1.36 kcal/mol，k 变化 10 倍。

DFT 能垒的典型误差在 1–3 kcal/mol 量级，不同泛函之间差异可能更大。周环反应的系统基准可参考 Houk 等人的工作（参考资料），做正式研究时应在小模型上比较几种泛函。

### 8.12 支撑材料应包含什么

- 所有驻点的笛卡尔坐标（Å）、电子能、单点能、ZPE、焓、qh-G，以及虚频（TS）；
- 完整的计算方案：程序与版本（如 Gaussian 16 Rev. C.01）、泛函、基组、色散、溶剂模型、积分格点、温度、标准态、频率缩放因子、qRRHO 方法与截断频率、使用的后处理程序及版本；
- 构象搜索的方法与能量窗口；
- 引用：Gaussian、所用泛函与基组、D3(BJ)、SMD、CREST、GoodVibes/Shermo 等。

GoodVibes 的 `--si` 选项可以输出包含坐标附录的 SI 表格；仓库脚本 [`g16_log_summary.py`](../scripts/g16_log_summary.py) 可以把一批 log 汇总为 CSV/JSON。

---

## 9. 结果自查清单

| 检查项 | 怎么查 | 合格标准 |
|---|---|---|
| 正常结束 | `grep -c "Normal termination" x.log` | 等于作业步数（`Opt Freq` 为 2） |
| 优化真正收敛 | `grep -A 5 "Converged?" x.log \| tail -n 6` | 4 个判据为 `YES`，或 Gaussian 明确接受收敛 |
| 虚频个数 | `grep "imaginary frequencies" x.log` | 极小点：无输出；TS：恰好 1 个 |
| 虚频模式 | GaussView 播放虚频 | 对应目标反应坐标 |
| IRC 连通性 | IRC 两端结构 + 端点 `Opt Freq` | 一端是反应物（或其复合物），另一端是产物 |
| 自旋污染（开壳层） | `grep "S\*\*2" x.log` | `<S**2>` 与 S(S+1) 相差不超过约 10% |
| 波函数稳定性 | `Stable` 计算 | `The wavefunction is stable under the perturbations considered.` |
| 方案一致性 | 比较所有 log 的 route | 方法、基组、溶剂、格点、温度完全一致 |
| Opt 与 Freq 同一水平 | 比较 route | 频率与优化的水平和溶剂完全相同 |
| 构象 | CREST + DFT 重排序 | 使用的是最低能构象（或已做 Boltzmann 加权） |
| 标准态 | 检查热化学处理 | 溶液反应用 1 M；双分子步骤的校正不能遗漏 |
| 低频处理 | 检查 GoodVibes/Shermo 设置 | 说明 qRRHO 方法和截断频率 |
| 单位换算 | 重算一遍 | 1 E<sub>h</sub> = 627.5095 kcal/mol = 2625.50 kJ/mol |
| 电荷与多重度 | 检查输入 | 与化学结构一致，电子数奇偶正确 |
| 版本与可重复性 | 日志开头或 archive 中的版本号 | 同一项目使用同一 Gaussian 版本；比较旧文献时注意 G09 与 G16 默认值不同（格点、SCRF 等） |

---

## 参数选择与收敛建议

| 场景 | 推荐设置 | 说明 |
|---|---|---|
| SCF 不收敛的第一步 | `SCF=XQC` | 常规 DIIS 失败后自动切换二次收敛 |
| 难收敛体系的初猜 | 小基组或闭壳层离子收敛后 `Guess=Read` | 往往比任何 SCF 选项都有效 |
| 开壳层 | `Stable=Opt`，检查 `<S**2>` | 防止收敛到不稳定或错误的态 |
| 优化振荡 | `Opt=CalcFC`、`MaxStep=10`、`Int=SuperFine` | 先诊断再选择 |
| 优化缓慢下降 | `Opt=(Restart,MaxCycles=200)` | chk 必须保留 |
| TS 优化 | `Opt=(TS,CalcFC,NoEigenTest)` | 之后必须 Freq + IRC |
| IRC | `IRC=(RCFC,MaxPoints=50,MaxCycle=50)`，失败时 `LQA` | 先确认 TS 质量 |
| 热化学 | qRRHO（截断 100 cm⁻¹）+ 1 M | GoodVibes/Shermo |
| 能量组合 | 高水平单点 + 低水平热校正 | 所有物种同一方案 |
| 构象 | CREST 初筛 + DFT 重排序 | 窗口内所有构象都做 DFT |

---

## 常见问题与报错

本节补充第 8 节完整流程中容易遇到的问题（通用报错见第 3–7 节）。

**1. TS 优化正常结束，但频率显示 0 个虚频**

- 原因：初猜离 TS 太远，优化滑到了极小点（通常是产物或复合物）。
- 处理：检查结构；用柔性扫描或 QST2 重新获得初猜；确认 route 中写了 `TS`。

**2. TS 有 2 个虚频**

- 原因：停在了二阶鞍点，或第二个虚频是很小的甲基转动等模式。
- 处理：播放两个虚频；沿非反应坐标的那个虚频方向位移结构后重新做 TS 优化；收紧 `Opt=Tight`、`Int=SuperFine`。

**3. IRC 两端连接的不是预期的物种**

- 原因：找到的是另一个反应（或另一构象转换）的 TS。
- 处理：这本身就是重要信息：说明你的反应可能经由其他路径；重新考虑 TS 初猜。

**4. GoodVibes 提示找不到单点能文件，或单点能没有被使用**

- 原因：单点 log 的文件名不是 `X_TZ.log`，或者单点作业没有正常结束。
- 处理：检查文件名配对；`--spc` 的参数与文件名后缀一致。

**5. 能垒为负或明显过低**

- 原因：反应物参考用错（用了 IRC 端点的复合物、非最低构象）；标准态校正加了两次；把 Hartree 当成 kcal/mol；TS 与反应物的方案不一致。
- 处理：按第 9 节清单逐项检查。

**6. 不同泛函给出的能垒相差 5 kcal/mol 以上**

- 原因：体系对泛函敏感（色散、自相互作用误差、多参考特征）。
- 处理：在小模型上与高精度方法（如 DLPNO-CCSD(T)，需要其他程序）比较，选择表现最好的泛函，并在文中说明。

---

## 练习

1. 运行 `18-error-cases` 中的所有错误文件，记录每个作业出错的 Link 和真正的错误信息，与第 7 节速查表对照。
2. 把 `err_l301_wrong_multiplicity_fixed.gjf` 改为 `ROB3LYP`，比较 `<S**2>` 和能量的变化。
3. 完成 Diels–Alder 案例，报告 ΔG‡、ΔG<sub>r</sub> 和 ΔH‡；对比加与不加 1 M 校正、qRRHO 与 RRHO 的差别。
4. 把单点泛函换成 DSD-PBEP86/def2-TZVP（`DSDPBEP86`），比较能垒变化。
5. 把乙烯换成丙烯腈（CH₂=CHCN），用 CREST 和 DFT 找出 endo 与 exo 两个 TS，计算选择性并用 GoodVibes 的 `--label` 选项做 Boltzmann 分析。
6. 由你计算的 ΔH‡ 换算出 E<sub>a</sub>，与实验值约 27.5 kcal/mol 比较，并讨论气相与溶液、计算误差各自的贡献。

---

## 参考资料

**报错大全与社区资源**

- [Solutions for Gaussian Error Messages（网页版）](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/) 与 [GitHub 仓库 liyuanhe211/Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)：按 Link 整理的报错大全（CC BY-NC-SA 4.0，本章只做链接引用）
- [计算化学公社 Gaussian 报错原帖](http://bbs.keinsci.com/thread-4829-1-1.html)
- [YuOZW/gaussian-tips](https://github.com/YuOZW/gaussian-tips)：日文 Gaussian 技巧合集，含 IRC 与报错部分

**Gaussian 官方文档**

- [Keyword List](https://gaussian.com/keywords/)、[Input 格式](https://gaussian.com/input/)、[SCF](https://gaussian.com/scf/)、[Guess](https://gaussian.com/guess/)、[Integral](https://gaussian.com/integral/)
- [Opt](https://gaussian.com/opt/)、[IRC](https://gaussian.com/irc/)、[Freq](https://gaussian.com/freq/)、[Stable](https://gaussian.com/stable/)
- [DFT 方法与 EmpiricalDispersion 参数表](https://gaussian.com/dft/)、[Gen/GenECP](https://gaussian.com/gen/)、[Pseudo](https://gaussian.com/pseudo/)
- [SCRF（SMD）](https://gaussian.com/scrf/)、[Thermochemistry in Gaussian](https://gaussian.com/thermo/)、[Utilities（testrt 等）](https://gaussian.com/utils/)

**工具文档**

- [CREST 文档](https://crest-lab.github.io/crest-docs/)、[GoodVibes 文档](https://goodvibespy.readthedocs.io/)、[Shermo](http://sobereva.com/soft/shermo/)、[Basis Set Exchange](https://www.basissetexchange.org/)

**相关 GitHub 仓库**（仅链接，未复制代码）

- [crest-lab/crest](https://github.com/crest-lab/crest)：GFN-xTB 构象搜索，DFT 之前的标准步骤
- [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)：qRRHO、1 M 标准态、单点组合、Boltzmann 加权、反应剖面（v4.4.0）
- [jvalegre/aqme](https://github.com/jvalegre/aqme)：从 SMILES 出发的构象搜索、输入生成、自动检查虚频与报错并重算（QCORR）
- [duartegroup/autodE](https://github.com/duartegroup/autodE)：从 SMILES 自动生成反应剖面（含 TS 搜索）
- [patonlab/pyQRC](https://github.com/patonlab/pyQRC)：沿虚频方向位移，消除残留小虚频或做“廉价 IRC”
- [jvalegre/IRC_to_com](https://github.com/jvalegre/IRC_to_com)：把 IRC 路径点（含端点）转成新的 Gaussian 输入
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：中文项目，含错误诊断与结果验证清单
- [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash)：含调用 `testrt` 检查 route 的工具
- [chauncey-garrett/gaussian-tools](https://github.com/chauncey-garrett/gaussian-tools)：逐项提取虚频、自旋污染、热校正等的脚本清单
- [lenhanpham/Gaussian-Extractor](https://github.com/lenhanpham/Gaussian-Extractor)：批量检查作业状态、汇总能量
- [lenhanpham/OpenThermo](https://github.com/lenhanpham/OpenThermo)：C++ 热化学程序，支持 qRRHO
- [kimariyb/easy-Shermo](https://github.com/kimariyb/easy-Shermo)：批量调用 Shermo，中文说明
- [patonlab/pyDFTD3](https://github.com/patonlab/pyDFTD3)：单独计算 D3(0)/D3(BJ) 色散能

**关键文献**

- S. Grimme, Supramolecular binding thermodynamics by dispersion-corrected density functional theory, *Chem. Eur. J.* **2012**, 18, 9955.（qRRHO）DOI: [10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)
- G. Luchini, J. V. Alegre-Requena, I. Funes-Ardoiz, R. S. Paton, GoodVibes: automated thermochemistry for heterogeneous computational chemistry data, *F1000Research* **2020**, 9, 291. DOI: [10.12688/f1000research.22758.1](https://doi.org/10.12688/f1000research.22758.1)
- P. Pracht, F. Bohle, S. Grimme, Automated exploration of the low-energy chemical space with fast quantum chemical methods, *Phys. Chem. Chem. Phys.* **2020**, 22, 7169. DOI: [10.1039/C9CP06869D](https://doi.org/10.1039/C9CP06869D)
- P. Pracht et al., CREST—A program for the exploration of low-energy molecular chemical space, *J. Chem. Phys.* **2024**, 160, 114110. DOI: [10.1063/5.0197592](https://doi.org/10.1063/5.0197592)
- A. V. Marenich, C. J. Cramer, D. G. Truhlar, Universal solvation model based on solute electron density and on a continuum model of the solvent defined by the bulk dielectric constant and atomic surface tensions, *J. Phys. Chem. B* **2009**, 113, 6378.（SMD）DOI: [10.1021/jp810292n](https://doi.org/10.1021/jp810292n)
- S. Grimme, S. Ehrlich, L. Goerigk, Effect of the damping function in dispersion corrected density functional theory, *J. Comput. Chem.* **2011**, 32, 1456.（D3(BJ)）DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- J.-D. Chai, M. Head-Gordon, Long-range corrected hybrid density functionals with damped atom–atom dispersion corrections, *Phys. Chem. Chem. Phys.* **2008**, 10, 6615.（ωB97X-D）DOI: [10.1039/b810189b](https://doi.org/10.1039/b810189b)
- F. Weigend, R. Ahlrichs, Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for H to Rn, *Phys. Chem. Chem. Phys.* **2005**, 7, 3297.（def2）DOI: [10.1039/b508541a](https://doi.org/10.1039/b508541a)
- V. Guner, K. S. Khuong, A. G. Leach, P. S. Lee, M. D. Bartberger, K. N. Houk, A standard set of pericyclic reactions of hydrocarbons for the benchmarking of computational methods, *J. Phys. Chem. A* **2003**, 107, 11445. DOI: [10.1021/jp035501w](https://doi.org/10.1021/jp035501w)
- D. Rowley, H. Steiner, Kinetics of diene reactions at high temperatures, *Discuss. Faraday Soc.* **1951**, 10, 198.（丁二烯 + 乙烯实验动力学）DOI: [10.1039/DF9511000198](https://doi.org/10.1039/DF9511000198)
- T. Lu, Q. Chen, Shermo: A general code for calculating molecular thermochemistry properties, *Comput. Theor. Chem.* **2021**, 1200, 113249. DOI: [10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)

---

[← 上一章](17-hpc-jobs.md) | [返回目录](../README.md)
