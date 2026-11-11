# 第4章 单点能计算与 SCF 收敛技巧（Single-Point Energy & SCF Convergence in Gaussian 16）

本章介绍 Gaussian 16 单点能（single-point energy, SP）计算的写法，以及 `SCF Done` 能量、Hartree 与 kcal/mol、kJ/mol、eV 之间的换算。重点是 SCF 不收敛（L502 Convergence failure）时怎样系统排查：`SCF=XQC/QC/YQC`、`SCF=VShift`、`SCF=Fermi`、`SCF=NoIncFock`、`Guess=Read/Fragment`、`Int=SuperFine` 和 `Stable` 稳定性检查各自适合什么情形。四个可直接运行的示例覆盖过渡金属配合物、开壳层体系和金属团簇。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 写出一个规范的单点能输入，从输出中准确取出能量并换算成 kcal/mol 或 kJ/mol；
- 用 `#p` 读懂 SCF 迭代输出，区分“正在收敛”“振荡”“停滞”三种情形；
- 按照“先查输入、再换初猜、后换算法”的顺序排查 L502 SCF 不收敛；
- 知道 SCF 收敛不等于得到正确的电子态，会用 `Stable` 和 `<S**2>` 检查；
- 理解积分格点、对称性和基组线性相关对能量与收敛的影响。

单点能计算的典型用途：

| 用途 | 说明 |
|---|---|
| 高精度单点（high//low） | 在低水平（如 B3LYP-D3(BJ)/def2-SVP）优化好的结构上，用更好的泛函或更大的基组算能量，例如 ωB97X-D/def2-TZVP//B3LYP-D3(BJ)/def2-SVP |
| 反应能、构象能排序 | ΔE = ΣE(产物) − ΣE(反应物)，相同水平下比较 |
| 生成波函数用于分析 | 输出 chk/fchk，交给 cubegen、Multiwfn 做轨道、ESP、键级分析（见[第16章](16-postprocess-multiwfn.md)） |
| 检查电子结构 | 比较不同自旋多重度的能量，做波函数稳定性分析 |
| 溶剂化能、BSSE 等组合计算 | 见[第9章](09-solvation.md)、[第14章](14-weak-bsse.md) |

> [!NOTE]
> 本章默认你已经会写基本输入（[第2章](02-input-anatomy.md)），并了解泛函与基组的选择（[第3章](03-method-basis.md)）。全部示例输入在 [../examples/04-sp-scf/README.md](../examples/04-sp-scf/README.md)。

---

## 2. 原理速览

### 2.1 SCF 在做什么

Hartree-Fock 和 Kohn-Sham DFT 都要求解 Roothaan-Hall 型方程：

$$
\mathbf{F}(\mathbf{P})\,\mathbf{C} = \mathbf{S}\,\mathbf{C}\,\boldsymbol{\varepsilon},
\qquad
P_{\mu\nu} = \sum_{i}^{\mathrm{occ}} C_{\mu i} C_{\nu i}
$$

Fock（或 Kohn-Sham）矩阵 $\mathbf{F}$ 依赖于密度矩阵 $\mathbf{P}$，而 $\mathbf{P}$ 又由 $\mathbf{F}$ 的本征矢决定，所以只能迭代求解：猜一个 $\mathbf{P}^{(0)}$ → 构造 $\mathbf{F}$ → 对角化得到新轨道 → 新 $\mathbf{P}^{(1)}$ → …，直到前后两轮的密度矩阵几乎不变，这就是自洽场（self-consistent field, SCF）。

单纯迭代往往振荡甚至发散，Gaussian 16 默认用 **EDIIS + CDIIS 外推**来加速（官方 SCF 文档：默认不开阻尼和 Fermi 展宽）。DIIS 利用前若干轮的误差向量 $\mathbf{e} = \mathbf{FPS} - \mathbf{SPF}$ 线性组合出更好的 Fock 矩阵，EDIIS 在远离收敛时更稳。

### 2.2 为什么会不收敛

| 原因 | 典型体系 | 表现 |
|---|---|---|
| HOMO-LUMO 能隙小、近简并 | 过渡金属、金属团簇、双自由基、拉长的键 | 能量上下振荡，占据轨道来回翻转 |
| 初猜太差 | 重元素、带电体系、非常规成键 | 前几轮能量变化很大，后面振荡 |
| 输入本身有错 | 电荷/多重度写错、结构不合理（原子重叠、单位错） | 怎么调都不收敛，或收敛到离谱的能量 |
| 数值噪声 | 弥散基组导致线性相关、格点太粗 | 收敛到 1e-6 附近就停滞 |
| 对称性约束 | 高对称结构但电子态简并 | 被迫保持初猜的不可约表示占据，可能振荡或收敛到激发态 |

### 2.3 Gaussian 16 的默认设置

| 项目 | G16 默认 | 说明 |
|---|---|---|
| 初猜 | `Guess=Harris` | 对 HF 和所有 DFT 方法都是默认（官方 Guess 文档） |
| 算法 | EDIIS + CDIIS | 无阻尼、无 Fermi 展宽、无能级移动 |
| 收敛判据 | `SCF=Tight`，即 `Conver=8` | 密度矩阵 RMS 变化 < 1e-8 **且**最大变化 < 1e-6 |
| 最大迭代 | `MaxCycle=64`（`SCF=QC`/`DM` 为 512） | 以官方文档为准；实际日志里 `within N cycles` 的 N 会随方法变化，DFT 作业常见 128，以自己的输出为准 |
| 积分方式 | Direct SCF + 增量 Fock（`IncFock`） + `VarAcc` | 前期用较低积分精度，后期切换为全精度 |
| DFT 积分格点 | `Int=UltraFine`，即 pruned (99,590) | G09 默认是 FineGrid (75,302) |
| 双电子积分精度 | `Acc2E=12`（1e-12） | |

官方文档还指出：收敛**只看密度矩阵变化，不看能量变化**；RMS 密度变化为 $10^{-N}$ 时，能量误差通常约为 $10^{-2N}$ Hartree。所以默认的 Conver=8 对能量来说已经非常严格，**不需要**为了“更准”再收紧。

---

## 3. 第一个单点能计算：H2O

### 3.1 输入文件

单点能是 Gaussian 的**默认任务类型**：route 里不写 `Opt`、`Freq` 等任务关键词，就是单点（也可以显式写 `SP`）。示例文件 [h2o_sp_b3lyp_def2tzvp.gjf](../examples/04-sp-scf/h2o_sp_b3lyp_def2tzvp.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=h2o_sp_b3lyp_def2tzvp.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ SP

H2O single point, B3LYP-D3(BJ)/def2-TZVP

0 1
O     0.000000     0.000000     0.117300
H     0.000000     0.757200    -0.469200
H     0.000000    -0.757200    -0.469200

```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `%nprocshared` | 与申请的核数一致 | 共享内存并行核数；也可以在作业脚本中用 `-c`/`GAUSS_CDEF` 指定，见[第17章](17-hpc-jobs.md) |
| `%mem` | 申请内存的 75–85% | G16 不写时默认只有 800 MB（100 MW） |
| `%chk` | 与输入同名 | 保存波函数，后续 `Guess=Read`、formchk 都要用 |
| `#p` | 推荐 | 详细打印，能看到每一轮 SCF 迭代 |
| `B3LYP/def2TZVP` | — | 泛函/基组；Gaussian 关键词里 def2 系列不写连字符 |
| `EmpiricalDispersion=GD3BJ` | B3LYP、PBE0 等可用 | D3(BJ) 色散校正，Minnesota 泛函只能用 GD3，见[第3章](03-method-basis.md) |
| `SP` | 可省略 | 单点能任务（默认） |

> [!IMPORTANT]
> 输入文件**必须以空行结尾**，title 行前后都要有空行。缺少空行是 L101 "End of file in ZSymb" 的常见原因。

### 3.2 运行

```bash
g16 < h2o_sp_b3lyp_def2tzvp.gjf > h2o_sp_b3lyp_def2tzvp.log
# or equivalently
g16 h2o_sp_b3lyp_def2tzvp.gjf          # writes h2o_sp_b3lyp_def2tzvp.log

# Did it finish?
tail -n 1 h2o_sp_b3lyp_def2tzvp.log    # expect: Normal termination of Gaussian 16 ...
```

在集群上用示例目录中的 Slurm 脚本：`sbatch run_g16.sh h2o_sp_b3lyp_def2tzvp.gjf`（脚本说明见示例 README，通用作业调度见 ../../Common/tutorials/02-job-scheduler.md）。

### 3.3 找到能量

```bash
grep "SCF Done" h2o_sp_b3lyp_def2tzvp.log
```

（示意输出，数值仅供参考）

```text
 SCF Done:  E(RB3LYP) =  -76.4660xxxxxx     A.U. after   10 cycles
```

这一行的含义：

- `E(RB3LYP)`：方法标签，R 表示闭壳层限制性计算，开壳层会显示 `E(UB3LYP)`；
- 单位 `A.U.`，即 Hartree（Eh）；
- `after 10 cycles`：迭代了 10 轮收敛；
- 下一行的 `Conv=` 是最终达到的密度收敛度，`-V/T` 是维里比，正常体系接近 2。

> [!NOTE]
> 加了 `EmpiricalDispersion` 时，Gaussian 把 D3 色散能并入核排斥能（输出里有 `Nuclear repulsion after empirical dispersion term` 一行），所以 `SCF Done` 的能量**已经包含色散校正**，不需要再手动加。

不同方法的“最终能量”不在同一行：

| 方法 | 去哪里找最终能量 |
|---|---|
| HF、DFT（含 D3） | `SCF Done:  E(...) =` |
| 双杂化（B2PLYPD3、DSDPBEP86 等） | `SCF Done` 只是参考行列式的能量，**还要加上 MP2 部分**，最终能量看其后的 `E2(...)`/`E(...)` 行 |
| MP2 | `EUMP2 =` |
| CCSD(T) | `CCSD(T)=` |
| 任意方法 | 文件末尾 archive 段（以 `1\1\` 开头、`@` 结尾）里的 `HF=`、`MP2=`、`CCSD(T)=` 等字段 |

批量提取多个 log 的能量，可以用本仓库的 [g16_log_summary.py](../scripts/g16_log_summary.py)（`python ../scripts/g16_log_summary.py -h` 查看用法），或用 cclib 交叉验证。

### 3.4 能量单位换算

Gaussian 输出的能量单位是 Hartree，相对能量通常换算成 kcal/mol 或 kJ/mol：

| 1 Hartree 等于 | 数值 |
|---|---|
| kcal/mol | 627.5095 |
| kJ/mol | 2625.500 |
| eV | 27.211386 |
| cm⁻¹ | 219474.63 |

例如两个构象 A、B 的能量差：

```python
# relative energy of conformer B vs A (values are illustrative)
HARTREE2KCAL = 627.5095
E_A = -232.123456   # Eh, from "SCF Done"
E_B = -232.120987   # Eh
print(f"dE = {(E_B - E_A) * HARTREE2KCAL:.2f} kcal/mol")   # -> 1.55 kcal/mol
```

> [!WARNING]
> 只有**同一方法、同一基组、同一溶剂模型、同一格点、同一程序版本**算出的能量才能直接相减。G09（FineGrid）和 G16（UltraFine）、不同 SCRF 默认设置之间的差别可能达到零点几 kcal/mol，足以改变构象排序。

---

## 4. 读懂 SCF 迭代输出（#p）

### 4.1 收敛判据行

`#p` 输出中，SCF 开始前会打印本次使用的判据（示意输出，数值仅供参考）：

```text
 Requested convergence on RMS density matrix=1.00D-08 within 128 cycles.
 Requested convergence on MAX density matrix=1.00D-06.
 Requested convergence on             energy=1.00D-06.
```

前两行就是 `Conver=8` 对应的两个判据；`within N cycles` 给出本次作业的最大迭代数。

### 4.2 每一轮迭代

（示意输出，数值仅供参考，字段与 G16 #p 输出一致，中间若干行省略）

```text
 Cycle   1  Pass 1  IDiag  1:
 E= -76.3958412345678
 DIIS: error= 4.1D-02 at cycle   1 NSaved=   1.
 ...
 Gap=     0.312 Goal=   None    Shift=    0.000
 RMSDP=9.8D-03 MaxDP=1.2D-01              OVMax= 1.1D-01

 Cycle   2  Pass 1  IDiag  1:
 E= -76.4421987654321     Delta-E=       -0.046357530864 Rises=F Damp=F
 DIIS: error= 1.6D-02 at cycle   2 NSaved=   2.
 ...
 RMSDP=3.1D-03 MaxDP=4.5D-02 DE=-4.64D-02 OVMax= 4.0D-02
```

| 字段 | 含义 | 怎么看 |
|---|---|---|
| `E=` | 本轮总能量 | 正常情况下单调下降，后期只在小数点后很多位变化 |
| `Delta-E=` | 与上一轮的能量差 | 符号来回变、幅度不减小，说明在振荡 |
| `DIIS: error=` | DIIS 误差向量的大小 | 应逐渐减小 |
| `Gap=` | 当前 HOMO-LUMO 能隙（Hartree） | 很小（< 0.05 Eh）时容易振荡；`Shift=` 显示当前的能级移动量 |
| `RMSDP` / `MaxDP` | 密度矩阵 RMS / 最大变化 | 收敛判据就是它们分别小于 1e-8 / 1e-6 |
| `Rises=T` | 能量比上一轮升高 | 偶尔一次没关系，频繁出现说明 DIIS 不稳定 |

监控一个正在运行的作业：

```bash
# energy of every SCF cycle (current or last SCF)
grep -E "^ E= " job.log | tail -n 20
# density changes
grep "RMSDP=" job.log | tail -n 10
# HOMO-LUMO gap and level shift during SCF
grep "Gap=" job.log | tail -n 5
```

### 4.3 三种典型迭代行为

| 现象 | 判断 | 对策 |
|---|---|---|
| 能量单调下降，RMSDP 稳定减小，只是没在 MaxCycle 内收敛 | **缓慢收敛** | 可以适当增加 `MaxCycle`，或者换 `SCF=XQC` |
| 能量在几个值之间来回跳，`Delta-E` 正负交替且幅度不减 | **振荡** | 加 `MaxCycle` 没用；换初猜，或用 `VShift`、`Fermi`、`QC` |
| RMSDP 降到 1e-6~1e-7 后不再下降 | **数值噪声导致停滞** | `SCF=NoIncFock`、`Int=SuperFine`，检查弥散基组的线性相关 |

### 4.4 失败时的输出

SCF 达到最大迭代仍未收敛时（示意输出，数值仅供参考）：

```text
 >>>>>>>>>> Convergence criterion not met.
 SCF Done:  E(UB3LYP) =  -1724.81xxxxxx     A.U. after  129 cycles
            NFock=128  Conv=0.21D-04     -V/T= 2.0021
 ...
 Convergence failure -- run terminated.
 Error termination via Lnk1e in /opt/g16/l502.exe at ...
```

注意这里**也有 SCF Done 行**，但前面有 `Convergence criterion not met`，这个能量**不能用**。写批处理脚本时要先检查 `Normal termination`，再取能量。g16_error_doctor.py 可以自动识别这类错误并给出中文建议（`python ../scripts/g16_error_doctor.py -h`）。

---

## 5. SCF 不收敛的系统排查流程

> [!IMPORTANT]
> 最常见的错误做法是一失败就把 `MaxCycle` 改成 512 或 1000。如果迭代在振荡，再多迭代也只是浪费机时。正确顺序是：**先确认输入没错 → 再换更好的初猜 → 最后才换算法、调参数**。

### 第 0 步：检查输入本身

1. **电荷与多重度**：电子数 = 各原子序数之和 − 总电荷（用 ECP 时 Gaussian 会自动扣除芯电子，但奇偶性不变）。电子数为奇数时多重度只能是 2、4、6…，为偶数时只能是 1、3、5…。写错时 L301 会直接报 `The combination of multiplicity 1 and 17 electrons is impossible.`。
2. **自旋态是否合理**：例如 Fe(II) 高自旋是五重态，低自旋是单重态；Cu(II) 是二重态；O2 基态是三重态。
3. **几何结构**：在 GaussView 中检查有没有原子重叠、键长异常（Bohr 与 Å 混用、从晶体结构复制时少了原子），以及多余或缺失的 H。
4. **方法与基组**：重元素是否带了 ECP（[第3章](03-method-basis.md)）；大分子用 aug- 系列弥散基组时容易出现线性相关。

### 第 1 步：看迭代行为

用 `#p` 重跑几轮，或者直接看已有日志，按 4.3 节判断是缓慢收敛、振荡还是停滞。

### 第 2 步：给一个更好的初猜（第 7、9 节）

- 先用小基组（def2-SVP、6-31G(d)）或 HF 收敛，再 `Guess=Read` 读入大基组计算；
- 先算对应的闭壳层离子（如阳离子或阴离子），再读入开壳层体系；
- 从附近已收敛的几何结构读入（扫描、优化的前一步）；
- 双自由基、反铁磁耦合体系用 `Guess=Mix` 或 `Guess=Fragment`（[第13章](13-openshell-stability.md)）。

### 第 3 步：换算法或加稳定化手段（第 6、8 节）

| 选项 | 作用 | 适用情形 | 代价与限制 |
|---|---|---|---|
| `SCF=XQC` | 先做常规 DIIS（默认最多 `MaxConventionalCycles=32` 轮），不收敛再自动切换到二次收敛 QC | **首选**的通用补救手段，容易收敛的体系没有额外开销 | QC 阶段每轮较慢 |
| `SCF=QC` | 全程用二次收敛算法（线搜索 + Newton-Raphson） | 振荡严重的体系 | 慢；**不支持 ROHF/ROKS** |
| `SCF=YQC` | 最速下降 → 缩放最速下降 → 常规 SCF，失败时才用 QC | 很大的分子 | — |
| `SCF=VShift[=N]` | 虚轨道能量上移 N×0.001 Eh（不写 N 时为 100，即 0.1 Eh） | 能隙小导致占据轨道翻转、振荡 | 收敛变慢；会关闭自动 archive |
| `SCF=Fermi` | 前期迭代使用温度展宽（分数占据），并隐含阻尼和能级移动 | 金属团簇、近简并体系 | — |
| `SCF=NoIncFock` | 每轮完整重建 Fock 矩阵，不用增量 Fock | 后期停滞在 1e-6~1e-7 | 每轮稍慢 |
| `SCF=(MaxCycle=N)` | 增加最大迭代数 | **只有**能量稳定下降时才有意义 | — |
| `SCF=CDIIS` / `SCF=Damp` | 只用 CDIIS（隐含阻尼）/ 前期动态阻尼 | 前几轮剧烈振荡 | 官方提示阻尼与 EDIIS 配合不好 |
| `SCF=NoVarAcc` | 全程使用全精度积分 | 怀疑前期低精度积分引起问题 | 稍慢 |

多个选项可以组合，例如 `SCF=(XQC,VShift=300,MaxCycle=300)`。

> [!TIP]
> 对 ROHF/ROKS 这类不能用 QC 的情况，官方 SCF 文档的建议是在 route 中加 `Use=L506`。

### 第 4 步：收敛后检查是不是正确的电子态（第 10 节）

- 开壳层计算看 `<S**2>` 是否接近 S(S+1)；
- 用 `Stable` 或 `Stable=Opt` 检查波函数稳定性；
- 对过渡金属体系，比较不同多重度的能量。

---

## 6. 实战 1：过渡金属配合物 [Fe(H2O)6]2+ 与 SCF=XQC

### 6.1 体系分析

六水合铁(II)：Fe²⁺ 是 d⁶，弱场配体水使其处于**高自旋**态（t2g⁴eg²，4 个未成对电子），S = 2，多重度 2S+1 = 5。电子数 = 26 + 6×10 − 2 = 84（偶数），与五重态一致。

难点：

- d 轨道近简并，HOMO-LUMO 能隙小，常规 DIIS 容易振荡；
- 起始结构按理想的 Th 对称性搭建，而高自旋 d⁶ 在高对称下是轨道简并态（会发生 Jahn-Teller 畸变），对称性约束可能让 SCF 停在错误的占据上，所以加 `NoSymm`。

### 6.2 输入文件

[fe_complex_sp_scf_xqc.gjf](../examples/04-sp-scf/fe_complex_sp_scf_xqc.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=fe_complex_sp_scf_xqc.chk
#p UB3LYP/def2TZVP EmpiricalDispersion=GD3BJ SCF=(XQC,MaxConventionalCycles=64) NoSymm

[Fe(H2O)6]2+ high-spin quintet, single point with SCF=XQC

2 5
Fe    0.000000    0.000000    0.000000
O     2.130000    0.000000    0.000000
H     2.713761    0.774676    0.000000
H     2.713761   -0.774676    0.000000
O    -2.130000    0.000000    0.000000
H    -2.713761    0.774676    0.000000
H    -2.713761   -0.774676    0.000000
O     0.000000    2.130000    0.000000
H     0.000000    2.713761    0.774676
H     0.000000    2.713761   -0.774676
O     0.000000   -2.130000    0.000000
H     0.000000   -2.713761    0.774676
H     0.000000   -2.713761   -0.774676
O     0.000000    0.000000    2.130000
H     0.774676    0.000000    2.713761
H    -0.774676    0.000000    2.713761
O     0.000000    0.000000   -2.130000
H     0.774676    0.000000   -2.713761
H    -0.774676    0.000000   -2.713761

```

| 关键词 | 说明 |
|---|---|
| `UB3LYP` | 开壳层用非限制性 Kohn-Sham；多重度大于 1 时不写 U，Gaussian 也会自动用 U，但显式写出更清楚 |
| `def2TZVP` | Fe 原子序数 < 37，def2-TZVP 对 Fe 是全电子基组，不需要 ECP |
| `SCF=XQC` | 常规 SCF 失败后自动切换到 QC |
| `MaxConventionalCycles=64` | 把切换前允许的常规迭代数从默认的 32 放宽到 64：常规 SCF 多数情况下更快，给它多一点机会 |
| `NoSymm` | 不使用对称性，避免轨道对称性约束导致占据错误；同时不再把分子转到标准取向 |
| `2 5` | 总电荷 +2，多重度 5 |

### 6.3 运行与检查

```bash
g16 < fe_complex_sp_scf_xqc.gjf > fe_complex_sp_scf_xqc.log

grep "Normal termination" fe_complex_sp_scf_xqc.log
grep "SCF Done" fe_complex_sp_scf_xqc.log
grep "S\*\*2 before annihilation" fe_complex_sp_scf_xqc.log
# did XQC switch to the quadratic algorithm? (l508 is the QC SCF link)
grep -c "l508" fe_complex_sp_scf_xqc.log
```

判断标准：

1. 有 `Normal termination`，`SCF Done` 前没有 `Convergence criterion not met`；
2. `<S**2>` 应接近 S(S+1) = 6.00（示意：`S**2 before annihilation 6.0040, after 6.0000`，数值仅供参考）。明显偏大说明混入了更高自旋态；
3. 如果常规 SCF 在 64 轮内收敛，日志里只有 L502 的迭代；如果切换到了 QC，会看到 L508 打印的二次收敛迭代。

> [!WARNING]
> 自旋态能量差（高自旋与低自旋）对泛函中 HF 交换的比例非常敏感：HF 交换越多越偏向高自旋，纯 GGA 往往偏向低自旋。研究自旋交叉时应该比较多种泛函，并参考实验或高精度基准，不能只看一个 B3LYP 结果。

---

## 7. 实战 2：用小基组初猜读入大基组（Guess=Read + Link1）

### 7.1 为什么有效

大基组的 SCF 更难收敛，每轮也更贵。用小基组收敛的轨道作为初猜，Gaussian 会把它**投影到新基组**上，通常只需要很少几轮就能收敛。同样的思路还可以用于：

- 先算闭壳层离子，再读入开壳层自由基；
- 先用 HF 或纯泛函收敛，再读入杂化泛函；
- 读入相邻几何结构（扫描、优化的前一点）的波函数。

### 7.2 输入文件

[guess_read_small_to_big_basis.gjf](../examples/04-sp-scf/guess_read_small_to_big_basis.gjf)，体系是平面四方形的 [Cu(NH3)4]²⁺（Cu(II) d⁹，电子数 29 + 4×10 − 2 = 67，二重态）：

```text
%nprocshared=8
%mem=12GB
%chk=cu_nh3_4_guess_read.chk
#p UB3LYP/def2SVP EmpiricalDispersion=GD3BJ

[Cu(NH3)4]2+ step 1: small basis (def2-SVP) to generate a guess

2 2
Cu    0.000000    0.000000    0.000000
N     2.030000    0.000000    0.000000
H     2.395535    0.000000    0.952252
H     2.395535   -0.824674   -0.476126
H     2.395535    0.824674   -0.476126
N     0.000000    2.030000    0.000000
H     0.000000    2.395535    0.952252
H     0.824674    2.395535   -0.476126
H    -0.824674    2.395535   -0.476126
N    -2.030000    0.000000    0.000000
H    -2.395535    0.000000    0.952252
H    -2.395535    0.824674   -0.476126
H    -2.395535   -0.824674   -0.476126
N     0.000000   -2.030000    0.000000
H     0.000000   -2.395535    0.952252
H    -0.824674   -2.395535   -0.476126
H     0.824674   -2.395535   -0.476126

--Link1--
%nprocshared=8
%mem=12GB
%chk=cu_nh3_4_guess_read.chk
#p UB3LYP/def2TZVP EmpiricalDispersion=GD3BJ Geom=Check Guess=Read

[Cu(NH3)4]2+ step 2: large basis (def2-TZVP), geometry and guess from chk

2 2

```

要点：

| 要点 | 说明 |
|---|---|
| `--Link1--` | 把两个作业串成一个输入，按顺序执行 |
| 两步都写 `%chk`，而且同名 | 第二步从这个 chk 读几何结构和轨道，结束后再把新波函数写回去 |
| `Geom=Check` | 从 chk 读几何结构；**仍然需要** title 段和电荷多重度行 |
| `Guess=Read` | 从 chk 读轨道作为初猜，基组不同时自动投影 |
| `Geom=AllCheck` | 连 title 和电荷多重度也从 chk 读，此时这两段都**不能**再写 |

如果想保留第一步的 chk 不被覆盖，第二步改成：

```text
%oldchk=cu_nh3_4_guess_read.chk
%chk=cu_nh3_4_tzvp.chk
```

`%oldchk` 只读，Gaussian 会把其中的内容复制到新的 `%chk` 后再开始计算。

### 7.3 检查

```bash
grep -E "Initial guess from the checkpoint file|SCF Done" guess_read_small_to_big_basis.log
grep "S\*\*2 before annihilation" guess_read_small_to_big_basis.log
```

第二步应出现 `Initial guess from the checkpoint file`，而且 `after N cycles` 应明显少于第一步从 Harris 初猜开始的轮数。Cu(II) 二重态的 `<S**2>` 应接近 0.75。

> [!TIP]
> 读入初猜时**电荷和多重度可以与 chk 中不同**：Gaussian 会按新的电子数重新填充轨道。这正是“先算闭壳层阳离子，再读入中性自由基”能行得通的原因。

---

## 8. 实战 3：小能隙开壳层体系与 SCF=VShift

### 8.1 能级移动的原理

level shifting（Saunders-Hillier）给虚轨道的能量人为加上一个正值 $b$：

$$
\mathbf{F}' = \mathbf{F} + b\,\mathbf{S}\,(\mathbf{S}^{-1} - \mathbf{P}_{\mathrm{occ}})\,\mathbf{S}
$$

效果相当于把 HOMO-LUMO 能隙人为拉大，抑制占据轨道与虚轨道之间的剧烈混合和“轨道翻转”，让迭代变得平稳但更慢。**收敛后的解不受 shift 影响**，shift 只改变到达解的路径。

Gaussian 写法：`SCF=VShift=N`，移动量为 N×0.001 Hartree（N=500 即 0.5 Eh）。`#p` 输出中 `Shift=` 字段会显示实际使用的值。

### 8.2 输入文件

金属小团簇是典型的小能隙体系。[radical_sp_vshift.gjf](../examples/04-sp-scf/radical_sp_vshift.gjf) 以 Cu3 团簇（29×3 = 87 个电子，二重态自由基）为例：

```text
%nprocshared=8
%mem=12GB
%chk=cu3_radical_vshift.chk
#p UPBE1PBE/def2TZVP SCF=(VShift=500,MaxCycle=300)

Cu3 doublet radical cluster, single point with level shifting

0 2
Cu    0.000000    0.000000    0.000000
Cu    2.300000    0.000000    0.000000
Cu    1.150000    2.163331    0.000000

```

| 关键词 | 说明 |
|---|---|
| `UPBE1PBE` | PBE0 杂化泛函（Gaussian 关键词 PBE1PBE），开壳层 |
| `VShift=500` | 虚轨道上移 0.5 Eh；0.1–0.5 Eh 是常用范围，太大收敛会非常慢 |
| `MaxCycle=300` | level shifting 使收敛变慢，所以同时放宽最大迭代数 |

> [!NOTE]
> 起始结构是一个等腰三角形（Cu3 有 Jahn-Teller 效应，等边三角形不是能量最低构型）。这里只演示单点和 SCF 选项；真正研究团簇时需要先做结构优化和多个自旋态比较。

### 8.3 其他两个常用选项

```text
#p UPBE1PBE/def2TZVP SCF=(Fermi,MaxCycle=300)
```

`SCF=Fermi` 在前期迭代中用 Fermi-Dirac 展宽让轨道分数占据，平滑占据数的突变；按官方文档，Fermi 默认还**隐含阻尼和能级移动**。它对金属团簇、近简并体系往往比单纯 VShift 更有效。最终收敛时占据数会回到整数 0/1。

```text
#p UB3LYP/def2TZVP SCF=(XQC,NoIncFock) Int=SuperFine
```

`SCF=NoIncFock` 让每一轮都完整重建 Fock 矩阵。Direct SCF 默认用增量 Fock（只计算密度变化部分的贡献），速度快但会累积数值误差，表现为 RMSDP 降到 1e-6~1e-7 就不动了。这种“最后一公里”停滞，用 NoIncFock 配合更细的格点通常可以解决。

### 8.4 检查

```bash
grep "Shift=" radical_sp_vshift.log | head -n 3     # level shift actually applied
grep -E "SCF Done|S\*\*2 before" radical_sp_vshift.log
```

二重态的 `<S**2>` 理论值是 0.75。另外注意，使用 VShift 会关闭自动 archive，所以日志末尾可能没有 `1\1\` 开头的归档块，这属于正常现象，提取能量要用 `SCF Done` 行。

---

## 9. 更多初猜手段

| 写法 | 作用 | 何时使用 |
|---|---|---|
| `Guess=Harris` | Harris 泛函初猜（HF/DFT 的默认值） | 通常不用写 |
| `Guess=Huckel` | 扩展 Hückel 初猜 | Harris 初猜给出明显错误的占据时可以尝试 |
| `Guess=Core` | 对角化核哈密顿量 | 很差，几乎只用于测试 |
| `Guess=Read` | 从 `%chk`/`%oldchk` 读轨道，自动投影到当前基组 | 小基组→大基组，相邻几何结构，换泛函 |
| `Guess=Mix` | 混合 HOMO 与 LUMO，打破 α/β 对称性 | 开壳层单重态（双自由基），见[第13章](13-openshell-stability.md) |
| `Guess=Fragment=N` | 用 N 个碎片的波函数拼成初猜 | 反铁磁耦合、broken-symmetry 计算、弱相互作用复合物 |
| `Guess=Only` | 只生成并打印初猜，然后结束 | 配合 `Pop=Full` 检查初猜的轨道占据和对称性 |
| `Guess=Alter` | 交换指定的占据/虚轨道 | 想得到特定电子态，例如某个对称性的激发态 |
| `Guess=Always` | 优化的每一步都重新生成初猜 | 优化过程中波函数特征发生突变时 |

### 9.1 Guess=Fragment：拉长 H2 的 broken-symmetry 初猜

`Guess=Fragment=N` 的规则：原子后面用 `(Fragment=n)` 标出所属碎片；电荷多重度行先写**总体**，再依次写**每个碎片**；碎片多重度写成负数表示该碎片的未成对电子取 β 自旋。

最小示例 [h2_stretched_guess_fragment.gjf](../examples/04-sp-scf/h2_stretched_guess_fragment.gjf)：H–H 拉到 2.0 Å，两个 H 原子各带一个电子，一个取 α、一个取 β，组成 broken-symmetry 单重态：

```text
%nprocshared=8
%mem=12GB
%chk=h2_stretched_guess_fragment.chk
#p UB3LYP/def2SVP Guess=Fragment=2

Stretched H2 (2.0 A), broken-symmetry singlet from fragment guess

0 1 0 2 0 -2
H(Fragment=1)    0.000000    0.000000    0.000000
H(Fragment=2)    0.000000    0.000000    2.000000

```

`0 1 0 2 0 -2` 表示：总体中性单重态；碎片 1 中性二重态（α）；碎片 2 中性二重态，负号表示其未成对电子取 β 自旋。收敛后 `<S**2>` 接近 1（单重态与三重态的等量混合），这是 broken-symmetry 解的特征，不是“出错”。Broken-symmetry 方法和磁耦合常数 J 的计算见[第13章](13-openshell-stability.md)。

---

## 10. 收敛到错误态的风险：Stable 稳定性检查

SCF 收敛只说明找到了能量泛函的一个**驻点**，不保证是最低能量的解。常见的情况有：

- 拉长的键、双自由基：闭壳层 RHF/RKS 解存在 RHF→UHF 不稳定性，真正更低的是 UHF/UKS broken-symmetry 解；
- 过渡金属：收敛到了占据方式错误的激发组态；
- 高对称性约束：被迫保持初猜的轨道对称性。

`Stable` 关键词检查当前波函数是否是“稳定的局部极小”（官方文档：默认检查所有实数不稳定性但不重新优化；`Stable=Opt` 会在发现不稳定性时自动降低约束重新优化，必要时变成 UHF）。

示例 [h2_stretched_stable_opt.gjf](../examples/04-sp-scf/h2_stretched_stable_opt.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=h2_stretched_stable_opt.chk
#p B3LYP/def2SVP Stable=Opt

Stretched H2 (2.0 A): RB3LYP is unstable, Stable=Opt finds the UB3LYP solution

0 1
H    0.000000    0.000000    0.000000
H    0.000000    0.000000    2.000000

```

```bash
grep -E "wavefunction is stable|instability|SCF Done|S\*\*2 before" h2_stretched_stable_opt.log
```

预期：先得到 RB3LYP 能量，稳定性分析报告存在 RHF→UHF 不稳定性，随后重新优化得到更低的 UB3LYP 能量，`<S**2>` 接近 1；最后的稳定性检查报告波函数稳定。拿到稳定波函数后，后续计算（优化、频率、TD-DFT）用 `Guess=Read` 读入它。

> [!IMPORTANT]
> 官方 Freq 文档明确指出：解析频率依赖的 CPHF 只有在同一自旋多重度下不存在更低能量的波函数时才有物理意义。对开壳层、过渡金属或可能有双自由基特征的体系，**在优化和频率之前先做一次 Stable 检查**是好习惯。详细讨论见[第13章](13-openshell-stability.md)。

---

## 11. 数值精度：积分格点、对称性与线性相关

### 11.1 积分格点

DFT 的交换相关能是在原子中心格点上数值积分得到的。G16 的命名格点（官方 Integral 文档）：

| 关键词 | 格点 | 说明 |
|---|---|---|
| `Int=SG1` | pruned (50,194) | 不推荐用于正式计算 |
| `Int=FineGrid` | pruned (75,302) | G09 默认 |
| `Int=UltraFine` | pruned (99,590) | **G16 默认** |
| `Int=SuperFine` | pruned (175,974)，第二周期以后的元素为 (250,974) | 数值噪声问题、低频模式、Minnesota 泛函 |

需要 `Int=SuperFine` 的情形：

- SCF 停滞在 1e-6~1e-7、不同轮之间能量有 1e-6 量级的抖动；
- M06 系列等对格点敏感的泛函，特别是计算弱相互作用和低频振动时；
- 计算非常小的能量差（如构象能 < 0.5 kcal/mol）时，要确认格点不影响结论。

> [!WARNING]
> 格点是“方法定义”的一部分：同一项研究里，所有需要相减的能量都必须用同一套格点。不要用 SuperFine 的能量减 UltraFine 的能量。

### 11.2 对称性：NoSymm 的取舍

| 写法 | 效果 |
|---|---|
| 默认 | 自动识别点群，按标准取向重新放置分子，利用对称性加速，并**保持初猜的轨道对称性占据** |
| `SCF=NoSymm`（等价于 `Guess=NoSymm`） | 只解除轨道对称性约束，积分仍利用对称性 |
| `NoSymm`（即 `Symmetry=None`） | 完全不用对称性，也不旋转分子（保持输入取向） |

需要 `NoSymm` 的情形：电子态简并或对称性破缺（如 6.1 节）、需要保持输入取向（做 cube 差分、外加电场、按碎片分析）、ModRedundant 约束或冻结原子让对称性不再成立。代价是失去对称性带来的加速，对高对称的小分子比较明显。

### 11.3 基组线性相关

在大分子上用很多弥散函数（aug-cc-pVTZ、6-311++G(2df,2pd)）时，基函数之间可能接近线性相关，重叠矩阵出现极小的本征值，SCF 数值不稳定。Gaussian 会自动剔除这些近线性相关的组合。输出中对比两行：

```bash
grep -E "NBasis=|NBsUse=" job.log | head -n 4
```

如果 `NBsUse` 小于 `NBasis`，说明有基函数被剔除。剔除少量是正常的；数量很多且 SCF 难收敛时，考虑只给需要的原子（阴离子中心、氢键受体）加弥散函数，或换成 def2-TZVP 这类本身较“干净”的基组。

---

## 12. 单点能的典型用法：high//low 组合

最常见的做法是：低水平优化加频率（得到结构和热校正），再在其结构上做高水平单点。

```text
%nprocshared=8
%mem=12GB
%chk=mol.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

Step 1: optimization + frequencies at the low level

0 1
...coordinates...

--Link1--
%nprocshared=8
%mem=12GB
%chk=mol.chk
#p wB97XD/def2TZVP Geom=Check Guess=Read

Step 2: high-level single point on the optimized geometry

0 1

```

（坐标处为占位符；完整可运行的版本见[第6章](06-freq-thermo.md)的 `highlevel_sp_plus_thermo.gjf`。）

- 记号 ωB97X-D/def2-TZVP//B3LYP-D3(BJ)/def2-SVP：`//` 前面是单点方法，后面是优化方法；
- 不同泛函之间 `Guess=Read` 照样可用，只是初猜不如同一泛函时那么好；
- 自由能组合 G = E(高水平单点) + G_corr(低水平频率)，详见[第6章](06-freq-thermo.md)第11节；
- 溶液相单点要在单点步骤中加同样的 SCRF 设置，见[第9章](09-solvation.md)。

---

## 参数选择与收敛建议

| 情形 | 推荐设置 | 说明 |
|---|---|---|
| 普通有机分子 | 默认（不写 SCF 选项） | G16 默认 Tight 判据已足够 |
| 第一次遇到 L502 | 先查电荷/多重度/结构，再加 `SCF=XQC` | 最省事、对易收敛体系无副作用 |
| 过渡金属、开壳层 | `SCF=XQC`，必要时 `NoSymm`；收敛后 `Stable` 检查 | 并比较不同自旋多重度 |
| 大基组难收敛 | 小基组先收敛，再 `Geom=Check Guess=Read` | 7 节 Link1 写法 |
| 能量振荡 | `SCF=VShift=300~500` 或 `SCF=Fermi`，或 `SCF=QC` | 加 MaxCycle 无效 |
| 能量稳定下降但没收敛完 | `SCF=(MaxCycle=256~512)` | 只有这种情况才增加迭代数 |
| 停滞在 1e-6~1e-7 | `SCF=NoIncFock` + `Int=SuperFine` | 数值噪声 |
| 很大的分子 | `SCF=YQC` | 官方针对大分子设计的算法 |
| ROHF/ROKS 难收敛 | `Use=L506` | QC 不支持 RO |
| 双自由基、反铁磁 | `Guess=Mix` 或 `Guess=Fragment=N` + `Stable=Opt` | 见第13章 |
| 弥散基组、大分子 | 检查 `NBsUse` 与 `NBasis` | 必要时减少弥散函数 |
| 需要比较小能量差 | 全部用相同格点；必要时 `Int=SuperFine` 复核 | |

---

## 常见问题与报错

**1. `Convergence failure -- run terminated.`（L502）**

- 原因：在最大迭代数内没有达到收敛判据。
- 处理：按第 5 节的流程依次检查输入、换初猜（`Guess=Read`）、加 `SCF=XQC`、再考虑 `VShift`/`Fermi`/`QC`。不要直接把 MaxCycle 改成 1000。

**2. `The combination of multiplicity 1 and 17 electrons is impossible.`（L301）**

- 原因：电子数奇偶性与多重度不匹配（奇数个电子却写了单重态）。
- 处理：重新数电子，改电荷或多重度。本仓库的 [gjf_builder.py](../scripts/gjf_builder.py) 生成输入时会自动做这项检查。

**3. `Inaccurate quadrature in CalDSu`（L502）**

- 原因：数值积分得到的电子数与实际电子数相差过大，常见于结构不合理（原子过近）或格点相对基组太粗。
- 处理：先检查结构，再用 `Int=SuperFine`；如果是在优化过程中出现，检查前几步结构是否已经跑偏。

**4. SCF 收敛了，但能量比预期高很多（相差零点几 Hartree）**

- 原因：电子态错误（错误的多重度或占据），或者电荷写错。
- 处理：核对电荷和多重度；做 `Stable=Opt`；比较不同多重度；用 `Guess=Only Pop=Full` 检查初猜占据。

**5. 开壳层计算 `<S**2>` 远大于 S(S+1)**

- 原因：自旋污染，混入了更高自旋态。偏差超过约 10% 时要警惕。
- 处理：做稳定性分析；考虑 ROKS、broken-symmetry 处理，或者体系本身需要多参考方法。见[第13章](13-openshell-stability.md)。

**6. 用了 `SCF=QC` 的 ROHF/ROKS 计算报错退出**

- 原因：官方文档说明 SCF=QC 不适用于 restricted open-shell。
- 处理：改用 UHF/UKS，或保留 RO 方法并在 route 中加 `Use=L506`。

**7. `Guess=Read` 读不到轨道，又从 Harris 初猜开始**

- 原因：`%chk` 文件名写错、路径不对，或 chk 来自不同版本或平台的 Gaussian。
- 处理：检查文件名；用 `%oldchk` 显式指定；跨机器先 `formchk` 再在新机器上 `unfchk`；C.01 之前生成的 chk 要用 `c8616` 升级（见[第16章](16-postprocess-multiwfn.md)）。

**8. 同一体系与文献能量对不上**

- 原因：Gaussian 版本（G09 默认 FineGrid，G16 默认 UltraFine）、SCRF 默认设置、物理常数版本或基组定义不同。
- 处理：比较相对能量而不是绝对能量；复现 G09 结果时可以用 `G09Defaults` 关键词。

更完整的报错索引见[第18章](18-errors-capstone.md)和 [Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)。

---

## 练习

1. 运行 `h2o_sp_b3lyp_def2tzvp.gjf`，把 `SCF Done` 能量换算成 kJ/mol 和 eV；然后把基组依次换成 def2SVP、def2QZVP，观察能量变化，体会基组收敛。
2. 对 `fe_complex_sp_scf_xqc.gjf` 去掉 `NoSymm` 和 `XQC` 各跑一次，比较迭代轮数和最终能量。再把多重度改成 1 和 3（电荷不变）计算，比较三个自旋态的能量顺序。
3. 运行 `guess_read_small_to_big_basis.gjf`，记录两步的 `after N cycles`；然后把第二步单独拿出来、不加 `Guess=Read` 重跑，比较轮数。
4. 分别运行 `h2_stretched_stable_opt.gjf` 和 `h2_stretched_guess_fragment.gjf`，比较两者的最终能量与 `<S**2>`，并与 R–H 距离为 0.74 Å 时的 RB3LYP 结果对照，思考为什么在拉长的键上 RKS 不稳定。

---

## 参考资料

**官方文档**

- [Gaussian 16 SCF keyword](https://gaussian.com/scf/)：XQC、QC、YQC、VShift、Fermi、NoIncFock、MaxCycle、Conver 的定义与默认值
- [Guess keyword](https://gaussian.com/guess/)：Harris、Huckel、Read、Mix、Fragment、Alter、Only
- [Stable keyword](https://gaussian.com/stable/)：波函数稳定性分析与 Stable=Opt
- [Integral keyword](https://gaussian.com/integral/)：UltraFine（G16 默认）、SuperFine 格点定义
- [Geom keyword](https://gaussian.com/geom/)：Check、AllCheck 读取 chk
- [Link 0 Commands](https://gaussian.com/link0/)：%chk、%oldchk、%mem、%nprocshared
- [Gaussian 16 Release Notes](https://gaussian.com/relnotes/)：G09 → G16 默认值变化、G09Defaults

**相关 GitHub 仓库（仅外链，未复制其代码或文本）**

- [liyuanhe211/Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)：按 Link 编号整理的 Gaussian 报错大全，L502 一节值得通读（CC BY-NC-SA 4.0）
- [cclib/cclib](https://github.com/cclib/cclib)：解析 Gaussian log 中的 SCF 能量、轨道能级等，可用来交叉验证自己的脚本
- [maserasgroup-repo/pyssian](https://github.com/maserasgroup-repo/pyssian)：按 Link 分块解析 Gaussian 输出的 Python 库
- [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash)：g16 bash 工具集，含 getenergy 和调用 testrt 的 route 语法检查
- [1234zou/MOKIT](https://github.com/1234zou/MOKIT)：在 Gaussian 与多参考程序之间传递轨道，适合 broken-symmetry 和强关联体系的后续处理

**关键文献**

- P. Pulay, Improved SCF convergence acceleration, *J. Comput. Chem.* **1982**, 3, 556. DOI: [10.1002/jcc.540030413](https://doi.org/10.1002/jcc.540030413)
- K. N. Kudin, G. E. Scuseria, E. Cancès, A black-box self-consistent field convergence algorithm: One step closer, *J. Chem. Phys.* **2002**, 116, 8255. DOI: [10.1063/1.1470195](https://doi.org/10.1063/1.1470195)
- G. B. Bacskay, A quadratically convergent Hartree-Fock (QC-SCF) method, *Chem. Phys.* **1981**, 61, 385. DOI: [10.1016/0301-0104(81)85156-7](https://doi.org/10.1016/0301-0104(81)85156-7)
- A. D. Rabuck, G. E. Scuseria, Improving self-consistent field convergence by varying occupation numbers, *J. Chem. Phys.* **1999**, 110, 695. DOI: [10.1063/1.478177](https://doi.org/10.1063/1.478177)
- V. R. Saunders, I. H. Hillier, A "Level-Shifting" method for converging closed shell Hartree-Fock wave functions, *Int. J. Quantum Chem.* **1973**, 7, 699. DOI: [10.1002/qua.560070407](https://doi.org/10.1002/qua.560070407)
- R. Bauernschmitt, R. Ahlrichs, Stability analysis for solutions of the closed shell Kohn-Sham equation, *J. Chem. Phys.* **1996**, 104, 9047. DOI: [10.1063/1.471637](https://doi.org/10.1063/1.471637)
- S. Lehtola, An overview of self-consistent field calculations within finite basis sets, *Molecules* **2020**, 25, 1218. DOI: [10.3390/molecules25051218](https://doi.org/10.3390/molecules25051218)

---

[← 上一章](03-method-basis.md) | [返回目录](../README.md) | [下一章 →](05-opt.md)
