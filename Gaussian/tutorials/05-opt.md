# 第5章 几何优化（Geometry Optimization with Gaussian 16）

本章讲解 Gaussian 16 几何优化（geometry optimization）的原理和实际操作。内容包括：读懂四项收敛判据（Max/RMS Force、Max/RMS Displacement），`Opt=(CalcFC, ReadFC, CalcAll, Tight, MaxCycles, MaxStep)` 等常用选项，坐标前加 `-1` 冻结原子，以及用 `Opt=ModRedundant` 约束键长和二面角。还会介绍如何诊断 L9999 "Number of steps exceeded" 优化不收敛、怎样用 `Opt=Restart` 或 `Geom=Check Guess=Read` 续算，以及构象和对称性带来的陷阱。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 写出 `Opt Freq` 联合作业，读懂每一步的收敛表，并找到优化后的结构和能量；
- 知道默认、`Tight`、`VeryTight` 判据的具体数值，以及什么时候需要收紧；
- 根据体系特点选择 `CalcFC`、`ReadFC`、`CalcAll`、`MaxStep` 等选项；
- 冻结原子、冻结键长、键角或二面角，做受限优化；
- 优化失败时，从能量和力的变化判断原因，并选择正确的续算方式；
- 理解优化只能找到“最近的”驻点，知道构象搜索和对称性的影响。

适用场景：几乎所有 Gaussian 研究都从优化开始，包括反应物、产物、中间体和构象的极小点，表面簇模型的受限优化，以及作为过渡态搜索（[第7章](07-ts-irc.md)）和柔性扫描（[第8章](08-scan.md)）的基础。

> [!NOTE]
> 全部示例输入在 [../examples/05-opt/README.md](../examples/05-opt/README.md)。单点能与 SCF 收敛见[第4章](04-sp-scf.md)，频率验证见[第6章](06-freq-thermo.md)。

---

## 2. 原理速览

### 2.1 势能面与驻点

在 Born-Oppenheimer 近似下，能量 $E(\mathbf{x})$ 是核坐标的函数，即势能面（PES）。优化的目标是找到**驻点**，也就是梯度为零的点：

$$
\mathbf{g} = \frac{\partial E}{\partial \mathbf{x}} = \mathbf{0}
$$

驻点的类型由 Hessian 矩阵 $\mathbf{H} = \partial^2 E / \partial \mathbf{x}\,\partial \mathbf{x}^{\mathsf T}$ 的本征值决定：全部为正（扣除平动和转动）是**极小点**，恰有一个负本征值是**过渡态**（一阶鞍点）。优化本身只保证梯度足够小，**是否为极小点必须用频率计算验证**（[第6章](06-freq-thermo.md)）。

### 2.2 准牛顿法与 Berny 算法

在当前点附近把能量展开到二阶，得到牛顿步：

$$
\Delta\mathbf{x} = -\mathbf{H}^{-1}\mathbf{g}
$$

精确 Hessian 很贵，所以 Gaussian 的默认 **Berny 算法**用近似 Hessian：

1. 第一步的 Hessian 由经验价力场估计（`NewEstmFC`，默认），或用 `CalcFC` 精确计算、用 `ReadFC` 从 chk 读入；
2. 之后每步用梯度变化更新 Hessian（极小化用 BFGS，过渡态用 Bofill）；
3. 极小化默认用 **GEDIIS** 外推，过渡态默认用 **RFO** 步；
4. 用信赖半径（trust radius，`MaxStep`）限制每一步的步长；
5. 默认在**冗余内坐标**（键长、键角、二面角，`Opt=Redundant`）中进行，通常比笛卡尔坐标收敛快得多。

### 2.3 收敛判据

每一步之后 Gaussian 检查四项（单位：力为 Hartree/Bohr 或 Hartree/rad，位移为 Bohr 或 rad）：

| 判据 | 默认 | `Opt=Tight` | `Opt=VeryTight` |
|---|---|---|---|
| Maximum Force | 0.000450 | 0.000015 | 0.000002 |
| RMS Force | 0.000300 | 0.000010 | 0.000001 |
| Maximum Displacement | 0.001800 | 0.000060 | 0.000006 |
| RMS Displacement | 0.001200 | 0.000040 | 0.000004 |

收敛规则：

- 四项全部为 `YES` 即判定收敛；
- 例外：如果**力已经比阈值小两个数量级**（小于阈值的 1/100），即使位移没达标，也判定收敛。这时输出 `Optimization completed on the basis of negligible forces.`，通常出现在势能面很平坦、步长估计偏大的体系上。
- `Opt=Loose`（官方文档：最大步长 0.01 au、RMS 力 0.0017 au）只适合大分子的粗优化，之后必须用默认判据再优化。

> [!TIP]
> 官方 Opt 文档建议：使用 `Tight`/`VeryTight` 的 DFT 计算应同时使用 `Int=UltraFine`。G16 默认已经是 UltraFine，所以只有显式改成更粗格点时才需要注意这一点。

---

## 3. 第一个优化：乙醛 Opt Freq

### 3.1 输入文件

[acetaldehyde_opt_freq.gjf](../examples/05-opt/acetaldehyde_opt_freq.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=acetaldehyde_opt_freq.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

Acetaldehyde optimization + frequencies, B3LYP-D3(BJ)/def2-SVP

0 1
C    -0.648568    -0.051139    -0.015523
C     0.843397     0.072740     0.022052
O     1.605373    -0.848266    -0.252785
H    -1.049168     0.669064    -0.732641
H    -0.933409    -1.061969    -0.317599
H    -1.049822     0.157790     0.978848
H     1.232196     1.061779     0.317648

```

| 关键词 | 说明 |
|---|---|
| `Opt` | 几何优化到极小点（默认 Berny + GEDIIS，冗余内坐标） |
| `Freq` | 在优化好的结构上自动接着做频率计算。`Opt Freq` 会被 Gaussian 自动拆成两步，第二步从 chk 读取结构和波函数 |
| `B3LYP/def2SVP EmpiricalDispersion=GD3BJ` | 优化和频率常用的“性价比”水平，见[第3章](03-method-basis.md) |
| `%chk` | 每一步优化后都会更新，是续算的依据 |

> [!IMPORTANT]
> 频率必须与优化**完全相同**的方法、基组、溶剂模型和格点。`Opt Freq` 写在同一行可以自动保证这一点。

### 3.2 运行

```bash
g16 < acetaldehyde_opt_freq.gjf > acetaldehyde_opt_freq.log
grep -c "Normal termination" acetaldehyde_opt_freq.log   # expect 2 (Opt step + Freq step)
```

`Opt Freq` 是两步作业，所以成功时有**两个** `Normal termination`。

### 3.3 读懂每一步的输出

每一步优化以 `GradGradGrad...` 分隔，并给出步数（示意输出，数值仅供参考）：

```text
 Berny optimization.
 Search for a local minimum.
 Step number   3 out of a maximum of   38
 ...
         Item               Value     Threshold  Converged?
 Maximum Force            0.002871     0.000450     NO
 RMS     Force            0.000958     0.000300     NO
 Maximum Displacement     0.010245     0.001800     NO
 RMS     Displacement     0.004012     0.001200     NO
 Predicted change in Energy=-5.123456D-06
```

`out of a maximum of 38`：默认最大步数取 20 与“2 倍冗余内坐标数”中的较大者，所以分子越大，允许的步数越多。

收敛时（示意输出）：

```text
         Item               Value     Threshold  Converged?
 Maximum Force            0.000031     0.000450     YES
 RMS     Force            0.000012     0.000300     YES
 Maximum Displacement     0.000420     0.001800     YES
 RMS     Displacement     0.000188     0.001200     YES
 Predicted change in Energy=-1.234567D-08
 Optimization completed.
    -- Stationary point found.
```

随后是 `Optimized Parameters` 表，列出最终的键长、键角和二面角，`-DE/DX` 一列是该坐标上残余的力。

常用 grep：

```bash
# convergence table of every step (4 lines per step); last one is the final step
grep -E "^ (Maximum|RMS) +(Force|Displacement)" acetaldehyde_opt_freq.log | tail -n 4
# energy of every step
grep "SCF Done" acetaldehyde_opt_freq.log
# confirm convergence
grep -E "Stationary point found|negligible forces" acetaldehyde_opt_freq.log
# number of optimization steps
grep -c "Step number" acetaldehyde_opt_freq.log
```

**优化后的能量**是 `Stationary point found` 之前最后一个 `SCF Done`（Freq 步会在同一结构上重算 SCF，能量相同）。**优化后的结构**是最后一个 `Standard orientation`（或 `Input orientation`）块；也可以用 GaussView 打开 log 查看，或用本仓库的 [log2xyz.py](../scripts/log2xyz.py) 提取（`python ../scripts/log2xyz.py -h`）。

### 3.4 Freq 步里的“再检查一次”

官方 Freq 文档指出：频率作业的开头和结尾都会运行 L103（优化模块），用**精确 Hessian** 重新判断一次收敛。偶尔会出现这种情况：优化用近似 Hessian 判定收敛了，但用精确 Hessian 算出的下一步仍然超过阈值。所以要看 Freq 步末尾那张收敛表：

```bash
# the LAST convergence table in an Opt Freq log comes from the Freq step
grep -E "^ (Maximum|RMS) +(Force|Displacement)" acetaldehyde_opt_freq.log | tail -n 4
```

如果这里出现 `NO`，说明结构其实没有完全收敛，处理办法见第 8 节和[第6章](06-freq-thermo.md)。最稳妥的做法是用 `Opt=CalcAll`（每步都算精确 Hessian，收敛后自动做频率分析），或者 `Geom=Check Guess=Read Opt=(ReadFC,Tight)` 再优化一次。

---

## 4. 常用 Opt 选项一览

| 选项 | 默认/取值 | 作用与适用场景 |
|---|---|---|
| `MaxCycles=N` | max(20, 2×冗余内坐标数) | 最大优化步数。只有在“缓慢但稳定下降”时增大才有意义 |
| `MaxStep=N` | 30（即 0.30 Bohr 或 rad） | 初始信赖半径为 0.01N。振荡时减小到 5–10 |
| `CalcFC` | — | 第一步精确计算 Hessian。适合起始结构较差、软模式多、氢键体系，以及过渡态 |
| `ReadFC` | — | 从 chk 读入力常数，最好来自低水平的 **Freq** 作业，其次是低水平优化结束时的近似 Hessian |
| `RCFC` | — | 读入 chk 中的**笛卡尔**力常数，适合频率计算发现需要降低对称性时 |
| `CalcAll` | — | 每一步都算精确 Hessian，很贵但最稳；收敛后自动做频率分析 |
| `RecalcFC=N` | — | 第 1 步算一次精确 Hessian，之后每 N 步重算一次，是 CalcFC 与 CalcAll 的折中。网上也常见 `Recalc=N` 的写法，建议写官方文档中的全称 |
| `Tight` / `VeryTight` | — | 收紧判据（第 2.3 节表格），用于软模式多的体系或高精度频率 |
| `Loose` | — | 放宽判据，只用于大分子的预优化 |
| `GEDIIS` | 极小化默认 | — |
| `RFO` | 过渡态默认 | 极小化时，如果 GEDIIS 在平坦势能面上来回振荡，可以改用 RFO（G03 的默认算法） |
| `Cartesian` | — | 在笛卡尔坐标中优化。内坐标出问题（线性角、L103 报错）时使用；纯笛卡尔优化不能冻结变量 |
| `Redundant` | 默认 | 冗余内坐标 |
| `ModRedundant` | — | 增加、删除、冻结或扫描内坐标（第 6.2 节） |
| `ReadOptimize` | — | 读入原子列表，指定哪些原子参与优化（第 6.1 节） |
| `NoFreeze` | — | 解冻全部原子，例如从带冻结信息的 chk 读结构时 |
| `Restart` | — | 从 chk 续算中断的优化（第 9.1 节） |
| `NGoUp=N` | 1 | **C.01 新增**：允许能量上升 N 次后才改为只做线搜索；-1 表示一旦能量上升就只做线搜索 |
| `NGoDown=N` | 3 | **C.01 新增**：在鞍点附近最多混合 N 个负本征值方向，用于离开鞍点；-1 关闭（等价于 `NoDownHill`） |
| `MaxEStep=N` | 600（0.6）；ONIOM QuadMac 为 100 | **C.01 新增**：离开鞍点时的步长，为 N/1000 Bohr 或 rad |
| `EigenTest` / `NoEigenTest` | 对 TS 开启 | 检查 Hessian 曲率，主要用于过渡态，见[第7章](07-ts-irc.md) |

> [!NOTE]
> 一些旧教程会写 `Opt=GDIIS`。这是 G03/G09 时代的选项，G16 官方 Opt 文档已不再列出，极小化的默认算法是 GEDIIS。如需改变算法，用 `Opt=RFO` 等文档中列出的选项。

选项可以组合，写在括号中用逗号分隔：`Opt=(CalcFC,MaxCycles=200,MaxStep=10)`。

---

## 5. 实战 1：起始结构较差的分子，用 CalcFC 和 ReadFC

### 5.1 CalcFC + MaxCycles

甘氨酸有分子内氢键、可旋转的 NH2 和 COOH 基团，势能面相对平坦，经验价力场估计的初始 Hessian 不太好。[opt_maxcycles_calcfc.gjf](../examples/05-opt/opt_maxcycles_calcfc.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=opt_maxcycles_calcfc.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt=(CalcFC,MaxCycles=200) Freq

Glycine: optimization with an exact initial Hessian (CalcFC) and more steps

0 1
N    -1.634151     0.328147     0.190861
C    -0.512051    -0.572745    -0.133908
C     0.822684     0.147488    -0.006423
O     0.996235     1.312429     0.317342
O     1.870662    -0.648107    -0.296614
H    -1.476858     0.704566     1.128418
H    -1.568371     1.147693    -0.416974
H    -0.625980    -0.929357    -1.161529
H    -0.524307    -1.421701     0.555496
H     2.652137    -0.068413    -0.176668

```

| 选项 | 为什么这样设 |
|---|---|
| `CalcFC` | 第一步精确 Hessian，相当于多花一次“小频率计算”的时间，换来更准确的步长和方向；对小分子很划算 |
| `MaxCycles=200` | 给足步数，防止因为平坦势能面在默认步数内没走完 |
| `Freq` | 确认极小点（0 个虚频）。起始结构来自 RDKit + MMFF，可能不是最低构象，见第 10 节 |

检查 `Derivative Info` 一列：用了 CalcFC 时，`Initial Parameters` 表中显示为 `calculate D2E/DX2 analytically`，默认则显示 `estimate D2E/DX2`。

```bash
grep -m 3 "D2E/DX2" opt_maxcycles_calcfc.log
```

### 5.2 ReadFC：读入低水平频率的 Hessian

对较大的分子，在高水平上做 CalcFC 很贵。更聪明的办法是：先在低水平做 Opt Freq，再把那一步的**精确 Hessian** 读给高水平优化。[opt_readfc_from_lowlevel.gjf](../examples/05-opt/opt_readfc_from_lowlevel.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=glycine_readfc.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

Glycine step 1: low-level opt + freq (provides an exact Hessian in the chk)

0 1
N    -1.634151     0.328147     0.190861
C    -0.512051    -0.572745    -0.133908
C     0.822684     0.147488    -0.006423
O     0.996235     1.312429     0.317342
O     1.870662    -0.648107    -0.296614
H    -1.476858     0.704566     1.128418
H    -1.568371     1.147693    -0.416974
H    -0.625980    -0.929357    -1.161529
H    -0.524307    -1.421701     0.555496
H     2.652137    -0.068413    -0.176668

--Link1--
%nprocshared=8
%mem=12GB
%chk=glycine_readfc.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt=ReadFC Freq Geom=Check Guess=Read

Glycine step 2: def2-TZVP optimization starting from the def2-SVP Hessian

0 1

```

`ReadFC` 读取的是当前 `%chk` 中保存的力常数。低水平 Freq 算出的精确 Hessian 远好于低水平优化结束时的近似 Hessian（官方文档原意）。用 `Initial Parameters` 表的 `Derivative Info` 列确认读入成功（应显示从 chk 读入的力常数，而不是 `estimate`）。

### 5.3 什么时候用 CalcAll

`Opt=CalcAll` 每一步都算精确 Hessian，比 CalcFC 贵得多，只在以下情况使用：

- CalcFC + 小 MaxStep 仍然振荡、不收敛；
- 势能面非常平坦或形状特殊，近似 Hessian 的更新总是失败；
- 过渡态优化非常困难（[第7章](07-ts-irc.md)）。

CalcAll 收敛后会自动做频率分析，所以不需要再写 Freq。折中方案是 `Opt=(CalcFC,RecalcFC=5)`。

---

## 6. 冻结原子与约束优化

### 6.1 坐标前加 -1：冻结原子

在原子符号和坐标之间加一个**冻结码**：`0` 表示参与优化（可省略），`-1` 表示冻结。这是表面簇模型最常用的做法：冻结底层或整个金属簇，只优化吸附分子。

[opt_frozen_atoms_surface_cluster.gjf](../examples/05-opt/opt_frozen_atoms_surface_cluster.gjf)：Cu(111) 两层 10 原子簇（7+3，最近邻距离按 fcc Cu 晶格常数 3.615 Å 取 2.556 Å）顶位吸附 CO，冻结全部 Cu，只优化 CO：

```text
%nprocshared=8
%mem=12GB
%chk=opt_frozen_atoms_surface_cluster.chk
#p PBE1PBE/def2SVP EmpiricalDispersion=GD3BJ Opt NoSymm

CO on a frozen Cu10(111) two-layer cluster (Cu frozen with flag -1)

0 1
Cu  -1    0.000000    0.000000    0.000000
Cu  -1    2.556191    0.000000    0.000000
Cu  -1    1.278096    2.213726    0.000000
Cu  -1   -1.278096    2.213726    0.000000
Cu  -1   -2.556191    0.000000    0.000000
Cu  -1   -1.278096   -2.213726    0.000000
Cu  -1    1.278096   -2.213726    0.000000
Cu  -1    1.278096    0.737909   -2.087121
Cu  -1   -1.278096    0.737909   -2.087121
Cu  -1    0.000000   -1.475818   -2.087121
C    0    0.000000    0.000000    1.850000
O    0    0.000000    0.000000    3.000000

```

| 要点 | 说明 |
|---|---|
| `Cu  -1  x y z` | 冻结码 -1：该原子在优化中固定不动 |
| `C  0  x y z` | 冻结码 0：参与优化 |
| `NoSymm` | 官方文档提示冻结原子的优化可能需要 NoSymm；同时保证坐标保持输入取向，便于和表面模型对应 |
| 电子数 | 10×29 + 6 + 8 = 304（偶数），单重态 |
| 基组 | Cu（Z = 29）用 def2-SVP 是全电子基组；也可以 Cu 用 LANL2DZ + ECP、C/O 用 6-31G(d)，写成 GenECP（[第3章](03-method-basis.md)） |

> [!NOTE]
> G16 的 Freq 会把冻结原子排除在振动分析之外（输出中它们的位移为零，默认不打印，`Freq=PrintFrozenAtoms` 可以打印）。这种部分 Hessian 的频率可以用来确认吸附分子的局部极小，但由此得到的热化学量只是近似值。

**另一种写法：`Opt=ReadOptimize`。** 不改坐标，在分子坐标后的单独一节（以空行结束）中列出要优化的原子：

```text
noatoms atoms=11-12
```

意思是“先清空列表，再加入第 11、12 号原子”，其余原子冻结。按元素选择也可以，例如 `atoms=C O`。

### 6.2 Opt=ModRedundant：冻结键长、键角、二面角

在分子坐标后（空行分隔）加一节 ModRedundant 输入，语法是：

```text
[Type] N1 [N2 [N3 [N4]]] [Action [params]]
```

| 动作代码 | 含义 |
|---|---|
| （不写） | 添加该坐标（例如让某个氢键距离出现在每一步的参数表中） |
| `F` | 冻结该坐标 |
| `A` | 激活（解冻）该坐标 |
| `B` | 添加该坐标并生成所有相关坐标 |
| `K` | 删除该坐标及所有包含它的相关坐标 |
| `R` | 只删除该坐标本身 |
| `D` | 对该坐标数值计算初始 Hessian 的行与列 |
| `H val` | 把该坐标初始 Hessian 的对角元设为 val（Hartree/au²） |
| `S n step` | 柔性扫描：坐标增加 step，共 n 步，每点做受限优化（[第8章](08-scan.md)） |

坐标类型由原子个数决定：1 个原子为笛卡尔坐标（X），2 个为键长（B），3 个为键角（A），4 个为二面角（D）。也可以显式写 `B 1 4 F`、`D 1 2 3 4 F`；`L` 表示线性弯曲。`*` 是通配符，例如 `* F` 冻结所有笛卡尔坐标，`3 *` 表示所有与 3 号原子相连的键。

**示例 A：冻结 O···O 距离的水二聚体。** [opt_modredundant_freeze_bond.gjf](../examples/05-opt/opt_modredundant_freeze_bond.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=opt_modredundant_freeze_bond.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt=ModRedundant

Water dimer: O1...O4 distance frozen at its input value (3.10 A)

0 1
O     0.000000     0.000000     0.000000
H     0.970000     0.000000     0.000000
H    -0.240400     0.929400     0.000000
O     3.100000     0.000000     0.000000
H     3.420100    -0.492900     0.759100
H     3.420100    -0.492900    -0.759100

B 1 4 F

```

`B 1 4 F`：1 号和 4 号原子（两个 O）之间的距离在输入值 3.10 Å 处冻结，其余坐标全部优化。原子 1–4 之间本来没有化学键，ModRedundant 会自动**添加**这个坐标再冻结。这是构造氢键距离势能曲线、做受限采样的基本操作。

**示例 B：冻结丁烷的 C–C–C–C 二面角。** [opt_modredundant_freeze_dihedral.gjf](../examples/05-opt/opt_modredundant_freeze_dihedral.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=opt_modredundant_freeze_dihedral.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt=ModRedundant

n-Butane: C1-C2-C3-C4 dihedral frozen at its input value (65 deg)

0 1
C     1.864934    -0.513946     0.360610
C     0.436665    -0.062823     0.621340
C    -0.436665    -0.223603    -0.621339
C    -0.593960    -1.685979    -1.006253
H     2.322519     0.075508    -0.440330
H     2.472162    -0.390283     1.262776
H     1.896475    -1.569410     0.071905
H     0.443168     0.986397     0.938951
H     0.019523    -0.649258     1.448303
H     0.006390     0.323818    -1.461466
H    -1.427354     0.208907    -0.438992
H    -1.064595    -2.255793    -0.198687
H    -1.222319    -1.777825    -1.897685
H     0.377274    -2.140234    -1.226947

D 1 2 3 4 F

```

```bash
# the frozen coordinate is flagged in the parameter tables
grep -A 40 "Optimized Parameters" opt_modredundant_freeze_dihedral.log | grep "D(1,2,3,4)"
```

最终表中该二面角的值应保持 65°，并且 `-DE/DX` 不为零（约束坐标上存在残余力，属正常现象）。

> [!WARNING]
> `F` 冻结的是**输入结构中的当前值**。想把坐标冻结在某个特定数值（例如二面角 90°），要先在 GaussView 中把结构调到目标值再导出。另外，ModRedundant 中的原子编号就是输入中原子的顺序（不计虚原子）；在 GaussView 中重新排序原子后，约束编号必须同步修改。

> [!IMPORTANT]
> 带约束的优化得到的不是驻点（约束方向上仍有力），所以**不要**在受限优化的结构上做频率和热化学分析，除非你明确知道自己在做什么（例如受限的模型簇）。

---

## 7. 用 GaussView 和脚本监控优化过程

优化失败时，第一件事是**看能量和力随步数的变化**，而不是马上加 MaxCycles。

**GaussView：** 打开 log 文件，`Results → Optimization...` 会画出每一步的能量和 RMS 梯度，可以逐帧查看结构。

**命令行：**

```bash
# step energies (one per optimization step; the Freq step adds one more at the end)
grep "SCF Done" job.log | awk '{print NR, $5}'
# Max Force of every step
grep "^ Maximum Force" job.log | awk '{print NR, $3, $5}'
```

**本仓库脚本：** [opt_monitor.py](../scripts/opt_monitor.py) 会解析每一步的能量和四项判据，并画出它们与阈值的对比图（需要 matplotlib；`python ../scripts/opt_monitor.py -h` 查看用法）。

---

## 8. 优化不收敛的诊断（L9999）

超过最大步数时，日志末尾是（示意输出）：

```text
 Optimization stopped.
    -- Number of steps exceeded,  NStep= 100
    -- Flag reset to prevent archiving.
 ...
 Error termination request processed by link 9999.
 Error termination via Lnk1e in /opt/g16/l9999.exe at ...
```

按能量和力的曲线形状判断原因：

| 曲线形状 | 判断 | 处理 |
|---|---|---|
| 能量单调下降，力稳定减小，只是没走完 | 起始结构离极小点远，或分子很大、很软 | 从最后一步续算（第 9 节），增大 `MaxCycles` |
| 能量在 1e-5 ~ 1e-6 Eh 范围内上下振荡，力不再减小 | 近似 Hessian 差，或数值噪声 | `Opt=CalcFC`（或 `ReadFC`），`MaxStep=5~10`；软模式多时用 `Int=SuperFine`；可试 `Opt=RFO` |
| 只有位移判据一直是 NO，力已经很小 | 势能面非常平坦（甲基转动、弱复合物） | 通常结构已经可用；需要高质量频率时用 `Opt=(Tight,CalcFC)` + `Int=SuperFine` |
| 能量突然跳变（升高或大幅下降） | SCF 收敛到了不同的电子态，或结构发生了构象跳变 | 检查每一步的 `<S**2>` 和结构；用 `Guess=Read`、`Stable` 检查（[第4章](04-sp-scf.md)） |
| 一开始就很奇怪：键断开、原子乱飞 | 起始结构不合理、电荷或多重度错误 | 回到建模阶段，先用 xtb 或小基组预优化 |
| 接近线性角时 L103 报错 | 冗余内坐标在 180° 附近奇异 | `Opt=Cartesian`，或从最新结构重启（第 9 节、常见问题） |

> [!TIP]
> 一个简单的经验：**振荡时加步数没用，缓慢下降时加步数才有用。** 振荡的优化给 500 步也只是在原地打转。

---

## 9. 从中断或失败的作业续算

### 9.1 Opt=Restart：完全接着原来的优化走

[opt_restart_from_chk.gjf](../examples/05-opt/opt_restart_from_chk.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=opt_maxcycles_calcfc.chk
#p Opt=(Restart,MaxCycles=300)

```

| 要点 | 说明 |
|---|---|
| 整个输入只有 Link0 和 route | 官方文档：route 中只写 Opt 和原作业的选项，再加 Restart；方法、基组、结构和 Hessian 都从 chk 读取，**不要**写 title 和坐标 |
| `%chk` | 必须是原作业的 chk（原作业必须写了 `%chk`），而且是同一版本 Gaussian 生成的 |
| `MaxCycles=300` | 这类参数选项可以在续算时修改 |
| `CalcFC` | 想在续算的第一步重新计算精确 Hessian，就必须在续算的 route 中写上 CalcFC，不论原作业是否写过 |
| Freq | 续算时不会自动接着做 Freq；优化收敛后另外用 `Freq Geom=AllCheck Guess=Read` 计算 |

> [!CAUTION]
> 续算会继续写入同一个 chk。先备份：`cp opt_maxcycles_calcfc.chk opt_maxcycles_calcfc.chk.bak`。如果原作业是因为结构跑偏而失败，Restart 只会沿着错误的方向继续走，这时应该用 9.2 节的方法从合适的结构重新开始。

### 9.2 Geom=Check Guess=Read：从最后一步结构重新开始

[opt_continue_geom_check.gjf](../examples/05-opt/opt_continue_geom_check.gjf)：

```text
%nprocshared=8
%mem=12GB
%oldchk=opt_maxcycles_calcfc.chk
%chk=glycine_opt_continue.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt=(CalcFC,MaxCycles=200) Freq Geom=Check Guess=Read

Glycine: continue the optimization from the last geometry stored in the old chk

0 1

```

和 Opt=Restart 的区别：

| | `Opt=Restart` | `Geom=Check Guess=Read` |
|---|---|---|
| 结构 | chk 中最后一步 | chk 中最后一步（也可以用 `Geom=(Check,Step=N)` 取第 N 步） |
| Hessian | 继承原作业的近似 Hessian | **丢失**，需要 CalcFC 或 ReadFC |
| 方法和选项 | 不能改 | 可以改（例如加 `Tight`、换基组、加 Freq） |
| 原 chk | 继续写入 | 用 `%oldchk` 只读，新结果写入新 chk |
| 适用 | 作业被意外中断（超时、节点故障） | 想换选项、换水平，或要从某个中间步重新开始 |

`Geom=(Check,Step=N)`：官方文档说明它可以取失败或部分完成的优化的第 N 步结构（不适用于成功完成的优化）；`Step=Original` 取初始结构。不是每一步都一定保存在 chk 中，日志中有 `Hessian updated` 信息的步骤才可以取。

### 9.3 没有 chk 时

从 log 中取最后一个结构（或能量最低的结构）重新写输入：

- GaussView 打开 log，选择目标帧，`File → Save` 另存为 gjf；
- 用 [log2xyz.py](../scripts/log2xyz.py) 提取最后一步或能量最低一步，它也能直接写出续算用的 gjf（`python ../scripts/log2xyz.py -h`）；
- 振荡体系建议取**能量最低**的那一步，而不是最后一步。

---

## 10. 低能构象、预优化与对称性

### 10.1 优化只找最近的极小点

Berny 优化是局部算法，只会滑到起始结构附近的那个极小点。柔性分子（多个可旋转键）有很多构象，**你优化的那个不一定是最低的**。标准流程：

1. 构象搜索：[CREST](https://github.com/crest-lab/crest)（GFN2-xTB 元动力学，推荐）、RDKit ETKDG、molclus 等；
2. 用 xtb 或小基组 DFT 对候选构象粗优化、去重；
3. 能量窗口内（例如 3–5 kcal/mol）的构象做 DFT `Opt Freq`；
4. 按自由能做 Boltzmann 加权（[第6章](06-freq-thermo.md)、[第18章](18-errors-capstone.md)）。

### 10.2 预优化

对较大或起始结构很差的分子，先用便宜的方法预优化能省很多 DFT 步数：

```bash
# GFN2-xTB pre-optimization (xtb writes xtbopt.xyz)
xtb start.xyz --opt
```

然后把 `xtbopt.xyz` 的坐标放进 Gaussian 输入。也可以用 [xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian) 通过 Gaussian 的 `External` 接口直接调用 xtb，用 Gaussian 自己的优化器在 GFN2-xTB 水平上预优化。另一种常见做法是先用 def2-SVP 优化，再用 `ReadFC` 在大基组上精修（第 5.2 节）。

### 10.3 对称性：保持还是放弃

Gaussian 默认识别分子点群并在优化中**保持对称性**：对称结构上的力也是对称的，所以优化不会自发破缺对称性。这带来两种结果：

- 好处：高对称分子计算更快，结构更“干净”；
- 风险：如果起始结构的对称性高于真实极小点，优化会收敛到**高对称的鞍点**。例如从完全平面的结构开始优化 NH3，会停在平面过渡态上。

应对方法：

| 情形 | 做法 |
|---|---|
| 希望得到对称的极小点 | 优化前在 GaussView 中用 `Edit → Symmetrize` 把结构对称化，让 Gaussian 识别出正确的点群 |
| 怀疑真实结构对称性更低 | 频率计算中出现虚频，沿虚频方向位移（GaussView 振动动画中手动位移，或用 pyQRC），破坏对称性后重新优化；必要时加 `NoSymm` |
| 冻结原子、ModRedundant 约束导致对称性不成立 | 加 `NoSymm` |
| 需要保持输入取向（后续做 cube 差分、表面模型） | 加 `NoSymm` |

`NoSymm` 的代价是失去对称性加速，并且数值噪声可能让本来对称的结构有微小的不对称。没有特殊需要时，保持默认即可。

---

## 参数选择与收敛建议

| 情形 | 推荐设置 | 说明 |
|---|---|---|
| 普通有机小分子 | `Opt Freq` | 默认足够 |
| 起始结构较差、软模式多 | `Opt=(CalcFC,MaxCycles=200)` | 小分子 CalcFC 很划算 |
| 较大分子的高水平优化 | 低水平 `Opt Freq` → 高水平 `Opt=ReadFC` | 第 5.2 节 |
| 优化振荡 | `Opt=(CalcFC,MaxStep=10)`；必要时 `Int=SuperFine` 或 `Opt=RFO` | 加 MaxCycles 无效 |
| 缓慢下降 | 续算并增大 `MaxCycles` | `Opt=Restart` 或 `Geom=Check Guess=Read` |
| 弱相互作用复合物、低频模式要准 | `Opt=Tight` + `Int=SuperFine`（频率用同样的格点） | 第6章 |
| 非常困难的情形 | `Opt=CalcAll` 或 `Opt=(CalcFC,RecalcFC=5)` | 贵但稳 |
| 线性或接近线性的结构 | `Opt=Cartesian` | 避免 L103 线性角问题 |
| 表面簇模型 | 冻结码 `-1`，或 `Opt=ReadOptimize`；加 `NoSymm` | 第 6.1 节 |
| 固定某个几何参数 | `Opt=ModRedundant` + `B/A/D ... F` | 第 6.2 节 |
| 大分子粗优化 | `Opt=Loose`（之后必须用默认判据再优化） | 不可单独使用 |
| 柔性分子 | CREST 构象搜索 → xtb 粗筛 → DFT `Opt Freq` | 第 10 节 |

---

## 常见问题与报错

**1. `Optimization stopped. -- Number of steps exceeded, NStep= 100`（L9999）**

- 原因：在最大步数内未收敛。
- 处理：先按第 8 节判断是振荡还是缓慢下降。缓慢下降就续算并增大 MaxCycles；振荡就用 CalcFC 和较小的 MaxStep，并检查连接性和起始结构。

**2. `FormBX had a problem` / `Linear angle in Bend` / `Linear angle in Tors`（L103）**

- 原因：键角接近 180°（炔烃、M–C≡O、线性三原子片段），冗余内坐标出现奇异。
- 处理：从最新结构重启（`Geom=Check Guess=Read`）让 Gaussian 重新生成内坐标；或用 `Opt=Cartesian`；或加虚原子；或在 ModRedundant 中用 `L` 类型定义线性弯曲。

**3. `Error in internal coordinate system.`（L103）**

- 原因：优化中结构变化很大，原来的内坐标定义不再适用。
- 处理：取最后一个合理结构，用 `Geom=Check Guess=Read` 重新开始（内坐标会重建），或改用 `Opt=Cartesian`。

**4. `Small interatomic distances encountered` / `Problem with the distance matrix`（L202）**

- 原因：原子重叠，常见于复制粘贴时单位错误（Bohr 与 Å）、重复原子，或者优化中结构崩坏。
- 处理：在 GaussView 中检查结构；如果是优化中途出现，换起始结构或减小 MaxStep。

**5. 优化收敛了，但 Freq 步的收敛表出现 `NO`**

- 原因：近似 Hessian 判定收敛，但用精确 Hessian 判断还没收敛（第 3.4 节）。
- 处理：`Geom=Check Guess=Read Opt=(ReadFC,Tight) Freq` 再优化一次，或直接用 `Opt=CalcAll`。

**6. 优化收敛，但出现一个小虚频（几十 i cm⁻¹ 以内）**

- 原因：收敛不够紧、格点噪声，或对称性把结构限制在鞍点上。
- 处理：`Opt=(Tight,CalcFC)` + `Int=SuperFine`；沿虚频方向位移后重新优化。详见[第6章](06-freq-thermo.md)。

**7. 冻结了原子却报错，或冻结没有生效**

- 原因：纯笛卡尔优化（`Opt=Cartesian`）不支持冻结变量；或者从 chk 读入的结构带着旧的冻结信息。
- 处理：冻结原子时使用默认的冗余内坐标；需要全部解冻时加 `Opt=NoFreeze`。

**8. 用了 ModRedundant，约束的却是错误的键或角**

- 原因：原子编号和输入顺序不一致（GaussView 中重新排序过，或者编号时把虚原子算了进去）。
- 处理：在 GaussView 中打开输入，用标签显示原子编号后再核对。

更多报错见[第18章](18-errors-capstone.md)和 [Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)。

---

## 练习

1. 运行 `acetaldehyde_opt_freq.gjf`，记录优化步数和最终能量；把 `Opt` 改成 `Opt=CalcFC` 和 `Opt=Tight` 各跑一次，比较步数、耗时和最终能量的差别（应在 1e-6 Eh 量级）。
2. 把 `opt_modredundant_freeze_bond.gjf` 中 O···O 的初始距离分别改成 2.7、2.9、3.1、3.3、3.5 Å（移动第二个水分子），得到受限优化的能量，画出势能曲线，并与[第8章](08-scan.md)的柔性扫描结果比较。
3. 故意把 `opt_maxcycles_calcfc.gjf` 改成 `Opt=(MaxCycles=5)` 运行，让它因 L9999 失败，然后分别用 `opt_restart_from_chk.gjf` 和 `opt_continue_geom_check.gjf` 续算，比较总步数。
4. 用 RDKit 或 CREST 生成甘氨酸的多个构象，分别优化，比较能量，看看示例中的起始结构是不是最低构象。

---

## 参考资料

**官方文档**

- [Gaussian 16 Opt keyword](https://gaussian.com/opt/)：Berny 算法说明、全部选项（CalcFC、ReadFC、CalcAll、RecalcFC、Tight、MaxCycles、MaxStep、ModRedundant、ReadOptimize、Restart、NGoUp/NGoDown/MaxEStep）
- [Geom keyword](https://gaussian.com/geom/)：Check、AllCheck、Step=N、NGeom、ModRedundant
- [Freq keyword](https://gaussian.com/freq/)：Freq 作业中 L103 再次检查收敛的说明
- [Integral keyword](https://gaussian.com/integral/)：UltraFine 与 SuperFine 格点
- [Gaussian 16 Release Notes](https://gaussian.com/relnotes/)：C.01 新增的 Opt=NGoUp/NGoDown/MaxEStep
- [CREST documentation](https://crest-lab.github.io/crest-docs/)：构象搜索

**相关 GitHub 仓库（仅外链，未复制其代码或文本）**

- [crest-lab/crest](https://github.com/crest-lab/crest)：GFN2-xTB 构象搜索，DFT 优化前的标准步骤
- [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian)：通过 External 让 Gaussian 调用 xtb，在 GFN2-xTB 水平上预优化
- [jvalegre/aqme](https://github.com/jvalegre/aqme)：从 SMILES 做构象搜索、生成 Gaussian 输入，并自动检查虚频和报错后修复重算
- [chunxiangzheng/gaussian_log_file_converter](https://github.com/chunxiangzheng/gaussian_log_file_converter)：把优化 log 中的结构转成 xyz/gjf
- [patonlab/pyQRC](https://github.com/patonlab/pyQRC)：沿虚频方向位移生成新输入，用于消除残留虚频
- [thisisntnathan/dftCourse](https://github.com/thisisntnathan/dftCourse)：Cornell 的 Gaussian 课程，含几何优化和效率章节

**关键文献**

- H. B. Schlegel, Optimization of equilibrium geometries and transition structures, *J. Comput. Chem.* **1982**, 3, 214. DOI: [10.1002/jcc.540030212](https://doi.org/10.1002/jcc.540030212)
- C. Peng, P. Y. Ayala, H. B. Schlegel, M. J. Frisch, Using redundant internal coordinates to optimize equilibrium geometries and transition states, *J. Comput. Chem.* **1996**, 17, 49. DOI: [10.1002/(SICI)1096-987X(19960115)17:1&lt;49::AID-JCC5&gt;3.0.CO;2-0](https://doi.org/10.1002/%28SICI%291096-987X%2819960115%2917:1%3C49::AID-JCC5%3E3.0.CO;2-0)
- X. Li, M. J. Frisch, Energy-represented direct inversion in the iterative subspace within a hybrid geometry optimization method, *J. Chem. Theory Comput.* **2006**, 2, 835. DOI: [10.1021/ct050275a](https://doi.org/10.1021/ct050275a)
- H. B. Schlegel, Geometry optimization, *WIREs Comput. Mol. Sci.* **2011**, 1, 790. DOI: [10.1002/wcms.34](https://doi.org/10.1002/wcms.34)
- P. Pracht, F. Bohle, S. Grimme, Automated exploration of the low-energy chemical space with fast quantum chemical methods, *Phys. Chem. Chem. Phys.* **2020**, 22, 7169. DOI: [10.1039/C9CP06869D](https://doi.org/10.1039/C9CP06869D)
- C. Bannwarth, S. Ehlert, S. Grimme, GFN2-xTB, *J. Chem. Theory Comput.* **2019**, 15, 1652. DOI: [10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176)

---

[← 上一章](04-sp-scf.md) | [返回目录](../README.md) | [下一章 →](06-freq-thermo.md)
