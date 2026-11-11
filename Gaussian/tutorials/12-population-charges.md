# 第12章 布居分析、NBO 与静电势/RESP 电荷（Population Analysis, NBO & ESP/RESP Charges with Gaussian）

本章系统讲解 Gaussian 16 中的原子电荷与成键分析：Mulliken 布居（`Pop=Regular/Full`）及其基组依赖性、内置 NBO 3.1 与外部 NBO 7（NPA 电荷、Wiberg 键级、二阶微扰 E(2) 供体-受体分析）、静电势拟合电荷（`Pop=MK`、`Pop=CHelpG`）、AMBER 力场所需的 RESP 电荷（antechamber 或 Multiwfn 拟合），以及 Hirshfeld/CM5 电荷和 Multiwfn 中的 ADCH、Mayer 键级、轨道成分分析。每种方法都给出完整输入、输出定位命令和“该用哪种电荷”的选择建议。

## 1. 本章目标与适用场景

学完本章你应该能够：

- 说清楚 Mulliken、NPA、Hirshfeld/CM5、ADCH、MK/CHelpG、RESP 各自的定义和适用范围，并且不再用 Mulliken 电荷讨论定量问题；
- 运行 NBO 分析，读懂自然布居、Lewis 结构、杂化轨道和 E(2) 表，用 `$NBO BNDIDX $END` 得到 Wiberg 键级；
- 按 AMBER/GAFF 惯例计算 HF/6-31G(d) 静电势并拟合两步 RESP 电荷；
- 用 `Density=Check` 在同一个密度上快速比较多种电荷，避免重复做 SCF；
- 为 ONIOM/MM 保存 CM5 或 Hirshfeld 电荷（`Pop=SaveCM5`/`SaveHirshfeld`，C.01 新增）。

| 任务 | 本章小节 | 配套示例 |
|---|---|---|
| Mulliken 与基组依赖性 | 第 4 节 | [12-nbo-pyridine-nma](../examples/12-nbo-pyridine-nma/README.md) |
| NBO、Wiberg 键级、E(2) | 第 5、6 节 | 同上 |
| ESP 拟合与 RESP | 第 7 节 | [12-esp-resp-paracetamol](../examples/12-esp-resp-paracetamol/README.md) |
| Hirshfeld/CM5、多种电荷对比 | 第 8 节 | [12-hirshfeld-cm5-acetate](../examples/12-hirshfeld-cm5-acetate/README.md) |

前置知识：[第4章](04-sp-scf.md) 单点计算、[第5章](05-opt.md) 几何优化。波函数后处理的通用流程（formchk、cubegen、Multiwfn）见 [第16章](16-postprocess-multiwfn.md)。

## 2. 原理速览

### 2.1 原子电荷不是可观测量

“原子电荷”取决于你怎样把分子的电子密度分给各个原子，不同定义之间没有对错之分，只有是否适合某个用途。常见方法分三类：

| 类别 | 方法 | 核心思想 | 主要特点 |
|---|---|---|---|
| 基于轨道/基函数 | Mulliken、Löwdin | 按基函数所属原子划分密度矩阵 | 便宜，但强烈依赖基组，弥散函数下可能失去意义 |
| 基于轨道/基函数 | NPA（NBO） | 先构造自然原子轨道（NAO）再计数 | 基组依赖性小，化学直观，适合讨论电荷转移 |
| 基于实空间密度划分 | Hirshfeld、CM5、ADCH、AIM | 用参考原子密度或零通量面划分 ρ(r) | Hirshfeld 绝对值偏小；CM5、ADCH 修正后能较好重现偶极矩 |
| 静电势拟合 | MK、CHelpG、RESP | 让点电荷尽量重现分子外部的静电势 | 适合分子力学和静电相互作用；埋藏原子的电荷不稳定 |

几个关键公式：

- Mulliken：$q_A=Z_A-\sum_{\mu\in A}(\mathbf{PS})_{\mu\mu}$，P 为密度矩阵，S 为重叠矩阵；
- Hirshfeld：权重函数 $w_A(\mathbf r)=\rho_A^0(\mathbf r)/\sum_B\rho_B^0(\mathbf r)$，$q_A=Z_A-\int w_A(\mathbf r)\rho(\mathbf r)\,d\mathbf r$；CM5 在 Hirshfeld 电荷上加经验修正，使其更好地重现偶极矩（Marenich 等 2012）；
- ESP 拟合：在分子范德华表面外的格点 $\mathbf r_k$ 上最小化

$$\chi^2=\sum_k\Big[V_{\rm QM}(\mathbf r_k)-\sum_A\frac{q_A}{|\mathbf r_k-\mathbf R_A|}\Big]^2,\qquad \sum_A q_A=Q_{\rm tot}$$

- RESP：在 χ² 上增加双曲型约束 $a\sum_A\big(\sqrt{q_A^2+b^2}-b\big)$，把埋藏原子（如 sp³ 碳）的电荷拉向 0，并分两步拟合（第二步只重新拟合 CH3/CH2 等需要等价约束的基团）。

### 2.2 NBO 的思路

NBO 把一阶约化密度矩阵依次变换为自然原子轨道（NAO）→ 自然杂化轨道（NHO）→ 自然键轨道（NBO），得到最接近 Lewis 结构的定域描述：成键轨道 BD、孤对 LP、芯轨道 CR，以及对应的反键 BD*、Rydberg RY*。Lewis 型轨道向非 Lewis 型轨道的离域用二阶微扰能估计：

$$E^{(2)}_{i\to j^*}=-q_i\,\frac{F_{ij}^2}{\varepsilon_{j^*}-\varepsilon_i}$$

$q_i$ 为供体轨道占据数，$F_{ij}$ 为 NBO 基下的 Fock 矩阵元。E(2) 越大，供体-受体作用（共轭、超共轭、氢键中的 n→σ*）越强。

## 3. Pop 关键词速查

`Population`（简写 `Pop`）是性质关键词，**必须带选项**。官方 [Population 页](https://gaussian.com/population/) 的要点：单点作业只做一次布居分析，几何优化在第一步和最后一步各做一次；分析使用的密度由 `Density` 关键词控制；一个作业步只能用一种密度和一种电荷拟合方法；结果在标准取向下给出。

| 选项 | 说明 |
|---|---|
| `Minimal` | 默认：只打印原子总电荷和轨道能 |
| `Regular` | 5 个最高占据和 5 个最低空轨道的系数、密度矩阵、完整 Mulliken 分析 |
| `Full` | 同 Regular，但打印全部轨道 |
| `Orbitals=N` / `AllOrbitals` | 每个 MO 的原子贡献分析（默认 N=10） |
| `Hirshfeld`（同义 `CM5`） | Hirshfeld 电荷，同时给出 CM5 电荷 |
| `MK`（同义 `ESP`、`MerzKollman`） | Merz-Singh-Kollman 静电势拟合 |
| `CHelpG` / `CHelp` / `HLY` | 其他 ESP 拟合格点方案 |
| `MKUFF` | MK 拟合但使用 UFF 半径（覆盖整个周期表） |
| `ReadRadii` / `ReadAtRadii` | 按元素 / 按原子读入拟合用半径（Å），空行结束 |
| `Dipole` | 拟合时约束电荷重现偶极矩 |
| `Resp=N` | C.01 起：在 ESP 拟合中加入 RESP 约束，权重 N×10⁻⁶，默认 N=2（见 [C.01 发行说明](https://gaussian.com/relnotes/)） |
| `NBO` / `NPA` / `NBORead` / `NBODel` | 内置 NBO 3.1：完整 NBO / 仅 NPA / 从输入读 $NBO 关键词 / 删除分析 |
| `NBO7` / `NPA7` / `NBO7Read` / `NBO7Delete` | C.01 起：调用外部 NBO 7（需另行获得并安装） |
| `SaveNBOs` / `SaveNLMOs` / `SaveMixed` | 把 NBO/NLMO 写入 chk 供可视化 |
| `SaveMulliken` / `SaveESP` / `SaveNPA` / `SaveCM5` / `SaveHirshfeld` | 把电荷存入 chk，供后续 MM/ONIOM 计算使用 |
| `NCS` | NBO 自然化学屏蔽分析（配合 GIAO NMR） |
| `NTO` / `SaveNTO` / `DCT` | 激发态分析，见 [第10章](10-tddft.md) |

> [!NOTE]
> 截至 2026-11，官方 Population 页面上 RESP 选项的说明处于隐藏状态，但 C.01 发行说明明确列出了 `Pop=(MK,Resp=N)`。它是在单步 ESP 拟合中加入对电荷平方的约束，与 AMBER 的标准两步双曲型 RESP 并不完全相同；为 GAFF/AMBER 力场拟合电荷时，仍建议使用第 7 节的 antechamber 或 Multiwfn 流程。

## 4. 实战一：Mulliken 电荷及其基组依赖性

示例输入 `pyridine_mulliken_basis.gjf` 在同一个吡啶结构上，用 5 个基组（STO-3G、6-31G(d)、6-311++G(d,p)、def2-TZVP、aug-cc-pVTZ）各算一次单点，并同时输出 Mulliken 与 Hirshfeld/CM5 电荷。其中一步如下（其余步骤只有基组不同）：

```text
%oldchk=pyridine_opt.chk
%chk=pyridine_631gd.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-31G(d) Pop=Hirshfeld Geom=Check

Pyridine Mulliken vs Hirshfeld/CM5 charges, B3LYP/6-31G(d)

0 1

```

`pyridine_opt.chk` 由同目录的 `pyridine_pop_nbo.gjf` 第 1 步生成，所以要先运行它。Mulliken 电荷总会打印（`Pop=Minimal` 也包含），无需单独请求：

```bash
grep -A12 "Mulliken charges:" pyridine_mulliken_basis.log | grep -E "^ +4 +N"
grep -A12 "Hirshfeld charges, spin densities" pyridine_mulliken_basis.log | grep -E "^ +4 +N"
```

吡啶中 N 为 4 号原子。比较不同基组下 N 的 Mulliken 电荷和 CM5 电荷，你会看到 Mulliken 值随基组变化明显，加入弥散函数后甚至可能出现不合理的数值，而 CM5 相对稳定。这就是 Mulliken 电荷只适合在**同一基组内做定性比较**的原因。

Mulliken 输出格式（示意输出，数值仅供参考）：

```text
 Mulliken charges:
               1
     1  C   -0.102345
     ...
     4  N   -0.431212
 Sum of Mulliken charges =   0.00000
 Mulliken charges with hydrogens summed into heavy atoms:
```

> [!TIP]
> `Pop=Orbitals=3` 会打印最高 3 个占据轨道和最低 3 个空轨道的原子贡献（格式示意：`Alpha occ 21 OE=-0.xxx is N4-p=0.xx C3-p=0.xx ...`，即轨道序号、轨道能和各原子各角动量的贡献），是快速判断 HOMO/LUMO 定域在哪里的便捷方法。更精细的轨道成分分析见第 9 节的 Multiwfn。

## 5. 实战二：吡啶的 NBO 分析

### 5.1 完整输入

```text
! Chapter 12 example: pyridine, full Mulliken population and built-in NBO 3.1 analysis
! Step 1: opt+freq   Step 2: single point with Pop=(Full,NBO)
%chk=pyridine_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

Pyridine opt+freq, B3LYP-D3(BJ)/def2-SVP

0 1
 C      -0.260785     1.152220     0.077153
 C       1.064933     0.732616     0.042780
 C       1.320285    -0.626437    -0.047956
 N       0.355378    -1.570155    -0.105138
 C      -0.922176    -1.132515    -0.069824
 C      -1.276943     0.204101     0.019943
 H      -0.500038     2.209301     0.147935
 H       1.879687     1.446327     0.085592
 H       2.338357    -1.003716    -0.077975
 H      -1.678817    -1.910312    -0.117148
 H      -2.319881     0.498569     0.044638

--Link1--
%oldchk=pyridine_opt.chk
%chk=pyridine_nbo.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2TZVP Pop=(Full,NBO) Geom=Check Guess=Read

Pyridine NBO analysis, B3LYP/def2-TZVP

0 1

```

| 关键词 | 说明 |
|---|---|
| `Pop=(Full,NBO)` | 打印全部 MO 和完整 Mulliken 分析，并调用内置 NBO 3.1 |
| `B3LYP/def2TZVP` | NBO 对基组不敏感，但建议使用合理的三重 ζ 基组 |
| 两步分开 | 优化用较小基组；分析单点放在 Link1 第 2 步，避免 NBO 在优化首末步各输出一次 |

### 5.2 输出导览

NBO 部分很长，按标题定位：

```bash
grep -n "Summary of Natural Population Analysis" pyridine_pop_nbo.log
grep -A20 "Summary of Natural Population Analysis" pyridine_pop_nbo.log
grep -n "Natural Electron Configuration" pyridine_pop_nbo.log
grep -n "NATURAL BOND ORBITAL ANALYSIS" pyridine_pop_nbo.log
grep -A40 "Natural Bond Orbitals (Summary)" pyridine_pop_nbo.log | head -60
grep -n "SECOND ORDER PERTURBATION THEORY ANALYSIS" pyridine_pop_nbo.log
```

（示意输出，数值仅供参考）

```text
 Summary of Natural Population Analysis:

                                       Natural Population
                Natural  -----------------------------------------------
    Atom  No    Charge         Core      Valence    Rydberg      Total
 -----------------------------------------------------------------------
      C    1   -0.2xxxx      1.99xxx     4.2xxxx    0.0xxxx     6.2xxxx
      ...
      N    4   -0.4xxxx      1.99xxx     5.4xxxx    0.0xxxx     7.4xxxx
```

需要关注的内容：

- **NPA 电荷**（Natural Charge 列）：比 Mulliken 更稳定，可用于讨论取代基效应、电荷转移；
- **Natural Electron Configuration**：如 N 的 `[core]2S( 1.3)2p( 4.1)` 形式，看价层电子分布；
- **Lewis 结构质量**：`Total Lewis` 占总电子数的百分比，有机分子通常 > 98%，明显偏低说明单一 Lewis 结构描述不充分（共振、离域）；
- **NBO 列表**：`LP ( 1) N 4` 是吡啶 N 的孤对，其杂化成分（如 `s( 3x%)p 1.xx`）说明它是 sp² 型、位于环平面内，不参与 π 共轭——这正是吡啶碱性的来源；`BD ( 2) C - N` 等为 π 键；
- 想直观看 NBO 形状：用 `Pop=(NBO,SaveNBOs)` 把 NBO 写入 chk，formchk 后在 GaussView 中按轨道编号显示。

## 6. 实战三：N-甲基乙酰胺的 Wiberg 键级与 E(2)

### 6.1 NBORead 输入

`Pop=NBORead` 在完整 NBO 分析的基础上，从输入流读取 `$NBO ... $END` 关键词。它放在分子说明（`Geom=Check` 时为空）之后的空行之后：

```text
! Chapter 12 example: N-methylacetamide (trans), NBO with extra $NBO keywords
! BNDIDX prints Wiberg bond indices; the E(2) donor-acceptor table is printed by default
%chk=nma_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

N-methylacetamide (trans) opt+freq, B3LYP-D3(BJ)/def2-SVP

0 1
 C      -1.883360     0.045696     0.103455
 C      -0.469513     0.551075    -0.002799
 O      -0.202773     1.747909     0.018523
 N       0.476371    -0.443387    -0.108308
 C       1.882259    -0.127160    -0.122652
 H      -2.224595     0.151443     1.136483
 H      -1.959161    -1.004416    -0.192281
 H      -2.524152     0.633729    -0.559475
 H       0.203337    -1.397028     0.088851
 H       2.428027    -0.987259    -0.517200
 H       2.204837     0.077985     0.901529
 H       2.068723     0.751413    -0.746128

--Link1--
%oldchk=nma_opt.chk
%chk=nma_nbo.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2TZVP Pop=NBORead Geom=Check Guess=Read

N-methylacetamide NBO, Wiberg bond index and second-order perturbation analysis

0 1

$NBO BNDIDX $END

```

常用的 NBO 3.1 关键词：`BNDIDX`（Wiberg 键级等键指数）；`ARCHIVE FILE=name`（写出 `name.47` 归档文件，可交给独立的 GenNBO 6/7 程序继续分析，例如 `$NBO ARCHIVE FILE=nma $END`）。更多关键词与 NBO 6/7 的新功能请查阅 [NBO 7 手册](https://nbo7.chem.wisc.edu/)。

### 6.2 读取 Wiberg 键级

```bash
grep -A20 "Wiberg bond index matrix in the NAO basis" nboread_bndidx_e2.log | head -24
```

Wiberg 键级矩阵按原子对列出。酰胺的 C-N 键级明显大于 1（部分双键性），C=O 键级小于 2，体现 N 孤对向羰基的离域（酰胺共振）。

### 6.3 读取二阶微扰 E(2) 表

```bash
grep -A60 "SECOND ORDER PERTURBATION THEORY ANALYSIS" nboread_bndidx_e2.log | grep -E "LP|Threshold" | head -30
```

（示意输出，数值仅供参考）

```text
 SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS

     Threshold for printing:   0.50 kcal/mol
                                                          E(2)  E(j)-E(i) F(i,j)
         Donor NBO (i)                     Acceptor NBO (j)       kcal/mol   a.u.    a.u.
 ===================================================================================================
   ...
  LP ( 1) N   4                 / 62.xx  BD*( 2) C   2 - O   3       6x.xx    0.2x    0.11x
  LP ( 2) O   3                 / 50.xx  BD*( 1) C   1 - C   2       2x.xx    0.6x    0.10x
```

解读：

- `LP(1) N4 → BD*(2) C2-O3` 即 $n_N\to\pi^*_{C=O}$，这是酰胺共振的 NBO 表述，E(2) 通常达数十 kcal/mol；
- `LP(2) O3 → BD*(1) C1-C2` 等为羰基氧孤对对相邻 σ* 的负超共轭；
- 只打印 ≥ 0.50 kcal/mol 的作用（阈值见表头）；氢键复合物中分子间的 $n\to\sigma^*_{\rm X-H}$ 作用也在这张表中，是判断氢键强弱的常用依据。

> [!WARNING]
> E(2) 是二阶微扰估计，各项**不能简单相加**得到“总共轭能”，其数值也依赖基组和泛函。需要定量的离域能时，用 `Pop=NBODel`（删除指定的 Fock 矩阵元后重新计算能量，仅限 SCF 方法）或能量分解方法，并只做同一水平下的相对比较。

### 6.4 外部 NBO 7（C.01 起）

- `Pop=NBO7`、`NPA7`、`NBO7Read`、`NBO7Delete` 通过外部接口调用 NBO 7，需要单独从 [NBO 7 官网](https://nbo7.chem.wisc.edu/) 获得授权并按其说明安装与 Gaussian 的接口脚本；
- NBO 7 提供 NRT（自然共振理论）、STERIC、CMO 等 NBO 3.1 没有的分析，也支持 NEDA 能量分解。C.01 发行说明：NEDA 通过 `NEDA=n` 请求，碎片定义使用与 Counterpoise 计算相同的语法（原子后标注 `(Fragment=n)`，电荷多重度行写总体加各碎片，见 [第14章](14-weak-bsse.md)），具体写法请以 NBO 7 手册为准；
- 没有 NBO 7 授权时，可以用 `$NBO ARCHIVE FILE=... $END` 写出 .47 文件，在有授权的机器上用 GenNBO 继续分析；[Molden2AIM](https://github.com/zorkzou/Molden2AIM) 也能把其他程序的 molden 文件转成 NBO .47 格式。

## 7. 实战四：静电势拟合电荷与 AMBER 风格 RESP

### 7.1 为什么是 HF/6-31G(d)

AMBER 力场（ff94 及以后）和 GAFF/GAFF2 的电荷协议采用**气相 HF/6-31G\***（即 6-31G(d)）静电势拟合 RESP 电荷。HF/6-31G* 会系统地高估分子偶极矩，这恰好近似地模拟了凝聚相中的极化效应，与力场其余参数是配套的。因此为 AMBER/GAFF 体系拟合电荷时，**不要**随意改成 B3LYP/def2-TZVP 等“更高”的水平；用于其他力场或其他目的（如 RESP2）时另当别论。

### 7.2 完整输入

示例目录：[12-esp-resp-paracetamol](../examples/12-esp-resp-paracetamol/README.md)，配体为对乙酰氨基酚（paracetamol）。

```text
! Chapter 12 example: AMBER-style ESP for RESP fitting of a small ligand (paracetamol)
! Step 1: HF/6-31G(d) opt (classic AMBER/GAFF protocol)
! Step 2: HF/6-31G(d) MK ESP. 6/33=2 prints ESP points, 6/41=10 layers, 6/42=17 point density,
!         6/50=1 writes an Antechamber .gesp file; its name is the last line of the step
%chk=lig_opt.chk
%mem=8GB
%nprocshared=8
#p HF/6-31G(d) Opt SCF=Tight

Paracetamol HF/6-31G(d) geometry optimization

0 1
 C       3.691373    -0.360571     0.117591
 C       2.324738    -0.235562    -0.509254
 O       2.189776    -0.267602    -1.727077
 N       1.324123    -0.118925     0.438503
 C      -0.061788     0.011747     0.213336
 C      -0.886809     0.103826     1.342915
 C      -2.270808     0.236126     1.205105
 C      -2.833598     0.277054    -0.062634
 O      -4.186799     0.406769    -0.155051
 C      -2.033565     0.188000    -1.195089
 C      -0.647108     0.055564    -1.054623
 H       3.775184     0.274994     1.003877
 H       3.868912    -1.403154     0.393946
 H       4.451223    -0.040538    -0.601250
 H       1.612526    -0.136313     1.408452
 H      -0.468575     0.074044     2.345249
 H      -2.905656     0.306769     2.083621
 H      -4.431592     0.419296    -1.095143
 H      -2.460147     0.218961    -2.192773
 H      -0.051411    -0.010485    -1.959704

--Link1--
%oldchk=lig_opt.chk
%chk=lig_esp.chk
%mem=8GB
%nprocshared=8
#p HF/6-31G(d) SCF=Tight Geom=Check Guess=Read NoSymm
   Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)

Paracetamol HF/6-31G(d) Merz-Kollman ESP points for RESP

0 1

paracetamol.gesp

```

| 关键词 / IOp | 含义 | 出处 |
|---|---|---|
| `Pop=MK` | Merz-Singh-Kollman 格点方案拟合 ESP 电荷 | 官方 Population 页 |
| `IOp(6/33=2)` | 在输出中打印拟合格点坐标与静电势值（antechamber/espgen 读取 log 时需要） | AMBER 社区惯用写法（当前官方 IOp 列表中未列出） |
| `IOp(6/41=10)` | ESP 拟合格点层数为 10（默认 4） | 官方 Overlay 6 IOp 说明 |
| `IOp(6/42=17)` | 每层格点密度（默认 1），17 对应每个原子约 2500 个点，R.E.D. 工具采用的设置 | 官方 Overlay 6 IOp 说明；R.E.D. FAQ |
| `IOp(6/50=1)` | 拟合时写出 Antechamber 数据文件，文件名在本步输入末尾单独一行给出 | 官方 Population 页 |
| `NoSymm` | 保持输入取向，避免对称化带来的坐标变化 | — |
| `SCF=Tight` | 更严格的 SCF 收敛，保证 ESP 精度 | — |

antechamber 自动生成的输入常用 `IOp(6/33=2,6/42=6)`（每个原子较少的点），也完全可用；重要的是同一项目中保持一致。

### 7.3 用 antechamber 拟合两步 RESP

需要 AmberTools（`antechamber` 和 `resp` 程序）。仓库示例里的 `run_antechamber_resp.sh` 封装了下面的命令：

```bash
# 方式 1：读取 Gaussian 写出的 .gesp 文件（推荐）
antechamber -i paracetamol.gesp -fi gesp -o paracetamol_resp.mol2 -fo mol2 \
            -c resp -nc 0 -rn LIG -at gaff2
# 方式 2：直接读取 Gaussian 输出（需要 IOp(6/33=2) 打印的 ESP 点）
antechamber -i lig_esp.log -fi gout -o paracetamol_resp.mol2 -fo mol2 \
            -c resp -nc 0 -rn LIG -at gaff2
```

- `-c resp`：两步 RESP 拟合；`-nc 0`：分子净电荷；`-at gaff2`：原子类型；
- 得到的 mol2 中最后一列即 RESP 电荷，可再用 `parmchk2` 检查缺失参数（属于 AMBER 建模流程，不在本章展开）；
- 对照项目 [make_resp](https://github.com/alexamist/make_resp) 演示了从 Gaussian MK 静电势到 RESP 的完整小流程。

### 7.4 用 Multiwfn 拟合 RESP

Multiwfn 的 RESP 模块（主功能 7 → 子功能 18）可直接读取 fchk，自己计算格点上的静电势：

```bash
formchk lig_esp.chk lig_esp.fchk
Multiwfn lig_esp.fchk
#  7    -> 布居分析与原子电荷
#  18   -> RESP 模块
#  1    -> 标准两步 RESP 拟合（Kollman 等的定义）
```

Multiwfn 的 RESP 模块还支持：多构象联合拟合（选项 −1，提供构象列表与权重）、自定义等价约束和电荷约束、改用 CHELPG 格点等，详见 [Multiwfn 手册 3.9.16 节](http://sobereva.com/multiwfn/manual/)。

### 7.5 内置 MK/CHelpG 与特殊情况

- 普通 ESP 电荷：直接 `Pop=MK` 或 `Pop=CHelpG`，输出在 `ESP charges:` 之后：

```bash
grep -A25 "ESP charges:" lig_esp.log | head -27
grep "Charges from ESP fit" lig_esp.log        # 含拟合 RMS，越小越好
```

- `Pop=(MK,Dipole)`：约束电荷重现量子化学偶极矩；
- **含金属等缺少默认半径的元素**：MK 默认半径不覆盖全周期表，可用 `Pop=MKUFF`（UFF 半径），或 `Pop=(MK,ReadRadii)` 在输入末尾按“元素 半径(Å)”成对给出，空行结束，例如：

```text
Zn 1.40

```

这里的 1.40 只是格式示例，半径取值请依据文献（如所用力场的约定）确定；
- **构象依赖性**：ESP 电荷依赖构象，柔性分子最好对多个构象联合拟合（Multiwfn 选项 −1 或 R.E.D. 工具）；
- **埋藏原子**：分子内部的原子离格点远，其 ESP 电荷数值不稳定，这正是 RESP 约束要解决的问题。

## 8. 实战五：Hirshfeld/CM5、MK 与 NPA 的对比（醋酸根）

### 8.1 用 Density=Check 复用同一个密度

官方 Population 页建议：一个作业步只能用一种密度和一种拟合方法，需要多种分析时，追加 `Density=Check` 的作业步（隐含 `Guess=Only ChkBasis`），不必重做 SCF。示例 `acetate_hirshfeld_cm5.gjf`：

```text
! Chapter 12 example: acetate anion in water, Hirshfeld/CM5 vs MK vs NPA charges
! Step 1: opt+freq (SMD water) with Pop=Hirshfeld
! Steps 2-3: re-analyse the SAME density read from the chk (Density=Check implies Guess=Only)
%chk=acetate.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-311+G(d,p) EmpiricalDispersion=GD3BJ Opt Freq SCRF=(SMD,Solvent=Water)
   Pop=Hirshfeld

Acetate anion opt+freq, B3LYP-D3(BJ)/6-311+G(d,p), SMD(water), Hirshfeld/CM5

-1 1
 C      -0.639020    -0.032922    -0.033364
 C       0.875403     0.051535     0.044772
 O       1.457987    -1.051972     0.239774
 O       1.358682     1.209230    -0.094474
 H      -0.931450    -0.691546    -0.856605
 H      -1.085127     0.951767    -0.202960
 H      -1.036475    -0.436090     0.902858

--Link1--
%oldchk=acetate.chk
%chk=acetate_mk.chk
%mem=8GB
%nprocshared=8
#p B3LYP Geom=AllCheck Density=Check Pop=MK

--Link1--
%oldchk=acetate.chk
%chk=acetate_npa.chk
%mem=8GB
%nprocshared=8
#p B3LYP Geom=AllCheck Density=Check Pop=NPA

```

- 醋酸根带一个负电荷，阴离子需要弥散函数（6-311+G(d,p)），并用 SMD(水) 模拟溶液环境；
- `Geom=AllCheck` 从 chk 读取标题、电荷多重度和结构，因此第 2、3 步没有标题和分子说明；
- `Density=Check` 读取第 1 步保存的密度（溶液中的 SCF 密度），基组也从 chk 读取，route 中不写基组。

### 8.2 读取 Hirshfeld/CM5 输出

```bash
grep -A10 "Hirshfeld charges, spin densities, dipoles, and CM5 charges" acetate_hirshfeld_cm5.log | tail -11
grep -A9  "ESP charges:" acetate_hirshfeld_cm5.log | tail -9
grep -A12 "Summary of Natural Population Analysis" acetate_hirshfeld_cm5.log | tail -9
```

（示意输出，数值仅供参考）

```text
 Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=      4:
                Q-H        S-H        Dx         Dy         Dz        Q-CM5
     1  C   -0.0xxxxx   0.000000  ...                                -0.2xxxxx
     2  C    0.1xxxxx   0.000000  ...                                 0.1xxxxx
     3  O   -0.5xxxxx   0.000000  ...                                -0.5xxxxx
```

`Q-H` 为 Hirshfeld 电荷，`Q-CM5` 为 CM5 电荷，`S-H` 为 Hirshfeld 自旋布居（闭壳层为 0）。注意 Opt+Freq 作业中布居分析会出现多次（优化首步、末步和 Freq 步），取最后一次。

比较时你通常会看到：Hirshfeld 电荷绝对值最小；CM5 与 MK 的量级相近；NPA 电荷绝对值较大；两个羧酸 O 在各种方法下都应接近等价。把结果整理成表，体会“电荷依定义而异”。

### 8.3 为 ONIOM/MM 保存电荷

C.01 起，`Pop=SaveCM5` 与 `Pop=SaveHirshfeld` 把对应电荷作为 MM 电荷存入 chk，供后续带 `Geom=Check` 的 MM 或 ONIOM 计算使用（类似的还有 `SaveMulliken`、`SaveESP`、`SaveNPA`）。为确保电荷被计算，可以与计算选项一起写，例如 `Pop=(Hirshfeld,SaveCM5)`、`Pop=(MK,SaveESP)`。多层 ONIOM 中默认只保存 QM 层中显式计算的电荷；`Uncharged` 选项可保留输入中已有的电荷。详见 [第15章](15-oniom.md)。

## 9. Multiwfn 中的电荷、键级与轨道成分

先把 chk 转成 fchk：`formchk pyridine_nbo.chk pyridine_nbo.fchk`，再 `Multiwfn pyridine_nbo.fchk`。常用菜单（编号见 [Multiwfn 手册](http://sobereva.com/multiwfn/manual/) 第 3.9–3.11 节）：

| 主功能 → 子功能 | 内容 |
|---|---|
| 7 → 1 | Hirshfeld 电荷 |
| 7 → 11 | ADCH 电荷（原子偶极矩校正的 Hirshfeld，能较好重现分子偶极和 ESP，Lu & Chen 2012） |
| 7 → 12 / 13 | CHELPG / MK 静电势拟合电荷 |
| 7 → 16 | CM5 电荷 |
| 7 → 18 | RESP 电荷 |
| 9 → 1 | Mayer 键级 |
| 9 → 3 | Löwdin 正交化基下的 Wiberg 键级 |
| 9 → 7 | 模糊键级（fuzzy bond order） |
| 8 → 7 / 8 / 9 | 轨道成分分析：NAO / Hirshfeld / Becke 方法 |

- 轨道成分分析能回答“HOMO 中金属 d 轨道占多少”“LUMO 在哪个片段上”这类问题，比 `Pop=Orbitals` 的粗略输出更可靠；
- 批量处理：仓库脚本 `batch_cube_multiwfn.sh` 可对目录内所有 chk 执行 formchk/cubegen，并可用命令流文件非交互运行 Multiwfn（`bash ../scripts/batch_cube_multiwfn.sh -h`）；[BaneWfn](https://github.com/bane-dysta/banewfn) 把固定菜单序列模板化，适合大批量分析；
- 分子表面静电势（ESP 着色、极值点）分析见 [第16章](16-postprocess-multiwfn.md)。

## 10. 该用哪种电荷？

| 用途 | 推荐 | 不推荐 |
|---|---|---|
| AMBER/GAFF 力场参数化 | HF/6-31G(d) 两步 RESP | Mulliken、高水平 DFT 的 MK（与力场不配套） |
| CHARMM/OPLS 等其他力场 | 按该力场自身协议 | 跨力场照搬 RESP |
| 讨论取代基效应、电荷转移趋势 | NPA、CM5、ADCH | 弥散基组下的 Mulliken |
| 估算偶极矩、静电相互作用 | CM5、ADCH、MK/CHelpG | Hirshfeld 原始值（偏小） |
| 反应活性位点（亲核/亲电） | ESP 表面极值（第16章）、Fukui 函数 | 只看单个原子电荷 |
| ONIOM 电子嵌入的 MM 电荷 | RESP 或 CM5（`SaveCM5`） | 未经检查的 Mulliken |
| 快速定性比较（同一基组） | Mulliken 也可 | 跨基组比较 Mulliken |

## 参数选择与收敛建议

| 项目 | 推荐 | 说明 |
|---|---|---|
| 结构 | 先优化并确认极小点 | 电荷依赖结构与构象 |
| 分析水平 | 与研究主体一致的泛函，三重 ζ 基组 | RESP 例外：AMBER 协议固定为 HF/6-31G(d) |
| 阴离子 | 含弥散函数 | 但 Mulliken 在弥散基组下不可靠，改用 NPA/CM5 |
| 溶剂 | 讨论溶液性质时用 SMD/PCM 下的密度 | RESP（AMBER）使用气相 |
| ESP 格点 | MK + `IOp(6/41=10,6/42=17)` 或 antechamber 默认 | 同一项目保持一致 |
| 多种电荷 | 追加 `Geom=AllCheck Density=Check Pop=...` 步骤 | 避免重复 SCF |
| NBO | 内置 3.1 满足大多数需求；NRT/NEDA 等需要 NBO 7 | E(2) 只做相对比较 |
| 键级 | NBO Wiberg（`BNDIDX`）或 Multiwfn Mayer | 不同定义的数值不可混比 |

## 常见问题与报错

**1. 写了 `Pop` 但没有选项，作业在输入解析阶段就报错退出**
原因：官方说明 Population 关键词必须带选项。
解决：写成 `Pop=Regular`、`Pop=NBO` 等。

**2. `$NBO ... $END` 似乎没有生效，或者作业在读取输入时出错**
原因：NBO 关键词行前缺少空行，或在 `Pop=NBO`（而不是 `Pop=NBORead`）下写了 $NBO 段。
解决：使用 `Pop=NBORead`；$NBO 行放在分子说明（或 `Geom=Check` 时的电荷多重度行）后的空行之后，并以空行结束输入。

**3. 使用 `Pop=NBO7` 等选项时作业无法调用 NBO 7**
原因：没有安装外部 NBO 7 或 Gaussian 找不到其接口脚本。
解决：按 NBO 7 官方说明安装并配置；没有授权时改用内置 `Pop=NBO`，或写出 .47 文件在其他机器上分析。

**4. Mulliken 电荷出现大于 1 或符号反常的数值**
原因：弥散函数使基函数归属失去物理意义。
解决：改用 NPA、CM5 或 ADCH，不要在含弥散函数的基组下解释 Mulliken 电荷。

**5. 含过渡金属体系做 `Pop=MK` 时因缺少原子半径而失败**
原因：MK 默认半径没有覆盖所有元素。
解决：`Pop=MKUFF`，或 `Pop=(MK,ReadRadii)` 并在输入末尾按“元素 半径”给出半径（空行结束），在论文中注明半径来源。

**6. antechamber 读取 log 时找不到 ESP 数据或 RESP 电荷异常**
原因：没有加 `IOp(6/33=2)`（log 中无格点数据）；多步 Link1 的 log 中包含多个结构；净电荷 `-nc` 与 Gaussian 中的电荷不一致。
解决：优先使用 `IOp(6/50=1)` 生成的 .gesp 文件；确认 `-nc` 与输入一致；检查拟合 RMS。

**7. 同一分子两次计算得到的 RESP 电荷差别较大**
原因：构象不同（例如 OH 或甲基取向不同），ESP 电荷对构象敏感。
解决：固定一组有代表性的构象，或做多构象联合拟合。

**8. 优化作业中读到的电荷与单点结果不同**
原因：Opt 作业的布居分析出现在第一步和最后一步，Opt Freq 中 Freq 步还会再打印一次，读到的是初始结构的结果。
解决：用 `grep ... | tail` 读取最后一次输出，或者像本章示例一样把分析放在单独的单点步骤。

**9. `Erroneous write` 或磁盘相关错误出现在 NBO/ESP 步骤**
原因：scratch 空间不足（大体系的格点和 NBO 中间文件较大）。
解决：为 `GAUSS_SCRDIR` 指定足够大的本地盘，参见 [第17章](17-hpc-jobs.md) 与 [第18章](18-errors-capstone.md)。

## 练习

1. 运行 `pyridine_mulliken_basis.gjf`，把 5 个基组下 N 原子的 Mulliken 与 CM5 电荷画成折线图，讨论二者的基组依赖性。
2. 对 N-甲基乙酰胺的顺式（cis）构象重复 NBO 计算，比较 $n_N\to\pi^*_{C=O}$ 的 E(2) 和 C-N Wiberg 键级。
3. 分别用 antechamber 和 Multiwfn 拟合 paracetamol 的 RESP 电荷，比较两组结果（注意等价约束的设置）。
4. 对醋酸根分别在气相和 SMD(水) 中计算 CM5 电荷，观察溶剂极化对氧原子电荷的影响。
5. 用 Multiwfn 计算吡啶的 ADCH 电荷和 Mayer 键级，与 NPA 电荷、NBO Wiberg 键级对比。

## 参考资料

官方与工具文档：

- [Population 关键词](https://gaussian.com/population/)、[Density 关键词](https://gaussian.com/density/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes（NBO7、Resp=N、SaveCM5/SaveHirshfeld）](https://gaussian.com/relnotes/)
- [IOp Overlay 6 说明（6/41、6/42、6/50）](https://gaussian.com/overlay6/)
- [NBO 7 官网与手册](https://nbo7.chem.wisc.edu/)
- [AmberTools / antechamber 文档](https://ambermd.org/AmberTools.php)
- [R.E.D. FAQ（ESP 格点设置的讨论）](https://upjv.q4md-forcefieldtools.org/RED/FAQ-I.php)
- [Multiwfn 主页](http://sobereva.com/multiwfn/) 与 [在线手册](http://sobereva.com/multiwfn/manual/)

相关 GitHub 仓库（仅链接）：

- [alexamist/make_resp](https://github.com/alexamist/make_resp)：Gaussian MK 静电势 + antechamber 拟合 RESP 的小流程
- [zorkzou/Molden2AIM](https://github.com/zorkzou/Molden2AIM)：molden 转 wfn/wfx/NBO .47
- [cclib/cclib](https://github.com/cclib/cclib)：解析 Mulliken、NPA 等原子电荷
- [theochem/iodata](https://github.com/theochem/iodata)：fchk、wfx、molden、cube 等格式读写
- [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn)：Multiwfn 批处理工作流
- [1234zou/MOKIT](https://github.com/1234zou/MOKIT)：在 Gaussian 与其他程序之间传递轨道

关键文献：

- R. S. Mulliken, *J. Chem. Phys.* 23, 1833 (1955). DOI: [10.1063/1.1740588](https://doi.org/10.1063/1.1740588)
- A. E. Reed, L. A. Curtiss, F. Weinhold, *Chem. Rev.* 88, 899 (1988). DOI: [10.1021/cr00088a005](https://doi.org/10.1021/cr00088a005)（NBO）
- E. D. Glendening, C. R. Landis, F. Weinhold, *J. Comput. Chem.* 40, 2234 (2019). DOI: [10.1002/jcc.25873](https://doi.org/10.1002/jcc.25873)（NBO 7.0）
- B. H. Besler, K. M. Merz, P. A. Kollman, *J. Comput. Chem.* 11, 431 (1990). DOI: [10.1002/jcc.540110404](https://doi.org/10.1002/jcc.540110404)（MK）
- C. M. Breneman, K. B. Wiberg, *J. Comput. Chem.* 11, 361 (1990). DOI: [10.1002/jcc.540110311](https://doi.org/10.1002/jcc.540110311)（CHelpG）
- C. I. Bayly, P. Cieplak, W. Cornell, P. A. Kollman, *J. Phys. Chem.* 97, 10269 (1993). DOI: [10.1021/j100142a004](https://doi.org/10.1021/j100142a004)（RESP）
- J. Wang, W. Wang, P. A. Kollman, D. A. Case, *J. Mol. Graph. Model.* 25, 247 (2006). DOI: [10.1016/j.jmgm.2005.12.005](https://doi.org/10.1016/j.jmgm.2005.12.005)（antechamber）
- F. L. Hirshfeld, *Theor. Chim. Acta* 44, 129 (1977). DOI: [10.1007/BF00549096](https://doi.org/10.1007/BF00549096)
- A. V. Marenich, S. V. Jerome, C. J. Cramer, D. G. Truhlar, *J. Chem. Theory Comput.* 8, 527 (2012). DOI: [10.1021/ct200866d](https://doi.org/10.1021/ct200866d)（CM5）
- T. Lu, F. Chen, *J. Theor. Comput. Chem.* 11, 163 (2012). DOI: [10.1142/S0219633612500113](https://doi.org/10.1142/S0219633612500113)（ADCH）
- T. Lu, F. Chen, *J. Comput. Chem.* 33, 580 (2012). DOI: [10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)（Multiwfn）

[← 上一章](11-nmr-props.md) | [返回目录](../README.md) | [下一章 →](13-openshell-stability.md)
