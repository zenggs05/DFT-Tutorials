# 第9章 隐式溶剂模型 PCM/SMD（Implicit Solvation with PCM & SMD in Gaussian 16）

本章讲解 Gaussian 16 中隐式溶剂模型的用法：`SCRF=(IEFPCM,Solvent=...)`、`SCRF=CPCM` 和 `SCRF=(SMD,Solvent=...)`，包括溶剂化自由能 ΔG_solv 的计算、自定义溶剂 `Solvent=Generic,Read`、溶液相 Opt Freq 与“气相结构 + 溶剂单点”两种策略的取舍、激发态的非平衡溶剂化与态特定 `ExternalIteration`，以及 pKa 和氧化还原电位两个热力学循环实例。关键词：Gaussian SMD solvation free energy、IEFPCM、CPCM、implicit solvent、pKa calculation、redox potential。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 写出 `SCRF=(SMD,Solvent=Water)`、`SCRF=(IEFPCM,Solvent=Toluene)`、`SCRF=(CPCM,Solvent=Acetonitrile)`，知道 G16 的默认 PCM 是什么、与 G09 有何不同；
2. 用“气相 + SMD 单点”计算溶剂化自由能 ΔG_solv，并知道 SMD 原始参数化使用了哪些方法；
3. 在“溶液中 Opt Freq”与“气相 Opt Freq + 溶剂单点”两种策略之间做出合理选择，并正确处理 1 atm → 1 M 标准态校正；
4. 用 `Solvent=Generic,Read` 定义列表中没有的溶剂，并会调整 PCM 空腔参数（Radii、Alpha、Surface）；
5. 理解平衡与非平衡溶剂化，知道 TD-DFT 中默认用哪一种，以及 `ExternalIteration`、`NonEquilibrium=Save/Read` 的用法；
6. 用热力学循环计算 pKa（苯酚）和氧化还原电位（四硫富瓦烯 TTF，乙腈中）；
7. 处理 `Inv3 failed in PCMMkU`、`GetVDW: no radius for atom` 等溶剂相关报错。

**适用场景**：溶液中的反应能垒与反应自由能（离子、极性 TS 尤其重要）、溶剂化自由能与分配系数、pKa、氧化还原电位、溶剂中的 UV-Vis/NMR 等光谱性质。

**局限**：隐式溶剂不描述**特定的**溶质–溶剂相互作用（强氢键、配位、质子穿梭）。这类情况要采用“显式溶剂分子 + 连续介质”的簇–连续介质（cluster–continuum）模型，或者做分子动力学。

**前置知识**：几何优化（[第5章](05-opt.md)）、频率与热化学、标准态校正（[第6章](06-freq-thermo.md)）。

本章配套示例：[../examples/09-solvation-smd/README.md](../examples/09-solvation-smd/README.md)

| 文件 | 内容 |
|---|---|
| `phenol_smd_water_opt_freq.gjf` | 苯酚在 SMD(水) 中 Opt Freq |
| `phenoxide_smd_water_opt_freq.gjf` | 苯氧负离子在 SMD(水) 中 Opt Freq（与上一个文件配合计算 pKa） |
| `solvation_free_energy_gas_and_smd.gjf` | Link1：气相 Opt Freq → 同一结构的 SMD(水) 单点，得到 ΔG_solv |
| `custom_solvent_generic_read.gjf` | Link1：内置乙醇 vs 手动输入乙醇参数的 `Generic,Read`，验证自定义溶剂写法 |
| `redox_potential_cycle_neutral_and_cation.gjf` | Link1：TTF 与 TTF•⁺ 在 SMD(乙腈) 中 Opt Freq + 高水平单点 |

---

## 2. 原理速览

### 2.1 连续介质模型的基本图像

隐式溶剂模型把溶质放进一个由原子球叠加而成的**分子空腔**中，空腔外是介电常数为 ε 的连续介质。溶质的电荷分布使介质极化，在空腔表面感应出**表观表面电荷**（apparent surface charge），这些电荷产生的**反应场**又反过来极化溶质波函数，二者自洽求解。

最简单的情形是半径为 $a$ 的球形空腔中的点电荷 $q$（Born 模型）：

$$
\Delta G_{\text{el}} = -\frac{q^2}{2a}\left(1-\frac{1}{\varepsilon}\right)
$$

它说明了两件事：离子的溶剂化能远大于中性分子（与 $q^2$ 成正比）；ε 从 1 增大到约 10 时，大部分静电溶剂化效应已经出现，ε 再增大，变化就不大了。

### 2.2 IEFPCM、CPCM 与 SMD

| 模型 | Gaussian 写法 | 包含的项 | 适合 |
|---|---|---|---|
| IEFPCM | `SCRF=PCM`（`IEFPCM` 为同义词），**SCRF 的默认方法** | 只有静电项（默认不含空穴形成、色散、排斥） | 溶剂对结构、光谱、相对能量的影响 |
| CPCM | `SCRF=CPCM` | 导体近似的静电项 | 与 IEFPCM 类似；高 ε 溶剂中两者结果非常接近 |
| SMD | `SCRF=(SMD,Solvent=...)` | IEFPCM 静电项（使用 SMD 自己的原子半径）+ CDS 项（空穴形成、色散、溶剂结构） | **溶剂化自由能 ΔG_solv**（官方文档推荐） |

SMD（Marenich, Cramer, Truhlar 2009）把溶剂化自由能写成

$$
\Delta G^{*}_{\text{solv}} = \Delta G_{\text{ENP}} + G_{\text{CDS}}
$$

其中 ENP 为电子–核–极化（即静电）项，CDS 项由原子表面张力参数与溶剂可及表面积计算。SMD 用同一套参数（原子库仑半径 + 原子表面张力系数）覆盖了 6 种电子结构方法：M05-2X/MIDI!6D、M05-2X/6-31G\*、M05-2X/6-31+G\*\*、M05-2X/cc-pVTZ、B3LYP/6-31G\* 和 HF/6-31G\*。按原文，用 6-31G\* 基组时，中性分子溶剂化自由能的平均无符号误差为 0.6–1.0 kcal/mol，离子约 4 kcal/mol。

> [!IMPORTANT]
> 默认的 IEFPCM 只有静电项，`E(IEFPCM) − E(gas)` **不是**溶剂化自由能。凡是要算 ΔG_solv、pKa、氧化还原电位、分配系数，都应该用 SMD。

### 2.3 Gaussian 16 的默认设置

| 项目 | G16 默认 | 说明 |
|---|---|---|
| SCRF 方法 | IEFPCM（对称形式） | G09 默认是非对称 IEFPCM；`SCRF=AsymmetricIEFPCM`（或 `AIEFPCM`）恢复 G09 行为 |
| 溶剂 | 水 | 不写 `Solvent=` 时默认为水 |
| 空腔原子半径（PCM/CPCM） | UFF 半径，乘以 1.1（`Alpha=1.1`），H 原子有独立的球 | SMD 使用自己的原子半径 |
| 表面类型 | `Surface=VDW`（不生成额外的平滑球） | G03 默认是 SES |
| 表面电荷 | 连续表面电荷（continuous surface charge）形式 | 保证能量和梯度对原子位置连续变化，优化更稳定 |
| 兼容旧版 | `SCRF=G09Defaults` / `G03Defaults` | 复现旧文献时使用 |

> [!WARNING]
> 拿 G16 的 PCM 结果与 G09 时代的文献直接比较时，默认积分格点（G16 为 UltraFine）和 PCM 形式（对称 vs 非对称 IEFPCM）都不同，能量差可能达到零点几 kcal/mol 甚至更多。整个项目要固定同一个版本和协议。

### 2.4 标准态

SMD 给出的 ΔG*_solv 对应**气相 1 mol/L → 溶液 1 mol/L**（Ben-Naim 标准态，带星号）。Gaussian 的气相热化学对应 **1 atm**。把气相自由能换到 1 mol/L：

$$
\Delta G^{\circ \rightarrow *} = RT\ln\frac{RT}{p^\circ V^*} = RT \ln(24.46) = 1.89\ \text{kcal/mol}\quad(298.15\ \text{K})
$$

如果溶剂本身是反应物或产物（如水），它的标准态是纯液体（55.34 mol/L），还要再加 $RT\ln(55.34) = 2.38$ kcal/mol。

---

## 3. 基本用法与输出

### 3.1 写法

```text
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(SMD,Solvent=Water) Opt Freq
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(IEFPCM,Solvent=Toluene) Opt Freq
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(CPCM,Solvent=Acetonitrile) Opt Freq
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(Solvent=DiMethylSulfoxide) Opt Freq
```

最后一行没有写模型，使用默认的 IEFPCM。

### 3.2 常用溶剂名（Gaussian 关键词，ε 取自官方溶剂表）

| 溶剂 | 关键词 | ε | 溶剂 | 关键词 | ε |
|---|---|---|---|---|---|
| 水 | `Water` | 78.3553 | 二氯甲烷 | `Dichloromethane` | 8.93 |
| 二甲亚砜 | `DiMethylSulfoxide` | 46.826 | 四氢呋喃 | `TetraHydroFuran` | 7.4257 |
| DMF | `n,n-DiMethylFormamide` | 37.219 | 氯仿 | `Chloroform` | 4.7113 |
| 乙腈 | `Acetonitrile` | 35.688 | 乙醚 | `DiethylEther` | 4.2400 |
| 甲醇 | `Methanol` | 32.613 | 甲苯 | `Toluene` | 2.3741 |
| 乙醇 | `Ethanol` | 24.852 | 苯 | `Benzene` | 2.2706 |
| 丙酮 | `Acetone` | 20.493 | 1,4-二氧六环 | `1,4-Dioxane` | 2.2099 |
| 吡啶 | `Pyridine` | 12.978 | 正己烷 | `n-Hexane` | 1.8819 |

完整列表（约 180 种溶剂）见 [SCRF 官方页面](https://gaussian.com/scrf/) 的 Solvents 部分。

> [!NOTE]
> 每种内置溶剂除 ε 外还有一整套参数（光学介电常数、溶剂半径；SMD 还包括氢键酸性、碱性、表面张力等）。官方文档特别提醒：只改 ε 并不能正确定义一个新溶剂。

### 3.3 输出中要检查什么

```bash
grep -i "solvent" job.log | head          # 确认溶剂和模型确实被使用
grep "SMD-CDS" job.log | tail -n 1        # SMD 的非静电项（已包含在总能量中）
grep "SCF Done" job.log | tail -n 1       # 溶液中的总电子能量
```

日志开头的 PCM 参数块大致如下（示意输出，具体格式以你的日志为准）：

```text
 Polarizable Continuum Model (PCM)
 =================================
 Model                : PCM.
 Atomic radii         : SMD-Coulomb.
 Solvent              : Water, Eps=  78.355300 Eps(inf)=   1.777849
```

SMD 计算会额外打印一行 CDS 能量（示意输出，数值仅供参考）：

```text
 SMD-CDS (non-electrostatic) energy       (kcal/mol) =       0.85
 (included in total energy above)
```

G16 中 `SCF Done` 的能量**已经包含了所有溶剂化贡献**（G03 的输出与此不同），不需要再手动加 CDS 项。

---

## 4. 实战一：溶液中的 Opt Freq（苯酚，SMD 水）

### 4.1 输入

完整输入 [`phenol_smd_water_opt_freq.gjf`](../examples/09-solvation-smd/phenol_smd_water_opt_freq.gjf)：

```text
%chk=phenol_smd_water.chk
%nprocshared=8
%mem=8GB
#p M062X/6-311+G(d,p) SCRF=(SMD,Solvent=Water) Opt Freq

Phenol in water (SMD): Opt Freq in solution

0 1
O      2.438305   -0.394243   -0.341688
C      1.101868   -0.192003   -0.160619
C      0.234248   -1.204037   -0.555917
C     -1.140527   -1.038919   -0.388906
C     -1.637634    0.137825    0.172537
C     -0.759877    1.149160    0.566837
C      0.616508    0.985891    0.400548
H      2.914002    0.389600   -0.022251
H      0.627566   -2.117527   -0.992261
H     -1.822487   -1.827084   -0.696335
H     -2.709139    0.266424    0.302663
H     -1.149961    2.064699    1.003880
H      1.287128    1.780214    0.711512

```

| 关键词 | 说明 |
|---|---|
| `M062X/6-311+G(d,p)` | M06-2X 配合带弥散函数的三重 ζ 基组，是 SMD 计算 pKa 时常见的组合（苯氧负离子需要弥散函数）。Minnesota 泛函对积分格点敏感，G16 默认的 UltraFine 格点已经足够 |
| `SCRF=(SMD,Solvent=Water)` | SMD 模型，溶剂为水。优化和频率都在溶剂中进行，属于平衡溶剂化 |
| `Opt Freq` | 溶液中的结构优化和频率。得到的 `Sum of electronic and thermal Free Energies` 是**溶液中**的自由能（1 atm 理想气体热校正） |

### 4.2 运行与检查

```bash
g16 < phenol_smd_water_opt_freq.gjf > phenol_smd_water_opt_freq.log
grep -c "Normal termination" phenol_smd_water_opt_freq.log           # 2
grep "imaginary frequencies" phenol_smd_water_opt_freq.log           # 无输出 = 0 个虚频
grep "Sum of electronic and thermal Free Energies" phenol_smd_water_opt_freq.log
```

溶液中的自由能（1 M 标准态）：

$$
G^{*}_{\text{sol}} = G_{\text{log}} + 1.89\ \text{kcal/mol}
$$

其中 $G_{\text{log}}$ 就是上面 grep 得到的那一行。

> [!TIP]
> 溶液中的频率计算偶尔会出现很小的虚频（几 i 到几十 i cm⁻¹），常见于甲基、羟基等转子，原因是空腔离散化带来的数值噪声。先用 `Opt=Tight` 收紧收敛判据，必要时加 `Int=SuperFine`。如果虚频仍然存在且模式是真实的转动，就沿该模式位移后重新优化（[第6章](06-freq-thermo.md)）。

---

## 5. 实战二：溶剂化自由能 ΔG_solv

### 5.1 输入（Link1：气相 → SMD 单点）

完整输入 [`solvation_free_energy_gas_and_smd.gjf`](../examples/09-solvation-smd/solvation_free_energy_gas_and_smd.gjf)：

```text
%chk=phenol_dgsolv.chk
%nprocshared=8
%mem=8GB
#p M052X/6-31G(d) Opt Freq

Phenol gas-phase Opt Freq at M05-2X/6-31G(d), one of the SMD parametrization levels

0 1
... (phenol Cartesian coordinates) ...

--Link1--
%chk=phenol_dgsolv.chk
%nprocshared=8
%mem=8GB
#p M052X/6-31G(d) SCRF=(SMD,Solvent=Water) Geom=Check Guess=Read

Phenol SMD(water) single point at the gas-phase optimized geometry

0 1

```

这里选用 **M05-2X/6-31G(d)**，因为它是 SMD 原始参数化所用的方法之一。G16 中 6-31G 系列默认使用 Cartesian d 函数（6D），与参数化时的 6-31G\* 定义一致。

### 5.2 计算

```bash
grep "SCF Done" solvation_free_energy_gas_and_smd.log
```

取第一步（气相优化）**最后一个** `SCF Done` 作为 $E_{\text{gas}}$，第二步（SMD 单点）的 `SCF Done` 作为 $E_{\text{SMD}}$：

$$
\Delta G^{*}_{\text{solv}} \approx E_{\text{SMD}} - E_{\text{gas}}\quad(\times 627.509\ \text{kcal/mol per Hartree})
$$

第一步是 `Opt Freq` 复合作业，会打印两次 `Normal termination`（Opt 和 Freq 各一次），第二步的 SMD 单点再打印一次。下面的命令记录每次 `Normal termination` 之前最后一个 `SCF Done`，取第一个（气相优化结构）和最后一个（SMD 单点）计算 ΔG_solv：

```bash
awk '/SCF Done/{e=$5} /Normal termination/{n++; E[n]=e}
     END{printf "E_gas = %.8f  E_SMD = %.8f  dG_solv = %.2f kcal/mol\n", E[1], E[n], (E[n]-E[1])*627.509}' \
     solvation_free_energy_gas_and_smd.log
```

苯酚在水中的实验溶剂化（水合）自由能约为 −6.6 kcal/mol（Minnesota Solvation Database 等数据集中的常用值），可以用来检验你的结果。

> [!NOTE]
> 严格地说，$E_{\text{SMD}} - E_{\text{gas}}$ 是“固定气相结构”的溶剂化自由能。如果溶质在溶液中的结构变化明显，应分别在两相中优化，ΔG_solv = E(SMD，溶液结构) − E(gas，气相结构)。对苯酚这样的刚性中性分子，两种做法的差别很小。

### 5.3 两种溶液自由能策略的取舍

| 策略 | 计算 | 优点 | 缺点 |
|---|---|---|---|
| **A：溶液中 Opt Freq** | $G^{*}_{\text{sol}} = G_{\text{log}}(\text{SMD}) + 1.89$ | 结构、频率都在溶剂中得到；能描述只在溶液中稳定的物种（两性离子、多价阴离子、离子对） | 溶液中的频率数值噪声稍大；比气相贵 |
| **B：气相 Opt Freq + 溶剂单点** | $G^{*}_{\text{sol}} = G^{\circ}_{\text{gas}} + \Delta G^{*}_{\text{solv}} + 1.89$ | 与 SMD 参数化时的做法一致（多用气相结构）；气相频率更稳健；溶剂单点可以单独换更大的基组 | 气相中不稳定的物种（两性离子等）会给出错误的结构 |
| **A + 高水平单点** | $G^{*}_{\text{sol}} = E_{\text{high}}(\text{SMD}) + [G_{\text{log}} - E](\text{SMD, low}) + 1.89$ | 目前文献中最常见的做法 | 要保证低水平与高水平使用同一溶剂模型 |

**选择建议**：

- 带电体系、强极性 TS、两性离子：用策略 A；
- 中性、刚性分子，或者需要严格对照 SMD 原始定义：用策略 B；
- 一个研究项目中所有物种必须使用**同一种**策略，不能混用。

---

## 6. 自定义溶剂：Solvent=Generic,Read

### 6.1 SMD 自定义溶剂需要 7 个参数

官方文档规定，用 SMD 定义新溶剂时 7 个参数**全部必需**，而且都必须写成实数：

| 关键词 | 物理量 | 说明与来源 |
|---|---|---|
| `Eps` | 静态介电常数 ε | 溶剂手册 |
| `EpsInf` | 光学介电常数 | 等于 293 K 下折射率的平方 $n^2$ |
| `HbondAcidity` | Abraham 氢键酸性 α | Minnesota 溶剂描述符数据库 / Abraham 参数表 |
| `HbondBasicity` | Abraham 氢键碱性 β | 同上 |
| `SurfaceTensionAtInterface` | 宏观表面张力 γ | 单位 cal mol⁻¹ Å⁻²，**1 dyn/cm = 1.4393 cal mol⁻¹ Å⁻²** |
| `CarbonAromaticity` | 芳香碳原子分数 φ | 溶剂分子中芳香碳占重原子的比例 |
| `ElectronegativeHalogenicity` | 卤素原子分数 ψ | F、Cl、Br 占重原子的比例 |

### 6.2 输入：先用内置溶剂验证写法

示例 [`custom_solvent_generic_read.gjf`](../examples/09-solvation-smd/custom_solvent_generic_read.gjf) 做了一个对照实验：第一步用内置的 `Solvent=Ethanol`，第二步用 `Solvent=Generic,Read` 手动输入乙醇的 7 个参数。两步能量应当一致，这样就验证了自定义溶剂的输入格式。第二步为：

```text
--Link1--
%chk=phenol_custom_solvent.chk
%nprocshared=8
%mem=8GB
#p M052X/6-31G(d) SCRF=(SMD,Solvent=Generic,Read) Geom=Check Guess=Read

Phenol in a user-defined SMD solvent (ethanol descriptors typed by hand)

0 1

Eps=24.852
EpsInf=1.852593
HbondAcidity=0.37
HbondBasicity=0.48
SurfaceTensionAtInterface=31.62
CarbonAromaticity=0.0
ElectronegativeHalogenicity=0.0

```

注意空行：电荷多重度行（`Geom=Check` 时分子说明只剩这一行）之后空一行，然后是 PCM 输入段，段末再空一行。

检查：

```bash
grep "SCF Done" custom_solvent_generic_read.log      # 两个能量应当一致
```

例如乙醇的表面张力约为 21.97 dyn/cm，换算后 21.97 × 1.4393 ≈ 31.62 cal mol⁻¹ Å⁻²，正好与上面的取值对应。用同样的方法就可以为其他溶剂准备参数。

### 6.3 只改介电常数（PCM/CPCM）

对于不含 CDS 项的 IEFPCM/CPCM，只需在 PCM 输入段中覆盖 `Eps` 和 `EpsInf`。没有指定的参数沿用 `Solvent=` 所给溶剂（默认是水）的值：

```text
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(IEFPCM,Solvent=Water,Read)

IEFPCM with a user-defined dielectric (eps = 10, n^2 = 2.0); other parameters from water

0 1
... (coordinates) ...

Eps=10.0
EpsInf=2.0

```

### 6.4 其他常用 PCM 输入段关键词

这些关键词写在 `SCRF=(...,Read)` 的附加输入段中。官方文档提醒：大部分参数是给专家准备的，正式计算中不建议随意修改默认值。

| 关键词 | 默认 | 作用 |
|---|---|---|
| `Radii=` | `UFF` | 原子半径集合：`UFF`、`UA0`、`UAHF`、`UAKS`、`Pauling`、`Bondi` |
| `Alpha=` | 1.1 | 静电空腔的半径缩放因子 |
| `Surface=` | `VDW` | 表面类型：`VDW`、`SES`、`SAS` |
| `ModifySph` | — | 修改某些原子球的半径（格式见官方文档） |
| `ExtraSph=N` / `NSph=N` | — | 增加球 / 完全自定义空腔 |
| `PDens=` | 5.0 Å⁻² | 表面积分点密度 |
| `QConv=` | `Tight`（10⁻⁹） | 迭代求解极化电荷的收敛阈值 |
| `Dis` / `Rep` / `Cav` | 不计算 | 在 PCM 中加入色散、排斥、空穴项（Floris–Tomasi / Pierotti 模型），**不能与 SMD 同用** |
| `PrintSpheres` | — | 把空腔信息写入 chk，供 GaussView 显示 |

---

## 7. 激发态：平衡与非平衡溶剂化

溶剂对溶质状态变化的响应分为两部分：**快**的电子极化（光学介电常数 $\varepsilon_\infty = n^2$）和**慢**的取向极化（分子转动）。

- **平衡溶剂化**：两部分都已经适应溶质的当前状态，适用于几何优化；
- **非平衡溶剂化**：只有快的部分跟上了溶质的变化，适用于垂直吸收、垂直发射。

Gaussian 16 的默认行为：

| 计算类型 | 默认溶剂化 |
|---|---|
| CIS/TD-DFT 单点能（垂直激发） | **非平衡**（线性响应） |
| CIS/TD-DFT 激发态几何优化 | 平衡 |
| `SCRF=ExternalIteration`（态特定） | 平衡，需要配合 `NonEquilibrium=Read` 才是非平衡 |

线性响应（LR）方法对电荷转移态的溶剂效应往往描述不足。这时可以用态特定方法，有两种：

1. `SCRF=(Solvent=...,CorrectedLR)`：对 LR 激发能做态特定校正，一步完成；
2. `SCRF=(...,ExternalIteration)`：外迭代，让溶剂反应场与激发态密度自洽。**只能用于能量计算**，不能用于 Opt/Freq。

用外迭代计算非平衡的垂直吸收，需要两步：

```text
%oldchk=mol_gs_opt.chk
%chk=mol_ss.chk
%nprocshared=8
%mem=8GB
#p CAM-B3LYP/def2TZVP SCRF=(Solvent=Acetonitrile,NonEquilibrium=Save) Geom=Check Guess=Read

Step 1: ground state at its optimized geometry, save the slow solvent polarization

0 1

--Link1--
%chk=mol_ss.chk
%nprocshared=8
%mem=8GB
#p CAM-B3LYP/def2TZVP TD=(NStates=6,Root=1) SCRF=(Solvent=Acetonitrile,ExternalIteration,NonEquilibrium=Read) Geom=Check Guess=Read

Step 2: state-specific (external iteration) energy of excited state 1, non-equilibrium

0 1

```

第二步日志末尾的 `After PCM corrections, the energy is ...` 减去第一步的基态能量，就是态特定非平衡溶剂化下的垂直激发能。荧光发射的完整六步流程（激发态优化、平衡溶剂化保存、基态非平衡读取）见 [第10章](10-tddft.md)。

---

## 8. 热力学循环一：pKa（苯酚）

### 8.1 直接法

$$
\mathrm{HA(aq)} \longrightarrow \mathrm{A^-(aq)} + \mathrm{H^+(aq)}
$$

$$
\Delta G^{*}_{\text{aq}} = G^{*}_{\text{aq}}(\mathrm{A^-}) + G^{*}_{\text{aq}}(\mathrm{H^+}) - G^{*}_{\text{aq}}(\mathrm{HA}),\qquad
\mathrm{p}K_a = \frac{\Delta G^{*}_{\text{aq}}}{RT\ln 10} = \frac{\Delta G^{*}_{\text{aq}}}{1.364\ \text{kcal/mol}}
$$

- $G^{*}_{\text{aq}}(\mathrm{HA})$、$G^{*}_{\text{aq}}(\mathrm{A^-})$：取自 SMD(水) 中的 Opt Freq（策略 A）。HA 和 A⁻ 各自的 +1.89 kcal/mol 相互抵消，可以都不加；
- $G^{*}_{\text{aq}}(\mathrm{H^+})$ 无法用连续介质模型计算，取文献值：

$$
G^{*}_{\text{aq}}(\mathrm{H^+}) = G^{\circ}_{\text{gas}}(\mathrm{H^+}) + \Delta G^{\circ\rightarrow *} + \Delta G^{*}_{\text{solv}}(\mathrm{H^+}) = -6.28 + 1.89 - 265.9 = -270.3\ \text{kcal/mol}
$$

其中 −6.28 kcal/mol 是 298.15 K、1 atm 下气相质子的自由能（Sackur–Tetrode 平动熵），−265.9 kcal/mol 是 Kelly、Cramer、Truhlar 推荐的质子绝对水合自由能（不确定度不小于 2 kcal/mol）。

示例文件 `phenol_smd_water_opt_freq.gjf` 和 `phenoxide_smd_water_opt_freq.gjf` 分别计算 HA 和 A⁻。取出自由能：

```bash
grep "Sum of electronic and thermal Free Energies" phenol_smd_water_opt_freq.log phenoxide_smd_water_opt_freq.log
```

然后用下面这段小程序计算（把两个自由能作为命令行参数传入，单位 Hartree）：

```python
# pka_direct.py  -- usage: python pka_direct.py G_HA G_A   (Hartree, SMD Opt Freq)
import sys

HARTREE_TO_KCAL = 627.509
G_PROTON_AQ = -6.28 + 1.89 - 265.9   # kcal/mol, 1 M standard state, 298.15 K
RT_LN10 = 1.3642                     # kcal/mol at 298.15 K

g_ha, g_a = float(sys.argv[1]), float(sys.argv[2])
dg = (g_a - g_ha) * HARTREE_TO_KCAL + G_PROTON_AQ
print(f"dG_aq = {dg:.2f} kcal/mol   pKa(direct) = {dg / RT_LN10:.1f}")
```

苯酚的实验 pKa 约为 10.0。

> [!WARNING]
> pKa 对误差极其敏感：**1.36 kcal/mol 的自由能误差就对应 1 个 pKa 单位**。SMD 对离子溶剂化自由能的平均误差约 4 kcal/mol，所以直接法的误差常达到 2–5 个 pKa 单位。不要把直接法的绝对值当作定量预测。

### 8.2 更可靠的做法

1. **质子交换（相对）法**：选一个结构相近、pKa 已知的参比酸 HRef：

   $$
   \mathrm{HA + Ref^- \longrightarrow A^- + HRef},\qquad
   \mathrm{p}K_a(\mathrm{HA}) = \mathrm{p}K_a(\mathrm{HRef}) + \frac{\Delta G^{*}_{\text{exch}}}{1.364}
   $$

   质子不出现在方程中，阴离子溶剂化误差大部分可以抵消。例如计算取代苯酚时，以苯酚（pKa ≈ 10.0）作为参比。
2. **簇–连续介质模型**：给阴离子加 1–3 个显式水分子，再放入 SMD 中，改善局部氢键的描述；
3. **线性校正**：对一组同类化合物，把计算的 ΔG 与实验 pKa 做线性拟合，然后用拟合关系预测新化合物。

方法综述见 Ho & Coote 2010（参考资料）。

---

## 9. 热力学循环二：氧化还原电位（TTF，乙腈）

### 9.1 输入

四硫富瓦烯（TTF）是经典的有机电子给体，其氧化 TTF → TTF•⁺ + e⁻ 是可逆的单电子过程。示例 [`redox_potential_cycle_neutral_and_cation.gjf`](../examples/09-solvation-smd/redox_potential_cycle_neutral_and_cation.gjf) 用 4 个 Link1 步骤完成计算：

| 步骤 | 物种 | 计算 | 作用 |
|---|---|---|---|
| 1 | TTF（0 1） | B3LYP-D3(BJ)/def2-SVP，SMD(乙腈)，Opt Freq | 结构 + 热校正 |
| 2 | TTF•⁺（1 2） | UB3LYP-D3(BJ)/def2-SVP，SMD(乙腈)，Opt Freq，从中性结构出发 | 结构 + 热校正 |
| 3 | TTF | M06-2X/def2-TZVP，SMD(乙腈)，单点 | 高水平电子能 |
| 4 | TTF•⁺ | UM06-2X/def2-TZVP，SMD(乙腈)，单点 | 高水平电子能 |

第 2 步的关键行：

```text
--Link1--
%oldchk=ttf_neutral.chk
%chk=ttf_cation.chk
%nprocshared=8
%mem=8GB
#p UB3LYP/def2SVP EmpiricalDispersion=GD3BJ SCRF=(SMD,Solvent=Acetonitrile) Opt Freq Geom=Check

Step 2: TTF radical cation in acetonitrile (SMD), Opt Freq starting from the neutral geometry

1 2

```

### 9.2 计算公式

$$
G_{\text{sol}}(X) = E_{\text{high}}(X) + \left[G_{\text{log}}(X) - E_{\text{low}}(X)\right]
$$

其中方括号内是低水平的热校正（即 `Thermal correction to Gibbs Free Energy` 行）。氧化自由能：

$$
\Delta G_{\text{ox}} = G_{\text{sol}}(\mathrm{TTF^{\bullet+}}) - G_{\text{sol}}(\mathrm{TTF})
$$

电子不在体系内显式出现，它的贡献包含在参比电极的绝对电位中。氧化过程两边溶质个数相同，所以 1 M 标准态校正相互抵消。相对某参比电极的标准电位（n = 1，ΔG 以 eV 为单位时，数值上就等于伏特）：

$$
E^\circ(\mathrm{TTF^{\bullet+}/TTF}) = \frac{\Delta G_{\text{ox}}}{nF} - E_{\text{abs}}(\text{ref})
$$

| 参比 | 绝对电位或换算 | 来源与说明 |
|---|---|---|
| SHE（水溶液） | 4.28 V | Isse & Gennaro 2010，与 ΔG_solv(H⁺) = −265.9 kcal/mol 自洽 |
| SHE（IUPAC 推荐） | 4.44 V | Trasatti 1986。两个数值相差 0.16 V，文献中都有人用，**必须注明** |
| Fc⁺/Fc（乙腈） | 约 +0.38 V vs SCE | Pavlishchuk & Addison 2000，乙腈中实测换算常数 |

单位换算：1 Hartree = 27.2114 eV，1 eV = 23.0605 kcal/mol = 96.485 kJ/mol。

> [!TIP]
> 有机溶剂中最稳妥的做法是**在同一理论水平下同时计算二茂铁**，报告相对 Fc⁺/Fc 的电位：
> $E(\text{vs Fc}^+/\text{Fc}) = [\Delta G_{\text{ox}}(\text{M}) - \Delta G_{\text{ox}}(\text{Fc})]/F$。
> 这样可以避开绝对电位的数值选择、液接电势等问题，系统误差也会大部分抵消。实验电化学数据通常也以 Fc⁺/Fc 为内标。

### 9.3 注意事项

- 自由基阳离子要检查自旋污染：`grep "S\*\*2" ...log`，`<S**2>` 应接近 0.75（[第13章](13-openshell-stability.md)）；
- 多电子过程（n = 2）或伴随质子转移的过程（PCET），需要构建更完整的热力学方格（square scheme）；
- 也可以用 Born–Haber 循环：$\Delta G_{\text{ox}}(\text{sol}) = \Delta G_{\text{ox}}(\text{gas}) + \Delta G_{\text{solv}}(\mathrm{M^+}) - \Delta G_{\text{solv}}(\mathrm{M})$，即“气相 Opt Freq + SMD 单点”（策略 B）。两种循环不要混用。

---

## 10. 与其他计算类型组合

| 组合 | 写法要点 | 详见 |
|---|---|---|
| 溶液中的 TS 与 IRC | Opt=TS、Freq、IRC、端点优化都用同一个 `SCRF=(...)` | [第7章](07-ts-irc.md) |
| 溶液中的扫描 | `Opt=ModRedundant SCRF=(SMD,Solvent=...)` | [第8章](08-scan.md) |
| UV-Vis / 荧光 | `TD=(NStates=...) SCRF=(...)`，默认非平衡；态特定用 `ExternalIteration` 或 `CorrectedLR` | [第10章](10-tddft.md) |
| NMR | `NMR SCRF=(SMD,Solvent=Chloroform)`；参比 TMS 也必须在同一溶剂中计算 | [第11章](11-nmr-props.md) |
| ESP 电荷 | RESP 电荷通常在气相 HF/6-31G(d) 下拟合；溶剂中的电荷要单独说明 | [第12章](12-population-charges.md) |
| ONIOM | `ONIOM(...) SCRF=(...)`，默认方案为 `ONIOMPCM=X`，另有 A/B/C 三种 | [第15章](15-oniom.md) |

---

## 参数选择与收敛建议

| 项目 | 推荐 | 理由 |
|---|---|---|
| 计算 ΔG_solv、pKa、电位 | `SCRF=(SMD,Solvent=...)` | IEFPCM 只有静电项 |
| 溶剂对结构/光谱的影响 | IEFPCM（默认）或 SMD 均可 | 静电项占主导 |
| 溶剂化能的理论水平 | 与 SMD 参数化一致更稳妥：M05-2X/6-31G(d) 或 B3LYP/6-31G(d)；常用的 M06-2X/6-311+G(d,p) 也可以 | SMD 参数在多种方法下拟合，换方法后误差会有一些变化 |
| 阴离子 | 基组中加弥散函数（6-311+G(d,p)、def2-TZVPD、aug-cc-pVTZ 等） | 外层电子分布弥散 |
| 溶液中的优化 | 默认判据；可疑小虚频时 `Opt=Tight`、`Int=SuperFine` | 空腔离散化带来数值噪声 |
| 标准态 | 每个溶质 +1.89 kcal/mol；作为反应物的溶剂再 +2.38 kcal/mol | 1 atm → 1 M；纯液体 55.34 M |
| 方案一致性 | 同一项目中所有物种使用相同的模型、半径、版本和策略（A 或 B） | 避免人为误差 |
| 版本 | 复现 G09 文献时用 `SCRF=G09Defaults` | G09 默认非对称 IEFPCM |
| 激发态 | 垂直激发默认非平衡 LR；CT 态考虑 `ExternalIteration` 或 `CorrectedLR` | LR 对 CT 态的溶剂效应描述不足 |

---

## 常见问题与报错

### 1. `Inv3 failed in PCMMkU`（L502）

**含义**：构建 PCM 方程时矩阵求逆失败，问题出在空腔或表面离散化上。

**可尝试的处理**（按顺序）：

1. 检查几何结构：有没有原子重叠、键长单位错误、断裂成碎片。不合理的结构会产生病态空腔；
2. 先在气相中预优化，再在溶剂中继续优化；
3. 稍微改变空腔，例如用 `SCRF=(...,Read)` 加 `Alpha=1.2`，或换一种半径集合（只限 IEFPCM/CPCM；SMD 的半径是模型参数的一部分，不应随意修改）；
4. 改用 `SCRF=CPCM` 做对照；
5. 旧版本（G09）中出现时，换到 G16（连续表面电荷形式更稳定）。

### 2. `GetVDW: no radius for atom`

**含义**：所用的半径集合中没有该元素的半径，常见于过渡金属、镧系元素，以及 `Radii=Bondi`、`UAHF` 等不覆盖全部元素的半径集合。

**处理**：

- IEFPCM/CPCM 中改回默认的 `Radii=UFF`（UFF 覆盖全部元素）；
- 用 PCM 输入段的 `ModifySph` 为该元素手动指定半径（格式见官方 SCRF 页面），并在论文中说明取值和依据；
- SMD 计算中如果遇到这种情况，同样只为缺失的元素补半径，其余元素保持 SMD 默认值。

### 3. 带 `ExternalIteration` 的优化或频率作业报错

`ExternalIteration` 只能用于能量计算。激发态优化用 LR 平衡溶剂化（`TD=(Root=N) SCRF=(...) Opt`），得到结构后再做外迭代单点。

### 4. 离子的“溶剂化能”有几十甚至上百 kcal/mol，是不是算错了？

不是。一价离子在水中的溶剂化自由能通常是 −60 到 −100 kcal/mol 量级（Born 公式中 $q^2$ 项的结果）。真正需要警惕的是：不同电荷物种之间相减时，误差不会像中性分子那样抵消，所以涉及电荷变化的反应误差更大。

### 5. 同一体系，G09 与 G16 的溶剂能量不同

默认 IEFPCM 形式（非对称 → 对称）、积分格点（FineGrid → UltraFine）、物理常数都变了。复现旧结果时用 `SCRF=G09Defaults`，或者干脆在 route 中加 `G09Defaults`。

### 6. 忘了标准态校正

溶液中由 Gaussian 输出得到的自由能是 1 atm 标准态。计算反应自由能时，如果反应前后分子数不同，必须对每个溶质加 1.89 kcal/mol（298.15 K），否则双分子反应的 ΔG 和 ΔG‡ 会有系统偏差。

### 7. 溶剂名报语法错误（QPErr）

溶剂名必须是官方列表中的关键词，例如 `DiMethylSulfoxide`、`n,n-DiMethylFormamide`、`TetraHydroFuran`。写成 `DMSO`、`DMF`、`THF` 等缩写之前，请先在 [SCRF 页面](https://gaussian.com/scrf/) 的溶剂表中确认它是否被接受。

### 8. 溶液中优化振荡、不收敛

溶剂中的势能面通常比气相更平。可以先在气相中优化，再以该结构为起点在溶剂中优化；或者用 `Opt=CalcFC`、减小 `MaxStep`。诊断方法见 [第5章](05-opt.md) 和仓库脚本 `opt_monitor.py`（`python ../scripts/opt_monitor.py -h`）。也可以用 [`g16_error_doctor.py`](../scripts/g16_error_doctor.py) 自动识别报错类型（`python ../scripts/g16_error_doctor.py -h` 查看用法）。

---

## 练习

1. 用 `solvation_free_energy_gas_and_smd.gjf` 的写法，分别计算苯酚在水、正辛醇（`n-Octanol`）中的 ΔG_solv，由两者之差估算 log P（$\log P = -[\Delta G_{\text{solv}}(\text{oct}) - \Delta G_{\text{solv}}(\text{water})]/(RT\ln 10)$），并与实验值比较。
2. 用 IEFPCM 和 SMD 分别计算苯酚在水中的 `E(SCRF) − E(gas)`，解释二者为什么不同。
3. 用直接法和质子交换法（参比酸选乙酸或苯酚）分别计算对硝基苯酚的 pKa，比较两种方法的误差。
4. 在同一水平下计算二茂铁的氧化自由能，报告 TTF 相对 Fc⁺/Fc 的电位。
5. 把 `custom_solvent_generic_read.gjf` 第二步中的 `HbondBasicity` 改为 0.0，看 ΔG_solv 变化多少，体会各个 SMD 描述符的作用。

---

## 参考资料

**官方文档**

- [SCRF 关键词（PCM/IEFPCM/CPCM/SMD、Solvent 列表、Generic,Read、ExternalIteration、NonEquilibrium、PCM 输入段关键词）](https://gaussian.com/scrf/)
- [TD 关键词（激发态与溶剂）](https://gaussian.com/td/)
- [Thermochemistry in Gaussian（白皮书，标准态与热校正）](https://gaussian.com/thermo/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes（G09 → G16 默认值变化：SCRF 默认对称 IEFPCM 等）](https://gaussian.com/relnotes/)
- [Basis Sets（6-31G 系列的元素范围与弥散函数）](https://gaussian.com/basissets/)

**相关 GitHub 仓库**（只链接，不复制代码）

- [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)：qRRHO、1 M 标准态校正、单点能与热校正组合
- [lenhanpham/OpenThermo](https://github.com/lenhanpham/OpenThermo)：热化学计算，支持 qRRHO 和温度/压力扫描
- [kimariyb/easy-Shermo](https://github.com/kimariyb/easy-Shermo)：批量调用 Shermo，自动配对 opt/freq 与高精度单点（中文 README）
- [jvalegre/aqme](https://github.com/jvalegre/aqme)：从 SMILES 到 Gaussian 输入的自动化流程，可统一设置溶剂
- [cclib/cclib](https://github.com/cclib/cclib)：解析日志中的能量与热化学数据
- 更多工具见 [../../Resources/github-repos.md](../../Resources/github-repos.md)

**本仓库脚本**：[`../scripts/thermo_qrrho.py`](../scripts/thermo_qrrho.py)（qRRHO、1 M 校正、高水平单点组合）、[`../scripts/g16_log_summary.py`](../scripts/g16_log_summary.py)（批量汇总能量与热化学）、[`../scripts/gjf_builder.py`](../scripts/gjf_builder.py)（生成带溶剂的输入）。各脚本用 `-h` 查看用法。

**关键文献**

- A. V. Marenich, C. J. Cramer, D. G. Truhlar, Universal solvation model based on solute electron density and on a continuum model of the solvent defined by the bulk dielectric constant and atomic surface tensions, *J. Phys. Chem. B* **2009**, 113, 6378–6396. DOI: [10.1021/jp810292n](https://doi.org/10.1021/jp810292n)
- J. Tomasi, B. Mennucci, R. Cammi, Quantum mechanical continuum solvation models, *Chem. Rev.* **2005**, 105, 2999–3094. DOI: [10.1021/cr9904009](https://doi.org/10.1021/cr9904009)
- G. Scalmani, M. J. Frisch, Continuous surface charge polarizable continuum models of solvation. I. General formalism, *J. Chem. Phys.* **2010**, 132, 114110. DOI: [10.1063/1.3359469](https://doi.org/10.1063/1.3359469)
- C. P. Kelly, C. J. Cramer, D. G. Truhlar, Aqueous solvation free energies of ions and ion–water clusters based on an accurate value for the absolute aqueous solvation free energy of the proton, *J. Phys. Chem. B* **2006**, 110, 16066–16081. DOI: [10.1021/jp063552y](https://doi.org/10.1021/jp063552y)
- J. Ho, M. L. Coote, A universal approach for continuum solvent pKa calculations: are we there yet?, *Theor. Chem. Acc.* **2010**, 125, 3–21. DOI: [10.1007/s00214-009-0667-0](https://doi.org/10.1007/s00214-009-0667-0)
- A. A. Isse, A. Gennaro, Absolute potential of the standard hydrogen electrode and the problem of interconversion of potentials in different solvents, *J. Phys. Chem. B* **2010**, 114, 7894–7899. DOI: [10.1021/jp100402x](https://doi.org/10.1021/jp100402x)
- V. V. Pavlishchuk, A. W. Addison, Conversion constants for redox potentials measured versus different reference electrodes in acetonitrile solutions at 25 °C, *Inorg. Chim. Acta* **2000**, 298, 97–102. DOI: [10.1016/S0020-1693(99)00407-7](https://doi.org/10.1016/S0020-1693%2899%2900407-7)

---

[← 上一章](08-scan.md) | [返回目录](../README.md) | [下一章 →](10-tddft.md)
