# 第6章 频率计算、红外拉曼与热化学（Frequencies, IR/Raman & Thermochemistry in Gaussian 16）

本章讲解 Gaussian 16 频率计算（frequency analysis）的正确用法：如何用虚频判断极小点和过渡态，用 `Freq=Raman` 得到 IR 和 Raman 光谱，用 `Freq=Anharmonic` 做非谐（VPT2）频率，以及用 `Freq=ReadIsotopes`/`Iso=` 做同位素替换。热化学部分讲解零点能（ZPE）、焓 H 和吉布斯自由能 G 的输出含义，频率缩放因子，Grimme qRRHO 准谐振校正（GoodVibes/Shermo），1 atm → 1 mol/L 标准态校正，以及“高精度单点 + 低水平热校正”的自由能组合公式。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 写出正确的 `Opt Freq` 作业，理解为什么频率必须与优化使用完全相同的理论水平；
- 根据虚频个数判断驻点类型，并处理优化后残留的小虚频；
- 得到 IR/Raman 强度、非谐频率和同位素位移，画出可与实验对比的光谱；
- 读懂热化学输出中的每一行，按需改变温度、压力和缩放因子；
- 对低频模式做 qRRHO 校正，对溶液反应做标准态校正；
- 正确组合高水平单点能与低水平热校正，得到可发表的自由能。

| 用途 | 本章相关小节 |
|---|---|
| 确认极小点或过渡态 | 第 4 节 |
| 模拟 IR/Raman 光谱，指认实验谱峰 | 第 5、6 节 |
| 反应焓、反应自由能、能垒 ΔG‡ | 第 7–11 节 |
| 同位素效应、动力学同位素效应的前置计算 | 第 6 节 |

> [!NOTE]
> 全部示例输入在 [../examples/06-freq-thermo/README.md](../examples/06-freq-thermo/README.md)。几何优化见[第5章](05-opt.md)，过渡态与 IRC 见[第7章](07-ts-irc.md)，溶剂模型见[第9章](09-solvation.md)。

---

## 2. 原理速览

### 2.1 谐振近似与简正模式

在驻点附近把能量展开到二阶，Gaussian 计算笛卡尔坐标下的 Hessian（力常数矩阵）$\mathbf{H}$，转换成质量加权坐标后对角化：

$$
\mathbf{H}^{\mathrm{mw}}_{ij} = \frac{H_{ij}}{\sqrt{m_i m_j}},\qquad
\mathbf{H}^{\mathrm{mw}}\,\mathbf{L} = \mathbf{L}\,\boldsymbol{\Lambda},\qquad
\tilde{\nu}_k = \frac{\sqrt{\lambda_k}}{2\pi c}
$$

投影掉 3 个平动和 3 个转动（线性分子是 2 个）后，剩下 3N−6（线性分子 3N−5）个振动模式。本征值 $\lambda_k < 0$ 时频率是虚数，Gaussian 以**负数**打印，称为虚频。

> [!IMPORTANT]
> 官方 Freq 文档强调：笛卡尔 Hessian 到质量加权简正坐标的变换**只在驻点上成立**。在 6-31G(d) 优化的结构上算 6-311G(d) 频率，结果没有意义。频率必须与优化使用相同的方法、基组、溶剂模型、积分格点（以及冻芯设置）。最简单的保证方式是 `Opt Freq` 写在同一个 route 里。

### 2.2 强度

- **IR 强度**正比于偶极矩对简正坐标导数的平方 $|\partial\boldsymbol{\mu}/\partial Q_k|^2$，单位 km/mol；Freq 默认计算。
- **Raman 活性**由极化率导数给出：$S_k = 45\,\alpha_k'^2 + 7\,\gamma_k'^2$（单位 Å⁴/amu）。DFT 默认**不算**，需要 `Freq=Raman`；HF 默认计算。

### 2.3 热化学（理想气体 + 刚性转子 + 谐振子，RRHO）

Gaussian 按理想气体正则系综计算，配分函数分解为 $q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}$：

$$
\mathrm{ZPE} = \frac{1}{2}\sum_k h c\,\tilde{\nu}_k
$$

$$
S_{\mathrm{vib}} = R\sum_k\left[\frac{x_k}{e^{x_k}-1} - \ln\left(1-e^{-x_k}\right)\right],\qquad x_k = \frac{h c\,\tilde{\nu}_k}{k_{\mathrm B}T}
$$

$$
S_{\mathrm{trans}} = R\left[\ln\left(\left(\frac{2\pi m k_{\mathrm B}T}{h^2}\right)^{3/2}\frac{k_{\mathrm B}T}{P}\right) + \frac{5}{2}\right]
$$

$$
H_{\mathrm{corr}} = E_{\mathrm{corr}} + k_{\mathrm B}T,\qquad G_{\mathrm{corr}} = H_{\mathrm{corr}} - T\,S_{\mathrm{tot}}
$$

几个要点：

- 压力只通过 $S_{\mathrm{trans}}$ 影响结果，这是后面做标准态校正的依据；
- 转动熵中含有转动对称数 σ（水为 2，NH3 为 3，CH4 为 12），由点群自动确定；
- $\tilde{\nu}\to 0$ 时 $S_{\mathrm{vib}}$ 发散，所以低频模式会让 RRHO 熵严重偏大，这是需要 qRRHO 校正的原因；
- 虚频不计入热化学，所以对过渡态可以直接使用 Gaussian 给出的热校正。

完整公式见 Gaussian 白皮书 [Thermochemistry in Gaussian](https://gaussian.com/thermo/) 和 [Vibrational Analysis in Gaussian](https://gaussian.com/vib/)。

---

## 3. 第一个频率计算：H2O 的 Opt Freq=Raman

### 3.1 输入文件

[h2o_opt_freq_raman.gjf](../examples/06-freq-thermo/h2o_opt_freq_raman.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=h2o_opt_freq_raman.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt Freq=Raman

H2O opt + harmonic frequencies with IR and Raman intensities

0 1
O     0.000000     0.000000     0.117300
H     0.000000     0.757200    -0.469200
H     0.000000    -0.757200    -0.469200

```

| 关键词 | 说明 |
|---|---|
| `Opt Freq=Raman` | 先优化，再在同一水平下算频率，并额外计算 Raman 活性 |
| `#p` | 详细输出；`#p` 下频率作业还会打印极化率张量，小分子还会打印简单的光谱示意图 |
| `%chk` | 保存力常数和简正模式，后面的 ReadFC、ReadIsotopes 示例都要读它 |

### 3.2 运行与快速检查

```bash
g16 < h2o_opt_freq_raman.gjf > h2o_opt_freq_raman.log
grep -c "Normal termination" h2o_opt_freq_raman.log          # expect 2
grep "Frequencies --" h2o_opt_freq_raman.log                  # all positive -> minimum
grep "imaginary frequencies" h2o_opt_freq_raman.log           # no output -> 0 imaginary modes
grep -A1 "Low frequencies" h2o_opt_freq_raman.log
```

### 3.3 频率输出结构

（示意输出，数值仅供参考）

```text
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:
                      1                      2                      3
                     A1                     A1                     B2
 Frequencies --   16xx.xxxx              38xx.xxxx              39xx.xxxx
 Red. masses --      1.08xx                 1.04xx                 1.08xx
 Frc consts --       1.7xxx                 8.8xxx                 9.6xxx
 IR Inten    --     7x.xxxx                 x.xxxx                5x.xxxx
 Raman Activ --      x.xxxx               10x.xxxx                3x.xxxx
 Depolar (P) --      0.xxxx                 0.xxxx                 0.7500
 Depolar (U) --      0.xxxx                 0.xxxx                 0.8571
  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
     1   8     0.00   0.00   0.07     0.00   0.00  -0.05     0.00   0.07   0.00
   ...
```

| 行 | 含义 |
|---|---|
| `Frequencies --` | 谐振频率（cm⁻¹），负数为虚频 |
| `Red. masses --` | 约化质量（amu） |
| `Frc consts --` | 力常数（mDyne/Å） |
| `IR Inten --` | IR 强度（km/mol） |
| `Raman Activ --` | Raman 活性（Å⁴/amu），只有 `Freq=Raman` 时出现 |
| `Depolar (P)/(U)` | 平面偏振、非偏振入射光下的退偏比；非全对称模式的理论值分别为 0.75 和 0.857 |
| 下方的 X Y Z | 简正模式中各原子的位移（归一化），三位有效数字；`Freq=HPModes` 额外输出五位有效数字的版本 |

作为参照，气相水的实验基频为 1595（弯曲）、3657（对称伸缩）、3756（反对称伸缩）cm⁻¹（NIST CCCBDB 实验振动数据）。谐振计算值通常比实验基频高约 3–5%，原因是非谐性和方法误差，这就是要做频率缩放（第 8 节）或非谐计算（第 5.3 节）的原因。

### 3.4 "Low frequencies" 两行

```text
 Low frequencies ---   -x.xxxx   -x.xxxx   -0.00xx    0.00xx    x.xxxx   xx.xxxx
 Low frequencies ---  16xx.xxxx 38xx.xxxx 39xx.xxxx
```

（示意输出，数值仅供参考）

第一行来自 `DiagFull` 分析（官方文档：对角化包含平动和转动的完整 Hessian，用来检验数值稳定性），前 6 个值对应平动和转动，理想情况下应接近 0。如果这几个值达到几十 cm⁻¹，说明结构没有充分收敛，或者格点数值噪声较大，此时低频振动模式也不可靠，应该收紧优化（`Opt=Tight`）、加密格点（`Int=SuperFine`）后重算。第二行是最低的几个真实振动频率。

---

## 4. 虚频判断与处理

### 4.1 判断规则

| 虚频个数 | 驻点类型 | 下一步 |
|---|---|---|
| 0 | 极小点（反应物、产物、中间体、构象） | 可以做热化学 |
| 1 | 一阶鞍点（过渡态候选） | 用 GaussView 动画确认虚频模式对应反应坐标，再做 IRC（[第7章](07-ts-irc.md)） |
| ≥ 2 | 高阶鞍点 | 不是想要的结构，沿虚频方向位移后重新优化 |

有虚频时日志中会出现（示意输出）：

```text
 ******    1 imaginary frequencies (negative Signs) ******
```

```bash
grep "imaginary frequencies" job.log      # number of imaginary modes
grep -m 1 "Frequencies --" job.log        # the first value is the lowest mode
```

### 4.2 极小点上残留的小虚频

优化“收敛”了，却出现一个绝对值很小的虚频（例如 −15 cm⁻¹，通常是甲基转动或分子间的软模式），常见原因和对策：

| 原因 | 检查方法 | 对策 |
|---|---|---|
| 优化收敛得不够紧 | Freq 步末尾的收敛表出现 `NO`（[第5章](05-opt.md) 3.4 节） | `Geom=Check Guess=Read Opt=(ReadFC,Tight) Freq` 重新优化 |
| 积分格点噪声 | "Low frequencies" 第一行的值偏大 | 优化和频率都用 `Int=SuperFine` |
| 对称性把结构限制在鞍点上 | 虚频模式破坏分子对称性 | 沿虚频位移后加 `NoSymm` 重新优化 |
| 真正的鞍点 | 虚频较大（几十到几百 i cm⁻¹）、模式明确 | 沿虚频方向位移，重新优化到真正的极小点 |

**沿虚频方向位移**的做法：

- GaussView：打开 log，`Results → Vibrations`，选中虚频模式，用 `Manual Displacement` 把结构沿该模式推一小段距离，然后另存为新的输入文件；
- [pyQRC](https://github.com/patonlab/pyQRC)：读取 Gaussian 频率输出，自动沿虚频模式位移并写出新输入（用法见其 README）；
- 位移后用 `Opt=CalcFC` 重新优化，再做频率验证。

> [!WARNING]
> 不要因为虚频“很小”就直接忽略。小虚频说明结构不是真正的驻点，热化学（特别是熵）会因此不可靠。如果反复处理后仍残留 < 10i cm⁻¹ 的虚频，并且确认来自格点噪声，至少要在论文中说明。

---

## 5. IR 与 Raman 光谱

### 5.1 从频率到光谱

Gaussian 给出的是离散的谱线（频率加强度），画光谱需要两步：

1. **缩放**频率：乘以对应方法基组的**基频**缩放因子（第 8 节）；
2. **展宽**：每条线用洛伦兹函数（IR/Raman 常用，FWHM 约 8–20 cm⁻¹）或高斯函数展宽后叠加。

可用的工具：

- GaussView：`Results → Vibrations → Spectrum`，可以设置峰宽和缩放；
- Multiwfn 主功能 11（[第16章](16-postprocess-multiwfn.md)）；
- 本仓库的 [spectrum_ir_raman.py](../scripts/spectrum_ir_raman.py)：解析频率、IR 强度和 Raman 活性，缩放、洛伦兹展宽并作图（`python ../scripts/spectrum_ir_raman.py -h` 查看用法）。

### 5.2 Raman 活性与 Raman 强度

Gaussian 输出的是 Raman **活性** $S_k$，与入射激光无关。实验测到的 Raman **强度**还取决于激光波数 $\tilde{\nu}_0$ 和温度，常用换算式为：

$$
I_k \propto \frac{(\tilde{\nu}_0 - \tilde{\nu}_k)^4\,S_k}{\tilde{\nu}_k\left[1 - \exp\left(-h c\,\tilde{\nu}_k / k_{\mathrm B}T\right)\right]}
$$

与实验 Raman 谱比较时，应先换算成强度（spectrum_ir_raman.py 提供此选项）。另外，预共振 Raman 可以用 `Freq=Raman` 配合 `CPHF=RdFreq` 并在输入中给出入射光频率来计算（见官方 Freq 文档）。

### 5.3 非谐频率：Freq=Anharmonic

`Freq=Anharmonic` 用二阶振动微扰理论（VPT2）计算非谐基频、倍频、组合频和非谐 ZPE。它通过沿各简正模式数值微分解析 Hessian 来得到三阶、四阶导数（官方文档：只适用于有解析二阶导数的方法，包括 HF、DFT、CIS 和 MP2），每个简正模式正负各位移一次，需要约 2×(3N−6) 次额外的 Hessian 计算，计算量随分子增大迅速上升，只适合小分子。

[formaldehyde_freq_anharmonic.gjf](../examples/06-freq-thermo/formaldehyde_freq_anharmonic.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=formaldehyde_freq_anharmonic.chk
#p B3LYP/def2TZVP Opt=VeryTight Int=SuperFine Freq=Anharmonic

Formaldehyde: VeryTight opt + anharmonic (VPT2) frequencies

0 1
C     0.000000     0.000000     0.000000
O     0.000000     0.000000     1.205000
H     0.000000     0.934900    -0.579700
H     0.000000    -0.934900    -0.579700

```

| 关键词 | 为什么 |
|---|---|
| `Opt=VeryTight` | 数值三阶、四阶导数对结构偏离驻点非常敏感，需要很紧的收敛 |
| `Int=SuperFine` | 降低格点噪声对数值微分的影响 |
| `Freq=Anharmonic` | VPT2 非谐分析；数值微分步长默认 0.025 Å（官方文档） |

结果解读（字段名与官方文档示例一致）：

```bash
grep -A 6 "Anharmonic Zero Point Energy" formaldehyde_freq_anharmonic.log
grep -A 12 "Fundamental Bands" formaldehyde_freq_anharmonic.log
```

- `Anharmonic Zero Point Energy` 块：谐振 ZPE、非谐修正和总非谐 ZPE；
- `Fundamental Bands` 表：`E(harm)` 为谐振频率，`E(anharm)` 为非谐基频，`I(anharm)` 为非谐 IR 强度；
- 之后是 `Overtones`（倍频）和 `Combination Bands`（组合频）。

甲醛气相实验基频约为 2782、1746、1500、1167、2843、1249 cm⁻¹（NIST CCCBDB 实验振动数据），可以与 `E(anharm)` 直接比较，而不需要缩放因子。

> [!NOTE]
> 存在 Fermi 共振（例如 C–H 伸缩与弯曲模式倍频接近）时，VPT2 对相关谱带的处理比较敏感，与实验比较时要重点核对这些谱带。高级用法（`ReadAnharm` 输入节、选模 `SelectAnharmonicModes`）见官方 Freq 文档。

---

## 6. 同位素替换：ReadIsotopes 与 Iso=

在 Born-Oppenheimer 近似下，同位素替换**不改变势能面**，所以结构和力常数都不变，只有质量变了。因此不需要重新优化，只需要用新质量重做简正分析。

### 6.1 方法一：Freq=(ReadFC,ReadIsotopes) 从 chk 重新分析

[d2o_freq_readisotopes.gjf](../examples/06-freq-thermo/d2o_freq_readisotopes.gjf)（先运行 3.1 节的 H2O 作业生成 chk）：

```text
%nprocshared=8
%mem=12GB
%oldchk=h2o_opt_freq_raman.chk
%chk=d2o_freq_readisotopes.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Freq=(ReadFC,ReadIsotopes) Geom=AllCheck Guess=Read

298.15 1.0
16
2
2

```

| 要点 | 说明 |
|---|---|
| `ReadFC` | 从 chk 读取已经算好的力常数，只重做简正和热化学分析，几乎不花时间。官方文档提醒：用 ReadFC 时基组从 chk 读取，不要再输入 Gen 基组；上一作业用过 Raman 的话，这里不要再写 Raman |
| `ReadIsotopes` | 读入一节“温度 压力 [缩放因子]”，后面每行一个原子的同位素质量，顺序与分子坐标相同 |
| `298.15 1.0` | 温度（K）和压力（atm），可选的第三个数是缩放因子 |
| `16` `2` `2` | 整数质量会自动换成对应同位素的精确质量（16O、D、D） |
| `Geom=AllCheck` | 结构、title、电荷和多重度全部来自 chk，所以 route 后面直接是 ReadIsotopes 输入节 |

### 6.2 方法二：在坐标中写 Iso=

[d2o_freq_iso_inline.gjf](../examples/06-freq-thermo/d2o_freq_iso_inline.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=d2o_freq_iso_inline.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt Freq

D2O opt + freq with inline isotope specification H(Iso=2)

0 1
O           0.000000     0.000000     0.117300
H(Iso=2)    0.000000     0.757200    -0.469200
H(Iso=2)    0.000000    -0.757200    -0.469200

```

这种写法直观，但需要从头计算。两种方法得到的 D2O 频率应当相同。

```bash
grep "Frequencies --" d2o_freq_readisotopes.log d2o_freq_iso_inline.log
```

D2O 的气相实验基频约为 1178、2671、2788 cm⁻¹（NIST CCCBDB），与 H2O 相比，伸缩模式的频率比接近 $1/\sqrt{2}$ 的谐振子极限（实际略高，因为 O 原子也参与振动）。

> [!TIP]
> 同样的 `ReadFC` 技巧可以用来**不重算 Hessian 就改温度和压力**（第 7.3 节），这对动力学同位素效应（KIE）的计算尤其方便：一次频率计算，多次换同位素。

---

## 7. 热化学输出解读

### 7.1 输出块

频率计算后，日志中有 `- Thermochemistry -` 一节（示意输出，H2O 平动项按理想气体公式可以严格算出；其余数值仅供参考）：

```text
 -------------------
 - Thermochemistry -
 -------------------
 Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
 ...
 Rotational symmetry number  2.
 ...
 Zero-point correction=                           0.02xxxx (Hartree/Particle)
 Thermal correction to Energy=                    0.02xxxx
 Thermal correction to Enthalpy=                  0.02xxxx
 Thermal correction to Gibbs Free Energy=         0.00xxxx
 Sum of electronic and zero-point Energies=          -76.4xxxxx
 Sum of electronic and thermal Energies=             -76.4xxxxx
 Sum of electronic and thermal Enthalpies=           -76.4xxxxx
 Sum of electronic and thermal Free Energies=        -76.4xxxxx

                     E (Thermal)             CV                S
                      KCal/Mol        Cal/Mol-Kelvin    Cal/Mol-Kelvin
 Total                   1x.xxx              6.0xx             45.xxx
 Electronic               0.000              0.000              0.000
 Translational            0.889              2.981             34.608
 Rotational               0.889              2.981             10.xxx
 Vibrational             1x.xxx              0.0xx              0.0xx
```

| 行 | 含义 | 关系 |
|---|---|---|
| `Zero-point correction` | ZPE（Hartree/分子） | $\frac{1}{2}\sum h c\tilde\nu$ |
| `Thermal correction to Energy` | 0 K → T 的内能校正（含 ZPE） | $E_{\mathrm{corr}}$ |
| `Thermal correction to Enthalpy` | 焓校正 | $E_{\mathrm{corr}} + k_{\mathrm B}T$ |
| `Thermal correction to Gibbs Free Energy` | 自由能校正 | $H_{\mathrm{corr}} - TS$ |
| `Sum of electronic and zero-point Energies` | $E_0 = E_{\mathrm{elec}} + \mathrm{ZPE}$ | |
| `Sum of electronic and thermal Energies` | $E = E_{\mathrm{elec}} + E_{\mathrm{corr}}$ | |
| `Sum of electronic and thermal Enthalpies` | $H = E_{\mathrm{elec}} + H_{\mathrm{corr}}$ | |
| `Sum of electronic and thermal Free Energies` | $G = E_{\mathrm{elec}} + G_{\mathrm{corr}}$ | **最常用的“自由能”** |

其中 $E_{\mathrm{elec}}$ 就是本作业的 `SCF Done` 能量。下方表格把内能、热容和熵分解为电子、平动、转动、振动贡献，单位分别是 kcal/mol 和 cal/(mol·K)。

```bash
grep -E "Zero-point correction|Thermal correction to|Sum of electronic and" job.log
grep -E "Temperature .* Kelvin|Rotational symmetry number" job.log
```

### 7.2 两个容易忽略的检查

**转动对称数。** `Rotational symmetry number` 行由识别出的点群决定。用了 `NoSymm`，或结构对称性没有被识别出来时，该值可能变成 1，转动熵就会多出 $R\ln\sigma$。以水为例，σ 从 2 变成 1，熵增加 $R\ln 2 = 1.38$ cal/(mol·K)，298.15 K 下 G 降低约 0.41 kcal/mol。高对称分子（CH4 σ = 12、苯 σ = 12）的误差更大。

**低频模式警告。** 存在很低的振动频率时，热化学部分会出现形如 `Warning -- explicit consideration of N degrees of freedom as vibrations may cause significant error` 的提示，说明有模式更像受阻转动，RRHO 熵不可靠，需要第 9 节的 qRRHO 校正。

### 7.3 改变温度、压力与缩放因子

**直接在 route 中指定。** [freq_temperature_pressure.gjf](../examples/06-freq-thermo/freq_temperature_pressure.gjf)：甲醇在 373.15 K、2.0 atm 下的热化学：

```text
%nprocshared=8
%mem=12GB
%chk=freq_temperature_pressure.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq Temperature=373.15 Pressure=2.0

Methanol opt + freq, thermochemistry at 373.15 K and 2.0 atm

0 1
C    -0.370577     0.019773    -0.009266
O     0.900681    -0.531378    -0.299847
H    -0.521270     0.039918     1.072865
H    -0.431098     1.032498    -0.415279
H    -1.141363    -0.602000    -0.470688
H     1.563626     0.041188     0.122216

```

`Temperature=` 和 `Pressure=` 只影响热化学分析，不影响结构、能量和频率本身。检查日志中的 `Temperature   373.150 Kelvin.  Pressure   2.00000 Atm.` 一行。

**不重算 Hessian，从 chk 重新分析：`Freq=ReadFC`。** [freq_readfc_rethermo.gjf](../examples/06-freq-thermo/freq_readfc_rethermo.gjf)：

```text
%nprocshared=8
%mem=12GB
%oldchk=h2o_opt_freq_raman.chk
%chk=h2o_readfc_rethermo.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Freq=ReadFC Geom=AllCheck Guess=Read Temperature=298.15 Pressure=24.4654

```

`Geom=AllCheck` 时 title、电荷多重度和结构都来自 chk，所以 route 后面只有一个空行，**不能**再写 title。压力取 24.4654 atm，是因为 298.15 K 下理想气体浓度为 1 mol/L 时的压力 $P = cRT$ 正好是这个值，所以这个作业直接给出 **1 M 标准态**的自由能（第 10 节）。注意：换温度后对应的压力也要随之改变（$P = cRT$）。

`Scale=` 关键词可以给热化学分析指定频率缩放因子（官方 Freq 文档：Temperature、Pressure、Scale 都可以写在 route 中），例如 B3LYP/6-31G(d) 的 ZPE 缩放：`Scale=0.9806`（Scott & Radom 1996）。Gaussian 只接受一个缩放因子，用于 ZPE 和所有热项。

**freqchk 工具。** Gaussian 自带的 freqchk 可以从 chk 中交互式地重做振动和热化学分析：

```bash
freqchk h2o_opt_freq_raman.chk
```

程序依次询问温度、压力、缩放因子、同位素等（直接回车为默认值），结果打印到屏幕。适合快速查看不同温度下的结果，不需要提交作业。工具列表见 [Gaussian Utilities](https://gaussian.com/utils/)。

---

## 8. 频率缩放因子

谐振近似忽略了非谐性，加上方法本身的误差，计算频率系统性偏高。经验做法是乘以缩放因子，但**不同用途的缩放因子不同**：

| 用途 | 缩放因子类型 | 说明 |
|---|---|---|
| 与实验 IR/Raman 基频比较、画光谱 | 基频（fundamental, F） | 使谐振频率接近实验基频 |
| ZPE（热化学、动力学） | ZPE 因子 | 一般比基频因子更接近 1 |
| 与谐振频率比较 | 谐振（H）因子 | 主要用于方法基准测试 |

同一方法基组在不同数据库中的数值示例：

| 方法/基组 | 基频因子 | ZPE 因子 | 来源 |
|---|---|---|---|
| HF/6-31G(d) | 0.8953 | 0.9135 | Scott & Radom 1996 |
| MP2(full)/6-31G(d) | 0.9427 | 0.9661 | Scott & Radom 1996 |
| B3LYP/6-31G(d) | 0.9614 | 0.9806 | Scott & Radom 1996 |
| B3LYP/6-31G(d) | — | 0.977 | Truhlar 组数据库，引自 GoodVibes 文档示例 |

可以看到同一水平下 ZPE 因子与基频因子不同，不同数据库给出的数值也略有差别。其他方法基组请查：

- [Truhlar 组 Database of Frequency Scale Factors（v5）](https://comp.chem.umn.edu/freqscale/)：按方法基组给出 ZPE、F、H 三类因子；
- [NIST CCCBDB 振动缩放因子](https://cccbdb.nist.gov/vibscalejust.asp)；
- [liyuanhe211/Collection_of_Frequency_Scale_Factors](https://github.com/liyuanhe211/Collection_of_Frequency_Scale_Factors)：多来源缩放因子的汇总（仅外链）。

> [!WARNING]
> 数据库里查不到你用的组合时，不要套用别的基组的数值。可以选用最接近的组合并在论文中注明，或者不缩放直接报告谐振结果。热化学中 ZPE 缩放带来的变化通常只有零点几 kcal/mol，但对含很多 X–H 键的大分子会累积。

---

## 9. 低频模式与准谐振（qRRHO）校正

### 9.1 问题

RRHO 中振动熵随频率降低而发散（2.3 节公式），频率低于约 100 cm⁻¹ 的模式（甲基转动、分子间平动和转动、柔性链的扭转）实际上更接近受阻转子。按谐振子处理会**高估熵**，使 $-TS$ 过负，对以下情形影响特别大：

- 两分子结合成一个分子的步骤（A + B → AB、双分子过渡态）：ΔG 和 ΔG‡ 偏差可达数 kcal/mol；
- 柔性大分子、弱相互作用复合物的结合自由能。

### 9.2 Grimme 的 qRRHO 熵

Grimme（2012）用阻尼函数在谐振子熵 $S_V$ 与自由转子熵 $S_R$ 之间插值：

$$
S = w(\tilde\nu)\,S_V + \left[1 - w(\tilde\nu)\right] S_R,
\qquad
w(\tilde\nu) = \frac{1}{1 + (\tilde\nu_0/\tilde\nu)^{\alpha}},\quad \tilde\nu_0 = 100\ \mathrm{cm^{-1}},\ \alpha = 4
$$

$$
S_R = R\left[\frac{1}{2} + \ln\left(\frac{8\pi^3\mu' k_{\mathrm B}T}{h^2}\right)^{1/2}\right],
\qquad
\mu' = \frac{\mu B_{\mathrm{av}}}{\mu + B_{\mathrm{av}}},\quad \mu = \frac{h}{8\pi^2\nu}
$$

其中 $B_{\mathrm{av}} = 10^{-44}$ kg·m² 是平均转动惯量上限。高频模式 $w\to 1$，仍按谐振子处理；低频模式逐渐过渡到自由转子。Head-Gordon 等人对焓（振动内能）提出了类似的插值。另一种常用做法是 **Truhlar 准谐振**：把低于阈值（通常 100 cm⁻¹）的频率直接提升到 100 cm⁻¹ 再计算。

### 9.3 工具

**GoodVibes**（v4.4.0，2026 年 9 月，MIT 许可）：

```bash
pip install goodvibes            # or: conda install -c conda-forge goodvibes

# Grimme qRRHO entropy + Head-Gordon enthalpy, cutoff 100 cm-1 (gas, 1 atm)
goodvibes h2o_opt_freq_raman.log -q -f 100
# Truhlar-type quasi-harmonic treatment instead of Grimme
goodvibes h2o_opt_freq_raman.log -q --qs truhlar -f 100
```

输出列：`E`（电子能）、`ZPE`、`H`、`T.S`、`T.qh-S`（qRRHO 熵项）、`G(T)`（RRHO 自由能）、`qh-G(T)`（准谐振自由能），单位均为 Hartree。根据 GoodVibes 文档，v4 版本会按日志中识别到的理论水平**自动套用** Truhlar 数据库的缩放因子；加 `-v 1.0` 可以关闭缩放，便于和 Gaussian 自己的输出直接对照。全部选项用 `goodvibes -h` 查看，完整脚本见示例目录中的 [goodvibes_run.sh](../examples/06-freq-thermo/goodvibes_run.sh)。

**Shermo**（Sobereva 开发，中文用户常用）：通过 `settings.ini` 设置温度、压力或浓度、缩放因子和低频处理方式，支持批量处理与高精度单点能组合，见 [Shermo 主页](http://sobereva.com/soft/shermo/)。批量调用可以参考 [easy-Shermo](https://github.com/kimariyb/easy-Shermo)。

**OpenThermo**：C++ 实现，支持 qRRHO 和温度、压力扫描，见 [lenhanpham/OpenThermo](https://github.com/lenhanpham/OpenThermo)。

**本仓库脚本**：[thermo_qrrho.py](../scripts/thermo_qrrho.py) 是一个教学用的精简实现（Grimme qRRHO、1 M 校正、缩放、单点组合、构象 Boltzmann 加权），便于对照公式理解每一项（`python ../scripts/thermo_qrrho.py -h`）。正式研究推荐使用 GoodVibes 或 Shermo。

> [!TIP]
> 论文中应写明：是否做了 qRRHO 校正、用的哪种方法（Grimme/Truhlar）、阈值多少、是否缩放，以及使用的程序和版本。同一个反应剖面上的所有物种必须用同一套处理方式。

---

## 10. 标准态校正：1 atm → 1 mol/L

Gaussian 的热化学默认是 **1 atm 理想气体**标准态。溶液中的反应通常以 **1 mol/L** 为标准态。两者只差平动熵中的体积项：

$$
\Delta G^{\circ}_{1\,\mathrm{atm}\to 1\,\mathrm{M}} = RT\ln\frac{V_{1\,\mathrm{atm}}}{V_{1\,\mathrm{M}}} = RT\ln\frac{RT/P^{\circ}}{1/c^{\circ}} = RT\ln(24.46)
$$

298.15 K 下，1 mol 理想气体在 1 atm 下的体积为 24.465 L，所以每个物种的 G 增加：

$$
RT\ln(24.46) = 0.5925 \times 3.197 \approx 1.89\ \mathrm{kcal/mol}\ (7.93\ \mathrm{kJ/mol})
$$

| 情形 | 校正是否影响 ΔG |
|---|---|
| A → B（分子数不变） | 不影响，校正相互抵消 |
| A + B → C | ΔG 变化 −1.89 kcal/mol，即更有利（反应物两个物种各加 1.89，产物只加一次） |
| 双分子过渡态 A + B → [AB]‡ | ΔG‡ 同样受影响 |
| 溶剂分子作为反应物（水在水中约 55.3 M） | 还需额外 $RT\ln(55.3) \approx 2.38$ kcal/mol |

实现方式（三选一，**不要重复加**）：

1. GoodVibes 加 `--conc 1.0`（根据其文档，这一选项会给每个物种的 G 加上 298 K 下的 1.89 kcal/mol）；
2. Shermo 中设置浓度；
3. Gaussian 中写 `Pressure=24.4654`（仅对 298.15 K 成立，第 7.3 节），或者手动给每个物种的 G 加 1.89 kcal/mol（0.003019 Hartree）。

---

## 11. 高精度单点能 + 低水平热校正

### 11.1 组合公式

DFT 的电子能对泛函和基组敏感，而热校正（ZPE、H_corr、G_corr）对理论水平相对不敏感。所以通行做法是：

$$
G = E_{\mathrm{SP}}^{\mathrm{high}} + G_{\mathrm{corr}}^{\mathrm{low}}\ \ (+\ \Delta G^{\circ}_{1\,\mathrm{atm}\to 1\,\mathrm{M}})
$$

其中 $G_{\mathrm{corr}}^{\mathrm{low}}$ 是低水平 Freq 的 `Thermal correction to Gibbs Free Energy`（或 qRRHO 校正后的值），$E_{\mathrm{SP}}^{\mathrm{high}}$ 是同一结构上高水平单点的 `SCF Done`（双杂化泛函要用含 MP2 部分的最终能量，见[第4章](04-sp-scf.md)）。焓同理：$H = E_{\mathrm{SP}}^{\mathrm{high}} + H_{\mathrm{corr}}^{\mathrm{low}}$。

> [!CAUTION]
> **不能**用 `Sum of electronic and thermal Free Energies` 加高水平单点，那一项已经包含了低水平的电子能。要用的是 `Thermal correction to Gibbs Free Energy` 这一项。

### 11.2 示例：Link1 一次完成

[highlevel_sp_plus_thermo.gjf](../examples/06-freq-thermo/highlevel_sp_plus_thermo.gjf)：

```text
%nprocshared=8
%mem=12GB
%chk=highlevel_sp_plus_thermo.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

Ethanol step 1: B3LYP-D3(BJ)/def2-SVP opt + freq (thermal corrections)

0 1
C    -0.886284     0.164102    -0.072797
C     0.472700    -0.504171    -0.048378
O     1.385279     0.279434     0.706651
H    -0.821053     1.157324    -0.529101
H    -1.605682    -0.435781    -0.637136
H    -1.264380     0.305590     0.945012
H     0.868553    -0.625847    -1.061187
H     0.405268    -1.491959     0.416669
H     1.445597     1.151307     0.280268

--Link1--
%nprocshared=8
%mem=12GB
%chk=highlevel_sp_plus_thermo.chk
#p wB97XD/def2TZVP Geom=Check Guess=Read

Ethanol step 2: wB97X-D/def2-TZVP single point on the optimized geometry

0 1

```

提取并组合：

```bash
grep "Thermal correction to Gibbs Free Energy" highlevel_sp_plus_thermo.log   # G_corr (low level)
grep "SCF Done" highlevel_sp_plus_thermo.log | tail -n 1                      # E_SP (high level, last SCF)
```

```python
# combine high-level SP with low-level thermal correction (numbers are illustrative)
HARTREE2KCAL = 627.5095
E_sp_high = -155.012345   # Eh, "SCF Done" of the wB97XD/def2TZVP step
G_corr_low = 0.052345     # Eh, "Thermal correction to Gibbs Free Energy" of the B3LYP step
dG_conc = 1.894 / HARTREE2KCAL   # 1 atm -> 1 M at 298.15 K
G = E_sp_high + G_corr_low + dG_conc
print(f"G(1 M) = {G:.6f} Eh")
```

或者直接交给 GoodVibes，它能识别 Link1 中的单点：

```bash
goodvibes highlevel_sp_plus_thermo.log -q --spc link --conc 1.0
```

### 11.3 反应自由能与构象加权

- 反应自由能：$\Delta G = \sum G(\text{产物}) - \sum G(\text{反应物})$，乘以 627.5095 换成 kcal/mol；
- 能垒：$\Delta G^{\ddagger} = G(\mathrm{TS}) - G(\text{反应物或前驱复合物})$，TS 的热校正自动排除了虚频（[第7章](07-ts-irc.md)）；
- 多构象体系按 Boltzmann 分布加权：

$$
p_i = \frac{\exp(-\Delta G_i/RT)}{\sum_j \exp(-\Delta G_j/RT)}
$$

GoodVibes 的 `--boltz` 和 thermo_qrrho.py 都能做构象加权。完整研究流程见[第18章](18-errors-capstone.md)。

> [!IMPORTANT]
> 同一项研究中所有物种的热校正必须来自**同一个**低水平（同一方法、基组、溶剂模型、格点、qRRHO 设置），单点也必须来自同一个高水平。论文的计算方法部分要写清楚，例如“ωB97X-D/def2-TZVP//B3LYP-D3(BJ)/def2-SVP，热校正在 B3LYP-D3(BJ)/def2-SVP 水平下用 Grimme qRRHO（100 cm⁻¹）计算，自由能已校正到 1 mol/L 标准态”。

---

## 参数选择与收敛建议

| 情形 | 推荐设置 | 说明 |
|---|---|---|
| 常规极小点或过渡态验证 | `Opt Freq`（同一 route） | 自动保证理论水平一致 |
| 低频模式重要（弱复合物、柔性分子） | `Opt=Tight` + `Int=SuperFine`，频率用同样的格点 | 降低小虚频与低频噪声 |
| 需要 Raman | `Freq=Raman` | DFT 默认不算 |
| 需要高精度简正模式 | `Freq=HPModes` | 五位有效数字 |
| 与实验 IR 精确对比的小分子 | `Opt=VeryTight Int=SuperFine Freq=Anharmonic` | 不需要缩放 |
| 同位素、改温度或压力 | `Freq=(ReadFC,ReadIsotopes)` 或 `Freq=ReadFC Temperature=... Pressure=...` | 不重算 Hessian |
| 画光谱 | 基频缩放因子 + 洛伦兹展宽 | 第 5、8 节 |
| ZPE 或热化学缩放 | ZPE 缩放因子（`Scale=` 或 GoodVibes `-v`/`--zpe-vscal`） | 与基频因子不同 |
| 溶液反应自由能 | qRRHO + 1 M 标准态校正 | GoodVibes `-q --conc 1.0` |
| 高精度能量 | 高水平单点 + 低水平 G_corr | 第 11 节 |
| 用了 NoSymm | 检查 `Rotational symmetry number` | σ 错误会导致熵误差 |

---

## 常见问题与报错

**1. 优化后出现小虚频（例如 −20 cm⁻¹）**

- 原因：收敛不够紧、格点噪声、对称性约束。
- 处理：第 4.2 节，`Opt=(ReadFC,Tight)` + `Int=SuperFine`，或者沿虚频位移（GaussView、pyQRC）后重新优化。

**2. Freq 步末尾的收敛表出现 `NO`**

- 原因：优化用近似 Hessian 判断已收敛，但用精确 Hessian 判断尚未收敛（[第5章](05-opt.md) 3.4 节）。
- 处理：此时频率和热化学不可靠，用 `Opt=CalcAll` 或 `Opt=(ReadFC,Tight)` 重新优化。

**3. 频率输出里没有 `Raman Activ` 行**

- 原因：DFT 默认不计算 Raman 活性。
- 处理：加 `Freq=Raman`；已有频率的 chk 时不能靠 ReadFC 补算 Raman，需要重新计算。

**4. 不同方法、基组之间“频率对不上”，或者在别人的结构上直接算频率**

- 原因：频率只在当前方法的驻点上有意义（第 2.1 节）。
- 处理：先在同一水平下重新优化，再算频率。

**5. `Freq=ReadFC` 读不到力常数，或基组报错**

- 原因：chk 中没有力常数（来自 Opt 而不是 Freq），或者在 ReadFC 作业中又输入了 Gen 基组。
- 处理：确认 chk 来自 Freq 作业；ReadFC 时基组从 chk 读，不要再写 Gen 输入节。

**6. 用了 `NoSymm` 后自由能与预期相差约 0.4 kcal/mol（或更多）**

- 原因：转动对称数变成 1。
- 处理：去掉 NoSymm 重算频率，或者用 GoodVibes、Shermo 等工具指定正确的对称数后重新做热化学。

**7. 冻结原子的模型算出很多低频或虚频**

- 原因：冻结的原子不在驻点上，部分 Hessian 只对活动区有意义。
- 处理：只用活动区的频率做定性判断；热化学数值仅为近似，需要在论文中说明。

**8. 频率作业内存不足：`galloc: could not allocate memory` 或作业被调度系统杀掉**

- 原因：%mem 超过了作业实际申请的内存（被杀），或者频率计算需要比单点更多的内存。
- 处理：%mem 设为 Slurm `--mem` 的 75–85%；分子较大时适当增加申请的内存，见[第17章](17-hpc-jobs.md)。

更多报错见[第18章](18-errors-capstone.md)。

---

## 练习

1. 运行 `h2o_opt_freq_raman.gjf`，把三个谐振频率与实验基频比较，讨论 B3LYP/def2-TZVP 的谐振频率为什么偏高；用 spectrum_ir_raman.py 画出 IR 和 Raman 谱。
2. 先后运行 `d2o_freq_readisotopes.gjf` 和 `d2o_freq_iso_inline.gjf`，比较两者的频率和 ZPE，计算 H2O 与 D2O 的 ZPE 差。
3. 运行 `freq_readfc_rethermo.gjf`，比较 1 atm 和 24.4654 atm 两种情况下的 `Thermal correction to Gibbs Free Energy`，验证差值约为 1.89 kcal/mol。
4. 对 `highlevel_sp_plus_thermo.log` 分别运行 `goodvibes -q --spc link` 和不加 `-q`，比较 `G(T)` 与 `qh-G(T)` 的差别；再对一个含甲基转动的较大分子做同样比较，体会低频模式的影响。
5. 运行 `formaldehyde_freq_anharmonic.gjf`，把 `E(harm)` 与 `E(anharm)` 分别和实验基频比较，计算两种情况下的平均绝对误差。

---

## 参考资料

**官方文档**

- [Gaussian 16 Freq keyword](https://gaussian.com/freq/)：Raman、Anharmonic、ReadFC、ReadIsotopes、HPModes、DiagFull，以及 Temperature/Pressure/Scale
- [Thermochemistry in Gaussian（白皮书）](https://gaussian.com/thermo/)：热化学各项的计算公式
- [Vibrational Analysis in Gaussian（白皮书）](https://gaussian.com/vib/)：简正模式分析的实现
- [Opt keyword](https://gaussian.com/opt/)：Tight/VeryTight、CalcAll
- [Integral keyword](https://gaussian.com/integral/)：SuperFine 格点
- [Gaussian Utilities（freqchk 等）](https://gaussian.com/utils/)
- [GoodVibes documentation](https://goodvibespy.readthedocs.io/)
- [Shermo](http://sobereva.com/soft/shermo/)
- [Truhlar 组频率缩放因子数据库](https://comp.chem.umn.edu/freqscale/)，[NIST CCCBDB 缩放因子](https://cccbdb.nist.gov/vibscalejust.asp)

**相关 GitHub 仓库（仅外链，未复制其代码或文本）**

- [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)：qRRHO、1 M 标准态校正、缩放因子库、单点组合、Boltzmann 加权（v4.4.0）
- [patonlab/pyQRC](https://github.com/patonlab/pyQRC)：沿虚频方向位移生成新输入，消除残留虚频
- [lenhanpham/OpenThermo](https://github.com/lenhanpham/OpenThermo)：C++ 热化学程序，支持 qRRHO 与温度、压力扫描
- [kimariyb/easy-Shermo](https://github.com/kimariyb/easy-Shermo)：批量调用 Shermo，自动配对频率与单点文件（中文 README）
- [eljost/thermoanalysis](https://github.com/eljost/thermoanalysis)：简短的 Python 热化学实现，适合对照白皮书理解公式
- [chauncey-garrett/gaussian-tools](https://github.com/chauncey-garrett/gaussian-tools)：从 log 中提取虚频、热校正等单项数据的小脚本集合
- [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py)：含热化学、沿振动模式位移等功能的 Python 工具包

**关键文献**

- A. P. Scott, L. Radom, Harmonic vibrational frequencies: An evaluation of Hartree-Fock, Møller-Plesset, quadratic configuration interaction, density functional theory, and semiempirical scale factors, *J. Phys. Chem.* **1996**, 100, 16502. DOI: [10.1021/jp960976r](https://doi.org/10.1021/jp960976r)
- I. M. Alecu, J. Zheng, Y. Zhao, D. G. Truhlar, Computational thermochemistry: Scale factor databases and scale factors for vibrational frequencies obtained from electronic model chemistries, *J. Chem. Theory Comput.* **2010**, 6, 2872. DOI: [10.1021/ct100326h](https://doi.org/10.1021/ct100326h)
- S. Grimme, Supramolecular binding thermodynamics by dispersion-corrected density functional theory, *Chem. Eur. J.* **2012**, 18, 9955. DOI: [10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)
- R. F. Ribeiro, A. V. Marenich, C. J. Cramer, D. G. Truhlar, Use of solution-phase vibrational frequencies in continuum models for the free energy of solvation, *J. Phys. Chem. B* **2011**, 115, 14556. DOI: [10.1021/jp205508z](https://doi.org/10.1021/jp205508z)
- Y.-P. Li, J. Gomes, S. M. Sharada, A. T. Bell, M. Head-Gordon, Improved force-field parameters for QM/MM simulations of the energies of adsorption for molecules in zeolites and a free rotor correction to the rigid rotor harmonic oscillator model for adsorption enthalpies, *J. Phys. Chem. C* **2015**, 119, 1840. DOI: [10.1021/jp509921r](https://doi.org/10.1021/jp509921r)
- G. Luchini, J. V. Alegre-Requena, I. Funes-Ardoiz, R. S. Paton, GoodVibes: automated thermochemistry for heterogeneous computational chemistry data, *F1000Research* **2020**, 9, 291. DOI: [10.12688/f1000research.22758.1](https://doi.org/10.12688/f1000research.22758.1)
- T. Lu, Q. Chen, Shermo: A general code for calculating molecular thermochemistry properties, *Comput. Theor. Chem.* **2021**, 1200, 113249. DOI: [10.1016/j.comptc.2021.113249](https://doi.org/10.1016/j.comptc.2021.113249)
- V. Barone, Anharmonic vibrational properties by a fully automated second-order perturbative approach, *J. Chem. Phys.* **2005**, 122, 014108. DOI: [10.1063/1.1824881](https://doi.org/10.1063/1.1824881)

---

[← 上一章](05-opt.md) | [返回目录](../README.md) | [下一章 →](07-ts-irc.md)
