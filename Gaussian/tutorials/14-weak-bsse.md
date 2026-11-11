# 第14章 弱相互作用与 BSSE（Counterpoise）（Noncovalent Interactions & Counterpoise BSSE Correction in Gaussian）

本章讲如何在 Gaussian 16 中可靠地计算分子间弱相互作用能（氢键、π-π 堆积、主客体、阳离子-偶极等）：超分子法相互作用能的定义，基组重叠误差（BSSE）的来源，以及用 `Counterpoise=2` 做 Boys-Bernardi 均衡校正的输入写法与输出读取（Counterpoise corrected energy、BSSE energy）。本章还会讲 CP 校正优化、形变能与结合能分解、色散校正与弥散函数的选择，并介绍用 Multiwfn 的 NCI/IGMH 方法把弱相互作用画出来（noncovalent interaction, BSSE, counterpoise correction, IGMH）。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 写出 ΔE_int、ΔE_bind、形变能、BSSE 的严格定义，知道每一项在哪个几何、哪个基组下计算；
- 正确写 `Counterpoise=N` 输入：每个原子的 `(Fragment=n)` 标记和“整体 + 各碎片”的电荷多重度行；
- 从输出中读出 `Counterpoise corrected energy`、`BSSE energy`、`sum of fragments`，换算出 kcal/mol 并自己核对；
- 判断什么时候需要 CP 校正优化（`Opt Counterpoise=2`），什么时候只需单点校正；
- 为弱相互作用选择合适的泛函、色散校正（GD3BJ/GD3）和基组；
- 用 Multiwfn 的 IGMH/NCI 分析把相互作用区域可视化。

| 体系类型 | 例子 | 主导作用 | 方法要点 |
|---|---|---|---|
| 氢键二聚体 | 水二聚体、甲酸二聚体、碱基对 | 静电 + 诱导 | 弥散函数、CP 校正 |
| π-π 堆积 | 苯二聚体、石墨烯片段、药物-DNA 嵌入 | 色散 | 必须加色散校正 |
| 主客体 | 冠醚·K⁺、环糊精、葫芦脲 | 静电 + 色散 | 形变能不可忽略 |
| 阳离子-π / 阴离子-π | 苯·Na⁺、C6F6·Cl⁻ | 静电 + 极化 | 阴离子要弥散函数 |
| 卤键、硫族键 | C–I···N、S···O | σ-hole 静电 + 色散 | 重元素基组/赝势 |
| 吸附（簇模型） | 小分子在金属/氧化物簇上 | 多种 | 碎片划分要合理 |

> [!NOTE]
> 本章示例目录：[../examples/14-cp-water-dimer/README.md](../examples/14-cp-water-dimer/README.md)（CP 单点与基组收敛）、[../examples/14-cp-benzene-dimer/README.md](../examples/14-cp-benzene-dimer/README.md)（CP 校正优化）、[../examples/14-cp-crown-k/README.md](../examples/14-cp-crown-k/README.md)（主客体与形变能）。碎片与电荷多重度行的通用语法见 [第2章](02-input-anatomy.md)，泛函与色散关键词见 [第3章](03-method-basis.md)。

---

## 2. 原理速览

### 2.1 超分子法相互作用能

对复合物 AB，记号 $E_X^{\,Y}(Z)$ 表示：体系 X，在基组 Y（Y 中所有原子的基函数）下，几何取 Z 中的构型。最常用的几个量：

$$
\Delta E_{\rm int}^{\rm raw} = E_{AB}^{\,AB}(AB) - E_{A}^{\,A}(AB) - E_{B}^{\,B}(AB)
$$

$$
\Delta E_{\rm int}^{\rm CP} = E_{AB}^{\,AB}(AB) - E_{A}^{\,AB}(AB) - E_{B}^{\,AB}(AB)
$$

$$
{\rm BSSE} = \left[E_{A}^{\,A}(AB) - E_{A}^{\,AB}(AB)\right] + \left[E_{B}^{\,B}(AB) - E_{B}^{\,AB}(AB)\right] \;>\; 0
$$

- 相互作用能 ΔE_int：单体**保持复合物中的几何**，只反映“把两个已经变形的单体放在一起”的能量变化；
- 结合能 ΔE_bind：相对于各自**优化好的**自由单体，需要加上形变能（§2.3）。

### 2.2 BSSE 从哪里来

有限基组下，复合物中的 A 可以“借用” B 的基函数来改善自身的描述，使 $E_{AB}$ 被人为降低，相互作用被高估——这就是基组重叠误差（basis set superposition error）。Boys-Bernardi 均衡校正（counterpoise, CP；[Boys & Bernardi 1970](https://doi.org/10.1080/00268977000101561)）的做法是：单体也在整个复合物基组下计算，另一个单体的位置只放基函数（鬼原子，ghost atom），没有核电荷和电子。

BSSE 的特点：

- 基组越小越严重，越接近完备基组越小；弥散函数的基组 BSSE 往往更大；
- 对 HF/DFT 已经不小，对 MP2、双杂化等含相关能的方法更明显；
- 色散校正（D3 等）只依赖原子坐标、与基组无关，本身不产生 BSSE；
- CP 校正有时被认为“过度校正”，所以也有人用 CP 校正值与未校正值的平均（half-counterpoise）；但在中小基组下，标准 CP 校正仍是最常见、最易复现的做法。

### 2.3 形变能与结合能分解

$$
E_{\rm def} = \left[E_{A}^{\,A}(AB) - E_{A}^{\,A}(A)\right] + \left[E_{B}^{\,B}(AB) - E_{B}^{\,B}(B)\right] \;\ge\; 0
$$

$$
\Delta E_{\rm bind}^{\rm CP} = \Delta E_{\rm int}^{\rm CP} + E_{\rm def}
= E_{AB}^{\,AB}(AB) - E_{A}^{\,A}(A) - E_{B}^{\,B}(B) + {\rm BSSE}
$$

其中 $E_X^{\,X}(X)$ 是自由单体在自身优化几何下的能量。刚性小分子（水二聚体）的形变能很小；柔性主体（冠醚、环糊精、蛋白口袋模型）的形变能可达数到十几 kcal/mol，绝不能忽略。若主体的自由构象与结合构象完全不同，还要说明自由主体用的是哪个构象。

### 2.4 色散校正与弥散函数

- **色散**：B3LYP、PBE0 等常规泛函几乎不描述长程色散，π-π 堆积、烷基链聚集等体系必须加色散校正。B3LYP、PBE0（`PBE1PBE`）用 `EmpiricalDispersion=GD3BJ`；M06-2X 等 Minnesota 泛函只有零阻尼参数，用 `EmpiricalDispersion=GD3`；`wB97XD`、`B97D3`、`B2PLYPD3`、`DSDPBEP86` 本身已含色散，**不要再加**（官方 DFT 与 MP 方法页面）。
- **弥散函数**：阴离子、氢键受体的孤对电子、极化率大的体系需要弥散函数（如 `6-311+G(d,p)`、`aug-cc-pVTZ`）。def2-TZVP 对中性体系的弱相互作用通常已可接受，加 CP 校正后结果更稳。

---

## 3. Counterpoise 关键词与输入格式

### 3.1 语法要点（G16 官方文档）

| 项目 | 规则 |
|---|---|
| 关键词 | `Counterpoise=N`，N 为碎片数（必须给出） |
| 碎片标记 | 每个原子写成 `O(Fragment=1)`；编号从 1 开始连续 |
| 电荷多重度行 | 先写整体，再按碎片编号依次写各碎片：`整体电荷 整体多重度 碎片1电荷 碎片1多重度 碎片2电荷 碎片2多重度` |
| 可用作业类型 | 单点、`Opt`、`Freq`、BOMD |
| 选项 | `NewGhost`（默认，推荐；鬼原子带 DFT 积分格点）/`OldGhost`（仅用于与旧结果比较） |
| 限制 | **不能与 `ONIOM` 或 `SCRF` 同时使用**；CP 作业不产生分子轨道（不能用于后续 `Guess=Read` 或波函数分析） |

典型的电荷多重度行：

| 体系 | 行 | 说明 |
|---|---|---|
| 水二聚体 | `0 1 0 1 0 1` | 两个中性单重态 |
| 18-冠-6·K⁺ | `1 1 0 1 1 1` | 整体 +1；碎片 1（冠醚）0/1；碎片 2（K⁺）+1/1 |
| 苯·Cl⁻ | `-1 1 0 1 -1 1` | 阴离子碎片 |
| 自由基·分子 | `0 2 0 2 0 1` | 整体二重态；自由基碎片二重态 |

> [!IMPORTANT]
> 碎片电荷之和必须等于整体电荷，每个碎片的多重度必须与它自己的电子数奇偶一致。漏标某个原子的 `Fragment`、碎片编号不连续或碎片电荷写错，都会导致报错或得到错误的 BSSE。

### 3.2 Gaussian 实际做了哪些计算

对 `Counterpoise=2`，Gaussian 依次做 5 个 SCF（每个都会打印一行 `SCF Done`），通常顺序为：

| 序号 | 体系 | 基组 | 对应符号 |
|---|---|---|---|
| 1 | 复合物 AB | 全部基函数 | $E_{AB}^{\,AB}(AB)$ |
| 2 | 碎片 1 | 全部基函数（碎片 2 为鬼原子） | $E_{A}^{\,AB}(AB)$ |
| 3 | 碎片 2 | 全部基函数（碎片 1 为鬼原子） | $E_{B}^{\,AB}(AB)$ |
| 4 | 碎片 1 | 自身基函数 | $E_{A}^{\,A}(AB)$ |
| 5 | 碎片 2 | 自身基函数 | $E_{B}^{\,B}(AB)$ |

输出中每个碎片计算前会有类似 `Counterpoise: doing DCBS calculation for fragment   1`（全基组）或 `... MCBS calculation ...`（单体基组）的提示行，核对顺序时以这些提示行为准。

---

## 4. 实战一：水二聚体 CP 单点与基组收敛

### 4.1 输入

文件 [../examples/14-cp-water-dimer/water_dimer_counterpoise.gjf](../examples/14-cp-water-dimer/water_dimer_counterpoise.gjf) 在同一结构下用三个基组做 CP 单点，再优化一个自由水分子用于形变能。结构为手工搭建的 Cs 型氢键二聚体（O···O 2.91 Å，O–H···O 近似线性，受体平面与供体平面垂直），不是文献坐标。第 1 步如下，第 2、3 步只把基组换成 `def2TZVP` 和 `aug-cc-pVTZ`：

```text
%nprocshared=8
%mem=8GB
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Counterpoise=2

Water dimer step 1: CP-corrected interaction energy, B3LYP-D3(BJ)/def2-SVP

0 1 0 1 0 1
O(Fragment=1)       0.000000    0.000000    0.000000
H(Fragment=1)       0.960000    0.000000    0.000000
H(Fragment=1)      -0.240365    0.929422    0.000000
O(Fragment=2)       2.910000    0.000000    0.000000
H(Fragment=2)       3.230100   -0.492911    0.759062
H(Fragment=2)       3.230100   -0.492911   -0.759062

```

第 4 步（自由单体）：

```text
--Link1--
%nprocshared=8
%mem=8GB
%chk=h2o_monomer.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt Freq

Water monomer step 4: relaxed monomer for the deformation energy

0 1
O                   0.000000    0.000000    0.000000
H                   0.960000    0.000000    0.000000
H                  -0.240365    0.929422    0.000000

```

| 关键词 / 位置 | 说明 |
|---|---|
| `Counterpoise=2` | 两个碎片：碎片 1 为氢键供体，碎片 2 为受体 |
| `EmpiricalDispersion=GD3BJ` | B3LYP 有 BJ 阻尼 D3 参数；色散对氢键体系也有约 10%–20% 的贡献 |
| `0 1 0 1 0 1` | 整体、碎片 1、碎片 2 均为中性单重态 |
| CP 步骤不写 `%chk` | CP 作业不产生可用于后续读取的分子轨道 |
| 三个基组 | 观察 BSSE 随基组增大而减小，CP 校正后的 ΔE_int 更快收敛 |

### 4.2 运行与读取

```bash
cd Gaussian/examples/14-cp-water-dimer
g16 < water_dimer_counterpoise.gjf > water_dimer_counterpoise.log
grep -c "Normal termination" water_dimer_counterpoise.log        # expect 4
grep -A4 "Counterpoise corrected energy" water_dimer_counterpoise.log
grep "Counterpoise: doing" water_dimer_counterpoise.log | head -4
```

每个 CP 步骤结束时会打印一个汇总块（以下为官方文档给出的输出格式，数值来自官方示例，与本例无关）：

```text
 Counterpoise corrected energy =    -622.483714884838
                  BSSE energy =       0.005664641553
             sum of fragments =    -622.333463706103
          complexation energy =     -97.84 kcal/mole (raw)
          complexation energy =     -94.28 kcal/mole (corrected)
```

| 输出行 | 含义 | 与 §2 符号的关系 |
|---|---|---|
| `Counterpoise corrected energy` | 复合物能量加上 BSSE | $E_{AB}^{\,AB}(AB) + {\rm BSSE}$ |
| `BSSE energy` | BSSE，正值，单位 Hartree | ${\rm BSSE}$ |
| `sum of fragments` | 各碎片在复合物几何、自身基组下的能量之和 | $E_{A}^{\,A}(AB)+E_{B}^{\,B}(AB)$ |
| `complexation energy ... (raw)` | 未校正的相互作用能 | $\Delta E_{\rm int}^{\rm raw}$ |
| `complexation energy ... (corrected)` | CP 校正的相互作用能 | $\Delta E_{\rm int}^{\rm CP}$ |

用官方示例数字自己核对一遍：$E_{\rm raw} = -622.483714885 - 0.005664642 = -622.489379527$；$(E_{\rm raw} - E_{\rm sum})\times 627.5095 = -97.84$ kcal/mol；$(E_{\rm CP}-E_{\rm sum})\times 627.5095 = -94.28$ kcal/mol，与输出一致。可见 Gaussian 打印的 complexation energy **相对于复合物几何下的单体**，即 ΔE_int，不含形变能。

### 4.3 用脚本汇总

示例目录中的 [../examples/14-cp-water-dimer/cp_summary.py](../examples/14-cp-water-dimer/cp_summary.py)（仅用标准库）逐个 Link1 步骤提取汇总块；给出自由单体能量后还会算形变能和结合能：

```bash
python cp_summary.py water_dimer_counterpoise.log
# monomer energy (Hartree) from step 4, used twice because both fragments are water
EMONO=$(grep "SCF Done" water_dimer_counterpoise.log | tail -1 | awk '{print $5}')
python cp_summary.py water_dimer_counterpoise.log --mono $EMONO $EMONO
```

注意：第 4 步的单体在 def2-TZVP 下优化，形变能只对第 2 步（同为 def2-TZVP）有意义；其他基组的步骤要用同基组的单体能量。

### 4.4 预期趋势与参考值

（定性预期，具体数值以你的计算为准）

| 基组 | BSSE | ΔE_int(raw) 与 ΔE_int(CP) 的差别 |
|---|---|---|
| def2-SVP | 最大，可达 1–2 kcal/mol 量级 | 很明显 |
| def2-TZVP | 明显减小 | 较小 |
| aug-cc-pVTZ | 最小 | 很小，两者接近 |

水二聚体相互作用能的高精度参考值约为 −5.0 kcal/mol（CCSD(T)/CBS，S22 基准集，[Jurečka et al. 2006](https://doi.org/10.1039/B600027D)；近似值，且对应其参考几何）。本例的几何是手工搭建的，与参考几何不同，只能做量级比较；正式比较时应使用基准集提供的几何。

---

## 5. 实战二：苯二聚体的 CP 校正优化

### 5.1 为什么做 CP 校正优化

普通优化是在“未校正”的势能面上找极小点。小基组下 BSSE 使单体人为靠近，平衡距离偏短、相互作用偏强。`Opt Counterpoise=2` 在每一步都用 CP 校正后的能量和梯度，相当于在 CP 势能面上优化。代价是每步要做 5 次 SCF（2 碎片时），耗时约为普通优化的 3–5 倍。

经验原则：

- 基组为三 ζ 及以上、体系较大时，通常做普通优化 + CP 单点即可；
- 双 ζ 基组下研究 π-π 堆积距离、弱络合物的结构时，CP 校正优化更合理；
- 无论哪种方式，都应在方法部分写清楚。

### 5.2 输入

文件 [../examples/14-cp-benzene-dimer/benzene_dimer_cp_opt.gjf](../examples/14-cp-benzene-dimer/benzene_dimer_cp_opt.gjf)：第 1 步 CP 校正优化，第 2 步同水平普通优化作对比，第 3 步优化自由苯。起始结构为平行错位（parallel-displaced）型：两个理想苯环（C–C 1.39 Å，C–H 1.09 Å）层间距 3.50 Å、水平错位 1.70 Å，是手工给的起始猜测。

```text
%nprocshared=16
%mem=16GB
%chk=bzd_cpopt.chk
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Counterpoise=2 NoSymm

Benzene dimer (parallel-displaced start) step 1: CP-corrected optimization

0 1 0 1 0 1
C(Fragment=1)       1.390000    0.000000    0.000000
C(Fragment=1)       0.695000    1.203775    0.000000
C(Fragment=1)      -0.695000    1.203775    0.000000
C(Fragment=1)      -1.390000    0.000000    0.000000
C(Fragment=1)      -0.695000   -1.203775    0.000000
C(Fragment=1)       0.695000   -1.203775    0.000000
H(Fragment=1)       2.480000    0.000000    0.000000
H(Fragment=1)       1.240000    2.147743    0.000000
H(Fragment=1)      -1.240000    2.147743    0.000000
H(Fragment=1)      -2.480000    0.000000    0.000000
H(Fragment=1)      -1.240000   -2.147743    0.000000
H(Fragment=1)       1.240000   -2.147743    0.000000
C(Fragment=2)       3.090000    0.000000    3.500000
C(Fragment=2)       2.395000    1.203775    3.500000
C(Fragment=2)       1.005000    1.203775    3.500000
C(Fragment=2)       0.310000    0.000000    3.500000
C(Fragment=2)       1.005000   -1.203775    3.500000
C(Fragment=2)       2.395000   -1.203775    3.500000
H(Fragment=2)       4.180000    0.000000    3.500000
H(Fragment=2)       2.940000    2.147743    3.500000
H(Fragment=2)       0.460000    2.147743    3.500000
H(Fragment=2)      -0.780000    0.000000    3.500000
H(Fragment=2)       0.460000   -2.147743    3.500000
H(Fragment=2)       2.940000   -2.147743    3.500000

```

| 关键词 | 说明 |
|---|---|
| `Opt Counterpoise=2` | CP 校正优化；官方 Counterpoise 页面的优化示例也采用这种写法 |
| `NoSymm` | 与官方 CP 优化示例一致，避免标准取向和对称性处理带来的麻烦 |
| def2-SVP | 故意用双 ζ 基组，使 CP 校正对结构的影响足够明显 |
| GD3BJ | 苯二聚体几乎完全由色散维系，不加色散 B3LYP 基本不成键 |

### 5.3 结果分析

```bash
cd Gaussian/examples/14-cp-benzene-dimer
grep -c "Normal termination" benzene_dimer_cp_opt.log      # expect 3
grep "Counterpoise corrected energy" benzene_dimer_cp_opt.log | tail -1
grep -c "Optimization completed" benzene_dimer_cp_opt.log  # expect 3
```

用 [../scripts/log2xyz.py](../scripts/log2xyz.py) 分别取出第 1、2 步的最终结构（`python ../scripts/log2xyz.py -h` 查看用法），在 GaussView 或 VMD 中测量两苯环质心之间的垂直距离和水平错位：

- CP 校正优化的层间距一般**大于**普通优化的结果；
- 两步的 ΔE_int(CP) 可以直接比较，普通优化结构上的 CP 能量应略高（不在 CP 势能面极小点）；
- 苯二聚体平行错位构型的高精度参考相互作用能约为 −2.6 ~ −2.7 kcal/mol（S22 系列基准，近似值），双 ζ 基组的结果只适合定性比较，定量结论请用 §8 推荐的单点水平。

> [!TIP]
> CP 校正优化后如果还要做频率，`Freq Counterpoise=2` 也是允许的（官方文档列出 Counterpoise 可用于频率计算），但耗时更高。很多研究只在普通优化结构上做频率，用于热校正。

---

## 6. 实战三：主客体复合物 18-冠-6·K⁺ 与形变能

### 6.1 流程

主客体体系的完整流程分两部分：

```text
host_guest_opt.gjf:   complex Opt Freq -> free host Opt Freq -> K+ single point   (B3LYP-D3(BJ)/def2-SVP)
host_guest_cp_sp.gjf: CP single points on the optimized complex                 (def2-TZVP, two functionals)
```

文件 [../examples/14-cp-crown-k/host_guest_opt.gjf](../examples/14-cp-crown-k/host_guest_opt.gjf) 的起始结构用 MMFF 力场预优化（K···O 2.7–2.9 Å），只作为 DFT 优化的起点。

### 6.2 CP 单点输入

文件 [../examples/14-cp-crown-k/host_guest_cp_sp.gjf](../examples/14-cp-crown-k/host_guest_cp_sp.gjf)（节选；43 个原子，冠醚 42 个原子为碎片 1，K⁺ 为碎片 2）：

```text
%nprocshared=16
%mem=16GB
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Counterpoise=2

18-crown-6 + K+ step 1: CP-corrected interaction energy, B3LYP-D3(BJ)/def2-TZVP
(replace coordinates by the optimized ones from host_guest_opt.log)

1 1 0 1 1 1
C(Fragment=1)      -0.173997   -3.468921   -1.127460
C(Fragment=1)      -1.637121   -3.267450   -0.730055
O(Fragment=1)      -1.907877   -1.863000   -0.620212
...                (42 crown-ether atoms in total, all Fragment=1)
K(Fragment=2)       0.005181    0.006849   -0.253118

--Link1--
%nprocshared=16
%mem=16GB
#p wB97XD/def2TZVP Counterpoise=2
...
```

| 要点 | 说明 |
|---|---|
| `1 1 0 1 1 1` | 整体 +1 单重态；冠醚中性单重态；K⁺ 为 +1 单重态（18 个电子） |
| 坐标 | 示例中是 MMFF 预优化坐标；**正式计算请换成 `host_guest_opt.log` 第 1 步的优化坐标** |
| 为什么不用 `Geom=Check` | 官方文档没有说明碎片信息是否随 chk 保存；最稳妥的做法是把优化坐标连同 `(Fragment=n)` 完整写出 |
| 第二个泛函 | `wB97XD` 已含色散，不加 `EmpiricalDispersion` |
| K 的基组 | def2 系列对 K 是全电子基组（def2 赝势从 Rb 开始） |

### 6.3 计算结合能

```bash
cd Gaussian/examples/14-cp-crown-k
python ../14-cp-water-dimer/cp_summary.py host_guest_cp_sp.log
grep "SCF Done" host_guest_opt.log | tail -1                  # K+ energy (step 3)
```

完整的结合能需要同一水平下的所有项。示例中优化在 def2-SVP、CP 单点在 def2-TZVP，所以形变能要在 def2-TZVP 下重新算：在 def2-TZVP 下对自由冠醚（第 2 步优化结构）和复合物中的冠醚构型各做一次单点，二者之差即冠醚的形变能（K⁺ 是单原子，形变能为 0）。

| 项 | 来源 |
|---|---|
| $\Delta E_{\rm int}^{\rm CP}$ | `host_guest_cp_sp.log`，`complexation energy (corrected)` |
| $E_{\rm def}$(冠醚) | def2-TZVP 单点：复合物中的冠醚构型 − 自由冠醚优化构型 |
| $\Delta E_{\rm bind}^{\rm CP}$ | 两者相加 |
| 热校正 | `host_guest_opt.log` 三步的 `Thermal correction to Gibbs Free Energy` |

> [!WARNING]
> 自由 18-冠-6 的最低能量构象与结合 K⁺ 时的 D3d 型构象不同。从复合物构型出发优化，往往只得到一个邻近的局部极小，形变能可能被低估。严格做法是对自由主体做构象搜索（例如 CREST，见 [第18章](18-errors-capstone.md)），取最低能量构象。

### 6.4 从电子能到结合自由能

$$
\Delta G_{\rm bind} = \Delta E_{\rm bind}^{\rm CP} + \Delta G_{\rm corr} + \Delta G_{\rm 1atm\to1M}
$$

- $\Delta G_{\rm corr}$：复合物与各单体 `Thermal correction to Gibbs Free Energy` 之差，建议用准谐振（qRRHO）处理低频，见 [第6章](06-freq-thermo.md)；
- 对 A + B → AB（分子数减少 1），1 atm → 1 mol/L 标准态校正为 $-RT\ln 24.46 \approx -1.89$ kcal/mol（298.15 K）；
- 溶液中的结合还要加溶剂化自由能差。由于 `Counterpoise` **不能与 `SCRF` 同时使用**，常用做法是：气相 CP 校正的 ΔE_bind，加上用 SMD 分别计算复合物与单体得到的溶剂化自由能差 $\Delta\Delta G_{\rm solv}$（[第9章](09-solvation.md)）。

---

## 7. 手动鬼原子计算（Bq）

需要单独检查某一项（例如只算 $E_{A}^{\,AB}(AB)$），或碎片划分很特殊时，可以手动写鬼原子：在元素符号后加分子力学类型 `Bq`，如 `O-Bq`、`H-Bq`。鬼原子带有该元素的基函数和积分格点，但没有核电荷和电子（G16 分子说明文档）。

水二聚体中“水 1 在二聚体基组下”的能量：

```text
%nprocshared=8
%mem=8GB
#p B3LYP/def2TZVP

Water 1 in the full dimer basis, water 2 as ghost atoms

0 1
O                   0.000000    0.000000    0.000000
H                   0.960000    0.000000    0.000000
H                  -0.240365    0.929422    0.000000
O-Bq                2.910000    0.000000    0.000000
H-Bq                3.230100   -0.492911    0.759062
H-Bq                3.230100   -0.492911   -0.759062

```

要点：

- 电荷与多重度只针对真实原子（这里是一个中性水分子）；
- 再算一次不带鬼原子的水 1 单点，二者之差就是碎片 1 对 BSSE 的贡献；
- 这里没有加色散校正：BSSE 是纯粹的基组效应，与 D3 无关，手动核对 BSSE 时不需要色散项；
- 日常使用推荐直接用 `Counterpoise=N`，不容易出错。

---

## 8. 方法与基组推荐（Benchmark 视角）

弱相互作用的常用基准集有 S22（[Jurečka et al. 2006](https://doi.org/10.1039/B600027D)）、S66（[Řezáč et al. 2011](https://doi.org/10.1021/ct2002946)）以及更大的综合基准 GMTKN55。它们给出了 CCSD(T)/CBS 级别的参考相互作用能和参考几何，可以用来检验自己的计算方案。

| 用途 | 推荐方案 | 说明 |
|---|---|---|
| 几何优化 | B3LYP-D3(BJ)/def2-SVP 或 def2-TZVP；ωB97X-D/def2-SVP | 双 ζ 基组下可考虑 CP 校正优化 |
| 常规单点 | B3LYP-D3(BJ)/def2-TZVP + CP | 性价比高，文献使用广泛 |
| 常规单点 | ωB97X-D/def2-TZVP + CP | 范围分离 + 色散，对电荷转移型复合物更稳 |
| 常规单点 | M06-2X-D3(0)/def2-TZVP + CP | `M062X EmpiricalDispersion=GD3`，注意 M06-2X 对积分格点较敏感 |
| 高精度单点 | DSD-PBEP86/def2-TZVP 或 def2-QZVP + CP | Gaussian 关键词 `DSDPBEP86`（已含 D3BJ，不要再加色散）；含 MP2 型相关，BSSE 比纯 DFT 大，CP 校正更重要 |
| 阴离子、σ-hole 受体 | 上述方法换成带弥散的基组（如 aug-cc-pVTZ、6-311+G(d,p)） | 弥散函数对阴离子必不可少 |

对应的 route 写法：

```text
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Counterpoise=2
#p wB97XD/def2TZVP Counterpoise=2
#p M062X/def2TZVP EmpiricalDispersion=GD3 Counterpoise=2
#p DSDPBEP86/def2TZVP Counterpoise=2
```

> [!TIP]
> 需要把相互作用能分解为静电、交换、诱导、色散等物理成分时，Gaussian 本身没有 SAPT。可选途径：NBO 7 接口的 NEDA（官方 Population 页面说明用 `NEDA=n` 指定碎片数，碎片定义方式与 Counterpoise 相同；需另行安装 NBO 7，见 [第12章](12-population-charges.md)）；或借助 MOKIT 的 `frag_guess_wfn` 生成 GAMESS 的 LMO-EDA 输入；或使用 Psi4 等程序做 SAPT。

---

## 9. 可视化分析：IGMH 与 NCI

能量回答“有多强”，图形化分析回答“在哪里、是什么类型”。Multiwfn 主功能 20 集中了弱相互作用可视化方法：其中子功能 11 为 IGMH（基于 Hirshfeld 划分的独立梯度模型，[Lu & Chen 2022](https://doi.org/10.1002/jcc.26812)），子功能 10 为 IGM；NCI（RDG）方法（[Johnson et al. 2010](https://doi.org/10.1021/ja100936w)）也在主功能 20 中（具体子功能编号以所用 Multiwfn 版本的菜单为准）。Multiwfn 作者推荐分析片段间相互作用时优先用 IGMH，因为它可以自定义片段、只显示片段间的作用区域，等值面也更平滑（见 [Sobereva 的 IGMH 介绍](http://sobereva.com/621)）。

### 9.1 准备波函数

CP 作业不产生分子轨道，所以要单独做一次带 `%chk` 的单点（不加 `Counterpoise`）：

```text
%nprocshared=8
%mem=8GB
%chk=water_dimer_wfn.chk
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ

Water dimer wavefunction for IGMH analysis

0 1
...dimer coordinates (no Fragment labels needed)...

```

IGMH 对计算级别要求不高，与几何优化相当的 DFT 级别（带极化函数的双 ζ 及以上基组）即可。然后：

```bash
formchk water_dimer_wfn.chk water_dimer_wfn.fchk
Multiwfn water_dimer_wfn.fchk
```

### 9.2 IGMH 操作要点

在 Multiwfn 中进入主功能 20 → 子功能 11（IGMH），按提示：

1. 定义片段：水二聚体中片段 1 为原子 1–3，片段 2 为原子 4–6；主客体中冠醚和 K⁺ 各为一个片段；
2. 选择格点质量（格点越细越平滑，但耗时按格点间距三次方增加）；
3. 导出 δg_inter 和 sign(λ2)ρ 的 cube 文件，用 VMD 绘制 sign(λ2)ρ 着色的 δg_inter 等值面。

着色约定：**蓝色**表示较强的吸引（如氢键、阳离子-氧配位），**绿色**为范德华（色散主导）作用，**红色**为位阻/排斥。等值面取值按体系调整（Sobereva 的示例中蛋白-配体体系用 0.002–0.006 a.u.）。完整的图文教程见 IGMH 官方教程包（链接在上述博文中）和 Multiwfn 手册；批量生成 fchk 与 cube 的脚本见 ../scripts/batch_cube_multiwfn.sh（`bash ../scripts/batch_cube_multiwfn.sh -h` 查看用法），Multiwfn 命令流批处理见 [第16章](16-postprocess-multiwfn.md)。

### 9.3 预期图像

| 体系 | IGMH 图中应看到 |
|---|---|
| 水二聚体 | O–H···O 之间一个小而偏蓝的圆盘状等值面 |
| 苯二聚体（平行错位） | 两环之间大面积、平行于环平面的绿色薄片（π-π 堆积特征） |
| 18-冠-6·K⁺ | K⁺ 与 6 个 O 之间的蓝色区域，冠醚 CH2 与 K⁺ 之间可能有绿色区域 |

---

## 参数选择与收敛建议

| 项目 | 推荐 | 原因 |
|---|---|---|
| 色散 | B3LYP/PBE0 加 GD3BJ；M06-2X 加 GD3；ωB97X-D、B97-D3、B2PLYPD3、DSDPBEP86 不再加 | 参数表决定可用组合，重复加会重复计算色散 |
| 基组 | 优化 ≥ def2-SVP；单点 ≥ def2-TZVP；阴离子加弥散 | BSSE 与基组不完备误差随基组增大而减小 |
| CP 校正 | 双 ζ、三 ζ 单点都建议做；报告时同时给出 raw 与 CP 值 | 便于读者判断 BSSE 大小 |
| CP 优化 | 双 ζ 基组研究弱络合物结构时考虑 | 小基组下 BSSE 会压缩分子间距离 |
| 形变能 | 柔性单体必须算；自由单体用最低能量构象 | 否则高估结合强度 |
| 积分格点 | 默认 UltraFine；M06 系列或极平坦势能面可用 `Int=SuperFine` | 弱相互作用能量差小，格点噪声相对更大 |
| 优化收敛 | 弱络合物势能面平坦，可用 `Opt=Tight`，并检查低频 | 避免残留虚频 |
| 溶液相 | 气相 CP + SMD 溶剂化能差 | `Counterpoise` 不能与 `SCRF` 合用 |
| 自由能 | 低频较多，用 qRRHO；加 1 atm→1 M 校正 | 弱络合物的熵贡献对低频敏感 |

---

## 常见问题与报错

**1. 碎片电荷多重度写错**

- 现象：作业在开头报错，或者碎片计算的电子数不对，得到荒谬的 BSSE。
- 原因：电荷多重度行只写了整体（如 `0 1`），或碎片电荷之和不等于整体电荷，或碎片多重度与其电子数奇偶不符（例如 `The combination of multiplicity 1 and 17 electrons is impossible.`）。
- 处理：按“整体、碎片 1、碎片 2……”逐对检查；每个碎片单独数电子。

**2. 漏标 `(Fragment=n)` 或编号不连续**

- 现象：报错，或某些原子被当作不属于任何碎片，BSSE 错误。
- 处理：每个原子都要有碎片标记；编号从 1 开始连续。GaussView 6 可用 Tools → Atom Groups（Gaussian Fragment 类）可视化地分配碎片，减少遗漏。

**3. `Counterpoise` 与 `SCRF` 或 `ONIOM` 同时使用时报错**

- 原因：官方文档明确说明 Counterpoise 不能与 ONIOM 或 SCRF 一起使用。
- 处理：气相做 CP，溶剂效应用单独的 SMD 单点处理（§6.4）。

**4. 后续步骤 `Guess=Read` 读不到轨道**

- 原因：CP 作业不产生分子轨道。
- 处理：需要波函数（IGMH、NBO、轨道图）时，另做一次普通单点并保存 chk（§9.1）。

**5. BSSE 为负或者非常大**

- BSSE 按定义应为正；为负通常说明碎片计算收敛到不同的电子态（尤其是开壳层或过渡金属碎片）。检查每个 `SCF Done` 和碎片的 ⟨S²⟩，必要时对碎片单独做 `Stable=Opt`（[第13章](13-openshell-stability.md)）。
- BSSE 很大（几 kcal/mol 以上）说明基组太小或缺少弥散函数，增大基组。

**6. CP 校正优化非常慢或不收敛**

- 每步 5 次 SCF，耗时本来就高；弱络合物势能面平坦，步数也多。
- 处理：先普通优化到接近收敛，再用其结构做 CP 优化；必要时 `Opt=(CalcFC,MaxCycles=200)`；或改用较大基组的普通优化 + CP 单点。

**7. 色散关键词报错 `R6DS8: Unable to choose the S8 parameter`**

- 原因：所用泛函没有对应的色散参数（例如 M06-2X 加 GD3BJ）。
- 处理：M06-2X 改用 GD3；查看官方 DFT 页面的 GD3/GD3BJ 参数表；ωB97X-D 等自带色散的泛函不加。

**8. 含赝势（ECP）元素的 CP 计算**

- Counterpoise 可以与 ECP 基组一起使用（官方示例使用 LANL2DZ）。若输出中的碎片电子数提示异常，先确认 GenECP 输入格式正确（[第3章](03-method-basis.md)），并检查每个碎片的 `SCF Done` 是否合理。

---

## 练习

1. 对水二聚体示例，计算 def2-SVP、def2-TZVP、aug-cc-pVTZ 三个基组下的 BSSE 和 ΔE_int(raw)、ΔE_int(CP)，作图观察两条曲线随基组的收敛。
2. 把苯二聚体第 2 步（普通优化）去掉 `EmpiricalDispersion=GD3BJ` 重新优化，观察两个苯环是否会远离；与带色散的结果比较。
3. 对冠醚·K⁺，分别用 B3LYP-D3(BJ) 和 ωB97X-D 算 ΔE_int(CP)，再算形变能，讨论两者的差别主要来自哪一项。
4. 用 §7 的鬼原子写法手动算出水二聚体的 BSSE，与 `BSSE energy` 行核对。
5. 对水二聚体或苯二聚体做 IGMH 分析，截图并说明等值面颜色与相互作用类型的对应关系。

---

## 参考资料

**官方文档**

- [Counterpoise keyword](https://gaussian.com/counterpoise/)：语法、`NewGhost/OldGhost`、可用作业类型与限制、输出示例
- [Molecule Specification](https://gaussian.com/molspec/)：`(Fragment=n)`、碎片电荷多重度行、鬼原子 `O-Bq`
- [Density Functional (DFT) Methods](https://gaussian.com/dft/)：`EmpiricalDispersion` 的 GD2/GD3/GD3BJ 参数表
- [MP methods / double hybrids](https://gaussian.com/mp/)：`DSDPBEP86`、`B2PLYPD3` 等双杂化泛函
- [Multiwfn 主页与手册](http://sobereva.com/multiwfn/)（2026 年起按发布日期命名版本）及 [在线手册](http://sobereva.com/multiwfn/manual/)

**相关 GitHub 仓库**

- [patonlab/pyDFTD3](https://github.com/patonlab/pyDFTD3)：对已有结构单独计算 D3(0)/D3(BJ) 色散能，可用来理解 `EmpiricalDispersion` 贡献的大小
- [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn)：中文文档的 Multiwfn 批处理工具，适合批量做 IGMH 等分析
- [polyluxus/runMultiwfn.bash](https://github.com/polyluxus/runMultiwfn.bash)：非交互运行 Multiwfn 的 bash 封装
- [1234zou/MOKIT](https://github.com/1234zou/MOKIT)：其 `frag_guess_wfn` 可生成碎片初猜及 LMO-EDA 等能量分解输入
- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian 输出，适合批量整理能量

**关键文献**

- S. F. Boys, F. Bernardi, *Mol. Phys.* 1970, 19, 553. [DOI: 10.1080/00268977000101561](https://doi.org/10.1080/00268977000101561)（Counterpoise 方法）
- S. Simon, M. Duran, J. J. Dannenberg, *J. Chem. Phys.* 1996, 105, 11024. [DOI: 10.1063/1.472902](https://doi.org/10.1063/1.472902)（CP 校正的势能面与优化）
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 2010, 132, 154104. [DOI: 10.1063/1.3382344](https://doi.org/10.1063/1.3382344)（DFT-D3）
- S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* 2011, 32, 1456. [DOI: 10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)（BJ 阻尼）
- P. Jurečka, J. Šponer, J. Černý, P. Hobza, *Phys. Chem. Chem. Phys.* 2006, 8, 1985. [DOI: 10.1039/B600027D](https://doi.org/10.1039/B600027D)（S22 基准集）
- J. Řezáč, K. E. Riley, P. Hobza, *J. Chem. Theory Comput.* 2011, 7, 2427. [DOI: 10.1021/ct2002946](https://doi.org/10.1021/ct2002946)（S66 基准集）
- E. R. Johnson et al., *J. Am. Chem. Soc.* 2010, 132, 6498. [DOI: 10.1021/ja100936w](https://doi.org/10.1021/ja100936w)（NCI 方法）
- T. Lu, Q. Chen, *J. Comput. Chem.* 2022, 43, 539. [DOI: 10.1002/jcc.26812](https://doi.org/10.1002/jcc.26812)（IGMH 方法）
- T. Lu, F. Chen, *J. Comput. Chem.* 2012, 33, 580. [DOI: 10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)（Multiwfn）

延伸阅读：后处理与出版级作图的通用流程见 ../../Common/tutorials/05-postprocess-visualization.md；从 log 中批量提取结构用 [../scripts/log2xyz.py](../scripts/log2xyz.py)，生成带碎片标记的 CP 输入可用 [../scripts/gjf_builder.py](../scripts/gjf_builder.py)（均可用 `-h` 查看用法）。

---

[← 上一章](13-openshell-stability.md) | [返回目录](../README.md) | [下一章 →](15-oniom.md)
