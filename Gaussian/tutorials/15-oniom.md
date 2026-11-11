# 第15章 ONIOM 多层方法与 QM/MM（ONIOM Multilayer Methods & QM/MM in Gaussian）

本章讲 Gaussian 16 的 ONIOM 多层方法：用 `ONIOM(B3LYP/6-31G(d):UFF)`、`ONIOM(B3LYP/def2SVP:PM6)`、三层 `ONIOM(High:Medium:Low)` 把大体系分层计算，以及用 Amber 力场和电子嵌入 `ONIOM=EmbedCharge` 做 QM/MM。内容包括原子行中的 MM 原子类型、电荷、冻结标记、层标记和链接原子（link atom）的写法，电荷多重度行的顺序，ONIOM 优化（微迭代、`Opt=QuadMac`）与频率，GaussView 设置方法，以及用 `Pop=SaveCM5` 为 QM 区原子提供电荷（ONIOM, QM/MM, link atom, electronic embedding）。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 理解 ONIOM 外推公式，知道 real、model、intermediate 体系分别是什么；
- 合理划分层边界，并为每个跨层共价键写出链接原子；
- 正确写出二层（3 对）和三层（6 对）的电荷多重度行；
- 写出带 Amber 原子类型、部分电荷、冻结标记和连接表的 QM/MM 输入；
- 理解机械嵌入（ME）与电子嵌入（EE）的区别，按“先 ME 后 EE”的策略优化；
- 读懂 `ONIOM: gridpoint ...` 与 `ONIOM: extrapolated energy` 输出，并用全 QM 参考检验分层是否合理。

| 适用场景 | 典型分层 | 说明 |
|---|---|---|
| 大分子中局部反应中心 | 反应中心 DFT，其余 PM6/UFF | 例：带大位阻取代基的催化剂 |
| 酶活性位点 | 底物 + 关键残基 DFT，蛋白其余部分 Amber | 经典 QM/MM，常用电子嵌入 |
| 显式溶剂中的溶质 | 溶质 DFT，周围水分子 Amber/TIP3P | 研究特定氢键或溶剂化光谱 |
| 固体/表面簇模型 | 活性中心 DFT，外围原子半经验或力场 | 外围原子常冻结 |
| 三层模型 | 高层 DFT，中层 PM6 或小基组 DFT，低层力场 | 兼顾精度与成本 |

> [!WARNING]
> ONIOM 不是“越大越好”的万能方案。对蛋白-小分子相互作用等问题，从蛋白中截取合理的簇模型做全 QM 计算往往更可靠、更容易收敛和解释；分层不当（边界离反应中心太近、切断极性键或共轭体系）会引入比低层方法本身更大的误差。建议在中文社区参考 Sobereva 关于“善用簇模型”的讨论（[sobereva.com/597](http://sobereva.com/597)），并始终用全 QM 或 S 值检验分层。

> [!NOTE]
> 本章示例目录：[../examples/15-oniom-hexanal/README.md](../examples/15-oniom-hexanal/README.md)（二层 MO:MM 与 MO:MO，含全 QM 参考）、[../examples/15-qmmm-phenol-water/README.md](../examples/15-qmmm-phenol-water/README.md)（Amber + 电子嵌入 QM/MM）、[../examples/15-oniom3-cyclohexanone/README.md](../examples/15-oniom3-cyclohexanone/README.md)（三层 ONIOM）。

---

## 2. 原理速览

### 2.1 ONIOM 外推

二层 ONIOM 把整个体系（real system, R）分成需要高精度描述的模型体系（model system, SM）和其余部分：

$$
E^{\rm ONIOM2} = E^{\rm low}(R) + E^{\rm high}(SM) - E^{\rm low}(SM)
$$

三层 ONIOM 再引入中间体系（intermediate model, IM）：

$$
E^{\rm ONIOM3} = E^{\rm low}(R) + E^{\rm medium}(IM) + E^{\rm high}(SM) - E^{\rm low}(IM) - E^{\rm medium}(SM)
$$

（公式见官方 ONIOM 技术说明 [gaussian.com/oniom_technote](https://gaussian.com/oniom_technote/)，原始方法见 [Dapprich et al. 1999](https://doi.org/10.1016/S0166-1280(98)00475-8)。）

可以把它理解为：用低层方法算整个体系，再把模型区域的低层描述“替换”成高层描述。只有当高、低层方法对模型体系之外的影响描述得差不多时，外推才准确——这正是分层合理性的判据。

### 2.2 链接原子

当共价键跨越层边界时，模型体系会出现悬挂键，需要用链接原子（通常是 H）封端。设低层原子 Y 与高层原子 X 成键，模型体系中 Y 被链接原子 H 取代，H 位于 X–Y 键轴上，距离按比例缩放：

$$
R_{\rm X-H} = g\,R_{\rm X-Y}
$$

缩放因子 g 由程序自动确定，也可以在原子行中显式给出。**链接原子必须由用户指定**，Gaussian 不会自动添加（GaussView 会自动生成，但仍需人工检查）。

### 2.3 机械嵌入与电子嵌入

| 方式 | 关键词 | QM–MM 静电作用 | QM 波函数是否被 MM 电荷极化 | 适用 |
|---|---|---|---|---|
| 机械嵌入（ME） | 默认（`NoEmbedCharge`） | 在 MM 层用点电荷计算 | 否 | 外围非极性、ONIOM(MO:MO)、初步优化 |
| 电子嵌入（EE） | `ONIOM(...)=EmbedCharge` | MM 点电荷进入 QM 哈密顿量 | 是 | 极性环境（酶、水溶液）中的 QM/MM |

电子嵌入只适用于含 MM 层的 MO:MM 和 MO:MO:MM 计算。默认的 `ScaleCharge=500`（即 555500）会关闭距 QM 区两个键以内的 MM 电荷、其余不缩放，避免边界附近的点电荷过度极化 QM 区。

### 2.4 分层的基本规则

1. 边界键尽量是非极性的 C–C 单键（sp³–sp³ 最好），不要切断双键、芳环、共轭体系或强极性键；
2. 边界离反应中心至少 2–3 个键；
3. 带电基团、参与反应的氢键网络、配位原子都应放进高层；
4. 层的划分在比较能量的所有结构（反应物、过渡态、产物）中必须完全相同；
5. 用全 QM 计算（小体系）或 S 值检验（`ONIOM=SValue`）评估外推误差。

---

## 3. ONIOM 输入格式详解

### 3.1 Route

```text
#p ONIOM(B3LYP/6-31G(d):UFF) Opt Freq
#p ONIOM(B3LYP/def2SVP:PM6) Opt Freq
#p ONIOM(B3LYP/6-31G(d):Amber)=EmbedCharge Opt
#p ONIOM(B3LYP/6-31G(d):PM6:UFF) Opt
```

以上四行依次是：二层 MO:MM（默认机械嵌入）、二层 MO:MO、二层 QM/MM 电子嵌入、三层 High:Medium:Low。

方法按 High:Medium:Low 的顺序写，用冒号分隔。常用选项（官方 ONIOM 页面）：

| 选项 | 作用 |
|---|---|
| `EmbedCharge` | 电子嵌入；默认 `NoEmbedCharge`（机械嵌入） |
| `ScaleCharge=ijklmn` | 边界附近 MM 电荷的缩放（各位数字 × 0.2）；默认 500；使用它即隐含 `EmbedCharge` |
| `MKUFF`（默认）/`MK`/`Mulliken`/`HLYGAt` | 电子嵌入的微迭代中，QM 区用哪种近似电荷 |
| `SValue` | 做完整的“方块”计算，给出 S 值检验 |
| `InputFiles` / `OnlyInputFiles` | 打印各个子计算的输入文件；后者只生成文件不计算 |
| `Compress`（默认） | 机械嵌入二阶导数计算时压缩到活动原子 |

ONIOM 支持能量、梯度和频率；各层还可以做 TD/CIS、NMR，`Gen`、`GenECP` 可用于相应的层。ONIOM **不能**与 `Counterpoise` 一起使用（[第14章](14-weak-bsse.md)）。ONIOM 与隐式溶剂的组合见 [第9章](09-solvation.md)。

### 3.2 原子行

官方给出的通用格式：

```text
atom [freeze-code] coordinate-spec layer [link-atom [bonded-to [scale-fac1 [scale-fac2 [scale-fac3]]]]]
```

| 字段 | 例子 | 说明 |
|---|---|---|
| atom | `C`、`C-CT`、`C-CT--0.18` | 元素[-MM 原子类型[-部分电荷]]；负电荷写成两个连字符，如 `O-OW--0.834` |
| freeze-code | `0`、`-1`、`-2` | 0 或省略：优化；−1：冻结；其他负整数：同值原子作为刚性块一起移动 |
| coordinate-spec | `1.23 0.45 -0.67` | 笛卡尔坐标（或 Z-matrix） |
| layer | `H`/`M`/`L` | High、Medium、Low |
| link-atom | `H`、`H-HC`、`H-HC-0.1` | 当本原子与更高层原子成键时，在较高层计算中用它替换本原子 |
| bonded-to | `7` | 链接原子连到哪个原子（较高层中的原子序号）；省略时程序尝试自动判断 |
| scale-fac | `0.709` | 链接原子键长缩放因子 g；一般省略，由程序确定 |

几个例子（前两行取自示例 15-1，第四行取自示例 15-2，第三行为说明格式用的虚构 Amber 原子）：

```text
C                  0    0.836684   -0.214571    0.306494  H
C                  0   -0.506895    0.217918   -0.281422  L H 7
C-CT--0.18         0    1.234000    0.000000    0.000000  L H-HC-0.1 5
O-OW--0.834000    -1    2.222605   -2.812808    0.140442  L
```

| 行 | 含义 |
|---|---|
| 1 | 高层碳原子，参与优化 |
| 2 | 低层碳原子，与高层原子 7 成键；在模型体系中被 H 取代 |
| 3 | Amber 原子：类型 CT、电荷 −0.18；在模型体系中被类型 HC、电荷 0.1 的 H 取代，该 H 连到原子 5 |
| 4 | 冻结（−1）的低层 TIP3P 水氧原子，类型 OW、电荷 −0.834 |

> [!IMPORTANT]
> 链接原子写在**较低层**那个原子的行上，`bonded-to` 填与它相连的**较高层**原子的序号。用 Amber 时链接原子也要带原子类型（如 `H-HC`）；用 UFF 时可以只写 `H`。

### 3.3 MM 原子类型、电荷与连接表

- **UFF、Dreiding**：程序可自动分配原子类型，可以不写类型（复杂体系仍建议显式写出）；
- **Amber**：**所有原子都必须显式给出 Amber 原子类型**（如 `C-CT`、`H-HC`、`O-OH`、`H-HO`、`C-CA`、`O-OW`、`H-HW`），电荷可选；
- 默认不给 MM 原子分配电荷，需要在原子行中写出，或用 `Amber=QEq` 等选项按 QEq 方法估算；
- MM 计算依赖原子间的连接关系，推荐用 `Geom=Connectivity` 给出连接表（GaussView 生成的输入默认就有）；
- 内置参数缺失时（例如 Amber 中没有某个键角参数）会报错，用 `Amber=HardFirst` 或 `SoftFirst` 在输入末尾补充参数，`$g16root/g16` 目录下有 `amber.prm`、`uff.prm`、`oplsaa.prm` 等参数文件可供参考。

### 3.4 电荷与多重度行

二层 ONIOM 最多写 4 对（第 4 对只用于 `SValue`）：

```text
charge(real,low) mult(real,low)  charge(model,high) mult(model,high)  charge(model,low) mult(model,low)
```

例如中性闭壳层体系写 `0 1 0 1 0 1`。只写一对时所有层都用它；写两对时第三对默认等于第二对。

三层 ONIOM 最多写 6 对（`SValue` 时 9 对），顺序为：

```text
real-low   int-medium   int-low   model-high   model-medium   model-low
```

例如 `0 1 0 1 0 1 0 1 0 1 0 1`。缺省对的默认规则见官方 ONIOM 页面中的表格。

> [!TIP]
> 模型体系的电荷 = 高层原子的电荷 + 链接原子（中性 H）。例如把带负电的羧酸根放在高层、把与之相连的碳链放在低层时，模型体系仍带 −1 电荷，第二、三对要写 `-1 1`；反之若带电基团在低层，模型体系是中性的，应写 `0 1`。每一对都必须与对应子体系的实际电子数一致（不过带电的反应中心通常不应放在低层）。

---

## 4. 实战一：二层 ONIOM(B3LYP/6-31G(d):UFF)——己醛

### 4.1 分层思路

以己醛 CH3(CH2)4CHO 为例：醛基及其相邻的两个 CH2 放在高层（B3LYP/6-31G(d)），其余丙基放在低层（UFF）。切断的是 C3–C4 这个 sp³–sp³ 单键，距羰基 3 个键。模型体系 = 高层原子 + 替换 C4 的链接 H，正好是丙醛（CH3CH2CHO）。

### 4.2 输入

文件 [../examples/15-oniom-hexanal/oniom2_b3lyp_uff_organic.gjf](../examples/15-oniom-hexanal/oniom2_b3lyp_uff_organic.gjf)（起始结构为 MMFF 预优化的全反式构象）：

```text
%nprocshared=8
%mem=8GB
%chk=hexanal_oniom_uff.chk
#p ONIOM(B3LYP/6-31G(d):UFF) Opt Freq

Hexanal 2-layer ONIOM: CHO-CH2-CH2 (high, B3LYP) + propyl tail (low, UFF)

0 1 0 1 0 1
O                  0    3.503965   -0.717638    1.086565  H
C                  0    3.339650    0.097274    0.185045  H
H                  0    4.207065    0.595327   -0.281611  H
C                  0    1.997438    0.503526   -0.374203  H
H                  0    2.007696    0.283038   -1.447059  H
H                  0    1.903302    1.587577   -0.248857  H
C                  0    0.836684   -0.214571    0.306494  H
H                  0    0.957244   -1.299380    0.194811  H
H                  0    0.853696   -0.005313    1.383406  H
C                  0   -0.506895    0.217918   -0.281422  L H 7
H                  0   -0.626859    1.302000   -0.165362  L
H                  0   -0.523033    0.004421   -1.357178  L
C                  0   -1.669989   -0.501791    0.400832  L
H                  0   -1.660920   -0.288341    1.476037  L
H                  0   -1.557099   -1.585932    0.284201  L
C                  0   -3.007957   -0.074340   -0.181099  L
H                  0   -3.166500    1.001042   -0.051661  L
H                  0   -3.825385   -0.601076    0.321170  L
H                  0   -3.062103   -0.303744   -1.250107  L

```

| 位置 | 内容 | 说明 |
|---|---|---|
| route | `ONIOM(B3LYP/6-31G(d):UFF) Opt Freq` | 高层 DFT、低层 UFF；默认机械嵌入 |
| 电荷多重度 | `0 1 0 1 0 1` | real-low、model-high、model-low 都是中性单重态 |
| 原子 1–9 | `H` | 醛基、α-CH2、β-CH2 |
| 原子 10 | `L H 7` | 低层碳 C4，在模型体系中被 H 取代，该 H 连到原子 7（β-碳） |
| 第二列 `0` | 冻结标记 | 全部原子参与优化 |
| 无原子类型 | UFF 自动分配 | 也可显式写 UFF 类型 |

### 4.3 运行与输出

```bash
cd Gaussian/examples/15-oniom-hexanal
g16 < oniom2_b3lyp_uff_organic.gjf > oniom2_b3lyp_uff_organic.log
grep "ONIOM: gridpoint" oniom2_b3lyp_uff_organic.log | tail -3
grep "ONIOM: extrapolated energy" oniom2_b3lyp_uff_organic.log | tail -1
grep -c "Optimization completed" oniom2_b3lyp_uff_organic.log
```

每个几何步都会打印三个子计算的能量和外推能量（示意输出，数值仅供参考；格式与官方技术说明一致）：

```text
 ONIOM: gridpoint  1 method:  low   system:  model energy:     0.00xxxxxxxxxx
 ONIOM: gridpoint  2 method:  high  system:  model energy:  -193.xxxxxxxxxxxx
 ONIOM: gridpoint  3 method:  low   system:  real  energy:     0.0xxxxxxxxxxx
 ONIOM: extrapolated energy =    -193.xxxxxxxxxxxx
```

| gridpoint | 计算 | 说明 |
|---|---|---|
| 1 | 低层（UFF）算模型体系 | MM 能量是相对力场零点的小数值，单独看没有绝对意义 |
| 2 | 高层（B3LYP）算模型体系（丙醛） | 一个普通的 DFT 计算，也有自己的 `SCF Done` |
| 3 | 低层（UFF）算整个体系 | — |
| extrapolated | $E_3 + E_2 - E_1$ | ONIOM 总能量，GaussView 显示的也是它；比较能量时用这一行 |

> [!NOTE]
> `SCF Done` 只是高层模型体系的能量，**不是** ONIOM 总能量。比较不同结构的能量、计算反应能时，一律使用 `ONIOM: extrapolated energy`（频率作业中的热力学量也基于 ONIOM 总能量）。

### 4.4 用全 QM 参考检验分层

同目录的 [../examples/15-oniom-hexanal/hexanal_full_qm_reference.gjf](../examples/15-oniom-hexanal/hexanal_full_qm_reference.gjf) 用 B3LYP/6-31G(d) 优化整个分子（原子顺序与 ONIOM 输入相同）。比较两者优化结构中醛基附近的键长、键角和二面角：差别应在千分之几埃和零点几度的量级；若差别大，说明边界离醛基太近或低层方法不合适。对反应问题，更好的检验是比较 ONIOM 与全 QM 的**相对能量**（反应能、能垒），而不是绝对能量。

---

## 5. 实战二：ONIOM(MO:MO)——B3LYP/def2-SVP:PM6

把低层换成半经验方法 PM6，就是 ONIOM(MO:MO)。文件 [../examples/15-oniom-hexanal/oniom2_b3lyp_pm6_momo.gjf](../examples/15-oniom-hexanal/oniom2_b3lyp_pm6_momo.gjf) 与 §4 的原子行完全相同，只改 route：

```text
%nprocshared=8
%mem=8GB
%chk=hexanal_oniom_pm6.chk
#p ONIOM(B3LYP/def2SVP:PM6) Opt Freq

Hexanal 2-layer ONIOM(MO:MO): high = B3LYP/def2-SVP, low = PM6

0 1 0 1 0 1
...same 19 atom lines as in section 4.2...

```

MO:MO 的特点：

- 低层也是量子化学方法，能描述低层中的电子效应（如共轭、电荷转移），通常比 MO:MM 外推更平滑；
- 三个子计算都有 SCF，若某个子计算的 SCF 难收敛，可用 `ONIOM=OnlyInputFiles` 打印各子计算输入，单独收敛后用 `Guess=Input` 为每个子计算指定初猜来源（每个子计算一行：`generate`、`read` 或 chk 文件名）；
- 电子嵌入不适用于 MO:MO（只有 MO:MM 与 MO:MO:MM 支持 `EmbedCharge`）；
- `ONIOM(...)=SValue` 会额外计算“高层算整个体系”等组合并打印 S 值，用于评估外推误差（只适合小体系）。

---

## 6. 实战三：QM/MM 电子嵌入——苯酚与显式水（Amber）

### 6.1 体系与分层

苯酚（13 个原子，高层 B3LYP/6-31G(d)）周围 6 个 TIP3P 水分子（低层 Amber）。QM 与 MM 之间没有共价键，因此不需要链接原子——这是“显式溶剂化”类 QM/MM 的典型情况。水分子 1、2 分别与苯酚 O–H 形成氢键（一个作受体、一个作供体），水 3 位于苯环上方；水 4–6 冻结，模拟外层环境。起始结构用 MMFF 预优化。

| 原子 | Amber 类型 | 电荷 | 层 | 冻结 |
|---|---|---|---|---|
| 苯酚 O | `OH` | 不写（0） | H | 0 |
| 与 O 相连的芳碳 | `C` | 不写 | H | 0 |
| 其余芳碳 | `CA` | 不写 | H | 0 |
| 羟基 H | `HO` | 不写 | H | 0 |
| 芳环 H | `HA` | 不写 | H | 0 |
| 水 O（6 个） | `OW` | −0.834 | L | 水 1–3：0；水 4–6：−1 |
| 水 H（12 个） | `HW` | +0.417 | L | 同上 |

苯酚的原子类型沿用 Amber 中酪氨酸侧链的类型；水用 TIP3P 电荷（O −0.834，H +0.417）。苯酚原子故意不写电荷：第 1 步用 `Pop=(Uncharged,SaveCM5)` 为它们拟合 CM5 电荷（§11）。

### 6.2 完整输入

文件 [../examples/15-qmmm-phenol-water/oniom2_qmmm_amber_embedcharge.gjf](../examples/15-qmmm-phenol-water/oniom2_qmmm_amber_embedcharge.gjf)：

```text
%nprocshared=8
%mem=8GB
%chk=phenol_qmmm.chk
#p ONIOM(B3LYP/6-31G(d):Amber) Geom=Connectivity Pop=(Uncharged,SaveCM5)

Phenol (QM) + 6 TIP3P waters (MM) step 1: ME single point, save QM CM5 charges

0 1 0 1 0 1
O-OH              0    1.664315   -0.239700   -0.069140  H
C-C               0    0.329439    0.044900   -0.114619  H
C-CA              0   -0.627481   -0.961724   -0.195016  H
C-CA              0   -1.980436   -0.621574   -0.238867  H
C-CA              0   -2.366411    0.719110   -0.204882  H
C-CA              0   -1.399673    1.721469   -0.123177  H
C-CA              0   -0.047410    1.383230   -0.074152  H
H-HO              0    1.787134   -1.224099   -0.097593  H
H-HA              0   -0.340157   -2.009129   -0.218556  H
H-HA              0   -2.734608   -1.402750   -0.297273  H
H-HA              0   -3.420236    0.985049   -0.242786  H
H-HA              0   -1.696405    2.768373   -0.110176  H
H-HA              0    0.703616    2.166467   -0.017991  H
O-OW--0.834000    0    2.222605   -2.812808    0.140442  L
H-HW-0.417000     0    1.830004   -3.136224    0.993001  L
H-HW-0.417000     0    2.840236   -3.542079   -0.043672  L
O-OW--0.834000    0    2.523336    0.598853    2.510257  L
H-HW-0.417000     0    2.358816    0.429870    1.556440  L
H-HW-0.417000     0    3.319708    1.155046    2.475818  L
O-OW--0.834000    0    0.547407   -1.042113    3.366029  L
H-HW-0.417000     0   -0.134834   -0.424518    3.044204  L
H-HW-0.417000     0    1.345102   -0.463298    3.271673  L
O-OW--0.834000   -1   -0.216169    1.191685   -3.403354  L
H-HW-0.417000    -1   -0.854421    0.705018   -2.848413  L
H-HW-0.417000    -1    0.630330    0.847465   -3.064842  L
O-OW--0.834000   -1    1.149070   -3.510847    2.503238  L
H-HW-0.417000    -1    0.899661   -2.614974    2.850515  L
H-HW-0.417000    -1    0.913612   -4.053819    3.275423  L
O-OW--0.834000   -1   -0.808643    3.872852   -3.127423  L
H-HW-0.417000    -1   -0.578612    2.919612   -3.197951  L
H-HW-0.417000    -1   -0.871480    4.119106   -4.062687  L

 1 2 1.0 8 1.0
 2 3 1.5 7 1.5
 3 4 1.5 9 1.0
 4 5 1.5 10 1.0
 5 6 1.5 11 1.0
 6 7 1.5 12 1.0
 7 13 1.0
 8
 9
 10
 11
 12
 13
 14 15 1.0 16 1.0
 15
 16
 17 18 1.0 19 1.0
 18
 19
 20 21 1.0 22 1.0
 21
 22
 23 24 1.0 25 1.0
 24
 25
 26 27 1.0 28 1.0
 27
 28
 29 30 1.0 31 1.0
 30
 31

--Link1--
%nprocshared=8
%mem=8GB
%chk=phenol_qmmm.chk
#p ONIOM(B3LYP/6-31G(d):Amber) Opt Geom=AllCheck Guess=Read

--Link1--
%nprocshared=8
%mem=8GB
%chk=phenol_qmmm.chk
#p ONIOM(B3LYP/6-31G(d):Amber)=EmbedCharge Opt Geom=AllCheck Guess=Read

```

| 步骤 | route 要点 | 目的 |
|---|---|---|
| 1 | `Geom=Connectivity Pop=(Uncharged,SaveCM5)` | 机械嵌入单点；保留输入中水的 TIP3P 电荷，只为未带电荷的 QM 原子（苯酚）拟合 CM5 电荷并存入 chk |
| 2 | `Opt Geom=AllCheck Guess=Read` | 机械嵌入优化；结构、层、原子类型、电荷、冻结信息和连接表都从 chk 读取 |
| 3 | `=EmbedCharge Opt Geom=AllCheck Guess=Read` | 在 ME 结构基础上做电子嵌入优化 |

连接表的格式：每行先写原子序号，再成对写“与之成键的、序号更大的原子 + 键级”；没有新键的原子也要单独占一行。芳环键级写 1.5。

> [!TIP]
> “先机械嵌入优化、再电子嵌入优化”是官方 ONIOM 页面推荐的策略：电子嵌入的每一步更贵，从较好的结构开始能省很多步。

### 6.3 运行与检查

```bash
cd Gaussian/examples/15-qmmm-phenol-water
g16 < oniom2_qmmm_amber_embedcharge.gjf > oniom2_qmmm_amber_embedcharge.log
grep -c "Normal termination" oniom2_qmmm_amber_embedcharge.log      # expect 3
grep "ONIOM: extrapolated energy" oniom2_qmmm_amber_embedcharge.log | tail -1
grep -i "undefined\|MM function not complete" oniom2_qmmm_amber_embedcharge.log
```

检查要点：

- 若出现 `Angle bend undefined between atoms ...` 和 `MM function not complete`，说明内置 Amber 参数缺少某个键角或键长项，按 §3.3 和“常见问题”用 `Amber=HardFirst` 补参数；
- 第 3 步的外推能量与第 2 步不能直接比较（嵌入方式不同，能量定义不同）；
- 比较第 2、3 步优化结构中苯酚 O–H···O(水) 的氢键长度，电子嵌入通常使氢键略有变化，体现 QM 区被环境极化。

### 6.4 电子嵌入的细节

- 电子嵌入时 MM 点电荷进入高层模型体系的哈密顿量；默认 `ScaleCharge=500` 关闭距 QM 区两个键以内的 MM 电荷，对有共价边界的体系很重要（本例无共价边界，不受影响）；
- 低层优化（微迭代）中 QM 区的静电用近似电荷表示，默认 `MKUFF`（基于 UFF 半径的 MK 拟合）；
- QM 区原子的 MM 电荷在电子嵌入的能量表达式中基本会抵消（real-low 与 model-low 中对应的项相减），但在机械嵌入和微迭代中有用，因此仍建议为它们提供合理电荷（§11）。

---

## 7. 实战四：三层 ONIOM(B3LYP/6-31G(d):PM6:UFF)

### 7.1 分层

以 4-叔丁基环己酮为例：

| 层 | 原子 | 加链接原子后得到的体系 |
|---|---|---|
| High（B3LYP/6-31G(d)） | C=O、两个 α-CH2 | 模型体系 = 丙酮 |
| Medium（PM6） | β-CH2 × 2、γ-CH | 中间体系 = 环己酮（叔丁基处加链接 H） |
| Low（UFF） | 叔丁基 | 真实体系 = 4-叔丁基环己酮 |

需要三个链接原子：两个 β-碳（中层）连到 α-碳（高层），叔丁基季碳（低层）连到 γ-碳（中层）。

### 7.2 输入

文件 [../examples/15-oniom3-cyclohexanone/oniom3_high_medium_low.gjf](../examples/15-oniom3-cyclohexanone/oniom3_high_medium_low.gjf)（MMFF 预优化的椅式、叔丁基在平伏键）：

```text
%nprocshared=8
%mem=8GB
%chk=tbu_chx_oniom3.chk
#p ONIOM(B3LYP/6-31G(d):PM6:UFF) Opt

4-tert-butylcyclohexanone 3-layer ONIOM: acetone / ring / tBu

0 1 0 1 0 1 0 1 0 1 0 1
O                  0    4.122036    0.136957   -0.633173  H
C                  0    3.002752   -0.071696   -0.168427  H
C                  0    2.244262   -1.341272   -0.464201  H
H                  0    2.313250   -1.993003    0.413434  H
H                  0    2.719285   -1.853964   -1.307291  H
C                  0    2.369477    0.866271    0.828466  H
H                  0    2.927042    1.808789    0.837501  H
H                  0    2.449750    0.413488    1.822599  H
C                  0    0.789301   -1.033011   -0.800459  M H 3
H                  0    0.269371   -1.986579   -0.932479  M
H                  0    0.755153   -0.508573   -1.762794  M
C                  0    0.912087    1.131695    0.467124  M H 6
H                  0    0.881395    1.717067   -0.459530  M
H                  0    0.481519    1.753582    1.257640  M
C                  0    0.103579   -0.175323    0.289372  M
H                  0    0.182603   -0.733238    1.234494  M
C                  0   -1.431990    0.063166    0.030840  L H 15
C                  0   -2.177817   -1.287418   -0.051602  L
H                  0   -1.906086   -1.851740   -0.949175  L
H                  0   -1.962500   -1.911264    0.822944  L
H                  0   -3.263383   -1.138444   -0.091803  L
C                  0   -1.702405    0.840723   -1.270832  L
H                  0   -2.778676    0.989002   -1.419802  L
H                  0   -1.241487    1.833214   -1.255093  L
H                  0   -1.328107    0.306071   -2.149341  L
C                  0   -2.056363    0.853811    1.202237  L
H                  0   -1.694080    1.885944    1.239497  L
H                  0   -3.147311    0.907906    1.106478  L
H                  0   -1.832658    0.377838    2.163373  L

```

| 位置 | 内容 | 说明 |
|---|---|---|
| 电荷多重度行 | 6 对 `0 1` | 依次为 real-low、int-medium、int-low、model-high、model-medium、model-low |
| 原子 9 | `M H 3` | 中层 β-碳；在模型体系（高层、中层都要计算）中被连到原子 3（α-碳）的 H 取代 |
| 原子 12 | `M H 6` | 另一个 β-碳，连到原子 6 |
| 原子 17 | `L H 15` | 叔丁基季碳；在中间体系（中层、低层都要计算）中被连到原子 15（γ-碳）的 H 取代，模型体系中不出现 |

### 7.3 输出

三层 ONIOM 有更多 gridpoint（低层算真实体系与中间体系、中层算中间体系与模型体系、高层算模型体系），外推能量仍在 `ONIOM: extrapolated energy` 行：

```bash
cd Gaussian/examples/15-oniom3-cyclohexanone
grep "ONIOM: gridpoint" oniom3_high_medium_low.log | tail -5
grep "ONIOM: extrapolated energy" oniom3_high_medium_low.log | tail -1
```

---

## 8. ONIOM 优化、过渡态、IRC 与频率

| 任务 | 推荐写法 | 说明（据官方 ONIOM 页面） |
|---|---|---|
| 极小点（MO:MM） | `Opt` | 默认：模型区域用内坐标 RFO 步，仅属于 MM 区的原子用线性微迭代；通常最佳 |
| 难收敛的极小点 | `Opt=QuadMac`，再不行 `Opt=(QuadMac,CalcFC)` 或 `Opt=(QuadMac,CalcAll)` | QuadMac 对所有原子做二次耦合步，考虑模型区与 MM 区的耦合 |
| 过渡态 | `Opt=(QuadMac,TS,CalcFC)` | 帮助沿正确方向找到过渡态 |
| 电子嵌入 | 先 ME 优化，再 EE 优化 | 比一开始就用 EE 更高效 |
| IRC | 先 `Freq` 确认过渡态，再 `IRC=RCFC`；大体系 `Freq=SaveNM` | MO:MM 默认积分器为 EulerPC，大体系最合适；中小体系可考虑 `IRC=HPC` |
| 微迭代次数 | `Opt=MaxMicroiterations=N` | 默认按原子数确定，至少 5000 |
| 冻结/解冻原子 | 原子行冻结标记，或 `Opt=ReadOpt` | 同一负值（< −1）的原子作为刚性块 |
| 频率 | `Freq`（与 `Opt` 同一 ONIOM 方案） | 任一层需要数值频率时，所有层都用数值频率 |

> [!WARNING]
> 冻结了大量原子时，被冻结的部分并不处在真正的驻点上，频率计算可能出现小虚频或异常低频，由此得到的热力学校正只能作为近似。比较能量时，所有结构必须使用**相同的分层、相同的冻结原子集合和相同的嵌入方式**。

---

## 9. 输出解读与常用 grep

| 想要的信息 | 命令 |
|---|---|
| 每个几何步的 ONIOM 总能量 | `grep "ONIOM: extrapolated energy" job.log` |
| 各子计算能量 | `grep "ONIOM: gridpoint" job.log` |
| 高层 SCF 是否收敛 | `grep "SCF Done" job.log` |
| 优化是否完成 | `grep -c "Optimization completed" job.log` |
| 虚频 | `grep "imaginary frequencies" job.log` |
| MM 参数缺失 | `grep -i "undefined\|not complete" job.log` |
| 作业是否成功 | `tail -1 job.log` 应为 `Normal termination` |

GaussView 打开 ONIOM 的 log 时，优化曲线显示的就是外推能量。用 [../scripts/opt_monitor.py](../scripts/opt_monitor.py) 可以把各步收敛判据画出来，[../scripts/g16_log_summary.py](../scripts/g16_log_summary.py) 汇总作业状态（均可用 `-h` 查看用法；它们针对普通作业设计，ONIOM 能量请以 `extrapolated energy` 行为准）。

---

## 10. 用 GaussView 设置 ONIOM；从 PDB/MD 结构生成输入

### 10.1 GaussView 6 的一般流程

（菜单名称以 GaussView 6 为例，小版本间可能略有差异。）

1. 建模或导入结构（PDB 文件会保留残基信息）；
2. 用层选择工具（Edit 菜单中的 Select Layer 一类功能）按原子、残基或距离选择高/中/低层；
3. 若低层用 Amber，在 MM 原子类型编辑功能中为每个原子指定类型和电荷，可从 PDB 残基自动分配，再人工检查非标准残基；
4. 在 Calculate → Gaussian Calculation Setup 的 Method 页启用多层 ONIOM 模型，分别设置各层方法、嵌入方式和各层电荷多重度；
5. 生成的输入会自动包含链接原子和 `Geom=Connectivity` 连接表。**务必检查链接原子是否只出现在你期望的边界键上**，删改原子后连接表也要同步修改。

### 10.2 从 PDB 或 MD 轨迹生成输入

- 蛋白体系先补全氢原子、确定质子化状态（如 PDB2PQR/PROPKA、H++ 等工具），再用 AmberTools（`tleap`、`antechamber`）为标准残基和配体指定原子类型与电荷（配体 RESP 电荷的做法见 [第12章](12-population-charges.md)）；
- 从 MD 轨迹中抽取帧时，通常以活性位点为中心保留一定半径内的残基和水，外层原子冻结；
- [jcerezochem/traj2oniom](https://github.com/jcerezochem/traj2oniom) 可从 MD 轨迹和拓扑生成带原子类型、层标记和冻结标记的 ONIOM 输入，可参考其输入格式；
- 外部程序作低层：ONIOM 也可以通过 `External` 关键词调用外部程序计算某一层（例如 xtb，见 [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian)）。

### 10.3 某个子计算 SCF 难收敛

官方推荐流程：

1. 在 route 中加 `ONIOM=OnlyInputFiles` 运行，得到各子计算的输入文件；
2. 单独运行有问题的那个（通常是高层模型体系），用 `Stable=Opt` 等手段（[第13章](13-openshell-stability.md)）得到正确的态，保存其 chk（如 `highmod.chk`）；
3. 原 ONIOM 作业加 `Guess=Input`，并在分子说明之后为每个子计算写一行初猜来源：MO:MM 只需一行（`highmod.chk`）；MO:MO 需要三行，依次为低层-真实体系、高层-模型体系、低层-模型体系，可写 `generate`、`read` 或 chk 文件名。

---

## 11. 为 QM 区原子提供 MM 电荷：Pop=SaveCM5 等

机械嵌入和微迭代需要 QM 区原子的 MM 电荷。G16 C.01 起可以把计算得到的原子电荷存入 chk，供后续 `Geom=Check`/`Geom=AllCheck` 的 MM 或 ONIOM 计算使用（官方 Population 页面）：

| 写法 | 作用 |
|---|---|
| `Pop=SaveCM5`、`Pop=SaveHirshfeld` | 保存 CM5 / Hirshfeld 电荷（CM5 与 Hirshfeld 电荷一起计算） |
| `Pop=SaveESP`、`Pop=SaveMulliken`、`Pop=SaveNPA` | 保存相应的电荷 |
| `Pop=(Uncharged,SaveCM5)` | 保留输入中已有的原子电荷，只为未带电荷的 QM 层原子拟合并保存 |

要点：

- 多层 ONIOM 中默认只保存显式计算的电荷，即 QM 层原子的电荷；不加 `Uncharged` 时，输入中原有的电荷会被新拟合的电荷替换；
- 后续步骤用 `Geom=Check` 或 `Geom=AllCheck` 读取，§6 的示例就是这样做的；
- CM5 电荷（[Marenich et al. 2012](https://doi.org/10.1021/ct200866d)）对基组依赖较小，适合给 QM 区原子提供合理的静电描述；若低层是 Amber 且要与蛋白电荷模型保持一致，也可以用 RESP 电荷（[第12章](12-population-charges.md)）。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 原因 |
|---|---|---|
| 边界键 | 非极性 sp³–sp³ C–C 单键，距反应中心 ≥ 2–3 个键 | 链接原子近似误差最小 |
| 高层大小 | 包含所有成键/断键原子、带电基团、关键氢键 | 外推误差主要来自边界附近 |
| 低层方法 | 有机体系 PM6 或 UFF；生物体系 Amber；需要电子效应时用 MO:MO | MM 不能描述低层中的电子结构变化 |
| 嵌入方式 | 极性环境用电子嵌入；先 ME 后 EE | 极化效应对能垒、光谱影响大 |
| MM 原子类型 | Amber 全部显式给出；UFF 复杂体系也建议显式 | 自动判断在非常规成键时不可靠 |
| 连接表 | `Geom=Connectivity`，删改原子后同步修改 | MM 项依赖连接关系 |
| 冻结 | 外层原子冻结，所有比较的结构冻结集合相同 | 保证能量可比 |
| 优化算法 | 默认 → `Opt=QuadMac` → `QuadMac,CalcFC` | 官方推荐顺序 |
| 过渡态 | `Opt=(QuadMac,TS,CalcFC)`，再 `Freq`、`IRC=RCFC` | 同第7章流程 |
| 验证 | 小体系做全 QM 对照，或 `ONIOM=SValue` | 量化分层误差 |

---

## 常见问题与报错

**1. `Angle bend undefined between atoms ...` / `MM function not complete`**

- 原因：内置力场参数表中没有某种原子类型组合的键角（或键长）参数。官方 MM 页面说明 Gaussian 只检查伸缩和弯曲项是否齐全，缺少的二面角项会被直接忽略。
- 处理：在 route 中写 `Amber=HardFirst`（ONIOM 中写在低层方法处，如 `ONIOM(B3LYP/6-31G(d):Amber=HardFirst)`），在输入末尾补充缺失参数，格式如官方 MM 页面示例中的 `HrmBnd1 HA CM HA 40.0 120.0`（函数名、三个原子类型、力常数、平衡键角）；参数值应来自文献或拟合，不能随意编造。也可以检查原子类型是否写错。

**2. Amber 计算中某些原子没有原子类型**

- 原因：Amber 要求所有原子显式指定类型（UFF/Dreiding 才会自动分配）。
- 处理：为每个原子写 `元素-类型`；GaussView 的 MM 类型编辑功能可以批量指定。

**3. 模型体系电子数与多重度不匹配（如 `The combination of multiplicity 1 and N electrons is impossible.`）**

- 原因：跨层共价键没有写链接原子，模型体系出现悬挂键；或带电/开壳层中心放错了层；或电荷多重度对数不够、顺序写错。
- 处理：检查每个跨层键都在较低层原子行写了链接原子；按 §3.4 写全各对电荷多重度。

**4. 链接原子连错了原子**

- 现象：高层计算的结构里出现奇怪的 H 位置，或能量异常。
- 处理：`bonded-to` 必须是与该低层原子直接成键的较高层原子序号；增删原子后序号会变，要同步修改。用 `ONIOM=OnlyInputFiles` 生成模型体系输入，在 GaussView 中检查。

**5. 高层（或某个子计算）SCF 不收敛**

- 处理：按 §10.3 用 `ONIOM=OnlyInputFiles` + 单独收敛 + `Guess=Input`。

**6. 微迭代或宏迭代不收敛，优化步数耗尽（L9999）**

- 处理：`Opt=QuadMac`；仍不行用 `Opt=(QuadMac,CalcFC)`；检查是否有 MM 原子因参数缺失而结构异常；必要时增加 `MaxCycles` 或 `MaxMicroiterations`。

**7. 电子嵌入时 QM 区被过度极化、能量或结构异常**

- 原因：带电 MM 原子离 QM 区太近（尤其是共价边界附近）。
- 处理：保持默认 `ScaleCharge` 或进一步缩放边界附近的电荷；把边界外移；把带电残基放进 QM 区。

**8. 想对 ONIOM 体系做 Counterpoise**

- `Counterpoise` 不能与 `ONIOM` 同时使用（官方 Counterpoise 页面）。需要 BSSE 校正时，在截取的全 QM 簇模型上做（[第14章](14-weak-bsse.md)）。

**9. 删改原子后作业读入出错**

- 原因：`Geom=Connectivity` 的连接表行数或序号与原子数不一致。
- 处理：重新生成连接表（GaussView 重新保存即可），确认每个原子都有一行。

---

## 练习

1. 运行己醛的三个输入（UFF、PM6、全 QM），比较 C=O 键长、C–C–C=O 二面角和优化步数，讨论哪种低层方法与全 QM 更接近。
2. 把己醛的层边界移到 C2–C3 之间（只有醛基和 α-CH2 在高层），重复练习 1，观察边界离反应中心更近时误差的变化。
3. 在苯酚-水 QM/MM 示例中，比较第 2 步（ME）和第 3 步（EE）优化结构中两个氢键的长度。
4. 对三层示例加 `ONIOM=OnlyInputFiles` 运行，查看生成的各子计算输入，确认模型体系是丙酮、中间体系是环己酮。
5. 用 GaussView 打开任一示例输入，查看层的着色和链接原子标记，尝试改变分层后重新保存，比较连接表的变化。

---

## 参考资料

**官方文档**

- [ONIOM keyword](https://gaussian.com/oniom/)：原子行语法、链接原子、电荷多重度对的默认规则、`EmbedCharge`/`ScaleCharge`/`SValue`/`OnlyInputFiles`、优化与 IRC 建议、`Guess=Input`
- [Investigating the Reactivity and Spectra of Large Molecules with ONIOM（技术说明）](https://gaussian.com/oniom_technote/)：外推公式与输出格式
- [Molecular Mechanics Methods](https://gaussian.com/mm/)：Amber/Dreiding/UFF、原子类型与电荷写法、QEq、`HardFirst`/`SoftFirst`/`SoftOnly`、参数文件
- [Opt keyword](https://gaussian.com/opt/)：`QuadMacro`、`MaxMicroiterations`、`ReadOpt`、微迭代选项
- [Population keyword](https://gaussian.com/population/)：`SaveCM5`、`SaveHirshfeld`、`SaveESP`、`Uncharged`
- [Guess keyword](https://gaussian.com/guess/)：`Guess=Input` 在 ONIOM 中的用法

**相关 GitHub 仓库**

- [jcerezochem/traj2oniom](https://github.com/jcerezochem/traj2oniom)：从 MD 轨迹和拓扑生成带原子类型、层标记和冻结标记的 ONIOM 输入
- [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian)：通过 `External` 让 Gaussian 调用 xtb，可作为便宜的低层或预优化手段
- [chauncey-garrett/gaussian-tools](https://github.com/chauncey-garrett/gaussian-tools)：其中有提取 ONIOM 外推能量的小脚本，可作为解析参考
- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian 输出
- [nanotech-empa/aiida-gaussian](https://github.com/nanotech-empa/aiida-gaussian)：在 AiiDA 中运行 Gaussian 并记录数据溯源，适合大批量 QM/MM 计算的管理

**关键文献**

- S. Dapprich, I. Komáromi, K. S. Byun, K. Morokuma, M. J. Frisch, *J. Mol. Struct. (THEOCHEM)* 1999, 461–462, 1. [DOI: 10.1016/S0166-1280(98)00475-8](https://doi.org/10.1016/S0166-1280(98)00475-8)（ONIOM 方法）
- T. Vreven, K. S. Byun, I. Komáromi, S. Dapprich, J. A. Montgomery Jr., K. Morokuma, M. J. Frisch, *J. Chem. Theory Comput.* 2006, 2, 815. [DOI: 10.1021/ct050289g](https://doi.org/10.1021/ct050289g)（ONIOM(QM:MM) 的扩展：电子嵌入、优化算法）
- L. W. Chung et al., *Chem. Rev.* 2015, 115, 5678. [DOI: 10.1021/cr5004419](https://doi.org/10.1021/cr5004419)（ONIOM 方法与应用综述）
- W. D. Cornell et al., *J. Am. Chem. Soc.* 1995, 117, 5179. [DOI: 10.1021/ja00124a002](https://doi.org/10.1021/ja00124a002)（Amber 力场）
- A. K. Rappé et al., *J. Am. Chem. Soc.* 1992, 114, 10024. [DOI: 10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040)（UFF）
- A. V. Marenich, S. V. Jerome, C. J. Cramer, D. G. Truhlar, *J. Chem. Theory Comput.* 2012, 8, 527. [DOI: 10.1021/ct200866d](https://doi.org/10.1021/ct200866d)（CM5 电荷）

延伸阅读：结构准备与建模（截取簇模型、补氢）见 ../../Common/tutorials/03-structure-modeling.md；作业提交与资源设置见 [第17章](17-hpc-jobs.md) 和 ../scripts/sub_g16.sh（`bash ../scripts/sub_g16.sh -h` 查看用法）；生成 ONIOM 输入骨架可参考 [../scripts/gjf_builder.py](../scripts/gjf_builder.py)（`python ../scripts/gjf_builder.py -h`）。

---

[← 上一章](14-weak-bsse.md) | [返回目录](../README.md) | [下一章 →](16-postprocess-multiwfn.md)
