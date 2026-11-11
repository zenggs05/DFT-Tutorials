# 第3章 方法、泛函与基组选择（含 Gen/GenECP 与赝势）（Choosing DFT Functionals, Dispersion Corrections, Basis Sets and ECPs in Gaussian 16）

本章讲 Gaussian 16 中"方法/基组"（model chemistry）怎么选、怎么写：Jacob 阶梯上的各类 DFT 泛函（B3LYP、PBE0、M06-2X、ωB97X-D、CAM-B3LYP、双杂化 DSD-PBEP86 等）及其 Gaussian 关键词，EmpiricalDispersion=GD3BJ/GD3 色散校正的正确搭配，Pople、Dunning 和 def2 基组，过渡金属与重元素的赝势（LANL2DZ、SDD、def2-ECP），以及用 Gen/GenECP/Pseudo=Read/ExtraBasis 自定义混合基组、从 Basis Set Exchange 下载基组的完整流程。最后介绍积分格点（Int=UltraFine）和 G09Defaults，帮助你写出既可靠又可复现的计算方案。

---

## 1. 本章目标与适用场景

| 你将学会 | 对应小节 |
|---|---|
| 理解 HF、DFT、MP2/双杂化的精度与成本差别 | §2 |
| 按 Jacob 阶梯选择泛函，并写出正确的 Gaussian 关键词 | §3 |
| 正确添加色散校正（GD3BJ、GD3），避免重复加或加错 | §4 |
| 读懂基组记号，选择 Pople、Dunning、def2 基组；知道什么时候需要弥散函数 | §5 |
| 用常见的 "high//low" 组合设计计算方案 | §6 |
| 为过渡金属和重元素选择赝势（ECP） | §7 |
| 写 Gen/GenECP 输入、从 BSE 下载并粘贴基组、用 Pseudo=Read 和 ExtraBasis | §8 |
| 理解积分格点和 G16/G09 默认值的差别 | §9 |
| 使用双杂化泛函，并理解冻芯设置 | §10 |

前置知识：输入文件结构见 [第2章](02-input-anatomy.md)。DFT 理论背景和学习路线见公共篇 第7章 DFT 理论速览、学习路线与资源汇总。本章示例在 [examples/03-method-basis](../examples/03-method-basis/README.md)。

> [!NOTE]
> 不存在"万能"的方法和基组。本章给出的是被大量基准测试支持的**常用起点**。正式研究中应查阅与你的体系最接近的基准研究，必要时自己用小体系做测试（与实验或高精度 CCSD(T) 数据比较）。

---

## 2. 原理速览

### 2.1 HF、DFT 与相关方法

Hartree–Fock（HF）把电子看成在其他电子的平均场中运动，完全忽略电子相关（交换项是精确的）。Kohn–Sham DFT 把能量写成电子密度的泛函：

$$
E_{\text{KS}}[\rho]=T_s[\rho]+\int v_{\text{ext}}(\mathbf r)\rho(\mathbf r)\,d\mathbf r+J[\rho]+E_{\text{xc}}[\rho]
$$

所有近似都集中在交换-相关泛函 $E_{\text{xc}}$ 上。杂化泛函把一部分 HF 精确交换混入：

$$
E_{\text{xc}}^{\text{hybrid}}=a_x E_x^{\text{HF}}+(1-a_x)E_x^{\text{DFA}}+E_c^{\text{DFA}}
$$

双杂化泛函再加入一部分 MP2 型的二阶微扰相关能：

$$
E_{\text{xc}}^{\text{DH}}=a_x E_x^{\text{HF}}+(1-a_x)E_x^{\text{DFA}}+(1-a_c)E_c^{\text{DFA}}+a_c E_c^{\text{PT2}}
$$

### 2.2 计算量随体系增长的规律

设基函数数为 $N$：

| 方法 | 形式标度 | 说明 |
|---|---|---|
| 纯泛函（GGA、meta-GGA） | $N^3$–$N^4$ | 可用密度拟合（DenFit）加速，是 DFT 中最便宜的 |
| 杂化泛函、HF | $N^4$（积分筛选后更低） | 精确交换使其比纯泛函慢 |
| MP2、双杂化 | $N^5$ | 内存和磁盘需求明显增加 |
| CCSD(T) | $N^7$ | 小分子基准用的"金标准" |

基组从双 ζ 换成三 ζ，基函数数大约翻一倍，杂化 DFT 的耗时增加数倍。这就是 §6 中"小基组优化 + 大基组单点"策略的来由。

### 2.3 基组是什么

分子轨道用原子中心的高斯函数线性展开：

$$
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r),\qquad
\chi_\mu=\sum_k d_{k\mu}\,N_k\,x^l y^m z^n e^{-\alpha_k r^2}
$$

一组 $(\alpha_k,d_{k\mu})$ 构成一个收缩函数。基组质量取决于：价层用几个 ζ（劈裂价层，split valence）、有没有极化函数（更高角动量，描述轨道形变）、有没有弥散函数（小指数，描述远离核的电子）。

---

## 3. 泛函：Jacob 阶梯与 Gaussian 关键词

### 3.1 Jacob 阶梯总表

Perdew 把泛函按所用信息的复杂程度排成"Jacob 阶梯"，一般阶梯越高越准确，也越昂贵。

| 阶梯 | 用到的信息 | 代表泛函 | Gaussian 关键词 | HF 交换比例 |
|---|---|---|---|---|
| LDA | $\rho$ | SVWN | `SVWN`（同义词 `LSDA`）、`SVWN5` | 0 |
| GGA | $\rho,\nabla\rho$ | BLYP、PBE、BP86、B97-D3 | `BLYP`、`PBEPBE`、`BP86`、`B97D3` | 0 |
| meta-GGA | $+\tau$ 或 $\nabla^2\rho$ | TPSS、M06-L、MN15-L | `TPSSTPSS`、`M06L`、`MN15L` | 0 |
| 杂化 GGA | + HF 交换 | B3LYP、PBE0、B3PW91、APFD | `B3LYP`、`PBE1PBE`、`B3PW91`、`APFD` | B3LYP 20%，PBE0 25% |
| 杂化 meta-GGA | meta-GGA + HF 交换 | M06-2X、MN15、TPSSh、PW6B95、BMK | `M062X`、`MN15`、`TPSSh`、`PW6B95`、`BMK` | M06-2X 54%，MN15 44%，TPSSh 10%，PW6B95 28%，BMK 42% |
| 范围分离杂化 | 短程/长程分别混合 HF 交换 | CAM-B3LYP、ωB97X-D、LC-ωHPBE、HSE06 | `CAM-B3LYP`、`wB97XD`、`LC-wHPBE`、`HSEH1PBE` | CAM-B3LYP 19%→65%，ωB97X-D 约 22%→100%，HSE06 短程 25%、长程 0 |
| 双杂化 | + MP2 型相关 | B2PLYP-D3、DSD-PBEP86、PBE0-DH、PBE-QIDH | `B2PLYPD3`、`DSDPBEP86`、`PBE0DH`、`PBEQIDH` | B2PLYP 53% HF + 27% PT2；PBE0-DH 50% + 12.5%；PBE-QIDH 约 69% + 33% |

关键词书写要点：

- Gaussian 中 **PBE0 的关键词是 `PBE1PBE`**，纯 PBE 是 `PBEPBE`（交换、相关各写一次），写成 `PBE0` 会报语法错误。
- **M06-2X 写作 `M062X`**，ωB97X-D 写作 `wB97XD`，HSE06 写作 `HSEH1PBE`，范围分离的 ωPBE 推荐写 `LC-wHPBE`。CAM-B3LYP 的关键词带连字符：`CAM-B3LYP`。
- 纯泛函名由"交换 + 相关"拼成：`B`+`LYP` = `BLYP`，`PBE`+`PBE` = `PBEPBE`。很多 GGA 都可以加前缀 `LC-` 得到长程校正版本，例如 `LC-BLYP`。
- 方法名前加 `U`（如 `UB3LYP`）表示非限制性计算，加 `RO` 表示限制性开壳层（[第13章](13-openshell-stability.md)）。
- Gaussian 的 B3LYP 使用 VWN **III** 型局域相关（官方 DFT 文档有说明），有些程序默认用 VWN V。跨程序比较 B3LYP 结果时，能量会有小差别。
- 列表之外的组合可以用 `IOp(3/76)`、`IOp(3/77)`、`IOp(3/78)` 自定义杂化比例（官方 DFT 页面有示例），但要谨慎，见 [第2章](02-input-anatomy.md) 对 IOp 的提醒。

### 3.2 各类泛函适合做什么

| 泛函 | 优点 | 局限 | 典型用途 |
|---|---|---|---|
| B3LYP-D3(BJ) | 文献数据极多，几何和频率较好，计算稳健 | 反应能垒偏低，色散必须另加 | 有机分子的优化与频率 |
| PBE0-D3(BJ) | 几何好，常用于过渡金属，激发能也不错 | 能垒仍偏低 | 金属配合物的优化，TD-DFT |
| M06-2X(-D3) | 主族热化学、能垒、非共价作用表现好 | 对格点敏感，**不适合过渡金属**（HF 交换比例高） | 有机反应能垒的单点 |
| ωB97X-D | 范围分离，内置色散，能垒、非共价作用和 CT 激发都较可靠 | 比 B3LYP 稍贵 | 单点、TS、激发态 |
| CAM-B3LYP | 修正长程行为，适合电荷转移激发 | 基态热化学一般 | TD-DFT（[第10章](10-tddft.md)） |
| MN15 | 主族和过渡金属的综合表现都好（原文献基准） | 对格点敏感；Gaussian 中没有为它定义 D3 参数，也不需要再加 | 金属体系的单点 |
| TPSSh | HF 交换比例低（10%） | — | 过渡金属自旋态能量的常用选择之一 |
| DSD-PBEP86 | 双杂化，GMTKN55 等基准中整体表现最好的一类 | 成本与 MP2 相当，只建议做单点 | 高精度单点 |

> [!TIP]
> 主族有机化学的一个稳健起点：B3LYP-D3(BJ)/def2-SVP（或 6-31G(d)）做优化和频率，再用 M06-2X-D3、ωB97X-D 或 DSD-PBEP86 加 def2-TZVP 做单点。过渡金属体系：PBE0-D3(BJ) 或 B3LYP-D3(BJ) 加 def2-SVP 做优化，加 def2-TZVP 做单点；自旋态问题要用多种泛函交叉验证。基准依据可参考 Goerigk 等的 GMTKN55 研究，以及 Bursch 等的"最佳实践"综述（见文末参考资料）。

---

## 4. 色散校正：EmpiricalDispersion

### 4.1 为什么必须加

常见的杂化泛函（B3LYP、PBE0 等）几乎描述不了长程色散（范德华）作用。这会导致：

- π–π 堆积、烷基链之间、主客体之间的结合被严重低估，甚至变成排斥；
- 构象能、大分子的折叠能、反应中的分子内弱作用出错；
- 优化出来的分子间距离偏大。

Grimme 的 D3 校正按原子对加一个经验的 $-C_6/r^6-C_8/r^8$ 项，几乎不增加计算量。**除非泛函本身已经包含色散（见 §4.3），否则正式研究中都应当加上。**

$$
E_{\text{disp}}^{\text{D3}}=-\sum_{A<B}\sum_{n=6,8}s_n\frac{C_n^{AB}}{r_{AB}^n}f_{\text{damp},n}(r_{AB})
$$

零阻尼（D3(0)）和 Becke–Johnson 阻尼（D3(BJ)）只是阻尼函数 $f_{\text{damp}}$ 的形式不同，参数（$s_6$、$s_8$、$s_{r,6}$ 或 $a_1$、$a_2$）对每个泛函单独拟合。**一个泛函只能使用已经为它拟合过参数的阻尼形式。**

### 4.2 Gaussian 中的选项

```text
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq
#p M062X/def2TZVP EmpiricalDispersion=GD3
```

| 选项 | 模型 | Gaussian 16 中定义了参数的泛函（摘自官方 DFT 文档） |
|---|---|---|
| `GD3BJ` | D3 + BJ 阻尼 | B3LYP、BLYP、PBE1PBE、PBEPBE、BP86、BPBE、B3PW91、TPSSTPSS、BMK、CAM-B3LYP、LC-wPBE；以及自带它的 B2PLYPD3、B97D3、PW6B95D3 |
| `GD3` | D3 + 零阻尼 | B3LYP、BLYP、PBE1PBE、PBEPBE、BP86、BPBE、B3PW91、TPSSTPSS、BMK、CAM-B3LYP、LC-wPBE、**M05、M052X、M06L、M06、M062X、M06HF**、B97D3、PW6B95D3、B2PLYPD3 |
| `GD2` | 较旧的 D2 | B3LYP、BLYP、BP86、PBEPBE、TPSSTPSS；以及自带它的 B97D、B2PLYPD、mPW2PLYPD |
| `PFD` | Petersson–Frisch 色散 | APFD 泛函自带 |

几点说明：

- **M06 系列只能用 `GD3`**。它们的 GD3 参数中 $s_8=0$，即只加 $C_6$ 项。
- MN15、MN15L、M08HX、M11、TPSSh、HSEH1PBE 等不在上表中，Gaussian 没有为它们定义 D3 参数，不能直接加 `EmpiricalDispersion`（MN15 系列本身也不建议再加）。
- 文献中常写作 B3LYP-D3(BJ)、M06-2X-D3(0)、PBE0-D3(BJ)。在论文方法部分写清楚阻尼类型。

### 4.3 已经自带色散的泛函

| 关键词 | 自带的色散 | 能否再加 `EmpiricalDispersion` |
|---|---|---|
| `wB97XD` | Grimme D2 的一个变体 | **不能** |
| `B97D3` | D3(BJ) | 不能 |
| `B97D` | D2 | 不能 |
| `B2PLYPD3` | D3(BJ) | 不能 |
| `DSDPBEP86` | D3(BJ) | 不能 |
| `PW6B95D3` | D3(BJ) | 不能 |
| `APFD` | PFD | 不能 |

> [!WARNING]
> 给不支持的泛函加 `GD3BJ`（最常见的是 `M062X ... EmpiricalDispersion=GD3BJ`）会在 L301 报 `R6DS8: Unable to choose the S8 parameter`。给已经自带色散的泛函再加一次色散，有时不会报错，但色散被重复计算，结果是错的。

### 4.4 在输出中查看色散能

色散能会加到核排斥能中一起打印，大致形如（示意输出，数值仅供参考）：

```text
 R6Disp:  Grimme-D3(BJ) Dispersion energy=       -0.0082345678 Hartrees.
 Nuclear repulsion after empirical dispersion term =      203.0012345678 Hartrees.
```

```bash
grep -m1 "Dispersion energy" job.log
```

因此 `SCF Done` 的能量**已经包含色散校正**，不需要再手动相加。

---

## 5. 基组

### 5.1 记号速读

| 记号 | 含义 |
|---|---|
| STO-3G | 最小基：每个原子轨道一个收缩函数，只用于测试 |
| 3-21G、6-31G | 劈裂价层（双 ζ）：芯层 1 个收缩函数（6 个初始高斯），价层拆成 2 个（3+1） |
| 6-311G | 价层三 ζ（3+1+1） |
| `(d)` 或 `*` | 重原子（非 H）加一组 d 极化函数 |
| `(d,p)` 或 `**` | 重原子加 d，H 再加 p 极化函数 |
| `(2df,2pd)` 等 | 多组、更高角动量的极化函数 |
| `+` / `++` | 重原子加弥散 sp 函数 / 再给 H 加弥散 s 函数 |
| cc-pVnZ（n = D、T、Q、5） | Dunning 相关一致基组，按 ζ 数系统增大，便于外推到完备基组极限 |
| aug-cc-pVnZ | 每种角动量各加一组弥散函数；jun-、may- 等"月份"变体逐步删减弥散函数 |
| def2-SVP、def2-TZVP、def2-TZVPP、def2-QZVP | Karlsruhe 基组：全周期表设计一致，从 Rb 起自动配 def2 ECP |
| def2-SV(P) | def2-SVP 去掉 H 上的 p 函数，Gaussian 关键词为 `Def2SVPP` |

### 5.2 Gaussian 关键词与覆盖元素

Gaussian 中 **def2 基组的关键词没有连字符**：`Def2SVP`、`Def2TZVP`、`Def2TZVPP`、`Def2QZVP`（大小写不限）。

| 基组 | Gaussian 关键词 | 覆盖元素（官方 Basis Sets 文档） | d 函数默认类型 |
|---|---|---|---|
| STO-3G | `STO-3G` | H–Xe | — |
| 3-21G | `3-21G` | H–Xe | 6D |
| 6-31G 系列 | `6-31G(d)`、`6-31+G(d,p)` 等 | **H–Kr** | **6D（笛卡尔）** |
| 6-311G 系列 | `6-311G(d,p)`、`6-311++G(2d,2p)` 等 | H–Kr | 5D（球谐） |
| cc-pVnZ | `cc-pVDZ`、`cc-pVTZ`、`aug-cc-pVTZ` 等 | H–Ar、Ca–Kr（随 n 有差别） | 5D |
| def2 系列 | `Def2SVP`、`Def2TZVP`、`Def2TZVPP`、`Def2QZVP` | H–La、Hf–Rn（Rb 起自动带 ECP） | 5D |
| LANL2DZ | `LANL2DZ` | H、Li–La、Hf–Bi（**Na 起即用 ECP**） | 见注 |
| SDD | `SDD` | 除 Fr、Ra 外全部元素（Ar 以前用 D95 全电子，之后用 Stuttgart/Dresden ECP） | 见注 |

注：官方文档说明，内置基组中 3-21G、6-21G、4-31G、6-31G（含 6-31G†/‡）、CEP-31G、D95、D95V 默认使用笛卡尔 d 函数（6D），其余大多为球谐 d（5D），f 及更高角动量函数一律为球谐（7F）。同一个作业中所有 d 壳层必须是同一种类型。可以用 `5D`/`6D`、`7F`/`10F` 关键词覆盖默认值。

> [!IMPORTANT]
> Pople 基组只覆盖到 Kr。分子中含 I、Pd、Pt、Au 等元素时，写 `6-31G(d)` 会在 L301 报 `Atomic number out of range for 6-31G basis set.`。解决办法：全体系改用 def2 系列（自动带 ECP），或用 GenECP 给重元素单独指定 LANL2DZ/SDD/def2 基组和 ECP（§7、§8）。

### 5.3 基组大小与基函数数目

以苯（C6H6）为例，按各基组的标准收缩方案计算基函数数目（实际数目以输出中的 `basis functions` 行为准）：

| 基组 | C 的收缩 | H 的收缩 | 苯的基函数数 |
|---|---|---|---|
| 6-31G(d) | [3s2p1d]，6D 共 15 | [2s] = 2 | 102 |
| 6-31G(d,p) | 15 | [2s1p] = 5 | 120 |
| 6-311+G(d,p) | [5s4p1d] = 22 | [3s1p] = 6 | 168 |
| def2-SVP | [3s2p1d] = 14 | [2s1p] = 5 | 114 |
| def2-TZVP | [5s3p2d1f] = 31 | [3s1p] = 6 | 222 |
| cc-pVTZ | [4s3p2d1f] = 30 | [3s2p1d] = 14 | 264 |
| def2-TZVPP | 31 | [3s2p1d] = 14 | 270 |
| aug-cc-pVTZ | [5s4p3d2f] = 46 | [4s3p2d] = 23 | 414 |
| def2-QZVP | [7s4p3d2f1g] = 57 | [4s3p2d1f] = 30 | 522 |

```bash
grep -m1 "basis functions" job.log
# benzene/def2-SVP (illustrative):  114 basis functions, ... primitive gaussians, 120 cartesian basis functions
```

### 5.4 什么时候需要弥散函数

| 情形 | 是否需要弥散函数 | 推荐 |
|---|---|---|
| 阴离子（尤其是小阴离子，如 F⁻、OH⁻）、电子亲和能 | **需要** | 6-31+G(d)、6-311+G(d,p)、aug-cc-pVTZ、def2-TZVPD（需从 BSE 获取） |
| 氢键、弱相互作用的结合能 | 建议（三 ζ 以上时影响变小） | def2-TZVP(D)、aug-cc-pVTZ，并做 BSSE 校正（[第14章](14-weak-bsse.md)） |
| Rydberg 激发、极化率、超极化率 | 需要 | aug-cc-pVTZ、6-311++G(2d,2p) |
| 中性分子的几何优化和频率 | 一般不需要 | def2-SVP、6-31G(d) |
| pKa、溶剂化自由能 | 视情况 | 参照所用模型（如 SMD）原文的建议（[第9章](09-solvation.md)） |

弥散函数会带来近线性相关问题：Gaussian 会剔除重叠矩阵中特征值很小的函数组合，输出中 `NBsUse` 小于 `NBasis` 就是这种情况。偶尔还会引起 SCF 收敛困难（[第4章](04-sp-scf.md)）。

> [!NOTE]
> Gaussian 内置的 def2 关键词只有 `Def2SV`、`Def2SVP`、`Def2SVPP`、`Def2TZV`、`Def2TZVP`、`Def2TZVPP`、`Def2QZV`、`Def2QZVP`、`Def2QZVPP` 和 `QZVP`。带弥散的 def2-SVPD、def2-TZVPD，以及 ma-def2、cc-pVnZ-PP、x2c 等基组都**不是内置的**，要从 Basis Set Exchange 下载后用 Gen 读入（§8.3）。

---

## 6. 常用组合与 high//low 记号

### 6.1 high//low

"单点水平//优化水平"记号表示：用 low 水平优化结构（并在同一水平算频率和热校正），再在该结构上用 high 水平算单点能。例如：

$$
\text{DSD-PBEP86/def2-TZVP//B3LYP-D3(BJ)/def2-SVP}
$$

最终的自由能通常这样组合：

$$
G = E_{\text{elec}}^{\text{high}} + G_{\text{corr}}^{\text{low}}
$$

其中 $G_{\text{corr}}$ 是低水平频率计算给出的 `Thermal correction to Gibbs Free Energy`（[第6章](06-freq-thermo.md)）。这样做的依据是：几何和振动频率对方法、基组的敏感度远低于电子能量。

### 6.2 推荐组合速查

| 任务 | 优化 + 频率 | 单点 | 备注 |
|---|---|---|---|
| 有机分子构象、反应机理 | B3LYP-D3(BJ)/def2-SVP 或 6-31G(d) | M06-2X-D3(0)/def2-TZVP、ωB97X-D/def2-TZVP、DSD-PBEP86/def2-TZVP（或 def2-QZVP） | 溶液体系的优化和单点都要加 SMD/PCM（[第9章](09-solvation.md)） |
| 过渡金属催化 | PBE0-D3(BJ) 或 B3LYP-D3(BJ)/def2-SVP | 同一泛函或多种泛函 + def2-TZVP | def2 自带 ECP，无需 GenECP；自旋态要交叉验证 |
| 非共价作用、主客体 | B3LYP-D3(BJ)/def2-SVP 或 def2-TZVP | ωB97X-D 或 B3LYP-D3(BJ)/def2-TZVPP(D)、DSD-PBEP86 | 计算 BSSE（[第14章](14-weak-bsse.md)） |
| 阴离子、电子亲和能 | B3LYP-D3(BJ)/6-31+G(d) | ωB97X-D/aug-cc-pVTZ 或 def2-TZVPD | 弥散函数必不可少 |
| 激发态、UV-Vis | 基态用 PBE0 或 B3LYP-D3(BJ)/def2-SVP | TD-CAM-B3LYP、TD-ωB97X-D、TD-PBE0 加 def2-TZVP | [第10章](10-tddft.md) |
| NMR | B3LYP-D3(BJ)/6-31G(d) | mPW1PW91 或 B3LYP 加 6-311+G(2d,p)，GIAO | [第11章](11-nmr-props.md) |
| 快速筛选 | 半经验 PM7 或 GFN2-xTB（外部程序） | B3LYP-D3(BJ)/def2-SVP | 只作预优化 |

### 6.3 示例：B3LYP-D3(BJ)/def2-SVP 优化 + 频率

[`benzene_b3lyp_gd3bj_def2svp.gjf`](../examples/03-method-basis/benzene_b3lyp_gd3bj_def2svp.gjf)：

```text
%chk=benzene_b3lyp_gd3bj_def2svp.chk
%mem=4GB
%nprocshared=8
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt Freq

benzene opt and freq, B3LYP with D3BJ dispersion, def2SVP

0 1
C     0.000000    1.395248    0.000000
C     1.208320    0.697624    0.000000
C     1.208320   -0.697624    0.000000
C     0.000000   -1.395248    0.000000
C    -1.208320   -0.697624    0.000000
C    -1.208320    0.697624    0.000000
H     0.000000    2.482360    0.000000
H     2.149787    1.241180    0.000000
H     2.149787   -1.241180    0.000000
H     0.000000   -2.482360    0.000000
H    -2.149787   -1.241180    0.000000
H    -2.149787    1.241180    0.000000

```

| 关键词 | 说明 |
|---|---|
| `B3LYP/def2SVP` | 杂化泛函 + 双 ζ 极化基组；def2-SVP 对全周期表一致，便于以后加入金属 |
| `EmpiricalDispersion=GD3BJ` | D3(BJ) 色散校正，B3LYP 有对应参数 |
| `Opt Freq` | 优化后在同一水平做频率，确认是极小点并得到热校正 |

### 6.4 示例：M06-2X 搭配 GD3 做单点

[`ethylene_m062x_gd3_def2tzvp_sp.gjf`](../examples/03-method-basis/ethylene_m062x_gd3_def2tzvp_sp.gjf)：

```text
%chk=ethylene_m062x_gd3_def2tzvp_sp.chk
%mem=4GB
%nprocshared=8
#p M062X/def2TZVP EmpiricalDispersion=GD3 SP

ethylene single point, M062X with zero damping D3, def2TZVP

0 1
C     0.000000    0.000000    0.667480
...
```

把 `GD3` 改成 `GD3BJ` 重新运行，就能亲眼看到 L301 的 `R6DS8` 报错。

### 6.5 示例：双杂化单点（high//low）

[`dsdpbep86_def2tzvp_sp.gjf`](../examples/03-method-basis/dsdpbep86_def2tzvp_sp.gjf) 直接从苯的优化 chk 读结构：

```text
%oldchk=benzene_b3lyp_gd3bj_def2svp.chk
%chk=benzene_dsdpbep86_def2tzvp_sp.chk
%mem=8GB
%nprocshared=8
#p DSDPBEP86/def2TZVP Geom=AllCheck SP

```

`%oldchk` 读、`%chk` 写，`Geom=AllCheck` 从旧 chk 读取结构、电荷多重度和 title，所以输入里只有 route（[第2章](02-input-anatomy.md)）。`DSDPBEP86` 已自带 D3(BJ)，不要再加 `EmpiricalDispersion`。最终能量的读法见 §10.2。

---

## 7. 重元素与过渡金属：赝势（ECP）

### 7.1 为什么用 ECP

对于第四周期以后的元素：

1. 内层电子多，全电子基组很大，计算昂贵；
2. 内层电子速度接近光速，**标量相对论效应**显著（例如 Au 的 6s 轨道收缩），非相对论的全电子计算会出现系统误差。

有效核势（ECP，effective core potential，也称赝势）用一个解析势代替内层电子，只显式处理价层电子。主流 ECP 是通过拟合相对论原子计算得到的，因此**隐式包含了标量相对论效应**。

| ECP / 基组 | 芯层大小 | 特点 | Gaussian 中的写法 |
|---|---|---|---|
| LANL2DZ（Hay–Wadt） | 过渡金属为小芯（保留次外层 (n−1)s、(n−1)p 电子），主族元素为大芯；具体替换的电子数随元素变化 | 老牌，价层双 ζ，没有极化函数；**Na–Ar 也用 ECP** | 关键词 `LANL2DZ`；Gen 中写 `LANL2DZ` |
| SDD（Stuttgart/Dresden） | 多为小芯，如第二过渡系 28 电子芯（MWB28） | 质量较好；Ar 以前用 D95 全电子 | 关键词 `SDD`；Gen 中可写 `SDD`，或 `MWB28` 这类带芯电子数的名称 |
| def2-ECP | Rb–Xe 28 电子，Cs–La 46 电子，Hf–Rn 60 电子 | 与 def2 基组配套 | 用 def2 关键词时自动加上；GenECP 的 ECP 段可写 `def2` |

> [!NOTE]
> 每个元素被 ECP 替换了多少电子，可以在输出中核对：使用 ECP 后，`alpha electrons`/`beta electrons` 行给出的只是价层电子数，比全电子计算少了芯电子数。也可以在 BSE 网站上查询。不确定时不要凭记忆写进论文。

### 7.2 怎么选

| 情形 | 推荐 |
|---|---|
| 新项目，含过渡金属或重主族元素 | **全体系 def2-SVP（优化）/def2-TZVP（单点）**，内置关键词，自动配 def2-ECP，最省事也最一致 |
| 需要与大量文献直接对比 | 跟随文献，例如"LANL2DZ（金属）+ 6-31G(d)（其他）"，用 GenECP |
| 想给金属加极化函数 | LANL2DZ + f 函数（如 LANL2DZ(f)），或 LANL08(f)、def2 系列（BSE 获取） |
| 需要全电子 + 相对论 | `Int=DKH2`（二阶 Douglas–Kroll–Hess 标量相对论）配合为其设计的全电子基组（如 SARC、x2c 系列，从 BSE 获取） |

> [!WARNING]
> 直接对整个分子写 `B3LYP/LANL2DZ`，C、H、N、O 得到的是 D95V 基组（没有极化函数），P、S、Cl 等**第三周期元素也会被换成 ECP**。这通常不是你想要的。需要"金属用 LANL2DZ、其他用 Pople"时，请用 GenECP（§8.2）。

---

## 8. Gen、GenECP、Pseudo=Read 与 ExtraBasis

### 8.1 Gen 输入格式

route 中把基组写成 `Gen`，基组就从坐标后的附加段读入。每种元素（或原子编号）一个"中心块"：

```text
C H O 0          <- center line: element symbols and/or atom numbers, terminated by 0
6-31G(d)         <- a built-in basis name can replace explicit shells
****             <- four asterisks end the block
N 0
S    3   1.00    <- explicit shell: type, number of primitives, scale factor
  ...exponent  coefficient...
****
                 <- blank line ends the whole basis section
```

规则（摘自官方 Gen 文档）：

- 中心行列出元素符号和/或原子序号，以 `0` 结尾，可以混写，例如 `C H 0`、`1 5 0`。
- 每个中心块以 `****` 结束。
- 块内可以直接写内置基组名（如 `6-31G(d)`、`Def2TZVP`、`LANL2DZ`、`SDD`），也可以写显式的壳层：`类型 原始高斯数 标度因子`，后面每行一个"指数 系数"（SP 壳层为"指数 s 系数 p 系数"）。数字可以用 `D` 或 `E` 表示指数部分。
- 同一元素可以出现在多个块中，函数会累加。官方例子就是先给 `C H 0` 指定 `6-31G(d,p)`，再给 `1 0`（第 1 号原子）单独加一个弥散 SP 壳层。
- **分子中每一种元素都必须有基组**；中心行中**不能出现分子里没有的元素**，否则作业中止。写基组库文件时可以用 `-H` 这种带负号的写法，表示"有就用，没有就跳过"。
- Gen 输入默认 **5D、7F**（与 route 中写的方法无关）。要与内置 `6-31G(d)`（6D）完全一致，需要在 route 中加 `6D`。
- 可以用 `@文件名/N` 引入外部的 .gbs 文件；`$G16BASIS` 目录下有 Gaussian 自带的通用基组文件（官方文档提到 `631.gbs`、`631s.gbs` 等）。

### 8.2 GenECP：混合基组 + 赝势

只要有任何原子使用 ECP，就要写 **`GenECP`**（它等价于 `Gen Pseudo=Read`）。输入在基组段之后、隔一个空行，再写 ECP 段：

[`pd_complex_genecp_lanl2dz.gjf`](../examples/03-method-basis/pd_complex_genecp_lanl2dz.gjf)（trans-PdCl2(PMe3)2，坐标略）：

```text
%chk=pdcl2_pme3_2_genecp.chk
%mem=8GB
%nprocshared=8
#p B3LYP/GenECP 6D EmpiricalDispersion=GD3BJ Opt Freq

trans PdCl2(PMe3)2, 6-31G(d) on C H P Cl, LANL2DZ basis and ECP on Pd

0 1
 Pd      0.000000      0.000000      0.000000
 Cl      2.310000      0.000000      0.000000
 ... (29 atoms in total) ...

C H P Cl 0
6-31G(d)
****
Pd 0
LANL2DZ
****

Pd 0
LANL2DZ

```

| 部分 | 说明 |
|---|---|
| `GenECP` | 基组段 + ECP 段都从输入读 |
| `6D` | 让 6-31G(d) 与内置定义一致（笛卡尔 d）。Gen 默认 5D；加不加都可以，但**同一项目中要比较的计算必须统一** |
| `C H P Cl 0` / `6-31G(d)` / `****` | 轻原子的基组 |
| `Pd 0` / `LANL2DZ` / `****` | Pd 的价层基组 |
| 空行 | 结束基组段 |
| `Pd 0` / `LANL2DZ` | ECP 段：使用 LANL2DZ 基组配套的 ECP（简化写法：直接写内置基组名） |
| 末尾空行 | 结束 ECP 段，也是文件末尾空行 |

ECP 段的简化写法也可以用 Stuttgart 名称，例如 `Pd 0` 下一行写 `SDD`（同时基组段也写 `SDD`），或用带芯电子数的 `MWB28`；用 def2 基组时 ECP 段可以写 `def2`（官方文档说明 `def2` 和同义词 `QZV` 表示 def2/QZV 基组所用的 ECP）。

> [!CAUTION]
> 只写 `Gen`、不读 ECP，却给 Pd 用了 LANL2DZ 价层基组：程序会对 Pd 做"全电子"计算，但基组根本描述不了内层电子，结果是毫无意义的能量或 SCF 不收敛。这种错误**不一定报错**，所以每次都要检查输出中是否确实读入了 ECP：最可靠的办法是核对电子数（`grep -m1 "alpha electrons" job.log`）。trans-PdCl2(PMe3)2 全电子共 164 个电子，LANL2DZ 的 ECP 替换掉 Pd 的 28 个芯电子后应剩 136 个，即 alpha、beta 各 68 个。也可以 `grep -i pseudopotential job.log` 查看 ECP 参数段是否出现。

### 8.3 从 Basis Set Exchange 下载并粘贴

[Basis Set Exchange](https://www.basissetexchange.org/)（BSE）是获取各种基组和 ECP 的标准来源。

**方法一：网页**

1. 打开 https://www.basissetexchange.org/ ，在周期表中选择元素（如 Ir）；
2. 在左侧选择基组（如 def2-TZVP）；
3. 格式（Format）选 **Gaussian**，点 Get Basis Set，复制或下载文本。

**方法二：命令行（API 或 Python 包）**

```bash
# REST API: def2-TZVP for Ir (Z = 77) in Gaussian format, basis + ECP in one file
curl -s "https://www.basissetexchange.org/api/basis/def2-tzvp/format/gaussian94/?version=1&elements=77" -o ir_def2tzvp.gbs

# or the official Python package (pip install basis_set_exchange); see 'bse --help'
bse get-basis def2-TZVP gaussian94 --elements Ir > ir_def2tzvp.gbs
```

下载的文件开头是若干 `!` 注释行和**空行**，然后是基组块（以 `****` 结束），再隔一个空行是 ECP 块。粘贴时：

1. **删掉文件头的空行**（`!` 注释行可以保留，也可以删掉）。空行会被当成基组段的结束，这是从 BSE 粘贴时最常见的错误。
2. 基组块放在 GenECP 的基组段中，ECP 块放在 ECP 段中，两段之间保留**一个**空行。
3. 其他元素可以继续用内置基组名。

[`ir_complex_gen_def2tzvp_from_bse.gjf`](../examples/03-method-basis/ir_complex_gen_def2tzvp_from_bse.gjf) 的结构（数值行有删节）：

```text
#p PBE1PBE/GenECP EmpiricalDispersion=GD3BJ SP

trans IrCl(CO)(PH3)2, def2TZVP everywhere, Ir basis and ECP pasted from BSE

0 1
Ir     0.000000    0.000000    0.000000
... (12 atoms) ...

C O P H Cl 0
Def2TZVP
****
! Ir orbital basis: def2-TZVP from the Basis Set Exchange (BSE 0.12,
! data from Turbomole 7.3), Gaussian94 format. BSE header removed.
Ir     0
S    3   1.00
     30.000000000            0.30797903228
     ...
F    1   1.00
      0.6100800              1.0000000
****

IR     0
IR-ECP     3     60
f potential
  1
2      3.03407192            21.53103107
s-f potential
  ...

```

ECP 块的格式：中心行 `IR 0`；`名称 最高角动量 芯电子数`（这里 `IR-ECP 3 60` 表示 f 为最高投影、替换 60 个芯电子）；然后每个分量依次是一行描述（如 `f potential`、`s-f potential`）、项数、每项一行 `r 的幂次  指数  系数`。自己命名 ECP 时不要用 Gaussian 内置的名字（CEP、CHF、LANL1、LANL2、LP-31、SDD、SHC）。

**验证粘贴是否正确**：Ir 的 def2-TZVP 在 Gaussian 中也是内置的，所以示例目录同时提供了 [`ir_complex_def2tzvp_builtin.gjf`](../examples/03-method-basis/ir_complex_def2tzvp_builtin.gjf)（`PBE1PBE/Def2TZVP`，ECP 自动加上）。两者的 `basis functions` 数应当相同，能量应一致或只差数值噪声。用这种"对照组"的方法检查手工粘贴的基组是个好习惯。真正需要 BSE 的场合是 Gaussian 没有内置的基组，例如 def2-TZVPD、ma-def2-TZVP、cc-pVTZ-PP、x2c-TZVPall、LANL08(f)。

```bash
grep -m1 "basis functions" ir_complex_*.log
grep "SCF Done" ir_complex_*.log
```

> [!TIP]
> 想把输出中实际使用的基组打印成 Gen 格式，可以在 route 中加 `GFInput`；加 `GFPrint` 则以表格形式打印基函数。用它们可以核对内置基组的定义，或把一次计算的基组原样搬到另一个输入中。

### 8.4 Pseudo=Read 与 ExtraBasis

- **`Pseudo=Read`**（同义词 `Pseudo=Cards`）：从输入读 ECP。`Gen Pseudo=Read` 与 `GenECP` 等价。官方文档特别指出，在 ONIOM 中只想给某一层读入 ECP 时应使用 `GenECP`，因为 `Pseudo` 会作用于所有层。
- **`ExtraBasis`**：在内置基组的基础上**追加**函数，例如给某个原子加一个弥散函数，或给金属加极化函数。追加段的格式与 Gen 相同。注意 d/f 函数的类型（5D/6D）由 route 中写的那个基组决定。

给 Pd 的 LANL2DZ 加一组 f 极化函数（LANL2DZ(f)，指数 1.472 取自 Ehlers 等 1993 年的工作）。在 GenECP 中再加一个 Pd 块即可，函数会累加：

```text
C H P Cl 0
6-31G(d)
****
Pd 0
LANL2DZ
****
Pd 0
F   1   1.00
      1.4720000              1.0000000
****

Pd 0
LANL2DZ

```

> [!NOTE]
> 不同文献、不同版本的 f 指数可能略有差别。使用前请对照原文（Ehlers et al., *Chem. Phys. Lett.* 1993）或 BSE 上的 LANL2DZ(f)/LANL08(f) 数据，并在论文中注明来源。

### 8.5 Gen/GenECP 格式自查清单

- [ ] route 写的是 `GenECP`（有 ECP 时）还是 `Gen`（只有全电子基组时）？
- [ ] 分子中的**每一种元素**都在基组段中出现，且只出现在需要的块中？中心行里没有分子中不存在的元素？
- [ ] 每个中心块都以 `****` 结尾？
- [ ] 基组段和 ECP 段之间恰好一个空行，ECP 段之后还有空行（文件末尾空行）？
- [ ] 从 BSE 粘贴的内容中没有残留的空行？
- [ ] 有 ECP 的元素在基组段和 ECP 段**都**有定义？
- [ ] 5D/6D 设置与对比的计算一致？
- [ ] 输出中的电子数、基函数数与预期一致（`grep -m1 "alpha electrons"`、`grep -m1 "basis functions"`）？

---

## 9. 积分格点

### 9.1 DFT 数值积分

交换-相关能无法解析积分，Gaussian 在每个原子周围的径向 × 角向格点上做数值积分。格点越密越准确，也越慢。

| 格点 | 关键词 | 规模（官方 Integral 文档） | 说明 |
|---|---|---|---|
| SG1 | `Int=SG1Grid` | 剪枝 (50,194) | 不推荐用于生产计算 |
| Fine | `Int=FineGrid` | 剪枝 (75,302)，约 7000 点/原子 | **G09 的默认** |
| UltraFine | `Int=UltraFine` | 剪枝 (99,590) | **G16 的默认** |
| SuperFine | `Int=SuperFine` | 第一行原子剪枝 (175,974)，第二行及以后 (250,974) | 难收敛的优化、低频模式、残留小虚频时使用 |

剪枝格点只为 H–Kr 定义。CPHF（解析二阶导数等）默认使用比积分格点低两档的格点：UltraFine 对应 SG1。

### 9.2 何时需要更密的格点

- **M06 系列等 meta-GGA 泛函对格点敏感**，能量随分子在空间中的取向变化会出现"噪声"（Wheeler & Houk, *JCTC* 2010）。G16 默认的 UltraFine 已经大大改善了这一问题；做低频模式、构象能差很小的比较或优化不易收敛时，可以改用 `Int=SuperFine`。
- 优化后残留很小的虚频（几十 cm⁻¹ 以下）：用更密的格点重新优化并计算频率（[第6章](06-freq-thermo.md)）。
- **需要相互比较的能量必须用同一个格点算**，包括同一项目中的反应物、产物、过渡态，以及优化和频率。

### 9.3 G16 与 G09 默认值的差别与 G09Defaults

G16 相对 G09 改变了四项默认值：

| 项目 | G09 | G16 |
|---|---|---|
| DFT 积分格点 | FineGrid | UltraFine |
| 双电子积分精度 | $10^{-10}$ | $10^{-12}$ |
| SCRF 默认形式 | 非对称 IEFPCM | 对称 IEFPCM |
| 物理常数 | 2006 版 | 2010 版 |

`G09Defaults` 关键词一次性把这四项恢复成 G09 的设置。它只用于**复现旧工作**，新研究应使用 G16 默认值。用不同版本或不同格点算出的能量直接相减，可能引入零点几 kcal/mol 甚至更大的人为误差。

---

## 10. 双杂化泛函与 MP2：成本、冻芯与输出

### 10.1 成本

双杂化泛函先做一次杂化 DFT 的 SCF，再算 MP2 型的二阶微扰相关能，成本与 MP2 相同（$N^5$），而且需要较多的内存（`%mem`）和磁盘（`MaxDisk`）。实践中：

- **只用于单点**，几何来自更便宜的方法（high//low）。
- 基组至少三 ζ（def2-TZVP、def2-TZVPP 或 def2-QZVP），否则 PT2 部分误差较大。
- `%mem` 给得充足可以显著提速（MP2 算法会根据 `%Mem` 和 `MaxDisk` 自动选择）。

### 10.2 输出中的能量

双杂化计算中，`SCF Done` 行**不是**最终能量，它只是 SCF 部分（不含 PT2 相关能）。最终能量在后面的 PT2 结果行中，以 Fortran 的 D 指数格式打印。官方 MP2 文档中 B2PLYP 的例子形如：

```text
 E2(B2PLYP) =    -0.3262340664D-01 E(B2PLYP) =    -0.39113226645200D+02
```

`E2(...)` 是 PT2 相关能的贡献，`E(...)` 是总能量（上例为 −39.113226645 Hartree）。其他双杂化的输出与此类似，括号中是方法名，以你实际的输出为准：

```bash
grep "SCF Done" dsdpbep86_def2tzvp_sp.log     # NOT the final energy
grep "E2(" dsdpbep86_def2tzvp_sp.log          # last field = total double-hybrid energy
```

同理，MP2 计算的最终能量在 `EUMP2 =` 一行。用脚本批量提取能量时一定要区分这些情况（仓库脚本 [`g16_log_summary.py`](../scripts/g16_log_summary.py) 汇总的能量类型，用 `python ../scripts/g16_log_summary.py -h` 查看）。

### 10.3 冻芯（frozen core）

Gaussian 中 MP2 和双杂化**默认冻芯（FC）**：内层轨道不参与相关能计算。官方文档对 FC 的定义是：6-31G 和 6-311G 基组按 G2 约定冻结，其他基组冻结最大的稀有气体芯，第三周期起的碱金属和碱土金属的外层 s、p 芯轨道不冻结。

| 写法 | 含义 |
|---|---|
| 默认（`FC`） | 冻芯 |
| `MP2(Full)`、`MP2=Full` | 全电子相关 |
| `MP2=FreezeG2` 等 | 其他冻结方案 |

用 ECP 的元素，芯电子已被 ECP 替换，冻芯只作用于剩下的价层内部。与文献比较时，要确认原文用的是冻芯还是全电子；DSD 类泛函的参数是在特定冻芯约定下拟合的，**不要随意改用 Full**。

### 10.4 Gaussian 中的其他高精度选择

Gaussian 还提供 CCSD(T)（`CCSD(T)/aug-cc-pVTZ`）以及 G4、CBS-QB3 等复合方法。它们适合小分子的高精度热化学或作为基准。对于几十个原子以上的体系，Gaussian 本身没有提供局域相关（DLPNO 类）方法，可以考虑外部程序（如 ORCA）做 DLPNO-CCSD(T) 单点。

---

## 11. 本章示例一览

| 文件 | 演示内容 |
|---|---|
| [`benzene_b3lyp_gd3bj_def2svp.gjf`](../examples/03-method-basis/benzene_b3lyp_gd3bj_def2svp.gjf) | B3LYP-D3(BJ)/def2-SVP Opt Freq（low 水平） |
| [`dsdpbep86_def2tzvp_sp.gjf`](../examples/03-method-basis/dsdpbep86_def2tzvp_sp.gjf) | `%oldchk` + `Geom=AllCheck` 的 DSD-PBEP86/def2-TZVP 单点（high 水平） |
| [`ethylene_m062x_gd3_def2tzvp_sp.gjf`](../examples/03-method-basis/ethylene_m062x_gd3_def2tzvp_sp.gjf) | M06-2X 只能配 GD3 |
| [`pd_complex_genecp_lanl2dz.gjf`](../examples/03-method-basis/pd_complex_genecp_lanl2dz.gjf) | GenECP：6-31G(d) + LANL2DZ/ECP，`6D` |
| [`ir_complex_gen_def2tzvp_from_bse.gjf`](../examples/03-method-basis/ir_complex_gen_def2tzvp_from_bse.gjf) | 粘贴 BSE 下载的 Ir 基组与 ECP |
| [`ir_complex_def2tzvp_builtin.gjf`](../examples/03-method-basis/ir_complex_def2tzvp_builtin.gjf) | 内置 `Def2TZVP` 作对照 |

运行顺序和预期结果见 [示例 README](../examples/03-method-basis/README.md)。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 说明 |
|---|---|---|
| 有机分子优化/频率 | B3LYP-D3(BJ)/def2-SVP 或 6-31G(d) | 成本低、稳健 |
| 高精度单点 | ωB97X-D、M06-2X-D3(0) 或 DSD-PBEP86 加 def2-TZVP（或更大） | 至少三 ζ |
| 色散校正 | 按 §4.2 的表选择 GD3BJ 或 GD3；自带色散的泛函不加 | 写清阻尼类型 |
| 过渡金属 | 全体系 def2 系列（自动带 ECP） | 比混合基组更一致，也更省事 |
| 必须用 LANL2DZ 时 | GenECP：金属 LANL2DZ(+f)，其他原子用 Pople 或 def2 | 不要对整个分子写 LANL2DZ |
| 阴离子、弱相互作用、激发态 | 加弥散函数 | 注意近线性相关 |
| 积分格点 | 默认 UltraFine；M06 类做频率、出现小虚频时用 SuperFine | 同一项目统一 |
| 双杂化 | 只做单点；`%mem` 给足；保持默认冻芯 | 最终能量看 `E2(` 那一行 |
| 版本一致性 | 同一项目用同一个 Gaussian 版本；不要混用 `G09Defaults` | 避免人为误差 |
| 方法描述 | 写成"方法-色散/基组//方法-色散/基组"，注明 ECP、溶剂、格点 | 便于复现 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `R6DS8: Unable to choose the S8 parameter`（l301） | 该泛函没有所选色散形式的参数，最常见的是 M06-2X 配 GD3BJ | 按 §4.2 改用 GD3，或换泛函 |
| `Atomic number out of range for 6-31G basis set.`（l301） | Pople 基组只到 Kr，分子中有 I、Pd 等 | 全体系用 def2 系列，或 GenECP |
| `End of file reading basis center.`（l301） | Gen 段格式错误：缺 `****`、多余空行（常见于 BSE 粘贴）、文件提前结束 | 按 §8.5 逐项检查 |
| `EOF while reading ECP pointer card.`（l301） | 用了 GenECP 却没写 ECP 段；ECP 段前缺空行；文件末尾缺空行 | 补 ECP 段和空行 |
| 作业在 L301 中止，提示 Gen 中的某个原子不存在 | 中心行列出了分子里没有的元素 | 删掉该元素，或在库文件中用 `-X` 写法 |
| 某元素没有基函数，L301 报错 | Gen 段漏写了该元素 | 补上该元素的中心块 |
| `QPErr`，指向泛函名 | 关键词写错，如 `PBE0`、`M06-2X`、`wB97X-D` | 改成 `PBE1PBE`、`M062X`、`wB97XD` |
| 含 Pd 的计算能量异常或 SCF 不收敛，但没报错 | 用 `Gen` 没读 ECP，或 ECP 段写错 | 改用 `GenECP`，检查电子数 |
| 双杂化能量看起来"和杂化泛函差不多" | 取的是 `SCF Done`，不是 `E(...)` | 按 §10.2 读取 |
| 同一分子 G09 与 G16 结果差零点几 kcal/mol | 默认格点、积分精度、SCRF、物理常数不同 | 统一版本；只在复现旧工作时用 `G09Defaults` |
| `NBsUse` 小于 `NBasis` | 弥散函数引起近线性相关，程序自动剔除了部分组合 | 通常无需处理；SCF 有问题时见第4章 |
| 同一项目中不同计算的基函数数对不上 | 5D/6D 不一致（Gen 默认 5D，内置 6-31G(d) 为 6D） | 统一加 `6D` 或 `5D` |

---

## 练习

1. 对 `benzene_b3lyp_gd3bj_def2svp.gjf` 分别去掉色散校正、改为 `GD3`，比较三者的优化能量和 C–C 键长。苯是单体，色散影响很小；再对平行堆积的苯二聚体做同样的比较，体会色散在分子间作用中的作用（[第14章](14-weak-bsse.md)）。
2. 用 6-31G(d)、def2-SVP、def2-TZVP、def2-QZVP 计算乙烯的 B3LYP-D3(BJ) 单点能，核对每个基组的基函数数，画出能量随基函数数的变化曲线。
3. 把 `ethylene_m062x_gd3_def2tzvp_sp.gjf` 的 `GD3` 改为 `GD3BJ`，运行并找到 R6DS8 报错所在的 Link。
4. 运行两个 Ir 示例，比较基函数数和能量。然后故意在粘贴的 Ir 基组块中间插入一个空行，观察报错。
5. 在 `pd_complex_genecp_lanl2dz.gjf` 中去掉 `6D`，比较基函数数和能量的变化。再按 §8.4 给 Pd 加上 f 极化函数。
6. 从 BSE 下载 def2-TZVPD 的 O 和 H 基组，用 Gen 计算 OH⁻ 的单点，与 `Def2TZVP` 的结果比较。

---

## 参考资料

**官方文档**

- [Density Functional (DFT) Methods](https://gaussian.com/dft/)：泛函关键词、EmpiricalDispersion 的 GD2/GD3/GD3BJ 参数表、积分格点精度说明
- [MP Methods](https://gaussian.com/mp/)：MP2 与双杂化 B2PLYPD3、DSDPBEP86、PBE0DH、PBEQIDH
- [Window Keyword and Frozen Core Options](https://gaussian.com/frozencore/)
- [Basis Sets](https://gaussian.com/basissets/)：内置基组、覆盖元素、5D/6D、密度拟合
- [Gen / GenECP / ExtraBasis](https://gaussian.com/gen/)
- [Pseudo (ECP input)](https://gaussian.com/pseudo/)：ECP 输入格式、Stuttgart ECP 可用性表
- [Integral keyword](https://gaussian.com/integral/)：格点、DKH 相对论
- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)：G09 到 G16 的默认值变化与 G09Defaults
- [Basis Set Exchange](https://www.basissetexchange.org/)：下载 Gaussian 格式的基组与 ECP

**相关 GitHub 仓库**（仅外链，更多见 [相关仓库汇总](../../Resources/github-repos.md)）

- [patonlab/pyDFTD3](https://github.com/patonlab/pyDFTD3)：独立计算 D3(0)/D3(BJ) 色散能，可用来理解 EmpiricalDispersion 的数值贡献
- [igor-1982/xDH4Gau](https://github.com/igor-1982/xDH4Gau)：在 Gaussian 中计算 XYG3/XYGJ-OS 等 xDH 型双杂化泛函
- [thisisntnathan/dftCourse](https://github.com/thisisntnathan/dftCourse)：课程中有 Gen 基组和效率相关章节
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：覆盖 GenECP 等任务的中文工作流项目
- [cclib/cclib](https://github.com/cclib/cclib)：解析 log 中的能量、基组信息，可用于交叉验证

**论文**

- J. P. Perdew, K. Schmidt, *AIP Conf. Proc.* **2001**, 577, 1. Jacob 阶梯。DOI: [10.1063/1.1390175](https://doi.org/10.1063/1.1390175)
- A. D. Becke, *J. Chem. Phys.* **1993**, 98, 5648. DOI: [10.1063/1.464913](https://doi.org/10.1063/1.464913)
- Y. Zhao, D. G. Truhlar, *Theor. Chem. Acc.* **2008**, 120, 215. M06 系列。DOI: [10.1007/s00214-007-0310-x](https://doi.org/10.1007/s00214-007-0310-x)
- J.-D. Chai, M. Head-Gordon, *Phys. Chem. Chem. Phys.* **2008**, 10, 6615. ωB97X-D。DOI: [10.1039/b810189b](https://doi.org/10.1039/b810189b)
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* **2010**, 132, 154104. D3。DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* **2011**, 32, 1456. D3(BJ)。DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- S. Kozuch, J. M. L. Martin, *J. Comput. Chem.* **2013**, 34, 2327. DSD-PBEP86 等。DOI: [10.1002/jcc.23391](https://doi.org/10.1002/jcc.23391)
- F. Weigend, R. Ahlrichs, *Phys. Chem. Chem. Phys.* **2005**, 7, 3297. def2 基组。DOI: [10.1039/b508541a](https://doi.org/10.1039/b508541a)
- P. J. Hay, W. R. Wadt, *J. Chem. Phys.* **1985**, 82, 299. LANL2DZ ECP。DOI: [10.1063/1.448975](https://doi.org/10.1063/1.448975)
- A. W. Ehlers et al., *Chem. Phys. Lett.* **1993**, 208, 111. 过渡金属 ECP 的 f 极化函数。DOI: [10.1016/0009-2614(93)80086-5](https://doi.org/10.1016/0009-2614(93)80086-5)
- B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson, T. L. Windus, *J. Chem. Inf. Model.* **2019**, 59, 4814. Basis Set Exchange。DOI: [10.1021/acs.jcim.9b00725](https://doi.org/10.1021/acs.jcim.9b00725)
- S. E. Wheeler, K. N. Houk, *J. Chem. Theory Comput.* **2010**, 6, 395. M06 类泛函的积分格点误差。DOI: [10.1021/ct900639j](https://doi.org/10.1021/ct900639j)
- L. Goerigk et al., *Phys. Chem. Chem. Phys.* **2017**, 19, 32184. GMTKN55 基准。DOI: [10.1039/C7CP04913G](https://doi.org/10.1039/C7CP04913G)
- M. Bursch, J.-M. Mewes, A. Hansen, S. Grimme, *Angew. Chem. Int. Ed.* **2022**, 61, e202205735. DFT 最佳实践。DOI: [10.1002/anie.202205735](https://doi.org/10.1002/anie.202205735)

---

[← 上一章](02-input-anatomy.md) | [返回目录](../README.md) | [下一章 →](04-sp-scf.md)
