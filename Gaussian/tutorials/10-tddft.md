# 第10章 激发态 TD-DFT：UV-Vis、荧光发射与 ECD（Excited States with Gaussian TD-DFT: UV-Vis, Fluorescence & ECD）

本章用 Gaussian 16 的 TD-DFT（`TD` / `TDA` 关键词）完成几类最常见的激发态任务：UV-Vis 吸收光谱、S1 激发态优化与荧光发射、T1 与磷光能量、ECD（电子圆二色谱）定绝对构型，以及 NTO、DCT、空穴-电子分析等激发态性质分析。每一步都给出完整输入、关键词逐项解释和输出检查方法，同时讲清泛函选择（CAM-B3LYP、ωB97X-D、M06-2X、PBE0）、溶剂模型（线性响应、cLR、ExternalIteration）和谱带展宽中的常见陷阱。

## 1. 本章目标与适用场景

学完本章你应该能够：

- 为有机分子、染料或配合物选择合适的泛函和基组做 TD-DFT 计算，并知道什么时候 B3LYP 会失效（电荷转移态、Rydberg 态）；
- 写出 `TD=(NStates=..., Root=..., Singlets/Triplets/50-50)` 和 `TDA` 输入，读懂 `Excited State` 输出中的激发能、波长、振子强度和组态系数；
- 用高斯展宽把离散的激发态画成 UV-Vis 和 ECD 谱，做 Boltzmann 加权和对映体翻转；
- 优化 S1、T1 结构，计算垂直发射能、绝热激发能、0-0 能量和 Stokes 位移；
- 用 NTO、`Pop=DCT` 和 Multiwfn 空穴-电子分析判断激发是局域激发（LE）还是电荷转移（CT）；
- 在溶液中区分平衡/非平衡溶剂化，并会使用态特定方法 cLR 与 ExternalIteration。

前置知识：几何优化与频率（[第5章](05-opt.md)、[第6章](06-freq-thermo.md)）、隐式溶剂（[第9章](09-solvation.md)）。开壳层体系的激发态先读 [第13章](13-openshell-stability.md) 的稳定性检查。DFT 理论背景可参考 [Common 第7章](../../Common/tutorials/07-learning-path.md)。

| 任务 | 本章小节 | 配套示例 |
|---|---|---|
| UV-Vis 吸收光谱 | 第 4、5 节 | [10-uvvis-emission-coumarin](../examples/10-uvvis-emission-coumarin/README.md) |
| S1 优化、荧光、态特定溶剂化 | 第 6 节 | 同上 |
| T1 优化、磷光 | 第 7 节 | 同上 |
| ECD 与绝对构型 | 第 8 节 | [10-ecd-methyloxirane](../examples/10-ecd-methyloxirane/README.md) |
| NTO、DCT、空穴-电子分析 | 第 9 节 | [10-nto-dct-dmabn](../examples/10-nto-dct-dmabn/README.md) |

## 2. 原理速览

### 2.1 线性响应 TD-DFT 与 TDA

TD-DFT 在基态 Kohn-Sham 解的基础上求解线性响应方程（Casida 方程），得到激发能 $\omega$ 和激发向量 $(\mathbf X,\mathbf Y)$：

$$
\begin{pmatrix}\mathbf A & \mathbf B\\ \mathbf B^{*} & \mathbf A^{*}\end{pmatrix}
\begin{pmatrix}\mathbf X\\ \mathbf Y\end{pmatrix}
=\omega
\begin{pmatrix}\mathbf 1 & \mathbf 0\\ \mathbf 0 & -\mathbf 1\end{pmatrix}
\begin{pmatrix}\mathbf X\\ \mathbf Y\end{pmatrix}
$$

$\mathbf X$ 对应占据轨道 i 到空轨道 a 的激发（输出里写作 `i -> a`），$\mathbf Y$ 对应退激发（输出里写作 `i <- a`）。Tamm-Dancoff 近似（TDA）令 $\mathbf B=0$，方程退化为 $\mathbf A\mathbf X=\omega\mathbf X$。TDA 更便宜，对三重态不稳定性和接近锥形交叉的区域更稳健，但振子强度不再满足求和规则，绝对值会有偏差。Gaussian 用 Davidson 迭代只求最低的 `NStates` 个根，所以**算多少个态必须由你指定**。

### 2.2 垂直、绝热与 0-0 能量

以单重态吸收与荧光为例（$R_{S0}$、$R_{S1}$ 分别为 S0、S1 的平衡结构）：

| 量 | 定义 | Gaussian 中怎么得到 |
|---|---|---|
| 垂直吸收能 | $E_{S1}(R_{S0})-E_{S0}(R_{S0})$ | S0 结构上的 TD 单点，`Excited State 1` |
| 垂直发射能 | $E_{S1}(R_{S1})-E_{S0}(R_{S1})$ | S1 优化最后一步的 `Excited State 1` |
| 绝热激发能 | $E_{S1}(R_{S1})-E_{S0}(R_{S0})$ | 两个优化作业的总能量差 |
| 0-0 能量 | $E_{\rm adia}+\Delta{\rm ZPE}$ | 再加 S1 与 S0 的零点能差 |
| Stokes 位移 | $E_{\rm abs}^{\rm vert}-E_{\rm em}^{\rm vert}$ | 常换算成 cm⁻¹ 与实验比较 |

按 Kasha 规则，荧光通常来自最低激发单重态 S1，磷光来自 T1（存在例外，如薁）。单位换算：$\lambda({\rm nm}) = 1239.84/E({\rm eV})$，$1\ {\rm eV} = 8065.54\ {\rm cm^{-1}}$。

### 2.3 振子强度与旋转强度

振子强度（原子单位）：

$$f_{0n}=\frac{2}{3}\,\Delta E_{0n}\,\left|\langle 0|\hat{\boldsymbol\mu}|n\rangle\right|^2$$

它决定 UV-Vis 峰的强弱。ECD 的强度由旋转强度决定：

$$R_{0n}=\mathrm{Im}\left[\langle 0|\hat{\boldsymbol\mu}|n\rangle\cdot\langle n|\hat{\mathbf m}|0\rangle\right]$$

其中 $\hat{\mathbf m}$ 为磁偶极算符。有限基组下，长度规范（length gauge）的 $R$ 依赖坐标原点，速度规范（velocity gauge）的 $R$ 与原点无关，因此**做 ECD 时一般使用 R(velocity)**。

### 2.4 泛函与基组选择

| 泛函 | Gaussian 关键词 | 精确交换比例（约） | 适用 | 注意 |
|---|---|---|---|---|
| B3LYP | `B3LYP` | 20% | 局域 ππ* 激发的快速估计 | CT、Rydberg 激发能严重低估，可能出现虚假低能 CT 态 |
| PBE0 | `PBE1PBE` | 25% | 有机染料价层激发，λmax 常有误差抵消 | CT 态仍然偏低 |
| M06-2X | `M062X` | 54% | 价层激发、部分 CT 改善 | 对积分格点敏感（G16 默认 UltraFine 一般够用） |
| CAM-B3LYP | `CAM-B3LYP` | 短程 19% → 长程 65% | CT、Rydberg、ECD，激发态优化 | 局域价层激发通常略高估 |
| ωB97X-D | `wB97XD` | 短程约 22% → 长程 100% | CT 体系、长链共轭分子 | 同上；已内含色散，不要再加 GD3BJ |
| LC-ωhPBE | `LC-wHPBE` | 短程 0 → 长程 100% | 适合做 ω 最优调节（optimal tuning） | 默认 ω 对特定分子未必最优 |

为什么 CT 态需要范围分离泛函：CT 激发能在大分离距离下应近似为 $IP_D-EA_A-1/R$，常规杂化泛函的交换核衰减太快，无法给出 $-1/R$ 行为，CT 能量因此被系统性低估，而且低估幅度随 D-A 距离增大。长程含 100% HF 交换的泛函（ωB97X-D、LC-ωhPBE）能恢复正确的渐近行为。

基组：价层激发用 def2-TZVP 或 6-311+G(d,p) 通常足够，def2-SVP 适合激发态优化等耗时任务。Rydberg 态、阴离子和 ECD（低能区常有 Rydberg 成分）需要弥散函数，例如 aug-cc-pVDZ、6-311++G(d,p)。根据大规模基准测试（如 Laurent & Jacquemin 2013 的综述），合适泛函下价层垂直激发能的典型误差在 0.2–0.3 eV 量级，不要期望 TD-DFT 给出 0.05 eV 级别的精度。

> [!IMPORTANT]
> 计算得到的垂直激发能并不等于实验吸收峰位置 λmax：实验谱带包含振动结构、溶剂与温度效应。比较时建议同时看能量差（eV）和谱形，不要只比较 nm 数值。

## 3. TD 关键词速查

`TD` 是方法类关键词，和泛函/基组写在一起：`#p CAM-B3LYP/def2TZVP TD=(NStates=20)`。下表选项均来自官方 [TD 关键词页](https://gaussian.com/td/)：

| 选项 | 默认值 | 说明 |
|---|---|---|
| `NStates=M` | 3 | 求解的激发态数目；与 `50-50` 联用时是每种自旋各 M 个 |
| `Root=N` | 1 | “感兴趣的态”：激发态优化、频率、`Density=Current` 等都针对它 |
| `Singlets` | 闭壳层默认 | 只求单重态（仅对闭壳层单重态参考有效） |
| `Triplets` | — | 只求三重态（仅对闭壳层参考有效） |
| `50-50` | — | 单重态、三重态各一半 |
| `Read` | — | 从 chk 读取激发向量作为初猜（基组必须相同） |
| `Add=N` | — | 读取已收敛的态，再多算 N 个（隐含 `Read`） |
| `Restart` | — | 从 chk 续算中断的 TD 作业，重复原 route 加上该选项即可 |
| `EqSolv` / `NonEqSolv` | 单点为非平衡；激发态优化和请求激发态密度时为平衡 | PCM 平衡/非平衡溶剂化 |
| `Conver=N` | 单点 4，梯度 6 | 能量收敛到 $10^{-N}$ |
| `GOccSt`、`GOccEnd`、`GDEMin`、`DEMin` | — | 限定能量窗口，用于只关心高能区（如芯层激发）的情况 |

其他相关要点：

- `TDA` 关键词与 `TD` 用法相同：`TDA=(NStates=10)`；
- 开壳层参考（如双重态自由基）时 `Singlets/Triplets/50-50` 会被忽略，输出中每个态带 `<S**2>`，需检查自旋污染；
- `#p` 会打印更多信息；`IOp(9/40=N)` 把组态系数的打印阈值从默认 0.1 降到 $10^{-N}$（官方 IOp 说明：阈值 = $10^{-N}$），Multiwfn 激发态分析要求 `IOp(9/40=4)`，大体系可用 3；
- Gaussian 16 的 TD-DFT 频率默认为解析二阶导数（G09 只能数值差分），B.01 起支持 `TD Freq=Raman`（激发态静态 Raman，耗时约为普通频率的 7 倍）；
- 官方说明：G96、P86、PKZB、wPBEh、PBEh 交换泛函，PKZB 相关泛函，以及 OHSE1PBE、OHSE2PBE 杂化泛函**不能**做 TD 梯度和频率。

## 4. 实战一：香豆素的 UV-Vis 吸收光谱

示例目录：[10-uvvis-emission-coumarin](../examples/10-uvvis-emission-coumarin/README.md)。香豆素（coumarin，C9H6O2，17 个原子）是荧光染料的母核，体积小，适合练手。

### 4.1 完整输入：S0 优化 + 垂直激发

两步用 `--Link1--` 串联：第 1 步在 SMD(乙醇) 中优化基态并做频率，第 2 步在相同溶剂下做 TD 单点。

```text
! Chapter 10 example: coumarin UV-Vis absorption in ethanol (SMD)
! Step 1: S0 opt+freq   Step 2: vertical TD-DFT absorption, 20 singlets
%chk=coumarin_s0.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2SVP Opt Freq SCRF=(SMD,Solvent=Ethanol)

Coumarin S0 opt+freq, CAM-B3LYP/def2-SVP, SMD(ethanol)

0 1
 O       3.432361     1.735762    -0.161262
 C       2.414341     1.055850    -0.099047
 C       2.498112    -0.410001     0.028425
 C       1.413650    -1.187845     0.099364
 C       0.094269    -0.579979     0.050257
 C      -1.075754    -1.338826     0.119792
 C      -2.322908    -0.707303     0.068413
 C      -2.405880     0.682098    -0.052401
 C      -1.242447     1.446105    -0.122366
 C      -0.007031     0.804837    -0.070103
 O       1.144112     1.613186    -0.143894
 H       3.499281    -0.821949     0.061392
 H       1.498663    -2.264641     0.192933
 H      -1.032587    -2.421035     0.213956
 H      -3.232915    -1.300530     0.122758
 H      -3.377925     1.167635    -0.091865
 H      -1.297342     2.526637    -0.216351

--Link1--
%oldchk=coumarin_s0.chk
%chk=coumarin_td.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2TZVP TD=(NStates=20,Singlets) SCRF=(SMD,Solvent=Ethanol)
   Geom=Check Guess=Read IOp(9/40=3)

Coumarin vertical absorption, 20 singlets, CAM-B3LYP/def2-TZVP, SMD(ethanol)

0 1

```

| 关键词 | 取值 | 说明 |
|---|---|---|
| `CAM-B3LYP/def2SVP Opt Freq` | 第 1 步 | 基态用与后续 S1 优化**相同**的泛函/基组，这样绝热激发能和 0-0 能量才有意义 |
| `SCRF=(SMD,Solvent=Ethanol)` | 两步一致 | 溶剂模型、溶剂名在整个流程中保持一致 |
| `TD=(NStates=20,Singlets)` | 20 个单重态 | 覆盖约 190–350 nm；`Singlets` 对闭壳层是默认值，写出来更清楚 |
| `def2TZVP` | TD 单点 | 激发能对基组较敏感，单点用更大基组 |
| `Geom=Check Guess=Read` | 读取优化结构与轨道 | 配合 `%oldchk`，不覆盖第 1 步的 chk |
| `IOp(9/40=3)` | 系数阈值 0.001 | 输出更多组态系数，便于判断激发性质 |

> [!NOTE]
> TD 单点在 PCM/SMD 下默认使用**非平衡**溶剂化（只有溶剂电子极化跟上垂直激发），这正是吸收过程的物理图像，无需额外关键词。

运行与检查：

```bash
g16 < coumarin_td_camb3lyp_smd.gjf > coumarin_td_camb3lyp_smd.log
grep -c "Normal termination" coumarin_td_camb3lyp_smd.log      # Link1 两步应为 3（Opt、Freq、TD）
grep "NImag" coumarin_td_camb3lyp_smd.log                      # 归档段中应为 NImag=0
grep "Excited State" coumarin_td_camb3lyp_smd.log | head -20
```

> [!TIP]
> `Opt Freq` 写在同一行时，Gaussian 在内部自动追加一个 Freq 步，所以两步输入会出现 3 个 `Normal termination`。也可以用仓库的 `python ../scripts/g16_log_summary.py -h` 查看如何批量汇总终止状态、能量与虚频。

### 4.2 读懂 Excited State 输出

（示意输出，数值仅供参考）

```text
 Excitation energies and oscillator strengths:

 Excited State   1:      Singlet-A'     4.3521 eV  284.88 nm  f=0.1873  <S**2>=0.000
      42 -> 44        -0.10512
      43 -> 44         0.68234
      43 <- 44        -0.05231
 This state for optimization and/or second-order correction.
 Total Energy, E(TD-HF/TD-DFT) =  -496.xxxxxxxxx
 Copying the excited state density for this state as the 1-particle RhoCI density.

 Excited State   2:      Singlet-A"     4.6012 eV  269.46 nm  f=0.0004  <S**2>=0.000
      41 -> 44         0.69517
```

| 字段 | 含义 |
|---|---|
| `Excited State   1:` | 第 1 个激发态（按能量排序） |
| `Singlet-A'` | 自旋多重度与不可约表示；对称性无法确定时显示 `Singlet-?Sym` |
| `4.3521 eV  284.88 nm` | 垂直激发能及对应波长 |
| `f=0.1873` | 振子强度；f < 0.01 为弱/禁阻跃迁，> 0.1 为强吸收 |
| `<S**2>` | 激发态的自旋平方期望值，闭壳层参考时为 0 |
| `43 -> 44  0.68234` | 轨道 43（HOMO）到 44（LUMO）的激发系数 X |
| `43 <- 44` | 退激发系数 Y，只在 TD（非 TDA）中出现 |
| `Total Energy, E(TD-HF/TD-DFT)` | Root 指定态的总能量（Hartree），用于计算绝热激发能 |

闭壳层单重态中，组态 $i\to a$ 的贡献约为 $2c^2$（Gaussian 的归一化条件为 $\langle X+Y|X-Y\rangle=1$）。上面 `43 -> 44  0.68234` 的贡献约为 $2\times0.682^2\approx 93\%$，即 S1 基本是 HOMO→LUMO 激发。一个激发态若由多个系数相当的组态构成，用单一“HOMO→LUMO”描述就不合适，应改用第 9 节的 NTO。

快速整理成表（激发态序号、对称性、eV、nm、f）：

```bash
grep "Excited State" coumarin_td_camb3lyp_smd.log | awk '{print $3, $4, $5, $7, $9}'
```

HOMO、LUMO 编号可以从 `Alpha  occ. eigenvalues` 行数或 `NAE=`（α 电子数）读出：闭壳层下 HOMO = NAE。

### 4.3 NStates 够不够

- 最后一个激发态的波长应**短于**你想绘制的光谱下限，再留出展宽宽度（约 0.3–0.5 eV）的余量。例如要画到 200 nm（6.2 eV），最高激发态至少应到 6.5–6.7 eV；
- 不够就用 `TD=(Read,Add=10)` 在已有结果基础上追加 10 个态，不必从头算（基组、结构必须相同）；
- 只关心最低几个态时 `NStates=6~10` 足矣；激发态优化建议 `NStates ≥ Root+2`，以稳定 Davidson 迭代并监控态交叉。

## 5. 谱带展宽与绘图

### 5.1 UV-Vis 展宽公式

Gaussian 官方技术说明 [Creating UV/Visible Plots](https://gaussian.com/uvvisplot/) 给出的高斯线型（$\tilde\nu$ 为波数，σ 为 1/e 处的半宽，单位 cm⁻¹）：

$$
\varepsilon(\tilde\nu)=\sum_i 1.3062974\times10^{8}\,\frac{f_i}{\sigma}\exp\!\left[-\left(\frac{\tilde\nu-\tilde\nu_i}{\sigma}\right)^2\right]
\qquad (\varepsilon\ \text{单位 L·mol}^{-1}\text{·cm}^{-1})
$$

要点：

- 展宽要在**能量/波数空间**做，然后再换算到 nm 作图。直接在 nm 空间用等宽高斯会使长波区峰形失真；
- σ 与半高全宽的关系：${\rm FWHM}=2\sqrt{\ln 2}\,\sigma\approx1.665\sigma$。上述技术说明中 GaussView 的默认 σ 为 0.4 eV；文献中常用 FWHM 0.3–0.5 eV，必须在论文中写明；
- 也可用洛伦兹线型，但电子光谱习惯用高斯线型（溶液中非均匀展宽占主导）；
- 与实验比较时允许整体平移（例如 −0.2 eV），但平移量要在全体系中统一并注明。

### 5.2 三种绘谱途径

1. **GaussView 6**：打开 `.log`，菜单 Results → UV-VIS，可调半宽并导出数据；
2. **Multiwfn**：载入 `.log`，进入主功能 11（谱图绘制），选择 UV-Vis，可设置线型、FWHM、平移，导出 `spectrum_curve.txt`（Multiwfn 2026 起按日期命名版本，见 [Multiwfn 手册 3.13 节](http://sobereva.com/multiwfn/manual/)）；
3. **本仓库脚本**：`spectrum_uvvis_ecd.py` 解析 TD 输出、在 eV 空间做高斯展宽、输出 CSV 并画 UV-Vis/ECD 图，支持能量平移、对映体翻转和多构象 Boltzmann 加和。运行 `python ../scripts/spectrum_uvvis_ecd.py -h` 查看用法。

出版级作图（字体、配色、坐标轴）参见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md) 与 [第16章](16-postprocess-multiwfn.md)。也可以用 [cclib](https://github.com/cclib/cclib) 解析 `etenergies`、`etoscs`、`etrotats` 后自行绘图。

## 6. 实战二：S1 激发态优化与荧光发射

### 6.1 完整输入：S1 Opt + Freq

```text
! Chapter 10 example: optimize S1 of coumarin (TD-DFT gradient) + analytic TD frequencies
! Needs coumarin_s0.chk from coumarin_td_camb3lyp_smd.gjf (same functional/basis/solvent as the S0 opt)
%oldchk=coumarin_s0.chk
%chk=coumarin_s1.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2SVP TD=(NStates=6,Root=1) Opt Freq SCRF=(SMD,Solvent=Ethanol)
   Geom=Check Guess=Read NoSymm

Coumarin S1 opt+freq, CAM-B3LYP/def2-SVP, SMD(ethanol), equilibrium LR solvation

0 1

```

| 关键词 | 说明 |
|---|---|
| `TD=(NStates=6,Root=1)` | 优化第 1 激发态；多算几个态便于监控态交叉 |
| `Opt Freq` | 激发态优化后做解析 TD 频率，确认无虚频，并给出激发态 ZPE |
| `SCRF=(SMD,...)` | 激发态优化时 Gaussian 默认用**平衡**溶剂化（溶剂完全弛豫） |
| `NoSymm` | 激发态常使分子降低对称性（例如 nπ* 态让平面分子出平面）。关闭对称性，避免被约束在对称鞍点上 |
| `Geom=Check Guess=Read` | 从 S0 结构出发 |

可选项：`Opt=CalcFC` 在第一步算激发态 Hessian（官方 SCRF 页的乙醛示例就这么做），对难收敛或对称性破缺的情况有帮助，代价约为一次 TD 频率；`Opt=(MaxCycles=200)` 用于大体系。

运行与检查：

```bash
g16 < s1_opt_freq_emission.gjf > s1_opt_freq_emission.log
grep "Optimization completed" s1_opt_freq_emission.log
grep "NImag" s1_opt_freq_emission.log
grep "Total Energy, E(TD-HF/TD-DFT)" s1_opt_freq_emission.log | tail -1
grep "Excited State   1:" s1_opt_freq_emission.log | tail -1     # 最后一步即垂直发射能
grep "Zero-point correction" s1_opt_freq_emission.log
```

### 6.2 从输出得到发射能、0-0 能量与辐射速率

- **垂直发射能（线性响应、平衡溶剂化）**：优化收敛后最后一次出现的 `Excited State   1:` 的激发能；它的振子强度 f 决定荧光亮度；
- **绝热激发能**：$E_{\rm adia}=E_{\rm TD}(S1,R_{S1})-E_{\rm SCF}(S0,R_{S0})$，前者取 `Total Energy, E(TD-HF/TD-DFT)` 的最后一个值，后者取 S0 优化最后的 `SCF Done`，两者必须是同一泛函、基组和溶剂；
- **0-0 能量**：$E_{00}=E_{\rm adia}+[{\rm ZPE}(S1)-{\rm ZPE}(S0)]$。TD 频率的热化学部分基于激发态总能量，可核对 `Sum of electronic and zero-point Energies` 与 `Total Energy, E(TD-HF/TD-DFT)` 之差是否等于激发态 ZPE；
- **辐射速率估计**（Einstein 自发辐射系数，$\tilde\nu$ 单位 cm⁻¹，n 为折射率）：$k_r\approx \dfrac{n^2\,\tilde\nu^2 f}{1.499}\ {\rm s^{-1}}$。它只给出数量级，不能代替含振动耦合的速率理论。

> [!WARNING]
> **Root flipping（根翻转）**：`Root=1` 表示“按能量排序的第 1 个激发态”，不是“具有某种性质的那个态”。优化中如果 S1 与 S2 交叉，Gaussian 会继续跟随新的第 1 个态，你可能在不知情的情况下换到另一个电子态的势能面。检查方法：比较优化前后主要组态和振子强度是否突变，必要时用第 9 节的 NTO 对比；若确实想跟随某个态，可改变 `Root` 或从不同初始结构重新优化。

### 6.3 态特定溶剂化：cLR 与 ExternalIteration

线性响应（LR）PCM 用的是跃迁密度的响应，对**极性 CT 态**的溶剂稳定化描述不足。Gaussian 提供两种态特定（state-specific）修正：

| 方法 | 关键词 | 特点 |
|---|---|---|
| LR-PCM | `SCRF=(Solvent=...)` + `TD` | 默认；支持激发态梯度和频率 |
| 修正线性响应 cLR | `SCRF=(...,CorrectedLR)` | 用激发态弛豫密度对激发能做一阶修正，只需一步 |
| 外迭代 EI | `SCRF=(...,ExternalIteration)` | 激发态密度与反应场自洽迭代（Link 124），最严格；**只能用于能量计算** |

平衡/非平衡默认值（官方 [SCRF 页](https://gaussian.com/scrf/)）：TD 单点默认非平衡；TD 激发态优化默认平衡；ExternalIteration 默认平衡；`NonEquilibrium=Save` 把慢极化电荷存进 chk，下一步用 `NonEquilibrium=Read` 读入，从而得到“溶剂慢极化停留在上一个电子态”的非平衡能量。

官方 SCRF 页以乙醛 n→π* 为例给出 6 步完整流程（S0 优化 → LR 垂直吸收 → cLR/EI 态特定吸收 → S1 优化 → cLR/EI 平衡 S1 → 非平衡 S0），示例中 LR 给出的吸收为 283.28 nm，cLR 为 284.35 nm，EI 为 277.70 nm，cLR 发射为 400.79 nm（以上数值引自官方页面）。本仓库把发射部分写成了两步输入 `s1_emission_state_specific.gjf`：

```text
! Chapter 10 example: state-specific (corrected linear-response, cLR) vertical emission
! Follows the 2-step scheme of the official Gaussian SCRF page; uses IEFPCM, not SMD.
! Needs coumarin_s1.chk from s1_opt_freq_emission.gjf
! Emission energy = E(step 1, 'Total energy after correction') - E(step 2, SCF Done)
%oldchk=coumarin_s1.chk
%chk=coumarin_s1_clr.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2SVP TD=(NStates=6,Root=1) Geom=Check Guess=Read NoSymm
   SCRF=(IEFPCM,Solvent=Ethanol,CorrectedLR,NonEquilibrium=Save)

Coumarin S1 at S1 minimum, cLR equilibrium solvation, save slow solvent charges

0 1

--Link1--
%oldchk=coumarin_s1_clr.chk
%chk=coumarin_s0_at_s1.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2SVP Geom=Check Guess=Read NoSymm
   SCRF=(IEFPCM,Solvent=Ethanol,NonEquilibrium=Read)

Coumarin S0 at S1 geometry, non-equilibrium solvation read from the S1 step

0 1

```

- 第 1 步输出中查找 `PCM Corrected Linear-Response - Equilibrium solvation:` 下的 `Total energy after correction`；
- 第 2 步取 `SCF Done`；两者之差即态特定垂直发射能；
- 换成 EI：把第 1 步的 `CorrectedLR` 改为 `ExternalIteration`，能量从 `After PCM corrections, the energy is` 行读取。

```bash
grep -A3 "PCM Corrected Linear-Response" s1_emission_state_specific.log
grep "SCF Done" s1_emission_state_specific.log | tail -1
grep "After PCM corrections, the energy is" *.log      # 仅 ExternalIteration 作业
```

> [!NOTE]
> 官方示例使用默认的 IEFPCM，本仓库也照此写；S1 结构来自 SMD 优化，这一不一致对垂直发射能的影响通常很小，但应在论文方法部分说明。态特定方法对局域 ππ* 态的修正一般较小，对强 CT 态则可能明显改变能量顺序。

### 6.4 进阶：振动分辨的吸收/发射谱带

若要得到含振动精细结构的谱形，Gaussian 16 的 `Freq=(FC,ReadFC,ReadFCHT)` 可以基于 S0 与 S1 两个频率作业的 chk 做 Franck-Condon（`FC`）、Herzberg-Teller（`HT`）或 `FCHT` 计算；附加输入里用 `Transition=Emission` 模拟发射，用 `Spectroscopy=CD` 模拟振动分辨 ECD/CPL。输入格式和选项较多，请直接参考官方 [Freq 关键词页](https://gaussian.com/freq/) 的 ReadFCHT Input 部分，前提是两个态都是可靠的极小点（无虚频）。

## 7. 实战三：三重态 T1 与磷光

### 7.1 途径 A：UKS 基态优化（ΔSCF）

T1 是该自旋多重度下的最低态，可以直接用非限制性 DFT 当作“基态”优化，代价与普通基态优化相同：

```text
! Chapter 10 example: lowest triplet T1 of coumarin by unrestricted ground-state DFT (Delta-SCF)
! Step 1: UKS triplet opt+freq   Step 2: closed-shell S0 single point at the T1 geometry
! Vertical phosphorescence energy ~ E(T1, step 1 final) - E(S0, step 2)
%oldchk=coumarin_s0.chk
%chk=coumarin_t1.chk
%mem=8GB
%nprocshared=8
#p UCAM-B3LYP/def2SVP Opt Freq SCRF=(SMD,Solvent=Ethanol) Geom=Check

Coumarin T1 (UKS, multiplicity 3) opt+freq, CAM-B3LYP/def2-SVP, SMD(ethanol)

0 3

--Link1--
%oldchk=coumarin_t1.chk
%chk=coumarin_s0_at_t1.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/def2SVP SCRF=(SMD,Solvent=Ethanol) Geom=Check

Coumarin S0 single point at the T1 geometry

0 1

```

- `Geom=Check` 时电荷与多重度行仍然要写，可以改成 `0 3`；
- 检查自旋污染：`grep "S\*\*2" t1_uks_opt.log | tail -3`，三重态理论值 $S(S+1)=2.0$，明显偏离（如 > 2.1）时需警惕（见 [第13章](13-openshell-stability.md)）；
- 绝热 S-T 能隙 $\Delta E_{\rm ST}^{\rm adia}=E(S1,R_{S1})-E(T1,R_{T1})$，垂直磷光能 $\approx E(T1,R_{T1})-E(S0,R_{T1})$。

### 7.2 途径 B：TD 三重态

以闭壳层 S0 为参考，用 TD 求三重态并优化：

```text
#p CAM-B3LYP/def2SVP TDA=(Triplets,NStates=4,Root=1) Opt SCRF=(SMD,Solvent=Ethanol) Geom=Check Guess=Read NoSymm
```

- TD 求高阶三重态（T2、T3……）是 ΔSCF 做不到的，适合研究 ISC 通道（S1 → T_n 能量匹配）；
- 范围分离泛函配合 TD 处理三重态时容易出现“三重态不稳定性”导致 T1 能量偏低，**三重态推荐用 TDA**；
- 需要 T1 的零点能时，可在优化结构上用 `TD=(Triplets,NStates=4,Root=1) Freq` 做解析 TD 频率，或直接采用途径 A（UKS 频率）；
- 想一次得到垂直 S/T 能级：`TD=(50-50,NStates=5)` 给出 5 个单重态和 5 个三重态。

### 7.3 旋轨耦合（SOC）与磷光强度

Gaussian 16 的 TD-DFT **不计算旋轨耦合矩阵元**，因此三重态的振子强度恒为 0，无法直接得到磷光速率或 ISC 速率。常用做法：

- [PySOC](https://github.com/gaox-qd/pysoc)：读取 Gaussian TD 结果计算 S-T 间 SOC 矩阵元，常用于 TADF 与 ISC 研究（需按其说明准备特定格式的输出）；
- 换用支持 SOC-TDDFT 的程序（如 ORCA、BDF 等）做单点。

## 8. 实战四：ECD 与绝对构型

### 8.1 推荐流程

1. 构象搜索：柔性分子用 [CREST](https://github.com/crest-lab/crest)（GFN2-xTB）搜索，保留 3–5 kcal/mol 窗口内的构象；
2. DFT 优化与频率：如 B3LYP-D3(BJ)/def2-SVP 或 def2-TZVP（含溶剂），确认无虚频，得到每个构象的 G；
3. TD 计算：CAM-B3LYP 或 ωB97X-D，带弥散函数的基组，`NStates` 覆盖实验谱区；
4. 每个构象用 R(velocity) 展宽，按 Boltzmann 权重加和；
5. 与实验谱比较（必要时统一平移），符号与形状一致则归属该构型，相反则为其对映体；
6. 用 UV 谱交叉验证激发能平移量，或用旋光度（`Polar=OptRot`，见 [第11章](11-nmr-props.md)）、VCD 做独立佐证。

### 8.2 完整输入：(R)-甲基环氧乙烷

示例目录：[10-ecd-methyloxirane](../examples/10-ecd-methyloxirane/README.md)。甲基环氧乙烷（propylene oxide）只有 10 个原子、没有构象自由度，是 ECD/旋光计算的经典测试分子，其低能激发带有明显的 Rydberg 特征，必须使用弥散函数。

```text
! Chapter 10 example: ECD of (R)-methyloxirane (propylene oxide)
! Step 1: opt+freq   Step 2: TD-DFT with a diffuse basis (low-lying states are Rydberg-like)
! Rotatory strengths are printed automatically by every TD job; use R(velocity)
%chk=mox_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt Freq

(R)-methyloxirane opt+freq, B3LYP-D3(BJ)/def2-TZVP, gas phase

0 1
 C      -1.213067     0.014264    -0.254999
 C       0.044295    -0.211739     0.505216
 C       1.379921     0.199319    -0.052469
 O       0.932783    -1.165198    -0.086882
 H      -1.642969     0.984900     0.010350
 H      -1.048074    -0.000908    -1.337365
 H      -1.945555    -0.761414    -0.011385
 H      -0.096051    -0.241705     1.576914
 H       1.415290     0.737330    -0.991227
 H       2.173427     0.445152     0.641848

--Link1--
%oldchk=mox_opt.chk
%chk=mox_ecd.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/Aug-cc-pVDZ TD=(NStates=20) Geom=Check Guess=Read

(R)-methyloxirane ECD/UV, 20 singlets, CAM-B3LYP/aug-cc-pVDZ, gas phase

0 1

```

> [!IMPORTANT]
> 初始坐标的手性必须正确：本例坐标由 SMILES `C[C@@H]1CO1` 生成，经 RDKit 的 CIP 判定为 R 构型（手性中心为 2 号 C）。在 GaussView 中可以用“Symmetrize/Mirror”或直接把所有 z 坐标取反得到对映体。几何优化不会改变手性，但输入错一个原子位置就可能得到另一个构型。

### 8.3 读取旋转强度

每个 TD 作业都会自动打印旋转强度，在 `Excited State` 汇总之前（示意输出，数值仅供参考）：

```text
 <0|del|b> * <b|rxdel|0> (Au), Rotatory Strengths (R) in
 cgs (10**-40 erg-esu-cm/Gauss)
       state          XX          YY          ZZ    R(velocity)    E-M Angle
         1        -2.1034     -0.8712      1.0261     -0.6495      91.27
         2         3.5520      1.1203     -0.4410      1.4104      87.95
 <0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R) in
 cgs (10**-40 erg-esu-cm/Gauss)
       state          XX          YY          ZZ     R(length)
         1        -2.3011     -0.9020      1.1124     -0.6969
```

```bash
grep -A 25 "R(velocity)" mox_ecd.log | head -30
grep -A 25 "R(length)"   mox_ecd.log | head -30
```

R(velocity) 与 R(length) 差别很大时，说明基组不够完备，应增大基组（例如 aug-cc-pVTZ）再确认结论。R 的单位是 $10^{-40}$ erg·esu·cm/Gauss（cgs）。

### 8.4 ECD 展宽、Boltzmann 加权与对映体

常用的高斯展宽公式（Stephens & Harada 2010，E、$\Delta E_i$、σ 均以 eV 为单位，$R_i$ 以 cgs 单位代入）：

$$
\Delta\varepsilon(E)=\frac{1}{2.296\times10^{-39}\sqrt{\pi}\,\sigma}\sum_i \Delta E_i\,R_i\,\exp\!\left[-\left(\frac{E-\Delta E_i}{\sigma}\right)^2\right]
$$

多构象 Boltzmann 权重（通常用含热校正的 G，必要时用高水平单点修正电子能）：

$$
w_k=\frac{\exp(-\Delta G_k/RT)}{\sum_j \exp(-\Delta G_j/RT)},\qquad \Delta\varepsilon_{\rm tot}(E)=\sum_k w_k\,\Delta\varepsilon_k(E)
$$

298.15 K 时 $RT=0.5925$ kcal/mol，1 kcal/mol 的能量差对应约 5.4 倍的权重比，因此构象能量的误差会直接影响谱形。

- 对映体的 ECD 谱严格镜像：把所有 R 取负即可，无需重算（`spectrum_uvvis_ecd.py` 支持对映体翻转，见 `python ../scripts/spectrum_uvvis_ecd.py -h`）；
- [SpecDis](https://github.com/TBruhn74/SpecDis) 可做 Boltzmann 加权、UV 平移和与实验谱的相似度因子计算，是 ECD 归属的常用工具；
- [AaronTools.py](https://github.com/QChASM/AaronTools.py) 也提供 UV-Vis/ECD 绘图与构象加权；
- Multiwfn 主功能 11 中选择 ECD 类型，同样可以对多个体系按权重叠加。

> [!CAUTION]
> 只用一个构象、或者只看一个峰的符号就下结论，是 ECD 定构型最常见的错误。Pescitelli 与 Bruhn 的“good computational practice”（Chirality 2016）建议：完整的构象搜索、溶剂模型、至少两种泛函交叉验证、整条谱与实验比较，并报告平移量和半宽。

## 9. 激发态性质分析：NTO、DCT 与空穴-电子

示例目录：[10-nto-dct-dmabn](../examples/10-nto-dct-dmabn/README.md)，体系为经典的推-拉分子 4-(二甲氨基)苯甲腈（DMABN）。

### 9.1 自然跃迁轨道 NTO

把跃迁密度矩阵 $\mathbf T$ 做奇异值分解，得到成对的“空穴” NTO 与“粒子” NTO 及本征值 $\lambda_k$（$\sum_k\lambda_k=1$）。通常一对 NTO 的 λ > 0.9，就能用一张“空穴 → 电子”图描述整个激发，比列出十几个 MO 组态直观得多。

Gaussian 的做法（官方 [Population 页](https://gaussian.com/population/) 的 `NTO`/`SaveNTO` 选项）：先做 TD 计算，再用 `Guess=Only` 步骤读取第 N 个态的跃迁密度生成 NTO，并写入**单独的 chk 副本**（NTO 会替换 chk 中的正则轨道）：

```text
--Link1--
%oldchk=dmabn_td.chk
%chk=dmabn_nto_s1.chk
%mem=8GB
%nprocshared=8
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=1) Pop=(Minimal,NTO,SaveNTO)

```

- `Density=(Check,Transition=1)`：使用 chk 中基态到第 1 激发态的跃迁密度；`Density=Check` 隐含 `Guess=Only ChkBasis`，不会重新做 SCF；
- `Geom=AllCheck`：标题、电荷多重度和结构全部从 chk 读取，所以这一步只有 Link0、route 和一个空行；
- 每个态各写一个 chk（`dmabn_nto_s1.chk`、`dmabn_nto_s2.chk`），然后 `formchk` 转成 fchk。在 GaussView 或 Multiwfn 中，最高的“占据”NTO 为空穴，最低的“空”NTO 为粒子，λ 存放在原本的轨道能量位置。

```bash
formchk dmabn_nto_s2.chk dmabn_nto_s2.fchk
cubegen 0 MO=HOMO dmabn_nto_s2.fchk nto_s2_hole.cub 0 h
cubegen 0 MO=LUMO dmabn_nto_s2.fchk nto_s2_particle.cub 0 h
```

cubegen 与 Multiwfn 的批量用法见 [第16章](16-postprocess-multiwfn.md)；仓库脚本 [`batch_cube_multiwfn.sh`](../scripts/batch_cube_multiwfn.sh) 可批量生成 cube（`bash ../scripts/batch_cube_multiwfn.sh -h`）。

### 9.2 DCT 电荷转移指标

`Pop=DCT`（Gaussian 16 新增功能）计算 Le Bahers–Adamo–Ciofini 的 $D_{\rm CT}$ 指标：比较激发态与基态密度差中电子增加区和减少区的重心距离 $D_{\rm CT}$（Å）、转移电荷量 $q_{\rm CT}$ 以及对应的偶极变化 $\mu_{\rm CT}$。按官方说明，它对作业中可用的未弛豫激发态密度计算，如果请求了弛豫密度（`Density=Current`），也会对弛豫密度计算。示例的第 5 步仿照官方 Density 页的写法，对 S2 计算弛豫密度并给出 DCT：

```text
--Link1--
%oldchk=dmabn_td.chk
%chk=dmabn_dct_s2.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/6-31+G(d) TD=(Read,NStates=6,Root=2) SCRF=(SMD,Solvent=Acetonitrile)
   Geom=AllCheck Guess=Read Density=Current Pop=DCT

```

$D_{\rm CT}$ 没有普适的“CT 阈值”，应与分子尺寸（如 D、A 中心距离）对比，或在同一体系的不同态之间比较。$D_{\rm CT}$ 很大时，必须用范围分离泛函，并考虑态特定溶剂化（第 6.3 节）。

> [!NOTE]
> DCT 输出块的具体格式随版本略有不同，第一次使用时建议用 `grep -i -A12 "DCT\|charge transfer" dmabn.log` 定位后再写解析脚本。

### 9.3 Multiwfn 空穴-电子分析

Multiwfn 主功能 18（电子激发分析）读取 **fchk + TD 输出 log**：

1. TD 作业的 route 加 `IOp(9/40=4)`（Multiwfn 手册要求，以打印绝对值 ≥ $10^{-4}$ 的组态系数；大体系可用 3）；
2. `formchk dmabn_td.chk dmabn_td.fchk`；
3. `Multiwfn dmabn_td.fchk` → 主功能 18 → 子功能 1（空穴-电子分析），输入 log 路径和激发态序号；
4. 读取 D 指数（空穴与电子重心距离，Å）、Sr（重叠程度）、H、t 指数（t > 0 通常表示空穴与电子明显分离）；
5. 子功能 4 给出 Δr 指数，子功能 6 也可生成 NTO。

批量处理可用 [BaneWfn](https://github.com/bane-dysta/banewfn) 等工具把菜单操作模板化。Δr 与 Sr 的经验判据见 Multiwfn 手册 3.21 节与对应原始文献。

### 9.4 判断 CT 的组合拳

- NTO：空穴和粒子是否位于分子的不同片段；
- 定量指标：$D_{\rm CT}$、D 指数、Δr、t 指数；
- 泛函依赖性：同一个态在 B3LYP 与 CAM-B3LYP 下能量相差 > 0.5 eV，往往说明它有明显 CT 成分；
- 溶剂依赖性：极性溶剂中明显红移的态通常具有 CT 特征。

[ELOTE](https://github.com/joaquinbarroso/ELOTE) 也可以分析 Gaussian TD 激发的电子定域特征，辅助判断 CT。

## 10. 溶剂中的激发态：小结

| 场景 | 推荐写法 | 默认溶剂化 | 备注 |
|---|---|---|---|
| 垂直吸收 | `TD SCRF=(SMD 或 IEFPCM,Solvent=...)` | 非平衡 LR | 最常用 |
| 吸收，CT 态 | 加 `CorrectedLR`；或 `NonEquilibrium=Save` 基态 + `ExternalIteration,NonEquilibrium=Read` | 非平衡 | 见官方 SCRF 页第 3 步 |
| 激发态优化 | `TD Opt SCRF=(...)` | 平衡 LR | 可加 Freq |
| 垂直发射 | cLR/EI 激发态（`NonEquilibrium=Save`）+ 基态 `NonEquilibrium=Read` | 激发态平衡、基态非平衡 | 第 6.3 节 |
| 溶剂化自由能类问题 | SMD | — | 见 [第9章](09-solvation.md) |

隐式溶剂无法描述特异性相互作用（氢键、π-π 堆积）。若溶质与溶剂有强氢键（如质子性溶剂中的羰基 n→π*），可加入 1–2 个显式溶剂分子再套 PCM，或用 [第15章](15-oniom.md) 的 ONIOM/QM-MM 方法。

## 11. 外部工具与扩展

- [cclib](https://github.com/cclib/cclib)：解析 TD 输出中的激发能、振子强度、旋转强度，适合批量处理；
- [LeoX](https://github.com/LeonardoESousa/LeoX)：基于 Gaussian TD-DFT 的核系综法（nuclear ensemble）吸收/发射谱，以及长程校正参数 ω 的最优调节；
- [PySOC](https://github.com/gaox-qd/pysoc)：Gaussian 之外补算 S-T 旋轨耦合；
- [SpecDis](https://github.com/TBruhn74/SpecDis)：ECD/UV 谱的 Boltzmann 加权与相似度因子；
- [AaronTools.py](https://github.com/QChASM/AaronTools.py)：UV-Vis/ECD 绘图、构象处理；
- [gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：中文项目，包含 TD-DFT 的输入模板与结果自查项，可对照本章的检查清单；
- [Multiwfn](http://sobereva.com/multiwfn/)：空穴-电子分析、谱图绘制、跃迁偶极分解。

## 参数选择与收敛建议

| 项目 | 推荐 | 原因 / 说明 |
|---|---|---|
| 泛函（价层 ππ*、nπ*） | PBE0、M06-2X、CAM-B3LYP | 至少比较两种泛函，看结论是否稳健 |
| 泛函（CT、Rydberg、ECD） | CAM-B3LYP、ωB97X-D、LC-ωhPBE | 长程交换校正 |
| 基组（单点） | def2-TZVP、6-311+G(d,p) | Rydberg/阴离子/ECD 用 aug-cc-pVDZ 及以上 |
| 基组（激发态优化） | def2-SVP、6-31+G(d) | TD 梯度和频率较贵 |
| NStates | 覆盖到谱区下限再多 0.5 eV；优化时 ≥ Root+2 | 防止谱带截断与根翻转 |
| TD 与 TDA | 光谱用 TD；三重态、优化困难时用 TDA | TDA 避免三重态不稳定性 |
| 展宽 | 高斯，FWHM 0.3–0.5 eV（UV），0.2–0.4 eV（ECD） | 必须注明；在能量空间展宽 |
| 溶剂 | 单点非平衡、优化平衡（默认） | CT 态用 cLR/EI |
| 组态系数打印 | `IOp(9/40=3)` 或 `IOp(9/40=4)` | Multiwfn 分析需要 4 |
| 激发态优化 | `NoSymm`，必要时 `Opt=CalcFC` | 防止被困在对称鞍点 |
| ECD 构象 | 3–5 kcal/mol 窗口内全部构象，用 G 加权 | 单构象结论不可靠 |

## 常见问题与报错

**1. 光谱在某个波长处“戛然而止”或缺少实验中的强峰**
原因：`NStates` 太少，高能激发态没有算到。
解决：检查最后一个 `Excited State` 的波长；用 `TD=(Read,Add=10)` 追加态，或增大 `NStates` 重算。

**2. B3LYP 算出一堆 f≈0 的低能态，实验中没有对应吸收**
原因：虚假的低能 CT 态（常规杂化泛函低估 CT 能量），常见于 D-A 分子、聚集体、配合物。
解决：换用 CAM-B3LYP、ωB97X-D 重算并比较；用 NTO 或 DCT 确认这些态的 CT 性质。

**3. 激发能出现负值或极小值（< 0.1 eV），或 TD 迭代难以收敛**
原因：基态波函数不稳定（实际存在更低的解，常见于双自由基、拉伸键、某些金属配合物），TD 的 A−B 矩阵出现负本征值。
解决：先做 `Stable=Opt` 得到稳定波函数（[第13章](13-openshell-stability.md)），或改用 TDA；对三重态尤其推荐 TDA。

**4. S1 优化后出现虚频，或收敛到对称结构但能量偏高**
原因：对称性约束使优化停留在鞍点上（例如平面分子的 nπ* 态）。
解决：加 `NoSymm`，沿虚频方向扰动结构（可用 [pyQRC](https://github.com/patonlab/pyQRC)）后重新优化，必要时加 `Opt=CalcFC`。

**5. S1 优化过程中振子强度或主要组态突然改变**
原因：根翻转，`Root=1` 跟随了交叉后的另一个态。
解决：见第 6.2 节 WARNING；增加 `NStates`，比较前后 NTO，必要时改用 `Root=2` 或从不同初始结构重新优化。

**6. 在 `Opt` 作业中加入 `SCRF=ExternalIteration` 后作业失败**
原因：官方 SCRF 页说明 ExternalIteration 只能用于能量计算。
解决：激发态结构用 LR-PCM 优化，再在优化结构上做 cLR/EI 单点（第 6.3 节）。

**7. `QPErr --- A syntax error was detected in the input line.`**
原因：route 拼写错误，例如 `TD=(NStaes=10)`、括号不配对、泛函名写成 `CAMB3LYP`。
解决：根据 Gaussian 在 route 下方打印的标记定位错误；常用正确写法为 `CAM-B3LYP`、`wB97XD`、`M062X`、`PBE1PBE`。

**8. 用 G96、P86、PKZB、wPBEh、PBEh 等泛函做 TD 优化或频率时报错退出**
原因：官方说明这些泛函缺少 TD 梯度/频率所需的高阶导数。
解决：换用文中表格里的泛函。

**9. `galloc: could not allocate memory.` 或作业被调度系统杀死**
原因：`%mem` 超过节点或作业申请的内存；TD 频率、大 `NStates` 时内存需求明显增加。
解决：`%mem` 设为作业申请内存的 75–85%（[第17章](17-hpc-jobs.md)）；反过来，`%mem` 太小会导致 TD 频率非常慢，官方 C.01 发行说明专门提到增大 `%Mem` 能显著提速。

**10. NTO 步骤出错或得到的轨道看起来是正则 MO**
原因：`%oldchk` 指错了文件（不是 TD 作业的 chk），或 `Transition=N` 大于 `NStates`，或把 NTO 直接写回了原 chk 导致后续步骤读到 NTO。
解决：严格按“`%oldchk`=TD 的 chk、`%chk`=新副本”书写；每个态各用一个 chk。

**11. 开壳层体系激发态的 `<S**2>` 远离理论值**
原因：UKS 参考的自旋污染在激发态中进一步放大，部分态实际是高自旋态的混合。
解决：剔除 `<S**2>` 严重偏离的态，或改用 ROKS 参考、多参考方法；参见 [第13章](13-openshell-stability.md)。

更多报错与排查思路见 [第18章](18-errors-capstone.md) 和 [Solutions for Gaussian Error Messages](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/)。

## 练习

1. 用 B3LYP、PBE0、CAM-B3LYP、ωB97X-D 四种泛函计算香豆素前 10 个单重态，比较 S1、S2 的能量与 f，并画出四条 UV-Vis 谱（FWHM 0.4 eV）。
2. 把 `coumarin_td_camb3lyp_smd.gjf` 第 2 步改成 `TDA=(NStates=20)`，比较激发能与振子强度的变化。
3. 用 S1 优化结果计算香豆素在乙醇中的 $E_{\rm adia}$、$E_{00}$ 和 Stokes 位移（cm⁻¹），并比较 LR 与 cLR 的垂直发射能差别。
4. 把甲基环氧乙烷的所有 z 坐标取反，重算 ECD，验证两条谱互为镜像；再把基组换成 def2-SVP（无弥散），观察 R(velocity) 与 R(length) 的差别。
5. 对 DMABN 的 S1、S2 分别做 NTO 和 Multiwfn 空穴-电子分析，判断哪个是 LE 态、哪个更具 CT 特征；改用气相重算，观察态顺序是否变化。

## 参考资料

官方文档：

- [TD 关键词](https://gaussian.com/td/)、[Density 关键词](https://gaussian.com/density/)、[Population 关键词（NTO、DCT）](https://gaussian.com/population/)
- [SCRF 关键词（非平衡溶剂化、CorrectedLR、ExternalIteration、荧光示例）](https://gaussian.com/scrf/)
- [Creating UV/Visible Plots from Gaussian results](https://gaussian.com/uvvisplot/)
- [Freq 关键词（FC/HT/FCHT 振动分辨光谱）](https://gaussian.com/freq/)
- [IOp 第 9 组说明（9/40 打印阈值）](https://gaussian.com/overlay9/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)
- [Multiwfn 主页](http://sobereva.com/multiwfn/) 与 [在线手册](http://sobereva.com/multiwfn/manual/)

相关 GitHub 仓库（仅链接）：

- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian 激发态数据
- [LeonardoESousa/LeoX](https://github.com/LeonardoESousa/LeoX)：核系综法光谱、ω 最优调节
- [joaquinbarroso/ELOTE](https://github.com/joaquinbarroso/ELOTE)：TD 激发的定域/CT 特征分析
- [gaox-qd/pysoc](https://github.com/gaox-qd/pysoc)：基于 Gaussian TD 结果计算 SOC
- [TBruhn74/SpecDis](https://github.com/TBruhn74/SpecDis)：ECD/UV 谱处理与相似度因子
- [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py)：UV-Vis/ECD 绘图与构象工具
- [crest-lab/crest](https://github.com/crest-lab/crest)：ECD 前的构象搜索
- [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn)：Multiwfn 批处理
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：含 TD-DFT 模板与自查清单

关键文献：

- M. E. Casida, in *Recent Advances in Density Functional Methods*, World Scientific, 1995（线性响应 TD-DFT）
- T. Yanai, D. P. Tew, N. C. Handy, *Chem. Phys. Lett.* 393, 51 (2004). DOI: [10.1016/j.cplett.2004.06.011](https://doi.org/10.1016/j.cplett.2004.06.011)（CAM-B3LYP）
- J.-D. Chai, M. Head-Gordon, *Phys. Chem. Chem. Phys.* 10, 6615 (2008). DOI: [10.1039/B810189B](https://doi.org/10.1039/B810189B)（ωB97X-D）
- C. Adamo, D. Jacquemin, *Chem. Soc. Rev.* 42, 845 (2013). DOI: [10.1039/C2CS35394F](https://doi.org/10.1039/C2CS35394F)
- A. D. Laurent, D. Jacquemin, *Int. J. Quantum Chem.* 113, 2019 (2013). DOI: [10.1002/qua.24438](https://doi.org/10.1002/qua.24438)
- M. Caricato et al., *J. Chem. Phys.* 124, 124520 (2006). DOI: [10.1063/1.2183309](https://doi.org/10.1063/1.2183309)（cLR）
- R. Improta et al., *J. Chem. Phys.* 125, 054103 (2006). DOI: [10.1063/1.2222364](https://doi.org/10.1063/1.2222364)（态特定 PCM）
- R. L. Martin, *J. Chem. Phys.* 118, 4775 (2003). DOI: [10.1063/1.1558471](https://doi.org/10.1063/1.1558471)（NTO）
- T. Le Bahers, C. Adamo, I. Ciofini, *J. Chem. Theory Comput.* 7, 2498 (2011). DOI: [10.1021/ct200308m](https://doi.org/10.1021/ct200308m)（DCT）
- P. J. Stephens, N. Harada, *Chirality* 22, 229 (2010). DOI: [10.1002/chir.20733](https://doi.org/10.1002/chir.20733)（ECD 计算与展宽）
- G. Pescitelli, T. Bruhn, *Chirality* 28, 466 (2016). DOI: [10.1002/chir.22600](https://doi.org/10.1002/chir.22600)（ECD 定构型的良好实践）
- T. Lu, F. Chen, *J. Comput. Chem.* 33, 580 (2012). DOI: [10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)（Multiwfn）
- X. Gao et al., *J. Chem. Theory Comput.* 13, 515 (2017). DOI: [10.1021/acs.jctc.6b00915](https://doi.org/10.1021/acs.jctc.6b00915)（PySOC）

[← 上一章](09-solvation.md) | [返回目录](../README.md) | [下一章 →](11-nmr-props.md)
