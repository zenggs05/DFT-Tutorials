# 第13章 开壳层体系、波函数稳定性与 Broken-Symmetry（Open-Shell Systems, Wavefunction Stability & Broken-Symmetry DFT in Gaussian）

本章讲 Gaussian 16 中开壳层体系的完整处理流程：UHF/UKS 与 ROHF/ROKS 怎么选（`UB3LYP`、`ROB3LYP`），如何读 `<S**2>` 判断自旋污染，以及如何用 `Stable`/`Stable=Opt` 检查并修复波函数不稳定性。在此基础上，用 `Guess=Mix` 处理单重态双自由基，用 `Guess=Fragment` 构造 broken-symmetry（BS）态，并按 Yamaguchi 公式计算磁耦合常数 J（open-shell singlet, spin contamination, wavefunction stability, broken-symmetry DFT, exchange coupling）。

---

## 1. 本章目标与适用场景

学完本章你应能够：

- 根据电子数和化学直觉写出正确的电荷与自旋多重度，并比较高自旋、低自旋态的能量；
- 区分 R/U/RO 三类波函数，知道 Gaussian 在不写前缀时默认用哪一种；
- 从输出中读出 `<S**2>`（湮灭前/后）并判断自旋污染是否可以接受；
- 用 `Stable` 发现 RHF→UHF、内部不稳定性，用 `Stable=Opt` 得到稳定解，再用 `Guess=Read` 接续后续计算；
- 用 `Guess=Mix` 或 `Stable=Opt` 得到单重态双自由基的 broken-symmetry 解，并用近似自旋投影修正能量；
- 用 `Guess=Fragment=N` 和负多重度构造反铁磁耦合的初猜，计算 J 并正确说明所用哈密顿量的约定；
- 知道单行列式 DFT 什么时候不再可靠，以及如何借助 MOKIT 把轨道交给 CASSCF/NEVPT2。

本章覆盖的典型体系：

| 体系类型 | 例子 | 关键问题 | 本章对应小节 |
|---|---|---|---|
| 二重态自由基 | CH3·、NO、苯氧自由基、自由基阳离子 | 自旋污染、ROKS 对照 | §5、§6 |
| 三重态基态分子 | O2、CH2（卡宾）、NH | 多重度设定、单重态-三重态能隙 | §4、§5 |
| 单重态双自由基 | 对苯炔（p-benzyne）、拉长的键、某些多环芳烃 | RHF→UHF 不稳定性、`Guess=Mix` | §7、§8 |
| 过渡金属多自旋态 | Fe(II)/Fe(III)、Co(II)、Mn(II) 配合物 | 高/低自旋能量排序、泛函依赖 | §5.4 |
| 反铁磁耦合多核配合物 | 双核 Cu(II)、Mn2O2、Fe2S2 | Broken-symmetry、J 的计算 | §9 |
| 键均裂/解离曲线 | H2、C–C 键拉长 | 单重态解离到两个自由基 | §7、§8 |

> [!NOTE]
> 本章三个示例目录：[../examples/13-spin-states-o2-ch2/README.md](../examples/13-spin-states-o2-ch2/README.md)（O2、CH2 多重度比较）、[../examples/13-diradical-pbenzyne/README.md](../examples/13-diradical-pbenzyne/README.md)（双自由基与稳定性）、[../examples/13-bs-cu2-dimer/README.md](../examples/13-bs-cu2-dimer/README.md)（broken-symmetry 与 J）。SCF 不收敛的通用对策见 [第4章](04-sp-scf.md)，输入中碎片与电荷多重度行的写法见 [第2章](02-input-anatomy.md)。

---

## 2. 原理速览

### 2.1 自旋多重度与电子数奇偶

总自旋量子数 S 与自旋多重度的关系是 $2S+1$。Gaussian 的单行列式计算取 $M_S = S$，即 α 电子比 β 电子多 $2S$ 个：

$$
N_\alpha - N_\beta = 2S = (2S+1) - 1,\qquad N_\alpha + N_\beta = N_e = \sum_A Z_A - Q
$$

所以电子数为奇数时只能是二重态、四重态……，为偶数时只能是单重态、三重态……。写错时 Gaussian 在 L301 直接报错（见“常见问题”）。

| 体系 | 电子数 | 常见多重度 | 说明 |
|---|---|---|---|
| CH3·（甲基自由基） | 9 | 2 | 一个未成对电子 |
| O2 | 16 | 3（基态）/ 1（激发态 a¹Δg） | 两个 π* 轨道各一个电子 |
| CH2（亚甲基卡宾） | 8 | 3（基态 ³B1）/ 1（¹A1） | 单重态-三重态能隙是经典基准 |
| Cu(II) 配合物（d⁹） | 奇数 | 2 | 双核时可为 1 或 3 |
| Fe(III) 配合物（d⁵） | 奇数 | 2 / 4 / 6 | 低、中、高自旋都可能 |
| Fe(II) 配合物（d⁶） | 偶数 | 1 / 3 / 5 | 自旋交叉体系的典型 |

### 2.2 R、U、RO 三类单行列式波函数

| 类型 | Gaussian 前缀 | 空间轨道 | ⟨S²⟩ | 优点 | 缺点 |
|---|---|---|---|---|---|
| 限制性（RHF/RKS） | `R`，如 `RB3LYP` | α、β 共用 | 严格为 0 | 闭壳层的标准做法，便宜 | 不能描述开壳层；键断裂时给出过高能量 |
| 非限制性（UHF/UKS） | `U`，如 `UB3LYP` | α、β 各自优化 | 可能偏离 $S(S+1)$ | 能量更低，可描述自旋极化和双自由基 | 自旋污染；不是 $\hat S^2$ 本征态 |
| 限制性开壳层（ROHF/ROKS） | `RO`，如 `ROB3LYP` | 双占据部分共用 | 严格等于 $S(S+1)$ | 无自旋污染 | 缺少自旋极化，自旋密度分布和超精细常数不可靠；后续性质支持较少 |

> [!IMPORTANT]
> Gaussian 的默认规则：**多重度为 1 时默认限制性（R），多重度大于 1 时默认非限制性（U）**。因此三重态写 `B3LYP` 就等于 `UB3LYP`；而开壳层**单重态**（双自由基、BS 态）必须显式写 `UB3LYP`，否则即使读入了 UHF 轨道也会按 RKS 计算。

### 2.3 ⟨S²⟩ 与自旋污染

对于单行列式 UHF，$\langle \hat S^2\rangle$ 可以由占据轨道的 α/β 重叠算出：

$$
\langle \hat S^2\rangle_{\rm UHF} = S_z(S_z+1) + N_\beta - \sum_{i}^{N_\alpha}\sum_{j}^{N_\beta}\left|\langle \phi_i^\alpha|\phi_j^\beta\rangle\right|^2
$$

α、β 空间轨道完全相同时，最后一项正好抵消 $N_\beta$，得到纯自旋态 $S(S+1)$；两者差别越大，混入的高自旋成分越多，这就是**自旋污染**。对 UKS，Gaussian 用 Kohn-Sham 行列式按同一公式计算，它不是真实体系的 ⟨S²⟩，但仍是很有用的诊断量。

| 多重度 | S | 理想 $S(S+1)$ | 约 10% 警戒线 |
|---|---|---|---|
| 1（闭壳层） | 0 | 0.000 | — |
| 2 | 1/2 | 0.750 | > 0.82 |
| 3 | 1 | 2.000 | > 2.2 |
| 4 | 3/2 | 3.750 | > 4.1 |
| 5 | 2 | 6.000 | > 6.6 |
| 6 | 5/2 | 8.750 | > 9.6 |

Gaussian 还会打印“湮灭第一个自旋污染成分”（annihilation of the first spin contaminant，即去掉 S+1 成分）后的值。**湮灭后的值只用于诊断，能量并没有随之修正。**湮灭前后差别很大，说明污染主要来自更高一级的自旋态。

### 2.4 波函数稳定性

SCF 收敛只说明找到了能量对轨道变化的**驻点**，不保证是极小点。稳定性分析计算轨道旋转的二阶导数矩阵（稳定性矩阵）的最低本征值：若有负本征值，沿该方向旋转轨道能量还会下降，即波函数**不稳定**。Gaussian 的 `Stable` 可以检验（[Seeger & Pople 1977](https://doi.org/10.1063/1.434318)；[Bauernschmitt & Ahlrichs 1996](https://doi.org/10.1063/1.471637)）：

| 不稳定性类型 | 含义 | 常见体系 |
|---|---|---|
| 内部不稳定（internal） | 在同样约束下存在能量更低的解，当前解是“错误的态” | 过渡金属、收敛到激发组态 |
| RHF→UHF（外部，triplet 型） | 放开 α/β 相同的约束能量下降 | 双自由基、拉长的键、O3、反铁磁体系 |
| 实→复（外部） | 允许复数轨道能量下降 | 少见；Gaussian 默认不检查 |
| 对称性降低 | 轨道对称性约束导致的高能解 | 高对称结构 |

### 2.5 开壳层单重态、broken-symmetry 与自旋投影

两个弱耦合的未成对电子可以构成单重态（S=0）或三重态（S=1）。真正的开壳层单重态至少需要两个行列式：

$$
{}^1\Psi = \tfrac{1}{\sqrt2}\left(|\cdots a^\alpha b^\beta| - |\cdots a^\beta b^\alpha|\right)
$$

单行列式 UKS 只能取其中一个 $|\cdots a^\alpha b^\beta|$，这就是 **broken-symmetry（BS）态**：它是单重态和三重态（$M_S=0$ 分量）约 1:1 的混合，⟨S²⟩ 接近 1。BS 能量通常已经比闭壳层 RKS 合理得多，但仍被三重态“污染”。Yamaguchi 的近似自旋投影（approximate projection, AP）假设污染只来自三重态，修正后的单重态能量为（[Yamaguchi et al. 1988](https://doi.org/10.1016/0009-2614(88)80378-6)）：

$$
E_{\rm S}^{\rm AP} = E_{\rm BS} + \alpha\,(E_{\rm BS} - E_{\rm T}),\qquad
\alpha = \frac{\langle S^2\rangle_{\rm BS}}{\langle S^2\rangle_{\rm T} - \langle S^2\rangle_{\rm BS}}
$$

当 $\langle S^2\rangle_{\rm BS}\approx 1$、$\langle S^2\rangle_{\rm T}\approx 2$ 时，它退化为 Ziegler-Rauk-Baerends 的“求和法”$E_{\rm S} \approx 2E_{\rm BS} - E_{\rm T}$（[Ziegler et al. 1977](https://doi.org/10.1007/BF00551551)）。

### 2.6 Heisenberg 模型与磁耦合常数 J

两个磁中心 A、B（自旋 $S_A$、$S_B$）之间的交换作用用 Heisenberg-Dirac-van Vleck 哈密顿量描述。本章统一采用

$$
\hat H = -2J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B
$$

此时 $E(S) = -J\,S(S+1) + {\rm const}$，**J < 0 为反铁磁（低自旋基态），J > 0 为铁磁**。两个 S=1/2 中心的单重态-三重态能隙为 $E_{\rm T} - E_{\rm S} = -2J$。

用高自旋（HS）和 BS 两次单行列式计算估算 J 的常用公式（$S_{\max}=S_A+S_B$）：

| 公式 | 表达式（$\hat H=-2J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$） | 适用情形 |
|---|---|---|
| Noodleman（弱重叠） | $J = \dfrac{E_{\rm BS}-E_{\rm HS}}{S_{\max}^2}$ | 磁性轨道重叠很小（[Noodleman 1981](https://doi.org/10.1063/1.440939)） |
| Ruiz/Bencini（强重叠） | $J = \dfrac{E_{\rm BS}-E_{\rm HS}}{S_{\max}(S_{\max}+1)}$ | 磁性轨道重叠大；也常被用作“不投影”的经验做法 |
| Yamaguchi | $J = \dfrac{E_{\rm BS}-E_{\rm HS}}{\langle S^2\rangle_{\rm HS}-\langle S^2\rangle_{\rm BS}}$ | 一般情形，自动在两种极限之间过渡（[Soda et al. 2000](https://doi.org/10.1016/S0009-2614(00)00166-4)） |

> [!WARNING]
> 文献中还广泛使用 $\hat H = -J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$ 和 $\hat H = +J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$ 等约定，数值相差 2 倍或符号相反。**报告 J 时必须写明哈密顿量**，与实验或文献比较前先统一约定。

---

## 3. Gaussian 开壳层相关关键词速查

| 关键词 / 写法 | 作用 | 备注 |
|---|---|---|
| `UB3LYP`、`UHF`、`UM062X` | 非限制性计算 | 多重度 > 1 时可省略 `U`；开壳层单重态必须写 |
| `ROB3LYP`、`ROHF` | 限制性开壳层 | 无自旋污染，用于对照或作为后续方法的参考 |
| `Stable` | 只做稳定性检查，不重新优化（默认 `NoOpt`） | 默认检查内部和实数外部不稳定性（`RExt`） |
| `Stable=Opt` | 发现不稳定就降低约束重新优化，直到稳定 | 从 RHF/RKS 出发时允许变为 UHF/UKS |
| `Guess=Mix` | 混合 HOMO 与 LUMO，打破 α-β 和空间对称性 | 用于得到单重态的 UHF/UKS（BS）解 |
| `Guess=(Fragment=N,Only)` | 由各碎片的初猜拼成整体初猜 | 碎片多重度为负表示其未成对电子取 β 自旋 |
| `Guess=Alter` | 交换指定轨道的占据 | 用于强制得到特定组态 |
| `Guess=Read` + `Geom=Check/AllCheck` | 从 chk 读入轨道/结构 | 接续 `Stable=Opt` 的标准做法 |
| `SCF=NoSymm` 或 `NoSymm` | 取消轨道对称性约束 | BS 解通常会破坏空间对称性 |
| `Pop=NaturalOrbitals`（`Pop=NO`） | 总密度的自然轨道分析（UHF/UKS 下即 UNO） | 占据数用来判断双自由基特征 |

---

## 4. 实战一：O2 的三重态与单重态（UKS、Stable、Guess=Mix）

### 4.1 为什么从 O2 开始

O2 基态是 $X\,{}^3\Sigma_g^-$，两个 π* 轨道各有一个 α 电子。它的最低单重态 $a\,{}^1\Delta_g$ 比基态高约 0.98 eV（约 22.6 kcal/mol；NIST Chemistry WebBook / Huber & Herzberg 光谱数据，近似值）。用三种单行列式描述同一几何下的“单重态”，可以直观看到 RKS、BS-UKS 和自旋投影的区别。

### 4.2 完整输入

文件 [../examples/13-spin-states-o2-ch2/o2_triplet_vs_singlet.gjf](../examples/13-spin-states-o2-ch2/o2_triplet_vs_singlet.gjf)：

```text
%nprocshared=8
%mem=8GB
%chk=o2_triplet.chk
#p UB3LYP/def2TZVP Opt Freq

O2 step 1: triplet ground state (X 3Sigma_g-), UKS opt + freq

0 3
O                   0.000000    0.000000    0.604000
O                   0.000000    0.000000   -0.604000

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=o2_triplet.chk
%chk=o2_singlet_rks.chk
#p RB3LYP/def2TZVP Stable Geom=Check

O2 step 2: closed-shell singlet (RKS) at the triplet geometry, stability test

0 1

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=o2_triplet.chk
%chk=o2_singlet_bs.chk
#p UB3LYP/def2TZVP Guess=Mix SCF=NoSymm Stable=Opt Geom=Check

O2 step 3: open-shell broken-symmetry singlet (Ms=0) from Guess=Mix

0 1

```

### 4.3 逐行解释

| 位置 | 写法 | 为什么这样写 |
|---|---|---|
| 第 1 步 | `UB3LYP ... 0 3` | 三重态，`U` 可省略但显式写出更清楚；`Opt Freq` 得到平衡键长并确认无虚频 |
| 第 2 步 | `RB3LYP ... Stable` + `0 1` | 闭壳层单重态 (π*)²，只做稳定性检查（不重新优化），看它是否存在 RHF→UHF 不稳定性 |
| 第 2、3 步 | `%oldchk` + `Geom=Check` | 读入第 1 步优化好的结构，各步写到不同的 chk，互不覆盖 |
| 第 3 步 | `UB3LYP Guess=Mix` + `0 1` | 显式 `U`；`Guess=Mix` 把 HOMO/LUMO 混合，得到 α、β 分别占据两个不同 π* 的 BS 解 |
| 第 3 步 | `SCF=NoSymm` | BS 解破坏 D∞h 的轨道对称性，取消约束避免 SCF 被“拉回”对称解 |
| 第 3 步 | `Stable=Opt` | 确认得到的 BS 解是稳定的；若不稳定会自动继续优化 |

### 4.4 运行与检查

```bash
cd Gaussian/examples/13-spin-states-o2-ch2
g16 < o2_triplet_vs_singlet.gjf > o2_triplet_vs_singlet.log
grep -c "Normal termination" o2_triplet_vs_singlet.log      # expect 3
grep "SCF Done" o2_triplet_vs_singlet.log | tail -3
grep "S\*\*2 before annihilation" o2_triplet_vs_singlet.log
grep -i "instability\|stable under" o2_triplet_vs_singlet.log
```

（示意输出，数值仅供参考）

```text
 SCF Done:  E(UB3LYP) =  -150.3xxxxxxx     a.u. after    8 cycles     <- step 1, last opt point
 S**2 before annihilation     2.0074,   after     2.0000
 The wavefunction has an RHF -> UHF instability.                      <- step 2
 SCF Done:  E(UB3LYP) =  -150.3xxxxxxx     a.u. after   19 cycles     <- step 3
 S**2 before annihilation     1.0068,   after     0.0xxx
 The wavefunction is stable under the perturbations considered.
```

判读要点：

1. 第 1 步 ⟨S²⟩ 应非常接近 2.00，说明三重态几乎没有污染；
2. 第 2 步出现 `RHF -> UHF instability`，说明闭壳层 (π*)² 不是这个几何下的最低单重态解；
3. 第 3 步 ⟨S²⟩ ≈ 1，正是“单重态 + 三重态 $M_S=0$ 分量”各半的 BS 特征，不是计算出错；
4. 能量顺序应为 $E_{\rm T} < E_{\rm BS} < E_{\rm RKS}$。

### 4.5 用自旋投影估算 a¹Δg 的激发能

$$
E({}^1\Delta_g) \approx 2E_{\rm BS} - E_{\rm T},\qquad
\Delta E_{\rm ST} = E({}^1\Delta_g) - E_{\rm T} \approx 2\,(E_{\rm BS}-E_{\rm T})
$$

更一般地用 §2.5 的 AP 公式（代入实际 ⟨S²⟩）。几点说明：

- 直接用 $E_{\rm BS}-E_{\rm T}$ 会把能隙**低估约一半**，因为 BS 态里一半是三重态；
- $E_{\rm RKS}-E_{\rm T}$ 通常**高估**，因为单个闭壳层行列式 (π_x*)² 不是纯的 ¹Δg（混有 ¹Σg⁺ 成分）；
- 比较实验值时注意：实验 0.98 eV 是绝热/光谱值，这里算的是三重态平衡几何下的近似值，两者键长差很小，可直接比较量级；不同泛函的结果相差可达零点几 eV。

---

## 5. 实战二：卡宾 CH2 的单重态-三重态能隙与 ROKS

### 5.1 输入

亚甲基 CH2 的基态是三重态 ³B1（H–C–H 角约 134°），最低单重态 ¹A1 是闭壳层（角约 102°），两者几何差别很大，必须**各自优化**。文件 [../examples/13-spin-states-o2-ch2/ch2_carbene_singlet_triplet.gjf](../examples/13-spin-states-o2-ch2/ch2_carbene_singlet_triplet.gjf)：

```text
%nprocshared=8
%mem=8GB
%chk=ch2_triplet.chk
#p UB3LYP/def2TZVP Opt Freq

CH2 step 1: triplet carbene (3B1), UKS opt + freq

0 3
C                   0.000000    0.000000    0.000000
H                   0.991000    0.000000   -0.421000
H                  -0.991000    0.000000   -0.421000

--Link1--
%nprocshared=8
%mem=8GB
%chk=ch2_singlet.chk
#p RB3LYP/def2TZVP Opt Freq

CH2 step 2: closed-shell singlet carbene (1A1), RKS opt + freq

0 1
C                   0.000000    0.000000    0.000000
H                   0.862600    0.000000   -0.698500
H                  -0.862600    0.000000   -0.698500

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=ch2_singlet.chk
%chk=ch2_singlet_stab.chk
#p RB3LYP/def2TZVP Stable Geom=AllCheck Guess=Read

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=ch2_triplet.chk
%chk=ch2_triplet_roks.chk
#p ROB3LYP/def2TZVP Geom=AllCheck Guess=Read

```

| 步骤 | 目的 | 要点 |
|---|---|---|
| 1 | 三重态 UKS 优化 + 频率 | 起始键角 134°，接近三重态平衡结构 |
| 2 | 单重态 RKS 优化 + 频率 | 起始键角 102°；¹A1 本质上是闭壳层，R 是合适的 |
| 3 | 单重态稳定性检查 | 确认闭壳层 ¹A1 在优化几何下稳定（若出现 RHF→UHF 不稳定，就要当双自由基处理） |
| 4 | 三重态 ROKS 单点 | 在 UKS 结构上比较 ROKS 与 UKS 的能量；`Geom=AllCheck` 省略标题和电荷多重度行 |

> [!TIP]
> `Geom=AllCheck` 从 chk 读入标题、电荷、多重度和结构，因此该步 route 之后**只需一个空行**。读入 UHF 轨道作为 ROHF 初猜时，Gaussian 会自动对 α、β 轨道做双正交化（见官方 `Guess` 关键词说明中的 `Biorthogonalize` 选项）。

### 5.2 结果处理

```bash
grep "SCF Done" ch2_carbene_singlet_triplet.log
grep "Sum of electronic and zero-point Energies" ch2_carbene_singlet_triplet.log
grep "S\*\*2 before annihilation" ch2_carbene_singlet_triplet.log | head -3
```

能隙的定义要说清楚：

$$
\Delta E_{\rm ST}^{\rm el} = E_{\rm el}({}^1A_1) - E_{\rm el}({}^3B_1),\qquad
\Delta E_{\rm ST}^{0} = \Delta E_{\rm ST}^{\rm el} + \Delta {\rm ZPE}
$$

实验推得的 $T_0 \approx 3147\ {\rm cm^{-1}} \approx 9.0$ kcal/mol（[Jensen & Bunker 1988](https://doi.org/10.1063/1.455184)），对应的是含零点能的 0–0 跃迁，所以应拿 `Sum of electronic and zero-point Energies` 之差去比较。

可以用 [../scripts/g16_log_summary.py](../scripts/g16_log_summary.py) 一次提取多步作业的能量、ZPE 和 ⟨S²⟩（用 `python ../scripts/g16_log_summary.py -h` 查看用法）。

### 5.3 UKS 与 ROKS 怎么选

| 比较项 | UKS（第 1 步） | ROKS（第 4 步） |
|---|---|---|
| ⟨S²⟩ | 略大于 2.000 | 严格 2.000 |
| 能量 | 通常略低（有自旋极化） | 略高 |
| 自旋密度 | 有负自旋密度（自旋极化），更接近实验 ESR | 只在单占据轨道上 |
| 解析频率、TD-DFT 等后续性质 | 支持最完整 | 支持较少 |

实践建议：有机自由基和轻元素三重态一般直接用 UKS；若 ⟨S²⟩ 超过警戒线（§2.3），用 ROKS 单点做对照，两者能量差很大时说明结果对自旋污染敏感，需要谨慎解读。

### 5.4 过渡金属的高自旋/低自旋比较

过渡金属配合物通常要把几种可能的多重度都算一遍：

```text
Fe(II) d6 octahedral:  0 1  (low spin)   0 3  (intermediate)   0 5  (high spin)
Fe(III) d5:            0 2               0 4                   0 6
```

（以上为中性配合物的示意；带电配合物请把第一个数字改为实际电荷。）

注意事项：

1. 每个多重度**各自优化结构**并做频率，比较 $G$ 或 $E+{\rm ZPE}$，不要在同一个结构上比单点；
2. 每个态都做 `Stable=Opt`，确认不是收敛到激发组态；
3. 自旋态能量差强烈依赖泛函中 HF 交换的比例：杂化成分越高越偏向高自旋。结论最好用两三个泛函（如 B3LYP、TPSSh、M06-L 等）交叉检验，并与实验磁性数据对照；
4. 第一过渡系金属的高自旋态有时 SCF 很难收敛，参考 [第4章](04-sp-scf.md) 的 `SCF=XQC`、`Guess=Read` 等方法。

---

## 6. 自旋污染的诊断与处理

### 6.1 在哪里看 ⟨S²⟩

每次 SCF 收敛后（`SCF Done` 之后）都会打印：

```text
 SCF Done:  E(UB3LYP) =  -xxx.xxxxxxxx     a.u. after   15 cycles
 <Sx>= 0.0000 <Sy>= 0.0000 <Sz>= 0.5000 <S**2>= 0.7534 S= 0.5011
 Annihilation of the first spin contaminant:
 S**2 before annihilation     0.7534,   after     0.7500
```

（示意输出，数值仅供参考；格式与 G16 输出一致。）

作业末尾的 archive 段也有汇总：`S2=` 为湮灭前，`S2A=` 为湮灭后。批量检查：

```bash
grep "S\*\*2 before annihilation" *.log | awk '{print $1, $5, $7}'
grep -h "S2=" *.log            # archive lines, may be wrapped across two lines
```

### 6.2 判据

| 现象 | 解读 | 建议 |
|---|---|---|
| 湮灭前 ⟨S²⟩ 偏离理想值 < 5% | 污染可忽略 | 正常使用 |
| 偏离 5–10% | 轻度污染 | 报告 ⟨S²⟩；用 ROKS 或不同泛函做对照 |
| 偏离 > 10%，湮灭后接近理想值 | 主要混入 S+1 态 | 检查是否收敛到错误的态；考虑 ROKS、改用 HF 成分较低的泛函，或该体系本身有多参考特征 |
| 湮灭后仍明显偏离 | 混入更高的自旋态 | 单行列式描述不可靠，考虑多参考方法（§10） |
| 单重态 ⟨S²⟩ 在 0.2–1.0 之间 | BS 解，有双自由基特征 | 这是预期行为，用 AP 投影修正能量（§8） |

> [!NOTE]
> UHF 的自旋污染通常远大于 UKS。用 UHF 作为 MP2/CCSD(T) 参考时问题更严重，可考虑 ROHF 参考（Gaussian 支持 ROMP2 等方法）。

### 6.3 自旋密度

开壳层体系一定要看未成对电子在哪里。G16 的默认布居分析会打印：

```bash
grep -A 30 "Mulliken charges and spin densities:" job.log | head -35
```

自旋密度的正负和大小能告诉你：自由基中心是否符合化学直觉；BS 态里两个磁中心是否一正一负；过渡金属上的未成对电子数是否与所设多重度一致。自旋密度等值面可以用 `cubegen` 生成 spin density cube 后在 GaussView 或 VMD 中查看，见 [第16章](16-postprocess-multiwfn.md)。

---

## 7. 波函数稳定性分析：Stable 与 Stable=Opt

### 7.1 Stable 的选项（G16 官方文档）

| 选项 | 含义 |
|---|---|
| `RExt`（默认） | 检查内部不稳定性和实数外部不稳定性 |
| `Int` | 只检查内部不稳定性 |
| `RRHF`（同 `Singlet`） | 测试/重新优化限制在实数、限制性波函数内 |
| `RUHF`（同 `Triplet`） | 限制在实数、非限制性波函数内 |
| `CRHF` / `CUHF` | 检查实→复不稳定性（L902，代码较旧，非必要不用） |
| `Opt`（同 `RepOpt`） | 发现不稳定就降低约束重新优化，并反复检查直到稳定；默认是 `NoOpt` |
| `1Opt` | 只重做一次 SCF |
| `QCOnly`（`Stable=Opt` 时默认） | 后续 SCF 只用二次收敛算法（L508） |
| `XQC` | 后续 SCF 先试常规 SCF（L502），失败再用二次收敛 |
| `Restart` | 从 chk 重启 |

### 7.2 输出解读

稳定性分析在 SCF 之后进行，输出格式如下（示意输出，数值仅供参考）：

```text
 ***********************************************************************
 Stability analysis using <AA,BB:AA,BB> singles matrix:
 ***********************************************************************
 Eigenvectors of the stability matrix:
 Eigenvector   1:      Triplet-?Sym  Eigenvalue=-0.0123456  <S**2>=2.000
      24 -> 25        0.70xxx
 The wavefunction has an RHF -> UHF instability.
```

| 关键信息 | 含义 |
|---|---|
| `Eigenvalue=` 为负 | 沿该轨道旋转能量会下降，存在不稳定性 |
| `Triplet-...` 本征向量为负 | RHF→UHF（外部）不稳定 |
| `The wavefunction has an internal instability.` | 同一约束下存在更低的解 |
| `The wavefunction is stable under the perturbations considered.` | 在所检查的扰动下稳定 |
| `The wavefunction is already stable.` | `Stable=Opt` 发现无需再优化 |

```bash
grep -i "instability\|stable under\|already stable" job.log
grep "Eigenvalue=" job.log | head
```

### 7.3 标准工作流

```text
Stable=Opt (at the start geometry)
   -> Opt Freq with Guess=Read Geom=Check, written as U<functional> explicitly
   -> Stable at the optimized geometry (repeat the loop if unstable)
   -> single points / properties, again with Guess=Read
```

> [!IMPORTANT]
> 1. 从 `RB3LYP Stable=Opt` 得到 UKS 解后，后续步骤必须写 `UB3LYP ... Guess=Read`，多重度为 1 时不写 `U` 会退回 RKS。
> 2. 几何优化过程中波函数可能“滑回”闭壳层或别的解，**优化结束后要在最终结构上再做一次 `Stable`**。
> 3. 官方文档指出：解析频率只有在波函数不存在内部不稳定性时才有效；限制性波函数需没有 singlet 不稳定性，非限制性波函数需没有任何（singlet 或 triplet）不稳定性。Møller-Plesset 计算同样要求参考波函数稳定。
> 4. `Stable` 要分析的是**收敛后**的波函数。不要在一个连 SCF 都没收敛的作业里加 `Stable` 试图“救火”。

---

## 8. 实战三：单重态双自由基——对苯炔（p-benzyne）

### 8.1 背景

对苯炔（1,4-didehydrobenzene）的两个自由基中心相距较远，单重态有明显双自由基特征；实验测得单重态比三重态低约 3.8 kcal/mol（[Wenthold, Squires & Lineberger 1998](https://doi.org/10.1021/ja9803355)）。闭壳层 RB3LYP 对它不稳定，是演示 `Stable=Opt` 和 `Guess=Mix` 的理想体系。

### 8.2 完整输入

文件 [../examples/13-diradical-pbenzyne/benzyne_stable_opt_guess_mix.gjf](../examples/13-diradical-pbenzyne/benzyne_stable_opt_guess_mix.gjf)（坐标为按 D2h 理想六边形搭建的起始结构）：

```text
%nprocshared=8
%mem=8GB
%chk=pbz_rks.chk
#p RB3LYP/6-31G(d) Stable=Opt NoSymm

p-benzyne step 1: start from RKS; Stable=Opt may switch to UKS if unstable

0 1
C                   1.390000    0.000000    0.000000
C                   0.695000    1.203775    0.000000
C                  -0.695000    1.203775    0.000000
C                  -1.390000    0.000000    0.000000
C                  -0.695000   -1.203775    0.000000
C                   0.695000   -1.203775    0.000000
H                   1.235000    2.139083    0.000000
H                  -1.235000    2.139083    0.000000
H                  -1.235000   -2.139083    0.000000
H                   1.235000   -2.139083    0.000000

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=pbz_rks.chk
%chk=pbz_mix.chk
#p UB3LYP/6-31G(d) Guess=Mix Stable=Opt NoSymm Geom=Check

p-benzyne step 2: same geometry, broken-symmetry guess from Guess=Mix

0 1

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=pbz_rks.chk
%chk=pbz_bs_opt.chk
#p UB3LYP/6-31G(d) Opt Freq NoSymm Geom=Check Guess=Read

p-benzyne step 3: optimize the broken-symmetry singlet (stable UKS guess)

0 1

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=pbz_bs_opt.chk
%chk=pbz_bs_stab.chk
#p UB3LYP/6-31G(d) Stable NoSymm Geom=AllCheck Guess=Read Pop=NaturalOrbitals

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=pbz_bs_opt.chk
%chk=pbz_uno.chk
#p Guess=(Save,Only,NaturalOrbitals) Geom=AllCheck ChkBasis

--Link1--
%nprocshared=8
%mem=8GB
%oldchk=pbz_rks.chk
%chk=pbz_triplet.chk
#p UB3LYP/6-31G(d) Opt Freq NoSymm Geom=Check

p-benzyne step 6: triplet state for the singlet-triplet gap

0 3

```

| 步骤 | 内容 | 期望看到 |
|---|---|---|
| 1 | `RB3LYP Stable=Opt` | 先报 `RHF -> UHF instability`，随后以 `E(UB3LYP)` 收敛到更低的 BS 解 |
| 2 | `UB3LYP Guess=Mix Stable=Opt` | 另一条路线，应得到与第 1 步相同（或非常接近）的能量和 ⟨S²⟩ |
| 3 | BS 单重态 `Opt Freq` | 从第 1 步的稳定 UKS 轨道出发，在 BS 势能面上优化 |
| 4 | 最终结构上的 `Stable` + `Pop=NaturalOrbitals` | 确认仍稳定；打印 UNO 占据数 |
| 5 | `Guess=(Save,Only,NaturalOrbitals)` | 把自然轨道存入 `pbz_uno.chk`（写法取自官方 Population 关键词页面），供可视化或交给 MOKIT |
| 6 | 三重态 `Opt Freq` | 用于单重态-三重态能隙和 AP 投影 |

> [!TIP]
> `Guess=Mix` 和 `Stable=Opt` 两条路线各有用处：`Guess=Mix` 便宜、直接；从 RKS 出发的 `Stable=Opt` 还能同时告诉你闭壳层解是否真的不稳定。对新体系建议两条都试，取能量更低、且稳定的那个解。

### 8.3 判读

1. **第 1 步与第 2 步的能量**应一致到约 $10^{-6}$ Hartree；若不一致，说明存在多个 UKS 解，取更低且稳定者；
2. **BS 单重态的 ⟨S²⟩** 应明显大于 0 而小于 1，数值越接近 1，双自由基特征越强；
3. **单重态-三重态能隙**（$E$ 取 `Sum of electronic and zero-point Energies`）：

$$
\Delta E_{\rm ST}^{\rm BS} = E_{\rm BS} - E_{\rm T},\qquad
\Delta E_{\rm ST}^{\rm AP} = E_{\rm S}^{\rm AP} - E_{\rm T} = (1+\alpha)\,(E_{\rm BS}-E_{\rm T})
$$

（AP 式中的 $E_{\rm BS}$、$E_{\rm T}$ 严格说应在同一几何下计算；常见做法是各自在优化结构上计算，近似地把几何松弛也包括进去，报告时说明即可。）

4. **UNO 占据数与双自由基指数**：UKS 总密度的自然轨道中，最高占据自然轨道（HONO）和最低未占自然轨道（LUNO）的占据数成对出现（$2-n$ 与 $n$）。LUNO 占据数 $n_{\rm LUNO}$ 越接近 1，双自由基特征越强。常用的 Yamaguchi 双自由基指数

$$
y = 1 - \frac{2T}{1+T^2},\qquad T = \frac{n_{\rm HONO}-n_{\rm LUNO}}{2}
$$

$y=0$ 为闭壳层，$y=1$ 为纯双自由基。MOKIT 的 Python 接口也提供从 fch 计算 Yamaguchi/Head-Gordon 未成对电子数的函数（`calc_unpaired_from_fch`），见 §10。

> [!NOTE]
> 邻苯炔（o-benzyne）的两个自由基中心距离近，形成较强的 π 键，闭壳层描述基本可行（实验 S–T 能隙约 37.5 kcal/mol，同一参考文献）；间位、对位依次增强双自由基特征。用这三个异构体做对比是很好的练习。

---

## 9. 实战四：Broken-symmetry 与磁耦合常数 J（双核 Cu(II) 模型）

### 9.1 体系与碎片设计

模型 $[{\rm Cu_2(\mu\text{-}OH)_2(NH_3)_4}]^{2+}$：两个 Cu(II)（d⁹，各一个未成对电子，$S_A=S_B=1/2$）由两个 OH⁻ 桥连，Cu–O 1.93 Å，Cu–O–Cu 约 102°，Cu···Cu 3.00 Å（手工搭建的平面 Cu2O2 核，不是晶体结构，也没有优化）。经验上，这类羟桥双核铜的 Cu–O–Cu 角大于约 97.5° 时倾向反铁磁（[Crawford et al. 1976](https://doi.org/10.1021/ic50163a019) 的磁-构关联），所以预期 J < 0。

| 碎片 | 包含原子 | 碎片电荷 | 碎片多重度 | 说明 |
|---|---|---|---|---|
| 1 | Cu1 + 两个 NH3 | +2 | 2 | 未成对电子取 α |
| 2 | Cu2 + 两个 NH3 | +2 | **−2** | 负号：未成对电子取 β |
| 3 | 两个桥连 OH⁻ | −2 | 1 | 闭壳层 |
| 整体 | — | +2 | 1 | BS 态 $M_S=0$ |

电荷多重度行依次写“整体、碎片 1、碎片 2、碎片 3”：`2 1 2 2 2 -2 -2 1`。碎片编号必须从 1 开始连续编号；负多重度只在 `Guess=Fragment` 作业中合法，在其他作业类型中会报错（官方分子说明文档）。

### 9.2 完整输入

文件 [../examples/13-bs-cu2-dimer/cu2_dimer_bs_guess_fragment.gjf](../examples/13-bs-cu2-dimer/cu2_dimer_bs_guess_fragment.gjf)（第 1 步与第 2 步坐标相同，此处第 1 步只列开头几行）：

```text
%nprocshared=16
%mem=16GB
%chk=cu2_hs.chk
#p UB3LYP/def2SVP Stable=Opt NoSymm

[Cu2(mu-OH)2(NH3)4]2+ step 1: high-spin triplet (Ms=1), Cu-O-Cu 102 deg

2 3
Cu                  1.500000    0.000000    0.000000
N                   2.903210    1.453066    0.000000
...                 (22 atoms, same coordinates as in step 2 without the labels)

--Link1--
%nprocshared=16
%mem=16GB
%chk=cu2_bs.chk
#p UB3LYP/def2SVP Guess=(Fragment=3,Only) NoSymm

[Cu2(mu-OH)2(NH3)4]2+ step 2: fragment guess, Cu1 alpha / Cu2 beta doublet

2 1 2 2 2 -2 -2 1
Cu(Fragment=1)      1.500000    0.000000    0.000000
N(Fragment=1)       2.903210    1.453066    0.000000
H(Fragment=1)       3.138570    1.696789   -0.956781
H(Fragment=1)       2.542527    2.272380    0.478391
H(Fragment=1)       3.734612    1.121197    0.478391
N(Fragment=1)       2.903210   -1.453066    0.000000
H(Fragment=1)       3.138570   -1.696789   -0.956781
H(Fragment=1)       3.734612   -1.121197    0.478391
H(Fragment=1)       2.542527   -2.272380    0.478391
Cu(Fragment=2)     -1.500000    0.000000    0.000000
N(Fragment=2)      -2.903210   -1.453066    0.000000
H(Fragment=2)      -3.138570   -1.696789   -0.956781
H(Fragment=2)      -2.542527   -2.272380    0.478391
H(Fragment=2)      -3.734612   -1.121197    0.478391
N(Fragment=2)      -2.903210    1.453066    0.000000
H(Fragment=2)      -3.138570    1.696789   -0.956781
H(Fragment=2)      -3.734612    1.121197    0.478391
H(Fragment=2)      -2.542527    2.272380    0.478391
O(Fragment=3)       0.000000    1.214455    0.000000
H(Fragment=3)       0.000000    2.184455    0.000000
O(Fragment=3)       0.000000   -1.214455    0.000000
H(Fragment=3)       0.000000   -2.184455    0.000000

--Link1--
%nprocshared=16
%mem=16GB
%chk=cu2_bs.chk
#p UB3LYP/def2SVP Guess=Read Geom=AllCheck Stable=Opt NoSymm

```

| 步骤 | 关键词 | 作用 |
|---|---|---|
| 1 | `2 3` + `Stable=Opt` | 高自旋三重态（$M_S=1$），先确认 HS 解稳定，它是 J 公式里的参考 |
| 2 | `Guess=(Fragment=3,Only)` | 只生成碎片初猜：各碎片先单独生成初猜，再拼成 Cu1 α、Cu2 β 的整体初猜并写入 chk，不做整体 SCF |
| 3 | `Guess=Read Geom=AllCheck Stable=Opt` | 从碎片初猜出发做 BS 的 SCF，并用 `Stable=Opt` 保证得到稳定解 |
| 全部 | `NoSymm` | BS 解不保持 Cu2O2 核的空间对称性 |

这种“高自旋 → 碎片初猜（`Only`）→ `Guess=Read` + `Stable=Opt`”的流程与官方反铁磁耦合说明（[gaussian.com/afc](https://gaussian.com/afc/)）的思路一致。在 GaussView 6 中也可以用 Tools → Atom Groups，把 Atom Group Class 设为 Gaussian Fragment，逐个定义碎片的电荷、多重度和 α/β 自旋，再在计算设置的 Guess 页面勾选使用碎片初猜。

### 9.3 检查 BS 解是否正确

```bash
cd Gaussian/examples/13-bs-cu2-dimer
grep -c "Normal termination" cu2_dimer_bs_guess_fragment.log        # expect 3
grep "S\*\*2 before annihilation" cu2_dimer_bs_guess_fragment.log
grep -i "instability\|stable under\|already stable" cu2_dimer_bs_guess_fragment.log
grep -A 25 "Mulliken charges and spin densities:" cu2_dimer_bs_guess_fragment.log | grep " Cu "
```

（示意输出，数值仅供参考）

```text
     1  Cu   0.9xxxxx   0.5xxxxx        <- HS step: both Cu positive
    10  Cu   0.9xxxxx   0.5xxxxx
     1  Cu   0.9xxxxx   0.5xxxxx        <- BS step: opposite signs
    10  Cu   0.9xxxxx  -0.5xxxxx
```

| 检查项 | 正确的 BS 解 | 常见错误 |
|---|---|---|
| ⟨S²⟩（湮灭前） | HS ≈ 2.0，BS ≈ 1.0 | BS ≈ 0：塌缩成闭壳层；BS ≈ 2：其实是 $M_S=1$ |
| 两个 Cu 的自旋密度 | 大小相近、符号相反 | 同号：读错了 chk 或多重度写成 3 |
| 自旋密度大小 | 明显小于 1（部分离域到 O、N 上） | 某个 Cu 接近 0：未成对电子跑到配体上，检查结构和电荷 |
| 稳定性 | `stable under the perturbations considered` | 仍不稳定：继续 `Stable=Opt` 或换碎片划分 |

### 9.4 计算 J

示例目录提供了 [../examples/13-bs-cu2-dimer/calc_j_yamaguchi.py](../examples/13-bs-cu2-dimer/calc_j_yamaguchi.py)（仅用 Python 标准库），它按 Link1 步骤切分 log，取指定步骤最后一次 `SCF Done` 能量和湮灭前的 ⟨S²⟩：

```bash
python calc_j_yamaguchi.py --hs cu2_dimer_bs_guess_fragment.log:1 \
                           --bs cu2_dimer_bs_guess_fragment.log:3 --smax 1
```

输出格式（示意输出，数值仅供参考）：

```text
HS : E(UB3LYP) = -xxxx.xxxxxxxx Eh   <S**2> = 2.00xx (ideal 2.0000)
BS : E(UB3LYP) = -xxxx.xxxxxxxx Eh   <S**2> = 1.00xx (weak-overlap ideal 1.0000)
E_BS - E_HS = -0.000xxxxx Eh = -xxx.x cm-1 = -x.xxx kcal/mol

formula                            J (H=-2J S.S)    J (H=-J S.S)   [cm-1]
Yamaguchi                                 -xxx.x          -xxx.x
Noodleman (weak overlap)                  -xxx.x          -xxx.x
Ruiz/Bencini (strong overlap)             -xxx.x          -xxx.x
```

手算也很简单：1 Hartree = 219474.63 cm⁻¹。例如（虚构数字，仅演示换算）$E_{\rm BS}-E_{\rm HS} = -0.00100$ Eh，$\langle S^2\rangle_{\rm HS}=2.004$，$\langle S^2\rangle_{\rm BS}=1.002$，则 Yamaguchi $J = -0.00100/1.002 \times 219474.6 \approx -219$ cm⁻¹（$\hat H=-2J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$），单重态-三重态能隙 $E_{\rm T}-E_{\rm S} = -2J \approx 438$ cm⁻¹。

### 9.5 让 J 更可信

- **结构**：用晶体结构或在 HS 态上优化的结构；J 对桥连角、桥上 H 的取向非常敏感；
- **泛函**：J 强烈依赖 HF 交换比例。B3LYP 是最常用的选择之一，但应说明并最好用另一个泛函对照；不同公式（Noodleman/Ruiz/Yamaguchi）与泛函的搭配在文献中有争论，选择一种并保持一致；
- **基组**：Cu 上至少用三 ζ 级别（如 def2-TZVP）做最终单点，示例为节省时间用了 def2-SVP；
- **多核体系**：n 个磁中心需要多个 BS 构型（不同的 α/β 分布），用最小二乘拟合多个 J；Fe2S2 这类体系的写法可参考官方 `Guess` 页面中的 Fe2S2 碎片初猜示例；
- **多个 BS 解**：重金属多核体系常存在多个 UKS 解，务必检查自旋密度分布，并取能量最低、稳定的解。

---

## 10. 交给多参考方法：UNO 与 MOKIT

当出现以下情况时，单行列式 DFT 的结论需要多参考方法验证：湮灭后 ⟨S²⟩ 仍明显偏离理想值；LUNO 占据数接近 1；不同泛函给出相反的自旋态排序；需要定量的单重态-三重态能隙或 J。

[MOKIT](https://github.com/1234zou/MOKIT)（邹靖翔开发，Apache-2.0）可以在 Gaussian、PySCF、ORCA、OpenMolcas 等程序之间传递轨道，其 `automr` 程序以 Gaussian 输入文件为输入，自动完成 UHF → UNO → GVB → CASSCF → NEVPT2/CASPT2 等流程。思路如下（具体关键词以 [MOKIT 文档](https://jeanwsr.gitlab.io/mokit-doc-mdbook/) 为准）：

1. 在 Gaussian 中得到**稳定的** UHF 解（用 `Stable=Opt` 检查），保存 chk 后 `formchk` 成 fch；MOKIT 文档建议 `readuhf` 读入 UHF 而不是 UDFT 轨道；
2. 写一个 `automr` 输入，route 写目标方法（如 `CASSCF/cc-pVDZ` 或 `NEVPT2/cc-pVDZ`），在 `mokit{}` 中用 `ist=1,readuhf='xxx.fch'` 读入轨道；
3. 若已有自然轨道（如 UNO、UDFT NO），文档还提供 `ist=5` 配合 `readno` 的路线；
4. 结果给出活性空间、自然轨道占据数和多参考能量，可与 BS-DFT 的结论对照。

```text
%mem=16GB
%nprocshared=8
#p CASSCF/cc-pVDZ
mokit{ist=1,readuhf='pbz_uhf.fch'}

```

（以上为 `automr` 输入的格式示意，`pbz_uhf.fch` 需要是你自己用 UHF + `Stable=Opt` 得到并 formchk 的文件；运行方式为 `automr pbz_cas.gjf > pbz_cas.out 2>&1`。）

> [!CAUTION]
> 多参考计算的活性空间选择会显著影响结果，`automr` 自动选择的活性空间也要人工检查（看自然轨道形状和占据数）。初学者先用小体系熟悉流程，再用于研究体系。

---

## 参数选择与收敛建议

| 项目 | 推荐做法 | 原因 |
|---|---|---|
| 开壳层单重态 | 显式写 `U` + `Guess=Mix` 或 `Stable=Opt` | 默认 R 永远得不到 BS 解 |
| 新体系首次计算 | 先 `Stable=Opt`，再优化 | 避免在错误的态上优化 |
| 优化之后 | 在最终结构上再做一次 `Stable` | 优化中波函数可能换解 |
| BS 计算 | `NoSymm`（或 `SCF=NoSymm`） | BS 解破坏空间对称性 |
| 读入轨道 | `%oldchk` + `Geom=Check/AllCheck` + `Guess=Read` | 每步写新 chk，可追溯 |
| 能量比较 | 各自旋态各自优化，比较 $E+{\rm ZPE}$ 或 $G$ | 不同自旋态几何差别可能很大 |
| 自旋污染判据 | 湮灭前 ⟨S²⟩ 偏离 > 10% 需警惕 | 湮灭后的值只是诊断 |
| J 的计算 | HS 与 BS 同一结构、同一方法基组；写明哈密顿量约定 | 约定不同数值差 2 倍或符号相反 |
| 泛函 | 自旋态能量差、J 用两种以上泛函交叉验证 | HF 交换比例影响大 |
| SCF 难收敛 | `SCF=XQC`、`Stable=(Opt,XQC)`、先小基组再 `Guess=Read` | 详见第4章 |

---

## 常见问题与报错

**1. `The combination of multiplicity 1 and 17 electrons is impossible.`（L301）**

- 原因：电荷与多重度不匹配，奇数电子不能是单重态。
- 处理：重新数电子 $N_e=\sum Z - Q$；奇数电子用 2、4、6……，偶数电子用 1、3、5……。碎片作业中每个碎片也要分别满足这一条件。

**2. 想算开壳层单重态，结果 ⟨S²⟩ = 0.0000，能量与 RKS 相同**

- 原因：route 写的是 `B3LYP`（单重态默认 R），或写了 `UB3LYP` 但初猜对称，SCF 收敛到闭壳层解。
- 处理：写 `UB3LYP Guess=Mix`，必要时加 `NoSymm`，并用 `Stable=Opt` 确认。

**3. `Stable=Opt` 之后的 `Opt` 又回到了闭壳层**

- 原因：后续步骤 route 没写 `U`，或没写 `Guess=Read`，或 `%oldchk` 指向了错误的文件。
- 处理：后续步骤统一写 `UB3LYP ... Guess=Read`，优化结束后再做一次 `Stable` 检查。

**4. BS 计算收敛到了高自旋态或闭壳层（⟨S²⟩ ≈ 2 或 ≈ 0）**

- 原因：碎片多重度写错（忘了负号）；碎片划分让未成对电子落在桥连配体上；或 SCF 过程中自旋翻转。
- 处理：检查电荷多重度行；看初猜步骤（`Guess=Only`）打印的初猜 ⟨S²⟩ 和自旋分布；必要时调整碎片（例如把桥连配体单独作为闭壳层碎片），并用 `Stable=Opt`。

**5. 负多重度导致作业报错**

- 原因：负多重度只在 `Guess=Fragment` 作业中合法，官方文档说明在其他作业类型中会产生错误。
- 处理：负多重度只用于生成初猜的那一步；后续步骤用 `Geom=AllCheck` 从 chk 读取，或只写整体的电荷多重度。

**6. `Convergence failure -- run terminated.`（L502，开壳层过渡金属常见）**

- 处理顺序：先检查电荷、多重度和结构；再用 `SCF=XQC`；用小基组或高自旋态的轨道 `Guess=Read`；`Stable=(Opt,XQC)`。不要只加大 `MaxCycle`。详见 [第4章](04-sp-scf.md)。

**7. 频率计算结果异常（出现不合理的虚频）但优化正常收敛**

- 原因之一：波函数不稳定，此时解析二阶导数没有物理意义。
- 处理：在该结构上 `Stable=Opt`，用稳定的波函数重新优化和计算频率。

**8. 湮灭后 ⟨S²⟩ 接近理想值，就认为没有自旋污染**

- 这是误用：湮灭只是诊断，打印的能量并没有投影。应以**湮灭前**的值判断，并在论文中报告。

**9. 用 ROKS 算频率或 TD-DFT 遇到问题**

- ROKS 在 Gaussian 中的后续性质支持不如 UKS 完整。常见做法：结构和频率用 UKS，ROKS 只做单点能对照。

---

## 练习

1. 用示例中的 O2 输入，把泛函换成 PBE0（`PBE1PBE`）和 M06-2X，比较三种泛函下 $2(E_{\rm BS}-E_{\rm T})$ 与实验 0.98 eV 的差别。
2. 对 CH2，把第 2 步改为 `UB3LYP Guess=Mix` 并加 `Stable=Opt`，验证 ¹A1 是否为闭壳层（⟨S²⟩ 是否为 0）。再用 `Sum of electronic and zero-point Energies` 计算 $T_0$ 并与 9.0 kcal/mol 比较。
3. 对邻苯炔、间苯炔、对苯炔分别做 §8 的流程，列表比较 RHF→UHF 不稳定性是否出现、BS 态 ⟨S²⟩、LUNO 占据数和 S–T 能隙。
4. 把 Cu2 模型的 Cu–O–Cu 角从 95° 改到 105°（每 2.5° 一个点，保持 Cu–O 键长），画出 J 随角度的变化，观察符号是否在 97–98° 附近改变。
5. 对 Cu2 模型分别用 Noodleman、Ruiz、Yamaguchi 公式计算 J，比较差异，并说明你选择的约定。

---

## 参考资料

**官方文档**

- [Stable keyword](https://gaussian.com/stable/)：稳定性分析选项与适用范围
- [Guess keyword](https://gaussian.com/guess/)：`Mix`、`Alter`、`Fragment`、`Read`、`NaturalOrbitals`，含 Fe2S2 碎片初猜示例
- [Antiferromagnetic Coupling](https://gaussian.com/afc/)：Mn2O2 体系的 HS → Stable=Opt → 碎片初猜 → BS 完整流程
- [Molecule Specification](https://gaussian.com/molspec/)：碎片定义、负多重度的含义
- [Population keyword](https://gaussian.com/population/)：`NaturalOrbitals` 与把自然轨道存入 chk 的写法
- [SCF keyword](https://gaussian.com/scf/)：`XQC`、`QC`、`NoSymm` 等

**相关 GitHub 仓库**

- [1234zou/MOKIT](https://github.com/1234zou/MOKIT)：在 Gaussian（fch）、PySCF、ORCA、OpenMolcas 等程序间传递轨道，`automr` 自动完成 UHF/UNO → GVB → CASSCF/NEVPT2，适合把 BS/UNO 结果交给多参考方法（中文作者，有中文资料）
- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian 输出，可批量读取能量和自旋相关信息
- [chauncey-garrett/gaussian-tools](https://github.com/chauncey-garrett/gaussian-tools)：小脚本集合，其中有从 log 提取湮灭前后自旋污染的工具，可作为“需要从 log 抓什么”的清单
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：中文项目，其结果自查清单包含自旋污染检查项

**关键文献**

- R. Seeger, J. A. Pople, *J. Chem. Phys.* 1977, 66, 3045. [DOI: 10.1063/1.434318](https://doi.org/10.1063/1.434318)（SCF 稳定性）
- R. Bauernschmitt, R. Ahlrichs, *J. Chem. Phys.* 1996, 104, 9047. [DOI: 10.1063/1.471637](https://doi.org/10.1063/1.471637)（DFT 稳定性分析）
- L. Noodleman, *J. Chem. Phys.* 1981, 74, 5737. [DOI: 10.1063/1.440939](https://doi.org/10.1063/1.440939)（broken-symmetry 方法）
- T. Ziegler, A. Rauk, E. J. Baerends, *Theor. Chim. Acta* 1977, 43, 261. [DOI: 10.1007/BF00551551](https://doi.org/10.1007/BF00551551)（求和法）
- K. Yamaguchi, F. Jensen, A. Dorigo, K. N. Houk, *Chem. Phys. Lett.* 1988, 149, 537. [DOI: 10.1016/0009-2614(88)80378-6](https://doi.org/10.1016/0009-2614(88)80378-6)（近似自旋投影）
- T. Soda et al., *Chem. Phys. Lett.* 2000, 319, 223. [DOI: 10.1016/S0009-2614(00)00166-4](https://doi.org/10.1016/S0009-2614(00)00166-4)（Yamaguchi 公式求 J 的比较）
- P. Jensen, P. R. Bunker, *J. Chem. Phys.* 1988, 89, 1327. [DOI: 10.1063/1.455184](https://doi.org/10.1063/1.455184)（CH2 单重态-三重态能隙）
- P. G. Wenthold, R. R. Squires, W. C. Lineberger, *J. Am. Chem. Soc.* 1998, 120, 5279. [DOI: 10.1021/ja9803355](https://doi.org/10.1021/ja9803355)（苯炔异构体 S–T 能隙）
- V. H. Crawford et al., *Inorg. Chem.* 1976, 15, 2107. [DOI: 10.1021/ic50163a019](https://doi.org/10.1021/ic50163a019)（羟桥双核铜的 J–角度关联）

延伸阅读：DFT 理论基础与学习路线见 ../../Common/tutorials/07-learning-path.md；批量提取 ⟨S²⟩、能量用 [../scripts/g16_log_summary.py](../scripts/g16_log_summary.py)，报错诊断用 ../scripts/g16_error_doctor.py（均可用 `-h` 查看用法）。

---

[← 上一章](12-population-charges.md) | [返回目录](../README.md) | [下一章 →](14-weak-bsse.md)
