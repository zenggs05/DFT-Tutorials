# 第3章 GPW/GAPW 方法、基组与赝势（GPW/GAPW Methods, Gaussian Basis Sets & GTH Pseudopotentials in CP2K）

本章讲清 CP2K Quickstep 的两种核心方法：GPW（Gaussian and Plane Waves，高斯基表示轨道、平面波表示密度、多重网格）与 GAPW（Gaussian and Augmented Plane Waves，全电子与近核性质），以及与之配套的 GTH 赝势（Goedecker–Teter–Hutter pseudopotentials）和高斯基组（MOLOPT、MOLOPT-SR、UZH、全电子 ae 基组）。重点是实用规则：`-qN` 价电子数必须匹配、`BASIS_MOLOPT_UZH` 与 `POTENTIAL_UZH` 成对使用、2026.2 新增的带版本号数据快照、辅助基组 `AUX_FIT`/`RI_AUX` 的写法，以及“DZVP 做结构与 MD、TZV2P 做能量精修”的基组选择流程（CP2K basis set and pseudopotential selection）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 用自己的话解释 GPW 中“轨道在高斯基、密度在网格”的分工，以及 `CUTOFF`、`REL_CUTOFF`、`NGRIDS` 分别控制什么；
- 判断什么时候需要 GAPW 全电子计算，并写出 `METHOD GAPW` + `POTENTIAL ALL` + 全电子基组的输入；
- 读懂基组名（`TZV2P-MOLOPT-PBE-GTH-q6`）和赝势名（`GTH-PBE-q6`），读懂数据文件中基组/赝势条目的格式；
- 为任意元素选出**互相匹配**的基组和赝势，并知道去哪个文件里找；
- 在输入中同时使用多个基组文件、为 ADMM/RI 指定辅助基组，用 `@INCLUDE` 复用 `&KIND` 模板；
- 理解 2026.2 数据目录中 `BASIS_MOLOPT_UZH_2026.2` 这类快照文件的意义，写出可复现的输入。

| 场景 | 推荐方法与基组 | 本章小节 |
|---|---|---|
| 常规凝聚相 DFT、AIMD、表面吸附 | GPW + UZH `DZVP-MOLOPT-PBE-GTH-qN`（或经典 `DZVP-MOLOPT-SR-GTH`）+ 匹配的 GTH 赝势 | §4、§8 |
| 分子性质、反应能、需要更高精度的能量 | GPW + `TZV2P` 级别基组精修 | §7、§11 |
| 与已发表的 CP2K 结果对比 | 沿用原文的经典 `BASIS_MOLOPT` + `GTH_POTENTIALS` | §4 |
| 全电子、内层能级（XAS）、NMR/超精细等近核性质 | GAPW + `POTENTIAL ALL` + `SVP/TZVPP-MOLOPT-PBE-ae` 或 Pople/Dunning 基组 | §2.3、§9 |
| 杂化泛函（ADMM） | 轨道基组 + `BASIS_SET AUX_FIT admm-dzp-qN` | §10、[第9章](09-hybrid-admm.md) |

> [!NOTE]
> 本章示例在 [`examples/03-basis-pseudo-h2o/`](../examples/03-basis-pseudo-h2o/README.md)。截断参数的系统收敛见 [第4章](04-cutoff-convergence.md)；本章只在需要时给出足够保守的取值。

---

## 2. 原理速览

### 2.1 GPW：轨道用高斯基，密度用平面波

Kohn–Sham 轨道展开在原子中心的收缩高斯函数上，密度矩阵 $P_{\mu\nu}$ 是 SCF 求解的对象：

$$
\psi_i(\mathbf r)=\sum_\mu c_{\mu i}\,\varphi_\mu(\mathbf r),\qquad
n(\mathbf r)=\sum_{\mu\nu}P_{\mu\nu}\,\varphi_\mu(\mathbf r)\,\varphi_\nu(\mathbf r)
$$

GPW 的总能量写成

$$
E[P]=E^{\mathrm T}[P]+E^{\mathrm V}[P]+E^{\mathrm H}[n]+E^{\mathrm{XC}}[n]+E^{\mathrm{II}}
$$

各项的计算方式不同：

| 能量项 | 在哪里计算 | 说明 |
|---|---|---|
| 动能 $E^{\mathrm T}$、赝势短程与非局域部分 $E^{\mathrm V}$ | 高斯基中**解析积分** | 高斯函数的积分有闭式解，矩阵稀疏 |
| Hartree 能 $E^{\mathrm H}$ | **实空间规则网格 + FFT**，在倒空间解 Poisson 方程 | $E^{\mathrm H}=2\pi\Omega\sum_{\mathbf G\neq0}\lvert\tilde n_{\mathrm{tot}}(\mathbf G)\rvert^2/G^2$，$\tilde n_{\mathrm{tot}}$ 包含电子密度与高斯型“核电荷”密度 |
| 交换关联能 $E^{\mathrm{XC}}$ | 在实空间网格上逐点求值 | 网格越细越准 |
| 离子–离子项 $E^{\mathrm{II}}$ | 解析 | 与核电荷的高斯展宽配合处理 |

把密度从高斯基“投影”到网格叫 **collocation**，把网格上的势“积分”回高斯基矩阵元叫 **integration**。这就是输出中会出现 `Core density on regular grids`、`Self energy of the core charge distribution`、`Overlap energy of the core charge distribution` 等条目的原因：赝势的长程局域部分被表示为高斯形状的核电荷，与电子密度一起在网格上求 Hartree 势。

**两个独立的收敛维度**：

1. **基组质量**（SZV → DZVP → TZV2P …）决定轨道能被表示得多好——这一点 `CUTOFF` 再高也无法弥补；
2. **网格精度**（`CUTOFF`、`REL_CUTOFF`）决定密度与势在网格上的表示误差。

网格截断与网格间距的关系：平面波截断 $E_{\mathrm{cut}}$（Ry）对应最大波矢 $G_{\max}=\sqrt{E_{\mathrm{cut}}}$（bohr⁻¹），网格间距约 $h\approx\pi/G_{\max}$。`CUTOFF 400` Ry 对应 $h\approx0.157$ bohr ≈ 0.083 Å。

### 2.2 多重网格（multigrid）与 REL_CUTOFF

高斯函数乘积 $\varphi_\mu\varphi_\nu$ 有的很尖锐（大指数），有的很平缓（小指数）。把平缓的函数也放到最细网格上纯属浪费，于是 Quickstep 准备了 `NGRIDS`（默认 4）级由细到粗的网格：

$$
E_{\mathrm{cut}}^{(i)}=\frac{E_{\mathrm{cut}}^{(1)}}{\alpha^{\,i-1}},\qquad i=1,\dots,\texttt{NGRIDS},\quad \alpha=\texttt{PROGRESSION\_FACTOR}=3
$$

其中 $E_{\mathrm{cut}}^{(1)}$ 就是 `CUTOFF`。每个高斯乘积被映射到“能让它覆盖足够多网格点的**最粗**那一级”；判据由 `REL_CUTOFF` 决定：它是指数为 1 的参考高斯所需的截断，指数为 $\eta$ 的高斯大致需要 $\eta\times$`REL_CUTOFF` 的截断。

以 `CUTOFF 400`、`REL_CUTOFF 50` 为例：

| 网格级 | 截断 (Ry) | 输出中 `cutoff [a.u.]`（Ha） | 大致承载的高斯指数范围 |
|---|---|---|---|
| 1（最细） | 400 | 200.00 | $\eta > 2.67$ |
| 2 | 133.3 | 66.67 | $0.89 < \eta \le 2.67$ |
| 3 | 44.4 | 22.22 | $0.30 < \eta \le 0.89$ |
| 4（最粗） | 14.8 | 7.41 | $\eta \le 0.30$ |

输出中的 `count for grid N` 行给出每级网格上实际分到的高斯乘积数。两个典型错误：

- `CUTOFF` 太低：所有网格都太粗，能量和力都不准；
- `REL_CUTOFF` 太低：即使 `CUTOFF` 很高，大部分高斯也被扔到粗网格上，等于白提高了 `CUTOFF`（默认值 20 Ry 对多数体系偏低，常用 50–60 Ry）。

> [!TIP]
> `BASIS_MOLOPT` 文件头部的说明指出：经典 MOLOPT 基组含较弥散的函数，使用 `NGRIDS 5` 可以更高效地处理这些平缓的高斯（文件注释中称约有 2 倍提速）。改变 `NGRIDS` 后需要重新确认收敛。完整的两步收敛流程见 [第4章](04-cutoff-convergence.md)。

### 2.3 GAPW：把近核的“硬”密度交给原子网格

赝势把内层电子去掉，使价电子密度在核附近变得平滑。全电子计算没有这个便利：1s 电子的密度在核附近极其尖锐，若全放到均匀网格上，需要的截断高得不现实。GAPW 的做法是把密度分解为

$$
n(\mathbf r)=\tilde n(\mathbf r)+\sum_A\Big[n_A(\mathbf r)-\tilde n_A(\mathbf r)\Big]
$$

- $\tilde n$：全空间的**软**密度，仍放在 GPW 的规则网格上；
- $n_A$、$\tilde n_A$：以原子 $A$ 为中心的**硬**和**软**单中心密度，只在原子周围的球内用原子网格（Lebedev 角向格点 × 径向格点）计算；
- Hartree 项中再引入高斯补偿电荷，保证多极矩正确。

| 参数 | 位置 | 默认值 | 作用 |
|---|---|---|---|
| `METHOD GAPW` | `&DFT/&QS` | `GPW` | 启用 GAPW；`GAPW_XC` 只对交换关联项使用 GAPW 分解 |
| `EPSFIT` | `&DFT/&QS` | 1.0E-4 | 指数“软/硬”划分阈值；调小会把更硬的函数放进软密度，需要更高 `CUTOFF` |
| `EPSRHO0` | `&DFT/&QS` | 1.0E-6 | 硬补偿电荷的作用范围 |
| `EPSSVD` | `&DFT/&QS` | 1.0E-8 | 投影子矩阵 SVD 的容差 |
| `LEBEDEV_GRID` | `&KIND` | 50 | 原子网格的角向格点数 |
| `RADIAL_GRID` | `&KIND` | 50 | 原子网格的径向格点数 |

适用场景：全电子计算、内层能级谱（XAS）、NMR 化学位移/超精细耦合等依赖近核密度的性质、部分小芯赝势。**常规的价电子 DFT 用 GPW 更简单、更快。**

### 2.4 GTH 赝势

CP2K 主要使用 Goedecker–Teter–Hutter（GTH）型模守恒、可分离、“双空间”赝势。局域部分：

$$
V_{\mathrm{loc}}(r)=-\frac{Z_{\mathrm{ion}}}{r}\,\mathrm{erf}\!\left(\frac{\bar r}{\sqrt2}\right)
+\exp\!\left(-\frac{\bar r^{2}}{2}\right)\left[C_1+C_2\bar r^{2}+C_3\bar r^{4}+C_4\bar r^{6}\right],
\qquad \bar r=\frac{r}{r_{\mathrm{loc}}}
$$

非局域部分是各角动量通道上的可分离投影子：

$$
V_{\mathrm{nl}}=\sum_{l}\sum_{m=-l}^{l}\sum_{i,j}\lvert p_i^{lm}\rangle\,h_{ij}^{l}\,\langle p_j^{lm}\rvert,
\qquad p_i^{l}(r)\propto r^{\,l+2(i-1)}\exp\!\left(-\frac{r^2}{2r_l^2}\right)
$$

所有部分都是高斯型函数，在实空间和倒空间都有解析形式——这正是它与 GPW 天然契合的原因。$Z_{\mathrm{ion}}$ 就是名字里的 `q`：`GTH-PBE-q6` 表示氧有 6 个价电子（2s²2p⁴）被显式处理。

**大芯与小芯（半芯）**：同一元素常有多个 `q`：

| 元素 | 可选 q（`POTENTIAL_UZH`） | 含义 | 建议 |
|---|---|---|---|
| Ti | q4、q12 | q4：3d²4s²；q12：额外包含 3s²3p⁶ 半芯 | 氧化物、磁性、高压体系用 q12（UZH 基组只提供 Ti-q12） |
| Cu | q1、q11、q19 | q1：仅 4s；q11：3d¹⁰4s¹；q19：再加 3s3p | 绝大多数情况用 q11；q1 只适合特殊用途 |

半芯赝势更准、更“硬”，但需要更大的基组和更高的 `CUTOFF`。

---

## 3. 读懂数据文件：基组与赝势条目的格式

### 3.1 基组条目

以 `BASIS_MOLOPT` 中氧的 `DZVP-MOLOPT-GTH` 为例（摘录头部和第一行数据）：

```text
 O  DZVP-MOLOPT-GTH DZVP-MOLOPT-GTH-q6
 1
 2 0 2 7 2 2 1
     12.015954705512 -0.060190841200  0.065738617900  0.036543638800 ...
     ...（共 7 行：每行 1 个指数 + 5 列收缩系数）
```

| 行 | 内容 | 本例解读 |
|---|---|---|
| 1 | 元素、基组名、别名（可多个） | 名字 `DZVP-MOLOPT-GTH`，别名 `DZVP-MOLOPT-GTH-q6` |
| 2 | 基组包含的“集合”（set）数 | 1 个集合 |
| 3 | `n lmin lmax nexp nshell(lmin) … nshell(lmax)` | 主量子数标记 2，角动量 0–2，7 个高斯指数，s/p/d 各收缩成 2/2/1 个函数 |
| 之后 | 每行一个指数及其在各收缩函数中的系数 | 7 行 × (1 + 5) 列 |

“2 个 s + 2 个 p + 1 个 d” 正是 **D**ouble-**Z**eta **V**alence + **P**olarization 的含义。CP2K 使用球谐高斯函数，氧原子的 DZVP 共 $2\times1+2\times3+1\times5=13$ 个基函数。几个对比：

| 基组（氧） | 文件 | 第 3 行 | 壳层 | 最小指数（越小越弥散） |
|---|---|---|---|---|
| `DZVP-MOLOPT-GTH` | `BASIS_MOLOPT` | `2 0 2 7 2 2 1` | 2s2p1d | 0.0468 |
| `DZVP-MOLOPT-SR-GTH` | `BASIS_MOLOPT` | `2 0 2 5 2 2 1` | 2s2p1d | 0.162 |
| `DZVP-MOLOPT-PBE-GTH-q6` | `BASIS_MOLOPT_UZH` | `2 0 2 5 2 2 1` | 2s2p1d | 0.155 |
| `TZV2P-MOLOPT-GTH` | `BASIS_MOLOPT` | `2 0 2 7 3 3 2` | 3s3p2d | — |
| `TZV2P-MOLOPT-PBE-GTH-q6` | `BASIS_MOLOPT_UZH` | `2 0 3 5 3 3 2 1` | 3s3p2d1f | — |

可以看出：经典 MOLOPT 的原始高斯更多、更弥散（适合气相分子、弱相互作用，但凝聚相中重叠矩阵更大、更病态）；SR（Short Range）与 UZH 版本更紧凑，更适合凝聚相；UZH 的 TZV2P 还多了 f 极化函数。

### 3.2 赝势条目

`GTH_POTENTIALS` 中的 `O GTH-PBE-q6`：

```text
O GTH-PBE-q6 GTH-PBE
    2    4
     0.24455430    2   -16.66721480     2.48731132
    2
     0.22095592    1    18.33745811
     0.21133247    0
```

| 行 | 内容 | 本例解读 |
|---|---|---|
| 1 | 元素、赝势名、别名 | 别名 `GTH-PBE` |
| 2 | 各角动量的价电子数 | 2 个 s + 4 个 p，$Z_{\mathrm{ion}}=6$ |
| 3 | $r_{\mathrm{loc}}$、局域系数个数、$C_1, C_2,\dots$ | $r_{\mathrm{loc}}=0.2446$ bohr，2 个系数 |
| 4 | 非局域通道数 | 2 个通道（s、p） |
| 5– | 每个通道：$r_l$、投影子个数、$h_{ij}$ | s 通道 1 个投影子；p 通道 0 个投影子 |

`POTENTIAL_UZH` 中同名的 `O GTH-PBE-q6`（别名 `GTH-GGA-q6`）：

```text
O  GTH-PBE-q6 GTH-GGA-q6
    2    4    0    0
    0.24446328480160       2  -16.67548222363837    2.48908598241780
       1
    0.22097110943471       1   18.33446866406285
```

> [!WARNING]
> **同名不同参**：两个文件里都叫 `GTH-PBE-q6`，但 UZH 版本是重新优化过的（参数略有不同，且只有 1 个非局域通道）。UZH 基组是针对 UZH 赝势优化的，经典 MOLOPT 基组是针对 `GTH_POTENTIALS` 优化的。**基组与赝势要来自同一个“家族”**，只看名字一致是不够的。

---

## 4. 基组家族与命名规则

### 4.1 名字的各个部分

以 `TZV2P-MOLOPT-PBE-GTH-q6` 为例：

| 片段 | 含义 |
|---|---|
| `SZV` / `DZVP` / `TZVP` / `TZV2P` / `TZV2PX` / `QZVPP` | 质量等级：单/双/三 zeta 价层，P 为极化函数（2P 为两组），X 为额外函数 |
| `MOLOPT` | 在分子集合上优化的收缩基组（VandeVondele & Hutter 2007），条件数好，气相/凝聚相通用 |
| `SR` | Short-Range：更少、更不弥散的原始高斯，凝聚相大体系更高效（BSSE 略增） |
| `PBE` / `SCAN` / `PBE0`（别名 `GGA` / `MGGA` / `HYB`） | UZH 基组针对的泛函家族 |
| `GTH` | 配合 GTH 赝势使用 |
| `-qN` | 对应赝势的价电子数 N |
| `ae` | all-electron，全电子基组（配 `POTENTIAL ALL` 与 GAPW） |

### 4.2 主要基组文件

| 文件 | 典型名字 | 覆盖范围 | 配对的赝势文件 |
|---|---|---|---|
| `BASIS_MOLOPT` | `SZV/DZVP/TZVP/TZV2P/TZV2PX-MOLOPT-GTH`；`SZV/DZVP-MOLOPT-SR-GTH` | 完整 MOLOPT 系列只有 H、C、N、O、F、Si、P、S、Cl；SR 系列（SZV、DZVP）覆盖周期表大部分元素 | `GTH_POTENTIALS` |
| `BASIS_MOLOPT_UZH` | `DZVP/TZVP/TZV2P-MOLOPT-{PBE,SCAN,PBE0}-GTH-qN`；全电子 `SVP/TZVPP/QZVPP-MOLOPT-PBE-ae` | 大部分元素，并按泛函家族分别优化 | `POTENTIAL_UZH`（全电子用 `POTENTIAL ALL`） |
| `BASIS_ADMM_UZH` | `admm-dz/dzp/tzp/tz2p-qN` | ADMM 辅助基组 | 与 UZH 轨道基组配合 |
| `BASIS_ADMM` | `cFIT3`、`cpFIT3`、`aug-cFIT3` 等 | 经典 ADMM 辅助基组 | 与经典 MOLOPT 配合 |
| `BASIS_RI_MOLOPT`、`BASIS_RI_AUG_MOLOPT` 等 | `RI_...` | RI-MP2/RPA/RI-HFX 辅助基组 | — |
| `BASIS_AUG_MOLOPT`（2026.1 起） | 带弥散函数的 MOLOPT | 阴离子、激发态、弱相互作用 | — |
| `EMSL_BASIS_SETS`、`ALL_BASIS_SETS` | `6-31Gx`（= `6-31G*`）、`6-311++Gxx` 等 | 传统全电子基组 | `POTENTIAL ALL` |
| `GTH_BASIS_SETS`、`BASIS_SET` | 老一代 GTH 基组 | 历史输入 | — |

> [!NOTE]
> 官方文档对新项目的建议是：**当元素和泛函家族有对应条目时，优先使用 UZH 协议的成对基组/赝势**（`BASIS_MOLOPT_UZH` + `POTENTIAL_UZH`）；经典 `BASIS_MOLOPT` + `GTH_POTENTIALS` 仍然重要，用于复现已发表结果或 UZH 尚未覆盖的情况。UZH 协议的系统验证见 Mirhosseini *et al.*, *J. Chem. Phys.* (2026)。

### 4.3 查找某元素有哪些基组和赝势

```bash
cd $CP2K_DATA_DIR
grep -E "^ *Ti " BASIS_MOLOPT                    # Ti: 只有 SZV/DZVP-MOLOPT-SR-GTH(-q12)
grep -E "^Ti " BASIS_MOLOPT_UZH | grep -v "\-ae" # Ti 的 UZH 基组（均为 q12）
grep -E "^Ti " POTENTIAL_UZH                     # Ti 的 UZH 赝势：q4 与 q12（PBE/SCAN/PBE0）
grep -E "^Ti GTH-PBE" GTH_POTENTIALS             # 经典：GTH-PBE-q12
grep -E "^ *O +6-31" EMSL_BASIS_SETS             # 氧的 Pople 基组名
```

图形化/网页工具：[cp2k-basis](https://www.cp2k.org/tools:cp2k-basis)（在线版 <https://cp2k-basis.pierrebeaujean.net>）可以按元素查出匹配的基组和赝势组合。

---

## 5. 赝势库与选择原则

| 文件 | 内容 | 用途 |
|---|---|---|
| `GTH_POTENTIALS` | 经典 GTH 赝势（`GTH-PBE-qN`、`GTH-BLYP-qN` 等，部分有短别名如 `GTH-PBE`） | 配合 `BASIS_MOLOPT` |
| `POTENTIAL` | 更大的合集（多种泛函的 GTH 赝势与 `ALL`） | `POTENTIAL_FILE_NAME` 的默认值 |
| `POTENTIAL_UZH`（+ `_2026.1`、`_2026.2` 快照） | UZH 协议赝势：`GTH-PBE-qN`、`GTH-SCAN-qN`、`GTH-PBE0-qN`（别名 `GGA`/`MGGA`/`HYB`）及 `ALL` | 配合 UZH 基组 |
| `NLCC_POTENTIALS` | 含非线性芯修正的赝势 | 专门用途 |
| `GTH_SOC_POTENTIALS` | 含自旋–轨道耦合参数的赝势 | SOC 计算 |
| `ECP_POTENTIALS` | 高斯积分型有效芯势（ECP） | 与 ECP 基组配合 |
| `ALL_POTENTIALS` | 全电子“势”定义（`ALLELECTRON ALL`） | GAPW 全电子 |

选择原则：

1. **赝势应与所用泛函属于同一家族**：PBE 计算用 `GTH-PBE-qN`；SCAN 用 `GTH-SCAN-qN`；PBE0/杂化用 `GTH-PBE0-qN`。混用（例如 PBE0 计算用 PBE 赝势）在探索性计算中有时可接受，但不是系统提高精度的途径。
2. **按性质决定芯大小**：涉及半芯态参与成键、磁性、高压、内层激发时选小芯（大 q）。
3. **与基组成对**：`-qN` 一致，且来自同一家族的文件（§3.2）。

> [!IMPORTANT]
> **`&DFT` 中的 `POTENTIAL_FILE_NAME` 只能写一次**（重复会报 `Keyword POTENTIAL_FILE_NAME in section DFT should not repeat.`），而 `BASIS_SET_FILE_NAME` 可以写多行。如果某个元素的赝势必须来自另一个文件，在该元素的 `&KIND` 中单独写 `POTENTIAL_FILE_NAME`：
>
> ```text
> &KIND Cu
>   BASIS_SET DZVP-MOLOPT-SR-GTH
>   POTENTIAL GTH-PBE-q11
>   POTENTIAL_FILE_NAME GTH_POTENTIALS      ! overrides the file given in &DFT for this kind only
> &END KIND
> ```

---

## 6. -qN 匹配规则与自检清单

| 检查项 | 怎么做 |
|---|---|
| 基组名真实存在于所列文件中 | `grep -E "^ *O .*DZVP-MOLOPT-PBE-GTH-q6" $CP2K_DATA_DIR/BASIS_MOLOPT_UZH` |
| 赝势名真实存在于赝势文件中 | `grep -E "^O +GTH-PBE-q6" $CP2K_DATA_DIR/POTENTIAL_UZH` |
| 基组后缀 `-qN` = 赝势 `-qN` | 经典 MOLOPT 基组名不带 q 时，看它的别名（如 `DZVP-MOLOPT-GTH-q6`） |
| 赝势与泛函家族一致 | PBE↔`GTH-PBE`，SCAN↔`GTH-SCAN`，PBE0↔`GTH-PBE0` |
| 基组与赝势来自同一家族 | UZH↔UZH，MOLOPT↔`GTH_POTENTIALS` |
| 运行后核对电子数 | 输出 `Number of electrons` / `Electronic density on regular grids` 与 Σq 一致 |

> [!CAUTION]
> CP2K 会在基组或赝势**找不到**时报错，但**不会替你检查** `-qN` 是否一致、家族是否匹配。用 q11 的基组配 q19 的赝势，程序照样会跑，结果却是错的（基组缺少描述半芯态的函数）。上表的检查要靠自己完成。

---

## 7. 实战一：MOLOPT 基组阶梯（h2o_szv_dzvp_tzv2p.inp）

### 7.1 输入要点

[`h2o_szv_dzvp_tzv2p.inp`](../examples/03-basis-pseudo-h2o/h2o_szv_dzvp_tzv2p.inp) 用预处理器变量拼出基组名，并打印偶极矩作为可观测量：

```text
&GLOBAL
  PROJECT h2o_${BASIS-DZVP}-${FAMILY-MOLOPT}
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD QS
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 500
      REL_CUTOFF 60
      NGRIDS ${NGRIDS-4}
    &END MGRID
    ...
    &PRINT
      &MOMENTS ON             ! dipole moment, printed in Debye
        PERIODIC .FALSE.      ! simple position operator; requires a non-periodic cell
        REFERENCE COM         ! centre of mass as origin
      &END MOMENTS
    &END PRINT
  &END DFT
  &SUBSYS
    ...
    &KIND O
      BASIS_SET ${BASIS-DZVP}-${FAMILY-MOLOPT}-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET ${BASIS-DZVP}-${FAMILY-MOLOPT}-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 变量/关键字 | 取值 | 说明 |
|---|---|---|
| `BASIS` | `SZV`、`DZVP`、`TZVP`、`TZV2P`、`TZV2PX` | 经典 MOLOPT 的全部等级（H、O 都有） |
| `FAMILY` | `MOLOPT`、`MOLOPT-SR` | SR 只有 SZV 与 DZVP；`TZV2P-MOLOPT-SR-GTH` 不存在，会报基组找不到 |
| `CUTOFF`/`REL_CUTOFF` | 500/60 Ry | 比[第2章](02-input-structure.md)更保守，使网格误差远小于基组差别 |
| `&MOMENTS PERIODIC .FALSE.` | — | 孤立分子用普通位置算符计算偶极矩；周期体系默认 `PERIODIC T`（Berry 相位公式） |
| `REFERENCE COM` | — | 以质心为原点（对中性分子偶极矩与原点无关，但带电体系有关） |

### 7.2 运行

```bash
cd CP2K/examples/03-basis-pseudo-h2o
bash run_basis_scan.sh
```

脚本对 7 个基组依次运行，并从输出中提取 `ENERGY| Total FORCE_EVAL`、`Dipole moment [Debye]` 下一行的 `Total=`、SCF 步数和 timing 表中 `CP2K` 行的总耗时。输出表头如下（示意输出，数值未给出，请以实际运行为准）：

```text
basis                           E_total[Ha]    dipole[D]      SCF    time[s]
SZV-MOLOPT-GTH                          ...          ...      ...        ...
DZVP-MOLOPT-GTH                         ...          ...      ...        ...
...
```

### 7.3 如何解读

- **总能量随基组单调下降**（变分原理）；SZV → DZVP 降幅最大，TZV2P → TZV2PX 的变化应很小，说明接近该泛函下的基组极限。
- **只比较同一家族、同一赝势下的总能量**。更有意义的是能量差（反应能、结合能）、结构和偶极矩等可观测量。水分子气相偶极矩实验值约 1.85 D；SZV 往往偏差明显，DZVP 以上逐渐稳定（PBE 本身也有误差，不要求与实验完全一致）。
- **耗时随基函数数增长很快**：氧的基函数数从 SZV 的 4 个、DZVP 的 13 个到 TZV2P 的 22 个。凝聚相体系中这种差别会被放大。
- **SR 版本**：`DZVP-MOLOPT-SR-GTH` 的能量略高于 `DZVP-MOLOPT-GTH`（少了弥散函数），但对凝聚相大体系通常更快、更稳定。
- 可以再用 `-E NGRIDS=5` 跑一次 `DZVP-MOLOPT-GTH`，比较 `count for grid` 分布与耗时（§2.2 的提示）。

> [!TIP]
> 高斯基组存在基组重叠误差（BSSE）：二聚体中的单体会“借用”另一个单体的基函数而人为降低能量。计算弱相互作用结合能时，基组越小 BSSE 越大；可用 `&KIND` 中的 `GHOST` 鬼原子做 counterpoise 校正，或换更大的基组检查结合能是否收敛。

---

## 8. 实战二：UZH 基组与赝势成对使用（h2o_uzh_pair.inp）

### 8.1 输入要点

[`h2o_uzh_pair.inp`](../examples/03-basis-pseudo-h2o/h2o_uzh_pair.inp)：

```text
  &DFT
    BASIS_SET_FILE_NAME ${BFILE-BASIS_MOLOPT_UZH}
    POTENTIAL_FILE_NAME ${PFILE-POTENTIAL_UZH}
    ...
  &SUBSYS
    &KIND O
      BASIS_SET ${Q-DZVP}-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET ${Q-DZVP}-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
```

UZH 命名中的泛函标签与别名完全等价，下面两种写法读到的是同一个基组和赝势：

```text
BASIS_SET DZVP-MOLOPT-PBE-GTH-q6      POTENTIAL GTH-PBE-q6
BASIS_SET DZVP-MOLOPT-GGA-GTH-q6      POTENTIAL GTH-GGA-q6
```

对应 SCAN 计算写 `DZVP-MOLOPT-SCAN-GTH-q6` + `GTH-SCAN-q6`（别名 `MGGA`），PBE0 计算写 `DZVP-MOLOPT-PBE0-GTH-q6` + `GTH-PBE0-q6`（别名 `HYB`）。

### 8.2 可重复性：带版本号的数据快照（2026.2 新增）

CP2K 2026.2 的 `data/` 目录中新增了：

| 文件 | 内容 |
|---|---|
| `BASIS_MOLOPT_UZH_2026.1`、`POTENTIAL_UZH_2026.1` | 与 2025.2/2026.1 发布的 `BASIS_MOLOPT_UZH`/`POTENTIAL_UZH` 内容一致 |
| `BASIS_MOLOPT_UZH_2026.2`、`POTENTIAL_UZH_2026.2` | 2026.2 的新版本 |
| `BASIS_MOLOPT_UZH`、`POTENTIAL_UZH` | 在 2026.2 源码树中是指向 `*_2026.2` 的链接，即“最新版” |

根据 2026.2 数据文件头部的更新记录，新快照替换了 Ca、Cr、Tc、Ru、Rh、Pd、Ba、Ir、Au 以及稀有气体（He、Ne、Ar、Kr、Xe、Rn）等元素的改进基组/赝势。这意味着：**同一个输入，用 2026.1 和 2026.2 运行，若体系含有这些元素，读到的基组/赝势并不相同，总能量不能直接比较。**

```bash
# 默认：读“最新版”（在 2026.2 中等于 _2026.2 快照）
cp2k.ssmp -i h2o_uzh_pair.inp -o h2o_uzh.out
# 显式冻结到 2026.2 快照：将来换新版本 CP2K 也读同一份数据（文件只在 >= 2026.2 中存在）
cp2k.ssmp -E BFILE=BASIS_MOLOPT_UZH_2026.2 -E PFILE=POTENTIAL_UZH_2026.2 \
          -i h2o_uzh_pair.inp -o h2o_uzh_2026.2.out
# 在 2026.2 中复现 2026.1 时期的结果
cp2k.ssmp -E BFILE=BASIS_MOLOPT_UZH_2026.1 -E PFILE=POTENTIAL_UZH_2026.1 \
          -i h2o_uzh_pair.inp -o h2o_uzh_2026.1.out
```

对 H、O 而言两个快照没有差别，三次计算的能量应一致；换成含 Au、Pd 的体系就能看到差别。

> [!TIP]
> 写论文或做长期项目时，建议在输入中直接写快照文件名（如 `BASIS_MOLOPT_UZH_2026.2`），并在方法部分注明。另一个保险的办法是：CP2K 写出的 `*-1.restart` 文件中，每个 `&KIND` 都附带了完整展开的 `&BASIS`/`&POTENTIAL` 数据，保留它就保留了实际使用的参数。

---

## 9. 实战三：GAPW 全电子计算（h2o_gapw_allelectron.inp）

### 9.1 输入要点

[`h2o_gapw_allelectron.inp`](../examples/03-basis-pseudo-h2o/h2o_gapw_allelectron.inp)：

```text
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    BASIS_SET_FILE_NAME EMSL_BASIS_SETS
    POTENTIAL_FILE_NAME ALL_POTENTIALS
    &MGRID
      CUTOFF 800              ! the soft density still needs a fine grid in GAPW
      REL_CUTOFF 80
    &END MGRID
    &QS
      METHOD GAPW             ! default is GPW
      EPS_DEFAULT 1.0E-12
      EPSFIT 1.0E-6           ! hard/soft split of the Gaussian exponents (default 1.0E-4)
      EPSRHO0 1.0E-8          ! range of the hard compensation density (default 1.0E-6)
      EPSSVD 1.0E-10          ! SVD tolerance for the projectors (default 1.0E-8)
    &END QS
    ...
  &SUBSYS
    &KIND O
      BASIS_SET ${AEB-SVP-MOLOPT-PBE-ae}
      POTENTIAL ALL           ! full nuclear potential, all 8 electrons
      LEBEDEV_GRID 110        ! angular points of the atom-centred grid (default 50)
      RADIAL_GRID 80          ! radial points of the atom-centred grid (default 50)
    &END KIND
    &KIND H
      BASIS_SET ${AEB-SVP-MOLOPT-PBE-ae}
      POTENTIAL ALL
      LEBEDEV_GRID 110
      RADIAL_GRID 80
    &END KIND
```

| 关键字 | 本例取值 | 说明 |
|---|---|---|
| `METHOD GAPW` | — | 启用 GAPW；全电子计算必需 |
| `POTENTIAL ALL` | — | 全电子“势”（裸核库仑势），从 `ALL_POTENTIALS`（或 `POTENTIAL`、`POTENTIAL_UZH` 中的 `ALLELECTRON ALL` 条目）读取 |
| `SVP-MOLOPT-PBE-ae` | 默认 | UZH 全电子基组（别名 `SVP-MOLOPT-GGA-ae`），基于 def2-SVP 构建；更大的有 `TZVPP-…-ae`、`QZVPP-…-ae` |
| `6-31Gx` | `-E AEB=6-31Gx` | 即 6-31G*，在 `EMSL_BASIS_SETS` 中 |
| `CUTOFF 800`、`EPSFIT` 等 | 与官方 GAPW 文档的测试示例同一量级 | 先保证正确，再按需放宽 |
| `LEBEDEV_GRID`/`RADIAL_GRID` | 110/80 | 原子网格越密，近核性质越准，开销越大 |

### 9.2 运行与检查

```bash
cd CP2K/examples/03-basis-pseudo-h2o
rm -f h2o_gapw*.out
cp2k.ssmp -i h2o_gapw_allelectron.inp -o h2o_gapw_svp.out
cp2k.ssmp -E AEB=6-31Gx -i h2o_gapw_allelectron.inp -o h2o_gapw_631gs.out

grep "Total charge density (r-space)" h2o_gapw_svp.out
grep "Total Rho_soft + Rho1_hard - Rho1_soft (r-space)" h2o_gapw_svp.out
grep "ENERGY| Total FORCE_EVAL" h2o_gapw_*.out
```

- **电子数检查**：官方文档建议在 SCF 收敛后检查 CP2K 打印的电子数，它是软/硬划分质量的诊断。水分子全电子共 10 个电子，上面两行电荷积分应非常接近 −10；偏差较大时提高原子网格或收紧 `EPSFIT`/`EPSRHO0`/`EPSSVD`，并相应提高 `CUTOFF`。
- **总能量参考点不同**：全电子总能量约在 −76 Ha 量级，赝势计算约 −17 Ha；两者绝不能相减或比较。
- **与量子化学程序对比**：CP2K 使用球谐 d 函数（5d），而 Gaussian 对 6-31G* 默认使用笛卡尔 6d，因此同名基组的能量不会与 Gaussian 逐位一致（Gaussian 中可用 `5d` 关键字对齐）。

> [!NOTE]
> GAPW 也可以搭配小芯 GTH 赝势使用，或者用 `METHOD GAPW_XC` 只对交换关联项做 GAPW 分解。内层激发（XAS）、NMR/EPR 等性质计算的具体输入不在本系列主线中，可参考官方方法文档 <https://manual.cp2k.org/trunk/methods/dft/gapw.html> 及其中链接的性质章节。

---

## 10. 多个基组文件、辅助基组与 KIND 模板

### 10.1 多个 BASIS_SET_FILE_NAME

`BASIS_SET_FILE_NAME` 可以重复，CP2K 会依次读取所有文件来查找每个 `&KIND` 要求的基组：

```text
&DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH     ! orbital basis sets
  BASIS_SET_FILE_NAME BASIS_ADMM_UZH       ! ADMM auxiliary basis sets
  BASIS_SET_FILE_NAME ./MY_BASIS           ! a private file in the working directory
  POTENTIAL_FILE_NAME POTENTIAL_UZH
  ...
&END DFT
```

找不到时，报错信息会列出所有被搜索的文件：

```text
The requested basis set <TZV2P-MOLOPT-SR-GTH> for element <O> was not found in the basis set files <BASIS_MOLOPT>
```

> [!TIP]
> 不要让同一个基组名出现在多个被读取的文件里（尤其是自己复制修改过的文件），以免搞不清实际用的是哪一份。用 `cp2k.psmp -e` 或输出中的基组信息块核对。

### 10.2 一个 KIND 中的多个基组：ORB、AUX_FIT、RI_AUX

`BASIS_SET` 关键字可重复，前面加“类型”指明用途；不写类型即 `ORB`（轨道基组）：

| 写法 | 用途 | 详见 |
|---|---|---|
| `BASIS_SET DZVP-MOLOPT-PBE0-GTH-q6`（= `BASIS_SET ORB ...`） | 轨道基组 | 本章 |
| `BASIS_SET AUX_FIT admm-dzp-q6` | ADMM 辅助基组（来自 `BASIS_ADMM_UZH`；经典组合可用 `BASIS_ADMM` 中的 `cFIT3` 等） | [第9章](09-hybrid-admm.md) |
| `BASIS_SET RI_AUX <名称>` | RI-MP2/RPA/RI-HFX 的辅助基组（`BASIS_RI_MOLOPT` 等文件） | 官方 RI 文档 |
| `BASIS_SET LRI <名称>` | LRI-GPW 的局域 RI 基组 | 官方文档 |

示例（PBE0 + ADMM 计算中的氧）：

```text
&KIND O
  BASIS_SET DZVP-MOLOPT-PBE0-GTH-q6
  BASIS_SET AUX_FIT admm-dzp-q6
  POTENTIAL GTH-PBE0-q6
&END KIND
```

### 10.3 用 @INCLUDE 复用 KIND 模板

[`kind_templates.inc`](../examples/03-basis-pseudo-h2o/kind_templates.inc) 把 H、C、N、O 的 `&KIND` 写成模板，用变量控制基组等级、泛函家族和是否加 ADMM 辅助基组：

```text
&KIND O
  BASIS_SET ${Q-DZVP}-MOLOPT-${XCTAG-PBE}-GTH-q6
  POTENTIAL GTH-${XCTAG-PBE}-q6
  @IF ${ADMM-0} == 1
  BASIS_SET AUX_FIT admm-dzp-q6
  @ENDIF
&END KIND
```

在主输入的 `&SUBSYS` 中写 `@INCLUDE 'kind_templates.inc'` 即可（见 [`ch3oh_kind_templates.inp`](../examples/03-basis-pseudo-h2o/ch3oh_kind_templates.inp)）。结构中不存在的元素对应的 `&KIND` 段不会被使用（CP2K 只为实际出现的原子类型查找匹配的 `&KIND`）。

```bash
cp2k.ssmp -e ch3oh_kind_templates.inp | grep -A3 "&KIND"            # 默认 DZVP + PBE
cp2k.ssmp -E Q=TZV2P -e ch3oh_kind_templates.inp | grep BASIS_SET    # 换成 TZV2P
cp2k.ssmp -E XCTAG=PBE0 -E ADMM=1 -e ch3oh_kind_templates.inp | grep -E "BASIS_SET|POTENTIAL"
```

最后一条命令只用于查看展开结果；真正做 PBE0 + ADMM 计算还需要第9章的 `&XC`/`&HF`/`&AUXILIARY_DENSITY_MATRIX_METHOD` 设置。

---

## 11. 基组选择流程

```text
1. 确定泛函（PBE / SCAN / PBE0 / HSE06 ...）
        |
2. 选家族：新项目 -> UZH 成对（BASIS_MOLOPT_UZH + POTENTIAL_UZH，标签与泛函一致）
           复现文献 -> 与原文相同的文件和名字
        |
3. 为每个元素选 q（芯大小）：涉及半芯态/磁性/高压时取小芯（大 q）
        |
4. 选等级：结构优化、AIMD、大体系 -> DZVP
           能量精修（反应能、吸附能、能垒）-> 在 DZVP 结构上做 TZV2P 单点，或检查 TZV2P 下结论不变
        |
5. 用 grep、cp2k-basis 工具或本节脚本 cp2k_inpgen.py 确认名字存在、-qN 一致
        |
6. 针对选定基组做 CUTOFF / REL_CUTOFF 收敛（第4章），基组变了要重新检查
```

| 任务 | 推荐起点 | 备注 |
|---|---|---|
| 水/溶液/界面 AIMD | `DZVP-MOLOPT-PBE-GTH-qN`（或 `DZVP-MOLOPT-SR-GTH`）+ D3 | 大体系耗时主要在 SCF，DZVP 性价比最高 |
| 固体、表面吸附 | DZVP（SR 或 UZH）；吸附能用 TZV2P 检验 | 注意 BSSE |
| 分子反应能、能垒 | DZVP 优化 + TZV2P 单点 | 与小分子量化方法对比时考虑全电子基组 |
| 弱相互作用、阴离子 | TZV2P 以上，或带弥散函数的基组（`BASIS_AUG_MOLOPT`，2026.1 起） | 弥散函数使重叠矩阵病态，凝聚相慎用 |
| 杂化泛函 | `…-PBE0-GTH-qN` + `GTH-PBE0-qN` + `admm-dzp-qN` | [第9章](09-hybrid-admm.md) |
| 全电子/近核性质 | GAPW + `TZVPP-MOLOPT-PBE-ae` 等 | §9 |

本节的输入生成脚本 [`cp2k_inpgen.py`](../scripts/cp2k_inpgen.py) 可以从 xyz/cif/POSCAR 结构出发，按元素查找 `BASIS_MOLOPT_UZH`/`POTENTIAL_UZH` 中的条目并核对 `-qN` 是否一致（用 `python ../scripts/cp2k_inpgen.py -h` 查看用法）；生成的输入仍应按本章 §6 的清单人工复核。

---

## 参数选择与收敛建议

| 参数 | 建议 | 说明 |
|---|---|---|
| 基组等级 | 结构/MD：DZVP；能量：TZV2P 检验 | 关键结论应在两个等级下一致 |
| 基组家族 | 新项目 UZH 成对；复现文献沿用原组合 | 不要跨家族混配基组与赝势 |
| 数据文件版本 | 长期项目写快照名（`*_2026.2`） | 升级 CP2K 时保持数据不变 |
| `REL_CUTOFF` | 50–60 Ry 起步 | 默认 20 Ry 偏低 |
| `CUTOFF` | 小芯赝势、硬元素（O、F、过渡金属半芯）需要更高 | 基组/赝势改变后重新收敛（[第4章](04-cutoff-convergence.md)） |
| `NGRIDS` | 4（默认）；弥散的经典 MOLOPT 可试 5 | 改动后核对能量 |
| GAPW 原子网格 | `LEBEDEV_GRID` ≥ 110、`RADIAL_GRID` ≥ 80 起步 | 以电子数积分与目标性质收敛为准 |
| GAPW `EPSFIT`/`EPSRHO0`/`EPSSVD` | 从官方测试示例的量级（1E-6/1E-8/1E-10）出发 | 收紧后提高 `CUTOFF` |

---

## 常见问题与报错

| 报错或现象 | 原因 | 解决 |
|---|---|---|
| `The requested basis set <DZVP-MOLOPT-GTH> for element <Cu> was not found in the basis set files <BASIS_MOLOPT>` | 该元素在此文件中没有这个基组（完整 MOLOPT 只覆盖少数元素） | 改用 `DZVP-MOLOPT-SR-GTH`（`BASIS_MOLOPT`）或 UZH 基组（`BASIS_MOLOPT_UZH`）；`grep` 确认名字 |
| `The requested atomic potential <GTH-PBE-q11> for element <Cu> was not found in the potential file <...>` | 赝势名或文件不对 | 核对 `POTENTIAL_FILE_NAME`；查该元素可用的 q（Cu 在 `POTENTIAL_UZH` 中有 q1/q11/q19） |
| `Keyword POTENTIAL_FILE_NAME in section DFT should not repeat.` | `&DFT` 中写了两行 `POTENTIAL_FILE_NAME` | 只保留一行；个别元素用 `&KIND` 中的 `POTENTIAL_FILE_NAME` |
| `The specified OLD file <BASIS_MOLOPT_UZH_2026.2> cannot be opened. It does not exist.` | 快照文件只在 CP2K ≥ 2026.2 的数据目录中 | 旧版本改用不带版本号的文件名，或指向 2026.2 的数据目录 |
| 能跑但结果异常（能量偏高很多、结构崩坏） | `-qN` 不一致，或基组/赝势来自不同家族 | 按 §6 清单逐项检查 |
| 升级 CP2K 后，含 Au/Pd 等元素的体系能量变了 | 2026.2 的 `BASIS_MOLOPT_UZH`/`POTENTIAL_UZH` 更新了部分元素 | 用 `*_2026.1` 快照复现旧结果，或统一用新数据重算 |
| 大基组 SCF 难收敛、出现线性相关警告 | 弥散函数在凝聚相中使重叠矩阵病态 | 换 SR/UZH 基组；必要时收紧 `EPS_DEFAULT`；SCF 策略见 [第5章](05-scf-strategy.md) |
| GAPW 电子数积分明显偏离总电子数 | 原子网格太粗或软/硬划分参数太松 | 提高 `LEBEDEV_GRID`/`RADIAL_GRID`，收紧 `EPSFIT` 等并提高 `CUTOFF` |
| 全电子计算忘了 `METHOD GAPW` | 用 GPW 处理全电子密度，需要极高截断，结果不收敛或极慢 | 在 `&QS` 中写 `METHOD GAPW` |
| 与 Gaussian 程序同名基组的能量对不上 | 球谐（5d）与笛卡尔（6d）d 函数的差别、泛函实现与积分网格的差别 | 对齐 5d/6d 设置；比较能量差而不是总能量 |

---

## 练习

1. 运行 `run_basis_scan.sh`，作出“总能量 vs 基函数数”和“偶极矩 vs 基组等级”两张图，判断 H2O 偶极矩在哪个等级收敛到 0.02 D 以内。
2. 对 `DZVP-MOLOPT-GTH` 分别用 `NGRIDS 4` 和 `NGRIDS 5` 运行，比较 `count for grid` 分布、总能量与耗时。
3. 用 `h2o_uzh_pair.inp` 比较 UZH `DZVP`/`TZVP`/`TZV2P` 的偶极矩，并与经典 MOLOPT 的结果对比；思考为什么两个家族的总能量不能直接比较。
4. 在 `$CP2K_DATA_DIR` 中查找 Fe、Cu、Pt 各有哪些 q 的 GTH 赝势（`GTH_POTENTIALS` 与 `POTENTIAL_UZH`），以及对应的 UZH 基组，写出它们的 `&KIND` 段。
5. 用 GAPW 输入分别跑 `SVP-MOLOPT-PBE-ae`、`TZVPP-MOLOPT-PBE-ae` 和 `6-31Gx`，比较总能量、偶极矩和电子数积分；再把 `LEBEDEV_GRID` 降到 50，观察电子数积分的变化。
6. （进阶）把 `kind_templates.inc` 扩展到 Si 和 Ti（UZH 中 Si 为 q4，Ti 的 UZH 基组只有 q12），并用 `cp2k.psmp -e` 确认展开结果。

---

## 参考资料

**官方文档**

- 基组：<https://manual.cp2k.org/trunk/methods/dft/basis_sets.html>
- 赝势（GTH、POTENTIAL_UZH、ALL）：<https://manual.cp2k.org/trunk/methods/dft/pseudopotentials.html>
- GPW 方法：<https://manual.cp2k.org/trunk/methods/dft/gpw.html>
- GAPW 方法（含测试过的全电子 H2O 示例）：<https://manual.cp2k.org/trunk/methods/dft/gapw.html>
- CUTOFF/REL_CUTOFF 收敛与多重网格：<https://manual.cp2k.org/trunk/methods/dft/cutoff.html>
- 杂化泛函与 ADMM：<https://manual.cp2k.org/trunk/methods/dft/hartree-fock/index.html>
- 基组与赝势查询工具 cp2k-basis：<https://www.cp2k.org/tools:cp2k-basis>；GTH 赝势数据库：<https://cp2k.org/static/potentials/>
- 输入参考 `&KIND`：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html>

**相关 GitHub 仓库**（只链接，不复制；更多见 [Resources/github-repos.md](../../Resources/github-repos.md)）

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`data/` 目录包含 `BASIS_MOLOPT`、`BASIS_MOLOPT_UZH`（及 2026.x 快照）、`POTENTIAL_UZH`、`GTH_POTENTIALS`、`BASIS_ADMM_UZH` 等全部数据文件。
- [juerghutter/BASIS](https://github.com/juerghutter/BASIS)：CP2K 核心开发者 Jürg Hutter 维护的基组集合，可了解 MOLOPT/UZH 基组的来源。
- [DCM-Uni-Paderborn/UZH-protocol](https://github.com/DCM-Uni-Paderborn/UZH-protocol)：UZH 基组与赝势的验证数据，与 `BASIS_MOLOPT_UZH_2026.2` 对齐，说明了新项目推荐 UZH 组合的依据。
- [pierre-24/cp2k-basis](https://github.com/pierre-24/cp2k-basis)：cp2k-basis 网页工具背后的项目，解析基组与 GTH 赝势并查找匹配组合。
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，其中基组章节可对照阅读（GPL-3.0，仅链接）。
- [Clementine-C/cp2kbuild](https://github.com/Clementine-C/cp2kbuild)：按 `-qN` 自动匹配基组与赝势并生成输入的中文工具，可借鉴其思路。
- [JunboLu/CP2K_kit](https://github.com/JunboLu/CP2K_kit)：含 GTH 赝势与基组优化功能的工具集。

**关键论文**

- G. Lippert, J. Hutter, M. Parrinello, *Mol. Phys.* 92, 477 (1997)（GPW 方法）. DOI: [10.1080/002689797170220](https://doi.org/10.1080/002689797170220)
- G. Lippert, J. Hutter, M. Parrinello, *Theor. Chem. Acc.* 103, 124 (1999)（GAPW 方法）. DOI: [10.1007/s002140050523](https://doi.org/10.1007/s002140050523)
- M. Krack, M. Parrinello, *Phys. Chem. Chem. Phys.* 2, 2105 (2000)（全电子 GAPW 分子动力学）. DOI: [10.1039/b001167n](https://doi.org/10.1039/b001167n)
- S. Goedecker, M. Teter, J. Hutter, *Phys. Rev. B* 54, 1703 (1996)（GTH 赝势）. DOI: [10.1103/PhysRevB.54.1703](https://doi.org/10.1103/PhysRevB.54.1703)
- C. Hartwigsen, S. Goedecker, J. Hutter, *Phys. Rev. B* 58, 3641 (1998)（相对论 GTH 赝势 H–Rn）. DOI: [10.1103/PhysRevB.58.3641](https://doi.org/10.1103/PhysRevB.58.3641)
- M. Krack, *Theor. Chem. Acc.* 114, 145 (2005)（针对 GGA 泛函优化的 GTH 赝势）. DOI: [10.1007/s00214-005-0655-y](https://doi.org/10.1007/s00214-005-0655-y)
- J. VandeVondele, J. Hutter, *J. Chem. Phys.* 127, 114105 (2007)（MOLOPT 基组）. DOI: [10.1063/1.2770708](https://doi.org/10.1063/1.2770708)
- H. Mirhosseini, T. M. A. Müller, M. Krack, T. D. Kühne, J. Hutter, The UZH protocol: Separating errors and constructing improved CP2K basis sets and pseudopotentials, *J. Chem. Phys.* 165, 104103 (2026). DOI: [10.1063/5.0347392](https://doi.org/10.1063/5.0347392)
- M. Iannuzzi *et al.*, The CP2K Program Package Made Simple, *J. Phys. Chem. B* 130, 1237 (2026). DOI: [10.1021/acs.jpcb.5c05851](https://doi.org/10.1021/acs.jpcb.5c05851)

---

[← 上一章](02-input-structure.md) | [返回目录](../README.md) | [下一章 →](04-cutoff-convergence.md)
