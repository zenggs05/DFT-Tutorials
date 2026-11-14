# 第11章 NMR 化学位移与其他分子性质（NMR Chemical Shifts, J-Couplings, Polarizability & Molecular Properties with Gaussian）

本章介绍如何用 Gaussian 16 计算 NMR 化学位移（GIAO 屏蔽常数、TMS 参比、线性回归缩放、等价核与多构象 Boltzmann 平均）和自旋-自旋耦合常数（`NMR=(SpinSpin,Mixed)`），以及偶极矩、多极矩、极化率、超极化率（`Polar`、`CPHF=RdFreq`、`Polar=Gamma`）、NICS 芳香性、旋光度、VCD/ROA 等常用分子性质。所有示例都给出完整的 Gaussian 输入、输出提取命令和结果检查方法，可直接用于结构鉴定、立体化学归属和非线性光学（NLO）材料筛选。

## 1. 本章目标与适用场景

学完本章你应该能够：

- 选择合适的几何结构水平与 NMR 计算水平，写出 GIAO 输入并读懂 `Isotropic` 屏蔽常数；
- 用同一水平计算的 TMS 或线性回归参数把屏蔽常数换算成化学位移，正确处理甲基等价氢、非对映 CH2 氢和多构象平均；
- 计算 H-H、C-H 自旋-自旋耦合常数，理解 Fermi 接触项、同位素和计算成本；
- 计算静态与频率相关的极化率 α、第一超极化率 β，以及第二超极化率 γ，并正确换算单位、注意约定差异；
- 了解 NICS、磁化率、旋光度（OptRot）、VCD、ROA 的 Gaussian 写法。

| 任务 | 本章小节 | 配套示例 |
|---|---|---|
| 1H/13C 化学位移 | 第 5、6 节 | [11-nmr-shift-camphor](../examples/11-nmr-shift-camphor/README.md) |
| J 耦合常数 | 第 7 节 | [11-nmr-jcoupling-ethanol](../examples/11-nmr-jcoupling-ethanol/README.md) |
| 极化率、超极化率 | 第 9 节 | [11-polar-hyperpol-pna](../examples/11-polar-hyperpol-pna/README.md) |

前置知识：几何优化与频率（[第5章](05-opt.md)、[第6章](06-freq-thermo.md)）、溶剂模型（[第9章](09-solvation.md)）。手性光谱 ECD 见 [第10章](10-tddft.md)。

## 2. 原理速览

### 2.1 屏蔽张量与化学位移

外磁场 $\mathbf B_0$ 在核处感应出局部场 $\mathbf B_{\rm loc}=(\mathbf 1-\boldsymbol\sigma)\mathbf B_0$，$\boldsymbol\sigma$ 为 3×3 的屏蔽张量。溶液 NMR 中分子快速翻滚，只观察到各向同性部分：

$$\sigma_{\rm iso}=\tfrac13(\sigma_{11}+\sigma_{22}+\sigma_{33})$$

Gaussian 输出中的 `Isotropic` 就是 $\sigma_{\rm iso}$（ppm），`Anisotropy` 描述张量的各向异性，主要用于固体 NMR。化学位移相对参比物（1H、13C、29Si 通常为 TMS）定义：

$$\delta=\frac{\sigma_{\rm ref}-\sigma}{1-\sigma_{\rm ref}}\approx\sigma_{\rm ref}-\sigma$$

因为 $\sigma_{\rm ref}\ll 10^6$ ppm，分母可以忽略。**$\sigma_{\rm ref}$ 必须在与样品完全相同的方法、基组、溶剂模型和几何水平下计算**，否则系统误差可达数 ppm。

### 2.2 规范原点问题：GIAO、CSGT、IGAIM

有限基组下，磁性质会依赖矢势的规范原点选择。Gaussian 提供的处理方法（官方 [NMR 关键词页](https://gaussian.com/nmr/)）：

- **GIAO**（默认）：每个基函数带有与磁场相关的相因子，结果与原点无关，收敛最快，是事实标准；
- **CSGT**：连续规范变换，需要较大基组才能达到相同精度；
- **IGAIM**：以原子中心作为规范原点，是 CSGT 的变体；
- **SingleOrigin**：单一原点，仅供比较，不推荐。

### 2.3 自旋-自旋耦合

间接自旋-自旋耦合（J 耦合）由四项组成：Fermi 接触项（FC）、自旋-偶极项（SD）、顺磁自旋-轨道项（PSO）和抗磁自旋-轨道项（DSO）。FC 项通常占主导，它取决于核位置处的自旋密度，需要描述核附近电子的紧致（tight）s 型基函数。约化耦合常数 $K$ 与同位素无关，$J$ 与核的磁旋比有关：

$$J_{AB}=\frac{h}{4\pi^2}\,\gamma_A\gamma_B\,K_{AB}$$

### 2.4 电学性质

分子在均匀电场 $\mathbf F$ 中的能量按 Taylor 级数展开（爱因斯坦求和）：

$$E(\mathbf F)=E_0-\mu_iF_i-\tfrac12\alpha_{ij}F_iF_j-\tfrac16\beta_{ijk}F_iF_jF_k-\tfrac1{24}\gamma_{ijkl}F_iF_jF_kF_l-\cdots$$

μ 为偶极矩，α 为极化率，β、γ 为第一、第二超极化率。常用原子单位换算（近似值）：

| 量 | 1 a.u. 等于 |
|---|---|
| 偶极矩 μ | 2.5417 Debye |
| 极化率 α | 0.14818 Å³（= $0.14818\times10^{-24}$ esu） |
| 第一超极化率 β | $8.6392\times10^{-33}$ esu |
| 第二超极化率 γ | $5.0367\times10^{-40}$ esu |

> [!WARNING]
> β 存在多种约定（Taylor 级数约定 T、微扰级数约定 B 等，相差系数 2 或 1/2），不同文献报道的 β 可能相差一倍。官方输出注明 “Values do not include the 1/n! factor”，与文献比较前务必确认约定。

## 3. NMR 关键词速查

`NMR` 是性质关键词，可用于 HF、所有 DFT 方法和 MP2，可与 `SCRF` 联用；HF 与 DFT 下 `NMR` 和 `Freq` 可以写在同一个 route 中。多个选项要加括号：`NMR=(SpinSpin,Mixed)`。

| 选项 | 说明 |
|---|---|
| `GIAO` | 默认，推荐 |
| `CSGT` / `IGAIM` / `SingleOrigin` / `All` | 其他规范处理；`All` 同时计算 SingleOrigin、IGAIM 和 CSGT |
| `SpinSpin` | 同时计算自旋-自旋耦合常数，成本约为频率计算的 2 倍，仅限 HF 与 DFT |
| `Mixed` | 两步法：第一步用解收缩并加紧致极化函数的基组算 FC 项，第二步用原基组算其余三项；对 6-311+G(d,p)、aug-cc-pVDZ 等价层基组显著提高 J 的精度 |
| `ReadAtoms` | 只对指定原子计算耦合，原子列表在单独的输入段中给出 |
| `FCOnly` / `ReadFC` | 只算 FC 项 / 从 chk 读取 FC 项再算其余项 |
| `PrintEigenvectors` | 打印屏蔽张量的本征向量 |
| `Susceptibility` | 同时计算磁化率 |

与 NMR 相关的另一个选项：`Pop=NCS`（基于 NBO 的自然化学屏蔽分析，把 GIAO 屏蔽分解为成键与孤对电子的贡献），用法见 [Population 关键词页](https://gaussian.com/population/) 与 [第12章](12-population-charges.md)。

## 4. 推荐计算方案

| 环节 | 常用选择 | 说明 |
|---|---|---|
| 几何结构 | B3LYP-D3(BJ)/6-31G(d) 或 def2-SVP，必要时 def2-TZVP，含溶剂模型 | 必须是频率确认过的极小点；结构误差直接传递到位移 |
| NMR 单点 | mPW1PW91/6-311+G(2d,p)、B3LYP/6-311+G(2d,p)、ωB97X-D/def2-TZVP | mPW1PW91 是 Tantillo 组 Cheshire 数据库常用的组合之一 |
| 溶剂 | `SCRF=(SMD 或 PCM,Solvent=...)`，与实验溶剂一致 | CDCl3 → `Chloroform`，DMSO-d6 → `DiMethylSulfoxide` |
| 参比 | 同水平 TMS，或同水平的线性回归参数 | 二选一，不要混用 |
| 构象 | CREST 搜索 + DFT 优化 + Boltzmann 平均 | 柔性分子必须做 |
| 专用基组 | pcSseg-n（从 [Basis Set Exchange](https://www.basissetexchange.org/) 下载，用 `Gen` 读入） | 针对屏蔽常数优化的基组，见 [第3章](03-method-basis.md) |

> [!TIP]
> [Cheshire NMR 数据库](http://cheshirenmr.info/) 给出了大量“几何水平 // NMR 水平 // 溶剂模型”组合的线性回归参数和测试误差。选方法时可以直接从中挑一个组合，严格按照其几何与 NMR 水平计算，再使用对应的 slope/intercept。

## 5. 实战一：TMS 参比

示例目录：[11-nmr-shift-camphor](../examples/11-nmr-shift-camphor/README.md)。

```text
! Chapter 11 example: TMS reference shieldings (1H, 13C, 29Si)
! The NMR step must use EXACTLY the same method/basis/solvent as the sample molecule
%chk=tms_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt Freq SCRF=(SMD,Solvent=Chloroform)

TMS opt+freq, B3LYP-D3(BJ)/6-31G(d), SMD(chloroform)

0 1
 C       1.461021     0.837646    -0.824830
 Si     -0.000000     0.000000     0.000000
 C      -0.957190     1.260495     1.005743
 C       0.617294    -1.363057     1.130311
 C      -1.121126    -0.735085    -1.311224
 H       2.034298     0.118915    -1.417873
 H       2.129200     1.274173    -0.076563
 H       1.122815     1.637710    -1.489988
 H      -1.815534     0.792097     1.496426
 H      -0.320707     1.702564     1.778003
 H      -1.327093     2.066102     0.364579
 H       1.274415    -0.955377     1.904204
 H       1.179510    -2.110635     0.562893
 H      -0.220412    -1.865843     1.622628
 H      -0.581696    -1.474432    -1.910642
 H      -1.981619    -1.229639    -0.850910
 H      -1.493176     0.044365    -1.982758

--Link1--
%oldchk=tms_opt.chk
%chk=tms_nmr.chk
%mem=8GB
%nprocshared=8
#p mPW1PW91/6-311+G(2d,p) NMR=GIAO SCRF=(SMD,Solvent=Chloroform)
   Geom=Check Guess=Read

TMS GIAO shieldings, mPW1PW91/6-311+G(2d,p), SMD(chloroform)

0 1

```

| 关键词 | 说明 |
|---|---|
| `B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt Freq` | 几何水平：B3LYP-D3(BJ)；Freq 确认无虚频 |
| `SCRF=(SMD,Solvent=Chloroform)` | 几何和 NMR 两步都用氯仿 |
| `mPW1PW91/6-311+G(2d,p) NMR=GIAO` | NMR 水平；`NMR` 不带选项时也是 GIAO |
| `Geom=Check Guess=Read` | 读取优化结构与轨道 |

NMR 输出格式（官方示例，节选）：

```text
Magnetic properties (GIAO method)

Magnetic shielding (ppm):
  1  C    Isotropic =    57.7345   Anisotropy =   194.4092
   XX=    48.4143   YX=      .0000   ZX=      .0000
   XY=      .0000   YY=   -62.5514   ZY=      .0000
   XZ=      .0000   YZ=      .0000   ZZ=   187.3406
```

提取并平均 TMS 中 12 个 H 和 4 个 C 的屏蔽常数：

```bash
grep "Isotropic" tms_reference_nmr.log | \
  awk '$2=="H"{h+=$5;nh++} $2=="C"{c+=$5;nc++} $2=="Si"{si=$5}
       END{printf "sigma_ref(1H)=%.4f (n=%d)  sigma_ref(13C)=%.4f (n=%d)  sigma(29Si)=%.4f\n", h/nh, nh, c/nc, nc, si}'
```

对 B3LYP 一类泛函配合三重 ζ 基组，TMS 的 $\sigma(^1{\rm H})$ 通常在 31–32 ppm 附近，$\sigma(^{13}{\rm C})$ 在 180–190 ppm 附近（近似范围，随方法变化，请以自己的计算为准）。如果你的值与此相差很远，先检查是否算错了分子或水平。

## 6. 实战二：樟脑的 1H/13C 化学位移

### 6.1 完整输入

樟脑（(1R,4R)-camphor，C10H16O）是刚性的双环萜类，常用作天然产物 NMR 计算的练习分子。输入与 TMS 采用完全相同的协议：

```text
! Chapter 11 example: 1H/13C shieldings of (1R,4R)-camphor (rigid natural-product-like molecule)
! Same protocol as tms_reference_nmr.gjf; convert with delta = sigma(TMS) - sigma
%chk=camphor_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt Freq SCRF=(SMD,Solvent=Chloroform)

Camphor opt+freq, B3LYP-D3(BJ)/6-31G(d), SMD(chloroform)

0 1
 C      -1.974491    -0.729766     0.060723
 C      -0.627665    -0.170949    -0.442057
 C      -0.478487    -0.635059    -1.904031
 C      -0.458552     1.362153    -0.250553
 C       1.016507     1.603848    -0.660990
 C       1.739555     0.341128    -0.141135
 C       0.587399    -0.486952     0.473105
 C       0.912644    -1.936917     0.739336
 C       0.199552     0.320333     1.716218
 O       0.461775     0.030355     2.873072
 C      -0.546773     1.566292     1.282145
 H      -2.031474    -1.811985    -0.103998
 H      -2.806299    -0.265872    -0.481328
 H      -2.152188    -0.572112     1.127935
 H      -1.283647    -0.217522    -2.519493
 H       0.464303    -0.347251    -2.374414
 H      -0.539721    -1.727233    -1.972101
 H      -1.164608     1.985532    -0.804537
 H       1.120057     1.690410    -1.747682
 H       1.429251     2.517179    -0.219593
 H       2.502595     0.606300     0.599957
 H       2.239139    -0.186092    -0.960877
 H       1.784208    -2.027092     1.397056
 H       0.081188    -2.448233     1.235619
 H       1.136985    -2.472405    -0.188449
 H      -0.035295     2.470573     1.619533
 H      -1.575962     1.541336     1.646541

--Link1--
%oldchk=camphor_opt.chk
%chk=camphor_nmr.chk
%mem=8GB
%nprocshared=8
#p mPW1PW91/6-311+G(2d,p) NMR=GIAO SCRF=(SMD,Solvent=Chloroform)
   Geom=Check Guess=Read

Camphor GIAO shieldings, mPW1PW91/6-311+G(2d,p), SMD(chloroform)

0 1

```

Gaussian 原子编号与樟脑常规编号的对应（由连接关系确定）：

| Gaussian 编号 | 樟脑编号 | 说明 |
|---|---|---|
| 7 | C1 | 季碳桥头，连 C10 甲基 |
| 9 / 10 | C2 / O | 羰基 |
| 11（H26, H27） | C3 | 与羰基相邻的 CH2，两个 H 为 exo/endo，不等价 |
| 4（H18） | C4 | 桥头 CH |
| 5（H19, H20） | C5 | CH2，exo/endo 不等价 |
| 6（H21, H22） | C6 | CH2，exo/endo 不等价 |
| 2 | C7 | 桥碳，带两个甲基 |
| 1（H12–H14）、3（H15–H17） | C8、C9 | C7 上的两个甲基，彼此**非对映**，不等价 |
| 8（H23–H25） | C10 | C1 上的甲基 |

### 6.2 屏蔽常数换算为化学位移

```bash
SH=31.80    # 用你自己 TMS 计算得到的 sigma_ref(1H) 替换（此处为占位值）
SC=184.00   # 用你自己 TMS 计算得到的 sigma_ref(13C) 替换（此处为占位值）
grep "Isotropic" camphor_nmr_giao_smd.log | \
  awk -v rh=$SH -v rc=$SC '$2=="H"{printf "%3d H  sigma=%8.3f  delta=%7.2f\n",$1,$5,rh-$5}
                           $2=="C"{printf "%3d C  sigma=%8.3f  delta=%7.2f\n",$1,$5,rc-$5}'
```

仓库脚本 [`nmr_shift.py`](../scripts/nmr_shift.py) 可以完成 TMS 参比或线性回归换算、等价核平均、多构象 Boltzmann 平均，并与实验值比较给出 MAE，运行 `python ../scripts/nmr_shift.py -h` 查看用法。GaussView 6 也能在 Results → NMR 中选择内置的 TMS 参比并显示谱图，但内置参比的计算水平必须与你的作业一致才有意义。

### 6.3 等价核平均

- **快速旋转的甲基**：三个 H 在 NMR 时间尺度上平均，计算值可能相差 1 ppm 以上，必须取算术平均（本例 H12–14、H15–17、H23–25 各自平均）；
- **对称等价核**：如对位取代苯环的 2,6 位 H，取平均；
- **非对映 H（diastereotopic）**：同一 CH2 上的两个 H 处在不同化学环境，实验上通常分开出峰，**不能**平均（本例 C3、C5、C6 上的 exo/endo H）；
- **可交换质子**（OH、NH、COOH）：受氢键、浓度和交换影响很大，一般不纳入统计。

### 6.4 线性回归缩放

TMS 参比法依赖单个分子的误差抵消。另一种做法是对一组已知分子拟合 $\sigma_{\rm calc}=m\,\delta_{\rm exp}+b$，然后

$$\delta_{\rm scaled}=\frac{b-\sigma_{\rm calc}}{-m}$$

例如 Cheshire 数据库 Table #1a 中，气相 B3LYP/6-31+G(d,p)//B3LYP/6-31G(d) 的 1H 参数为 slope = −1.0381、intercept = 31.7427（引自 Cheshire 网站）。这些参数**只适用于完全相同的协议**（几何水平、NMR 水平、溶剂模型、Gaussian 版本差异也可能有影响）；换了任何一项都应查对应的条目，或用 20–30 个可靠实验数据的分子自己拟合。

### 6.5 多构象 Boltzmann 平均

对柔性分子，每个构象 k 分别计算位移 $\delta_{i,k}$，再按自由能加权：

$$\bar\delta_i=\sum_k w_k\,\delta_{i,k},\qquad w_k=\frac{e^{-\Delta G_k/RT}}{\sum_j e^{-\Delta G_j/RT}}$$

- 构象搜索用 [CREST](https://github.com/crest-lab/crest)，DFT 重新优化后去重，保留约 3 kcal/mol 窗口内的构象；
- ΔG 最好用高水平单点能加上低水平热校正（[第6章](06-freq-thermo.md)），[GoodVibes](https://github.com/patonlab/GoodVibes) 可直接输出 Boltzmann 权重；
- 区分非对映异构体时可以使用 DP4+ 等统计方法（Grimblat, Zanardi, Sarotti, *J. Org. Chem.* 2015），它综合比较 1H 与 13C 的计算-实验偏差，给出各候选结构的概率。

## 7. 实战三：乙醇的自旋-自旋耦合常数

### 7.1 完整输入

示例目录：[11-nmr-jcoupling-ethanol](../examples/11-nmr-jcoupling-ethanol/README.md)。

```text
! Chapter 11 example: indirect spin-spin coupling constants of ethanol (anti conformer)
! NMR=(SpinSpin,Mixed): Fermi-contact term with a modified (uncontracted + tight) basis,
! the other three terms with the route basis. Default isotopes: use the H-H couplings;
! see README for couplings involving 13C.
%chk=etoh_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-311+G(d,p) EmpiricalDispersion=GD3BJ Opt Freq

Ethanol (anti) opt+freq, B3LYP-D3(BJ)/6-311+G(d,p)

0 1
 C      -1.012894     0.103762     0.116897
 C       0.426069    -0.316378    -0.093435
 O       1.261496     0.826907     0.009591
 H      -1.309983     0.849384    -0.628010
 H      -1.139273     0.568910     1.100139
 H      -1.685279    -0.755560     0.043851
 H       0.728655    -1.043556     0.665981
 H       0.555645    -0.759303    -1.085447
 H       2.175565     0.525833    -0.129567

--Link1--
%oldchk=etoh_opt.chk
%chk=etoh_nmr_jj.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-311+G(d,p) NMR=(SpinSpin,Mixed) Geom=Check Guess=Read

Ethanol shieldings + spin-spin couplings, B3LYP/6-311+G(d,p), Mixed FC scheme

0 1

```

`Mixed` 会自动运行两个子步骤，最终结果在第二步末尾给出，并已包含第一步的 FC 项。官方示例输出格式（节选，HF 分子）：

```text
 Total nuclear spin-spin coupling K (Hz):
                1             2
      1  0.000000D+00
      2  0.147308D+02  0.000000D+00
 Total nuclear spin-spin coupling J (Hz):
                1             2
      1  0.000000D+00
      2  0.432614D+03  0.000000D+00
```

K、J 都以下三角矩阵打印，原子较多时会按每 5 列一块分段输出。提取最终 J 矩阵：

```bash
grep -n "Total nuclear spin-spin coupling J (Hz)" ethanol_nmr_spinspin_mixed.log
# 取最后一次出现的位置，向后打印足够多的行（9 个原子约 15 行）
awk '/Total nuclear spin-spin coupling J \(Hz\)/{n=NR} {line[NR]=$0} END{for(i=n;i<=n+15;i++) print line[i]}' \
    ethanol_nmr_spinspin_mixed.log
```

注意 `D+02` 是 Fortran 的科学计数法（$10^{2}$），用 Python 解析时先把 `D` 替换成 `E`。

### 7.2 结果解读

- 乙醇中 H4/H5/H6（甲基）与 H7/H8（亚甲基）之间的 $^3J_{\rm HH}$ 取决于二面角（Karplus 关系）。由于甲基快速旋转，实验值是平均值：应对 3×2 = 6 个 H-H 耦合取平均后再与实验比较；
- OH 质子的耦合在多数溶剂中因快速交换而观察不到；
- 乙醇存在 anti 与 gauche 构象（C-C-O-H 二面角），严格比较时需对构象做 Boltzmann 平均；
- 自旋-自旋耦合计算的成本约为频率计算的 2 倍（官方说明），大分子可以用 `ReadAtoms` 只算关心的原子。

### 7.3 同位素与 13C 耦合

`J` 矩阵按作业中实际使用的同位素计算（官方说明：显式指定的或默认的同位素），而 `K` 与同位素无关。碳的默认同位素是天然丰度最高的 12C（核自旋为 0），所以想得到 $^1J_{\rm CH}$ 等 13C 耦合时，应在分子说明中用 `Iso=13` 显式指定（[分子说明文档](https://gaussian.com/molspec/)）：

```text
#p B3LYP/6-311+G(d,p) NMR=(SpinSpin,Mixed)

Ethanol J(C,H) with explicit 13C labels (paste optimized coordinates here)

0 1
 C(Iso=13)   x1  y1  z1
 C(Iso=13)   x2  y2  z2
 O           x3  y3  z3
 ...
```

`Geom=Check` 无法在读取结构时再修改同位素，因此先用仓库脚本 `log2xyz.py`（`python ../scripts/log2xyz.py -h`）从优化 log 中提取最终坐标，再写成上面的独立输入。也可以用 $J=\frac{h}{4\pi^2}\gamma_A\gamma_B K$ 从 K 手工换算。

### 7.4 ReadAtoms：只算部分原子

官方给出的原子列表语法（单独输入段，空行结束）：

```text
atoms=3-6,17 notatoms=5
```

含义：加入 3、4、6、17 号原子，排除 5 号原子。也可以按元素：`atoms=C N notatoms=5`；裸整数直接表示原子编号，如 `1,3,5 7`。route 写成 `NMR=(SpinSpin,ReadAtoms)`，该输入段放在分子说明之后。

## 8. 其他磁性质

### 8.1 NICS 芳香性

核独立化学位移（NICS）是在环中心（NICS(0)）或环平面上方 1 Å（NICS(1)）放置虚原子 Bq，计算其屏蔽：$\mathrm{NICS}=-\sigma_{\rm iso}({\rm Bq})$。显著为负表示芳香性（抗磁环流），为正表示反芳香性。下面是苯的示例输入（理想 D6h 结构，C-C 1.39 Å、C-H 1.09 Å，仅用于演示；正式工作应先优化）：

```text
%chk=benzene_nics.chk
%mem=4GB
%nprocshared=4
#p B3LYP/6-311+G(d,p) NMR=GIAO

Benzene NICS(0) and NICS(1) with ghost atoms

0 1
 C    1.390000    0.000000    0.000000
 C    0.695000    1.203775    0.000000
 C   -0.695000    1.203775    0.000000
 C   -1.390000    0.000000    0.000000
 C   -0.695000   -1.203775    0.000000
 C    0.695000   -1.203775    0.000000
 H    2.480000    0.000000    0.000000
 H    1.240000    2.147743    0.000000
 H   -1.240000    2.147743    0.000000
 H   -2.480000    0.000000    0.000000
 H   -1.240000   -2.147743    0.000000
 H    1.240000   -2.147743    0.000000
 Bq   0.000000    0.000000    0.000000
 Bq   0.000000    0.000000    1.000000

```

`Bq` 原子没有核电荷、电子和基函数，只作为计算屏蔽的探针。输出中 `13  Bq   Isotropic = ...` 一行的数值取负即为 NICS(0)。NICS(1)zz（只取垂直环平面的 zz 分量）比各向同性值更能反映 π 环流，读取 `ZZ=` 分量即可（前提是环平面为 xy 平面，Gaussian 标准取向下苯正好如此）。

### 8.2 磁化率与 GaussView 显示

- `NMR=Susceptibility` 或 `Polar=Susceptibility`（仅 HF、DFT）同时计算磁化率；
- GaussView 打开 NMR 作业的 log，Results → NMR 显示谱线，可选择元素与参比。

## 9. 实战四：p-硝基苯胺的极化率与超极化率

### 9.1 完整输入

示例目录：[11-polar-hyperpol-pna](../examples/11-polar-hyperpol-pna/README.md)。p-硝基苯胺（pNA）是推-拉型 NLO 分子的标准测试体系。

```text
! Chapter 11 example: polarizability and first hyperpolarizability of p-nitroaniline (pNA)
! Step 1: opt+freq   Step 2: static + dynamic (1064 nm) alpha and beta via CPHF=RdFreq
%chk=pna_opt.chk
%mem=8GB
%nprocshared=8
#p B3LYP/6-311+G(d,p) EmpiricalDispersion=GD3BJ Opt Freq

p-Nitroaniline opt+freq, B3LYP-D3(BJ)/6-311+G(d,p), gas phase

0 1
 N      -2.814048    -0.077446    -0.266096
 C      -1.432583    -0.013928    -0.064026
 C      -0.728316    -1.154152     0.336857
 C       0.667370    -1.141629     0.424927
 C       1.368118     0.015152     0.066332
 C       0.681616     1.145223    -0.391314
 C      -0.714183     1.114575    -0.472914
 N       2.833187     0.036722     0.152334
 O       3.415545     1.078876    -0.181038
 O       3.402674    -0.987357     0.556452
 H      -3.290597    -0.763289     0.311617
 H      -3.280706     0.824408    -0.255079
 H      -1.264427    -2.067415     0.584886
 H       1.185050    -2.038535     0.757846
 H       1.210486     2.044443    -0.699481
 H      -1.239186     1.984351    -0.861302

--Link1--
%oldchk=pna_opt.chk
%chk=pna_polar.chk
%mem=8GB
%nprocshared=8
#p CAM-B3LYP/Aug-cc-pVDZ Polar CPHF=RdFreq Geom=Check Guess=Read

p-Nitroaniline alpha/beta, static and 1064 nm, CAM-B3LYP/aug-cc-pVDZ

0 1

1064nm

```

| 关键词 | 说明 |
|---|---|
| `Polar` | 计算偶极极化率；HF 与大多数 DFT 有解析三阶导数，默认同时给出超极化率（官方 Polar 页） |
| `CPHF=RdFreq` | 频率相关（动态）CPHF，从输入读入入射光频率；默认单位为 Hartree，可加后缀 `nm`（波长）或 `cm`（cm⁻¹） |
| `1064nm` | Nd:YAG 激光波长；放在电荷多重度行后的空行之后（`Geom=Check` 时分子说明为空） |
| `CAM-B3LYP/Aug-cc-pVDZ` | 超极化率对弥散函数非常敏感；常规杂化泛函往往高估推-拉体系的 β，范围分离泛函更稳妥 |

按 CPHF 的 `Static` 默认行为，动态计算时会自动包含静态（ω=0）结果。

### 9.2 读取输出

（示意输出，只展示格式，`...` 处为具体数值）

```text
 Dipole polarizability, Alpha (input orientation).
 (esu units = cm**3 , SI units = C**2 m**2 J**-1)
 Alpha(0;0):
               (au)            (10**-24 esu)      (10**-40 SI)
   iso         ...
   aniso       ...
   xx          ...
 First dipole hyperpolarizability, Beta (input orientation).
 ||, _|_  parallel and perpendicular components, (z) with respect to z axis,
 vector components x,y,z.  Values do not include the 1/n! factor of 1/2.
 Beta(0;0,0):
               (au)            (10**-30 esu)      (10**-50 SI)
   || (z)      ...
   x           ...
   y           ...
   z           ...
```

```bash
grep -n "Alpha(" pna_polar_hyperpolarizability.log          # 列出静态与动态极化率块
grep -n "Beta("  pna_polar_hyperpolarizability.log          # 列出静态与动态超极化率块
grep -A12 "Beta(0;0,0)" pna_polar_hyperpolarizability.log | head -14
grep "Isotropic polarizability" pna_polar_hyperpolarizability.log
```

动态结果的标签中会出现频率 w（如 `Alpha(-w;w)`、`Beta(-w;w,0)` 等），具体打印哪些 β 分量（电光效应 EOPE、二次谐波 SHG 等）以你的输出标签为准。官方 Polar 页说明：经典格式中极化率以下三角顺序（xx, xy, yy, xz, yz, zz）、超极化率以下四面体顺序（xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz）打印，取向为标准取向；G16 的汇总块会在标题中注明所用取向（如 input orientation），读取分量前先看清标题，不同取向下的分量不能混用。

常用标量：

$$\bar\alpha=\tfrac13(\alpha_{xx}+\alpha_{yy}+\alpha_{zz}),\qquad
\beta_i=\beta_{iii}+\tfrac13\sum_{j\ne i}(\beta_{ijj}+\beta_{jij}+\beta_{jji}),\qquad
\beta_{\rm tot}=\sqrt{\beta_x^2+\beta_y^2+\beta_z^2}$$

静态情况下 β 张量完全对称，$\beta_i=\beta_{iii}+\sum_{j\ne i}\beta_{ijj}$。$\beta_\parallel$（`||`）是 β 在偶极方向上的投影，常用于与 EFISH 实验比较。Multiwfn 主功能 24 提供（超）极化率分析，包括单位换算和不同约定下的量。

### 9.3 第二超极化率、旋光度、VCD 与 ROA

第二超极化率 γ（`Polar=Gamma` 等价于 `Polar=(DCSHG,Cubic)`，隐含 `CPHF=RdFreq`，需要给出频率）：

```text
#p CAM-B3LYP/Aug-cc-pVDZ Polar=Gamma Geom=Check Guess=Read

pNA second hyperpolarizability

0 1

1064nm

```

旋光度（`Polar=OptRot`，仅 HF/DFT，不能与 NMR 同时使用）。钠 D 线 589.3 nm 时输出 `[Alpha] ( 5893.0 A) = ... deg.`：

```text
#p CAM-B3LYP/Aug-cc-pVDZ Polar=OptRot CPHF=RdFreq Geom=Check Guess=Read

Specific rotation at the sodium D line

0 1

589.3nm

```

VCD 强度：在 HF/DFT 频率计算中加 `Freq=VCD`，例如 `#p B3LYP/6-311+G(d,p) Opt Freq=VCD`；比较前同样需要构象搜索和 Boltzmann 平均，可与 ECD（[第10章](10-tddft.md)）互相印证绝对构型。

ROA：官方推荐两步法——先用较小基组做 `Opt Freq`，再用更大基组做 `Polar=ROA`（从 chk 读取力常数），并在输入末尾给出入射光波长（如 `532nm`）。官方 [Polar 页](https://gaussian.com/polar/) 给出的推荐组合是 APFD/6-311+G(2d,2p) 做 Opt Freq，APFD/spAug-cc-pVTZ 做 ROA。

### 9.4 偶极矩与多极矩

每个作业默认都会打印偶极矩、四极矩、八极矩和十六极矩：

```bash
grep -A1 "Dipole moment (field-independent basis, Debye)" job.log | tail -2
grep -A2 "Traceless Quadrupole moment" job.log | tail -3
```

注意：偶极矩在标准取向下给出，Gaussian 的符号约定为从负电中心指向正电中心（与部分化学教材的箭头方向相反）；带电体系的偶极矩依赖坐标原点，没有绝对意义。实验偶极矩多为溶液测量，比较时要考虑溶剂和构象平均。

## 参数选择与收敛建议

| 项目 | 推荐 | 说明 |
|---|---|---|
| NMR 规范方法 | GIAO | CSGT 需要更大基组 |
| NMR 泛函/基组 | mPW1PW91、B3LYP、ωB97X-D + 6-311+G(2d,p) 或 def2-TZVP、pcSseg-2 | 与参比/回归参数严格一致 |
| 几何结构 | 含色散的 DFT + 溶剂，频率确认 | 结构误差是位移误差的主要来源之一 |
| 参比 | 同水平 TMS 或同协议回归参数 | 不要用别人的 TMS 数值 |
| 等价核 | 甲基 H、对称等价核取平均 | 非对映 H 不平均 |
| 柔性分子 | 构象搜索 + ΔG 加权 | DP4+ 用于异构体判别 |
| J 耦合 | `NMR=(SpinSpin,Mixed)`，6-311+G(d,p) 或 aug-cc-pVDZ 以上 | 13C 耦合要 `Iso=13` |
| 极化率 | 至少含一组弥散函数（aug-cc-pVDZ、6-311++G(d,p)） | 无弥散函数时 α 明显偏小 |
| 超极化率 | 范围分离泛函 + 弥散基组；注明约定与频率 | β 对方法、基组极敏感 |
| 积分格点 | G16 默认 UltraFine，CPHF 默认 SG1 | 与 G09 结果比较时注意默认值不同 |

## 常见问题与报错

**1. 计算位移整体偏差数 ppm**
原因：TMS 参比与样品不是同一水平（基组、泛函、溶剂或几何水平不同），或者拿了文献中的 TMS 数值。
解决：在完全相同的协议下重算 TMS，或改用同协议的回归参数。

**2. 同一个甲基的三个 H 给出三个不同的位移**
原因：静态结构中三个 H 的环境不同，而实验观察到的是快速旋转后的平均值。
解决：对甲基 H 取平均；但不要把非对映的 CH2 氢也平均掉。

**3. 所有含碳的 J 都是 0**
原因：碳的默认同位素是 12C（自旋为 0），J 按实际同位素计算。
解决：在坐标中写 `C(Iso=13)`，或从同位素无关的 K 矩阵换算（第 7.3 节）。

**4. `QPErr --- A syntax error was detected in the input line.`**
原因：多个选项没有加括号，例如 `NMR=SpinSpin,Mixed`；或泛函名拼错（如把 `mPW1PW91` 写成 `mPW1PW92`）。
解决：写成 `NMR=(SpinSpin,Mixed)`；对照 [关键词列表](https://gaussian.com/keywords/)。

**5. 在未优化或有虚频的结构上算 NMR，与实验偏差很大**
原因：NMR 对键长、二面角敏感；构象错误时偏差可达数 ppm。
解决：先 Opt Freq 确认极小点；柔性分子做构象搜索。

**6. 作业在读取频率处异常终止，或得到的只有静态结果**
原因：`CPHF=RdFreq`/`Polar=OptRot`/`Polar=ROA` 需要的频率行缺失，或空行位置不对。
解决：在电荷多重度行（`Geom=Check` 时其后为空分子说明）后留一个空行，再写 `1064nm`，最后以空行结束；单位后缀写成 `nm` 或 `cm`。

**7. `Polar=OptRot` 与 `NMR` 写在同一 route 中报错**
原因：官方说明 OptRot 不能与 NMR 组合。
解决：拆成两个 Link1 步骤。

**8. 对 CCSD(T) 等方法同时写 `Freq Polar` 失败**
原因：官方说明缺少解析梯度的方法（MP4(SDTQ)、QCISD(T)、CCSD(T) 等）不能把 Freq 与 Polar 组合。
解决：分开计算，或改用有解析导数的方法。

**9. 自旋-自旋耦合作业耗时或内存远超预期**
原因：SpinSpin 的成本约为频率计算的 2 倍，`Mixed` 的第一步还使用了更大的修改基组。
解决：用 `ReadAtoms` 限定原子；增加 `%mem` 并确保与作业申请一致（[第17章](17-hpc-jobs.md)）。

更多通用报错见 [第18章](18-errors-capstone.md)。

## 练习

1. 用 B3LYP/6-311+G(2d,p) 替换 mPW1PW91 重算 TMS 与樟脑，比较两种泛函得到的 13C 位移 MAE（实验数据请查阅文献或 SDBS 等谱图数据库）。
2. 把樟脑的溶剂模型去掉（气相 NMR），观察 C2（羰基碳）位移的变化。
3. 对乙醇的 gauche 构象（C-C-O-H 约 60°）做同样的耦合计算，按 Boltzmann 权重平均后与 anti 构象比较 $^3J_{\rm HH}$。
4. 计算苯和环丁二烯（单重态矩形结构）的 NICS(1)zz，比较芳香与反芳香体系。
5. 对 pNA 分别用 B3LYP、CAM-B3LYP、ωB97X-D 计算静态 $\beta_{\rm tot}$，再把基组换成 6-31G(d)（无弥散函数），讨论方法与基组的影响。

## 参考资料

官方文档：

- [NMR 关键词](https://gaussian.com/nmr/)、[Polar 关键词](https://gaussian.com/polar/)、[CPHF 关键词（RdFreq）](https://gaussian.com/cphf/)
- [Freq 关键词（VCD、ROA）](https://gaussian.com/freq/)、[分子说明（Iso、Spin、Bq）](https://gaussian.com/molspec/)
- [SCRF 关键词与溶剂列表](https://gaussian.com/scrf/)、[Population 关键词（NCS）](https://gaussian.com/population/)
- [Cheshire NMR 数据库与线性回归参数](http://cheshirenmr.info/ScalingFactors.htm)
- [Basis Set Exchange（pcSseg-n 等 NMR 专用基组）](https://www.basissetexchange.org/)
- [Multiwfn 在线手册](http://sobereva.com/multiwfn/manual/)

相关 GitHub 仓库（仅链接）：

- [cclib/cclib](https://github.com/cclib/cclib)：解析 Gaussian 输出中的屏蔽常数、极化率等性质
- [crest-lab/crest](https://github.com/crest-lab/crest)：NMR 多构象平均前的构象搜索
- [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)：自由能校正与 Boltzmann 权重
- [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py)：输入生成、结果解析与构象处理
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：中文工作流项目，含 NMR 模板与自查项
- [lmmentel/awesome-python-chemistry](https://github.com/lmmentel/awesome-python-chemistry)：化学 Python 工具索引

关键文献：

- K. Wolinski, J. F. Hinton, P. Pulay, *J. Am. Chem. Soc.* 112, 8251 (1990). DOI: [10.1021/ja00179a005](https://doi.org/10.1021/ja00179a005)（GIAO）
- J. R. Cheeseman, G. W. Trucks, T. A. Keith, M. J. Frisch, *J. Chem. Phys.* 104, 5497 (1996). DOI: [10.1063/1.471789](https://doi.org/10.1063/1.471789)（DFT 下 GIAO 与 CSGT 的比较）
- M. W. Lodewyk, M. R. Siebert, D. J. Tantillo, *Chem. Rev.* 112, 1839 (2012). DOI: [10.1021/cr200106v](https://doi.org/10.1021/cr200106v)（NMR 计算与结构鉴定综述）
- N. Grimblat, M. M. Zanardi, A. M. Sarotti, *J. Org. Chem.* 80, 12526 (2015). DOI: [10.1021/acs.joc.5b02396](https://doi.org/10.1021/acs.joc.5b02396)（DP4+）
- W. Deng, J. R. Cheeseman, M. J. Frisch, *J. Chem. Theory Comput.* 2, 1028 (2006). DOI: [10.1021/ct600110u](https://doi.org/10.1021/ct600110u)（`Mixed` 自旋-自旋耦合方案）
- T. Helgaker, M. Jaszuński, K. Ruud, *Chem. Rev.* 99, 293 (1999). DOI: [10.1021/cr960017t](https://doi.org/10.1021/cr960017t)（NMR 参数的从头计算综述）
- P. v. R. Schleyer et al., *J. Am. Chem. Soc.* 118, 6317 (1996). DOI: [10.1021/ja960582d](https://doi.org/10.1021/ja960582d)（NICS）
- J. R. Cheeseman, M. J. Frisch, *J. Chem. Theory Comput.* 7, 3323 (2011). DOI: [10.1021/ct200507e](https://doi.org/10.1021/ct200507e)（Raman/ROA 的基组依赖性）

[← 上一章](10-tddft.md) | [返回目录](../README.md) | [下一章 →](12-population-charges.md)
