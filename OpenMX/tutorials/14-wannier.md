# 第14章 Wannier 函数：MLWF 与 CWF（Wannier Functions in OpenMX: MLWF and CWF）

本章介绍 OpenMX 中两种构造 Wannier 函数的方法：最大局域化 Wannier 函数（MLWF，Marzari–Vanderbilt–Souza 方法，3.5 起提供）和 4.0 新增的最接近 Wannier 函数（CWF，Closest Wannier Functions）。全章用三个算例串起来：Si 的 sp3 型与键中心 MLWF、石墨烯 pz 模型的解纠缠，以及 GaAs 的 AO/HO 引导 CWF。内容包括内外能量窗口怎么设、投影子怎么写、多项式窗口函数、Wannier 插值能带的检验、spread 与距离函数（DM_func）的判读，以及如何输出 Wannier90 格式文件交给 Wannier90、postw90、WannierTools 继续计算。关键词：OpenMX Wannier、MLWF、CWF、`Wannier.Func.Calc`、`CWF.Calc`、disentanglement、Wannier interpolation、tight-binding model、Wannier90 interface。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清楚 Wannier 函数的规范自由度，以及 MLWF（最小化 spread）与 CWF（最接近给定引导函数）在思路上的区别；
- 写出完整的 MLWF 输入：`Wannier.Func.Calc`、`Wannier.Func.Num`、外/内能量窗口、`<Wannier.Initial.Projectors`（含 `Wannier.Initial.Projectors.Unit FRAC|ANG|AU`）、`Wannier.Kgrid`、`Wannier.MaxShells`；
- 控制两步优化：解纠缠（`Wannier.Dis.*`）和规范相关部分的最小化（`Wannier.Minimizing.*`），并会用 `grep DISE/CONV/SPRD/CENT` 检查收敛；
- 写出 4.0 的 CWF 输入：`CWF.Calc`、`CWF.Guiding.Orbitals AO|HO|MO`、`<CWF.Guiding.AOs`、`CWF.window.function Poly|Fermi`、`CWF.disentangling.Erange`、`CWF.disentangling.smearing.kBT`、`CWF.fileout.cube`、`CWF.Wannier90.fileout`；
- 用插值能带与原始能带的对比、spread、DM_func、跳跃积分随距离的衰减来判断 Wannier 函数“好不好”；
- 生成 Wannier90 接口文件，在 Wannier90/postw90/WannierTools 中继续计算。

| 典型需求 | 推荐方法 | 说明 |
|---|---|---|
| 化学键分析、轨道图像（键轨道、杂化轨道、d 轨道） | CWF（AO/HO 引导） | 形状接近引导函数，物理意义直观；一次 SVD 完成，无迭代 |
| 少带紧束缚模型（如石墨烯 pz、过渡金属 $t_{2g}$） | MLWF 或 CWF | MLWF 可以只选 pz 这类单个 m 分量；CWF 的 AO 模式只能整壳层选择（见第 7 节） |
| 分子晶体费米面附近少数几条带 | CWF（MO 引导） | 以分子轨道为引导，见第 8.4 节 |
| 极化、Wannier 电荷中心、局域性定量比较 | MLWF | spread 有明确定义，便于与文献比较 |
| Berry 曲率、AHC、光电导、表面态 | 任一方法 + Wannier90/WannierTools | OpenMX 负责生成 `.amn/.mmn/.eig/.win` |
| 原子电荷（第四种布居分析） | CWF | `.out` 中给出 CWF 布居，见 [第10章](10-charge-mo-cube.md) |

配套示例：

- Si 的 8 个 sp3 型 MLWF、4 个键中心 MLWF、Wannier90 接口：[../examples/14-si-mlwf-sp3/README.md](../examples/14-si-mlwf-sp3/README.md)
- 石墨烯 pz 双带模型（解纠缠）：[../examples/14-graphene-pz-mlwf/README.md](../examples/14-graphene-pz-mlwf/README.md)
- GaAs 的 AO/HO 引导 CWF 与窗口扫描：[../examples/14-gaas-cwf/README.md](../examples/14-gaas-cwf/README.md)

> [!NOTE]
> 版本说明：MLWF 从 OpenMX 3.5 起提供（H. Weng 实现），Wannier90 接口（`Wannier90.fileout`）从 3.8 起提供，CWF 是 4.0（2026-04）的新功能，3.9 没有。4.0 手册第 47 章 [Generation of Wannier Functions](https://www.openmx-square.org/openmx_man4.0/s47_wannier.html) 分为 47.1（CWF）和 47.2（MLWF）两部分。本章的 MLWF 输入在 3.9 下也能用（删掉 `scf.Generation.Kpoint` 一行即可）；CWF 输入只能在 4.0 下运行。

---

## 2. 原理速览

### 2.1 Bloch 函数、Wannier 函数与规范自由度

对一组 $J$ 条能带的 Bloch 函数 $\psi_{n\mathbf k}$，Wannier 函数定义为

$$
w_{n\mathbf R}(\mathbf r)=\frac{V}{(2\pi)^3}\int_{\mathrm{BZ}} d\mathbf k\, e^{-i\mathbf k\cdot\mathbf R}\sum_{m=1}^{J} U^{(\mathbf k)}_{mn}\,\psi_{m\mathbf k}(\mathbf r)
$$

其中 $U^{(\mathbf k)}$ 是任意 $J\times J$ 幺正矩阵，也就是“规范”。不同规范给出形状、局域性完全不同的 Wannier 函数，但它们张成同一个子空间。实际计算中积分换成 $N_k$ 个 k 点的离散求和，得到的 Wannier 函数在 $N_1\times N_2\times N_3$ 的超胞内周期重复。

两种方法的区别就在于怎样选 $U^{(\mathbf k)}$：

- **MLWF**：选使 spread 最小的规范；
- **CWF**：选使 Wannier 函数与给定引导函数（PAO、杂化轨道、分子轨道）“距离”最小的规范。

### 2.2 MLWF：spread 泛函与两步优化

Marzari–Vanderbilt（1997）定义 spread

$$
\Omega=\sum_{n}\left[\langle w_{n\mathbf 0}|r^2|w_{n\mathbf 0}\rangle-|\langle w_{n\mathbf 0}|\mathbf r|w_{n\mathbf 0}\rangle|^2\right]=\Omega_{\mathrm I}+\Omega_{\mathrm D}+\Omega_{\mathrm{OD}}
$$

$\Omega_{\mathrm I}$ 与规范无关，$\Omega_{\mathrm D}$、$\Omega_{\mathrm{OD}}$（对角、非对角部分）依赖规范。在倒空间里，$\Omega$ 由相邻 k 点 Bloch 函数周期部分的重叠矩阵

$$
M^{(\mathbf k,\mathbf b)}_{mn}=\langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle
$$

的有限差分表示，$\mathbf b$ 是连接近邻 k 点的矢量。OpenMX 按距离一壳一壳地搜索 $\mathbf b$，搜索的最大壳数由 `Wannier.MaxShells` 控制。

能带分两种情况：

1. **孤立能带组**（如 Si 的 4 条价带）：能带数 = Wannier 函数数，只需最小化 $\Omega_{\mathrm D}+\Omega_{\mathrm{OD}}$；
2. **纠缠能带**（如 Si 的导带、石墨烯的 π*）：先在外窗口（outer window）内的 $N_k^{\mathrm{win}}\ge J$ 个态中挑出一个 $J$ 维子空间，使 $\Omega_{\mathrm I}$ 最小，即 Souza–Marzari–Vanderbilt（2001）的解纠缠（disentanglement）；内窗口（inner/frozen window）内的态被强制保留在子空间中，因此这一能量范围内的能带会被精确复现。第二步再在这个子空间内最小化规范相关部分。

OpenMX 的标准输出用 `DISE` 标记第一步（$\Omega_{\mathrm I}$ 的自洽迭代），用 `CONV`、`SPRD` 标记第二步（梯度模与 spread 各分量），用 `CENT` 标记 Wannier 中心与单个 spread（单位 Å、Å²）。

### 2.3 CWF：距离函数与一次 SVD

Ozaki（2024）提出的 CWF 选择与引导函数 $Q_p$ 在 Hilbert 空间中距离最小的 Wannier 函数：

$$
F=\sum_{p}\left\| w_{p\mathbf 0}-Q_{\mathbf 0 p}\right\|^2
$$

先用窗口函数 $w(\varepsilon)$ 给每个 Bloch 态加权，构造投影矩阵

$$
a_{\mathbf k\mu,p}=w(\varepsilon_{\mathbf k\mu})\,\langle\phi_{\mathbf k\mu}|Q_{\mathbf 0p}\rangle
$$

对它做奇异值分解 $A=V\Sigma W^\dagger$，最优规范就是 $VW^\dagger$（即 $A$ 的 Löwdin 正交化）。整个过程是**确定性的**，不需要迭代，所以计算量远小于 MLWF；解纠缠通过平滑窗口函数自动完成，并且总是在进行。手册给出的经验判据是：每个 CWF 的 $F$ 约 0.16 算足够小，$F>0.7$ 说明引导函数与所选 Bloch 态的对称性明显不符。OpenMX 4.0 在标准输出和 `.out` 中以 `DM_func` 打印这个量，同时给出“每个 CWF、每个自旋”的值。

多项式窗口 `Poly` 由四个能量 $\varepsilon_0<\varepsilon_1<\varepsilon_2<\varepsilon_3$ 定义：$\varepsilon\le\varepsilon_0$ 与 $\varepsilon>\varepsilon_3$ 处权重严格为 0，$\varepsilon_1$ 到 $\varepsilon_2$ 之间为 1，两侧用四次多项式平滑过渡（值与一阶导连续）。Fermi 窗口 `Fermi` 是两个 Fermi 函数的乘积再加一个小常数 $\delta$，权重处处非零。能量都以化学势为零点。

### 2.4 Wannier 插值

在 Wannier 基下，实空间哈密顿量

$$
H_{mn}(\mathbf R)=\frac{1}{N_k}\sum_{\mathbf k}e^{-i\mathbf k\cdot\mathbf R}\left[U^{(\mathbf k)\dagger}H(\mathbf k)U^{(\mathbf k)}\right]_{mn}
$$

只有少数几个近邻 $\mathbf R$ 的矩阵元明显不为零。于是任意 $\mathbf k'$ 处的能带可以通过 $H(\mathbf k')=\sum_{\mathbf R}e^{i\mathbf k'\cdot\mathbf R}H(\mathbf R)$ 对一个 $J\times J$ 小矩阵对角化得到。这就是 Wannier 插值，也是 Wannier90 计算 AHC、光电导、费米面等需要稠密 k 网格的量的基础。插值是否可靠，取决于 $H(\mathbf R)$ 是否足够短程，所以检查跳跃积分随距离的衰减很重要。

### 2.5 MLWF 与 CWF 对比

| 项目 | MLWF | CWF |
|---|---|---|
| 变分目标 | spread $\Omega$ 最小 | 与引导函数的距离 $F$ 最小 |
| 求解 | 解纠缠 SCF + SD/CG 迭代，可能陷入局部极小 | 一次 SVD，结果确定 |
| 初猜/引导 | `<Wannier.Initial.Projectors`：PAO 及其杂化（sp、sp2、sp3、sp3dz2、sp3deg），可放在任意位置、任意取向 | PAO（AO）、在位杂化轨道（HO）、分子轨道（MO） |
| 解纠缠 | 外窗口 + 内（冻结）窗口 | 平滑窗口函数（Poly/Fermi） |
| k 点 | `Wannier.Kgrid`（独立于 SCF） | 直接使用 `scf.Kgrid`，目前只支持奇数 |
| 主要输出 | `.Wannier_Band`、`.HWR`、`.mlwf*_r.cube`、`.mmn/.amn/.eigen` | `.CWF.Band`、`.CWF.Hop`、`.CWF.Dis_vs_H`、`.cwf_*.cube`、CWF 布居 |
| Wannier90 接口 | `Wannier90.fileout on` | `CWF.Wannier90.fileout on` |
| 非共线/SOC | 支持（官方示例 `GaAs_NC.dat`、`NC_CVO.dat`） | 支持（官方示例 `Bi2Se3-wSO-2.dat`，输出 `.nccube`） |

---

## 3. 工作流程总览

1. 普通 SCF：先收敛 k 点和截断能；
2. 能带 + PDOS：确定目标能带的能量范围、主要轨道成分，以及是否与其他能带纠缠；
3. 选方法并写关键词：MLWF 需要 `Wannier.Func.Num`、外/内窗口、投影子、`Wannier.Kgrid`；CWF 需要 `CWF.Guiding.Orbitals`、`<CWF.Guiding.AOs`、窗口函数；
4. 运行：SCF 与 Wannier 可以在同一次计算中完成，反复调参时用重启文件；
5. 检查：插值能带 vs 原始能带、spread 或 DM_func、中心位置、t 随距离的衰减、cube 形状；
6. 输出 Wannier90 文件，交给 wannier90.x、postw90.x 或 WannierTools。

| 输出文件 | 方法 | 内容 |
|---|---|---|
| `System.Name.Band` | 两者 | 原始能带（`Band.dispersion on`），Hartree，未平移 |
| `System.Name.Wannier_Band` | MLWF | 插值能带，格式同 `.Band`（`Wannier.Interpolated.Bands on`）。注意大写 `B`：手册 47.2.2 节写作 `.Wannier_band`，但 4.0 源码 `Generate_Wannier.c` 实际写出的是 `.Wannier_Band`，运行时标准输出也会打印 `Band data file ...` 给出完整文件名 |
| `System.Name.CWF.Band` | CWF | 插值能带，格式同 `.Band` |
| `System.Name.HWR` | MLWF | Wigner–Seitz 超胞内的 $\langle m,\mathbf 0|H|n,\mathbf R\rangle$（Hartree） |
| `System.Name.CWF.Hop` / `.CWF.Dis_vs_H` | CWF | 二进制跳跃积分 / 文本格式（距离 Å、跳跃积分 eV） |
| `System.Name.mmn` / `.amn` / `.eigen` / `.eig` / `.win` | MLWF | 重叠矩阵、投影矩阵、本征值与本征态（`.mmn` 与 `.eigen` 供重启读入）；`.eig`（eV，相对费米能级）与 `.win` 是给 Wannier90 用的 |
| `System.Name.mlwf{spin}_{n}_r.cube`（及 `_i`） | MLWF | MLWF 的实部/虚部 |
| `System.Name.cwf_{spin}_{atom}_{n}_r.cube` / `.cwf_{atom}_{n}.nccube` | CWF | 共线 / 非共线 CWF |
| `System.Name.CWF_MOs_info` | CWF（MO） | 各原子组分子轨道的在位能与占据数 |

`.Band` 系列文件可以用 bandgnu13 转成 gnuplot 数据（见 [第8章](08-band.md)），也可以用本章示例中的 [../examples/14-si-mlwf-sp3/band_compare.py](../examples/14-si-mlwf-sp3/band_compare.py) 直接定量比较两套能带，或者用仓库脚本 [../scripts/omx_band_plot.py](../scripts/omx_band_plot.py) 作图（它的 `--compare` 选项可以叠加第二个 `.Band` 格式文件；用 `python ../scripts/omx_band_plot.py -h` 查看用法）。

---

## 4. MLWF 关键词详解

下表的默认值来自 4.0 源码 `Input_std.c`。手册中有几处写着“无默认值”或与源码不同，建议在输入中**全部显式写出**。

| 关键词 | 推荐值 / 源码默认 | 说明 |
|---|---|---|
| `Wannier.Func.Calc` | `on` / `off` | 总开关 |
| `Wannier.Func.Num` | 目标 WF 数 / 0（手册：无默认值） | `Wannier.Initial.Guess on` 时必须等于投影子总数（sp3 计 4 个，p 计 3 个……），否则报错 |
| `Wannier.Outer.Window.Bottom/Top` | 按能带图设 / −10.0、10.0（手册：无默认值） | eV，相对费米能级；每个 k 点窗口内能带数必须 ≥ `Wannier.Func.Num` |
| `Wannier.Inner.Window.Bottom/Top` | 按需要 / 0.0、0.0 | eV；上下限相等表示不用内窗口；必须落在外窗口内；每个 k 点内窗口能带数必须 ≤ `Wannier.Func.Num` |
| `Wannier.Initial.Guess` | `on` / `on` | 用投影子做初猜，强烈建议开 |
| `Wannier.Initial.Projectors.Unit` | `FRAC` 或 `ANG` / `Ang` | 投影子中心坐标的单位：`FRAC`、`ANG`、`AU` |
| `<Wannier.Initial.Projectors` | 见 4.1 | 每行：`物种-轨道  中心(3)  局部 z 轴(3)  局部 x 轴(3)` |
| `Wannier.Kgrid` | 8 8 8（体相）/ 4 4 4 | MP 网格；三个方向的 k 间距应大致相同 |
| `Wannier.MaxShells` | 12 / 12 | 搜索 $\mathbf b$ 矢量的最大壳数，找不到时增大 |
| `Wannier.Dis.SCF.Max.Steps` | 1000–2000 / 200 | 解纠缠最大步数 |
| `Wannier.Dis.Conv.Criterion` | 1e-10 / 1e-8 | 作用于 `Delta_I`（相邻两步 $\Omega_{\mathrm I}$ 之差，Å²） |
| `Wannier.Dis.Mixing.Para` | 0.5 / 0.5 | 解纠缠中输入/输出子空间投影的混合系数 |
| `Wannier.Minimizing.Scheme` | 2 / 0 | 0 = 最速下降 SD，1 = 共轭梯度 CG，2 = 先 SD 后 CG |
| `Wannier.Minimizing.StepLength` | 2.0 / 2.0 | SD 步长；振荡时减小（官方 `wf_example/Si.dat` 用 0.1，`Secant.StepLength` 用 0.5） |
| `Wannier.Minimizing.Secant.Steps` | 5 / 5 | CG 中割线法线搜索的最大步数 |
| `Wannier.Minimizing.Secant.StepLength` | 2.0 / 2.0 | 割线法初始步长 |
| `Wannier.Minimizing.Max.Steps` | 200–500 / 200 | 最大优化步数；混合方案中 SD、CG 各自使用这个上限 |
| `Wannier.Minimizing.Conv.Criterion` | 1e-10 / 1e-8 | 作用于梯度模 `Mode of Gradient` |
| `Wannier.Interpolated.Bands` | `on` / `off` | 输出 `.Wannier_Band`；必须同时 `Band.dispersion on` |
| `Wannier.Function.Plot` | `on` / `off` | 输出 MLWF 的 cube 文件 |
| `Wannier.Function.Plot.SuperCells` | 1 1 1 / 0 0 0 | 作图超胞 $(2l_1+1)\times(2l_2+1)\times(2l_3+1)$；要求 $2l_i+1\le$ `Wannier.Kgrid` |
| `Wannier.Readin.Overlap.Matrix` | 首次 `off`，重算 `on` / `off` | `on` 时读入 `.mmn` 与 `.eigen`，跳过重叠矩阵计算 |
| `Wannier90.fileout` | 需要时 `on` / `off` | 为 Wannier90 准备 `.amn/.mmn/.eig/.win`。4.0 源码中打开它会把 `Wannier.Dis.SCF.Max.Steps` 与 `Wannier.Minimizing.Max.Steps` 强制设为 0，即 OpenMX 不再做自己的局域化（见 9.1 节） |

> [!WARNING]
> 手册 47.2.1 节的示例注释把 `Wannier.Readin.Overlap.Matrix` 写成“default is on”，但 4.0 源码与官方示例 `work/wf_example/Si.dat` 的默认值都是 `off`。另外，源码在该关键词为 `on` 时会把 `scf.maxIter` 强制设为 1。Wannier 部分的本征值、本征态和费米能级从 `.eigen` 读入，不受影响；但若不同时设 `scf.restart on`，这一步 SCF 是从初始猜测只迭代 1 次，同一次运行输出的 `.Band` 等结果就没有意义。所以请同时设 `scf.restart on`，从已收敛的重启文件出发（见 5.5 节）。

### 4.1 投影子：物种、轨道、位置与取向

投影子用 PAO 构造，需要先在 `<Definition.of.Atomic.Species` 中额外定义一个“投影子物种”。它不出现在原子坐标中，所写的 VPS 也不会参与计算，只是为了让输入格式一致：

```text
Species.Number              2
<Definition.of.Atomic.Species
  Si     Si7.0-s2p2d1       Si_PBE19
  proj   Si7.0-s1p1d1f1     Si_PBE19
Definition.of.Atomic.Species>
```

手册要求每个角动量只用一个径向函数，并推荐一律写 `s1p1d1f1`。投影子行的格式为：

```text
<Wannier.Initial.Projectors
  proj-sp3   0.00 0.00 0.00    0.0 0.0  1.0    1.0 0.0  0.0
  proj-sp3   0.25 0.25 0.25   -1.0 0.0  0.0    0.0 0.0 -1.0
Wannier.Initial.Projectors>
```

- `proj-sp3`：物种名与轨道名用连字符连接；
- 接下来三个数是中心坐标，单位由 `Wannier.Initial.Projectors.Unit` 决定，投影子可以放在晶胞内任意位置（不必在原子上）；
- 再三个数是局部 z 轴，最后三个数是局部 x 轴（均为笛卡尔分量），用来旋转轨道取向。

可用的轨道名（手册表 10，表外名称不能用）：

| 轨道名 | 个数 | 说明 |
|---|---|---|
| `s`、`px`、`py`、`pz`、`dz2`、`dx2-y2`、`dxy`、`dxz`、`dyz` | 1 | 单个 PAO |
| `p` / `d` / `f` | 3 / 5 / 7 | 整个壳层 |
| `fz3`、`fxz2`、`fyz2`、`fzx2`、`fxyz`、`fx3-3xy2`、`f3yx2-y3` | 1 | 单个 f 轨道 |
| `sp` / `sp2` / `sp3` | 2 / 3 / 4 | s 与 p 的杂化 |
| `sp3dz2` / `sp3deg` | 5 / 6 | 三角双锥 / 八面体型杂化 |

> [!TIP]
> 金刚石结构中，`sp3` 的四个瓣在默认坐标系下指向 (1,1,1)、(1,−1,−1)、(−1,1,−1)、(−1,−1,1)，正好对着原点处原子的四个近邻。位于 $a/4(1,1,1)$ 的原子的键方向是它们的反向，所以第二行把局部 z 轴设为 $-\hat x$、局部 x 轴设为 $-\hat z$：这个旋转把 (1,1,1) 映射到 (−1,−1,−1)。取向写错不会报错，只会让初始投影很差、spread 难以收敛。

---

## 5. 实战一：Si 的 sp3 型 MLWF（Si_MLWF_sp3.dat）

目标：从 4 条价带和最低的 4 条导带构造 8 个 sp3 型 MLWF。导带与更高的能带纠缠，所以需要解纠缠。完整输入如下（示例目录中的同名文件）：

```text
System.CurrrentDirectory    ./
System.Name                 si_mlwf_sp3
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            1

Species.Number              2
<Definition.of.Atomic.Species
  Si     Si7.0-s2p2d1       Si_PBE19
  proj   Si7.0-s1p1d1f1     Si_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Si  0.00  0.00  0.00   2.0 2.0
  2  Si  0.25  0.25  0.25   2.0 2.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  0.000  2.715  2.715
  2.715  0.000  2.715
  2.715  2.715  0.000
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            220.0
scf.maxIter                 100
scf.EigenvalueSolver        Band
scf.Kgrid                   10 10 10
scf.Generation.Kpoint       GCenter
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.criterion               1.0e-9

MD.Type                     Nomd

Band.dispersion             on
<Band.KPath.UnitCell
  5.430  0.000  0.000
  0.000  5.430  0.000
  0.000  0.000  5.430
Band.KPath.UnitCell>
Band.Nkpath                 5
<Band.kpath
  26   0.50 0.50 0.50    0.00 0.00 0.00    L  G
  30   0.00 0.00 0.00    1.00 0.00 0.00    G  X
  15   1.00 0.00 0.00    1.00 0.50 0.00    X  W
  11   1.00 0.50 0.00    0.75 0.75 0.00    W  K
  32   0.75 0.75 0.00    0.00 0.00 0.00    K  G
Band.kpath>

Wannier.Func.Calc                     on
Wannier.Func.Num                      8
Wannier.Outer.Window.Bottom         -14.0
Wannier.Outer.Window.Top             12.0
Wannier.Inner.Window.Bottom         -14.0
Wannier.Inner.Window.Top              1.0
Wannier.Initial.Guess                 on
Wannier.Initial.Projectors.Unit       FRAC
<Wannier.Initial.Projectors
  proj-sp3   0.00 0.00 0.00    0.0 0.0  1.0    1.0 0.0  0.0
  proj-sp3   0.25 0.25 0.25   -1.0 0.0  0.0    0.0 0.0 -1.0
Wannier.Initial.Projectors>
Wannier.Kgrid                         8 8 8
Wannier.MaxShells                     12
Wannier.Dis.SCF.Max.Steps             2000
Wannier.Dis.Conv.Criterion            1.0e-10
Wannier.Dis.Mixing.Para               0.5
Wannier.Minimizing.Scheme             2
Wannier.Minimizing.StepLength         2.0
Wannier.Minimizing.Secant.Steps       5
Wannier.Minimizing.Secant.StepLength  2.0
Wannier.Minimizing.Max.Steps          300
Wannier.Minimizing.Conv.Criterion     1.0e-10
Wannier.Interpolated.Bands            on
Wannier.Function.Plot                 on
Wannier.Function.Plot.SuperCells      1 1 1
Wannier.Readin.Overlap.Matrix         off
```

### 5.1 关键参数解释

| 关键词 | 本例取值 | 为什么这样设 |
|---|---|---|
| `Wannier.Func.Num` | 8 | 两个原子各 4 个 sp3 投影子 |
| 外窗口 | −14 ~ 12 eV | 下限在价带底（约 −12 eV）之下；上限要保证每个 k 点至少有 8 条能带。手册的同类示例用 −14 ~ 11 eV |
| 内窗口 | −14 ~ 1 eV | 冻结全部价带和导带底附近的态，使带隙附近的插值能带与原始能带严格一致；每个 k 点内窗口中的能带不超过 6 条，满足 ≤ 8 的要求 |
| `Wannier.Kgrid` | 8 8 8 | 512 个 k 点。网格越密，插值越准，但 `.mmn` 越大 |
| `Wannier.Minimizing.Scheme` | 2 | 先 SD 稳定下降，再 CG 快速收敛 |
| `Wannier.Function.Plot.SuperCells` | 1 1 1 | 在 3×3×3 超胞中画 MLWF；要求 3 ≤ 8 |
| 投影子物种 `proj` | `Si7.0-s1p1d1f1` | 手册推荐写法；VPS 不参与计算 |

> [!IMPORTANT]
> OpenMX 中能量窗口的零点是费米能级（化学势）。对 300 K 下的半导体，化学势落在带隙中间附近，而不是价带顶。所以 Si 的 “0 eV” 大约位于价带顶之上 0.3 eV 处。写窗口前，务必先看一遍 `.Band` 或 DOS（bandgnu13 输出的 `GNUBAND` 已经以化学势为零点，见第 8 章）。

### 5.2 运行

```bash
cd OpenMX/examples/14-si-mlwf-sp3
sed -i "s#^DATA.PATH .*#DATA.PATH                   $HOME/openmx4.0/DFT_DATA19#" *.dat
export OMP_NUM_THREADS=1
mpirun -np 8 openmx Si_MLWF_sp3.dat -nt 1 > si_mlwf_sp3.std 2>&1
# or: sbatch job.slurm   (generic template: ../scripts/job_omx_slurm.sh)
```

SCF 收敛后，OpenMX 在同一次运行中继续计算重叠矩阵、解纠缠和最小化。

### 5.3 检查输出

先确认每个 k 点上外/内窗口内的能带数：

```bash
grep -A 12 "Selected bands within the Outer and Inner Windows" si_mlwf_sp3.std
# the table has one line per Wannier k point (512 here); summarize all of them:
grep -E '^ *[0-9] +\| *[0-9]+\| +\(' si_mlwf_sp3.std | \
  awk -F'|' '{split($3,o,"]"); split($4,w,"]"); n=o[2]+0; m=w[2]+0;
              if (NR==1 || n<nmin) nmin=n; if (m>mmax) mmax=m}
              END {print "min bands in OUTER window:", nmin, "  max bands in INNER window:", mmax}'
```

表头为 `spin| kpt | Outer Window (Nk) |  Inner Window (Mk)  |`，每行给出 `(下标, 上标] 带数 | (下标, 上标] 带数`。外窗口的带数必须 ≥ 8，内窗口的带数必须 ≤ 8，否则程序直接停止（见“常见问题”）。表格按 `Wannier.Kgrid` 的每个 k 点各打印一行，只看前几行不够，用上面的 awk 汇总最小/最大值。然后按手册 47.2.3 节的标记检查优化过程：

```bash
grep "DISE" si_mlwf_sp3.std | tail -n 3        # Omega_I converged, Delta_I below criterion
grep "CONV" si_mlwf_sp3.std | tail -n 4        # gradient norm, "CONVERGENCE ACHIEVED !"
grep "SPRD" si_mlwf_sp3.std | tail -n 4        # Omega_I, Omega_D, Omega_OD, Tot_Omega (Ang^2)
grep "CENT" si_mlwf_sp3.std | tail -n 10       # final centres (Ang) and spreads (Ang^2)
python3 wf_centers.py si_mlwf_sp3.std          # all of the above in one table
```

`wf_centers.py` 打印的表格形如（示意输出，数值仅供参考）：

```text
last SPRD: Omega_I = 16.9xxxx  Omega_D = 0.1xxxx  Omega_OD = 5.7xxxx  Total = 22.8xxxx Ang^2
'CONVERGENCE ACHIEVED' (SPRD) found: yes
 WF        x (Ang)     y (Ang)     z (Ang)   spread (Ang^2)
   1      1.14xxx     1.14xxx     1.14xxx      2.85xxx
 ...
```

判断标准：

- `CONVERGENCE ACHIEVED !` 出现，`Delta_I` 和梯度模都低于设定的判据；
- 由于对称性，8 个 MLWF 的 spread 应两两相等（4 + 4，或全部相等），中心应位于 Si 原子与键中点的连线上；手册用 LDA、$a=5.40$ Å 的同类计算得到每个 WF 约 2.85 Å² 的 spread（手册 47.2.3 节），本例参数不同，数值会有差异，但量级应相同；
- 插值能带与原始能带在内窗口中几乎重合：

```bash
python3 band_compare.py si_mlwf_sp3.Band si_mlwf_sp3.Wannier_Band --emin -14 --emax 1 --plot si_sp3.png
```

内窗口内的 `max |dE|` 应在 meV 量级。窗口之外的高能导带本来就不在 8 维子空间中，不必重合。

### 5.4 可视化 MLWF

`Wannier.Function.Plot on` 会为每个 MLWF 写出 `si_mlwf_sp3.mlwf0_*_r.cube` 与 `_i.cube`（第一个数字是自旋，第二个是 WF 序号，`r`/`i` 表示实部/虚部）。用 VESTA 或 OpenMX Viewer 打开，等值面取 ±0.05 左右，应能看到沿 Si–Si 键伸展、在键中点附近有大振幅的 sp3 型瓣。虚部应接近 0（收敛良好的 MLWF 可以选为实函数）。截面和积分可以用 [`omx_cube_tools.py`](../scripts/omx_cube_tools.py) 处理（用 `python ../scripts/omx_cube_tools.py -h` 查看用法）。

### 5.5 只重做 Wannier 部分：重用重叠矩阵

调整内窗口、初猜或收敛参数时，不必重算 SCF 和重叠矩阵：

```bash
sed -e 's/^Wannier.Readin.Overlap.Matrix .*/Wannier.Readin.Overlap.Matrix  on/' \
    -e 's/^Wannier.Inner.Window.Top .*/Wannier.Inner.Window.Top  0.0/' \
    Si_MLWF_sp3.dat > Si_MLWF_sp3_redo.dat
echo "scf.restart  on" >> Si_MLWF_sp3_redo.dat
mpirun -np 8 openmx Si_MLWF_sp3_redo.dat -nt 1 > si_mlwf_sp3_redo.std 2>&1
```

规则（手册 47.2.1 节）：外窗口和 `Wannier.Kgrid` 必须与生成 `.mmn` 时相同，读入时会做一致性检查；内窗口、投影子和收敛参数可以修改。`System.Name` 不能改，因为 `.mmn`、`.eigen` 和 `si_mlwf_sp3_rst/` 都按它命名。

### 5.6 变体：4 个键中心 MLWF（Si_MLWF_bond4.dat）

只取 4 条价带（孤立能带组）。把外窗口上限设在带隙中（0 eV），不用内窗口，初猜改成放在 4 个键中点的 s 型投影子，并用 `ANG` 单位写笛卡尔坐标：

```text
Wannier.Func.Num                      4
Wannier.Outer.Window.Bottom         -14.0
Wannier.Outer.Window.Top              0.0
Wannier.Inner.Window.Bottom           0.0
Wannier.Inner.Window.Top              0.0
Wannier.Initial.Projectors.Unit       ANG
<Wannier.Initial.Projectors
  proj-s    0.67875  0.67875  0.67875    0.0 0.0 1.0    1.0 0.0 0.0
  proj-s    0.67875 -0.67875 -0.67875    0.0 0.0 1.0    1.0 0.0 0.0
  proj-s   -0.67875  0.67875 -0.67875    0.0 0.0 1.0    1.0 0.0 0.0
  proj-s   -0.67875 -0.67875  0.67875    0.0 0.0 1.0    1.0 0.0 0.0
Wannier.Initial.Projectors>
```

$a/8=0.67875$ Å，四个点是原点处 Si 与四个近邻连线的中点。运行时会出现提示 `Message:WF  The top and bottom of INNER window is the same.`，表示不使用内窗口，这是正常的。检查：

```bash
mpirun -np 8 openmx Si_MLWF_bond4.dat -nt 1 > si_mlwf_bond4.std 2>&1
python3 wf_centers.py si_mlwf_bond4.std \
   --cell 0 2.715 2.715  2.715 0 2.715  2.715 2.715 0 \
   --points "0.67875,0.67875,0.67875;0.67875,-0.67875,-0.67875;-0.67875,0.67875,-0.67875;-0.67875,-0.67875,0.67875"
python3 band_compare.py si_mlwf_bond4.Band si_mlwf_bond4.Wannier_Band --emin -14 --emax 0
```

对于 Si 这种有反演中心、同核成键的体系，4 个键中心 WF 的中心应落在键中点上（`d_nearest_ref` 接近 0），4 个 spread 相等，价带插值与原始能带重合。换成 GaAs（练习 3），中心会向电负性更大的 As 一侧移动，这是 Wannier 中心描述键极性的经典例子。

---

## 6. 实战二：石墨烯 pz 双带模型（Graphene_pz_MLWF.dat）

目标：每个 C 原子一个 pz 型 MLWF，得到描述 π/π* 能带和 Dirac 锥的最小紧束缚模型。π 带与 σ 带在 −8 ~ −3 eV 范围内能量重叠，π* 带在 Γ 附近升得很高并与其他态纠缠，所以这是一个典型的解纠缠问题。

```text
System.CurrrentDirectory    ./
System.Name                 graphene_pz
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            1

Species.Number              2
<Definition.of.Atomic.Species
  C       C6.0-s2p2d1     C_PBE19
  projC   C6.0-s1p1       C_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  C   0.333333333333  0.666666666667  0.50   2.0 2.0
  2  C   0.666666666667  0.333333333333  0.50   2.0 2.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
   2.460000   0.000000   0.000000
  -1.230000   2.130422   0.000000
   0.000000   0.000000  15.000000
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            220.0
scf.maxIter                 100
scf.EigenvalueSolver        Band
scf.Kgrid                   15 15 1
scf.Generation.Kpoint       GCenter
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.criterion               1.0e-9

MD.Type                     Nomd

Band.dispersion             on
Band.Nkpath                 3
<Band.kpath
  40   0.000000 0.000000 0.0    0.500000 0.000000 0.0    G  M
  23   0.500000 0.000000 0.0    0.333333 0.333333 0.0    M  K
  46   0.333333 0.333333 0.0    0.000000 0.000000 0.0    K  G
Band.kpath>

Wannier.Func.Calc                     on
Wannier.Func.Num                      2
Wannier.Outer.Window.Bottom         -10.0
Wannier.Outer.Window.Top             16.0
Wannier.Inner.Window.Bottom          -2.0
Wannier.Inner.Window.Top              1.5
Wannier.Initial.Guess                 on
Wannier.Initial.Projectors.Unit       FRAC
<Wannier.Initial.Projectors
  projC-pz   0.333333333333  0.666666666667  0.50    0.0 0.0 1.0    1.0 0.0 0.0
  projC-pz   0.666666666667  0.333333333333  0.50    0.0 0.0 1.0    1.0 0.0 0.0
Wannier.Initial.Projectors>
Wannier.Kgrid                         12 12 1
Wannier.MaxShells                     12
Wannier.Dis.SCF.Max.Steps             2000
Wannier.Dis.Conv.Criterion            1.0e-10
Wannier.Dis.Mixing.Para               0.5
Wannier.Minimizing.Scheme             2
Wannier.Minimizing.StepLength         2.0
Wannier.Minimizing.Secant.Steps       5
Wannier.Minimizing.Secant.StepLength  2.0
Wannier.Minimizing.Max.Steps          300
Wannier.Minimizing.Conv.Criterion     1.0e-10
Wannier.Interpolated.Bands            on
Wannier.Function.Plot                 on
Wannier.Function.Plot.SuperCells      1 1 0
Wannier.Readin.Overlap.Matrix         off
```

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| 晶格 | $a=2.46$ Å，120° 约定，$c=15$ Å | 这个约定下 K = (1/3, 1/3, 0)，M = (1/2, 0, 0)（倒格矢单位） |
| `projC-pz` | 每个 C 一个 | 只选 pz 一个 m 分量，这是 MLWF 投影子相对 CWF AO 模式的优势 |
| 外窗口 | −10 ~ 16 eV | 下限包住 π 带底（Γ 点约 −8 eV）；上限要把 Γ 附近的 π* 也包进来，否则 Γ 点窗口内没有 pz 成分的态，初始投影接近奇异 |
| 内窗口 | −2 ~ 1.5 eV | 只包含 Dirac 点附近的 π/π*（σ 带顶在 Γ 点约 −3 eV 以下），保证 Dirac 锥被精确复现 |
| `Wannier.Kgrid` | 12 12 1 | 二维体系 c 方向只取 1 个 k 点 |
| `Wannier.Function.Plot.SuperCells` | 1 1 0 | 第三个方向必须为 0，因为 $2l_3+1\le 1$ |

> [!CAUTION]
> 上面的窗口是根据石墨烯能带的一般特征给出的起始值，不是计算结果。用你的基组和泛函先算一遍 `graphene_pz.Band`，确认 σ 带顶和 Γ 点 π* 的位置，再调整内外窗口。如果外窗口内某个 k 点只有 1 条带，或内窗口内某个 k 点超过 2 条带，程序会直接报错（见“常见问题”）。

运行与检查：

```bash
cd OpenMX/examples/14-graphene-pz-mlwf
mpirun -np 8 openmx Graphene_pz_MLWF.dat -nt 1 > graphene_pz.std 2>&1
grep "CONVERGENCE ACHIEVED" graphene_pz.std | tail -n 2
python3 ../14-si-mlwf-sp3/wf_centers.py graphene_pz.std
python3 hwr_hopping.py graphene_pz.HWR --top 12
python3 ../14-si-mlwf-sp3/band_compare.py graphene_pz.Band graphene_pz.Wannier_Band \
        --emin -2.0 --emax 1.5 --plot graphene_pz_bands.png --ylim -10 10
```

看什么：

- 两个 WF 的中心应在 C 原子上（z 方向不偏离原子平面），spread 相等；
- `hwr_hopping.py` 列出的最大跳跃积分应当是 R = (0,0,0) 中的 $H_{12}$（同一晶胞内的两个 C 互为最近邻）以及与之等价的另外两个最近邻，数值为负，量级为 −2.5 ~ −3 eV。紧束缚文献中常用的最近邻 $|t|\approx2.7$–2.8 eV（Castro Neto 等 2009 综述）；
- 两个 pz WF 的在位能 $\varepsilon_0$ 相同（A、B 子格等价）。它与 Dirac 点相差不大但不严格相等：含次近邻跳跃 $t_2$ 时，K 点能量为 $\varepsilon_0-3t_2$（见练习 5）；
- 插值能带在 K 点形成 Dirac 锥，内窗口内与原始能带重合，π 带底附近允许有小偏差。

> [!TIP]
> 如果需要包含 σ 键的完整模型，可以参考官方示例 `work/wf_example/symGra.dat`：每个 C 用 sp2 + pz 投影子，共 8 个 WF，内窗口覆盖全部占据态。在这种情况下比较 spread 和能带复现程度，可以体会“WF 个数越多，模型越完整，但越不简约”。

---

## 7. CWF 关键词详解（OpenMX 4.0）

| 关键词 | 推荐值 / 源码默认 | 说明 |
|---|---|---|
| `CWF.Calc` | `on` / `off` | 总开关 |
| `CWF.Guiding.Orbitals` | `AO` 或 `HO` / `HO` | `AO` = PAO 本身；`HO` = 每个原子上 PAO 的杂化轨道（由在位密度矩阵 SVD 得到）；`MO` = 原子组的分子轨道 |
| `<CWF.Guiding.AOs` | 见下 | AO/HO 模式下选择 PAO：每个物种一行，**每个径向函数都要写 0 或 1** |
| `CWF.window.function` | `Poly` / `Poly` | `Poly`（多项式，窗口外权重严格为 0）或 `Fermi` |
| `CWF.disentangling.Erange` | `Poly` 写 4 个数，`Fermi` 写 2 个数 | eV，以化学势为零点。同名关键词，参数个数随窗口类型而变 |
| `CWF.disentangling.smearing.kBT` | 仅 `Fermi` / 0.3 0.3 | 两侧 Fermi 函数的 $k_BT_0$、$k_BT_1$（eV） |
| `CWF.disentangling.smearing.bound` | 仅 `Fermi` | Fermi 窗口的下限常数 $\delta$，通常不必修改（手册写默认 1.0e-12；4.0 源码 `Input_std.c` 在 `Fermi` 分支中的默认值为 1.0e-10） |
| `CWF.Dis.vs.H` | `on` / `off` | 输出文本格式跳跃积分 `.CWF.Dis_vs_H` |
| `CWF.fileout.flag` | 需要时 `on` / `off` | 输出 CWF 的 cube（要求 `CWF.Calc on`） |
| `CWF.Plot.SuperCells` | 1 1 1 / 0 0 0 | cube 的超胞 $(2l_i+1)$；要求 $2l_i+1\le$ `scf.Kgrid` |
| `CWF.fileout.Num` | 原子数或组数 | AO/HO 时是原子数，MO 时是组数 |
| `<CWF.fileout.cube` | 原子序号或组序号 | 个数必须等于 `CWF.fileout.Num` |
| `<CWF.MO.Grouped.Atoms` | MO 模式 | 每个原子一行：`原子序号  组号` |
| `<CWF.MO.Selection` | MO 模式 | 每组一行：`组号  选取的分子轨道个数` |
| `CWF.Wannier90.fileout` | 先 `off`，确认能带复现后再 `on` / `off` | 写出 `System.Name_{spin}.amn/.mmn/.eig/.win` |

两条硬性限制：

1. **CWF 使用 `scf.Kgrid` 作为 Wannier 的 k 网格，并且目前只支持奇数 k 点**（手册 47.1.7 节，这是实现上的限制）；
2. `CWF.Guiding.AOs` 按径向函数整体选择，`1p 1` 表示 px、py、pz 三个都选上，没有办法只选 pz。需要单个 m 分量时，用 MLWF 的投影子（第 6 节），或者考虑 MO/HO 模式。

`CWF.Guiding.AOs` 的写法必须与 `Definition.of.Atomic.Species` 中的基组逐项对应。例如 `Ga7.0-s3p2d2` 有 `1s 2s 3s 1p 2p 1d 2d` 共 7 项：

```text
<CWF.Guiding.AOs
  Ga  1s 1 2s 0 3s 0 1p 1 2p 0 1d 0 2d 0
  As  1s 1 2s 0 3s 0 1p 1 2p 0 1d 0 2d 0
CWF.Guiding.AOs>
```

对于把 3d 当作半芯价态的 `Ga_PBE19`，`1s`、`1p`、`1d` 分别对应 Ga-4s、4p、3d（最低的那个径向函数）。CWF 个数等于选中的轨道数之和：$(1+3)\times 2=8$。

> [!NOTE]
> 源码中，`<CWF.Guiding.AOs` 会按 `Species.Number` 的数目读取行数，所以每个物种都要写一行，包括全选 0 的物种。块内不要插入注释行，结束标记要紧跟最后一行数据，否则会报 `Format error for CWF.Guiding.AOs`。

---

## 8. 实战三：GaAs 的 CWF（GaAs_CWF_AO.dat / GaAs_CWF_HO.dat）

目标：以 Ga-4s、4p 和 As-4s、4p 为引导函数，构造张成 4 条价带（含 As-4s 带）和最低导带的 8 个 CWF；再改用杂化轨道引导，得到指向近邻原子的“键方向”CWF。GaAs 中 As-3d（约 −36 eV）、Ga-3d（约 −15 eV）、As-4s（约 −11 eV）都是孤立能带（手册 47.1.1 节），本例用多项式窗口排除前两者。

### 8.1 完整输入

```text
System.CurrrentDirectory    ./
System.Name                 gaas_cwf_ao
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            1

Species.Number              2
<Definition.of.Atomic.Species
  Ga   Ga7.0-s3p2d2   Ga_PBE19
  As   As7.0-s3p2d2   As_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Ga   0.00  0.00  0.00   6.5  6.5
  2  As   0.25  0.25  0.25   7.5  7.5
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  2.8265   2.8265   0.0000
  2.8265   0.0000   2.8265
  0.0000   2.8265   2.8265
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0
scf.energycutoff            220.0
scf.maxIter                 100
scf.EigenvalueSolver        Band
scf.Kgrid                   11 11 11
scf.Generation.Kpoint       GCenter
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.criterion               1.0e-9
scf.restart                 off

MD.Type                     Nomd

Band.dispersion             on
Band.Nkpath                 4
<Band.kpath
  30  0.500 0.500 0.500   0.000 0.000 0.000   L  G
  35  0.000 0.000 0.000   0.500 0.500 0.000   G  X
  12  0.500 0.500 0.000   0.625 0.625 0.250   X  U
  37  0.750 0.375 0.375   0.000 0.000 0.000   K  G
Band.kpath>

CWF.Calc                    on
CWF.Guiding.Orbitals        AO
<CWF.Guiding.AOs
  Ga  1s 1 2s 0 3s 0 1p 1 2p 0 1d 0 2d 0
  As  1s 1 2s 0 3s 0 1p 1 2p 0 1d 0 2d 0
CWF.Guiding.AOs>
CWF.window.function         Poly
CWF.disentangling.Erange    -14.5  -13.8  1.0  8.0
CWF.Dis.vs.H                on
CWF.Wannier90.fileout       off
CWF.fileout.flag            on
CWF.Plot.SuperCells         1 1 1
CWF.fileout.Num             2
<CWF.fileout.cube
  1
  2
CWF.fileout.cube>
```

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `scf.Kgrid` | 11 11 11 | 同时用作 CWF 的 k 网格，必须是奇数 |
| `CWF.Guiding.AOs` | Ga/As 的 `1s`、`1p` | 4s + 4p，每个原子 4 个，共 8 个 CWF |
| `CWF.disentangling.Erange` | −14.5 −13.8 1.0 8.0 | ≤ −14.5 eV 权重为 0，排除 Ga-3d 与 As-3d；−13.8 ~ 1.0 eV 权重为 1，完整包含 As-4s 带和全部价带；1 ~ 8 eV 平滑衰减，导带底附近的态以较大权重参与 |
| `CWF.Dis.vs.H` | `on` | 便于分析 t 随距离的衰减 |
| `CWF.Plot.SuperCells` | 1 1 1 | 3 ≤ 11，满足要求 |

> [!WARNING]
> Ga-3d 带的位置取决于基组和泛函。正式计算前先看一遍 PDOS（[第9章](09-dos-pdos-cohp-optics.md)）或 `.Band`，确认 Ga-3d 整条带都在 $\varepsilon_0$ 之下、As-4s 带底在 $\varepsilon_1$ 之上。如果 $\varepsilon_0$ 切进了 Ga-3d 带，部分 3d 态会以小权重混进来，DM_func 会变大。
>
> 参考数据：官方参考输出 `work/cwf_example/GaAs-CWF1.out`（PBE，`Ga7.0-s3p2d2f1`/`As7.0-s3p2d2f1`，$a=5.65$ Å，13×13×13 网格）中列出了全部 SCF k 点的本征值。按其中的化学势换算后，Ga-3d 带的上沿约为 −15.2 eV，As-4s 带约在 −13.1 ~ −10.5 eV（相对化学势）。本例的 $\varepsilon_0=-14.5$、$\varepsilon_1=-13.8$ eV 与两者各留约 0.7 eV 余量；基组或晶格常数不同时请自行核对。

### 8.2 运行与检查

```bash
cd OpenMX/examples/14-gaas-cwf
mpirun -np 8 openmx GaAs_CWF_AO.dat -nt 1 > gaas_cwf_ao.std 2>&1
grep "DM_func" gaas_cwf_ao.std gaas_cwf_ao.out | tail -n 2
grep -A 12 "Populations evaluated by CWF" gaas_cwf_ao.out
python3 cwf_analyze.py gaas_cwf_ao --plot gaas_cwf_ao_t.png
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ao.Band gaas_cwf_ao.CWF.Band --emin -13.8 --emax 1.0
ls gaas_cwf_ao.cwf_0_*_r.cube | head
```

判断标准：

- `DM_func` 第二个数（每个 CWF、每个自旋）越小越好。手册中以全部价态为引导的 GaAs 例子为 0.16，可作为“足够小”的参照；明显大于 0.7 时，说明引导函数与窗口内态的对称性不匹配；
- 插值能带 `.CWF.Band` 在权重为 1 的区间（−13.8 ~ 1 eV）与 `.Band` 重合，导带越往上偏差越大，这是正常的；
- `.CWF.Dis_vs_H` 中的 |t| 随距离迅速衰减（`cwf_analyze.py` 按距离壳层列出最大 |t|），说明 CWF 在实空间是局域的；
- cube 文件 `gaas_cwf_ao.cwf_0_{atom}_{n}_r.cube` 中，各 CWF 保持了 s、px、py、pz 的形状。共线计算中 CWF 总是实函数。

CWF 布居会写在 `.out` 的 `Populations evaluated by CWF` 段中，格式与 Mulliken 布居类似。本例只选了 4s/4p，3d 半芯电子不在 CWF 子空间里，所以 `total` 会明显小于 `ideal(neutral)=28`，这不是错误。若想做完整的电荷分析，应像手册的 GaAs-CWF1 那样把所有价态（包括 3d）都选为引导函数，并把窗口下限放到 As-3d 之下。

### 8.3 杂化轨道引导（GaAs_CWF_HO.dat）与键轨道图像

`GaAs_CWF_HO.dat` 只改了 `CWF.Guiding.Orbitals HO`（以及 `System.Name`）。HO 模式仍然用 `<CWF.Guiding.AOs` 选出每个原子上参与杂化的 PAO，再由在位密度矩阵的 SVD 得到杂化轨道。对 GaAs 来说，得到的是 4 个 sp3 型杂化轨道，分别指向 4 个近邻：

```bash
mpirun -np 8 openmx GaAs_CWF_HO.dat -nt 1 > gaas_cwf_ho.std 2>&1
python3 cwf_analyze.py gaas_cwf_ho
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ho.Band gaas_cwf_ho.CWF.Band --emin -13.8 --emax 1.0
```

预期（手册 47.1.2 节的结论）：选用相同的 PAO 时，AO 与 HO 张成**同一个子空间**，所以插值能带完全相同；不同的是 CWF 的形状，HO 型 CWF 沿键方向伸展。Ga 与 As 上指向同一根键的两个杂化 CWF，其成键/反键组合就对应 sp3 键轨道图像：成键组合主要落在价带，反键组合落在导带。在 `.CWF.Dis_vs_H` 的非在位项（距离 ≠ 0）中，“同一根键上两个杂化轨道之间的跳跃积分”预期是绝对值最大的一项，可以理解为键的共价强度；请用 `cwf_analyze.py` 的距离壳层表自行确认，本仓库不附参考数值。

### 8.4 用 scf.restart pp 扫描窗口

窗口参数通常要试几组。4.0 的 `scf.restart pp` 读入 `System.Name_rst/` 后只做一次对角化，然后进行后处理：程序自动设置 `scf.maxIter=1`、`MD.Type=nomd`，并且不覆盖重启文件（手册 17.1 节）。示例脚本 [scan_cwf_window.sh](../examples/14-gaas-cwf/scan_cwf_window.sh) 就按这个思路，对几组 $(\varepsilon_0,\varepsilon_1,\varepsilon_2,\varepsilon_3)$ 依次计算，记录 `DM_func`。因为 `System.Name` 不能改，每一步都会覆盖原始计算的 `.out`、`.CWF.Band` 等文件，所以脚本先把它们备份为 `gaas_cwf_ao_orig.*`，扫描结束后再恢复，每个窗口的结果另存为 `gaas_cwf_ao_w<i>.*`：

```bash
mpirun -np 8 openmx GaAs_CWF_AO.dat -nt 1 > gaas_cwf_ao.std 2>&1   # creates gaas_cwf_ao_rst/
NP=8 bash scan_cwf_window.sh
cat cwf_window_scan.txt
```

脚本中最后一组窗口（−10 −9 …）把 As-4s 带排除在外，但引导函数里仍然有 As-4s，预期 DM_func 会明显变大。这正好用来体会“引导函数个数和窗口内能带要匹配”。

> [!NOTE]
> 手册没有逐条列出 `pp` 模式支持哪些后处理。从 4.0 源码（`Input_std.c`）看，`pp` 等价于 `scf.restart on` 再加上 `scf.maxIter=1`、`MD.Type=nomd` 和“不写回重启文件”，SCF 之后的后处理（包括 CWF）照常执行。如果你的版本在 `pp` 下没有输出 CWF 结果，就改用 `scf.restart on`：从收敛的重启文件出发，SCF 一两步就会结束，然后照常进行 CWF 计算（注意这样会覆盖重启文件）。

### 8.5 分子轨道引导（MO）与非共线 CWF

分子晶体的单胞原子很多，而费米面附近往往只有几条带。MO 模式先把原子分组，用窗口内的态为每组构造密度矩阵，再通过 SVD 得到“类分子轨道”作为引导函数：

```text
CWF.Guiding.Orbitals        MO
CWF.window.function         Poly
CWF.disentangling.Erange    -0.6 -0.5 0.5 0.6
<CWF.MO.Grouped.Atoms
  1  1
  2  1
  3  2
  4  2
CWF.MO.Grouped.Atoms>
<CWF.MO.Selection
  1  1
  2  1
CWF.MO.Selection>
```

上面是一个 4 原子、分成 2 组、每组取 1 个分子轨道的格式示意。`<CWF.MO.Grouped.Atoms` 必须按原子序号 1, 2, 3, … 逐行列出全部原子；`<CWF.MO.Selection` 必须按组号逐行列出全部组。每组取几个分子轨道，可以先看 `System.Name.CWF_MOs_info` 中各分子轨道的占据数（第 3 列）来决定。手册强调：选取的引导函数总数应等于窗口内的目标能带数，引导函数多于能带数时 SVD 问题不确定。完整的例子见官方 `work/cwf_example/Rb3C60-CWF3.dat`（C60 的 3 条 t1u 带）和 `TTF-TCNQ-3.dat`。

非共线与 SOC 计算（`scf.SpinPolarization nc`、`scf.SpinOrbit.Coupling on`）同样支持 CWF。此时 `.CWF.Dis_vs_H` 没有自旋列，cube 输出为 `System.Name.cwf_{atom}_{n}.nccube`，依次存放 α 分量的实部、虚部和 β 分量的实部、虚部，可用 OpenMX Viewer 查看。官方示例 `Bi2Se3-wSO-2.dat` 只选 Se-4p 与 Bi-6p，窗口为 (−7, −6, 3, 6) eV。SOC 与拓扑分析见 [第12章](12-soc-topology.md)。

### 8.6 自己写 Wannier 插值程序：analysis_cwf_col

源码目录提供了两个辅助程序 `analysis_cwf_col.c`（共线）和 `analysis_cwf_noncol.c`（非共线），它们读取二进制的 `System.Name.CWF.Hop`，构造 $H(\mathbf k)$ 并对角化：

```bash
cd openmx4.0/source
make analysis_cwf_col
./analysis_cwf_col /path/to/gaas_cwf_ao.CWF.Hop     # default: eigenvalues at one k point
```

代码写得比较直白，适合作为自己实现 Berry 曲率、输运等 Wannier 插值计算的起点（手册 47.1.9 节）。

---

## 9. 与 Wannier90 / postw90 / WannierTools 对接

### 9.1 MLWF 路线：Wannier90.fileout

在 MLWF 输入中加上 `Wannier90.fileout on`（`Wannier.Func.Calc` 也必须为 `on`），计算结束时会提示 `The input files for Wannier90, ... are successfully generated.`，并生成 `System.Name.amn`、`.mmn`、`.eig`、`.win`。示例 `Si_MLWF_W90.dat` 就是在 `Si_MLWF_sp3.dat` 的基础上做了以下修改：

```text
System.Name                           si_w90
Wannier.Dis.SCF.Max.Steps             0        # minimization left to wannier90.x
Wannier.Minimizing.Max.Steps          0        # minimization left to wannier90.x
Wannier.Interpolated.Bands            off
Wannier.Function.Plot                 off
Wannier.Function.Plot.SuperCells      0 0 0
Wannier90.fileout                     on
```

把 OpenMX 自己的优化步数设为 0 是官方示例 `Si-Wannier90.dat` 的做法：OpenMX 只负责生成矩阵，局域化交给 wannier90.x。实际上 4.0 源码（`Generate_Wannier.c`）在 `Wannier90.fileout on` 时会自动把这两个步数置 0，示例中显式写出只是为了让意图一目了然；所以打开这个开关后，OpenMX 自己不再做解纠缠和最小化，此时的 `.HWR` 与中心/spread 也不代表最大局域化结果。4.0 源码写出的 `.win` 包含 `num_bands`、`num_wann`、`num_iter 15000`、晶胞与原子分数坐标、`dis_win_min/max`（等于外窗口）、`dis_froz_min/max`（等于内窗口）、`mp_grid` 与 k 点列表、`bands_plot T` 与 k 路径，以及 postw90 的 `berry T`、`berry_task kubo`（光电导；改成 `ahc` 即反常霍尔电导）。能量零点是费米能级（`fermi_energy = 0.0`）。

```bash
mpirun -np 8 openmx Si_MLWF_W90.dat -nt 1 > si_w90.std 2>&1
cat w90_extra.win >> si_w90.win      # adds write_hr / write_xyz
wannier90.x si_w90                    # -> si_w90.wout, si_w90_band.dat, si_w90_hr.dat
grep -A 10 "Final State" si_w90.wout  # Wannier centres and spreads
postw90.x si_w90                      # optional: berry_task in si_w90.win
```

> [!TIP]
> 内窗口上下限相等（不用冻结窗口）时，`.win` 中 `dis_froz_min` 与 `dis_froz_max` 也会相等。若想在 Wannier90 中改用其他冻结窗口，可以直接编辑 `.win`；但 `dis_win_min/max` 只能在 OpenMX 生成矩阵时的外窗口范围内缩小，不能扩大。

### 9.2 CWF 路线：CWF.Wannier90.fileout

```text
CWF.Wannier90.fileout       on
```

会生成 `System.Name_0.amn/.eig/.mmn/.win`。`_0` 是自旋序号；共线自旋极化时还有 `_1`；非共线计算的文件名不带自旋后缀。4.0 源码写出的 `.win` 中设置了 `dis_num_iter 0` 和 `num_iter 0`：Wannier90 **直接沿用 CWF 的规范**，不再做解纠缠和局域化。如果想以 CWF 为起点继续做最大局域化，可以把 `num_iter` 改为正数。默认的 postw90 任务是 `berry_task = sc`（shift current），可按需改为 `kubo`、`ahc` 等。

手册特别提醒：`.mmn`、`.amn`、`.eig` 的大小分别随窗口内 Bloch 态数 $N_B$ 按 $O(N_B^2)$、$O(N_B)$、$O(N_B)$ 增长。Fermi 窗口处处权重非零，会把所有 Bloch 态都包含进来，`.mmn` 可能达到几十 GB。所以对接 Wannier90 时要用 `Poly` 窗口，并且先在 `CWF.Wannier90.fileout off` 的情况下确认插值能带合格，再打开它。

### 9.3 WannierTools 与其他下游工具

WannierTools 读取 Wannier90 的 `seedname_hr.dat`：在 `.win` 中加 `write_hr = .true.`，运行 wannier90.x 后，在 WannierTools 的 `wt.in` 里用 `Hrfile` 指定这个文件，就可以计算表面态、Wilson loop、Weyl 点等（输入格式见 WannierTools 文档）。其他常用路线：

| 工具 | 输入 | 用途 |
|---|---|---|
| Wannier90 / postw90 | `.amn/.mmn/.eig/.win` | 插值能带、费米面、AHC、光电导、shift current、BoltzWann 输运 |
| WannierTools | `seedname_hr.dat` | 表面态（迭代格林函数）、Wilson loop、Weyl/节线 |
| [HopTB.jl](https://github.com/HopTB/HopTB.jl) | OpenMX 的 PAO 哈密顿量 | 不经过 Wannier 化，直接在 PAO 紧束缚模型上算 Berry 曲率、光学响应 |
| [HWC](https://github.com/Ncmexp2717/HWC) | `.scfout` | 混合 Wannier 电荷中心（Z2、层极化），见第 12 章 |

---

## 10. 安装后的自动测试

手册提供了 MLWF 与 CWF 的自动测试，在 `work/` 目录下运行：

```bash
mpirun -np 8 openmx -runtestWF  -nt 2     # 8 MLWF cases in work/wf_example  -> runtestWF.result
mpirun -np 8 openmx -runtestCWF -nt 2     # 7 CWF cases in work/cwf_example  -> runtestCWF.result
```

结果文件给出与参考值（Xeon 集群上算得）之差：MLWF 比较 spread，CWF 比较距离函数。小数点后 7 位一致即可认为安装正确。正式计算前跑一次，可以排除编译问题（见 [第2章](02-install.md)）。

---

## 参数选择与收敛建议

| 参数 | 推荐做法 | 说明 |
|---|---|---|
| SCF k 点 / 截断能 | 先按 [第5章](05-scf-kpoints-convergence.md) 收敛 | Wannier 化不会修正 SCF 的误差 |
| 方法选择 | 要轨道图像、快速建模用 CWF；要 spread 定量、文献对比用 MLWF | 两者都可以接 Wannier90 |
| WF 个数 | 等于目标子空间维数；引导函数个数与窗口内能带数匹配 | MLWF 要求投影子总数 = `Wannier.Func.Num` |
| 外窗口 / `Poly` 的 $\varepsilon_0,\varepsilon_3$ | 刚好包住目标能带和与之纠缠的部分；排除孤立的半芯带 | 太宽会使 `.mmn` 变大、解纠缠变慢 |
| 内窗口 / `Poly` 的 $\varepsilon_1,\varepsilon_2$ | 覆盖最关心的能量区间（通常是 $E_F$ 附近） | MLWF 内窗口每个 k 点带数 ≤ WF 数 |
| `Wannier.Kgrid` | 体相 6×6×6 ~ 10×10×10；二维 $N\times N\times1$；三个方向 k 间距相近 | 收敛判据：插值能带与 spread 不再随网格变化 |
| CWF 的 k 网格 | `scf.Kgrid` 取奇数，并检查插值能带对网格的收敛 | 网格太粗时 $H(\mathbf R)$ 截断在小超胞内，插值会出现振荡 |
| 解纠缠收敛 | `Wannier.Dis.Conv.Criterion` 1e-10，`Max.Steps` 1000–2000 | 看 `DISE` 的 `Delta_I` |
| 最小化 | `Scheme 2`；振荡时减小 `StepLength` | 看 `CONV` 的梯度模和 `SPRD` 的 `Tot_Omega` |
| 结果检验 | 插值能带、spread 或 DM_func、中心、t 衰减、cube 形状五项都看 | 任何一项异常都要回到窗口和初猜 |
| 重复调参 | MLWF 用 `Wannier.Readin.Overlap.Matrix on` + `scf.restart on`；CWF 用 `scf.restart pp` | `System.Name` 保持不变 |

---

## 常见问题与报错

下列报错原文取自 4.0 源码（`Input_std.c`、`Generate_Wannier.c`）。

| 现象 / 报错 | 原因 | 解决办法 |
|---|---|---|
| `Error:WF  The number (4) of projectors is not same as that (8) of Wannier functions` | 投影子总数（sp3 = 4、p = 3、d = 5……）与 `Wannier.Func.Num` 不等 | 按第 4.1 节的表数清楚投影子个数 |
| `Error:WF  INNER window (...) is not inside of OUTER window (...).` | 内窗口超出外窗口 | 让内窗口落在外窗口内（边界可以相等） |
| `Error:WF  The top of the OUTER window should be higher than the bottom.` | 外窗口上下限写反或相等 | 改正上下限 |
| `Message:WF  The top and bottom of INNER window is the same.` | 内窗口上下限相等 | 只是提示：表示不使用冻结窗口 |
| `* Bands number within OUTER window [...] * is less than Wannier Function number 8.` | 某些 k 点上外窗口内的能带数少于 WF 数 | 扩大外窗口（通常是提高上限）；用 `grep -A 20 "Selected bands within"` 找出是哪些 k 点 |
| `* Bands number within INNER window [...] * is larger than wannier function number 2.` | 某些 k 点上内窗口内的能带数多于 WF 数 | 缩小内窗口，或者增加 WF 个数 |
| `Error:WF You are requiring for plotting the interpolated bands but without turn on key word Band.dispersion.` | `Wannier.Interpolated.Bands on`，但没开能带 | 加 `Band.dispersion on` 与 `<Band.kpath` |
| `Error:WF, the following condition should be satisfied (2*Wannier.Function.Plot.SuperCells+1)<=Wannier.Kgrid` | 作图超胞大于 Wannier k 网格对应的超胞 | 减小 `Wannier.Function.Plot.SuperCells`；二维体系第三个数写 0 |
| `Error:WF  not found xxx.mmn` | `Wannier.Readin.Overlap.Matrix on`，但没有对应的 `.mmn`（或改了 `System.Name`） | 先用 `off` 算一次；保持 `System.Name` 不变 |
| `Outer window size in Mmn(k,b) file is smaller than that in present calculation.` / `b vectors in Mmn(k,b) file are not the same ...` | 重用 `.mmn` 时改了外窗口或 `Wannier.Kgrid` | 恢复原设置，或者设 `off` 重新计算重叠矩阵 |
| `Can not find proper b vectors, please increase parameter` / `MAXSHELL OR change Wannier.Kgrids.`，或 `Weights for b vectors (totally N) are not found.` | 晶胞长宽比大、k 网格各方向间距相差太多 | 增大 `Wannier.MaxShells`；按倒格矢长度之比设置 `Wannier.Kgrid` |
| `Check Wannier.Kgrids, the equivalent k points for k+b not found.` | k 网格与 b 矢量不匹配 | 按输出中提示的倒格矢长度比例调整 `Wannier.Kgrid` |
| `Format error for Wannier.Initial.Projectors` | 块的结束标记缺失或拼错，或者某行字段不足 10 个 | 每行写满“名称 + 9 个数”；块内不要写注释行 |
| `Could not assign the multiple projector for the same angular part of 'proj'` | 投影子物种的 PAO 在同一角动量上有多个径向函数（例如写成 `Si7.0-s2p2d1`） | 投影子物种只能用 `s1`、`s1p1`、`s1p1d1` 或 `s1p1d1f1`，推荐 `s1p1d1f1` |
| `Error:WF  The top of the INNER window should be higher than the bottom.` | 内窗口上下限写反 | 改正上下限；不想用内窗口时让两者相等 |
| spread 很大、`CONV` 不收敛或来回振荡 | 初猜取向或位置错误；窗口不合理 | 检查局部 z/x 轴；减小 `Wannier.Minimizing.StepLength`；先用 SD（`Scheme 0`） |
| 插值能带在内窗口外大幅振荡 | 外窗口过宽、k 网格太粗、WF 不够局域 | 加密 `Wannier.Kgrid`；收窄外窗口；增加 WF 个数 |
| `Format error for CWF.Guiding.AOs` | 某物种的项数与 PAO 定义不符、漏写物种行、块内有注释 | 按基组逐项写全 0/1，每个物种一行 |
| `Error for CWF: the following condition should be satisfied (2*CWF.Plot.SuperCells+1)<=Kspace_grid` | CWF cube 超胞大于 `scf.Kgrid` 所允许的大小 | 减小 `CWF.Plot.SuperCells` 或加密 `scf.Kgrid` |
| `CWF.Calc should be switched on when CWF.fileout.flag=on.` | 打开了 CWF cube 输出，但没开 CWF | 加 `CWF.Calc on` |
| `Format error for CWF.fileout.cube` | 序号超出原子数（或组数），或个数与 `CWF.fileout.Num` 不一致 | 核对序号和个数 |
| `Format error for CWF.MO.Grouped.Atoms` / `Format error for CWF.MO.Selection` | 原子或组没有按 1, 2, 3 … 顺序逐行列全 | 每个原子、每个组各一行，序号连续 |
| CWF 结果异常或运行出错，`scf.Kgrid` 为偶数 | CWF 目前只支持奇数 k 点 | 把 `scf.Kgrid` 改成奇数 |
| DM_func 每 CWF 大于 0.7 | 引导函数与窗口内态的对称性不符，或个数不匹配 | 用 PDOS 重新确认能带成分；调整 `CWF.Guiding.AOs` 和窗口 |
| `.mmn` 有几十 GB | CWF 用了 `Fermi` 窗口去对接 Wannier90 | 改用 `Poly` 窗口 |
| 3.9 中 CWF 关键词不起作用 | CWF 是 4.0 新功能 | 升级到 4.0（+patch 4.0.1） |

---

## 练习

1. 在 `Si_MLWF_sp3.dat` 中把 `Wannier.Kgrid` 依次设为 4×4×4、6×6×6、8×8×8，记录 `Tot_Omega` 和 `band_compare.py` 在 −14 ~ 1 eV 内的 `max |dE|`，判断网格何时收敛。
2. 用 5.5 节的重算方式，把内窗口上限从 1.0 eV 改为 0.0 eV、−14.0 eV（相当于不冻结），比较带隙附近插值能带的误差和 spread 的变化。
3. 把 `Si_MLWF_bond4.dat` 改成 GaAs（$a=5.653$ Å，换成 `Ga7.0-s3p2d2`/`As7.0-s3p2d2` 和对应电子数，外窗口下限要高于 Ga-3d），计算 4 个价带 MLWF，看它们的中心向哪个原子偏移、偏移多少。
4. 在 `Graphene_pz_MLWF.dat` 中把内窗口改为 −1.0 ~ 1.0 eV，再关闭内窗口（上下限都写 0.0），比较 Dirac 点附近费米速度（K 点附近能带斜率）的复现情况。
5. 用 `hwr_hopping.py` 读出石墨烯的 $t_1$（最近邻）和 $t_2$（次近邻），代入 $E(\mathbf k)=\varepsilon_0\pm t_1|f(\mathbf k)|+t_2(|f(\mathbf k)|^2-3)$ 画出解析能带，与 `.Wannier_Band` 对比。
6. 对 `GaAs_CWF_AO.dat` 运行 `scan_cwf_window.sh`，解释为什么窗口 (−10, −9, 1, 8) 的 DM_func 明显变大；然后把 As 的 `1s` 改为 0（7 个 CWF），再用这个窗口计算，比较结果。
7. 在 `GaAs_CWF_AO.dat` 中打开 `CWF.Wannier90.fileout`，确认 `gaas_cwf_ao_0.win` 中的 `num_iter 0`；运行 wannier90.x 后，比较 `gaas_cwf_ao_0_band.dat` 与 `.CWF.Band`。

---

## 参考资料

官方文档：

- OpenMX 4.0 手册第 47 章 Generation of Wannier Functions：<https://www.openmx-square.org/openmx_man4.0/s47_wannier.html>
- 47.1 CWF（PAO/杂化/分子轨道引导、解纠缠、可视化、k 点、电荷、Wannier90 接口）：<https://www.openmx-square.org/openmx_man4.0/s47_1_cwf.html>
- 47.1.1 PAO 引导函数与窗口函数：<https://www.openmx-square.org/openmx_man4.0/s47_1_1_pao-guiding.html>；47.1.10 CWF–Wannier90 接口：<https://www.openmx-square.org/openmx_man4.0/s47_1_10_interface-w90.html>
- 47.2 MLWF（生成步骤、分析、优化过程检查、输出文件、Wannier90 接口）：<https://www.openmx-square.org/openmx_man4.0/s47_2_mlwf.html>
- 47.2.1 MLWF 生成步骤与关键词：<https://www.openmx-square.org/openmx_man4.0/s47_2_1_generation.html>；47.2.5 MLWF 输出文件格式：<https://www.openmx-square.org/openmx_man4.0/s47_2_5_output-files.html>
- 第 17 章 重启（`scf.restart on|off|pp|c2n`）：<https://www.openmx-square.org/openmx_man4.0/s17_restart.html>
- CWF 计算示例网站：<https://www.openmx-square.org/cwf/>
- What's new（4.0 新功能：CWF）：<https://www.openmx-square.org/whatsnew/whatsnew.html>
- 官方示例：`openmx4.0/work/wf_example/`（`Si.dat`、`GaAs.dat`、`symGra.dat`、`Benzene.dat`、`pmCVO.dat`、`NC_CVO.dat`、`GaAs_NC.dat`、`VBz.dat`），`work/cwf_example/`（`GaAs-CWF1/2/3.dat`、`Cufcc-CWF0.dat`、`Rb3C60-CWF3.dat`、`TTF-TCNQ-3.dat`、`Bi2Se3-wSO-2.dat`），`work/Si-Wannier90.dat`、`work/SrVO3-Wannier90.dat`
- Wannier90 官网与用户手册：<https://wannier.org/>
- WannierTools 文档：<https://www.wanniertools.org/>

相关 GitHub 仓库（只作参考链接，详见 [../../Resources/github-repos.md](../../Resources/github-repos.md)）：

- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90)：Wannier90 源码（MLWF、postw90 的 Berry 模块等）
- [quanshengwu/wannier_tools](https://github.com/quanshengwu/wannier_tools)：WannierTools，读取 `_hr.dat` 计算表面态、Wilson loop
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：从 OpenMX 哈密顿量构建紧束缚模型做 Berry 曲率、光学与非线性响应
- [knhr0727/pymx](https://github.com/knhr0727/pymx)：读 scfout 构建 H(k)，含 Bi2Se3 的 Wannier 电荷中心（Wilson loop）示例，适合理解思路
- [Ncmexp2717/HWC](https://github.com/Ncmexp2717/HWC)：读 scfout 计算混合 Wannier 中心与层极化
- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：官方下载页链接的 4.0 GPU 版，源码树中包含 `wf_example/`、`cwf_example/` 输入
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站（无许可证，仅供阅读）

关键文献：

- N. Marzari, D. Vanderbilt, Phys. Rev. B 56, 12847 (1997)（MLWF）. DOI: [10.1103/PhysRevB.56.12847](https://doi.org/10.1103/PhysRevB.56.12847)
- I. Souza, N. Marzari, D. Vanderbilt, Phys. Rev. B 65, 035109 (2001)（纠缠能带的解纠缠）. DOI: [10.1103/PhysRevB.65.035109](https://doi.org/10.1103/PhysRevB.65.035109)
- H. Weng, T. Ozaki, K. Terakura, Phys. Rev. B 79, 235118 (2009)（OpenMX 中基于 PAO 的 MLWF 实现；手册要求使用该功能时引用）. DOI: [10.1103/PhysRevB.79.235118](https://doi.org/10.1103/PhysRevB.79.235118)
- T. Ozaki, Phys. Rev. B 110, 125115 (2024)（Closest Wannier functions）. DOI: [10.1103/PhysRevB.110.125115](https://doi.org/10.1103/PhysRevB.110.125115)
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, D. Vanderbilt, Rev. Mod. Phys. 84, 1419 (2012)（MLWF 综述）. DOI: [10.1103/RevModPhys.84.1419](https://doi.org/10.1103/RevModPhys.84.1419)
- G. Pizzi et al., J. Phys.: Condens. Matter 32, 165902 (2020)（Wannier90 v3）. DOI: [10.1088/1361-648X/ab51ff](https://doi.org/10.1088/1361-648X/ab51ff)
- Q. Wu, S. Zhang, H.-F. Song, M. Troyer, A. A. Soluyanov, Comput. Phys. Commun. 224, 405 (2018)（WannierTools）. DOI: [10.1016/j.cpc.2017.09.033](https://doi.org/10.1016/j.cpc.2017.09.033)
- A. H. Castro Neto, F. Guinea, N. M. R. Peres, K. S. Novoselov, A. K. Geim, Rev. Mod. Phys. 81, 109 (2009)（石墨烯紧束缚模型与跳跃积分）. DOI: [10.1103/RevModPhys.81.109](https://doi.org/10.1103/RevModPhys.81.109)

相关章节：[第4章 基组与赝势](04-pao-vps.md)、[第8章 能带](08-band.md)、[第9章 DOS/PDOS](09-dos-pdos-cohp-optics.md)、[第10章 电荷分析](10-charge-mo-cube.md)、[第12章 SOC 与拓扑](12-soc-topology.md)；VASP 中的 Wannier90 用法见 [../../VASP/tutorials/17-advanced-electronic.md](../../VASP/tutorials/17-advanced-electronic.md)；作业提交见 [../../Common/tutorials/02-job-scheduler.md](../../Common/tutorials/02-job-scheduler.md)。

[← 上一章](13-dft-u.md) | [返回目录](../README.md) | [下一章 →](15-berry-efield.md)
