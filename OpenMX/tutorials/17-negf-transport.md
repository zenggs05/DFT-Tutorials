# 第17章 NEGF 电子输运（Electron Transport with NEGF in OpenMX）

本章讲 OpenMX 中基于非平衡格林函数（Non-Equilibrium Green's Function，NEGF）的第一性原理电子输运计算：电极（lead）计算、零偏压和有限偏压下的自洽 NEGF 计算，以及透射谱、电流、电导、本征通道（eigenchannel）和电流密度的后处理。我们以一维碳链和 7-扶手椅型石墨烯纳米带（7-AGNR）为例，完整演示三步法（Step 1/2/3）、只重算 Step 3、I–V 曲线扫描、偏压插值、非共线/SOC 输运，以及如何从透射谱计算 Seebeck 系数等热电量。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 OpenMX NEGF 的器件模型（L∣C∣R，输运方向固定为 **a 轴**），能按“基函数只与相邻胞重叠”的规则确定电极胞和中心区的长度；
- 写出 Step 1（电极能带 + `.hks`）与 Step 2（`scf.EigenvalueSolver NEGF`）的输入，看懂 `<TRAN_Calc_GridBound>` 的区域标签；
- 读懂 Step 3 的 `.tran#_%`、`.current`、`.conductance` 文件以及本征通道、电流密度的输出；
- 用 `NEGF.tran.SCF.skip on` 只重算 Step 3，用 `NEGF.Output.for.TranMain` 处理零偏压周期体系；
- 完成 I–V 曲线计算，了解偏压插值方法，知道非共线/SOC 输运与热电后处理该从哪里入手。

| 场景 | 是否适合 OpenMX NEGF | 说明 |
|---|---|---|
| 分子结、原子链、纳米带、隧道结（Fe\|MgO\|Fe） | 适合 | 典型两电极器件，可做有限偏压自洽 |
| 二维材料面内输运、超晶格 | 适合 | bc 面周期 + k 分辨透射 |
| 栅压调控 | 部分 | `NEGF.gate.voltage` 给出简化的高斯型栅势 |
| 自旋输运、磁性隧道结、畴壁 | 适合 | 共线与非共线都支持，可与 SOC、DFT+U、自旋约束联用 |
| 声子输运、电–声耦合 | 不适合 | OpenMX 不包含；声子相关工作可考虑 alamode 等外部工具 |
| 体材料的扩散输运（迁移率、电导率张量） | 不适合 | 改用 Boltzmann 输运：BoltzTraP/BoltzTraP2 接口（[第9章](09-dos-pdos-cohp-optics.md)） |

示例目录：[17-negf-carbon-chain](../examples/17-negf-carbon-chain/README.md)、[17-negf-7agnr](../examples/17-negf-7agnr/README.md)。前置章节：[第5章 SCF](05-scf-kpoints-convergence.md)、[第8章 能带](08-band.md)、[第10章 cube 可视化](10-charge-mo-cube.md)、[第11章 非共线磁性](11-spin-noncollinear.md)。

---

## 2. 原理速览

### 2.1 器件模型

OpenMX 4.0 处理的体系由左、右两个半无限电极和连接它们的中心区组成，bc 面上为二维周期。对 bc 面做 Bloch 变换后，问题化为一系列依赖 $\mathbf k_\parallel=(k_2,k_3)$ 的一维问题：

```text
 ... | L2 | L1 | L0 |   C0   | R0 | R1 | R2 | ...      (transport along the a-axis)
            \___ extended central region C = L0|C0|R0 ___/
```

为了把电极/中心区界面附近的电子结构弛豫考虑进去，程序自洽求解的是**扩展中心区** $C\equiv L_0|C_0|R_0$ 的格林函数。由于 OpenMX 的 PAO 基函数在实空间严格局域，模型必须满足两个条件（手册第 46.3 节）：

1. C0 中的基函数与 L0、R0 有重叠，但与 L1、R1 **没有**重叠；
2. 电极胞 Li（Ri）中的基函数只与相邻的 L(i±1) 重叠，不与更远的胞重叠。

实际操作中，令 $r_c^\mathrm{max}$ 为体系中最大的 PAO 截断半径（例如 C6.0 对应 6.0 Bohr = 3.175 Å），检查

$$d_{\min}(L_i, L_{i+2}) > 2\,r_c^\mathrm{max},\qquad d_{\min}(C_0, L_1) > 2\,r_c^\mathrm{max}$$

其中 $d_{\min}$ 是两个区域原子之间的最近距离。条件不满足时，就加长电极胞（取 2 倍、3 倍原胞）或者加长 L0/R0。

### 2.2 格林函数、透射与电流

对每个 $\mathbf k_\parallel$ 和自旋 σ，中心区的推迟格林函数为

$$G^r(E)=\left[(E+i\eta)S_C-H_C-\Sigma_L(E)-\Sigma_R(E)\right]^{-1}$$

其中电极自能 $\Sigma_{L/R}$ 由电极表面格林函数得到，展宽矩阵为 $\Gamma_{L/R}=i(\Sigma_{L/R}-\Sigma_{L/R}^\dagger)$。透射系数与 Landauer 电流分别为

$$T_\sigma(E)=\mathrm{Tr}\left[\Gamma_L G^r\Gamma_R G^{r\dagger}\right],\qquad I=\frac{e}{h}\sum_\sigma\int T_\sigma(E,V)\left[f_L(E)-f_R(E)\right]\mathrm dE$$

线性响应（低偏压）下的电导为 $G=\frac{e^2}{h}\sum_\sigma T_\sigma(\mu_L)$。OpenMX 输出的电导以 $G_0=e^2/h$ 为单位，**每个自旋通道单独列出**（`.conductance` 文件头中写明 `G0 = e^2/h`），所以非自旋极化体系的总电导是两列之和。

### 2.3 密度矩阵的计算

自洽计算需要中心区的密度矩阵，它分为两部分：

- **平衡部分**：在复平面上用 Fermi 函数连分式展开（Ozaki 2007）做围道积分，极点数由 `NEGF.Num.Poles` 控制；
- **非平衡部分**（仅有限偏压）：非平衡格林函数在复平面上一般不解析，只能在实轴上加一个小虚部 `NEGF.bias.neq.im.energy`，按步长 `NEGF.bias.neq.energy.step` 做矩形积分，积分范围覆盖偏压窗口及其两侧的 Fermi 尾巴。

总计算量 ∝（极点数 + 实轴能量点数）×（$\mathbf k_\parallel$ 点数）×（自旋数）。手册给出的例子是：碳链在 0.5 V 下实轴积分用 120 个能量点。

### 2.4 OpenMX NEGF 的功能清单（4.0）

零偏压/有限偏压 SCF、栅压、与 DFT+U 联用、自旋相关透射与电流、bc 面 k 分辨透射、I–V 特性、高精度围道积分、偏压插值、三维周期体系的快速零偏压透射、本征通道分析（Paulsson & Brandbyge 方法）、实空间电流与自旋流密度，以及非共线 DFT（含 SOC、约束）下的 NEGF。实现细节见 Ozaki, Nishio & Kino, *PRB* **81**, 035116 (2010)。

---

## 3. 三步法总览

| 步骤 | 程序与输入 | 关键设置 | 产物 |
|---|---|---|---|
| Step 1 电极 | `openmx Lead-*.dat` | 普通 `Band` 计算 + `NEGF.output_hks on` | `*.hks`（电极哈密顿量、电子密度、差分 Hartree 势） |
| Step 2 NEGF SCF | `openmx NEGF-*.dat` | `scf.EigenvalueSolver NEGF`，L0/C0/R0 三个坐标块 | `System.Name.tranb`、`.out` 等 |
| Step 3 输运分析 | 同一次 `openmx` 运行，Step 2 结束后自动进行 | `NEGF.tran.*`、`NEGF.Channel.*` | `.tran#_%`、`.current`、`.conductance`、本征通道 cube、电流密度 cube/xsf |

> [!NOTE]
> OpenMX 早期版本的 Step 3 由单独的程序 `TranMain` 完成。至少从 3.9 起，Step 3 已并入 `openmx`，在 Step 2 之后自动运行，标准输出中会出现 `Welcome to TRAN_Main_Analysis` 字样。仍然沿用的关键词 `NEGF.Output.for.TranMain` 名字里带有历史痕迹（见第 7.2 节）。

安装后先跑一遍官方的 NEGF 自动测试，确认编译正确（手册第 46.11 节）：

```bash
cd ~/openmx4.0/work
mpirun -np 16 ./openmx -runtestNEGF            # or: mpirun -np 8 ./openmx -runtestNEGF -nt 2
cat runtestNEGF.result                         # diff Utot / Force / EigenChannel / CurrentDensity
```

参考结果用 16 个 MPI 进程算出；各项差值在小数点后第 7 位以内一致，即可认为安装正确。

---

## 4. Step 1：电极计算

### 4.1 几何约定

- **输运方向固定为 a 轴（第一条晶格矢量）**。不必旋转结构，只要调整晶格矢量的顺序。例如原本沿 z 的方向要作为输运方向时，把 `0 0 L` 写到第一行即可；但 `scf.Kgrid` 要随之调整，第一个数对应输运方向。官方的 8-ZGNR 示例就是 a=(0, 0, 4.88) Å 的写法。
- 左、右电极在 bc 面内的晶格矢量必须完全相同，否则 Step 2 会以 `Warning: The b- or c-axis of the unit cell for the left lead is not same as that for the right lead.` 退出。中心区的 bc 面晶胞默认与电极相同。
- 电极胞沿 a 轴要足够长，满足 2.1 节的条件 2。

### 4.2 关键词与输入

| 关键词 | 值 | 说明 |
|---|---|---|
| `NEGF.output_hks` | `on` | 输出 Step 2 所需的电极信息（源码读取 `NEGF.Output_HKS`，关键词大小写不敏感） |
| `NEGF.filename.hks` | 例如 `lead-chain.hks` | 输出文件名，默认 `NEGF.hks` |
| `scf.EigenvalueSolver` | `Band` | 普通能带计算 |
| `scf.Kgrid` | `n_a n_b n_c`，$n_a$ 取大 | Step 2 用表面格林函数处理半无限电极，Step 1 沿 a 轴的 k 点要足够密，两步才能自洽一致。官方碳链示例（4.5 Å 胞）用 200，8-ZGNR（4.88 Å）用 201；$n_a\times a\gtrsim 900$ Å 是一个经验起点 |
| `scf.Kgrid` 的 $n_b, n_c$ | 与 Step 2 的 `NEGF.scf.Kgrid` **相同** | 手册明确要求 |
| `scf.criterion` | 1e-9 – 1e-10 | 电极的密度要收敛得很好 |

碳链电极的完整输入 `Lead-Chain.dat`：4 个 C 原子，C–C 1.5 Å，周期 6.0 Å；用 C6.0 基组时满足 $d_{\min}(L_i,L_{i+2})=7.5$ Å $>6.35$ Å。

```text
System.CurrrentDirectory      ./
System.Name                   lead-chain
DATA.PATH                     /path/to/openmx4.0/DFT_DATA19   # EDIT
level.of.stdout               1
level.of.fileout              0

NEGF.output_hks               on
NEGF.filename.hks             lead-chain.hks

Species.Number                1
<Definition.of.Atomic.Species
  C    C6.0-s2p2d1    C_PBE19
Definition.of.Atomic.Species>

Atoms.Number                  4
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
     1  C       0.000000     5.000000     5.000000    2.0 2.0
     2  C       1.500000     5.000000     5.000000    2.0 2.0
     3  C       3.000000     5.000000     5.000000    2.0 2.0
     4  C       4.500000     5.000000     5.000000    2.0 2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit        Ang
<Atoms.UnitVectors
     6.000000     0.000000     0.000000
     0.000000    10.000000     0.000000
     0.000000     0.000000    10.000000
Atoms.UnitVectors>

scf.XcType                    GGA-PBE
scf.SpinPolarization          off
scf.ElectronicTemperature     300.0
scf.energycutoff              220.0
scf.maxIter                   200
scf.EigenvalueSolver          Band
scf.Kgrid                     160 1 1
scf.Mixing.Type               RMM-DIISK
scf.Init.Mixing.Weight        0.05
scf.Min.Mixing.Weight         0.001
scf.Max.Mixing.Weight         0.20
scf.Mixing.History            30
scf.Mixing.StartPulay         10
scf.criterion                 1.0e-10

Band.dispersion               on
Band.Nkpath                   1
<Band.kpath
  60   0.0 0.0 0.0   0.5 0.0 0.0   G X
Band.kpath>

MD.Type                       Nomd
```

```bash
mpirun -np 8 openmx Lead-Chain.dat -nt 1 > lead-chain.std 2>&1
ls -l lead-chain.hks
grep "Chemical potential (Hartree)" lead-chain.out
```

同时输出电极能带（`lead-chain.Band`），可以确认电极的金属性、带隙和通道数：在某个能量 E 处，与 E 相交的能带条数（每个自旋）就是理想透射的上限。

---

## 5. Step 2：NEGF 自洽计算

### 5.1 几何：三个坐标块

| 块 | 内容 | 要求 |
|---|---|---|
| `Atoms.Number` + `<Atoms.SpeciesAndCoordinates` | 中心区 C0 | 满足 2.1 节条件 1 |
| `LeftLeadAtoms.Number` + `<LeftLeadAtoms.SpeciesAndCoordinates` | L0 | 必须是 Step 1 左电极胞的整体平移（原子数、顺序、相对位置相同），b、c 坐标与 Step 1 一致 |
| `RightLeadAtoms.Number` + `<RightLeadAtoms.SpeciesAndCoordinates` | R0 | 同上，对应右电极 |

- `Atoms.SpeciesAndCoordinates.Unit` 只能用 `Ang` 或 `AU`，**不支持 FRAC**；
- 不需要写 `<Atoms.UnitVectors`，半无限电极由 Step 1 的晶胞（存放在 `.hks` 中）自动构造；
- 所有 L0、C0、R0 中出现的物种名，都要在 `<Definition.of.Atomic.Species` 里定义。

程序在标准输出中打印每个原子所属的区域，这是检查几何最有效的办法：

```bash
grep -A12 "<TRAN_Calc_GridBound>" negf-chain.std
```

```text
   Atom1  = 12 Atom2  = 12 Atom3  = 12 Atom4  = 12 Atom5  =  1 ...
   ... Atom13 = 13 Atom14 = 13 Atom15 = 13 Atom16 = 13
```

（示意输出：本例 L0 4 原子、C0 8 原子、R0 4 原子时应得到的标签模式）

标签含义：`12`/`2` 表示属于 L0（`12` 还表示与 L1 有重叠），`13`/`3` 表示属于 R0（`13` 还表示与 R1 有重叠），`1` 表示属于 C0。若 C0 的原子被标成 12 或 13，或者标签模式与你的设计不符，就说明几何不满足重叠条件。

### 5.2 关键词（默认值按 4.0 源码 `TRAN_Input_std.c` 核对）

| 关键词 | 默认 | 推荐 / 说明 |
|---|---|---|
| `scf.EigenvalueSolver` | — | `NEGF` |
| `NEGF.filename.hks.l` / `.r` | `NEGF.hks.l` / `NEGF.hks.r` | 左/右电极的 `.hks`；对称结可以用同一个文件 |
| `NEGF.Num.Poles` | 150 | 平衡部分的围道积分极点数；手册示例取 100。低电子温度或宽能谱体系应作收敛测试 |
| `NEGF.scf.Kgrid` | 1 1 | b、c 方向的 k 点，必须与 Step 1 的 $n_b, n_c$ 相同 |
| `NEGF.SCF.Iter.Band` | **源码默认 10**（手册写 6） | 前若干步沿 a 轴也按周期体系做普通对角化，可以显著改善 NEGF 的 SCF 收敛；难收敛体系可加到 20–50（官方碳链示例用 50） |
| `NEGF.Poisson.Solver` | `FD` | `FD`：2D FFT + 1D 有限差分，严格满足电极/中心区界面的静电边界条件，理论上更合理；`FFT`：3D FFT，SCF 通常更容易收敛。非极性体系两者结果接近，极性体系差别可能很大 |
| `NEGF.bias.voltage` | 0.0 | 源漏偏压（eV，数值上等于 V） |
| `NEGF.bias.neq.im.energy` | 0.01 eV | 非平衡部分实轴积分的虚部 |
| `NEGF.bias.neq.energy.step` | 0.02 eV | 实轴积分步长；计算时间与能量点数成正比 |
| `NEGF.gate.voltage` | 0.0 | 栅压 $V_g^{(0)}$（eV），见 5.5 节 |
| `scf.criterion` | 1e-6 | NEGF 中作用于 **NormRD**（其他求解器作用于 dUele）。建议 1e-7 |

### 5.3 完整输入（零偏压碳链）

```text
System.CurrrentDirectory      ./
System.Name                   negf-chain
DATA.PATH                     /path/to/openmx4.0/DFT_DATA19   # EDIT
level.of.stdout               1
level.of.fileout              1

Species.Number                1
<Definition.of.Atomic.Species
  C    C6.0-s2p2d1    C_PBE19
Definition.of.Atomic.Species>

Atoms.Number                  8
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
     1  C       6.000000     5.000000     5.000000    2.0 2.0
     2  C       7.500000     5.000000     5.000000    2.0 2.0
     3  C       9.000000     5.000000     5.000000    2.0 2.0
     4  C      10.500000     5.000000     5.000000    2.0 2.0
     5  C      12.000000     5.000000     5.000000    2.0 2.0
     6  C      13.500000     5.000000     5.000000    2.0 2.0
     7  C      15.000000     5.000000     5.000000    2.0 2.0
     8  C      16.500000     5.000000     5.000000    2.0 2.0
Atoms.SpeciesAndCoordinates>

LeftLeadAtoms.Number          4
<LeftLeadAtoms.SpeciesAndCoordinates
     1  C       0.000000     5.000000     5.000000    2.0 2.0
     2  C       1.500000     5.000000     5.000000    2.0 2.0
     3  C       3.000000     5.000000     5.000000    2.0 2.0
     4  C       4.500000     5.000000     5.000000    2.0 2.0
LeftLeadAtoms.SpeciesAndCoordinates>

RightLeadAtoms.Number         4
<RightLeadAtoms.SpeciesAndCoordinates
     1  C      18.000000     5.000000     5.000000    2.0 2.0
     2  C      19.500000     5.000000     5.000000    2.0 2.0
     3  C      21.000000     5.000000     5.000000    2.0 2.0
     4  C      22.500000     5.000000     5.000000    2.0 2.0
RightLeadAtoms.SpeciesAndCoordinates>

scf.XcType                    GGA-PBE
scf.SpinPolarization          off
scf.ElectronicTemperature     300.0
scf.energycutoff              220.0
scf.maxIter                   300
scf.EigenvalueSolver          NEGF
scf.Mixing.Type               RMM-DIISK
scf.Init.Mixing.Weight        0.02
scf.Min.Mixing.Weight         0.001
scf.Max.Mixing.Weight         0.05
scf.Mixing.History            30
scf.Mixing.StartPulay         20
scf.criterion                 1.0e-7

NEGF.filename.hks.l           lead-chain.hks
NEGF.filename.hks.r           lead-chain.hks
NEGF.Num.Poles                150
NEGF.scf.Kgrid                1 1
NEGF.SCF.Iter.Band            10
NEGF.Poisson.Solver           FD
NEGF.bias.voltage             0.0
NEGF.bias.neq.im.energy       0.01
NEGF.bias.neq.energy.step     0.02

NEGF.tran.Analysis            on
NEGF.tran.energyrange         -5.0 5.0 1.0e-3
NEGF.tran.energydiv           1001
NEGF.tran.Kgrid               1 1
NEGF.tran.CurrentDensity      on
NEGF.tran.Channel             on
NEGF.Channel.Nkpoint          1
<NEGF.Channel.kpoint
  0.0  0.0
NEGF.Channel.kpoint>
NEGF.Channel.Nenergy          1
<NEGF.Channel.energy
  0.0
NEGF.Channel.energy>
NEGF.Channel.Num              4

MD.Type                       Nomd
```

> [!IMPORTANT]
> Step 1 与 Step 2 的 `scf.XcType`、`scf.SpinPolarization`、`scf.ElectronicTemperature`、`scf.energycutoff`、基组/赝势必须一致。源码只检查自旋设置（不一致时报 `scf.SpinPolarization conflicts between leads or lead and center.`），其余几项不一致**不会报错**，但电极与中心区的势无法衔接，会在界面处产生虚假散射。

### 5.4 运行与 SCF 监控

```bash
mpirun -np 8 openmx NEGF-Chain.dat -nt 1 > negf-chain.std 2>&1
grep -A2 "Intrinsic chemical potential" negf-chain.std
grep "NormRD" negf-chain.out | tail -n 3
grep "normally finished" negf-chain.std
```

NEGF 的 SCF 通常比普通计算难收敛。按以下顺序调整：

1. 加大 `NEGF.SCF.Iter.Band`（10 → 30 → 50）；
2. 降低混合权重（`scf.Max.Mixing.Weight` 0.02–0.05），`scf.Mixing.StartPulay` 设为 20–30；
3. 换成 `NEGF.Poisson.Solver FFT`。对非极性体系影响很小，但要在报告中写明；
4. 加长 L0/R0，让界面离 C0 更远；
5. 先在零偏压下收敛，再以零偏压结果为初猜计算有限偏压（第 8 节）。

### 5.5 偏压与栅压的约定

- **偏压的施加方式（4.0）**：读源码 `TRAN_Input_std.c` 可知，偏压对称地加在两侧，$\mu_L=\mu_0-V/2$、$\mu_R=\mu_0+V/2$。官方示例 `negf-chain.out`（0.1 V）中两侧化学势分别为 −0.190201 和 −0.186526 Hartree，相差正好 0.1 eV，并关于电极本征化学势 −0.188364 Hartree 对称。手册第 46.3 节示例输出里“只给右电极加 0.5 eV”的写法出自旧版。核对偏压的方法：

  ```bash
  grep "Chemical potential of left lead\|Chemical potential of right lead" negf-chain-0.5.out
  ```

- **透射谱的能量零点**始终是左电极化学势 $\mu_L$（源码注释 *The zero energy is the chemical potential of the left lead*）。有限偏压下，偏压窗口在 tran 文件的能量轴上对应 $[0,\,\mu_R-\mu_L]=[0,V]$。
- **栅压**：$V_g(x)=V_g^{(0)}\exp\left[-\left(\frac{x-x_c}{d}\right)^8\right]$，其中 $x_c$ 是 C0 的中心，$d$ 是 C0 沿 a 轴的长度。这个势主要作用在 C0 上，形状类似镜像电荷产生的势。设置方式是 `NEGF.gate.voltage 1.0`（eV）。

### 5.6 NEGF 中的 DOS

```text
Dos.fileout                   on
NEGF.Dos.energyrange          -3.0 3.0 5.0e-3   # lower, upper (eV), imaginary smearing
NEGF.Dos.energy.div           601
NEGF.Dos.Kgrid                1 1               # may be denser than NEGF.scf.Kgrid
```

运行后生成 `System.Name.Dos.val/.Dos.vec`，用 DosMain 处理的方法与普通计算相同（[第9章](09-dos-pdos-cohp-optics.md)），能量零点同样是左电极化学势。

---

## 6. Step 3：透射、电流、电导、本征通道与电流密度

### 6.1 关键词

| 关键词 | 默认 | 说明 |
|---|---|---|
| `NEGF.tran.Analysis` | on | 计算透射、电流和电导 |
| `NEGF.tran.energyrange` | −10 10 1e-3 | 下限、上限（eV，相对 $\mu_L$）、虚部（平滑用） |
| `NEGF.tran.energydiv` | 200 | 能量点数，例如 −5…5 eV 取 1001 点，间隔 0.01 eV |
| `NEGF.tran.Kgrid` | 等于 `NEGF.scf.Kgrid` | b、c 方向的 k 点，可以比 SCF 用的更密，建议更密 |
| `NEGF.tran.CurrentDensity` | on | 实空间电流（自旋流）密度。打开它会自动打开 `NEGF.tran.Analysis` |
| `NEGF.tran.Channel` | on | 本征通道。打开它会自动打开 `NEGF.tran.Analysis` |
| `NEGF.Channel.Nkpoint` + `<NEGF.Channel.kpoint` | 1 / (0,0) | 计算本征通道的 $(k_2,k_3)$，以 b、c 倒格矢为单位 |
| `NEGF.Channel.Nenergy` + `<NEGF.Channel.energy` | 1 / 0.0 | 计算本征通道的能量（eV，相对 $\mu_L$） |
| `NEGF.Channel.Num` | 5（共线），10（非共线） | 每个 (k, E, σ) 输出的实空间通道数，按透射本征值从大到小 |
| `NEGF.OffDiagonalCurrent` | off | 非共线体系的非对角自旋流分量（4.0） |

> [!TIP]
> 本征通道和电流密度默认都是打开的，每次运行都会写出一批 cube/xsf 文件。做 I–V 扫描或收敛测试时，建议显式关闭（`NEGF.tran.Channel off`、`NEGF.tran.CurrentDensity off`），`run_iv.sh` 就是这样处理的。需要分析时再用 `NEGF.tran.SCF.skip on` 单独重算 Step 3。

### 6.2 输出文件与列定义

| 文件 | 内容 |
|---|---|
| `System.Name.tran#_%` | 每个 $(k_2,k_3)$ 点一个文件（`#`、`%` 为 k 网格序号），9 列：`iw`、Re E (Hartree)、Im E (Hartree)、**Re E (eV)**、Im E (eV)、**Re T↑**、Im T↑、**Re T↓**、Im T↓。画透射谱用第 4 列对第 6（8）列 |
| `System.Name.current` | k 分辨电流与平均电流（A）；文件头含 `# average current for up and down spins (ampere)=` |
| `System.Name.conductance` | k 分辨电导与平均电导（$G_0=e^2/h$）；文件头含 `# average conductance for up and down spins (G0)=` |
| `System.Name.traneval#k_#E_#s` / `tranevec…` | 本征通道的透射本征值和 LCAO 系数 |
| `System.Name.tranec#k_#E_#s_#c_r.cube` / `_i.cube` | 本征通道波函数的实部和虚部（VESTA/XCrySDen 可视化） |
| `System.Name.curden1.cube`、`scurden1.cube` | 电流（自旋流）密度的 a 分量 |
| `System.Name.curden.xsf`、`scurden.xsf` | 网格上的电流矢量场，在 XCrySDen 中用 Display→Forces 显示 |
| `System.Name.curden_atom.xsf` | 按 Voronoi 多面体平均到各原子的电流 |

一次运行的快速检查：

```bash
grep "average conductance" negf-chain.conductance
grep "average current" negf-chain-0.5.current
awk '!/^#/ && NF>=8 && $4>-0.005 && $4<0.005 {print "T_up(E~0) =", $6}' negf-chain.tran0_0
```

### 6.3 作图与平均

示例里的 [plot_tran.py](../examples/17-negf-carbon-chain/plot_tran.py) 对多个 k 文件做等权平均，画出透射谱，并给出零偏压电导和 Seebeck 系数。OpenMX 部分的公共脚本 omx_negf_plot.py 还可以画 I–V 曲线和 k 分辨透射图，用法见 `python ../scripts/omx_negf_plot.py -h`。

```bash
python3 plot_tran.py negf-chain.tran0_0 --plot negf-chain_T.png
python3 plot_tran.py "negf-graphene.tran*_*" --csv T_avg.csv      # many k files
```

**结果判断（理想碳链）**：电极与中心区完全相同时，$T(E)$ 必须等于电极在该能量上的能带条数（每个自旋）。碳链在 $E=0$ 附近有两条简并的 π 带，所以 $T\approx2$，`.conductance` 给出约 2、2（单位 $e^2/h$）。偏离整数通常说明：几何不满足重叠条件、Step 1 的 a 方向 k 点不够、Step 1 与 Step 2 的设置不一致，或者 `NEGF.tran.energyrange` 的虚部太大导致台阶被抹平。

### 6.4 本征通道

本征通道把总透射分解为互不干涉的独立通道，$T=\sum_n\tau_n$，适合分析“电流走哪条路”。`traneval` 文件中的本征值之和应等于该能量处的 $T$。可视化时，对 `*_r.cube` 和 `*_i.cube` 分别画等值面（或者画 $|\psi|^2$）。散射体附近振幅突然衰减的地方，就是发生背散射的位置。

---

## 7. 只重算 Step 3 与零偏压周期体系

### 7.1 `NEGF.tran.SCF.skip`

Step 2 是最耗时的部分。如果只想换一个能量窗口、加密 k 点或增加本征通道能量，可以只重算 Step 3：

```text
NEGF.tran.SCF.skip            on          # read System.Name.tranb, skip Step 2
```

`System.Name` 必须与 Step 2 的计算相同，因为程序会读取 `System.Name.tranb`。示例中的 [NEGF-Chain_step3.dat](../examples/17-negf-carbon-chain/NEGF-Chain_step3.dat) 把能量网格加密到 −2…2 eV、2001 点，虚部减小到 1e-4，并在 −1、0、1 eV 三个能量处计算本征通道。

反过来，若只想做完 Step 2 就停止：

```text
NEGF.tran.SCF.skip            off
NEGF.tran.Analysis            off
NEGF.tran.Channel             off
```

### 7.2 三维周期体系的零偏压透射：`NEGF.Output.for.TranMain`

如果体系沿 a 轴本身就是周期的（超晶格、多层膜、周期磁结构），零偏压透射可以直接用普通能带计算的哈密顿量求出，**省去 Step 2**（手册第 46.6 节）：

```text
NEGF.Output.for.TranMain      on          # in the Step 1-type band calculation
```

计算正常结束后会生成可供 Step 3 使用的 `*.tranb`，用来快速扫描大量几何或磁构型的透射。

> [!NOTE]
> 4.0 源码中 `NEGF.tran.SCF.skip` 还接受第三个值 `Periodic`，另有一个 `NEGF.tranb.file` 关键词（默认等于 System.Name），二者都没有写进手册。周期体系的 Step 3 具体怎么组织输入，请对照手册第 46.6 节和 `work/negf_example`，本教程不展开未记载的用法。

---

## 8. 有限偏压与 I–V 曲线

### 8.1 逐个偏压做 SCF

非线性 I–V 需要在每个偏压下分别做 NEGF 自洽。示例脚本 [run_iv.sh](../examples/17-negf-carbon-chain/run_iv.sh) 为每个偏压生成一个输入（`System.Name negf-chain-V<V>`），运行后把电流汇总到 `iv.csv`：

```bash
cd OpenMX/examples/17-negf-carbon-chain
NP=8 bash run_iv.sh                                  # 0.0 0.2 ... 1.0 V
BIASES="0.0 0.1 0.2 0.3" NP=8 bash run_iv.sh         # custom list
RESTART=1 NP=8 bash run_iv.sh                        # start from previous bias
cat iv.csv
```

**加速技巧**：用前一个偏压收敛的密度作初猜。OpenMX 的 `scf.restart on` 默认读取同名 `System.Name_rst/` 中的文件；4.0 源码读取重启文件时使用 `scf.restart.filename`（默认等于 System.Name，手册只在 `c2n` 场景下提到它），因此可以写：

```text
scf.restart                   on
scf.restart.filename          negf-chain-V0.4     # System.Name of the previous bias
```

> [!WARNING]
> 重启目录中的密度文件 `*.crst<spin>_<rank>_<n>` 按 MPI 进程分块保存，程序只校验网格数是否一致。不同偏压之间请保持**相同的 MPI 进程数和相同的网格**（`scf.energycutoff`、几何都不变），否则读入的密度是错乱的。这一用法请先在小体系上验证。

偏压越大，实轴积分的能量点越多（约为 (V + 若干 $k_BT$)/步长），计算时间随之增加。标准输出会打印实际使用的点数：

```bash
grep -A4 "Parameters for the integration of the non-equilibrium part" negf-chain-V1.0.std
```

Step 3 计算电流时，若 `NEGF.current.energy.step`（默认 0.01 eV，按 4.0 源码）对当前偏压窗口偏大，程序会在标准输出中给出提示和推荐值，按提示调整即可。

### 8.2 偏压插值

大体系对每个偏压都做 SCF 太贵，可以用插值方法（手册第 46.7 节）：先在两三个偏压下做 SCF，然后在 Step 3 中对中心区和右电极的哈密顿量块以及 $\mu_R$ 做线性插值：

$$H_{\sigma,C}(\mathbf k)=\lambda H_{\sigma,C}(\mathbf k,1)+(1-\lambda)H_{\sigma,C}(\mathbf k,2),\qquad \mu_R=\lambda\mu_R^{(1)}+(1-\lambda)\mu_R^{(2)}$$

```text
NEGF.tran.SCF.skip            on
NEGF.tran.interpolate         on                   # default off
NEGF.tran.interpolate.file1   negf-chain-V0.4.tranb
NEGF.tran.interpolate.file2   negf-chain-V0.6.tranb
NEGF.tran.interpolate.coes    0.5 0.5              # V = 0.5*0.4 + 0.5*0.6 = 0.5 V
```

λ 应取在 0 与 1 之间，只做内插。手册用碳链做了对比：用 0、0.5、1.0 V 三个 SCF 结果插值，得到的电流和透射与完整 SCF 非常接近。

---

## 9. 实战二：7-AGNR，纯净 vs N 掺杂

示例目录：[17-negf-7agnr](../examples/17-negf-7agnr/README.md)。要点：

| 项目 | 设置 | 理由 |
|---|---|---|
| 电极胞 | 2 个周期，8.52 Å，36 原子 | 1 个周期（4.26 Å）时 $d_{\min}(L_i,L_{i+2})=5.68$ Å < 6.35 Å，不满足条件；2 个周期时为 9.94 Å |
| L0/C0/R0 | 各 2 个周期 | C0 到 L1 的最近距离同样 > 6.35 Å |
| 基组 | `C6.0-s2p2`、`H5.0-s2`、`N6.0-s2p2`（Quick） | NEGF 计算量随轨道数增长很快，教学示例取 Quick |
| k 点 | Step 1 `100 1 1`，Step 2 `NEGF.scf.Kgrid 1 1` | 孤立纳米带，b、c 方向只需 1 个 k 点 |
| 掺杂 | 中心区第 12 号原子 C→N（坐标行 2.5 2.5） | 取代式施主 |
| 输出 | 透射、E=±1.2 eV 处的本征通道、NEGF DOS | 对比纯净与掺杂 |

```bash
cd OpenMX/examples/17-negf-7agnr
mpirun -np 32 openmx Lead-7AGNR.dat  -nt 1 > lead-7agnr.std 2>&1
mpirun -np 32 openmx NEGF-7AGNR.dat  -nt 1 > negf-7agnr.std 2>&1
mpirun -np 32 openmx NEGF-7AGNR-N.dat -nt 1 > negf-7agnr-n.std 2>&1
python3 ../17-negf-carbon-chain/plot_tran.py negf-7agnr.tran0_0   --csv T_pristine.csv
python3 ../17-negf-carbon-chain/plot_tran.py negf-7agnr-n.tran0_0 --csv T_Ndoped.csv
```

**结果判断**：纯净纳米带的 $T(E)$ 是整数台阶，费米能级附近有一段 $T=0$ 的带隙区，台阶位置与 `lead-7agnr.Band` 的子带边一一对应；N 掺杂后，导带一侧的台阶出现凹陷（背散射），E=+1.2 eV 处本征通道的振幅在 N 附近被削弱。

---

## 10. 非共线与 SOC 输运

OpenMX 4.0 的 NEGF 与非共线 DFT 完全兼容，包括 SOC、DFT+U 以及对自旋/轨道磁矩方向的约束（手册第 46.9 节）。操作流程与共线情况相同：

- Step 1 和 Step 2 都设 `scf.SpinPolarization nc`，需要时加上 `scf.SpinOrbit.Coupling on`；
- 左右电极可以有不同的磁化方向，例如官方示例 8-ZGNR 的左电极自旋朝上、右电极朝右。电极在 Step 1 中用自旋约束固定磁化方向，Step 2 中不施加约束，中心区的自旋会自洽地逐渐转向；
- `NEGF.Channel.Num` 的默认值变为 10；`NEGF.OffDiagonalCurrent on` 输出非对角自旋流（实部、虚部分开存放）；
- 官方输入在 `work/negf_example`：`Lead-L-8ZGNR-NC.dat`、`Lead-R-8ZGNR-NC.dat`、`NEGF-8ZGNR-NC.dat`，以及 Au 原子链 `Lead-Au-Chain-NC.dat`、`NEGF-Au-Chain-NC.dat`。

共线自旋极化的磁隧道结（例如 Fe|MgO|Fe 的平行/反平行构型），用两组 Step 1（左右电极的磁矩取相同或相反的初始占据）配合 `scf.SpinPolarization on` 即可。隧穿磁阻 TMR = $(G_P-G_{AP})/G_{AP}$，k 分辨电导图需要非常密的 `NEGF.tran.Kgrid`；手册图 46 中的 Fe|MgO|Fe 用了 120×120 个 k 点。初始磁矩的写法见[第11章](11-spin-noncollinear.md)。

---

## 11. 热电：从透射谱到 Seebeck 系数与 ZT

在线性响应、弹道输运的框架下，由零偏压透射谱 $T(E)$（自旋求和）可以得到电子的热电输运系数：

$$L_n=\int T(E)\,(E-\mu)^n\left(-\frac{\partial f}{\partial E}\right)\mathrm dE$$

$$G=\frac{e^2}{h}L_0,\qquad S=-\frac{1}{eT}\frac{L_1}{L_0},\qquad \kappa_e=\frac{1}{hT}\left(L_2-\frac{L_1^2}{L_0}\right),\qquad ZT=\frac{S^2GT}{\kappa_e+\kappa_\mathrm{ph}}$$

`plot_tran.py --mu ... --temp 300` 可以对一组化学势计算 $G$ 和 $S$。功率因子 $S^2G$ 和 $\kappa_e$ 的计算方法相同，自己补几行即可。声子热导 $\kappa_\mathrm{ph}$ 需要另外计算（OpenMX 本身没有声子模块）。思路可参考 [TnakaTnakaTnaka/TEP_tools](https://github.com/TnakaTnakaTnaka/TEP_tools)，它从 NEGF 透射谱出发，计算 BZ 平均透射、电导、Seebeck 系数、功率因子和 ZT（无许可证，只作思路参考，不要复制代码）。

> [!NOTE]
> 半导体纳米带的 Seebeck 系数对 μ 的位置非常敏感。不要只报告 $\mu=\mu_L$ 一个点，应当画出 $S(\mu)$ 曲线，并说明对应的掺杂水平。对于体材料的扩散输运，请改用 Boltzmann 方法（BoltzTraP/BoltzTraP2 接口）。

---

## 12. 并行与计算量

- MPI 并行分配在 (能量点 × k 点 × 自旋) 的循环上。手册中 8-ZGNR 的 0.5 V 计算共有 197 个能量点、只用 Γ 点、自旋极化，共 394 个组合可以分配；4.0 中 MPI 进程数超过原子数时也能高效并行；
- 链接 MKL 时，格林函数中的矩阵乘法和求逆还可以用 OpenMP 并行，即 `-nt` 混合并行。手册的基准测试（CRAY-XT5）中，64 进程 × 4 线程的加速比约为 140；
- 估算计算量：围道积分点数（`NEGF.Num.Poles`）+ 实轴点数（偏压窗口 / 步长），再乘以 $N_{k_\parallel}$ 和自旋数，每个点需要对扩展中心区的矩阵求逆一次（维度 = 基函数总数）。

```bash
mpirun -np 64 openmx NEGF-7AGNR.dat -nt 2 > negf-7agnr.std 2>&1
grep "Total Computational Time\|Diagonalization" negf-7agnr.out | head
```

更多并行、作业脚本和排错内容见第18章。

---

## 参数选择与收敛建议

| 参数 | 起点 | 收敛判据 |
|---|---|---|
| 电极胞长度 | $d_{\min}(L_i,L_{i+2})>2r_c^\mathrm{max}$ | 硬性条件 |
| L0/R0、C0 长度 | 至少一个电极胞；界面附近留出缓冲 | 加长后 $T(E_F)$ 变化 < 1% |
| Step 1 `scf.Kgrid` 的 $n_a$ | $n_a\times a\gtrsim900$ Å | 理想体系的 $T$ 为整数 |
| `NEGF.scf.Kgrid` | 孤立体系 1 1；二维体系与 Step 1 的 $n_b,n_c$ 相同 | 与 Step 1 一致是硬性条件 |
| `NEGF.tran.Kgrid` | ≥ `NEGF.scf.Kgrid`，二维体系往往要几十×几十 | 平均电导变化 < 1% |
| `NEGF.Num.Poles` | 100–150 | 总能量 / NormRD 不随极点数变化 |
| `NEGF.bias.neq.energy.step` | 0.02 eV（默认） | 电流变化 < 1% |
| `NEGF.tran.energyrange` 虚部 | 1e-3 eV；精细谱用 1e-4 | 台阶不被抹平 |
| `NEGF.SCF.Iter.Band` | 10（难收敛时 20–50） | — |
| `NEGF.Poisson.Solver` | FD；不收敛时改 FFT | 两者结果对比（极性体系尤其要比） |
| `scf.criterion` | 1e-7（作用于 NormRD） | — |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `Warning: The b- or c-axis of the unit cell for the left lead is not same as that for the right lead.` | 左右电极的 bc 面晶格矢量不同 | 两个 Step 1 输入使用完全相同的 b、c |
| `scf.SpinPolarization conflicts between leads or lead and center.` | Step 1 与 Step 2 的自旋设置不一致 | 统一为 off/on/nc |
| `LeftLeadAtoms.Number may be wrong.` / `RightLeadAtoms.Number may be wrong.` | 缺少电极块，或原子数 ≤ 0 | 补全 `LeftLeadAtoms.*` 和 `RightLeadAtoms.*` |
| `NEGF.scf.Kgrid should be over 1` / `NEGF.tran.Kgrid should be over 1` | 写了 0 或负数 | 至少写 1 1 |
| `NEGF.bias.neq.im.energy should be positive.` | 虚部为负 | 改为正值（默认 0.01） |
| 找不到 `.hks` 文件、程序在读取电极时退出 | `NEGF.filename.hks.l/r` 写错，或 Step 1 没有在当前目录运行 | 检查文件名与路径 |
| 理想体系的 $T(E)$ 不是整数 | Step 1 沿 a 的 k 点不够、几何不满足重叠条件、Step 1/2 设置不一致 | 加密 $n_a$；检查 `<TRAN_Calc_GridBound>` 标签；统一 XcType/温度/截断能 |
| `<TRAN_Calc_GridBound>` 中 C0 原子被标成 12/13 | C0 太短，与 L1/R1 有重叠 | 加长 L0/R0 或 C0 |
| NEGF 的 SCF 振荡、不收敛 | 界面电荷重排、极性体系 | 见 5.4 节：加大 Iter.Band、降低权重、改用 FFT、加长 L0/R0 |
| 有限偏压下电流噪声大 | 实轴步长太大，或 Step 3 电流积分步长偏大 | 减小 `NEGF.bias.neq.energy.step`；按标准输出的提示调整 `NEGF.current.energy.step` |
| 只重算 Step 3 时找不到 `.tranb` | `System.Name` 与 Step 2 不同 | 保持同名 |
| 用 FRAC 坐标后结构错乱 | NEGF 不支持 FRAC | 改用 Ang/AU |

---

## 练习

1. 把碳链电极改成 3 原子、4.5 Å 的胞，其余不变，观察 `<TRAN_Calc_GridBound>` 的输出和 $T(E)$ 发生了什么变化；然后把基组换成 C5.0-s2p2d1 再试一次。
2. 在碳链 C0 中间把一根键拉长到 2.0 Å，计算 $T(E)$ 和 0–1 V 的 I–V 曲线，比较与理想链的差别，并用本征通道定位散射。
3. 用 0、0.5、1.0 V 三个 SCF 结果，以插值方法计算 0.25 V 和 0.75 V 的电流，与直接 SCF 的结果比较误差。
4. 对 7-AGNR-N 做 `NEGF.Poisson.Solver FD` 与 `FFT` 的对比，看看透射谱的差别。
5. 由 7-AGNR 的透射谱计算 300 K 下 μ 在 −1.5…+1.5 eV 范围内的 $S(\mu)$ 曲线，找出 $|S|$ 的最大值及其位置。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 46 章 NEGF：<https://www.openmx-square.org/openmx_man4.0/s46_negf.html>
  - 46.1 概述：<https://www.openmx-square.org/openmx_man4.0/s46_1_overview.html>
  - 46.2 Step 1：<https://www.openmx-square.org/openmx_man4.0/s46_2_step1-leads.html>
  - 46.3 Step 2：<https://www.openmx-square.org/openmx_man4.0/s46_3_step2-negf.html>
  - 46.4 Step 3：<https://www.openmx-square.org/openmx_man4.0/s46_4_step3-trans.html>
  - 46.5 只重算 Step 3：<https://www.openmx-square.org/openmx_man4.0/s46_5_recalc-step3.html>
  - 46.6 零偏压周期体系：<https://www.openmx-square.org/openmx_man4.0/s46_6_zero-bias.html>
  - 46.7 偏压插值：<https://www.openmx-square.org/openmx_man4.0/s46_7_bias-interp.html>
  - 46.8 并行：<https://www.openmx-square.org/openmx_man4.0/s46_8_parallel.html>
  - 46.9 非共线 NEGF：<https://www.openmx-square.org/openmx_man4.0/s46_9_ncl.html>
  - 46.10 示例与 46.11 自动测试：<https://www.openmx-square.org/openmx_man4.0/s46_10_examples.html>
- 第 17 章 重启（`scf.restart`、`scf.restart.filename`）：<https://www.openmx-square.org/openmx_man4.0/s17_restart.html>
- 第 63 章 BoltzTraP/BoltzTraP2 接口（体材料输运）：<https://www.openmx-square.org/openmx_man4.0/s63_boltztrap.html>
- OpenMX 论坛：<https://www.openmx-square.org/forum/patio.cgi>
- Ozaki 讲座幻灯片（NEGF，ISSP 2018）：<https://t-ozaki.issp.u-tokyo.ac.jp/ISS18/slides-1st/7-NEGF-Ozaki.pdf>

**相关 GitHub 仓库（只链接，不复制）**

- [pulkin/openmx-hks](https://github.com/pulkin/openmx-hks)：读取 NEGF Step 1 生成的 `.hks`，把哈密顿量/重叠矩阵导出为 hdf5/json/MATLAB 格式，附能带与弹道输运示例
- [TnakaTnakaTnaka/TEP_tools](https://github.com/TnakaTnakaTnaka/TEP_tools)：由 NEGF 透射谱计算 BZ 平均透射、电导、Seebeck 系数、功率因子与 ZT（热电后处理思路）
- [samlipton/pyfesm](https://github.com/samlipton/pyfesm)：轻量 Python 解析器，可读取 NEGF 透射输出，并处理单位换算
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：从 OpenMX 哈密顿量构建紧束缚模型，做 Berry 曲率、光学和非线性响应计算，可与输运结果互补
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，基础计算的对照

**关键文献**

- T. Ozaki, K. Nishio, H. Kino, *Phys. Rev. B* **81**, 035116 (2010)，OpenMX 的 NEGF 实现：<https://doi.org/10.1103/PhysRevB.81.035116>
- T. Ozaki, *Phys. Rev. B* **75**, 035123 (2007)，Fermi 函数连分式展开与围道积分：<https://doi.org/10.1103/PhysRevB.75.035123>
- M. Brandbyge, J.-L. Mozos, P. Ordejón, J. Taylor, K. Stokbro, *Phys. Rev. B* **65**, 165401 (2002)，TranSIESTA 的 NEGF-DFT 方法：<https://doi.org/10.1103/PhysRevB.65.165401>
- M. Paulsson, M. Brandbyge, *Phys. Rev. B* **76**, 115117 (2007)，透射本征通道：<https://doi.org/10.1103/PhysRevB.76.115117>
- Y.-W. Son, M. L. Cohen, S. G. Louie, *Phys. Rev. Lett.* **97**, 216803 (2006)，石墨烯纳米带的能隙：<https://doi.org/10.1103/PhysRevLett.97.216803>

---

[← 上一章](16-surface-esm-workfunction.md) | [返回目录](../README.md)
