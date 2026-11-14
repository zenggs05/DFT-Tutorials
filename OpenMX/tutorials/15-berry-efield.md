# 第15章 Berry 相宏观极化与外加电场（Berry-Phase Polarization & External Electric Fields in OpenMX）

本章讲 OpenMX 中与“电极化”和“电场”有关的两类功能。第一类是零场下的 **Berry 相宏观极化**：用 `HS.fileout on` 生成 scfout 文件，交给后处理程序 `polB` 计算电偶极（Debye）和极化（μC/cm²），在此基础上用有限位移法求 **Born 有效电荷 Z\***，并以 BaTiO3 为例求铁电体中心对称结构与极化结构之间的自发极化差 ΔP。第二类是**有限外电场**：用 `scf.Electric.Field` 施加锯齿势表示的均匀电场（单位 GV/m），以硝基苯为例求偶极响应、静态极化率和 `diff_gcube` 差分密度，同时说明体相无真空时锯齿势不连续面的问题，以及何时应改用 ESM on4（OpenMX Berry phase polarization, polB, Born effective charge, ferroelectric BaTiO3, sawtooth electric field, polarizability tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清楚为什么周期性晶体的极化不能用“电荷 × 位置”直接积分，以及 Berry 相公式给出的极化为什么是**多值**的（极化量子与分支）；
- 完成 OpenMX 的 polB 工作流：`HS.fileout on` → `System.Name.scfout` → `make polB` → `polB x.scfout < in > out`，并读懂输出中 Background / Core / Electron / Total 四列的含义；
- 用 ±u 位移的中心差分求 Born 有效电荷张量，并用芯电荷/电子分解和声学求和规则检查结果；
- 用“同一晶胞、λ 路径插值”的办法求铁电体自发极化，避开极化量子带来的分支错误；
- 正确使用 `scf.Electric.Field`：单位、方向（沿 a、b、c 晶格矢量）、符号约定、不连续面的位置，以及在分子和板模型中的适用范围；
- 用有限差分求分子静态极化率，用 `diff_gcube` 画场致电荷密度差；
- 知道 OpenMX 4.0 中**没有**体相 Berry 相有限电场方法，需要体相介电响应时有哪些替代路线。

| 任务 | OpenMX 中的做法 | 本章示例 |
|---|---|---|
| 体相宏观极化 P | SCF + `HS.fileout on`，后处理 `polB` | [15-nacl-polb-born](../examples/15-nacl-polb-born/README.md) |
| Born 有效电荷 Z\* | 原子 ±u 位移，各算一次 P，中心差分 | 同上，`born_charge.py` |
| 铁电自发极化 ΔP | 中心对称 → 极化结构的 λ 路径，逐点 polB，按量子展开 | [15-batio3-polarization](../examples/15-batio3-polarization/README.md) |
| 分子在均匀电场中 | `scf.Electric.Field`（锯齿势） | [15-nitrobenzene-efield](../examples/15-nitrobenzene-efield/README.md) |
| 场致电荷重排 | 两个 `tden.cube` 相减：`diff_gcube` | 同上 |
| 板模型 / 电容器中的均匀场 | ESM `on4`（或有真空的 slab + 锯齿势） | [第16章](16-surface-esm-workfunction.md) |

前置知识：输入文件结构见[第3章](03-input-first-calc.md)，k 点与 SCF 收敛见[第5章](05-scf-kpoints-convergence.md)，结构弛豫见[第6章](06-optimization-eos-neb.md)，cube 文件与 `diff_gcube` 见[第10章](10-charge-mo-cube.md)，Berry 曲率与 Chern 数见[第12章](12-soc-topology.md)，Wannier 中心见[第14章](14-wannier.md)。

---

## 2. 原理速览

### 2.1 为什么不能直接积分偶极

对孤立分子，偶极矩 $\boldsymbol\mu=\int\mathbf r\,\rho(\mathbf r)\,\mathrm d^3r$ 定义明确。OpenMX 在每次计算的 `.out` 里都会打印 `Dipole moment (Debye)`，分子计算可以直接使用。

对周期晶体，把同样的积分放到一个晶胞里，结果依赖于晶胞边界怎么取：平移晶胞原点、或者换一种等价的晶胞选法，“偶极”都会改变。根源在于位置算符 $\hat{\mathbf r}$ 在 Bloch 态之间不是良定义的。可观测的物理量其实是极化的**变化**，即绝热变化过程中流过晶胞的电流积分：

$$\Delta\mathbf P=\int_0^1\frac{\mathrm d\mathbf P}{\mathrm d\lambda}\,\mathrm d\lambda$$

现代极化理论（King-Smith–Vanderbilt, 1993；Resta, 1994）证明，这个差值可以写成两个端点上各自的 Berry 相之差。

### 2.2 King-Smith–Vanderbilt 公式

极化分成离子（芯）和电子两部分（$e>0$ 为元电荷，$\Omega$ 为晶胞体积）：

$$\mathbf P=\mathbf P_\mathrm{ion}+\mathbf P_\mathrm{el},\qquad
\mathbf P_\mathrm{ion}=\frac{e}{\Omega}\sum_I Z_I^\mathrm{ion}\,\mathbf r_I$$

$$\mathbf P_\mathrm{el}=-\frac{f\,e}{(2\pi)^3}\sum_{n}^{\mathrm{occ}}\int_\mathrm{BZ}\mathbf A_n(\mathbf k)\,\mathrm d^3k,\qquad
\mathbf A_n(\mathbf k)=i\langle u_{n\mathbf k}|\nabla_{\mathbf k}u_{n\mathbf k}\rangle$$

其中 $Z_I^\mathrm{ion}$ 是赝势的价电荷（OpenMX 中就是 `.vps` 的价电子数，半芯态也计入），$f$ 是占据数（自旋非极化时 $f=2$，共线自旋极化或非共线时按每个自旋/旋量分别计 $f=1$），$\mathbf A_n$ 是 Berry 联络。

在离散 k 网格上，沿倒格矢 $\mathbf b_i$ 方向取一串（string）k 点 $\mathbf k_0,\dots,\mathbf k_{J-1}$（$\mathbf k_J=\mathbf k_0+\mathbf b_i$），每串的 Berry 相为

$$\phi_i(\mathbf k_\perp)=-\mathrm{Im}\,\ln\prod_{j=0}^{J-1}\det M^{(j,j+1)},\qquad
M^{(j,j+1)}_{mn}=\langle u_{m\mathbf k_j}|u_{n\mathbf k_{j+1}}\rangle$$

对垂直方向的所有串取平均得到 $\bar\phi_i$，电子极化为

$$\mathbf P_\mathrm{el}=-\frac{f\,e}{2\pi\,\Omega}\sum_{i=1}^{3}\bar\phi_i\,\mathbf a_i$$

$\mathbf a_i$ 是与 $\mathbf b_i$ 对偶的实空间晶格矢量。polB 输入的第一行（如 `9 9 9`）就是三个倒格矢方向上的分割数 $J$；沿某个方向计算时，另外两个方向的网格点构成串的数目（手册示例中 9 × 9 = 81 串）。OpenMX 使用非正交的局域原子轨道基，重叠矩阵 $M$ 要用 scfout 中的哈密顿量和重叠矩阵在每个 k 点求解本征矢后构造，所以 polB 的计算量主要在逐 k 点对角化。

### 2.3 极化量子与分支选择

$\bar\phi_i$ 只确定到 $2\pi$ 的整数倍，离子项也一样：把一个离子平移一个晶格矢量，$\mathbf P_\mathrm{ion}$ 会改变 $eZ_I\mathbf R/\Omega$。因此 P 是一个**格点集合**而不是一个数：

$$\mathbf P\;\equiv\;\mathbf P+\frac{e}{\Omega}\,\mathbf R,\qquad \mathbf R=n_1\mathbf a_1+n_2\mathbf a_2+n_3\mathbf a_3$$

量子 $\mathbf P_q=e\mathbf R/\Omega$。自旋简并、且所有离子价电荷都是偶数时，可能出现的差别只有 $2e\mathbf R/\Omega$ 的整数倍；只要有一个奇数价电荷的离子（例如 Na_CA19 价电子 9），或者体系自旋极化，就要按 $e\mathbf R/\Omega$ 来看。保险的做法是统一按 $e\mathbf R/\Omega$ 检查。

由此得到三条实用规则：

1. **单个结构的 P 绝对值没有物理意义**。polB 打印的 Core 列就是 $\sum_I Z_I\mathbf r_I$，随原点平移而改变；Electron 列的取值落在哪个分支也是任意的。
2. **只有差值 ΔP 可观测**，而且要沿一条“绝热路径”把 P(λ) 连成连续曲线，保证相邻两点之差远小于量子，才能确定正确的分支。
3. **中心对称结构的 P 只能是 0 或 $\mathbf P_q/2$**（模量子）。如果参考结构算出的 P 不在这些值上，说明结构的对称性有问题或计算没收敛。

以四方 BaTiO3 为例（a = 3.992 Å），沿 c 的量子为 $e c/\Omega=e/a^2=1602.18/3.992^2\approx100.5$ μC/cm²，而实验自发极化约 26 μC/cm²，两者同一量级，所以分支错一个就会得到完全错误的结论。

单位换算（本章脚本里使用）：

| 量 | 换算 |
|---|---|
| 1 e·Å | 4.80320 Debye |
| 1 e·bohr | 2.541746 Debye |
| 1 e/Å² | 1602.18 μC/cm² |
| 1 e/bohr² | 5721.48 μC/cm² |
| P（μC/cm²） | $=\;D(\text{e·bohr})\,/\,\Omega(\text{bohr}^3)\times5721.48$ |

### 2.4 Born 有效电荷

Born 有效电荷（动力学电荷）定义为宏观电场为零时，原子位移引起的极化变化：

$$Z^*_{I,\alpha\beta}=\frac{\Omega}{e}\,\frac{\partial P_\alpha}{\partial u_{I\beta}}\Big|_{\mathbf E=0}
=\frac{1}{e}\,\frac{\partial D_\alpha}{\partial u_{I\beta}}$$

$D_\alpha=\Omega P_\alpha$ 是一个晶胞的偶极。polB 同时输出 D（Debye）和 P（μC/cm²），用 D 计算时不需要体积。等价的定义是 $Z^*_{I,\alpha\beta}=\partial F_{I\beta}/\partial(eE_\alpha)$，即电场引起的原子受力，这在分子或板模型中可以用 `scf.Electric.Field` 算，体相则不行（见 2.5 节）。

几个有用的性质：

- **声学求和规则**：$\sum_I Z^*_{I,\alpha\beta}=0$（整体平移不改变极化）。这是检验计算质量最简单的办法。
- **芯电荷 + 电子屏蔽**：$Z^*=Z^\mathrm{ion}+Z^\mathrm{el}$。离子项恰好等于赝势价电荷，电子项是负的，抵消大部分。例如手册中 NaCl 的 Na 就是 9 − 7.95 ≈ 1.05（见 4.4 节）。
- **反常 Born 电荷**：钙钛矿铁电体中 B 位阳离子和沿 B–O 键方向的 O 的 Z\* 远大于名义电荷。立方 BaTiO3 的 Berry 相（LDA）计算值为 Z\*(Ba) = 2.75、Z\*(Ti) = 7.16、Z\*(O∥) = −5.69、Z\*(O⊥) = −2.11（Zhong、King-Smith、Vanderbilt, PRL 72, 3618 (1994)；Ghosez 等 1998 年的 DFPT 结果与之接近），Ti 的名义电荷只有 +4。可以顺手验证声学求和规则：2.75 + 7.16 − 5.69 − 2 × 2.11 = 0.00。巨大的 Z\* 正是 LO–TO 劈裂大、铁电不稳定性强的原因。

### 2.5 锯齿势表示的均匀外电场

均匀电场 $\mathbf E$ 对应线性势 $\varphi(\mathbf r)=-\mathbf E\cdot\mathbf r$，它与周期边界不相容。常用办法是锯齿势：在晶胞内线性变化，到晶胞边界处跳回原值。跳变的大小为 $e|\mathbf E|L$（L 为场方向上的晶胞长度）。

| 体系 | 不连续面位置 | 是否可用 |
|---|---|---|
| 分子（四周有真空） | 真空中 | 可用，手册硝基苯示例 |
| 有真空层的板模型 | 应在真空中，需用 `v0.cube` 确认（6.5 节） | 可用，但两侧真空的势被人为改变 |
| 无真空的体相 | 必然切过材料 | 手册第 30 章：可能引起数值不稳定（numerical instability may be caused），不应使用 |

另外两个要注意的问题：

- **场致电子溢出**：跳变 $e|\mathbf E|L$ 很容易达到几十 eV（10 GV/m × 24 Å = 24 eV）。如果真空低势一侧的势降到分子 HOMO 以下，电子会“漏”进真空，得到异常大的偶极和异常低的能量。场强和盒子都不宜过大。
- **线性响应范围**：分子极化率 $\alpha_{\alpha\beta}=\partial\mu_\alpha/\partial E_\beta$ 用 ±F 中心差分计算，F 取 1–3 GV/m（约 0.002–0.006 a.u.）；10 GV/m = 1 V/Å 已明显进入非线性区。

OpenMX 4.0 源码（`Input_std.c`）中，输入单位 GV/m 乘以 $0.1944688\times10^{-2}$ 换算为 Hartree/(e·bohr)（即 1 a.u. ≈ 514.2 GV/m）。三个分量分别是**沿 a、b、c 晶格矢量方向**的分量，程序会把它们换算成笛卡尔方向，并在标准输出打印方向单位矢量和大小。对正交晶胞，三个分量就是 x、y、z。

### 2.6 OpenMX 中与电场有关的工具一览

| 工具 | 物理量 | 适用体系 | 本书位置 |
|---|---|---|---|
| `polB`（Berry 相） | 零场下的 P、ΔP、Z\* | 体相绝缘体（共线、非共线均支持） | 本章 3–5 节 |
| `scf.Electric.Field`（锯齿势） | 有限场下的密度、偶极、力 | 分子、有真空的板 | 本章 6–7 节 |
| ESM `on4` | 两块理想金属间的均匀场 | 板模型 | [第16章](16-surface-esm-workfunction.md) 7.3 节 |
| 光电导/介电函数（Kubo–Greenwood） | $\varepsilon(\omega)$，ω→0 给出 $\varepsilon_\infty$ | 体相（共线） | [第9章](09-dos-pdos-cohp-optics.md) |
| `Dipole moment (Debye)`（`.out`） | 孤立体系的偶极 | 分子 | 本章 6 节 |

> [!IMPORTANT]
> 关于“OpenMX 4.0 新增基于 Berry 相的有限电场方法”的说法：笔者核对了 OpenMX 4.0 手册目录（第 30 章 Electric Field 与第 44 章 Berry phase polarization）、4.0 关键词总表和 What's new 中 4.0 的 24 项更新，**均未找到**体相有限电场（Souza–Íñiguez–Vanderbilt 或 Umari–Pasquarello 类方法）的关键词或说明。4.0 中能施加电场的只有锯齿势 `scf.Electric.Field` 和 ESM `on4`，polB 是零场后处理。若你需要体相的有限电场计算（如 $\varepsilon_\infty$、压电或电场下的晶格响应），可以 (1) 在 OpenMX 内用 polB 求 Z\*、再结合声子计算得到离子贡献；(2) 改用支持体相有限电场的程序，如 VASP 的 `LCALCEPS`/`EFIELD_PEAD` 或 Quantum ESPRESSO 的 `lelfield`。以后版本若加入此功能，以官方手册为准。

---

## 3. polB 的编译与工作流程

### 3.1 编译

polB 和 openmx 共用 `source/makefile` 中的 `CC`、`LIB` 设置（含 MPI），在源码目录编译：

```bash
cd ~/openmx4.0/source
make polB              # or: make all  (builds openmx and the post-processing codes)
find ~/openmx4.0 -name polB -type f    # the manual says the binary is placed in work/
```

| 版本 | 说明 |
|---|---|
| 3.9（patch 3.9.9） | polB 已存在，手册章节名为 “Macroscopic polarization by Berry's phase”，3.9 版对效率和稳定性做过改进 |
| 4.0（patch 4.0.1） | 手册第 44 章；polB 读入时打印 `The file format of the SCFOUT file:  3`，并列出读该格式的后处理程序（jx、polB、kSpin、Z2FH、calB 等） |

> [!WARNING]
> 后处理程序要与生成 scfout 的 openmx 来自同一版本、同一次编译。用 3.9 编译的 polB 去读 4.0 的 scfout（或反过来）可能读错数据，也可能直接崩溃。

### 3.2 三步流程

```bash
# (1) SCF with HS.fileout on  ->  nacl.scfout
mpirun -np 8 openmx NaCl.dat -nt 1 > nacl.std
# (2) Berry phase, interactive
mpirun -np 8 polB nacl.scfout
# (2') or batch mode with a two-line input file
printf "9 9 9\n1 1 1\n" > in
mpirun -np 8 polB nacl.scfout < in > out
```

SCF 输入里与 polB 有关的关键词：

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `HS.fileout` | `on` | 写出 `System.Name.scfout`（实空间 H、S、密度矩阵等），polB 的唯一输入。默认 off |
| `scf.EigenvalueSolver` | `band` | 体相用 band；polB 会打印 `Previous eigenvalue solver = Band` |
| `scf.Kgrid` | 与普通 SCF 相同，做收敛测试 | 只决定 SCF 密度的质量；polB 的 k 网格由它自己的输入决定，与此无关 |
| `scf.criterion` | `1.0e-9`～`1.0e-10` | 极化差和 Z\* 是小量之差，SCF 要收紧 |
| `scf.ElectronicTemperature` | 300 K 或更低 | Berry 相只对占据态有定义，要求绝缘体；温度过高会出现分数占据 |
| `scf.SpinPolarization` | `off` / `on` / `nc` | polB 对共线和非共线都支持 |
| `level.of.fileout` | `0` | polB 不需要 cube 文件，关掉可以省磁盘 |

polB 交互模式下依次询问两项（手册第 44 章）：

| 提示 | 输入示例 | 含义 |
|---|---|---|
| `Specify the number of grids to discretize reciprocal a-, b-, and c-vectors` | `9 9 9` | 三个倒格矢方向的 k 分割数 |
| `Specify the direction of polarization as reciprocal a-, b-, and c-vectors` | `1 1 1` | 是否计算沿 a、b、c 的 Berry 相（1 算，0 不算） |

> [!TIP]
> polB 源码（`polB.c`）只用两次 `scanf` 读这两行，与自旋设置无关，所以共线自旋极化和非共线体系的批处理输入同样是两行。要注意第一行：三个数都必须 ≥ 1。如果第一行含 0 或负数，程序会打印 `The number of grids should be one or more` 并要求重新输入；批处理模式下标准输入已经读完，它会不停重复这句提示，作业看起来像“卡住”，输出文件也越写越大。`in` 文件为空或少一行时，读入的值没有定义，结果同样不可信。

### 3.3 读懂输出

polB 先回显 scfout 的信息（原子数、化学势 `ChemP`、自旋处理方式、`tv1`–`tv3` 实空间晶格矢量（Bohr）、`rtv1`–`rtv3` 倒格矢、`Cell_Volume`（Bohr³）），然后逐串计算，最后给出两个表。下面是手册第 44 章 NaCl 示例（Na 沿 x 位移 +0.05 Å）的结果，引自官方手册：

```text
*******************************************************
              Electric dipole  (Debye) : Berry phase
*******************************************************
 Absolute dipole moment      164.07947170
               Background        Core             Electron          Total
 Dx           -0.00000000       96.80863302       -1.90924436       94.89938866
 Dy           -0.00000000       94.64718996       -0.00000658       94.64718338
 Dz           -0.00000000       94.64718996       -0.00000436       94.64718560
***************************************************************
              Electric polarization (muC/cm^2) : Berry phase
***************************************************************
               Background        Core             Electron          Total
 Px           -0.00000000      723.82242623      -14.27511002      709.54731620
 ...
```

| 列 | 含义 |
|---|---|
| Background | 带电体系（`scf.system.charge` ≠ 0）的均匀背景电荷贡献；中性体系为 0 |
| Core | $\sum_I Z_I\mathbf r_I$，离子贡献，**随原点选择而变** |
| Electron | Berry 相给出的电子贡献，**落在哪个分支是任意的** |
| Total | 前三列之和；单独看没有意义，只用于求差 |

可以验证两个表的一致性：94.89938866 D ÷ 2.541746 = 37.336 e·bohr，除以体积 301.066 bohr³ 再乘 5721.48，得到 709.55 μC/cm²，与 Px 一致。

手册的示例输出到 P 表为止。按 polB 源码（这里核对的是基于 4.0 的 [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu) 中的 `source/polB.c`），还有两项输出值得留意（量子表在 P 表之后，金属警告在 D 表与 P 表之间）：

- 极化量子表：标题为 `Note: Bulk electric polarization has arbitrariness of eRi/Omega.`，下面三行 `eRix/Omega`、`eRiy/Omega`、`eRiz/Omega` 依次给出三个晶格矢量 R1、R2、R3 对应的 $e\mathbf R_i/\Omega$ 的 x、y、z 分量（μC/cm²）。这就是 2.3 节的量子 $\mathbf P_q$（f = 1），可以直接拿来做分支展开，不必自己算。
- 金属警告：polB 在每个 k 点数一数化学势 `ChemP` 以下的本征态个数，若与价电子数给出的占据带数（自旋非极化为 N/2，非共线为 N）不一致，计算过程中会打印 `Metallic !!`，D 表与 P 表之间会出现 `##!! Warning!! The system may be metallic! Result is not correct!`。出现这句时结果没有意义，先检查带隙和电子温度。注意：共线自旋极化（`scf.SpinPolarization on`）时，源码直接用每个 k 点实际的占据数，**不做这项检查**，需要自己从 `.out` 的本征值或 DOS 确认两个自旋通道都有带隙。

```bash
grep -A5 "arbitrariness of eRi/Omega" out      # polarization quanta (muC/cm^2)
grep "may be metallic" out                     # must print nothing
```

---

## 4. 实战一：NaCl 中 Na 的 Born 有效电荷

示例目录：[15-nacl-polb-born](../examples/15-nacl-polb-born/README.md)。

### 4.1 输入文件 `NaCl_polB.dat`

```text
System.CurrrentDirectory    ./           # note: official spelling has three r
System.Name                 nacl_x0
level.of.stdout             1
level.of.fileout            0            # no cube files needed here
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/

Species.Number      2
<Definition.of.Atomic.Species
  Na   Na9.0-s3p2d1    Na_CA19
  Cl   Cl7.0-s2p2d1    Cl_CA19
Definition.of.Atomic.Species>

Atoms.Number        2
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1  Na   0.00000   0.00000   0.00000   4.5  4.5
  2  Cl   2.82000   2.82000   2.82000   3.5  3.5
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  0.00000   2.82000   2.82000
  2.82000   0.00000   2.82000
  2.82000   2.82000   0.00000
Atoms.UnitVectors>

scf.XcType                  LDA          # must match *_CA19 pseudopotentials
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0        # K
scf.energycutoff            220.0        # Ry, set explicitly
scf.maxIter                 200
scf.EigenvalueSolver        band
scf.Kgrid                   9 9 9
scf.Mixing.Type             rmm-diisk
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          7
scf.Mixing.StartPulay       8
scf.criterion               1.0e-10      # Hartree

MD.Type                     nomd
HS.fileout                  on
```

| 关键词 | 设置 | 说明 |
|---|---|---|
| `Na9.0-s3p2d1` / `Na_CA19` | DFT_DATA19 Standard 基组 | Na_CA19 价电子 9（2s² 2p⁶ 3s¹），坐标行写 4.5 4.5 |
| `Cl7.0-s2p2d1` / `Cl_CA19` | Quick 基组 | Cl 价电子 7；Standard 为 `Cl7.0-s2p2d1f1`，收敛测试时可以换上 |
| `scf.XcType LDA` | 与 `*_CA19` 配套 | 与手册示例一致，便于对照；也可整体换成 GGA-PBE + `*_PBE19` |
| 晶格 | a0 = 5.64 Å | 接近室温实验值；手册 `NaCl.dat` 用 5.63 Å |
| `scf.Kgrid 9 9 9` | — | 4.0 默认 `scf.Generation.Kpoint GCenter`；用 3.9 时为 Regular 网格，结果会有微小差异 |
| `scf.criterion 1.0e-10` | 紧 | Z\* 来自约 0.5 D 量级的差值 |
| `HS.fileout on` | 必需 | 生成 `nacl_x0.scfout` |

> [!NOTE]
> 手册第 44 章那组输出对应的是 Na **已经沿 x 位移 +0.05 Å** 的结构，不是平衡结构：Core 列 Dx 比 Dy 大 2.16144 D = 0.45 e·Å = 9 × 0.05 Å，而手册后文列出的 x = +0.05 Å 的 Px 也正是 94.89938866 D。不要把那组输出误当作平衡结构的结果。

### 4.2 生成位移结构并运行

`make_disp.sh` 只改 Na 的 x 坐标和 `System.Name`，生成 u = −0.05、0、+0.05 Å 三个输入：

```bash
cd OpenMX/examples/15-nacl-polb-born
bash make_disp.sh                       # NaCl_dx_m0.05.dat NaCl_dx_p0.00.dat NaCl_dx_p0.05.dat
for u in m0.05 p0.00 p0.05; do
  mpirun -np 8 openmx NaCl_dx_$u.dat -nt 1 > NaCl_dx_$u.std 2>&1
  grep -q "The calculation was normally finished." NaCl_dx_$u.std || echo "$u FAILED"
  mpirun -np 8 polB nacl_dx_$u.scfout < polB.in > polB_NaCl_dx_$u.out
done
```

集群上直接 `sbatch job_slurm.sh`，脚本会依次完成上述步骤并调用 `born_charge.py`。polB 的 MPI 并行效率不受原子数限制（手册：进程数超过原子数时仍能高效并行）。手册中 NaCl 的 polB 一次运行约 76 s，作为量级参考。

### 4.3 检查输出

```bash
grep -A3 "SCF history" nacl_dx_p0.00.out | head          # converged?
grep -i "Chemical Potential" nacl_dx_p0.00.out
grep -A6 "Electric dipole" polB_NaCl_dx_p0.00.out
grep "may be metallic" polB_NaCl_dx_*.out                # must print nothing
```

- `p0.00` 是立方结构，Total 列 Dx = Dy = Dz（到小数点后多位）。
- `p0.05` 中只有 Dx 变化，Dy、Dz 与 `p0.00` 相同；`m0.05` 的 Dx 变化与 `p0.05` 大小相等、方向相反。
- 若三个点的 Dx 不在一条直线上（`born_charge.py` 会给出警告），先怀疑分支跳变，再怀疑 SCF 或 k 网格。

### 4.4 计算 Z\*，并做芯/电子分解

```bash
python3 born_charge.py --atom Na --dir x \
  -0.05:polB_NaCl_dx_m0.05.out 0.00:polB_NaCl_dx_p0.00.out 0.05:polB_NaCl_dx_p0.05.out
```

脚本对 $D_\alpha(u)$ 做线性拟合，输出 $Z^*_{xx}$、$Z^*_{yx}$、$Z^*_{zx}$ 的 total/core/electron 三部分，并用 P 表和 `Cell_Volume` 交叉核对。用 `python born_charge.py -h` 查看全部选项。

用手册第 44 章公布的数字手算一遍，可以看清各部分的意义（u 从 0 到 +0.05 Å）：

| 部分 | ΔDx（Debye） | ΔDx（e·Å） | ΔDx / Δu |
|---|---|---|---|
| Core | 96.80863 − 94.64719 = 2.16144 | 0.45000 | **9.000** = Na_CA19 价电荷 |
| Electron | −1.90924 − (≈0) | −0.39749 | **−7.950** |
| Total | 0.25220 | 0.05251 | **1.050** |

手册给出的 Z\*(Na) = 1.05，表中另一计算值为 1.09，实验值为 1.12。手册用的是 ±0.05 Å 的中心差分：Px = 94.39498187、94.64718713、94.89938866 D（x = −0.05、0、+0.05 Å），

$$Z^*_{xx}=\frac{(94.89938866-94.39498187)/2.54174776}{0.1/0.529177}=1.050$$

NaCl 是立方结构，非对角元为 0，且 $Z^*_{xx}=Z^*_{yy}=Z^*_{zz}$。

> [!TIP]
> 完整的 Z\* 张量需要对每个不等价原子分别沿 x、y、z 位移（3 × 2 次 SCF + polB），每次位移读出 Px、Py、Pz，得到张量的一列。高对称结构可以用对称性减少计算量。

### 4.5 收敛与检验

| 检验项 | 做法 | 判据 |
|---|---|---|
| polB k 网格 | `polB.in` 第 1 行取 6/9/12 | Z\* 变化 < 0.01 |
| SCF k 网格 | `scf.Kgrid` 取 7/9/11 | 同上 |
| 位移大小 | ±0.02 与 ±0.05 Å | 结果一致（线性区） |
| 基组 | Cl 换 `s2p2d1f1`，Na 换 `s3p2d2` | 变化在百分之几以内 |
| 声学求和规则 | 再算 Cl 的 Z\* | Z\*(Na) + Z\*(Cl) ≈ 0 |

---

## 5. 实战二：四方 BaTiO3 的自发极化

示例目录：[15-batio3-polarization](../examples/15-batio3-polarization/README.md)。

### 5.1 思路：同一晶胞中的 λ 路径

立方 BaTiO3 是中心对称的顺电相；室温四方相（P4mm）中 Ti 和 O 沿 c 反向偏离中心，产生沿 c 的自发极化。自发极化定义为

$$\mathbf P_s=\mathbf P(\lambda=1)-\mathbf P(\lambda=0)$$

λ = 0 是中心对称参考结构，λ = 1 是极化结构。本例做了两个简化，让结果只取决于离子位移：

- 参考结构使用**和极化结构相同的四方晶胞**，只把原子放回中心对称位置（P4/mmm）。若用真正的立方晶胞，路径上还会同时改变晶胞，引入额外的应变贡献，也不便插值；
- 在 λ = 0, 0.25, 0.5, 0.75, 1 上逐点算 P，再把 P(λ) 连成连续曲线，自然地选出正确分支。

### 5.2 输入文件

中心对称参考结构 `BaTiO3_cubic_ref.dat`（完整）：

```text
System.CurrrentDirectory    ./           # note: official spelling has three r
System.Name                 bto_centro
level.of.stdout             1
level.of.fileout            0
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/

Species.Number      3
<Definition.of.Atomic.Species
  Ba   Ba10.0-s3p2d2   Ba_CA19
  Ti   Ti7.0-s3p2d1    Ti_CA19
  O    O6.0-s2p2d1     O_CA19
Definition.of.Atomic.Species>

Atoms.Number        5
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Ba   0.000000   0.000000   0.000000    5.0  5.0
  2  Ti   0.500000   0.500000   0.500000    6.0  6.0
  3  O    0.500000   0.500000   0.000000    3.0  3.0
  4  O    0.500000   0.000000   0.500000    3.0  3.0
  5  O    0.000000   0.500000   0.500000    3.0  3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  3.992000   0.000000   0.000000
  0.000000   3.992000   0.000000
  0.000000   0.000000   4.036000
Atoms.UnitVectors>

scf.XcType                  LDA          # must match *_CA19 pseudopotentials
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0        # K
scf.energycutoff            250.0        # Ry, set explicitly
scf.maxIter                 200
scf.EigenvalueSolver        band
scf.Kgrid                   8 8 8
scf.Mixing.Type             rmm-diisk
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          7
scf.Mixing.StartPulay       8
scf.criterion               1.0e-10      # Hartree

MD.Type                     nomd
HS.fileout                  on
```

极化结构 `BaTiO3_tetra_polB.dat` 与之相比只有坐标块和 `System.Name` 不同：

```text
System.Name                 bto_polar
<Atoms.SpeciesAndCoordinates
  1  Ba   0.000000   0.000000   0.000000    5.0  5.0
  2  Ti   0.500000   0.500000   0.515000    6.0  6.0
  3  O    0.500000   0.500000   0.975000    3.0  3.0
  4  O    0.500000   0.000000   0.485000    3.0  3.0
  5  O    0.000000   0.500000   0.485000    3.0  3.0
Atoms.SpeciesAndCoordinates>
```

这组 z 坐标只是量级合理的**初猜**（Ti 上移，O 下移），不是实验或计算结果。`BaTiO3_tetra_relax.dat` 从这组初猜出发，用 `MD.Type RF` 在固定晶胞下弛豫内坐标，并用 `<MD.Fixed.XYZ` 固定 Ba 以消除整体平移：

```text
MD.Type                     RF
MD.maxIter                  100
MD.Opt.criterion            1.0e-4       # Hartree/Bohr

<MD.Fixed.XYZ
  1   1 1 1
  2   0 0 0
  3   0 0 0
  4   0 0 0
  5   0 0 0
MD.Fixed.XYZ>
```

| 关键词 / 设置 | 值 | 说明 |
|---|---|---|
| 基组 | Ba10.0-s3p2d2、Ti7.0-s3p2d1、O6.0-s2p2d1 | Ti、O 为 DFT_DATA19 推荐表（手册 12.3 节 Table 1）的 Standard 档；Ba 用 `Ba10.0.pao` 中的 s3p2d2 子集，正式计算前请对照 Table 2 的推荐档并做基组收敛测试 |
| 价电子 | Ba 10、Ti 12、O 6（合计 40） | 坐标行上下自旋各一半；全为偶数，所以自旋简并时差别只会是 2eR/Ω 的倍数 |
| 泛函 | LDA + `*_CA19` | GGA-PBE 会显著高估 BaTiO3 的四方性，固定实验晶胞时 LDA 更常用 |
| 晶格 | a = 3.992、c = 4.036 Å | 室温实验值的近似；所有 λ 共用 |
| `scf.Kgrid` | 8 8 8 | 需收敛测试 |
| `polB.in` | `6 6 8` / `0 0 1` | 只算沿 c 的 Berry 相，c 方向 k 点更密 |

### 5.3 运行

```bash
cd OpenMX/examples/15-batio3-polarization
sbatch --export=ALL,STAGE=relax job_slurm.sh          # optional but recommended
grep -A8 "Fractional coordinates of the final structure" bto_relax.out

sbatch --export=ALL,STAGE=path job_slurm.sh
#   = python3 make_path.py BaTiO3_cubic_ref.dat BaTiO3_tetra_polB.dat \
#         --lambdas 0 0.25 0.5 0.75 1.0 --from-out bto_relax.out
#     + SCF and polB on every BaTiO3_lam*.dat
#     + python3 polar_path.py 0.00:polB_bto_lam0.00.out ... 1.00:polB_bto_lam1.00.out
```

`make_path.py` 按最小像插值（所以 0.975 与 0.000 之间按 −0.025 处理），`polar_path.py` 读出每个 λ 的 Pz，按 $Q=f\,e\,c/\Omega$ 把相邻点之差折到最小，并打印 ΔP。两个脚本都可以用 `python <脚本名> -h` 查看用法。

### 5.4 结果判读

下表演示 `polar_path.py` 的输出格式（示意输出，数值仅供参考，省略号处为你的计算结果）：

```text
Polarization quantum along z (f=1, R=tv3): 100.540 muC/cm^2
 lambda         core     electron   total(raw) total(unwrapped)  file
  0.000          ...          ...          ...              ...  polB_bto_lam0.00.out
  0.250          ...          ...          ...              ...  polB_bto_lam0.25.out
  ...
  1.000          ...          ...          ...              ...  polB_bto_lam1.00.out

Delta P = P(lambda=1) - P(lambda=0) = ... muC/cm^2
```

第一行的量子是确定的：$e\,c/\Omega=e/a^2=1602.18/3.992^2\approx100.5$ μC/cm²，应与 polB 输出末尾 `eRiz/Omega` 行 R3 列的值一致（`grep -A5 "arbitrariness" polB_bto_lam0.00.out`）。判读要点：

1. 中心对称的 λ = 0 结构，P 只能是 0 或半个“真正的”量子。本例自旋简并、价电荷全为偶数（Ba 10、Ti 12、O 6），真正的量子是 $Q_2=2ec/\Omega\approx201$ μC/cm²，λ = 0 的 Pz 只能 ≡ 0 或 100.5 μC/cm²（模 201）；换成脚本默认的 f = 1 量子 Q ≈ 100.5 μC/cm² 来看，就是 **Pz ≡ 0（模 100.5）**。这一点可以分开看：Core 列在 λ = 0 时为 $\sum_I Z_I z_I=(12+6+6)\times c/2=12c$，对应 12 个 Q，本身 ≡ 0；电子部分每个 Berry 相只能是 0 或 π，乘以 f = 2 后给出 0 或一个完整的 Q，也 ≡ 0。若 λ = 0 的 Pz 明显偏离 Q 的整数倍，说明参考结构不是严格中心对称，或 SCF/k 网格没收敛。
2. unwrapped 列应随 λ 单调、近似线性地变化；`WARNING: a step along the path exceeds 0.3 quantum` 表示路径点太稀。
3. ΔP 与实验值（约 26 μC/cm²，Wieder 1955）比较时，注意计算结果强烈依赖于弛豫后的离子位移、晶胞和泛函。**未弛豫的初猜坐标只能用来演示流程**。
4. 用 Z\* 做线性估算可以交叉检验：$\Delta P_z\approx\frac{e}{\Omega}\sum_I Z^*_{I,zz}\,\Delta z_I$。把 2.4 节的立方 Z\* 和弛豫后的位移代入，得到的数应与 Berry 相结果同量级（Z\* 在四方相中会变小，所以线性估算一般偏大）。

### 5.5 推广到 PbTiO3 和其他铁电体

- PbTiO3：把 Ba 换成 Pb（赝势与基组以 DFT_DATA19 实际文件和推荐表为准），换成 PbTiO3 的四方晶格常数，流程不变。PbTiO3 的位移和 c/a 都更大，λ 路径要取得更密。
- 极化方向不沿某一晶格矢量时（如菱方相），把 `polB.in` 第 2 行改为 `1 1 1`，三个方向都要做分支展开（`polar_path.py --comp x/y/z --R 1/2/3`）。
- 层状、二维铁电或需要分层极化时，用混合 Wannier 中心（见第 8 节 HWC）更直观。

---

## 6. 实战三：硝基苯在外加电场中

示例目录：[15-nitrobenzene-efield](../examples/15-nitrobenzene-efield/README.md)。思路与手册第 30 章的 `Nitro_Benzene.dat` 一致，输入为本书重新编写。

### 6.1 输入文件 `Nitro_Benzene_EF.dat`

```text
System.CurrrentDirectory    ./           # note: official spelling has three r
System.Name                 nben_f0
level.of.stdout             1
level.of.fileout            1            # writes tden/v0/... cube files
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/

Species.Number      4
<Definition.of.Atomic.Species
  H    H6.0-s2p1       H_PBE19
  C    C6.0-s2p2d1     C_PBE19
  N    N6.0-s2p2d1     N_PBE19
  O    O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>

Atoms.Number        14
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
   1  C    10.61000   10.00000    7.00000   2.0 2.0
   2  C    11.30500   11.20378    7.00000   2.0 2.0
   3  C    12.69500   11.20378    7.00000   2.0 2.0
   4  C    13.39000   10.00000    7.00000   2.0 2.0
   5  C    12.69500    8.79622    7.00000   2.0 2.0
   6  C    11.30500    8.79622    7.00000   2.0 2.0
   7  H    10.76500   12.13908    7.00000   0.5 0.5
   8  H    13.23500   12.13908    7.00000   0.5 0.5
   9  H    14.47000   10.00000    7.00000   0.5 0.5
  10  H    13.23500    7.86092    7.00000   0.5 0.5
  11  H    10.76500    7.86092    7.00000   0.5 0.5
  12  N     9.14000   10.00000    7.00000   2.5 2.5
  13  O     8.56724   11.07720    7.00000   3.0 3.0
  14  O     8.56724    8.92280    7.00000   3.0 3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  24.0   0.0   0.0
   0.0  20.0   0.0
   0.0   0.0  14.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0        # K
scf.energycutoff            220.0        # Ry, set explicitly
scf.maxIter                 200
scf.EigenvalueSolver        cluster
scf.Mixing.Type             rmm-diis
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.40
scf.Mixing.History          7
scf.Mixing.StartPulay       6
scf.criterion               1.0e-9       # Hartree

scf.Electric.Field          0.0 0.0 0.0

MD.Type                     nomd
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `scf.Electric.Field` | `F 0.0 0.0`（GV/m） | 默认 `0.0 0.0 0.0`。三个分量沿 a、b、c 晶格矢量；1 V/Å = 10 GV/m。手册：符号按作用于电子的场定义 |
| `scf.EigenvalueSolver` | `cluster` | 孤立分子；不需要 `scf.Kgrid` |
| 盒子 | 24 × 20 × 14 Å，分子居中 | 场沿 a，每侧留出约 8.5 Å 以上真空，让不连续面落在真空中 |
| `level.of.fileout` | `1` | 写出 `tden.cube`、`v0.cube` 等，供 `diff_gcube` 使用 |
| `scf.criterion` | `1.0e-9` | 偶极变化量小，SCF 要收紧 |
| 几何 | 理想化平面结构 | C–C 1.39、C–H 1.08、C–N 1.47、N–O 1.22 Å，O–N–O 124°；定量计算前先在零场下弛豫 |
| 电子数 | 46 | 坐标行上下自旋之和：C 4、H 1、N 5、O 6 |

### 6.2 运行一组场强

```bash
cd OpenMX/examples/15-nitrobenzene-efield
bash run_fields.sh                     # writes Nitro_Benzene_EF_{m10,m2,p2,p10}.dat
RUN=1 NP=8 bash run_fields.sh          # write and run, or: sbatch job_slurm.sh
grep -A2 "Applied External Electric Field" Nitro_Benzene_EF_p10.std
```

`run_fields.sh` 只改两行：`System.Name`（如 `nben_p10`）和 `scf.Electric.Field`（如 `10 0.0 0.0`）。加场且 `level.of.stdout` ≥ 1 时，标准输出中会出现（格式取自 `Input_std.c`，与 3.9.9 源码一致）：

```text
<Applied External Electric Field>
  direction (x, y, z)  =    1.00000    0.00000    0.00000
  magnitude (10^9 V/m) =   10.00000
```

没有这一段，说明关键词没有生效（拼写错误，或三个分量都是 0）。

### 6.3 偶极与极化率

每个 `.out` 中都有偶极块，格式如下（数值为示意）：

```text
                    Dipole moment (Debye)
 Absolute D        ...
                      Dx                Dy                Dz
 Total              ...               ...               ...
 Core               ...
 Electron           ...
 Back ground        ...
```

```bash
grep -A12 "Dipole moment (Debye)" nben_f0.out | grep -E "Absolute|Total"
python3 dipole_alpha.py 0:nben_f0.out -2:nben_m2.out 2:nben_p2.out -10:nben_m10.out 10:nben_p10.out
```

`dipole_alpha.py` 读 Total 行，只用 |F| ≤ 3 GV/m 的点做线性拟合：

$$\alpha_{xx}=\frac{\mu_x(+F)-\mu_x(-F)}{2F},\qquad
\alpha\,[\text{Å}^3]=0.148185\times\alpha\,[\text{bohr}^3]$$

±10 GV/m 的点单独列出，用来观察非线性（超极化率）。极化率必须为正：如果脚本得到负的 $\alpha_{xx}$，说明在这套约定下，正的关键词值对应 −x 方向的物理电场。手册只说“符号按作用于电子的场定义”，用偶极变化的方向来实际确认，是最不容易出错的做法。

### 6.4 场致电荷密度差

```bash
gcc $(find ~/openmx4.0 -name diff_gcube.c | head -1) -lm -o diff_gcube
./diff_gcube nben_p10.tden.cube nben_f0.tden.cube drho_p10.cube
./diff_gcube nben_m10.tden.cube nben_f0.tden.cube drho_m10.cube
```

`diff_gcube` 的用法是 `diff_gcube input1.cube input2.cube output.cube`，输出 input1 − input2（手册第 66 章；手册说 `diff_gcube.c` 在 `work/` 目录，上面用 `find` 是为了不依赖具体位置）。只改场强、不改结构和晶胞时，两次计算的网格完全一致，可以逐点相减。如果屏幕上出现 `Found a difference in ...` 之类的提示，说明两个 cube 的网格不一致；源码中询问是否继续的那段已被注释掉，程序不会停下来，往往**照样写出** output.cube（网格数不同时也可能直接崩溃），但这个差分没有意义，必须查明原因后重算。在 VESTA 中用正负两种颜色显示等值面（如 ±0.002 e/bohr³）：

- 手册第 30 章图 33 的结论是：电子转移主要发生在 −NO2 的 O 原子与对位 C、对位 H 之间；
- `drho_p10` 与 `drho_m10` 的增减区域大致互换；如果两者明显不对称，说明 10 GV/m 已进入强非线性区，或电子溢出到了真空中（见 2.5 节）。

做平面平均可以把 Δρ 化成一维曲线。可以用第16章的 `gcube2oned`，也可以用仓库脚本 [omx_cube_tools.py](../scripts/omx_cube_tools.py)（`python3 ../../scripts/omx_cube_tools.py -h` 查看全部子命令）：

```bash
python3 ../../scripts/omx_cube_tools.py avg drho_p10.cube --axis a -o drho_p10_avg_a.dat   # planar average along a
python3 ../../scripts/omx_cube_tools.py integrate drho_p10.cube                           # should be ~0 (charge conserved)
```

### 6.5 检查不连续面的位置

```bash
gcc $(find ~/openmx4.0 -name gcube2oned.c | head -1) -lm -o gcube2oned   # see chapter 16, sec. 5.3
./gcube2oned nben_p10.v0.cube 1 > v0_p10_1d.txt   # average over b and c, profile along a
./gcube2oned nben_f0.v0.cube  1 > v0_f0_1d.txt
paste v0_p10_1d.txt v0_f0_1d.txt | awk '{print $2, $3-$6}' > dv_1d.txt   # position (Ang), dV (Hartree)
```

`gcube2oned` 输出三列：网格序号、从 cube 原点量起的位置（Å）、平面平均值（Hartree），所以 `paste` 之后第 3 列与第 6 列相减就是场致的势变化。

如果外场项包含在 `v0` 中，`dv_1d.txt` 应是一段线性斜坡加一个跳变，跳变位置必须落在真空中、远离分子。若这两个文件的差几乎为零，说明当前版本没有把外场写进 `v0`，此时只能根据盒子几何来保证分子远离晶胞边界（本例分子居中，两侧各约 8.5 Å 真空）。

> [!TIP]
> 用同样的有限场计算，还可以求分子的原子极性张量 $\partial F_{I\beta}/\partial E_\alpha$（读 `.out` 中 `<coordinates.forces` 块的原子受力，对 ±F 做差分），它与红外强度直接相关，也是分子版本的 Born 电荷。
>
> 手册第 30 章还说明，锯齿势电场既可用于 SCF，也可用于结构优化：把 `MD.Type` 改成 `RF` 等优化方法，就能得到分子在外场下弛豫后的几何（场强同样要保证不发生电子溢出）。

---

## 7. 板模型中的电场：锯齿势还是 ESM on4

对有真空层的 slab，可以用两种办法加均匀场：

| | `scf.Electric.Field` | ESM `on4` |
|---|---|---|
| 边界条件 | 三维周期 + 锯齿势 | 法向两侧为理想金属，电势差 `ESM.potential.diff`（eV） |
| 法向 | 任意晶格矢量方向 | a 轴（`<Atoms.UnitVectors` 第一行）作法向，且须平行于某个笛卡尔轴；`ESM.direction x\|y\|z`（默认 x）指明是哪个轴（见第16章 3.1 节） |
| 不连续面 | 有，必须落在真空中 | 无 |
| 板自身的偶极 | 不修正（4.0 关键词表中没有偶极修正关键词） | 由边界条件处理 |
| 带电体系 | 用均匀背景 | 可以 |
| 推荐用途 | 分子、对称板的小场响应 | 电容器、栅压、非对称板 |

> [!WARNING]
> 对无真空的体相晶体加 `scf.Electric.Field`，不连续面必然切过材料，手册第 30 章指出这可能引起数值不稳定；即使 SCF 收敛，体相中的锯齿势也不代表真实的均匀宏观电场。体相的介电/压电响应请用 polB + Born 电荷（零场），或改用支持体相有限电场的程序（见 2.6 节）。

slab + 锯齿势的几点建议：

1. 真空至少 15 Å，场强 ≤ 5 GV/m，避免真空低势一侧的势降到表面态以下；
2. 用 `v0.cube` 的平面平均（6.5 节方法）确认两侧真空各有一段平坦或线性区；
3. 非对称板、吸附体系、极性表面一律优先用 ESM（第16章）。

---

## 8. 与其他工具联动

| 工具 | 输入 | 用途 |
|---|---|---|
| [Ncmexp2717/HWC](https://github.com/Ncmexp2717/HWC) | scfout（`HS.fileout on`） | 复制进 OpenMX `source/` 后 `make HWC`；计算混合 Wannier 中心，得到逐层（layer-resolved）极化，适合超晶格、二维铁电和拓扑分析 |
| [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl) | OpenMX 哈密顿量 | 构建紧束缚模型，计算 Berry 曲率、光学和非线性响应（如 shift current） |
| OpenMX MLWF/CWF | SCF 输入中的 Wannier 关键词 | 所有占据 Wannier 函数中心之和给出 $\mathbf P_\mathrm{el}=-\frac{fe}{\Omega}\sum_n\bar{\mathbf r}_n$，与 Berry 相等价（模量子），见[第14章](14-wannier.md) |
| calB / Z2FH | scfout | Berry 曲率、Chern 数、Z2，见[第12章](12-soc-topology.md) |

除本章示例自带的脚本外，OpenMX 部分的通用脚本也可配合使用：[omx_parse_out.py](../scripts/omx_parse_out.py) 汇总 `.out` 中的能量与收敛信息，[omx_cube_tools.py](../scripts/omx_cube_tools.py) 做 cube 的差分与平面平均，[job_omx_slurm.sh](../scripts/job_omx_slurm.sh) 是通用 Slurm 作业模板。各脚本都可以用 `-h` 查看用法（如 `python3 omx_parse_out.py -h`）。

---

## 参数选择与收敛建议

| 项目 | 建议起点 | 收敛判据 / 说明 |
|---|---|---|
| 体系类型 | 体相绝缘体 | 金属的宏观极化无定义；小带隙体系需更密的 polB k 网格 |
| SCF `scf.Kgrid` | 与总能收敛测试相同 | ΔP 变化 < 0.5 μC/cm²，Z\* 变化 < 0.01 |
| polB k 网格 | 每个方向 6–9；极化方向 8–12 | 沿极化方向的分割数影响最大 |
| `scf.criterion` | 1e-9～1e-10 Hartree | 极化差是小量 |
| 电子温度 | ≤ 300 K | 不能出现分数占据 |
| Born 电荷位移 | ±0.02～0.05 Å，中心差分 | 两种位移结果一致；声学求和规则误差 < 0.05 |
| λ 路径 | 5 个点起步，畸变大时 9 个以上 | 相邻点 ΔP 远小于量子 |
| 外电场（分子） | 线性响应 1–3 GV/m，可视化 5–10 GV/m | ±F 结果对称；无电子溢出 |
| 盒子（场方向） | 分子两侧各 ≥ 8 Å 真空 | 跳变 eEL 不宜过大；α 随盒子增大不变 |
| 基组 | DFT_DATA19 Standard | 极化对极化函数（d、f）较敏感，收敛测试时升到 Precise 档 |
| 版本 | 4.0 + patch 4.0.1 | 后处理程序与 openmx 同一版本编译；3.9 的 k 网格为 Regular，结果略有不同 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `bash: polB: command not found` | 没有编译 polB，或可执行文件不在 PATH 中 | `cd source && make polB`；用 `find ~/openmx4.0 -name polB -type f` 找到后写完整路径（作业脚本中的 `POLB=` 变量） |
| 目录中没有 `*.scfout`，polB 无法读入 | SCF 输入里没开 `HS.fileout on`；或 `System.Name` 与预期不同 | 打开 `HS.fileout`；scfout 文件名是 `System.Name.scfout`，不是输入文件名 |
| polB 读入后崩溃或打印的晶格、原子数明显不对 | polB 与 openmx 版本或编译选项不一致 | 用同一 makefile、同一版本重新编译两者 |
| 批处理时 polB 不停打印 `The number of grids should be one or more` | `in` 第一行的 k 分割数含 0 或负数，程序反复要求重新输入 | 立即取消作业；`in` 必须恰好两行：三个 ≥ 1 的整数 + 三个 0/1 开关（3.2 节） |
| 相邻两个结构的 P 相差几十到上百 μC/cm² | 分支跳变（差了整数个极化量子） | 沿 λ 路径加密取点，用 `polar_path.py` 展开；比较时只看差值 |
| 中心对称结构的 P 既不是 0 也不是半量子 | 结构并不严格中心对称（坐标位数不够、原子顺序错）；SCF 未收敛 | 检查坐标；收紧 `scf.criterion`；加密 polB k 网格 |
| Z\* 的 ±u 单侧斜率不同 | 位移过大进入非线性区；某点分支跳变；SCF 未收敛 | 减小位移；检查各点的 Electron 列是否连续；复查 SCF |
| 声学求和规则严重不满足 | k 网格不足、基组不足，或各原子计算的设置不一致 | 统一设置后加密 k 网格、升级基组 |
| polB 打印 `Metallic !!` 和 `##!! Warning!! The system may be metallic! Result is not correct!` | 某些 k 点上 `ChemP` 以下的本征态数与占据带数不符：体系是金属，或电子温度过高、带隙过小 | 确认带隙（`.out` 中的本征值或 DOS）；降低电子温度；金属的宏观极化无定义 |
| 金属或小带隙体系的 P 随 k 网格剧烈变化 | 极化对金属无定义；分数占据（共线自旋极化时 polB 不会给出金属警告） | 确认两个自旋通道都有带隙；降低电子温度 |
| 体相加 `scf.Electric.Field` 后 SCF 振荡不收敛 | 锯齿势的不连续面切过材料 | 改用有真空的 slab、ESM `on4`，或用 polB 做零场计算 |
| 加场后偶极异常大、`Utot` 骤降 | 场强 × 盒长太大，电子溢出到真空低势一侧 | 减小场强或盒子长度；画 Δρ 检查真空中是否有电子 |
| 标准输出中没有 `<Applied External Electric Field>` | 关键词拼写错误、三个分量都为 0，或 `level.of.stdout 0`（该段只在 ≥ 1 时打印） | 核对 `scf.Electric.Field`（大小写不敏感，但拼写要准确） |
| 场的方向与预期不符 | 三个分量是沿 a、b、c 晶格矢量的分量，不是笛卡尔分量 | 非正交晶胞中按晶格矢量换算；看标准输出打印的方向单位矢量 |
| 求得的极化率为负 | 正的关键词值对应 −x 方向的物理电场 | 取绝对值，并以偶极变化方向为准记住该约定 |
| `<ESM:Warning> Non-zero ESM.ElectricField is not valid except for on4.`，随后程序退出 | 在 ESM on1/on2/on3 下设置了非零 `ESM.potential.diff`（4.0 `Input_std.c` 检查后直接 `exit`） | 改用 `on4` 或删除该行（第16章） |
| `diff_gcube` 打印 `Found a difference in the number of grid on a-axis`、`Found a difference in x-coordinate of the origin` 等 | 两次计算网格不同（结构、晶胞、截断能或网格原点不同）；程序不会停止，写出的差分无效 | 场响应只改场强；结构不同时按第10章用 `scf.fixed.grid` 固定网格原点后重算 |

---

## 练习

1. 用 `15-nacl-polb-born` 计算 Na 的 Z\*，再把 `polB.in` 第 1 行改为 `6 6 6` 和 `12 12 12`，画出 Z\* 随 k 网格的收敛曲线。
2. 修改 `make_disp.sh`，对 Cl 做 ±0.05 Å 位移，验证声学求和规则 $Z^*(\mathrm{Na})+Z^*(\mathrm{Cl})\approx0$。
3. 在 BaTiO3 示例中，用弛豫后的结构计算 Ti 沿 z 的 $Z^*_{zz}$，与立方相的 7.16 比较，并解释为什么四方相中会变小。
4. 把 BaTiO3 的 λ 路径从 5 个点减到 2 个点（只算 λ = 0 和 1），用 `polar_path.py` 观察是否还能确定正确分支；再用 `--factor 2` 看量子格点的变化。
5. 用硝基苯示例的 ±1、±2、±3 GV/m 三组场强分别求 α_xx，判断线性区范围；再把盒子 a 从 24 Å 改为 30 Å，检查 α 是否变化。
6. 把场改加在 y 方向（`0.0 F 0.0`），求 α_yy，与 α_xx 比较分子极化率的各向异性。
7. 选做：用 Z\* 与第 5 节弛豫得到的位移，按 $\Delta P_z\approx\frac{e}{\Omega}\sum_I Z^*_{I,zz}\Delta z_I$ 估算 BaTiO3 的自发极化，与 Berry 相结果对比。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 44 章 Berry 相宏观极化（polB）：<https://www.openmx-square.org/openmx_man4.0/s44_berry-polar.html>
- 第 30 章 外加电场（锯齿势、硝基苯示例）：<https://www.openmx-square.org/openmx_man4.0/s30_electric-field.html>
- 第 66 章 两个 cube 文件之差（`diff_gcube`）：<https://www.openmx-square.org/openmx_man4.0/s66_cube-diff.html>
- 第 49 章 ESM（含 on4）：<https://www.openmx-square.org/openmx_man4.0/s49_esm.html>
- 第 70.1 节 scfout 文件格式：<https://www.openmx-square.org/openmx_man4.0/s70_1_scfout.html>
- 第 58 章 Chern 数与 Berry 曲率（calB）：<https://www.openmx-square.org/openmx_man4.0/s58_chern-berry-fhs.html>
- 关键词总表：<https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html>
- 官方输入示例列表（`NaCl.dat`、`Nitro_Benzene.dat`）：<https://www.openmx-square.org/openmx_man4.0/s78_input-examples.html>
- 推荐基组表（DFT_DATA19）：<https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html>
- What's new（4.0 更新列表、3.9 中 polB 的改进）：<https://www.openmx-square.org/whatsnew/whatsnew.html>

**相关 GitHub 仓库（只链接，不复制）**

- [Ncmexp2717/HWC](https://github.com/Ncmexp2717/HWC)：读 scfout 计算混合 Wannier 中心与分层极化，可作为 Berry 相极化的补充
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：从 OpenMX 哈密顿量构建紧束缚模型，计算 Berry 曲率、光学与非线性响应
- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：官方下载页链接的 4.0 GPU 版，`source/` 中含 `polB.c`、`diff_gcube.c` 等源码，本章用来核对 polB 的输出格式与金属判据
- [hikaruri/Z2-support-OpenMX](https://github.com/hikaruri/Z2-support-OpenMX)：OpenMX 自带 Z2FH 的使用说明，Berry 联络/曲率类后处理的另一参考
- [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-)：印尼语教程，含 Berry 曲率、反常霍尔等分支，可作功能对照
- 更多仓库见 [../../Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**

- R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993)，Berry 相极化理论：<https://doi.org/10.1103/PhysRevB.47.1651>
- D. Vanderbilt, R. D. King-Smith, *Phys. Rev. B* **48**, 4442 (1993)，极化作为体相性质及其与表面电荷的关系：<https://doi.org/10.1103/PhysRevB.48.4442>
- R. Resta, *Rev. Mod. Phys.* **66**, 899 (1994)，宏观极化综述：<https://doi.org/10.1103/RevModPhys.66.899>
- N. A. Spaldin, *J. Solid State Chem.* **195**, 2 (2012)，现代极化理论入门：<https://doi.org/10.1016/j.jssc.2012.05.010>
- W. Zhong, R. D. King-Smith, D. Vanderbilt, *Phys. Rev. Lett.* **72**, 3618 (1994)，钙钛矿的反常 Born 电荷与巨 LO–TO 劈裂：<https://doi.org/10.1103/PhysRevLett.72.3618>
- Ph. Ghosez, J.-P. Michenaud, X. Gonze, *Phys. Rev. B* **58**, 6224 (1998)，ABO3 化合物的动力学电荷（含 BaTiO3 Z\* 对照表）：<https://doi.org/10.1103/PhysRevB.58.6224>
- H. H. Wieder, *Phys. Rev.* **99**, 1161 (1955)，BaTiO3 单晶自发极化的实验测量：<https://doi.org/10.1103/PhysRev.99.1161>
- I. Souza, J. Íñiguez, D. Vanderbilt, *Phys. Rev. Lett.* **89**, 117602 (2002)，体相有限电场的 Berry 相方法（对照）：<https://doi.org/10.1103/PhysRevLett.89.117602>
- P. Umari, A. Pasquarello, *Phys. Rev. Lett.* **89**, 157602 (2002)，体相有限电场方法（对照）：<https://doi.org/10.1103/PhysRevLett.89.157602>
- M. Otani, O. Sugino, *Phys. Rev. B* **73**, 115407 (2006)，ESM 方法：<https://doi.org/10.1103/PhysRevB.73.115407>

相关章节：[第5章 SCF 与 k 点](05-scf-kpoints-convergence.md)、[第6章 结构优化](06-optimization-eos-neb.md)、[第9章 光学性质](09-dos-pdos-cohp-optics.md)、[第10章 cube 与 diff_gcube](10-charge-mo-cube.md)、[第12章 Berry 曲率与拓扑](12-soc-topology.md)、[第16章 ESM](16-surface-esm-workfunction.md)。

[← 上一章](14-wannier.md) | [返回目录](../README.md) | [下一章 →](16-surface-esm-workfunction.md)
