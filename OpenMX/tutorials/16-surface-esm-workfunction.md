# 第16章 表面与界面：ESM、功函数、电荷掺杂（Surfaces & Interfaces in OpenMX: ESM, Work Function & Charge Doping）

本章讲在 OpenMX 中如何搭建表面板模型（slab），并用有效屏蔽介质法（Effective Screening Medium，ESM）处理真空、带电电极和外加电场这三类边界条件。在此基础上，演示怎样从 `v0.cube` 的平面平均势计算金属表面功函数（work function），以及怎样用 `scf.system.charge` 给表面或二维材料做电子/空穴掺杂（charge doping）。配套示例包括 Al(111) 功函数（ESM on1 + 空原子层）、带电 Pt(111) 电极（ESM on3/on4）和石墨烯掺杂；最后简要介绍 STM 模拟和 OpenMX 4.0 新增的 XPS + ESM 组合计算。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 按 OpenMX 的约定搭建 slab：表面法向必须是**第一条晶格矢量（a 轴）**；会选层数、真空厚度和 k 点，会用 `<MD.Fixed.XYZ` 固定底层；
- 理解 `ESM.switch` 的四种模式（on1/on2/on3/on4）分别对应的物理边界，知道 4.0 源码对 ESM 几何的检查规则以及那些不会报错的坑；
- 用 `gcube2oned` 或 Python 求平面平均势，得到 $\Phi=\bar V_\mathrm{vac}-\mu$，并做空原子层、层数、k 点的收敛测试；
- 处理非对称和极性表面的偶极：OpenMX 没有 VASP `LDIPOL` 那样的偶极修正关键词，要用 ESM；
- 用 `scf.system.charge` 掺杂，并区分“均匀背景补偿”和“ESM 理想金属补偿”两种物理图像；
- 知道 STM（Tersoff–Hamann）、XPS 芯能级与 ESM 联用、带势垒墙的 ESM-MD 分别在哪里查、怎么起步。

| 任务 | 关键设置 | 本章示例 |
|---|---|---|
| 金属表面功函数 | `ESM.switch on1` + 空原子层 + `v0.cube` 平面平均 | [16-al111-esm-workfunction](../examples/16-al111-esm-workfunction/README.md) |
| 带电电极 / 电容器 | `ESM.switch on3` + `scf.system.charge` | [16-pt111-esm-charging](../examples/16-pt111-esm-charging/README.md) |
| 均匀外电场（双金属板） | `ESM.switch on4` + `ESM.potential.diff` | 同上，`Pt111_ESM_on4.dat` |
| 表面弛豫（固定底层） | `MD.Type RF` + `<MD.Fixed.XYZ` | 同上，`Pt111_relax_fixbottom.dat` |
| 二维材料掺杂 / 栅压 | `scf.system.charge`，背景补偿或 `on3` + `ESM.direction z` | [16-graphene-charge-doping](../examples/16-graphene-charge-doping/README.md) |

前置知识：输入文件结构见[第3章](03-input-first-calc.md)，基组与空原子见[第4章](04-pao-vps.md)，SCF 与 k 点见[第5章](05-scf-kpoints-convergence.md)，结构优化见[第6章](06-optimization-eos-neb.md)，cube 文件与差分密度见[第10章](10-charge-mo-cube.md)。锯齿势外电场 `scf.Electric.Field` 见第15章。通用的表面建模方法（超胞、切面、异质结）见 Common 第3章。

---

## 2. 原理速览

### 2.1 为什么 slab 需要特殊的静电处理

OpenMX 在实空间网格上用 FFT 解 Poisson 方程，因此在三个方向上都默认周期性。把一块板放进带真空的超胞，相当于沿法向无限重复“板–真空–板–真空……”。会出现两个问题：

1. **偶极相互作用**：如果板两侧不对称（单面吸附、极性表面、外加电荷或电场），板会带一个法向偶极矩 $p_\perp$（单位面积 $p_\perp/A$）。周期边界要求势在超胞两端相等，结果真空中被强加了一个人为的反向电场。两侧真空能级之差的物理值为

   $$\Delta V=\frac{e\,p_\perp}{\varepsilon_0 A}$$

   如果不处理，这个差值就会被周期性“抹平”，功函数和吸附能都会出错。
2. **带电体系**：净电荷 $q\neq0$ 时，周期体系的静电能发散，FFT 求解器只能引入均匀背景电荷来中和（见 2.4 节）。

平面波程序通常用偶极修正（dipole correction，Neugebauer–Scheffler / Bengtsson 方法）来处理第一个问题。OpenMX 4.0 的 `Input_std.c` 中没有这类关键词，对应的工具是 ESM。

### 2.2 ESM 方法

ESM（Otani & Sugino, 2006）在 slab 法向上放置两个**半无限介质**，并用 Green 函数法求解带介电函数的 Poisson 方程（原子单位）：

$$\nabla\cdot\left[\varepsilon(\mathbf r)\,\nabla V(\mathbf r)\right]=-4\pi\rho(\mathbf r)$$

面内（b、c）方向保持周期性，法向变成开放边界。真空介质取 $\varepsilon=1$，理想金属取 $\varepsilon=\infty$（介质内电势为常数，可以容纳镜像电荷）。OpenMX 把两个介质放在晶胞边界 $x=0$ 和 $x=|\mathbf a|$ 处，组合如下：

| `ESM.switch` | ESM(I) @ x=0 | ESM(II) @ x=\|a\| | 物理图像 | 允许带电？ |
|---|---|---|---|---|
| `off` | — | — | 普通 3D 周期 | 用背景电荷 |
| `on1` | 真空 | 真空 | 孤立板（功函数、偶极板） | **否**，`scf.system.charge` 必须为 0 |
| `on2` | 理想金属 | 理想金属 | 两块金属之间的板 | 是 |
| `on3` | 真空 | 理想金属 | 电极表面 + 对电极（电化学界面、栅极） | 是 |
| `on4` | 理想金属 | 理想金属 + 电势差 | 平行板电容器中的板（均匀外场） | 是 |

“板”不一定是金属表面，也可以是二维材料或二维排列的分子层。ESM 还能与 MD 联用，用来模拟固/液界面，这时需要势垒墙把液体分子限制在指定区域，见第 9.3 节。

### 2.3 功函数

功函数定义为把一个电子从费米能级移到表面外“无穷远”所需的最小能量。借助 Janak 定理，可以写成（OpenMX 手册第 50 章式 (58)）：

$$\Phi=\left(\phi_\infty+E[N-1]\right)-E[N]=\phi_\infty-\mu$$

其中 $\phi_\infty$ 是真空中远离表面处的势，$\mu$ 是化学势（`.out` 中的 `Chemical potential (Hartree)`）。在实际计算中，$\phi_\infty$ 取 KS 势在 bc 面上平面平均后、真空平台区的值：

$$\bar V(x)=\frac{1}{A}\iint_A V_\mathrm{KS}(x,y,z)\,\mathrm dy\,\mathrm dz$$

所用的势文件是 `System.Name.v0.cube`：自旋向上的 KS 势，已扣除非局域赝势部分；`scf.ProExpn.VNA on`（默认）时还扣除了中性原子势 VNA。这两部分都是短程的，在真空区为零，所以不影响真空能级，但**板内部**的 `v0` 并不是完整势，不要拿它解释体内能带排列。两个量的单位都是 Hartree：$\Phi(\mathrm{eV})=27.2114\times(\bar V_\mathrm{vac}-\mu)$。

对有带隙的体系，$\mu$ 落在带隙中，位置受电子温度影响，式 (58) 只能作近似（手册原话）。更稳妥的做法是分别给出电离势 $\mathrm{IP}=\bar V_\mathrm{vac}-E_\mathrm{VBM}$ 和电子亲和能 $\mathrm{EA}=\bar V_\mathrm{vac}-E_\mathrm{CBM}$。

### 2.4 电荷掺杂的两种补偿方式

`scf.system.charge q` 表示从体系中**移走** q 个电子：正值为空穴掺杂（体系带正电），负值为电子掺杂，可以取分数。电荷由谁来补偿，物理图像完全不同：

| 方式 | 设置 | 补偿电荷 | 适合 | 注意 |
|---|---|---|---|---|
| 均匀背景 | `ESM.switch off` | 充满整个晶胞（包括真空）的均匀反号背景 | 体相掺杂、粗略的二维掺杂 | 总能量含背景静电项，不同电荷态的能量不能直接比较；真空区的势会弯曲 |
| ESM 理想金属 | `on2`/`on3`/`on4` | 金属介质中的镜像电荷 | 电极充电、栅压、电化学界面 | 体系等价于电容器，过剩电荷会偏聚到朝向金属的一侧 |
| 孤立带电体系 | `scf.coulomb.cutoff on` | 无（精确 Coulomb 截断） | 分子的 IP/EA（手册第 61 章） | 截断半径为最短晶格矢量的一半，盒子要足够大 |

面电荷密度与电场的换算（A 为表面胞面积，单位 Å$^2$）：

$$\sigma\,[\mu\mathrm{C/cm^2}]=1602.18\times\frac{q}{A},\qquad E=\frac{\sigma}{\varepsilon_0}\;\Rightarrow\;q=0.01,\ A=6.65\ \text{Å}^2\ \text{时}\ E\approx2.7\ \mathrm{GV/m}$$

---

## 3. 板模型构建

### 3.1 OpenMX 中 ESM 的几何约定（务必先读）

手册第 49.1 节和 4.0 源码（`Input_std.c`）共同规定：

1. **a 轴（`<Atoms.UnitVectors` 第一行）就是表面法向**，而且必须与某个笛卡尔轴平行，分量为正。`ESM.direction x|y|z`（默认 x）告诉程序 a 轴沿哪个笛卡尔方向；它**不会**把别的晶格矢量当作法向。
2. b、c 必须与 a 正交，也就是说在 ESM 方向上的分量为 0。
3. ESM 介质位于 a 轴的两端 0 与 |a|。原子沿 a 轴的分数坐标必须落在 0–1 之间（用 FRAC 时程序会检查）。
4. ESM 计算必须显式给出 `<Atoms.UnitVectors`，不能让程序自动定胞。
5. k 点网格的第一个数对应 a 轴，所以 slab 总写成 `scf.Kgrid 1 n n`。

违反第 1、2 条时，程序会打印下面的信息并退出：

```text
<ESM:Warning>
The a-axis of the unit cell used for the ESM calculation
must be parallel to surface normal axis,
where its component must be POSITIVE.
```

```text
<ESM:Warning>
The b- and c-axes of the unit cell used for the ESM calculation
must be orthogonal to the a-axis.
```

两种常见写法如下：

| 写法 | `ESM.direction` | `<Atoms.UnitVectors` | 例子 |
|---|---|---|---|
| 法向沿 x（默认） | `x` | `L 0 0` / `0 b_y b_z` / `0 c_y c_z` | 手册 Al111_WorkFunc_*.dat；本章 Al、Pt 示例 |
| 法向沿 z（多数建模软件的习惯） | `z` | `0 0 L` / `b_x b_y 0` / `c_x c_y 0` | 手册 N1Be2_IS.dat（XPS+ESM）；本章石墨烯示例 |

> [!IMPORTANT]
> 从 VASP/ASE/OpenMX Viewer 导出的 slab，通常是第三条矢量 c 沿 z 作法向。要用 ESM，必须把晶格矢量**重排**成 (c, a, b)（循环置换，晶胞手性不变），并设 `ESM.direction z`。坐标若用 Ang 写，不必改动；若用 FRAC，分数坐标的列也要做同样的置换。只设 `ESM.direction z` 而不重排矢量，会触发上面的 `<ESM:Warning>`。

### 3.2 层数、真空与 k 点

| 参数 | 推荐起点 | 收敛判据 / 说明 |
|---|---|---|
| 金属层数 | fcc(111) 6–8 层；fcc(100) 7–9 层 | 功函数变化 < 0.02 eV；中心层要像体相 |
| 半导体层数 | 足以让中心层恢复体相 DOS | 极性面需要钝化或 ESM |
| 真空（ESM on1） | 每侧 ≥ 10 Å（原子到边界） | ESM 下势收敛很快，真空不必像周期计算那样取 20–30 Å；但基函数尾巴不能碰到边界 |
| 原子到 ESM 边界的距离 | ≥ 基函数截断半径 + 余量（Al7.0 约 3.7 Å，建议 ≥ 8–10 Å） | 4.0 源码**不检查**，见第 4 节 |
| 面内 k 点 | 金属 1×1 胞约 1×15×15 以上；2×2 胞相应减半 | 用功函数或表面能收敛 |
| 法向 k 点 | 1 | — |
| `scf.energycutoff` | 220–300 Ry，显式写出 | 平面平均势的平滑程度 |

slab 是“厚金属 + 大真空”，SCF 容易出现电荷振荡（charge sloshing）。经验做法：`RMM-DIISK`，`scf.Init.Mixing.Weight` 0.01–0.05，`scf.Max.Mixing.Weight` 0.05–0.1，`scf.Mixing.StartPulay` 15–25，`scf.Mixing.History` 30。手册的 Al(111) 示例用的就是 0.01 量级的混合权重。

### 3.3 用 ASE 生成符合约定的坐标

示例目录里的 [make_slab_omx.py](../examples/16-al111-esm-workfunction/make_slab_omx.py) 用 ASE 的 `fcc111`/`fcc100`/`bcc110`/`hcp0001` 生成板，做循环置换让法向成为 a 轴，把板居中，并可在两侧追加空原子层：

```bash
cd OpenMX/examples/16-al111-esm-workfunction
python3 make_slab_omx.py --surface fcc111 --element Al --a 4.05 \
        --layers 6 --length 50 --empty 2 --direction x > al111_blocks.txt
python3 make_slab_omx.py -h        # all options
```

输出的第二行会给出最外层原子到两个 ESM 边界的距离，例如：

```text
# distance of outermost atoms to the ESM boundaries: 14.48 / 14.48 Ang
```

不想写代码的话，可以用 OpenMX Viewer（2026-05 版）的 **Slab Builder**：按 Miller 指数、板厚、真空厚度和终止面生成 slab，再导出 OpenMX 输入。导出后务必按 3.1 节检查法向是否为第一条矢量。其他建模工具（cif2omx、cif2x、acoord 等）见[第3章](03-input-first-calc.md)和 Common 第3章。

### 3.4 固定底层做表面弛豫

表面弛豫的常用做法是固定底部 1–2 层，让它们模拟体相，只弛豫上部各层。`<MD.Fixed.XYZ` 的每一行依次是原子序号和 x、y、z 三个方向的标志，1 表示固定，0 表示放开（第6章）：

```text
MD.Type                       RF
MD.maxIter                    100
MD.Opt.criterion              3.0e-4      # Hartree/Bohr
<MD.Fixed.XYZ
   1   1 1 1
   2   1 1 1
   3   0 0 0
   4   0 0 0
MD.Fixed.XYZ>
```

完整输入见 [Pt111_relax_fixbottom.dat](../examples/16-pt111-esm-charging/Pt111_relax_fixbottom.dat)。ESM 下也能正常计算力和做优化。如果要算吸附，建议采用“对称板”（两面都吸附）或者“单面吸附 + ESM on1”，否则偶极会污染能量。

> [!WARNING]
> `<...>` 块里的行**不能写注释**，也不能插入注释行。OpenMX 按位置读取块中的每一行，例如物种块用 `sscanf("%s %s %s %lf")` 读取，第 4 列被当作原子质量；坐标行多出来的列会被解释为非共线角度等参数。注释请写在块外。

### 3.5 空原子层：补全真空一侧的基组

PAO 基函数在截断半径之外严格为零。表面外侧几 Å 内的电子密度尾巴（溢出电子，spill-out）决定了表面偶极层，而这部分只能靠表面原子的轨道来描述，往往不够。办法是在表面外放置**空原子**（empty atom）：使用零赝势 `E.vps`、不带电子，只提供基函数。

```text
Species.Number                2
<Definition.of.Atomic.Species
  Al   Al7.0-s2p2d1   Al_PBE19
  E    Al7.0-s2p2d1   E
Definition.of.Atomic.Species>
```

坐标行写 `0.0 0.0`，位置一般放在“延续晶格”的格点上，即把 slab 向外多堆几层。手册第 50 章的测试显示，Al(111) 两侧各加 2 层空原子后功函数已收敛。DFT_DATA19 推荐表中给 `E` 列出的基组是 Kr10.0-s1p1/s2p1d1/s2p2d1f1（Quick/Standard/Precise），而手册的 Al 示例直接把 Al 的 PAO 放在空原子上。空原子用什么基组、加几层，本身就是需要做收敛测试的参数。

---

## 4. ESM 关键词详解（按 4.0 源码核对）

| 关键词 | 取值 / 默认 | 说明 |
|---|---|---|
| `ESM.switch` | `off`\|`on1`\|`on2`\|`on3`\|`on4`，默认 `off` | 见 2.2 节表格 |
| `ESM.direction` | `x`\|`y`\|`z`，默认 `x` | a 轴（法向）所沿的笛卡尔方向 |
| `ESM.potential.diff` | eV，默认 0.0 | 仅 `on4`：两块理想金属之间的电势差，以底部金属为参考；电场 ≈ 电势差 / a 轴长度。在其他模式下设成非零值会报错退出 |
| `ESM.wall.switch` | `on`\|`off`，默认 `off` | ESM-MD 的立方势垒墙开关（手册正文没有单独介绍，但官方示例 Al-Si111_ESM.dat 用到了它） |
| `ESM.wall.position` | Å；手册写默认 10.0，**源码默认 1.0** | 势垒原点到晶胞上边缘的距离 $a-x_b$。用到时务必显式设置 |
| `ESM.wall.height` | eV，默认 100.0 | 势垒高度，定义为 $x_b+1$ Å 处的值 |
| `ESM.buffer.range` | Å，默认 10.0 | 手册说明为“原子与 ESM 之间的禁区宽度” |

> [!CAUTION]
> 在 4.0 源码中，`ESM.buffer.range` 和 `ESM.wall.height` 只有在 `ESM.wall.switch on` 时才会被读取，而检查“原子离边界太近”的代码块处于 `#if 0` 停用状态。也就是说，原子贴着 ESM 边界时**程序不会报错**，只会悄悄给出错误结果。要自己保证所有原子（包括空原子）离 x=0 和 x=|a| 有足够距离：至少大于基函数截断半径，本章示例都在 12 Å 以上。另外，关键词拼错同样不会报错，而是被静默忽略。例如官方示例 Al-Si111_ESM.dat 里写的 `ESM.potential.difference` 并不是源码读取的 `ESM.potential.diff`。好在那个示例中它的值本来就是 0，没有造成影响。

确认 ESM 已经生效的方法是看标准输出中的横幅：

```bash
grep -A1 "The following calc. is implemented with ESM method" x.std
#     The following calc. is implemented with ESM method.
#     Boundary condition = Vacuum|Vacuum|Vacuum           (on1)
#     Boundary condition = Metal|Vacuum|Metal             (on2; on4 adds "plus Uniform electric field")
#     Boundary condition = Vacuum|Vacuum|Metal            (on3)
```

---

## 5. 实战一：Al(111) 功函数（ESM on1）

示例目录：[16-al111-esm-workfunction](../examples/16-al111-esm-workfunction/README.md)。

### 5.1 输入文件

下面是 `Al111_ESM_on1_2E.dat` 的完整内容：6 层 Al，两侧各 2 层空原子，a 轴长 50 Å，板居中。

```text
System.CurrrentDirectory      ./          # official spelling: three 'r'
System.Name                   al111_esm_on1_2e
DATA.PATH                     /path/to/openmx4.0/DFT_DATA19   # EDIT
level.of.stdout               1
level.of.fileout              1           # >= 1 is needed to write *.v0.cube

Species.Number                2
<Definition.of.Atomic.Species
  Al   Al7.0-s2p2d1   Al_PBE19
  E    Al7.0-s2p2d1   E
Definition.of.Atomic.Species>

Atoms.Number                  10
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
     1  E      14.477791     0.000000     1.653406    0.0  0.0
     2  E      16.816060     1.431891     0.826703    0.0  0.0
     3  Al     19.154329     0.000000     0.000000    1.5  1.5
     4  Al     21.492597     0.000000     1.653406    1.5  1.5
     5  Al     23.830866     1.431891     0.826703    1.5  1.5
     6  Al     26.169134     0.000000     0.000000    1.5  1.5
     7  Al     28.507403     0.000000     1.653406    1.5  1.5
     8  Al     30.845671     1.431891     0.826703    1.5  1.5
     9  E      33.183940     0.000000     0.000000    0.0  0.0
    10  E      35.522209     0.000000     1.653406    0.0  0.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit        Ang
<Atoms.UnitVectors
     50.000000     0.000000     0.000000
      0.000000     2.863782     0.000000
      0.000000    -1.431891     2.480108
Atoms.UnitVectors>

scf.XcType                    GGA-PBE
scf.SpinPolarization          off
scf.ElectronicTemperature     300.0
scf.energycutoff              250.0
scf.maxIter                   300
scf.EigenvalueSolver          Band
scf.Kgrid                     1 15 15
scf.Generation.Kpoint         GCenter     # 4.0 default; Regular = 3.9 behaviour
scf.Mixing.Type               RMM-DIISK
scf.Init.Mixing.Weight        0.02
scf.Min.Mixing.Weight         0.001
scf.Max.Mixing.Weight         0.05
scf.Mixing.History            30
scf.Mixing.StartPulay         20
scf.criterion                 1.0e-8

ESM.switch                    on1
ESM.direction                 x

MD.Type                       Nomd
```

| 关键词 | 本例取值 | 为什么 |
|---|---|---|
| 面内晶格 | b=(0, 2.8638, 0)，c=(0, −1.4319, 2.4801) | 实验 a₀=4.05 Å，$a_0/\sqrt2$，夹角 120° |
| 层序 | A(0,0)→B(1/3,2/3)→C(2/3,1/3)（以 b、c 分数计） | fcc(111) 的 ABC 堆垛，层间距 $a_0/\sqrt3$ |
| `level.of.fileout` | 1 | 0 不输出 cube 文件，就没有 `v0.cube` 可用 |
| `scf.Kgrid` | 1 15 15 | a 轴取 1；1×1 金属表面胞需要较密的面内网格 |
| `scf.ElectronicTemperature` | 300 K | 金属；k 点足够密时不必取得很高 |
| 混合参数 | 0.02 / 0.05，StartPulay 20 | 抑制大真空 slab 的电荷振荡 |
| `scf.criterion` | 1e-8 | 功函数取决于真空能级与 μ 之差，需要较紧的收敛 |
| `ESM.switch` | on1 | 中性孤立板，真空\|真空 |

### 5.2 运行与检查

```bash
cd OpenMX/examples/16-al111-esm-workfunction
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat
mpirun -np 8 openmx Al111_ESM_on1_2E.dat -nt 1 > al111_esm_on1_2e.std 2>&1

grep "Boundary condition" al111_esm_on1_2e.std
grep "The calculation was normally finished." al111_esm_on1_2e.std
grep "Required cutoff energy\|Num. of grids" al111_esm_on1_2e.out
grep "Chemical potential (Hartree)" al111_esm_on1_2e.out
ls al111_esm_on1_2e.v0.cube
```

SCF 是否收敛，看 `.out` 中 `SCF history` 一段最后几行的 `NormRD` 是否持续减小、相邻 `Uele` 之差是否小于 `scf.criterion`。集群上的完整作业脚本见示例目录的 `job_slurm.sh`。

### 5.3 平面平均：gcube2oned

`gcube2oned` 是 OpenMX 自带的小工具：把 cube 数据在其余两个方向上平均，输出沿指定轴（1/2/3 分别对应 a/b/c）的一维曲线。

```bash
# 4.0 源码包中 gcube2oned.c 在 source/ 目录（手册写的是 work/）
gcc ~/openmx4.0/source/gcube2oned.c -lm -o gcube2oned      # or: cd source && make gcube2oned
./gcube2oned al111_esm_on1_2e.v0.cube 1 > 1d_pot.txt
head -3 1d_pot.txt
```

输出三列：网格序号、沿该轴的位置（Å，从 cube 原点量起）、平面平均势（Hartree）。

> [!NOTE]
> 读源码可以发现，轴参数取 2 或 3 时，`gcube2oned` 的第一列打印的是一个固定值而不是序号，所以只应使用第 2、3 列。另外第二列是从 cube 原点量起的距离，不是绝对坐标。

### 5.4 计算功函数

用示例中的 [workfunc.py](../examples/16-al111-esm-workfunction/workfunc.py) 自动选取两侧真空窗口（默认为晶胞的 5%–15% 和 85%–95%），并分别给出左右两侧的功函数：

```bash
python3 workfunc.py 1d_pot.txt al111_esm_on1_2e.out --plot al111_wf.png
python3 workfunc.py 1d_pot.txt al111_esm_on1_2e.out --left 3 9 --right 41 47   # custom windows
```

```text
# cell length along the averaging axis :    50.0000 Ang
# chemical potential mu               : -0.15xxxxxx Hartree = -4.1xxx eV
left  vacuum [  2.50,  7.50] Ang: V = ...  ->  Phi =   4.1xx eV
right vacuum [ 42.50, 47.50] Ang: V = ...  ->  Phi =   4.1xx eV
vacuum-level difference V_right - V_left =   0.000 eV (non-zero for asymmetric/polar slabs)
```

（示意输出，数值仅供参考）

想自己写脚本的话，核心只有几行：读入 cube，在 b、c 两个方向求平均，取真空平台减去 μ。这一步也可以交给 OpenMX 部分的公共脚本 omx_cube_tools.py，它直接对 cube 做平面平均，不依赖 `gcube2oned`，用法见 `python ../scripts/omx_cube_tools.py -h`：

```python
import numpy as np
def read_cube(fn):
    with open(fn) as f:
        f.readline(); f.readline()
        nat = abs(int(f.readline().split()[0]))
        n = [int(f.readline().split()[0]) for _ in range(3)]
        for _ in range(nat):
            f.readline()
        data = np.array(f.read().split(), dtype=float)
    return data.reshape(n)               # cube data: outer loop over a, inner loop over c
v = read_cube("al111_esm_on1_2e.v0.cube")
vbar = v.mean(axis=(1, 2))               # average over b and c -> profile along a (Hartree)
```

**结果判断**：

- 真空平台平坦（波动 < 1 meV），左右功函数相同（对称板）；
- 手册第 50 章报告：用 7 层空原子时 Al(111) 计算值为 4.19 eV，实验值为 4.26±0.03 eV（Eastment & Mee 1973）。本例的层数、基组和截断能与手册不同，结果相差零点几 eV 以内都算合理；
- 再跑一遍 `Al111_ESM_on1.dat`（无空原子）对比：两者之差就是真空侧基组不完备造成的误差。

### 5.5 收敛测试清单

| 参数 | 测试序列 | 判据 |
|---|---|---|
| 空原子层数（每侧） | 0, 1, 2, 3 | ΔΦ < 0.02 eV（手册结论：Al(111) 2 层即可） |
| Al 层数 | 6, 8, 10 | ΔΦ < 0.02–0.03 eV（量子尺寸效应会让 Φ 随层数振荡） |
| 面内 k 点 | 1×12×12, 1×15×15, 1×18×18 | ΔΦ < 0.01 eV |
| `scf.energycutoff` | 200, 250, 300 Ry | ΔΦ < 0.01 eV |
| 基组 | Al7.0-s2p1d1 → s2p2d1 → s3p2d2 | 与空原子层的测试一起做 |
| 电子温度 | 100, 300, 600 K | μ 的漂移 |

---

## 6. 非对称与极性表面

单面吸附（例如 Pt(111) 顶面吸附 O）、两种终止面不同的板、极性面（如 ZnO(0001)/(000-1)），两侧真空能级不同。ESM on1 在两侧都是真空、法向开放的条件下求解 Poisson 方程，偶极造成的真空能级台阶 $\Delta V=e p_\perp/(\varepsilon_0A)$ 会自然出现，**不需要**额外的偶极修正：

1. 用 `on1` 计算；
2. 用 `gcube2oned … 1` 得到 $\bar V(x)$，会看到两个高度不同的平台；
3. `workfunc.py` 分别给出 $\Phi_\mathrm{left}$ 和 $\Phi_\mathrm{right}$，并打印 $\Delta V$。

功函数的变化 $\Delta\Phi$ 应取吸附一侧的功函数减去清洁表面的值。由 $\Delta V$ 还可以反推出表面偶极（单位面积）：

$$\frac{p_\perp}{A}=\frac{\varepsilon_0\,\Delta V}{e}$$

> [!TIP]
> 对极性半导体板（手册第 50 章也提醒“有带隙的极性表面必须用 ESM”），除了 ESM 之外，还应该钝化底面（例如用赝氢原子）或者构造对称板，否则内建电场可能把带隙“压塌”，导致非物理的电荷转移。判断方法：看板内部的平面平均势是否在一个常数附近周期振荡，而不是带有整体斜坡。

---

## 7. 实战二：带电 Pt(111) 电极（ESM on3）与电容器（on4）

示例目录：[16-pt111-esm-charging](../examples/16-pt111-esm-charging/README.md)。

### 7.1 on3 输入要点

`Pt111_ESM_on3.dat` 使用 4 层 Pt(111)（实验 a₀ = 3.92 Å，1×1 胞），板位于 x = 12.00–18.79 Å，a 轴长 40 Å。x=0 一侧是真空，x=40 Å 一侧是理想金属。关键部分如下：

```text
Species.Number                1
<Definition.of.Atomic.Species
  Pt   Pt7.0-s3p2d2f1   Pt_PBE19
Definition.of.Atomic.Species>
...
scf.ElectronicTemperature     500.0
scf.energycutoff              250.0
scf.EigenvalueSolver          Band
scf.Kgrid                     1 15 15
scf.system.charge             0.000       # >0 hole doping (slab positive), <0 electron doping
ESM.switch                    on3
ESM.direction                 x
```

`Pt_PBE19` 有 16 个价电子，坐标行写 `8.0 8.0`。完整文件见示例目录。

### 7.2 电荷扫描与微分电容

```bash
cd OpenMX/examples/16-pt111-esm-charging
NP=16 NT=1 bash run_charge_scan.sh          # q = -0.02 ... +0.02
cat charge_scan.csv
```

脚本对每个 q 生成 `pt111_on3_q{m|p}0.0xx.dat`，运行后从 `.out` 中提取 `Chemical potential (Hartree)` 和 `Utot.`，换算出 σ 写入 CSV。分析思路：

- **μ–σ 关系**：在 on3 下，板与理想金属之间是真空电容器。μ 随 σ 近似线性变化，斜率给出单位面积微分电容的倒数 $1/C=\mathrm d\phi/\mathrm d\sigma$；按平行板电容粗估，$C\approx\varepsilon_0/d_\mathrm{eff}$（$d_\mathrm{eff}$ 约为表面到 ESM(II) 的距离，本例约 21 Å，对应约 0.4 μF/cm²，这是数量级估算）；
- **过剩电荷分布**：用 `diff_gcube`（源码 `source/diff_gcube.c`，用法为 `diff_gcube input1.cube input2.cube output.cube`）对带电与中性的 `tden.cube` 作差，再用 `gcube2oned` 做平面平均。过剩电荷会集中在朝向金属的那一面，这与手册第 49.2 节 Al–Si(111) 的结果相同；
- **能量比较**：不同电荷态的 `Utot` 对应的电子数不同，不能直接相减。若要讨论恒电势下的稳定性，应使用巨正则能 $\Omega=E-\mu_e N$，其中 $\mu_e$ 是电子库（电极电势）的化学势，相关方法见 ESM 文献（Otani & Sugino 2006；Sugino *et al.* 2007）。

### 7.3 on4：两块理想金属之间的均匀电场

```text
ESM.switch                    on4
ESM.direction                 x
ESM.potential.diff            1.0         # eV; only valid with on4
```

两块金属分别位于 x=0 和 x=|a|，电势差为 1.0 V（以底部金属为参考）。真空区的电场约为 ΔV / 有效间距，板内部的电场被金属屏蔽。检查方法是对 `pt111_on4_dv1.v0.cube` 做平面平均，真空区应当出现线性斜坡。

on4 与 `scf.Electric.Field`（锯齿势，第15章）的区别：

| | `scf.Electric.Field` | ESM `on4` |
|---|---|---|
| 实现 | 周期锯齿势，在晶胞中某处有不连续面 | 两侧理想金属边界，势连续 |
| 体相无真空 | 不稳定（不连续面落在材料内） | 不适用（本来就是板） |
| 单位 | GV/m（符号按作用于电子定义） | eV（电势差） |
| 带电体系 | 用背景电荷 | 可以 |
| 典型用途 | 分子、有真空的板的介电响应 | 电容器中的板、电极/栅极 |

---

## 8. 实战三：石墨烯电荷掺杂

示例目录：[16-graphene-charge-doping](../examples/16-graphene-charge-doping/README.md)。三个输入的几何相同：a 轴是沿 z 的法向（20 Å），b、c 在 xy 面内、夹角 120°，用 FRAC 坐标：

```text
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
     1  C     0.500000   0.333333   0.666667    2.0  2.0
     2  C     0.500000   0.666667   0.333333    2.0  2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit        Ang
<Atoms.UnitVectors
      0.000000     0.000000    20.000000
      2.467000     0.000000     0.000000
     -1.233500     2.136485     0.000000
Atoms.UnitVectors>
scf.Kgrid                     1 24 24
```

三个输入的差别只在电荷与 ESM 设置：

| 文件 | `scf.system.charge` | ESM |
|---|---|---|
| `Graphene_q0.dat` | 0.0 | `off` |
| `Graphene_qm0.02_bg.dat` | −0.02 | `off`（背景电荷） |
| `Graphene_qm0.02_esm_on3.dat` | −0.02 | `on3` + `ESM.direction z` |

三个计算都输出 Γ–M–K–Γ 能带，K=(0, 1/3, 1/3)（按 a、b、c 倒格矢计）。用示例里的 `dirac_point.py` 读取 `.Band`，求 K 点狄拉克点相对 μ 的位置：

```bash
python3 dirac_point.py graphene_q0.Band graphene_qm002_bg.Band graphene_qm002_on3.Band
```

`.Band` 中的本征值和 μ 都以 Hartree 为单位、未平移；`bandgnu13` 输出的 `GNUBAND`/`BANDDAT1` 已换算为 eV，并以 μ 为零点（第8章）。也可以用 OpenMX 部分的公共脚本 omx_band_plot.py 画图（`python ../scripts/omx_band_plot.py -h`）。

**如何判断结果**：中性时 $E_D-\mu\approx0$；电子掺杂 0.02 e/胞（$n\approx3.8\times10^{13}$ cm$^{-2}$）时，按线性色散 $E_F-E_D=\hbar v_F\sqrt{\pi n}$ 估算，费米能级应比狄拉克点高 0.6–0.7 eV（取 $v_F\approx(0.85\text{–}1.0)\times10^6$ m/s 的估算值）。背景电荷和 ESM 两种做法下，能带相对 μ 的位移应当接近，区别在静电势：ESM on3 下石墨烯与“栅极”之间有线性电势坡（栅电场），背景电荷做法里真空中的势则呈抛物线弯曲。

> [!NOTE]
> 背景电荷模型下，掺杂后的总能量和真空能级都受背景电荷影响，不适合用来比较不同电荷态的能量，也不适合计算掺杂后的功函数。需要这些量时，请改用 ESM `on2`/`on3`。

---

## 9. 进阶

### 9.1 STM 图像（Tersoff–Hamann）

STM 模拟就是计算以 μ 为参考的能量窗口内的部分电荷密度（手册第 53 章）。它在 DOS 计算中顺带完成，所以必须同时打开 DOS：

```text
partial.charge                 on         # default off
partial.charge.energy.window   -1.0       # eV from mu; <0 occupied (negative sample bias), >0 empty states
Dos.fileout                    on
Dos.Erange                     -20.0 20.0
Dos.Kgrid                      1 24 24
```

运行后得到 `System.Name.pden.cube`，可以在 VESTA 里画等值面（恒流模式），也可以取板上方某一高度的截面（恒高模式）。手册的示例是 work/Graphene_STM.dat（能量窗口 2 eV，图像用 WSxM 显示）。表面外侧的部分电荷密度同样受 3.5 节所说的基组尾巴问题影响，STM 计算也建议加空原子层。

### 9.2 XPS 芯能级结合能 + ESM（4.0 新增）

4.0 把绝对芯能级结合能计算（Ozaki & Lee, PRL 2017，基于惩罚泛函和精确 Coulomb 截断）扩展到可以与 ESM 联用，适合表面和二维体系（手册第 60 章）：

- 初态和终态是两个计算：终态使用 `vps_pao_core2019` 数据库中带芯空穴的 VPS/PAO（例如 `N6.0_1s-s3p2d1` + `N_PBE19_1s`）；
- ESM 只作用于**初态**电子密度，所以即使终态带电（半导体的终态取 `scf.system.charge 1.0`），也用 `ESM.switch on1`；
- 官方示例 work/N1Be2_IS.dat、N1Be2_FS.dat 采用 `ESM.direction z`，a 轴为 (0, 0, 10) Å。手册表 14 显示：不用 ESM 时，结合能随 a 轴长度（10→100 Å）缓慢收敛；用 ESM 后，在 10 Å 时就已收敛。

具体的初态/终态输入写法和结合能公式，请对照手册第 60 章与 [vps_pao_core2019](https://www.openmx-square.org/vps_pao_core2019/) 数据库页。

### 9.3 ESM + MD：固/液界面与势垒墙

固/液界面 MD 的布置是：板贴近一侧，液体放在板与另一侧 ESM 之间，再用立方势垒墙防止分子跑进 ESM 区域（手册第 49.1 节第 5 点）：

```text
ESM.switch                    on3
ESM.wall.switch               on
ESM.wall.position             6.0         # Ang from the upper cell edge (set explicitly: source default 1.0)
ESM.wall.height               100.0       # eV at x_b + 1 Ang
```

手册建议 MD 中固定板底部的原子（`<MD.Fixed.XYZ`）。MD 的系综与控温方法见[第7章](07-md.md)。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 备注 |
|---|---|---|
| 法向 | 第一条矢量，`ESM.direction` 与之一致 | 源码强制检查 |
| 原子到 ESM 边界 | ≥ 8–10 Å（≥ 基组截断半径 + 余量） | 源码不检查，靠自己保证 |
| 金属 slab 层数 | 6–10 层，做收敛 | 功函数会随层数振荡 |
| 空原子层 | 每侧 2 层起步 | Al(111) 2 层即可（手册） |
| 面内 k 点 | 1×1 金属胞 ≥ 15×15 | 第一维写 1 |
| `scf.energycutoff` | 220–300 Ry | 显式写出 |
| 混合 | RMM-DIISK，权重 0.01–0.05 | 大真空 slab 防振荡 |
| `scf.criterion` | 1e-8（功函数）、1e-7（一般） | |
| 带电体系 | `on2`/`on3`/`on4`，或背景电荷 | `on1` 必须中性 |
| 功函数 | `v0.cube` 平面平均，真空平台减去 μ | 有带隙体系改用 IP/EA |
| 非对称板 | ESM on1，两侧分别读平台 | OpenMX 无偶极修正关键词 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `<ESM:Warning> The a-axis of the unit cell used for the ESM calculation must be parallel to surface normal axis, where its component must be POSITIVE.` | 法向不是第一条矢量，或第一条矢量与 `ESM.direction` 指定的轴不平行，或分量为负 | 把法向矢量放在 `<Atoms.UnitVectors` 第一行，`ESM.direction` 与之对应 |
| `<ESM:Warning> The b- and c-axes ... must be orthogonal to the a-axis.` | b 或 c 在法向上有分量（斜切晶胞） | 构造正交于法向的表面胞 |
| `<ESM:Warning> A unit cell must be specified in the ESM calculation.` | 没写 `<Atoms.UnitVectors` | 写出晶胞 |
| `<ESM:Warning> The fractional coordinate of a-axis should be kept within 0 to 1 in the ESM calculation.` | FRAC 坐标超出 0–1 | 平移到晶胞内（板居中） |
| `<ESM:Warning> Non-zero ESM.ElectricField is not valid except for on4.` | 在 on1/on2/on3 下设了非零 `ESM.potential.diff` | 改用 on4 或删掉该行 |
| 功函数明显偏离、两侧平台不平 | 真空太小、原子离 ESM 边界太近、SCF 未收敛 | 加大 a 轴长度；检查 `NormRD`；把平台窗口移到远离板的位置 |
| 加空原子前后功函数差很多 | 真空侧基组不完备 | 每侧加 2–3 层空原子并测试收敛 |
| SCF 在 slab 中来回振荡 | 大真空 + 金属造成电荷振荡 | 降低混合权重，加大 `scf.Mixing.StartPulay`，适当提高电子温度 |
| 带电计算没有按预期补偿 | 用了 on1，或者 `ESM.switch` 拼错（被静默忽略） | 用 on2/on3；`grep "Boundary condition" x.std` 确认 |
| `Found an undefined species name E ...` | 用了空原子 `E`，但物种块里没有定义 | 在 `<Definition.of.Atomic.Species` 中加一行 `E  <PAO>  E` |
| `Could not find .../VPS/E.vps` | `DATA.PATH` 错误 | 指向含 PAO/ 和 VPS/ 的 DFT_DATA19 |
| 块内一行写了注释后，物种/坐标读错 | 块行按位置解析 | 注释移到块外 |
| 掺杂后不同电荷态的能量差异常 | 背景电荷的静电能 | 改用 ESM（巨正则比较）或 `scf.coulomb.cutoff on`（孤立分子） |

---

## 练习

1. 把 `Al111_ESM_on1_2E.dat` 改成每侧 1 层和 3 层空原子（`make_slab_omx.py --empty 1/3`），画出 Φ 随空原子层数的曲线，与手册图 68 的结论对比。
2. 在 Pt(111) 顶面 fcc 空位上放一个 O 原子（单面吸附），用 on1 计算两侧功函数，由 $\Delta V$ 求表面偶极 $p_\perp/A$。
3. 用 `run_charge_scan.sh` 的结果画 μ–σ 曲线，求斜率得到微分电容，与 $\varepsilon_0/d_\mathrm{eff}$ 比较。然后把 a 轴从 40 Å 加长到 50 Å（板位置不变），观察电容如何变化。
4. 石墨烯掺杂扫描 q = −0.05…+0.05，画 $E_D-\mu$ 与 $\mathrm{sign}(n)\sqrt{|n|}$ 的关系，拟合费米速度 $v_F$。
5. 对 Al(111) 打开 `partial.charge`（窗口 −0.5 eV），用 VESTA 查看 `pden.cube` 中表面上方的电子密度，比较有无空原子层时的差别。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 49 章 ESM：<https://www.openmx-square.org/openmx_man4.0/s49_esm.html>（49.1 概述、49.2 Al–Si(111) 示例）
- 第 50 章 功函数：<https://www.openmx-square.org/openmx_man4.0/s50_work-function.html>
- 第 31 章 电荷掺杂：<https://www.openmx-square.org/openmx_man4.0/s31_charge-doping.html>
- 第 30 章 外电场：<https://www.openmx-square.org/openmx_man4.0/s30_electric-field.html>
- 第 12.5 节 空原子：<https://www.openmx-square.org/openmx_man4.0/s12_5_placement-empty-atoms.html>
- 第 18.4 节 约束优化 `MD.Fixed.XYZ`：<https://www.openmx-square.org/openmx_man4.0/s18_4_constraints.html>
- 第 53 章 STM：<https://www.openmx-square.org/openmx_man4.0/s53_stm.html>
- 第 60 章 XPS 芯能级（含 60.4 与 ESM 联用）：<https://www.openmx-square.org/openmx_man4.0/s60_xps-core-binding.html>
- 第 61 章 孤立体系 IP/EA 与 `scf.coulomb.cutoff`：<https://www.openmx-square.org/openmx_man4.0/s61_ip-ea.html>
- 第 9 章 输出文件（`v0.cube` 的定义）：<https://www.openmx-square.org/openmx_man4.0/s9_output-files.html>
- 关键词总表：<https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html>
- What's new（4.0 的 XPS+ESM、Viewer Slab Builder）：<https://www.openmx-square.org/whatsnew/whatsnew.html>
- OpenMX Viewer：<https://www.openmx-square.org/viewer/>
- 芯能级 VPS/PAO 数据库：<https://www.openmx-square.org/vps_pao_core2019/>

**相关 GitHub 仓库（只链接，不复制）**

- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站（安装、能带、DOS、+U、SOC），可作为基础操作的对照
- [wxia529/acoord](https://github.com/wxia529/acoord)：在 VS Code 内编辑、可视化 OpenMX `.dat` 结构，便于检查 slab 和空原子位置
- [issp-center-dev/cif2x](https://github.com/issp-center-dev/cif2x)：从 CIF 和模板批量生成 OpenMX 输入（基于 pymatgen）
- [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx)：POSCAR→`.dat`，可以按晶格长度自动设定 k 点；注意导出后要按 ESM 约定重排法向
- [mfukudaQED/FLPQ](https://github.com/mfukudaQED/FLPQ)：4.0 区域化学势（RCP）分析的配套程序，可用来讨论表面吸附位点

**关键文献**

- M. Otani, O. Sugino, *Phys. Rev. B* **73**, 115407 (2006)，ESM 方法：<https://doi.org/10.1103/PhysRevB.73.115407>
- O. Sugino, I. Hamada, M. Otani, Y. Morikawa, T. Ikeshoji, Y. Okamoto, *Surf. Sci.* **601**, 5237 (2007)，ESM 用于电极界面
- T. Ohwaki, M. Otani, T. Ikeshoji, T. Ozaki, *J. Chem. Phys.* **136**, 134101 (2012)，OpenMX 中的 ESM 实现与大规模 MD：<https://doi.org/10.1063/1.3698583>
- J. F. Janak, *Phys. Rev. B* **18**, 7165 (1978)，Janak 定理：<https://doi.org/10.1103/PhysRevB.18.7165>
- J. Tersoff, D. R. Hamann, *Phys. Rev. B* **31**, 805 (1985)，STM 理论：<https://doi.org/10.1103/PhysRevB.31.805>
- T. Ozaki, C.-C. Lee, *Phys. Rev. Lett.* **118**, 026401 (2017)，绝对芯能级结合能：<https://doi.org/10.1103/PhysRevLett.118.026401>
- J. Neugebauer, M. Scheffler, *Phys. Rev. B* **46**, 16067 (1992)，偶极修正（对照）：<https://doi.org/10.1103/PhysRevB.46.16067>
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999)，slab 偶极修正（对照）：<https://doi.org/10.1103/PhysRevB.59.12301>
- R. M. Eastment, C. H. B. Mee, *J. Phys. F: Met. Phys.* **3**, 1738 (1973)，Al(111) 功函数实验值（手册表 12 引用）

---

[← 上一章](13-dft-u.md) | [返回目录](../README.md) | [下一章 →](17-negf-transport.md)
