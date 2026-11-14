# 第13章 DFT+U（Hubbard U 与 Hund J）（DFT+U with Hubbard U and Hund's J in OpenMX）

本章讲解如何在 OpenMX 4.0 中做 DFT+U（LDA+U / GGA+U）计算：从 Dudarev 简化形式（只含有效 U）到含 Hund J 的一般形式，再到 sFLL/sAMF/cFLL/cAMF 四种双计数项和 onsite/full/dual 三种占据数算符。示例体系是 AFM-II 型 NiO 和 MnO，我们会对比 GGA 与 GGA+U 的带隙、Ni/Mn 磁矩和 d 轨道占据矩阵，写一个 U 值扫描脚本，并说明线性响应法确定 U 的思路。关键词：OpenMX DFT+U、Hubbard U、Hund J、`Hubbard.U.values`、NiO 反铁磁、transition-metal oxide band gap。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 打开 DFT+U：`scf.Hubbard.U on`，并且知道它**必须**配合 `scf.SpinPolarization on` 或 `nc`；
- 区分 `scf.DFTU.Type 1`（Dudarev 简化形式，只读 U）和 `scf.DFTU.Type 2`（一般形式，同时读 U 与 J，需要选 `scf.dc.Type`）；
- 按 PAO 基组定义逐项写出 `<Hubbard.U.values` 与 `<Hund.J.values`，知道为什么 U 一般只加在 `1d`（或 `1f`）上；
- 理解 `scf.Hubbard.Occupation onsite|full|dual` 的区别，以及为什么默认 `dual`；
- 用 `RMM-DIISH` 混合让 +U 计算稳定收敛，必要时用坐标行的轨道极化开关寻找轨道有序的基态；
- 读懂 `.std` 中的 `DFT+U Type and DC` 信息、`.out` 中的 `Uhub` 能量项与 `Occupation Number in LDA+U and Constraint DFT` 占据矩阵；
- 用 U 值扫描、实验数据或线性响应思路给出有依据的 U。

| 体系类型 | 典型例子 | 是否需要 +U | 备注 |
|---|---|---|---|
| 3d 过渡金属氧化物/氟化物（Mott 或电荷转移绝缘体） | NiO、MnO、CoO、FeO、Cr2O3 | 几乎总是需要 | GGA 带隙严重偏小，甚至给出金属 |
| 4f 稀土化合物 | CeO2 的 Ce 4f、NdCo5 的 Nd 4f | 常用 | U 加在 `1f` 上 |
| 过渡金属掺杂/缺陷的局域态 | ZnO:Co、TiO2 中的 Ti3+ 极化子 | 视问题而定 | 结论对 U 敏感，必须做 U 扫描 |
| 普通 sp 半导体、简单金属 | Si、GaAs、Al、Cu | 不需要 | 带隙问题应换更高阶方法，不是加 U |

配套示例：[../examples/13-dft-u-nio-mno/README.md](../examples/13-dft-u-nio-mno/README.md)。

> [!NOTE]
> 版本说明：Dudarev 简化形式在 OpenMX 3.8 已有；一般形式（`scf.DFTU.Type 2`、`scf.dc.Type`、`scf.Slater.Ratio`、`scf.Yukawa`）从 3.9 开始提供。4.0 沿用同一套关键词，对应手册第 40 章 [DFT+U method](https://www.openmx-square.org/openmx_man4.0/s40_dft-plus-u.html)。本章输入在 4.0 下编写；用 3.9 时只需删掉 `scf.Generation.Kpoint` 一行（3.9 不认识这个关键词，会忽略它，但 k 网格默认就是 Regular，见第 5 章）。

---

## 2. 原理速览

### 2.1 为什么需要 +U

LDA/GGA 的交换关联泛函存在自相互作用误差（也叫离域化误差）。对 3d、4f 这类局域轨道，这个误差会让电子倾向于“平摊”到多个轨道上，得到分数占据，带隙和轨道极化因此被低估。NiO 是教科书例子：实验上它是带隙约 4 eV 的反铁磁绝缘体（光电子能谱/逆光电子能谱给出约 4.3 eV，见 Sawatzky & Allen 1984），而 LSDA/GGA 只给出不到 1 eV 的小带隙，Ni 磁矩也偏小。

DFT+U 的思路是：选定一组局域轨道（OpenMX 中就是某个原子上某个径向函数的 d 或 f 赝原子轨道），构造它们的占据数矩阵 $n^{I\sigma}_{mm'}$，再加上一个“鼓励整数占据”的修正能量。

### 2.2 Dudarev 简化形式（`scf.DFTU.Type 1`）

Dudarev 等（1998）给出的旋转不变简化形式为

$$
E_{\mathrm{DFT}+U} = E_{\mathrm{DFT}} + \frac{U_{\mathrm{eff}}}{2}\sum_{I,\sigma}\mathrm{Tr}\left[\mathbf n^{I\sigma} - \mathbf n^{I\sigma}\mathbf n^{I\sigma}\right]
$$

对应的轨道相关势为

$$
V^{I\sigma}_{mm'} = U_{\mathrm{eff}}\left(\frac{1}{2}\delta_{mm'} - n^{I\sigma}_{mm'}\right)
$$

占据数接近 1 的轨道被压低约 $U_{\mathrm{eff}}/2$，接近 0 的轨道被抬高约 $U_{\mathrm{eff}}/2$，因此 d 轨道的占据态和空态之间被拉开约 $U_{\mathrm{eff}}$。这里只出现一个参数，文献中常写成 $U_{\mathrm{eff}} = U - J$。OpenMX 的 `scf.DFTU.Type 1` 只读 `Hubbard.U.values` 中的数值，它就是这里的 $U_{\mathrm{eff}}$；手册说明此时隐含使用 sFLL 双计数。

### 2.3 一般形式：U、J 与 Slater 积分（`scf.DFTU.Type 2`）

一般形式（Liechtenstein 等 1995 的旋转不变形式及其推广）用完整的在位库仑相互作用张量 $\langle m_1 m_3|V_{ee}|m_2 m_4\rangle$ 构造修正能。对 d 壳层，这个张量由 Slater 积分 $F^0, F^2, F^4$ 和 Racah–Wigner 系数决定，输入参数 U、J 与 Slater 积分的关系是

$$
U = F^0,\qquad J = \frac{F^2 + F^4}{14},\qquad \frac{F^4}{F^2}\approx 0.625
$$

$F^4/F^2 = 0.625$ 是球对称原子环境下的近似值，OpenMX 用 `scf.Slater.Ratio` 调整它（默认 0.625）；也可以用 `scf.Yukawa on` 假设 Yukawa 型屏蔽库仑势，由输入的 U 自动估计 J 与 $F^4/F^2$（见 7.3 节）。

与简化形式相比，一般形式多了 J 的各向异性交换效应，对轨道有序、Hund 规则耦合、高低自旋竞争等问题更合适。

### 2.4 双计数项：FLL 与 AMF，s 与 c

DFT 本身已经以平均场方式包含了局域电子之间的部分相互作用，加 U 时必须减掉这部分，即双计数（double counting, DC）修正。OpenMX 4.0 提供四种组合：

| `scf.dc.Type` | 含义 | 适用的直观图像 |
|---|---|---|
| `sFLL`（默认） | spin-density DFT + fully localized limit | 局域性强、接近整数占据的体系；Liechtenstein 等的原始形式 |
| `sAMF` | spin-density DFT + around mean field | 占据数离整数较远、偏金属性的体系 |
| `cFLL` | charge-only DFT + FLL | 自旋极化完全交给 +U 部分处理 |
| `cAMF` | charge-only DFT + AMF | 同上，AMF 版本 |

两个字母的含义：

- **FLL / AMF**：FLL 假设局域电子完全局域、占据为整数时 DC 最合理；AMF 假设占据接近平均值（Czyżyk & Sawatzky 1994 讨论了二者差异）。以 FLL 为例，自旋版本（s）的双计数能量为 $E_{\mathrm{dc}} = \frac{U}{2}N(N-1) - \frac{J}{2}\sum_\sigma N^\sigma(N^\sigma-1)$；
- **s / c**：s 表示底层 DFT 用自旋密度泛函（LSDA/σGGA）；c 表示底层只用电荷密度泛函（LDA/GGA 的非自旋部分）。手册明确说明：选 `cFLL` 或 `cAMF` 时，SCF 循环中**不再计入**LDA/GGA 的自旋密度交换关联势，自旋劈裂全部来自 U、J 项。

Ryee & Han（2018）系统比较了这几种组合，发现 J 在 sDFT+U 与 cDFT+U 中的作用并不相同。OpenMX 手册的 NiO 图例（U = 5 eV，J = 0.5 或 1.0 eV）也显示：加 U 后带隙变大，而 J 的效果取决于 `scf.dc.Type`。所以使用 Type 2 时，**必须在文章中写明 DC 类型**，并且不要把 sFLL 与 cFLL 的结果混在一起比较。

### 2.5 占据数算符：`onsite` / `full` / `dual`

OpenMX 的基组是非正交的赝原子轨道（PAO），“某原子 d 轨道上有多少电子”没有唯一定义。Han、Ozaki、Yu（2006）讨论了三种占据数算符，可以粗略理解为：

| `scf.Hubbard.Occupation` | 直观理解 | 求和规则 |
|---|---|---|
| `onsite` | 只用密度矩阵的在位块，忽略与邻近原子轨道的重叠 | 不满足 |
| `full` | 用完整的投影算符，邻近原子的重叠贡献也投影进来 | 不满足 |
| `dual`（默认） | 用对偶基组构造，类似 Mulliken 布居的对称化形式 | **满足**：占据矩阵的迹之和等于总电子数 |

手册明确指出三者中只有 `dual` 满足求和规则，因此除非要复现某篇用 `onsite`/`full` 的文献，一律用默认的 `dual`。另外，非共线约束和 Zeeman 项只支持 `dual`；`full` 与应力计算不兼容（见常见问题）。

---

## 3. DFT+U 关键词总览

| 关键词 | 取值 / 默认值 | 说明 |
|---|---|---|
| `scf.Hubbard.U` | `on` \| `off`，默认 `off` | 总开关 |
| `scf.SpinPolarization` | 必须为 `on` 或 `nc` | `off` 时 +U 不生效（不要做非自旋极化 +U） |
| `scf.DFTU.Type` | `1` \| `2`，默认 `1` | 1 = Dudarev 简化形式（只读 U）；2 = 一般形式（读 U 和 J） |
| `scf.Hubbard.Occupation` | `onsite` \| `full` \| `dual`，默认 `dual` | 占据数算符 |
| `<Hubbard.U.values` | 块，单位 eV | 每个物种、每个径向轨道的 U（Type 1 时即 $U_{\mathrm{eff}}$） |
| `<Hund.J.values` | 块，单位 eV | 每个径向轨道的 J，仅 Type 2 读取（`scf.Yukawa on` 时忽略） |
| `scf.dc.Type` | `sFLL` \| `sAMF` \| `cFLL` \| `cAMF`，默认 `sFLL` | 仅 Type 2 有效；Type 1 隐含 sFLL |
| `scf.Slater.Ratio` | 默认 `0.625` | $F^4/F^2$，仅 Type 2 有效 |
| `scf.Yukawa` | `on` \| `off`，默认 `off` | 仅 Type 2：由 U 自动估计 J 和 $F^4/F^2$ |
| 坐标行最后一列 | `on` \| `off`，默认 `off` | 轨道极化促进开关：共线为第 8 列，非共线为第 13 列 |
| `scf.Mixing.Type` | 推荐 `RMM-DIISH` | 混合 KS 哈密顿量，手册表 4 推荐用于 +U 与约束计算 |

---

## 4. 写对 `Hubbard.U.values`：与 PAO 定义逐项对应

这是 OpenMX DFT+U 最容易出错的地方。规则只有一条：**`Definition.of.Atomic.Species` 中定义的每一个物种、每一个径向函数，都要按 s → p → d → f、多重度从 1 开始的顺序列出一个数值**。

例如本章使用的物种定义：

```text
<Definition.of.Atomic.Species
  Ni   Ni6.0S-s3p2d2   Ni_PBE19S
  O    O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>
```

`Ni6.0S-s3p2d2` 表示 3 组 s、2 组 p、2 组 d 径向函数，于是 Ni 行要写 `1s 2s 3s 1p 2p 1d 2d` 共 7 项；`O6.0-s2p2d1` 要写 `1s 2s 1p 2p 1d` 共 5 项：

```text
<Hubbard.U.values           # eV
  Ni  1s 0.0 2s 0.0 3s 0.0 1p 0.0 2p 0.0 1d 5.0 2d 0.0
  O   1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0
Hubbard.U.values>
```

更多对照：

| PAO 定义 | 需要写的项（按顺序） | 典型 U 位置 |
|---|---|---|
| `Ni6.0S-s3p2d1`（Standard） | `1s 2s 3s 1p 2p 1d` | `1d` |
| `Ni6.0S-s3p2d2` | `1s 2s 3s 1p 2p 1d 2d` | `1d` |
| `Ni6.0S-s2p2d2f1`（手册示例） | `1s 2s 1p 2p 1d 2d 1f` | `1d` |
| `Mn6.0-s3p2d2f1`（Precise） | `1s 2s 3s 1p 2p 1d 2d 1f` | `1d` |
| `O6.0-s2p2d1` | `1s 2s 1p 2p 1d` | 全 0（O 的 `1d` 是极化函数，不是 3d 电子） |
| `Nd8.0-s3p2d2f1` 一类 4f 元素 | `1s 2s 3s 1p 2p 1d 2d 1f` | `1f` |

> [!IMPORTANT]
> 阅读 4.0 源码 `Input_std.c` 可以看到：程序读取这个块时，对每个物种按 l = 0,1,2,3 和多重度循环，依次读入“一个标签 + 一个数值”，**标签本身只是占位，不参与匹配**。因此，项数或顺序写错时程序不一定报错，而是把数值错位地分配到别的轨道上。写完后请务必用第 6.3 节的方法核对占据矩阵和 `Uhub` 能量项。另外，没有 U 的物种（如 O、空原子 `E`、Wannier 投影用的辅助物种）**也必须列出**，全部写 0.0。

为什么只在 `1d` 上加 U？OpenMX 的优化 PAO 中，第一个 d 径向函数最接近原子的 3d 轨道，局域性最强；第二个 d 径向函数（`2d`）主要提供径向变分自由度，相当于极化/弛豫函数。给 `2d` 加 U 会让“局域轨道”的定义变得含糊，结果也难以和文献比较。4.0 的 `work/` 目录中，面向物理问题的 +U 示例（`NiO-cFLL.dat`、`NiO-sFLL.dat`、`CoO_NC.dat`、`FeO_NC.dat`、`MnO_NC.dat`、`SrVO3-Wannier90.dat`）都只在 `1d` 上加 U，`NdCo5_4f+U.dat` 加在 Nd 的 `1f` 上；`NiO-Yukawa.dat` 给 `2d` 加了较小的 U，用来演示 Yukawa 功能。`SrVO3-Wannier90.dat` 还表明，Wannier 投影用的辅助物种（`projV` 等）同样要在 U 块中列出。

基组较大时手写容易漏项，可以用下面的小函数生成（只处理 `s3p2d2f1` 这种不含收缩符号的写法）：

```python
import re

def hubbard_line(species, pao, values=None):
    """Return one line of <Hubbard.U.values for a PAO name such as 'Ni6.0S-s3p2d2'.
    values: dict like {'1d': 5.0}; all other orbitals get 0.0."""
    values = values or {}
    spec = pao.split("-")[-1]
    items = []
    for l, n in re.findall(r"([spdf])(\d+)", spec):
        for mul in range(1, int(n) + 1):
            key = f"{mul}{l}"
            items.append(f"{key} {values.get(key, 0.0):.2f}")
    return f"  {species:<4s}" + " ".join(items)

print(hubbard_line("Ni", "Ni6.0S-s3p2d2", {"1d": 5.0}))
print(hubbard_line("O", "O6.0-s2p2d1"))
```

输出：

```text
  Ni  1s 0.00 2s 0.00 3s 0.00 1p 0.00 2p 0.00 1d 5.00 2d 0.00
  O   1s 0.00 2s 0.00 1p 0.00 2p 0.00 1d 0.00
```

---

## 5. 实战一：NiO AFM-II 的 GGA 基准

### 5.1 结构：为什么用 4 原子菱方胞

NiO 是岩盐结构（实验晶格常数约 4.18 Å），奈尔温度以下形成 AFM-II 序：沿 [111] 方向，Ni 的 (111) 铁磁层一层自旋向上、一层自旋向下交替排列。描述这种序至少需要沿 [111] 加倍的菱方胞，含 2 个 Ni、2 个 O。取 $h = a/2$：

$$
\mathbf a_1 = (2h, h, h),\quad \mathbf a_2 = (h, 2h, h),\quad \mathbf a_3 = (h, h, 2h)
$$

Ni1 位于 (0,0,0)（满足 $x+y+z = 0$ 的 (111) 面），Ni2 位于 $(h,h,0)$（$x+y+z = 2h$，相邻的 (111) 面），O 位于 $(h,0,0)$ 与 $(h,h,h)$。每个 Ni 有 6 个次近邻 Ni（经 O 的 180° 超交换），其中全部与它反平行——这正是 AFM-II 能量最低的原因。

### 5.2 输入文件 `NiO_AFM_GGA.dat`

```text
#
# NiO AFM-II (rhombohedral 4-atom cell), GGA-PBE reference without +U
# Target: OpenMX 4.0 (3.9: delete scf.Generation.Kpoint)
#
System.CurrrentDirectory    ./
System.Name                 nio_gga
level.of.stdout             1
level.of.fileout            1
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: dir that contains PAO/ and VPS/

Species.Number      2
<Definition.of.Atomic.Species
  Ni   Ni6.0S-s3p2d2   Ni_PBE19S
  O    O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>

Atoms.Number        4
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1  Ni   0.00000   0.00000   0.00000   10.0   6.0
  2  Ni   2.08850   2.08850   0.00000    6.0  10.0
  3  O    2.08850   0.00000   0.00000    3.0   3.0
  4  O    2.08850   2.08850   2.08850    3.0   3.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit   Ang
<Atoms.UnitVectors
  4.17700   2.08850   2.08850
  2.08850   4.17700   2.08850
  2.08850   2.08850   4.17700
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.ElectronicTemperature   300.0
scf.energycutoff            250.0
scf.EigenvalueSolver        Band
scf.Kgrid                   7 7 7
scf.Generation.Kpoint       GCenter
scf.maxIter                 300
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.05
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15
scf.criterion               1.0e-8

MD.Type                     nomd

Dos.fileout                 on
Dos.Erange                 -10.0  10.0
Dos.Kgrid                   11 11 11
```

逐项说明：

| 关键词 | 本例取值 | 为什么 |
|---|---|---|
| `Ni6.0S-s3p2d2` + `Ni_PBE19S` | 价电子 16（3p⁶3d⁸4s²） | `S`（soft）赝势不含 3s 半芯；比 Standard（`s3p2d1`）多一组 d，给 3d 留径向弛豫自由度。想更省可改回 Standard，但 U 行要同步删去 `2d` |
| `O6.0-s2p2d1` + `O_PBE19` | 价电子 6 | 第 4 章推荐的 Standard |
| 上/下自旋初始电子 | Ni1 = 10.0/6.0，Ni2 = 6.0/10.0，O = 3.0/3.0 | 两列之和等于价电子数（16、6）；Ni1、Ni2 符号相反，给出 AFM 初猜 |
| `scf.SpinPolarization` | `on` | 共线自旋极化；+U 也要求它 |
| `scf.Kgrid` | `7 7 7` | 菱方胞 $\lvert\mathbf a_i\rvert\approx 5.12$ Å，7 × 5.12 ≈ 36 Å，绝缘体足够；正式计算请做收敛测试（第 5 章） |
| `scf.energycutoff` | 250 Ry | 4.0 手册对默认值的写法不一致，务必显式设置 |
| `scf.Mixing.Type` | `RMM-DIISK` | 无 +U 时最稳健的选择 |
| `scf.criterion` | 1e-8 Hartree | 比较不同 U 的磁矩/带隙时，宁可收得紧一点 |
| `Dos.*` | 开 | 后面用 PDOS 比较 Ni-3d 与 O-2p（第 11 节） |

### 5.3 运行与检查

```bash
cd OpenMX/examples/13-dft-u-nio-mno
export OMP_NUM_THREADS=1
mpirun -np 8 openmx NiO_AFM_GGA.dat -nt 1 | tee nio_gga.std

grep "Utot\." nio_gga.out                     # total energy (Hartree)
grep -A6 "Total spin moment" nio_gga.out | head -8   # Mulliken spin moments
python3 gap_from_out.py nio_gga.out           # gap on the SCF k-mesh + moments
```

判断 SCF 是否真的收敛：`.std` 最后几步的 `NormRD` 应降到 1e-8 以下、`dUele` 小于 `scf.criterion`，且没有达到 `scf.maxIter` 才停（第 5 章）。

Mulliken 表中 `Diff` 列是上下自旋电子数之差（单位 μB）。AFM-II 收敛正确时应看到：Ni1 与 Ni2 的 `Diff` 大小相等、符号相反；两个 O 的 `Diff` 为 0（对称性保证）；`Total spin moment` 接近 0。

```text
（示意输出，数值仅供参考）
  Total spin moment (muB)   0.000000000

                    Up spin      Down spin     Sum           Diff
      1   Ni      9.3xxxxxxx   7.9xxxxxxx   17.2xxxxxxx   1.3xxxxxxx
      2   Ni      7.9xxxxxxx   9.3xxxxxxx   17.2xxxxxxx  -1.3xxxxxxx
      3    O      ...
```

`gap_from_out.py` 是示例目录中的小脚本，读取 `.out` 里在 SCF k 网格上打印的全部本征值和化学势，给出（各自旋及总的）VBM、CBM 与带隙估计。它只用 SCF 网格上的点，因此只是估计；最终带隙请用第 8 章的能带或第 9 章的 DOS 确认。

> [!TIP]
> 若 GGA 收敛到了铁磁态（Ni1、Ni2 的 `Diff` 同号）或非磁态（`Diff` ≈ 0），先检查坐标行的初始占据是否写成了相反符号，再把 `scf.Init.Mixing.Weight` 降到 0.02、`scf.Mixing.StartPulay` 提到 20–30。

---

## 6. 实战二：Dudarev GGA+U（U = 5 eV）

### 6.1 与 GGA 输入的差异

`NiO_AFM_U_Dudarev.dat` 在 GGA 输入基础上只改动以下几处（完整文件见示例目录）：

```text
System.Name                 nio_dudarev_u5

#
# DFT+U part (simplified rotationally invariant form of Dudarev et al.)
#
scf.Hubbard.U               on
scf.DFTU.Type               1          # 1: simplified (Dudarev), only U is read
scf.Hubbard.Occupation      dual       # onsite|full|dual, default=dual
<Hubbard.U.values           # eV, U_eff on the first d radial function of Ni
  Ni  1s 0.0 2s 0.0 3s 0.0 1p 0.0 2p 0.0 1d 5.0 2d 0.0
  O   1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0
Hubbard.U.values>

<Atoms.SpeciesAndCoordinates
  1  Ni   0.00000   0.00000   0.00000   10.0   6.0  off
  2  Ni   2.08850   2.08850   0.00000    6.0  10.0  off
  3  O    2.08850   0.00000   0.00000    3.0   3.0  off
  4  O    2.08850   2.08850   2.08850    3.0   3.0  off
Atoms.SpeciesAndCoordinates>

scf.Mixing.Type             RMM-DIISH  # mixes the KS Hamiltonian; recommended for +U
scf.restart                 off
```

要点：

- 坐标行第 8 列 `off` 是轨道极化促进开关（第 8 节）。打开 +U 后程序才读取这一列；写了就要**每个原子都写**，不写时默认 `off`；
- `RMM-DIISH` 直接混合 KS 哈密顿量。+U 项依赖占据数矩阵，哈密顿量不再只是总密度的泛函，密度看起来已经变化很小时哈密顿量仍可能大幅振荡，混合哈密顿量本身更稳定（手册 16.1 节）；
- U = 5 eV 与手册 NiO 示例一致。文献中 NiO 的 $U_{\mathrm{eff}}$ 常在 5–6.5 eV 之间（例如 Wang、Maxisch、Ceder 2006 拟合氧化焓得到 Ni 的 U ≈ 6.4 eV，但那是 VASP PAW 投影下的值，不能直接照搬，见 10.3 节）。

### 6.2 运行

```bash
mpirun -np 8 openmx NiO_AFM_U_Dudarev.dat -nt 1 | tee nio_dudarev_u5.std
```

### 6.3 确认 +U 真的生效

第一步，看 SCF 开始前 `.std` 中打印的 DFT+U 类型：

```bash
grep -A3 "DFT+U Type and DC" nio_dudarev_u5.std
```

Type 1 时应看到（格式摘自 4.0 手册 40.1.2 节）：

```text
*******************************************************
                   DFT+U Type and DC
*******************************************************
              scf.DFTU.Type: 1(Simplified)
```

第二步，看能量分解中的 `Uhub`（LDA+U 能量）不为 0：

```bash
grep -E "Uhub\.|Utot\." nio_dudarev_u5.out
```

第三步，看占据数矩阵。`.out` 中标题为 `Occupation Number in LDA+U and Constraint DFT` 的部分，对每个原子、每个自旋给出占据数矩阵的本征值（`Individual` 一行）和本征矢（按 s、px、py、pz、d3z²−r²、dx²−y²、dxy、dxz、dyz … 的轨道顺序）：

```bash
grep -n "Occupation Number in LDA+U" nio_dudarev_u5.out
sed -n "$(grep -n 'Occupation Number in LDA+U' nio_dudarev_u5.out | head -1 | cut -d: -f1),+40p" nio_dudarev_u5.out
```

对 Ni²⁺（d⁸，八面体场中 $t_{2g}^6 e_g^2$，高自旋）：多数自旋的 5 个 d 本征值都应接近 1；少数自旋 3 个（$t_{2g}$）接近 1、2 个（$e_g$）明显小于 1。加 U 后，这种“接近整数”的分布比 GGA 更清楚。因为 `dual` 算符满足求和规则，`Sum` 一行给出该原子该自旋的总占据。

### 6.4 GGA 与 GGA+U 对比

```bash
python3 gap_from_out.py nio_gga.out nio_dudarev_u5.out
```

请把你自己的结果填进下表。“文献常见范围”只是帮助判断数量级的近似值，不是 OpenMX 的参考输出：

| 量 | GGA（PBE） | GGA+U（$U_{\mathrm{eff}}$ = 5 eV） | 文献常见范围（近似） | 实验 |
|---|---|---|---|---|
| 带隙 (eV) | 你的结果 | 你的结果 | PBE 通常 < 1；PBE+U（5–6 eV）约 3 | 约 4.3（Sawatzky & Allen 1984） |
| Ni 自旋磁矩 (μB) | 你的结果 | 你的结果 | PBE 约 1.2–1.4；PBE+U 约 1.6–1.8 | 约 1.6–1.9（中子衍射，例如 Cheetham & Hope 1983 给出 1.64） |
| 价带顶主要成分 | Ni-3d | O-2p 与 Ni-3d 混合 | — | 电荷转移型绝缘体 |

> [!NOTE]
> Mulliken 磁矩依赖布居分析方式，与 VASP 在 PAW 球内积分的磁矩不能逐位比较；比较时用同一程序、同一方法的相对变化。OpenMX 输出的自旋磁矩符号约定见第 11 章。

> [!WARNING]
> 不同 U 值对应不同的能量泛函。**不要**比较不同 U 下的 `Utot`，也不要用不同 U 的总能计算反应能或形成能。同一组比较（例如 AFM 与 FM 的能量差、缺陷形成能）必须使用同一个 U。

---

## 7. 实战三：一般形式 U + J 与 cFLL 双计数

### 7.1 输入差异

`NiO_AFM_UJ_cFLL.dat` 与 Dudarev 输入的差别：

```text
System.Name                 nio_cfll_u6_j1

scf.Hubbard.U               on
scf.DFTU.Type               2          # 2: general form, U and J are both read
scf.dc.Type                 cFLL       # sFLL|sAMF|cFLL|cAMF, default=sFLL
scf.Slater.Ratio            0.625      # F4/F2, default=0.625
scf.Yukawa                  off        # on: estimate J and F4/F2 from U
scf.Hubbard.Occupation      dual

<Hubbard.U.values           # eV
  Ni  1s 0.0 2s 0.0 3s 0.0 1p 0.0 2p 0.0 1d 6.0 2d 0.0
  O   1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0
Hubbard.U.values>

<Hund.J.values              # eV
  Ni  1s 0.0 2s 0.0 3s 0.0 1p 0.0 2p 0.0 1d 1.0 2d 0.0
  O   1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0
Hund.J.values>
```

- `<Hund.J.values` 的书写规则与 `<Hubbard.U.values` 完全相同，项数和顺序必须与 PAO 定义一致；
- 本例取 U = 6 eV、J = 1 eV，使 $U - J = 5$ eV 与上一节的 $U_{\mathrm{eff}}$ 相同，方便讨论“同样的 U−J，换成一般形式和 cFLL 后有何差别”。NiO 的 J 在文献中常取 0.9–1.0 eV 量级（例如 Anisimov 等 1991 使用 J ≈ 0.95 eV）；
- 运行后 `.std` 中应出现 `scf.DFTU.Type: 2(General)   scf.dc.Type: cFLL`。

```bash
mpirun -np 8 openmx NiO_AFM_UJ_cFLL.dat -nt 1 | tee nio_cfll_u6_j1.std
grep -A3 "DFT+U Type and DC" nio_cfll_u6_j1.std
python3 gap_from_out.py nio_gga.out nio_dudarev_u5.out nio_cfll_u6_j1.out
```

### 7.2 双计数怎么选

- 想和大量 VASP/QE 文献的 “GGA+U（Dudarev）” 结果对照：用 `scf.DFTU.Type 1`；
- 需要 J 的物理（Hund 耦合、轨道有序、自旋态竞争），且希望 J 的作用“干净”：可以考虑 `cFLL`。此时自旋劈裂全部由 U、J 产生，底层 GGA 不再提供自旋交换关联势；
- `sFLL` + J ≠ 0 时，J 的效应会和自旋密度泛函中已有的交换劈裂叠加，解读要小心（Ryee & Han 2018）；
- 占据明显非整数、偏金属的体系可以试 AMF 版本，但要和 FLL 对比后再下结论。

一个很好的自检：`scf.DFTU.Type 2` + `scf.dc.Type sFLL` + J = 0 在理论上与 Dudarev 形式（U_eff = U）是同一个泛函，两次计算的 `Utot`、带隙和磁矩应在数值精度内一致（练习 2）。如果差别明显，通常是 U/J 块项数写错。

### 7.3 调整 F4/F2 与用 Yukawa 势估计 J

`scf.Slater.Ratio` 只在 Type 2 下有效。固体中 $F^4/F^2$ 可能偏离 0.625，一般不需要改动；若要研究它的影响，可以在 0.6–0.8 之间扫描。

`scf.Yukawa on` 时，程序只读 `Hubbard.U.values`，忽略 `Hund.J.values` 与 `scf.Slater.Ratio`，通过估计与输入 U 对应的 Thomas–Fermi 屏蔽长度，自动生成 J 和 $F^4/F^2$，并在 SCF 前打印：

```bash
grep -A10 "Calculating Thomas-Fermi screening length" nio_yukawa.std
```

手册 40.2.2 节给出的输出格式是：

```text
  <species: Ni, angular momentum= 2, multiplicity number= 0>
   TF-screening-length lambda= ... 1/au
   Hubbard U= ... eV
   Hund J= ... eV
   Slater F0= ... eV
   Slater F2= ... eV
   Slater F4= ... eV
   F4/F2= ...
```

想尝试时，把 cFLL 输入复制一份，把 `scf.Yukawa` 改成 `on`、`System.Name` 改成 `nio_yukawa` 即可（官方 `work/NiO-Yukawa.dat` 也是这样做的）。注意 Yukawa 估计的是 **J 与 F4/F2**，并不能替你确定 U。

---

## 8. 实战四：MnO 与轨道极化开关

### 8.1 MnO AFM-II

MnO 与 NiO 同为岩盐结构、AFM-II 序（实验晶格常数约 4.445 Å）。Mn²⁺ 为 d⁵ 高自旋，多数自旋 5 个 d 轨道全满、少数自旋全空，没有轨道简并问题，是检验 +U 是否正确设置的好例子：磁矩应接近 5 μB（Mulliken 值因杂化略小），带隙随 U 明显增大。`MnO_AFM_U.dat` 的关键部分：

```text
<Definition.of.Atomic.Species
  Mn   Mn6.0-s3p2d2    Mn_PBE19
  O    O6.0-s2p2d1     O_PBE19
Definition.of.Atomic.Species>

<Hubbard.U.values           # eV
  Mn  1s 0.0 2s 0.0 3s 0.0 1p 0.0 2p 0.0 1d 4.0 2d 0.0
  O   1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0
Hubbard.U.values>

<Atoms.SpeciesAndCoordinates
  1  Mn   0.00000   0.00000   0.00000   10.0   5.0  off
  2  Mn   2.22250   2.22250   0.00000    5.0  10.0  off
  3  O    2.22250   0.00000   0.00000    3.0   3.0  off
  4  O    2.22250   2.22250   2.22250    3.0   3.0  off
Atoms.SpeciesAndCoordinates>
```

`Mn_PBE19` 的价电子数为 15（3s²3p⁶3d⁵4s²），初始占据 10.0/5.0 给出 5 μB 的初猜。实验数据可参考：带隙约 3.6–4 eV（van Elp 等 1991），Mn 磁矩约 4.6 μB（Cheetham & Hope 1983）。

### 8.2 轨道极化促进开关

DFT+U 泛函在轨道自由度上有多个极值，SCF 往往收敛到离初猜最近的局部解，不一定是基态。对 CoO（Co²⁺ d⁷）、FeO（Fe²⁺ d⁶）这类少数自旋 $t_{2g}$ 部分占据、存在轨道简并的体系，可能存在“有轨道极化”和“无轨道极化”两类解，二者带隙、磁各向异性都不同。

OpenMX 在坐标行的最后一列提供轨道极化促进开关：

```text
# collinear: column 8
<Atoms.SpeciesAndCoordinates
  1  Co   0.00000   0.00000   0.00000    9.0   6.0  on
  2  Co   2.13000   2.13000   0.00000    6.0   9.0  on
  3  O    2.13000   0.00000   0.00000    3.0   3.0  on
  4  O    2.13000   2.13000   2.13000    3.0   3.0  on
Atoms.SpeciesAndCoordinates>
```

（以上为说明格式的片段，晶格取 CoO 的近似值 a ≈ 4.26 Å；Co 的初始占据按 `Co_PBE19S` 的 15 个价电子设置，9.0/6.0 对应 Co²⁺ 高自旋约 3 μB 的初猜。）

- `on` 时，程序在最初几步 SCF 中主动促进轨道极化，之后转为普通 SCF；
- 必须对**所有**原子写 `on` 或 `off`；
- 手册提醒：这种处理常能找到基态，但也可能把收敛引向不稳定方向，**不是万能方法**；`scf.dc.Type cFLL/cAMF` 时不推荐使用；
- 正确做法是同时跑 `off` 和 `on` 两种初值，比较同一 U 下的 `Utot`，取较低者，并检查占据矩阵本征值是否出现明显的整数化（例如少数自旋 $t_{2g}$ 中两个接近 1、一个接近 0）。

MnO 是 d⁵ 半满壳层，没有轨道简并，本例保持 `off`；做 CoO/FeO 时再练习 `on`（官方 `work/CoO_NC.dat`、`work/FeO_NC.dat` 是非共线 + SOC + U 的例子，可作参考）。

---

## 9. 非共线 DFT+U（`scf.SpinPolarization nc`）

DFT+U 同样支持非共线计算（含 SOC）。此时占据数矩阵是包含自旋非对角元的 $2(2l+1)\times 2(2l+1)$ 矩阵，修正在自旋空间保持旋转不变。设置要点：

```text
scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on          # optional; needs j-dependent VPS (DFT_DATA19 is fine)
scf.Hubbard.U               on
scf.DFTU.Type               1
scf.Hubbard.Occupation      dual

<Atoms.SpeciesAndCoordinates
  1  Ni  0.00000  0.00000  0.00000  10.0  6.0  0.0   0.0   0.0   0.0   0      off
  2  Ni  2.08850  2.08850  0.00000   6.0 10.0  0.0   0.0   0.0   0.0   0      off
  3  O   2.08850  0.00000  0.00000   3.0  3.0  0.0   0.0   0.0   0.0   0      off
  4  O   2.08850  2.08850  2.08850   3.0  3.0  0.0   0.0   0.0   0.0   0      off
Atoms.SpeciesAndCoordinates>
```

非共线坐标行共 13 列：第 6–7 列为上/下自旋初始电子，第 8–9 列为自旋方向的 Euler 角 θ、φ（度），第 10–11 列为轨道磁矩方向的 θ、φ，第 12 列为非共线约束开关（0/1，配合 `scf.Constraint.NC.Spin`），**第 13 列为 +U 轨道极化开关**。各列的详细含义与约束 DFT 的用法见 [第11章](11-spin-noncollinear.md)，SOC 的设置见 [第12章](12-soc-topology.md)。

> [!TIP]
> nc + SOC + U 收敛很慢。可以先从较高的电子温度（如 1000 K）起步再降温重启，混合仍首选 `RMM-DIISH`。4.0 的 `scf.restart c2n` 能把共线计算的重启文件转换给非共线计算读取（手册第 17 章），它主要为第二变分法 MAE 设计；拿它给非共线 SCF 当初值时，请确认后续 SCF 确实充分收敛。

---

## 10. U 值扫描与 U 的确定

### 10.1 扫描 U：看带隙和磁矩怎么变

U 不是“算出来就对”的量，至少要知道结论对 U 有多敏感。示例目录的 `scan_U.sh` 以 `NiO_AFM_U_Dudarev.dat` 为模板，用 `sed` 改写 Ni 的 `1d` 值和 `System.Name`，依次计算并汇总：

```bash
NP=8 bash scan_U.sh "0 2 4 5 6 7 8"
# chain restarts (same System.Name, scf.restart on after the first U):
RESTART=1 NP=8 bash scan_U.sh "0 2 4 5 6 7 8"
```

脚本最后调用 `gap_from_out.py` 输出一张表（带隙、各 Ni 的 `Diff`、`Uhub`）。把带隙和磁矩对 U 作图：

```python
import csv
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("scan_U_summary.csv")))
U = [float(r["U"]) for r in rows]
gap = [float(r["gap_eV"]) for r in rows]
mom = [abs(float(r["moment_atom1"])) for r in rows]

fig, ax1 = plt.subplots(figsize=(4.5, 3.2))
ax1.plot(U, gap, "o-", label="gap")
ax1.set_xlabel("U_eff (eV)")
ax1.set_ylabel("band gap (eV)")
ax2 = ax1.twinx()
ax2.plot(U, mom, "s--", color="tab:red", label="|m(Ni1)|")
ax2.set_ylabel("Mulliken moment (muB)")
fig.tight_layout()
fig.savefig("nio_scan_U.png", dpi=200)
```

典型趋势（定性）：带隙随 U 近似线性增大，价带顶从 Ni-3d 为主逐渐转为 O-2p 为主；Ni 磁矩随 U 增大并向 2 μB 靠近。U = 0 且 `scf.Hubbard.U on` 时应与 GGA 结果一致，这是检查脚本和输入的好办法。

> [!NOTE]
> `RESTART=1` 模式在同一个 `System.Name` 下依次改 U，并用 `scf.restart on` 读入上一个 U 的收敛密度/哈密顿量，通常能显著减少 SCF 步数。但相邻 U 的解可能落在不同的局部极小，若发现某个 U 处磁矩或带隙突变，要用独立初值（非重启）复算该点。

### 10.2 线性响应法确定 U 的思路

Cococcioni 与 de Gironcoli（2005）提出用线性响应确定 U：在某个位点 I 的局域轨道上加一个小的势移 $\alpha_I$，测量占据数 $n_I$ 的变化，得到响应矩阵

$$
\chi_{IJ} = \frac{\partial n_I}{\partial \alpha_J}\ (\text{自洽}),\qquad
\chi^0_{IJ} = \frac{\partial n_I}{\partial \alpha_J}\ (\text{非自洽，第一步 KS 响应})
$$

$$
U = \left(\chi_0^{-1} - \chi^{-1}\right)_{II}
$$

其物理含义是：去掉“非相互作用”的 KS 响应后，剩下的就是电子–电子相互作用对局域占据的屏蔽效应。实际操作需要超胞（让被扰动位点彼此隔离）、多组 $\alpha$（例如 ±0.05、±0.1 eV）和线性拟合。

> [!IMPORTANT]
> OpenMX 4.0 手册中没有内置的线性响应 U 计算关键词。若要在 OpenMX 中实践上述思路，需要自己实现对局域轨道的势移扰动（或借助其他代码），这已超出本章范围。更重要的是：U 值与局域轨道（投影子）的定义绑定。OpenMX 用 PAO + `dual` 占据算符，VASP 用 PAW 投影球，QE 用原子波函数投影——同一体系的“最佳 U”在不同代码间不能直接照搬。

### 10.3 其他确定 U 的途径

| 方法 | 做法 | 优点 | 注意 |
|---|---|---|---|
| 文献经验值 | 采用同类体系、同类投影定义的 U | 便于与文献比较 | 必须注明来源和投影方式 |
| 拟合实验 | 让带隙、磁矩、晶格常数或氧化焓与实验一致（如 Wang、Maxisch、Ceder 2006 拟合氧化焓） | 有明确物理依据 | 不同性质拟合出的 U 可能不同 |
| 线性响应 / cRPA | 第一性原理估计 | 不依赖实验 | 依赖投影子定义，计算量大 |
| U 扫描 + 结论稳健性 | 在合理范围内（如 3–7 eV）报告结论是否改变 | 最诚实 | 文章中应给出扫描结果 |

---

## 11. 用 DOS/PDOS 验证 +U 效果

三个 NiO 输入都打开了 `Dos.fileout`，SCF 结束后会生成 `System.Name.Dos.val` 与 `System.Name.Dos.vec`。用 DosMain 的四面体法计算 Ni1 和 O3 的 PDOS（交互流程详见 [第9章](09-dos-pdos-cohp-optics.md)）：

```bash
# answers: 1 = tetrahedron, 2 = PDOS, then atom indices
DosMain nio_dudarev_u5.Dos.val nio_dudarev_u5.Dos.vec <<EOF
1
2
1 3
EOF
ls nio_dudarev_u5.PDOS.Tetrahedron.atom1.*
```

也可以用仓库脚本 [`omx_dos_auto.py`](../scripts/omx_dos_auto.py) 非交互地驱动 DosMain 并按元素和角动量求和（用 `python ../scripts/omx_dos_auto.py -h` 查看用法）。

看什么：

- GGA：价带顶附近以 Ni-3d（$e_g$ 多数自旋与 $t_{2g}$ 少数自旋）为主，带隙很小；
- GGA+U：占据的 Ni-3d 被压到价带深处（约 −6 ~ −8 eV），价带顶变为 O-2p 与 Ni-3d 的杂化态，导带底是少数自旋 Ni-$e_g$——这正是“电荷转移型绝缘体”的图像；
- 自旋极化时 PDOS 文件第 2、3 列分别为上/下自旋，Ni1 与 Ni2 的上下自旋应互换。

---

## 参数选择与收敛建议

| 参数 | 推荐做法 | 说明 |
|---|---|---|
| `scf.DFTU.Type` | 与文献对照用 1；需要 J 物理用 2 | 文章中写明 Type 与 DC |
| `scf.Hubbard.Occupation` | `dual`（默认） | 唯一满足求和规则 |
| U 的位置 | 只加在 `1d`（3d/4d/5d）或 `1f`（4f/5f） | 其他项写 0.0，但必须写全 |
| U 的数值 | 3d 氧化物常见 3–7 eV；做 U 扫描 | U 与投影定义绑定 |
| J | 3d 元素常见 0.7–1.0 eV | 只在 Type 2 下有意义 |
| `scf.Mixing.Type` | `RMM-DIISH` | `Init.Mixing.Weight` 0.02–0.1，`History` 20–30，`StartPulay` 10–30 |
| `scf.ElectronicTemperature` | 绝缘体 300 K；难收敛时先 1000 K 再降温重启 | 高温会抹平小带隙 |
| `scf.criterion` | 1e-8 Hartree（比较磁矩/能量差时） | 默认 1e-6 偏松 |
| k 点 / 截断能 | 先在 GGA 下收敛，再加 U | +U 一般不改变收敛参数的量级 |
| 初始磁构型 | 按目标磁序给对称/反对称初值 | 不同磁序要分别计算并比较 `Utot`（同一 U） |
| 轨道极化 | 简并体系同时试 `off` 与 `on` | 取同一 U 下能量较低的解 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决办法 |
|---|---|---|
| `Found an undefined species name ... in Atoms.SpeciesAndCoordinates, Hubbard.U.values, or Number.of.Electrons.Species` | U 块中的物种名与 `Definition.of.Atomic.Species` 不一致；或文件被 Windows 编辑器存成 UTF-16/带 BOM/含全角空格 | 核对物种名（区分大小写）；`iconv -f UTF-16 -t ASCII`、`dos2unix` 转码 |
| `Format error for Hubbard.U.values`（或 `Hund.J.values`） | 项数与 PAO 定义不符；缺少结束标记；或最后一行数据与结束标记之间夹了空行/注释行（4.0 源码 `Inputtools.c` 读完数据后只再读一行来找结束标记） | 按第 4 节逐项核对；块内不要写注释行，结束标记紧跟最后一行数据 |
| 加了 U 但带隙/磁矩几乎不变 | `scf.Hubbard.U` 未打开；`scf.SpinPolarization off`；U 写到了 `2d` 或因项数错位落到别的轨道 | 检查 `.std` 中 `DFT+U Type and DC`、`.out` 中 `Uhub` ≠ 0、占据矩阵是否整数化 |
| Type 2 结果与预期不符 | 忘记写 `<Hund.J.values`，或 `scf.Yukawa on` 时 J 块被忽略 | Type 2 按手册要求给出 J；检查 `.std` 打印的 Type 与 DC |
| SCF 振荡、`NormRD` 长期不降 | 仍在用密度混合；混合权重太大 | 改 `RMM-DIISH`，降低 `scf.Init.Mixing.Weight`，增大 `scf.Mixing.StartPulay`；运行中可写 `System.Name_SCF_keywords` 文件即时调整（第 5 章） |
| AFM 计算收敛到铁磁或非磁 | 初始占据设置错误；初猜磁矩太小 | 检查坐标行上下自旋；增大初猜磁矩；先在 GGA 下收敛 AFM 再加 U 重启 |
| 同一 U 下多次计算得到不同带隙 | 收敛到了不同的轨道/磁局部极小 | 比较 `Utot` 取最低；对简并体系试轨道极化开关 |
| `'full' local propjector in DFT+U is not compatible with stress calculation.`（源码原文拼写） | `scf.Hubbard.Occupation full` 与应力计算（如晶胞优化）同时使用 | 改用 `dual` |
| 非共线约束或 Zeeman 项报 `... is supported in case of scf.Hubbard.Occupation=dual.` | 非共线约束/Zeeman 只支持 `dual` | 改用 `dual` |
| 与 VASP 的 GGA+U 结果对不上 | U 投影定义不同、DC 不同、磁矩计算方式不同 | 不追求逐位一致；比较趋势，并在同一代码内做 U 扫描 |

---

## 练习

1. 用 `scan_U.sh` 对 NiO 做 U = 0–8 eV 扫描，画出带隙与 Ni 磁矩随 U 的变化；确认 U = 0（`scf.Hubbard.U on`）与 `NiO_AFM_GGA.dat` 的 `Utot` 一致。
2. 复制 `NiO_AFM_UJ_cFLL.dat`，改为 `scf.dc.Type sFLL`、U = 5、J = 0，与 `NiO_AFM_U_Dudarev.dat` 比较 `Utot`、带隙和磁矩，验证二者在数值精度内一致。
3. 固定 U − J = 5 eV，取 (U, J) = (5, 0)、(5.5, 0.5)、(6, 1) 并分别用 sFLL 和 cFLL 计算，观察 J 对带隙和磁矩的影响如何依赖 DC 类型。
4. 对 NiO 再构造一个铁磁输入（两个 Ni 初始占据同号），在 GGA 与 GGA+U（5 eV）下计算 $E_{\mathrm{FM}} - E_{\mathrm{AFM}}$，体会 U 如何减弱超交换。
5. 用 `MnO_AFM_U.dat` 对 U = 2、4、6 eV 计算，比较 Mn 磁矩是否基本不变而带隙明显变化，并解释原因（d⁵ 半满壳层）。

---

## 参考资料

官方文档：

- OpenMX 4.0 手册第 40 章 DFT+U：<https://www.openmx-square.org/openmx_man4.0/s40_dft-plus-u.html>
- 40.1 标准设置（Type、DC、占据算符、轨道极化）：<https://www.openmx-square.org/openmx_man4.0/s40_1_standard.html>
- 40.2 附加功能（`scf.Slater.Ratio`、`scf.Yukawa`）：<https://www.openmx-square.org/openmx_man4.0/s40_2_extras.html>
- 第 16 章 SCF 收敛（RMM-DIISH 的适用场景）：<https://www.openmx-square.org/openmx_man4.0/s16_scf.html>
- 关键词列表：<https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html>
- 推荐基组表（DFT_DATA19）：<https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html>
- 官方示例：`openmx4.0/work/NiO-cFLL.dat`、`NiO-sFLL.dat`、`NiO-Yukawa.dat`、`CoO_NC.dat`、`MnO_NC.dat`

相关 GitHub 仓库（只作参考链接，详见 [../../Resources/github-repos.md](../../Resources/github-repos.md)）：

- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，含 NiO DFT+U 示例与 Python 后处理笔记本（无许可证，仅供阅读）
- [by-student-2017/cif2cell-informal](https://github.com/by-student-2017/cif2cell-informal)：cif2cell 非官方版，可输出 OpenMX 输入并带 LDA+U 选项
- [pohao82/vasp-omx-moment-setter](https://github.com/pohao82/vasp-omx-moment-setter)：交互式设置 FM/AFM/非共线初始磁矩并写入 OpenMX 坐标块
- [mailhexu/TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX)：从 OpenMX 的 scfout 计算交换常数，可用于比较不同 U 下的超交换

关键文献：

- V. I. Anisimov, J. Zaanen, O. K. Andersen, Phys. Rev. B 44, 943 (1991). DOI: [10.1103/PhysRevB.44.943](https://doi.org/10.1103/PhysRevB.44.943)
- A. I. Liechtenstein, V. I. Anisimov, J. Zaanen, Phys. Rev. B 52, R5467 (1995). DOI: [10.1103/PhysRevB.52.R5467](https://doi.org/10.1103/PhysRevB.52.R5467)
- S. L. Dudarev et al., Phys. Rev. B 57, 1505 (1998). DOI: [10.1103/PhysRevB.57.1505](https://doi.org/10.1103/PhysRevB.57.1505)
- M. T. Czyżyk, G. A. Sawatzky, Phys. Rev. B 49, 14211 (1994)（FLL 与 AMF 双计数）. DOI: [10.1103/PhysRevB.49.14211](https://doi.org/10.1103/PhysRevB.49.14211)
- M. J. Han, T. Ozaki, J. Yu, Phys. Rev. B 73, 045110 (2006)（OpenMX 中基于 PAO 的 LDA+U 与占据算符）. DOI: [10.1103/PhysRevB.73.045110](https://doi.org/10.1103/PhysRevB.73.045110)
- S. Ryee, M. J. Han, Sci. Rep. 8, 9559 (2018)（双计数、自旋密度与 Hund J 在不同 DFT+U 泛函中的作用）. DOI: [10.1038/s41598-018-27731-4](https://doi.org/10.1038/s41598-018-27731-4)
- M. Cococcioni, S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)（线性响应 U）. DOI: [10.1103/PhysRevB.71.035105](https://doi.org/10.1103/PhysRevB.71.035105)
- L. Wang, T. Maxisch, G. Ceder, Phys. Rev. B 73, 195107 (2006)（拟合氧化焓确定 U）. DOI: [10.1103/PhysRevB.73.195107](https://doi.org/10.1103/PhysRevB.73.195107)
- G. A. Sawatzky, J. W. Allen, Phys. Rev. Lett. 53, 2339 (1984)（NiO 带隙）. DOI: [10.1103/PhysRevLett.53.2339](https://doi.org/10.1103/PhysRevLett.53.2339)
- A. K. Cheetham, D. A. O. Hope, Phys. Rev. B 27, 6964 (1983)（MnO、NiO 等的磁矩）. DOI: [10.1103/PhysRevB.27.6964](https://doi.org/10.1103/PhysRevB.27.6964)
- J. van Elp et al., Phys. Rev. B 44, 1530 (1991)（MnO 电子结构与带隙）. DOI: [10.1103/PhysRevB.44.1530](https://doi.org/10.1103/PhysRevB.44.1530)

相关章节：[第4章 基组与赝势](04-pao-vps.md)、[第5章 SCF 与收敛](05-scf-kpoints-convergence.md)、[第9章 DOS/PDOS](09-dos-pdos-cohp-optics.md)、[第11章 自旋与非共线](11-spin-noncollinear.md)、[第18章 并行与报错](18-parallel-on-errors.md)、[通用后处理与可视化](../../Common/tutorials/05-postprocess-visualization.md)。

[← 上一章](12-soc-topology.md) | [返回目录](../README.md) | [下一章 →](14-wannier.md)
