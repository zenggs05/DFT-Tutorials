# 第11章 自旋极化、非共线磁性与交换耦合（Spin Polarization, Noncollinear Magnetism & Exchange Coupling with OpenMX）

本章讲解 OpenMX 中磁性体系的完整计算流程：共线自旋极化（`scf.SpinPolarization on`）与初始磁矩的设置，铁磁/反铁磁构型与 Mulliken/Voronoi 磁矩，非共线磁性（`scf.SpinPolarization nc`）的 Euler 角输入与矢量可视化，约束 DFT（`scf.Constraint.NC.Spin`）与 Zeeman 项，磁各向异性能 MAE 的二次变分法（`scf.restart c2n`），基于广义 Bloch 定理的自旋螺旋，以及用 `jx`、TB2J-OpenMX 计算海森堡交换耦合 J_ij 并估算居里温度（OpenMX magnetism tutorial: spin polarization, noncollinear magnetism, constrained DFT, MAE, spin spiral, exchange coupling Jij）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据 VPS 的价电子数，正确写出坐标行第 6、7 列的上/下自旋初始电子数，从而设定铁磁（FM）、反铁磁（AFM）或非磁（NM）初态；
- 从 `.out` 中读出总自旋磁矩、各原子 Mulliken / Voronoi 磁矩，并理解 OpenMX 的"自旋磁矩符号约定"与物理磁矩相反；
- 写出非共线计算坐标行的第 8–13 列（自旋 Euler 角、轨道磁矩 Euler 角、约束开关、+U 轨道极化开关），用 XCrySDen 显示 `nc.xsf` 自旋矢量；
- 用约束 DFT 计算能量随自旋夹角的变化曲线，了解 Zeeman 项的用途与限制；
- 用二次变分法（共线 SCF → `scf.restart c2n` 单次非共线对角化）计算 MAE；
- 用 `scf.Generalized.Bloch` + `Spin.Spiral.Vector` 做自旋螺旋计算，并知道 4.0 中必须改用 `scf.Generation.Kpoint Regular`；
- 用 `HS.fileout on` + `jx` 计算 J_ij，用平均场公式估算 T_C；知道何时改用 TB2J-OpenMX、DFTforge.jl。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| Fe/Co/Ni/Mn/Cr 金属、合金的磁基态 | `scf.SpinPolarization on`，比较 FM/AFM/NM | §3、§4 |
| 自旋夹角、复杂磁结构 | `scf.SpinPolarization nc` + Euler 角 | §5 |
| 固定自旋方向求能量曲线 | 约束 DFT `scf.Constraint.NC.Spin on` | §6 |
| 外加磁场下的磁化 | Zeeman 项 `scf.NC.Zeeman.Spin` | §7 |
| 轨道磁矩 | nc + SOC，读 "Orbital moments" | §8 |
| 易磁化轴、MAE | 二次变分 `scf.restart c2n` 或全自洽约束 | §9 |
| 非公度螺旋磁序 | 广义 Bloch `scf.Generalized.Bloch on` | §10 |
| 交换常数、居里温度 | `jx`（共线）/ TB2J-OpenMX（含 SOC、DMI） | §11、§12 |
| 过渡金属氧化物（NiO、MnO） | 自旋极化 + DFT+U | [第13章](13-dft-u.md) |

> [!NOTE]
> 本章需要的基础：输入文件结构（[第3章](03-input-first-calc.md)）、PAO/VPS 与价电子数（[第4章](04-pao-vps.md)）、SCF 混合与 k 点收敛（[第5章](05-scf-kpoints-convergence.md)）。磁性体系对 k 点和电子温度比非磁体系敏感得多，务必先做收敛测试。VASP 中同类内容可对照 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)。

---

## 2. 原理速览

### 2.1 共线与非共线自旋密度

共线自旋极化把电子分成上、下两种自旋，分别求解 Kohn–Sham 方程，磁化密度 $m(\mathbf r)=\rho_\uparrow(\mathbf r)-\rho_\downarrow(\mathbf r)$ 只有"大小与正负"，没有空间方向。非共线 DFT 中每个 Kohn–Sham 轨道是二分量旋量，密度是 2×2 矩阵：

$$
\rho_{\alpha\beta}(\mathbf r)=\sum_{n} f_n\,\psi_{n\alpha}(\mathbf r)\psi^{*}_{n\beta}(\mathbf r),\qquad
\rho=\mathrm{Tr}\,\rho_{\alpha\beta},\qquad
\mathbf m=\mathrm{Tr}\,(\boldsymbol\sigma\,\rho_{\alpha\beta}),
$$

每个实空间格点上的磁化方向都可以不同。OpenMX 实现的是这种最一般的形式（不限制每个格点的自旋量子化轴），并可自洽地加入自旋轨道耦合（SOC，[第12章](12-soc-topology.md)）。

> [!IMPORTANT]
> **没有 SOC 时，自旋方向相对晶格的取向是简并的**：把所有自旋整体转一个角度，能量不变；只有自旋之间的**相对**夹角（FM、AFM、螺旋）有能量意义。要谈"易磁化轴"、MAE、轨道磁矩，必须打开 SOC。

### 2.2 OpenMX 的磁矩符号约定

OpenMX 手册明确说明：输出的自旋磁矩与轨道磁矩都是"OpenMX convention"，**与电子的物理磁矩符号相反**。原因是电子的磁矩 $\boldsymbol\mu_s=-g\mu_B\mathbf s/\hbar$ 与自旋反平行，而 OpenMX 报告的是"自旋"意义下的量（上自旋多则为正）。这一约定对共线计算的上/下自旋同样成立。实际影响：

- 报告"某原子磁矩 2.3 μB"时只关心大小或相对符号（FM/AFM），不受影响；
- 非共线计算中 Euler 角 (θ, φ) 描述的是**自旋（OpenMX 约定）**方向，物理磁矩方向与之相反；与中子衍射给出的磁矩方向比较时要翻转；
- MAE 只依赖 $\sin^2\theta$ 之类的偶函数，不受影响；DMI 手性、Zeeman 场方向等与符号有关的量，要按约定仔细换算，最好先用小体系测试一次。

### 2.3 海森堡模型与 Liechtenstein 公式

把磁性体系的低能激发映射到经典海森堡模型

$$
H=-\sum_{i\neq j}J_{ij}\,\mathbf e_i\cdot\mathbf e_j ,
$$

$J_{ij}>0$ 表示铁磁耦合。Liechtenstein 等人用磁力定理把 $J_{ij}$ 写成单粒子格林函数与在位交换劈裂 $\hat P_i=\hat H_{i\uparrow}-\hat H_{i\downarrow}$ 的乘积，只需一次共线自洽计算即可得到任意距离的 $J_{ij}$。OpenMX 的后处理程序 `jx` 在局域基组上实现了这一公式，并用 Ozaki 的有限极点近似处理费米函数（约 60 个极点即可达到 0.05 meV 精度，见 Terasawa et al. 2019）。在平均场近似下（手册式 (48)），单一磁性子格的居里温度为

$$
k_B T_C^{\rm MFA}=\frac{2}{3}\sum_{\mathbf R\neq 0}J_{0\mathbf R}.
$$

> [!NOTE]
> 不同程序/论文的海森堡哈密顿量前因子（是否含 1/2、是否把 $S^2$ 吸收进 J、求和是否重复计数）各不相同。比较不同来源的 J 之前，必须先统一约定；本章所有公式采用 OpenMX 手册的约定。

### 2.4 约束 DFT、二次变分与广义 Bloch 定理

- **约束 DFT**：在能量泛函中加一个罚函数，使各原子自旋方向保持在指定的 Euler 角附近（`on` 只约束方向，`on2` 同时约束大小），用来计算能量随自旋夹角的曲线。罚函数能量记为 `Ucs`。
- **二次变分法**：SOC 相对较弱时，先做无 SOC 的共线自洽，再用其密度构造含 SOC 的非共线哈密顿量只对角化一次；依据 Harris 泛函，能量差的误差是密度误差的二阶小量，所以比较不同自旋取向的**能带能量 Uele** 就能得到 MAE。
- **广义 Bloch 定理**：无 SOC 时，自旋螺旋 $\mathbf m_i=m\,(\cos(\mathbf q\cdot\mathbf R_i+\varphi)\sin\theta,\ \sin(\mathbf q\cdot\mathbf R_i+\varphi)\sin\theta,\ \cos\theta)$ 可以在一个原胞内计算，无需超胞。SOC 破坏这种螺旋对称性，所以两者不能同时使用。

---

## 3. 共线自旋极化：bcc Fe 铁磁（Fe_bcc_FM.dat）

完整文件见 [../examples/11-fe-bcc-fm/README.md](../examples/11-fe-bcc-fm/README.md)。

### 3.1 初始占据：磁性计算的第一步

坐标行第 6、7 列是该原子上、下自旋的**初始电子数**，二者之和必须等于 VPS 的价电子数，二者之差就是初始磁矩（μB）。常用磁性元素（DFT_DATA19，GGA-PBE）：

| VPS | 价电子数 | 典型 FM 初值（上 / 下） | 说明 |
|---|---|---|---|
| `Cr_PBE19` | 14 | 8.0 / 6.0 | 含 3s3p 半芯态 |
| `Mn_PBE19` | 15 | 10.0 / 5.0 | Mn²⁺ 类高自旋可给到 5 μB |
| `Fe_PBE19H` | 16 | 10.0 / 6.0 | "Hard"，含 3s3p 半芯 |
| `Fe_PBE19S` | 14 | 9.0 / 5.0 | "Soft"，只含 3p |
| `Co_PBE19H` / `Co_PBE19S` | 17 / 15 | 10.0 / 7.0、9.0 / 6.0 | |
| `Ni_PBE19H` / `Ni_PBE19S` | 18 / 16 | 9.5 / 8.5、8.5 / 7.5 | Ni 磁矩小，初值不宜过大 |
| `O_PBE19` | 6 | 3.0 / 3.0 | 配体一般不给磁矩 |

> [!TIP]
> 初始磁矩宁大勿小：给 2–4 μB 的初值通常比 0.5 μB 更容易落到正确的高自旋解；但不要超过物理上可能的最大值（例如 Fe 的 d 壳层最多 5 μB）。仓库脚本 `../scripts/omx_magconf.py` 可以根据价电子数和目标磁矩自动写出上/下占据（FM/AFM/自定义），并检查电子数守恒，用 `python ../scripts/omx_magconf.py -h` 查看用法。

### 3.2 输入文件

```text
#
# bcc Fe, ferromagnetic, collinear spin-polarized GGA
#
System.CurrrentDirectory    ./
System.Name                 fe_fm
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            1          # writes fe_fm.sden.cube

Species.Number              1
<Definition.of.Atomic.Species
  Fe   Fe5.5H-s3p2d1   Fe_PBE19H
Definition.of.Atomic.Species>

Atoms.Number                1
Atoms.SpeciesAndCoordinates.Unit   FRAC
# Fe_PBE19H: 16 valence electrons -> 10.0 up + 6.0 down = initial moment 4 muB
<Atoms.SpeciesAndCoordinates
  1  Fe   0.0  0.0  0.0   10.0  6.0
Atoms.SpeciesAndCoordinates>

# bcc primitive cell, a = 2.866 Ang (expt.)
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  -1.4330   1.4330   1.4330
   1.4330  -1.4330   1.4330
   1.4330   1.4330  -1.4330
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.EigenvalueSolver        Band
scf.Kgrid                   16 16 16
scf.energycutoff            250.0
scf.ElectronicTemperature   300.0
scf.maxIter                 200
scf.criterion               1.0e-8
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.02
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15

MD.Type                     Nomd
Voronoi.charge              on         # Voronoi spin moment for comparison
```

### 3.3 关键词说明

| 关键词 | 取值 | 说明 |
|---|---|---|
| `scf.XcType` | GGA-PBE | 自旋极化必须用 LSDA-CA、LSDA-PW 或 GGA-PBE；若写 `LDA` 又打开自旋，OpenMX 报 `SpinP_switch should be OFF for this exchange functional.` |
| `scf.SpinPolarization` | on | 共线自旋极化；`nc` 为非共线；`off` 为非磁 |
| `Fe5.5H-s3p2d1` / `Fe_PBE19H` | Standard | DFT_DATA19 推荐表中 Fe（Hard 赝势）的 Standard 组合 |
| `scf.Kgrid` | 16×16×16 | 磁性金属的磁矩与能量对 k 点敏感；需做收敛测试（[第5章](05-scf-kpoints-convergence.md)） |
| `scf.ElectronicTemperature` | 300 K | 磁性金属不要盲目调高：过高的电子温度会压低磁矩 |
| `scf.Mixing.Type` | RMM-DIISK | 金属与磁性体系最稳的方案；初始混合权重取小 |
| `Voronoi.charge` | on | 额外输出 Voronoi 自旋磁矩（多组基函数时 Mulliken 自旋矩可能偏大） |

### 3.4 运行与读结果

```bash
mpirun -np 8 openmx Fe_bcc_FM.dat -nt 1 > fe_fm.std 2>&1
grep "Utot\." fe_fm.out
grep "Total spin moment" fe_fm.out
grep -A 6 "Mulliken populations" fe_fm.out
grep "by Voronoi charges" fe_fm.out
```

"Mulliken populations" 表（4.0 的格式）：

```text
  Total spin moment (muB)   2.4xxxxxxxx

                    Up spin      Down spin     Sum           Diff
      1   Fe       9.2xxxxxxx   6.7xxxxxxx   16.000000000   2.4xxxxxxxx
```

（示意输出，数值仅供参考）

- `Diff` 列 = 上自旋布居 − 下自旋布居 = 该原子的 Mulliken 自旋磁矩（μB，OpenMX 约定）；
- 实验上 bcc Fe 的饱和磁矩约 2.2 μB/atom，GGA 一般给出 2.2–2.3 μB 左右；Mulliken 值（多组 d 函数时）可能略偏大，Voronoi 值可作对照；
- 若磁矩塌缩到 0，说明初值不对或混合太激进（见 §13）。

**非磁对照**：把 `scf.SpinPolarization` 改为 `off`，初始占据写 8.0 8.0，System.Name 改为 `fe_nm`，重算后比较 `Utot`：E(FM) 应明显低于 E(NM)（每原子数百 meV 量级）。两次计算必须使用相同的基组、截断能、k 点与电子温度。

### 3.5 自旋密度可视化

`level.of.fileout 1` 时自旋极化计算会写出 `fe_fm.sden.cube`（ρ↑ − ρ↓，e/Bohr³）。在 VESTA 中取 ±0.005 左右的等值面，可以看到以 Fe 为中心、具有 t2g/eg 角向特征的自旋密度，以及间隙区很小的负自旋极化。

---

## 4. 反铁磁：bcc Cr 的 FM / AFM / NM 对比（Cr_bcc_AFM.dat）

完整文件见 [../examples/11-cr-bcc-afm/README.md](../examples/11-cr-bcc-afm/README.md)。

### 4.1 两种写法

**写法 A：同一 species，相反的初始占据**（最简单，本例采用）：

```text
<Atoms.SpeciesAndCoordinates
  1  Cr   0.0  0.0  0.0    8.5  5.5
  2  Cr   0.5  0.5  0.5    5.5  8.5
Atoms.SpeciesAndCoordinates>
```

**写法 B：同一元素定义两个 species**，便于按 species 统计 PDOS、给不同子格加不同的 U 或不同基组：

```text
Species.Number              2
<Definition.of.Atomic.Species
  Cr1   Cr6.0-s3p2d1   Cr_PBE19
  Cr2   Cr6.0-s3p2d1   Cr_PBE19
Definition.of.Atomic.Species>
<Atoms.SpeciesAndCoordinates
  1  Cr1   0.0  0.0  0.0    8.5  5.5
  2  Cr2   0.5  0.5  0.5    5.5  8.5
Atoms.SpeciesAndCoordinates>
```

两种写法在物理上等价——OpenMX 的 SCF 不靠"species 相同"来强加对称性，决定磁构型的是初始占据。

### 4.2 输入文件（AFM）

```text
#
# bcc Cr, simple commensurate antiferromagnet (CsCl-type) in the 2-atom cubic cell
#
System.CurrrentDirectory    ./
System.Name                 cr_afm
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            1

Species.Number              1
<Definition.of.Atomic.Species
  Cr   Cr6.0-s3p2d1   Cr_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
# Cr_PBE19: 14 valence electrons; atom 1 starts with +3 muB, atom 2 with -3 muB
<Atoms.SpeciesAndCoordinates
  1  Cr   0.0  0.0  0.0    8.5  5.5
  2  Cr   0.5  0.5  0.5    5.5  8.5
Atoms.SpeciesAndCoordinates>

# simple cubic 2-atom cell, a ~ 2.88 Ang (close to expt.)
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  2.8800   0.0000   0.0000
  0.0000   2.8800   0.0000
  0.0000   0.0000   2.8800
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.EigenvalueSolver        Band
scf.Kgrid                   14 14 14
scf.energycutoff            250.0
scf.ElectronicTemperature   300.0
scf.maxIter                 300
scf.criterion               1.0e-8
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.02
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15

MD.Type                     Nomd
Voronoi.charge              on
```

### 4.3 运行三种磁态

示例目录中的 `run_fm_afm_nm.sh` 用 sed 生成 FM（两原子都写 8.5 5.5）与 NM（`scf.SpinPolarization off`，7.0 7.0）两个输入，依次运行并汇总：

```bash
bash run_fm_afm_nm.sh           # runs cr_afm, cr_fm, cr_nm
grep "Utot\." cr_afm.out cr_fm.out cr_nm.out
grep -A 5 "Mulliken populations" cr_afm.out
```

判断标准：

- AFM 计算中两个 Cr 的 `Diff` 大小相近、符号相反，`Total spin moment` ≈ 0；
- 能量顺序预期 E(AFM) < E(NM)，而 FM 往往不稳定（可能塌缩成非磁或很小磁矩）；
- 实际 Cr 的基态是非公度自旋密度波（SDW），简单 AFM 只是近似；GGA 通常高估 Cr 的局域磁矩，这是已知的泛函问题，不是输入错误。

> [!WARNING]
> AFM 计算"跑成"FM 或 NM 是最常见的问题。先检查初始占据是否真的相反；再减小 `scf.Init.Mixing.Weight`、`scf.Max.Mixing.Weight`；或先用较高电子温度收敛，再用 `scf.restart on` 降温续算。

---

## 5. 非共线磁性（Mol_MnO_NC.dat）

完整文件见 [../examples/11-mno-noncollinear/README.md](../examples/11-mno-noncollinear/README.md)。

### 5.1 坐标行的 13 列

`scf.SpinPolarization nc` 时，坐标行可写到 13 列：

| 列 | 含义 |
|---|---|
| 1 | 原子序号 |
| 2 | species 名 |
| 3–5 | x, y, z（或分数坐标） |
| 6, 7 | 上、下自旋初始电子数（二者之差 = 初始磁矩大小） |
| 8, 9 | 自旋磁矩（OpenMX 约定）方向的 Euler 角 θ, φ（度）；同时用于构造初始自旋密度 |
| 10, 11 | 轨道磁矩方向的 Euler 角 θ, φ（度） |
| 12 | 约束开关：当使用 `scf.Constraint.NC.Spin`、`scf.NC.Zeeman.Spin` 或 `scf.NC.Zeeman.Orbital` 时，1 = 对该原子施加约束/场，0 = 不施加 |
| 13 | DFT+U 轨道极化开关：`on` 促进轨道极化，`off` 为常规计算（[第13章](13-dft-u.md)） |

Euler 角按常规球坐标理解：θ 是与 z 轴的夹角，φ 是在 xy 平面内从 x 轴起算的方位角，方向矢量为 (sinθ cosφ, sinθ sinφ, cosθ)。非磁原子或不需要约束时，第 8 列以后可以省略（官方 Bi2Se3 的非共线输入就只写 7 列）。

### 5.2 输入文件

MnO 分子，键长取 1.65 Å（近似值，未优化），基态为高自旋（总自旋磁矩约 5 μB）。初始占据 Mn 10.0/5.0（Mn_PBE19 价电子 15），O 3.0/3.0；自旋方向取 θ = 45°、φ = 0°（在 xz 平面内与 z 轴成 45°），并打开 SOC 与约束：

```text
#
# MnO molecule: noncollinear DFT with spin-orbit coupling and constrained spin direction
#
System.CurrrentDirectory    ./
System.Name                 mno_nc
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            1          # writes mno_nc.nc.xsf, mno_nc.ncsden.xsf, mno_nc.nco.xsf

Species.Number              2
<Definition.of.Atomic.Species
  Mn   Mn6.0-s3p2d1   Mn_PBE19
  O    O6.0-s2p2d1    O_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   Ang
# columns: no sp x y z up dn theta_s phi_s theta_o phi_o constraint +U-orbital-polarization
#          (never put comment lines INSIDE a block: OpenMX reads block lines one by one)
<Atoms.SpeciesAndCoordinates
  1  Mn   0.00000  0.00000  0.00000  10.0  5.0  45.0  0.0   45.0  0.0   1    off
  2  O    1.65000  0.00000  0.00000   3.0  3.0  45.0  0.0   45.0  0.0   1    off
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  14.0   0.0   0.0
   0.0  14.0   0.0
   0.0   0.0  14.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on
scf.Constraint.NC.Spin      on         # constrain the spin direction (not the magnitude)
scf.Constraint.NC.Spin.v    0.5        # prefactor of the penalty functional (eV)
scf.EigenvalueSolver        Cluster
scf.Kgrid                   1 1 1
scf.energycutoff            220.0
scf.ElectronicTemperature   300.0
scf.maxIter                 300
scf.criterion               1.0e-9
scf.Mixing.Type             RMM-DIIS
scf.Init.Mixing.Weight      0.005
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       10

MD.Type                     Nomd
scf.restart                 off
```

### 5.3 读非共线结果

```bash
mpirun -np 4 openmx Mol_MnO_NC.dat -nt 1 > mno_nc.std 2>&1
grep -A 8 "Mulliken populations" mno_nc.out
grep -A 8 "Orbital moments" mno_nc.out
grep "Ucs\." mno_nc.out
```

非共线的 Mulliken 表多了 theta、phi 两列：

```text
   Total spin moment (muB)   4.99xxxxxxx   Angles (Deg) 45.0xxxxxx  -0.0xxxxxx

               Up       Down      Sum      Diff        theta      phi
    1   Mn   9.6xxxx  4.7xxxx  14.3xxxx   4.8xxxx    45.0xxxx   -0.0xxxx
    2    O   3.4xxxx  3.2xxxx   6.6xxxx   0.1xxxx    45.0xxxx   -0.0xxxx
```

（示意输出，数值仅供参考）

这里的 Up/Down 是沿该原子局域自旋量子化轴的布居，`Diff` 是局域自旋磁矩大小，theta/phi 是其方向（OpenMX 约定）。检查要点：

- 总自旋磁矩约 5 μB（高自旋 MnO）；
- 有约束时，各原子 theta/phi 应接近输入的 45°/0°；
- 能量分解中的 `Ucs`（自旋取向约束能）应远小于你要比较的能量差，否则约束本身在"污染"能量。

### 5.4 矢量可视化（XCrySDen）

非共线且 `level.of.fileout` ≥ 1 时生成：

| 文件 | 内容 |
|---|---|
| `mno_nc.nc.xsf` | 每个原子的 Mulliken 自旋磁矩矢量（OpenMX 约定） |
| `mno_nc.ncsden.xsf` | 实空间格点上的自旋磁化矢量场 |
| `mno_nc.nco.xsf` | 每个原子的轨道磁矩矢量 |

在 XCrySDen 中 File → Open Structure → Open XSF，然后 Display → Forces 打开箭头显示，Modify → Force Settings 调整箭头长度。OpenMX Viewer 也能读取 XSF 中的原子矢量。

> [!TIP]
> 非共线 SCF 比共线难收敛，因为改变自旋方向的能量代价很小。可行的策略：(1) 先做共线计算，确认磁矩大小合理；(2) 非共线计算先用约束 DFT 把方向"钉住"收敛，再关掉约束、`scf.restart on` 续算；(3) 混合权重取小、`scf.Mixing.History` 取大（20–40）。

---

## 6. 约束 DFT：能量随自旋夹角的曲线

### 6.1 关键词

```text
scf.Constraint.NC.Spin       on       # off | on | on2
scf.Constraint.NC.Spin.v     0.5      # eV, strength of the penalty functional
```

| 取值 | 约束内容 | 适用 |
|---|---|---|
| `on` | 只约束自旋**方向**到坐标行第 8、9 列给定的 Euler 角，大小自由优化 | 绝大多数夹角扫描、MAE 全自洽对照 |
| `on2` | 方向与大小都约束到初始值（第 6、7 列之差） | 需要固定磁矩大小的特殊研究 |

- 坐标行第 12 列为 1 的原子才受约束；
- 约束 DFT 与 DFT+U 可以同时使用；
- 只能用于非共线计算，共线时报 `The constraint scheme is not supported for a collinear DFT calculation.`

### 6.2 推荐流程（手册第 41 章的做法）

局域 d/f 电子的部分占据态在占据数空间中有很多局域解，直接从大角度起步容易落入错误的解。手册推荐"从易轴出发、小步旋转、逐步续算"：

1. 先在易轴（或 FM 共线方向）上做一次自洽，得到重启文件；
2. 下一个角度只转约 15°–30°，`scf.restart on` 读取上一步的重启文件（`System.Name` 必须相同）；
3. 每步都把 `.out` 另存，读 `Utot` 与 `Ucs`；
4. 对每个角度重复，得到平滑曲线。

下面是一个通用的扫描骨架（以两原子二聚体的相对夹角为例：原子 1 转 −θ/2，原子 2 转 +θ/2）：

```bash
#!/bin/bash
# relative-angle scan with constrained DFT; template.dat contains TH1/TH2 placeholders
# in columns 8 and 10 of the two coordinate lines and a line "scf.restart  RESTART"
first=1
for th in 0 30 60 90 120 150 180; do
  h=$(echo "$th/2" | bc -l)
  rs=on; [ $first -eq 1 ] && rs=off
  sed -e "s/TH1/-$h/g" -e "s/TH2/$h/g" -e "s/RESTART/$rs/" template.dat > run.dat
  mpirun -np 8 openmx run.dat -nt 1 > run.std 2>&1
  cp dimer.out dimer_theta${th}.out
  printf "%6s %s %s\n" "$th" "$(grep 'Utot\.' dimer.out | awk '{print $2}')" \
                             "$(grep 'Ucs\.'  dimer.out | awk '{print $2}')"
  first=0
done
```

官方 work 目录中的 `Cr2_CNC.dat`（Cr 二聚体，键长 2.0 Å）就是这一类计算的示例，手册图 40 给出了其能量与磁矩随夹角的变化。本仓库的 MnO 示例目录提供了对整个分子做"整体转动"的扫描脚本 `scan_theta.sh`（有 SOC 时得到分子的磁各向异性曲线，能量差很小，需要很紧的收敛标准）。

> [!WARNING]
> 罚函数参数 `scf.Constraint.NC.Spin.v` 太小，方向约束不住；太大，SCF 难收敛且 `Ucs` 偏大。建议用两个不同的 v（如 0.5 和 1.0 eV）各算一遍关键角度，确认能量差不随 v 明显变化。

---

## 7. Zeeman 项：外加磁场

Zeeman 项把均匀磁场与磁矩的相互作用直接加到哈密顿量中。与约束 DFT 不同，它不仅能改变方向，也能改变磁矩大小，且总能量中包含 Zeeman 能本身（能量分解中的 `Uzs`、`Uzo`）。

```text
scf.NC.Zeeman.Spin          on        # Zeeman term for the spin moment
scf.NC.Mag.Field.Spin       100.0     # Tesla
scf.Constraint.NC.Spin      off       # must be off together with the spin Zeeman term

scf.NC.Zeeman.Orbital       on        # Zeeman term for the orbital moment (needs SOC)
scf.NC.Mag.Field.Orbital    100.0     # Tesla
```

- 场的方向逐原子给定：作用于自旋的场方向由第 8、9 列的 Euler 角给出，作用于轨道磁矩的场方向由第 10、11 列给出，第 12 列为 1 的原子才受场；
- 自旋 Zeeman 与自旋约束互斥：同时打开会报 `For spin magnetic moment, the constraint scheme and the Zeeman term are mutually exclusive.`；
- 轨道 Zeeman 需要 SOC，否则报 `The Zeeman term for orbital moment is not supported without the SO term.`；
- Zeeman 项要求 `scf.Hubbard.Occupation dual`（默认值即为 dual）；
- 只能用于非共线计算。

> [!NOTE]
> 100 T 在实验上已是极强的场，但对应的 Zeeman 能仅约 $\mu_B B\approx 5.8$ meV/μB。手册中把它用作"按原子指定方向的约束手段"来研究各种非共线结构。解释结果时注意总能量包含 Zeeman 能，不能直接与零场能量曲线比较。

---

## 8. 轨道磁矩

非共线计算默认就会计算各原子的轨道磁矩，但只有打开 SOC 时它才可能非零。`.out` 中的 "Orbital moments" 段给出每个原子的轨道磁矩大小与方向（OpenMX 约定），以及按 s、p、d、f 和多重度分解的值（"Decomposed Orbital Moments"）：

```bash
grep -A 10 "Orbital moments" xxx.out
grep -A 60 "Decomposed Orbital Moments" xxx.out | less
```

要点：

- 3d 过渡金属中轨道磁矩被晶体场部分猝灭，只有 SOC（以及 DFT+U 的轨道极化）才能恢复有限值；
- 轨道磁矩对基组（需要包含极化函数的较大基组）、U 值、k 点与电子温度都很敏感，自旋磁矩则相对稳定；
- 手册表 7 用 LDA+U（Mn 3 eV、Fe 5 eV、Co 7 eV、Ni 7 eV）给出了 MnO、FeO、CoO、NiO 的自旋与轨道磁矩，例如 FeO 的轨道磁矩约 0.80 μB，可作为复现练习的参考（输入文件 `FeO_NC.dat` 等在官方 work 目录）。

---

## 9. 磁各向异性能 MAE：二次变分法

### 9.1 两种方法

| 方法 | 做法 | 比较的能量 | 优缺点 |
|---|---|---|---|
| 全自洽（约束 DFT） | nc + SOC + `scf.Constraint.NC.Spin on`，每个角度完全自洽 | `Utot` | 最严格；昂贵，局域 d/f 体系可能陷入不同局域解 |
| 二次变分（推荐） | 共线自洽 → `scf.restart c2n` 读入 → nc + SOC 只对角化一次 | `Uele` | 快、稳定；只适用于铁磁体 |

手册以 L1₀ FePt 为例：全自洽 MAE 为 2.78 meV/f.u.，二次变分在相同 k 点/电子温度下为 2.89 meV/f.u.；若第一步用 800 K、9³ k 点，第二步用 300 K、17³ k 点，得到 2.90 meV/f.u.——说明第一步可以用较粗的参数，第二步必须用低电子温度、密 k 网格。

### 9.2 第一步：共线自洽（L1₀ FePt 原胞）

以下用 2 原子的简单四方原胞（a = a_fcc/√2 ≈ 2.72 Å，c ≈ 3.71 Å，近似实验值）：

```text
System.CurrrentDirectory    ./
System.Name                 fept
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            0

Species.Number              2
<Definition.of.Atomic.Species
  Fe   Fe5.5H-s3p2d1     Fe_PBE19H
  Pt   Pt7.0-s3p2d2f1    Pt_PBE19
Definition.of.Atomic.Species>

Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Fe   0.0  0.0  0.0   10.0  6.0
  2  Pt   0.5  0.5  0.5    8.0  8.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  2.7230   0.0000   0.0000
  0.0000   2.7230   0.0000
  0.0000   0.0000   3.7100
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.EigenvalueSolver        Band
scf.Kgrid                   12 12 9
scf.energycutoff            250.0
scf.ElectronicTemperature   600.0
scf.maxIter                 200
scf.criterion               1.0e-9
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.02
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15
MD.Type                     Nomd
```

### 9.3 第二步：读入共线密度，含 SOC 对角化一次

对每个自旋方向各准备一个输入（只列出与第一步不同的行；晶胞、基组、截断能必须与第一步**完全相同**）：

```text
System.Name                 fept_nc_z              # different name for each direction
scf.SpinPolarization        nc
scf.SpinOrbit.Coupling      on
scf.restart                 c2n                    # convert the collinear restart files
scf.restart.filename        fept                   # System.Name of the first step
scf.Restart.Spin.Angle.Theta  0.0                  # 0 = along c; 90 = in-plane
scf.Restart.Spin.Angle.Phi    0.0
scf.maxIter                 1                      # one diagonalization only
scf.Init.Mixing.Weight      0.0
scf.Min.Mixing.Weight       0.0
scf.Max.Mixing.Weight       0.0
scf.Kgrid                   24 24 18               # denser than in step 1
scf.ElectronicTemperature   300.0                  # lower than in step 1
```

坐标行保持 7 列即可：`c2n` 会把所有格点上的自旋方向统一设为 `scf.Restart.Spin.Angle.Theta/Phi` 给定的方向（这也是二次变分法只适用于铁磁体的原因）。

把第一步的输入保存为 `FePt_col.dat`，第二步两个方向分别保存为 `FePt_nc_z.dat`（Theta = 0）与 `FePt_nc_x.dat`（Theta = 90，System.Name 改为 `fept_nc_x`），三个文件放在同一目录、依次运行：

```bash
mpirun -np 16 openmx FePt_col.dat  -nt 1 > fept.std 2>&1
mpirun -np 16 openmx FePt_nc_z.dat -nt 1 > fept_nc_z.std 2>&1
mpirun -np 16 openmx FePt_nc_x.dat -nt 1 > fept_nc_x.std 2>&1
grep "Uele\." fept_nc_z.out fept_nc_x.out
```

$$
\mathrm{MAE}=U_{\rm ele}(\theta=90^\circ)-U_{\rm ele}(\theta=0^\circ),
$$

MAE > 0 表示易轴沿 c（θ = 0）。本例原胞含 1 个 FePt 化学式单元，直接乘 27211.4 换算为 meV/f.u.。

> [!IMPORTANT]
> - 二次变分法比较的是 **Uele**（能带能量），不是 Utot。输出中的双计数项看起来会随角度变化，这是因为 OpenMX 用输出密度重新计算各项，不要被误导；
> - MAE 通常只有 μeV–meV 量级，第二步的 k 点必须加密到 MAE 收敛（逐步加密直到变化 < 10%），电子温度要低；
> - `scf.restart c2n` 只能用于非共线计算（否则报 `scf.restart=c2n is supported only for non-collinear calculations.`）；3.9 手册的 MAE 一节已使用同样的写法，4.0 将 `c2n` 与 `pp` 一起正式列为 `scf.restart` 的选项。

官方 work 目录的 `FePt.dat`、`FePt-NC.dat`、`FePt-NC-SCF.dat` 分别对应手册中的第一步、第二步与全自洽对照，可用来核对你的设置。

---

## 10. 自旋螺旋（广义 Bloch 定理）

### 10.1 关键词

```text
scf.SpinPolarization        nc
scf.Generalized.Bloch       on         # default off
scf.SpinOrbit.Coupling      off        # SOC is NOT allowed with the generalized Bloch theorem
Spin.Spiral.Vector          0.30 0.30 0.00   # q in units of the reciprocal lattice vectors
scf.Generation.Kpoint       Regular    # OpenMX 4.0: required for spin spirals (see below)
```

- 锥角 θ 由坐标行第 8 列给出：θ = 90° 为平面螺旋，0 < θ < 90° 为锥形螺旋；第 9 列给初始相位 φ；可配合约束 DFT 固定锥角；
- `Spin.Spiral.Vector` 以 `Atoms.UnitVectors` 的倒格矢为基（分数坐标）；
- 手册推荐自旋螺旋计算用 `scf.Mixing.Type RMM-DIISH`（基准测试中收敛最快）；
- 共线计算中打开会报 `Generalized Bloch theorem is not supported for the collinear DFT calculation.`；与 SOC 同时打开会报 `Spin-orbit coupling is not supported for the generalized Bloch theorem.`

> [!WARNING]
> **OpenMX 4.0 的已知问题**：4.0 默认 `scf.Generation.Kpoint GCenter`，与广义 Bloch 计算同时使用时，论坛报告会出现 MPI 错误 `Fatal error in internal_Comm_size: Invalid communicator`。解决办法是显式写 `scf.Generation.Kpoint Regular`（3.9 的行为，论坛帖 no=3557）。

### 10.2 示例：fcc Fe 沿 Γ–X 的螺旋

下面只列出与普通输入不同或关键的部分（`System.CurrrentDirectory`、`System.Name`、`DATA.PATH`、`level.of.*` 等同前）。fcc 原胞取 a1 = a/2(0,1,1)、a2 = a/2(1,0,1)、a3 = a/2(1,1,0)，则 `Spin.Spiral.Vector  t/2 t/2 0` 对应笛卡尔 q = (2π/a)(0, 0, t)。手册图 85 显示 fcc Fe 在 q = (0, 0, 0.6)·2π/a 与 (0.2, 0, 1)·2π/a 附近出现稳定的螺旋态（LSDA，官方输入 `Fefcc-SpinSpiral.dat`）。

```text
Species.Number              1
<Definition.of.Atomic.Species
  Fe   Fe5.5H-s3p2d1   Fe_PBE19H
Definition.of.Atomic.Species>

Atoms.Number                1
Atoms.SpeciesAndCoordinates.Unit   FRAC
# theta = 90 deg in column 8: flat spiral
<Atoms.SpeciesAndCoordinates
  1  Fe   0.0  0.0  0.0   10.0  6.0   90.0  0.0   0.0  0.0   1  off
Atoms.SpeciesAndCoordinates>

# fcc Fe primitive cell, a ~ 3.59 Ang (approximate)
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  0.000  1.795  1.795
  1.795  0.000  1.795
  1.795  1.795  0.000
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        nc
scf.Generalized.Bloch       on
scf.SpinOrbit.Coupling      off
scf.Constraint.NC.Spin      on
scf.Constraint.NC.Spin.v    1.0
Spin.Spiral.Vector          0.30 0.30 0.00     # q = (0,0,0.6) 2pi/a
scf.Generation.Kpoint       Regular
scf.EigenvalueSolver        Band
scf.Kgrid                   15 15 15
scf.energycutoff            250.0
scf.ElectronicTemperature   300.0
scf.maxIter                 800
scf.criterion               1.0e-8
scf.Mixing.Type             RMM-DIISH
scf.Init.Mixing.Weight      0.01
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.05
scf.Mixing.History          40
scf.Mixing.StartPulay       30
MD.Type                     Nomd
```

扫描 t = 0, 0.1, …, 1.0（每个 q 一个输入，`Utot` 对 t 作图）即得 E(q) 曲线；t = 0 对应 FM，t = 1 对应沿 z 的单层 AFM。fcc Fe 的螺旋基态对晶格常数非常敏感（手册引文中有专门讨论），这里的晶格常数仅作示例。

---

## 11. 交换耦合 J_ij：jx 程序（Fe_bcc_jx.dat）

完整文件见 [../examples/11-fe-bcc-jx/README.md](../examples/11-fe-bcc-jx/README.md)。

### 11.1 流程

1. 共线自旋极化 SCF，打开 `HS.fileout on`，生成 `System.Name.scfout`；
2. 写 jx 配置文件（格式与 OpenMX 输入相同）；
3. `jx System.Name.scfout jx.config`（周期体系可 MPI 并行，按 k 点分配）；
4. 把 J 按距离分壳，估算 T_C。

限制（OpenMX 4.0 手册第 45 章）：只支持**共线**计算（Cluster 与 Band 都可）；手册建议用较小的基组（如 Fe 用 s2p2d2），因为实现假设赝原子基函数足够局域，大基组会削弱这一假设。

### 11.2 SCF 输入

```text
System.CurrrentDirectory    ./
System.Name                 fe_jx
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # edit
level.of.stdout             1
level.of.fileout            0
HS.fileout                  on                          # -> fe_jx.scfout

Species.Number              1
# small basis (s2p2d2) as recommended for jx
<Definition.of.Atomic.Species
  Fe   Fe6.0H-s2p2d2   Fe_PBE19H
Definition.of.Atomic.Species>

Atoms.Number                1
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1  Fe   0.0  0.0  0.0   10.0  6.0
Atoms.SpeciesAndCoordinates>

Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  -1.4330   1.4330   1.4330
   1.4330  -1.4330   1.4330
   1.4330   1.4330  -1.4330
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.EigenvalueSolver        Band
scf.Kgrid                   16 16 16
scf.energycutoff            300.0
scf.ElectronicTemperature   300.0
scf.maxIter                 200
scf.criterion               1.0e-8
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.02
scf.Min.Mixing.Weight       0.001
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15
MD.Type                     Nomd
```

### 11.3 jx 配置文件

bcc 原胞基矢为 a1 = a/2(−1,1,1)、a2 = a/2(1,−1,1)、a3 = a/2(1,1,−1)，格矢 R = l1·a1 + l2·a2 + l3·a3。第一近邻 a/2(±1,±1,±1) 对应 (l1,l2,l3) = (1,0,0)、(0,1,0)、(0,0,1)、(1,1,1) 及其负值；第二近邻 a(±1,0,0) 等对应 (0,1,1)、(1,0,1)、(1,1,0) 及其负值；第三近邻 a(±1,±1,0) 等 12 个。配置文件（节选，完整 26 对见示例目录）：

```text
Flag.PeriodicSum     off          # off: J(i0,jR) for each listed R (Eq. 44 of the manual)
Num.Poles            60           # finite-pole approximation of the Fermi function
Num.Kgrid            16 16 16     # >= scf.Kgrid of the SCF run
Num.ij.pairs         26           # must equal the number of lines in the block
Bunch.ij.pairs       26           # reduce if memory is insufficient

<ijpairs.cellid
1 1  1  0  0
1 1  0  1  0
1 1  0  0  1
1 1  1  1  1
1 1 -1  0  0
1 1  0 -1  0
1 1  0  0 -1
1 1 -1 -1 -1
1 1  0  1  1
...
ijpairs.cellid>
```

| 关键词 | 含义与建议 |
|---|---|
| `Flag.PeriodicSum` | off：对每个列出的 (i, j, R) 计算 $J_{i\mathbf 0,j\mathbf R}$；on：计算对 R 求和的 $J_{ij}=\sum_{\mathbf R}J_{i\mathbf 0,j\mathbf R}$，此时用 `<ijpairs.nocellid` 只写 i j |
| `Num.Poles` | 费米函数有限极点数；bcc Fe 约 60 个极点可达 0.05 meV 精度，300 K 下多数体系够用，建议自己检查一次收敛 |
| `Num.Kgrid` | 布里渊区离散化，应 ≥ SCF 的 k 网格；Cluster 计算中无作用 |
| `Num.ij.pairs` | 必须等于块中的行数 |
| `Bunch.ij.pairs` | 一次同时计算的对数；内存不足时调小（会变慢） |
| `<ijpairs.cellid` | 每行 `i j l1 l2 l3`：原子 i（原胞 0）与原子 j（原胞 R）之间的 J |

### 11.4 编译与运行

```bash
cd ~/openmx4.0/source && make jx            # jx is copied to ../work
mpirun -np 8 openmx Fe_bcc_jx.dat -nt 1 > fe_jx.std 2>&1
ls -l fe_jx.scfout
mpirun -np 16 ~/openmx4.0/work/jx fe_jx.scfout Fe_bcc_jx.config | tee jx.log
```

jx 只对 k 点做 MPI 并行，进程数超过 k 点总数是浪费；Cluster 体系不并行。输出表头为 `i  j  c1  c2  c3  J [meV]  J [mRy] ...`。由于对称性，同一壳层内各 R 的 J 应相等（数值噪声范围内），这是检查 k 点与极点数是否足够的好方法。

### 11.5 平均场居里温度

示例目录的 `jx_tc_mf.py` 读取 jx 输出与 `.dat` 中的晶胞，按距离分壳列出平均 J，并按手册式 (48) 在**所列壳层截断**的近似下给出平均场 T_C：

```bash
python jx_tc_mf.py jx.log Fe_bcc_jx.dat
```

参考：手册表 9 用 Fe6.0H-s2p2d2、GGA-PBE 得到 bcc Fe 的平均场 T_C 为 1321 K（实验 1040 K），hcp Co 1640 K（实验 1131 K），fcc Ni 445 K（实验 627 K）。平均场近似通常高估 T_C；金属中 J 是长程振荡（RKKY 型）的，截断到 3 个壳层的求和只能给出粗略估计——要得到手册那样的结果，需要像官方 `Fe_Bulk_jx.config` 那样列出数百对 (i, j, R)，或使用 `Flag.PeriodicSum on`（手册式 (47)–(49)）。

> [!TIP]
> 用能量映射交叉验证：在 2 原子立方胞中比较 FM 与 AFM 的总能量，$E_{\rm AFM}-E_{\rm FM}$ 与按壳层求和的 J 应在量级和符号上一致。注意这种映射对 Fe 这类巡游磁体只是粗略检验（AFM 态的磁矩大小与 FM 不同）。

---

## 12. TB2J-OpenMX 与 DFTforge.jl

| 工具 | 输入 | 能算什么 | 备注 |
|---|---|---|---|
| `jx`（OpenMX 自带） | 共线 scfout | 各向同性 J_ij、平均场 T_C | 官方支持，MPI 并行 |
| [TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX)（[TB2J](https://github.com/mailhexu/TB2J) 插件） | scfout + `System.Name.xyz` | 各向同性 J、DMI、各向异性交换、磁振子能带，输出可供自旋动力学使用 | GPL 插件需单独安装；仓库示例包含 SOC 情形 |
| [DFTforge.jl](https://github.com/KAIST-ELST/DFTforge.jl) | scfout | J(q)、轨道分辨的交换 | Julia 实现，是 OpenMX Jx 方法的扩展 |

TB2J-OpenMX 的基本用法（据 TB2J 文档）：

```bash
pip install TB2J-OpenMX
# SCF with HS.fileout on; then in the calculation directory:
openmx2J.py --prefix fe_jx --elements Fe --kmesh 7 7 7
ls TB2J_results/
```

`--prefix` 是 OpenMX 的 `System.Name`（程序读取 `fe_jx.xyz` 与 `fe_jx.scfout`），`--kmesh` 需要按收敛测试加密（文档示例的 7×7×7 只用于演示）。其余参数用 `openmx2J.py --help` 查看。

> [!WARNING]
> TB2J-OpenMX 内部使用 OpenMX 的 scfout 读取代码，与 OpenMX 版本相关。用 4.0 生成的 scfout 时请先确认你安装的插件版本支持该格式（可用简单的 bcc Fe 与 jx 结果对照）；J 的符号与前因子约定请以 TB2J 文档 "Conventions of Heisenberg Model" 为准，与 jx 比较前先换算。

---

## 参数选择与收敛建议

| 参数 | 推荐 | 说明 |
|---|---|---|
| 初始磁矩 | 比预期稍大（Fe 4 μB、Mn 5 μB、Ni 1 μB） | 多试几种初值（FM/AFM/NM），取能量最低者 |
| `scf.Kgrid`（磁性金属） | 原胞 14³–20³ 起步 | 磁矩变化 < 0.01 μB、能量差变化 < 1 meV/atom |
| `scf.ElectronicTemperature` | 100–300 K（金属可到 500 K 帮助收敛） | 过高会压低磁矩；比较能量用同一温度 |
| `scf.Mixing.Type` | 共线/非共线 RMM-DIISK；+U、约束、螺旋 RMM-DIISH | 权重 0.01–0.05 起步，`History` 30–40 |
| `scf.criterion` | 磁构型能量差 1e-8；MAE、夹角扫描 1e-9～1e-10 | 能量差越小，标准越严 |
| 基组 | 磁矩/能量：Standard；轨道磁矩：含极化函数的大基组；jx：小基组（s2p2d2） | jx 对大基组反而不可靠 |
| `scf.Constraint.NC.Spin.v` | 0.5–1.0 eV | 检查 `Ucs` 小且结果不随 v 显著变化 |
| MAE 第二步 k 点 | 远密于第一步（如 24³ 以上） | 逐步加密至 MAE 收敛 |
| jx `Num.Poles` | 60 | 改为 100 检查 J 变化 < 0.05 meV |
| jx `Num.Kgrid` | ≥ SCF k 网格 | 同一壳层 J 的离散度可作噪声指标 |
| 自旋螺旋 | `scf.Generation.Kpoint Regular`（4.0），RMM-DIISH | 每个 q 用相同参数 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 磁矩收敛到 0 | 初始占据上下相等，或初值太小、混合太激进 | 给明显不对称的初值；减小混合权重；先高电子温度再降温续算 |
| AFM 跑成 FM / NM | 初始占据没写反，或 SCF 中途翻转 | 核对两原子的第 6、7 列；减小 `scf.Max.Mixing.Weight`；用 `scf.restart on` 从部分收敛的 AFM 续算 |
| `Found an undefined species name ... Please check your input file` | species 名拼写不一致，或文件含全角字符/BOM/CRLF | 统一名称；`dos2unix`，确认文件为 ASCII |
| `Format error of the sequential number ... in <Atoms.SpeciesAndCoordinates` | 坐标块内写了整行注释或空行、序号不连续 | 块内只写数据行（说明性注释放到块外）；序号从 1 连续编号 |
| 初始电子数之和与 VPS 价电子数不符 | 第 6、7 列之和写错（如 Fe_PBE19H 写成 14） | 按 VPS 表核对（Fe_PBE19H = 16，Fe_PBE19S = 14） |
| `SpinP_switch should be OFF for this exchange functional.` | `scf.XcType LDA` 与自旋极化同时使用 | 改用 LSDA-CA / LSDA-PW / GGA-PBE |
| `Spin-orbit coupling is not supported for collinear DFT calculations.` | 共线 + SOC | 改 `scf.SpinPolarization nc` |
| `The constraint scheme is not supported for a collinear DFT calculation.` | 共线计算打开了约束 | 约束只用于 nc |
| `For spin magnetic moment, the constraint scheme and the Zeeman term are mutually exclusive.` | 自旋约束与自旋 Zeeman 同时打开 | 二选一 |
| `The Zeeman term for orbital moment is not supported without the SO term.` | 轨道 Zeeman 未开 SOC | `scf.SpinOrbit.Coupling on` |
| `Generalized Bloch theorem is not supported for the collinear DFT calculation.` / `Spin-orbit coupling is not supported for the generalized Bloch theorem.` | 螺旋计算设置错误 | nc + SOC off |
| 自旋螺旋报 `Fatal error in internal_Comm_size: Invalid communicator` | 4.0 默认 GCenter k 点与广义 Bloch 冲突 | `scf.Generation.Kpoint Regular` |
| `scf.restart=c2n is supported only for non-collinear calculations.` | 第二步忘了改为 nc | `scf.SpinPolarization nc` |
| c2n 第二步读不到重启文件 | `scf.restart.filename` 不是第一步的 System.Name，或目录不对 | 文件名写第一步的 System.Name（不带 `_rst`），在同一目录运行 |
| MAE 符号/大小随 k 点剧烈跳动 | k 点不足、电子温度太高 | 加密第二步 k 网格、降低电子温度；比较 Uele 而不是 Utot |
| 非共线 SCF 长时间振荡 | 自旋方向自由度能量尺度小 | 先约束后放开；RMM-DIISK 小权重 + 大 History；适当提高电子温度 |
| 非共线 + 仅 Γ（Cluster）时 HOMO 序号错误 | 4.0 原版 `Mulliken_Charge.c` 的问题 | 打 patch 4.0.1 |
| jx 结果同一壳层数值差异大 | `Num.Kgrid` 或 SCF k 点太少、基组过大 | 加密 k 点；按手册用较小基组 |
| jx 内存不足 | `Bunch.ij.pairs` 默认等于全部对数 | 调小 `Bunch.ij.pairs` |
| `jx` 用于非共线 scfout 得不到结果 | 4.0 的 jx 只支持共线 | 用共线 SCF；需要 DMI/各向异性时用 TB2J-OpenMX |

---

## 练习

1. **初值敏感性**：bcc Fe 分别用 (8.5, 7.5)、(10.0, 6.0)、(11.0, 5.0) 三种初值，比较收敛步数与最终磁矩；再换成 `Fe_PBE19S`（价电子 14），初值应如何改？
2. **体积–磁矩**：对 bcc Fe 做 a = 2.75–2.95 Å 的扫描，画磁矩随体积的变化（固定 `scf.Ngrid`，见 [第6章](06-optimization-eos-neb.md)）。
3. **Cr 的 AFM 稳定性**：用 `run_fm_afm_nm.sh` 比较三种磁态，并把电子温度从 300 K 提高到 1000 K，观察 AFM 磁矩的变化。
4. **MAE 收敛**：用 §9 的 FePt 输入，把第二步 k 网格从 12×12×9 逐步加密到 30×30×22，画 MAE 收敛曲线。
5. **J 的壳层收敛**：在 jx 配置中加入第 4、5 近邻壳层（自己推导 (l1,l2,l3)），观察 `jx_tc_mf.py` 给出的 T_C 如何变化。
6. **约束扫描**：用 `scan_theta.sh` 对 MnO 分子做 θ = 0°–90° 的整体转动扫描，分别在 `scf.Constraint.NC.Spin.v` = 0.5 和 1.0 eV 下比较能量曲线。

---

## 参考资料

**官方文档**

- OpenMX 4.0 手册 第 8.2 章 关键词（`Atoms.SpeciesAndCoordinates`、`scf.SpinPolarization`、约束）：<https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html>
- 第 37 章 Non-collinear DFT：<https://www.openmx-square.org/openmx_man4.0/s37_ncl-dft.html>
- 第 39 章 Orbital magnetic moment：<https://www.openmx-square.org/openmx_man4.0/s39_orb-magmom.html>
- 第 41 章 约束 DFT：<https://www.openmx-square.org/openmx_man4.0/s41_constraint-ncl.html>
- 第 42 章 二次变分法与 MAE：<https://www.openmx-square.org/openmx_man4.0/s42_mae.html>
- 第 43 章 Zeeman 项：<https://www.openmx-square.org/openmx_man4.0/s43_zeeman.html>
- 第 45 章 交换耦合参数（jx）：<https://www.openmx-square.org/openmx_man4.0/s45_jx.html>
- 第 57 章 自旋螺旋：<https://www.openmx-square.org/openmx_man4.0/s57_spin-spiral.html>
- 第 17 章 重启（`scf.restart on|pp|c2n`）：<https://www.openmx-square.org/openmx_man4.0/s17_restart.html>
- 第 15 章 k 点生成（GCenter/Regular）：<https://www.openmx-square.org/openmx_man4.0/s15_k-sampling.html>
- OpenMX 论坛：<https://www.openmx-square.org/forum/patio.cgi>
- TB2J 的 OpenMX 接口文档：<https://tb2j.readthedocs.io/en/latest/src/openmx.html>

**相关 GitHub 仓库（仅链接）**

- [mailhexu/TB2J](https://github.com/mailhexu/TB2J)：TB2J 主仓库，基于磁力定理与格林函数计算海森堡参数。
- [mailhexu/TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX)：TB2J 的 OpenMX 插件，从 scfout 计算 J、DMI、各向异性交换与磁振子能带。
- [KAIST-ELST/DFTforge.jl](https://github.com/KAIST-ELST/DFTforge.jl)：Julia 框架，读 OpenMX scfout 计算 J(q) 等，是 OpenMX Jx 方法的扩展实现。
- [pohao82/vasp-omx-moment-setter](https://github.com/pohao82/vasp-omx-moment-setter)：交互式设置 FM/AFM/非共线初始磁矩并写入 OpenMX 坐标块。
- [zeleznyj-htp/aiida-plugin-openmx](https://github.com/zeleznyj-htp/aiida-plugin-openmx)：较新的 AiiDA 插件，含 jx 交换常数工作链。
- [pranabdas/openmx](https://github.com/pranabdas/openmx)：英文 OpenMX 教程站，含磁矩可视化示例。
- [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-)：印尼语教程，含磁性、交换参数与居里温度、MAE 主题分支，可作覆盖面对照。

**关键文献**

- A. I. Liechtenstein, M. I. Katsnelson, V. P. Antropov, V. A. Gubanov, J. Magn. Magn. Mater. 67, 65 (1987). DOI: [10.1016/0304-8853(87)90721-9](https://doi.org/10.1016/0304-8853(87)90721-9)
- M. J. Han, T. Ozaki, J. Yu, Phys. Rev. B 70, 184421 (2004)（OpenMX 中的交换耦合计算）. DOI: [10.1103/PhysRevB.70.184421](https://doi.org/10.1103/PhysRevB.70.184421)
- A. Terasawa, M. Matsumoto, T. Ozaki, Y. Gohda, J. Phys. Soc. Jpn. 88, 114706 (2019)（jx 的有限极点算法）. DOI: [10.7566/JPSJ.88.114706](https://doi.org/10.7566/JPSJ.88.114706)
- T. Ozaki, Phys. Rev. B 75, 035123 (2007)（费米函数连分式/有限极点展开）. DOI: [10.1103/PhysRevB.75.035123](https://doi.org/10.1103/PhysRevB.75.035123)
- M. Pajda, J. Kudrnovský, I. Turek, V. Drchal, P. Bruno, Phys. Rev. B 64, 174402 (2001)（Fe、Co、Ni 的交换常数与 T_C）. DOI: [10.1103/PhysRevB.64.174402](https://doi.org/10.1103/PhysRevB.64.174402)
- X. He, N. Helbig, M. J. Verstraete, E. Bousquet, Comput. Phys. Commun. 264, 107938 (2021)（TB2J）. DOI: [10.1016/j.cpc.2021.107938](https://doi.org/10.1016/j.cpc.2021.107938)
- L. M. Sandratskii, Adv. Phys. 47, 91 (1998)（非共线磁性与广义 Bloch 定理综述）. DOI: [10.1080/000187398243573](https://doi.org/10.1080/000187398243573)
- J. Harris, Phys. Rev. B 31, 1770 (1985)（Harris 泛函）. DOI: [10.1103/PhysRevB.31.1770](https://doi.org/10.1103/PhysRevB.31.1770)

[← 上一章](10-charge-mo-cube.md) | [返回目录](../README.md) | [下一章 →](12-soc-topology.md)
