# 第8章 势能面扫描（刚性与柔性扫描）（Rigid & Relaxed PES Scans in Gaussian 16）

本章讲解如何用 Gaussian 16 做势能面扫描（potential energy surface scan, PES scan）：用 `Scan` 关键词配合 Z-matrix 做刚性扫描，用 `Opt=ModRedundant` 的 `S` 动作做柔性（relaxed）扫描，以及二维扫描和 G16 的广义内坐标（GIC）扫描。还会介绍如何提取扫描曲线、如何把能量最高点交给 [第7章](07-ts-irc.md) 的过渡态优化，并处理线性角报错、结构跳变、原子编号错位等常见问题。关键词：Gaussian relaxed scan、ModRedundant scan、rigid scan Z-matrix、dihedral scan、2D PES scan、GIC。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 区分刚性扫描和柔性扫描，知道各自的物理含义、代价和适用场合；
2. 用 Z-matrix + `Scan` 做刚性扫描（H₂O₂ 的 H–O–O–H 二面角）；
3. 用 `Opt=ModRedundant` 写出 `B 1 2 S 10 0.1`、`A 1 2 3 S 18 5.0`、`D 1 2 3 4 S 36 10.0` 形式的柔性扫描（正丁烷 C–C–C–C 二面角）；
4. 用键长扫描为过渡态搜索提供初猜（Cl⁻ + CH₃Br 的 C···Cl 距离）；
5. 写出二维扫描，理解计算量随网格点数增长的方式；
6. 用 `Opt=AddGIC` 定义“键长差”这类复合坐标并沿它扫描；
7. 从日志中提取能量曲线，在 GaussView 和 Python 中作图，并取出最高点结构。

**适用场景**：

| 目的 | 典型扫描 |
|---|---|
| 构象分析、转动能垒 | 二面角柔性扫描 |
| 为 TS 搜索提供初猜 | 沿“形成/断裂键长”做柔性扫描，取最高点 |
| 力场参数化（二面角项拟合） | 二面角柔性扫描（必要时冻结其他二面角） |
| 解离曲线、键能估计 | 键长扫描（共价键均裂需要开壳层处理，见 [第13章](13-openshell-stability.md)） |
| 光谱/振动的一维模型势 | 沿简正模式或单一内坐标的刚性扫描 |
| 二维反应势能面（如 More O'Ferrall–Jencks 图） | 两个键长的二维柔性扫描 |

**前置知识**：输入文件结构与 Z-matrix（[第2章](02-input-anatomy.md)）、几何优化与 ModRedundant 冻结（[第5章](05-opt.md)）。

本章配套示例：[../examples/08-pes-scan/README.md](../examples/08-pes-scan/README.md)

| 文件 | 内容 |
|---|---|
| `h2o2_rigid_zmatrix_scan.gjf` | 刚性扫描：H₂O₂ 二面角 0→360°，步长 10° |
| `butane_dihedral_relaxed_scan.gjf` | 柔性扫描：正丁烷 C1–C2–C3–C4 二面角，36 步 × 10° |
| `sn2_bond_scan_ts_guess.gjf` | 柔性扫描：Cl⁻ + CH₃Br 的 C···Cl 距离 3.2→1.8 Å，为 TS 提供初猜 |
| `two_dimensional_scan.gjf` | 二维柔性扫描：C···Cl 与 C–Br 两个距离，7 × 7 = 49 点 |
| `gic_bond_difference_scan.gjf` | GIC 扫描：反应坐标 RC = r(C–Br) − r(C–Cl) |

---

## 2. 原理速览

### 2.1 刚性扫描与柔性扫描

设扫描坐标为 $q$（某个键长、键角或二面角），其余自由度记为 $\mathbf{y}$：

- **刚性扫描（rigid scan）**：其余坐标固定在初始值 $\mathbf{y}_0$，每个点只算一次单点能：

$$
E_{\text{rigid}}(q) = E(q,\ \mathbf{y}_0)
$$

- **柔性扫描（relaxed scan）**：每个点上把 $q$ 固定在目标值，其余坐标做约束优化：

$$
E_{\text{relaxed}}(q) = \min_{\mathbf{y}} E(q,\ \mathbf{y})
$$

因为柔性扫描在每个点都取了最小值，所以总有 $E_{\text{relaxed}}(q) \le E_{\text{rigid}}(q)$。刚性扫描得到的能垒一般**偏高**，偏离平衡结构越远，误差越大。

| 对比项 | 刚性扫描 `Scan` | 柔性扫描 `Opt=ModRedundant` + `S` |
|---|---|---|
| 每点计算 | 1 次单点 | 1 次约束优化（通常 5–20 个 SCF+梯度） |
| 代价 | 低 | 高（约为刚性扫描的 5–20 倍） |
| 结构输入 | **必须用 Z-matrix**，被扫描的量是 Z-matrix 变量 | 任意（笛卡尔坐标最常用），被扫描的量是冗余内坐标 |
| 物理意义 | 其余结构冻结的截面，能垒偏高 | 最小能量路径的近似（绝热势能曲线） |
| 典型用途 | 快速预览、很小的分子、谱学模型势 | 转动能垒、TS 初猜、力场拟合 |

### 2.2 网格点数

每个扫描坐标写成“步数 $n$ + 步长 $\Delta$”，**点数是 $n+1$**（包括起点）。多个坐标同时扫描时形成矩形网格：

$$
N_{\text{points}} = \prod_{i}(n_i + 1)
$$

例如 `D 1 2 3 4 S 36 10.0` 共 37 个点；二维扫描 `S 6 ...` × `S 6 ...` 共 7 × 7 = 49 个点。二维柔性扫描很容易达到几百次约束优化，务必先估算计算量。

### 2.3 扫描曲线的最高点 ≠ 过渡态

柔性扫描把一个坐标强行固定，最高点只是“沿该坐标的约束最大值”。真实反应坐标往往是多个坐标的组合（例如 S<sub>N</sub>2 中 C–Nu 缩短的同时 C–LG 伸长），所以最高点只是 TS 的**初猜**，必须再做 `Opt=(TS,CalcFC)` + Freq + IRC（[第7章](07-ts-irc.md)）。如果反应坐标与扫描坐标耦合很弱，扫描曲线还可能出现**滞后**（正向和反向扫描曲线不重合）或突然跳变，这时需要换一个更合适的扫描坐标，例如第 7 节的键长差 GIC。

---

## 3. 刚性扫描：Scan + Z-matrix（H₂O₂ 二面角）

### 3.1 输入

完整输入 [`h2o2_rigid_zmatrix_scan.gjf`](../examples/08-pes-scan/h2o2_rigid_zmatrix_scan.gjf)：

```text
%chk=h2o2_rigid_scan.chk
%nprocshared=4
%mem=4GB
#p B3LYP/6-311+G(d,p) Scan NoSymm

H2O2 rigid scan of the H-O-O-H dihedral, 0 to 360 deg in 10 deg steps

0 1
H
O 1 ROH
O 2 ROO 1 AOOH
H 3 ROH 2 AOOH 1 DHOOH

ROH 0.968
ROO 1.452
AOOH 100.5
DHOOH 0.0 36 10.0

```

### 3.2 逐项解释

| 部分 | 说明 |
|---|---|
| `Scan` | 任务类型为刚性扫描，每个点只算单点能 |
| `NoSymm` | 二面角从 0°（C<sub>2v</sub>）经过 C<sub>2</sub> 到 180°（C<sub>2h</sub>），对称性会改变。官方文档要求扫描变量破坏对称性时必须加 `NoSymm`，否则作业会报错 |
| Z-matrix 部分 | 第 1 行 H；第 2 行 O 连到原子 1，键长 ROH；第 3 行 O 连到 2，键长 ROO，键角 ∠(O3–O2–H1)=AOOH；第 4 行 H 连到 3，二面角 H4–O3–O2–H1 = DHOOH |
| 空行 + 变量段 | 变量值写在坐标后面，二者之间用空行隔开 |
| `DHOOH 0.0 36 10.0` | **变量名 初值 步数 步长**：从 0.0° 开始，走 36 步，每步 +10°，共 37 个点（0°、10°、…、360°） |
| 其余变量 | 只有初值，扫描中保持不变 |

> [!IMPORTANT]
> 刚性扫描中其余几何参数是**固定**的，所以初值应取自同一理论水平的优化结构。示例中的 ROH、ROO、AOOH 是近似值，正式计算前应先在 B3LYP/6-311+G(d,p) 下优化 H₂O₂，再把优化值填进来。两个 O–H 键长共用变量 ROH、两个键角共用 AOOH，可以在扫描中保持分子的 C<sub>2</sub> 对称性。

步长单位由 `Units` 关键词控制，默认是 Å 和度。如果有多个变量都写了“步数 步长”，程序会对所有组合逐一计算，例如 R 走 3 步、A 走 2 步，共 4 × 3 = 12 个点。中断的刚性扫描可以用 `Scan=Restart` 从 chk 续算。

### 3.3 运行与结果

```bash
g16 < h2o2_rigid_zmatrix_scan.gjf > h2o2_rigid_zmatrix_scan.log
grep "SCF Done" h2o2_rigid_zmatrix_scan.log | wc -l        # 应为 37
grep -A 45 "Summary of the potential surface scan" h2o2_rigid_zmatrix_scan.log
```

刚性扫描结束时打印 `Scan completed.` 和汇总表，每行一个点：序号、扫描变量值和能量（示意输出，数值仅供参考）：

```text
 Summary of the potential surface scan:
  N      DHOOH        SCF
 ----  ---------  -----------
    1     0.0000  -151.53xxxx
    2    10.0000  -151.53xxxx
  ...
```

**如何判断结果合理**：H₂O₂ 的平衡二面角约为 110–120°，trans（180°）势垒约 1 kcal/mol，cis（0°）势垒约 7 kcal/mol（文献中高精度计算和光谱拟合的近似值）。刚性扫描的势垒会比柔性扫描略高，因为 O–O 键长和 OOH 键角在 cis/trans 构型下实际上会有所变化。

---

## 4. 柔性扫描：Opt=ModRedundant（正丁烷二面角）

### 4.1 输入

完整输入 [`butane_dihedral_relaxed_scan.gjf`](../examples/08-pes-scan/butane_dihedral_relaxed_scan.gjf)（坐标为 anti 构象，C1–C2–C3–C4 = 180°）：

```text
%chk=butane_scan.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ Opt=ModRedundant NoSymm

n-butane C1-C2-C3-C4 relaxed dihedral scan: 36 steps x 10 deg starting from anti

0 1
C      1.513585   -0.457577    1.164336
C      0.234199    0.182604    0.649763
C     -0.234194   -0.469009   -0.649781
C     -1.513580    0.171173   -1.164354
H      2.324021   -0.355442    0.435564
H      1.830723    0.023855    2.094695
H      1.365247   -1.523138    1.366635
H      0.407239    1.253047    0.488550
H     -0.546135    0.091927    1.414378
H     -0.407234   -1.539451   -0.488568
H      0.546140   -0.378331   -1.414396
H     -1.365242    1.236733   -1.366652
H     -1.830718   -0.310259   -2.094713
H     -2.324016    0.069037   -0.435581

D 1 2 3 4 S 36 10.0

```

坐标后面空一行，然后是 ModRedundant 输入段，段末再空一行。

### 4.2 ModRedundant 语法

每一行的格式为：

```text
[Type] N1 [N2 [N3 [N4]]] [Action [Params]]
```

**坐标类型**：缺省时按原子个数自动判断：1 个原子 = 笛卡尔坐标，2 个 = 键长，3 个 = 键角，4 个 = 二面角。也可以显式写类型字母：`X`（笛卡尔）、`B`（键长）、`A`（键角）、`D`（二面角）、`L`（线性弯曲）。

**动作代码**：

| 代码 | 含义 | 示例 |
|---|---|---|
| （无） | 添加该坐标（不改变其状态） | `3 8`：加入 3–8 之间的键（例如氢键），以便在输出中跟踪 |
| `A` | 激活之前被冻结的坐标 | `B 1 2 A` |
| `F` | 冻结该坐标 | `B 1 7 F`、`D 1 2 3 4 F` |
| `B` | 添加坐标并生成所有相关的键角、二面角 | `3 8 B` |
| `K` | 删除该坐标及所有包含它的坐标 | `* * * K` |
| `R` | 只删除该坐标本身 | `* * * * R` |
| `D` | 对该坐标对应的初始 Hessian 行列做数值微分 | `1 2 D` |
| `H diag` | 把该坐标的初始 Hessian 对角元改为给定值 | `1 2 H 0.55` |
| `S nsteps stepsize` | **柔性扫描**：每次把该坐标增加 stepsize，共 nsteps 次，每个点做约束优化 | `D 1 2 3 4 S 36 10.0` |

原子编号从 1 开始，**哑原子不计数**。星号 `*` 是通配符：`* *` 表示所有已定义的键，`3 *` 表示所有与原子 3 相连的键，`* 4 *` 表示以 4 为顶点的所有键角，`* 3 4 *` 表示绕 3–4 键的所有二面角。

本章最常用的三种写法（与本书其他章节一致）：

```text
B 1 2 S 10 0.1
A 1 2 3 S 18 5.0
D 1 2 3 4 S 36 10.0
```

| 行 | 含义 | 点数 |
|---|---|---|
| `B 1 2 S 10 0.1` | 键长 1–2 每步 +0.1 Å，走 10 步 | 11 |
| `A 1 2 3 S 18 5.0` | 键角 1–2–3 每步 +5°，走 18 步 | 19 |
| `D 1 2 3 4 S 36 10.0` | 二面角 1–2–3–4 每步 +10°，走 36 步（转一整圈） | 37 |

步长可以为负数，例如 `B 1 2 S 14 -0.1` 表示键长每步缩短 0.1 Å。扫描坐标以外的约束可以同时写在同一段中，例如再加一行 `D 5 1 2 3 F` 冻结甲基的二面角。

> [!NOTE]
> 柔性扫描只能用 Berny 算法（默认）。官方文档同样要求：扫描变量会破坏对称性时必须加 `NoSymm`。二面角扫描几乎总会破坏对称性，所以本章所有柔性扫描都加了 `NoSymm`。

### 4.3 运行与监控

```bash
g16 < butane_dihedral_relaxed_scan.gjf > butane_dihedral_relaxed_scan.log

# 已经完成的扫描点数（每个点收敛时打印一次）
grep -c "Stationary point found" butane_dihedral_relaxed_scan.log
# 当前进度：G16 在每一步的 "Step number" 行中给出所处的扫描点
grep "Step number" butane_dihedral_relaxed_scan.log | tail -n 1
```

`Step number` 行通常写成类似 `Step number   4 out of a maximum of  52 on scan point     7 out of    37` 的形式（示意），其中 maximum 是**每个扫描点**允许的最大优化步数，可以用 `Opt=(ModRedundant,MaxCycles=N)` 调整。

全部完成后，日志末尾会有柔性扫描的汇总表，按列列出每个点的扫描坐标值和能量（示意输出，数值仅供参考）：

```bash
grep -A 30 "Summary of Optimized Potential Surface Scan" butane_dihedral_relaxed_scan.log
```

```text
 Summary of Optimized Potential Surface Scan (add -158.0 to energies):
                           1         2         3         4         5
     Eigenvalues --    -0.4xxxx  -0.4xxxx  -0.4xxxx  -0.4xxxx  -0.4xxxx
           D1        -180.00000 -170.00000 -160.00000 -150.00000 -140.00000
```

在常见的 G16 输出中，这张表里标为 `Eigenvalues` 的一行其实是各扫描点的能量（Hartree，需要加上标题中括号里的常数），下面各行是扫描坐标和其他冗余内坐标的值。这只是输出格式的历史写法。如果拿不准，就以第 8 节的 awk 方法逐点提取。

> [!WARNING]
> 一定要确认**每个点**都收敛了：`grep -c "Stationary point found"` 的结果应等于点数（本例 37）。如果某个点出现 `Number of steps exceeded`，该点的能量不可信，需要加大 `MaxCycles`、减小步长，或对该点单独重算。

### 4.4 结果解读

正丁烷的扫描曲线应当呈现：anti（180°）为全局极小；gauche（约 ±65°）为局部极小，比 anti 高约 0.6–0.9 kcal/mol（气相实验 ΔH 约 0.67 kcal/mol）；±120° 附近是 anti/gauche 之间的重叠式能垒；0°（syn）是最高能垒，约 4.5–6 kcal/mol，具体数值依赖方法。这些参考值都是文献中的近似范围，用来检查曲线形状是否合理，不能代替你自己的计算结果。

曲线应当关于 180° 大致对称（−60° 与 +60° 两个 gauche 等价）。如果明显不对称，常见原因是某个点的甲基转到了另一个构象，或者某个点没有收敛。

---

## 5. 扫描找 TS 初猜：S<sub>N</sub>2 的 C···Cl 距离

### 5.1 输入

起点是 Cl⁻···CH₃Br 离子–偶极复合物（C···Cl 3.2 Å），每步把 C···Cl 缩短 0.1 Å，共 14 步，到 1.8 Å（接近 CH₃Cl 的平衡键长）。Br 不受约束，会随着 Cl 靠近自动离开。完整输入 [`sn2_bond_scan_ts_guess.gjf`](../examples/08-pes-scan/sn2_bond_scan_ts_guess.gjf)：

```text
%chk=sn2_scan.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ Opt=ModRedundant NoSymm

Cl- + CH3Br relaxed scan of the C...Cl distance, 3.2 -> 1.8 A, to obtain a TS guess

-1 1
C     0.000000    0.000000    0.000000
Cl    0.000000    0.000000   -3.200000
Br    0.000000    0.000000    1.990000
H     1.030000    0.000000   -0.325000
H    -0.515000    0.892006   -0.325000
H    -0.515000   -0.892006   -0.325000

B 1 2 S 14 -0.1

```

原子编号与 [第7章](07-ts-irc.md) 的 S<sub>N</sub>2 示例完全一致：1 = C，2 = Cl，3 = Br，4–6 = H。

### 5.2 怎么用扫描结果

1. 提取能量曲线（第 8 节），找到最高点，一般在 C···Cl ≈ 2.3–2.4 Å 附近；
2. 取出该点的结构：在 GaussView 中打开 `.log`，Results → Scan，点击最高点，再 File → Save 另存为 `.gjf`；或用仓库脚本 [`scan_irc_profile.py`](../scripts/scan_irc_profile.py) 报告最高点，并用 `--write-max` 把最高点结构直接写成 xyz（`python ../scripts/scan_irc_profile.py -h` 查看用法），也可以配合 [`log2xyz.py`](../scripts/log2xyz.py) 导出结构；
3. 把这个结构作为初猜，运行 `Opt=(TS,CalcFC,NoEigenTest) Freq`，然后做 IRC（[第7章](07-ts-irc.md) 第 4–7 节）。

> [!TIP]
> 曲线最高点附近如果很平，可以在最高点两侧各取一个点，用更小的步长（0.02–0.05 Å）再细扫一次。TS 初猜越接近真实鞍点，`Opt=TS` 越容易一次成功。

### 5.3 扫描方向与滞后

从反应物一侧扫（C···Cl 缩短）和从产物一侧扫（C–Br 缩短）得到的曲线不一定重合。当扫描坐标不是好的反应坐标时，另一个“没被约束”的坐标可能在某一点突然跳变，曲线上出现台阶或尖峰。处理办法：

- 两个方向都扫一遍，比较最高点结构；
- 改用二维扫描（第 6 节），或者用两个键长之差作为扫描坐标（第 7 节 GIC）；
- 减小步长，让每一步的结构变化更连续。

---

## 6. 二维扫描

### 6.1 输入

同时写两行 `S`，就得到二维网格。完整输入 [`two_dimensional_scan.gjf`](../examples/08-pes-scan/two_dimensional_scan.gjf) 的 ModRedundant 段为：

```text
B 1 2 S 6 -0.2
B 1 3 S 6 0.2

```

起始结构中 C···Cl = 3.0 Å、C–Br = 2.0 Å。C···Cl 从 3.0 缩短到 1.8 Å，C–Br 从 2.0 伸长到 3.2 Å，共 7 × 7 = 49 个约束优化。每个点上两个距离都固定，只有三个 H（CH₃ 伞形）弛豫，所以每个点收敛很快。

### 6.2 二维扫描的注意事项

| 问题 | 建议 |
|---|---|
| 计算量 | 点数为 (n₁+1)(n₂+1)。先用粗网格（0.2 Å 或 30°）看全貌，再在关注区域加密 |
| 点的顺序 | 不同坐标谁先变化、谁后变化，以日志汇总表中每个点列出的坐标值为准，不要凭假设对号入座 |
| 网格角落的高能区 | 例如 C···Cl 和 C–Br 都很短的角落，能量很高，SCF 可能难收敛。若只关心反应通道，可以缩小范围 |
| 可视化 | Gaussian 16 C.01 起会把多维柔性扫描信息写入 fchk，可在 GaussView 中显示；仓库脚本 [`scan_irc_profile.py`](../scripts/scan_irc_profile.py) 能识别二维柔性扫描并画等高线图（`python ../scripts/scan_irc_profile.py -h` 查看用法） |
| 二面角二维扫描 | 例如丙氨酸二肽的 φ/ψ（Ramachandran 图），30° 网格就有 12 × 12 = 144 个点，建议先用小基组或 xtb 预扫 |

二维 S<sub>N</sub>2 势能面上，左上角（C···Cl 长、C–Br 短）是反应物谷，右下角（C···Cl 短、C–Br 长）是产物谷，鞍点位于两谷之间的“山口”处。这就是 More O'Ferrall–Jencks 图的思路，可以直观看出反应是协同的还是分步的。

---

## 7. GIC：沿“键长差”扫描（G16 B.01 起）

### 7.1 为什么需要 GIC

S<sub>N</sub>2 反应的真实反应坐标近似是“C–Br 伸长 + C–Cl 缩短”的反对称组合。普通 ModRedundant 只能扫描单个键长、键角或二面角，而 Gaussian 16 的**广义内坐标**（generalized internal coordinates, GIC）允许用数学表达式定义任意坐标，并用于约束和扫描。GIC 在 G16 Rev. B.01 中引入，C.01/C.02 文档中有完整的语法说明（Opt 和 Geom 关键词页面的 “GIC Info” 部分）。

### 7.2 输入

完整输入 [`gic_bond_difference_scan.gjf`](../examples/08-pes-scan/gic_bond_difference_scan.gjf)：

```text
%chk=sn2_gic_scan.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ Opt=AddGIC NoSymm

Cl- + CH3Br relaxed scan along the GIC RC = r(C-Br) - r(C-Cl), in Angstrom

-1 1
C     0.000000    0.000000    0.000000
Cl    0.000000    0.000000   -3.200000
Br    0.000000    0.000000    1.990000
H     1.030000    0.000000   -0.325000
H    -0.515000    0.892006   -0.325000
H    -0.515000   -0.892006   -0.325000

RC(NSteps=12,StepSize=0.2)=(R(1,3)-R(1,2))*0.529177

```

### 7.3 语法要点

| 写法 | 含义 |
|---|---|
| `Opt=AddGIC` | 先自动生成默认的 GIC 坐标，再读取坐标后面的 GIC 输入段进行增删改（等价于 `Opt=(ModRedundant,GIC)`） |
| `Opt=ReadAllGIC` | 不自动生成任何坐标，全部由用户定义（高级用法） |
| `Label(Options)=Expression` | 定义一个带标签的坐标。标签必须以字母开头，不区分大小写 |
| `R(i,j)`、`A(i,j,k)`、`D(i,j,k,l)` | 键长、键角、二面角函数（也有 `B(i,j)`/`Bond(i,j)` 等同义写法） |
| `NSteps=N,StepSize=x` | 柔性扫描：把该坐标增加 x，共 N 次。写在表达式后面（Expression Options 形式）时，逗号改为空格 |
| `Freeze` / `Frozen` | 冻结坐标，例如 `R(5,9) Freeze` 或 `HOH(Freeze)=A(2,1,3)` |
| `Inactive` / `Remove` | 不把该坐标用作优化变量（可作为中间量） |
| `Value=x` | 设定坐标初值（程序会调整笛卡尔坐标去匹配，官方建议谨慎少用） |
| `Min=a,Max=b` | 条件动作：坐标值在区间内时才激活/冻结/移除 |

> [!WARNING]
> **单位问题**：纯键长、键角、二面角在输出中以 Å 和度显示。但只要坐标是一个“表达式”（例如两个键长之差），它就是 *generic GIC*，程序内部用 Bohr 和弧度计算，输出注明 `GIC-type coordinates are in arbitrary units`。示例中乘以 0.529177（1 Bohr = 0.529177 Å），把 RC 换算成 Å，这与官方文档中 `F1F2=...*0.529177` 的写法相同，StepSize=0.2 也按这个坐标自身的单位理解。需要说明的是，官方 Opt 页面称 StepSize 采用 `Units` 关键词定义的单位（默认 Å/度），而 Geom 页面的 GIC Info 写的是“原子单位”，两处说法并不完全一致。因此第一次使用时，建议先把 `NSteps` 设为 2 试跑，核对输出 `Initial Parameters` 表中 RC 的初值（本例应约为 1.99 − 3.20 = −1.21）以及前两个扫描点的 RC 是否按预期递增，再正式提交。

> [!NOTE]
> 官方文档指出：由多个二面角组合成的复合坐标在作为**活动**优化变量时可能表现不佳（周期性处理问题），但用于冻结或扫描是可以的；键长、键角组合成的坐标一般表现良好。另外，旧式 ModRedundant 写法和新的 GIC 写法不能混在同一个输入段中。

### 7.4 用途扩展

- **约束优化**：`RC(Freeze)=(R(1,3)-R(1,2))*0.529177`，在固定反应坐标的同时弛豫其余自由度，可以替代第 7 章 3.5 节的双键冻结；
- **碎片质心距离**：官方 GIC 文档演示了用 `XCntr/YCntr/ZCntr` 定义两个碎片的质心，再扫描质心距离，适合弱相互作用复合物的解离曲线；
- **线性弯曲**：`L(i,j,k,L,M)` 处理接近 180° 的角，避免线性角问题。

---

## 8. 扫描结果的提取与作图

### 8.1 GaussView

打开 `.log` 文件，Results → Scan，会画出能量随扫描坐标变化的曲线。点击曲线上的点可以查看对应结构；在结构窗口中 File → Save 可以把它另存为新的输入。二维扫描建议先 `formchk` 生成 `.fchk` 再打开（`formchk` 的用法见 [第16章](16-postprocess-multiwfn.md)）。

### 8.2 命令行快速提取

柔性扫描中，每个点收敛前的最后一个 `SCF Done` 就是该点的能量。下面的 awk 命令按顺序列出每个收敛点的能量，并换算为相对最低点的 kcal/mol：

```bash
LOG=butane_dihedral_relaxed_scan.log
awk '/SCF Done/{e=$5} /Stationary point found/{n++; printf "%3d %18.8f\n", n, e}' "$LOG" > scan_E.dat
awk 'NR==FNR{if(min==""||$2<min)min=$2; next}{printf "%3d %10.3f\n", $1, ($2-min)*627.509}' scan_E.dat scan_E.dat
```

如果要同时取出每个点的二面角数值，可以在 `Stationary point found` 之后的 `Optimized Parameters` 表中查找 `D(1,2,3,4)` 所在的行：

```bash
awk '/Stationary point found/{f=1} f && /D\(1,2,3,4\)/{print $4; f=0}' "$LOG"
```

对于刚性扫描，每个点只有一个 `SCF Done`，直接 `grep "SCF Done"` 即可，或者读取 `Summary of the potential surface scan` 表。

### 8.3 Python 作图

仓库脚本 [`scan_irc_profile.py`](../scripts/scan_irc_profile.py) 可以解析一维/二维柔性扫描和 IRC，以 kcal/mol 相对最低点作图，并报告最高点作为 TS 初猜：

```bash
python ../scripts/scan_irc_profile.py -h
```

如果只想快速看一眼曲线，也可以把上面的 `scan_E.dat` 直接交给 gnuplot、Excel 或 matplotlib。出版级作图规范见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md)。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 说明 |
|---|---|---|
| 扫描类型 | 默认用柔性扫描；刚性扫描只用于快速预览或小分子 | 刚性能垒偏高 |
| 起始结构 | 先在同一水平优化好的极小点 | 第一个点就不收敛会拖慢整个扫描 |
| 二面角步长 | 10–15°（预览可用 30°） | 步长过大容易漏掉极值，也容易让其他转子跳变 |
| 键长步长 | 0.05–0.1 Å；最高点附近细扫 0.02 Å | 键长变化对能量影响大 |
| 键角步长 | 2–5° | 避免扫过 180° |
| 每点最大步数 | `Opt=(ModRedundant,MaxCycles=50)` 左右 | 默认值通常够用，平坦体系可适当加大 |
| 对称性 | 一律加 `NoSymm` | 防止扫描中对称性改变导致报错 |
| 理论水平 | 预扫可用小基组/半经验/xtb，重要区域再用目标水平 | 二维扫描尤其要控制成本 |
| 解离曲线 | 均裂用 UKS + `Guess=Mix`，并检查稳定性 | 闭壳层 RKS 在长键长处给出错误的渐近线 |
| 力场拟合 | 冻结其他软二面角，或报告它们的变化 | 避免多个扭转耦合导致曲线难以拟合 |

---

## 常见问题与报错

### 1. `FormBX had a problem` / `Linear angle in Bend` / `Linear angle in Tors`（L103）

**原因**：冗余内坐标在接近 180° 的键角处变得奇异。扫描中某个角被推过（或接近）180°，或者体系中本来就有近线性片段（炔、腈、M–C≡O，S<sub>N</sub>2 的 Nu–C–LG）。

**处理**：

- 避免让扫描坐标穿过 180°，例如把键角扫描终点设在 175°；
- 从出错前最后一个点的结构重新开始，并缩小步长；
- 改用 `Opt=(ModRedundant,Cartesian)` 在笛卡尔坐标下优化（注意纯笛卡尔优化不能冻结内坐标，约束需求复杂时不适用）；
- 用 GIC 的线性弯曲坐标 `L(i,j,k,L,M)` 代替普通键角；
- 必要时加哑原子（Z-matrix 输入）。

更多细节见 [第18章](18-errors-capstone.md)。

### 2. 扫描曲线中途突然跳变

**现象**：曲线在某一点突然下降或上升几 kcal/mol，前后结构差别明显（例如甲基转了 120°、羟基氢键换了方向、另一个二面角翻转）。

**原因**：其他软自由度在约束优化中落到了另一个局部极小；或者扫描坐标与真实反应坐标耦合很弱（滞后）。

**处理**：减小步长；从跳变后的结构反向再扫一遍，取两条曲线在每个点上的较低值；冻结不相关的软二面角；或者换用更好的扫描坐标（二维扫描、GIC 组合坐标）。

### 3. 约束的原子编号对不上

- ModRedundant 和 GIC 中的编号是**分子说明部分的原子顺序**，从 1 开始，哑原子不计，与 GaussView 中显示的编号一致；
- 在 GaussView 中增删原子、或者用别的程序转换格式后，编号可能整体改变。提交前在 GaussView 中逐一核对扫描坐标涉及的原子；
- 带 `Geom=Connectivity` 的 GaussView 输入在原子数改变后，连接表也必须同步更新，否则读入时报错（见 [第2章](02-input-anatomy.md)）。

### 4. 扫描破坏对称性后报错

官方文档明确说明：扫描变量破坏分子对称性时，作业会出错，除非 route 中加了 `NoSymm`。给所有扫描作业都加上 `NoSymm` 即可避免。

### 5. 刚性扫描报错或扫描变量没有变化

`Scan` 刚性扫描要求分子用 **Z-matrix** 定义，被扫描的必须是 Z-matrix 中的变量（写在变量段中，并带“步数 步长”）。用笛卡尔坐标输入再加 `Scan` 达不到刚性扫描的效果，这种情况应改用柔性扫描。

### 6. 某些点没有收敛（`Number of steps exceeded`）

用 `grep -c "Stationary point found"` 与点数对比。没有收敛的点能量不可信：加大 `MaxCycles`；从该点前一个收敛结构出发，用 `Opt=ModRedundant` 单独重算该点；或者在该区域减小步长。判断优化是振荡还是缓慢下降，可以用 `python ../scripts/opt_monitor.py -h` 查看用法。

### 7. 键解离曲线在长键长处能量不对

闭壳层 RKS/RHF 描述均裂时会给出过高的解离极限。使用 `UB3LYP ... Guess=Mix` 并配合 `Stable=Opt` 检查波函数稳定性，得到 broken-symmetry 解，详见 [第13章](13-openshell-stability.md)。异裂（得到离子对）则通常需要溶剂模型才能得到有意义的曲线。

---

## 练习

1. 对 H₂O₂ 分别做刚性扫描和柔性扫描（`D 1 2 3 4 S 36 10.0`，用笛卡尔坐标输入），把两条曲线画在一张图上，比较 cis 和 trans 势垒的差别。
2. 把正丁烷扫描的步长改为 30°，看能否正确找到 gauche 极小的位置和能量；再在 gauche 附近用 5° 细扫。
3. 从产物侧（CH₃Cl···Br⁻ 复合物）反向扫描 C···Br 距离，与本章 C···Cl 扫描的最高点结构比较。
4. 用 `two_dimensional_scan.gjf` 的结果画出二维等高线图，标出反应物谷、产物谷和鞍点区域。
5. 用 GIC 的 `RC(Freeze)=...` 写法，把 RC 固定在二维扫描的鞍点附近做约束优化，再交给 `Opt=(TS,CalcFC)`。

---

## 参考资料

**官方文档**

- [Scan 关键词（刚性扫描，Z-matrix 变量步数与步长）](https://gaussian.com/scan/)
- [Opt 关键词（ModRedundant 语法、动作代码、通配符、柔性扫描、GIC 选项）](https://gaussian.com/opt/)
- [Geom 关键词（ModRedundant、Checkpoint、GIC Info：GIC 的完整语法与单位说明）](https://gaussian.com/geom/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes（多维柔性扫描写入 fchk 等）](https://gaussian.com/relnotes/)
- [GaussView 6](https://gaussian.com/gaussview6/)

**相关 GitHub 仓库**（只链接，不复制代码）

- [cclib/cclib](https://github.com/cclib/cclib)：可解析 Gaussian 日志中的扫描能量与几何，适合交叉验证自写脚本
- [chunxiangzheng/gaussian_log_file_converter](https://github.com/chunxiangzheng/gaussian_log_file_converter)：从 log 提取结构转 xyz/gjf
- [maserasgroup-repo/pyssian](https://github.com/maserasgroup-repo/pyssian)：按 Link 分块解析 Gaussian 输出的 Python 库
- [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian)：用 xtb 做低成本预扫描
- [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py)：结构操作与输出解析工具包
- 更多工具见 [../../Resources/github-repos.md](../../Resources/github-repos.md)

**本仓库脚本**：[`../scripts/scan_irc_profile.py`](../scripts/scan_irc_profile.py)（扫描与 IRC 能量曲线，报告最高点）、[`../scripts/log2xyz.py`](../scripts/log2xyz.py)（导出所有扫描点结构为多帧 xyz）、[`../scripts/opt_monitor.py`](../scripts/opt_monitor.py)（单点优化收敛诊断）、[`../scripts/gjf_builder.py`](../scripts/gjf_builder.py)（生成 scan 模板）。各脚本用 `-h` 查看用法。

**关键文献**

- C. Peng, P. Y. Ayala, H. B. Schlegel, M. J. Frisch, Using redundant internal coordinates to optimize equilibrium geometries and transition states, *J. Comput. Chem.* **1996**, 17, 49–56.
- J. B. Foresman, Æ. Frisch, *Exploring Chemistry with Electronic Structure Methods*, 3rd ed., Gaussian, Inc., 2015（第 6 章讨论势能面扫描）。
- R. A. More O'Ferrall, Relationships between E2 and E1cB mechanisms of β-elimination, *J. Chem. Soc. B* **1970**, 274–277（二维反应势能面 More O'Ferrall–Jencks 图的出处之一）。

---

[← 上一章](07-ts-irc.md) | [返回目录](../README.md) | [下一章 →](09-solvation.md)
