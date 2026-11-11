# 第5章 SCF 参数、k 点与收敛测试（SCF Settings, k-Point Sampling & Convergence Tests in OpenMX）

本章系统讲解 OpenMX 自洽场（SCF）计算的各个旋钮：交换关联泛函 `scf.XcType`、自旋 `scf.SpinPolarization`、电子温度、实空间网格截断能 `scf.energycutoff` 与 `scf.Ngrid`、本征值求解器与 k 点网格（含 4.0 新增的 `scf.Generation.Kpoint GCenter|Regular`）、七种混合方案（RMM-DIISK 等）、运行中调参文件与 `scf.restart on|pp`。最后用 Si 和 bcc Fe 两个可运行算例，完成截断能与 k 点收敛测试（OpenMX cutoff energy and k-point convergence test, SCF mixing troubleshooting）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据体系选择 `scf.XcType`、`scf.SpinPolarization`、`scf.ElectronicTemperature`，知道 OpenMX 4.0 不提供杂化泛函；
- 理解 `scf.energycutoff` 控制的是实空间积分网格而不是基组，会读 `.out` 中的网格信息，知道何时必须用 `scf.Ngrid` 固定网格；
- 选择 `scf.EigenvalueSolver` 和 `scf.Kgrid`，理解 GCenter 与 Regular 两种 k 点生成方式的差别，以及它对复现 3.9 结果的影响；
- 按"先 RMM-DIISK、再调历史步数和起始步、金属提温度、+U 换 RMM-DIISH"的顺序处理难收敛体系；
- 在计算运行中通过 `System.Name_SCF_keywords` 文件修改混合参数，或优雅地终止计算；
- 用 `scf.restart on` 续算，用 `scf.restart pp` 只做一次对角化来算能带和 DOS；
- 完成截断能和 k 点收敛测试，并按 1 meV/atom 的判据给出推荐值。

| 场景 | 本章对应小节 |
|---|---|
| 新体系第一次计算，参数怎么填 | §3–§6、"参数选择与收敛建议" |
| E–V 曲线有锯齿、能量随结构跳变 | §4.3–§4.4 |
| 4.0 与 3.9 的结果对不上 | §5.3 |
| 金属、磁性、表面体系 SCF 振荡 | §6、§10 |
| 跑了很久还没收敛，不想重启 | §7 |
| 已有自洽电荷，只想补算能带/DOS | §8 |
| 写论文前的收敛性测试 | §9 |

> [!NOTE]
> 输入文件的整体结构和第一个甲烷计算见 [第3章](03-input-first-calc.md)；基组（PAO/VPS）怎么选见 [第4章](04-pao-vps.md)。本章做收敛测试时固定基组，只改截断能和 k 点。

---

## 2. 原理速览

### 2.1 SCF 循环在做什么

Kohn–Sham 方程的哈密顿量依赖于电子密度，而密度又由哈密顿量的本征态决定。OpenMX 的每个 SCF 步：

1. 由输入的密度矩阵 $\rho^{\text{in}}_n$（或电子密度、势、哈密顿量，取决于混合方案）构造哈密顿量；
2. 在每个 k 点上求解广义本征值问题 $H\,c=\varepsilon\,S\,c$（LCPAO 基组不正交，$S$ 为重叠矩阵）；
3. 按费米–狄拉克分布填充占据数，得到输出 $\rho^{\text{out}}_n$；
4. 用混合方案由历史上的输入、输出构造下一步的输入。

残差 $R_n=X^{\text{out}}_n-X^{\text{in}}_n$ 的范数在 `.out` 中记为 `NormRD`；能带能量（band energy）记为 `Uele`。OpenMX 的收敛判据是相邻两步 `Uele` 之差的绝对值 `dUele` 小于 `scf.criterion`（单位 Hartree）。

### 2.2 截断能不是基组截断

OpenMX 的基组是 PAO（第4章）。`scf.energycutoff` 决定的是**实空间积分网格**的疏密，用于计算差分 Hartree 势、交换关联势的矩阵元，以及用 FFT 解 Poisson 方程。因此：

- 截断能提高时，总能**不会**像平面波那样从上方单调收敛，而是在收敛值附近小幅振荡；
- 网格固定在晶胞上，原子相对网格的位置不同，能量会有微小差别（egg-box effect）；
- 网格点数随晶格常数离散变化，做 E–V 曲线时会产生跳变（§4.3）。

网格间距 $\Delta$（bohr）与截断能 $E_{\text{cut}}$（Ry）满足 $E_{\text{cut}}=(\pi/\Delta)^2$。沿晶格矢量 $\mathbf a_i$ 取 $N_i$ 个网格点时，对应的截断能为

$$
E_{\text{cut},i}=\left(\frac{\pi N_i}{d_i}\right)^2,\qquad d_i=\frac{\Omega}{|\mathbf a_j\times\mathbf a_k|},
$$

$d_i$ 是与 $\mathbf a_i$ 对应的晶面间距（bohr）。程序会选取满足 $E_{\text{cut},i}\ge$ `scf.energycutoff` 的 $N_i$，并把实际值打印出来（§4.2）。

### 2.3 布里渊区采样与电子温度

周期体系的总能、密度要对布里渊区积分，用均匀网格 $N_1\times N_2\times N_3$ 近似。绝缘体的被积函数光滑，较稀的网格就能收敛；金属的占据数在费米面处突变，需要更密的网格，并配合费米–狄拉克展宽（`scf.ElectronicTemperature`）把突变抹平：

$$
f(\varepsilon)=\frac{1}{1+\exp[(\varepsilon-\mu)/k_BT_e]},\qquad k_B\times 300\,\text{K}\approx 25.9\ \text{meV}.
$$

电子温度越高，SCF 越稳定，k 点收敛越快，但会改变占据数，影响磁矩和细小的能量差。

---

## 3. 交换关联、自旋与电子温度

### 3.1 `scf.XcType`

| 取值 | 含义 | 配套 VPS |
|---|---|---|
| `LDA` | Ceperley–Alder 局域密度近似（非自旋） | `*_CA19` |
| `LSDA-CA` | Ceperley–Alder 的自旋极化版本 | `*_CA19` |
| `LSDA-PW` | Perdew–Wang 参数化的局域自旋密度近似 | `*_CA19` |
| `GGA-PBE` | Perdew–Burke–Ernzerhof 广义梯度近似 | `*_PBE19` |

OpenMX 4.0 只提供以上四种，**没有杂化泛函和 meta-GGA**。需要更好的带隙时，常用 DFT+U（[第13章](13-dft-u.md)）；范德华作用用 DFT-D2/D3 补充（[第6章 §6](06-optimization-eos-neb.md)）。手册说明 GGA 的梯度用实空间 27 点多项式拟合计算。

> [!IMPORTANT]
> 选 `LDA` 时 `scf.SpinPolarization` 必须为 `off`；要做自旋极化的局域密度计算，请用 `LSDA-CA` 或 `LSDA-PW`。

### 3.2 `scf.SpinPolarization`

| 取值 | 含义 | 坐标行第 6、7 列 |
|---|---|---|
| `off` | 非自旋极化 | 上下各一半，如 Si `2.0 2.0` |
| `on` | 共线自旋极化 | 用上下差设定初始磁矩，如 Fe（16 e）`9.0 7.0` |
| `nc` | 非共线（SOC 也必须用它） | 还要加 Euler 角等列，见 [第11章](11-spin-noncollinear.md) |

磁性体系的初始磁矩对 SCF 收敛和最终磁态影响很大，详见第11章。

### 3.3 `scf.ElectronicTemperature`

默认 300 K。建议：

| 体系 | 推荐值 | 说明 |
|---|---|---|
| 分子、宽带隙绝缘体 | 300 K（默认） | 对结果几乎无影响 |
| 半导体 | 300 K | — |
| 普通金属 | 300–1000 K | 手册指出金属体系温度太低常出现数值不稳定 |
| 难收敛的磁性金属、大表面 | 先 1000 K 收敛，再用 `scf.restart on` 降温 | 降温后磁矩、能量要用最终温度的结果 |
| 磁各向异性能（MAE）等微小能量差 | 尽量低，并加密 k 点 | 手册第42章强调低温和密 k 点 |

> [!WARNING]
> 比较能量差时，所有参与比较的计算必须用同一个电子温度。

---

## 4. 实空间网格：`scf.energycutoff` 与 `scf.Ngrid`

### 4.1 默认值说法不一，务必显式设置

| 来源 | `scf.energycutoff` 默认值 |
|---|---|
| OpenMX 3.9 手册关键词说明 | 150 Ry |
| OpenMX 4.0 手册关键词说明 | 220 Ry |
| 官方示例输入（如 `Methane.dat`）中的注释 | 仍写 150 Ry |

因为不同版本、不同文档的说法不一致，**每个输入都应显式写出 `scf.energycutoff`**，不要依赖默认值。

手册给出的经验：多数情况下 150–200 Ry 就是合适的选择，很多体系在 200 Ry 附近收敛；但在势能面很平的结构优化、不同自旋序之间的微小能量差等问题上，可能需要 300 Ry 以上。Ozaki 的实践讲义还提到，含较深半芯态的赝势（如 `Fe_PBE19H` 这类 Hard 赝势）可能需要 300–400 Ry。

### 4.2 读懂网格信息

每次计算开始时，`.out`（以及标准输出）都会打印实际使用的网格：

```bash
grep -E "Required cutoff energy|Used cutoff energy|Num. of grids of" si_scf.out
```

输出形如（**示意输出，数值仅供参考**）：

```text
  Required cutoff energy (Ryd) for 3D-grids = 220.0000
      Used cutoff energy (Ryd) for 3D-grids = 2xx.xxxx, 2xx.xxxx, 2xx.xxxx
  Num. of grids of a-, b-, and c-axes = NN, NN, NN
```

`Used cutoff energy` 是三个方向实际达到的截断能，总是不低于 `Required`；`Num. of grids` 是 $N_1, N_2, N_3$。比较两个计算时，先看网格是否一致。

### 4.3 `scf.Ngrid`：固定网格点数

```text
scf.Ngrid     32 32 32      # N1 N2 N3 along a, b, c; overrides scf.energycutoff
```

一旦写了 `scf.Ngrid`，程序直接按给定点数划分三个晶格矢量，不再参考 `scf.energycutoff`。以下情况必须用它：

- **E–V 曲线、EvsLC 扫描**：晶格常数变化时网格点数会跳变，E–V 曲线出现锯齿（手册以 bcc Fe 为例说明了这一现象）。固定 `scf.Ngrid` 后曲线变光滑；
- **一组需要相减的计算晶胞形状不同**（例如应变扫描）；
- **超胞**：让 $N_i$ 能被超胞在该方向的扩胞倍数整除，网格才与原胞的平移对称性一致。论坛中也有自动选取的网格在三个方向严重不均衡、改用 `scf.Ngrid` 后结果恢复正常的报告。

选 `scf.Ngrid` 的简单办法：先按目标截断能跑一次，读出 `Num. of grids`，再在 E–V 扫描中把它固定（或取稍大一点的数），这样扫描中最"拉伸"的结构也不低于目标截断能。

### 4.4 `scf.fixed.grid`：固定网格原点

总能依赖于原子相对网格的位置。计算一系列只有少数原子移动的结构（结合能曲线、CP 修正中的单体、表面吸附位点比较）时，手册建议固定网格原点：

```bash
grep Grid_Origin first_run.std          # Grid_Origin is printed to standard output
```

把读到的三个数写进其余所有输入：

```text
scf.fixed.grid    xxx  yyy  zzz         # same Grid_Origin for every calculation in the series
```

结构优化、分子动力学生成的续算文件 `*.dat#` 中会自动写入 `scf.fixed.grid`（见第6章）。

### 4.5 `1DFFT.EnergyCutoff`

这是 PAO 和非局域投影子径向函数做一维傅里叶变换用的截断，默认 3600 Ry，手册说明默认值已充分收敛、代价也小，一般不用改。

---

## 5. 本征值求解器与 k 点

### 5.1 `scf.EigenvalueSolver`

| 取值 | 适用体系 | 说明 |
|---|---|---|
| `Cluster` | 分子、团簇 | 不用 k 点；盒子要足够大（第4章 §5.3） |
| `Band` | 周期体系（体相、表面、一维链） | 需要 `scf.Kgrid`；大超胞也可用 `1 1 1` |
| `DC` / `DC-LNO` / `Krylov` | 上千原子的大体系 | O(N) 方法，见 第18章 |
| `ON2` | 中大体系 | 数值精确的低阶标度方法 |
| `NEGF` | 输运 | 见 [第17章](17-negf-transport.md) |

大规模计算时还可以用 `scf.eigen.lib` 选择本征值求解库（ELPA1/ELPA2 等），见第18章。

### 5.2 `scf.Kgrid`：网格密度怎么定

`scf.Kgrid n1 n2 n3` 沿三个倒格矢均匀划分。常用经验是让 $n_i\cdot|\mathbf a_i|$（单位 Å）大致相等：

| 体系 | n_i × a_i 参考值（a_i 为晶格矢量长度） | Si 原胞（a_i = 3.84 Å）对应 |
|---|---|---|
| 宽带隙绝缘体、分子晶体 | ≈ 20–30 Å | 6×6×6 – 8×8×8 |
| 半导体 | ≈ 30–40 Å | 8×8×8 – 10×10×10 |
| 金属、磁性金属 | ≥ 50–60 Å | 14×14×14 以上 |
| 表面 / 二维材料 | 面内同上，真空方向取 1 | — |
| 能带、DOS 后处理 | 用 `Band.kpath` 和 `Dos.Kgrid`，不靠 `scf.Kgrid` | 第8、9章 |

（`pohao82/poscar2openmx` 等工具的默认规则是 $n_i\cdot d_i\ge 35$ Å，可作为起点。）这些只是起点，最终必须做 §9 的收敛测试。计算量大致与 k 点数成正比。

### 5.3 4.0 新增：`scf.Generation.Kpoint GCenter | Regular`

```text
scf.Generation.Kpoint    GCenter     # GCenter|Regular, default=GCenter (OpenMX >= 4.0)
```

| 方式 | 版本 | 生成规则（手册第15章） | 是否含 Γ |
|---|---|---|---|
| `GCenter` | 4.0 默认 | 以 Γ 为中心对称地布点；手册说明它无论 $N_i$ 奇偶都包含 Γ，并按构造包含时间反演不变动量（TRIM）点 | 总是包含 |
| `Regular` | 3.9 的唯一方式 | $q_{il}=-\tfrac12+\tfrac{2l+1}{2N_i}$，$l=0\ldots N_i-1$，即 Monkhorst–Pack 型网格 | 仅当 $N_i$ 全为奇数时；不含布里渊区边界点 |

实用要点：

- 日常计算用默认 `GCenter` 即可，手册也推荐它；
- **复现 3.9 的结果**（或与 3.9 时代的文献比较）时，显式写 `scf.Generation.Kpoint Regular`。同样的 `scf.Kgrid 8 8 8`，两种方式的 k 点集合不同，总能会有小差别，k 点不够密时差别更明显；
- 自旋螺旋（广义 Bloch 定理）计算在 GCenter 下有报错的报告（见"常见问题"），此时改用 `Regular`；
- 3.9 不认识这个关键词，3.9 输入里不要写它。

### 5.4 分子和大超胞

- 孤立分子用 `Cluster`，不需要 `scf.Kgrid`；
- 几百原子以上的超胞可以用 `Band` + `scf.Kgrid 1 1 1`（只算 Γ），但仍要检查 `2 2 2` 是否改变结果；
- 表面板模型真空方向取 1；一维链在两个真空方向取 1。

---

## 6. 混合方案与 SCF 收敛控制

### 6.1 七种混合方案

| `scf.Mixing.Type` | 混合的量 | 表示空间 | 主要用途（据手册第16章） |
|---|---|---|---|
| `Simple` | 密度矩阵 | 实空间 | 小分子、绝缘体；会参考历史自动调节权重 |
| `GR-Pulay` | 密度矩阵 | 实空间 | 用历史稳定收敛（保证残差下降的 Pulay） |
| `RMM-DIIS` | 密度矩阵 | 实空间 | 通用 Pulay 型混合 |
| `Kerker` | 电子密度 | 傅里叶空间 | 抑制金属中的长波电荷振荡（charge sloshing） |
| `RMM-DIISK` | 电子密度 | 傅里叶空间 | **多数体系的首选**，金属也首选 |
| `RMM-DIISV` | Kohn–Sham 势 | 傅里叶空间 | 金属和疑难体系 |
| `RMM-DIISH` | Kohn–Sham 哈密顿量 | 矩阵表示 | DFT+U、约束计算 |

手册的总结是：RMM-DIISK 对很多体系最稳健；DFT+U 和约束 DFT 中哈密顿量依赖占据数矩阵，只混合电子密度可能不足以稳定含 U 项的哈密顿量，所以首选 RMM-DIISH。

### 6.2 相关关键词

| 关键词 | 默认值 | 推荐值 | 说明 |
|---|---|---|---|
| `scf.Init.Mixing.Weight` | 0.3 | 0.01–0.10 | 初始混合比例，太大会发散 |
| `scf.Min.Mixing.Weight` | 0.001 | 0.001–0.01 | Simple/Kerker 混合权重下限 |
| `scf.Max.Mixing.Weight` | 0.4 | 0.1–0.3 | Simple/Kerker 混合权重上限；振荡时调小 |
| `scf.Mixing.History` | 手册未明确给出 | 20–40 | Pulay 型混合参考的历史步数，手册建议 30 左右 |
| `scf.Mixing.StartPulay` | 6 | 10–30 | 从第几步开始 Pulay 型混合，之前用 Simple/Kerker |
| `scf.Mixing.EveryPulay` | 1 | 1 | 每隔几步做一次 Pulay 混合（用于 RMM-DIISK；手册第16章也列出 RMM-DIISV/H），1 即每步 |
| `scf.Kerker.factor` | 自动估计 | 先不写 | Kerker、RMM-DIISK、RMM-DIISV 的长波抑制因子，只有自动值不行时才手调 |
| `scf.criterion` | 1.0e-6 Ha | 1e-7 – 1e-10 | `dUele` 的收敛阈值 |
| `scf.maxIter` | 40 | 100–300 | 最大 SCF 步数，到了就停，**不管是否收敛** |

`scf.criterion` 的选择：

| 用途 | 推荐 |
|---|---|
| 一般单点、结构优化 | 1e-7 – 1e-8 |
| 能带、DOS、磁性能量差 | 1e-8 – 1e-9 |
| 弱相互作用、BSSE、声子有限位移 | 1e-9 – 1e-10 |

### 6.3 三套起步参数

绝缘体、半导体、分子：

```text
scf.Mixing.Type            RMM-DIISK
scf.Init.Mixing.Weight     0.10
scf.Min.Mixing.Weight      0.01
scf.Max.Mixing.Weight      0.30
scf.Mixing.History         20
scf.Mixing.StartPulay      8
scf.Mixing.EveryPulay      1
scf.ElectronicTemperature  300.0
scf.criterion              1.0e-8
scf.maxIter                100
```

金属、磁性金属、表面：

```text
scf.Mixing.Type            RMM-DIISK
scf.Init.Mixing.Weight     0.02
scf.Min.Mixing.Weight      0.01
scf.Max.Mixing.Weight      0.20
scf.Mixing.History         30
scf.Mixing.StartPulay      15
scf.Mixing.EveryPulay      1
scf.ElectronicTemperature  600.0       # 300-1000 K
scf.criterion              1.0e-8
scf.maxIter                200
```

DFT+U、约束 DFT：

```text
scf.Mixing.Type            RMM-DIISH
scf.Init.Mixing.Weight     0.02
scf.Mixing.History         20
scf.Mixing.StartPulay      10
scf.criterion              1.0e-8
scf.maxIter                200
```

### 6.4 调参顺序

手册推荐初学者按下面的顺序一次只改一个参数，否则很难判断到底是哪个改动起了作用：

1. 用 `RMM-DIISK`，不写 `scf.Kerker.factor`（自动估计）；
2. 收敛慢：把 `scf.Mixing.History` 加到 30–40；
3. 发散或大幅振荡：把 `scf.Mixing.StartPulay` 从 10 加到 30，并减小 `scf.Max.Mixing.Weight`；
4. 金属或疑难体系：提高 `scf.ElectronicTemperature`；
5. DFT+U 或约束计算：换成 `RMM-DIISH`；
6. 仍不稳定：再手动设置 `scf.Kerker.factor`（设得太大，长波分量被过度抑制，收敛反而变慢）。

### 6.5 读懂 SCF 过程

`.out` 中有每一步的历史：

```bash
grep -A60 "SCF history at MD" si_scf.out | head -40
```

格式如下（**示意输出，数值仅供参考**）：

```text
   SCF=   1  NormRD=  1.000000000000  Uele= -x.xxxxxxxxxxxx
   SCF=   2  NormRD=  0.2xxxxxxxxxxx  Uele= -x.xxxxxxxxxxxx
   ...
   SCF=  18  NormRD=  0.000000000xxx  Uele= -x.xxxxxxxxxxxx
```

标准输出（`.std`）每步也打印 `NormRD`、`Uele` 和 `dUele`，适合用 `tail -f` 实时观察：

```bash
mpirun -np 8 openmx Si_scf.dat -nt 1 | tee si_scf.std
grep -E "NormRD|dUele" si_scf.std | tail -6
```

判断方法：

- `NormRD` 总体下降几个数量级，`Uele` 的变化越来越小，直到 `dUele` 小于 `scf.criterion`；
- Pulay 型混合中 `NormRD` 偶尔会短暂上升，不要立即判为失败，看几步的趋势；
- **最可靠的收敛标志**：`.md` 文件每个结构块的第二行带有 `scf_conv=1`（收敛）或 `scf_conv=0`（未收敛）：

```bash
grep -o "scf_conv=[01]" si_scf.md
```

- `scf.maxIter` 用完而没收敛时，程序会照常输出能量和力，但这些数是不可信的。只要看到 `scf_conv=0`，就不要使用该结果。

---

## 7. 运行中调参：`System.Name_SCF_keywords`

大体系一次 SCF 可能要几小时。OpenMX 允许在运行中修改以下参数：`scf.maxIter`、`scf.Min.Mixing.Weight`、`scf.Max.Mixing.Weight`、`scf.Kerker.factor`、`scf.Mixing.StartPulay`。

做法：在运行目录（`System.CurrrentDirectory`，默认 `./`）中新建一个名为 `<System.Name>_SCF_keywords` 的文件。例如 `System.Name fe_bcc` 时：

```bash
cat > fe_bcc_SCF_keywords << 'EOF'
scf.maxIter                300
scf.Min.Mixing.Weight      0.005
scf.Max.Mixing.Weight      0.10
scf.Mixing.StartPulay      30
EOF
```

OpenMX 每个 SCF 步都会检查这个文件，读取成功后在标准输出打印：

```text
The keywords for SCF iteration are renewed by ./fe_bcc_SCF_keywords.
```

两个技巧：

- **优雅终止**：把 `scf.maxIter` 写成负数（例如 `scf.maxIter -1`），OpenMX 会停止计算，比直接 `scancel` 更干净；
- 改完要**删掉或改名**这个文件，否则下次用同一个 `System.Name` 运行时它会再次生效。

> [!TIP]
> 调参后若想从当前密度继续，可以结束本次运行，然后在输入中加 `scf.restart on` 重新提交（§8），新的混合参数直接写进输入即可。

---

## 8. 续算与后处理：`scf.restart`

```text
scf.restart     on      # on|off|pp|c2n, default=off ; pp and c2n need OpenMX >= 4.0
```

| 取值 | 作用 | 是否覆盖 `System.Name_rst/` |
|---|---|---|
| `off` | 从原子密度叠加开始 | 写入新的重启文件 |
| `on` | 读取 `System.Name_rst/` 中的哈密顿量和密度作为初值，继续 SCF 直到收敛 | 覆盖 |
| `pp` | 读入后只做**一次对角化**，然后做能带、DOS 等后处理；程序自动设 `scf.maxIter=1`、`MD.Type=nomd` | **不覆盖** |
| `c2n` | 把共线计算的重启文件转换给非共线计算使用（磁各向异性能的二阶变分法） | — |

要点：

- 续算时 `System.Name` 必须与第一次相同，程序据此找到 `System.Name_rst/` 目录；
- 重启文件同时包含密度矩阵混合和傅里叶空间混合所需的信息，所以续算时可以**换混合方案**；
- 结构不能改（改了结构就不是同一个体系）；k 点网格、电子温度、后处理选项可以改；
- 3.9 只有 `on` 和 `off`；`pp`、`c2n` 是 4.0 的功能；`c2n` 的用法见 [第11章](11-spin-noncollinear.md) 和手册第42章。

### 8.1 实战：先 SCF，再用 `pp` 算能带和 DOS

第一步 [`Si_scf.dat`](../examples/05-scf-conv-si/README.md) 做自洽（`System.Name si_scf`）。第二步 `Si_restart_pp.dat` 中与后处理相关的部分如下（结构、基组、SCF 参数与第一步完全相同）：

```text
System.Name                 si_scf      # MUST be identical to the SCF run (reads si_scf_rst/)
scf.restart                 pp          # one diagonalization, restart files are not modified

Band.dispersion             on
<Band.KPath.UnitCell
  5.43  0.00  0.00
  0.00  5.43  0.00
  0.00  0.00  5.43
Band.KPath.UnitCell>
Band.Nkpath                 5
<Band.kpath
  40   0.5  0.5  0.5    0.0  0.0  0.0    L  G
  40   0.0  0.0  0.0    1.0  0.0  0.0    G  X
  20   1.0  0.0  0.0    1.0  0.5  0.0    X  W
  15   1.0  0.5  0.0    0.75 0.75 0.0    W  K
  40   0.75 0.75 0.0    0.0  0.0  0.0    K  G
Band.kpath>

Dos.fileout                 on
Dos.Erange                  -15.0  10.0     # eV, relative to the chemical potential
Dos.Kgrid                   16 16 16        # denser than scf.Kgrid
```

```bash
mpirun -np 8 openmx Si_scf.dat -nt 1 > si_scf.std           # step 1: SCF
ls si_scf_rst/                                               # restart files exist
mpirun -np 8 openmx Si_restart_pp.dat -nt 1 > si_pp.std     # step 2: band + DOS
ls si_scf.Band si_scf.Dos.val si_scf.Dos.vec
grep -E "NormRD" si_pp.std | wc -l                           # expect only one SCF step
```

注意第二步会用同一个 `System.Name` 覆盖 `si_scf.out`，如需保留第一步的输出，先 `cp si_scf.out si_scf_step1.out`。能带和 DOS 的后处理见 [第8章](08-band.md)、[第9章](09-dos-pdos-cohp-optics.md)。

> [!NOTE]
> 结构优化和分子动力学中，程序会自动用前几步的密度外推当前步的初始密度，外推所用的历史步数由 `scf.ExtCharge.History` 控制（默认 2）。手册的测试表明 2 最好，更大反而不稳定，保持默认即可。

---

## 9. 收敛测试实战：Si 的截断能与 k 点

示例目录 [`examples/05-scf-conv-si`](../examples/05-scf-conv-si/README.md) 包含模板和驱动脚本。

### 9.1 测试流程与判据

1. **固定基组和结构**（Standard 基组 `Si7.0-s2p2d1`，a = 5.43 Å）；
2. **截断能测试**：k 点固定为 8×8×8，`scf.energycutoff` 取 120–400 Ry；
3. **k 点测试**：截断能固定为上一步的推荐值，`scf.Kgrid` 取 4–16；
4. **判据**：相对最高精度计算，总能差 |ΔE| < 1 meV/atom；同时检查你真正关心的量（晶格常数、力、带隙、磁矩）；
5. **交叉检查**：用最终选定的组合再与 "更高截断能 + 更密 k 点" 比较一次。

> [!IMPORTANT]
> 收敛判据看的是**能量差**是否稳定，而不只是绝对总能。例如计算吸附能时，应该测试吸附能本身随截断能的变化。另外，OpenMX 的总能随截断能是振荡收敛的，不要用"能量单调下降"来判断。

### 9.2 截断能扫描

模板 `conv_ecut/Si_ecut_template.dat` 中有两处占位符：

```text
System.Name                 si_ecut_@ECUT@
scf.energycutoff            @ECUT@
scf.Kgrid                   8 8 8
```

驱动脚本 `run_conv.sh` 用 sed 替换占位符、依次运行，再从 `.out` 中提取 `Utot.`、`Elapsed.Time.`、网格点数和 `.md` 中的 `scf_conv`：

```bash
cd OpenMX/examples/05-scf-conv-si
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat conv_*/*.dat   # once
MPI_RUN="mpirun -np 8" bash run_conv.sh ecut              # run + summary
ECUTS="150 200 250 300" bash run_conv.sh ecut collect     # re-print a subset
```

汇总表的格式（**示意输出，数值仅供参考**）：

```text
# ecut  Utot(Ha/cell)  dE_vs_last(meV/atom)  time(s)  grid  scf_conv  ok(<1)
   120    -7.8xxxxxxxxx        x.xxx        xx.x     NN,NN,NN    1   no
   150    -7.8xxxxxxxxx       -x.xxx        xx.x     NN,NN,NN    1   no
   ...
   400    -7.8xxxxxxxxx        0.000        xx.x     NN,NN,NN    1  yes
```

读表要点：

- `dE_vs_last` 以列表中最后一个（最高截断能）为参考；
- 低截断能处能量可能偏高也可能偏低，这是积分网格的振荡，属正常；
- 从某个截断能开始，后面所有点都在 ±1 meV/atom 以内，这个截断能就是推荐值；
- `time(s)` 列让你看到代价：网格点数约随 $E_{\text{cut}}^{3/2}$ 增长。

### 9.3 k 点扫描

```text
System.Name                 si_k@K@
scf.energycutoff            220.0       # use the value chosen in the ecut test
scf.Kgrid                   @K@ @K@ @K@
```

```bash
MPI_RUN="mpirun -np 8" bash run_conv.sh kgrid
KGRIDS="6 8 10 12 14 16 20" bash run_conv.sh kgrid       # custom list
```

对 Si 这样的半导体，ΔE 通常在 8×8×8 附近就进入 1 meV/atom 以内；金属需要密得多的网格（§10）。

### 9.4 用仓库脚本做同样的事

本仓库的 `omx_conv_test.py` 可以基于模板批量生成输入、调用 mpirun、解析 `Utot` 并输出 CSV 和 ΔE 曲线，用 `python ../scripts/omx_conv_test.py -h` 查看用法；`omx_parse_out.py` 可以把 `.out` 中的能量、SCF 历史、Mulliken 电荷导出为 JSON/CSV。

### 9.5 把结果写进记录

建议在项目目录中保存一份 `convergence.md`，至少记录：OpenMX 版本、数据库版本、基组、截断能测试表、k 点测试表、最终选择和理由。审稿人问到时可以直接回答。

---

## 10. 金属体系：bcc Fe 的 SCF 与 k 点

示例目录 [`examples/05-metal-scf-fe`](../examples/05-metal-scf-fe/README.md) 中的 `Fe_bcc_metal_scf.dat`：

```text
#
# bcc Fe, ferromagnetic, collinear GGA-PBE: metallic SCF settings
# OpenMX 4.0 (+patch 4.0.1), database DFT_DATA19
#
System.CurrrentDirectory    ./          # official spelling has three r's
System.Name                 fe_bcc
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            0

Species.Number              1
<Definition.of.Atomic.Species
  Fe   Fe5.5H-s3p2d1   Fe_PBE19H    # Standard basis; Fe_PBE19H has 16 valence electrons
Definition.of.Atomic.Species>

Atoms.Number                1
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1   Fe   0.0  0.0  0.0   9.0  7.0      # up + down = 16; initial moment 2 muB
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  -1.4350   1.4350   1.4350             # bcc primitive vectors, a = 2.87 Ang
   1.4350  -1.4350   1.4350
   1.4350   1.4350  -1.4350
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        on
scf.ElectronicTemperature   600.0       # K; metals: 300-1000 K
scf.energycutoff            300.0       # Ry; hard VPS with 3s/3p semicore states
scf.EigenvalueSolver        Band
scf.Kgrid                   16 16 16
scf.Generation.Kpoint       GCenter     # OpenMX 4.0 default; remove for OpenMX 3.9
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.02
scf.Min.Mixing.Weight       0.01
scf.Max.Mixing.Weight       0.20
scf.Mixing.History          30
scf.Mixing.StartPulay       15
scf.Mixing.EveryPulay       1
scf.maxIter                 200
scf.criterion               1.0e-8      # Hartree

MD.Type                     Nomd
```

参数解释：

| 关键词 | 取值 | 理由 |
|---|---|---|
| VPS / PAO | `Fe_PBE19H` + `Fe5.5H-s3p2d1` | 推荐表 Standard；致密 bcc 用 5.5 bohr 小半径避免过完备（第4章 §10） |
| 初始占据 | `9.0 7.0` | 16 个价电子，初始磁矩 2 μB，接近铁磁基态，SCF 更快 |
| `scf.energycutoff` | 300 Ry | Hard 赝势含 3s/3p 半芯态，比 Si 需要更细的网格 |
| `scf.ElectronicTemperature` | 600 K | 金属费米面附近的占据平滑化，SCF 更稳定 |
| `scf.Kgrid` | 16×16×16 | bcc 原胞晶格矢量长 2.49 Å，16 × 2.49 ≈ 40 Å；还需测试 |
| 混合参数 | 小权重、长历史、晚开始 Pulay | 金属/磁性体系的稳妥组合（§6.3） |

运行和检查：

```bash
cd OpenMX/examples/05-metal-scf-fe
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
mpirun -np 8 openmx Fe_bcc_metal_scf.dat -nt 1 > fe_bcc.std
grep -o "scf_conv=[01]" fe_bcc.md                     # 1 = converged
grep -E "Utot\.|Chemical potential" fe_bcc.out
grep "Total spin moment" fe_bcc.out                   # spin moment per cell in muB
grep "Total spin S" fe_bcc.out                        # S; moment = 2S muB
```

bcc Fe 的实验磁矩约 2.2 μB/atom，PBE 计算一般给出 2.2–2.3 μB 左右，作为量级参考。

### 10.1 电子温度与 k 点的耦合测试

脚本 `run_scan.sh` 可以扫描电子温度（`temp`）或 k 点（`kgrid`），并汇总 Utot、磁矩和收敛标志：

```bash
TEMPS="300 600 1000" MPI_RUN="mpirun -np 8" bash run_scan.sh temp
KGRIDS="10 12 14 16 18 20 24" MPI_RUN="mpirun -np 8" bash run_scan.sh kgrid
```

观察要点：

- 电子温度越低，达到同样精度所需的 k 点越多；
- 磁矩对电子温度比总能更敏感，报告磁矩时注明所用温度；
- 能量差（例如 FM 与 AFM、不同晶格常数之间）应在**同一温度、同一 k 网格**下比较。

### 10.2 金属 SCF 的排障清单

1. 确认初始磁矩合理（铁磁 Fe 不要写成 `8.0 8.0`，否则可能收敛到非磁态）；
2. `scf.Init.Mixing.Weight` 降到 0.01–0.02，`scf.Mixing.StartPulay` 推后到 20–30；
3. 电子温度提高到 1000 K，收敛后再用 `scf.restart on` 降回目标温度；
4. 仍振荡时减小 `scf.Max.Mixing.Weight`，最后才手动设 `scf.Kerker.factor`；
5. 运行中可以随时用 `fe_bcc_SCF_keywords` 文件（示例目录中有 `.example` 模板）调整，不必重新开始。

---

## 参数选择与收敛建议

| 参数 | 起步推荐 | 何时调整 | 验证 |
|---|---|---|---|
| `scf.XcType` | `GGA-PBE` + `*_PBE19` | 与 LDA 文献对比时用 `LDA`/`LSDA-*` + `*_CA19` | VPS 的 `xc.type` |
| `scf.energycutoff` | 200–250 Ry（显式写出） | 半芯 Hard 赝势、弱相互作用、微小磁能差：300–400 Ry | §9.2，1 meV/atom |
| `scf.Ngrid` | 不写 | E–V、应变扫描、超胞：必须写 | 网格一致 |
| `scf.Kgrid` | n_i × a_i ≈ 30–40 Å（金属 ≥ 50 Å） | 金属、小胞、DOS | §9.3 |
| `scf.Generation.Kpoint` | `GCenter`（4.0 默认） | 复现 3.9、自旋螺旋：`Regular` | — |
| `scf.ElectronicTemperature` | 300 K | 金属 300–1000 K；MAE 等微小能量差要低 | 磁矩、能量差 |
| `scf.Mixing.Type` | `RMM-DIISK` | +U/约束：`RMM-DIISH` | `NormRD` 曲线 |
| `scf.Mixing.History` / `StartPulay` | 20–30 / 8–15 | 振荡：StartPulay 20–30 | — |
| `scf.criterion` | 1e-8 Ha | 弱相互作用 1e-9–1e-10 | `scf_conv=1` |
| `scf.maxIter` | 100–200 | 大体系、金属更多 | `scf_conv=1` |

---

## 常见问题与报错

**1. SCF 不收敛：`NormRD` 长时间振荡或上升，`.md` 中出现 `scf_conv=0`**

- 原因：混合太激进；金属电子温度太低；初始磁矩不合理；+U 体系只混合了密度；基组过完备。
- 解决：按 §6.4 的顺序调整；金属提高电子温度；+U 改 `RMM-DIISH`；运行中用 `System.Name_SCF_keywords` 降低混合权重；排查过完备见 [第4章 §10](04-pao-vps.md)。

**2. 程序"正常结束"，但结果不对**

- 原因：`scf.maxIter` 用完，SCF 未收敛，程序仍照常输出能量和力。
- 解决：每次都检查 `grep -o "scf_conv=[01]" *.md`；把 `scf.maxIter` 设得足够大。

**3. E–V 曲线或应变–能量曲线出现锯齿、跳变**

- 原因：晶格变化导致实空间网格点数跳变。
- 解决：用 `scf.Ngrid` 固定网格（§4.3），或用 `MD.Type EvsLC` 并固定 `scf.Ngrid`（第6章）。

**4. 同一个输入在 4.0 和 3.9 下总能不同**

- 原因：4.0 默认 `scf.Generation.Kpoint GCenter`，3.9 只有 Regular 网格；另外两者对 `scf.energycutoff` 默认值的说明也不同。
- 解决：显式写 `scf.Generation.Kpoint Regular` 和 `scf.energycutoff`，再比较。

**5. 自旋螺旋（广义 Bloch）计算报 `Fatal error in internal_Comm_size: Invalid communicator`**

- 原因：OpenMX 论坛报告该功能与 4.0 默认的 GCenter k 点生成方式组合时会出现此 MPI 错误。
- 解决：设 `scf.Generation.Kpoint Regular`。

**6. `scf.restart on` 后像是从头开始，或报找不到重启文件**

- 原因：第二次运行的 `System.Name` 与第一次不同，或 `System.Name_rst/` 目录被删除、不在 `System.CurrrentDirectory` 下。
- 解决：保持 `System.Name` 一致；确认 `ls <System.Name>_rst/` 有文件。

**7. 改了混合参数文件却没有生效**

- 原因：文件名不是 `<System.Name>_SCF_keywords`，或不在运行目录；或者改的关键词不在可在线修改的列表中。
- 解决：看标准输出是否出现 `The keywords for SCF iteration are renewed by ...`。可在线修改的只有 §7 列出的几个参数。

**8. 大体系或开 OpenMP 后出现 `Segmentation fault (signal 11)`**

- 原因：栈空间不足或内存不够。
- 解决：作业脚本中加 `ulimit -s unlimited` 和 `export OMP_STACKSIZE=512M`，减少 `-nt` 线程数或增加节点，详见 第18章。

---

## 练习

1. 对 Si 分别用 `scf.Generation.Kpoint GCenter` 和 `Regular` 做 4–12 的 k 点扫描，比较两条收敛曲线。哪种在偶数网格下收敛更快？
2. 用 `run_conv.sh ecut` 的结果画 ΔE–E_cut 曲线，并同时画出计算时间，选择性价比最高的截断能。
3. 把 Fe 例子的初始占据改成 `8.0 8.0`，观察 SCF 收敛到什么磁态；再改回 `9.0 7.0`。
4. 把 Fe 例子的 `scf.Mixing.Type` 换成 `Simple`，记录 SCF 步数，并与 RMM-DIISK 比较。
5. 在 Fe 计算运行期间创建 `fe_bcc_SCF_keywords`，写入 `scf.maxIter -1`，确认程序会停止，再用 `scf.restart on` 续算。

---

## 参考资料

**官方文档**

- [OpenMX 4.0 手册 8.2 Keywords](https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html) — 全部 SCF 关键词与默认值
- [手册 第14章 Cutoff energy / grid](https://www.openmx-square.org/openmx_man4.0/s14_cutoff-grid.html) — 截断能收敛、E–V 曲线跳变与 `scf.Ngrid`、`scf.fixed.grid`
- [手册 第15章 k-point sampling](https://www.openmx-square.org/openmx_man4.0/s15_k-sampling.html) — GCenter 与 Regular
- [手册 第16章 SCF convergence](https://www.openmx-square.org/openmx_man4.0/s16_scf.html) — 七种混合方案、调参顺序、Kerker 因子自动估计、在线调参
- [手册 第17章 Restart](https://www.openmx-square.org/openmx_man4.0/s17_restart.html) — `scf.restart on|pp|c2n`、`scf.ExtCharge.History`
- [手册 第11章 Functionals](https://www.openmx-square.org/openmx_man4.0/s11_functionals.html) — `scf.XcType` 与 `scf.SpinPolarization`
- [OpenMX 3.9 手册](https://www.openmx-square.org/openmx_man3.9/index.html) — 与 4.0 对照（`scf.energycutoff` 默认值、无 GCenter）
- [What's new](https://www.openmx-square.org/whatsnew/whatsnew.html) — 4.0 的 GCenter k 点、`scf.restart pp`
- [Practical guide to OpenMX calculations（PDF）](https://t-ozaki.issp.u-tokyo.ac.jp/lectures-2023/OpenMX_Practical_Guide.pdf) — 截断能与混合方案的经验建议
- [OpenMX Forum](https://www.openmx-square.org/forum/patio.cgi) — 网格、GCenter 等问题的讨论

**相关 GitHub 仓库（仅链接）**

- [pranabdas/openmx](https://github.com/pranabdas/openmx) — 英文教程，含截断能和 k 点收敛测试示例（无许可证，仅参考）
- [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx) — 按 $n_i\cdot d_i\ge 35$ Å 自动设置 `scf.Kgrid`
- [samlipton/pyfesm](https://github.com/samlipton/pyfesm) — 轻量 OpenMX 输出解析器，可借鉴其解析思路
- [hirtatsu/Calculations](https://github.com/hirtatsu/Calculations) — 日文实操笔记，含运行与日志检查流程
- [ASE OpenMX calculator](https://gitlab.com/ase/ase/-/tree/master/ase/calculators/openmx) — 用 Python 批量生成和运行 OpenMX 计算（见 Common 第4章）

**关键文献**

- T. Ozaki, H. Kino, *Phys. Rev. B* 72, 045121 (2005)（OpenMX 的网格积分与投影展开方法）. DOI: [10.1103/PhysRevB.72.045121](https://doi.org/10.1103/PhysRevB.72.045121)
- G. Kresse, J. Furthmüller, *Phys. Rev. B* 54, 11169 (1996)（RMM-DIIS 与 Kerker 度量的 Pulay 混合）. DOI: [10.1103/PhysRevB.54.11169](https://doi.org/10.1103/PhysRevB.54.11169)
- G. P. Kerker, *Phys. Rev. B* 23, 3082 (1981)（Kerker 混合）. DOI: [10.1103/PhysRevB.23.3082](https://doi.org/10.1103/PhysRevB.23.3082)
- P. Pulay, *Chem. Phys. Lett.* 73, 393 (1980)（DIIS）. DOI: [10.1016/0009-2614(80)80396-4](https://doi.org/10.1016/0009-2614(80)80396-4)
- D. R. Bowler, M. J. Gillan, *Chem. Phys. Lett.* 325, 473 (2000)（GR-Pulay）. DOI: [10.1016/S0009-2614(00)00750-8](https://doi.org/10.1016/S0009-2614(00)00750-8)
- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* 13, 5188 (1976)（均匀 k 点网格）. DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)
- N. D. Mermin, *Phys. Rev.* 137, A1441 (1965)（有限温度密度泛函）. DOI: [10.1103/PhysRev.137.A1441](https://doi.org/10.1103/PhysRev.137.A1441)

---

[← 上一章](04-pao-vps.md) | [返回目录](../README.md) | [下一章 →](06-optimization-eos-neb.md)
