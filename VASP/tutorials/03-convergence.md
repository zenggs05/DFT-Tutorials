# 第3章 收敛测试与计算精度控制（Convergence Tests & Accuracy Control in VASP: ENCUT, k-points, SIGMA）

本章系统讲解 VASP 计算的收敛测试（convergence test）与精度控制：如何扫描 ENCUT 截断能和 k 点网格（KPOINTS/KSPACING）、如何为半导体和金属选择 ISMEAR/SIGMA 并检查熵项 `T*S`，以及 PREC、EDIFF、LREAL、LASPH、ADDGRID 等精度标签的取值依据。除能量外，本章还强调力与应力（Pulay stress）的收敛，并提供 bash + Python 的批量测试脚本（Si 与金属 Al 两个算例），可直接改造用于你自己的体系。

> [!NOTE]
> 本章的全部输入与脚本在示例目录 [03-convergence-si-al](../examples/03-convergence-si-al/README.md) 中。POTCAR 需自备（见 [第2章](02-input-output-files.md)）。文中带“示意输出”标注的数据仅用于说明格式和趋势，不是真实计算结果。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| ENCUT | 从约 1.0× 到 1.6×（及以上）max(ENMAX) 扫描，以 ΔE/atom < 1 meV 为判据；变胞计算 ≥ 1.3×ENMAX |
| k 点 | 按网格 N×N×N 或 KSPACING（0.15–0.25 Å⁻¹ 量级）扫描；金属更密；表面 z 方向取 1 |
| 展宽 | 半导体 ISMEAR=0、SIGMA=0.01–0.05；金属优化 ISMEAR=1/2、SIGMA=0.1–0.2；DOS 与精确总能 ISMEAR=−5；检查 \|`T*S`\| < 1 meV/atom |
| 精度标签 | PREC=Accurate、EDIFF（1E-5→1E-6，声子/弹性 1E-8）、LREAL、LASPH、ADDGRID |
| 力与应力 | 不只看能量；力常数类计算需要更严格的设置 |
| 自动化 | 批量建目录、本地运行或提交 Slurm、汇总、作图 |

适用场景：**每一个新体系、新赝势组合或新泛函都应该做收敛测试**。测试结果决定了后续结构优化（[第4章](04-relax-eos.md)）、能带（[第5章](05-band-dos.md)）、声子（[第12章](12-phonon.md)）等所有计算的基本参数。

---

## 2. 原理速览

### 2.1 平面波截断能

VASP 用平面波展开 Kohn–Sham 轨道，只保留动能低于截断能的平面波：

$$
\frac{\hbar^2}{2m_e}|\mathbf{G}+\mathbf{k}|^2 < E_{\text{cut}} = \texttt{ENCUT}
$$

基组越大，对给定结构的总能越低（变分原理），并逐渐收敛。要点：

- **能量差比总能收敛得快**：结构相近的两个构型（如优化前后），截断误差大部分相互抵消；差别很大的构型（体相 vs 孤立原子的结合能）收敛则慢得多（官方 Energy cutoff and FFT meshes 页面）；
- POTCAR 中的 ENMAX 是“推荐”截断能，官方说明用默认最小值时结合能误差一般小于约 10 meV；**低于 ENMAX 的值不应使用**；
- **ENCUT 必须在 INCAR 中显式写出**，并在所有需要比较能量的计算中保持一致（默认值随元素组合变化，例如化合物与单质的默认 ENCUT 不同）。

### 2.2 Pulay 应力

VASP 的变胞优化在**固定平面波数目**（固定基组）下进行。体积变化时，固定数目的平面波对应的有效截断能随之改变，基组不完备会产生一个虚假的应力——Pulay stress。它几乎总是使晶胞**偏小**（表现为负的虚假压强）。对策：

1. 涉及应力的计算（ISIF ≥ 3、状态方程、弹性常数）使用 **ENCUT ≥ 1.3×max(ENMAX)**；
2. 变胞优化结束后用 CONTCAR 重新开始再优化一两次，直到体积不再变化（[第4章](04-relax-eos.md)）；
3. 在 ENCUT 测试中同时记录压强，观察它随 ENCUT 的收敛（本章第 4 节）。

### 2.3 k 点积分

布里渊区积分被离散为加权求和 $\frac{1}{\Omega_{\text{BZ}}}\int d^3k \to \sum_{\mathbf{k}} w_{\mathbf{k}}$。对有带隙的体系，被积函数光滑，误差随网格加密迅速下降；对金属，费米面处占据数从 1 跳到 0，直接求和收敛很慢，因此需要**展宽（smearing）**或**四面体方法**。k 点网格的密度与晶胞大小成反比：晶胞越大，所需的 $N_i$ 越小。

### 2.4 展宽与自由能

引入展宽后，VASP 最小化的是广义自由能

$$
F = E - \sum_{n\mathbf{k}} w_{\mathbf{k}}\,\sigma\,S[f_\sigma(\epsilon_{n\mathbf{k}})] = E - TS
$$

官方给出的 σ→0 外推为 $E_0 = \frac{1}{2}(F+E)$，即 OUTCAR 中的 `energy(sigma->0)`。需要记住：

- **力和应力是 $F$ 的导数**，不是 $E_0$ 的导数，因此力和应力也要对 SIGMA 收敛；
- Methfessel–Paxton（ISMEAR ≥ 1）能给出很准的金属总能，但占据数非单调，**不能用于半导体和绝缘体**（官方警告：声子频率误差可超过 20%）；
- 带 Blöchl 修正的四面体方法（ISMEAR = −5）给出精确的总能与 DOS，但对占据数不是变分的，**金属中的力和应力可能偏差 5–10%**，所以不用于金属的结构优化。

### 2.5 “收敛”到什么程度

| 物理量 | 常用收敛目标 | 备注 |
|---|---|---|
| 能量差（每原子） | < 1 meV/atom | 形成能、吸附能等更关心差值 |
| 力 | 变化 < 0.005–0.01 eV/Å | 声子/力常数要求更高 |
| 应力 / 压强 | 变化 < 1 kB 左右 | 晶格常数、弹性常数 |
| 熵项 \|`T*S`\| | < 1 meV/atom | 金属展宽检查 |
| 力的漂移 `total drift` | 远小于关心的力（官方：一般应小于你关心的力，常见量级 0.1 eV/Å 以下） | 反映 FFT 网格是否足够 |

> [!IMPORTANT]
> “收敛”永远是针对**你关心的物理量**而言的。只测试单个体系的总能，并不能保证吸附能、能垒或声子频率收敛。条件允许时，直接对目标量（能量差、力、频率、带隙）做收敛测试。

---

## 3. 测试策略总览

```text
1. 固定泛函、POTCAR、PREC=Accurate、EDIFF 比生产计算严一个量级
2. 中等密度 k 点 (如 Si 用 8x8x8) 下扫描 ENCUT       -> 选 ENCUT
3. 用选定的 ENCUT 扫描 k 点网格                        -> 选 k 网格
4. 金属: 在若干 k 网格上扫描 ISMEAR/SIGMA              -> 选 SIGMA 与 k 网格
5. 用最终参数复核: 力、应力、(必要时) 目标物理量
```

几条原则：

- **一次只改一个参数**，其余全部固定（包括 ISTART=0，避免读入旧 WAVECAR 引入差异）；
- 测试体系尽量接近生产体系：同样的元素、POTCAR、泛函、自旋设置；大超胞的 k 点可由原胞结果按晶胞尺寸换算；
- 测试时 EDIFF 要足够小（如 1E-7），否则 SCF 噪声会掩盖收敛趋势；
- 选取参数时留余量：收敛点再往上取一档，生产计算中各种体系变化不会立即让你掉出收敛区。

---

## 4. ENCUT 收敛测试（Si）

### 4.1 输入

`encut_test/INCAR.template`（`__ENCUT__` 由脚本替换）：

```text
SYSTEM = Si ENCUT convergence test

ISTART = 0
ICHARG = 2

PREC   = Accurate
ENCUT  = __ENCUT__     # replaced by run_encut.sh
LREAL  = .FALSE.

ALGO   = Normal
EDIFF  = 1E-7          # tighter than production: SCF noise << 1 meV/atom
NELM   = 120

ISMEAR = 0
SIGMA  = 0.05

IBRION = -1            # single point
NSW    = 0
ISIF   = 2             # forces and stress tensor are computed

LWAVE  = .FALSE.
LCHARG = .FALSE.
```

| 关键词 | 取值 | 原因 |
|---|---|---|
| `PREC` | Accurate | 与生产计算一致，避免 FFT 网格引入额外误差 |
| `EDIFF` | 1E-7 | 让电子收敛误差远小于 1 meV/atom 的判据 |
| `LREAL` | .FALSE. | 消除实空间投影误差，使趋势只反映 ENCUT |
| `ISIF` | 2 | 单点计算时也输出应力张量，用于观察 Pulay 应力 |
| `LWAVE`/`LCHARG` | .FALSE. | 每个目录只保留小文件 |

KPOINTS 固定为 Γ 中心 8×8×8。POSCAR 为理想金刚石 Si（a = 5.43 Å）；另有 `POSCAR_displaced`，把第二个 Si 原子的分数坐标从 0.25 改为 0.27（沿 $\mathbf{a}_1$ 约 0.08 Å），使原子受力不为零，用于力的收敛测试。

Si 的 ENMAX = 245.3 eV，脚本默认扫描 250–500 eV：

| ENCUT (eV) | 250 | 275 | 300 | 325 | 350 | 375 | 400 | 450 | 500 |
|---|---|---|---|---|---|---|---|---|---|
| ENCUT / ENMAX | 1.02 | 1.12 | 1.22 | 1.32 | 1.43 | 1.53 | 1.63 | 1.83 | 2.04 |

最高的 450、500 eV 作为“参考值”，用于计算其他点的偏差。

### 4.2 运行

```bash
cd VASP/examples/03-convergence-si-al/encut_test
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

# 本机/交互节点直接运行（默认 mpirun -np 4 vasp_std）
NP=4 bash run_encut.sh                    # 理想结构
NP=4 bash run_encut.sh POSCAR_displaced   # 有位移的结构 -> 力的收敛

# 或者在集群上为每个 ENCUT 提交一个作业，完成后再汇总
MODE=slurm bash run_encut.sh
MODE=collect bash run_encut.sh            # 作业全部结束后执行
```

`MODE=slurm` 时脚本会把 [示例 01 的 job_cpu.slurm](../examples/01-setup-compile/job_cpu.slurm) 复制到每个子目录并 `sbatch`（可用 `JOB_TEMPLATE=/path/to/your.slurm` 指定自己的模板）。

脚本的核心逻辑非常简单（完整文件见 [run_encut.sh](../examples/03-convergence-si-al/encut_test/run_encut.sh)，公共函数见 [common.sh](../examples/03-convergence-si-al/common.sh)）：

```bash
for e in $ENCUTS; do
    d=${TAG}_ENCUT_$e
    mkdir -p "$d"
    sed "s/__ENCUT__/$e/" INCAR.template > "$d/INCAR"
    cp "$POSCAR_FILE" "$d/POSCAR"
    cp KPOINTS POTCAR "$d/"
    run_or_submit "$d"            # local: (cd $d && mpirun ... vasp_std) ; slurm: sbatch
done
```

汇总时，`outcar_row` 函数从每个 OUTCAR 中提取：`energy(sigma->0)`、`TOTEN`、`EENTRO`、`external pressure`、最后一步的最大原子受力、`NKPTS` 和耗时，写入 `encut_POSCAR.dat`。

### 4.3 分析

```bash
python ../plot_conv.py encut_POSCAR.dat
python ../plot_conv.py encut_POSCAR_displaced.dat --thr-f 5
```

[plot_conv.py](../examples/03-convergence-si-al/plot_conv.py) 以最高 ENCUT 为参考，打印每个点的能量偏差（meV/atom）、压强偏差（kB）和最大力偏差（meV/Å），给出“该点及之后所有点都满足阈值”的最小 ENCUT，并生成 PNG 图。输出格式如下（**示意输出，由合成数据生成，仅演示格式，不是 Si 的真实结果**）：

```text
reference: ENCUT = 500
     ENCUT  nkpts   E0/atom (eV)  dE (meV/at)    P (kB)  dP (kB)  Fmax (eV/A)  dF (meV/A)
       250     29      -5.347893       71.838    128.68   131.17      0.56435       14.30
       ...
       400     29      -5.418799        0.933      1.70     4.19      0.55036        0.31
       450     29      -5.419551        0.180     -1.45     1.04      0.55012        0.07
       500     29      -5.419732        0.000     -2.49     0.00      0.55005        0.00
    energy converged within 1 meV/atom: ENCUT = 400
  pressure converged within 1 kB: ENCUT = 500
max. force converged within 5 meV/A: ENCUT = 300
```

如何解读真实结果：

1. **能量**：ΔE/atom 随 ENCUT 增大而单调减小（若出现明显非单调，检查 PREC、LREAL、EDIFF 是否一致）。取 ΔE < 1 meV/atom 的最小值，再结合 2.2 节的 1.3×ENMAX 规则；
2. **压强**：常见情况是压强比能量收敛得慢——能量已经收敛时压强仍在变化，这正是 Pulay 应力。做变胞优化、状态方程时以压强的收敛为准；
3. **力**（`POSCAR_displaced`）：力的偏差通常比能量更早满足阈值，但声子等需要高精度力的计算要单独确认；
4. 最终写进生产 INCAR 的 ENCUT 一般取收敛点并向上取整（例如 Si 常用 400–520 eV 量级做高精度变胞计算，具体以你的测试为准）。

> [!TIP]
> 也可以用 Python 直接解析各目录（ASE 会把 `energy(sigma->0)` 作为势能返回）：
>
> ```python
> import glob
> from ase.io import read
> for d in sorted(glob.glob("POSCAR_ENCUT_*"), key=lambda s: int(s.split("_")[-1])):
>     at = read(f"{d}/OUTCAR")
>     print(d.split("_")[-1], at.get_potential_energy() / len(at))
> ```

### 4.4 多元素体系与特殊情况

- 多元素体系的 ENMAX 取所有元素中的**最大值**（例如含 O 的氧化物，O 的 ENMAX = 400 eV，1.3 倍就是 520 eV）；
- 含 `_h`、`_GW` 或第一行元素（C、N、O、F）的体系截断能偏高，测试范围相应上移；
- 杂化泛函、GW 等方法的计算量对 ENCUT 非常敏感，但截断能仍需与 PBE 测试结果一致或更高；
- 机器学习力场（MLFF）训练与 NPT 模拟也需要较高的 ENCUT（[第15章](15-md-mlff.md)）。

---

## 5. k 点收敛测试

### 5.1 网格扫描（N×N×N）

使用上一节确定的 ENCUT（示例 INCAR 中为 320 eV），扫描 Γ 中心网格：

```bash
cd ../kpoint_test
cp ../encut_test/POTCAR .
bash run_kpts.sh                     # MESHES="2 4 6 8 10 12 14 16"
python ../plot_conv.py kpts_mesh.dat
```

脚本为每个 N 写入：

```text
Gamma-centred 8x8x8
0
Gamma
  8 8 8
  0 0 0
```

金刚石 Si 原胞（Oh 点群）中，Γ 中心 N×N×N 网格约化后的不可约 k 点数（按 48 个点群操作计算，以 OUTCAR 中的 `NKPTS` 为准）：

| N | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|
| NKPTS | 3 | 8 | 16 | 29 | 47 | 72 | 104 | 145 |

计算量大致与 NKPTS 成正比，因此从 8×8×8 到 16×16×16，成本约增加 5 倍。

### 5.2 KSPACING 扫描

不写 KPOINTS 文件、在 INCAR 中设置 KSPACING 也是常用的方法，它对不同大小的晶胞自动给出“密度一致”的网格：

```bash
SCAN=kspacing bash run_kpts.sh                 # SPACINGS="0.5 0.4 0.3 0.25 0.2 0.15"
python ../plot_conv.py kpts_kspacing.dat --ref min
```

按 $N_i = \max(1, \lceil 2\pi|\mathbf{b}_i| / \text{KSPACING} \rceil)$ 计算得到的网格：

| KSPACING (Å⁻¹) | 0.5 | 0.4 | 0.3 | 0.25 | 0.2 | 0.15 | 0.12 |
|---|---|---|---|---|---|---|---|
| Si（a = 5.43 Å，原胞） | 5³ | 6³ | 7³ | 9³ | 11³ | 14³ | 17³ |
| Al（a = 4.05 Å，原胞） | 6³ | 7³ | 9³ | 11³ | 14³ | 18³ | 23³ |
| Cu（a = 3.61 Å，原胞） | 7³ | 8³ | 11³ | 13³ | 16³ | 21³ | 26³ |

> [!NOTE]
> KSPACING 由上取整决定网格，晶格常数稍有变化就可能让某个方向跳一格（例如 Si 在 KSPACING = 0.2 时 $2\pi|\mathbf{b}|/0.2 \approx 10.02$，刚好越过 10）。状态方程、弹性常数等需要比较不同晶格的计算，最好固定显式网格。

### 5.3 经验与特殊体系

| 体系 | 经验 |
|---|---|
| 宽带隙绝缘体、半导体 | 收敛快；KSPACING 0.25–0.3 Å⁻¹ 常已足够（需测试） |
| 金属 | 收敛慢，需要更密的网格（KSPACING 0.15 Å⁻¹ 或更小），并与 SIGMA 联合测试（第 6 节） |
| 表面板层 / 二维材料 | 面内按体相密度换算，**真空方向取 1**（如 `12 12 1`），Γ 中心 |
| 超胞 | $N_i$ 与晶胞长度成反比：原胞用 12×12×12，则 2×2×2 超胞用 6×6×6 即密度相同 |
| 分子、团簇（大盒子） | 只用 Γ 点（`Gamma 1 1 1` + `vasp_gam`） |
| 六方、fcc | Γ 中心（偶数 MP 网格会破坏对称性） |
| DOS、光学、费米面 | 比总能收敛需要更密的网格，单独测试 |

官方 KPOINTS 页面给出的长度参数经验：$R_k = N_i|\mathbf{a}_i|$ 约 10 Å（宽带隙绝缘体）到约 100 Å（d 金属）之间，与 $R_k = 2\pi/\text{KSPACING}$ 的关系对应 KSPACING 约 0.06–0.6 Å⁻¹。这个范围很宽，所以**测试不可省略**。

> [!TIP]
> 选网格时顺便考虑并行：KPAR 最好能整除不可约 k 点数（[第1章](01-setup-compile.md) 第 13 节）。例如 29 个 k 点是素数，KPAR 只能取 1 或 29。

---

## 6. ISMEAR / SIGMA 的选择与检查（金属 Al）

### 6.1 选择表（依据官方 Smearing technique 页面）

| 体系 / 任务 | ISMEAR | SIGMA (eV) | 说明 |
|---|---|---|---|
| 不清楚是否有带隙、高通量 | 0 | 0.03–0.1 | 最安全的默认；VASP ≥ 6.4 可加 `EFERMI = MIDGAP` |
| 半导体、绝缘体 | 0 | 0.01–0.05 | 或 ISMEAR = −5（≥ 4 个 k 点）；**不要用 ISMEAR > 0** |
| 金属结构优化、受力、声子 | 1 或 2 | 0.1–0.2（官方：默认 0.2 对金属常是合理选择） | 保证 \|`T*S`\| < 1 meV/atom |
| 精确总能、DOS（不做优化） | −5 | — | 需 Γ 中心网格、≥ 4 个 k 点；金属的力不可靠 |
| 需要真实电子温度 | −1 | $k_BT$ | Fermi–Dirac |
| 分子、Γ 点计算、表面 Γ 计算 | 0 | 0.01–0.05 | 四面体方法不可用 |

> [!WARNING]
> VASP 的默认值是 `ISMEAR = 1`、`SIGMA = 0.2`——这是为金属准备的。**半导体、绝缘体、分子如果忘了写 ISMEAR，就会用上 Methfessel–Paxton 展宽**，可能得到不合理的占据数、错误的能量和声子频率。所以每个 INCAR 都应显式写出 ISMEAR 和 SIGMA。

### 6.2 检查熵项

```bash
grep "EENTRO" OUTCAR | tail -1                     # entropy T*S (eV, 整个晶胞)
n=$(sed -n 7p POSCAR | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}')
grep "EENTRO" OUTCAR | tail -1 | awk -v n=$n '{printf "T*S per atom = %.3f meV\n", 1000*$NF/n}'
grep -E "free +energy +TOTEN|energy\(sigma->0\)" OUTCAR | tail -2
```

判据：每原子 `|T*S|` < 1 meV（官方对 Methfessel–Paxton 给出的标准）。若超出，减小 SIGMA 并相应加密 k 点；对 Gaussian 展宽，比较能量时用 `energy(sigma->0)`。

### 6.3 Al 的展宽 + k 点联合测试

fcc Al（a = 4.05 Å，1 原子原胞，ENCUT = 320 eV）是经典的简单金属。[run_sigma.sh](../examples/03-convergence-si-al/sigma_check_Al/run_sigma.sh) 默认扫描：

- ISMEAR ∈ {1, 0}，SIGMA ∈ {0.05, 0.1, 0.2, 0.3}，k 网格 N ∈ {12, 16, 20}；
- 外加一个参考计算：ISMEAR = −5、24×24×24 Γ 网格。

```bash
cd ../sigma_check_Al
cat $VASP_PP_PATH/potpaw_PBE/Al/POTCAR > POTCAR
bash run_sigma.sh               # 25 个 1 原子小计算，单机几分钟到十几分钟
```

输出 `sigma_Al.dat` 的每一行包含 ISMEAR、SIGMA、N、E0、TOTEN、EENTRO、每原子 `T*S`（meV）和相对四面体参考的 ΔE0（meV/atom），脚本最后逐行标出同时满足 `|T*S|` < 1 meV/atom 与 |ΔE0| < 1 meV/atom 的组合。解读要点（定性趋势）：

1. 同一 SIGMA 下，k 网格加密使 E0 趋于稳定；SIGMA 越小，需要的 k 网格越密；
2. 同一 SIGMA 下，Methfessel–Paxton（ISMEAR=1）的 `|T*S|` 通常明显小于 Gaussian（ISMEAR=0），这就是金属优化推荐 MP 的原因；
3. 选择“满足判据且 SIGMA 尽量大”的组合：SIGMA 大则 SCF 稳定、所需 k 点少；
4. 用最终选定的 SIGMA，再回到第 5 节的 k 点测试复核一次。

> [!NOTE]
> 较新版本的 VASP 提供 `ISMEAR = -3` 与 `SMEARINGS` 标签，可以在一次运行中依次计算多组展宽参数（借用离子循环实现：NSW 设为展宽组数，因此不能与结构优化同时使用）。使用前请在官方 SMEARINGS 页面确认你的版本是否支持；本章脚本采用“一组参数一个目录”的方式，适用于所有版本。

---

## 7. 精度控制标签详解

### 7.1 PREC

PREC 设定 ENCUT 的默认值、FFT 网格与实空间投影精度（ROPT）。官方表格的要点：

| PREC | 默认 ENCUT | 波函数网格 NGX,Y,Z | 密度细网格 NGXF,Y,Z | ROPT（LREAL=Auto） | 适用 |
|---|---|---|---|---|---|
| Normal（默认） | max(ENMAX) | 覆盖到 $\frac{3}{2}G_{\text{cut}}$ | 2×NGX | −5×10⁻⁴ | 常规计算、大体系 |
| Accurate | max(ENMAX) | 覆盖到 $2G_{\text{cut}}$（无混叠误差） | 2×NGX | −2.5×10⁻⁴ | 高精度：力、应力、声子、弹性 |
| Single（VASP 6） | max(ENMAX) | 同 Accurate | = NGX（不用双网格） | −5×10⁻⁴ | 减小 CHGCAR 体积、STM 等 |
| Low / Medium / High | — | — | 由 ENAUG 启发式决定 | — | **已弃用**，仅为兼容旧输入 |

> [!IMPORTANT]
> PREC = Accurate **并不意味着结果已收敛**，它只是消除了 FFT 网格的混叠误差；ENCUT 仍需单独测试。另外，PREC 设定的默认值会被 INCAR 中显式给出的 ENCUT、NGX 等覆盖。

### 7.2 EDIFF 与 EDIFFG

| 场景 | EDIFF | EDIFFG | 说明 |
|---|---|---|---|
| 粗略筛选 | 1E-5 | −0.05 | — |
| 常规生产（优化、能量） | 1E-6 | −0.01 ~ −0.02 | 官方认为 1E-6 是最佳折中 |
| 收敛测试 | 1E-7 | — | 让 SCF 误差远小于测试判据 |
| 声子（有限位移/DFPT）、弹性常数 | 1E-7 ~ 1E-8 | 优化到 −1E-3 量级 | 力常数对噪声敏感 |
| 大体系 meta-GGA | 1E-6 | — | 官方提醒 1E-7/1E-8 对大体系 meta-GGA 很难达到 |

EDIFF 的含义：相邻两步的自由能变化与能带结构能变化**都**小于 EDIFF 时电子步停止。EDIFFG 为负值时按力判断离子收敛（所有原子受力的模都小于 |EDIFFG|），正值时按能量变化判断；默认值为 EDIFF×10。结构优化中电子收敛不够会导致力有噪声，进而引发 `ZBRENT` 一类线搜索错误（[第4章](04-relax-eos.md)）。

### 7.3 LREAL

| 取值 | 含义 | 建议 |
|---|---|---|
| `.FALSE.`（默认） | 投影在倒空间计算，精确 | 小晶胞（< 约 30 原子）；最终精确能量；声子、弹性 |
| `Auto`（或 `A`） | 实空间投影，自动优化投影算符 | 大体系（官方：> 约 30 原子）提速明显 |
| `On`、`.TRUE.` | 旧方案 | 官方不推荐 |

官方说明：实空间投影总会带来一个小的、通常是“每原子近似常数”的能量偏移。因此计算能量差时，**所有相关计算（如表面能中的板层与体相）必须使用相同的 ENCUT、PREC、LREAL、ROPT**；或者用 `LREAL = Auto` 优化结构，最后用 `LREAL = .FALSE.` 做一次单点能。显著提高 ENCUT 时，实空间投影可能不可靠，需要减小 ROPT 的绝对值（官方提到可能需要小到 1E-4 量级）。

### 7.4 LASPH

默认 `.FALSE.`，即 PAW 球内只计算梯度修正的球对称部分。官方指出以下情况**必须**设 `LASPH = .TRUE.`：f 元素（如 CeO₂）、所有 3d 元素（过渡金属氧化物）、第二周期的磁性原子（B–F），特别是使用 DFT+U、杂化泛函、meta-GGA（SCAN/r2SCAN）或 vdW-DF 时。它的额外开销很小，对含过渡金属的体系建议默认打开。

### 7.5 ADDGRID

`ADDGRID = .TRUE.` 在计算增广电荷时使用一个比细网格再密 2 倍（总点数 8 倍）的辅助网格，通常能减小力的噪声，常用于声子和高精度力的计算。官方提醒它并非总是有益、应谨慎使用；建议在声子计算前用小测试比较有无 ADDGRID 时的力和频率。

### 7.6 FFT 网格、混叠与力的漂移

检查 FFT 网格是否足够的两个简便方法：

```bash
grep "dimension x,y,z" OUTCAR          # NGX/NGY/NGZ 与 NGXF/NGYF/NGZF
grep -i "wrap" OUTCAR                   # 手动设置网格过小时会出现 wrap-around 警告
grep "total drift" OUTCAR | tail -1     # 力的漂移（总力应为零，漂移是数值误差）
```

官方建议：力的漂移绝不能超过力本身的大小，一般应小于你关心的力；漂移偏大时，增加 PREC/细网格或打开 ADDGRID。

### 7.7 精度档位速查

| 档位 | PREC | ENCUT | EDIFF | LREAL | 其他 | 典型用途 |
|---|---|---|---|---|---|---|
| 快速筛选 | Normal | max(ENMAX) | 1E-5 | Auto（大体系） | — | 初始结构、构型筛选 |
| 常规生产 | Accurate | 收敛值，变胞 ≥ 1.3×ENMAX | 1E-6 | 小体系 .FALSE.，大体系 Auto | LASPH（d/f 元素） | 优化、能量、能带 |
| 高精度 | Accurate | 收敛值 + 余量 | 1E-7 ~ 1E-8 | .FALSE. | ADDGRID（按需）、EDIFFG ≈ −1E-3 | 声子、弹性常数、介电 |

---

## 8. 力与应力的收敛

能量收敛不代表力和应力收敛。建议至少做两件事：

1. **位移结构的力测试**：如 4.2 节的 `POSCAR_displaced`，观察最大力随 ENCUT 和 k 点的变化，阈值可取 5 meV/Å（`plot_conv.py --thr-f 5`）；
2. **应力随 ENCUT 的变化**：`grep "external pressure" OUTCAR`，Pulay 应力在此直接体现。

对力常数类计算（有限位移声子、IBRION = 5/6、弹性常数），还要额外注意：

- 位移只有 0.01–0.015 Å 量级，产生的力很小，因此对数值噪声极其敏感：EDIFF 1E-7–1E-8、PREC = Accurate、LREAL = .FALSE.，必要时 ADDGRID；
- **蛋箱效应（egg-box effect）测试**：把整个结构刚性平移一个小量（如所有原子的 z 坐标加 0.1 Å），能量和力应保持不变；若变化明显，说明 FFT 网格（或 ENCUT）不足；
- 金属的声子计算要同时收敛 SIGMA 与 k 点；不要用 ISMEAR = −5 计算力。

详见 [第12章 声子](12-phonon.md) 与 [第13章 弹性常数](13-elastic.md)。

---

## 9. 自动化：bash + Python 批量测试

示例目录的组织方式可以直接套用到你的体系：

```text
03-convergence-si-al/
  common.sh            # natoms / outcar_row / run_or_submit / check_potcar
  plot_conv.py         # 读取 .dat，计算偏差、判断收敛、画图
  encut_test/          # INCAR.template + run_encut.sh
  kpoint_test/         # INCAR + run_kpts.sh (SCAN=mesh | kspacing)
  sigma_check_Al/      # INCAR.template + run_sigma.sh
```

三种运行模式（环境变量 `MODE`）：

| MODE | 行为 | 适用 |
|---|---|---|
| `local`（默认） | 在每个子目录中执行 `$VASP_CMD`（默认 `mpirun -np $NP vasp_std`），然后汇总 | 工作站、交互式节点 |
| `slurm` | 把 `JOB_TEMPLATE`（默认示例 01 的 job_cpu.slurm）复制到每个子目录并 `sbatch` | 集群 |
| `collect` | 不运行，只从已有 OUTCAR 汇总 | 作业结束后 |

改造成你自己的体系只需三步：替换 POSCAR、POTCAR；在 INCAR.template 中改成你的泛函和自旋设置；按需修改 `ENCUTS`、`MESHES`、`SIGMAS` 环境变量。例如：

```bash
ENCUTS="400 450 500 550 600 650 700" MODE=slurm bash run_encut.sh
MESHES="3 5 7 9 11" bash run_kpts.sh
```

> [!TIP]
> 仓库还提供了更通用的测试驱动 [`conv_test.sh`](../scripts/conv_test.sh) 与分析脚本 [`conv_analyze.py`](../scripts/conv_analyze.py)（分别用 `bash ../scripts/conv_test.sh -h` 和 `python ../scripts/conv_analyze.py -h` 查看用法）。GitHub 上的 [kavanase/vaspup2.0](https://github.com/kavanase/vaspup2.0) 也是一个流行的 ENCUT/k 点收敛测试流程，可对照参考。

对于更大规模的自动化（数十个体系同时测试），可以考虑 pymatgen 输入集 + atomate2/jobflow 工作流，见公共章节 [第6章 工作流自动化与高通量](../../Common/tutorials/06-workflow-automation.md)。

---

## 参数选择与收敛建议

| 参数 | 建议起点 | 收敛判据 | 备注 |
|---|---|---|---|
| ENCUT | 扫描 1.0–1.6×（及以上）max(ENMAX) | ΔE < 1 meV/atom；压强变化 < 1 kB | 变胞/EOS/弹性 ≥ 1.3×ENMAX；所有对比计算统一 |
| k 点（半导体） | KSPACING 0.25–0.3 Å⁻¹ 起扫 | ΔE < 1 meV/atom | Γ 中心；fcc/六方禁用偶数 MP |
| k 点（金属） | KSPACING 0.15–0.2 Å⁻¹ 起扫 | ΔE < 1 meV/atom，与 SIGMA 联合 | 越小的 SIGMA 需要越密的网格 |
| k 点（表面/2D） | 面内按体相换算，z = 1 | 同上 | 真空 ≥ 15 Å |
| ISMEAR/SIGMA（半导体） | 0 / 0.05 | \|`T*S`\| ≈ 0 | 或 −5（≥ 4 个 k 点） |
| ISMEAR/SIGMA（金属优化） | 1 / 0.2 | \|`T*S`\| < 1 meV/atom | 不用 −5 做优化 |
| ISMEAR（DOS/精确能量） | −5 | — | Γ 中心，≥ 4 个 k 点 |
| PREC | Accurate | — | Normal 仅用于粗筛 |
| EDIFF | 1E-6（生产）/ 1E-7（测试）/ 1E-8（声子、弹性） | — | 大体系 meta-GGA 用 1E-6 |
| LREAL | 小体系 .FALSE.；> 30 原子 Auto | — | 能量差计算中保持一致 |
| LASPH | d/f 元素、DFT+U、meta-GGA、vdW-DF、杂化：.TRUE. | — | 开销小 |
| ADDGRID | 声子、高精度力：按需 .TRUE. | 力的噪声/漂移减小 | 谨慎使用，先测试 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| 能量随 ENCUT 不单调、来回跳 | EDIFF 太松；不同 ENCUT 用了不同 PREC/LREAL；读了旧 WAVECAR | EDIFF = 1E-7；固定其他参数；`ISTART = 0` 并删除旧文件 |
| 能量已收敛但优化出的晶格常数随 ENCUT 明显变化 | Pulay 应力 | ENCUT ≥ 1.3×ENMAX；看压强的收敛；变胞后重新优化 |
| 金属的能量随 k 点振荡、迟迟不收敛 | SIGMA 太小；或用了 ISMEAR = 0 且 SIGMA 很小 | 改用 ISMEAR = 1、SIGMA ≈ 0.2 并检查 `T*S`；加密 k 点 |
| `EENTRO` 每原子超过 1 meV | SIGMA 太大 | 减小 SIGMA，同时加密 k 点；比较能量用 `energy(sigma->0)` |
| 半导体的带隙、占据数或声子异常 | 用了默认的 `ISMEAR = 1` | 显式设置 `ISMEAR = 0`、`SIGMA = 0.05` 或 `ISMEAR = -5` |
| 报错中含 `Tetrahedron method fails` | ISMEAR = −5 但 k 点 < 4（如 Γ 点、表面 1×1×1） | 改为 ISMEAR = 0 |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | fcc/六方用了偶数 MP 网格或不相容的网格 | 改为 Γ 中心网格 |
| `Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accurate calculation.` | 手动设置的网格或 PREC 太低 | 使用 PREC = Accurate，删除手动 NGX 设置 |
| `total drift` 与力同量级 | FFT 网格/ENCUT 不足、SCF 未充分收敛 | PREC = Accurate、提高 ENCUT、减小 EDIFF，必要时 ADDGRID |
| 改变进程数后能量有细微差别 | NBANDS 被自动调整、SCF 收敛容差 | 差别应在 EDIFF 量级；需要严格一致时固定 NBANDS 与并行设置 |
| 收敛测试的 KSPACING 网格“跳格” | 上取整规则导致临界值附近网格突变 | 查看 OUTCAR 中实际网格；生产计算固定显式网格 |

---

## 练习

1. 运行 `encut_test`，分别以能量（1 meV/atom）和压强（1 kB）为判据确定 Si 的 ENCUT，比较两者差别并解释原因。
2. 用 `POSCAR_displaced` 做 ENCUT 与 k 点测试，确定使最大力误差 < 5 meV/Å 的设置。
3. 在 `kpoint_test` 中比较 SCAN=mesh 与 SCAN=kspacing 两种方式：画出 ΔE 随 NKPTS 的变化，二者是否落在同一条曲线上？
4. 对 Al 运行 `run_sigma.sh`，找出满足 `|T*S|` < 1 meV/atom 且 |ΔE0| < 1 meV/atom、计算量最小的 (ISMEAR, SIGMA, N) 组合。
5. （进阶）对 Si 的 8 原子立方超胞做蛋箱效应测试：把所有原子沿 z 平移 0.05、0.10、0.15 Å，比较 `PREC = Normal` 与 `PREC = Accurate` 下能量和力的变化。

---

## 参考资料

**官方文档**

- ENCUT：https://www.vasp.at/wiki/index.php/ENCUT ；Energy cutoff and FFT meshes：https://www.vasp.at/wiki/index.php/Energy_cutoff_and_FFT_meshes
- PREC：https://www.vasp.at/wiki/index.php/PREC ；EDIFF：https://www.vasp.at/wiki/index.php/EDIFF ；EDIFFG：https://www.vasp.at/wiki/index.php/EDIFFG
- LREAL：https://www.vasp.at/wiki/index.php/LREAL ；LASPH：https://www.vasp.at/wiki/index.php/LASPH ；ADDGRID：https://www.vasp.at/wiki/index.php/ADDGRID
- ISMEAR：https://www.vasp.at/wiki/index.php/ISMEAR ；SIGMA：https://www.vasp.at/wiki/index.php/SIGMA ；Smearing technique：https://www.vasp.at/wiki/index.php/Smearing_technique ；EFERMI：https://www.vasp.at/wiki/index.php/EFERMI ；SMEARINGS：https://www.vasp.at/wiki/index.php/SMEARINGS
- KPOINTS：https://www.vasp.at/wiki/index.php/KPOINTS ；KSPACING：https://www.vasp.at/wiki/index.php/KSPACING
- VASP 官方教程（含收敛测试练习）：https://www.vasp.at/tutorials/latest/

**相关 GitHub 仓库**（仅链接）

- [kavanase/vaspup2.0](https://github.com/kavanase/vaspup2.0)：ENCUT/k 点能量与力的收敛测试 bash 流程。
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：输入集（MPRelaxSet 等）中的 ENCUT/k 点密度约定可作参考。
- [materialsproject/custodian](https://github.com/materialsproject/custodian)：收录大量真实报错及自动修复策略。
- [SMTG-Bham/surfaxe](https://github.com/SMTG-Bham/surfaxe)：表面板层层数与真空的收敛测试工具。
- [cndaqiang/DFT-EXERCISES](https://github.com/cndaqiang/DFT-EXERCISES)：Sholl《DFT: A Practical Introduction》习题的中文解答，含收敛测试练习。

**关键文献**

- H. J. Monkhorst, J. D. Pack, Phys. Rev. B 13, 5188 (1976). DOI: 10.1103/PhysRevB.13.5188
- M. Methfessel, A. T. Paxton, Phys. Rev. B 40, 3616 (1989). DOI: 10.1103/PhysRevB.40.3616
- P. E. Blöchl, O. Jepsen, O. K. Andersen, Phys. Rev. B 49, 16223 (1994). DOI: 10.1103/PhysRevB.49.16223
- G. P. Francis, M. C. Payne, J. Phys.: Condens. Matter 2, 4395 (1990)（Pulay 应力的有限基组修正）. DOI: 10.1088/0953-8984/2/19/007
- K. Lejaeghere et al., Science 351, aad3000 (2016)（DFT 程序间的可重复性比较）. DOI: 10.1126/science.aad3000
- P. Wisesa, K. A. McGill, T. Mueller, Phys. Rev. B 93, 155109 (2016)（广义 k 点网格）. DOI: 10.1103/PhysRevB.93.155109

[← 上一章](02-input-output-files.md) | [返回目录](../README.md) | [下一章 →](04-relax-eos.md)
