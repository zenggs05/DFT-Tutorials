# 第14章 振动分析、红外与热化学（Vibrational Analysis, IR Intensities & Thermochemistry with CP2K）

本章讲解 CP2K 的 `RUN_TYPE VIBRATIONAL_ANALYSIS`：用有限差分构建 Hessian，得到简谐振动频率、红外（IR）强度、零点能（ZPE）和热化学量（熵、焓、Gibbs 自由能校正）。以 H2O（气相分子：频率 + IR + 热化学）、Pt(111) 表面吸附的 CO（固定衬底的部分 Hessian：吸附物的 ZPE 与振动熵）和[第13章](13-neb-dimer.md) H/Cu(111) 的过渡态（只有一个虚频的验证）三个例子，逐一说明 `DX`、`NPROC_REP`、`INTENSITIES`、`THERMOCHEMISTRY`、`FULLY_PERIODIC`、`&MOLDEN_VIB`、`&MODE_SELECTIVE` 等关键字，以及虚频的诊断与处理、CP2K 热化学输出的局限（对称数恒为 1、默认温度 273.15 K）和 2026.2 新增的未质量加权 Hessian 输出。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 理解有限差分 Hessian、质量加权与简谐频率之间的关系，知道计算量随可动原子数如何增长；
- 为振动分析准备“严格收敛的几何 + 高精度 SCF”，并用同一个 `&FORCE_EVAL` 串起优化和频率计算；
- 写出 `VIBRATIONAL_ANALYSIS` 的完整输入，用 `NPROC_REP` 并行处理各个位移结构；
- 读懂 `VIB|` 输出：频率、IR 强度、约化质量、力常数、热化学数据；用 molden 文件可视化振动模式；
- 用固定原子做表面吸附物的部分 Hessian，计算 ZPE、振动熵和振动自由能；
- 用频率计算验证过渡态，诊断并处理虚频；
- 知道 CP2K 热化学模块的适用范围，以及何时需要自己后处理。

| 任务 | 做法 | 示例 |
|---|---|---|
| 气相分子：频率、IR、ZPE、G | 全部原子可动，`INTENSITIES T` + `&MOMENTS`，`THERMOCHEMISTRY T` | `h2o_opt_tight.inp` → `h2o_vib_ir.inp` |
| 表面吸附物：ZPE、振动熵 | 固定衬底原子（部分 Hessian），简谐极限后处理 | `co_pt111_opt.inp` → `co_on_pt_partial_hessian.inp` |
| 过渡态验证 | 与 NEB/Dimer 同一 `FORCE_EVAL`，检查只有一个虚频 | `ts_freq_check.inp` |
| 大体系里只关心某几个模式 | `&MODE_SELECTIVE`（Davidson 迭代跟踪指定模式） | §5.6 |
| 晶体 Γ 点声子 | `FULLY_PERIODIC T`；色散关系请用 phonopy 等 | §6 |
| 含非谐效应的红外光谱 | AIMD 偶极自相关函数 | [第12章](12-aimd.md) |

示例输入集：[../examples/14-vibration-thermo/README.md](../examples/14-vibration-thermo/README.md)。

---

## 2. 原理速览

### 2.1 简谐近似与质量加权 Hessian

在势能极小点（或鞍点）$\mathbf R_0$ 附近把能量展开到二阶：

$$
E(\mathbf R) \approx E(\mathbf R_0) + \frac12\sum_{ij} H_{ij}\,\Delta x_i\,\Delta x_j,\qquad H_{ij} = \frac{\partial^2 E}{\partial x_i\,\partial x_j}
$$

（一阶项因 $\nabla E=0$ 而消失，这正是“必须先严格优化”的原因。）引入质量加权坐标 $q_i = \sqrt{m_i}\,\Delta x_i$，质量加权 Hessian 为 $\tilde H_{ij} = H_{ij}/\sqrt{m_i m_j}$。它的本征值 $\lambda_k = \omega_k^2$ 给出角频率，波数为

$$
\tilde\nu_k = \frac{\omega_k}{2\pi c} = \frac{\sqrt{\lambda_k}}{2\pi c}
$$

$\lambda_k<0$ 时 $\omega_k$ 是虚数，即**虚频**。CP2K 把它打印成负的波数。

### 2.2 有限差分

CP2K 用解析力做中心差分来构造 Hessian：

$$
H_{ij} \approx -\frac{F_i(x_j+\Delta) - F_i(x_j-\Delta)}{2\Delta}
$$

其中 $\Delta$ 就是 `DX`（默认 0.01 bohr）。每个可动坐标正、负各位移一次，再加一次参考结构，共 $6N_\mathrm{mob}+1$ 次能量 + 力计算（$N_\mathrm{mob}$ 为可动原子数），算完后再把矩阵对称化。误差有两个来源：

- 截断误差 $\propto \Delta^2$：`DX` 越大越不准，尤其是软模式和非谐性强的模式；
- 数值噪声 $\propto \delta F/\Delta$：SCF 和积分网格的力噪声 $\delta F$ 被 $1/\Delta$ 放大，`DX` 越小越吵。

所以 `DX` 取 0.005–0.02 bohr，同时把力噪声压下去（`EPS_SCF` 1E-7 至 1E-8、较紧的 `EPS_DEFAULT`、足够高的 CUTOFF），比单纯调 `DX` 更有效。

### 2.3 平动与转动的去除

孤立分子的 Hessian 有 6 个（线性分子 5 个）零频模式，对应整体平动和转动。数值误差会让它们变成几十 cm⁻¹ 的“假频率”，并和真实的低频模式混在一起。CP2K 先由惯性主轴构造平动/转动子空间，在其正交补空间里对角化 Hessian（输出中的 `Internal Low frequencies`）：

- 默认：去掉平动和转动；
- `FULLY_PERIODIC T`：不去掉转动（周期性晶体中整体转动不是对称操作）；
- **有固定原子时**：CP2K 不做任何平动/转动投影（固定原子本身已破坏这两种对称性），直接对角化可动原子的 Hessian 块。

### 2.4 红外强度

第 $k$ 个简正模式的 IR 吸收强度正比于偶极矩对简正坐标导数的平方：

$$
I_k \propto \left|\frac{\partial \boldsymbol\mu}{\partial Q_k}\right|^2
$$

CP2K 在每个位移结构上计算偶极矩（由 `&DFT/&PRINT/&MOMENTS` 提供），再与振动模式组合，输出单位为 km/mol。孤立体系用普通偶极算符（`PERIODIC F`），周期体系用 Berry 相公式（`PERIODIC T`，默认）。Raman 强度需要极化率导数，要求额外设置 `&PROPERTIES/&LINRES/&POLAR`（见 §4.8）。

### 2.5 热化学：RRHO 与简谐极限

理想气体分子在刚性转子–简谐振子（RRHO）模型下，自由能校正由平动、转动、振动和电子四部分配分函数给出。振动部分（$x_k = h\nu_k/k_\mathrm BT$）：

$$
\mathrm{ZPE} = \sum_k \tfrac12 h\nu_k,\qquad
U_\mathrm{vib} = \sum_k\left[\tfrac12 h\nu_k + \frac{h\nu_k}{e^{x_k}-1}\right],\qquad
S_\mathrm{vib} = k_\mathrm B\sum_k\left[\frac{x_k}{e^{x_k}-1} - \ln\!\left(1-e^{-x_k}\right)\right]
$$

转动熵含对称数 $\sigma$：$S_\mathrm{rot} \propto R\ln(1/\sigma)$；电子熵 $S_\mathrm{el} = R\ln g$（$g$ 为基态简并度，如三重态 O2 的 $g=3$）。

**表面吸附物**通常采用**简谐极限**：吸附后的所有自由度（包括受阻平动、受阻转动）都当作振动，只算 $\mathrm{ZPE}$、$U_\mathrm{vib}$、$S_\mathrm{vib}$，不计平动和转动配分函数。

> [!WARNING]
> 低于 50–100 cm⁻¹ 的软模式会让 $S_\mathrm{vib}$ 发散式增大（$\nu\to0$ 时 $-\ln(1-e^{-x})\to\infty$），而这类模式的简谐近似本身就不可靠。CP2K 2026.2 的输入参考也明确说明其热化学基于 RRHO，在存在很低频率的柔性体系中会失效。常见补救办法：把低频抬高到一个下限（如 50 cm⁻¹）；对分子用 quasi-RRHO（如 Grimme 的方案，可借助 Shermo 等程序）；对表面吸附物用受阻平动/受阻转动模型（Sprowl, Campbell & Árnadóttir, 2016）。

---

## 3. 前提：严格优化的几何与高精度 SCF

在非驻点上做振动分析，梯度项不为零，得到的“频率”会有系统误差，最常见的症状是出现几十 cm⁻¹ 的小虚频。建议：

| 设置 | 推荐值 | 说明 |
|---|---|---|
| 优化 `MAX_FORCE` | 分子 ≤ 1E-5，表面/大体系 ≤ 1E-4 Ha/bohr | 默认 4.5E-4 对频率计算太松 |
| `EPS_SCF` | 1E-7 至 1E-8（OT 时内外层一致） | 力噪声直接进入 Hessian |
| `&QS/EPS_DEFAULT` | 1E-12 至 1E-14 | 积分筛选误差也是力噪声的来源 |
| CUTOFF | 比单点能收敛值再提高 100–200 Ry | 网格造成的“蛋盒效应”会产生假的平动频率 |
| `DX` | 0.01 bohr（默认）；软模式多时试 0.005 与 0.02 对比 | §2.2 |
| 一致性 | 优化与频率用**同一个** `&FORCE_EVAL`（基组、泛函、色散、CUTOFF、固定原子） | 否则优化得到的“极小点”在新设置下并不是极小点 |

本章示例都用 `@INCLUDE` 共享 `&FORCE_EVAL` 来保证最后一条。H2O 的共享文件 `h2o_force_eval.inc`：

```fortran
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-14
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-8
      MAX_SCF 30
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-8
        MAX_SCF 20
      &END OUTER_SCF
      &PRINT
        &RESTART
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &MOMENTS ON
        PERIODIC F
        REFERENCE COM
      &END MOMENTS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.0 12.0 12.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
      ! MASS 2.014102
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

- TZV2P 基组和 600 Ry 比结构优化常用的 DZVP/400 Ry 更稳妥，H2O 很小，代价可以忽略；
- 孤立分子用 `PERIODIC NONE` + MT 泊松求解器，盒子是分子尺寸的 2 倍以上（[第6章](06-periodic-surface-kpoints.md)）；
- `&MOMENTS` 为 IR 强度提供偶极矩；孤立体系用 `PERIODIC F`、以质心为参考点（`REFERENCE COM`）。

先做紧优化（`h2o_opt_tight.inp`，`&MOTION/&GEO_OPT` 中 `MAX_FORCE 1.0E-5`、`RMS_FORCE 5.0E-6`、`MAX_DR 1.0E-4`、`RMS_DR 5.0E-5`）：

```bash
cd CP2K/examples/14-vibration-thermo
mpirun -np 8 cp2k.psmp -i h2o_opt_tight.inp -o h2o_opt_tight.out
grep "GEOMETRY OPTIMIZATION COMPLETED" h2o_opt_tight.out
grep "OPT| Maximum gradient " h2o_opt_tight.out | tail -1
bash last_frame.sh h2o_opt-pos-1.xyz h2o_opt.xyz
```

---

## 4. 气相分子：H2O 的频率、红外强度与热化学

### 4.1 完整输入

`h2o_vib_ir.inp`：

```fortran
@SET COORD ${COORD-h2o_opt.xyz}

&GLOBAL
  PROJECT h2o_vib
  RUN_TYPE VIBRATIONAL_ANALYSIS
  PRINT_LEVEL MEDIUM
&END GLOBAL

@INCLUDE 'h2o_force_eval.inc'

&VIBRATIONAL_ANALYSIS
  DX 0.01
  NPROC_REP ${NPROC_REP-2}
  INTENSITIES T
  THERMOCHEMISTRY T
  TC_TEMPERATURE [K] 298.15
  TC_PRESSURE [Pa] 101325.0
  &PRINT
    &MOLDEN_VIB ON
    &END MOLDEN_VIB
  &END PRINT
&END VIBRATIONAL_ANALYSIS
```

`&VIBRATIONAL_ANALYSIS` 是与 `&FORCE_EVAL`、`&MOTION` 并列的**顶层段**，不要写进 `&MOTION`。

| 关键词 | 本例 | 默认值 | 说明 |
|---|---|---|---|
| `RUN_TYPE`（`&GLOBAL`） | `VIBRATIONAL_ANALYSIS` | — | 别名 `NORMAL_MODES` |
| `DX` | 0.01 | 0.01 bohr | 有限差分步长 |
| `NPROC_REP` | 2 | 1 | 每个“副本”（同时计算的一个位移结构）使用的 MPI 进程数，见 §4.2 |
| `PROC_DIST_TYPE` | 不写 | `BLOCKED` | 进程到副本的映射方式 |
| `INTENSITIES` | `T` | `F` | 计算 IR（和 Raman）强度。IR 需要 `&DFT/&PRINT/&MOMENTS`，Raman 需要 `&PROPERTIES/&LINRES/&POLAR` |
| `THERMOCHEMISTRY` | `T` | `F` | 按 RRHO 模型输出热化学数据，**只适用于气相孤立分子**，且只考虑 Γ 点 |
| `TC_TEMPERATURE` | 298.15 K | **273.15 K** | 热化学温度。默认值不是室温，务必显式设置 |
| `TC_PRESSURE` | 101325 Pa | 101325 Pa | 平动熵用的压力 |
| `FULLY_PERIODIC` | 不写 | `F` | `T` 时不从 Hessian 中去掉转动，用于周期性晶体 |
| `&PRINT/&MOLDEN_VIB` | `ON` | 打印级别 `LOW` | 写出 `<PROJECT>-VIBRATIONS-1.mol`，可视化振动模式 |
| `&PRINT/&HESSIAN` | 不写 | `LOW` | 把（质量加权的）Hessian 写成二进制文件 `*.hess` |
| `&PRINT/&CARTESIAN_EIGS` | 不写 | — | 把笛卡尔 Hessian 的本征值和本征矢写成二进制 `*.eig`，可用于 MD 初始化 |
| `&PRINT/&PROGRAM_RUN_INFO` | 不写 | `MEDIUM` | 控制各副本能量/力、Hessian 矩阵等详细信息的打印 |
| `&MODE_SELECTIVE` | 不写 | — | 只跟踪少数模式，见 §5.6 |

### 4.2 并行：`NPROC_REP`

全部 MPI 进程被分成 `N_MPI / NPROC_REP` 个副本组，每组同时计算一个位移结构。H2O 有 9 个坐标，正负位移各 9 个结构。用 18 个进程、`NPROC_REP 2` 时有 9 个副本组，正向位移一轮、负向位移一轮、加上参考结构，共约 3 轮单点计算的墙钟时间。

| 体系 | 可动坐标数 | 建议 |
|---|---|---|
| 小分子 | 9–30 | `NPROC_REP` 取 1–4，副本组数尽量等于坐标数 |
| 吸附物（部分 Hessian） | 6–30 | 每个副本用单点计算的“甜点”进程数，组数取坐标数的约数 |
| 大体系全 Hessian | 数百以上 | 计算量为 $6N$ 次力计算，考虑部分 Hessian 或 `&MODE_SELECTIVE` |

示例的 [job_slurm.sh](../examples/14-vibration-thermo/job_slurm.sh) 会根据 `SLURM_NTASKS` 自动设置 `NPROC_REP`。

### 4.3 运行与输出

```bash
mpirun -np 18 cp2k.psmp -i h2o_vib_ir.inp -o h2o_vib_ir.out

grep -A4 "Minimum Structure - Energy and Forces" h2o_vib_ir.out   # 参考结构的力应接近 0
grep "VIB|Frequency"      h2o_vib_ir.out
grep "VIB|IR int"         h2o_vib_ir.out
grep -A12 "THERMOCHEMICAL DATA" h2o_vib_ir.out
ls h2o_vib-VIBRATIONS-1.mol
```

输出中与振动分析相关的内容依次是：

1. `VIB| Vibrational Analysis Info` 和每个位移结构的 `VIB| REPLICA Nr. ... - Energy and Forces for particle: ... coordinate: ...`（`PRINT_LEVEL MEDIUM` 时）；
2. `VIB| Minimum Structure - Energy and Forces:`：参考结构的总能量和各原子受力。**先看这里的力是否足够小**；
3. Hessian 矩阵。2026.2 起先打印未质量加权的 `VIB| Hessian (before multiplying by 1E6/(sqrt(mass_i)*sqrt(mass_j)))`（原子单位，Ha/bohr²），再打印 `VIB| Hessian in cartesian coordinates (mass weighted)`；更早的版本只打印后者；
4. `VIB| Cartesian Low frequencies ---`（投影前的本征值）和 `VIB| Internal  Low frequencies ---`（去掉平动/转动后的本征值，原子单位）；
5. `NORMAL MODES - CARTESIAN DISPLACEMENTS` 表：每三个模式一组；
6. `NORMAL MODES - THERMOCHEMICAL DATA`（`THERMOCHEMISTRY T` 时）。

简正模式表的样子（示意输出，数值仅供参考）：

```text
 VIB|                        NORMAL MODES - CARTESIAN DISPLACEMENTS
 VIB|
 VIB|                     1                     2                     3
 VIB|Frequency (cm^-1)  1.5874312201E+03     3.6795403311E+03     3.7841098725E+03
 VIB|IR int (KM/Mole)   6.5120391137E+01     2.1450018720E+00     4.6631729012E+01
 VIB|Red.Masses (a.u.)  1.0830016720E+00     1.0453119002E+00     1.0816610734E+00
 VIB|Frc consts (a.u.)  2.4987311009E-02     1.0221570013E-01     1.1088320761E-01
  ATOM  EL           X     Y     Z         X     Y     Z         X     Y     Z
     1  O         0.00  0.00 -0.07      0.00  0.00  0.05     -0.07  0.00  0.00
     2  H         0.42  0.00  0.56     -0.58  0.00 -0.40      0.56  0.00 -0.43
     3  H        -0.42  0.00  0.56      0.58  0.00 -0.40      0.56  0.00  0.43
```

各行含义：`Frequency` 为波数（负值即虚频）；`IR int` 为红外强度（km/mol，无偶极信息时这一行不出现，并打印 `VIB| No IR intensities available. Check input`）；`Red.Masses` 为约化质量（原子质量单位）；`Frc consts` 为力常数（原子单位）。下面的矩阵是每个模式的笛卡尔位移（归一化）。

> [!NOTE]
> CP2K 2026.1 的发布说明中提到修复了振动分析中的力常数（`Frc consts` 一行）。如果你用更早的版本，请不要依赖这一行的数值；频率本身不受影响。

**判断是否成功**：H2O 应有 3 个实频模式（弯曲约 1600 cm⁻¹，对称和反对称 O–H 伸缩在 3600–3800 cm⁻¹ 附近），没有虚频。参考值：气相 H2O 的实验基频为 1595、3657、3756 cm⁻¹（NIST Chemistry WebBook）。注意计算得到的是**简谐**频率，而实验值是包含非谐效应的**基频**，O–H 伸缩的非谐修正约 100–200 cm⁻¹；PBE 本身也有系统误差。二者的差别不能单纯归咎于数值设置。

**可视化**：`h2o_vib-VIBRATIONS-1.mol` 是标准 molden 格式（`[FREQ]`、`[FR-COORD]`、`[FR-NORM-COORD]`、启用强度时还有 `[INT]`），可以用 Molden、Jmol、Avogadro 等程序打开并播放振动动画。

### 4.4 画红外光谱

把每个模式的 IR 强度展宽成洛伦兹线形，就得到可以和实验对比的“理论 IR 谱”。下面的脚本从输出中读取频率和 IR 强度（只依赖 numpy 和 matplotlib）：

```python
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

num = r"[-+]?\d+\.\d*(?:[Ee][-+]?\d+)?"
freq, inten = [], []
for line in open("h2o_vib_ir.out"):
    if line.strip().startswith("VIB|Frequency (cm^-1)"):
        freq += [float(x) for x in re.findall(num, line.split(")", 1)[1])]
    elif line.strip().startswith("VIB|IR int (KM/Mole)"):
        inten += [float(x) for x in re.findall(num, line.split(")", 1)[1])]
freq, inten = np.array(freq), np.array(inten)

x = np.linspace(400, 4200, 4000)
gamma = 10.0                               # HWHM in cm^-1
spec = sum(I * gamma / np.pi / ((x - f) ** 2 + gamma ** 2) for f, I in zip(freq, inten) if f > 0)
plt.plot(x, spec)
plt.vlines(freq[freq > 0], 0, inten[freq > 0] / inten.max() * spec.max(), color="gray", lw=0.8)
plt.xlabel("wavenumber (cm$^{-1}$)"); plt.ylabel("IR intensity (arb. units)")
plt.gca().invert_xaxis(); plt.tight_layout(); plt.savefig("h2o_ir.png", dpi=200)
```

（若输出中有多个振动分析块，请先截取最后一块。）如需更真实的谱形（温度展宽、非谐、溶剂效应），请用 AIMD 偶极自相关函数计算，见 [第12章](12-aimd.md)。

### 4.5 热化学数据及其局限

`THERMOCHEMISTRY T` 时的输出（示意输出，数值仅供参考；2026.2 格式）：

```text
 VIB|                        NORMAL MODES - THERMOCHEMICAL DATA
 VIB|              [q = gamma only, rigid-rotor harmonic oscillator (RRHO) model]

 VIB|       Symmetry number:                                                   1
 VIB|       Temperature [K]:                                              298.15
 VIB|       Pressure [Pa]:                                             101325.00

 VIB|       Electronic energy (U) [kJ/mol]:                      -45137.12345678
 VIB|       Zero-point correction [kJ/mol]:                          54.81234567
 VIB|       Entropy [kJ/(mol K)]:                                     0.19234567
 VIB|       Enthalpy correction (H-U) [kJ/mol]:                      64.73456789
 VIB|       Gibbs energy correction [kJ/mol]:                         7.38345678
 VIB|       Heat capacity [kJ/(mol*K)]:                               0.03345678
```

各项含义（按 CP2K 2026.2 源码 `vibrational_analysis.F` 中的实现）：

| 输出项 | 含义 |
|---|---|
| `Electronic energy (U)` | 参考结构的电子总能量，换算成 kJ/mol |
| `Zero-point correction` | $\mathrm{ZPE}=\sum\tfrac12 h\nu_k$ |
| `Entropy` | $S = S_\mathrm{trans}+S_\mathrm{rot}+S_\mathrm{vib}+S_\mathrm{el}$ |
| `Enthalpy correction (H-U)` | $E_\mathrm{rot}+H_\mathrm{trans}+U_\mathrm{vib}$（$U_\mathrm{vib}$ 已含 ZPE，$H_\mathrm{trans}=\tfrac52 RT$ 含 $pV$） |
| `Gibbs energy correction` | $(H-U) - TS$ |

于是 $G(T) = E_\mathrm{elec} + G_\mathrm{corr}$，$H(T) = E_\mathrm{elec} + (H-U)$。使用时注意三点局限：

1. **对称数恒为 1**（源码中写死 `sym_num = 1`）。H2O（$C_{2v}$）的 $\sigma=2$，CH4 为 12，NH3 为 3，苯为 12。真实转动熵应减去 $R\ln\sigma$，即 $G$ 要加上 $RT\ln\sigma$。H2O 在 298.15 K：$RT\ln 2 = 1.72$ kJ/mol（$S$ 减少 5.76 J mol⁻¹ K⁻¹）；
2. **电子简并度恒为 1**（调用时 `spin = 1`）。开壳层分子（如三重态 O2）应再加 $S_\mathrm{el}=R\ln(2S+1)$，即 $G$ 减去 $RT\ln 3$；
3. **只对孤立分子有意义**：只考虑 Γ 点，基于 RRHO；有固定原子时 CP2K 不输出这一块（§5）。若存在虚频，CP2K 也会把它的负值代入求和，所以过渡态的热化学必须自己去掉虚频后重算。

> [!TIP]
> 需要 quasi-RRHO、正确对称数、不同温度/压力扫描时，可以把 molden 文件或频率列表交给专门的程序处理（例如 Shermo），或用示例目录中的 `harmonic_thermo.py` 算振动部分、再手工加上平动/转动项。本仓库的 [vib_thermo.py](../scripts/vib_thermo.py) 也能解析振动输出并计算 ZPE、U、S、G（用 `python ../scripts/vib_thermo.py -h` 查看用法）。

### 4.6 同位素效应

改变质量只改变质量加权，不改变势能面。两种做法：

- 在 `&KIND` 中设置 `MASS`（如 `&KIND H` 中写 `MASS 2.014102` 就把 H 换成 D），重新运行振动分析；
- **2026.2 起**，输出中先打印未质量加权的 Hessian（`VIB| Hessian (before multiplying by 1E6/(sqrt(mass_i)*sqrt(mass_j)))`），可以把它取出来，离线用任意质量重新加权、对角化，一次计算得到所有同位素体的频率。

### 4.7 计算量与重启

振动分析没有“中途重启”：每个位移结构都要从头做 SCF（CP2K 会用参考结构的波函数作为初猜）。因此：

- 优化后保留 `<PROJECT>-RESTART.wfn`，在振动分析的 `&DFT` 中设 `WFN_RESTART_FILE_NAME` 并在 `&SCF` 中设 `SCF_GUESS RESTART`，参考结构的 SCF 会快很多（[第5章](05-scf-strategy.md)）；
- 大体系优先考虑部分 Hessian（§5）或 `&MODE_SELECTIVE`（§5.6）。

### 4.8 Raman（简述）

`INTENSITIES T` 时，如果每个位移结构还算了极化率（需要 `&FORCE_EVAL/&PROPERTIES/&LINRES/&POLAR`，线性响应计算，代价明显高于单点），输出还会多出 `VIB|Raman (A^4/amu)` 和两行退偏比 `VIB|Depol Ratio (P)`、`VIB|Depol Ratio (U)`；没有极化率时打印 `VIB| No Raman intensities available. Check input`。具体写法请参考官方回归测试 `tests/QS/regtest-gpw-2-1/N2_raman.inp`，并先用小分子确认设置无误。

---

## 5. 表面吸附物：部分 Hessian、ZPE 与振动熵

### 5.1 为什么只算吸附物

对 27 个 Pt 原子 + CO 的 slab 做全 Hessian 需要 $6\times29=174$ 次力计算，而真正关心的 CO 振动与衬底声子耦合很弱。通常的做法是**固定衬底、只位移吸附物**（部分 Hessian 近似）：只需 $6\times2=12$ 次力计算，得到 6 个模式。

在 CP2K 中，这不需要特殊关键字：`&MOTION/&CONSTRAINT/&FIXED_ATOMS` 中列出的原子不会被位移，Hessian 只包含可动原子。源码层面的两个后果需要知道：

- 有固定原子时**不做平动/转动投影**，可动原子的全部 $3N_\mathrm{mob}$ 个模式都输出。吸附物的受阻平动/受阻转动会以低频模式出现，这正是简谐极限所需要的；
- 有固定原子时**不输出 `THERMOCHEMISTRY` 块**（即使设了 `THERMOCHEMISTRY T`）。吸附物的 ZPE 和振动熵需要自己从频率计算，见 §5.4。

### 5.2 输入

共享的 `co_pt111_force_eval.inc` 与[第13章](13-neb-dimer.md)的 Cu 例子类似：PBE + D3(BJ)，`DZVP-MOLOPT-PBE-GTH-q18/q4/q6`，对角化 + Broyden + Fermi–Dirac 展宽，`EPS_SCF 1.0E-8`，CUTOFF 500 Ry，偶极校正。先优化（固定底部两层 1–18）：

```bash
mpirun -np 64 cp2k.psmp -i co_pt111_opt.inp -o co_pt111_opt.out
bash last_frame.sh co_pt111_opt-pos-1.xyz co_pt111_opt.xyz
```

`co_pt111_opt.inp` 的 `&GEO_OPT` 用 `MAX_FORCE 1.0E-4`、`RMS_FORCE 5.0E-5`。然后是部分 Hessian 输入 `co_on_pt_partial_hessian.inp`：

```fortran
@SET COORD ${COORD-co_pt111_opt.xyz}

&GLOBAL
  PROJECT co_pt111_vib
  RUN_TYPE VIBRATIONAL_ANALYSIS
  PRINT_LEVEL MEDIUM
&END GLOBAL

@INCLUDE 'co_pt111_force_eval.inc'

&VIBRATIONAL_ANALYSIS
  DX 0.01
  NPROC_REP ${NPROC_REP-16}
  &PRINT
    &MOLDEN_VIB ON
    &END MOLDEN_VIB
  &END PRINT
&END VIBRATIONAL_ANALYSIS

&MOTION
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..27
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

注意：优化时固定 1–18、振动分析时固定 1–27（全部 Pt）。这是允许的，因为振动分析只要求“可动原子处在驻点”，而 C、O 在优化中是放开的。反过来不行：振动分析中可动的原子必须在优化中也可动。

```bash
mpirun -np 96 cp2k.psmp -i co_on_pt_partial_hessian.inp -o co_on_pt_partial_hessian.out
grep "VIB|Frequency" co_on_pt_partial_hessian.out
```

### 5.3 结果解读

应得到 6 个模式：一个高频的 C–O 伸缩（2000 cm⁻¹ 上下），一个 Pt–C 伸缩（几百 cm⁻¹），以及两对简并或近简并的受阻转动、受阻平动（几十到几百 cm⁻¹）。参考：气相 CO 的实验基频为 2143 cm⁻¹（NIST）；吸附到金属顶位后 C–O 伸缩红移，Pt(111) 顶位 CO 的实验 IR 峰在 2100 cm⁻¹ 附近（随覆盖度变化，近似值）。如果某个低频模式是小虚频，通常是优化不够紧或数值噪声，见 §7。

### 5.4 吸附物的 ZPE、振动熵与自由能

示例目录中的 `harmonic_thermo.py` 按 §2.5 的公式在简谐极限下计算（只用标准库），可以读主输出，也可以读 molden 文件：

```bash
python3 harmonic_thermo.py co_on_pt_partial_hessian.out                  # 298.15 K
python3 harmonic_thermo.py co_pt111_vib-VIBRATIONS-1.mol -T 500 --floor 50
```

它会列出全部频率、标出并剔除虚频，`--floor 50` 把低于 50 cm⁻¹ 的实频抬高到 50 cm⁻¹，最后给出 ZPE、$U_\mathrm{vib}$、$S_\mathrm{vib}$、$TS_\mathrm{vib}$ 和 $F_\mathrm{vib}=U_\mathrm{vib}-TS_\mathrm{vib}$（eV 与 kJ/mol）。

吸附自由能的常见写法（以 CO 为例）：

$$
\Delta G_\mathrm{ads}(T) = \underbrace{E_\mathrm{CO^*} - E_\mathrm{slab} - E_\mathrm{CO(g)}}_{\Delta E_\mathrm{ads}}
+ \left[F_\mathrm{vib}^\mathrm{CO^*}(T)\right] - \left[G^\mathrm{CO(g)}(T) - E_\mathrm{CO(g)}\right]
$$

其中气相分子的 $G - E$ 可以用 §4 的 `THERMOCHEMISTRY` 输出（CO 的对称数为 1，正好不用修正），或者查热力学表。表面上的其他物种同样用简谐极限的 $F_\mathrm{vib}$。为了误差相消，各物种的可动原子选择应一致（只放开吸附物，或者都放开吸附物加表层原子），并且所有能量用同一套 `&FORCE_EVAL`。

### 5.5 选择可动原子：收敛性检查

部分 Hessian 是一种近似。检查办法：在可动集合中逐步加入与吸附物成键的衬底原子（如 CO 下方的 Pt 原子，再加其最近邻），比较 C–O 伸缩频率和 ZPE 的变化。对 ZPE 和 C–O 伸缩，通常只放开吸附物就已足够；低频的受阻模式对衬底原子更敏感，若要精确计算熵，应当加大可动区域。

### 5.6 `&MODE_SELECTIVE`：只跟踪感兴趣的模式

体系很大、又只关心一两个特征模式（例如某个 C=O 或 O–H 伸缩）时，可以用 mode-selective 振动分析：它用 Davidson 迭代只收敛与初猜相近的若干个模式，而不构造完整的 Hessian。下面的写法参考了官方回归测试 `tests/QS/regtest-ot-1-vib/H2O-VIB-MS-INT-1.inp`：

```fortran
&VIBRATIONAL_ANALYSIS
  DX 0.01
  NPROC_REP 4
  &MODE_SELECTIVE
    ! atoms displaced to build the initial guess vector
    ATOMS 28 29
    INITIAL_GUESS ATOMIC
    ! look for the mode closest to this wavenumber (cm^-1)
    FREQUENCY 2050
    EPS_MAX_VAL 1.0E-5
    EPS_NORM 1.0E-4
  &END MODE_SELECTIVE
&END VIBRATIONAL_ANALYSIS
```

| 关键词 | 说明 |
|---|---|
| `ATOMS` | 构造初猜时被位移的原子 |
| `INITIAL_GUESS` | `ATOMIC`（对 `ATOMS` 做随机位移，默认）、`BFGS_HESS`（取 BFGS 近似 Hessian）、`RESTART`/`RESTART_VEC`/`MOLDEN_RESTART`（从之前的结果继续） |
| `FREQUENCY` / `RANGE` / `&INVOLVED_ATOMS` | 三选一（互斥）：跟踪最接近某个波数的模式；跟踪某个波数区间内的模式；跟踪主要由某些原子贡献的模式（`&INVOLVED_ATOMS` 中可再设 `RANGE`，未设时只看 400 cm⁻¹ 以上） |
| `EPS_MAX_VAL`、`EPS_NORM` | Davidson 收敛判据（残差向量的最大分量和范数） |
| `LOWEST_FREQUENCY` | 写输出时包含的最低频率，负值可以把虚频也写出来 |
| `NPROC_REP` | 副本数大于 1 时使用块 Davidson，同时跟踪多个最接近的模式 |

Davidson 迭代的次数无法事先确定，对小体系不一定比全 Hessian 便宜；它的优势在于成本不随原子数线性增长。

---

## 6. 周期性晶体：Γ 点声子与 `FULLY_PERIODIC`

对体相晶体做 `VIBRATIONAL_ANALYSIS`，得到的是**当前晶胞的 Γ 点声子**。此时应设 `FULLY_PERIODIC T`：周期性晶体没有整体转动对称性，去掉“转动”会误删真实的模式；三个整体平动仍然是零频声学模式。要点：

- 晶胞必须先优化好（[第7章](07-geo-cell-opt.md)），原子受力要非常小；
- 一个 $n\times n\times n$ 超胞的 Γ 点频率等价于原胞布里渊区中对应 $q$ 点上的频率，但得不到完整的色散关系；
- 声子色散和态密度请用 phonopy 一类的工具：它支持 CP2K，在超胞中生成位移结构，用 CP2K 计算力，再构造力常数矩阵。CP2K 2026.2 的输入参考也注明，`THERMOCHEMISTRY` 只考虑 Γ 点（$q=0$），不支持一般 $q$ 点上的声子频率。

---

## 7. 过渡态验证与虚频处理

### 7.1 验证第13章的过渡态

`ts_freq_check.inp` 通过 `@INCLUDE '../13-neb-dimer-h-cu111/cu111_h_force_eval.inc'` 使用与 NEB/Dimer **完全相同**的 `&FORCE_EVAL`，并在包含之前用 `@SET EPS_SCF 1.0E-8` 把 SCF 收紧：

```fortran
@SET COORD ${COORD-ts_dimer.xyz}
@SET EPS_SCF 1.0E-8

&GLOBAL
  PROJECT ts_freq
  RUN_TYPE VIBRATIONAL_ANALYSIS
  PRINT_LEVEL MEDIUM
&END GLOBAL

@INCLUDE '../13-neb-dimer-h-cu111/cu111_h_force_eval.inc'

&VIBRATIONAL_ANALYSIS
  DX 0.01
  NPROC_REP ${NPROC_REP-8}
  &PRINT
    &MOLDEN_VIB ON
    &END MOLDEN_VIB
  &END PRINT
&END VIBRATIONAL_ANALYSIS

&MOTION
  &CONSTRAINT
    &FIXED_ATOMS
      LIST 1..${NFIX-18}
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

（`@INCLUDE` 的路径相对于运行 CP2K 时的当前目录。）默认固定底部两层（1–18），可动的是表层 9 个 Cu 和 H，共 30 个坐标、60 次力计算。`-E NFIX=27` 则只放开 H，3 个模式，适合快速检查。

```bash
cp ../13-neb-dimer-h-cu111/ts_dimer.xyz .
mpirun -np 64 cp2k.psmp -i ts_freq_check.inp -o ts_freq_check.out
grep "VIB|Frequency" ts_freq_check.out
python3 harmonic_thermo.py ts_freq_check.out     # 列出并剔除虚频，给出 TS 的 ZPE
```

判据：**恰好一个**负频率（虚频），用 molden 文件可视化时它把 H 从 fcc 推向 hcp。这个虚频的绝对值通常有几百 cm⁻¹（H 原子很轻）。

### 7.2 虚频的诊断与处理

| 现象 | 可能原因 | 处理 |
|---|---|---|
| 极小点出现 1 个 < 50i cm⁻¹ 的小虚频 | 数值噪声；软模式（甲基转动、吸附物的受阻平动） | 收紧 `EPS_SCF`、`EPS_DEFAULT`，提高 CUTOFF；改 `DX`（0.005/0.02）看是否稳定；检查优化是否真的收敛 |
| 极小点出现明显虚频（> 100i cm⁻¹） | 结构不是极小点（停在鞍点或山脊上，常见于对称性被人为保持的结构） | 沿该模式扰动后重新优化（下面的 `displace_mode.py`） |
| TS 有 2 个以上明显虚频 | 二阶鞍点 | 沿第二个虚频模式扰动后重新做 Dimer（[第13章 §8](13-neb-dimer.md)） |
| TS 没有虚频 | Dimer/CI-NEB 实际滑到了极小点，或可动原子选得不对 | 检查 Dimer 曲率和能量；确认虚频模式涉及的原子都在可动集合中 |
| 周期体系出现几十 cm⁻¹ 的“平动”频率 | 网格蛋盒效应、平动未能完全去除 | 提高 CUTOFF；`FULLY_PERIODIC` 设置是否合理 |

**沿模式扰动。** 示例目录中的 `displace_mode.py` 从 molden 文件中读出指定模式，把结构沿正、负方向各推一段距离（默认最大原子位移 0.1 Å），写出两个 xyz 文件：

```bash
python3 displace_mode.py h2o_opt.xyz h2o_vib-VIBRATIONS-1.mol --mode 1 --amp 0.10
# -> h2o_opt_mode1_plus.xyz, h2o_opt_mode1_minus.xyz ；分别重新 GEO_OPT，保留能量较低者
```

对过渡态，同样的操作就是 [第13章 §9](13-neb-dimer.md) 中的“伪 IRC”：两个扰动结构分别优化，应当落到初态和末态。

---

## 参数选择与收敛建议

| 参数 | 推荐 | 说明 |
|---|---|---|
| 优化 `MAX_FORCE` | 分子 1E-5，表面 1E-4 Ha/bohr | 频率计算前必须收紧 |
| `EPS_SCF` | 1E-7 至 1E-8 | 力噪声 → Hessian 噪声 |
| `EPS_DEFAULT` | 1E-12 至 1E-14 | 同上 |
| CUTOFF | 单点收敛值 + 100–200 Ry | 抑制蛋盒效应 |
| `DX` | 0.01 bohr；软模式多时与 0.005、0.02 对比 | 截断误差与噪声的折中 |
| `NPROC_REP` | 使副本组数 ≈ 可动坐标数（或其约数） | 位移结构并行计算 |
| `TC_TEMPERATURE` | 显式写 298.15 K（或目标温度） | 默认 273.15 K |
| 对称数/自旋 | 自行加 $RT\ln\sigma$、$-RT\ln g$ 修正 | CP2K 固定 $\sigma=1$、$g=1$ |
| 吸附物 | 固定衬底的部分 Hessian + 简谐极限后处理 | 有固定原子时不输出热化学 |
| 低频模式 | 抬高到 50 cm⁻¹，或 quasi-RRHO / 受阻平动模型 | RRHO 在低频处失效 |
| 晶体 | `FULLY_PERIODIC T`，色散用 phonopy | 仅 Γ 点 |

---

## 常见问题与报错

**1. `VIB| No IR intensities available. Check input`**
原因：设置了 `INTENSITIES T`，但 `&DFT/&PRINT` 中没有 `&MOMENTS`，CP2K 拿不到各位移结构的偶极矩。解决：加入 `&MOMENTS ON`（孤立体系 `PERIODIC F`，周期体系保持默认的 Berry 相）。

**2. 热化学数值明显不对，比如熵偏小、$G$ 校正偏大。**
原因：`TC_TEMPERATURE` 默认 273.15 K，不是 298.15 K。解决：显式设置 `TC_TEMPERATURE [K] 298.15`。另外别忘了对称数（§4.5）。

**3. 设了 `THERMOCHEMISTRY T`，输出中却没有 `THERMOCHEMICAL DATA`。**
原因：体系中有固定原子（部分 Hessian），CP2K 此时不计算热化学。解决：这是设计如此；用 `harmonic_thermo.py` 或 [vib_thermo.py](../scripts/vib_thermo.py) 在简谐极限下计算。

**4. 极小点出现几十 cm⁻¹ 的虚频。**
原因：几何没有严格收敛、SCF/网格噪声。解决：看 `Minimum Structure - Energy and Forces` 中的残余力；收紧优化和 SCF 判据；提高 CUTOFF；比较不同 `DX`；必要时沿该模式扰动后重新优化。

**5. 输出中“频率”的个数比预期多 5–6 个，且都接近 0 或几十 cm⁻¹。**
原因：有固定原子时不去除平动/转动（部分 Hessian 是这样设计的）；或孤立分子误设了 `FULLY_PERIODIC T`。解决：前者无需处理（这些就是受阻平动/转动）；后者删掉 `FULLY_PERIODIC`。

**6. 某个位移结构 SCF 不收敛，报 `SCF run NOT converged` 并中止。**
原因：金属体系的 SCF 设置偏激进，或参考结构本身就难收敛。解决：按 [第5章](05-scf-strategy.md) 调整混合与展宽，用优化得到的波函数重启（`SCF_GUESS RESTART`）。不要用 `IGNORE_CONVERGENCE_FAILURE`，否则 Hessian 中会混入未收敛的力。

**7. `Frc consts` 一行的数值与预期不符。**
原因：2026.1 之前的版本中力常数的计算有误（2026.1 发布说明中提到已修复）。解决：升级，或只用频率和约化质量。

**8. 振动分析计算量太大。**
原因：全 Hessian 需要 $6N+1$ 次力计算。解决：部分 Hessian（固定远离反应中心的原子）、`&MODE_SELECTIVE`、或者用 xTB 先做预筛选（[第16章](16-qmmm-xtb.md)）。

**9. 过渡态的 ZPE/热化学包含了虚频。**
原因：CP2K 的热化学把负频率也代入求和。解决：用 `harmonic_thermo.py` 等工具剔除虚频后重新计算。

---

## 练习

1. 对 H2O 分别用 `DX` = 0.002、0.005、0.01、0.02、0.05 bohr 计算频率，画出三个频率随 `DX` 的变化，找出既不受噪声也不受截断误差影响的区间。
2. 把 `EPS_SCF` 放宽到 1E-5、`EPS_DEFAULT` 放宽到 1E-8 重做 H2O，观察平动/转动“假频率”（`Cartesian Low frequencies`）和真实频率的变化。
3. 在 `h2o_force_eval.inc` 中启用 `MASS 2.014102`，计算 D2O 的频率，验证伸缩模式近似按 $1/\sqrt2$ 缩放；若使用 2026.2，再尝试用未质量加权 Hessian 离线重新加权。
4. 对 CO/Pt(111) 依次放开：只有 CO；CO + 下方 Pt；CO + 表层全部 Pt。比较 C–O 伸缩频率、ZPE 和 298 K 下的 $TS_\mathrm{vib}$。
5. 用 `THERMOCHEMISTRY T` 计算气相 CO 的 $G-E$，再与 `harmonic_thermo.py` 对 CO\* 的 $F_\mathrm{vib}$ 组合，估算 298 K、1 bar 下 CO 在 Pt(111) 顶位的吸附自由能。
6. 在 `ts_freq_check.inp` 中分别用 `NFIX=27` 和 `NFIX=18`，比较虚频和 TS 的 ZPE，评估部分 Hessian 的误差。

---

## 参考资料

**官方文档**

- 输入参考 `&VIBRATIONAL_ANALYSIS`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/VIBRATIONAL_ANALYSIS.html>（子段 `MODE_SELECTIVE`、`PRINT/MOLDEN_VIB`、`PRINT/HESSIAN` 在同目录下）
- 输入参考 `&DFT/&PRINT/&MOMENTS`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/MOMENTS.html>
- Geometry and Cell Optimization（优化收敛判据、输出文件）：<https://manual.cp2k.org/trunk/methods/optimization/geometry_and_cell_opt.html>
- CP2K 2026.2 release notes（振动输出修订、未质量加权 Hessian）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>；PR #5395：<https://github.com/cp2k/cp2k/pull/5395>
- CP2K 2026.1 release notes（振动分析力常数修复）：<https://github.com/cp2k/cp2k/releases/tag/v2026.1>
- Infrared Spectroscopy（官方页面目前只给出外部链接）：<https://manual.cp2k.org/trunk/methods/properties/infrared.html>
- NIST Chemistry WebBook，H2O：<https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185>；CO：<https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080>
- phonopy（含 CP2K 接口）：<https://phonopy.github.io/phonopy/>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-gpw-9/h2o_GRRM.inp`（IR 强度 + `&MOMENTS`）、`tests/QS/regtest-ot-1-vib/`（mode-selective）、`tests/QS/regtest-gpw-2-1/N2_raman.inp`（Raman）、`src/motion/vibrational_analysis.F`（热化学的具体实现）
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：中文 Wiki，含 CP2K 计算 ZPE 的实践
- [csc-training/sscc-cp2k](https://github.com/csc-training/sscc-cp2k)：NEB → 过渡态振动分析的完整练习
- [arvpinto/enzyme_neb_cp2k](https://github.com/arvpinto/enzyme_neb_cp2k)：QM/MM 过渡态的频率验证与分析脚本（CC0）
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含振动分析一章（2025.2 验证）

**本仓库脚本**：[vib_thermo.py](../scripts/vib_thermo.py)（解析频率与 IR 强度、检查虚频、ZPE 与 U/S/G）、[cp2k_outparse.py](../scripts/cp2k_outparse.py)（检查优化是否收敛）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations*（McGraw-Hill, 1955；Dover 重印），简正模式分析的经典教材
- D. A. McQuarrie, *Statistical Mechanics*（University Science Books, 2000），RRHO 配分函数
- S. Grimme, *Chem. Eur. J.* **18**, 9955 (2012)，quasi-RRHO 处理低频模式：<https://doi.org/10.1002/chem.201200497>
- L. H. Sprowl, C. T. Campbell, L. Árnadóttir, *J. Phys. Chem. C* **120**, 9719 (2016)，吸附物的受阻平动/受阻转动模型：<https://doi.org/10.1021/acs.jpcc.5b11616>
- J. K. Nørskov *et al.*, *J. Phys. Chem. B* **108**, 17886 (2004)，计算氢电极中 ZPE 与熵校正的典型用法：<https://doi.org/10.1021/jp047349j>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](13-neb-dimer.md) | [返回目录](../README.md) | [下一章 →](15-metadynamics-plumed.md)
