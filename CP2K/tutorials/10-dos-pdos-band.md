# 第10章 电子结构分析：DOS、PDOS 与能带（DOS, PDOS & Band Structure with CP2K）

本章讲如何用 CP2K Quickstep 计算并解读态密度（DOS）、投影态密度（PDOS）、原子列表局域 PDOS（LDOS）和能带结构（band structure）。重点有三项：一是 CP2K 2026.2 重构后的 `&DFT/&PRINT/&DOS` 新语法，包括 `&PDOS`、`&CURVE`、`&BROADEN`；二是 2026.1 及以前旧 `&PRINT/&PDOS` 写法的对照与迁移；三是 `&BAND_STRUCTURE` + `&KPOINT_SET` 的能带路径设置。我们用 H2O 分子、金红石 TiO2、石墨烯和体相 Si 四个例子走完“SCF → 输出文件 → 展宽 → 对齐能量零点 → 作图”的流程。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 按所用版本选对语法：2026.2 用 `&DOS` 套 `&PDOS`/`&CURVE`，2026.1 及以前用 `&PDOS`，并能把旧输入迁移到新版本；
- 说清 OT 与对角化两种 SCF 下“空轨道从哪里来”（`NLUMO` 与 `ADDED_MOS` 的分工）；
- 读懂 `.dos`、`.pdos`、`*_curve-*.pdos` 和 `.bs` 文件的每一列，知道哪些是原子单位、哪些是 eV、零点在哪里；
- 为周期体系设置 k 点 DOS（含 `MP_GRID` 加密）与能带路径，判断带隙类型；
- 正确选择能量零点（费米能、HOCO/VBM、真空能级），避免跨体系比较时出错；
- 用本仓库脚本或自己写的几十行 Python 完成展宽与作图。

| 想回答的问题 | 用什么 | 本章小节 |
|---|---|---|
| 体系是金属、半导体还是绝缘体，带隙多大 | 总 DOS + 能带 | 4、7、8 |
| 价带顶/导带底、HOMO/LUMO 来自哪种元素、哪个轨道 | `&PDOS`（`COMPONENTS T` 看 px/py/pz、d 分量） | 4、6 |
| 缺陷、掺杂、吸附物在带隙中有没有局域态 | `&LDOS LIST` / `&R_LDOS` | 9 |
| 能带色散、直接/间接带隙、有效质量的定性判断 | `&BAND_STRUCTURE` | 8 |
| 自旋极化体系（磁性氧化物、自由基）的上下自旋态 | ALPHA/BETA 两套文件 | 10 |

> [!NOTE]
> 本章默认你已经掌握 [第5章 SCF 策略](05-scf-strategy.md)（OT、对角化、展宽、`ADDED_MOS`）和 [第6章 周期性体系与 k 点](06-periodic-surface-kpoints.md)。示例输入集在 [../examples/10-dos-pdos-band/README.md](../examples/10-dos-pdos-band/README.md)。实空间性质（轨道 cube、电荷密度、功函数）见 [第11章](11-cube-charge-workfunction.md)。

---

## 2. 原理速览

### 2.1 DOS 与展宽

Kohn–Sham 本征值 $\varepsilon_{n\mathbf k}$ 给出的 DOS 是一串 δ 函数：

$$
g(E)=\sum_{n}\sum_{\mathbf k} w_{\mathbf k}\,\delta(E-\varepsilon_{n\mathbf k}),\qquad \sum_{\mathbf k} w_{\mathbf k}=1
$$

要画成连续曲线，需要把每个 δ 函数换成归一化的线型。CP2K 2026.2 的 `&BROADEN` 支持三种线型，`WIDTH` 一律是**半高全宽（FWHM）**：

| `TYPE` | 线型 | FWHM 与参数的关系 |
|---|---|---|
| `GAUSSIAN`（默认） | $\frac{1}{\sigma\sqrt{2\pi}}e^{-(E-\varepsilon)^2/2\sigma^2}$ | $\mathrm{FWHM}=2\sqrt{2\ln2}\,\sigma\approx2.3548\,\sigma$ |
| `LORENTZIAN` | $\frac{1}{\pi}\frac{\gamma}{(E-\varepsilon)^2+\gamma^2}$ | $\mathrm{FWHM}=2\gamma$ |
| `PSEUDO_VOIGT` | $\eta L+(1-\eta)G$ | 同一 FWHM；$\eta$ = `VOIGT_MIXING`（默认 0.5） |

展宽只是为了画图，**不能代替 k 点收敛**：稀疏 k 网格配大展宽，曲线很平滑，但物理上仍可能没收敛（官方 DOS 文档专门强调了这一点）。分子和大超胞只有 Γ 点，展宽曲线本质上是离散能级的“可视化”。

### 2.2 PDOS：CP2K 用的是 Löwdin 型投影

CP2K 的标准 PDOS 不是 Mulliken 分析，而是在原子轨道（AO）基下的 Löwdin 型投影：先用重叠矩阵的平方根把 MO 系数正交化，

$$
\widetilde{C}=S^{1/2}C,\qquad
p_{i}^{(A,l)}=\sum_{\mu\in A,\;l}\big|\widetilde{C}_{\mu i}\big|^2,\qquad
\sum_{A,l}p_{i}^{(A,l)}=1
$$

再把第 $i$ 个态的权重 $p_i$ 按“原子种类（KIND）× 角动量 $l$”累加。由于 $\widetilde C$ 列正交归一，每个态在所有通道上的权重之和为 1，且都是非负数。k 点计算对每个 $\mathbf k$ 做同样的投影，再按 $w_{\mathbf k}$ 加权。

由此有两点推论：

1. PDOS 依赖基组：弥散函数越多，“属于谁”越模糊。用不同基组得到的 PDOS 只能定性比较；
2. 对 PDOS 积分得到的“轨道布居”与 Mulliken、Löwdin、Hirshfeld 电荷都不同，定量讨论电荷请看 [第11章](11-cube-charge-workfunction.md)。

### 2.3 能带结构：在给定 k 路径上做非自洽对角化

CP2K 用局域高斯基，Kohn–Sham 矩阵先在实空间按晶格平移 $\mathbf R$ 存储，任意 $\mathbf k$ 点的哈密顿量由 Bloch 求和得到：

$$
H_{\mu\nu}(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}H_{\mu\nu}(\mathbf R),\qquad
S_{\mu\nu}(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}S_{\mu\nu}(\mathbf R)
$$

`&BAND_STRUCTURE` 在 SCF 收敛后，固定密度（即固定 $H(\mathbf R)$），对路径上每个 k 点求解 $H(\mathbf k)c=\varepsilon S(\mathbf k)c$。这一步是**非自洽**的，所以 SCF 用的 k 网格必须先收敛；能带路径本身不参与积分。

### 2.4 KS 带隙与能量零点

- GGA（PBE）的 KS 带隙系统性低估实验光学/基本带隙，常低估 30–50% 甚至更多；要定量比较可用杂化泛函（[第9章](09-hybrid-admm.md)）或 GW；
- KS 本征值的绝对数值依赖于静电势零点。周期体系中，Poisson 方程的 $\mathbf G=0$ 分量被设为零，即平均静电势为零，所以**不同晶胞之间不能直接比较本征值的绝对位置**。比较能带排列（band alignment）时，要借助真空能级或芯能级做静电势对齐（见第 11 节）。

---

## 3. 版本差异总览：2026.2 新语法与旧语法

CP2K 2026.2（2026-07-15）重构了 DOS/PDOS 输入（PR #5326）。这是一个**破坏性变更**：两套写法互不兼容。

| 功能 | 2026.1 及以前（旧） | 2026.2 及以后（新） |
|---|---|---|
| 总 DOS | `&PRINT/&DOS`（`DELTA_E`、`APPEND`、`MP_GRID`），只用 SCF 已有的态 | `&PRINT/&DOS`，新增 `NLUMO`、`NDIGITS`（默认 6） |
| 按 KIND 的 PDOS | `&PRINT/&PDOS`（`NLUMO`、`COMPONENTS`） | `&PRINT/&DOS/&PDOS`（只剩 `COMPONENTS`；`NLUMO` 移到 `&DOS`） |
| 原子列表 / 空间区域 | `&PDOS/&LDOS`、`&PDOS/&R_LDOS` | `&DOS/&LDOS`、`&DOS/&R_LDOS` |
| 展宽曲线 | 无，需自己展宽 | `&DOS/&CURVE`（`ENERGY_UNIT`、`ENERGY_ZERO`）+ `&CURVE/&BROADEN`（`TYPE`、`WIDTH`、`VOIGT_MIXING`） |
| k 点 PDOS | 不支持（输出警告 `Projected density of states (pDOS) is not implemented for k points`） | 支持，但只输出展宽曲线 |
| `.pdos` 文件列 | `MO  Eigenvalue[a.u.]  Occupation  s  p  d ...` | `MO  Energy[a.u.]  Occupation  Total  s  p  d ...`（多了 `Total` 列，表头多一行 `E(HOCO)`） |
| XAS / XAS_TDP 的 `&PDOS` | 原样 | 原样（不受此次重构影响） |

在 2026.2 上运行旧输入，会在解析阶段中止：

```text
unknown subsection PDOS of section PRINT
```

反过来，把新输入拿到 2026.1 及以前版本运行，会报 `&DOS` 段里不认识的子段或关键字（`CURVE`、`PDOS`、`NLUMO`）。

> [!IMPORTANT]
> **`DELTA_E` 没有单位。** 在 2026.2 的输入定义里（`src/input_cp2k_print_dft.F`），`&DOS/DELTA_E` 是没有关联单位的实数，按 Hartree（原子单位）解释，默认 `1.0E-3`（≈ 0.027 eV）。官方 DOS 方法页的示例写成了 `DELTA_E [eV] 0.01`；但根据 2026.2 的源码，这个关键字不接受单位标记，解析器会把 `[eV]` 当成数字读而失败（报错形如 `A floating point type object was expected, found <[eV]>`；这是我们读源码得出的结论，未在所有版本上实测）。请直接写数值，例如 `DELTA_E 0.0005`（≈ 0.0136 eV）。`&BROADEN/WIDTH` 则带单位（默认 eV），`WIDTH [eV] 0.1` 是合法的。

用下面的命令先检查输入是否与可执行文件版本匹配：

```bash
cp2k.psmp --version | head -n 3      # 看 CP2K version 2026.2 / 2025.2 ...
cp2k.psmp -c h2o_pdos_2026.inp       # 只做语法检查，不计算
```

---

## 4. 实战一：孤立 H2O 分子的 DOS/PDOS（2026.2 新语法，OT）

### 4.1 输入文件

`h2o_pdos_2026.inp`（完整文件见示例目录）：

```text
&GLOBAL
  PROJECT     h2o_pdos
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     400
      REL_CUTOFF 60
    &END MGRID
    &POISSON
      PERIODIC       NONE         ! must match &CELL PERIODIC NONE
      POISSON_SOLVER MT           ! box >= 2x the size of the charge distribution
    &END POISSON
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   30
      EPS_LUMO      1.0E-6        ! accuracy of the post-SCF virtual orbitals
      MAX_ITER_LUMO 400
      &OT
        MINIMIZER      DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &DOS
        NLUMO   10                ! virtual states included in DOS/PDOS
        DELTA_E 0.0005            ! grid spacing in HARTREE (no unit allowed), 0.0136 eV
        &PDOS
          COMPONENTS T            ! split p into py pz px, d into 5 components
        &END PDOS
        &CURVE
          ENERGY_UNIT EV          ! energy axis in eV, density in states/eV
          ENERGY_ZERO HOCO        ! HOMO at 0 eV (molecule, no smearing)
          &BROADEN
            TYPE  GAUSSIAN
            WIDTH [eV] 0.2        ! FWHM
          &END BROADEN
        &END CURVE
      &END DOS
      &MO
        ENERGIES           T
        OCCUPATION_NUMBERS T
        &EACH
          QS_SCF 0                ! only after the converged SCF
        &END EACH
      &END MO
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.0 12.0 12.0
      PERIODIC NONE
    &END CELL
    &COORD
      ! experimental gas-phase geometry, molecule in the yz plane
      O   0.000000   0.000000   0.117300
      H   0.000000   0.757000  -0.468600
      H   0.000000  -0.757000  -0.468600
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### 4.2 关键字逐项说明

| 位置 / 关键字 | 本例取值 | 说明 |
|---|---|---|
| `&PRINT/&DOS`（段本身） | 写出即开启 | DOS 打印键的默认级别是 `DEBUG`，**不写不输出**。写出该段即开启总 DOS（`.dos` 文件） |
| `&DOS/NLUMO` | 10 | 纳入 DOS/PDOS 的空轨道数，`-1` 表示全部。OT 计算在 SCF 后用 OT 本征求解器补算这些空轨道；对角化计算中，若 `ADDED_MOS` 不够，CP2K 会自动增加 `ADDED_MOS`（只增不减） |
| `&DOS/DELTA_E` | 0.0005 | 能量网格间距，**Hartree**。对离散输出是直方图的格宽，对曲线输出是采样间隔；应明显小于展宽宽度（经验上 ≤ FWHM/5） |
| `&DOS/&PDOS` | 写出即开启 | 每种 KIND 输出一个 `.pdos` 文件（自旋极化时每种 KIND×每个自旋一个） |
| `&PDOS/COMPONENTS` | T | 把 p、d 拆成各个实球谐分量，列名为 `s py pz px d-2 d-1 d0 d+1 d+2 ...`；F 时只给 `s p d` 之和 |
| `&DOS/&CURVE` | 写出即开启 | 额外输出可直接作图的展宽曲线文件；普通 `.dos`/`.pdos` 文件仍保留未展宽的数据 |
| `&CURVE/ENERGY_UNIT` | `EV` | 曲线能量轴的单位（`HARTREE` 默认 / `EV`），态密度随之换算成 states/eV |
| `&CURVE/ENERGY_ZERO` | `HOCO` | 曲线的能量零点：`AUTO`（默认，有展宽或分数占据时取 `FERMI`，否则取 `HOCO`）、`FERMI`、`HOCO`（最高占据晶体轨道）、`ABSOLUTE`（不平移）。**只作用于曲线文件** |
| `&BROADEN/TYPE`、`WIDTH` | `GAUSSIAN`、0.2 eV | 线型与 FWHM（默认 0.1 eV）。分子画图常用 0.1–0.3 eV |
| `&SCF/EPS_LUMO`、`MAX_ITER_LUMO` | 1E-6、400 | 控制 OT 补算空轨道的精度和最大迭代数（默认 1E-5、299） |
| `&PRINT/&MO` | `ENERGIES T` | 在主输出里列出本征值和占据数，便于核对 |
| `&POISSON PERIODIC NONE` + `MT` | — | 孤立分子必须与 `&CELL PERIODIC NONE` 一致，见第6章 |

### 4.3 运行与输出文件

```bash
cd CP2K/examples/10-dos-pdos-band
mpirun -np 4 cp2k.psmp -i h2o_pdos_2026.inp -o h2o_pdos_2026.out
ls -1 h2o_pdos*
```

主输出中与 DOS 相关的行：

```bash
grep -n "Calculate DOS at iteration step\|Calculate PDOS at iteration step\|additional unoccupied KS orbitals" h2o_pdos_2026.out
grep -n "MO| " h2o_pdos_2026.out | head -20      # eigenvalues printed by &MO
```

按 2026.2 源码（`qs_dos.F`、`qs_pdos.F`）的命名规则，会生成下列文件（项目名 `h2o_pdos`，单点计算的迭代后缀为 `-1`；请以 `ls` 结果为准）：

| 文件 | 内容 | 能量单位与零点 |
|---|---|---|
| `h2o_pdos-1.dos` | 未展宽总 DOS：`DELTA_E > 0` 时是直方图（只列非零格），列为 `Energy[a.u.]  Density  Occupation` | a.u.，绝对能量 |
| `h2o_pdos-curve-1.dos` | 展宽总 DOS 曲线，列为 `Energy[eV]  Density  Occupation`（自旋极化时为 `Alpha_Density Occupation Beta_Density Occupation`） | 按 `ENERGY_UNIT`/`ENERGY_ZERO` |
| `h2o_pdos-k1-1.pdos`、`h2o_pdos-k2-1.pdos` | 逐态投影权重（未展宽），每种 KIND 一个 | a.u.，绝对能量 |
| `h2o_pdos-k1_curve-1.pdos`、`h2o_pdos-k2_curve-1.pdos` | 该 KIND 的展宽 PDOS 曲线，第 2 列为 `total` | 按 `ENERGY_UNIT`/`ENERGY_ZERO` |

`kN` 中的 N 是原子种类编号，与输出里 `ATOMIC KIND INFORMATION` 的顺序一致；每个文件第一行都写明了元素，例如 `# Projected DOS for atomic kind O at iteration step i = ...`，不要靠猜。

> [!NOTE]
> 官方 DOS 文档的 gnuplot 示例里写的是 `project-k1-curve-1.pdos`，而 2026.2 源码中的中间名是 `k1_curve`（下划线）。不同小版本或补丁可能有出入，脚本里请用通配符（如 `*k2*curve*.pdos`）匹配，不要把文件名写死。

一个曲线文件的表头大致如下（示意输出，数值仅供参考）：

```text
# Projected DOS for atomic kind O at iteration step i = 1
# E(Fermi) =    -0.2xxxxx a.u. =    -7.xxxxx eV
# E(HOCO)  =    -0.2xxxxx a.u. =    -7.xxxxx eV
# Energy zero: HOCO
# Gaussian broadening, FWHM = 0.00734987 a.u. = 0.200000 eV
# Energy[eV]  total  s  py  pz  px  d-2  d-1  d0  d+1  d+2
```

### 4.4 如何判断结果是否合理

1. **SCF 收敛**：`grep "SCF run converged" h2o_pdos_2026.out`。2024.1 起不收敛会直接中止，不会产生 DOS 文件。
2. **占据态数目**：H2O 有 8 个价电子（O 6 + H 2×1），应有 4 个占据态，`h2o_pdos-k*-1.pdos` 中 `Occupation` 列为 2 的行恰好 4 行，后面跟着 10 个占据数为 0 的空轨道。
3. **轨道特征**：分子在 yz 平面内，HOMO（1b1）垂直于分子平面，在 O 的 PDOS 里应几乎全部落在 `px` 列；最低的价轨道（2a1）以 O 的 `s` 为主。这是检验 `COMPONENTS` 列顺序是否理解正确的好办法。
4. **HOMO–LUMO 间隙**：在 `.pdos` 文件中，取最后一个占据态和第一个空态的能量差（a.u. × 27.2114 = eV）。PBE 的 KS 能隙与实验光学能隙、电离能都不是一回事；例如 GGA 的 HOMO 本征值远浅于实验电离能（H2O 约 12.6 eV，NIST），这是 GGA 自相互作用误差的表现，并非计算出错。

快速作图（gnuplot 会忽略 `#` 开头的行）：

```bash
gnuplot -persist -e 'set xrange [-25:8]; set xlabel "E - E_{HOMO} (eV)";
  plot for [f in system("ls h2o_pdos-k*_curve-1.pdos")] f u 1:2 w l t f'
```

用本仓库脚本作图（按元素或角动量求和、分别画自旋通道等功能见帮助）：

```bash
python ../../scripts/pdos_plot.py -h
```

---

## 5. 空轨道从哪里来：OT 的 NLUMO 与对角化的 ADDED_MOS

DOS/PDOS 的上半部分（空态）是否完整，取决于 SCF 方法：

| SCF 方法 | 空轨道来源 | 2026.2 中该怎么写 | 2026.1 及以前 |
|---|---|---|---|
| OT（`&OT`） | SCF 中只有占据态；空轨道在 SCF 结束后用 OT 本征求解器单独计算 | `&DOS/NLUMO n`；精度由 `&SCF/EPS_LUMO`、`MAX_ITER_LUMO` 控制 | `&PDOS/NLUMO n`（作用相同） |
| 对角化（`&DIAGONALIZATION`） | 对角化直接给出 `ADDED_MOS` 个空轨道 | 在 `&SCF` 中设 `ADDED_MOS`；`NLUMO` 只作为下限，不足时 CP2K 自动加大 `ADDED_MOS` | `ADDED_MOS`；`NLUMO` 超过已有空轨道数时会额外补算 |
| k 点（必然是对角化） | 同上，每个 k 点都有 `ADDED_MOS` 个空带 | `ADDED_MOS`；能带计算另有 `&BAND_STRUCTURE/ADDED_MOS` | 同左 |

几点实用建议：

- 只关心价带和带边时，`NLUMO`/`ADDED_MOS` 取占据态数的 10–20% 就够；要看到导带以上数 eV，就要更多，可以先用小值试算看空态覆盖到多高的能量；
- `NLUMO -1`（全部空轨道）对大体系代价很高，旧版本的手册还提醒过：请求的态数超过基函数能提供的上限可能导致 Cholesky 分解出错；
- 金属体系用展宽时，`ADDED_MOS` 还必须覆盖费米能以上几个 $k_BT$，这是 SCF 正确性的要求，不只是为了画 DOS（见第5章）；
- OT 补算空轨道若 `MAX_ITER_LUMO` 用完仍未收敛，高处的空态能量不可靠，可以提高 `MAX_ITER_LUMO`，或改用对角化做一次单点（从 OT 收敛的 `RESTART.wfn` 起步，通常几步就收敛）。

---

## 6. 实战二：金红石 TiO2 超胞的 PDOS（2026.1 及以前的旧语法）

很多集群仍在用 2023.x–2025.x。旧语法里，PDOS 是 `&PRINT` 下的独立段，而且**不支持 k 点**，因此周期体系只能用 Γ 点超胞。

### 6.1 输入文件

`tio2_pdos_legacy_2025.inp` 的核心部分（完整文件见示例目录，文件末尾以注释形式给出了 2026.2 的等价写法）：

```text
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   200
      ADDED_MOS 100               ! virtual states come from the diagonalization
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.2
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC        ! always give METHOD explicitly (default changed between versions)
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
    &END SCF
    ...
    &PRINT
      ! LEGACY PDOS section (<= 2026.1)
      &PDOS
        NLUMO 100                 ! unoccupied states written to the .pdos files (-1 = all)
        COMPONENTS F              ! s, p, d sums per kind
        &LDOS
          LIST 1                  ! PDOS projected on atom 1 only (a Ti atom)
        &END LDOS
      &END PDOS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ! rutile, approximate experimental lattice parameters
      ABC 4.594 4.594 2.959
      MULTIPLE_UNIT_CELL 2 2 3
    &END CELL
    &TOPOLOGY
      MULTIPLE_UNIT_CELL 2 2 3    ! must match &CELL MULTIPLE_UNIT_CELL
    &END TOPOLOGY
    &COORD
      SCALED
      ! internal parameter u = 0.3048
      Ti  0.0000  0.0000  0.0000
      Ti  0.5000  0.5000  0.5000
      O   0.3048  0.3048  0.0000
      O   0.6952  0.6952  0.0000
      O   0.8048  0.1952  0.5000
      O   0.1952  0.8048  0.5000
    &END COORD
    &KIND Ti
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q12
      POTENTIAL GTH-PBE-q12
    &END KIND
```

| 关键字 | 说明 |
|---|---|
| `MULTIPLE_UNIT_CELL 2 2 3`（`&CELL` 与 `&TOPOLOGY` 各写一次） | 由 6 原子原胞复制出 72 原子超胞（9.19 × 9.19 × 8.88 Å）。两处必须一致，否则报晶胞或原子数错误 |
| `ADDED_MOS 100` | 对角化给出 100 个空轨道。体系有 576 个价电子（Ti 24×12 + O 48×6），288 个占据态 |
| `&SMEAR ... METHOD FERMI_DIRAC` | TiO2 是绝缘体，300 K 的 Fermi–Dirac 展宽对占据数几乎没有影响，只是让对角化 SCF 更稳。**务必写 `METHOD`**：2026.1 及以前的默认值是 `ENERGY_WINDOW`，2026.2 改为 `GAUSSIAN` |
| `&PDOS/NLUMO 100` | 写进 `.pdos` 的空轨道数（`-1` = 全部）。对角化时若超过 `ADDED_MOS`，会再用迭代求解器补算 |
| `&LDOS/LIST 1` | 额外输出投影到第 1 号原子（一个 Ti）的 PDOS，文件名含 `list1`。复制超胞后原子编号以输出中的坐标表为准 |
| `GTH-PBE-q12` + `...-q12` 基组 | Ti 含 3s3p 半芯态，基组与赝势的 `-q12` 必须一致 |

### 6.2 旧格式 `.pdos` 文件

旧版本的 `tio2_pdos-k1-1.pdos` 形如（示意输出，数值仅供参考）：

```text
# Projected DOS for atomic kind Ti at iteration step i = 0, E(Fermi) =    -0.0xxxxx a.u.
#     MO Eigenvalue [a.u.]      Occupation                 s                 p                 d
       1       -2.xxxxxx        2.000000        0.xxxxxxxx        0.xxxxxxxx        0.xxxxxxxx
```

- 第 1 列：MO 序号；第 2 列：本征值（Hartree，绝对值）；第 3 列：占据数（闭壳层最大为 2）；之后是各角动量通道的权重；
- `E(Fermi)` 写在第一行末尾；2026.2 的 `.pdos` 把 `E(Fermi)`、`E(HOCO)` 各放一行，并在 `Occupation` 后多了 `Total` 列。**自己写的解析脚本要按表头的列名定位，不要按固定列号取数**。

### 6.3 自己做展宽：一个兼容新旧格式的 Python 片段

旧格式只有离散态，需要自己展宽。下面这段代码按表头识别列名，因此对 2026.2 的非曲线 `.pdos` 同样适用（依赖 numpy、matplotlib）：

```python
import glob, re
import numpy as np
import matplotlib.pyplot as plt

HA2EV = 27.211386

def read_pdos(fname):
    """Return element, E_F (eV), energies (eV), {channel: weights} of a .pdos file."""
    with open(fname) as f:
        head = [f.readline() for _ in range(4)]
    elem = re.search(r"atomic kind\s+(\S+)", head[0]).group(1)
    ef = float(re.search(r"E\(Fermi\)\s*=\s*(-?\d+\.\d+)", "".join(head)).group(1)) * HA2EV
    labels = next(h for h in head if "Occupation" in h).split()
    labels = labels[labels.index("Occupation") + 1:]          # s p d ... (or Total s p d ...)
    data = np.loadtxt(fname)                                   # '#' lines are skipped
    ener = data[:, 1] * HA2EV
    chans = {lab: data[:, 3 + i] for i, lab in enumerate(labels) if lab.lower() != "total"}
    return elem, ef, ener, chans

def gauss(grid, ener, weights, fwhm=0.2):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    x = (grid[:, None] - ener[None, :]) / sigma
    return np.exp(-0.5 * x * x) @ weights / (sigma * np.sqrt(2.0 * np.pi))

grid = np.arange(-20.0, 10.0, 0.01)
for fn in sorted(glob.glob("*-k[0-9]*-1.pdos")):
    if "curve" in fn:
        continue
    elem, ef, ener, chans = read_pdos(fn)
    for lab, w in chans.items():
        plt.plot(grid, gauss(grid, ener - ef, w), label=f"{elem} {lab}")
plt.axvline(0.0, ls="--", c="k", lw=0.8)
plt.xlabel("E - E$_F$ (eV)"); plt.ylabel("PDOS (states/eV)"); plt.legend()
plt.savefig("pdos.png", dpi=200)
```

对绝缘体，把零点换成价带顶（最后一个占据态的能量）通常更直观；E(Fermi) 在未展宽或展宽很小时对绝缘体的定义并不统一。

### 6.4 迁移到 2026.2

把整个 `&PDOS ... &END PDOS` 换成：

```text
    &PRINT
      &DOS
        NLUMO   100
        DELTA_E 0.0005
        &PDOS
          COMPONENTS F
        &END PDOS
        &LDOS
          LIST 1
        &END LDOS
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO AUTO
          &BROADEN
            WIDTH [eV] 0.1
          &END BROADEN
        &END CURVE
      &END DOS
    &END PRINT
```

要点：`NLUMO` 从 `&PDOS` 移到 `&DOS`；`&LDOS`/`&R_LDOS` 从 `&PDOS` 下移到 `&DOS` 下，与 `&PDOS` 并列；需要展宽曲线时再加 `&CURVE`。2026.2 起金红石也可以直接用原胞加 k 点算 PDOS（见第 7 节），不必再用超胞。

---

## 7. 实战三：石墨烯的 k 点 DOS 与 PDOS（2026.2）

### 7.1 输入要点

`graphene_dos_kpoints.inp` 的关键部分：

```text
    &KPOINTS
      ! one k-point along the vacuum direction
      SCHEME MONKHORST-PACK 21 21 1
    &END KPOINTS
    &SCF
      ADDED_MOS 10
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.3
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
    &END SCF
    ...
    &PRINT
      &DOS
        DELTA_E 0.0002            ! Hartree, 0.0054 eV
        MP_GRID 60 60 1           ! dense mesh for the TOTAL DOS only
        &PDOS
          COMPONENTS T            ! columns: s py pz px (d components follow)
        &END PDOS
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO FERMI       ! metal / semimetal with smearing
          &BROADEN
            TYPE  GAUSSIAN
            WIDTH [eV] 0.1
          &END BROADEN
        &END CURVE
      &END DOS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 2.46 2.46 20.0
      ALPHA_BETA_GAMMA 90.0 90.0 120.0
    &END CELL
    &COORD
      SCALED
      C  0.333333333  0.666666667  0.500000000
      C  0.666666667  0.333333333  0.500000000
    &END COORD
```

| 设置 | 为什么这样选 |
|---|---|
| `MONKHORST-PACK 21 21 1` | 二维体系，真空方向只取 1 个 k 点。对 Monkhorst–Pack 网格，奇数 n 包含 Γ；n 又能被 3 整除时，K 点 $(1/3,1/3)$ 恰在网格上，能抓住石墨烯的狄拉克点 |
| 对角化 + `&SMEAR METHOD FERMI_DIRAC` | 2026.2 中 k 点计算只能用对角化（OT+k 点在 trunk 才支持）；半金属需要展宽稳定占据数 |
| `&DOS/MP_GRID 60 60 1` | 只对**总 DOS** 生效：SCF 收敛后，在更密的 60×60×1 网格上非自洽地重新求本征值来统计 DOS，代价很低。PDOS 仍用 SCF 的 21×21×1 网格 |
| `ENERGY_ZERO FERMI` | 有展宽时 `AUTO` 也会选 `FERMI`，这里显式写出 |
| `COMPONENTS T` | 区分 σ（s、px、py）与 π（pz）态 |

### 7.2 k 点 PDOS 的限制（2026.2）

根据 2026.2 源码和官方文档：

- k 点下**没有逐态的 `.pdos` 文件**，只输出展宽曲线（每个 KIND 一个 `*_curve-*.pdos`），因此 **k 点 PDOS 必须配 `&CURVE`**，否则什么也不写。输出里会出现警告 `State-resolved k-point PDOS output is not implemented yet`，属于正常现象；
- `&LDOS`/`&R_LDOS` 在 k 点下不实现，会被忽略并给出警告；
- 打开 k 点原子对称约化（`&KPOINTS SYMMETRY T`，2026.2 中仍是实验性功能）时，约化后的 k 点会“旋转”轨道分量，`COMPONENTS T` 的 px/py/pz 分量不可靠。官方建议分析分量时关闭对称性。默认的 k/−k 时间反演配对不受影响；
- 总 DOS 用 `MP_GRID` 加密时，各 KIND 的 PDOS 之和与总 DOS 来自不同的网格，两者不会严格相等。只有在同一网格、同一展宽下，PDOS 的 `total` 列对所有 KIND 求和才等于总 DOS。

### 7.3 结果检查

```bash
grep -n "Fermi energy:" graphene_dos_kpoints.out | tail -n 1     # Hartree
grep -n "Using a" graphene_dos_kpoints.out                         # k-point grid used for the DOS
head -n 8 graphene_dos-curve-1.dos
```

- 总 DOS 在 0 eV（费米能）附近应呈“V”形、在狄拉克点处趋于零（展宽和有限网格会让底部有一点残余）；
- 在 C 的 PDOS 中，−2 到 +2 eV 范围几乎全部是 `pz`；σ 带（s、px、py）在费米能以下约 3 eV 以外才出现；
- 把 `MP_GRID` 从 60 改到 90，若曲线形状基本不变，说明总 DOS 已收敛到展宽宽度的分辨率。

> [!TIP]
> 二维材料、表面 slab 的 DOS 可以完全照此设置：面内 k 点要密，真空方向取 1。能量零点若要与其他体系对齐，请用真空能级而不是费米能，见 [第11章 功函数](11-cube-charge-workfunction.md)。

---

## 8. 实战四：体相 Si 的能带结构（&BAND_STRUCTURE）

### 8.1 输入文件

`si_band_structure.inp` 的核心部分：

```text
    &KPOINTS
      SCHEME MONKHORST-PACK 8 8 8
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-8
      MAX_SCF   100
      ADDED_MOS 8                 ! 4 occupied + 8 empty bands in the SCF
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.4
        NBUFFER 8
      &END MIXING
    &END SCF
    ...
    &PRINT
      &BAND_STRUCTURE
        FILE_NAME si_band.bs
        ADDED_MOS 8               ! extra bands for the band-structure step
        &KPOINT_SET
          ! coordinates in fractions of the reciprocal lattice vectors
          UNITS   B_VECTOR
          NPOINTS 30              ! points per segment
          SPECIAL_POINT L      0.500 0.500 0.500
          SPECIAL_POINT GAMMA  0.000 0.000 0.000
          SPECIAL_POINT X      0.500 0.000 0.500
          SPECIAL_POINT U      0.625 0.250 0.625
        &END KPOINT_SET
        &KPOINT_SET
          UNITS   B_VECTOR
          NPOINTS 30
          SPECIAL_POINT K      0.375 0.375 0.750
          SPECIAL_POINT GAMMA  0.000 0.000 0.000
        &END KPOINT_SET
      &END BAND_STRUCTURE
      &DOS
        DELTA_E 0.0005            ! Hartree
        &PDOS
          COMPONENTS F
        &END PDOS
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO HOCO        ! VBM at 0 eV, same zero as the band plot
          &BROADEN
            WIDTH [eV] 0.1
          &END BROADEN
        &END CURVE
      &END DOS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ! primitive fcc vectors, a = 5.431 Angstrom (experimental)
      A  0.00000  2.71550  2.71550
      B  2.71550  0.00000  2.71550
      C  2.71550  2.71550  0.00000
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Si  0.00  0.00  0.00
      Si  0.25  0.25  0.25
    &END COORD
```

| 关键字 | 说明 |
|---|---|
| `&KPOINTS SCHEME MONKHORST-PACK 8 8 8` | **自洽**积分用的网格，先按第6章收敛。能带路径与它无关 |
| `&BAND_STRUCTURE/FILE_NAME` | 能带数据文件名，按原样使用（不加项目名前缀）。**以追加方式打开**，重复运行前要删除旧文件；不写则数据打印到主输出 |
| `&BAND_STRUCTURE/ADDED_MOS`（别名 `ADDED_BANDS`） | 在能带计算中额外增加的带数，与 `&SCF/ADDED_MOS` 是两回事 |
| `&KPOINT_SET`（可重复） | 一段连续路径。多个 `&KPOINT_SET` 之间不连续，可以用来表示路径中的“跳跃”（如 U\|K） |
| `UNITS` | `B_VECTOR`（默认，倒格矢分数坐标）、`CART_ANGSTROM`（2π/Å）、`CART_BOHR`（2π/bohr）。推荐 `B_VECTOR`：它与晶胞在空间中的朝向无关 |
| `SPECIAL_POINT 名称 kx ky kz` | 名称可省略（只写 3 个数）；名称会写进 `.bs` 文件，方便作图标注 |
| `NPOINTS` | 每两个相邻特殊点之间的分段数；一个 `KPOINT_SET` 的总点数为 $(N_\mathrm{special}-1)\times\mathrm{NPOINTS}+1$，本例两段分别为 91 和 31 个 k 点 |

### 8.2 高对称点从哪里来

高对称点坐标**依赖于你选的原胞矢量**。本例用的 fcc 原胞矢量 $\mathbf a_1=(0,\tfrac a2,\tfrac a2)$、$\mathbf a_2=(\tfrac a2,0,\tfrac a2)$、$\mathbf a_3=(\tfrac a2,\tfrac a2,0)$ 与 Setyawan–Curtarolo（2010）的约定一致，对应的倒格矢分数坐标为：

| 点 | 坐标（B_VECTOR） | 点 | 坐标（B_VECTOR） |
|---|---|---|---|
| Γ | (0, 0, 0) | L | (1/2, 1/2, 1/2) |
| X | (1/2, 0, 1/2) | W | (1/2, 1/4, 3/4) |
| K | (3/8, 3/8, 3/4) | U | (5/8, 1/4, 5/8) |

换了晶胞（例如用常规立方胞，或让程序重新取向）时，这张表就不再适用。可靠的做法是用 SeeK-path（Hinuma 等，2017）或 pymatgen 的 `HighSymmKpath` 生成路径（见 Common 第4章），再按 CP2K 的 `SPECIAL_POINT` 格式写入；注意这些工具可能会先把结构标准化成它们约定的原胞，此时要把**标准化后的晶胞**写回 CP2K 输入，或者把 k 点坐标转换回你的晶胞。

> [!TIP]
> 想确认 CP2K 理解的 k 点与你预期一致，看主输出中的 `KPOINTS| Special point` 行，它以 `KPOINTS| In units of b-vector [2pi/Bohr]` 的形式回显了每个特殊点。

### 8.3 运行与 `.bs` 文件格式

```bash
rm -f si_band.bs
mpirun -np 4 cp2k.psmp -i si_band_structure.inp -o si_band_structure.out
grep "KPOINTS|" si_band_structure.out
```

`si_band.bs` 的结构（按 2026.2 源码 `qs_band_structure.F` 的写法；示意输出，数值仅供参考）：

```text
# Set 1: 4 special points, 91 k-points, 20 bands
#  Special point 1      0.50000000     0.50000000     0.50000000  L
#  Special point 2      0.00000000     0.00000000     0.00000000  GAMMA
...
#  Point 1  Spin 1:     0.50000000     0.50000000     0.50000000      0.01098901
#   Band    Energy [eV]     Occupation
       1     -x.xxxxxxxx     2.00000000
       2     -x.xxxxxxxx     2.00000000
...
# Set 2: 2 special points, 31 k-points, 20 bands
```

- 每个 `# Set` 对应一个 `&KPOINT_SET`；`Point` 行给出 k 点的倒格矢分数坐标和权重；
- 能量单位是 **eV**，为**绝对值**（未减去费米能或 VBM），作图时要自己平移；
- 自旋极化计算中，每个 k 点会先后出现 `Spin 1` 和 `Spin 2` 两块。

### 8.4 作图与带隙判断

推荐用本仓库脚本，它会解析 Set/Point/Band 块、生成高对称点标签，并导出 Origin 可用的 CSV：

```bash
python ../../scripts/band_plot.py -h
```

也可以用官方 cp2k-output-tools 中的 `cp2k_bs2csv` 把 `.bs` 转成 CSV。若想自己处理，下面这段代码展示了格式要点（只处理非自旋极化，依赖 numpy、matplotlib）：

```python
import numpy as np
import matplotlib.pyplot as plt

a = 5.431
cell = np.array([[0, a/2, a/2], [a/2, 0, a/2], [a/2, a/2, 0]])     # rows = A, B, C
bmat = 2 * np.pi * np.linalg.inv(cell).T                              # rows = b1, b2, b3

sets, cur = [], None
for line in open("si_band.bs"):
    t = line.replace(":", " ").split()
    if line.startswith("# Set"):
        cur = {"labels": [], "k": [], "e": []}
        sets.append(cur)
    elif line.startswith("#  Special point"):
        cur["labels"].append(" ".join(t[7:]) or "?")
    elif line.startswith("#  Point"):
        cur["k"].append([float(x) for x in t[5:8]])
        cur["e"].append([])
    elif not line.startswith("#") and t:
        cur["e"][-1].append(float(t[1]))

vbm = max(max(ek[:4]) for s in sets for ek in s["e"])     # Si: 4 occupied bands per k-point
x0, ticks = 0.0, []
for s in sets:
    k = np.array(s["k"]) @ bmat
    x = x0 + np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(k, axis=0), axis=1))])
    e = np.array(s["e"]) - vbm
    plt.plot(x, e, "b-", lw=1)
    n = (len(x) - 1) // (len(s["labels"]) - 1)
    ticks += [(x[i * n], lab) for i, lab in enumerate(s["labels"])]
    x0 = x[-1]
plt.xticks([p for p, _ in ticks], [l.replace("GAMMA", r"$\Gamma$") for _, l in ticks])
plt.axhline(0, ls="--", c="k", lw=0.8); plt.ylim(-13, 6); plt.ylabel("E - E$_{VBM}$ (eV)")
plt.savefig("si_band.png", dpi=200)
```

判断要点：

- Si 是间接带隙：VBM 在 Γ，CBM 在 Γ–X 线上靠近 X 处。能带路径上 CBM 的位置取决于采样密度，`NPOINTS` 太小会高估带隙；
- 实验带隙约 1.17 eV（低温），PBE 的 KS 带隙明显偏小（平面波 PBE 文献值约 0.5–0.6 eV）。CP2K 中的数值还受基组影响：DZVP 对导带的描述不如 TZV2P 完备，导带可能偏高；定量比较时请换 `TZV2P-MOLOPT-PBE-GTH-q4` 复查；
- 同一个输入里的 `&DOS` 用 `ENERGY_ZERO HOCO`，曲线零点就是 VBM，可与能带图并排画在同一纵轴上。

> [!NOTE]
> 能带只对 k 点自洽计算有意义：先有收敛的 `&KPOINTS` SCF，再在路径上做非自洽对角化（官方 K-Points 文档的建议）。`&KPOINTS SCHEME GENERAL` 是给 SCF 积分用的加权 k 点集合，不是画能带路径的接口。杂化泛函的能带要用 RI-HFXk，见 [第9章](09-hybrid-admm.md)。

---

## 9. LDOS 与 R_LDOS：投影到指定原子或空间区域

CP2K 中的 “LDOS” 是**指定原子列表上的 PDOS**，不是 STM 分析中常说的实空间局域态密度 $\rho(\mathbf r,E)$（官方文档特别提醒过这一区别）。

```text
      &DOS
        NLUMO 50
        DELTA_E 0.0005
        &LDOS
          LIST 1 2 3            ! adsorbate atoms
          COMPONENTS F
        &END LDOS
        &LDOS
          LIST 25 26 27 28      ! surface atoms below the adsorbate
        &END LDOS
        &R_LDOS
          LIST 1                ! volume around atom 1
          XRANGE -1.5 1.5       ! Angstrom, relative to the listed atoms
          YRANGE -1.5 1.5
          ZRANGE  0.0 3.0
          ERANGE -0.5 0.1       ! Hartree, only project states in this window
        &END R_LDOS
        &CURVE
          ENERGY_UNIT EV
          &BROADEN
            WIDTH [eV] 0.1
          &END BROADEN
        &END CURVE
      &END DOS
```

| 关键字 | 单位 | 说明 |
|---|---|---|
| `&LDOS/LIST` | — | 原子序号（从 1 开始），可写多行累加；每个 `&LDOS` 段输出一个文件（`...-list1-1.pdos`、`...-list2-1.pdos`） |
| `&LDOS/COMPONENTS` | — | 同 `&PDOS` |
| `&R_LDOS/LIST` + `XRANGE/YRANGE/ZRANGE` | Å | 以列表中原子为参照的空间范围：距列表中至少一个原子在该范围内的网格点都计入 |
| `&R_LDOS/ERANGE` | Hartree | 只投影能量在该区间内的态（默认全部） |

典型用途：判断缺陷能级是否局域在缺陷原子上；比较吸附物与表面原子在同一能量窗口内的态（杂化程度）；对 slab 区分表层与“体相层”。LDOS/R_LDOS **不支持 k 点**，需要用 Γ 点大超胞。在 2026.1 及以前，这两个段写在 `&PRINT/&PDOS` 下面。

---

## 10. 自旋极化与 DFT+U 体系的 DOS

开启 `UKS` 后（见 [第8章](08-dispersion-spin-dftu.md)）：

- PDOS 文件按自旋分开：`project-ALPHA_k1-1.pdos`、`project-BETA_k1-1.pdos`，曲线文件同理（`ALPHA_k1_curve`、`BETA_k1_curve`）；
- 总 DOS 曲线文件中同时给出 `Alpha_Density Occupation Beta_Density Occupation` 四列；
- 两个自旋通道使用**同一个**能量零点（取两者中较高的费米能或 HOCO），所以上下自旋曲线可以直接画在同一能量轴上；习惯上把 β 通道乘以 −1 画在下方；
- 反铁磁体系（如 NiO）要用两个 KIND（`Ni1`、`Ni2`，`ELEMENT Ni`）设置相反的初始磁矩。这样 PDOS 也按 `Ni1`/`Ni2` 分开输出，正好用来检查两套子晶格的自旋是否相反；
- DFT+U 会显著改变 d 态的位置，比较 PBE 与 PBE+U 的 PDOS 是判断 U 值影响最直观的方法。2026.2 起 Mulliken 型 DFT+U 支持 k 点，原胞加 k 点即可算 PDOS。

---

## 11. 能量零点与跨体系对齐

| 零点 | 适用 | CP2K 中怎么得到 |
|---|---|---|
| 费米能 $E_F$ | 金属、带展宽的计算 | 主输出 `Fermi energy:`（Hartree，仅在用展宽时打印）；曲线文件表头 `E(Fermi)` |
| HOCO / VBM / HOMO | 分子、半导体、绝缘体 | `ENERGY_ZERO HOCO`；`.pdos` 中最后一个占据态 |
| 真空能级 $E_\mathrm{vac}$ | 比较不同表面、分子与表面的能级排列、功函数 | `V_HARTREE_CUBE` 平面平均，见 [第11章](11-cube-charge-workfunction.md) |
| 芯能级或平均静电势 | 体相异质结的能带偏移 | 需单独的界面计算做静电势对齐 |

`ENERGY_ZERO` 只是给曲线文件统一平移，**不是**静电势对齐。官方 DOS 文档也明确说明：比较不同晶胞或表面的绝对能级时，仅改变 DOS 零点不能代替静电势对齐。对同一结构的不同电荷态、不同自旋态，至少要使用同一种零点约定。

---

## 12. 后处理工具一览

| 工具 | 作用 | 备注 |
|---|---|---|
| `CP2K/scripts/pdos_plot.py`（本仓库） | 读旧版 `.pdos` 与 2026.2 曲线文件，高斯/洛伦兹展宽，按元素或角动量求和，分自旋作图 | `python ../scripts/pdos_plot.py -h` 查看用法 |
| `CP2K/scripts/band_plot.py`（本仓库） | 解析 `.bs`，生成横轴与高对称点标签，以 VBM 或费米能为零点，导出 CSV | `python ../scripts/band_plot.py -h` |
| gnuplot | 直接画 2026.2 的曲线文件（`#` 行自动忽略） | 官方 DOS 文档给出了示例 |
| [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools) | `cp2k_bs2csv`（能带转 CSV）、`cp2k_pdos`（对旧 `.pdos` 做展宽） | 最新 tag 仍是 v0.4.0，新输出格式请先小样本验证 |
| [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata) | 读 pdos、cube 和输出文件的 Python 库 | 2026.2 的 `.pdos` 多了 `Total` 列，使用前先确认解析正确 |
| Multiwfn | 读 CP2K 的 molden 文件画 TDOS/PDOS/OPDOS，做轨道成分分析 | molden 的生成见 [第11章](11-cube-charge-workfunction.md)；k 点计算无法输出 molden |

通用的出版级作图（配色、字体、能带与 DOS 并排）见 Common 第5章 后处理与可视化。

---

## 参数选择与收敛建议

| 参数 | 推荐起点 | 判据 / 说明 |
|---|---|---|
| SCF k 网格（DOS 用） | 比总能收敛所需更密，约 1.5–2 倍 | 带边与费米能附近的 DOS 形状不随网格加密而变化 |
| `&DOS/MP_GRID` | SCF 网格的 2–4 倍 | 仅总 DOS；k 点计算专用 |
| `DELTA_E`（Hartree） | 0.0002–0.0005 | 采样间隔 ≤ FWHM/5；不要写单位 |
| `&BROADEN/WIDTH` | 0.05–0.1 eV（固体），0.1–0.3 eV（分子） | 小于要分辨的特征（带隙、自旋劈裂、缺陷峰） |
| `NLUMO` / `ADDED_MOS` | 占据态数的 10–30% | 空态覆盖到你关心的最高能量；金属展宽时还要覆盖 $E_F$ 以上几个 $k_BT$ |
| `&BAND_STRUCTURE/NPOINTS` | 20–50 / 段 | 带边位置稳定；CBM 不在高对称点时要更密 |
| 基组 | DZVP 看价带；TZV2P 复查导带和带隙 | 空态对基组更敏感 |
| `EPS_SCF` | 1E-7（OT）～1E-8（对角化） | 非自洽能带依赖收敛的密度 |
| `&KPOINTS SYMMETRY` | 默认关闭 | 分析 `COMPONENTS` 时保持关闭 |

---

## 常见问题与报错

**1. 2026.2 运行旧输入报 `unknown subsection PDOS of section PRINT`。**
原因：2026.2 重构了 DOS/PDOS 输入。解决：按第 6.4 节把 `&PRINT/&PDOS` 改为 `&PRINT/&DOS/&PDOS`，`NLUMO` 移到 `&DOS`，`&LDOS` 移到 `&DOS` 下。

**2. 旧版本运行新输入报 `&DOS` 段中的未知关键字或子段（如 `NLUMO`、`CURVE`、`PDOS`）。**
原因：新语法需要 CP2K ≥ 2026.2。解决：升级，或改回旧语法；输入开头用注释写明适用版本。

**3. 写了 `DELTA_E [eV] 0.01`，解析报错（形如 `A floating point type object was expected ...`）。**
原因：2026.2 中 `DELTA_E` 没有定义单位。解决：直接写 Hartree 数值，如 `DELTA_E 0.0005`。

**4. 有 `.dos` 文件但没有任何 `*_curve*` 文件，输出里有警告 `Broadened DOS output requires DELTA_E > 0 and will be skipped`。**
原因：`DELTA_E` ≤ 0（或 `WIDTH` ≤ 0）。解决：设一个正的 `DELTA_E`。

**5. k 点计算只得到总 DOS，没有 PDOS，输出警告 `State-resolved k-point PDOS output is not implemented yet`。**
原因：k 点 PDOS 只以展宽曲线形式输出。解决：在 `&DOS` 中加入 `&CURVE`（含 `&BROADEN`）。若同时设了 `&LDOS`，还会看到 `LDOS/R_LDOS are not implemented for k-point PDOS and will be ignored`，需要改用 Γ 点超胞。

**6. 旧版本 k 点计算报警告 `Projected density of states (pDOS) is not implemented for k points`，没有 PDOS 文件。**
原因：2026.1 及以前不支持 k 点 PDOS。解决：升级到 2026.2，或用 Γ 点超胞。

**7. PDOS 的空态部分是空的，或只到 LUMO 附近。**
原因：OT 计算没设 `NLUMO`，或对角化的 `ADDED_MOS` 太小。解决：见第 5 节；OT 下检查补算空轨道是否在 `MAX_ITER_LUMO` 内收敛。

**8. `.bs` 文件里出现了两倍、三倍的 `# Set` 块。**
原因：`&BAND_STRUCTURE` 以追加方式写文件，旧文件没删。解决：每次运行前 `rm -f xxx.bs`（本示例的 `job_slurm.sh` 会自动删除）。

**9. 能带图的带隙比预期大，CBM 不在应在的位置。**
原因：路径采样太稀，或高对称点坐标与晶胞约定不一致。解决：增大 `NPOINTS`；核对输出中的 `KPOINTS| Special point`；用 SeeK-path 重新生成路径。

**10. PDOS 中 `px`/`py`/`pz` 的相对大小与对称性不符。**
原因：打开了 k 点对称约化，或者体系在晶胞中的朝向与你设想的不同。解决：关闭 `SYMMETRY`；必要时检查 CP2K 是否对晶胞做了规范化（`&CELL/CANONICALIZE`）。

**11. 2026.2 中 `&SMEAR` 没写 `METHOD`，结果和老版本不一样。**
原因：2026.2 新增了 Gaussian、Methfessel–Paxton、Marzari–Vanderbilt 展宽，`METHOD` 默认值变为 `GAUSSIAN`（宽度由 `SIGMA` 控制），旧版本默认是 `ENERGY_WINDOW`。解决：始终显式写 `METHOD FERMI_DIRAC`（或你需要的方法）。

**12. 用 `cp2kdata` 或旧脚本读 2026.2 的 `.pdos`，各通道错位一列。**
原因：2026.2 在 `Occupation` 后增加了 `Total` 列。解决：按表头列名解析（第 6.3 节代码），或升级/修改解析工具。

---

## 练习

1. 在 H2O 示例中把 `COMPONENTS` 改为 F，比较 `p` 列与 `py+pz+px` 之和；再把分子旋转到 xz 平面，观察 HOMO 的主要分量如何变化。
2. 把石墨烯的 SCF 网格从 21×21×1 改成 12×12×1，保持 `MP_GRID 60 60 1`，比较总 DOS 与 `pz` PDOS 的变化，体会“总 DOS 用密网格、PDOS 用 SCF 网格”的区别。
3. 在 2026.2 上用金红石**原胞** + 6×6×10 k 点重做 PDOS（需要 `&CURVE`），与 72 原子 Γ 点超胞的结果比较带隙和 O 2p/Ti 3d 峰位。
4. Si 能带：分别用 DZVP 和 TZV2P 基组计算，比较 Γ、X、L 处导带能量和带隙。
5. 给 Si 加上 `&KPOINTS SYMMETRY T`，比较 SCF 用到的 k 点数和总能量、DOS 是否改变（对称约化在 2026.2 中是实验性功能，结果有差异时以关闭对称性的结果为准）。

---

## 参考资料

**官方文档**

- Density of States（新的 `&DOS/&PDOS/&CURVE` 接口，2026.2+）：<https://manual.cp2k.org/trunk/methods/electronic_structure/dos.html>
- 输入参考 `&DFT/&PRINT/&DOS`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/DOS.html>
- 输入参考 `&BAND_STRUCTURE`：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/PRINT/BAND_STRUCTURE.html>
- K-Points（含能带路径与 SCF 网格的区别、对称约化）：<https://manual.cp2k.org/trunk/methods/dft/k-points.html>
- How to make a SCF run converge：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- CP2K v2026.2 发布说明（DOS/PDOS 重构 #5326、k 点 PDOS #5287/#5299、新展宽方式 #4958）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- 各版本手册入口（查旧语法）：<https://manual.cp2k.org/trunk/versions.html>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-kp-1`、`regtest-kp-2`、`regtest-ot` 中有 DOS、PDOS 和能带的 regtest 输入，是查关键字真实写法的权威来源
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：`cp2k_bs2csv` 把 `.bs` 转 CSV，`cp2k_pdos` 对 PDOS 做高斯展宽
- [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata)：读 pdos、cube、输出文件；新版本输出格式请先验证
- [mmm-group/dos-broadening](https://github.com/mmm-group/dos-broadening)：对旧格式 `.pdos` 做展宽的小脚本（2026.2 起 `&CURVE` 可直接输出展宽曲线）
- [RobinLi-95/CP2K_tools](https://github.com/RobinLi-95/CP2K_tools)：能带作图等脚本
- [Doria58/DFT_Script](https://github.com/Doria58/DFT_Script)：中文脚本，CP2K 能带转 Origin 数据、生成能带路径
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含能带章节（2025.2 验证，旧语法）
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，CP2K 专题
- [nanotech-empa/cp2k-spm-tools](https://github.com/nanotech-empa/cp2k-spm-tools)：读 `.wfn` 做 STM/STS 与轨道杂化分析（表面体系进阶）

**关键文献**

- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* **49**, 299 (2010)，高对称 k 点路径约定：<https://doi.org/10.1016/j.commatsci.2010.05.010>
- Y. Hinuma *et al.*, *Comput. Mater. Sci.* **128**, 140 (2017)，基于晶体学的能带路径（SeeK-path）：<https://doi.org/10.1016/j.commatsci.2016.10.015>
- J. VandeVondele *et al.*, *Comput. Phys. Commun.* **167**, 103 (2005)，Quickstep：<https://doi.org/10.1016/j.cpc.2004.12.014>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>
- The CP2K Program Package Made Simple, *J. Phys. Chem. B* (2025)：<https://doi.org/10.1021/acs.jpcb.5c05851>

---

[← 上一章](09-hybrid-admm.md) | [返回目录](../README.md) | [下一章 →](11-cube-charge-workfunction.md)
