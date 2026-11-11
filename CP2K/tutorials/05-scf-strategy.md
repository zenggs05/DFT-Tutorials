# 第5章 SCF 策略：OT、对角化、展宽与重启（SCF Convergence in CP2K: OT, Diagonalization, Smearing & Restart）

本章系统讲解 CP2K Quickstep 的自洽场（SCF）求解策略：OT（orbital transformation）直接极小化、对角化加 Broyden/Pulay 密度混合、金属与小带隙体系的 Fermi–Dirac 展宽与 `ADDED_MOS`，以及用 `.wfn`/`.kp` 波函数文件搭建“低精度 → 高精度”的重启链。我们用 H2O（OT）、体相 Si（对角化 + k 点）、fcc Al（金属展宽）三个可运行示例逐个演示，同时说明 2024.1 起 SCF 不收敛即中止的默认行为、2026.2 与 trunk 在 OT + k 点/展宽上的版本差异，最后给出一份 SCF 不收敛的排查流程。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据体系类型（分子/绝缘体、半导体 + k 点、金属）选择 OT 或对角化；
- 正确设置 `&SCF` 的通用关键字（`EPS_SCF`、`MAX_SCF`、`SCF_GUESS`、`ADDED_MOS`）以及 `&DFT/WFN_RESTART_FILE_NAME`；
- 为 OT 选择 `MINIMIZER`、`PRECONDITIONER`，并用 `&OUTER_SCF` 稳定收敛；
- 为对角化设置 `&MIXING`（`BROYDEN_MIXING`/`PULAY_MIXING`、`ALPHA`、`BETA`、`NBUFFER`）和 `&SMEAR`；
- 读懂 SCF 迭代表、`outer SCF` 信息、费米能和熵项，判断收敛是否可信；
- 用重启文件把便宜计算的结果传给昂贵计算，并避开 `.wfn` 与 `.kp` 混用等陷阱。

| 体系 | 推荐求解器 | 示例 |
|---|---|---|
| 分子、团簇、有能隙的大体系（Γ 点） | OT + `&OUTER_SCF` | `h2o_ot.inp` |
| 半导体/绝缘体 + k 点 | 对角化 + Broyden/Pulay 混合 | `si_diag_broyden.inp` |
| 金属、小带隙体系（Γ 或 k 点） | 对角化 + 混合 + `&SMEAR` + 足够的 `ADDED_MOS` | `al_metal_smear.inp` |
| 从已有波函数继续（换基组/泛函/任务） | 任意，`SCF_GUESS RESTART` | `restart_from_wfn.inp` |

示例输入集：[../examples/05-scf-strategy/README.md](../examples/05-scf-strategy/README.md)。

> [!NOTE]
> 本章使用的截断（H2O 500 Ry、Si 350 Ry、Al 300 Ry）来自 [第4章](04-cutoff-convergence.md) 的测试思路，只是示例值。正式计算请先完成你自己体系的截断与 k 点收敛（k 点见 [第6章](06-periodic-surface-kpoints.md)）。

---

## 2. 原理速览

### 2.1 SCF 是一个不动点问题

Kohn–Sham 方程 $\hat H_\mathrm{KS}[\rho]\,\psi_i = \varepsilon_i \psi_i$ 中，哈密顿量依赖密度 $\rho$，而密度又由轨道 $\psi_i$ 决定。SCF 要找的是不动点 $\rho_\mathrm{out}[\rho_\mathrm{in}] = \rho_\mathrm{in}$。CP2K 提供两类求法：

1. **对角化 + 密度混合**：每步构建 KS 矩阵，对角化得到本征态，按占据数生成新的密度矩阵 $P_\mathrm{out}$，再与历史密度混合得到下一步的 $P_\mathrm{in}$。适用面最广，支持展宽、k 点和任意数量的空轨道；
2. **OT（orbital transformation）直接极小化**：把总能量看作占据轨道系数的泛函，在保持正交归一的约束下直接极小化（VandeVondele & Hutter, 2003）。只处理占据轨道，不需要完整对角化；对有能隙的体系，只要预条件合适就能稳定收敛到极小值。

### 2.2 OT 的要点

OT 用一个辅助变量 $X$（满足 $X^\mathrm{T} S C_0 = 0$）参数化占据轨道：

$$
C(X) = C_0 \cos U + X\,U^{-1}\sin U,\qquad U = \left(X^\mathrm{T} S X\right)^{1/2}
$$

这样，无论 $X$ 取什么值，$C(X)$ 都自动满足正交归一，于是带约束的极小化变成无约束问题，可以用 CG、DIIS 等方法求解。每步的主要开销是稀疏矩阵乘法，规模约为 $N_\mathrm{basis}^2 N_\mathrm{occ}$，没有 $N_\mathrm{basis}^3$ 的对角化，所以大体系上明显更快。收敛快慢取决于**预条件子（preconditioner）**：好的预条件子能让迭代次数下降一个数量级以上。

OT 的局限（2026.2 及以前）：它只优化固定数目的占据轨道，不能附加空轨道（`ADDED_MOS` 必须为 0），因此**做不了真正意义上的展宽，也不能用于 k 点**；对金属或几乎无能隙的体系收敛很差。

### 2.3 密度混合

最简单的线性混合是 $P_\mathrm{in}^{(n+1)} = (1-\alpha) P_\mathrm{in}^{(n)} + \alpha P_\mathrm{out}^{(n)}$，`ALPHA` 就是这里的 $\alpha$（默认 0.4）。Broyden 和 Pulay 方法利用前 `NBUFFER` 步的残差历史构造拟牛顿更新，收敛快得多。CP2K 的 `BROYDEN_MIXING`、`PULAY_MIXING`、`KERKER_MIXING` 在倒空间对密度做混合，并带 Kerker 阻尼

$$
\rho_\mathrm{mix}(\mathbf G) = \rho_\mathrm{in}(\mathbf G) + \alpha\,\frac{G^2}{G^2+\beta^2}\left[\rho_\mathrm{out}(\mathbf G)-\rho_\mathrm{in}(\mathbf G)\right]
$$

`BETA`（$\beta$，默认 0.5 bohr⁻¹）压低长波分量的更新，用来抑制大晶胞、金属和 slab 中常见的“电荷晃荡”（charge sloshing）。

### 2.4 展宽（smearing）与电子熵

金属在费米能附近有连续的能级，整数占据会让占据态在 SCF 迭代之间来回跳动，导致不收敛。Fermi–Dirac 展宽给每个态分数占据

$$
f_i = \frac{1}{1+\exp\!\left[(\varepsilon_i-\mu)/k_\mathrm{B}T_\mathrm{el}\right]}
$$

此时被极小化的是自由能 $F = E - T_\mathrm{el}S$。CP2K 输出中的 `Electronic entropic energy` 就是 $-T_\mathrm{el}S$（通常为负），且**已经计入** `Total energy`。零温能量可以用 Gillan 外推近似：$E_0 \approx F - \tfrac12(-T_\mathrm{el}S)$。2026.2 起，Fermi–Dirac 展宽时 CP2K 会直接打印 `Total energy (extrapolated to T->0)` 一行。注意：**力和应力对应的是自由能 $F$**，而不是外推能量。

展宽需要空轨道来容纳费米面以上的分数占据，这就是 `ADDED_MOS` 的作用。

---

## 3. 先说清楚：2024.1 起 SCF 不收敛会中止

从 CP2K 2024.1 开始，Quickstep 的 SCF 达到 `MAX_SCF`（默认 50）仍未满足 `EPS_SCF` 时，程序**默认中止**，输出类似下面的框：

```text
 *******************************************************************************
 *   ___                                                                       *
 *  /   \                                                                      *
 * [ABORT]                                                                     *
 *  \___/     SCF run NOT converged. To continue the calculation regardless,   *
 *    |             please set the keyword IGNORE_CONVERGENCE_FAILURE.         *
 *  O/|                                                                        *
 * /| |                                                                        *
 * / \                                                            qs_scf.F:702 *
 *******************************************************************************
```

（文件名和行号随版本变化。）在 `&SCF` 中设 `IGNORE_CONVERGENCE_FAILURE TRUE` 后，这里只会打印一条 `*** WARNING ... SCF run NOT converged ***`，计算继续往下走。

> [!CAUTION]
> `IGNORE_CONVERGENCE_FAILURE` 只能作为排查问题时的**临时手段**，绝不能放进默认模板。没有收敛的 SCF 给出的能量、力、占据数、能带都不可信；在 MD 中会表现为温度失控，在几何优化中会把结构推向错误的方向。正确做法是按第 10 节的流程让 SCF 真正收敛。

这一行为变化还会影响旧的工作流：

- 依赖“MAX_SCF 1 + 非自洽能量”的旧截断测试脚本，必须显式加上 `IGNORE_CONVERGENCE_FAILURE`（见 [第4章](04-cutoff-convergence.md)）；
- 用 OT 时若开启了 `&OUTER_SCF`，内层循环达到 `MAX_SCF` 只会进入下一个外层循环，不算失败；只有外层也用完 `MAX_SCF` 仍未收敛才会中止。

---

## 4. `&SCF` 通用关键字与重启文件

| 关键字 | 位置 | 默认 | 推荐 | 说明 |
|---|---|---|---|---|
| `EPS_SCF` | `&SCF` | 1.0E-5 | 1.0E-6（常规）；1.0E-7 至 1.0E-8（振动、截断测试） | 收敛阈值。OT 判据是能量梯度的范数，对角化判据是密度矩阵的最大变化，两者不能直接比较 |
| `MAX_SCF` | `&SCF` | 50 | OT 内层 20–30；对角化 50–200 | 用 OT 加 `&OUTER_SCF` 时，这是**内层**最大步数 |
| `SCF_GUESS` | `&SCF` | `ATOMIC` | `ATOMIC` / `RESTART` | `RESTART` 读取波函数文件，找不到文件时回退到 `ATOMIC` 并打印警告 |
| `ADDED_MOS` | `&SCF` | 0 | 金属：占据数的 10–30%，小胞取更多；`-1` 表示全部 | 每个自旋通道额外计算的空轨道数；2026.2 及以前与 OT 不兼容 |
| `IGNORE_CONVERGENCE_FAILURE` | `&SCF` | F | 不要打开 | 见第 3 节 |
| `LEVEL_SHIFT` | `&SCF` | 0 [hartree] | 0.1–0.3（主要用于对角化，作为补救） | 抬高虚轨道能级以抑制振荡；会减慢收敛，收敛后应去掉再确认一次 |
| `EPS_DIIS` / `MAX_DIIS` | `&SCF` | 0.1 / 4 | 一般默认 | 对角化中残差小于 `EPS_DIIS` 后启用 DIIS；OT 默认不切换 |
| `WFN_RESTART_FILE_NAME` | **`&DFT`** | 由 PROJECT 推断 | 显式写出 | 注意它不在 `&SCF` 里 |
| `&PRINT/&RESTART` | `&SCF` | 打开 | `BACKUP_COPIES 0` | 控制 `<project>-RESTART.wfn`（Γ）或 `<project>-RESTART.kp`（k 点）的写出 |

> [!TIP]
> `EPS_SCF` 的选择要与任务匹配。能量对波函数误差是二阶的，所以能量通常比 `EPS_SCF` 本身精确得多；但力、应力和振动频率对 SCF 误差更敏感。经验上，单点能和几何优化用 1.0E-6，振动分析用 1.0E-7 到 1.0E-8，AIMD 用 1.0E-6 配合波函数外推（见 第12章）。

---

## 5. OT 实战：H2O 分子（h2o_ot.inp）

### 5.1 完整输入

```text
# H2O single point + forces with OT (orbital transformation) and an outer SCF loop
&GLOBAL
  PROJECT     h2o_ot
  RUN_TYPE    ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON                       ! isolated molecule, see chapter 6
      PERIODIC       NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   30                 ! inner loop: steps per preconditioner rebuild
      &OT ON
        MINIMIZER      DIIS        ! CG is the default and more robust
        PRECONDITIONER FULL_ALL    ! most effective for small/medium systems
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6             ! should not be looser than the inner EPS_SCF
        MAX_SCF 10                 ! up to 10 x 30 = 300 inner steps in total
      &END OUTER_SCF
      &PRINT
        &RESTART ON                ! writes h2o_ot-RESTART.wfn
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.0 12.0 12.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   h2o.xyz
      COORD_FILE_FORMAT XYZ
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
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

### 5.2 OT 关键字详解

| 关键字（`&SCF/&OT`） | 默认 | 本例 | 说明 |
|---|---|---|---|
| 段参数 `&OT ON` | F | ON | 打开 OT。不能与 `&DIAGONALIZATION ON` 同时出现 |
| `MINIMIZER` | `CG` | `DIIS` | `CG` 最稳健，每步能量都应下降；`DIIS` 对表现良好的体系可快约 50%，但较不可靠；`SD` 不推荐；还有 `BROYDEN` |
| `PRECONDITIONER` | `FULL_KINETIC` | `FULL_ALL` | 见下表 |
| `ENERGY_GAP` | -1（由程序决定） | 默认 | 对 `FULL_ALL` 应**低估** HOMO-LUMO 能隙（如 0.002 Ha）；对 `FULL_SINGLE_INVERSE` 是下限，低于 0.05 可能不稳定 |
| `LINESEARCH` | `2PNT` | 默认 | `CG` 不稳时可改 `3PNT`；`GOLD` 最稳健也最贵 |
| `ALGORITHM` | `STRICT` | 默认 | `IRAC` 用于特定的大体系优化 |
| `STEPSIZE` | -1 | 默认 | 初猜很差时调小（如 0.05）可以稳住前几步 |

预条件子怎么选（据 2026.2 输入参考的说明整理）：

| `PRECONDITIONER` | 构造成本 | 效果 | 适用 |
|---|---|---|---|
| `FULL_ALL` | 高（基于对角化） | 最好 | 中小体系的首选；需要合适的 `ENERGY_GAP` |
| `FULL_SINGLE_INVERSE` | 中（Cholesky 求逆） | 接近 `FULL_SINGLE`，稍不稳健 | 大体系（数百原子以上）推荐 |
| `FULL_KINETIC` | 低 | 稳健、较好 | 默认；超大体系 |
| `FULL_S_INVERSE` | 低 | 不如 `FULL_KINETIC` | 很少用 |
| `NONE` | 无 | 差 | 仅用于测试 |

`&OUTER_SCF` 的作用：内层循环每走 `SCF/MAX_SCF` 步，就退出一次，**重建预条件子**，再开始新的内层循环。官方文档建议内层 `MAX_SCF` 取 16–32 左右、外层取 8–16。总迭代上限是“内层 × 外层”。外层的 `EPS_SCF` 不应比内层宽松，通常两者相同。

### 5.3 运行与检查

```bash
cd CP2K/examples/05-scf-strategy
cp2k.psmp -c h2o_ot.inp                       # syntax check
mpirun -np 4 cp2k.psmp -i h2o_ot.inp -o h2o_ot.out

grep -E "OT (DIIS|CG|LS|SD)" h2o_ot.out | head    # inner OT iterations
grep "outer SCF" h2o_ot.out                   # outer loop summary
grep "SCF run converged" h2o_ot.out
grep "ENERGY| Total FORCE_EVAL" h2o_ot.out
ls -l h2o_ot-RESTART.wfn                      # needed by restart_from_wfn.inp
```

SCF 迭代表大致是这样（示意输出，数值仅供参考）：

```text
 Step     Update method      Time    Convergence         Total energy    Change
 ------------------------------------------------------------------------------
    1 OT DIIS     0.15E+00    0.1     0.03178120       -17.0896475152 -1.71E+01
    2 OT DIIS     0.15E+00    0.1     0.01633541       -17.1531062216 -6.35E-02
    ...
   11 OT DIIS     0.15E+00    0.1     0.00000088       -17.2207634721 -3.15E-11

  *** SCF run converged in    11 steps ***
  ...
  outer SCF iter =    1 RMS gradient =   0.88E-06 energy =        -17.2207634721
  outer SCF loop converged in   1 iterations or   11 steps
```

各列含义：`Update method` 是本步所用的方法（`OT DIIS`、`OT CG`、`OT LS` 表示线搜索步）；`Convergence` 是与 `EPS_SCF` 比较的量（OT 下为梯度）；`Change` 是相对上一步的能量变化。判断成功的依据：

1. 出现 `SCF run converged` 和 `outer SCF loop converged`；
2. `Convergence` 在最后几步单调下降到 `EPS_SCF` 以下，没有反复跳动；
3. 能量变化在最后几步只剩 10⁻⁹ Ha 以下的量级。

> [!TIP]
> 如果 `OT DIIS` 的能量出现上升或大幅振荡，先把 `MINIMIZER` 改回 `CG`；仍然困难时加 `LINESEARCH 3PNT`，并保持 `PRECONDITIONER FULL_ALL`。大体系（上千原子）构造 `FULL_ALL` 太慢时，换 `FULL_SINGLE_INVERSE`。

---

## 6. 对角化 + 密度混合：体相 Si + k 点（si_diag_broyden.inp）

2026.2 及以前的版本中，**k 点计算只能用对角化**。把 `&OT` 和 `&KPOINTS` 放在一起会直接中止，报 `OT not possible with kpoint calculations`。下面用 8 原子 Si 胞加 4×4×4 Monkhorst–Pack 网格演示。

### 6.1 完整输入

```text
# Bulk Si (8-atom cubic cell) with k-points: diagonalization + Broyden density mixing
&GLOBAL
  PROJECT     si8_diag
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     350
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 4
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   80
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD  BROYDEN_MIXING
        ALPHA   0.4                ! fraction of the new density
        BETA    0.5                ! Kerker damping parameter [bohr^-1]
        NBUFFER 8                  ! history length (default 4)
      &END MIXING
      &PRINT
        &RESTART ON                ! k-point runs write si8_diag-RESTART.kp
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 5.431 5.431 5.431
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   si_bulk8.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### 6.2 对角化与混合关键字

| 关键字 | 默认 | 本例 | 说明 |
|---|---|---|---|
| `&DIAGONALIZATION ALGORITHM` | `STANDARD` | `STANDARD` | LAPACK/ScaLAPACK（或 ELPA）完整对角化；k 点只能用它 |
| `&MIXING METHOD` | `DIRECT_P_MIXING` | `BROYDEN_MIXING` | 还可选 `PULAY_MIXING`、`KERKER_MIXING`、`MULTISECANT_MIXING`；2026.2 新增 `MODIFIED_BROYDEN_MIXING`、`NEW_PULAY_MIXING` |
| `ALPHA` | 0.4 | 0.4 | 新密度所占比例。振荡时降到 0.1–0.2 |
| `BETA` | 0.5 bohr⁻¹ | 0.5 | Kerker 阻尼；倒空间混合方法不允许设为 0 |
| `NBUFFER` | 4 | 8 | 历史步数（别名 `NBROYDEN`、`NPULAY`）；难收敛体系可取 8–20 |
| `NSKIP` | 0 | 默认 | 前几步不混合 |
| `ALPHA_MAG` / `BETA_MAG` | -1（同 ALPHA/BETA） | — | **2026.2 新增**：自旋极化时磁化密度通道单独的混合参数，过渡金属磁性体系可把 `ALPHA_MAG` 设得比 `ALPHA` 大 |

> [!NOTE]
> `&MIXING` 只对对角化（以及线性标度 SCF）起作用，OT 会忽略它。另一方面，用对角化时**必须**有混合：若设 `METHOD NONE`，CP2K 会中止并提示 `Diagonalization procedures without mixing are not recommendable`。

### 6.3 运行与检查

```bash
mpirun -np 8 cp2k.psmp -i si_diag_broyden.inp -o si_diag_broyden.out
grep "BRILLOUIN| List of Kpoints" si_diag_broyden.out   # number of irreducible k-points
grep -E "Broy\./Diag\.|Pulay/Diag\.|DIIS/Diag\.|Mix/Diag\." si_diag_broyden.out | head
grep "SCF run converged" si_diag_broyden.out
ls -l si8_diag-RESTART.kp
```

对角化的迭代行以 `NoMix/Diag.`（第一步没有可混合的历史）开头，随后是 `Broy./Diag.`（Broyden）或 `Pulay/Diag.` 等。`Convergence` 列是密度矩阵变化，通常比 OT 的梯度下降得更“台阶式”。如果它在 10⁻³–10⁻⁴ 之间来回跳而不再下降，说明出现了电荷晃荡或占据数翻转，处理方法见第 9 节。

---

## 7. 金属与小带隙体系：fcc Al + Fermi–Dirac 展宽（al_metal_smear.inp）

### 7.1 完整输入

```text
# fcc Al (conventional 4-atom cell): metal -> diagonalization + mixing + smearing + ADDED_MOS
&GLOBAL
  PROJECT     al_fcc_smear
  RUN_TYPE    ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     300
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 8 8 8   ! converge this for production (chapter 6)
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   100
      ADDED_MOS 10                  ! 6 occupied MOs per k-point -> 16 MOs in total
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD  BROYDEN_MIXING
        ALPHA   0.2                 ! metals: smaller ALPHA is safer
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC          ! always give METHOD explicitly (default changed in 2026.2)
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 4.05 4.05 4.05            ! experimental lattice constant ~4.05 Angstrom
      PERIODIC XYZ
    &END CELL
    &COORD
      SCALED
      Al  0.0  0.0  0.0
      Al  0.0  0.5  0.5
      Al  0.5  0.0  0.5
      Al  0.5  0.5  0.0
    &END COORD
    &KIND Al
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q3
      POTENTIAL GTH-PBE-q3
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

### 7.2 展宽相关关键字

| 关键字（`&SCF/&SMEAR`） | 默认 | 本例 | 说明 |
|---|---|---|---|
| 段参数 `&SMEAR ON` | F | ON | 打开分数占据 |
| `METHOD` | ≤2026.1 文档默认 `ENERGY_WINDOW`；2026.2 默认 `GAUSSIAN` | `FERMI_DIRAC` | **务必显式写出**。各版本都有 `FERMI_DIRAC`、`ENERGY_WINDOW`、`LIST`；2026.2 新增 `GAUSSIAN`、`METHFESSEL_PAXTON`、`MARZARI_VANDERBILT` |
| `ELECTRONIC_TEMPERATURE` | 300 K | 500 K | Fermi–Dirac 的电子温度。常用 300–1000 K；越高越好收敛，但熵误差越大 |
| `SIGMA` | 0.002 Ha | — | **2026.2 新增**，用于 Gaussian/MP/MV 展宽的宽度，如 `SIGMA [eV] 0.1` |
| `FIXED_MAGNETIC_MOMENT` | -100（不固定） | — | 自旋极化时固定 $n_\uparrow - n_\downarrow$ |
| `ADDED_MOS`（`&SCF`） | 0 | 10 | 见下文 |

2026.2 中可以用冷展宽（Marzari–Vanderbilt）写成：

```text
&SMEAR ON
  METHOD MARZARI_VANDERBILT        ! CP2K >= 2026.2 only
  SIGMA [eV] 0.1
&END SMEAR
```

> [!WARNING]
> `METHFESSEL_PAXTON` 可能给出负的或大于 1 的占据数，**不要用于半导体和绝缘体**。Gaussian 展宽同样需要把 SIGMA 外推到 0 才有意义。如果输入要在 2026.2 和旧版本之间共用，只用 `FERMI_DIRAC`。

**ADDED_MOS 取多少？** 默认情况下 CP2K 只计算占据轨道。展宽时费米面以上几个 $k_\mathrm{B}T$ 范围内的态都会有分数占据，必须计算出来。原则如下：

- 最高的那条轨道在每个 k 点上的占据数应当接近 0（低于 `EPS_FERMI_DIRAC`）。不满足时，CP2K 会打印警告 `Fermi-Dirac smearing includes the last MO => Add more MOs for proper smearing.`；
- 大体系（上百原子）通常取占据轨道数的 10–30%；本例胞很小（每个 k 点只有 6 条占据轨道），能带色散又大，所以取 10 条；
- `ADDED_MOS -1` 表示计算全部轨道（即基函数个数），最保险但最贵，可用来验证某个 ADDED_MOS 是否足够；
- 不能超过基函数总数（本例 4 × 13 = 52）。

### 7.3 运行与检查

```bash
mpirun -np 8 cp2k.psmp -i al_metal_smear.inp -o al_metal_smear.out
grep "Fermi energy"                        al_metal_smear.out | tail -1
grep "Electronic entropic energy"          al_metal_smear.out | tail -1
grep "Total energy (extrapolated to T->0)" al_metal_smear.out | tail -1   # 2026.2+
grep -i "smearing includes the last MO"    al_metal_smear.out            # must be empty
grep "SCF run converged"                   al_metal_smear.out
```

能量部分大致是这样（示意输出，数值仅供参考）：

```text
  Electronic entropic energy:                                -0.00041826317455
  Fermi energy:                                               0.27631402155712

  Total energy:                                              -8.38512716294132
  Total energy (extrapolated to T->0):                       -8.38491803135404
```

判断要点：

1. `Electronic entropic energy` 的绝对值相对能量差应当很小（例如小于 1 meV/原子）；若很大，说明电子温度过高；
2. 计算能量差时，各结构要使用**同一个**电子温度，并统一采用自由能或外推能量中的一种；
3. 力对应自由能 $F$。几何优化和 MD 中，用这个温度得到的力是自洽的，但需要确认展宽对结构的影响可以忽略（例如对比 300 K 与 500 K 的晶格常数）。

---

## 8. 版本差异：OT、k 点与展宽

| 功能 | ≤ 2025.2 | 2026.1 | 2026.2（当前稳定版） | trunk（2026.2 之后） |
|---|---|---|---|---|
| OT + Γ 点（整数占据） | 支持 | 支持 | 支持 | 支持 |
| OT + k 点 | 不支持，中止：`OT not possible with kpoint calculations` | 同左 | 同左（k 点兼容表标为 Unsupported） | 支持（复数 k 点轨道） |
| OT + `ADDED_MOS` / 展宽 | 不支持，中止：`OT with ADDED_MOS/=0 not implemented` | 同左 | 同左 | 支持：`ROTATION` + `ENERGIES`，可选 `OCCUPATION_PRECONDITIONER`，`ADDED_MOS AUTO` |
| `&SMEAR METHOD` 可选值 | FERMI_DIRAC / ENERGY_WINDOW / LIST | 同左 | 新增 GAUSSIAN / METHFESSEL_PAXTON / MARZARI_VANDERBILT 和 `SIGMA` | 同 2026.2 |
| `&MIXING` 新选项 | — | — | `MODIFIED_BROYDEN_MIXING`、`NEW_PULAY_MIXING`、`ALPHA_MAG`、`BETA_MAG` | 同 2026.2 |
| k 点波函数外推（GEO_OPT/MD） | 不支持 | 不支持 | 支持 | 支持 |
| 输出 `Total energy (extrapolated to T->0)` | 无 | 无 | 有 | 有 |

trunk 版本中 OT 处理金属的写法如下，仅作了解；**2026.2 及以前会直接报错**：

```text
# trunk only (after 2026.2) -- NOT accepted by CP2K <= 2026.2
&SCF
  ADDED_MOS AUTO                   ! grows the virtual buffer automatically
  &SMEAR ON
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
  &OT ON
    MINIMIZER      CG
    PRECONDITIONER FULL_ALL
    ROTATION                       ! variational rotations within the occupied space
    ENERGIES                       ! auxiliary orbital energies -> occupations
    OCCUPATION_PRECONDITIONER      ! optional, system dependent
  &END OT
&END SCF
```

官方 OT 文档还建议：用这条路径计算新的金属体系时，至少拿一个代表性结果与“对角化 + 密度混合”对比，并检查最高能带没有被意外占据。在下一个正式版本发布、你的集群升级之前，金属和 k 点计算请继续使用第 6、7 节的对角化方案。

---

## 9. 重启链：.wfn 与 .kp

### 9.1 从低精度到高精度：restart_from_wfn.inp

先运行 `h2o_ot.inp`（DZVP，500 Ry），得到 `h2o_ot-RESTART.wfn`；再用它作为初猜，计算 TZV2P 基组、600 Ry、更紧 EPS_SCF 的版本：

```text
# Restart a TZV2P calculation of H2O from the DZVP wavefunction of h2o_ot.inp
# Run h2o_ot.inp first: it writes h2o_ot-RESTART.wfn
&GLOBAL
  PROJECT     h2o_tzv2p            ! new project name: the source .wfn is not overwritten
  RUN_TYPE    ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME   BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME   POTENTIAL_UZH
    WFN_RESTART_FILE_NAME h2o_ot-RESTART.wfn   ! keyword of &DFT, not of &SCF
    &MGRID
      CUTOFF     600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC       NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS RESTART            ! falls back to ATOMIC (with a warning) if the file is missing
      EPS_SCF   1.0E-7
      MAX_SCF   30
      &OT ON
        MINIMIZER      DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.0 12.0 12.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   h2o.xyz
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND H
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q1   ! same pseudopotential as in h2o_ot.inp
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

运行与检查：

```bash
mpirun -np 4 cp2k.psmp -i h2o_ot.inp           -o h2o_ot.out
mpirun -np 4 cp2k.psmp -i restart_from_wfn.inp -o restart_from_wfn.out
grep -i "does not exist" restart_from_wfn.out   # must be empty
grep "SCF run converged" h2o_ot.out restart_from_wfn.out
```

如果重启文件不存在，CP2K 会打印类似下面的警告，并改用 ATOMIC 初猜（计算不会中止，所以一定要检查）：

```text
 *** WARNING in qs_initial_guess.F:312 :: User requested to restart the    ***
 *** wavefunction from the file named: project_name-RESTART.wfn. This file ***
 *** does not exist. Please check the existence of the file or change      ***
 *** properly the value of the keyword WFN_RESTART_FILE_NAME. Calculation  ***
 *** continues using ATOMIC GUESS.                                         ***
```

（摘自官方 Restarting 文档中的示例；行号随版本变化。）

重启成功时，第二个计算的 SCF 步数应明显少于从 ATOMIC 起步（可以把 `SCF_GUESS` 改成 `ATOMIC` 跑一次作对比）。

### 9.2 常用的重启链

| 链条 | 做法 | 注意 |
|---|---|---|
| 小基组 → 大基组 | DZVP 的 `.wfn` 作为 TZV2P 的初猜 | **赝势必须相同**；Γ 点格式允许基函数数不同 |
| 低截断/松 EPS → 高精度 | 直接 `SCF_GUESS RESTART` | 原子数和元素顺序不能变 |
| PBE → PBE0/HSE06 | 先 PBE 收敛，再以其 `.wfn` 起步做杂化计算，并可开 `SCREEN_ON_INITIAL_P` | 见 [第9章](09-hybrid-admm.md) |
| ENERGY → GEO_OPT → VIBRATIONAL_ANALYSIS | 每一步都读上一步的 `.wfn` | 几何变化后 `.wfn` 只是初猜；优化与 MD 中的初猜主要靠 `&QS EXTRAPOLATION`（默认 `ASPC`） |
| 闭壳层 → 开壳层（UKS） | Γ 点格式只在“重启文件的自旋数少于输入、且输入 `MULTIPLICITY` 为 1”时接受 | 自旋态改变时通常应重新初猜，见 [第8章](08-dispersion-spin-dftu.md) |
| 中断的长任务续算 | 运动部分用 `*-1.restart` 加 `&EXT_RESTART`，电子部分用 `.wfn` | 见 [第7章](07-geo-cell-opt.md)、第12章 |

### 9.3 .wfn 与 .kp 不能互换

| | Γ 点（不写 `&KPOINTS`） | k 点（写了 `&KPOINTS`，哪怕 `SCHEME GAMMA`） |
|---|---|---|
| 重启文件 | `<project>-RESTART.wfn` | `<project>-RESTART.kp` |
| 内容 | 原子基组信息、每个自旋的 MO 系数、本征值、占据数 | 每个自旋、每个晶胞像的密度矩阵 |
| 原子数不同 | 拒绝，回退初猜 | 拒绝 |
| 基函数数不同 | 可以 | 拒绝 |
| 自旋数不同 | 有限条件下可以 | 拒绝 |

`.wfn` 和 `.kp` 是两套独立的读写机制，**不能互相重启**。此外，这两种文件都是二进制（unformatted）格式，不保证跨机器、跨版本兼容。不要把它们当作后处理的数据来源；轨道分析请用 molden/cube 输出（见 [第11章](11-cube-charge-workfunction.md)）。

> [!TIP]
> 想把重启后的结构写回输入，或者从 `-pos-1.xyz` 取最后一帧，可以用本仓库的 ../scripts/restart_update.py（`python ../scripts/restart_update.py -h` 查看用法）。

---

## 10. SCF 不收敛的排查流程

建议按“从物理到数值”的顺序逐项排查，每次只改一处：

1. **结构**：原子是否重叠、缺 H 或多 H、晶胞是否写错（xyz 文件没有晶胞信息，CIF 可能有无序位点）。用可视化软件连同周期像一起检查。
2. **电子态**：`CHARGE`、`MULTIPLICITY` 是否正确；奇数电子必须开 `UKS`；磁性体系用 `&KIND/MAGNETIZATION` 给初猜（见 [第8章](08-dispersion-spin-dftu.md)）。
3. **基组与赝势**：q 值是否一致；体系里有没有过于弥散的基函数（凝聚相优先用 SR 或 UZH 基组）。
4. **精度参数**：`CUTOFF`/`REL_CUTOFF`/`EPS_DEFAULT` 太低时 SCF 也会难收敛。适当提高反而可能减少总步数。
5. **算法**：
   - 用 OT 的有能隙体系：`MINIMIZER CG`，`PRECONDITIONER FULL_ALL`（大体系用 `FULL_SINGLE_INVERSE`），加 `&OUTER_SCF`（内层 20–30，外层 8–16）；仍振荡时加 `LINESEARCH 3PNT` 或减小 `STEPSIZE`；
   - 金属或小带隙体系：改用对角化，`BROYDEN_MIXING` 或 `PULAY_MIXING`，`ALPHA` 降到 0.1–0.2，`NBUFFER` 增到 8–15，加 `&SMEAR` 和足够的 `ADDED_MOS`；大晶胞或 slab 中出现电荷晃荡时可增大 `BETA`；
   - 对角化中仍振荡：临时加 `LEVEL_SHIFT 0.1`（Ha），或者先用较高的电子温度收敛，再以 RESTART 降温。
6. **初猜**：从便宜计算（小基组、低截断、纯 GGA、较高电子温度）的 `RESTART` 文件起步。
7. **最后手段**：只有在确认问题出在某个瞬时构型（比如 MD 中的个别步）且后果可控时，才考虑 `IGNORE_CONVERGENCE_FAILURE`，并且要事后逐步检查结果。

本仓库的 ../scripts/cp2k_outparse.py 能汇总输出中的 SCF 收敛状态、步数和 WARNING（`python ../scripts/cp2k_outparse.py -h`），便于批量检查。更完整的报错目录见 第18章。

---

## 参数选择与收敛建议

| 场景 | 求解器 | 关键设置 | EPS_SCF |
|---|---|---|---|
| 分子 / 团簇（Γ） | OT | `MINIMIZER DIIS`（不稳时 `CG`），`PRECONDITIONER FULL_ALL`，`&OUTER_SCF` 10 × 内层 30 | 1.0E-6（振动 1.0E-7 至 1.0E-8） |
| 大体系绝缘体（数百原子以上，Γ） | OT | `PRECONDITIONER FULL_SINGLE_INVERSE` 或 `FULL_KINETIC`，`MINIMIZER CG`，`&OUTER_SCF` | 1.0E-6 |
| 半导体 + k 点 | 对角化 | `BROYDEN_MIXING`，`ALPHA 0.3–0.4`，`NBUFFER 8` | 1.0E-6 至 1.0E-7 |
| 金属 + k 点 | 对角化 + 展宽 | `ALPHA 0.1–0.2`，`NBUFFER 8–15`，`FERMI_DIRAC` 300–1000 K，`ADDED_MOS` 足够大 | 1.0E-6 至 1.0E-7 |
| 金属大超胞（Γ） | 对角化 + 展宽 | 同上，`ADDED_MOS` 取占据数的 10–30% | 1.0E-6 |
| 杂化泛函 | 先 PBE 再 RESTART | OT + `FULL_ALL`（Γ）；见第9章 | 1.0E-6 |
| AIMD | OT（绝缘体） | `&QS EXTRAPOLATION ASPC`，`EXTRAPOLATION_ORDER 3` | ≤ 1.0E-6 |

---

## 常见问题与报错

**1. `SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`**
原因：2024.1 起的默认中止行为，SCF 在 `MAX_SCF`（或外层 `MAX_SCF`）内没有收敛。解决：按第 10 节排查；不要把 `IGNORE_CONVERGENCE_FAILURE` 当作解决方案。

**2. `Use the LSD option for an odd number of electrons`**
原因：电子数为奇数却没开自旋极化（不加展宽时才检查）。解决：在 `&DFT` 中加 `UKS`（等价于 `LSD`）并设正确的 `MULTIPLICITY`。如果体系本应是闭壳层，多半是结构中少了或多了 H 原子，或者 `CHARGE` 写错了。

**3. `Don't activate OT and Diagonaliztion together`**（原文拼写如此）
原因：`&OT ON` 与 `&DIAGONALIZATION ON` 同时出现。解决：二选一。两个段都不写时，默认使用对角化。

**4. `OT not possible with kpoint calculations`**
原因：2026.2 及以前的版本中，OT 不支持 k 点。解决：删掉 `&OT`，改用 `&DIAGONALIZATION` 加 `&MIXING`（第 6 节）；或者大超胞只用 Γ 点（删掉 `&KPOINTS`）。

**5. `OT with ADDED_MOS/=0 not implemented`**
原因：2026.2 及以前的版本中，OT 不能与 `ADDED_MOS`（也就不能与展宽）同时使用。解决：金属改用对角化；如果只是想得到 LUMO 等空轨道信息，请用相应的 `&PRINT` 选项在 SCF 之后计算（见 [第10章](10-dos-pdos-band.md)）。

**6. `Fermi-Dirac smearing includes the last MO => Add more MOs for proper smearing.`**
原因：`ADDED_MOS` 太少，费米分布的尾巴被截断，能量和力都会出错，SCF 也可能振荡。解决：增大 `ADDED_MOS`，或用 `ADDED_MOS -1` 验证；同时检查电子温度是否过高。

**7. `Diagonalization procedures without mixing are not recommendable`**
原因：对角化时把 `&MIXING METHOD` 设成了 `NONE`。解决：使用 `DIRECT_P_MIXING`（默认）、`BROYDEN_MIXING` 或 `PULAY_MIXING`。

**8. `Mixing employing the Kerker damping factor needs BETA /= 0.0`**
原因：倒空间混合方法（Broyden、Pulay、Kerker 等）中 `BETA` 被设为 0。解决：恢复默认 0.5 或设一个正值。

**9. 重启没有生效，SCF 步数和从头算一样多。**
原因：`WFN_RESTART_FILE_NAME` 写在了 `&SCF` 中（应在 `&DFT`），或者文件名、路径不对；原子数或元素顺序变了；在 Γ 与 k 点格式之间互相重启。解决：查 `grep -i "does not exist"`，核对文件和体系是否一致。

**10. 对角化中 `Convergence` 卡在 10⁻³–10⁻⁴ 来回跳。**
原因：电荷晃荡、费米面附近占据数翻转，或者混合过激进。解决：`ALPHA` 降到 0.1，`NBUFFER` 增大，加展宽或提高电子温度，大晶胞或 slab 增大 `BETA`；必要时换 `PULAY_MIXING` 对比。

---

## 练习

1. 把 `h2o_ot.inp` 中的 `PRECONDITIONER` 依次换成 `FULL_ALL`、`FULL_SINGLE_INVERSE`、`FULL_KINETIC`，`MINIMIZER` 分别用 `DIIS` 和 `CG`，记录 SCF 步数和耗时，填一张 3 × 2 的表。
2. 在 `si_diag_broyden.inp` 中把 `METHOD` 换成 `PULAY_MIXING`、`DIRECT_P_MIXING`，比较收敛步数；再把 `ALPHA` 改成 0.1 和 0.8，观察迭代行为的变化。
3. 对 `al_metal_smear.inp` 分别用 300、500、1000、2000 K 计算，画出 `Total energy`（自由能）和外推能量随温度的变化，验证外推能量对温度更不敏感。若你使用 2026.2，再用 `MARZARI_VANDERBILT` 加 `SIGMA` 做同样的比较。
4. 把 `al_metal_smear.inp` 的 `ADDED_MOS` 改成 1，看看是否出现 `smearing includes the last MO` 警告，并比较能量的变化。
5. 先运行 `h2o_ot.inp`，再分别以 `SCF_GUESS RESTART` 和 `SCF_GUESS ATOMIC` 运行 `restart_from_wfn.inp`，比较两者的 SCF 步数。

---

## 参考资料

**官方文档**

- How to make a SCF run converge：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- Orbital Transformation（trunk，含 OT + 展宽/k 点的新功能）：<https://manual.cp2k.org/trunk/methods/dft/orbital_transformation.html>
- Restarting CP2K Calculations（`.wfn`/`.kp`、`EXT_RESTART`）：<https://manual.cp2k.org/trunk/methods/restarting.html>
- 输入参考 `&SCF`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/SCF.html>
- 输入参考 `&OT`、`&OUTER_SCF`、`&MIXING`、`&SMEAR`、`&DIAGONALIZATION`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/SCF/OT.html>（同目录下替换段名即可）
- K-Points（功能兼容表）：<https://manual.cp2k.org/trunk/methods/dft/k-points.html>
- Troubleshooting：<https://manual.cp2k.org/trunk/getting-started/troubleshooting.html>
- CP2K 2026.2 release notes：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/` 下的 regtest 覆盖 OT、对角化、展宽、k 点重启（如 `tests/QS/regtest-kp-1/`）的真实写法
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：中文 Wiki，有 SCF 设置与 DFT+U、slab 等实践经验
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程（2025.2 验证）
- [ChiahsinChu/cp2k-user-tutorial](https://github.com/ChiahsinChu/cp2k-user-tutorial)：零基础教程，含密度混合的讲解
- [szt2006/cp2k-aimd-skill](https://github.com/szt2006/cp2k-aimd-skill)：中文“症状到处方”对照表，可作为排查清单参考（license 不明，仅链接）
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：OT/DIAG 一键切换、restart 回写的思路（无 license，仅参考）
- [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools)：`cp2kparse` 把输出解析为 JSON/YAML，便于批量检查收敛

**关键文献**

- J. VandeVondele, J. Hutter, *J. Chem. Phys.* **118**, 4365 (2003)，OT 方法：<https://doi.org/10.1063/1.1543154>
- V. Weber, J. VandeVondele, J. Hutter, A. M. N. Niklasson, *J. Chem. Phys.* **128**, 084113 (2008)，正交约束下的直接极小化：<https://doi.org/10.1063/1.2841077>
- P. Pulay, *Chem. Phys. Lett.* **73**, 393 (1980)，DIIS：<https://doi.org/10.1016/0009-2614(80)80396-4>
- D. D. Johnson, *Phys. Rev. B* **38**, 12807 (1988)，修正 Broyden 混合：<https://doi.org/10.1103/PhysRevB.38.12807>
- G. P. Kerker, *Phys. Rev. B* **23**, 3082 (1981)，Kerker 阻尼：<https://doi.org/10.1103/PhysRevB.23.3082>
- N. D. Mermin, *Phys. Rev.* **137**, A1441 (1965)，有限温度 DFT：<https://doi.org/10.1103/PhysRev.137.A1441>
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* **82**, 3296 (1999)，冷展宽：<https://doi.org/10.1103/PhysRevLett.82.3296>
- M. Methfessel, A. T. Paxton, *Phys. Rev. B* **40**, 3616 (1989)：<https://doi.org/10.1103/PhysRevB.40.3616>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](04-cutoff-convergence.md) | [返回目录](../README.md) | [下一章 →](06-periodic-surface-kpoints.md)
