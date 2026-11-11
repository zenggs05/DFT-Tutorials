# 第6章 杂化泛函与 meta-GGA（HSE06 / r2SCAN / mBJ）（Hybrid Functionals & meta-GGA in VASP: HSE06, r2SCAN, mBJ）

PBE 等半局域泛函会系统性低估带隙。本章讲 VASP 中三类"超越 GGA"的常用方法：HSE06 杂化泛函（LHFCALC、HFSCREEN、AEXX、ALGO=Damped、PRECFOCK、NKRED）、r2SCAN meta-GGA（METAGGA=R2SCAN、LASPH、LMIXTAU，可配 rVV10 或 DFT-D4）以及只用于带隙的 mBJ 势（METAGGA=MBJ、CMBJ）。重点讲清楚为什么杂化泛函和 meta-GGA 的能带**不能用 ICHARG=11**，以及零权重 k 点法和 KPOINTS_OPT 法的完整操作，并比较 Si、单层 MoS2、ZnO 的 PBE/HSE06/mBJ 带隙与计算代价。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 从 PBE 波函数出发，稳定地收敛一个 HSE06 自洽计算，并合理使用 `PRECFOCK`、`NKRED`、`KPAR` 控制计算量；
- 用**零权重 k 点**或 **KPOINTS_OPT** 计算 HSE06 能带，知道何时需要 `HFRCUT = -1`；
- 用 r2SCAN 做结构优化（可加 rVV10 或 D4 色散修正），知道 POTCAR 和版本上的坑；
- 用 mBJ 便宜地估算带隙，并清楚它**不能**用于能量和受力；
- 根据研究目标在 PBE、r2SCAN、HSE06、mBJ 之间做出选择。

| 目标 | 推荐 | 说明 |
|---|---|---|
| 半导体/绝缘体带隙、能带排列、缺陷能级 | HSE06（或调过 AEXX 的 HSE） | 最常用的"可信带隙"方法，代价高 |
| 大量材料快速估带隙 | mBJ | 代价接近 meta-GGA，带隙常接近实验，但只是势，没有能量 |
| 结构、形成能、多种成键类型兼顾 | r2SCAN（+ rVV10/D4） | 比 PBE 更均衡，代价约为 PBE 的数倍 |
| 金属、磁性体系的常规计算 | PBE / PBE+U | 杂化泛函对金属代价高且未必更好 |
| 准粒子能带、激子 | GW / BSE | 超出本章，见 [第14章](14-optics-dielectric.md) 的指引 |

---

## 2. 原理速览

### 2.1 Jacob 阶梯

LDA → GGA（依赖 $\nabla\rho$）→ meta-GGA（再依赖动能密度 $\tau=\tfrac12\sum_i|\nabla\psi_i|^2$）→ 杂化泛函（混入非局域 Hartree–Fock 交换）→ RPA 等。越往上越准、越贵。

### 2.2 HSE06

HSE 把库仑相互作用按 $\mathrm{erfc}(\mu r)/r + \mathrm{erf}(\mu r)/r$ 拆成短程与长程两部分，只在短程部分混入精确交换：

$$
E_{xc}^{\mathrm{HSE}} = a\,E_{x}^{\mathrm{HF,SR}}(\mu) + (1-a)\,E_{x}^{\mathrm{PBE,SR}}(\mu) + E_{x}^{\mathrm{PBE,LR}}(\mu) + E_{c}^{\mathrm{PBE}},
$$

HSE06 取 $a = 0.25$、$\mu = 0.2\ \mathrm{\AA^{-1}}$（VASP 中 `AEXX = 0.25`、`HFSCREEN = 0.2`；HSE03 对应 `HFSCREEN = 0.3`；不设 `HFSCREEN` 则是 PBE0）。短程化使它对 k 点收敛更快，也让 Fock 交换可以在更粗的 q 网格上计算（`NKRED`）。

计算代价：Fock 交换需要对 k 和 q 双重求和，大约随 $N_k \times N_q \times N_{\mathrm{occ}}^2$ 增长，通常比 PBE 慢 1–2 个数量级，内存也大得多。

### 2.3 为什么杂化与 meta-GGA 的能带不能用 ICHARG=11

GGA 的哈密顿量只由电荷密度决定，固定 CHGCAR 就固定了哈密顿量。杂化泛函的 Fock 算符依赖**规则 k 网格上的所有占据轨道**；meta-GGA 依赖由规则网格轨道算出的动能密度 $\tau$。只读 CHGCAR 无法重建它们，因此：

- 杂化泛函**禁止**使用 `ICHARG = 11`（VASP 官方明确警告）；
- meta-GGA 的能带也按与杂化泛函相同的方法计算（规则网格必须参与）。

可行的两种做法：**零权重 k 点**（规则网格 + 权重为 0 的路径点一起自洽）和 **KPOINTS_OPT**（自洽后在路径上做一次性计算）。

### 2.4 r2SCAN 与 mBJ

- **SCAN**（2015）满足全部 17 个已知的精确约束，但数值上对积分网格敏感；**r2SCAN**（2020）在保持大部分约束的同时做了正则化，数值稳定得多，是目前 VASP 中推荐的 SCAN 类泛函（`METAGGA = R2SCAN`，VASP ≥ 6.2.0）。
- **mBJ**（Tran–Blaha 修正的 Becke–Johnson 交换势）只是一个**交换势**：

$$
v_x^{\mathrm{mBJ}}(\mathbf r) = c\,v_x^{\mathrm{BR}}(\mathbf r) + (3c-2)\frac{1}{\pi}\sqrt{\frac{5}{12}}\sqrt{\frac{2\,\tau(\mathbf r)}{\rho(\mathbf r)}},\qquad
c = \alpha + \beta\left(\frac{1}{V_{\mathrm{cell}}}\int_{\mathrm{cell}}\frac{|\nabla\rho|}{\rho}\,d^3r\right)^{1/2}
$$

VASP 默认 $\alpha=-0.012$、$\beta=1.023\ \mathrm{bohr}^{1/2}$（`CMBJA`、`CMBJB`），$c$ 每个电子步自洽更新，也可以用 `CMBJ` 固定。没有对应的能量泛函（VASP 用 LDA 能量表达式占位），所以**不能做结构优化、不能比较总能**，只用来算带隙等电子结构。

---

## 3. 通用准备：从 PBE 出发

1. **结构**：先在 PBE（或 r2SCAN）下优化好结构（[第4章](04-relax-eos.md)）。HSE 计算最常用"PBE 结构 + HSE 单点"，记作 HSE06@PBE；发表时写清楚用的是谁的结构。
2. **POTCAR**：meta-GGA（r2SCAN、mBJ）需要 POTCAR 含芯电子动能密度信息：

   ```bash
   grep kinetic POTCAR
   # each element should show lines such as "kinetic energy-density" and "mkinetic energy-density pseudized"
   ```

   少数 POTCAR（例如部分 `_GW` 变体）缺少这些信息；H、He、Li_sv、Be_sv 这类无芯电子的势不显示 pseudized 一行，但仍可使用。mBJ 这类与 PBE 差别很大的势，VASP 建议测试更硬、半芯态更多的 PAW（如 `_sv`、`_GW`）。
3. **先跑 PBE 得到 WAVECAR**：杂化和 meta-GGA 从 PBE 波函数起步，收敛快得多也稳得多。PBE 与后续计算的 `ENCUT`、`NBANDS`、k 网格要一致（k 点集合不同时 WAVECAR 无法复用）。
4. **k 网格**：杂化泛函的规则网格往往可以比 PBE 稍粗（Si 用 6×6×6 已足够做示范），但必须做收敛测试；若想用 `NKRED = 2`，网格每个方向都要能被 2 整除。

---

## 4. HSE06 自洽计算

示例目录：[../examples/06-hse06-band-si/README.md](../examples/06-hse06-band-si/README.md)（`0-pbe-scf`、`1-hse-scf`）。

### 4.1 第 0 步：PBE 自洽（产生 WAVECAR）

```text
SYSTEM  = Si PBE SCF (WAVECAR for the HSE06 runs)

PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
NBANDS  = 16         # keep NBANDS identical in all following runs

NSW     = 0
IBRION  = -1
LWAVE   = .TRUE.     # WAVECAR is the starting point of the hybrid runs
LCHARG  = .FALSE.
```

KPOINTS（之后的 HSE 自洽和 KPOINTS_OPT 计算都用同一个网格）：

```text
Gamma-centered 6x6x6
0
Gamma
6 6 6
0 0 0
```

### 4.2 第 1 步：HSE06 自洽

```text
SYSTEM   = Si HSE06 SCF (starts from the PBE WAVECAR)

ISTART   = 1         # read WAVECAR copied from 0-pbe-scf
PREC     = Accurate
ENCUT    = 400
EDIFF    = 1E-6
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
NBANDS   = 16

# HSE06
LHFCALC  = .TRUE.
GGA      = PE
HFSCREEN = 0.2       # 0.2 1/A -> HSE06 (0.3 -> HSE03, unset -> PBE0)
AEXX     = 0.25      # default for LHFCALC=.TRUE., written for clarity
ALGO     = Damped    # robust; ALGO = Normal is a fast alternative; never ALGO = Fast
TIME     = 0.4       # damping time step for ALGO = Damped; reduce (0.2-0.3) if SCF oscillates
PRECFOCK = Normal    # FFT grid for the exchange; Fast for pre-convergence, Accurate for phonons
HFRCUT   = -1        # truncated Coulomb kernel: recommended for gapped systems and band structures
# NKRED  = 2         # optional: Fock operator on a 3x3x3 q-mesh; test the effect on the gap first

NSW      = 0
IBRION   = -1
LWAVE    = .TRUE.
LCHARG   = .FALSE.
```

逐项说明：

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `LHFCALC` | .TRUE. | 打开 Hartree–Fock/杂化泛函 |
| `GGA` | PE | 半局域部分用 PBE；`HFSCREEN` 只能与 `GGA = PE/PS/CA` 组合 |
| `HFSCREEN` | 0.2 | HSE06 的屏蔽参数 μ（Å⁻¹） |
| `AEXX` | 0.25 | 精确交换比例；可调（见 §7.2） |
| `ALGO` | Damped（或 All / Normal） | 官方推荐直接最小化类算法（Damped/All）；`Normal` 配合默认的 `LFOCKACE = .TRUE.` 往往更快；`Fast` **不被正确支持且不会报警告** |
| `TIME` | 0.4（默认） | Damped/All 的步长；不收敛时减小 |
| `PRECFOCK` | Normal | Fock 交换的 FFT 网格；`Fast` 快 2–3 倍但力有噪声，可用于预收敛；比较能量时保持一致 |
| `NKRED` | 1（默认）或 2 | q 网格降采样，计算量约降为 1/NKRED³；**金属和掺杂（有部分占据）体系慎用** |
| `HFRCUT` | −1 | 截断库仑核；有带隙体系和能带计算推荐；金属可能 0（默认）收敛更快 |
| `ISTART` | 1 | 读 PBE 的 WAVECAR |

> [!WARNING]
> 杂化泛函下 `ISIF` 默认为 0，不计算应力。用 HSE 做变胞优化必须显式写 `ISIF = 3`（同时注意代价）。另外，`vasp_ncl` 中带 SOC 的杂化泛函计算在 6.6.1 及更早版本存在 PAW 单中心交换项的已知问题（known issue 116），做 HSE+SOC 前请查看 [Known issues](https://www.vasp.at/wiki/index.php/Known_issues) 的最新状态。

### 4.3 运行与监控

```bash
cd 1-hse-scf
cp ../POSCAR ../POTCAR ../0-pbe-scf/WAVECAR .
mpirun -np 16 vasp_std > vasp.out

grep -E "^[A-Z]+ *:" OSZICAR | tail -5             # last electronic steps (line prefix depends on ALGO)
grep "FFT grid for exact exchange" -A 2 OUTCAR       # PRECFOCK grid actually used
grep -A12 "Band structure" OUTCAR                    # gap on the SCF mesh (VASP >= 6.4.3)
```

（OSZICAR 中电子步行以算法缩写加冒号开头，具体缩写随 `ALGO` 不同而不同，上面的正则对各种 ALGO 都适用。）

一个 HSE 自洽通常需要 15–40 个电子步。每一步都比 PBE 慢得多：Si 原胞在几十个核上几分钟内就能完成，而 100 原子的超胞可能需要几百甚至上千核时。

### 4.4 并行与计算量控制

| 手段 | 做法 | 注意 |
|---|---|---|
| k 点并行 | `KPAR` = 节点数或不可约 k 点数的因子 | 杂化泛函通常从 KPAR 中获益明显；内存随 KPAR 增加 |
| NCORE | 按节点核数设置后实测 | 与 PBE 的最佳值可能不同 |
| 降采样 | `NKRED = 2`（或 `NKREDX/Y/Z`） | 先在小体系上比较带隙变化（通常 < 0.05 eV 可接受） |
| 预收敛 | 先 `PRECFOCK = Fast`，再 `Normal` 续算 | 两步之间保留 WAVECAR |
| 减少空带 | NBANDS 只比占据带多需要的数量 | 空带越多 Fock 越贵 |
| GPU | OpenACC 版对杂化泛函加速明显 | VASP 6.6.x 的 NVIDIA GPU 版在杂化泛函 + `ALGO = Normal` 时有随机崩溃的已知问题（known issue 99），改用 `ALGO = Damped` |

---

## 5. HSE06 能带：方法 1——零权重 k 点

示例目录：`2-hse-band-zero-weight`。

### 5.1 思路

把**规则网格的不可约 k 点（带正常权重）**和**能带路径上的 k 点（权重为 0）**写进同一个显式 KPOINTS 文件，做一次完整的 HSE 自洽。零权重点不参与电荷密度和 Fock 算符的构造，但它们的本征值会被正确地对角化出来。

### 5.2 生成 KPOINTS

第一部分取自同一网格的 PBE 计算产生的 `IBZKPT`（k 点 + 权重），第二部分追加路径点（权重 0）。显式 KPOINTS 的结构如下（示意，k 点数与坐标省略）：

```text
Explicit k-points: IBZKPT of the 6x6x6 mesh + zero-weight band path
   N_ibz + N_path
Reciprocal lattice
    0.00000000    0.00000000    0.00000000     1
    0.16666667    0.00000000    0.00000000     8
    ...                                         (weighted k-points copied from IBZKPT)
    0.50000000    0.50000000    0.50000000     0   L
    0.46428571    0.46428571    0.46428571     0
    ...                                         (path points, weight 0)
    0.00000000    0.00000000    0.00000000     0   GAMMA
```

几种生成方法：

```bash
# (a) this repository's example helper (reads IBZKPT + a line-mode path file)
python make_zero_weight_kpoints.py --ibzkpt ../0-pbe-scf/IBZKPT --kpath KPATH --out KPOINTS

# (b) VASPKIT 1.5.x: task 251 generates the KPOINTS for hybrid band structures
#     (asks for k-mesh resolution and path density); task 252 post-processes the bands
vaspkit -task 251

# (c) sumo: an IBZKPT must be present in the current directory
sumo-kgen --hybrid            # option --split N writes several smaller KPOINTS files
```

> [!TIP]
> 路径点越多越贵：所有零权重点在**每一个电子步**都要对角化，而且参与 Fock 交换的计算。每段 10–20 个点通常足够看清能带；需要更细时可以把路径拆成几份分别计算（sumo-kgen `--split`、VASP 官方也建议这样做）。

### 5.3 计算步骤

零权重计算的 k 点集合与 `0-pbe-scf` 不同，不能直接读那里的 WAVECAR，所以在本目录先用**同一个 KPOINTS** 跑一次 PBE，再跑 HSE：

```bash
cd 2-hse-band-zero-weight
cp ../POSCAR ../POTCAR .
python make_zero_weight_kpoints.py --ibzkpt ../0-pbe-scf/IBZKPT --kpath KPATH --out KPOINTS

cp INCAR.pbe INCAR && mpirun -np 16 vasp_std > vasp_pbe.out     # PBE with the same k-point list -> WAVECAR
cp INCAR.hse INCAR && mpirun -np 16 vasp_std > vasp_hse.out     # HSE06 from that WAVECAR
```

`INCAR.hse` 与 §4.2 基本相同，另加 `LORBIT = 11`（fat band）；`HFRCUT = -1` 在这里尤其重要——默认的辅助函数法（`HFRCUT = 0`）在能带路径上会造成 Γ 点附近的不连续。

> [!NOTE]
> VASP 官方提示：零权重法**不要**从已收敛的 HSE WAVECAR 重新开始，因为自洽判据只看总能（零权重点不影响总能），可能在零权重点尚未收敛时就提前停止；从 PBE 波函数开始，零权重点会在整个自洽过程中逐步收敛。如果一定要从 HSE WAVECAR 续算，请把 `NELMIN` 设大一些并仔细检查能带。

### 5.4 后处理

画图时必须**剔除有权重的规则网格点**，只保留零权重的路径点：

- VASPKIT：`vaspkit -task 252`；
- sumo：`sumo-bandplot`（可识别并拼接 `--split` 产生的 split-01、split-02… 目录）；
- pymatgen：`Vasprun("vasprun.xml").get_band_structure(kpoints_filename="KPOINTS", line_mode=True)`，显式 KPOINTS 中高对称点要带标签；
- 本仓库 `../scripts/plot_band.py` 会自动剔除非零权重点（用法见 `python ../scripts/plot_band.py -h`）。

---

## 6. HSE06 能带：方法 2——KPOINTS_OPT（VASP ≥ 6.3）

示例目录：`3-hse-band-kpoints-opt`。

### 6.1 思路

和第5章 PBE 的用法完全一样：`KPOINTS` 写规则网格，`KPOINTS_OPT` 写 Line-mode 路径。VASP 先做 HSE 自洽，再在路径点上做一次性计算；路径点的收敛单独判断（相邻两次迭代的占据带能量之和不再变化）。

```text
SYSTEM   = Si HSE06 SCF + band path (KPOINTS_OPT, VASP >= 6.3)

ISTART   = 1         # PBE WAVECAR of the same uniform mesh (0-pbe-scf)
PREC     = Accurate
ENCUT    = 400
EDIFF    = 1E-6
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
NBANDS   = 16
LORBIT   = 11

LHFCALC  = .TRUE.
GGA      = PE
HFSCREEN = 0.2
AEXX     = 0.25
ALGO     = Damped
TIME     = 0.4
PRECFOCK = Normal
HFRCUT   = -1        # required for unscreened hybrids with KPOINTS_OPT, recommended for HSE as well
# KPOINTS_OPT_NKBATCH = 10   # fewer k-points per batch if memory is tight

NSW      = 0
IBRION   = -1
LWAVE    = .FALSE.
LCHARG   = .FALSE.
```

```bash
cd 3-hse-band-kpoints-opt
cp ../POSCAR ../POTCAR ../0-pbe-scf/WAVECAR .
mpirun -np 16 vasp_std > vasp.out
```

作图（py4vasp）：

```python
import py4vasp
calc = py4vasp.Calculation.from_path("3-hse-band-kpoints-opt")
calc.band.plot("kpoints_opt")
print(calc.bandgap.fundamental())
```

### 6.2 两种方法对比

| | 零权重 k 点 | KPOINTS_OPT |
|---|---|---|
| 版本 | 所有版本 | ≥ 6.3（hybrid 能带的 Coulomb 发散修正问题 issue 80 在 6.3.1 修复，建议 ≥ 6.5） |
| 准备工作 | 需要 IBZKPT + 拼接脚本；另跑一次同 k 集的 PBE | 两个普通 KPOINTS 文件 |
| 路径点的代价 | 每个电子步都算，且参与 Fock | 只在自洽后算一次（分批，`KPOINTS_OPT_NKBATCH`） |
| 内存 | 零权重点多时可能爆内存，需拆分 | 用 NKBATCH 控制 |
| 收敛判断 | 只看总能，零权重点可能未充分收敛 | 路径点单独判断收敛 |
| 后处理 | 需剔除加权点 | py4vasp/pymatgen 直接识别 |

**建议**：VASP ≥ 6.5 时优先用 KPOINTS_OPT；只有旧版本或需要与旧工作流兼容时用零权重法。

> [!IMPORTANT]
> 无论哪种方法，都**不要**设置 `ICHARG = 11`。如果忘了设 `HFRCUT`，能带在 Γ 点附近会出现不连续（半芯态上最明显）；官方说明可以通过扣除半芯态的伪色散部分补救，但最省事的做法是设好 `HFRCUT = -1` 重算。

---

## 7. HSE 的其他常见用法

### 7.1 HSE 结构优化

```text
LHFCALC = .TRUE. ; HFSCREEN = 0.2 ; ALGO = Damped ; TIME = 0.4
IBRION  = 2      ; ISIF = 3        ; NSW = 50     ; EDIFFG = -0.02
PRECFOCK = Normal
```

代价是 PBE 优化的几十倍以上。常见折中：PBE（或 r2SCAN）优化 → HSE 单点；对晶格常数敏感的性质（如带隙随应变变化）再考虑 HSE 优化。

### 7.2 调整 AEXX

对某些材料（如 ZnO、GaN），HSE06 仍明显低估带隙。文献中常把 `AEXX` 调到 0.3–0.4 去拟合实验带隙。这样做时要说明：调参后的泛函不再是 HSE06，其他性质（形成能、缺陷能级）的准确性未必随之改善。更系统的做法是介电依赖杂化泛函（DDH，`LMODELHF`；6.6.0 起支持 RS-DDH），见 [VASP Wiki: List of hybrid functionals](https://www.vasp.at/wiki/index.php/List_of_hybrid_functionals)。

### 7.3 其他杂化泛函

| 泛函 | INCAR |
|---|---|
| PBE0 | `LHFCALC = .TRUE.`、`GGA = PE`（不设 HFSCREEN）；未屏蔽，k 点收敛更慢，能带建议 `HFRCUT = -1` |
| HSE03 | `LHFCALC = .TRUE.`、`GGA = PE`、`HFSCREEN = 0.3` |
| HSE06 | `LHFCALC = .TRUE.`、`GGA = PE`、`HFSCREEN = 0.2` |
| meta-GGA 杂化（如 r2SCAN0 一类） | VASP ≥ 6.4.0 支持 meta-GGA 与 HF 交换混合（`AEXX`、`AMGGAX`、`AMGGAC`），具体写法以官方 List of hybrid functionals 为准 |

---

## 8. r2SCAN：更均衡的 meta-GGA

示例目录：[../examples/06-r2scan-relax-zno/README.md](../examples/06-r2scan-relax-zno/README.md)。

### 8.1 INCAR（纤锌矿 ZnO 变胞优化）

```text
SYSTEM  = ZnO wurtzite, r2SCAN full relaxation

PREC    = Accurate
ENCUT   = 520        # >= 1.3 x ENMAX(O ~400 eV); cell relaxation
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.05
LREAL   = .FALSE.

METAGGA = R2SCAN     # VASP >= 6.2.0; no GGA tag is needed
LASPH   = .TRUE.     # aspherical one-centre terms: strongly recommended for meta-GGA
LMIXTAU = .TRUE.     # pass the kinetic-energy density through the mixer (helps convergence)
ALGO    = Normal     # density mixing; if SCF struggles try ALGO = All (+ ISEARCH = 1 on >= 6.5.1)

IBRION  = 2
ISIF    = 3
NSW     = 60
EDIFFG  = -0.01

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `METAGGA` | R2SCAN | 也可用 SCAN、RSCAN；SCAN 数值上更不稳定 |
| `LASPH` | .TRUE. | VASP 强烈建议 meta-GGA 打开；d/f 元素、强方向性成键尤其重要 |
| `LMIXTAU` | .TRUE. | 默认 .FALSE.；密度混合类算法（ALGO = Normal/Fast）收敛困难时打开 |
| `ALGO` | Normal 或 All | 含真空的体系（slab、分子）中 `ALGO = All` 往往更稳 |
| `ENCUT` | ≥ 1.3×ENMAX | meta-GGA 对基组略敏感，变胞时更要提高 |

POTCAR 检查：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Zn,O}/POTCAR > POTCAR
grep kinetic POTCAR | sort | uniq -c       # every element must provide kinetic-energy-density data
```

### 8.2 色散修正

r2SCAN 只描述了部分中程 vdW 作用，层状/分子晶体仍需色散修正：

**r2SCAN + rVV10**（Ning 等 2022 的参数，非局域 vdW，VASP Wiki 给出的写法）：

```text
METAGGA  = R2SCAN
LUSE_VDW = .TRUE.
BPARAM   = 11.95     # r2SCAN+rVV10 value (default 6.3 is for rVV10 with GGA=ML)
CPARAM   = 0.0093
LASPH    = .TRUE.
```

rVV10 不需要 `vdw_kernel.bindat`。VASP 6.2 及更早版本 rVV10 的应力张量有误，**变胞优化请用 ≥ 6.3.0**。

**r2SCAN + D4**（需要编译时链接 DFT-D4 库，见 [第1章](01-setup-compile.md)）：

```text
METAGGA = R2SCAN
IVDW    = 13         # DFT-D4; VASP passes the functional name to the D4 library
# DFTD4_XC = r2scan # VASP >= 6.6.0: set the functional explicitly if needed (names: param.f90 of dftd4)
```

是否需要 `DFTD4_XC` 取决于版本和泛函：官方说明对 PBE、SCAN、HSE06 等常见泛函只设 `IVDW = 13` 即可；对其他泛函可用 `DFTD4_XC` 指定（6.6.0+）。运行后请在 OUTCAR 中确认 D4 使用的参数确实对应 r2SCAN。更多 vdW 方法见 [第8章](08-vdw-2d.md)。

### 8.3 r2SCAN 的注意事项

- **VASP 6.6.0 的已知 bug（issue 94）**：自旋极化（`ISPIN = 2`）meta-GGA 的受力错误。做磁性体系的 r2SCAN 优化或 MD 请使用 **6.6.1**，或给 6.6.0 打官方 `pot.F.patch` 后重新编译。
- 收敛困难时：先用 PBE 收敛得到 WAVECAR，再以 `ISTART = 1` 启动 r2SCAN；或 `ALGO = All`（6.5.1 起可加 `ISEARCH = 1` 使用改进的线搜索）。
- 能带：与杂化泛函相同，使用 KPOINTS_OPT 或零权重法，**不能** `ICHARG = 11`（VASP Wiki: Band-structure calculation using meta-GGA functionals）。
- 读取 CHGCAR 续算时（`ICHARG = 1`），若目录中有 `TAUCAR`，VASP 也会读入动能密度；否则从空的 τ 开始，前几步可能不稳。

参考量级：纤锌矿 ZnO 的实验晶格常数约 a = 3.25 Å、c = 5.21 Å；PBE 一般高估约 1–2%，r2SCAN 通常更接近实验（具体数值随设置而变，请以你自己的计算为准）。

---

## 9. mBJ：便宜的带隙

示例目录：[../examples/06-mbj-gap-si/README.md](../examples/06-mbj-gap-si/README.md)。

### 9.1 INCAR

```text
SYSTEM  = Si mBJ band gap (potential-only: no total energies, no forces)

ISTART  = 1          # start from the PBE WAVECAR of the previous step
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
NBANDS  = 16

METAGGA = MBJ
LASPH   = .TRUE.
ALGO    = All
TIME    = 0.1        # small step: mBJ converges slowly
NELM    = 200
# CMBJ  = 1.20       # fix c: needed for slabs/interfaces and for restarts (value from a bulk run)
LORBIT  = 11

NSW     = 0
IBRION  = -1
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

示例中 KPOINTS 为规则网格、KPOINTS_OPT 为能带路径，一次运行同时得到 mBJ 的带隙和能带。

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `METAGGA` | MBJ | Tran–Blaha 修正的 BJ 势；`LMBJ` 为局域版本，适合含真空/界面的体系 |
| `LASPH` | .TRUE. | 必须 |
| `ALGO`/`TIME` | All / 0.1 | 官方故障排查页给出的稳妥组合 |
| `NELM` | 100–200 | mBJ 收敛慢 |
| `CMBJ` | 不设（自洽）或固定值 | 不设时 c 每个电子步自洽更新；**slab、界面体系必须固定**（真空区的 $|\nabla\rho|/\rho$ 无意义） |
| `CMBJA`/`CMBJB` | 默认 −0.012 / 1.023 | 一般不改 |

### 9.2 读取 c 与续算

```bash
grep CMBJ OUTCAR | tail -1        # converged self-consistent c (printed every electronic step)
```

mBJ 在 c 自洽时无法从 WAVECAR "接着算"而不重复自洽过程；需要续算（例如另算一条能带路径）时，把上面得到的 c 写进 `CMBJ` 固定下来。官方推荐的难收敛体系流程：PBE 收敛 → 固定 `CMBJ` 的 mBJ（`ALGO = All`、`TIME = 0.1`）→ 不固定 `CMBJ` 的 mBJ，每一步都从上一步的 WAVECAR 开始。

### 9.3 mBJ 能做与不能做的

| 能做 | 不能做 |
|---|---|
| 带隙、能带、DOS、有效质量（定性）、磁矩 | 结构优化、MD（没有对应的能量泛函，受力无意义） |
| 与 SOC 结合（`vasp_ncl`） | 比较总能、形成能、吸附能 |
| 作为高通量带隙筛选的工具 | 未固定 c 时用于 slab/界面 |

---

## 10. 方法比较：带隙与计算代价

下表是**文献/实验量级的参考值**，用于判断你的结果是否落在合理范围，不是本仓库的计算输出；具体数值依赖 POTCAR、结构、k 点与实现细节。

| 体系 | PBE | r2SCAN | HSE06 | mBJ | 实验 |
|---|---|---|---|---|---|
| Si（间接） | ~0.5–0.6 eV | 略大于 PBE | ~1.1–1.3 eV | ~1.1–1.2 eV（Tran & Blaha 2009） | ~1.17 eV（低温） |
| 单层 MoS2（K 点直接） | ~1.6–1.7 eV | 略大于 PBE | ~2.1–2.3 eV | 文献差异较大 | 光学带隙 ~1.8–1.9 eV（Mak 等 2010；准粒子带隙更大，约 2.4 eV 量级，差值为激子结合能） |
| ZnO（纤锌矿，直接） | ~0.7–0.8 eV | 略大于 PBE | ~2.4–2.5 eV | ~2.7 eV（Tran & Blaha 2009） | ~3.4 eV |

经验结论：

- PBE 系统性低估带隙（常低估 40–60%）；r2SCAN 只有小幅改善；
- HSE06 对 sp 半导体（Si、GaAs）通常很好，对宽带隙离子性氧化物（ZnO、MgO）仍偏低；
- mBJ 以接近 meta-GGA 的代价给出接近实验的带隙，但能带宽度、有效质量等细节不一定可靠；
- 二维材料要区分**准粒子带隙**（GW/STS）与**光学带隙**（PL/吸收），HSE06 的值介于两者之间，不能简单与 PL 对比。

相对代价（同一 Si 原胞、同一 k 网格的数量级估计）：

| 方法 | 相对 PBE 的耗时 | 内存 | 能否优化结构 |
|---|---|---|---|
| PBE | 1 | 1 | 能 |
| r2SCAN | 约 2–5 倍 | 略高 | 能 |
| mBJ | 约 2–10 倍（主要因为电子步多） | 略高 | **不能** |
| HSE06 | 约 10–100 倍以上（随 k 点与体系增大） | 高得多 | 能，但很贵 |

---

## 参数选择与收敛建议

| 参数 | HSE06 | r2SCAN | mBJ |
|---|---|---|---|
| 起点 | PBE WAVECAR（同 k 网格、同 NBANDS） | 可从头；困难时用 PBE WAVECAR | PBE WAVECAR |
| `ALGO` | Damped（TIME 0.4）/ All / Normal；**不用 Fast** | Normal + `LMIXTAU`；困难时 All | All（TIME 0.1） |
| `PRECFOCK` | Normal（预收敛 Fast，声子 Accurate） | — | — |
| `NKRED` | 2（半导体，测试后使用）；金属/掺杂慎用 | — | — |
| `HFRCUT` | −1（有带隙体系、能带计算） | — | — |
| `LASPH` | d/f 元素建议 .TRUE. | .TRUE. | .TRUE. |
| `NELM` | 60–100 | 60–100 | 100–200 |
| k 网格 | 收敛测试；可比 PBE 稍粗；用 NKRED 时需可被整除 | 同 PBE | 同 PBE |
| 能带 | KPOINTS_OPT（≥ 6.5 推荐）或零权重 | 同 HSE | 同 HSE（c 自洽时用 KPOINTS_OPT 最方便） |
| 版本 | ≥ 6.3.1（hybrid 能带修正）；HSE+SOC 关注 issue 116 | ≥ 6.2.0；自旋极化避开 6.6.0（issue 94）；rVV10 变胞 ≥ 6.3.0 | 任意 6.x |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| HSE 自洽振荡、不收敛 | 从随机波函数开始；ALGO 不合适；TIME 太大 | 先跑 PBE 得到 WAVECAR；`ALGO = Damped` 并把 `TIME` 降到 0.2–0.3；或 `ALGO = All` |
| 用了 `ALGO = Fast` 的 HSE 结果奇怪或不收敛 | 杂化泛函下 Fast 不被正确支持，且**不会给警告** | 改用 Damped / All / Normal |
| HSE 计算中晶格完全没变 | 杂化泛函默认 `ISIF = 0`，不计算应力 | 显式写 `ISIF = 3` |
| 对 HSE/SCAN 用了 `ICHARG = 11`，得到错误或无法运行 | 杂化/meta-GGA 不能固定电荷密度做非自洽 | 改用零权重 k 点或 KPOINTS_OPT |
| HSE 能带在 Γ 点附近有"台阶"或不连续 | 默认 `HFRCUT = 0`（辅助函数法）不适合路径上的 k 点 | 设 `HFRCUT = -1` 重算 |
| 零权重法读 WAVECAR 后很快"收敛"但能带乱 | 从已收敛的 HSE WAVECAR 开始，零权重点未收敛就停止 | 从同 k 集的 PBE WAVECAR 开始；或增大 `NELMIN` |
| 内存不足、作业被系统杀掉 | 零权重点太多 / NBANDS 太多 / KPAR 太大 | 拆分路径（sumo-kgen `--split`）；KPOINTS_OPT 用 `KPOINTS_OPT_NKBATCH`；`NKRED = 2`；减小 KPAR；多申请节点 |
| `NKRED = 2` 后带隙明显变化（> 0.1 eV）或金属体系占据异常 | 降采样不适用于金属/部分占据体系 | 改回 `NKRED = 1` |
| WAVECAR 读入后从头开始或报错 | k 点集合、NBANDS 或 ENCUT 与生成 WAVECAR 的计算不一致 | 保证 PBE 与 HSE 的 KPOINTS、NBANDS、ENCUT 一致 |
| r2SCAN 第一步就发散 / 报 POTCAR 缺少动能密度信息 | POTCAR 不含芯电子动能密度（如部分 `_GW` 势） | `grep kinetic POTCAR` 检查，换用合适的 PAW |
| r2SCAN 收敛很慢 | 密度混合中未混合 τ | `LMIXTAU = .TRUE.`；或 `ALGO = All`；或从 PBE WAVECAR 开始 |
| VASP 6.6.0 中自旋极化 r2SCAN 的优化/MD 结果异常 | known issue 94 | 升级到 6.6.1 或打 `pot.F.patch` |
| mBJ 在 slab/界面体系中 c 发散 | 真空区 $|\nabla\rho|/\rho$ 无意义 | 固定 `CMBJ`（取体相的值），或改用 `METAGGA = LMBJ` |
| mBJ 优化后结构奇怪 | mBJ 没有能量泛函，受力无意义 | 不要用 mBJ 做优化；结构用 PBE/r2SCAN 优化 |

---

## 练习

1. 对 Si 分别用 `NKRED = 1` 和 `NKRED = 2` 做 HSE06 自洽，比较带隙与耗时。
2. 用零权重法和 KPOINTS_OPT 法分别计算 Si 的 HSE06 能带，叠图比较；再把 KPOINTS_OPT 计算中的 `HFRCUT = -1` 去掉，观察 Γ 点附近的变化。
3. 对 ZnO 计算 PBE、r2SCAN 优化后的晶格常数和 HSE06@PBE 带隙；把 `AEXX` 调到 0.375 重算，观察带隙变化，并讨论这种"调参"的局限。
4. 对 Si 做 mBJ 计算，用 `grep CMBJ OUTCAR | tail -1` 读出 c；固定这个 c 再算一次，确认带隙一致。
5. 估算：一个 64 原子 Si 超胞、2×2×2 网格的 HSE06 单点，相对 PBE 大约要多少核时？用小体系的实测耗时外推，并说明 NKRED、KPAR 的作用。

---

## 参考资料

**官方文档**

- [VASP Wiki: LHFCALC](https://www.vasp.at/wiki/index.php/LHFCALC)、[HFSCREEN](https://www.vasp.at/wiki/index.php/HFSCREEN)、[AEXX](https://www.vasp.at/wiki/index.php/AEXX)、[PRECFOCK](https://www.vasp.at/wiki/index.php/PRECFOCK)、[NKRED](https://www.vasp.at/wiki/index.php/NKRED)、[HFRCUT](https://www.vasp.at/wiki/index.php/HFRCUT)、[TIME](https://www.vasp.at/wiki/index.php/TIME)
- [VASP Wiki: List of hybrid functionals](https://www.vasp.at/wiki/index.php/List_of_hybrid_functionals)
- [VASP Wiki: Downsampling of the Hartree-Fock operator](https://www.vasp.at/wiki/index.php/Downsampling_of_the_Hartree-Fock_operator) — NKRED 的原理与不适用场景
- [VASP Wiki: Band-structure calculation using hybrid functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals) — 零权重法与 KPOINTS_OPT 法
- [VASP Wiki: Band-structure calculation using meta-GGA functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_meta-GGA_functionals)
- [VASP Wiki: KPOINTS_OPT](https://www.vasp.at/wiki/index.php/KPOINTS_OPT)
- [VASP Wiki: METAGGA](https://www.vasp.at/wiki/index.php/METAGGA)、[LMIXTAU](https://www.vasp.at/wiki/index.php/LMIXTAU)、[LASPH](https://www.vasp.at/wiki/index.php/LASPH)、[CMBJ](https://www.vasp.at/wiki/index.php/CMBJ)
- [VASP Wiki: Nonlocal vdW-DF functionals](https://www.vasp.at/wiki/index.php/Nonlocal_vdW-DF_functionals)（r2SCAN+rVV10 参数）、[DFT-D4](https://www.vasp.at/wiki/index.php/DFT-D4)
- [VASP Wiki: Troubleshooting electronic convergence](https://www.vasp.at/wiki/index.php/Troubleshooting_electronic_convergence)（mBJ 分步收敛方案）
- [VASP Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — issue 80、94、99、116
- [VASPKIT 教程](https://vaspkit.com/tutorials.html)（任务 251/252）、[sumo-kgen 文档](https://smtg-bham.github.io/sumo/sumo-kgen.html)、[py4vasp 文档](https://vasp.at/py4vasp/latest/)

**相关 GitHub 仓库（仅链接）**

- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) — `sumo-kgen --hybrid` 生成零权重 KPOINTS，`sumo-bandplot` 自动拼接
- [JMSkelton/VASP-Examples](https://github.com/JMSkelton/VASP-Examples) — Si 的 Wannier90 插值 HSE06 能带、SnS2 的 PBE/HSE06/GW 对比输入
- [chouyoudou/mBJ-VASP-tutorial](https://github.com/chouyoudou/mBJ-VASP-tutorial) — mBJ 带隙计算的小教程
- [QijingZheng/VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding) — 含 HSE k 点生成等脚本
- [dftd4/dftd4](https://github.com/dftd4/dftd4) — IVDW = 13 的外部库（param.f90 中列出泛函名）
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) — 读取 KPOINTS_OPT 能带与带隙
- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90) — 用 Wannier 插值得到平滑的杂化泛函能带（见 [第17章](17-advanced-electronic.md)）

**关键文献**

- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* 118, 8207 (2003). DOI: [10.1063/1.1564060](https://doi.org/10.1063/1.1564060)；Erratum *J. Chem. Phys.* 124, 219906 (2006). DOI: [10.1063/1.2204597](https://doi.org/10.1063/1.2204597)
- A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, G. E. Scuseria, *J. Chem. Phys.* 125, 224106 (2006)（HSE06）. DOI: [10.1063/1.2404663](https://doi.org/10.1063/1.2404663)
- J. Paier et al., *J. Chem. Phys.* 124, 154709 (2006)（VASP 中的屏蔽杂化泛函）. DOI: [10.1063/1.2187006](https://doi.org/10.1063/1.2187006)
- J. Spencer, A. Alavi, *Phys. Rev. B* 77, 193110 (2008)（截断库仑核，HFRCUT = -1）. DOI: [10.1103/PhysRevB.77.193110](https://doi.org/10.1103/PhysRevB.77.193110)
- J. Sun, A. Ruzsinszky, J. P. Perdew, *Phys. Rev. Lett.* 115, 036402 (2015)（SCAN）. DOI: [10.1103/PhysRevLett.115.036402](https://doi.org/10.1103/PhysRevLett.115.036402)
- J. W. Furness et al., *J. Phys. Chem. Lett.* 11, 8208 (2020)（r2SCAN）. DOI: [10.1021/acs.jpclett.0c02405](https://doi.org/10.1021/acs.jpclett.0c02405)
- J. Ning et al., *Phys. Rev. B* 106, 075422 (2022)（r2SCAN+rVV10）. DOI: [10.1103/PhysRevB.106.075422](https://doi.org/10.1103/PhysRevB.106.075422)
- E. Caldeweyher et al., *J. Chem. Phys.* 150, 154122 (2019)（DFT-D4）. DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
- A. D. Becke, E. R. Johnson, *J. Chem. Phys.* 124, 221101 (2006). DOI: [10.1063/1.2213970](https://doi.org/10.1063/1.2213970)
- F. Tran, P. Blaha, *Phys. Rev. Lett.* 102, 226401 (2009)（mBJ）. DOI: [10.1103/PhysRevLett.102.226401](https://doi.org/10.1103/PhysRevLett.102.226401)
- K. F. Mak, C. Lee, J. Hone, J. Shan, T. F. Heinz, *Phys. Rev. Lett.* 105, 136805 (2010)（单层 MoS2 光学带隙）. DOI: [10.1103/PhysRevLett.105.136805](https://doi.org/10.1103/PhysRevLett.105.136805)

---

[← 上一章](05-band-dos.md) | [返回目录](../README.md) | [下一章 →](07-magnetism-u-soc.md)
