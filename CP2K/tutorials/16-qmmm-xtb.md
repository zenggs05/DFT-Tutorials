# 第16章 QM/MM 与半经验方法（xTB/DFTB）（QM/MM and Semi-empirical Methods: GFN-xTB, DFTB）

本章介绍 CP2K 中两类比全 DFT 便宜得多的方法，以及它们的组合。第一类是紧束缚半经验方法：原生 GFN0/GFN1-xTB、通过 tblite 接入的 GFN2-xTB（支持 k 点，UKS 下即 spGFN2-xTB），以及 SCC-DFTB。第二类是 QM/MM：用 `&FORCE_EVAL METHOD QMMM` 把反应中心交给 DFT 或 xTB，其余环境用 AMBER、CHARMM 或输入文件里自定义的力场描述。文中依次讲解 `&QMMM`（`&QM_KIND MM_INDEX`、`E_COUPL GAUSS/COULOMB`、`&CELL`、`CENTER`、`&LINK`、`&PERIODIC`）、`&MM &FORCEFIELD PARMTYPE AMBER` 与 `&TOPOLOGY CONN_FILE_FORMAT AMBER` 的写法，以及 QM 区选择、link atom 电荷修正、周期性 QM/MM、QM/MM 下的 NEB/Dimer 过渡态搜索和 GROMACS–CP2K 接口。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 DFTB、GFN0/GFN1/GFN2-xTB 的能量表达式差异，知道它们各自适合做什么、不适合做什么；
- 按 CP2K 版本写出正确的 xTB 输入：2026.2 用 `GFN_TYPE TBLITE`，2025.2/2026.1 只需写 `&TBLITE`，周期性原生 xTB/DFTB 在 2026.1 及以前必须写 `DO_EWALD T`；
- 用 tblite GFN2-xTB 做带 k 点的周期性计算，用 UKS 做开壳层 spGFN2-xTB 计算；
- 用 CP2K 自带的 Slater–Koster 参数做 SCC-DFTB 分子动力学；
- 搭出完整的 QM/MM 输入：QM 原子选择、QM 盒子、静电耦合方式、MM 力场（输入文件内定义或 AMBER prmtop）；
- 在切断共价键时正确设置 `&LINK`，并处理边界电荷；
- 理解周期性 QM/MM（`&QMMM &PERIODIC`）的作用与代价；
- 按“MM 平衡 → QM/MM 优化 → CI-NEB → 部分 Hessian → Dimer”的流程在酶或溶液体系中搜索过渡态；
- 知道什么时候改用 GROMACS–CP2K 接口，以及该接口的使用限制。

| 方法 | 相对成本（同一体系） | 适用场景 | 不适用 |
|---|---|---|---|
| GFN1/GFN2-xTB | DFT 的 1/100–1/1000 | 构象搜索、预优化、长时间 MD 预采样、大体系筛选、给 DFT 提供初始结构 | 精确反应能垒、带隙、金属电子结构 |
| GFN0-xTB | 比 GFN1 更低（非自洽） | 超大体系的快速几何预处理 | 电荷转移、极化明显的体系 |
| SCC-DFTB | 与 xTB 相近 | 已有可靠 Slater–Koster 参数的元素组合（CP2K 自带 C H N O P S Zn） | 参数缺失的元素组合 |
| DFT/MM（GPW + 力场） | 只与 QM 区大小有关 | 溶液或酶中的反应、光谱、局域电子结构 | QM 区需要几百个原子以上时，单步成本很高 |
| xTB/MM、DFTB/MM | 非常低 | 长时间 QM/MM 采样、DFT/MM 之前的预平衡 | 定量能垒（应当用 DFT/MM 精修） |

本章示例都在 [../examples/16-qmmm-xtb/README.md](../examples/16-qmmm-xtb/README.md)：

| 文件 | 内容 | 对应小节 |
|---|---|---|
| `ethanol_gfn1_native.inp` | 原生 GFN1-xTB 优化孤立乙醇 | §4 |
| `ice_gfn2_tblite_kpoints.inp` | tblite GFN2-xTB + 4×4×2 k 点，冰 Ih | §5 |
| `o2_triplet_spgfn2.inp` | 三重态 O2，UKS 下的 spGFN2-xTB | §5.4 |
| `dftb_md.inp` | 125 个水分子的 SCC-DFTB NVT MD | §6 |
| `water_qmmm_tip3p.inp` | 1 个 PBE 水分子 + 124 个 MM 水分子，周期性 QM/MM MD | §8 |
| `water_qmmm_amber.inp` + `tleap_water.in` | 同一体系改读 AMBER prmtop | §9 |
| `xtb_qmmm_water.inp` | GFN1-xTB/MM，`E_COUPL COULOMB` | §10 |
| `gromacs_qmmm_fragment.mdp` | GROMACS–CP2K 接口的 mdp 选项片段 | §14 |

> [!NOTE]
> 本章所有示例都按 CP2K 2026.2 编写，并在输入文件头部注明了与 2025.x/2026.1 的差异。tblite 相关示例要求编译时启用 tblite，用 `cp2k.psmp -v | grep -o tblite` 检查。仓库不附带参考输出，文中出现的输出片段都标注为示意。

---

## 2. 原理速览

### 2.1 紧束缚 DFT：从 Kohn–Sham 到 DFTB

DFTB 把 Kohn–Sham 总能量在一个参考密度 $\rho_0$（中性原子密度的叠加）附近做泰勒展开。保留到二阶就是 SCC-DFTB（也叫 DFTB2）：

$$
E_{\mathrm{DFTB2}}=\sum_i^{\mathrm{occ}} n_i\langle\psi_i|\hat H^0|\psi_i\rangle
+\frac12\sum_{A,B}\gamma_{AB}\,\Delta q_A\,\Delta q_B
+E_{\mathrm{rep}}
$$

- 第一项：零阶哈密顿量 $\hat H^0$ 在最小价层基组中的矩阵元，事先按原子对和距离制成 Slater–Koster 表（`.spl` 文件）；
- 第二项：Mulliken 电荷涨落 $\Delta q$ 之间的屏蔽库仑相互作用，需要自洽求解，"SCC" 指的就是这一步；
- 第三项：原子对排斥势，同样逐原子对拟合。

由于积分全部预先制表，DFTB 的计算量集中在对角化一个很小的矩阵上，单步比 GGA-DFT 快两到三个数量级。代价是可移植性取决于参数：没有 Slater–Koster 表的元素组合根本算不了，参数也只在拟合时覆盖的化学环境内可靠。

### 2.2 GFN-xTB 家族

Grimme 组的 GFN-xTB 用一套覆盖 Z = 1–86 的全局参数取代逐原子对的表格。官方文档给出的 GFN1 与 GFN2 能量表达式为：

$$
E_{\mathrm{GFN1}}=E_{\mathrm{EL}}+E_{\mathrm{IES}}+E_{\mathrm{REP}}+E_{\mathrm{DISP}}+E_{\mathrm{XB}},\qquad
E_{\mathrm{GFN2}}=E_{\mathrm{EL}}+E_{\mathrm{AES}}+E_{\mathrm{REP}}+E_{\mathrm{DISP}}
$$

其中 IES 是各向同性静电（含二阶自洽项和三阶对角项 $\frac13\sum_A\Gamma_A q_A^3$），AES 是 GFN2 增加的原子偶极、四极各向异性静电，XB 是 GFN1 的卤键校正。

| 方法 | 自洽 | 静电 | 色散 | CP2K 实现 |
|---|---|---|---|---|
| GFN0-xTB | 否（EEQ 电荷） | 各向同性 | D4 | 原生，`GFN_TYPE 0`（≥ 2025.1） |
| GFN1-xTB | 是 | 各向同性 + 三阶 | D3(BJ) | 原生（默认 `GFN_TYPE 1`），也可通过 tblite |
| IPEA1-xTB | 是 | 同 GFN1 | D3(BJ) | tblite，`METHOD IPEA1` |
| GFN2-xTB | 是 | 加偶极/四极 | 自洽 D4 | 只能通过 tblite，`METHOD GFN2` |
| spGFN2-xTB | 是 | 同 GFN2，另加自旋极化项 | 自洽 D4 | tblite + `UKS`（`LSD`） |

> [!TIP]
> xTB 的设计目标是**几何结构、振动频率和非共价相互作用**（GFN 即 Geometry, Frequency, Noncovalent）。它非常适合给 DFT 提供起点，但不适合直接给出反应能、带隙或磁性能量差。常见用法是：xTB 构象搜索或 MD 预采样 → DFT 重新优化 → DFT 单点或自由能计算（见 [第15章](15-metadynamics-plumed.md) 的“xTB 预采样，DFT 精修”）。

### 2.3 QM/MM 的能量分解

CP2K 采用加和式（additive）QM/MM：

$$
E_{\mathrm{tot}}=E_{\mathrm{QM}}(\mathrm{QM})+E_{\mathrm{MM}}(\mathrm{MM}+\mathrm{QM})+E^{\mathrm{el}}_{\mathrm{QM/MM}}
$$

- $E_{\mathrm{QM}}$：Quickstep（GPW、xTB、DFTB 等）只对 QM 原子（加上 link 氢原子）计算；
- $E_{\mathrm{MM}}$：力场对**全体**原子计算，但 QM 原子的 MM 电荷被置零，纯 QM 内部的成键项被去掉。QM–MM 之间的范德华和跨边界的成键项仍由力场负责；
- $E^{\mathrm{el}}_{\mathrm{QM/MM}}$：MM 点电荷与 QM 电子密度、原子核之间的静电作用，由 `E_COUPL` 决定：

| `E_COUPL` | 嵌入方式 | 适用的 QM 方法 |
|---|---|---|
| `NONE`（**默认**） | 机械嵌入：QM 区感受不到 MM 电荷，不会被极化 | 任意 |
| `GAUSS` | 静电嵌入，GEEP 高斯展开 | GPW/GAPW（2026.2 文档注明不适用于 DFTB） |
| `S-WAVE` | GEEP 的 s 波变体 | GPW/GAPW |
| `COULOMB` | 静电嵌入，解析 $1/r$ 势 | 半经验、DFTB、xTB（文档注明不适用于 GPW/GAPW） |
| `POINT_CHARGE` | QM 区用推导出的点电荷与 MM 电荷作用 | 紧束缚方法（见官方回归测试） |

> [!WARNING]
> `E_COUPL` 的默认值是 `NONE`。忘了写这一行，程序照样能跑，但结果只是机械嵌入，溶剂和蛋白的静电环境对 QM 区完全没有作用。

### 2.4 GEEP：为什么 GPW 能高效地做静电嵌入

在 GPW 中，QM 电子密度定义在实空间格点上。若直接用 $q_J/r$ 计算每个 MM 电荷对每个格点的势，代价是 $N_{\mathrm{MM}}\times N_{\mathrm{grid}}$，而且 $r\to0$ 处发散。GEEP（Gaussian Expansion of the Electrostatic Potential）先把 MM 电荷“抹开”，再用一组宽度不同的高斯函数展开其势：

$$
v_J(\mathbf r)=q_J\,\frac{\mathrm{erf}\!\left(|\mathbf r-\mathbf R_J|/r_{c,J}\right)}{|\mathbf r-\mathbf R_J|}
\approx q_J\sum_{g}A_g\exp\!\left[-\left(\frac{|\mathbf r-\mathbf R_J|}{G_g}\right)^2\right]+R_{\mathrm{low}}(|\mathbf r-\mathbf R_J|)
$$

$r_{c,J}$ 就是 `&MM_KIND` 中的 `RADIUS`，用于消除短程发散，防止 QM 电子“溢出”到正电荷 MM 原子上。窄的高斯函数在细格点上收集，宽的高斯函数在粗格点上收集，再插值回来，因此总成本近似线性。`USE_GEEP_LIB N` 表示使用内置的 N 项展开（2–15），不需要外部的 `MM_POTENTIAL` 文件。

### 2.5 切断共价键：link atom（IMOMM）

QM/MM 边界穿过共价键时，QM 区出现悬挂键。CP2K 默认的 `LINK_TYPE IMOMM` 在 QM 原子与 MM 原子的连线上放一个封端氢原子：

$$
\mathbf r_{\mathrm H}=\mathbf r_{\mathrm{QM}}+\frac{1}{\alpha}\left(\mathbf r_{\mathrm{MM}}-\mathbf r_{\mathrm{QM}}\right),\qquad
\alpha=\frac{r_{\mathrm{eq}}(\mathrm{QM{-}MM})}{r_{\mathrm{eq}}(\mathrm{QM{-}H})}
$$

$\alpha$ 对应 `ALPHA_IMOMM`（默认 1.38）。作用在封端氢上的力按杠杆关系分配给 QM 原子和 MM 原子，所以封端氢不是独立的自由度。另外两种 `LINK_TYPE` 是 `GHO`（广义杂化轨道）和 `PSEUDO`（单价赝势）。

### 2.6 周期性 QM/MM

MM 体系通常是周期性的（溶剂盒子、晶体），而 QM 区只占一个小盒子。Laino 等人（2006）的做法是：QM 密度与**所有** MM 周期映像的长程静电作用用 Ewald 方法处理（`&QMMM &PERIODIC`）；QM 盒子本身的周期映像之间的虚假相互作用，则用密度拟合得到的原子点电荷做多极去耦（`&MULTIPOLE`，默认开启）。这样 QM 盒子不必和 MM 盒子一样大。

---

## 3. CP2K 中的实现与版本差异

| 功能 | 2024.x 及更早 | 2025.1 | 2025.2 / 2026.1 | 2026.2 |
|---|---|---|---|---|
| 原生 GFN1-xTB | 有（不写 `GFN_TYPE`） | 有 | 有 | 有 |
| 原生 GFN0-xTB | 无 | `GFN_TYPE 0`（整数） | `GFN_TYPE 0` | `GFN_TYPE 0` |
| tblite（GFN1/GFN2/IPEA1） | 无 | 无 | 只写 `&TBLITE` 段即启用 | 必须 `GFN_TYPE TBLITE` + `&TBLITE` |
| `&TBLITE` 中的 `ACCURACY` | — | — | 无 | 有（默认 1.0） |
| `SCC_MIXER`（AUTO/TBLITE/CP2K/NONE） | — | — | 无 | 有（默认 AUTO） |
| 周期性 xTB/DFTB 的 Ewald | 必须写 `DO_EWALD T` | 同左 | 同左 | 由 `&CELL PERIODIC` 自动决定，`DO_EWALD` 标为 deprecated |
| QM/MM：xTB/DFTB 用 `GAUSS` 耦合 | 不可用 | 不可用 | 不可用 | 不建议：文档注明 `GAUSS` 不适用于 DFTB，xTB/DFTB 的 QM/MM 回归测试只用 `COULOMB`/`POINT_CHARGE`（trunk 文档已扩展到 DFTB/xTB，需配合 `NOCOMPATIBILITY`） |

> [!IMPORTANT]
> 2026.2 中 `GFN_TYPE` 的合法值是 `0 | 1 | TBLITE`；2025.1–2026.1 中它是整数。因此 `GFN_TYPE TBLITE` 拿到旧版本会在解析阶段报错。反过来，2026.2 的 `&TBLITE` 段只是“选项集合”，真正切换后端的是 `GFN_TYPE TBLITE`。示例 `ice_gfn2_tblite_kpoints.inp` 用预处理变量同时兼容两种写法（见 §5）。

**参数文件**：原生 xTB 读 `xTB0_parameters`、`xTB1_parameters`，DFTB 读 `DFTB/scc/scc_parameter` 及其中列出的 `.spl` 表、`DFTB/uff_table`，GFN1 的 D3 色散读 `dftd3.dat`，这些都在 CP2K 源码的 `data/` 目录下。tblite 的参数内置在库中。运行前设置好 `CP2K_DATA_DIR`：

```bash
export CP2K_DATA_DIR=/path/to/cp2k/data
ls "$CP2K_DATA_DIR"/xTB1_parameters "$CP2K_DATA_DIR"/DFTB/scc/scc_parameter "$CP2K_DATA_DIR"/DFTB/uff_table
cp2k.psmp -v | grep -i -o tblite         # 无输出说明这个可执行文件没有编译 tblite
```

conda-forge 和官方 Docker 镜像是否带 tblite 取决于具体构建，请以 `cp2k.psmp -v` 的输出为准；从源码构建时 tblite 是可选依赖，见 [第1章](01-install.md)。

---

## 4. 原生 GFN1-xTB：孤立乙醇的几何优化

先从最简单的情况入手：原生 GFN1-xTB 不需要基组、赝势和格点截断，`&DFT` 中只有 `&QS METHOD XTB` 和 SCF 设置。

```fortran
&GLOBAL
  PROJECT     ethanol_gfn1
  RUN_TYPE    GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER  200
    MAX_FORCE 4.5E-4               ! hartree/bohr (CP2K default)
  &END GEO_OPT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    CHARGE       0
    MULTIPLICITY 1
    &QS
      METHOD      XTB
      EPS_DEFAULT 1.0E-10
      &XTB
        CHECK_ATOMIC_CHARGES   T   ! abort if Mulliken charges leave the chemical range
        USE_HALOGEN_CORRECTION T   ! XB term (default T, no effect without halogens)
      &END XTB
    &END QS
    &SCF
      SCF_GUESS MOPAC              ! cheap semi-empirical starting guess
      EPS_SCF   1.0E-7
      MAX_SCF   50
      &OT ON
        MINIMIZER      DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC 15.0 15.0 15.0
      PERIODIC NONE                ! isolated molecule: no Ewald sum for xTB
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   ethanol.xyz
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `&QS METHOD` | `XTB` | 选择 xTB；不写 `GFN_TYPE` 即原生 GFN1，兼容 2023.x–2026.2 |
| `GFN_TYPE` | `1`（默认）/ `0` / `TBLITE` | `0` 为原生 GFN0（≥ 2025.1），`TBLITE` 为 tblite 后端（2026.2） |
| `CHECK_ATOMIC_CHARGES` | `T`（默认） | Mulliken 电荷超出化学合理范围时中止，可防止三阶项数值失稳 |
| `USE_HALOGEN_CORRECTION` | `T`（默认） | GFN1 的卤键项 $E_{\mathrm{XB}}$ |
| `DO_NONBONDED` | `F`（默认） | 打开后可加 CP2K 特有的非键校正势，结果不再是标准 GFN1 |
| `EPS_DEFAULT` | 1E-10 – 1E-12 | 积分筛选阈值；做振动分析时取 1E-12 |
| `SCF_GUESS` | `MOPAC` | 半经验初猜，对 xTB 比 `ATOMIC` 更稳 |
| `&OT` | `DIIS` + `FULL_SINGLE_INVERSE` | 有带隙的分子体系用 OT 最快；金属或 k 点体系改用对角化 |
| `&CELL PERIODIC` | `NONE` | 孤立分子；2026.1 及以前的周期性原生 xTB 必须另写 `DO_EWALD T` |

运行与检查：

```bash
cd CP2K/examples/16-qmmm-xtb
cp2k.psmp -i ethanol_gfn1_native.inp -o ethanol_gfn1_native.out
grep "GEOMETRY OPTIMIZATION COMPLETED" ethanol_gfn1_native.out
grep "ENERGY| Total FORCE_EVAL" ethanol_gfn1_native.out | tail -1
grep -c "SCF run converged" ethanol_gfn1_native.out          # 每个优化步一次
grep -i "xtb" ethanol_gfn1_native.out | head -20             # 方法与参数信息
```

成功标准：几十步内出现 `GEOMETRY OPTIMIZATION COMPLETED`；`ethanol_gfn1-pos-1.xyz` 的最后一帧中 C–O、C–C、O–H 键长与常识一致（约 1.4、1.5、0.97 Å）。优化轨迹可以用 `python ../scripts/plot_opt.py -h` 查看用法后作图。

> [!TIP]
> 想用 GFN0 只需在 `&XTB` 中加一行 `GFN_TYPE 0`（≥ 2025.1）。GFN0 不自洽，`SCC_MIXER AUTO` 会被当作 `NONE`。周期性的原生 GFN1，在 2026.1 及以前要写 `DO_EWALD T`，同时 `&DFT &POISSON` 使用周期性求解器；2026.2 起只需 `&CELL PERIODIC XYZ`。官方文档的原生 GFN1 示例就是一个 `PERIODIC XYZ` 加 `DO_EWALD T` 的写法。

---

## 5. tblite GFN2-xTB 与 k 点：冰 Ih

### 5.1 为什么需要 tblite

GFN2-xTB 只能通过 tblite 库调用。官方文档说明 tblite 后端完全支持 k 点，所以小晶胞的分子晶体、层状材料也能用 GFN2 计算，不必为了 Γ 点近似去扩大超胞。示例是 4 个水分子的六方冰晶胞（质子有序模型，O 位置按理想晶格放置，H 按冰规则放置，只是一个起始模型）。

```fortran
&GLOBAL
  PROJECT     ice_gfn2_k442
  RUN_TYPE    ${RUN-ENERGY_FORCE}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&MOTION
  &GEO_OPT
    OPTIMIZER LBFGS
    MAX_ITER  200
  &END GEO_OPT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &QS
      METHOD      XTB
      EPS_DEFAULT 1.0E-12
      &XTB
@IF ${NEWTB-1} == 1
        GFN_TYPE TBLITE
@ENDIF
        &TBLITE
          METHOD   GFN2            ! GFN1 | GFN2 | IPEA1
@IF ${NEWTB-1} == 1
          ACCURACY 0.1             ! 2026.2+: smaller = tighter tblite thresholds (default 1.0)
@ENDIF
        &END TBLITE
      &END XTB
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 2
    &END KPOINTS
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF   1.0E-8
      MAX_SCF   200
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA  0.2
      &END MIXING
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC              4.50 4.50 7.33
      ALPHA_BETA_GAMMA 90.0 90.0 120.0
      PERIODIC         XYZ
    &END CELL
    &COORD
      SCALED
      O    0.33333333   0.66666667   0.06250000
      H    0.09618077   0.54809038   0.01803389
      H    0.45190962   0.90381923   0.01803389
      O    0.66666667   0.33333333   0.56250000
      H    0.54809038   0.09618077   0.51803389
      H    0.90381923   0.45190962   0.51803389
      O    0.66666667   0.33333333   0.93750000
      H    0.54809038   0.45190962   0.98196611
      H    0.66666667   0.33333333   0.80380286
      O    0.33333333   0.66666667   0.43750000
      H    0.33333333   0.66666667   0.30380286
      H    0.45190962   0.54809038   0.48196611
    &END COORD
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

### 5.2 关键字说明

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `GFN_TYPE` | `TBLITE` | 2026.2 中切换到 tblite 后端；此时必须有 `&TBLITE` 段 |
| `&TBLITE METHOD` | `GFN2` | 可选 `GFN1`、`GFN2`（默认）、`IPEA1` |
| `ACCURACY` | 1.0（默认）；需要精确的力、应力时取 0.1 | 缩放 tblite 内部的各种数值阈值，越小越严；使用 tblite 自带的 SCC 混合器时，它也是 SCC 的收敛判据（对应 tblite 命令行的 `--acc`） |
| `PARAM` | 不写 | 读入自定义的 tblite 参数文件（TOML 格式），只在重新参数化时使用 |
| `SCC_MIXER` | `AUTO`（默认） | 2026.2 新增。`AUTO` 对 tblite 使用 tblite 自己的 Broyden 混合器；`CP2K` 改用 CP2K 的电荷混合，这时 SCC 收敛跟随 `EPS_SCF` |
| `&KPOINTS SCHEME` | `MONKHORST-PACK n1 n2 n3` | 与 DFT 相同；k 点计算一律用对角化（2026.2 中 OT 不支持 k 点） |
| `EPS_SCF`、`&MIXING` | 1E-8；`DIRECT_P_MIXING` | 2026.2 默认的 tblite 混合器下作用有限；`SCC_MIXER CP2K` 或 2025.2/2026.1 下起作用 |
| `${RUN-ENERGY_FORCE}` | 预处理变量 | `-E RUN=GEO_OPT` 即可改为结构优化 |
| `${NEWTB-1}` | 预处理变量 | `-E NEWTB=0` 去掉 `GFN_TYPE TBLITE` 和 `ACCURACY`，适用于 2025.2/2026.1 |

### 5.3 运行与检查

```bash
cp2k.psmp -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out                    # 2026.2
cp2k.psmp -E NEWTB=0 -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out         # 2025.2 / 2026.1
cp2k.psmp -E RUN=GEO_OPT -i ice_gfn2_tblite_kpoints.inp -o ice_opt.out    # 结构优化（建议 >= 2026.1）

grep "BRILLOUIN|" ice_sp.out | head                 # k 点方案与 k 点数，确认 k 点生效
grep "ENERGY| Total FORCE_EVAL" ice_sp.out
grep -A 16 "ATOMIC FORCES in" ice_sp.out            # 力：H 原子通常最大，因为 H 位置是按规则放的
```

k 点收敛的做法与 DFT 相同：依次取 2×2×1、3×3×2、4×4×2、6×6×3，看每个水分子的能量变化是否小于 1 meV。可以用 `python ../scripts/conv_kpoints.py -h` 查看自动化脚本的用法。

> [!NOTE]
> 周期性 tblite 的梯度和应力修复进入了 2026.1，所以用 tblite 做 `GEO_OPT` 或 `CELL_OPT` 请用 2026.1 或更新的版本（见示例文件头部注释）。2026.2 的晶胞优化写法（已删除 `&CELL_OPT TYPE`）见 [第7章](07-geo-cell-opt.md)。

### 5.4 开壳层：UKS 下的 spGFN2-xTB

官方文档说明：在 `&DFT` 中打开 `LSD`（即 `UKS`）后，tblite 会自动使用带自旋极化项的 spGFN2-xTB（Neugebauer 等，2023）。三重态 O2 示例：

```fortran
&GLOBAL
  PROJECT     o2_triplet_spgfn2
  RUN_TYPE    GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER  100
  &END GEO_OPT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    UKS                            ! = LSD: spin-polarised (spGFN2-xTB in tblite)
    MULTIPLICITY 3                 ! triplet ground state of O2
    CHARGE       0
    &QS
      METHOD      XTB
      EPS_DEFAULT 1.0E-12
      &XTB
        GFN_TYPE TBLITE
        &TBLITE
          METHOD   GFN2
          ACCURACY 0.1
        &END TBLITE
      &END XTB
    &END QS
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF   1.0E-8
      MAX_SCF   200
      ADDED_MOS 2 2                ! a few virtual orbitals per spin channel
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC      15.0 15.0 15.0
      PERIODIC NONE
    &END CELL
    &COORD
      O   0.000000   0.000000   0.000000
      O   1.230000   0.000000   0.000000
    &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `UKS`（别名 `LSD`、`SPIN_POLARIZED`） | 开壳层时必写 | tblite 使用自旋极化哈密顿量，即 spGFN2-xTB |
| `MULTIPLICITY` | 2S+1 | 决定未配对电子数；O2 基态为 3 |
| `ADDED_MOS` | 每个自旋通道几个 | 对角化时额外求解的空轨道，便于检查 HOMO–LUMO |

```bash
cp2k.psmp -i o2_triplet_spgfn2.inp -o o2_triplet_spgfn2.out
grep "GEOMETRY OPTIMIZATION COMPLETED" o2_triplet_spgfn2.out
grep -i "multiplicity" o2_triplet_spgfn2.out | head -3
```

> [!CAUTION]
> spGFN2 适合快速筛选自旋态的相对顺序（原文献针对过渡金属配合物做了高通量测试），但自旋态能量差的定量结论仍应以 DFT 或更高级的方法为准。该示例只按 2026.2 编写，没有为旧版本提供兼容开关。

---

## 6. SCC-DFTB：125 个水分子的 NVT MD

CP2K 自带一套 Slater–Koster 参数（`data/DFTB/scc/`，覆盖 C H N O P S Zn 之间的原子对），可以直接用于有机分子和水。

```fortran
&GLOBAL
  PROJECT     water125_dftb_nvt
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE    NVT
    STEPS       2000
    TIMESTEP    [fs] 0.5
    TEMPERATURE [K] 300.0
    &THERMOSTAT
      TYPE   CSVR
      REGION GLOBAL
      &CSVR
        TIMECON [fs] 50.0          ! short coupling for equilibration; use ~100-500 fs for production
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &VELOCITIES OFF
    &END VELOCITIES
    &RESTART
      BACKUP_COPIES 0
      &EACH
        MD 100
      &END EACH
    &END RESTART
    &RESTART_HISTORY OFF
    &END RESTART_HISTORY
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &QS
      METHOD DFTB
      &DFTB
        SELF_CONSISTENT T          ! SCC-DFTB (DFTB2)
        DISPERSION      T          ! UFF pair dispersion; needs UFF_FORCE_FIELD below
        DO_EWALD        T          ! periodic electrostatics, see version notes above
        &PARAMETER
          PARAM_FILE_PATH DFTB/scc
          PARAM_FILE_NAME scc_parameter
          UFF_FORCE_FIELD ../uff_table
        &END PARAMETER
      &END DFTB
    &END QS
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA      0.35
        GMAX       32
        O_SPLINE   5
      &END EWALD
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   50
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA  0.2
      &END MIXING
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC      15.54 15.54 15.54
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   water125.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `&QS METHOD` | `DFTB` | 选择 DFTB |
| `SELF_CONSISTENT` | `T`（默认） | SCC-DFTB；`F` 为非自洽 DFTB，此时应改用 `DFTB/nonscc` 中的参数 |
| `DISPERSION` | `T` | 打开色散校正，类型由 `&PARAMETER DISPERSION_TYPE` 决定（`UFF` 默认，另有 `D2`、`D3`、`D3(BJ)`） |
| `DO_EWALD` | ≤ 2026.1：周期体系必写 `T` | 2026.2 起由 `&CELL PERIODIC` 自动决定，此关键字被标为 deprecated（仍可作调试开关）；全部机器升级到 2026.2 后可删去 |
| `HB_SR_GAMMA` | 视参数而定 | 针对氢键修正的 γ 函数，指数由 `HB_SR_PARAM`（默认 4.0）给出 |
| `DIAGONAL_DFTB3` | 默认 `F` | 对角形式的三阶校正（DFTB3）；需要与之匹配的参数 |
| `PARAM_FILE_PATH` | `DFTB/scc` | 参数目录（相对 `CP2K_DATA_DIR`） |
| `PARAM_FILE_NAME` | `scc_parameter` | 文本文件，每行格式为 `ATOM1 ATOM2 filename.spl` |
| `SK_FILE` | 可重复 | 不用列表文件，直接写 `SK_FILE O H oh.spl` 这样的原子对 |
| `UFF_FORCE_FIELD` | `../uff_table` | `DISPERSION T` 且类型为 UFF 时**必须**给出，否则手册明确写到会段错误 |
| `&POISSON &EWALD` | `SPME` | 周期性 SCC 静电项的 Ewald 参数 |

运行与分析：

```bash
mpirun -np 8 cp2k.psmp -i dftb_md.inp -o dftb_md.out
tail -3 water125_dftb_nvt-1.ener        # 列：Step  Time[fs]  Kin[a.u.]  Temp[K]  Pot[a.u.]  Cons Qty[a.u.]  UsedTime[s]
grep -c "SCF run converged" dftb_md.out
python ../scripts/md_analyze.py -h      # 查看 RDF、温度/能量曲线等分析脚本的用法
```

成功标准：温度在 300 K 附近波动；守恒量（Cons Qty）在平衡后没有系统性漂移；O–O 径向分布函数的第一峰位于 2.7–2.9 Å 附近（具体峰位和峰高取决于 DFTB 参数，标准 SCC-DFTB 已知会给出过于结构化的液态水）。`water125.xyz` 是理想化的格点初始结构，前约 1 ps 必须舍弃。AIMD 的一般分析方法见 [第12章](12-aimd.md)。

---

## 7. QM/MM 输入结构总览

QM/MM 计算在一个 `&FORCE_EVAL` 里同时包含四部分：

```text
&FORCE_EVAL
  METHOD QMMM
  &DFT ... &END DFT          QM part: GPW / xTB / DFTB, CHARGE and MULTIPLICITY of the QM region only
  &MM                        MM part: force field for ALL atoms + MM Ewald
    &FORCEFIELD ... &END FORCEFIELD
    &POISSON ... &END POISSON
  &END MM
  &QMMM                      coupling
    E_COUPL GAUSS | COULOMB | ...
    &CELL ... &END CELL      QM box (smaller than the MM box)
    &QM_KIND O               which atoms are QM, grouped by element
      MM_INDEX 1
    &END QM_KIND
    &MM_KIND O               Gaussian width of the MM charge (GEEP)
      RADIUS 0.78
    &END MM_KIND
    &LINK ... &END LINK      only if a covalent bond is cut
    &PERIODIC ... &END PERIODIC   periodic QM/MM electrostatics
  &END QMMM
  &SUBSYS                    whole system: MM cell, coordinates, topology, &KIND for QM basis/pseudo
  &END SUBSYS
&END FORCE_EVAL
```

几个容易混淆的点：

- **原子编号**：`MM_INDEX` 是原子在整个体系（坐标文件或拓扑）中的序号，从 1 开始；`&QM_KIND` 的名字是元素符号，同一元素的所有 QM 原子写在同一段里，可以分多行写多个 `MM_INDEX`；
- **两个 `&CELL`**：`&SUBSYS &CELL` 是 MM 盒子（真实周期）；`&QMMM &CELL` 是 QM 盒子，GPW 的实空间格点只铺在这个盒子里，它决定了 QM 部分的成本；
- **两个电荷**：`&DFT CHARGE` 只是 QM 区的净电荷；MM 部分的电荷来自力场。切断共价键后，两部分加起来不一定等于体系总电荷（见 §11）；
- **`CENTER`**：默认 `EVERY_STEP`，每一步把 QM 盒子重新居中到 QM 原子上；`SETUP_ONLY` 只在开始时居中，`NEVER` 不居中。`CENTER_TYPE PBC_AWARE_MAX_MINUS_MIN` 在 QM 区跨越 MM 盒子边界时有用；
- **QM 原子的 `&KIND`**：GPW 时要在 `&SUBSYS` 中为 QM 元素给出 `BASIS_SET`、`POTENTIAL`；xTB/DFTB 不需要。

---

## 8. 实例：DFT/MM 水溶液（输入文件内定义 TIP3P 类力场）

第一个 QM/MM 例子尽量不依赖外部文件：MM 参数全部写在输入文件里（`PARMTYPE OFF`），QM 区是 1 个水分子（原子 1–3），用 PBE/DZVP 描述；其余 124 个水分子是柔性的 TIP3P 类模型。

```fortran
&GLOBAL
  PROJECT     water_qmmm_tip3p
  RUN_TYPE    MD
  PRINT_LEVEL LOW
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE    NVT
    STEPS       200
    TIMESTEP    [fs] 0.5
    TEMPERATURE [K] 300.0
    &THERMOSTAT
      TYPE   CSVR
      REGION GLOBAL
      &CSVR
        TIMECON [fs] 50.0
      &END CSVR
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 5
      &END EACH
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 0
      &EACH
        MD 20
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD QMMM
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    CHARGE       0                 ! total charge of the QM region only
    MULTIPLICITY 1
    &MGRID
      CUTOFF     400
      REL_CUTOFF 50
      COMMENSURATE                 ! mandatory for E_COUPL GAUSS / S-WAVE (run aborts otherwise)
    &END MGRID
    &QS
      EPS_DEFAULT         1.0E-12
      EXTRAPOLATION       ASPC
      EXTRAPOLATION_ORDER 3
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-6
      MAX_SCF   30
      &OT ON
        MINIMIZER      DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &MM
    &FORCEFIELD
      PARMTYPE OFF                 ! all parameters are given below
      &CHARGE
        ATOM   O
        CHARGE -0.834
      &END CHARGE
      &CHARGE
        ATOM   H
        CHARGE 0.417
      &END CHARGE
      &BOND
        ATOMS O H
        KIND  AMBER
        K     0.24678
        R0    [angstrom] 0.9572
      &END BOND
      &BEND
        ATOMS  H O H
        KIND   AMBER
        K      0.15936
        THETA0 [deg] 104.52
      &END BEND
      &NONBONDED
        &LENNARD-JONES
          ATOMS   O O
          EPSILON [K_e] 76.54        ! = 0.1521 kcal/mol (TIP3P O-O)
          SIGMA   [angstrom] 3.1507
          RCUT    [angstrom] 7.5
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS   O H
          EPSILON [K_e] 0.0
          SIGMA   [angstrom] 1.0
          RCUT    [angstrom] 7.5
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS   H H
          EPSILON [K_e] 0.0
          SIGMA   [angstrom] 1.0
          RCUT    [angstrom] 7.5
        &END LENNARD-JONES
      &END NONBONDED
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA      0.44
        GMAX       32
      &END EWALD
    &END POISSON
  &END MM
  &QMMM
    E_COUPL      GAUSS             ! GEEP electrostatic embedding (GPW/GAPW only)
    USE_GEEP_LIB 15                ! internal Gaussian expansion, no MM_POTENTIAL file needed
    CENTER       EVERY_STEP        ! re-centre the QM box on the QM atoms (default)
    &CELL
      ABC      10.0 10.0 10.0      ! QM box: must contain the QM density, smaller than the MM box
      PERIODIC XYZ
    &END CELL
    &PERIODIC                      ! periodic QM/MM electrostatics (Laino et al. 2006)
      &MULTIPOLE ON                ! decouple the periodic images of the QM box
      &END MULTIPOLE
    &END PERIODIC
    &QM_KIND O
      MM_INDEX 1
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 2 3
    &END QM_KIND
    &MM_KIND O
      RADIUS 0.78                  ! Gaussian width of the smeared MM charge, angstrom
    &END MM_KIND
    &MM_KIND H
      RADIUS 0.44
    &END MM_KIND
    &PRINT                         ! these print keys default to MEDIUM; switch them on for LOW
      &PROGRAM_RUN_INFO ON         ! "Number of QM atoms", "Number of MM atoms"
      &END PROGRAM_RUN_INFO
      &QMMM_CHARGES ON             ! list of MM point charges and their TOTAL CHARGE
      &END QMMM_CHARGES
      &QMMM_LINK_INFO ON           ! QM/MM link summary ("NO QM/MM LINKS DETECTED" here)
      &END QMMM_LINK_INFO
    &END PRINT
  &END QMMM
  &SUBSYS
    &CELL
      ABC      15.54 15.54 15.54
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME   water125.xyz
      COORD_FILE_FORMAT XYZ
      CONN_FILE_FORMAT  GENERATE   ! default: bonds, angles, molecules guessed from distances
    &END TOPOLOGY
    &KIND O                        ! basis/pseudo used for the QM atoms of kind O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `&FORCE_EVAL METHOD` | `QMMM` | 打开 QM/MM 驱动 |
| `PARMTYPE` | `OFF` | 不读参数文件，所有参数写在 `&CHARGE`、`&BOND`、`&BEND`、`&NONBONDED` 中 |
| `&BOND K` / `&BEND K` | 原子单位 | `K` 关键字不做单位换算，示例头部注释给出了 kcal/mol 到 Hartree 的换算 |
| `EPSILON [K_e]` | 76.54 | Lennard-Jones 阱深默认单位是 K；0.1521 kcal/mol ≈ 76.54 K（TIP3P O–O，Jorgensen 1983） |
| `E_COUPL` | `GAUSS` | GPW 的静电嵌入；不写则为机械嵌入 |
| `USE_GEEP_LIB` | 10–15 | 内置高斯展开的项数（2–15），越多越准 |
| `&QMMM &CELL ABC` | QM 原子包络 + 每侧 ≥ 3–4 Å 真空 | QM 盒子决定格点数；太小会截断电子密度，太大浪费 |
| `&MULTIPOLE` | 默认开启 | 去除 QM 盒子周期映像间的相互作用；QM 盒子与 MM 盒子同样大时可关掉 |
| `&MM_KIND RADIUS` | O 0.78、H 0.44 Å | 官方回归测试和教程中的常用值；`CORR_RADIUS` 可另行指定修正半径 |
| `COMMENSURATE` | **必须写** | `E_COUPL GAUSS` 或 `S-WAVE` 时强制要求（GEEP 的势要在多重网格之间插值），漏写会中止：`QM/MM with coupling GAUSS or S-WAVE requires keyword FORCE_EVAL/DFT/MGRID/COMMENSURATE to be enabled.`（2026.2 源码 `qs_environment.F`）。打开后相邻格点的截断比固定为 4 |
| `EXTRAPOLATION ASPC` | 阶数 3 | MD 中用前几步波函数外推初猜，见 [第12章](12-aimd.md) |
| `CHARGE` / `MULTIPLICITY` | QM 区 | 只针对 QM 原子（加上 link 氢） |
| `&QMMM &PRINT` | 打开 `PROGRAM_RUN_INFO`、`QMMM_CHARGES`、`QMMM_LINK_INFO` | 这三个打印键的默认级别是 `MEDIUM`，`PRINT_LEVEL LOW` 时不写就看不到下面要 grep 的内容 |

运行与检查（这些输出字符串可以在 CP2K 源码 `src/qmmm_create.F` 和 `src/qmmm_init.F` 中查到）：

```bash
mpirun -np 8 cp2k.psmp -i water_qmmm_tip3p.inp -o water_qmmm_tip3p.out
grep "Number of QM atoms" water_qmmm_tip3p.out           # 应为 3
grep -A 3 "MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL" water_qmmm_tip3p.out
grep "TOTAL CHARGE:" water_qmmm_tip3p.out                # MM 点电荷总和
grep "NO QM/MM LINKS DETECTED" water_qmmm_tip3p.out      # 本例没有切断共价键
tail -3 water_qmmm_tip3p-1.ener
```

`TOTAL CHARGE:` 应为 0：QM 水分子的 3 个原子从 MM 点电荷中去掉之后，剩下 124 个中性水分子的电荷之和仍为 0。如果你的体系在边界处切断了带电残基，这里会出现非整数电荷，处理方法见 §11。

---

## 9. 使用 AMBER 拓扑：tleap → prmtop → CP2K

实际的生物体系不会把力场写在输入文件里，而是用 AmberTools 的 tleap 生成 `prmtop` + `inpcrd`。示例中的 tleap 脚本：

```bash
cat tleap_water.in
#   source leaprc.water.tip3p
#   sys = loadpdb water125.pdb
#   set sys box { 15.54 15.54 15.54 }
#   check sys
#   saveamberparm sys water125.prmtop water125.inpcrd
tleap -f tleap_water.in            # 生成 water125.prmtop / water125.inpcrd
```

> [!NOTE]
> AmberTools 默认的 TIP3P 残基是为刚性水（SHAKE/SETTLE）准备的：prmtop 中每个水分子有两根 O–H 键和一根 H–H“键”（三者构成刚性三角形），而不是 §8 那样的 O–H 键 + H–O–H 键角。CP2K 不会自动把它们变成约束，而是当作普通谐振键计算，所以这个拓扑描述的水模型与 §8 的柔性水（O–H 键 + 键角）不同。本例只做单点能和力，影响不大；如果要用 prmtop 做 MD，可以在 tleap 中加 `set default FlexibleWater on` 生成柔性水拓扑，或者在 CP2K 中给水分子加 `&CONSTRAINT` 约束。

与 §8 相比，`water_qmmm_amber.inp` 只改了 `&MM`、`&QMMM` 中的 `&MM_KIND` 和 `&SUBSYS &TOPOLOGY`，`RUN_TYPE` 改为 `ENERGY_FORCE`（完整文件见 [示例目录](../examples/16-qmmm-xtb/README.md)）：

```fortran
  &MM
    &FORCEFIELD
      PARMTYPE       AMBER
      PARM_FILE_NAME water125.prmtop
      EI_SCALE14     0.8333333     ! AMBER 1-4 electrostatics scaling (1/1.2)
      VDW_SCALE14    0.5           ! AMBER 1-4 van der Waals scaling (1/2)
      &SPLINE
        RCUT_NB [angstrom] 7.5     ! must stay below half of the 15.54 A box
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA      0.44
        GMAX       32
      &END EWALD
    &END POISSON
  &END MM
  ...
    &MM_KIND OW                    ! kinds are named after the AMBER atom types
      RADIUS 0.78
    &END MM_KIND
    &MM_KIND HW
      RADIUS 0.44
    &END MM_KIND
  ...
    &TOPOLOGY
      CONN_FILE_FORMAT  AMBER
      CONN_FILE_NAME    water125.prmtop
      COORD_FILE_FORMAT CRD
      COORD_FILE_NAME   water125.inpcrd
    &END TOPOLOGY
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `PARMTYPE` | `AMBER` | 从 prmtop 读取力场参数；CHARMM 用 `CHM` 加 `.prm` 文件 |
| `PARM_FILE_NAME` | prmtop 文件 | 与 `CONN_FILE_NAME` 通常是同一个文件 |
| `EI_SCALE14` / `VDW_SCALE14` | 1/1.2、0.5 | AMBER 力场的 1-4 缩放；CP2K 的默认值与 AMBER 不同，必须显式写出 |
| `&SPLINE RCUT_NB` | ≤ MM 盒子最短边的一半 | 非键截断；小盒子时必须调小 |
| `CONN_FILE_FORMAT` | `AMBER` | 从 prmtop 读取成键关系与分子划分；CHARMM 用 `PSF` |
| `COORD_FILE_FORMAT` | `CRD` | 读 `inpcrd`/`rst7`；也可以用 PDB |
| `&MM_KIND` 名称 | AMBER 原子类型 | 读 prmtop 时 MM 原子的 kind 名取自原子类型（TIP3P 水为 `OW`、`HW`）；同一类型的原子若电荷不同，CP2K 会在类型名后加数字后缀区分（如 `CT1`、`CT2`，见 `src/topology_amber.F`），写 `&MM_KIND` 前先在输出中核对 kind 名；QM 原子的 kind 名则取自 `&QM_KIND`（这里是 `O`、`H`），所以 `&SUBSYS &KIND` 仍写 `O`、`H` |

> [!WARNING]
> AMBER 力场中 TIP3P 水的 H 以及丝氨酸、苏氨酸、酪氨酸羟基的 H（HO 类型）没有 Lennard-Jones 参数。纯 MM 时没有问题，但在静电嵌入下，这些裸露的正电荷会被 QM 电子密度强烈吸引，严重时导致 MD 崩溃。官方 QM/MM 教程的做法是用 ParmEd 给这些类型补上 LJ 参数，例如 `changeLJSingleType :WAT@H1 0.3019 0.047`（Rmin/2 = 0.3019 Å，ε = 0.047 kcal/mol，即 GAFF2 醇羟基 H 的参数），然后 `outparm` 写出新的 prmtop。

在正式的生物体系研究中，QM/MM 之前必须先完成纯 MM 的准备：能量最小化 → NVT → NPT 平衡（可以用 CP2K 的 `METHOD FIST`，也可以用 AMBER/GROMACS），再从平衡结构开始 QM/MM。官方教程（分支酸变位酶 + AM1/MM 元动力学）完整演示了这一流程，链接见参考资料。

---

## 10. 半经验 QM/MM：GFN1-xTB 嵌入水溶液

把 §8 的 QM 部分换成原生 GFN1-xTB，每个 MD 步的成本可以下降到秒级以下，适合测试 QM 区选择、采样长度和分析脚本，或者在 DFT/MM 之前做预平衡。`xtb_qmmm_water.inp` 与 §8 的区别集中在 `&DFT` 和 `&QMMM`（完整文件见示例目录）：

```fortran
  &DFT
    CHARGE       0                 ! charge of the QM region only
    MULTIPLICITY 1
    &QS
      METHOD      XTB
      EPS_DEFAULT 1.0E-10
      &XTB
        CHECK_ATOMIC_CHARGES T
      &END XTB
    &END QS
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF   1.0E-7
      MAX_SCF   100
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA  0.3
      &END MIXING
    &END SCF
  &END DFT
  ...
  &QMMM
    E_COUPL COULOMB                ! analytical 1/r coupling for tight-binding QM regions
    CENTER  EVERY_STEP
    &CELL
      ABC      10.0 10.0 10.0      ! QM box, must stay smaller than the MM box
      PERIODIC XYZ
    &END CELL
    &QM_KIND O
      MM_INDEX 1
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 2 3
    &END QM_KIND
    &PRINT                         ! default level MEDIUM; switched on for PRINT_LEVEL LOW
      &PROGRAM_RUN_INFO ON
      &END PROGRAM_RUN_INFO
      &QMMM_LINK_INFO ON
      &END QMMM_LINK_INFO
    &END PRINT
  &END QMMM
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `E_COUPL` | `COULOMB` | 紧束缚 QM 区的静电嵌入：MM 电荷与 xTB/DFTB 的 Mulliken 电荷通过解析 $1/r$ 作用 |
| `&PERIODIC` | 本例不写 | 与官方回归测试 `tests/QMMM/xTB/regtest/wdim_scc.inp` 一致：MM–MM 是周期性的（SPME），QM–MM 静电不做周期处理 |
| QM 原子的 `&KIND` | 不需要 | xTB 不使用 GTO 基组文件 |
| `&MIXING` | `DIRECT_P_MIXING` | 原生 GFN1 用 CP2K 的电荷混合（2026.2 中 `SCC_MIXER AUTO` 对原生 GFN1 即保持此行为） |

换成 DFTB/MM 只需把 `&QS` 改为 §6 的 `METHOD DFTB` 写法（QM 区为团簇时去掉 `DO_EWALD`），`E_COUPL` 同样用 `COULOMB`。

> [!NOTE]
> 关于紧束缚 QM/MM 的周期性和耦合方式：CP2K 2026.2 的源码中，xTB 和 DFTB 的周期性 QM/MM 走 Ewald 路径（`&QMMM &PERIODIC &POISSON &EWALD`），回归测试 `tests/QMMM/xTB/regtest/fdeb_7.inp` 使用 `E_COUPL POINT_CHARGE` 加 `&PERIODIC`；经典半经验方法（AM1、PM6 等）则会报 `QM/MM periodic calculations not implemented for semi empirical methods`。trunk 的输入参考中 `GAUSS` 已扩展到 DFTB/xTB（“Gaussian MM charges for DFTB/xTB with NOCOMPATIBILITY”），预计进入下一个正式版本。在自己的版本上使用前，先用 `cp2k.psmp -c` 检查语法，再用短 MD 检查能量守恒。

```bash
cp2k.psmp -i xtb_qmmm_water.inp -o xtb_qmmm_water.out
grep "Number of QM atoms" xtb_qmmm_water.out
tail -3 water_qmmm_xtb-1.ener
```

---

## 11. QM 区选择与 link atom

### 11.1 QM 区怎么选

- **必须包含**：所有成键或断键的原子、发生质子转移的供体和受体、金属离子及其第一配位层，以及直接参与电荷转移的 π 体系；
- **切键位置**：选择极性最小、远离反应中心的 C–C 单键（蛋白残基通常切 Cα–Cβ），绝不切极性键、双键或芳香环；切点离反应中心至少 3 个键；
- **规模**：GPW 下 50–200 个 QM 原子比较常见；先用 xTB/MM 比较不同 QM 区大小下的能垒或几何，再决定 DFT/MM 的 QM 区；
- **电荷**：统计 QM 区（含 link 氢）的形式电荷，写进 `&DFT CHARGE`；奇数电子数说明漏了或多了 H，或者需要 `UKS`；
- **收敛性检查**：能垒随 QM 区增大的变化小于 1–2 kcal/mol，才说明 QM 区选够了。

### 11.2 `&LINK` 写法（模板）

下面是一个**模板**，原子编号仅作示意：假设蛋白某残基在 Cα（MM，原子 1411）与 Cβ（QM，原子 1413）之间切断。

```fortran
  &QMMM
    E_COUPL      GAUSS
    USE_GEEP_LIB 15
    &CELL
      ABC 20.0 20.0 20.0
    &END CELL
    &QM_KIND C
      MM_INDEX 1413 1416 1419
    &END QM_KIND
    # ... further &QM_KIND sections for the other QM atoms ...
    &LINK
      QM_INDEX          1413       ! QM atom of the cut bond (C-beta)
      MM_INDEX          1411       ! MM atom of the cut bond (C-alpha)
      QM_KIND           H          ! capping atom element (default H)
      LINK_TYPE         IMOMM      ! default
      ALPHA_IMOMM       1.38       ! r_eq(C-C) / r_eq(C-H), default 1.38
      QMMM_SCALE_FACTOR 0.0        ! remove the C-alpha charge from the QM/MM potential only
    &END LINK
  &END QMMM
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `QM_INDEX` / `MM_INDEX` | 被切断键两端的原子 | 每个被切断的键写一个 `&LINK` 段 |
| `QM_KIND` | `H`（默认） | 封端原子的元素；该元素在 `&SUBSYS &KIND` 中要有基组和赝势 |
| `LINK_TYPE` | `IMOMM`（默认） | 另有 `GHO`、`PSEUDO` |
| `ALPHA_IMOMM` | 1.38（默认） | 按力场键长估算：$\alpha=r_{\mathrm{eq}}(\mathrm{QM{-}MM})/r_{\mathrm{eq}}(\mathrm{QM{-}H})$，C–C/C–H 约 1.53/1.09 ≈ 1.40 |
| `QMMM_SCALE_FACTOR` | 0.0–1.0（默认 1.0） | 只缩放 MM 边界原子电荷在 QM/MM 静电势中的贡献，不影响纯 MM 静电；常设 0 以免它过度极化附近的封端氢 |
| `FIST_SCALE_FACTOR` | 默认 1.0 | 同时修改 FIST（MM）中的电荷，并由此进入 QM/MM 势 |
| `RADIUS` / `CORR_RADIUS` | 默认同 `&MM_KIND` | 只为该 MM 边界原子单独设置 GEEP 半径 |

### 11.3 边界电荷修正

把残基的侧链移入 QM 区后，这些原子的 MM 电荷从 MM 体系中移除，残留在主链上的 MM 电荷之和通常不是整数。常见的修正方式：

1. **重新分配到 MM 主链原子上**（官方教程的做法）：计算残留电荷（例如官方教程中 Arg 侧链移出后为 −0.0362 e），把等量反号的电荷平均分配到该残基剩下的主链原子上，用 ParmEd 的 `change charge @<原子号> <新电荷>` 写进新的 prmtop，再让 `PARM_FILE_NAME` 和 `CONN_FILE_NAME` 指向它；
2. **在输入中处理**：`&LINK` 内的 `&ADD_MM_CHARGE`（在 `ATOM_INDEX_1` 与 `ATOM_INDEX_2` 的连线上，按 `ALPHA` 确定的位置添加一个电荷 `CHARGE`）或 `&MOVE_MM_CHARGE`（把边界 MM 原子的电荷沿连线移动）；
3. **只关闭边界原子的 QM/MM 电荷**：`QMMM_SCALE_FACTOR 0.0`，是最简单但最粗糙的做法。

无论用哪种方式，都要检查两点：QM 区电荷 + MM 电荷总和等于体系的真实总电荷（输出中的 `TOTAL CHARGE:` 是 MM 点电荷之和）；不同方案下关键能量差的变化在可接受范围内。

运行后检查 link 是否被识别（`PRINT_LEVEL LOW` 时先在 `&QMMM &PRINT` 中打开 `&QMMM_LINK_INFO` 和 `&QMMM_CHARGES`，见 §8）：

```bash
grep -A 6 "QM/MM LINKS" my_qmmm.out          # 有切键时应看到 "IMOMM TYPE LINK" 和每个 link 的信息
grep "TOTAL CHARGE:" my_qmmm.out
```

> [!CAUTION]
> 一个 MM 原子只能连接一个 QM 原子。两个 `&LINK` 段共用同一个 `MM_INDEX` 时，程序报 `In the IMOMM scheme no more than one QM atom can be bounded to the same MM atom. Multiple link MM atom not allowed. Check your link sections.`，此时应把切点移到别处。

---

## 12. 周期性 QM/MM 与 QM 盒子

`&QMMM &PERIODIC` 让 QM 密度感受到所有 MM 周期映像的长程静电。对溶液、晶体和带电体系，这比只用最近映像更物理，代价是额外的 Ewald 计算。

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `&PERIODIC GMAX` | 1.0 | QM/MM Ewald 倒空间求和的最大 G |
| `&PERIODIC NGRIDS` | 50 50 50 | G 空间项插值格点 |
| `&PERIODIC REPLICA` | −1（自动） | 实空间求和考虑的映像数 |
| `&MULTIPOLE` | 默认开启（即使不写） | 用密度拟合得到的原子电荷去除 QM 盒子周期映像间的耦合；QM 盒子与 MM 盒子一样大时可写 `&MULTIPOLE OFF` |
| `&MULTIPOLE EWALD_PRECISION` | 1E-6 Hartree | 多极去耦 Ewald 的精度 |
| `&PERIODIC &POISSON &EWALD` | — | 紧束缚 QM 区的周期性 QM/MM 使用的 Ewald 参数（见 §10 的回归测试） |

QM 盒子相关的建议：

- QM 盒子在每个方向上至少比 QM 原子的包络大 6–8 Å（每侧 3–4 Å），保证电子密度在盒子边界处衰减到可以忽略；
- QM 盒子不能大于 MM 盒子，通常明显小于 MM 盒子。若两者一样大（例如整个周期体系都是 QM，只有少量 MM），请关闭 `&MULTIPOLE`；
- MD 中 QM 原子可能漂移到 QM 盒子边缘，程序会警告 `One or few QM atoms are within the SKIN of the quantum box. Check your run and you may possibly consider: the activation of the QMMM WALLS around the QM box, switching ON the centering of the QM box or increase ...`。处理办法：保持 `CENTER EVERY_STEP`，增大 QM 盒子，或加 `&WALLS`（`TYPE REFLECTIVE` 默认 / `QUADRATIC`，`WALL_SKIN` 默认 0.5 Å）；
- 机械嵌入（`E_COUPL NONE`）时写 `&PERIODIC` 没有意义，程序会提示 `QMMM periodic calculation with coupling NONE was requested! Switching off the periodic keyword ...` 并自动关闭；
- 使用 `&IMAGE_CHARGE`（金属表面镜像电荷）时，QM 盒子与 MM 盒子必须一样大，否则报 `QM and MM box need to have the same size when using image charges`。

---

## 13. QM/MM 下的 NEB 与 Dimer 过渡态搜索

酶或溶液中的反应能垒通常用 QM/MM 势能面上的 CI-NEB 加 Dimer 求得。完整的实战教程是 [arvpinto/enzyme_neb_cp2k](https://github.com/arvpinto/enzyme_neb_cp2k)（磷脂酶 A2 水解反应，PBE/AMBER 静电嵌入，CC0 协议），本节总结其流程，方法本身见 [第13章](13-neb-dimer.md) 和 [第14章](14-vibration-thermo.md)。

**流程**

1. **MM 平衡**：完成溶剂化、最小化、NVT/NPT 平衡，取一帧代表性结构；
2. **QM/MM 优化反应物**：把 `&FORCE_EVAL`（QM/MM 设置）写成一个 include 文件（enzyme_neb_cp2k 中叫 `forceeval_qmmm.cp2k.inc`），供后续所有输入共用。离 QM 区较远的原子（例如 6–8 Å 以外）全部固定，只优化活性中心附近的区域；
3. **构造初始路径**：沿集体变量（例如两个键长之差）做受约束的势能面扫描，把扫描得到的结构作为 NEB 的初始副本；
4. **CI-NEB**：所有副本使用同一组固定原子，`ROTATE_FRAMES F`、`ALIGN_FRAMES F`（大部分原子被固定时，副本之间不能再整体转动或平移）；
5. **部分 Hessian 频率**：用 `&MODE_SELECTIVE` 只求与反应原子相关的模式，确认最高副本有一个沿反应方向的虚频；
6. **Dimer 精修**：用该虚频模式初始化 `&DIMER_VECTOR`，收敛到鞍点后再做一次部分 Hessian 验证。

**固定原子列表**：大体系的 `&FIXED_ATOMS` 列表可能有上万个原子，不能写在一行里。enzyme_neb_cp2k 的做法是每 100 个原子写一行 `LIST`，生成一个 include 文件：

```bash
# fixed.txt: one atom index per line (e.g. exported from VMD: atoms beyond 8 A of the QM region)
awk 'NR % 100 == 1 {if (NR > 1) print ""; printf "LIST"} {printf " %s", $0} END {print ""}' fixed.txt > fixed_atoms.inc
```

**CI-NEB 输入模板**（`forceeval_qmmm.inc` 中是完整的 `&FORCE_EVAL METHOD QMMM ... &END FORCE_EVAL`，副本坐标来自上一步扫描）：

```fortran
&GLOBAL
  PROJECT  enzyme_cineb
  RUN_TYPE BAND
  PRINT_LEVEL LOW
&END GLOBAL

@INCLUDE 'forceeval_qmmm.inc'

&MOTION
  &BAND
    BAND_TYPE         CI-NEB
    NUMBER_OF_REPLICA 12
    NPROC_REP         ${NPROC_REP-32}
    K_SPRING          0.05
    ROTATE_FRAMES     F
    ALIGN_FRAMES      F
    &CI_NEB
      NSTEPS_IT 5
    &END CI_NEB
    &OPTIMIZE_BAND
      OPT_TYPE            DIIS
      OPTIMIZE_END_POINTS F        ! end points are already optimised minima
    &END OPTIMIZE_BAND
    &REPLICA
      COORD_FILE_NAME scan_01.xyz
    &END REPLICA
    # ... one &REPLICA per scan structure; missing ones are interpolated by CP2K ...
    &REPLICA
      COORD_FILE_NAME scan_12.xyz
    &END REPLICA
  &END BAND
  &CONSTRAINT
    &FIXED_ATOMS
      @INCLUDE 'fixed_atoms.inc'
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `NPROC_REP` × `NUMBER_OF_REPLICA` | = 总 MPI 进程数 | enzyme_neb_cp2k 中 16 个副本 × 64 进程 = 1024 核；各副本的报错在 `-BAND*.out` 中 |
| `ROTATE_FRAMES` / `ALIGN_FRAMES` | `F` | 大部分原子固定时必须关闭 |
| `OPTIMIZE_END_POINTS` | `F` | 端点已在同一 QM/MM 设置下严格优化 |
| `K_SPRING` | 0.05 | 副本分布不均时适当增大 |
| `&FIXED_ATOMS` | 所有输入一致 | 优化、NEB、频率、Dimer 必须使用同一组固定原子，否则能量不可比 |

Dimer 与部分 Hessian 的写法分别见 [第13章 §8](13-neb-dimer.md) 和 [第14章](14-vibration-thermo.md)。QM/MM 下需要注意：`&MODE_SELECTIVE` 的 `ATOMS` 写反应原子，`&INVOLVED_ATOMS`（例如 `RANGE -5000 0`）用来追踪虚频模式；`&DIMER_VECTOR` 需要给出**全部原子**的分量（每行一个原子，归一化），enzyme_neb_cp2k 用自带脚本把部分 Hessian 得到的虚频模式写成这种格式。dimer 向量最好在 Dimer 起始的那个结构上计算，否则初始方向可能是正曲率，搜索会多走很多步。NEB 能量曲线可以用 `python ../scripts/neb_plot.py -h` 查看作图脚本的用法，频率与零点能用 `python ../scripts/vib_thermo.py -h`。

> [!TIP]
> GGA 泛函通常低估反应能垒。enzyme_neb_cp2k 的作者在 PBE/MM 路径上找到鞍点后，再用更高精度的方法对关键结构做单点修正。比较省时的组合是：xTB/MM 上做扫描和 NEB 预收敛 → DFT/MM 上 CI-NEB + Dimer → 杂化泛函（ADMM，见 [第9章](09-hybrid-admm.md)）单点。

---

## 14. GROMACS–CP2K 接口

当 MM 部分很大（蛋白 + 膜 + 溶剂，几十万原子），或者需要 GROMACS 的约束算法、自由能工具、增强采样生态时，可以让 GROMACS 负责 MM 与积分，CP2K 只作为 QM 引擎计算 QM 区的能量和力。GROMACS 自 2022 版起提供这一接口。

**工作方式**（GROMACS 手册）：QM–MM 静电使用 GEEP 静电嵌入；QM 原子的 MM 电荷被置零，QM–QM 的 LJ 被排除，QM–MM 的 LJ 仍由 GROMACS 计算；被切断的键由 CP2K 自动加封端氢；所有拓扑修改在 `gmx grompp` 预处理阶段自动完成。

**编译**：GROMACS 需要静态链接 `libcp2k`：CMake 选项 `-DGMX_CP2K=ON`、`-DCP2K_DIR=<libcp2k 所在目录>`，以及 `-DBUILD_SHARED_LIBS=OFF -DGMXAPI=OFF -DGMX_INSTALL_NBLIB_API=OFF`，FFT、BLAS、LAPACK 要与 CP2K 使用同一套。GROMACS 安装指南中生成 `libcp2k` 的步骤仍是旧的 Makefile 命令 `make ... libcp2k`，而 CP2K 从 2026.1 起只支持 CMake 构建，请按 CP2K 的构建文档用 CMake 生成库，必要时参考 BioExcel 教程中的说明。

**mdp 设置**（示例 `gromacs_qmmm_fragment.mdp`，追加到一个已平衡的普通 MD 的 mdp 中）：

```ini
qmmm-cp2k-active         = true
qmmm-cp2k-qmgroup        = QMatoms   ; index group with the QM atoms
qmmm-cp2k-qmmethod       = PBE       ; PBE, PBE-D3, BLYP, BLYP-D3, B3LYP, ... or INPUT
qmmm-cp2k-qmcharge       = 0         ; total charge of the QM region
qmmm-cp2k-qmmultiplicity = 1         ; spin multiplicity of the QM region
; qmmm-cp2k-qmfilenames  =           ; default: <tpr name>_cp2k
```

| mdp 选项 | 说明 |
|---|---|
| `qmmm-cp2k-active` | `true` 打开 QM/MM，需要链接了 CP2K 的 GROMACS |
| `qmmm-cp2k-qmgroup` | index 文件中 QM 原子组的名字（默认 `System`，即全部原子） |
| `qmmm-cp2k-qmmethod` | 内置方法（PBE、BLYP 及其 `-D3` 版本；较新版本另有 B3LYP、PBE0、CAM-B3LYP、WB97X 等，配 6-31G*）或 `INPUT`；请以所用 GROMACS 版本的 mdp 文档为准 |
| `qmmm-cp2k-qmcharge` / `qmmm-cp2k-qmmultiplicity` | QM 区电荷与自旋多重度 |
| `qmmm-cp2k-qmfilenames` | CP2K 的 `.inp`、`.out`、`.pdb` 文件名前缀，默认是 tpr 名加 `_cp2k` |

```bash
gmx make_ndx -f conf.gro -o index.ndx                               # 建立 QMatoms 组
gmx grompp -f md_qmmm.mdp -c conf.gro -p topol.top -n index.ndx -o qmmm.tpr
gmx mdrun -deffnm qmmm                                               # 生成 qmmm_cp2k.inp / .out / .pdb
echo "Quantum" | gmx energy -f qmmm.edr -o quantum.xvg              # 能量文件中新增的 "Quantum En." 项
```

**自定义 CP2K 输入**：`qmmm-cp2k-qmmethod = INPUT`，并在 `gmx grompp` 中用 `-qmi my_cp2k.inp` 传入。该文件有以下限制：`RUN_TYPE` 必须为 `ENERGY_FORCE`；必须有 `CHARGE`、`MULTIPLICITY`；`COORD_FILE_NAME` 指向 pdb 文件，并同时写 `COORD_FILE_FORMAT PDB` 和 `CHARGE_EXTENDED TRUE`；不允许 `@INCLUDE`。最稳妥的做法是先用内置方法跑几步，把自动生成的 `qmmm_cp2k.inp` 拿来修改。

**限制**（GROMACS 手册）：模拟过程中 QM 区拓扑不能改变；忽略拓扑中的 B 态参数，不能做双拓扑自由能微扰；QM/MM 力对维里没有贡献，NPT 下压力耦合可能出现伪影，因此生产模拟建议用 NVT。

学习路线建议：GROMACS 手册的 QM/MM 页面 → BioExcel 的 [gromacs-2022-cp2k-tutorial](https://github.com/bioexcel/gromacs-2022-cp2k-tutorial)（附 Handbook，包含 NMA、stilbene、eGFP 等示例）→ BioExcel QM/MM 最佳实践指南 → 华东师大 [bwang-ecnu/gromacs_cp2k_tutorial](https://github.com/bwang-ecnu/gromacs_cp2k_tutorial) 的示例体系。

---

## 参数选择与收敛建议

| 项目 | 建议 | 理由 |
|---|---|---|
| xTB 方法 | 分子结构、构象用 GFN2；GFN2 不可用时用 GFN1；超大体系预处理用 GFN0 | GFN2 的静电与色散更完整 |
| `ACCURACY`（tblite） | 单点、MD 用 1.0；优化、频率、应力用 0.1 | 越小越严，也越慢 |
| `EPS_DEFAULT`（xTB/DFTB） | 1E-10；频率计算 1E-12 | 力的数值噪声 |
| 周期性 xTB/DFTB | ≤ 2026.1 写 `DO_EWALD T`；2026.2 只写 `&CELL PERIODIC` | 漏写时旧版本会按非周期处理静电 |
| xTB 的 k 点 | 按 DFT 的做法做收敛测试；对角化 + 混合 | 2026.2 中 OT 不支持 k 点 |
| QM 区 | 包含所有成键/断键原子，切点离反应中心 ≥ 3 个键 | 边界误差随距离迅速衰减 |
| `E_COUPL` | GPW：`GAUSS` + `USE_GEEP_LIB 10–15`；xTB/DFTB：`COULOMB` | 默认 `NONE` 只是机械嵌入 |
| `&MM_KIND RADIUS` | O 0.78、H 0.44 Å，其他元素取共价半径量级 | 防止 QM 电子溢出到 MM 正电荷上 |
| QM 盒子 | QM 原子包络 + 每侧 3–4 Å；保持 `CENTER EVERY_STEP` | 截断密度会造成能量与力的误差 |
| `CUTOFF`（QM/MM） | 与纯 QM 一样做收敛测试，见 [第4章](04-cutoff-convergence.md)；`GAUSS`/`S-WAVE` 耦合必须写 `COMMENSURATE` | 不写会直接中止（见常见问题 16） |
| AMBER 拓扑 | `EI_SCALE14 0.8333333`、`VDW_SCALE14 0.5`；为没有 LJ 的 H 补参数 | CP2K 默认值与 AMBER 不同；裸露的正电荷会导致 QM 区坍缩 |
| 时间步长 | 含 H 时 0.5 fs | 与 AIMD 相同 |
| NEB/Dimer | 所有步骤使用相同的固定原子；`ROTATE_FRAMES F`、`ALIGN_FRAMES F` | 能量可比性 |
| 准确度验证 | xTB/MM 只用于采样与筛选；能垒用 DFT/MM 加更高级单点修正 | 半经验方法的能垒误差可达数个 kcal/mol 以上 |

---

## 常见问题与报错

**1. `[ABORT] SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`（xTB/DFTB 或 QM/MM 都可能出现）**
原因：初猜不好、混合太激进、QM 电荷或多重度写错、QM 区被附近裸露的 MM 正电荷过度极化。解决：xTB 用 `SCF_GUESS MOPAC`；对角化时减小 `ALPHA`，或在 2026.2 中尝试 `SCC_MIXER TBLITE`/`CP2K`；检查 `CHARGE`/`MULTIPLICITY`；为没有 LJ 参数的 H 补参数（§9）；其他 SCF 策略见 [第5章](05-scf-strategy.md)。

**2. `The XTB/TBLITE section requires XTB/GFN_TYPE TBLITE` 或 `XTB/GFN_TYPE TBLITE requires an XTB/TBLITE section`（2026.2）**
原因：2026.2 中 `GFN_TYPE TBLITE` 与 `&TBLITE` 段必须同时出现。把 2025.2/2026.1 的输入（只有 `&TBLITE`）直接拿到 2026.2 运行，会触发第一条；只写了 `GFN_TYPE TBLITE` 而没有 `&TBLITE` 段，会触发第二条。解决：两者都写上。反过来，2025.1–2026.1 中 `GFN_TYPE` 是整数关键字，写 `GFN_TYPE TBLITE` 会在解析阶段报错，应去掉这一行（示例中用 `-E NEWTB=0`）。

**3. `Built without TBLITE`**
原因：可执行文件编译时没有启用 tblite。解决：`cp2k.psmp -v | grep -o tblite` 检查；换用带 tblite 的构建，或改用原生 GFN1（去掉 `GFN_TYPE TBLITE` 和 `&TBLITE`）。

**4. 找不到 xTB/DFTB 参数文件（`xTB1_parameters`、`scc_parameter`、`.spl`）。**
原因：`CP2K_DATA_DIR` 未设置或指向了错误的目录；`PARAM_FILE_PATH` 写成了绝对路径以外的错误相对路径。解决：`export CP2K_DATA_DIR=<cp2k>/data`，并用 §3 的 `ls` 命令确认文件存在。

**5. DFTB 打开 `DISPERSION T` 后程序直接段错误（Segmentation fault）。**
原因：没有给出 `UFF_FORCE_FIELD`，手册明确说明这种情况会崩溃。解决：在 `&PARAMETER` 中写 `UFF_FORCE_FIELD ../uff_table`（相对 `PARAM_FILE_PATH`）。若改用 `DISPERSION_TYPE D3` 或 `D3(BJ)`，还要用 `DISPERSION_PARAMETER_FILE dftd3.dat` 指定 D3 参数文件。

**6. 周期性 xTB/DFTB 的能量与 Γ 点期望不符、随盒子大小变化，或在 ≤ 2026.1 中能量明显不对。**
原因：2026.1 及以前没写 `DO_EWALD T`，周期体系的长程静电被当作非周期处理。解决：旧版本写 `DO_EWALD T` 并配置 `&POISSON &EWALD`；2026.2 中确认 `&CELL PERIODIC` 正确即可。

**7. `QM/MM periodic calculations not implemented for semi empirical methods`**
原因：AM1、PM6 等经典半经验方法（`METHOD AM1/PM6/...`）不支持 `&QMMM &PERIODIC`。解决：去掉 `&PERIODIC`，或改用 xTB/DFTB 或 GPW。

**8. `One or few QM atoms are within the SKIN of the quantum box. ...`**
原因：MD 中 QM 原子移动到 QM 盒子边缘。解决：`CENTER EVERY_STEP`（默认），增大 QM 盒子，或加 `&WALLS`（§12）。

**9. `In the IMOMM scheme no more than one QM atom can be bounded to the same MM atom. Multiple link MM atom not allowed. Check your link sections.`**
原因：两个 `&LINK` 共用一个 MM 原子。解决：换一个切点，让每个 MM 边界原子只连接一个 QM 原子。

**10. `Use the LSD option for an odd number of electrons`（QM/MM 中）**
原因：QM 区电荷写错，或选择 QM 区时漏选、多选了 H。注意 link 氢也计入 QM 电子数。解决：重新统计 QM 区（含封端氢）的形式电荷；确实是开壳层时加 `UKS` 和正确的 `MULTIPLICITY`。

**11. QM/MM MD 几十步后 QM 区结构坍缩、出现 `GEOMETRY wrong or EMAX_SPLINE too small`，或者能量突然大幅下降。**
原因：MM 中没有 LJ 参数的 H（TIP3P、羟基 H）与 QM 区电子发生“静电坍缩”；或者初始结构中原子重叠。解决：用 ParmEd 补 LJ 参数（§9）；先做纯 MM 最小化和平衡；检查 `&MM_KIND RADIUS`。

**12. 忘了写 `E_COUPL`，结果与气相计算几乎一样。**
原因：默认 `E_COUPL NONE`，即机械嵌入。解决：GPW 写 `E_COUPL GAUSS` + `USE_GEEP_LIB`，xTB/DFTB 写 `COULOMB`；输出中应出现 `MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL` 一节。

**13. `QMMM periodic calculation with coupling NONE was requested! Switching off the periodic keyword ...`**
原因：在机械嵌入下写了 `&PERIODIC`。这只是警告，程序会自动关闭周期性；若本意是静电嵌入，请补上 `E_COUPL`。

**14. NEB 某个副本在 QM/MM 下 SCF 不收敛，或者副本结构整体“漂移”。**
原因：线性插值使 QM 原子靠得过近；固定原子在各输入中不一致；`ROTATE_FRAMES`/`ALIGN_FRAMES` 没有关闭。解决：用扫描结构代替线性插值；所有输入共用 `fixed_atoms.inc`；两个关键字都设为 `F`。报错在对应的 `-BAND*.out` 中（见 [第13章](13-neb-dimer.md)）。

**15. GROMACS–CP2K 的 NPT 模拟中盒子尺寸异常变化。**
原因：GROMACS 手册说明 QM/MM 力不计入维里。解决：生产模拟用 NVT；盒子先用纯 MM 的 NPT 平衡好。

**16. `QM/MM with coupling GAUSS or S-WAVE requires keyword FORCE_EVAL/DFT/MGRID/COMMENSURATE to be enabled.`**
原因：GPW 的 GEEP 耦合要求粗细格点互为整数倍，`&MGRID` 中漏写了 `COMMENSURATE`。解决：在 `&DFT &MGRID` 中加一行 `COMMENSURATE`（此时不能再用 `MULTIGRID_SET` + `MULTIGRID_CUTOFF` 手动指定各级格点的截断，`PROGRESSION_FACTOR` 也被固定为 4）；之后按 [第4章](04-cutoff-convergence.md) 重新确认 `CUTOFF` 收敛。

---

## 练习

1. **GFN0 / GFN1 / GFN2 对比**：用 `ethanol_gfn1_native.inp` 分别以 `GFN_TYPE 0`、默认 GFN1、tblite GFN2 优化乙醇，比较 C–O 键长和 C–C–O 键角，并与 [第7章](07-geo-cell-opt.md) 的 PBE 优化结果对照。
2. **k 点收敛**：把冰示例的 k 点网格依次改为 1×1×1、2×2×1、3×3×2、4×4×2、6×6×3，作出每个水分子能量对 k 点数的曲线，判断多大的网格足够。
3. **自旋态**：对 O2 分别计算 `MULTIPLICITY 1`（去掉 `UKS` 并改为闭壳层）和 `MULTIPLICITY 3`，比较能量与键长，并说明为什么定量结论仍需要 DFT 验证。
4. **DFTB 液态水结构**：用 `dftb_md.inp` 跑 10 ps（丢弃前 1–2 ps），用 `python ../scripts/md_analyze.py -h` 查看用法后计算 O–O RDF，与实验或 [第12章](12-aimd.md) 的 AIMD 结果比较。
5. **嵌入方式的影响**：在 `water_qmmm_amber.inp` 中把 `E_COUPL` 依次设为 `NONE`、`GAUSS`，比较 QM 水分子的偶极矩（在 `&DFT &PRINT &MOMENTS` 中打开）和 QM 原子受力。
6. **QM 盒子与 GEEP**：把 QM 盒子从 8 Å 增加到 14 Å，`USE_GEEP_LIB` 从 6 增加到 15，观察 QM/MM 能量的变化，确定收敛的设置。
7. **xTB/MM 与 DFT/MM**：分别运行 `xtb_qmmm_water.inp` 和 `water_qmmm_tip3p.inp` 各 1 ps，比较每步耗时和 QM 水分子 O–H 键长的分布。
8. **设计一个 link**：从 PDB 下载一个小肽，用 tleap 生成 prmtop，把某个 Lys 侧链设为 QM 区（在 Cα–Cβ 之间切断），写出 `&QM_KIND`、`&LINK`，用 ParmEd 修正边界电荷，并用 `grep "TOTAL CHARGE:"` 检查电荷是否正确。

---

## 参考资料

**官方文档**

- Extended Tight Binding（原生 GFN0/GFN1 与 tblite GFN2，含 k 点与 spGFN2 示例）：<https://manual.cp2k.org/trunk/methods/semiempiricals/xtb.html>
- QM/MM with Built-in Force Field（分支酸变位酶：AMBER 拓扑、ParmEd 补 LJ、link atom、电荷修正、元动力学）：<https://manual.cp2k.org/trunk/methods/qm_mm/builtin.html>
- QM/MM with GROMACS：<https://manual.cp2k.org/trunk/methods/qm_mm/gromacs.html>
- QM/MM 方法总览：<https://manual.cp2k.org/trunk/methods/qm_mm/index.html>
- 输入参考（2026.2）`&XTB`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/XTB.html>；`&TBLITE`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/XTB/TBLITE.html>
- 输入参考（2026.2）`&DFTB`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/DFTB.html>；`&DFTB/&PARAMETER`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS/DFTB/PARAMETER.html>
- 输入参考（2026.2）`&QMMM`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/QMMM.html>（子段 `LINK`、`PERIODIC`、`WALLS`、`MM_KIND`、`QM_KIND` 在同目录下）
- 旧版 HowTo：QM/MM for biochemical systems：<https://www.cp2k.org/howto:biochem_qmmm>
- GROMACS 手册：QM/MM with CP2K interface：<https://manual.gromacs.org/current/reference-manual/special/qmmm.html>；mdp 选项：<https://manual.gromacs.org/current/user-guide/mdp-options.html>；安装指南：<https://manual.gromacs.org/current/install-guide/index.html>
- BioExcel：Best Practices in QM/MM Simulation of Biomolecular Systems：<https://docs.bioexcel.eu/qmmm_bpg/en/main/index.html>
- tblite 文档（参数文件格式）：<https://tblite.readthedocs.io/en/latest/spec/parameter.html>
- CP2K 2026.2 release notes：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QMMM/QS/`、`tests/QMMM/xTB/`、`tests/QMMM/DFTB/`、`tests/xTB/`、`tests/DFTB/` 是查关键字真实写法最可靠的来源
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：`qm_mm/` 目录是官方 QM/MM 教程的输入文件
- [bioexcel/gromacs-2022-cp2k-tutorial](https://github.com/bioexcel/gromacs-2022-cp2k-tutorial)：GROMACS–CP2K 接口教程与 Handbook
- [bioexcel/cp2K_qmmm_tutorials_for_biological_simulations](https://github.com/bioexcel/cp2K_qmmm_tutorials_for_biological_simulations)：CP2K 原生 `&QMMM` 教程（AMBER prmtop、link atom、电荷修正），版本较老
- [bioexcel/2020-10-cp2k_qmmm_intro](https://github.com/bioexcel/2020-10-cp2k_qmmm_intro)：ARCHER2 QM/MM 实践课材料
- [arvpinto/enzyme_neb_cp2k](https://github.com/arvpinto/enzyme_neb_cp2k)：酶反应 QM/MM 的 CI-NEB、Dimer 与频率验证（CC0）
- [bwang-ecnu/gromacs_cp2k_tutorial](https://github.com/bwang-ecnu/gromacs_cp2k_tutorial)：GROMACS–CP2K QM/MM 示例体系
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含 xTB 与 QM/MM 章节（GPL-3.0，只链接）
- [tblite/tblite](https://github.com/tblite/tblite)：CP2K 中 GFN2-xTB 的后端库

**本仓库脚本**：[md_analyze.py](../scripts/md_analyze.py)（MD 轨迹与能量分析）、[neb_plot.py](../scripts/neb_plot.py)（NEB 能量曲线）、[vib_thermo.py](../scripts/vib_thermo.py)（频率与热化学）、[cp2k_outparse.py](../scripts/cp2k_outparse.py)（输出解析）、[conv_kpoints.py](../scripts/conv_kpoints.py)（k 点收敛）、[job_cp2k_slurm.sh](../scripts/job_cp2k_slurm.sh)（Slurm 作业模板）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.* **13**, 1989 (2017)，GFN1-xTB：<https://doi.org/10.1021/acs.jctc.7b00118>
- C. Bannwarth, S. Ehlert, S. Grimme, *J. Chem. Theory Comput.* **15**, 1652 (2019)，GFN2-xTB：<https://doi.org/10.1021/acs.jctc.8b01176>
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, ChemRxiv (2019)，GFN0-xTB：<https://doi.org/10.26434/chemrxiv.8326202.v1>
- C. Bannwarth *et al.*, *WIREs Comput. Mol. Sci.* **11**, e1493 (2021)，xTB 方法综述：<https://doi.org/10.1002/wcms.1493>
- H. Neugebauer *et al.*, *J. Comput. Chem.* **44**, 2120 (2023)，spGFN-xTB：<https://doi.org/10.1002/jcc.27185>
- M. Elstner *et al.*, *Phys. Rev. B* **58**, 7260 (1998)，SCC-DFTB：<https://doi.org/10.1103/PhysRevB.58.7260>
- L. Zhechkov *et al.*, *J. Chem. Theory Comput.* **1**, 841 (2005)，DFTB 的 UFF 色散校正：<https://doi.org/10.1021/ct050065y>
- T. Laino, F. Mohamed, A. Laio, M. Parrinello, *J. Chem. Theory Comput.* **1**, 1176 (2005)，GEEP 实空间多重网格 QM/MM 静电耦合：<https://doi.org/10.1021/ct050123f>
- T. Laino, F. Mohamed, A. Curioni, M. Parrinello, *J. Chem. Theory Comput.* **2**, 1370 (2006)，周期性 QM/MM：<https://doi.org/10.1021/ct6001169>
- F. Maseras, K. Morokuma, *J. Comput. Chem.* **16**, 1170 (1995)，IMOMM：<https://doi.org/10.1002/jcc.540160911>
- H. M. Senn, W. Thiel, *Angew. Chem. Int. Ed.* **48**, 1198 (2009)，生物体系 QM/MM 综述：<https://doi.org/10.1002/anie.200802019>
- W. L. Jorgensen *et al.*, *J. Chem. Phys.* **79**, 926 (1983)，TIP3P：<https://doi.org/10.1063/1.445869>
- Revisiting the reaction pathways for phospholipid hydrolysis catalyzed by phospholipase A2 with QM/MM methods, *Chem. Sci.* **15**, 9793 (2024)（enzyme_neb_cp2k 的配套论文）：<https://doi.org/10.1039/d4sc02315c>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>

---

[← 上一章](15-metadynamics-plumed.md) | [返回目录](../README.md) | [下一章 →](17-tools-ecosystem.md)
