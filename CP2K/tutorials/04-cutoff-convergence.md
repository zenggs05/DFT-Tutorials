# 第4章 CUTOFF / REL_CUTOFF 收敛测试（Converging CUTOFF & REL_CUTOFF in CP2K）

本章讲 CP2K Quickstep（GPW）多重网格（multigrid）的工作原理，以及 `&MGRID` 中 `CUTOFF`、`REL_CUTOFF`、`NGRIDS`、`PROGRESSION_FACTOR` 各自控制什么。我们以 8 原子体相 Si 和孤立 H2O 分子为例，完整走一遍 CP2K 截断能收敛测试（cutoff convergence test）：先扫 CUTOFF，再扫 REL_CUTOFF。内容还包括如何读 `MULTIGRID INFO` 中的 `count for grid`、`EPS_DEFAULT` 对精度的影响、为什么 CELL_OPT/NPT 需要更高的截断，以及怎样用 `@SET` 模板加脚本批量生成、运行和汇总作图。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清楚 CP2K 的 `CUTOFF` 与平面波程序（VASP 的 `ENCUT`、QE 的 `ecutwfc`）在含义上的区别；
- 理解多重网格的层级划分，以及高斯函数乘积被分配到各级网格的规则；
- 用“两步法”为新体系确定 `CUTOFF` 和 `REL_CUTOFF`，并根据总能量和 `count for grid` 两类判据判断是否收敛；
- 知道 `EPS_DEFAULT`、基组（SR 与非 SR）、赝势硬度对所需截断的影响；
- 为 CELL_OPT / NPT 单独做应力（压力）收敛，并会使用 `CELL_REF`；
- 用模板加 shell 脚本（或本仓库的 `conv_cutoff.py`）批量完成扫描，在本地或 Slurm 集群上运行。

什么时候需要做本章的测试？

| 场景 | 是否需要重新做截断收敛 |
|---|---|
| 第一次计算某种元素组合（尤其含过渡金属、含半芯态赝势，如 Ti-q12、Cu-q11、Pt-q18） | 必须 |
| 更换基组族（MOLOPT → MOLOPT-SR → UZH）或赝势（q 值不同） | 必须 |
| 从单点能/几何优化转向 CELL_OPT、NPT-MD、声子、弹性常数（依赖应力或高阶导数） | 必须，且判据要看应力 |
| 只改变体系大小（同一材料从 2×2×2 扩到 3×3×3 超胞） | 一般不需要，截断是“局域”性质 |
| 换泛函（PBE → PBE0、SCAN） | 建议复查（meta-GGA 对网格更敏感） |

> [!NOTE]
> 本章默认你已经读过 [第2章 输入文件结构](02-input-structure.md) 和 [第3章 GPW/GAPW、基组与赝势](03-gpw-basis-pseudo.md)。示例输入集在 [../examples/04-cutoff-convergence-si/README.md](../examples/04-cutoff-convergence-si/README.md)。

---

## 2. 原理速览

### 2.1 GPW 为什么需要“平面波截断”

GPW（Gaussian and Plane Waves）方法用原子中心的高斯函数展开 Kohn–Sham 轨道，但电子密度 $\rho(\mathbf r)$ 要映射（collocate）到均匀的实空间网格上。Hartree 势在倒空间用 FFT 求解，交换关联势也在网格上计算。网格的疏密由平面波截断 $E_\mathrm{cut}$ 决定：

$$
\frac{1}{2}|\mathbf G|^2 \le E_\mathrm{cut}\ (\text{Ha})
\quad\Longleftrightarrow\quad
|\mathbf G|_\mathrm{max} = \sqrt{E_\mathrm{cut}\,[\text{Ry}]}\ \ \text{bohr}^{-1}
$$

对应的实空间网格间距约为 $h \approx \pi / |\mathbf G|_\mathrm{max}$。例如 `CUTOFF 400`（Ry）给出 $|\mathbf G|_\mathrm{max}=20\ \text{bohr}^{-1}$，网格间距约 0.157 bohr（0.083 Å）。对边长 5.431 Å（10.26 bohr）的 Si 立方胞，每个方向约 65 个点，CP2K 会再向上取到便于 FFT 的整数。网格点总数大致随 $E_\mathrm{cut}^{3/2}$ 增长，所以截断每提高一倍，FFT 和积分的开销约增加 2.8 倍。

> [!IMPORTANT]
> CP2K 的 `CUTOFF` 是**密度网格**的截断，对应 QE 的 `ecutrho`，而不是波函数截断（QE 的 `ecutwfc`、VASP 的 `ENCUT`）。轨道本身由高斯基组描述，基组是否完备与 `CUTOFF` 无关。所以 CP2K 的 400 Ry 不能拿去和平面波程序的 400 eV 或 30 Ry 比较。截断收敛不能代替基组收敛（DZVP → TZV2P），两者要分别测试。

### 2.2 多重网格：NGRIDS 与 PROGRESSION_FACTOR

高斯函数宽窄差别很大。很陡的高斯（如 O 的 1s 型收缩、或指数 $\alpha \sim 10\ \text{bohr}^{-2}$ 的原函数）需要很细的网格；很宽的弥散函数放在细网格上则是浪费。Quickstep 因此同时维护 `NGRIDS` 套网格（默认 4 套），第 $i$ 级的截断为

$$
E_\mathrm{cut}^{(i)} = \frac{E_\mathrm{cut}^{(1)}}{\alpha^{\,i-1}},\qquad i = 1,\dots,N_\mathrm{grids}
$$

其中 $E_\mathrm{cut}^{(1)}$ 就是 `CUTOFF`，$\alpha$ 是 `PROGRESSION_FACTOR`（默认 3.0）。以 `CUTOFF 400` 为例：

| 网格级别 | 截断 / Ry | 截断 / Ha（输出中的 `[a.u.]`） | 用途 |
|---|---|---|---|
| 1（最细） | 400.0 | 200.00 | 最陡的高斯乘积；电子密度最终总是汇总到最细网格 |
| 2 | 133.3 | 66.67 | 中等宽度 |
| 3 | 44.4 | 22.22 | 较宽 |
| 4（最粗） | 14.8 | 7.41 | 最弥散的函数 |

> [!NOTE]
> 输出中 `MULTIGRID INFO` 和 `QS| Density cutoff [a.u.]` 的单位都是 Hartree，数值是你在输入里写的 Ry 值的一半。别把 `200.00` 误读成“截断设成了 200”。

### 2.3 REL_CUTOFF：高斯乘积被分配到哪一级网格

两个原函数 $e^{-\alpha r_A^2}$ 和 $e^{-\beta r_B^2}$ 的乘积仍是高斯，指数为 $\zeta=\alpha+\beta$。CP2K 2026.2 源码（`gaussian_gridlevels.F`）中的分配规则可以写成：

$$
E_\mathrm{need} = \zeta \times \texttt{REL\_CUTOFF},\qquad
\text{取满足 } E_\mathrm{cut}^{(i)} \ge E_\mathrm{need} \text{ 的最粗一级 } i
$$

若连最细网格都不满足，就放在第 1 级。直观理解是：`REL_CUTOFF` 是“指数为 1 的高斯所需的参考截断”。它越大，乘积高斯就越倾向于落到更细的网格，每个高斯覆盖的网格点数越多，积分越精确，代价也越高。

用本仓库 `CP2K/data/BASIS_MOLOPT_UZH` 中 O 的 `DZVP-MOLOPT-PBE-GTH-q6` 估算一下：最大指数约 10.28，最小指数约 0.155（单位 bohr⁻²）。在 `CUTOFF 400`、`REL_CUTOFF 60` 时：

- 最陡乘积 $\zeta \approx 20.6$，$E_\mathrm{need}\approx 1233$ Ry > 400 Ry，放到第 1 级；
- 最宽乘积 $\zeta \approx 0.31$，$E_\mathrm{need}\approx 18.6$ Ry。第 4 级（14.8 Ry）不够，第 3 级（44.4 Ry）满足，所以放到第 3 级。

由此得到两条结论，对理解扫描结果很重要：

1. **REL_CUTOFF 太小**：即使 CUTOFF 很高，大部分乘积也会被“降级”到粗网格，有效积分精度并不高；
2. **只提高 CUTOFF、不提高 REL_CUTOFF**：各级网格整体变细，原本需要第 1 级的乘积会逐渐满足第 2、3 级的条件，最细网格上的乘积数（`count for grid 1`）反而减少，甚至降到 0。此后继续提高 CUTOFF，能量几乎不再变化，但你也没有得到更高的精度。

### 2.4 基组与赝势决定了“需要多高”

所需截断主要由两样东西决定：**最陡的基函数指数**和**赝势局域部分的硬度**。从本仓库 `CP2K/data` 中读出的最大原函数指数（bohr⁻²）可以说明 SR 基组为什么对截断更友好：

| 元素 | 基组 | 最大指数 | 最小指数 | 基函数数 |
|---|---|---|---|---|
| O | DZVP-MOLOPT-GTH（BASIS_MOLOPT） | 12.02 | 0.047 | 13 |
| O | DZVP-MOLOPT-SR-GTH（BASIS_MOLOPT） | 10.39 | 0.162 | 13 |
| O | DZVP-MOLOPT-PBE-GTH-q6（BASIS_MOLOPT_UZH） | 10.28 | 0.155 | 13 |
| Si | DZVP-MOLOPT-GTH | 2.69 | 0.068 | 13 |
| Si | DZVP-MOLOPT-SR-GTH | 1.26 | 0.087 | 13 |
| Si | DZVP-MOLOPT-PBE-GTH-q4（UZH） | 1.23 | 0.048 | 13 |
| Pt | DZVP-MOLOPT-PBE-GTH-q18（UZH） | 5.87 | 0.078 | 26 |

（数据由脚本解析数据文件得到，指数只取各收缩中出现的原函数最大/最小值。）

可以看到：

- **SR（short-range）基组**的最大指数更小，所需 CUTOFF 更低；最小指数更大，即弥散函数更少，重叠矩阵条件数更好，更适合凝聚相；
- 第二周期元素（C、N、O、F）本身指数大，赝势也“硬”，有机分子和水体系往往需要 400–600 Ry；
- 带半芯态的过渡金属（q11、q12、q18）的 GTH 赝势局域部分较硬，所需 CUTOFF 通常也偏高。

> [!TIP]
> 下面给出一个**经验起点**，只用来确定扫描范围，不能代替测试：Si/Al 等 sp 元素 250–350 Ry；含 H/C/N/O 的分子与水 400–600 Ry；带半芯态的过渡金属 500–800 Ry；REL_CUTOFF 绝大多数体系在 50–70 Ry 收敛。

### 2.5 其它影响网格精度的参数：EPS_DEFAULT

`&DFT/&QS/EPS_DEFAULT`（默认 `1.0E-10`）是一个“总开关”，它按比例设定一系列阈值：

| 派生关键字 | 默认值（随 EPS_DEFAULT） | 控制内容 |
|---|---|---|
| `EPS_PGF_ORB` | $\sqrt{\text{EPS\_DEFAULT}}$ | 重叠矩阵元截断、近邻列表范围 |
| `EPS_RHO`（`EPS_RHO_RSPACE`/`EPS_RHO_GSPACE`） | EPS_DEFAULT | 密度映射到网格时高斯的截断半径 |
| `EPS_GVG_RSPACE` | $\sqrt{\text{EPS\_DEFAULT}}$ | 实空间 KS 矩阵元积分精度 |
| `EPS_CORE_CHARGE` | EPS_DEFAULT / 100 | 核电荷映射精度 |

EPS_DEFAULT 太松时，总能量会带一个与网格无关的“噪声底”。这时无论 CUTOFF 设多高，能量差都低不过这个噪声；MD 中还会表现为守恒量漂移。所以做截断收敛时，要把 EPS_DEFAULT 设得比生产计算**至少同样严格**（推荐 `1.0E-12`），否则会把 EPS 噪声误判成网格误差。第 6 节给出单独扫描 EPS_DEFAULT 的方法。

---

## 3. 收敛判据：收敛到多少才算够

CP2K 中能量随 CUTOFF 的变化**不是单调的**：GPW 对网格不满足变分原理，曲线可能先降后升或小幅振荡。因此判据应当是“与参考值（扫描中最高的那个截断）的偏差小于阈值”，而不是“能量一直下降”。

| 判据 | 常用阈值 | 适用 |
|---|---|---|
| 总能量偏差 | < 1×10⁻⁶ Ha/原子（≈0.03 meV/原子） | 严格设置，官方教程的量级 |
| 总能量偏差 | < 1 meV/原子（≈3.7×10⁻⁵ Ha/原子） | 常规结构优化、吸附能、反应能 |
| 原子力偏差 | < 1×10⁻⁴ Ha/bohr（≈5 meV/Å） | 几何优化、振动分析、MD |
| 压力（应力迹/3）偏差 | < 1 kbar（0.1 GPa），精细工作 < 0.5 kbar | CELL_OPT、NPT |
| 网格分布 | 第 1 级 `count for grid` 不为 0，多数乘积在中间级别 | 所有体系 |
| 电子数积分误差 | `Electronic density on regular grids` 第二列 ≲ 1×10⁻⁸ | 辅助判断 |

单位换算：1 Ha = 2 Ry = 27.211386 eV；1 mHa ≈ 27.2 meV；1×10⁻⁶ Ha ≈ 0.0272 meV。

> [!TIP]
> 能量**差**（吸附能、反应能、能垒）因误差抵消通常收敛得更快。但用总能量做判据更保守，换来的是一次测试可用于多种性质。若计算成本吃紧，可以对目标能量差直接做收敛测试。

---

## 4. 实战：8 原子体相 Si 的两步法

示例文件位于 [../examples/04-cutoff-convergence-si/README.md](../examples/04-cutoff-convergence-si/README.md)：

| 文件 | 作用 |
|---|---|
| `template.inp` | 参数化模板，顶部两行 `@SET` 由脚本替换 |
| `si_bulk8.xyz` | 8 原子立方胞坐标（a = 5.431 Å，实验晶格常数） |
| `run_cutoff_scan.sh` | 第一步：固定 REL_CUTOFF 扫 CUTOFF |
| `run_relcutoff_scan.sh` | 第二步：固定 CUTOFF 扫 REL_CUTOFF |
| `collect_scan.sh` | 提取能量、SCF 收敛标志和各级网格计数，写成 CSV |
| `plot_scan.py` | 画 \|ΔE\|（meV/原子，对数坐标）随截断的变化 |
| `h2o_template.inp`、`h2o.xyz` | 分子体系（H2O + MT 求解器）的对应模板，见第 5 节 |
| `job_scan_slurm.sh` | Slurm 作业脚本，在计算节点上依次跑完两步扫描 |

### 4.1 模板输入

```text
# CP2K template for the CUTOFF / REL_CUTOFF convergence test
# System: bulk Si, 8-atom cubic cell (a = 5.431 Angstrom, experimental value)
# The two @SET lines are rewritten by run_cutoff_scan.sh / run_relcutoff_scan.sh.
@SET CUTOFF 300
@SET REL_CUTOFF 60

&GLOBAL
  PROJECT     si8_conv
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM            ! MEDIUM is needed to print the MULTIGRID INFO block
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      NGRIDS     4              ! default, written out for clarity
      CUTOFF     ${CUTOFF}      ! Ry, finest grid
      REL_CUTOFF ${REL_CUTOFF}  ! Ry, reference grid for a Gaussian with exponent 1
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12       ! tighter than the default 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-8          ! SCF noise must be far below the target accuracy
      MAX_SCF   100
      ADDED_MOS 10
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.4
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &PRINT
        &RESTART OFF            ! no .wfn files needed for a scan
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

逐项说明：

| 关键字 | 取值 | 为什么这样设 |
|---|---|---|
| `@SET CUTOFF` / `@SET REL_CUTOFF` | 300 / 60 | 预处理器变量，下文用 `${CUTOFF}` 引用；脚本只改这两行，生成的每个输入都自成一体、可复现 |
| `RUN_TYPE` | `ENERGY` | 收敛测试只需单点能，不算力，成本最低（测力的收敛时改为 `ENERGY_FORCE`） |
| `PRINT_LEVEL` | `MEDIUM` | `LOW` 不打印 `MULTIGRID INFO`，就无法检查网格分布 |
| `NGRIDS` | 4 | 默认值。一般不需要改，见第 8 节 |
| `EPS_DEFAULT` | 1.0E-12 | 比默认值严格，避免 EPS 噪声干扰判断 |
| `EPS_SCF` | 1.0E-8 | SCF 残差带来的能量误差必须远小于 10⁻⁶ Ha/原子 |
| `ADDED_MOS` + `&SMEAR` | 10 + Fermi–Dirac 300 K | 只算 Γ 点的 8 原子小胞，折叠后的 HOMO–LUMO 间隙可能偏小，加少量展宽和空轨道作为保险，让对角化稳定收敛；所有截断点用同样的设置，不影响比较 |
| `&MIXING` | Broyden, α=0.4 | 对角化必须配密度混合，详见 [第5章](05-scf-strategy.md) |
| `&PRINT &RESTART OFF` | — | 扫描几十个点时不写 `.wfn` 文件，节省磁盘 |
| `BASIS_MOLOPT_UZH` + `POTENTIAL_UZH` | DZVP-…-q4 + GTH-PBE-q4 | 新项目推荐的成对组合，q 值必须一致（见[第3章](03-gpw-basis-pseudo.md)） |

> [!NOTE]
> 官方 HowTo（[How to Converge the CUTOFF and REL_CUTOFF](https://manual.cp2k.org/trunk/methods/dft/cutoff.html)）为了省时间使用了 `MAX_SCF 1`，即只做一次对角化、不迭代到自洽，并用 `IGNORE_CONVERGENCE_FAILURE` 让程序继续。理由是网格足够细时，非自洽能量同样应当稳定。自 2024.1 起 SCF 不收敛默认中止，这个技巧必须配合 `IGNORE_CONVERGENCE_FAILURE TRUE` 才能运行。本教程对 Si8 这样的小体系直接做完整 SCF：结果更直观，也不会在模板里留下 `IGNORE_CONVERGENCE_FAILURE`，避免被误抄到生产输入中。体系很大、扫描成本高时，可以改用官方的 `MAX_SCF 1` 方案。

先检查语法再开始扫描：

```bash
cd CP2K/examples/04-cutoff-convergence-si
export CP2K_DATA_DIR=/path/to/cp2k/data     # conda/容器版本通常已内置，可省略
cp2k.psmp -c template.inp                    # 只做语法检查，不计算
cp2k.psmp -e template.inp | grep -A3 "&MGRID"   # 回显展开后的输入，确认 ${CUTOFF} 已被替换
```

### 4.2 第一步：固定 REL_CUTOFF = 60 Ry，扫描 CUTOFF

`run_cutoff_scan.sh` 的核心逻辑如下（完整脚本见示例目录）：

```bash
#!/usr/bin/env bash
# Step 1: scan CUTOFF at fixed REL_CUTOFF (core logic, see the example directory)
set -euo pipefail
TEMPLATE=${TEMPLATE:-template.inp}
EXTRA_FILES=${EXTRA_FILES:-si_bulk8.xyz}
REL_CUTOFF=${REL_CUTOFF:-60}
CUTOFFS=${CUTOFFS:-"150 200 250 300 350 400 450 500 600"}
CP2K=${CP2K:-cp2k.psmp}

for c in $CUTOFFS; do
  d=cutoff_${c}Ry
  mkdir -p "$d"
  sed -E -e "s/^@SET[[:space:]]+CUTOFF[[:space:]].*/@SET CUTOFF ${c}/" \
         -e "s/^@SET[[:space:]]+REL_CUTOFF[[:space:]].*/@SET REL_CUTOFF ${REL_CUTOFF}/" \
         "$TEMPLATE" > "$d/input.inp"
  for f in $EXTRA_FILES; do cp "$f" "$d/"; done
  (cd "$d" && $CP2K -i input.inp -o output.out)
done
bash collect_scan.sh cutoff 8
```

运行方式：

```bash
# serial, one MPI rank per calculation
bash run_cutoff_scan.sh

# 4 MPI ranks per calculation, 2 calculations at the same time (8 cores in total)
CP2K="mpirun -np 4 cp2k.psmp" NPAR=2 bash run_cutoff_scan.sh

# only generate inputs (e.g. to inspect them or submit them separately)
DRYRUN=1 bash run_cutoff_scan.sh
```

扫描范围选 150–600 Ry。最高点（600 Ry）作为参考值，它应明显高于预期的收敛值。如果最后几个点之间仍有可见变化，就继续向上加点（如 700、800）。

### 4.3 读懂输出：总能量与 MULTIGRID INFO

每个目录的 `output.out` 中需要关注以下几处：

```bash
cd cutoff_300Ry
grep "QS| Density cutoff"             output.out   # 确认实际使用的截断（Ha）
grep "QS| Relative density cutoff"    output.out
grep -A8 "MULTIGRID INFO"             output.out   # 各级网格的高斯乘积计数
grep "Electronic density on regular"  output.out   # 网格上积分得到的电子数
grep "SCF run converged"              output.out   # 必须出现
grep "ENERGY| Total FORCE_EVAL"       output.out   # 最终总能量（Ha）
```

`MULTIGRID INFO` 块的格式如下（示意输出，数值仅供参考）：

```text
 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:            912          cutoff [a.u.]          150.00
 count for grid        2:           4416          cutoff [a.u.]           50.00
 count for grid        3:           3528          cutoff [a.u.]           16.67
 count for grid        4:            200          cutoff [a.u.]            5.56
 total gridlevel count  :           9056
```

- 第 5 列是落在该级网格上的高斯乘积个数，最后一列是该级截断（Ha）；
- 第 1 级计数为 0，说明在当前 REL_CUTOFF 下最细网格没有用上，CUTOFF 再高也只是浪费；
- 绝大多数乘积挤在第 3、4 级，说明 REL_CUTOFF 偏小（第二步会测试这一点）。

`Electronic density on regular grids:` 一行有两个数：第一个是网格上积分的电子密度（为负值，CP2K 中电子带负号），第二个是它与精确电子数之差。差值越接近 0，网格对密度的表示越好；若在 10⁻⁵ 量级，就说明截断明显不足。

### 4.4 第一步结果的判读

`collect_scan.sh` 运行后会生成 `cutoff_scan.csv` 并在屏幕上打印表格。下面是一组**示意输出，数值仅供参考**，用来说明怎样读：

```text
cutoff_Ry  energy_Ha        dE_meV_per_atom  scf_converged  grid_counts
150        -31.6573812269   0.8428           yes            3096;4536;1224;200
200        -31.6575932852   0.1215           yes            2048;4192;2616;200
250        -31.6576168929   0.0412           yes            1360;4024;3272;400
300        -31.6576255364   0.0118           yes            912;4416;3528;200
350        -31.6576266535   0.0080           yes            648;3744;4272;392
400        -31.6576279471   0.0036           yes            400;3632;4224;800
450        -31.6576286821   0.0011           yes            304;3104;4640;1008
500        -31.6576288291   0.0006           yes            200;2600;5024;1232
600        -31.6576290055   0.0000           yes            96;2032;5216;1712
```

（示意输出，数值仅供参考；`grid_counts` 列按 `grid1;grid2;grid3;grid4` 排列，`dE` 以表中最高截断为参考。）

判读要点：

1. 看 `dE_meV_per_atom` 从哪个截断开始**稳定**低于你的阈值（如 0.03 meV/原子 ≈ 10⁻⁶ Ha/原子），而不是只看某一个点；
2. 取满足阈值的**最低**截断，再上浮 50 Ry 左右作为安全余量（生产计算中原子位置变化会让网格误差略有波动）；
3. 检查该截断下第 1 级网格计数不为 0。

按上面这组示意数据：在 0.03 meV/原子的阈值下，300 Ry 是第一个稳定满足阈值的点，加安全余量后取 350 Ry；若只要求 1 meV/原子，200 Ry 就已足够。同时注意 `grid_counts` 的变化：随 CUTOFF 升高，第 1 级的计数从 3096 一路降到 96，乘积不断被“推”到较粗的网格，这正是 2.3 节说的现象。

用示例目录里的脚本作图：

```bash
python plot_scan.py cutoff_scan.csv --threshold 0.03 -o cutoff_scan.png
```

（`plot_scan.py` 只依赖 numpy 和 matplotlib，用 `python plot_scan.py -h` 查看参数。本仓库的通用驱动脚本 ../scripts/conv_cutoff.py 也能完成“生成、运行、汇总、作图”全流程，用法见 `python ../scripts/conv_cutoff.py -h`。）

### 4.5 第二步：固定 CUTOFF，扫描 REL_CUTOFF

假设第一步选出 CUTOFF = 350 Ry，把它固定下来扫 REL_CUTOFF：

```bash
CUTOFF=350 bash run_relcutoff_scan.sh
# 默认扫描 20 30 40 50 60 70 80 90 100；可用 RELCUTOFFS="..." 覆盖
```

这时生成的目录名为 `relcutoff_20Ry`、`relcutoff_30Ry` …，结果写入 `relcutoff_scan.csv`。随 REL_CUTOFF 增大，你会看到各级计数**向细网格移动**（第 1、2 级增多，第 4 级减少），能量在某个值之后进入平台。取进入平台的最低 REL_CUTOFF；平台通常出现在 50–70 Ry。

> [!WARNING]
> 如果第二步选出的 REL_CUTOFF 明显高于第一步固定的 60 Ry（如 80 Ry），就用新的 REL_CUTOFF **重做一遍第一步**。两个参数会互相影响：REL_CUTOFF 变了，同一 CUTOFF 下的网格分布和能量都会变。一般迭代一次就足够。

### 4.6 换一种参数注入方式：命令行 `-E/--set`

除了用 sed 改写 `@SET` 行，CP2K 还能在命令行给预处理器变量设**初始值**：

```bash
# template_noset.inp: identical to template.inp but WITHOUT the two @SET lines
cp2k.psmp -i template_noset.inp -o si_400.out -E CUTOFF=400 -E REL_CUTOFF=60
```

注意以下两点：

- `-E name=value`（长写法 `--set`）只设定变量的初始值。如果输入文件里又出现同名的 `@SET`，会被文件中的值覆盖。所以用这种方式时，模板中**不能**再写这两行 `@SET`；
- 该选项在 2023.1 的源码中已经存在，7.1 中还没有。老版本请用 sed 方案。

批量生成参数扫描输入的其它途径：

- `cp2k-input-tools` 的 `cp2kgen` 命令按表达式生成参数扫描输入，见 第17章；
- 本仓库的 ../scripts/conv_cutoff.py（`python ../scripts/conv_cutoff.py -h`）；
- 官方示例 [cp2k/cp2k-examples converging_grid](https://github.com/cp2k/cp2k-examples/tree/master/converging_grid) 中的脚本套件。

### 4.7 在 Slurm 集群上跑扫描

示例目录里的 `job_scan_slurm.sh` 把两步扫描放进一个作业，在同一节点上**依次**运行，每个计算使用作业分到的全部核：

```bash
sbatch job_scan_slurm.sh                 # 用默认的 CUTOFF=350 做第二步
sbatch --export=ALL,CUTOFF2=400 job_scan_slurm.sh   # 第二步改用 400 Ry
```

体系较大、单点要几十分钟时，更好的做法是把每个截断作为 Slurm 作业数组（job array）的一个元素并行提交。思路如下（片段，需要按集群修改）：

```bash
#SBATCH --array=0-8
CUTOFFS=(150 200 250 300 350 400 450 500 600)
c=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}
cd cutoff_${c}Ry && srun cp2k.psmp -i input.inp -o output.out
```

先用 `DRYRUN=1 bash run_cutoff_scan.sh` 生成全部目录，再提交上述数组作业；全部完成后执行 `bash collect_scan.sh cutoff 8` 汇总。Slurm/PBS 的更多写法见 Common 第2章 HPC 作业调度 和本仓库的 ../scripts/job_cp2k_slurm.sh。

---

## 5. 分子体系：H2O 的截断收敛

分子体系与体相 Si 的差别在于三点：

1. O、H 的 GTH 赝势和基组指数都更大，所需 CUTOFF 明显更高；
2. 要用非周期 Poisson 求解器（MT），盒子边长至少是分子电荷分布尺寸的 2 倍（详见 [第6章](06-periodic-surface-kpoints.md)）；
3. 有 HOMO-LUMO 能隙，直接用 OT，无需展宽。

示例目录中的 `h2o_template.inp`：

```text
# CP2K template for the CUTOFF / REL_CUTOFF convergence test of an isolated H2O molecule
@SET CUTOFF 500
@SET REL_CUTOFF 60

&GLOBAL
  PROJECT     h2o_conv
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF     ${CUTOFF}
      REL_CUTOFF ${REL_CUTOFF}
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC       NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   30
      &OT
        MINIMIZER      DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
      &PRINT
        &RESTART OFF
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
&END FORCE_EVAL
```

为了不和 Si 的扫描目录混在一起，在单独的子目录中运行：

```bash
mkdir -p h2o && cp h2o_template.inp h2o.xyz run_cutoff_scan.sh run_relcutoff_scan.sh collect_scan.sh plot_scan.py h2o/
cd h2o
TEMPLATE=h2o_template.inp EXTRA_FILES=h2o.xyz \
  CUTOFFS="300 400 500 600 700 800 1000 1200" bash run_cutoff_scan.sh
```

脚本从 xyz 文件第一行读取原子数（H2O 为 3），用于折算“每原子”能量差。分子体系的总能量差除以 3 个原子后数值偏小，建议同时看**总能量差**是否低于 1×10⁻⁵ Ha（约 0.27 meV）这样的绝对阈值。

> [!NOTE]
> 对含水体系（水盒子 AIMD、水/界面），CP2K 社区常用 400 Ry（DZVP-MOLOPT-SR + GTH-PBE），追求精度时用 500–600 Ry。这些是经验值，**请以本节的测试结果为准**，并在论文中写明测试过程。

---

## 6. EPS_DEFAULT 的单独测试

固定已经收敛的 CUTOFF/REL_CUTOFF，只改变 EPS_DEFAULT，就能看出阈值带来的“噪声底”。在 `template.inp` 中把 `EPS_DEFAULT 1.0E-12` 改成 `EPS_DEFAULT ${EPS}`，并在顶部加 `@SET EPS 1.0E-12`，然后：

```bash
for e in 1.0E-8 1.0E-10 1.0E-12 1.0E-14; do
  d=eps_${e}; mkdir -p "$d"
  sed -E -e "s/^@SET[[:space:]]+CUTOFF[[:space:]].*/@SET CUTOFF 350/" \
         -e "s/^@SET[[:space:]]+EPS[[:space:]].*/@SET EPS ${e}/" template.inp > "$d/input.inp"
  cp si_bulk8.xyz "$d/"
  (cd "$d" && cp2k.psmp -i input.inp -o output.out)
  printf "%-8s %s\n" "$e" "$(grep 'ENERGY| Total FORCE_EVAL' $d/output.out | awk '{print $NF}')"
done
```

推荐值：

| 任务 | EPS_DEFAULT |
|---|---|
| 一般单点、几何优化 | 1.0E-10（默认）至 1.0E-12 |
| AIMD（守恒量要求高）、振动分析 | 1.0E-12 |
| 杂化泛函、应力 / CELL_OPT、需要高精度能量差 | 1.0E-12 至 1.0E-14 |

---

## 7. 力和应力的收敛：CELL_OPT / NPT 需要更高的 CUTOFF

能量收敛并不保证力和应力收敛。应力是能量对晶胞应变的导数，而 GPW 的实空间网格跟着晶胞一起伸缩：

- 截断有限时，网格表示的误差会随晶胞形变而改变，给出一个虚假的附加应力。它与平面波程序中的 Pulay 应力类似，常见的表现是压力偏负、优化后的晶胞偏小；
- CELL_OPT 或 NPT 过程中，晶胞变化到一定程度时网格点数会**跳变**，能量、力和应力随之出现不连续。

因此做 CELL_OPT / NPT 之前，要专门对**压力**做一次截断收敛。在模板中加两处：

```text
&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL      ! compute the analytical stress tensor
  &PRINT
    &STRESS_TENSOR ON           ! print it (default unit: bar)
    &END STRESS_TENSOR
  &END PRINT
  ...
&END FORCE_EVAL
```

同时把 `RUN_TYPE` 改为 `ENERGY_FORCE`：`RUN_TYPE ENERGY` 只算能量、不算力，应力也就不会计算。然后提取压力：

```bash
for d in cutoff_*Ry; do
  p=$(grep "STRESS| 1/3 Trace" $d/output.out | tail -1 | awk '{print $NF}')
  echo "${d#cutoff_} $p"
done | sort -n
```

输出中 `STRESS| Analytical stress tensor [bar]` 一行标明了应力是解析计算的。如果显示的是 `Numerical stress tensor`，说明当前方法组合不支持解析应力，CP2K 改用了有限差分，此时 CELL_OPT 会慢很多。

经验上，应力收敛需要的 CUTOFF 往往是能量收敛值的 1.5–2 倍。对含 O/H 的体系，常见做法是 CELL_OPT/NPT 用 ≥ 600 Ry，REL_CUTOFF 60–80 Ry。这只是经验范围，最终请以“压力随 CUTOFF 的变化 < 1 kbar”为准。

另一个关键设置是 `&CELL_REF`：在其中给出一个**比优化过程中可能出现的最大晶胞还要大一些**的参考晶胞，网格点数就按参考晶胞固定，CELL_OPT 中不会再跳变：

```text
&SUBSYS
  &CELL
    ABC 5.431 5.431 5.431
    &CELL_REF
      ABC 5.70 5.70 5.70        ! ~5% larger than the expected maximum cell
    &END CELL_REF
  &END CELL
  ...
&END SUBSYS
```

CELL_OPT 的完整流程（`EXTERNAL_PRESSURE`、`PRESSURE_TOLERANCE`，以及 2026.2 删除 `&CELL_OPT/TYPE` 的破坏性变更）见 [第7章](07-geo-cell-opt.md)；NPT-MD 见 第12章。

---

## 8. NGRIDS、PROGRESSION_FACTOR 与手动多重网格

绝大多数情况下保持默认（`NGRIDS 4`、`PROGRESSION_FACTOR 3.0`）即可，只调 CUTOFF 和 REL_CUTOFF。下面这些选项了解即可：

| 关键字（`&MGRID`） | 默认 | 说明 |
|---|---|---|
| `NGRIDS` | 4 | 网格级数。级数少，乘积更多地落在细网格上，更贵；级数多，最粗一级可能过粗 |
| `PROGRESSION_FACTOR` | 3.0 | 相邻两级的截断比 |
| `MULTIGRID_SET` | F | 设为 T 后可用 `MULTIGRID_CUTOFF`（别名 `CUTOFF_LIST`）手动给出每一级的截断 |
| `COMMENSURATE` | F | 让各级网格可公度；为 T 时会覆盖 PROGRESSION_FACTOR 和子网格截断 |
| `&RS_GRID` 子段 | — | 实空间网格的并行分布方式，影响性能不影响精度 |

> [!TIP]
> 如果第一步扫描中，`count for grid 1` 在你选的 CUTOFF 下已经是 0，这不是去调 NGRIDS 的理由。正确的做法是提高 REL_CUTOFF，或者接受一个更低、但第 1 级仍有乘积的 CUTOFF。

---

## 参数选择与收敛建议

| 参数 | 推荐起点 | 收敛判据 | 备注 |
|---|---|---|---|
| `CUTOFF` | sp 元素 250–350 Ry；H/C/N/O 400–600 Ry；半芯过渡金属 500–800 Ry | ΔE < 10⁻⁶ Ha/原子（严格）或 < 1 meV/原子（常规） | 确定后再上浮约 50 Ry |
| `REL_CUTOFF` | 60 Ry | 能量进入平台，`count for grid 1` 不为 0 | 若与 60 差别大，回头重做第一步 |
| `NGRIDS` | 4 | — | 一般不改 |
| `EPS_DEFAULT` | 1.0E-12（测试与 MD） | 噪声底远低于截断误差 | 杂化、应力用 1.0E-12 至 1.0E-14 |
| `EPS_SCF`（测试用） | 1.0E-8（对角化）/ 1.0E-7（OT） | SCF 噪声 ≪ 目标阈值 | 测试结束后恢复生产值 |
| CELL_OPT / NPT | 能量收敛值的 1.5–2 倍，常见 ≥ 600 Ry | 压力变化 < 1 kbar | 加 `STRESS_TENSOR ANALYTICAL` 与 `&CELL_REF` |
| 基组 | 凝聚相优先 DZVP-MOLOPT-SR 或 UZH | — | 换基组后重做截断测试 |

---

## 常见问题与报错

**1. 输出中找不到 `MULTIGRID INFO`。**
原因：`PRINT_LEVEL LOW`（或 `SILENT`）。解决：在 `&GLOBAL` 中设 `PRINT_LEVEL MEDIUM`。

**2. 扫描中某些点报 `SCF run NOT converged` 并中止。**
原文（2024.1 起的默认行为）：

```text
 *  \___/     SCF run NOT converged. To continue the calculation regardless,   *
 *    |             please set the keyword IGNORE_CONVERGENCE_FAILURE.         *
```

原因：低截断（如 100–150 Ry）时网格误差大，SCF 难收敛；或 MAX_SCF 太小。解决：提高 `MAX_SCF`，金属或小带隙体系加展宽（见 [第5章](05-scf-strategy.md)），或者直接去掉过低的截断点。`collect_scan.sh` 会跳过没有最终能量的目录（并在屏幕上提示）；若某个点有最终能量却没有 `SCF run converged`（例如用了 `IGNORE_CONVERGENCE_FAILURE`），该行的 `scf_converged` 列会标为 `NO`，这样的点不能用于判断。

**3. 能量随 CUTOFF 不单调，甚至在高截断处“回升”。**
原因：GPW 对网格不满足变分原理，这是正常现象。解决：按“与参考值的偏差”判断，不要按“是否单调下降”判断；同时确认 EPS_DEFAULT 足够严格，不要把 EPS 噪声误当成网格误差。

**4. 子目录中运行时报基组或赝势找不到：**

```text
The requested basis set <DZVP-MOLOPT-PBE-GTH-q4> for element <Si> was not found in the basis set files BASIS_MOLOPT_UZH
```

原因：没有设置 `CP2K_DATA_DIR`，或者基组名拼错、与文件不对应。解决：`export CP2K_DATA_DIR=<cp2k>/data`；或者在 `BASIS_SET_FILE_NAME` 中写绝对路径；并核对基组名是否确实在 `BASIS_MOLOPT_UZH` 中（`grep "DZVP-MOLOPT-PBE-GTH-q4" $CP2K_DATA_DIR/BASIS_MOLOPT_UZH`）。

**5. `count for grid 1` 在所有截断下都是 0。**
原因：REL_CUTOFF 过小，或者体系只有很软的赝势和基组。解决：先做第二步（扫 REL_CUTOFF），再回到第一步。

**6. 预处理器报变量未定义：**
`${CUTOFF}` 没有对应的 `@SET`，也没有在命令行用 `-E CUTOFF=...` 给出。检查模板顶部的 `@SET` 行有没有被删掉或改名。另外，`-E` 只设定初始值，文件中的 `@SET` 会覆盖它（见 4.6 节）。

**7. CELL_OPT 中晶胞持续收缩，压力始终偏负。**
原因：截断不足引起的虚假应力。解决：按第 7 节对压力做截断收敛，提高 CUTOFF，加 `&CELL_REF`，并确认输出中是 `Analytical stress tensor`。

**8. 不同截断下算出的能量拿来做差（如吸附能）。**
这是错误做法：能量差中的各项必须使用**完全相同**的 CUTOFF、REL_CUTOFF、基组、赝势和 EPS 设置。截断测试的结论只用来确定统一的参数，不能混用。

---

## 练习

1. 用 `run_cutoff_scan.sh` 在 REL_CUTOFF = 40、60、80 Ry 下各扫一次 CUTOFF，画在同一张图上，观察 REL_CUTOFF 对收敛速度的影响。
2. 把 `template.inp` 中的基组换成 `BASIS_MOLOPT` 文件里的 `DZVP-MOLOPT-GTH`（赝势用 `GTH_POTENTIALS` 中的 `GTH-PBE-q4`），重复第一步。和 UZH 基组比较，收敛截断有何不同？对照 2.4 节表中的最大指数解释原因。
3. 把 `h2o_template.inp` 的 `RUN_TYPE` 改为 `ENERGY_FORCE`，并在 `&FORCE_EVAL` 中加 `&PRINT` / `&FORCES ON` 打印段，扫描 CUTOFF；提取 O 原子的力（`grep -A5 "ATOMIC FORCES"`），比较力与能量各自收敛所需的截断。
4. 对 Si8 打开 `STRESS_TENSOR ANALYTICAL`（同时把 `RUN_TYPE` 改为 `ENERGY_FORCE`），扫描压力随 CUTOFF 的变化，求出压力收敛到 1 kbar 以内所需的截断，并与能量收敛的截断比较。

---

## 参考资料

**官方文档**

- How to Converge the CUTOFF and REL_CUTOFF：<https://manual.cp2k.org/trunk/methods/dft/cutoff.html>
- 输入参考 `&MGRID`（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/MGRID.html>
- 输入参考 `&QS`（EPS_DEFAULT 等）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/QS.html>
- Basis Sets：<https://manual.cp2k.org/trunk/methods/dft/basis_sets.html>
- Pseudopotentials：<https://manual.cp2k.org/trunk/methods/dft/pseudopotentials.html>
- Geometry and Cell Optimization（CELL_REF、应力）：<https://manual.cp2k.org/trunk/methods/optimization/geometry_and_cell_opt.html>
- 旧版 HowTo 合集：<https://www.cp2k.org/howto>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库，`converging_grid` 是截断收敛脚本的标准范例
- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/` 中的 regtest 是查关键字真实写法的权威来源，`data/` 中有基组和赝势文件
- [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools)：`cp2kgen` 批量生成参数扫描输入，`cp2klint` 检查语法
- [EPCCed/intro_to_cp2k](https://github.com/EPCCed/intro_to_cp2k)：ARCHER2 入门课，含截断收敛练习
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含截断与周期 Si 章节（2025.2 验证）
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，CP2K 入门与收敛测试
- [compchem-cybertraining/Tutorials_CP2K](https://github.com/compchem-cybertraining/Tutorials_CP2K)：以钙钛矿为例的收敛测试（版本较老，部分关键字过时）
- [Clementine-C/cp2kbuild](https://github.com/Clementine-C/cp2kbuild)：中文工具，可生成 CUTOFF/REL_CUTOFF 收敛脚本
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：cutoff/relcutoff/k 点收敛脚本思路（无 license，仅供参考）
- [pierre-24/cp2k-basis](https://github.com/pierre-24/cp2k-basis)：按元素查询匹配的基组和赝势

**关键文献**

- G. Lippert, J. Hutter, M. Parrinello, *Mol. Phys.* **92**, 477 (1997)，GPW 方法：<https://doi.org/10.1080/002689797170220>
- J. VandeVondele *et al.*, *Comput. Phys. Commun.* **167**, 103 (2005)，Quickstep：<https://doi.org/10.1016/j.cpc.2004.12.014>
- J. VandeVondele, J. Hutter, *J. Chem. Phys.* **127**, 114105 (2007)，MOLOPT 基组：<https://doi.org/10.1063/1.2770708>
- S. Goedecker, M. Teter, J. Hutter, *Phys. Rev. B* **54**, 1703 (1996)，GTH 赝势：<https://doi.org/10.1103/PhysRevB.54.1703>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>
- The CP2K Program Package Made Simple, *J. Phys. Chem. B* (2025)：<https://doi.org/10.1021/acs.jpcb.5c05851>

---

[← 上一章](03-gpw-basis-pseudo.md) | [返回目录](../README.md) | [下一章 →](05-scf-strategy.md)
