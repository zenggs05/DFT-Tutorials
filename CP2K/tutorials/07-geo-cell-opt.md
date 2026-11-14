# 第7章 几何优化与晶胞优化（Geometry & Cell Optimization with CP2K）

本章讲 CP2K 的几何优化（`RUN_TYPE GEO_OPT`）和晶胞优化（`RUN_TYPE CELL_OPT`）：BFGS、LBFGS、CG 三种优化器怎么选，`MAX_FORCE`、`RMS_FORCE`、`MAX_DR`、`RMS_DR` 收敛判据的单位和含义，表面 slab 计算中用 `&FIXED_ATOMS` 固定原子，以及晶胞优化中的 `STRESS_TENSOR ANALYTICAL`、`EXTERNAL_PRESSURE`、`KEEP_ANGLES`、`KEEP_SYMMETRY`、`KEEP_SPACE_GROUP`、`KEEP_VOLUME` 和 `CONSTRAINT`。本章同时整理了 CP2K 2026.2 的破坏性变更（删除 `&CELL_OPT/TYPE`）和新功能（定体积晶胞优化、优化结束后输出 CIF/EXTXYZ），并演示怎样监控、续算和检查一次 geometry optimization / cell optimization。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出孤立分子、表面 slab、体相晶体三类体系的 GEO_OPT / CELL_OPT 输入，并说明每个关键字的作用；
- 看懂输出中的 `OPT|` 收敛块，判断优化是否真正收敛；
- 按体系大小和初始结构质量选择 BFGS、LBFGS 或 CG，知道振荡或不收敛时先调哪个参数；
- 用 `&MOTION/&CONSTRAINT/&FIXED_ATOMS` 固定 slab 底层原子，并能验证约束确实生效；
- 为晶胞优化准备足够高的 `CUTOFF`、解析应力和 `CELL_REF`，理解 Pulay 应力的来源；
- 用 `KEEP_ANGLES`、`KEEP_SYMMETRY`、`CONSTRAINT` 等控制晶胞自由度，正确处理层状材料和二维材料；
- 把旧版（≤ 2026.1）输入迁移到 2026.2，并用 `-1.restart`、`-RESTART.wfn` 续算中断的优化。

GEO_OPT 与 CELL_OPT 的选择：

| 体系 / 问题 | 推荐 RUN_TYPE | 说明 |
|---|---|---|
| 孤立分子、团簇（大盒子、`PERIODIC NONE`） | `GEO_OPT` | 盒子不是物理量，不能优化 |
| 表面 slab、吸附构型、界面（晶格常数已由体相确定） | `GEO_OPT` + `&FIXED_ATOMS` | 面内晶格取体相优化值，真空方向不优化 |
| 体相晶体的平衡晶格常数、体积、c/a | `CELL_OPT` | 同时优化原子坐标和晶胞 |
| 层状材料（石墨、MoS2、h-BN）层间距 | `CELL_OPT` + 色散校正 | 没有色散校正时层间距严重偏大，见 [第8章](08-dispersion-spin-dftu.md) |
| 二维单层（真空方向不能动） | `CELL_OPT` + `CONSTRAINT Z` 或 `KEEP_*` | 只优化面内晶格 |
| 有限温度下的平均晶格（热膨胀） | NPT 分子动力学 | 优化得到的是势能面极小，不含零点能和热效应，见 [第12章](12-aimd.md) |
| 过渡态 | `GEO_OPT` + `TYPE TRANSITION_STATE` 或 NEB | 见 [第13章](13-neb-dimer.md) |

本章示例输入集：[../examples/07-geo-cell-opt/README.md](../examples/07-geo-cell-opt/README.md)。

> [!NOTE]
> 本章默认你已经掌握 [第5章 SCF 策略](05-scf-strategy.md)（OT、对角化、展宽）和 [第6章 周期性体系与 k 点](06-periodic-surface-kpoints.md)（`&CELL`、`PERIODIC`、`&KPOINTS`、偶极校正）。优化中的每一步都是一次完整的 SCF，SCF 设置不稳，优化一定不稳。

---

## 2. 原理速览

### 2.1 几何优化在做什么

在 Born–Oppenheimer 近似下，原子核坐标 $\mathbf R$ 对应一个势能 $E(\mathbf R)$，原子受力为 $\mathbf F_I = -\partial E/\partial \mathbf R_I$。几何优化寻找满足 $\mathbf F_I \approx 0$ 的驻点（通常是局部极小）。准牛顿法每一步按

$$
\Delta \mathbf R = -\mathbf H^{-1}\,\mathbf g ,\qquad \mathbf g = -\mathbf F
$$

更新坐标，其中 $\mathbf H$ 是 Hessian（能量二阶导数矩阵）的近似。BFGS 用相邻两步的位移 $\mathbf s_k$ 和梯度差 $\mathbf y_k$ 迭代修正 $\mathbf H$：

$$
\mathbf H_{k+1} = \mathbf H_k + \frac{\mathbf y_k \mathbf y_k^{\mathsf T}}{\mathbf y_k^{\mathsf T}\mathbf s_k} - \frac{\mathbf H_k \mathbf s_k \mathbf s_k^{\mathsf T}\mathbf H_k}{\mathbf s_k^{\mathsf T}\mathbf H_k \mathbf s_k}
$$

CP2K 的 BFGS 默认用 Lindh 模型 Hessian 作初值（`USE_MODEL_HESSIAN T`），并用信赖半径（`TRUST_RADIUS`）限制单步位移。LBFGS 只保存最近若干步的 $(\mathbf s,\mathbf y)$ 对，不存完整 Hessian，内存随原子数线性增长。CG（共轭梯度）沿搜索方向做一维线搜索，每一步可能要好几次能量和力计算，但对噪声大、初始结构差的情况更稳。

> [!NOTE]
> “0 K 下的优化”这个说法并不严格：优化的是电子总能量（势能面），不包含零点振动能和热修正。需要这些量时请看 [第14章 振动分析与热化学](14-vibration-thermo.md)。

### 2.2 晶胞优化与应力

晶胞优化把晶胞矩阵 $\mathbf h=[\mathbf a,\mathbf b,\mathbf c]$ 也作为变量。驱动力是应力张量

$$
\sigma_{\alpha\beta} = \frac{1}{V}\frac{\partial E}{\partial \varepsilon_{\alpha\beta}},\qquad
P_\mathrm{int} = -\frac{1}{3}\mathrm{Tr}\,\boldsymbol\sigma
$$

（符号约定：压缩状态下内压为正）。CELL_OPT 的目标是让内压（以及完整应力张量）与外压 `EXTERNAL_PRESSURE` 平衡，同时原子受力趋于零。CP2K 2026.2 的晶胞优化总是“直接法”（direct cell optimization）：每一步同时更新原子坐标和晶胞，每一步都计算应力。

### 2.3 Pulay 应力：为什么 CELL_OPT 要提高 CUTOFF

CP2K 的轨道用原子中心的高斯基组展开，密度和势放在实空间网格上。晶胞变化时，网格点数是离散的，网格间距跟着晶胞一起变化；基组也不是完备的。两者都会在应力中引入与“体积变化方向”相关的系统误差，统称 Pulay 应力，表现为晶胞系统性地偏小或偏大。对策有三条：

1. `CUTOFF` 比单点能计算明显提高（经验上 ≥ 500–600 Ry，用应力而不是总能量做收敛判据，见 [第4章](04-cutoff-convergence.md)）；
2. 设置 `&CELL_REF`，让整个优化过程使用固定的参考晶胞来确定网格点数；
3. 必要时换更大的基组（DZVP → TZV2P）复核晶格常数。

### 2.4 收敛判据与单位

CP2K 的判据用原子单位，换算成常用单位后更直观：

| 关键字 | 默认值 | 单位 | 换算 | 含义 |
|---|---|---|---|---|
| `MAX_FORCE` | 4.5E-4 | Ha/bohr | ≈ 0.023 eV/Å | 最大力分量 |
| `RMS_FORCE` | 3.0E-4 | Ha/bohr | ≈ 0.015 eV/Å | 力的均方根 |
| `MAX_DR` | 3.0E-3 | bohr | ≈ 0.0016 Å | 相邻两步最大位移 |
| `RMS_DR` | 1.5E-3 | bohr | ≈ 0.0008 Å | 位移均方根 |
| `PRESSURE_TOLERANCE`（仅 CELL_OPT） | 100 | bar | 0.01 GPa | 内压与外压之差 |

1 Ha/bohr ≈ 51.42 eV/Å。作为对照，VASP 常用的 `EDIFFG = -0.02`（eV/Å）约等于 3.9E-4 Ha/bohr，与 CP2K 默认的 `MAX_FORCE` 接近。BFGS 和 CG 要求**全部**判据同时满足才算收敛；LBFGS 另有自己的判据，见第 4.2 节。

---

## 3. 第一个几何优化：孤立水分子（BFGS）

### 3.1 输入文件

示例从一个故意扭曲的水分子开始（O–H 为 1.05/1.00 Å，键角 120°），文件 [../examples/07-geo-cell-opt/h2o_geo_opt_bfgs.inp](../examples/07-geo-cell-opt/h2o_geo_opt_bfgs.inp)：

```text
@SET PROJ h2o_opt
@SET OPT  BFGS

&GLOBAL
  PROJECT ${PROJ}
  RUN_TYPE GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
      ! re-use the wavefunctions of previous steps as SCF guess
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 30
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
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
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME h2o_distorted.xyz
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER ${OPT}
    MAX_ITER 200
    MAX_FORCE 4.5E-4
    RMS_FORCE 3.0E-4
    MAX_DR 3.0E-3
    RMS_DR 1.5E-3
    &BFGS
      TRUST_RADIUS [angstrom] 0.25
      USE_MODEL_HESSIAN T
    &END BFGS
  &END GEO_OPT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        GEO_OPT 1
      &END EACH
    &END TRAJECTORY
    &RESTART
      BACKUP_COPIES 1
      &EACH
        GEO_OPT 1
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
```

（示例目录中的文件还额外包含 `&LBFGS`、`&CG` 子段，方便切换优化器时直接使用。）

### 3.2 关键字逐项说明

| 关键字 | 本例取值 | 说明 |
|---|---|---|
| `RUN_TYPE` | `GEO_OPT` | 固定晶胞，只优化原子坐标 |
| `&QS/EXTRAPOLATION` | `ASPC`（默认） | 用前几步的波函数外推出本步的 SCF 初猜，可显著减少 SCF 步数；`EXTRAPOLATION_ORDER` 默认 3 |
| `EPS_SCF` | 1.0E-7 | 比单点能更紧。力的误差大致与 $\sqrt{\text{EPS\_SCF}}$ 同量级，SCF 太松会让优化在终点附近来回振荡 |
| `&OUTER_SCF` | 10 次外循环 | OT 的标准配置，见[第5章](05-scf-strategy.md) |
| `TYPE` | `MINIMIZATION`（默认） | 找极小；`TRANSITION_STATE` 用于 Dimer 等过渡态搜索 |
| `OPTIMIZER` | `BFGS`（默认） | 别名 `MINIMIZER`；可选 `BFGS`、`LBFGS`、`CG` |
| `MAX_ITER` | 200（默认） | 优化步数上限；CG 和 LBFGS 的一步可能包含多次力计算 |
| `MAX_FORCE` 等 | 默认值 | 见 2.4 节；要接振动分析时把 `MAX_FORCE` 收紧到 1.0E-4 |
| `&BFGS/TRUST_RADIUS` | 0.25 Å（默认） | 单步最大位移；振荡或“越走越远”时首先减小它 |
| `&BFGS/USE_MODEL_HESSIAN` | `T`（默认） | 用 Lindh 模型 Hessian 作为初值，通常比单位矩阵收敛快 |
| `&MOTION/&PRINT/&TRAJECTORY` | 每步输出 | 写 `<PROJECT>-pos-1.xyz`；`FORMAT` 可选 `XYZ`、`PDB`、`DCD`、`EXTXYZ` 等 |
| `&MOTION/&PRINT/&RESTART` | 每步输出 | 写 `<PROJECT>-1.restart`；`BACKUP_COPIES` 控制保留几份 `.bak-*` |

> [!TIP]
> 孤立分子用 `PERIODIC NONE` 加 `POISSON_SOLVER MT`，并用 `&CENTER_COORDINATES` 把分子放到盒子中心。MT 求解器要求盒子边长约为电荷分布尺寸的 2 倍，水分子用 10 Å 足够。细节见 [第6章](06-periodic-surface-kpoints.md)。

### 3.3 运行

```bash
cd CP2K/examples/07-geo-cell-opt
cp2k.psmp -c h2o_geo_opt_bfgs.inp              # syntax check only
mpirun -np 4 cp2k.psmp -i h2o_geo_opt_bfgs.inp -o h2o_geo_opt_bfgs.out
# or on a Slurm cluster:
sbatch --export=ALL,INP=h2o_geo_opt_bfgs.inp job_slurm.sh
```

> [!WARNING]
> 命令行的 `-E VAR=value`（`--set`）只是给预处理变量设**初值**，输入文件里的 `@SET` 会把它覆盖掉（2026.2 源码中 `@SET` 对已有变量直接重新赋值）。想切换优化器，请修改文件中的 `@SET` 行，例如：
>
> ```bash
> sed -e 's/^@SET OPT .*/@SET OPT LBFGS/' -e 's/^@SET PROJ .*/@SET PROJ h2o_lbfgs/' \
>     h2o_geo_opt_bfgs.inp > h2o_lbfgs.inp
> ```

### 3.4 读懂 OPT| 收敛块

每一步优化后，CP2K 打印一个以 `OPT|` 开头的信息块（字段名称来自 2026.2 源码，数值为示意）：

```text
 OPT| **************************************************************************
 OPT| Step number                                                              3
 OPT| Optimization method                                                   BFGS
 OPT| Total energy [hartree]                                      -17.2034517782
 OPT| Effective energy change [hartree]                            -0.0004231158
 OPT| Predicted energy change [hartree]                            -0.0003987413
 OPT| Step size                                                     0.0412776510
 OPT| Trust radius                                                  0.4724315332
 OPT| Decrease in energy                                                     YES
 OPT|
 OPT| Maximum step size                                             0.0302114781
 OPT| Convergence limit for maximum step size                       0.0030000000
 OPT| Maximum step size is converged                                          NO
 OPT|
 OPT| RMS step size                                                 0.0178562304
 OPT| Convergence limit for RMS step size                           0.0015000000
 OPT| RMS step size is converged                                              NO
 OPT|
 OPT| Maximum gradient                                              0.0061873517
 OPT| Convergence limit for maximum gradient                        0.0004500000
 OPT| Maximum gradient is converged                                           NO
 OPT|
 OPT| RMS gradient                                                  0.0035718864
 OPT| Convergence limit for RMS gradient                            0.0003000000
 OPT| RMS gradient is converged                                               NO
 OPT| **************************************************************************
```

（示意输出，数值仅供参考）

要点：

- `Maximum gradient` 就是最大力分量（梯度 = 负的力），与 `MAX_FORCE` 比较；`Maximum step size` 与 `MAX_DR` 比较，单位都是原子单位；
- `Decrease in energy` 偶尔为 `NO` 是正常的，判断收敛看的是力和位移，不是能量；
- `Trust radius` 一行在输出中以 bohr 为单位，所以和输入里写的 0.25 Å 数值不同；
- 四个 `is converged` 都为 `YES`（CELL_OPT 再加上 `Pressure is converged`）时，CP2K 打印：

```text
 *******************************************************************************
 ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
 *******************************************************************************

                    Reevaluating energy at the minimum
```

随后在最终结构上再做一次能量计算。如果到达 `MAX_ITER` 仍未收敛，则打印 `*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***`。

常用检查命令：

```bash
grep -c "OPT| Step number" h2o_geo_opt_bfgs.out            # number of optimization steps
grep "OPT| Maximum gradient  " h2o_geo_opt_bfgs.out | tail -3
grep "GEOMETRY OPTIMIZATION COMPLETED" h2o_geo_opt_bfgs.out
grep "ENERGY| Total FORCE_EVAL" h2o_geo_opt_bfgs.out | tail -1
grep -i "SCF run NOT converged" h2o_geo_opt_bfgs.out       # must be empty
```

也可以用本仓库的脚本解析并作图：`python ../../scripts/cp2k_outparse.py -h` 和 `python ../../scripts/plot_opt.py -h` 查看用法（从示例目录调用时的相对路径）。

### 3.5 如何判断结果合理

- 输出中出现 `GEOMETRY OPTIMIZATION COMPLETED`，且最后几步的 SCF 都正常收敛；
- 用可视化软件（VESTA、Ovito、ASE GUI：`ase gui h2o_opt-pos-1.xyz`）播放轨迹，结构变化平滑，没有原子“飞走”；
- PBE 下水分子的 O–H 键长约 0.97 Å、键角约 104°（近似值；实验值 0.957 Å、104.5°，PBE 通常略高估键长）。可以用 ASE 读最后一帧：

```bash
python - <<'EOF'
from ase.io import read
a = read("h2o_opt-pos-1.xyz", index=-1)
print("d(O-H) =", a.get_distance(0, 1), a.get_distance(0, 2), "A")
print("angle  =", a.get_angle(1, 0, 2), "deg")
EOF
```

---

## 4. 优化器选择：BFGS、LBFGS 与 CG

### 4.1 对比与推荐

| 优化器 | 适用体系 | 优点 | 缺点 | 关键子段 |
|---|---|---|---|---|
| `BFGS`（默认） | 小到中等体系（经验上几百个原子以内） | 收敛快，步数少 | 需要存储并对角化完整 Hessian（$3N\times 3N$），大体系内存和时间开销高 | `&BFGS`：`TRUST_RADIUS`、`USE_MODEL_HESSIAN`、`RESTART_HESSIAN` |
| `LBFGS` | 大体系（上千原子）、凝聚相 | 内存小 | 有自己的收敛判据，可能在一般判据未满足时提前结束 | `&LBFGS`：`MAX_H_RANK`、`TRUST_RADIUS`、`WANTED_PROJ_GRADIENT`、`WANTED_REL_F_ERROR` |
| `CG` | 初始结构很差、力有噪声、准牛顿法不稳定 | 稳健 | 每步要做线搜索，力计算次数多 | `&CG`：`MAX_STEEP_STEPS`、`&LINE_SEARCH TYPE 2PNT/GOLD/FIT` |

经验流程：

```text
初始结构合理？ ── 否 ──> 先用 CG（或 xTB/力场预优化），力降到 ~1e-2 Ha/bohr 以下
      │是
原子数 < 几百？ ── 是 ──> BFGS（默认设置）
      │否
      └──────────────> LBFGS（MAX_H_RANK 适当调大）
振荡 / 能量反复上升 ──> 减小 TRUST_RADIUS（BFGS 0.25 -> 0.1 A），检查 SCF 精度
```

### 4.2 常用子段写法

BFGS（默认已开启模型 Hessian）：

```text
&GEO_OPT
  OPTIMIZER BFGS
  &BFGS
    TRUST_RADIUS [angstrom] 0.10   ! smaller steps for floppy / oscillating systems
    USE_MODEL_HESSIAN T            ! default; F = start from a unit matrix
  &END BFGS
&END GEO_OPT
```

BFGS 会把 Hessian 写到 `<PROJECT>-BFGS.Hessian`，续算时可以读回来，省去重新“学习”Hessian 的步数：

```text
&BFGS
  RESTART_HESSIAN T
  RESTART_FILE_NAME h2o_opt-BFGS.Hessian
&END BFGS
```

LBFGS：

```text
&GEO_OPT
  OPTIMIZER LBFGS
  &LBFGS
    MAX_H_RANK 20        ! default 5; larger = better Hessian approximation, more memory
    MAX_F_PER_ITER 20    ! max force evaluations per line search (default 20)
  &END LBFGS
&END GEO_OPT
```

> [!WARNING]
> LBFGS 有两个专用判据 `WANTED_PROJ_GRADIENT` 和 `WANTED_REL_F_ERROR`，它们可以“覆盖”一般判据。官方文档指出，LBFGS 可能打印 `Specific L-BFGS convergence criteria ... satisfied .... run CONVERGED!` 而一般判据（`MAX_FORCE` 等）并未全部满足。看到这段提示时，请检查最后一个 `OPT|` 块；不满足就收紧判据或换 BFGS/CG 续算。LBFGS 还会在工作目录写一个 `iterate.dat`，其中 `projg` 和 `f` 两列对应这两个专用判据。2026.2 新增了控制 LBFGS 打印量的 `&LBFGS/PRINT_LEVEL`。

CG：

```text
&GEO_OPT
  OPTIMIZER CG
  &CG
    MAX_STEEP_STEPS 5        ! a few steepest-descent steps first (default 0)
    &LINE_SEARCH
      TYPE 2PNT              ! cheapest; GOLD (default) is more robust; FIT is most robust
    &END LINE_SEARCH
  &END CG
&END GEO_OPT
```

| `LINE_SEARCH/TYPE` | 特点 | 何时用 |
|---|---|---|
| `2PNT` | 两点外推，最便宜 | 势能面光滑、力精度高 |
| `GOLD`（默认） | 黄金分割搜索，较贵但稳 | 一般困难情形 |
| `FIT` | 多点拟合抛物线，最贵，对数值噪声最不敏感 | 力或能量有噪声（如 SCF 精度有限） |

### 4.3 其它控制手段

- **紧急停止**：官方文档说明，在工作目录创建名为 `EXIT` 或 `EXIT_GEO` 的文件可以让优化在当前步之后体面退出（例如 `touch EXIT_GEO`），重启文件会正常写出；
- **墙钟时间**：`&GLOBAL/WALLTIME` 到时后 CP2K 主动结束并写重启文件，建议设得比作业队列时限短 10–20 分钟；
- **步数估算**：先跑一次 `RUN_TYPE ENERGY_FORCE` 计时，乘以预计步数（一般几十到两三百步）即可估算总耗时。

---

## 5. 约束：固定原子与表面 slab 优化

### 5.1 `&FIXED_ATOMS` 语法

```text
&MOTION
  &CONSTRAINT
    &FIXED_ATOMS
      COMPONENTS_TO_FIX XYZ      ! X, Y, Z, XY, XZ, YZ or XYZ (default)
      LIST 1..16                 ! atom indices in the order of &COORD / the xyz file
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

| 关键字 | 说明 |
|---|---|
| `LIST` | 原子序号（从 1 开始，按 `&COORD` 或坐标文件的顺序）。支持区间写法 `1..16`，可写多行 `LIST` |
| `COMPONENTS_TO_FIX` | 固定哪些分量，默认 `XYZ`。例如只允许原子在 z 方向弛豫时写 `XY` |
| `MOLNAME` | 按分子名固定（配合 `&TOPOLOGY` 的分子定义，常用于 QM/MM） |
| `&FIXED_ATOMS` 可重复 | 不同原子组可以固定不同分量 |

> [!NOTE]
> 输入手册说明 `&FIXED_ATOMS` 约束的是**分数坐标**分量，而不是笛卡尔分量。对正交晶胞两者等价；对斜晶胞做 CELL_OPT 时，`COMPONENTS_TO_FIX X` 固定的是沿 $\mathbf a$ 方向的分数坐标。

### 5.2 实例：MgO(001) 表面上的水分子

示例 [../examples/07-geo-cell-opt/slab_geo_opt_fixed.inp](../examples/07-geo-cell-opt/slab_geo_opt_fixed.inp) 的模型：2×2 面内超胞（8.42 Å × 8.42 Å）、3 层 MgO（48 个原子），一个水分子以 O 端吸附在表层 Mg 上方约 2.2 Å，z 方向晶胞 25 Å（真空约 17 Å）。坐标文件 `mgo001_h2o.xyz` 中第 1–16 号原子是最底层，固定不动。核心片段：

```text
&DFT
  SURFACE_DIPOLE_CORRECTION T     ! molecule on one side only -> asymmetric slab
  SURF_DIP_DIR Z
  ...
  &XC
    &XC_FUNCTIONAL PBE
    &END XC_FUNCTIONAL
    &VDW_POTENTIAL                ! adsorption: include dispersion (chapter 8)
      DISPERSION_FUNCTIONAL PAIR_POTENTIAL
      &PAIR_POTENTIAL
        TYPE DFTD3(BJ)
        PARAMETER_FILE_NAME dftd3.dat
        REFERENCE_FUNCTIONAL PBE
      &END PAIR_POTENTIAL
    &END VDW_POTENTIAL
  &END XC
&END DFT
&SUBSYS
  &CELL
    ABC 8.42 8.42 25.0
    PERIODIC XYZ
  &END CELL
  ...
&END SUBSYS

&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 300
    &BFGS
      TRUST_RADIUS [angstrom] 0.20
    &END BFGS
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      COMPONENTS_TO_FIX XYZ
      LIST 1..16
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

设计要点：

- **面内晶格**必须取自同一套设置（泛函、基组、色散）下的体相 CELL_OPT 结果，本例为演示直接用实验值 4.21 Å；生产计算请先按第 6 节优化体相 MgO 再建 slab；
- **层数与固定层**：示例只有 3 层、固定 1 层，便于快速运行。生产计算通常至少 4 层，固定底部 1–2 层，并测试吸附能随层数的收敛；
- **偶极校正**：分子只吸附在一侧，slab 两侧不对称，应打开 `SURFACE_DIPOLE_CORRECTION`；
- **k 点**：MgO 是宽带隙绝缘体，8.4 Å 的面内超胞用 Γ 点 + OT 已基本够用；金属表面必须收敛 k 网格或用更大的超胞，并改用对角化和展宽（第5、6章）。

运行与检查：

```bash
mpirun -np 32 cp2k.psmp -i slab_geo_opt_fixed.inp -o slab_geo_opt_fixed.out
grep "OPT| Maximum gradient  " slab_geo_opt_fixed.out | tail -3
grep "GEOMETRY OPTIMIZATION COMPLETED" slab_geo_opt_fixed.out
```

验证底层确实没动（比较第一帧和最后一帧）：

```bash
python - <<'EOF'
import numpy as np
from ase.io import read
traj = read("mgo_h2o_opt-pos-1.xyz", index=":")
d = np.abs(traj[-1].positions[:16] - traj[0].positions[:16]).max()
print(f"max displacement of fixed atoms: {d:.2e} A")       # should be ~0
d2 = np.linalg.norm(traj[-1].positions[16:] - traj[0].positions[16:], axis=1).max()
print(f"max displacement of free atoms : {d2:.3f} A")
EOF
```

> [!CAUTION]
> 约束写错是“看似收敛、实则错误”的常见原因：序号从 0 开始数、坐标文件里底层原子不在最前面、换了建模软件后原子顺序被重排……每次改动结构文件后都要重新核对 `LIST`。另外，吸附计算中不要把整个衬底都固定，表层原子应允许弛豫。

---

## 6. 晶胞优化 CELL_OPT：体相 Si

### 6.1 输入文件

示例 [../examples/07-geo-cell-opt/si_cell_opt.inp](../examples/07-geo-cell-opt/si_cell_opt.inp)：8 原子立方 Si，初始晶格常数故意取 5.55 Å（实验值 5.431 Å）。

```text
@SET PROJ  si_cellopt
@SET ALAT  5.55

&GLOBAL
  PROJECT ${PROJ}
  RUN_TYPE CELL_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL
  &PRINT
    &STRESS_TENSOR ON
    &END STRESS_TENSOR
  &END PRINT
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 600
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
      EPS_SCF 1.0E-7
      MAX_SCF 100
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.3
        NBROYDEN 8
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC ${ALAT} ${ALAT} ${ALAT}
      SYMMETRY CUBIC
      &CELL_REF
        ABC 5.70 5.70 5.70
      &END CELL_REF
    &END CELL
    &COORD
      SCALED
      Si 0.00 0.00 0.00
      Si 0.00 0.50 0.50
      Si 0.50 0.00 0.50
      Si 0.50 0.50 0.00
      Si 0.25 0.25 0.25
      Si 0.25 0.75 0.75
      Si 0.75 0.25 0.75
      Si 0.75 0.75 0.25
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &CELL_OPT
    OPTIMIZER BFGS
    MAX_ITER 100
    EXTERNAL_PRESSURE [bar] 1.01325
    PRESSURE_TOLERANCE [bar] 100.0
    KEEP_SYMMETRY T
    MAX_FORCE 4.5E-4
    RMS_FORCE 3.0E-4
    MAX_DR 3.0E-3
    RMS_DR 1.5E-3
  &END CELL_OPT
  &PRINT
    &TRAJECTORY
      FORMAT XYZ
    &END TRAJECTORY
    &CELL ON
    &END CELL
    &STRESS ON
    &END STRESS
    &RESTART
      BACKUP_COPIES 1
    &END RESTART
  &END PRINT
&END MOTION
```

### 6.2 关键字说明

| 关键字 | 本例取值 | 说明 |
|---|---|---|
| `STRESS_TENSOR`（`&FORCE_EVAL` 下） | `ANALYTICAL` | 解析应力。可选 `NUMERICAL`（有限差分，非常慢）、`DIAGONAL_ANALYTICAL` 等。不写时 CELL_OPT 可能退回数值应力 |
| `&FORCE_EVAL/&PRINT/&STRESS_TENSOR` | `ON` | 在输出中打印以 STRESS 开头的应力表，默认只在 `PRINT_LEVEL HIGH` 时打印 |
| `CUTOFF` | 600 Ry | 比单点能更高，降低 Pulay 应力 |
| `&KPOINTS` | 4×4×4 | 8 原子 Si 立方胞需要 k 点；2026.2 及以前 k 点只能配对角化 |
| `SYMMETRY`（`&CELL` 下） | `CUBIC` | 声明初始晶胞的对称性，配合 `KEEP_SYMMETRY` 使用 |
| `&CELL_REF` | 5.70 Å | 参考晶胞，用来固定网格点数；应**大于**优化过程中可能出现的最大晶胞 |
| `EXTERNAL_PRESSURE` | 1.01325 bar | 外压（1 atm）。默认值是 100 bar，务必显式写出。可给 1 个数（静水压）或 9 个分量 |
| `PRESSURE_TOLERANCE` | 100 bar（默认） | 内压与外压之差的收敛阈值；精确晶格常数可收紧到 10–50 bar |
| `KEEP_SYMMETRY` | `T` | 保持 `&CELL/SYMMETRY` 声明的晶胞对称性（只作用于晶胞矢量） |
| `&MOTION/&PRINT/&CELL`、`&STRESS` | `ON` | 写 `<PROJECT>-1.cell` 和 `<PROJECT>-1.stress`，默认只在 HIGH 级别输出 |

> [!IMPORTANT]
> CP2K 约定晶胞矢量 $\mathbf a$ 沿 x 轴、$\mathbf b$ 在 xy 平面内（晶胞矩阵为下三角）。CELL_OPT 的手册明确说明，对上三角三个分量的更新没有完整的算法支持和测试，所以输入时这三个分量必须严格为 0。用 `ABC` + `ALPHA_BETA_GAMMA` 定义晶胞时 CP2K 会自动按此约定放置；直接写 `A`/`B`/`C` 矢量时需自己保证。2026.2 新增 `&CELL/CANONICALIZE`（`AUTO`/`TRUE`/`FALSE`），可以把一般取向的晶胞变换到这一约定。FCC 原胞等“非常规”取向的晶胞，请先做坐标变换再输入。

### 6.3 运行与输出

```bash
mpirun -np 16 cp2k.psmp -i si_cell_opt.inp -o si_cell_opt.out
grep "OPT| Internal pressure"  si_cell_opt.out
grep "OPT| Pressure deviation" si_cell_opt.out | tail -3
grep "CELL| Volume"           si_cell_opt.out | tail -1
grep "CELL| Vector a"         si_cell_opt.out | tail -1
grep -A5 "STRESS| Analytical stress tensor" si_cell_opt.out | tail -6
```

CELL_OPT 的 `OPT|` 块在 GEO_OPT 的基础上多了压力相关的几行（示意输出，数值仅供参考）：

```text
 OPT| Internal pressure [bar]                                  -22417.5838210771
 ...
 OPT| Pressure deviation [bar]                                  22418.5970710771
 OPT| Pressure tolerance [bar]                                    100.0000000000
 OPT| Pressure is converged                                                   NO
```

初始晶格偏大时内压为负（体系“想”收缩），随优化逐步趋近外压。`<PROJECT>-1.cell` 每行包含步数、时间、三个晶胞矢量的 9 个分量和体积，可以直接画体积–步数曲线；[`plot_opt.py`](../scripts/plot_opt.py) 也能读取（`python ../../scripts/plot_opt.py -h`）。

### 6.4 结果判断

- `GEOMETRY OPTIMIZATION COMPLETED` 出现，且最后一个 `OPT|` 块中压力、力、位移全部 `YES`；
- 最终晶格常数 = `CELL| Vector a` 一行中的 `|a|`。PBE 对 Si 的文献值一般在 5.46–5.47 Å 左右（近似范围，具体数值随基组和赝势变化；实验值 5.431 Å），偏离太多时检查 CUTOFF、k 点和基组；
- **用 E–V 曲线复核**：在最终晶格常数 ±2% 范围内取 5–7 个点做固定晶胞的单点计算，拟合状态方程。E–V 曲线的极小与 CELL_OPT 结果相差明显时，说明应力有系统误差（多半是 Pulay 应力），应提高 CUTOFF 或换基组。

> [!TIP]
> CELL_OPT 的应力对 CUTOFF 的收敛远慢于总能量。用 400 Ry 算单点能够用，不代表 400 Ry 下的晶格常数可靠。做法是：固定结构，在 400/500/600/800 Ry 下比较 `STRESS| 1/3 Trace`，变化小于目标精度（如 0.1 GPa）后再开始 CELL_OPT。批量扫描可用 `python ../../scripts/conv_cutoff.py -h`。

---

## 7. 控制晶胞自由度：KEEP_* 与 CONSTRAINT

### 7.1 关键字一览

| 关键字（`&MOTION/&CELL_OPT`） | 作用 | 版本 / 依赖 | 典型用途 |
|---|---|---|---|
| `KEEP_ANGLES T` | 晶胞夹角不变，三个边长独立变化 | 各版本 | 六方/层状材料只优化 a 和 c；三斜晶胞保持角度 |
| `KEEP_SYMMETRY T` | 保持 `&CELL/SYMMETRY` 声明的晶胞对称性 | 各版本 | 立方、四方、六方晶体 |
| `KEEP_SPACE_GROUP T` | 用 spglib 识别空间群并在优化中保持（作用于坐标、力、应力） | 需编译链接 spglib（`cp2k.psmp -v` 的 cp2kflags 中有 `spglib`） | 防止高对称结构因数值噪声破缺对称性 |
| `KEEP_VOLUME T` | 体积不变，只优化晶胞形状 | **2026.2 新增** | 定体积下的 c/a、应变研究 |
| `CONSTRAINT` | 固定部分晶胞分量 | 写法随版本变化，见 7.3 | 二维材料、slab 只优化面内 |
| `EPS_SYMMETRY` | 空间群识别容差（分数坐标） | 默认 1.0E-4 | 配合 `KEEP_SPACE_GROUP`、`SHOW_SPACE_GROUP` |
| `SHOW_SPACE_GROUP T` | 优化结束后打印识别到的空间群 | 需 spglib | 检查对称性是否被破坏 |

`KEEP_VOLUME` 可以与 `KEEP_ANGLES` 或 `KEEP_SYMMETRY` 组合使用。`KEEP_SPACE_GROUP` 通常与 `KEEP_SYMMETRY` 一起用。固定晶胞的 GEO_OPT 也可以用 `&GEO_OPT/KEEP_SPACE_GROUP` 保持原子坐标的空间群对称性。

> [!WARNING]
> 对称性约束只在你**确知**极小结构保持该对称性时使用，否则会把结构锁在高对称的鞍点上（例如 BaTiO3 的立方相、有 Jahn–Teller 畸变的体系）。不确定时，先不加约束优化，再用 `SHOW_SPACE_GROUP T` 或 spglib 检查结果的对称性。

### 7.2 层状材料：2H-MoS2 + D3(BJ)

示例 [../examples/07-geo-cell-opt/layered_cell_opt_keep_angles.inp](../examples/07-geo-cell-opt/layered_cell_opt_keep_angles.inp) 的初始 c = 13.00 Å（实验值约 12.29 Å），用 `KEEP_ANGLES T` 保持 90°/90°/120°，a 和 c 独立弛豫：

```text
&CELL
  ABC ${A0} ${A0} ${C0}
  ALPHA_BETA_GAMMA 90.0 90.0 120.0
&END CELL
...
&XC
  &XC_FUNCTIONAL PBE
  &END XC_FUNCTIONAL
  &VDW_POTENTIAL
    DISPERSION_FUNCTIONAL PAIR_POTENTIAL
    &PAIR_POTENTIAL
      TYPE DFTD3(BJ)
      PARAMETER_FILE_NAME dftd3.dat
      REFERENCE_FUNCTIONAL PBE
    &END PAIR_POTENTIAL
  &END VDW_POTENTIAL
&END XC
...
&CELL_OPT
  OPTIMIZER BFGS
  MAX_ITER 150
  EXTERNAL_PRESSURE [bar] 1.01325
  PRESSURE_TOLERANCE [bar] 50.0
  KEEP_ANGLES T
&END CELL_OPT
```

为什么一定要加色散校正：层间结合主要来自范德华作用，纯 PBE 几乎不给层间结合，c 会弛豫到明显偏大的值，甚至层间“解离”。把 `&VDW_POTENTIAL` 段删掉重跑一次，对比两次的最终 c，就能直观看到差别。

> [!TIP]
> 六方晶体也可以用 `&CELL/SYMMETRY HEXAGONAL` 配 `KEEP_SYMMETRY T`，效果是 a = b 且角度不变。`KEEP_ANGLES` 只保角，不保证 a = b，适合对称性较低的层状结构。

### 7.3 二维材料与 slab：CONSTRAINT

对真空方向不能变的体系（单层 MoS2、石墨烯、slab），要么直接用 GEO_OPT（面内晶格取体相值），要么在 CELL_OPT 中固定 c 方向：

```text
&CELL_OPT
  CONSTRAINT Z          ! 2026.2: c vector fixed (only c has a z component)
  KEEP_ANGLES T
  EXTERNAL_PRESSURE [bar] 1.01325
&END CELL_OPT
```

`CONSTRAINT` 的取值随版本有差别：

| 版本 | 可选值 | 含义 |
|---|---|---|
| ≤ 2026.2 | `NONE`、`X`、`Y`、`Z`、`XY`、`XZ`、`YZ` | 固定压力张量/晶胞的对应分量（如 `Z` 固定 z 分量） |
| trunk（2026.2 之后的开发版手册） | 在上述基础上增加 `A`、`B`、`C`、`AX`、`BX`、`BY`、`CX`、`CY`、`CZ`、`ALPHA`、`BETA`、`GAMMA`，一行可写多个 | 更细的分量控制；`Z` 等价于固定整个 c 矢量 |

只用 `Z`、`XY` 这类旧写法可以在各版本间通用。

> [!NOTE]
> 二维材料面内晶格常数的优化对真空层中的“假压力”很敏感：真空方向的应力分量不为零也不应被优化，所以必须约束。另外，官方文档提醒：外压只作用在被优化的方向上，slab 和含真空的体系除非有物理理由，否则不要让真空方向参与弛豫。

---

## 8. CP2K 2026.2 的变化与旧输入迁移

### 8.1 删除 `&CELL_OPT/TYPE`（破坏性变更）

2025.x 及以前，`&CELL_OPT` 有一个 `TYPE` 关键字：

| 旧 `TYPE` 值 | 含义（2025.2 手册） |
|---|---|
| `DIRECT_CELL_OPT`（默认） | 每一步同时优化原子和晶胞，每步计算应力 |
| `GEO_OPT` | 每个晶胞步之间先做一次完整的几何优化，再在优化好的结构上算应力 |
| `MD` | 用一段 MD 计算应力 |

2026.2（PR #5257）删除了这个关键字，晶胞优化固定使用直接法。旧输入中的 `TYPE DIRECT_CELL_OPT` 在 2026.2 下会直接中止。根据 2026.2 源码（输入解析器对“已删除关键字”的处理），报错信息形如：

```text
The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE
has been removed because cell optimizations now always use DIRECT_CELL_OPT.
```

迁移方法：**删掉这一行即可**。因为旧版本的默认值本来就是 `DIRECT_CELL_OPT`，删掉后的输入在 2023.x–2026.2 都能运行，本章示例都按这种方式书写。原来依赖 `TYPE GEO_OPT` 的工作流，可以改成“GEO_OPT ↔ 固定晶胞 E–V 扫描”或直接用 `DIRECT_CELL_OPT` 并收紧判据。

### 8.2 新功能速查

| 功能 | 关键字 / 输出 | 版本 |
|---|---|---|
| 定体积晶胞优化 | `&CELL_OPT/KEEP_VOLUME T` | 2026.2 新增 |
| 优化结束输出最终结构 | `&MOTION/&PRINT/&FINAL_STRUCTURE`，写 `<PROJECT>-FINAL-1_<步数>.cif` 和同名 `.xyz`（扩展 XYZ，注释行含 `Lattice=`） | 2026.2 新增，默认打印级别 LOW，即默认开启 |
| 晶胞规范化 | `&CELL/CANONICALIZE AUTO`（默认）/`TRUE`/`FALSE` | 2026.2 新增 |
| GExt 波函数外推 | `&QS/EXTRAPOLATION GEXT_PROJ`、`GEXT_PROJ_QTR` | 2026.2 新增 |
| k 点计算的波函数外推 | 与 Γ 点相同的 `EXTRAPOLATION` 设置 | 2026.2 新增 |
| LBFGS 打印控制 | `&LBFGS/PRINT_LEVEL` | 2026.2 新增 |

`FINAL_STRUCTURE` 的 CIF 一律以 P1 空间群写出，分数坐标落在 [0,1)；`PRINT_CIF`、`PRINT_XYZ` 可分别关闭；`PRINT_ATOM_KIND T` 会在扩展 XYZ 中写 KIND 名（如 `Ni1`）而不是元素符号。显式写法：

```text
&MOTION
  &PRINT
    &FINAL_STRUCTURE ON
      PRINT_CIF T
      PRINT_XYZ T
    &END FINAL_STRUCTURE
  &END PRINT
&END MOTION
```

> [!CAUTION]
> 上面的 `&FINAL_STRUCTURE` 和 `KEEP_VOLUME`、`CANONICALIZE`、`GEXT_PROJ` 在 2026.1 及更早版本中不存在，写进去会报 unknown keyword/subsection。需要跨版本使用的输入，请把这些行放进 `@IF` 块或另存一个 2026.2 专用版本，运行前用 `cp2k.psmp -c x.inp` 检查。

---

## 9. 加速优化中的 SCF：外推与重启

每个优化步的结构变化很小，用上一步的电子结构做初猜可以把 SCF 步数从几十降到几步。`&QS/EXTRAPOLATION` 的常用选项：

| 选项 | 说明 | 建议 |
|---|---|---|
| `ASPC`（默认） | Always Stable Predictor-Corrector，稳定 | 通用首选，`EXTRAPOLATION_ORDER 3` |
| `PS` | 密度矩阵×重叠矩阵的高阶外推，初猜精度高 | GEO_OPT 中常比 ASPC 少几步 SCF |
| `GEXT_PROJ` | GExt 外推（2026.2 新增），官方文档称通常是很好的首选 | 2026.2 起可试用 |
| `USE_PREV_WF` | 直接用上一步波函数 | 保守的后备选项 |
| `USE_GUESS` | 不外推，每步按 `SCF_GUESS` 重新初猜 | 结构变化很大或怀疑外推出错时 |
| `USE_PREV_P`、`LINEAR_P` | 基于密度矩阵 | CELL_OPT 中慎用，见下 |

> [!WARNING]
> 晶胞或近邻列表变化时，基于密度矩阵的外推（`USE_PREV_P`、`LINEAR_P`）可能给出很差的初猜，CP2K 会警告 `Change in cell neighborlist: might affect quality of initial guess`。如果同时出现能量、力或 SCF 行为异常，改用 `ASPC` 或 `USE_PREV_WF`。另外，k 点对称约化（2026.2 实验性功能）下，若对称性随结构更新而改变，外推可能自动退回 `USE_GUESS`。

其它要点：

- 不要在优化中使用 `IGNORE_CONVERGENCE_FAILURE`。SCF 不收敛的那一步力是错的，下一步的外推也会被污染，结构“炸开”并不奇怪；
- 从一次便宜计算（低 CUTOFF 单点、或 PBE 预优化）的 `<PROJECT>-RESTART.wfn` 起步，可以省掉第一步的大量 SCF；
- 2026.2 及以前，k 点计算只能用对角化（OT 在 k 点支持表中标为 Unsupported），金属体系配合 `&SMEAR` 和足够的 `ADDED_MOS`，见 [第5章](05-scf-strategy.md)。

---

## 10. 输出文件、监控与续算

### 10.1 输出文件

| 文件 | 内容 | 控制 |
|---|---|---|
| `<PROJECT>-pos-1.xyz` | 每一步的结构（轨迹） | `&MOTION/&PRINT/&TRAJECTORY`，默认 XMOL/XYZ 格式，不含晶胞 |
| `<PROJECT>-1.restart` | 最新结构、晶胞和完整输入，可直接作为输入续算 | `&MOTION/&PRINT/&RESTART` |
| `<PROJECT>-1.restart.bak-1` … | 前几步的备份 | `BACKUP_COPIES` |
| `<PROJECT>-RESTART.wfn` | 最新波函数（Γ 点）；k 点计算为 `<PROJECT>-RESTART.kp` | `&DFT/&SCF/&PRINT/&RESTART` |
| `<PROJECT>-BFGS.Hessian` | BFGS 近似 Hessian | BFGS 自动写出 |
| `<PROJECT>-1.cell` | 每步晶胞矢量与体积 | `&MOTION/&PRINT/&CELL`（默认 HIGH，需手动打开） |
| `<PROJECT>-1.stress` | 每步应力张量 | `&MOTION/&PRINT/&STRESS`（默认 HIGH，需手动打开） |
| `<PROJECT>-FINAL-1_<N>.cif/.xyz` | 最终结构 | 2026.2 新增 `&FINAL_STRUCTURE` |

> [!NOTE]
> 普通 XYZ 不带晶胞信息。CELL_OPT 的轨迹想在 Ovito/ASE 中带晶胞显示，可以把 `&TRAJECTORY` 的 `FORMAT` 改成 `EXTXYZ`（2026.2 输入手册中已列出该选项）或 `PDB`，或者结合 `-1.cell` 文件自己组装。

### 10.2 监控长时间优化

```bash
# live view of the convergence block
tail -f si_cell_opt.out | grep --line-buffered -E "OPT\| (Step number|Maximum gradient |Pressure deviation)"
# how many SCF iterations per optimization step (should drop to a few after the first steps)
grep -c "SCF run converged" si_cell_opt.out
```

除了数字，也要定期下载轨迹在可视化软件里看一眼：原子是否成键异常、晶胞是否畸变。发现问题及时 `touch EXIT_GEO` 或 `scancel`。

### 10.3 续算

方法一：直接把 `-1.restart` 当输入（它是一个完整的 CP2K 输入，含最新坐标、晶胞和计步器）：

```bash
cp si_cellopt-1.restart si_cellopt_cont.inp
mpirun -np 16 cp2k.psmp -i si_cellopt_cont.inp -o si_cellopt_cont.out
```

方法二：保留原输入，加 `&EXT_RESTART`：

```text
&EXT_RESTART
  RESTART_FILE_NAME si_cellopt-1.restart
  ! RESTART_DEFAULT T (default) restores positions, cell, counters ...
&END EXT_RESTART
```

同时在 `&DFT` 中设置 `WFN_RESTART_FILE_NAME si_cellopt-RESTART.kp`（Γ 点计算为 `.wfn`）和 `&SCF/SCF_GUESS RESTART`，省掉第一步的 SCF。BFGS 还可以加 `RESTART_HESSIAN T` 读回 Hessian。

> [!WARNING]
> Γ 点的 `.wfn` 和 k 点的 `.kp` 不能互换；换了原子数、元素顺序、基组或赝势后，波函数重启会失败或给出错误初猜。另外，续算时若沿用同一个 `PROJECT`，新的 `-1.restart` 会覆盖旧文件，重要节点请先备份。

把最终坐标写回新输入或 xyz（例如为后续的能带、振动计算准备结构），可以用 `python ../../scripts/restart_update.py -h` 查看用法；2026.2 也可以直接使用 `FINAL` 开头的 CIF/扩展 XYZ 文件。

---

## 参数选择与收敛建议

| 参数 | 一般 GEO_OPT | 后接振动分析 / 过渡态验证 | CELL_OPT | 说明 |
|---|---|---|---|---|
| `MAX_FORCE`（Ha/bohr） | 4.5E-4（默认） | 1.0E-4 或更小 | 4.5E-4 → 精修时 1.0E-4 | 1E-4 Ha/bohr ≈ 0.005 eV/Å |
| `RMS_FORCE` | 3.0E-4 | 约 `MAX_FORCE` 的 2/3 | 同左 | |
| `MAX_DR` / `RMS_DR`（bohr） | 3.0E-3 / 1.5E-3 | 1.0E-3 / 5.0E-4 | 默认 | 松软体系位移判据可能很难满足，可适当放宽 |
| `EPS_SCF` | 1E-6 ~ 1E-7 | 1E-7 ~ 1E-8 | 1E-7 | 比单点能更紧 |
| `EPS_DEFAULT` | 1E-12 | 1E-12 ~ 1E-14 | 1E-12 | |
| `CUTOFF`（Ry） | 单点收敛值 | 单点收敛值 + 100 | ≥ 500–600，用应力判收敛 | 配 `CELL_REF` |
| `PRESSURE_TOLERANCE`（bar） | — | — | 100（默认），精修 10–50 | 1 GPa = 10 kbar = 10000 bar |
| `TRUST_RADIUS`（Å） | 0.25（默认） | 0.1–0.2 | 0.1–0.25 | 振荡时减小 |
| 优化器 | BFGS | BFGS | BFGS；困难时 CG | 大体系用 LBFGS |
| `EXTRAPOLATION` | ASPC（默认）/ PS | ASPC | ASPC / USE_PREV_WF | 2026.2 可试 GEXT_PROJ |
| `MAX_ITER` | 200（默认） | 200–400 | 100–200 | 柔性大分子可放宽到 300–400 |

---

## 常见问题与报错

**1. `The specified keyword 'TYPE' is not available anymore: The keyword MOTION/CELL_OPT/TYPE has been removed ...`**
原因：在 2026.2 中使用了旧输入的 `TYPE DIRECT_CELL_OPT`（或 Multiwfn、旧脚本生成的输入）。
解决：删除 `&CELL_OPT` 中的 `TYPE` 行。删掉后的输入在旧版本中行为不变（默认就是直接法）。

**2. `*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***`**
原因：`MAX_ITER` 用完仍未收敛。
解决：先看 `OPT|` 块中哪个判据卡住。若力已很小、只有位移判据不满足（松软体系常见），可适当放宽 `MAX_DR`/`RMS_DR`；若力在一个范围内来回振荡，减小 `TRUST_RADIUS`、收紧 `EPS_SCF`，或换 CG；然后从 `-1.restart` 续算，而不是从头开始。

**3. 能量/力在几个值之间振荡，始终达不到 `MAX_FORCE`**
原因：力的数值噪声大于收敛阈值，常见来源有 SCF 不够紧（`EPS_SCF` 1E-5）、`EPS_DEFAULT` 太松、CUTOFF 太低导致的网格“蛋箱效应”。
解决：`EPS_SCF` 收紧到 1E-7，`EPS_DEFAULT 1E-12`，提高 CUTOFF；CG 的线搜索换成 `FIT`。

**4. `[ABORT] SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`（出现在优化中途）**
原因：自 2024.1 起 SCF 不收敛默认中止。优化中途出现，往往是结构走到了不合理的位置（原子太近、键断裂）或外推给出了很差的初猜。
解决：查看最后几帧结构；把 `EXTRAPOLATION` 改为 `USE_PREV_WF` 或 `USE_GUESS`；减小 `TRUST_RADIUS`；按[第5章](05-scf-strategy.md)增强 SCF（`OUTER_SCF`、换预条件子、金属用对角化 + 展宽）。不要用 `IGNORE_CONVERGENCE_FAILURE` 硬跑优化。

**5. LBFGS 提示 `Specific L-BFGS convergence criteria ... run CONVERGED!`，但最大力仍大于 `MAX_FORCE`**
原因：LBFGS 的专用判据先满足了。
解决：从 `-1.restart` 改用 BFGS（或 CG）续算到一般判据满足，或收紧 `WANTED_PROJ_GRADIENT`/`WANTED_REL_F_ERROR`。

**6. CELL_OPT 得到的晶格常数明显偏小（或偏大），且随 CUTOFF 变化**
原因：Pulay 应力；或者忘了写 `STRESS_TENSOR ANALYTICAL`，走了数值应力；或者 `EXTERNAL_PRESSURE` 用了默认的 100 bar。
解决：提高 CUTOFF 并以应力为判据做收敛，设置 `CELL_REF`，检查输出里 `STRESS| Analytical stress tensor` 是否出现，显式写 `EXTERNAL_PRESSURE [bar] 1.01325`。

**7. `*** WARNING ... Change in cell neighborlist: might affect quality of initial guess ***`**
原因：CELL_OPT 中晶胞变化，基于密度矩阵的外推初猜质量下降。
解决：若 SCF 和能量正常可以忽略；若伴随异常，改用 `ASPC` 或 `USE_PREV_WF`。

**8. `Symmetry library SPGLIB not available`**
原因：用了 `KEEP_SPACE_GROUP` 或 `SHOW_SPACE_GROUP`，但 CP2K 编译时没有链接 spglib。
解决：去掉这些关键字改用 `KEEP_SYMMETRY`，或重新编译带 spglib 的版本（toolchain 中 `--with-spglib`，或用 conda-forge/容器版本，先用 `cp2k.psmp -v` 确认 cp2kflags 中有 `spglib`）。

**9. 二维材料/slab 在 CELL_OPT 中真空层被压缩或拉长**
原因：没有约束真空方向。
解决：`CONSTRAINT Z`（c 沿 z 时），或干脆用 GEO_OPT、面内晶格取体相值。

**10. 晶胞矢量的上三角分量不为零导致 CELL_OPT 行为异常**
原因：直接用 `A`/`B`/`C` 写了“非常规”取向的晶胞（如 FCC 原胞的标准写法）。
解决：改用 `ABC` + `ALPHA_BETA_GAMMA`，或先把晶胞和坐标一起旋转到 a∥x、b 在 xy 面内的取向（ASE 中 `atoms.cell.standard_form()` 可以给出这种下三角形式）；2026.2 也可用 `&CELL/CANONICALIZE TRUE`。

**11. `GEOMETRY wrong or EMAX_SPLINE too small`**
原因：原子重叠，常见于建模时周期边界处产生的重复原子，或 xyz 文件缺少晶胞信息导致原子被错误折回。
解决：用 VESTA/Ovito 显示周期像检查结构，删除重复原子，确认 `&CELL` 与坐标匹配。

---

## 练习

1. 用示例中的 `h2o_geo_opt_bfgs.inp` 分别以 BFGS、LBFGS、CG（`2PNT` 与 `GOLD`）运行，统计每种方法的优化步数和总 SCF 迭代次数（`grep -c "OPT| Step number"`、`grep "outer SCF loop converged" | wc -l` 或 [`cp2k_outparse.py`](../scripts/cp2k_outparse.py)），比较哪种最省。
2. 把 `MAX_FORCE` 收紧到 1.0E-4 继续优化（从 `-1.restart` 续算），观察键长和键角还变化多少；思考为什么振动分析之前需要这一步。
3. 对 Si 分别用 CUTOFF 300、400、600、800 Ry 运行 `si_cell_opt.inp`，画出晶格常数随 CUTOFF 的变化；再去掉 `&CELL_REF` 重复一次，比较曲线的平滑程度。
4. 在 Si 的最终晶格常数附近做 7 个点的 E–V 扫描（固定晶胞 `RUN_TYPE ENERGY`），用 ASE 的 `EquationOfState` 拟合，比较 E–V 极小与 CELL_OPT 结果之差。
5. 去掉 MoS2 示例中的 `&VDW_POTENTIAL` 重跑，比较 c 轴长度；再把 `KEEP_ANGLES T` 换成 `&CELL/SYMMETRY HEXAGONAL` + `KEEP_SYMMETRY T`，检查结果是否一致。
6. （2026.2）在 Si 示例中加入 `KEEP_VOLUME T`，同时把初始晶胞改成略微四方畸变（如 `ABC 5.40 5.40 5.62`、去掉 `SYMMETRY CUBIC` 和 `KEEP_SYMMETRY`），观察体积不变时晶胞形状如何回到立方。

---

## 参考资料

**官方文档**

- Geometry and Cell Optimization（方法文档，含输出说明与检查清单）：<https://manual.cp2k.org/trunk/methods/optimization/geometry_and_cell_opt.html>
- `MOTION/GEO_OPT` 输入参考：<https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/GEO_OPT.html>
- `MOTION/CELL_OPT` 输入参考（2026.2 版本，无 `TYPE`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/CELL_OPT.html>
- `MOTION/CELL_OPT` 输入参考（2025.2 版本，含 `TYPE`）：<https://manual.cp2k.org/cp2k-2025_2-branch/CP2K_INPUT/MOTION/CELL_OPT.html>
- `MOTION/CONSTRAINT/FIXED_ATOMS`：<https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/CONSTRAINT/FIXED_ATOMS.html>
- `MOTION/PRINT/FINAL_STRUCTURE`（2026.2 新增）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/MOTION/PRINT/FINAL_STRUCTURE.html>
- Restarting CP2K Calculations：<https://manual.cp2k.org/trunk/methods/restarting.html>
- How to make a SCF run converge：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- CP2K v2026.2 release notes（`&CELL_OPT/TYPE` 删除、定体积优化、CIF/EXTXYZ 输出）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- 旧版 HowTo（geometry_optimisation 等）：<https://www.cp2k.org/howto>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：官方源码；`tests/` 下的 regtest 是查 GEO_OPT/CELL_OPT 真实写法的权威来源
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库，含 `cp2k-made-simple` 配套输入
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含几何优化、slab 与约束的实践经验
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程，含 GEO_OPT/CELL_OPT 章节（2025.2 验证，注意 2026.2 的 `TYPE` 变化）
- [EPCCed/intro_to_cp2k](https://github.com/EPCCed/intro_to_cp2k)：ARCHER2 入门课，含几何优化练习
- [ishikawa-group/cp2k-tutorial](https://github.com/ishikawa-group/cp2k-tutorial)：日文入门教程，Pt 表面加苯分子的优化例子可作 slab 参考
- [ybx20021218/CP2Kview--Visualization-of-CP2K-Output-Files](https://github.com/ybx20021218/CP2Kview--Visualization-of-CP2K-Output-Files)：在浏览器里查看优化轨迹和收敛曲线
- [Doria58/DFT_Script](https://github.com/Doria58/DFT_Script)：含把 CP2K 优化结构转成 VASP POSCAR 的脚本
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：OT/DIAG/GEO/CELL/MD 一键切换、restart 回写的思路（无 license，仅供参考）
- [RobinLi-95/CP2K_tools](https://github.com/RobinLi-95/CP2K_tools)：含 EOS 拟合脚本，可用于 E–V 复核

**本仓库脚本**：`../scripts/plot_opt.py`（优化收敛图）、`../scripts/cp2k_outparse.py`（解析 OPT| 块）、`../scripts/restart_update.py`（从 restart 回写坐标/晶胞）、`../scripts/job_cp2k_slurm.sh`（作业模板）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- R. H. Byrd, P. Lu, J. Nocedal, C. Zhu, *SIAM J. Sci. Comput.* **16**, 1190 (1995)，L-BFGS-B：<https://doi.org/10.1137/0916069>
- R. Lindh, A. Bernhardsson, G. Karlström, P.-Å. Malmqvist, *Chem. Phys. Lett.* **241**, 423 (1995)，模型 Hessian：<https://doi.org/10.1016/0009-2614(95)00646-L>
- J. Kolafa, *J. Comput. Chem.* **25**, 335 (2004)，ASPC 外推：<https://doi.org/10.1002/jcc.10385>
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999)，表面偶极校正：<https://doi.org/10.1103/PhysRevB.59.12301>
- A. Togo, I. Tanaka, Spglib: a software library for crystal symmetry search, arXiv:1808.01590 (2018)：<https://arxiv.org/abs/1808.01590>
- T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020)，CP2K 综述：<https://doi.org/10.1063/5.0007045>
- M. Iannuzzi *et al.*, The CP2K Program Package Made Simple, *J. Phys. Chem. B* **130**, 1237 (2026)：<https://doi.org/10.1021/acs.jpcb.5c05851>

---

[← 上一章](06-periodic-surface-kpoints.md) | [返回目录](../README.md) | [下一章 →](08-dispersion-spin-dftu.md)
