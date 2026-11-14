# 第6章 结构优化、晶胞优化、状态方程与 NEB（Geometry & Cell Optimization, Equation of State and NEB with OpenMX）

本章讲 OpenMX 中与"找结构"有关的全部功能：原子位置优化（最速下降 Opt、DIIS、BFGS、RF、EF）及约束 `MD.Fixed.XYZ`、用 `*.dat#` 续算、晶胞优化（OptC1–OptC7、RFC5–RFC7，推荐 RFC5 + `MD.Fixed.Cell.Vectors`）与焓优化、`MD.Type EvsLC` 状态方程扫描与 Birch–Murnaghan 拟合、DFT-D2/D3 范德华修正，以及用 NEB（nudged elastic band）求最小能量路径和扩散势垒。每一部分都配有可直接运行的输入（OpenMX geometry optimization, variable-cell relaxation, E–V curve, DFT-D3, NEB barrier）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 选择合适的 `MD.Type` 和收敛判据，读懂 `.out` 中的优化历史，判断优化是否真正收敛；
- 用 `MD.Fixed.XYZ` 固定原子坐标，用 `*.dat#` 文件续算，处理收敛到鞍点或长时间不收敛的情况；
- 用 RFC5 做晶胞优化，用 `MD.Fixed.Cell.Vectors` 约束晶格分量，用 `MD.applied.pressure` 做加压下的焓优化；
- 用 `MD.Type EvsLC` 扫描晶格常数，固定 `scf.Ngrid` 得到光滑的 E–V 曲线，并拟合 V0、B0、B0'；
- 正确开启 DFT-D3（BJ 阻尼），处理层状材料和表面的 vdW 方向设置；
- 完成一次 NEB 计算：端点优化、构造路径、并行设置、监控各 image 的 SCF、从 `.neb.ene` 读出势垒。

| 任务 | `MD.Type` | 本章小节 | 示例目录 |
|---|---|---|---|
| 分子、表面、缺陷的原子弛豫 | `RF`（首选）/ `EF` / `BFGS` | §3 | [06-geomopt-h2o-ch4](../examples/06-geomopt-h2o-ch4/README.md) |
| 初始结构很差，先粗调 | `Opt`（最速下降） | §3.4 | 同上 |
| 晶格常数和原子同时优化 | `RFC5` | §4 | [06-cellopt-eos-si-c](../examples/06-cellopt-eos-si-c/README.md) |
| 体模量、E–V 曲线 | `EvsLC` | §5 | 同上 |
| 层状、分子晶体、物理吸附 | 上述类型 + `scf.dftD on` | §6 | [06-vdw-d3-graphite](../examples/06-vdw-d3-graphite/README.md) |
| 反应路径、扩散势垒 | `NEB` | §7 | [06-neb-h-in-si](../examples/06-neb-h-in-si/README.md) |
| 分子动力学 | `NVE`、`NVT_*`、`NPT_*` | — | [第7章](07-md.md) |

> [!NOTE]
> 结构优化用的截断能、k 点和基组应来自 [第4章](04-pao-vps.md)、[第5章](05-scf-kpoints-convergence.md) 的收敛测试。力和应力对这些参数比总能更敏感，必要时用比单点计算更严格的设置。

---

## 2. 原理速览

### 2.1 力、应力与单位

每个优化步先做 SCF，再计算每个原子受的力 $\mathbf F_I=-\partial E/\partial\mathbf R_I$。PAO 基组以原子为中心、随原子移动，因此力中除了 Hellmann–Feynman 项还包含基组移动带来的项，OpenMX 解析地计算全部贡献。晶胞优化还需要应力张量：手册说明 OpenMX 根据总能对应变的梯度解析计算应力，并把它变换成对晶格矢量的梯度（这种处理中基函数中心随应变移动、形状保持不变）。

OpenMX 中力的单位是 Hartree/Bohr。常用换算：

| Hartree/Bohr | eV/Å | 含义 |
|---|---|---|
| 1 | 51.422 | — |
| 3.0×10⁻⁴ | 0.0154 | `MD.Opt.criterion` 默认值 |
| 1.0×10⁻⁴ | 0.0051 | 常用的较严格判据 |
| 5.0×10⁻⁴ | 0.0257 | 手册对弱相互作用体系建议的折中判据 |

`MD.Opt.criterion` 比较的是**所有原子受力分量绝对值的最大值**（晶胞优化时还包括对晶格矢量的梯度）。

### 2.2 优化算法

| 方法 | 思路 | 特点 |
|---|---|---|
| `Opt` | 最速下降，步长因子随收敛行为自动调节 | 稳健但慢，适合离极小点很远时粗调 |
| `DIIS` | 用历史几何和梯度外推，近似 Hessian 取单位矩阵 | 便宜，接近极小点时有效 |
| `BFGS` | 逐步更新近似逆 Hessian 的经典拟牛顿法 | 与 EF 表现接近 |
| `RF` | 有理函数法，用增广 Hessian 稳定每一步，Hessian 用 BFGS 更新 | 手册基准测试中多数情况下**最稳健、最高效** |
| `EF` | 本征矢跟踪，监控并修正近似 Hessian 的本征值 | 稳健，Hessian 不合适时也较稳定 |

拟牛顿方法在开始前先跑 `MD.Opt.StartDIIS` 步最速下降。若在越过拐点之前就切换到拟牛顿法，结构可能收敛到**鞍点**而非极小点，这时应增大 `MD.Opt.StartDIIS`。

### 2.3 状态方程

三阶 Birch–Murnaghan 方程：

$$
E(V)=E_0+\frac{9V_0B_0}{16}\left\{\left[\left(\tfrac{V_0}{V}\right)^{2/3}-1\right]^3B_0'+\left[\left(\tfrac{V_0}{V}\right)^{2/3}-1\right]^2\left[6-4\left(\tfrac{V_0}{V}\right)^{2/3}\right]\right\}
$$

它等价于 $E$ 对 $x=V^{-2/3}$ 的三次多项式，所以可以用线性最小二乘拟合，再由多项式的导数求 $V_0$、$B_0=V\,\partial^2E/\partial V^2|_{V_0}$ 和 $B_0'$。本仓库示例中的 `fit_evslc.py` 就是这样做的。

### 2.4 NEB

NEB 在反应物和产物之间放置若干 image，用虚拟弹簧连成"弹性带"。每个 image 受到的力为

$$
\mathbf F_i^{\text{NEB}}=\mathbf F_i^{\perp}+\mathbf F_i^{s,\parallel},\qquad
\mathbf F_i^{\perp}=-\nabla V(\mathbf R_i)+\left[\nabla V(\mathbf R_i)\cdot\hat\tau_i\right]\hat\tau_i,
$$

$$
\mathbf F_i^{s,\parallel}=k\left(|\mathbf R_{i+1}-\mathbf R_i|-|\mathbf R_i-\mathbf R_{i-1}|\right)\hat\tau_i .
$$

真实力只保留垂直于路径的分量，使 image 落到最小能量路径（MEP）上；弹簧力只保留平行分量，使 image 沿路径均匀分布。OpenMX 4.0 用 Henkelman–Jónsson 改进切线定义 $\hat\tau_i$，并用 DIIS + BFGS 混合方法优化整条路径，支持 MPI/OpenMP 混合并行。

---

## 3. 原子位置优化

### 3.1 关键词

| 关键词 | 默认值 | 推荐 | 说明 |
|---|---|---|---|
| `MD.Type` | `Nomd` | `RF` 或 `EF` | `Opt`、`DIIS`、`BFGS`、`RF`、`EF` |
| `MD.maxIter` | **1** | 100–300 | 最大优化步数；**默认只走 1 步**，忘写是最常见的错误 |
| `MD.Opt.criterion` | 3.0×10⁻⁴ Hartree/Bohr | 1×10⁻⁴（体相、要算声子时更严） | 最大力的收敛阈值 |
| `MD.Opt.StartDIIS` | 5 | 5–20 | 前几步用最速下降，再切换到 DIIS/BFGS/RF/EF |
| `MD.Opt.DIIS.History` | 3 | 3–5 | 拟牛顿法更新 Hessian 时参考的历史步数 |
| `MD.Opt.EveryDIIS` | 200 | 保持默认 | 手册未展开说明，按默认即可 |
| `MD.Opt.Init.Hessian` | `Schlegel` | `Schlegel` | 初始 Hessian：`Schlegel`（基于键伸缩的模型 Hessian）或 `iden`；**`Schlegel` 不能用于 `BFGS` 和 `DIIS`** |
| `MD.Fixed.XYZ` | 全部放开 | — | 固定坐标分量，见 §3.5 |

> [!NOTE]
> 手册的关键词总表把最速下降写作 `SD`，但第18章正文和官方示例 `Methane2.dat` 使用的都是 `MD.Type Opt`，本教程也用 `Opt`。

### 3.2 实战：用 RF 优化水分子

[`examples/06-geomopt-h2o-ch4/H2O_opt_RF.dat`](../examples/06-geomopt-h2o-ch4/README.md) 的完整内容如下，初始结构故意扭曲（r(OH) = 1.012 Å，∠HOH = 114.2°）：

```text
#
# H2O molecule: geometry optimization with the rational-function (RF) method
# Start geometry is distorted on purpose: r(OH) = 1.012 Ang, angle(HOH) = 114.2 deg
# OpenMX 4.0 (+patch 4.0.1), database DFT_DATA19
#
System.CurrrentDirectory    ./          # official spelling has three r's
System.Name                 h2o_rf
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            0           # no cube files during optimization

Species.Number              2
<Definition.of.Atomic.Species
  O    O6.0-s2p2d1   O_PBE19
  H    H6.0-s2p1     H_PBE19
Definition.of.Atomic.Species>

Atoms.Number                3
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1   O    5.00000   5.00000   5.30000   3.0  3.0
  2   H    5.85000   5.00000   4.75000   0.5  0.5
  3   H    4.15000   5.00000   4.75000   0.5  0.5
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  10.0   0.0   0.0
   0.0  10.0   0.0
   0.0   0.0  10.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0       # K
scf.energycutoff            250.0       # Ry
scf.EigenvalueSolver        Cluster
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.01
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.maxIter                 100
scf.criterion               1.0e-9      # Hartree; tight SCF gives clean forces

#
# Geometry optimization
#
MD.Type                     RF          # Opt|DIIS|BFGS|RF|EF
MD.maxIter                  100         # default is 1 (= a single step!)
MD.Opt.criterion            1.0e-4      # Hartree/Bohr, max |force|; default 3.0e-4
MD.Opt.StartDIIS            5           # steepest-descent steps before RF starts (default 5)
MD.Opt.DIIS.History         3           # default 3
MD.Opt.EveryDIIS            200         # default 200
MD.Opt.Init.Hessian         Schlegel    # Schlegel|iden (Schlegel only for RF and EF)
```

要点：

- **SCF 收得紧**（`scf.criterion 1e-9`）：力是能量的导数，SCF 不充分时力有噪声，优化会在极小点附近来回跳；
- **截断能 250 Ry**：手册指出更高的 `scf.energycutoff` 使能量面更平滑、可加速优化，弱相互作用体系建议 300–400 Ry；
- `level.of.fileout 0`：优化过程中不写 cube 文件，节省时间和磁盘；
- 盒子 10 Å：分子尺寸约 1.7 Å，加上 2 × 3.2 Å 的 PAO 半径仍小于 10 Å。

运行与检查：

```bash
cd OpenMX/examples/06-geomopt-h2o-ch4
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
mpirun -np 4 openmx H2O_opt_RF.dat -nt 1 > h2o_rf.std
grep -A30 "History of" h2o_rf.out          # optimization history table
grep -o "scf_conv=[01]" h2o_rf.md | sort | uniq -c
ls H2O_opt_RF.dat# h2o_rf.xyz h2o_rf.md
```

`.out` 中的优化历史表格式如下（**示意输出，数值仅供参考**）：

```text
  MD_iter   SD_scaling     |Maximum force|   Maximum step        Utot
                           (Hartree/Bohr)        (Ang)         (Hartree)
      1     0.xxxxxxxx      0.0xxxxxxx        0.0xxxxxxx     -17.xxxxxxxx
      2     0.xxxxxxxx      0.0xxxxxxx        0.0xxxxxxx     -17.xxxxxxxx
    ...
     11     0.xxxxxxxx      0.0000xxxx        0.000xxxxx     -17.xxxxxxxx
```

判断优化是否真正收敛：

1. 最后一行的 `|Maximum force|` 小于 `MD.Opt.criterion`，而且步数小于 `MD.maxIter`（如果恰好停在 `MD.maxIter`，说明是被截断的）；
2. `Utot` 在最后几步基本不变；
3. `.md` 中所有步都是 `scf_conv=1`；
4. 几何合理：用 OpenMX Viewer 打开 `h2o_rf.md` 观察轨迹，或用本仓库的 [`omx_traj.py`](../scripts/omx_traj.py)（`python ../scripts/omx_traj.py -h`）画能量和最大力随步数的变化。

作为量级参考，官方手册表 3（同为 H2O，基组不同）中 PBE 的 r(OH) ≈ 0.971 Å、∠HOH ≈ 103.6–103.7°，实验值为 0.957 Å、104.5°。

### 3.3 输出文件

| 文件 | 内容 |
|---|---|
| `System.Name.out` | 输入回显、每步 SCF 历史、优化历史、最终能量、力、Mulliken 电荷 |
| `System.Name.md` | 每一步的结构（扩展 xyz 格式：晶胞、坐标、力、速度、净电荷、磁矩、`scf_conv`），OpenMX Viewer 可直接播放；dpdata 也支持这种格式 |
| `System.Name.xyz` | 最终结构 |
| `System.Name.ene` | 每步的计算量（具体列见源码 `iterout_md.c`） |
| `输入文件名#`（如 `H2O_opt_RF.dat#`） | 每步更新的续算输入：插入了最新坐标和 `scf.fixed.grid`（网格原点） |
| `System.Name_rst/` | SCF 重启文件，优化中自动用于密度外推 |

手册把续算文件写作 `System.Name.dat#`；从官方 NEB 示例的输出（`C2H4_NEB.dat#`）看，它实际以**输入文件名**加 `#` 命名，用 `ls *.dat#` 找到即可。

### 3.4 最速下降：甲烷例子

`Methane2_SD.dat` 把 C 原子从正四面体中心沿 x 方向移开 0.3 Å，用 `MD.Type Opt`：

```text
MD.Type                     Opt         # steepest descent
MD.maxIter                  200
MD.Opt.criterion            3.0e-4      # Hartree/Bohr (default)
```

```bash
mpirun -np 4 openmx Methane2_SD.dat -nt 1 > ch4_sd.std
grep -A60 "History of" ch4_sd.out | head -70
```

手册中同类例子的收敛曲线显示：最速下降会根据收敛情况放大步长因子，因此结构可能先**越过**驻点，之后再收敛回来，步数明显多于 RF。练习：把 `MD.Type` 改成 `RF`（`sed 's/^MD.Type .*/MD.Type RF/'`）比较步数。

### 3.5 约束：`MD.Fixed.XYZ`

```text
<MD.Fixed.XYZ
  1   1 1 1        # atom 1: x, y, z fixed
  2   1 0 0        # atom 2: only x fixed
  3   0 0 0
  4   0 0 0
MD.Fixed.XYZ>
```

- 第 1 列是原子序号（与坐标块一致），第 2–4 列是 x、y、z 的开关：**1 = 固定，0 = 放开**（注意与 VASP 的 `T/F` 含义相反）；
- 体系有 N 个原子就写 N 行；
- 对所有优化方法和分子动力学都有效，可以与 `MD.Init.Velocity` 同时使用；
- 典型用法：板模型固定底部 1–2 层（见 [第16章](16-surface-esm-workfunction.md)）；大体系只弛豫局部；固定一个远离活性区的原子以消除整体平移。

### 3.6 续算与疑难情况

**续算。** 优化在 `MD.maxIter` 或作业时限内没收敛时，直接用 `#` 文件继续：

```bash
cp H2O_opt_RF.dat# H2O_opt_RF_r2.dat
sed -i 's/^System.Name .*/System.Name                 h2o_rf_r2/' H2O_opt_RF_r2.dat   # optional
mpirun -np 4 openmx H2O_opt_RF_r2.dat -nt 1 > h2o_rf_r2.std
```

手册说明，OpenMX 每步都保存近似 Hessian，用 `#` 文件续算时会自动重用，所以分多次提交与一次跑完所需的总步数相同（3.9 起即如此）。这对有时限的共享集群很有用。

**收敛到鞍点。** 若怀疑优化停在鞍点（例如分子构型对称性异常高、频率计算出现虚频），增大 `MD.Opt.StartDIIS`（如 20–30），让最速下降先把结构带到极小点附近；或者给初始结构加一点小的随机扰动后重新优化。

**弱相互作用体系很难收敛。** 手册"Known issues"指出，分子体系等弱相互作用体系完全优化需要很多步。办法：

- `scf.energycutoff` 提高到 300–400 Ry，使能量面更平滑；
- `scf.criterion` 收紧到 1e-9 或 1e-10；
- 作为折中，把 `MD.Opt.criterion` 从 3×10⁻⁴ 放宽到 5×10⁻⁴。

**SCF 在优化中途不收敛。** 优化和 MD 会用前几步的密度外推初始密度（`scf.ExtCharge.History`，默认 2，手册建议保持默认）。若某一步 SCF 失败（`.md` 中出现 `scf_conv=0`），按 [第5章 §6](05-scf-kpoints-convergence.md) 调混合参数，并从 `#` 文件续算。

---

## 4. 晶胞优化

### 4.1 `MD.Type` 选项

| 取值 | 晶格 | 内坐标 | 约束 | 算法 |
|---|---|---|---|---|
| `OptC1` | 优化 | 固定 | 无 | 最速下降 |
| `OptC2` | 只优化长度 | 固定 | 晶格夹角固定 | 最速下降 |
| `OptC3` | 只优化长度 | 固定 | 夹角固定，三个矢量等长 | 最速下降 |
| `OptC4` | 只优化长度 | 固定 | 夹角固定，a1、a2 等长且与 a3 不同 | 最速下降 |
| `OptC5` | 优化 | 优化 | 无 | 最速下降 |
| `OptC6` | 优化（a3 固定） | 优化 | a3 固定 | 最速下降 |
| `OptC7` | 优化（a2、a3 固定） | 优化 | a2、a3 固定 | 最速下降 |
| `RFC5` | 优化 | 优化 | 无（可配合 `MD.Fixed.Cell.Vectors`） | RF + DIIS 混合，BFGS 更新 Hessian |
| `RFC6` | 优化（a3 固定） | 优化 | a3 固定 | 同 RFC5 |
| `RFC7` | 优化（a2、a3 固定） | 优化 | a2、a3 固定 | 同 RFC5 |

其余控制关键词（`MD.maxIter`、`MD.Opt.criterion`、`MD.Opt.StartDIIS`、`MD.Opt.DIIS.History`、`MD.Opt.Init.Hessian`）与原子优化完全相同。手册推荐用 **RFC5**：它已在大量晶体的高通量变胞优化中使用过；配合 `MD.Fixed.Cell.Vectors` 可以实现 RFC6、RFC7 的约束以及更一般的约束。

> [!IMPORTANT]
> **适用范围**：晶胞优化只支持**共线**计算（含 DFT+U），不能用于非共线/SOC。4.0 手册第19.1节写明晶胞优化不支持 DFT-D2/D3；而 3.9.7 补丁的发布说明曾称已支持含 DFT-D 的晶胞优化。两处说法不一致，保守起见，需要 vdW 的层状材料请用 §6 的 E(c) 扫描来确定层间晶格常数。

### 4.2 实战：Si 的 RFC5 优化

[`Si_RFC5_cellopt.dat`](../examples/06-cellopt-eos-si-c/README.md) 从一个故意扭曲的晶胞出发（三个晶格矢量长度和夹角都偏离 fcc 原胞，原子 2 偏离 (¼,¼,¼)），核心设置：

```text
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1   Si   0.00   0.00   0.00   2.0  2.0
  2   Si   0.27   0.24   0.26   2.0  2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  0.0000  2.6500  2.7800
  2.7400  0.0000  2.6900
  2.7600  2.7000  0.0500
Atoms.UnitVectors>

scf.energycutoff            250.0       # Ry
scf.Kgrid                   8 8 8
scf.criterion               1.0e-9      # Hartree

MD.Type                     RFC5        # OptC1..OptC7 | RFC5 | RFC6 | RFC7
MD.maxIter                  100         # default 1
MD.Opt.criterion            1.0e-4      # Hartree/Bohr, max |gradient| (atoms and cell)
MD.Opt.StartDIIS            5
MD.Opt.DIIS.History         3
MD.Opt.Init.Hessian         Schlegel
```

```bash
cd OpenMX/examples/06-cellopt-eos-si-c
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
mpirun -np 8 openmx Si_RFC5_cellopt.dat -nt 1 > si_rfc5.std
grep -A25 "History of cell optimization" si_rfc5.out
```

历史表的列为 `MD_iter`、`SD_scaling`、`|Maximum force|`（对原子坐标或晶格矢量的最大梯度，Hartree/Bohr）、`Maximum step`（Å）、`Utot`、`Enpy`（焓，Hartree）、`Volume`（Å³）。不加压时 `Enpy` 与 `Utot` 相同。收敛后：

- 最终晶格应恢复为 fcc 原胞（三个矢量等长、夹角 60°），原子 2 回到 (¼,¼,¼)；
- 体积对应的立方晶格常数 $a=(4V)^{1/3}$ 应接近 PBE 的约 5.47 Å；
- 用 §5 的 EvsLC 结果交叉验证：两种方法给出的 a0 应基本一致（相差 < 0.3%）。

最终晶格从 `si_rfc5.out` 末尾、`.md` 最后一帧的 `Cell_Vectors=` 字段或 `Si_RFC5_cellopt.dat#` 中读取：

```bash
tail -n 3 si_rfc5.md | head -1     # comment line of the last frame contains Cell_Vectors=
grep -A4 "<Atoms.UnitVectors" Si_RFC5_cellopt.dat#
```

> [!TIP]
> 晶胞变化较大时，用优化后的结构（`#` 文件）再优化一次，直到晶格不再变化。因为截断能固定时，实空间网格点数是按初始晶胞确定的，晶胞形状变化较大后，最好在新结构上重新确定网格再确认一次。

### 4.3 约束晶格分量：`MD.Fixed.Cell.Vectors`

```text
<MD.Fixed.Cell.Vectors
  0 0 1          # ax ay az
  0 0 1          # bx by bz
  1 1 1          # cx cy cz
MD.Fixed.Cell.Vectors>
```

三行分别对应 **a、b、c 三个矢量的笛卡尔分量**，1 = 固定在初始值，0 = 可变。上面的例子固定了整个 c 矢量和 a、b 的 z 分量，适合 c 沿真空方向的二维材料或板模型：只优化面内晶格，不让真空层塌缩。

> [!NOTE]
> **版本差异**：OpenMX 3.9 中 `MD.Fixed.Cell.Vectors` 只对 `OptC1`、`OptC2`、`OptC5` 有效；**`RFC5` 与 `MD.Fixed.Cell.Vectors` 的组合是 4.0 新增的**，也是 4.0 手册推荐的统一约束方式。用 3.9 时只能选 `RFC6`/`RFC7` 或最速下降类的 `OptC*`。

### 4.4 焓优化与应力输出

加压优化（最小化 $H=E+pV$）只需再加：

```text
MD.applied.pressure         10.0        # GPa, positive = compression (default 0)
MD.applied.pressure.flag    1 1 1       # default 1 1 1 (isotropic)
```

适用于 OptC1–OptC7 和 RFC5–RFC7。`MD.applied.pressure.flag` 可以只沿某个轴加压：正交晶系用对角形式的晶格矢量；六方晶系 `0 0 1` 表示只沿 c 方向（垂直于另外两个矢量所在平面）加压。历史表中 `Enpy` 列就是焓。

只想看当前结构的应力（对晶格矢量的梯度），不做优化：

```text
scf.stress.tensor           on          # default off; gradient w.r.t. lattice vectors is written to .out
```

---

## 5. 状态方程：EvsLC 扫描与 Birch–Murnaghan 拟合

### 5.1 关键词

```text
MD.Type                     EvsLC
MD.EvsLC.Step               0.5         # percent per step (default 0.4)
MD.EvsLC.flag               1 1 1       # per axis: 1 = expand, 0 = fixed, -1 = contract
MD.maxIter                  15          # number of lattice constants (default 1)
scf.Ngrid                   24 24 24    # REQUIRED in practice: fixed real-space grid
```

- 程序从输入的初始晶格出发，每一步按 `MD.EvsLC.Step`（百分比）均匀放大（flag = 1）或缩小（flag = −1）晶格矢量，共 `MD.maxIter` 步；flag = 0 的轴保持不变（层状材料可用 `0 0 1` 只扫 c）；
- 所以**初始晶格要取在平衡值的一侧**（例如压缩 3–5%），让扫描跨过极小点；
- 结果写在 `System.Name.EvsLC`，计算过程中就会逐步写入。每行前 9 列是 a、b、c 三个晶格矢量的 x、y、z 分量（Å），最后一列是总能（Hartree）。每一步实际用的晶格以这个文件为准；
- EvsLC 只均匀缩放晶格，**不弛豫内坐标**。对金刚石、fcc/bcc 金属等内坐标由对称性固定的结构没有问题；内坐标可变的结构见 §5.4。

> [!WARNING]
> 一定要用 `scf.Ngrid` 固定实空间网格。否则晶格常数变化时网格点数会离散跳变，E–V 曲线出现锯齿，拟合出的 B0 和 B0' 会严重失真（手册以 bcc Fe 为例展示了这一现象）。`scf.Ngrid` 的取法见 [第5章 §4.3](05-scf-kpoints-convergence.md)。

### 5.2 实战：金刚石的 E–V 曲线

`Cdia_EvsLC.dat` 从 a = 3.45 Å（比实验值 3.567 Å 小约 3.3%）开始，每步 0.5%，共 15 步，`scf.Ngrid 24 24 24`（在 a = 3.567 Å 时约相当于 375 Ry）：

```bash
mpirun -np 8 openmx Cdia_EvsLC.dat -nt 1 > cdia_evslc.std
cat cdia_evslc.EvsLC
python fit_evslc.py cdia_evslc.EvsLC --lattice fcc --natoms 2 --plot cdia_eos.png
```

`fit_evslc.py`（只依赖 numpy，`--plot` 需要 matplotlib）的输出格式如下（**示意输出，数值仅供参考**）：

```text
# 3rd-order Birch-Murnaghan fit of 15 points
V0  = 11.xxxx A^3/cell   (5.xxxx A^3/atom)
a0  = 3.5xxx A          (fcc, from V0)
E0  = -xxx.xxxxxx eV/cell   (-xxx.xxxxxx eV/atom)
B0  = 4xx.xx GPa
B0' = 3.xx
RMS residual = 0.0xx meV/cell
```

判断标准：

- V0 落在扫描范围中部（脚本在靠近边缘或超出时给出 WARNING）；
- RMS 残差远小于 1 meV/cell；若残差大且数据点呈锯齿，先检查是否忘了 `scf.Ngrid`；
- 量级参考：金刚石实验 a = 3.567 Å、B0 ≈ 443 GPa；PBE 文献值约 3.57–3.58 Å、约 430 GPa（Haas et al. 2009）。

也可以用 ASE 交叉验证：

```python
import numpy as np
from ase.eos import EquationOfState

d = np.loadtxt("cdia_evslc.EvsLC", comments="#")
cells = d[:, :9].reshape(-1, 3, 3)
vol = np.abs(np.linalg.det(cells))             # A^3 per cell
ene = d[:, 9] * 27.211386                       # Hartree -> eV
v0, e0, b0 = EquationOfState(vol, ene, eos="birchmurnaghan").fit()
print("a0 = %.4f A, B0 = %.1f GPa" % ((4 * v0) ** (1 / 3), b0 * 160.21766))
```

本仓库的 [`omx_eos_fit.py`](../scripts/omx_eos_fit.py) 提供"按晶格常数缩放生成输入（固定 `scf.Ngrid`）→ 运行 → 拟合"的一体化流程，用 `python ../scripts/omx_eos_fit.py -h` 查看用法。

### 5.3 `MD.Type DF`：Δ 因子用的 7 点扫描

```text
MD.Type      DF
```

OpenMX 会自动计算相对初始结构体积 −6、−4、−2、0、2、4、6% 的 7 个点，结果写入 `System.Name.DF`，可以直接用于 Δ 因子（delta gauge）的计算。同样建议固定 `scf.Ngrid`。

### 5.4 内坐标可变的晶体：手动扫描

六方 ZnO、钙钛矿、层状硫化物等结构，晶格改变时内坐标也要随之弛豫。做法是在每个固定体积（或固定 c/a）下用 `MD.Type RF` 优化原子，再收集能量拟合：

```bash
for s in 0.97 0.98 0.99 1.00 1.01 1.02 1.03; do
  # scale the three lattice vectors (FRAC coordinates follow automatically)
  awk -v s=$s 'BEGIN{inb=0}
       /^<Atoms.UnitVectors/{print; inb=1; next}
       /^Atoms.UnitVectors>/{inb=0}
       inb{printf "  %.6f %.6f %.6f\n", $1*s, $2*s, $3*s; next}
       {print}' base.dat | sed "s/^System.Name .*/System.Name   v_${s}/" > v_${s}.dat
  mpirun -np 8 openmx v_${s}.dat -nt 1 > v_${s}.std
done
```

`base.dat` 中设 `MD.Type RF`、`MD.maxIter 100`、固定的 `scf.Ngrid`，坐标用 `FRAC`。取每个结构优化最后一步的 `Utot` 和晶格（`#` 文件或 `.md` 最后一帧）组成 10 列的表，就能用同一个拟合脚本处理。

---

## 6. 范德华修正：DFT-D2 / DFT-D3

### 6.1 关键词

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `scf.dftD` | `off` | `on` 打开经验色散修正 |
| `version.dftD` | 2 | 2 = DFT-D2，3 = DFT-D3 |
| `DFTD3.damp` | `bj` | D3 的阻尼函数：`zero` 或 `bj`（Becke–Johnson） |
| `DFTD.Unit` | — | `Ang` 或 `AU`，决定截断半径的单位；手册建议 D3 用 `AU` |
| `DFTD.rcut_dftD` | 100 | 两体相互作用截断半径（`DFTD.Unit`） |
| `DFTD.cncut_dftD` | 40 | D3 配位数计算的截断（`DFTD.Unit`） |
| `DFTD.IntDirection` | `1 1 1` | 是否计入沿 a、b、c 方向周期镜像的 vdW 作用（1 计入，0 不计） |
| `<DFTD.periodicity` | 全部 1 | 逐原子指定是否周期：第 1 列原子序号，第 2 列 1（周期）或 0（非周期） |
| `DFTD.d`、`DFTD.scale6` | 20、0.75 | D2 参数（PBE） |
| `DFTD.scale6`、`DFTD.scale8`、`DFTD.sr6`、`DFTD.a1`、`DFTD.a2` | D3 的 PBE 默认值 | 用其他泛函时才需修改 |

默认参数对应 GGA-PBE。修正能量计入总能，在 `.out` 能量分解中显示为 `UvdW.`：

```bash
grep -E "UvdW\.|Utot\." graphite_d3.out
```

### 6.2 实战一：石墨的层间晶格常数（E(c) 扫描）

由于晶胞优化与 DFT-D 的兼容性在手册中写明不支持（§4.1），这里用 EvsLC 只扫描 c 轴。[`Graphite_D3_opt.dat`](../examples/06-vdw-d3-graphite/README.md) 的关键部分：

```text
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1   C   0.000000000   0.000000000   0.25   2.0  2.0
  2   C   0.333333333   0.666666667   0.25   2.0  2.0
  3   C   0.000000000   0.000000000   0.75   2.0  2.0
  4   C   0.666666667   0.333333333   0.75   2.0  2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
   2.464000   0.000000   0.000000
  -1.232000   2.133886   0.000000
   0.000000   0.000000   6.200000
Atoms.UnitVectors>

scf.Ngrid                   24 24 72    # fixed grid while c changes
scf.Kgrid                   15 15 5

scf.dftD                    on
version.dftD                3           # 2|3, default 2
DFTD3.damp                  bj          # zero|bj, default bj
DFTD.Unit                   AU          # the manual recommends AU for D3
DFTD.rcut_dftD              100.0       # in DFTD.Unit
DFTD.cncut_dftD             40          # in DFTD.Unit
DFTD.IntDirection           1 1 1       # bulk: vdW images along a, b and c

MD.Type                     EvsLC
MD.EvsLC.Step               1.0         # percent per step
MD.EvsLC.flag               0 0 1       # keep a and b, scale c only
MD.maxIter                  16
```

```bash
cd OpenMX/examples/06-vdw-d3-graphite
mpirun -np 8 openmx Graphite_D3_opt.dat -nt 1 > graphite_d3.std
python ../06-cellopt-eos-si-c/fit_evslc.py graphite_d3.EvsLC --mode c --natoms 4
```

对照实验：示例中的 `scan_c.sh` 在同一组 c 值上做单点计算，可以用 `DFTD=off` 关掉色散修正：

```bash
DFTD=on  MPI_RUN="mpirun -np 8" bash scan_c.sh
DFTD=off MPI_RUN="mpirun -np 8" bash scan_c.sh
python ../06-cellopt-eos-si-c/fit_evslc.py graphite_off_scan.EvsLC --mode c --natoms 4
```

预期：PBE-D3(BJ) 在实验值 c ≈ 6.71 Å 附近有明确的极小；纯 PBE 的能量随 c 增大几乎单调下降或只有极浅的极小，`fit_evslc.py` 会给出"最低点在扫描边缘"的警告。这正说明 GGA 本身描述不了层间色散作用。

### 6.3 实战二：双层石墨烯的层间距（原子优化 + D3）

板模型中原子位置可以自由弛豫，D3 可以直接和 RF 原子优化一起用。`Graphene_bilayer_D3_RF.dat` 从 3.70 Å 的层间距出发，c = 20 Å（真空约 16 Å）：

```text
scf.dftD                    on
version.dftD                3
DFTD3.damp                  bj
DFTD.Unit                   AU
DFTD.IntDirection           1 1 0       # no vdW images along the vacuum direction
MD.Type                     RF
MD.maxIter                  100
MD.Opt.criterion            1.0e-4
```

```bash
mpirun -np 8 openmx Graphene_bilayer_D3_RF.dat -nt 1 > blg_d3.std
grep -A30 "History of" blg_d3.out
tail -n 6 blg_d3.xyz                 # final structure; interlayer distance = z(layer 2) - z(layer 1)
```

`DFTD.IntDirection 1 1 0` 去掉了沿真空方向的镜像色散作用，这对表面、吸附、二维材料都应这样设置。优化后的层间距应收缩到 3.3–3.4 Å 量级（AB 堆垛石墨的层间距实验值约 3.35 Å）。

> [!CAUTION]
> 局域基组有 BSSE（[第4章 §8](04-pao-vps.md)）。对层间结合能、物理吸附能这类小量，应同时用空原子做 counterpoise 修正，手册第54章的苯二聚体例子也显示 BSSE 不可忽略。

---

## 7. NEB：最小能量路径与势垒

### 7.1 三步流程

手册把 NEB 计算分为三步，**三次计算必须使用完全相同的计算参数**（晶胞、截断能、基组、赝势、电子温度、k 点等），这样数值误差才能尽量抵消：

1. 优化反应物（初态）；
2. 优化产物（终态）；
3. 用两个优化结果构造 NEB 输入，优化连接二者的路径。

### 7.2 关键词

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `MD.Type` | — | `NEB` |
| `MD.NEB.Number.Images` | 10 | 两端点之间的 image 数，**不含**两个端点 |
| `MD.NEB.Spring.Const` | 0.1 Hartree/Bohr² | 弹簧常数，手册指出多数情况下路径对它不敏感 |
| `MD.Opt.DIIS.History` | 7（NEB） | DIIS+BFGS 混合优化参考的历史步数 |
| `MD.Opt.StartDIIS` | 5 | 先走几步最速下降 |
| `MD.maxIter` | 1 | 路径优化的最大步数 |
| `MD.Opt.criterion` | 3.0×10⁻⁴ Hartree/Bohr | 收敛判据 |
| `<Atoms.SpeciesAndCoordinates` | — | 反应物（image 0） |
| `<NEB.Atoms.SpeciesAndCoordinates` | — | 产物（image N+1） |
| `<NEB1.Atoms.SpeciesAndCoordinates` … `<NEBn.…` | — | 用户自定义初始路径时，给出每个 image 的坐标（同时设 `scf.restart on`） |
| `MD.NEB.Parallel.Number` | 自动 | 把 image 分成几组并行，节省内存（§7.5） |
| `MD.Fixed.XYZ` | — | 可以和 NEB 一起用 |

> [!WARNING]
> NEB 中体相体系的坐标**必须用 `Ang` 或 `AU`，不能用 `FRAC`**。FRAC 坐标会先被归一化到 0–1 之间，跨越边界的原子会被"折回"，image 在几何空间中不再连续，NEB 计算失败（手册51.8节）。

### 7.3 实战：Si 中间隙 H 在两个四面体位之间的扩散

[`examples/06-neb-h-in-si`](../examples/06-neb-h-in-si/README.md) 用 8 原子立方 Si 胞（a = 5.43 Å），H 从四面体间隙位 T1 = (½,½,½)a 扩散到相邻的 T2 = (¾,¾,¾)a，直线路径经过六角间隙位 (⅝,⅝,⅝)a。中性 H 共有 33 个电子，所以用 `scf.SpinPolarization on`，H 的初始占据写 `1.0 0.0`。

**第 1、2 步：优化端点。** `H_Si8_T1_opt.dat` 与 `H_Si8_T2_opt.dat` 除 H 的位置和 `System.Name` 外完全相同，都用 `MD.Type RF`、`MD.Opt.criterion 1.0e-4`：

```bash
cd OpenMX/examples/06-neb-h-in-si
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
mpirun -np 8 openmx H_Si8_T1_opt.dat -nt 1 > h_si8_t1.std
mpirun -np 8 openmx H_Si8_T2_opt.dat -nt 1 > h_si8_t2.std
grep "Utot\." h_si8_t1.out | tail -1; grep "Utot\." h_si8_t2.out | tail -1   # should be equal by symmetry
```

T1 与 T2 通过六角位处的反演中心互相对称，两个端点的能量应几乎相同（差值反映了数值误差，例如原子相对网格的位置不同）。

**第 3 步：构造 NEB 输入并运行。** `build_neb.py` 把两个 `#` 文件中的坐标分别填入模板 `H_diffusion_NEB.dat` 的 `<Atoms.SpeciesAndCoordinates` 和 `<NEB.Atoms.SpeciesAndCoordinates`，并检查单位、原子顺序，以及是否有原子在端点之间"穿过周期边界"（若有则自动按最小镜像平移回来）：

```bash
python build_neb.py                    # -> H_diffusion_NEB_relaxed.dat
mpirun -np 10 openmx H_diffusion_NEB_relaxed.dat -nt 1 > h_si8_neb.std
```

NEB 部分的关键词：

```text
MD.Type                     NEB
MD.NEB.Number.Images        5           # images between the end points (default 10)
MD.NEB.Spring.Const         0.1         # Hartree/Bohr^2 (default 0.1)
MD.maxIter                  100         # default 1
MD.Opt.criterion            3.0e-4      # Hartree/Bohr
MD.Opt.StartDIIS            10          # steepest descent before DIIS+BFGS
MD.Opt.DIIS.History         4           # NEB default is 7
```

### 7.4 输出与结果分析

| 文件 | 内容 |
|---|---|
| `h_si8_neb.neb.opt` | 路径优化历史：最大力（Maximum force）、最大步长、力的范数（Norm）、各 image 总能之和 |
| `h_si8_neb.neb.ene` | 每个 image 的总能（Hartree）、与相邻 image 的距离、与 image 0 的累计距离（Bohr） |
| `h_si8_neb.neb.xyz` | 各 image 的结构（xyz），可做动画 |
| `H_diffusion_NEB_relaxed.dat#` | 续算输入，含所有 image 的最新坐标 |
| `H_diffusion_NEB_relaxed.dat_0` … `_6` | 程序自动生成的各 image 输入（0 和 N+1 为端点） |
| `h_si8_neb_0.out` … `h_si8_neb_6.out` | 各 image 的输出，可分析电子结构沿路径的变化 |
| `h_si8_neb_#.DFTSCF` | 各 image 的 SCF 过程（临时文件） |

`.neb.ene` 的注释行说明了各列含义（第 1 列 image 序号，0 与 `MD.NEB.Number.Images+1` 为端点；第 2 列总能；第 3、4 列为距离）。从中计算势垒：

```bash
awk '!/^#/ && NF>=4 { if (n==0) e0=$2; n++; if (n==1 || $2>emax) emax=$2;
      printf "%3d  s=%8.4f A  dE=%8.4f eV\n", $1, $4*0.529177, ($2-e0)*27.211386 }
     END { printf "forward barrier = %.3f eV\n", (emax-e0)*27.211386 }' h_si8_neb.neb.ene
```

判断标准：

- `.neb.opt` 最后的 `|Maximum force|` 小于 `MD.Opt.criterion`；
- 势垒曲线平滑，image 在路径上大致均匀分布（`.neb.ene` 第 3 列近似相等）；
- 本例两端对称，能量曲线应关于中点（六角位）对称，最高点就在中点附近；
- **每个 image 的 SCF 都要收敛**。手册特别提醒：标准输出只显示 image 1 的信息，其余 image 是否收敛要自己检查：

```bash
for f in h_si8_neb_*.DFTSCF; do echo "== $f"; tail -n 3 "$f"; done
grep -L "scf_conv=1" h_si8_neb_*.md 2>/dev/null   # files without any converged step, if .md files exist
```

> [!NOTE]
> 本例是 8 原子的教学小胞，H 与其周期镜像只相距 5.43 Å，势垒数值不能直接与文献比较。实际研究要用 64 或 216 原子超胞，并考虑 H 的电荷态（H⁺、H⁰、H⁻ 的稳定位置不同）。OpenMX 手册中未提供爬坡图像（climbing image）开关；需要精确的鞍点能量时，可以增加 image 数，或在最高点附近再做一次更密的 NEB。

### 7.5 续算、自定义初始路径与并行

**续算**：路径优化达到 `MD.maxIter` 仍未收敛时，用每步更新的 `#` 文件（含所有 image 的坐标）作为新的输入继续：

```bash
cp H_diffusion_NEB_relaxed.dat# H_diffusion_NEB_r2.dat
mpirun -np 10 openmx H_diffusion_NEB_r2.dat -nt 1 > h_si8_neb_r2.std
```

**自定义初始路径**：默认初始路径是两端点间的直线插值。如果直线路径上原子靠得太近（例如分子穿过环、原子交换），就要显式给出每个 image 的坐标，方法与续算相同：

```text
<NEB1.Atoms.SpeciesAndCoordinates
  1   Si   ...
  ...
NEB1.Atoms.SpeciesAndCoordinates>
<NEB2.Atoms.SpeciesAndCoordinates
  ...
NEB2.Atoms.SpeciesAndCoordinates>
scf.restart     on
```

必须为 `MD.NEB.Number.Images` 指定的每个 image 都给出坐标，并设 `scf.restart on`。

**并行**：OpenMX 根据进程数和线程数自动分配并行方式，默认优先在 image 之间并行，进程数多于 image 数时再在每个 image 内部并行。手册建议 **MPI 进程数取 image 数的整数倍**（本例 5 个 image，用 10 或 20 个进程），MPI 进程数可以多于原子数。大体系在默认设置下若进程不够多，可能因内存不足异常退出，这时用：

```text
MD.NEB.Parallel.Number      3           # split the images into 3 groups computed one after another
```

image 被分成 3 组，组内并行，各组依次计算；`MD.NEB.Parallel.Number 1` 表示不在 image 之间并行，所有进程计算同一个 image，内存占用最小。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 何时调整 | 检查方法 |
|---|---|---|---|
| 原子优化方法 | `MD.Type RF` | 离极小点很远：先 `Opt` 或增大 `MD.Opt.StartDIIS` | `History of` 表 |
| `MD.maxIter` | 100–300 | 大体系分段续算 | 是否停在上限 |
| `MD.Opt.criterion` | 1×10⁻⁴ Ha/Bohr（默认 3×10⁻⁴） | 声子、弹性：更严；弱相互作用：可放宽到 5×10⁻⁴ | 最大力 |
| `scf.criterion`（优化中） | 1e-8 – 1e-9 | 力有噪声时收紧 | `scf_conv=1` |
| `scf.energycutoff`（优化中） | 比单点高 20–30% | 弱相互作用 300–400 Ry | 力的稳定性 |
| 晶胞优化 | `RFC5`（+ `MD.Fixed.Cell.Vectors`，4.0） | 二维材料固定 c | `History of cell optimization` |
| E–V 曲线 | `EvsLC` + `scf.Ngrid` | 内坐标可变：逐点 RF 优化 | 拟合残差 |
| vdW | `scf.dftD on`、`version.dftD 3`、`DFTD3.damp bj` | 板模型 `DFTD.IntDirection` 真空方向设 0 | `UvdW.` |
| NEB image 数 | 5–8 | 路径长或曲折时增加 | `.neb.ene` 平滑 |
| NEB 并行 | 进程数 = image 数的整数倍 | 内存不足：`MD.NEB.Parallel.Number` | — |

---

## 常见问题与报错

**1. 优化只走了一步就结束**

- 原因：`MD.maxIter` 默认是 1。
- 解决：显式写 `MD.maxIter 100`（或更大）。

**2. 优化在极小点附近来回振荡，`|Maximum force|` 降不下去**

- 原因：SCF 收敛不够导致力有噪声；截断能偏低使能量面不平滑；弱相互作用体系能量面太平。
- 解决：`scf.criterion` 收紧到 1e-9；`scf.energycutoff` 提高到 300–400 Ry；必要时 `MD.Opt.criterion` 放宽到 5×10⁻⁴。

**3. 优化"收敛"了，但结构对称性异常或频率有虚频**

- 原因：拟牛顿法开始得太早，收敛到了鞍点。
- 解决：增大 `MD.Opt.StartDIIS`，或加小扰动后用 `#` 文件重新优化。

**4. 作业到时限被杀，优化没完成**

- 原因：大体系步数多。
- 解决：用 `输入文件名#` 续算，近似 Hessian 会自动重用，不损失收敛速度。

**5. 晶胞优化不能用于非共线 / SOC 计算，或与 DFT-D 同用时结果可疑**

- 原因：手册第19.1节写明晶胞优化只支持共线计算（含 +U），并写明不支持 DFT-D2/D3。
- 解决：非共线体系先用共线计算优化晶格；需要 vdW 的层状材料用 E(c) 扫描（§6.2）。

**6. `MD.Fixed.Cell.Vectors` 加了却没有效果（3.9）**

- 原因：3.9 中该关键词只对 `OptC1`、`OptC2`、`OptC5` 有效，RFC5 的组合是 4.0 新增。
- 解决：升级到 4.0，或在 3.9 中改用 `RFC6`/`RFC7`/`OptC5`。

**7. E–V 曲线有锯齿、拟合出的 B0 不合理**

- 原因：没有固定 `scf.Ngrid`，网格点数随晶格跳变；或扫描范围没有跨过极小点。
- 解决：加 `scf.Ngrid`；调整初始晶格和 `MD.EvsLC.Step`，让 V0 落在扫描范围中部。

**8. NEB 路径中原子"飞"过整个晶胞，能量异常高**

- 原因：坐标用了 `FRAC`，或端点之间有原子穿过周期边界（一端在 0.01，另一端在 0.99）。
- 解决：坐标用 `Ang`/`AU`；端点坐标按最小镜像对齐（示例中的 `build_neb.py` 会自动检查）。

**9. NEB 中某些 image 的能量明显偏离光滑曲线**

- 原因：该 image 的 SCF 未收敛（标准输出只显示 image 1）。
- 解决：检查各 `System.Name_#.DFTSCF`；调整混合参数后用 `#` 文件续算。

**10. NEB 大体系计算中途因内存不足退出**

- 原因：默认优先在 image 之间并行，每个 image 分到的进程太少，内存未被分摊。
- 解决：设 `MD.NEB.Parallel.Number`（如 1–3），或增加进程数，并保持进程数是 image 数的整数倍。

---

## 练习

1. 用 `EF`、`BFGS`、`RF`、`Opt` 四种方法优化 `H2O_opt_RF.dat`，记录步数。注意 `BFGS` 时要把 `MD.Opt.Init.Hessian` 改成 `iden`（或删掉该行）。
2. 把 `Si_RFC5_cellopt.dat` 加上 `MD.applied.pressure 10.0`，比较体积和 `Enpy`。
3. 用 `MD.Type DF` 对 Si 做 7 点扫描，与 `EvsLC` 的拟合结果比较。
4. 在 `Cdia_EvsLC.dat` 中删掉 `scf.Ngrid`、改用 `scf.energycutoff 300`，重新扫描，观察 E–V 曲线和拟合残差的变化。
5. 把 NEB 的 image 数从 5 改为 7，势垒变化多少？进程数应怎样调整？
6. 对双层石墨烯，比较 `DFTD.IntDirection 1 1 1` 与 `1 1 0` 的总能和优化后的层间距。

---

## 参考资料

**官方文档**

- [OpenMX 4.0 手册 第18章 Geometry optimization](https://www.openmx-square.org/openmx_man4.0/s18_geom-opt.html) — Opt/DIIS/BFGS/RF/EF、初始 Hessian、约束、续算、md 文件格式
- [手册 第19章 Variable-cell optimization](https://www.openmx-square.org/openmx_man4.0/s19_vcell-opt.html) — OptC1–7、RFC5–7、`MD.Fixed.Cell.Vectors`、焓优化、应力
- [手册 第14章 Cutoff energy](https://www.openmx-square.org/openmx_man4.0/s14_cutoff-grid.html) — E–V 曲线跳变与 `scf.Ngrid`
- [手册 第17章 Restart](https://www.openmx-square.org/openmx_man4.0/s17_restart.html) — `#` 续算文件、密度外推
- [手册 第51章 NEB](https://www.openmx-square.org/openmx_man4.0/s51_neb.html) — 关键词、输出、自定义路径、并行与注意事项
- [手册 第54章 DFT-D2/D3](https://www.openmx-square.org/openmx_man4.0/s54_vdw-d2-d3.html) — vdW 关键词与苯二聚体例子
- [手册 第64章 Lattice constant vs. energy](https://www.openmx-square.org/openmx_man4.0/s64_lattice-energy.html) — `MD.Type EvsLC` 与 `DF`
- [手册 第79章 Known issues](https://www.openmx-square.org/openmx_man4.0/s79_known-issues.html) — 弱相互作用体系的优化
- [What's new](https://www.openmx-square.org/whatsnew/whatsnew.html) — 4.0 的 RFC5 + `MD.Fixed.Cell.Vectors`、3.9.7 补丁说明
- [OpenMX Viewer](https://www.openmx-square.org/viewer/) — 播放 `.md` 轨迹、查看 NEB 路径
- [ASE: Equation of state](https://docs.ase-lib.org/ase/eos.html)

**相关 GitHub 仓库（仅链接）**

- [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata) — 读取 OpenMX 的 `.dat` + `.md`，可把优化/MD 轨迹转为机器学习势训练数据
- [Tomoki-YAMASHITA/CrySPY](https://github.com/Tomoki-YAMASHITA/CrySPY) — 晶体结构预测，可调用 OpenMX 做结构弛豫
- [alamode-code/alamode](https://github.com/alamode-code/alamode) — 带 OpenMX 接口的声子/非谐工具，需要先得到充分优化的结构
- [issp-center-dev/abICS](https://github.com/issp-center-dev/abICS) — 构型采样，可生成 OpenMX 输入
- [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-) — 印尼语教程，含晶格优化与 NEB 分支
- [pranabdas/openmx](https://github.com/pranabdas/openmx) — 英文 OpenMX 教程站（无许可证，仅参考）
- ASE 的 OpenMX 计算器与 `ase.mep` NEB 工具见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)，更多仓库见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**

- G. Henkelman, H. Jónsson, *J. Chem. Phys.* 113, 9978 (2000)（改进切线 NEB）. DOI: [10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* 113, 9901 (2000)（climbing-image NEB）. DOI: [10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
- A. Banerjee, N. Adams, J. Simons, R. Shepard, *J. Phys. Chem.* 89, 52 (1985)（RF/EF 优化）. DOI: [10.1021/j100247a015](https://doi.org/10.1021/j100247a015)
- H. B. Schlegel, *Theor. Chim. Acta* 66, 333 (1984)（模型 Hessian）. DOI: [10.1007/BF00554788](https://doi.org/10.1007/BF00554788)
- F. Birch, *Phys. Rev.* 71, 809 (1947). DOI: [10.1103/PhysRev.71.809](https://doi.org/10.1103/PhysRev.71.809)
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 132, 154104 (2010)（DFT-D3）. DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* 32, 1456 (2011)（BJ 阻尼）. DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- S. Grimme, *J. Comput. Chem.* 27, 1787 (2006)（DFT-D2）. DOI: [10.1002/jcc.20495](https://doi.org/10.1002/jcc.20495)
- P. Haas, F. Tran, P. Blaha, *Phys. Rev. B* 79, 085104 (2009)（各泛函晶格常数与体模量基准）. DOI: [10.1103/PhysRevB.79.085104](https://doi.org/10.1103/PhysRevB.79.085104)

---

[← 上一章](05-scf-kpoints-convergence.md) | [返回目录](../README.md) | [下一章 →](07-md.md)
