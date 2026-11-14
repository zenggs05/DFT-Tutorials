# 第11章 过渡态搜索：NEB / CI-NEB / Dimer（VTST）（Transition-State Search in VASP: NEB, CI-NEB and Dimer with VTST）

本章讲解如何用 VASP + VTST（Henkelman 组的 Transition State Tools）搜索反应过渡态、计算反应能垒：从初末态准备、线性/IDPP 插值、NEB 与爬坡图像 CI-NEB（climbing-image nudged elastic band）的 INCAR 设置和两阶段收敛策略，到 nebresults.pl 能垒拟合、Dimer 方法精修鞍点，以及用频率计算验证"有且仅有一个虚频"。配套 Pt 吸附原子表面扩散和 Li 在 Si 中扩散两个完整案例，覆盖 VASP NEB 计算（transition state, energy barrier, diffusion）中最常见的坑。

---

## 1. 本章目标与适用场景

| 场景 | 推荐方法 | 本章小节 | 配套示例 |
|---|---|---|---|
| 初态和末态都已知（扩散、表面反应、解离） | NEB → CI-NEB | §4–§8 | [11-pt-adatom-neb](../examples/11-pt-adatom-neb/README.md)、[11-li-si-neb](../examples/11-li-si-neb/README.md) |
| 已有粗略 NEB，想精确收敛鞍点或换更高精度 | Dimer（neb2dim.pl） | §10 | [11-dimer-from-neb](../examples/11-dimer-from-neb/README.md) |
| 只知道初态，想找附近未知的鞍点 | Dimer（随机或指定方向） | §10 | — |
| 确认得到的是一阶鞍点、计算指前因子 | 频率分析（IBRION=5 / VTST 动力学矩阵） | §9 | [11-ts-freq-check](../examples/11-ts-freq-check/README.md) |
| 没有 VTST 版 VASP | VASP 内置 NEB（无爬坡）或 ASE NEB | §11 | — |

**前置知识**：结构优化（[第4章](04-relax-eos.md)）、slab 模型（[第9章](09-surface-workfunction.md)）、吸附与频率（[第10章](10-adsorption-catalysis.md)）。VTST 补丁的编译在[第1章](01-setup-compile.md)。作业调度与多节点并行见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)。

> [!NOTE]
> 本章以 VASP 6.4–6.6 + 对应版本的 vtstcode 为准（[henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode) 按 VASP 版本分目录，必须选与你的 VASP 版本匹配的那一份）。文中的输出片段除特别说明外都是**示意输出，数值仅供参考**。

---

## 2. 原理速览

### 2.1 势能面、最小能量路径与过渡态理论

反应的**最小能量路径**（MEP）连接初态（IS）与末态（FS）两个极小点，路径上的能量最高点是**一阶鞍点**——过渡态（TS）：沿反应坐标是极大，其余所有方向是极小，Hessian 恰好有一个负本征值（一个虚频）。

谐振子过渡态理论（hTST）给出速率常数：

$$
k = \nu^* \exp\left(-\frac{E_a}{k_BT}\right),\qquad
\nu^* = \frac{\prod_{i=1}^{3N}\nu_i^\mathrm{IS}}{\prod_{j=1}^{3N-1}\nu_j^\mathrm{TS}}
$$

其中 $E_a=E_\mathrm{TS}-E_\mathrm{IS}$，$\nu^*$ 为 Vineyard 指前因子（TS 的乘积中去掉虚频）。典型 $\nu^*\sim10^{12}$–$10^{13}$ s⁻¹。

### 2.2 NEB

在 IS 与 FS 之间放 $N$ 个"图像"（image），相邻图像之间用弹簧相连。第 $i$ 个图像受到的 NEB 力为：

$$
\mathbf F_i^\mathrm{NEB} = -\nabla V(\mathbf R_i)\big|_\perp + \mathbf F_i^\mathrm{s}\big|_\parallel
$$

- 只保留真实力**垂直于路径切线**的分量（让图像落到 MEP 上）；
- 只保留弹簧力**平行于切线**的分量（让图像沿路径等间距分布）。

这个"投影"（nudging）就是 NEB 与普通弹性带的区别。Henkelman 和 Jónsson 改进了切线定义（相邻能量较高图像方向），避免路径弯曲处出现"扭结"。

### 2.3 CI-NEB

普通 NEB 的图像未必正好落在鞍点上，能垒需要插值估计。**爬坡图像**（climbing image）让能量最高的图像不受弹簧力，并把真实力沿切线的分量反向：

$$
\mathbf F_\mathrm{CI} = -\nabla V + 2\left(\nabla V\cdot\hat\tau\right)\hat\tau
$$

于是该图像沿路径"爬"到能量极大、在其他方向上极小，收敛后**精确位于鞍点**。

### 2.4 Dimer

Dimer 方法只用一个构型：在中心点 $\mathbf R$ 两侧沿单位向量 $\hat N$ 各放一个相距 $\Delta R$ 的图像。先**旋转** dimer 使其沿最低曲率方向（相当于找 Hessian 最小本征向量，不需要算整个 Hessian），再**平移**：把沿 $\hat N$ 的力分量反向，使中心点沿最低曲率方向上坡、其他方向下坡，最终收敛到鞍点。它不需要末态，适合精修 NEB 结果或探索未知反应。

### 2.5 方法对比

| 方法 | 需要末态 | 图像数 | 精确鞍点 | 优点 | 局限 |
|---|---|---|---|---|---|
| NEB（VTST 或内置） | 是 | 3–8 | 否（需插值） | 给出整条 MEP | 鞍点能量有插值误差 |
| CI-NEB（VTST） | 是 | 3–8 | 是 | MEP + 精确 TS，最常用 | 两侧图像间距不等 |
| Dimer（VTST） | 否 | 1（2 次力计算/步） | 是 | 便宜、可探索未知路径 | 依赖初始方向，可能收敛到无关鞍点 |
| 频率分析 | — | — | 验证 | 判定鞍点阶数、求 $\nu^*$ | 计算量 ∝ 放开原子数 |

---

## 3. 准备工作：VTST 版 VASP 与 VTST 脚本

### 3.1 VTST 版 VASP

CI-NEB、Dimer 以及 `IOPT` 力优化器都来自 VTST 补丁，需要把 vtstcode 的源文件加入 VASP 源码并修改 `.objects`、`main.F` 等后重新编译（步骤见[第1章](01-setup-compile.md)与 [VTST 安装说明](https://theory.cm.utexas.edu/vtsttools/installation.html)）。

快速判断你的 VASP 是否带 VTST：

> [!TIP]
> 拿一个小体系（如稍微偏离平衡键长的 2 原子分子）设置 `IBRION = 3`、`POTIM = 0`、`IOPT = 3` 跑几步离子步。普通 VASP 中 `POTIM = 0` 意味着时间步长为零，原子不会移动（不同版本也可能直接报错退出）；VTST 版本中 `IOPT` 优化器接管，原子会移动、能量会下降。对比 OSZICAR 中各步的 `E0` 即可判断。

### 3.2 VTST 脚本

```bash
git clone https://github.com/henkelmangroup/vtstscripts.git ~/soft/vtstscripts
echo 'export PATH=$HOME/soft/vtstscripts:$PATH' >> ~/.bashrc && source ~/.bashrc
which nebmake.pl nebef.pl nebresults.pl neb2dim.pl dist.pl
```

脚本主要用 Perl（依赖同目录的 `Vasp.pm`），部分绘图脚本调用 gnuplot。本章用到的脚本：

| 脚本 | 用法 | 作用 |
|---|---|---|
| `nebmake.pl` | `nebmake.pl POSCAR1 POSCAR2 N` | 线性插值生成 00 … N+1 目录 |
| `dist.pl` | `dist.pl POSCAR1 POSCAR2` | 两个构型之间的距离（各原子位移平方和开方） |
| `nebmovie.pl` | `nebmovie.pl 0`（POSCAR）/ `1`（CONTCAR） | 生成 movie.xyz 动画 |
| `nebef.pl` | `nebef.pl` | 打印各图像的受力与能量（运行中监控） |
| `nebbarrier.pl` | `nebbarrier.pl` | 生成 neb.dat（反应坐标、能量、切向力） |
| `nebspline.pl` | `nebspline.pl` | 由 neb.dat 做基于力的三次样条，生成 spline.dat、exts.dat、mep.eps |
| `nebresults.pl` | `nebresults.pl` | 一次性运行 nebef/nebspline/nebmovie/nebconverge |
| `vfin.pl` | `vfin.pl 目录名` | 整理一次运行的结果到指定目录，并把 CONTCAR 覆盖到 POSCAR（便于续算） |
| `neb2dim.pl` | `neb2dim.pl` | 由 NEB 结果生成 Dimer 输入目录 saddle_dimer |
| `modemake.pl` | `modemake.pl POSCAR1 POSCAR2` | 两个构型之差生成 MODECAR |

> [!IMPORTANT]
> VTST 的 NEB 脚本要求 `00/` 和最后一个目录中**放有初、末态的 OUTCAR**，否则端点能量缺失，能垒曲线会出错。

仓库的 Python 辅助脚本 [`neb_tools.py`](../scripts/neb_tools.py) 也提供插值、检查原子顺序、收集能量和画能垒曲线等功能（`python ../scripts/neb_tools.py -h` 查看用法）。

---

## 4. 第一步：充分优化初态与末态

NEB 把弹性带"挂"在两个固定端点上，端点的任何误差都会传递到能垒。原则：

1. **相同设置**：IS、FS、NEB、Dimer、频率计算的 POTCAR、ENCUT、PREC、ISMEAR/SIGMA、k 点、ISPIN、LDIPOL、IVDW 完全一致；
2. **原子顺序一一对应**：FS 中第 $k$ 个原子必须是 IS 中第 $k$ 个原子移动后的位置（插值是逐原子进行的）；
3. **相同晶胞**：VASP/VTST 的普通 NEB 不变胞；
4. **相同的 Selective dynamics**；
5. **端点收敛更严**：`EDIFFG = -0.02`（比 NEB 的 −0.03～−0.05 更严）。

示例 [11-pt-adatom-neb](../examples/11-pt-adatom-neb/README.md) 的端点 INCAR：

```bash
SYSTEM  = Pt adatom on Pt(100) p(3x3) 4L - endpoint relaxation (IS or FS)

# ---- electronic part: IDENTICAL in IS, FS, NEB, dimer and frequency runs ----
PREC    = Normal
ENCUT   = 350          # ENMAX(Pt) is about 230 eV; check your POTCAR
EDIFF   = 1E-6
ALGO    = Fast
NELM    = 150
ISMEAR  = 1            # Methfessel-Paxton for the metal
SIGMA   = 0.20
ISPIN   = 1
LREAL   = Auto         # 37 atoms: real-space projectors are fine

# ---- ionic part ----
IBRION  = 2
ISIF    = 2
NSW     = 200
EDIFFG  = -0.02        # endpoints must be tighter than the NEB criterion

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

```bash
cd VASP/examples/11-pt-adatom-neb
export VASP_PP_PATH=/path/to/your/paw/library
mkdir -p IS FS
cp POSCAR_IS IS/POSCAR; cp POSCAR_FS FS/POSCAR
for d in IS FS; do cp INCAR_relax $d/INCAR; cp KPOINTS $d/; done
bash ../10-co-pt111-adsorption/make_potcar.sh IS FS
cp IS/POTCAR .                                   # NEB 根目录也需要 POTCAR
# 端点用单节点即可，例如：
(cd IS && mpirun -np 32 vasp_std > vasp.out)
(cd FS && mpirun -np 32 vasp_std > vasp.out)
grep "reached required accuracy" IS/OUTCAR FS/OUTCAR
```

### 4.1 检查端点

```bash
dist.pl IS/CONTCAR FS/CONTCAR                    # 总位移；吸附原子跳一个格点约 2.8 A
python idpp_interpolate.py IS/CONTCAR FS/CONTCAR --check-only
```

`idpp_interpolate.py --check-only` 会检查原子数、元素顺序、晶胞和 Selective dynamics 是否一致，并列出位移最大的几个原子（**示意输出**）：

```text
total IS->FS distance (sqrt of sum of squares): 2.8xx A
atoms that move most (1-based index, element, |d| in A):
     37 Pt    2.8xx
     22 Pt    0.0xx
     ...
```

> [!WARNING]
> 如果"位移最大"的原子不止一个且都在 3 Å 以上，而你的反应只涉及一个原子迁移，几乎可以肯定是**原子顺序不一致**（例如 FS 是用另一个脚本重新生成的，原子被重新排序了）。必须手工调整 FS 的原子顺序，否则插值路径会让多个原子"穿过"彼此。

**周期性边界**：如果某个原子在 IS 中分数坐标为 0.98、在 FS 中为 0.02（实际只移动了 0.04），简单线性插值会让它横穿整个晶胞。`idpp_interpolate.py` 用最小镜像约定处理这种情况；使用其他工具时要先手工把坐标"展开"。

---

## 5. 第二步：生成中间图像

### 5.1 图像数怎么选

- 经验上相邻图像间距 ≲ 0.8–1 Å；单原子跳跃 2.8 Å 用 3 个图像（间距约 0.7 Å）即可；
- 官方 Wiki 建议**图像数越少越好**：可以先用 1–3 个图像收敛，再增加；
- 路径中如果存在中间极小（两个鞍点），应拆成两段 NEB 分别计算；
- 对称路径（如本例 hollow→bridge→hollow）用**奇数**个图像，中间图像恰好从对称鞍点附近开始。

### 5.2 线性插值：nebmake.pl

```bash
nebmake.pl IS/CONTCAR FS/CONTCAR 3     # 生成 00 01 02 03 04
cp IS/OUTCAR 00/ ; cp FS/OUTCAR 04/
nebmovie.pl 0                          # movie.xyz，用 VESTA/OVITO/Jmol 查看
```

### 5.3 IDPP 插值（推荐用于多原子同时运动）

线性插值在多个原子协同运动、或分子转动时容易让原子靠得过近，导致前几步受力巨大甚至 SCF 失败。**IDPP**（image dependent pair potential，Smidstrup 等，[J. Chem. Phys. 2014](https://doi.org/10.1063/1.4878664)）先在线性插值的原子间距上构造一个目标势，优化后得到更合理的初始路径。示例脚本 [idpp_interpolate.py](../examples/11-pt-adatom-neb/idpp_interpolate.py) 基于 ASE：

```bash
python idpp_interpolate.py IS/CONTCAR FS/CONTCAR -n 3              # IDPP（默认）
python idpp_interpolate.py IS/CONTCAR FS/CONTCAR -n 3 --method linear
```

脚本会写出 `00/POSCAR … 04/POSCAR`，保留 Selective dynamics，打印相邻图像间距，并警告原子间距 < 1 Å 的图像（**示意输出**）：

```text
wrote 00 ... 04 (idpp interpolation, IMAGES = 3)
distance between neighbouring images (A):
   00 -> 01 :  0.702   (largest single-atom move 0.702)
   01 -> 02 :  0.702   (largest single-atom move 0.702)
   ...
```

核心代码只有几行（ASE ≥ 3.23 中 NEB 位于 `ase.mep`，旧版为 `ase.neb`）：

```python
from ase.io import read
from ase.mep import NEB          # ASE < 3.23: from ase.neb import NEB
ini, fin = read("IS/CONTCAR"), read("FS/CONTCAR")
images = [ini] + [ini.copy() for _ in range(3)] + [fin]
NEB(images).interpolate(method="idpp", mic=True)
```

示例中的 `make_images.sh` 把这几步合在一起（插值 + 复制端点 OUTCAR + dist.pl 检查）：

```bash
bash make_images.sh IS FS 3            # IDPP
bash make_images.sh IS FS 3 linear     # 线性（优先调用 nebmake.pl）
```

最终的目录结构：

```text
11-pt-adatom-neb/
├── INCAR            # NEB 的 INCAR（含 IMAGES = 3）
├── KPOINTS  POTCAR  # 所有图像共用
├── 00/  POSCAR OUTCAR   # 初态（不参与计算，OUTCAR 供脚本读取能量）
├── 01/  POSCAR
├── 02/  POSCAR
├── 03/  POSCAR
└── 04/  POSCAR OUTCAR   # 末态
```

---

## 6. 第三步：NEB 的 INCAR 详解

### 6.1 第一阶段：普通 NEB，粗收敛

```bash
SYSTEM  = Pt adatom hop on Pt(100) - NEB stage 1: plain NEB, coarse convergence (VTST)

# ---- electronic part: identical to INCAR_relax ----
PREC    = Normal
ENCUT   = 350
EDIFF   = 1E-6
ALGO    = Fast
NELM    = 150
ISMEAR  = 1
SIGMA   = 0.20
ISPIN   = 1
LREAL   = Auto

# ---- nudged elastic band ----
IMAGES  = 3            # intermediate images = directories 01 02 03; MPI ranks must be a multiple of 3
SPRING  = -5           # spring constant (eV/A^2); negative value switches on nudging
ICHAIN  = 0            # 0 = NEB (VTST)
LCLIMB  = .FALSE.      # VTST default is .TRUE.; keep the climbing image off in stage 1

# ---- VTST force-based optimiser ----
IBRION  = 3            # with POTIM = 0 the VASP optimiser is disabled ...
POTIM   = 0            # ... and IOPT selects the VTST optimiser
IOPT    = 3            # 3 = quick-min: robust when the band is far from the MEP
ISIF    = 2
NSW     = 150
EDIFFG  = -0.10        # loose: only a rough path is needed before climbing

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `IMAGES` | 3–8 | 中间图像数（不含端点）；**MPI 进程总数必须是 IMAGES 的整数倍** |
| `SPRING` | −5 | 弹簧常数（eV/Å²），负值表示开启 nudging；VTST 与官方默认都是 −5，过大反而难收敛 |
| `ICHAIN` | 0 | VTST 方法选择：0 = NEB，1 = 动力学矩阵，2 = Dimer，3 = Lanczos |
| `LCLIMB` | 第一阶段 `.FALSE.`，第二阶段 `.TRUE.` | **VTST 中默认值是 `.TRUE.`**，想跑普通 NEB 必须显式关掉 |
| `IBRION` | 3 | 与 `POTIM = 0` 一起关闭 VASP 自带优化器 |
| `POTIM` | 0 | 必须为 0，否则 VASP 的阻尼 MD 与 VTST 优化器同时移动原子 |
| `IOPT` | 3（QM）/ 7（FIRE）→ 1（LBFGS）/ 2（CG） | 见下表 |
| `EDIFFG` | 阶段一 −0.1～−0.05；阶段二 −0.03 | 必须为负（按力判据）；NEB 的力是投影后的力 |
| `EDIFF` | 1E-6（TS 精修/Dimer 用 1E-7） | 力不准会导致优化器振荡 |
| `NSW` | 100–300 | NEB 收敛通常需要几十到几百步 |

**VTST 优化器**（来自 [VTST optimizers 文档](https://theory.cm.utexas.edu/vtsttools/optimizers.html)）：

| `IOPT` | 方法 | 特点 / 何时用 |
|---|---|---|
| 0 | 使用 IBRION 指定的 VASP 优化器（默认） | 不用 VTST 优化器时 |
| 1 | LBFGS | 可对整条带做全局优化（`LGLOBAL = .TRUE.` 默认）；力准确时快 |
| 2 | CG（基于力的牛顿线搜索共轭梯度） | 力准确时快，Dimer 推荐 |
| 3 | QM（quick-min） | 最稳健、对新手友好，远离 MEP 时用 |
| 4 | SD（最速下降） | 主要用于测试 |
| 7 | FIRE | 与 QM 类似但通常更快，适合受力大或力有噪声的阶段 |
| 8 | ML-PyAMFF（vtstcode6.3 起，实验性） | 机器学习加速优化 |

常用优化器参数：`MAXMOVE`（每步最大位移，默认 0.2 Å）、`TIMESTEP`（QM/FIRE 的时间步，默认 0.1）、`ILBFGSMEM`（LBFGS 记忆步数，默认 20）。

> [!CAUTION]
> 不要在 NEB 中使用 `IBRION = 2`（共轭梯度）：NEB 投影后的力与能量不再一致，基于能量的线搜索会失败。VASP 自带优化器中只有 `IBRION = 1`（准牛顿）和 `IBRION = 3`（阻尼 MD）是纯力驱动的；用 VTST 时统一用 `IBRION = 3, POTIM = 0` + `IOPT`。

### 6.2 第二阶段：打开爬坡图像

```bash
# ---- climbing-image NEB ----
IMAGES  = 3
SPRING  = -5
ICHAIN  = 0
LCLIMB  = .TRUE.       # highest image climbs to the saddle point

# ---- VTST optimiser ----
IBRION  = 3
POTIM   = 0
IOPT    = 1            # 1 = LBFGS (global); use 7 (FIRE) or 3 (QM) if it oscillates
ISIF    = 2
NSW     = 300
EDIFFG  = -0.03        # max force on every image < 0.03 eV/A
```

（完整文件见示例中的 `INCAR_neb_stage2`，电子部分与第一阶段相同。）

**为什么要分两阶段？** VTST 文档明确指出：如果一开始最高图像离鞍点很远就打开 climbing，鞍点两侧的图像间距会变得非常不均匀，路径质量变差。先用普通 NEB 把整条带拉到 MEP 附近，再让最高图像爬坡，收敛更快也更可靠。

### 6.3 并行设置

```bash
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32          # 96 ranks = 3 images x 32 ranks
mpirun -np 96 vasp_std > vasp.out
```

- VASP 把 MPI 进程平均分给各图像（每个图像 96/3 = 32 个进程），图像之间几乎不通信，并行效率很高；
- `NCORE`、`KPAR` 作用于**每个图像内部**：本例每图像 32 进程，`NCORE = 4` 即可；
- 示例的 [job.slurm](../examples/11-pt-adatom-neb/job.slurm) 会在提交前检查进程数能否被 IMAGES 整除。

---

## 7. 第四步：运行、监控与两阶段切换

### 7.1 运行与监控

```bash
cp INCAR_neb_stage1 INCAR
sbatch job.slurm

# 运行中监控（在 NEB 根目录）
nebef.pl                                  # 每个图像的受力与能量
for i in 01 02 03; do echo -n "$i "; grep "energy(sigma->0)" $i/OUTCAR | tail -1 | awk '{print $NF}'; done
tail -3 01/OSZICAR                        # 每个图像有自己的 OSZICAR、OUTCAR、CONTCAR
```

`nebef.pl` 的输出列依次为图像号、最大受力、能量、相对 00 的能量（**示意输出，数值仅供参考**）：

```text
  0    0.012   -224.53xxx    0.000000
  1    0.08x   -224.3xxxx    0.2xxxxx
  2    0.09x   -224.0xxxx    0.5xxxxx
  3    0.07x   -224.3xxxx    0.2xxxxx
  4    0.013   -224.53xxx    0.00xxxx
```

收敛的标志：所有中间图像的受力都低于 `|EDIFFG|`，每个图像的 OUTCAR 中出现 `reached required accuracy`。

### 7.2 切换到 CI-NEB

第一阶段达到 −0.1 eV/Å 量级后（不必完全收敛），用最后的结构作为第二阶段起点：

```bash
vfin.pl stage1          # 结果整理到 stage1/，并把各图像 CONTCAR 覆盖到 POSCAR
# 不用 vfin.pl 时：for i in 01 02 03; do cp $i/CONTCAR $i/POSCAR; done
cp INCAR_neb_stage2 INCAR
sbatch job.slurm
```

> [!TIP]
> 第二阶段如果 LBFGS 出现能量/受力大幅振荡，换 `IOPT = 7`（FIRE）或 `IOPT = 3`，或减小 `MAXMOVE`（如 0.1）。如果是 SCF 噪声导致，把 `EDIFF` 收紧到 1E-7。

---

## 8. 第五步：结果分析与能垒

### 8.1 nebresults.pl

计算结束后在 NEB 根目录运行：

```bash
nebresults.pl
ls neb.dat spline.dat exts.dat mep.eps movie.xyz
cat exts.dat
```

| 文件 | 内容 |
|---|---|
| `neb.dat` | 每个图像的累计反应坐标（Å）、相对能量（eV）、沿带的切向力 |
| `spline.dat` | 用能量和切向力做的三次样条插值曲线 |
| `exts.dat` | 样条曲线上所有极值点的位置与能量（鞍点能量即能垒） |
| `mep.eps` | MEP 图 |
| `movie.xyz` | 路径动画 |

正向能垒 $E_a^\mathrm{f}=E_\mathrm{TS}-E_\mathrm{IS}$，逆向能垒 $E_a^\mathrm{r}=E_\mathrm{TS}-E_\mathrm{FS}$，反应能 $\Delta E=E_\mathrm{FS}-E_\mathrm{IS}$。CI-NEB 收敛后，爬坡图像本身的能量就是鞍点能量，与样条极值应基本一致；差别较大时说明图像太少或未收敛。

### 8.2 用 Python 画能垒（可选）

```python
from ase.io import read
from ase.mep import NEBTools        # ASE < 3.23: from ase.neb import NEBTools
import matplotlib.pyplot as plt

images = [read(f"{i:02d}/OUTCAR") for i in range(5)]   # 00 与 04 的 OUTCAR 必须存在
nt = NEBTools(images)
barrier, dE = nt.get_barrier()                         # 基于能量+力的样条拟合
print("barrier = %.3f eV, reaction energy = %.3f eV" % (barrier, dE))
fig, ax = plt.subplots()
nt.plot_band(ax)
fig.savefig("neb_barrier.png", dpi=200)
```

也可以使用仓库脚本 [`neb_tools.py`](../scripts/neb_tools.py)（用法见 `-h`），或 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md) 中的出版级作图方法。

### 8.3 检查清单

- [ ] 爬坡图像是能量最高的图像，且只有一个能量极大；
- [ ] 路径动画中原子运动合理，没有原子穿越或异常靠近；
- [ ] 端点能量与单独优化的 IS/FS 一致（00/04 的 OUTCAR 是对的）；
- [ ] 能垒对图像数不敏感（3 → 5 个图像时能垒变化 < 0.02–0.05 eV）；
- [ ] 下一节的频率分析显示 TS 只有一个虚频。

---

## 9. 第六步：验证过渡态——频率分析

### 9.1 IBRION = 5 有限差分

把爬坡图像的 CONTCAR 拿出来，放开参与反应的原子（及其近邻），做有限差分频率计算（示例 [11-ts-freq-check](../examples/11-ts-freq-check/README.md)）：

```bash
mkdir -p ts_freq && cd ts_freq
# 本例中第 37 号原子是吸附原子；先只放开它（3 个模，最快的检查）
python ../../10-co-pt111-adsorption/set_freq_selective.py ../02/CONTCAR -i 37 -o POSCAR
# 更完整：放开吸附原子 + 表面层（z 高于约 13 A 的原子，按你的结构调整）
# python ../../10-co-pt111-adsorption/set_freq_selective.py ../02/CONTCAR -i 37 --zmin 13.0 -o POSCAR
cp ../../11-ts-freq-check/INCAR ../KPOINTS ../POTCAR .
mpirun -np 32 vasp_std > vasp.out
grep -E "f/i|f  =" OUTCAR | head -12
```

INCAR 要点：

```bash
EDIFF   = 1E-7         # tighter than the NEB
IBRION  = 5            # selective dynamics is honoured only by IBRION = 5
POTIM   = 0.015
NFREE   = 2
NSW     = 1
```

判断标准（**示意输出，数值仅供参考**）：

```text
   1 f  =    x.xxxxxx THz ...  1xx.xxxxxx cm-1 ...
   2 f  =    x.xxxxxx THz ...   xx.xxxxxx cm-1 ...
   3 f/i=    x.xxxxxx THz ...  1xx.xxxxxx cm-1 ...
```

- **一阶鞍点：有且仅有一个明显的虚频**（`f/i`），其本征向量指向反应方向（本例为吸附原子沿跳跃方向的运动）；
- 两个或更多明显虚频：是高阶鞍点，沿第二个虚频方向位移后重新搜索；
- 没有虚频：并非鞍点（可能是 CI-NEB 未收敛或爬坡图像选错）；
- 很小的虚频（< ~50 cm⁻¹）通常是数值噪声或放开原子太少导致，需要提高精度或放开更多原子再确认。

振动模式可视化：[VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF) 可把 OUTCAR 中的模式转换成 XSF 或动画，直观检查虚频方向。

### 9.2 VTST 动力学矩阵方法与指前因子

VTST 还提供 `ICHAIN = 1` 的动力学矩阵计算：用 `dymselsph.pl`（以某原子为中心按半径选原子）或 `dymseldsp.pl`（选 IS→TS 位移最大的原子）生成 DISPLACECAR，计算后用 `dymmatrix.pl` 得到频率（freq.dat），再用 `dymprefactor.pl freq.dat(IS) freq.dat(TS)` 计算 Vineyard 指前因子。必要参数见 [VTST Dynamical Matrix 文档](https://theory.cm.utexas.edu/vtsttools/dymmat.html)（`IBRION = 3`、`POTIM = 0`、`EDIFF = 1E-8`、`EDIFFG` 设得很小以防提前停止、`NSW` = 自由度数 + 1）。

用 IBRION=5 的结果手算也很简单：

$$
\nu^* = \frac{\prod_{i}\nu_i^\mathrm{IS}}{\prod_{j\neq \mathrm{imag}}\nu_j^\mathrm{TS}}
$$

IS 与 TS 必须放开**相同的原子**，且 IS 有 $3N_\mathrm{free}$ 个实频、TS 有 $3N_\mathrm{free}-1$ 个实频。

---

## 10. Dimer 方法

### 10.1 从 NEB 生成 Dimer 输入（neb2dim.pl）

Dimer 只需一个构型，计算量约为每步 2 次力计算，非常适合：①把一个粗收敛的 NEB 精修成精确鞍点；②在更高 ENCUT/更密 k 点下重新收敛鞍点（不必重跑整条 NEB）。

```bash
# 在 NEB 根目录（已完成一个 NEB 运行）
vfin.pl neb_final        # 构型写回各图像的 POSCAR
nebresults.pl            # 生成 exts.dat（neb2dim.pl 用它找能量最高点）
neb2dim.pl               # -> saddle_dimer/ （含 POSCAR 与 MODECAR）
# 或指定图像：neb2dim.pl 2
cd saddle_dimer
cp ../../11-dimer-from-neb/INCAR ../KPOINTS ../POTCAR .
```

`neb2dim.pl` 在能量最高点附近插值得到中心构型（POSCAR），并利用 MEP 的切线方向给出初始 dimer 方向（MODECAR）。没有 NEB 时也可以用 `modemake.pl POSCAR1 POSCAR2` 由两个构型之差生成 MODECAR。MODECAR 的格式是**每个原子一行、三个笛卡尔分量**（原子顺序与 POSCAR 一致），表示一个单位向量；其中只要反应相关原子的分量足够大即可，通常无需手写。

### 10.2 Dimer 的 INCAR

```bash
SYSTEM  = Pt adatom on Pt(100) - dimer saddle-point refinement (VTST, ICHAIN = 2)
# Start from saddle_dimer/ produced by neb2dim.pl (POSCAR + MODECAR).
# Do NOT set IMAGES here: the dimer works on a single configuration.

# ---- electronic part: same as the NEB, but tighter SCF (curvature needs accurate forces) ----
PREC    = Normal
ENCUT   = 350
EDIFF   = 1E-7
ALGO    = Fast
NELM    = 150
ISMEAR  = 1
SIGMA   = 0.20
ISPIN   = 1
LREAL   = Auto

# ---- dimer method ----
ICHAIN  = 2            # 2 = dimer
DdR     = 0.005        # dimer separation (A), default; up to 0.01 if the torque is noisy
DRotMax = 4            # max rotations per translation step (default 1)
DFNMin  = 0.01         # rotational force below which the dimer is not rotated (default)
DFNMax  = 1.0          # rotational force below which rotation stops (default)

# ---- VTST optimiser (required) ----
IBRION  = 3
POTIM   = 0
IOPT    = 2            # 2 = CG, recommended for the dimer by the VTST documentation
ISIF    = 2
NSW     = 300
EDIFFG  = -0.01        # converged when max force < 0.01 eV/A AND curvature < 0

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

| 关键词 | 默认值 | 说明 |
|---|---|---|
| `ICHAIN` | — | 必须为 2（新版 VTST 要求显式设置） |
| `DdR` | 5E-3 Å | dimer 两像间距（中心到每个像的两倍）；力噪声大时可增到 0.01 |
| `DRotMax` | 1 | 每次平移前最多旋转次数；适当增大（3–5）可让方向更准 |
| `DFNMin` | 0.01 | 旋转力低于此值不旋转 |
| `DFNMax` | 1.0 | 旋转力低于此值停止旋转；介于两者之间至少旋转一次 |
| `IOPT` | — | Dimer 必须设置 IOPT（官方推荐 2），否则作业报错退出 |

### 10.3 输出与收敛判断

`DIMCAR` 每行包含：Step、Force（最大受力）、Torque（旋转力）、Energy（中心点能量，**不是**外推到零展宽的鞍点能量）、Curvature（沿 dimer 的曲率）、Angle（旋转角）。

- 收敛条件：Force < |EDIFFG| **且** Curvature < 0；
- Torque 应随旋转系统性下降；若不下降，说明力不够准——把 `EDIFF` 降到 1E-7～1E-8，或把 `DdR` 增大到 0.01；
- Angle 一开始可能较大，应逐渐变小并保持很小；持续接近 1 说明有问题；
- Curvature 为正表示还没进入鞍点区域，继续运行即可；
- 鞍点能量从 OUTCAR 中读取（`energy(sigma->0)`），而不是 DIMCAR 的 Energy 列。

续算：

```bash
cp CENTCAR POSCAR          # 当前 dimer 中心
cp NEWMODECAR MODECAR      # 当前最低曲率方向
```

收敛后对中心构型做 §9 的频率分析确认。`dimmins.pl` 可以沿鞍点两侧各推一步生成两个极小构型，用于确认鞍点连接的是哪两个极小点。Lanczos 方法（`ICHAIN = 3`，`neb2lan.pl`）与 Dimer 类似，见 VTST 文档。

---

## 11. 没有 VTST 时：VASP 内置 NEB 与 ASE NEB

### 11.1 VASP 内置 NEB

官方 VASP 本身支持普通 NEB（[Nudged elastic bands](https://www.vasp.at/wiki/index.php/Nudged_elastic_bands)）：

```bash
IMAGES  = 3
SPRING  = -5           # negative: nudged elastic band
IBRION  = 3            # or 1; do not use 2 (CG) with NEB
POTIM   = 0.1          # damped-MD step for IBRION = 3 (test 0.05-0.2)
SMASS   = 2            # official advice: IBRION=3 + SMASS=2 to pre-converge
EDIFFG  = -0.05
NSW     = 200
```

- 目录结构、进程数为 IMAGES 整数倍等要求与 VTST 相同；
- **没有爬坡图像**（VTST 文档说明官方实现不含 climbing image 和改进切线，`LCLIMB` 只在 VTST 中可用），鞍点能量需靠插值估计，或在最高点附近加密图像；
- 官方 Wiki 也指出内置 NEB **不能处理变胞的结构相变**；VTST 提供 SS-NEB（`LNEBCELL = .TRUE.`，配合 `ISIF = 3`、`IOPT = 3`）。

### 11.2 ASE NEB + VASP 计算器

不能重新编译 VASP 时，也可以在 ASE 中运行 NEB/CI-NEB，由 ASE 的 `Vasp` 计算器逐图像调用 VASP（每个力计算都是一个独立的 VASP 单点）。优点是不需要补丁、可以用 ASE 的各种优化器和 AutoNEB；缺点是效率较低（每步重启 VASP，可读 WAVECAR 部分缓解）。ASE 驱动 VASP 的方法见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。

---

## 12. 案例一：Pt 吸附原子在 Pt(100) 上的跳跃扩散

**体系**：Pt(100) p(3×3)、4 层（底部 2 层固定），一个 Pt 吸附原子从 4 重空位（hollow）跳到相邻的 hollow，经过 bridge 位。路径对称，3 个图像时中间图像（02）就从 bridge 位开始。

**完整流程**（示例 [11-pt-adatom-neb](../examples/11-pt-adatom-neb/README.md)）：

```bash
cd VASP/examples/11-pt-adatom-neb
# 1) 端点优化（§4）
mkdir -p IS FS && cp POSCAR_IS IS/POSCAR && cp POSCAR_FS FS/POSCAR
for d in IS FS; do cp INCAR_relax $d/INCAR; cp KPOINTS $d/; done
bash ../10-co-pt111-adsorption/make_potcar.sh IS FS && cp IS/POTCAR .
(cd IS && sbatch ../../10-co-pt111-adsorption/job.slurm)
(cd FS && sbatch ../../10-co-pt111-adsorption/job.slurm)
# 2) 插值（§5）
bash make_images.sh IS FS 3
# 3) NEB 第一阶段（§6–7）
cp INCAR_neb_stage1 INCAR && sbatch job.slurm
# 4) CI-NEB 第二阶段
vfin.pl stage1 && cp INCAR_neb_stage2 INCAR && sbatch job.slurm
# 5) 结果
nebresults.pl && cat exts.dat
```

**预期**：能量最高的是 02（bridge 位），能垒曲线对称，正逆向能垒相同。用 Dimer（§10）和频率分析（§9）可进一步确认。

> [!NOTE]
> 场离子显微镜实验表明 Pt 在 Pt(001) 上的自扩散主要通过**交换机制**（吸附原子挤入表面层、把一个表面原子顶出来）完成，而非直接跳跃（Kellogg 与 Feibelman，[Phys. Rev. Lett. 1990](https://doi.org/10.1103/PhysRevLett.64.3143)）。本例的跳跃路径只是教学用的最简单路径；练习 3 中你可以构造交换机制的末态，比较两种机制的能垒。VASP 官方教程中也有 Pt 吸附原子在 fcc-Pt(001) 上"集体跳跃"的 NEB 例子可供对照。

---

## 13. 案例二：Li 在晶体 Si 中的扩散（体相 NEB）

**体系**：2×2×2 立方超胞 Si64 + 1 个 Li 间隙原子。Li 在金刚石 Si 中占据四面体间隙位（Td，4 个 Si 近邻，Li–Si ≈ 2.37 Å），沿 ⟨111⟩ 经六角间隙位（Hex）跳到相邻 Td 位，跳跃距离 = Si–Si 键长 ≈ 2.37 Å。

与表面 NEB 的不同点：

| 项目 | 表面吸附原子 | 体相间隙原子 |
|---|---|---|
| 固定原子 | 底部两层固定 | 通常全部放开（可固定一个远处原子防止整体漂移） |
| k 点 | z 方向 1 | 三个方向都要（2×2×2） |
| 晶胞 | 不变 | 不变（`ISIF = 2`，用纯 Si 的 PBE 晶格常数） |
| 偶极校正 | 需要 | 不需要 |
| 超胞收敛 | 表面胞大小 | 超胞大小（64 → 216 原子检查能垒） |

INCAR 要点（示例 [11-li-si-neb](../examples/11-li-si-neb/README.md)，Li 用 `Li_sv` POTCAR）：

```bash
ENCUT   = 500          # Li_sv has ENMAX ~ 500 eV
ISMEAR  = 0
SIGMA   = 0.05
IMAGES  = 3            # the middle image (02) starts at the hexagonal site
LCLIMB  = .TRUE.       # symmetric path with a good guess: climbing from the start is acceptable here
IBRION  = 3
POTIM   = 0
IOPT    = 7            # FIRE
EDIFFG  = -0.03
```

```bash
cd VASP/examples/11-li-si-neb
mkdir -p IS FS && cp POSCAR_IS IS/POSCAR && cp POSCAR_FS FS/POSCAR
for d in IS FS; do cp INCAR_relax $d/INCAR; cp KPOINTS $d/; done
POT_Li=Li_sv bash ../10-co-pt111-adsorption/make_potcar.sh IS FS && cp IS/POTCAR .
# ... 端点优化后：
bash ../11-pt-adatom-neb/make_images.sh IS FS 3
cp INCAR_neb INCAR && sbatch ../11-pt-adatom-neb/job.slurm
```

**参考量级**：文献中 DFT 计算的 Td→Hex→Td 能垒多在 0.5–0.65 eV 范围（近似，依赖超胞与泛函）；经典实验扩散激活能约 0.66 eV（近似值）。由于路径对称且初始猜测很好，本例直接打开了 CI；对于不对称或复杂路径仍应采用两阶段策略。

---

## 14. 案例三（思路）：N2 在金属表面的解离

N2 解离是氨合成的速控步，Ru(0001) 台阶位（B5 位）显著降低解离能垒（Dahl 等，[Phys. Rev. Lett. 1999](https://doi.org/10.1103/PhysRevLett.83.1814)）。搭建这类表面反应 NEB 的要点：

1. **初态**：N2 以端基（垂直）或侧基（平躺）吸附，比较两者，选能量更低且与末态"几何上相连"的构型；
2. **末态**：两个 N 原子分别位于相邻的 hollow 位（hcp/fcc），N–N 距离约 2.5–3 Å；要尝试多种末态（相邻 vs 次近邻 hollow），末态应是局部极小；
3. **插值**：N–N 键断裂伴随分子转动，**必须用 IDPP**，图像数取 5–7；
4. **中间极小**：常见"先从端基转为侧基吸附（前驱态），再解离"，若 NEB 路径上出现中间极小，应把它优化出来，拆成两段 NEB；
5. **自旋与偶极**：Ru 非磁，但含 Fe/Co/Ni 表面需要 `ISPIN = 2`；不对称 slab 打开偶极校正；
6. **验证**：TS 频率中唯一虚频对应 N–N 伸缩/分离方向。

---

## 15. 参数选择与收敛建议

| 参数 | 推荐 | 说明 |
|---|---|---|
| 端点 EDIFFG | −0.01～−0.02 | 比 NEB 更严 |
| NEB EDIFFG | 阶段一 −0.1～−0.05，阶段二 −0.03（精确工作 −0.02） | 投影力 |
| EDIFF | 1E-6（NEB），1E-7（Dimer、频率） | 力噪声是优化器振荡的常见原因 |
| IMAGES | 图像间距 ≲ 0.8–1 Å；对称路径用奇数 | 进程数为 IMAGES 整数倍 |
| SPRING | −5 | 一般不必改 |
| IOPT | 3/7 → 1/2 | 先稳健再快速 |
| LCLIMB | 阶段一 `.FALSE.`，阶段二 `.TRUE.` | VTST 默认 `.TRUE.` |
| 插值 | 单原子迁移：线性即可；多原子/分子转动：IDPP | 检查最小原子间距 |
| 超胞/表面胞 | 迁移原子与其周期像之间 ≥ 8–10 Å | 能垒变化 < 0.05 eV |
| k 点、ENCUT | 与端点完全相同；最终可用 Dimer 在更高精度下复核 | |
| 频率验证 | IBRION=5、POTIM=0.015、NFREE=2，放开反应原子及近邻 | 只有一个虚频 |

---

## 16. 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 作业一开始就退出，stdout 提示无法按图像划分进程（如 `M_divide: can not subdivide ...`） | MPI 进程数不是 IMAGES 的整数倍 | 进程数改为 IMAGES 的整数倍 |
| 图像"乱飞"、原子重叠、前几步受力几百 eV/Å | IS/FS 原子顺序不一致；跨周期边界的原子被线性插值穿过晶胞；线性插值不合理 | `idpp_interpolate.py --check-only` 检查顺序；用最小镜像插值或 IDPP；`nebmovie.pl 0` 看初始路径 |
| 设置了 IOPT 但原子不动 | VASP 不是 VTST 版本（`POTIM = 0` 时普通 VASP 不移动原子） | 用 VTST 编译的可执行文件 |
| 原子移动但不收敛、能量振荡 | 一开始就开 CI；LBFGS 步子太大；SCF 精度不足 | 两阶段策略；改 FIRE/QM、减小 `MAXMOVE`；`EDIFF = 1E-7` |
| 用了 `IBRION = 2` 后 NEB 收敛很差 | CG 线搜索依赖能量，与 NEB 投影力不一致 | VTST：`IBRION = 3, POTIM = 0` + `IOPT`；内置 NEB：`IBRION = 1` 或 3 |
| 以为跑的是普通 NEB，结果最高图像在爬坡 | VTST 中 `LCLIMB` 默认为 `.TRUE.` | 阶段一显式写 `LCLIMB = .FALSE.` |
| nebresults.pl 曲线端点能量为 0 或异常 | 00/ 和 NN/ 中没有 OUTCAR | 把 IS/FS 的 OUTCAR 复制进去 |
| 能垒曲线有两个极大值 | 路径上存在中间极小 | 优化中间极小，拆成两段 NEB |
| CI 图像不是能量最高的图像，或爬坡图像在两个图像间"跳" | 两个图像能量接近、未收敛 | 继续优化；加密图像 |
| 频率分析出现 2 个明显虚频 | 高阶鞍点，或 TS 未收敛 | 沿第二个虚频方向扰动后用 Dimer 重新收敛 |
| 频率分析没有虚频 | 爬坡图像未到鞍点，或放开原子太少 | 继续 CI-NEB/Dimer；放开更多原子 |
| Dimer 报错要求 IOPT | Dimer 只能用 VTST 优化器 | 设置 `IOPT = 2`（推荐）、`IBRION = 3`、`POTIM = 0` |
| Dimer 曲率始终为正、Torque 不下降 | 力不准或初始方向差 | `EDIFF = 1E-8`、`DdR = 0.01`、增大 `DRotMax`；用 neb2dim.pl 给出更好的 MODECAR |
| Dimer 收敛到一个与目标反应无关的鞍点 | 初始方向不对 | 从更接近 TS 的 NEB 图像出发，或检查 MODECAR 中反应原子的分量 |
| `ZBRENT: fatal error in bracketing` | 出现在端点优化中（IBRION=2） | 见[第18章](18-postprocess-troubleshooting.md)：续算、减小 POTIM、改 IBRION=1 或 FIRE |

---

## 练习

1. 用 1、3、5 个图像分别计算 Pt 吸附原子的跳跃能垒，观察能垒随图像数的变化；比较普通 NEB（只跑第一阶段并收敛到 −0.03）与 CI-NEB 的能垒差别。
2. 对爬坡图像做频率分析：①只放开吸附原子；②放开吸附原子 + 表面层。比较虚频大小，并用两种结果分别估算指前因子。
3. 构造交换机制的末态（吸附原子进入表面层、原来的表面原子被顶到相邻 hollow 位），用 IDPP + CI-NEB 计算能垒，与跳跃机制比较。
4. 把 Li-in-Si 的超胞换成 3×3×3（216 原子），检查能垒对超胞大小的收敛。
5. 在 Li-in-Si 的 NEB 结果上运行 `neb2dim.pl`，用 Dimer 在 ENCUT = 600 eV 下重新收敛鞍点，比较能垒变化。

---

## 参考资料

**官方文档**
- [VTST：Nudged Elastic Band](https://theory.cm.utexas.edu/vtsttools/neb.html)、[Optimizers](https://theory.cm.utexas.edu/vtsttools/optimizers.html)、[Dimer](https://theory.cm.utexas.edu/vtsttools/dimer.html)、[Dynamical Matrix](https://theory.cm.utexas.edu/vtsttools/dymmat.html)、[Scripts](https://theory.cm.utexas.edu/vtsttools/scripts.html)、[Installation](https://theory.cm.utexas.edu/vtsttools/installation.html)
- VASP Wiki：[Nudged elastic bands](https://www.vasp.at/wiki/index.php/Nudged_elastic_bands)、[IMAGES](https://www.vasp.at/wiki/index.php/IMAGES)、[SPRING](https://www.vasp.at/wiki/index.php/SPRING)、[IBRION](https://www.vasp.at/wiki/index.php/IBRION)
- [VASP 官方 Tutorials（含过渡态部分）](https://www.vasp.at/tutorials/latest/)
- [ASE NEB 文档](https://docs.ase-lib.org/ase/neb.html)

**相关 GitHub 仓库**（只链接，不转载代码）
- [henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode) —— VTST Fortran 补丁（CI-NEB、Dimer、Lanczos、IOPT 优化器），按 VASP 版本分目录
- [henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts) —— nebmake.pl、nebef.pl、nebresults.pl、neb2dim.pl、dymmatrix.pl 等脚本
- [QijingZheng/VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF) —— 把 OUTCAR 振动模式转成 XSF/动画，确认虚频方向
- [drinwater/Nudged-Elastic-Band-Tutorial](https://github.com/drinwater/Nudged-Elastic-Band-Tutorial) —— VTST CI-NEB 入门模板与 TS 频率验证示例
- [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script) —— 中文社区脚本集（含 NEB 插点等），仅供参考
- [northword/wiki](https://github.com/northword/wiki) —— 中文笔记，含 VTST 版 VASP 编译
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**
- G. Henkelman, B. P. Uberuaga, H. Jónsson, "A climbing image nudged elastic band method for finding saddle points and minimum energy paths", *J. Chem. Phys.* 113, 9901 (2000). [doi:10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
- G. Henkelman, H. Jónsson, "Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points", *J. Chem. Phys.* 113, 9978 (2000). [doi:10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
- G. Henkelman, H. Jónsson, "A dimer method for finding saddle points on high dimensional potential surfaces using only first derivatives", *J. Chem. Phys.* 111, 7010 (1999). [doi:10.1063/1.480097](https://doi.org/10.1063/1.480097)
- A. Heyden, A. T. Bell, F. J. Keil, "Efficient methods for finding transition states in chemical reactions: Comparison of improved dimer method and partitioned rational function optimization method", *J. Chem. Phys.* 123, 224101 (2005). [doi:10.1063/1.2104507](https://doi.org/10.1063/1.2104507)
- D. Sheppard, R. Terrell, G. Henkelman, "Optimization methods for finding minimum energy paths", *J. Chem. Phys.* 128, 134106 (2008). [doi:10.1063/1.2841941](https://doi.org/10.1063/1.2841941)
- S. Smidstrup, A. Pedersen, K. Stokbro, H. Jónsson, "Improved initial guess for minimum energy path calculations", *J. Chem. Phys.* 140, 214106 (2014). [doi:10.1063/1.4878664](https://doi.org/10.1063/1.4878664)
- G. H. Vineyard, "Frequency factors and isotope effects in solid state rate processes", *J. Phys. Chem. Solids* 3, 121 (1957). [doi:10.1016/0022-3697(57)90059-8](https://doi.org/10.1016/0022-3697(57)90059-8)
- G. L. Kellogg, P. J. Feibelman, "Surface self-diffusion on Pt(001) by an atomic exchange mechanism", *Phys. Rev. Lett.* 64, 3143 (1990). [doi:10.1103/PhysRevLett.64.3143](https://doi.org/10.1103/PhysRevLett.64.3143)
- S. Dahl et al., "Role of steps in N2 activation on Ru(0001)", *Phys. Rev. Lett.* 83, 1814 (1999). [doi:10.1103/PhysRevLett.83.1814](https://doi.org/10.1103/PhysRevLett.83.1814)

---

[← 上一章](10-adsorption-catalysis.md) | [返回目录](../README.md) | [下一章 →](12-phonon.md)
