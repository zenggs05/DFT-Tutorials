# 第4章 结构优化与状态方程（Structure Relaxation & Equation of State with VASP）

本章讲 VASP 结构优化（geometry optimization / structure relaxation）的完整流程：如何选 IBRION、ISIF、EDIFFG、POTIM，如何处理 Pulay 应力，二维材料怎样用 LATTICE_CONSTRAINTS 只优化面内晶格，分子和孤立体系怎样放进真空盒子用 Gamma 点计算。后半部分用 Birch–Murnaghan 状态方程（equation of state, E–V curve）拟合平衡体积 V0、体模量 B0 和 B0'，并给出可直接运行的 Si、MoS2、CO 算例和拟合脚本。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 根据体系的大小和离极小点的远近，选择 `IBRION = 1 / 2 / 3`（或 VTST 的 `IOPT` 优化器）；
- 用 `ISIF` 控制哪些自由度参与优化（原子位置、晶胞形状、晶胞体积），知道 `ISIF=3` 什么时候必须提高 `ENCUT`；
- 正确设置 `EDIFFG`、`NSW`、`EDIFF`，读懂 OSZICAR/OUTCAR，判断一次优化是否真的收敛（不只是"跑完了"）；
- 处理 Pulay 应力：变胞优化后用 CONTCAR 重新优化，直到体积不再变化；
- 用 `LATTICE_CONSTRAINTS`（VASP ≥ 6.4.3，推荐 ≥ 6.5.0）只优化二维材料的面内晶格常数；
- 在真空盒子中优化分子（`vasp_gam`、Gamma 点、`ISMEAR = 0`）；
- 做 E–V 扫描，用 Birch–Murnaghan 方程拟合得到 V0、B0、B0'。

| 场景 | 推荐做法 | 本章对应小节 |
|---|---|---|
| 体相晶体，晶格常数未知 | `ISIF=3` 变胞优化 + 重新优化消除 Pulay 应力；或 E–V 拟合 | §4、§8 |
| 已知晶格（如表面、吸附、缺陷超胞） | `ISIF=2` 只动原子 | §3.2 |
| 二维材料单层 | `ISIF=3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` | §6 |
| 分子/团簇 | 大盒子 + Gamma 点 + `ISIF=2` | §7 |
| 需要体模量 B0 | E–V 状态方程拟合 | §8 |
| 后续要算声子、弹性常数 | 更严格的 `EDIFFG`（≈ −1E-3 eV/Å）与 `EDIFF` | §9、[第12章](12-phonon.md)、[第13章](13-elastic.md) |

> [!NOTE]
> 本章默认你已经完成 [第2章](02-input-output-files.md) 的输入文件基础和 [第3章](03-convergence.md) 的 ENCUT / k 点收敛测试。结构优化的 ENCUT 和 k 点应当取自收敛测试的结果，而不是凭感觉填写。

---

## 2. 原理速览

### 2.1 结构优化在做什么

结构优化就是在势能面（potential energy surface）上找离初始结构最近的局部极小点。每一个离子步，VASP 先把电子自洽算到 `EDIFF` 精度，得到能量 $E$、每个原子的 Hellmann–Feynman 受力 $\mathbf{F}_I=-\partial E/\partial \mathbf{R}_I$，以及（`ISIF ≥ 2` 时）应力张量

$$
\sigma_{ij} = -\frac{1}{\Omega}\frac{\partial E}{\partial \varepsilon_{ij}},
$$

然后由优化算法更新原子坐标和晶格，直到受力（和应力）足够小。

> [!IMPORTANT]
> 优化算法只能找到**局部**极小。初始结构的对称性默认会被保留（`ISYM` 默认开启），所以如果你怀疑真实结构对称性更低（比如八面体倾转、Jahn–Teller 畸变），需要**主动给初始结构加一点小扰动**，而不是指望优化器自己破坏对称性。VASP 官方也建议优先扰动初始结构，而不是直接关掉对称性（`ISYM = 0`）。

### 2.2 Pulay 应力从哪里来

VASP 用固定的截断能 `ENCUT` 来确定平面波基组。平面波数目大约是

$$
N_{\mathrm{PW}} \propto \Omega \, E_{\mathrm{cut}}^{3/2},
$$

在一次变胞优化中，VASP **保持平面波数目（G 矢量集合）不变**，晶胞变化时等效截断能就随之改变：体积变小，等效截断能升高；形状改变，截断球变成椭球。基组不完备带来的这种人为应力称为 Pulay 应力，它通常使对角应力偏负，导致优化出来的体积**偏小**。

解决办法有三条，后面都会用到：

1. 变胞优化时把 `ENCUT` 提高到 ≥ 1.3 × max(ENMAX)，同时 `PREC = Accurate`；
2. 优化结束后把 CONTCAR 拷成 POSCAR **重新开始一次优化**（新的一次运行会按新晶胞重建基组），重复到体积不再变化；
3. 改用固定体积的一系列计算做 E–V 拟合（每个体积都重建基组，Pulay 误差在体积变化 5–10% 量级时基本被平均掉）。

> [!TIP]
> 只做"体积不变"的优化（`ISIF=2` 或 `ISIF=4`）时，Pulay 应力基本是各向同性的常数，对结构影响很小。真正需要警惕的是 `ISIF=3/6/7/8` 这类会改变体积的优化。

### 2.3 状态方程

固定若干体积 $V$ 算出能量 $E(V)$，用三阶 Birch–Murnaghan 方程拟合：

$$
E(V) = E_0 + \frac{9V_0B_0}{16}\left\{\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^3 B_0' + \left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]^2\left[6-4\left(\frac{V_0}{V}\right)^{2/3}\right]\right\}
$$

其中 $V_0$ 为平衡体积，$B_0 = V\,\partial^2E/\partial V^2\big|_{V_0}$ 为体模量，$B_0' = \partial B/\partial P\big|_{P=0}$。单位换算：$1\ \mathrm{eV/\AA^3} = 160.2177\ \mathrm{GPa}$。

一个有用的事实：三阶 BM 方程恰好是 $x = V^{-2/3}$ 的三次多项式，所以用线性最小二乘就能精确拟合，不需要非线性优化的初值。本章示例脚本 `fit_bm.py` 就采用这种做法。

---

## 3. 结构优化的关键标签

### 3.1 IBRION：选择优化算法

| IBRION | 算法 | 特点 | 适用 |
|---|---|---|---|
| 2 | 共轭梯度（CG） | 每个搜索方向做线搜索（trial step + corrector step），稳健；对 `POTIM` 不太敏感 | **默认首选**，初始结构离极小较远、不确定怎么设 POTIM 时 |
| 1 | RMM-DIIS 准牛顿 | 用历史步构造近似 Hessian，接近极小时很快；只看力和应力，不看能量；要求受力精确 | 已接近极小、自由度较多（>20）的体系；第二阶段精修 |
| 3 | 阻尼分子动力学（damped MD） | 由 `POTIM`（时间步）与 `SMASS`（阻尼因子 μ）控制；`SMASS < 0` 时用速度淬火（quick-min） | 离极小很远的大体系做"粗优化"；配合 VTST 优化器 |
| 3 + `POTIM = 0` + `IOPT` | VTST 力基优化器 | 需要 VTST 补丁编译的 VASP | NEB/Dimer 必需；普通优化也可用 FIRE、LBFGS |

几点经验：

- **IBRION = 2** 的 stdout 中每个 corrector 步会打印 `trialstep`，把当前 `POTIM` 乘以它就是更合适的步长（VASP 官方建议用 CG 先找出合理的 POTIM）。
- **IBRION = 1** 依赖精确受力，建议同时设 `NELMIN = 4`（体相）到 `8`（表面），保证每个离子步至少做几次电子迭代；可用 `NFREE`（如 10）限制历史步数，低维体系 `NFREE` 不宜太大。
- **IBRION = 3** 的 `SMASS` 在这里是阻尼因子 μ：μ = 2 等价于不带线搜索的最速下降，μ = 0 无阻尼，官方建议的初值约 0.4；时间步 `POTIM` 必须给（太大发散、太小很慢）。

**VTST 优化器**（需要按 [第1章](01-setup-compile.md) 打 VTST 补丁）：设 `IBRION = 3`、`POTIM = 0` 关闭 VASP 自带的离子更新，再用 `IOPT` 选择算法：

| IOPT | 算法 | 说明 |
|---|---|---|
| 0 | 使用 IBRION 指定的 VASP 自带算法 | 默认 |
| 1 | LBFGS | 受力精确时效率高，`MAXMOVE` 默认 0.2 Å |
| 2 | CG（基于力的牛顿线搜索） | 需要精确受力 |
| 3 | Quick-Min | 最稳妥，VTST 推荐新手使用 |
| 7 | FIRE | 类似 quick-min，通常更快，适合受力大或噪声大的情况 |

> [!NOTE]
> VTST 优化器是"纯力"算法（不使用能量），因此能用于 NEB 等能量与力不一致的情形；对普通体相优化，它只移动原子，变胞优化仍建议用 VASP 自带的 `IBRION = 1/2`。详细用法见 [第11章](11-neb-dimer.md)。

### 3.2 ISIF：哪些自由度参与优化

`ISIF` 同时决定"是否计算应力"和"哪些自由度可变"（下表依据 VASP Wiki 的 ISIF 页面整理）：

| ISIF | 计算受力 | 计算应力 | 原子位置 | 晶胞形状 | 晶胞体积 | 典型用途 |
|---|---|---|---|---|---|---|
| 0 | 是 | 否 | 是 | 否 | 否 | MD（`IBRION=0`）和杂化泛函的默认值 |
| 1 | 是 | 仅迹（压强） | 是 | 否 | 否 | 很少用 |
| 2 | 是 | 是 | 是 | 否 | 否 | 表面、吸附、缺陷等固定晶胞的优化（非 MD、非杂化时的默认值） |
| 3 | 是 | 是 | 是 | 是 | 是 | 体相完全优化 |
| 4 | 是 | 是 | 是 | 是 | 否 | 固定体积优化（E–V 扫描中的非立方体系） |
| 5 | 是 | 是 | 否 | 是 | 否 | 只变形状 |
| 6 | 是 | 是 | 否 | 是 | 是 | 只变晶胞 |
| 7 | 是 | 是 | 否 | 否 | 是 | 只变体积（各向同性缩放） |
| 8 | 是 | 是 | 是 | 否 | 是 | 原子 + 体积，形状不变（**VASP ≥ 6.4.1**） |

> [!WARNING]
> 杂化泛函（`LHFCALC = .TRUE.`）下 `ISIF` 的默认值是 0，**不计算应力**。如果你要用 HSE06 做变胞优化，必须显式写 `ISIF = 3`，否则晶格根本不会动。详见 [第6章](06-hybrid-metagga.md)。

`POSCAR` 中的 `Selective dynamics` 可以进一步固定个别原子的个别方向（例如 slab 底部几层），见 [第2章](02-input-output-files.md) 和 [第9章](09-surface-workfunction.md)。注意 Selective dynamics 的 T/F 标志是沿**晶格矢量（分数坐标）方向**的，不是笛卡尔方向。

### 3.3 收敛判据：EDIFFG、EDIFF、NSW

| 标签 | 含义 | 推荐值 | 说明 |
|---|---|---|---|
| `EDIFFG` | 离子步收敛判据 | `-0.01` ~ `-0.02` eV/Å（常规）；`-1E-3`（声子/弹性前的预优化） | **负值**：所有原子受力的模都小于 \|EDIFFG\| 时停止；正值：相邻两个离子步能量差小于 EDIFFG 时停止；默认是 `EDIFF × 10`（正值，按能量判断，偏松） |
| `EDIFF` | 电子自洽收敛判据 | `1E-6`（优化）；`1E-7`~`1E-8`（声子/弹性） | 受力精度受它限制。`EDIFF` 太松时 CG 线搜索容易失败（ZBRENT 报错） |
| `NSW` | 最大离子步数 | 50–200 | 到达 NSW 仍未收敛，VASP 会正常结束但**不会**打印 "reached required accuracy" |
| `NELMIN` | 每个离子步最少电子步数 | 4–8 | 对 `IBRION=1` 特别有用 |
| `POTIM` | 步长缩放（IBRION=1/2/3） | 默认 0.5 | `IBRION=1/3` 对它敏感；`IBRION=2` 会自动调整 |

> [!NOTE]
> 对于 `ISIF ≥ 3`，VASP 开发者在论坛中说明负的 `EDIFFG` 也会同时检查应力分量（以内部单位比较，对用户并不直观）。因此变胞优化结束后，请**自己**到 OUTCAR 中确认 `external pressure` 和各应力分量足够小（通常 |P| < 1–2 kB 即可，具体取决于体系硬度），不要只看 "reached required accuracy"。

### 3.4 如何选择：一个实用流程

```text
初始结构离极小很远（手搭结构、替换了元素、受力 > 1 eV/A）
  └─> 第一阶段：IBRION=2（或 IBRION=3 + SMASS 阻尼 / VTST FIRE），PREC=Normal，EDIFFG=-0.05
        └─> 第二阶段：cp CONTCAR POSCAR，IBRION=1 或 2，PREC=Accurate，EDIFFG=-0.01
              └─> 若 ISIF=3：再 cp CONTCAR POSCAR 重新优化，直到体积和能量不再变化
```

---

## 4. 实战一：Si 体相变胞优化（ISIF = 3）与 Pulay 应力处理

完整输入在 [../examples/04-relax-si-isif3/README.md](../examples/04-relax-si-isif3/README.md)。

### 4.1 输入文件

故意把初始晶格常数设为 5.55 Å（比实验值 5.431 Å 大约 2%），观察优化过程。

POSCAR（fcc 原胞，金刚石结构）：

```text
Si diamond primitive cell, start from a = 5.55 A
5.55
  0.0  0.5  0.5
  0.5  0.0  0.5
  0.5  0.5  0.0
Si
2
Direct
  0.00  0.00  0.00
  0.25  0.25  0.25
```

KPOINTS（fcc 晶格推荐 Γ 中心网格）：

```text
Gamma-centered 8x8x8
0
Gamma
8 8 8
0 0 0
```

INCAR：

```text
SYSTEM  = Si bulk relaxation (ISIF=3)

# Electronic part
PREC    = Accurate
ENCUT   = 400        # >= 1.3 x ENMAX(Si ~245 eV); required for stress/volume relaxations
EDIFF   = 1E-6
ISMEAR  = 0          # Gaussian smearing, semiconductor
SIGMA   = 0.05
LREAL   = .FALSE.    # small cell: reciprocal-space projection is more accurate

# Ionic part
IBRION  = 2          # conjugate gradient, robust default
ISIF    = 3          # relax positions + cell shape + cell volume
NSW     = 60
EDIFFG  = -0.005     # force criterion (stress is also checked for ISIF>=3)

# Output
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

POTCAR（仓库中永远不提供 POTCAR，请从你自己的赝势库拼接）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR     # 按本地目录名调整，如 potpaw_PBE.64
grep -E "TITEL|ENMAX" POTCAR                         # 确认元素与 ENMAX
```

关键词逐一说明：

| 关键词 | 取值 | 为什么 |
|---|---|---|
| `PREC` | Accurate | FFT 网格足够精确地计算 $H|\psi\rangle$，减少应力噪声和 E–V 曲线的锯齿 |
| `ENCUT` | 400 | Si 的 ENMAX ≈ 245 eV，变胞优化要 ≥ 1.3×ENMAX ≈ 320 eV，再留一些余量 |
| `ISMEAR`/`SIGMA` | 0 / 0.05 | 半导体用 Gaussian 展宽；金属可改为 `ISMEAR=1`、`SIGMA=0.1~0.2` |
| `IBRION` | 2 | CG 最稳健；Si 只有晶格一个自由度，任何算法都能很快收敛 |
| `ISIF` | 3 | 位置 + 形状 + 体积全部放开 |
| `EDIFFG` | −0.005 | 金刚石 Si 原子位置被对称性固定（受力恒为 0），收敛实际由应力决定 |

### 4.2 运行

```bash
cd VASP/examples/04-relax-si-isif3
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
sbatch job.slurm                    # 或交互测试：mpirun -np 4 vasp_std > vasp.out
```

作业脚本的写法见 [Common 第2章](../../Common/tutorials/02-job-scheduler.md) 和本仓库的 [`job_vasp_slurm.sh`](../scripts/job_vasp_slurm.sh)（用 `bash ../scripts/job_vasp_slurm.sh -h` 或打开文件查看说明）。

### 4.3 监控与检查

```bash
# 每个离子步一行：F 为自由能 TOTEN，E0 为 sigma->0 外推能量
grep "F=" OSZICAR

# 是否正常收敛（没有这一行说明 NSW 用完了或中途出错）
grep "reached required accuracy" OUTCAR

# 最后一个离子步的压强与应力（单位 kB）
grep "external pressure" OUTCAR | tail -1
grep "in kB" OUTCAR | tail -1

# 体积随离子步的变化
grep "volume of cell" OUTCAR

# 优化后的晶格矢量
head -5 CONTCAR
```

OSZICAR 的典型形式（示意输出，数值仅供参考）：

```text
       N       E                     dE             d eps       ncg     rms          rms(c)
DAV:   1    -0.108094E+02   -0.10809E+02   -0.48322E+02   864   0.126E+02    0.110E+01
...
   1 F= -.10790000E+02 E0= -.10790000E+02  d E =-.107900E+02
...
   4 F= -.10842000E+02 E0= -.10842000E+02  d E =-.520000E-01
```

OUTCAR 中压强行的典型形式（示意输出，数值仅供参考）：

```text
  external pressure =        0.42 kB  Pullay stress =        0.00 kB
```

> [!NOTE]
> 在 VASP 的约定下，`external pressure` 为**正**表示晶胞想**膨胀**（处于压缩状态），为**负**表示晶胞想收缩。"Pullay stress" 一栏显示的是你在 INCAR 中通过 `PSTRESS` 施加的值，而不是 VASP 自动估算的 Pulay 误差。

### 4.4 消除 Pulay 应力：重新优化

```bash
mkdir -p run1 && cp INCAR KPOINTS POTCAR POSCAR OUTCAR OSZICAR CONTCAR run1/
cp CONTCAR POSCAR
sbatch job.slurm                    # 第二次优化，新运行会按新晶胞重建平面波基组
```

判断标准：

- 第二次运行在 1–3 个离子步内结束；
- 两次的体积差 < 0.1–0.2%，能量差 < 1 meV/atom；
- 最后的 `external pressure` 的绝对值在 1 kB 左右以内。

> [!WARNING]
> 续算时不要设置 `ISTART = 2`。`ISTART = 2` 表示"保持基组不变"的精确续算，会把 Pulay 应力也原样带过去。默认的 `ISTART = 1`（存在 WAVECAR 时）或干脆不读 WAVECAR 都会在新晶胞上重建基组。

本仓库示例的 `job.slurm` 已经把"优化→拷贝 CONTCAR→再优化"写成循环，直到某一轮离子步数 ≤ 3 为止。

### 4.5 结果怎么看

- PBE 对 Si 晶格常数通常给出约 5.46–5.47 Å（略高于实验值 5.431 Å，GGA 普遍高估约 1%）。系统性的基准可见 Haas, Tran & Blaha, *PRB* 79, 085104 (2009)。你得到的值取决于 ENCUT、k 点和 POTCAR 版本，偏差 0.005 Å 量级属正常。
- 从 CONTCAR 读晶格常数：fcc 原胞的晶格矢量长度 $|\mathbf a_1| = a/\sqrt{2}$，所以 $a = \sqrt{2}\,|\mathbf a_1|$。

```python
# Python one-off check (needs ase)
from ase.io import read
atoms = read("CONTCAR")
a1 = atoms.cell.lengths()[0]
print("conventional a = %.4f A" % (a1 * 2 ** 0.5))
print("volume per atom = %.4f A^3" % (atoms.get_volume() / len(atoms)))
```

---

## 5. 输出检查：怎样才算"真的收敛了"

一次合格的结构优化至少满足以下全部条件：

1. OUTCAR 中出现 `reached required accuracy - stopping structural energy minimisation`；
2. 最后几个离子步的电子步数没有达到 `NELM`（否则说明电子没收敛，受力不可信）；
3. 自由原子上的最大受力 < |EDIFFG|（要排除 Selective dynamics 固定的原子）；
4. 对 `ISIF=3`，压强/应力足够小，并且重新优化后体积不再变化；
5. 结构合理：键长、配位没有异常（用 VESTA 或 ASE 看一眼 CONTCAR）。

### 5.1 每个离子步用了多少电子步

```bash
# 统计 OSZICAR 中每个离子步之前的电子步数（DAV/RMM/CG 开头的行）
awk '/^[A-Z]+ *:/{n++} / F= /{printf "ionic step %3d : %3d electronic steps\n", $1, n; n=0}' OSZICAR
grep NELM INCAR OUTCAR | head -3          # 对照 NELM 上限（默认 60）
```

如果很多离子步的电子步数都等于 NELM，请先解决电子收敛问题（`ALGO`、`AMIX`、初始磁矩等，见 [第18章](18-postprocess-troubleshooting.md)），再谈结构优化。

### 5.2 自由原子上的最大受力

OUTCAR 里每个离子步都有一个 `TOTAL-FORCE` 块。没有固定原子时，最大受力可以直接这样求：

```bash
# 最后一个 TOTAL-FORCE 块中最大的 |F|（不考虑 Selective dynamics）
NIONS=$(grep -m1 "NIONS" OUTCAR | awk '{print $NF}')
grep -A $((NIONS+1)) "TOTAL-FORCE" OUTCAR | tail -n $NIONS | \
  awk '{f=sqrt($4^2+$5^2+$6^2); if(f>m){m=f; i=NR}} END{printf "max |F| = %.4f eV/A (atom %d)\n", m, i}'
```

有 Selective dynamics 时，需要把固定方向的受力排除。下面这个小脚本读取 POSCAR 的 T/F 标志和 OUTCAR 最后一个受力块：

```python
#!/usr/bin/env python3
"""Max force on free degrees of freedom (last ionic step). Usage: python maxforce.py [POSCAR] [OUTCAR]"""
import sys
import numpy as np

def free_mask(poscar):
    lines = open(poscar).read().splitlines()
    counts = [int(x) for x in lines[6].split()]          # VASP5+ format: line 7 = atom counts
    nat = sum(counts)
    if not lines[7].strip()[:1] in ("s", "S"):           # no Selective dynamics
        return np.ones((nat, 3), dtype=bool)
    rows = [l.split()[3:6] for l in lines[9:9 + nat]]     # skip 'Selective' + 'Direct/Cartesian'
    return np.array([[t.upper().startswith("T") for t in r] for r in rows])

def last_forces(outcar, nat):
    lines = open(outcar).read().splitlines()
    k = max(i for i, l in enumerate(lines) if "TOTAL-FORCE" in l)
    return np.array([[float(x) for x in l.split()[3:6]] for l in lines[k + 2:k + 2 + nat]])

poscar = sys.argv[1] if len(sys.argv) > 1 else "POSCAR"
outcar = sys.argv[2] if len(sys.argv) > 2 else "OUTCAR"
mask = free_mask(poscar)
f = np.where(mask, last_forces(outcar, len(mask)), 0.0)
norm = np.linalg.norm(f, axis=1)
i = int(np.argmax(norm))
print("max |F| on free dof = %.4f eV/A (atom %d)" % (norm[i], i + 1))
```

> [!NOTE]
> T/F 标志作用在分数坐标方向上，脚本把它直接套到笛卡尔受力分量上，这对"c 轴垂直于 ab 面、只固定整个原子或只固定 z"的常见 slab 是精确的，对斜交晶胞的部分固定只是近似。本仓库的 [`check_vasp.py`](../scripts/check_vasp.py) 会把电子步、离子步、是否收敛、最大受力、常见报错一起检查，用法见 `python ../scripts/check_vasp.py -h`。

### 5.3 用 Python 快速看优化轨迹

```python
from ase.io import read
traj = read("OUTCAR", index=":")          # every ionic step
for i, a in enumerate(traj):
    fmax = (a.get_forces() ** 2).sum(axis=1).max() ** 0.5
    print("%3d  E = %.6f eV  V = %.3f A^3  fmax = %.4f eV/A" % (i + 1, a.get_potential_energy(), a.get_volume(), fmax))
```

更多 ASE/pymatgen 解析方法见 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。

---

## 6. 实战二：二维 MoS2 只优化面内晶格（LATTICE_CONSTRAINTS）

完整输入在 [../examples/04-relax-mos2-2d/README.md](../examples/04-relax-mos2-2d/README.md)。

### 6.1 问题：为什么不能对 slab 直接 ISIF=3

单层材料的 POSCAR 在 c 方向有十几埃的真空。直接 `ISIF = 3` 时，真空区几乎没有应力"撑着"，c 轴会被压缩（或在数值噪声下乱动），真空层可能塌缩到层间相互作用的距离，结果毫无意义。我们真正想要的是：**a、b 放开，c 固定**。

### 6.2 LATTICE_CONSTRAINTS 的用法与版本要求

```text
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.     # a1, a2 free; a3 (vacuum direction) fixed
```

- 三个逻辑值分别对应晶格矢量 a1、a2、a3。对正交晶胞，`.FALSE.` 的方向晶格常数不变。
- 对非正交晶胞（如六角），它的作用是把应力张量对应的**行和列**置零，从而相应的晶格矩阵分量不更新。对 "a1、a2 在 xy 平面内、a3 沿 z" 的二维材料，`T T F` 恰好保证 a3 不变、a1/a2 只在面内变化。
- **版本**：与 `IBRION = 1/2` 一起使用从 **VASP 6.4.3** 起才生效（6.4.3 之前该标签在 IBRION=1/2 下根本不会被读取，等于没设）；VASP Wiki 还提示 **< 6.5.0** 的版本存在相关 bug（6.5.0 修复了"约束在检查收敛判据之前未被施加"的问题）；与 `ISIF = 4/5` 组合则需要 **≥ 6.5.1**（known issue 34）。**推荐使用 6.5.1 及以上版本**，并在 OUTCAR 中核对晶格是否真的按预期变化。

> [!WARNING]
> 如果你的集群只有 6.4.2 或更早的版本，写上 `LATTICE_CONSTRAINTS` 不会报错，但会被**静默忽略**，c 轴照样会动。旧版本的替代方案见 §6.5。

### 6.3 输入文件

POSCAR（1H-MoS2，120° 六角晶胞，真空方向 c = 20 Å，初始 a = 3.20 Å 略大于实验值）：

```text
MoS2 monolayer 1H, a = 3.20 A start, c = 20 A
1.0
   3.2000000000   0.0000000000   0.0000000000
  -1.6000000000   2.7712812921   0.0000000000
   0.0000000000   0.0000000000  20.0000000000
Mo S
1 2
Direct
  0.3333333333  0.6666666667  0.5000000000
  0.6666666667  0.3333333333  0.5792500000
  0.6666666667  0.3333333333  0.4207500000
```

KPOINTS（二维体系：Γ 中心，z 方向取 1）：

```text
Gamma-centered 12x12x1 for 2D
0
Gamma
12 12 1
0 0 0
```

INCAR：

```text
SYSTEM  = MoS2 monolayer, in-plane relaxation only

PREC    = Accurate
ENCUT   = 450        # >= 1.3 x max(ENMAX) of Mo_pv/S
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.05
LREAL   = .FALSE.
LASPH   = .TRUE.     # recommended for d elements

IBRION  = 2
ISIF    = 3
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.   # VASP >= 6.4.3 (>= 6.5.1 recommended)
NSW     = 100
EDIFFG  = -0.01

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

POTCAR：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
grep TITEL POTCAR     # order must be Mo then S, same as POSCAR line 6
```

### 6.4 检查结果

```bash
grep "reached required accuracy" OUTCAR
head -5 CONTCAR                      # a1, a2 changed; a3 must still be 0 0 20
grep "in kB" OUTCAR | tail -1        # XX, YY, XY should be ~0; ZZ is not relaxed and can stay nonzero
```

```python
from ase.io import read
a = read("CONTCAR")
print("a = %.4f A, c = %.4f A" % (a.cell.lengths()[0], a.cell.lengths()[2]))
z = sorted(a.positions[:, 2])
print("S-S vertical distance = %.4f A" % (z[-1] - z[0]))
```

参考量级：单层 MoS2 实验晶格常数约 3.16 Å（与体相 2H-MoS2 接近），PBE 一般给出约 3.18–3.19 Å（近似值，随 POTCAR 和 ENCUT 略有不同）。

> [!TIP]
> 单层材料的面内晶格基本不受层间 vdW 影响，可以不加 vdW 修正；但如果后面要和双层、异质结比较能量，请从一开始就用同一套 vdW 设置，见 [第8章](08-vdw-2d.md)。

### 6.5 旧版本 VASP 的替代方案

1. **二维 E–a 扫描**：固定 c 和真空，缩放面内晶格常数 a（例如 ±4%，7 个点），每个点 `ISIF = 2` 只优化原子，再对 E(a) 做多项式拟合求最小值。思路与 §8 相同，只是把"体积"换成"面积"。
2. **重新编译打补丁**：社区补丁 [Chengcheng-Xiao/VASP_OPT_AXIS](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS) 可以在旧版 VASP 中固定晶格分量（仅作链接，不在本仓库转载）。6.4.3 之后建议直接用 `LATTICE_CONSTRAINTS`。

---

## 7. 实战三：分子与孤立体系（CO 放在真空盒子中）

完整输入在 [../examples/04-relax-co-molecule/README.md](../examples/04-relax-co-molecule/README.md)。

### 7.1 要点

| 要点 | 做法 | 原因 |
|---|---|---|
| 盒子大小 | 分子到相邻镜像的最近距离 ≥ 10 Å（小分子常用 12–15 Å 的盒子） | 平面波方法天然周期，需要真空隔开镜像 |
| k 点 | 只用 Γ 点：`Gamma 1 1 1` | 孤立体系没有色散 |
| 可执行程序 | `vasp_gam` | 利用 Γ 点波函数为实数，速度约快 1 倍、内存减半 |
| 展宽 | `ISMEAR = 0`，`SIGMA = 0.01~0.05` | 分子能级离散，**不能**用 `ISMEAR = -5`（四面体法要求 ≥ 4 个 k 点） |
| 自旋 | 闭壳层（CO、H2O、CH4）`ISPIN = 1`；开壳层（O2 三重态、自由基）`ISPIN = 2` 并给 `MAGMOM` 或 `NUPDOWN` | 自旋态决定能量，吸附能计算中尤其重要 |
| ISIF | 2 | 盒子是人为的，不能优化晶胞 |
| 极性分子 | 可加偶极校正（`LDIPOL`、`IDIPOL = 4`） | 见 [第9章](09-surface-workfunction.md) 和 [第10章](10-adsorption-catalysis.md) |

### 7.2 输入文件

POSCAR（笛卡尔坐标，CO 沿 z 放在盒子中心附近，初始键长 1.16 Å，略长于实验值 1.128 Å）：

```text
CO molecule in a 15 A cubic box
1.0
  15.0  0.0  0.0
   0.0 15.0  0.0
   0.0  0.0 15.0
C O
1 1
Cartesian
  7.5  7.5  6.92
  7.5  7.5  8.08
```

KPOINTS：

```text
Gamma point only
0
Gamma
1 1 1
0 0 0
```

INCAR：

```text
SYSTEM  = CO molecule in a box

PREC    = Accurate
ENCUT   = 520
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.01
ISPIN   = 1          # closed shell; use ISPIN=2 (+MAGMOM) for O2, radicals
LREAL   = .FALSE.

IBRION  = 2
ISIF    = 2          # never relax the box of an isolated molecule
NSW     = 50
EDIFFG  = -0.01

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

运行（注意用 `vasp_gam`）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{C,O}/POTCAR > POTCAR
mpirun -np 8 vasp_gam > vasp.out       # 或 sbatch job.slurm（脚本里已写 vasp_gam）
grep "reached required accuracy" OUTCAR
python -c "from ase.io import read; a=read('CONTCAR'); print('d(C-O) = %.4f A' % a.get_distance(0, 1, mic=True))"
```

参考量级：CO 实验键长 1.128 Å，PBE 一般给出约 1.13–1.14 Å（近似值）。

> [!WARNING]
> `vasp_gam` 只能用于 KPOINTS 中**只有 Γ 点**的计算；对体相用 `vasp_gam` 会直接报错。反过来，用 `vasp_std` 算 Γ 点也能得到正确结果，只是慢一些。

### 7.3 盒子大小收敛

对要比较能量的分子（如吸附能中的气相分子），建议做一次盒子测试：12、15、18 Å 三个盒子，能量差 < 1–2 meV 即可。极性分子（如 H2O、CO）受镜像偶极作用影响较大，盒子要更大或加偶极校正。

---

## 8. 状态方程（E–V 曲线）拟合

完整脚本在 [../examples/04-eos-si/README.md](../examples/04-eos-si/README.md)。

### 8.1 什么时候用 E–V 拟合，而不是 ISIF=3

- 需要体模量 B0、B0'（例如与实验对比、准谐近似的输入）；
- 想在**默认截断能**下得到可靠的平衡体积（每个体积都重建基组，Pulay 误差被平均掉；官方建议体积变化取 5–10% 量级）；
- 比较不同相（如 fcc/bcc/hcp、相变压力）的能量—体积关系；
- `ISIF=3` 收敛困难或对称性复杂时作为交叉验证。

### 8.2 计算设置

| 体系 | 每个体积点做什么 | INCAR 要点 |
|---|---|---|
| 立方且原子位置被对称性固定（Si、Cu、NaCl） | 单点能 | `NSW = 0`、`IBRION = -1`；半导体可用 `ISMEAR = -5` 得到更精确的总能 |
| 立方但有内坐标（如钙钛矿畸变相） | 固定体积优化原子 | `ISIF = 2`、`IBRION = 2` |
| 非立方（六角、四方、正交……） | 固定体积优化原子 + 形状 | `ISIF = 4`、`IBRION = 2`，并提高 `ENCUT` |

体积范围：一般取 $V/V_{\mathrm{ref}} = 0.94 \sim 1.06$（7 个点），或 $0.90 \sim 1.10$（9 个点），参考体积选在估计的极小附近，让极小点落在扫描区间中部。

### 8.3 示例：Si 的 E–V 扫描

示例目录中的 `run_eos.sh` 分三步完成：

```bash
cd VASP/examples/04-eos-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

bash run_eos.sh prep        # 生成 V_0.94 ... V_1.06 七个目录（缩放 POSCAR 第 2 行）
bash run_eos.sh submit      # 每个目录 sbatch ../job.slurm
bash run_eos.sh collect     # 读取每个 OUTCAR 的体积和 energy(sigma->0)，写出 eos.dat
python fit_bm.py eos.dat --lattice fcc --plot
```

`run_eos.sh` 中缩放 POSCAR 的核心逻辑是：体积缩放 $s$ 对应长度缩放 $s^{1/3}$，所以把第 2 行的缩放因子乘以 $s^{1/3}$（如果缩放因子写成负数，它表示体积，此时直接乘以 $s$）。

单点计算的 INCAR：

```text
SYSTEM  = Si E-V scan (fixed volume, static)
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
ISMEAR  = -5         # tetrahedron with Bloechl corrections: accurate total energy for a semiconductor
LREAL   = .FALSE.
NSW     = 0
IBRION  = -1
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

`eos.dat` 的格式是两列：体积（Å³/晶胞）与能量（eV/晶胞）。能量取 OUTCAR 中最后一个 `energy(sigma->0)`：

```bash
for d in V_*; do
  V=$(grep "volume of cell" $d/OUTCAR | tail -1 | awk '{print $NF}')
  E=$(grep "energy(sigma->0)" $d/OUTCAR | tail -1 | awk '{print $NF}')
  echo "$V $E"
done | sort -n > eos.dat
```

### 8.4 拟合：三种方法

**方法 1：本章示例脚本 `fit_bm.py`**（只依赖 numpy，matplotlib 可选），输出 V0、E0、B0（GPa）、B0'，`--lattice fcc` 时还会换算出立方晶格常数 $a_0 = (4V_0)^{1/3}$。

```bash
python fit_bm.py eos.dat --lattice fcc --plot
```

**方法 2：ASE**

```python
import numpy as np
from ase.eos import EquationOfState
from ase.units import kJ

V, E = np.loadtxt("eos.dat", unpack=True)
eos = EquationOfState(V, E, eos="birchmurnaghan")
v0, e0, B = eos.fit()
print("V0 = %.4f A^3, E0 = %.6f eV, B0 = %.1f GPa" % (v0, e0, B / kJ * 1.0e24))
eos.plot("eos_ase.png")
```

**方法 3：pymatgen**

```python
import numpy as np
from pymatgen.analysis.eos import EOS

V, E = np.loadtxt("eos.dat", unpack=True)
fit = EOS(eos_name="birch_murnaghan").fit(V, E)
print(fit.v0, fit.e0, fit.b0_GPa, fit.b1)     # b1 is B0'
```

### 8.5 怎样判断拟合可信

- 极小点 $V_0$ 必须落在扫描区间**内部**，最好在中间三分之一；否则平移参考体积重扫；
- 拟合残差（RMS）应在 meV 量级以下；某个点明显偏离曲线，多半是该点电子没收敛或 SCF 落到了别的磁态；
- 用 BM、Murnaghan、Vinet 三种形式拟合，V0 差别应 < 0.1%，B0 差别在几个 GPa 内；
- 与 `ISIF=3`（重新优化后）得到的体积对比，差别应 < 0.5%。

参考量级：Si 实验体模量约 98–99 GPa；PBE 一般给出约 88–90 GPa（GGA 高估体积、低估 B0 是普遍趋势），具体数值随设置不同。

> [!NOTE]
> E–V 拟合得到的 $B_0$ 是等温体模量在 0 K、不含零点振动的近似。要和室温实验严格比较，需要准谐近似（QHA）等更高级的处理，见 [第12章](12-phonon.md)。

---

## 9. 进阶话题

### 9.1 外加压强：PSTRESS

`PSTRESS`（单位 kB）在优化时给能量加上 $PV$ 项，并从应力对角元中减去该值，于是可以在指定压强下做结构优化。例如 10 GPa = 100 kB：

```text
ISIF    = 3
PSTRESS = 100        # 100 kB = 10 GPa external pressure
```

此时比较不同相的稳定性要看焓 $H = E + PV$（OUTCAR 中会额外给出包含 $PV$ 项的焓）。`PSTRESS` 也可以用来近似修正 Pulay 应力（取低/高截断能下 `external pressure` 之差），但 VASP 官方强调这只是 ENCUT 无法提高时的权宜之计。

### 9.2 声子、弹性计算前的"紧"优化

有限位移声子和弹性常数对残余受力非常敏感。建议：

```text
PREC    = Accurate
EDIFF   = 1E-8
EDIFFG  = -1E-3
ADDGRID = .TRUE.     # optional, smoother forces
LREAL   = .FALSE.
```

并且用与后续声子/弹性计算**完全相同**的 ENCUT、k 点和泛函。见 [第12章](12-phonon.md)、[第13章](13-elastic.md)。

### 9.3 磁性、DFT+U、vdW、杂化泛函下的优化

- 磁性体系：`ISPIN = 2` 并给合理的 `MAGMOM`，否则可能优化到非磁态；不同磁构型分别优化，见 [第7章](07-magnetism-u-soc.md)。
- 层状材料：不加 vdW 修正时层间距会严重高估，见 [第8章](08-vdw-2d.md)。
- 杂化泛函优化：显式写 `ISIF`，计算量大，常见做法是 PBE 结构上做 HSE 单点（记为 HSE06@PBE），见 [第6章](06-hybrid-metagga.md)。
- meta-GGA（SCAN/r2SCAN）+ 自旋极化：**VASP 6.6.0** 的受力有误（known issue 94），请用 6.6.1 或打官方补丁后再做优化/MD。

### 9.4 续算与断点

VASP 没有"断点续算"一说，优化中断后的标准做法是：

```bash
cp CONTCAR POSCAR          # CONTCAR is updated after every ionic step
# optionally keep WAVECAR/CHGCAR to speed up the first SCF
sbatch job.slurm
```

如果中断时 CONTCAR 是空的（作业在写文件时被杀），就用上一份备份或 XDATCAR 最后一帧（`ase.io.read("XDATCAR", index=-1)`）。

### 9.5 用脚本批量生成输入

本仓库的 [`vasp_inputgen.py`](../scripts/vasp_inputgen.py) 可以按任务类型（relax/scf/band/...）生成带注释的 INCAR、KPOINTS，并从 `VASP_PP_PATH` 拼接 POTCAR；具体参数用 `python ../scripts/vasp_inputgen.py -h` 查看。高通量场景下也可以用 pymatgen 的 `MPRelaxSet` 或 atomate2，见 [Common 第6章](../../Common/tutorials/06-workflow-automation.md)。

---

## 参数选择与收敛建议

| 参数 | 常规优化 | 高精度（声子/弹性前） | 说明 |
|---|---|---|---|
| `IBRION` | 2（远离极小）→ 1（接近极小） | 1 或 2 | VTST 用户可用 `IBRION=3, POTIM=0, IOPT=1/2/3/7` |
| `ISIF` | 2（固定晶胞）/ 3（体相） | 3 | 二维用 3 + `LATTICE_CONSTRAINTS` |
| `ENCUT` | ISIF=2：收敛测试值；ISIF=3：≥ 1.3×ENMAX | ≥ 1.3×ENMAX | 变胞必须提高 |
| `PREC` | Accurate | Accurate | 减少应力噪声 |
| `EDIFF` | 1E-5 ~ 1E-6 | 1E-7 ~ 1E-8 | 影响受力精度 |
| `EDIFFG` | −0.01 ~ −0.02 | −1E-3 ~ −5E-3 | 负值按受力判断 |
| `NSW` | 50–200 | 100–200 | 不够就续算 |
| `POTIM` | 0.5（默认） | 0.1–0.3（IBRION=1 振荡时） | IBRION=2 会打印建议的 trialstep |
| `ISMEAR`/`SIGMA` | 半导体 0/0.05；金属 1/0.1–0.2 | 同左 | 金属优化不要用 −5（受力不准） |
| `NELMIN` | 4–8（IBRION=1 时） | 4–8 | 保证每步受力可靠 |
| `LREAL` | 大体系 Auto；小体系 .FALSE. | .FALSE. | Auto 有实空间投影误差 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `ZBRENT: fatal error in bracketing`（或 `ZBRENT: fatal error: bracketing interval incorrect`） | CG 线搜索失败：结构离极小很远、受力噪声大（EDIFF 太松）、能量面很平 | `cp CONTCAR POSCAR` 续算；`EDIFF` 收紧到 1E-6 或更小；减小 `POTIM`；改用 `IBRION = 1` 或 VTST FIRE（`IBRION=3, POTIM=0, IOPT=7`） |
| 跑满 NSW 但没有 "reached required accuracy" | 步数不够或优化在振荡 | 看 `grep "F=" OSZICAR` 是否持续下降；下降则续算，振荡则减小 POTIM / 换算法 |
| `WARNING: Sub-Space-Matrix is not hermitian in DAV` / `Error EDDDAV: Call to ZHEGV failed` | 原子靠得太近、步长过大使结构崩坏、或 ALGO 不合适 | 检查 CONTCAR 的原子间距；减小 POTIM；`ALGO = Normal` 改 `Fast` 或 `All`；删除旧 WAVECAR 从头开始 |
| `ISIF = 3` 后晶格明显偏小 | Pulay 应力（ENCUT 太低） | ENCUT ≥ 1.3×ENMAX；CONTCAR → POSCAR 重新优化直到体积不变；或改用 E–V 拟合 |
| 二维材料真空层被压缩 | 对 slab 直接用了 `ISIF = 3` | 用 `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`（≥ 6.4.3）；旧版本改用 E–a 扫描 |
| 写了 `LATTICE_CONSTRAINTS` 但 c 轴仍然在变 | VASP < 6.4.3 时 IBRION=1/2 下不读取该标签 | 升级到 6.5.1+；或用 §6.5 的替代方案 |
| `VERY BAD NEWS! internal error in subroutine IBZKPT` | 优化后晶格对称性改变，k 网格（如六角体系用了偶数 MP 网格）与晶格类型不匹配 | 改用 Γ 中心网格；重新生成 KPOINTS 后续算 |
| `internal error in subroutine PRICEL` / `Found some non-integer element in rotation matrix` | 结构有"近似对称"的数值噪声 | 用 spglib/vaspkit 对称化结构；调整 `SYMPREC`（如 1E-4 或 1E-6）；必要时 `ISYM = 0` |
| 分子用 `ISMEAR = -5` 报 `Tetrahedron method fails (number of k-points < 4)` | 四面体方法需要至少 4 个 k 点 | 孤立体系用 `ISMEAR = 0`、`SIGMA = 0.01~0.05` |
| 优化后对称性没有降低（明知应该畸变） | 对称性约束保持了初始高对称 | 给初始结构加小扰动（0.01–0.05 Å）；或降低 `SYMPREC`/关闭对称性做对照 |
| E–V 曲线锯齿状 | 基组/FFT 网格随体积不连续变化，k 点太稀 | 提高 ENCUT、`PREC = Accurate`、加密 k 点 |
| VASP 6.6.0 下自旋极化 r2SCAN 优化结果异常 | known issue 94：自旋极化 meta-GGA 受力错误 | 升级到 6.6.1，或给 6.6.0 打官方 `pot.F.patch` 后重编译 |

---

## 练习

1. 把 Si 示例中的 `ENCUT` 改为 250 eV（约等于 ENMAX），重复 `ISIF=3` 优化与重新优化，记录两轮的体积差，体会 Pulay 应力的影响。
2. 对 Si 分别用 `IBRION = 1、2、3`（IBRION=3 时设 `SMASS = 0.4`、`POTIM = 0.1`）优化同一个偏离的初始结构，比较离子步数。
3. 用 `run_eos.sh` + `fit_bm.py` 得到 Si 的 V0 和 B0，再与 §4 中 `ISIF=3` 的体积比较；把体积区间缩小到 ±2%，观察 B0' 的拟合稳定性。
4. 对 O2 分子分别用 `ISPIN = 1` 和 `ISPIN = 2, MAGMOM = 2*1` 优化，比较键长和能量，判断哪一个是基态。
5. 用 §6.5 的 E–a 扫描方法求单层 MoS2 晶格常数，并与 `LATTICE_CONSTRAINTS` 的结果对比。

---

## 参考资料

**官方文档**

- [VASP Wiki: IBRION](https://www.vasp.at/wiki/index.php/IBRION) — 优化、MD、有限差分等离子更新方式
- [VASP Wiki: ISIF](https://www.vasp.at/wiki/index.php/ISIF) — 应力计算与自由度表（ISIF=8 自 6.4.1 起）
- [VASP Wiki: Structure optimization](https://www.vasp.at/wiki/index.php/Structure_optimization) — RMM-DIIS / CG / 阻尼 MD 的原理与输出解读
- [VASP Wiki: Pulay stress](https://www.vasp.at/wiki/index.php/Pulay_stress) 与 [Volume relaxation](https://www.vasp.at/wiki/index.php/Volume_relaxation) — Pulay 应力的来源和修正、E–V 拟合流程
- [VASP Wiki: LATTICE_CONSTRAINTS](https://www.vasp.at/wiki/index.php/LATTICE_CONSTRAINTS) — 固定晶格分量及版本说明
- [VASP Wiki: EDIFFG](https://www.vasp.at/wiki/index.php/EDIFFG)、[POTIM](https://www.vasp.at/wiki/index.php/POTIM)、[PSTRESS](https://www.vasp.at/wiki/index.php/PSTRESS)、[NFREE](https://www.vasp.at/wiki/index.php/NFREE)
- [VASP Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — issue 34（LATTICE_CONSTRAINTS + ISIF=4/5）、issue 94（6.6.0 自旋极化 meta-GGA 受力）
- [VTST: Force based optimizers](https://theory.cm.utexas.edu/vtsttools/optimizers.html) — IOPT 各优化器参数
- [ASE: Equation of state](https://docs.ase-lib.org/ase/eos.html)、[pymatgen.analysis.eos](https://pymatgen.org/pymatgen.analysis.html)

**相关 GitHub 仓库（仅链接）**

- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — `MPRelaxSet` 输入集、`EOS` 拟合、Vasprun/Outcar 解析
- [henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode) — VTST 优化器（IOPT）源码补丁
- [Chengcheng-Xiao/VASP_OPT_AXIS](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS) — 旧版 VASP 固定晶格分量的补丁；6.4.3 起可用 LATTICE_CONSTRAINTS 代替
- [Ionizing/rsgrad](https://github.com/Ionizing/rsgrad) — 快速查看 OUTCAR 收敛过程、各步受力与磁矩
- [materialsproject/custodian](https://github.com/materialsproject/custodian) — 自动处理 ZBRENT、EDDDAV 等报错的规则库
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book) — VASP + ASE 经典教程，含 E–V 拟合与优化的例题
- [cndaqiang/DFT-EXERCISES](https://github.com/cndaqiang/DFT-EXERCISES) — Sholl《DFT: A Practical Introduction》习题的中文解答

**关键文献**

- F. Birch, *Phys. Rev.* 71, 809 (1947). DOI: [10.1103/PhysRev.71.809](https://doi.org/10.1103/PhysRev.71.809)
- F. D. Murnaghan, *PNAS* 30, 244 (1944). DOI: [10.1073/pnas.30.9.244](https://doi.org/10.1073/pnas.30.9.244)
- G. P. Francis, M. C. Payne, *J. Phys.: Condens. Matter* 2, 4395 (1990)（有限基组修正与 Pulay 应力）. DOI: [10.1088/0953-8984/2/19/007](https://doi.org/10.1088/0953-8984/2/19/007)
- E. Bitzek et al., *Phys. Rev. Lett.* 97, 170201 (2006)（FIRE 优化器）. DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)
- D. Sheppard, R. Terrell, G. Henkelman, *J. Chem. Phys.* 128, 134106 (2008)（NEB 优化器比较）. DOI: [10.1063/1.2841941](https://doi.org/10.1063/1.2841941)
- P. Haas, F. Tran, P. Blaha, *Phys. Rev. B* 79, 085104 (2009)（各泛函晶格常数与体模量基准）. DOI: [10.1103/PhysRevB.79.085104](https://doi.org/10.1103/PhysRevB.79.085104)

---

[← 上一章](03-convergence.md) | [返回目录](../README.md) | [下一章 →](05-band-dos.md)
