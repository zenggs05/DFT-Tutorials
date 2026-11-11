# 第13章 弹性常数与力学性质（Elastic Constants & Mechanical Properties with VASP）

本章讲解如何用 VASP 计算晶体的弹性常数（elastic constants / elastic tensor C_ij）：内置的应力-应变法（stress-strain，`IBRION = 6` + `ISIF = 3`）与能量-应变法（energy-strain，VASPKIT 201），并从 C_ij 出发判断 Born 力学稳定性、计算 Voigt-Reuss-Hill 平均的体模量、剪切模量、杨氏模量、泊松比、Pugh 比和德拜温度。最后给出二维材料面内弹性常数（单位 N/m）的完整流程，以及用 ElATools 等工具做弹性各向异性可视化的方法。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 `IBRION = 6`、`ISIF = 3` 的弹性常数 INCAR，知道 `NFREE`、`POTIM`、`ENCUT`、`EDIFF` 各自为什么这样取；
- 从 OUTCAR 中读出 `TOTAL ELASTIC MODULI (kBar)`，**除以 10 换算为 GPa**，并把 VASP 的分量顺序 XX YY ZZ **XY YZ ZX** 正确重排成 Voigt 顺序 xx yy zz **yz zx xy**；
- 区分钳位离子（clamped-ion / rigid-ion）与弛豫离子（relaxed-ion）弹性常数；
- 用 Born 判据判断力学稳定性，计算 B、G、E、ν、B/G、各向异性指数、声速和德拜温度；
- 用能量-应变法（含 VASPKIT 201）交叉验证，并计算二维材料的面内弹性常数；
- 把 C_ij 交给 ElATools / ELATE 等工具做方向依赖的杨氏模量、泊松比可视化。

| 方法 | VASP 设置 | 优点 | 局限 | 本章 |
|---|---|---|---|---|
| 应力-应变（VASP 内置有限差分） | `IBRION = 6`，`ISIF = 3` | 一次运行得到完整 6×6 张量，自动包含离子弛豫贡献，利用对称性 | 应力对平面波基组敏感，ENCUT 要显著提高并系统收敛 | §3 |
| 能量-应变 | 一组应变结构，`ISIF = 2` 弛豫原子 | 只用总能，对 Pulay 应力不敏感；VASPKIT 201 自动化 | 计算次数随独立常数个数增加；需要拟合 | §6 |
| 应力-应变（外部工作流） | 一组应变结构，`ISIF = 2`，读应力 | 应变幅度可控；pymatgen / atomate2 全自动 | 需要自己组织拟合 | §6.4 |
| 二维面内弹性常数 | 面内应变，固定 z 方向 | 单位 N/m，与实验可比 | 3D 张量含真空，不能直接用 | §7 |
| DFPT（`IBRION = 7/8`） | — | — | 官方说明：DFPT 不计算钳位离子弹性张量（未实现对应变的微扰） | — |

> [!NOTE]
> 本章默认你已经完成 [第3章](03-convergence.md) 的 ENCUT / k 点收敛测试和 [第4章](04-relax-eos.md) 的结构优化。`IBRION = 6` 的有限位移部分与 [第12章 声子](12-phonon.md) 共用同一套机制，建议先读第12章中关于有限差分的部分。

---

## 2. 原理速览

### 2.1 胡克定律与 Voigt 记号

在线弹性范围内，应力张量与应变张量通过四阶弹性张量相联系：

$$
\sigma_{ij} = \sum_{kl} C_{ijkl}\,\varepsilon_{kl}
$$

由于应力、应变都是对称张量，且存在应变能函数，C_ijkl 最多只有 21 个独立分量，通常用 6×6 的 Voigt 矩阵 C_IJ 表示：

| Voigt 下标 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| 张量分量 | xx | yy | zz | **yz** | **zx (xz)** | **xy** |

Voigt 记号中剪切应变使用**工程剪切应变** γ：ε_4 = 2ε_yz，ε_5 = 2ε_zx，ε_6 = 2ε_xy。于是

$$
\sigma_I = \sum_J C_{IJ}\,\varepsilon_J ,\qquad
\frac{\Delta E}{V_0} = \frac{1}{2}\sum_{IJ} C_{IJ}\,\varepsilon_I\,\varepsilon_J + O(\varepsilon^3)
$$

第二个式子是能量-应变法的出发点（§6）；第一个式子是应力-应变法的出发点（§3）。柔度矩阵 S = C⁻¹。

### 2.2 独立弹性常数的个数

| 晶系 | 独立 C_ij 个数 | 独立分量 |
|---|---|---|
| 立方 | 3 | C11, C12, C44 |
| 六方 | 5 | C11, C12, C13, C33, C44（C66 = (C11−C12)/2） |
| 三方（I/II） | 6 / 7 | C11, C12, C13, C14, C33, C44（II 类再加 C15） |
| 四方（I/II） | 6 / 7 | C11, C12, C13, C33, C44, C66（II 类再加 C16） |
| 正交 | 9 | C11, C22, C33, C12, C13, C23, C44, C55, C66 |
| 单斜 | 13 | — |
| 三斜 | 21 | 全部 |

> [!IMPORTANT]
> 这些"标准形式"都假设晶体按 IEEE / 标准取向放置（例如六方 c 轴沿 z，立方晶轴沿 x、y、z）。如果你的 POSCAR 把晶格旋转了，VASP 给出的张量也是旋转后笛卡尔坐标系中的张量，看起来会有很多"不该有"的非零分量。做弹性计算前，先用 `vaspkit` 603/604 或 pymatgen 把结构标准化。

### 2.3 钳位离子与弛豫离子

对原胞中原子不全处在对称性固定位置上的晶体（如金刚石结构的 Si、纤锌矿 ZnO），宏观应变会引起原子之间的内部位移（internal strain）。因此有两种弹性常数：

- **钳位离子（clamped-ion，刚性离子）**：只按应变仿射地移动原子，不再弛豫；
- **弛豫离子（relaxed-ion）**：在每个应变下原子重新弛豫到受力为零，这是与实验（静态/低频）可比的量。

两者的差别来自离子 Hessian（力常数矩阵）Φ 与内应变张量 Λ（应变引起的原子受力）：

$$
C^{\mathrm{relax}}_{IJ} = C^{\mathrm{clamp}}_{IJ} - \frac{1}{\Omega}\sum_{mn}\Lambda_{mI}\,\big(\Phi^{-1}\big)_{mn}\,\Lambda_{nJ}
$$

（Ω 为晶胞体积，求逆时去掉三个声学平移模；此式只给出形式，具体符号约定见 Wu, Vanderbilt & Hamann 2005。）离子弛豫总是使弹性常数**减小**（或不变）。fcc Al 原胞只有一个原子，这一修正为零；对 Si 的 C44，这一修正就很显著。

### 2.4 VASP 的应力符号与单位

VASP 把应力定义为能量对应变的**负**导数（见官方 `ISIF` 页面），OUTCAR 中的应力以 kBar（kB）为单位：

- 对角元为**正**：体系处于压缩状态，"想要膨胀"；为负：处于拉伸状态。
- 1 kBar = 0.1 GPa，所以 **GPa = kBar / 10**。

OUTCAR 中打印的弹性模量同样以 kBar 为单位，已经按通常的物理符号给出（稳定晶体的 C11 为正值）。

### 2.5 `IBRION = 6` + `ISIF ≥ 3` 做了什么

`IBRION = 6` 是"利用对称性的有限差分"（`IBRION = 5` 不用对称性，官方已不推荐）。当 `ISIF ≥ 3` 时，除了逐个移动原子（得到力常数和内应变张量）之外，VASP 还对晶格施加六种有限畸变，由应力-应变关系得到弹性张量（Le Page & Saxe 2002 的方法）。OUTCAR 中依次给出：

| OUTCAR 标题 | 含义 |
|---|---|
| `SYMMETRIZED ELASTIC MODULI (kBar)` | 钳位离子（刚性离子）弹性张量，已对称化 |
| `ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)` | 离子弛豫贡献（由 Hessian 求逆乘以内应变张量得到） |
| `TOTAL ELASTIC MODULI (kBar)` | 两者之和 = 弛豫离子弹性张量，**通常引用的就是这一组** |

官方特别强调：应力张量需要足够大的 ENCUT 才能收敛，默认截断能通常要提高约 30%，并建议以约 15% 的步长继续增大 ENCUT，直到弹性常数完全收敛。

---

## 3. 实战一：fcc Al 的应力-应变法（IBRION = 6）

完整输入在 [../examples/13-elastic-al-ibrion6/README.md](../examples/13-elastic-al-ibrion6/README.md)。

### 3.1 计算前的收敛测试

弹性常数是能量对应变的二阶导数，对数值精度的要求明显高于总能和结构优化：

| 参数 | 建议 | 原因 |
|---|---|---|
| `ENCUT` | ≥ 1.3 × max(ENMAX)，再以约 15% 步长测试 | 应力受 Pulay 误差影响大 |
| k 点 | 比总能收敛再加密；金属更密（Al 原胞可用 21×21×21 或更高） | 应变降低对称性，C44 对 k 点最敏感 |
| `ISMEAR`/`SIGMA` | 半导体 `ISMEAR = 0`、`SIGMA = 0.01–0.05`；金属 `ISMEAR = 1`、`SIGMA = 0.05–0.2` 并检查 entropy T*S | 展宽会改变费米面附近的能量 |
| `EDIFF` | 1E-7 ~ 1E-8 | 有限差分需要准确的应力和力 |
| `PREC` | `Accurate` | 力与应力的 FFT 网格更精确 |
| `LREAL` | `.FALSE.` | 实空间投影引入的噪声会进入二阶导数 |

收敛判据建议：ENCUT 或 k 点再提高一档，各 C_ij 变化小于 1–2 GPa（或 1–2%）。可以用 [../scripts/conv_test.sh](../scripts/conv_test.sh) 批量生成测试目录（用 `bash ../scripts/conv_test.sh -h` 查看用法），只需把测试的输出量从总能换成 C_ij。

### 3.2 第 1 步：用同一套参数充分弛豫

弹性张量是在**零应力的平衡结构**上定义的。若参考结构残余应力较大，有限差分得到的 C_ij 会混入应力项，立方晶体的 C12 与 C44 甚至会出现明显偏差。先做 `ISIF = 3` 变胞优化：

```bash
# INCAR.relax —— 与第 2 步保持相同的 ENCUT / k 点 / SIGMA
SYSTEM  = Al fcc ISIF=3 relaxation before elastic constants

# electronic
PREC    = Accurate
ENCUT   = 400          # ~1.66 x ENMAX(Al)=240 eV; stresses need a large basis
EDIFF   = 1E-8
ISMEAR  = 1            # Methfessel-Paxton for a metal
SIGMA   = 0.10         # check entropy T*S < 1 meV/atom in OUTCAR
LREAL   = .FALSE.

# ionic
IBRION  = 2
ISIF    = 3
NSW     = 60
EDIFFG  = -1E-3        # eV/A; for one atom only the stress matters

# output
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

```bash
cat $VASP_PP_PATH/potpaw_PBE/Al/POTCAR > POTCAR        # 仓库不提供 POTCAR
cp INCAR.relax INCAR
mpirun -np 16 vasp_std > relax1.out
cp CONTCAR POSCAR && mpirun -np 16 vasp_std > relax2.out   # 再弛豫一轮消除 Pulay 应力
grep "external pressure" OUTCAR | tail -1                   # 应接近 0（|P| < 1 kB）
grep "volume of cell" OUTCAR | tail -2                      # 两轮体积几乎不变
```

> [!TIP]
> 变胞优化的 PAW 基组在优化过程中不更新（官方 `ISIF` 页面的警告），所以至少再用 CONTCAR 重启一次，直到体积不变。这一点在 [第4章](04-relax-eos.md) 中有详细说明。

### 3.3 第 2 步：弹性常数计算的 INCAR

```bash
# INCAR.elastic —— 在弛豫后的结构上运行：cp CONTCAR POSCAR
SYSTEM  = Al fcc elastic constants IBRION=6 ISIF=3

# electronic: identical to INCAR.relax
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-8
ISMEAR  = 1
SIGMA   = 0.10
LREAL   = .FALSE.

# finite differences
IBRION  = 6            # finite differences WITH symmetry
ISIF    = 3            # >= 3: also strain the lattice -> elastic tensor
NFREE   = 2            # central differences (+/-); 4 is more accurate and 2x cost
POTIM   = 0.015        # ionic displacement in A (internal-strain part)
NSW     = 1            # only needed by very old versions; harmless

# output
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `IBRION` | 6 | 有限差分 + 对称性；5 为无对称性版本（已弃用，但它是唯一支持 Selective dynamics 的选项） |
| `ISIF` | 3 | `ISIF ≥ 3` 才计算弹性张量与内应变张量；`ISIF = 2` 只算 Γ 点声子 |
| `NFREE` | 2 或 4 | 每个方向的位移数：2 为中心差分 ±POTIM；4 再加 ±2×POTIM，误差更小、代价加倍；1 强烈不推荐，3、5 等无效值会被重置为 2 |
| `POTIM` | 0.015 | 离子位移幅度（Å），对 `IBRION = 5/6` 若给得过大，VASP 5.1 起会自动重置为 0.015 Å |
| `ENCUT` | ≥ 1.3×ENMAX，系统收敛 | 本例 400 eV；正式研究应测试 400/460/530 eV |
| `EDIFF` | 1E-8 | 官方提示有限差分的默认 EDIFF（1E-4）误差不可接受 |
| `ISYM` | 默认 | `IBRION = 6` 依赖对称性减少位移数；不要关掉 |
| `ADDGRID` | 不设 | 官方建议在有限差分计算中不要未经测试就打开 |
| `NCORE`/`KPAR` | 按常规设置 | 与 `LEPSILON` 不同，`IBRION = 6` 没有 NCORE 限制；小体系 k 点多，`KPAR` 效果好 |

> [!NOTE]
> **混合泛函**：`LHFCALC = .TRUE.` 时 `ISIF` 的默认值是 0（不计算应力），必须显式写 `ISIF = 3`。杂化泛函的应力计算很贵，通常先用 PBE / PBEsol / r2SCAN 计算弹性常数。
>
> **VASP 6.6 断点续算**：`IBRION = 6` 在 VASP ≥ 6.6.0（需 HDF5）会默认把位移数据写入 `vaspcheckfd.h5`，作业中断后可设 `CHECKPOINT_FD = CONTINUE` 续算，也可以 `PREPARE` + `SINGLE` 拆分成多个独立作业。官方文档以声子为例；用于 `ISIF = 3` 的应变部分前，建议先用小体系验证。

### 3.4 运行与监控

```bash
cp INCAR.elastic INCAR
mpirun -np 16 vasp_std > elastic.out
grep -c "free  energy   TOTEN" OUTCAR      # 已完成的构型数（含平衡结构）
grep "Iteration" OUTCAR | tail -1          # 当前所在的离子步/电子步
```

对 1 个原子的原胞，位移部分很少，大部分时间花在六种晶格畸变（以及中心差分的正负两侧）上。超胞越大、对称性越低，位移数越多，计算量近似正比于"位移数 × 单点计算时间"。

### 3.5 读取 OUTCAR 中的弹性张量

```bash
grep -A 10 "SYMMETRIZED ELASTIC MODULI" OUTCAR   # 钳位离子
grep -A 10 "ELASTIC MODULI CONTR FROM IONIC" OUTCAR
grep -A 10 "TOTAL ELASTIC MODULI" OUTCAR         # 总弹性张量（最常用）
```

输出块的结构如下（示意输出，数值仅供参考）：

```text
 TOTAL ELASTIC MODULI (kBar)
 Direction    XX          YY          ZZ          XY          YZ          ZX
 --------------------------------------------------------------------------------
 XX        1050.0000    590.0000    590.0000      0.0000      0.0000      0.0000
 YY         590.0000   1050.0000    590.0000      0.0000      0.0000      0.0000
 ZZ         590.0000    590.0000   1050.0000      0.0000      0.0000      0.0000
 XY           0.0000      0.0000      0.0000    300.0000      0.0000      0.0000
 YZ           0.0000      0.0000      0.0000      0.0000    300.0000      0.0000
 ZX           0.0000      0.0000      0.0000      0.0000      0.0000    300.0000
 --------------------------------------------------------------------------------
```

判断是否成功：

- 计算正常结束（OUTCAR 末尾有 `General timing and accounting informations`）；
- 张量对称（C_IJ ≈ C_JI），立方晶体满足 C11=C22=C33、C12=C13=C23、C44=C55=C66，其余元素只有小于 0.5 GPa 的数值噪声；
- 三个块中 `CONTR FROM IONIC RELAXATION` 对单原子原胞应为 0。

### 3.6 单位换算与 Voigt 重排（最常见的错误）

VASP 的行列顺序是 **XX YY ZZ XY YZ ZX**，而 Voigt 顺序是 **xx yy zz yz zx xy**。对应关系：

| VASP 第 n 行/列 | 1 (XX) | 2 (YY) | 3 (ZZ) | 4 (XY) | 5 (YZ) | 6 (ZX) |
|---|---|---|---|---|---|---|
| Voigt 下标 | 1 | 2 | 3 | **6** | **4** | **5** |

立方晶体三个剪切分量相等，重排与否看不出区别；但对六方晶体，VASP 第 4 行（XY）给出的是 **C66 = (C11−C12)/2**，第 5、6 行（YZ、ZX）才是 C44 = C55。直接把 VASP 第 4 行当作 C44，就会把 C66 错当成 C44，后续 Born 判据和 VRH 平均全部出错。

重排的核心只有几行 numpy（完整脚本见算例目录中的 `voigt_vrh.py`）：

```python
import numpy as np

VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]   # Voigt index i takes VASP row VASP_TO_VOIGT[i]

def read_total_moduli(outcar="OUTCAR", header="TOTAL ELASTIC MODULI (kBar)"):
    lines = open(outcar, errors="replace").readlines()
    start = max(i for i, l in enumerate(lines) if header in l)
    rows = [l.split()[1:7] for l in lines[start + 1:start + 12]
            if l.split() and l.split()[0] in ("XX", "YY", "ZZ", "XY", "YZ", "ZX")]
    return np.array(rows[:6], dtype=float)          # kBar, VASP order

c_vasp = read_total_moduli()
c = c_vasp[np.ix_(VASP_TO_VOIGT, VASP_TO_VOIGT)] / 10.0   # -> GPa, Voigt order
c = 0.5 * (c + c.T)                                        # remove numerical asymmetry
print(np.round(c, 1))
```

在算例目录中直接运行：

```bash
python3 voigt_vrh.py -f OUTCAR --poscar CONTCAR    # 总张量 + VRH + 德拜温度
python3 voigt_vrh.py --clamped                     # 钳位离子张量
```

### 3.7 其他读取方式

- **py4vasp**（需要编译时开启 HDF5，生成 `vaspout.h5`）：

  ```python
  from py4vasp import Calculation
  calc = Calculation.from_path(".")
  print(calc.elastic_modulus)              # 钳位离子与弛豫离子两组张量（kBar）
  data = calc.elastic_modulus.read()       # 字典；不同 py4vasp 版本的返回格式略有差别，用 help() 确认
  ```

- **pymatgen**：`Outcar("OUTCAR").read_elastic_tensor()` 把 `TOTAL ELASTIC MODULI` 读入 `outcar.data["elastic_tensor"]`（仍是 kBar、VASP 顺序），重排并除以 10 后交给 `ElasticTensor.from_voigt(c)`，即可使用 `k_vrh`、`g_vrh`、`homogeneous_poisson`、`universal_anisotropy`、`debye_temperature(structure)` 等属性/方法（注意 `y_mod` 返回的是 SI 单位 Pa）。
- **VASPKIT 203**：直接读取 OUTCAR 中的弹性常数并输出力学性质与稳定性判据；202 读取用户整理好的 `ELASTIC_TENSOR.in`。
- **仓库脚本**：[../scripts/elastic_post.py](../scripts/elastic_post.py) 完成单位换算、Voigt 重排、各晶系 Born 判据、VRH、德拜温度和声速，用 `python ../scripts/elastic_post.py -h` 查看用法。

---

## 4. 力学稳定性判据（Born 判据）

无应力晶体力学稳定的充要条件是：弹性矩阵 C 正定，即**六个本征值全部为正**（Mouhat & Coudert 2014）。这一通用判据对任何晶系都适用，`voigt_vrh.py` 就是这样判断的。各晶系也有等价的显式形式（下标为 Voigt 记号）：

| 晶系 | 充要条件 |
|---|---|
| 立方 | C11 − C12 > 0；C11 + 2C12 > 0；C44 > 0 |
| 六方 | C11 > \|C12\|；2C13² < C33(C11 + C12)；C44 > 0；C66 = (C11 − C12)/2 > 0 |
| 四方 I（4mm, −42m, 422, 4/mmm） | C11 > \|C12\|；2C13² < C33(C11 + C12)；C44 > 0；C66 > 0 |
| 三方 I（3m, −3m, 32） | C11 > \|C12\|；C44 > 0；C13² < ½C33(C11 + C12)；C14² < ½C44(C11 − C12) |
| 正交 | C11 > 0；C11C22 > C12²；C11C22C33 + 2C12C13C23 − C11C23² − C22C13² − C33C12² > 0；C44 > 0；C55 > 0；C66 > 0 |
| 单斜、三斜及其他 | 用本征值判据 |

> [!WARNING]
> 这些判据针对**零应力**状态。在有外压的结构上（例如 `PSTRESS ≠ 0` 优化得到的高压相），需要使用包含应力项的修正判据（Born 判据在应力下的推广），不能直接套用上表。

力学不稳定（某个本征值为负）可能有两种原因：真实的力学不稳定性（结构会自发变形），或数值问题（参考结构没弛豫好、ENCUT/k 点不够、Voigt 顺序弄错）。先排除后者。

---

## 5. 多晶平均与派生力学量

### 5.1 Voigt–Reuss–Hill 平均

实验常测量多晶样品的体模量 B 和剪切模量 G。Voigt 平均假设均匀应变（上限），Reuss 平均假设均匀应力（下限），Hill 取两者的算术平均：

$$
9B_V = (C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{23}+C_{31})
$$

$$
15G_V = (C_{11}+C_{22}+C_{33}) - (C_{12}+C_{23}+C_{31}) + 3(C_{44}+C_{55}+C_{66})
$$

$$
\frac{1}{B_R} = (S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{23}+S_{31}),\qquad
\frac{15}{G_R} = 4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{23}+S_{31}) + 3(S_{44}+S_{55}+S_{66})
$$

$$
B_H = \frac{B_V+B_R}{2},\qquad G_H = \frac{G_V+G_R}{2}
$$

### 5.2 杨氏模量、泊松比与经验指标

$$
E = \frac{9BG}{3B+G},\qquad \nu = \frac{3B-2G}{2(3B+G)}
$$

| 指标 | 公式 | 常见解读（经验性，需谨慎） |
|---|---|---|
| Pugh 比 | B/G | > 1.75 倾向延性，< 1.75 倾向脆性（Pugh 1954） |
| Cauchy 压力（立方） | C12 − C44 | 正值常与金属键/延性相关，负值与方向性共价键相关（Pettifor 1992） |
| Zener 各向异性（立方） | A = 2C44/(C11 − C12) | A = 1 为各向同性 |
| 普适各向异性指数 | A^U = 5G_V/G_R + B_V/B_R − 6 | A^U = 0 为各向同性（Ranganathan & Ostoja-Starzewski 2008） |

### 5.3 声速与德拜温度

由 Hill 平均的 B、G 和密度 ρ 可估算各向同性声速与德拜温度（Anderson 1963）：

$$
v_l = \sqrt{\frac{B+\frac{4}{3}G}{\rho}},\qquad v_t = \sqrt{\frac{G}{\rho}},\qquad
v_m = \left[\frac{1}{3}\left(\frac{2}{v_t^3}+\frac{1}{v_l^3}\right)\right]^{-1/3}
$$

$$
\Theta_D = \frac{h}{k_B}\left(\frac{3n}{4\pi}\right)^{1/3} v_m
$$

其中 n 为单位体积内的原子数（原子数 / 晶胞体积）。

### 5.4 算例：用 Al 的实验弹性常数检验公式

把 Al 的室温实验值 C11 = 106.8、C12 = 60.7、C44 = 28.2 GPa（Kittel《固体物理导论》表格数据，近似）代入上面的公式（即 `voigt_vrh.py` 的计算），得到：

| 量 | 结果 | 对照 |
|---|---|---|
| B (V = R = H) | 76.1 GPa | 立方晶体 B = (C11 + 2C12)/3 |
| G_V / G_R / G_H | 26.1 / 25.9 / 26.0 GPa | Al 接近各向同性，三者很接近 |
| E | 70.1 GPa | 多晶 Al 的常用工程值约 70 GPa |
| ν | 0.347 | 工程值约 0.33–0.35 |
| B/G | 2.92 | 延性金属 |
| A（Zener） | 1.22 | — |
| v_l / v_t | 6406 / 3105 m/s | — |
| Θ_D | 约 407 K | 量热法测得约 428 K（Kittel 表格，近似） |

这组数说明脚本与公式是自洽的；你的 PBE 结果与实验通常相差 5–15%（PBE 一般略微高估晶格常数、低估弹性常数）。

---

## 6. 实战二：能量-应变法

### 6.1 原理

对平衡结构施加一组小应变 δ（典型 ±0.5%、±1%、±1.5%），每个应变下固定晶胞、弛豫原子（`ISIF = 2`，得到弛豫离子常数），把 ΔE(δ) 拟合成二次多项式，二次项系数就是若干 C_ij 的线性组合。新晶格矢量（行矢量）为

$$
\mathbf{A}' = \mathbf{A}\,(\mathbf{I} + \boldsymbol{\varepsilon}),\qquad
\boldsymbol{\varepsilon} =
\begin{pmatrix}
\varepsilon_1 & \varepsilon_6/2 & \varepsilon_5/2\\
\varepsilon_6/2 & \varepsilon_2 & \varepsilon_4/2\\
\varepsilon_5/2 & \varepsilon_4/2 & \varepsilon_3
\end{pmatrix}
$$

立方晶体只需三组应变：

| 应变 (ε1…ε6) | ΔE / V0 | 得到 |
|---|---|---|
| (δ, δ, δ, 0, 0, 0) | (3/2)(C11 + 2C12) δ² | 体模量组合 C11 + 2C12 |
| (δ, −δ, 0, 0, 0, 0) | (C11 − C12) δ² | 四方剪切 C11 − C12 |
| (0, 0, 0, δ, δ, δ) | (3/2) C44 δ² | C44 |

低对称晶系需要更多组应变（VASPKIT 会根据空间群自动选择）。

> [!IMPORTANT]
> 若参考结构存在残余应力 σ⁰，ΔE 中会出现线性项 V0 Σ σ⁰_I ε_I。拟合时保留一次项可以吸收它，但更好的做法是先把结构弛豫到零应力；一次项系数可以作为"参考结构是否弛豫好"的检查。

### 6.2 VASPKIT 201 自动化流程

VASPKIT（1.x 版本）用任务 201 实现能量-应变法，流程分三步：

1. 准备**完全弛豫且标准化的惯用胞** POSCAR（可用 vaspkit 603/604 生成标准惯用胞）、较密的 KPOINTS（vaspkit 102）以及 INCAR。INCAR 中 `ISIF = 2`、`IBRION = 2`、`NSW` 适当、`EDIFF = 1E-7~1E-8`、`PREC = Accurate`、`LREAL = .FALSE.`。
2. 写 `VPKIT.in`（第 1 行 1 = 前处理 / 2 = 后处理；第 2 行 `3D` 或 `2D`；第 3 行应变点数；第 4 行应变值），例如：

   ```text
   1                                                 ! 1: pre-processing, 2: post-processing
   3D                                                ! 3D bulk or 2D sheet
   7                                                 ! number of strain values
   -0.015 -0.010 -0.005 0.000 0.005 0.010 0.015      ! strain values
   ```

   运行 `vaspkit -task 201`（或交互菜单输入 201），VASPKIT 按晶系自动建立 `C11_C12_I/strain_-0.015/` 之类的目录。
3. 批量提交所有目录中的 VASP 计算；全部完成后把 `VPKIT.in` 第 1 行改为 2，再运行一次 201，得到刚度矩阵、柔度矩阵、VRH 平均、Pugh 比、Cauchy 压力、声速、德拜温度、本征值稳定性判据等。

> [!TIP]
> VASPKIT 的输出很全面，但务必自己检查拟合精度（它会打印每组拟合的残差）和单个目录是否都已收敛（`grep -L "reached required accuracy" */strain_*/OUTCAR`）。

### 6.3 应力-应变法 vs 能量-应变法

| 对比项 | 应力-应变（`IBRION = 6`） | 能量-应变（VASPKIT 201） |
|---|---|---|
| 需要的量 | 应力（一阶导数） | 总能 |
| 计算次数 | 少（一次运行内部完成） | 多（每组应变 5–9 个点 × 组数） |
| 对 ENCUT 的敏感度 | 高（Pulay 应力） | 较低（固定晶胞比较能量） |
| 离子弛豫 | 由内应变张量解析处理 | 每个应变点显式弛豫（`ISIF = 2`） |
| 适合 | 高对称晶体、快速得到完整张量 | 交叉验证、低对称或二维材料、对应力精度没把握时 |

两种方法得到的 C_ij 相差 < 3–5% 时，可以认为结果可靠。

### 6.4 用 pymatgen / atomate2 自建工作流（简介）

pymatgen 的 `DeformedStructureSet` 可以按指定应变生成变形结构，`ElasticTensor.from_independent_strains(strains, stresses, eq_stress=...)` 由多组（应变，应力）拟合弹性张量；atomate2 的 `ElasticMaker` 把"弛豫 → 生成变形 → VASP 计算 → 拟合"做成了完整工作流，适合高通量。Materials Project 的弹性数据库（de Jong et al. 2015）就是用这类应力-应变工作流计算的。

> [!NOTE]
> 从 OUTCAR 读应力时同样要注意符号与单位：pymatgen 的 `Vasprun` / `Outcar` 读出的应力已经按其约定处理过，混用自己写的解析代码时要统一约定（VASP 中压缩为正）。

---

## 7. 实战三：二维材料的面内弹性常数（单层 MoS2）

完整输入在 [../examples/13-elastic-mos2-2d/README.md](../examples/13-elastic-mos2-2d/README.md)。

### 7.1 为什么二维材料要用 N/m

二维材料的晶胞包含人为设定的真空层，按"晶胞体积"归一化的 C_ij（GPa）会随真空厚度变化，没有物理意义。应改用"单位面积"的二维弹性常数：

$$
\frac{\Delta E}{A_0} = \frac{1}{2}C_{11}\varepsilon_1^2 + \frac{1}{2}C_{22}\varepsilon_2^2 + C_{12}\,\varepsilon_1\varepsilon_2 + \frac{1}{2}C_{66}\,\varepsilon_6^2
$$

单位为 N/m（1 eV/Å² = 16.0218 N/m）。若已用 3D 方法得到 C_ij（GPa），可换算为 C^2D(N/m) = C^3D(GPa) × L_z(Å) / 10，其中 L_z 是晶胞 z 方向长度；只有面内分量（11、22、12、66）有意义。

### 7.2 应变模式

| 模式 | 面内应变 (ε1, ε2, ε6) | ΔE / A0 | 得到 |
|---|---|---|---|
| xx 单轴 | (δ, 0, 0) | ½ C11 δ² | C11 |
| yy 单轴 | (0, δ, 0) | ½ C22 δ² | C22 |
| 双轴 | (δ, δ, 0) | ½ (C11 + C22 + 2C12) δ² | C12 |
| 剪切 | (0, 0, δ) | ½ C66 δ² | C66 |

对六方材料（MoS2、石墨烯、h-BN），C11 = C22，且 C66 = (C11 − C12)/2，可以作为自洽检查。

### 7.3 操作步骤

```bash
cd VASP/examples/13-elastic-mos2-2d
# 0) 先求面内平衡晶格常数（第4章：ISIF=3 + LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.，VASP >= 6.4.3）
# 1) 生成 4 种模式 x 7 个应变 = 28 个目录，并提交
bash run_strain.sh            # 只建目录
bash run_strain.sh submit     # 建目录并 sbatch
# 2) 收集能量并拟合
grep -L "reached required accuracy" strain_*/*/OUTCAR    # 应无输出
python3 fit_2d_elastic.py --plot
```

每个应变点的 INCAR 关键设置：`ISIF = 2`（固定应变后的晶胞，只弛豫原子）、`EDIFF = 1E-7`、`EDIFFG = -1E-3`、`PREC = Accurate`、`LREAL = .FALSE.`；KPOINTS 用 Γ 中心 18×18×1。能量差只有 meV 量级，所以电子收敛必须很严格。

`fit_2d_elastic.py` 输出的格式（示意输出，数值仅供参考）：

```text
    mode  npts    c1 (eV)   2c2/A0 (N/m)  max|res| (meV)
      xx     7     0.0002         13x.xx           0.004
      yy     7    -0.0001         13x.xx           0.004
 biaxial     7     0.0001         3xx.xx           0.003
   shear     7     0.0000          4x.xx           0.005
C11 = ... N/m   C22 = ... N/m   C12 = ... N/m   C66 = ... N/m
```

### 7.4 二维力学量与稳定性

$$
Y_x = \frac{C_{11}C_{22}-C_{12}^2}{C_{22}},\qquad
Y_y = \frac{C_{11}C_{22}-C_{12}^2}{C_{11}},\qquad
\nu_{xy} = \frac{C_{12}}{C_{22}}
$$

六方材料中 Y = (C11² − C12²)/C11，ν = C12/C11。二维 Born 判据：C11 > 0、C22 > 0、C11C22 > C12²、C66 > 0（六方：C11 > |C12|、C66 > 0）。

量级参考：文献中单层 MoS2 的 DFT 结果大约为 C11 ≈ 130–150 N/m、C12 ≈ 30–50 N/m，Y ≈ 120–140 N/m，ν ≈ 0.25（不同泛函与参数差异明显，仅供参考）；AFM 纳米压痕实验给出 180 ± 60 N/m（Bertolazzi et al. 2011）。

### 7.5 二维材料能否直接用 IBRION = 6？

可以，但要注意：

- 只取面内分量（VASP 顺序中的 XX、YY、XY 行列）乘以 L_z/10 换成 N/m；
- 真空方向的分量（ZZ、YZ、ZX）接近 0，整个 6×6 矩阵近乎奇异，**不要**对它做 VRH 平均或求逆，只对 3×3 面内块（11、22、66）做分析；
- 参考结构需先用 `LATTICE_CONSTRAINTS` 只弛豫面内晶格；
- VASPKIT 201 的 `VPKIT.in` 第 2 行写 `2D` 时会自动按二维处理。

---

## 8. 弹性各向异性与可视化

单晶的杨氏模量、剪切模量、泊松比、线压缩率都随方向变化。任意单位方向 **n** 上的杨氏模量为

$$
\frac{1}{E(\mathbf{n})} = \sum_{ijkl} n_i\,n_j\,n_k\,n_l\,S_{ijkl}
$$

其中 S_ijkl 是四阶柔度张量（由 Voigt 柔度矩阵转换时，剪切分量要乘以 1/2 或 1/4 的因子）。立方晶体可化简为 1/E(**n**) = S11 − 2(S11 − S12 − ½S44)(n1²n2² + n2²n3² + n3²n1²)。用 Al 的实验常数计算得到 E[100] ≈ 62.8 GPa、E[110] ≈ 71.7 GPa、E[111] ≈ 75.3 GPa，可见即使"近似各向同性"的 Al 也有约 20% 的方向差异。

常用工具：

| 工具 | 输入 | 能做什么 |
|---|---|---|
| [ElATools](https://github.com/shahramyalameha/ElATools) | 6×6 刚度矩阵（Voigt，GPa），支持 2D 与 3D | 本征值稳定性检查、E/G/ν/线压缩率的 3D 曲面与截面、求解 Christoffel 方程得到声速、硬度等，配有终端交互界面 |
| [ELATE](https://progs.coudert.name/elate) | 6×6 刚度矩阵 | 在线绘制方向依赖的力学性质（Gaillac, Pullumbi & Coudert 2016） |
| VASPKIT 201/203 | OUTCAR 或 VPKIT 流程 | 各向异性最小/最大值、VRH、稳定性判据 |
| pymatgen `ElasticTensor` | 6×6 矩阵 | `directional_elastic_mod`、`universal_anisotropy` 等程序化接口 |

> [!TIP]
> 把 C_ij 交给这些工具之前，确认已经**换算成 GPa 并重排为 Voigt 顺序**。很多"奇怪的各向异性图"都源于把 VASP 的 XY 行当成了第 4 行。

---

## 参数选择与收敛建议

| 参数 | 推荐取值 | 说明 |
|---|---|---|
| `ENCUT` | ≥ 1.3×ENMAX，并以约 15% 步长测试 | 应力-应变法的首要误差来源 |
| k 点 | 比总能收敛再加密 1–2 档；金属更密 | C44 / C66 最敏感 |
| `EDIFF` | 1E-7 ~ 1E-8 | 有限差分、能量-应变法都需要 |
| `PREC` / `LREAL` | `Accurate` / `.FALSE.` | 减小力和应力的噪声 |
| `IBRION` / `ISIF` | 6 / 3 | 弹性张量 + 内应变；ISIF=2 只得 Γ 声子 |
| `NFREE` / `POTIM` | 2（或 4）/ 0.015 Å | 原子位移部分；4 更准、代价翻倍 |
| 参考结构 | `ISIF = 3` 弛豫至 \|P\| < 1 kB，并重启一次 | 残余应力会污染 C_ij |
| `SIGMA`（金属） | 0.05–0.2，测试 | 同时检查 entropy T*S < 1 meV/atom |
| 能量-应变幅度 | ±0.5% ~ ±1.5%（≤ 2%） | 太小噪声大，太大非谐 |
| 2D 材料 | 面内应变 + 固定 L_z；单位 N/m | 不要对含真空的 6×6 求逆 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决办法 |
|---|---|---|
| `grep "TOTAL ELASTIC MODULI" OUTCAR` 无输出 | `ISIF < 3`，或作业未正常结束（超时） | 设 `ISIF = 3`；检查 OUTCAR 末尾；6.6 可用 `CHECKPOINT_FD = CONTINUE` 续算 |
| 立方晶体出现 C11 ≠ C22 或非零的 C14、C16 | 晶格没有按标准取向放置，或结构有对称性噪声 | 用 vaspkit 603/604、pymatgen 标准化结构；检查 `SYMPREC` |
| 六方晶体"C44"与文献相差很大，且恰好约等于 (C11−C12)/2 | 把 VASP 的 XY 行当成了 Voigt 第 4 行 | 按 §3.6 重排：VASP (XX,YY,ZZ,XY,YZ,ZX) → Voigt (1,2,3,6,4,5) |
| 数值大 10 倍 | 忘了 kBar → GPa | 除以 10 |
| 出现负本征值，Born 判据不满足 | 参考结构未弛豫、ENCUT/k 点不足；或真实的力学不稳定 | 先检查残余压力、提高精度重算；确认后再讨论物理不稳定性 |
| 弛豫离子与钳位离子结果差别很大 | 内应变贡献本身较大（如 Si 的 C44），或离子 Hessian 有软模/虚频 | 检查 Γ 点频率（`grep THz OUTCAR`）中是否有 `f/i`；虚频会使 Hessian 求逆失真 |
| VASP 报错找不到原网格中的某些 k 点（三方晶系 + `IBRION = 6` + `ISIF = 3`） | 已知问题 7：畸变后低对称网格与原网格不匹配，6.4.3 修复 | 升级到 ≥ 6.4.3；旧版本可设 `KBLOWUP = .FALSE.` 绕过 |
| 结果随 POTIM 明显变化 | 位移太大超出谐振区或太小噪声大 | 保持 0.01–0.02 Å，用 `NFREE = 4` 比较 |
| 二维材料的 C33、C44 接近 0，VRH 结果荒谬 | 真空方向没有刚度，矩阵近奇异 | 只分析面内 3×3 块并换算为 N/m（§7.5） |
| 在 `IBRION = 6` 中使用 Selective dynamics 无效 | 只有 `IBRION = 5` 支持选择性位移 | 去掉 Selective dynamics，或改用 `IBRION = 5`（不利用对称性，代价更高） |
| `IBRION = 8` 算不出完整弹性张量 | DFPT 未实现对应变的微扰，只给离子弛豫贡献 | 弹性常数用 `IBRION = 6` |
| 杂化泛函计算中没有应力输出 | `LHFCALC = .TRUE.` 时 `ISIF` 默认为 0 | 显式设置 `ISIF = 3` |

---

## 练习

1. 对 fcc Al 分别用 ENCUT = 300、400、500 eV 和 k 点 15³、21³、27³ 计算弹性常数，画出 C11、C12、C44 的收敛曲线，判断哪个参数更关键。
2. 用金刚石结构 Si（2 原子原胞）做 `IBRION = 6` + `ISIF = 3`，比较 `SYMMETRIZED`（钳位离子）与 `TOTAL`（弛豫离子）的 C44，解释差别的来源。
3. 对纤锌矿 ZnO 计算弹性张量，故意不重排直接读取，看看"C44"变成了什么；再正确重排，验证 C66 = (C11 − C12)/2。
4. 用 VASPKIT 201 对 Al 做能量-应变法，与 `IBRION = 6` 的结果比较，差异是否在 5% 以内？
5. 对单层 MoS2 运行 §7 的流程，验证 C66 ≈ (C11 − C12)/2，并把结果交给 ElATools 画出 Y(θ) 的极坐标图。

---

## 参考资料

**官方文档**

- [VASP Wiki: IBRION](https://www.vasp.at/wiki/index.php/IBRION) — 有限差分（5/6）与 DFPT（7/8），`ISIF ≥ 3` 时计算弹性与内应变张量
- [VASP Wiki: Phonons from finite differences](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences) — `SYMMETRIZED` / `CONTR FROM IONIC RELAXATION` / `TOTAL ELASTIC MODULI` 三组输出与 ENCUT 收敛建议
- [VASP Wiki: Phonons from density-functional-perturbation theory](https://www.vasp.at/wiki/index.php/Phonons_from_density-functional-perturbation_theory) — DFPT 不计算钳位离子弹性张量的说明
- [VASP Wiki: ISIF](https://www.vasp.at/wiki/index.php/ISIF)、[NFREE](https://www.vasp.at/wiki/index.php/NFREE)、[POTIM](https://www.vasp.at/wiki/index.php/POTIM) — 应力符号约定、位移数与步长
- [VASP Wiki: CHECKPOINT_FD](https://www.vasp.at/wiki/index.php/CHECKPOINT_FD) — 6.6 起有限差分计算的续算与拆分
- [VASP Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — issue 7（三方晶系 `IBRION = 6` + `ISIF = 3`）
- [py4vasp: elastic_modulus](https://vasp.at/py4vasp/latest/calculation/elastic_modulus/) — 从 vaspout.h5 读取弹性张量
- [VASPKIT 教程（弹性常数：能量-应变法）](https://vaspkit.com/tutorials.html) — 任务 201/202/203 与 VPKIT.in 格式
- [pymatgen.analysis.elasticity](https://pymatgen.org/pymatgen.analysis.elasticity.html) — `ElasticTensor`、`DeformedStructureSet`

**相关 GitHub 仓库（仅链接）**

- [keeeto/VASP-Elastic](https://github.com/keeeto/VASP-Elastic) — 从 OUTCAR 提取弹性张量并计算 B、G，演示了 VASP 分量的重排
- [shahramyalameha/ElATools](https://github.com/shahramyalameha/ElATools) — 2D/3D 弹性各向异性分析与 3D 可视化（文档：<https://yalameha.gitlab.io/elastictools/>）
- [materialsproject/atomate2](https://github.com/materialsproject/atomate2) — 含弹性常数工作流（`ElasticMaker`）的高通量框架
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — `Outcar.read_elastic_tensor`、`ElasticTensor` 等
- [vaspkit/VASPKIT_NatureProtocols](https://github.com/vaspkit/VASPKIT_NatureProtocols) — VASPKIT 论文配套示例（含弹性常数输入输出）
- [Bit-Part-Young/scripts](https://github.com/Bit-Part-Young/scripts) — 中文注释的脚本合集，含弹性常数相关后处理脚本
- [haiw201/dft-tutorial](https://github.com/haiw201/dft-tutorial) — 中文 DFT+VASP 课程，可对照其力学性质章节

**关键文献**

- Y. Le Page, P. Saxe, *Phys. Rev. B* 65, 104104 (2002)（应力-应变法最小二乘提取弹性常数）. DOI: [10.1103/PhysRevB.65.104104](https://doi.org/10.1103/PhysRevB.65.104104)
- X. Wu, D. Vanderbilt, D. R. Hamann, *Phys. Rev. B* 72, 035105 (2005)（内应变与离子弛豫贡献）. DOI: [10.1103/PhysRevB.72.035105](https://doi.org/10.1103/PhysRevB.72.035105)
- F. Mouhat, F.-X. Coudert, *Phys. Rev. B* 90, 224104 (2014)（各晶系 Born 稳定性判据）. DOI: [10.1103/PhysRevB.90.224104](https://doi.org/10.1103/PhysRevB.90.224104)
- R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952)（Voigt-Reuss-Hill 平均）. DOI: [10.1088/0370-1298/65/5/307](https://doi.org/10.1088/0370-1298/65/5/307)
- S. F. Pugh, *Phil. Mag.* 45, 823 (1954)（B/G 与延脆性）. DOI: [10.1080/14786440808520496](https://doi.org/10.1080/14786440808520496)
- D. G. Pettifor, *Mater. Sci. Technol.* 8, 345 (1992)（Cauchy 压力与成键）. DOI: [10.1179/mst.1992.8.4.345](https://doi.org/10.1179/mst.1992.8.4.345)
- O. L. Anderson, *J. Phys. Chem. Solids* 24, 909 (1963)（由弹性常数估算德拜温度）. DOI: [10.1016/0022-3697(63)90067-2](https://doi.org/10.1016/0022-3697(63)90067-2)
- S. I. Ranganathan, M. Ostoja-Starzewski, *Phys. Rev. Lett.* 101, 055504 (2008)（普适各向异性指数）. DOI: [10.1103/PhysRevLett.101.055504](https://doi.org/10.1103/PhysRevLett.101.055504)
- M. de Jong et al., *Sci. Data* 2, 150009 (2015)（Materials Project 弹性数据库与工作流）. DOI: [10.1038/sdata.2015.9](https://doi.org/10.1038/sdata.2015.9)
- R. Gaillac, P. Pullumbi, F.-X. Coudert, *J. Phys.: Condens. Matter* 28, 275201 (2016)（ELATE）. DOI: [10.1088/0953-8984/28/27/275201](https://doi.org/10.1088/0953-8984/28/27/275201)
- S. Yalameha, Z. Nourbakhsh, D. Vashaee, *Comput. Phys. Commun.* 271, 108195 (2022)（ElATools）. DOI: [10.1016/j.cpc.2021.108195](https://doi.org/10.1016/j.cpc.2021.108195)
- V. Wang et al., *Comput. Phys. Commun.* 267, 108033 (2021)（VASPKIT）. DOI: [10.1016/j.cpc.2021.108033](https://doi.org/10.1016/j.cpc.2021.108033)
- S. Bertolazzi, J. Brivio, A. Kis, *ACS Nano* 5, 9703 (2011)（单层 MoS2 力学测量）. DOI: [10.1021/nn203879f](https://doi.org/10.1021/nn203879f)

---

[← 上一章](12-phonon.md) | [返回目录](../README.md) | [下一章 →](14-optics-dielectric.md)
