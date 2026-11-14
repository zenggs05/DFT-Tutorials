# 第2章 四大输入文件与输出文件全解析（VASP INCAR, POSCAR, KPOINTS, POTCAR & Output Files Explained）

本章逐行讲解 VASP 的四个输入文件 INCAR、POSCAR、KPOINTS、POTCAR 的语法与常见陷阱，并系统介绍 OUTCAR、OSZICAR、CONTCAR、vasprun.xml、vaspout.h5、CHGCAR、WAVECAR 等输出文件的内容与读取方法（grep / py4vasp / pymatgen / ASE）。最后用 fcc Si 的单点能计算（Si SCF single-point with VASP）把整个流程走一遍，并讲清 ISTART/ICHARG 续算逻辑以及 TOTEN、energy(sigma->0) 等能量的区别。

> [!NOTE]
> 本章所有输入文件都在示例目录 [02-input-files-si](../examples/02-input-files-si/README.md) 中。仓库**不包含 POTCAR**，请按第 6 节用自己的赝势库生成。

---

## 1. 本章目标与适用场景

| 目标 | 你将学会 |
|---|---|
| INCAR | `TAG = value` 语法、注释、`N*value`、续行、Tab/CRLF 陷阱、常用标签分组 |
| POSCAR | 每一行的含义、缩放因子（含负值=体积）、Selective dynamics、Direct/Cartesian、精度与对称性 |
| KPOINTS | Γ 中心 / Monkhorst-Pack 网格、显式列表、Line-mode、KSPACING/KGAMMA |
| POTCAR | 推荐后缀（_pv/_sv/_d/_GW）、头部字段 TITEL/ZVAL/ENMAX、按 POSCAR 顺序拼接 |
| 输出文件 | 每个文件写什么、什么时候写、用哪些命令检查 |
| 续算 | ISTART、ICHARG 的组合与典型场景 |
| 第一个算例 | Si 单点能：准备 → 运行 → 判断是否成功 → 读取结果 |

适用读者：已完成 [第1章](01-setup-compile.md) 的安装，准备开始第一次计算的新手；也适合作为日常查阅的“输入输出速查”。

---

## 2. 原理速览：VASP 需要知道什么

一次 Kohn–Sham DFT 计算需要四类信息，正好对应四个输入文件：

| 需要的信息 | 文件 | 物理含义 |
|---|---|---|
| 原子在哪里、晶胞多大 | POSCAR | 外势 $V_{\text{ext}}(\mathbf{r})$ 的几何 |
| 每种原子的芯区怎么处理 | POTCAR | PAW 赝势：价电子数、投影子、参考能量 |
| 布里渊区怎么采样 | KPOINTS | 把 $\frac{1}{\Omega_{\text{BZ}}}\int d^3k$ 换成 $\sum_{\mathbf{k}} w_{\mathbf{k}}$ |
| 算什么、算多准 | INCAR | 截断能、泛函、收敛判据、离子运动、输出控制 |

平面波基组由截断能 $E_{\text{cut}}$（INCAR 中的 `ENCUT`）决定：所有满足 $\frac{\hbar^2}{2m}|\mathbf{G}+\mathbf{k}|^2 < E_{\text{cut}}$ 的平面波都被纳入。PAW 方法把芯区附近剧烈振荡的波函数交给 POTCAR 中预先构造的投影子与分波处理，因此平面波只需描述较平滑的部分。

计算流程可以概括为：

```text
POSCAR + POTCAR + KPOINTS + INCAR
        |
        v
  初始电荷/波函数 (ISTART, ICHARG)
        |
        v
  电子自洽循环 (ALGO, EDIFF, NELM)  --> OSZICAR 中的 DAV:/RMM: 行
        |
        v
  能量、力、应力  --> 离子步 (IBRION, NSW, EDIFFG) --> CONTCAR, XDATCAR
        |
        v
  OUTCAR, vasprun.xml, vaspout.h5, CHGCAR, WAVECAR, DOSCAR, EIGENVAL ...
```

---

## 3. INCAR：控制参数文件

### 3.1 语法规则

INCAR 是“标签 = 值”格式的纯文本文件。官方规则与实践要点：

| 规则 | 示例 | 说明 |
|---|---|---|
| 每行一个语句 | `ENCUT = 520` | 标签不区分大小写，`=` 两侧空格随意 |
| 一行多个语句用分号 | `ISMEAR = 0 ; SIGMA = 0.05` | — |
| 注释用 `#` 或 `!` | `EDIFF = 1E-6  # SCF criterion` | 注释后的内容被忽略 |
| 重复值简写 `N*value` | `MAGMOM = 2*3.0 4*0` | 等价于 `3.0 3.0 0 0 0 0` |
| 逻辑值 | `.TRUE.` / `.FALSE.`（也接受 `T` / `F`） | — |
| 长行续行用反斜杠 | `MAGMOM = 0 0 1 \`（下一行继续） | 反斜杠后不要留空格 |
| 用双引号包住多行内容 | `WANNIER90_WIN = "...多行..."` | 用于嵌入其他程序的输入 |
| 嵌套标签（较新功能） | `KERNEL_TRUNCATION/LTRUNCATE = T` 或 `PLUGINS { ... }` | 6.x 新功能使用 |

> [!WARNING]
> **拼错的标签会被静默忽略**，VASP 不会报错。例如把 `ISMEAR` 写成 `ISMAER`，计算会用默认的 `ISMEAR = 1`（Methfessel–Paxton），对半导体会给出错误结果。每次计算后都应检查 OUTCAR 中 VASP 实际读到的参数（见 3.3 节）。

> [!CAUTION]
> **Tab 字符和 Windows 换行（CRLF）**是 INCAR 出问题的头号原因：在 Windows 记事本里编辑后上传的 INCAR，可能出现 `Error reading item ... from file INCAR` 或者某些标签被忽略。上传后执行：
>
> ```bash
> dos2unix INCAR POSCAR KPOINTS        # 转换为 Unix 换行
> cat -A INCAR | grep -n -E '\^I|\^M'   # ^I = Tab, ^M = CR；应无输出
> ```

### 3.2 常用标签分组速查

下表列出新手最常用的标签与官方默认值；“推荐起点”是常见计算的初始设置，具体数值需要用第3章的收敛测试确认。

| 类别 | 标签 | 默认值 | 推荐起点 | 作用 |
|---|---|---|---|---|
| 一般 | `SYSTEM` | — | 体系描述 | 写入 OUTCAR 的说明文字 |
| 一般 | `ISTART` | WAVECAR 存在时 1，否则 0 | 0 | 是否读入 WAVECAR |
| 一般 | `ICHARG` | ISTART=0 时 2，否则 0 | 2 | 初始电荷密度来源 |
| 精度 | `PREC` | Normal | Accurate | 决定 FFT 网格密度 |
| 精度 | `ENCUT` | POTCAR 中最大 ENMAX | ≥ 1.3×max(ENMAX)（变胞） | 平面波截断能（eV），**务必手动设置** |
| 精度 | `LREAL` | .FALSE. | 小体系 .FALSE.；> 约 30 原子 Auto | 投影算符在实空间还是倒空间计算 |
| 电子 | `ALGO` | Normal | Normal / Fast | 电子优化算法 |
| 电子 | `EDIFF` | 1E-4 | 1E-6 | 电子自洽收敛判据（eV） |
| 电子 | `NELM` | 60 | 100 | 最大电子步数 |
| 展宽 | `ISMEAR` | 1 | 0 | 占据数展宽方法 |
| 展宽 | `SIGMA` | 0.2 | 0.05（半导体） | 展宽宽度（eV） |
| 自旋 | `ISPIN` | 1 | 磁性体系 2 | 是否自旋极化 |
| 自旋 | `MAGMOM` | 每原子 1.0（ISPIN=2） | 按磁性设置 | 初始磁矩，个数 = 原子数（vasp_ncl 为 3×原子数） |
| 离子 | `IBRION` | NSW=0 或 −1 时为 −1，否则 0 | 优化用 2 | 离子运动算法 |
| 离子 | `NSW` | 0 | 优化 100–200 | 最大离子步数 |
| 离子 | `ISIF` | 2（MD 与杂化泛函为 0） | 2 或 3 | 是否算应力、哪些自由度可变 |
| 离子 | `EDIFFG` | EDIFF×10 | −0.02 | 离子收敛判据；负值表示按力（eV/Å） |
| 输出 | `LWAVE` | .TRUE. | 不需要时 .FALSE. | 写 WAVECAR |
| 输出 | `LCHARG` | .TRUE. | 按需 | 写 CHGCAR 与 CHG |
| 输出 | `LORBIT` | — | DOS/投影 11 | 写 PROCAR 与分轨道 DOSCAR |
| 输出 | `NWRITE` | 2 | 2 | OUTCAR 详细程度 |
| 并行 | `NCORE`、`KPAR`、`NSIM` | 1、1、4 | 见第1章 | 并行与性能 |

两个“小”标签的用法补充：

- `SYSTEM`：任意说明文字，只写入 OUTCAR 等输出文件，便于日后辨认计算内容，不影响结果；
- `NWRITE`：控制 OUTCAR 的详细程度（0–3，默认 2）。官方建议：长时间 MD 用 0 或 1 以控制文件大小；短计算用默认 2；出问题排查时可用 3（每个电子步都写本征值等信息）；4 仅供调试。

示例目录中的 [INCAR_annotated_template](../examples/02-input-files-si/INCAR_annotated_template) 是一份按上述分组整理、带英文注释的 INCAR 模板（输入文件保持纯 ASCII，中文说明都放在本教程里）。它的主体结构如下：

```text
SYSTEM  = annotated_template     # free text, echoed in OUTCAR

#==================== start / restart ====================
ISTART  = 0          # 0 new start; 1 read WAVECAR (default 1 if WAVECAR exists)
ICHARG  = 2          # 2 atomic densities; 1 read CHGCAR; 11 non-SCF with fixed CHGCAR

#==================== basis set and grids ====================
PREC    = Accurate   # Normal (default) or Accurate (denser FFT grids)
ENCUT   = 520        # eV; ALWAYS set explicitly; >= 1.3 x max(ENMAX) for cell relaxations
LREAL   = .FALSE.    # .FALSE. for small cells / accurate energies; Auto for > ~30 atoms

#==================== electronic minimisation ====================
ALGO    = Normal     # Normal (Davidson, default) | Fast (Davidson + RMM-DIIS) | All | Damped
NELM    = 100        # max. SCF steps (default 60)
EDIFF   = 1E-6       # SCF stopping criterion in eV (default 1E-4)

#==================== occupations / smearing ====================
ISMEAR  = 0          # 0 Gaussian | 1,2 Methfessel-Paxton (metals) | -5 tetrahedron (DOS)
SIGMA   = 0.05       # eV; keep entropy T*S below ~1 meV/atom

#==================== output ====================
LWAVE   = .FALSE.    # WAVECAR (default .TRUE.); large file
LCHARG  = .TRUE.     # CHGCAR and CHG (default .TRUE.)
```

（完整文件还包括自旋、离子、KSPACING、输出与并行等被注释掉的可选项。）

> [!TIP]
> 官方建议：如果不确定体系是绝缘体还是金属，用 `ISMEAR = 0` 加较小的 `SIGMA`（0.03–0.1）；VASP ≥ 6.4 可加 `EFERMI = MIDGAP`，让费米能级落在带隙中间，结果更稳定。展宽与收敛的详细讨论见 [第3章](03-convergence.md)。

### 3.3 检查 VASP 实际读到了什么

VASP 会把它对 INCAR 的“理解”写进 OUTCAR，这是发现拼写错误和格式问题的最可靠方法：

```bash
grep -E "ENCUT|PREC|ISMEAR|SIGMA|EDIFF |NELM|LREAL|ISPIN" OUTCAR | head -20
grep -A 3 "Startparameter for this run" OUTCAR     # 参数回显段落的开头
```

此外，`vasprun.xml` 的 `<incar>` 段与 `vaspout.h5` 的 `input`、`original` 组也记录了解析后的参数和原始文件。

---

## 4. POSCAR：结构文件

### 4.1 逐行解读

以示例 [Si_scf/POSCAR](../examples/02-input-files-si/Si_scf/POSCAR)（金刚石 Si 原胞）为例：

```text
Si2 diamond, a = 5.43 A
   5.43000000000000
     0.0000000000    0.5000000000    0.5000000000
     0.5000000000    0.0000000000    0.5000000000
     0.5000000000    0.5000000000    0.0000000000
   Si
     2
Direct
     0.0000000000    0.0000000000    0.0000000000
     0.2500000000    0.2500000000    0.2500000000
```

| 行 | 内容 | 说明 |
|---|---|---|
| 1 | 注释 | 任意文字；官方说明最长 40 字符，超出部分被截断 |
| 2 | 缩放因子 $s$ | 晶格矢量与 Cartesian 坐标都乘以它；**负数表示目标晶胞体积（Å³）**；也可以写 3 个数（用法见官方 POSCAR 页面） |
| 3–5 | 晶格矢量 $\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3$ | 每行一个矢量的笛卡尔分量（乘以 $s$ 后单位为 Å）；这里是 fcc 原胞，$a = 5.43$ Å |
| 6 | 元素符号（VASP 5 起的格式） | **顺序必须与 POTCAR 一致**；元素名最多 2 个字符 |
| 7 | 每种元素的原子数 | 总数即 NIONS |
| 8 | 可选：`Selective dynamics` | 只看首字母 S/s |
| 8/9 | 坐标模式 | 首字母为 C/c/K/k 表示 Cartesian，其余一律按 Direct（分数坐标）处理 |
| 之后 | 原子坐标 | 行数 = NIONS；按第 6、7 行的元素顺序排列 |

金刚石结构的两个 Si 原子分别位于 (0,0,0) 和 (1/4,1/4,1/4)（分数坐标），对应 Fd-3m 空间群。

> [!WARNING]
> 第 6 行的元素符号**只是注释**：真正决定每种原子是什么的是 POTCAR 的顺序与第 7 行的原子数。把第 6 行写成 `As Ga` 而 POTCAR 是 Ga+As，VASP 往往不会报错，却会把 Ga 原子当作 As 来算。务必用第 6 节的方法核对。

### 4.2 Direct 与 Cartesian

- **Direct**：$\mathbf{R} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + x_3\mathbf{a}_3$，与晶胞缩放无关，最常用；
- **Cartesian**：$\mathbf{R} = s\,(x_1, x_2, x_3)$，注意坐标也会乘以第 2 行的缩放因子。

用 ASE 互相转换：

```python
from ase.io import read, write
atoms = read("POSCAR")
write("POSCAR_cart", atoms, format="vasp", direct=False)   # Cartesian
write("POSCAR_frac", atoms, format="vasp", direct=True)    # Direct
print(atoms.get_scaled_positions())                         # fractional coordinates
```

### 4.3 Selective dynamics：固定部分原子

表面计算通常固定底部几层。示例 [POSCAR_selective_dynamics_example](../examples/02-input-files-si/POSCAR_selective_dynamics_example) 是 4 层 Cu(111) 板层（fcc Cu，$a$ = 3.61 Å，接近实验值约 3.615 Å；真空层约 16 Å），底部两层固定：

```text
Cu(111) 1x1 4-layer slab, 16 A vacuum
   1.00000000000000
     2.5526554800    0.0000000000    0.0000000000
     1.2763277400    2.2106644900    0.0000000000
     0.0000000000    0.0000000000   22.2527034200
   Cu
     4
Selective dynamics
Direct
     0.0000000000    0.0000000000    0.3595068800   F   F   F
     0.3333333333    0.3333333333    0.4531689600   F   F   F
     0.6666666667    0.6666666667    0.5468310400   T   T   T
     0.0000000000    0.0000000000    0.6404931200   T   T   T
```

- 每个坐标后跟 3 个 `T`/`F`，分别表示该原子沿 3 个方向能否移动；
- 官方说明：**这 3 个标志对应的是沿晶格矢量 $\mathbf{a}_1,\mathbf{a}_2,\mathbf{a}_3$ 的分数坐标方向**，与坐标用 Direct 还是 Cartesian 输入无关。对于非正交晶胞，`T F T` 并不等于“只能沿 x 和 z 移动”；
- 固定原子在 OUTCAR 中仍会显示受力，判断收敛时要排除它们（仓库脚本 [`check_vasp.py`](../scripts/check_vasp.py) 会考虑这一点，用 `python ../scripts/check_vasp.py -h` 查看用法）。

板层建模与层数/真空收敛见 [第9章](09-surface-workfunction.md)；建模工具见公共章节 [第3章 结构准备](../../Common/tutorials/03-structure-modeling.md)。

### 4.4 精度、对称性与 CONTCAR

- 官方建议坐标与晶格参数至少写 **7 位有效数字**。VASP 用 `SYMPREC`（默认 $10^{-5}$）判断对称性，只写 5 位有效数字是临界情况，可能只找到部分对称操作，进而导致 k 点约化错误或 `PRICEL`、`SGRCON` 一类报错；
- 低精度的结构可以用 spglib/pymatgen 对称化后重新写出（见公共章节 [第4章 ASE 与 pymatgen 实战](../../Common/tutorials/04-ase-pymatgen.md)）；
- **CONTCAR 与 POSCAR 格式相同**，优化中断或需要续算时 `cp CONTCAR POSCAR` 即可；MD 的 CONTCAR 还包含速度和预测-校正坐标，若换系综（如 NVT→NPT）续算，要删掉预测-校正坐标块，否则会崩溃。

---

## 5. KPOINTS：k 点采样

### 5.1 自动网格（最常用）

```text
Gamma-centred mesh (use for fcc, hexagonal, 2D and tetrahedron runs)
0              ! 0 = generate the mesh automatically
Gamma          ! first character G/g = Gamma-centred
  8  8  8      ! N1 N2 N3, roughly proportional to |b1| |b2| |b3|
  0  0  0      ! optional extra shift (fractions of a mesh step)
```

| 行 | 含义 |
|---|---|
| 1 | 注释 |
| 2 | `0` 表示自动生成 |
| 3 | 首字母 `G`/`g`：Γ 中心；`M`/`m`：Monkhorst–Pack |
| 4 | 沿三个倒格矢的分割数 $N_1, N_2, N_3$ |
| 5 | 可选的额外平移（一般写 0 0 0 或省略） |

两种网格的差别（官方 KPOINTS 页面）：

$$
\Gamma\text{ 中心: } \mathbf{k} = \sum_i \frac{n_i + s_i}{N_i}\mathbf{b}_i, \qquad
\text{MP: } \mathbf{k} = \sum_i \frac{n_i + s_i + \frac{1 - N_i}{2}}{N_i}\mathbf{b}_i
$$

$N_i$ 为奇数时两者完全相同；$N_i$ 为偶数时 MP 网格相对 Γ 中心网格平移了半个格点。偶数 MP 网格可能破坏晶格对称性：

| 晶格 | 建议 |
|---|---|
| 面心立方（fcc，含金刚石、闪锌矿、岩盐） | **只用 Γ 中心网格** |
| 六方 / 三方（含石墨烯、MoS₂、ZnO、hcp 金属） | Γ 中心网格（或官方说明的面内 1/3 平移网格）；**避免偶数 MP** |
| 二维材料、表面板层 | Γ 中心，真空方向取 1（如 `12 12 1`） |
| 四面体方法（ISMEAR = −5） | Γ 中心，且至少 4 个 k 点 |
| 简单立方、四方、正交等 | Γ 中心与 MP 都可以 |

> [!TIP]
> 分割数的比例应近似为 $N_1:N_2:N_3 \approx |\mathbf{b}_1|:|\mathbf{b}_2|:|\mathbf{b}_3|$，对正交晶胞等价于与晶格常数成反比。官方给出的经验：长度参数 $R_k = N_i |\mathbf{a}_i|$ 从约 10 Å（宽带隙绝缘体）到约 100 Å（d 金属）不等，具体必须做收敛测试（第3章）。

示例目录中的 [KPOINTS_MP](../examples/02-input-files-si/KPOINTS_MP) 演示了 MP 格式（用的是奇数 7×7×7，因此与 Γ 中心网格等价）。

### 5.2 KSPACING：不写 KPOINTS 文件

如果目录中**没有 KPOINTS 文件**，VASP 按 INCAR 中的 `KSPACING`（单位 Å⁻¹）自动生成网格：

$$
N_i = \max\left(1,\ \left\lceil \frac{2\pi |\mathbf{b}_i|}{\text{KSPACING}} \right\rceil \right)
$$

其中 $\mathbf{b}_i$ 满足 $\mathbf{b}_i\cdot\mathbf{a}_j = \delta_{ij}$（不含 2π）。`KGAMMA = .TRUE.`（默认）生成 Γ 中心网格。以 a = 5.43 Å 的 Si 原胞为例（按上式计算，$2\pi|\mathbf{b}_i| = 2.004$ Å⁻¹）：

| KSPACING (Å⁻¹) | 0.5 | 0.4 | 0.3 | 0.25 | 0.2 | 0.15 |
|---|---|---|---|---|---|---|
| 网格 | 5×5×5 | 6×6×6 | 7×7×7 | 9×9×9 | 11×11×11 | 14×14×14 |

> [!WARNING]
> - **KSPACING 的默认值是 0.5**：忘记放 KPOINTS 又没写 KSPACING 时，VASP 会悄悄用 0.5 Å⁻¹，对金属和小晶胞往往太粗；
> - KPOINTS 文件存在时，`KSPACING`、`KGAMMA` 被忽略；
> - 官方 Known issues 中记录过：支持 KPOINTS_OPT/QPOINTS 的版本在没有 KPOINTS 文件时可能出现 HDF5 内部错误（issue 46），遇到时补一个 KPOINTS 文件即可；
> - KSPACING 生成的网格随晶格常数变化，做状态方程或比较不同体积时，最好固定显式网格。

### 5.3 显式列表

```text
Explicit list: Chadi-Cohen 2-point set for fcc (format demo)
2              ! number of k points
Cartesian      ! units of 2*pi/a, a = scaling factor (line 2) of POSCAR
  0.25  0.25  0.25   1     ! kx ky kz weight (relative weights)
  0.75  0.25  0.25   3
```

- 第 2 行为 k 点个数，第 3 行首字母 C/c/K/k 表示 Cartesian（单位 $2\pi/a$，$a$ 为 POSCAR 缩放因子），其余字符表示倒格矢分数坐标；
- 每行 3 个坐标 + 1 个权重，权重只看相对值，VASP 会自动归一化；
- 上例是文献中常用的 fcc Chadi–Cohen 双点集（Phys. Rev. B 8, 5747 (1973)），这里**只用来演示格式**，实际计算请用收敛的自动网格；
- 显式列表的主要用途：有效质量的有限差分点、杂化泛函能带的“零权重 k 点”方法（[第6章](06-hybrid-metagga.md)）。自动网格生成的 `IBZKPT` 文件就是这种格式，可以直接复制为 KPOINTS 复用。

### 5.4 Line-mode：能带路径

```text
fcc path G-X-W-K-G-L-U-W-L-K|U-X (primitive cell of Si_scf/POSCAR)
 40            ! points per segment (end points included)
Line-mode      ! first character L/l = line mode
Reciprocal     ! coordinates in units of b1 b2 b3
  0.000  0.000  0.000   GAMMA
  0.500  0.000  0.500   X

  0.500  0.000  0.500   X
  0.500  0.250  0.750   W
...
```

- 每两行定义一段路径，段与段之间的空行和点的标签都不是必需的，但 py4vasp 等工具会读取标签；
- 上例的坐标对应 [Si_scf/POSCAR](../examples/02-input-files-si/Si_scf/POSCAR) 的晶格矢量取法（$\mathbf{a}_1 = \frac{a}{2}(0,1,1)$ 等）。例如 X = (0.5, 0, 0.5) 换算到笛卡尔坐标是 $\frac{2\pi}{a}(0,1,0)$，K = (0.375, 0.375, 0.75) 对应 $\frac{2\pi}{a}(\frac34,\frac34,0)$。**换了原胞取法，分数坐标也要跟着换**；
- Line-mode **不能用于自洽计算**：必须先用均匀网格自洽，再 `ICHARG = 11` 固定电荷做非自洽计算；VASP ≥ 6.3 还可以用 `KPOINTS_OPT` 一步完成（详见 [第5章](05-band-dos.md)）；meta-GGA 和杂化泛函的能带需要特殊处理（[第6章](06-hybrid-metagga.md)）。

完整路径见示例 [KPOINTS_line](../examples/02-input-files-si/KPOINTS_line)。

> [!NOTE]
> 旧的“全自动”模式（第 3 行写 `Auto`，第 4 行给长度 $R_k$）已被官方标为不推荐，请改用 `KSPACING` 或显式网格。

---

## 6. POTCAR：PAW 赝势文件

### 6.1 赝势库版本与后缀

VASP 的 PAW 赝势按“泛函 + 发布版本”组织，常见目录名如 `potpaw_PBE.64`（最新，推荐）、`potpaw_PBE.54`（VASP 5.4 时代）、`potpaw_LDA.64`。PBE 赝势用于 PBE 以及其他 GGA、meta-GGA、杂化泛函；只有 LDA 计算才用 LDA 赝势。同一元素有多种变体，以后缀区分：

| 后缀 | 含义 | 典型用法 |
|---|---|---|
| 无 | 标准价电子划分 | 多数主族元素 |
| `_pv` / `_sv` | 半芯 p（/s+p）态作为价电子 | 过渡金属、碱金属、碱土金属；磁性、高压 |
| `_d` | 半芯 d 态作为价电子 | Ga、Ge、In、Sn、Pb、Bi 等 |
| `_s` / `_h` | 软 / 硬 | `_s` 只做粗筛，**不要**用于杂化泛函或换泛函的计算 |
| `_GW` | 优化了高能散射性质 | GW、BSE、光学 |
| `_2`、`_3` | 镧系 f 电子冻结在芯内 | 不关心 4f 的镧系化合物 |

官方推荐列表在 **2026 年 10 月**有更新（标准 DFT）：`Be→Be_sv`、`Mg→Mg_pv`、`Fe→Fe_pv`、`Co→Co_pv`、`Ni→Ni_pv`、`Cu→Cu_pv`、`Eu_2→Eu_3`、`Yb_2→Yb_3`、`Re→Re_pv`、`Os→Os_pv`。完整的元素–后缀–ZVAL–ENMAX 速查表见 [potcar_recommended_list.md](../examples/02-input-files-si/potcar_recommended_list.md)。

> [!IMPORTANT]
> 复现旧文献时要用**当时**的赝势选择；同一项研究中所有需要比较能量的计算必须使用同一组 POTCAR 和同一个 ENCUT。论文中写明 `grep TITEL POTCAR` 的输出。

### 6.2 头部关键字段

POTCAR 是只读文件，**不要编辑**（换泛函请用 INCAR 的 `GGA`、`METAGGA`、`XC` 标签）。用 grep 查看头部信息即可：

```bash
grep -E "TITEL|ZVAL|ENMAX|VRHFIN|LEXCH" POTCAR
```

| 字段 | 含义 | 用途 |
|---|---|---|
| `TITEL` | 赝势名称与生成日期，如 `PAW_PBE Si 05Jan2001` | 核对元素顺序、写进论文 |
| `VRHFIN` | 价电子组态 | 判断是否包含半芯态 |
| `LEXCH` | 生成赝势所用泛函（`PE` = PBE） | 确认赝势集 |
| `ZVAL` | 价电子数 | 计算 NELECT、Bader 净电荷 |
| `POMASS` | 原子质量 | MD |
| `ENMAX` / `ENMIN` | 推荐 / 最低截断能（eV） | 设定 ENCUT 的依据 |

以 Si 为例：`PAW_PBE Si`，ZVAL = 4，ENMAX = 245.345 eV（官方 Available pseudopotentials 页面）。

### 6.3 按 POSCAR 顺序拼接

POTCAR 由各元素的 POTCAR 按 **POSCAR 第 6 行的顺序**首尾相接而成：

```bash
# 假设赝势库目录结构为 $VASP_PP_PATH/potpaw_PBE/<变体名>/POTCAR（ASE 也使用这一约定）
export VASP_PP_PATH=$HOME/vasp/potentials
cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR    # GaAs：Ga 在前
grep TITEL POTCAR                                         # 应依次为 Ga_d、As
sed -n '6p' POSCAR                                        # 应为 Ga As
```

示例目录提供了 [make_potcar.sh](../examples/02-input-files-si/make_potcar.sh)：读取 POSCAR 第 6 行，按官方 2026-10 推荐后缀自动拼接，并打印每个元素的 TITEL、ZVAL、ENMAX 以及 1.3×max(ENMAX)：

```bash
bash make_potcar.sh                     # 当前目录的 POSCAR
bash make_potcar.sh POSCAR Ti=Ti_pv     # 覆盖默认选择
PP_SET=potpaw_PBE.54 bash make_potcar.sh
```

其他工具：VASPKIT 的 `103`（按默认设置生成 POTCAR，需要在 `~/.vaspkit` 中配置赝势路径）；pymatgen 的 `Potcar` 类（需先用 `pmg config` 设置 `PMG_VASP_PSP_DIR`）；仓库脚本 [`vasp_inputgen.py`](../scripts/vasp_inputgen.py) 也能按推荐列表拼接 POTCAR（用 `python ../scripts/vasp_inputgen.py -h` 查看用法）。

> [!CAUTION]
> 仓库中**永远不要**提交 POTCAR（本仓库的 `.gitignore` 已排除 POTCAR、WAVECAR、CHGCAR 等）。`vaspout.h5` 中也含有 POTCAR 文本，同样不能公开分享。

---

## 7. 输出文件

### 7.1 总览

| 文件 | 内容 | 何时写出 / 控制标签 | 常用查看方式 |
|---|---|---|---|
| `OUTCAR` | 最完整的文本输出：参数回显、对称性、k 点、每步能量、力、应力、磁矩、计时 | 总是；详细程度由 `NWRITE` 控制 | `grep` |
| `OSZICAR` | 每个电子步、离子步的能量摘要 | 总是 | `cat`、`tail` |
| 标准输出（如 `vasp.out`） | 与 OSZICAR 类似，另含警告和报错 | 重定向得到 | `less` |
| `CONTCAR` | 最新结构（POSCAR 格式） | 每个离子步更新 | VESTA、ASE |
| `XDATCAR` | 离子步轨迹 | 有离子步时 | ASE、OVITO、py4vasp |
| `vasprun.xml` | XML 格式的完整记录（参数、能量、本征值、DOS、介电函数等） | 总是 | pymatgen、ASE、sumo |
| `vaspout.h5` | HDF5 格式的输入与结果 | 编译时开启 `-DVASP_HDF5` | py4vasp |
| `CHGCAR` / `CHG` | 电荷密度（CHGCAR 含 PAW 单中心占据，可用于续算） | `LCHARG`（默认 .TRUE.） | VESTA、`ICHARG=1/11` |
| `WAVECAR` | 波函数（二进制，通常很大） | `LWAVE`（默认 .TRUE.） | `ISTART=1`；后处理工具 |
| `vaspwave.h5` | HDF5 格式的波函数/电荷 | `LH5 = .TRUE.`（≥ 6.0），此时不再写 WAVECAR/CHGCAR | py4vasp |
| `DOSCAR` | 态密度（`LORBIT=11` 时含分轨道 PDOS） | 总是写，内容取决于 LORBIT | 绘图脚本、sumo、py4vasp |
| `EIGENVAL` | 每个 k 点的本征值与占据 | 总是 | 能带绘图 |
| `PROCAR` | 波函数在原子/轨道上的投影 | `LORBIT ≥ 10` | 胖能带（fat band） |
| `LOCPOT` | 局域静电势 | `LVTOT` 或 `LVHAR` | 功函数（[第9章](09-surface-workfunction.md)） |
| `IBZKPT` | 不可约 k 点及权重 | 使用自动网格时 | 可复制为 KPOINTS |

VASP 6.4.2 及以前的版本中，`LCHARG` 的默认值写作 `.NOT. LH5`；由于 `LH5` 默认为 .FALSE.，效果都是“写 CHGCAR”。

### 7.2 OSZICAR：一眼看懂收敛过程

（示意输出，数值仅供参考）

```text
       N       E                     dE             d eps       ncg     rms          rms(c)
DAV:   1     0.157386453128E+02    0.15739E+02   -0.31012E+03   432   0.794E+02
DAV:   2    -0.103217806325E+02   -0.26060E+02   -0.25226E+02   592   0.118E+02
DAV:   3    -0.108386402911E+02   -0.51686E+00   -0.51122E+00   480   0.174E+01    0.312E+00
...
DAV:  11    -0.108446715391E+02   -0.12345E-06   -0.10012E-07   480   0.121E-03    0.315E-04
   1 F= -.10844672E+02 E0= -.10844672E+02  d E =-.108447E+02
```

| 列 | 含义 |
|---|---|
| `DAV:` / `RMM:` | 当前电子步使用的算法（Davidson / RMM-DIIS） |
| `N` | 电子步序号 |
| `E` | 当前自由能（eV） |
| `dE` | 与上一步自由能之差；与 `d eps` 都小于 EDIFF 时电子步收敛 |
| `d eps` | 能带结构能的变化 |
| `ncg` | 哈密顿量作用次数 |
| `rms`、`rms(c)` | 波函数残差、输入输出电荷密度之差 |
| 离子步行 `F=` | 自由能 TOTEN |
| 离子步行 `E0=` | 外推到 σ→0 的能量 |
| 离子步行 `d E` | 与上一离子步的能量差 |

自旋极化计算的离子步行末尾还有 `mag=`（总磁矩）。判断要点：电子步数应明显小于 `NELM`；若某个离子步的电子步数恰好等于 NELM，说明该步**没有收敛**。

### 7.3 OUTCAR：最常用的 grep 命令

```bash
grep "General timing and accounting" OUTCAR      # 出现 = 正常结束
grep -E "free +energy +TOTEN" OUTCAR | tail -1   # 最终自由能 F
grep "energy(sigma->0)" OUTCAR | tail -1         # E 与 E(sigma->0)
grep "EENTRO" OUTCAR | tail -1                   # 熵项 T*S
grep "E-fermi" OUTCAR | tail -1                  # 费米能级
grep "NKPTS" OUTCAR                              # 不可约 k 点数与 NBANDS
grep "NELECT" OUTCAR                             # 电子数
grep "reached required accuracy" OUTCAR          # 离子优化达到 EDIFFG
N=2; grep -A $((N+1)) "TOTAL-FORCE" OUTCAR | tail -n $N   # 最后一步每个原子的坐标与受力（N = 原子数）
grep "total drift" OUTCAR | tail -1              # 力的漂移，应远小于关心的力
grep "external pressure" OUTCAR | tail -1        # 外压与 Pulay stress（kB）
grep -A 3 "in kB" OUTCAR | tail -4               # 应力张量
grep "magnetization" OUTCAR | tail -3            # 磁矩（ISPIN=2）
grep -E "Elapsed time|Maximum memory used" OUTCAR
grep "LOOP+" OUTCAR | tail -3                    # 每个离子步的耗时
grep -n -i -E "warning|advice" OUTCAR | head     # VASP 的提示（6.x 中会用字符画框出）
```

OUTCAR 中的最终能量块（示意输出，数值仅供参考）：

```text
  free  energy   TOTEN  =       -10.84467200 eV

  energy  without entropy=      -10.84467200  energy(sigma->0) =      -10.84467200
```

### 7.4 三种能量的区别

使用展宽（ISMEAR ≥ −1）时，占据数不再是 0/1，VASP 实际最小化的是广义自由能 $F = E - TS$：

| 名称 | 出现位置 | 定义 | 何时使用 |
|---|---|---|---|
| `TOTEN`（free energy，F） | OUTCAR、OSZICAR 的 `F=` | $F = E - TS$ | 与受力、应力严格对应；MP 展宽下比较能量时常用 |
| `energy without entropy`（E） | OUTCAR | $E = F + TS$ | 一般不直接用 |
| `energy(sigma->0)`（E₀） | OUTCAR、OSZICAR 的 `E0=` | 官方给出的外推式 $E_0 = \frac{1}{2}(F + E)$ | Gaussian 展宽下比较能量的首选 |

- 半导体用小 σ 的 Gaussian 展宽或四面体方法时，$TS \approx 0$，三者几乎相等；
- 金属中 $TS$ 不可忽略，应检查 `EENTRO`：**每个原子的 `|T*S|` 应小于约 1 meV**，否则减小 SIGMA（[第3章](03-convergence.md)）；
- 力和应力对应的是 $F$，而不是 $E_0$。

### 7.5 用 Python 读取结果

```python
# py4vasp (needs a VASP build with -DVASP_HDF5)
from py4vasp import Calculation
calc = Calculation.from_path(".")
calc.energy.print()                 # energies of the last ionic step
print(calc.energy.read())           # same data as a dict
print(calc.structure.read())        # final structure

# pymatgen
from pymatgen.io.vasp import Vasprun
vr = Vasprun("vasprun.xml")
print(vr.final_energy, vr.converged)

# ASE
from ase.io import read
atoms = read("OUTCAR")              # or read("vasprun.xml")
print(atoms.get_potential_energy())                         # energy(sigma->0)
print(atoms.get_potential_energy(force_consistent=True))    # free energy TOTEN
print(atoms.get_forces())
```

更多后处理工具（VASPKIT、sumo、pyprocar）见公共章节 [第5章 后处理与可视化](../../Common/tutorials/05-postprocess-visualization.md) 与 [第18章](18-postprocess-troubleshooting.md)。

---

## 8. 续算逻辑：ISTART 与 ICHARG

### 8.1 选项含义

| ISTART | 含义 |
|---|---|
| 0 | 从头开始，不读 WAVECAR |
| 1 | 读 WAVECAR，按新的截断能/晶胞重新建立平面波基组（“截断能不变”续算）；适用于 ENCUT 测试、换体积的计算 |
| 2 | 读 WAVECAR，保持原来的平面波基组不变（“基组不变”续算）；变胞优化需要与前一次基组一致时使用 |
| 3 | 同 2，另读 TMPCAR（MD 的轨道/电荷预测），一般不推荐 |

| ICHARG | 含义 |
|---|---|
| 0 | 由初始波函数计算电荷密度（ISTART ≠ 0 时的默认值） |
| 1 | 从 CHGCAR 读电荷密度，并按新位置外推 |
| 2 | 原子电荷密度叠加（ISTART = 0 时的默认值） |
| 4 | 从 POT 文件读局域势 |
| 11 | 从 CHGCAR 读电荷且**整个计算中保持不变**（非自洽），用于能带和 DOS |
| 12 | 原子电荷叠加且保持不变（Harris–Foulkes 非自洽） |

> [!NOTE]
> 若 WAVECAR 缺失或与当前体系不兼容（能带数、k 点数不符），VASP 会把 ISTART 自动退回 0，此时 ICHARG 也会变为 2。因此“我明明设了 ISTART=1”并不保证真的读入了波函数——查看 OUTCAR 开头的提示或 OSZICAR 第一步的能量即可判断。

### 8.2 典型场景

| 场景 | 设置 | 注意 |
|---|---|---|
| 新计算 | `ISTART = 0`、`ICHARG = 2` | 目录里不要留旧的 WAVECAR/CHGCAR |
| 结构优化中断后续算 | `cp CONTCAR POSCAR`；保留 WAVECAR 则 `ISTART = 1` | 变胞优化的 Pulay 应力问题见 [第4章](04-relax-eos.md) |
| 能带 / DOS 非自洽 | `ICHARG = 11`，读自洽的 CHGCAR | ENCUT、ISPIN、POTCAR、LMAXMIX 与自洽计算一致；DFT+U 需 `LMAXMIX = 4`（d）/ `6`（f） |
| HSE 从 PBE 波函数出发 | 保留 PBE 的 WAVECAR，`ISTART = 1` | 见 [第6章](06-hybrid-metagga.md) |
| 系列小改动（如换 U 值） | `ICHARG = 1` 读上一个 CHGCAR | 官方建议，可减少电子步 |

> [!WARNING]
> 不同体系、不同原子数的计算**不要**共用目录。旧 WAVECAR/CHGCAR 被意外读入是 `BRMIX: very serious problems` 等报错的常见来源；有疑问时删掉它们从头算。

---

## 9. 第一个算例：fcc Si 单点能

### 9.1 输入文件

示例目录 [Si_scf/](../examples/02-input-files-si/README.md) 中的 INCAR：

```text
SYSTEM = Si bulk - first SCF calculation

# start from scratch
ISTART = 0            # do not read WAVECAR
ICHARG = 2            # superposition of atomic charge densities

# basis set and accuracy
PREC   = Accurate     # dense FFT grids, no wrap-around errors
ENCUT  = 320          # 1.3 x ENMAX of PAW_PBE Si (245.3 eV), rounded
LREAL  = .FALSE.      # projection in reciprocal space (small cell)

# electronic self-consistency
ALGO   = Normal       # blocked Davidson
EDIFF  = 1E-6         # eV
NELM   = 100

# occupations of a semiconductor
ISMEAR = 0            # Gaussian smearing
SIGMA  = 0.05         # eV

# static calculation (no ionic steps)
IBRION = -1
NSW    = 0

# output control
LWAVE  = .FALSE.      # no WAVECAR
LCHARG = .TRUE.       # write CHGCAR (input for a later band structure)
LORBIT = 11           # write PROCAR and lm-projected DOSCAR
```

| 关键词 | 取值 | 为什么这样选 |
|---|---|---|
| `PREC` | Accurate | FFT 网格更密，避免混叠（wrap-around）误差；内存略增 |
| `ENCUT` | 320 | Si 的 ENMAX = 245.3 eV，1.3 倍约 319 eV；这一档在后续需要计算应力/变胞时也够用（第3章验证） |
| `LREAL` | .FALSE. | 2 原子小晶胞，倒空间投影既精确又不慢 |
| `EDIFF` | 1E-6 | 官方推荐的“收敛良好”电子判据 |
| `ISMEAR`/`SIGMA` | 0 / 0.05 | 半导体用 Gaussian 小展宽；**不要用默认的 ISMEAR = 1** |
| `IBRION`/`NSW` | −1 / 0 | 单点能，不移动原子 |
| `LCHARG` | .TRUE. | 保留 CHGCAR，第5章做能带时可直接用 |
| `LORBIT` | 11 | 顺便输出 PROCAR 与分轨道 DOSCAR，便于熟悉文件格式 |

KPOINTS 为 Γ 中心 8×8×8（fcc 必须用 Γ 中心），POSCAR 见 4.1 节。

### 9.2 运行

```bash
EX=$HOME/DFT/VASP/examples/02-input-files-si      # 本仓库示例目录，按实际路径修改
mkdir -p ~/calc && cp -r $EX/Si_scf ~/calc/Si_scf && cd ~/calc/Si_scf
bash $EX/make_potcar.sh                           # 或 cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
bash $EX/check_inputs.sh                          # 运行前自检（只读）
mpirun -np 4 vasp_std > vasp.out 2>&1
```

[check_inputs.sh](../examples/02-input-files-si/check_inputs.sh) 会检查：文件是否齐全、Tab/CRLF、POSCAR 坐标行数与 Selective dynamics 标志、POTCAR 元素顺序是否与 POSCAR 一致、ENCUT 与 max(ENMAX) 之比、MAGMOM 个数、四面体方法的 k 点数与偶数 MP 网格等。

### 9.3 检查结果

```bash
tail -3 OSZICAR                                   # 电子步数应远小于 100
grep "General timing and accounting" OUTCAR       # 正常结束
grep "NKPTS" OUTCAR                               # 8x8x8 Gamma 网格在 Oh 对称性下约化为 29 个 k 点
grep "NELECT" OUTCAR                              # 2 x ZVAL(4) = 8
grep "energy(sigma->0)" OUTCAR | tail -1          # 最终能量
grep "EENTRO" OUTCAR | tail -1                    # 半导体应接近 0
grep "E-fermi" OUTCAR | tail -1
ls -lh CHGCAR DOSCAR EIGENVAL PROCAR vasprun.xml  # 输出文件
```

| 检查项 | 预期 |
|---|---|
| 不可约 k 点数 | 29（按 fcc 的 48 个点群操作约化 8×8×8 Γ 网格的结果；以 OUTCAR 为准） |
| NELECT | 8 |
| NBANDS | 默认公式 $\max\left(\frac{1}{2}\text{nint}(N_e+2) + \max(\frac{N_{\text{ions}}}{2}, 3),\ \text{int}(\frac{3}{5}N_e)\right) = 8$；并行进程数不能整除时会被自动增大 |
| 电子步数 | 十几步（示意，取决于版本与并行） |
| 总能 | 约 −10.8 eV（2 原子，经验量级，仅供参考） |
| EENTRO | 接近 0（Si 有带隙） |

如果 `energy(sigma->0)` 与 `TOTEN` 相差明显、或 EENTRO 不接近零，通常说明 ISMEAR/SIGMA 没设对，或者体系实际上是金属。

### 9.4 用 py4vasp 快速看一眼

```python
from py4vasp import Calculation
calc = Calculation.from_path(".")
calc.energy.print()
calc.dos.plot()          # LORBIT = 11 gives projected DOS; opens an interactive figure
```

能带与 DOS 的正规计算流程见 [第5章](05-band-dos.md)；ENCUT 和 k 点是否收敛见 [第3章](03-convergence.md)。

---

## 参数选择与收敛建议

| 项目 | 建议 | 说明 |
|---|---|---|
| INCAR 格式 | 纯 ASCII、LF 换行、无 Tab | `dos2unix`、`cat -A` 检查 |
| ENCUT | 始终手动设置；单点 ≥ max(ENMAX)，变胞 ≥ 1.3×max(ENMAX) | 同一研究统一 |
| PREC | Accurate（小体系/高精度），Normal（大体系初筛） | — |
| EDIFF | 1E-6（一般），1E-7–1E-8（声子、弹性） | — |
| ISMEAR/SIGMA | 半导体 0/0.05；金属优化 1/0.1–0.2；DOS −5 | 检查 EENTRO |
| POSCAR 精度 | ≥ 7 位有效数字 | 对称性识别 |
| KPOINTS | fcc/六方/二维/四面体用 Γ 中心 | 避免偶数 MP |
| KSPACING | 仅在无 KPOINTS 时生效，默认 0.5 | 生产计算建议显式网格 |
| POTCAR | 官方推荐后缀（2026-10 起）；顺序与 POSCAR 一致 | `grep TITEL POTCAR` |
| 输出控制 | 不需要时 `LWAVE = .FALSE.` | 节省磁盘 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决办法 |
|---|---|---|
| `Error reading item 'MAGMOM' from file INCAR.` | MAGMOM 个数与原子数不符（非共线时每原子 3 个分量），或含非法字符 | 用 `N*value` 写法核对总数；vasp_ncl 写 3×NIONS 个数；`check_inputs.sh` 会检查 |
| 标签“没起作用”，结果像用了默认值 | 标签拼错、Tab、CRLF | `grep` OUTCAR 中的回显；`dos2unix INCAR` |
| 计算正常结束但结果明显不对（晶格常数、能量离谱） | POTCAR 元素顺序与 POSCAR 不一致 | `grep TITEL POTCAR` 对照 `sed -n 6p POSCAR` |
| 能量、费米能级和预期差很多，且 k 点很少 | 没有 KPOINTS 文件，VASP 用了默认 KSPACING = 0.5 | 放入 KPOINTS 或显式设置 KSPACING |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六方/fcc 等晶格用了偶数 Monkhorst–Pack 网格 | 改为 Γ 中心网格 |
| 报错中含 `Tetrahedron method fails` | ISMEAR = −5 但 k 点少于 4 个（如 Γ 点计算） | 改用 `ISMEAR = 0`、`SIGMA = 0.05` |
| `internal error in subroutine PRICEL` / `SGRCON` 相关报错 | 坐标精度不足导致对称性判断不一致 | 对称化结构、写足有效数字，必要时调整 `SYMPREC` 或 `ISYM = 0` |
| `The distance between some ions is very small` | 原子重叠：坐标单位/模式写错（Direct 与 Cartesian 混用）、缩放因子错误 | 用 VESTA/ASE 可视化 POSCAR，检查第 2 行与坐标模式行 |
| `WARNING: dimensions on CHGCAR file are different` | 读入的 CHGCAR 来自 FFT 网格不同的计算 | ICHARG=11 时保持 ENCUT、PREC 与自洽计算一致 |
| `BRMIX: very serious problems the old and the new charge density differ` | 读入了不匹配的 CHGCAR/WAVECAR，或对称性判断不一致 | 删除旧 CHGCAR/WAVECAR 重算；必要时调整 SYMPREC 或 `ISYM = 0` |
| 磁盘被占满 | 每个计算都写了 WAVECAR/CHG | 不需要时 `LWAVE = .FALSE.`；定期清理 |

---

## 练习

1. 把 Si_scf 的 POSCAR 改写为 Cartesian 坐标，再改为“缩放因子 = 1、晶格矢量直接写 Å”的形式，确认三种写法得到相同能量。
2. 把 POSCAR 第 2 行改成负数（等于原胞体积的相反数，$V = a^3/4$），确认结果不变。
3. 删除 KPOINTS，在 INCAR 中分别设置 `KSPACING = 0.5, 0.3, 0.2`，用 `grep NKPTS OUTCAR` 与 `head IBZKPT` 查看网格，与 5.2 节的表格对照。
4. 故意把 INCAR 中的 `ISMEAR` 拼成 `ISMAER`，运行后在 OUTCAR 中找到证据说明 VASP 使用了默认值。
5. 用 GaAs（闪锌矿，a ≈ 5.65 Å）建立 POSCAR，用 `make_potcar.sh` 生成 POTCAR，再故意把 POTCAR 顺序调反，用 `check_inputs.sh` 检查。

---

## 参考资料

**官方文档**

- INCAR：https://www.vasp.at/wiki/index.php/INCAR
- POSCAR：https://www.vasp.at/wiki/index.php/POSCAR
- KPOINTS：https://www.vasp.at/wiki/index.php/KPOINTS ；KSPACING：https://www.vasp.at/wiki/index.php/KSPACING
- POTCAR：https://www.vasp.at/wiki/index.php/POTCAR ；Choosing pseudopotentials：https://www.vasp.at/wiki/index.php/Choosing_pseudopotentials ；Available pseudopotentials：https://www.vasp.at/wiki/index.php/Available_PAW_potentials
- OUTCAR / OSZICAR / CONTCAR / vaspout.h5：https://www.vasp.at/wiki/index.php/OUTCAR ，https://www.vasp.at/wiki/index.php/OSZICAR ，https://www.vasp.at/wiki/index.php/CONTCAR ，https://www.vasp.at/wiki/index.php/Vaspout.h5
- ISTART / ICHARG：https://www.vasp.at/wiki/index.php/ISTART ，https://www.vasp.at/wiki/index.php/ICHARG
- Smearing technique 与 K-point integration：https://www.vasp.at/wiki/index.php/Smearing_technique
- py4vasp：https://vasp.at/py4vasp/latest/ ；pymatgen.io.vasp：https://pymatgen.org/pymatgen.io.vasp.html ；ASE VASP 接口：https://docs.ase-lib.org/ase/calculators/vasp.html

**相关 GitHub 仓库**（仅链接）

- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：读写 INCAR/POSCAR/KPOINTS/Potcar，解析 Vasprun/Outcar。
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：官方 HDF5 后处理库。
- [PytLab/VASPy](https://github.com/PytLab/VASPy)：轻量级 VASP 文件读写库，可参考其解析思路。
- [bzkarimi/VASP](https://github.com/bzkarimi/VASP)：INCAR 常用参数说明，可用来交叉核对参数速查表。
- [cpfpengfei/VASP-Examples](https://github.com/cpfpengfei/VASP-Examples)：初学者笔记，含 Si 的输入文件要点。
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book)：经典的 VASP + ASE 英文教程书。
- [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script)：中文社区常用脚本合集（含 INCAR 转 Unix 格式等）。

**关键文献**

- H. J. Monkhorst, J. D. Pack, Phys. Rev. B 13, 5188 (1976). DOI: 10.1103/PhysRevB.13.5188
- D. J. Chadi, M. L. Cohen, Phys. Rev. B 8, 5747 (1973). DOI: 10.1103/PhysRevB.8.5747
- G. Kresse, D. Joubert, Phys. Rev. B 59, 1758 (1999). DOI: 10.1103/PhysRevB.59.1758
- P. E. Blöchl, Phys. Rev. B 50, 17953 (1994). DOI: 10.1103/PhysRevB.50.17953

[← 上一章](01-setup-compile.md) | [返回目录](../README.md) | [下一章 →](03-convergence.md)
