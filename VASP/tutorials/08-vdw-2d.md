# 第8章 范德华修正、二维与层状材料（van der Waals Corrections, 2D & Layered Materials with VASP）

本章讲解 VASP 中范德华（vdW，色散）修正的全部主流方案——DFT-D2/D3/D4、simple-DFT-D3（IVDW=15，VASP 6.6+）、Tkatchenko–Scheffler 与 MBD/libMBD、以及 optB86b-vdW、rVV10、BEEF-vdW 等非局域 vdW-DF 泛函——并以石墨、单层/双层 MoS2 为例，演示层间距、层间结合能、二维材料真空层与 k 点设置、真空能级对齐和异质结晶格匹配（VASP vdW correction / DFT-D3 / 2D materials / layered materials tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说明为什么 PBE 等半局域泛函描述不了层状材料的层间结合；
- 根据 VASP 版本和编译选项，选择合适的 `IVDW`（11/12/13/15/20/21/202/14 …）或非局域泛函（`LUSE_VDW`）；
- 正确写出 optB86b-vdW、rVV10、r2SCAN+rVV10、BEEF-vdW 的 INCAR，知道 `vdw_kernel.bindat` 在不同版本中的处理方式；
- 计算石墨的平衡层间距与层间结合能，并与 RPA/实验参考值比较；
- 为二维材料设置真空层、k 点（z 方向取 1）、`LATTICE_CONSTRAINTS`，做真空收敛并把带边对齐到真空能级；
- 扫描双层 MoS2 的结合能曲线，初步掌握异质结的晶格失配与超胞匹配。

| 体系 | vdW 是否关键 | 推荐起点 | 本章小节 |
|---|---|---|---|
| 石墨、h-BN、MoS2 等层状体相 | 决定层间距与结合能 | D3(BJ)、optB86b-vdW、rev-vdW-DF2 | §6、§7 |
| 二维单层（面内性质） | 影响较小，但应与多层保持一致 | 与比较对象相同的 vdW 设置 | §8 |
| 双层、异质结、堆垛能 | 关键 | D3(BJ)/D4，或 vdW-DF 类 | §9、§10 |
| 分子晶体、物理吸附 | 关键 | D3(BJ)、D4、MBD | §3、§5 |
| 金属表面吸附（催化） | 显著 | D3(BJ)、optB86b-vdW、BEEF-vdW | §5，[第10章](10-adsorption-catalysis.md) |

---

## 2. 原理速览

### 2.1 为什么需要 vdW 修正

色散作用来自电子密度的长程关联涨落，两个中性、不重叠的电荷分布之间的吸引能在远距离按 $-C_6/r^6$ 衰减。LDA/GGA 的交换关联能只依赖局部（或局部梯度）的电子密度，**在密度几乎不重叠的区域不可能产生这种 $r^{-6}$ 吸引**。结果：PBE 计算的石墨层间几乎不结合（c 轴过长），分子晶体体积偏大，物理吸附能严重偏小。（LDA 往往"碰巧"给出看似合理的层间距，但那是误差抵消，不可依赖。）

### 2.2 原子对加和型修正（DFT-D 一族、TS）

在 DFT 总能上加一个只依赖原子坐标的经验项：

$$
E_{\rm disp}=-\frac12\sum_{i,j}\sum_{\mathbf L}{}^{\prime}\left(f_{6}(r_{ij,\mathbf L})\frac{C_{6,ij}}{r_{ij,\mathbf L}^{6}}+f_{8}(r_{ij,\mathbf L})\frac{C_{8,ij}}{r_{ij,\mathbf L}^{8}}\right)
$$

阻尼函数 $f_n$ 在短程把修正关掉，避免与 DFT 本身已描述的部分重复。各方法的区别在于 $C_6$ 怎么来：D2 用固定原子参数；D3 依赖配位数；D4 还依赖原子部分电荷；TS 由 Hirshfeld 分区的电子密度得到；MBD 进一步包含多体极化屏蔽。由于修正项是坐标的解析函数，受力和应力都可以直接得到，可用于结构优化、MD 和有限位移声子。

### 2.3 非局域 vdW-DF 泛函

vdW-DF 一族在交换关联泛函中加入一个非局域关联项：

$$
E_c^{\rm nl}=\frac12\iint n(\mathbf r)\,\Phi(\mathbf r,\mathbf r')\,n(\mathbf r')\,d^3r\,d^3r'
$$

核函数 $\Phi$ 依赖两点的密度与梯度。VASP 用 Román-Pérez–Soler 的 FFT 算法计算它，额外耗时主要取决于 FFT 网格大小（即晶胞大小，包括真空），与原子数关系不大。不同 vdW-DF 变体（optB86b、optB88、rev-vdW-DF2 …）的区别主要在于搭配的交换泛函和核函数参数。

### 2.4 层间结合能的定义

对含 $N_L$ 层的体相晶胞（如石墨的 2 层 AB 胞）：

$$
E_b=\frac{E_{\rm mono}}{N_{\rm mono}}-\frac{E_{\rm bulk}}{N_{\rm bulk}}\quad(\text{meV/atom}),\qquad
E_b^{\rm area}=\frac{N_L E_{\rm mono}-E_{\rm bulk}}{N_L A}\quad(\text{meV/Å}^2)
$$

注意区分**层间结合能**（binding，所有层完全分开的平均）与**剥离能**（exfoliation，从体相表面剥下一层）；两者数值接近但定义不同。单位换算：1 meV/Å² ≈ 0.016 J/m²。

---

## 3. IVDW：原子对加和与多体色散修正

### 3.1 方法总览与版本要求

| `IVDW` | 方法 | 可用版本 / 编译要求 |
|---|---|---|
| 1 或 10 | DFT-D2（Grimme） | VASP ≥ 5.2.11 |
| 11 | DFT-D3，零阻尼 D3(0) | VASP ≥ 5.3.4 |
| 12 | DFT-D3，Becke–Johnson 阻尼 D3(BJ) | VASP ≥ 5.3.4 |
| 13 | DFT-D4 | VASP ≥ 6.2，外部库 dftd4，编译 `-DDFTD4`；dftd4 ≥ 4.0 的新 API 需 VASP ≥ 6.6.0（旧 API ≤ 3.7 用 `-DDFTD4_API_V3`） |
| 15 | simple-DFT-D3 库中的 D3 变体 | VASP ≥ 6.6.0，外部库 simple-dftd3，编译 `-DSDFTD3` |
| 3 | DFT-ulg | VASP ≥ 5.3.5 |
| 4 | dDsC | VASP ≥ 5.4.1 |
| 2 或 20 | Tkatchenko–Scheffler（TS） | VASP ≥ 5.3.3 |
| 21 | TS + 迭代 Hirshfeld 分区（TS-HI） | VASP ≥ 5.3.5 |
| 202 | MBD@rsSCS | VASP ≥ 5.4.1 |
| 263 | MBD@rSC/FI（分数离子极化率模型） | VASP ≥ 6.1.0 |
| 14 | libMBD 库中的 TS/MBD 方法（`LIBMBD_METHOD`） | VASP ≥ 6.4.3，外部库 libMBD，编译 `-DLIBMBD` |

旧标签 `LVDW = .TRUE.` 已废弃，等价于 `IVDW = 1`。外部库的编译方法（在 `makefile.include` 的 `CPP_OPTIONS` 中加对应宏、在 `LLIBS`/`INCS` 中加库和头文件路径）见 [第1章](01-setup-compile.md) 与官方 [Makefile.include](https://www.vasp.at/wiki/index.php/Makefile.include) 页面。

### 3.2 DFT-D3（IVDW = 11 / 12）

对 PBE、RPBE、revPBE、PBEsol、TPSS、M06L、SCAN、B3LYP、PBE0 和 HF，VASP 内置了 D3 阻尼参数，**只写 `IVDW` 即可**。其他泛函（例如 r2SCAN、optB86b）没有内置参数时，要么手动给出参数，要么改用 D4 / simple-DFT-D3（它们的参数表更全）。可调参数：

| 标签 | 含义 | 适用 |
|---|---|---|
| `VDW_S6` | 偶极–偶极项缩放 s6 | D3(0)/D3(BJ)：VASP ≥ 6.6.0 才可设置 |
| `VDW_S8` | 偶极–四极项缩放 s8 | D3(0)/D3(BJ) |
| `VDW_SR` / `VDW_SR8` | 零阻尼的半径缩放（`VDW_SR8` 需 ≥ 6.6.0） | D3(0) |
| `VDW_A1` / `VDW_A2` | BJ 阻尼参数 a1、a2 | D3(BJ) |
| `VDW_RADIUS` | 两体作用截断（Å） | 全部；6.6.0 修复了大于默认截断时不生效的问题 |
| `VDW_CNRADIUS` | 配位数截断（Å） | 全部 |

> [!TIP]
> 除非你很清楚自己在做什么（例如复现某篇文章的自定义参数），不要改动默认阻尼参数。D3(BJ) 通常比 D3(0) 更好地描述短程区域，是大多数固体/表面计算的默认推荐。

### 3.3 DFT-D4（IVDW = 13）与 simple-DFT-D3（IVDW = 15）

- **D4**：对 PBE、SCAN、HSE06 等常见泛函，`IVDW = 13` 即可，VASP 会自动把泛函名传给 dftd4；其他泛函用 `DFTD4_XC` 指定名称（6.6.0+）。6.6.0 还新增 `DFTD4_MODEL`、`VDW_S9`（三体项）、`VDW_RADIUS`、`VDW_CNRADIUS` 的设置。
- **simple-DFT-D3**（6.6.0+）：`IVDW = 15` 并用 `SDFTD3_DAMPING` 选阻尼函数：`zero`、`rational`（即 BJ）、`mzero`、`mrational`、`optimizedpower`；`SDFTD3_XC` 指定泛函名；`VDW_S9` 打开 Axilrod–Teller–Muto 三体项。与 VASP 内置 D3 相比，它的阻尼选项更多、元素覆盖到 Lr、速度更快。

```text
# PBE-D4 (VASP >= 6.2 compiled with dftd4)
IVDW = 13

# PBE-D3(BJ) through the simple-DFT-D3 library (VASP >= 6.6.0, -DSDFTD3)
IVDW           = 15
SDFTD3_DAMPING = rational
```

> [!NOTE]
> dftd4 与 simple-dftd3 库的新版本可能默认采用平滑截断，周期体系的结果会与旧版本、以及 VASP 内置 D3 略有差别。项目内比较能量时要固定库的版本。

### 3.4 TS、MBD 与 libMBD

TS 类方法（`IVDW = 20/21`）和 MBD（`IVDW = 202/263`）的 C6 来自电子密度的 Hirshfeld 分区，适合离子性较强或极化率随环境变化大的体系。libMBD 实现（`IVDW = 14`）用解析力，比 VASP 自带实现（数值力）快得多，官方推荐有条件时使用；但要注意：

- VASP ≤ 6.5.1 的 libMBD 接口在结构优化中**不会把新的原子坐标和晶胞传给 libMBD**（known issue 55，6.6.0 修复，6.5.1 有官方补丁）；
- 传给 libMBD 的 k 网格是按默认 `KSPACING = 0.5` 生成的，而不是 KPOINTS 中的网格（known issue 91，截至本文仍未修复）；
- libMBD 编译时不要开 ScaLAPACK/MPI，否则 NEB 计算结果错误。

### 3.5 检查色散修正是否生效

```bash
grep -i "IVDW" OUTCAR | head -3                  # the tag was read
grep -i -E "edisp|dispersion|vdw" OUTCAR | head  # dispersion energy / parameter lines
```

不同方法打印的内容不同（D3 会列出所用参数，D4 自 6.4.3 起输出更多信息），至少应确认 `IVDW` 被正确读入、没有 "not compiled" 一类的报错。

---

## 4. 非局域 vdW-DF 泛函（LUSE_VDW）

### 4.1 常用泛函的 INCAR 写法

所有 vdW-DF 都用 `LUSE_VDW = .TRUE.` 打开非局域关联，并建议 `LASPH = .TRUE.`。下表的设置取自 VASP Wiki（`ZAB_VDW` 默认 −0.8491）：

| 泛函 | INCAR 设置 |
|---|---|
| vdW-DF（Dion） | `GGA = RE`、`AGGAC = 0.0` |
| vdW-DF2 | `GGA = ML`、`AGGAC = 0.0`、`ZAB_VDW = -1.8867` |
| optPBE-vdW | `GGA = OR`、`AGGAC = 0.0` |
| optB88-vdW | `GGA = BO`、`PARAM1 = 0.1833333333`、`PARAM2 = 0.22`、`AGGAC = 0.0` |
| **optB86b-vdW** | `GGA = MK`、`PARAM1 = 0.1234`、`PARAM2 = 1.0`、`AGGAC = 0.0` |
| rev-vdW-DF2（vdW-DF2-B86R） | `GGA = MK`、`PARAM1 = 0.1234568`、`PARAM2 = 0.7114`、`AGGAC = 0.0`、`ZAB_VDW = -1.8867` |
| vdW-DF-cx | `GGA = CX`、`AGGAC = 0.0` |
| rVV10 | `GGA = ML`、`IVDW_NL = 2`（`BPARAM = 6.3`、`CPARAM = 0.0093` 为默认） |
| SCAN+rVV10 | `METAGGA = SCAN`、`BPARAM = 15.7` |
| r2SCAN+rVV10 | `METAGGA = R2SCAN`、`BPARAM = 11.95` |
| PBE+rVV10L | `GGA = PE`、`BPARAM = 10` |
| BEEF-vdW | `GGA = BF`、`ZAB_VDW = -1.8867`；需 libbeef（`-Dlibbeef`）或 Libxc 编译 |

optB86b-vdW 的完整写法：

```text
GGA      = MK
PARAM1   = 0.1234
PARAM2   = 1.0
AGGAC    = 0.0           # remove PBE gradient correction to correlation
LUSE_VDW = .TRUE.        # non-local correlation (Dion kernel)
LASPH    = .TRUE.        # recommended for all vdW-DF functionals
```

r2SCAN+rVV10 的写法：

```text
METAGGA  = R2SCAN
LUSE_VDW = .TRUE.
BPARAM   = 11.95
CPARAM   = 0.0093
LASPH    = .TRUE.
LMIXTAU  = .TRUE.        # usually helps SCF convergence of meta-GGAs (see chapter 6)
```

### 4.2 核函数文件 vdw_kernel.bindat 与版本差异

| 版本 | 行为 |
|---|---|
| ≤ 6.4.2 | Dion 核（`IVDW_NL = 1`，即除 rVV10 外的所有 vdW-DF）需要预先计算的 `vdw_kernel.bindat`（随 VASP 在下载门户分发）。找不到时 VASP 会现场生成，**耗时数小时**。请把它复制到计算目录 |
| ≥ 6.4.3 | 核函数生成被大幅加速并行化（8 个 MPI 进程约 2 分钟），不复制也可以；不同核（`IVDW_NL = 1/3/4`）的文件不兼容时会自动重新生成 |
| rVV10 | 任何版本都不需要核文件 |
| ≥ 6.4.0 | 新增 `IVDW_NL`，可任意组合 GGA/meta-GGA 与 Dion/rVV10 核；新增 `LSPIN_VDW`（自旋极化的非局域项，仅 Dion 核 + GGA） |
| ≤ 6.2 | rVV10 的应力张量有误，变胞优化请用 ≥ 6.3.0 |

### 4.3 使用注意

- **POTCAR**：没有专门的 vdW-DF POTCAR，直接用 PBE 的 PAW 势；非局域项用赝价电子密度 + 部分芯密度近似计算，用更硬或更多价电子的 POTCAR 时精度更高。
- **数值噪声**：部分 vdW-DF 的能量噪声较大，`ALGO = All`（共轭梯度）尤其容易出现能量突然上升而不收敛；可以加密 FFT 网格（`PREC = Accurate`）。
- **不要叠加**：`LUSE_VDW` 已经包含色散，不要再同时设置 `IVDW`。
- **GPU**：用 Cray 编译器在 AMD GPU（OpenMP offload）上运行 `LUSE_VDW` 时，应力张量可能因竞争条件出错（known issue），变胞优化前请与 `OMP_TARGET_OFFLOAD=DISABLED` 的结果对比。

---

## 5. 如何选择 vdW 方法，以及它们的限制

| 场景 | 推荐 | 说明 |
|---|---|---|
| 层状材料层间距、堆垛 | PBE-D3(BJ)、optB86b-vdW、rev-vdW-DF2、r2SCAN+rVV10 | 层间距通常在实验值的百分之几以内；同时关注面内晶格常数 |
| 分子晶体、有机/无机界面 | D3(BJ)、D4、MBD | 多体效应强时考虑 MBD |
| 金属表面上的分子吸附 | D3(BJ)、optB86b-vdW、BEEF-vdW | 不同方法吸附能可差 0.1–0.3 eV；BEEF 还能给误差估计 |
| meta-GGA（r2SCAN） | r2SCAN+rVV10、r2SCAN-D4 | 内置 D3 没有 r2SCAN 参数 |
| 杂化泛函（HSE06） | HSE06-D4（`IVDW = 13`）或 simple-DFT-D3 | VASP 内置 D3 官方列出的杂化泛函参数只有 B3LYP 与 PBE0；用内置 D3 搭配 HSE06 时请核对/手动设置参数 |
| 有限位移声子、MD | 任意 IVDW 或 vdW-DF | 受力、应力都可用 |
| DFPT 声子（`IBRION = 7/8`） | **不支持 IVDW 修正** | 需要 vdW 时改用有限位移（`IBRION = 5/6` 或 phonopy），见 [第12章](12-phonon.md) |

> [!IMPORTANT]
> 一个项目中参与比较的所有能量（体相、单层、分子参考态、吸附体系……）必须使用同一种 vdW 方法、同一套参数和同一版本的外部库。论文中请写明方法、阻尼形式和 VASP/库的版本。

---

## 6. 实战一：石墨层间距（PBE vs D3(BJ) vs optB86b-vdW）

完整输入见 [../examples/08-graphite-vdw/README.md](../examples/08-graphite-vdw/README.md)。

### 6.1 输入文件

POSCAR（Bernal AB 堆垛，P6₃/mmc，实验初始结构）：

```text
graphite 2H (Bernal AB), experimental start a = 2.464 A, c = 6.711 A
1.0
   2.4640000000   0.0000000000   0.0000000000
  -1.2320000000   2.1338865949   0.0000000000
   0.0000000000   0.0000000000   6.7110000000
C
4
Direct
   0.0000000000   0.0000000000   0.2500000000
   0.0000000000   0.0000000000   0.7500000000
   0.3333333333   0.6666666667   0.2500000000
   0.6666666667   0.3333333333   0.7500000000
```

KPOINTS（六角晶胞必须用 Γ 中心网格；用偶数 Monkhorst–Pack 网格会触发 IBZKPT 报错）：

```text
Gamma-centred mesh (hexagonal cell -> always Gamma-centred)
0
Gamma
15 15 5
0 0 0
```

INCAR_D3BJ：

```text
SYSTEM  = graphite, PBE-D3(BJ)

PREC    = Accurate
ENCUT   = 520            # 1.3 x ENMAX(C = 400 eV): needed for ISIF = 3 (Pulay stress)
EDIFF   = 1E-6
LREAL   = .FALSE.
ISMEAR  = 0              # graphite is a semimetal; small Gaussian smearing is safe
SIGMA   = 0.05

# ---- dispersion correction ----
IVDW    = 12             # DFT-D3 with Becke-Johnson damping (parameters set for PBE automatically)

IBRION  = 2
ISIF    = 3              # relax a and c
NSW     = 80
EDIFFG  = -0.005

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

INCAR_optB86b 只把 `IVDW = 12` 换成 §4.1 的 optB86b-vdW 六行；INCAR_PBE 去掉 `IVDW`（对照组）。

| 关键词 | 值 | 为什么 |
|---|---|---|
| `ENCUT` | 520 | 变胞优化需要 1.3×ENMAX 以上，否则 Pulay 应力会使 c 偏小 |
| `ISIF` | 3 | 同时优化 a 和 c；层状材料的 c 对应力非常敏感 |
| `EDIFFG` | −0.005 | 层间作用弱，势能面平坦，受力判据要比常规更严 |
| `ISMEAR` | 0 | 石墨是半金属，小 SIGMA 的 Gaussian 展宽稳妥 |

### 6.2 运行与检查

```bash
cd ../examples/08-graphite-vdw
cat $VASP_PP_PATH/potpaw_PBE/C/POTCAR > POTCAR
bash run_compare.sh relax
# ... after the three jobs finish:
grep "reached required accuracy" */relax/OUTCAR
for m in PBE D3BJ optB86b; do echo $m; sed -n 3,5p $m/relax/CONTCAR; done
```

预期趋势（定性）：

- **PBE**：c 迅速增大，远超实验值 6.711 Å，可能在 `NSW` 用完前都不收敛——层间几乎没有结合；
- **PBE-D3(BJ)、optB86b-vdW**：c 回到 6.6–6.8 Å 附近，a ≈ 2.46–2.47 Å。

> [!TIP]
> 对层状材料，变胞优化的结果建议再用"固定 a、扫描 c"的 E–c 曲线核对：取 c = 6.2–7.4 Å 的 7–9 个点，每点 `ISIF = 2`（或固定原子的单点），对 E(c) 拟合多项式求极小。这能排除优化器在平坦势能面上"提前停下"的问题。

生成一组缩放 c 的 POSCAR：

```python
import numpy as np
from ase.io import read, write

atoms = read("D3BJ/relax/CONTCAR")
for c in np.arange(6.2, 7.41, 0.15):
    a = atoms.copy()
    cell = a.cell.copy()
    cell[2] = cell[2] / np.linalg.norm(cell[2]) * c
    a.set_cell(cell, scale_atoms=True)        # keep fractional coordinates
    write(f"POSCAR_c{c:.2f}", a, format="vasp", direct=True)
```

---

## 7. 实战二：石墨的层间结合能

在各自的优化结构上做单点能：石墨（4 原子）和相同面内晶格、c = 20 Å 的单层石墨烯（2 原子）。

```bash
bash run_compare.sh static     # graphite + graphene for each method
bash run_compare.sh collect    # -> analyze_graphite.py
```

`make_graphene.py` 从石墨的 CONTCAR 读出 a1、a2，生成石墨烯 POSCAR；石墨烯的 KPOINTS 自动设为 15×15×1。`analyze_graphite.py` 计算：

$$
E_b=\frac{E_{\rm graphene}}{2}-\frac{E_{\rm graphite}}{4}
$$

示意输出（数值仅供参考，不是实际计算结果）：

```text
method        a (A)    c (A)  E_b (meV/atom)  E_b (meV/A^2)
PBE          2.4670   8.9xxx             1.x           0.xx
D3BJ         2.4640   6.7xxx            5x.x          2x.xx
optB86b      2.4650   6.6xxx            6x.x          2x.xx
```

参考值：RPA 给出石墨层间结合能约 48 meV/atom（Lebègue et al., PRL 105, 196401 (2010)）；Björkman 等对大量层状化合物的系统研究表明，层间结合能普遍约 20 meV/Å²（PRL 108, 235502 (2012)）。不同 vdW 方法的偏差可达 10–30%，这本身就是方法选择需要考虑的因素。

> [!NOTE]
> 石墨烯与石墨必须使用相同的面内 k 网格、ENCUT 和 smearing，石墨烯的 k 点在 z 方向只取 1。结合能是两个大数之差，k 点不足会带来几 meV/atom 的误差。

---

## 8. 二维材料建模要点：以单层 MoS2 为例

完整输入见 [../examples/08-mos2-monolayer/README.md](../examples/08-mos2-monolayer/README.md)；只做面内晶格优化的基础版见 [第4章](04-relax-eos.md) §6。

### 8.1 必须注意的设置

| 项目 | 建议 | 原因 |
|---|---|---|
| 真空层 | ≥ 15–20 Å，并做收敛测试 | 防止层与周期镜像相互作用；功函数、带边对齐对真空更敏感 |
| 原子位置 | 放在 c 方向中央（z ≈ 0.5） | 避免原子跨越周期边界；偶极校正时 DIPOL 取质心 |
| k 点 | Γ 中心 n×n×1 | 真空方向没有色散；六角晶格必须 Γ 中心 |
| 晶格优化 | `ISIF = 3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` | 只放开面内晶格；VASP ≥ 6.4.3（推荐 ≥ 6.5.1） |
| vdW | 与双层/异质结/体相保持一致 | 单层面内晶格受 vdW 影响小，但能量比较需一致 |
| 偶极 | 对称单层不需要；Janus（如 MoSSe）、吸附一侧时需要 `LDIPOL` | 见 [第9章](09-surface-workfunction.md) |
| 应力/弹性 | VASP 的应力按含真空的晶胞体积归一化 | 面内量要乘以 c 换算成 N/m，见 [第13章](13-elastic.md) |
| k 路径 | Γ–M–K–Γ | 六角二维布里渊区 |

### 8.2 优化输入

```text
SYSTEM  = MoS2 monolayer, PBE-D3(BJ), in-plane relaxation only

PREC    = Accurate
ENCUT   = 450            # >= 1.3 x max(ENMAX) of Mo_pv / S; ISIF = 3 needs the margin
EDIFF   = 1E-6
LREAL   = .FALSE.
LASPH   = .TRUE.
ISMEAR  = 0
SIGMA   = 0.05

IVDW    = 12             # same vdW setting as the bilayer / heterostructure you compare to

IBRION  = 2
ISIF    = 3
# relax a1, a2 only and keep c (vacuum) fixed; VASP >= 6.4.3 (6.5.1+ recommended)
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.
NSW     = 100
EDIFFG  = -0.01

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

POSCAR 中的 MoS2 层位于 z = c/2，c = 23.13 Å（约 20 Å 真空），KPOINTS 为 Γ 中心 12×12×1。

```bash
cd ../examples/08-mos2-monolayer
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
mkdir -p relax && cp INCAR_relax relax/INCAR && cp POSCAR KPOINTS POTCAR job.slurm relax/
(cd relax && sbatch job.slurm)
head -5 relax/CONTCAR                 # a changed; the third lattice vector must stay 0 0 23.13
```

> [!WARNING]
> VASP ≤ 6.4.2 中 `LATTICE_CONSTRAINTS` 与 `IBRION = 1/2` 一起使用时会被静默忽略——c 轴照样会变，真空层可能被压缩。没有新版本时，改用"固定 c、扫描 a"的 E–a 拟合，或社区补丁 [VASP_OPT_AXIS](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS)（仅链接）。

### 8.3 真空收敛测试

`vacuum_test.sh` 用 `set_vacuum.py` 保持原子的笛卡尔坐标不变、把 c 改为 15/20/25/30 Å 并重新居中：

```bash
bash vacuum_test.sh prepare
bash vacuum_test.sh collect
```

判据：相邻两个 c 的总能量差 < 1 meV/atom。注意每个晶胞的 E_F 绝对值都不同——VASP 的势能零点是晶胞平均静电势，随真空大小改变，**只有相对真空能级的能量才有物理意义**。

### 8.4 带边对齐到真空能级

静态计算（`INCAR_scf`）中设置 `LVHAR = .TRUE.`，LOCPOT 只包含离子势 + Hartree 势（不含交换关联势），真空区更快变平：

```bash
python ../09-pt111-workfunction-dipole/plot_locpot_planar.py c_20/LOCPOT --outcar c_20/OUTCAR \
       --vbm <VBM_eV> --cbm <CBM_eV>
```

得到真空能级 E_vac 后，电离能 IP = E_vac − E_VBM，电子亲和能 EA = E_vac − E_CBM。VBM/CBM 可取自 OUTCAR 的带隙报告（VASP ≥ 6.4.3）或能带计算。VASPKIT 也提供二维材料工具（921 居中、922 调真空、923 标准化二维晶胞、927 以真空能级为参考的带边），功能编号以 [vaspkit.com](https://vaspkit.com/) 的 1.5.x 文档为准。

参考量级（近似值）：单层 MoS2 的 PBE 面内晶格常数约 3.18–3.19 Å（实验约 3.16 Å）；PBE 带隙为 K 点直接带隙约 1.6–1.7 eV（不含 SOC），SOC 使 K 点价带顶劈裂约 0.15 eV（SOC 计算见 [第7章](07-magnetism-u-soc.md) §8）。

### 8.5 二维库仑截断（VASP ≥ 6.5.0，选读）

VASP 6.5.0 引入了库仑核截断（Coulomb kernel truncation），可去除二维体系沿法向与周期镜像之间的静电相互作用，对带电或有偶极的二维体系、表面尤其有用：

```text
KERNEL_TRUNCATION {
     LTRUNCATE       = T
     IDIMENSIONALITY = 2
     ISURFACE        = 3
     LCOARSEN        = T
}
```

`LTRUNCATE` 是"总开关"，不打开时其他子标签被忽略。使用填充真空（padding）时，原子不能位于法向的晶胞边界，应把结构居中。6.6.0 的二维截断会产生噪声势（known issue 72），6.6.1 已修复；6.5.x 中还有法向 Ewald 力的符号错误等问题（6.6.0 修复）。常规中性、对称的二维材料用足够的真空即可，不必启用。

---

## 9. 实战三：双层 MoS2 层间结合能曲线

完整输入见 [../examples/08-mos2-bilayer-binding/README.md](../examples/08-mos2-bilayer-binding/README.md)。

### 9.1 结构

体相 2H-MoS2 中相邻层是 AA' 堆垛：上层相对下层旋转 180°，上层 Mo 正好位于下层 S 的上方。`make_bilayer.py` 用两个参数描述单层（a 和 S–S 高度 t），层间距 d 定义为上下两层 Mo 原子的竖直距离（体相实验值 c/2 ≈ 6.15 Å）：

```bash
python make_bilayer.py --d 6.2 --a 3.19 --t 3.13 --c 30 > POSCAR     # bilayer
python make_bilayer.py --mono   --a 3.19 --t 3.13 --c 30 > POSCAR     # monolayer reference
```

生成的双层 POSCAR（d = 6.15 Å 时）：

```text
MoS2 bilayer 2H (AA') d(Mo-Mo)=6.150 a=3.19 t=3.13 c=30.0
1.0
    3.1900000000   0.0000000000   0.0000000000
   -1.5950000000   2.7626210381   0.0000000000
    0.0000000000   0.0000000000  30.0000000000
Mo S
2 4
Direct
    0.3333333333   0.6666666667   0.3975000000
    0.6666666667   0.3333333333   0.6025000000
    0.6666666667   0.3333333333   0.3453333333
    0.6666666667   0.3333333333   0.4496666667
    0.3333333333   0.6666666667   0.5503333333
    0.3333333333   0.6666666667   0.6546666667
```

（Mo–S 键长 2.42 Å，层间 S–S 竖直间隙约 3.02 Å。）

### 9.2 扫描与分析

INCAR 为 PBE-D3(BJ) 单点（`IVDW = 12`、`ENCUT = 450`、`NSW = 0`），KPOINTS 为 Γ 中心 12×12×1，单层参考在**同一晶胞**中计算：

```bash
cd ../examples/08-mos2-bilayer-binding
cat $VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR
export A=3.19 T=3.13 C=30.0          # take a and t from your relaxed monolayer
bash run.sh prepare                  # mono/ + d_5.6 ... d_10.0
bash run.sh collect                  # E_b(d) table, cubic fit, binding_curve.png
```

$$
E_b(d)=\frac{E_{\rm bilayer}(d)-2E_{\rm mono}}{A},\qquad A=\frac{\sqrt3}{2}a^2
$$

示意输出（数值仅供参考）：

```text
 d(Mo-Mo) A   E_b meV/A^2  E_b meV/f.u.
       6.00        -1x.xx        -1xx.x
       6.20        -2x.xx        -1xx.x
      10.00         -x.xx         -xx.x
cubic fit around the minimum: d_eq ~ 6.2 A, E_b ~ -2x meV/A^2
```

检查：

- 极小值应在 d ≈ 6.1–6.3 Å 附近，大 d 处 E_b → 0；
- 去掉 `IVDW` 的对照计算几乎没有结合；
- 换成 `IVDW = 11`、`13` 或 optB86b-vdW 比较 d_eq 与 E_b，体会方法依赖性。

> [!TIP]
> 刚性扫描只是第一步。在 d_eq 附近做一次完整弛豫，或用 `ISIF = 3` + `LATTICE_CONSTRAINTS` 优化双层后，再与单层比较。研究不同堆垛（AA、AB、AA'……）时，只需修改上层原子的面内位置；VASPKIT 也提供堆垛能面相关功能。

---

## 10. 双层与异质结的晶格失配和超胞匹配（简介）

两种二维材料（或体相与二维材料）堆叠时，面内晶格常数一般不相等。常见处理流程：

1. **计算失配度**：$\delta = (L_{\rm film}-L_{\rm sub})/L_{\rm sub}$，其中 L 是两边超胞的边长。一般希望 |δ| ≲ 1–3%，并在论文中说明应变加在哪一层（或平均分摊）。
2. **寻找公度超胞**：对两套晶格搜索整数倍（必要时加旋转）使边长接近。下面的小脚本对两种六角晶格搜索 n×n 与 m×m（R0°）或 m(√3×√3)R30° 的组合：

```python
import math

a_sub, a_film = 3.19, 2.467        # MoS2 and graphene in-plane lattice constants (A)
nat_sub, nat_film = 3, 2           # atoms per primitive cell
cands = []
for n in range(1, 9):              # substrate n x n
    for m in range(1, 13):         # film m x m (R0) or m(sqrt3 x sqrt3)R30
        for rot, f in (("R0", 1.0), ("R30", math.sqrt(3.0))):
            L_sub, L_film = n * a_sub, m * f * a_film
            mis = (L_film - L_sub) / L_sub
            if abs(mis) < 0.03:
                nat = n * n * nat_sub + round((m * f) ** 2) * nat_film
                cands.append((nat, abs(mis), n, m, rot, L_sub))
for nat, mis, n, m, rot, L in sorted(cands)[:6]:
    print(f"MoS2 {n}x{n} / graphene {m}x{m} {rot}: L = {L:6.2f} A, "
          f"strain on graphene = {100 * mis:4.2f} %, atoms = {nat}")
```

运行结果（这是脚本的真实输出，只依赖输入的两个晶格常数）：

```text
MoS2 4x4 / graphene 3x3 R30: L =  12.76 A, strain on graphene = 0.46 %, atoms = 102
MoS2 7x7 / graphene 9x9 R0: L =  22.33 A, strain on graphene = 0.57 %, atoms = 309
MoS2 8x8 / graphene 6x6 R30: L =  25.52 A, strain on graphene = 0.46 %, atoms = 408
```

即 MoS2 的 4×4 超胞与石墨烯的 3(√3×√3)R30° 超胞只差约 0.5%，102 个原子，是常用的计算模型。

3. **更通用的工具**：pymatgen 的 `ZSLGenerator`（`pymatgen.analysis.interfaces.zsl`，Zur–McGill 晶格匹配算法）与 `CoherentInterfaceBuilder`（`pymatgen.analysis.interfaces.coherent_interfaces`）可以对任意晶面自动搜索匹配并生成界面结构；VASPKIT 804（由两个 slab 建异质结）、825（六角莫尔超晶格）也能完成类似工作。建模的通用知识见 Common 第3章。
4. **计算与分析**：异质结同样要固定真空、居中、使用一致的 vdW；不对称异质结必须加偶极校正。能带排列可用两种方式：分别计算两个单层、以各自的真空能级为参考（Anderson 规则）；或在异质结中用投影能带/PDOS 直接读带边，两者的差别反映了界面偶极与杂化。

---

## 参数选择与收敛建议

| 项目 | 建议值 / 做法 | 备注 |
|---|---|---|
| vdW 方法 | 固体/层状：D3(BJ) 或 optB86b-vdW；有 D4 库时可用 D4 | 一个项目内统一 |
| ENCUT | 变胞优化 ≥ 1.3×ENMAX | 层状材料 c 轴对 Pulay 应力敏感 |
| EDIFFG | −0.005 ~ −0.01 eV/Å | 层间势能面平坦 |
| 层间距 | 变胞优化 + E–c 扫描核对 | 防止优化器提前停止 |
| 二维真空 | ≥ 15–20 Å，测试到 ΔE < 1 meV/atom | 带边/功函数需要更大真空 |
| 二维 k 点 | Γ 中心 n×n×1 | 单层与多层面内网格一致 |
| 结合能 | 同晶胞（同 a、c）、同 k、同 smearing 计算参考态 | 两个大数相减，误差要小 |
| vdW-DF 核 | ≤ 6.4.2 复制 `vdw_kernel.bindat`；≥ 6.4.3 可现场生成 | rVV10 不需要 |
| 外部库 | 记录 dftd4 / simple-dftd3 / libMBD 版本 | 平滑截断等默认行为可能随版本变 |
| libMBD | ≥ 6.6.0（或打补丁的 6.5.1） | known issue 55、91 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| 纯 PBE 优化石墨/MoS2 体相时 c 一直变大、优化不收敛 | GGA 没有层间色散 | 加 `IVDW` 或用 vdW-DF |
| vdW-DF 计算开始后长时间卡在第一步 | ≤ 6.4.2 找不到 `vdw_kernel.bindat`，正在现场生成核函数（可能数小时） | 从 VASP 门户获取并复制到计算目录；或升级到 ≥ 6.4.3 |
| vdW-DF 计算 `ALGO = All` 时能量突然升高、SCF 不收敛 | 非局域项数值噪声 | `PREC = Accurate`、加密 FFT；改用 `ALGO = Normal/Fast` |
| `IVDW = 13` 或 `15` 时程序报错退出 / 提示未编译 | VASP 编译时没有链接 dftd4 或 simple-dftd3 | 按 [第1章](01-setup-compile.md) 重新编译；或改用内置 `IVDW = 12` |
| 换了 dftd4 版本后能量变了几 meV | 库版本的默认截断（如平滑截断）改变 | 固定库版本；重新计算全部参考态 |
| D3 用于 r2SCAN、optB86b 等泛函时结果异常 | 内置 D3 没有该泛函的参数 | 用 D4/simple-DFT-D3（自动识别泛函），或 vdW-DF 类 |
| DFPT（`IBRION = 7/8`）声子未包含 vdW | DFPT 不支持 IVDW 修正 | 改用 `IBRION = 5/6` 或 phonopy 有限位移 |
| 二维材料变胞优化后真空层塌缩/变化 | 没有约束 c 轴；或 ≤ 6.4.2 中 `LATTICE_CONSTRAINTS` 被忽略 | 用 ≥ 6.4.3（推荐 6.5.1+）；或 E–a 扫描 |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角/二维晶胞用了 Monkhorst–Pack 偶数网格 | 改为 Gamma 中心网格 |
| `Tetrahedron method fails (number of k-points < 4)` | 二维/分子体系用了 `ISMEAR = -5` 且 k 点太少 | 用 `ISMEAR = 0`、`SIGMA = 0.05` |
| libMBD 结构优化中色散能不随结构变化 | known issue 55（≤ 6.5.1） | 升级到 6.6.0+ 或打官方补丁 |
| 不同 c 的二维计算 E_F 差很多 | 势能零点随晶胞变化，属正常 | 以真空能级为参考比较 |

---

## 练习

1. 对石墨用 D3(0)（`IVDW = 11`）、D3(BJ)、optB86b-vdW、rVV10 优化 c，画出"c 与结合能"的方法对比图。
2. 用"固定 a、扫描 c"的方法得到 PBE-D3(BJ) 石墨的 E–c 曲线，与 `ISIF = 3` 的结果比较。
3. 对单层 MoS2 做 c = 15–35 Å 的真空收敛，分别以总能、E_F 和"VBM 相对真空能级"为指标，比较哪个收敛最慢。
4. 修改 `make_bilayer.py`，比较 2H（AA'）与 AA 堆垛双层 MoS2 的结合能和平衡层间距。
5. 用 §10 的脚本为 MoS2/WS2、graphene/h-BN 寻找失配 < 1% 的超胞，估算原子数。

---

## 参考资料

**VASP 官方文档**

- IVDW 总览：<https://www.vasp.at/wiki/index.php/IVDW>
- DFT-D3：<https://www.vasp.at/wiki/index.php/DFT-D3>；DFT-D4：<https://www.vasp.at/wiki/index.php/DFT-D4>；simple-DFT-D3：<https://www.vasp.at/wiki/index.php/Simple-DFT-D3>
- libMBD：<https://www.vasp.at/wiki/index.php/LIBMBD_METHOD>
- 非局域 vdW-DF：<https://www.vasp.at/wiki/index.php/Nonlocal_vdW-DF_functionals>、<https://www.vasp.at/wiki/index.php/LUSE_VDW>、<https://www.vasp.at/wiki/index.php/IVDW_NL>
- 编译外部库：<https://www.vasp.at/wiki/index.php/Makefile.include>
- LATTICE_CONSTRAINTS：<https://www.vasp.at/wiki/index.php/LATTICE_CONSTRAINTS>
- 库仑核截断：<https://www.vasp.at/wiki/index.php/KERNEL_TRUNCATION/LTRUNCATE>
- Known issues：<https://www.vasp.at/wiki/index.php/Known_issues>；Changelog：<https://www.vasp.at/wiki/index.php/Changelog>

**相关 GitHub 仓库（仅链接）**

- [dftd4/dftd4](https://github.com/dftd4/dftd4)：IVDW=13 的外部库；VASP 6.6.0 适配 4.x 新 API
- [dftd3/simple-dftd3](https://github.com/dftd3/simple-dftd3)：IVDW=15（VASP 6.6.0+）的外部库
- [Chengcheng-Xiao/VASP_OPT_AXIS](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS)：旧版 VASP 固定晶格分量的补丁（6.4.3+ 用 LATTICE_CONSTRAINTS 代替）
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：`ZSLGenerator`、`CoherentInterfaceBuilder` 界面匹配
- [skelton-group/VASP-Workshop-Chester-2019](https://github.com/skelton-group/VASP-Workshop-Chester-2019)：含分子晶体 vdW 案例的 CC-BY 工作坊材料
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book)：VASP + ASE 经典教程
- [vaspkit/VASPKIT_NatureProtocols](https://github.com/vaspkit/VASPKIT_NatureProtocols)：VASPKIT 二维材料工具等功能的示例文件

**关键文献**

- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys. 132, 154104 (2010)（DFT-D3）. DOI: [10.1063/1.3382344](https://doi.org/10.1063/1.3382344)
- S. Grimme, S. Ehrlich, L. Goerigk, J. Comput. Chem. 32, 1456 (2011)（BJ 阻尼）. DOI: [10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759)
- E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)（DFT-D4）. DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
- A. Tkatchenko, M. Scheffler, PRL 102, 073005 (2009)（TS）. DOI: [10.1103/PhysRevLett.102.073005](https://doi.org/10.1103/PhysRevLett.102.073005)
- M. Dion et al., PRL 92, 246401 (2004)（vdW-DF）. DOI: [10.1103/PhysRevLett.92.246401](https://doi.org/10.1103/PhysRevLett.92.246401)
- J. Klimeš, D. R. Bowler, A. Michaelides, PRB 83, 195131 (2011)（optB86b-vdW 等在固体中的应用）. DOI: [10.1103/PhysRevB.83.195131](https://doi.org/10.1103/PhysRevB.83.195131)
- R. Sabatini, T. Gorni, S. de Gironcoli, PRB 87, 041108 (2013)（rVV10）. DOI: [10.1103/PhysRevB.87.041108](https://doi.org/10.1103/PhysRevB.87.041108)
- H. Peng et al., PRX 6, 041005 (2016)（SCAN+rVV10）. DOI: [10.1103/PhysRevX.6.041005](https://doi.org/10.1103/PhysRevX.6.041005)
- J. Wellendorff et al., PRB 85, 235149 (2012)（BEEF-vdW）. DOI: [10.1103/PhysRevB.85.235149](https://doi.org/10.1103/PhysRevB.85.235149)
- S. Lebègue et al., PRL 105, 196401 (2010)（石墨 RPA 结合能）. DOI: [10.1103/PhysRevLett.105.196401](https://doi.org/10.1103/PhysRevLett.105.196401)
- T. Björkman et al., PRL 108, 235502 (2012)（层状化合物的 vdW 结合）. DOI: [10.1103/PhysRevLett.108.235502](https://doi.org/10.1103/PhysRevLett.108.235502)

[← 上一章](07-magnetism-u-soc.md) | [返回目录](../README.md) | [下一章 →](09-surface-workfunction.md)
