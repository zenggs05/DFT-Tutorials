# 第10章 吸附、催化与自由能（含隐式溶剂）（Adsorption, Catalysis, Free Energy and Implicit Solvation with VASP）

本章以 CO/Pt(111) 和 H/Pt(111) 为例，手把手完成 VASP 表面吸附计算的完整流程：吸附能（adsorption energy）、吸附位点与覆盖度、气相分子参考、频率计算与 ZPE/熵校正、计算氢电极（CHE）模型下的 HER/OER/ORR 自由能台阶图与过电位、d 带中心（d-band center），以及 VASPsol 隐式溶剂（implicit solvation）和恒电势方法简介。学完后你可以独立搭建一套"VASP 催化计算（DFT catalysis, free energy diagram）"工作流，并知道每个参数为什么这样选。

---

## 1. 本章目标与适用场景

| 你想做的事 | 本章对应小节 | 配套示例 |
|---|---|---|
| 计算分子在金属/氧化物表面的吸附能 | §4–§7 | [10-co-pt111-adsorption](../examples/10-co-pt111-adsorption/README.md) |
| 比较 top/bridge/hollow 位点、研究覆盖度 | §5 | 同上（`3_co_top`、`4_co_fcc`） |
| 得到吸附物的 ZPE、熵与吉布斯自由能 | §8 | 同上（`5_freq_co_top`、`6_freq_co_gas`） |
| 画 HER/OER/ORR 自由能台阶图、算过电位 | §9 | [10-her-oer-free-energy](../examples/10-her-oer-free-energy/README.md) |
| 用 d 带中心解释吸附强弱趋势 | §10 | `7_slab_dos` + `dband_center.py` |
| 在隐式水溶剂中计算、了解恒电势方法 | §11 | [10-vaspsol-implicit-solvent](../examples/10-vaspsol-implicit-solvent/README.md) |

**前置知识**：结构优化（[第4章](04-relax-eos.md)）、收敛测试（[第3章](03-convergence.md)）、slab 模型与偶极校正（[第9章](09-surface-workfunction.md)）、vdW 修正（[第8章](08-vdw-2d.md)）、DOS（[第5章](05-band-dos.md)）。建模工具见 [Common 第3章](../../Common/tutorials/03-structure-modeling.md) 和 [Common 第4章](../../Common/tutorials/04-ase-pymatgen.md)。过渡态（反应能垒）放在下一章 [第11章](11-neb-dimer.md)。

> [!NOTE]
> 本章所有输入针对 VASP 6.4–6.6（推荐 6.6.1），对 VASP 5.4.4 也基本适用。VASPsol 需要单独打补丁编译（§11）。文中出现的数值输出除特别注明来源外，均为**示意输出，数值仅供参考**，请以你自己的计算为准。

---

## 2. 原理速览

### 2.1 吸附能

最常用的定义（以气相分子为参考）：

$$
E_\mathrm{ads} = E_\mathrm{slab+ads} - E_\mathrm{slab} - E_\mathrm{gas}
$$

- $E_\mathrm{ads}<0$：吸附放热、吸附态更稳定；数值越负吸附越强。也有文献把吸附能定义成正值（结合能 $E_b=-E_\mathrm{ads}$），**写论文时一定要声明符号约定**。
- 解离吸附（如 $\mathrm{H_2}\to 2\mathrm{H^*}$）通常按每个吸附原子计：$E_\mathrm{ads}(\mathrm{H})=E_\mathrm{slab+H}-E_\mathrm{slab}-\tfrac12E_\mathrm{H_2}$。
- 多个吸附物时区分**平均吸附能** $\bar E_\mathrm{ads}=\frac{1}{n}\left[E_{\mathrm{slab}+n\,\mathrm{ads}}-E_\mathrm{slab}-nE_\mathrm{gas}\right]$ 与**微分吸附能** $E_{\mathrm{slab}+n}-E_{\mathrm{slab}+(n-1)}-E_\mathrm{gas}$。
- 三个能量都取 OUTCAR 中最后一个离子步的 `energy(sigma->0)`（外推到零展宽的能量，见[第2章](02-input-output-files.md)），这样金属 slab（大 SIGMA）与分子（小 SIGMA）的能量才能直接相减。

### 2.2 从电子能到吉布斯自由能

DFT 给出 0 K、无零点振动的电子能 $E_\mathrm{DFT}$。在温度 $T$ 下：

$$
G = E_\mathrm{DFT} + E_\mathrm{ZPE} + \int_0^T C_p\,\mathrm{d}T - TS
$$

- **吸附物**：平动、转动被表面"冻结"成受阻振动，通常用**谐振子近似**——把吸附物的 $3N$ 个振动模全部当作谐振动（此时 $\int C_v\,\mathrm{d}T$ 与 $\int C_p\,\mathrm{d}T$ 相同，不含 $pV$ 项）：

$$
E_\mathrm{ZPE}=\sum_i \frac{h\nu_i}{2},\qquad
S_\mathrm{vib}=k_B\sum_i\left[\frac{x_i}{e^{x_i}-1}-\ln\left(1-e^{-x_i}\right)\right],\quad x_i=\frac{h\nu_i}{k_BT}
$$

- **气相分子**：用理想气体模型（平动 + 转动 + 振动 + 电子），熵与压力有关：$S(p)=S^\circ-k_B\ln(p/p^\circ)$。
- 低频模（< 50–100 cm⁻¹）对熵贡献极大且数值噪声大，常见处理是把低于 50 cm⁻¹ 的频率抬高到 50 cm⁻¹（VASPKIT 501 就是这样做的），或改用受阻平动/转动模型。

### 2.3 计算氢电极（CHE）模型

Nørskov 等提出的 CHE 模型（[J. Phys. Chem. B 2004](https://doi.org/10.1021/jp047349j)）利用平衡 $\mathrm{H^+ + e^-} \rightleftharpoons \tfrac12\mathrm{H_2}$ 在 0 V（RHE）下成立，绕开了对溶剂化质子的直接计算：

$$
\mu(\mathrm{H^+}+\mathrm{e^-}) = \tfrac12\mu(\mathrm{H_2}) - eU_\mathrm{RHE}
$$

于是每一步质子-电子转移反应的自由能变化随电势线性移动：氧化方向（OER）每步 $\Delta G_i(U)=\Delta G_i(0)-eU$，还原方向（ORR、HER）每步 $\Delta G_i(U)=\Delta G_i(0)+eU$。用 SHE 标度时再加 pH 校正 $-k_BT\ln 10\cdot\mathrm{pH}\approx -0.059\,\mathrm{pH}$ eV。**极限电势** $U_L$ 是使所有步骤都"下坡"的电势，过电位 $\eta=|U_L-U_\mathrm{eq}|$。

### 2.4 d 带中心模型

Hammer–Nørskov 的 d 带模型认为过渡金属表面的吸附强度主要由表面原子 d 态相对费米能级的位置决定。d 带中心是 d 投影态密度的一阶矩：

$$
\varepsilon_d=\frac{\int E\,n_d(E)\,\mathrm{d}E}{\int n_d(E)\,\mathrm{d}E}\quad(\text{以 }E_F=0\text{ 为零点})
$$

$\varepsilon_d$ 越靠近 $E_F$，反键态越难被填充，吸附通常越强。它是**趋势描述符**，绝对数值依赖积分窗口与未占据带的处理。

### 2.5 隐式溶剂

隐式溶剂把溶剂看成介电连续介质，介电常数在溶质电子密度附近从 1 平滑过渡到 $\varepsilon_b$（水为 78.4），求解广义 Poisson 方程 $\nabla\cdot[\varepsilon(\mathbf r)\nabla\phi]=-4\pi\rho$。VASPsol 还包含空腔形成能（cavitation）项，并可加入线性化 Poisson–Boltzmann 电解质（Debye 屏蔽）。

---

## 3. 整体工作流与示例目录

```text
体相优化得到 a0 (第4章) ──> 切 slab、真空/层数收敛 (第9章)
        │
        ├─ 1_slab        清洁 slab 优化            ──┐
        ├─ 2_co_gas      气相分子（盒子 + Gamma 点） ├─> E_ads (§7)
        ├─ 3_co_top/4_co_fcc  不同位点吸附优化     ──┘
        ├─ 5_freq_co_top / 6_freq_co_gas  频率      ──> ZPE、TS、G (§8)
        ├─ 7_slab_dos    静态 DOS (LORBIT=11)       ──> d 带中心 (§10)
        └─ HER/OER 中间体 + CHE 脚本                ──> 台阶图、过电位 (§9)
                         └─ VASPsol 溶剂化单点/优化  ──> 溶剂化校正 (§11)
```

配套示例 [VASP/examples/10-co-pt111-adsorption](../examples/10-co-pt111-adsorption/README.md) 的文件：

| 路径 | 内容 |
|---|---|
| `1_slab/` | Pt(111) p(2×2) 4 层 slab，底部 2 层固定，a = 3.97 Å |
| `2_co_gas/` | CO 分子，15.0×15.5×16.0 Å 盒子，Gamma 点 |
| `3_co_top/`、`4_co_fcc/` | CO 吸附在 top / fcc hollow 位（1/4 ML） |
| `5_freq_co_top/`、`6_freq_co_gas/` | IBRION=5 频率计算 INCAR/KPOINTS |
| `7_slab_dos/` | 12×12×1 静态 DOS（d 带中心） |
| `ads_energy.py` | 从三个 OUTCAR 计算吸附能 |
| `set_freq_selective.py` | 把 CONTCAR 改成只放开吸附质的 POSCAR |
| `dband_center.py` | 从 DOSCAR 计算 d 带中心、宽度、填充度 |
| `make_potcar.sh`、`job.slurm` | 按 POSCAR 拼接 POTCAR；Slurm 作业脚本 |

> [!IMPORTANT]
> **同一套对比计算中，slab、吸附体系、气相分子必须使用完全相同的 POTCAR、ENCUT、PREC、LASPH、LREAL、GGA/IVDW 设置**，slab 与吸附体系还要用相同的 k 点和晶胞。任何一项不同，能量差都会带入系统误差，而且通常比你想研究的效应还大。

---

## 4. 第一步：清洁 Pt(111) slab

### 4.1 模型

- 晶格常数用**同一泛函**下体相优化得到的值（PBE 下 Pt 约 3.97 Å，近似值；实验值约 3.92 Å）。如果你打开 D3，就要用 PBE-D3 重新优化 a0。
- p(2×2) 表面胞、4 层原子、底部 2 层固定（模拟体相）、上部 2 层放开；z 方向约 18 Å 真空。层数、真空的收敛方法见[第9章](09-surface-workfunction.md)。
- 六角表面晶胞**必须用 Γ 中心 k 网格**，否则可能遇到 `IBZKPT` 报错（见§13）。

POSCAR 的开头（完整文件见示例）：

```text
Pt(111) p(2x2) 4L slab, bottom 2 layers fixed, a=3.97
   1.0
      5.6144278426     0.0000000000     0.0000000000
      2.8072139213     4.8622371394     0.0000000000
      0.0000000000     0.0000000000    24.8762417060
   Pt
   16
Selective dynamics
Direct
    0.0000000000   0.0000000000   0.3617909854   F F F
    ...
    0.0000000000   0.0000000000   0.6382090146   T T T
```

### 4.2 INCAR

```bash
SYSTEM  = Pt(111) p(2x2) 4-layer slab, clean surface relaxation

# ---- electronic part: keep IDENTICAL in slab, slab+CO and CO-gas runs ----
PREC    = Accurate
ENCUT   = 450          # C and O POTCARs have ENMAX ~400 eV; same value in every run
EDIFF   = 1E-6
ALGO    = Fast
NELM    = 150
ISMEAR  = 1            # Methfessel-Paxton smearing for the metal slab
SIGMA   = 0.10         # check: entropy T*S in OUTCAR < ~1 meV per atom
ISPIN   = 1            # Pt and CO are non-magnetic
LREAL   = .FALSE.      # small cell (16-18 atoms): reciprocal-space projectors
LASPH   = .TRUE.       # aspherical gradient corrections inside PAW spheres
# IVDW  = 12           # D3(BJ); if switched on, use it in ALL runs (and for bulk a0)

# ---- dipole correction along the surface normal (asymmetric slab) ----
LDIPOL  = .TRUE.
IDIPOL  = 3
DIPOL   = 0.5 0.5 0.5  # centre of the slab (direct coordinates)

# ---- ionic part ----
IBRION  = 2            # conjugate gradient
ISIF    = 2            # cell fixed: in-plane lattice from the bulk a0
NSW     = 200
EDIFFG  = -0.02        # stop when all free forces < 0.02 eV/A

# ---- output / parallelisation ----
LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4            # adjust to your node (divisor of cores per node)
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ENCUT` | 450 eV | 取体系中最大 ENMAX（C、O 约 400 eV）的 1.1–1.3 倍；吸附能对 ENCUT 的收敛要单独测（§12） |
| `ISMEAR`/`SIGMA` | 1 / 0.1–0.2 | 金属用 Methfessel–Paxton；检查 `grep "entropy T" OUTCAR` 每原子 < 1 meV |
| `LREAL` | `.FALSE.` | 小胞更准；>50 原子的大 slab 可用 `Auto`（但所有对比计算保持一致） |
| `LASPH` | `.TRUE.` | 含 d 元素、使用 vdW-DF/meta-GGA 时尤其需要 |
| `LDIPOL`/`IDIPOL` | `.TRUE.`/3 | 底部固定、顶部弛豫（以及吸附后）的 slab 两侧不对称，有净偶极 |
| `DIPOL` | slab 质心 | 不设时 VASP 自动寻找电荷密度最小的平面；显式给出更稳妥 |
| `EDIFFG` | −0.02 | 吸附能精度要求 ~0.01 eV 时足够；频率计算前建议收紧到 −0.01 |
| `ISIF` | 2 | slab 计算**绝不能**变胞 |

> [!TIP]
> 偶极校正会让 SCF 变慢。如果 OSZICAR 中能量来回振荡，可加 `AMIN = 0.01`（抑制电荷在 slab 两侧来回晃动）并增大 `NELM`。这也是官方 LDIPOL 页面给出的建议。

### 4.3 运行与检查

```bash
cd VASP/examples/10-co-pt111-adsorption
export VASP_PP_PATH=/path/to/your/paw/library     # 你的赝势库，仓库不提供 POTCAR
bash make_potcar.sh 1_slab 2_co_gas 3_co_top 4_co_fcc
cd 1_slab && sbatch ../job.slurm && cd ..

# 结束后
grep "reached required accuracy" 1_slab/OUTCAR
grep "energy(sigma->0)" 1_slab/OUTCAR | tail -1
grep "entropy T" 1_slab/OUTCAR | tail -1          # T*S / 18 原子 应 << 1 meV/原子
grep "dipolmoment" 1_slab/OUTCAR | tail -1        # 偶极矩（LDIPOL 开启时输出）
```

通用的作业健康检查可以用仓库脚本 [`check_vasp.py`](../scripts/check_vasp.py)（`python ../scripts/check_vasp.py -h` 查看用法）。

---

## 5. 吸附位点与覆盖度

### 5.1 fcc(111) 表面的四种高对称位

| 位点 | 位置 | 配位 | 典型初始高度（C 原子到表面层） |
|---|---|---|---|
| top（ontop） | 正上方一个表面原子 | 1 | ~1.85 Å |
| bridge | 两个表面原子中点 | 2 | ~1.5 Å |
| fcc hollow | 三个表面原子中心，下方第二层**无**原子 | 3 | ~1.3 Å |
| hcp hollow | 三个表面原子中心，下方第二层**有**原子 | 3 | ~1.3 Å |

用 ASE 生成（这就是示例中 POSCAR 的生成方式）：

```python
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.io import write

slab = fcc111("Pt", size=(2, 2, 4), a=3.97, vacuum=9.0)       # p(2x2), 4 layers
slab.set_constraint(FixAtoms(mask=[t > 2 for t in slab.get_tags()]))  # tags: 1 = top layer
for site, h in [("ontop", 1.85), ("bridge", 1.50), ("fcc", 1.30), ("hcp", 1.30)]:
    s = slab.copy()
    add_adsorbate(s, "C", h, site)            # C down
    add_adsorbate(s, "O", h + 1.15, site)     # C-O ~1.15 A
    write(f"POSCAR_CO_{site}", s, format="vasp", direct=True, sort=False)
```

对于低对称表面（台阶、合金、氧化物），手工找位点容易遗漏，可以用 pymatgen 的 `AdsorbateSiteFinder`（Delaunay 三角剖分找 ontop/bridge/hollow）：

```python
from pymatgen.core import Molecule, Structure
from pymatgen.analysis.adsorption import AdsorbateSiteFinder   # 新版实际位于 pymatgen.core.adsorption
from pymatgen.io.vasp import Poscar

slab = Structure.from_file("1_slab/CONTCAR")
asf = AdsorbateSiteFinder(slab, height=0.9)          # 距最高原子 0.9 A 内算"表面原子"
sites = asf.find_adsorption_sites(distance=1.8)      # dict: ontop / bridge / hollow / all
for k in ("ontop", "bridge", "hollow"):
    print(k, len(sites[k]))

co = Molecule(["C", "O"], [[0, 0, 0], [0, 0, 1.15]])
structs = asf.generate_adsorption_structures(co, repeat=[1, 1, 1],
                                             find_args={"distance": 1.8})
for i, s in enumerate(structs):
    Poscar(s).write_file(f"POSCAR_ads_{i:02d}")
```

> [!WARNING]
> `generate_adsorption_structures` 会对**对称不等价**位点逐一放分子，但不会替你检查分子取向（端基 C-down/O-down、平躺）和吸附后的固定层设置。写出的结构需要逐一在 VESTA 里检查，并确认 Selective dynamics 与清洁 slab 一致。

### 5.2 覆盖度与横向相互作用

覆盖度 $\theta$ = 吸附物数 / 表面原子数：p(2×2) 中 1 个 CO 是 1/4 ML，p(3×3) 是 1/9 ML，($\sqrt3\times\sqrt3$)R30° 是 1/3 ML。吸附物间存在**横向相互作用**（偶极–偶极、通过金属的电子相互作用），所以：

1. 比较不同表面/催化剂时必须在**相同覆盖度**下；
2. 研究低覆盖极限时，逐步增大表面胞（2×2 → 3×3 → 4×4），直到 $E_\mathrm{ads}$ 变化 < 0.02–0.05 eV；
3. 真实反应条件下的覆盖度可能很高（如 CO 氧化），这时要构建覆盖度-吸附能关系或做微观动力学。

### 5.3 一个经典教训：CO/Pt(111) 之谜

实验上低覆盖 CO 在 Pt(111) 上吸附于 **top 位**，而 LDA/GGA（包括 PBE）通常给出 **fcc hollow 更稳定**——这就是著名的"CO/Pt(111) puzzle"（Feibelman 等，[J. Phys. Chem. B 2001](https://doi.org/10.1021/jp002302t)）。原因与 GGA 低估 CO 的 HOMO–LUMO 间隙、高估 2π\* 反馈键有关。示例中同时给出 `3_co_top` 和 `4_co_fcc`，你可以亲自验证 PBE 的位点偏好。

> [!NOTE]
> 这提醒我们：GGA 的**绝对吸附能**常有 0.2–0.5 eV 的误差（通常高估化学吸附强度），但同一泛函下的**趋势和差值**更可靠。对位点偏好敏感的问题，可以比较 RPBE、BEEF-vdW（需 libbeef，见[第8章](08-vdw-2d.md)）或杂化泛函（[第6章](06-hybrid-metagga.md)）的结果。

### 5.4 vdW 修正要不要加？

- 物理吸附、芳香分子平躺吸附、大分子：**必须**加（`IVDW = 12` D3(BJ) 或 optB86b-vdW 等非局域泛函，[第8章](08-vdw-2d.md)）。
- 小分子化学吸附（CO、O、H）：D3 会使吸附更强 0.1–0.3 eV 量级，常被用于更接近实验的绝对值。
- 一旦使用，**体相 a0、slab、吸附体系、气相分子全部都要用**。示例中 `IVDW` 被注释掉，你可以打开它做对比练习。

---

## 6. 第二步：气相分子参考能量

### 6.1 设置要点

| 要点 | 做法 | 原因 |
|---|---|---|
| 盒子 | ≥ 15 Å，最好三个边长略有不同（如 15.0×15.5×16.0） | 消除周期像间相互作用；不等边打破立方对称，避免 O2 等简并轨道占据出错 |
| k 点 | 只用 Γ 点，配合 `vasp_gam` | 孤立分子无色散；`vasp_gam` 快约 2 倍 |
| 展宽 | `ISMEAR = 0`，`SIGMA = 0.01`（≤ 0.05） | 分子有 HOMO–LUMO 能隙，大展宽会导致分数占据 |
| 自旋 | O2 三重态 `ISPIN = 2`，`MAGMOM = 2*1`；自由基（OH、NO）也要自旋极化 | 检查 `grep "number of electron" OUTCAR` 中 magnetization 为 2.0 |
| 偶极 | 极性分子可加 `IDIPOL = 4`（只做能量校正） | 15 Å 盒子中 CO 的偶极作用很小，通常可忽略 |
| 优化 | `ISIF = 2`，绝不能变胞 | |

CO 分子的 INCAR（`2_co_gas/INCAR`）：

```bash
SYSTEM  = CO molecule in a 15.0 x 15.5 x 16.0 A box - relaxation

# ---- electronic part: same PREC/ENCUT/LASPH/LREAL (and IVDW) as the slab ----
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-6
ALGO    = Fast
NELM    = 150
ISMEAR  = 0            # Gaussian smearing for a molecule (HOMO-LUMO gap)
SIGMA   = 0.01
ISPIN   = 1            # CO is closed shell (O2 would need ISPIN = 2, triplet)
LREAL   = .FALSE.
LASPH   = .TRUE.
# IVDW  = 12

# ---- ionic part ----
IBRION  = 2
ISIF    = 2            # never change the box of an isolated molecule
NSW     = 50
EDIFFG  = -0.02

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

KPOINTS：

```text
Gamma point only (isolated molecule; run with vasp_gam)
0
Gamma
1 1 1
0 0 0
```

运行与检查：

```bash
cd 2_co_gas && VASP_EXE=vasp_gam sbatch --export=ALL ../job.slurm && cd ..
grep "energy(sigma->0)" 2_co_gas/OUTCAR | tail -1
python -c "from ase.io import read; a=read('2_co_gas/CONTCAR'); print('d(C-O) = %.3f A' % a.get_distance(0,1))"
```

PBE 优化的 C–O 键长约 1.14 Å（实验 1.128 Å，近似比较）。

### 6.2 O2 的已知误差与"绕开 O2"

GGA 对 O2 的结合能误差较大（PBE 明显高估 O=O 结合）。因此 OER/ORR 的 CHE 计算中**一般不直接用 DFT 的 O2 能量**，而是用水的实验生成自由能：

$$
\mathrm{2H_2O(l)\to O_2(g) + 2H_2(g)},\qquad \Delta G = 4.92\ \mathrm{eV}\ (298\ \mathrm{K})
$$

即 $G(\mathrm{O_2}) = 2G(\mathrm{H_2O}) - 2G(\mathrm{H_2}) + 4.92$ eV。示例 [10-her-oer-free-energy](../examples/10-her-oer-free-energy/README.md) 中保留了 `O2_gas/`，用于让你自己比较"DFT O2"与"水参考 O2"的差别。

---

## 7. 第三步：吸附结构优化与吸附能

### 7.1 优化

`3_co_top/INCAR` 与 slab 完全相同，只改了 `SYSTEM`。POSCAR 中 C、O 放在 top 位上方，Pt 的固定方式与清洁 slab 一致。

```bash
for d in 3_co_top 4_co_fcc; do (cd $d && sbatch ../job.slurm); done
```

> [!TIP]
> 更省时的做法：先把清洁 slab 优化好，再把吸附物放到**优化后的 slab（1_slab/CONTCAR）**上，而不是放到理想体相截断的 slab 上。示例为了让各目录可以并行提交，直接用了理想 slab；生产计算中建议改用 CONTCAR。

优化结束后检查：

```bash
grep "reached required accuracy" 3_co_top/OUTCAR 4_co_fcc/OUTCAR
# 吸附物是否跑位：比较 C 原子与最近 Pt 的距离和配位数
python - <<'EOF'
from ase.io import read
for d in ("3_co_top", "4_co_fcc"):
    a = read(f"{d}/CONTCAR")
    iC = [i for i, s in enumerate(a.get_chemical_symbols()) if s == "C"][0]
    dPt = sorted(a.get_distance(iC, j, mic=True) for j, s in enumerate(a.get_chemical_symbols()) if s == "Pt")
    print(d, "C-O = %.3f" % a.get_distance(iC, iC + 1), " shortest C-Pt:", " ".join("%.2f" % x for x in dPt[:3]))
EOF
```

- top 位：1 个 Pt–C ≈ 1.85 Å，其余 > 3 Å；fcc 位：3 个 Pt–C ≈ 2.1 Å。
- 如果 top 位的 CO 优化后有 2–3 个近邻 Pt，说明它**滑到了 bridge/hollow**——这本身是有用的信息（top 位在该泛函下不是极小点），但不能再把它当作"top 位吸附能"。需要固定 C 的 x、y（Selective dynamics `F F T`）才能得到约束的 top 位能量，并在论文中注明。

### 7.2 计算吸附能

```bash
python ads_energy.py --slab 1_slab --ads 3_co_top 4_co_fcc --gas 2_co_gas
```

输出形如（**示意输出，数值仅供参考**）：

```text
E(slab)  1_slab                -95.xxxxxx eV  (relaxed)
E(gas)   2_co_gas              -14.8xxxxx eV  (relaxed)
------------------------------------------------------------------
E(ads)   3_co_top             -112.xxxxxx eV  (relaxed)
  -> E_ads                             -1.7xx eV
E(ads)   4_co_fcc             -112.xxxxxx eV  (relaxed)
  -> E_ads                             -1.8xx eV
```

也可以手工：

```bash
for d in 1_slab 2_co_gas 3_co_top; do
  printf "%-10s %s\n" $d "$(grep 'energy(sigma->0)' $d/OUTCAR | tail -1 | awk '{print $NF}')"
done
```

> [!NOTE]
> 文献中 PBE 计算的 1/4 ML CO/Pt(111) 吸附能大多落在约 −1.6 到 −1.9 eV（不同文献的层数、k 点、POTCAR 不同，这只是近似范围）；实验吸附热（量热/TPD，低覆盖）通常报道在约 1.3–1.7 eV（方法与覆盖度不同而有差异）。PBE 结果整体偏强，正是 §5.3 讨论的 GGA 过度结合问题。

### 7.3 吸附诱导的电荷与功函数变化

吸附后可以顺带分析：
- 功函数变化 $\Delta\Phi$（`LVHAR = .TRUE.` 输出 LOCPOT，用 [`locpot_workfunction.py`](../scripts/locpot_workfunction.py) 处理，见[第9章](09-surface-workfunction.md)）；
- 差分电荷 $\Delta\rho=\rho(\mathrm{slab+CO})-\rho(\mathrm{slab})-\rho(\mathrm{CO})$ 与 Bader 电荷（[第16章](16-charge-bonding.md)，脚本 [`chgdiff.py`](../scripts/chgdiff.py)）。

---

## 8. 第四步：频率计算——ZPE、熵与自由能校正

### 8.1 为什么用 IBRION = 5 并且只放开吸附质

- `IBRION = 5`：对每个**可动**原子的每个笛卡尔方向做 ±POTIM 位移（NFREE=2，中心差分），由力的差分构造 Hessian 并对角化。
- 官方文档明确说明 **Selective dynamics 目前只被 `IBRION = 5` 支持**（`IBRION = 6` 利用对称性，但不支持选择性动力学）。
- 吸附质与金属振动频率相差很大（金属声子 < 250 cm⁻¹），把 slab 固定是常用且合理的近似；计算量只与放开原子数 $N_\mathrm{free}$ 有关：$2\times3N_\mathrm{free}$ 次 SCF。

### 8.2 准备 POSCAR

```bash
mkdir -p 5_freq_co_top
python set_freq_selective.py 3_co_top/CONTCAR -e C O -o 5_freq_co_top/POSCAR
#  -> wrote 5_freq_co_top/POSCAR: 18 atoms, 2 free atoms = 6 degrees of freedom
#     IBRION=5, NFREE=2 -> 12 displaced structures, 6 vibrational modes in OUTCAR
bash make_potcar.sh 5_freq_co_top
```

`set_freq_selective.py` 只修改 Selective dynamics 标志，坐标原样保留；也可以按序号（`-i 17-18`）或高度（`--zmin 15.0`）选择原子，`-h` 查看全部选项。

### 8.3 INCAR

```bash
SYSTEM  = CO on Pt(111) top - vibrational frequencies of the adsorbate (IBRION=5)

# ---- electronic part: identical to the relaxation, but tighter EDIFF ----
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-7         # accurate forces are needed for finite differences
ALGO    = Fast
NELM    = 150
ISMEAR  = 1
SIGMA   = 0.10
ISPIN   = 1
LREAL   = .FALSE.
LASPH   = .TRUE.
LDIPOL  = .TRUE.
IDIPOL  = 3
DIPOL   = 0.5 0.5 0.5

# ---- finite-difference Hessian ----
IBRION  = 5            # all free Cartesian directions; honours selective dynamics
POTIM   = 0.015        # displacement step in A
NFREE   = 2            # central differences: +POTIM and -POTIM
NSW     = 1

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `IBRION` | 5 | 有限差分 Hessian；6 会利用对称性但忽略 Selective dynamics |
| `POTIM` | 0.015 Å | 官方经验值；对 H 等轻原子可试 0.01–0.02 检查稳定性 |
| `NFREE` | 2 | 中心差分；4 更准但贵一倍；1（单侧）官方强烈不推荐 |
| `EDIFF` | 1E-7（至少 1E-6） | 力的差分对 SCF 噪声很敏感 |
| `NSW` | 1 | VASP ≥ 5.1 中可省略，保留无害 |

气相 CO 的频率（`6_freq_co_gas`）：把 `2_co_gas/CONTCAR` 复制为 POSCAR，所有原子放开（不加 Selective dynamics）即可。

### 8.4 读取频率

```bash
grep "THz" 5_freq_co_top/OUTCAR | head -6
```

（**示意输出，数值仅供参考**）

```text
   1 f  =   61.xxxxxx THz   38x.xxxxxx 2PiTHz 20xx.xxxxxx cm-1   25x.xxxxxx meV
   2 f  =   1x.xxxxxx THz    xx.xxxxxx 2PiTHz  4xx.xxxxxx cm-1    5x.xxxxxx meV
   ...
   6 f  =    1.xxxxxx THz     x.xxxxxx 2PiTHz   5x.xxxxxx cm-1     x.xxxxxx meV
```

- `f` 为实频，`f/i` 为虚频；吸附物的 6 个模包括 C–O 伸缩（~2000–2100 cm⁻¹）、Pt–C 伸缩、两个受阻转动、两个受阻平动（几十 cm⁻¹）。
- OUTCAR 中频率**会打印两遍**（"Eigenvectors and eigenvalues of the dynamical matrix" 与 "Eigenvectors after division by SQRT(mass)"），用 `grep` 统计时注意只取前 $3N_\mathrm{free}$ 行。
- 极小点不应有明显虚频；出现 < 50 cm⁻¹ 的小虚频通常是结构没优化到位或 SCF 精度不足：收紧 `EDIFFG = -0.01`、`EDIFF = 1E-7` 重新优化后再算。
- 气相 CO 的 6 个模中只有 1 个是真振动（线形分子 3N−5 = 1），其余 5 个接近 0（可能带 f/i），属于平动/转动，热力学计算时要剔除。

振动模式可视化：[VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF) 把 OUTCAR 振动模式转成 XSF 动画，在 VESTA/XCrySDen 中查看。

### 8.5 计算 ZPE、TS 与 G 校正

**方法 1：VASPKIT**（在频率计算目录中运行，交互输入温度等）：

| 功能编号 | 用途 | 模型 |
|---|---|---|
| `501` | 吸附物热校正 | 只含振动，低于 50 cm⁻¹ 的频率按 50 cm⁻¹ 计，忽略 pV |
| `502` | 气相分子热校正 | 平动 + 转动 + 振动 + 电子（与 Gaussian 的处理等价），需输入温度、压力、自旋多重度 |
| `503` | d 带中心（见§10） | |

**方法 2：示例脚本 `calc_dG.py thermo`**（不依赖交互，便于批处理）：

```bash
cd ../10-her-oer-free-energy
python calc_dG.py thermo ../10-co-pt111-adsorption/5_freq_co_top -T 298.15          # 吸附物：谐振子
python calc_dG.py thermo ../10-co-pt111-adsorption/6_freq_co_gas --gas \
       --geometry linear --symmetry 1 --spin 0 \
       --structure ../10-co-pt111-adsorption/6_freq_co_gas/POSCAR                    # 气相：理想气体（需要 ASE）
```

`--gas` 模式内部调用 `ase.thermochemistry.IdealGasThermo`；吸附物模式自己实现谐振子公式（§2.2），默认把低于 50 cm⁻¹ 的频率抬高到 50 cm⁻¹（`--cutoff 0` 关闭）。

**方法 3：直接用 ASE**：

```python
from ase.thermochemistry import HarmonicThermo
import numpy as np
cm2eV = 1.239841984e-4
freqs_cm = np.array([2050., 480., 410., 405., 60., 55.])   # 示意数值：换成你 OUTCAR 中的实频
thermo = HarmonicThermo(vib_energies=freqs_cm * cm2eV, potentialenergy=0.0)
print("G_corr(298 K) = %.3f eV" % thermo.get_helmholtz_energy(298.15, verbose=False))
```

**方法 4**：[VaspGibbs](https://github.com/ftherrien/VaspGibbs) 等脚本也能从 OUTCAR 一键得到 ZPE、熵和吉布斯自由能。

### 8.6 常用的气相校正值

自己算是首选。作为数量级参考，CHE 文献中常见的 298 K 近似值（如 Nørskov 2004、Valdés 2008 等采用的数值，属于近似量）：

| 物种 | ZPE (eV) | TS (eV) | 备注 |
|---|---|---|---|
| H2 (g, 1 bar) | ~0.27 | ~0.40–0.41 | |
| H2O (l) | ~0.56 | ~0.67 | 用 0.035 bar（约 300 K 时液态水的饱和蒸气压）下的气相水代替液态水 |
| OH\* / O\* / OOH\* | ~0.35 / ~0.05–0.07 / ~0.40–0.45 | 很小 | 依赖表面，必须自己算 |

> [!WARNING]
> 不要把不同来源、不同温度/压力约定的校正值混用。尤其是液态水：要么用 0.035 bar 的气相水（与液态水化学势相等），要么直接用实验液态水的标准熵，二者不能叠加。

---

## 9. 计算氢电极：HER / OER / ORR 自由能台阶图

### 9.1 HER：ΔG_H\*

$$
\Delta G_\mathrm{H^*} = \Delta E_\mathrm{H^*} + \Delta E_\mathrm{ZPE} - T\Delta S,\qquad
\Delta E_\mathrm{H^*} = E_\mathrm{slab+H}-E_\mathrm{slab}-\tfrac12E_\mathrm{H_2}
$$

Nørskov 等（[J. Electrochem. Soc. 2005](https://doi.org/10.1149/1.1856988)）对金属表面给出 $\Delta E_\mathrm{ZPE}-T\Delta S\approx 0.24$ eV，所以常写 $\Delta G_\mathrm{H^*}\approx\Delta E_\mathrm{H^*}+0.24$ eV。$|\Delta G_\mathrm{H^*}|\to 0$ 的表面处于火山图顶点（Pt 接近这一点）。

示例流程（H 在 Pt(111) fcc 位，1/4 ML）：

```bash
cd VASP/examples/10-her-oer-free-energy
bash ../10-co-pt111-adsorption/make_potcar.sh H_fcc H2_gas
(cd H_fcc  && sbatch ../../10-co-pt111-adsorption/job.slurm)
(cd H2_gas && VASP_EXE=vasp_gam sbatch --export=ALL ../../10-co-pt111-adsorption/job.slurm)
# 清洁 slab 直接复用 ../10-co-pt111-adsorption/1_slab（设置完全相同）
python calc_dG.py her --slab ../10-co-pt111-adsorption/1_slab --slabH H_fcc --h2 H2_gas
```

更严格的做法是用自己的频率替代 0.24 eV：

```bash
python calc_dG.py thermo freq_H -T 298.15                     # 得到 G_corr(H*)
python calc_dG.py thermo freq_H2 --gas --geometry linear --symmetry 2 --spin 0 --structure freq_H2/POSCAR
python calc_dG.py her --slab ... --slabH ... --h2 ... --gcorr-h <G_corr(H*)> --gcorr-h2 <G_corr(H2)>
```

（`freq_H`、`freq_H2` 目录按 §8 的方法用 `INCAR_freq_ads`、`INCAR_freq_gas` 自行创建。）

> [!NOTE]
> ΔG_H\* 依赖覆盖度。Pt(111) 在 HER 条件下 H 覆盖度较高，文献常比较 1/4 ML 与 1 ML 的结果；做趋势比较时要统一覆盖度定义（平均还是微分）。

### 9.2 OER 四电子机理（酸性/RHE 标度）

$$
\begin{aligned}
&(1)\quad \mathrm{* + H_2O \to OH^* + H^+ + e^-} \\
&(2)\quad \mathrm{OH^* \to O^* + H^+ + e^-} \\
&(3)\quad \mathrm{O^* + H_2O \to OOH^* + H^+ + e^-} \\
&(4)\quad \mathrm{OOH^* \to * + O_2 + H^+ + e^-}
\end{aligned}
$$

在 $U=0$ V（RHE）：

$$
\begin{aligned}
\Delta G_1 &= G_\mathrm{OH^*}-G_\mathrm{*}-G_\mathrm{H_2O}+\tfrac12G_\mathrm{H_2}\\
\Delta G_2 &= G_\mathrm{O^*}-G_\mathrm{OH^*}+\tfrac12G_\mathrm{H_2}\\
\Delta G_3 &= G_\mathrm{OOH^*}-G_\mathrm{O^*}-G_\mathrm{H_2O}+\tfrac12G_\mathrm{H_2}\\
\Delta G_4 &= 4.92\ \mathrm{eV}-(\Delta G_1+\Delta G_2+\Delta G_3)
\end{aligned}
$$

- OER 极限电势 $U_L^\mathrm{OER}=\max_i\Delta G_i/e$，过电位 $\eta_\mathrm{OER}=U_L-1.23$ V；
- ORR 是逆过程：$U_L^\mathrm{ORR}=\min_i\Delta G_i/e$，$\eta_\mathrm{ORR}=1.23-U_L$ V；
- 理想催化剂四步都等于 1.23 eV。

**标度关系**：在许多表面上 $\Delta G_\mathrm{OOH^*}-\Delta G_\mathrm{OH^*}\approx 3.2$ eV（Man 等，[ChemCatChem 2011](https://doi.org/10.1002/cctc.201000397)），而理想值是 2.46 eV，由此得到 OER 的理论最小过电位约 $(3.2-2.46)/2\approx0.37$ V。打破标度关系（双位点、限域、质子-电子解耦步骤等）是 OER 催化剂设计的核心思路之一。

### 9.3 用 `calc_dG.py oer` 计算并作图

把各中间体的能量（数字或计算目录）和 G 校正写进 JSON（格式见 [oer_input_demo.json](../examples/10-her-oer-free-energy/oer_input_demo.json)，其中数值是**虚构的演示数字**）：

```json
{
  "E":      {"slab": "slab", "OH": "OH_top", "O": "O_fcc", "OOH": "OOH_top",
             "H2O": "H2O_gas", "H2": "H2_gas"},
  "G_corr": {"slab": 0.0, "OH": 0.30, "O": 0.05, "OOH": 0.35, "H2O": -0.01, "H2": -0.04}
}
```

```bash
python calc_dG.py oer my_oer.json --plot oer.png            # U = 0 V
python calc_dG.py oer my_oer.json --U 1.23                  # 平衡电势下各步
```

用演示 JSON 运行的输出（**虚构数据，仅演示脚本功能**）：

```text
Descriptors (vs H2O/H2, U = 0 V_RHE):
  dG_OH* =  0.900   dG_O* =  2.450   dG_OOH* =  4.100 eV
  dG_OOH* - dG_OH* =  3.200 eV   (scaling relation: ~3.2 eV on many surfaces)
OER steps at U = 0.00 V_RHE:
  dG1 =  0.900 eV   * + H2O -> OH* + (H+ + e-)
  dG2 =  1.550 eV   OH* -> O* + (H+ + e-)
  dG3 =  1.650 eV   O* + H2O -> OOH* + (H+ + e-)
  dG4 =  0.820 eV   OOH* -> * + O2 + (H+ + e-)
OER: limiting potential U_L = 1.650 V, overpotential eta = 0.420 V (step 3)
ORR: limiting potential U_L = 0.820 V, overpotential eta = 0.410 V (step 4, reverse)
```

`--plot` 会画出 U = 0、1.23 V 和 U_L 三条台阶线。

### 9.4 OER/ORR 计算的实务要点

1. **中间体逐一做位点搜索**：OH\*、O\*、OOH\* 的最稳定位点往往不同（Pt(111) 上 O\* 偏好 fcc，OH\* 偏好 top/bridge）。
2. **OOH\* 容易解离**成 O\* + OH\*，优化后检查 O–O 键（~1.4–1.5 Å 为完好）。
3. **溶剂化**：水对 OH\*、OOH\* 有显著氢键稳定作用。常见做法是加显式水层，或用 VASPsol（§11），或采用文献中的经验溶剂化校正（不同文献数值不同，务必注明出处）。
4. **自旋**：氧化物表面（如 IrO2、RuO2、含过渡金属的钙钛矿）要考虑自旋极化与 DFT+U（[第7章](07-magnetism-u-soc.md)）；U 值不同的能量不能互相比较。
5. **碱性条件**：在 RHE 标度下，CHE 推导出的 ΔG 与 pH 无关（只是参考电极换了）；用 SHE 标度时每步加 −0.059·pH eV。

---

## 10. d 带中心

### 10.1 静态 DOS 计算

在优化好的 slab 上做一次 `LORBIT = 11` 的静态计算（`7_slab_dos`），k 点加密到 12×12×1：

```bash
SYSTEM  = Pt(111) relaxed slab - static run for the projected DOS (d-band centre)

PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-6
ALGO    = Normal
NELM    = 150
ISMEAR  = -5           # tetrahedron method with Bloechl corrections (Gamma mesh, >= 4 k-points)
ISPIN   = 1
LREAL   = .FALSE.
LASPH   = .TRUE.
LDIPOL  = .TRUE.
IDIPOL  = 3
DIPOL   = 0.5 0.5 0.5

IBRION  = -1           # no ionic update
NSW     = 0

LORBIT  = 11           # site- and lm-projected DOS in DOSCAR / PROCAR
NEDOS   = 3001         # fine energy grid
EMIN    = -15
EMAX    = 10

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

```bash
cp 1_slab/CONTCAR 7_slab_dos/POSCAR && bash make_potcar.sh 7_slab_dos
(cd 7_slab_dos && sbatch ../job.slurm)
```

### 10.2 计算 ε_d

示例 slab 中原子 13–16 是表面层（最上层），9–12 是次表面层：

```bash
python dband_center.py 7_slab_dos/DOSCAR -a 13-16                # 全能量窗口
python dband_center.py 7_slab_dos/DOSCAR -a 13-16 -w -10 2 --dat pdos_d_surface.dat
```

输出 ε_d（一阶矩）、W_d（二阶矩开方，带宽）、填充度与窗口内 d 电子数；自旋极化体系分别给出 up/down/sum。也可以用：

- VASPKIT `503`：逐原子输出 d 带中心，可自定义积分窗口；
- 仓库脚本 [`plot_dos.py`](../scripts/plot_dos.py)：画 PDOS 并计算 d 带中心（`python ../scripts/plot_dos.py -h` 查看用法）；
- VTST 脚本 `split_dos` + `dosanalyze.pl`。

> [!IMPORTANT]
> d 带中心的数值依赖：①积分窗口（只到 E_F 还是包含未占据部分）；②NBANDS（未占据 d 态是否完整）；③NEDOS 与展宽方法；④选取哪些原子。比较不同表面时这四项必须统一，并在文中写明。它适合解释**同一类过渡金属/合金表面**的吸附趋势，对氧化物、单原子催化剂等未必适用。

---

## 11. 隐式溶剂：VASPsol

### 11.1 安装要点（详细编译见[第1章](01-setup-compile.md)）

VASPsol 由 Hennig 组开发（[henniggroup/VASPsol](https://github.com/henniggroup/VASPsol)），**不是**官方 VASP 的一部分。对 VASP 6.x，其 README 给出的步骤概括为：

1. 把 VASPsol 的 `src/solvation.F` 复制到 VASP 的 `src/`；
2. 确认 `src/.objects` 中 `solvation.o` 位于 `pot.o` 之前；
3. 在 `makefile.include` 的 `CPP_OPTIONS` 中加入 `-Dsol_compat`；
4. 按 README 应用仓库中针对 VASP 6 的补丁（提供 FERMI_SHIFT 等功能），然后正常编译。

> [!NOTE]
> 原 VASPsol 仓库已声明不再开发，新功能（非线性响应、非局域空腔等）转到 VASPsol++（[VASPsol/VASPsol](https://github.com/VASPsol/VASPsol)，论文 J. Chem. Phys. 159, 234117 (2023)）。VASPsol++ 的关键词与本节不完全相同，使用前请阅读其文档。

### 11.2 两步计算法

VASPsol 官方用法要求：**先做真空计算并保存 WAVECAR，再以 `ISTART = 1` 读入，打开 `LSOL`**。

```bash
SYSTEM  = CO/Pt(111) top - step 2: VASPsol implicit water, single point
# Requires a VASP binary compiled with the VASPsol patch (see chapter 1 / 10).

PREC    = Accurate
ENCUT   = 450          # must equal the vacuum run
EDIFF   = 1E-6
ALGO    = Fast
NELM    = 300          # solvent runs often need more SCF steps
ISMEAR  = 1
SIGMA   = 0.10
ISPIN   = 1
LREAL   = .FALSE.
LASPH   = .TRUE.

ISTART  = 1            # start from the vacuum WAVECAR

# ---- VASPsol (values below are the VASPsol defaults for water) ----
LSOL    = .TRUE.       # switch on the implicit solvent
EB_K    = 78.4         # relative permittivity of the bulk solvent (water)
TAU     = 0.000525     # effective surface tension of the cavity; TAU = 0 drops cavitation
# NC_K    = 0.0025     # cut-off charge density that defines the cavity (default, water)
# SIGMA_K = 0.6        # width of the dielectric transition (default, water)
# LRHOB   = .TRUE.     # write the bound charge to RHOB (slow, every SCF step)

IBRION  = -1
NSW     = 0

LWAVE   = .TRUE.
LCHARG  = .FALSE.
NCORE   = 4
```

| 关键词 | 默认值（VASPsol 源码） | 说明 |
|---|---|---|
| `LSOL` | `.FALSE.` | 打开隐式溶剂 |
| `EB_K` | 78.4 | 溶剂相对介电常数（水） |
| `TAU` | 5.25E-4 | 空腔表面张力参数；设 0 则不计空腔能。空腔项需要更高 ENCUT 才能收敛 |
| `NC_K` | 0.0025 | 决定空腔边界的电子密度阈值（对水拟合） |
| `SIGMA_K` | 0.6 | 介电函数过渡宽度（对水拟合） |
| `LAMBDA_D_K` | 不设（无电解质） | Debye 长度（Å），设置后启用线性化 Poisson–Boltzmann 电解质模型 |
| `LRHOB` | `.FALSE.` | 输出束缚电荷 RHOB（每个 SCF 步都写，慢） |

> [!WARNING]
> `NC_K`、`SIGMA_K`、`TAU` 是**针对水**拟合的参数。换成其他溶剂只改 `EB_K` 是不够的，需要参考文献重新拟合这几个参数。

运行示例（[10-vaspsol-implicit-solvent](../examples/10-vaspsol-implicit-solvent/README.md)）：

```bash
cd VASP/examples/10-vaspsol-implicit-solvent
bash ../10-co-pt111-adsorption/make_potcar.sh .
bash run_vaspsol.sh            # vac/ -> sol/，最后打印 E(solvent) - E(vacuum)
bash run_vaspsol.sh lpb        # 额外进行 LAMBDA_D_K = 3.0 的电解质计算
```

### 11.3 溶剂化对吸附能的影响

$$
\Delta E_\mathrm{ads}^\mathrm{sol}=E^\mathrm{sol}_\mathrm{slab+ads}-E^\mathrm{sol}_\mathrm{slab}-E^\mathrm{sol(or\ gas)}_\mathrm{molecule}
$$

- 气相反应物（CO、O2、H2）通常仍取**真空**能量；液相物种（H2O、离子）可取溶剂中的能量——参考态选择要与实验条件一致，并在文中说明。
- 也可以在溶剂中做结构优化（`IBRION = 2` + `LSOL = .TRUE.`），VASPsol 提供溶剂化力；但 SCF 更难收敛，建议从真空优化好的结构与 WAVECAR 开始。
- 对 OH\*、OOH\* 等形成氢键的吸附物，隐式溶剂往往低估稳定化（无方向性氢键）；显式水 + 隐式溶剂的混合模型更可靠。
- 社区常见做法是溶剂计算中不再叠加 `LDIPOL`（溶剂已经屏蔽了真空中的偶极场）；若必须同时使用，请先做对照测试。

### 11.4 恒电势方法简介

CHE 模型中电极电势只通过 $-eU$ 进入，表面本身是电中性的。真实电极在给定电势下带有表面电荷，**恒电势（巨正则）DFT** 尝试直接描述这一点：

1. 在隐式电解质（VASPsol 的 LPB，`LAMBDA_D_K`）中用 `NELECT` 改变体系电子数 $n$（电解质中的反离子会补偿电荷）；
2. 计算功函数/费米能级 $\mu_e(n)$，换算为电极电势 $U=-\mu_e/e-\Phi_\mathrm{SHE}$（$\Phi_\mathrm{SHE}$ 常取 IUPAC 推荐的 4.44 V，部分 VASPsol 文献会做自己的标定）；VASPsol 输出的 `FERMI_SHIFT` 要加到 E_F 上，使电解质深处的电势为零；
3. 对每个 $n$ 计算巨正则能 $\Omega(n)=E(n)-n\mu_e(n)$（注意 VASPsol 文档中关于电势参考移动时能量要加 $QV$ 修正的说明），拟合 $\Omega(U)$（通常近似为抛物线），得到任意电势下的反应能/能垒。

这类方法的原理与算法可参考 Sundararaman 等的综述（[J. Chem. Phys. 2017](https://doi.org/10.1063/1.4978411)）。实现细节（参考电势、电荷补偿方式、电容拟合）差异很大，入门阶段建议先用 CHE 得到趋势，再针对关键步骤做恒电势校验。

---

## 12. 参数选择与收敛建议

| 参数 | 推荐起点 | 收敛判据 / 说明 |
|---|---|---|
| ENCUT | 最大 ENMAX × 1.1–1.3（C/O/N 体系 450–520 eV） | $E_\mathrm{ads}$ 变化 < 0.01–0.02 eV |
| k 点（slab） | k 间距 ≈ 0.15–0.2 Å⁻¹，z 方向 1，Γ 中心 | $E_\mathrm{ads}$ 变化 < 0.02 eV；表面胞增大时按比例减少 |
| k 点（分子） | Γ 点（`vasp_gam`） | — |
| 层数 | 金属 4–5 层（底部 2 层固定），氧化物按化学计量单元 | 增加一层后 $E_\mathrm{ads}$ 变化 < 0.03 eV |
| 真空 | 吸附物顶端到上方周期像 ≥ 10 Å，总真空 15–20 Å | 配合偶极校正 |
| 表面胞（覆盖度） | p(2×2) 起，低覆盖用 p(3×3)/p(4×4) | 研究低覆盖极限时 $E_\mathrm{ads}$ 变化 < 0.05 eV |
| SIGMA | 金属 0.1–0.2（MP），分子 0.01 | 每原子 T\*S < 1 meV；比较用 energy(sigma->0) |
| EDIFF / EDIFFG | 1E-6 / −0.02（优化）；1E-7 / −0.01（频率前） | 频率计算无明显虚频 |
| 频率 | IBRION=5、POTIM=0.015、NFREE=2，只放开吸附质 | 低频 < 50 cm⁻¹ 统一处理 |
| 自旋 | O2、自由基、磁性表面 ISPIN=2 | 比较自旋/非自旋能量，确认基态 |
| vdW | 统一使用或统一不用 | 体相 a0 用同一泛函重新优化 |

收敛测试的自动化可参考 [`conv_test.sh`](../scripts/conv_test.sh) 与 [`conv_analyze.py`](../scripts/conv_analyze.py)（用法见各自 `-h`），高通量吸附计算的工作流见 [Common 第6章](../../Common/tutorials/06-workflow-automation.md)。

---

## 13. 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角表面胞用了 Monkhorst-Pack 偶数网格 | KPOINTS 改为 `Gamma` |
| `Tetrahedron method fails (number of k-points < 4)` | 分子/Γ 点计算用了 `ISMEAR = -5` | 分子用 `ISMEAR = 0`，`SIGMA = 0.01` |
| 吸附能为正且很大（> +1 eV） | 三个计算的 POTCAR/ENCUT 不一致；气相分子没优化或自旋错误；用了 TOTEN 而非 energy(sigma->0) | `grep TITEL POTCAR`、`grep ENCUT OUTCAR` 逐一对照；重算气相 |
| O2 的 magnetization 为 0 | 未设 `ISPIN = 2` 或初始 MAGMOM 太小 | `ISPIN = 2`、`MAGMOM = 2*1`，必要时 `NUPDOWN = 2` |
| top 位 CO 优化后跑到 hollow | 该泛函下 top 位不是极小点 | 如实报告；需要约束能量时固定 C 的 xy（`F F T`） |
| SCF 不收敛、能量在两个值之间振荡（开了 LDIPOL） | 电荷在 slab 两侧来回晃动 | `AMIN = 0.01`、`NELM = 200`；先不加偶极校正得到 WAVECAR 再续算 |
| 频率出现 2–3 个 50–150 cm⁻¹ 的虚频 | 结构没优化到位、EDIFF 太松、POTIM 不合适 | `EDIFFG = -0.01` 重新优化，`EDIFF = 1E-7`，试 `POTIM = 0.01/0.02` |
| `IBRION = 6` 时 Selective dynamics 被忽略、算了所有原子 | 只有 IBRION=5 支持选择性动力学 | 改用 `IBRION = 5` |
| 频率结果中出现 6 个以上近零模（吸附物） | slab 原子也被放开了 | 用 `set_freq_selective.py` 只放开吸附质 |
| VASPsol 计算能量与真空完全相同 | 使用的不是打过 VASPsol 补丁的可执行文件 | 检查编译；`grep -i solv vasp.out` 应能看到 VASPsol 相关输出 |
| VASPsol SCF 发散 | 直接从头算、ENCUT 太低、真空层过薄 | 先真空 WAVECAR，`ISTART = 1`；提高 ENCUT；真空 ≥ 15 Å |
| d 带中心与文献差 0.5 eV 以上 | 积分窗口、原子选择、NBANDS 不同 | 统一窗口（如 [−10, E_F] 或全范围）并注明 |

更多通用报错（ZBRENT、EDDDAV、BRMIX 等）见[第18章](18-postprocess-troubleshooting.md)。

---

## 练习

1. 用示例输入比较 CO 在 top 与 fcc 位的吸附能，验证 PBE 下的"CO/Pt(111) 之谜"；再打开 `IVDW = 12`（体相 a0 一并重新优化）看位点顺序是否改变。
2. 把表面胞从 p(2×2) 扩大到 p(3×3)，计算 1/9 ML 的 CO 吸附能，估计横向相互作用的大小。
3. 用 `set_freq_selective.py` 分别放开"仅 CO"和"CO + 最近的 Pt"，比较 C–O 伸缩频率与 ZPE 的变化。
4. 计算 H 在 Pt(111) 的 fcc、hcp、top 位的 ΔG_H\*（0.24 eV 近似与自算频率两种方式），讨论哪种更接近 0。
5. 在 Pt(111) 上计算 OH\*、O\*、OOH\*，用 `calc_dG.py oer` 画台阶图，检查 $\Delta G_\mathrm{OOH^*}-\Delta G_\mathrm{OH^*}$ 是否接近 3.2 eV。
6. 用 VASPsol 计算 OH\* 的溶剂化能，与 CO\* 的溶剂化能对比，解释差异。

---

## 参考资料

**官方文档**
- VASP Wiki：[LDIPOL](https://www.vasp.at/wiki/index.php/LDIPOL)、[IDIPOL](https://www.vasp.at/wiki/index.php/IDIPOL)、[DIPOL](https://www.vasp.at/wiki/index.php/DIPOL)、[IBRION](https://www.vasp.at/wiki/index.php/IBRION)、[NFREE](https://www.vasp.at/wiki/index.php/NFREE)、[Phonons from finite differences](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences)（IBRION=5 与 Selective dynamics）
- [VASP 官方 Tutorials（含表面与吸附部分）](https://www.vasp.at/tutorials/latest/)
- [VASPKIT 文档：Thermo Energy Correction（501/502）与 d-band Center（503）](https://vaspkit.com/tutorials.html)
- [ASE thermochemistry](https://docs.ase-lib.org/ase/thermochemistry/thermochemistry.html)、[pymatgen.analysis.adsorption](https://pymatgen.org/pymatgen.analysis.html)
- [VASPsol 使用说明（USAGE.md）](https://github.com/henniggroup/VASPsol/blob/master/docs/USAGE.md)

**相关 GitHub 仓库**（只链接，不转载代码）
- [henniggroup/VASPsol](https://github.com/henniggroup/VASPsol) —— 隐式溶剂补丁（LSOL、EB_K、TAU、LAMBDA_D_K），需重新编译 VASP
- [VASPsol/VASPsol](https://github.com/VASPsol/VASPsol) —— VASPsol++，原 VASPsol 的后续开发
- [ftherrien/VaspGibbs](https://github.com/ftherrien/VaspGibbs) —— 从频率计算得到 ZPE、熵与吉布斯自由能
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) —— `AdsorbateSiteFinder`、slab 生成、VASP 输入输出解析
- [QijingZheng/VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF) —— OUTCAR 振动模式转 XSF/动画
- [skelton-group/VASP-Workshop-Chester-2019](https://github.com/skelton-group/VASP-Workshop-Chester-2019) —— 含 NaCl(100) 表面吸附案例的 CC-BY 工作坊材料
- [FedeDat/VASP-tutorial](https://github.com/FedeDat/VASP-tutorial) —— 面向电化学体系的 VASP 入门课
- [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script) —— 中文社区常用脚本（含 ZPE/熵校正等），仅供参考
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**
- J. K. Nørskov et al., "Origin of the overpotential for oxygen reduction at a fuel-cell cathode", *J. Phys. Chem. B* 108, 17886 (2004). [doi:10.1021/jp047349j](https://doi.org/10.1021/jp047349j)
- J. K. Nørskov et al., "Trends in the exchange current for hydrogen evolution", *J. Electrochem. Soc.* 152, J23 (2005). [doi:10.1149/1.1856988](https://doi.org/10.1149/1.1856988)
- J. Rossmeisl, Z.-W. Qu, H. Zhu, G.-J. Kroes, J. K. Nørskov, "Electrolysis of water on oxide surfaces", *J. Electroanal. Chem.* 607, 83 (2007). [doi:10.1016/j.jelechem.2006.11.008](https://doi.org/10.1016/j.jelechem.2006.11.008)
- I. C. Man et al., "Universality in oxygen evolution electrocatalysis on oxide surfaces", *ChemCatChem* 3, 1159 (2011). [doi:10.1002/cctc.201000397](https://doi.org/10.1002/cctc.201000397)
- B. Hammer, J. K. Nørskov, "Why gold is the noblest of all the metals", *Nature* 376, 238 (1995). [doi:10.1038/376238a0](https://doi.org/10.1038/376238a0)
- P. J. Feibelman et al., "The CO/Pt(111) puzzle", *J. Phys. Chem. B* 105, 4018 (2001). [doi:10.1021/jp002302t](https://doi.org/10.1021/jp002302t)
- K. Mathew et al., "Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways", *J. Chem. Phys.* 140, 084106 (2014). [doi:10.1063/1.4865107](https://doi.org/10.1063/1.4865107)
- K. Mathew et al., "Implicit self-consistent electrolyte model in plane-wave density-functional theory", *J. Chem. Phys.* 151, 234101 (2019). [doi:10.1063/1.5132354](https://doi.org/10.1063/1.5132354)
- R. Sundararaman, W. A. Goddard III, T. A. Arias, "Grand canonical electronic density-functional theory: Algorithms and applications to electrochemistry", *J. Chem. Phys.* 146, 114104 (2017). [doi:10.1063/1.4978411](https://doi.org/10.1063/1.4978411)

---

[← 上一章](09-surface-workfunction.md) | [返回目录](../README.md) | [下一章 →](11-neb-dimer.md)
