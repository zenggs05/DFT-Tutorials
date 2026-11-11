# 第14章 介电与光学性质（Dielectric & Optical Properties with VASP: LOPTICS / LEPSILON / LCALCEPS）

本章讲解如何用 VASP 计算介电函数与光学性质（dielectric function & optical properties）：独立粒子近似下的频率相关介电函数 ε(ω)（`LOPTICS`），由 ε1、ε2 得到吸收系数 α、折射率 n、消光系数 k、反射率 R、能量损失谱和 Tauc 图；静态介电常数 ε∞ / ε0 与 Born 有效电荷（DFPT `LEPSILON`、有限电场 `LCALCEPS`，以及 `IBRION = 8` 的离子贡献）；杂化泛函 HSE06 + LOPTICS / LCALCEPS，并简要介绍 GW（`ALGO = EVGW0`）与 BSE（`ALGO = BSE`）。后处理覆盖 py4vasp、VASPKIT 711、sumo-optplot 与自写 Python 脚本。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 写出 `LOPTICS` 计算的两步 INCAR，知道 `NBANDS`、`NEDOS`、`CSHIFT`、`ISMEAR`/`SIGMA`、k 点各自如何影响光谱；
- 从 vasprun.xml 读取 density-density 与 current-current 两种介电函数，计算 α、n、k、R、EELS 并画 Tauc 图；
- 用 `LEPSILON` 得到 ε∞（含/不含局域场效应）、Born 有效电荷和压电张量，用 `IBRION = 8` 得到离子贡献和 ε0；
- 在杂化泛函下用 `LCALCEPS` 计算 ε∞ 和 Born 电荷，用 HSE06 + `LOPTICS` 改善光谱的带隙问题；
- 了解 GW、BSE 的计算链条和各自需要的关键标签，知道什么时候必须超越独立粒子近似。

| 想要的量 | 推荐方法 | 关键标签 | 本章 |
|---|---|---|---|
| ε(ω)、吸收谱、n、k、R（IPA） | Green–Kubo 求和 + Kramers–Kronig | `LOPTICS`、`NBANDS`、`CSHIFT`、`NEDOS` | §3、§4 |
| 电子静态介电常数 ε∞、Born 电荷、压电张量（LDA/GGA） | DFPT（Sternheimer 方程） | `LEPSILON`、`LRPA` | §5 |
| 同上，但用杂化泛函 | 有限电场（PEAD） | `LCALCEPS`、`EFIELD_PEAD` | §6 |
| 离子贡献 ε_ion、静态 ε0 | 有限差分 / DFPT 声子 + 场响应 | `IBRION = 5/6/7/8` + `LEPSILON`/`LCALCEPS` | §5.3 |
| 准粒子带隙、含激子的吸收谱 | GW + BSE / 模型 BSE / TDHF | `ALGO = EVGW0`、`ALGO = BSE`、`ALGO = TDHF` | §7 |

> [!NOTE]
> 本章用到 [第5章](05-band-dos.md) 的能带/DOS 基础和 [第6章](06-hybrid-metagga.md) 的杂化泛函设置；离子贡献与非解析修正部分与 [第12章 声子](12-phonon.md) 衔接。

---

## 2. 原理速览

### 2.1 复介电函数与光学常数

线性响应下，介电函数是复数张量 ε_αβ(ω) = ε1(ω) + iε2(ω)。ε2 描述吸收，ε1 描述色散与极化，两者由 Kramers–Kronig（KK）关系联系。由 ε 可以得到所有常用光学常数：

| 量 | 公式 | 备注 |
|---|---|---|
| 复折射率 | ñ = n + ik = √ε | n = √[(\|ε\| + ε1)/2]，k = √[(\|ε\| − ε1)/2] |
| 吸收系数 | α(ω) = 2ωk/c = 2Ek/(ħc) | ħc = 1.9733×10⁻⁵ eV·cm，α 单位 cm⁻¹ |
| 正入射反射率 | R = [(n−1)² + k²] / [(n+1)² + k²] | — |
| 能量损失函数（EELS） | L(ω) = −Im[1/ε] = ε2 / (ε1² + ε2²) | 峰位对应等离激元 |
| 光电导（实部，SI） | σ1(ω) = ε0 ω ε2(ω) | — |

各向异性晶体要分别看 xx、yy、zz 分量；"各向同性平均"通常取 (ε_xx + ε_yy + ε_zz)/3 后再计算 n、k、α。

### 2.2 独立粒子近似（`LOPTICS`）

`LOPTICS = .TRUE.` 在基态收敛后，对所有占据态 v 与空态 c 的带间跃迁求和得到 ε2：

$$
\varepsilon^{(2)}_{\alpha\beta}(\omega) = \frac{4\pi^2 e^2}{\Omega}\lim_{q\to 0}\frac{1}{q^2}
\sum_{c,v,\mathbf{k}} 2w_{\mathbf{k}}\,\delta(\epsilon_{c\mathbf{k}}-\epsilon_{v\mathbf{k}}-\omega)\,
\langle u_{c\mathbf{k}+\mathbf{e}_\alpha q}|u_{v\mathbf{k}}\rangle
\langle u_{v\mathbf{k}}|u_{c\mathbf{k}+\mathbf{e}_\beta q}\rangle
$$

再用带复数位移 η 的 KK 变换得到 ε1：

$$
\varepsilon^{(1)}_{\alpha\beta}(\omega) = 1 + \frac{2}{\pi}\,\mathcal{P}\!\int_0^\infty
\frac{\varepsilon^{(2)}_{\alpha\beta}(\omega')\,\omega'}{\omega'^2-\omega^2+i\eta}\,d\omega'
$$

η 就是 `CSHIFT`（洛伦兹展宽，eV）。要点：

- 求和需要**大量空带**：官方建议 `NBANDS` 约为默认值的 2–3 倍；
- 频率网格点数由 `NEDOS` 控制，官方强烈建议约 2000；
- 这一近似**不含局域场效应**，也**不含激子效应**（电子-空穴相互作用）；
- 对杂化泛函同样适用（对 k 的导数改用有限差分），所以 HSE06 + LOPTICS 是改善带隙误差的常用办法；
- VASP ≥ 6.6.0 支持 LOPTICS 与依赖动能密度的 meta-GGA（如 r2SCAN）组合。

### 2.3 局域场效应、RPA 与 "in DFT"

宏观介电函数要对微观介电矩阵求逆后取 G = G' = 0 分量，包含了晶胞内电势的微观变化（局域场效应，local field effects）。VASP 的静态方法会同时给出几组结果：

| OUTCAR 中的描述 | 含义 |
|---|---|
| excluding local field effects | 不含局域场，应与 `LOPTICS` 的 ε1(ω→0) 一致 |
| including local field effects in RPA (Hartree) | 只包含 Hartree 势的局域变化 |
| including local field effects in DFT | Hartree + 交换关联势的局域变化（`LRPA = .FALSE.`，默认） |

官方 `LRPA` 说明指出，加入交换关联项通常会使介电常数增大；Born 电荷和压电张量只在 `LRPA = .FALSE.` 时有意义。

### 2.4 ε∞、ε0 与离子贡献

- **ε∞（电子 / 钳位离子介电常数）**：频率远高于声子、远低于带隙时的介电常数，只有电子响应。
- **ε0（静态介电常数）**：再加上离子位移的贡献，ε0 = ε∞ + ε_ion。离子贡献由 Born 有效电荷 Z* 与 Γ 点力常数决定：

$$
\varepsilon^{\mathrm{ion}}_{\alpha\beta} = \frac{4\pi}{\Omega}\sum_{\nu}\frac{S_{\alpha\beta,\nu}}{\omega_\nu^2},\qquad
S_{\alpha\beta,\nu}=\Big(\sum_{I\gamma} Z^*_{I\alpha\gamma}\,e_{I\gamma,\nu}\Big)\Big(\sum_{J\delta} Z^*_{J\beta\delta}\,e_{J\delta,\nu}\Big)
$$

（e 为本征位移，ω_ν 为 Γ 点 TO 频率；完整的频率相关形式见官方"Dielectric properties"页面与 Gonze & Lee 1997。）对立方双原子晶体，Lyddane–Sachs–Teller（LST）关系 ε0/ε∞ = (ω_LO/ω_TO)² 是很好的自洽检查。

### 2.5 density-density 与 current-current

电磁场的规范选择不同，会得到两种响应函数：密度-密度（χ_ρρ，标量势耦合电荷密度）与电流-电流（χ_jj，矢势耦合电流）。在有限 q 与 ω 下两者给出相同的介电函数；VASP 6 的 vasprun.xml 中通常同时写出两组：

```xml
<dielectricfunction comment="density-density"> ... </dielectricfunction>
<dielectricfunction comment="current-current"> ... </dielectricfunction>
```

官方说明：current-current 在 q→0 与 q = 0 都严格成立，对金属能给出 ω=0 附近正确的 Drude 尾，但数值上更不稳定。半导体/绝缘体的吸收谱通常用 density-density；比较不同程序或与实验对比时要明确用的是哪一种（Sangalli et al. 2017 讨论了这一常见陷阱）。

---

## 3. 实战一：Si 的频率相关介电函数（LOPTICS 两步法）

完整输入在 [../examples/14-optics-si-loptics/README.md](../examples/14-optics-si-loptics/README.md)。

### 3.1 第 1 步：自洽基态

```bash
# INCAR.scf
SYSTEM  = Si LOPTICS step 1 (SCF)
PREC    = Accurate
ENCUT   = 350
EDIFF   = 1E-8
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
LWAVE   = .TRUE.
LCHARG  = .TRUE.
```

KPOINTS 用 Γ 中心 16×16×16（光谱比总能更依赖 k 点密度），POSCAR 为 Si 原胞（a = 5.47 Å）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
cp INCAR.scf INCAR
mpirun -np 16 vasp_std > scf.out
```

### 3.2 第 2 步：大量空带 + LOPTICS

```bash
# INCAR.optics
SYSTEM  = Si LOPTICS step 2 (optics)
ISTART  = 1            # read WAVECAR from step 1
ALGO    = Exact        # exact diagonalisation: all NBANDS are converged in one step
NELM    = 1
NBANDS  = 32           # default would be 8 for Si2; LOPTICS needs 2-3x (test 24/32/48)

PREC    = Accurate     # keep identical to step 1
ENCUT   = 350
EDIFF   = 1E-8
ISMEAR  = 0
SIGMA   = 0.01         # Gaussian width << CSHIFT, so only the Lorentzian remains
LREAL   = .FALSE.

LOPTICS = .TRUE.
CSHIFT  = 0.100        # complex shift (Lorentzian broadening, eV) in Kramers-Kronig
NEDOS   = 2000         # frequency grid points

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `LOPTICS` | `.TRUE.` | 计算频率相关介电矩阵（IPA），写入 vasprun.xml 与 OUTCAR |
| `NBANDS` | 默认值的 2–3 倍起，做收敛测试 | 空带不够会使 ε1(0) 偏小、高能部分谱形错误 |
| `ALGO` + `NELM` | `Exact` + 1（读入 WAVECAR） | 与官方 GW 教程中准备空带的做法相同；小体系也可直接一步 `ALGO = Normal` |
| `CSHIFT` | 0.1（默认） | KK 变换的洛伦兹展宽；带隙很小（接近 2×CSHIFT）时应减小，并同时增大 NEDOS |
| `NEDOS` | 2000 左右 | 频率网格点数（同时也是 DOS 网格点数） |
| `ISMEAR`/`SIGMA` | 0 / 0.01–0.05 | 高斯展宽直接作用于 ε2；官方建议 SIGMA 远小于 CSHIFT，避免两种展宽叠加。`ISMEAR = -5`（四面体）与 LOPTICS 组合需 VASP ≥ 6.3；`ISMEAR = -2` 不支持 |
| k 点 | 远密于总能收敛值 | 光谱是 k 空间积分，粗网格会产生虚假的尖峰 |
| `LNABLA` | 不设 | 横向表达式，官方说明精度更低，仅用于测试 |

```bash
cp INCAR.optics INCAR
mpirun -np 16 vasp_std > optics.out
```

> [!WARNING]
> 不要把 `LOPTICS` 与 `LEPSILON` 放在同一次计算中（官方不推荐，即使某些版本能运行）。`LEPSILON` 不需要空带，增加 `NBANDS` 反而让它变慢。

### 3.3 输出检查

```bash
grep "NBANDS=" OUTCAR                        # 实际能带数（并行时可能被自动调大）
grep -n "DIELECTRIC FUNCTION" OUTCAR | head  # IMAGINARY / REAL，density-density / current-current
grep "<dielectricfunction" vasprun.xml       # vasprun.xml 中的数据块
```

OUTCAR 中的数据块标题形如 `frequency dependent IMAGINARY DIELECTRIC FUNCTION (independent particle, no local field effects) density-density`，下面每行依次是能量（eV）和 xx、yy、zz、xy、yz、zx 六个分量。vasprun.xml 中 `<imag>` 与 `<real>` 两个数组的每一行 `<r>` 也是"能量 + 6 个分量"。

成功的标志：

- ε2 在带隙以下（ω < E_g）接近 0，在带隙以上出现吸收峰；
- ε1(ω→0) 为正且与 `LEPSILON` 中"excluding local field effects"的值一致（§3.4）；
- 立方晶体 xx = yy = zz，非对角元为 0。

用示例脚本画图并导出数据：

```bash
python3 plot_optics.py                        # 各向同性平均，density-density
python3 plot_optics.py --kind current         # current-current
python3 plot_optics.py --component zz         # 单个分量（各向异性晶体）
```

脚本输出（示意输出，数值仅供参考）：

```text
read 2000 frequency points (density-density), E = 0.000 ... 3x.xxx eV
eps1(omega->0) = 1x.xxx   (compare with LEPSILON, no local fields)
n(omega->0)    = 3.xxx
wrote optics.dat
wrote optics.png
```

量级参考（近似值）：Si 的实验 ε∞ 约 11.9，ε2 主峰（E1、E2 跃迁）约在 3.4 eV 与 4.2 eV。PBE 带隙偏小，因此 ε∞ 偏大、吸收峰整体红移约 0.5 eV——这是泛函误差，不是输入错误。

### 3.4 收敛测试：NBANDS、k 点与 CSHIFT

最可靠的 NBANDS 判据来自官方建议：先做一次 `LEPSILON = .TRUE.`（DFPT，不需要空带），记下"excluding local field effects"的 ε∞；再逐步增加 LOPTICS 的 `NBANDS`，直到 ε1(ω→0) 与之一致。

| 测试 | 做法 | 收敛判据 |
|---|---|---|
| NBANDS | 24 → 32 → 48 → 64 | ε1(0) 变化 < 1%，高能谱形不再变 |
| k 点 | 12³ → 16³ → 20³ → 24³ | 主峰位置与高度稳定，没有锯齿状尖峰 |
| CSHIFT | 0.2 → 0.1 → 0.05（同时增大 NEDOS） | 只影响展宽；报告结果时注明所用展宽 |

```bash
# 例：扫描 NBANDS（每个目录都需要第 1 步的 WAVECAR）
for nb in 24 32 48 64; do
  mkdir -p nb_$nb && cp INCAR.optics POSCAR KPOINTS POTCAR WAVECAR nb_$nb/
  sed -i "s/^NBANDS.*/NBANDS  = $nb/" nb_$nb/INCAR.optics && mv nb_$nb/INCAR.optics nb_$nb/INCAR
  (cd nb_$nb && mpirun -np 16 vasp_std > out && python3 ../plot_optics.py | grep "eps1")
done
```

### 3.5 HSE06 + LOPTICS

PBE 低估带隙会使吸收边红移、ε∞ 偏大。杂化泛函可以直接与 LOPTICS 组合：

1. HSE06 自洽（参考 [第6章](06-hybrid-metagga.md)），保存 WAVECAR；
2. 读入 WAVECAR，增加空带，打开 LOPTICS：

```bash
# INCAR (HSE06 + LOPTICS, step 2) -- fragment
ISTART   = 1
LHFCALC  = .TRUE.
HFSCREEN = 0.2
ALGO     = Normal        # or All / Damped if the SCF struggles
NBANDS   = 32            # increase as in the PBE case
LOPTICS  = .TRUE.
CSHIFT   = 0.1
NEDOS    = 2000
ISMEAR   = 0
SIGMA    = 0.01
```

杂化泛函计算量大，k 点往往无法加到 PBE 那么密；常见的折中是用 PBE 的光谱形状 + HSE 带隙做比较，或用 Wannier 插值（见 [第17章](17-advanced-electronic.md)）。

> [!WARNING]
> 非共线（SOC）+ 依赖动能密度的 meta-GGA + LOPTICS 在 VASP 6.6 中存在已知问题（known issue 100：光学矩阵元的动能密度势项作用在错误的旋量分量上，介电函数偏差约百分之几）。使用这一组合前请查看官方 Known issues 页面的最新状态与规避方法。

### 3.6 金属：带内（Drude）贡献

`LOPTICS` 的 density-density 结果只包含带间跃迁。金属的低频光学响应由带内跃迁（自由载流子）主导，需要额外加上 Drude 项：

$$
\varepsilon^{\mathrm{intra}}(\omega) = 1 - \frac{\omega_p^2}{\omega(\omega + i\gamma)}
$$

VASP 在 OUTCAR 中给出由带内跃迁计算的等离子体频率平方（`grep -i "plasma frequency squared" OUTCAR`），γ 需要取自实验或经验值。current-current 响应在 ω→0 附近包含 Drude 尾，但数值上更敏感。金属需要非常密的 k 点。

---

## 4. 吸收系数、折射率、反射率、EELS 与 Tauc 图

### 4.1 由 ε1、ε2 计算光学常数

示例脚本 `plot_optics.py`（只依赖 numpy、matplotlib）读取 vasprun.xml，输出 `optics.dat`，各列为：

| 列 | 量 | 单位 |
|---|---|---|
| 1 | 光子能量 E | eV |
| 2, 3 | ε1, ε2 | — |
| 4, 5 | n, k | — |
| 6 | α | cm⁻¹ |
| 7 | R | — |
| 8 | 能量损失函数 L | — |

核心计算只有几行：

```python
import numpy as np
HBAR_C = 1.973269804e-5                      # eV*cm
mod = np.sqrt(e1**2 + e2**2)
n = np.sqrt((mod + e1) / 2)
k = np.sqrt((mod - e1) / 2)
alpha = 2 * E * k / HBAR_C                   # cm^-1
R = ((n - 1)**2 + k**2) / ((n + 1)**2 + k**2)
loss = e2 / (e1**2 + e2**2)
```

波长换算：λ(nm) = 1239.84 / E(eV)。半导体吸收边附近 α 通常在 10⁴–10⁶ cm⁻¹ 量级，适合用对数坐标作图。

### 4.2 Tauc 图

Tauc 关系 (αhν)^{1/r} = A(hν − E_g) 常用于从吸收谱外推光学带隙，r = 1/2 对应直接允许跃迁，r = 2 对应间接允许跃迁：

```bash
python3 plot_optics.py --tauc direct --fit 3.0 3.6     # 在 3.0-3.6 eV 的线性段拟合并外推
```

> [!CAUTION]
> 独立粒子近似只包含竖直（直接）跃迁，**不包含声子辅助的间接跃迁**，也不包含激子。因此：
> - 对间接带隙材料（如 Si），用 IPA 光谱做"间接 Tauc 图"没有物理意义，吸收边对应的是直接带隙；
> - 外推得到的"光学带隙"继承了所用泛函的带隙误差，并受 CSHIFT/SIGMA 展宽影响，拟合窗口的选择也带有主观性。
> 报告 Tauc 带隙时请同时给出泛函、展宽和拟合区间，并与电子带隙（`BANDGAP` 标签或能带图）对照。

### 4.3 二维材料与含真空的体系

对含真空层的 slab 或单层材料，ε 按整个超胞体积归一化，因而 α、n、R 等"体材料量"随真空厚度变化，不能直接与实验比较。常用的处理：

- 报告与真空无关的量，例如薄层极限下的吸收率（absorbance）A(ω) ≈ (ω L / c) ε2(ω)，其中 L 为超胞在垂直方向的长度（正入射、弱吸收近似）；或二维极化率 α_2D(ω) = L[ε(ω) − 1]/(4π)；
- 区分面内（xx、yy）与面外（zz）偏振；
- 官方 VASPKIT 文档说明其 711 功能不适用于低维材料，二维体系请自行处理。

---

## 5. 实战二：GaAs 的静态介电常数与 Born 有效电荷（LEPSILON）

完整输入在 [../examples/14-dielectric-gaas-lepsilon/README.md](../examples/14-dielectric-gaas-lepsilon/README.md)。

### 5.1 INCAR

```bash
# INCAR (Part A)
SYSTEM   = GaAs LEPSILON (DFPT, electronic part)
PREC     = Accurate
ENCUT    = 400         # Ga_d/As: ~1.4 x max(ENMAX)
EDIFF    = 1E-8
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.

LEPSILON = .TRUE.      # linear response to an electric field (Sternheimer equations)
LRPA     = .FALSE.     # default: local fields at Hartree + xc level; Born charges printed

IBRION   = -1
NSW      = 0
LWAVE    = .FALSE.
LCHARG   = .FALSE.

# Parallelisation: do NOT set NCORE or NPAR for LEPSILON runs (keep NCORE = 1).
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `LEPSILON` | `.TRUE.` | DFPT 计算 ε∞（含/不含局域场）、Born 有效电荷、钳位离子压电张量 |
| `LRPA` | `.FALSE.`（默认） | `.TRUE.` 只含 Hartree 局域场；Born 电荷与压电张量只在 `.FALSE.` 时有意义 |
| `NBANDS` | 默认 | 只需占据态；增大反而变慢 |
| `NCORE` / `NPAR` | 不设置 | 设置后会报错退出；可用 `KPAR` |
| `EDIFF` | 1E-8 | 线性响应迭代需要精确的基态 |
| k 点 | Γ 中心，比总能更密 | Born 电荷与 ε∞ 收敛较慢，小带隙材料尤甚 |
| 泛函 | LDA / GGA | 不支持杂化泛函（改用 §6 的 `LCALCEPS`）；`IBRION = 7/8` 也不支持 `IVDW` |

POSCAR 采用实验晶格常数 a = 5.653 Å：PBE 在其自身平衡晶格下 GaAs 带隙接近 0，线性响应会变得很难收敛。KPOINTS 为 Γ 中心 12×12×12。

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR
mpirun -np 16 vasp_std > lepsilon.out
```

运行时 OSZICAR / stdout 中会出现以 `RMM:` 开头、带 `d eps` 列的迭代信息——这是在逐个方向求解线性 Sternheimer 方程，属于正常现象。

### 5.2 读取输出

```bash
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR" OUTCAR   # 各种近似下的 eps_inf
grep -A 10 "BORN EFFECTIVE CHARGES" OUTCAR                # 每个离子一个 3x3 张量
grep -A 6 "PIEZOELECTRIC TENSOR" OUTCAR                   # 单位 C/m^2
```

介电张量块形如（示意输出，数值仅供参考）：

```text
 MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
 ------------------------------------------------------
          1x.xxxxxx     0.000000     0.000000
           0.000000    1x.xxxxxx     0.000000
           0.000000     0.000000    1x.xxxxxx
 ------------------------------------------------------
```

判断标准：

- 立方对称：张量对角、三个对角元相等；
- Born 电荷满足电荷中性求和规则：所有原子的 Z* 之和接近零（Ga 为正、As 为负，大小相等）。偏差大说明 k 点不够；
- 实验参考（近似值）：ε∞ ≈ 10.9、ε0 ≈ 12.9、Z* ≈ 2.1。PBE 低估带隙 → ε∞ 往往明显偏大。

> [!NOTE]
> 已知问题：带电晶胞中 VASP 会把 Born 电荷之和强制移到零（known issue 92），这只对中性晶胞成立；非共线计算在 VASP < 6.5.0 中 `LEPSILON` + 对称性结果错误（known issue 20，旧版本用 `ISYM = -1`）；GPU 上的非共线线性响应也有问题（known issue 114）。

### 5.3 离子贡献：`IBRION = 8` + `LEPSILON`

在**同一套参数、已充分弛豫**的结构上，把 `IBRION` 改为 8（DFPT 声子，利用对称性），VASP 会同时计算 Γ 点力常数和 Born 电荷，进而给出离子对介电张量的贡献：

```bash
# INCAR.ionic (Part B) -- differences to Part A
LEPSILON = .TRUE.
IBRION   = 8
NSW      = 1          # only needed by very old versions; harmless
```

```bash
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION" OUTCAR
grep "THz" OUTCAR      # Gamma frequencies: 3 acoustic (~0) + 3 TO modes for GaAs
```

静态介电常数 ε0 = ε∞（Part A，含局域场）+ ε_ion（Part B）。用 LST 关系检查：ε0/ε∞ 应等于 (ω_LO/ω_TO)²，其中 ω_TO 来自 `grep THz`，ω_LO 可由非解析修正得到或与实验比较。

对需要 vdW 修正（`IVDW`）或杂化泛函的体系，DFPT 声子不可用，改用 `IBRION = 6` + `LEPSILON`（GGA）或 `IBRION = 6` + `LCALCEPS`（杂化）。

### 5.4 与声子非解析修正的衔接

极性材料的声子谱在 Γ 点有 LO–TO 劈裂，需要 Born 电荷和 ε∞ 做非解析修正（NAC）。phonopy 的 `phonopy-vasp-born` 可以从本节 `LEPSILON` 计算的 vasprun.xml 生成 BORN 文件，具体流程见 [第12章](12-phonon.md)。

---

## 6. 实战三：ZnO 的 HSE06 + LCALCEPS（有限电场法）

完整输入在 [../examples/14-dielectric-zno-hse-lcalceps/README.md](../examples/14-dielectric-zno-hse-lcalceps/README.md)。

### 6.1 为什么需要 LCALCEPS

DFPT（`LEPSILON`）目前只支持不显式依赖轨道的泛函，杂化泛函不能用。`LCALCEPS = .TRUE.` 改为对体系施加有限的均匀电场（基于 PEAD 的电焓泛函），由极化 P、受力 F、应力 σ 的变化做有限差分：

$$
\varepsilon^\infty_{ij}=\delta_{ij}+\frac{4\pi}{\epsilon_0}\frac{\partial P_i}{\partial \mathcal{E}_j},\qquad
Z^*_{ij}=\frac{1}{e}\frac{\partial F_j}{\partial \mathcal{E}_i}
$$

VASP 会依次做零场和沿 x、y、z 加场（默认 0.01 eV/Å）的自洽计算。

### 6.2 INCAR

```bash
# INCAR.hse_lcalceps -- start from a PBE WAVECAR (INCAR.pbe)
SYSTEM   = ZnO HSE06 LCALCEPS
ISTART   = 1
PREC     = Accurate
ENCUT    = 520
EDIFF    = 1E-7        # the field-polarised states are converged to ~EDIFF/100
ISMEAR   = 0
SIGMA    = 0.01
LREAL    = .FALSE.
NBANDS   = 24          # a few empty bands so that VASP can check the field strength vs. gap

# HSE06
LHFCALC  = .TRUE.
HFSCREEN = 0.2
AEXX     = 0.25
ALGO     = Normal      # if the SCF struggles: ALGO = All
PRECFOCK = Normal

# finite electric field
LCALCEPS = .TRUE.
# EFIELD_PEAD = 0.01 0.01 0.01   # default field strength (eV/A); reduce if VASP warns

IBRION   = -1
NSW      = 0
LWAVE    = .FALSE.
LCHARG   = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `LCALCEPS` | `.TRUE.` | 有限电场法计算 ε∞、Z*、钳位离子压电张量 |
| `EFIELD_PEAD` | 默认 3×0.01 eV/Å | 电场太大时电焓泛函失去极小值，VASP 会在场强超过带隙相关的安全阈值时给出警告；与 `LCALCEPS` 联用时小于 1E-5 的值会被重置为默认值 |
| KPOINTS | Γ 中心规则网格 | 官方说明 PEAD 例程只适用于包含 Γ 点的规则网格（Γ 中心或奇数网格） |
| `EDIFF` | 1E-6 ~ 1E-7 | 加场后的自洽默认收敛到 EDIFF/100，代价较高 |
| `NBANDS` | 占据带 + 少量空带 | 没有空带时 VASP 无法估计带隙、无法检查电场大小，会给出警告 |
| `IBRION = 6` | 可选 | 与 `LCALCEPS` 组合可再得到离子对介电与压电张量的贡献 |

### 6.3 输出与检查

```bash
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)" OUTCAR
grep -A 14 "BORN EFFECTIVE CHARGES (including local field effects)" OUTCAR
grep -A 8 "PIEZOELECTRIC TENSOR (including local field effects)" OUTCAR
```

- 纤锌矿是单轴晶体：εxx = εyy ≠ εzz；
- 实验 ε∞ 约 3.7（近似值）；PBE 带隙约 0.7–0.8 eV（实验约 3.4 eV），ε∞ 明显偏大；HSE06 带隙改善到 2.5 eV 左右（仍偏低），ε∞ 更接近实验（数值仅供参考）。
- 同一体系用 PBE 的 `LEPSILON` 与 `LCALCEPS` 应给出一致的结果，这是检验 LCALCEPS 设置（k 点、场强、EDIFF）的好办法。

> [!WARNING]
> **版本与平台**：GPU（OpenACC）版本在 6.6.1 之前用 `LCALCEPS` 得到的 Born 有效电荷有误，介电张量不受影响（known issue 93）；VASP 5.4.4 及更早版本的压电张量符号错误（需打官方补丁）。用 GPU 做 LCALCEPS 时请使用 ≥ 6.6.1，或用 CPU 复核 Born 电荷。

---

## 7. 超越独立粒子近似：GW 与 BSE 简介

独立粒子近似有两个根本缺陷：带隙来自 Kohn–Sham 本征值（GGA 明显偏小），且没有电子-空穴相互作用（激子）。对宽带隙绝缘体、二维材料和分子晶体，激子效应可以使吸收峰红移零点几到 1 eV 以上，并改变谱形。VASP 的多体微扰理论模块提供了完整的解决方案，这里只给出计算链条和关键标签，具体请跟随官方教程。

### 7.1 GW 准粒子能量

VASP ≥ 6.3 可以一步完成 GW（内部先做精确对角化得到所有可用空带），此时 INCAR 中**不能**设置 `NBANDS`：

```bash
# one-step GW (VASP >= 6.3) -- fragment
ALGO    = EVGW0      # G0W0 / partially self-consistent GW0 (VASP 5.x: ALGO = GW0)
NELMGW  = 1          # number of GW iterations (VASP <= 6.2: use NELM)
LOPTICS = .TRUE.     # insulators: WAVEDER for the head/wings of the dielectric matrix
ISMEAR  = 0
SIGMA   = 0.05
# do NOT set NBANDS for the one-step procedure
```

传统的两步法：先做 DFT（或 HSE）基态，再读 WAVECAR 用 `ALGO = Exact`、`NELM = 1`、很大的 `NBANDS` 和 `LOPTICS = .TRUE.` 生成空带与 WAVEDER；然后 `ALGO = EVGW0`、相同 `NBANDS`、`NOMEGA`（频率点数）做 GW。官方提醒：

- 必须仔细收敛 `NBANDS` 与 `NOMEGA`（以及 `ENCUTGW`）；
- 金属一般应去掉 `LOPTICS` 并删除 WAVEDER；
- 准粒子能量写在 OUTCAR 中（`grep -A 20 "QP shifts" OUTCAR` 附近的表格），可配合 Wannier90 插值得到 GW 能带。

### 7.2 BSE 吸收谱

BSE 需要前一步 GW 计算写出的屏蔽库仑核文件（`Wxxxx.tmp` 只含对角项，`WFULLxxxx.tmp` 为完整矩阵；分子、表面等必须用完整矩阵），以及 WAVECAR 与 WAVEDER：

```bash
# BSE step -- fragment (after a GW run with LWAVE = .TRUE.)
ALGO     = BSE
NBANDS   = same as in the GW run
NBANDSO  = 4          # occupied bands included in the BSE Hamiltonian
NBANDSV  = 8          # unoccupied bands included
OMEGAMAX = 20         # maximum electron-hole pair energy (eV) kept in the Hamiltonian
ISMEAR   = 0
SIGMA    = 0.05
# ANTIRES = 0         # Tamm-Dancoff approximation; 2 = beyond TDA (use LORBITALREAL = .TRUE.)
```

- `NBANDSO`/`NBANDSV` 只需包含带边附近的少数能带，但内存随能带数迅速增加；
- 收敛 `NBANDSV` 与 `OMEGAMAX` 时要特别注意 ε1（实部）——官方指出实部比虚部敏感得多；
- BSE 的介电函数与激子能量写入 vasprun.xml；VASP 6.5 起 BSE 可在 GPU 上运行，6.6 起可用 BSE 计算 X 射线吸收谱（XAS）；
- 官方建议的一致性检查：关闭 `LADDER` 与 `LHARTREE` 时 BSE 应重现 LOPTICS 的 IPA 谱（两者设同样的 `CSHIFT`）。

### 7.3 更便宜的替代：模型 BSE、TDHF 与时间演化

- **模型 BSE（mBSE）**：用模型介电函数代替 GW 屏蔽（`LMODELHF = .TRUE.` + `AEXX`、`HFSCREEN` 参数，`ALGO = TDHF`、`IBSE = 0`），无需 GW 步骤；
- **TDHF / TDDFT**（`ALGO = TDHF`）：用杂化泛函核求解 Casida 方程；
- **时间演化**（`ALGO = TIMEEV`）：施加 δ 脉冲电场跟踪偶极矩，可在合适的核设置下重现 BSE 谱（VASP < 6.6.0 时须 `ISYM = 0`，见 known issue 44）。

| 方法 | 带隙 | 激子 | 相对代价 |
|---|---|---|---|
| PBE + LOPTICS | 偏小 | 无 | 1 |
| HSE06 + LOPTICS | 改善 | 无 | 10–100 |
| 模型 BSE / TDHF | 取决于参数 | 有 | 10–100 |
| G0W0 + BSE | 准粒子水平 | 有 | 100–1000 |

---

## 8. 后处理工具

| 工具 | 用法 | 说明 |
|---|---|---|
| 示例脚本 `plot_optics.py` | `python3 plot_optics.py --tauc direct` | 读 vasprun.xml，输出 optics.dat 与图；区分 density / current |
| py4vasp（vaspout.h5） | 见下 | 官方 Python 库，需 HDF5 编译 |
| VASPKIT 711 | `vaspkit -task 711` | 读 vasprun.xml（或 REAL.in/IMAG.in），可选能量单位 eV / nm / THz，输出吸收、折射率、消光、反射、损失谱等数据文件；官方说明不适用于低维材料 |
| sumo-optplot | `sumo-optplot`、`sumo-optplot absorption loss --anisotropic` | 出版级图；可选 absorption、loss、eps_real、eps_imag、n_real、n_imag，`--bandgap` 标出带隙，`--gaussian` 额外展宽 |
| pymatgen | `Vasprun("vasprun.xml").dielectric`、`.optical_absorption_coeff` | 程序化读取能量、实部、虚部和吸收系数（cm⁻¹） |
| 仓库脚本 | ../scripts/optics_post.py | 介电函数、光学常数、Tauc 图以及 LEPSILON/LCALCEPS 张量与 Born 电荷的读取；用 `python ../scripts/optics_post.py -h` 查看用法 |

py4vasp 的基本用法：

```python
from py4vasp import Calculation
calc = Calculation.from_path(".")
calc.dielectric_function.plot()                 # 默认画各向同性平均的实部与虚部
print(calc.dielectric_function.selections())    # 查看可选的分量、方向、实部/虚部等关键字
data = calc.dielectric_function.read()          # 字典：能量 + 3x3 张量
calc.dielectric_tensor.print()                  # LEPSILON / LCALCEPS 的静态介电张量
calc.born_effective_charge.print()              # Born 有效电荷
```

---

## 参数选择与收敛建议

| 参数 | 推荐取值 | 说明 |
|---|---|---|
| `NBANDS`（LOPTICS） | 默认值的 2–3 倍起，用 LEPSILON 的 ε∞ 标定 | 不够时 ε1(0) 偏小 |
| `NEDOS` | ≈ 2000 | 减小 CSHIFT 时同步增大 |
| `CSHIFT` | 0.1（默认），小带隙时减小 | 报告光谱时注明展宽 |
| `SIGMA` | 远小于 CSHIFT（0.01–0.05） | 避免高斯与洛伦兹展宽叠加 |
| k 点（光谱） | 明显密于总能；Γ 中心 | 粗网格产生尖峰；用 12³→24³ 测试（Si 原胞） |
| `ISMEAR` | 0；`-5` 需 VASP ≥ 6.3；不能用 `-2` | — |
| `LEPSILON` 并行 | 不设 NCORE/NPAR；可用 KPAR | — |
| `LCALCEPS` k 网格 | 包含 Γ 的规则网格 | PEAD 例程要求 |
| `EFIELD_PEAD` | 默认 0.01 eV/Å；有警告时减半 | 场强过大失去极小值 |
| `EDIFF` | 1E-7 ~ 1E-8 | 线性响应与有限场都需要 |
| 泛函 | PBE 做趋势，HSE06/GW 做定量 | PBE 低估带隙 → 光谱红移、ε∞ 偏大 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决办法 |
|---|---|---|
| LEPSILON 计算报错并提示 PEAD 例程不支持 `NCORE`，或要求去掉 `NPAR` | 线性响应 / PEAD 不能做能带并行 | 删除 `NCORE`、`NPAR`（默认 NCORE = 1），改用 `KPAR` |
| 杂化泛函 + `LEPSILON` 无法运行或结果无意义 | DFPT 不支持杂化泛函 | 改用 `LCALCEPS = .TRUE.` |
| `WARNING: The number of bands has been changed from X to Y` | 并行设置要求 NBANDS 为某数的整数倍 | 光学计算中多几条带无害；需要固定能带数时按并行组数取整 |
| ε2 谱布满尖锐的锯齿 | k 点太稀、CSHIFT/SIGMA 太小 | 加密 k 点，CSHIFT 回到 0.1 |
| ε1(0) 明显小于 LEPSILON 的 ε∞（不含局域场） | 空带不足 | 增大 NBANDS 直到一致 |
| VASP 警告 CSHIFT 太小 | 洛伦兹展宽小于频率网格可分辨的宽度 | 增大 CSHIFT 或 NEDOS |
| `LOPTICS` + `ISMEAR = -2` 失败 | 官方不支持该组合 | 改用 `ISMEAR = 0` |
| `ISMEAR = -5` + LOPTICS 出错（旧版本） | 四面体法与 LOPTICS 组合需 ≥ 6.3 | 升级或改用 `ISMEAR = 0` |
| LCALCEPS 输出电场过大的警告，或加场自洽不收敛 | 场强相对 E_g/N_k 太大，电焓泛函无极小值 | 减小 `EFIELD_PEAD`；加少量空带让 VASP 能估计带隙 |
| Born 电荷之和明显不为零 | k 点不够、EDIFF 太松；带电晶胞的已知问题 | 加密 k 点、收紧 EDIFF；带电体系注意 known issue 92 |
| 金属的低频介电函数与实验完全不符 | 缺少 Drude 带内项 | 用 OUTCAR 的等离子体频率加 Drude 项，或看 current-current |
| 二维材料的 α、n 随真空层变化 | 体材料量按超胞体积归一化 | 改报吸收率 A(ω) 或二维极化率（§4.3） |
| `LOPTICS` 与 `LEPSILON` 同时开启后结果异常 | 官方不建议在同一次运行中组合 | 分两次计算 |

---

## 练习

1. 对 Si 做 NBANDS = 16、24、32、48 的 LOPTICS 计算，画出 ε1(0) 随 NBANDS 的变化，并与 LEPSILON 的"excluding local field effects"值比较。
2. 比较 density-density 与 current-current 两组 Si 介电函数，在哪个能量范围差别最大？
3. 对 GaAs 计算 ε∞、ε0、Z*，用 LST 关系由 ω_TO 估算 ω_LO，并与实验的 LO 频率比较。
4. 对 ZnO 分别用 PBE（LEPSILON）、PBE（LCALCEPS）、HSE06（LCALCEPS）计算 ε∞，讨论带隙误差对 ε∞ 的影响。
5. 计算单层 MoS2 的 LOPTICS，改变真空层厚度（15、20、25 Å），验证 ε2 随真空变化而吸收率 A(ω) 基本不变。

---

## 参考资料

**官方文档**

- [VASP Wiki: LOPTICS](https://www.vasp.at/wiki/index.php/LOPTICS) — IPA 介电函数公式、NBANDS/NEDOS 建议、展宽说明与 ISMEAR 限制
- [VASP Wiki: CSHIFT](https://www.vasp.at/wiki/index.php/CSHIFT)、[NEDOS](https://www.vasp.at/wiki/index.php/NEDOS)、[LNABLA](https://www.vasp.at/wiki/index.php/LNABLA)
- [VASP Wiki: LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON)、[LRPA](https://www.vasp.at/wiki/index.php/LRPA) — DFPT 静态介电常数、Born 电荷、压电张量与局域场
- [VASP Wiki: LCALCEPS](https://www.vasp.at/wiki/index.php/LCALCEPS)、[EFIELD_PEAD](https://www.vasp.at/wiki/index.php/EFIELD_PEAD) — 有限电场法与 Γ 中心网格要求
- [VASP Wiki: Category Dielectric properties](https://www.vasp.at/wiki/index.php/Category:Dielectric_properties) — 各种方法的总览、离子贡献、density-density vs current-current、EELS、吸收、反射率
- [VASP Wiki: Practical guide to GW calculations](https://www.vasp.at/wiki/index.php/Practical_guide_to_GW_calculations)、[Bethe-Salpeter-equations calculations](https://www.vasp.at/wiki/index.php/Bethe-Salpeter-equations_calculations)
- [VASP Known issues](https://www.vasp.at/wiki/index.php/Known_issues) — issue 20、44、92、93、100、114
- [VASP 官方教程：线性响应 / GW / BSE](https://www.vasp.at/tutorials/latest/)
- [py4vasp: dielectric_function](https://vasp.at/py4vasp/latest/calculation/dielectric_function/)
- [sumo-optplot 文档](https://smtg-bham.github.io/sumo/sumo-optplot.html)、[VASPKIT 教程（Linear Optical Properties）](https://vaspkit.com/tutorials.html)

**相关 GitHub 仓库（仅链接）**

- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) — `sumo-optplot` 绘制吸收谱、损失谱、n/k
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) — 官方 Python 后处理库（介电函数、介电张量、Born 电荷）
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — `Vasprun.dielectric`、`optical_absorption_coeff`、`Outcar.read_lepsilon`
- [raman-sc/VASP](https://github.com/raman-sc/VASP) — 有限位移极化率法计算非共振 Raman 活性（光学/声子扩展）
- [JMSkelton/VASP-Examples](https://github.com/JMSkelton/VASP-Examples) — SnS2 的 PBE/HSE06/G0W0 对比输入，可作 GW 入门参考
- [kavanase/vaspup2.0](https://github.com/kavanase/vaspup2.0) — 含介电常数收敛测试的 bash 流程
- [haiw201/dft-tutorial](https://github.com/haiw201/dft-tutorial) — 中文 DFT+VASP 课程，含光学性质章节

**关键文献**

- M. Gajdoš, K. Hummer, G. Kresse, J. Furthmüller, F. Bechstedt, *Phys. Rev. B* 73, 045112 (2006)（VASP 中 PAW 线性光学性质的实现）. DOI: [10.1103/PhysRevB.73.045112](https://doi.org/10.1103/PhysRevB.73.045112)
- S. Baroni, R. Resta, *Phys. Rev. B* 33, 7017 (1986)（DFPT 介电常数）. DOI: [10.1103/PhysRevB.33.7017](https://doi.org/10.1103/PhysRevB.33.7017)
- X. Gonze, C. Lee, *Phys. Rev. B* 55, 10355 (1997)（Born 电荷、离子介电贡献、LO-TO）. DOI: [10.1103/PhysRevB.55.10355](https://doi.org/10.1103/PhysRevB.55.10355)
- R. W. Nunes, X. Gonze, *Phys. Rev. B* 63, 155107 (2001)（有限电场与 PEAD）. DOI: [10.1103/PhysRevB.63.155107](https://doi.org/10.1103/PhysRevB.63.155107)
- I. Souza, J. Íñiguez, D. Vanderbilt, *Phys. Rev. Lett.* 89, 117602 (2002)（有限电场下的电焓泛函）. DOI: [10.1103/PhysRevLett.89.117602](https://doi.org/10.1103/PhysRevLett.89.117602)
- D. Sangalli et al., *Phys. Rev. B* 95, 155203 (2017)（current-current 与 density-density 响应的陷阱）. DOI: [10.1103/PhysRevB.95.155203](https://doi.org/10.1103/PhysRevB.95.155203)
- M. Shishkin, G. Kresse, *Phys. Rev. B* 74, 035101 (2006)（VASP 中的 GW 实现）. DOI: [10.1103/PhysRevB.74.035101](https://doi.org/10.1103/PhysRevB.74.035101)
- J. Tauc, R. Grigorovici, A. Vancu, *phys. stat. sol.* 15, 627 (1966)（Tauc 图）. DOI: [10.1002/pssb.19660150224](https://doi.org/10.1002/pssb.19660150224)
- A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, G. E. Scuseria, *J. Chem. Phys.* 125, 224106 (2006)（HSE06）. DOI: [10.1063/1.2404663](https://doi.org/10.1063/1.2404663)

---

[← 上一章](13-elastic.md) | [返回目录](../README.md) | [下一章 →](15-md-mlff.md)
