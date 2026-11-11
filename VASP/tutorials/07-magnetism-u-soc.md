# 第7章 磁性、DFT+U 与自旋轨道耦合（Magnetism, DFT+U & Spin-Orbit Coupling with VASP）

本章系统讲解 VASP 中的磁性计算：自旋极化（ISPIN、MAGMOM）、铁磁/反铁磁构型与超胞、DFT+U（Dudarev 方案、LDAUL/LDAUU/LMAXMIX、线性响应 LDAUTYPE=3 求 U）、由磁构型能量差映射 Heisenberg 交换常数 J 并估算居里/奈尔温度，以及自旋轨道耦合（SOC，vasp_ncl、LSORBIT、SAXIS、GGA_COMPAT）和磁晶各向异性能（MAE）的力定理计算。所有步骤都配有可运行的 bcc Fe、NiO 算例（VASP magnetism / DFT+U / spin-orbit coupling / magnetic anisotropy tutorial）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 正确书写 `MAGMOM`（共线每原子 1 个数、非共线每原子 3 个数），理解"初始磁矩决定最终磁态"，并从 OSZICAR/OUTCAR 读出总磁矩与原子磁矩；
- 为反铁磁（AFM）构型建立磁性超胞，比较 FM / AFM / NM 能量；
- 用 Dudarev DFT+U 处理过渡金属氧化物（NiO），知道 `LMAXMIX`、`LDAUPRINT` 的作用和"不同 U 的能量不能比较"等陷阱；
- 用线性响应（`LDAUTYPE = 3`）从第一性原理求 U；
- 把多个磁构型的能量映射到 Heisenberg 模型，求 J1/J2 并粗估 T_C / T_N；
- 用 `vasp_ncl` 做 SOC 计算，用力定理（`ICHARG = 11` + 不同 `SAXIS`）计算 MAE；
- 了解约束磁矩 `I_CONSTRAINED_M` 的用途和基本写法。

| 场景 | 推荐做法 | 本章小节 |
|---|---|---|
| 含 Fe/Co/Ni/Mn/Cr 等的金属、合金 | `ISPIN=2` + 合理的 `MAGMOM`，比较多种初始磁态 | §3、§4 |
| 过渡金属氧化物、氟化物、f 电子体系 | DFT+U（`LDAUTYPE=2`），`LMAXMIX=4/6`，`LASPH` | §5、§6 |
| 需要第一性原理的 U | 线性响应 `LDAUTYPE=3`（或 cRPA） | §6 |
| 交换常数、居里温度 | 多磁构型能量映射 + 平均场/蒙特卡洛 | §7 |
| 重元素能带劈裂、Rashba、拓扑 | SOC：`vasp_ncl` + `LSORBIT` | §8 |
| 易磁化轴、MAE | 共线自洽 → 非自洽 SOC（力定理） | §9 |
| 固定磁矩方向/大小求能量曲线 | `I_CONSTRAINED_M` | §10 |

> [!NOTE]
> 本章假设你已经掌握 [第2章](02-input-output-files.md) 的输入文件、[第3章](03-convergence.md) 的收敛测试和 [第4章](04-relax-eos.md) 的结构优化。磁性计算对 ENCUT、k 点、smearing 的要求通常比非磁计算更苛刻，务必先做收敛测试。

---

## 2. 原理速览

### 2.1 共线自旋极化

自旋极化 DFT 把电子密度分成自旋向上和向下两部分：

$$
n(\mathbf r)=n_\uparrow(\mathbf r)+n_\downarrow(\mathbf r),\qquad
m(\mathbf r)=n_\uparrow(\mathbf r)-n_\downarrow(\mathbf r)
$$

交换关联能 $E_{xc}[n_\uparrow,n_\downarrow]$ 依赖两种自旋密度，两套 Kohn–Sham 方程分别求解（`ISPIN = 2`）。共线计算中磁矩只有"大小和正负"，没有空间方向——它与晶格的相对取向只有在加入 SOC 后才有意义。

关键事实：**自旋极化 DFT 的能量面常有多个自洽解**（FM、各种 AFM、低自旋/高自旋、非磁）。计算会落入哪个解，很大程度上取决于初始磁矩 `MAGMOM` 和对称性。这就是为什么磁性计算要"试多个初始构型、取能量最低者"。

### 2.2 DFT+U

半局域泛函对局域的 d/f 电子存在自相互作用误差，倾向于让 d/f 轨道"部分占据、离域"，导致过渡金属氧化物带隙严重偏小甚至成为金属。DFT+U 对选定原子的 d（或 f）壳层加一个 Hubbard 型修正。VASP 默认的 Dudarev 方案（`LDAUTYPE = 2`）：

$$
E_{\rm DFT+U}=E_{\rm DFT}+\frac{U-J}{2}\sum_{\sigma}\mathrm{Tr}\left[\rho^{\sigma}-\rho^{\sigma}\rho^{\sigma}\right]
$$

其中 $\rho^{\sigma}$ 是该原子 d 壳层的占据矩阵（由 PAW 投影子定义）。惩罚项在占据数为 0 或 1 时为零，在分数占据时为正，因此 U 会把占据态往下推、空态往上推（大约各 $U_{\rm eff}/2$），打开带隙、增大局域磁矩。Dudarev 方案只依赖 $U_{\rm eff}=U-J$；Liechtenstein 方案（`LDAUTYPE = 1`）则分别使用 U 和 J。

### 2.3 非共线磁性与自旋轨道耦合

非共线计算把每个 k 点的波函数写成两分量旋量，处理完整的 2×2 自旋密度矩阵，磁化密度 $\mathbf m(\mathbf r)$ 成为矢量。自旋轨道耦合（SOC）作为相对论修正：

$$
H_{\rm SOC}\propto \xi(r)\,\mathbf L\cdot\boldsymbol\sigma
$$

VASP 只在 PAW 球内计算 SOC 矩阵元（SOC 主要来自原子核附近）。SOC 把自旋自由度与晶格耦合起来，于是磁矩沿不同晶体方向时能量不同——这就是磁晶各向异性能（MAE）。

### 2.4 Heisenberg 模型与力定理

本章采用如下约定（经典单位矢量、每对原子只计一次、J>0 表示铁磁耦合）：

$$
E=E_0-\sum_{\langle ij\rangle}J_{ij}\,\mathbf e_i\cdot\mathbf e_j
$$

把几种共线磁构型的 DFT 总能代入这个式子，就能解出 J。MAE 则用"力定理"：固定共线计算得到的电荷/磁化密度，只改变磁化方向做一次非自洽 SOC 计算，两方向的能量差（近似等于占据态本征值之和的差）就是 MAE。

---

## 3. 自旋极化计算基础：bcc Fe 铁磁

完整输入见 [../examples/07-fe-bcc-fm/README.md](../examples/07-fe-bcc-fm/README.md)。

### 3.1 关键标签

| 标签 | 推荐值 | 说明 |
|---|---|---|
| `ISPIN` | 2 | 打开共线自旋极化（默认 1，非磁） |
| `MAGMOM` | 每原子 1 个数（μB） | 初始磁矩；默认对 `ISPIN=2` 为 NIONS*1.0。建议取预期值的 1.2–1.5 倍 |
| `LORBIT` | 11 | 在 OUTCAR 输出 `magnetization (x)` 表（按原子、按 s/p/d/f 分解），并写 lm 分解的 PROCAR/DOSCAR |
| `NUPDOWN` | 不设（默认全松弛） | 固定自旋向上与向下电子数之差，如 O₂ 三重态 `NUPDOWN = 2` |
| `AMIX_MAG` / `BMIX_MAG` | 默认 1.6 / 1.0 | 磁化密度的混合参数；难收敛时可改用线性混合（见 §3.5） |
| `ISYM` | 默认；必要时 0 | `MAGMOM` 会自动降低对称性；怀疑对称性把磁矩"抹掉"时关掉对称性 |
| `LASPH` | `.TRUE.` | 3d 元素、DFT+U、杂化泛函时推荐打开 |

### 3.2 MAGMOM 的写法

- 个数必须等于原子数（共线）或 3 × 原子数（非共线），**顺序与 POSCAR 中原子顺序一致**；
- 支持 `N*value` 简写和反斜杠续行：

```text
MAGMOM = 2*3.0                       # 2 Fe atoms, 3 muB each
MAGMOM = 2*3.0 2*-3.0 4*0            # 2 up, 2 down, 4 non-magnetic atoms (8 in total)
MAGMOM = 3.0 -3.0 3.0 -3.0 \
         16*0.0                      # 4 magnetic + 16 non-magnetic atoms
```

`N*value` 中 N 是重复次数（正整数），负磁矩写成 `2*-3.0`（星号后面是数值本身）。原子数不符时 VASP 报：

```text
Error reading item 'MAGMOM' from file INCAR.
```

> [!TIP]
> 原子多时手写 MAGMOM 很容易数错。可以用下面的 ASE 小脚本读 POSCAR、按元素自动生成（同号磁矩会被压缩成 `N*value`；AFM 的正负号仍需按子晶格自己指定）。

```python
from ase.io import read

atoms = read("POSCAR")
moments = {"Fe": 3.0, "Ni": 2.0, "Co": 2.0, "Mn": 4.0}      # initial guesses (muB)
mag = [moments.get(s, 0.0) for s in atoms.get_chemical_symbols()]

# compress consecutive equal values into N*value
out, prev, n = [], None, 0
for m in mag + [None]:
    if m == prev:
        n += 1
        continue
    if prev is not None:
        out.append(f"{n}*{prev}" if n > 1 else f"{prev}")
    prev, n = m, 1
print("MAGMOM = " + " ".join(out))
```

仓库中的 [../scripts/vasp_incar_writer.py](../scripts/vasp_incar_writer.py)（作者 Gaoshan Zeng）可以根据 POSCAR 为非磁、磁性和 DFT+U 三种静态计算写出 `INCAR.static`，磁矩以 `2*0_2_-2_6*0` 这样的下划线格式传入；具体参数请运行 `python ../scripts/vasp_incar_writer.py -h` 查看。

### 3.3 输入文件

POSCAR（2 原子 bcc 惯用胞，便于之后改成 AFM）：

```text
bcc Fe, 2-atom conventional cubic cell (a = 2.83 A, close to PBE)
2.83
   1.0000000000   0.0000000000   0.0000000000
   0.0000000000   1.0000000000   0.0000000000
   0.0000000000   0.0000000000   1.0000000000
Fe
2
Direct
   0.0000000000   0.0000000000   0.0000000000
   0.5000000000   0.5000000000   0.5000000000
```

KPOINTS：

```text
Gamma-centred 12x12x12 mesh for the 2.83 A cubic cell
0
Gamma
12 12 12
0 0 0
```

INCAR：

```text
SYSTEM  = bcc Fe ferromagnetic (FM), static run

# ---- basis / accuracy ----
PREC    = Accurate
ENCUT   = 400          # ~1.35 x ENMAX of Fe_pv (293 eV)
EDIFF   = 1E-6
NELM    = 100
LASPH   = .TRUE.       # non-spherical PAW gradient terms (recommended for 3d)
LREAL   = .FALSE.

# ---- smearing for a metal ----
ISMEAR  = 1            # Methfessel-Paxton order 1
SIGMA   = 0.10         # check T*S per atom < 1 meV in OUTCAR

# ---- spin polarisation ----
ISPIN   = 2
MAGMOM  = 2*3.0        # initial moments (muB); slightly above the expected ~2.2
LORBIT  = 11           # site/orbital projected moments -> "magnetization (x)"

# ---- static calculation ----
IBRION  = -1
NSW     = 0

LWAVE   = .FALSE.
LCHARG  = .TRUE.
```

| 关键词 | 值 | 为什么 |
|---|---|---|
| `ENCUT` | 400 | `Fe_pv` 的 ENMAX≈293 eV；静态计算 1.3 倍已足够，要做变胞优化同样适用 |
| `ISMEAR/SIGMA` | 1 / 0.10 | 磁性金属用 Methfessel–Paxton；SIGMA 过大会人为压低磁矩 |
| `MAGMOM` | 3.0 | 实验 Fe 磁矩约 2.2 μB，初值略大更容易落到高自旋铁磁解 |
| `LORBIT` | 11 | 不设 `LORBIT` 时 OUTCAR 中没有按原子分解的磁矩表 |

POTCAR 与运行：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR
mpirun -np 16 vasp_std > vasp.out 2>&1      # or: sbatch job.slurm
```

### 3.4 输出检查

OSZICAR 最后一行的 `mag=` 是晶胞总磁矩（μB）：

```bash
grep "F=" OSZICAR | tail -1
```

示意输出（数值仅供参考）：

```text
   1 F= -.16xxxxxxE+02 E0= -.16xxxxxxE+02  d E =0.000000E+00  mag=     4.4xxx
```

OUTCAR 中每个电子步都会打印总电子数与总磁化，最后一次即收敛值；`LORBIT = 11` 时末尾有按原子分解的磁矩表：

```bash
grep "number of electron" OUTCAR | tail -1
grep -A6 "magnetization (x)" OUTCAR | tail -7
```

示意输出（数值仅供参考）：

```text
 magnetization (x)

# of ion       s       p       d       tot
------------------------------------------
    1       -0.0xx  -0.0xx   2.3xx   2.2xx
    2       -0.0xx  -0.0xx   2.3xx   2.2xx
```

判断标准：

- 两个 Fe 的 `tot` 相同且约 2.2–2.3 μB（PBE 量级；实验约 2.22 μB/Fe）；
- `magnetization (x)` 的原子磁矩是 PAW 球内的积分，球外的磁化不计入，所以原子磁矩之和不严格等于 OSZICAR 的 `mag=`，这是正常的；
- 共线计算中只有 `(x)` 一张表；非共线计算会有 `(x)`、`(y)`、`(z)` 三张。

用 py4vasp（需要 HDF5 版 VASP 输出的 `vaspout.h5`）也能直接读局域磁矩：

```python
import py4vasp
calc = py4vasp.Calculation.from_path(".")
print(calc.local_moment.read())          # charges and magnetic moments per atom/orbital
```

### 3.5 让计算收敛到"正确"的磁态

1. **初始磁矩取大一点**：官方建议取实验值的 1.2–1.5 倍；不知道实验值时，可参考 Bilbao 磁结构数据库（MAGNDATA）或 Huebsch 等人的初始磁矩生成方案（PRX 11, 011031 (2021)）。
2. **从非磁解起步**：先跑 `ISPIN = 1` 得到 WAVECAR/CHGCAR，再以 `ICHARG = 1`、`ISPIN = 2` 和目标 `MAGMOM` 续算——此时 CHGCAR 里没有磁化密度，VASP 会用 MAGMOM 初始化磁矩。
3. **续算时不要删掉 MAGMOM**：从带磁化密度的 CHGCAR/WAVECAR 续算时，MAGMOM 不再决定磁矩大小，但仍参与对称性判断；删掉它，VASP 可能用更高的对称性把 AFM 磁矩"对称化"为零。
4. **调整混合参数**：对 slab、磁性体系和绝缘体，VASP Wiki 建议可尝试先用线性混合：

   ```text
   AMIX     = 0.2
   BMIX     = 0.0001    # almost zero, but 0 may crash some versions
   AMIX_MAG = 0.8
   BMIX_MAG = 0.0001
   ```

5. **固定总磁矩再放开**：对分子或自旋态明确的体系，先用 `NUPDOWN` 固定总磁矩收敛，再去掉 `NUPDOWN` 续算，检查是否保持在该自旋态。
6. **版本问题**：VASP 5.4.4–6.5.x 中，若 `MAGMOM` 含负值，默认 `NBANDS` 按带符号磁矩计算而偏少（known issue 51，6.6.0 修复），可能导致奇怪的崩溃。遇到时手动设置 `NBANDS`。

### 3.6 FM / AFM / NM 能量对比

算例目录中的 `run_fm_afm_nm.sh` 用同一 POSCAR 生成三种计算：

```bash
bash run_fm_afm_nm.sh            # FM: MAGMOM = 3 3 ; AFM: 3 -3 ; NM: ISPIN = 1
bash run_fm_afm_nm.sh collect
```

对 bcc Fe 应得到 E(FM) 最低。AFM（B2 型：体心原子与顶角原子反平行）能量更高，且局域磁矩往往明显减小甚至塌缩——这提醒我们：对巡游磁体，不同构型的磁矩大小可能差别很大，后面用 Heisenberg 模型映射 J 时要特别小心（§7）。

---

## 4. 反铁磁与磁构型：超胞的构建

反铁磁序的周期往往大于晶体学原胞，必须在**磁性晶胞**中计算。以岩盐结构 NiO 为例，实验基态是 AFM-II：沿 [111] 方向的 Ni 原子面内铁磁、相邻面反平行。

| 构型 | 磁性晶胞 | 晶格矢量（a 为立方晶格常数） | Ni 磁矩排列 |
|---|---|---|---|
| FM | 任意（原胞即可） | — | 全部平行 |
| AFM-I | 四方胞，2 个 NiO | a(½,½,0)、a(−½,½,0)、a(0,0,1) | 沿 [001] 逐层交替 |
| AFM-II | 菱方胞，2 个 NiO | a(1,½,½)、a(½,1,½)、a(½,½,1) | 沿 [111] 逐层交替 |

AFM-II 的 POSCAR（算例 [../examples/07-nio-afm-dftu/README.md](../examples/07-nio-afm-dftu/README.md)）：

```text
NiO AFM-II (rhombohedral 2 f.u. cell, moments alternate along [111]), a_cubic = 4.17 A
4.17
   1.0000000000   0.5000000000   0.5000000000
   0.5000000000   1.0000000000   0.5000000000
   0.5000000000   0.5000000000   1.0000000000
Ni O
2 2
Direct
   0.0000000000   0.0000000000   0.0000000000
   0.5000000000   0.5000000000   0.5000000000
   0.2500000000   0.2500000000   0.2500000000
   0.7500000000   0.7500000000   0.7500000000
```

配合 `MAGMOM = 2.0 -2.0 2*0.0`：第一个 Ni 向上、第二个 Ni 向下、两个 O 无磁矩。

构建磁性超胞的几点经验：

- **在 POSCAR 中把同一元素的原子按"自旋子晶格"排在一起**，MAGMOM 才好写；但不要为了区分子晶格而把同一元素写成两个物种——那样需要重复拼接 POTCAR，且 DFT+U 的 `LDAUL/LDAUU` 也要多写一列（线性响应求 U 时恰恰需要这样做，见 §6）。
- 用 ASE 的 `atoms.repeat()` 或 `make_supercell()` 扩胞后，按子晶格给 `set_initial_magnetic_moments()`，再用 §3.2 的脚本输出 MAGMOM 行。
- 需要系统枚举大量磁构型时，可用 pymatgen 的 `MagneticStructureEnumerator`（`pymatgen.analysis.magnetism.analyzer`，反铁磁枚举依赖 enumlib）。
- VASP 会根据 MAGMOM 自动去掉被磁序破坏的对称操作，一般无需手动设置 `ISYM`。

---

## 5. DFT+U：以 AFM NiO 为例

### 5.1 标签详解

| 标签 | 推荐值 | 说明 |
|---|---|---|
| `LDAU` | `.TRUE.` | 打开 DFT+U |
| `LDAUTYPE` | 2 | 1 = Liechtenstein（U、J 分别起作用）；2 = Dudarev（只有 U−J 起作用，默认）；3 = 线性响应用的附加势（§6）；4 = 同 1 但不含交换劈裂 |
| `LDAUL` | 每个物种一个数 | 加 U 的角动量：d = 2、f = 3、不加 = −1 |
| `LDAUU` | 每个物种一个数（eV） | U；Dudarev 方案中实际起作用的是 `LDAUU − LDAUJ` |
| `LDAUJ` | 每个物种一个数（eV） | J；Dudarev 方案通常写 0，把 U_eff 全部写进 `LDAUU` |
| `LMAXMIX` | 4（d）/ 6（f） | 让 d/f 的单中心占据矩阵参与密度混合并写入 CHGCAR；**DFT+U 必设**，否则收敛慢，且 `ICHARG = 11` 能带与自洽结果不一致 |
| `LDAUPRINT` | 1 | 每个电子步在 OUTCAR 输出占据矩阵（`onsite density matrix`）及总占据矩阵的本征值/本征矢（`occupancies and eigenvectors`） |
| `LASPH` | `.TRUE.` | 3d 氧化物与 DFT+U 计算的非球形梯度修正，官方强烈建议 |

> [!IMPORTANT]
> `LDAUL`、`LDAUU`、`LDAUJ` 是"每个**物种**（POSCAR 第 6 行的每一项）一个数"，不是每个原子一个数。POSCAR 写 `Ni O` 时三行各写两个数，顺序与 POSCAR、POTCAR 一致。

### 5.2 NiO 的 INCAR

```text
SYSTEM    = NiO AFM-II, PBE+U (Dudarev), static

# ---- basis / accuracy ----
PREC      = Accurate
ENCUT     = 520          # 1.3 x ENMAX(O = 400 eV); keep identical for all runs you compare
EDIFF     = 1E-6
NELM      = 150
LASPH     = .TRUE.       # important for 3d oxides and DFT+U
LREAL     = .FALSE.
ALGO      = Normal

# ---- insulator smearing ----
ISMEAR    = 0
SIGMA     = 0.05

# ---- spin: AFM-II, one up and one down Ni, O non-magnetic ----
ISPIN     = 2
MAGMOM    = 2.0 -2.0 2*0.0
LORBIT    = 11

# ---- DFT+U (Dudarev, LDAUTYPE = 2: only U_eff = U - J matters) ----
LDAU      = .TRUE.
LDAUTYPE  = 2
LDAUL     = 2 -1         # one value per species in POSCAR order: Ni (d) , O (none)
LDAUU     = 6.2 0.0      # eV; 6.2 eV = value used by Materials Project for Ni in oxides
LDAUJ     = 0.0 0.0
LDAUPRINT = 1            # write on-site occupation matrices to OUTCAR
LMAXMIX   = 4            # d electrons: needed for mixing and for later ICHARG=11 runs

# ---- static ----
IBRION    = -1
NSW       = 0
LWAVE     = .FALSE.
LCHARG    = .TRUE.
```

KPOINTS 用 Γ 中心 10×10×10（菱方胞），POTCAR 用 VASP 推荐的 `Ni_pv` 与 `O`：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Ni_pv,O}/POTCAR > POTCAR
grep TITEL POTCAR
mpirun -np 32 vasp_std > vasp.out 2>&1
```

> [!NOTE]
> POTCAR 的选择对 DFT+U 影响不小。VASP Wiki 的 NiO 示例表明：不加 U 时 `Ni`、`Ni_pv`、`Ni_sv_GW` 给出几乎相同的晶格常数；加较大的 U 后，没有半芯 3p 态的 `Ni` 势会明显低估体积，`Ni_pv` 好得多。DFT+U 计算请优先使用带半芯态的推荐势。

### 5.3 结果检查

```bash
grep "F=" OSZICAR | tail -1                     # mag= ~ 0 for AFM-II
grep -A8 "magnetization (x)" OUTCAR | tail -9   # Ni +-1.7 muB (approx.), O ~ 0
grep -A12 "occupancies and eigenvectors" OUTCAR | tail -14
```

- 两个 Ni 的磁矩大小相同、符号相反；PBE+U(6.2 eV) 一般在 1.6–1.8 μB 量级（示意；实验约 1.6–1.9 μB）。
- 占据矩阵本征值应接近 0 或 1（高度极化），Ni²⁺ 的 d⁸ 组态：自旋多数通道 5 个 d 轨道全满，少数通道 t2g 满、eg 空。
- 带隙：PBE 只有不到 1 eV；PBE+U 一般在 3 eV 左右（示意；实验光学带隙约 4 eV）。VASP ≥ 6.4.3 会在 OUTCAR 中报告带隙（`BANDGAP` 标签控制详略），也可按 [第5章](05-band-dos.md) 画 DOS。

### 5.4 U 值从哪里来

| 来源 | 做法 | 特点 |
|---|---|---|
| 经验/文献值 | 例如 Materials Project 对氧化物和氟化物的 U：Co 3.32、Cr 3.7、Fe 5.3、Mn 3.9、Mo 4.38、Ni 6.2、V 3.25、W 6.2 eV（pymatgen `MPRelaxSet`） | 便于与数据库对比；拟合自氧化反应焓（Wang, Maxisch, Ceder, PRB 73, 195107 (2006) 的思路） |
| 拟合实验 | 调 U 使带隙、磁矩、晶格常数或形成焓与实验一致 | 简单，但依赖所选的目标性质 |
| 线性响应 | `LDAUTYPE = 3`，见 §6 | 第一性原理；需要超胞收敛 |
| cRPA | 受限随机相近似计算屏蔽库仑相互作用（VASP 6.6.0 增加了谱 cRPA 等新功能） | 最严格，计算量大 |

### 5.5 必须牢记的陷阱

> [!WARNING]
> **不同 U 值下的总能量不能直接比较。** 能量随 U（Dudarev 方案下随 U−J）变化，比较反应能、吸附能、相稳定性时，所有参与比较的计算必须用完全相同的 U（同一元素在所有体系中用同一个 U），包括单质参考态。

- **GGA 与 GGA+U 混用**：例如计算 Fe₂O₃ 的形成焓时，Fe 金属通常用 GGA、Fe₂O₃ 用 GGA+U，两者能量尺度不同。Materials Project 用经验修正方案处理（pymatgen 的 `MaterialsProject2020Compatibility`，参见 Wang et al., Sci. Rep. 11, 15496 (2021)）。自己做时一定要说明处理方式。
- **能带计算**：自洽步骤就要写 `LMAXMIX = 4`（f 元素 6），再做 `ICHARG = 11`；或者直接用 `KPOINTS_OPT` 一步完成（VASP ≥ 6.3，见 [第5章](05-band-dos.md)）。
- **非自旋极化 + U**：`LDAUTYPE = 1/2/3` 也能用于 `ISPIN = 1`，但对磁性氧化物通常没有意义。

### 5.6 亚稳态与 U ramping

DFT+U 的能量面比纯 DFT 更"崎岖"：不同的 d 轨道占据方式（轨道序）对应不同的局部极小，SCF 很容易停在亚稳态，表现为同一结构、不同起点得到不同能量。常用对策：

1. **多起点**：改变 `MAGMOM`（高自旋/低自旋、不同 AFM 序）、改变 `ISYM`，取能量最低的解。
2. **U ramping**：先在 U = 0 收敛，然后逐步增大 U（如 0 → 2 → 4 → 6.2 eV），每一步读入上一步的 WAVECAR（`ISTART = 1`）。

   ```bash
   for U in 0.0 2.0 4.0 6.2; do
       sed "s/^LDAUU .*/LDAUU     = $U 0.0/" INCAR.template > INCAR
       mpirun -np 32 vasp_std > vasp_U$U.out 2>&1
       grep "F=" OSZICAR | tail -1
       cp OUTCAR OUTCAR_U$U          # WAVECAR is kept for the next step (LWAVE = .TRUE.)
   done
   ```

3. **检查占据矩阵**：`LDAUPRINT = 1` 输出的 `occupancies and eigenvectors` 能告诉你哪些轨道被占据，是判断是否落在预期轨道序的直接证据。
4. **占据矩阵控制（OMC）**：社区提供了强制初始占据矩阵的补丁（Allen & Watson, PCCP 16, 21016 (2014)），需要重新编译 VASP，适合系统地搜索轨道序。

---

## 6. 线性响应法求 U（LDAUTYPE = 3）

完整流程见 [../examples/07-linear-response-u/README.md](../examples/07-linear-response-u/README.md)。

### 6.1 原理

在超胞中只对一个原子 I 的 d 壳层施加一个小的球对称势 V，测量它的 d 占据数 N_I 的响应：

$$
\chi_0=\frac{\partial N_I^{\rm NSCF}}{\partial V},\qquad
\chi=\frac{\partial N_I^{\rm SCF}}{\partial V},\qquad
U=\chi^{-1}-\chi_0^{-1}
$$

非自洽响应 χ0（电荷密度固定，`ICHARG = 11`）反映 Kohn–Sham 轨道的"裸"重排；自洽响应 χ 包含电子屏蔽。两者倒数之差就是 DFT 本身"缺失"的那部分在位相互作用。

### 6.2 输入要点

1. **超胞 + 单独的物种**：在 2×2×2 的 AFM-II NiO 超胞（32 原子）中，把第 1 个 Ni 写成单独的物种（POSCAR 第 6、7 行为 `Ni Ni O` / `1 15 16`），POTCAR 相应拼接成 `Ni_pv, Ni_pv, O`。
2. **基态**（无 U）：`LORBIT = 11`（OUTCAR 末尾 `total charge` 表给出每个原子的 s/p/d 占据）、`LMAXMIX = 4`、保存 CHGCAR 和 WAVECAR。
3. **扰动计算**：

   ```text
   LDAU     = .TRUE.
   LDAUTYPE = 3
   LDAUL    = 2 -1 -1
   LDAUU    = 0.10 0.00 0.00     # strength V (eV) on the spin-up d manifold of species 1
   LDAUJ    = 0.10 0.00 0.00     # strength V (eV) on the spin-down d manifold
   ICHARG   = 11                 # only for the NON-self-consistent response
   ```

   `LDAUTYPE = 3` 时 `LDAUU`/`LDAUJ` 不再是 U/J，而是作用在自旋向上/向下 d 壳层上的附加势强度；两者取相同值即扰动总 d 占据。
4. **扫描 V**：例如 V = ±0.05、±0.10、±0.20 eV，每个 V 各做一次 NSCF（`ICHARG = 11`，读基态 CHGCAR）和一次 SCF（不设 `ICHARG = 11`）。
5. **线性拟合** N_d(V) 的斜率，得到 χ0、χ 和 U。

### 6.3 运行与分析

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Ni_pv,Ni_pv,O}/POTCAR > POTCAR
export NP=32 VASP_EXE=vasp_std
bash run_lr.sh gs
bash run_lr.sh response
python fit_u.py --gs gs/OUTCAR --dirs nscf_* scf_* --plot
```

示意输出（数值仅供参考）：

```text
nscf : points = 7, dN/dV =  0.49xxx 1/eV, max |residual| = 2.xxe-04
scf  : points = 7, dN/dV =  0.13xxx 1/eV, max |residual| = 1.xxe-04
U = 1/chi - 1/chi0 = 5.6 eV
```

VASP 官方 NiO 教程用这一方法、在相同大小的超胞中得到约 5.6 eV（线性拟合结果）。

> [!CAUTION]
> 线性响应 U 依赖于：超胞大小（扰动原子与其镜像的距离）、k 点、POTCAR（投影子定义了"d 占据"）、以及基态是 DFT 还是 DFT+U。报告 U 时请同时报告这些设置，并做超胞收敛。社区脚本 [sylearn/Ueff_VASP](https://github.com/sylearn/Ueff_VASP) 提供了另一套自动化实现（仅链接）。

---

## 7. 交换耦合常数 J 与居里/奈尔温度

### 7.1 先统一约定

文献中 Heisenberg 模型至少有三种常见写法，数值可以差 2 倍甚至 S² 倍：

| 写法 | 与本章 J 的关系 |
|---|---|
| $E=-\sum_{\langle ij\rangle}J_{ij}\mathbf e_i\cdot\mathbf e_j$（每对一次，单位矢量） | 本章约定 |
| $E=-\sum_{i\neq j}J'_{ij}\mathbf e_i\cdot\mathbf e_j$（每对计两次） | $J'=J/2$ |
| $E=-\sum_{\langle ij\rangle}\tilde J_{ij}\mathbf S_i\cdot\mathbf S_j$（含自旋长度 S） | $\tilde J=J/S^2$ |
| 符号相反的写法 $E=+\sum J\,\mathbf S_i\cdot\mathbf S_j$ | J 变号（J>0 为反铁磁） |

**比较文献数值之前，先把约定换算一致。**

### 7.2 能量映射的一般做法

1. 选 n 个磁构型（FM 加若干 AFM），**全部在同一结构、同一 ENCUT/k 点密度/U 下计算**；
2. 对每个构型数出每个磁性原子各近邻壳层中平行与反平行邻居的数目，写出 E = E0 − Σ c_k J_k；
3. 解线性方程（构型数多于未知数时用最小二乘）。

**例 1：bcc Fe 的有效 J（2 原子立方胞）。** FM 与 B2 型 AFM 中，每个 Fe 的 8 个最近邻全部平行/全部反平行，而次近邻（6 个）和第三近邻在两种构型中都平行、相互抵消，于是：

$$
J_{1}^{\rm eff}=\frac{E_{\rm AFM}-E_{\rm FM}}{8}\quad(\text{能量取每个 Fe 原子})
$$

> [!WARNING]
> bcc Fe 是巡游磁体，AFM 构型中磁矩明显变小，这种"大角度"映射得到的 J 只能作为粗略估计。对金属，更可靠的方法是基于格林函数的小角度方法（Liechtenstein 公式），例如 [TB2J](https://github.com/mailhexu/TB2J)（Wannier 函数或 LCAO 哈密顿量 → J_ij，He et al., Comput. Phys. Commun. 264, 107938 (2021)）。

**例 2：NiO 的 J1 和 J2。** Ni 有 12 个最近邻（J1，约 90° Ni–O–Ni 超交换）和 6 个次近邻（J2，180° 超交换）。按 §2.4 的约定，每个 Ni 的能量为：

| 构型 | 最近邻 Σ e·e | 次近邻 Σ e·e | 每个 Ni 的能量 |
|---|---|---|---|
| FM | +12 | +6 | E0 − 6J1 − 3J2 |
| AFM-I | −4 | +6 | E0 + 2J1 − 3J2 |
| AFM-II | 0 | −6 | E0 + 3J2 |

由此：

$$
J_1=\frac{E_{\rm AFM\text{-}I}-E_{\rm FM}}{8},\qquad
J_2=\frac{E_{\rm AFM\text{-}II}-E_{\rm FM}}{6}-J_1
$$

（能量均为每个 Ni 原子。）算例中的 `run_exchange.sh` 自动建立 FM（菱方胞）、AFM-II、AFM-I 以及"四方胞中的 FM"四个计算；最后一个用于检查两种晶胞的 k 点是否一致收敛——同一个 FM 态在两种晶胞中算出的每 Ni 能量应在 1 meV 以内相同：

```bash
cd ../examples/07-nio-afm-dftu
bash run_exchange.sh
bash run_exchange.sh collect
```

邻居计数可以自己用 ASE 验证（这也是检查 POSCAR/MAGMOM 是否写对的好办法）：

```python
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list

atoms = read("POSCAR")                         # AFM-II cell
m = np.array([1, -1, 0, 0])                    # same pattern as MAGMOM
i, j, d = neighbor_list("ijd", atoms, 4.2)
for name, lo, hi in (("J1", 2.9, 3.0), ("J2", 4.1, 4.2)):
    sel = (i == 0) & (atoms.numbers[j] == 28) & (d > lo) & (d < hi)
    print(name, "neighbours:", sel.sum(), " sum e.e:", (m[0] * m[j[sel]]).sum())
```

### 7.3 平均场估算 T_C / T_N

在本章约定下，经典自旋的平均场（MFA）结果为：

$$
k_BT_C^{\rm MFA}=\frac{1}{3}\sum_{j}J_{0j}\quad(\text{铁磁}),\qquad
k_BT_N^{\rm MFA}=\frac{1}{3}\Big|\sum_j J_{0j}\,s_j\Big|\quad(\text{反铁磁，}s_j=\pm1)
$$

对 NiO 的 AFM-II，最近邻贡献相互抵消，$k_BT_N^{\rm MFA}=-2J_2$。平均场通常**高估**转变温度（几十个百分点）；要得到可与实验（NiO 的 T_N ≈ 523 K）比较的数值，请把 J 交给经典蒙特卡洛程序（例如 [VAMPIRE](https://vampire.york.ac.uk/)）或用 RPA 近似。

---

## 8. 自旋轨道耦合（SOC）

### 8.1 程序与标签

SOC 计算必须使用 **`vasp_ncl`**（编译时不带 `-DNGZhalf`/`-DNGXhalf` 的版本，见 [第1章](01-setup-compile.md)）。用 `vasp_std` 运行非共线 INCAR 会直接报错退出，报错大意是"非共线计算要求程序在编译时不带 -DNGXhalf/-DNGZhalf"。

| 标签 | 推荐值 | 说明 |
|---|---|---|
| `LSORBIT` | `.TRUE.` | 打开 SOC，并自动设置 `LNONCOLLINEAR = .TRUE.`；只支持 PAW |
| `LNONCOLLINEAR` | （自动） | 非共线磁性；不加 SOC 的非共线磁结构单独设它即可 |
| `MAGMOM` | 每原子 3 个数 | 在 `SAXIS` 定义的自旋空间坐标系中给出；非磁体系写 `3*N*0` 形式，如 `MAGMOM = 24*0` 对 8 原子 |
| `SAXIS` | 默认 `0 0 1` | 全局自旋量子化轴（相对笛卡尔坐标）；VASP 读写的所有磁矩都在这个坐标系中 |
| `GGA_COMPAT` | `.FALSE.` | 恢复 GGA 梯度场的完整晶格对称性；对 MAE **强烈建议** |
| `LASPH` | `.TRUE.` | 官方推荐非共线计算打开 |
| `ISYM` | −1 或 0（视情况） | SOC 下建议测试关闭对称性是否改变结果；MAE 用 −1 保证各方向 k 点集合相同 |
| `LMAXMIX` | 4（d）/ 6（f） | 非共线 d/f 体系常需要，便于收敛和 `ICHARG = 11` 续算 |
| `NBANDS` | 共线计算的 2 倍 | 每个能带现在是旋量，能带数翻倍 |
| `LORBMOM` | `.TRUE.` | 输出轨道磁矩 `orbital moment (x)`（笛卡尔坐标系，仅 `LSORBIT = .TRUE.` 时可用） |
| `AMIX_MAG`/`BMIX_MAG` | 视情况 | 非共线体系难收敛时调整磁化密度混合 |

### 8.2 版本差异（务必核对你的 VASP 版本）

| 版本 | 行为 |
|---|---|
| < 6.5.0 | 非共线计算中 `ISPIN` 被内部忽略（设成 1），写 `ISPIN = 2` 不会出错 |
| 6.5.0 | `ISPIN = 2` 与 `LNONCOLLINEAR`/`LSORBIT` 同时出现时，因新的 MAGMOM 个数检查而报错停止（known issue 26） |
| ≥ 6.5.1 | `LNONCOLLINEAR` 优先，`ISPIN` 自动回退为 1 |
| 6.6.0 | 官方称大幅改善了非共线 GGA / meta-GGA 计算的稳定性 |
| ≤ 6.6.1 | `vasp_ncl` 中 Hartree–Fock/杂化泛函/GW 与 SOC（或不满足 k→−k 闭合的 k 点集）同时使用时，部分 PAW 单中心交换项有误（known issue 116，截至 2026-10 仍未修复） |

> [!TIP]
> 最省心的写法：**SOC/非共线 INCAR 里永远不写 `ISPIN`**，只写 `LSORBIT = .TRUE.`（或 `LNONCOLLINEAR = .TRUE.`）。这在所有版本中都正确。

### 8.3 MAGMOM 与 SAXIS：两种指定磁矩方向的方法

- **方法 A**：保持默认 `SAXIS = 0 0 1`，直接在笛卡尔坐标中写磁矩，例如沿 x 轴：`MAGMOM = 3 0 0`。
- **方法 B**：令 `SAXIS` 指向想要的方向，MAGMOM 只写第三个分量：`SAXIS = 1 1 1`、`MAGMOM = 0 0 3`。

两种方法原则上能量相同；VASP Wiki 指出方法 B 由于实现原因可能更精确，MAE 计算一般用方法 B。VASP 根据 SAXIS 计算两个欧拉角（OUTCAR 中可搜索 `Euler`），OUTCAR/PROCAR 中的磁矩都在旋转后的自旋坐标系中输出，而 `LORBMOM` 的轨道矩在笛卡尔坐标系中——比较两者时要做坐标变换。

> [!WARNING]
> 杂化泛函（`LHFCALC = .TRUE.`）会设置 `ISYM = 3`，而这只支持默认 `SAXIS`。杂化泛函 + SOC 时请保持 `SAXIS = 0 0 1`、改为旋转 MAGMOM（方法 A），或设置 `ISYM = -1`。

### 8.4 推荐的两步流程

非共线计算直接从头开始常常收敛很慢。官方推荐：

1. 先算非磁（或共线自旋极化）基态，得到 WAVECAR 与 CHGCAR；
2. 再用 `vasp_ncl` 读入，给出非共线 `MAGMOM`。共线计算的磁化被解释为沿 `SAXIS` 方向；只能整体旋转磁结构，不能在续算时单独旋转某些原子。

非磁半导体（如 GaAs、单层 MoS2）的 SOC 能带：

```text
SYSTEM     = GaAs with SOC, self-consistent step (vasp_ncl)
PREC       = Accurate
ENCUT      = 400
EDIFF      = 1E-7
ISMEAR     = 0
SIGMA      = 0.05
LSORBIT    = .TRUE.
MAGMOM     = 6*0          # 2 atoms x 3 components, non-magnetic
GGA_COMPAT = .FALSE.
LASPH      = .TRUE.
LMAXMIX    = 4            # Ga_d POTCAR contains 3d states
LORBIT     = 11
LCHARG     = .TRUE.
```

之后的能带步骤与 [第5章](05-band-dos.md) 相同（`ICHARG = 11` + Line-mode，或在本步直接加 `KPOINTS_OPT`），只是全程都用 `vasp_ncl`。自旋纹理、按自旋分量投影的能带可用 [pyprocar](https://github.com/romerogroup/pyprocar) 从 PROCAR 中绘制。

> [!NOTE]
> 对含反演对称性的非磁体系，SOC 不会打开能带的自旋简并；单层 MoS2 这类没有反演中心的体系，在 K 点会出现明显的价带自旋劈裂（约 0.15 eV 量级）。

---

## 9. 磁晶各向异性能（MAE）

完整输入见 [../examples/07-fe-bcc-soc-mae/README.md](../examples/07-fe-bcc-soc-mae/README.md)。

### 9.1 力定理流程

1. **共线高精度自洽**（`vasp_std`）：`PREC = Accurate`、`LREAL = .FALSE.`、`EDIFF = 1E-8`、`LMAXMIX = 4`、`GGA_COMPAT = .FALSE.`、`LCHARG = .TRUE.`，记下 NBANDS。
2. **非自洽 SOC**（`vasp_ncl`）：复制 CHGCAR，`ICHARG = 11`、`LSORBIT = .TRUE.`、`NBANDS` 翻倍、`ISYM = -1`，`MAGMOM = 0 0 m`，`SAXIS` 依次指向各个方向。
3. **MAE** = E(方向) − E(易轴)；两个方向的 SOC 计算除 `SAXIS` 外必须完全相同。

第二步的 INCAR（以 [111] 为例）：

```text
SYSTEM     = bcc Fe + SOC, spin quantisation axis [111] (run with vasp_ncl)

PREC       = Accurate
ENCUT      = 450
EDIFF      = 1E-8
LREAL      = .FALSE.
LASPH      = .TRUE.
GGA_COMPAT = .FALSE.
ISMEAR     = -5

LSORBIT    = .TRUE.          # switches on SOC and LNONCOLLINEAR; needs vasp_ncl
# Do NOT set ISPIN = 2 together with LSORBIT (VASP 6.5.0 stops; ignored otherwise)
SAXIS      = 1 1 1
MAGMOM     = 0 0 3.0         # 3 components per atom in the SAXIS frame -> moment along SAXIS
ICHARG     = 11              # keep the collinear CHGCAR fixed (force theorem)
ISYM       = -1              # same full k-point set for every spin direction
LMAXMIX    = 4
NBANDS     = 32              # 2 x NBANDS of the collinear run
LORBIT     = 11
LORBMOM    = .TRUE.          # orbital moments -> "orbital moment (x)" in OUTCAR

NSW        = 0
LCHARG     = .FALSE.
LWAVE      = .FALSE.
```

运行：

```bash
cd ../examples/07-fe-bcc-soc-mae
cat $VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR
export NP=8 VASP_STD=vasp_std VASP_NCL=vasp_ncl
bash run_mae.sh col
bash run_mae.sh soc
bash run_mae.sh collect
```

检查：

```bash
grep "NBANDS=" col/OUTCAR soc_001/OUTCAR | head -2      # doubled
grep -i "euler" soc_111/OUTCAR | head -2
grep -A4 "magnetization (z)" soc_111/OUTCAR | tail -2   # moment along SAXIS
grep "free  energy   TOTEN" soc_*/OUTCAR | awk '{print $1, $(NF-1)}'
```

### 9.2 为什么本例同时算 [001] 和 [100]

立方 bcc Fe 中 [001] 与 [100] 对称等价，E[100] − E[001] 理论上严格为 0，它的数值大小就是你当前 k 点、smearing、收敛标准下的**噪声水平**。bcc Fe 真正的各向异性 E[111] − E[001] 只有约 1 μeV/atom 量级（由实验各向异性常数估算），远小于很多常规设置下的噪声——这正是 Fe/Ni 的 MAE 被称为"计算难题"的原因（Halilov et al., PRB 57, 9557 (1998)）。

实用建议：

- 逐步加密 k 网格（16³ → 24³ → 32³ …），画 MAE–k 曲线，直到符号和大小稳定；
- 用 `ISMEAR = -5`（四面体）或很小 SIGMA 的 smearing，并在所有方向上保持一致；
- 练习时换成 MAE 大得多的体系：L1₀-FePt、CoPt（meV/f.u. 量级），或 CrI₃ 等二维铁磁体（面内 vs 面外）；
- 结果换算：MAE（J/m³）= MAE（eV/atom）× 1.602×10⁻¹⁹ / 原子体积（m³）。

### 9.3 自洽 SOC 的替代做法

也可以对每个方向做完全自洽的 SOC 计算（`ICHARG = 1` 或从头算），再比较总能。此时自旋会缓慢地向易轴转动，需要较松的收敛标准或配合约束磁矩（§10）把方向"钉住"。官方 Wiki 中两种做法都有示例，对大 MAE 体系两者结果通常一致。

---

## 10. 约束磁矩 I_CONSTRAINED_M（简介）

约束局域磁矩方法在总能中加一个罚函数，把每个原子在半径 `RWIGS` 球内积分得到的磁矩"拉"向 `M_CONSTR` 指定的方向（或方向+大小）。它只能用于非共线计算（`vasp_ncl`）。

| 标签 | 说明 |
|---|---|
| `I_CONSTRAINED_M` | 1 = 约束方向；2 = 约束方向和大小；4 = 约束方向与符号（VASP ≥ 6.4.0） |
| `M_CONSTR` | 每原子 3 个数（笛卡尔坐标）；某原子写 `0 0 0` 表示不约束 |
| `LAMBDA` | 罚函数权重；逐步增大（每次读上一次的 WAVECAR）可使罚能趋于零 |
| `RWIGS` | 每个物种的积分球半径（Å），**必须设置**；球之间尽量不重叠又尽量大 |

示例：bcc Fe 2 原子胞中，让两个磁矩在 xz 平面内夹角为 θ = 30°（用于由小角度能量求 J）：

```text
SYSTEM          = bcc Fe, two moments canted by 30 deg (vasp_ncl)
PREC            = Accurate
ENCUT           = 400
EDIFF           = 1E-7
ISMEAR          = 1
SIGMA           = 0.10
LNONCOLLINEAR   = .TRUE.           # no SOC needed for exchange energies
MAGMOM          = 0 0 2.2   1.1 0 1.905
I_CONSTRAINED_M = 1                # constrain directions only
M_CONSTR        = 0 0 1     0.5 0 0.866
LAMBDA          = 10               # increase stepwise (10 -> 20 -> 50) restarting from WAVECAR
RWIGS           = 1.20             # Fe sphere radius (A): < half the Fe-Fe distance (2.45 A)
LORBIT          = 11
GGA_COMPAT      = .FALSE.
LASPH           = .TRUE.
```

OSZICAR 中会多出 `E_p`（罚能）和每个原子的 `MW_int`/`M_int` 积分磁矩；`E_p` 应远小于你关心的能量差，否则继续增大 `LAMBDA`。注意：VASP 6.6.0 之前，GPU 版的 `I_CONSTRAINED_M` 计算可能崩溃或给出错误结果（6.6.0 修复），旧版本请在 CPU 上运行。

---

## 参数选择与收敛建议

| 任务 | 关键设置 | 建议 |
|---|---|---|
| 磁性金属（Fe/Co/Ni） | `ISMEAR = 1`、`SIGMA ≤ 0.1`、k 点比非磁计算更密 | 检查 T*S < 1 meV/atom；比较能量用 `energy(sigma->0)` |
| 初始磁矩 | 预期值 × 1.2–1.5 | 多起点比较，取能量最低者 |
| DFT+U | `LDAUTYPE = 2`、`LMAXMIX = 4/6`、`LASPH = .TRUE.`、带半芯态 POTCAR | 所有比较体系用同一 U；报告 U 的来源 |
| 线性响应 U | V = ±0.05…±0.2 eV、超胞 ≥ 2×2×2 | 检查线性度与超胞收敛 |
| 交换常数 J | 同一结构/ENCUT/k 点密度/U | 用"同一构型在不同晶胞"检查 k 点一致性（误差 < 1 meV/磁性原子） |
| SOC 能带 | `vasp_ncl`、`LSORBIT`、`MAGMOM` 3N 个数、`LMAXMIX` | 不写 `ISPIN`；NBANDS 约为共线的 2 倍 |
| MAE | `EDIFF = 1E-8`、`GGA_COMPAT = .FALSE.`、`ISYM = -1`、NBANDS×2、极密 k 点 | 用等价方向估计噪声；做 k 点收敛曲线 |
| 约束磁矩 | `RWIGS` 必设、`LAMBDA` 逐步增大 | `E_p` 远小于目标能量差 |

---

## 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `Error reading item 'MAGMOM' from file INCAR.` | MAGMOM 个数与原子数不符；非共线时没写 3×NIONS 个数；`N*value` 写法有误 | 逐个数清楚；用 §3.2 的脚本生成 |
| 设置 `ISPIN = 2` 与 `LSORBIT` 后报 MAGMOM 个数错误（仅 6.5.0） | known issue 26 | 删除 `ISPIN`，或升级到 ≥ 6.5.1 |
| 用 `vasp_std` 跑 SOC 立即退出，提示与 `-DNGXhalf/-DNGZhalf` 有关 | 非共线计算需要 `vasp_ncl` | 改用 `vasp_ncl` |
| 计算结束磁矩为 0（想要的是磁性解） | 初始 MAGMOM 太小；续算时删除了 MAGMOM 导致对称化；SIGMA 太大 | 增大初值；保留 MAGMOM；减小 SIGMA；尝试从非磁 CHGCAR 起步 |
| AFM 计算变成 FM（或反之） | 对称性/初值不合适；该 AFM 态本身不稳定 | 检查 MAGMOM 顺序；必要时 `ISYM = 0`；比较能量判断稳定性 |
| DFT+U 计算 SCF 振荡、不收敛 | 未设 `LMAXMIX`；亚稳态间来回跳 | `LMAXMIX = 4/6`；`ALGO = All` 或调小 `AMIX`；U ramping |
| DFT+U 的 `ICHARG = 11` 能带与自洽 DOS 对不上 | 自洽步骤没设 `LMAXMIX = 4/6` | 重做自洽（设置 LMAXMIX），或改用 KPOINTS_OPT |
| 同一结构 DFT+U 不同次计算能量相差几十 meV 以上 | 落入不同的轨道序（亚稳态） | 多起点 + 检查 `occupancies and eigenvectors`，取最低能量解 |
| MAE 结果随 k 点剧烈跳动、甚至变号 | k 点不足；不同方向 k 点集不同 | 加密 k 网格；`ISYM = -1`；用等价方向估计噪声 |
| SOC 计算中 NBANDS 被自动修改，且与共线不再是 2 倍关系 | 并行设置要求 NBANDS 为某数的整数倍 | 手动设 NBANDS 为合适的倍数，并在所有方向保持一致 |
| 自旋极化 + 负 MAGMOM 时奇怪的崩溃（< 6.6.0） | 默认 NBANDS 偏少（known issue 51） | 手动设置 NBANDS |
| 使用 6.6.0 做自旋极化 meta-GGA（如 r2SCAN）优化，受力异常 | known issue 94 | 升级到 6.6.1 或打官方补丁 |

---

## 练习

1. 用 [07-fe-bcc-fm](../examples/07-fe-bcc-fm/README.md) 分别对 FM 和 NM 态做 E–V 扫描（[第4章](04-relax-eos.md)），比较两者的平衡晶格常数，体会"磁体积效应"。
2. 在 NiO 中把 U 从 0 扫到 8 eV（U ramping），画出带隙和 Ni 磁矩随 U 的变化；思考为什么不能用这些计算的总能量画"能量–U"曲线来找最优 U。
3. 运行线性响应算例，比较 1×1×1 与 2×2×2 k 点、以及更大超胞得到的 U。
4. 用 `run_exchange.sh` 求 NiO 的 J1、J2 和平均场 T_N，并换算成"每对计两次"和"含 S = 1"的约定，与文献比较。
5. 构建 L1₀-FePt（c 轴为四方轴），计算 E[100] − E[001]，验证易轴沿 c 方向，并检查 MAE 随 k 网格的收敛。

---

## 参考资料

**VASP 官方文档**

- ISPIN / MAGMOM / NUPDOWN：<https://www.vasp.at/wiki/index.php/ISPIN>、<https://www.vasp.at/wiki/index.php/MAGMOM>、<https://www.vasp.at/wiki/index.php/NUPDOWN>
- DFT+U：<https://www.vasp.at/wiki/index.php/LDAU>、<https://www.vasp.at/wiki/index.php/LDAUTYPE>、<https://www.vasp.at/wiki/index.php/LMAXMIX>、<https://www.vasp.at/wiki/index.php/LDAUPRINT>
- 线性响应求 U：<https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA%2BU>
- POTCAR 选择（含 NiO + U 的赝势比较）：<https://www.vasp.at/wiki/index.php/Choosing_pseudopotentials>
- SOC 与非共线：<https://www.vasp.at/wiki/index.php/LSORBIT>、<https://www.vasp.at/wiki/index.php/LNONCOLLINEAR>、<https://www.vasp.at/wiki/index.php/SAXIS>、<https://www.vasp.at/wiki/index.php/GGA_COMPAT>、<https://www.vasp.at/wiki/index.php/LORBMOM>
- MAE 示例：<https://www.vasp.at/wiki/index.php/Determining_the_Magnetic_Anisotropy>
- 约束磁矩：<https://www.vasp.at/wiki/index.php/I_CONSTRAINED_M>
- 已知问题（issue 26、51、94、116 等）：<https://www.vasp.at/wiki/index.php/Known_issues>
- 官方磁性教程：<https://www.vasp.at/tutorials/latest/>

**相关 GitHub 仓库（仅链接）**

- [sylearn/Ueff_VASP](https://github.com/sylearn/Ueff_VASP)：线性响应（LDAUTYPE=3）自洽求 U 的流程脚本
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：`MPRelaxSet` 中的 U 值、`MagneticStructureEnumerator` 磁构型枚举、`MaterialsProject2020Compatibility` 能量修正
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar)：SOC 能带的自旋投影、自旋纹理
- [Ionizing/rsgrad](https://github.com/Ionizing/rsgrad)：快速查看 OUTCAR 中各步磁矩与收敛过程
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：`calc.local_moment` 读取局域磁矩
- [materialsproject/atomate2](https://github.com/materialsproject/atomate2)：含磁性构型枚举与排序的高通量工作流
- [mailhexu/TB2J](https://github.com/mailhexu/TB2J)：基于格林函数的交换参数计算

**关键文献**

- S. L. Dudarev et al., PRB 57, 1505 (1998). DOI: [10.1103/PhysRevB.57.1505](https://doi.org/10.1103/PhysRevB.57.1505)
- A. I. Liechtenstein, V. I. Anisimov, J. Zaanen, PRB 52, R5467 (1995). DOI: [10.1103/PhysRevB.52.R5467](https://doi.org/10.1103/PhysRevB.52.R5467)
- M. Cococcioni, S. de Gironcoli, PRB 71, 035105 (2005). DOI: [10.1103/PhysRevB.71.035105](https://doi.org/10.1103/PhysRevB.71.035105)
- L. Wang, T. Maxisch, G. Ceder, PRB 73, 195107 (2006). DOI: [10.1103/PhysRevB.73.195107](https://doi.org/10.1103/PhysRevB.73.195107)
- D. Hobbs, G. Kresse, J. Hafner, PRB 62, 11556 (2000)（PAW 非共线磁性）. DOI: [10.1103/PhysRevB.62.11556](https://doi.org/10.1103/PhysRevB.62.11556)
- S. Steiner et al., PRB 93, 224425 (2016)（PAW 中的 MAE 与 SOC）. DOI: [10.1103/PhysRevB.93.224425](https://doi.org/10.1103/PhysRevB.93.224425)
- S. V. Halilov et al., PRB 57, 9557 (1998)（Fe/Co/Ni 的 MAE）. DOI: [10.1103/PhysRevB.57.9557](https://doi.org/10.1103/PhysRevB.57.9557)
- P.-W. Ma, S. L. Dudarev, PRB 91, 054420 (2015)（约束磁矩泛函）. DOI: [10.1103/PhysRevB.91.054420](https://doi.org/10.1103/PhysRevB.91.054420)
- M. Pajda et al., PRB 64, 174402 (2001)（交换参数与居里温度）. DOI: [10.1103/PhysRevB.64.174402](https://doi.org/10.1103/PhysRevB.64.174402)
- M.-T. Huebsch et al., PRX 11, 011031 (2021)（磁结构预测基准与初始磁矩）. DOI: [10.1103/PhysRevX.11.011031](https://doi.org/10.1103/PhysRevX.11.011031)
- J. P. Allen, G. W. Watson, PCCP 16, 21016 (2014)（占据矩阵控制）. DOI: [10.1039/C4CP01083C](https://doi.org/10.1039/C4CP01083C)
- A. Wang et al., Sci. Rep. 11, 15496 (2021)（MP2020 能量修正）. DOI: [10.1038/s41598-021-94550-5](https://doi.org/10.1038/s41598-021-94550-5)
- X. He et al., Comput. Phys. Commun. 264, 107938 (2021)（TB2J）. DOI: [10.1016/j.cpc.2021.107938](https://doi.org/10.1016/j.cpc.2021.107938)

[← 上一章](06-hybrid-metagga.md) | [返回目录](../README.md) | [下一章 →](08-vdw-2d.md)
