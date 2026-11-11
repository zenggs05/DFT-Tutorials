# 第16章 电荷与成键分析：Bader、差分电荷、ELF、STM、COHP（Charge & Bonding Analysis with VASP）

本章讲怎样从 VASP 的电荷密度和波函数里提取化学信息：Bader 电荷（LAECHG + AECCAR + bader 程序）、差分电荷密度（charge density difference, CDD）、电子局域函数（ELF, LELF）、部分电荷与 STM 模拟（LPARD、PARCHG、py4vasp STM），以及基于 LOBSTER 的晶体轨道哈密顿布居分析（COHP / ICOHP）。每一部分都给出完整的 INCAR、命令行、输出文件的读法和"结果是否可信"的判据，配套算例为 CO/Pt(111)、Si、单层 MoS2 和 GaAs。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清楚 CHGCAR、CHG、AECCAR0/1/2、PARCHG、ELFCAR 里分别存的是什么、单位是什么、怎样读；
- 用 `LAECHG = .TRUE.` + `chgsum.pl` + `bader CHGCAR -ref CHGCAR_sum` 做正确的 Bader 电荷分析，并用 POTCAR 中的 `ZVAL` 换算成净电荷；
- 设计"三次单点、同一晶胞、同一网格"的差分电荷计算，画出三维等值面和沿 z 的平面平均 Δρ(z)、电荷位移曲线 ΔQ(z)；
- 计算和解读 ELF，知道它在 VASP 中的限制（`KPAR = 1`、不支持非共线）；
- 用 `LPARD`、`IBAND`/`EINT`、`NBMOD`、`KPUSE`、`LSEPB`/`LSEPK` 得到能带分解或能量窗口内的部分电荷，按 Tersoff–Hamann 近似模拟恒高/恒流 STM 图像（VASP ≥ 6.5 可直接用 `LPARDH5` + py4vasp）；
- 为 LOBSTER 准备 VASP 计算（`ISYM`、`NBANDS`、`LWAVE`、`NSW = 0`），写 `lobsterin`，读 `ICOHPLIST.lobster`、检查 spilling，并用 LobsterPy 自动分析。

| 你想回答的问题 | 推荐方法 | VASP 关键设置 | 后处理工具 | 本章小节 |
|---|---|---|---|---|
| 某个原子得失了多少电子？ | Bader 电荷 | `LAECHG`、`LCHARG`、细 FFT 网格 | `bader`、`chgsum.pl`、baderkit | §4 |
| 吸附/掺杂/界面处电子往哪里流动？ | 差分电荷密度 | 三次单点，结构与网格完全一致 | VESTA、vaspkit 314、`chgdiff.pl`、[chgdiff.py](../scripts/chgdiff.py) | §5 |
| 键是共价、离子还是金属性？孤对电子在哪？ | ELF | `LELF = .TRUE.`、`KPAR = 1` | VESTA、vaspkit 428/429 | §6 |
| 某条能带/缺陷态/HOMO 的空间分布？STM 像什么样？ | 部分电荷 PARCHG | `LPARD`、`IBAND`/`EINT`、`NBMOD` | VESTA、py4vasp、vaspkit 325 | §7 |
| 两个原子之间成键还是反键？键强对比？ | COHP / ICOHP | `ISYM`、`NBANDS`、`LWAVE`、`NSW = 0` | LOBSTER、LobsterPy、wxDragon | §8 |

> [!NOTE]
> 本章默认你已经掌握 [第2章](02-input-output-files.md) 的输入输出文件、[第3章](03-convergence.md) 的收敛测试，以及 [第10章](10-adsorption-catalysis.md) 中 CO/Pt(111) 吸附结构的优化。电荷分析几乎都是在**已经优化好的结构**上做单点计算（`NSW = 0`），结构没优化好时任何电荷数字都没有意义。

---

## 2. 原理速览

### 2.1 VASP 的"电荷密度"到底是什么

PAW 方法把真实波函数写成平滑的赝波函数加上原子球内的修正。VASP 输出的 CHGCAR 中存放的是**赝价电子密度加补偿电荷**，定义在细 FFT 网格（`NGXF × NGYF × NGZF`）上。几点必须记住：

1. **只有价电子**。芯电子不在 CHGCAR 里，所以 CHGCAR 积分得到的是 POTCAR 中 `ZVAL` 之和（即 OUTCAR 里的 `NELECT`），不是原子序数之和。
2. **存储的是 ρ(r)·V_cell**。CHGCAR 每个格点的数值等于电荷密度乘以晶胞体积，所以
   $$N_\mathrm{el}=\frac{1}{N_\mathrm{grid}}\sum_{\mathbf r} \mathrm{data}(\mathbf r),\qquad \rho(\mathbf r)=\frac{\mathrm{data}(\mathbf r)}{V_\mathrm{cell}}\;(e/\text{Å}^3)$$
   VESTA 打开 CHGCAR 时已经自动除以了体积。
3. **原子核附近不是真实密度**。赝化后某些原子（例如过渡金属、带 `_sv` 半芯的元素）的密度极大值可能不在原子核上，甚至出现"甜甜圈"形分布。这正是 Bader 分析必须用全电子参考密度的原因（§4）。
4. **自旋极化时有两块数据**。`ISPIN = 2` 时 CHGCAR 中第一块是总密度 ρ↑+ρ↓，第二块是磁化密度 ρ↑−ρ↓（不是 ρ↓）。

### 2.2 Bader（QTAIM）划分

Bader 的"分子中的原子"理论用电荷密度的**零通量面**（zero-flux surface，$\nabla\rho\cdot\mathbf n=0$）把空间分成若干盆地，每个盆地包含一个密度极大值（通常就是一个原子）。对盆地内的价电子密度积分，就得到该原子的 Bader 布居 $N_A$，净电荷为

$$q_A = Z_\mathrm{val}^{A} - N_A$$

其中 $Z_\mathrm{val}^{A}$ 是该元素 POTCAR 中的 `ZVAL`。Henkelman 组的 `bader` 程序用格点上的最陡上升路径划分盆地，结果依赖于网格密度，因此要做网格收敛测试。

### 2.3 差分电荷密度

把体系 AB 拆成 A、B 两个片段（例如 slab + 吸附分子），在**完全相同的晶胞和原子坐标**下各算一次：

$$\Delta\rho(\mathbf r)=\rho_{AB}(\mathbf r)-\rho_{A}(\mathbf r)-\rho_{B}(\mathbf r)$$

Δρ > 0 表示成键后电子聚集，Δρ < 0 表示电子流失。对平行于表面的平面取平均得到 $\Delta\lambda(z)$，再沿 z 积分得到电荷位移曲线

$$\Delta Q(z)=\int_0^{z}\Delta\lambda(z')\,A\,\mathrm dz'$$

ΔQ(z) 的极值点和符号直观地给出"有多少电子从哪一侧流向哪一侧"。

### 2.4 电子局域函数（ELF）

ELF 衡量在参考电子附近找到同自旋电子的概率被"压低"的程度（Becke–Edgecombe 定义）：

$$\mathrm{ELF}=\frac{1}{1+\chi^2},\quad \chi=\frac{D}{D_h},\quad D=\tau-\frac{1}{8}\frac{|\nabla\rho|^2}{\rho},\quad D_h=\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}$$

其中 $\tau=\frac12\sum_i|\nabla\psi_i|^2$ 是动能密度。ELF 取值 0–1：接近 1 表示电子高度局域（共价键、孤对电子），0.5 对应均匀电子气（金属键），很小的值出现在原子壳层之间或离子之间。

### 2.5 部分电荷与 Tersoff–Hamann STM

部分电荷只对选定的能带 n 和 k 点求和：

$$\rho_\mathrm{part}(\mathbf r)=\sum_{n\in\text{bands}}\sum_{\mathbf k\in\text{kpts}} w_{\mathbf k}\,|\psi_{n\mathbf k}(\mathbf r)|^2$$

Tersoff–Hamann 近似认为隧穿电流正比于针尖位置处、费米能级到偏压窗口内的局域态密度：

$$I(\mathbf r_0,V)\propto\int_{E_F}^{E_F+eV}\rho(\mathbf r_0,E)\,\mathrm dE$$

样品偏压为负（V < 0）时电子从样品的占据态流向针尖，成像的是**占据态**；偏压为正时成像的是**空态**。恒高模式取某一高度的 ρ(x, y)，恒流模式取某一等值面 ρ = ρ_iso 的高度 z(x, y)。

### 2.6 COHP / ICOHP

把能带结构能量按原子轨道对分解（Dronskowski–Blöchl）：

$$\mathrm{COHP}_{\mu\nu}(E)=H_{\mu\nu}\sum_{n,\mathbf k}w_{\mathbf k}\,c^{*}_{\mu,n\mathbf k}c_{\nu,n\mathbf k}\,\delta(E-\varepsilon_{n\mathbf k}),\qquad \mathrm{ICOHP}=\int^{E_F}\mathrm{COHP}(E)\,\mathrm dE$$

成键态的 $H_{\mu\nu}$ 为负，所以习惯上画 **−COHP**：正值（右侧）为成键，负值（左侧）为反键。−ICOHP 越大，这一对原子的共价作用越强（只是"提示"键强，不等于键能）。平面波没有局域轨道，LOBSTER 把 PAW 波函数投影到局域基组上再计算（projected COHP），投影损失用 spilling 衡量。

### 2.7 各种"原子电荷"的比较

| 方法 | 依据 | 在 VASP 中如何得到 | 特点 |
|---|---|---|---|
| Bader | 电荷密度拓扑 | `LAECHG` + bader 程序 | 只依赖密度，可比较不同体系；对网格敏感 |
| DDEC6 | 原子密度迭代划分 | AECCAR0/AECCAR2/CHGCAR + Chargemol | 对多孔材料、MOF 更"化学合理" |
| Mulliken / Löwdin | 局域基组布居 | LOBSTER 输出 `CHARGE.lobster` | 依赖基组，适合同一基组下的相对比较 |
| OUTCAR 中 `total charge` | Wigner-Seitz 球（RWIGS）内积分 | `LORBIT = 10/11` | **不是**原子电荷：球之间的电荷没被计入，总和不等于 NELECT |

> [!WARNING]
> 不同方法给出的数值可以相差很大（同一个氧原子，Bader 可能是 −1.2 e，Löwdin 可能是 −0.6 e）。论文中比较电荷时必须使用同一种方法、同一套 POTCAR 和同等精度，只比较**趋势和差值**。

---

## 3. 体积数据文件速查

| 文件 | 产生条件 | 内容 | 典型用途 |
|---|---|---|---|
| `CHGCAR` | `LCHARG = .TRUE.`（默认） | 赝价电子密度 ×V（含 PAW 补偿电荷的平面波部分），后附 augmentation occupancies；ISPIN=2 时第二块为磁化密度 | 续算、差分电荷、Bader 被积密度、能带非自洽 |
| `CHG` | 默认写出 | 与 CHGCAR 同格式但精度低、不含 augmentation 信息；MD/优化时每隔若干步写一次 | 快速可视化 |
| `AECCAR0` | `LAECHG = .TRUE.` | 芯电子密度（PAW 重构的全电子芯密度） | Bader 参考密度的一半 |
| `AECCAR1` | `LAECHG = .TRUE.` | 原子价电子密度叠加（proto-atomic） | 一般用不到 |
| `AECCAR2` | `LAECHG = .TRUE.` | 自洽后的"全电子"价电子密度 | Bader 参考密度的另一半；DDEC6 |
| `ELFCAR` | `LELF = .TRUE.` | ELF 值（0–1，不乘体积）；ISPIN=2 时先 ELF↑ 后 ELF↓ | 成键类型分析 |
| `PARCHG*` | `LPARD = .TRUE.` | 部分电荷 ×V，格式同 CHG（无 augmentation 信息） | 能带分解电荷、STM |
| `LOCPOT` | `LVTOT` / `LVHAR` | 局域势（eV） | 功函数，见 [第9章](09-surface-workfunction.md) |
| `vaspout.h5` | 编译时 `-DVASP_HDF5` | 6.5+ 可写入部分电荷（`LPARDH5`），py4vasp 读取 | 现代后处理 |

文件结构（以 CHGCAR 为例）：前面是 POSCAR 格式的结构块 → 空行 → 一行 `NGXF NGYF NGZF` → 共 NGXF·NGYF·NGZF 个数，x 变化最快、z 最慢 → augmentation occupancies →（ISPIN=2 时）第二组网格行与磁化密度。用 Python 读取最简单的办法是 ASE：

```python
from ase.calculators.vasp import VaspChargeDensity

vcd = VaspChargeDensity("CHGCAR")
atoms = vcd.atoms[-1]          # 结构
rho = vcd.chg[-1]              # 已除以体积的密度，形状 (NGXF, NGYF, NGZF)，单位 e/Å^3
print(rho.shape, rho.sum() * atoms.get_volume() / rho.size)   # 最后一个数应等于 NELECT
if vcd.is_spin_polarized():
    mag = vcd.chgdiff[-1]      # 磁化密度 rho_up - rho_down
```

> [!TIP]
> 想把 CHGCAR 拆成自旋向上和向下两部分（例如做自旋分辨的 Bader 分析），可以用上面的 `rho` 和 `mag` 计算 `(rho+mag)/2`、`(rho-mag)/2`，把结果赋给 `vcd.chg[-1]`，再令 `vcd.chgdiff = []`、`vcd.aug = ""`、`vcd.augdiff = ""`（去掉与总密度对应的 augmentation 块），最后用 `vcd.write("CHGCAR_up", format="chgcar")` 写出。这样得到的文件只能用于可视化和 Bader，不能用于续算。

---

## 4. Bader 电荷分析（CO/Pt(111)）

配套算例：[examples/16-bader-co-pt111](../examples/16-bader-co-pt111/README.md)。

### 4.1 为什么必须用 AECCAR0 + AECCAR2 作参考

`bader` 程序假设每个原子核位置有一个密度极大值。CHGCAR 中的赝价电子密度在原子核附近被"挖空"，零通量面会被放错位置，导致电荷分给错误的原子，甚至出现不属于任何原子的"伪盆地"。正确做法是：

1. 用 `AECCAR0 + AECCAR2`（芯密度 + 全电子价密度）这个**在原子核处有尖峰**的密度来确定盆地边界；
2. 在这些盆地里积分 **CHGCAR**（价电子），得到每个原子的价电子数；
3. 净电荷 = `ZVAL` − Bader 布居。

### 4.2 INCAR

在第10章优化好的 CO/Pt(111) 结构上做一次单点：

```fortran
SYSTEM  = CO/Pt(111) p(2x2) Bader single point

# electronic
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-6
ISMEAR  = 1          ! Methfessel-Paxton for the metal slab
SIGMA   = 0.1
LREAL   = Auto
ALGO    = Fast

# single point on the relaxed geometry
NSW     = 0
IBRION  = -1

# charge output for Bader
LCHARG  = .TRUE.     ! CHGCAR: valence density to be integrated
LAECHG  = .TRUE.     ! AECCAR0 (core) + AECCAR2 (valence, all-electron) -> reference
LWAVE   = .FALSE.

# finer FFT grid for the charge density (about 2x the default of your system);
# read the default values from a previous OUTCAR: grep "NGXF" OUTCAR
# NGXF = 160 ; NGYF = 160 ; NGZF = 768
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `LAECHG` | `.TRUE.` | 重构全电子密度，写出 AECCAR0/1/2。AECCAR0、AECCAR1 在开始时写出，AECCAR2 在自洽结束后写出 |
| `LCHARG` | `.TRUE.` | 被积分的价电子密度 |
| `PREC` | `Accurate` | 默认 `NGXF = 2 × NGX`，避免 wrap-around 误差 |
| `NGXF/NGYF/NGZF` | 默认值的约 1.5–2 倍 | 芯密度极尖锐，网格越细 Bader 越准；用两次不同网格的结果做收敛判断 |
| `NSW` / `IBRION` | `0` / `-1` | 单点 |
| `ISPIN` | 与研究体系一致 | 自旋极化时 bader 默认只分析第一块（总密度） |
| `LDIPOL` 等 | 与第10章保持一致 | 电荷分析不改变物理设置，只加输出关键词 |

> [!IMPORTANT]
> 细网格会显著增加内存和 CHGCAR/AECCAR 的文件大小（一个 160×160×768 的网格单个文件就有数百 MB）。不要把这些文件提交到 Git 仓库。

### 4.3 运行与后处理

```bash
# 1) 生成 POTCAR（元素顺序必须与 POSCAR 第 6 行一致：Pt C O）
cat $VASP_PP_PATH/potpaw_PBE/{Pt,C,O}/POTCAR > POTCAR
grep -E "TITEL|ZVAL" POTCAR

# 2) 运行 VASP
mpirun -np 32 vasp_std > vasp.out

# 3) 检查文件是否齐全
ls -lh CHGCAR AECCAR0 AECCAR2
grep -E "NGXF|NELECT" OUTCAR | head -3

# 4) 参考密度 = 芯 + 全电子价密度（VTST 脚本，输出 CHGCAR_sum）
chgsum.pl AECCAR0 AECCAR2

# 5) Bader 分析：用 CHGCAR_sum 定边界，积分 CHGCAR
bader CHGCAR -ref CHGCAR_sum

# 6) 输出：ACF.dat（每个原子）、BCF.dat（每个 Bader 极大值）、AVF.dat（原子与盆地的对应关系；部分版本文件名为 AtomVolumes.dat）
cat ACF.dat
```

`chgsum.pl` 来自 [VTST scripts](https://theory.cm.utexas.edu/vtsttools/scripts.html)（GitHub：[henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts)），`bader` 可执行文件从 [Henkelman 组官网](https://theory.cm.utexas.edu/henkelman/code/bader/) 下载（当前 v1.05）。不想用 Perl 的话，[chgdiff.py](../scripts/chgdiff.py) 也可以做 AECCAR0 + AECCAR2 的求和（`python ../scripts/chgdiff.py -h` 查看用法），或者使用纯 Python 的 [baderkit](https://github.com/SWeav02/baderkit)。

### 4.4 读懂 ACF.dat 并换算净电荷

ACF.dat 的格式如下（示意输出，数值仅供参考）：

```text
    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL
 --------------------------------------------------------------------------------
    1      0.0000      0.0000     10.0000    10.0213      1.2501     15.3302
  ...
   17      0.0000      0.0000     18.7260     2.9xxx      0.7xxx      5.xxxx
   18      0.0000      0.0000     19.8760     7.1xxx      0.9xxx     19.xxxx
 --------------------------------------------------------------------------------
    VACUUM CHARGE:               0.0000
    VACUUM VOLUME:               0.0000
    NUMBER OF ELECTRONS:       170.0000
```

- `CHARGE` 是 Bader 布居 $N_A$（价电子数），不是净电荷；
- `NUMBER OF ELECTRONS` 应等于 OUTCAR 中的 `NELECT`（本例 16×10 + 4 + 6 = 170，若用 `Pt_pv` 等其他 POTCAR 会不同），误差应小于约 0.01 e；
- `MIN DIST` 是原子到其盆地边界的最小距离，应大于 POTCAR 的芯半径 RCORE 的量级，太小说明网格或参考密度有问题。

配套脚本 `bader_net_charge.py` 从 POTCAR 读取每种元素的 `ZVAL`、从 POSCAR 读取原子数，输出每个原子的净电荷和分组求和（如 CO 分子整体带电多少）：

```bash
python bader_net_charge.py --acf ACF.dat --poscar POSCAR --potcar POTCAR --group 17-18
```

净电荷为负表示得到电子。对 CO/Pt(111)，常见的结论是 CO 整体从金属得到少量电子（π* 反馈），但具体数值依赖于覆盖度、泛函和网格，请以自己的收敛结果为准。

### 4.5 网格收敛与其他注意事项

```bash
# 网格收敛：同一结构用两套 NGXF/NGYF/NGZF 各算一次，比较目标原子的电荷
paste <(awk 'NR>2 && $1~/^[0-9]+$/ {print $1,$5}' grid1/ACF.dat) \
      <(awk 'NR>2 && $1~/^[0-9]+$/ {print $5}'    grid2/ACF.dat) | \
      awk '{printf "%4d %9.4f %9.4f %8.4f\n", $1, $2, $3, $3-$2}'
```

- 差值小于约 0.01–0.02 e 即可认为收敛；
- 带真空的 slab：bader 1.03 起 `-vac` 默认为 `off`，真空中的低密度会被分配给最近的原子盆地，一般不影响结论；
- 自旋极化：bader 只读第一块（总密度），需要自旋分辨时先按 §3 的方法拆分；
- `-b weight` 选项使用 Yu–Trinkle 权重法，对粗网格更稳健，可作为交叉检验；
- 想把每个原子的 Bader 体积导出可视化，用 `bader CHGCAR -ref CHGCAR_sum -p all_atom`（生成 `BvAtxxxx.dat`，可用 VESTA 打开）；
- vaspkit 508 可以把 Bader 结果转换为 PQR 格式，在 VMD 中按电荷着色。

---

## 5. 差分电荷密度（CO/Pt(111)）

配套算例：[examples/16-cdd-co-pt111](../examples/16-cdd-co-pt111/README.md)。

### 5.1 三次计算的设计原则

| 计算 | 目录 | 结构来源 | 要求 |
|---|---|---|---|
| AB：CO/Pt(111) | `AB/` | 第10章优化后的 CONTCAR | 单点 |
| A：Pt(111) slab | `slab/` | 从 AB 中**删掉 C、O**，其余坐标不动 | 不再优化 |
| B：CO 分子 | `co/` | 从 AB 中**只保留 C、O**，放在同一个晶胞中 | 不再优化 |

三次计算必须使用：

- **完全相同的晶胞**（一个字节都不能改）和原子坐标；
- 相同的 `ENCUT`、`PREC`，从而相同的 `NGXF/NGYF/NGZF`（最保险的做法是把 AB 的 OUTCAR 中的数值显式写进三个 INCAR）；
- 相同的 k 点网格、POTCAR、泛函、`ISMEAR/SIGMA`；
- 一般 `ISPIN` 也保持一致；若片段本身是磁性的（如 O2、单个过渡金属原子），可以分别做自旋极化计算，但要在报告中说明，并注意 CHGCAR 第二块是磁化密度，求差时只用第一块。

> [!CAUTION]
> 片段千万不要重新优化。一旦片段原子移动了哪怕 0.01 Å，差分图里就会在原子核附近出现大块正负交替的伪影，掩盖真正的成键信息。

### 5.2 INCAR（三次计算共用）

```fortran
SYSTEM  = CDD single point (same file for AB, slab and CO)
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-6
ISMEAR  = 1
SIGMA   = 0.1
LREAL   = Auto
ALGO    = Fast
NSW     = 0
IBRION  = -1
ISPIN   = 1
LCHARG  = .TRUE.
LWAVE   = .FALSE.
# Force identical fine grids in all three runs (copy from AB/OUTCAR):
# NGXF = 120 ; NGYF = 120 ; NGZF = 576
```

### 5.3 运行

配套的 `run_cdd.sh` 做了这几件事：

```bash
# 1) 用 split_fragments.py 从 AB 的 POSCAR 拆出 slab 与 CO 两个片段（坐标、晶胞完全不变）
python split_fragments.py POSCAR --keep Pt --out slab/POSCAR
python split_fragments.py POSCAR --keep C O --out co/POSCAR

# 2) 每个目录按自己的元素拼 POTCAR，复制 INCAR 与 KPOINTS
cat $VASP_PP_PATH/potpaw_PBE/{Pt,C,O}/POTCAR > AB/POTCAR
cat $VASP_PP_PATH/potpaw_PBE/Pt/POTCAR       > slab/POTCAR
cat $VASP_PP_PATH/potpaw_PBE/{C,O}/POTCAR    > co/POTCAR

# 3) 依次运行（或分别提交作业）
for d in AB slab co; do (cd $d && mpirun -np 32 vasp_std > vasp.out); done

# 4) 检查三者网格一致
grep "NGXF" AB/OUTCAR slab/OUTCAR co/OUTCAR | grep dimension
```

### 5.4 求差的四种方法

```bash
# 方法 1：VTST 脚本（CHGCAR1 - CHGCAR2 写入 CHGCAR_diff），两次相减
chgdiff.pl AB/CHGCAR slab/CHGCAR && mv CHGCAR_diff CHGCAR_AB-slab
chgdiff.pl CHGCAR_AB-slab co/CHGCAR && mv CHGCAR_diff CHGCAR_cdd

# 方法 2：chgsum.pl 带系数：CHGCAR1*fact1 + CHGCAR2*fact2
chgsum.pl AB/CHGCAR slab/CHGCAR 1.0 -1.0 && mv CHGCAR_sum tmp
chgsum.pl tmp co/CHGCAR 1.0 -1.0 && mv CHGCAR_sum CHGCAR_cdd

# 方法 3：vaspkit 314，按提示输入三个文件路径：AB/CHGCAR slab/CHGCAR co/CHGCAR
vaspkit      # -> 31 -> 314

# 方法 4：本例自带的 cdd_planar.py（只用 numpy），同时输出平面平均与 ΔQ(z)
python cdd_planar.py AB/CHGCAR slab/CHGCAR co/CHGCAR --axis z --out CHGCAR_cdd
```

仓库的 [chgdiff.py](../scripts/chgdiff.py) 提供同样的功能并会先检查网格与晶胞是否一致（`python ../scripts/chgdiff.py -h`）。

> [!WARNING]
> `chgdiff.pl` / `chgsum.pl` 对自旋极化的 CHGCAR 只适合处理第一块（总密度）。如果三个文件的 ISPIN 不同，先只保留第一块再求差。

### 5.5 可视化与定量

**三维等值面（VESTA）**：File → Open 打开 `CHGCAR_cdd`；Properties → Isosurfaces 设置等值面水平（例如 ±0.002 e/Å³，需要根据体系调整），正值（电子聚集）常用黄色，负值（电子流失）用青色。VESTA 的 "Edit → Edit Data → Volumetric data" 还可以做两组数据的加减，但要求网格一致。

**平面平均与电荷位移曲线**：`cdd_planar.py` 输出 `cdd_planar.dat`，三列为 z（Å）、Δλ(z)（e/Å）、ΔQ(z)（e）。在 CHGCAR 约定下，

$$\Delta Q(z_k)=\sum_{k'\le k}\frac{1}{N_z}\,\overline{\Delta\mathrm{data}}(z_{k'})$$

其中 $\overline{\Delta\mathrm{data}}$ 是该层所有格点的平均值。读图要点：

- 在 Pt 表面与 C 原子之间的区域 Δλ 有正有负，说明 σ 给予（CO → Pt）与 π 反馈（Pt → CO）同时存在；
- ΔQ(z) 在 CO 与表面之间某一位置的值，表示穿过该平面的净电荷转移量（示意意义大于绝对数值）；
- 积分到真空区 ΔQ 应回到约 0（总电荷守恒），否则说明三个计算的 NELECT 或网格不一致。

**二维切片**：vaspkit 310（2D Slice of Charge Density）或 py4vasp 的 `calc.density.to_contour(...)`（需要 vaspout.h5 中有密度数据）。

### 5.6 自旋密度

自旋极化体系常用"自旋密度"代替差分电荷：直接把 CHGCAR 第二块（ρ↑−ρ↓）作为体积数据画等值面即可。vaspkit 312（Spin Density）可直接提取；VESTA 打开自旋极化 CHGCAR 时会把两块都读进来，可以在等值面设置中切换。

---

## 6. 电子局域函数（ELF）

配套算例：[examples/16-elf-si](../examples/16-elf-si/README.md)。

### 6.1 INCAR

```fortran
SYSTEM  = Si ELF
PREC    = Accurate      ! avoid wrap-around errors in ELFCAR
ENCUT   = 400
EDIFF   = 1E-6
ISMEAR  = 0
SIGMA   = 0.05
NSW     = 0
LELF    = .TRUE.        ! write ELFCAR
KPAR    = 1             ! ELF is not implemented for KPAR > 1
NPAR    = 1             ! the LELF page of the VASP wiki asks for NPAR = 1 explicitly
LCHARG  = .TRUE.
LWAVE   = .FALSE.
```

| 关键词 | 值 | 说明 |
|---|---|---|
| `LELF` | `.TRUE.` | 自洽结束后写出 ELFCAR |
| `KPAR` | `1` | 否则报 `ELF: KPAR>1 not implemented` |
| `NPAR` | `1` | VASP Wiki 的 LELF 页面要求显式设置 `NPAR = 1`（等价于所有进程按平面波/FFT 并行），小体系影响不大，大体系会变慢 |
| 非共线 | 不支持 | `vasp_ncl` 中会给出 `WARNING: ELF not implemented for non collinear case`；SOC 体系用共线计算的 ELF 近似 |
| `PREC` | `Accurate` | ELFCAR 文件页建议避免 wrap-around（旧写法为 `PREC = High`） |

### 6.2 运行与可视化

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
mpirun -np 8 vasp_std > vasp.out
ls -lh ELFCAR
```

- **VESTA**：直接打开 ELFCAR，等值面设为 0.75–0.85，可以看到 Si–Si 键中点处的局域"哑铃"；用 (1 1 0) 晶面切片看二维分布更清楚。
- **vaspkit**：428/429 把 LOCPOT/ELFCAR 转换为 XSF 或 Gaussian cube 格式，便于在其他软件中处理。
- **Python**：ELFCAR 与 CHGCAR 格式相同，但数值**没有乘体积**，用 `VaspChargeDensity("ELFCAR")` 读入后要乘回体积才是原始 ELF 值（ASE 读入时统一除以了体积）。

### 6.3 解读与局限

- 共价晶体（Si、金刚石）：键中点 ELF ≈ 0.9 以上；
- 离子晶体（NaCl、MgO）：ELF 高值集中在阴离子周围，阴阳离子之间很低；
- 金属：ELF 在间隙区接近 0.5 的平台；
- ELF 依赖于赝势：VASP 的 ELF 在平面波网格上计算，原子核附近的值受 PAW 赝化影响，**只解读原子之间（价层）的区域**。VASP Known issues 中"ELF is POTCAR dependent"一条标注在 6.4 中已处理，旧版本更要注意这一点；
- ELF 不是电荷，不能积分得到电子数。需要定量时，结合 Bader 或 COHP。

---

## 7. 部分电荷与 STM 模拟（单层 MoS2）

配套算例：[examples/16-stm-mos2](../examples/16-stm-mos2/README.md)。

### 7.1 两步工作流

`LPARD` 是一个**后处理**步骤：它读入已经收敛的 WAVECAR，不做电子自洽，几乎瞬间完成。

```text
步骤 1：自洽计算  LWAVE = .TRUE.，NBANDS 足够（要成像空态时尤其重要）
步骤 2：部分电荷  复制 POSCAR / KPOINTS / POTCAR / WAVECAR，INCAR 加 LPARD 相关标签
                  ENCUT、ISPIN、ISYM、NBANDS 必须与步骤 1 相同
```

### 7.2 选择能带与 k 点的标签

| 标签 | 默认 | 作用 |
|---|---|---|
| `LPARD` | `.FALSE.` | 打开部分电荷计算；不支持非共线（`LNONCOLLINEAR = .TRUE.`） |
| `IBAND` | 未设置 | 能带编号列表，如 `IBAND = 20 21 22 23 45`；设置后 `NBMOD` 自动等于列表长度 |
| `EINT` | 未设置 | 能量窗口（eV）。两个数为区间；只写一个数时另一端取 E_F；`EINT` 设了而 `NBMOD` 没设时，按 `NBMOD = -2`（绝对能量）处理 |
| `NBMOD` | 见下 | `n`（>0）：用 IBAND 中的 n 条带；`0`：所有能带（含空带）；`-1`：所有占据带，写入 CHGCAR（不含 augmentation）而不是 PARCHG；`-2`：EINT 为绝对能量；`-3`：EINT 相对于 E_F |
| `KPUSE` | 未设置 | 只用列表中的 k 点；设置后所有 k 点权重被重置为 1，一般需要自洽和后处理都用 `ISYM = -1` 才正确 |
| `LSEPB` | `.FALSE.` | `.TRUE.` 时每条带单独输出 `PARCHG.nb.ALLK` |
| `LSEPK` | `.FALSE.` | `.TRUE.` 时每个 k 点单独输出 `PARCHG.ALLB.nk`；与 LSEPB 都为 `.TRUE.` 时输出 `PARCHG.nb.nk`。分开 k 点时权重同样被设为 1 |
| `LPARDH5` | `.FALSE.` | **VASP ≥ 6.5.0**，且编译时开启 HDF5：输出重定向到 vaspout.h5，可用 py4vasp 分析与模拟 STM |

只设置 `LPARD = .TRUE.` 而不给任何选择标签时，`NBMOD` 默认为 −1，VASP 只把所有占据态的价电子密度写入 CHGCAR，**不会产生 PARCHG**——这是新手最常见的"为什么没有 PARCHG"的原因。

### 7.3 自洽计算 INCAR（步骤 1）

```fortran
SYSTEM  = MoS2 monolayer SCF for partial charges
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.05
NBANDS  = 24          ! default is ~16 for Mo_pv+2S; extra empty bands for empty-state images
NSW     = 0
LWAVE   = .TRUE.      ! WAVECAR is the input of the LPARD step
LCHARG  = .TRUE.
LORBIT  = 11
```

### 7.4 部分电荷 INCAR（步骤 2）

占据态（模拟负偏压约 −1 V 的 STM）：

```fortran
SYSTEM  = MoS2 partial charge, occupied states E_F-1.0 ... E_F
ENCUT   = 400         ! identical to the SCF run
ISMEAR  = 0
SIGMA   = 0.05
NBANDS  = 24          ! identical to the SCF run
ISTART  = 1           ! read WAVECAR
LPARD   = .TRUE.
NBMOD   = -3          ! EINT relative to the Fermi energy
EINT    = -1.0 0.0
LSEPB   = .FALSE.
LSEPK   = .FALSE.
# LPARDH5 = .TRUE.    ! VASP >= 6.5 with HDF5: write to vaspout.h5 for py4vasp
```

空态（正偏压）：MoS2 是半导体，VASP 给出的 E_F 位于价带顶附近，所以空态窗口应从导带底开始。先查带隙：

```bash
grep E-fermi OUTCAR
# VASP >= 6.4.3 默认（BANDGAP = COMPACT）在 OUTCAR 中报告 VBM、CBM 和带隙；也可用 vaspkit 911
grep -E "val. band max|cond. band min|fundamental gap" OUTCAR | tail -3
```

再把窗口设为导带底以上约 1 eV，例如（PBE 单层 MoS2 带隙约 1.6–1.8 eV，仅为近似值，以自己的计算为准）：

```fortran
NBMOD   = -3
EINT    = 1.0 3.0     ! relative to E_F; must cover CBM ... CBM+~1 eV of YOUR calculation
```

也可以用 `NBMOD = -2` 直接给绝对能量，或用 `IBAND` 指定导带底的能带编号（本例用 `Mo_pv` + `S` 时 NELECT = 12 + 2×6 = 24，占据 12 条带，导带底为第 13 条；换用 `Mo_sv` 时 NELECT = 26。务必用 `grep NELECT OUTCAR` 核对）。

### 7.5 输出文件与检查

```bash
mpirun -np 8 vasp_std > vasp.out     # LPARD 步骤很快，也可以用更少的核
ls PARCHG*
head -12 PARCHG                      # 结构块 + 网格行
```

- `LSEPB = LSEPK = .FALSE.` → 单个 `PARCHG`；
- `IBAND = 9 10 11`、`KPUSE = 1 34`、`LSEPB = LSEPK = .TRUE.` → `PARCHG.0009.0001`、`PARCHG.0009.0034`……
- 自旋极化时，PARCHG 第一块为选中态的总密度，第二块为其磁化密度。要单独看某一自旋通道，可以用 `EINT` 给足够多的小数位只框住一个自旋的本征值，或者对两块做加减。
- 如果 PARCHG 全是 0 或非常小：窗口里没有态（检查 E_F、NBMOD 的绝对/相对含义）、NBANDS 不够（空态）、或者 `EINT` 两端写反（VASP 会自动交换，但窗口仍可能落在带隙里）。

### 7.6 模拟 STM 图像

**方法 A：VASP ≥ 6.5 + py4vasp（推荐）**。在步骤 2 中加入 `LPARDH5 = .TRUE.`，然后：

```python
import py4vasp

calc = py4vasp.Calculation.from_path(".")
print(calc.partial_density.selections())
# constant-height mode: tip 3 A above the topmost atom, 6x6 tiling of the surface cell
stm = calc.partial_density.to_stm("constant_height", tip_height=3.0, supercell=6)
stm.to_plotly().write_image("stm_constant_height.png")   # needs the 'kaleido' package
# constant-current mode: isosurface defined by a tunnelling current (nA)
calc.partial_density.to_stm("constant_current", current=1.0, supercell=6).show()
```

`to_stm` 的 `selection` 还可以加自旋：`"constant_height(up)"`；平滑参数通过 `calc.partial_density.stm_settings`（`sigma_xy`、`sigma_z`、`truncate`、`enhancement_factor`、`interpolation_factor`）调节。数据组织在 vaspout.h5 的 `/results/partial_charges` 下，可以用 `h5ls -r vaspout.h5` 查看。

**方法 B：任意版本 + PARCHG**。本例的 `stm_from_parchg.py` 只依赖 numpy/matplotlib：

```bash
# 恒高：针尖在最高原子上方 3 Å
python stm_from_parchg.py PARCHG --mode height --height 3.0 --repeat 6 -o stm_h.png
# 恒流：取 rho = 1e-4 e/A^3 的等值面高度
python stm_from_parchg.py PARCHG --mode current --iso 1e-4 --repeat 6 -o stm_c.png
```

其他选择：vaspkit 325（STM Simulation）、VESTA 中对 PARCHG 做平行于表面的切片（恒高）或画等值面后从上方观察（恒流）。

> [!NOTE]
> 实验的偏压和针尖高度不能一一对应到模拟中。Tersoff–Hamann 近似忽略了针尖电子结构，真空中的密度衰减也受 DFT 交换关联势渐近行为的影响。模拟 STM 适合比较"亮点在哪个原子上、对称性和缺陷衬度"，不适合定量比较电流大小。

### 7.7 部分电荷的其他用途

- **缺陷态 / 杂质能级**：在 DOS 中找到带隙内的态，用 `IBAND` 画出它的空间分布，判断它局域在哪个原子；
- **分子轨道**：在大盒子中算分子（`vasp_gam`），用 `IBAND` 分别画 HOMO、LUMO；
- **能带边**：用 `KPUSE` 指定 VBM 所在的 k 点和 `IBAND` 指定价带顶，画出 VBM 的实空间轨道（记得 `ISYM = -1`）；
- 已知问题：VASP 6 之前的版本在 `LORBIT ≥ 11` 且 `ISYM = 2` 时部分电荷对称化不正确（Known issue 49），旧版本建议先对称计算自洽，再用 `ISYM = 0` 重新算部分电荷。

---

## 8. 晶体轨道哈密顿布居：LOBSTER COHP（GaAs）

配套算例：[examples/16-lobster-gaas](../examples/16-lobster-gaas/README.md)。

[LOBSTER](http://www.cohp.de/) 是 Dronskowski 组开发的免费学术软件（不允许用于商业研究；当前版本 6.0.0，2026-09-16；上一版 5.1.1），它读取 VASP 的 WAVECAR 等文件，把平面波投影到局域原子轨道基上，输出 pCOHP、pCOOP、COBI、Mulliken/Löwdin 电荷和投影 DOS。本节的 `lobsterin` 语法按 LOBSTER 5.x 编写，并与 pymatgen 的 `Lobsterin` 关键词表一致；升级到 6.0 后请对照新手册确认是否有变化。

### 8.1 VASP 计算的要求

| 要求 | 设置 | 原因 |
|---|---|---|
| 单点 | `NSW = 0`、`IBRION = -1` | LOBSTER 只分析一个结构的波函数 |
| 保存波函数 | `LWAVE = .TRUE.` | LOBSTER 读取 WAVECAR |
| 对称性 | `ISYM = 0`（VASP 6 常用）或 `-1` | LOBSTER 需要未经空间群对称约化的 k 点集合（`ISYM = 0` 时只按时间反演约化），pymatgen/LobsterPy 默认写 `ISYM = 0` |
| 能带数 | 手动设 `NBANDS` ≥ 局域基函数总数 | 投影需要足够多的能带；同时 NBANDS 要能被并行的能带组数整除，否则 VASP 会自动改动 |
| 赝势 | PAW（所有 VASP 5/6 的 POTCAR） | LOBSTER 的投影基于 PAW |
| 可执行程序 | `vasp_std` 或 `vasp_ncl`，**不能用 vasp_gam** | Gamma-only 版本的 WAVECAR 格式不受支持 |
| 输出 | 默认的 vasprun.xml、OUTCAR，以及 POSCAR、POTCAR、INCAR、KPOINTS 放在同一目录 | LOBSTER 从这些文件中读取结构、费米能和基本参数 |

**怎样确定 NBANDS？** 以 pbeVaspFit2015 基组为例，`Ga_d` 的推荐基函数是 3d 4s 4p（5 + 1 + 3 = 9 个），`As` 是 4s 4p（4 个），GaAs 原胞共 13 个。VASP 对 GaAs 的默认 NBANDS 约为 13，刚好卡在边界上，所以显式设 `NBANDS = 16`。如果用 32 个进程、`KPAR = 4`、`NCORE = 4`，能带组数为 32/4/4 = 2，16 能被 2 整除，VASP 不会再改动。

### 8.2 INCAR

```fortran
SYSTEM  = GaAs for LOBSTER
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.05         ! use the same width as gaussianSmearingWidth in lobsterin
NSW     = 0
IBRION  = -1
ISYM    = 0            ! full k-point set (time reversal only)
NBANDS  = 16           ! >= 13 local basis functions (Ga 3d4s4p + As 4s4p)
LWAVE   = .TRUE.       ! LOBSTER reads WAVECAR
LCHARG  = .TRUE.
NEDOS   = 2001
LORBIT  = 11           ! optional: compare VASP PDOS with LOBSTER DOS
```

### 8.3 lobsterin

```text
COHPstartEnergy -20.0
COHPendEnergy    5.0
gaussianSmearingWidth 0.05
basisSet pbeVaspFit2015
basisfunctions Ga 3d 4s 4p
basisfunctions As 4s 4p
cohpGenerator from 2.0 to 2.6 type Ga type As orbitalwise
saveProjectionToFile
```

| 关键词 | 说明 |
|---|---|
| `COHPstartEnergy` / `COHPendEnergy` | 输出 COHP 的能量范围（eV，相对 E_F）。Ga 3d 在价带顶以下约 15 eV 处，起点取 −20 才能包含它 |
| `gaussianSmearingWidth` | VASP 用 `ISMEAR = 0` 时设为与 SIGMA 相同；用 `ISMEAR = -5` 时 LOBSTER 使用四面体法，不写此行 |
| `basisSet pbeVaspFit2015` | 针对 PBE POTCAR 拟合的局域基组（覆盖到 Lr） |
| `basisfunctions` | 每种元素使用的轨道；要与 POTCAR 的价电子组态对应（`Ga_d` 包含 3d） |
| `cohpGenerator from a to b type X type Y` | 自动寻找距离在 a–b Å 之间的 X–Y 原子对；`orbitalwise` 额外输出轨道分辨的 COHP。GaAs 的 Ga–As 键长 = a√3/4 ≈ 2.45 Å（a = 5.653 Å） |
| `cohpbetween atom i atom j` | 手动指定原子对（编号按 POSCAR 顺序），可多行 |
| `saveProjectionToFile` | 保存投影，之后改分析参数时用 `loadProjectionFromFile` 跳过投影步骤 |

> [!TIP]
> 不想手写 lobsterin？在 VASP 输入目录中运行 `lobsterpy create-inputs`，LobsterPy 会根据 POSCAR/POTCAR/INCAR 生成 lobsterin 和修改后的 INCAR（目前仅支持 PBE_54 系列 POTCAR 的推荐基组）。pymatgen 的 `Lobsterin.standard_calculations_from_vasp_files(...)` 提供同样的功能，默认 `COHPstartEnergy = -35`、`COHPendEnergy = 5`、`cohpGenerator from 0.1 to 6.0 orbitalwise`。

### 8.4 运行 LOBSTER

```bash
# VASP
cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR
mpirun -np 16 vasp_std > vasp.out
grep -E "NBANDS|number of bands has been changed" OUTCAR | head

# LOBSTER（OpenMP 并行；可执行文件名以下载包为准，例如 lobster-5.1.1）
export OMP_NUM_THREADS=16
lobster > lobster.out
ls *.lobster lobsterout
```

主要输出：

| 文件 | 内容 |
|---|---|
| `lobsterout` | 运行日志，包含 spilling、警告 |
| `COHPCAR.lobster` / `ICOHPLIST.lobster` | COHP 曲线 / 每个原子对积分到 E_F 的 ICOHP |
| `COOPCAR.lobster` / `ICOOPLIST.lobster` | COOP（按重叠布居加权） |
| `COBICAR.lobster` / `ICOBILIST.lobster` | COBI（键指数，LOBSTER 4.1 起） |
| `DOSCAR.lobster` | 投影得到的 DOS（格式同 VASP DOSCAR） |
| `CHARGE.lobster` | Mulliken 与 Löwdin 电荷 |
| `GROSSPOP.lobster` | 轨道总布居 |
| `MadelungEnergies.lobster` | Madelung 能 |

### 8.5 检查投影质量

```bash
grep -i "spilling" lobsterout
```

示意输出（数值仅供参考）：

```text
abs. total spilling:   1.xx%
abs. charge spilling:  0.xx%
```

charge spilling 表示投影后丢失的价电子比例。经验上希望它在 1–2% 以内，越小越好；超过约 5% 时 COHP 的可靠性下降，需要换更大的基组（例如加入半芯或极化轨道）、检查 NBANDS 是否足够、确认 POTCAR 与基组匹配。`lobsterout` 中的 WARNING（例如能带数、k 点相关的提示）也要逐条看。

### 8.6 读 ICOHPLIST 与画 COHP

ICOHPLIST.lobster 每行给出一个原子对：编号、两个原子、距离、（新版本中的）晶格平移、积分到 E_F 的 ICOHP（eV）。示意输出（数值仅供参考）：

```text
  COHP#  atomMU  atomNU  distance   translation   ICOHP(ef)  for spin  1
      1     Ga1     As2   2.44785    0  0  0      -x.xxxxx
      2     Ga1     As2   2.44785   -1  0  0      -x.xxxxx
  ...
```

- ICOHP 为负表示净成键；比较键强时看 −ICOHP 的大小；
- 同一种化学键的多个对称等价原子对数值应相同，可以用来检查计算是否合理；
- `orbitalwise` 输出中还有 Ga 4s–As 4p 等轨道对的贡献。

画图：

```bash
pip install lobsterpy
lobsterpy description                 # 自动识别重要的阳离子-阴离子键并给出文字描述
lobsterpy plot-automatic              # 自动画出最重要的 COHP（别名 automatic-plot）
lobsterpy plot 1 2 --summed           # 画 COHPCAR.lobster 中第 1、2 个键的加和
lobsterpy plot-dos                    # 画 DOSCAR.lobster
lobsterpy description --help          # 查看全部参数
```

也可以用 Dronskowski 组推荐的 wxDragon，或用 pymatgen 的 `pymatgen.io.lobster.Cohpcar` / `Icohplist` 读入后自己作图。

**解读 −COHP 图**：E_F 以下的正峰为成键态，负峰为反键态。若 E_F 附近存在明显的占据反键态，说明该键被"过度填充"，常与结构不稳定、晶格软化或倾向于掺杂空穴有关；若 E_F 以上还有大量空的成键态，说明该键还有增强空间。

---

## 参数选择与收敛建议

| 分析 | 关键参数 | 推荐值 / 做法 | 收敛判据 |
|---|---|---|---|
| Bader | `NGXF/NGYF/NGZF` | 默认值的 1.5–2 倍；`PREC = Accurate` | 两套网格下目标原子电荷差 < 0.01–0.02 e；`NUMBER OF ELECTRONS` 等于 NELECT |
| Bader | 参考密度 | `AECCAR0 + AECCAR2` | `MIN DIST` 合理，无孤立的空盆地 |
| 差分电荷 | 结构 / 网格 | 片段坐标与晶胞完全相同，显式写 NGXF 等 | 真空区 Δρ ≈ 0；ΔQ 在真空区回到 0 |
| 差分电荷 | 等值面 | ±0.001–0.005 e/Å³ 之间试几次 | 图中特征不随等值面小幅变化而消失 |
| ELF | 并行 | `KPAR = 1`、`NPAR = 1` | 不报错；ELFCAR 存在 |
| ELF | 精度 | `PREC = Accurate`，k 点与能量计算一致 | 键区 ELF 形状稳定 |
| 部分电荷 | `NBANDS` | 成像空态时比默认多 30–50% | 窗口内能带完整 |
| 部分电荷 | k 点 | STM 需要收敛的 k 网格；`KPUSE` 时 `ISYM = -1` | 改变 k 网格图像不变 |
| STM | 针尖高度 / 等值面 | 恒高 2–4 Å；恒流 ρ_iso ≈ 1E-5–1E-3 e/Å³ | 主要亮点位置不随参数小幅变化 |
| LOBSTER | `NBANDS` | ≥ 基函数数，能被能带组数整除 | OUTCAR 中无 "number of bands has been changed" |
| LOBSTER | 基组 | pbeVaspFit2015 推荐基组，必要时加轨道 | charge spilling 约 < 1–2% |
| LOBSTER | k 点 | 与 DOS 计算相当的密度 | ICOHP 随 k 点变化 < 约 1% |

---

## 常见问题与报错

| 现象 / 报错信息 | 原因 | 解决办法 |
|---|---|---|
| Bader 电荷明显不合理（例如金属原子带 −3 e，或某个原子电荷远大于 ZVAL） | 没用 AECCAR0+AECCAR2 作参考，直接 `bader CHGCAR` | `LAECHG = .TRUE.`，`chgsum.pl AECCAR0 AECCAR2`，`bader CHGCAR -ref CHGCAR_sum` |
| 把 ACF.dat 的 `CHARGE` 当成了净电荷 | `CHARGE` 是价电子布居 | 净电荷 = POTCAR 中的 ZVAL − CHARGE（ZVAL 不是原子序数） |
| ACF.dat 的 `NUMBER OF ELECTRONS` ≠ NELECT | 网格太粗、真空区设了 `-vac`、读错了文件（用 AECCAR 作被积密度） | 加密 NGXF；被积密度必须是 CHGCAR |
| 加密网格后 Bader 结果仍在变 | 芯密度尖锐，网格还不够 | 继续加密，或用 `-b weight` 交叉检验 |
| 差分电荷在原子核附近出现大块正负伪影 | 片段坐标与 AB 不一致（重新优化过、用错了 CONTCAR） | 片段从 AB 的同一个 POSCAR 中删原子得到，不再优化 |
| `chgdiff.pl` / 自写脚本提示网格不一致（grid mismatch） | 三次计算的晶胞、ENCUT 或 PREC 不同 | 统一 INCAR，并显式写 `NGXF/NGYF/NGZF` |
| `WARNING: dimensions on CHGCAR file are different` | 续算或 ICHARG=1/11 时读入的 CHGCAR 网格与当前设置不同 | 保持 ENCUT/PREC/NGXF 一致，或删除旧 CHGCAR 重新自洽 |
| `ELF: KPAR>1 not implemented` | ELF 不支持 k 点并行 | `KPAR = 1`（并按 Wiki 设 `NPAR = 1`） |
| `WARNING: ELF not implemented for non collinear case` | `vasp_ncl` / `LNONCOLLINEAR` | 改用共线计算求 ELF |
| `LPARD = .TRUE.` 后没有 PARCHG，只有一个新的 CHGCAR | 没设 IBAND/EINT，`NBMOD` 默认为 −1 | 加 `IBAND` 或 `NBMOD = -3` + `EINT` |
| PARCHG 全为 0 或极小 | 窗口落在带隙中；`EINT` 被当作绝对能量（`NBMOD` 未设为 −3）；空态 NBANDS 不够 | 查 E-fermi 与带边；设 `NBMOD = -3`；增大 NBANDS 并重新自洽 |
| LPARD 步骤报 `while reading WAVECAR` 一类错误 | 后处理的 ENCUT、NBANDS、晶胞或 k 点与生成 WAVECAR 的计算不同 | 完全复制自洽计算的 INCAR，只增加 LPARD 相关标签 |
| 用 `KPUSE` 或 `LSEPK` 得到的密度对称性/大小不对 | k 点权重被重置为 1，而自洽用了对称约化 | 自洽和部分电荷都设 `ISYM = -1` |
| 非共线（SOC）体系的部分电荷算不出来 | LPARD 不支持非共线 | 用共线计算的部分电荷近似，或用 WAVECAR 解析工具（如 vaspkit 51x）处理 |
| LOBSTER 的 charge spilling 很大（> 5%） | 基组太小、NBANDS 不够、POTCAR 与基组不匹配 | 用推荐基组并适当增加轨道；增大 NBANDS；使用 PBE_54 推荐 POTCAR |
| OUTCAR 中出现 `The number of bands has been changed` | NBANDS 不能被能带组数整除，VASP 自动加大 | 调整 NBANDS 或并行参数（KPAR、NCORE）使之整除，保证与 lobsterin 的假设一致 |
| LOBSTER 无法读取 WAVECAR | 用了 `vasp_gam` | 改用 `vasp_std` 重新计算 |
| LOBSTER 提示 k 点/对称性相关问题 | VASP 使用了完整的空间群对称约化 | `ISYM = 0` 或 `-1` 重新计算 |

---

## 练习

1. **Bader 网格收敛**：对 CO/Pt(111) 分别用默认 NGXF 和 2 倍 NGXF 计算，比较 C、O 和与 C 相连 Pt 原子的 Bader 电荷，确定你需要的网格。
2. **参考密度的影响**：同一组文件分别运行 `bader CHGCAR` 和 `bader CHGCAR -ref CHGCAR_sum`，对比 Pt 原子的电荷，体会参考密度的作用。
3. **覆盖度效应**：把 p(2×2) 换成 p(3×3) 超胞，比较 CO 的总电荷和 ΔQ(z) 曲线的变化。
4. **ELF 比较**：分别计算 Si、MgO、Al 的 ELF，用 (110) 切片比较共价、离子、金属三类成键的特征。
5. **STM 衬度**：分别用 `EINT = -1.0 0.0`（占据态）和导带窗口（空态）模拟单层 MoS2 的 STM，观察亮点分别位于 S 还是 Mo 的位置，并与能带的轨道投影（[第5章](05-band-dos.md)）对照。
6. **COHP 与成键**：对 GaAs 画出 Ga–As 的 −COHP，确认 E_F 以下以成键态为主；再对比 GaN 或 InAs 的 −ICOHP，讨论键强趋势。

---

## 参考资料

**官方文档**

- VASP Wiki：[LAECHG](https://www.vasp.at/wiki/index.php/LAECHG)、[LELF](https://www.vasp.at/wiki/index.php/LELF)、[ELFCAR](https://www.vasp.at/wiki/index.php/ELFCAR)、[LPARD](https://www.vasp.at/wiki/index.php/LPARD)、[NBMOD](https://www.vasp.at/wiki/index.php/NBMOD)、[EINT](https://www.vasp.at/wiki/index.php/EINT)、[IBAND](https://www.vasp.at/wiki/index.php/IBAND)、[KPUSE](https://www.vasp.at/wiki/index.php/KPUSE)、[LSEPB](https://www.vasp.at/wiki/index.php/LSEPB)、[LSEPK](https://www.vasp.at/wiki/index.php/LSEPK)、[LPARDH5](https://www.vasp.at/wiki/index.php/LPARDH5)、[PARCHG](https://www.vasp.at/wiki/index.php/PARCHG)
- VASP Wiki：[Partial charge densities and STM simulations](https://www.vasp.at/wiki/index.php/Partial_charge_densities_and_STM_simulations)、[Known issues](https://www.vasp.at/wiki/index.php/Known_issues)
- VASP 官方表面教程（含恒高/恒流 STM）：<https://www.vasp.at/tutorials/latest/surface/part3/>
- py4vasp `partial_density`：<https://vasp.at/py4vasp/latest/calculation/partial_density/>
- Bader 程序（Henkelman 组）：<https://theory.cm.utexas.edu/henkelman/code/bader/>；VTST 电荷脚本说明：<https://theory.cm.utexas.edu/vtsttools/scripts.html>
- LOBSTER 官网与下载：<http://www.cohp.de/>；LobsterPy 文档：<https://digimatchem.github.io/LobsterPy/>
- VASPKIT 功能列表（31x 电荷、325 STM、42x 势/ELF 转换、508 Bader2PQR）：<https://vaspkit.com/features.html>

**相关 GitHub 仓库**（只链接，不转载代码）

- [henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts)：VTST 官方脚本，`chgsum.pl`、`chgdiff.pl`、`chg2cube.pl` 等
- [SWeav02/baderkit](https://github.com/SWeav02/baderkit)：Bader 分析的 Python 实现，便于在脚本中调用
- [DigiMatChem/LobsterPy](https://github.com/DigiMatChem/LobsterPy)：LOBSTER 输入生成、自动成键分析与作图
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp)：VASP 官方 Python 后处理，读取 vaspout.h5 模拟 STM
- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen)：`pymatgen.io.lobster` 读写 lobsterin、COHPCAR、ICOHPLIST
- [WMD-group/MacroDensity](https://github.com/WMD-group/MacroDensity)：体积数据的平面平均与宏观平均（与 §5 的平面平均思路相同）
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**

- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* 36, 354 (2006). DOI: [10.1016/j.commatsci.2005.04.010](https://doi.org/10.1016/j.commatsci.2005.04.010)
- E. Sanville, S. D. Kenny, R. Smith, G. Henkelman, *J. Comput. Chem.* 28, 899 (2007). DOI: [10.1002/jcc.20575](https://doi.org/10.1002/jcc.20575)
- W. Tang, E. Sanville, G. Henkelman, *J. Phys.: Condens. Matter* 21, 084204 (2009). DOI: [10.1088/0953-8984/21/8/084204](https://doi.org/10.1088/0953-8984/21/8/084204)
- M. Yu, D. R. Trinkle, *J. Chem. Phys.* 134, 064111 (2011)（weight 方法）. DOI: [10.1063/1.3553716](https://doi.org/10.1063/1.3553716)
- T. A. Manz, N. G. Limas, *RSC Adv.* 6, 47771 (2016)（DDEC6）. DOI: [10.1039/C6RA04656H](https://doi.org/10.1039/C6RA04656H)
- A. D. Becke, K. E. Edgecombe, *J. Chem. Phys.* 92, 5397 (1990). DOI: [10.1063/1.458517](https://doi.org/10.1063/1.458517)
- B. Silvi, A. Savin, *Nature* 371, 683 (1994). DOI: [10.1038/371683a0](https://doi.org/10.1038/371683a0)
- J. Tersoff, D. R. Hamann, *Phys. Rev. B* 31, 805 (1985). DOI: [10.1103/PhysRevB.31.805](https://doi.org/10.1103/PhysRevB.31.805)
- R. Dronskowski, P. E. Blöchl, *J. Phys. Chem.* 97, 8617 (1993). DOI: [10.1021/j100135a014](https://doi.org/10.1021/j100135a014)
- V. L. Deringer, A. L. Tchougréeff, R. Dronskowski, *J. Phys. Chem. A* 115, 5461 (2011). DOI: [10.1021/jp202489s](https://doi.org/10.1021/jp202489s)
- S. Maintz, V. L. Deringer, A. L. Tchougréeff, R. Dronskowski, *J. Comput. Chem.* 34, 2557 (2013). DOI: [10.1002/jcc.23424](https://doi.org/10.1002/jcc.23424)；*J. Comput. Chem.* 37, 1030 (2016). DOI: [10.1002/jcc.24300](https://doi.org/10.1002/jcc.24300)
- R. Nelson et al., *J. Comput. Chem.* 41, 1931 (2020)（LOBSTER 4）. DOI: [10.1002/jcc.26353](https://doi.org/10.1002/jcc.26353)
- J. George et al., *ChemPlusChem* 87, e202200123 (2022)（LobsterPy 自动分析方法）. DOI: [10.1002/cplu.202200123](https://doi.org/10.1002/cplu.202200123)
- K. Momma, F. Izumi, *J. Appl. Cryst.* 44, 1272 (2011)（VESTA）. DOI: [10.1107/S0021889811038970](https://doi.org/10.1107/S0021889811038970)

---

[← 上一章](15-md-mlff.md) | [返回目录](../README.md) | [下一章 →](17-advanced-electronic.md)
