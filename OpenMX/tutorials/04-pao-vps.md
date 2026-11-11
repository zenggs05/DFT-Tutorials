# 第4章 基组 PAO 与赝势 VPS 的选择（Choosing PAO Basis Sets & VPS Pseudopotentials in OpenMX）

OpenMX 用数值赝原子轨道（PAO, pseudo-atomic orbital）做基组，用模守恒赝势（VPS）描述离子实。基组需要用户自己挑，这一点和 VASP、QE 等平面波程序很不一样，也是 OpenMX 结果准不准、算得快不快的第一个决定因素。本章讲清楚 PAO 和 VPS 的命名规则、DFT_DATA19 数据库的 Quick/Standard/Precise 推荐表、半芯态与价电子数、同元素多 species、空原子（ghost atom）与 BSSE 修正，并用 Si 体相和水二聚体两个可运行算例演示基组收敛测试（OpenMX basis set convergence, PAO/VPS selection）。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 读懂 `Si7.0-s2p2d1`、`Fe5.5H-s3p2d1`、`C5.0-s4>1p4>1` 这类写法，并算出每个原子的基函数个数；
- 区分 `Fe_PBE19H` 与 `Fe_PBE19S`、`*_PBE19` 与 `*_CA19`，从 `.vps` 文件读出价电子数，正确填写坐标行的上/下自旋初始电子数；
- 按官方 DFT_DATA19 推荐表在 Quick / Standard / Precise 三档之间做取舍，知道 OpenMX 4.0 改了 Sn 的推荐；
- 理解为什么 PAO 和 VPS 必须来自同一版本数据库；
- 用同一元素定义多个 species，给不同区域分配不同精度的基组；
- 用空原子 `E` 做 counterpoise（CP）修正，估计基组重叠误差（BSSE）；
- 做基组收敛测试（能量、晶格常数、能带），并识别致密金属中的基组过完备（overcompleteness）问题。

| 场景 | 推荐做法 | 本章对应小节 |
|---|---|---|
| 第一次算某个材料 | 直接用推荐表 Standard 列 | §5 |
| 高通量筛选、大体系预弛豫 | Quick 列，最后用 Standard 复核 | §5、§9 |
| 发表用的能量差、弱相互作用、精细能带 | Standard 起步，用 Precise 验证 | §9 |
| 分子、表面、真空区域重要 | 选截断半径更大的 PAO 文件，必要时加空原子 | §5.3、§8 |
| 分子间结合能、吸附能 | CP 修正 BSSE | §8 |
| 致密 fcc/bcc/hcp 金属出现异常本征值 | 换小的优化基组或 `scf.ProExpn.VNA off` | §10 |
| 大体系只关心局部 | 同一元素定义多个 species，局部用大基组 | §7 |

> [!NOTE]
> 本章默认你已按 [第2章](02-install.md) 装好 OpenMX 4.0（含 patch 4.0.1），并在 [第3章](03-input-first-calc.md) 跑通过甲烷分子。截断能 `scf.energycutoff` 和 k 点的收敛测试放在 [第5章](05-scf-kpoints-convergence.md)，本章做基组测试时把它们固定住。

---

## 2. 原理速览

### 2.1 LCPAO：用局域数值轨道展开 Kohn–Sham 波函数

OpenMX 把 Kohn–Sham 轨道展开成以各原子为中心的数值轨道的线性组合（LCPAO 方法）。每个基函数写成径向部分与实球谐函数的乘积：

$$
\chi_{i\alpha}(\mathbf r)=R_{i\,l\,\zeta}(|\mathbf r-\mathbf R_i|)\,Y_{lm}(\widehat{\mathbf r-\mathbf R_i}),
\qquad R(r)=0\ \ (r>r_c)
$$

径向函数 $R$ 以数值表的形式存放在 `.pao` 文件里，超过截断半径 $r_c$ 严格为零。基组质量由两个旋钮控制：

1. **截断半径 $r_c$**：写在文件名里（如 `Si7.0` 表示 7.0 bohr）。$r_c$ 越大，轨道越"软"、越能描述长程尾巴，但每个原子的近邻数和矩阵非零元也越多；
2. **每个角动量的径向函数个数 $\zeta$**：写在 `-s2p2d1` 这一段。数字越大，描述极化和成键形变的自由度越多。

基函数总数为

$$
N_{\text{basis}}=\sum_{\text{atoms}}\ \sum_{l}\ n_l\,(2l+1),
$$

对角化代价约随 $N_{\text{basis}}^3$ 增长（O(N) 方法除外，见 第18章），所以基组多一档，时间和内存往往翻倍。

**原始轨道与优化轨道。** ADPACK 在限制势（confinement potential）下求解赝原子的基态和激发态，得到的径向函数称为原始 PAO（primitive）：基态无节点，第一激发态一个节点，依次类推。OpenMX 的轨道变分优化（orbital optimization）把若干原始函数线性组合成优化函数，组合系数叫收缩系数（contraction coefficient）。DFT_DATA19 中的 PAO 已经在 2–3 种化学环境（孤立原子、二聚体、体相或水分子等）下做过变分优化，再合并到同一个文件里，所以可移植性较好。数据库网站说明，每个 `.pao` 文件为 l = 0–3 的每个角动量各存了 15 个径向函数，按本征值从低到高排列；`s2` 就是取前两个 s 径向函数。

### 2.2 VPS：MBK 模守恒赝势，全相对论 j 依赖

DFT_DATA19 的 VPS 采用 Morrison–Bylander–Kleinman（MBK）方案生成的模守恒赝势，带非线性部分芯修正（PCC），并且是**全相对论、j 依赖**的：

- 非共线计算中打开 `scf.SpinOrbit.Coupling on`，就能直接自洽地计入自旋轨道耦合（见 [第12章](12-soc-topology.md)）；
- 共线或非自旋计算时，程序读入 j 依赖赝势后会按 j 的简并度自动做加权平均，变成标量相对论赝势。所以**同一个 VPS 文件既能算 SOC，也能算普通计算**，不需要另备一套。

### 2.3 和平面波程序的三点关键区别

| 方面 | 平面波（VASP/QE） | OpenMX（PAO） |
|---|---|---|
| 基组完备性 | 由 ENCUT 一个参数系统提高 | 由 $r_c$ 和 $\zeta$ 两个参数控制，不能无限加 |
| 基组随原子移动 | 不随原子移动，无 BSSE | 随原子移动，有 BSSE，需要 CP 修正（§8） |
| 截断能的含义 | 平面波基组截断 | `scf.energycutoff` 只决定实空间积分网格，**不是基组截断**（第5章） |
| 真空区 | 真空越厚越贵 | 真空区没有原子就没有基函数，几乎不增加对角化代价 |
| 过大基组 | 只会更慢 | 致密体系可能过完备，出现数值不稳定（§10） |

### 2.4 数据库目录 DFT_DATA19 与 DATA.PATH

```text
openmx4.0/
├── source/
├── work/                      # 官方示例输入，默认在这里运行
└── DFT_DATA19/
    ├── PAO/                   # H5.0.pao  H6.0.pao ... Si7.0.pao ... Fe5.5H.pao ...
    └── VPS/                   # H_PBE19.vps  H_CA19.vps ... Fe_PBE19H.vps  Fe_PBE19S.vps  E.vps ...
```

`DATA.PATH` 的默认值是相对路径 `../DFT_DATA19`，所以官方示例只能在 `work/` 里直接跑。在自己的目录里计算时，请写成**绝对路径**，指向同时包含 `PAO/` 和 `VPS/` 的那一层：

```text
DATA.PATH   /home/user/openmx4.0/DFT_DATA19     # 正确：指向 DFT_DATA19 本身
# DATA.PATH /home/user/openmx4.0/DFT_DATA19/PAO # 错误：不要指到 PAO 子目录
```

用官方 deb 包安装时，数据库在 `/opt/openmx4.0/DFT_DATA19`（本仓库示例的默认值就是它）。

---

## 3. PAO 命名规则详解

### 3.1 逐段拆解

以 `<Definition.of.Atomic.Species` 中的一行为例：

```text
<Definition.of.Atomic.Species
  Si   Si7.0-s2p2d1   Si_PBE19
Definition.of.Atomic.Species>
```

| 片段 | 含义 |
|---|---|
| 第1列 `Si` | species 名，任意取（如 `Si1`、`Si_bulk`），坐标块里用它引用 |
| `Si7.0` | PAO 文件名去掉扩展名，对应 `DATA.PATH/PAO/Si7.0.pao` |
| `7.0` | 生成该 PAO 时限制势的截断半径，单位 bohr（7.0 bohr ≈ 3.70 Å） |
| `-s2p2d1` | 从文件中取 2 个 s、2 个 p、1 个 d 径向函数 |
| 第3列 `Si_PBE19` | VPS 文件名去掉扩展名，对应 `DATA.PATH/VPS/Si_PBE19.vps` |

每个径向函数乘上 $2l+1$ 个实球谐函数：s 为 1 个，p 为 3 个，d 为 5 个，f 为 7 个。几个常见例子：

| 写法 | 计算 | 每原子基函数数 |
|---|---|---|
| `H6.0-s2p1` | 2×1 + 1×3 | 5 |
| `O6.0-s2p2d1`（或 `O7.0-s2p2d1`） | 2×1 + 2×3 + 1×5 | 13 |
| `Si7.0-s2p1d1` | 2 + 3 + 5 | 10 |
| `Si7.0-s2p2d1` | 2 + 6 + 5 | 13 |
| `Si7.0-s3p3d2` | 3 + 9 + 10 | 22 |
| `Fe5.5H-s3p2d1` | 3 + 6 + 5 | 14 |
| `Mo7.0-s3p2d2` | 3 + 6 + 10 | 19 |
| `Bi8.0-s3p2d2f1` | 3 + 6 + 10 + 7 | 26 |

> [!TIP]
> 截断半径不改变基函数个数：`O6.0-s2p2d1` 与 `O7.0-s2p2d1` 都是 13 个，只是轨道的空间范围不同。

### 3.2 收缩写法 `s4>1`

完整写法是"原始函数个数 `>` 收缩后函数个数"。官方手册的关键词说明里，`H5.0-s1>1p1>1` 与 `H5.0-s1p1` 等价（不收缩时可省略 `>1`）。轨道优化（`orbitalOpt.Method`，手册第26章）里会用到真正的收缩，例如：

```text
<Definition.of.Atomic.Species
  H   H5.0-s4>1        H_CA19      # 1 optimized s function built from 4 primitive s functions
  C   C5.0-s4>1p4>1    C_CA19
Definition.of.Atomic.Species>
```

这表示用 4 个原始 s 函数线性组合出 1 个优化 s 函数，基函数个数按 `>` 右边的数计算（H 为 1 个，C 为 1 + 3 = 4 个）。DFT_DATA19 的 PAO 已经优化过，日常计算**不需要**再写收缩。只有在给特殊体系自己优化基组时才用，优化结果写成 `C_1.pao` 之类的新文件（`CntOrb.fileout on`）。

### 3.3 用脚本数基函数、估计代价

下面的小函数解析 PAO 写法，并统计整个体系的基函数总数，可以在提交大计算前估算内存和时间：

```python
import re

L_DEG = {"s": 1, "p": 3, "d": 5, "f": 7}

def n_basis(spec: str) -> int:
    """'Si7.0-s2p2d1' -> 13; 'C5.0-s4>1p4>1' -> 4"""
    tail = spec.split("-", 1)[1]
    total = 0
    for l, n_prim, n_cont in re.findall(r"([spdf])(\d+)(?:>(\d+))?", tail):
        total += int(n_cont or n_prim) * L_DEG[l]
    return total

species = {"Fe": "Fe5.5H-s3p2d1", "O": "O6.0-s2p2d1"}
counts = {"Fe": 32, "O": 32}          # e.g. a 64-atom FeO supercell
nb = sum(n_basis(species[s]) * counts[s] for s in counts)
print("basis functions per cell:", nb)  # 32*14 + 32*13 = 864
```

本仓库的 `omx_make_input.py` 内置了 DFT_DATA19 推荐表和价电子数，可以从 POSCAR/CIF 直接生成带 `Definition.of.Atomic.Species` 的输入，用 `python ../scripts/omx_make_input.py -h` 查看用法。

### 3.4 查看 PAO 文件里有什么

`.pao` 是纯文本，文件末尾附近记录了最大角动量和每个 l 的径向函数个数，前部还保留了生成它的 ADPACK 参数：

```bash
DATA=/opt/openmx4.0/DFT_DATA19          # edit to your installation
ls $DATA/PAO | grep '^Si'                # Si7.0.pao  Si8.0.pao  Si9.0.pao
grep -E '^PAO\.Lmax|^PAO\.Mul|radial\.cutoff\.pao' $DATA/PAO/Si7.0.pao
```

对 DFT_DATA19 的 `Si7.0.pao`，可以看到 `PAO.Lmax 3`、`PAO.Mul 15`、`radial.cutoff.pao 7.0`，即 l = 0–3 每个角动量 15 个径向函数，截断 7.0 bohr。所以 `s2p2d2f1` 合法，但写 `g1` 会出错（文件里没有 l = 4）。

---

## 4. VPS 命名与价电子数

### 4.1 命名规则

| 写法 | 含义 |
|---|---|
| `C_PBE19` | GGA-PBE 赝势，2019 版数据库 |
| `C_CA19` | LDA（Ceperley–Alder）赝势，配 `scf.XcType LDA` 或 `LSDA-CA`/`LSDA-PW` |
| `Fe_PBE19H` | Hard：把更多半芯态当价电子（Fe：3s、3p、3d、4s，共 16 个电子） |
| `Fe_PBE19S` | Soft：较少半芯态（Fe：3p、3d、4s，共 14 个电子） |
| `Nd_PBE19_OC` | open-core：4f 电子冻结在芯内（部分镧系元素提供） |
| `E` | 空原子（empty atom）：零赝势、0 个价电子（§8） |

赝势和交换关联泛函要配套：`scf.XcType GGA-PBE` 用 `*_PBE19`，`LDA` 用 `*_CA19`，不要混用。

### 4.2 从 VPS 文件读价电子数

坐标行第 6、7 列是该原子上、下自旋的初始电子数，**两者之和必须等于 VPS 的价电子数**。不要凭元素周期表猜，直接查文件：

```bash
grep -E 'valence\.electron|xc\.type|j\.dependent' $DATA/VPS/Si_PBE19.vps
grep -E 'valence\.electron' $DATA/VPS/Fe_PBE19H.vps $DATA/VPS/Fe_PBE19S.vps
```

`Si_PBE19.vps` 中有 `valence.electron 4.0`、`xc.type GGA` 和 `j.dependent.pseudo.potentials on` 三行。常用元素的价电子数（取自官方手册表 1、表 2）：

| VPS | 价电子 | 包含的态 | 坐标行示例（非磁性起始） |
|---|---|---|---|
| `H_PBE19` | 1 | 1s | `0.5 0.5` |
| `C_PBE19` | 4 | 2s 2p | `2.0 2.0` |
| `O_PBE19` | 6 | 2s 2p | `3.0 3.0` |
| `Na_PBE19` | 9 | 2s 2p 3s | `4.5 4.5` |
| `Mg_PBE19` | 8 | 2p 3s | `4.0 4.0` |
| `Si_PBE19` | 4 | 3s 3p | `2.0 2.0` |
| `Ti_PBE19` | 12 | 3s 3p 3d 4s | `6.0 6.0` |
| `Fe_PBE19H` | 16 | 3s 3p 3d 4s | `8.0 8.0`（铁磁可用 `9.0 7.0`） |
| `Fe_PBE19S` | 14 | 3p 3d 4s | `7.0 7.0` |
| `Ni_PBE19H` / `Ni_PBE19S` | 18 / 16 | 含 / 不含 3s | `9.0 9.0` / `8.0 8.0` |
| `Cu_PBE19H` / `Cu_PBE19S` | 19 / 11 | 3s3p3d4s / 3d4s | `9.5 9.5` / `5.5 5.5` |
| `Ga_PBE19` | 13 | 3d 4s 4p | `6.5 6.5` |
| `Ge_PBE19` | 4 | 4s 4p | `2.0 2.0` |
| `As_PBE19` | 15 | 3d 4s 4p | `7.5 7.5` |
| `Mo_PBE19` | 14 | 4s 4p 4d 5s | `7.0 7.0` |
| `Pt_PBE19` | 16 | 5p 5d 6s | `8.0 8.0` |
| `Bi_PBE19` | 15 | 5d 6s 6p | `7.5 7.5` |

> [!IMPORTANT]
> 价电子里包含的半芯态必须有对应的基函数。例如 `Fe_PBE19H` 把 3s 和 4s 都当价电子，数据库说明它的最小基组是 `Fe*.*H-s2p1d1`（至少两个 s 函数）；`Fe_PBE19S` 的最小基组是 `Fe*.*S-s1p1d1`。`Ga_PBE19` 含 3d，所以推荐表里 Ga 的基组至少有 `d1`。

### 4.3 Hard 还是 Soft？

| 考虑 | 选 H（更多半芯态） | 选 S（较少半芯态） |
|---|---|---|
| 精度 | 半芯极化、高压、短键、磁性更可靠 | 常规结构和能带通常够用 |
| 代价 | 电子数和基函数更多，常需更高 `scf.energycutoff`（300–400 Ry，见第5章） | 更便宜 |
| 与文献对比 | 需与对方使用的价电子构型一致 | 同左 |
| DFT+U | `Hubbard.U.values` 要按所选 PAO 逐个列轨道（见 [第13章](13-dft-u.md)） | 同左 |

拿不准时，先用推荐表 Standard 列的那一个，并在论文方法部分写清楚 VPS 和 PAO 的全名。

---

## 5. DFT_DATA19 推荐表：Quick / Standard / Precise

### 5.1 三档的含义

官方手册第12.3节（Database Ver. 2019）为每个 VPS 给出三档基组：

- **Quick**：快速计算，精度尚可，适合预弛豫、高通量初筛；
- **Standard**：大多数用户的合适选择，兼顾精度和效率；
- **Precise**：多数情况下接近收敛的结果，用来验证 Standard。

数据库用 Δ 因子（delta gauge，与全电子 WIEN2k 状态方程的偏差）做过验证：71 个元素平均 1.774 meV/atom，标准差 1.702 meV/atom。各元素页面还给出了二聚体、体相 E–V 曲线和能带与全电子结果的对比，可以直接作为选基组的依据。

### 5.2 常用元素推荐表（摘自 OpenMX 4.0 手册表 1、表 2）

| VPS | 价电子 | Quick | Standard | Precise |
|---|---|---|---|---|
| `H_PBE19` | 1 | H5.0-s2 | H6.0-s2p1 | H7.0-s2p2d1 |
| `C_PBE19` | 4 | C6.0-s2p2 | C6.0-s2p2d1 | C6.0-s3p2d2 |
| `N_PBE19` | 5 | N6.0-s2p2 | N6.0-s2p2d1 | N6.0-s3p2d2 |
| `O_PBE19` | 6 | O6.0-s2p2 | O6.0-s2p2d1 | O6.0-s3p2d2 |
| `F_PBE19` | 7 | F6.0-s2p2 | F6.0-s2p2d1 | F6.0-s3p3d2f1 |
| `Na_PBE19` | 9 | Na9.0-s3p2 | Na9.0-s3p2d1 | Na9.0-s3p2d2 |
| `Mg_PBE19` | 8 | Mg9.0-s2p2 | Mg9.0-s3p2d1 | Mg9.0-s3p2d2 |
| `Al_PBE19` | 3 | Al7.0-s2p1d1 | Al7.0-s2p2d1 | Al7.0-s3p2d2 |
| `Si_PBE19` | 4 | Si7.0-s2p1d1 | Si7.0-s2p2d1 | Si7.0-s3p3d2 |
| `S_PBE19` | 6 | S7.0-s2p2d1 | S7.0-s2p2d1f1 | S7.0-s3p2d2f1 |
| `Cl_PBE19` | 7 | Cl7.0-s2p2d1 | Cl7.0-s2p2d1f1 | Cl7.0-s3p2d2f1 |
| `Ti_PBE19` | 12 | Ti7.0-s2p2d1 | Ti7.0-s3p2d1 | Ti7.0-s3p2d2f1 |
| `Mn_PBE19` | 15 | Mn6.0-s2p2d1 | Mn6.0-s3p2d1 | Mn6.0-s3p2d2f1 |
| `Fe_PBE19H` | 16 | Fe5.5H-s2p2d1 | Fe5.5H-s3p2d1 | Fe5.5H-s3p2d2f1 |
| `Fe_PBE19S` | 14 | Fe6.0S-s2p2d1 | Fe6.0S-s3p2d1 | Fe6.0S-s3p2d2f1 |
| `Co_PBE19H` | 17 | Co6.0H-s2p2d1 | Co6.0H-s3p2d1 | Co6.0H-s3p2d2f1 |
| `Ni_PBE19H` | 18 | Ni6.0H-s2p2d1 | Ni6.0H-s3p2d1 | Ni6.0H-s3p2d2f1 |
| `Ni_PBE19S` | 16 | Ni6.0S-s2p2d1 | Ni6.0S-s3p2d1 | Ni6.0S-s3p2d2f1 |
| `Cu_PBE19H` | 19 | Cu6.0H-s2p2d1 | Cu6.0H-s3p2d1 | Cu6.0H-s3p2d2f1 |
| `Cu_PBE19S` | 11 | Cu6.0S-s2p1d1 | Cu6.0S-s3p2d1 | Cu6.0S-s3p2d2f1 |
| `Ga_PBE19` | 13 | Ga7.0-s2p2d1 | Ga7.0-s3p2d2 | Ga7.0-s3p2d2f1 |
| `As_PBE19` | 15 | As7.0-s3p2d1 | As7.0-s3p2d2 | As7.0-s3p2d2f1 |
| `Se_PBE19` | 6 | Se7.0-s3p2d1 | Se7.0-s3p2d2 | Se7.0-s3p2d2f1 |
| `Mo_PBE19` | 14 | Mo7.0-s3p2d1 | Mo7.0-s3p2d2 | Mo7.0-s3p2d2f1 |
| `Sn_PBE19` | 14 | Sn7.0-s2p2d2 | Sn7.0-s3p2d3 | Sn7.0-s3p2d3f1 |
| `Ba_PBE19` | 10 | Ba10.0-s3p2d1 | Ba10.0-s3p2d2 | Ba10.0-s3p2d2f1 |
| `W_PBE19` | 12 | W7.0-s3p2d2 | W7.0-s3p2d2f1 | W7.0-s3p3d2f1 |
| `Pt_PBE19` | 16 | Pt7.0-s3p2d2 | Pt7.0-s3p2d2f1 | Pt7.0-s3p3d2f1 |
| `Au_PBE19` | 17 | Au7.0-s3p2d2 | Au7.0-s3p2d2f1 | Au7.0-s3p3d2f1 |
| `Pb_PBE19` | 14 | Pb8.0-s3p2d2 | Pb8.0-s3p2d2f1 | Pb8.0-s3p3d2f1 |
| `Bi_PBE19` | 15 | Bi8.0-s3p2d2 | Bi8.0-s3p2d2f1 | Bi8.0-s3p3d2f1 |
| `E`（空原子） | 0 | Kr10.0-s1p1 | Kr10.0-s2p1d1 | Kr10.0-s2p2d1f1 |

完整表格（H–Bi，含 Co/Ni/Cu/Zn 的 H/S 两种和部分镧系元素）见手册 [12.3 节](https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html)。

> [!NOTE]
> **版本差异**：OpenMX 4.0 修改了 Sn 的推荐基组，现为 `Sn7.0-s2p2d2`（Quick）、`Sn7.0-s3p2d3`（Standard）、`Sn7.0-s3p2d3f1`（Precise）。用 3.9 时代的输入（例如 Standard 写 `Sn7.0-s3p2d2`）复算时，结果会有差别。4.0 数据库还新增了 Tb、Yb 的 VPS/PAO 以及 XPS 芯能级计算用的 VPS/PAO。

### 5.3 截断半径怎么选

推荐表里每个元素只给一个半径，但数据库通常提供好几个 PAO 文件（例如 `H5.0/H6.0/H7.0`、`O5.0/O6.0/O7.0`、`Si7.0/Si8.0/Si9.0`、`Fe5.0H/5.5H/6.0H/8.0H/10.0H`）。数据库各元素页面的原话是：较小半径（如 `H5.0`、`O5.0`、`C5.0`、`Si7.0`）对体相已经足够，**分子体系**更适合较大半径（如 `H6.0`/`H7.0`、`O6.0`/`O7.0`、`Si8.0`/`Si9.0`）。规律可以概括为：

- **体相致密晶体**：用推荐表的半径即可，过大反而可能过完备（§10）；
- **分子、团簇、弱相互作用、表面外的电子密度**：选大一号的半径，让轨道尾巴更长；
- **周期表左侧的元素（Li、Na、K、Cs 等）**：手册指出它们需要较大的截断半径；**右侧的元素（F、Cl 等）**则需要更高角动量的极化函数；
- **团簇计算的盒子**：盒子边长至少要大于"分子尺寸 + 2 r_c"，否则相邻镜像的基函数会重叠。这个尺寸和截断能无关，只和 PAO 半径有关。

### 5.4 推荐流程

1. 先用 Standard 把结构、SCF 参数（第5章）调通；
2. 大规模扫描可以降到 Quick，但最后几个关键结构必须回到 Standard；
3. 对要发表的关键数值（能量差、带隙、吸附能），用 Precise 复核一次，差别小于你关心的精度才算收敛（§9）；
4. 所有要相减的能量（反应能、形成能、吸附能），**每个元素在每个体系里都必须用同一套 PAO + VPS**。

---

## 6. PAO 与 VPS 必须来自同一版本数据库

OpenMX 历史上发布过多个数据库版本（如 DFT_DATA13 与 DFT_DATA19）。手册第12.6节特别提醒：**有些元素在不同版本中当作价电子的半芯态数目不同**，PAO 是针对同一版本的 VPS 生成和优化的。用 DFT_DATA13 的 PAO 配 `*_PBE19` 的 VPS，半芯态的描述会对不上，结果不可信，而且程序未必会报错。

正确做法：

- 一个计算只用一个 `DATA.PATH`，所有 PAO 和 VPS 都从这一个目录取；
- 需要保留旧版本数据时，放在不同目录（如 `DFT_DATA13/`、`DFT_DATA19/`），用 `DATA.PATH` 切换；
- 在输入文件开头用注释写明数据库版本，便于以后复现。

> [!WARNING]
> 网上流传的旧输入常写 `C_PBE13`、`Fe_PBE13H` 等 2013 版名称。拿来用之前，把 PAO 和 VPS **一起**换成 DFT_DATA19 中对应的名称，并重新核对价电子数和坐标行的初始电子数。

---

## 7. 同一元素定义多个 species

`Definition.of.Atomic.Species` 的第一列只是标签，同一元素可以定义多个 species，配不同的基组甚至不同的赝势：

```text
Species.Number   3
<Definition.of.Atomic.Species
  Pt_s   Pt7.0-s3p2d2f1   Pt_PBE19     # surface layers near the adsorbate: Standard
  Pt_b   Pt7.0-s3p2d2     Pt_PBE19     # deeper bulk-like layers: Quick
  O      O6.0-s2p2d1      O_PBE19
Definition.of.Atomic.Species>
```

坐标块里给表层原子写 `Pt_s`、底层原子写 `Pt_b`，就能在保证吸附区精度的同时降低代价。常见用途：

| 用途 | 做法 | 相关章节 |
|---|---|---|
| 大体系分区精度 | 活性区 Standard/Precise，远处 Quick | 本节 |
| 反铁磁初始构型 | 同元素两个 species（如 `Ni1`、`Ni2`），初始自旋相反 | [第11章](11-spin-noncollinear.md) |
| 只给部分原子加 U | 不同 species 设不同 `Hubbard.U.values` | [第13章](13-dft-u.md) |
| 芯空穴 / XPS | 某个原子用芯空穴赝势 | 手册第13.3节、第60节 |
| 空原子 | 第3列写 `E` | §8 |

> [!CAUTION]
> 分区基组得到的总能不能直接和"全体系统一基组"的总能相减。做能量差时，比较的两个体系必须用完全相同的 species 定义；吸附能计算中，干净表面、吸附构型、孤立分子三者的基组设置要一致。

---

## 8. 空原子 E：BSSE 修正、空位与表面浮动态

### 8.1 什么是空原子

数据库里有一个特殊赝势 `E.vps`：势为零、价电子数为零。把任意 PAO 配给它，就得到一个"只有基函数、没有原子核和电子"的空原子（empty atom，也叫 ghost atom）：

```text
<Definition.of.Atomic.Species
  O    O6.0-s2p2d1   O_PBE19
  H    H6.0-s2p1     H_PBE19
  EO   O6.0-s2p2d1   E          # ghost O: basis functions only
  EH   H6.0-s2p1     E          # ghost H
Definition.of.Atomic.Species>
```

坐标行里空原子的初始电子数写 `0.0 0.0`。空原子能把基函数放到空间任意位置，常见用途有三类：

1. **counterpoise（CP）修正 BSSE**：本节重点；
2. **空位与间隙态**：在空位处放空原子，描述被俘获的电子（手册示例 `NaCl_FC.dat`：NaCl 中 Cl 空位的 F 色心，配合 `scf.system.charge -1`）；
3. **表面外的浮动态（floating states）、像势态**：在金属表面外的真空区放几层空原子。手册推荐表对 `E` 给出的是大半径的 `Kr10.0` 系列 PAO，正是为了覆盖这类弥散态。

### 8.2 BSSE 与 counterpoise 修正

基组随原子移动。在二聚体 AB 中，A 可以"借用" B 的基函数降低自己的能量，于是直接相减得到的结合能偏大（过度结合），这就是基组重叠误差（BSSE）。Boys–Bernardi 的 CP 方法让单体也在完整二聚体基组下计算：

$$
\Delta E_{\text{raw}} = E_{AB}^{AB} - E_{A}^{A} - E_{B}^{B},
\qquad
\Delta E_{\text{CP}} = E_{AB}^{AB} - E_{A}^{AB} - E_{B}^{AB},
$$

上标表示所用基组，$E_A^{AB}$ 就是"A 的真实原子 + B 位置上的空原子"。BSSE 的大小为 $\Delta E_{\text{CP}}-\Delta E_{\text{raw}}$。Ozaki 的 OpenMX 实践讲义给出的经验值是：用数据库基组时，分子体系的 BSSE 约为 0.5 kcal/mol 量级。对氢键、范德华复合物、物理吸附，这个量级不能忽略。

> [!NOTE]
> 这里单体用的是二聚体中的几何（刚性单体）。如果单体先单独优化过，还要加上形变能（deformation energy），即单体在二聚体几何与自身平衡几何之间的能量差。

### 8.3 实战：水二聚体的 CP 修正

完整输入在 [`examples/04-bsse-ghost-h2o-dimer`](../examples/04-bsse-ghost-h2o-dimer/README.md)，一共 5 个计算：

| 文件 | 体系 | 基组 | 对应能量 |
|---|---|---|---|
| `H2O_dimer.dat` | 两个水分子 | AB | $E_{AB}^{AB}$ |
| `H2O_dimer_BSSE_ghost.dat` | 水 A + 空原子 B | AB | $E_{A}^{AB}$ |
| `H2O_dimer_BSSE_ghost2.dat` | 空原子 A + 水 B | AB | $E_{B}^{AB}$ |
| `H2O_monoA.dat` | 只有水 A | A | $E_{A}^{A}$ |
| `H2O_monoB.dat` | 只有水 B | B | $E_{B}^{B}$ |

`H2O_dimer_BSSE_ghost.dat` 的完整内容如下（原子 1–3 为真实的水 A，原子 4–6 为水 B 位置上的空原子）：

```text
#
# Water dimer counterpoise (BSSE) set: E_A(AB): real monomer A + ghost atoms of B
# Rigid monomers (r_OH = 0.9572 Ang, HOH = 104.52 deg), linear H-bond, R(O...O) = 2.91 Ang.
# Monomer A = atoms 1-3 (donor), monomer B = atoms 4-6 (acceptor).
# All five inputs share the same box, cutoff, basis and SCF settings.
# OpenMX 4.0 (+patch 4.0.1), database DFT_DATA19
#
System.CurrrentDirectory    ./          # official spelling has three r's
System.Name                 h2o_A_ghostB
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            0

Species.Number              4
<Definition.of.Atomic.Species
  O    O6.0-s2p2d1   O_PBE19
  H    H6.0-s2p1     H_PBE19
  EO   O6.0-s2p2d1   E          # ghost O: basis functions only, zero potential
  EH   H6.0-s2p1     E          # ghost H
Definition.of.Atomic.Species>

Atoms.Number                6
Atoms.SpeciesAndCoordinates.Unit   Ang
<Atoms.SpeciesAndCoordinates
  1   O     4.31910   6.00935   6.00000   3.0  3.0
  2   H     5.27630   6.00935   6.00000   0.5  0.5
  3   H     4.07911   6.93598   6.00000   0.5  0.5
  4   EO    7.22910   6.00935   6.00000   0.0  0.0
  5   EH    7.54819   5.51799   6.75695   0.0  0.0
  6   EH    7.54819   5.51799   5.24305   0.0  0.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  12.0   0.0   0.0
   0.0  12.0   0.0
   0.0   0.0  12.0
Atoms.UnitVectors>

scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0       # K
scf.energycutoff            250.0       # Ry
scf.EigenvalueSolver        Cluster     # isolated molecule, no k points
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.01
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.maxIter                 100
scf.criterion               1.0e-9      # Hartree; binding energies are only a few kcal/mol

MD.Type                     Nomd
```

逐项说明：

| 关键词 | 取值 | 说明 |
|---|---|---|
| `EO`/`EH` 的 VPS | `E` | 空原子，零势、零电子 |
| 空原子 PAO | 与真实原子相同 | CP 修正要求"同一套基函数"，所以空原子用和真实 O、H 完全一样的 PAO |
| 空原子初始电子 | `0.0 0.0` | 空原子不贡献电子；体系总电子数 = 6 + 1 + 1 = 8 |
| 盒子 | 12 Å 立方 | 二聚体在 x 方向约 3.5 Å，加上 2 × 3.2 Å（6.0 bohr）的 PAO 半径仍小于 12 Å |
| `scf.EigenvalueSolver` | `Cluster` | 孤立分子用团簇模式，不需要 k 点 |
| `scf.criterion` | `1e-9` | 结合能只有几 kcal/mol，SCF 要收得紧一些 |
| 所有 5 个文件 | 同一盒子、同一截断能、同一基组 | 这样网格误差才能最大程度抵消 |

运行与汇总（`run_bsse.sh` 依次跑 5 个输入，并用 awk 计算 ΔE_raw、ΔE_CP 和 BSSE，单位 kcal/mol）：

```bash
cd OpenMX/examples/04-bsse-ghost-h2o-dimer
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
MPI_RUN="mpirun -np 4" bash run_bsse.sh
grep "Utot\." *.out                  # total energies in Hartree
```

判断标准：

- 5 个计算都要 SCF 收敛（`.out` 中 SCF history 的最后一个 `NormRD` 很小，见第5章）；
- `h2o_A_ghostB.out` 中 Mulliken 布居的总电子数应为 8，空原子上会有少量布居（它们承载了 A 借用的那部分电子密度），这是正常的；
- ΔE_CP 应比 ΔE_raw 弱（绝对值更小），两者之差就是 BSSE；
- 作为量级参考，高精度 CCSD(T) 基准（S22 数据集）给出的水二聚体结合能约 −5.0 kcal/mol。本例的几何是手工搭建的刚性结构，并非 PBE 极小，数值只用于对比 ΔE_raw 与 ΔE_CP，不应期望与基准逐位吻合。

> [!TIP]
> 单体只放在自己的基组里计算时（`H2O_monoA.dat`），原子相对积分网格的位置和二聚体计算不同，会引入额外的网格误差（egg-box effect）。要做到 0.01 kcal/mol 级别的比较，可以先从二聚体的标准输出（`.std`）中读出 `Grid_Origin`，再在单体输入里用 `scf.fixed.grid` 固定成同一个原点（第5章 §4.4）。

---

## 9. 基组收敛测试：Si 的 Quick → Standard → Precise

### 9.1 测试思路

平面波程序只需加大 ENCUT。OpenMX 的基组收敛要同时看三类性质，而且测试时其他参数（截断能、k 点、几何）必须固定：

| 性质 | 怎么比 | 合格的参考标准 |
|---|---|---|
| 总能 | 相同结构下比较 Utot/atom | 能量差（而非绝对值）变化 < 你关心的精度，常用 1–5 meV/atom |
| 晶格常数、体模量 | 各基组做 E–V 拟合 | a0 变化 < 0.2–0.5%，B0 变化 < 2–3% |
| 能带 | 同一 k 路径比较带隙与色散 | 价带色散基本重合，带隙变化 < 0.05 eV |

> [!NOTE]
> LCPAO 对基组是变分的：在**同一积分网格**下，基组越大，总能越低。如果 Precise 的能量反而高于 Standard，先检查两次计算的网格、k 点和结构是否真的相同，再考虑过完备（§10）。

### 9.2 输入文件

[`examples/04-basis-test-si`](../examples/04-basis-test-si/README.md) 中有三个只在基组一行不同的输入。以 Standard 为例：

```text
#
# Si (diamond) basis set test: STANDARD  (Si7.0-s2p2d1, 13 functions/atom)
# OpenMX 4.0 (+patch 4.0.1), database DFT_DATA19
# Only the basis line and System.Name differ between the quick/standard/precise files.
#

#
# File names
#
System.CurrrentDirectory    ./          # official spelling has three r's
System.Name                 si_standard
DATA.PATH                   /opt/openmx4.0/DFT_DATA19   # EDIT: directory that contains PAO/ and VPS/
level.of.stdout             1
level.of.fileout            0           # 0: no cube files (energies and bands only)

#
# Atomic species: name  PAO  VPS
#
Species.Number              1
<Definition.of.Atomic.Species
  Si   Si7.0-s2p2d1   Si_PBE19
Definition.of.Atomic.Species>

#
# Structure: fcc primitive cell, a = 5.43 Ang (lattice vector components = a/2)
#
Atoms.Number                2
Atoms.SpeciesAndCoordinates.Unit   FRAC
<Atoms.SpeciesAndCoordinates
  1   Si   0.00   0.00   0.00   2.0  2.0
  2   Si   0.25   0.25   0.25   2.0  2.0
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit      Ang
<Atoms.UnitVectors
  0.0000  2.7150  2.7150
  2.7150  0.0000  2.7150
  2.7150  2.7150  0.0000
Atoms.UnitVectors>

#
# SCF (kept identical for all three basis sets)
#
scf.XcType                  GGA-PBE
scf.SpinPolarization        off
scf.ElectronicTemperature   300.0       # K
scf.energycutoff            220.0       # Ry, real-space grid (not a basis cutoff)
scf.EigenvalueSolver        Band
scf.Kgrid                   8 8 8
scf.Generation.Kpoint       GCenter     # OpenMX 4.0 default; remove this line for OpenMX 3.9
scf.Mixing.Type             RMM-DIISK
scf.Init.Mixing.Weight      0.10
scf.Min.Mixing.Weight       0.01
scf.Max.Mixing.Weight       0.30
scf.Mixing.History          20
scf.Mixing.StartPulay       8
scf.maxIter                 100
scf.criterion               1.0e-8      # Hartree

MD.Type                     Nomd

#
# Band structure along L-G-X-W-K-G (conventional cubic cell used for the k-path)
#
Band.dispersion             on
<Band.KPath.UnitCell
  5.43  0.00  0.00
  0.00  5.43  0.00
  0.00  0.00  5.43
Band.KPath.UnitCell>
Band.Nkpath                 5
<Band.kpath
  40   0.5  0.5  0.5    0.0  0.0  0.0    L  G
  40   0.0  0.0  0.0    1.0  0.0  0.0    G  X
  20   1.0  0.0  0.0    1.0  0.5  0.0    X  W
  15   1.0  0.5  0.0    0.75 0.75 0.0    W  K
  40   0.75 0.75 0.0    0.0  0.0  0.0    K  G
Band.kpath>
```

另外两个文件只改两行：

```text
# Si_basis_quick.dat
System.Name                 si_quick
  Si   Si7.0-s2p1d1   Si_PBE19          # 10 functions/atom

# Si_basis_precise.dat
System.Name                 si_precise
  Si   Si7.0-s3p3d2   Si_PBE19          # 22 functions/atom
```

关键设置的理由：

| 关键词 | 取值 | 为什么 |
|---|---|---|
| `Atoms.SpeciesAndCoordinates.Unit` | `FRAC` | 改晶格常数时只需改 `Atoms.UnitVectors`，原子跟着缩放 |
| `scf.energycutoff` | 220 Ry | 三个基组用同一个值；截断能本身的收敛在第5章测试 |
| `scf.Kgrid` | 8 8 8 | 对 Si 原胞足够（第5章会验证）；三次计算完全一致 |
| `<Band.KPath.UnitCell` | 惯用立方胞 5.43 Å | 用惯用胞的倒格矢写 fcc 高对称点 L、X、W、K，坐标才和教科书一致（见 [第8章](08-band.md)） |
| `level.of.fileout` | 0 | 只比较能量和能带，不输出 cube，节省磁盘 |

### 9.3 运行与比较总能

```bash
cd OpenMX/examples/04-basis-test-si
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once
MPI_RUN="mpirun -np 4" bash run_basis_test.sh scf
```

脚本依次运行三个输入，然后从 `.out` 中取 `Utot.`，输出每原子能量及相对 Precise 的差（meV/atom）。手动检查可以用：

```bash
grep "Utot\." si_quick.out si_standard.out si_precise.out
grep -A3 "Required cutoff energy" si_standard.out    # same grid for all three?
```

下面是结果表的格式（**示意输出，数值仅供参考**）：

```text
# basis      nbasis/atom   Utot(Ha/cell)       dE vs precise (meV/atom)
quick           10         -7.8xxxxxxxxx           +xx.x
standard        13         -7.8xxxxxxxxx            +x.x
precise         22         -7.8xxxxxxxxx             0.0
```

预期的趋势是 Quick > Standard > Precise（基组越大能量越低）。作为参照，数据库 Si 页面给出：在 E–V 曲线的极小点，`Si7.0-s2p2d1` 与 `Si7.0-s3p3d2f1` 的总能相差 0.025 eV/atom。绝对总能的这点差别对大多数**能量差**影响很小，所以 Standard 是性价比最高的选择。

### 9.4 比较晶格常数：用 EvsLC 扫描

总能本身不是物理量，真正要看的是基组对平衡结构的影响。`run_basis_test.sh eos` 会把三个输入改写成 `MD.Type EvsLC` 扫描（从 a = 5.24 Å 开始，每步放大 0.7%，共 12 步，并用 `scf.Ngrid` 固定实空间网格），再用第6章示例中的拟合脚本求 a0 和 B0：

```bash
MPI_RUN="mpirun -np 4" bash run_basis_test.sh eos
for b in quick standard precise; do
  python ../06-cellopt-eos-si-c/fit_evslc.py si_${b}_eos.EvsLC --lattice fcc --natoms 2
done
```

`EvsLC` 的原理、为什么必须固定 `scf.Ngrid`，在 [第6章 §5](06-optimization-eos-neb.md) 详细介绍。判断：Standard 与 Precise 的 a0 差应在 0.01 Å 量级；PBE 下 Si 的晶格常数文献值约 5.47 Å（实验 5.431 Å），GGA 略微高估属正常。

### 9.5 比较能带

三个计算都打开了 `Band.dispersion`，会生成 `si_quick.Band` 等文件。用 bandgnu13 或本仓库的 `omx_band_plot.py`（`python ../scripts/omx_band_plot.py -h` 查看用法）画图并读出带隙，详见 [第8章](08-band.md)。一般会看到：

- 价带三个基组几乎重合；
- 导带对基组更敏感，Quick 的导带可能偏高；
- 数据库 Si 页面的能带对比显示，`Si7.0-s2p2d1` 已接近 `Si7.0-s3p3d2f1` 和全电子结果。

---

## 10. 基组过完备：识别与处理

### 10.1 现象

官方手册"Known issues"一节写明：在 fcc、hcp、bcc 等**高密度体相**中，如果给每个原子分配很多基函数（特别是大截断半径 + 多 ζ），基组可能变成过完备（overcomplete），重叠矩阵接近奇异，数值不稳定，更容易出现**伪本征值（spurious eigenvalues）**。常见症状：

- 能带图或本征值列表里出现远低于价带、不该存在的平带，或能量异常低的态；
- SCF 难以收敛，`NormRD` 长时间振荡，换混合方案也没用；
- 换更大的基组后总能"降低"得不合理，晶格常数或磁矩反而变差；
- E–V 曲线在某些晶格常数处出现跳变（也可能是网格问题，先按第6章固定 `scf.Ngrid` 排除）。

### 10.2 检查方法

```bash
# lowest Kohn-Sham eigenvalues; compare Quick / Standard / Precise
grep -A12 "Eigenvalues (Hartree)" fe_precise.out | head -20
# is the SCF really converged? (scf_conv=1 in every block of the .md file)
grep -o "scf_conv=[01]" fe_precise.md | sort | uniq -c
```

把不同基组得到的最低几个本征值放在一起看：半芯态（例如 Fe 3s、3p）的能量应基本不随基组变化。如果大基组多出一个能量很低、在小基组里不存在的态，就要怀疑过完备。

### 10.3 处理办法（按手册建议）

1. **换小一些的优化基组**：例如致密 bcc Fe 用推荐表里的 `Fe5.5H-s3p2d1`（推荐表对 `Fe_PBE19H` 三档都用 5.5 bohr 的 PAO），不要上来就用 `Fe8.0H-s3p3d2f1` 这类大半径、多 ζ 的组合；
2. **关闭中性原子势的投影展开**：

```text
scf.ProExpn.VNA      off      # default = on
scf.energycutoff     300.0    # raise the grid cutoff when VNA is evaluated on the grid
```

`scf.ProExpn.VNA off` 后，中性原子势 VNA 的矩阵元改在实空间网格上计算，通常要相应提高 `scf.energycutoff`（第5章）。数据库网站上的基准输入（如 Si 的 E–V 测试）本身就设置了 `scf.ProExpn.VNA off`，可以作为参考。

---

## 11. 生成输入时的基组设置工具

- **本仓库** `omx_make_input.py`：从 POSCAR/CIF 生成 `.dat`，按 Quick/Standard/Precise 自动填基组和价电子数（`python ../scripts/omx_make_input.py -h`）；
- **OpenMX Viewer**（浏览器端，2026-05 大更新）：读入 CIF/POSCAR 后导出 OpenMX 输入，可以在界面里选基组；
- **cif2omx**（4.0 起手册第8.3节正式收录）：从 CIF 生成 `.dat`；
- 结构准备与格式转换的通用方法见 Common 第3章，ASE 的 OpenMX 接口见 Common 第4章。

> [!WARNING]
> 自动工具写出的基组只是起点。生成后务必打开 `.dat` 核对三件事：PAO 与 VPS 的版本一致、坐标行上下自旋电子数之和等于 `valence.electron`、`DATA.PATH` 指向正确的 DFT_DATA19。

---

## 参数选择与收敛建议

| 项目 | 起步推荐 | 何时升级 | 验证方式 |
|---|---|---|---|
| 基组档次 | Standard | 弱相互作用、能量差 < 10 meV、精细能带 | Standard 与 Precise 比较（§9） |
| PAO 截断半径 | 推荐表的值 | 分子、表面、真空态：大一号；致密金属：别加大 | 数据库元素页面的基准 |
| VPS（H/S） | 推荐表任选其一，常用 H | 高压、强磁、短键、需要半芯极化 | 比较晶格常数与磁矩 |
| 泛函与 VPS | GGA-PBE + `*_PBE19` | 只在需要与 LDA 文献对比时用 `*_CA19` | 检查 `xc.type` |
| 能量差计算 | 所有体系同一 species 定义 | — | 输入 diff 检查 |
| 分子间作用、吸附 | CP 修正 | BSSE > 你要的精度 | §8 |
| 过完备 | 推荐表小半径 | 出现伪本征值 | `scf.ProExpn.VNA off` + 提高截断能 |

---

## 常见问题与报错

**1. `Found an undefined species name ... in Atoms.SpeciesAndCoordinates or Hubbard.U.values` / `Please check your input file`**

- 原因 A：坐标块或 `Hubbard.U.values` 里用的 species 名在 `Definition.of.Atomic.Species` 中没有定义（大小写、拼写不一致，或 `Species.Number` 写小了）。
- 原因 B：输入文件在 Windows 编辑器里存成了 UTF-16、带 BOM、含全角空格或 CRLF 换行，程序读到的名字里混入了不可见字符。
- 解决：核对拼写；用 `file Si.dat` 检查编码，必要时 `iconv -f UTF-16 -t ASCII Si.dat > Si_fixed.dat`，再 `dos2unix Si_fixed.dat`。

**2. 找不到 PAO/VPS 文件（程序在读基组阶段退出，标准输出提示无法打开对应文件）**

- 原因：`DATA.PATH` 是默认的 `../DFT_DATA19`，而你不在 `work/` 下运行；或者路径写到了 `PAO/` 子目录；或者 PAO/VPS 名称拼错（如 `Si7-s2p2d1` 少写了 `.0`，`Fe5.5-s3p2d1` 漏了 `H`）。
- 解决：`DATA.PATH` 写绝对路径并指向 `DFT_DATA19` 本身；用 `ls $DATA/PAO/Si7.0.pao $DATA/VPS/Si_PBE19.vps` 确认文件存在。

**3. Mulliken 总电子数与预期不符，或 SCF 一开始就发散**

- 原因：坐标行上下自旋初始电子数之和不等于 VPS 价电子数（例如 Fe 写成 `4.0 4.0`，以为价电子是 8）。
- 解决：`grep valence.electron` 查 VPS；带电体系用 `scf.system.charge`，不要靠改初始占据数。

**4. 基组写法报错或基函数个数不对**

- 原因：要求的角动量或径向函数个数超出了 `.pao` 文件的内容（DFT_DATA19 每个文件 l ≤ 3，每个 l 15 个径向函数）；或者 VPS 的半芯态没有对应基函数（例如 `Fe5.5H-s1p1d1`，少了一个 s）。
- 解决：用 §3.4 的 `grep PAO.Lmax` 查看文件；半芯态 VPS 按数据库给的最小基组以上设置。

**5. 换了大基组后能量、磁矩或能带反而异常**

- 原因：致密体系基组过完备（§10）。
- 解决：回到推荐表的半径；或 `scf.ProExpn.VNA off` 并提高 `scf.energycutoff`。

**6. 复现文献或旧输入时结果对不上**

- 原因：数据库版本不同（DFT_DATA13 与 DFT_DATA19 的半芯态不同），Sn 推荐基组在 4.0 中更改，或者 4.0 默认的 k 点生成方式不同（第5章 GCenter 与 Regular）。
- 解决：确认对方使用的 PAO/VPS 全名与版本，统一后再比较。

---

## 练习

1. 写出 `Mo7.0-s3p2d2f1`、`Fe5.5H-s3p2d2f1`、`Sn7.0-s3p2d3` 每原子的基函数个数，再估算 MoS2 单层 3×3 超胞（9 个 Mo、18 个 S，用推荐表 Standard）的基函数总数。
2. 对 bcc Fe 分别用 `Fe_PBE19H`（`Fe5.5H-s3p2d1`）和 `Fe_PBE19S`（`Fe6.0S-s3p2d1`）计算，比较磁矩和计算时间。坐标行的初始电子数要怎么改？
3. 在水二聚体例子中把 H、O 的 PAO 换成 Precise（`H7.0-s2p2d1`、`O6.0-s3p2d2`），BSSE 变大还是变小？为什么？
4. 用 `run_basis_test.sh eos` 求三档基组的 a0，画出 a0 随每原子基函数个数的变化。
5. 把 Si 的 PAO 换成 `Si9.0-s3p3d2`，观察总能、SCF 步数和最低本征值，判断是否出现过完备迹象。

---

## 参考资料

**官方文档**

- [OpenMX 4.0 手册 第12章 Basis functions](https://www.openmx-square.org/openmx_man4.0/s12_basis.html) — PAO 概念、原始与优化轨道、空原子、`DATA.PATH`
- [手册 12.3 Database Ver. 2019 推荐基组表](https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html) — Quick/Standard/Precise 与价电子数
- [手册 12.5 Placement of empty atoms](https://www.openmx-square.org/openmx_man4.0/s12_5_placement-empty-atoms.html) — 空原子、CP 方法、F 色心示例
- [手册 第13章 Pseudopotentials](https://www.openmx-square.org/openmx_man4.0/s13_pseudo.html) — MBK 赝势、价电子数、open-core 赝势
- [手册 38.3 Scalar-relativistic treatment](https://www.openmx-square.org/openmx_man4.0/s38_3_scalar-rel.html) — j 依赖赝势在共线计算中的自动平均
- [手册 第26章 Optimization of basis functions](https://www.openmx-square.org/openmx_man4.0/s26_orb-opt.html) — 收缩写法 `s4>1` 与轨道优化
- [手册 第79章 Known issues](https://www.openmx-square.org/openmx_man4.0/s79_known-issues.html) — 基组过完备
- [Database (2019) of optimized VPS and PAO](https://www.openmx-square.org/vps_pao2019/) — 各元素页面的基准计算与 Δ 因子
- [What's new（4.0 的 Sn 推荐基组变更、Tb/Yb 数据）](https://www.openmx-square.org/whatsnew/whatsnew.html)
- [Practical guide to OpenMX calculations（Ozaki 2023 讲义，PDF）](https://t-ozaki.issp.u-tokyo.ac.jp/lectures-2023/OpenMX_Practical_Guide.pdf) — 基组选择与 BSSE 的经验值

**相关 GitHub 仓库（仅链接）**

- [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx) — POSCAR 转 `.dat`，可选 Quick/Standard/Precise 基组
- [issp-center-dev/cif2x](https://github.com/issp-center-dev/cif2x) — CIF + 模板批量生成 OpenMX 输入
- [mitsuaki1987/cif2input](https://github.com/mitsuaki1987/cif2input) — CIF 转 OpenMX 输入，可与官方 cif2omx 对照
- [by-student-2017/cif2cell-informal](https://github.com/by-student-2017/cif2cell-informal) — cif2cell 非官方版，支持输出 OpenMX 格式并指定 PAO 目录
- [wxia529/acoord](https://github.com/wxia529/acoord) — VS Code 中编辑 `.dat` 结构并保留参数块
- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu) — 4.0 GPU 版，扩展了部分镧系元素的 VPS/PAO
- [pranabdas/openmx](https://github.com/pranabdas/openmx) — 英文 OpenMX 教程站（无许可证，仅供参考）
- 更多见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**

- T. Ozaki, *Phys. Rev. B* 67, 155108 (2003)（变分优化原子轨道）. DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)
- T. Ozaki, H. Kino, *Phys. Rev. B* 69, 195113 (2004)（H–Kr 数值原子基组）. DOI: [10.1103/PhysRevB.69.195113](https://doi.org/10.1103/PhysRevB.69.195113)
- T. Ozaki, H. Kino, *Phys. Rev. B* 72, 045121 (2005)（中性原子势的投影展开 `scf.ProExpn.VNA`）. DOI: [10.1103/PhysRevB.72.045121](https://doi.org/10.1103/PhysRevB.72.045121)
- I. Morrison, D. M. Bylander, L. Kleinman, *Phys. Rev. B* 47, 6728 (1993)（MBK 赝势）. DOI: [10.1103/PhysRevB.47.6728](https://doi.org/10.1103/PhysRevB.47.6728)
- S. F. Boys, F. Bernardi, *Mol. Phys.* 19, 553 (1970)（counterpoise 方法）. DOI: [10.1080/00268977000101561](https://doi.org/10.1080/00268977000101561)
- K. Lejaeghere et al., *Science* 351, aad3000 (2016)（Δ 因子与代码间精度比较）. DOI: [10.1126/science.aad3000](https://doi.org/10.1126/science.aad3000)
- P. Jurečka, J. Šponer, J. Černý, P. Hobza, *Phys. Chem. Chem. Phys.* 8, 1985 (2006)（S22 非共价相互作用基准）. DOI: [10.1039/B600027D](https://doi.org/10.1039/B600027D)

---

[← 上一章](03-input-first-calc.md) | [返回目录](../README.md) | [下一章 →](05-scf-kpoints-convergence.md)
