# 第12章 声子：phonopy 有限位移与 DFPT（Phonons with VASP and phonopy: Finite Displacement and DFPT）

本章系统讲解用 VASP + phonopy 计算声子谱、声子态密度和热力学性质：有限位移法（finite displacement / frozen phonon）的完整流程、VASP 密度泛函微扰理论（DFPT，IBRION=7/8）+ phonopy 的力常数路线、极性材料的 Born 有效电荷与非解析项修正（NAC，LO–TO 劈裂）、band.conf/mesh.conf 的写法与虚频诊断，以及 VASP 6.x 原生声子功能（LPHON_DISPERSION、VASP 6.6 的 CHECKPOINT_FD 断点续算）、准谐近似（QHA）热膨胀和 phono3py 晶格热导率入门。phonopy 部分同时给出 4.x（phonopy-init）新写法与 2.x/3.x 旧写法。

---

## 1. 本章目标与适用场景

| 目标 | 方法 | 本章小节 | 配套示例 |
|---|---|---|---|
| 半导体/金属的声子谱、DOS、热力学 | phonopy 有限位移 | §4–§7 | [12-si-phonopy-fd](../examples/12-si-phonopy-fd/README.md) |
| 用 DFPT 一次算出超胞力常数 | VASP IBRION=8 + phonopy | §8 | [12-mgo-dfpt-nac](../examples/12-mgo-dfpt-nac/README.md) |
| 极性材料的 LO–TO 劈裂 | LEPSILON + BORN + NAC | §9 | 同上 |
| 判断结构动力学稳定性、处理虚频 | 虚频诊断 | §10 | — |
| 不用 phonopy，直接在 VASP 中出声子谱；大体系拆分/续算 | IBRION=6 + LPHON_DISPERSION + CHECKPOINT_FD | §11 | — |
| 二维材料声子 | 有限位移 + 2D 注意事项 | §12 | [12-graphene-phonon](../examples/12-graphene-phonon/README.md) |
| 热膨胀、晶格热导率 | QHA、phono3py | §13 | — |

**前置知识**：结构优化与 Pulay 应力（[第4章](04-relax-eos.md)）、收敛测试（[第3章](03-convergence.md)）、介电性质（[第14章](14-optics-dielectric.md)）。弹性常数（同样基于有限差分 IBRION=6）在[第13章](13-elastic.md)。吸附物频率（Γ 点、Selective dynamics）已在[第10章](10-adsorption-catalysis.md)介绍。

> [!NOTE]
> 版本：VASP 6.4–6.6（推荐 6.6.1）；phonopy 4.10.0（2026-11）为写作时最新版。phonopy 4.0（2026-05）做了命令行拆分等破坏性更新，本章对新旧写法分别标注。文中输出片段除注明来源外均为**示意输出，数值仅供参考**。

---

## 2. 原理速览

### 2.1 简谐近似与力常数

把原子位移 $u$ 展开到二阶，二阶力常数为

$$
\Phi_{\alpha\beta}(l\kappa,l'\kappa')=\frac{\partial^2E}{\partial u_\alpha(l\kappa)\,\partial u_\beta(l'\kappa')}
\approx-\frac{F_\beta(l'\kappa')}{u_\alpha(l\kappa)}
$$

后一式就是**有限位移法**：在超胞中让一个原子位移 $u$（phonopy 默认 0.01 Å），计算所有原子受力 $F$，由"力/位移"得到力常数（实际用正负位移和对称性求解）。

动力学矩阵与声子频率：

$$
D_{\alpha\beta}(\kappa\kappa',\mathbf q)=\sum_{l'}\frac{\Phi_{\alpha\beta}(0\kappa,l'\kappa')}{\sqrt{m_\kappa m_{\kappa'}}}e^{i\mathbf q\cdot[\mathbf r(l'\kappa')-\mathbf r(0\kappa)]},\qquad
D(\mathbf q)\,\mathbf e=\omega^2\mathbf e
$$

- $\omega^2<0$ 对应**虚频**（phonopy/VASP 中以负频率或 `f/i` 显示）；
- 超胞必须足够大，使力常数在超胞边界衰减到可忽略，否则傅里叶插值得到的非 Γ 点频率不准；
- 声学求和规则（ASR）要求 Γ 点三个声学支频率为 0，数值误差会使其偏离。

### 2.2 DFPT

DFPT 通过求解线性 Sternheimer 方程直接得到能量对原子位移的二阶导数。VASP 的实现（IBRION=7/8）只计算与超胞公度的 **q = 0 微扰**，因此同样需要超胞；官方文档指出它仅支持 LDA/GGA，不计算刚性离子弹性张量，相对有限差分的主要优点是不必选择位移幅度。

### 2.3 极性材料与非解析项

在离子晶体中，长程偶极–偶极相互作用使 $\mathbf q\to0$ 时纵光学（LO）模与横光学（TO）模劈裂。非解析项修正需要高频介电张量 $\varepsilon^\infty$ 与 Born 有效电荷 $Z^*$：

$$
D^\mathrm{NA}_{\alpha\beta}(\kappa\kappa',\mathbf q\to0)=\frac{4\pi}{\Omega\sqrt{m_\kappa m_{\kappa'}}}
\frac{(\mathbf q\cdot Z^*_\kappa)_\alpha(\mathbf q\cdot Z^*_{\kappa'})_\beta}{\mathbf q\cdot\varepsilon^\infty\cdot\mathbf q}
$$

phonopy 默认使用 Gonze–Lee 方法（`NAC_METHOD = GONZE`），也可选 Wang 等的方法。

### 2.4 热力学性质

由声子态密度可得（每个原胞）：

$$
F=\frac12\sum_{\mathbf q\nu}\hbar\omega_{\mathbf q\nu}+k_BT\sum_{\mathbf q\nu}\ln\left[1-e^{-\hbar\omega_{\mathbf q\nu}/k_BT}\right],\qquad
C_V=\sum_{\mathbf q\nu}k_B\left(\frac{\hbar\omega}{k_BT}\right)^2\frac{e^{\hbar\omega/k_BT}}{(e^{\hbar\omega/k_BT}-1)^2}
$$

以及熵 $S=-\partial F/\partial T$。在多个体积下重复，就得到准谐近似（QHA）下的 $G(T,p)$、热膨胀系数和 $C_p$。

---

## 3. 工具与版本说明

### 3.1 phonopy 新旧写法对照（重要）

| 步骤 | phonopy ≥ 4.0（当前 4.10） | phonopy 2.x / 3.x |
|---|---|---|
| 生成位移超胞 | `phonopy-init -d --dim 2 2 2 -c POSCAR` | `phonopy -d --dim 2 2 2 -c POSCAR` |
| 收集力 → FORCE_SETS | `phonopy-init -f disp-*/vasprun.xml` | `phonopy -f disp-*/vasprun.xml` |
| DFPT 力常数 → FORCE_CONSTANTS | `phonopy-init --fc vasprun.xml` | `phonopy --fc vasprun.xml` |
| 读配置文件做后处理 | `phonopy --config band.conf -p -s` | `phonopy -p -s band.conf` |
| 开启 NAC | 有 BORN 文件即自动开启（`--nonac` 关闭） | `NAC = .TRUE.` 或 `--nac` |
| 原胞（PRIMITIVE_AXES） | 默认 `AUTO` | 默认单位矩阵（输入胞即原胞） |

> [!WARNING]
> phonopy 4.x 的 `phonopy` 命令会拒绝 `-d`、`-f`、`--fc` 等前处理选项并提示改用 `phonopy-init`。在 4.10.0 中，`phonopy band.conf` 这种把配置文件作为位置参数的旧写法会报错：`"band.conf" could not be read as a phonopy.yaml-like file. If "band.conf" is a configuration file, give it by --config.`——必须写成 `phonopy --config band.conf`。本章示例的 conf 文件同时写了 `DIM` 与显式 `PRIMITIVE_AXES`，所以在新旧版本下结果一致。

其他与版本相关的事项：

- phonopy 4.8.2 起需要 **Python ≥ 3.11**；2.46.0 起不再支持 Python 3.9。集群 Python 较旧时，可在 conda 环境中安装新版，或使用对应的旧版本 phonopy。
- `phonopy-vasp-born` 读取 vaspout.h5 时 Born 电荷被转置的 bug 存在于 4.4.0–4.7.x，**4.8.0 已修复**；默认读取 vasprun.xml 不受影响。
- `phonopy-qha` 命令计划在未来版本移除，官方建议新计算改用 Python API `run_qha`。

安装：

```bash
conda install -c conda-forge phonopy seekpath     # 或 pip install phonopy seekpath
phonopy --version
```

### 3.2 VASP 侧的声子相关功能

| 功能 | 关键词 | 版本 |
|---|---|---|
| 有限差分（全部位移 / 利用对称性） | `IBRION = 5 / 6`，`POTIM`，`NFREE` | 5.1+ |
| DFPT（全部 / 利用对称性） | `IBRION = 7 / 8` | 5.1+ |
| 介电张量与 Born 有效电荷 | `LEPSILON`（DFPT）/ `LCALCEPS`（有限电场） | — |
| VASP 内部插值声子色散 | `LPHON_DISPERSION` + QPOINTS 文件 | ≥ 6.3.2 |
| 读入已有力常数（vaspin.h5） | `LPHON_READ_FORCE_CONSTANTS` | ≥ 6.4.0 |
| 有限差分断点续算/拆分 | `CHECKPOINT_FD`（仅 IBRION=6，需 HDF5） | ≥ 6.6.0 |
| 声子频率写入 vaspout.h5 | — | 6.6 |

> [!IMPORTANT]
> VASP 自身的 IBRION=5–8 只给出**超胞的 Γ 点**频率；要得到整个布里渊区的色散，需要 phonopy（或 VASP ≥ 6.3.2 的 LPHON_DISPERSION）对超胞力常数做傅里叶插值。

### 3.3 整体流程

```text
1. 体相紧致优化 (ISIF=3, EDIFFG ~ -1E-3)  ──────────────────────────────┐
2. phonopy-init -d  →  SPOSCAR + POSCAR-001 ...                         │
3a. 有限位移：每个 POSCAR-XXX 做单点 (IBRION=-1) → phonopy-init -f      │
3b. DFPT：SPOSCAR 做 IBRION=8 → phonopy-init --fc                        │
4. (极性) 原胞 LEPSILON → phonopy-vasp-born → BORN                       │
5. phonopy --config band.conf / mesh.conf → 声子谱、DOS、热力学          │
6. 虚频诊断；(可选) QHA、phono3py  ←──────────────────────────────────────┘
```

仓库脚本 `phonopy_workflow.sh` 把上述步骤（生成位移、批量建目录与提交、收集、由 seekpath 生成 band.conf、作图与虚频检测）做成了自动化流程，用法见 `bash ../scripts/phonopy_workflow.sh -h`。下面先手动走一遍，理解每一步。

---

## 4. 第一步：高精度结构优化

声子对结构是否处于能量极小点非常敏感：残余受力或残余应力会直接表现为虚频。

```bash
SYSTEM  = Si conventional cell - tight relaxation before the phonon calculation

PREC    = Accurate
ENCUT   = 400          # ~1.6 x ENMAX(Si ~245 eV): converged stress, no Pulay error
EDIFF   = 1E-8
ALGO    = Normal
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.

IBRION  = 2
ISIF    = 3            # relax cell and positions; repeat (CONTCAR -> POSCAR) until stable
NSW     = 100
EDIFFG  = -1E-3        # phonons need a structure very close to the minimum

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ENCUT` | ≥ 1.3 × ENMAX | 变胞优化必须避免 Pulay 应力；声子计算中**沿用同一 ENCUT** |
| `EDIFF` | 1E-8 | 力常数需要 ~1E-4 eV/Å 量级的力精度 |
| `EDIFFG` | −1E-3 | 比普通优化严一个数量级 |
| `ISIF` | 3 | 优化后 `cp CONTCAR POSCAR` 再跑一次，直到体积不再变化（见[第4章](04-relax-eos.md)） |

```bash
cd VASP/examples/12-si-phonopy-fd
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR        # 仓库不提供 POTCAR
mkdir -p relax
cp INCAR_relax relax/INCAR; cp KPOINTS_relax relax/KPOINTS; cp POSCAR POTCAR relax/
(cd relax && mpirun -np 16 vasp_std > vasp.out)
grep "reached required accuracy" relax/OUTCAR
grep "external pressure" relax/OUTCAR | tail -1    # |P| 应 < ~0.5 kB
cp relax/CONTCAR POSCAR                            # 用优化后的结构生成超胞
```

> [!TIP]
> 示例 POSCAR 的晶格常数 5.469 Å 只是 PBE 的近似值（实验 5.431 Å）。务必用**你自己的** ENCUT/POTCAR 优化后的 CONTCAR 来生成位移超胞。

---

## 5. 第二步：phonopy 生成位移超胞

```bash
# phonopy >= 4.0
phonopy-init -d --dim 2 2 2 -c POSCAR
# phonopy 2.x / 3.x
# phonopy -d --dim 2 2 2 -c POSCAR
ls
# POSCAR  POSCAR-001  SPOSCAR  phonopy_disp.yaml
```

| 输出文件 | 含义 |
|---|---|
| `SPOSCAR` | 无位移的完美超胞（DFPT 用它） |
| `POSCAR-001 …` | 带位移的超胞，每个做一次单点计算；金刚石 Si 由于对称性高只需 1 个 |
| `phonopy_disp.yaml` | 晶胞、超胞矩阵、位移信息；后处理要读它（不要删） |

**超胞大小怎么选**：

- 经验上超胞各方向 ≥ 10–12 Å，复杂/软模体系更大；
- 最终要做收敛测试：2×2×2 → 3×3×3（或 `DIM = -1 1 1 1 -1 1 1 1 -1` 这类非对角超胞）比较声子谱；
- 本例输入是 8 原子惯用胞，`--dim 2 2 2` 得到 64 原子超胞（边长约 10.9 Å）。

**其他常用选项**：

| 选项 / 标签 | 说明 |
|---|---|
| `--amplitude 0.01`（`DISPLACEMENT_DISTANCE`） | 位移幅度，默认 0.01 Å |
| `--pm`（`PM = .TRUE.`） | 同时生成正负位移（默认 `AUTO`：只在正负不等价时生成负位移） |
| `--nodiag`（`DIAG = .FALSE.`） | 只沿晶轴位移，晶格常数差别很大时推荐 |
| `--pa auto` / `PRIMITIVE_AXES` | 原胞变换；v4 默认 auto |
| `MAGMOM` | 磁性体系，phonopy 会写出超胞的 MAGMOM 文件 |

---

## 6. 第三步：VASP 单点力计算

### 6.1 INCAR

```bash
SYSTEM  = Si 2x2x2 supercell (64 atoms) - single-point forces for phonopy

# ---- accuracy: forces of ~1E-4 eV/A matter for force constants ----
PREC    = Accurate
ENCUT   = 400          # exactly the value used for the relaxation
EDIFF   = 1E-8
NELMIN  = 5            # avoid a premature SCF stop
ALGO    = Normal
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.      # real-space projection adds force noise
ISYM    = 0            # optional: no symmetrisation of the displaced cell
# ADDGRID = .TRUE.     # common in phonopy inputs; the VASP wiki advises testing it first

# ---- single point: the displaced atom must NOT relax ----
IBRION  = -1
NSW     = 0

LWAVE   = .FALSE.
LCHARG  = .FALSE.
NCORE   = 4
```

| 关键词 | 推荐值 | 原因 |
|---|---|---|
| `IBRION`/`NSW` | −1 / 0 | **绝对不能**让位移超胞弛豫，否则位移被"抹平" |
| `PREC` | Accurate | 减少混叠误差，力更准 |
| `EDIFF` | 1E-8 | 位移 0.01 Å 引起的力只有 ~0.01–0.1 eV/Å，需要高 SCF 精度 |
| `NELMIN` | 5 | 防止 SCF 在前几步因能量变化小而提前结束 |
| `LREAL` | `.FALSE.` | 实空间投影引入的力噪声会产生虚频 |
| `ISYM` | 0（可选） | 关闭对称性，避免位移结构对称性判断的噪声；代价是更多不可约 k 点 |
| `ADDGRID` | 慎用 | 很多 phonopy 示例使用 `ADDGRID = .TRUE.`；VASP 官方声子页面则提醒不经测试不要设置 |
| `ISMEAR`/`SIGMA` | 半导体 0 / 0.01；金属 1 / 0.1–0.2 | 金属需要测试 SIGMA 对声子的影响 |
| `ENCUT`、POTCAR | 与优化完全相同 | 保持同一势能面 |

### 6.2 KPOINTS：随超胞缩放

超胞的 k 点密度应与原胞收敛时的密度**相当**。VASP Wiki 给出的例子：Si 原胞需要约 12×12×12 才能得到准确的 Γ 点频率，2×2×2 超胞用 6×6×6、4×4×4 超胞用 3×3×3 即可。本例 8 原子惯用胞（约相当于原胞 2 倍边长）扩成 2×2×2，用 Γ 中心 3×3×3：

```text
Gamma-centred 3x3x3 mesh for the 2x2x2 supercell (10.9 A): same density as ~6x6x6 in the unit cell
0
Gamma
3 3 3
0 0 0
```

### 6.3 批量建目录与提交

示例的 [run_disp.sh](../examples/12-si-phonopy-fd/run_disp.sh) 自动识别 phonopy 版本（有 `phonopy-init` 就用新写法）：

```bash
bash run_disp.sh prepare     # phonopy-init -d ... 并建立 disp-001/{POSCAR,INCAR,KPOINTS,POTCAR}
bash run_disp.sh submit      # 每个 disp-XXX 中 sbatch ../job.slurm
bash run_disp.sh check       # 检查所有 vasprun.xml 是否完整（以 </modeling> 结尾）
bash run_disp.sh collect     # phonopy-init -f disp-*/vasprun.xml -> FORCE_SETS
```

位移超胞很多（低对称、多原子）时，用 Slurm job array 提交更方便，见 Common 第2章。

### 6.4 检查单点计算

```bash
for d in disp-*/; do
  n=$(grep -c "^RMM\|^DAV" $d/OSZICAR)
  echo "$d  SCF steps: $n  $(grep -c 'F=' $d/OSZICAR) ionic step(s)"
done
grep "TOTAL-FORCE" -A 3 disp-001/OUTCAR | head -5   # 第一个原子（被位移的原子）受力明显
```

- 每个目录只应有 1 个离子步（OSZICAR 中一行 `F=`）；
- SCF 步数不应等于 NELM（否则没收敛）；
- `run_disp.sh check` 确认 vasprun.xml 完整。

---

## 7. 第四步：收集力与后处理

### 7.1 生成 FORCE_SETS

```bash
phonopy-init -f disp-*/vasprun.xml          # phonopy >= 4.0
# phonopy -f disp-*/vasprun.xml             # phonopy 2.x / 3.x
```

`disp-*/` 的通配展开顺序必须与 POSCAR-001、002 … 的顺序一致（三位数编号按字典序即可）。后处理时 phonopy 会打印 `Max drift of force constants`——力常数的平移不变性偏差，数值大（> 0.01 eV/Å² 量级）说明力不够准。

### 7.2 band.conf：声子谱

```text
# phonopy settings: phonon band structure of diamond Si
# The input POSCAR is the 8-atom conventional cell; phonons are plotted in the
# 2-atom fcc primitive cell defined by PRIMITIVE_AXES.
ATOM_NAME = Si
DIM = 2 2 2
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
# Gamma-X-W-K-Gamma-L in reduced coordinates of the primitive reciprocal lattice
BAND = 0 0 0  1/2 0 1/2  1/2 1/4 3/4  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2
BAND_LABELS = $\Gamma$ X W K $\Gamma$ L
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
```

| 标签 | 说明 |
|---|---|
| `DIM` | 超胞矩阵，必须与生成位移时一致（3 个数为对角，9 个数为一般矩阵） |
| `PRIMITIVE_AXES` | 从输入胞到原胞的变换；fcc 惯用胞 → 原胞用上面的矩阵；`AUTO` 自动判断（v4 默认） |
| `BAND` | 高对称点（**原胞倒格子**的分数坐标），依次连线；逗号表示断开 |
| `BAND_LABELS` | 标签数 = 路径段数 + 1；LaTeX 写在 `$...$` 中 |
| `BAND_POINTS` | 每段的采样点数（默认 51） |
| `BAND_CONNECTION` | 按本征矢判断能带连接，正确处理交叉 |
| `BAND = AUTO` | 用 seekpath 自动生成路径（需安装 seekpath；坐标基于 seekpath 的标准原胞约定） |

运行：

```bash
phonopy --config band.conf -p -s           # phonopy >= 4.x
# phonopy -p -s band.conf                  # phonopy 2.x / 3.x
phonopy-bandplot --gnuplot band.yaml > band.dat   # 导出纯文本，便于 Origin/matplotlib
```

输出：`band.yaml`（全部数据）、`band.pdf`（`-s` 保存）。phonopy 还会把本次计算的完整设置写入 `phonopy.yaml`。

**结果检查**：Si 应无虚频，Γ 点三支声学模为 0，光学模在 Γ 点约 15 THz 量级（PBE 近似值；实验 Raman 峰约 520 cm⁻¹ ≈ 15.6 THz）。

### 7.3 mesh.conf：声子 DOS 与热力学

```text
# phonopy settings: phonon DOS and thermal properties of diamond Si
ATOM_NAME = Si
DIM = 2 2 2
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
MESH = 24 24 24
GAMMA_CENTER = .TRUE.
TPROP = .TRUE.
TMIN = 0
TMAX = 1000
TSTEP = 10
```

```bash
phonopy --config mesh.conf -p -s           # 新版；旧版: phonopy -p -s mesh.conf
phonopy --config pdos.conf -p -s           # 原子投影 DOS（PDOS = 1, 2）
ls total_dos.dat projected_dos.dat thermal_properties.yaml
```

| 标签 | 说明 |
|---|---|
| `MESH`（同 `MP`） | q 点网格；DOS 与热力学要对网格收敛（24³ → 32³ 比较） |
| `GAMMA_CENTER` | Γ 中心网格；fcc 等晶格下偶数 MP 网格会破坏对称性、增加不可约点数 |
| `TPROP` | 计算 F、S、C_V（`thermal_properties.yaml`） |
| `TMIN/TMAX/TSTEP` | 温度范围（默认 0、1000、10 K） |
| `PDOS` | 投影 DOS，`PDOS = 1 2, 3 4` 表示原子 1+2 一条曲线、3+4 一条 |
| `EXCLUDE_GAMMA_ACOUSTIC` | （phonopy ≥ 4.7.1）在热力学求和中排除 Γ 点声学模，避免近零频率引入的噪声 |

`thermal_properties.yaml` 中的单位：自由能 kJ/mol、熵与热容 J/K/mol，均**以原胞（PRIMITIVE_AXES 定义的胞）为单位**。换算成 eV/原胞：1 kJ/mol ≈ 0.010364 eV。

### 7.4 用 Python 读取结果

```python
import phonopy
ph = phonopy.load("phonopy_disp.yaml", force_sets_filename="FORCE_SETS",
                  primitive_matrix=[[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
ph.run_qpoints([[0, 0, 0]])
print("Gamma frequencies (THz):", ph.qpoints.frequencies[0])
ph.run_mesh([24, 24, 24], is_gamma_center=True)
ph.run_thermal_properties(t_min=0, t_max=1000, t_step=100)
tp = ph.thermal_properties
for T, F, S, Cv in zip(tp.temperatures, tp.free_energy, tp.entropy, tp.heat_capacity):
    print(f"{T:6.0f} K  F={F:8.3f} kJ/mol  S={S:7.3f} J/K/mol  Cv={Cv:7.3f} J/K/mol")
```

（以上属性写法在 phonopy 4.10 中验证可用；旧代码里常见的 `get_qpoints_dict()`、`get_thermal_properties_dict()` 在 4.10 中仍可运行，但会给出 DeprecationWarning。）

更多出版级作图方法（sumo 的声子绘图等）见 Common 第5章。

---

## 8. DFPT 路线：IBRION = 7/8 + phonopy

### 8.1 流程

```bash
cd VASP/examples/12-mgo-dfpt-nac
# 0) POSCAR-unitcell 先用 INCAR_relax 紧致优化（同 §4），并准备 POTCAR (Mg_pv 或 Mg, O)
bash setup_mgo.sh          # phonopy-init -d 生成 SPOSCAR，建立 dfpt/ 与 born/
(cd dfpt && sbatch ../../12-si-phonopy-fd/job.slurm)
(cd born && sbatch ../../12-si-phonopy-fd/job.slurm)
bash setup_mgo.sh post     # FORCE_CONSTANTS + BORN + 声子谱
```

DFPT 只需要**一次**超胞计算（SPOSCAR），VASP 内部完成所有微扰。`setup_mgo.sh` 中的关键命令：

```bash
phonopy-init -d --dim 2 2 2 -c POSCAR-unitcell   # 只用到 SPOSCAR 和 phonopy_disp.yaml
# ... dfpt/ 中 IBRION = 8 计算完成后
grep -c 'name="hessian"' dfpt/vasprun.xml       # vasprun.xml 中必须有 hessian
phonopy-init --fc dfpt/vasprun.xml               # -> FORCE_CONSTANTS（旧版: phonopy --fc）
```

### 8.2 INCAR_dfpt

```bash
SYSTEM  = MgO 2x2x2 supercell (64 atoms) - DFPT force constants (IBRION = 8)

PREC    = Accurate
ENCUT   = 520          # same as the relaxation
EDIFF   = 1E-8
ALGO    = Normal
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
# ADDGRID = .TRUE.     # used in the phonopy DFPT example; test before use

IBRION  = 8            # DFPT, symmetry-reduced displacements (7 = no symmetry)
# Do NOT set NCORE/NPAR (DFPT linear response needs NCORE = 1)
# Do NOT set IVDW (DFPT does not support vdW corrections; use finite differences)

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

band.conf 中加 `FORCE_CONSTANTS = READ`（phonopy ≥ 4 在 FORCE_CONSTANTS 文件存在时默认就会读取，写上兼容旧版）。

### 8.3 有限位移与 DFPT 怎么选

| 对比项 | 有限位移（phonopy + IBRION=−1） | DFPT（IBRION=8） |
|---|---|---|
| 泛函 | 任意（LDA/GGA/meta-GGA/杂化/DFT+U/vdW） | 仅 LDA/GGA（官方文档） |
| vdW 修正（IVDW） | 支持 | **不支持** |
| 并行 | 各位移独立，天然高通量并行 | 单个大作业；`NCORE` 必须为 1 |
| 精度控制 | 需选位移幅度（默认 0.01 Å 一般可靠） | 无位移幅度问题 |
| 中断后 | 只需重算失败的位移 | 整个重算（有限差分可用 CHECKPOINT_FD 续算，§11） |
| 适用 | **默认推荐** | 快速对照、检验位移幅度 |

---

## 9. 极性材料：Born 有效电荷与非解析项修正

### 9.1 LEPSILON 计算

在**单胞**（不是超胞）上用 DFPT 对电场的响应求 $\varepsilon^\infty$ 和 $Z^*$，k 点要比超胞计算密（约两倍）：

```bash
SYSTEM  = MgO conventional cell - Born effective charges and eps_inf (LEPSILON)

PREC    = Accurate
ENCUT   = 520
EDIFF   = 1E-8
ALGO    = Normal
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.

IBRION  = -1
NSW     = 0
LEPSILON = .TRUE.      # DFPT response to an electric field: eps_inf and Z*
LRPA    = .FALSE.      # include exchange-correlation local fields (needed for Z*)
# Do NOT set NCORE/NPAR here either.

LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

OUTCAR 中查看：

```bash
grep -A 4 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)" born/OUTCAR
grep -A 12 "BORN EFFECTIVE CHARGES" born/OUTCAR | head -14
```

VASP Wiki 的 MgO 示例给出 $\varepsilon^\infty\approx3.13$、$Z^*_\mathrm{Mg}\approx+1.97$、$Z^*_\mathrm{O}\approx-1.97$（来自官方页面，作为量级参考）。检查要点：所有原子 $Z^*$ 之和应接近 0（电荷中性求和规则），偏差大说明 k 点不够。

> [!NOTE]
> 杂化泛函不能用 LEPSILON（DFPT 不支持 HF/杂化），可改用 `LCALCEPS = .TRUE.`（有限电场法），详见[第14章](14-optics-dielectric.md)。VASP 6.6.1 修复了 OpenACC GPU 版 LCALCEPS 计算 Born 有效电荷的错误，GPU 用户应升级。

### 9.2 生成 BORN 文件

```bash
cd born
phonopy-vasp-born --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0" > ../BORN
cat ../BORN
```

BORN 文件的第一行是注释（如 `# epsilon and Z* of atoms 1 5`，列出对称不等价原子的序号），第二行是 $\varepsilon^\infty$ 的 9 个分量，之后每个不等价原子一行 $Z^*$。

- `--pa` 必须与 band.conf 中的 `PRIMITIVE_AXES` 一致；在 phonopy 4.10 中实测 `phonopy-vasp-born` 的 `--pa` 只接受数值矩阵（不接受 `auto` 或 `F` 这类写法）；
- 默认会按晶体对称性对张量做对称化，`--nost` 关闭；`--outcar` 改为读取 OUTCAR + POSCAR。

### 9.3 打开 NAC 后处理

```text
# phonopy settings: MgO phonon bands from VASP-DFPT force constants with NAC
ATOM_NAME = Mg O
DIM = 2 2 2
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
FORCE_CONSTANTS = READ
NAC = .TRUE.
BAND = 0 0 0  1/2 0 1/2  1/2 1/4 3/4  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2
BAND_LABELS = $\Gamma$ X W K $\Gamma$ L
BAND_POINTS = 101
```

```bash
phonopy --config band.conf -p -s           # phonopy >= 4: BORN 存在即自动 NAC（日志显示 NAC parameters were read from "BORN"）
# phonopy -c POSCAR-unitcell -p -s band.conf   # 旧版，NAC = .TRUE. 生效
```

- 打开 NAC 后，Γ 点附近光学支出现 **LO–TO 劈裂**，且 LO 频率依赖 q→0 的方向（声子谱在 Γ 两侧不连续是正常的）；
- 需要某一方向的 Γ 点 LO 频率时，用 `QPOINTS = 0 0 0` 加 `Q_DIRECTION = 1 0 0`；
- MgO 的 TO/LO 实验值约为 12 THz / 21 THz 量级（近似），可用来粗略检验。

---

## 10. 虚频诊断

| 现象 | 可能原因 | 处理 |
|---|---|---|
| Γ 点附近声学支出现很小的虚频（< 0.3 THz 量级） | 力精度不足、ASR 未满足、超胞偏小 | `EDIFF = 1E-8`、`PREC = Accurate`、`LREAL = .FALSE.`；增大超胞；新版 phonopy 默认用 symfc 投影对力常数做对称化 |
| 整条声学/光学支在布里渊区边界变虚 | 结构未充分优化（残余应力/力）、或真实的动力学不稳定 | 先收紧优化（ISIF=3、EDIFFG=−1E-3）重新算；若仍存在，可能是真实软模 |
| 只在某个特定 q 点（如 X、M、R）出现明显虚频 | 真实的结构不稳定（如钙钛矿八面体倾转、CDW） | 沿软模本征矢位移、在对应超胞中重新优化，得到更稳定的低对称结构 |
| 金属体系出现奇异虚频 | SIGMA/k 点不够 | 测试 SIGMA 与 k 点；金属需要更密的 k 点 |
| 二维材料 ZA 支在 Γ 附近小虚频 | 数值误差 + 旋转不变性未满足 | 见 §12 |
| 磁性材料虚频 | 磁构型错误或超胞破坏了磁有序 | 用正确的磁胞；在 phonopy 中设置 `MAGMOM` |
| 只在某些位移计算中力异常大 | SCF 未收敛、某个 disp 作业中断 | `run_disp.sh check` 重查；重算该位移 |

> [!TIP]
> VASP Wiki 的 *How to handle imaginary phonon modes* 页面给出了处理真实软模的系统流程：先在超胞中定位软模所在的 q 点，构造能把该 q 折叠到 Γ 的最小超胞，再用 IRC（`IBRION = 40`）或沿本征矢的冻结声子扫描找到真正的能量极小，最后重新优化并验证。官方同时提醒：未充分弛豫的结构是"假虚频"最常见的来源。

phonopy 端可以用 `MODULATION` 标签沿某个模构造调制结构，帮助沿软模方向搜索新相（用法见 phonopy 文档）。

---

## 11. VASP 原生声子功能：IBRION = 6、CHECKPOINT_FD 与 LPHON_DISPERSION

### 11.1 在超胞中计算力常数（IBRION = 6）

不借助 phonopy，VASP 也可以在超胞中直接做有限差分（POSCAR 用超胞，例如 phonopy 生成的 SPOSCAR）：

```bash
SYSTEM  = Si 2x2x2 supercell - finite differences inside VASP
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-8
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
IBRION  = 6            # symmetry-reduced finite differences (5 = all displacements)
POTIM   = 0.015        # displacement step (A); too-large values are reset to 0.015
NFREE   = 2            # central differences
ISIF    = 2            # ISIF >= 3 would also compute elastic constants (chapter 13)
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

```bash
grep "THz" OUTCAR | head        # 超胞 Gamma 点的 3N 个模；f/i 为虚频；最后 3 个是平移模
```

### 11.2 VASP 6.6：CHECKPOINT_FD 断点续算与拆分

VASP 6.6.0 新增 `CHECKPOINT_FD`（**仅适用于 IBRION = 6，且需要编译时打开 HDF5**），把每个位移的结果写入 `vaspcheckfd.h5`：

| 取值 | 作用 |
|---|---|
| `RESET`（默认） | 新建/覆盖 vaspcheckfd.h5，每完成一个位移更新一次 |
| `CONTINUE` | 读取已有的 vaspcheckfd.h5，从最后完成的位移继续；也用于汇总拆分后的计算 |
| `NONE` | 不写 vaspcheckfd.h5（6.6.0 之前的行为） |
| `PREPARE` | 只做未位移结构的 SCF，写出 `CONTCAR_disp-N`（每个位移一个结构）后停止 |
| `SINGLE` | 运行 PREPARE 生成的单个位移计算 |

**作业中断后续算**：在 INCAR 中加 `CHECKPOINT_FD = CONTINUE`，在原目录重新提交，stdout 中出现 `Continuing from previous run` 即表示从断点继续。

**拆分成独立作业**（大体系可同时排队多个小作业）：

```bash
# 1) 准备：INCAR 中 CHECKPOINT_FD = PREPARE，运行一次 -> CONTCAR_disp-1 ... CONTCAR_disp-N
# 2) 每个位移一个目录，CHECKPOINT_FD 改为 SINGLE
for f in CONTCAR_disp-*; do
    n=${f#CONTCAR_disp-}
    mkdir -p disp-$n
    cp $f disp-$n/POSCAR
    cp KPOINTS POTCAR disp-$n/
    sed 's/^ *CHECKPOINT_FD.*/CHECKPOINT_FD = SINGLE/' INCAR > disp-$n/INCAR
done
# 各 disp-N 分别提交运行（官方示例中还复用父目录的 CHGCAR 以加速 SCF）
# 3) 全部完成后，在父目录把 CHECKPOINT_FD 改为 CONTINUE 再运行一次，汇总得到声子频率
```

> [!NOTE]
> 拆分模式要求所有子计算都已完成，否则汇总步骤会报错。具体的目录命名与文件要求以官方页面 [Restarting finite differences calculations](https://www.vasp.at/wiki/index.php/Restarting_finite_differences_calculations) 为准。

### 11.3 LPHON_DISPERSION：VASP 内部插值声子色散

VASP ≥ 6.3.2 可以把超胞力常数傅里叶插值到 VASP 自动确定的原胞上：

```bash
# 在力常数计算（IBRION = 5-8）的同一 INCAR 中加：
LPHON_DISPERSION = .TRUE.     # 读 QPOINTS（格式同 KPOINTS，常用 line mode）
PHON_NWRITE      = 1          # OUTCAR 输出详细程度（6.3.2+）
```

更省时的做法是**复用已经算好的力常数**（VASP ≥ 6.4.0）：

```bash
cp vaspout.h5 vaspin.h5
# INCAR:
LPHON_READ_FORCE_CONSTANTS = .TRUE.
LPHON_DISPERSION           = .TRUE.
# 修改 QPOINTS 即可换路径，VASP 只做色散计算后退出
```

声子 DOS：在 QPOINTS 中给均匀网格，设置 `PHON_DOS = 1`（高斯展宽，宽度 `PHON_SIGMA`）或 `2`（四面体），`PHON_NEDOS` 控制能量点数。极性材料用 `LPHON_POLAR = .TRUE.`，并在超胞 INCAR 中手工写入 `PHON_DIELECTRIC`（9 个数）和 `PHON_BORN_CHARGES`（每原子 9 个数），数值取自单胞 LEPSILON/LCALCEPS 计算。

> [!WARNING]
> QPOINTS 中 q 点的分数坐标是相对于 **VASP 自动识别的原胞**的倒格子而言的。对复杂结构，请参考官方 [Computing the phonon dispersion and DOS](https://www.vasp.at/wiki/index.php/Computing_the_phonon_dispersion_and_DOS) 页面及其推荐的 k 路径工具确认坐标约定；不确定时，phonopy 路线更透明。

### 11.4 用 py4vasp 画图

VASP 编译时打开 HDF5 后，声子数据写入 vaspout.h5（VASP 6.6 起声子频率也写入）。py4vasp 提供：

```python
import py4vasp
calc = py4vasp.Calculation.from_path(".")
calc.phonon.band.plot()          # 需要 LPHON_DISPERSION 计算
calc.phonon.dos.plot()           # 需要 PHON_DOS 计算
modes = calc.phonon.mode.read()  # Gamma 点频率与本征矢（IBRION = 5-8）
```

---

## 12. 二维材料：石墨烯声子

二维材料声子计算的特殊之处：

1. **真空层**：POSCAR 中 c ≥ 15–20 Å；KPOINTS 中 z 方向取 1，六角晶格用 Γ 中心网格；
2. **超胞只在面内扩展**：`--dim 5 5 1`（不要在真空方向扩胞）；
3. **晶格优化**：只优化面内晶格常数（VASP ≥ 6.4 用 `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` 配合 `ISIF = 3`，旧版本用能量–晶格常数扫描）；
4. **ZA 支**（面外弯曲模）在 Γ 附近呈二次色散 $\omega\propto q^2$，数值误差很容易让它在 Γ 附近出现很小的虚频"口袋"。处理：增大超胞、`EDIFF = 1E-8`、`PREC = Accurate`、对力常数施加平移/旋转不变性（新版 phonopy 默认的 symfc 投影对称化有帮助；更严格的旋转求和规则可借助 hiphive 等工具）。只出现在 Γ 点极近处、幅度很小的 ZA 虚频通常不代表真实不稳定。

示例 [12-graphene-phonon](../examples/12-graphene-phonon/README.md) 的 band.conf：

```text
# phonopy settings: phonon bands of graphene (POSCAR = primitive 2-atom cell)
ATOM_NAME = C
DIM = 5 5 1
PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1
# Gamma-M-K-Gamma for a hexagonal cell with a 120 degree angle between a1 and a2
BAND = 0 0 0  1/2 0 0  1/3 1/3 0  0 0 0
BAND_LABELS = $\Gamma$ M K $\Gamma$
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
```

```bash
cd VASP/examples/12-graphene-phonon
# 优化后：
DIM="5 5 1" bash ../12-si-phonopy-fd/run_disp.sh prepare     # 50 原子超胞，只需 1 个位移
(cd disp-001 && sbatch ../../12-si-phonopy-fd/job.slurm)
DIM="5 5 1" bash ../12-si-phonopy-fd/run_disp.sh collect
phonopy --config band.conf -p -s
```

石墨烯声子谱的特征：K 点附近最高光学支有 Kohn 异常（对 k 点和 SIGMA 敏感），Γ 点 E2g 光学模约 1580 cm⁻¹ 量级（对应 Raman G 峰，近似比较）。二维材料的介电屏蔽与 NAC 处理比体相复杂（2D 的 LO–TO 行为不同），这里不展开。更多二维材料设置见[第8章](08-vdw-2d.md)。

---

## 13. 进阶：准谐近似（QHA）与 phono3py 晶格热导率

### 13.1 QHA 热膨胀

步骤：

1. 在平衡体积附近取 ≥ 5 个（推荐 7–11 个）体积，例如 −5% 到 +5%，每个体积固定晶胞、只弛豫原子（`ISIF = 2`；非立方体系要在每个体积下优化形状，可用 `ISIF = 4`）；
2. 每个体积计算静态能量 $U(V)$（写入 e-v.dat：体积 Å³、能量 eV，**按 phonopy 的原胞/单胞计**，与 thermal_properties.yaml 一致）；
3. 每个体积做一次完整的声子计算，得到 `thermal_properties.yaml`（所有体积用相同的温度范围和步长，TMAX 要比目标温度高）；
4. 运行：

```bash
# 把各体积的 thermal_properties.yaml 改名为 thermal_properties.yaml-{-5..5}
phonopy-qha -p -s e-v.dat thermal_properties.yaml-{-{5..1},{0..5}}
phonopy-qha -b e-v.dat                                  # 只拟合 E-V 得到 B0
```

输出包括 $V(T)$、体积热膨胀系数、$G(T)$、$B_T(T)$、$C_p(T)$ 等（`-s` 写文件）。金属可用 `phonopy-vasp-efe` 加入电子自由能（`--efe`）。phonopy-qha 计划在未来移除，新项目建议直接用 Python API `phonopy.qha` / `run_qha`。

### 13.2 phono3py 晶格热导率（简介）

晶格热导率需要三阶力常数（声子–声子散射）。phono3py 4.x 同样拆分出 `phono3py-init`：

```bash
phono3py-init -d --dim 2 2 2 --pa F -c POSCAR-unitcell     # phono3py_disp.yaml + POSCAR-xxxxx（数量可达数百）
# 每个 POSCAR-xxxxx 做单点（与本章 INCAR 相同的高精度设置）
phono3py-init --cf3 disp-{00001..00755}/vasprun.xml       # -> FORCES_FC3（编号范围按实际修改）
phono3py --mesh 11 11 11 --br                              # RTA 热导率
```

位移数随原子数和超胞快速增长，通常需要较小的三阶力常数超胞（或借助机器学习势/随机位移拟合）。可以用 `--dim-fc2` 为二阶力常数单独指定更大的超胞。详细流程见 [phono3py 文档](https://phonopy.github.io/phono3py/)。

---

## 14. 参数选择与收敛建议

| 参数 | 推荐 | 收敛判据 / 说明 |
|---|---|---|
| 结构优化 | ISIF=3，EDIFFG=−1E-3，EDIFF=1E-8，ENCUT ≥ 1.3×ENMAX | 外压 < 0.5 kB，力 < 1E-3 eV/Å |
| ENCUT | 与优化相同 | Γ 点光学模频率变化 < 0.05 THz |
| EDIFF | 1E-8 | 更松会引入虚频 |
| PREC / LREAL | Accurate / .FALSE. | — |
| 超胞 | 每方向 ≥ 10–12 Å，再增大一次比较 | 声子谱变化 < 0.1 THz（关注边界点） |
| k 点（超胞） | 与原胞收敛密度相当（原胞 N×N×N → n 倍超胞约 N/n） | Γ 点频率稳定 |
| 位移幅度 | 0.01 Å（phonopy 默认） | 0.005–0.02 Å 结果应一致；可与 DFPT 对比 |
| MESH（DOS、热力学） | 20³–40³，Γ 中心 | 热容/自由能变化 < 0.1% |
| BORN 计算的 k 点 | 约为超胞等效密度的 2 倍 | ε∞、Z* 变化 < 1%，ΣZ* ≈ 0 |
| 金属 SIGMA | 0.1–0.2（MP），需测试 | 声子谱对 SIGMA 不敏感 |

---

## 15. 常见问题与报错

| 现象 / 报错 | 原因 | 解决 |
|---|---|---|
| `"band.conf" could not be read as a phonopy.yaml-like file. If "band.conf" is a configuration file, give it by --config.` | phonopy 4.x 不再接受配置文件作位置参数 | `phonopy --config band.conf -p -s` |
| phonopy 4.x 运行 `phonopy -d ...` / `phonopy -f ...` 报迁移错误 | 前处理选项已移到 phonopy-init | 改用 `phonopy-init -d` / `phonopy-init -f` / `phonopy-init --fc` |
| `phonopy-vasp-born --pa auto` 报 `could not convert string to float` | 该工具的 `--pa` 只接受数值矩阵（4.10 实测） | 写出完整矩阵，如 `--pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"` |
| 升级 phonopy 后声子谱的 q 路径/形状变了 | v4 默认 `PRIMITIVE_AXES = AUTO` | 在 conf 中显式写 `PRIMITIVE_AXES`（或 `--pa P` 保持 v3 行为） |
| DFPT/LEPSILON 报 `PEAD routines do not work for NCORE` 或提示 `remove the tag NPAR from the INCAR file` | 线性响应不支持 NCORE/NPAR > 1 | 删除 INCAR 中的 NCORE/NPAR |
| DFPT 计算加了 IVDW | DFPT 不支持 vdW 修正 | 改用有限位移（phonopy 或 IBRION=5/6） |
| `phonopy-init --fc` 报错 / vasprun.xml 中没有 hessian | 不是 IBRION=7/8 计算或作业未完成 | 检查 INCAR 与 vasprun.xml 是否完整 |
| FORCE_SETS 生成失败，提示文件数与位移数不符 | disp 目录数量或顺序与 phonopy_disp.yaml 不一致 | 用 `disp-*/vasprun.xml` 按编号顺序给出；检查是否漏算 |
| 声子谱在 Γ 附近有小虚频 | 精度/超胞/ASR | 见 §10 |
| 整体出现大虚频 | 位移超胞被弛豫了（NSW>0）或结构未优化 | `IBRION = -1`、`NSW = 0`；重新紧致优化 |
| `thermal_properties` 出现 GammaAcousticWarning | Γ 点声学模近零频率进入求和 | phonopy ≥ 4.7.1 可设 `EXCLUDE_GAMMA_ACOUSTIC = .TRUE.` |
| CHECKPOINT_FD 无效或报错 | VASP < 6.6.0、未编译 HDF5、或 IBRION ≠ 6 | 升级并打开 `-DVASP_HDF5`；只用于 IBRION = 6 |
| LPHON_READ_FORCE_CONSTANTS 找不到力常数 | 没把 vaspout.h5 复制为 vaspin.h5，或 VASP < 6.4.0 | `cp vaspout.h5 vaspin.h5` |

更多通用报错见[第18章](18-postprocess-troubleshooting.md)。

---

## 练习

1. 用 2×2×2 与 3×3×3 超胞分别计算 Si 声子谱，比较 X、L 点频率，判断超胞收敛情况。
2. 把位移幅度改为 0.005 Å 和 0.03 Å（`--amplitude`），观察声子谱变化。
3. 对 MgO 分别画"无 NAC"（`--nonac`）和"有 NAC"的声子谱，读出 Γ 点 TO 与 LO 频率。
4. 用 IBRION=6 在 Si 超胞中直接计算 Γ 点频率，与 phonopy 结果比较；若使用 VASP ≥ 6.6，人为中断作业后用 `CHECKPOINT_FD = CONTINUE` 续算。
5. 计算 Si 在 −4%…+4% 共 9 个体积下的声子，用 phonopy-qha 得到 300 K 的热膨胀系数（提示：Si 在低温有负热膨胀，PBE 能定性复现）。
6. 计算石墨烯 4×4×1、5×5×1、6×6×1 超胞的声子谱，观察 ZA 支在 Γ 附近的行为。

---

## 参考资料

**官方文档**
- phonopy：[VASP & phonopy calculation](https://phonopy.github.io/phonopy/vasp.html)、[VASP-DFPT & phonopy calculation](https://phonopy.github.io/phonopy/vasp-dfpt.html)、[phonopy-init command](https://phonopy.github.io/phonopy/phonopy-init.html)、[Setting tags](https://phonopy.github.io/phonopy/setting-tags.html)、[Migrating from phonopy v3 to v4](https://phonopy.github.io/phonopy/migration-v4.html)、[Quasi harmonic approximation](https://phonopy.github.io/phonopy/qha.html)、[Auxiliary tools（phonopy-vasp-born 等）](https://phonopy.github.io/phonopy/auxiliary-tools.html)
- [phono3py 文档](https://phonopy.github.io/phono3py/)
- VASP Wiki：[Phonons from finite differences](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences)、[Phonons from density-functional-perturbation theory](https://www.vasp.at/wiki/index.php/Phonons_from_density-functional-perturbation_theory)、[Computing the phonon dispersion and DOS](https://www.vasp.at/wiki/index.php/Computing_the_phonon_dispersion_and_DOS)、[CHECKPOINT_FD](https://www.vasp.at/wiki/index.php/CHECKPOINT_FD)、[LPHON_DISPERSION](https://www.vasp.at/wiki/index.php/LPHON_DISPERSION)、[LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON)、[How to handle imaginary phonon modes](https://www.vasp.at/wiki/index.php/How_to_handle_imaginary_phonon_modes)
- [VASP 官方 Tutorials（声子、电声部分）](https://www.vasp.at/tutorials/latest/)、[py4vasp 文档](https://vasp.at/py4vasp/latest/)

**相关 GitHub 仓库**（只链接，不转载代码）
- [phonopy/phonopy](https://github.com/phonopy/phonopy) —— 声子计算标准工具（4.10.0），有限位移、VASP-DFPT 接口、热力学、QHA
- [phonopy/phono3py](https://github.com/phonopy/phono3py) —— 三阶力常数与晶格热导率
- [gabkrenzer/PhononFlow](https://github.com/gabkrenzer/PhononFlow) —— VASP + phonopy 流程：声子谱、热力学、QHA、Grüneisen
- [Stanford-MCTG/VASP-plot-modes](https://github.com/Stanford-MCTG/VASP-plot-modes) —— 把 DFPT 振动模式导出为 VESTA 图像
- [raman-sc/VASP](https://github.com/raman-sc/VASP) —— 有限位移极化率法计算 Raman 活性（声子扩展）
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) —— 读取 vaspout.h5，绘制 VASP 原生声子结果
- [materialsproject/atomate2](https://github.com/materialsproject/atomate2) —— 含声子工作流的高通量框架
- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) —— 出版级绘图（含声子谱）
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**
- A. Togo, I. Tanaka, "First principles phonon calculations in materials science", *Scr. Mater.* 108, 1 (2015). [doi:10.1016/j.scriptamat.2015.07.021](https://doi.org/10.1016/j.scriptamat.2015.07.021)
- A. Togo, L. Chaput, T. Tadano, I. Tanaka, "Implementation strategies in phonopy and phono3py", *J. Phys.: Condens. Matter* 35, 353001 (2023). [doi:10.1088/1361-648X/acd831](https://doi.org/10.1088/1361-648X/acd831)
- K. Parlinski, Z. Q. Li, Y. Kawazoe, "First-principles determination of the soft mode in cubic ZrO2", *Phys. Rev. Lett.* 78, 4063 (1997). [doi:10.1103/PhysRevLett.78.4063](https://doi.org/10.1103/PhysRevLett.78.4063)
- S. Baroni, S. de Gironcoli, A. Dal Corso, P. Giannozzi, "Phonons and related crystal properties from density-functional perturbation theory", *Rev. Mod. Phys.* 73, 515 (2001). [doi:10.1103/RevModPhys.73.515](https://doi.org/10.1103/RevModPhys.73.515)
- M. Gajdoš et al., "Linear optical properties in the projector-augmented wave methodology", *Phys. Rev. B* 73, 045112 (2006). [doi:10.1103/PhysRevB.73.045112](https://doi.org/10.1103/PhysRevB.73.045112)
- X. Gonze, C. Lee, "Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory", *Phys. Rev. B* 55, 10355 (1997). [doi:10.1103/PhysRevB.55.10355](https://doi.org/10.1103/PhysRevB.55.10355)
- Y. Wang et al., "A mixed-space approach to first-principles calculations of phonon frequencies for polar materials", *J. Phys.: Condens. Matter* 22, 202201 (2010). [doi:10.1088/0953-8984/22/20/202201](https://doi.org/10.1088/0953-8984/22/20/202201)

---

[← 上一章](11-neb-dimer.md) | [返回目录](../README.md) | [下一章 →](13-elastic.md)
