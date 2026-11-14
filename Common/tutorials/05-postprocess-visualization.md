# 第5章 后处理与可视化：VASPKIT、sumo、pyprocar、VESTA、Multiwfn 与出版级作图（Post-processing & Visualization: VASPKIT, sumo, PyProcar, py4vasp, VESTA, OVITO, Multiwfn and Publication-quality Figures）

本章讲的是算完之后的工作：怎样把 VASP、Gaussian、CP2K、OpenMX 的输出变成可信、好看、能直接放进论文的图。前半部分以 GaAs 为例，用 VASPKIT 1.5.1（交互菜单与非交互批处理、能带 211/213、DOS 111/113、光学 711、功函数 426、热校正 501/502、差分电荷 314）、sumo 3.0（sumo-kgen、sumo-bandplot、sumo-dosplot、sumo-bandstats）、PyProcar 6（新的 `dirname` API、fat band、费米面）、py4vasp（vaspout.h5）和 matplotlib 样式文件处理能带与态密度，并说明自旋极化与 SOC 能带怎么画。后半部分介绍三维可视化：VESTA 的等值面、截面和差分电荷，OVITO 渲染 AIMD 轨迹，VMD + Multiwfn 画轨道与弱相互作用等值面，Avogadro/Jmol 与 phonon website 播放振动动画，以及在 Jupyter 中用 nglview、weas-widget，在浏览器中用 matterviz 交互查看结构（band structure plotting, PDOS, fat bands, effective mass, isosurface, charge density difference, trajectory rendering, matplotlib mplstyle）。

> [!NOTE]
> 本章是跨软件的工具总览与作图规范。各软件“怎么算”分别见：VASP [第5章 能带与态密度](../../VASP/tutorials/05-band-dos.md)、[第6章 杂化泛函](../../VASP/tutorials/06-hybrid-metagga.md)、[第14章 光学](../../VASP/tutorials/14-optics-dielectric.md)、[第16章 电荷与成键](../../VASP/tutorials/16-charge-bonding.md)、[第17章 反折叠与费米面](../../VASP/tutorials/17-advanced-electronic.md)、[第18章 后处理生态](../../VASP/tutorials/18-postprocess-troubleshooting.md)；Gaussian [第16章 Multiwfn](../../Gaussian/tutorials/16-postprocess-multiwfn.md)；CP2K [第10章 DOS/能带](../../CP2K/tutorials/10-dos-pdos-band.md)、[第11章 cube 与功函数](../../CP2K/tutorials/11-cube-charge-workfunction.md)；OpenMX [第8章 能带](../../OpenMX/tutorials/08-band.md)、[第9章 DOS](../../OpenMX/tutorials/09-dos-pdos-cohp-optics.md)、[第10章 cube](../../OpenMX/tutorials/10-charge-mo-cube.md)。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 正确准备后处理所需的计算 | 知道能带要 `ICHARG = 11` + line-mode、DOS 要 `ISMEAR = -5` + `NEDOS` + `LORBIT = 11`，杂化泛函能带要零权重 k 点 |
| 用 VASPKIT 批量处理 | 配置 `~/.vaspkit`，用 `echo -e` / 应答文件把交互菜单变成脚本，对几十个目录一次性提取能带、DOS、带隙、功函数、热校正 |
| 出版级能带/DOS | 用 sumo 一条命令出 PDF，用 PyProcar/py4vasp 画投影能带，用 matplotlib + mplstyle 完全控制版式 |
| 自旋与 SOC | 正确画自旋向上/向下能带和镜像 DOS，理解 SOC 计算中投影数据的变化 |
| 三维体数据 | 在 VESTA 中设置等值面、截面，做 AB − A − B 差分电荷，导出高分辨率图 |
| 轨迹与振动 | 用 OVITO 在集群上无界面渲染 XDATCAR，用 phonon website / Avogadro / Jmol 播放振动 |
| 波函数分析 | 用 Multiwfn 命令流批量处理 fchk/molden，用 VMD 渲染等值面 |
| 作图规范 | 掌握图幅尺寸、字体嵌入、矢量格式、色盲友好配色和能量零点标注 |

先根据“要画什么”选工具：

| 要画的东西 | 首选工具 | 需要的文件 | 备选 |
|---|---|---|---|
| 能带、带隙、有效质量 | sumo（PDF）、VASPKIT 211/911 | vasprun.xml 或 EIGENVAL + KPOINTS | py4vasp、PyProcar、pymatgen |
| fat band（投影能带） | sumo `--project`、PyProcar | vasprun.xml / PROCAR（`LORBIT = 11`） | VASPKIT 213/214、py4vasp |
| DOS / PDOS | sumo-dosplot、VASPKIT 111/113/115 | vasprun.xml / DOSCAR | py4vasp、PyProcar、pymatgen |
| HSE 能带 | sumo（自动拼接 split）、VASPKIT 252 | 零权重 KPOINTS 的 vasprun.xml | py4vasp（KPOINTS_OPT） |
| 光学谱 | VASPKIT 711、sumo-optplot | LOPTICS 的 vasprun.xml | py4vasp |
| 超胞能带反折叠 | easyunfold | WAVECAR | VASPKIT 28x |
| 费米面 | IFermi、PyProcar FermiHandler | 均匀密网格 vasprun.xml | VASPKIT 26x + FermiSurfer |
| 电荷密度、ELF、差分电荷 | VESTA | CHGCAR / ELFCAR / PARCHG / cube | VASPKIT 314、matterviz |
| AIMD 轨迹 | OVITO | XDATCAR / extxyz | VMD、nglview |
| 分子轨道、ESP、IGMH | Multiwfn + VMD | fchk / molden / wfx / cube | GaussView、Avogadro |
| 振动与声子动画 | Avogadro、Jmol（分子）；phonon website（声子） | Gaussian .log；phonopy band.yaml | GaussView |

配套示例：

- [05-band-dos-plotting](../examples/05-band-dos-plotting/README.md)：GaAs 的 SCF/能带/DOS/HSE 输入，sumo、PyProcar、py4vasp、matplotlib 四种作图方式；
- [05-vaspkit-batch](../examples/05-vaspkit-batch/README.md)：VASPKIT 批处理脚本、各功能的应答文件、光学/功函数/频率 INCAR；
- [05-visualization-3d](../examples/05-visualization-3d/README.md)：VESTA 操作清单、OVITO 渲染脚本、nglview/weas-widget 笔记本、Multiwfn 批处理、phonon website 用的 band.conf。

---

## 2. 原理速览

### 2.1 能带：沿路径的本征值

Kohn–Sham 方程在每个 k 点给出一组本征值 $\varepsilon_{n\mathbf{k}}$。能带图就是沿布里渊区高对称路径（如 fcc 的 L–Γ–X–U|K–Γ）把 $\varepsilon_{n\mathbf{k}}$ 连成曲线。路径上的点不参与构造电荷密度，所以通常先做自洽计算得到密度，再在固定密度下做非自洽计算（VASP 的 `ICHARG = 11`）。杂化泛函的哈密顿量依赖轨道而不只是密度，不能这样做，而要把路径点以**零权重**加到自洽网格中（或用 VASP ≥ 6.3 的 `KPOINTS_OPT`）。

路径约定不唯一：Setyawan–Curtarolo（pymatgen 默认）、Hinuma/SeeK-path、Bradley–Cracknell（sumo 默认）等对同一晶体可能给出不同的原胞取向和点标签。**同一张图的计算和标签必须来自同一约定**，生成路径的工具若改写了原胞（sumo 的 `POSCAR_prim`、VASPKIT 的 `PRIMCELL.vasp`），能带计算必须使用这个原胞。

### 2.2 态密度与投影

总态密度

$$g(E)=\sum_{n}\sum_{\mathbf{k}} w_{\mathbf{k}}\,\delta(E-\varepsilon_{n\mathbf{k}})$$

需要在整个布里渊区上积分，所以要用**均匀且较密**的网格，不能用能带路径。δ 函数要么用高斯/Methfessel–Paxton 展宽（`ISMEAR ≥ 0`，展宽宽度 `SIGMA`），要么用四面体积分（`ISMEAR = -5`，Blöchl 修正），后者对半导体和绝缘体的 DOS 形状更准确。

投影态密度（PDOS）把每个态投影到以原子 $i$ 为中心、角动量为 $lm$ 的球谐函数上：

$$g_{i,lm}(E)=\sum_{n\mathbf{k}} w_{\mathbf{k}}\,\bigl|\langle Y^{i}_{lm}\,|\,\psi_{n\mathbf{k}}\rangle\bigr|^{2}\,\delta(E-\varepsilon_{n\mathbf{k}})$$

VASP 中 `LORBIT = 11` 输出按 $lm$ 分解的投影（PROCAR、DOSCAR、vasprun.xml）。投影是在以原子为中心的球内做的；`LORBIT = 10/11` 时球半径取自 POTCAR（INCAR 中的 `RWIGS` 只对 `LORBIT = 0/1/2` 起作用），球之间的间隙区没有计入，所以各原子 PDOS 之和一般小于总 DOS——这是正常现象，不是作图错误。fat band 就是把同样的权重 $|\langle Y^{i}_{lm}|\psi_{n\mathbf{k}}\rangle|^2$ 画在能带上（颜色或点的大小）。

### 2.3 有效质量

在能带极值附近做抛物线近似：

$$\frac{1}{m^{*}}=\frac{1}{\hbar^{2}}\,\frac{\partial^{2}\varepsilon(\mathbf{k})}{\partial k^{2}}$$

sumo-bandstats、VASPKIT 912 等工具在 VBM/CBM 附近取几个 k 点拟合。结果对路径点密度和拟合点数敏感，非抛物线能带（如 GaAs 导带离 Γ 稍远处）尤其如此。

### 2.4 体数据与差分电荷

CHGCAR、LOCPOT、ELFCAR、cube 文件都是在实空间网格上存储的标量场。VASP 的 CHGCAR 中存的是 $\rho(\mathbf r)\times V_\text{cell}$，不同读取程序会自行换算。差分电荷常用

$$\Delta\rho(\mathbf r)=\rho_{AB}(\mathbf r)-\rho_{A}(\mathbf r)-\rho_{B}(\mathbf r)$$

三项必须在**同一晶胞、同一原子坐标、同一 FFT 网格**上计算，否则逐点相减没有意义。等值面就是 $\rho(\mathbf r)=c$ 的曲面，$c$ 的选择会改变图像观感，必须写在图注里。

### 2.5 能量零点

能带和 DOS 的横/纵坐标几乎总是相对某个参考能量：费米能级 $E_\mathrm F$（VASPKIT、py4vasp 默认）或价带顶 $E_\mathrm{VBM}$（sumo 对半导体默认）。PyProcar 6.5 对 VASP 数据**不会自动平移**：不传 `fermi=` 时只给出警告，纵轴是绝对能量。注意两点：

1. ICHARG = 11 的能带计算打印的 E-fermi **不可靠**（路径上的 k 点不能代表整个布里渊区），参考能量应取自洽或 DOS 计算的值；
2. 绝缘体的 $E_\mathrm F$ 落在带隙中的位置取决于展宽方法和程序实现，不同计算之间比较时用 VBM 对齐或真空能级对齐更稳妥。

---

## 3. 准备计算：能带、DOS、HSE 与光学的 VASP 输入

工作流（示例目录 [05-band-dos-plotting](../examples/05-band-dos-plotting/README.md)）：

```text
scf/   ISMEAR=0, 8x8x8 Gamma          -> CHGCAR, E-fermi(参考)
 ├── band/  ICHARG=11, line-mode KPOINTS, LORBIT=11  -> vasprun.xml, PROCAR, EIGENVAL
 └── dos/   ICHARG=11, 16x16x16, ISMEAR=-5, NEDOS=3001, LORBIT=11 -> DOSCAR, vasprun.xml
hse_band/   LHFCALC, IBZKPT(带权重) + 路径点(零权重)，自洽               -> vasprun.xml
```

### 3.1 能带：INCAR_band + KPOINTS_band_linemode

```text
SYSTEM  = GaAs PBE band structure along a line-mode path (non-SCF)

PREC    = Accurate
ENCUT   = 400        # MUST equal the SCF value
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.02
LREAL   = .FALSE.

ICHARG  = 11         # read CHGCAR, keep the density fixed (non-self-consistent)
NBANDS  = 24         # 18 valence electrons -> 9 occupied bands; add empty bands
LORBIT  = 11         # lm-decomposed PROCAR / vasprun.xml projections (fat bands)

NSW     = 0
IBRION  = -1
LCHARG  = .FALSE.    # never overwrite the SCF CHGCAR
LWAVE   = .FALSE.    # set .TRUE. only if you need WAVECAR (unfolding, wavefunctions)
```

```text
GaAs fcc path L-G-X-U|K-G (Setyawan-Curtarolo labels)
40              ! points per segment, both end points included
Line-mode
Reciprocal
  0.500  0.500  0.500   L
  0.000  0.000  0.000   GAMMA

  0.000  0.000  0.000   GAMMA
  0.500  0.000  0.500   X

  0.500  0.000  0.500   X
  0.625  0.250  0.625   U

  0.375  0.375  0.750   K
  0.000  0.000  0.000   GAMMA
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ICHARG` | 11 | 读入自洽 CHGCAR 并保持不变；只适用于 LDA/GGA/GGA+U |
| `ENCUT`、`PREC`、POTCAR | 与 SCF 相同 | 否则读入的密度与当前基组不匹配 |
| `NBANDS` | 占据带数的 1.5–3 倍 | 想看导带就要加空带；默认值通常只多一点点 |
| `LORBIT` | 11 | 输出 $lm$ 分解投影，fat band 和 PDOS 都需要 |
| `LCHARG` | .FALSE. | 防止非自洽计算覆盖 SCF 的 CHGCAR |
| KPOINTS 第 2 行 | 20–60 | 每段的点数（Line-mode），段内等分；`U|K` 这种断点写成两段 |
| 点标签 | 第 4 列文字 | VASPKIT、sumo、pymatgen 都从这里读标签；漏写时 VASPKIT 的 KLABELS 会出现 “Undefined” |

生成路径的工具：VASPKIT 303（体相）/302（二维）、`sumo-kgen`、SeeK-path 网站、本仓库的 [kpath_gen.py](../scripts/kpath_gen.py)（`python kpath_gen.py -h` 查看用法）。

### 3.2 DOS：INCAR_dos

```text
SYSTEM  = GaAs PBE DOS/PDOS on a dense mesh (non-SCF)

PREC    = Accurate
ENCUT   = 400        # MUST equal the SCF value
EDIFF   = 1E-7
LREAL   = .FALSE.

ICHARG  = 11         # fixed SCF density
ISMEAR  = -5         # tetrahedron method with Bloechl corrections
LORBIT  = 11         # site- and lm-projected DOS in DOSCAR / vasprun.xml
NEDOS   = 3001       # energy grid points between EMIN and EMAX
EMIN    = -20        # absolute energies in eV (not relative to E-fermi)
EMAX    = 15
NBANDS  = 24

NSW     = 0
IBRION  = -1
LCHARG  = .FALSE.
LWAVE   = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `ISMEAR` | −5（半导体/绝缘体 DOS） | 四面体法 + Blöchl 修正；网格要 Γ 中心且足够密，k 点太少（< 4）会报错停止；金属也可用 −5 画 DOS，但不要用它算力 |
| `NEDOS` | 2000–5000 | 能量格点数；默认 301 太粗，画出来是折线 |
| `EMIN`/`EMAX` | 覆盖感兴趣范围 | 绝对能量，不是相对 E-fermi；不设时 VASP 自动选择 |
| `LORBIT` | 11 | 10 只分到 s/p/d，11 分到 $p_x,p_y,p_z,d_{xy}\dots$ |
| KPOINTS | 比 SCF 密 1.5–2 倍 | 本例 16×16×16；DOS 的平滑程度首先取决于 k 网格 |

### 3.3 HSE06 能带：零权重 KPOINTS

杂化泛函不能用 `ICHARG = 11`。传统做法是：把自洽网格的不可约 k 点（带权重，来自 IBZKPT）和路径点（权重 0）写在同一个显式 KPOINTS 中，做一次自洽计算。示例 `KPOINTS_hse_zero_weight` 的开头与结尾：

```text
GaAs HSE: 6x6x6 Gamma IBZ (16 weighted points) + 21 zero-weight points L-G-X
37
Reciprocal lattice
    0.00000000    0.00000000    0.00000000     1
    0.00000000    0.00000000    0.16666667     8
    ...                                              (16 weighted points in total)
    0.50000000    0.50000000    0.50000000     0    L
    0.45000000    0.45000000    0.45000000     0
    ...
    0.50000000    0.00000000    0.50000000     0    X
```

| 项目 | 说明 |
|---|---|
| 第 2 行 | k 点总数 = 带权重点数 + 零权重点数 |
| 带权重部分 | 必须与该结构、该网格、该对称性下 VASP 自己写出的 IBZKPT 等价（权重之和 = 网格总点数）；自己的体系**不要手工拼凑**，用 VASPKIT 251、`sumo-kgen --hybrid`（需要当前目录中的 IBZKPT）或先跑一次 PBE 取 IBZKPT |
| 零权重部分 | 路径点；点数越多越贵（每个点都参与 Fock 交换的计算），HSE 一般取较稀的路径（VASPKIT 251 提示路径分辨率 0.04–0.06） |
| 太贵怎么办 | `sumo-kgen --hybrid --split 30 --folders` 把路径拆到多个 `split-0*` 目录分别计算，`sumo-bandplot` 在上级目录自动拼接 |
| VASP ≥ 6.3 的替代 | `KPOINTS_OPT` 一步法，见 [VASP 第6章](../../VASP/tutorials/06-hybrid-metagga.md) |

作图：sumo-bandplot（自动去掉带权重点）、py4vasp、VASPKIT 252（投影能带用 253–257）。本例的 matplotlib 脚本只处理 line-mode 运行。

> [!WARNING]
> VASPKIT 252–257 只支持由 **VASPKIT 251 生成**的 KPOINTS：251 在第一行写入 `Parameters to Generate KPOINTS (Don't Edit This Line): ...`，252 靠这一行识别哪些是零权重路径点（见 [VASPKIT 教程](https://vaspkit.com/tutorials.html) 的杂化泛函能带一节）。本例的 `KPOINTS_hse_zero_weight` 是按对称性手工写出的（16 个带权重点的权重已按 6×6×6 网格、$T_d$ 点群加时间反演核对过，和为 216），第一行不是这种格式，所以只能用 sumo 或 py4vasp 作图。想用 VASPKIT 后处理，就用 251 重新生成 KPOINTS（303 → 251）后再算。

### 3.4 线性光学：INCAR_optics_LOPTICS

VASPKIT 711 和 sumo-optplot 读取 LOPTICS 计算写进 vasprun.xml 的介电函数。两步做法：先普通 SCF（`LWAVE = .TRUE.`），再读 WAVECAR、加大 `NBANDS`、`ALGO = Exact`、`LOPTICS = .TRUE.`、`CSHIFT = 0.1`、`NEDOS = 2000`（完整文件见 [05-vaspkit-batch/INCAR_optics_LOPTICS](../examples/05-vaspkit-batch/INCAR_optics_LOPTICS)）。物理背景、k 点收敛和 GW/BSE 见 [VASP 第14章](../../VASP/tutorials/14-optics-dielectric.md)。

### 3.5 运行与检查

```bash
cd Common/examples/05-band-dos-plotting
cat $VASP_PP_PATH/potpaw_PBE/Ga_d/POTCAR $VASP_PP_PATH/potpaw_PBE/As/POTCAR > POTCAR
sbatch job_band_dos.slurm                 # RUN_HSE=1 追加 HSE06 能带

grep "E-fermi" scf/OUTCAR dos/OUTCAR      # 参考能量：用这两个值，不用 band/OUTCAR 的
grep "General timing" band/OUTCAR dos/OUTCAR   # 两个非自洽计算都正常结束
grep "NKPTS" band/OUTCAR                  # 路径点数 = 4 段 x 40 = 160
```

成功标准：band、dos 两个 OUTCAR 都有 timing 信息；`band/OUTCAR` 中 `NKPTS = 160`；`dos/OUTCAR` 中的 E-fermi 与 `scf/OUTCAR` 相差很小（四面体法与高斯展宽在带隙中的位置可以不同，绝缘体相差零点几 eV 也属正常）。

---

## 4. VASPKIT 1.5.1：交互菜单、命令行与批处理

### 4.1 版本与安装配置

VASPKIT 只提供编译好的可执行文件（Linux/macOS/Windows），1.5.1（2024-01）仍是最新的公开版本，另有收费的 VASPKIT Pro。老版中文手册（0.73，GitHub 上的 Markdown 版本）只能参考，**功能编号和提示顺序以 [vaspkit.com](https://vaspkit.com/features.html) 的 1.5 文档为准**。

```bash
tar -zxvf vaspkit.1.5.1.linux.x64.tar.gz       # 文件名以下载页为准
cd vaspkit.1.5.1
cp -f how_to_set_environment_variable ~/.vaspkit
vi ~/.vaspkit                                   # 至少改 PBE_PATH；作图要设 PYTHON_BIN
echo 'export PATH=$HOME/vaspkit.1.5.1/bin:$PATH' >> ~/.bashrc && source ~/.bashrc
```

| `~/.vaspkit` 条目 | 推荐值 | 说明 |
|---|---|---|
| `PBE_PATH` | 解压的 PBE 赝势库目录 | 103/104 生成 POTCAR 时使用，通常只需改这一项 |
| `RECOMMENDED_POTCAR` | `.TRUE.` | 自动选用 VASP 推荐的后缀（`_pv`、`_sv`、`_d`） |
| `SET_FERMI_ENERGY_ZERO` | `.TRUE.` | BAND.dat、DOS 数据中把 E-fermi 平移到 0 |
| `PYTHON_BIN` | `dft` 环境的 python 路径 | 需要 Python 3 + numpy、scipy、matplotlib；指向系统 Python 时作图会报找不到模块 |
| `AUTO_PLOT` | `.TRUE.`（≥ 1.3.0） | 后处理时自动出图；首次运行生成 `PLOT.in`，可修改后重跑 |

1.5.0 起，若 `~/.vaspkit` 不存在会读取 `/etc/vaspkit`（集群管理员可以提供全局配置）。完整模板见 [01-env-setup/vaspkit_rc.example](../examples/01-env-setup/vaspkit_rc.example)；环境自检可用 [env_check.sh](../scripts/env_check.sh)（`bash env_check.sh -h` 查看用法）。

### 4.2 四种运行方式

| 方式 | 示例 | 适用 |
|---|---|---|
| 交互菜单 | `vaspkit`，输入 `21` → `211` | 第一次使用、确认提示顺序 |
| 命令行参数 | `vaspkit -task 102 -kpr 0.04` | 官方说明只有部分功能支持，`vaspkit -help` 查看 |
| 管道 | `echo -e "102\n2\n0.04" \| vaspkit` | 所有功能，按提示顺序依次给出答案 |
| 应答文件 | `vaspkit < cmd_211_band.in` | 答案较多时清晰，可纳入版本管理 |

管道中可以直接输入三位功能号（如 `211`），不必先输一级菜单号。`echo -e '102\n2\n0.04'` 的含义是：功能 102 → 选 2（Γ 中心网格）→ k 点分辨率 0.04（单位 2π/Å）。每个方向的网格数约为 $N_i\approx|\mathbf b_i|/(2\pi\cdot\text{kpr})$（$|\mathbf b_i|$ 含 2π 因子，取整且不小于 1），分辨率数值越小网格越密。

> [!NOTE]
> 打开 `AUTO_PLOT` 后，部分版本在写完数据文件后还会多问一个作图选项（例如 211 的 “If you want use the default setting, type 0”）。应答文件里少了这一行时，程序读到文件结尾可能报错或停住。批处理前先交互跑一次：多出来的提示就在应答文件末尾补一行 `0`；也可以在批处理时临时关掉 `AUTO_PLOT`。

### 4.3 本章用到的功能编号

| 编号 | 功能（1.5 文档名称） | 读取 | 写出 / 屏幕输出 |
|---|---|---|---|
| 101 / 102 / 103 | INCAR 模板 / SCF 用 KPOINTS / 默认 POTCAR | POSCAR | INCAR、KPOINTS、POTCAR |
| 303 | 体相 k 路径（302 二维、301 一维） | POSCAR | KPATH.in、PRIMCELL.vasp、HIGH_SYMMETRY_POINTS |
| 305 | phonopy 用的 k 路径 | POSCAR | phonopy 可用的路径设置（以实际输出为准） |
| 251 / 252 | 杂化泛函能带的 KPOINTS / 后处理（253–257 为投影能带） | POSCAR + KPATH.in / 计算输出 + **251 写出的** KPOINTS | 零权重 KPOINTS / BAND.dat 等 |
| 211 / 213 | 能带 / 分元素投影能带 | EIGENVAL、KPOINTS、POSCAR、DOSCAR（读 E-fermi）；213 另需 PROCAR（`LORBIT = 10/11`） | BAND.dat、BAND_REFORMATTED.dat、KLINES.dat、KLABELS、BAND_GAP；213 写 PBAND_<El>.dat |
| 111 / 113 / 115 | 总 DOS / 分元素 PDOS / 指定原子和轨道求和 | DOSCAR、POSCAR | TDOS.dat、PDOS_<El>.dat、PDOS_USER.dat |
| 911 | 带隙 | EIGENVAL | 屏幕输出 VBM/CBM/带隙 |
| 711 | 体相线性光学谱（710 二维） | vasprun.xml | 吸收系数、折射率、消光系数、反射率、能量损失等 .dat |
| 201 | 能量–应变法弹性常数（200 为应力–应变法） | 各应变目录的计算 | 弹性张量与力学性质，见 [VASP 第13章](../../VASP/tutorials/13-elastic.md) |
| 426 | LOCPOT 一维平面平均势 | LOCPOT | 平面平均数据（1.5 教程中为 `POTPAVG.dat`，旧版为 `PLANAR_AVERAGE.dat`）与屏幕上的 `Vacuum-Level (eV)`，功函数 = E_vac − E_F |
| 501 / 502 | 吸附物 / 气体分子的热力学校正 | OUTCAR（频率）、CONTCAR | ZPE、U(T)、H(T)、G(T) 校正与熵 |
| 314 | 差分电荷 | 2 或 3 个 CHGCAR 路径 | CHGDIFF.vasp（可用 VESTA 打开） |

> [!WARNING]
> 有些教程把 305 写成“杂化泛函 k 路径”——在 1.5 文档中 305 是 phonopy 用的 k 路径，杂化泛函能带用 251/252。另外，501 是**吸附物**（只算振动，频率低于 50 cm⁻¹ 的按 50 cm⁻¹ 处理，不含 PV 项），502 是**气体分子**（含平动、转动、振动和电子贡献，需要输入温度、压强和自旋多重度）。编号若与你的版本不符，以 `vaspkit` 菜单为准。

### 4.4 批处理脚本

示例 [vaspkit_batch.sh](../examples/05-vaspkit-batch/vaspkit_batch.sh) 对一组目录执行同一功能，先检查必需文件，屏幕输出保存为 `vaspkit_<task>.log`，再从数据文件中提取结果：

```bash
cd Common/examples/05-vaspkit-batch
bash vaspkit_batch.sh band  ../05-band-dos-plotting/band            # 211
bash vaspkit_batch.sh dos   ../05-band-dos-plotting/dos             # 111 + 113
bash vaspkit_batch.sh gap   strain_*/band                           # 911
KPR=0.03 bash vaspkit_batch.sh kpoints relax_*/                     # 102, Gamma, 0.03
TEMP=298.15 PRES=1 MULT=1 bash vaspkit_batch.sh thermo-gas H2O_freq CO2_freq   # 502
MULT=3 bash vaspkit_batch.sh thermo-gas O2_freq                     # O2 为三重态
bash vaspkit_batch.sh workfunc slab_*/                              # 426, z 方向
```

脚本内部就是这样的管道：

```bash
( cd "$d" && printf '502\n%s\n%s\n%s\n' "$TEMP" "$PRES" "$MULT" | vaspkit > vaspkit_thermo-gas.log 2>&1 )
grep -i -m1 "Thermal correction to G" "$d/vaspkit_thermo-gas.log"
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| 每个计算一个目录 | 必须 | VASPKIT 写固定文件名，会互相覆盖 |
| 结果提取 | 读 BAND_GAP、PDOS_*.dat 等文件或按关键字 grep | 按行号截屏幕输出（`sed -n '40p'`）在升级后容易失效 |
| 211/213 的 E-fermi | 能带计算**结束后**把它自己的 DOSCAR 改名备份，再把 SCF（或 DOS）的 DOSCAR 复制进来；或关掉 `SET_FERMI_ENERGY_ZERO` 自己平移 | 211/213 从 DOSCAR 读 E-fermi；能带计算会写出自己的 DOSCAR，其中的 E-fermi 不可靠 |
| 501/502 的输入 | 充分弛豫后的频率计算（[INCAR_freq_thermo](../examples/05-vaspkit-batch/INCAR_freq_thermo)） | 吸附物计算中固定衬底原子；虚频要先消除 |

（批处理的输出格式示意见示例 README；没有真实运行结果时不要把示意数值当作参考值。）

### 4.5 输出文件给 Origin/Excel

`BAND_REFORMATTED.dat`：第一列为路径长度（Å⁻¹），其后每列一条能带，可直接导入 Origin；高对称点位置在 `KLABELS`。`TDOS.dat`、`PDOS_*.dat` 第一列为能量。自旋极化计算会分别给出自旋向上/向下的数据（文件名或列标题带 UP/DW 字样，以实际输出为准）。

---

## 5. sumo：命令行出版级能带与 DOS

[sumo](https://github.com/SMTG-Bham/sumo)（3.0.0，MIT）基于 pymatgen，一条命令输出排版良好的 PDF，适合直接放进论文。安装：`pip install sumo`（3.0.0 要求 Python ≥ 3.10）。

### 5.1 sumo-kgen：生成路径

```bash
sumo-kgen -p POSCAR                       # 默认 Bradley-Cracknell 约定 -> KPOINTS_band
sumo-kgen -p POSCAR --pymatgen            # Setyawan-Curtarolo（与 pymatgen/MP 一致）
sumo-kgen -p POSCAR --seekpath            # SeeK-path（Hinuma）约定
sumo-kgen --hybrid --split 30 --folders   # HSE：当前目录需有 IBZKPT，路径点零权重追加，按每 30 点拆分
```

| 选项 | 说明 |
|---|---|
| `--hybrid` / `-y` | 把路径点以零权重追加到 IBZKPT 后面；缺少 IBZKPT 时报 `ERROR: Hybrid specified but no IBZKPT file found.` |
| `--split N` | 每个文件 N 个路径点，用于 HSE 分批计算 |
| `--folders` / `-f` | 自动建 `split-01`、`split-02`…目录并复制输入文件 |
| `--density` | 路径点密度（默认 60） |
| `--symprec` | 对称性容差（默认 0.01） |

> [!IMPORTANT]
> 如果输入结构不是所选约定下的标准原胞，sumo 会把标准原胞写成 `POSCAR_prim` 并在屏幕提示——此时能带计算必须使用 `POSCAR_prim`（以及它对应的 CHGCAR，即自洽也要用这个原胞重算），否则路径坐标与晶胞不匹配，能带完全错误。

### 5.2 sumo-bandplot 与 sumo-dosplot

```bash
cd band
sumo-bandplot --band-edges --ymin -6 --ymax 4 --format pdf                     # band.pdf + band.dat
sumo-bandplot --project Ga.s,Ga.p,As.p --mode rgb --ymin -6 --ymax 4 --prefix proj
sumo-bandplot --dos ../dos/vasprun.xml --elements Ga.s.p.d,As.s.p --prefix banddos
cd ../dos
sumo-dosplot --elements Ga.s.p.d,As.s.p --gaussian 0.05 --xmin -6 --xmax 4 --prefix gaas
```

| 选项 | 说明 |
|---|---|
| `--project Ga.s,Ga.p,As` | 投影能带；元素后用 `.` 接轨道，同一元素可多次出现；`--mode rgb`（默认，最多 3 组，颜色混合）或 `stacked` |
| `--dos <vasprun.xml>` | 在右侧附 DOS 面板，DOS 相关选项（`--elements`、`--orbitals`、`--atoms`、`--gaussian`）与 sumo-dosplot 相同 |
| `--elements "C.s.p,O"` | 只画指定元素/轨道的 PDOS |
| `--orbitals Ru.d` | 把指定轨道拆成 $lm$ 分量 |
| `--spin up` / `down` | 自旋极化计算只画一个通道 |
| `--zero-energy E` | 手动指定零点；**默认：半导体取 VBM，金属取 E-fermi** |
| `--style`、`--config` | 自定义 matplotlib 样式与线条颜色 |
| `--format`、`--dpi`、`--width`、`--height` | 输出格式（默认 pdf）与尺寸 |

sumo-bandplot 会在当前目录自动寻找 vasprun.xml(.gz)；若存在 `split-0*` 子目录则自动拼接。同时写出 `band.dat`（或 `<prefix>_band.dat`），sumo-dosplot 写出 `total_dos.dat` 和 `<El>_dos.dat`（加前缀时为 `<prefix>_total_dos.dat` 等），便于在 Origin 中重画。

### 5.3 sumo-bandstats 与 sumo-optplot

```bash
cd band && sumo-bandstats            # VBM/CBM 位置、带隙、沿各方向的有效质量（写 sumo-bandstats.log）
sumo-bandstats -s 5                  # 拟合用 5 个采样点（默认 3）
sumo-bandstats -n                    # 改用非抛物线模型拟合
sumo-optplot                         # 在 LOPTICS 计算目录中画吸收谱（读 vasprun.xml）
```

sumo-bandstats **默认用抛物线模型**拟合，`-n`（`--nonparabolic`）改用非抛物线模型（sumo 3.0 源码中该选项是 `store_false`，传给 `parabolic` 参数）。结果写进当前目录的 `sumo-bandstats.log`，屏幕上的信息走的是 stderr，所以用 `tee` 保存时要写 `sumo-bandstats 2>&1 | tee ...`。有效质量要在路径点足够密的能带上算，结果与 VASPKIT 912、py4vasp 交叉验证。

完整脚本：[sumo_plot.sh](../examples/05-band-dos-plotting/sumo_plot.sh)。

---

## 6. PyProcar 6：投影能带、PDOS 与费米面

[PyProcar](https://github.com/romerogroup/pyprocar)（GPL-3.0；本文写作时 PyPI 版本为 6.5.0）擅长投影能带、自旋纹理和三维费米面。**v6 起 API 由“传 PROCAR 文件路径”改为“传代码名 + 目录”**，旧教程中的 `pyprocar.bandsplot('PROCAR', outcar='OUTCAR')` 在新版中不再可用。

```python
import pyprocar

ef = 3.05   # E-fermi of the SCF/DOS run (example value; use your own grep E-fermi)

# plain bands
pyprocar.bandsplot(code="vasp", dirname="band", mode="plain", fermi=ef,
                   elimit=[-6, 4], savefig="bands_plain.png", show=False)

# projected bands: colour = weight of As p (atom index 1, orbitals 1-3 = py, pz, px)
pyprocar.bandsplot(code="vasp", dirname="band", mode="parametric", fermi=ef,
                   atoms=[1], orbitals=[1, 2, 3], elimit=[-6, 4],
                   savefig="bands_As_p.png", show=False)

# PDOS stacked by species
pyprocar.dosplot(code="vasp", dirname="dos", mode="stack_species", fermi=ef,
                 elimit=[-6, 4], savefig="dos_species.png", show=False)

# bands + DOS in one figure
pyprocar.bandsdosplot(code="vasp",
                      bands_settings={"dirname": "band", "mode": "plain", "fermi": ef},
                      dos_settings={"dirname": "dos", "mode": "plain", "fermi": ef},
                      elimit=[-6, 4], savefig="bands_dos.png", show=False)
```

| 参数 | 说明 |
|---|---|
| `code`、`dirname` | 程序名（`"vasp"`、`"qe"` 等）与计算目录；VASP 目录中需要 PROCAR（`LORBIT = 11`）、OUTCAR、KPOINTS、POSCAR、vasprun.xml |
| `mode`（能带） | `plain`、`parametric`（颜色表示权重）、`scatter`（点大小）、`atomic`、`overlay`、`overlay_species`、`overlay_orbitals`、`ipr` |
| `mode`（DOS） | `plain`、`parametric`、`parametric_line`、`stack`、`stack_species`、`stack_orbitals`、`overlay*` |
| `atoms`、`orbitals` | **从 0 开始**的原子与轨道索引；VASP 轨道顺序 0 s，1–3 $p_y,p_z,p_x$，4–8 $d_{xy},d_{yz},d_{z^2},d_{xz},d_{x^2-y^2}$ |
| `spins` | 自旋极化：`[0]` 向上、`[1]` 向下；非共线（SOC）计算中含义变为总投影与各自旋分量，见官方文档 |
| `fermi` | **务必传自洽/DOS 计算的 E-fermi**；6.5 中 VASP 数据不传 `fermi` 时能带不做平移（只打印 “`fermi` is not set!” 警告），纵轴是绝对能量 |
| `elimit`、`savefig`、`show` | 能量窗口、保存文件、是否弹窗（集群上设 `show=False`） |

三维费米面（金属，需在均匀密网格上做非自洽计算；渲染依赖 pyvista）。v5 教程中的 `pyprocar.fermi3D(procar, outcar, ...)` 在 v6 中已被 `FermiHandler` 类取代：

```python
fs = pyprocar.FermiHandler(code="vasp", dirname="fermi_mesh", fermi=ef)
fs.plot_fermi_surface(mode="plain", show=False, save_2d="fermi3d.png")
```

完整脚本：[pyprocar_fatband.py](../examples/05-band-dos-plotting/pyprocar_fatband.py)（`python pyprocar_fatband.py -h` 查看选项）。

---

## 7. py4vasp：读取 vaspout.h5

[py4vasp](https://vasp.at/py4vasp/latest)（VASP 官方，Apache-2.0，0.11.x；拆分为 py4vasp-core、py4vasp 与 vasp-backend 三个分发包）读取 HDF5 格式的 `vaspout.h5`，能画能带、DOS、声子、介电函数、STM 等，也支持 VASP 6.3+ 的 `KPOINTS_OPT`。前提是 **VASP 编译时打开 HDF5**：在 `makefile.include` 的 `CPP_OPTIONS` 中加 `-DVASP_HDF5`，并按官方模板链接 HDF5 库；之后每次运行都会写出 `vaspout.h5`。VASP 6.5 起 HDF5 输出进一步增强（可在运行中访问，见 [vaspout.h5 Wiki](https://vasp.at/wiki/Vaspout.h5) 与 `LSYNCH5`）。

```python
import py4vasp

dos  = py4vasp.Calculation.from_path("dos")
band = py4vasp.Calculation.from_path("band")
ef = dos.dos.read()["fermi_energy"]          # reference from the dense-mesh run

band.band.plot(fermi_energy=ef)                         # plain bands (interactive plotly figure)
band.band.plot("Ga(s), As(p)", fermi_energy=ef)         # fat bands, needs LORBIT = 11
dos.dos.plot("Ga, As(s), As(p)")                        # projected DOS
band.band.to_csv(filename="band.csv")                   # data for Origin / matplotlib
# band.band.to_image(filename="band.pdf")               # static image, needs: pip install kaleido
```

| 写法 | 含义 |
|---|---|
| `"Ga(s), As(p)"` | 元素 + 轨道，逗号分隔多条曲线 |
| `"1(d)"` | 第 1 个原子（**从 1 开始**）的 d 轨道 |
| `"dxy + dxz + dyz"` | 多个轨道求和（如 t2g） |
| `"up"`、`"down"`、`"up(As(p))"` | 自旋极化计算的通道 |
| `"sigma_x"`、`"sigma_z"` | 非共线计算的自旋分量 |
| `"kpoints_opt"` | 画 KPOINTS_OPT 路径上的能带 |
| `.selections()` | 列出当前计算可用的选择 |

笔记本：[py4vasp_demo.ipynb](../examples/05-band-dos-plotting/py4vasp_demo.ipynb)。py4vasp 的图是 plotly 交互图，导出静态图需要 kaleido；论文用图建议导出 CSV 后用第 9 节的样式文件重画，保证全文风格统一。

---

## 8. 自旋极化、SOC、反折叠与费米面

### 8.1 自旋极化（ISPIN = 2）

| 要点 | 做法 |
|---|---|
| 参考能量 | 两个自旋通道**共用同一个 E-fermi**，不要分别平移 |
| 能带 | 向上实线、向下虚线（或两种颜色），图例写明；也可左右两个面板 |
| DOS | 向下通道取负值画在下方/左侧（镜像），横轴 0 处加细线 |
| sumo | 默认两个通道都画；`--spin up` / `--spin down` 只画一个 |
| PyProcar | `spins=[0]` / `spins=[1]` |
| py4vasp | 选择 `"up"`、`"down"`，或组合 `"up(Fe(d))"` |
| matplotlib 脚本 | `plot_band_dos_matplotlib.py` 自动识别：向下能带虚线、DOS 镜像 |
| 磁矩检查 | 作图前先 `grep -A10 "magnetization (x)" OUTCAR` 确认磁性状态是你想要的那个（铁磁/反铁磁） |

### 8.2 自旋轨道耦合（LSORBIT = .TRUE.，vasp_ncl）

- 非共线计算没有“上/下”两个通道，每条能带只容纳一个电子，**能带数和占据带数都翻倍**（`NBANDS` 需相应加倍）；
- PROCAR 中每个态给出总投影以及 $m_x, m_y, m_z$ 三个自旋分量，可以画自旋纹理（VASPKIT 651–653、PyProcar、py4vasp 的 `sigma_x/z`、IFermi 的 spin 性质）；
- SOC 与非 SOC 能带对比时，用同一路径、同一参考（VBM 或 E-fermi），并在图注写明；
- 能带反转或 Rashba 劈裂这类细节对路径点密度敏感，在关心的 k 点附近加密路径。

SOC 的计算细节见 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)。

### 8.3 能带反折叠（easyunfold）与费米面（IFermi）

缺陷、合金和表面超胞的能带折叠成一团，需要反折叠回原胞布里渊区：

```bash
pip install easyunfold
easyunfold generate primitive/POSCAR supercell/POSCAR primitive/KPOINTS_band --matrix "2 2 2"
# 用生成的 k 点做超胞非自洽计算（LWAVE = .TRUE.），然后：
easyunfold unfold calculate WAVECAR
easyunfold unfold plot
```

金属费米面：在均匀密网格上做非自洽计算，`ifermi info` 查看费米面信息、`ifermi plot` 画交互式三维图（可按费米速度、自旋着色）。两者的完整流程、HSE 反折叠和费米速度见 [VASP 第17章](../../VASP/tutorials/17-advanced-electronic.md)。

---

## 9. 自己用 matplotlib 画：mplstyle 与出版级规范

脚本工具很方便，但一篇论文里的所有图最好风格统一（字体、字号、线宽、配色）。做法是写一个样式文件，所有作图脚本都 `plt.style.use()` 它。

### 9.1 样式文件 dft_paper.mplstyle（节选）

```text
figure.figsize      : 3.35, 2.6        # single column ~ 85 mm
savefig.dpi         : 600              # raster formats only
savefig.bbox        : tight
pdf.fonttype        : 42               # embed TrueType fonts (editable text)
svg.fonttype        : none
font.family         : sans-serif
font.sans-serif     : Arial, Helvetica, Liberation Sans, DejaVu Sans
font.size           : 8
mathtext.fontset    : dejavusans
axes.linewidth      : 0.6
axes.prop_cycle     : cycler('color', ['0072B2', 'D55E00', '009E73', 'CC79A7', 'E69F00', '56B4E9', '000000'])
xtick.direction     : in
ytick.direction     : in
xtick.top           : True
ytick.right         : True
lines.linewidth     : 1.0
legend.frameon      : False
image.cmap          : viridis
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| `figure.figsize` | 单栏 3.3–3.5 in（85–90 mm），双栏 6.9–7.1 in（175–180 mm） | 按期刊要求的**最终尺寸**作图，避免排版时缩放导致字号不一 |
| `font.size` | 7–9 pt | 缩放到最终尺寸后不小于 6 pt |
| `pdf.fonttype` | 42 | 嵌入 TrueType 字体；默认的 Type 3 字体在部分投稿系统和 Illustrator 中显示异常 |
| `svg.fonttype` | none | SVG 中保留文字对象，便于后期编辑 |
| `savefig.dpi` | 600（线条图）| 只对 PNG/TIFF 等位图有效；矢量格式（PDF/SVG/EPS）与 dpi 无关 |
| 颜色循环 | Okabe–Ito 色盲友好配色 | mplstyle 中十六进制颜色**不写 `#`**（`#` 是注释符） |
| `image.cmap` | viridis / cividis（顺序）、RdBu（发散） | 不要用 jet/rainbow：亮度不单调、色盲不可读 |

### 9.2 脚本 plot_band_dos_matplotlib.py

[plot_band_dos_matplotlib.py](../examples/05-band-dos-plotting/plot_band_dos_matplotlib.py) 用 pymatgen 读 vasprun.xml，按分支（`bs.branches`）逐段画能带（避免在 U|K 这类断点处画出竖线），合并重合的高对称点标签，右侧画总 DOS 和分元素 PDOS，自动处理自旋极化：

```bash
cd Common/examples/05-band-dos-plotting
python plot_band_dos_matplotlib.py --band band --dos dos                       # E-fermi(DOS run) = 0
python plot_band_dos_matplotlib.py --band band --dos dos --project Ga,As --pdos-orbitals \
       --ylim -6 4 --formats pdf svg png --out gaas_band_dos
python plot_band_dos_matplotlib.py --band band --zero vbm --efermi 3.05        # 只画能带，VBM = 0
```

| 选项 | 说明 |
|---|---|
| `--efermi` | 手动给参考 E-fermi（取 SCF 的 OUTCAR）；不给时用 DOS 计算的值 |
| `--zero efermi\|vbm` | 能量零点 |
| `--project Ga,As` | 分元素 fat band（需要 `LORBIT = 11`，读取较慢） |
| `--pdos-orbitals` | PDOS 按 s/p/d/f 分开 |
| `--smear 0.05` | 对 DOS 额外做高斯展宽（eV） |
| `--formats pdf svg png` | 同时导出多种格式 |

本仓库 VASP 板块还有单独的能带/DOS 脚本：[plot_band.py](../../VASP/scripts/plot_band.py)、[plot_dos.py](../../VASP/scripts/plot_dos.py)。

### 9.3 出版级作图检查清单

1. **坐标轴**：能量轴写成 “$E-E_\mathrm F$ (eV)” 或 “$E-E_\mathrm{VBM}$ (eV)”，零点在图注中说明；DOS 单位写 states/eV（必要时注明是每晶胞还是每原子）。
2. **高对称点**：标签用正体，Γ 用希腊字母；断点写成 “U|K”。
3. **线宽与字号**：最终尺寸下线宽 ≥ 0.5 pt、字号 ≥ 6 pt；子图编号 (a)(b) 统一位置。
4. **颜色**：同一元素在全文所有图中颜色一致；不只靠颜色区分（线型、标记、直接标注）；灰度打印检查一次。
5. **格式**：期刊接受时优先 PDF/EPS/SVG 矢量图；位图线条图 ≥ 600 dpi，照片类（如 VESTA 渲染）≥ 300 dpi，以目标期刊作者指南为准。
6. **可复现**：把作图脚本、样式文件和导出的数据（band.dat、CSV）与原始计算一起存档；图注写明泛函、软件版本、展宽方法与参数。
7. **Origin 用户**：VASPKIT 的 BAND_REFORMATTED.dat/TDOS.dat、sumo 的 band.dat/total_dos.dat、py4vasp 的 CSV 都能直接导入；同样应遵守上面的尺寸、字体与配色规范。

---

## 10. 三维结构与体数据：VESTA

[VESTA](https://jp-minerals.org/vesta/en/)（免费）是晶体结构和体数据可视化最常用的工具。稳定版仍为 **3.5.8**（2021）；**3.90.x 开发版**（2026-08 的 3.90.9a）提供 Apple Silicon、Windows on Arm 和 Linux GTK3 原生构建，Mac M 系列芯片用户可以直接用开发版。

| 任务 | 操作 | 说明 |
|---|---|---|
| 打开结构/体数据 | File → Open：POSCAR/CONTCAR/cif、CHGCAR/PARCHG/ELFCAR/LOCPOT、cube、xsf | 自旋极化 CHGCAR 默认显示总密度；只看自旋密度先用 VASPKIT 312 提取 |
| 等值面 | Properties → Isosurfaces，设置 Isosurface level；有正负的量加两个等值面 | 电荷密度从 0.01–0.05 试起，差分电荷常用 0.001–0.005 量级；数值写进图注 |
| 截面 | Utilities → 2D Data Display，设定 (hkl) 与位置 | 画成键平面上的电荷或 ELF 等高线 |
| 差分电荷 | Edit → Edit Data → Volumetric Data → Import，选择减去 | AB、A、B 网格与坐标必须一致；等价方法 VASPKIT 314、[chgdiff.py](../../VASP/scripts/chgdiff.py) |
| 出图 | File → Export Raster Image（可放大）/ Export Vector Image | 平行投影、白底、去掉坐标轴 |

逐步操作见 [vesta_isosurface_howto.md](../examples/05-visualization-3d/vesta_isosurface_howto.md)。差分电荷的物理分析（Bader、平面平均电荷转移）见 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md)。

> [!TIP]
> 大体系的 CHGCAR 动辄几百 MB，不要在集群上开 VESTA 远程图形界面。先用 VASPKIT 310/316（截面、平面平均）或 pymatgen 在集群上处理，必要时只下载需要的 CHGCAR；结构和体数据的格式转换可用 [convert_any.py](../scripts/convert_any.py)（`python convert_any.py -h` 查看用法）。

---

## 11. 轨迹与动画：OVITO、VMD、Avogadro/Jmol、phonon website

### 11.1 OVITO：渲染 XDATCAR / AIMD 轨迹

[OVITO](https://www.ovito.org/) 3.17 的基础版免费，Python 模块可以 `pip install ovito`。注意两点：**3.16.1 起要求 Python ≥ 3.11**；官方提示**不要把 pip 版装进 Anaconda/conda 环境**（会与 conda 的 Qt 等库冲突），建议用普通 venv：

```bash
python3.12 -m venv ~/venv-ovito && source ~/venv-ovito/bin/activate
pip install ovito
python ovito_render_xdatcar.py XDATCAR --every-nth 5 --movie traj.gif --rdf-cutoff 6.0
```

[ovito_render_xdatcar.py](../examples/05-visualization-3d/ovito_render_xdatcar.py) 的核心几行：

```python
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

pipeline = import_file("XDATCAR")          # VASP format detected automatically
pipeline.add_to_scene()
vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(2, 1, -1))
vp.zoom_all(size=(800, 600))
vp.render_image(size=(800, 600), frame=0, filename="frame0.png", renderer=TachyonRenderer())
vp.render_anim(filename="traj.mp4", size=(800, 600), fps=10, every_nth=5, renderer=TachyonRenderer())
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| 渲染器 | `TachyonRenderer()` | 软件光线追踪，不需要 GPU 和显示器，适合集群 |
| 动画格式 | `.mp4`/`.avi`/`.mov`/`.gif`，或 `frame*.png` 图片序列 | 文件名扩展名决定格式 |
| `every_nth` | 5–20 | AIMD 动辄上万步，抽帧渲染 |
| RDF | `CoordinationAnalysisModifier(cutoff=6.0, number_of_bins=200)` | 脚本对所选帧取平均写出 `traj_rdf.dat` |

MSD、扩散系数等定量分析见 [VASP 第15章](../../VASP/tutorials/15-md-mlff.md)。

### 11.2 VMD + Multiwfn：轨道、ESP 表面与弱相互作用等值面

Multiwfn（2026 年起按发布日期命名版本，下载页最新为 2026.10.18，更新历史已有 2026.11.3 的条目，并改进了对 CP2K 2026.2 的兼容）负责计算格点数据（cube），[VMD](https://www.ks.uiuc.edu/Research/vmd/) 负责渲染：

1. Multiwfn 计算 → 导出 cube（轨道、密度、ESP、IGMH 的 `sl2r.cub` 与 `dg_inter.cub` 等）；
2. VMD 载入 cube，Representations 中加 Isosurface（Solid Surface），轨道的正负相位各一个；
3. 用 Multiwfn `examples` 目录中作者提供的 VMD 脚本（如 IGMH 的 `IGM_inter.vmd`）一键设置配色与等值面；
4. File → Render 选 Tachyon，输出高分辨率图片。

详细步骤与命令流见 [Gaussian 第16章](../../Gaussian/tutorials/16-postprocess-multiwfn.md) 第 8、9 节。

### 11.3 Avogadro / Jmol：分子振动动画

- **Avogadro 2**（免费、开源）：打开 Gaussian `.log`（或 `.fchk`）后在振动（Vibrational Modes）窗口选择频率、播放动画，可导出图片；也是 GaussView 的开源替代，可以生成 Gaussian/CP2K 等程序的输入；
- **Jmol**（Java，或浏览器中的 JSmol）：直接读取 Gaussian 输出中的振动模式，命令 `vibration on` 播放，适合嵌入网页；支持的其他程序格式见 Jmol 文件格式列表；
- **GaussView**：商业软件，Gaussian 用户的默认选择。

确认过渡态时，播放唯一虚频的振动方向是否对应反应坐标，是必做的一步（见 [Gaussian 第7章](../../Gaussian/tutorials/07-ts-irc.md)）。

### 11.4 phonon website：声子振动动画

[phonon website](https://henriquemiranda.github.io/phononwebsite/) 在浏览器中显示声子色散，点击任一模式即可看到原子振动动画，非常适合组会汇报和教学。用 phonopy 输出带本征矢的 band.yaml：

```text
DIM = 2 2 2
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
BAND = 0 0 0  1/2 0 1/2  5/8 1/4 5/8,  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2
BAND_LABELS = $\Gamma$ X U K $\Gamma$ L
BAND_POINTS = 51
EIGENVECTORS = .TRUE.
BAND_CONNECTION = .TRUE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `EIGENVECTORS` | `.TRUE.` | 必需：网站用本征矢做动画 |
| `BAND_CONNECTION` | `.TRUE.` | 按本征矢连接交叉的声子支 |
| `BAND_LABELS` | 与 `BAND` 的点一一对应 | 网站上显示的标签 |

```bash
phonopy --config band_phononwebsite.conf      # phonopy 4.x（2.x/3.x：phonopy band_phononwebsite.conf）
# 浏览器打开 phonon website，点 Choose files 上传 band.yaml
```

声子计算本身见 [VASP 第12章](../../VASP/tutorials/12-phonon.md)；文件：[band_phononwebsite.conf](../examples/05-visualization-3d/band_phononwebsite.conf)。

---

## 12. Multiwfn 批处理（跨程序）

Multiwfn 能读 Gaussian 的 fchk、各程序的 molden（ORCA、CP2K 等）、wfx/wfn、cube 等格式，因此是 Gaussian 与 CP2K 等程序共用的波函数分析工具（布居、键级、AIM、RDG/IGMH、分子表面静电势、空穴–电子分析等）。批处理的思路是：**交互做一遍，把每一步的输入按顺序写成命令流文件**，再对一批文件重放。

```bash
export OMP_STACKSIZE=200M      # 每个 OpenMP 线程的栈；大体系可以更大
ulimit -s unlimited            # 不设置时大体系常直接 Segmentation fault
formchk phenol.chk phenol.fchk # Multiwfn 不能读 .chk
Multiwfn phenol.fchk -silent < multiwfn_esp_surface_cmds.txt > phenol_esp.out
```

| 设置 | 推荐值 | 说明 |
|---|---|---|
| `-silent` | 批处理必加 | 不弹图形窗口，否则程序停在窗口等待 |
| `settings.ini` 的 `nthreads` | 与作业申请核数一致 | 当前目录有 settings.ini 时通常优先读取；在输出开头确认实际线程数 |
| `OMP_STACKSIZE` | 200M–1000M | 线程栈不足会崩溃 |
| 输入格式 | fchk 信息最全 | Gaussian 用 formchk 转换；CP2K 可输出 molden |
| 输出文件 | 每个输入单独目录 | Multiwfn 使用固定文件名（surfanalysis.txt、*.cub），批量时会互相覆盖 |

示例 [run_multiwfn_batch.sh](../examples/05-visualization-3d/run_multiwfn_batch.sh) 自动完成 chk→fchk 转换、私有 settings.ini（线程数）、每个文件一个 `<name>_mwfn/` 目录：

```bash
NTHREADS=8 bash run_multiwfn_batch.sh \
    ../../../Gaussian/examples/16-postprocess-multiwfn/multiwfn_esp_surface_cmds.txt  mols/*.fchk
```

经过测试的命令流（ESP 表面、IGMH、UV-Vis）见 [Gaussian 第16章示例](../../Gaussian/examples/16-postprocess-multiwfn/README.md)。命令流对 Multiwfn 版本敏感，升级后先交互核对一次菜单编号。中文社区的 Multiwfn 批处理工具（如 [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn)）可作参考。

---

## 13. Jupyter 与浏览器：nglview、weas-widget、matterviz

在 Jupyter 中交互查看结构，适合检查建模结果、在笔记本中记录分析过程：

```python
import nglview as nv
from ase.io import read

atoms = read("CONTCAR")
v = nv.show_ase(atoms.repeat((2, 2, 2)))
v.add_unitcell()
v                                         # rotate / zoom in the notebook

traj = read("XDATCAR", index="::10")      # every 10th frame
nv.show_asetraj(traj)                     # trajectory with a frame slider

from weas_widget import WeasWidget
viewer = WeasWidget()
viewer.from_ase(atoms)                    # ball / ball+stick / polyhedra, editable
viewer
edited = viewer.to_ase()                  # structure after editing in the widget
```

| 工具 | 安装 | 特点 |
|---|---|---|
| [nglview](https://github.com/nglviewer/nglview) | `conda install -c conda-forge nglview` | 成熟稳定，支持 ASE/pymatgen/MDAnalysis 对象和轨迹 |
| [weas-widget](https://github.com/superstar54/weas-widget) | `pip install weas-widget` | 可在界面中移动、删除、替换原子，并 `to_ase()` 取回；支持多面体、等值面、矢量（磁矩、声子） |
| [matterviz](https://github.com/janosh/matterviz) | 浏览器打开 <https://matterviz.janosh.dev>；或装 VS Code 扩展（支持 Remote SSH）、`pip install matterviz-jupyterlab` | 无需安装 Python；拖入 CIF/POSCAR/XDATCAR/CHGCAR/cube 即可查看，还能显示等值面与截面、能带+DOS、布里渊区、费米面、相图 |

笔记本：[nglview_demo.ipynb](../examples/05-visualization-3d/nglview_demo.ipynb)（没有 DFT 输出时用 ASE 自建结构演示）。

---

## 参数选择与收敛建议

| 项目 | 推荐 | 判据 / 说明 |
|---|---|---|
| 能带路径点数 | PBE 每段 20–60 点；HSE 路径分辨率 0.04–0.06（VASPKIT 251 提示值） | 有效质量和能带交叉处要更密 |
| 能带 `NBANDS` | 占据带数的 1.5–3 倍 | 导带完整显示到需要的能量 |
| DOS 网格 | 比 SCF 密 1.5–2 倍，Γ 中心 | 加密后 DOS 形状不再变化 |
| DOS 展宽 | 半导体 `ISMEAR = -5`；作图时可再加 0.02–0.1 eV 高斯展宽 | 展宽值写进图注 |
| `NEDOS` | 2000–5000 | 能量间隔 ≈ (EMAX − EMIN)/NEDOS ≤ 0.01–0.02 eV |
| 参考能量 | SCF/DOS 的 E-fermi 或 VBM | 不用能带计算的 E-fermi；不同体系比较时用真空能级或 VBM 对齐 |
| 有效质量拟合 | 3–5 个点，路径步长 ≤ 0.01 Å⁻¹ 量级 | 与第二种工具交叉验证 |
| 等值面数值 | 电荷 0.01–0.05；差分电荷 0.001–0.005；轨道 0.02–0.05 a.u. | 同一组对比图用同一数值 |
| 轨迹渲染 | 每 5–20 帧取 1 帧 | 动画 10–20 fps |
| 图幅 | 单栏 85–90 mm，双栏 175–180 mm | 最终尺寸下字号 ≥ 6 pt |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| 能带整体平移了零点几 eV，带隙中间穿过一条能带 | 用了 ICHARG = 11 能带计算自己打印的 E-fermi | 用 SCF/DOS 计算的 E-fermi（`--fermi`、`fermi_energy=`、`--efermi`） |
| sumo 和 VASPKIT 的图零点不同 | sumo 对半导体默认以 VBM 为零，VASPKIT 默认 E-fermi | 统一用 `--zero-energy` 或在 VASPKIT 中关掉 `SET_FERMI_ENERGY_ZERO` 后自行平移，图注写明 |
| 能带图在某个高对称点出现竖直连线 | 路径断点（如 U\|K）被当成连续路径画线 | 按分支分段画（`bs.branches`）；sumo/PyProcar 已自动处理 |
| VASPKIT 的 KLABELS 中出现 `Undefined` | KPOINTS 中没有写高对称点标签 | 在 line-mode KPOINTS 每行第 4 列写标签，或手工编辑 KLABELS |
| `vaspkit: /lib64/libc.so.6: version 'GLIBC_2.14' not found (required by vaspkit)` | 系统 glibc 太旧 | 换较新的系统/容器，或联系开发者索取低版本 glibc 编译的版本（官方 FAQ） |
| VASPKIT 自动作图报找不到 matplotlib | `PYTHON_BIN` 指向系统 Python | 改成 `dft` 环境中 `which python` 的路径 |
| 管道批处理结果错位或卡住 | 版本更新后提示顺序变了 | 交互运行一次核对，修改应答文件 |
| `ERROR: Hybrid specified but no IBZKPT file found.`（sumo-kgen） | `--hybrid` 需要当前目录中的 IBZKPT | 先跑一次（可很短的）自洽得到 IBZKPT |
| `ERROR: IBZKPT contains tetrahedron information.`（sumo-kgen） | 该 IBZKPT 来自 `ISMEAR = -5` 的计算 | 用 `ISMEAR = 0` 的计算重新生成 IBZKPT |
| HSE 能带乱跳或与 PBE 形状差很多 | 带权重部分与真实 IBZKPT 不一致，或误用了 ICHARG = 11 | 用工具重新生成零权重 KPOINTS；杂化泛函必须自洽 |
| PyProcar 旧代码 `bandsplot('PROCAR', outcar='OUTCAR')` 报参数或文件错误 | v6 改为 `code` + `dirname` 接口 | `pyprocar.bandsplot(code='vasp', dirname='band', mode='plain')`；PROCAR 需由 `LORBIT = 11` 生成 |
| PyProcar 打印 “`fermi` is not set! ... The plot did not shift the bands by the Fermi energy.”，图的零点在几 eV 之外 | 6.5 对 VASP 数据不再自动读取 E-fermi 平移 | 显式传 `fermi=`（SCF/DOS 计算的 E-fermi） |
| VASPKIT 252 读手写或 sumo 生成的零权重 KPOINTS 时出错、能带错乱 | 252–257 只认 251 写在 KPOINTS 第一行的生成参数 | 用 303 → 251 重新生成 KPOINTS 后再算，或改用 sumo-bandplot / py4vasp |
| fat band 全是零 / PyProcar 报找不到投影 | 能带计算没设 `LORBIT = 11` | 加 `LORBIT = 11` 重新做非自洽 |
| 各原子 PDOS 之和明显小于总 DOS | 投影只在 Wigner–Seitz 球内，间隙区未计入 | 正常现象；比较的是相对贡献 |
| DOS 是锯齿/尖峰 | k 网格太稀、`NEDOS` 太小 | 加密网格、`ISMEAR = -5`、`NEDOS ≥ 2000`，必要时作图时加小展宽 |
| pymatgen 读 vasprun.xml 报 `ParseError` | 作业中断，XML 不完整 | 确认计算正常结束；`Vasprun(..., exception_on_bad_xml=False)` 只用于抢救数据 |
| py4vasp 报找不到 `vaspout.h5` | VASP 编译时没有 `-DVASP_HDF5` | 重新编译 VASP 并链接 HDF5，或改用 vasprun.xml 系工具 |
| `Image export using the "kaleido" engine requires the kaleido package` | plotly 导出静态图需要 kaleido | `pip install -U kaleido`，或导出 CSV 后用 matplotlib 画 |
| `findfont: Font family ['Arial'] not found. Falling back to DejaVu Sans` | Linux 上没有 Arial | 装 Liberation Sans 等度量兼容字体后清除 matplotlib 字体缓存，或接受 DejaVu Sans |
| PDF 中文字在 Illustrator 中变成轮廓/乱码 | 默认 Type 3 字体 | `pdf.fonttype : 42` |
| `pip install ovito` 装上了旧版本或找不到版本 | 环境的 Python < 3.11 | 用 Python 3.11/3.12 的 venv（不要装进 conda 环境） |
| nglview 单元格输出为空白 | JupyterLab 前端扩展未加载、版本不匹配 | 在同一环境用 conda-forge 安装 nglview 与 jupyterlab，重启 JupyterLab 和浏览器 |
| VESTA 差分电荷满屏条纹 | AB、A、B 网格或坐标不同 | 显式统一 `NGXF/NGYF/NGZF`、ENCUT 与坐标，片段不弛豫 |
| Multiwfn 运行大体系时 `Segmentation fault` | 栈空间不足 | `ulimit -s unlimited`、`export OMP_STACKSIZE=200M`（或更大），settings.ini 中设 `nthreads` |
| Multiwfn 批处理一直不结束 | 没加 `-silent`，程序在等图形窗口 | 加 `-silent` 或在 settings.ini 中设 `isilent=1` |

---

## 练习

1. 运行 [05-band-dos-plotting](../examples/05-band-dos-plotting/README.md)，分别用 sumo、VASPKIT 211 和 `plot_band_dos_matplotlib.py` 画 GaAs 能带，比较三者的零点与带隙；把 VASPKIT 的零点改成 VBM，让三张图一致。
2. 用 `sumo-bandplot --project Ga.s,Ga.p,As.p` 与 PyProcar 的 `parametric` 模式分别画 fat band，判断 GaAs 价带顶和导带底主要由哪个元素、哪个轨道贡献。
3. 把 DOS 计算的 `ISMEAR` 改为 0（`SIGMA = 0.05`）、`NEDOS` 改为 301，与 `ISMEAR = -5`、`NEDOS = 3001` 的结果画在同一张图上，说明差别来自哪里。
4. 用 `sumo-bandstats` 计算 GaAs 导带底的电子有效质量，改变路径点数（20、40、80）和 `-s` 采样点数，观察结果的变化，并与文献实验值（约 0.067 $m_e$，Vurgaftman 等 2001）比较，讨论 PBE 带隙偏小对有效质量的影响。
5. 选一个吸附体系（如 [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md) 的 CO/Pt(111)），分别用 VASPKIT 314 和 VESTA 的 Edit Data 做差分电荷，确认两者一致，并导出一张带正负等值面、注明等值面数值的图。
6. 用 `ovito_render_xdatcar.py` 渲染一段 AIMD 轨迹的 GIF，并比较 `--every-nth 1` 与 `--every-nth 10` 计算的 RDF。
7. 修改 `dft_paper.mplstyle` 为双栏尺寸（7.0 × 3.0 in），把能带、PDOS 和光学吸收谱排成 1×3 子图，导出 PDF 并在 Inkscape/Illustrator 中确认文字可编辑。

---

## 参考资料

官方文档：

- VASPKIT 1.5 文档：功能列表 <https://vaspkit.com/features.html>，安装与 `~/.vaspkit` <https://vaspkit.com/installation.html>，教程 <https://vaspkit.com/tutorials.html>，更新日志 <https://vaspkit.com/changelog.html>
- sumo 文档与教程：<https://smtg-bham.github.io/sumo/>，<https://smtg-bham.github.io/sumo/tutorials.html>
- PyProcar 文档（v6）：<https://pyprocar.readthedocs.io/en/stable/>
- py4vasp 文档：<https://vasp.at/py4vasp/latest>；VASP Wiki：[vaspout.h5](https://vasp.at/wiki/Vaspout.h5)、[LORBIT](https://vasp.at/wiki/LORBIT)、[ISMEAR](https://vasp.at/wiki/ISMEAR)、[LOPTICS](https://vasp.at/wiki/LOPTICS)
- easyunfold 文档：<https://smtg-bham.github.io/easyunfold/>；IFermi 文档：<https://fermisurfaces.github.io/IFermi/>
- VESTA 下载（含 3.90.x 开发版）：<https://jp-minerals.org/vesta/en/download.html>
- OVITO 文档与更新说明：<https://docs.ovito.org/>，<https://docs.ovito.org/new_features.html>；Python 模块 `ovito.vis`：<https://docs.ovito.org/python/modules/ovito_vis.html>
- Multiwfn 官网与手册：<http://sobereva.com/multiwfn/>，<http://sobereva.com/multiwfn/Multiwfn_manual.html>
- VMD：<https://www.ks.uiuc.edu/Research/vmd/>；Avogadro：<https://two.avogadro.cc/>；Jmol：<https://jmol.sourceforge.net/>
- phonon website：<https://henriquemiranda.github.io/phononwebsite/>；phonopy 设置项（EIGENVECTORS、BAND_CONNECTION）：<https://phonopy.github.io/phonopy/setting-tags.html>
- matplotlib 样式表与 rcParams：<https://matplotlib.org/stable/users/explain/customizing.html>
- SeeK-path 在线工具：<https://www.materialscloud.org/work/tools/seekpath>

相关 GitHub 仓库（只链接，不复制代码）：

- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) — 命令行出版级能带/DOS/光学/声子图，sumo-kgen 生成 k 路径（含杂化泛函拆分），sumo-bandstats 算有效质量。
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar) — 投影能带、费米面、能带反折叠；v6 起 API 改为 `bandsplot(code='vasp', dirname=...)`。
- [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) — VASP 官方 Python 后处理接口，读取 vaspout.h5。
- [SMTG-Bham/easyunfold](https://github.com/SMTG-Bham/easyunfold) — 超胞能带反折叠，支持 VASP 与 CASTEP。
- [fermisurfaces/IFermi](https://github.com/fermisurfaces/IFermi) — 费米面生成与可视化（费米速度、自旋纹理）。
- [tamaswells/VASPKIT_manual](https://github.com/tamaswells/VASPKIT_manual) — VASPKIT 中文手册（v0.73）源码，编号以 1.5 文档为准；无许可证，仅供链接。
- [nglviewer/nglview](https://github.com/nglviewer/nglview) — Jupyter 中交互显示 ASE/pymatgen 结构与轨迹。
- [superstar54/weas-widget](https://github.com/superstar54/weas-widget) — Jupyter 中可编辑的结构查看器。
- [janosh/matterviz](https://github.com/janosh/matterviz) — 浏览器中交互展示晶体、轨迹、体数据、能带+DOS 与相图。
- [OpenChemistry/avogadroapp](https://github.com/OpenChemistry/avogadroapp) — 开源分子编辑器，可显示轨迹与振动。
- [digital-chemistry-laboratory/multiwfn-mac-build](https://github.com/digital-chemistry-laboratory/multiwfn-mac-build) — macOS 上编译 Multiwfn 的方案。
- [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn) — Multiwfn 批处理/工作流工具（中文开发者）。

关键文献：

- V. Wang, N. Xu, J.-C. Liu, G. Tang, W.-T. Geng, VASPKIT: A user-friendly interface facilitating high-throughput computing and analysis using VASP code, *Comput. Phys. Commun.* **267**, 108033 (2021). DOI: [10.1016/j.cpc.2021.108033](https://doi.org/10.1016/j.cpc.2021.108033)
- A. M. Ganose, A. J. Jackson, D. O. Scanlon, sumo: Command-line tools for plotting and analysis of periodic ab initio calculations, *J. Open Source Softw.* **3**, 717 (2018). DOI: [10.21105/joss.00717](https://doi.org/10.21105/joss.00717)
- U. Herath et al., PyProcar: A Python library for electronic structure pre/post-processing, *Comput. Phys. Commun.* **251**, 107080 (2020). DOI: [10.1016/j.cpc.2019.107080](https://doi.org/10.1016/j.cpc.2019.107080)
- B. Zhu, S. R. Kavanagh, D. O. Scanlon, easyunfold: A Python package for unfolding electronic band structures, *J. Open Source Softw.* **9**, 5974 (2024). DOI: [10.21105/joss.05974](https://doi.org/10.21105/joss.05974)
- A. M. Ganose, A. Searle, A. Jain, S. M. Griffin, IFermi: A python library for Fermi surface generation and analysis, *J. Open Source Softw.* **6**, 3089 (2021). DOI: [10.21105/joss.03089](https://doi.org/10.21105/joss.03089)
- K. Momma, F. Izumi, VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data, *J. Appl. Crystallogr.* **44**, 1272 (2011). DOI: [10.1107/S0021889811038970](https://doi.org/10.1107/S0021889811038970)
- A. Stukowski, Visualization and analysis of atomistic simulation data with OVITO – the Open Visualization Tool, *Modelling Simul. Mater. Sci. Eng.* **18**, 015012 (2010). DOI: [10.1088/0965-0393/18/1/015012](https://doi.org/10.1088/0965-0393/18/1/015012)
- T. Lu, F. Chen, Multiwfn: A multifunctional wavefunction analyzer, *J. Comput. Chem.* **33**, 580 (2012). DOI: [10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)
- W. Humphrey, A. Dalke, K. Schulten, VMD: Visual molecular dynamics, *J. Mol. Graph.* **14**, 33 (1996). DOI: [10.1016/0263-7855(96)00018-5](https://doi.org/10.1016/0263-7855(96)00018-5)
- H. Nguyen, D. A. Case, A. S. Rose, NGLview – interactive molecular graphics for Jupyter notebooks, *Bioinformatics* **34**, 1241 (2018). DOI: [10.1093/bioinformatics/btx789](https://doi.org/10.1093/bioinformatics/btx789)
- P. E. Blöchl, O. Jepsen, O. K. Andersen, Improved tetrahedron method for Brillouin-zone integrations, *Phys. Rev. B* **49**, 16223 (1994). DOI: [10.1103/PhysRevB.49.16223](https://doi.org/10.1103/PhysRevB.49.16223)
- W. Setyawan, S. Curtarolo, High-throughput electronic band structure calculations: Challenges and tools, *Comput. Mater. Sci.* **49**, 299 (2010). DOI: [10.1016/j.commatsci.2010.05.010](https://doi.org/10.1016/j.commatsci.2010.05.010)
- Y. Hinuma et al., Band structure diagram paths based on crystallography, *Comput. Mater. Sci.* **128**, 140 (2017). DOI: [10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015)
- I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, Band parameters for III–V compound semiconductors and their alloys, *J. Appl. Phys.* **89**, 5815 (2001). DOI: [10.1063/1.1368156](https://doi.org/10.1063/1.1368156)（GaAs 带隙与有效质量的参考值）
- B. Wong, Points of view: Color blindness, *Nat. Methods* **8**, 441 (2011). DOI: [10.1038/nmeth.1618](https://doi.org/10.1038/nmeth.1618)（Okabe–Ito 色盲友好配色）
- J. D. Hunter, Matplotlib: A 2D graphics environment, *Comput. Sci. Eng.* **9**, 90 (2007). DOI: [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)

[← 上一章](04-ase-pymatgen.md) | [返回目录](../README.md) | [下一章 →](06-workflow-automation.md)
