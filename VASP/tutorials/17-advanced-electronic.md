# 第17章 高级电子结构：Wannier90、有效质量、能带反折叠、费米面与拓扑（Advanced Electronic Structure with VASP）

本章介绍在普通能带和态密度之上的进阶电子结构分析：通过 VASP2WANNIER90 接口构造最大局域 Wannier 函数（MLWF）并做能带插值和紧束缚模型，用有限差分/抛物线拟合计算载流子有效质量（effective mass），用 easyunfold、vaspkit 做超胞能带反折叠（band unfolding），用 IFermi、vaspkit + FermiSurfer 画费米面（Fermi surface），以及 Berry 曲率、Z2、不可约表示、Rashba 自旋劈裂等拓扑与自旋相关分析。最后简要介绍点缺陷形成能与电荷修正（doped、pydefect）。配套算例为 Si、GaAs、Si 空位超胞和 fcc Cu。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 编译带 Wannier90 接口的 VASP，用 `LWANNIER90`、`NUM_WANN`、`WANNIER90_WIN` 等标签生成 Wannier90 输入，设置投影和解纠缠窗口，检查 spread 收敛，并把插值能带与 DFT 能带对比；
- 知道 VASP 内部的 Wannier 化方法（`LOCPROJ`、`LSCDM`）和 `KPOINTS_WAN`（VASP ≥ 6.3）及其已知问题；
- 计算有效质量张量的主值，做正确的单位换算，判断结果是否受非抛物线性、能带简并和带隙误差影响；
- 完成一次超胞能带反折叠，理解谱权重的含义；
- 在稠密网格上画费米面并提取费米速度；
- 了解 Berry 曲率、陈数、Z2、不可约表示、Rashba 劈裂的计算路线和常用工具；
- 了解带电点缺陷形成能的公式、化学势和有限尺寸修正的处理思路。

| 任务 | VASP 关键设置 | 主要工具 | 本章小节 |
|---|---|---|---|
| Wannier 插值能带、紧束缚哈密顿量 | `-DVASP2WANNIER90`、`LWANNIER90`、`ISYM = -1`、`NBANDS`、`NUM_WANN` | Wannier90、postw90 | §3 |
| HSE/GW 能带的廉价插值 | 杂化泛函自洽 + `LWANNIER90` | Wannier90 | §3.7 |
| 载流子有效质量 | 带边附近的 k 点（非自洽或 KPOINTS_OPT） | 本章脚本、emc、effmass、vaspkit 912、sumo-bandstats | §4 |
| 超胞（缺陷、合金、表面）能带 | 非自洽 + `LWAVE = .TRUE.` | easyunfold、vaspkit 28x、VaspBandUnfolding | §5 |
| 金属费米面、费米速度 | 稠密均匀网格非自洽 | IFermi、vaspkit 26x、FermiSurfer、XCrySDen | §6 |
| Berry 曲率、陈数、Z2、表面态 | SOC（`vasp_ncl`）、Wannier90 | VASPBERRY、postw90、WannierTools、Z2Pack | §7 |
| 能带不可约表示、拓扑量子化学 | 标准化晶胞、WAVECAR | IRVSP | §7.4 |
| Rashba 劈裂、自旋纹理 | `LSORBIT = .TRUE.` | vaspkit 65x、pyprocar、IFermi | §7.1 |
| 点缺陷形成能、转变能级 | 带电超胞（`NELECT`）、LOCPOT | doped、pydefect | §8 |

> [!NOTE]
> 本章建立在 [第5章](05-band-dos.md)（能带与 DOS、KPOINTS_OPT）、[第6章](06-hybrid-metagga.md)（HSE06）和 [第7章](07-magnetism-u-soc.md)（SOC）的基础上。所有分析都应在已经收敛的结构和参数上进行。

---

## 2. 原理速览

### 2.1 Wannier 函数

Wannier 函数是 Bloch 态的傅里叶变换，并允许在每个 k 点对一组能带做幺正混合：

$$|w_{m\mathbf R}\rangle=\frac{V}{(2\pi)^3}\int_\mathrm{BZ}\mathrm d\mathbf k\,e^{-i\mathbf k\cdot\mathbf R}\sum_{n}U_{nm}^{(\mathbf k)}|\psi_{n\mathbf k}\rangle$$

Marzari–Vanderbilt 方法选择 $U^{(\mathbf k)}$ 使总扩展（spread）最小：

$$\Omega=\sum_m\left[\langle w_{m\mathbf 0}|r^2|w_{m\mathbf 0}\rangle-|\langle w_{m\mathbf 0}|\mathbf r|w_{m\mathbf 0}\rangle|^2\right]$$

当目标能带与其他能带纠缠在一起时（例如 Si 的导带），需要先**解纠缠**（Souza–Marzari–Vanderbilt）：在外窗口（`dis_win_min`–`dis_win_max`）内的能带中挑选出最光滑的 `num_wann` 维子空间，冻结窗口（`dis_froz_min`–`dis_froz_max`）内的态被完整保留，因而插值能带在冻结窗口内与 DFT 精确一致。

得到的 Wannier 函数构成一个最小基组，哈密顿量在实空间的矩阵元 $H_{mn}(\mathbf R)=\langle w_{m\mathbf 0}|H|w_{n\mathbf R}\rangle$（`wannier90_hr.dat`）可以在任意 k 点廉价地对角化——这就是 Wannier 插值，也是 WannierTools 等拓扑计算的输入。

### 2.2 有效质量

在带边 $\mathbf k_0$ 附近展开：

$$E(\mathbf k)\approx E(\mathbf k_0)+\frac{\hbar^2}{2}\sum_{ij}\Delta k_i\,(M^{-1})_{ij}\,\Delta k_j,\qquad (M^{-1})_{ij}=\frac{1}{\hbar^2}\frac{\partial^2E}{\partial k_i\partial k_j}$$

对角化逆有效质量张量得到主值与主轴。沿某个方向 $\hat{\mathbf u}$：

$$\frac{m^*}{m_0}=\frac{\hbar^2}{m_0}\left(\frac{\mathrm d^2E}{\mathrm dk^2}\right)^{-1}=\frac{7.61996\ \mathrm{eV\,\AA^2}}{\mathrm d^2E/\mathrm dk^2\ [\mathrm{eV\,\AA^2}]}$$

若拟合 $E=E_0+a\,\Delta k^2$，则 $m^*/m_0 = 3.80998/a$（$a$ 的单位为 eV·Å²，k 的单位为 Å⁻¹ 且**包含 2π**）。原子单位下（k 用 bohr⁻¹、E 用 Hartree）则直接有 $m^*=1/(\mathrm d^2E/\mathrm dk^2)$。

三点有限差分：$\dfrac{\mathrm d^2E}{\mathrm dk^2}\approx\dfrac{E(k_0+h)-2E(k_0)+E(k_0-h)}{h^2}$。

### 2.3 能带反折叠

超胞的布里渊区更小，原胞能带被"折叠"进去，变得难以辨认。对超胞本征态 $|\mathbf K m\rangle$，把它投影到原胞的 k 点 $\mathbf k_i$（满足 $\mathbf k_i=\mathbf K+\mathbf G$）上得到谱权重：

$$P_{\mathbf Km}(\mathbf k_i)=\sum_n|\langle\mathbf K m|\mathbf k_i n\rangle|^2,\qquad A(\mathbf k_i,E)=\sum_m P_{\mathbf Km}(\mathbf k_i)\,\delta(E-E_m)$$

在平面波基下，$P$ 只需对满足 $\mathbf G=\mathbf g+\mathbf k_i-\mathbf K$ 的平面波系数平方求和（Popescu–Zunger）。完美超胞的 $P$ 为 0 或 1；有缺陷、无序时 $P$ 变为分数，能带"模糊"。

### 2.4 费米面

金属的费米面是 $E_n(\mathbf k)=E_F$ 的等能面。费米速度 $\mathbf v_n=\frac{1}{\hbar}\nabla_\mathbf kE_n$。由于等能面需要精细的 k 网格，通常在中等密度网格上算本征值，再用傅里叶插值（BoltzTraP2）或线性插值加密。

### 2.5 Berry 曲率、陈数与 Z2

对第 n 条带，Berry 联络与曲率为

$$\mathbf A_n(\mathbf k)=i\langle u_{n\mathbf k}|\nabla_\mathbf k|u_{n\mathbf k}\rangle,\qquad \boldsymbol\Omega_n(\mathbf k)=\nabla_\mathbf k\times\mathbf A_n(\mathbf k)$$

二维体系占据带的陈数 $C=\frac{1}{2\pi}\sum_{n\in\mathrm{occ}}\int_\mathrm{BZ}\Omega_n^z\,\mathrm d^2k$ 为整数。时间反演不变体系的 Z2 拓扑不变量可以通过 Wannier 电荷中心（Wilson loop）随 k 的演化来判断，有反演对称性时也可以用时间反演不变动量处的宇称乘积（Fu–Kane 判据）。

### 2.6 Rashba 自旋劈裂

缺少反演对称的表面、界面或极性晶体中，SOC 产生动量依赖的有效磁场。最简单的二维 Rashba 模型：

$$H_R=\alpha_R(\boldsymbol\sigma\times\mathbf k)\cdot\hat z,\qquad E_\pm(k)=\frac{\hbar^2k^2}{2m^*}\pm\alpha_R|k|,\qquad \alpha_R=\frac{2E_R}{k_R}$$

其中 $E_R$ 是能带交叉点到能带极值的能量差，$k_R$ 是极值相对交叉点的动量偏移，可以直接从 SOC 能带上读出。

---

## 3. Wannier90 接口

配套算例：[examples/17-wannier90-si](../examples/17-wannier90-si/README.md)。

### 3.1 编译与版本

VASP 侧（详细编译流程见 [第1章](01-setup-compile.md)）：

```make
# makefile.include 中需要增加的内容（只列改动，不复制官方模板）
CPP_OPTIONS    += -DVASP2WANNIER90
WANNIER90_ROOT ?= /path/to/wannier90-3.1.0      # directory that contains lib/libwannier.a
LLIBS          += -L$(WANNIER90_ROOT)/lib -lwannier
```

Wannier90 侧（3.x）：

```bash
cd wannier90-3.1.0
cp config/make.inc.ifort make.inc        # 或 make.inc.gfort，按编译器选择
# 按 VASP Wiki 的要求，链接到 VASP 的库应编译为串行版本：注释掉 make.inc 中的 COMMS = mpi
make lib                                 # 生成 libwannier.a
make wannier                             # 同时编译独立程序 wannier90.x（以及 make post 得到 postw90.x）
```

版本说明：

| 组合 | 说明 |
|---|---|
| VASP ≤ 6.1 + Wannier90 1.x / 2.x | `-DVASP2WANNIER90`（1.x）或 `-DVASP2WANNIER90v2`（2.x），已过时 |
| VASP ≥ 6.2 + Wannier90 2.x / 3.x | 统一用 `-DVASP2WANNIER90`；VASP Wiki 目前只写明支持 2.x/3.x |
| Wannier90 4.0.x（2026 年发布，4.0.3 起库名统一为 `libwannier90`、提供 MPI 并行库接口） | VASP 官方文档尚未说明与 4.x 库的兼容性；建议链接时仍用 3.1.0。独立运行的 `wannier90.x` 最好也用同一版本，避免关键词差异（4.0.3 删除了一些过时关键词） |

> [!TIP]
> 判断 VASP 是否带有接口：设 `LWANNIER90 = .TRUE.` 跑一个小算例，若运行后没有 `wannier90.amn`、`wannier90.mmn`、`wannier90.eig` 等文件，说明可执行程序没有编译接口。

### 3.2 INCAR 标签

| 标签 | 默认 | 说明 |
|---|---|---|
| `LWANNIER90` | `.FALSE.` | 打开接口：调用 `wannier_setup`，写出 `wannier90.win`（不存在时写模板，存在时只补充缺少的晶胞、原子、`mp_grid`、k 点块）、`.amn`、`.mmn`、`.eig` |
| `LWANNIER90_RUN` | `.FALSE.` | 在 VASP 内以库模式继续调用 `wannier_run`，一次完成 Wannier 化；此时默认不写 `.amn/.mmn` |
| `LWRITE_MMN_AMN` | `LWANNIER90` 时为 `.TRUE.` | 强制写出 `.mmn/.amn`（库模式下想保留文件以便日后独立运行时设置） |
| `LWRITE_UNK` | `.FALSE.` | 写出 `UNKp.s` 周期部分波函数，用于 Wannier90 画 Wannier 函数（`wannier_plot`）；非共线时写出的 UNK 文件在 6.5.1 之前不能被 Wannier90 读取 |
| `NUM_WANN` | `NBANDS` | Wannier 函数个数（VASP ≥ 6.2 也用于 Wannier90，写入 win 文件的 `num_wann` 总是这个值） |
| `WANNIER90_WIN` | 无 | 多行字符串，直接在 INCAR 中给出 win 文件内容（VASP ≥ 6.2）；总长度上限为 8192 字符（Known issue 64） |
| `ISYM` | — | 设 `-1`：Wannier90 需要完整的均匀 Monkhorst–Pack 网格，不能对称约化 |
| `NBANDS` | — | 外窗口内要有足够的能带；必须 ≥ `NUM_WANN` |

`WANNIER90_WIN` 写法示例（与示例目录中独立 win 文件等价，VASP 会自动补充 `num_wann`、`num_bands`、晶胞、原子、`mp_grid` 和 k 点）：

```fortran
NUM_WANN = 8
WANNIER90_WIN = "
num_iter      = 200
dis_num_iter  = 500
dis_win_max   = 17.0
dis_froz_max  = 8.0
begin projections
Si:sp3
end projections
bands_plot    = true
begin kpoint_path
L 0.5 0.5 0.5  G 0.0 0.0 0.0
G 0.0 0.0 0.0  X 0.5 0.0 0.5
end kpoint_path
"
```

### 3.3 Si 的完整流程

**第一步：自洽**（可以使用对称性）：

```fortran
SYSTEM  = Si SCF before Wannierisation
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-8
ISMEAR  = 0
SIGMA   = 0.05
LCHARG  = .TRUE.
```

**第二步：非自洽 + 接口**（完整 6×6×6 Γ 中心网格）：

```fortran
SYSTEM     = Si non-SCF + VASP2WANNIER90 interface
PREC       = Accurate
ENCUT      = 400
ISMEAR     = 0
SIGMA      = 0.05
ICHARG     = 11
ISYM       = -1
NBANDS     = 16
NUM_WANN   = 8
LWANNIER90 = .TRUE.
LCHARG     = .FALSE.
LWAVE      = .FALSE.
```

**第三步：wannier90.win**（关键部分）：

```text
num_wann        = 8
num_iter        = 200
dis_win_max     = 17.0      ! absolute energy: ~E_F + 11 eV
dis_froz_max    = 8.0       ! absolute energy: ~E_F + 2 eV
dis_num_iter    = 500
begin projections
Si:sp3
end projections
guiding_centres = true
write_hr        = true
bands_plot      = true
begin kpoint_path
L 0.500 0.500 0.500  G 0.000 0.000 0.000
G 0.000 0.000 0.000  X 0.500 0.000 0.500
X 0.500 0.000 0.500  U 0.625 0.250 0.625
K 0.375 0.375 0.750  G 0.000 0.000 0.000
end kpoint_path
```

**第四步：运行 Wannier90**：

```bash
cd wannier
mpirun -np 16 vasp_std > vasp.out           # 写出 wannier90.amn/.mmn/.eig，并补全 wannier90.win
ls wannier90.*
wannier90.x wannier90                        # 独立运行（参数为 seedname）
```

### 3.4 wannier90.win 关键词

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `num_wann` | 8 | 两个 Si 原子各 4 个 sp3 杂化轨道；与 INCAR 的 `NUM_WANN` 一致 |
| `num_bands` | 由 VASP 写入（= NBANDS） | 参与 Wannier 化的能带数；也可用 `exclude_bands` 排除半芯等能带 |
| `begin projections` | `Si:sp3` | 初始猜测。常见写法：`Ga:d`、`As:p`、`f=0.25,0.25,0.25:s`（指定位置的 s 轨道）、`C:pz` |
| `dis_win_min/max` | 默认最低 / 17.0 | 外窗口（绝对能量，eV） |
| `dis_froz_min/max` | 默认 / 8.0 | 冻结窗口：窗口内每个 k 点的态数不能超过 `num_wann` |
| `dis_num_iter` | 500 | 解纠缠迭代步数 |
| `num_iter` | 200 | 局域化迭代步数 |
| `guiding_centres` | true | 迭代初期用投影中心引导，避免 Wannier 中心漂移 |
| `bands_plot`、`kpoint_path`、`bands_num_points` | true / 路径 / 100 | 输出 `wannier90_band.dat`、`wannier90_band.gnu` |
| `write_hr` | true | 输出实空间哈密顿量 `wannier90_hr.dat` |
| `write_xyz` | true | 输出 Wannier 中心 `wannier90_centres.xyz` |
| `wannier_plot` | （画图时）true | 输出 Wannier 函数的实空间格点文件，需要 `LWRITE_UNK = .TRUE.` |
| `spinors` | 由 VASP 写入 | 非共线计算时为 true |

### 3.5 怎样选窗口

1. 先算能带和投影 DOS（第5章），确定目标轨道（这里是 Si 3s3p）主要分布的能量范围；
2. `grep E-fermi OUTCAR` 得到费米能（绝对能量），窗口必须用同一零点；
3. 冻结窗口上限选在"目标子空间仍然完整、且任一 k 点上窗口内态数 ≤ num_wann"的位置。Si 的 8 个 sp3 函数包含 4 条价带 + 4 条导带，冻结到导带底以上约 1–2 eV 一般可行；
4. 外窗口要足够宽，使每个 k 点窗口内至少有 `num_wann` 个态；
5. 窗口不合适时 Wannier90 会直接报错（见常见问题），按提示调整即可。

示例的作业脚本会读取自洽的 E-fermi，自动把窗口设为 E_F+2 eV 和 E_F+11 eV。

### 3.6 检查 Wannier 化质量

```bash
grep "<-- DIS"  wannier90.wout | tail -3     # 解纠缠收敛：Omega_I 不再变化
grep "<-- CONV" wannier90.wout | tail -3     # 局域化收敛：总 spread 不再下降
grep -A 12 "Final State" wannier90.wout      # 每个 WF 的中心坐标与 spread
gnuplot -persist wannier90_band.gnu          # 插值能带
```

合格的 Wannier 化应满足：

- spread 单调下降并收敛，各个对称等价 WF 的 spread 接近；
- Wannier 中心在化学上合理的位置（Si 的 sp3 中心在键上）；
- 冻结窗口内插值能带与 DFT 能带重合到 meV 级，没有"毛刺"；
- 若插值能带在冻结窗口内出现尖峰振荡，通常是 k 网格太稀（实空间截断不够）或投影选择不当。

### 3.7 用 Wannier 插值得到 HSE/GW 能带

杂化泛函不能做 `ICHARG = 11` 非自洽能带（第6章），零权重 k 点法又很贵。Wannier 插值提供了第三种途径：

```fortran
# HSE06 self-consistent run on a uniform mesh, Wannier90 files written at the end
LHFCALC    = .TRUE.
HFSCREEN   = 0.2
ALGO       = Damped
TIME       = 0.4
ISYM       = -1          ! full mesh for Wannier90 (more expensive than with symmetry)
NBANDS     = 16
NUM_WANN   = 8
LWANNIER90 = .TRUE.
```

自洽结束后运行 `wannier90.x`，得到的插值能带即为 HSE06 能带。GW 计算也可以用同样的接口把准粒子能量插值到路径上（参考 [JMSkelton/VASP-Examples](https://github.com/JMSkelton/VASP-Examples) 中 Si 的 HSE06 Wannier 插值算例）。注意：`ISYM = -1` 的杂化泛函计算 k 点更多，代价可能是对称约化的数倍。

### 3.8 VASP 内部的 Wannier 化与 KPOINTS_WAN

VASP 6 起可以不调用 Wannier90 而在内部构造（非最大局域的）Wannier 函数：

| 方法 | 标签 | 特点 |
|---|---|---|
| 局域投影 + 一次 SVD | `LOCPROJ` | 用户指定投影目标（如 `LOCPROJ = 1 : sp3 : Hy 1 4.`），控制力强 |
| SCDM | `LSCDM = .TRUE.`、`NUM_WANN`、`CUTOFF_TYPE/CUTOFF_MU/CUTOFF_SIGMA` | 从密度矩阵自动挑选局域函数，无需投影 |
| 交给 Wannier90 做最大局域化 | 以上两者 + `LWANNIER90` | 用 SCDM 时建议 `dis_num_iter = 0` |

VASP ≥ 6.3.0 读取 `KPOINTS_WAN`（格式同 KPOINTS）时，会用 Wannier 函数插值得到这些 k 点上的本征值，写入 vasprun.xml 和 vaspout.h5 的 `kpoints_wan` 段，可用 py4vasp 的 `calc.band.plot("kpoints_wan")` 画出，并与 `KPOINTS_OPT` 的直接计算结果对比。

> [!WARNING]
> VASP Known issues 中有两条与此相关、截至 Wiki 最新状态仍未解决的问题：ID 105——对高对称晶格，`KPOINTS_WAN` 插值能带是错的（等距晶格矢量只保留了最后一个），KPOINTS 网格越稀误差越明显；ID 106——`LSCDM` 不经过 Wannier90 库时能量窗口被应用了两次。因此对 Si、Cu 这类立方晶体，建议用 Wannier90 自己的 `bands_plot`，或至少把 `KPOINTS_WAN` 的结果与 `KPOINTS_OPT` 对比后再使用。

### 3.9 SOC 与磁性体系

- 用 `vasp_ncl` + `LSORBIT = .TRUE.`，VASP 会在 win 中写 `spinors = true`，投影数（`num_wann`）加倍；
- 必须 `ISYM = -1`：6.4.3 之前，非共线 + 对称性时写出的 `.amn/.mmn` 缺少旋量部分的转动（Known issue 19）；
- VASP 6.2.0–6.3.0 中非共线时 win 文件里的自旋量子化轴读取有误（Known issue 79，6.3.1 修复）；使用 VASP 5.4.4 + Wannier90 2.x 的用户可参考 [Chengcheng-Xiao/VASP2WAN90_v2_fix](https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix)；
- 共线自旋极化（`ISPIN = 2`）时，两个自旋通道要分别做 Wannier 化：VASP 会为每个通道生成一套 Wannier90 输入（seedname 中带自旋标记，以实际生成的文件名为准），再对每一套分别运行 Wannier90。

---

## 4. 有效质量

配套算例：[examples/17-effmass-gaas](../examples/17-effmass-gaas/README.md)。

### 4.1 流程

1. 先算能带（第5章），确定 VBM、CBM 所在的 k 点（GaAs 都在 Γ）；
2. 自洽得到 CHGCAR；
3. 在带边 k0 附近沿若干方向取小步长 k 点，做非自洽（`ICHARG = 11`），或者在自洽目录中放 `KPOINTS_OPT`（VASP ≥ 6.3）；
4. 对每个方向拟合 $E-E_0=a\,\Delta k^2$，换算 $m^*/m_0=3.80998/a$；
5. 需要完整张量时，补充非对角方向（如 [110] 与 [1-10]），由二阶导数构造 3×3 矩阵后对角化。

### 4.2 k 点文件

`KPOINTS_fd` 使用 Cartesian 坐标，单位为 2π/a（a 为 POSCAR 第二行的缩放因子），这样方向和步长一目了然：

```text
k-points around Gamma for finite-difference effective masses (Cartesian, units of 2pi/a)
16
Cartesian
    0.00000000   0.00000000   0.00000000   1   ! Gamma
    0.00500000   0.00000000   0.00000000   1   ! [100]  |k|=0.005
    ...（[100]、[110]、[111] 各 5 个点，|k| = 0.005 ... 0.025）
```

对 a = 5.653 Å，步长 0.005·2π/a ≈ 0.0056 Å⁻¹，最大 |Δk| ≈ 0.028 Å⁻¹。步长选择的原则：太小时能量差只有零点几 meV，被数值噪声淹没；太大时偏离抛物线。用 `--npoints` 比较只用内侧 3 个点与全部 5 个点的结果，差别小说明抛物线近似成立。

### 4.3 INCAR

```fortran
SYSTEM  = GaAs non-SCF on k-points around Gamma
PREC    = Accurate
ENCUT   = 400
EDIFF   = 1E-8         ! meV-scale energy differences
ISMEAR  = 0
SIGMA   = 0.01
ICHARG  = 11
NBANDS  = 16
```

### 4.4 拟合

```bash
python3 fit_effmass.py --eigenval EIGENVAL --poscar POSCAR          # 默认 CBM、VBM
python3 fit_effmass.py --band 7 --band 8 --npoints 3                 # 价带顶简并的其他能带
```

输出格式（示意输出，数值仅供参考）：

```text
band 10 (cbm), E(k0) = x.xxxx eV, spin 1
  direction (cartesian)        npts    m*_fit/m0    m*_2pt/m0    max|dk|
  ( 0.5770  0.5770  0.5770)         5       0.0xxx       0.0xxx     0.0278
  ( 0.7070  0.7070  0.0000)         5       0.0xxx       0.0xxx     0.0278
  ( 1.0000  0.0000  0.0000)         5       0.0xxx       0.0xxx     0.0278
```

物理上需要注意：

- **带隙误差**：k·p 理论给出 Γ 导带质量大致与带隙成正比，PBE 严重低估 GaAs 带隙（甚至接近 0），电子有效质量也会被明显低估。实验值约 0.067 m0（Vurgaftman 等, 2001）。定量计算请用 HSE06（可加 SOC）或 mBJ；
- **价带简并**：不含 SOC 时 GaAs 价带顶在 Γ 三重简并；含 SOC 后劈裂出自旋–轨道分裂带（split-off），重空穴、轻空穴的质量都会改变。III-V 半导体的空穴质量应在 SOC 下计算；
- **翘曲**：重空穴带在 [100]、[110]、[111] 方向质量不同，报告时要写明方向，或给出方向平均；
- **非抛物线性**：窄带隙材料的导带非抛物线性强，"有效质量"依赖于能量范围。effmass 工具区分了曲率质量、输运质量和光学质量等定义。

### 4.5 杂化泛函下的有效质量

HSE 不能用 `ICHARG = 11`，两种做法：

1. **零权重 k 点**：把 `IBZKPT`（自洽的不可约 k 点及权重）复制为 KPOINTS，再在末尾追加 `KPOINTS_fd` 中的点并把权重设为 0，k 点总数相应修改，然后做杂化泛函自洽；
2. **KPOINTS_OPT**（VASP ≥ 6.3）：KPOINTS 用均匀网格，把 `KPOINTS_fd` 改名为 `KPOINTS_OPT`；对未屏蔽的杂化泛函（如 PBE0）需要 `HFRCUT = -1`（见第6章）。结果写入 `EIGENVAL_OPT`/vaspout.h5。

### 4.6 现成工具

| 工具 | 方法 | 备注 |
|---|---|---|
| [afonari/emc](https://github.com/afonari/emc) | 在指定 k 点周围生成差分 KPOINTS，读 EIGENVAL 后对角化有效质量张量 | 经典工具，给出主值与主轴 |
| [lucydot/effmass](https://github.com/lucydot/effmass) | 读 VASP 输出，自动找极值，支持抛物线、非抛物线（Kane）、光学质量等多种定义 | Python 包，`pip install effmass` |
| vaspkit 912 / 913 / 914 | 912 读 `VPKIT.in`（预处理生成 KPOINTS、后处理拟合），914 角分辨有效质量 | 官方说明目前只支持非磁、不带电的半导体 |
| `sumo-bandstats` | 从能带计算的 vasprun.xml 对带边做抛物线拟合 | sumo 3.0.0 |
| `easyunfold unfold effective-mass` | 从反折叠能带中拟合有效质量 | 适合合金与缺陷体系 |
| 本仓库 [plot_band.py](../scripts/plot_band.py) | 报告 VBM/CBM 位置与带隙，帮助确定 k0 | `python ../scripts/plot_band.py -h` |

vaspkit 912 的 `VPKIT.in` 结构（按 vaspkit 官方教程）：第一行 1 = 预处理（生成 KPOINTS）/ 2 = 后处理；第二行拟合点数；第三行 k 截断半径（Å⁻¹）；第四行任务数；之后每行给出起点、终点两个 k 点的分数坐标和方向标签。

---

## 5. 能带反折叠

配套算例：[examples/17-unfold-si-vacancy](../examples/17-unfold-si-vacancy/README.md)。

### 5.1 easyunfold 工作流

```bash
pip install easyunfold

# 1) 原胞的能带路径（建议用 sumo-kgen 生成，并确认原胞是"标准原胞"）
sumo-kgen -p POSCAR_prim            # 写出 KPOINTS_band；若提示非标准原胞，则使用 POSCAR_prim 的标准化版本

# 2) 超胞需要计算的 k 点：superccell = M @ primitive
easyunfold generate POSCAR_prim POSCAR_sc KPOINTS_band --matrix "3 3 3"
#   -> easyunfold.json（反折叠信息）与 KPOINTS_easyunfold

# 3) 超胞非自洽：KPOINTS_easyunfold 作为 KPOINTS，ICHARG = 11，LWAVE = .TRUE.
# 4) 计算谱权重并作图
easyunfold unfold calculate WAVECAR
easyunfold unfold plot --eref <E-fermi> --emin -4 --emax 4
```

超胞非自洽 INCAR：

```fortran
SYSTEM  = Si vacancy supercell - non-SCF for unfolding
PREC    = Accurate
ENCUT   = 320
ISMEAR  = 0
SIGMA   = 0.05
ICHARG  = 11
NBANDS  = 144          ! occupied (~106) + conduction bands
LORBIT  = 11           ! PROCAR for atom-projected unfolding (optional)
LWAVE   = .TRUE.       ! easyunfold reads WAVECAR
```

### 5.2 要点

- `--matrix` 是超胞相对原胞的变换矩阵（行向量约定 cell_super = M @ cell_primitive），对角时写 3 个数，否则按行写 9 个数；超胞没有做晶胞弛豫时可以省略，easyunfold 会自动识别；
- 超胞做过**晶胞**弛豫时，反折叠对应的是"同样形变后的原胞"布里渊区，与理想原胞略有差别；
- easyunfold 默认利用原胞对称性展开并平均等价路径（`--no-symm-average` 关闭），这是缺陷超胞对称性降低时得到正确 EBS 的关键；
- `vasp_gam` 或 `vasp_ncl` 生成的 WAVECAR 需要在 `unfold calculate` 时加 `--gamma` 或 `--ncl`；
- 杂化泛函：`easyunfold generate ... --scf-kpoints IBZKPT --nk-per-split 60` 生成多个带零权重 k 点的 KPOINTS 分批计算，再把所有 WAVECAR 一起传给 `unfold calculate`；
- 原子/轨道投影：`easyunfold unfold plot-projections`（需要 `LORBIT = 11` 的 PROCAR），可以区分缺陷近邻原子与体相原子对 EBS 的贡献。

### 5.3 其他工具

| 工具 | 说明 |
|---|---|
| vaspkit 281–285 | 281 生成 KPOINTS（需要原胞的 `KPATH.in` 和描述超胞矩阵的 `TRANSMAT.in`），282 计算 EBS，283–285 为轨道投影的 EBS |
| [QijingZheng/VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding) | 直接解析 WAVECAR 的 Python 脚本集，除反折叠外还能算跃迁偶极矩、生成 HSE 能带 k 点 |
| [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar) | 也提供基于 PROCAR/WAVECAR 的反折叠功能 |

---

## 6. 费米面

配套算例：[examples/17-fermi-surface-cu](../examples/17-fermi-surface-cu/README.md)。

### 6.1 计算

```fortran
# non-SCF on a dense uniform mesh (Gamma-centred 25x25x25 for fcc Cu)
ICHARG  = 11
ISMEAR  = 1
SIGMA   = 0.1
NBANDS  = 12
LORBIT  = 11
```

- 金属的费米面对 k 点非常敏感：自洽可以用中等网格（15×15×15），画费米面的非自洽用更密的网格；
- 用原胞计算（超胞会把费米面折叠，难以解读）；
- 半金属、半导体掺杂研究时，可以把等能面能量 `--mu` 移到导带或价带中。

### 6.2 IFermi

```bash
ifermi info                                   # 读当前目录的 vasprun.xml：费米面个数、面积、维度、取向
ifermi info --property velocity               # 平均费米速度
ifermi plot                                   # plotly 交互式三维图
ifermi plot --interpolation-factor 10         # 插值倍数（默认 8）
ifermi plot --property velocity --property-colormap viridis
ifermi plot --slice 1 1 0 0                   # (110) 截面，过原点
ifermi plot --type matplotlib --output fs.png # 静态图
```

- IFermi 通过 pymatgen 读取 vasprun.xml，用 BoltzTraP2 做傅里叶插值，所以对称约化的网格也可以使用；
- `--property spin` 投影自旋纹理需要非共线计算，并且网格必须覆盖整个布里渊区（`ISYM = -1`）；
- 默认把费米面放在 Wigner–Seitz 胞中，`--no-wigner` 改用倒格子平行六面体。

### 6.3 vaspkit + FermiSurfer / XCrySDen

```bash
vaspkit     # 261：按分辨率（如 0.01×2π/Å）生成费米面用的 KPOINTS
            # 262：生成 FERMISURFACE.bxsf（XCrySDen）
            # 263/264：FermiSurfer 格式（264 带轨道投影，需要 PROCAR）
            # 266/267：二维材料的费米面
fermisurfer FERMISURFACE.frmsf   # 输出文件名以 vaspkit 实际生成的为准
```

pyprocar 也提供三维费米面绘制，适合需要轨道投影着色的场合。

---

## 7. 拓扑与自旋相关性质

### 7.1 SOC 能带、Rashba 劈裂与自旋纹理

```fortran
# vasp_ncl, read the collinear CHGCAR (chapter 7)
LSORBIT  = .TRUE.
ICHARG   = 11
ISYM     = -1
LORBIT   = 11
GGA_COMPAT = .FALSE.
```

- 在 Γ 附近沿两个相反方向（如 −M → Γ → M）取足够密的 k 点，SOC 劈裂的两条能带会形成"双抛物线"。从图中读出交叉点与极值的能量差 $E_R$ 和动量差 $k_R$，计算 $\alpha_R=2E_R/k_R$；
- 自旋纹理：PROCAR 在非共线计算中给出每个态的 $\langle\sigma_x\rangle,\langle\sigma_y\rangle,\langle\sigma_z\rangle$；vaspkit 651（生成二维 k 网格）+ 652/653（二维/体相自旋纹理，仅 SOC）、pyprocar 和 IFermi（`--property spin`）都可以画；
- Rashba 效应的经典体系是极性半导体 BiTeI 和重金属表面，参见 Ishizaka 等的工作。

### 7.2 Berry 曲率与陈数

| 路线 | 工具 | 输入 |
|---|---|---|
| 直接用平面波波函数 | [Infant83/VASPBERRY](https://github.com/Infant83/VASPBERRY) | WAVECAR（通常需要关闭对称性的均匀网格），用 Fukui–Hatsugai–Suzuki 格点方法求 Berry 曲率与陈数，也可计算圆二色性；具体输入以 README 为准 |
| Wannier 插值 | Wannier90 的 `postw90.x`（`berry = true`、`berry_task = ahc` 等） | `.chk` 与 `.win`，可计算反常霍尔电导、Berry 曲率分布 |
| Wannier 紧束缚 | WannierTools | `wannier90_hr.dat` |

Wannier 插值方法对 k 点收敛非常有效：一次 DFT + Wannier 化后，可以在 100×100×100 的插值网格上积分 Berry 曲率。

### 7.3 Z2 与表面态

- [WannierTools](https://github.com/quanshengwu/wannier_tools) 读入 `wannier90_hr.dat`，通过 `wt.in` 中的开关计算体能带、Wannier 电荷中心演化（Z2）、表面格林函数表面态、费米弧、Weyl 点搜索、Berry 曲率等。构造 Wannier 函数时要包含 SOC，并保证插值能带在费米能附近与 DFT 精确一致；
- [Z2Pack](https://github.com/Z2PackDev/Z2Pack) 通过混合 Wannier 电荷中心计算 Z2 与陈数，可以直接驱动 VASP 或使用紧束缚模型；
- 有反演对称性的体系可以用 Fu–Kane 宇称判据：在 8 个时间反演不变动量处统计占据 Kramers 对的宇称乘积。

### 7.4 能带的不可约表示（IRVSP）

[IRVSP](https://github.com/zjwang11/IRVSP) 读取 VASP 的 WAVECAR 与 OUTCAR（对称操作），给出每条能带在高对称点的不可约表示（支持 230 个空间群和 SOC 双群）：

```bash
# 1) 标准化晶胞：phonopy --symmetry --tolerance 0.01 -c POSCAR，再用 IRVSP 自带的 pos2aBR 转成 POSCAR_std
# 2) 用 POSCAR_std 做 VASP 计算（高对称 k 点，保留 WAVECAR）
# 3) 运行 IRVSP：空间群编号与能带范围
irvsp -sg 227 -nb 1 8 > outir
```

得到的不可约表示可用于能带标注、检查简并是否由对称性保护，以及在拓扑量子化学框架下判断拓扑性（与 Bilbao 晶体学服务器的相关工具配合）。

---

## 8. 扩展：点缺陷形成能与电荷修正

缺陷 X 在电荷态 q 下的形成能：

$$E_f[X^q]=E_\mathrm{tot}[X^q]-E_\mathrm{tot}[\text{bulk}]-\sum_i n_i\mu_i+q\,(E_F+E_\mathrm{VBM})+E_\mathrm{corr}$$

- $n_i$：加入（+）或移除（−）的第 i 种原子数；$\mu_i$ 为其化学势，受热力学稳定性约束（需要计算竞争相）；
- $E_F$ 以价带顶 $E_\mathrm{VBM}$ 为零点，在带隙内变化；
- 带电超胞：VASP 用 `NELECT` 指定电子总数，电荷 q = 中性电子数 − NELECT（例如 q = +1 时 `NELECT` 比默认值少 1）；
- $E_\mathrm{corr}$：周期性镜像电荷相互作用与势对齐修正，常用 Freysoldt–Neugebauer–Van de Walle（FNV）方法或 Kumagai–Oba 的扩展 FNV（eFNV，用原子位点势，适合各向异性体系），需要体相介电张量（[第14章](14-optics-dielectric.md)）和缺陷、体相的静电势（LOCPOT 或 OUTCAR 中的原子位点平均势）；
- 转变能级：$\varepsilon(q/q')=\dfrac{E_f[X^q;E_F=0]-E_f[X^{q'};E_F=0]}{q'-q}$；
- 带隙误差会直接影响转变能级，定量研究常用 HSE06（可调混合比例）。

工具：

- [SMTG-Bham/doped](https://github.com/SMTG-Bham/doped)：生成缺陷与电荷态、竞争相与化学势、解析 VASP 输出并自动施加 FNV/eFNV 修正、画形成能图和转变能级，与 ShakeNBreak（缺陷结构搜索）配合；
- pydefect：Kumagai 等开发的高通量缺陷计算工具，内置 eFNV 修正；
- 缺陷态的空间分布可以用 [第16章](16-charge-bonding.md) 的部分电荷（`IBAND`）可视化，缺陷超胞能带可以用 §5 的反折叠分析。

---

## 参数选择与收敛建议

| 任务 | 参数 | 推荐 | 判据 |
|---|---|---|---|
| Wannier90 | k 网格 | 半导体 6×6×6 起步，金属更密 | 冻结窗口内插值能带与 DFT 一致，无振荡 |
| Wannier90 | `NBANDS` | 外窗口上限以上仍有若干能带 | 每个 k 点外窗口内态数 ≥ `num_wann` |
| Wannier90 | `dis_froz_max` | E_F 或带隙中 + 1–2 eV（视目标子空间） | 不报 "more states in frozen window" |
| Wannier90 | `num_iter` / `dis_num_iter` | 200 / 500 起，按需增加 | `<-- DIS`、`<-- CONV` 收敛 |
| 有效质量 | 步长 | 0.005–0.02·2π/a | 两点差分与多点拟合一致 |
| 有效质量 | `EDIFF` | 1E-8 | 改为 1E-9 结果不变 |
| 有效质量 | 泛函 | 定性用 PBE，定量用 HSE/mBJ（+SOC） | 与实验或文献同等级计算对比 |
| 反折叠 | 超胞 k 点 | easyunfold 自动生成 | 完美超胞反折叠后与原胞能带重合 |
| 反折叠 | `NBANDS` | 覆盖感兴趣的导带范围 | 图中导带完整 |
| 费米面 | 非自洽网格 | Cu 原胞 20³–30³；IFermi 插值 8–10 倍 | 费米面形状与面积随网格不再变化 |
| 拓扑 | SOC + Wannier | 冻结窗口覆盖 E_F ± 1 eV | Wannier 能带在 E_F 附近与 DFT 误差 < 几 meV |
| 缺陷 | 超胞尺寸 | ≥ 10 Å 量级，并做尺寸外推或修正 | 修正后的形成能随超胞变化 < 约 0.1 eV |

---

## 常见问题与报错

| 现象 / 报错信息 | 原因 | 解决办法 |
|---|---|---|
| 设了 `LWANNIER90 = .TRUE.` 但没有 `wannier90.amn/.mmn/.eig` | VASP 未用 `-DVASP2WANNIER90` 编译 | 按 §3.1 重新编译并链接 libwannier.a |
| 链接时报找不到 `wannier_setup` 等符号 | Wannier90 库版本或路径不对；Wannier90 3.x 库编译成了 MPI 版本 | 使用 3.1.0、按串行方式 `make lib`，检查 `WANNIER90_ROOT` |
| Wannier90 报 `dis_windows: More states in the frozen window than target WFs` | 冻结窗口内的态比 `num_wann` 多 | 降低 `dis_froz_max`，或增加 Wannier 函数个数 |
| Wannier90 报外窗口内态数少于目标 Wannier 函数数 | `dis_win_max` 太低或 `NBANDS` 不够 | 提高 `dis_win_max`，增大 `NBANDS` 后重算 VASP |
| Wannier90 报 k 点/`mp_grid` 不一致 | VASP 用了对称约化（`ISYM` 未设 -1），或 KPOINTS 不是完整均匀网格 | `ISYM = -1`，KPOINTS 用均匀网格（Γ 中心） |
| 插值能带在冻结窗口内有尖峰或振荡 | k 网格太稀；投影不合适；spread 未收敛 | 加密 k 网格；修改投影；增大迭代步数 |
| 非共线计算的 Wannier 结果不对称或错误 | 旧版本非共线 + 对称性的 amn/mmn 缺陷（Known issue 19，6.4.3 修复） | `ISYM = -1`，升级 VASP |
| `KPOINTS_WAN` 插值能带在立方晶体中明显错误 | Known issue 105 | 改用 Wannier90 的 `bands_plot`，或加密 KPOINTS 并与 `KPOINTS_OPT` 对比 |
| `WANNIER90_WIN` 内容被截断 | 超过 8192 字符（Known issue 64） | 改用独立的 wannier90.win 文件 |
| 有效质量为极大值或符号异常 | 选错能带（简并交叉）；步长太小被噪声主导；k0 不是极值点 | 用 `--band` 指定能带，比较多个 `--npoints`；先确认 k0 处确是极值 |
| 有效质量与文献差一个数量级 | 单位换算漏了 2π；k 用了分数坐标 | k 必须换成 Å⁻¹（含 2π），用 m*/m0 = 3.81/a |
| 反折叠图几乎没有权重 | `--matrix` 与实际超胞不符；KPOINTS 没换成 `KPOINTS_easyunfold` | 检查矩阵（`easyunfold unfold status`）；重新做非自洽 |
| `unfold calculate` 读 WAVECAR 失败 | Gamma-only 或非共线 WAVECAR 未加选项；WAVECAR 不完整 | 加 `--gamma` / `--ncl`；确认 VASP 正常结束 |
| IFermi 提示没有费米面 | 半导体，E_F 在带隙中 | 用 `--mu` 移到导带或价带 |
| 费米面锯齿、颈部断开 | 网格太粗 | 加密非自洽网格或增大 `--interpolation-factor` |
| 带电缺陷形成能随超胞剧烈变化 | 没有做镜像电荷修正或势对齐 | 使用 FNV/eFNV 修正（doped、pydefect） |

---

## 练习

1. **Wannier 窗口**：把 Si 示例的 `dis_froz_max` 分别设为 E_F+1、E_F+2、E_F+4 eV，观察 Wannier90 的报错或插值能带的变化。
2. **只构造价带**：只用 4 条价带（`NBANDS = 4`、`NUM_WANN = 4`、键中心 s 型投影 `f=0.125,0.125,0.125:s` 等），比较键中心 Wannier 函数与 sp3 函数的中心和 spread。
3. **HSE 插值**：用 §3.7 的设置得到 Si 的 HSE06 插值能带，并与第6章零权重 k 点法的结果对比。
4. **有效质量与泛函**：分别用 PBE 和 HSE06（零权重 k 点或 KPOINTS_OPT）计算 GaAs 的 Γ 点电子有效质量，与实验值 0.067 m0 比较；再加入 SOC，观察空穴质量的变化。
5. **反折叠检验**：用完美的 3×3×3 Si 超胞做一次反折叠，确认得到的 EBS 与原胞能带完全重合，再与空位超胞对比。
6. **费米面收敛**：对 Cu 分别用 11³、17³、25³ 的非自洽网格，比较 `ifermi info` 给出的费米面面积。

---

## 参考资料

**官方文档**

- VASP Wiki：[LWANNIER90](https://www.vasp.at/wiki/index.php/LWANNIER90)、[LWANNIER90_RUN](https://www.vasp.at/wiki/index.php/LWANNIER90_RUN)、[WANNIER90_WIN](https://www.vasp.at/wiki/index.php/WANNIER90_WIN)、[NUM_WANN](https://www.vasp.at/wiki/index.php/NUM_WANN)、[LWRITE_MMN_AMN](https://www.vasp.at/wiki/index.php/LWRITE_MMN_AMN)、[LWRITE_UNK](https://www.vasp.at/wiki/index.php/LWRITE_UNK)、[LSCDM](https://www.vasp.at/wiki/index.php/LSCDM)、[LOCPROJ](https://www.vasp.at/wiki/index.php/LOCPROJ)、[KPOINTS_WAN](https://www.vasp.at/wiki/index.php/KPOINTS_WAN)
- VASP Wiki：[Constructing Wannier orbitals](https://www.vasp.at/wiki/index.php/Constructing_Wannier_orbitals)、[Precompiler options](https://www.vasp.at/wiki/index.php/Precompiler_options)、[Makefile.include](https://www.vasp.at/wiki/index.php/Makefile.include)、[KPOINTS](https://www.vasp.at/wiki/index.php/KPOINTS)、[Known issues](https://www.vasp.at/wiki/index.php/Known_issues)
- Wannier90 官网与文档：<https://wannier.org/>
- easyunfold 文档：<https://smtg-bham.github.io/easyunfold/>
- IFermi 文档：<https://fermisurfaces.github.io/IFermi/>
- effmass 文档：<https://effmass.readthedocs.io/en/latest/>
- VASPKIT 教程（有效质量、反折叠、费米面）：<https://vaspkit.com/tutorials.html>

**相关 GitHub 仓库**（只链接，不转载代码）

- [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90)：Wannier90 源码（VASP ≥ 6.2 链接 2.x/3.x）
- [JMSkelton/VASP-Examples](https://github.com/JMSkelton/VASP-Examples)：Si 的 Wannier90 插值 HSE06 能带等输入示例
- [Chengcheng-Xiao/VASP2WAN90_v2_fix](https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix)：VASP 5.4.4 与 Wannier90 2.x 接口问题的修复
- [afonari/emc](https://github.com/afonari/emc)、[lucydot/effmass](https://github.com/lucydot/effmass)：有效质量
- [SMTG-Bham/easyunfold](https://github.com/SMTG-Bham/easyunfold)、[QijingZheng/VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding)：能带反折叠
- [fermisurfaces/IFermi](https://github.com/fermisurfaces/IFermi)：费米面
- [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar)：投影能带、自旋纹理、费米面
- [Infant83/VASPBERRY](https://github.com/Infant83/VASPBERRY)：从 WAVECAR 计算 Berry 曲率与陈数
- [zjwang11/IRVSP](https://github.com/zjwang11/IRVSP)：能带不可约表示
- [SMTG-Bham/doped](https://github.com/SMTG-Bham/doped)：点缺陷计算与修正
- [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo)：k 路径生成与 `sumo-bandstats`
- 更多工具见 [Resources/github-repos.md](../../Resources/github-repos.md)

**关键文献**

- N. Marzari, D. Vanderbilt, *Phys. Rev. B* 56, 12847 (1997). DOI: [10.1103/PhysRevB.56.12847](https://doi.org/10.1103/PhysRevB.56.12847)
- I. Souza, N. Marzari, D. Vanderbilt, *Phys. Rev. B* 65, 035109 (2001). DOI: [10.1103/PhysRevB.65.035109](https://doi.org/10.1103/PhysRevB.65.035109)
- N. Marzari et al., *Rev. Mod. Phys.* 84, 1419 (2012). DOI: [10.1103/RevModPhys.84.1419](https://doi.org/10.1103/RevModPhys.84.1419)
- G. Pizzi et al., *J. Phys.: Condens. Matter* 32, 165902 (2020)（Wannier90 v3）. DOI: [10.1088/1361-648X/ab51ff](https://doi.org/10.1088/1361-648X/ab51ff)
- C. Franchini et al., *J. Phys.: Condens. Matter* 24, 235602 (2012)（VASP2WANNIER90）. DOI: [10.1088/0953-8984/24/23/235602](https://doi.org/10.1088/0953-8984/24/23/235602)
- A. Damle, L. Lin, *Multiscale Model. Simul.* 16, 1392 (2018)（SCDM）. DOI: [10.1137/17M1129696](https://doi.org/10.1137/17M1129696)
- I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, *J. Appl. Phys.* 89, 5815 (2001)（III-V 能带参数）. DOI: [10.1063/1.1368156](https://doi.org/10.1063/1.1368156)
- L. D. Whalley et al., *J. Open Source Softw.* 3, 797 (2018)（effmass）. DOI: [10.21105/joss.00797](https://doi.org/10.21105/joss.00797)
- V. Popescu, A. Zunger, *Phys. Rev. Lett.* 104, 236403 (2010); *Phys. Rev. B* 85, 085201 (2012). DOI: [10.1103/PhysRevLett.104.236403](https://doi.org/10.1103/PhysRevLett.104.236403)、[10.1103/PhysRevB.85.085201](https://doi.org/10.1103/PhysRevB.85.085201)
- B. Zhu et al., *J. Open Source Softw.* 9, 5974 (2024)（easyunfold）. DOI: [10.21105/joss.05974](https://doi.org/10.21105/joss.05974)
- A. M. Ganose et al., *J. Open Source Softw.* 6, 3089 (2021)（IFermi）. DOI: [10.21105/joss.03089](https://doi.org/10.21105/joss.03089)
- M. Kawamura, *Comput. Phys. Commun.* 239, 197 (2019)（FermiSurfer）. DOI: [10.1016/j.cpc.2019.01.017](https://doi.org/10.1016/j.cpc.2019.01.017)
- T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* 74, 1674 (2005). DOI: [10.1143/JPSJ.74.1674](https://doi.org/10.1143/JPSJ.74.1674)
- Q. Wu et al., *Comput. Phys. Commun.* 224, 405 (2018)（WannierTools）. DOI: [10.1016/j.cpc.2017.09.033](https://doi.org/10.1016/j.cpc.2017.09.033)
- D. Gresch et al., *Phys. Rev. B* 95, 075146 (2017)（Z2Pack）. DOI: [10.1103/PhysRevB.95.075146](https://doi.org/10.1103/PhysRevB.95.075146)
- J. Gao, Q. Wu, C. Persson, Z. Wang, *Comput. Phys. Commun.* 261, 107760 (2021)（IRVSP）. DOI: [10.1016/j.cpc.2020.107760](https://doi.org/10.1016/j.cpc.2020.107760)
- K. Ishizaka et al., *Nat. Mater.* 10, 521 (2011)（BiTeI 巨 Rashba 劈裂）. DOI: [10.1038/nmat3051](https://doi.org/10.1038/nmat3051)
- C. Freysoldt, J. Neugebauer, C. G. Van de Walle, *Phys. Rev. Lett.* 102, 016402 (2009). DOI: [10.1103/PhysRevLett.102.016402](https://doi.org/10.1103/PhysRevLett.102.016402)
- Y. Kumagai, F. Oba, *Phys. Rev. B* 89, 195205 (2014). DOI: [10.1103/PhysRevB.89.195205](https://doi.org/10.1103/PhysRevB.89.195205)
- C. Freysoldt et al., *Rev. Mod. Phys.* 86, 253 (2014)（点缺陷计算综述）. DOI: [10.1103/RevModPhys.86.253](https://doi.org/10.1103/RevModPhys.86.253)
- S. R. Kavanagh et al., *J. Open Source Softw.* 9, 6433 (2024)（doped）. DOI: [10.21105/joss.06433](https://doi.org/10.21105/joss.06433)

---

[← 上一章](16-charge-bonding.md) | [返回目录](../README.md) | [下一章 →](18-postprocess-troubleshooting.md)
