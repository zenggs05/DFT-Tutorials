# 第8章 色散校正、自旋极化与 DFT+U（Dispersion, Spin Polarization & DFT+U in CP2K）

本章讲 CP2K 中三类常一起出现的“修正”：范德华色散校正（DFT-D3、DFT-D3(BJ)、DFT-D4，以及 rVV10、vdW-DF 等非局域泛函），自旋极化计算（`UKS`/`LSD`、`MULTIPLICITY`、`CHARGE`、`MAGNETIZATION`、`&BS`），以及过渡金属氧化物常用的 DFT+U（`PLUS_U_METHOD`、`&DFT_PLUS_U`、`U_MINUS_J`）。示例包括苯二聚体相互作用能（含 BSSE 校正）、O2 三重态、NiO 反铁磁 PBE+U（2026.2 起支持 k 点）以及石墨的 PBE-D4 晶胞优化，并说明各版本之间的默认值变化。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 在 `&XC/&VDW_POTENTIAL` 中正确开启 D3、D3(BJ)、D4，知道 `TYPE`、`PARAMETER_FILE_NAME`、`REFERENCE_FUNCTIONAL` 的作用以及默认值在 2025.2 → 2026.2 之间的变化；
- 计算分子间相互作用能，并用 `RUN_TYPE BSSE` 做 counterpoise 校正；
- 会写 rVV10、vdW-DF 类非局域色散泛函的输入；
- 为自由基、三重态分子、磁性固体设置 `UKS`、`MULTIPLICITY`、`CHARGE`，读懂自旋布居和 $\langle S^2\rangle$；
- 用多个 `&KIND`（`Fe1`/`Fe2` + `ELEMENT Fe`）和 `MAGNETIZATION` 构造铁磁与反铁磁初猜，必要时用 `&BS` 精细控制轨道占据；
- 写出 DFT+U 输入，理解 `U_MINUS_J` 的单位陷阱、Mulliken 与 Löwdin 投影的区别、U 斜坡（`U_RAMPING`），以及 2026.2 新增的 k 点 DFT+U。

| 体系 | 需要色散校正 | 需要自旋极化 | 需要 DFT+U |
|---|---|---|---|
| 分子晶体、层状材料、物理吸附、液态水 | 是 | 一般否 | 否 |
| 自由基、O2、卡宾、含奇数电子的离子 | 视情况 | 是 | 否 |
| Fe、Co、Ni 金属 | 通常不必 | 是（铁磁） | 一般不用 |
| NiO、MnO、FeO、CoO 等 Mott 绝缘体 | 视情况 | 是（多为反铁磁） | 通常需要 |
| TiO2 中的极化子、CeO2 的 Ce 4f | 视情况 | 是 | 常用 |

示例输入集：[../examples/08-dispersion-spin-dftu/README.md](../examples/08-dispersion-spin-dftu/README.md)。

---

## 2. 原理速览

### 2.1 原子对色散校正：D2 / D3 / D3(BJ) / D4

半局域泛函（LDA、GGA）描述不了长程电子关联引起的色散作用。Grimme 型校正在 DFT 总能量上加一个原子对（必要时含三体）能量项：

$$
E_\mathrm{disp} = -\sum_{A<B}\ \sum_{n=6,8} s_n\,\frac{C_n^{AB}}{r_{AB}^{\,n}}\ f_{\mathrm{damp},n}(r_{AB})
$$

- **D2**：$C_6$ 为固定的原子参数，现已很少使用；
- **D3**：$C_6^{AB}$ 依赖配位数（coordination number），零阻尼函数（zero damping），参数为 $(s_6, s_{r,6}, s_8)$；
- **D3(BJ)**：Becke–Johnson 有理阻尼 $f_n = r^n/\big(r^n + (a_1 R_0^{AB} + a_2)^n\big)$，短程行为更合理，参数为 $(s_6, a_1, s_8, a_2)$，是目前最常用的选择；
- **D4**：$C_6$ 进一步依赖原子部分电荷（由 EEQ 模型得到），对离子性体系和金属有机体系更可靠。

阻尼参数与泛函一一对应，这就是 `REFERENCE_FUNCTIONAL` 存在的原因：PBE 和 B3LYP 用的是不同的 $(s_6,a_1,s_8,a_2)$。

### 2.2 非局域 vdW 泛函

vdW-DF、vdW-DF2、optB88-vdW、rVV10 等把非局域关联写成密度的双重积分：

$$
E_c^\mathrm{nl} = \frac12 \iint \rho(\mathbf r)\,\Phi(\mathbf r,\mathbf r')\,\rho(\mathbf r')\,\mathrm d\mathbf r\,\mathrm d\mathbf r'
$$

核函数 $\Phi$ 预先制成表格（CP2K data 目录下的 `vdW_kernel_table.dat`、`rVV10_kernel_table.dat`）。这类方法不需要原子类型参数，但要与特定的交换泛函搭配，计算量也比 D3 大。

### 2.3 自旋极化

自旋非限制 Kohn–Sham（UKS，CP2K 中别名 `LSD`、`SPIN_POLARIZED`）分别求解 α、β 两套轨道。多重度 $M = 2S+1$，α、β 电子数为

$$
N_\alpha = \frac{N + M - 1}{2},\qquad N_\beta = \frac{N - M + 1}{2}
$$

这正是 CP2K 源码中分配电子数的方式。单行列式 UKS 波函数的 $\langle S^2\rangle$ 与理想值 $S(S+1)$ 的偏差用来衡量自旋污染。

### 2.4 DFT+U（Dudarev 形式）

对局域的 d、f 电子，在 DFT 能量上加 Hubbard 型修正，CP2K 采用 Dudarev 的旋转不变简化形式：

$$
E_{U} = \frac{U_\mathrm{eff}}{2}\sum_{I}\sum_{\sigma}\mathrm{Tr}\!\left[\mathbf n^{I\sigma} - \mathbf n^{I\sigma}\mathbf n^{I\sigma}\right],\qquad U_\mathrm{eff} = U - J
$$

其中 $\mathbf n^{I\sigma}$ 是原子 $I$ 上指定角动量壳层（如 d 壳层，$l=2$）的占据矩阵。CP2K 用原子轨道的 Mulliken 或 Löwdin 布居构造 $\mathbf n$。投影方式不同，$\mathbf n$ 的数值就不同，所以**同样的 U 值在 CP2K（高斯基 + Mulliken 投影）和 VASP（PAW 投影）中并不等价**，不能直接照搬。

---

## 3. DFT-D3(BJ)：苯二聚体相互作用能

### 3.1 输入文件

示例 [../examples/08-dispersion-spin-dftu/benzene_dimer_d3bj.inp](../examples/08-dispersion-spin-dftu/benzene_dimer_d3bj.inp)。用 `@SET SYS` 选择计算二聚体（`dimer`）还是单体（`a`、`b`），用 `@SET VDW` 开关色散校正：

```text
@SET SYS dimer
@SET VDW D3BJ

&GLOBAL
  PROJECT bz_${SYS}_${VDW}
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 30
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
@IF ${VDW} == D3BJ
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF [angstrom] 15.0
          &PRINT_DFTD ON
          &END PRINT_DFTD
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
@ENDIF
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 20.0 20.0 20.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME benzene_${SYS}.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND C
      BASIS_SET TZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      BASIS_SET TZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### 3.2 关键字说明

| 关键字 | 本例取值 | 说明 |
|---|---|---|
| `DISPERSION_FUNCTIONAL` | `PAIR_POTENTIAL` | 原子对校正（D2/D3/D4）；`NON_LOCAL` 为非局域泛函（第 5 节） |
| `TYPE` | `DFTD3(BJ)` | 可选 `DFTD2`、`DFTD3`、`DFTD3(BJ)`、`DFTD4`（trunk 还有 `DFTD2(CHG)`）。**默认值随版本变化**：2025.2 为 `DFTD3`（零阻尼），2026.2 为 `DFTD3(BJ)`，务必显式写出 |
| `PARAMETER_FILE_NAME` | `dftd3.dat`（默认） | D2/D3 参考 $C_6$ 数据文件，在 `CP2K_DATA_DIR` 下查找，也可写完整路径；D4 不使用此文件 |
| `REFERENCE_FUNCTIONAL` | `PBE` | 选择该泛函的阻尼参数。2026.2 起（PR #5296）对 PBE、PBE0、B3LYP、BLYP、BP86、TPSS 等常见泛函可自动识别；为兼容旧版本和组合泛函，建议始终显式写出 |
| `R_CUTOFF` | 15 Å | 对势截断范围（别名 `D3_CUTOFF`），默认约 10.58 Å；周期体系可适当增大并检查收敛 |
| `CALCULATE_C9_TERM` | 未开（默认 F） | D3 的三体（ATM）项，对致密分子晶体有影响 |
| `D3BJ_SCALING` | 未写 | 手动给出 $(s_6, a_1, s_8, a_2)$，覆盖 `REFERENCE_FUNCTIONAL` 的选择；例如 PBE-D3(BJ) 的 Grimme 参数为 `1.0 0.4289 0.7875 4.4407` |
| `D3_EXCLUDE_KIND` | 未写 | 排除某些 KIND（按 KIND 编号）不计色散 |
| `&PRINT_DFTD` | `ON` | 打印色散计算的详细信息（默认 HIGH 级别才打印） |

> [!NOTE]
> 坐标文件已经把二聚体放在 20 Å 盒子的中心，单体文件使用与二聚体**完全相同**的坐标。这里不能加 `&CENTER_COORDINATES`，否则单体会被各自平移到盒子中心，虽然对孤立体系的能量没有影响，但会让结构对比和后续 BSSE 计算变得混乱。

### 3.3 运行与计算相互作用能

```bash
cd CP2K/examples/08-dispersion-spin-dftu
for s in dimer a b; do
  for v in D3BJ NONE; do
    sed -e "s/^@SET SYS .*/@SET SYS $s/" -e "s/^@SET VDW .*/@SET VDW $v/" \
        benzene_dimer_d3bj.inp > bz_${s}_${v}.inp
    mpirun -np 8 cp2k.psmp -i bz_${s}_${v}.inp -o bz_${s}_${v}.out
  done
done
grep "Dispersion energy:" bz_*_D3BJ.out
for v in D3BJ NONE; do
  python - "$v" <<'EOF'
import re, sys
v = sys.argv[1]
def etot(f):
    return float(re.findall(r"ENERGY\| Total FORCE_EVAL.*?:\s+(-?\d+\.\d+)", open(f).read())[-1])
e = etot(f"bz_dimer_{v}.out") - etot(f"bz_a_{v}.out") - etot(f"bz_b_{v}.out")
print(f"{v:5s}  E_int = {e:.6f} Ha = {e*627.509:.2f} kcal/mol")
EOF
done
```

能量分解中会出现一行 `Dispersion energy:`（单位 hartree）。期望现象：纯 PBE 的相互作用能很弱，甚至为正（排斥）；加上 D3(BJ) 后变为明显的负值（吸引）。作为参照，S22 基准集中苯二聚体（平行错位、T 型构型）的 CCSD(T)/CBS 相互作用能约为 −2.7 kcal/mol（Jurečka 等，DOI 见参考资料）。本例的几何是按理想键长构造的平行错位构型，不是 S22 的参考几何，只能做定性比较。

### 3.4 BSSE：counterpoise 校正

高斯基组计算分子间相互作用能会有基组重叠误差（BSSE）：二聚体中每个单体都“借用”了对方的基函数，使结合能偏大。CP2K 提供 `RUN_TYPE BSSE` 自动完成 Boys–Bernardi counterpoise 计算，示例 [../examples/08-dispersion-spin-dftu/benzene_dimer_bsse.inp](../examples/08-dispersion-spin-dftu/benzene_dimer_bsse.inp)：

```text
&GLOBAL
  PROJECT bz_dimer_bsse
  RUN_TYPE BSSE
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &BSSE
    &FRAGMENT
      LIST 1..12
    &END FRAGMENT
    &FRAGMENT
      LIST 13..24
    &END FRAGMENT
  &END BSSE
  ...
&END FORCE_EVAL
```

运行后在输出中查找结果块：

```bash
grep -A12 "BSSE RESULTS" bz_dimer_bsse.out
grep "BSSE-free interaction energy" bz_dimer_bsse.out     # hartree
```

TZVP 基组的 BSSE 通常比 DZVP 小得多；如果 BSSE 校正量与相互作用能本身相当，说明基组不够，应换更大的基组。

---

## 4. DFT-D4 与版本差异

D4 在 CP2K 中通过外部库 `libdftd4` 实现。确认你的版本是否可用：

```bash
cp2k.psmp -v | grep -i -o "libdftd4"      # must print libdftd4
```

输入写法（示例 [../examples/08-dispersion-spin-dftu/graphite_d4_cell_opt.inp](../examples/08-dispersion-spin-dftu/graphite_d4_cell_opt.inp) 的核心部分）：

```text
&VDW_POTENTIAL
  DISPERSION_FUNCTIONAL PAIR_POTENTIAL
  &PAIR_POTENTIAL
    TYPE DFTD4
    REFERENCE_FUNCTIONAL PBE
  &END PAIR_POTENTIAL
&END VDW_POTENTIAL
```

| 关键字 | 说明 |
|---|---|
| `TYPE DFTD4` | 参数取自 libdftd4，不读 `dftd3.dat` |
| `REFERENCE_FUNCTIONAL` | D4 的可用泛函名见 dftd4 源码中的参数表；2026.2 起常见泛函可自动识别（含 SCAN、r2SCAN） |
| `D4_SCALING` | 手动给出 $(s_6,a_1,s_8,a_2)$ |
| `D4_CUTOFF`、`D4_CN_CUTOFF`、`R_CUTOFF` | 二体项、配位数和三体项（D4 中 `R_CUTOFF` 用于三体项）的截断 |
| `D4_REFERENCE_CODE` | 是否完全交给外部库计算。默认值随版本变化：2025.2 为 `T`；2026.2 为 `F`（使用 CP2K 内置实现，2026.2 手册说明外部库不支持 MPI 并行，因此 MPI 版本不建议打开）；trunk 手册又改为 `T`。无论哪种设置都需要链接 libdftd4 |

> [!TIP]
> 编译时没有 libdftd4 的话，把示例中的 `@SET VDW D4` 改成 `D3BJ`（文件中已用 `@IF` 写好两种分支）。D3(BJ) 对石墨层间距的描述同样不错，可作为对照。

### 4.1 石墨的 PBE-D4 晶胞优化

AB 堆垛石墨（4 原子），初始 c = 7.20 Å（实验值 6.711 Å，层间距 3.355 Å），k 点 12×12×4，对角化 + Fermi–Dirac 展宽（石墨是半金属）：

```text
&SCF
  ADDED_MOS 20
  &DIAGONALIZATION
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
    ALPHA 0.2
    NBROYDEN 8
  &END MIXING
  &SMEAR ON
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
&END SCF
...
&CELL_OPT
  OPTIMIZER BFGS
  EXTERNAL_PRESSURE [bar] 1.01325
  PRESSURE_TOLERANCE [bar] 50.0
  KEEP_ANGLES T
&END CELL_OPT
```

检查：

```bash
grep "CELL| Vector c" graphite_cellopt_D4.out | tail -1
grep "Dispersion energy:" graphite_cellopt_D4.out | tail -1
```

最终 c 应明显小于初始的 7.20 Å，并接近实验值（具体数值请以计算为准）。作为对比，去掉色散校正（纯 PBE）时层间几乎没有结合，c 会弛豫到远大于实验值的位置，甚至无法收敛到合理极小。晶胞优化的一般设置见 [第7章](07-geo-cell-opt.md)。

---

## 5. 非局域 vdW 泛函：rVV10 与 vdW-DF

### 5.1 rVV10

rVV10 一般与 rPW86 交换 + PBE 关联搭配（需要 CP2K 链接 LibXC）：

```text
&XC
  &XC_FUNCTIONAL
    &GGA_X_RPW86
    &END GGA_X_RPW86
    &GGA_C_PBE
    &END GGA_C_PBE
  &END XC_FUNCTIONAL
  &VDW_POTENTIAL
    DISPERSION_FUNCTIONAL NON_LOCAL
    &NON_LOCAL
      TYPE RVV10
      PARAMETERS 6.3 0.0093
      KERNEL_FILE_NAME rVV10_kernel_table.dat
      CUTOFF 100
    &END NON_LOCAL
  &END VDW_POTENTIAL
&END XC
```

### 5.2 vdW-DF 家族

vdW-DF、vdW-DF2 及 optB88-vdW、optPBE-vdW 等使用 `TYPE DRSLL`（vdW-DF 型核函数）或 `TYPE LMKLL`（vdW-DF2 型核函数），核函数文件为 `vdW_kernel_table.dat`。以 optB88-vdW 为例（需要 LibXC）：

```text
&XC
  &XC_FUNCTIONAL
    &GGA_X_OPTB88_VDW
    &END GGA_X_OPTB88_VDW
    &PW92
    &END PW92
  &END XC_FUNCTIONAL
  &VDW_POTENTIAL
    DISPERSION_FUNCTIONAL NON_LOCAL
    &NON_LOCAL
      TYPE DRSLL
      KERNEL_FILE_NAME vdW_kernel_table.dat
      CUTOFF 100
    &END NON_LOCAL
  &END VDW_POTENTIAL
&END XC
```

CP2K 源码的 `data/xc_section/` 目录提供了 `vdW-DF1.sec`、`vdW-DF2.sec`、`vdW-DF-optB88.sec`、`vdW-DF-optPBE.sec`、`vdW-DF-cx.sec` 等预设的 `&XC` 段，可以用预处理指令 `@XCTYPE vdW-DF-optB88.sec` 直接引入，避免手写出错。该指令在当前目录找不到文件时会改为查找 `xc_section/vdW-DF-optB88.sec`；如果你的安装找不到它，把对应的 `.sec` 文件从 CP2K 的 `data/xc_section/` 复制到工作目录即可。

| 关键字（`&NON_LOCAL`） | 说明 |
|---|---|
| `TYPE` | `DRSLL`（默认）、`LMKLL`、`RVV10` |
| `KERNEL_FILE_NAME` | DRSLL/LMKLL 用 `vdW_kernel_table.dat`，RVV10 用 `rVV10_kernel_table.dat` |
| `PARAMETERS` | rVV10 的 $b$ 和 $C$，默认 `6.3 0.0093` |
| `CUTOFF` | 计算非局域项所用 FFT 网格的截断（Ry），独立于 `&MGRID`，需要单独测试收敛 |
| `SCALE` | 非局域能量的缩放因子 |

> [!NOTE]
> 非局域泛函的交换部分必须与核函数配套：vdW-DF 用 revPBE 交换，vdW-DF2 用 rPW86 交换，optB88-vdW 用 optB88 交换，且这一家族的局域关联部分只用 LDA。随意拿 PBE 交换加一个 vdW 核函数得到的不是任何已发表的泛函。优先参考 `data/xc_section/` 中的预设文件。

---

## 6. 自旋极化计算：UKS、MULTIPLICITY 与 CHARGE

### 6.1 O2 三重态

O2 的基态是三重态 $X\,^3\Sigma_g^-$（两个未配对电子）。示例 [../examples/08-dispersion-spin-dftu/o2_triplet_uks.inp](../examples/08-dispersion-spin-dftu/o2_triplet_uks.inp) 的核心设置：

```text
@SET MULT 3
@SET PROJ o2_mult${MULT}
...
&DFT
  BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
  POTENTIAL_FILE_NAME POTENTIAL_UZH
  UKS
  CHARGE 0
  MULTIPLICITY ${MULT}
  ...
&END DFT
...
&GEO_OPT
  OPTIMIZER BFGS
  MAX_ITER 100
  MAX_FORCE 1.0E-4
&END GEO_OPT
```

| 关键字（`&DFT` 下） | 说明 |
|---|---|
| `UKS` | 自旋非限制计算，别名 `LSD`、`SPIN_POLARIZED`、`UNRESTRICTED_KOHN_SHAM` |
| `MULTIPLICITY` | $2S+1$，别名 `MULTIP`。不写时：偶数电子取 1，奇数电子取 2 |
| `CHARGE` | 体系总电荷，例如 `CHARGE -1` 表示多一个电子 |
| `ROKS` | 限制性开壳层 Kohn–Sham，可替代 UKS 处理高自旋态 |
| `RELAX_MULTIPLICITY` | 不强制按 `MULTIPLICITY` 分配 α/β 占据，而按能量（Aufbau）占据，单位 Ha；需要 `ADDED_MOS` > 0，因此不适用于 OT |

运行三重态和闭壳层（`MULT 1`）两次并比较：

```bash
mpirun -np 4 cp2k.psmp -i o2_triplet_uks.inp -o o2_mult3.out
sed 's/^@SET MULT .*/@SET MULT 1/' o2_triplet_uks.inp > o2_mult1.inp
mpirun -np 4 cp2k.psmp -i o2_mult1.inp -o o2_mult1.out
grep "ENERGY| Total FORCE_EVAL" o2_mult*.out | awk '{print $1, $NF}'
grep "Ideal and single determinant S\*\*2" o2_mult3.out | tail -1
grep "Integrated absolute spin density" o2_mult3.out | tail -1
```

检查要点：

- 三重态能量应低于 `MULT 1` 的解。注意 UKS 的 `MULT 1` 得到的是闭壳层（或自旋对称性破缺）解，**不是**真正的 $a\,^1\Delta_g$ 激发态（实验上 $a\,^1\Delta_g$ 比基态高约 0.98 eV）；单重态需要多组态或专门的方法处理；
- `Ideal and single determinant S**2 :` 一行给出理想值 2.0 和实际的 $\langle S^2\rangle$，两者接近说明自旋污染小；
- Mulliken 布居分析（`PRINT_LEVEL MEDIUM` 时输出）最后一列 `Spin moment`，两个 O 各约 1；
- 优化后 O–O 键长：实验值 1.2075 Å，PBE 通常略长（约 1.22–1.23 Å，近似范围）。

### 6.2 奇数电子与 LSD 报错

体系电子数为奇数而没有开 UKS 时，CP2K 中止并报：

```text
Use the LSD option for an odd number of electrons
```

解决：在 `&DFT` 中加 `UKS`，`MULTIPLICITY` 设为 2（或更高的正确值）。如果体系本应是闭壳层（如中性有机分子），出现此错多半是结构中少了或多了 H 原子，或 `CHARGE` 写错。另一个相关报错 `LSD: too few electrons for this multiplicity` 表示 `MULTIPLICITY` 超过了电子数允许的最大值。

> [!TIP]
> 自旋态不确定时（如过渡金属配合物），对每个可能的多重度分别计算，比较能量。不要只算一个多重度就下结论。

---

## 7. 磁性初猜：MAGNETIZATION、多 KIND 与 &BS

### 7.1 MAGNETIZATION

`&KIND` 中的 `MAGNETIZATION` 只影响原子初猜（`SCF_GUESS ATOMIC`）：给该原子加 `MAGNETIZATION/2` 个 α 电子、减 `MAGNETIZATION/2` 个 β 电子。SCF 收敛后的磁矩由体系本身决定，但好的初猜决定了收敛到哪个磁性解。

铁磁 bcc Fe（示意片段）：

```text
&DFT
  UKS
  ! a ferromagnetic metal: let the moment relax with smearing (see 7.3)
&END DFT
&KIND Fe
  ELEMENT Fe
  BASIS_SET DZVP-MOLOPT-PBE-GTH-q16
  POTENTIAL GTH-PBE-q16
  MAGNETIZATION 2.5
&END KIND
```

### 7.2 反铁磁：两个 KIND

同一元素自旋方向不同的原子必须属于不同的 KIND，`ELEMENT` 指明元素，`MAGNETIZATION` 取相反符号：

```text
&KIND Fe1
  ELEMENT Fe
  BASIS_SET DZVP-MOLOPT-PBE-GTH-q16
  POTENTIAL GTH-PBE-q16
  MAGNETIZATION  4.0
&END KIND
&KIND Fe2
  ELEMENT Fe
  BASIS_SET DZVP-MOLOPT-PBE-GTH-q16
  POTENTIAL GTH-PBE-q16
  MAGNETIZATION -4.0
&END KIND
```

坐标中（`&COORD` 或 xyz 文件第一列）对应原子写成 `Fe1`、`Fe2`。总磁矩为零的反铁磁态在 `&DFT` 中设 `UKS` 和 `MULTIPLICITY 1`；亚铁磁或铁磁态则把 `MULTIPLICITY` 设为 $2S_\mathrm{tot}+1$。

> [!WARNING]
> 只写一个 `&KIND Fe` 并期望 SCF 自己“找到”反铁磁解，几乎总会得到铁磁或非磁解。xyz 文件里用 `Fe1`/`Fe2` 这类标签后，部分可视化软件识别不了元素，可以另存一份只含元素符号的文件用于作图。

### 7.3 金属磁体与展宽

金属磁性体系要用对角化 + 展宽。展宽时 α/β 电子数可以在 SCF 中调整，`&SMEAR` 中的 `FIXED_MAGNETIC_MOMENT` 控制总磁矩：

| `FIXED_MAGNETIC_MOMENT` | 含义 |
|---|---|
| 负值（默认 −100） | 允许磁矩变化 |
| −1 | 保持 α、β 电子数为整数 |
| ≥ 0 | 固定 $n_\uparrow - n_\downarrow$ 为该值 |

### 7.4 &BS：精细控制初猜占据

`&KIND/&BS` 按角动量通道增减初猜中的电子（只对 `SCF_GUESS ATOMIC` 有效），常用于给出离子性、高自旋或破缺对称性（broken symmetry）初猜。`&ALPHA`、`&BETA` 子段中 `N` 为主量子数、`L` 为角量子数、`NEL` 为增减的电子数，三者按位置一一对应。以 NiO 中的 Ni²⁺（d⁸，高自旋 S = 1）为例，假设原子初猜为 3d⁸4s²（α、β 各 4s¹3d⁴），自旋向上的 Ni 可以写成（思路示意）：

```text
&KIND Ni1
  ELEMENT Ni
  BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
  POTENTIAL GTH-PBE-q18
  &BS
    &ALPHA
      N   4  3
      L   0  2
      NEL -1 1       ! alpha: 4s1 -> 4s0, 3d4 -> 3d5
    &END ALPHA
    &BETA
      N   4  3
      L   0  2
      NEL -1 -1      ! beta : 4s1 -> 4s0, 3d4 -> 3d3
    &END BETA
  &END BS
&END KIND
```

自旋向下的 `Ni2` 把 α、β 两段对调；O 可以相应地在 2p（`N 2`、`L 1`）上各加电子，使初猜接近 Ni²⁺O²⁻。`&BS` 只是初猜，最终电子结构仍由 SCF 决定，计算后要检查 Mulliken 自旋矩。

---

## 8. DFT+U：NiO 反铁磁

### 8.1 输入文件（2026.2，k 点）

NiO 是典型的 Mott/电荷转移绝缘体，低温下为 AFM-II 序：Ni 的 (111) 面内铁磁排列、相邻 (111) 面反向。最小磁胞是含 2 个 NiO 的菱方胞。示例 [../examples/08-dispersion-spin-dftu/nio_afm_dftu_kp.inp](../examples/08-dispersion-spin-dftu/nio_afm_dftu_kp.inp)：

```text
@SET PROJ  nio_afm_kp
@SET UEFF  6.0

&GLOBAL
  PROJECT ${PROJ}
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    UKS
    MULTIPLICITY 1
    PLUS_U_METHOD MULLIKEN
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 4
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 300
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 10
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &PLUS_U ON
      &END PLUS_U
      &MULLIKEN ON
      &END MULLIKEN
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 5.107186 5.107186 5.107186
      ALPHA_BETA_GAMMA 33.557310 33.557310 33.557310
    &END CELL
    &COORD
      SCALED
      Ni1 0.00 0.00 0.00
      Ni2 0.50 0.50 0.50
      O   0.25 0.25 0.25
      O   0.75 0.75 0.75
    &END COORD
    &KIND Ni1
      ELEMENT Ni
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
      POTENTIAL GTH-PBE-q18
      MAGNETIZATION 2.0
      &DFT_PLUS_U
        L 2
        U_MINUS_J [eV] ${UEFF}
      &END DFT_PLUS_U
    &END KIND
    &KIND Ni2
      ELEMENT Ni
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
      POTENTIAL GTH-PBE-q18
      MAGNETIZATION -2.0
      &DFT_PLUS_U
        L 2
        U_MINUS_J [eV] ${UEFF}
      &END DFT_PLUS_U
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

晶胞由岩盐相晶格常数 a = 4.17 Å（实验值）构造：三个菱方基矢长度为 $a\sqrt{3/2}=5.1072$ Å，夹角 33.557°；Ni1 和 Ni2 位于体对角线上相隔一半的位置，正好分属相邻的 (111) 面。

### 8.2 关键字说明

| 关键字 | 位置 | 说明 |
|---|---|---|
| `PLUS_U_METHOD` | `&DFT` | `MULLIKEN`（默认，基于 net AO 与重叠布居，便宜）、`MULLIKEN_CHARGES`（基于 gross orbital population）、`LOWDIN`（需要对角化重叠矩阵，较贵但可能更稳健） |
| `&DFT_PLUS_U` | `&KIND` | 对该 KIND 开启 +U |
| `L` | `&DFT_PLUS_U` | 施加修正的角动量：2 = d，3 = f |
| `U_MINUS_J` | `&DFT_PLUS_U` | $U_\mathrm{eff}=U-J$，别名 `U_EFF`。**默认单位是 hartree**，必须写 `[eV]` |
| `U_RAMPING`、`EPS_U_RAMPING` | `&DFT_PLUS_U` | 从 0 开始按步长逐渐增大 U，每当 SCF 收敛到 `EPS_U_RAMPING` 以下就增加一次，直到 `U_MINUS_J` |
| `INIT_U_RAMPING_EACH_SCF` | `&DFT_PLUS_U` | 每次波函数优化（如每个 GEO_OPT 步）都从 U = 0 重新开始斜坡；默认只在第一次 |
| `&ENFORCE_OCCUPATION` | `&DFT_PLUS_U` | 强制指定轨道占据（高级用法） |
| `U`、`J`、`N`、`OCCUPATION`、`ALPHA`、`BETA` | `&DFT_PLUS_U` | 手册注明只在使用 SIRIUS 平面波（pwdft）时生效，Quickstep 计算中忽略 |
| `&PRINT/&PLUS_U` | `&DFT` | 打印 +U 壳层的占据矩阵（默认 HIGH 级别） |

> [!CAUTION]
> `U_MINUS_J 6.0` 不写单位时表示 6.0 **hartree**（约 163 eV），结果会完全失去物理意义而且不一定报错。请始终写成 `U_MINUS_J [eV] 6.0`。

### 8.3 版本差异：k 点 DFT+U

| 版本 | DFT+U 与 k 点 | 建议 |
|---|---|---|
| ≤ 2026.1 | 不支持（release notes 将其列为 2026.2 新功能） | 用 Γ 点超胞，见下面的 `nio_afm_dftu_gamma.inp` |
| 2026.2 | 支持 Mulliken 类方法（PR #4855）；2026.2 源码中 Löwdin + k 点计算力时会中止（`Lowdin forces with k-points NYA in DFT+U`） | 用 `PLUS_U_METHOD MULLIKEN` |

另外，2025.2 与 2026.2 源码都有 `The stress tensor is not implemented for the DFT+U method LOWDIN`：需要应力的 CELL_OPT / NPT 不能用 `LOWDIN`。

旧版本可以用示例 [../examples/08-dispersion-spin-dftu/nio_afm_dftu_gamma.inp](../examples/08-dispersion-spin-dftu/nio_afm_dftu_gamma.inp)：2×2×2 常规岩盐超胞（8.34 Å，32 Ni + 32 O），Γ 点 + OT，同时演示 `U_RAMPING`：

```text
&KIND Ni1
  ELEMENT Ni
  BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
  POTENTIAL GTH-PBE-q18
  MAGNETIZATION 2.0
  &DFT_PLUS_U
    L 2
    U_MINUS_J [eV] ${UEFF}
    U_RAMPING [eV] 0.5
    EPS_U_RAMPING 1.0E-4
  &END DFT_PLUS_U
&END KIND
```

### 8.4 运行与检查

```bash
mpirun -np 32 cp2k.psmp -i nio_afm_dftu_kp.inp -o nio_afm_dftu_kp.out
grep "DFT+U energy:" nio_afm_dftu_kp.out | tail -1
grep -A8 "Mulliken Population Analysis" nio_afm_dftu_kp.out | tail -8
grep "Total energy:" nio_afm_dftu_kp.out | tail -1
```

Mulliken 表的列为 `#  Atom  Element  Kind  Atomic population (alpha,beta) Net charge  Spin moment`。判断标准：

- Ni1 与 Ni2 的 `Spin moment` 大小相近、符号相反，O 接近 0，总自旋接近 0——这说明得到了 AFM 解；如果两个 Ni 同号，说明落入了铁磁解，检查 KIND 标签与 MAGNETIZATION；
- Ni 自旋矩的实验值约 1.6–1.9 μB（文献范围），Mulliken 数值依赖基组，作半定量对比；
- 带隙：PBE 对 NiO 给出很小的带隙，加 U 后显著打开（实验光学带隙约 4 eV 左右）。带隙要从 DOS 中读取，见 [第10章](10-dos-pdos-band.md)。

### 8.5 U 值怎么取

1. **文献/数据库值只能作起点**：例如 Materials Project 的 GGA+U 设置中 Ni 取 U = 6.2 eV，但那是基于 VASP PAW 投影拟合的，与 CP2K 的 Mulliken 投影不等价；
2. **按物性标定**：对一组 U 值（如 3、4、5、6、7 eV）计算带隙、磁矩、晶格常数或氧化反应能，选择与实验或高精度方法最接近的值，并在文中注明；
3. **线性响应**：理论上可以按 Cococcioni–de Gironcoli 的线性响应方法计算 U，但需要在 CP2K 中自行实现施加微扰势的流程，超出本章范围；
4. 无论选哪种，比较不同体系的能量时必须使用**同一个 U** 和同一投影方式。

---

## 9. 组合使用的建议

磁性氧化物表面吸附、层状磁性材料等体系经常需要同时开启色散、自旋和 +U。推荐的工作顺序：

1. 先用 PBE（无 U、无色散）把磁性初猜调通，确认得到目标磁序；
2. 加 U（必要时用 `U_RAMPING`），从上一步的波函数重启（`SCF_GUESS RESTART`）；
3. 加色散校正（对 SCF 几乎没有影响，可以直接加）；
4. 体相 CELL_OPT（Mulliken 类 +U 支持应力）→ 建 slab → GEO_OPT；
5. 每次改变 U、磁序或色散方案，都重新检查自旋矩和能量。

> [!NOTE]
> 色散校正对能量和力是“加性”的，不改变电子结构；+U 和自旋设置会改变电子结构。所以 SCF 难收敛时问题几乎总在后两者，而不在色散。

---

## 参数选择与收敛建议

| 项目 | 推荐设置 | 说明 |
|---|---|---|
| 色散方案（GGA） | PBE-D3(BJ) 或 PBE-D4 | 显式写 `TYPE` 和 `REFERENCE_FUNCTIONAL` |
| 色散方案（杂化） | PBE0-D3(BJ)、B3LYP-D3(BJ) | `REFERENCE_FUNCTIONAL` 改成对应泛函，见[第9章](09-hybrid-admm.md) |
| `R_CUTOFF` | 默认约 10.6 Å；周期体系可增大到 15–20 Å 测试 | 能量变化 < 0.1 meV/原子即可 |
| 相互作用能基组 | TZVP 以上，或做 BSSE 校正 | DZVP 的 BSSE 可能与色散能量同量级 |
| 非局域 vdW `CUTOFF` | 从 100 Ry 起测试 | 独立于 `&MGRID` |
| `MULTIPLICITY` | 逐个尝试可能的值 | 比较能量后再定基态 |
| 磁性初猜 | `MAGNETIZATION` ≈ 预期磁矩（μB） | 反铁磁必须分 KIND |
| `U_MINUS_J` | 必须写 `[eV]`；按物性标定 | 常见 3d 氧化物在 3–7 eV 范围内试算 |
| `PLUS_U_METHOD` | `MULLIKEN`（默认） | k 点和应力都支持；`LOWDIN` 不支持应力 |
| SCF（磁性/+U） | 对角化：Broyden，`ALPHA` 0.05–0.2；OT：`FULL_SINGLE_INVERSE` + `OUTER_SCF` | `MAX_SCF` 适当放大 |

---

## 常见问题与报错

**1. `Use the LSD option for an odd number of electrons`**
原因：奇数电子却没有开自旋极化。
解决：`&DFT` 中加 `UKS` 并设 `MULTIPLICITY 2`（或正确值）。若体系本应闭壳层，检查结构中的 H 原子数和 `CHARGE`。

**2. `LSD: too few electrons for this multiplicity`**
原因：`MULTIPLICITY` 太大（例如 `CHARGE` 写错导致电子数变少）。
解决：核对总电子数（输出中的 `Number of electrons:`）和 `CHARGE`。

**3. 反铁磁体系收敛成铁磁或非磁解**
原因：只用了一个 KIND，或 `MAGNETIZATION` 符号没有区分，或 `MULTIPLICITY` 与目标磁序不符。
解决：按 7.2 节拆分 `Ni1`/`Ni2` 等 KIND，给相反的 `MAGNETIZATION`，`MULTIPLICITY 1`；必要时用 `&BS` 构造离子性初猜，或先用较小的 `ALPHA` 混合。

**4. DFT-D3 参数文件找不到 / 泛函识别失败**
原因：`dftd3.dat` 不在 `CP2K_DATA_DIR` 下；或使用了自定义组合泛函（如手写 `&PBE SCALE_X 0.75` 的 PBE0），CP2K 无法推断参考泛函。
解决：设置 `CP2K_DATA_DIR` 或给 `PARAMETER_FILE_NAME` 写完整路径；显式写 `REFERENCE_FUNCTIONAL`（如 `PBE0`），或用 `D3BJ_SCALING` 直接给参数。

**5. 同一个旧输入在 2025.2 与 2026.2 下色散能量不同**
原因：`&PAIR_POTENTIAL/TYPE` 的默认值从 `DFTD3`（零阻尼）变成了 `DFTD3(BJ)`。
解决：显式写 `TYPE`。比较新旧结果之前先核对输出中的色散类型。

**6. `TYPE DFTD4` 报错或 D4 能量为零**
原因：CP2K 编译时没有链接 libdftd4。
解决：`cp2k.psmp -v` 检查 cp2kflags 是否含 `libdftd4`；没有就改用 `DFTD3(BJ)` 或换带 D4 的构建。

**7. DFT+U 结果离谱（带隙巨大、能量异常、SCF 发散）**
原因：`U_MINUS_J` 没有写 `[eV]`，被当作 hartree。
解决：写成 `U_MINUS_J [eV] 6.0`，并在输出中检查 DFT+U 设置的打印值。

**8. `The stress tensor is not implemented for the DFT+U method LOWDIN`**
原因：CELL_OPT、NPT 或打开了 `STRESS_TENSOR` 时使用了 `PLUS_U_METHOD LOWDIN`。
解决：改用 `MULLIKEN`。

**9. `Lowdin forces with k-points NYA in DFT+U`（2026.2）**
原因：k 点 + Löwdin DFT+U 计算力（GEO_OPT、ENERGY_FORCE）。
解决：k 点计算使用 `PLUS_U_METHOD MULLIKEN`。

**10. +U 体系 SCF 振荡、多个解之间跳动**
原因：DFT+U 势能面上有多个局部极小（不同的轨道占据）。
解决：先不加 U 收敛，再 `U_RAMPING` 逐步加；减小混合 `ALPHA`；用 `&BS` 指定占据；对不同初猜分别计算并比较能量，选最低者。

---

## 练习

1. 用苯二聚体示例分别计算 PBE、PBE-D3（零阻尼，`TYPE DFTD3`）、PBE-D3(BJ) 的相互作用能，再用 `RUN_TYPE BSSE` 得到 BSSE 校正后的值，整理成表格。
2. 把苯二聚体的层间距（坐标文件中第二个苯的 z 坐标）从 3.2 Å 扫到 4.4 Å（步长 0.2 Å，用 Python 改写 xyz），画出 PBE 与 PBE-D3(BJ) 的势能曲线，比较极小位置。
3. 对 O2 分别用 `MULT 3` 和 `MULT 1` 优化，比较键长和能量；再用 `ROKS` 计算三重态，与 UKS 结果比较 $\langle S^2\rangle$ 与能量。
4. 对 NiO 的 4 原子菱方胞，在 2026.2 下取 `UEFF` = 0、2、4、6、8 eV，记录 Ni 自旋矩和 `DFT+U energy`；结合[第10章](10-dos-pdos-band.md)的 DOS 读出带隙，画出带隙随 U 的变化。
5. 把 NiO 输入中两个 Ni 都改成 `MAGNETIZATION 2.0` 并设 `MULTIPLICITY 5`（铁磁，2 个 Ni 各 2 μB），与 AFM 能量比较，估算交换作用强弱。
6. 去掉石墨示例中的色散校正重跑 CELL_OPT，再分别用 D3(BJ) 和 D4 计算，比较三种情况下的 c 和层间结合能。

---

## 参考资料

**官方文档**

- `&XC/&VDW_POTENTIAL/&PAIR_POTENTIAL` 输入参考（2026.2）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/XC/VDW_POTENTIAL/PAIR_POTENTIAL.html>
- `&NON_LOCAL`（vdW-DF、rVV10）：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/XC/VDW_POTENTIAL/NON_LOCAL.html>
- `&DFT`（`UKS`、`MULTIPLICITY`、`PLUS_U_METHOD`）：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT.html>
- `&KIND/&DFT_PLUS_U`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND/DFT_PLUS_U.html>
- `&KIND/&BS`：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND/BS.html>
- `&FORCE_EVAL/&BSSE`：<https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/BSSE.html>
- How to make a SCF run converge（含磁性体系建议）：<https://manual.cp2k.org/trunk/methods/dft/convergence.html>
- CP2K v2026.2 release notes（k 点 DFT+U 等）：<https://github.com/cp2k/cp2k/releases/tag/v2026.2>
- 预设 vdW-DF 泛函段 `data/xc_section/`：<https://github.com/cp2k/cp2k/tree/master/data/xc_section>
- DFT-D4 参数表（`REFERENCE_FUNCTIONAL` 可用名称）：<https://github.com/dftd4/dftd4/blob/main/src/dftd4/param.f90>

**相关 GitHub 仓库（只链接，不复制）**

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/QS/regtest-dft-vdw-corr-*`、`regtest-dft-d4-auto-ref`、`regtest-plus_u`、`regtest-plus_u-kp`、`regtest-bs` 等 regtest 展示了色散、DFT+U 和 `&BS` 的标准写法
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含 CP2K DFT+U 与磁性体系设置经验
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：中文分级教程（2025.2 验证）
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库
- [JunboLu/CP2K_kit](https://github.com/JunboLu/CP2K_kit)：国内开发的 CP2K 工具集，可批量生成输入
- [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts)：含 U 值收敛扫描的脚本思路（无 license，仅供参考）
- [szt2006/cp2k-aimd-skill](https://github.com/szt2006/cp2k-aimd-skill)：中文“症状→处方”对照表，可作为排错选题参考（license 不明确，只链接）

**本仓库脚本**：`../scripts/cp2k_outparse.py`（提取总能量、Mulliken 电荷/自旋、警告汇总）、`../scripts/plot_opt.py`（CELL_OPT 收敛图）。用法请运行 `python <脚本> -h` 查看。

**关键文献**

- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* **132**, 154104 (2010)，DFT-D3：<https://doi.org/10.1063/1.3382344>
- S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* **32**, 1456 (2011)，BJ 阻尼：<https://doi.org/10.1002/jcc.21759>
- E. Caldeweyher *et al.*, *Phys. Chem. Chem. Phys.* **22**, 8499 (2020)，周期体系 D4：<https://doi.org/10.1039/d0cp00502a>
- M. Dion *et al.*, *Phys. Rev. Lett.* **92**, 246401 (2004)，vdW-DF：<https://doi.org/10.1103/PhysRevLett.92.246401>
- R. Sabatini, T. Gorni, S. de Gironcoli, *Phys. Rev. B* **87**, 041108 (2013)，rVV10：<https://doi.org/10.1103/PhysRevB.87.041108>
- F. Tran, J. Hutter, *J. Chem. Phys.* **138**, 204103 (2013)，CP2K 中的非局域 vdW 泛函：<https://doi.org/10.1063/1.4807332>
- P. Jurečka, J. Šponer, J. Černý, P. Hobza, *Phys. Chem. Chem. Phys.* **8**, 1985 (2006)，S22 基准：<https://doi.org/10.1039/B600027D>
- S. F. Boys, F. Bernardi, *Mol. Phys.* **19**, 553 (1970)，counterpoise 校正：<https://doi.org/10.1080/00268977000101561>
- S. L. Dudarev *et al.*, *Phys. Rev. B* **57**, 1505 (1998)，DFT+U（Dudarev 形式，以 NiO 为例）：<https://doi.org/10.1103/PhysRevB.57.1505>
- M. Cococcioni, S. de Gironcoli, *Phys. Rev. B* **71**, 035105 (2005)，线性响应求 U：<https://doi.org/10.1103/PhysRevB.71.035105>

---

[← 上一章](07-geo-cell-opt.md) | [返回目录](../README.md) | [下一章 →](09-hybrid-admm.md)
