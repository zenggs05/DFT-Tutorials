# 第2章 输入文件结构与第一个计算（CP2K Input File Structure & First Calculation）

本章系统讲解 CP2K 输入文件（`.inp`）的语法与层级结构：`&GLOBAL`、`&FORCE_EVAL`（`&DFT`、`&SUBSYS`）、`&MOTION`、`&EXT_RESTART` 四大顶层段，`&CELL`/`&COORD`/`&TOPOLOGY`/`&KIND` 的写法，方括号单位，以及 `@SET`、`${VAR}`、`@INCLUDE`、`@IF` 预处理器和 `-E/--set` 命令行变量。随后完成第一批 CP2K 单点计算（single-point energy & forces）：孤立 H2O 分子、周期性 Si 晶体和金红石 TiO2，逐段解读输出文件（SCF 迭代、总能量、原子力、Mulliken 电荷、timing 表、重启文件），并说明为什么小原胞只用 Gamma 点是不够的。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 看懂任何一个 CP2K 输入文件的骨架：哪一段管“算什么”、哪一段管“体系是什么”、哪一段管“原子怎么动”；
- 正确书写段（section）、关键字（keyword）、段参数、默认关键字（如 `&COORD` 中的原子行）和单位 `[angstrom]`；
- 用预处理器写“模板化”的输入，用 `cp2k.psmp -E VAR=value` 从命令行切换参数，并理解 `@SET` 与 `-E` 的优先级；
- 用 `cp2k.psmp -c`（语法检查）和 `-e`（展开并显式写出全部默认值）在提交作业前排错；
- 完成 `RUN_TYPE ENERGY` 与 `ENERGY_FORCE` 计算，用 `grep` 从输出中取出总能量、SCF 收敛信息、原子力、Mulliken 电荷与耗时；
- 知道各种 `RUN_TYPE` 分别在哪一章展开。

本章示例：

| 示例目录 | 内容 |
|---|---|
| [`examples/02-first-calc-h2o-si/`](../examples/02-first-calc-h2o-si/README.md) | `h2o_energy.inp`、`h2o_energy_force.inp`、`h2o.xyz`、`si_bulk8_energy.inp`、`preprocessor_demo.inp` + `h2o_kinds.inc` |
| [`examples/02-energy-force-rutile/`](../examples/02-energy-force-rutile/README.md) | 金红石 TiO2：原有 Gamma 点输入 `input_file.inp`，改进版 `rutile_kpoints.inp`、`rutile_supercell_gamma.inp`，Slurm/LSF 作业脚本 |

> [!NOTE]
> 本章假设你已按 [第1章](01-install.md) 装好 CP2K 并设置了 `CP2K_DATA_DIR`。示例输入在 CP2K 2023.x–2026.2 中均可运行；输出格式有差别的地方会注明版本。

---

## 2. 输入文件的基本语法

### 2.1 段与关键字

CP2K 输入是**层级结构**，由“段”（SECTION）和“关键字”（KEYWORD）组成：

```text
&SECTION_NAME [段参数]        ! 段以 & 开头
  KEYWORD_NAME  value          ! 关键字 + 值
  &SUBSECTION
    KEYWORD  value1 value2
  &END SUBSECTION
&END SECTION_NAME              ! 段以 &END 结束（建议写上段名，便于检查配对）
```

| 规则 | 说明 |
|---|---|
| 大小写 | 段名、关键字名、枚举值（如 `PBE`、`MONKHORST-PACK`）不区分大小写；**文件名、预处理器变量名区分大小写** |
| 缩进 | 不要求，习惯上每层缩进 2 个空格 |
| 注释 | `#` 或 `!` 之后到行尾为注释，可以单独成行，也可以写在行尾 |
| 续行 | 行尾写 `\` 表示下一行接着写（长列表时有用） |
| 段的顺序 | 顶层段之间、同级段之间的顺序任意（习惯把 `&GLOBAL` 写在最前） |
| 可重复的段 | 如 `&KIND`（每个元素/原子类型一个）、`&FORCE_EVAL`（多力场时） |
| 默认值 | 显式写出某个段就启用该功能；段内没写的关键字取默认值 |

**两种特殊写法**：

1. **段参数（SECTION_PARAMETERS）**：值直接写在 `&SECTION` 后面。例如 `&FORCES ON`、`&OT ON`、`&SMEAR ON`、`&KIND O`（`O` 是这个 KIND 的名字）、`&XC_FUNCTIONAL PBE`。
2. **默认关键字（DEFAULT_KEYWORD）**：不写关键字名，直接写数据行。最典型的是 `&COORD` 里的 `O 0.0 0.0 0.0`。

### 2.2 值的类型

| 类型 | 写法 |
|---|---|
| 逻辑 | 真：`T`、`.TRUE.`、`TRUE`、`Y`、`YES`、`ON`、`1`；假：`F`、`.FALSE.`、`FALSE`、`N`、`NO`、`OFF`、`0`。**只写关键字名不写值，表示真** |
| 实数 | `0.001`、`1E-3`、`1.0E-6` 均可；还可以写简单表达式，如 `2/3`、`-SQRT(3)/2`（支持大写的 `COS`、`EXP`、`LOG`、`LOG10`、`SIN`、`SQRT`、`TAN`） |
| 整数、字符串、枚举 | 按手册给出的类型；枚举值必须是允许列表中的一个 |
| 带单位 | `[单位]` 写在关键字名和数值之间，如 `ABC [angstrom] 10 10 10`（见 §5） |

> [!TIP]
> 不确定某个关键字的类型、默认值或允许的取值时，查 2026.2 输入参考 <https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT.html>（trunk 版：<https://manual.cp2k.org/trunk/CP2K_INPUT.html>），或者用 §7 的 `cp2k.psmp -e` 让程序把默认值全部写出来。

---

## 3. 顶层段一览

一个典型输入的骨架如下（`...` 表示省略）：

```text
&GLOBAL                 ! 全局：项目名、任务类型、输出详细程度
  PROJECT  h2o
  RUN_TYPE ENERGY_FORCE
&END GLOBAL

&FORCE_EVAL             ! 力的计算引擎：用什么方法、算什么体系
  METHOD QS
  &DFT                  ! 电子结构：基组/赝势文件、网格、SCF、泛函、泊松求解器 ...
    ...
  &END DFT
  &SUBSYS               ! 体系：晶胞、坐标、元素类型（基组、赝势）
    &CELL ... &END CELL
    &COORD ... &END COORD
    &KIND O ... &END KIND
  &END SUBSYS
  &PRINT ... &END PRINT ! FORCE_EVAL 级别的输出（原子力、应力）
&END FORCE_EVAL

&MOTION                 ! 原子怎么动：GEO_OPT / CELL_OPT / MD / BAND / 约束 / 轨迹输出
  ...
&END MOTION

&EXT_RESTART            ! 从 *-1.restart 文件续算
  RESTART_FILE_NAME h2o-1.restart
&END EXT_RESTART
```

此外还有一些专用顶层段：`&VIBRATIONAL_ANALYSIS`（振动分析，[第14章](14-vibration-thermo.md)）、`&MULTIPLE_FORCE_EVALS`（多个力场组合）、`&FARMING`（批量运行多个输入）、`&ATOM`（原子程序，生成/测试赝势和基组）等。

### 3.1 &GLOBAL

| 关键字 | 默认值 | 说明 |
|---|---|---|
| `PROJECT`（= `PROJECT_NAME`） | `PROJECT` | 所有输出文件的前缀，如 `h2o-RESTART.wfn`、`h2o-pos-1.xyz`。含空格时会被替换成下划线并给出警告 |
| `RUN_TYPE` | **`ENERGY_FORCE`** | 任务类型（§13）。不写时默认算能量和力 |
| `PRINT_LEVEL`（= `IOLEVEL`） | `MEDIUM` | `SILENT`/`LOW`/`MEDIUM`/`HIGH`/`DEBUG`；大体系 MD 常用 `LOW` 减少输出 |
| `WALLTIME` | 无 | 最长运行时间（秒或 `HH:MM:SS`）；设得比作业系统时限略短，让 CP2K 在被强制杀掉之前自行结束（优化/MD 可从最近的重启文件续算） |
| `ECHO_INPUT` | `F` | 在输出开头回显（默认值已显式化的）输入 |
| `PREFERRED_DIAG_LIBRARY` | 自动 | `ELPA`/`SCALAPACK` 等，排错和性能调优时用（第18章） |

### 3.2 &FORCE_EVAL

| 内容 | 说明 |
|---|---|
| `METHOD` | `QS`（= `QUICKSTEP`，DFT 等电子结构）、`FIST`（经典力场）、`QMMM`、`SIRIUS`、`NNP`、`MIXED` 等 |
| `STRESS_TENSOR` | `NONE`（默认）/`ANALYTICAL`/`NUMERICAL`；晶胞优化和 NPT 需要（[第7章](07-geo-cell-opt.md)） |
| `&DFT` | 电子结构设置：`BASIS_SET_FILE_NAME`、`POTENTIAL_FILE_NAME`、`CHARGE`、`MULTIPLICITY`、`UKS`，以及子段 `&MGRID`、`&QS`、`&SCF`、`&XC`、`&POISSON`、`&KPOINTS`、`&PRINT` |
| `&SUBSYS` | 体系定义（§4） |
| `&PRINT` | 如 `&FORCES ON`（打印原子力）、`&STRESS_TENSOR ON` |

`&DFT` 下最常用的子段：

| 子段 | 作用 | 详见 |
|---|---|---|
| `&MGRID` | 多重网格截断 `CUTOFF`、`REL_CUTOFF`、`NGRIDS` | [第3章](03-gpw-basis-pseudo.md)、[第4章](04-cutoff-convergence.md) |
| `&QS` | Quickstep 总体设置：`METHOD GPW/GAPW/xTB...`、`EPS_DEFAULT`、`EXTRAPOLATION` | [第3章](03-gpw-basis-pseudo.md) |
| `&SCF` | 自洽：`EPS_SCF`、`MAX_SCF`、`SCF_GUESS`，`&OT` 或 `&DIAGONALIZATION` + `&MIXING` + `&SMEAR` | [第5章](05-scf-strategy.md) |
| `&XC` | 交换关联泛函、色散校正 `&VDW_POTENTIAL`、杂化泛函 `&HF` | [第8章](08-dispersion-spin-dftu.md)、[第9章](09-hybrid-admm.md) |
| `&POISSON` | 静电边界条件：`PERIODIC`、`POISSON_SOLVER` | [第6章](06-periodic-surface-kpoints.md) |
| `&KPOINTS` | k 点采样；不写即 Gamma 点 | [第6章](06-periodic-surface-kpoints.md) |
| `&PRINT` | 电子结构输出：Mulliken/Hirshfeld 电荷、MO、cube、DOS | [第10章](10-dos-pdos-band.md)、[第11章](11-cube-charge-workfunction.md) |

### 3.3 &MOTION

`&MOTION` 只在原子需要“动”的任务中起作用（`GEO_OPT`、`CELL_OPT`、`MD`、`BAND` 等）。常用子段：`&GEO_OPT`、`&CELL_OPT`、`&MD`、`&BAND`、`&CONSTRAINT`（固定原子、约束）、`&FREE_ENERGY`（元动力学）、`&PRINT`（轨迹 `TRAJECTORY`、速度、力、重启文件 `RESTART`）。单点计算不需要写 `&MOTION`。

### 3.4 &EXT_RESTART

用于从上一次运行写出的 `<PROJECT>-1.restart` 续算：

```text
&EXT_RESTART
  RESTART_FILE_NAME h2o-1.restart
&END EXT_RESTART
```

`.restart` 文件本身就是一个完整的 CP2K 输入（包含最新坐标、晶胞、速度、计数器），也可以直接把它当输入运行。续算策略在 [第7章](07-geo-cell-opt.md)、第12章 和 第18章 展开。

---

## 4. &SUBSYS 详解：晶胞、坐标与元素类型

### 4.1 &CELL

| 关键字 | 默认单位 | 说明 |
|---|---|---|
| `ABC` | Å | 三个晶格矢量的长度 |
| `ALPHA_BETA_GAMMA` | deg | 与 `ABC` 配合的三个夹角，默认 `90 90 90`；六方晶系写 `ALPHA_BETA_GAMMA 90 90 120` |
| `A`、`B`、`C` | Å | 直接给出三个晶格矢量的笛卡尔分量（写了 `ABC` 或 `CELL_FILE_NAME` 时被忽略） |
| `PERIODIC` | — | `XYZ`（默认）、`XY`、`XZ`、`YZ`、`X`、`Y`、`Z`、`NONE`；必须与 `&POISSON PERIODIC` 一致 |
| `CELL_FILE_NAME` + `CELL_FILE_FORMAT` | — | 从文件读晶胞，格式 `CP2K`/`CIF`/`XSC`/`EXTXYZ`/`XYZ`/`PDB`（2026.2 输入参考列出的取值） |
| `MULTIPLE_UNIT_CELL` | — | 把晶胞重复 n1×n2×n3 次；**必须同时在 `&TOPOLOGY` 中写同样的 `MULTIPLE_UNIT_CELL`，坐标才会被复制** |
| `SYMMETRY` | — | 晶胞对称性约束（晶胞优化时使用，[第7章](07-geo-cell-opt.md)） |

三种等价写法示例（立方 Si，a = 5.431 Å）：

```text
&CELL
  ABC [angstrom] 5.431 5.431 5.431
&END CELL

&CELL
  A [angstrom] 5.431 0.000 0.000
  B [angstrom] 0.000 5.431 0.000
  C [angstrom] 0.000 0.000 5.431
&END CELL

&CELL
  CELL_FILE_NAME  si.cif
  CELL_FILE_FORMAT CIF
&END CELL
```

> [!IMPORTANT]
> **孤立分子也必须写 `&CELL`**：GPW 需要一个盒子来定义实空间网格。孤立体系写 `PERIODIC NONE`，并在 `&DFT/&POISSON` 中同样写 `PERIODIC NONE` 和非周期泊松求解器（`MT` 或 `WAVELET`）。两处 `PERIODIC` 不一致会导致错误或物理上错误的静电能。

### 4.2 坐标：&COORD 或 &TOPOLOGY

**方式一：直接写在 `&COORD` 中**（小体系、教学示例）：

```text
&COORD
  UNIT angstrom        ! 默认就是 angstrom
  SCALED F             ! T 表示分数坐标（以晶格矢量为单位）
  O   0.000000   0.000000   0.117300
  H   0.000000   0.756990  -0.468540
  H   0.000000  -0.756990  -0.468540
&END COORD
```

每行第一列是 **KIND 名**（通常就是元素符号），后面三列是坐标。

**方式二：从外部文件读取**（大体系、推荐）：

```text
&TOPOLOGY
  COORD_FILE_NAME   h2o.xyz
  COORD_FILE_FORMAT XYZ          ! XYZ / CIF / PDB / CP2K / CRD / G96 / XTL
  &CENTER_COORDINATES            ! 把分子平移到盒子中心（孤立分子常用）
  &END CENTER_COORDINATES
&END TOPOLOGY
```

| 格式 | 注意 |
|---|---|
| `XYZ` | 标准 xyz：第 1 行原子数，第 2 行注释，之后每行 `元素 x y z`（Å）。**xyz 不含晶胞信息**，晶胞要在 `&CELL` 中另写（2026.2 起可用 `CELL_FILE_FORMAT EXTXYZ` 从扩展 xyz 的注释行读晶胞） |
| `CIF` | 同时写 `&CELL CELL_FILE_NAME x.cif` + `CELL_FILE_FORMAT CIF`，晶胞与坐标都从 CIF 读；含无序占位的 CIF 必须先处理 |
| `PDB` | 生物/QM-MM 体系常用，可携带残基信息 |

结构建模（超胞、表面、缺陷）与格式转换见 公共篇第3章 和 公共篇第4章。

### 4.3 &KIND：给每种原子指定基组和赝势

```text
&KIND O
  BASIS_SET DZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q6
&END KIND

&KIND Fe1                ! KIND 名可以不是元素符号（例如区分两种自旋方向的 Fe）
  ELEMENT Fe             ! 这时必须用 ELEMENT 指明元素
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE-q16
  MAGNETIZATION 4        ! 初始磁矩（第8章）
&END KIND
```

| 关键字 | 说明 |
|---|---|
| `BASIS_SET [类型] 名称` | 类型默认 `ORB`（轨道基组）；`AUX_FIT`（ADMM）、`RI_AUX`（RI）等辅助基组也在这里指定（[第3章](03-gpw-basis-pseudo.md)） |
| `POTENTIAL [类型] 名称` | GTH 赝势名（如 `GTH-PBE-q6`）；全电子写 `POTENTIAL ALL` |
| `ELEMENT` | KIND 名不是元素符号时必须写 |
| `MAGNETIZATION` | 自旋极化计算的初始磁矩 |
| `GHOST` | 鬼原子（只有基组没有核和电子），用于 BSSE 校正 |
| `POTENTIAL_FILE_NAME` | 为这个 KIND 单独指定赝势文件（覆盖 `&DFT` 中的设置） |

基组名、赝势名必须在 `&DFT` 指定的文件中真实存在，且价电子数后缀 `-qN` 要一致，这是 [第3章](03-gpw-basis-pseudo.md) 的主题。

---

## 5. 单位写法

带量纲的关键字可以在名字后用方括号写单位，CP2K 在读入时立即换算成内部原子单位：

```text
ABC [angstrom] 10.0 10.0 10.0
ABC [bohr] 18.897 18.897 18.897
CUTOFF [Ry] 400
CUTOFF [eV] 5442              ! 约 400 Ry
TIMESTEP [fs] 0.5
TEMPERATURE [K] 300
U_MINUS_J [eV] 4.0
EXTERNAL_PRESSURE [bar] 1.01325
```

常用单位名：能量 `hartree`、`eV`、`Ry`、`kcalmol`、`kjmol`、`K_e`、`wavenumber_e`；长度 `angstrom`、`bohr`、`nm`、`pm`；时间 `fs`、`ps`；温度 `K`；压强 `bar`、`atm`、`GPa`、`Pa`；角度 `deg`、`rad`。完整列表见 <https://manual.cp2k.org/cp2k-2026_2-branch/units.html>。

**每个关键字的默认单位不同**，不写方括号时按默认单位理解。最容易出错的几个（取自 2026.2 源码中的关键字定义）：

| 关键字 | 默认单位 | 默认值 | 易错点 |
|---|---|---|---|
| `&CELL ABC`、`&COORD` 坐标 | Å | — | 与 VASP/QE 习惯一致 |
| `&MGRID CUTOFF` | Ry | 280 | 不是 eV，也不是 Hartree |
| `&MGRID REL_CUTOFF` | Ry | 20 | — |
| `&MD TIMESTEP` | fs | 0.5 | — |
| `&MD TEMPERATURE` | K | 300 | — |
| `&SMEAR ELECTRONIC_TEMPERATURE` | K | 300 | — |
| `&GEO_OPT MAX_FORCE` | hartree/bohr | 4.5E-4 | 换算成 eV/Å 约为 0.023 |
| `&GEO_OPT MAX_DR` | bohr | 3.0E-3 | — |
| `&CELL_OPT EXTERNAL_PRESSURE` | bar | **100** | 默认不是 1 atm，做常压优化要显式写 `1.01325` |
| `&DFT_PLUS_U U_MINUS_J` | **hartree** | 0 | 文献中的 U 多以 eV 给出，务必写 `[eV]` |
| `ALPHA_BETA_GAMMA` | deg | 90 90 90 | — |

> [!IMPORTANT]
> **输出文件中的数值，除非特别标注，一律是原子单位**：能量 Hartree，长度 bohr，力 Hartree/bohr。常用换算：1 Ha = 27.211386 eV，1 bohr = 0.529177 Å，1 Ha/bohr = 51.42208 eV/Å。在 `&GLOBAL/&PRINT` 中设 `PHYSCON T` 可让 CP2K 打印它使用的全部物理常数与换算因子。

---

## 6. 预处理器：@SET、${VAR}、@INCLUDE、@IF

CP2K 在解析输入前会先运行一个简单的预处理器。它让你把一个输入写成“模板”，用变量控制截断、基组、k 点等参数，是收敛测试和批量计算的基础（[第4章](04-cutoff-convergence.md)）。

### 6.1 指令一览

| 指令 | 作用 | 示例 |
|---|---|---|
| `@SET 变量 值` | 定义或重新赋值变量（行内其余部分都是值） | `@SET CUT 400` |
| `${变量}` 或 `$变量` | 展开变量；推荐带花括号的写法 | `CUTOFF ${CUT}` |
| `${变量-默认值}` | 变量未定义时使用默认值（变量名不能含 `-`） | `CUTOFF ${CUT-400}` |
| `@INCLUDE 文件` | 把文件内容原样插入此处；文件名可加单/双引号 | `@INCLUDE 'kinds.inc'` |
| `@IF 表达式` … `@ENDIF` | 条件为假时跳过中间各行 | `@IF ${USE_OT} == 1` |
| `@PRINT 文本` | 解析时把文本（变量已展开）打印到输出，用于调试 | `@PRINT cutoff=${CUT}` |
| `@XCTYPE 名称` | 插入数据目录 `xc_section/名称.sec` 中预定义的 `&XC` 段 | `@XCTYPE PBE` |

`@IF` 的判断规则（2026.2 源码 `cp_parser_inpp_methods.F`）：

- `@IF a == b` 与 `@IF a /= b`：**按字符串比较**（不是数值比较），`4.0` 与 `4` 不相等；
- `@IF 值`：值为空或为 `0` 时为假，其余为真；
- 没有 `@ELSE`；`@IF` 不能嵌套（跳过时遇到的第一个 `@ENDIF` 就结束跳过）。

### 6.2 命令行变量 -E / --set 与优先级

```bash
cp2k.psmp -E CUT=500 -E BASIS=TZV2P -i template.inp -o run.out
cp2k.psmp --set CUT=500 -i template.inp -o run.out      # 等价写法
```

`-E` 只是给变量一个**初值**。解析时文件里的 `@SET` 会**覆盖**这个初值，所以：

| 输入中的写法 | `-E VAR=x` 是否生效 |
|---|---|
| `@SET VAR 400` 然后用 `${VAR}` | **不生效**（被 `@SET` 覆盖） |
| 只用 `${VAR}`，文件中没有 `@SET VAR` | 生效；但不加 `-E` 时报 `Variable ${VAR} not defined` |
| 只用 `${VAR-400}` | 生效；不加 `-E` 时取 400（**推荐写法**） |

### 6.3 示例：preprocessor_demo.inp

示例 [`preprocessor_demo.inp`](../examples/02-first-calc-h2o-si/preprocessor_demo.inp) 把上面的指令都用了一遍（节选）：

```text
@SET SYSTEM h2o
@SET XC PBE
@SET BOX 10.0

@PRINT system=${SYSTEM} xc=${XC} basis=${BASIS-DZVP} cutoff=${CUTOFF-350} ot=${USE_OT-1}

&GLOBAL
  PROJECT ${SYSTEM}_${BASIS-DZVP}_c${CUTOFF-350}_ot${USE_OT-1}
  RUN_TYPE ENERGY
  PRINT_LEVEL ${PL-MEDIUM}
&END GLOBAL
...
    &MGRID
      CUTOFF ${CUTOFF-350}
      REL_CUTOFF 50
    &END MGRID
...
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 50
      @IF ${USE_OT-1} == 1
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      @ENDIF
    &END SCF
...
    &TOPOLOGY
      COORD_FILE_NAME ${SYSTEM}.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    @INCLUDE 'h2o_kinds.inc'
```

被包含的 [`h2o_kinds.inc`](../examples/02-first-calc-h2o-si/h2o_kinds.inc) 可以直接使用主文件中的变量：

```text
&KIND O
  BASIS_SET ${BASIS-DZVP}-MOLOPT-GTH
  POTENTIAL GTH-${XC}-q6
&END KIND
&KIND H
  BASIS_SET ${BASIS-DZVP}-MOLOPT-GTH
  POTENTIAL GTH-${XC}-q1
&END KIND
```

试运行：

```bash
cd CP2K/examples/02-first-calc-h2o-si
cp2k.ssmp -e preprocessor_demo.inp | less                       # 看展开后的完整输入
cp2k.ssmp -E BASIS=TZV2P -E CUTOFF=500 -i preprocessor_demo.inp -o demo_tzv2p.out
grep "GLOBAL| Project name" demo_tzv2p.out                       # h2o_TZV2P_c500_ot1
cp2k.ssmp -E USE_OT=0 -i preprocessor_demo.inp -o demo_diag.out  # 不写 &OT -> 默认对角化
cp2k.ssmp -E XC=BLYP -e preprocessor_demo.inp | grep XC_FUNCTIONAL   # 仍是 PBE：被 @SET 覆盖
```

> [!WARNING]
> 预处理器的几个陷阱：
> - **`@SET` 行末不要写注释**：该行变量名之后的全部内容都会成为变量的值。
> - 预处理器会处理**非注释行中出现的任何 `$` 和 `@`**，即使它们出现在行尾注释里。行尾注释中不要写 `$`、`@`（整行注释不受影响）。
> - `@INCLUDE` 的文件名先按当前工作目录（或绝对路径）查找，再到数据目录中找；在别的目录运行时要注意相对路径。
> - 变量名只能由字母、数字和下划线组成，以字母或下划线开头，**区分大小写**。

---

## 7. 提交前检查输入：-c、-e 与 cp2klint

| 命令 | 作用 | 典型用途 |
|---|---|---|
| `cp2k.psmp -c x.inp` | 语法检查（段/关键字是否存在、类型是否正确、预处理器是否正确） | 每次提交前都跑一下，几秒钟 |
| `cp2k.psmp -e x.inp` | 输出“预处理后 + 所有默认值显式写出”的完整输入 | 查默认值；确认 `-E` 变量是否生效 |
| `cp2k.psmp --xml` | 导出当前版本的输入结构 `cp2k_input.xml` | 给第三方工具用 |
| `cp2klint x.inp` | cp2k-input-tools 提供的 linter（`pip install cp2k-input-tools`） | 在没有 CP2K 的笔记本上检查 |

```bash
cp2k.psmp -c h2o_energy.inp
# SUCCESS, the input could be parsed correctly.
#          This does not guarantee that this input is meaningful
#          or will run successfully
```

`-c` 只检查**语法**，不会检查基组名是否存在、结构是否合理、PERIODIC 是否一致；这些错误要到真正运行时才暴露。

> [!NOTE]
> `cp2klint` 依赖 cp2k-input-tools 自带的 CP2K 输入结构描述（XML），其版本可能落后于 CP2K 2026.2，遇到新关键字（如 2026.2 的 `&DOS`）可能误报；以你所用版本的 `cp2k.psmp -c` 为准。cp2k-input-tools 的其他功能见 第17章。

---

## 8. 实战一：H2O 单点能量（RUN_TYPE ENERGY）

### 8.1 输入文件

[`h2o_energy.inp`](../examples/02-first-calc-h2o-si/h2o_energy.inp)（坐标在 [`h2o.xyz`](../examples/02-first-calc-h2o-si/h2o.xyz) 中，水分子取实验平衡构型 r(O–H) = 0.9572 Å、∠HOH = 104.52°）：

```text
&GLOBAL
  PROJECT h2o_energy          ! prefix of all files written by this run
  RUN_TYPE ENERGY             ! single point, no forces (default would be ENERGY_FORCE)
  PRINT_LEVEL MEDIUM          ! SILENT / LOW / MEDIUM / HIGH / DEBUG
&END GLOBAL

&FORCE_EVAL
  METHOD QS                   ! Quickstep (DFT)
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    CHARGE 0                  ! default 0
    MULTIPLICITY 1            ! default 1 for an even number of electrons
    &MGRID
      CUTOFF 400              ! Ry, plane-wave cutoff of the finest grid
      REL_CUTOFF 50           ! Ry, controls the mapping of Gaussians to grid levels
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12     ! overall integral/screening accuracy (default 1.0E-10)
    &END QS
    &POISSON
      PERIODIC NONE           ! isolated molecule: must match &CELL PERIODIC
      POISSON_SOLVER MT       ! Martyna-Tuckerman; box >= 2x the size of the charge density
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC        ! superposition of atomic densities (default)
      EPS_SCF 1.0E-6          ! convergence threshold of the inner SCF loop
      MAX_SCF 30              ! inner OT iterations per outer cycle
      &OT ON                  ! orbital transformation: fast and robust for gapped systems
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON           ! restarts the OT minimisation if not yet converged
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 10.0 10.0 10.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME h2o.xyz
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES     ! shift the molecule to the centre of the box
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6    ! 6 valence electrons: 2s2 2p4
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### 8.2 关键字逐项说明

| 关键字 | 本例取值 | 为什么这样选 |
|---|---|---|
| `RUN_TYPE` | `ENERGY` | 只要能量；默认值 `ENERGY_FORCE` 会多算一遍力 |
| `BASIS_SET_FILE_NAME` / `POTENTIAL_FILE_NAME` | `BASIS_MOLOPT` / `GTH_POTENTIALS` | 经典 MOLOPT 基组 + GTH 赝势组合，所有版本都有；新项目的 UZH 组合见 [第3章](03-gpw-basis-pseudo.md) |
| `CHARGE` / `MULTIPLICITY` | 0 / 1 | 中性闭壳层水分子；奇数电子体系必须加 `UKS` 并设正确多重度（[第8章](08-dispersion-spin-dftu.md)） |
| `CUTOFF` / `REL_CUTOFF` | 400 / 50 Ry | 对 MOLOPT 基组的小分子是常见的起点，需按 [第4章](04-cutoff-convergence.md) 收敛 |
| `EPS_DEFAULT` | 1E-12 | 比默认 1E-10 更严格，能量和力更稳定；代价很小 |
| `POISSON_SOLVER MT` | — | 孤立体系去掉周期性镜像的静电相互作用；MT 要求盒子至少是分子电荷分布尺寸的 2 倍，水分子约 2–3 Å，10 Å 的盒子足够 |
| `&OT` | DIIS + FULL_SINGLE_INVERSE | 有带隙体系的首选 SCF 算法；SCF 策略详见 [第5章](05-scf-strategy.md) |
| `&OUTER_SCF` | 内层 30 × 外层 10 | OT 的“重启”机制，提高难收敛体系的稳健性 |
| `&CENTER_COORDINATES` | — | 自动把分子放到盒子中央，避免分子贴着盒子边界 |

### 8.3 运行

```bash
cd CP2K/examples/02-first-calc-h2o-si
cp2k.psmp -c h2o_energy.inp                     # 语法检查
rm -f h2o_energy.out                            # -o 是追加写入
OMP_NUM_THREADS=2 cp2k.ssmp -i h2o_energy.inp -o h2o_energy.out
# 或 MPI：mpirun -np 2 -x OMP_NUM_THREADS=1 cp2k.psmp -i h2o_energy.inp -o h2o_energy.out
```

### 8.4 读懂输出文件

CP2K 的输出按时间顺序分为几个部分。建议第一次完整地翻一遍 `h2o_energy.out`，之后用 `grep` 定位。

**(1) 程序与环境信息**：版本、`cp2kflags`、数据目录、输入文件名、`GLOBAL|` 行（项目名、任务类型、FFT/对角化库）、并行信息：

```bash
grep -E "CP2K\| version string|CP2K\| Data directory path|GLOBAL\| (Project name|Run type)" h2o_energy.out
grep "DBCSR| MPI: Number of processes\|DBCSR| OMP: Current number of threads" h2o_energy.out
```

**(2) 体系与方法信息**：原子 KIND、基组与赝势的详细参数、电子数、多重网格信息：

（示意输出，数值仅供参考）

```text
 -------------------------------------------------------------------------------
 ----                             MULTIGRID INFO                            ----
 -------------------------------------------------------------------------------
 count for grid        1:            420          cutoff [a.u.]          200.00
 count for grid        2:           3010          cutoff [a.u.]           66.67
 count for grid        3:           1830          cutoff [a.u.]           22.22
 count for grid        4:             95          cutoff [a.u.]            7.41
 total gridlevel count  :           5355
```

`cutoff [a.u.]` 以 Hartree 为单位（400 Ry = 200 Ha），四级网格的截断依次除以 `PROGRESSION_FACTOR`（默认 3）。`count for grid` 是映射到各级网格的高斯函数乘积个数，收敛测试时会用到（[第3章](03-gpw-basis-pseudo.md)、[第4章](04-cutoff-convergence.md)）。

**(3) SCF 迭代**：

（示意输出，数值仅供参考）

```text
  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 OT DIIS     0.15E+00    0.1     0.03128500       -17.0802151522 -1.71E+01
     2 OT DIIS     0.15E+00    0.1     0.01601743       -17.1386102385 -5.84E-02
   ...
    10 OT DIIS     0.15E+00    0.1     0.00000081       -17.1625318470 -3.10E-11

  *** SCF run converged in    10 steps ***
```

| 列 | 含义 |
|---|---|
| `Update method` | 本步使用的算法（OT DIIS / OT CG / Diag. 等） |
| `Convergence` | 收敛指标，降到 `EPS_SCF` 以下即收敛 |
| `Total energy` / `Change` | 当前总能量（Ha）及与上一步的差 |

收敛后会打印能量分解（`Overlap energy of the core charge distribution`、`Self energy of the core charge distribution`、`Core Hamiltonian energy`、`Hartree energy`、`Exchange-correlation energy`、`Total energy`），以及 `Electronic density on regular grids` 等电荷积分。**电子密度积分应非常接近 −（价电子数）**，偏差大说明网格太粗。

> [!CAUTION]
> 自 CP2K 2024.1 起，SCF 不收敛时**默认中止计算**并打印 `SCF run NOT converged. To continue the calculation regardless, please set the keyword IGNORE_CONVERGENCE_FAILURE.`。不要一看到这个就加 `IGNORE_CONVERGENCE_FAILURE`，先按 [第5章](05-scf-strategy.md) 排查结构、电荷/多重度与 SCF 设置。

**(4) 布居分析**：默认打印 Mulliken 电荷：

（示意输出，数值仅供参考）

```text
                     Mulliken Population Analysis

 #  Atom  Element  Kind  Atomic population                           Net charge
       1     O        1          6.5xxxxx                             -0.5xxxxx
       2     H        2          0.7xxxxx                              0.2xxxxx
       3     H        2          0.7xxxxx                              0.2xxxxx
 # Total charge                              8.000000                 -0.000000
```

注意这里的“Atomic population”是**价电子**布居（O 为 6 个价电子），Mulliken 电荷强烈依赖基组，只能作定性比较；更合理的电荷分析见 [第11章](11-cube-charge-workfunction.md)。

**(5) 最终能量**：

```bash
grep "ENERGY| Total FORCE_EVAL" h2o_energy.out
```

| 版本 | 输出格式 |
|---|---|
| 2025.1 及以后 | `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]      -17.16...` |
| 2024.x 及以前 | `ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:          -17.16...` |

2025.1 起，这一行的能量单位可用 `&FORCE_EVAL/&PRINT/&PROGRAM_RUN_INFO` 中的 `ENERGY_UNIT` 修改（如 `eV`），默认 Hartree。

**(6) 结尾**：警告数、内存峰值、参考文献、timing 表、结束标记：

```bash
grep "The number of warnings for this run is" h2o_energy.out
grep "MEMORY| Estimated peak process memory" h2o_energy.out
grep -A12 "T I M I N G" h2o_energy.out | head -15      # timing 表开头
grep "PROGRAM ENDED AT" h2o_energy.out                  # 有这一行才说明正常结束
```

timing 表按子程序列出调用次数和耗时；表中名为 `CP2K` 的那一行给出总耗时。大体系调优时看哪个子程序最耗时（第18章）。输出末尾的 `REFERENCES` 列出了本次计算用到的方法文献，发表论文时应引用。

> [!TIP]
> 本节的解析脚本可以一次性提取 SCF 是否收敛、步数、总能量、Mulliken 电荷、警告与耗时：用 `python ../scripts/cp2k_outparse.py -h` 查看用法。

### 8.5 运行目录里多出的文件

| 文件 | 来源 | 用途 |
|---|---|---|
| `h2o_energy.out` | `-o` 指定 | 主输出 |
| `h2o_energy-RESTART.wfn` | `&SCF/&PRINT/&RESTART`（默认开启） | 收敛的波函数（二进制），可作下一次计算的初猜（`SCF_GUESS RESTART`） |
| `h2o_energy-RESTART.wfn.bak-1` | 同上 | 上一次写出的备份 |

---

## 9. 实战二：H2O 能量与原子力（RUN_TYPE ENERGY_FORCE）

[`h2o_energy_force.inp`](../examples/02-first-calc-h2o-si/h2o_energy_force.inp) 与上例的区别：

```text
&GLOBAL
  PROJECT h2o_ef
  RUN_TYPE ENERGY_FORCE       ! energy + analytic forces on all atoms
&END GLOBAL

&FORCE_EVAL
  METHOD QS
  &PRINT
    &FORCES ON                ! print forces (default unit: hartree/bohr)
    &END FORCES
  &END PRINT
  &DFT
    ...
    &SCF
      EPS_SCF 1.0E-7          ! tighter SCF: forces are more sensitive than energies
      ...
      &PRINT
        &RESTART ON           ! write h2o_ef-RESTART.wfn (on by default, shown for clarity)
        &END RESTART
      &END PRINT
    &END SCF
    ...
```

运行与查看：

```bash
rm -f h2o_ef.out
cp2k.ssmp -i h2o_energy_force.inp -o h2o_ef.out
grep -A6 "FORCES| Atomic forces" h2o_ef.out     # 2025.1 及以后
grep -A6 "ATOMIC FORCES" h2o_ef.out             # 2024.x 及以前
```

（示意输出，数值仅供参考；2026.2 格式）

```text
 FORCES| Atomic forces [hartree/bohr]
 FORCES| Atom x                 y                 z                 |f|
 FORCES|       1  0.00000000E+00  1.2xxxxxxxE-07 -1.1xxxxxxxE-02  1.1xxxxxxxE-02
 FORCES|       2  0.00000000E+00 -8.3xxxxxxxE-03  5.6xxxxxxxE-03  1.0xxxxxxxE-02
 FORCES|       3  0.00000000E+00  8.3xxxxxxxE-03  5.6xxxxxxxE-03  1.0xxxxxxxE-02
 FORCES| Sum            0.00000000E+00  1.xxxxxxxxE-07  2.xxxxxxxxE-05
 FORCES| Total atomic force                                         2.xxxxxxxxE-05
```

判断要点：

1. **力之和应接近 0**（孤立分子没有外力）。`Sum` 明显偏离 0 说明网格太粗（`CUTOFF`/`REL_CUTOFF` 太低）或 `EPS_DEFAULT` 太松——这是检查网格质量的一个廉价指标（“egg-box” 效应）。
2. 实验构型并不是 PBE/DZVP 的能量极小点，所以各原子受力不为零；要得到该方法下的平衡结构需做几何优化（[第7章](07-geo-cell-opt.md)），优化收敛判据 `MAX_FORCE` 默认 4.5E-4 Ha/bohr。
3. 单位换算：1 Ha/bohr = 51.422 eV/Å。例如 0.011 Ha/bohr ≈ 0.57 eV/Å。

> [!NOTE]
> 2025.1 起，`&FORCE_EVAL/&PRINT/&FORCES` 的屏幕输出改为带 `FORCES|` 前缀的新格式，并可用 `FORCE_UNIT` 关键字（2026.2 输入参考中存在）改单位；若设置了 `FILENAME`，写入单独文件时仍使用稳定的旧格式（`ATOMIC FORCES in [a.u.]`、`SUM OF ATOMIC FORCES`），便于脚本解析。

---

## 10. 实战三：周期性体系——Si 8 原子立方晶胞

[`si_bulk8_energy.inp`](../examples/02-first-calc-h2o-si/si_bulk8_energy.inp) 演示周期体系的写法，以及用预处理器 `@IF` 选择性地启用 k 点：

```text
&GLOBAL
  PROJECT si8_k${NK-0}
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD QS
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 300
      REL_CUTOFF 50
    &END MGRID
    ! The KPOINTS section is only read when NK is set to a non-zero value.
    @IF ${NK-0} /= 0
    &KPOINTS
      SCHEME MONKHORST-PACK ${NK-0} ${NK-0} ${NK-0}
    &END KPOINTS
    @ENDIF
    &SCF
      EPS_SCF 1.0E-7
      MAX_SCF 100
      ADDED_MOS 20            ! empty bands needed by the smearing
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 5.431 5.431 5.431
      PERIODIC XYZ            ! default
    &END CELL
    &COORD
      SCALED T                ! fractional coordinates
      Si  0.00  0.00  0.00
      Si  0.00  0.50  0.50
      ...（共 8 个原子，见示例文件）
    &END COORD
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

与分子例子的不同点：

| 项目 | 分子（H2O） | 周期体系（Si） |
|---|---|---|
| `&CELL PERIODIC` | `NONE` | `XYZ`（默认） |
| `&POISSON` | `PERIODIC NONE` + `MT` | 不写（默认周期性 FFT 求解器） |
| 坐标 | 笛卡尔（Å） | `SCALED T` 分数坐标 |
| 基组 | `DZVP-MOLOPT-GTH` | `DZVP-MOLOPT-SR-GTH`（短程版，凝聚相更高效，[第3章](03-gpw-basis-pseudo.md)） |
| SCF | OT | 对角化 + Broyden 混合 + Fermi–Dirac 展宽（同时兼容 k 点） |

运行与比较：

```bash
rm -f si8_k0.out si8_k4.out
cp2k.ssmp -i si_bulk8_energy.inp -o si8_k0.out              # Gamma 点
cp2k.ssmp -E NK=4 -i si_bulk8_energy.inp -o si8_k4.out      # 4x4x4 Monkhorst-Pack
grep "BRILLOUIN|" si8_k4.out                                 # 确认 k 点方案和网格
grep "ENERGY| Total FORCE_EVAL" si8_k0.out si8_k4.out
```

8 原子 Si 晶胞边长只有 5.43 Å，Gamma 点能量与 4×4×4 网格的结果会有明显差别（通常每原子几十 meV 以上）；继续加密到 6×6×6、8×8×8 可以看到收敛。k 点收敛方法见 [第6章](06-periodic-surface-kpoints.md)，也可用本节脚本 `conv_kpoints.py`（`python ../scripts/conv_kpoints.py -h` 查看用法）自动化。

> [!NOTE]
> 为什么这里不用 OT？**在 CP2K 2026.2 及以前，OT 不支持 k 点**（2026.2 手册 k 点兼容性表中标为 Unsupported；trunk 开发版已开始支持）。所以只要输入可能启用 k 点，就用对角化。`ADDED_MOS` 为展宽提供空轨道，其取值原则见 [第5章](05-scf-strategy.md)。

---

## 11. 实战四：金红石 TiO2 的能量与力（改造仓库原有示例）

[`examples/02-energy-force-rutile/`](../examples/02-energy-force-rutile/README.md) 中的 `input_file.inp` 是本仓库最早的 CP2K 示例：金红石 TiO2 原胞（2 Ti + 4 O），PBE，`DZVP-MOLOPT-SR-GTH` + `GTH-PBE-q6`/`GTH-PBE-q12`，OT，`RUN_TYPE ENERGY_FORCE`。它的结构写法值得一看：

```text
  &SUBSYS
    &CELL
      ABC [angstrom]    4.593 4.593 2.959
    &END CELL
    &COORD
@include rutile.xyz
    &END COORD
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND Ti
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q12
    &END KIND
  &END SUBSYS
```

- `@include rutile.xyz`（预处理器指令不区分大小写）把一个**没有表头**的坐标文件（6 行 `元素 x y z`）原样插入 `&COORD`；
- Ti 用 `GTH-PBE-q12`：价电子包括 3s²3p⁶3d²4s² 共 12 个（半芯态赝势），基组 `DZVP-MOLOPT-SR-GTH` 与之匹配；
- `&SCF` 中 `SCF_GUESS RESTART`：第一次运行时 `Universality-RESTART.wfn` 不存在，CP2K 会警告 `... This file does not exist ... Calculation continues using ATOMIC GUESS.` 并自动改用原子初猜。

### 11.1 问题：6 原子原胞只用 Gamma 点

原输入没有 `&KPOINTS`，即 **Gamma-only**。CP2K 以大超胞 + Gamma 点见长，但这个原胞只有 4.59 × 4.59 × 2.96 Å，倒格子很“大”，只取 k = 0 一个点对布里渊区积分的误差非常大——能量、力、电子结构都不可信。

> [!WARNING]
> 经验判断：晶胞**每个方向的长度都在 10–15 Å 以上**（且不是金属）时，Gamma 点通常才够用；更小的晶胞必须做 k 点收敛，或者扩成超胞。金属需要更密的采样。

两种改进方法都放在示例目录中：

**方法 A：原胞 + k 点**（[`rutile_kpoints.inp`](../examples/02-energy-force-rutile/rutile_kpoints.inp)）

```text
    &KPOINTS
      ! reciprocal lengths: 1/a = 0.218 1/A, 1/c = 0.338 1/A -> use more points along c
      SCHEME MONKHORST-PACK ${KX-4} ${KY-4} ${KZ-6}
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 100
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.3
        NBUFFER 8
      &END MIXING
    &END SCF
```

c 轴较短，所以沿 c 方向取更多 k 点（4×4×6）。由于 OT 不支持 k 点（≤ 2026.2），改用对角化 + Broyden 混合；TiO2 是绝缘体，不需要展宽。

**方法 B：超胞 + Gamma 点**（[`rutile_supercell_gamma.inp`](../examples/02-energy-force-rutile/rutile_supercell_gamma.inp)）

```text
    &CELL
      ABC [angstrom] 4.593 4.593 2.959        ! the UNIT cell
      MULTIPLE_UNIT_CELL ${NX-2} ${NY-2} ${NZ-3}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      MULTIPLE_UNIT_CELL ${NX-2} ${NY-2} ${NZ-3}
    &END TOPOLOGY
    &COORD
@INCLUDE rutile.xyz
    &END COORD
```

`&CELL` 与 `&TOPOLOGY` 中的 `MULTIPLE_UNIT_CELL` 必须一致：前者扩晶胞，后者复制坐标。2×2×3 超胞含 72 个原子（9.19 × 9.19 × 8.88 Å），仍可用 OT。**Gamma 点超胞 n1×n2×n3 与原胞上以 Gamma 为中心的 n1×n2×n3 k 网格采样的是同一组 k 点**，因此 `-E NX=3 -E NY=3 -E NZ=3`（162 原子）的每 TiO2 能量应与 `rutile_kpoints.inp` 用 `-E KX=3 -E KY=3 -E KZ=3`（奇数网格包含 Gamma 点）的结果非常接近——这是检验你是否理解 k 点的一个好练习。

### 11.2 运行与对比

```bash
cd CP2K/examples/02-energy-force-rutile
for f in input_file rutile_kpoints rutile_supercell_gamma; do cp2k.psmp -c $f.inp; done

rm -f input_file.out rutile_k446.out rutile_sc223.out
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i input_file.inp -o input_file.out
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i rutile_kpoints.inp -o rutile_k446.out
mpirun -np 8 -x OMP_NUM_THREADS=1 cp2k.psmp -i rutile_supercell_gamma.inp -o rutile_sc223.out

# 每 TiO2 单元的能量（原胞 2 个单元，2x2x3 超胞 24 个单元）
for f in input_file.out:2 rutile_k446.out:2 rutile_sc223.out:24; do
  out=${f%%:*}; n=${f##*:}
  awk -v n=$n -v f=$out '/ENERGY\| Total FORCE_EVAL/{e=$NF} END{printf "%-22s E/TiO2 = %.8f Ha\n", f, e/n}' $out
done
```

集群上用目录中的作业脚本：`sbatch --export=ALL,INP=rutile_kpoints.inp job_slurm.sh`（Slurm），或 `bsub < job_lsf.sh`（LSF）。原仓库的 `cp2k.slurm` 名为 Slurm、实为 LSF 脚本，且运行的是不存在的 `input.inp`，现已废弃并由这两个脚本替代。

关于原子力：Ti 位于反演中心，受力由对称性严格为零；O 位于 $(u,u,0)$ 一类位置，沿 [110] 方向的力取决于内坐标 $u$ 是否处于该方法下的平衡值。比较 Gamma 点原胞与 k 点结果中 O 的受力，可以直观看到采样不足带来的误差。

---

## 12. 输出文件与重启文件一览

| 文件名 | 何时产生 | 内容与用途 |
|---|---|---|
| `<PROJECT>.out`（`-o` 指定） | 总是 | 主输出 |
| `<PROJECT>-RESTART.wfn`（及 `.bak-N`） | Gamma 点 SCF 计算（默认开启） | 波函数重启文件；`SCF_GUESS RESTART` 读入 |
| `<PROJECT>-RESTART.kp` | k 点计算 | k 点波函数重启文件；**与 `.wfn` 不能互换** |
| `<PROJECT>-1.restart`（及 `.bak-N`） | `GEO_OPT`、`CELL_OPT`、`MD`、`BAND` 等有 `&MOTION` 的任务 | 完整输入格式的重启文件，配合 `&EXT_RESTART` 续算；**`ENERGY`/`ENERGY_FORCE` 单点计算不写这个文件** |
| `<PROJECT>-pos-1.xyz` | 优化、MD | 轨迹（每步坐标） |
| `<PROJECT>-1.ener` | MD | 每步动能、温度、势能、守恒量 |
| `<PROJECT>-frc-1.xyz`、`-vel-1.xyz`、`-1.cell` | MD/优化中按 `&MOTION/&PRINT` 请求 | 力、速度、晶胞 |
| `*.pdos`、`*.cube`、`*.molden`、`*.bs` | 按 `&DFT/&PRINT` 请求 | 态密度、实空间性质、轨道、能带 |

> [!TIP]
> 同一目录中用同一个 `PROJECT` 名跑两个不同的计算，会互相覆盖重启文件（甚至错误地读入别的体系的波函数）。**一个计算一个目录，或一个计算一个 PROJECT 名。**

---

## 13. RUN_TYPE 一览

| RUN_TYPE | 做什么 | 需要的段 | 详见 |
|---|---|---|---|
| `ENERGY` | 单点能量 | `&FORCE_EVAL` | 本章 |
| `ENERGY_FORCE`（**默认**） | 单点能量 + 原子力（可加应力） | `&FORCE_EVAL`（+ `&PRINT &FORCES`） | 本章 |
| `GEO_OPT`（= `GEOMETRY_OPTIMIZATION`） | 几何优化（也可做 Dimer 过渡态） | `&MOTION/&GEO_OPT` | [第7章](07-geo-cell-opt.md)、[第13章](13-neb-dimer.md) |
| `CELL_OPT` | 晶胞 + 原子位置同时优化 | `&MOTION/&CELL_OPT` + `STRESS_TENSOR` | [第7章](07-geo-cell-opt.md) |
| `MD`（= `MOLECULAR_DYNAMICS`） | 分子动力学（NVE/NVT/NPT…） | `&MOTION/&MD` | 第12章、第15章 |
| `BAND` | NEB / CI-NEB 等链方法 | `&MOTION/&BAND` | [第13章](13-neb-dimer.md) |
| `VIBRATIONAL_ANALYSIS`（= `NORMAL_MODES`） | 有限差分 Hessian、频率、IR 强度 | 顶层 `&VIBRATIONAL_ANALYSIS` | [第14章](14-vibration-thermo.md) |
| `LINEAR_RESPONSE`（= `LR`） | 线性响应性质（NMR、极化率等） | `&FORCE_EVAL/&PROPERTIES` | — |
| `RT_PROPAGATION`、`EHRENFEST_DYN` | 实时 TDDFT、Ehrenfest 动力学 | `&DFT/&REAL_TIME_PROPAGATION` | — |
| `MC`、`PINT`、`BSSE`、`DEBUG`、`NONE` 等 | 蒙特卡洛、路径积分、BSSE、数值验证力、什么都不做 | 各自的段 | 输入参考 |

`ENERGY` 还有别名 `WFN_OPT`/`WAVEFUNCTION_OPTIMIZATION`。`DEBUG` 用有限差分检验解析力是否正确，开发和验证新功能时很有用。

---

## 参数选择与收敛建议

第一次搭建计算时可以参考下表作为起点（之后必须针对具体体系做收敛测试）：

| 参数 | 起点值 | 说明 |
|---|---|---|
| `CUTOFF` | 400 Ry（MOLOPT 基组的有机/水体系）；400–600 Ry（含过渡金属、要算应力） | 按 [第4章](04-cutoff-convergence.md) 收敛，判据如每原子 < 1 meV |
| `REL_CUTOFF` | 50–60 Ry | 默认 20 Ry 偏低 |
| `EPS_DEFAULT` | 1E-12（要力/MD）；1E-10（默认，只要能量） | 越小越准越慢 |
| `EPS_SCF` | 1E-6（能量）；1E-7（力、优化）；更严（振动分析） | 力比能量更敏感 |
| 孤立分子盒子 | 分子尺寸 + 每侧 ≥ 4–5 Å，且满足 MT 的 2 倍要求 | 配合 `&CENTER_COORDINATES` |
| k 点 | 小原胞必须收敛；大超胞（各方向 ≳ 10–15 Å，绝缘体）可用 Gamma 点 | [第6章](06-periodic-surface-kpoints.md) |
| SCF 算法 | 有带隙 + Gamma 点：OT；金属或 k 点：对角化 + 混合 + 展宽 | [第5章](05-scf-strategy.md) |
| `PRINT_LEVEL` | 调试 `MEDIUM`；长时间 MD `LOW` | — |

---

## 常见问题与报错

| 报错信息 | 原因 | 解决 |
|---|---|---|
| `found an unknown keyword XXX in section YYY` | 关键字拼错、放错了段，或版本不支持（例如 2026.2 删除了 `&CELL_OPT` 中的 `TYPE`） | 报错前 CP2K 通常会打印 `Possible matches for unknown keyword`，列出拼写相近的候选；再查对应版本的输入参考，用 `cp2k.psmp -c` 定位 |
| `unknown subsection XXX of section YYY` | 子段名错误或位置不对（例如 2026.2 下写旧的 `&DFT/&PRINT/&PDOS`） | 同上，注意版本差异（[第10章](10-dos-pdos-band.md)） |
| `A (integer\|floating point\|string) type object was expected, found ...` | 值的类型或个数不对（例如 `ABC` 只写了两个数，或 xyz 文件格式不对） | 按手册给出正确类型和个数的值 |
| `Variable ${X} not defined in file: ...` | 用了未定义的预处理器变量 | 加 `@SET X ...`、用 `-E X=...`，或改写为 `${X-默认值}` |
| `Missing '}' in file: ...` | `${` 未闭合 | 补全花括号 |
| `Error while searching for matching @ENDIF directive ...` | `@IF` 没有对应的 `@ENDIF` | 补上 `@ENDIF`；注意不支持嵌套 |
| `The specified OLD file <h2o.xyz> cannot be opened. It does not exist.` | 坐标/包含文件不在当前工作目录 | 在输入所在目录运行，或写绝对路径 |
| `The requested basis set <...> for element <...> was not found in the basis set files ...` | 基组名不在所列文件中 | 检查拼写与文件（[第3章](03-gpw-basis-pseudo.md)） |
| `Keyword XXX in section YYY should not repeat.` | 不可重复的关键字写了两次（如 `POTENTIAL_FILE_NAME`） | 删去重复行；可重复的只有手册标明 repeats 的关键字（如 `BASIS_SET_FILE_NAME`） |
| `Use the LSD option for an odd number of electrons` | 奇数电子却没开自旋极化 | `&DFT` 中加 `UKS` 并设正确 `MULTIPLICITY`；闭壳层体系出现此错多半是结构少/多了原子 |
| `GEOMETRY wrong or EMAX_SPLINE too small` | 原子重叠或距离过近（常见于 xyz 没有晶胞信息、扩胞后重复原子） | 可视化检查结构与晶胞 |
| `An invalid set of cell vectors was specified. The cell volume is too small` | `&CELL` 写错（单位、漏写某个矢量） | 检查 `ABC`/`A B C` 的数值和单位 |
| SCF 第一步就很怪、能量异常 | `&CELL PERIODIC` 与 `&POISSON PERIODIC` 不一致，或孤立分子盒子太小 | 两处保持一致；增大盒子 |
| 输出中有两次运行的内容 | `-o` 追加写入 | 运行前删除旧输出 |
| 警告 `... RESTART.wfn. This file does not exist ... Calculation continues using ATOMIC GUESS.` | `SCF_GUESS RESTART` 但重启文件还不存在 | 第一次计算属正常；确认后续计算确实读到了想要的文件 |

更多报错及其排查流程见 第18章 与官方排错指南。

---

## 练习

1. 运行 `h2o_energy.inp`，找出：SCF 步数、总能量（换算成 eV）、Mulliken 电荷、四级网格的 `count for grid`、总耗时。把 `PRINT_LEVEL` 改为 `LOW` 与 `HIGH` 再运行，比较输出长度。
2. 用 `cp2k.ssmp -e h2o_energy.inp > full.inp` 导出完整输入，找出 `&SCF` 中 `MAX_SCF`、`&OT` 中 `ENERGY_GAP`、`&MGRID` 中 `NGRIDS` 的默认值。
3. 把 `h2o_energy.inp` 中的 `&POISSON` 段删掉但保留 `&CELL PERIODIC NONE`，或反过来，观察 CP2K 的反应；再把盒子缩小到 5 Å，比较能量变化。
4. 用 `preprocessor_demo.inp` 分别以 `BASIS=SZV/DZVP/TZVP/TZV2P` 运行，比较能量与耗时（基组的意义见 [第3章](03-gpw-basis-pseudo.md)）。试试 `-E XC=BLYP` 为什么无效，并修改输入使其有效（提示：`${XC-PBE}`）。
5. 用 `si_bulk8_energy.inp` 依次运行 `NK=0, 2, 3, 4, 6`，画出每原子能量随 k 网格的变化，判断多密的网格可以达到 1 meV/原子的精度。
6. 完成 §11 的等价性检验：比较 3×3×3 超胞（Gamma 点）与 3×3×3 k 网格的每 TiO2 能量，并解释两者为何应该一致、为何又不完全相同。

---

## 参考资料

**官方文档**

- 第一个计算（H2O 示例与输入语法规则）：<https://manual.cp2k.org/trunk/getting-started/first-calculation.html>
- 输入参考（全部段与关键字）：<https://manual.cp2k.org/trunk/CP2K_INPUT.html>；2026.2 版本：<https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT.html>
- 单位列表：<https://manual.cp2k.org/cp2k-2026_2-branch/units.html>
- k 点（含 2026.2 的功能兼容性表）：<https://manual.cp2k.org/cp2k-2026_2-branch/methods/dft/k-points.html>
- 重启：<https://manual.cp2k.org/trunk/methods/restarting.html>
- 官方排错指南：<https://manual.cp2k.org/trunk/getting-started/troubleshooting.html>
- 旧版 HowTo（单点能量与力 static_calculation 等）：<https://www.cp2k.org/howto>

**相关 GitHub 仓库**（只链接，不复制；更多见 [Resources/github-repos.md](../../Resources/github-repos.md)）

- [cp2k/cp2k](https://github.com/cp2k/cp2k)：`tests/` 下的回归测试输入是查关键字真实写法的最权威来源（例如 `tests/QS/regtest-gpw-1/H2-inpp.inp` 覆盖了预处理器的各种用法）。
- [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples)：官方示例库，含《The CP2K Program Package Made Simple》配套输入。
- [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools)：`cp2klint` 语法检查、`fromcp2k`/`tocp2k` 格式互转、`cp2kgen` 参数扫描。
- [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io)：厦大程俊组中文 Wiki，含 CP2K 入门与能量/力计算。
- [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn)：分级中文教程，单点能与周期 Si 章节可对照阅读（GPL-3.0，仅链接）。
- [ChiahsinChu/cp2k-user-tutorial](https://github.com/ChiahsinChu/cp2k-user-tutorial)：面向零基础用户的入门教程，汇总了中文学习资源。
- [ybx20021218/CP2Kview](https://github.com/ybx20021218/CP2Kview--Visualization-of-CP2K-Output-Files)：浏览器中查看 CP2K 输出与收敛曲线的小工具。

**关键论文**

- J. VandeVondele *et al.*, *Comput. Phys. Commun.* 167, 103 (2005)（Quickstep/GPW 实现）. DOI: [10.1016/j.cpc.2004.12.014](https://doi.org/10.1016/j.cpc.2004.12.014)
- T. D. Kühne *et al.*, *J. Chem. Phys.* 152, 194103 (2020)（CP2K 程序综述）. DOI: [10.1063/5.0007045](https://doi.org/10.1063/5.0007045)
- G. J. Martyna, M. E. Tuckerman, *J. Chem. Phys.* 110, 2810 (1999)（孤立体系的 MT 泊松求解器）. DOI: [10.1063/1.477923](https://doi.org/10.1063/1.477923)
- J. K. Burdett *et al.*, *J. Am. Chem. Soc.* 109, 3639 (1987)（金红石/锐钛矿 TiO2 中子衍射结构）. DOI: [10.1021/ja00246a021](https://doi.org/10.1021/ja00246a021)
- H. J. Monkhorst, J. D. Pack, *Phys. Rev. B* 13, 5188 (1976)（Monkhorst–Pack k 点网格）. DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)

---

[← 上一章](01-install.md) | [返回目录](../README.md) | [下一章 →](03-gpw-basis-pseudo.md)
