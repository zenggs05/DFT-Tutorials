# CP2K 辅助脚本（CP2K helper scripts）

本目录收录 CP2K 计算的前处理、批量收敛测试、作业提交和后处理脚本：从结构文件生成输入（[`cp2k_inpgen.py`](cp2k_inpgen.py)）、CUTOFF/REL_CUTOFF 与 k 点收敛测试（[`conv_cutoff.py`](conv_cutoff.py)、[`conv_kpoints.py`](conv_kpoints.py)）、输出汇总（[`cp2k_outparse.py`](cp2k_outparse.py)）、优化收敛曲线（[`plot_opt.py`](plot_opt.py)）、能带/PDOS 作图（[`band_plot.py`](band_plot.py)、[`pdos_plot.py`](pdos_plot.py)）、cube 文件与功函数（[`cube_tools.py`](cube_tools.py)）、AIMD 分析（[`md_analyze.py`](md_analyze.py)）、NEB 能垒（[`neb_plot.py`](neb_plot.py)）、振动热化学（[`vib_thermo.py`](vib_thermo.py)）、续算（[`restart_update.py`](restart_update.py)），以及 Slurm/PBS 作业模板。脚本解析的输出格式按 CP2K 2026.2 源码核对，同时兼容 2023.x–2025.x 的写法（差异见各节说明）。

> [!NOTE]
> 所有 Python 脚本都支持 `-h` 查看完整参数（带子命令的脚本用 `python <脚本> <子命令> -h`）；两个作业模板可用 `bash <脚本> -h` 显示文件头部的用法说明。本 README 中的参数均直接取自各脚本的 argparse 定义，命令行示例只演示用法，没有附带运行结果。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [cp2k_inpgen.py](cp2k_inpgen.py) | 从 xyz/cif/POSCAR 生成完整的 Quickstep（GPW, PBE）输入，自动查找并核对 MOLOPT 基组与 GTH 赝势 | ase（numpy） | [第2章](../tutorials/02-input-structure.md)、[第3章](../tutorials/03-gpw-basis-pseudo.md)、[第17章](../tutorials/17-tools-ecosystem.md) |
| [conv_cutoff.py](conv_cutoff.py) | CUTOFF / REL_CUTOFF 两步收敛测试：生成目录、本地运行或批量提交、汇总 CSV 并作图 | 标准库；作图需 numpy、matplotlib | [第4章](../tutorials/04-cutoff-convergence.md) |
| [conv_kpoints.py](conv_kpoints.py) | Monkhorst-Pack k 点收敛测试：按长度/k 间距/显式网格生成序列，运行、提交、汇总 | 标准库；作图需 numpy、matplotlib；`--structure` 需 ase | [第6章](../tutorials/06-periodic-surface-kpoints.md) |
| [cp2k_outparse.py](cp2k_outparse.py) | 输出文件一站式汇总：正常结束、SCF、能量、费米能/能隙、Mulliken、WARNING、ABORT、OPT\| 步骤，可导出 JSON/CSV | 仅标准库 | [第2章](../tutorials/02-input-structure.md)、[第5章](../tutorials/05-scf-strategy.md)、[第7章](../tutorials/07-geo-cell-opt.md)、[第18章](../tutorials/18-troubleshooting.md) |
| [plot_opt.py](plot_opt.py) | GEO_OPT / CELL_OPT 收敛过程作图（能量、梯度、位移、压力、晶胞体积） | numpy、matplotlib | [第7章](../tutorials/07-geo-cell-opt.md) |
| [restart_update.py](restart_update.py) | 把 `.restart` 或轨迹中的最新坐标/晶胞写回输入或 xyz，可附加 `&EXT_RESTART` | 仅标准库 | [第5章](../tutorials/05-scf-strategy.md)、[第7章](../tutorials/07-geo-cell-opt.md)、[第12章](../tutorials/12-aimd.md) |
| [band_plot.py](band_plot.py) | 解析 `&BAND_STRUCTURE` 输出（.bs），画能带、判断直接/间接带隙、导出 CSV | numpy、matplotlib | [第10章](../tutorials/10-dos-pdos-band.md) |
| [pdos_plot.py](pdos_plot.py) | 读取旧版 `&PDOS` 与 2026.2 `&DOS/&PDOS/&CURVE` 的 .pdos 文件，展宽、按元素/角动量求和、分自旋作图 | numpy、matplotlib | [第10章](../tutorials/10-dos-pdos-band.md) |
| [cube_tools.py](cube_tools.py) | cube 文件信息/检查、积分、平面平均（V_HARTREE → 真空能级与功函数）、差分电荷密度 | numpy；作图需 matplotlib | [第6章](../tutorials/06-periodic-surface-kpoints.md)、[第11章](../tutorials/11-cube-charge-workfunction.md) |
| [md_analyze.py](md_analyze.py) | AIMD 后处理：.ener 温度/守恒量漂移、RDF 与配位数、MSD 与扩散系数、VACF/VDOS、偶极 IR 谱 | numpy、matplotlib | [第12章](../tutorials/12-aimd.md)、[第15章](../tutorials/15-metadynamics-plumed.md)、[第16章](../tutorials/16-qmmm-xtb.md) |
| [neb_plot.py](neb_plot.py) | RUN_TYPE BAND（NEB/CI-NEB）能量曲线、正/逆能垒、导出最高能副本 xyz | numpy；作图需 matplotlib | [第13章](../tutorials/13-neb-dimer.md) |
| [vib_thermo.py](vib_thermo.py) | 振动分析输出的频率、IR 强度、虚频检查，ZPE 与 U/H/S/G（吸附物或理想气体模型，可选 quasi-RRHO） | numpy | [第14章](../tutorials/14-vibration-thermo.md) |
| [job_cp2k_slurm.sh](job_cp2k_slurm.sh) | Slurm 作业模板（cp2k.psmp，MPI × OpenMP），运行后自动检查结束状态并提示续算方法 | bash、Slurm | [第4章](../tutorials/04-cutoff-convergence.md)、[第12章](../tutorials/12-aimd.md)，以及 [作业调度通用章节](../../Common/tutorials/02-job-scheduler.md) |
| [job_cp2k_pbs.sh](job_cp2k_pbs.sh) | PBS Pro / OpenPBS / Torque 作业模板，功能同上 | bash、PBS/Torque | [第12章](../tutorials/12-aimd.md)，以及 [作业调度通用章节](../../Common/tutorials/02-job-scheduler.md) |

## 安装依赖与通用约定

```bash
# Python >= 3.9
pip install numpy matplotlib      # band_plot / pdos_plot / plot_opt / cube_tools / md_analyze / neb_plot / vib_thermo
pip install ase                   # cp2k_inpgen.py；conv_kpoints.py 的 --structure 选项
# cp2k_outparse.py 和 restart_update.py 只用标准库，无需安装任何包
```

通用约定：

- **作图**：matplotlib 只在需要作图时才导入，没有图形界面的集群上会自动用 Agg 后端保存图片，加 `--show` 时才弹出窗口。多数脚本的 `-o/--output` 可设为 `none`，这样只打印表格、写 CSV，不出图；[`conv_cutoff.py`](conv_cutoff.py)、[`conv_kpoints.py`](conv_kpoints.py)、[`plot_opt.py`](plot_opt.py) 改用 `--no-plot`。
- **单位**：CP2K 输出的能量是 hartree，长度是 bohr（cube 文件、梯度）或 Å（xyz、`CELL|` 行）。脚本内部按 CODATA 2018 的 1 Ha = 27.211386245988 eV 换算，各节会写明输入/输出单位。
- **文件命名**：示例中的 `<PROJECT>` 指 `&GLOBAL/PROJECT`。CP2K 生成的文件名形如 `<PROJECT>-1.restart`、`<PROJECT>-pos-1.xyz`、`<PROJECT>-1.ener`、`<PROJECT>-1.cell`、`<PROJECT>-k1-1.pdos`。
- **版本**：第一次在新版本 CP2K 上使用时，建议先用 [`cp2k_outparse.py`](cp2k_outparse.py) 看一个小算例的解析结果是否完整。如有字段缺失，欢迎提 Issue 并附上对应的输出片段。

> [!IMPORTANT]
> 生成输入或批量提交前，先用 `cp2k.psmp -c input.inp`（只检查输入语法，不计算）检查一遍。[`cp2k_inpgen.py`](cp2k_inpgen.py) 的结果、`conv_*` 生成的目录也一样要检查。

---

## 1. cp2k_inpgen.py — 从结构文件生成输入

**用途**：读入任意 ASE 能读的结构（xyz、cif、POSCAR、extxyz……），写出一个完整、可直接 `cp2k.psmp -c` 检查的 Quickstep 输入。脚本会处理 `&CELL`/`PERIODIC`、基组与赝势配对、SCF 策略、k 点、色散、自旋、偶极校正，以及 GEO_OPT / CELL_OPT / MD 的 `&MOTION`。适合快速起草输入，再按第 2–12 章的说明手工细调。

**依赖**：`pip install ase`（numpy 随 ase 安装）。ase 在运行时才导入，缺失时会给出提示 `cp2k_inpgen.py needs ASE to read structures: pip install ase`。

**用法示例**：

```bash
python cp2k_inpgen.py h2o.xyz                                   # isolated molecule, OT, ENERGY
python cp2k_inpgen.py h2o.xyz --run-type GEO_OPT --vacuum 6
python cp2k_inpgen.py Si.cif --kpoints 4 4 4 --smear 300 --run-type CELL_OPT
python cp2k_inpgen.py POSCAR --kpoints 6 6 1 --d3bj --surface-dipole --run-type GEO_OPT
python cp2k_inpgen.py fe.cif --uks --multiplicity 5 --magmom Fe=4 --smear 500 --kpoints 8 8 8
python cp2k_inpgen.py water64.xyz --box 12.42 12.42 12.42 --periodic XYZ --run-type MD \
       --temperature 330 --steps 20000 --d3bj
python cp2k_inpgen.py TiO2.cif --pdos --version 2025.2          # legacy &PRINT/&PDOS syntax
python cp2k_inpgen.py Cu.cif --q Cu=11 --basis TZVP --data-dir /opt/cp2k/data
cp2k.psmp -c h2o.inp                                            # always syntax-check the result
```

**参数**（括号内为默认值）：

| 类别 | 参数 | 说明 |
|---|---|---|
| 通用 | `structure` | 结构文件（位置参数） |
| | `--format FMT` | ASE 格式名，猜不出格式时指定（如 `vasp`） |
| | `-o/--output OUT` | 输出输入文件名（`<project>.inp`） |
| | `--project NAME` | `PROJECT` 名（结构文件名去掉扩展名） |
| | `--force` | 覆盖已有文件 |
| | `--version` | 目标 CP2K 版本：`2026.2`（默认）、`2025.2`、`2024.1` 等 `YYYY.N` 形式，或 `trunk`（`master`/`dev` 同义）；决定 DOS/PDOS 与 CELL_OPT 的写法 |
| 体系 | `--run-type` | `ENERGY`（默认）/ `ENERGY_FORCE` / `GEO_OPT` / `CELL_OPT` / `MD` |
| | `--periodic XYZ\|NONE` | 周期性（默认由文件判断：有晶胞 → XYZ，无晶胞 → NONE） |
| | `--box A B C` | 无晶胞结构（xyz）使用的正交盒子（Å），分子会被居中；文件自带晶胞时报错 |
| | `--vacuum` | 孤立分子：每个方向盒长 = max(2×分子尺寸, 分子尺寸 + 2×vacuum)（5.0 Å） |
| | `--cubic-box` | 孤立分子使用立方盒 |
| | `--charge N` | 总电荷（0） |
| | `--uks` | 自旋极化 |
| | `--multiplicity M` | 2S+1，不等于 1 时自动开启 UKS |
| | `--magmom El=M ...` | 按元素写初始 `MAGNETIZATION`（配合 `--uks`） |
| | `--surface-dipole` | `SURFACE_DIPOLE_CORRECTION T`，`SURF_DIP_DIR Z`（真空沿 z 的 slab） |
| | `--coord-file XYZ` | 坐标写到单独的 xyz，输入中用 `&TOPOLOGY` 引用，而不是内联 `&COORD` |
| 基组/赝势 | `--basis` | `DZVP`（默认）/ `TZVP` / `TZV2P`，对应 `<BASIS>-MOLOPT-PBE-GTH-qN` |
| | `--q El=N ...` | 指定某元素的价电子数 q，如 `Cu=11 Na=9` |
| | `--kind-basis El=NAME ...` | 直接指定某元素的基组名 |
| | `--kind-potential El=NAME ...` | 直接指定某元素的赝势名 |
| | `--basis-file` | `BASIS_SET_FILE_NAME`（`BASIS_MOLOPT_UZH`） |
| | `--potential-file` | `POTENTIAL_FILE_NAME`（`POTENTIAL_UZH`） |
| | `--data-dir DIR` | CP2K data 目录（未给时依次找 `$CP2K_DATA_DIR`、脚本旁的 `../data`） |
| | `--no-check` | 名字核对失败时只警告、照样写出 |
| SCF/数值 | `--cutoff` / `--rel-cutoff` | `CUTOFF`（400 Ry）/ `REL_CUTOFF`（60 Ry），正式计算前应按[第4章](../tutorials/04-cutoff-convergence.md)做收敛测试 |
| | `--eps-scf` | `EPS_SCF`（1e-6） |
| | `--scf ot\|diag` | SCF 求解器（默认 OT；给了 `--kpoints` 或 `--smear` 时自动改为 diag） |
| | `--smear T_K` | Fermi-Dirac 展宽的电子温度（K），隐含 diag |
| | `--added-mos N` | 对角化时的 `ADDED_MOS`（默认 max(10, 占据态数的 20%)） |
| | `--kpoints N1 N2 N3` | Monkhorst-Pack 网格，隐含 diag |
| | `--d3bj` | DFT-D3(BJ) 色散（PBE 参数） |
| | `--stress` | 解析应力 + `STRESS_TENSOR` 打印（CELL_OPT 自动打开） |
| | `--pdos` | 添加 DOS/PDOS 打印块（语法随 `--version`） |
| 优化/MD | `--optimizer` | `BFGS`（默认）/ `LBFGS` / `CG` |
| | `--max-iter` | `MAX_ITER`（200） |
| | `--max-force` | `MAX_FORCE`，hartree/bohr（4.5e-4） |
| | `--pressure` | CELL_OPT 的 `EXTERNAL_PRESSURE`，bar（1.01325） |
| | `--ensemble` | `NVT`（默认，CSVR 恒温器）/ `NVE` |
| | `--temperature` / `--timestep` / `--steps` / `--timecon` | MD 温度 300 K / 步长 0.5 fs / 步数 1000 / CSVR 时间常数 100 fs |

**输出与行为要点**：

- 只生成 PBE 输入（基组和赝势名均为 `*-PBE-*`）。杂化泛函/ADMM、DFT+U、SCAN 等需在生成后手工修改，见 [第8章](../tutorials/08-dispersion-spin-dftu.md) 和 [第9章](../tutorials/09-hybrid-admm.md)。
- 晶胞写成 `ABC` + `ALPHA_BETA_GAMMA`，笛卡尔坐标转换到 CP2K 约定的下三角取向（a 沿 x，b 在 xy 平面内），因此坐标与晶胞取向一致。
- 孤立分子（`PERIODIC NONE`）会建一个居中的盒子，并写 `&POISSON` 的 `PERIODIC NONE` + `PSOLVER MT`。MT 求解器要求盒子足够大，见 [第6章](../tutorials/06-periodic-surface-kpoints.md)。
- **基组/赝势核对**：对每个元素，取 `BASIS_MOLOPT_UZH` 中 `<BASIS>-MOLOPT-PBE-GTH-q*` 与 `POTENTIAL_UZH` 中 `GTH-PBE-q*` 的 q 值交集。有多个可选时取最大值（含半芯态），并列出其他选择，可用 `--q` 改选。脚本还会核对两者的 `-qN` 后缀一致、名字确实存在于数据文件中。
- 找不到数据文件时，只有 H、C、N、O、F、Si、P、S、Cl、Br、I（价电子数无歧义）能自动处理，其他元素必须给 `--q El=N`。
- **SCF 策略**：默认 OT（+ `OUTER_SCF`）。给了 `--kpoints` 或 `--smear` 时自动改为对角化 + Broyden mixing，并写 `ADDED_MOS`。对 2026.2 及更早版本，`--scf ot` 与 k 点或展宽同时出现会直接报错退出，因为这些版本的 OT 不支持 k 点和展宽。
- 价电子总数为奇数却没开 `--uks` 时会报错（否则 CP2K 会停在 "Use the LSD option for an odd number of electrons"）。多重度与电子数矛盾时同样报错。
- **版本差异**：`--version` ≥ 2026.2 时写 `&DFT/&PRINT/&DOS` + `&PDOS`/`&CURVE`/`&BROADEN`（`DELTA_E` 按 2026.2 源码不带单位），CELL_OPT 不写 `TYPE`。更早的版本写旧版 `&PRINT/&PDOS`（使用 k 点时省略该块），CELL_OPT 写 `TYPE DIRECT_CELL_OPT`。把 2026.2 之前的 CELL_OPT 输入直接拿到 2026.2 运行，会报 "The specified keyword 'TYPE' is not available anymore"，详见 [第7章](../tutorials/07-geo-cell-opt.md) 和 [第18章](../tutorials/18-troubleshooting.md)。
- D3(BJ) 与 MD/CSVR 块与 [第8章](../tutorials/08-dispersion-spin-dftu.md)、[第12章](../tutorials/12-aimd.md) 的示例输入写法一致。

**局限**：不做对称性分析，不生成能带路径（`&BAND_STRUCTURE` 请参考[第10章](../tutorials/10-dos-pdos-band.md)手写）；`--box` 只支持正交盒；生成的 CUTOFF/REL_CUTOFF 只是起点，不是收敛值。

---

## 2. conv_cutoff.py — CUTOFF / REL_CUTOFF 收敛测试

**用途**：实现 [第4章](../tutorials/04-cutoff-convergence.md) 的两步法。先固定 REL_CUTOFF 扫 CUTOFF，再固定 CUTOFF 扫 REL_CUTOFF。脚本用 `@SET` 变量替换模板中的数值，生成 `cutoff_<值>Ry/` 或 `relcutoff_<值>Ry/` 目录，可以本地运行或逐个提交作业，最后汇总能量差、SCF 收敛情况、各级网格的 `count for grid` 和（如有打印）压力。

**依赖**：标准库即可完成 gen/run/submit/collect；作图需 `pip install numpy matplotlib`（没有时用 `--no-plot`）。

**子命令与参数**：

| 子命令 | 作用 | 参数 |
|---|---|---|
| 公共 | — | `--scan cutoff\|relcutoff`（必需）、`--inp-name`（`input.inp`）、`--out-name`（`output.out`）、`--force`（重新生成/重新提交已有输出的点） |
| `gen` | 只生成目录 | `-t/--template`（必需）、`--values V ...`（Ry，必需）、`--cutoff`（扫 relcutoff 时固定的 CUTOFF）、`--rel-cutoff`（扫 cutoff 时固定的 REL_CUTOFF）、`--extra FILES`（复制到每个目录的文件，如 xyz）、`--var-cutoff`（`CUTOFF`）、`--var-rel`（`REL_CUTOFF`） |
| `run` | gen + 本地运行 + collect | gen 的全部参数 + `--cp2k`（`"cp2k.psmp"`，可含 MPI 启动器）、`--npar`（同时运行的点数，1）+ collect 参数 |
| `submit` | 逐个目录提交作业 | `--job SCRIPT`（必需）、`--submit-cmd`（`sbatch`；或 `qsub`） |
| `collect` | 解析输出 → CSV + 图 | `--natoms`（默认从输出的 `- Atoms:` 读取）、`--threshold`（0.03 meV/atom，约 1e-6 Ha/atom）、`--csv`、`-o/--outfile`、`--dpi`（200）、`--no-plot`、`--show` |

**模板要求**：模板中必须有两行预处理指令（数值任意）

```text
@SET CUTOFF 400
@SET REL_CUTOFF 60
```

并在 `&MGRID` 中写 `CUTOFF ${CUTOFF}`、`REL_CUTOFF ${REL_CUTOFF}`。变量名不同时用 `--var-cutoff/--var-rel` 指定。现成模板见 [examples/04-cutoff-convergence-si/template.inp](../examples/04-cutoff-convergence-si/template.inp)。模板还需满足：`PRINT_LEVEL MEDIUM` 及以上（否则 MULTIGRID INFO 里的 `count for grid` 不打印）、`EPS_SCF` 1e-8、`EPS_DEFAULT` 1e-12，扫描时关闭 .wfn 重启。

**用法示例**：

```bash
# step 1: CUTOFF scan at fixed REL_CUTOFF = 60 Ry, 4 MPI ranks, run locally then collect
python conv_cutoff.py run -t template.inp --scan cutoff --values 150 200 250 300 400 500 600 \
       --rel-cutoff 60 --extra si_bulk8.xyz --cp2k "mpirun -np 4 cp2k.psmp"

# step 2: REL_CUTOFF scan at the chosen CUTOFF, one batch job per point
python conv_cutoff.py gen -t template.inp --scan relcutoff --values 40 50 60 70 80 100 \
       --cutoff 400 --extra si_bulk8.xyz
python conv_cutoff.py submit --scan relcutoff --job ../../scripts/job_cp2k_slurm.sh
python conv_cutoff.py collect --scan relcutoff --threshold 0.03     # after all jobs finished

# rebuild CSV + figure from existing directories only
python conv_cutoff.py collect --scan cutoff --natoms 8 -o cutoff_scan.png
```

**输出**：

- `<scan>_scan.csv`：列为 `<scan>_Ry`（`cutoff_Ry` 或 `relcutoff_Ry`）、`energy_Ha`、`dE_meV_per_atom`、`scf_converged`、`pressure_bar`、`grid1_count`、`grid2_count`……
- `<scan>_scan.png`：|ΔE| 的对数坐标曲线；输出中有 `STRESS| 1/3 Trace` 时加压力面板，解析到 `count for grid` 时再加各级网格的 count 面板。
- 屏幕上打印表格，并给出“从哪个值起 |ΔE| < 阈值”。参考值是**扫描到的最大值**（最后一行）；如果判定收敛的值正好是倒数第二个点，说明扫描范围不够，应继续往上扩展。

**实现细节**：同时识别 2025.1+ 的 `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]` 和 2024.x 及更早的 `[a.u.]:` 格式。某个点 `count for grid 1` 为 0（最粗网格没有分到高斯函数，通常说明 REL_CUTOFF 过低）或 SCF 未收敛时会给出警告。`submit` 在每个目录内执行 `sbatch --export=ALL,INP=...,OUT=... <job>` 或 `qsub -v INP=...,OUT=... <job>`，本目录的两个作业模板正好读取 `INP`/`OUT` 变量。

**局限**：只比较总能量（及可选的压力），不检查力；`run` 的并行只是同时启动多个进程，请自行保证核数不超订。

---

## 3. conv_kpoints.py — k 点收敛测试

**用途**：按晶格自动给出一组 Monkhorst-Pack 网格（非周期方向固定为 1），把模板中 `&FORCE_EVAL/&DFT/&KPOINTS` 的 `SCHEME` 行替换为 `SCHEME MONKHORST-PACK n1 n2 n3`，生成 `kp_<n1>x<n2>x<n3>/` 目录。之后可本地运行或批量提交，最后汇总每原子能量、相对最密网格的能量差、费米能、压力和实际使用的 k 点数。对应 [第6章](../tutorials/06-periodic-surface-kpoints.md)。

**依赖**：标准库；作图需 `numpy matplotlib`；只有 `--structure` 需要 `ase`。

**子命令与参数**：

| 子命令 | 参数 |
|---|---|
| 网格选项（`mesh`/`gen`/`run`） | `-t/--template`（必需）；三选一：`--lengths R ...`（Å，默认 `10 15 20 25 30 40 50`，n_i = max(1, ceil(R·\|b_i\|))，\|b_i\| 不含 2π）、`--kspacing K ...`（1/Å，含 2π，n_i = max(1, ceil(2π\|b_i\|/K))，与 VASP KSPACING 同义）、`--grids 4x4x4 ...`（显式网格）；`--cell a b c [alpha beta gamma]`、`--structure FILE`（用 ase 读晶胞）、`--nonperiodic AXES`（如 `Z`，默认按 `&CELL PERIODIC` 判断）、`--even-only`（把大于 1 的奇数 n 向上取为偶数） |
| `mesh` | 只打印将要生成的网格 |
| `gen` | 网格选项 + `--extra FILES`；公共 `--inp-name`、`--out-name`、`--force` |
| `run` | gen 的全部参数 + `--cp2k`（`"cp2k.psmp"`）、`--npar`（1）+ collect 参数 |
| `submit` | `--job SCRIPT`（必需）、`--submit-cmd`（`sbatch` 或 `qsub`）；公共 `--inp-name`、`--out-name`、`--force` |
| `collect` | `--natoms`、`--threshold`（1.0 meV/atom）、`--csv`（`kpoint_scan.csv`）、`-o/--outfile`（`kpoint_scan.png`）、`--dpi`（200）、`--no-plot`、`--show` |

**模板要求**：

- `&FORCE_EVAL/&DFT` 内有 `&KPOINTS` 段和一行 `SCHEME`。该行整行替换，原来写成 `${NK}` 之类的变量也可以，例如 [examples/06-periodic-surface-kpoints/kpoint_convergence_template.inp](../examples/06-periodic-surface-kpoints/kpoint_convergence_template.inp)。
- 晶胞从 `&SUBSYS/&CELL` 读取（`ABC` [+ `ALPHA_BETA_GAMMA`] 或 `A/B/C` 矢量，支持 `[angstrom]`/`[bohr]` 单位）。晶胞写成 `${变量}` 或用 `CELL_FILE_NAME` 时，请给 `--cell` 或 `--structure`。
- CP2K 2026.2 及更早版本中，OT 不能与 k 点同时使用（只有 trunk 支持），模板必须用 `&DIAGONALIZATION` + `&MIXING`；金属还需要 `&SMEAR` 和 `ADDED_MOS`。
- 需要压力列时，打开 `STRESS_TENSOR ANALYTICAL` 和 `&FORCE_EVAL/&PRINT/&STRESS_TENSOR`。

**用法示例**：

```bash
python conv_kpoints.py mesh -t si_primitive_kpoints.inp --lengths 10 20 30 40 50     # preview only
python conv_kpoints.py run -t kpoint_template.inp --lengths 10 15 20 25 30 40 \
       --cp2k "mpirun -np 8 cp2k.psmp"
python conv_kpoints.py gen -t slab.inp --kspacing 0.5 0.4 0.3 0.25 0.2 --extra slab.xyz
python conv_kpoints.py submit --job ../../scripts/job_cp2k_slurm.sh
python conv_kpoints.py collect --threshold 1.0
python conv_kpoints.py gen -t t.inp --grids 2x2x2 4x4x4 6x6x6 8x8x8 --structure POSCAR
```

**输出**：`kpoint_scan.csv` 的列为 `mesh, n1, n2, n3, nkpts, energy_Ha, energy_per_atom_eV, dE_meV_per_atom, fermi_eV, pressure_bar, scf_converged`。`nkpts` 来自 `BRILLOUIN| List of Kpoints`，`fermi_eV` 只在使用 `&SMEAR` 时有值。网格按 n1·n2·n3 排序，参考值是**最密的网格**。图中显示 |ΔE| 随网格的变化。

**局限**：ΔE 是相对最密网格而言的，本身不能证明最密网格已经收敛，所以最后一两个点之间的差值最好也远小于阈值。只做 Monkhorst-Pack，不支持 Gamma-centred 偏移和显式 k 点列表。

---

## 4. cp2k_outparse.py — 输出文件汇总

**用途**：一次性汇总 CP2K 主输出的关键信息，适合计算结束后快速体检和批量检查。汇总内容包括：是否正常结束（`PROGRAM ENDED AT`）、SCF 是否收敛与步数、内外层 SCF、总能量（含展宽时的 free energy）、费米能、HOMO-LUMO/带隙、Mulliken 电荷与自旋、WARNING 去重汇总、`[ABORT]` 内容、峰值内存与总耗时；GEO_OPT/CELL_OPT 还会列出每一步 `OPT|` 块中的 Max/RMS step、Max/RMS grad、内压与压力偏差。

**依赖**：仅标准库（Python ≥ 3.9）。

**参数**：`outputs ...`（一个或多个输出文件，支持 shell 通配）、`--json`（输出 JSON）、`--brief`（每个文件一行）、`--opt`（打印全部 OPT| 步骤表）、`--opt-csv FILE`（OPT| 步骤写入 CSV）、`--mulliken`（打印最后一个 Mulliken 表）、`--warnings N`（列出的不同 WARNING 条数，默认 10）。

```bash
python cp2k_outparse.py h2o_energy.out                  # human-readable summary
python cp2k_outparse.py h2o_geo_opt_bfgs.out --opt      # + table of every OPT| step
python cp2k_outparse.py fe_afm.out --mulliken           # + Mulliken charges / spin moments
python cp2k_outparse.py run.out --warnings 20
python cp2k_outparse.py */output.out --brief            # one line per file (batch check)
python cp2k_outparse.py si_cell_opt.out --json > si.json
python cp2k_outparse.py si_cell_opt.out --opt-csv opt_steps.csv
python cp2k_outparse.py x.out --brief || echo "check x.out"
```

**退出码**：所有文件都正常结束且没有 ABORT 时返回 0，否则返回 1，可直接用在 shell 循环和作业脚本里。

**解析的输出行**（按 2026.2 源码核对，2023.x–2025.x 相同，除非注明）：`CP2K| version string`、`GLOBAL| Project name / Run type / Total number of message passing processes / Number of threads`、`- Atoms:`、`*** SCF run converged in N steps ***`、`outer SCF loop converged / FAILED to converge`、`SCF run NOT converged`；`ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]`（2025.1+）与 `[a.u.]:`（2024.x 及更早），单位从行内读取，因此改过 `ENERGY_UNIT` 也能正确识别；`Fermi energy:`（hartree，展宽时打印）、`MO| E(Fermi):`、`MO| Band gap:`、`HOMO - LUMO gap [eV] :`；`*** WARNING in <file>:<line> :: ... ***`（多行消息会拼接）；`MEMORY| Estimated peak process memory [MiB]`；计时表 `CP2K` 行；`GEOMETRY OPTIMIZATION COMPLETED`、`MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED`、`STRESS| 1/3 Trace`、`CELL| Volume`。

**局限**：`OPT|` 块是 2023.1 起的格式，更早的 `Informations at step` 格式只能尽力解析；不解析 `.ener`、`.restart` 等附属文件（分别用 [`md_analyze.py`](md_analyze.py)、[`restart_update.py`](restart_update.py)）。

---

## 5. plot_opt.py — 优化收敛曲线

**用途**：读取 GEO_OPT / CELL_OPT 输出中每一步的 `OPT|` 块，画出能量（相对最低值）、最大/RMS 梯度、最大/RMS 位移、压力偏差随步数的变化，并用虚线标出收敛阈值（实心点表示该判据已满足）。如果找到 `<PROJECT>-1.cell`（`&MOTION/&PRINT/&CELL` 打开时生成），会再加一幅晶胞体积曲线。优化还没结束时也可以用，便于判断是否卡住。对应 [第7章](../tutorials/07-geo-cell-opt.md)。

**依赖**：`pip install numpy matplotlib`。

**参数**：`output`（GEO_OPT/CELL_OPT 主输出）、`--cell FILE`（指定 .cell 文件，默认在输出旁自动查找）、`--no-cell`、`--energy-unit hartree|eV|meV|kcal/mol|kJ/mol`（默认 meV）、`-o/--outfile`（默认 `<输出文件名>_opt.png`）、`--dpi`（150）、`--csv FILE`、`--no-plot`、`--show`。

```bash
python plot_opt.py h2o_geo_opt_bfgs.out                      # -> h2o_geo_opt_bfgs_opt.png
python plot_opt.py si_cell_opt.out                           # finds si_cell-1.cell automatically
python plot_opt.py si_cell_opt.out --cell si_cell-1.cell -o si_cell.png --dpi 200
python plot_opt.py opt.out --energy-unit kcal/mol --csv opt_steps.csv
```

**说明**：梯度和位移保持 CP2K 原始单位（hartree/bohr、bohr），可直接与 `MAX_FORCE`、`MAX_DR` 等阈值比较。第 0 步只有初始能量。`.cell` 文件的列为 `Step Time[fs] Ax Ay Az Bx By Bz Cx Cy Cz Volume[Angstrom^3]`。旧版（≤ 2022）`Informations at step` 格式只能尽力解析，图中不画阈值线。

---

## 6. restart_update.py — 把最新结构写回输入

**用途**：从 `<PROJECT>-1.restart`（或 `-pos-1.xyz` 轨迹的最后一帧/指定帧，可配合 `-1.cell`）读取最新坐标与晶胞，然后：

1. 写成 xyz 或带 `Lattice=` 的扩展 xyz，用于可视化或后续计算；
2. 写回一个已有的输入：替换 `&SUBSYS` 中的 `&COORD` 和 `&CELL` 的 A/B/C。如果输入用的是 `&TOPOLOGY/COORD_FILE_NAME`，则写出新的 xyz 并修改文件名。其他设置原样保留。

适合优化中断后续算，或把优化好的结构交给能带、振动等后续计算。对应 [第5章](../tutorials/05-scf-strategy.md)、[第7章](../tutorials/07-geo-cell-opt.md)、[第12章](../tutorials/12-aimd.md)。

**依赖**：仅标准库。

**参数**：`source`（`.restart` 或 xyz/扩展 xyz 轨迹）、`-i/--input INP`（要更新的输入，原文件不改动）、`-o/--output NEW_INP`（默认 `<输入名>_updated.inp`）、`--xyz FILE`（写出 xyz；与 `-i` 和 `COORD_FILE_NAME` 同用时，这就是新输入引用的文件，默认 `<新输入名>.xyz`）、`--extxyz`、`--frame N`（0 起算，默认 -1 = 最后一个完整帧）、`--cell-file CELL`（为轨迹帧提供晶胞，按相同步号匹配，找不到则用最后一行）、`--no-cell`（不改输入的晶胞）、`--ext-restart`（追加 `&EXT_RESTART` / `RESTART_FILE_NAME <source>`，source 必须是 .restart）、`--force`（允许原子数不同、覆盖文件、非 XYZ 的 `COORD_FILE_FORMAT`）。

```bash
python restart_update.py h2o_opt-1.restart --xyz h2o_opt_final.xyz
python restart_update.py si_cell-1.restart -i si_band.inp -o si_band_opt.inp
python restart_update.py si_cell-pos-1.xyz --cell-file si_cell-1.cell -i si.inp -o si_new.inp
python restart_update.py md-pos-1.xyz --frame 10 --cell-file md-1.cell --xyz frame10.xyz --extxyz
python restart_update.py h2o_opt-1.restart -i h2o_opt.inp -o h2o_opt_cont.inp --ext-restart
```

**说明**：

- 坐标以 Å（笛卡尔）写出，原子标签保持 restart/轨迹中的 kind 名（如 AFM 的 `Fe1`/`Fe2`），`&COORD` 行多出的分子名、残基号等列也会保留。
- `&CELL` 顶层的 `ABC`、`ALPHA_BETA_GAMMA`、`A/B/C`、`CELL_FILE_NAME`、`CELL_FILE_FORMAT` 会被替换为 A/B/C 矢量（Å）；`PERIODIC`、`SYMMETRY`、`MULTIPLE_UNIT_CELL` 和 `&CELL_REF` 保持不变。
- `.restart` 里除结构外还有速度、恒温器状态、优化器历史、波函数文件名等。要**原样续算同一个任务**，最干净的做法是 `&EXT_RESTART`（`--ext-restart`）；本脚本主要用于把结构搬到新的或修改过的输入里。
- 普通 xyz 轨迹不含晶胞，CELL_OPT/NPT 请配合 `--cell-file` 或使用扩展 xyz。
- 用了 `SPLIT_RESTART_FILE` 时，`.restart` 文本里没有坐标，请改用 `-pos-1.xyz` 轨迹。

---

## 7. band_plot.py — 能带图与带隙

**用途**：解析 `&DFT/&PRINT/&BAND_STRUCTURE` 写出的能带文件（`FILE_NAME` 指定的 .bs，内含 `# Set / Special point / Point / Band` 块），按各 `KPOINT_SET` 生成横轴与高对称点标签，以 VBM 或费米能为零点画图，判断直接/间接带隙，可导出 CSV。对应 [第10章](../tutorials/10-dos-pdos-band.md)，示例输入见 [examples/10-dos-pdos-band/si_band_structure.inp](../examples/10-dos-pdos-band/si_band_structure.inp)。

**依赖**：`pip install numpy matplotlib`。

**参数**：`bsfile`（.bs 文件，或未设 `FILE_NAME` 时含能带数据的主输出）、`--cell-from OUT`（从主输出的 `CELL| Vector a/b/c [angstrom]` 读晶胞）、`--lattice 9 个数`（a1 a2 a3 b1 b2 b3 c1 c2 c3，Å）、`--zero vbm|fermi|none`（默认 vbm）、`--efermi EV`（绝对费米能，配合 `--zero fermi`）、`--fermi-from OUT`（读最后一个 `Fermi energy:`）、`--nocc N`（每个自旋通道的占据带数，代替按占据数判断）、`--erange EMIN EMAX`（相对零点，默认 -6 6 eV）、`--labels L ...`（覆盖高对称点标签，断点写 `U|K`）、`--csv FILE`、`-o/--output`（`band.png`，`none` 不出图）、`--title`、`--dpi`（200）、`--show`。

```bash
python band_plot.py si_band.bs --cell-from si_band_structure.out --csv si_band.csv
python band_plot.py si_band.bs --lattice 0 2.7155 2.7155 2.7155 0 2.7155 2.7155 2.7155 0 --erange -13 6
python band_plot.py graphene.bs --cell-from graphene.out --zero fermi --fermi-from graphene.out
python band_plot.py si_band_structure.out --nocc 4 -o si_band.png
```

**注意**：

- CP2K 以**追加模式**写 .bs 文件。同一文件里有多次运行（重复出现 `# Set 1`）时，只使用最后一次。
- k 点以倒格矢分数坐标给出，必须提供晶胞（`--cell-from` 或 `--lattice`）才能得到真实的笛卡尔路径长度（1/Å）。没有晶胞时脚本按分数坐标计算并给出警告，非立方晶胞的能带图横轴会变形。
- 相邻两个 KPOINT_SET 首尾点不同时，标签合并为 `U|K`，表示路径断开。
- `--zero vbm` 时以占据数 > 最大占据数一半的最高能级为 VBM；金属请用 `--zero fermi`。CSV 写出 x、kx、ky、kz 和各能带（按自旋分列），刻度位置写到 `FILE_ticks.csv`。

---

## 8. pdos_plot.py — DOS/PDOS 作图

**用途**：读取 CP2K 的 .pdos 文件，展宽、对齐能量零点、按元素/角动量求和，分自旋作图并导出 CSV。支持两种格式（同一张图里不能混用）：

1. **逐态原始文件**（未展宽）：`<PROJECT>-k1-1.pdos`、`<PROJECT>-ALPHA_k1-1.pdos`、`<PROJECT>-list1-1.pdos` 等。CP2K ≤ 2026.1 的 `&PRINT/&PDOS` 把 `E(Fermi)` 写在首行；2026.2 的 `&PRINT/&DOS/&PDOS` 把 `E(Fermi)`、`E(HOCO)` 单独成行，并在 `Occupation` 后多一列 `Total`。脚本按表头列名定位各列，不依赖固定列号。能量（hartree）换算为 eV 后再做高斯或洛伦兹展宽。
2. **2026.2 的 `&CURVE` 展宽曲线**：`<PROJECT>-k1_curve-1.pdos` 等，CP2K 已按 `ENERGY_ZERO` 平移并展宽，默认原样绘制。

**依赖**：`pip install numpy matplotlib`。

**参数**：`files ...`（文件或通配；不给时读 `-d` 目录下所有 `*.pdos`）、`-d/--dir`（`.`）、`--mode auto|raw|curve`（两类文件都存在时选哪类；auto：有原始文件就用原始文件）、`--group total|element|l|element-l|orbital`（默认 element）、`--by-kind`（保留 `Fe1`/`Fe2` 等 kind 名，不合并为元素）、`--elements EL ...`、`--zero fermi|hoco|none|file`（默认：原始文件 fermi，曲线文件 file）、`--efermi EV`（绝对参考能量，优先于 `--zero`）、`--fwhm`（原始文件展宽 FWHM，默认 0.2 eV）、`--lorentzian`、`--de`（能量网格步长，0.01 eV）、`--erange EMIN EMAX`（默认 -10 5）、`--spin both|up|down|sum`（默认 both，beta 画在零线以下）、`--total`（另画所有已载入文件之和）、`--scale`、`--fill`、`-o/--output`（`pdos.png`）、`--csv FILE`、`--title`、`--ymax`、`--dpi`（200）、`--show`。

```bash
python pdos_plot.py                                         # all *.pdos in the current directory
python pdos_plot.py tio2_pdos-k*-1.pdos --group element-l --erange -8 6
python pdos_plot.py *-k*-1.pdos --zero hoco --fwhm 0.3 --lorentzian -o h2o_pdos.png
python pdos_plot.py h2o_pdos-k*_curve-1.pdos --group orbital --elements O
python pdos_plot.py fe-ALPHA_k1-1.pdos fe-BETA_k1-1.pdos --group l --csv fe_pdos.csv
python pdos_plot.py *.pdos --group total --total --show
```

**注意**：`--zero fermi` 取各自旋文件头中 E(Fermi) 的最大值（与 CP2K 一致）；`hoco` 优先用 2026.2 的 `E(HOCO)` 行，旧文件则取占据数 > 最大占据数一半的最高态。文件中有多个追加块（如 GEO_OPT 中按 `&EACH` 输出）时，只用最后一块。两种 PDOS 输入写法的对照见 [第10章](../tutorials/10-dos-pdos-band.md)，示例见 [examples/10-dos-pdos-band/README.md](../examples/10-dos-pdos-band/README.md)。

---

## 9. cube_tools.py — cube 文件工具箱

**用途**：处理 CP2K 写出的 Gaussian 格式 cube 文件（电子密度、V_HARTREE、MO、差分密度等）。对应 [第11章](../tutorials/11-cube-charge-workfunction.md)，第6章的表面偶极校正也用它检查真空电势平台。

**依赖**：`pip install numpy`；`planar` 作图需 matplotlib。

**子命令与参数**：

| 子命令 | 作用 | 参数 |
|---|---|---|
| `info CUBE ...` | 头信息、网格、晶胞、数值范围、积分 | — |
| `check CUBE ...` | 格式检查；多个文件时检查网格是否一致（做差分的前提），有问题返回 1 | — |
| `integrate CUBE` | 积分 = Σ值 × 体素体积；对 ELECTRON_DENSITY 即（价）电子数，另给正、负部分 | `--axis a\|b\|c\|x\|y\|z`（c）、`--range Z1 Z2`（只积分该区间的平面，Å，从 cube 原点量起）、`--abs` |
| `planar CUBE` | 沿某轴的平面平均。V_HARTREE（hartree→eV）配合费米能时给出真空能级与功函数 Φ = V_vac − E_F；密度/差分密度给出单位长度电荷与累积 ΔQ(z) | `--axis`（c）、`--kind auto\|potential\|density\|other`、`--fermi EV` 或 `--out CP2K.out`（二选一，读最后一个 `Fermi energy:`）、`--window Z1 Z2`（手动真空窗口，Å）、`--vac-offset`（4 Å）、`--vac-width`（2 Å）、`--csv`、`-o/--output`（`planar.png`）、`--dpi`、`--show` |
| `diff OUT A B [C ...]` | OUT = A − B − C …（或 `--coeffs` 指定的线性组合），要求同一网格 | `--coeffs W ...`（默认 1 -1 -1 …） |

```bash
python cube_tools.py info h2o_cube-ELECTRON_DENSITY-1_0.cube
python cube_tools.py integrate h2o_cube-ELECTRON_DENSITY-1_0.cube
python cube_tools.py planar pt111_wf-v_hartree-1_0.cube --out pt111_workfunction.out -o vh.png --csv vh.csv
python cube_tools.py planar pt111_wf-v_hartree-1_0.cube --fermi -5.21 --window 20 22
python cube_tools.py check co_cu-ELECTRON_DENSITY-1_0.cube cu-ELECTRON_DENSITY-1_0.cube co-ELECTRON_DENSITY-1_0.cube
python cube_tools.py diff drho.cube co_cu-ELECTRON_DENSITY-1_0.cube cu-ELECTRON_DENSITY-1_0.cube co-ELECTRON_DENSITY-1_0.cube
python cube_tools.py planar drho.cube --kind density -o drho_z.png
```

**注意**：

- CP2K 的 V_HARTREE cube 单位是 hartree，表示电子的静电势能（与本征值同号），因此真空平台位于 E_F 之上，Φ = V_vac − E_F > 0。
- 自动真空窗口的取法：沿该轴原子间最大的空隙视为真空，在两侧最外层原子外 `--vac-offset` 处各取宽 `--vac-width` 的窗口求平均。两侧数值不同说明存在偶极，参见 [第6章](../tutorials/06-periodic-surface-kpoints.md) 的偶极校正。
- `Fermi energy:` 只在使用展宽（`&SMEAR`）时才会打印；绝缘体请用 `--fermi` 手动给出参考能量，并在文中说明取的是哪一个。
- 差分电荷的三个 cube 必须使用相同的晶胞、网格和原子位置，先用 `check` 确认。示例见 [examples/11-cube-charge-workfunction/README.md](../examples/11-cube-charge-workfunction/README.md)。

---

## 10. md_analyze.py — AIMD 后处理

**用途**：处理 `&MOTION/&MD` 的各类输出，无需 ase 或 MDAnalysis。对应 [第12章](../tutorials/12-aimd.md)，[第15章](../tutorials/15-metadynamics-plumed.md)（元动力学）和[第16章](../tutorials/16-qmmm-xtb.md)（QM/MM、xTB MD）也使用。

**依赖**：`pip install numpy matplotlib`。

**子命令**：

| 子命令 | 输入 | 输出与主要参数 |
|---|---|---|
| `ener` | `<PROJECT>-1.ener`（列：Step Nr.、Time[fs]、Kin.[a.u.]、Temp[K]、Pot.[a.u.]、Cons Qty[a.u.]、UsedTime[s]） | T 与 E_pot 的均值/标准差、守恒量漂移（线性拟合与首末差，给 `--natoms` 时按每原子计）；`--natoms`、`--skip-fs`；默认图 `ener.png` |
| `rdf` | `<PROJECT>-pos-1.xyz` + 晶胞 | 偏径向分布函数 g(r)（最小镜像，支持三斜和 NPT）与配位数 n(r)，报告第一峰与第一极小；`--pairs A-B ...`（O-O）、`--rmax`（默认最小晶胞宽度的一半）、`--dr`（0.02 Å）；默认图 `rdf.png` |
| `msd` | `<PROJECT>-pos-1.xyz` + 晶胞 | 按元素的 MSD（多时间原点，先展开跨边界跳跃并扣除质心漂移），由 MSD = 6Dt + c 拟合 D；`--elements`、`--max-lag-fs`、`--origin-stride`（1）、`--fit T1 T2`（ps，默认最大延迟的 20%–80%）、`--keep-com`；默认图 `msd.png` |
| `vacf` | `<PROJECT>-vel-1.xyz`（或位置轨迹 + `--from-positions`） | VACF（默认质量加权，`--no-mass` 关闭）与振动态密度（Hann 窗 FFT）；`--elements`、`--max-lag-fs`、`--max-freq`（4500 cm⁻¹）、`--pad`（4）；默认图 `vdos.png` |
| `ir` | `&DFT/&PRINT/&MOMENTS` 写出的偶极文件（或含该数据的主输出） | 偶极导数自相关得到的 IR 线形，Berry 相偶极按打印的晶胞量子展开；`--dt-fs`（必需，= MD TIMESTEP × `&EACH/MD`）、`--skip-fs`、`--max-lag-fs`、`--max-freq`、`--pad`、`--qcf none\|harmonic`（量子校正因子）、`--temperature`（300 K）；默认图 `ir.png` |

轨迹类子命令（rdf/msd/vacf）的公共参数：`--stride N`（每 N 帧取一帧）、`--skip-fs T`（丢弃前 T fs 平衡段）、`--dt-fs`（xyz 标题行没有 `time =` 时指定帧间隔）。晶胞三选一：`--cell <PROJECT>-1.cell`（逐帧晶胞，NPT）、`--abc A B C [--angles α β γ]`（固定晶胞，角度默认 90 90 90）、`--lattice 9 个数`。所有子命令都有 `-o/--output`、`--csv`、`--dpi`（200）、`--show`。

```bash
python md_analyze.py ener water32_nve-1.ener --natoms 96 -o nve_ener.png
python md_analyze.py rdf water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --pairs O-O O-H --skip-fs 1000
python md_analyze.py rdf water32_npt-pos-1.xyz --cell water32_npt-1.cell --pairs O-O
python md_analyze.py msd water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --elements O --fit 1 4
python md_analyze.py vacf water32_nve-vel-1.xyz --elements H O --max-lag-fs 1000
python md_analyze.py ir water32_nve_dipole.dat --dt-fs 1.0 --temperature 300 --qcf harmonic
```

**局限**：假定帧间隔相等。32 个水分子、几个 ps 的轨迹统计噪声很大，扩散系数还需要有限尺寸校正（Yeh & Hummer, *J. Phys. Chem. B* 108, 15873 (2004), https://doi.org/10.1021/jp0477147）。示例体系见 [examples/12-aimd-water/README.md](../examples/12-aimd-water/README.md)。

---

## 11. neb_plot.py — NEB 能量曲线

**用途**：收集 `RUN_TYPE BAND`（NEB / CI-NEB）各副本的能量，以累计位移为反应坐标画最小能量路径（MEP），标注正/逆反应能垒和反应能，并导出能量最高的中间副本（过渡态近似）为 xyz，供 Dimer 或振动分析继续使用。对应 [第13章](../tutorials/13-neb-dimer.md)，示例见 [examples/13-neb-dimer-h-cu111/README.md](../examples/13-neb-dimer-h-cu111/README.md)。

**依赖**：`pip install numpy`；作图需 matplotlib。

**可识别的输入**（自动判断）：

- 主输出中 `&BAND/&PROGRAM_RUN_INFO` 的汇总块（`DISTANCES REP = ...`，bohr；`ENERGIES [au] = ...`，hartree），或 2026.2 起 `PLOT_REL_ENERGY T` 时的 `REPLICA  ENERGY [au] ...` 表格。该打印键在 2026.2 中默认级别为 LOW，更早版本为 MEDIUM。
- `<PROJECT>-*.ener`（`&BAND/&ENERGY` 打印键），每个 band 步一行。
- `<PROJECT>-FINAL-1.xyz`（2026.2 起），标题行含 `Replica=`、`Energy=`、`Ener_loc_max=` 的扩展 xyz。
- 全部副本轨迹 `<PROJECT>-pos-Replica_nr_*-1.xyz`，取每个文件的最后一帧，能量取自标题行 `E = ...`。

**参数**：`inputs ...`（支持通配）、`--step N`（从 .out/.ener 中取第几个 band 步，1 起算，负数从末尾数，默认 -1）、`--xyz FILE ...`（能量来自 .out/.ener 时用于 `--ts-xyz` 的结构）、`--ts-xyz FILE`、`--x dist|index`（默认累计路径长度，Å）、`--mic`（副本间距用最小镜像，需要 FINAL-1.xyz 中的 `Lattice=`）、`--ref first|min`（能量零点：第一个副本或最低副本）、`--units eV|kJ/mol|kcal/mol`（默认 eV）、`--spline`（自然三次样条）、`--csv FILE`、`-o/--output`（`mep.png`）、`--title`、`--dpi`、`--show`。

```bash
python neb_plot.py h_diffusion_cineb.out -o mep.png
python neb_plot.py h_cu111_neb-1.ener --step -1 --units kJ/mol --csv mep.csv
python neb_plot.py h_cu111_neb-FINAL-1.xyz --ts-xyz ts_from_neb.xyz --spline
python neb_plot.py h_cu111_neb-pos-Replica_nr_*-1.xyz --ts-xyz ts_guess.xyz
python neb_plot.py h_diffusion_cineb.out --xyz "h_cu111_neb-pos-Replica_nr_*-1.xyz" --ts-xyz ts.xyz
```

**注意**：每个副本组的日志 `<PROJECT>-BAND*.out` 对应的是副本**组**（一个组可含多个副本），只用于排查错误，不要拿来读能量。样条只用于显示，能垒数值按离散副本计算。最高能副本只是过渡态的近似，应按第13–14章用 Dimer 或振动分析（恰好一个虚频）确认。

---

## 12. vib_thermo.py — 振动频率与热化学

**用途**：解析 `RUN_TYPE VIBRATIONAL_ANALYSIS` 的主输出和/或 `<PROJECT>-VIBRATIONS-1.mol`（molden），列出频率、IR 强度与约化质量，检查虚频，计算 ZPE 及谐振子近似下的 U、H、S、G。对应 [第14章](../tutorials/14-vibration-thermo.md)，示例见 [examples/14-vibration-thermo/README.md](../examples/14-vibration-thermo/README.md)。

**依赖**：`pip install numpy`。

**输入**：主输出中**最后一个** `NORMAL MODES - CARTESIAN DISPLACEMENTS` 表（`VIB|Frequency (cm^-1)`，负值表示虚频；`VIB|IR int (KM/Mole)` 只在 `INTENSITIES T` 时出现；`VIB|Red.Masses (a.u.)`）、电子能量（`Minimum Structure - Energy and Forces` 后的 `Total energy:`）以及 CP2K 自带的热化学块（用于对比）。molden 文件提供 `[FREQ]`、`[FR-COORD]`（bohr）和 `[INT]`，其中 `[INT]` 不是 km/mol，只能作相对强度。

**模型**：

- `--mode adsorbate`（默认）：全部模式按振动处理，没有平动和转动，适用于吸附物和部分 Hessian（`&MODE_SELECTIVE`、固定原子）。
- `--mode gas`：理想气体 RRHO 分子，加上平动（压力 p 下的 Sackur-Tetrode）、刚性转子（对称数 `--sigma`，线性与否由转动惯量判断）和电子熵 R ln(2S+1)（`--mult`）；H = U + RT，G = H − TS。需要几何结构（molden 或 `--xyz`）。
- 虚频会列出，且**始终不计入**求和。CP2K 自带的热化学会把虚频算进去，并固定使用对称数 1，所以两者结果可能不同。`--floor F` 把低于 F cm⁻¹ 的实频抬到 F；`--qrrho` 使用 Grimme quasi-RRHO 振动熵（*Chem. Eur. J.* 18, 9955 (2012), https://doi.org/10.1002/chem.201200497；w0 = `--qrrho-w0`，B_av = 1e-44 kg m²）。

**参数**：`files ...`、`-T/--temperature K ...`（可多个，默认 298.15）、`-p/--pressure`（Pa，101325，用于 gas）、`--mode`、`--sigma`（1）、`--mult`（1）、`--linear auto|yes|no`、`--xyz FILE`、`--mass EL=M ...`（同位素质量，如 `H=2.014`）、`--floor CM`（0 = 关闭）、`--qrrho`、`--qrrho-w0`（100 cm⁻¹）、`--soft`（低于此值的实频标记为软模，50 cm⁻¹）、`--ts`（过渡态：应恰有一个虚频）、`--e0 HA`（电子能量，默认读输出）、`--csv FILE`。

```bash
python vib_thermo.py h2o_vib.out h2o_vib-VIBRATIONS-1.mol --mode gas --sigma 2 -T 298.15
python vib_thermo.py co_on_pt_partial_hessian.out -T 300 500 --floor 50
python vib_thermo.py ts_freq_check.out --ts
python vib_thermo.py co_pt111_vib-VIBRATIONS-1.mol --qrrho -T 298.15 --csv modes.csv
```

**局限**：只做谐振子/RRHO 近似，不处理受阻转子与非谐性；频率的可信度取决于振动分析前的优化精度和 `DX` 步长，见第14章。

---

## 13. job_cp2k_slurm.sh / job_cp2k_pbs.sh — 作业模板

**用途**：提交 `cp2k.psmp`（MPI × OpenMP）的作业模板。脚本会做以下事情：加载环境（文件中有可编辑的环境块，toolchain 的 `cp2k_env`、`module load`、conda 三选一）；设置 `OMP_NUM_THREADS`/`OMP_PLACES`/`OMP_PROC_BIND` 和 `ulimit -s unlimited`；已有的输出文件改名为 `<out>.<n>.bak`，不会被覆盖；运行结束后检查 `PROGRAM ENDED AT`、`[ABORT]`、`SCF run NOT converged`、优化是否收敛、WARNING 条数，异常时打印 `&EXT_RESTART` 续算片段和 [`restart_update.py`](restart_update.py) 的命令。正常结束退出码为 0，否则为 1。

**依赖**：bash；分别需要 Slurm 或 PBS Pro / OpenPBS / Torque。使用前先编辑头部的 `#SBATCH` / `#PBS` 资源行和环境块。

| 变量 | Slurm 模板 | PBS 模板 |
|---|---|---|
| `INP` | 第一个命令行参数，或环境变量 `INP`（默认 `input.inp`） | 只能用 `qsub -v INP=...`（PBS 不传脚本参数） |
| `OUT` | 默认 `${INP%.*}.out` | 同左 |
| `CP2K_EXE` | `cp2k.psmp` | `cp2k.psmp` |
| 启动方式 | `LAUNCHER=srun`（默认，`srun --cpu-bind=cores`）或 `mpirun` | `MPIRUN=mpirun`，`MPI_EXTRA` 附加选项 |
| 线程数 | 取 `--cpus-per-task`，并设置 `SRUN_CPUS_PER_TASK`（Slurm ≥ 22.05 的 srun 不再继承该值） | `NTHREADS`（默认取 `OMP_NUM_THREADS`，否则 1） |
| MPI 进程数 | Slurm 分配的 ntasks | `$PBS_NODEFILE` 行数；Torque `nodes=N:ppn=P` 时设 `NODEFILE_LINES=cores`，进程数 = 行数 / NTHREADS |

```bash
sbatch job_cp2k_slurm.sh my.inp
sbatch --export=ALL,INP=my.inp,OUT=my.out job_cp2k_slurm.sh
sbatch -N 2 --ntasks-per-node=16 --cpus-per-task=4 job_cp2k_slurm.sh my.inp
qsub -v INP=my.inp job_cp2k_pbs.sh
qsub -v INP=my.inp,OUT=my.out,NTHREADS=2 job_cp2k_pbs.sh
bash job_cp2k_slurm.sh -h                  # print the usage header
```

> [!TIP]
> `conv_cutoff.py submit` 和 `conv_kpoints.py submit` 会通过 `--export=ALL,INP=...,OUT=...`（sbatch）或 `-v INP=...,OUT=...`（qsub）把文件名传给这两个模板，可直接用作 `--job`。NEB 的总进程数必须等于 `NUMBER_OF_REPLICA × NPROC_REP`。

MPI × OpenMP 的搭配建议（杂化泛函、大内存任务适合少进程、多线程）见 [第9章](../tutorials/09-hybrid-admm.md) 和 [第18章](../tutorials/18-troubleshooting.md)，调度系统基础见 [Common/tutorials/02-job-scheduler.md](../../Common/tutorials/02-job-scheduler.md)。

---

## 典型工作流

```bash
# 1) draft an input, check the syntax
python cp2k_inpgen.py Si.cif --kpoints 4 4 4 --smear 300 -o si.inp
cp2k.psmp -c si.inp
# 2) converge CUTOFF / REL_CUTOFF and the k-mesh (templates derived from si.inp)
python conv_cutoff.py run -t template.inp --scan cutoff --values 200 300 400 500 600 --rel-cutoff 60
python conv_kpoints.py run -t kpoint_template.inp --lengths 10 20 30 40
# 3) production run + health check
sbatch job_cp2k_slurm.sh si_cell_opt.inp
python cp2k_outparse.py si_cell_opt.out --opt
python plot_opt.py si_cell_opt.out
# 4) move the optimized structure into the band / DOS input and plot
python restart_update.py si_cell-1.restart -i si_band.inp -o si_band_opt.inp
python band_plot.py si_band.bs --cell-from si_band_opt.out
python pdos_plot.py *-k*-1.pdos --group element-l
```

## 相关资源

- CP2K 官方手册（输入参考）：https://manual.cp2k.org/
- CP2K 官方教程与练习：https://www.cp2k.org/howto
- 其他社区工具（cp2k-input-tools、cp2kdata、ASE 的 CP2K calculator、Multiwfn 等）的介绍见 [第17章](../tutorials/17-tools-ecosystem.md)，自动化示例见 [examples/17-automation-tools/README.md](../examples/17-automation-tools/README.md)。

> [!CAUTION]
> 脚本按 CP2K 输出的真实格式编写，但没有覆盖所有打印级别和版本组合。用在新体系或新版本上时，先拿一个小算例手工 `grep` 对照脚本输出；如果发现解析错误，欢迎在 GitHub Issues 中附上输出片段反馈。

[返回目录](../README.md)
