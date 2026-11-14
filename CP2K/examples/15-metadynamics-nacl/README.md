# 示例 15：NaCl 在水中解离的自由能（元动力学、PLUMED、blue moon 约束 MD）

配套教程：[第15章 增强采样：元动力学、PLUMED 与约束 MD](../../tutorials/15-metadynamics-plumed.md)

## 目的

以 Na⁺Cl⁻ + 62 H2O（12.42 Å 立方盒子）为例，用三种方法计算离子对解离的自由能曲线 $F(d_\mathrm{Na-Cl})$，并演示“xTB 预采样 → DFT 精修”的流程：

| 方法 | 输入 | 说明 |
|---|---|---|
| 平衡 + CV 涨落 | `nacl_equil.inp` → `cv_monitor.inp` | 带谐振 `RESTRAINT` 的 NVT 平衡；无偏 MD（`DO_HILLS F`）测量 CV 涨落，用来选高斯宽度 |
| Well-tempered 元动力学（CP2K 原生） | `nacl_dissociation_wtmetad.inp` | 1D：偏置 $d$(Na–Cl)，`DELTA_T 1500`（γ = 6），上限墙 5.61 Å |
| 2D WT-MetaD | `nacl_wtmetad_2d.inp` | $d$(Na–Cl) + Na–O 配位数，`WTGAMMA 8` |
| PLUMED 驱动的 WT-MetaD | `plumed_metad.inp` + `plumed.dat` | 与 1D 原生输入物理上等价，需要带 PLUMED 编译的 CP2K |
| Blue moon 约束 MD | `constrained_md_distance.inp` + `run_bluemoon.sh` | 每个窗口固定 $d$，对 SHAKE 拉格朗日乘子求平均后积分 |

所有输入都通过 `@INCLUDE 'nacl_${LEVEL}_force_eval.inc'` 共用同一个 `&FORCE_EVAL`：`-E LEVEL=xtb` 选 GFN1-xTB，`-E LEVEL=dft` 选 PBE-D3(BJ)。两个 include 文件里定义的 `&COLVAR` 完全相同，所以同一套 CV 编号在两个层次上都适用。

## 文件说明

| 文件 | 说明 |
|---|---|
| `nacl_water.xyz` | 起始结构：1 = Na，2 = Cl，之后按 O H H 排列（O 为 3, 6, …, 186），共 188 个原子；由格点生成，**未经平衡** |
| `nacl_xtb_force_eval.inc` | GFN1-xTB（CP2K 原生实现），OT + ASPC；COLVAR 1 = $d$(Na–Cl)，COLVAR 2 = Na 被 O 配位的配位数（R0 = 3.0 Å，NN 12，ND 24） |
| `nacl_dft_force_eval.inc` | PBE-D3(BJ)，DZVP-MOLOPT-PBE-GTH（UZH 文件），400/50 Ry，OT + ASPC，共 512 个价电子；COLVAR 定义同上 |
| `nacl_equil.inp` | 步骤 1：5 ps NVT 平衡（CSVR，`TIMECON 50`），用 `&COLLECTIVE/&RESTRAINT` 把 $d$ 拉在 2.80 Å 附近 |
| `cv_monitor.inp` | 步骤 1b：无偏 NVT，`DO_HILLS F`，只输出 COLVAR，用于测量 CV 涨落 |
| `nacl_dissociation_wtmetad.inp` | 步骤 2：1D WT-MetaD（`NT_HILLS 40`，`WW 4.0E-4`，`SCALE 0.35` bohr，`DELTA_T 1500`，二次墙 10.6 bohr） |
| `nacl_wtmetad_2d.inp` | 步骤 2c：2D WT-MetaD（`WTGAMMA 8`，CN 方向带下限墙） |
| `plumed_metad.inp` | 步骤 2b：`USE_PLUMED T` + `PLUMED_INPUT_FILE plumed.dat` |
| `plumed.dat` | PLUMED 输入：`DISTANCE`、`COORDINATION`、`METAD`（BIASFACTOR 6）、`UPPER_WALLS`，单位 Å / kJ/mol / fs |
| `constrained_md_distance.inp` | 步骤 3：一个 blue-moon 窗口，`CONSTRAINT_INIT T`，`&LAGRANGE_MULTIPLIERS` 写到 `bluemoon.LagrangeMultLog` |
| `run_bluemoon.sh` | 依次（`CHAIN=1`）或独立（`CHAIN=0`）运行各窗口 |
| `pick_windows.py` | 从 xTB 元动力学轨迹中为每个目标距离挑一帧，写成 `window_<t>/start.xyz`（只用标准库） |
| `fes_from_hills.py` | 由 `HILLS.metadynLog`（或全精度的 `-1.restart`）重构 1D 自由能曲线或 2D 自由能面并投影到一维；`--blocks` 做收敛检查（需要 numpy） |
| `bluemoon_integrate.py` | 读取各窗口的乘子，做块平均误差估计和梯形积分，`--pmf` 再加 $2k_BT\ln r$ 项（需要 numpy） |
| `last_frame.sh` | 取 xyz 轨迹的最后一帧 |
| `job_slurm.sh` | Slurm 脚本，用 `STEP=` 选择步骤：`equil`、`cvmon`、`wtmetad`、`wt2d`、`plumed`、`bluemoon`；可选 `STEPS=` 覆盖步数 |

## 运行步骤

```bash
cd CP2K/examples/15-metadynamics-nacl
export CP2K_DATA_DIR=/path/to/cp2k/data          # 若已内置可省略

# 1. 平衡（xTB，5 ps），取最后一帧作为后续计算的起点
mpirun -np 16 cp2k.psmp -i nacl_equil.inp -o nacl_equil_xtb.out
bash last_frame.sh nacl_equil_xtb-pos-1.xyz nacl_eq.xyz

# 1b. 无偏 MD 看 CV 涨落（xTB，10 ps）
mpirun -np 16 cp2k.psmp -i cv_monitor.inp -o cv_monitor_xtb.out

# 2. xTB 预采样：1D WT-MetaD，100 ps
mpirun -np 16 cp2k.psmp -E LEVEL=xtb -E STEPS=200000 -i nacl_dissociation_wtmetad.inp -o wtmetad_xtb.out
python3 fes_from_hills.py nacl_wtmetad_xtb-HILLS.metadynLog --bohr2ang --blocks 4 -o fes_xtb.dat

# 2'. DFT 层次的 WT-MetaD（默认 20 ps，只够看出定性形状）
mpirun -np 64 cp2k.psmp -i nacl_dissociation_wtmetad.inp -o wtmetad_dft.out
python3 fes_from_hills.py nacl_wtmetad_dft-HILLS.metadynLog --bohr2ang --blocks 4 -o fes_dft.dat
# HILLS 中的高度只有 5 位小数（hartree），定量结果请用重启文件里的全精度高斯
python3 fes_from_hills.py nacl_wtmetad_dft-1.restart --bohr2ang --blocks 4 -o fes_dft_restart.dat

# 2b. PLUMED 版本（先确认 cp2kflags 里有 plumed2）
cp2k.psmp -v | grep plumed2
mpirun -np 64 cp2k.psmp -i plumed_metad.inp -o plumed_metad.out
plumed sum_hills --hills HILLS --mintozero --outfile fes_plumed.dat

# 2c. 2D WT-MetaD（xTB，200 ps）
mpirun -np 16 cp2k.psmp -i nacl_wtmetad_2d.inp -o wt2d_xtb.out
python3 fes_from_hills.py nacl_wt2d_xtb-HILLS.metadynLog --bohr2ang 1 -T 300 -o fes_d_2d.dat

# 3. DFT blue moon：从 xTB 轨迹挑起点，各窗口相互独立
python3 pick_windows.py nacl_wtmetad_xtb-pos-1.xyz --box 12.42 \
    --targets 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20
CHAIN=0 bash run_bluemoon.sh                      # 或每个窗口单独提交一个作业
python3 bluemoon_integrate.py window_*/ --skip 2000 --pmf -T 300 -o bluemoon_fes.dat

# 也可以用 Slurm（STEPS 可选，覆盖输入中的默认步数）
sbatch --export=ALL,STEP=equil,LEVEL=xtb job_slurm.sh
sbatch --export=ALL,STEP=wtmetad,LEVEL=xtb,STEPS=200000 job_slurm.sh
```

各脚本的全部选项请用 `python3 <脚本> -h` 查看。

## 预期结果 / 如何检查

```bash
grep "PROGRAM ENDED AT" *.out                                   # 每个计算都应正常结束
tail -3 nacl_wtmetad_xtb-COLVAR.metadynLog                       # t、d(bohr)、高斯力、墙力、V_bias、V_wall
wc -l nacl_wtmetad_xtb-HILLS.metadynLog                          # 高斯数 = STEPS / NT_HILLS
awk '{print $1, $NF}' nacl_wtmetad_xtb-HILLS.metadynLog | tail   # WT：高度应随时间明显下降
grep -c "Shake  Lagrangian Multipliers" window_3.20/bluemoon.LagrangeMultLog
grep -i "not converged" window_*/bm.out                          # SHAKE/SCF 不应出现未收敛
```

- **COLVAR 曲线**：WT-MetaD 中 $d$ 应在接触离子对（约 2.7–2.9 Å）和溶剂分隔离子对（约 4.5–5.2 Å）之间来回多次；只跑过去一次说明远没有收敛。
- **自由能曲线**：在定性上应能看到接触离子对极小、3.5–4 Å 附近的能垒和溶剂分隔离子对极小（文献参考见教程“结果解读”一节）；`--blocks 4` 输出的后两条曲线在相关区间内应基本重合。具体数值依赖泛函、水模型质量和采样长度，**本仓库没有提供计算结果**。
- **xTB 与 DFT**：GFN1-xTB 对液态水只是定性描述，自由能曲线的形状和能垒与 PBE-D3 会有明显差别。xTB 只用来摸清 CV 范围、参数和起点，定量结论请用 DFT 结果。
- **Blue moon**：`bluemoon_fes.dat` 中每个窗口的误差（块平均）应远小于相邻窗口平均力之差；误差大的窗口要延长。

## 注意事项

- `nacl_water.xyz` 是格点起始结构，一定要先跑 `nacl_equil.inp`（≥ 5 ps）；
- 所有 CP2K 内部的 CV、`SCALE`、墙位置和 `RESTRAINT K` 都使用**内部单位**（距离 bohr，能量 hartree），只有写了 `[angstrom]` 的 `TARGET` 例外；`plumed.dat` 中按 `UNITS` 行使用 Å 和 kJ/mol；
- 墙必须放在盒长的一半以内（这里 5.61 Å < 6.21 Å），否则最小镜像距离会“折返”；
- 续算元动力学：取消输入里 `&EXT_RESTART` 的注释，`STEPS` 写**本次还要走的**步数（步号从重启文件恢复，需要总步数上限时在 `&MD` 加 `MAX_STEPS`），已放置的高斯保存在 `.restart` 文件里；
- CP2K 以追加方式写 `*.metadynLog` 和 `*.LagrangeMultLog`：同一目录从头重跑前请先删除或改名旧文件，否则新旧高斯会叠加（`run_bluemoon.sh` 会自动删除各窗口的旧乘子文件）；
- PLUMED 版本会把 `HILLS`、`COLVAR` 和日志 `PLUMED.OUT` 写到当前目录；CP2K 不把温度传给 PLUMED，所以 `plumed.dat` 里的 `TEMP=300` 不能删；续算时在 `plumed.dat` 的 `UNITS` 行之后加一行 `RESTART`，否则 PLUMED 会备份旧文件并从零开始；
- 示例的步数和核数只是起点，请按自己的机器调整。没有运行结果的数值都不要当作参考值。

返回：[第15章教程](../../tutorials/15-metadynamics-plumed.md) | [CP2K 目录](../../README.md)
