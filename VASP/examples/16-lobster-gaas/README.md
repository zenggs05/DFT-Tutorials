# 示例 16-5：GaAs 的 LOBSTER COHP / ICOHP 分析

## 目的

- 按 LOBSTER 的要求准备 VASP 单点（`ISYM = 0`、`NBANDS`、`LWAVE`、`NSW = 0`，不用 `vasp_gam`）；
- 写 `lobsterin`，计算 Ga–As 的 −COHP、−ICOHP、Mulliken/Löwdin 电荷；
- 检查 charge spilling，并用 LobsterPy 自动分析和作图。

对应教程：[第16章 电荷与成键分析](../../tutorials/16-charge-bonding.md)（§8）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 闪锌矿 GaAs 原胞，a = 5.653 Å（实验值），Ga–As 键长 ≈ 2.448 Å |
| `KPOINTS` | Γ 中心 8×8×8 |
| `INCAR` | `ISYM = 0`、`NBANDS = 16`（≥ 13 个基函数）、`LWAVE = .TRUE.`、`ISMEAR = 0`、`SIGMA = 0.05` |
| `lobsterin` | pbeVaspFit2015 基组，Ga 3d 4s 4p / As 4s 4p，`cohpGenerator from 2.0 to 2.6 type Ga type As orbitalwise` |
| `job.slurm` | 先运行 VASP，再运行 LOBSTER（OpenMP） |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/16-lobster-gaas
cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR
grep -E "TITEL|ENMAX|ZVAL" POTCAR

# LOBSTER 可执行文件名以你下载的版本为准
export LOBSTER_EXE=/path/to/lobster-5.1.1
sbatch job.slurm
```

不想手写 `lobsterin` 时，也可以在只有 POSCAR/POTCAR/INCAR 的目录中运行 `lobsterpy create-inputs`（目前针对 PBE_54 系列 POTCAR），比较它生成的文件与本例的差别。

## 预期结果 / 如何检查

```bash
grep -E "NBANDS" OUTCAR | head -2                       # 应为 16，且没有 "number of bands has been changed"
grep -i "spilling" lobsterout                            # charge spilling 越小越好，经验上约 < 1-2%
head -6 ICOHPLIST.lobster                                # Ga-As 原子对及其 ICOHP（负值 = 净成键）
cat CHARGE.lobster                                       # Mulliken 与 Loewdin 电荷

pip install lobsterpy
lobsterpy description                                    # 自动文字描述
lobsterpy plot-automatic                                 # 自动画最重要的 COHP
lobsterpy plot 1 2 --summed                              # 画第 1、2 个键的加和
```

- 闪锌矿中每个 Ga 有 4 个等价的 As 近邻，ICOHPLIST 中这些原子对的 ICOHP 应基本相同；
- −COHP 图中费米能级以下主要为正值（成键），费米能级以上为反键，符合半导体中价带成键、导带反键的图像；
- 具体数值请以自己的计算为准，本仓库不提供参考数值。

## 注意事项

- `vasp_gam` 产生的 WAVECAR 不能被 LOBSTER 读取；
- 若 VASP 自动改动了 NBANDS，调整 `NBANDS` 或并行参数（`KPAR`、`NCORE`）使其整除后重算；
- `lobsterin` 中的 `gaussianSmearingWidth` 应与 INCAR 的 `SIGMA` 一致（`ISMEAR = -5` 时删掉这一行）；
- LOBSTER 是学术免费软件，不允许用于商业研究；本例语法按 LOBSTER 5.x 编写，使用 6.0 时请对照新版手册；
- WAVECAR 体积较大，不要提交到仓库。
