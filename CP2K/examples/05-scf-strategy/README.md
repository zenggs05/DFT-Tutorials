# 示例 05：SCF 策略——OT、对角化 + 混合、金属展宽与波函数重启

配套教程：[第5章 SCF 策略：OT、对角化、展宽与重启](../../tutorials/05-scf-strategy.md)

## 目的

用四个小体系演示 CP2K Quickstep 的几种 SCF 求解方式：

| 输入 | 体系 | 演示内容 |
|---|---|---|
| `h2o_ot.inp` | 孤立 H2O（12 Å 盒子，MT 求解器） | OT（`DIIS` + `FULL_ALL`）+ `&OUTER_SCF`，写出 `h2o_ot-RESTART.wfn` |
| `restart_from_wfn.inp` | 同上，TZV2P 基组，600 Ry | `SCF_GUESS RESTART` + `&DFT/WFN_RESTART_FILE_NAME`，从 DZVP 波函数起步 |
| `si_diag_broyden.inp` | 体相 Si（8 原子立方胞），MP 4×4×4 | k 点下必须用对角化；`BROYDEN_MIXING`；写出 `.kp` 重启文件 |
| `al_metal_smear.inp` | fcc Al（4 原子惯用胞），MP 8×8×8 | 金属：对角化 + Broyden（ALPHA 0.2）+ Fermi–Dirac 500 K + `ADDED_MOS 10` |

所有输入都用 PBE、`BASIS_MOLOPT_UZH` 加 `POTENTIAL_UZH`（DZVP/TZV2P-MOLOPT-PBE-GTH-qN 配 GTH-PBE-qN），只使用 CP2K 2024.1–2026.2 都支持的关键字。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_ot.inp` | OT 示例（RUN_TYPE ENERGY_FORCE） |
| `restart_from_wfn.inp` | 波函数重启示例，**必须在 `h2o_ot.inp` 之后运行** |
| `si_diag_broyden.inp` | 对角化 + Broyden + k 点示例（RUN_TYPE ENERGY） |
| `al_metal_smear.inp` | 金属展宽示例（RUN_TYPE ENERGY_FORCE），坐标直接写在 `&COORD SCALED` 中 |
| `h2o.xyz` | H2O 坐标（实验气相几何），用 `&CENTER_COORDINATES` 放到盒子中心 |
| `si_bulk8.xyz` | Si8 坐标（a = 5.431 Å） |
| `job_slurm.sh` | Slurm 作业脚本，默认按正确顺序运行全部四个输入，也可只运行指定的输入 |

## 运行步骤

```bash
cd CP2K/examples/05-scf-strategy
export CP2K_DATA_DIR=/path/to/cp2k/data      # 若已内置可省略

# 语法检查（逐个）
for f in *.inp; do cp2k.psmp -c "$f" > /dev/null && echo "$f OK"; done

# 本地运行
mpirun -np 4 cp2k.psmp -i h2o_ot.inp           -o h2o_ot.out
mpirun -np 4 cp2k.psmp -i restart_from_wfn.inp -o restart_from_wfn.out
mpirun -np 8 cp2k.psmp -i si_diag_broyden.inp  -o si_diag_broyden.out
mpirun -np 8 cp2k.psmp -i al_metal_smear.inp   -o al_metal_smear.out

# 或提交到 Slurm
sbatch job_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep "SCF run converged"   *.out            # 每个输出都应出现
grep "outer SCF loop"      h2o_ot.out       # OT 外层循环收敛信息
grep -E "OT (DIIS|CG|LS)"  h2o_ot.out | head
grep -i "does not exist"   restart_from_wfn.out   # 应为空：说明重启文件被找到
grep "BRILLOUIN| List of Kpoints" si_diag_broyden.out al_metal_smear.out
grep -E "Broy\./Diag\."    si_diag_broyden.out | head
grep -E "Fermi energy|Electronic entropic energy" al_metal_smear.out | tail -2
grep -i "smearing includes the last MO" al_metal_smear.out   # 应为空
ls *-RESTART.wfn *-RESTART.kp
```

- `h2o_ot.out`：SCF 迭代行的方法列为 `OT DIIS`（或 `OT LS`），最后出现 `outer SCF loop converged in ... iterations`；
- `restart_from_wfn.out`：SCF 步数应少于从 `ATOMIC` 起步的同一计算（可以把 `SCF_GUESS` 改成 `ATOMIC` 对比）；
- `si_diag_broyden.out`：第一步是 `NoMix/Diag.`，之后是 `Broy./Diag.`，并生成 `si8_diag-RESTART.kp`（k 点计算写 `.kp` 而不是 `.wfn`）；
- `al_metal_smear.out`：打印 `Electronic entropic energy`（即 −TS，为负的小量）和 `Fermi energy`；用 CP2K ≥ 2026.2 时还会打印 `Total energy (extrapolated to T->0)`。

本仓库没有附带参考输出；能量、步数和耗时取决于 CP2K 版本、编译方式和并行设置。

## 注意事项

- 截断和 k 点网格（H2O 500/600 Ry，Si 350 Ry + 4×4×4，Al 300 Ry + 8×8×8）只是演示值，正式计算请按[第4章](../../tutorials/04-cutoff-convergence.md)、[第6章](../../tutorials/06-periodic-surface-kpoints.md)做收敛测试；
- 2026.2 及以前，`&OT` 不能与 `&KPOINTS` 或 `ADDED_MOS > 0` 同时使用，分别报 `OT not possible with kpoint calculations` 和 `OT with ADDED_MOS/=0 not implemented`；
- `&SMEAR` 中务必显式写 `METHOD`：≤ 2026.1 的输入参考把默认值标为 `ENERGY_WINDOW`，2026.2 改为 `GAUSSIAN`；
- `WFN_RESTART_FILE_NAME` 属于 `&DFT`，不属于 `&SCF`；从 DZVP 重启到 TZV2P 时赝势必须相同；
- 不要在这些输入中加 `IGNORE_CONVERGENCE_FAILURE`；若某个计算不收敛，按教程第 10 节逐项排查。
