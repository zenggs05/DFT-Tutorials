# 示例 02a：第一个 CP2K 计算——H2O 分子与 Si 晶体（First Calculations: H2O & Bulk Si）

## 目的

熟悉 CP2K 输入文件的结构和输出文件的阅读方法：

- 孤立 H2O 分子的单点能量（`RUN_TYPE ENERGY`）与能量 + 原子力（`RUN_TYPE ENERGY_FORCE`）；
- 周期性 Si 8 原子立方晶胞的单点能量，用预处理器从命令行切换 Gamma 点 / Monkhorst–Pack k 点；
- 预处理器 `@SET`、`${VAR}`、`${VAR-默认值}`、`@INCLUDE`、`@IF/@ENDIF`、`@PRINT` 与命令行 `-E` 的用法。

对应教程 [第2章 输入文件结构与第一个计算](../../tutorials/02-input-structure.md) 第 6、8、9、10 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o.xyz` | 水分子坐标（标准 xyz：原子数、注释行、3 行坐标，Å），实验平衡构型 r(OH) = 0.9572 Å，∠HOH = 104.52° |
| `h2o_energy.inp` | H2O 单点能量：PBE，`DZVP-MOLOPT-GTH` + `GTH-PBE`，`CUTOFF 400`/`REL_CUTOFF 50`，非周期 MT 泊松求解器，OT |
| `h2o_energy_force.inp` | 同上，`RUN_TYPE ENERGY_FORCE` 并打印原子力，SCF 更严格 |
| `si_bulk8_energy.inp` | Si 8 原子立方晶胞（a = 5.431 Å），对角化 + Broyden + Fermi–Dirac 展宽；`-E NK=n` 启用 n×n×n k 点 |
| `preprocessor_demo.inp` | 预处理器演示；变量 `BASIS`、`CUTOFF`、`USE_OT`、`PL` 可用 `-E` 修改 |
| `h2o_kinds.inc` | 被 `preprocessor_demo.inp` 用 `@INCLUDE` 引入的 `&KIND` 段 |

所有输入适用于 CP2K 2023.x–2026.2。

## 运行步骤

```bash
cd CP2K/examples/02-first-calc-h2o-si

# 0) 语法检查
for f in h2o_energy h2o_energy_force si_bulk8_energy preprocessor_demo; do cp2k.psmp -c $f.inp; done

# 1) H2O 单点能量与力（-o 会追加写入，先删除旧输出）
rm -f h2o_energy.out h2o_ef.out
OMP_NUM_THREADS=2 cp2k.ssmp -i h2o_energy.inp -o h2o_energy.out
OMP_NUM_THREADS=2 cp2k.ssmp -i h2o_energy_force.inp -o h2o_ef.out

# 2) Si：Gamma 点与 4x4x4 k 点
rm -f si8_k0.out si8_k4.out
cp2k.ssmp -i si_bulk8_energy.inp -o si8_k0.out
cp2k.ssmp -E NK=4 -i si_bulk8_energy.inp -o si8_k4.out

# 3) 预处理器
cp2k.ssmp -e preprocessor_demo.inp | less          # 展开后的完整输入（含全部默认值）
rm -f demo*.out
cp2k.ssmp -i preprocessor_demo.inp -o demo.out
cp2k.ssmp -E BASIS=TZV2P -E CUTOFF=500 -i preprocessor_demo.inp -o demo_tzv2p.out
cp2k.ssmp -E USE_OT=0 -i preprocessor_demo.inp -o demo_diag.out
```

MPI 版本：把 `cp2k.ssmp` 换成 `mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp`（OpenMPI）或 `mpiexec -n 4 -genv OMP_NUM_THREADS=1 cp2k.psmp`（MPICH）。

## 预期结果 / 如何检查

```bash
grep "SCF run converged" *.out
grep "ENERGY| Total FORCE_EVAL" *.out
grep -A6 "FORCES| Atomic forces" h2o_ef.out        # CP2K >= 2025.1
grep -A6 "ATOMIC FORCES" h2o_ef.out                # CP2K <= 2024.x
grep -A6 "Mulliken Population Analysis" h2o_energy.out
grep "count for grid" h2o_energy.out
grep "BRILLOUIN|" si8_k4.out
grep "GLOBAL| Project name" demo*.out
grep "PROGRAM ENDED AT" *.out
```

- H2O：SCF 应在十几步内收敛，总能量约 −17.2 Ha（示意；具体数值取决于版本与网格）。`h2o_ef.out` 中原子力之和（`Sum` / `SUM OF ATOMIC FORCES`）应接近 0；各原子受力不为零，因为实验构型不是 PBE/DZVP 的能量极小点。
- Si：Gamma 点与 4×4×4 的每原子能量差别明显（8 原子晶胞太小，Gamma 点不够），说明需要 k 点收敛。
- 预处理器：`demo_tzv2p.out` 的项目名为 `h2o_TZV2P_c500_ot1`；`demo_diag.out` 中 SCF 使用对角化（`Update method` 列不是 `OT`）。`-E XC=BLYP` 无效，因为文件中的 `@SET XC PBE` 会覆盖命令行初值。

## 注意事项

- 本目录的输入需要在**本目录**中运行：`h2o.xyz` 与 `h2o_kinds.inc` 按当前工作目录查找。
- 这里的 `CUTOFF`/`REL_CUTOFF`、k 点数只是教学用的起点，正式计算前按 [第4章](../../tutorials/04-cutoff-convergence.md)、[第6章](../../tutorials/06-periodic-surface-kpoints.md) 收敛。
- Si 输入使用对角化而不是 OT，因为 CP2K ≤ 2026.2 的 OT 不支持 k 点。
- 预处理器会处理非注释行中出现的任何 `$` 与 `@`，不要在行尾注释中写这两个字符；`@SET` 行末不要加注释。

返回教程：[第2章 输入文件结构与第一个计算](../../tutorials/02-input-structure.md)
