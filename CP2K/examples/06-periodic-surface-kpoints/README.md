# 示例 06：孤立分子、表面 slab 偶极校正与 k 点收敛

配套教程：[第6章 周期性体系、表面模型与 k 点](../../tutorials/06-periodic-surface-kpoints.md)

## 目的

| 输入 | 体系 | 演示内容 |
|---|---|---|
| `h2o_isolated_mt.inp` | 孤立 H2O | `PERIODIC NONE` 在 `&CELL` 和 `&POISSON` 两处一致，`POISSON_SOLVER MT`，`&CENTER_COORDINATES`；`@SET BOX` 用于盒子尺寸测试 |
| `pt111_slab_dipole.inp` | Pt(111) p(2×2) 4 层 + O（fcc 位） | 3D 周期 + `SURFACE_DIPOLE_CORRECTION`/`SURF_DIP_DIR Z`，k 点 4×4×1，金属展宽；`@SET RUNTYPE GEO_OPT` 时固定底部两层 |
| `si_primitive_kpoints.inp` | 体相 Si 原胞（2 原子） | `ABC` + `ALPHA_BETA_GAMMA 60 60 60` + 分数坐标；`&KPOINTS MONKHORST-PACK 8 8 8`；k 点下用对角化 |
| `kpoint_convergence_template.inp` + `run_kpoint_scan.sh` | 同上 | `@SET NK` 模板 + 扫描脚本，汇总 k 点数、能量和压力 |

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_isolated_mt.inp` | 孤立分子示例（RUN_TYPE ENERGY，OT） |
| `h2o.xyz` | H2O 坐标（实验气相几何） |
| `pt111_slab_dipole.inp` | slab 示例（默认 RUN_TYPE ENERGY_FORCE；把顶部的 `@SET RUNTYPE` 改为 `GEO_OPT` 即做结构优化） |
| `pt111_o_p2x2.xyz` | 用 ASE `fcc111('Pt', size=(2,2,4), a=3.97)` 生成并加 O 的初始结构（未优化），晶胞：ABC 5.614428 5.614428 26.0，角度 90 90 60，原子 1–8 为底部两层 |
| `si_primitive_kpoints.inp` | Si 原胞 + 8×8×8 MP 网格（RUN_TYPE ENERGY_FORCE） |
| `kpoint_convergence_template.inp` | k 点收敛模板，`@SET NK` 由脚本改写；开启 `STRESS_TENSOR ANALYTICAL` 以同时检查压力 |
| `run_kpoint_scan.sh` | 生成 `kp_<NK>/` 目录、运行、汇总到 `kpoint_scan.csv`；环境变量 `NKS`、`CP2K`、`NATOMS`、`DRYRUN`、`COLLECT_ONLY` |
| `job_slurm.sh` | Slurm 作业脚本（默认 32 核），可传入要运行的输入；`KSCAN=1` 时附带运行 k 点扫描 |

## 运行步骤

```bash
cd CP2K/examples/06-periodic-surface-kpoints
export CP2K_DATA_DIR=/path/to/cp2k/data          # 若已内置可省略
for f in h2o_isolated_mt.inp pt111_slab_dipole.inp si_primitive_kpoints.inp kpoint_convergence_template.inp; do
  cp2k.psmp -c "$f" > /dev/null && echo "$f OK"
done

# 1) isolated molecule + box-size test
mpirun -np 4 cp2k.psmp -i h2o_isolated_mt.inp -o h2o_isolated_mt.out
for L in 10 12 14 16; do
  sed -E "s/^@SET[[:space:]]+BOX[[:space:]].*/@SET BOX ${L}.0/" h2o_isolated_mt.inp > box_${L}.inp
  mpirun -np 4 cp2k.psmp -i box_${L}.inp -o box_${L}.out
done
grep "ENERGY| Total FORCE_EVAL" box_*.out

# 2) Si primitive cell with k-points, then the k-point convergence scan
mpirun -np 8 cp2k.psmp -i si_primitive_kpoints.inp -o si_primitive_kpoints.out
CP2K="mpirun -np 8 cp2k.psmp" bash run_kpoint_scan.sh

# 3) Pt(111)+O slab (17 atoms, 4x4x1 k-points; use a cluster node)
sbatch job_slurm.sh pt111_slab_dipole.inp
```

## 预期结果 / 如何检查

```bash
grep "SCF run converged" *.out
grep "BRILLOUIN| List of Kpoints" si_primitive_kpoints.out pt111_slab_dipole.out
grep -A1 "Total dipole moment perpendicular" pt111_slab_dipole.out | tail -2
grep -i "smearing includes the last MO" pt111_slab_dipole.out    # should be empty
column -s, -t kpoint_scan.csv
```

- `h2o_isolated_mt`：盒长从 10 Å 增加到 16 Å 时，总能量应趋于平台；若 10 Å 与 12 Å 差别明显，说明 MT 盒子不够大；
- `si_primitive_kpoints`：`BRILLOUIN| List of Kpoints` 为 256（8³/2，默认只做时间反演配对）；
- `kpoint_scan.csv`：`dE_meV_per_atom` 随 NK 增大趋于 0，压力收敛得比能量慢；
- `pt111_slab_dipole`：每个 SCF 步都打印 `Total dipole moment perpendicular to the slab [electrons-Angstroem]`。若出现 `Dipole correction needs more vacuum space above the surface`，说明真空不够或 slab 跨过了晶胞边界。

本仓库没有附带参考输出；数值取决于 CP2K 版本、编译方式和并行设置。

## 注意事项

- 截断（350/500 Ry）、电子温度（500 K）和 k 网格都是演示值，正式计算请按[第4章](../../tutorials/04-cutoff-convergence.md)、[第6章](../../tutorials/06-periodic-surface-kpoints.md)的方法收敛；Pt 的晶格常数 3.97 Å 是近似的 PBE 值，正式计算应先用同一套参数优化体相晶格常数；
- 2026.2 及以前，k 点计算（包括 `SCHEME GAMMA`）不能用 OT；`h2o_isolated_mt.inp` 没有 `&KPOINTS`，因此可以用 OT；
- `SURFACE_DIPOLE_CORRECTION` 要求表面法向严格沿 `SURF_DIP_DIR` 指定的笛卡尔轴（本例 C 沿 z，A、B 在 xy 平面）；
- slab 示例的 `&MOTION` 段只在 `RUN_TYPE GEO_OPT`/`MD` 时生效；固定原子的编号与 `pt111_o_p2x2.xyz` 中的顺序对应，改结构后要同步修改；
- 2026.2 新增的 `GAMMA_CENTERED`、`INVERSION_SYMMETRY_ONLY`、`SYMMETRY_BACKEND`、`EPS_SYMMETRY` 等关键字在旧版本中不存在，本目录的输入没有使用它们，因此可在 2024.1–2026.2 上运行。
