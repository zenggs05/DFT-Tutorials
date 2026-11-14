# 示例 06-1：分子结构优化——H2O（RF）与 CH4（最速下降）

## 目的

- 用 `MD.Type RF`（有理函数法）优化一个故意扭曲的水分子，学会读 `.out` 中的优化历史和 `.md` 轨迹；
- 用 `MD.Type Opt`（最速下降）优化 C 原子偏离中心 0.3 Å 的甲烷，体会最速下降"稳但慢、会越过驻点"的特点；
- 练习用 `输入文件名#` 续算。

对应教程：[第6章 结构优化、晶胞优化、状态方程与 NEB](../../tutorials/06-optimization-eos-neb.md)（§3）。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `H2O_opt_RF.dat` | `h2o_rf` | 初始 r(OH) = 1.012 Å、∠HOH = 114.2°；RF，`MD.Opt.criterion 1.0e-4`，`MD.maxIter 100` |
| `Methane2_SD.dat` | `ch4_sd` | C 沿 x 偏离 0.3 Å（正四面体中 C–H = 1.09 Å）；`MD.Type Opt`，默认判据 3.0e-4 |
| `job.slurm` | — | Slurm 作业脚本 |

两个输入都是团簇计算（`scf.EigenvalueSolver Cluster`），10 Å 立方盒子，250 Ry，Standard 基组。

## 运行步骤

```bash
cd OpenMX/examples/06-geomopt-h2o-ch4
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

mpirun -np 4 openmx H2O_opt_RF.dat -nt 1 > h2o_rf.std
mpirun -np 4 openmx Methane2_SD.dat -nt 1 > ch4_sd.std
# or: sbatch job.slurm H2O_opt_RF.dat ; sbatch job.slurm Methane2_SD.dat

# compare with RF for methane
sed -e 's/^MD.Type .*/MD.Type                     RF/' \
    -e 's/^System.Name .*/System.Name                 ch4_rf/' Methane2_SD.dat > Methane2_RF.dat
mpirun -np 4 openmx Methane2_RF.dat -nt 1 > ch4_rf.std
```

## 预期结果 / 如何检查

```bash
grep -A30 "History of" h2o_rf.out              # MD_iter, SD_scaling, |Maximum force|, Maximum step, Utot
grep -o "scf_conv=[01]" h2o_rf.md | sort | uniq -c
cat h2o_rf.xyz                                 # final geometry
ls *.dat#                                      # restart inputs, updated every step
```

- 收敛时最后一步的最大力 < `MD.Opt.criterion`，且总步数小于 `MD.maxIter`；
- 水分子应收敛到 r(OH) ≈ 0.97 Å、∠HOH ≈ 104° 量级（官方手册表 3 的 PBE 结果为 0.971 Å、103.6–103.7°，所用基组与本例不同；实验 0.957 Å、104.5°）；
- 甲烷应恢复为正四面体，C–H 键长约 1.09–1.10 Å；
- 最速下降的步数明显多于 RF，最大力曲线可能先越过驻点再回落；
- 用 OpenMX Viewer 打开 `h2o_rf.md` 可以播放优化过程；本仓库脚本 [`omx_traj.py`](../../scripts/omx_traj.py)（`python ../../scripts/omx_traj.py -h` 查看用法）可画能量和最大力随步数的变化。

## 注意事项

- `MD.maxIter` 默认是 1，一定要显式设置。
- `MD.Opt.Init.Hessian Schlegel` 只能用于 RF 和 EF；改用 `BFGS` 或 `DIIS` 时请改成 `iden` 或删除该行。
- 没收敛时：`cp H2O_opt_RF.dat# H2O_opt_RF_r2.dat`，然后运行新文件即可续算，近似 Hessian 会被自动重用。
- 输入文件只含 ASCII 字符和 LF 换行。
