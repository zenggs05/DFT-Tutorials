# 示例 07：CP2K 几何优化与晶胞优化（GEO_OPT / CELL_OPT）

配套教程：[第7章 几何优化与晶胞优化](../../tutorials/07-geo-cell-opt.md)

## 目的

用四个小体系演示 CP2K 中最常用的结构优化场景：

1. 孤立水分子的 GEO_OPT（BFGS，可切换 LBFGS/CG）；
2. MgO(001) 表面吸附水分子，固定底层原子的 slab GEO_OPT；
3. 体相 Si 的 CELL_OPT（解析应力、`CELL_REF`、`KEEP_SYMMETRY`、k 点 + 对角化）；
4. 层状 2H-MoS2 的 CELL_OPT（`KEEP_ANGLES`、DFT-D3(BJ) 色散校正）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_geo_opt_bfgs.inp` | 水分子几何优化，`@SET OPT` 可选 BFGS/LBFGS/CG |
| `h2o_distorted.xyz` | 故意扭曲的初始水分子（O–H 1.05/1.00 Å，键角 120°） |
| `slab_geo_opt_fixed.inp` | MgO(001) 2×2 三层 slab + H2O，固定第 1–16 号原子，含偶极校正和 D3(BJ) |
| `mgo001_h2o.xyz` | slab 坐标（51 个原子；晶胞 8.42 × 8.42 × 25.0 Å，坐标文件不含晶胞，晶胞写在输入中） |
| `si_cell_opt.inp` | 8 原子立方 Si 的 CELL_OPT，初始 a = 5.55 Å |
| `layered_cell_opt_keep_angles.inp` | 2H-MoS2 的 CELL_OPT，初始 c = 13.00 Å，`KEEP_ANGLES T` |
| `job_slurm.sh` | 通用 Slurm 作业脚本，用 `INP` 变量指定输入 |

所有输入使用 `BASIS_MOLOPT_UZH` + `POTENTIAL_UZH`，PBE 泛函。输入文件不含 `&CELL_OPT/TYPE`，因此在 CP2K 2023.x–2026.2 中都能运行（2026.2 删除了该关键字，旧版本默认值就是 `DIRECT_CELL_OPT`）。

## 运行步骤

```bash
# 0) environment (pick one) and data directory
source /path/to/cp2k/install/cp2k_env        # or: module load cp2k / conda activate cp2k_env
cp2k.psmp -v | head -20                      # version and cp2kflags
cp2k.psmp -c si_cell_opt.inp                 # syntax check before submitting

# 1) H2O, a few minutes on 4 cores
mpirun -np 4 cp2k.psmp -i h2o_geo_opt_bfgs.inp -o h2o_geo_opt_bfgs.out

#    same input with LBFGS / CG (edit the @SET lines; "-E" cannot override an @SET in the file)
sed -e 's/^@SET OPT .*/@SET OPT LBFGS/' -e 's/^@SET PROJ .*/@SET PROJ h2o_lbfgs/' \
    h2o_geo_opt_bfgs.inp > h2o_lbfgs.inp
sed -e 's/^@SET OPT .*/@SET OPT CG/' -e 's/^@SET PROJ .*/@SET PROJ h2o_cg/' \
    h2o_geo_opt_bfgs.inp > h2o_cg.inp

# 2)-4) larger jobs: submit to Slurm
sbatch --export=ALL,INP=slab_geo_opt_fixed.inp job_slurm.sh
sbatch --export=ALL,INP=si_cell_opt.inp job_slurm.sh
sbatch --export=ALL,INP=layered_cell_opt_keep_angles.inp job_slurm.sh
```

计算量参考：水分子几分钟；Si 和 MoS2 的 CELL_OPT 在 16–32 核上通常在一小时量级；MgO slab（51 原子、392 个价电子）的 GEO_OPT 需要更长时间（以上均为粗略估计，取决于机器和收敛步数）。

## 预期结果 / 如何检查

```bash
grep "GEOMETRY OPTIMIZATION COMPLETED" *.out            # every run should print this
grep -c "OPT| Step number" h2o_geo_opt_bfgs.out          # number of optimization steps
grep "OPT| Maximum gradient  " h2o_geo_opt_bfgs.out | tail -2
grep "CELL| Vector a" si_cell_opt.out | tail -1           # final lattice vector a
grep "OPT| Pressure deviation" si_cell_opt.out | tail -1
grep "CELL| Vector c" layered_cell_opt_keep_angles.out | tail -1
grep -i "SCF run NOT converged" *.out                    # should print nothing
```

- 水分子：PBE 下 O–H 约 0.97 Å、键角约 104°（近似值；实验 0.957 Å、104.5°）。
- Si：PBE 晶格常数一般在 5.46–5.47 Å 附近（近似范围，与基组、赝势、CUTOFF 有关；实验 5.431 Å）。初始的 `Internal pressure` 为负值，最后应在外压 ±100 bar 以内。
- MoS2：c 应从 13.00 Å 明显收缩，接近实验值 12.29 Å 附近（PBE-D3(BJ) 的具体数值请以自己的计算为准）；去掉 `&VDW_POTENTIAL` 再算一次，c 会显著偏大。
- slab：第 1–16 号原子的坐标在整个轨迹中保持不变；水分子 O 与表层 Mg 的距离在 2.1–2.3 Å 左右（示意范围）。检查方法见教程第 5.2 节。

可视化：`ase gui h2o_opt-pos-1.xyz`，或用 VESTA/Ovito 打开 `*-pos-1.xyz`。CELL_OPT 的晶胞随步数变化记录在 `*-1.cell` 中；CP2K 2026.2 还会写出 `*-FINAL-1_<步数>.cif` 和扩展 XYZ。作图可用 `python ../../scripts/plot_opt.py -h` 查看用法。

## 注意事项

- **CP2K 2026.2**：不要在 `&CELL_OPT` 中写 `TYPE`，否则直接中止；`KEEP_VOLUME`、`&FINAL_STRUCTURE`、`&CELL/CANONICALIZE`、`EXTRAPOLATION GEXT_PROJ` 只有 2026.2 才有。
- **k 点与 OT**：2026.2 及以前，k 点计算只能用对角化，因此 Si 和 MoS2 示例用 `&DIAGONALIZATION` + Broyden 混合。
- **CUTOFF**：晶胞优化的 600 Ry 只是示例值，正式计算请按[第4章](../../tutorials/04-cutoff-convergence.md)用应力做收敛测试。
- **MgO 晶格常数**：slab 示例直接用了实验值 4.21 Å，正式计算应改为同一设置下体相 CELL_OPT 的结果，并增加层数。
- **数据文件**：`BASIS_MOLOPT_UZH`、`POTENTIAL_UZH`、`dftd3.dat` 需在 `CP2K_DATA_DIR` 中可找到（toolchain/conda/容器版本一般已自动设置）。
- **续算**：中断后可直接 `cp <PROJECT>-1.restart new.inp` 再运行，或在原输入中加入 `&EXT_RESTART`，详见教程第 10.3 节。
