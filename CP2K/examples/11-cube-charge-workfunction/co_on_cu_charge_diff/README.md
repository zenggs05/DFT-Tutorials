# 子示例：CO/Cu(111) 差分电荷密度

配套教程：[第11章 第7节 差分电荷密度](../../../tutorials/11-cube-charge-workfunction.md)；上级示例说明见 [../README.md](../README.md)。

## 目的

计算 CO 顶位吸附在 Cu(111) p(2×2) 三层 slab 上的差分电荷密度

$$\Delta\rho(\mathbf r)=\rho_{AB}(\mathbf r)-\rho_{A}(\mathbf r)-\rho_{B}(\mathbf r)$$

其中 AB 为吸附体系，A 为干净 slab，B 为单独的 CO。三次计算的原子位置完全相同（从同一套坐标中删去另一部分），晶胞、`CUTOFF`、`REL_CUTOFF`、k 点和 cube 的 `STRIDE` 也完全相同，这样三个 cube 文件位于同一个实空间网格上，可以逐点相减。

## 文件说明

| 文件 | 说明 |
|---|---|
| `dft_common.inc` | 三个输入共用的 `&DFT` 设置（`@INCLUDE` 引入）；`ADDED_MOS ${NADD}` 由各输入用 `@SET` 指定 |
| `full.inp` / `cu111_co_full.xyz` | 体系 AB：12 Cu + C + O |
| `slab.inp` / `cu111_slab.xyz` | 体系 A：12 Cu（坐标与 AB 中的 Cu 完全相同） |
| `co.inp` / `co_only.xyz` | 体系 B：CO（坐标与 AB 中的 C、O 完全相同） |
| `chgdiff.py` | 读三个 cube、检查网格一致、写 `diff.cube` 与平面平均 `diff_planar.dat`（只需 numpy） |
| `run_all.sh` | 依次运行三个计算，再调用 `chgdiff.py` |

结构说明：Cu 晶格常数取实验值 3.615 Å，层间距 $a/\sqrt3$；C 在顶层 Cu 正上方 1.85 Å，C–O 1.15 Å（近似值，未优化）；c = 28 Å，真空约 21 Å。正式研究应先用同一套参数优化 AB 结构，再从优化后的结构中拆出 A 和 B。

## 运行步骤

```bash
cd CP2K/examples/11-cube-charge-workfunction/co_on_cu_charge_diff
export CP2K_DATA_DIR=/path/to/cp2k/data
cp2k.psmp -c full.inp            # 语法检查（会展开 @INCLUDE）
bash run_all.sh                  # 本地：默认 mpirun -np 4 cp2k.psmp
# Slurm 分配内：CP2K="srun cp2k.psmp" bash run_all.sh
```

## 预期结果 / 如何检查

- 三个输出都有 `SCF run converged` 和 `PROGRAM ENDED AT`；
- 生成 `full-ELECTRON_DENSITY-1_0.cube`、`slab-ELECTRON_DENSITY-1_0.cube`、`co-ELECTRON_DENSITY-1_0.cube`（文件名以 `ls` 为准）；
- `chgdiff.py` 打印的 `integral of drho` 应接近 0（三个体系的价电子数满足 $N_{AB}=N_A+N_B$，残差来自数值积分误差）；
- 在 VESTA 中打开 `diff.cube`，取 ±0.002–0.005 e/bohr³ 的等值面：CO 的 5σ 区域和 Cu–C 之间应能看到电荷重新分布（σ 给予、π 反馈的图像），具体形状以计算结果为准；
- `diff_planar.dat` 的第 3 列 $\Delta Q(z)$ 描述沿表面法向的电荷转移，在真空区应回到约 0。

## 注意事项

- **不要**在三个输入中使用 `&TOPOLOGY/&CENTER_COORDINATES`，它会把每个子体系分别居中，导致三套坐标不再重合；
- 三次计算的晶胞、`CUTOFF`/`REL_CUTOFF`、`STRIDE` 必须完全一致，`chgdiff.py` 会检查网格；
- CP2K 的 cube 只含价电子（赝势），差分密度在原子核附近的细节不具物理意义；
- k 点 4×4×1 只用于演示，金属表面的定量结果需要收敛 k 点；
- 体系 B 的基函数只有 26 个、占据态 5 个，所以 `ADDED_MOS` 取 10（`@SET NADD 10`）。
