# 示例 02b：金红石 TiO2 的能量与原子力（Rutile TiO2 Energy & Forces: Gamma vs k-points vs Supercell）

## 目的

用金红石（rutile）TiO2 演示一个周期性体系的 `RUN_TYPE ENERGY_FORCE` 计算，并说明一个常被忽视的问题：**6 原子原胞（4.593 × 4.593 × 2.959 Å）太小，只用 Gamma 点采样远未收敛**。本目录保留仓库原有的输入 `input_file.inp`（Gamma 点 + OT），并给出两种改进做法：

1. `rutile_kpoints.inp`：原胞 + Monkhorst–Pack k 点网格（对角化 + Broyden 混合）；
2. `rutile_supercell_gamma.inp`：用 `MULTIPLE_UNIT_CELL` 扩成超胞，仍用 Gamma 点 + OT。

对应教程 [第2章 输入文件结构与第一个计算](../../tutorials/02-input-structure.md) 第 11 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| `input_file.inp` | **原有输入（保持不变）**：PBE，`DZVP-MOLOPT-SR-GTH` + `GTH-PBE-q6/q12`，`CUTOFF 400`/`REL_CUTOFF 60`，OT + OUTER_SCF，**仅 Gamma 点**，`PROJECT Universality` |
| `rutile.xyz` | 原有坐标文件：6 行 `元素 x y z`（Å），**没有** xyz 文件通常的“原子数 + 注释”两行表头，专供 `&COORD` 中 `@INCLUDE` 使用 |
| `rutile_kpoints.inp` | 改进版 1：同样的结构、基组、赝势和泛函，加 `&KPOINTS SCHEME MONKHORST-PACK 4 4 6`；k 点数可用 `-E KX=… -E KY=… -E KZ=…` 修改 |
| `rutile_supercell_gamma.inp` | 改进版 2：`MULTIPLE_UNIT_CELL 2 2 3`（72 原子）的 Gamma 点超胞；超胞大小可用 `-E NX=… -E NY=… -E NZ=…` 修改 |
| `job_slurm.sh` | Slurm 作业脚本（默认运行 `input_file.inp`，可用 `--export=ALL,INP=...` 换输入） |
| `job_lsf.sh` | LSF 作业脚本（原 `cp2k.slurm` 的修正版） |
| `cp2k.slurm` | 已废弃：原文件名为 slurm 实为 LSF 脚本，且运行的是不存在的 `input.inp`；现在只打印一行提示，请改用上面两个脚本 |

## 运行步骤

```bash
cd CP2K/examples/02-energy-force-rutile

# 0) 语法检查（几秒钟）
cp2k.psmp -c input_file.inp
cp2k.psmp -c rutile_kpoints.inp
cp2k.psmp -c rutile_supercell_gamma.inp

# 1) 原始输入：Gamma 点（交互运行示例；集群上请用作业脚本）
rm -f input_file.out
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i input_file.inp -o input_file.out

# 2) k 点版本（默认 4x4x6）
rm -f rutile_k446.out
mpirun -np 4 -x OMP_NUM_THREADS=1 cp2k.psmp -i rutile_kpoints.inp -o rutile_k446.out

# 3) 超胞版本（默认 2x2x3，72 原子）
rm -f rutile_sc223.out
mpirun -np 8 -x OMP_NUM_THREADS=1 cp2k.psmp -i rutile_supercell_gamma.inp -o rutile_sc223.out

# 集群提交
sbatch job_slurm.sh                                        # 运行 input_file.inp
sbatch --export=ALL,INP=rutile_kpoints.inp job_slurm.sh
bsub < job_lsf.sh                                          # LSF 集群
```

MPICH/Intel MPI 用 `mpiexec -n 4 -genv OMP_NUM_THREADS=1 ...` 代替 `mpirun -np 4 -x ...`。

## 预期结果 / 如何检查

```bash
grep "SCF run converged" *.out                       # Gamma/OT 版本；k 点对角化版本同样打印
grep "ENERGY| Total FORCE_EVAL" *.out                # 总能量（Hartree）
grep -A10 "FORCES| Atomic forces" rutile_k446.out    # 2025.1 起的力输出格式
grep -A10 "ATOMIC FORCES" rutile_k446.out            # 2024.x 及以前
grep "BRILLOUIN|" rutile_k446.out                    # 确认 k 点方案与网格（如 K-Point grid 4 4 6）
```

判断要点：

- **每 TiO2 单元的能量**：原胞含 2 个 TiO2，`E/2`；超胞 `E/(2·NX·NY·NZ)`。Gamma 点原胞结果与 k 点结果相差很大（远大于 1 meV/原子），说明 Gamma 点不可用；k 点网格从 3×3×3 → 4×4×6 → 5×5×8 的变化应迅速减小。
- **原子力**：Ti 位于反演中心，受力由对称性严格为零；O 位于 $(u,u,0)$ 等位置，沿 [110] 方向的力不受对称性约束，大小取决于 $u$ 是否等于该方法下的平衡值。k 点收敛后，O 的残余力应很小（$u=0.305$ 接近实验值，也接近 PBE 平衡值）；Gamma 点原胞计算的 O 受力会明显偏离，这正是采样不足的表现。各原子力之和（`Sum`/`SUM OF ATOMIC FORCES`）应接近 0。
- **等价性检验**：`-E NX=3 -E NY=3 -E NZ=3` 的超胞（162 原子，Gamma 点）与 `-E KX=3 -E KY=3 -E KZ=3` 的 k 点计算（奇数网格包含 Gamma 点）采样的是同一组 k 点，每 TiO2 单元能量应非常接近（差异来自实空间网格的细微不同，通常远小于 1 meV/单元）。

（以上为定性预期，未附实际数值；请以你的运行结果为准。）

## 注意事项

- `input_file.inp` 中 `SCF_GUESS RESTART`：第一次运行时 `Universality-RESTART.wfn` 不存在，CP2K 会给出警告 “... This file does not exist ... Calculation continues using ATOMIC GUESS.” 并改用原子初猜，属正常现象；重复运行时会读入上一次的波函数。改进版输入改为明确的 `SCF_GUESS ATOMIC`。
- `@INCLUDE rutile.xyz` 把文件内容原样插入 `&COORD` 段，所以 `rutile.xyz` 不能带 xyz 表头。若要用标准 xyz 文件，请改用 `&TOPOLOGY COORD_FILE_NAME ... COORD_FILE_FORMAT XYZ`。
- k 点计算在 CP2K ≤ 2026.2 中**只能用对角化**（OT 不支持 k 点）；k 点重启文件是 `<PROJECT>-RESTART.kp`，与 Gamma 点的 `-RESTART.wfn` 不能互换。
- 晶格常数 4.593/2.959 Å 接近室温实验值（a ≈ 4.594 Å，c ≈ 2.959 Å，Burdett et al., *J. Am. Chem. Soc.* 109, 3639 (1987)，DOI: [10.1021/ja00246a021](https://doi.org/10.1021/ja00246a021)），并非 PBE 平衡值；要得到 PBE 晶格常数需做晶胞优化，见 [第7章](../../tutorials/07-geo-cell-opt.md)。
- `CUTOFF 400 Ry` 与 k 点网格都需要收敛测试，见 [第4章](../../tutorials/04-cutoff-convergence.md) 和 [第6章](../../tutorials/06-periodic-surface-kpoints.md)。

返回教程：[第2章 输入文件结构与第一个计算](../../tutorials/02-input-structure.md)
