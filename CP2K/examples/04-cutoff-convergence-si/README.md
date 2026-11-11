# 示例 04：体相 Si 与 H2O 的 CUTOFF / REL_CUTOFF 收敛测试

配套教程：[第4章 CUTOFF / REL_CUTOFF 收敛测试](../../tutorials/04-cutoff-convergence.md)

## 目的

用“两步法”为 PBE + `DZVP-MOLOPT-PBE-GTH`（BASIS_MOLOPT_UZH）+ `GTH-PBE`（POTENTIAL_UZH）组合确定 `&MGRID` 的 `CUTOFF` 和 `REL_CUTOFF`：

1. 固定 `REL_CUTOFF = 60 Ry`，扫描 `CUTOFF`（默认 150–600 Ry）；
2. 固定第一步选出的 `CUTOFF`，扫描 `REL_CUTOFF`（默认 20–100 Ry）。

主体系是 8 原子立方胞体相 Si（a = 5.431 Å，实验值），附带一个孤立 H2O 分子模板，用于比较分子体系所需的更高截断。

## 文件说明

| 文件 | 说明 |
|---|---|
| `template.inp` | Si8 的参数化输入：顶部 `@SET CUTOFF` / `@SET REL_CUTOFF` 由脚本改写；对角化 + Broyden + Fermi–Dirac 300 K，`PRINT_LEVEL MEDIUM` 用于打印 `MULTIGRID INFO` |
| `si_bulk8.xyz` | Si8 坐标（Å） |
| `h2o_template.inp` | H2O（12 Å 盒子，`PERIODIC NONE` + MT 求解器，OT）的参数化输入 |
| `h2o.xyz` | H2O 坐标（实验气相几何） |
| `run_cutoff_scan.sh` | 第一步：生成 `cutoff_<值>Ry/` 目录、运行、调用 `collect_scan.sh` |
| `run_relcutoff_scan.sh` | 第二步：生成 `relcutoff_<值>Ry/` 目录、运行、汇总 |
| `collect_scan.sh` | 从 `output.out` 提取 `ENERGY| Total FORCE_EVAL`、`SCF run converged` 和 `count for grid`，写出 `cutoff_scan.csv` / `relcutoff_scan.csv` |
| `plot_scan.py` | 画 \|ΔE\|（meV/原子，对数坐标）随截断的变化，并报告从哪个值起稳定低于阈值（需要 numpy、matplotlib） |
| `job_scan_slurm.sh` | Slurm 作业：在一个节点上依次跑完两步扫描 |

脚本的可调参数都通过环境变量传入：`TEMPLATE`、`EXTRA_FILES`、`CUTOFFS`、`RELCUTOFFS`、`REL_CUTOFF`、`CUTOFF`、`CP2K`（可包含 `mpirun -np N`）、`NPAR`（同时运行的计算个数）、`DRYRUN=1`（只生成输入）、`NATOMS`（默认取 xyz 文件第一行）。

## 运行步骤

```bash
cd CP2K/examples/04-cutoff-convergence-si
export CP2K_DATA_DIR=/path/to/cp2k/data      # 若你的 CP2K 已内置数据路径可省略

# 0) 语法检查
cp2k.psmp -c template.inp

# 1) 第一步：扫 CUTOFF（REL_CUTOFF = 60）
bash run_cutoff_scan.sh
#    多核并行示例：每个计算 4 个 MPI 进程，同时跑 2 个
#    CP2K="mpirun -np 4 cp2k.psmp" NPAR=2 bash run_cutoff_scan.sh
python plot_scan.py cutoff_scan.csv --threshold 0.03 -o cutoff_scan.png

# 2) 第二步：用第一步选出的值（下面以 350 Ry 为例）扫 REL_CUTOFF
CUTOFF=350 bash run_relcutoff_scan.sh
python plot_scan.py relcutoff_scan.csv --threshold 0.03 -o relcutoff_scan.png

# 集群上：一个作业跑完两步
sbatch job_scan_slurm.sh
```

H2O 分子请放在单独的子目录中运行，以免与 Si 的结果目录混在一起：

```bash
mkdir -p h2o && cp h2o_template.inp h2o.xyz run_*.sh collect_scan.sh plot_scan.py h2o/ && cd h2o
TEMPLATE=h2o_template.inp EXTRA_FILES=h2o.xyz \
  CUTOFFS="300 400 500 600 700 800 1000 1200" bash run_cutoff_scan.sh
```

## 预期结果 / 如何检查

- 每个子目录的 `output.out` 中都应出现 `*** SCF run converged in ... steps ***` 和 `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]`；
- `cutoff_scan.csv` 中 `dE_meV_per_atom` 随截断升高逐渐趋于 0，但**不一定单调**（GPW 对网格不满足变分原理）；
- `grid_counts` 列（`grid1;grid2;grid3;grid4`）随 CUTOFF 升高时，第 1 级计数逐渐减少；随 REL_CUTOFF 升高时，计数向细网格移动；
- 选择标准：从某个值起所有点都低于阈值（严格阈值 0.03 meV/原子 ≈ 10⁻⁶ Ha/原子；常规 1 meV/原子），取该值再加约 50 Ry 余量，并确认第 1 级计数不为 0；
- H2O 所需的 CUTOFF 明显高于 Si（O 的基组指数和赝势更硬）。

本仓库没有附带参考输出；具体数值取决于 CP2K 版本、编译选项和并行设置，请以你自己的运行结果为准。

## 注意事项

- 模板中 `EPS_DEFAULT 1.0E-12`、`EPS_SCF 1.0E-8` 比生产计算更严格，这是为了让 EPS 和 SCF 噪声远低于网格误差；确定截断后，生产输入可以恢复常规值；
- 输入中的基组名与赝势名必须与 `BASIS_MOLOPT_UZH`、`POTENTIAL_UZH` 中的条目完全一致，`-q4`、`-q6`、`-q1` 必须成对匹配；
- 2024.1 起 SCF 不收敛会直接中止，`collect_scan.sh` 会跳过没有最终能量的目录；如果低截断点普遍不收敛，从扫描范围中去掉这些点即可，不要在模板里加 `IGNORE_CONVERGENCE_FAILURE`；
- 若要测试应力收敛（CELL_OPT/NPT 用），把 `RUN_TYPE` 改为 `ENERGY_FORCE`（`ENERGY` 不计算力和应力），在 `&FORCE_EVAL` 中加入 `STRESS_TENSOR ANALYTICAL` 和 `&PRINT` / `&STRESS_TENSOR ON` 打印段，再用 `grep "STRESS| 1/3 Trace"` 提取压力，方法见教程第 7 节；
- 不同截断下得到的能量不能互相做差；正式计算中所有用来求能量差的结构必须使用同一组 `&MGRID` 参数。
