# 示例 18：报错排查练习、版本迁移与并行性能测试

配套教程：第18章 常见错误排查、性能调优与版本迁移

## 目的

1. **排查练习**：`broken_scf.inp` 故意放了 6 个常见错误（A–F），按报错顺序逐个修正，最终得到 `broken_scf_fixed.inp`；
2. **版本迁移**：一对可 `diff` 的输入（2025 写法 / 2026.2 写法）、一个按文本扫描旧语法的脚本，以及一份迁移清单；
3. **性能调优**：在一个整节点上比较多种 MPI×OpenMP 组合的耗时。

## 文件说明

| 文件 | 说明 |
|---|---|
| `broken_scf.inp` | 故意出错的 fcc Cu（4 原子惯用胞，k 点 8×8×8）输入，错误用 (A)–(F) 标注 |
| `broken_scf_fixed.inp` | 修正版：UZH q11 基组/赝势、对角化 + Broyden + Fermi–Dirac 展宽 + `ADDED_MOS`；DOS 输出用 `NEWDOS` 变量在新旧语法间切换 |
| `legacy_si_cellopt_2025.inp` | Si 8 原子晶胞优化，2025 写法（`CELL_OPT TYPE`、`&PRINT &PDOS`），在 2026.2 上会报错 |
| `si_cellopt_2026_2.inp` | 同一计算的 2026.2 写法 |
| `check_migration.sh` | 文本扫描：找出 2026.2 删除/改动的关键字（不展开 `@INCLUDE`，`@IF` 块内的发现只作 INFO） |
| `migrate_2025_to_2026_2.md` | 2025.x → 2026.2 迁移清单（输入、构建、Python 工具） |
| `job_mpi_omp_hybrid.sh` | Slurm 整节点基准：依次运行若干 MPI×OpenMP 组合，汇总到 `bench_summary.dat` |

## 运行步骤

### 1. 排查练习

```bash
cd CP2K/examples/18-troubleshooting
export CP2K_DATA_DIR=/path/to/cp2k/data

cp2k.psmp -c broken_scf.inp                               # 第 1 轮：语法检查
bash check_migration.sh broken_scf.inp                    # 文本扫描（与 CP2K 版本无关）
# 修一处、跑一次：
cp2k.psmp -i broken_scf.inp -o broken_scf.out; tail -25 broken_scf.out
```

在 CP2K 2026.2 上，错误按下面的顺序出现（每修好一个，下一个才会暴露）：

| 顺序 | 标记 | 阶段 | 典型报错（节选） | 修正 |
|---|---|---|---|---|
| 1 | A | 解析输入 | `found an unknown keyword TYPE in section CELL_OPT` | 删除 `TYPE`；本例是单点计算，整个 `&MOTION` 都可以删 |
| 2 | B | 解析输入 | `unknown subsection PDOS of section PRINT` | 改用 `&PRINT &DOS &PDOS`（2026.2）或删除 |
| 3 | C | 读取基组 | `The requested basis set <DZVP-MOLOPT-SR-GTH> for element <Cu> was not found in the basis set files BASIS_MOLOPT_UZH` | 改成文件中真实存在的名字 `DZVP-MOLOPT-PBE-GTH-q11` |
| 4 | D | SCF 初始化 | `OT not possible with kpoint calculations` | 删除 `&OT`，改用对角化 |
| 5 | E | 无报错 | ——（CP2K 不会按名字里的 `-qN` 去核对基组与赝势） | 把 `GTH-PBE-q1` 改成 `GTH-PBE-q11`，并在输出中核对电子数 |
| 6 | F | SCF 迭代 | 很可能出现 `SCF run NOT converged`（或收敛到错误的整数占据态） | 加 `&SMEAR`（显式 `METHOD`）、`ADDED_MOS`、Broyden 混合、`MAX_SCF 100` |

在 2026.1 及更早版本上，A 和 B 不是错误，练习从 C 开始。

修正版的运行：

```bash
cp2k.psmp -i broken_scf_fixed.inp -o fixed.out                     # CP2K >= 2026.2
cp2k.psmp -E NEWDOS=0 -i broken_scf_fixed.inp -o fixed.out         # CP2K <= 2026.1
grep "Number of electrons" fixed.out | head -2                       # q11: 4 x 11 = 44
grep -E "SCF run converged|Fermi energy" fixed.out | tail -2
ls cu_fixed*.pdos cu_fixed-RESTART.kp
```

### 2. 版本迁移

```bash
diff legacy_si_cellopt_2025.inp si_cellopt_2026_2.inp
bash check_migration.sh legacy_si_cellopt_2025.inp si_cellopt_2026_2.inp
cp2k.psmp -c si_cellopt_2026_2.inp                     # 2026.2
cp2k.psmp -c legacy_si_cellopt_2025.inp                # 2025.x / 2026.1
```

### 3. MPI × OpenMP 基准

```bash
sbatch job_mpi_omp_hybrid.sh                                   # 默认用 broken_scf_fixed.inp，自动按节点核数生成组合
sbatch job_mpi_omp_hybrid.sh si_cellopt_2026_2.inp "32x1 16x2 8x4"
cat bench_summary.dat
```

## 预期结果 / 如何检查

- 修正版：输出中 `Number of electrons` 为 44，SCF 在 `MAX_SCF` 以内收敛，打印 `Fermi energy` 与 `Electronic entropic energy`；生成 `cu_fixed-RESTART.kp`（k 点计算写 `.kp`）；
- `check_migration.sh legacy_si_cellopt_2025.inp` 报告两个 ERROR（第 14 行 `TYPE`、第 57 行 `&PDOS`），对 `si_cellopt_2026_2.inp` 报告 0 个；
- `bench_summary.dat` 每行一个组合；`status` 一列若出现 `MPI_ranks_reported=1` 之类的标记，说明 MPI 环境错配（多个串行进程各自运行），见教程第 4 节。

本仓库没有附带参考输出；不同组合的快慢取决于硬件、MPI 库、编译选项和体系大小，没有通用的最优值。

## 注意事项

- `broken_scf.inp` 只用于练习，不要拿它的结果做任何分析；
- `check_migration.sh` 是基于文本的启发式检查，可能漏报或误报；以目标版本的 `cp2k.psmp -c` 为准；
- 基准测试请用与生产计算同体系、同参数的短任务（几个 SCF 或 MD 步），不要用玩具分子；记得把输入引用的坐标文件一起复制到 `bench_*` 目录（脚本中有注释说明）；
- Cu 的晶格常数 3.615 Å 是室温实验值的近似，仅作示例；本例的 CUTOFF、k 点和展宽温度未经收敛测试。
