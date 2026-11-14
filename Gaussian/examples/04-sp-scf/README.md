# 示例 04：单点能计算与 SCF 收敛技巧（Gaussian 16）

配套教程：[第4章 单点能计算与 SCF 收敛技巧](../../tutorials/04-sp-scf.md)

## 目的

通过 6 个小体系练习 Gaussian 16 单点能计算，以及常用的 SCF 收敛手段：

- 读懂 `#p` 下的 SCF 迭代输出和 `SCF Done` 能量；
- 过渡金属高自旋配合物用 `SCF=XQC` + `NoSymm`；
- 先用小基组收敛、再 `Guess=Read` 读入大基组（`--Link1--` 两步作业）；
- 小能隙开壳层团簇用 `SCF=VShift`；
- 用 `Guess=Fragment` 构造 broken-symmetry 初猜；
- 用 `Stable=Opt` 发现并修复 RHF→UHF 不稳定性。

## 文件说明

| 文件 | 体系 | 电荷/多重度 | 演示内容 |
|---|---|---|---|
| `h2o_sp_b3lyp_def2tzvp.gjf` | H2O | 0 / 1 | 最基本的单点能，`#p` 查看 SCF 迭代 |
| `fe_complex_sp_scf_xqc.gjf` | [Fe(H2O)6]²⁺（高自旋 d⁶） | +2 / 5 | `SCF=(XQC,MaxConventionalCycles=64)` + `NoSymm` |
| `guess_read_small_to_big_basis.gjf` | [Cu(NH3)4]²⁺（d⁹） | +2 / 2 | 第 1 步 def2-SVP，第 2 步 def2-TZVP `Geom=Check Guess=Read` |
| `radical_sp_vshift.gjf` | Cu3 团簇 | 0 / 2 | `SCF=(VShift=500,MaxCycle=300)` |
| `h2_stretched_guess_fragment.gjf` | 拉长的 H2（2.0 Å） | 0 / 1 | `Guess=Fragment=2`，碎片多重度 `0 2 0 -2` |
| `h2_stretched_stable_opt.gjf` | 拉长的 H2（2.0 Å） | 0 / 1 | `Stable=Opt`：RB3LYP → UB3LYP |
| `run_g16.sh` | — | — | 最简 Slurm 作业脚本（8 核、16 GB，与输入中的 `%nprocshared=8`、`%mem=12GB` 对应） |

结构说明：H2O 为常用实验近似结构；配合物和 Cu3 的坐标按典型键长手工搭建（Fe–O 2.13 Å、Cu–N 2.03 Å、Cu–Cu 2.30/2.45 Å），只用于演示单点和 SCF 选项，**不是**优化后的结构。

## 运行步骤

本地（已配置好 `g16`）：

```bash
cd Gaussian/examples/04-sp-scf
for f in h2o_sp_b3lyp_def2tzvp h2_stretched_stable_opt h2_stretched_guess_fragment; do
    g16 < $f.gjf > $f.log
done
g16 < guess_read_small_to_big_basis.gjf > guess_read_small_to_big_basis.log
g16 < fe_complex_sp_scf_xqc.gjf > fe_complex_sp_scf_xqc.log
g16 < radical_sp_vshift.gjf > radical_sp_vshift.log
```

Slurm 集群：

```bash
sbatch run_g16.sh fe_complex_sp_scf_xqc.gjf
sbatch run_g16.sh guess_read_small_to_big_basis.gjf
```

`run_g16.sh` 中的 `module load` 或 `g16.profile` 两行需要按本集群的实际情况修改；scratch 目录默认放在 `/tmp/$USER` 下，可以用环境变量 `SCR_BASE` 改到更大的本地盘。更完整的提交脚本是 [`sub_g16.sh`](../../scripts/sub_g16.sh)（`bash ../../scripts/sub_g16.sh -h` 查看用法），作业调度的基础知识见 [Common 第2章 HPC 作业调度](../../../Common/tutorials/02-job-scheduler.md)。

H2O 和 H2 只需几秒；[Cu(NH3)4]²⁺ 和 [Fe(H2O)6]²⁺ 在 8 核上通常为几分钟到几十分钟，取决于机器。

## 预期结果 / 如何检查

```bash
# 1) every job terminated normally (the Link1 job prints it twice)
grep -c "Normal termination" *.log

# 2) energies; reject any SCF Done that follows "Convergence criterion not met"
grep "SCF Done" *.log
grep -l "Convergence criterion not met" *.log      # should print nothing

# 3) spin contamination of the open-shell jobs
grep "S\*\*2 before annihilation" fe_complex_sp_scf_xqc.log guess_read_small_to_big_basis.log radical_sp_vshift.log

# 4) Guess=Read really used in step 2
grep "Initial guess from the checkpoint file" guess_read_small_to_big_basis.log

# 5) stability analysis of stretched H2
grep -E "instability|wavefunction is stable" h2_stretched_stable_opt.log
```

判断标准：

| 作业 | 应该看到 |
|---|---|
| H2O | 1 个 `Normal termination`；`SCF Done: E(RB3LYP)` 约为 −76.4 Eh；`-V/T` 接近 2 |
| [Fe(H2O)6]²⁺ | `<S**2>` 接近 S(S+1) = 6.00 |
| [Cu(NH3)4]²⁺ | 2 个 `Normal termination`；第 2 步出现 `Initial guess from the checkpoint file`，并且迭代轮数少于第 1 步；`<S**2>` 接近 0.75 |
| Cu3 | `<S**2>` 接近 0.75；`#p` 输出中 `Shift=` 字段不为 0 |
| H2（Stable=Opt） | 先报告存在 RHF→UHF 不稳定性，重新优化后能量更低，`<S**2>` 接近 1 |
| H2（Guess=Fragment） | 直接收敛到 broken-symmetry 解，能量和 `<S**2>` 与上一行的最终结果一致（差别在数值误差范围内） |

这里不给出具体能量数值，请以自己的输出为准；不同版本、不同机器之间末几位小数可能不同。

## 注意事项

- 所有输入都以空行结尾，修改时不要删掉最后的空行。
- 内存和核数：`%mem=12GB` 对应 Slurm `--mem=16G`（约 75%）；修改其中一个时要同步修改另一个。
- `radical_sp_vshift.gjf` 用了 VShift，Gaussian 会关闭自动 archive，所以日志末尾没有 `1\1\` 开头的归档块，属于正常现象。
- 配合物和团簇的自旋态能量对泛函很敏感，本目录只演示 SCF 技巧，不代表这些体系的推荐计算方案。
- 如果作业仍然 L502 不收敛，按教程第 5 节的流程排查，也可以用 [`g16_error_doctor.py`](../../scripts/g16_error_doctor.py) 自动识别错误（`python ../../scripts/g16_error_doctor.py -h`）。
