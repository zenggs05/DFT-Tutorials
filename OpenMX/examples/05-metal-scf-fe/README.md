# 示例 05-2：bcc Fe 的金属 SCF 设置、电子温度与 k 点扫描

## 目的

- 给出一套适合磁性金属的 SCF 起步参数（RMM-DIISK、小混合权重、长历史、较高电子温度）；
- 扫描电子温度和 k 点，观察总能与磁矩的收敛行为；
- 练习运行中通过 `fe_bcc_SCF_keywords` 文件在线修改混合参数或终止计算。

对应教程：[第5章 SCF 参数、k 点与收敛测试](../../tutorials/05-scf-kpoints-convergence.md)（§7、§10）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Fe_bcc_metal_scf.dat` | bcc Fe 原胞（a = 2.87 Å），铁磁共线 GGA-PBE；`Fe_PBE19H` + `Fe5.5H-s3p2d1`，初始占据 `9.0 7.0`，300 Ry，16×16×16，600 K |
| `fe_bcc_SCF_keywords.example` | 在线调参文件模板；运行中复制为 `fe_bcc_SCF_keywords` 即生效 |
| `run_scan.sh` | `temp` 模式扫描电子温度，`kgrid` 模式扫描 k 网格，并汇总 Utot、磁矩、耗时、`scf_conv` |
| `job.slurm` | Slurm 作业脚本 |

## 运行步骤

```bash
cd OpenMX/examples/05-metal-scf-fe
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat   # once

# single SCF
mpirun -np 8 openmx Fe_bcc_metal_scf.dat -nt 1 | tee fe_bcc.std

# while it is running (other terminal): slow down the mixing on the fly
cp fe_bcc_SCF_keywords.example fe_bcc_SCF_keywords
grep "renewed" fe_bcc.std          # confirms that the file was read
rm fe_bcc_SCF_keywords             # remove it afterwards!

# scans
TEMPS="300 600 1000" MPI_RUN="mpirun -np 8" bash run_scan.sh temp
KGRIDS="10 12 14 16 18 20 24" MPI_RUN="mpirun -np 8" bash run_scan.sh kgrid
sbatch job.slurm run_scan.sh kgrid   # or on a Slurm cluster
```

## 预期结果 / 如何检查

```bash
grep -o "scf_conv=[01]" fe_bcc.md                # 1 = converged
grep -E "Utot\.|Chemical potential" fe_bcc.out
grep "Total spin moment" fe_bcc.out              # spin moment of the cell (muB)
grep "Total spin S" fe_bcc.out                   # S; moment = 2S muB
grep -A60 "SCF history at MD" fe_bcc.out | head -40
```

- 磁矩应收敛到铁磁态。bcc Fe 的实验磁矩约 2.2 μB/atom，PBE 一般给出 2.2–2.3 μB 量级（作为参考，不是本例的实际输出）；
- `NormRD` 应总体下降多个数量级，金属体系偶有反弹属正常；
- `run_scan.sh` 输出的表中，ΔE 相对列表最后一项（最高温度或最密 k 网格）；温度越低，需要越密的 k 网格才能使能量和磁矩稳定；
- 在线调参生效时，标准输出会出现 `The keywords for SCF iteration are renewed by ./fe_bcc_SCF_keywords.`。

## 注意事项

- 初始占据之和必须等于 `Fe_PBE19H` 的价电子数 16；写成 `8.0 8.0` 可能收敛到非磁态。
- `Fe_PBE19H` 含 3s/3p 半芯态，所以截断能取 300 Ry；用 `Fe_PBE19S`（14 个价电子，PAO `Fe6.0S-s3p2d1`）时要同时改坐标行的初始占据（如 `8.0 6.0`）。
- 在线调参文件中可以修改的只有 `scf.maxIter`、`scf.Min.Mixing.Weight`、`scf.Max.Mixing.Weight`、`scf.Kerker.factor`、`scf.Mixing.StartPulay`；把 `scf.maxIter` 写成负数会让程序停止。
- 能量差（例如不同晶格常数、FM 与 AFM）必须在相同电子温度和相同 k 网格下比较。
- 输入文件只含 ASCII 字符和 LF 换行。
