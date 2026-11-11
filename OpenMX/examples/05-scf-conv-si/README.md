# 示例 05-1：Si 的 SCF 参数、截断能与 k 点收敛测试，以及 `scf.restart pp`

## 目的

- 给出一个把所有 SCF 关键词都显式写出的参考输入 `Si_scf.dat`；
- 用模板 + 驱动脚本完成 `scf.energycutoff` 和 `scf.Kgrid` 的收敛测试，按 1 meV/atom 判据选参数；
- 演示 OpenMX 4.0 的 `scf.restart pp`：读入自洽结果，只做一次对角化来计算能带和 DOS 数据。

对应教程：[第5章 SCF 参数、k 点与收敛测试](../../tutorials/05-scf-kpoints-convergence.md)（§8、§9）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Si_scf.dat` | Si 原胞（a = 5.43 Å）自洽计算，`System.Name si_scf`，Standard 基组，220 Ry，8×8×8，RMM-DIISK |
| `Si_restart_pp.dat` | 同一 `System.Name` 的后处理输入：`scf.restart pp` + 能带 + `Dos.fileout on`（`Dos.Kgrid 16 16 16`） |
| `conv_ecut/Si_ecut_template.dat` | 截断能扫描模板，占位符 `@ECUT@`（k 点固定 8×8×8） |
| `conv_kgrid/Si_k_template.dat` | k 点扫描模板，占位符 `@K@`（截断能固定 220 Ry，可按截断能测试结果修改） |
| `run_conv.sh` | 驱动脚本：`ecut` 或 `kgrid` 模式，`run`（默认）或 `collect` |
| `job.slurm` | Slurm 作业脚本 |

## 运行步骤

```bash
cd OpenMX/examples/05-scf-conv-si
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat conv_*/*.dat   # once

# A) cutoff-energy test, then k-point test
MPI_RUN="mpirun -np 8" bash run_conv.sh ecut
MPI_RUN="mpirun -np 8" bash run_conv.sh kgrid
#    on a Slurm cluster:
sbatch job.slurm run_conv.sh ecut
#    custom lists / tolerance:
ECUTS="150 200 250 300 400" THRESH=0.5 bash run_conv.sh ecut collect

# B) SCF, then post-processing restart
mpirun -np 8 openmx Si_scf.dat -nt 1 > si_scf.std
cp si_scf.out si_scf_step1.out                  # pp run will overwrite si_scf.out
mpirun -np 8 openmx Si_restart_pp.dat -nt 1 > si_pp.std
```

## 预期结果 / 如何检查

**收敛测试**：`run_conv.sh` 在 `conv_ecut/ecut_conv.txt`（或 `conv_kgrid/kgrid_conv.txt`）中给出表格：截断能（或 k 网格）、Utot、相对列表最后一项的 ΔE（meV/atom）、耗时、网格点数、`scf_conv`、是否在阈值内。

- 截断能扫描中，ΔE 会在 0 附近振荡而不是单调下降，这是实空间积分网格的特性；
- 找到"从这一项开始，后面所有项都 `ok=yes`"的最小值作为推荐值；
- 对 Si，k 点扫描通常在 8×8×8 附近进入 1 meV/atom 以内（以你的实际结果为准）；
- `scf_conv` 列必须全部为 1。

**restart pp**：

```bash
grep -c "NormRD" si_pp.std                        # only one SCF step is expected
ls si_scf.Band si_scf.Dos.val si_scf.Dos.vec      # band and DOS data
ls -l si_scf_rst/                                 # timestamps unchanged by the pp run
```

能带后处理（bandgnu13 或 `python ../../scripts/omx_band_plot.py -h`）见第8章，DOS（DosMain 或 `python ../../scripts/omx_dos_auto.py -h`）见第9章。

## 注意事项

- `run_conv.sh` 需要 GNU sed 和 awk；它从 `.out` 读取 `Utot.`、`Elapsed.Time.` 和 `Num. of grids of a-, b-, and c-axes` 三行，从 `.md` 读取 `scf_conv=`。
- 截断能扫描时 k 点固定，k 点扫描时截断能固定。k 点模板中的 `scf.energycutoff 220.0` 请改成你在截断能测试中选定的值。
- `scf.restart pp` 是 OpenMX 4.0 的功能。3.9 用户可以用 `scf.restart on` 加较小的 `scf.maxIter` 代替，但重启文件会被覆盖。
- `pp` 运行时结构、基组、`scf.XcType`、`scf.SpinPolarization` 必须与自洽计算一致，`System.Name` 必须相同。
- 输入文件只含 ASCII 字符和 LF 换行。
