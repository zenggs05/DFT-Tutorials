# 示例 18-1：OpenMX 作业脚本（Slurm / PBS）、运行检查与 MPI×OpenMP 基准

配套教程：[第18章 并行、大体系 O(N)、作业脚本与常见错误排查](../../tutorials/18-parallel-on-errors.md)

## 目的

提供一套可以直接用于生产计算的 OpenMX 4.0 作业模板：

- 在 Slurm 和 PBS Pro/OpenPBS 上以 MPI × OpenMP 混合方式运行 `openmx`，保证 `-nt` 与 `OMP_NUM_THREADS` 一致，并设置 `ulimit -s unlimited` 和 `OMP_STACKSIZE`；
- 兼容官方 deb 包：设置 `LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib`；
- 计算结束后按 `.std` 末尾的 `The calculation was normally finished.` 判断是否成功，不依赖退出码；
- 批量检查已完成（或崩溃）的计算：是否正常结束、有无已知报错、SCF 步数与最后的 `NormRD`、`Utot`、耗时；
- 对同一个输入扫描不同的“进程数 × 线程数”组合，找出本机最快的配置。

## 文件说明

| 文件 | 说明 |
|---|---|
| `job_slurm.sh` | Slurm 模板：`sbatch job_slurm.sh input.dat`。默认 2 节点 × 32 进程 × 2 线程；`LAUNCH=srun` 改用 `srun`（显式传 `--cpus-per-task` 并设置 `SRUN_CPUS_PER_TASK`，兼容不继承该选项的 Slurm 版本），`USE_DEB=1` 用于 deb 包，`KEEP_SCRATCH=0` 时成功后删除 `_rst/` |
| `job_pbs.sh` | PBS Pro / OpenPBS 模板：`qsub -v INPUT=input.dat job_pbs.sh`。用 `select=...:mpiprocs=...:ompthreads=...` 申请资源；文件头注明了 Torque 的写法 |
| `check_run.sh` | 运行检查：`bash check_run.sh`（当前目录所有 `*.std`）或 `bash check_run.sh a.std b.std`；有问题时退出码为 1，方便放进流水线 |
| `bench_np_nt.sh` | 并行配置基准：`bash bench_np_nt.sh input.dat "8x1 4x2 2x4"`，每个组合用不同的 `System.Name` 运行，汇总到 `bench_<name>.csv` |

这几个脚本不依赖仓库中的其他文件，可以直接复制到自己的计算目录或放进 `PATH` 长期使用。Slurm/PBS 资源写法的通用说明见 [公共篇 作业调度](../../../Common/tutorials/02-job-scheduler.md)；仓库的通用模板 [job_omx_slurm.sh](../../scripts/job_omx_slurm.sh)、[job_omx_pbs.sh](../../scripts/job_omx_pbs.sh) 思路相同，用 `bash ../../scripts/job_omx_slurm.sh -h` 查看用法。

## 运行步骤

```bash
# 1. adapt the #SBATCH / #PBS header (partition, account, time) and the module lines
vi job_slurm.sh

# 2. submit
sbatch job_slurm.sh si512.dat                    # Slurm
qsub -v INPUT=si512.dat job_pbs.sh               # PBS Pro / OpenPBS
USE_DEB=1 sbatch job_slurm.sh si512.dat          # official .deb installation in /opt/openmx4.0

# 3. after the job
bash check_run.sh si512_dclno.std

# optional: find the fastest NP x NT split on one node (use a short test input, e.g. scf.maxIter 5)
bash bench_np_nt.sh test.dat "32x1 16x2 8x4 4x8"
```

## 预期结果 / 如何检查

- 作业日志（`openmx-<jobid>.log`）开头打印 `NP=... NT=...`，结尾打印 `finished normally` 以及最后的 `Utot.` 行和 `Total Computational Time` 行；失败时打印 `.std` 最后 30 行，并以退出码 1 结束。
- 用 `top -H` 或 `htop` 登录计算节点检查：每个 `openmx` 进程应占用约 `NT×100%` 的 CPU，节点上的总线程数不超过物理核数。
- `check_run.sh` 的输出（示意输出，数值仅供参考）：

```text
================ si512_dclno.std
status      : normally finished
SCF         : 23 steps in last block, last NormRD = 0.000000012345, |dUele| = 4.2e-09 Hartree
...
```

- `bench_np_nt.sh` 生成的 CSV 中，所有 `ok` 行的 `Utot_Hartree` 应一致（相差约 1e-8 Hartree 以内）；比较同一核数下 `total_time_s` 最小的组合。

## 注意事项

- **`-nt` 必须与申请的 `--cpus-per-task`（或 `ompthreads`）一致**，并同时导出 `OMP_NUM_THREADS`，否则 MKL/OpenBLAS 可能在每个进程里再开满线程，造成严重超订（oversubscription）。
- `ulimit -s unlimited` 只对当前 shell 及其子进程生效。多节点运行时，远端节点上的进程是否继承，取决于 MPI 的启动方式；若仍出现段错误，请咨询管理员把栈限制写进 `/etc/security/limits.conf` 或 Slurm 的 `PropagateResourceLimits` 配置。
- Intel MPI 用 `I_MPI_PIN_DOMAIN=omp` 绑定；Open MPI 用脚本注释中的 `--map-by ppr:N:node:PE=NT --bind-to core`。两者的选项不能混用。
- `bench_np_nt.sh` 解析的 `Diagonalization` 行来自 `.out` 末尾的耗时表，不同版本的行名可能略有差异，解析失败时该列为空，不影响总时间。
- 脚本里的 `module load` 行只是示例，请换成你所在集群的模块名。

返回：[第18章](../../tutorials/18-parallel-on-errors.md) | [OpenMX 目录](../../README.md)
