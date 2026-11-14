# 示例 02-job-templates：Slurm/PBS/LSF 作业脚本模板（纯 MPI、MPI×OpenMP、共享内存、作业数组、续算）

## 目的

- 提供一组**与具体软件解耦**的作业脚本模板，覆盖三种并行方式：纯 MPI（VASP、OpenMX、CP2K）、MPI×OpenMP 混合（cp2k.psmp、openmx -nt）、单节点共享内存（Gaussian 16）；
- 给出同一个 MPI 作业在 Slurm、PBS Pro/Torque、LSF 三种调度器下的写法，方便对照迁移；
- 演示几个实用技巧：提交前检查 CRLF 换行和 NCORE/KPAR 是否与核数匹配、每个作业独立的本地 scratch、临近时限时用信号让程序体面退出并自动续投、用 `sacct` 汇总作业状态。

对应教程：[第2章 HPC 作业调度：Slurm/PBS/LSF 模板与并行设置](../../tutorials/02-job-scheduler.md)。

## 文件说明

| 文件 | 调度器 | 并行方式 | 典型用途 |
|---|---|---|---|
| [slurm_mpi.sh](slurm_mpi.sh) | Slurm | 纯 MPI | 通用模板：`EXE` 可设为 `vasp_std`/`cp2k.psmp`/`openmx`；含 CRLF 检查、`--signal` + `trap` 示例 |
| [slurm_vasp.sh](slurm_vasp.sh) | Slurm | 纯 MPI | VASP 专用：`VARIANT=std/gam/ncl`，检查四个输入文件，检查 NCORE 是否整除每节点进程数、KPAR 是否整除总进程数，结束后判断 OUTCAR 是否完整 |
| [slurm_resubmit_vasp.sh](slurm_resubmit_vasp.sh) | Slurm | 纯 MPI | 长时间弛豫/MD：时限前 15 分钟写 `STOPCAR`（`LSTOP = .TRUE.`），CONTCAR→POSCAR 后自动续投，最多 `MAX_ROUND` 轮 |
| [slurm_hybrid_mpi_omp.sh](slurm_hybrid_mpi_omp.sh) | Slurm | MPI×OpenMP | CP2K（`cp2k.psmp`）、OpenMX（`-nt`）；含 `OMP_PLACES`/`OMP_PROC_BIND`/`OMP_STACKSIZE` 和 Open MPI、Intel MPI 的等价写法 |
| [slurm_shared_memory.sh](slurm_shared_memory.sh) | Slurm | 单节点多线程 | Gaussian 16：`GAUSS_SCRDIR` 指向本地盘上每个作业独立的目录，退出时自动清理；检查 `%NProcShared` 与 `--cpus-per-task` 是否一致 |
| [slurm_array.sh](slurm_array.sh) | Slurm | 作业数组 | 通用参数扫描：`template/` 中含占位符 `@VALUE@`，每个数组任务建一个 `run_<值>/` |
| [pbs_mpi.sh](pbs_mpi.sh) | PBS Pro / Torque | 纯 MPI | 默认 `EXE=vasp_std`（即 PBS 版 VASP 脚本）；`select=...:mpiprocs=...` 语法，注释中给出 Torque 的 `nodes=:ppn=` |
| [lsf_mpi.sh](lsf_mpi.sh) | LSF | 纯 MPI | 默认 `EXE=vasp_std`（即 LSF 版 VASP 脚本）；`-n` 总槽数 + `span[ptile=]`，注意用 `bsub < lsf_mpi.sh` 提交 |
| [sacct_report.sh](sacct_report.sh) | Slurm | — | 在登录节点汇总最近 N 天作业的状态、耗时、内存峰值，并列出失败/超时/OOM 的作业 |

各软件更完整的专用脚本（本仓库其他板块）：

| 软件 | 脚本 |
|---|---|
| VASP | [VASP/scripts/job_vasp_slurm.sh](../../../VASP/scripts/job_vasp_slurm.sh) |
| Gaussian 16 | [Gaussian/scripts/sub_g16.sh](../../../Gaussian/scripts/sub_g16.sh) |
| CP2K | [CP2K/scripts/job_cp2k_slurm.sh](../../../CP2K/scripts/job_cp2k_slurm.sh) |
| OpenMX | [OpenMX/scripts/job_omx_slurm.sh](../../../OpenMX/scripts/job_omx_slurm.sh) |

另外，Common 板块的脚本 [mkjob.py](../../scripts/mkjob.py) 可以按调度器和软件生成作业脚本，[job_status.py](../../scripts/job_status.py) 可以批量检查目录中各计算是否正常结束，用 `python ../../scripts/mkjob.py -h`、`python ../../scripts/job_status.py -h` 查看用法。

## 运行步骤

```bash
# 1) 复制到集群的工作盘（不要放在 /home 下跑大计算）
rsync -av Common/examples/02-job-templates/ hpc:~/work/templates/

# 2) 在集群上查清楚本机的分区、每节点核数、内存
sinfo -s                                  # 分区列表
sinfo -N -o "%N %c %m %P" | head          # 节点名、CPU 数、内存(MB)、分区
lscpu | grep -E "^CPU\(s\)|Thread|Socket|NUMA"

# 3) 按上面的结果修改模板头部的 --partition、--ntasks-per-node、--mem、module load
vi slurm_vasp.sh

# 4) 静态检查：语法、换行符、调度器能否接受（--test-only 不会真正提交）
bash -n slurm_vasp.sh
file slurm_vasp.sh                        # 不应出现 "with CRLF line terminators"
sbatch --test-only slurm_vasp.sh

# 5) 放到算例目录中提交
cp slurm_vasp.sh ~/work/si_scf/ && cd ~/work/si_scf && sbatch slurm_vasp.sh
squeue -u $USER
```

PBS 和 LSF 的提交命令分别是 `qsub pbs_mpi.sh` 和 `bsub < lsf_mpi.sh`。

## 预期结果 / 如何检查

- `sbatch --test-only` 输出形如 `sbatch: Job 123456 to start at ... using 64 processors on nodes ... in partition cpu`，说明资源请求合法；若报 `Requested node configuration is not available`，说明每节点核数或内存超过了该分区节点的实际配置。
- 作业运行后，`%x-%j.out` 的开头会打印节点名、进程数、线程数；VASP 模板还会打印 `ranks total=... per node=... NCORE=... KPAR=...`，若出现 `WARNING: NCORE=... does not divide ...`，应调整 INCAR。
- `slurm_vasp.sh` 正常结束时打印最后一个 `free  energy   TOTEN` 行；OUTCAR 没有 `General timing and accounting` 时脚本以非零退出码结束，`sacct` 中显示为 `FAILED`，可以让依赖链（`afterok`）自动停下。
- `slurm_resubmit_vasp.sh`：每轮的输出保存为 `vasp.out.1`、`vasp.out.2`……，上一轮的 OUTCAR 改名为 `OUTCAR.1` 等；收敛后打印 `converged in round N`。
- `bash sacct_report.sh 7` 列出最近 7 天的作业和各状态的计数。

## 注意事项

- `--partition`、`--account`、`--qos`、每节点核数、`--mem=0` 是否允许、本地 scratch 的路径（`$TMPDIR`、`/scratch/$USER`、`$SLURM_TMPDIR` 等）都因集群而异，以本集群的用户手册为准。
- 模板中所有 `module load` 都是注释掉的示例；必须换成**编译该程序时**所用的编译器、MPI、数学库模块，否则常见 `error while loading shared libraries` 或运行时段错误。
- 用 `srun` 启动 MPI 程序要求 MPI 库带有 Slurm 的 PMI/PMIx 支持；若 `srun` 启动后出现“每个进程都以为自己是 rank 0”（输出重复 N 次），改用 `mpirun` 或咨询管理员应使用的 `--mpi=` 类型。
- 在 Windows 上编辑过的脚本先 `dos2unix`，否则 `sbatch: error: Batch script contains DOS line breaks (\r\n)`。
- `slurm_resubmit_vasp.sh` 需要计算节点能执行 `sbatch`；有些集群禁止，此时去掉自动续投，改为手动 `sbatch`。
- VASP 的 POTCAR 受许可证保护，不在本仓库中，需用自己的赝势库生成。
