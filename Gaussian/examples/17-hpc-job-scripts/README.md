# 示例 17：Gaussian 16 作业脚本、并行、批量提交与断点续算

## 目的

提供一套可以直接修改使用的集群脚本，并用苯（C₆H₆）这个小分子演示：

1. Slurm、PBS、GPU、Linda 四种作业模板，资源只在调度指令中写一次，由脚本换算成 `GAUSS_CDEF`/`GAUSS_PDEF` 和 `GAUSS_MDEF`；
2. 批量提交整个目录、作业数组、依赖作业（优化成功后自动做高精度单点）；
3. 缩放测试：同一输入用不同核数运行，计算加速比和并行效率；
4. 断点续算：`Opt=Restart`，以及解析频率的 `%RWF` + `Restart`；
5. `Default.Route` 的写法。

对应教程：[第17章 作业脚本、并行与批量计算](../../tutorials/17-hpc-jobs.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `g16_slurm.sh` | **主模板**（Slurm，单节点 SMP）。读取 `Cpus_allowed_list` 设置 `GAUSS_CDEF`，按 `--mem` 的 80% 设置 `GAUSS_MDEF`；检查 CRLF、末尾空行和输入中的 Link 0 资源行；独立 scratch 并在退出时清理；结束后检查 `Normal termination`、对所有 `%chk` 做 formchk；失败时退出码为 1 |
| `g16_pbs.sh` | PBS/Torque 模板（注释中给出 PBS Pro 的 `select` 写法），用 `qsub -v INPUT=job.gjf` 传入输入文件 |
| `g16_gpu_slurm.sh` | GPU 模板：`--gres=gpu:2`，自动组成 `GAUSS_GDEF="0,1=<控制核>"`，打印 `nvidia-smi topo -m` 供核对 |
| `g16_linda_slurm.sh` | 两节点 Linda 模板：`GAUSS_WDEF` 取自 `$SLURM_JOB_NODELIST`，每个节点建立 scratch |
| `batch_submit_all_gjf.sh` | 批量提交目录中所有 .gjf，跳过已完成和已在队列中的作业，并限制排队数量 |
| `g16_array_slurm.sh` | 作业数组：第 i 个任务运行 `joblist.txt` 第 i 行的输入（调用 `g16_slurm.sh`） |
| `chain_opt_sp_dependency.sh` | 依赖作业：`benzene_opt_freq.gjf` 成功后才运行 `benzene_sp_tz.gjf` |
| `scaling_test.sh` | 缩放测试：在一次分配中依次用 1/2/4/8/16/32 核运行，输出 `scaling_<name>.dat` |
| `Default.Route.example` | `Default.Route` 示例：默认内存 4 GB、`MaxDisk=100GB`。复制为 `Default.Route` 放到作业目录才生效 |
| `benzene_opt_freq.gjf` | B3LYP-D3(BJ)/def2-SVP 优化 + 频率。**没有** `%nprocshared`/`%mem`，资源由脚本设置 |
| `benzene_sp_tz.gjf` | `%oldchk` 读取优化结构，ωB97X-D/def2-TZVP 单点（依赖作业的第二步） |
| `benzene_opt_restart.gjf` | `Opt=(Restart,MaxCycles=200)` 续算优化，再用 `--Link1--` 补做频率 |
| `benzene_freq_rwf.gjf` | 解析频率，`%rwf` 命名 read-write 文件 + `%nosave`，使作业可以续算 |
| `benzene_freq_restart.gjf` | 上一个作业中断后的续算输入：只有 Link 0 和 `#p Restart` |
| `benzene_scaling_freq.gjf` | 缩放测试用的 def2-TZVP 频率计算（只用于计时，结构未优化，热化学结果没有意义） |

## 运行步骤

先按你的集群修改每个脚本开头的 site specific 部分：`module load` 或 `g16root`、分区名（`--partition`）、scratch 位置（环境变量 `SCR_BASE`，默认 `/tmp`）。

```bash
cd Gaussian/examples/17-hpc-job-scripts

# 1) one job with the main template (resources can be overridden on the command line)
sbatch g16_slurm.sh benzene_opt_freq.gjf
sbatch --cpus-per-task=8 --mem=16G g16_slurm.sh benzene_opt_freq.gjf

# 2) PBS
qsub -v INPUT=benzene_opt_freq.gjf g16_pbs.sh

# 3) opt+freq, then the single point only if the first job succeeded
bash chain_opt_sp_dependency.sh

# 4) submit everything in a folder / as a job array
bash batch_submit_all_gjf.sh . 20
ls benzene_opt_freq.gjf benzene_scaling_freq.gjf > joblist.txt
sbatch --array=1-2%2 g16_array_slurm.sh

# 5) scaling test (one exclusive node, 32 cores)
sbatch scaling_test.sh benzene_scaling_freq.gjf "1 2 4 8 16 32"

# 6) restart demos
sbatch --time=00:02:00 g16_slurm.sh benzene_opt_freq.gjf   # let it time out on purpose
sbatch g16_slurm.sh benzene_opt_restart.gjf                # continue from benzene_opt_freq.chk
```

`benzene_freq_rwf.gjf` 的 rwf 写在作业目录中（相对路径），所以 `g16_slurm.sh` 清理 scratch 时不会删除它；作业中断后提交 `benzene_freq_restart.gjf` 即可续算。

## 预期结果 / 如何检查

| 检查项 | 命令 | 期望 |
|---|---|---|
| 资源是否按预期传给 Gaussian | `grep "GAUSS_.DEF" *.out` | `GAUSS_CDEF` 中的核数等于 `--cpus-per-task`，`GAUSS_MDEF` 约为 `--mem` 的 80% |
| 没有 Link 0 覆盖 | `grep WARNING *.out` | 无输出（苯的输入中没有资源行） |
| 作业成功 | `tail -n 1 benzene_opt_freq.log` | `Normal termination of Gaussian 16 ...` |
| fchk 自动生成 | `ls benzene_opt_freq.fchk` | 存在 |
| CPU 效率 | `seff <jobid>` | 小分子在 8–16 核时 CPU 效率明显低于 100% 是正常的 |
| 依赖作业 | `squeue -u $USER` | 单点作业状态为 `PD`，原因 `(Dependency)`，优化完成后才开始 |
| 缩放测试 | `cat scaling_benzene_scaling_freq.dat` | 苯这种小分子在少量核之后加速比明显变平 |
| 优化续算 | `grep -c "Normal termination" benzene_opt_restart.log` | 续算的优化和频率都正常结束 |

本仓库没有附带运行结果；墙钟时间和加速比取决于你的硬件。

## 注意事项

- 本目录的输入都**故意不写** `%nprocshared`/`%mem`。如果你的输入里有这些行，它们会覆盖脚本设置的环境变量，`g16_slurm.sh` 会打印警告。
- 共享节点上不要把 `%CPU=0-15` 写死在输入中：Slurm 可能分配给你的是其他编号的核。
- `--dependency=afterok` 依赖作业脚本的退出码；这里的模板只在 `Normal termination` 时返回 0。
- `scaling_test.sh` 使用 `--exclusive`，在不允许独占节点的集群上删除这一行，并把核数列表限制在申请的范围内。
- GPU 模板默认让前两个分配到的核控制两块 GPU，这只是最简单的选择。请根据日志中 `nvidia-smi topo -m` 的 CPU Affinity 把控制核改到与 GPU 相连的那颗 CPU 上。A100 需要 G16 Rev. C.02。
- Linda 模板需要 Linda 版本的 Gaussian 和节点间免密码 ssh；普通（SMP 版）Gaussian 不能跨节点。
- `Default.Route.example` 只有改名为 `Default.Route` 并放在作业目录（或 `$g16root/g16`）时才生效，它的优先级最低。
- 想要功能更完整的提交封装，可以使用仓库脚本 [`sub_g16.sh`](../../scripts/sub_g16.sh)（`bash ../../scripts/sub_g16.sh -h` 查看用法）。

[← 返回第17章](../../tutorials/17-hpc-jobs.md)
