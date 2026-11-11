# 示例 02-job-templates（Common example）

本目录是 Common 板块第 2 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [lsf_mpi.sh](lsf_mpi.sh) | .sh |
| [pbs_mpi.sh](pbs_mpi.sh) | .sh |
| [slurm_array.sh](slurm_array.sh) | .sh |
| [slurm_hybrid_mpi_omp.sh](slurm_hybrid_mpi_omp.sh) | .sh |
| [slurm_mpi.sh](slurm_mpi.sh) | .sh |
| [slurm_shared_memory.sh](slurm_shared_memory.sh) | .sh |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[Common 目录](../../README.md)
