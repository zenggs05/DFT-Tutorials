# 示例 18-job-scripts（OpenMX example）

本目录是 OpenMX 板块第 18 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [check_run.sh](check_run.sh) | .sh |
| [job_pbs.sh](job_pbs.sh) | .sh |
| [job_slurm.sh](job_slurm.sh) | .sh |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[OpenMX 目录](../../README.md)
