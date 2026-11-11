# 示例 06-custodian-vasp-si（Common example）

本目录是 Common 板块第 6 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [custodian_spec.yaml](custodian_spec.yaml) | .yaml |
| [custodian_vasp_run.py](custodian_vasp_run.py) | .py |
| [INCAR](INCAR) | 输入文件 |
| [job_custodian_slurm.sh](job_custodian_slurm.sh) | .sh |
| [KPOINTS](KPOINTS) | 输入文件 |
| [POSCAR](POSCAR) | 输入文件 |
| [show_custodian_log.py](show_custodian_log.py) | .py |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[Common 目录](../../README.md)
