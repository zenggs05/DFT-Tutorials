# 示例 02-array-chain-si（Common example）

本目录是 Common 板块第 2 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [chain_relax_static_band.sh](chain_relax_static_band.sh) | .sh |
| [collect_encut.sh](collect_encut.sh) | .sh |
| [dpdispatcher_encut.py](dpdispatcher_encut.py) | .py |
| [INCAR.band](INCAR.band) | .band |
| [INCAR.encut.template](INCAR.encut.template) | .template |
| [INCAR.relax](INCAR.relax) | .relax |
| [INCAR.static](INCAR.static) | .static |
| [job_step.sh](job_step.sh) | .sh |
| [KPOINTS](KPOINTS) | 输入文件 |
| [KPOINTS.band](KPOINTS.band) | .band |
| [POSCAR](POSCAR) | 输入文件 |
| [slurm_array_encut.sh](slurm_array_encut.sh) | .sh |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[Common 目录](../../README.md)
