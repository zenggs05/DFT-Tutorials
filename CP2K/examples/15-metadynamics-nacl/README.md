# 示例 15-metadynamics-nacl（CP2K example）

本目录是 CP2K 板块第 15 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [bluemoon_integrate.py](bluemoon_integrate.py) | .py |
| [constrained_md_distance.inp](constrained_md_distance.inp) | .inp |
| [fes_from_hills.py](fes_from_hills.py) | .py |
| [job_slurm.sh](job_slurm.sh) | .sh |
| [last_frame.sh](last_frame.sh) | .sh |
| [nacl_dft_force_eval.inc](nacl_dft_force_eval.inc) | .inc |
| [nacl_dissociation_wtmetad.inp](nacl_dissociation_wtmetad.inp) | .inp |
| [nacl_equil.inp](nacl_equil.inp) | .inp |
| [nacl_water.xyz](nacl_water.xyz) | .xyz |
| [nacl_xtb_force_eval.inc](nacl_xtb_force_eval.inc) | .inc |
| [plumed.dat](plumed.dat) | .dat |
| [plumed_metad.inp](plumed_metad.inp) | .inp |
| [run_bluemoon.sh](run_bluemoon.sh) | .sh |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[CP2K 目录](../../README.md)
