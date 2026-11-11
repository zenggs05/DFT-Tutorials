# 示例 06-aiida-cp2k-vasp（Common example）

本目录是 Common 板块第 6 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [aiida_cp2k_sp.py](aiida_cp2k_sp.py) | .py |
| [aiida_vasp_relax.py](aiida_vasp_relax.py) | .py |
| [computer_slurm.yaml](computer_slurm.yaml) | .yaml |
| [computer_ssh.yaml](computer_ssh.yaml) | .yaml |
| [cp2k_code.yaml](cp2k_code.yaml) | .yaml |
| [vasp_code.yaml](vasp_code.yaml) | .yaml |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[Common 目录](../../README.md)
