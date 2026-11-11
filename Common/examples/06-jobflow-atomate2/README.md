# 示例 06-jobflow-atomate2（Common example）

本目录是 Common 板块第 6 章的配套示例。各文件的用途与参数说明见文件头部注释。

## 文件列表

| 文件 | 类型 |
|---|---|
| [atomate2.yaml.example](atomate2.yaml.example) | .example |
| [atomate2_bandstructure_local.py](atomate2_bandstructure_local.py) | .py |
| [atomate2_more_flows.py](atomate2_more_flows.py) | .py |
| [jfremote_project.yaml](jfremote_project.yaml) | .yaml |
| [jobflow.yaml.example](jobflow.yaml.example) | .example |
| [jobflow_emt_eos.py](jobflow_emt_eos.py) | .py |
| [submit_jfremote.py](submit_jfremote.py) | .py |

## 使用方法

1. 将整个目录复制到计算集群的工作目录；
2. 阅读各文件头部注释，按需修改路径、核数与计算参数；
3. 脚本类文件可用 `python <脚本> -h` 或 `bash -n <脚本>` 先行检查。

返回：[Common 目录](../../README.md)
