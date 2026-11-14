# 06-quacc-dpdispatcher：quacc 配方与 dpdispatcher 轻量批量提交

## 目的

- **quacc**：配方（recipe）就是带 `@job`/`@flow` 装饰器的普通 Python 函数，内部调用 ASE 计算器。不配置工作流引擎时调用即执行，不需要数据库和守护进程。这里用 EMT 势演示 Cu 体相弛豫和 Cu(111) 表面能，把 import 换成 `quacc.recipes.vasp.core` 或 `quacc.recipes.gaussian.core` 就能改用 VASP 或 Gaussian。
- **dpdispatcher**：不建工作流，只负责“把一批目录变成作业、提交、等待、取回文件”。这里用它批量提交 VASP 的 ENCUT 扫描（300–600 eV，7 个点），结束后汇总 E0。

## 文件说明

| 文件 | 作用 |
|---|---|
| `quacc_emt_demo.py` | quacc EMT 配方：`relax_job(relax_cell=True)` 弛豫体相，`relax_job` + `static_job` 计算 Cu(111) 表面能（仅作示意） |
| `quacc.yaml.example` | quacc 设置模板，复制为 `~/.quacc.yaml`；VASP 的启动命令由 `VASP_PARALLEL_CMD` 和 `VASP_CMD` 拼接而成 |
| `dpdisp_encut_scan.py` | 由模板目录生成 `work/encut_*`，调用 dpdispatcher 提交并汇总能量；支持 `--dry-run` 和 `--local` |
| `dpdispatcher_machine.json` | 机器：Slurm + SSHContext（从自己电脑经 SSH 提交），`remote_root` 为集群上的工作目录 |
| `dpdispatcher_resources.json` | 资源：1 节点 32 核、分区 `cpu`、`group_size`、模块、环境变量、`ulimit` |

## 运行步骤

```bash
# quacc（quacc 1.8.x 要求 Python >= 3.11）
pip install quacc
python quacc_emt_demo.py                       # writes one quacc-* results folder per job

# 用 quacc 跑 VASP 前
cp quacc.yaml.example ~/.quacc.yaml            # edit VASP_PARALLEL_CMD, VASP_PP_PATH, RESULTS_DIR

# dpdispatcher
pip install dpdispatcher
# 编辑 dpdispatcher_machine.json 中的 hostname/username/key_filename/remote_root
# key_filename 必须写绝对路径：dpdispatcher 1.1.0 不展开 "~"
# 以及 dpdispatcher_resources.json 中的分区、核数、模块
python dpdisp_encut_scan.py --template ../06-custodian-vasp-si --dry-run   # upload + write scripts only
python dpdisp_encut_scan.py --template ../06-custodian-vasp-si             # submit and wait (use tmux)
```

`--template` 目录中必须有 POTCAR（仓库中不提供）。`run_submission()` 会一直等到所有作业结束，请在 tmux/screen 中运行；中途断开后重新运行同一命令，dpdispatcher 会接着跟踪已提交的作业，而不是重复提交。

## 预期结果 / 如何检查

- `quacc_emt_demo.py` 打印体相晶格常数、每原子能量和 Cu(111) 表面能（EMT 结果只用来演示流程），并给出每个 job 的结果目录 `dir_name`。返回值是普通字典，`atoms`、`results`、`dir_name` 等键可以直接传给下一个配方。
- `dpdisp_encut_scan.py` 结束时打印 ENCUT–E0 表，格式如下（示意输出，数值仅供参考）：

```text
 ENCUT          E0 (eV)
   300     -10.8xxxxxxx
   350     -10.8xxxxxxx
   ...
   600     -10.8xxxxxxx
```

检查要点：相邻 ENCUT 的每原子能量差小于 1 meV/atom 即可认为收敛，判据见 [VASP 第3章](../../../VASP/tutorials/03-convergence.md)。每个 `work/encut_*` 中应有取回的 `OUTCAR`、`OSZICAR`、`vasp.out`；失败的任务可以用 `grep -L "General timing" work/*/OUTCAR` 找出来。

## 注意事项

- dpdispatcher **不做纠错**。它只负责提交和取回，作业失败后要自己检查并重提，或者把 `command` 换成调用 custodian 的脚本（见 [06-custodian-vasp-si](../06-custodian-vasp-si/README.md)）。
- `group_size` > 1 时，多个任务会写进同一个 Slurm 脚本里依次运行，适合大量小任务，可以减少排队作业数。
- `backward_files` 只列需要的小文件，不要取回 WAVECAR、CHGCAR。
- quacc 的 `WORKFLOW_ENGINE` 只能通过配置文件或环境变量 `QUACC_WORKFLOW_ENGINE` 设置，不能在脚本运行中途切换。
- quacc 调用 VASP 时默认经由 custodian（`VASP_USE_CUSTODIAN: true`），所以同样会生成 `custodian.json`。

返回教程：[第6章 工作流自动化与高通量](../../tutorials/06-workflow-automation.md)
