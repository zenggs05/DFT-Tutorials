# 06-aiida-cp2k-vasp：AiiDA 驱动 CP2K、VASP、Gaussian（带完整溯源）

## 目的

用 [AiiDA](https://aiida.readthedocs.io/) 运行三种程序的最小示例，并把每次计算的输入、输出、代码版本和调用关系存进数据库（provenance graph）：

- **aiida-cp2k**：H2O 分子 PBE 单点能，CP2K 输入用嵌套 Python 字典写出；
- **aiida-vasp 5.x**：用 `VaspRelaxWorkChain.get_builder_from_protocol` 和 `balanced` 协议弛豫体相 Si（协议默认是 PBEsol，脚本改成了 PBE），WorkChain 会反复弛豫直到体积变化足够小，最后做一次静态计算；
- **aiida-gaussian**：用 `GaussianBaseWorkChain` 优化 H2O，SCF 不收敛时自动依次改用 `SCF=YQC`、`SCF=XQC` 重算。

同时给出 Computer（集群）和 Code（可执行程序）的 YAML 配置模板。

## 文件说明

| 文件 | 作用 |
|---|---|
| `computer_slurm.yaml` | `verdi computer setup --config`：SSH + Slurm 集群 `hpc`，工作目录、`mpirun_command`、预置命令 |
| `computer_ssh.yaml` | `verdi computer configure core.ssh hpc --config`：用户名、私钥、超时等 SSH 参数 |
| `cp2k_code.yaml` | 在 `localhost` 上注册 `cp2k.psmp`（插件 `cp2k`） |
| `vasp_code.yaml` | 在 `localhost` 上注册 `vasp_std`（插件 `vasp.vasp`） |
| `gaussian_code.yaml` | 在 `localhost` 上注册 `g16`（插件 `gaussian`，`with_mpi: false`） |
| `aiida_cp2k_sp.py` | aiida-cp2k 单点能，打印 `output_parameters` 中的能量（Hartree 与 eV） |
| `aiida_vasp_relax.py` | aiida-vasp 5.x 弛豫，`overrides` 覆盖 ENCUT 等 INCAR 参数，`relax_settings` 设置力收敛标准 |
| `aiida_gaussian_base.py` | aiida-gaussian 的 `GaussianBaseWorkChain`，B3LYP/6-31G(d) 优化 |

## 运行步骤

```bash
# 0) 安装（aiida-core 2.9.x 支持 Python >= 3.10）
conda create -n aiida python=3.12 -y && conda activate aiida
pip install aiida-core aiida-cp2k aiida-vasp aiida-gaussian

# 1) 建 profile：SQLite 存储，自动配置 localhost 计算机
verdi presto
verdi status                                # all lines should be green / ok

# 2) 注册程序（先修改 YAML 中的可执行文件路径）
verdi code create core.code.installed --config cp2k_code.yaml
verdi code create core.code.installed --config vasp_code.yaml
verdi code create core.code.installed --config gaussian_code.yaml
verdi code list

# 3) VASP 赝势只需上传一次，族名要与协议中的 potential_family 一致（默认 PBE.54）
aiida-vasp potcar uploadfamily --path /path/to/potpaw_PBE.54 --name PBE.54 --description "PBE.54"

# 4) 运行
verdi run aiida_cp2k_sp.py
verdi run aiida_vasp_relax.py
verdi run aiida_gaussian_base.py

# 5) 远程集群（可选）
verdi computer setup --config computer_slurm.yaml
verdi computer configure core.ssh hpc --config computer_ssh.yaml
verdi computer test hpc
#    然后把 *_code.yaml 中的 computer 改为 hpc、label 改名后重新创建 code，用 --code xxx@hpc 运行
verdi daemon start 2
python aiida_vasp_relax.py --code vasp_std@hpc --submit
```

## 预期结果 / 如何检查

```bash
verdi process list -a -p 1              # processes of the last day, with exit status
verdi process show <PK>                 # inputs, outputs, called sub-processes
verdi process report <PK>               # WorkChain log: which handler restarted what
verdi calcjob outputcat <PK>            # main output file of a CalcJob (aiida.out / OUTCAR ...)
verdi calcjob gotocomputer <PK>         # open a shell in the remote working directory
verdi node graph generate <PK>          # provenance graph as a PDF
```

成功的标志：`verdi process list` 中状态为 `Finished [0]`。`Finished [非零]` 表示进程结束但计算失败，用 `verdi process report` 查看原因。

脚本最后的打印格式如下（示意输出，数值仅供参考）：

```text
CalcJob pk=123  exit_status=0
E(PBE) = -17.2xxxxxxx Ha = -46.9xxxx eV
provenance graph: verdi node graph generate 123
```

## 注意事项

- `verdi presto` 建的是 SQLite profile，适合学习和单人小规模使用。课题组生产环境建议用 `verdi presto --use-postgres` 或 `verdi profile setup core.psql_dos` 配合 PostgreSQL。
- 没有可用的消息中间件（broker）时只能用 `verdi run` / `run_get_node` 阻塞运行，不能 `submit` 给 daemon。aiida-core 2.9 的 `verdi presto` 检测不到 RabbitMQ 时会改用 ZeroMQ broker；更早的版本会建立无 broker 的 profile，需要另装 RabbitMQ。
- AiiDA 会把 POTCAR 内容存进数据库。导出（`verdi archive create`）和公开分享数据库前，必须确认不包含受许可证保护的文件。
- aiida-cp2k 自动生成 `aiida.inp`，并根据 `StructureData` 写出坐标文件和晶胞。不要在参数字典里再写 `CELL/A`、`B`、`C`，否则会冲突。
- aiida-gaussian 用 cclib（当前稳定版 1.8.1）解析输出，`scfenergies` 的单位是 eV。
- aiida-vasp 5.1.0 的 `balanced` 协议默认 `GGA = PS`（PBEsol）、`ISPIN = 2`、`ENCUT = 550`。`aiida_vasp_relax.py` 在 `overrides` 中显式改成 `gga = pe`（PBE）、`ispin = 1`、`encut = 400`；用协议时一定要用 `verdi process show` 或 `verdi calcjob inputcat <PK> INCAR` 确认最终 INCAR。

返回教程：[第6章 工作流自动化与高通量](../../tutorials/06-workflow-automation.md)
