# 06-jobflow-atomate2：jobflow 入门、atomate2 VASP 工作流与 jobflow-remote 集群提交

## 目的

本目录按由浅入深的顺序准备了三层示例：

1. **jobflow**：用 ASE 的 EMT 势代替 DFT，搭一个“7 个单点 → 1 个 EOS 拟合”的扇入工作流，几秒钟就能在笔记本上跑完，用来理解 `@job`、`Flow`、`OutputReference` 和 `run_locally`；
2. **atomate2**：用 `RelaxBandStructureMaker` 对 Si 做双弛豫 → 静态 → 均匀网格与高对称线非自洽计算；用 `PhononMaker`、`VaspLobsterMaker` 构建声子和 COHP 工作流。每个 VASP 步骤都由 custodian 包装；
3. **jobflow-remote 1.x**：把 atomate2 工作流写入 MongoDB 队列，由 runner 守护进程自动上传、提交到 Slurm、下载并解析结果。

## 文件说明

| 文件 | 作用 | 依赖 |
|---|---|---|
| `jobflow_emt_eos.py` | EMT 版 Cu 状态方程工作流，`--folders` 时每个 job 建一个目录 | jobflow、ase |
| `atomate2_bandstructure_local.py` | Si 弛豫 + 能带（`RelaxBandStructureMaker`），`--dry-run` 只打印作业图 | atomate2（Python ≥ 3.12）、VASP、POTCAR |
| `atomate2_more_flows.py` | `PhononMaker`（Si）和 `VaspLobsterMaker`（GaAs），默认只列出作业，`--run` 才真正运行 | `atomate2[phonons,lobster]` |
| `atomate2.yaml.example` | atomate2 设置模板，复制为 `~/.atomate2.yaml`，最关键的是 `VASP_CMD` | — |
| `jobflow.yaml.example` | jobflow 设置模板，复制为 `~/.jobflow.yaml`，把结果存到 MongoDB（`JOB_STORE`） | MongoDB |
| `jfremote_project.yaml` | jobflow-remote 项目文件：一个远程 Slurm worker、一个本地 shell worker、队列库、结果库、执行配置 `vasp6` | jobflow-remote 1.x、MongoDB |
| `submit_jfremote.py` | 把 Si 能带（或 `--phonon` 声子）工作流提交到 jobflow-remote 队列 | jobflow-remote、atomate2 |

## 运行步骤

```bash
# 0) 环境（atomate2 0.1.x 要求 Python >= 3.12）
conda create -n wf python=3.12 -y && conda activate wf
pip install jobflow ase                     # step 1 only
pip install "atomate2[phonons,lobster]"     # steps 2-3
pip install jobflow-remote                  # step 3

# 1) jobflow 入门（任何电脑都能运行）
python jobflow_emt_eos.py
python jobflow_emt_eos.py --folders         # creates job_<UTC timestamp>-<5 random digits> folders

# 2) atomate2 + VASP（必须在已分配计算资源的作业中运行）
cp atomate2.yaml.example ~/.atomate2.yaml   # edit VASP_CMD, VASP_INCAR_UPDATES
pmg config --add PMG_VASP_PSP_DIR /path/to/pymatgen_potcars
python atomate2_bandstructure_local.py --dry-run   # see the job graph first
sbatch my_atomate2_job.sh                          # a job script that runs the line below
#   python atomate2_bandstructure_local.py
python atomate2_more_flows.py phonon               # list the top-level phonon jobs
python atomate2_more_flows.py lobster              # list the LOBSTER jobs

# 3) jobflow-remote（runner 所在机器需能连上 MongoDB 和集群）
mkdir -p ~/.jfremote && cp jfremote_project.yaml ~/.jfremote/dft_tutorial.yaml
jf project check --errors
jf admin reset                              # first time only: wipes the queue database
jf runner start
python submit_jfremote.py POSCAR
jf job list
jf flow list
```

作业脚本可以用本板块的 [mkjob.py](../../scripts/mkjob.py) 生成（`python ../../scripts/mkjob.py -h` 查看用法），也可以参考 [第2章](../../tutorials/02-job-scheduler.md) 的 Slurm 模板。

## 预期结果 / 如何检查

- `jobflow_emt_eos.py` 最后打印一行 `a0 = ... A   B = ... GPa   (EMT, not DFT)`。EMT 是经验势，结果只用来验证流程能跑通，不能当作物性数据。
- `atomate2_bandstructure_local.py --dry-run` 列出 5 个作业，名称类似 `relax 1`、`relax 2`、`static`、`non-scf uniform`、`non-scf line`（具体名称以所装版本为准）。
- 真正运行后，每个作业有一个 `job_*` 目录，其中包含 gzip 压缩的 VASP 输出和 `custodian.json.gz`（atomate2 默认 `VASP_ZIP_FILES: atomate`，会把 custodian 日志一起压缩）；脚本打印每一步的能量和带隙。检查方法：

```bash
ls -d job_*                                   # one folder per VASP step
zgrep -c "reached required accuracy" job_*/OUTCAR*.gz
python ../06-custodian-vasp-si/show_custodian_log.py job_XXXX/custodian.json.gz   # handles .gz
```

- jobflow-remote：`jf job list` 中的状态依次为 `READY → CHECKED_OUT → UPLOADED → SUBMITTED → RUNNING → RUN_FINISHED → DOWNLOADED → COMPLETED`。出错的作业显示为 `FAILED` 或 `REMOTE_ERROR`，用 `jf job info <id>` 查看原因。

## 注意事项

- 声子工作流的位移计算由 `run_phonon_displacements` 在运行时动态生成，所以 `atomate2_more_flows.py phonon` 只列出顶层作业，看不到位移计算的个数。
- `VASP_CMD` 是完整的启动命令（如 `srun vasp_std`）。atomate2 先展开其中的环境变量，再按 shell 规则拆分后交给 custodian。
- 不配置 `~/.jobflow.yaml` 时，`run_locally` 使用内存数据库（MemoryStore），脚本结束后结果只保留在磁盘上的 `job_*` 目录里。
- jobflow-remote 会忽略 `~/.jobflow.yaml`，所有数据库设置都写在项目文件里。修改项目文件后要 `jf runner stop` 再 `jf runner start`。
- jobflow-remote 用 `jf job rerun` 在**同一目录**重跑作业。重跑前先清理或移走旧的 VASP 输出，避免旧 CONTCAR/WAVECAR 被误读。
- `jf admin reset` 会清空队列数据库，只在第一次使用时运行。
- 项目文件中写有数据库密码，请设置 `chmod 600 ~/.jfremote/*.yaml`，不要提交到公开仓库。

返回教程：[第6章 工作流自动化与高通量](../../tutorials/06-workflow-automation.md)
