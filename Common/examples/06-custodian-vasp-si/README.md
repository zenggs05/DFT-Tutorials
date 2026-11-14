# 06-custodian-vasp-si：用 custodian 自动纠错运行 VASP（Si 弛豫 → 弛豫 → 静态）

## 目的

用 [custodian](https://materialsproject.github.io/custodian/) 包装一个最小的 VASP 计算序列：体相 Si 两次全弛豫（`ISIF = 3`），然后在弛豫后的结构上做一次静态计算。VASP 运行期间，custodian 每隔一段时间检查一次 `vasp.out`、`OSZICAR`、`vasprun.xml`。一旦发现 ZBRENT、EDDDAV、EDDRMM、BRMIX、电子步不收敛等已知错误，就把现场打包备份成 `error.N.tar.gz`，按内置规则修改 INCAR/POSCAR/KPOINTS 后重启，所有修正记录写入 `custodian.json`。

本例同时给出 Python 和 YAML 两种写法，并附一个只用标准库的日志查看脚本。

## 文件说明

| 文件 | 作用 |
|---|---|
| `INCAR` | Si 全弛豫参数：`ENCUT = 400`（约 1.6 × ENMAX，避免 Pulay 应力）、`ALGO = Fast`、`ISIF = 3`、`EDIFFG = -0.01` |
| `KPOINTS` | Γ 中心 8×8×8 网格 |
| `POSCAR` | 金刚石结构 Si 原胞，a = 5.43 Å |
| `custodian_vasp_run.py` | Python 写法：3 个 `VaspJob`（`.relax1`/`.relax2`/`.static`）+ 8 个 handler + 2 个 validator；`max_errors = 10` |
| `custodian_spec.yaml` | YAML 写法，用 `cstdn run custodian_spec.yaml` 运行（弛豫 → 静态两步） |
| `job_custodian_slurm.sh` | Slurm 作业脚本：加载模块和 conda 环境、生成 POTCAR、去除 CRLF、设置 `VASP_CMD` 和 `CUSTODIAN_WALLTIME` 后运行 custodian |
| `show_custodian_log.py` | 只用标准库，把 `custodian.json`（或 atomate2 目录中的 `custodian.json.gz`）整理成可读列表（哪个 job、哪个 handler、什么错误、改了什么） |

> [!IMPORTANT]
> 仓库中不包含 POTCAR（VASP 许可证限制）。作业脚本会在 POTCAR 不存在时从 `$VASP_PP_PATH/potpaw_PBE/Si/POTCAR` 生成，请确认该路径指向你有授权的赝势库。

## 运行步骤

```bash
# 0) 环境：Python >= 3.11 的 conda 环境中安装 custodian
conda create -n wf python=3.12 -y && conda activate wf
pip install custodian            # 2026.10.21 at the time of writing

# 1) 复制到集群工作目录
rsync -av 06-custodian-vasp-si/ hpc:/work/$USER/si_custodian/
ssh hpc && cd /work/$USER/si_custodian

# 2) 按集群修改 job_custodian_slurm.sh 中的 module、conda 环境名、分区和 --time
#    （CUSTODIAN_WALLTIME 必须和 --time 一致，单位为秒）
bash -n job_custodian_slurm.sh   # syntax check only
sbatch job_custodian_slurm.sh

# 3) 结束后查看纠错记录（本地电脑上也能运行，不需要装 custodian）
python show_custodian_log.py custodian.json
```

YAML 写法：在作业脚本中把 `python custodian_vasp_run.py` 换成 `cstdn run custodian_spec.yaml`。YAML 中 `$vasp_cmd` 里的 `$SLURM_NTASKS` 会在运行时按环境变量展开。

## 预期结果 / 如何检查

正常结束后，目录中应当有：

```text
INCAR.orig  KPOINTS.orig  POSCAR.orig             <- backups made by the first VaspJob
INCAR.relax1  OUTCAR.relax1  vasprun.xml.relax1  CONTCAR.relax1 ...   (copies)
INCAR.relax2  OUTCAR.relax2  vasprun.xml.relax2  CONTCAR.relax2 ...   (copies)
INCAR.static  OUTCAR.static  vasprun.xml.static  CONTCAR.static ...   (final=True: moved)
custodian.json                                     <- list of jobs and corrections
error.1.tar.gz ...                                 <- only if a handler fired
```

非最终步（`final=False`）的输出被**复制**成带后缀的文件；最终步（`final=True`）的输出被**移动**成带后缀的文件，所以结束后目录里没有不带后缀的 OUTCAR。

检查要点：

```bash
grep -c "reached required accuracy" OUTCAR.relax1 OUTCAR.relax2   # 1 for each relaxation
grep "free  energy   TOTEN" OUTCAR.static | tail -1                # static total energy
python show_custodian_log.py                                       # "no corrections needed" for a clean run
ls error.*.tar.gz 2>/dev/null                                      # empty for a clean run
```

`show_custodian_log.py` 的输出格式如下（示意输出，数值仅供参考）：

```text
=== job 1: VaspJob .relax1
    no corrections needed
=== job 2: VaspJob .relax2
  [1] handler = VaspErrorHandler
      errors  = zbrent
      action  = {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}}
      action  = {"dict": "INCAR", "action": {"_set": {"EDIFF": 1e-07}}}
      action  = {"dict": "INCAR", "action": {"_set": {"NELMIN": 8}}}
=== job 3: VaspJob .static
    no corrections needed
```

Si 原胞很小，正常情况下不会触发任何 handler。要观察纠错过程，可以人为制造问题，例如把 `NELM` 改成 6、`ALGO` 改成 `VeryFast`，此时 `UnconvergedErrorHandler` 会按 VeryFast → Fast → Normal → All 的顺序逐级修改 `ALGO`，最后改混合参数。handler 不会修改 NELM，所以这个作业最终会以 `MaxCorrectionsPerJobError` 或 `NonRecoverableError` 结束，这正说明自动纠错救不了输入本身的错误。

## 注意事项

- custodian 在**计算节点**上运行，由它启动 VASP。不要在登录节点运行 `custodian_vasp_run.py`。
- `vasp_cmd` 必须是参数列表（脚本中用 `shlex.split` 转换），不能是一整个字符串。
- `VaspErrorHandler` 只在 `vasp.out`（VaspJob 默认把标准输出重定向到这里）中找报错，不要自己再重定向 VASP 的输出。
- custodian 会改写 INCAR，原始文件保存在 `*.orig`。重新提交前先确认 INCAR 是否已被修改过。
- `max_errors` 是整个序列允许的修正次数上限，超过后抛出 `MaxCorrectionsError`；`max_errors_per_job` 是单个 VaspJob 的上限。
- Slurm 不会把 `--time` 传给作业环境，因此脚本中显式设置了 `CUSTODIAN_WALLTIME`，供 `WalltimeHandler` 在到时前写入 STOPCAR。

返回教程：[第6章 工作流自动化与高通量](../../tutorials/06-workflow-automation.md)
