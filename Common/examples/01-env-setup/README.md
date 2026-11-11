# 示例 01：DFT 计算环境配置模板（conda 环境、bashrc、pymatgen/ASE/VASPKIT/SSH 配置）

## 目的

- 用一份 `environment.yml` 建好 Python 3.12 的 `dft` 环境（ASE、pymatgen、phonopy、sumo、pyprocar、py4vasp、cclib、GoodVibes、LobsterPy、mp-api、dpdispatcher、RDKit）；
- 提供 `~/.bashrc`、`~/.config/.pmgrc.yaml`、`~/.config/ase/config.ini`、`~/.vaspkit`、`~/.ssh/config`、`~/.condarc` 的可直接修改的模板；
- 用 `print_versions.py` 和仓库脚本 `env_check.sh` 自检。

对应教程：[第1章 计算环境搭建](../../tutorials/01-linux-hpc.md)。

## 文件说明

| 文件 | 复制到 | 说明 |
|---|---|---|
| `environment.yml` | （直接使用） | conda-forge + pip 的环境定义，Python 3.12 |
| `bashrc_dft_snippet.sh` | 追加到 `~/.bashrc` | `VASP_PP_PATH`、`ASE_VASP_COMMAND`、MP key、Gaussian/CP2K/OpenMX/Multiwfn 变量、别名 |
| `pmgrc.yaml.example` | `~/.config/.pmgrc.yaml` | pymatgen 的 POTCAR 目录、默认泛函、API key |
| `ase_config.ini.example` | `~/.config/ase/config.ini` | ASE 新式配置文件（替代环境变量） |
| `vaspkit_rc.example` | `~/.vaspkit` | VASPKIT 1.5.x 配置，按官方模板键名整理 |
| `ssh_config.example` | `~/.ssh/config` | 集群别名、心跳、连接复用、跳板机 |
| `condarc_mirror.example` | `~/.condarc` | 只用 conda-forge，可选清华镜像 |
| `print_versions.py` | （直接运行） | 打印 Python 与各包版本，缺包不会报错 |

完整的环境检查脚本在仓库的 Common/scripts/env_check.sh（检查包版本、环境变量、`~/.vaspkit`、`ulimit -s`、`OMP_NUM_THREADS`、`OMP_STACKSIZE`，并给出修复建议）。

## 运行步骤

```bash
cd Common/examples/01-env-setup

# 1) Python 环境
mamba env create -f environment.yml        # 或 conda env create -f environment.yml
conda activate dft
python print_versions.py | tee versions.txt
pip check

# 2) Shell 环境变量（先编辑其中的路径）
cp bashrc_dft_snippet.sh ~/dft_env.sh && vi ~/dft_env.sh
echo 'source ~/dft_env.sh' >> ~/.bashrc && source ~/.bashrc

# 3) pymatgen（推荐用命令生成配置，而不是手抄模板）
pmg config -p ~/vasp_pp ~/pmg_potcars
pmg config --add PMG_VASP_PSP_DIR ~/pmg_potcars PMG_DEFAULT_FUNCTIONAL PBE_54

# 4) ASE 配置文件（可选，替代环境变量）
mkdir -p ~/.config/ase && cp ase_config.ini.example ~/.config/ase/config.ini && vi ~/.config/ase/config.ini
ase info --calculators | head -40

# 5) VASPKIT
cp vaspkit_rc.example ~/.vaspkit && vi ~/.vaspkit      # 改 PBE_PATH、PYTHON_BIN

# 6) SSH（在你自己的电脑上）
cp ssh_config.example ~/.ssh/config && chmod 600 ~/.ssh/config && vi ~/.ssh/config

# 7) 全面自检
bash ../../scripts/env_check.sh
```

## 预期结果 / 如何检查

- `python print_versions.py` 列出各包版本（2026-11 时的参考：ase 3.29.0、pymatgen 2026.10.24、phonopy 4.10.0、spglib 2.8.0、mp-api 0.47.0；你的版本可能更新），且没有 “too old” 警告。
- `echo $VASP_PP_PATH` 指向一个**包含** `potpaw_PBE/` 的目录：`ls $VASP_PP_PATH/potpaw_PBE/Si/POTCAR` 存在。
- `python -c "from pymatgen.core import SETTINGS; print(SETTINGS.get('PMG_VASP_PSP_DIR'))"` 输出你的 pymatgen POTCAR 目录。
- `ase info --calculators` 的 Configuration 部分显示 `Loaded: /home/.../.config/ase/config.ini`。
- `ssh hpc hostname` 不需要密码即可返回集群登录节点名。

## 注意事项

- 所有路径（`/home/your_name/...`）都是示例，必须改成你自己的。
- POTCAR、API key、私钥都**不能**提交到 git 仓库或公开分享；含 key 的文件请 `chmod 600`。
- `pmg config --add` 一次可以写多个“键 值”对；它会把旧文件备份为 `.pmgrc.yaml.bak`。
- `environment.yml` 用的是版本下限（`>=`），便于长期使用；需要精确复现时用 `conda env export > environment-full.yml` 另存一份锁定版本。
- `vaspkit_rc.example` 中的键名来自 VASPKIT 官方安装文档；如果你下载的版本自带的 `how_to_set_environment_variable` 有新增条目，以自带文件为准。
- 本目录中除 README 外的文件均为 ASCII、LF 换行；在 Windows 下编辑后请确认没有变成 CRLF。

[← 返回第1章](../../tutorials/01-linux-hpc.md)
