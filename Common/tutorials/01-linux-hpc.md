# 第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输（Setting up a DFT Computing Environment: Linux/WSL, conda, SSH & File Transfer）

本章从零搭建第一性原理计算的日常工作环境：在 Windows 上用 WSL2（或直接用 Linux/macOS）获得 Linux 命令行，用 SSH 密钥登录 HPC 集群，用 scp/rsync 传输文件、用 tmux/screen 保持会话，用 conda/mamba 建一个 Python ≥ 3.11 的 `dft` 环境，安装 ASE、pymatgen、phonopy、sumo、pyprocar、py4vasp、cclib、GoodVibes、LobsterPy 等工具，并正确设置 `VASP_PP_PATH`、`PMG_VASP_PSP_DIR`、`MP_API_KEY` 等环境变量。最后用 `environment.yml` 锁定版本、用自检脚本确认一切就绪（Linux / WSL2 / conda / SSH / rsync setup for VASP, Gaussian, CP2K and OpenMX users）。

> [!NOTE]
> 本章讲的是“周边环境”：Python 前后处理工具、登录与传输方法、环境变量。VASP、Gaussian、CP2K、OpenMX 本体的编译与安装分别见 [VASP 第1章](../../VASP/tutorials/01-setup-compile.md)、[Gaussian 第1章](../../Gaussian/tutorials/01-setup-first-job.md)、[CP2K 第1章](../../CP2K/tutorials/01-install.md) 和 [OpenMX 第2章](../../OpenMX/tutorials/02-install.md)；作业脚本与并行设置见 第2章。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 获得 Linux 环境 | 在 Windows 上安装 WSL2 + Ubuntu，或在 macOS/Linux 上配好终端 |
| 登录集群 | 生成 ed25519 密钥、免密登录、写 `~/.ssh/config`、通过跳板机登录 |
| 传输文件 | 用 `scp`、`rsync` 上传输入、只下载需要的输出，断点续传 |
| 长时间操作不掉线 | 用 `tmux`/`screen` 保持会话 |
| Python 环境 | 用 Miniforge（conda/mamba）建 Python 3.12 的 `dft` 环境，装齐本教程用到的包 |
| 环境变量 | 设置 `VASP_PP_PATH`、`ASE_VASP_COMMAND`、`PMG_VASP_PSP_DIR`、`MP_API_KEY`、`CP2K_DATA_DIR` 等 |
| 软件栈一致 | 理解 `module load` 与编译器/MPI 版本的对应关系 |
| 避开换行符陷阱 | 识别并修复 Windows CRLF 换行 |
| 可复现 | 用 `environment.yml` 锁定版本，并用自检脚本确认环境 |

适用读者：刚进课题组、第一次接触集群的同学；也适合已经会用 VASP 但 Python 环境“装一个坏一个”的用户。配套示例在 [examples/01-env-setup](../examples/01-env-setup/README.md)。

---

## 2. 原理速览：你的电脑、登录节点与计算节点

一次典型的 DFT 计算会经过三类机器，弄清楚“什么事在哪里做”可以避免大部分新手问题：

```text
 你的电脑 (Windows+WSL2 / macOS / Linux)
   |  建模、可视化（VESTA、Avogadro、GaussView）、写输入、作图
   |  ssh / scp / rsync（端口 22，可能需要 VPN 或跳板机）
   v
 登录节点 login01, login02 ...           <- 所有用户共享，只做轻量操作
   |  编辑文件、module load、提交作业(sbatch/qsub/bsub)、简单后处理
   |  调度器 Slurm / PBS / LSF 分配资源
   v
 计算节点 cn001 ... cn512                <- 真正运行 vasp_std / g16 / cp2k.psmp / openmx
      共享文件系统：/home（小、有备份）、/work 或 /data（大、无备份）、
      本地盘 /tmp 或 /scratch（最快，作业结束后可能被清空）
```

几条基本规则：

1. **登录节点不跑计算**。在登录节点上直接 `mpirun -np 64 vasp_std` 会影响所有人，通常会被管理员杀掉甚至封号。几秒钟的小测试和 Python 后处理一般可以，具体以集群规定为准。
2. **计算放在大容量工作盘**。`/home` 通常有几十 GB 的配额，一个 WAVECAR 就可能把它塞满；用 `quota -s`、`df -h ~`、`du -sh *` 查看用量（有些集群用自己的命令，如 `lfs quota -h -u $USER /work`）。
3. **软件环境要和作业一致**。登录节点上 `module load` 过的东西不会自动进入作业，作业脚本里必须再写一遍（见 第2章）。
4. **大文件不下载**。WAVECAR、CHGCAR、`*.wfn`、`*.rwf` 动辄 GB 级，后处理尽量在集群上完成，只下载图和小文件。

---

## 3. 准备 Linux 命令行

### 3.1 Windows：安装 WSL2 + Ubuntu

WSL2（Windows Subsystem for Linux 2）在 Windows 10（2004 及以上）和 Windows 11 中运行真正的 Linux 内核，是 Windows 用户做 DFT 前后处理的首选。以管理员身份打开 PowerShell：

```powershell
wsl --install                      # 启用 WSL 并安装默认发行版（Ubuntu），需要重启
wsl --list --online                # 查看可安装的发行版
wsl --install -d Ubuntu-24.04      # 指定版本安装
wsl -l -v                          # 确认 VERSION 一列为 2
wsl --update                       # 更新 WSL 内核
wsl --shutdown                     # 修改 .wslconfig 后需要重启 WSL
```

首次启动会要求设置 Linux 用户名和密码（与 Windows 账户无关）。进入 Ubuntu 后先更新并安装常用工具：

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential gfortran git wget curl vim \
                    rsync tmux htop dos2unix tree unzip \
                    openmpi-bin libopenmpi-dev     # 本地编译/测试小算例时需要
```

WSL2 默认最多使用一半物理内存。若要在本地跑小算例（例如 OpenMX、CP2K 的教学算例），在 Windows 用户目录下创建 `%UserProfile%\.wslconfig`：

```ini
[wsl2]
memory=16GB        # 分配给 WSL2 的内存上限
processors=8       # 可用的逻辑核数
swap=8GB
```

> [!WARNING]
> **不要在 `/mnt/c/...` 下做计算。** WSL2 访问 Windows 盘（`/mnt/c`、`/mnt/d`）要跨文件系统，I/O 比 Linux 自己的文件系统（`~`）慢一个数量级，而且 Windows 编辑器保存的文件常带 CRLF 换行（见第 10 节）。计算放在 `~/work`，需要用 Windows 软件打开时，在资源管理器地址栏输入 `\\wsl$\Ubuntu-24.04\home\<用户名>`（新版也可用 `\\wsl.localhost\...`）访问。

推荐搭配：**Windows Terminal**（多标签终端）+ **VS Code** 的 WSL 与 Remote-SSH 扩展（在 Windows 界面里直接编辑 WSL 和集群上的文件，并且可以把换行符固定为 LF）。另一个常见选择是 MobaXterm（自带 SSH、SFTP 浏览器和 X11 显示）。

### 3.2 macOS

macOS 自带的“终端”就是类 Unix 环境，`ssh`、`scp`、`rsync` 都已可用。建议：

```bash
xcode-select --install                 # 命令行编译工具
# 安装 Homebrew 后（见 https://brew.sh）：
brew install tmux rsync dos2unix wget htop
```

Apple Silicon（M 系列）上 conda-forge 提供原生 arm64 包，下文的 Miniforge 安装方式相同；VESTA 3.90.x 开发版已有原生 Apple Silicon 构建（见 第5章）。

### 3.3 原生 Linux 或虚拟机

Ubuntu 22.04/24.04、Rocky/Alma Linux 9 都可以。包管理器分别是 `apt` 和 `dnf`，其余操作与 WSL 相同。

### 3.4 DFT 用户最常用的 Linux 命令

| 任务 | 命令 | 说明 |
|---|---|---|
| 看目录 | `ls -lhtr` | 按修改时间排序，最新的在最后 |
| 建目录 | `mkdir -p proj/si/relax` | `-p` 自动创建上级目录 |
| 复制 | `cp -r relax static`，`cp relax/CONTCAR static/POSCAR` | 续算时最常用的一行 |
| 删除 | `rm -i file`，`rm -r dir` | 集群上没有回收站，删除前三思 |
| 看大文件 | `less OUTCAR`，`tail -n 30 OUTCAR`，`tail -f vasp.out` | `less` 中 `/关键字` 搜索，`G` 跳到结尾 |
| 搜索 | `grep -n "TOTEN" OUTCAR \| tail -1` | 最后一个离子步的能量 |
| 查找文件 | `find . -name OUTCAR -newer INCAR` | 找出比 INCAR 新的 OUTCAR |
| 列处理 | `awk '{print $5}'`，`sed -n '5,10p' POSCAR` | 提取列、打印行 |
| 磁盘 | `du -sh *`，`df -h .`，`quota -s` | 找出占空间的目录 |
| 打包 | `tar -czf run.tar.gz run/`，`tar -xzf run.tar.gz` | 传输前打包小文件 |
| 权限 | `chmod +x job.sh`，`chmod 600 ~/.ssh/id_ed25519` | 脚本可执行、私钥只读 |
| 链接 | `ln -s $HOME/vasp_pp/potpaw_PBE potpaw_PBE` | 软链接，节省空间 |
| 硬件 | `nproc`，`lscpu`，`free -h`，`nvidia-smi` | 核数、NUMA、内存、GPU |
| 进程 | `top`/`htop`，`ps -u $USER`，`kill <PID>` | 登录节点上检查自己有没有残留进程 |
| 路径 | `which vasp_std`，`echo $PATH`，`type python` | 当前到底用的是哪个程序 |

几个与 DFT 直接相关的组合命令：

```bash
grep "free  energy   TOTEN" OUTCAR | tail -n 1                       # VASP 最终能量
grep -c "reached required accuracy" OUTCAR                           # 弛豫是否收敛（1 表示是）
grep "SCF Done" h2o.log | tail -n 1                                  # Gaussian 最后一次 SCF 能量
grep "ENERGY| Total FORCE_EVAL" cp2k.out | tail -n 1                 # CP2K 总能
grep "Utot\." Si.out                                                 # OpenMX 总能
for d in */; do echo -n "$d "; grep TOTEN $d/OUTCAR | tail -1 | awk '{print $5}'; done   # 批量汇总
```

---

## 4. 登录集群：SSH 密钥与 `~/.ssh/config`

### 4.1 第一次登录

```bash
ssh your_name@login.hpc.example.edu          # 默认端口 22
ssh -p 2222 your_name@login.hpc.example.edu  # 非标准端口
```

第一次连接会提示核对主机指纹（`Are you sure you want to continue connecting (yes/no/[fingerprint])?`），确认与集群文档公布的一致后输入 `yes`。很多学校的集群要求先连校园 VPN，或者经过跳板机（堡垒机）、使用动态口令（OTP）；这些以集群使用手册为准。

### 4.2 生成密钥并免密登录

在**你自己的电脑**（WSL、macOS 或 Windows PowerShell）上生成 ed25519 密钥：

```bash
ssh-keygen -t ed25519 -C "your_name@lab-laptop"
# 一路回车即可：私钥 ~/.ssh/id_ed25519，公钥 ~/.ssh/id_ed25519.pub
# 建议设置口令（passphrase），再用 ssh-agent 免去重复输入
```

把**公钥**（`.pub`）追加到集群的 `~/.ssh/authorized_keys`：

```bash
# Linux / macOS / WSL
ssh-copy-id -i ~/.ssh/id_ed25519.pub your_name@login.hpc.example.edu

# Windows PowerShell 没有 ssh-copy-id，用管道代替：
type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh your_name@login.hpc.example.edu "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys"
```

> [!CAUTION]
> 私钥 `id_ed25519`（没有 `.pub` 后缀的那个）等同于密码，**绝不能**复制到集群、发给别人或放进 git 仓库。有的集群不允许密钥登录、或要求密钥 + 动态口令，请遵守集群规定。

权限不对时 SSH 会拒绝使用密钥：

```bash
chmod 700 ~/.ssh
chmod 600 ~/.ssh/id_ed25519 ~/.ssh/config ~/.ssh/authorized_keys
chmod 644 ~/.ssh/id_ed25519.pub
chmod go-w ~                     # 集群端：家目录不能对组/其他人可写，否则 sshd 会忽略 authorized_keys
```

### 4.3 用 `~/.ssh/config` 起别名

把常用参数写进 `~/.ssh/config`，以后 `ssh hpc`、`scp file hpc:~/`、`rsync ... hpc:...` 以及 VS Code Remote-SSH 都能直接用别名。完整示例见 [ssh_config.example](../examples/01-env-setup/ssh_config.example)，核心部分：

```text
Host hpc
    HostName login.hpc.example.edu
    User your_name
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
    ControlMaster auto
    ControlPath ~/.ssh/cm-%r@%h:%p
    ControlPersist 10m

Host hpc-inner
    HostName inner-login.hpc.example.edu
    User your_name
    ProxyJump your_name@gateway.hpc.example.edu
```

| 关键字 | 作用 | 说明 |
|---|---|---|
| `HostName` | 真实主机名或 IP | `Host` 后面是你起的别名 |
| `IdentityFile` | 使用哪把私钥 | 多把密钥时必须写 |
| `ServerAliveInterval 60` | 每 60 s 发心跳 | 防止空闲时被防火墙断开 |
| `ControlMaster/ControlPersist` | 复用一条连接 | 需要动态口令的集群上，第一次认证后 10 分钟内的 scp/rsync 不再要求输入；Windows 原生 OpenSSH 不支持，WSL/macOS/Linux 可用 |
| `ProxyJump` | 经跳板机登录 | OpenSSH ≥ 7.3；等价于命令行 `ssh -J gateway inner` |

### 4.4 图形界面与端口转发

```bash
ssh -X hpc            # X11 转发（Linux/WSLg/macOS+XQuartz），可以远程打开 ase gui、xmgrace
ssh -Y hpc            # “可信”X11 转发，部分程序需要
ssh -N -L 8888:localhost:8888 hpc   # 把集群上 8888 端口的 Jupyter 映射到本地浏览器
```

> [!TIP]
> X11 转发显示 VESTA、GaussView 这类重图形程序非常慢。更高效的做法是：在集群上把结构/电荷密度处理好，用 `rsync` 拉回小文件（POSCAR、CHGCAR 截取后的 `.vasp`、`.cube`）在本地打开；或者在集群的计算节点上开 Jupyter，用 nglview 交互查看（见 第5章）。

---

## 5. 文件传输：scp、rsync、sftp

### 5.1 scp：简单拷贝

```bash
scp POSCAR INCAR KPOINTS hpc:~/work/si/relax/        # 上传多个文件
scp -r si_project hpc:~/work/                         # 上传整个目录
scp hpc:~/work/si/relax/OUTCAR ./                     # 下载
scp -P 2222 file user@host:~/                         # 注意 scp 的端口参数是大写 -P
```

### 5.2 rsync：增量同步（推荐）

`rsync` 只传输变化的部分，支持断点续传、排除大文件，是日常同步计算目录的首选：

```bash
# 上传：本地 si_project/ 的“内容”同步到集群 ~/work/si_project/
rsync -avzP si_project/ hpc:~/work/si_project/

# 下载：排除波函数、电荷密度等大文件
rsync -avzP --exclude='WAVECAR' --exclude='CHG*' --exclude='vaspout.h5' \
      --exclude='*.wfn' --exclude='*.rwf' hpc:~/work/si_project/ ./si_project/

# 只下载每个子目录里的 OUTCAR、vasprun.xml、CONTCAR（保留目录结构，删除空目录）
rsync -avzPm --include='*/' --include='OUTCAR' --include='vasprun.xml' \
      --include='CONTCAR' --exclude='*' hpc:~/work/si_project/ ./si_project/

# 先演习：-n（--dry-run）只列出将要传输的文件，不真正传输
rsync -avzn --delete si_project/ hpc:~/work/si_project/
```

| 选项 | 含义 |
|---|---|
| `-a` | 归档模式：递归、保留权限和时间戳 |
| `-v` / `-z` / `-P` | 详细输出 / 传输时压缩 / 显示进度并保留部分传输的文件（断点续传） |
| `-m` | 删除同步后为空的目录（配合 include/exclude） |
| `--exclude` / `--include` | 过滤规则，**按顺序**匹配，先匹配到的生效 |
| `--delete` | 删除目标端多余的文件，危险，先配合 `-n` 演习 |
| `-n` | dry run |

> [!IMPORTANT]
> **结尾斜杠有区别**：`rsync -a dir/ host:target/` 同步的是 `dir` 里面的内容；`rsync -a dir host:target/` 会在 target 下再建一个 `dir` 子目录。不确定时先加 `-n`。

### 5.3 其他方式

- `sftp hpc`：交互式传输（`put`、`get`、`ls`、`lcd`）。
- 图形工具：WinSCP、FileZilla、MobaXterm 左侧的 SFTP 浏览器、VS Code Remote-SSH 的资源管理器（可拖拽上传下载）。
- 大量小文件：先打包 `tar -czf results.tar.gz */OUTCAR */CONTCAR`，传一个文件远比传几千个小文件快。
- 校验完整性：`sha256sum big.tar.gz` 在两端各算一次，结果一致才删除原文件。

---

## 6. 会话保持：tmux 与 screen

网络一断，普通 SSH 会话里正在运行的命令（比如编译 VASP、`rsync` 一个大目录、交互式后处理）都会被终止。`tmux` 和 `screen` 让会话在服务器上“挂起”，断线后可以重新接上。

```bash
tmux new -s dft          # 新建名为 dft 的会话
# ……执行耗时操作……
# 按 Ctrl-b 然后按 d：脱离（detach），会话继续在后台运行
tmux ls                  # 列出会话
tmux attach -t dft       # 重新接上（可简写 tmux a -t dft）
tmux kill-session -t dft # 结束会话
```

| tmux 快捷键（先按 Ctrl-b） | 作用 |
|---|---|
| `d` | 脱离会话 |
| `c` / `n` / `p` | 新窗口 / 下一个 / 上一个 |
| `%` / `"` | 左右 / 上下分屏 |
| 方向键 | 在分屏之间切换 |
| `[` | 进入滚动模式（`q` 退出） |

`screen` 的对应用法：`screen -S dft` 新建，`Ctrl-a d` 脱离，`screen -ls` 列出，`screen -r dft` 恢复。

> [!WARNING]
> 集群通常有多个登录节点（login01、login02……），`ssh hpc` 可能被轮询分配到不同节点。**tmux 会话只存在于创建它的那个节点**，接不上时先 `hostname` 看看自己在哪个节点，再 `ssh login01` 回到原节点。另外，tmux 只是保持会话，**不是**在登录节点跑长时间计算的理由——计算一律通过作业调度系统提交。

只想让一条命令在退出登录后继续运行（例如本地 WSL 里跑一个小测试）时，也可以用 `nohup cmd > log 2>&1 &`。

---

## 7. Python 环境：Miniforge + conda/mamba

### 7.1 为什么用 Miniforge 和独立环境

- **Miniforge** 是 conda-forge 社区维护的最小安装包，默认只用 conda-forge 频道，并自带 `mamba`。Anaconda 的 defaults 频道受其服务条款约束（大型机构可能需要商业许可），conda-forge 没有这一限制，而且材料/化学类软件包（phonopy、spglib、rdkit、pymatgen 等）更新也更及时。
- **不要往 base 环境里装东西**。为 DFT 工具单独建一个 `dft` 环境；需要不同版本时（例如复现老论文）再建 `dft-old`，互不干扰。
- 新版 conda（23.10 起）默认使用 libmamba 求解器，`conda` 和 `mamba` 的速度已经相近，下文两者可以互换。

### 7.2 安装 Miniforge（不需要 root）

```bash
cd ~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3     # -b 静默安装，-p 安装路径
$HOME/miniforge3/bin/conda init bash                                  # 往 ~/.bashrc 写入初始化代码
source ~/.bashrc
conda --version && mamba --version
```

> [!TIP]
> 集群家目录配额小时，把 Miniforge 装到工作盘（`-p /work/$USER/miniforge3`），或用 `conda config --add pkgs_dirs /work/$USER/conda_pkgs` 把包缓存挪走；定期 `conda clean -a` 清理缓存。国内网络慢可以参考 [condarc_mirror.example](../examples/01-env-setup/condarc_mirror.example) 配置清华 TUNA 镜像，pip 镜像用 `pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple`。

### 7.3 为什么是 Python ≥ 3.11（本教程用 3.12）

2026 年发布的几个核心包都提高了 Python 下限（数据取自 PyPI，2026-11-10）：

| 包 | 最新版本 | 要求的 Python | 说明 |
|---|---|---|---|
| ase | 3.29.0 | ≥ 3.10 | 3.27 起要求 3.10；3.28 起 NEB 默认 `improvedtangent` |
| pymatgen | 2026.10.24 | ≥ 3.11 | 2026.3.23 起核心代码拆为 `pymatgen-core`（2026.10.31），`pip install pymatgen` 会自动装上 |
| phonopy | 4.10.0 | ≥ 3.11 | 4.8.2 起要求 3.11；旧 Python 会被 pip 静默装成旧版 |
| mp-api | 0.47.0 | ≥ 3.12 | Materials Project 新 API 客户端 |
| lobsterpy | 0.6.1 | ≥ 3.11, < 3.15 | LOBSTER 成键分析 |
| sumo | 3.0.0 | ≥ 3.10 | 能带/DOS 出版级作图 |
| py4vasp | 0.11.3 | ≥ 3.10 | 读取 vaspout.h5 |
| pyprocar | 6.5.0 | ≥ 3.8 | v6 起 API 改为 `dirname=` 参数 |
| cclib | 1.8.1 | ≥ 3.7 | 2.0a1 为 alpha 版，解析 Gaussian 有已知 bug |
| goodvibes | 4.4.0 | ≥ 3.9 | 命令改为 `goodvibes`（不再用 `python -m goodvibes`） |
| ovito（pip 版） | 3.17.2 | ≥ 3.11 | 官方提示不要 pip 装进 Anaconda 环境，见第5章 |

结论：**Python 3.12 能同时满足以上所有包**；3.13 也可以，3.10 及以下会导致 pymatgen、phonopy、mp-api 无法升级到当前版本。

### 7.4 创建 `dft` 环境

推荐直接用配套的 [environment.yml](../examples/01-env-setup/environment.yml)：

```bash
cd Common/examples/01-env-setup
mamba env create -f environment.yml          # 约需 5-15 分钟，取决于网络
conda activate dft
python print_versions.py                      # 打印各包版本
```

`environment.yml` 的思路：能从 conda-forge 装的（numpy、scipy、ase、pymatgen、spglib、seekpath、phonopy、rdkit、jupyterlab、dos2unix）用 conda 装，只在 PyPI 上更新及时的（mp-api、sumo、pyprocar、py4vasp、cclib、goodvibes、lobsterpy、dpdispatcher）放在 `pip:` 段。也可以手动一步步来：

```bash
mamba create -n dft -c conda-forge python=3.12 numpy scipy matplotlib-base pandas h5py \
      ase pymatgen spglib seekpath phonopy rdkit jupyterlab dos2unix
conda activate dft
pip install mp-api sumo pyprocar py4vasp "cclib==1.8.1" goodvibes lobsterpy dpdispatcher
```

> [!IMPORTANT]
> **先 conda、后 pip，而且 pip 只在激活的环境里用。** 先用 conda 装好带编译扩展的基础包，再用 pip 补齐剩下的；反过来容易让 pip 覆盖 conda 的 numpy，出现 `numpy.dtype size changed` 之类的二进制不兼容错误。永远不要用 `sudo pip`，也不要在 base 环境里 `pip install`。装完后用 `pip check` 检查依赖冲突。

各包在本教程中的用途：

| 包 | 用途 | 相关章节 |
|---|---|---|
| ase | 结构读写/建模、计算器驱动四款软件、NEB、ase db | 第3章、第4章 |
| pymatgen + mp-api | 结构对象、VASP 输入集、输出解析、相图、MP 数据库 | 第3章、第4章 |
| spglib / seekpath | 对称性、原胞/标准胞、高对称 k 路径 | 第3章、第5章 |
| phonopy | 声子谱、热力学性质 | [VASP 第12章](../../VASP/tutorials/12-phonon.md) |
| sumo / pyprocar / py4vasp | 能带、DOS、投影能带作图 | 第5章 |
| cclib / goodvibes | Gaussian 输出解析、准谐自由能 | 第4章、[Gaussian 第6章](../../Gaussian/tutorials/06-freq-thermo.md) |
| lobsterpy | COHP 成键分析 | [VASP 第16章](../../VASP/tutorials/16-charge-bonding.md) |
| rdkit | SMILES → 3D 分子结构 | 第3章 |
| dpdispatcher | 批量生成/提交作业脚本 | 第2章 |

### 7.5 无外网的集群：conda-pack 打包迁移

很多集群的计算节点甚至登录节点不能访问外网。做法是在能联网的 Linux 机器（或 WSL）上建好同名环境，用 `conda-pack` 打包后上传：

```bash
# 在能联网的 Linux/WSL 上（CPU 架构和 glibc 版本要与集群兼容）
mamba install -n base -c conda-forge conda-pack
conda pack -n dft -o dft.tar.gz
rsync -avP dft.tar.gz hpc:/work/your_name/

# 在集群上
mkdir -p /work/your_name/envs/dft
tar -xzf /work/your_name/dft.tar.gz -C /work/your_name/envs/dft
source /work/your_name/envs/dft/bin/activate
conda-unpack                         # 修正环境中写死的路径，只需执行一次
python -c "import pymatgen.core, ase; print('ok')"
```

### 7.6 在集群上用 Jupyter

不要在登录节点长时间开 Jupyter。正确做法是交互式申请一个计算节点（`salloc`/`srun --pty`，见第2章），在节点上 `jupyter lab --no-browser --port=8888`，再在本地做两级端口转发（具体命令因集群而异，很多中心提供 Open OnDemand 网页门户，可以一键启动 Jupyter）。

---

## 8. 环境变量：让工具找到赝势、程序和 API key

### 8.1 总表

| 变量 | 谁使用 | 示例值 | 说明 |
|---|---|---|---|
| `VASP_PP_PATH` | ASE | `$HOME/vasp_pp` | **包含** `potpaw_PBE/`、`potpaw/`（LDA）的上一级目录 |
| `ASE_VASP_COMMAND` | ASE | `mpirun -np 8 vasp_std` | ASE 运行 VASP 的命令；也可用 `VASP_COMMAND`、`VASP_SCRIPT` |
| `PMG_VASP_PSP_DIR` | pymatgen | `$HOME/pmg_potcars` | 由 `pmg config -p` 生成的目录 |
| `PMG_DEFAULT_FUNCTIONAL` | pymatgen | `PBE_54` | 默认 POTCAR 版本 |
| `MP_API_KEY` | mp-api | 32 位字符串 | Materials Project 新 API key |
| `PMG_MAPI_KEY` | pymatgen 的 MPRester | 同上 | pymatgen 自带客户端读这个名字 |
| `ASE_GAUSSIAN_COMMAND` | ASE | `g16 < PREFIX.com > PREFIX.log` | `PREFIX` 会被替换为 label |
| `ASE_CP2K_COMMAND` | ASE | `mpirun -np 4 cp2k.psmp -s` | CP2K 以 shell 模式运行 |
| `CP2K_DATA_DIR` | CP2K、ASE | `$HOME/cp2k/data` | 基组、赝势、dftd3.dat 所在目录 |
| `OPENMX_DFT_DATA_PATH` | ASE | `$HOME/openmx4.0/DFT_DATA19` | OpenMX 的 PAO/VPS 数据库 |
| `g16root`、`GAUSS_SCRDIR` | Gaussian | `/opt/gaussian`、`/tmp/$USER/g16` | 需 `source $g16root/g16/bsd/g16.profile` |
| `OMP_NUM_THREADS` | 所有 OpenMP 程序 | 纯 MPI 时为 `1` | 见第2章 |
| `OMP_STACKSIZE` | Multiwfn、OpenMP 程序 | `200M` | 太小会段错误 |

把这些写进 `~/.bashrc` 的模板见 [bashrc_dft_snippet.sh](../examples/01-env-setup/bashrc_dft_snippet.sh)。修改后执行 `source ~/.bashrc`，并用 `echo $VASP_PP_PATH`、`env | grep -E "VASP|PMG|MP_"` 检查。

### 8.2 VASP 赝势目录（ASE 用）

VASP 赝势库（POTCAR）随 VASP 许可证分发，**不能**放进公开仓库或分享给没有许可证的人。解压后整理成如下结构：

```text
$HOME/vasp_pp/                   <- VASP_PP_PATH 指向这里
├── potpaw_PBE/                  <- PBE（可以是 potpaw_PBE.64 的软链接）
│   ├── Si/POTCAR
│   ├── Fe_pv/POTCAR
│   └── ...
└── potpaw/                      <- LDA（ASE 约定的名字）
```

ASE 按 `$VASP_PP_PATH/potpaw_PBE/<元素+后缀>/POTCAR` 拼接 POTCAR。若你的目录名带版本号（如 `potpaw_PBE.64`），可以建软链接 `ln -s potpaw_PBE.64 potpaw_PBE`，或者设置 `VASP_PP_VERSION=64`（ASE 会去找 `potpaw_PBE.64`）。验证方法（只写输入，不运行 VASP）：

```bash
python - <<'EOF'
from ase.build import bulk
from ase.calculators.vasp import Vasp
si = bulk("Si", "diamond", a=5.43)
calc = Vasp(xc="PBE", encut=400, kpts=(8, 8, 8), gamma=True, directory="test_ase_vasp")
calc.write_input(si)            # 生成 INCAR/POSCAR/KPOINTS/POTCAR
EOF
head -1 test_ase_vasp/POTCAR     # 应显示 PAW_PBE Si ...
```

如果没设置 `VASP_PP_PATH`，会报找不到赝势；运行计算时若没设置运行命令，会报 `Please set either command in calculator or one of the following environment variables (prioritized as follows): ASE_VASP_COMMAND, VASP_COMMAND, VASP_SCRIPT`。

> [!NOTE]
> 较新的 ASE 推荐把这些设置写进配置文件 `~/.config/ase/config.ini` 的 `[environment]` 段（模板见 [ase_config.ini.example](../examples/01-env-setup/ase_config.ini.example)，路径可用 `ASE_CONFIG_PATH` 改）。ASE 3.29 仍然读取同名环境变量，但会给出 `ASEEnvDeprecationWarning: Loaded ASE_VASP_COMMAND from environment. Please use configfile.` 的警告；用 `ase info --calculators` 可以查看 ASE 实际读到的配置。

### 8.3 pymatgen 的 POTCAR 与 API key：`pmg config`

pymatgen 需要把 POTCAR 重新整理成自己的目录结构，只需做一次：

```bash
# 第一个参数：包含 potpaw_PBE、potpaw_PBE.54、potpaw_PBE.64 等子目录的原始目录
# 第二个参数：pymatgen 格式的输出目录（会生成 POT_GGA_PAW_PBE、POT_GGA_PAW_PBE_54、POT_PAW_PBE_64 ...）
pmg config -p ~/vasp_pp ~/pmg_potcars
pmg config --add PMG_VASP_PSP_DIR ~/pmg_potcars
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54       # 可选
pmg config --add PMG_MAPI_KEY <你的32位key>           # 可选：pymatgen 自带的 MPRester 用
```

`pmg config --add` 把设置写入 `~/.config/.pmgrc.yaml`（新版位置；旧位置 `~/.pmgrc.yaml` 仍作为后备读取），同名环境变量（`PMG_VASP_PSP_DIR` 等）优先级更高。模板见 [pmgrc.yaml.example](../examples/01-env-setup/pmgrc.yaml.example)。检查：

```bash
python -c "from pymatgen.core import SETTINGS; print(SETTINGS.get('PMG_VASP_PSP_DIR'))"
python -c "from pymatgen.io.vasp import PotcarSingle; p = PotcarSingle.from_symbol_and_functional('Si', 'PBE_54'); print(p.header, p.nelectrons)"
```

> [!WARNING]
> 网上老教程中的 `from pymatgen import Structure` 早已失效（报 `ImportError: cannot import name 'Structure' from 'pymatgen'`），一律写 `from pymatgen.core import Structure`；pymatgen-core 拆分后写法不变。

### 8.4 Materials Project 新 API key

Materials Project 的旧版 REST API（`materialsproject.org/rest/v2`）已于 **2025 年 9 月下线**，旧的 16 位 key 全部失效。新 key 为 32 位，在登录后的 dashboard 页面查看：

```bash
export MP_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # mp-api 读取
export PMG_MAPI_KEY="$MP_API_KEY"                      # pymatgen.ext.matproj.MPRester 读取
unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT                # 删除指向旧 v2 地址的残留设置
python -c "from mp_api.client import MPRester; m=MPRester(); print(m.get_structure_by_material_id('mp-149'))"
```

`pymatgen.ext.matproj.MPRester` 在 key 长度不是 32 时会直接报 `ValueError: Invalid or old API key. Please obtain an updated API key at https://materialsproject.org/dashboard.`。用法详见 第3章。

> [!CAUTION]
> API key 是个人凭证。不要写进会上传到 GitHub 的脚本或 notebook；放在 `~/.bashrc` 或 `~/.config/.pmgrc.yaml` 并 `chmod 600`。

### 8.5 VASPKIT：`~/.vaspkit`

VASPKIT 1.5.1（2024-01）仍是最新公开版本。解压后把自带的 `how_to_set_environment_variable` 复制为 `~/.vaspkit`，通常只需要改 `PBE_PATH`（指向 `potpaw_PBE`）和 `PYTHON_BIN`（指向 `dft` 环境里的 python，作图需要 numpy/scipy/matplotlib）：

```bash
tar -zxvf vaspkit.1.5.1.linux.x64.tar.gz -C ~/software/      # 包名以实际下载为准
cp ~/software/vaspkit.1.5.1/how_to_set_environment_variable ~/.vaspkit
vi ~/.vaspkit            # 修改 PBE_PATH、PYTHON_BIN（可参考本仓库的 vaspkit_rc.example）
echo 'export PATH=$HOME/software/vaspkit.1.5.1/bin:$PATH' >> ~/.bashrc
which python             # conda activate dft 之后，把这个路径填到 PYTHON_BIN
```

[vaspkit_rc.example](../examples/01-env-setup/vaspkit_rc.example) 按官方模板的键名整理。VASPKIT 1.5.0 起，如果 `~/.vaspkit` 不存在会读取 `/etc/vaspkit`（管理员统一配置）；新版本若引入了新变量，需要手动把新模板中的条目补进你的 `~/.vaspkit`。VASPKIT 的具体功能见 第5章。

### 8.6 Gaussian、CP2K、OpenMX、Multiwfn

```bash
# Gaussian 16：g16root 是包含 g16/ 的目录
export g16root=/opt/gaussian
source $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR=/tmp/$USER/g16scr && mkdir -p $GAUSS_SCRDIR   # 必须存在且空间充足

# CP2K：数据目录（BASIS_MOLOPT_UZH、POTENTIAL_UZH、dftd3.dat ...）
export CP2K_DATA_DIR=$HOME/cp2k/data

# OpenMX：输入中的 DATA.PATH 要指向 DFT_DATA19（含 PAO/ 和 VPS/ 的那一层）
export OPENMX_DFT_DATA_PATH=$HOME/openmx4.0/DFT_DATA19

# Multiwfn：大体系需要大栈空间
export Multiwfnpath=$HOME/software/Multiwfn
export PATH=$Multiwfnpath:$PATH
export OMP_STACKSIZE=200M
ulimit -s unlimited
```

---

## 9. module 系统：编译器、MPI 与程序版本要对应

### 9.1 常用命令

集群上的软件通常用 Environment Modules 或 Lmod 管理：

```bash
module avail                     # 列出可用模块（Lmod 下 module spider vasp 可搜索所有层级）
module spider vasp/6.5.1         # Lmod：显示加载 vasp/6.5.1 之前需要先加载哪些模块
module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl vasp/6.5.1
module list                      # 当前已加载的模块
module show vasp/6.5.1           # 该模块修改了哪些 PATH/LD_LIBRARY_PATH
module purge                     # 卸载全部（作业脚本开头建议先 purge）
```

模块名完全取决于集群，以上只是示例。

### 9.2 为什么版本必须对应

`vasp_std`、`cp2k.psmp`、`openmx` 这类 MPI 程序在编译时链接了某个特定的编译器运行库、MPI 库和数学库（MKL/OpenBLAS/ScaLAPACK）。运行时：

- **`mpirun` 必须来自编译时用的那套 MPI**（Intel MPI 编译的程序用 Open MPI 的 `mpirun` 启动，或者反过来，都会出问题）；
- 动态库要能找到（`LD_LIBRARY_PATH` 由 module 设置）。

检查方法：

```bash
which mpirun && mpirun --version                 # 当前 mpirun 来自哪套 MPI
ldd $(which vasp_std) | grep -i -E "mpi|mkl|scalapack|fftw|not found"
```

典型症状对照：

| 现象 | 原因 | 解决 |
|---|---|---|
| `error while loading shared libraries: libmkl_intel_lp64.so.2: cannot open shared object file` | 没加载编译时的 MKL/oneAPI 模块 | `module load` 对应模块；`ldd` 中不应有 `not found` |
| 程序输出里同样的启动信息重复 N 遍、每份都显示只有 1 个进程（CP2K 中表现为多个 `PROGRAM STARTED AT` 和 `DBCSR\| MPI: Number of processes 1`） | `mpirun` 与程序链接的 MPI 不是同一套，N 个进程各自独立串行运行 | 用与编译时相同的 MPI 模块；检查 `which mpirun` |
| `conda activate dft` 之后 MPI 程序启动失败或变慢 | conda 环境里装了 `openmpi`/`mpich`（例如作为 mpi4py 的依赖），它的 `mpirun` 排到了 PATH 最前面 | 作业脚本中不激活 conda，或在 `module load` 之后再确认 `which mpirun` |
| `forrtl: severe (174): SIGSEGV, segmentation fault occurred` | 栈空间不足，或 MPI/OpenMP 混用不当，或编译器运行库不匹配 | `ulimit -s unlimited`；纯 MPI 时 `export OMP_NUM_THREADS=1`；核对模块版本 |

> [!TIP]
> 把“能跑通”的模块组合原样写进作业脚本（见第2章模板），并在计算记录里保存 `module list` 的输出。半年后复现结果时，这比任何记忆都可靠。

---

## 10. 换行符与编码：dos2unix

Windows 文本文件每行以 `\r\n`（CRLF）结尾，Linux 用 `\n`（LF）。在 Windows 记事本、Word、某些编辑器里写的作业脚本和输入文件传到集群后，会出现下面这些看似奇怪的错误：

| 报错 / 现象 | 原因 |
|---|---|
| `sbatch: error: Batch script contains DOS line breaks (\r\n)` 和 `sbatch: error: instead of expected UNIX line breaks (\n).` | 作业脚本是 CRLF |
| `/bin/bash^M: bad interpreter: No such file or directory` | 第一行 `#!/bin/bash\r`，系统去找名为 `bash\r` 的程序 |
| `$'\r': command not found` | 脚本中的空行实际是 `\r` |
| INCAR 中某个参数“没生效”（OUTCAR 中显示的仍是默认值） | 该行带 `\r` 或 Tab，VASP 没能正确解析 |
| Gaussian 读坐标或基组时出错、CP2K 关键字无法识别 | 行尾 `\r` 被当作值的一部分 |

检查与修复：

```bash
file job.sh INCAR                 # 显示 "with CRLF line terminators" 即有问题
cat -A INCAR | head               # 行尾显示 ^M$ 的就是 CRLF；^I 是 Tab
dos2unix job.sh INCAR KPOINTS POSCAR *.gjf *.inp *.dat   # 原地转换为 LF
sed -i 's/\r$//' job.sh           # 没有 dos2unix 时的替代
grep -nP '[^\x00-\x7F]' INCAR     # 找出非 ASCII 字符（中文、全角等号/逗号、不间断空格）
```

> [!IMPORTANT]
> 输入文件（INCAR、KPOINTS、POSCAR、`*.gjf`、`*.inp`、`*.dat`）**只用 ASCII 字符和英文注释**。从网页或 Word 复制过来的“＝”“，”、中文引号、不间断空格都可能让程序静默忽略参数或直接报错。

预防胜于修复：

- VS Code：右下角状态栏点击 `CRLF` 改为 `LF`；在设置中把 `files.eol` 设为 `\n`。Notepad++：编辑 → 文档格式转换 → 转为 Unix (LF)。
- Git：在仓库根目录 `.gitattributes` 写 `* text=auto eol=lf`（本仓库已经这样配置），或 `git config --global core.autocrlf input`。
- 本仓库第2章的作业模板在运行前会自动检查输入文件里有没有 CRLF。

---

## 11. 可复现：用 environment.yml 锁定版本并自检

### 11.1 三种导出方式

```bash
conda env export --from-history > environment.yml   # 只记录你显式安装的包（跨平台、易读）
conda env export > environment-full.yml              # 记录所有包的精确版本和 build（同平台精确复现）
conda list --explicit > spec-file.txt                # 带下载 URL 的显式清单（同平台、可离线复现）
pip freeze > requirements-pip.txt                    # pip 安装部分的精确版本
```

| 方式 | 可读性 | 跨平台 | 精确程度 | 用途 |
|---|---|---|---|---|
| 手写 `environment.yml`（带 `>=` 下限） | 高 | 是 | 低 | 教程、课题组共享“推荐环境” |
| `--from-history` | 高 | 是 | 中（只含你指定的版本约束） | 日常备份 |
| `conda env export` 全量 | 低 | 否 | 高 | 论文附带的精确环境 |
| `conda list --explicit` | 低 | 否 | 最高 | 归档、离线重建 `conda create -n x --file spec-file.txt` |

> [!TIP]
> 写论文时，在 SI 里注明关键包的版本（例如 “pymatgen 2026.10.24, phonopy 4.10.0, ASE 3.29.0”），并把 `environment-full.yml` 和计算输入一起归档。版本变化会影响结果：例如 ASE 3.28 起 NEB 默认方法改为 `improvedtangent`，pymatgen 的输入集默认参数也会随版本调整。

更新环境：修改 `environment.yml` 后执行 `mamba env update -f environment.yml --prune`（`--prune` 删除文件中已不存在的包）。

### 11.2 自检

最简单的检查：

```bash
conda activate dft
python --version
python -c 'import ase; print(ase.__version__)'
python -c 'import pymatgen.core as p; print(p.__version__)'
python -c 'import phonopy; print(phonopy.__version__)'
pip list | grep -i -E "^(ase|pymatgen|phonopy|sumo|pyprocar|py4vasp|cclib|goodvibes|lobsterpy|mp-api) "
pip check                                     # 依赖是否一致
python Common/examples/01-env-setup/print_versions.py
```

更完整的检查用本仓库的 env_check.sh：它会逐项检查 Python 版本和关键包版本、`VASP_PP_PATH`/`PMG_VASP_PSP_DIR`/`MP_API_KEY` 是否设置、`~/.vaspkit` 是否存在、`ulimit -s`、`OMP_NUM_THREADS` 和 Multiwfn 所需的 `OMP_STACKSIZE`，并给出修复建议。运行 `bash Common/scripts/env_check.sh`，可用选项见脚本开头的说明或 `-h`。

（示意输出，数值仅供参考）

```text
Python   3.12.7 (/home/your_name/miniforge3/envs/dft/bin/python)
--------------------------------------------
ase              3.29.0
pymatgen         2026.10.24
pymatgen-core    2026.10.31
spglib           2.8.0
phonopy          4.10.0
mp-api           0.47.0
...
```

---

## 推荐配置速查

| 项目 | 推荐 | 理由 |
|---|---|---|
| Windows 用户 | WSL2 + Ubuntu 24.04，计算放在 `~` 而非 `/mnt/c` | 真 Linux 内核，I/O 快，避免 CRLF |
| SSH 密钥 | ed25519，设置口令，私钥 600 权限 | 安全、短、现代 OpenSSH 默认支持 |
| 登录配置 | `~/.ssh/config` 起别名 + `ServerAliveInterval 60` | 少打字、不掉线 |
| 文件同步 | `rsync -avzP`，排除 WAVECAR/CHG*/`*.wfn` | 增量、可续传、省流量 |
| 会话 | tmux | 断线不丢操作；记住所在登录节点 |
| conda 发行版 | Miniforge（conda-forge 频道） | 无商业条款限制，科学包更新快，自带 mamba |
| Python | 3.12 | 满足 pymatgen/phonopy（≥3.11）、mp-api（≥3.12）、lobsterpy（<3.15） |
| 安装顺序 | 先 conda 后 pip，不用 base 环境 | 避免二进制不兼容 |
| cclib | 固定 `cclib==1.8.1` | 2.0a1 解析 Gaussian 有 bug |
| 环境变量 | 写入 `~/.bashrc`（作业需要的放在交互判断之前） | 批处理作业也能读到 |
| 版本记录 | `environment.yml` + `conda env export` + `module list` | 结果可复现 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| `Permission denied (publickey)` | 公钥没加到 `authorized_keys`；集群端 `~/.ssh` 或家目录权限过宽；用错私钥 | `ssh -v hpc` 看尝试了哪把钥匙；修正权限（第 4.2 节）；在 config 中写 `IdentityFile` |
| `WARNING: UNPROTECTED PRIVATE KEY FILE!` / `Permissions 0644 for '.../id_ed25519' are too open.` | 私钥权限过宽（WSL 中把 `/mnt/c` 下的钥匙拿来用时常见） | `chmod 600 ~/.ssh/id_ed25519`；钥匙放在 WSL 自己的 `~/.ssh` |
| `WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED!` | 集群重装/换了主机密钥，或存在中间人风险 | 先向管理员确认，再 `ssh-keygen -R login.hpc.example.edu` |
| `ssh: connect to host ... port 22: Connection timed out` | 不在校园网/VPN、端口不对、防火墙 | 连 VPN；确认端口；通过跳板机 `ProxyJump` |
| 连接空闲一会儿就断 | 防火墙回收空闲连接 | `ServerAliveInterval 60`；用 tmux 防止丢失工作 |
| `conda: command not found` | 没执行 `conda init` 或没 `source ~/.bashrc` | `~/miniforge3/bin/conda init bash && source ~/.bashrc` |
| `CondaHTTPError` / 下载极慢 | 网络问题 | 配置镜像（condarc_mirror.example）；或在其他机器上 conda-pack |
| `pip install phonopy` 装上的是 2.x 旧版 / 报 `Requires-Python >=3.11` | 环境的 Python 太旧 | 新建 Python 3.12 环境 |
| `ImportError: cannot import name 'Structure' from 'pymatgen'` | 旧写法 | `from pymatgen.core import Structure` |
| `No POTCAR for Fe_pv with functional PBE found. Please set the PMG_VASP_PSP_DIR ...`（新版为 `OSError: You do not have the right POTCAR ...`） | 未配置 pymatgen 的 POTCAR 目录或缺少该版本 | `pmg config -p` + `pmg config --add PMG_VASP_PSP_DIR`；需要 PBE_54 时设置 `PMG_DEFAULT_FUNCTIONAL` |
| `Please set either command in calculator or one of the following environment variables ... ASE_VASP_COMMAND, VASP_COMMAND, VASP_SCRIPT` | ASE 不知道如何运行 VASP | `export ASE_VASP_COMMAND='mpirun -np 8 vasp_std'` 或写进 config.ini |
| ASE 报 `VASP_PP_PATH does not exist` 或找不到 POTCAR | 路径写错，或指向了 `potpaw_PBE` 本身 | `VASP_PP_PATH` 应指向**包含** `potpaw_PBE` 的上一级目录 |
| `ValueError: Invalid or old API key` / 403 | 用了旧版 MP key 或旧地址 | 到 dashboard 取 32 位新 key；`unset MP_API_ENDPOINT` |
| `sbatch: error: Batch script contains DOS line breaks (\r\n)` | 脚本在 Windows 下编辑 | `dos2unix job.sh` |
| `/bin/bash^M: bad interpreter` | 同上 | 同上 |
| VASPKIT 作图功能报错找不到 matplotlib | `PYTHON_BIN` 指向了系统 Python | 改为 `dft` 环境中 `which python` 的路径 |
| `error while loading shared libraries: ...so: cannot open shared object file` | 运行时没加载编译时的模块 | `module load` 对应模块，`ldd` 检查 |
| `forrtl: severe (174): SIGSEGV` | 栈空间不足 / MPI 与 OpenMP 混用不当 | `ulimit -s unlimited`，`OMP_NUM_THREADS=1`（纯 MPI） |
| `Disk quota exceeded` | 家目录配额满（常见于 conda 缓存或 WAVECAR） | `du -sh ~/* ~/.conda ~/miniforge3/pkgs`；`conda clean -a`；计算移到工作盘 |

---

## 练习

1. 在 WSL2（或你的 Linux/macOS）上生成 ed25519 密钥，配置 `~/.ssh/config` 中名为 `hpc` 的别名，实现 `ssh hpc` 免密登录。
2. 用 `environment.yml` 建好 `dft` 环境，运行 `print_versions.py`，把输出保存为 `versions.txt`。
3. 按第 8.2 节组织一个 VASP 赝势目录（有许可证的前提下），用 ASE 为 Si 写出输入文件并检查 POTCAR 第一行；再用 `pmg config` 配置 pymatgen，比较两者生成的 POTCAR 是否一致（`diff`）。
4. 在 Windows 记事本里写一个两行的 `test.sh`，传到 Linux 后执行 `bash test.sh` 和 `sbatch test.sh`（若有集群），观察报错；再用 `file`、`cat -A`、`dos2unix` 诊断并修复。
5. 在集群上建一个包含若干子目录的测试目录，每个子目录放一个大文件 `WAVECAR`（可用 `dd if=/dev/zero of=WAVECAR bs=1M count=100` 生成）和一个小文件 `OUTCAR`，用 `rsync --include/--exclude` 只把 OUTCAR 拉回本地。

---

## 参考资料

官方文档：

- WSL 安装（Microsoft）：<https://learn.microsoft.com/windows/wsl/install>
- WSL 高级设置 `.wslconfig`：<https://learn.microsoft.com/windows/wsl/wsl-config>
- OpenSSH 手册（`ssh_config(5)`）：<https://man.openbsd.org/ssh_config>
- rsync 手册：<https://download.samba.org/pub/rsync/rsync.1>
- tmux Wiki：<https://github.com/tmux/tmux/wiki>
- Miniforge：<https://github.com/conda-forge/miniforge>
- conda 环境管理：<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>
- conda-pack：<https://conda.github.io/conda-pack/>
- ASE 文档（2026-06 起迁到 docs.ase-lib.org）：<https://docs.ase-lib.org/>；计算器配置文件说明：<https://docs.ase-lib.org/ase/calculators/calculators.html>
- ASE 版本说明：<https://docs.ase-lib.org/releasenotes.html>
- pymatgen 文档与更新日志：<https://pymatgen.org/>，<https://pymatgen.org/CHANGES.html>
- Materials Project API 入门（新 key）：<https://docs.materialsproject.org/downloading-data/using-the-api/getting-started>
- VASPKIT 安装与 `~/.vaspkit` 说明：<https://vaspkit.com/installation.html>；更新日志：<https://vaspkit.com/changelog.html>
- phonopy 文档：<https://phonopy.github.io/phonopy/>
- Gaussian 16 Link 0 命令与环境变量：<https://gaussian.com/link0/>

相关 GitHub 仓库（只链接，不复制代码）：

- [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) — 本教程的 Python 基础库，2026 年核心拆分为 pymatgen-core。
- [ase/ase（GitLab）](https://gitlab.com/ase/ase) — ASE 官方仓库在 GitLab，统一驱动四款软件。
- [materialsproject/api](https://github.com/materialsproject/api) — mp-api 客户端，旧 API 已于 2025-09 下线。
- [phonopy/phonopy](https://github.com/phonopy/phonopy) — 4.x 起要求 Python ≥ 3.11。
- [dstdev/awesome-hpc](https://github.com/dstdev/awesome-hpc) — HPC 工具与资料的 awesome 列表，延伸阅读。
- [northword/wiki](https://github.com/northword/wiki) — 中文学习笔记，含 Linux 基础与 PBS 作业管理入门。
- [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script) — 中文社区常用 VASP 小脚本（含 CRLF 修正），无许可证，仅供参考链接。

[返回目录](../README.md)
