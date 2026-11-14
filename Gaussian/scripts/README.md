# Gaussian 辅助脚本（Gaussian helper scripts）

本目录收录 Gaussian 16 计算的前处理、作业提交和后处理脚本，覆盖一次计算的全过程：生成输入（`gjf_builder.py`），提交作业（`sub_g16.sh`），监控优化（`opt_monitor.py`），诊断报错（`g16_error_doctor.py`），批量汇总结果（`g16_log_summary.py`），提取结构（`log2xyz.py`），画扫描/IRC 能量曲线（`scan_irc_profile.py`），计算 quasi-RRHO 热化学（`thermo_qrrho.py`），画 IR/Raman 与 UV-Vis/ECD 光谱（`spectrum_ir_raman.py`、`spectrum_uvvis_ecd.py`），换算 NMR 化学位移（`nmr_shift.py`），以及批量生成 cube 文件、调用 Multiwfn（`batch_cube_multiwfn.sh`）。解析规则按 Gaussian 16 Rev. C.01/C.02 的输出格式编写。

> [!NOTE]
> 所有 Python 脚本都支持 `python <脚本> -h` 查看完整参数，两个 Bash 脚本支持 `bash <脚本> -h`。本 README 的参数都取自各脚本的 argparse 定义或 usage 说明。命令行示例只演示用法，不附运行结果；文中出现的文件名（`opt.log`、`ts.gjf` 等）都是占位符。

## 总览

| 脚本 | 功能 | 依赖 | 对应章节 |
|---|---|---|---|
| [gjf_builder.py](gjf_builder.py) | 用任务模板（sp/opt/freq/optfreq/ts/qst2/qst3/irc/scan/td/tdopt/nmr/cp/oniom），从 xyz、log 或 gjf 批量生成输入，可加溶剂、色散、--Link1--、GenECP，并检查电荷、多重度与电子数是否匹配 | 仅标准库 | [第2章](../tutorials/02-input-anatomy.md)、[第4章](../tutorials/04-sp-scf.md)、[第7章](../tutorials/07-ts-irc.md)、[第8章](../tutorials/08-scan.md)、[第9章](../tutorials/09-solvation.md)、[第14章](../tutorials/14-weak-bsse.md)、[第15章](../tutorials/15-oniom.md)、[第17章](../tutorials/17-hpc-jobs.md) |
| [sub_g16.sh](sub_g16.sh) | 向 Slurm/PBS 提交 g16 作业：核数与内存保持一致，每个作业独立 scratch，自动 formchk，失败时以非零码退出 | bash、Slurm 或 PBS、Gaussian 16 | [第1章](../tutorials/01-setup-first-job.md)、[第7章](../tutorials/07-ts-irc.md)、[第15章](../tutorials/15-oniom.md)、[第17章](../tutorials/17-hpc-jobs.md) |
| [g16_error_doctor.py](g16_error_doctor.py) | 诊断失败或可疑的作业：识别出错的 Link，匹配约 25 类常见报错，给出中文原因与处理建议，可起草续算 gjf | 仅标准库（`--draft` 需同目录的 log2xyz.py） | [第1章](../tutorials/01-setup-first-job.md)、[第4章](../tutorials/04-sp-scf.md)、[第7章](../tutorials/07-ts-irc.md)、[第9章](../tutorials/09-solvation.md)、[第13章](../tutorials/13-openshell-stability.md)、[第18章](../tutorials/18-errors-capstone.md) |
| [g16_log_summary.py](g16_log_summary.py) | 批量汇总 .log/.out：结束状态、能量、优化收敛、虚频、ZPE/H/G、<S**2>、耗时；输出为表格、CSV 或 JSON | 仅标准库；cclib 可选 | [第1章](../tutorials/01-setup-first-job.md)、[第3章](../tutorials/03-method-basis.md)、[第4章](../tutorials/04-sp-scf.md)、[第7章](../tutorials/07-ts-irc.md)、[第10章](../tutorials/10-tddft.md)、[第13章](../tutorials/13-openshell-stability.md)、[第17章](../tutorials/17-hpc-jobs.md)、[第18章](../tutorials/18-errors-capstone.md) |
| [log2xyz.py](log2xyz.py) | 从 opt/scan/IRC 输出提取结构（最后一个、能量最低、全部、扫描收敛点、IRC 路径），写成 xyz，可生成续算 gjf | 仅标准库 | [第1章](../tutorials/01-setup-first-job.md)、[第5章](../tutorials/05-opt.md)、[第7章](../tutorials/07-ts-irc.md)、[第8章](../tutorials/08-scan.md)、[第11章](../tutorials/11-nmr-props.md)、[第16章](../tutorials/16-postprocess-multiwfn.md)、[第18章](../tutorials/18-errors-capstone.md) |
| [opt_monitor.py](opt_monitor.py) | 逐步解析几何优化的能量与四项收敛判据，诊断振荡、缓慢下降或平坦势能面，并画判据与阈值对比图 | 标准库；作图需 matplotlib | [第5章](../tutorials/05-opt.md)、[第7章](../tutorials/07-ts-irc.md)、[第8章](../tutorials/08-scan.md)、[第18章](../tutorials/18-errors-capstone.md) |
| [scan_irc_profile.py](scan_irc_profile.py) | 一维/二维柔性扫描、刚性扫描和 IRC 的能量曲线；正反向 IRC 自动合并；可导出最高点结构作为 TS 初猜 | 标准库；作图需 matplotlib | [第7章](../tutorials/07-ts-irc.md)、[第8章](../tutorials/08-scan.md)、[第16章](../tutorials/16-postprocess-multiwfn.md) |
| [thermo_qrrho.py](thermo_qrrho.py) | 教学版热化学重算：Grimme qRRHO 熵、Head-Gordon 焓校正、Truhlar 准谐振、1 M 标准态校正、频率缩放、高精度单点组合、构象 Boltzmann 加权 | numpy | [第6章](../tutorials/06-freq-thermo.md)、[第7章](../tutorials/07-ts-irc.md)、[第9章](../tutorials/09-solvation.md)、[第18章](../tutorials/18-errors-capstone.md) |
| [spectrum_ir_raman.py](spectrum_ir_raman.py) | 从频率计算画 IR/Raman 谱：频率缩放、洛伦兹或高斯展宽、Raman 活性换算为强度、多个体系叠加对比 | numpy、matplotlib | [第6章](../tutorials/06-freq-thermo.md)、[第16章](../tutorials/16-postprocess-multiwfn.md) |
| [spectrum_uvvis_ecd.py](spectrum_uvvis_ecd.py) | 从 TD-DFT/TDA/CIS 结果画 UV-Vis 与 ECD 谱：高斯展宽、能量平移、对映体翻转、构象 Boltzmann 加和、叠加实验谱 | numpy、matplotlib | [第10章](../tutorials/10-tddft.md)、[第16章](../tutorials/16-postprocess-multiwfn.md) |
| [nmr_shift.py](nmr_shift.py) | 把 GIAO/CSGT 各向同性屏蔽常数换算成化学位移（参比法或线性回归缩放），做等价核平均与构象 Boltzmann 平均，统计与实验值的误差 | 仅标准库 | [第11章](../tutorials/11-nmr-props.md) |
| [batch_cube_multiwfn.sh](batch_cube_multiwfn.sh) | 批量 formchk、读取 HOMO/LUMO 能量、用 cubegen 生成轨道/密度/ESP/自旋 cube，并按命令流文件批量运行 Multiwfn | bash、formchk/cubegen；Multiwfn 可选 | [第10章](../tutorials/10-tddft.md)、[第12章](../tutorials/12-population-charges.md)、[第14章](../tutorials/14-weak-bsse.md)、[第16章](../tutorials/16-postprocess-multiwfn.md) |

## 安装依赖与通用约定

```bash
# Python >= 3.9
pip install numpy matplotlib   # thermo_qrrho（numpy）、两个光谱脚本（numpy + matplotlib）
                               # opt_monitor / scan_irc_profile 只有作图时才需要 matplotlib
pip install cclib              # 可选：g16_log_summary.py --cclib 交叉验证
# gjf_builder / log2xyz / g16_log_summary / g16_error_doctor / nmr_shift 只用标准库
```

通用约定：

- **能量单位**：脚本读入的能量均为 Hartree；相对能量默认用 kcal/mol（1 Eh = 627.509474 kcal/mol），部分脚本可以改用 kJ/mol 或 eV（见各节说明）。
- **能量取哪一个**：多数脚本按"最具体的能量"取值，优先级为 Counterpoise 校正 > ONIOM 外推 > CCSD(T) > MP2/双杂化 > TD 总能量 > SCF Done，具体顺序见各脚本文件头部的说明。
- **多步作业**：--Link1-- 多步作业都能处理。一般取最后一个打印了相应数据的作业步，例如最后一个频率块、最后一个激发态块、最后一个屏蔽常数块。
- **作图**：matplotlib 在用到时才导入；没有图形界面时自动用 Agg 后端保存图片，加 `--show` 才会弹出窗口。`opt_monitor.py` 和 `scan_irc_profile.py` 加 `--no-plot` 后可以完全不装 matplotlib。
- **不会覆盖输入**：会写文件的脚本都写新文件；`gjf_builder.py`、`g16_error_doctor.py --draft` 遇到同名文件时拒绝覆盖，需加 `-f/--force`（g16_error_doctor 为 `--force`）。`sub_g16.sh --rewrite` 修改输入前会先备份为 `.bak`。
- **Windows 与 Linux**：Python 脚本两边都能用；Bash 脚本要在集群或 Linux 环境中运行。在 Windows 上编辑过的 gjf 如果带 CRLF 换行，`sub_g16.sh` 会在提交前自动转换。

> [!WARNING]
> 脚本给出的诊断、续算草稿和热化学结果只是辅助判断，不能代替自己阅读输出。发表前请用标准工具（GoodVibes、Shermo、Multiwfn 等）复核，并在文中写明所用方法和参数。

---

## gjf_builder.py：批量生成输入文件

**用途**：用任务模板从 xyz、Gaussian 输出（取最后一个或指定帧的结构）或已有 gjf（只读笛卡尔坐标）生成 `.gjf`。脚本会自动写 `%chk/%mem/%nprocshared`、溶剂、色散和多步 --Link1-- 作业，可以从本地基组库插入 GenECP 块，生成前还会检查电荷、多重度与电子数是否匹配。

**依赖**：Python >= 3.9 标准库。

**任务模板**（`-t/--task`，默认 `optfreq`）：`sp`、`opt`、`freq`、`optfreq`、`ts`、`qst2`、`qst3`、`irc`、`scan`、`td`、`tdopt`、`nmr`、`cp`、`oniom`。

**常用参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-m/--method`、`-b/--basis` | `B3LYP`、`def2SVP` | 方法与基组 |
| `--disp` | `none` | `GD3BJ`、`GD3`、`GD2` 或 `none` |
| `--solvent`、`--solvent-model` | 气相、`smd` | 溶剂名；模型可选 `smd`、`pcm`、`iefpcm`、`cpcm` |
| `-c/--charge`、`-s/--mult` | 来自源文件，或电荷 0 + 最低多重度 | 电荷与自旋多重度 |
| `--route-prefix` | `#p` | `#`、`#p`、`#n`、`#t`（大小写均可） |
| `--extra` | — | 额外的 route 关键词，如 `'Int=SuperFine NoSymm'` |
| `--no-check` | — | 跳过电荷/多重度/电子数检查 |
| `--mem`、`--nproc`、`--cpu` | `8GB`、8 | `%mem`（`none` 表示不写）、`%nprocshared`（0 表示不写）；`--cpu 0-15` 改写 `%cpu` |
| `--chk`、`--no-chk`、`--oldchk` | 输出文件名 | `--oldchk` 同时加 `Guess=Read` |
| `--frame`、`--all-frames` | 最后一帧 | 选择源文件中的某一帧，或每帧生成一个输入（构象批量） |
| `--freeze`、`--modred`、`--scan` | — | 冻结原子（如 `1-10,15`），ModRedundant 行，扫描行（可重复） |
| `--maxcycles`、`--calcfc`、`--tight`、`--no-freq` | — | 修改 Opt 选项；`--no-freq` 让 ts/qst2/qst3 不附带 Freq |
| `--product`、`--ts-guess` | — | qst2/qst3 的产物与 TS 初猜（原子顺序必须一致） |
| `--irc-points`、`--irc-step`、`--irc-direction`、`--irc-integrator` | 30、10、`both`、`default` | IRC 的 MaxPoints、StepSize、方向与积分器（`lqa`/`gs2`/`hpc`） |
| `--nstates`、`--root`、`--states`、`--tda` | 10、—、`singlets` | TD-DFT 设置；`--states` 可选 `singlets`/`triplets`/`50-50` |
| `--nmr-method` | `giao` | `giao` 或 `csgt` |
| `--fragments`、`--frag-cm` | — | Counterpoise 片段（`'1-3/4-6'` 或 `auto`）以及各片段的电荷和多重度 |
| `--high`、`--medium`、`--low-method`、`--medium-method`、`--oniom-job` | —、—、`UFF`、`PM6`、`opt` | ONIOM 分层与方法 |
| `--heavy-basis`、`--heavy-elements`、`--basis-lib` | —、Kr 以后的元素、`$GJF_BASIS_LIB` | GenECP：`LANL2DZ`/`LANL2MB`/`SDD`、基组文件路径，或基组库中的名称 |
| `--link1` | — | 追加 --Link1-- 步：`'task[:method/basis[:extra words]]'`，可重复 |
| `-o/--output`、`--outdir`、`--suffix`、`--stdout`、`-f/--force` | `<stem>_<task>.gjf`，写在当前目录 | 输出文件（`-o` 只能用于单个结构）、输出目录、后缀、打印到屏幕、覆盖 |

**用法示例**：

```bash
python gjf_builder.py water.xyz -t optfreq                       # B3LYP/def2SVP Opt Freq
python gjf_builder.py water.xyz -t optfreq -m B3LYP -b "6-31G(d)" --disp GD3BJ --solvent Water
python gjf_builder.py conf.xyz --all-frames -t opt --mem 16GB --nproc 16
python gjf_builder.py opt.log -t sp -m M062X -b def2TZVP --disp GD3     # last geometry of a log
python gjf_builder.py pd.xyz -t opt -b "6-31G(d)" --heavy-basis LANL2DZ  # GenECP for Pd
python gjf_builder.py r.xyz -t optfreq --link1 "sp:wB97XD/def2TZVP"      # Opt Freq + high-level SP
python gjf_builder.py ts_guess.xyz -t ts --link1 irc
python gjf_builder.py react.xyz -t qst2 --product prod.xyz
python gjf_builder.py butane.xyz -t scan --scan "D 1 2 3 4 S 36 10.0"
python gjf_builder.py dye.xyz -t td --nstates 20 -m CAMB3LYP --solvent Ethanol
python gjf_builder.py dimer.xyz -t cp --fragments auto -b "aug-cc-pVDZ"
python gjf_builder.py enzyme_model.xyz -t oniom --high 1-12 --low-method UFF
```

**基组库**：从 [Basis Set Exchange](https://www.basissetexchange.org) 以 "Gaussian" 格式下载基组文件，放到 `--basis-lib` 指定的目录（或环境变量 `GJF_BASIS_LIB` 指向的目录），例如 `DIR/def2-TZVP.gbs`。之后 `--heavy-basis def2-TZVP` 会从该文件中取出各重元素的轨道基组，文件里有 ECP 时一并取出。

**输入/输出**：输入为 `.xyz`（可多帧）、`.gjf/.com`、`.log/.out`；输出为 ASCII、LF 换行、以空行结尾的 `.gjf`。

**局限**：不读 Z-matrix 输入；ONIOM 只生成框架，连接关系和 link atom 要自己检查；检查只是提示，不能保证输入完全正确。提交前务必通读生成的文件。

## sub_g16.sh：提交 Gaussian 作业（Slurm / PBS）

**用途**：为每个输入写一个作业脚本 `<stem>.sub_g16.sh`，放在输入旁边，然后用 `sbatch --chdir` 或 `qsub` 提交。它主要解决三个常见问题：`%nprocshared/%mem` 与调度器申请的资源不一致，scratch 文件残留，以及失败作业没有被及时发现。

**依赖**：bash；Slurm（`sbatch`）或 PBS/Torque/PBS Pro（`qsub`）；计算节点上可用的 `g16` 和 `formchk`。

**站点设置**：直接修改脚本开头标有 `USER SETTINGS` 的区块，或者设置环境变量 `SUB_G16_MODULE`（module 名）、`g16root`（默认 `/opt`）、`SCR_BASE`（scratch 父目录，默认 `/tmp`）、`SUB_G16_PARTITION`、`SUB_G16_ACCOUNT`。

**资源的两种模式**：

- **默认（从输入读取）**：读取每个作业步的 `%nprocshared`（或 `%nproc`、`%cpu`）和 `%mem`。申请的核数等于 nprocshared，内存为 `%mem / 0.8` 并向上取整到整 GB。输入里没写的 Link 0 行，用 `GAUSS_PDEF/GAUSS_MDEF` 补上（MDEF 取申请内存的 80%）。如果输入与 `-c/-m` 不一致，或 `%mem` 超过申请量的 90%，脚本会给出警告。
- **`--rewrite`（必须同时给 `-c` 和 `-m`）**：把每个 Link1 步的 `%nprocshared` 改为 `-c`，`%mem` 改为 `-m` 的 80%，并替换 `%cpu`；原文件备份为 `<input>.bak`。

**参数**：

| 参数 | 说明 |
|---|---|
| `-c/--cpus N`、`-m/--mem SIZE` | 核数；调度器内存（如 `32G`、`64000M`） |
| `-t/--time HH:MM:SS` | 墙钟时间（默认 `24:00:00`） |
| `-q/--partition`（`--queue`）、`-A/--account`、`-J/--name` | 分区/队列、账户、作业名（`-J` 只能用于单个输入） |
| `--rewrite`、`--pin` | 改写 Link 0；`--pin` 在 Slurm 下用 `GAUSS_CDEF` 把线程绑定到分配的 CPU |
| `-s/--scheduler slurm\|pbs\|pbspro` | 默认自动检测：有 `sbatch` 用 Slurm，否则有 `qsub` 用 PBS。`pbs` 写 `nodes=1:ppn=N`，`pbspro` 写 `select=1:ncpus=N:mem=...` |
| `-d/--dependency DEP` | 例如 `afterok:123456` |
| `-x/--extra "ARGS"` | 传给 sbatch/qsub 的其他参数，如 `"--exclusive"` |
| `-a/--all [DIR...]` | 提交目录中所有 `*.gjf`、`*.com`（默认当前目录） |
| `--scr-base DIR`、`--keep-scr-on-fail` | scratch 父目录；作业失败时保留 scratch |
| `--no-formchk`、`--module NAME`、`--g16root DIR` | 不运行 formchk；指定 module；指定 Gaussian 根目录 |
| `-n/--dry-run`、`-h/--help` | 只生成作业脚本、不提交；显示帮助 |

**用法示例**：

```bash
bash sub_g16.sh water.gjf                          # resources read from the input
bash sub_g16.sh -c 16 -m 32G -t 48:00:00 ts.gjf     # explicit request (warns on mismatch)
bash sub_g16.sh -c 32 -m 120G --rewrite big.gjf    # make %nprocshared/%mem match the request
bash sub_g16.sh -a -c 8 -m 16G                     # every *.gjf / *.com in this directory
bash sub_g16.sh -s pbs -q normal job.gjf           # PBS/Torque
bash sub_g16.sh -s pbspro job.gjf                  # PBS Pro
bash sub_g16.sh -d afterok:123456 freq.gjf         # start after job 123456 succeeded
bash sub_g16.sh -n job.gjf                         # dry run: inspect job.sub_g16.sh first
```

**生成的作业脚本会做这些事**：

1. 加载 module 或 source `g16.profile`，执行 `ulimit -s unlimited`；
2. 为每个作业创建独立的 `GAUSS_SCRDIR`，在退出、`scancel/qdel` 或超时（EXIT/TERM/INT）时删除；
3. 在后台运行 g16 并 `wait`，这样能及时响应信号；
4. 检查是否有 `Normal termination`，然后对输入中的每个 `%chk` 运行 formchk；
5. 失败时打印 Error termination 前 25 行，提示用 `g16_error_doctor.py` 诊断，并以非零码退出，使 `afterok` 依赖的后续作业不再运行。

提交前，脚本还会把 CRLF 换行转成 LF（原文件备份为 `.crlf.bak`），并补上输入末尾缺少的空行。Gaussian 的输出写在输入旁边的 `<stem>.log`。

**局限**：只支持单节点共享内存并行（不处理 Linda 多节点）。不同集群的 module 名、scratch 位置和分区策略各不相同，第一次使用时先用 `-n` 检查生成的脚本。相关说明见[第17章](../tutorials/17-hpc-jobs.md)和 [Common 作业调度章节](../../Common/tutorials/02-job-scheduler.md)。

## g16_error_doctor.py：报错诊断

**用途**：读取 `.log/.out`，判断作业是正常结束、报错结束还是被中断，并给出中文的原因解释和处理建议。

- 优先从 `Error termination via Lnk1e in .../lNNN.exe` 识别出错的 Link（兼容 Linux 与 Windows 路径）；没有这一行时，改用 `Error termination request processed by link NNN`。
- 只在**出错的那个作业步**内匹配约 25 类报错原文，因此前面 --Link1-- 步中无害的警告不会被误判。例如：L502 `Convergence failure`、L508 QC 失败、L9999 `Number of steps exceeded`、L123 corrector、`QPErr`、`Illegal IType`、`End of file in ZSymb`、多重度与电子数不匹配、`Atomic number out of range`、`R6DS8`、基组中心/ECP 读到文件尾、原子间距过小、FormBX/内坐标失败、`Wrong number of Negative eigenvalues`、`Inaccurate quadrature`、冻核混合、`galloc`、`Erroneous write/read`、`NtrErr`、segmentation violation、`Link 1 not found`。`Linear angle in Bend`、`Bend failed` 这类在正常作业中也会出现的信息，只有出现在报错前最后 80 行内时才会被采用。
- 没有匹配到已知报错时，说明出错 Link 的功能。
- 文件中没有任何结束行时，提示作业可能被调度器杀掉，建议用 `sacct` 查看；如果文件在最近 5 分钟内还有改动，则提示作业可能仍在运行。
- 对正常结束的作业也做检查：是否有 `Convergence criterion not met`；Opt 作业是否缺少 `Stationary point found`；虚频个数（极小点应为 0，TS 应为 1）；`<S**2>` 自旋污染；`Stable` 检查是否发现不稳定。

**依赖**：Python >= 3.9 标准库；`--draft` 需要同目录下的 [log2xyz.py](log2xyz.py)。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `logs...` | — | 输出文件或目录 |
| `-r/--recursive` | — | 递归搜索目录 |
| `-c/--context N` | 12 | 显示 Error termination 上方的行数 |
| `-q/--quiet` | — | 每个文件只输出一行状态 |
| `--json FILE` | — | 另外把诊断写成 JSON（`-` 表示输出到 stdout） |
| `--draft [GJF]` | `<stem>_restart.gjf` | 起草续算输入（只能用于单个 log） |
| `--force` | — | 允许 `--draft` 覆盖已有文件 |
| `--list-errors` | — | 列出所有可识别的报错并退出 |

**用法示例**：

```bash
python g16_error_doctor.py job.log                    # diagnosis in Chinese
python g16_error_doctor.py *.log -q                   # one status line per file
python g16_error_doctor.py runs/ -r --json report.json
python g16_error_doctor.py job.log --context 30
python g16_error_doctor.py job.log --draft            # writes job_restart.gjf for review
python g16_error_doctor.py job.log --draft fix.gjf --force
python g16_error_doctor.py --list-errors
```

**`--draft` 如何修改 route**：SCF 不收敛时加 `SCF=XQC`；QC 失败时去掉 QC、加 `VShift=500`；积分精度问题加 `Int=SuperFine`；步数用完时 `MaxCycles` 取 max(200, 2 × 原值)；TS 问题加 `CalcFC,NoEigenTest`。同时去掉 `Geom=Check` 和 `Opt Restart`，log 中已有收敛的 SCF 时加 `Guess=Read`。结构取自 log 中最后一个结构（调用 log2xyz.py）。输入错误、IRC 错误和内存错误不生成草稿，需要人工处理。

**退出码**：只要有一个 log 是报错结束或未完成，返回 1，可在 shell 脚本中直接判断。

**局限**：报错原文按 G16 C.01/C.02 编写，其他版本的措辞可能略有不同，请用 `--context` 自己看报错前的内容。Gen/GenECP 基组块不在 log 中，草稿里要手工补上。草稿只是起点，提交前必须检查。报错的详细解释见[第18章](../tutorials/18-errors-capstone.md)，SCF 问题见[第4章](../tutorials/04-sp-scf.md)，TS/IRC 问题见[第7章](../tutorials/07-ts-irc.md)。

## g16_log_summary.py：批量汇总输出

**用途**：每个文件只读一遍（几 GB 的 log 也能处理），汇总以下内容：结束状态与出错的 Link、route、最终能量及其类型、优化是否收敛、虚频数、ZPE/H/G 热校正、`<S**2>`（湮灭前/后）、基函数数、CPU 与墙钟时间。多步作业的能量取最后一个打印了能量的作业步，热化学取最后一个频率步。

**依赖**：Python >= 3.9 标准库；`--cclib` 需要 `pip install cclib`。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `paths...` | — | 文件、通配符或目录 |
| `-r/--recursive` | — | 递归搜索子目录 |
| `--ext` | `.log,.out` | 在目录中搜索的扩展名 |
| `--csv FILE`、`--json FILE` | — | 把全部字段导出为 CSV/JSON |
| `--columns` | `file,status,error_link,e_final,e_final_type,opt_status,n_imag,g_sum,s2_before,wall_time_s` | 屏幕表格中显示的列（逗号分隔） |
| `--sort` | `name` | `name`、`energy`、`G`、`status` |
| `--relative` | — | 增加 `dE_kcal`、`dG_kcal` 两列（相对最低值，kcal/mol） |
| `--cclib` | — | 用 cclib 交叉验证 SCF 能量和虚频数 |
| `-q/--quiet` | — | 不打印表格（配合 `--csv/--json` 使用） |

**用法示例**：

```bash
python g16_log_summary.py water.log
python g16_log_summary.py *.log --csv summary.csv --json summary.json
python g16_log_summary.py calc_dir/ -r --sort G --relative
python g16_log_summary.py ts.log --cclib
```

**说明**：`e_final` 取最具体的能量，优先级为 Counterpoise 校正 > ONIOM 外推 > 复合方法（0 K）> CCSD(T) > MP2/双杂化 > TD-DFT 总能量 > SCF Done，类型记录在 `e_final_type` 中。能量单位为 Hartree。

**局限**：只汇总，不判断方法是否合理。比较不同体系的能量前，先确认它们的计算水平、溶剂和标准态一致。

## log2xyz.py：提取结构与生成续算输入

**用途**：从 Gaussian 输出提取结构并写成 xyz（单帧或多帧），能量以 `key=value` 形式写在注释行（ASE extxyz 可以读取）。也可以为中断或失败的作业生成续算 `.gjf`。

**依赖**：Python >= 3.9 标准库。

**参数**：

| 参数 | 说明 |
|---|---|
| `--last`（默认）/`--lowest`/`--all`/`--converged`（别名 `--scan`）/`--irc`/`--index N` | 选择哪些帧，几个选项互斥：最后一个、能量最低、全部、每个 `Stationary point found` 的收敛结构（柔性扫描的每个点）、按反应坐标排序的 IRC 路径（reverse 端 → TS → forward 端）、第 N 帧（从 1 开始，负数表示从末尾倒数） |
| `--orient auto\|input\|standard` | 写出哪种取向；`auto` 时有 Input orientation 就用它，否则用 Standard |
| `-o/--output FILE` | 输出文件（`-` 表示 stdout），默认 `<log stem>[_mode].xyz` |
| `--list` | 只打印所有帧的列表 |
| `--gjf FILE` | 用选中的那一帧另外生成续算 gjf |
| `--route TEXT`、`--maxcycles N`、`--guess-read`、`--add-route TEXT` | 续算 route：替换整行、设 Opt MaxCycles、加 `Guess=Read` + `%oldchk`、追加关键词（如 `'SCF=XQC'`） |

**用法示例**：

```bash
python log2xyz.py opt.log                      # last geometry -> opt.xyz
python log2xyz.py opt.log --lowest             # -> opt_lowest.xyz
python log2xyz.py opt.log --all -o traj.xyz    # whole optimization trajectory
python log2xyz.py scan.log --converged         # every relaxed-scan point
python log2xyz.py irc.log --irc                # reverse end -> TS -> forward end
python log2xyz.py opt.log --list
python log2xyz.py opt.log --index 12 -o -
python log2xyz.py failed.log --gjf restart.gjf --maxcycles 200 --guess-read
```

**续算 gjf 的内容**：沿用 log 中第一个用户 route，但去掉 `Geom=Check/AllCheck/Connectivity`，因为续算文件写的是显式坐标。另外沿用回显的 Link 0 行、电荷与多重度、"Symbolic Z-matrix" 回显中的原子标签/冻结标记/ONIOM 层列，以及 ModRedundant 行。

**局限**：Gen/GenECP 基组块不在 log 中，必须手工补上。IRC 每个点写出的是 "Point Number: n" 之前打印的最后一个结构；用默认的 HPC 积分器时，这是预测步的结构，与汇总表中校正后的点非常接近，但并不完全相同。

## opt_monitor.py：优化过程监控

**用途**：每个 Berny 步（即每个 `Converged?` 表）记录一次能量和四项收敛判据（Maximum/RMS Force、Maximum/RMS Displacement），打印表格，诊断属于"振荡""缓慢但稳定下降"还是"平坦势能面"，并给出建议。遇到 L9999（步数用完）时，可以用它决定下一步怎么做。正在运行的作业也能分析。

**依赖**：Python >= 3.9；作图需 matplotlib（`--no-plot` 时不需要）。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-o/--output` | `<log stem>_opt.png` | 图片文件（单个 log） |
| `--no-plot`、`--show`、`--dpi` | —、—、150 | 不作图 / 弹出窗口 / 分辨率 |
| `--csv FILE` | — | 把逐步表格写成 CSV（单个 log） |
| `--scan-point K` | — | 只看柔性扫描的第 K 个点 |
| `--job N` | — | 只看第 N 个作业步（Link1）中的优化步 |
| `--last N` | — | 只看最后 N 步 |
| `-q/--quiet` | — | 不打印表格，只输出诊断 |

**用法示例**：

```bash
python opt_monitor.py opt.log                    # table + diagnosis + opt_opt.png
python opt_monitor.py opt.log -o conv.png --dpi 200
python opt_monitor.py opt.log --no-plot --csv steps.csv
python opt_monitor.py scan.log --scan-point 7
python opt_monitor.py ts.log --last 30 --show
```

**说明**：只分析 route 中含 Opt 的作业步。Opt Freq 中频率步打印的那张附加收敛表（用精确 Hessian 重新检验）单独报告，不算作优化步。Gaussian 判定收敛有两种情况：四项全部为 YES；或者力已比阈值小 100 倍，此时即使位移未达标也算收敛（"力可以忽略"）。图中第一栏是相对最低点的能量（kcal/mol），另外四栏为对数坐标，虚线表示阈值；实心点表示该判据已满足，空心点表示未满足。

**局限**：诊断是基于规则的提示。具体处理方法（CalcFC、减小步长、改用 Cartesian 坐标等）见[第5章](../tutorials/05-opt.md)和[第18章](../tutorials/18-errors-capstone.md)。

## scan_irc_profile.py：扫描与 IRC 能量曲线

**用途**：

- **柔性扫描**（`Opt=ModRedundant`，一维或二维）：取每个扫描点收敛结构的能量；扫描坐标值根据 ModRedundant 中 `S` 行的原子编号，由结构直接计算（键长、键角或二面角）。
- **刚性扫描**（`Scan` 关键词）：读取 "Summary of the potential surface scan" 汇总表；没有汇总表时逐点读取 SCF Done。
- **IRC**：读取 "Summary of reaction path following"；作业未完成时逐点读取。正反两个方向（同一个 log，或 forward/reverse 两个 log）自动合并成一条曲线。

结果以相对能量作图并写入 CSV，同时报告最高点（可作 TS 初猜）和各个局部极大值。二维扫描画等高线图。

**依赖**：Python >= 3.9 标准库；作图需 matplotlib（`--no-plot` 时不需要）。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--mode auto\|scan\|irc` | `auto` | 作业类型（默认自动识别） |
| `--coord SPEC` | — | 手动指定坐标（可重复），如 `'B 1 2'`、`'A 1 2 3'`、`'D 1 2 3 4'`；适用于 GIC 或 Z-matrix 变量扫描 |
| `--no-unwrap` | — | 二面角保留在 (-180, 180] 区间，不展开成连续曲线 |
| `--ref` | `min` | 能量零点：`min`、`first`、`ts`（IRC），或一个 Hartree 数值 |
| `--unit` | `kcal` | `kcal`、`kj`、`ev` |
| `--flip`、`--separate` | — | IRC：反应坐标取反号；多个 log 不合并 |
| `-o/--output` | `<first log stem>_profile.png` | 图片文件 |
| `--csv FILE`、`--write-max XYZ` | — | 导出曲线；把扫描最高点的结构写成 xyz |
| `--no-plot`、`--show`、`--dpi` | —、—、150 | 作图控制 |

**用法示例**：

```bash
python scan_irc_profile.py butane_scan.log                        # 1D relaxed scan
python scan_irc_profile.py sn2_scan.log --write-max ts_guess.xyz  # TS guess from the highest point
python scan_irc_profile.py scan2d.log --csv scan2d.csv            # 2D scan -> contour plot
python scan_irc_profile.py irc.log                                # both IRC directions in one log
python scan_irc_profile.py irc_fwd.log irc_rev.log -o irc.png     # two logs merged
python scan_irc_profile.py irc.log --ref ts --unit kj
python scan_irc_profile.py gic_scan.log --coord "B 1 5" --coord "B 1 6"
```

**局限**：未收敛的柔性扫描点（没有 `Stationary point found`）会被跳过并单独列出，因为它们的能量不可信。IRC 反应坐标的单位是 Gaussian 打印的 Bohr amu^1/2，reverse 方向为负值。IRC 两端只是近似的极小点，报告能量前要先单独优化。没有 ModRedundant 行、也没有给 `--coord` 时，横坐标用点的序号。

## thermo_qrrho.py：quasi-RRHO 热化学（教学版）

**用途**：读取最后一个频率块和热化学部分（温度、压力、分子量、转动常数、转动对称数、多重度、频率），从头重算理想气体 RRHO 热化学。设置与 log 相同时，会和 Gaussian 自己输出的数值对比，作为自检。在此基础上可以选择：

- `--qs grimme`：Grimme qRRHO 熵，$w = 1/(1+(\nu_0/\nu)^\alpha)$，$\nu_0 = 100\ \mathrm{cm^{-1}}$，$\alpha = 4$（Grimme, *Chem. Eur. J.* **2012**, 18, 9955, [doi:10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)）；
- `--qs truhlar`：把低于截断值的频率提高到截断值，再计算振动熵（Ribeiro et al., *J. Phys. Chem. B* **2011**, 115, 14556, [doi:10.1021/jp205508z](https://doi.org/10.1021/jp205508z)）；
- `--qh`：Head-Gordon 型焓插值（Li et al., *J. Phys. Chem. C* **2015**, 119, 1840, [doi:10.1021/jp509921r](https://doi.org/10.1021/jp509921r)），本脚本对每个模式的整个振动能（ZPE + 热）做插值，与其他程序的实现细节可能不同；
- `--conc C`：标准态校正，$\Delta G = RT\ln(CRT/P)$（1 M、298.15 K 时为 1.89 kcal/mol）；
- `--sp`：用高精度单点能替换电子能量（可以给 Hartree 数值或 log 文件）；
- `--boltz`：把所有 log 当作同一物种的构象，做 Boltzmann 加权。

**依赖**：`pip install numpy`。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--qs grimme\|truhlar\|none` | `grimme` | 低频熵的处理方法 |
| `--qh` | — | 同时做焓的准谐振插值 |
| `--cutoff`、`--alpha` | 100、4 | 截断频率 $\nu_0$（cm⁻¹）与阻尼指数 |
| `-t/--temp`、`-p/--pressure` | 取 log 中的值 | 温度（K）、压力（atm） |
| `-c/--conc` | — | 标准态浓度（mol/L），如 `1.0` |
| `-f/--freq-scale`、`--zpe-scale` | 1.0、与 `--freq-scale` 相同 | 频率与 ZPE 缩放因子 |
| `--mult` | 取 log 中的值 | 覆盖计算 $S_\mathrm{elec} = R\ln(\mathrm{mult})$ 所用的多重度 |
| `--sp ...`、`--sp-key` | —、`auto` | 每个 log 对应一个单点能；`--sp-key` 指定从 archive 中取哪个能量（`HF`、`MP2`、`CCSD(T)` 等） |
| `--boltz` | — | 构象 Boltzmann 加权 |
| `--csv`、`--json` | — | 导出结果 |
| `-v/--verbose` | — | 打印熵的各项分量和逐模式表格 |

**用法示例**：

```bash
python thermo_qrrho.py h2o_opt_freq.log
python thermo_qrrho.py ts.log --qs grimme --qh --conc 1.0
python thermo_qrrho.py react.log --sp react_dsdpbep86.log --freq-scale 0.986
python thermo_qrrho.py conf1.log conf2.log conf3.log --qs grimme --boltz --csv confs.csv
python thermo_qrrho.py a.log b.log --sp -382.31234 -382.30998 --temp 333.15
python thermo_qrrho.py mol.log --qs grimme -v
```

> [!CAUTION]
> 本脚本用于教学，**不能替代** [GoodVibes](https://github.com/patonlab/GoodVibes) 或 [Shermo](http://sobereva.com/soft/shermo/)。科研计算请使用这些工具并引用它们。不同 qRRHO 实现的截断函数、$B_\mathrm{av}$ 和焓处理方式略有差别，同一项研究中不要混用。原理见[第6章](../tutorials/06-freq-thermo.md)。

## spectrum_ir_raman.py：IR / Raman 光谱

**用途**：读取最后一个普通格式（非 HPModes）的 "Harmonic frequencies" 块，包括 `Frequencies`、`IR Inten`、`Raman Activ`、`Depolar`；频率依赖（预共振）Raman 作业还读取 `RamAct Fr= n`。虚频会报告后丢弃，`*****` 溢出值记为 NaN。频率缩放后做面积归一化的洛伦兹或高斯展宽。IR 纵轴为十进制摩尔吸光系数 $\varepsilon$（L mol⁻¹ cm⁻¹），由 Gaussian 输出的积分强度（km/mol）换算得到。加 `--raman-intensity` 时，按激光波长和温度把 Raman 活性（Å⁴/amu）换算为相对强度：

$$I_i \propto \frac{(\tilde\nu_0-\tilde\nu_i)^4 S_i}{\tilde\nu_i\left[1-\exp(-hc\tilde\nu_i/k_BT)\right]}$$

**依赖**：`pip install numpy matplotlib`。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--labels` | — | 图例标签（每个 log 一个） |
| `-s/--scale` | 1.0 | 频率缩放因子 |
| `-w/--fwhm`、`--shape` | 10 cm⁻¹、`lorentzian` | 半高全宽、线型（`lorentzian`/`gaussian`） |
| `--mode auto\|ir\|raman\|both` | `auto` | 画哪些谱（`auto`：log 中有 Raman 数据时一起画） |
| `--raman-intensity`、`--laser`、`-t/--temp` | —、532 nm、298.15 K | Raman 活性换算为强度 |
| `--raman-fr N` | — | 使用 `RamAct Fr= N` 块（`CPHF=RdFreq`） |
| `--xmin`、`--xmax`、`--step` | 0、自动、1 | 波数范围与网格间距（cm⁻¹） |
| `--normalize`、`--sticks`、`--ir-ascending` | — | 归一化到最大值 1；叠加多个体系时也画棒状谱；IR 横轴改为从小到大（默认从大到小） |
| `-o/--out`、`--dpi`、`--show` | 第一个 log 名、300 | 输出前缀、分辨率、弹出窗口 |

**输出**：`<prefix>_modes.csv`（各模式数据）、`<prefix>_spectrum.csv`（展宽后的曲线）、`<prefix>.png`。

**用法示例**：

```bash
python spectrum_ir_raman.py h2o_freq.log
python spectrum_ir_raman.py mol.log --scale 0.967 --fwhm 8 --xmin 400 --xmax 4000
python spectrum_ir_raman.py mol.log --raman-intensity --laser 785 --temp 298.15
python spectrum_ir_raman.py h2o.log d2o.log --labels H2O D2O --mode ir -o isotopes
python spectrum_ir_raman.py prer.log --raman-fr 2 --raman-intensity --laser 532
```

**局限**：只处理谐振近似（不读 `Freq=Anharmonic` 的非谐结果）。缩放因子要与所用方法和基组匹配，可查 [NIST CCCBDB](https://cccbdb.nist.gov/)。不读 HPModes 格式的输出。原理见[第6章](../tutorials/06-freq-thermo.md)。

## spectrum_uvvis_ecd.py：UV-Vis / ECD 光谱

**用途**：读取最后一个 "Excitation energies and oscillator strengths" 块（Opt TD(Root=n)、--Link1-- 和 TD=(50-50) 作业都能处理），以及它前面的旋转强度表（默认用 R(velocity)，`--rot length` 改用 R(length)）。按 Gaussian 白皮书 "Creating UV/Visible Plots" 的公式，在能量（eV）空间做高斯展宽：默认 $\sigma = 0.4$ eV（GaussView 的默认值，$\sigma$ 为 1/e 处的半宽，FWHM = 1.665σ）。支持能量平移、对映体翻转、多构象 Boltzmann 加和，以及叠加实验谱。

**依赖**：`pip install numpy matplotlib`。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--sigma` / `--fwhm` | 0.4 eV / — | 展宽（二选一，单位 eV） |
| `--sigma-ecd` | 与 UV 相同 | ECD 单独使用的展宽 |
| `--shift` | 0 | 所有激发能平移的量（eV），如 `-0.2` |
| `--flip` | — | R 取反号（对映体） |
| `--rot velocity\|length` | `velocity` | 旋转强度的规范 |
| `--weights ...` / `--energies ...` | — | 构象权重（自动归一化），或构象能量（数值或 log；log 中有频率时取 G），二选一 |
| `--energy-unit`、`-t/--temp` | `hartree`、298.15 K | 数值形式 `--energies` 的单位（`hartree`/`kcal`/`kj`/`ev`），以及 Boltzmann 温度 |
| `--xunit`、`--xmin`、`--xmax`、`--npts` | `nm`、自动、自动、2000 | 横轴单位（`nm`/`eV`/`cm-1`）、范围与网格点数 |
| `--no-ecd` | — | 不输出 ECD |
| `--exp-uv`、`--exp-ecd`、`--exp-unit` | —、—、`nm` | 实验谱文件（两列 x y，按最大值缩放后叠加） |
| `-o/--out`、`--dpi`、`--show` | 第一个 log 名、300 | 输出控制 |

**输出**：`<prefix>_states.csv`（激发态列表）、`<prefix>_spectrum.csv`、`<prefix>.png`。

**用法示例**：

```bash
python spectrum_uvvis_ecd.py coumarin_td.log
python spectrum_uvvis_ecd.py td.log --sigma 0.3 --shift -0.2 --xmin 200 --xmax 500
python spectrum_uvvis_ecd.py conf1_td.log conf2_td.log --energies c1_freq.log c2_freq.log \
       --flip --exp-ecd exp_ecd.txt -o methyloxirane
python spectrum_uvvis_ecd.py a.log b.log --weights 0.7 0.3 --xunit eV --no-ecd
```

**局限**：展宽只能近似表示振动精细结构，不能代替振动分辨光谱（`FC`/`FCHT`）。`--shift` 是经验修正，使用时要在文中说明。多构象计算的原子编号与计算水平必须一致。原理与 ECD 判断绝对构型的流程见[第10章](../tutorials/10-tddft.md)。

## nmr_shift.py：NMR 化学位移换算

**用途**：读取 NMR 计算（`NMR=GIAO/CSGT/IGAIM`）中的 `Isotropic =` 屏蔽常数，换算为化学位移：

- **参比法**：$\delta = \sigma_\mathrm{ref} - \sigma$，$\sigma_\mathrm{ref}$ 来自**同一水平**计算的 TMS（`--ref-log`，脚本自动对同种原子取平均），也可以用 `--sigma-ref` 手工给出；
- **线性回归缩放**：$\delta = (\mathrm{intercept} - \sigma)/(-\mathrm{slope})$（`--scale`）。参数必须与你的计算协议完全一致，例如 [CHESHIRE](http://cheshirenmr.info/) 数据库中对应的条目；
- 等价核平均（甲基 H 等），多构象 Boltzmann 平均；
- 与实验值比较：给出每个元素的 MAE、RMSE 和最大偏差；`--fit` 对 $\sigma$ 与 $\delta_\mathrm{exp}$ 做线性拟合。

**依赖**：Python >= 3.9 标准库（不需要 numpy）。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `logs...` | — | NMR 输出；多个 log 视为同一分子的不同构象（原子编号必须一致） |
| `--ref-log LOG` | — | 同一水平计算的参比分子（如 TMS）输出 |
| `--sigma-ref EL=SIGMA` | — | 手工给出参比屏蔽常数，可重复，如 `H=31.80` |
| `--scale EL=SLOPE,INTERCEPT` | — | 线性回归缩放参数，可重复 |
| `--assign FILE` | — | 归属文件，每行为"标签 原子编号 [实验位移]" |
| `--equiv ATOMS` | — | 不用归属文件时对指定原子取平均，如 `12-14`，可重复 |
| `--elements EL...` | 全部 | 只列出这些元素 |
| `--weight scf\|gibbs\|file\|equal` | `scf` | Boltzmann 权重的来源：最后的 SCF 能量、同一 log 中的自由能、`--energies` 文件，或等权重 |
| `--energies FILE`、`--temp` | —、298.15 K | 每行为"log 名 能量（Hartree）"（配合 `--weight file`）；温度 |
| `--fit`、`--csv FILE`、`--show-sigma` | — | 线性拟合；导出 CSV；打印每个构象的屏蔽常数 |

**归属文件示例**（`#` 开头为注释，原子编号中不要有空格）：

```text
# label   atoms (Gaussian numbering)   [experimental shift / ppm]
C1        7                            57.7
H8        12-14                        0.84
H3exo     26                           2.35
```

**用法示例**：

```bash
python nmr_shift.py camphor_nmr.log --ref-log tms_nmr.log
python nmr_shift.py camphor_nmr.log --sigma-ref H=31.80 --sigma-ref C=184.00 --elements H C
python nmr_shift.py conf*.log --ref-log tms.log --assign assign.txt --weight gibbs --temp 298.15
python nmr_shift.py conf*.log --ref-log tms.log --weight file --energies conf_G.csv
python nmr_shift.py mol.log --scale H=-1.0381,31.7427 --equiv 12-14 --equiv 15-17
python nmr_shift.py mol.log --ref-log tms.log --assign assign.txt --fit --csv shifts.csv
```

**局限**：每个 log 只取最后一个屏蔽常数块。参比分子必须在完全相同的水平下计算（几何优化水平、NMR 泛函/基组、溶剂模型都相同）；脚本会比较两者的 route，不一致时给出警告。上面示例中的数值只是演示格式的占位值，不是推荐参数。可交换质子（OH/NH/COOH）通常不参与误差统计。详见[第11章](../tutorials/11-nmr-props.md)。

## batch_cube_multiwfn.sh：批量 cube 与 Multiwfn

**用途**：对目录中的每个 `*.chk`（或 `*.fchk`），依次执行：

1. `formchk` 把 chk 转为 fchk（fchk 比 chk 新时跳过，加 `-f` 强制重做）；
2. 从 fchk 读取 HOMO/LUMO 能量（eV），每个文件写一行到 `batch_cube_summary.tsv`；
3. 用 `cubegen` 生成所需的 cube。开壳层波函数的轨道分别输出 alpha 和 beta（`AMO=n`/`BMO=n`）；ESP 在密度 cube 的同一网格上计算；
4. 对每个 `--mwfn FILE[:TAG]` 命令流，在独立目录 `<stem>_<TAG>/` 中执行 `Multiwfn file.fchk -silent`，避免 Multiwfn 固定名称的输出文件互相覆盖。

**依赖**：bash；Gaussian 16 的 `formchk`、`cubegen`（先 source `g16.profile`）；Multiwfn（可选，只在使用 `--mwfn` 时需要，`settings.ini` 从 `$Multiwfnpath` 读取）。

**参数**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-k/--kinds LIST` | `homo,lumo,density` | 逗号分隔：`homo`、`lumo`、`homo-N`、`lumo+N`、`mo=N`、`density`、`spin`（仅开壳层）、`esp`；也可以写任意 cubegen kind，如 `Density=Current` |
| `-n/--npts N` | -3 | cubegen 的 npts：0 = 80³；-2/-3/-4 = 每 Bohr 3/6/12 个点 |
| `-p/--nproc N` | `$SLURM_CPUS_PER_TASK` 或 4 | 计算 ESP cube 时用的核数 |
| `-M/--mem SIZE` | — | formchk/cubegen 的内存（以 `-m=SIZE` 传入，需 G16 C.01 及以后版本） |
| `--no-cube` | — | 只做 formchk（以及 Multiwfn） |
| `-w/--mwfn FILE[:TAG]`、`--multiwfn PATH` | —、`$MULTIWFN` 或 PATH 中的 `Multiwfn` | Multiwfn 命令流（可重复）；可执行文件路径 |
| `-o/--outdir DIR`、`-r/--recursive`、`-f/--force`、`--dry-run` | — | 输出目录、递归搜索、强制重做、只打印命令 |

**用法示例**：

```bash
bash batch_cube_multiwfn.sh                                  # *.chk in ., homo,lumo,density
bash batch_cube_multiwfn.sh -k homo-1,homo,lumo,lumo+1 -n -4 calc/
bash batch_cube_multiwfn.sh -k density,esp -p 16 phenol.chk  # ESP on the density grid
bash batch_cube_multiwfn.sh -k spin radical.fchk
bash batch_cube_multiwfn.sh --no-cube --mwfn esp_cmds.txt:esp --mwfn adch_cmds.txt:adch runs/
bash batch_cube_multiwfn.sh -r -o cubes/ project/
```

**输出**：`<stem>_homo.cub`（开壳层为 `_alpha/_beta`）、`<stem>_dens.cub`、`<stem>_esp.cub`、`<stem>_spin.cub`、`batch_cube_summary.tsv`，以及 Multiwfn 结果目录 `<stem>_<TAG>/`（内含 `<TAG>.out`）。

**Multiwfn 命令流**：每行一个回答，和交互运行时输入的内容完全一样（空行表示回车）。先手动跑一遍确认流程，Multiwfn 每次更新后都要重新核对。示例见 [multiwfn_esp_surface_cmds.txt](../examples/16-postprocess-multiwfn/multiwfn_esp_surface_cmds.txt) 和 [multiwfn_igmh_dimer_cmds.txt](../examples/16-postprocess-multiwfn/multiwfn_igmh_dimer_cmds.txt)。`.out` 文件最后出现 "end-of-file during read" 只表示命令流已经读完。

**局限**：Multiwfn 的菜单编号会随版本变化，命令流不能保证通用。等值面推荐值：轨道 0.02–0.05 a.u.；映射 ESP 用的密度面取 0.001 a.u.。详见[第16章](../tutorials/16-postprocess-multiwfn.md)。

---

## 典型工作流

```bash
# 1) build inputs for a set of conformers and submit them
python gjf_builder.py confs.xyz --all-frames -t optfreq --disp GD3BJ --solvent Water --nproc 16 --mem 32GB
bash sub_g16.sh -a -t 48:00:00

# 2) watch / diagnose
python opt_monitor.py confs_f001_optfreq.log --no-plot
python g16_error_doctor.py . -q

# 3) collect and post-process
python g16_log_summary.py . --sort G --relative --csv summary.csv
python thermo_qrrho.py confs_f*_optfreq.log --qs grimme --conc 1.0 --boltz
```

`--all-frames` 时每帧的文件名为 `<stem>_fNNN_<task>.gjf`（如 `confs_f001_optfreq.gjf`），`sub_g16.sh` 把输出写成同名的 `.log`。完整的研究流程案例见[第18章](../tutorials/18-errors-capstone.md)，作业脚本与并行设置见[第17章](../tutorials/17-hpc-jobs.md)。

## 参考资料

- Gaussian 16 关键词手册：https://gaussian.com/keywords/
- Gaussian 16 运行环境与 Link 0：https://gaussian.com/running/ 、https://gaussian.com/link0/
- formchk / cubegen 说明：https://gaussian.com/formchk/ 、https://gaussian.com/cubegen/
- Gaussian 白皮书 "Creating UV/Visible Plots from the Results of Excited States Calculations"：https://gaussian.com/uvvisplot/
- Multiwfn：http://sobereva.com/multiwfn/
- GoodVibes：https://github.com/patonlab/GoodVibes
- cclib：https://github.com/cclib/cclib

> [!NOTE]
> 脚本按 Gaussian 16 C.01/C.02 的真实输出格式编写，但不同版本的输出可能略有差异。正式使用前，建议先用小体系核对结果；发现解析错误请到 GitHub Issues 反馈。

[返回目录](../README.md)
