# 第16章 检查点文件、cubegen 与 Multiwfn 后处理（Checkpoint Files, cubegen & Multiwfn Post-processing）

Gaussian 的计算结果大部分存放在二进制检查点文件（.chk）里，而不是 .log 里。本章介绍 Gaussian 16 自带的 formchk、unfchk、c8616、chkchk、freqchk、newzmat、cubegen、cubman 等工具，讲解如何生成轨道、电子密度和静电势（ESP）的 cube 格点文件，以及如何在 GaussView 6 中查看轨道、ESP 映射表面、振动、IRC 与扫描曲线。后半部分用 Multiwfn（2026 年起按发布日期命名版本）完成电荷布居、键级、轨道成分、光谱绘制、分子表面 ESP 分析、空穴-电子分析和 IGMH 弱相互作用可视化，并说明如何写 Multiwfn 命令流做批处理、如何用 VMD 渲染出版级图片、如何导出 wfx/molden，以及如何用 MOKIT 把 Gaussian 轨道交给 PySCF/ORCA 等程序。关键词：Gaussian formchk cubegen、fchk、cube file、Multiwfn ESP surface analysis、IGMH、hole-electron analysis、Multiwfn batch mode、VMD rendering、MOKIT。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 说清 .chk、.fchk、.log、.cub、.wfx 各自存什么，什么时候用哪个；
2. 熟练使用 `formchk`（含 C.01 起的 `-m` 选项）、`unfchk`、`c8616`、`chkchk`、`freqchk`、`newzmat` 处理检查点文件，能把 chk 安全地搬到另一台机器或另一个 Gaussian 版本；
3. 用 `cubegen` 生成 HOMO/LUMO、总电子密度、自旋密度和 ESP 的 cube 文件，知道 `npts` 各种取值的含义，并能让多个 cube 共用同一套格点；
4. 在 GaussView 6 中查看轨道、ESP 映射表面、振动动画、优化/IRC/扫描曲线和光谱；
5. 安装并配置 Multiwfn，知道主功能 0、7、8、9、11、12、18、20 分别做什么；
6. 写出 Multiwfn 命令流文件，用 `Multiwfn job.fchk < cmds.txt > out.txt` 非交互批量运行；
7. 用 VMD 渲染 ESP 表面和 IGMH 等值面；用 `Output=WFX` 导出波函数；用 MOKIT 把轨道转给其他程序。

| 你想做的事 | 本章对应小节 |
|---|---|
| 把作业结果拿到 GaussView/Multiwfn 里看 | 3.2 formchk、5 GaussView |
| 把 chk 从集群拷到自己电脑，或从旧版本升级 | 3.3 unfchk、3.4 c8616 |
| 忘了某个 chk 是什么作业算的 | 3.5 chkchk |
| 换温度、压力、缩放因子重算热化学，不重跑频率 | 3.6 freqchk |
| 从 chk 中取出结构（某一优化步、某个扫描点） | 3.7 newzmat |
| 画 HOMO/LUMO、电子密度、ESP、差分密度 | 4 cubegen、cubman |
| 原子电荷、键级、轨道成分 | 7.1–7.3 |
| 画 IR/Raman/UV-Vis/ECD/NMR 谱 | 7.4 |
| 分子表面静电势定量分析（极值点、面积、方差） | 7.5 |
| 激发态空穴-电子、NTO、电荷转移分析 | 7.6 |
| 弱相互作用可视化（NCI/IGMH） | 7.7 |
| 几十上百个 fchk 自动跑同一套分析 | 8 命令流批处理 |

**前置知识**：输入文件结构与 `%chk`/`%oldchk`（[第2章](02-input-anatomy.md)）、单点能（[第4章](04-sp-scf.md)）。具体性质的计算方法分散在各章：TD-DFT 见 [第10章](10-tddft.md)，布居与 ESP 电荷见 [第12章](12-population-charges.md)，弱相互作用见 [第14章](14-weak-bsse.md)。通用的可视化工具（VESTA、Multiwfn 跨程序用法、出版级作图规范）见 [Common 第5章](../../Common/tutorials/05-postprocess-visualization.md)。

本章配套示例：[../examples/16-postprocess-multiwfn/README.md](../examples/16-postprocess-multiwfn/README.md)（苯酚的 HOMO/LUMO、密度与 ESP cube，ESP 表面分析，TD-DFT UV-Vis 绘谱；水二聚体 IGMH；`Output=WFX` 示例）。

---

## 2. 原理速览

### 2.1 Gaussian 的几种结果文件

| 文件 | 格式 | 内容 | 可移植性 |
|---|---|---|---|
| `.log`/`.out` | 文本 | 人读的输出：能量、优化步、频率、激发态、布居分析 | 任意平台 |
| `.chk` | 二进制 | 几何、基组、分子轨道系数、密度矩阵、力常数（Hessian）、优化/扫描/IRC 历史等 | **依赖平台和版本**；C.01 起带版本号 |
| `.fchk`（Windows 下为 `.fch`） | 文本 | chk 的格式化版本，按“标签 + 数组”存放 | 任意平台，几乎所有后处理程序都能读 |
| `.cub`/`.cube` | 文本 | 三维格点上的标量场（轨道、密度、ESP 等） | 任意平台 |
| `.wfn`/`.wfx` | 文本 | AIM 风格波函数：原始高斯函数 + 自然轨道/分子轨道 | 任意平台 |
| `.rwf` | 二进制 | 运行时读写文件，作业结束一般自动删除 | 仅用于断点续算（见 [第17章](17-hpc-jobs.md)） |

> [!IMPORTANT]
> .log 里**没有**完整的波函数。想做轨道图、ESP、Multiwfn 分析，计算时就必须写 `%chk=xxx.chk`，并在作业结束后保留它。没写 `%chk` 的作业，检查点存放在 scratch 目录，结束后会被删掉，只能重算。

### 2.2 cube 文件里存的是什么

cube 文件把空间划分成 $N_1\times N_2\times N_3$ 个格点，每个格点上存一个（或几个）数值。常见的几种标量场：

- 分子轨道值：$\psi_i(\mathbf r)=\sum_\mu c_{\mu i}\,\phi_\mu(\mathbf r)$，有正有负，作图时用 ±等值面两种颜色表示相位；
- 电子密度：$\rho(\mathbf r)=\sum_i n_i\,|\psi_i(\mathbf r)|^2$，其中 $n_i$ 为占据数；
- 自旋密度：$\rho^\alpha(\mathbf r)-\rho^\beta(\mathbf r)$，开壳层体系才有意义；
- 静电势：

$$
V(\mathbf r)=\sum_A \frac{Z_A}{|\mathbf r-\mathbf R_A|}-\int \frac{\rho(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d\mathbf r'
$$

前一项是原子核的贡献，后一项是电子的贡献。每个格点都要做一次积分，所以 ESP cube 比密度 cube 慢得多，这也是 `cubegen` 的并行参数只对 ESP 有意义的原因。

格点越密，图越平滑，文件也越大：80×80×80 = 512 000 个点，按每个数 13 个字符估算约 6–7 MB；150³ 就超过 40 MB。所以要按用途选格点（见 4.3 节）。

### 2.3 分子表面 ESP 与弱相互作用指标

分子表面 ESP 分析一般采用 Bader 的定义：取电子密度 $\rho=0.001$ a.u. 的等值面作为范德华表面，把 ESP 映射到这个面上，再统计极大/极小点、正负区域面积、平均值和方差。正 ESP 区域（如 O–H 上的 H、卤键的 σ-hole）倾向与负 ESP 区域（孤对电子、π 电子）接触，因此常用来解释氢键、卤键、反应位点。

IGMH（independent gradient model based on Hirshfeld partition）用 Hirshfeld 划分得到各碎片的密度梯度，定义 $\delta g$ 描述密度“重叠”的程度。$\delta g^{\mathrm{inter}}$ 只反映碎片之间的相互作用。作图时把 $\mathrm{sign}(\lambda_2)\rho$ 映射到 $\delta g^{\mathrm{inter}}$ 等值面上：蓝色代表强吸引（氢键），绿色代表范德华作用，红色代表位阻排斥。原理详见 Multiwfn 手册 3.23 节和 IGMH 原始论文（参考资料）。

---

## 3. 检查点文件工具链

### 3.1 先把 chk 写对

```text
%oldchk=phenol_opt.chk
%chk=phenol.chk
%mem=8GB
%nprocshared=8
#p B3LYP/def2TZVP Geom=Check Guess=Read
```

- `%chk` 指定本步作业读写的检查点；`%oldchk`（可选）在本步开始时把旧 chk 的内容**复制**进新 chk，旧文件不被改动。读取旧结构做性质计算时推荐这种写法（[第2章](02-input-anatomy.md)）。
- 多步 `--Link1--` 作业的每一步都要写 `%chk`，否则那一步的 chk 不会保留。
- 如果后续要做 Multiwfn 分析，最好让**最后一步就是你要分析的那个波函数**（例如优化完成后再在大基组上做单点），这样 fchk 里的轨道就是你想要的。

`g16` 命令还有一个 `-fchk` 选项：在生成二进制 chk 的同时生成一个 fchk，chk 每次更新时 fchk 也同步更新，省去事后手动 `formchk`：

```bash
g16 -fchk=phenol.fchk phenol_opt_sp.gjf
```

### 3.2 formchk：二进制 chk → 文本 fchk

```bash
formchk phenol.chk                 # writes phenol.fchk (base name + .fchk)
formchk phenol.chk phenol_tz.fchk  # explicit output name
formchk -3 phenol.chk              # version-3 fchk (G09/G16 format)
formchk -m=1gb big.chk             # Rev. C.01+: memory for the utility; -m must come first
```

| 选项 | 作用 | 什么时候用 |
|---|---|---|
| `-m=1gb` | 指定工具可用内存（C.01 起所有工具都支持） | 大体系 formchk 报内存不足时；`-m` 必须写在所有其他选项和文件名之前 |
| `-3` | 生成 version 3 的 fchk，包含 G16 的全部特性 | 官方建议在不同类型的计算机之间搬运文件时加上 |
| `-2` | 生成 G03 时代的 version 2 格式 | 只在旧程序只认旧格式时 |
| `-c` | MM 原子类型以字符串而不是整数写出 | ONIOM(QM:MM) 结果要导入其他程序时 |
| `-r` | 从中途失败的 ONIOM chk 中恢复 real 体系（低层）的结构等信息 | ONIOM 作业中断后想在 GaussView 中修改结构 |
| `-matrix` / `-rawmatrix` | 输出矩阵元文件（`.dat`），供 GauOpen 接口读取 | 自己写程序读取 Gaussian 的积分/矩阵时 |

工具的默认内存也可以用环境变量 `GAUSS_MEMDEF` 设置（单位默认为 words，可加 `GB` 等后缀）；formchk 的默认选项可以写在 `Default.Route` 的 `-F-` 行或环境变量 `GAUSS_FDEF` 中（[第17章](17-hpc-jobs.md)）。

检查 fchk 是否完整：

```bash
head -n 3 phenol.fchk                          # title, job type/method/basis, ...
grep -E "^(Number of atoms|Charge|Multiplicity|Number of basis functions)" phenol.fchk
grep -c "Alpha MO coefficients" phenol.fchk    # should be 1
```

fchk 的前两行是标题和“作业类型 / 方法 / 基组”，之后每一项都是 `标签  类型  N=  数组长度` 的形式，例如 `Alpha Orbital Energies  R  N=  ...`，下面跟着数值。这种格式很容易用脚本解析（示例目录中的 `formchk_cubegen_homo_lumo.sh` 就是直接从 fchk 读取 HOMO/LUMO 能量）。

> [!TIP]
> 集群上算完后，建议**在集群上** formchk，再把 .fchk 拷回本地。fchk 是文本文件，压缩率很高（`gzip` 后通常只有原来的几分之一）。

### 3.3 unfchk：fchk → 本机二进制 chk

```bash
unfchk phenol.fchk phenol.chk     # rebuild a binary chk for THIS machine/version
unfchk -fixver old.fchk old.chk   # mark the resulting chk as current version (C.01+)
```

典型用途：

- **跨机器续算**：在 A 机器上 `formchk`，把 fchk 拷到 B 机器，再 `unfchk`，得到 B 机器可用的 chk，然后 `Geom=Check Guess=Read` 或 `Opt=Restart`；
- **Windows ↔ Linux**：G16W 生成的 chk 不能直接给 Linux 版 Gaussian 用，反之亦然，必须走 formchk/unfchk；
- **修复版本号**：C.01 起 chk 带版本号，读取没有版本号的旧 chk 会被拒绝。`unfchk -fixver` 会把新生成的 chk 标记为当前版本，即使输入的 fchk 里没有版本信息。

### 3.4 c8616：升级旧版本 chk

```bash
cp taxol.chk taxol.chk.bak    # c8616 converts IN PLACE, keep a backup
c8616 taxol.chk
```

`c8616` 把 Gaussian 86 到 Gaussian 09 的 chk 原地转换为 Gaussian 16 格式（文件名不变）。如果你用 G16 C.01/C.02 读取 G09 或 G16 A.03/B.01 时期的 chk 出现奇怪的错误，先试 `c8616`，或者在旧机器上 `formchk`、在新机器上 `unfchk -fixver`。

### 3.5 chkchk：这个 chk 是什么作业？

```bash
chkchk phenol           # .chk extension may be omitted
chkchk -p phenol > phenol_from_chk.txt
```

`chkchk` 打印 chk 中保存的标题、route、Gaussian 版本、原子数、电荷与多重度，并列出 chk 里有哪些数据（坐标、基组、MO 系数、力常数等）。B.01 起还会报告作业状态：正常完成、失败还是仍在运行。输出大致如下（示意输出，内容仅供参考）：

```text
 Checkpoint file phenol.chk:
 Title:  Phenol: B3LYP-D3(BJ)/def2-SVP opt+freq, then B3LYP/def2-TZVP SP
 Route:  #p B3LYP/def2TZVP Geom=Check Guess=Read
 ...
 MO coefficients present.
 This checkpoint is from a job which completed successfully.
```

`-p` 选项还会打印坐标、Gen 格式的基组和 `Guess=Cards` 格式的轨道，适合在丢失输入文件时“复原”一个作业。

### 3.6 freqchk：不重跑频率，换条件重算热化学

`freqchk` 从 chk（或 fchk）中读出 Hessian，按你给的温度、压力、频率缩放因子和同位素重新做振动分析与热化学：

```bash
freqchk phenol.chk -o phenol_thermo_373K.txt
```

运行后会依次询问（以 Gaussian 官方文档中的会话为准）：

| 提示 | 含义 | 说明 |
|---|---|---|
| `Write Hyperchem files?` | 是否写 HyperChem 文件 | 一般答 `n` |
| `Temperature (K)? [0=>298.15]` | 温度 | 输入 `0` 取默认值；直接回车不行 |
| `Pressure (Atm)? [0=>1 atm]` | 压力 | 同上 |
| `Scale factor for frequencies during thermochemistry?` | 热化学中使用的频率缩放因子 | 输入 `0` 用方括号中的默认值 |
| `Do you want the principal isotope masses? [Y]` | 是否使用最丰同位素 | 答 `n` 后可逐个原子输入同位素质量 |

官方文档给出的完整语法是 `freqchk checkpoint-file [options] [answers to prompts]`，即也可以把答案按顺序写在命令行上，适合在脚本中批量运行。由于提示的顺序和默认值可能随修订版本略有不同，建议先交互运行一次确认。常用选项：`-o file` 把结果写到文件；`-np=N` 多核做频率分析；`-sel`/`-read`/`-save` 选择、读取或保存简正模式。

> [!NOTE]
> 只想得到另一温度下的热化学，也可以在原作业中加 `Temperature=`、`Pressure=`，或者用 GoodVibes/Shermo 等程序直接读 .log（[第6章](06-freq-thermo.md)）。freqchk 的优势是能改同位素，而且不需要保留 .log。

### 3.7 newzmat：从 chk 中取结构、转格式

```bash
newzmat -ichk -oxyz  phenol.chk phenol.xyz      # last geometry -> xyz
newzmat -ichk -opdb  phenol.chk phenol.pdb      # -> PDB (with connectivity)
newzmat -ichk -ocart phenol.chk phenol_new.com  # -> Gaussian input, Cartesian coordinates
newzmat -ichk -oxyz -ot 2-5 scan.chk scan_2to5.xyz   # structures 2-5 stored on the chk
newzmat -ipdb -prompt ligand.pdb ligand.com      # PDB -> Gaussian input, prompts for route/title/charge
```

| 选项 | 含义 |
|---|---|
| `-i<fmt>` / `-o<fmt>` | 输入/输出格式：`chk`、`zmat`、`cart`、`xyz`、`pdb`、`mopac` 等 |
| `-ngeom N` | 取 chk 中第 N 个结构（等同 `Geom=(NGeom=N)`） |
| `-ot list` | 取多个结构，例如 `-ot 3,7` 或 `-ot 2-5` |
| `-step N` | 取优化第 N 步的结构；**对默认的冗余内坐标优化不可用**，此时改用 `Geom=(Check,Step=N)` |
| `-prompt` | 生成 Gaussian 输入时询问 route、标题、电荷与多重度（默认生成 HF/6-31G(d) 单点） |

从 .log 中提取结构（每一步、能量最低的一步、IRC 路径）更方便的是仓库脚本 [`log2xyz.py`](../scripts/log2xyz.py)，用 `python ../scripts/log2xyz.py -h` 查看用法。

### 3.8 工具速查

| 工具 | 一句话功能 | G16W 菜单中有 |
|---|---|---|
| `formchk` | chk → fchk | 是 |
| `unfchk` | fchk → chk | 是 |
| `c8616` | 旧 chk 升级到 G16 格式 | 否 |
| `chkchk` | 显示 chk 的 route、标题、内容与作业状态 | 是 |
| `cubegen` | 从 fchk 生成 cube | 是 |
| `cubman` | cube 加减、缩放、求性质 | 是 |
| `freqchk` | 从 chk 重做频率分析与热化学 | 是 |
| `freqmem` | 估计频率计算所需内存 | 否 |
| `newzmat` | 几何格式转换、从 chk 取结构 | 是 |
| `testrt` | 检查 route 语法 | 是 |

---

## 4. cubegen：生成格点文件

### 4.1 语法

```text
cubegen  nprocs  kind  fchkfile  cubefile  npts  format  [cubefile2]
```

| 参数 | 含义 | 常用值 |
|---|---|---|
| `nprocs` | 计算 ESP 时用的共享内存核数；0 等同于 1。**只要后面还有参数，这一项就必须写** | `0`（轨道、密度）；`8`（ESP） |
| `kind` | 要生成的量，见 4.2 | `MO=Homo`、`Density=SCF`、`Potential=SCF` |
| `fchkfile` | 输入的 fchk | `phenol.fchk` |
| `cubefile` | 输出文件名，默认 `test.cube` | `phenol_homo.cub` |
| `npts` | 格点设置，见 4.3 | `0`、`-3`、`-1` |
| `format` | `h` 写文件头（默认），`n` 不写 | `h` |
| `cubefile2` | 从这个 cube 读取格点范围（此时 `npts` 必须为 `-1`） | 做差分密度、让 ESP 与密度共用格点 |

参数不区分大小写。不带任何参数运行时，cubegen 会询问 fchk 文件名，然后相当于执行 `cubegen 1 density=scf test.cube 80 h`。内存不够时用 `cubegen -m=2gb ...`（`-m` 写在最前面）或设置 `GAUSS_MEMDEF`。

### 4.2 kind 的取值

| kind | 含义 | 备注 |
|---|---|---|
| `MO=n` | 第 n 个分子轨道 | 也可写 `MO=Homo`、`MO=Lumo`、`MO=All`、`MO=OccA`、`MO=OccB`、`MO=Valence`、`MO=Virtuals` |
| `AMO=n` / `BMO=n` | 只取 α 或 β 轨道 | 开壳层体系中 `MO=Homo` 会同时输出 α 和 β 的 HOMO，想单独看时用 `AMO=Homo` |
| `Density=SCF` | 总电子密度 | 也可以是 `MP2`、`CC`、`CI` 等（需要计算时保存了对应的密度）；不支持 `Current` |
| `Spin=SCF` | 自旋密度 ρα − ρβ | 自由基、开壳层过渡金属 |
| `Alpha=SCF` / `Beta=SCF` | α 或 β 密度 | |
| `Potential=SCF` | 静电势 | 最耗时，用 `nprocs` 并行 |
| `Gradient` | 密度及其梯度（每点 4 个值） | |
| `Laplacian` | 密度的拉普拉斯 ∇²ρ | |
| `NormGradient` | 密度梯度的模 | |
| `CurrentDensity=I` | 磁感应电流密度（GIAO），I 为外场方向 X/Y/Z | 需要 NMR 计算的 chk |
| `ShieldingDensity=IJN` | 第 N 个原子的屏蔽密度 | 同上 |

### 4.3 npts：格点怎么选

| npts | 含义 | 适用 |
|---|---|---|
| `0` | 默认：程序自动确定的矩形盒子，每边 80 个点（不一定是正方体） | 快速看轨道 |
| 正整数 N | 每边 N 个点，例如 `100` 共 10⁶ 个点 | 想控制总点数时 |
| `-2` / `-3` / `-4` | Coarse / Medium / Fine：每 Bohr 3 / 6 / 12 个点 | `-3` 适合大多数作图；`-4` 用于定量积分或出版级细节 |
| ≤ `-5` | 点间距为 \|npts\|×10⁻³ Å，例如 `-80` 表示 0.08 Å | 想精确指定间距时 |
| `-1` | 从标准输入或 `cubefile2` 读入格点定义 | 多个 cube 共用格点、只算某个区域 |
| `-5` | 特殊：从标准输入读任意点列表（Å），输出每行“坐标 + 数值” | 只要少数几个点的值 |

> [!WARNING]
> `-5` 这个值有两种含义：作为“间距”它是 0.005 Å（极密，文件巨大）；官方文档同时把 `-5` 定义为“读入任意点列表”。不要用 `-5` 去指定间距，用 `-10` 及以下即可。

`npts=-1` 时，从标准输入读入的格点定义格式为：

```text
IFlag, X0, Y0, Z0      ! output unit (negative = formatted file) and origin
N1, X1, Y1, Z1         ! number of points and step vector along direction 1
N2, X2, Y2, Z2
N3, X3, Y3, Z3
```

`N1 < 0` 表示坐标单位为 Bohr，否则为 Å；三个方向按给定向量使用，不做正交化。实际中更常用的是 `cubefile2` 方式：先生成一个带文件头的密度 cube，再让 ESP 按它的格点生成。

### 4.4 常用命令

```bash
# orbitals: Medium grid (6 points/Bohr)
cubegen 0 MO=Homo phenol.fchk phenol_homo.cub -3 h
cubegen 0 MO=Lumo phenol.fchk phenol_lumo.cub -3 h
cubegen 0 MO=25   phenol.fchk phenol_mo25.cub  -3 h     # a specific orbital index

# total density, then ESP on exactly the same grid (needed for ESP-mapped surfaces)
cubegen 0 Density=SCF   phenol.fchk phenol_dens.cub -3 h
cubegen 8 Potential=SCF phenol.fchk phenol_esp.cub  -1 h phenol_dens.cub

# spin density of a radical (open-shell fchk)
cubegen 0 Spin=SCF radical.fchk radical_spin.cub -3 h
```

HOMO 的轨道序号可以在 fchk 中查 `Number of alpha electrons`（闭壳层时 HOMO 序号就等于它），也可以在 .log 的 `Alpha  occ. eigenvalues` 行数里数。示例脚本 `formchk_cubegen_homo_lumo.sh` 会自动打印 HOMO、LUMO 能量和能隙。

### 4.5 cube 的头部与检查

cube 文件的结构（示意，数值仅供参考）：

```text
 <title line 1>
 <title line 2: what was computed>
   -13   -6.123456   -7.234567   -5.345678    1      <- NAtoms (negative for MO cubes), origin (Bohr), NVal
    85    0.166667    0.000000    0.000000           <- N1 and step vector 1 (Bohr)
    90    0.000000    0.166667    0.000000           <- N2
    70    0.000000    0.000000    0.166667           <- N3
     8    8.000000    4.607680   -0.745012   -0.645706   <- Z, nuclear charge, x, y, z (Bohr)
   ...
     1   25                                          <- MO cubes only: number of MOs and their indices
  1.23456E-05  2.34567E-05 ...                       <- data, 6 values per line
```

要点：

- 文件中的长度单位都是 Bohr，数值都是原子单位；
- 轨道 cube 的原子数为负，并在原子坐标后多一行轨道编号；
- 第 3–6 行可以直接看出格点数和间距，用来核对两个 cube 是否“同格点”（做加减前必须一致）。

### 4.6 cubman：cube 的加减与统计

```bash
cubman      # interactive; actions: add, copy, diff, prop, subtract, scale, square, ...
```

典型用途是**差分密度** Δρ = ρ(AB) − ρ(A) − ρ(B)：三个 cube 必须用同一套格点（第一个用 `-3`，后两个用 `npts=-1` 加 `cubefile2`），然后 `subtract` 两次。cubman 会询问输入文件、是否为格式化文件（回答 `y`）和输出文件名。用 `prop` 可以得到 cube 的积分、最大/最小值等统计量。Multiwfn 主功能 13 和 Python 工具 Cube-Toolz 也能做同样的运算。

---

## 5. GaussView 6 中的可视化

GaussView 能直接打开 .log、.chk（同一平台）和 .fchk。跨平台时一律用 .fchk。下面的菜单路径以 GaussView 6.x 为准，不同小版本的按钮名称可能略有差别。

### 5.1 分子轨道

1. File → Open 打开 `phenol.fchk`；
2. Edit → MOs 打开 MO Editor，左侧列表显示每个轨道的编号、能量和占据情况；
3. 选中 HOMO、LUMO（可多选），在 Visualize 标签页设置等值面（Isovalue），点 Update 生成表面；
4. 在主窗口用 View → Display Format 调整表面透明度和颜色。

轨道等值面常用 0.02–0.05 a.u.。太小时表面会包住整个分子，看不出节面；太大时只剩零碎的小块。

### 5.2 电子密度与 ESP 映射表面

Results → Surfaces/Contours 打开对话框：

1. Cube Actions → New Cube，Type 选 Total Density，Grid 选 Medium 或 Fine，生成密度 cube；
2. Surface Actions → New Mapped Surface，选择用 ESP 着色，GaussView 会在密度等值面上计算 ESP；
3. 在映射表面的设置中把密度等值面改为 0.001 a.u.（与 Multiwfn 的默认范德华表面一致，便于比较），并设定对称的色标范围，例如 −0.05 到 +0.05 a.u.；
4. 色标两端的颜色要在图注中写清楚（GaussView 默认红色为负、蓝色为正）。

GaussView 也可以直接读入外部 cube（File → Open 选择 .cub），因此可以用 cubegen 生成高质量格点后再在 GaussView 中作图。

> [!NOTE]
> 论文中定量讨论 ESP（极值点位置和数值、正负面积）时，推荐用 Multiwfn 主功能 12（7.5 节）。GaussView 的映射表面更适合快速定性查看。

### 5.3 振动动画与红外谱

打开频率作业的 .log 或 .fchk，Results → Vibrations：

- 列表显示每个简正模式的频率、IR 强度（有 Raman 时还有 Raman 活性），**负频率即虚频**；
- 选中一个模式，点 Start 播放动画；勾选 Show Displacement Vectors 显示位移矢量；
- Spectrum 按钮显示 IR/Raman 谱，可以设置峰宽和频率缩放因子；
- 对过渡态，务必播放虚频动画，确认它对应你关心的反应坐标（[第7章](07-ts-irc.md)）。

### 5.4 优化、扫描与 IRC 曲线

| 菜单 | 用途 | 说明 |
|---|---|---|
| Results → Optimization | 优化过程中能量和梯度随步数的变化 | 判断是振荡还是缓慢下降（[第5章](05-opt.md)），也可配合仓库脚本 [`opt_monitor.py`](../scripts/opt_monitor.py) |
| Results → Scan | 扫描曲线 | 点击曲线上的点查看结构；二维扫描请打开 fchk（C.01 起保存多维信息，见 [第8章](08-scan.md)） |
| Results → IRC/Path | IRC 能量曲线 | 两端结构可另存为新输入，再做端点优化 |
| Results → UV-VIS / NMR 等 | 激发态和 NMR 谱 | 可设置展宽；NMR 可选择参比物 |

曲线窗口中右键菜单可以把数据导出为文本或保存图片。批量提取扫描和 IRC 能量可用仓库脚本 [`scan_irc_profile.py`](../scripts/scan_irc_profile.py)。

### 5.5 导出图片

File → Save Image 保存当前视图，可设置分辨率（出版级建议 ≥ 300 dpi，或直接按期刊要求的像素宽度）和背景色。GaussView 的渲染质量有限，追求更好效果时用 VMD（第 9 节）。

---

## 6. Multiwfn 入门

### 6.1 版本与安装

Multiwfn 是卢天（Tian Lu）开发的免费波函数分析程序。最后一个数字编号版本是 3.8（2026 年 1 月），**此后改为按发布日期命名**：本文写作时（2026-11）官网最新版为 2026.11.3，2026.10.18 等也都是 2026 年发布的版本。下载页：<http://sobereva.com/multiwfn/>；在线手册：<http://sobereva.com/multiwfn/manual/>。

| 平台 | 安装要点 |
|---|---|
| Windows | 解压即可用；编辑同目录的 `settings.ini`：`nthreads=` 设为物理核数，`cubegenpath=` 填 cubegen 的完整路径（用 cubegen 加速 ESP 计算） |
| Linux（完整版） | 需要 motif 库（CentOS/Rocky 用 `yum install motif`，Ubuntu 用 `sudo apt-get install libxm4 libgl1`）；在 `~/.bashrc` 中加入下列环境变量 |
| Linux（noGUI 版） | 无图形依赖，适合集群；不能弹出图形窗口，但所有计算和文件导出功能都能用 |
| macOS | 官方不提供，社区维护编译方案 [multiwfn-mac-build](https://github.com/digital-chemistry-laboratory/multiwfn-mac-build) |

Linux 下 `~/.bashrc` 的典型设置（路径换成你的解压目录）：

```bash
export Multiwfnpath=$HOME/soft/Multiwfn_2026.11.3_bin_Linux_noGUI
export PATH=$PATH:$Multiwfnpath
export OMP_STACKSIZE=1000M     # stack size per OpenMP thread
ulimit -s unlimited            # remove the main stack limit
```

然后 `chmod +x $Multiwfnpath/Multiwfn`。通过无图形界面的远程终端使用完整版时，如果载入文件后卡顿约两秒，手册建议加 `export DISPLAY=":0"`。

> [!TIP]
> 让 Multiwfn 调用 Gaussian 的 cubegen 计算 ESP，在核数不多（< 10 核）时明显更快：在 `settings.ini` 里把 `cubegenpath` 设为 cubegen 的路径即可，此时输入文件必须是 fch/fchk。G16 的 cubegen 会使用 `nthreads` 个核。大 fchk 调用 cubegen 失败时，先 `export GAUSS_MEMDEF=4GB`。

### 6.2 喂给 Multiwfn 什么文件

| 分析 | 推荐输入 | Gaussian 端需要的设置 |
|---|---|---|
| 轨道、密度、ESP、电荷、键级、IGMH、表面分析 | `.fchk` | `%chk` + formchk |
| IR/Raman/VCD/ROA 谱 | 频率作业的 `.log` | `Freq`（Raman 需要 `Freq=Raman`） |
| UV-Vis/ECD 谱 | TD 作业的 `.log` | `TD(NStates=...)` |
| NMR 谱 | NMR 作业的 `.log` | `NMR` |
| 激发态分析（空穴-电子、NTO、IFCT 等） | 先载入 TD 作业的 `.fchk`，进入功能后再给 `.log` | route 中加 `IOp(9/40=4)`，让 Gaussian 打印所有绝对值大于 10⁻⁴ 的组态系数 |
| 其他程序的波函数 | `.molden`、`.wfx`、`.mwfn` | — |

> [!IMPORTANT]
> 激发态分析时，fchk 和 log 必须对应**同一几何、同一方法和基组**。例如 fchk 来自 B3LYP/6-31G* 基态计算，而 TD 是在 PBE0/6-31G* 下做的，结果就毫无意义。最稳妥的做法是直接用 TD 作业自己的 chk 生成 fchk。

使用赝势（ECP）的体系，波函数里缺少内层电子的密度。Multiwfn 在读取 fchk 等文件时可以自动补充内核电子密度（EDF），细节见手册 6.4 节。对 ESP、电荷等受内层电子影响的分析要特别注意这一点。

### 6.3 主菜单速查

Multiwfn 启动后先要求输入文件路径，然后显示主菜单，输入数字进入对应功能：

| 主功能 | 内容 | 本章例子 |
|---|---|---|
| 0 | 显示分子结构、查看轨道等值面（GUI） | 交互检查轨道 |
| 5 | 计算任意实空间函数的格点数据并导出 cube | ELF、LOL、密度梯度等 |
| 7 | 布居分析与原子电荷（Hirshfeld=1、VDD=2、Mulliken=5、Löwdin=6、ADCH=11、CHELPG=12、MK=13、AIM=14、Hirshfeld-I=15、CM5=16、RESP=18、MBIS=20 等） | 7.1 |
| 8 | 轨道成分分析（Mulliken/SCPA、NAO、Hirshfeld、Becke 等） | 7.2 |
| 9 | 键级分析（Mayer=1、多中心键级=2、Wiberg（Löwdin 正交基）=3、Fuzzy=7、Laplacian=8 等） | 7.3 |
| 11 | 绘制光谱：1 IR、2 Raman、3 UV-Vis、4 ECD、5 VCD、6 ROA、7 NMR | 7.4 |
| 12 | 分子表面定量分析（默认映射 ESP） | 7.5 |
| 13 | 处理格点数据（cube 加减、截取、积分等） | 差分密度 |
| 18 | 电子激发分析（1 空穴-电子、4 Δr、6 NTO、8 IFCT、14 Λ 指数等） | 7.6 |
| 20 | 弱相互作用可视化（1 NCI、4 IRI、10 IGM、11 IGMH、−10 mIGM 等） | 7.7 |
| 100 | 其他功能，其中子功能 2 可导出 wfx、molden、fch、mwfn 等格式，或生成量化程序输入 | 9.3 |

括号中的子功能编号来自 Multiwfn 在线手册对应章节的标题。Multiwfn 版本更新频繁，**以你所用版本屏幕上的菜单为准**。

---

## 7. 典型分析实战

下面以示例目录中的苯酚（`phenol.fchk`，B3LYP/def2-TZVP 单点）和水二聚体为例。启动方式统一为：

```bash
Multiwfn phenol.fchk          # interactive
```

### 7.1 原子电荷（主功能 7）

```text
7        <- population analysis and atomic charges
11       <- ADCH charges
```

之后程序会询问原子密度的来源等问题，按屏幕提示选择即可（一般选内置的球形化自由原子密度）。完成后屏幕打印每个原子的电荷，并可选择导出 `.chg` 文件（元素、坐标、电荷），该文件可以被 VMD 等程序读取，也可以供主功能 12 用“基于原子电荷的 ESP”快速估算。

选择建议：

- **ADCH、CM5**：基组依赖小，适合比较不同分子的电荷；
- **CHELPG、MK、RESP**：重现 ESP，适合给分子力场提供电荷（RESP 拟合的完整流程见 [第12章](12-population-charges.md)）；
- **Mulliken**：基组依赖严重，加弥散函数后经常不合理，不推荐用于定量讨论。

### 7.2 轨道成分（主功能 8）

```text
8        <- orbital composition analysis
8        <- Hirshfeld method (atom/fragment contributions)
```

随后输入轨道序号即可得到各原子的贡献百分比（输入的轨道号可在主功能 0 或 fchk 中查到）。也可以先定义片段（例如苯环的 6 个 C，或配合物中的金属），得到片段贡献，用于判断“HOMO 主要在哪个部分”。对弥散基组，Mulliken 类方法可能给出负贡献，此时改用 Hirshfeld 或 Becke 方法。

### 7.3 键级（主功能 9）

```text
9        <- bond order analysis
1        <- Mayer bond order
```

Mayer 键级对单/双/三键分别接近 1/2/3，适合比较同一基组下的一系列分子。`3`（Löwdin 正交基下的 Wiberg 键级）与 NBO 中的 Wiberg 键级概念相同但数值不完全一致；`7`（Fuzzy 键级）对基组更稳定。定性上，苯环 C–C 的 Mayer 键级介于单键和双键之间（芳香性），苯酚的 C–O 键级略高于典型单键（O 孤对电子与苯环共轭）；具体数值随方法和基组变化，请以自己的计算为准。

### 7.4 绘制光谱（主功能 11）

载入**频率或 TD 作业的 .log**：

```bash
Multiwfn phenol_td.log
```

```text
11       <- plot spectrum
3        <- UV-Vis
0        <- show the spectrum (GUI version)
```

谱图界面中的常用选项（编号来自手册 3.13 节）：

| 选项 | 作用 |
|---|---|
| `0` | 画图（GUI 窗口），同时在终端打印谱峰位置 |
| `1` | 保存图片（格式由 `-4` 选择，推荐矢量格式 pdf/svg） |
| `2` | 导出离散线和展宽曲线数据：`spectrum_line.txt`、`spectrum_curve.txt` |
| `3` / `4` | 设置 X 轴、左 Y 轴范围和刻度 |
| `6` | 选择展宽函数（Gaussian、Lorentzian、Pseudo-Voigt） |
| `8` | 设置半峰全宽 FWHM |
| `10` | UV-Vis/ECD 切换能量单位（eV、nm、1000 cm⁻¹） |
| `12` | X 方向平移（与实验对齐时用，必须在文中说明） |
| `14` | 频率/激发能缩放因子 |
| `15` | 输出各个跃迁对谱的贡献（找出每个峰由哪些激发态构成） |
| `19` | Raman 活性转换为强度（需输入激发光波数和温度） |
| `20` | 修改跃迁强度（例如只保留 S₁ 画荧光谱） |

IR 谱选 `1`，展宽默认用 Lorentzian；ECD 选 `4`，会询问用长度还是速度表示的转动强度（推荐速度表示，见 [第10章](10-tddft.md)）；NMR 选 `7`，在其中可以设置参比物的屏蔽值把横轴换成化学位移。多构象加权谱的做法是写一个列出各构象文件和权重的文本文件作为输入，手册 4.11.4 节有完整例子。

仓库也提供了不依赖 Multiwfn 的 Python 脚本：[`spectrum_uvvis_ecd.py`](../scripts/spectrum_uvvis_ecd.py)（UV-Vis/ECD）和 [`spectrum_ir_raman.py`](../scripts/spectrum_ir_raman.py)（IR/Raman），用 `-h` 查看用法。

### 7.5 分子表面静电势分析（主功能 12）

```bash
Multiwfn phenol.fchk
```

```text
12       <- quantitative analysis of molecular surface
0        <- start analysis with default settings (rho = 0.001 a.u. surface, mapped function = ESP)
```

计算完成后屏幕打印表面极小点、极大点和汇总。下面是 Multiwfn 手册中苯酚（B3PW91/6-31G** 波函数）的结果摘录，用于说明各项的含义；你的数值会因方法和基组不同而有差异：

```text
Global surface minimum: -0.041203 a.u. at   1.455097   3.343708  -0.007902 Ang.
Global surface maximum:  0.085761 a.u. at  -1.936645   3.093464   0.021360 Ang.
...
Overall surface area:         476.05682 Bohr^2  ( 133.30951 Angstrom^2)
Positive surface area:        231.09186 Bohr^2  (  64.71232 Angstrom^2)
Negative surface area:        244.96497 Bohr^2  (  68.59719 Angstrom^2)
Overall average value:   -0.00020233 a.u. (  -0.12695332 kcal/mol)
...
Overall variance (sigma^2_tot):  0.00041488 a.u.^2 ( 163.34165024 (kcal/mol)^2)
Balance of charges (nu):   0.18762182
Molecular polarity index (MPI):   0.50154872 eV (     11.56600 kcal/mol)
```

| 输出项 | 含义 | 怎么用 |
|---|---|---|
| Global surface minimum / maximum | 表面 ESP 全局极小、极大值及位置 | 苯酚的极小在 O 孤对电子附近，极大在羟基 H 外侧：分别是氢键受体和供体位点 |
| Number of surface minima/maxima | 所有局部极值点（带 `*` 的是全局极值） | 讨论 π 电子区域、σ-hole 等 |
| Positive/Negative surface area | ESP 为正/负的表面积 | 比较分子极性分布 |
| Overall variance σ²tot、ν | ESP 方差与正负平衡度 | Politzer 等人的 GIPF 方法中用于关联沸点、溶解度等性质 |
| MPI | 分子极性指数 | 不同分子的整体极性比较 |
| Polar/Nonpolar surface area | \|ESP\| 大于/小于 10 kcal/mol 的面积 | 极性表面积 |

后处理菜单中的常用选项（手册 3.15.4 节）：

| 选项 | 作用 |
|---|---|
| `0` | GUI 中查看结构和极值点（红球为极大，蓝球为极小） |
| `1` | 极值点导出为 `surfanalysis.txt` |
| `2` | 极值点导出为 `surfanalysis.pdb`（B 因子列为 ESP 值），配合 VMD 作图 |
| `6` | 所有表面顶点导出为 `vtx.pdb` |
| `7` | 所有顶点及 ESP 值导出为 `vtx.txt`，下次分析可直接读入，省去重算 |
| `-2` | 用于生成表面的密度格点导出为 `surf.cub` |
| `3` / `4` | 删除指定数值范围内的极小/极大点（例如中性分子中 ESP 为正的极小点通常没有化学意义） |

> [!TIP]
> - 推荐用 B3LYP/def2-TZVP 等水平的波函数做 ESP 分析；6-31G** 也能得到定性正确的结果。
> - 阴离子或改用更小的密度等值面时，表面可能被格点盒子截断，此时在主功能 12 的选项 4（高级选项）中加大格点的延伸范围。
> - 对称分子中对称等价的极值点数值会有细微差别，这是表面格点不满足分子对称性造成的；需要时用选项 3 减小格点间距。

### 7.6 激发态分析（主功能 18）

TD 作业需要在 route 中加 `IOp(9/40=4)`，并保留它的 chk。示例目录中的 `phenol_td.gjf` 已经这样写好。

```bash
formchk phenol_td.chk phenol_td.fchk
Multiwfn phenol_td.fchk
```

```text
18       <- electron excitation analysis
1        <- hole-electron analysis (then input the path of phenol_td.log and the state index)
```

进入子功能后，Multiwfn 会询问 TD 作业 .log 的路径，然后询问要分析第几个激发态；载入后会打印组态系数归一化的偏差（偏差大说明没有加 `IOp(9/40=4)`）。之后在子菜单中可以计算空穴和电子的分布、重叠指数 Sr、质心距离 D 和 t 指数，并导出 cube 用 VMD 作图。其他常用子功能：

| 子功能 | 内容 |
|---|---|
| 4 | Δr 指数，衡量电荷转移距离（Δr 大于约 2 Å 时通常认为是 CT 激发） |
| 6 | 生成自然跃迁轨道（NTO），并导出 fch/molden 用于作图 |
| 8 | IFCT 方法计算片段间电荷转移量 |
| 14 | Λ 指数 |
| 2 | 原子/片段跃迁矩阵热图 |

Gaussian 本身也能生成 NTO（`Pop=SaveNTO Density=(Transition=n)`，见 [第10章](10-tddft.md)），两者可以互相对照。

### 7.7 弱相互作用可视化：IGMH（主功能 20）

水二聚体的原子顺序是：1–3 为氢键供体水分子，4–6 为受体水分子。

```bash
Multiwfn water_dimer.fchk
```

```text
20       <- visual study of weak interaction
11       <- IGMH analysis
2        <- define two fragments
1-3      <- atoms of fragment 1
4-6      <- atoms of fragment 2
2        <- medium-quality grid
3        <- export cube files to the current folder
```

当前目录会生成若干 cube，其中作图需要的是 `sl2r.cub`（sign(λ₂)ρ）和 `dg_inter.cub`（δg^inter）。把这两个文件和 Multiwfn `examples` 目录下的 `IGM_inter.vmd` 脚本放到 VMD 的工作目录，在 VMD 控制台运行 `source IGM_inter.vmd`，就会得到用 sign(λ₂)ρ 着色的 δg^inter 等值面。O–H···O 氢键区域应显示为蓝色的小圆盘。

在后处理菜单中还可以计算原子和原子对的 δG 指数（手册中的选项 6），用来定量比较各原子对相互作用的贡献。

> [!TIP]
> - 大体系只关心片段间作用时，可以在 `settings.ini` 中把 `IGMvdwscl` 设为 2.0 左右，只计算碎片重叠区域的格点，耗时能减少数倍（此时只能画 δg^inter）。
> - 只有结构、没有波函数的大体系（例如 MD 快照），可以用 mIGM（子功能 −10），它基于预设的原子密度，计算成本与 IGM 相当，但图像质量明显更好。
> - NCI（子功能 1）是更早的方法，作图时用 RDG 等值面；它分不出片段内和片段间作用，所以研究分子间作用时 IGMH 通常更清楚。

---

## 8. Multiwfn 命令流批处理

### 8.1 原理

Multiwfn 是交互式程序，每次“输入一行再回车”。把这些输入按顺序写进一个文本文件，再用重定向喂给程序，就实现了非交互运行：

```bash
Multiwfn phenol.fchk -silent < multiwfn_esp_surface_cmds.txt > phenol_esp_surface.out
```

- `<` 把命令流文件作为标准输入，`>` 把屏幕输出存到文件，这是操作系统的功能，与 Multiwfn 无关；
- 输入文件名写在命令行上，**不要**写进命令流文件；
- `-silent`（或在 `settings.ini` 中设 `isilent=1`）禁止弹出图形窗口，否则程序会停在图形窗口等你手动关闭；
- 命令流结束后，如果程序仍在等待输入，会打印 `forrtl: severe (24): end-of-file during read` 之类的信息后退出。手册说明这不是真正的错误；想“优雅退出”就在最后回到主菜单并输入 `q`。

> [!WARNING]
> 命令流文件中**不要写注释**。Multiwfn 会把整行（包括你以为的注释）当作输入读入，读文件名、原子范围（如 `1-3`）这类字符串时尤其会出错。说明文字请写在 README 中。

### 8.2 写命令流的正确方法

1. 先**交互运行一遍**，把每一步输入记在纸上或文本里；
2. 按顺序写成命令流文件，每行一个输入；需要“直接回车”的地方写一个空行；
3. 用一个文件测试，检查输出文件中每一步的提示和你的输入是否对得上；
4. Multiwfn 升级后重新检查一遍（菜单编号偶尔会变）。

Linux 下也可以不写文件，直接用 `echo -e` 或 here-document：

```bash
# ADCH charges for every fchk in the folder
for f in *.fchk; do
  echo -e "7\n11\n1\ny\n0\nq" | Multiwfn "$f" -silent > "${f%.fchk}_adch.out"
done
```

上面 `7\n11` 之后的几项（`1`、`y`、`0`、`q`）对应 ADCH 子功能中的后续提示，**需要你先交互运行一次确认**它们在你的版本中的含义，再放进脚本。

Windows 命令提示符同样支持 `<` 和 `>`；PowerShell 不支持 `<`，要改写为：

```powershell
Get-Content multiwfn_esp_surface_cmds.txt | .\Multiwfn.exe phenol.fchk > phenol_esp_surface.out
```

### 8.3 批处理脚本与工具

```bash
#!/bin/bash
# run the same Multiwfn command stream on all fchk files, keep per-file outputs
cmd=multiwfn_esp_surface_cmds.txt
for f in *.fchk; do
  base=${f%.fchk}
  echo "== $base"
  Multiwfn "$f" -silent < "$cmd" > "${base}_esp.out" 2>&1
  for x in surfanalysis.txt surfanalysis.pdb vtx.pdb; do
    [ -f "$x" ] && mv "$x" "${base}_$x"     # outputs have fixed names: rename them
  done
done
```

Multiwfn 导出的文件名大多是固定的（`surfanalysis.txt`、`spectrum_curve.txt`、`sl2r.cub` 等），批处理时一定要在每次运行后立即改名，否则会被下一个文件覆盖。

现成工具：

- 仓库脚本 [`batch_cube_multiwfn.sh`](../scripts/batch_cube_multiwfn.sh)：批量 formchk + cubegen（HOMO/LUMO/密度/ESP），并可按命令流文件调用 Multiwfn（`bash ../scripts/batch_cube_multiwfn.sh -h` 查看用法）；
- [BaneWfn](https://github.com/bane-dysta/banewfn)：中文文档的 Multiwfn 批处理工作流工具，把固定的菜单序列模板化；
- [runMultiwfn.bash](https://github.com/polyluxus/runMultiwfn.bash)：Linux 下非交互运行 Multiwfn 的 bash 封装。

在集群上运行 Multiwfn 时，用 noGUI 版本，在作业脚本中设置 `settings.ini` 的 `nthreads` 与申请的核数一致（可以把 settings.ini 复制到作业目录，用 `sed -i 's/nthreads= *[0-9]*/nthreads= 16/' settings.ini` 修改这份私有副本；当前目录中存在 settings.ini 时 Multiwfn 通常优先读取它，请在输出开头确认实际使用的线程数）。作业脚本模板见 [第17章](17-hpc-jobs.md)。

---

## 9. VMD 渲染与波函数格式导出

### 9.1 用 VMD 画 cube 等值面

[VMD](https://www.ks.uiuc.edu/Research/vmd/) 是免费的分子可视化程序，渲染质量远好于 GaussView。基本步骤：

1. File → New Molecule 载入 cube 文件（文件类型选 Gaussian Cube）；
2. Graphics → Representations：保留一个 CPK 或 Licorice 表示结构，新建一个 Isosurface 表示，Draw 选 Solid Surface，Isovalue 设为 0.02（轨道）或 0.001（密度）；
3. 轨道要正负两个相位：再建一个 Isosurface，Isovalue 取负值，换一种颜色；
4. ESP 映射表面：载入密度 cube 和 ESP cube（作为同一分子的两个 volume），Isosurface 的 Coloring Method 选 Volume，指定用 ESP 那个 volume 着色，并在 Trajectory 标签页设置色标范围（例如 −0.03 到 0.03 a.u.）；
5. Display 中设置正交投影（Orthographic）、关闭坐标轴、白色背景；File → Render 选择 Tachyon（或 TachyonInternal）渲染，可把输出分辨率设大一些。

Multiwfn 作者提供了大量现成的 VMD 脚本和图文教程（例如 IGMH 用的 `IGM_inter.vmd`，以及 ESP 着色分子表面的教程），放在 Multiwfn 的 `examples` 目录和手册 4.A 节中，强烈建议直接使用。

### 9.2 Output=WFX：写出 AIM 风格波函数

有些程序（AIMAll、部分 QTAIM 工具）只认 wfn/wfx。Gaussian 用 `Output=WFX`（或 `Output=WFN`）写出，**文件名作为单独一行写在分子坐标之后的输入段中**：

```text
%chk=water_wfx.chk
#p B3LYP/def2TZVP Output=WFX

Water SP, write a .wfx file

0 1
O    0.000000    0.000000    0.117300
H    0.000000    0.757200   -0.469200
H    0.000000   -0.757200   -0.469200

water.wfx

```

注意坐标后的空行、文件名行、文件末尾的空行都不能少。wfx 相比 wfn 是较新的格式，推荐优先使用；Multiwfn 读 fchk 时信息最全，能读 fchk 时就不必转换为 wfx。

### 9.3 其他格式转换

- **Multiwfn 主功能 100 → 子功能 2**：把当前载入的波函数导出为 wfx、wfn、molden、fch、mwfn 等格式，或生成其他量化程序的输入；
- [Molden2AIM](https://github.com/zorkzou/Molden2AIM)：把 molden 转成 wfn/wfx/NBO 的 .47 文件；
- [iodata](https://github.com/theochem/iodata)：Python 库，可读写 fchk、log、wfn/wfx、molden、cube；
- [Cube-Toolz](https://github.com/funkymunkycool/Cube-Toolz)：cube 加减、积分、平面平均（无许可证，只作链接参考）；
- [cclib](https://github.com/cclib/cclib)：解析 log/fchk 中的能量、几何、频率、激发态和轨道，适合在 Python 中批量提取数据。

---

## 10. MOKIT：把 Gaussian 轨道交给其他程序

[MOKIT](https://github.com/1234zou/MOKIT)（邹靖翔开发，GitHub 与 [GitLab](https://gitlab.com/jxzou/mokit) 上均有仓库）可以在 Gaussian 和 PySCF、ORCA、GAMESS、OpenMolcas、PSI4、Molpro、BDF 等程序之间传递分子轨道。最常见的场景是：在 Gaussian 中得到稳定的 broken-symmetry UKS 或 UNO 轨道（[第13章](13-openshell-stability.md)），再交给多参考方法（CASSCF/NEVPT2）作为初猜。

| 命令 | 方向 |
|---|---|
| `fch2mkl` | Gaussian → ORCA（生成 .mkl，再用 ORCA 的 `orca_2mkl` 转为 .gbw） |
| `bas_fch2py` | Gaussian → PySCF（生成可运行的 PySCF 脚本） |
| `fch2inp` | Gaussian → GAMESS |
| `fch2inporb` | Gaussian → (Open)Molcas |
| `fch2psi` | Gaussian → PSI4 |
| `mkl2fch`、`molden2fch` | ORCA/其他程序 → Gaussian |

使用流程示例：

```bash
formchk bs_uks.chk bs_uks.fch
fch2mkl bs_uks.fch          # writes ORCA .mkl and a template .inp; see MOKIT docs for options
```

MOKIT 文档建议在 Gaussian 输入中加入 `nosymm int=nobasistransform`，以便轨道顺利传递并帮助后续 SCF 收敛。各工具的完整选项见 [MOKIT 文档](https://jeanwsr.gitlab.io/mokit-doc-mdbook/)。

---

## 参数选择与收敛建议

| 项目 | 推荐值 | 说明 |
|---|---|---|
| 轨道 cube 格点 | `npts=-3`（6 点/Bohr） | 默认 `0`（80³）对小分子够用；大分子用 `-3` 更均匀 |
| 定量积分用的 cube | `npts=-4` 或更密 | 差分密度积分、电荷转移量计算 |
| ESP cube | 与密度 cube 共用格点（`-1` + `cubefile2`），`nprocs` 设为可用核数 | 做 ESP 映射表面必须同格点 |
| 轨道等值面 | 0.02–0.05 a.u. | 同一篇文章内保持一致 |
| 密度等值面（范德华表面） | 0.001 a.u. | Bader 定义，Multiwfn 主功能 12 默认值 |
| ESP 色标 | 对称范围，如 ±0.03–0.05 a.u.（约 ±20–30 kcal/mol） | 图注中写明数值和颜色含义 |
| ESP 分析的波函数 | B3LYP/def2-TZVP（6-31G** 可定性） | 阴离子需要弥散函数 |
| 激发态分析 | TD 作业加 `IOp(9/40=4)`，用 TD 作业自己的 fchk | 体系大时可改 `IOp(9/40=3)` |
| IGMH 格点 | Medium（选项 2）作图；High 用于 δG 指数 | 大体系配合 `IGMvdwscl=2.0` |
| IGMH 等值面 | δg^inter = 0.005–0.01 a.u. | 想只显示最强作用时适当调大 |
| Multiwfn 并行 | `nthreads` = 物理核数 | 集群上与作业申请核数一致 |

---

## 常见问题与报错

**1. GaussView 或 Multiwfn 打不开 .chk**

- 原因：chk 是平台相关的二进制文件，从 Linux 集群拷到 Windows 后无法读取；或者 chk 由其他 Gaussian 版本生成。
- 处理：在生成它的机器上 `formchk`，拷贝 .fchk；如需在另一台机器上续算，再在那台机器上 `unfchk`。

**2. 读取旧 chk 时报版本或格式相关的错误（`Geom=Check`、`Guess=Read` 失败）**

- 原因：C.01 起 chk 带版本号，C.01 之前生成的 chk 会被拒绝或读错。
- 处理：`c8616 old.chk` 原地升级（先备份），或在旧环境 `formchk`、新环境 `unfchk -fixver`。

**3. `cubegen MO=Homo` 生成的文件包含两个轨道**

- 原因：开壳层体系中 `Homo` 同时选中 α 和 β 的 HOMO。
- 处理：用 `AMO=Homo` 或 `BMO=Homo`，或用具体的轨道序号。

**4. ESP cube 生成极慢**

- 原因：ESP 每个格点都要做积分，且 `nprocs` 写成了 0 或 1。
- 处理：`cubegen 16 Potential=SCF ...` 使用多核；格点不要过密；只需要表面上的 ESP 时直接用 Multiwfn 主功能 12（只在表面顶点上计算）。

**5. cubegen 或 formchk 报内存不足**

- 处理：`formchk -m=2gb ...`、`cubegen -m=2gb ...`（`-m` 必须写在最前面，C.01 起支持），或 `export GAUSS_MEMDEF=2GB`。

**6. 差分密度 cube 相减后全是噪声或报格点不一致**

- 原因：三个 cube 的格点不同。
- 处理：第一个 cube 用固定的 `npts`，其余用 `npts=-1` 并指定第一个 cube 为 `cubefile2`；三个计算的原子坐标也必须完全相同（片段计算要用复合物中的几何）。

**7. Multiwfn 终端最后打印 `forrtl: severe (24): end-of-file during read`**

- 原因：命令流结束时程序仍在等待输入。
- 处理：无害，可以忽略；或在命令流最后回到主菜单并输入 `q`。

**8. Linux 下 Multiwfn 计算大体系时崩溃（segmentation fault）**

- 原因：栈空间不足。
- 处理：确认 `ulimit -s unlimited` 和 `export OMP_STACKSIZE=1000M` 已生效（作业脚本中也要写）；检查 `/proc/sys/kernel/shmmax` 是否过小。

**9. 空穴-电子分析给出“Deviation to expected normalization value”很大**

- 原因：Gaussian 默认只打印绝对值大于 0.1 的组态系数。
- 处理：TD 作业的 route 中加 `IOp(9/40=4)` 重新计算；并确认 fchk 与 log 来自同一几何和同一水平。

**10. 批处理后发现所有分子的结果都一样**

- 原因：Multiwfn 的输出文件名固定，后一个分子覆盖了前一个；或命令流中写了文件名，导致每次都读同一个文件。
- 处理：每次运行后立刻改名（见 8.3 节）；输入文件名只写在命令行上。

**11. `Output=WFX` 作业正常结束但没有 wfx 文件**

- 原因：没有在坐标后的输入段中写文件名行，或缺少结尾空行。
- 处理：按 9.2 节的格式补全。

---

## 练习

1. 对示例中的苯酚分别用 `npts=0`、`-3`、`-4` 生成 HOMO cube，比较文件大小和 GaussView 中的显示效果。
2. 用 `cubegen` 生成苯酚的密度和同格点 ESP，在 VMD 中画 ESP 映射表面；再用 Multiwfn 主功能 12 找到极值点，比较两者给出的最负区域是否一致。
3. 在 B3LYP/def2-SVP 和 B3LYP/def2-TZVP 两个水平下计算苯酚的 ADCH 和 Mulliken 电荷，比较哪种电荷对基组更敏感。
4. 写一个命令流，对 3 个不同取代苯酚批量运行 ESP 表面分析，提取全局极小值，讨论取代基对 O 孤对电子 ESP 的影响。
5. 把水二聚体换成甲酸二聚体，用 IGMH 比较两个 O–H···O 氢键，并计算原子对 δG 指数。
6. 用 `freqchk` 把苯酚的热化学从 298.15 K 改为 373.15 K，与直接在输入中加 `Temperature=373.15` 的结果比较。

---

## 参考资料

**Gaussian 官方文档**

- [Gaussian Utilities 总览](https://gaussian.com/utils/)：formchk、unfchk、c8616、chkchk、cubegen、cubman、freqchk、newzmat、testrt 等
- [formchk](https://gaussian.com/formchk/)、[unfchk](https://gaussian.com/unfchk/)、[c8616](https://gaussian.com/c8616/)、[chkchk](https://gaussian.com/chkchk/)
- [cubegen](https://gaussian.com/cubegen/)：kind、npts、cube 文件格式
- [cubman](https://gaussian.com/cubman/)、[freqchk](https://gaussian.com/freqchk/)、[newzmat](https://gaussian.com/newzmat/)
- [Output 关键词（WFN/WFX）](https://gaussian.com/output/)
- [g16 命令行选项（含 -fchk）](https://gaussian.com/options/)、[Link 0 命令](https://gaussian.com/link0/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes](https://gaussian.com/relnotes/)：工具的 `-m` 选项、chk 版本号、多维扫描信息写入 fchk
- [GaussView 6](https://gaussian.com/gaussview6/)、[GaussView 6.1.1 Release Notes](https://gaussian.com/gv611rn/)

**Multiwfn 与可视化**

- [Multiwfn 主页与下载](http://sobereva.com/multiwfn/)、[在线手册](http://sobereva.com/multiwfn/manual/)：重点阅读 2.1（安装）、3.13（光谱）、3.15（表面分析）、3.21（激发态分析）、3.23（弱相互作用）、5.2–5.3（静默与批处理）、5.7（调用 cubegen）
- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [MOKIT 文档](https://jeanwsr.gitlab.io/mokit-doc-mdbook/)

**相关 GitHub 仓库**（仅链接，未复制代码）

- [1234zou/MOKIT](https://github.com/1234zou/MOKIT)：Gaussian 与 PySCF/ORCA/Molcas 等程序间的轨道传递，适合把 BS/UNO 轨道交给多参考方法
- [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn)：中文文档的 Multiwfn 批处理工作流工具
- [polyluxus/runMultiwfn.bash](https://github.com/polyluxus/runMultiwfn.bash)：非交互运行 Multiwfn 的 bash 封装
- [digital-chemistry-laboratory/multiwfn-mac-build](https://github.com/digital-chemistry-laboratory/multiwfn-mac-build)：macOS 上编译 Multiwfn
- [zorkzou/Molden2AIM](https://github.com/zorkzou/Molden2AIM)：molden → wfn/wfx/NBO .47
- [theochem/iodata](https://github.com/theochem/iodata)：fchk、wfx、molden、cube 等格式的 Python 读写库
- [funkymunkycool/Cube-Toolz](https://github.com/funkymunkycool/Cube-Toolz)：cube 加减、积分、平面平均（无许可证，只外链）
- [cclib/cclib](https://github.com/cclib/cclib)：Python 解析 Gaussian log/fchk
- [aligfellow/xyzrender](https://github.com/aligfellow/xyzrender)：从 xyz 生成出版级分子结构图
- [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py)：含 cube 输出、光谱作图等功能的 Python 工具包

**关键文献**

- T. Lu, F. Chen, Multiwfn: A multifunctional wavefunction analyzer, *J. Comput. Chem.* **2012**, 33, 580. DOI: [10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885)
- T. Lu, A comprehensive electron wavefunction analysis toolbox for chemists, Multiwfn, *J. Chem. Phys.* **2024**, 161, 082503. DOI: [10.1063/5.0216272](https://doi.org/10.1063/5.0216272)
- T. Lu, F. Chen, Quantitative analysis of molecular surface based on improved Marching Tetrahedra algorithm, *J. Mol. Graph. Model.* **2012**, 38, 314. DOI: [10.1016/j.jmgm.2012.07.004](https://doi.org/10.1016/j.jmgm.2012.07.004)
- T. Lu, Q. Chen, Independent gradient model based on Hirshfeld partition: A new method for visual study of interactions in chemical systems, *J. Comput. Chem.* **2022**, 43, 539. DOI: [10.1002/jcc.26812](https://doi.org/10.1002/jcc.26812)
- E. R. Johnson, S. Keinan, P. Mori-Sánchez, J. Contreras-García, A. J. Cohen, W. Yang, Revealing noncovalent interactions, *J. Am. Chem. Soc.* **2010**, 132, 6498. DOI: [10.1021/ja100936w](https://doi.org/10.1021/ja100936w)
- W. Humphrey, A. Dalke, K. Schulten, VMD: Visual molecular dynamics, *J. Mol. Graph.* **1996**, 14, 33. DOI: [10.1016/0263-7855(96)00018-5](https://doi.org/10.1016/0263-7855(96)00018-5)

---

[← 上一章](15-oniom.md) | [返回目录](../README.md) | [下一章 →](17-hpc-jobs.md)
