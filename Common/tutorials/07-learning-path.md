# 第7章 DFT 理论速览、学习路线与资源汇总（含机器学习势前沿）（DFT Primer, Learning Path, Convergence Checklist, Resources and Machine-Learned Potentials）

本章是整个仓库的“导航页”：先用最少的公式讲清密度泛函理论（DFT）、Kohn–Sham 方程、交换关联泛函（Jacob 阶梯）、平面波/高斯/数值原子轨道基组与赝势，再给出覆盖 VASP、Gaussian、CP2K、OpenMX 四个板块的学习顺序和 12 周学习计划。之后把四款软件的收敛测试（ENCUT、k 点、SIGMA、CUTOFF/REL_CUTOFF、scf.energycutoff、基组与积分格点）归纳成一张通用清单，并给出自动建目录、提交、收集、画图的模板。最后介绍机器学习势前沿（VASP MLFF、NEP/GPUMD、DeePMD-kit、MACE 与 MACE-MP 基础模型），汇总教材、课程、awesome 列表、官方论坛与中文社区，并说明引用规范和许可证注意事项（DFT tutorial roadmap, convergence tests, MLIP, learning resources）。

> [!NOTE]
> 本章不追求理论推导的完整性，目标是让你知道“每个输入参数背后对应哪一个近似”，以及“下一步该读哪一章”。严格的推导请看第 6.1 节列出的教材。配套示例：[07-convergence-templates](../examples/07-convergence-templates/README.md)（收敛测试模板）、[07-mlip-nep-mace](../examples/07-mlip-nep-mace/README.md)（机器学习势）、[07-learning-resources](../examples/07-learning-resources/README.md)（资源书签与学习打卡表）。

---

## 1. 本章目标与适用场景

| 目标 | 学完后你应当能够 |
|---|---|
| 理论框架 | 说清 Born–Oppenheimer 近似、Hohenberg–Kohn 定理、Kohn–Sham 方程、自洽场（SCF）各是什么 |
| 近似层级 | 知道 LDA/GGA/meta-GGA/杂化/双杂化各自的代价与典型用途，知道色散校正和 DFT+U 解决什么问题 |
| 基组与赝势 | 理解 ENCUT、CUTOFF/REL_CUTOFF、scf.energycutoff、Gaussian 基组为什么都需要做收敛测试 |
| 学习路线 | 按自己的课题选软件，按本仓库给出的顺序学习，不在不需要的章节上浪费时间 |
| 收敛测试 | 用统一的“1 meV/atom”思路为四款软件确定截断能、k 点、展宽、基组和格点 |
| 自动化 | 用脚本批量建目录、提交作业、收集能量并作图 |
| 机器学习势 | 知道 MLIP 能做什么、不能做什么；会用 MACE-MP 预优化结构和生成 NEB 初猜，会准备 NEP/MACE 训练集 |
| 资源与社区 | 找到合适的教材、课程、官方论坛和中文社区，学会正确提问 |
| 学术规范 | 正确引用软件与工具论文，理解“只链接、不复制”的许可证要求 |

适用读者：刚进组、需要从零上手的本科生和研究生；会用一款软件、想迁移到另一款的同学；想把机器学习势引入自己 DFT 工作流的读者。

---

## 2. 原理速览

### 2.1 从多体薛定谔方程到 Born–Oppenheimer 近似

含 $N$ 个电子、$M$ 个原子核的体系，非相对论哈密顿量为（原子单位）：

$$
\hat H = -\sum_{i}\frac{1}{2}\nabla_i^2 - \sum_{I}\frac{1}{2M_I}\nabla_I^2 - \sum_{i,I}\frac{Z_I}{|\mathbf r_i-\mathbf R_I|} + \sum_{i<j}\frac{1}{|\mathbf r_i-\mathbf r_j|} + \sum_{I<J}\frac{Z_IZ_J}{|\mathbf R_I-\mathbf R_J|}
$$

原子核比电子重三个数量级以上，**Born–Oppenheimer 近似**把核看作固定，只解电子问题，得到势能面 $E(\{\mathbf R_I\})$。之后所有“结构优化、过渡态、声子、分子动力学”都是在这个势能面上移动原子核；**机器学习势要拟合的也正是这张势能面**（第 5 节）。

### 2.2 Hohenberg–Kohn 定理与 Kohn–Sham 方程

Hohenberg–Kohn 定理（1964）指出：基态电子密度 $n(\mathbf r)$ 唯一决定外势，从而决定体系的所有基态性质；基态能量是密度的泛函，并在真实密度处取最小值。Kohn–Sham 方法（1965）引入一组无相互作用的辅助电子，使其密度等于真实密度：

$$
E[n] = T_s[n] + \int v_{\mathrm{ext}}(\mathbf r)\,n(\mathbf r)\,d\mathbf r + E_{\mathrm H}[n] + E_{\mathrm{xc}}[n]
$$

$$
\Big[-\tfrac{1}{2}\nabla^2 + v_{\mathrm{ext}}(\mathbf r) + v_{\mathrm H}[n](\mathbf r) + v_{\mathrm{xc}}[n](\mathbf r)\Big]\psi_i(\mathbf r) = \varepsilon_i\,\psi_i(\mathbf r),\qquad n(\mathbf r)=\sum_i f_i\,|\psi_i(\mathbf r)|^2
$$

其中 $T_s$ 是无相互作用动能，$E_{\mathrm H}$ 是经典库仑（Hartree）能，所有“不知道怎么写”的多体效应都塞进交换关联能 $E_{\mathrm{xc}}$。由于 $v_{\mathrm H}$ 和 $v_{\mathrm{xc}}$ 依赖于 $n$，方程必须**自洽求解**：猜密度 → 构造势 → 解本征方程 → 得新密度 → 混合 → 重复，直到能量/密度变化小于阈值。这就是各软件里 SCF 参数（VASP `EDIFF`、Gaussian `SCF=Tight`、CP2K `EPS_SCF`、OpenMX `scf.criterion`）的来源。

> [!IMPORTANT]
> Kohn–Sham 本征值 $\varepsilon_i$ 严格来说只有最高占据轨道（与电离能相关）有明确物理意义。用 GGA 本征值算出的带隙普遍偏小，这是方法本身的误差，不是你算错了。需要准确带隙时，用杂化泛函、meta-GGA（如 mBJ）或 GW，见 [VASP 第6章](../../VASP/tutorials/06-hybrid-metagga.md)。

### 2.3 交换关联泛函：Jacob 阶梯

Perdew 把泛函按所用信息的多少排成“Jacob 阶梯”，越往上越贵，一般也越准（但不保证）：

| 阶梯 | 类型 | 依赖的量 | 代表泛函 | 典型用途 |
|---|---|---|---|---|
| 1 | LDA | $n$ | PZ81、PW92 | 历史基准；均匀体系 |
| 2 | GGA | $n,\nabla n$ | PBE、PBEsol、BLYP、RPBE | 固体与表面的默认选择（PBE），晶格常数可选 PBEsol |
| 3 | meta-GGA | $+\ \tau$（动能密度） | SCAN、r2SCAN、TPSS、M06-L | 更好的结构与能量；mBJ 用于带隙 |
| 4 | 杂化 | $+$ 精确交换 | B3LYP、PBE0、HSE06、M06-2X、ωB97X-D | 分子化学的主力（Gaussian）；固体带隙（HSE06） |
| 5 | 双杂化 / RPA | $+$ 虚轨道 | B2PLYP、DSD-PBEP86、RPA | 高精度热化学，代价最高 |

两类常见“补丁”：

- **色散校正**：GGA/杂化泛函描述不了长程范德华作用。Grimme D3/D3(BJ)/D4 在能量上加经验项，VASP 用 `IVDW`，Gaussian 用 `EmpiricalDispersion=GD3BJ`，CP2K 用 `&VDW_POTENTIAL`，OpenMX 用 `scf.dftD`。层状材料、吸附、分子晶体基本都要加。详见 [VASP 第8章](../../VASP/tutorials/08-vdw-2d.md)、[CP2K 第8章](../../CP2K/tutorials/08-dispersion-spin-dftu.md)。
- **DFT+U**：给局域的 d/f 电子加 Hubbard 修正，缓解自相互作用误差（过渡金属氧化物、稀土化合物）。U 值要说明来源（线性响应或文献），见 [VASP 第7章](../../VASP/tutorials/07-magnetism-u-soc.md)、[OpenMX 第13章](../../OpenMX/tutorials/13-dft-u.md)。

### 2.4 基组与赝势：四款软件的根本差别

Kohn–Sham 轨道要用一组基函数展开。基组和核区处理方式不同，决定了每款软件“要收敛哪个参数”：

| 软件 | 基组 | 核区处理 | 需要收敛的离散化参数 | 详细章节 |
|---|---|---|---|---|
| VASP | 平面波 $e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}$ | PAW | `ENCUT`、k 点（KPOINTS/`KSPACING`） | [VASP 第3章](../../VASP/tutorials/03-convergence.md) |
| Gaussian | 原子中心高斯函数（Pople、Dunning、def2） | 全电子；重元素用 ECP | 基组大小、DFT 积分格点 | [Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md) |
| CP2K | 高斯轨道 + 辅助平面波表示密度（GPW） | GTH 赝势（GAPW 可全电子） | `CUTOFF`、`REL_CUTOFF`、基组、k 点 | [CP2K 第4章](../../CP2K/tutorials/04-cutoff-convergence.md) |
| OpenMX | 数值赝原子轨道（PAO） | 模守恒赝势（VPS） | PAO 截断半径与轨道数、`scf.energycutoff`、`scf.Kgrid` | [OpenMX 第4章](../../OpenMX/tutorials/04-pao-vps.md)、[第5章](../../OpenMX/tutorials/05-scf-kpoints-convergence.md) |

平面波基组只由一个数控制：保留满足

$$
\tfrac{1}{2}|\mathbf k+\mathbf G|^2 < E_{\mathrm{cut}}
$$

的所有平面波。$E_{\mathrm{cut}}$ 越大基组越完备，能量单调下降（变分）。基组不完备还会引起 **Pulay 应力**：晶胞体积改变时平面波数目固定、有效截断能随之改变，使能量–体积曲线不光滑、平衡体积偏小（第 4.3 节）。原子中心基组（Gaussian、CP2K、OpenMX）还存在**基组重叠误差（BSSE）**，算结合能时要做 Counterpoise 校正，见 [Gaussian 第14章](../../Gaussian/tutorials/14-weak-bsse.md)。

### 2.5 周期性、k 点与展宽

周期体系的波函数满足 Bloch 定理 $\psi_{n\mathbf k}(\mathbf r+\mathbf R)=e^{i\mathbf k\cdot\mathbf R}\psi_{n\mathbf k}(\mathbf r)$，总能量是对布里渊区的积分，实际用有限的 k 点网格求和。Monkhorst–Pack 网格的 k 点为

$$
\mathbf k_{n_1n_2n_3} = \sum_{i=1}^{3}\frac{2n_i-N_i-1}{2N_i}\,\mathbf b_i,\qquad n_i = 1,\dots,N_i
$$

$N_i$ 为偶数时网格不含 Γ 点；Γ 中心网格则总包含 Γ。六角晶系一般用 Γ 中心网格。实空间越大，倒空间越小，需要的 k 点越少：大超胞、表面 slab 的真空方向取 1。

金属在费米面附近占据数突变，k 点求和收敛很慢。**展宽（smearing）** 把占据数抹平，代价是引入一个虚拟的电子熵 $S$，体系的变分量变成 Mermin 自由能 $F = E - TS$（$T$ 即展宽宽度，VASP 中是 `SIGMA`）。展宽越大 k 点收敛越快，但 $TS$ 项越大、能量偏离越多，所以要检查每原子 $|TS|$ 是否足够小（第 4.2 节）。

### 2.6 DFT 能做什么、要警惕什么

| 适合 | 需要额外处理 | DFT 本身不擅长 |
|---|---|---|
| 结构、晶格常数、相对能量、反应热 | 弱相互作用（加色散校正） | 强关联体系（Mott 绝缘体、重费米子，需 DFT+U/DMFT） |
| 能带、态密度、电荷分布 | 带隙（杂化泛函/mBJ/GW） | 激发态的精确描述（TD-DFT 有局限，电荷转移态需长程校正泛函） |
| 声子、弹性、振动光谱 | 过渡金属氧化物（DFT+U） | 多参考体系（双自由基、键断裂，见 [Gaussian 第13章](../../Gaussian/tutorials/13-openshell-stability.md)） |
| 吸附能、反应能垒（定性到半定量） | 溶剂效应（隐式溶剂模型） | 需要化学精度（约 1 kcal/mol）的热化学，应改用高级波函数方法 |

> [!TIP]
> 判断“结果可信吗”的三个问题：(1) 计算参数是否收敛（第 4 节）；(2) 方法近似是否适合这个问题（本节）；(3) 是否与实验或更高精度的计算对比过。第一条是你能完全控制的，也是审稿人最常问的。

### 2.7 从 DFT 到机器学习势

DFT 计算量大致随电子数的三次方增长，几百个原子、皮秒量级的分子动力学已接近常规上限。机器学习原子间势（MLIP）利用相互作用的**局域性**，把总能量写成原子能量之和：

$$
E \approx \sum_{i=1}^{N}\varepsilon\big(\{\mathbf r_{ij}\}_{r_{ij}<r_c}\big),\qquad \mathbf F_i = -\nabla_{\mathbf r_i}E,\qquad \sigma = \frac{1}{V}\frac{\partial E}{\partial \boldsymbol\epsilon}
$$

$\varepsilon$ 由描述符（对称函数、ACE、Chebyshev 多项式等）加神经网络/线性回归给出，或由等变消息传递网络（如 MACE）直接学习。训练时最小化与 DFT 能量、力、应力（位力）的误差：

$$
\mathcal L = \lambda_e\,\mathrm{RMSE}(E) + \lambda_f\,\mathrm{RMSE}(\mathbf F) + \lambda_v\,\mathrm{RMSE}(W) + \text{正则项}
$$

关键认识：**MLIP 的精度上限就是训练数据的 DFT 精度**，而且只在训练数据覆盖的构型范围内可靠。所以 DFT 的收敛测试和设置一致性（第 4 节）同样是训练 MLIP 的前提。

---

## 3. 全仓库学习路线

### 3.1 路线总览

```text
                 +--------------------------------------------------+
   第 0 阶段     | Common 1-6：Linux/HPC、作业脚本、建模、ASE/pymatgen、 |
   （所有人）    | 后处理与可视化、工作流                               |
                 +-------------------------+------------------------+
                                           |
           +---------------+---------------+---------------+----------------+
           v               v                               v                v
   VASP（周期固体、   Gaussian（分子、         CP2K（大体系 AIMD、     OpenMX（磁性、
   表面、催化）       反应、光谱）             液体/界面、QM/MM）      输运、O(N)）
   1-5 入门           1-6 入门                 1-5 入门                1-5 入门
   6-10 进阶          7-12 进阶                6-11 进阶               6-10 进阶
   11-17 专题         13-16 专题               12-17 专题              11-17 专题
   18 排错            17-18 作业与排错          18 排错                 18 并行与排错
           |               |                               |                |
           +---------------+---------------+---------------+----------------+
                                           v
                 本章第 4 节收敛清单 + 第 5 节机器学习势（VASP 第15章、Scripts/）
```

各板块的完整章节列表见 [VASP 板块目录](../../VASP/README.md)、[Gaussian 板块目录](../../Gaussian/README.md)、[CP2K 板块目录](../../CP2K/README.md)、[OpenMX 板块目录](../../OpenMX/README.md)；机器学习势训练集与热输运脚本见 [Scripts/README.md](../../Scripts/README.md)；230 多个相关 GitHub 仓库的分类索引见 [Resources/github-repos.md](../../Resources/github-repos.md)。

### 3.2 第 0 阶段：通用基础（所有人，约 1–2 周）

| 顺序 | 章节 | 学完的标志 |
|---|---|---|
| 1 | [第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输](01-linux-hpc.md) | `ssh hpc` 免密登录，`dft` 环境装好并通过自检 |
| 2 | [第2章 HPC 作业调度：Slurm/PBS/LSF 模板与并行设置](02-job-scheduler.md) | 能写作业脚本、设置 MPI/OpenMP、用作业数组 |
| 3 | [第3章 结构准备：数据库下载、建模与格式转换](03-structure-modeling.md) | 能从 Materials Project 下载结构，搭超胞/表面/缺陷 |
| 4 | [第4章 ASE 与 pymatgen 实战](04-ase-pymatgen.md) | 能用 Python 读写四款软件的输入输出 |
| 5 | [第5章 后处理与可视化](05-postprocess-visualization.md) | 能画出版级能带/DOS 图 |
| 6 | [第6章 工作流自动化与高通量](06-workflow-automation.md) | （可选）知道什么时候该用 atomate2/AiiDA |

本板块的辅助脚本：[mkjob.py](../scripts/mkjob.py)（生成作业脚本）、[kpath_gen.py](../scripts/kpath_gen.py)（高对称路径）、[convert_any.py](../scripts/convert_any.py)（结构格式互转）、[job_status.py](../scripts/job_status.py)（批量查看作业状态）、[env_check.sh](../scripts/env_check.sh)（环境自检），都在 `Common/scripts/` 中，在 `Common/` 目录下用 `python scripts/<name> -h` 查看用法（`env_check.sh` 用 `bash scripts/env_check.sh -h`）。

### 3.3 选哪款软件

| 你的问题 | 首选 | 备选 | 理由 |
|---|---|---|---|
| 晶体结构、能带、态密度、缺陷、表面催化 | VASP | CP2K、OpenMX | PAW 精度高，后处理生态最完整 |
| 有机/无机分子、反应机理、过渡态、热化学 | Gaussian | CP2K（气相分子也可） | 杂化/双杂化泛函、频率、IRC、溶剂模型成熟 |
| UV-Vis、ECD、NMR 等分子光谱 | Gaussian | — | TD-DFT、GIAO 等功能完整 |
| 几百到上千原子的液体、界面、AIMD | CP2K | VASP + MLFF | GPW + OT 对大体系效率高 |
| 非共线磁性、交换耦合、电子输运（NEGF） | OpenMX | VASP（磁性） | 数值原子轨道 + NEGF 一体化，开源 |
| 大规模长时间 MD、热导率 | MLIP（NEP/GPUMD、DeePMD、MACE） | VASP MLFF | DFT 生成训练集，MLIP 跑生产 |

> [!TIP]
> 许可证也是选型因素：VASP 和 Gaussian 是商业软件，需要课题组购买许可证；CP2K（GPL）和 OpenMX（GPL）免费开源。没有 VASP 许可证时，周期体系可以先用 CP2K 或 OpenMX 学会全部流程。

### 3.4 VASP 学习顺序

| 阶段 | 章节 | 说明 |
|---|---|---|
| 入门 | [第1章 编译安装与作业提交](../../VASP/tutorials/01-setup-compile.md) → [第2章 四大输入文件](../../VASP/tutorials/02-input-output-files.md) → [第3章 收敛测试](../../VASP/tutorials/03-convergence.md) → [第4章 结构优化与状态方程](../../VASP/tutorials/04-relax-eos.md) → [第5章 能带与态密度](../../VASP/tutorials/05-band-dos.md) | 必修，按顺序做完每章的 Si 算例 |
| 进阶 | [第6章 杂化泛函与 meta-GGA](../../VASP/tutorials/06-hybrid-metagga.md)、[第7章 磁性/DFT+U/SOC](../../VASP/tutorials/07-magnetism-u-soc.md)、[第8章 范德华与二维材料](../../VASP/tutorials/08-vdw-2d.md)、[第9章 表面与功函数](../../VASP/tutorials/09-surface-workfunction.md)、[第10章 吸附与催化](../../VASP/tutorials/10-adsorption-catalysis.md) | 按课题选读 |
| 专题 | [第11章 NEB/Dimer](../../VASP/tutorials/11-neb-dimer.md)、[第12章 声子](../../VASP/tutorials/12-phonon.md)、[第13章 弹性](../../VASP/tutorials/13-elastic.md)、[第14章 光学](../../VASP/tutorials/14-optics-dielectric.md)、[第15章 AIMD 与 MLFF](../../VASP/tutorials/15-md-mlff.md)、[第16章 电荷与成键](../../VASP/tutorials/16-charge-bonding.md)、[第17章 高级电子结构](../../VASP/tutorials/17-advanced-electronic.md) | 用到再学 |
| 排错 | [第18章 后处理生态与报错排查](../../VASP/tutorials/18-postprocess-troubleshooting.md) | 遇到报错时查 |

### 3.5 Gaussian 学习顺序

| 阶段 | 章节 |
|---|---|
| 入门 | [第1章 第一个计算](../../Gaussian/tutorials/01-setup-first-job.md) → [第2章 输入文件结构](../../Gaussian/tutorials/02-input-anatomy.md) → [第3章 方法与基组](../../Gaussian/tutorials/03-method-basis.md) → [第4章 单点与 SCF](../../Gaussian/tutorials/04-sp-scf.md) → [第5章 几何优化](../../Gaussian/tutorials/05-opt.md) → [第6章 频率与热化学](../../Gaussian/tutorials/06-freq-thermo.md) |
| 进阶 | [第7章 过渡态与 IRC](../../Gaussian/tutorials/07-ts-irc.md)、[第8章 势能面扫描](../../Gaussian/tutorials/08-scan.md)、[第9章 溶剂模型](../../Gaussian/tutorials/09-solvation.md)、[第10章 TD-DFT](../../Gaussian/tutorials/10-tddft.md)、[第11章 NMR](../../Gaussian/tutorials/11-nmr-props.md)、[第12章 布居与电荷](../../Gaussian/tutorials/12-population-charges.md) |
| 专题 | [第13章 开壳层与 Broken-Symmetry](../../Gaussian/tutorials/13-openshell-stability.md)、[第14章 BSSE](../../Gaussian/tutorials/14-weak-bsse.md)、[第15章 ONIOM](../../Gaussian/tutorials/15-oniom.md)、[第16章 Multiwfn 后处理](../../Gaussian/tutorials/16-postprocess-multiwfn.md) |
| 作业与排错 | [第17章 作业脚本与并行](../../Gaussian/tutorials/17-hpc-jobs.md)、[第18章 报错大全与完整案例](../../Gaussian/tutorials/18-errors-capstone.md) |

### 3.6 CP2K 学习顺序

| 阶段 | 章节 |
|---|---|
| 入门 | [第1章 安装](../../CP2K/tutorials/01-install.md) → [第2章 输入结构](../../CP2K/tutorials/02-input-structure.md) → [第3章 基组与赝势](../../CP2K/tutorials/03-gpw-basis-pseudo.md) → [第4章 CUTOFF 收敛](../../CP2K/tutorials/04-cutoff-convergence.md) → [第5章 SCF 策略](../../CP2K/tutorials/05-scf-strategy.md) |
| 进阶 | [第6章 周期性与 k 点](../../CP2K/tutorials/06-periodic-surface-kpoints.md)、[第7章 几何与晶胞优化](../../CP2K/tutorials/07-geo-cell-opt.md)、[第8章 色散/自旋/DFT+U](../../CP2K/tutorials/08-dispersion-spin-dftu.md)、[第9章 杂化泛函与 ADMM](../../CP2K/tutorials/09-hybrid-admm.md)、[第10章 DOS/PDOS/能带](../../CP2K/tutorials/10-dos-pdos-band.md)、[第11章 cube 与功函数](../../CP2K/tutorials/11-cube-charge-workfunction.md) |
| 专题 | [第12章 AIMD](../../CP2K/tutorials/12-aimd.md)、[第13章 NEB 与 Dimer](../../CP2K/tutorials/13-neb-dimer.md)、[第14章 振动与热化学](../../CP2K/tutorials/14-vibration-thermo.md)、[第15章 元动力学与 PLUMED](../../CP2K/tutorials/15-metadynamics-plumed.md)、[第16章 QM/MM 与 xTB/DFTB](../../CP2K/tutorials/16-qmmm-xtb.md)、[第17章 自动化与工具生态](../../CP2K/tutorials/17-tools-ecosystem.md) |
| 排错 | [第18章 报错排查、性能调优与版本迁移](../../CP2K/tutorials/18-troubleshooting.md) |

> [!WARNING]
> CP2K 2026.2 改动了部分输入语法：`&MOTION/&CELL_OPT/TYPE` 被移除（旧输入会报 `The specified keyword 'TYPE' is not available anymore`），`&DFT/&PRINT/&PDOS` 改为 `&DFT/&PRINT/&DOS`（内含 `&PDOS`/`&CURVE`/`&BROADEN`）。网上 2025 年以前的教程照抄会出错，以 CP2K 板块中的新旧对照为准。

### 3.7 OpenMX 学习顺序

| 阶段 | 章节 |
|---|---|
| 入门 | [第1章 概述](../../OpenMX/tutorials/01-overview.md) → [第2章 编译](../../OpenMX/tutorials/02-install.md) → [第3章 第一个计算](../../OpenMX/tutorials/03-input-first-calc.md) → [第4章 PAO 与 VPS](../../OpenMX/tutorials/04-pao-vps.md) → [第5章 SCF 与 k 点](../../OpenMX/tutorials/05-scf-kpoints-convergence.md) |
| 进阶 | [第6章 优化与 NEB](../../OpenMX/tutorials/06-optimization-eos-neb.md)、[第7章 分子动力学](../../OpenMX/tutorials/07-md.md)、[第8章 能带](../../OpenMX/tutorials/08-band.md)、[第9章 DOS/COHP](../../OpenMX/tutorials/09-dos-pdos-cohp-optics.md)、[第10章 电荷与分子轨道](../../OpenMX/tutorials/10-charge-mo-cube.md) |
| 专题 | [第11章 非共线磁性](../../OpenMX/tutorials/11-spin-noncollinear.md)、[第12章 SOC 与拓扑](../../OpenMX/tutorials/12-soc-topology.md)、[第13章 DFT+U](../../OpenMX/tutorials/13-dft-u.md)、[第14章 Wannier 函数](../../OpenMX/tutorials/14-wannier.md)、[第15章 Berry 相与外加电场](../../OpenMX/tutorials/15-berry-efield.md)、[第16章 ESM 与功函数](../../OpenMX/tutorials/16-surface-esm-workfunction.md)、[第17章 NEGF 输运](../../OpenMX/tutorials/17-negf-transport.md) |
| 并行与排错 | [第18章 并行、O(N) 方法、作业脚本与常见错误](../../OpenMX/tutorials/18-parallel-on-errors.md) |

### 3.8 机器学习势路线

1. 先掌握一款 DFT 软件的 AIMD：[VASP 第15章](../../VASP/tutorials/15-md-mlff.md) 或 [OpenMX 第7章](../../OpenMX/tutorials/07-md.md)；
2. 体验“开箱即用”的 MLIP：本章第 5.2、5.3 节（MACE-MP 预优化、NEB 初猜）；
3. 训练专用势：VASP 内置 MLFF（[VASP 第15章](../../VASP/tutorials/15-md-mlff.md)），或 NEP/GPUMD 全流程（[Scripts/README.md](../../Scripts/README.md)，含训练集生成、主动学习与热输运分析），或 MACE/DeePMD（本章第 5.4 节）；
4. 验证：与 DFT 对比声子、RDF、扩散系数等（第 5.6 节）。

### 3.9 12 周学习计划示例

以“周期体系 + VASP”为主线，每周约 10 小时；Gaussian 方向把第 3–8 周换成 Gaussian 第1–12章即可。

| 周 | 内容 | 产出 |
|---|---|---|
| 1 | Common 第1–2章；读 Sholl & Steckel 第1–3章 | 能登录集群、提交 Si 单点 |
| 2 | Common 第3–4章；VASP 第1–2章 | 自己从 MP 下载结构并写出四大输入 |
| 3 | VASP 第3章 + 本章第 4 节；用 [07-convergence-templates](../examples/07-convergence-templates/README.md) 做一次完整收敛测试 | 收敛测试报告（图 + 选定参数） |
| 4 | VASP 第4章：EOS 与结构优化 | 晶格常数与实验/文献对比 |
| 5 | VASP 第5章 + Common 第5章 | 出版级能带 + PDOS 图 |
| 6 | VASP 第6–7章（按课题选） | HSE 带隙或磁性基态 |
| 7 | VASP 第8–9章 | 二维材料或表面 slab 模型 |
| 8 | VASP 第10章 | 一组吸附能 |
| 9 | VASP 第11章（NEB）+ 本章第 5.3 节 | 一个扩散/反应能垒 |
| 10 | VASP 第12章或第15章 | 声子谱或 AIMD 轨迹分析 |
| 11 | Common 第6章 + 本章第 5 节 | 批量计算脚本；MACE-MP 预优化体验 |
| 12 | 复现一篇文献中的关键图，整理成组会报告 | 可复现的计算目录（输入、脚本、版本记录） |

可打印的打卡表见 [07-learning-resources/study_checklist.md](../examples/07-learning-resources/study_checklist.md)。

---

## 4. 收敛测试通用清单

### 4.1 通用原则

1. **一次只改一个参数**，其余（泛函、PREC、展宽、赝势/基组、SCF 阈值）全部固定；
2. 以**每原子能量差**为主要判据，常用 1 meV/atom；SCF 阈值要比它严格约两个数量级；
3. 你真正关心的是**物理量的差值**（吸附能、能垒、形成能、晶格常数），它通常比总能量收敛得快，但不总是如此——关键量最好直接测试；
4. 最大值不是“真值”：参考点取扫描中最精确的一点，并确认最后两点彼此一致；
5. 不同参数相互影响：先 ENCUT/CUTOFF，再 k 点，再 SIGMA；金属体系 k 点和 SIGMA 要一起看；
6. 结构优化（特别是变胞）、声子、弹性常数对力和应力敏感，比单点能需要更高的精度；
7. 所有要相减的计算（体相、表面、分子、吸附构型）必须使用**同一套**参数。

| 判据 | 常用标准 | 适用 |
|---|---|---|
| 总能量（每原子） | < 1 meV/atom | 一般计算 |
| 能量差 | < 1–5 meV（或 < 0.1 kcal/mol） | 吸附能、反应能、构象能 |
| 力 | 收敛后变化 < 0.01 eV/Å | 优化、声子 |
| 应力/压力 | 变化 < 1 kbar 量级 | 变胞优化、EOS、弹性 |
| 熵项 \|T*S\| | < 1 meV/atom | VASP 展宽检查 |

### 4.2 VASP：ENCUT、k 点、ISMEAR/SIGMA

**ENCUT**：POTCAR 中每个元素给出 `ENMAX`，多元素体系取最大值。从 `max(ENMAX)` 起扫描到约 1.6–2.0 倍；**涉及体积/晶胞变化的计算，官方建议至少 1.3×ENMAX**（VASP Wiki “Volume relaxation”）。下面是本章模板中的单点 INCAR（[encut_scan/INCAR.template](../examples/07-convergence-templates/encut_scan/INCAR.template)）：

```bash
# ENCUT convergence test - single-point energy (VASP 6.4-6.6)
SYSTEM  = encut-scan
ISTART  = 0
ICHARG  = 2
PREC    = Accurate
ENCUT   = @ENCUT@         # replaced by the scan script
EDIFF   = 1E-7
ALGO    = Normal
LREAL   = .FALSE.
LASPH   = .TRUE.
ISMEAR  = 0
SIGMA   = 0.05
NSW     = 0
IBRION  = -1
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `PREC` | Accurate | 同时决定 FFT 网格；扫描中保持不变 |
| `ENCUT` | max(ENMAX) 到约 2×ENMAX | 扫描对象；变胞计算 ≥ 1.3×ENMAX |
| `EDIFF` | 1E-7 | SCF 噪声远小于 1 meV/atom |
| `LREAL` | .FALSE. | 倒空间投影，避免实空间投影误差混入 |
| `LASPH` | .TRUE. | d/f 元素和 meta-GGA 需要；对 Si 无害 |
| `ISMEAR`/`SIGMA` | 0 / 0.05（半导体）；1 / 0.2（金属） | 扫描 ENCUT 时固定 |
| `NSW`/`IBRION` | 0 / −1 | 单点，不移动原子 |

运行与收集（脚本：[run_encut_scan.sh](../examples/07-convergence-templates/encut_scan/run_encut_scan.sh)）：

```bash
cd Common/examples/07-convergence-templates/encut_scan
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR        # 自己生成，仓库不提供
grep ENMAX POTCAR                                       # ENMAX = 245.345; ...（以你的 POTCAR 为准）
bash run_encut_scan.sh prepare                          # 按 1.0-2.0 x ENMAX 建目录
bash run_encut_scan.sh run                              # 或 submit（每个点一个作业）
bash run_encut_scan.sh collect                          # -> encut.dat
grep "energy  without entropy" encut_*/OUTCAR | awk '{print $1, $NF}'
python ../plot_convergence.py encut.dat --xlabel "ENCUT (eV)"
```

**k 点**：两种写法二选一。

- KPOINTS 文件中的 Monkhorst–Pack 或 Γ 中心网格（六角晶系、含 Γ 的需求用 Gamma）；
- INCAR 中的 `KSPACING`（单位 Å⁻¹，**包含 2π**，默认 0.5）。VASP 按 $N_i=\max\big(1,\lceil |\mathbf b_i|\,2\pi/\mathrm{KSPACING}\rceil\big)$ 生成网格，`KGAMMA = .TRUE.`（默认）时以 Γ 为中心。只有目录中**没有 KPOINTS 文件**时才使用 `KSPACING`。

`KSPACING` 的好处是不同晶胞得到相近的 k 点密度，适合高通量；缺点是相邻取值可能给出相同网格。[kpoint_scan.py](../examples/07-convergence-templates/kpoint_scan.py) 会按 POSCAR 计算每个 KSPACING 对应的网格并跳过重复点：

```bash
cd ..                                   # back to 07-convergence-templates/
python kpoint_scan.py prepare --template encut_scan --encut 400 --kspacing 0.5 0.4 0.3 0.25 0.2 0.15
python kpoint_scan.py prepare --template encut_scan --encut 400 --mesh 4 6 8 10 12 14 16 --outdir kscan_mesh
python kpoint_scan.py collect --outdir kscan
python plot_convergence.py kpoints.dat --xlabel "N1*N2*N3" --logx
```

**ISMEAR/SIGMA**：半导体/绝缘体用 `ISMEAR = 0`、`SIGMA = 0.01–0.05`；金属优化用 `ISMEAR = 1` 或 2、`SIGMA` 约 0.1–0.2；精确总能和 DOS 用四面体法 `ISMEAR = -5`（需足够 k 点）。判据是**每原子 \|T*S\| < 1 meV**，OUTCAR 中对应：

```bash
grep "EENTRO" OUTCAR | tail -1          # entropy T*S    EENTRO = ...（整个晶胞，eV）
python ../vasp_tsigma_check.py sigma_*/OUTCAR
```

（示意输出，数值仅供参考）

```text
file                                     ISMEAR  SIGMA  T*S/atom(meV)         E0(eV)  ok?
sigma_0.20/OUTCAR                             1    0.2         -0.412     -3.745012  yes
sigma_0.40/OUTCAR                             1    0.4         -2.871     -3.746530  NO -> lower SIGMA (and densify k-points)
```

完整讲解（含 Al 金属算例、PREC/LREAL/ADDGRID 等精度标签）见 [VASP 第3章](../../VASP/tutorials/03-convergence.md)。

### 4.3 力与应力收敛、Pulay 应力

- 结构优化的终止条件 `EDIFFG`（负值表示力判据，如 −0.01 eV/Å）只有在 ENCUT、k 点足够时才有意义：截断不足时力本身就不准；
- **Pulay 应力**来自平面波基组不完备，总是使计算得到的平衡体积偏小。检验方法：在优化后的结构上用同一 ENCUT 再算一次单点，看 OUTCAR 中 `external pressure` 是否接近 0；或比较 ENCUT 和 1.3×ENCUT 下的压力差；
- 变胞优化（`ISIF = 3`）的标准做法：ENCUT 取 ≥ 1.3×ENMAX；优化结束后**从 CONTCAR 重新开始再优化一次**，因为优化过程中平面波数目固定、实际截断随体积变化；
- 声子和弹性常数对力的精度要求更高（`EDIFF` 1E-8、`ADDGRID`），见 [VASP 第12章](../../VASP/tutorials/12-phonon.md)、[第13章](../../VASP/tutorials/13-elastic.md)。

```bash
grep "external pressure" OUTCAR | tail -1    # external pressure = 1.23 kB  Pullay stress = 0.00 kB
NIONS=$(grep -m1 NIONS OUTCAR | awk '{print $NF}')
grep -A $((NIONS+1)) "TOTAL-FORCE" OUTCAR | tail -$NIONS | awk '{print sqrt($4^2+$5^2+$6^2)}' | sort -g | tail -1
```

（OUTCAR 中原文就拼作 `Pullay stress`，对应 INCAR 标签 `PSTRESS`，不是自动修正量。）

### 4.4 CP2K：CUTOFF 与 REL_CUTOFF

CP2K 的密度在多级实空间网格上表示：第 $i$ 级网格截断为 $E^{(i)} = \texttt{CUTOFF}/\alpha^{\,i-1}$（$\alpha$ 为 `PROGRESSION_FACTOR`，默认 3）。高斯函数乘积（指数 $\zeta$）被分配到满足 $E^{(i)} \ge \zeta\cdot\texttt{REL\_CUTOFF}$ 的最粗网格上。因此：

1. **先固定 `REL_CUTOFF`（如 60 Ry）扫 `CUTOFF`**，再**固定收敛的 `CUTOFF` 扫 `REL_CUTOFF`**；
2. 除能量外，还要看输出 `MULTIGRID INFO` 中各级网格上的高斯个数（`count for grid`）：只升高 CUTOFF 而不升高 REL_CUTOFF 时，高斯会逐渐“移到”粗网格，最细网格可能变空，能量不再变化，但这并不代表精度提高；
3. 再看网格上的**电子数积分**：`Electronic density on regular grids:` 一行的第一个数是网格上积分得到的电子密度（CP2K 中电子带负号，应非常接近 −价电子数），第二个数是它与精确电子数之差。8 原子 Si 有 32 个价电子（GTH-PBE-q4），差值越接近 0，密度在网格上的表示越好；差值明显偏大（如 10⁻⁵ 量级）说明截断不足，详见 [CP2K 第4章](../../CP2K/tutorials/04-cutoff-convergence.md)。

模板 [cp2k_cutoff_scan/template.inp](../examples/07-convergence-templates/cp2k_cutoff_scan/template.inp) 的关键部分：

```text
@SET CUTOFF 400
@SET REL_CUTOFF 60
&GLOBAL
  RUN_TYPE    ENERGY
  PRINT_LEVEL MEDIUM            ! prints MULTIGRID INFO
&END GLOBAL
&FORCE_EVAL
  &DFT
    &MGRID
      NGRIDS     4
      CUTOFF     ${CUTOFF}
      REL_CUTOFF ${REL_CUTOFF}
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      EPS_SCF 1.0E-7
    &END SCF
  &END DFT
&END FORCE_EVAL
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `CUTOFF` | 扫描 200–800 Ry | 最细网格的平面波截断（Ry） |
| `REL_CUTOFF` | 扫描 30–100 Ry | 决定高斯被映射到哪一级网格 |
| `NGRIDS` | 4（默认） | 网格级数 |
| `EPS_DEFAULT` | 1.0E-12 | 积分筛选阈值；太松会掩盖截断的影响 |
| `PRINT_LEVEL` | MEDIUM | 才会打印 `count for grid` |

```bash
cd Common/examples/07-convergence-templates/cp2k_cutoff_scan
bash run_scan.sh cutoff run && bash run_scan.sh cutoff collect
CUT=500 bash run_scan.sh rel run && CUT=500 bash run_scan.sh rel collect
grep "ENERGY| Total FORCE_EVAL" cut_500_rel_60/output.out
grep "count for grid" cut_500_rel_60/output.out
grep "Electronic density on regular grids" cut_500_rel_60/output.out
python ../plot_convergence.py cutoff.dat rel.dat --unit hartree --labels CUTOFF REL_CUTOFF
```

能量单位是 Hartree，`plot_convergence.py --unit hartree` 会换算成 meV/atom。CELL_OPT 和 NPT 模拟需要比单点更高的截断（应力对网格更敏感），详见 [CP2K 第4章](../../CP2K/tutorials/04-cutoff-convergence.md)。

### 4.5 OpenMX：scf.energycutoff 与 scf.Kgrid

OpenMX 的基组（PAO）本身由 `Definition.of.Atomic.Species` 中的轨道名决定（如 `Si7.0-s2p2d1`：截断半径 7.0 bohr，s、p 各 2 个、d 1 个）；`scf.energycutoff`（Ry，默认 150）只控制实空间积分网格的密度，`scf.Kgrid` 控制 k 点。两者都要做收敛测试，基组本身的测试见 [OpenMX 第4章](../../OpenMX/tutorials/04-pao-vps.md)。

```text
scf.energycutoff            @ECUT@      # Ry, scanned 120-400
scf.Kgrid                   @KGRID@     # e.g. 8 8 8 while scanning the cutoff
scf.criterion               1.0e-9      # Hartree
```

```bash
cd Common/examples/07-convergence-templates/openmx_cutoff_scan
# 先编辑 template.dat 中的 DATA.PATH
bash run_scan.sh ecut run  && bash run_scan.sh ecut collect
EFIX=220 bash run_scan.sh kgrid run && EFIX=220 bash run_scan.sh kgrid collect
grep "Utot\." scan_ecut/ecut_220.out
python ../plot_convergence.py ecut.dat --unit hartree --xlabel "scf.energycutoff (Ry)"
```

由于实空间网格与原子位置的相对关系会造成能量随平移轻微振荡（“egg-box”效应），截断偏低时力也会有噪声；结构优化和声子要比单点用更高的 `scf.energycutoff`。完整算例见 [OpenMX 第5章](../../OpenMX/tutorials/05-scf-kpoints-convergence.md)。

### 4.6 Gaussian：基组与积分格点

分子计算没有 k 点和平面波截断，“收敛测试”变成：

- **基组**：沿同一系列逐级增大（def2-SVP → def2-TZVP → def2-QZVP；cc-pVDZ → cc-pVTZ → cc-pVQZ），看你关心的相对能量（反应能、能垒、构象能）是否稳定。阴离子、弱相互作用、激发态需要弥散函数（aug-、def2-…D）；
- **DFT 积分格点**：Gaussian 16 默认 `Int=UltraFine`（剪枝 99,590），Gaussian 09 默认 `FineGrid`（75,302）。meta-GGA（如 M06 系列）对格点敏感，低频振动、构象能差很小或优化不易收敛时用 `Int=SuperFine`。比较能量时，所有计算的格点必须一致。

模板 [gaussian_basis_scan/template.gjf](../examples/07-convergence-templates/gaussian_basis_scan/template.gjf) 用两步 Link1 分别优化 HCN 和 HNC，扫描它们的异构化能：

```text
%nprocshared=8
%mem=8GB
%chk=@TAG@_hcn.chk
#p PBE1PBE/@BASIS@ Opt Int=@GRID@ SCF=Tight

HCN, basis/grid scan (@BASIS@, @GRID@)

0 1
H    0.000000    0.000000   -1.066000
C    0.000000    0.000000    0.000000
N    0.000000    0.000000    1.156000

--Link1--
(second step: HNC, same route section)
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `PBE1PBE` | — | Gaussian 中 PBE0 的关键词 |
| `@BASIS@` | def2SVP/def2TZVP/def2QZVP、cc-pVnZ | 扫描对象 |
| `Int=@GRID@` | FineGrid / UltraFine / SuperFine | 扫描对象；G16 默认 UltraFine |
| `SCF=Tight` | — | 收紧 SCF，避免噪声 |
| `Opt` | — | 每个基组各自优化几何，得到自洽的比较 |

```bash
cd Common/examples/07-convergence-templates/gaussian_basis_scan
bash run_scan.sh run          # 或 prepare 后自行提交，作业脚本见 Gaussian 第17章
bash run_scan.sh collect      # -> basis_scan.dat
grep "SCF Done" def2TZVP_UltraFine.log | tail -2
grep -c "Normal termination" def2TZVP_UltraFine.log     # 应为 2
```

判据由问题决定：例如反应能随基组变化 < 0.5 kJ/mol、格点变化 < 0.1 kJ/mol。详见 [Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md) 和 [第4章](../../Gaussian/tutorials/04-sp-scf.md)。

### 4.7 自动化：建目录 → 提交 → 收集 → 画图

[vaspup2.0](https://github.com/kavanase/vaspup2.0) 的思路很值得借鉴：把模板输入和一个 `CONFIG` 文件放在 `input/` 中，脚本生成每个扫描点的目录、提交作业，结束后再用一条命令汇总能量和力。本章模板沿用同一思路，四款软件的脚本都分为 `prepare / run / collect` 三步（VASP ENCUT 脚本另有 `submit`）。VASP、CP2K、OpenMX 的汇总表前三列统一为“参数 能量 原子数”（`encut.dat` 第 4 列附带压力，`kpoints.dat` 附带 k 点数和网格），可直接交给 `plot_convergence.py`；Gaussian 的 `basis_scan.dat` 则直接给出两个异构体的能量和异构化能（kJ/mol），不经过 `plot_convergence.py`：

```text
# x        E_total        NATOMS
250        -10.81234567   2
...
```

提交可以用作业数组，一次排队、限制并发（[slurm_array_dirs.sh](../examples/07-convergence-templates/slurm_array_dirs.sh)）：

```bash
ls -d encut_scan/encut_* > dirs.txt
sbatch --array=1-$(wc -l < dirs.txt)%4 --export=ALL,RUN="mpirun vasp_std" slurm_array_dirs.sh
squeue -u $USER
python ../../scripts/job_status.py -h       # 批量检查各目录是否正常结束
```

最后用 [plot_convergence.py](../examples/07-convergence-templates/plot_convergence.py) 统一作图，它会打印“从这一点往后都在阈值以内”的最小参数值（示意输出，数值仅供参考）：

```text
== encut_scan/encut.dat  (reference x = 490)
           x    dE (meV/atom)
         250            9.812  *
         270            3.204  *
         290            0.951
         ...
converged (|dE| < 1.0 meV/atom from here on): x = 290
```

更完整的高通量方案（custodian 自动纠错、atomate2 收敛工作流、AiiDA-VASP 的 converge WorkChain）见 [第6章](06-workflow-automation.md)。

---

## 5. 机器学习势前沿

### 5.1 概览：四类常用方案

| 方案 | 类型 | 训练/运行环境 | 适合 | 入口 |
|---|---|---|---|---|
| VASP MLFF | 核方法 + 贝叶斯误差估计，在线（on-the-fly）训练 | VASP ≥ 6.3（`ML_MODE` 自 6.4 起） | 已有 VASP 许可证、想在 AIMD 中直接加速 | [VASP 第15章](../../VASP/tutorials/15-md-mlff.md) |
| NEP / GPUMD | 神经演化势（Chebyshev 描述符 + 单隐层网络，SNES 进化算法训练） | GPU；`nep` 训练，`gpumd` 跑 MD | 热输运、大体系长时间 MD，速度快 | [GPUMD 文档](https://gpumd.org/)、[Scripts/README.md](../../Scripts/README.md) |
| DeePMD-kit / DP-GEN | 深度势（Deep Potential，DPA 系列预训练模型） | TensorFlow/PyTorch 等多后端；LAMMPS 接口 | 并发学习（DP-GEN）、生态成熟 | [DeePMD-kit 文档](https://docs.deepmodeling.com/projects/deepmd/) |
| MACE | 高阶等变消息传递网络；MACE-MP、MACE-OFF 等基础模型 | PyTorch；ASE、LAMMPS、OpenMM 接口 | 开箱即用的通用势、微调（fine-tuning） | [MACE 文档](https://mace-docs.readthedocs.io/) |

其他常见名字：NequIP/Allegro、CHGNet、M3GNet（matgl）、SevenNet、GRACE（VASP 6.6.0 起可通过 VASPml 以仅预测模式实验性调用，需 TensorFlow/cppflow，见 VASP Changelog）等。通用势的横向比较可参考 [Matbench Discovery](https://matbench-discovery.materialsproject.org/) 排行榜，最新的项目汇总见 [best-of-atomistic-machine-learning](https://github.com/JuDFTteam/best-of-atomistic-machine-learning)。

> [!CAUTION]
> 通用 MLIP（基础模型）只是“很好的初猜”，不是 DFT 的替代品：它在训练集未覆盖的化学环境（特殊表面、缺陷、电荷态、高压、磁性）中可能给出看似合理但错误的结果。**正确用法是 MLIP 预筛选/预优化，DFT 精修和报告最终数值**，论文中说明 MLIP 的用途和模型版本。

### 5.2 用法一：MACE-MP 预优化结构，再交给 DFT

手搭或替换得到的结构往往离平衡位置很远，直接上 DFT 会花很多离子步，甚至 SCF 发散。用 MACE-MP 先优化几十秒到几分钟，可以显著减少 DFT 的离子步数。脚本 [mace_ase_preopt.py](../examples/07-mlip-nep-mace/mace_ase_preopt.py) 的核心逻辑：

```python
from ase.io import read, write
from ase.filters import FrechetCellFilter          # ASE >= 3.23
from ase.optimize import FIRE
from mace.calculators import mace_mp

atoms = read("POSCAR_guess")
atoms.calc = mace_mp(model="medium", device="cpu", default_dtype="float64")
FIRE(FrechetCellFilter(atoms), trajectory="preopt.traj").run(fmax=0.05, steps=500)
atoms.calc = None
write("POSCAR_preopt", atoms, format="vasp", direct=True, sort=False)
```

```bash
conda activate mlip                      # pip install ase mace-torch（单独环境）
cd Common/examples/07-mlip-nep-mace
python mace_ase_preopt.py POSCAR_guess --cell --fmax 0.05 -o POSCAR_preopt
python mace_ase_preopt.py slab.vasp --d3 --device cuda          # 吸附/层状体系加 D3
python mace_ase_preopt.py mol.xyz --calc mace-off -o mol_preopt.xyz
python mace_ase_preopt.py POSCAR --calc deepmd --model frozen_model.pth
python mace_ase_preopt.py POSCAR --calc nep --model nep.txt      # 需要 calorine
```

| 选项 | 推荐值 | 说明 |
|---|---|---|
| `--calc` | mace-mp | 无机材料用 MACE-MP；有机分子用 MACE-OFF；有专用模型时用 `mace`/`deepmd`/`nep` |
| `--model` | medium | MACE-MP-0 的 small/medium/large、`medium-mpa-0` 等模型名或模型文件路径；脚本显式默认 `medium`（MACE-MP-0），而 mace-torch ≥ 0.3.10 中 `mace_mp()` 不传参数时默认是 `medium-mpa-0` |
| `--cell` | 体相需要 | 同时优化晶胞（FrechetCellFilter）；slab/分子不要加 |
| `--d3` | 层状、吸附体系 | MACE-MP 训练数据是 PBE(+U)，不含色散，需要时叠加 D3 |
| `--fmax` | 0.05 eV/Å | 预优化无需太紧，DFT 还会继续优化 |
| `--dtype` | float64 | 优化用双精度更稳；float32 更快 |

检查：输出的最短原子间距是否合理（原子“塌缩”说明模型外推失败）、体积变化是否离谱；然后**以 POSCAR_preopt 为起点做正常的 DFT 优化**（[VASP 第4章](../../VASP/tutorials/04-relax-eos.md)），并比较两者的能量和结构差异，积累对模型可靠性的判断。

> [!NOTE]
> MACE-MP 是在 MPtrj（Materials Project 弛豫轨迹，PBE/PBE+U）上训练的，其能量是原始 VASP 能量，不能与 Materials Project 数据库中经过 MP2020 兼容性修正的能量直接比较（MACE 官方说明）。此外 MACE-OFF、MACE-OMAT 等部分模型采用 ASL 学术许可，商用前需确认。

### 5.3 用法二：NEB 初猜

NEB 的初始路径越接近最小能量路径，DFT 收敛越快。用 IDPP 插值（ASE）加 MLIP CI-NEB 先松弛一遍路径，再交给 VASP+VTST（[VASP 第11章](../../VASP/tutorials/11-neb-dimer.md)）：

```bash
python mlip_neb_guess.py IS/CONTCAR FS/CONTCAR --nimages 5 --outdir neb_init
ls neb_init                       # 00 01 02 03 04 05 06
cp IS/OUTCAR neb_init/00/ ; cp FS/OUTCAR neb_init/06/
```

| 选项 | 推荐值 | 说明 |
|---|---|---|
| `--nimages` | 3–7 | 中间像数，对应 VASP 的 `IMAGES`；总核数需为其整数倍 |
| `--method` | improvedtangent | ASE 3.28 起 NEB 默认切线方法改为 improvedtangent，脚本显式写出以便复现 |
| `--fmax` | 0.1 eV/Å | 初猜即可，不必收紧 |
| `--no-neb` | — | 只输出 IDPP 路径（不需要安装 MLIP） |

初末态必须原子顺序一致，并且已经用 DFT 优化好。MLIP 给出的能垒只用于判断路径是否合理（例如是否有不合理的中间极小），最终能垒以 DFT CI-NEB 为准。CP2K 和 OpenMX 的 NEB 见 [CP2K 第13章](../../CP2K/tutorials/13-neb-dimer.md)、[OpenMX 第6章](../../OpenMX/tutorials/06-optimization-eos-neb.md)，可用 [convert_any.py](../scripts/convert_any.py) 把各像的 POSCAR 转为对应格式（在 `07-mlip-nep-mace/` 目录下运行 `python ../../scripts/convert_any.py -h` 查看用法）。

### 5.4 用法三：用 DFT 数据训练专用势

```text
 DFT 数据            转格式                 训练                     验证 / 主动学习
 AIMD 抽帧  ──┐
 微扰/应变   ──┼─> dft_to_extxyz.py ─> train.xyz / test.xyz ─> nep / mace_run_train / dp train
 单点能      ──┘   (dpdata 也可)                                   │
     ^                                                              v
     └──── 不确定度高或误差大的构型回算 DFT <──── MLIP MD 采样（GPUMD / LAMMPS / ASE）
```

**数据准备**（[dft_to_extxyz.py](../examples/07-mlip-nep-mace/dft_to_extxyz.py)）：

```bash
python dft_to_extxyz.py run_*/vasprun.xml --format nep  -o train.xyz
python dft_to_extxyz.py aimd/vasprun.xml --every 10 --format mace -o mace_train.xyz
python -c "from ase.io import read; print(len(read('train.xyz', ':')))"
```

NEP 读取的 extxyz 必须包含 `lattice`、`energy`（整个晶胞，eV）和 `properties=species:S:1:pos:R:3:force:R:3`（或 `forces`），位力 `virial` 可选，单位 eV（正值表示压缩），详见 [train.xyz 格式说明](https://gpumd.org/nep/input_files/train_test_xyz.html)。MACE 默认读 `energy`/`forces` 键；为避免与 ASE 计算器结果混淆，示例统一写成 `REF_energy`、`REF_forces`、`REF_stress`，并在 YAML 中用 `energy_key` 等显式指定。

**NEP 训练输入**（[nep.in](../examples/07-mlip-nep-mace/nep.in)）：

```text
type          2 Mo S
version       4
cutoff        6 4
n_max         6 6
basis_size    6 6
l_max         4 1
neuron        30
lambda_e      1.0
lambda_f      1.0
lambda_v      0.1
batch         1000
population    50
generation    20000
```

| 关键词 | 推荐值 | 说明 |
|---|---|---|
| `type` | 元素数 + 元素符号 | 唯一必填项，顺序会写入 nep.txt |
| `version` | 4 | 当前只支持 NEP4 |
| `cutoff` | 默认 8 4（Å） | 径向/角向截断，要求 3 ≤ 角向 ≤ 径向 ≤ 10 Å；较小的径向截断可提速 |
| `n_max`、`basis_size` | 6 6、6 6（默认） | 描述符规模 |
| `l_max` | 4 1（默认） | 三体角向阶数与四体描述符开关 |
| `neuron` | 30（默认） | 隐藏层神经元数 |
| `lambda_e/f/v` | 1.0 / 1.0 / 0.1（默认） | 能量、力、位力权重；没有可靠应力数据时 `lambda_v 0` |
| `batch` | 1000（默认） | 官方建议多数情况下 100–1000 |
| `generation` | 默认 100000 | 简单体系 1e4–1e5，复杂体系可到 1e6；示例取小值只为检查流程 |
| `fine_tune` | — | 可从官方 NEP89 基础模型微调，架构参数必须与基础模型一致（见 GPUMD 文档） |

```bash
sbatch --export=ALL,TASK=nep job_mlip_gpu.sh
tail -3 loss.out                          # 各列：generation total L1 L2 E_train F_train V_train E_test ...
python ../../../Scripts/nep/plot_nep_results.py    # 先复制到工作目录、按说明修改文件头变量
```

作者整理的 NEP/GPUMD 全流程脚本（AIMD 抽帧、微扰/应变/rattle 数据集、离群结构筛选、最远点采样、声子对比、NEMD/HNEMD/EMD 热导率分析）在 [Scripts/README.md](../../Scripts/README.md) 中逐个说明。

**MACE 训练**（[mace_config.yaml](../examples/07-mlip-nep-mace/mace_config.yaml)）：

```bash
mace_run_train --config mace_config.yaml            # 命令行参数会覆盖 YAML 中的同名项
mace_run_train --config mace_config.yaml --foundation_model small --E0s estimated   # 从 MACE-MP 微调
mace_eval_configs --configs test.xyz --model MoS2_mace.model --output test_pred.xyz
```

| 选项 | 推荐值 | 说明 |
|---|---|---|
| `r_max` | 5.0 Å | 截断半径；感受野 = 层数 × r_max |
| `num_channels` / `max_L` | 128 / 1 | 对应官方推荐的 `128x0e + 128x1o` |
| `E0s` | 孤立原子能量（`config_type=IsolatedAtom`）优先；`average` 为快速替代 | 官方建议优先把孤立原子放进训练集；`average` 用最小二乘回归估计原子能量；微调基础模型时官方示例用 `estimated` |
| `batch_size` / `max_num_epochs` | 按数据量调整 | 官方经验：总梯度更新次数约 20 万 |
| `stage_two` | yes | 训练后段（官方说明为最后约 25% 轮次，或从 `start_stage_two` 指定的轮次起）加大能量权重（旧名 swa/start_swa） |
| `default_dtype` | float64（默认） | float32 训练约快 2 倍 |

**DeePMD-kit**：示例 [deepmd_input.json](../examples/07-mlip-nep-mace/deepmd_input.json) 是 `se_e2_a` 描述符的最小能量模型，`sel` 设为 `"auto"` 由程序统计近邻数。训练与冻结：`dp --pt train deepmd_input.json`、`dp --pt freeze`（PyTorch 后端；TensorFlow 后端为 `dp train`/`dp freeze`，模型文件为 `.pb`）。数据需要用 dpdata 转为 `deepmd/npy` 格式，参见 [04-output-parsing/dpdata_convert.py](../examples/04-output-parsing/dpdata_convert.py)。并发学习流程（探索–筛选–标注循环）可用 [DP-GEN](https://github.com/deepmodeling/dpgen)。

### 5.5 VASP 内置 MLFF（简述）

VASP 的 MLFF 在 AIMD 过程中边算边学：贝叶斯误差估计超过阈值时调用 DFT 并把构型加入训练集，否则直接用力场预测。关键流程（详见 [VASP 第15章](../../VASP/tutorials/15-md-mlff.md)）：

```bash
# 1) on-the-fly training              INCAR: ML_LMLFF = .TRUE. ; ML_MODE = train
# 2) continue training                cp ML_ABN ML_AB   (never reuse the old ML_AB)
# 3) refit for fast prediction        ML_MODE = refit   -> ML_FFN
# 4) production                       cp ML_FFN ML_FF ; ML_MODE = run
grep "ERR" ML_LOGFILE | tail          # training errors per step
```

它的优点是与 VASP 设置天然一致、无需另装软件；缺点是只能在 VASP 中运行。需要在 LAMMPS/GPUMD 中跑更大体系时，可以把 VASP AIMD 数据导出后训练 NEP/MACE/DeePMD。

### 5.6 如何验证一个机器学习势

| 检查 | 方法 | 合格的表现 |
|---|---|---|
| 测试集误差 | 独立的 test.xyz（不同 AIMD 片段、不同温度），看能量 RMSE（meV/atom）、力 RMSE（meV/Å） | 测试集误差与训练集接近，没有明显过拟合 |
| 对角图 | MLIP vs DFT 的能量/力散点 | 点集中在对角线，无离群簇 |
| 物理性质 | 晶格常数、弹性常数、声子谱（与 phonopy+DFT 对比）、RDF、扩散系数 | 与 DFT 一致到你的研究需要的精度 |
| 稳定性 | 在目标温度/压力下跑足够长的 MD | 不出现原子重叠、能量漂移或结构崩溃 |
| 外推检测 | 委员会模型（多个模型的偏差）、VASP MLFF 的贝叶斯误差、描述符空间距离 | 生产 MD 中的构型都在训练分布内 |

误差阈值没有通用标准，应与所研究的物理量对应：例如比较相对稳定性时，能量误差要明显小于候选结构的能量差。

---

## 6. 资源汇总

### 6.1 教材

| 书 | 完整书目信息 | 适合 |
|---|---|---|
| Sholl & Steckel | D. S. Sholl, J. A. Steckel, *Density Functional Theory: A Practical Introduction*, 2nd ed., John Wiley & Sons, Hoboken, NJ, 2023, ISBN 978-1-119-84086-2。第 1 版：2009, ISBN 978-0-470-37317-0, DOI [10.1002/9780470447710](https://doi.org/10.1002/9780470447710) | **首选入门书**，偏实战，每章有习题 |
| Martin | R. M. Martin, *Electronic Structure: Basic Theory and Practical Methods*, 2nd ed., Cambridge University Press, Cambridge, 2020, ISBN 978-1-108-42990-0, DOI [10.1017/9781108555586](https://doi.org/10.1017/9781108555586) | 固体电子结构理论的标准参考书 |
| Kohanoff | J. Kohanoff, *Electronic Structure Calculations for Solids and Molecules: Theory and Computational Methods*, Cambridge University Press, 2006, ISBN 978-0-521-81591-8, DOI [10.1017/CBO9780511755613](https://doi.org/10.1017/CBO9780511755613) | 赝势、基组、MD 的方法细节 |
| Giustino | F. Giustino, *Materials Modelling Using Density Functional Theory: Properties and Predictions*, Oxford University Press, Oxford, 2014, ISBN 978-0-19-966244-9 | 面向本科生的材料 DFT 教材 |
| Koch & Holthausen | W. Koch, M. C. Holthausen, *A Chemist's Guide to Density Functional Theory*, 2nd ed., Wiley-VCH, Weinheim, 2001, ISBN 978-3-527-30372-4, DOI [10.1002/3527600043](https://doi.org/10.1002/3527600043) | 化学背景读者的 DFT 入门 |
| Parr & Yang | R. G. Parr, W. Yang, *Density-Functional Theory of Atoms and Molecules*, Oxford University Press, New York, 1989（平装本 1994, ISBN 978-0-19-509276-9） | DFT 理论经典 |
| Szabo & Ostlund | A. Szabo, N. S. Ostlund, *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*, Dover Publications, Mineola, NY, 1996, ISBN 978-0-486-69186-2 | Hartree–Fock 与基组的基础 |
| Jensen | F. Jensen, *Introduction to Computational Chemistry*, 3rd ed., John Wiley & Sons, Chichester, 2017, ISBN 978-1-118-82599-0 | 计算化学方法全景 |
| Cramer | C. J. Cramer, *Essentials of Computational Chemistry: Theories and Models*, 2nd ed., John Wiley & Sons, Chichester, 2004, ISBN 978-0-470-09182-1 | 同上，溶剂模型讲得好 |
| Foresman & Frisch | J. B. Foresman, Æ. Frisch, *Exploring Chemistry with Electronic Structure Methods*, 3rd ed., Gaussian, Inc., Wallingford, CT, 2015, ISBN 978-1-935522-03-4 | Gaussian 官方实战书 |
| Marx & Hutter | D. Marx, J. Hutter, *Ab Initio Molecular Dynamics: Basic Theory and Advanced Methods*, Cambridge University Press, 2009, ISBN 978-0-521-89863-8, DOI [10.1017/CBO9780511609633](https://doi.org/10.1017/CBO9780511609633) | AIMD 理论（CP2K 作者之一） |
| Frenkel & Smit | D. Frenkel, B. Smit, *Understanding Molecular Simulation: From Algorithms to Applications*, 3rd ed., Academic Press (Elsevier), London, 2023, ISBN 978-0-323-90292-2 | MD/MC 与统计力学，做 MLIP-MD 必读 |

### 6.2 在线书、课程与工作坊

| 资源 | 链接 | 说明 |
|---|---|---|
| dft-book（J. Kitchin） | [jkitchin/dft-book](https://github.com/jkitchin/dft-book)，在线版 <https://kitchingroup.cheme.cmu.edu/dft-book/> | ASE + VASP 实战：分子、体相、表面、原子热力学 |
| Sholl & Steckel 习题中文 VASP 解答 | [cndaqiang/DFT-EXERCISES](https://github.com/cndaqiang/DFT-EXERCISES) | 第 2–9 章习题；无许可证，只供参考 |
| Materials Project workshop | <https://workshop.materialsproject.org/>，源码 [materialsproject/workshop](https://github.com/materialsproject/workshop) | pymatgen、MP API、atomate2、机器学习 |
| pymatgen 视频教程 | [computron/pymatgen_tutorials](https://github.com/computron/pymatgen_tutorials) | 结构、VASP I/O、MP API、对称性（视频 + notebook） |
| VASP 官方教程 | <https://www.vasp.at/tutorials/latest/> | 分子、体相、表面、MD、MLFF 等官方练习 |
| CP2K HOWTOs | <https://www.cp2k.org/howto> | 官方入门与专题，含 [截断收敛](https://www.cp2k.org/howto:converging_cutoff) |
| OpenMX workshop | <https://www.openmx-square.org/workshop/workshop.html> | 历届 OpenMX 讲习班讲义 |
| GPUMD/NEP 教程 | [brucefan1983/GPUMD-Tutorials](https://github.com/brucefan1983/GPUMD-Tutorials) | NEP 训练与 GPUMD 模拟的官方示例 |
| MACE 教程 | <https://mace-docs.readthedocs.io/>（含 Colab 训练/微调教程） | 训练、微调、基础模型使用 |
| DeePMD-kit 文档 | <https://docs.deepmodeling.com/projects/deepmd/> | 快速上手与 DPA 预训练模型 |
| CompMatBook | [stanfordbshan/CompMatBook](https://github.com/stanfordbshan/CompMatBook) | 开源计算材料学教材 |
| 中文 DFT+VASP 20 课 | [haiw201/dft-tutorial](https://github.com/haiw201/dft-tutorial) | 可与本仓库对照查漏补缺（无许可证，只链接） |

### 6.3 awesome 列表

GitHub 上没有高星的 “awesome-dft” 或 “awesome-computational-chemistry”，下面五个是最接近的替代：

- [tilde-lab/awesome-materials-informatics](https://github.com/tilde-lab/awesome-materials-informatics) — 材料信息学：AiiDA、ASE、atomate、pymatgen、数据库；
- [lmmentel/awesome-python-chemistry](https://github.com/lmmentel/awesome-python-chemistry) — 化学相关 Python 包（模拟、可视化、机器学习）；
- [JuDFTteam/best-of-atomistic-machine-learning](https://github.com/JuDFTteam/best-of-atomistic-machine-learning) — 原子尺度机器学习项目按热度排名（MLIP、描述符、数据集）；
- [hsiaoyi0504/awesome-cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics) — 化学信息学（RDKit、Open Babel 等），分子建模时有用；
- [dstdev/awesome-hpc](https://github.com/dstdev/awesome-hpc) — HPC 软件与工具（调度器、MPI、性能分析）。

本仓库自己整理的 230 多个相关仓库按软件和类别分类在 [Resources/github-repos.md](../../Resources/github-repos.md)。

### 6.4 官方论坛与社区

| 社区 | 链接 | 适合问什么 |
|---|---|---|
| VASP Forum | <https://www.vasp.at/forum/> | VASP 输入、报错、编译（需注册，开发者会回复） |
| matsci.org | <https://matsci.org/> | pymatgen、atomate2、Materials Project、ASE 等材料软件的讨论区 |
| CP2K Google Group | <https://groups.google.com/g/cp2k> | CP2K 用户与开发者邮件组 |
| OpenMX Forum | <https://www.openmx-square.org/forum/patio.cgi> | OpenMX 官方论坛 |
| Computational Chemistry List（CCL） | <https://ccl.net/> | 计算化学综合邮件列表，Gaussian 相关问题也常见 |
| 计算化学公社 | <http://bbs.keinsci.com/> | 中文最活跃的计算化学/材料论坛，Gaussian、VASP、CP2K 均有版块 |
| sobereva 博客 | <http://sobereva.com/> | Multiwfn 作者的中文博客，量子化学方法与 Gaussian 实践文章 |
| VASPKIT | <https://vaspkit.com/>；旧版中文手册 [tamaswells/VASPKIT_manual](https://github.com/tamaswells/VASPKIT_manual) | VASP 前后处理；中文手册版本较旧，以官网 1.5 文档为准 |
| GPUMD / MACE / DeePMD | 各自 GitHub 仓库的 Issues 与 Discussions | 训练报错、接口问题 |

> [!TIP]
> 提问前先搜索论坛历史。提问时给出：软件与版本、完整输入（INCAR/KPOINTS/POSCAR 或 .gjf/.inp/.dat）、报错前后几十行原文、作业脚本、已经尝试过的办法，最好附一个能复现问题的最小算例。不要只贴一张截图问“为什么报错”。

### 6.5 本仓库内的资源入口

- 四个软件板块：[VASP](../../VASP/README.md)、[Gaussian](../../Gaussian/README.md)、[CP2K](../../CP2K/README.md)、[OpenMX](../../OpenMX/README.md)；
- NEP/GPUMD 机器学习势与数据集脚本：[Scripts/README.md](../../Scripts/README.md)；
- GitHub 仓库与官方文档精选：[Resources/github-repos.md](../../Resources/github-repos.md)；
- 本板块目录：[Common/README.md](../README.md)；
- 本章的书签清单：[07-learning-resources/resources.md](../examples/07-learning-resources/resources.md)。

---

## 7. 引用规范与许可证注意事项

### 7.1 该引用什么

发表文章时，除了引用所用软件的主论文，还要引用**你实际用到的方法和工具**：泛函、色散校正、赝势/基组、后处理工具、机器学习势和基础模型。许多程序会在输出或文档中说明引用方式（Gaussian 16 的引用格式见其官网，VASP、CP2K 在官网列出推荐引用）。

| 对象 | 推荐引用（示例） |
|---|---|
| VASP | G. Kresse, J. Furthmüller, *Phys. Rev. B* **54**, 11169 (1996), DOI [10.1103/PhysRevB.54.11169](https://doi.org/10.1103/PhysRevB.54.11169)；PAW：P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), DOI [10.1103/PhysRevB.50.17953](https://doi.org/10.1103/PhysRevB.50.17953)；G. Kresse, D. Joubert, *Phys. Rev. B* **59**, 1758 (1999), DOI [10.1103/PhysRevB.59.1758](https://doi.org/10.1103/PhysRevB.59.1758) |
| CP2K | T. D. Kühne *et al.*, *J. Chem. Phys.* **152**, 194103 (2020), DOI [10.1063/5.0007045](https://doi.org/10.1063/5.0007045) |
| OpenMX | T. Ozaki, *Phys. Rev. B* **67**, 155108 (2003), DOI [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108) |
| Gaussian 16 | 按 Gaussian 官网给出的完整引用格式（含版本号，如 Revision C.01/C.02） |
| ASE | A. H. Larsen *et al.*, *J. Phys.: Condens. Matter* **29**, 273002 (2017), DOI [10.1088/1361-648X/aa680e](https://doi.org/10.1088/1361-648X/aa680e) |
| pymatgen | S. P. Ong *et al.*, *Comput. Mater. Sci.* **68**, 314 (2013), DOI [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028) |
| NEP / GPUMD | Z. Fan *et al.*, *Phys. Rev. B* **104**, 104309 (2021)；Z. Fan *et al.*, *J. Chem. Phys.* **157**, 114801 (2022)（见参考资料） |
| MACE / MACE-MP | I. Batatia *et al.*, arXiv:2206.07697；I. Batatia *et al.*, arXiv:2401.00096（基础模型） |
| DeePMD-kit | H. Wang *et al.*, *Comput. Phys. Commun.* **228**, 178 (2018)；J. Zeng *et al.*, *J. Chem. Phys.* **159**, 054801 (2023) |

### 7.2 许可证：只链接，不复制

- **没有许可证的仓库默认保留全部权利**（如 tamaswells/VASP_script、haiw201/dft-tutorial、QijingZheng/VaspBandUnfolding）：可以阅读和链接，不能把代码或文字复制进自己的公开仓库；
- GPL 类代码（如 pyprocar、WannierTools、Z2Pack）可以使用，但修改后再发布需要遵守 GPL 条款；
- **VASP 的 POTCAR、源码和可执行文件受许可证保护**，绝不能上传到 GitHub 或分享给组外人员；Gaussian 的程序和许可证文件同理；
- 机器学习势的**模型权重也有许可证**：MACE-MP-0 为 MIT，MACE-OFF、MACE-OMAT 等为 ASL（学术软件许可）；使用和再分发前查看模型说明；
- 本仓库采用 MIT 协议，引用的外部资源一律给出链接，不复制原文。

### 7.3 可复现性

在文章“计算方法”部分或补充材料中写明：软件与版本号、泛函与色散校正、赝势/基组（含版本，如 PBE_54 POTCAR、DFT_DATA19）、截断能和 k 点网格（以及收敛测试的判据）、展宽方式、收敛阈值、MLIP 模型名称与版本。DFT 不同实现之间的一致性研究（Δ 测试、ACWF 验证集）表明，现代代码在参数收敛时可以给出高度一致的结果——前提正是参数收敛且记录完整（Lejaeghere 2016；Bosoni 2024）。

---

## 参数选择与收敛建议

| 项目 | 建议起点 | 判据 / 说明 |
|---|---|---|
| VASP ENCUT | 从 max(ENMAX) 扫到 ~2×ENMAX；变胞 ≥ 1.3×ENMAX | ΔE < 1 meV/atom；检查 `external pressure` |
| VASP k 点 | KSPACING 0.5 → 0.15 Å⁻¹ 或 n×n×n 逐级加密 | ΔE < 1 meV/atom；金属需要更密 |
| VASP SIGMA | 半导体 ISMEAR=0、0.01–0.05；金属 ISMEAR=1/2、0.1–0.2 | 每原子 \|T*S\| < 1 meV |
| CP2K | 先 CUTOFF（REL_CUTOFF 60），再 REL_CUTOFF | ΔE < 1 meV/atom + `count for grid` 分布合理 + 网格电子数积分误差小 |
| OpenMX | `scf.energycutoff` 120–400 Ry；`scf.Kgrid` 逐级加密 | ΔE < 1 meV/atom；优化/声子用更高截断 |
| Gaussian | 基组系列逐级增大；格点 UltraFine（G16 默认）/ SuperFine | 关心的相对能量稳定 |
| SCF 阈值 | VASP `EDIFF` 1E-6（测试时 1E-7）；CP2K `EPS_SCF` 1E-6–1E-7 | 比目标精度严格约两个数量级 |
| MLIP 预优化 | MACE-MP medium，`fmax` 0.05 | 之后必须 DFT 精修 |
| MLIP 训练 | 同一 DFT 设置；训练/测试集分开；孤立原子能量作参考 | 测试集误差、物理性质、长时间 MD 稳定性 |

---

## 常见问题与报错

| 报错 / 现象 | 原因 | 解决 |
|---|---|---|
| ENCUT 扫描中能量不单调或有跳变 | `PREC`、`LREAL`、`EDIFF` 在各点不同；或 SCF 未收敛 | 只改 ENCUT；`LREAL=.FALSE.`、`EDIFF=1E-7`；检查每个 OUTCAR 是否正常结束 |
| 设置了 `KSPACING` 但 k 点没变 | 目录中仍有 KPOINTS 文件，VASP 优先读取它 | 删除 KPOINTS（`kpoint_scan.py` 的 KSPACING 模式会自动删除） |
| 相邻两个 KSPACING 得到完全相同的能量 | 生成的网格相同 | 用 `kpoint_scan.py` 打印的网格去重 |
| 变胞优化后再算单点，`external pressure` 明显不为 0 | Pulay 应力 / 平面波数固定 | ENCUT ≥ 1.3×ENMAX，从 CONTCAR 重新优化 |
| CP2K 提高 CUTOFF 后能量“不再变化”，但结果与文献差异大 | REL_CUTOFF 太小，高斯都落在粗网格，`count for grid 1` 接近 0 | 再做 REL_CUTOFF 扫描 |
| CP2K 扫描某点出现 `SCF run NOT converged` | 展宽/混合参数不合适，或截断过低 | 修正 SCF 设置后重算该点，未收敛的点不能用于判断 |
| Gaussian `Atomic number out of range for ... basis set.` | 该基组未定义此元素（如 6-31G 用于重元素） | 换用 def2 系列或用 `Gen`/`GenECP` 自定义基组（[Gaussian 第3章](../../Gaussian/tutorials/03-method-basis.md)） |
| Gaussian 不同格点能量差很大（M06 类泛函） | meta-GGA 对积分格点敏感 | 使用 `Int=SuperFine` 并全项目统一 |
| `ModuleNotFoundError: No module named 'mace'` | 没装或装错了包，PyPI 包名是 `mace-torch` | `pip install mace-torch`（建议单独的 conda 环境） |
| `RuntimeError: Please install torch-dftd to use dispersion corrections (see https://github.com/pfnet-research/torch-dftd)` | 使用 `--d3`（`mace_mp(dispersion=True)`）需要 torch-dftd | `pip install torch-dftd` |
| `ModuleNotFoundError: No module named 'ase.filters'` 或无法导入 `FrechetCellFilter` | ASE 版本低于 3.23 | 升级 ASE；脚本会自动回退到旧的 `ExpCellFilter` |
| 首次调用 `mace_mp()` 很慢或在计算节点上失败 | 模型需要从 GitHub 下载到 `~/.cache/mace`，计算节点通常不能联网 | 在登录节点先运行一次下载，或用 `--model` 指向本地模型文件 |
| MLIP 预优化后原子重叠、能量异常低 | 结构超出模型训练分布 | 丢弃该结果；换更大的模型、加 D3，或直接用 DFT |
| `The model in nep.in is inconsistent with nep.txt:` | 目录中残留的 nep.txt/nep.restart 与当前 nep.in 的架构参数不同 | 续训时保持参数一致；重新训练前删除 `nep.restart`（只残留 nep.txt 时仅警告） |
| NEP 提示参考能量低于 −100 eV/atom | NEP 用单精度训练，绝对能量太大损失精度 | 减去孤立原子能量后再训练（全电子 Gaussian 能量尤其需要） |
| MLIP 测试误差很小，但 MD 中结构崩溃 | 训练集缺少高温/大形变构型 | 主动学习：从 MD 轨迹抽帧回算 DFT 并加入训练集 |
| VASP MLFF 续训后力场反而变差 | 续训时没有执行 `cp ML_ABN ML_AB`，丢失了上一轮数据 | 每轮训练后先复制 ML_ABN 为 ML_AB |
| `sbatch: error: Batch script contains DOS line breaks (\r\n)` | 脚本在 Windows 下编辑 | `dos2unix *.sh` |

---

## 练习

1. 用 [07-convergence-templates](../examples/07-convergence-templates/README.md) 对 Si 做 ENCUT 和 k 点（KSPACING 0.5–0.15）测试，给出满足 1 meV/atom 的最小参数；再对 ENCUT = ENMAX 与 1.3×ENMAX 各做一次 `ISIF = 3` 优化，比较晶格常数和单点的 `external pressure`。
2. 对金属 Al（或 Cu）比较 `ISMEAR = 0`、1 和 −5 在不同 SIGMA 下的能量和每原子 \|T*S\|，用 `vasp_tsigma_check.py` 汇总，说明为什么金属优化推荐 Methfessel–Paxton。
3. CP2K：用模板分别扫描 CUTOFF（REL_CUTOFF 固定 40 和 60 两组）和 REL_CUTOFF，记录 `count for grid`，解释“只升 CUTOFF 不升 REL_CUTOFF”时最细网格上的高斯个数如何变化。
4. Gaussian：用 `gaussian_basis_scan` 计算 HCN→HNC 异构化能随基组的变化，画出 def2 与 cc-pVnZ 两条曲线；比较 FineGrid 与 UltraFine 的差别。
5. 用 MACE-MP 预优化一个你手搭的结构（如 MoS2 中一个 S 空位的 4×4 超胞），再做 DFT 优化，比较两者的最大位移、所需离子步数和能量变化。
6. 从 [VASP 第15章](../../VASP/tutorials/15-md-mlff.md) 的 AIMD 算例中每 10 步抽一帧，用 `dft_to_extxyz.py` 生成 NEP 训练集，按 9:1 分为训练集和测试集，训练一个小 NEP 模型，用 [Scripts/nep/plot_nep_results.py](../../Scripts/nep/plot_nep_results.py) 画对角图。
7. 选一篇与你课题相关的文章，按第 7.3 节列出它在“计算方法”中缺失的信息，并写出你自己的计算方法段落。

---

## 参考资料

官方文档：

- VASP Wiki：[ENCUT](https://www.vasp.at/wiki/index.php/ENCUT)、[KSPACING](https://www.vasp.at/wiki/index.php/KSPACING)、[ISMEAR](https://www.vasp.at/wiki/index.php/ISMEAR)、[Pulay stress](https://www.vasp.at/wiki/index.php/Pulay_stress)、[Volume relaxation](https://www.vasp.at/wiki/index.php/Volume_relaxation)、[ML_MODE](https://www.vasp.at/wiki/index.php/ML_MODE)、[Best practices for machine-learned force fields](https://www.vasp.at/wiki/index.php/Best_practices_for_machine-learned_force_fields)
- CP2K：[Multigrid 与截断说明](https://manual.cp2k.org/trunk/methods/dft/cutoff.html)、[Converging the CUTOFF HOWTO](https://www.cp2k.org/howto:converging_cutoff)
- OpenMX 手册与论坛：<https://www.openmx-square.org/>
- Gaussian 16 Integral 关键词：<https://gaussian.com/integral/>
- GPUMD/NEP：[nep.in](https://gpumd.org/nep/input_files/nep_in.html)、[train.xyz 格式](https://gpumd.org/nep/input_files/train_test_xyz.html)、[fine_tune](https://gpumd.org/nep/input_parameters/fine_tune.html)
- MACE：<https://mace-docs.readthedocs.io/>；基础模型发布：[ACEsuit/mace-foundations](https://github.com/ACEsuit/mace-foundations)
- DeePMD-kit：<https://docs.deepmodeling.com/projects/deepmd/>
- calorine（NEP 的 ASE 接口）：<https://calorine.materialsmodeling.org/>
- ASE NEB（3.28 起默认 improvedtangent）：<https://docs.ase-lib.org/releasenotes.html>

相关 GitHub 仓库（只链接，不复制）：

- [ACEsuit/mace](https://github.com/ACEsuit/mace) — MACE 与 MACE-MP，可作为 ASE 计算器做预优化与 NEB 初猜。
- [brucefan1983/GPUMD](https://github.com/brucefan1983/GPUMD) — GPUMD 与 NEP 训练程序。
- [deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit)、[deepmodeling/dpgen](https://github.com/deepmodeling/dpgen)、[deepmodeling/dpdata](https://github.com/deepmodeling/dpdata) — 深度势、并发学习与数据格式转换。
- [JuDFTteam/best-of-atomistic-machine-learning](https://github.com/JuDFTteam/best-of-atomistic-machine-learning) — 原子尺度机器学习项目排行。
- [kavanase/vaspup2.0](https://github.com/kavanase/vaspup2.0) — ENCUT/k 点收敛测试脚本（MIT）。
- [jkitchin/dft-book](https://github.com/jkitchin/dft-book) — ASE + VASP 实战书。
- [materialsproject/workshop](https://github.com/materialsproject/workshop)、[computron/pymatgen_tutorials](https://github.com/computron/pymatgen_tutorials) — 课程与视频教程。
- [cndaqiang/DFT-EXERCISES](https://github.com/cndaqiang/DFT-EXERCISES)、[haiw201/dft-tutorial](https://github.com/haiw201/dft-tutorial)、[northword/wiki](https://github.com/northword/wiki) — 中文学习资料（无许可证或 CC 协议，只链接）。
- [tilde-lab/awesome-materials-informatics](https://github.com/tilde-lab/awesome-materials-informatics)、[lmmentel/awesome-python-chemistry](https://github.com/lmmentel/awesome-python-chemistry)、[hsiaoyi0504/awesome-cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics)、[dstdev/awesome-hpc](https://github.com/dstdev/awesome-hpc) — awesome 列表。
- [stanfordbshan/CompMatBook](https://github.com/stanfordbshan/CompMatBook) — 开源计算材料学教材。

关键文献：

- P. Hohenberg, W. Kohn, Inhomogeneous electron gas, *Phys. Rev.* **136**, B864 (1964). DOI: [10.1103/PhysRev.136.B864](https://doi.org/10.1103/PhysRev.136.B864)
- W. Kohn, L. J. Sham, Self-consistent equations including exchange and correlation effects, *Phys. Rev.* **140**, A1133 (1965). DOI: [10.1103/PhysRev.140.A1133](https://doi.org/10.1103/PhysRev.140.A1133)
- J. P. Perdew, K. Burke, M. Ernzerhof, Generalized gradient approximation made simple, *Phys. Rev. Lett.* **77**, 3865 (1996). DOI: [10.1103/PhysRevLett.77.3865](https://doi.org/10.1103/PhysRevLett.77.3865)
- J. P. Perdew, K. Schmidt, Jacob's ladder of density functional approximations for the exchange-correlation energy, *AIP Conf. Proc.* **577**, 1 (2001). DOI: [10.1063/1.1390175](https://doi.org/10.1063/1.1390175)
- H. J. Monkhorst, J. D. Pack, Special points for Brillouin-zone integrations, *Phys. Rev. B* **13**, 5188 (1976). DOI: [10.1103/PhysRevB.13.5188](https://doi.org/10.1103/PhysRevB.13.5188)
- G. P. Francis, M. C. Payne, Finite basis set corrections to total energy pseudopotential calculations, *J. Phys.: Condens. Matter* **2**, 4395 (1990). DOI: [10.1088/0953-8984/2/19/007](https://doi.org/10.1088/0953-8984/2/19/007)
- K. Lejaeghere *et al.*, Reproducibility in density functional theory calculations of solids, *Science* **351**, aad3000 (2016). DOI: [10.1126/science.aad3000](https://doi.org/10.1126/science.aad3000)
- E. Bosoni *et al.*, How to verify the precision of density-functional-theory implementations via reproducible and universal workflows, *Nat. Rev. Phys.* **6**, 45 (2024). DOI: [10.1038/s42254-023-00655-3](https://doi.org/10.1038/s42254-023-00655-3)
- J. Behler, M. Parrinello, Generalized neural-network representation of high-dimensional potential-energy surfaces, *Phys. Rev. Lett.* **98**, 146401 (2007). DOI: [10.1103/PhysRevLett.98.146401](https://doi.org/10.1103/PhysRevLett.98.146401)
- A. P. Bartók, M. C. Payne, R. Kondor, G. Csányi, Gaussian approximation potentials, *Phys. Rev. Lett.* **104**, 136403 (2010). DOI: [10.1103/PhysRevLett.104.136403](https://doi.org/10.1103/PhysRevLett.104.136403)
- R. Jinnouchi, J. Lahnsteiner, F. Karsai, G. Kresse, M. Bokdam, Phase transitions of hybrid perovskites simulated by machine-learning force fields trained on the fly with Bayesian inference, *Phys. Rev. Lett.* **122**, 225701 (2019). DOI: [10.1103/PhysRevLett.122.225701](https://doi.org/10.1103/PhysRevLett.122.225701)
- R. Jinnouchi, F. Karsai, G. Kresse, On-the-fly machine learning force field generation: Application to melting points, *Phys. Rev. B* **100**, 014105 (2019). DOI: [10.1103/PhysRevB.100.014105](https://doi.org/10.1103/PhysRevB.100.014105)
- Z. Fan *et al.*, Neuroevolution machine learning potentials: Combining high accuracy and low cost in atomistic simulations and application to heat transport, *Phys. Rev. B* **104**, 104309 (2021). DOI: [10.1103/PhysRevB.104.104309](https://doi.org/10.1103/PhysRevB.104.104309)
- Z. Fan *et al.*, GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations, *J. Chem. Phys.* **157**, 114801 (2022). DOI: [10.1063/5.0106617](https://doi.org/10.1063/5.0106617)
- K. Song *et al.*, General-purpose machine-learned potential for 16 elemental metals and their alloys, *Nat. Commun.* **15**, 10208 (2024)（NEP4）. DOI: [10.1038/s41467-024-54554-x](https://doi.org/10.1038/s41467-024-54554-x)
- H. Wang, L. Zhang, J. Han, W. E, DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics, *Comput. Phys. Commun.* **228**, 178 (2018). DOI: [10.1016/j.cpc.2018.03.016](https://doi.org/10.1016/j.cpc.2018.03.016)
- J. Zeng *et al.*, DeePMD-kit v2: A software package for deep potential models, *J. Chem. Phys.* **159**, 054801 (2023). DOI: [10.1063/5.0155600](https://doi.org/10.1063/5.0155600)
- Y. Zhang *et al.*, DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models, *Comput. Phys. Commun.* **253**, 107206 (2020). DOI: [10.1016/j.cpc.2020.107206](https://doi.org/10.1016/j.cpc.2020.107206)
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, G. Csányi, MACE: Higher order equivariant message passing neural networks for fast and accurate force fields, *NeurIPS* 2022, arXiv: [2206.07697](https://arxiv.org/abs/2206.07697)
- I. Batatia *et al.*, A foundation model for atomistic materials chemistry, arXiv: [2401.00096](https://arxiv.org/abs/2401.00096)（MACE-MP-0）

[← 上一章](06-workflow-automation.md) | [返回目录](../README.md)
