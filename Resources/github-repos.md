# DFT 相关 GitHub 仓库精选（Awesome DFT GitHub Repositories）

本页汇总与 **VASP、Gaussian、CP2K、OpenMX** 以及通用第一性原理（DFT）计算相关的优秀 GitHub 仓库：教程、课程、脚本合集、后处理工具与工作流框架，便于按需查找。清单由检索 GitHub 与官方文档整理而来，Star 数与最近更新时间统计于 **2026-11-10**。

> [!IMPORTANT]
> 本仓库只**链接**这些项目，不复制它们的代码或文本。使用前请阅读各项目自己的许可证；标注“未声明”的仓库默认保留所有权利，不能直接转载或再分发。

**目录**：[VASP](#vasp) · [Gaussian](#gaussian) · [CP2K](#cp2k) · [OpenMX](#openmx) · [通用工具与学习资源](#general) · [官方文档](#official-docs)

<a id="vasp"></a>

## VASP

### 教程 / 书籍

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [jkitchin/dft-book](https://github.com/jkitchin/dft-book) | 470 | other | 2026-10-28 | 经典英文 VASP+ASE 实战教程书，按分子、体相、表面、原子热力学、高级电子结构分章，适合系统学习用 ASE 驱动 VASP。 |
| [cndaqiang/DFT-EXERCISES](https://github.com/cndaqiang/DFT-EXERCISES) | 69 | 未声明 | 2021-08-27 | Sholl & Steckel《DFT: A Practical Introduction》第 2-9 章习题的中文 VASP 解答，附输入文件（不含 POTCAR），可作练习参考。 |
| [gabkrenzer/PhononFlow](https://github.com/gabkrenzer/PhononFlow) | 35 | 未声明 | 2022-08-19 | VASP+phonopy 声子计算流程的技巧与脚本：声子谱、热力学、QHA、Grüneisen 参数、ModeMap 势能面。 |
| [gabkrenzer/vaspMD](https://github.com/gabkrenzer/vaspMD) | 29 | 未声明 | 2023-05-25 | VASP AIMD 与 MLFF 分子动力学流程的技巧与脚本，含 MLFF 基准测试和热容计算。 |
| [drinwater/Nudged-Elastic-Band-Tutorial](https://github.com/drinwater/Nudged-Elastic-Band-Tutorial) | 19 | 未声明 | 2022-01-28 | VASP+VTST CI-NEB 入门演练与模板，并演示用频率计算确认过渡态（唯一虚频）。 |
| [ishikawa-group/vasp-tutorial](https://github.com/ishikawa-group/vasp-tutorial) | 10 | 未声明 | 2025-09-15 | 日文 VASP 循序渐进教程：分子能量→体相/表面→优化→吸附→MD→空位→DOS，附编译笔记。 |
| [cpfpengfei/VASP-Examples](https://github.com/cpfpengfei/VASP-Examples) | 8 | 未声明 | 2020-06-18 | VASP 初学者笔记与示例：Si fcc、DOS/能带、INCAR/KPOINTS/POSCAR 要点与作业提交。 |
| [chouyoudou/mBJ-VASP-tutorial](https://github.com/chouyoudou/mBJ-VASP-tutorial) | 6 | 未声明 | 2022-04-08 | VASP 中 mBJ（METAGGA=MBJ）带隙计算的简明教程。 |

### 课程 / 工作坊

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [FedeDat/VASP-tutorial](https://github.com/FedeDat/VASP-tutorial) | 34 | 未声明 | 2022-10-30 | 面向电化学体系的 5 讲 VASP 入门课，配 YouTube 视频、上机练习和 VASP 官方讲义。 |
| [haiw201/dft-tutorial](https://github.com/haiw201/dft-tutorial) | 20 | 未声明 | 2026-05-07 | 中文 20 课 DFT+VASP 课程，从理论到能带/DOS/电荷/光学/磁性、缺陷/AIMD/声子/热导，附 INCAR 参数速查表和报错速查表。 |
| [skelton-group/VASP-Workshop-Chester-2019](https://github.com/skelton-group/VASP-Workshop-Chester-2019) | 8 | cc-by-4.0 | 2019-10-29 | 2019 年 Chester VASP 工作坊讲义与 6 个案例：Cu2ZnSnS4、琥珀酸分子晶体（vdW）、ZnO、NaCl(100) 表面吸附、PbS、BaTiO3。 |

### 文档 / 笔记

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [tamaswells/VASPKIT_manual](https://github.com/tamaswells/VASPKIT_manual) | 86 | 未声明 | 2019-05-29 | 旧版 VASPKIT 中文手册（v0.73）的 Markdown 源码，功能编号已有变动，请以 vaspkit.com 的 1.5.x 文档为准。 |
| [bzkarimi/VASP](https://github.com/bzkarimi/VASP) | 28 | 未声明 | 2020-07-14 | VASP 实用指南，解释 ISTART/ICHARG/ISMEAR/IBRION/ISIF 等常用 INCAR 参数，可用于核对参数含义。 |
| [northword/wiki](https://github.com/northword/wiki) | 25 | other | 2026-08-12 | 中文学习笔记：Linux/PBS、VASP 输入输出与结构优化、vaspkit 与 VTST 安装（含 VTST 版 VASP 编译）、Materials Studio，适合新手入门。 |

### 输入文件合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [JMSkelton/VASP-Examples](https://github.com/JMSkelton/VASP-Examples) | 36 | gpl-3.0 | 2016-10-30 | VASP 计算示例输入：Wannier90 插值的 HSE06 能带（Si），以及 SnS2 的 PBE/HSE06/G0W0/GW 对比。 |
| [vaspkit/VASPKIT_NatureProtocols](https://github.com/vaspkit/VASPKIT_NatureProtocols) | 4 | 未声明 | 2025-02-05 | VASPKIT 2025 年 Nature Protocols 论文配套示例文件，覆盖主要功能的输入与输出。 |

### 脚本合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [QijingZheng/VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding) | 295 | 未声明 | 2026-04-27 | VASP 输出进阶脚本集：解析 WAVECAR 得到实空间波函数、超胞能带反折叠、跃迁偶极矩、HSE k 点生成、BSE 激子分析等。 |
| [tamaswells/VASP_script](https://github.com/tamaswells/VASP_script) | 207 | 未声明 | 2021-12-28 | 中文社区常用 VASP 小脚本：换行修正、SIGMA 检查、能隙、ZPE/熵校正、固定原子、NEB 插点等，README 注明了已被 vaspkit 取代的功能。 |
| [raman-sc/VASP](https://github.com/raman-sc/VASP) | 191 | mit | 2026-01-07 | 以 VASP 为后端，用有限位移极化率法计算非共振 Raman 活性。 |
| [Chengcheng-Xiao/VASP_OPT_AXIS](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS) | 148 | 未声明 | 2025-11-27 | 给旧版 VASP 打补丁以在弛豫时固定部分晶格分量（如二维材料只优化面内）；VASP 6.4 起可直接用 LATTICE_CONSTRAINTS 标签代替。 |
| [Chengcheng-Xiao/VASP2WAN90_v2_fix](https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix) | 111 | 未声明 | 2023-09-19 | 更新版 VASP2WANNIER90v2 接口，修复 VASP 5 与 Wannier90 2.x 的接口问题（如 SOC 自旋轴），适用于仍在用 VASP 5.4.4 的情况。 |
| [ftherrien/VaspGibbs](https://github.com/ftherrien/VaspGibbs) | 64 | 未声明 | 2025-03-06 | 由 VASP 频率计算得到 ZPE、熵和吉布斯自由能，用于吸附/催化体系的自由能校正。 |
| [kavanase/vaspup2.0](https://github.com/kavanase/vaspup2.0) | 60 | mit | 2026-05-23 | bash 实现的 VASP 收敛测试流程：ENCUT/k 点的能量与力收敛，以及介电常数收敛（基于 input/CONFIG 模板）。 |
| [QijingZheng/VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF) | 47 | 未声明 | 2026-09-29 | 把 OUTCAR 中的振动模式转成 XSF 或动画，便于查看过渡态虚频的振动方向。 |
| [sylearn/Ueff_VASP](https://github.com/sylearn/Ueff_VASP) | 35 | 未声明 | 2025-07-24 | 用线性响应法（LDAUTYPE=3）在 VASP 中自洽求 DFT+U 的 U 值的流程。 |
| [Trebonius91/utils4VASP](https://github.com/Trebonius91/utils4VASP) | 33 | 未声明 | 2026-11-06 | 体相/界面 VASP 计算的建模与分析脚本及程序（涵盖 ML-FF），附 Wiki 使用说明。 |
| [keeeto/VASP-Elastic](https://github.com/keeeto/VASP-Elastic) | 27 | 未声明 | 2015-10-23 | 从 OUTCAR 提取完整弹性张量（IBRION=6, ISIF=3），计算体模量和剪切模量，并演示分量重排。 |
| [Stanford-MCTG/VASP-plot-modes](https://github.com/Stanford-MCTG/VASP-plot-modes) | 27 | 未声明 | 2022-08-14 | 把 DFPT（IBRION=7/8）计算的 OUTCAR 振动模式导出为 VESTA 图像。 |
| [jeffwdoak/vasp_scripts](https://github.com/jeffwdoak/vasp_scripts) | 25 | mit | 2017-12-10 | 对 VASP 输出做简单后处理与分析的 bash/Python 脚本集。 |
| [Bit-Part-Young/scripts](https://github.com/Bit-Part-Young/scripts) | 24 | 未声明 | 2026-11-03 | 上海交大作者维护的近 200 个中文注释脚本，涵盖 ASE、pymatgen、phonopy、vaspkit、弹性、MLIP 与 HPC。 |
| [henkelmangroup/vtstscripts](https://github.com/henkelmangroup/vtstscripts) | 4 | 未声明 | 2026-06-22 | VTST 官方脚本（nebmake.pl、nebresults.pl、nebbarrier.pl、dimer 系列、chgsum.pl 等），用于 NEB/Dimer 过渡态和 Bader 电荷分析。 |

### 工具 / 软件包

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen) | 2010 | other | 2026-11-06 | 材料计算 Python 基础库：VASP 输入集与输出解析、slab/k 路径、能带/DOS 绘图、相图、MP API，含 Gaussian/CP2K I/O；2026 年核心拆为 pymatgen-core。 |
| [phonopy/phonopy](https://github.com/phonopy/phonopy) | 506 | bsd-3-clause | 2026-11-10 | 声子计算标准工具（最新 4.10.0）：有限位移/DFPT 声子谱、DOS、热力学性质与 QHA，原生支持 VASP 和 CP2K；4.x 起要求 Python≥3.11。 |
| [materialsproject/atomate2](https://github.com/materialsproject/atomate2) | 360 | other | 2026-11-09 | 基于 jobflow 的模块化高通量工作流（弛豫、能带、弹性、声子、LOBSTER、缺陷、MLFF），以 VASP 为主，也支持 CP2K 等。 |
| [wannier-developers/wannier90](https://github.com/wannier-developers/wannier90) | 349 | lgpl-2.1 | 2026-11-08 | 最大局域化 Wannier 函数标准程序（2026 年发布 v4.0.x），VASP、CP2K 均有接口；VASP 6.2 起以 -DVASP2WANNIER90 链接 2.x/3.x 的 libwannier.a。 |
| [PytLab/VASPy](https://github.com/PytLab/VASPy) | 305 | mit | 2026-05-15 | 国内开发者编写的 VASP 文件读写库（POSCAR/OUTCAR/CHGCAR/XDATCAR、差分电荷），比 pymatgen 更轻量。 |
| [SMTG-Bham/doped](https://github.com/SMTG-Bham/doped) | 285 | mit | 2026-11-05 | 点缺陷一体化工具：超胞生成、有限尺寸电荷修正、形成能与转变能级分析，可与 ShakeNBreak 配合使用。 |
| [SMTG-Bham/sumo](https://github.com/SMTG-Bham/sumo) | 262 | mit | 2026-09-03 | 命令行出版级绘图工具（3.0.0）：sumo-kgen 生成 k 路径（含 HSE 零权重点），可画能带、DOS、光学谱和声子谱，bandstats 算有效质量。 |
| [dftd4/dftd4](https://github.com/dftd4/dftd4) | 232 | lgpl-3.0 | 2026-10-28 | DFT-D4 色散校正库，供 VASP IVDW=13 调用；VASP 6.6.0 已适配 4.x 新 API，3.7 以下旧版需加 -DDFTD4_API_V3。 |
| [romerogroup/pyprocar](https://github.com/romerogroup/pyprocar) | 225 | gpl-3.0 | 2026-10-10 | PROCAR 投影能带（fat band）、自旋纹理、费米面、能带反折叠与能带+DOS 组合图；v6 起 API 改为 bandsplot(code='vasp', dirname=...)。 |
| [materialsproject/custodian](https://github.com/materialsproject/custodian) | 190 | mit | 2026-11-02 | VASP/QChem/CP2K 作业自动纠错重启框架，VaspErrorHandler 收录大量真实报错（ZBRENT、BRMIX、EDDDAV、Sub-Space-Matrix 等）及修复策略。 |
| [henniggroup/VASPsol](https://github.com/henniggroup/VASPsol) | 183 | other | 2025-10-01 | VASP 隐式溶剂模型补丁（LSOL、EB_K 等标签），用于电催化/电化学体系；需要重新编译 VASP。 |
| [phonopy/phono3py](https://github.com/phonopy/phono3py) | 171 | bsd-3-clause | 2026-11-10 | 计算三阶力常数与晶格热导率（RTA/LBTE），用于声子-声子相互作用等进阶声子计算。 |
| [orest-d/p4vasp](https://github.com/orest-d/p4vasp) | 159 | gpl-2.0 | 2022-05-21 | 老牌 vasprun.xml 可视化 GUI（Python 2 时代），已基本停更。 |
| [zjwang11/IRVSP](https://github.com/zjwang11/IRVSP) | 138 | gpl-3.0 | 2026-05-26 | 计算 VASP 能带电子态的不可约表示，用于拓扑分析和能带对称性标注。 |
| [DigiMatChem/LobsterPy](https://github.com/DigiMatChem/LobsterPy) | 129 | bsd-3-clause | 2026-10-28 | LOBSTER 自动成键分析：生成 INCAR 与 lobsterin，自动绘制 COHP/COBI/COOP 并解析 ICOHP，支持 LOBSTER 5.x。 |
| [afonari/emc](https://github.com/afonari/emc) | 111 | other | 2020-10-12 | 有限差分法计算有效质量张量的经典工具：先生成 KPOINTS，再读取 EIGENVAL 求解。 |
| [fermisurfaces/IFermi](https://github.com/fermisurfaces/IFermi) | 104 | mit | 2026-11-01 | 由 VASP 结果生成、分析和可视化费米面，可叠加费米速度与自旋纹理，适合金属体系后处理。 |
| [vasp-dev/py4vasp](https://github.com/vasp-dev/py4vasp) | 101 | apache-2.0 | 2026-11-03 | VASP 官方 Python 后处理库（0.11.3），读取 vaspout.h5（需编译时开启 HDF5）画能带/DOS/声子/STM 等，VASP 6.x 后处理首选。 |
| [Infant83/VASPBERRY](https://github.com/Infant83/VASPBERRY) | 99 | 未声明 | 2026-10-09 | 从 VASP 的 WAVECAR 计算 Berry 曲率和陈数，用于拓扑性质分析。 |
| [WMD-group/MacroDensity](https://github.com/WMD-group/MacroDensity) | 92 | mit | 2026-01-10 | LOCPOT/CHGCAR 平面平均与宏观平均、真空能级提取，用于功函数计算和表面/界面能带对齐。 |
| [dftd3/simple-dftd3](https://github.com/dftd3/simple-dftd3) | 90 | lgpl-3.0 | 2026-10-28 | D3 色散校正库（D3(BJ)/D3(op) 等），VASP 6.6.0 新增的 IVDW=15 即通过该外部库调用。 |
| [lucydot/effmass](https://github.com/lucydot/effmass) | 86 | mit | 2026-05-11 | 从半导体能带计算多种定义的有效质量（抛物线拟合、非抛物线、光学有效质量）。 |
| [SMTG-Bham/easyunfold](https://github.com/SMTG-Bham/easyunfold) | 71 | mit | 2026-11-02 | 超胞（缺陷、掺杂、合金）能带反折叠工具，支持 VASP 和 CASTEP，命令行使用简单。 |
| [SMTG-Bham/surfaxe](https://github.com/SMTG-Bham/surfaxe) | 71 | mit | 2023-09-17 | 表面计算工具：自动生成 slab 并做层数/真空收敛测试，计算表面能与静电势。 |
| [shahramyalameha/ElATools](https://github.com/shahramyalameha/ElATools) | 58 | gpl-3.0 | 2026-07-21 | 由弹性常数矩阵分析二维/三维材料的各向异性杨氏模量和泊松比，并做 3D 可视化。 |
| [utf/pymlff](https://github.com/utf/pymlff) | 47 | mit | 2025-02-24 | 轻量 Python 包，用于读写与合并 VASP MLFF 训练集文件 ML_AB。 |
| [Ionizing/rsgrad](https://github.com/Ionizing/rsgrad) | 32 | mit | 2026-03-30 | Rust 编写的 VASP 作业监控工具，可快速查看 OUTCAR 收敛过程、各步受力和磁矩，也能画能带/DOS。 |
| [SWeav02/baderkit](https://github.com/SWeav02/baderkit) | 13 | bsd-3-clause | 2026-09-01 | Bader 电荷分析的 Python 实现，可替代 Henkelman 的 Fortran 版 bader，便于在脚本中调用。 |
| [henkelmangroup/vtstcode](https://github.com/henkelmangroup/vtstcode) | 5 | apache-2.0 | 2026-06-22 | VTST 补丁源码（CI-NEB、Dimer、Lanczos、IOPT/ML 优化器），按 VASP 版本分目录（6.4.3/6.5.1/6.6.0），需按版本打补丁并重新编译 VASP。 |

<a id="gaussian"></a>

## Gaussian

### 教程 / 书籍

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [ishikawa-group/gaussian-tutorial](https://github.com/ishikawa-group/gaussian-tutorial) | 4 | 未声明 | 2024-11-18 | 结构清晰的 Gaussian 讲义：能量、优化与频率（maxcycle/readFC/calcFC 等）、过渡态（opt=TS/modredundant）、TD 激发态、MD 与超算使用。 |
| [YuOZW/gaussian-tips](https://github.com/YuOZW/gaussian-tips) | 0 | 未声明 | 2023-09-02 | 日文 Gaussian 技巧合集（输入、输出、IRC、报错，如 TS 负本征值个数报错时加 NoEigenTest），另含 ORCA/xTB/Psi4/GRRM 输入。 |

### 课程 / 工作坊

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [thisisntnathan/dftCourse](https://github.com/thisisntnathan/dftCourse) | 8 | CC-BY-4.0 (stated in README; the API reports NOASSERTION) | 2022-01-26 | Cornell 的 Gaussian 入门课程网站（Short/Long Course）：Linux 与 SLURM、输入输出解读、几何优化、计算效率、Gen 基组、过渡态搜索、溶剂模型。 |

### 文档 / 笔记

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [lmmentel/awesome-python-chemistry](https://github.com/lmmentel/awesome-python-chemistry) | 1453 | other (NOASSERTION) | 2025-09-21 | 化学相关 Python 包的 awesome 列表（模拟、可视化、机器学习、学习资源），收录 GoodVibes、pyiron 等。 |
| [liyuanhe211/Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message) | 98 | CC-BY-NC-SA-4.0 (from the LICENSE file; the API reports NOASSERTION). Non-commercial: link to it, do not copy. | 2023-06-01 | 目前最全的 Gaussian 报错解决大全，按 Link 编号（L1/L101/L502/L9999 等）整理，并覆盖内存、磁盘、CPU 类错误，附中文版。 |

### 输入文件合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [thisisntnathan/dftCourseCodeRepo](https://github.com/thisisntnathan/dftCourseCodeRepo) | 4 | CC0-1.0 | 2022-01-11 | dftCourse 课程配套的 Gaussian 输入/输出文件与脚本示例。 |

### 脚本合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [polyluxus/tools-for-g16.bash](https://github.com/polyluxus/tools-for-g16.bash) | 35 | GPL-3.0 | 2020-01-20 | 成熟的 Gaussian 16 bash 工具集：chk2xyz、getenergy、getfreq（热化学汇总）、排队提交、testroute 语法检查，及派生 IRC/freq 等输入。 |
| [chauncey-garrett/gaussian-tools](https://github.com/chauncey-garrett/gaussian-tools) | 33 | 未声明 | 2017-04-20 | 一组 zsh 小脚本，各从 Gaussian log 提取一项：最低虚频、虚频个数、自旋污染（湮灭前/后）、热校正能量、ONIOM 外推能、点群，并可生成 SI。 |
| [chunxiangzheng/gaussian_log_file_converter](https://github.com/chunxiangzheng/gaussian_log_file_converter) | 16 | Apache-2.0 | 2017-11-30 | 把 Gaussian 优化 log 中的结构转换为 xyz 或 gjf 文件。 |
| [igor-1982/xDH4Gau](https://github.com/igor-1982/xDH4Gau) | 15 | MIT | 2022-09-21 | 用标准 Gaussian 计算 XYG3/XYGJ-OS 等 xDH 型双杂化泛函（徐昕课题组方法）的 Python 脚本。 |
| [polyluxus/runMultiwfn.bash](https://github.com/polyluxus/runMultiwfn.bash) | 11 | GPL-3.0 | 2019-09-10 | Linux 下非交互运行 Multiwfn 的 bash 封装（设置环境、喂入命令流），演示了 Multiwfn 批处理的思路。 |
| [alexamist/make_resp](https://github.com/alexamist/make_resp) | 10 | MIT | 2019-03-26 | 仅从 .xyz 出发，用 Gaussian（HF/6-31G* MK 静电势）加 AmberTools antechamber 拟合 RESP 电荷并生成 .mol2/.frcmod 的完整流程。 |
| [kimariyb/easy-Shermo](https://github.com/kimariyb/easy-Shermo) | 8 | MIT | 2026-06-12 | 批量调用 Shermo（Sobereva 开发的热化学程序），自动配对 opt/freq 与高精度单点文件并输出自由能，有中文 README。 |
| [MauricioCafiero/Gaussian16-Scripts-and-Functions](https://github.com/MauricioCafiero/Gaussian16-Scripts-and-Functions) | 4 | CC0-1.0 | 2025-01-03 | Gaussian 16 的 Slurm 模板与 submit_g16.sh，自动使输入中的 %mem/%nprocshared 与作业申请一致，并提供批量提交 runall。 |
| [Chempower-ZZH/auto_gjf](https://github.com/Chempower-ZZH/auto_gjf) | 3 | MIT | 2026-11-04 | gjf2gjf/xyz2gjf/log2gjf 三个脚本，用任务模板与泛函简写（如 B3BJ）批量生成/改写 Gaussian 输入，面向 WSL 用户。 |
| [jvalegre/IRC_to_com](https://github.com/jvalegre/IRC_to_com) | 3 | MIT | 2020-01-19 | 把 IRC 路径上的点（含两端点）按指定关键词转成新的 Gaussian 输入，便于端点优化或沿路径做高精度单点。 |
| [TonksGroupUMN/Tonks-Gaussian-Scripts](https://github.com/TonksGroupUMN/Tonks-Gaussian-Scripts) | 2 | MIT | 2022-11-07 | Gaussian 16 脚本集：Slurm/PBS 提交脚本、xyz2com、manual-IRC.py、HOMO-LUMO gap 提取、shake 去虚频等。 |
| [jcerezochem/traj2oniom](https://github.com/jcerezochem/traj2oniom) | 1 | MIT | 2024-04-30 | 从 MD 轨迹和拓扑生成带原子类型、层标记与冻结标记的 Gaussian ONIOM（QM/MM）输入。 |
| [smleechem/gaussjob_slurm](https://github.com/smleechem/gaussjob_slurm) | 0 | MIT | 2025-10-06 | Slurm 提交 Gaussian 16 的脚本：每个作业使用独立的 GAUSS_SCRDIR=$SCR_BASE/$SLURM_JOB_ID，捕获信号后自动清理 scratch。 |

### 工具 / 软件包

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [aligfellow/xyzrender](https://github.com/aligfellow/xyzrender) | 624 | MIT | 2026-11-09 | 从 xyz 或输出结构生成出版级分子图，适合论文配图和结果展示。 |
| [cclib/cclib](https://github.com/cclib/cclib) | 434 | BSD-3-Clause | 2026-11-09 | 标准的计算化学日志解析库，统一解析 Gaussian/ORCA 等输出（log/fchk）中的能量、几何、频率、激发态、轨道，可用 ccread/ccwrite 转换格式；建议用稳定版 1.8.1。 |
| [crest-lab/crest](https://github.com/crest-lab/crest) | 365 | LGPL-3.0 | 2026-06-16 | 构象搜索标准工具：用 GFN2-xTB 自动搜索低能构象，再交给 Gaussian 等做 DFT 优化。 |
| [duartegroup/autodE](https://github.com/duartegroup/autodE) | 219 | MIT | 2026-11-02 | 从 SMILES 自动生成反应剖面（构象、TS 猜测与优化、IRC 式验证），后端支持 Gaussian/ORCA/xtb。 |
| [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes) | 181 | MIT | 2026-11-06 | 由 Gaussian 频率输出计算准谐自由能（Grimme/Truhlar），支持标准态校正、频率缩放、单点能组合（--spc）、Boltzmann 加权；v4.4.0 新增动力学/能量跨度模块。 |
| [theochem/iodata](https://github.com/theochem/iodata) | 173 | LGPL-3.0 | 2026-11-10 | 读写并转换 fchk、log、wfn/wfx、molden、cube 等计算化学文件格式，并可生成 Gaussian 输入。 |
| [jvalegre/aqme](https://github.com/jvalegre/aqme) | 135 | MIT | 2026-09-28 | 从 SMILES 出发的自动化流程（v2.0.0）：RDKit/CREST 构象搜索、生成 Gaussian/ORCA 输入、自动检查虚频与报错并重算（QCORR）、提取描述符。 |
| [1234zou/MOKIT](https://github.com/1234zou/MOKIT) | 81 | Apache-2.0 | 2026-11-01 | 在 Gaussian（fch）、PySCF、ORCA、OpenMolcas 等程序间传递轨道，常用于把 broken-symmetry/UNO 结果交给多参考计算；国内邹靖翔开发，有中文说明。 |
| [zorkzou/Molden2AIM](https://github.com/zorkzou/Molden2AIM) | 65 | MIT | 2025-07-17 | 把 molden 文件转为 AIM-WFN/WFX 和 NBO .47 文件，便于不同程序的波函数进入 Multiwfn/NBO 分析。 |
| [QChASM/AaronTools.py](https://github.com/QChASM/AaronTools.py) | 63 | GPL-3.0 | 2026-10-22 | 功能全面的 Python 工具包：生成 Gaussian/ORCA 输入、解析输出、热化学、沿振动模式位移、UV-Vis/ECD 绘图、cube 输出、取代基/配体替换、Sterimol/埋藏体积。 |
| [digital-chemistry-laboratory/multiwfn-mac-build](https://github.com/digital-chemistry-laboratory/multiwfn-mac-build) | 53 | MIT | 2026-07-26 | 在 macOS 上编译 Multiwfn 的方案（另有 homebrew-multiwfn 仓库），供 Mac 用户安装使用。 |
| [gaox-qd/pysoc](https://github.com/gaox-qd/pysoc) | 50 | Apache-2.0 | 2024-05-31 | 读取 Gaussian TD-DFT 结果计算单重态-三重态旋轨耦合矩阵元（Gaussian 本身不算），用于 TADF/ISC 研究。 |
| [funkymunkycool/Cube-Toolz](https://github.com/funkymunkycool/Cube-Toolz) | 49 | 未声明 | 2023-07-06 | cube 文件差分、积分、平面平均、平移等操作，适用于 Gaussian cubegen 与 CP2K 密度/V_HARTREE cube，可算差分电荷密度和功函数。 |
| [patonlab/pyQRC](https://github.com/patonlab/pyQRC) | 45 | MIT | 2026-09-27 | 读取 Gaussian 频率结果并沿虚频方向位移生成新输入，可做“廉价 IRC”，也是消除优化后残留小虚频的常用办法。 |
| [eljost/thermoanalysis](https://github.com/eljost/thermoanalysis) | 35 | GPL-3.0 | 2024-10-21 | 独立的 Python 热化学实现（支持 ORCA/Gaussian），代码短小，便于对照 Gaussian 热化学白皮书理解 ZPE/H/S/G 各项的计算。 |
| [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian) | 30 | MIT | 2020-08-18 | 借助 External 让 Gaussian 调用 xtb，用 Gaussian 优化器在 GFN2-xTB 水平预搜索过渡态（opt=(calcfc,ts,noeigentest)、QST2）。 |
| [lenhanpham/OpenThermo](https://github.com/lenhanpham/OpenThermo) | 30 | MIT | 2026-09-20 | C++ 编写的热化学计算程序，读取 Gaussian/ORCA 等输出，支持 qRRHO、低频处理和温度/压力扫描，是 GoodVibes/Shermo 之外的另一选择。 |
| [tommason14/autochem](https://github.com/tommason14/autochem) | 27 | MIT | 2025-04-27 | GAMESS、ORCA、PSI4、Gaussian 多程序自动化脚本集（输入生成、结果提取、作业管理）。 |
| [maserasgroup-repo/pyssian](https://github.com/maserasgroup-repo/pyssian) | 20 | MIT | 2026-02-10 | 面向对象的 Gaussian 输入/日志解析库，按 Link 分块解析，附命令行工具，便于编写批处理脚本。 |
| [patonlab/pyDFTD3](https://github.com/patonlab/pyDFTD3) | 19 | MIT | 2026-02-17 | Python 版 Grimme D3 色散校正，可对已有结构单独计算 D3(0)/D3(BJ) 色散能，理解 EmpiricalDispersion 的贡献，或为不支持 GD3BJ 的泛函补算。 |
| [bane-dysta/banewfn](https://github.com/bane-dysta/banewfn) | 13 | Apache-2.0 | 2026-10-18 | 中文文档的 Multiwfn 工作流工具，把重复的菜单操作序列模板化为脚本，适合批量做 fchk 波函数分析。 |
| [TBruhn74/SpecDis](https://github.com/TBruhn74/SpecDis) | 12 | other (NOASSERTION) | 2017-10-13 | ECD/UV、ORD、VCD/IR 计算谱处理的经典工具：Boltzmann 加权、UV shift、对映体翻转、与实验谱算相似度因子，常用于 ECD 定绝对构型。 |
| [nanotech-empa/aiida-gaussian](https://github.com/nanotech-empa/aiida-gaussian) | 12 | MIT | 2026-10-15 | AiiDA 的 Gaussian 插件，在工作流中运行 Gaussian 并记录数据溯源，适合高通量与可复现计算。 |
| [njzjz/gaussianrunner](https://github.com/njzjz/gaussianrunner) | 8 | LGPL-3.0 | 2026-11-09 | 国内开发者编写，从 SMILES 批量生成并并行运行 Gaussian、汇总自由能，支持 MPI 跨节点。 |
| [mizu-bai/VSCode-Gaussian-Input-File-Support](https://github.com/mizu-bai/VSCode-Gaussian-Input-File-Support) | 8 | BSD-2-Clause | 2023-12-01 | VSCode 的 Gaussian 输入文件（.gjf）语法高亮插件，可减少关键词拼写错误（QPErr）。 |
| [shin1koda/dmf-g16](https://github.com/shin1koda/dmf-g16) | 8 | LGPL-2.1 | 2026-08-16 | PyDMF 的 Gaussian 封装，用直接 MaxFlux 方法做双端过渡态搜索，是 QST2/QST3 之外的现代替代方案。 |
| [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot) | 8 | MIT | 2026-07-20 | 中文 Codex 工作流插件：生成 Gaussian 输入与 Slurm 脚本（优化、频率、TS/IRC、TD-DFT、NMR 等），提供报错诊断和结果自查清单（虚频、自旋污染等）。 |
| [LeonardoESousa/LeoX](https://github.com/LeonardoESousa/LeoX) | 6 | MIT | 2026-10-08 | 基于 Gaussian TD(A)-DFT 用核系综法模拟吸收/发射谱、估算激子性质，并做长程校正参数 ω 调优（optimal tuning）。 |
| [joaquinbarroso/ELOTE](https://github.com/joaquinbarroso/ELOTE) | 6 | MIT | 2026-08-05 | 分析 Gaussian TD-DFT 激发的电子定域与跃迁特征，辅助判断电荷转移（CT）激发。 |
| [lenhanpham/Gaussian-Extractor](https://github.com/lenhanpham/Gaussian-Extractor) | 5 | MIT | 2025-11-14 | C++ 编写的高性能并行 log 数据提取工具，可检查作业状态、汇总能量、生成新输入，适合大批量结果整理。 |

<a id="cp2k"></a>

## CP2K

### 教程 / 书籍

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [compchem-cybertraining/Tutorials_CP2K](https://github.com/compchem-cybertraining/Tutorials_CP2K) | 68 | CC0-1.0 | 2021-12-16 | 以二维钙钛矿 (BA)2PbI4 为例，覆盖收敛测试、几何/晶胞优化、TDDFT、杂化泛函、大体系、MD 与非绝热耦合（CP2K/Libra）；基于 v6.1/7.1，部分关键字已过时。 |
| [bioexcel/gromacs-2022-cp2k-tutorial](https://github.com/bioexcel/gromacs-2022-cp2k-tutorial) | 45 | 未声明 | 2026-01-29 | BioExcel 的 GROMACS–CP2K QM/MM 接口教程，附 Handbook PDF 及 eGFP、stilbene、NMA 示例。 |
| [bioexcel/cp2K_qmmm_tutorials_for_biological_simulations](https://github.com/bioexcel/cp2K_qmmm_tutorials_for_biological_simulations) | 26 | 未声明 | 2019-10-18 | CP2K 原生 &QMMM 生物体系教程（改编自 AMBER/GMX/NAMD/CPMD），讲解 AMBER prmtop 拓扑、link atom、电荷修正等输入搭建，版本较老。 |
| [arvpinto/enzyme_neb_cp2k](https://github.com/arvpinto/enzyme_neb_cp2k) | 10 | CC0-1.0 | 2026-10-25 | 酶反应 QM/MM 过渡态教程（2026 年仍在更新）：CI-NEB、Dimer、频率验证，附 cp2k_energy_analysis.sh 和 cp2k_freq_analysis.sh。 |
| [ishikawa-group/cp2k-tutorial](https://github.com/ishikawa-group/cp2k-tutorial) | 4 | 未声明 | 2024-09-18 | 日文 CP2K 入门教程：安装、能量、优化、MD、元动力学、泛函，含苯分子与 Pt 表面示例。 |
| [hyfaust/cp2k-learn](https://github.com/hyfaust/cp2k-learn) | 3 | GPL-3.0 | 2026-06-09 | 16 章分级 CP2K 实操教程（中文 README，30+ 输入经 v2025.2 验证），从单点、基组、截断到 MD、元动力学、xTB、QM/MM、TDDFT、能带，附交互网页。 |
| [svarner9/CP2K-Workshop](https://github.com/svarner9/CP2K-Workshop) | 2 | MIT | 2025-06-13 | 用 ASE 和 PySAGES 驱动 CP2K 做增强采样（偏置）AIMD 的 workshop 示例。 |
| [ChiahsinChu/cp2k-user-tutorial](https://github.com/ChiahsinChu/cp2k-user-tutorial) | 0 | GPL-3.0 | 2026-03-10 | 面向零基础的“CP2K 作为第一个 DFT 程序”教程：能量计算、密度混合、slurm 模板与基组文件，并汇总中文学习资源。 |
| [wkpahuwkp98/cp2k-nac-fssh](https://github.com/wkpahuwkp98/cp2k-nac-fssh) | 0 | CC-BY-4.0 | 2026-10-22 | 中文非绝热动力学进阶教程，从 Tully FSSH 推导到 PYXAID/CP2K 实操，附 CP2K 高斯基 NAC 实现。 |

### 课程 / 工作坊

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [bioexcel/2020-10-cp2k_qmmm_intro](https://github.com/bioexcel/2020-10-cp2k_qmmm_intro) | 5 | 未声明 | 2022-06-16 | ARCHER2 的 CP2K 生物分子 QM/MM 实践课材料，可与 BioExcel 其他 CP2K QM/MM 教程互补。 |
| [csc-training/sscc-cp2k](https://github.com/csc-training/sscc-cp2k) | 4 | CC-BY-4.0 | 2026-03-01 | 完整的反应研究案例（Jupyter）：PBE 下做 NEB，对过渡态做振动分析，用 PBE0 重新收敛 MEP，最后做溶液中 AIMD。 |
| [EPCCed/intro_to_cp2k](https://github.com/EPCCed/intro_to_cp2k) | 3 | 未声明 | 2024-04-08 | EPCC/ARCHER2 官方 CP2K 入门课，5 个练习：能量、截断收敛、几何优化、MD、并行性能。 |

### 文档 / 笔记

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [chenggroup/chenggroup.github.io](https://github.com/chenggroup/chenggroup.github.io) | 111 | MIT | 2026-11-10 | 厦大程俊组中文 Wiki，CP2K 专题含入门、能量/力、几何优化、DFT+U、NEB、slab、constrained MD、REFTRAJ、ZPE、杂化泛函、DeePMD 联用及工具介绍。 |
| [hfp/xconfigure](https://github.com/hfp/xconfigure) | 45 | BSD-3-Clause | 2026-11-06 | Intel 工具链下编译 CP2K 等软件的高性能配置配方；新版 CP2K 已改用 CMake，仅作历史参考。 |
| [DCM-Uni-Paderborn/CP2G-Decalogo](https://github.com/DCM-Uni-Paderborn/CP2G-Decalogo) | 3 | 未声明 | 2026-05-18 | Kühne 组 LiveCoMS 最佳实践论文源码，讲解 CP2K 第二代 Car-Parrinello AIMD 的步长、ASPC 外推阶数、Langevin 校正与能量守恒判据。 |
| [szt2006/cp2k-aimd-skill](https://github.com/szt2006/cp2k-aimd-skill) | 3 | 未声明 | 2026-09-27 | 中文 CP2K AIMD AI skill 知识库，references/playbook.md 提供“症状→处方”对照表（SCF 不收敛、几何振荡、报错原文），并摘录官方文档。 |
| [DCM-Uni-Paderborn/UZH-protocol](https://github.com/DCM-Uni-Paderborn/UZH-protocol) | 1 | 未声明 | 2026-10-03 | UZH 基组与赝势验证工作流的数据（对应 arXiv 2606.11064，与 BASIS_MOLOPT_UZH_2026.2 对齐），说明了新项目为何推荐 UZH 组合。 |

### 输入文件合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [cp2k/cp2k-examples](https://github.com/cp2k/cp2k-examples) | 22 | 未声明 | 2026-04-26 | CP2K 官方示例库：converging_grid 截断收敛范例、cp2k-made-simple 论文配套输入，以及元动力学、QM/MM、RI-HFX、GW/BSE、CDFT 等专题。 |
| [juerghutter/BASIS](https://github.com/juerghutter/BASIS) | 12 | 未声明 | 2026-10-16 | CP2K 核心开发者 Jürg Hutter 维护的基组集合（2026 年仍在更新），有助于理解 MOLOPT/UZH 基组的来源。 |
| [bwang-ecnu/gromacs_cp2k_tutorial](https://github.com/bwang-ecnu/gromacs_cp2k_tutorial) | 4 | 未声明 | 2021-02-01 | 华东师大王斌举组的 GROMACS–CP2K QM/MM 接口示例体系，方便中文读者对照学习。 |

### 脚本合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [LePingKYXK/Dipole-ACF](https://github.com/LePingKYXK/Dipole-ACF) | 21 | 未声明 | 2023-06-19 | 由 CP2K AIMD 总偶极自相关函数的 FFT 计算 IR 光谱的经典脚本（基于 Python 2.7）。 |
| [RobinLi-95/CP2K_tools](https://github.com/RobinLi-95/CP2K_tools) | 9 | 未声明 | 2021-09-28 | CP2K 脚本：EOS 拟合、能带作图、V_HARTREE cube 平面平均（vtotav）、功函数计算、ghost layer slab 与吸附物放置。 |
| [Doria58/DFT_Script](https://github.com/Doria58/DFT_Script) | 9 | 未声明 | 2023-11-01 | 中文 VASP/CP2K/QE 数据处理脚本：band.py 把 CP2K 能带转为 Origin 可用数据，ToPOSCAR 转为 VASP 结构，FindKPATHForCP2K 生成能带路径。 |
| [ma455173220/cp2k_scripts](https://github.com/ma455173220/cp2k_scripts) | 1 | 未声明 | 2026-10-28 | 中文 README 的 CP2K 脚本集：cutoff/rel_cutoff/U/k 点收敛、NEB 初始化与分析、restart 回写、任务类型一键切换、自动识别 PBS/Slurm 提交。 |
| [mmm-group/dos-broadening](https://github.com/mmm-group/dos-broadening) | 1 | MIT | 2019-04-05 | 对 CP2K 旧格式 .pdos 做展宽的小脚本；2026.2 起 CP2K 的 &DOS/&CURVE/&BROADEN 已可直接输出展宽曲线。 |

### 工具 / 软件包

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [cp2k/cp2k](https://github.com/cp2k/cp2k) | 1239 | GPL-2.0 | 2026-11-09 | CP2K 官方源码：docs/ 为官方手册源文件，tests/ 下数千个 regtest 输入可查关键字的真实写法，data/ 含 BASIS_MOLOPT_UZH 等基组与赝势。 |
| [molmod/psiflow](https://github.com/molmod/psiflow) | 144 | MIT | 2026-11-09 | Parsl 驱动的可扩展分子模拟工作流，以 CP2K 为参考 QM 引擎训练和验证机器学习势。 |
| [SINGROUP/pycp2k](https://github.com/SINGROUP/pycp2k) | 106 | LGPL-3.0 | 2022-06-07 | 以面向对象方式生成完整 CP2K 输入；自带类只对应老版本，新版本需用 cp2k --xml 导出 XML 重新生成，已基本停止维护。 |
| [chenggroup/ai2-kit](https://github.com/chenggroup/ai2-kit) | 100 | MIT | 2026-11-05 | 程俊组的 AI 与第一性原理工具链，以 CP2K 为标注引擎做 DP-GEN 式势函数训练、批量作业管理与质子转移分析等。 |
| [robinzyb/cp2kdata](https://github.com/robinzyb/cp2kdata) | 80 | LGPL-3.0 | 2026-02-21 | 国内作者维护的 CP2K 后处理库，读取输出、cube 和 pdos，带 dpdata 插件可把 AIMD 数据转为 DeePMD 训练集；输出格式变化时（如 2025.1）需先验证解析结果。 |
| [elcorto/pwtools](https://github.com/elcorto/pwtools) | 73 | BSD-3-Clause | 2026-09-12 | QE、CPMD、CP2K、LAMMPS 通用的前后处理库，提供输出解析器和数据结构。 |
| [cp2k/cp2k-input-tools](https://github.com/cp2k/cp2k-input-tools) | 63 | MIT | 2026-11-09 | CP2K 官方输入工具：cp2klint 语法检查、inp 与 JSON/YAML 互转、cp2kgen 批量生成扫描输入、cp2kget 读 restart，另有编辑器补全语言服务器。 |
| [avishart/CP2K_Editor](https://github.com/avishart/CP2K_Editor) | 60 | 未声明 | 2019-07-15 | 老牌 CP2K 输入文件 GUI 生成器，已不再维护。 |
| [cp2k/cp2k-output-tools](https://github.com/cp2k/cp2k-output-tools) | 55 | MIT | 2026-11-09 | CP2K 输出处理工具：cp2kparse 转 JSON/YAML，cp2k_bs2csv 转能带 CSV，cp2k_pdos 展宽 PDOS，xyz_restart_parser 去除续算重复帧。 |
| [jhaens/amaceing_toolkit](https://github.com/jhaens/amaceing_toolkit) | 27 | 未声明 | 2026-09-01 | 同时生成 CP2K 输入与 MACE、MatterSim、SevenNet、ORB、Grace 等机器学习势训练输入，并提供后处理。 |
| [aiidateam/aiida-cp2k](https://github.com/aiidateam/aiida-cp2k) | 26 | MIT | 2026-11-03 | AiiDA 的 CP2K 插件，以 Python 嵌套字典写出 CP2K section，适合高通量与可溯源计算；v2.2.1（2026-11）维护活跃。 |
| [cp2k/cp2k-containers](https://github.com/cp2k/cp2k-containers) | 17 | GPL-2.0 | 2026-08-20 | 用 Apptainer/Singularity 拉取运行 CP2K 镜像的说明（含 run_tests 自检、CUDA 镜像）；2026.2 起 Dockerfile 迁至主仓库 tools/docker。 |
| [JunboLu/CP2K_kit](https://github.com/JunboLu/CP2K_kit) | 15 | Apache-2.0 | 2022-11-13 | 国内作者开发的 CP2K 辅助工具：输入生成、轨迹分析、GTH 赝势与基组优化、热力学积分、DeePMD 联用。 |
| [nanotech-empa/cp2k-spm-tools](https://github.com/nanotech-empa/cp2k-spm-tools) | 15 | MIT | 2026-11-01 | 读取 CP2K .wfn 做 STM/STS 扫描探针模拟和轨道杂化分析，适合表面科学研究。 |
| [KMNitesh05/cp2k_2_deepmdkit](https://github.com/KMNitesh05/cp2k_2_deepmdkit) | 10 | MIT | 2024-08-19 | 把 CP2K AIMD 的 pos、frc、ener、cell 输出转为 DeepMD-kit 的 npy 格式，可与 cp2kdata/dpdata 对照使用。 |
| [ybx20021218/CP2Kview--Visualization-of-CP2K-Output-Files](https://github.com/ybx20021218/CP2Kview--Visualization-of-CP2K-Output-Files) | 9 | AGPL-3.0 | 2026-06-17 | 零依赖浏览器工具，可视化 CP2K 优化结构与轨迹、分析 OUT/LOG 收敛曲线、查看 INP/restart，适合初学者检查输出。 |
| [Shxnek/cp2k-windows-installer](https://github.com/Shxnek/cp2k-windows-installer) | 8 | MIT | 2026-04-24 | 中文 CP2K Windows 一键安装包，适合在笔记本上演示学习；正式计算仍建议用 Linux、WSL 或容器。 |
| [superstar54/xcp2k](https://github.com/superstar54/xcp2k) | 8 | 未声明 | 2021-11-06 | ASE 风格的 CP2K 接口，可直接访问全部输入段，可与 ASE 官方 CP2K 计算器对比选用。 |
| [pierre-24/cp2k-basis](https://github.com/pierre-24/cp2k-basis) | 6 | 未声明 | 2026-03-10 | CP2K 官网 cp2k-basis 工具背后的项目，解析基组与 GTH 赝势，可按元素查找匹配的基组+赝势组合。 |
| [conda-forge/cp2k-feedstock](https://github.com/conda-forge/cp2k-feedstock) | 4 | BSD-3-Clause | 2026-10-21 | CP2K 的 conda-forge 配方；2026.2 仅提供 linux-64 和 linux-aarch64，osx-64 停在 7.1.0，Mac 用户建议改用容器。 |
| [ljpyztd-cyber/CP2K_Forge](https://github.com/ljpyztd-cyber/CP2K_Forge) | 3 | 未声明 | 2026-11-05 | Windows 桌面 GUI，调用 Multiwfn 生成 CP2K 输入（目标 2026.1/2026.2），支持可视化、批量生成和 Slurm 脚本；使用时需引用 Multiwfn 的两篇论文。 |
| [CommonCE/cp2kkit](https://github.com/CommonCE/cp2kkit) | 2 | 未声明 | 2026-11-06 | 类似 VASPKIT 的 CP2K 前后处理工具箱：读取输出、cube 差分与平面平均、生成输入（目标版本 2026.2），目前仍为预发布版。 |
| [LouieSlocombe/ase-cp2k-noshell](https://github.com/LouieSlocombe/ase-cp2k-noshell) | 1 | MIT | 2025-01-27 | 不依赖 cp2k_shell 的 ASE CP2K 计算器，适用于部分发行版或 CMake 构建缺少 cp2k_shell 的情况。 |
| [Clementine-C/cp2kbuild](https://github.com/Clementine-C/cp2kbuild) | 0 | MIT | 2026-11-08 | 中文 CP2K 输入工具链（GUI+CLI）：由 xyz/pdb/POSCAR/cif 生成 inp，自动匹配基组赝势，生成截断收敛脚本，内置 2026.2 手册校验关键字。 |

<a id="openmx"></a>

## OpenMX

### 教程 / 书籍

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [pranabdas/openmx](https://github.com/pranabdas/openmx) | 9 | 未声明 | 2026-09-05 | 内容完整的英文 OpenMX 教程站：安装、能带/DOS、收敛测试、DFT+U、SOC、轨道投影、磁矩可视化，附 4.0+patch4.0.1 一键编译脚本与 Python notebook。 |
| [hikaruri/DFT_startup](https://github.com/hikaruri/DFT_startup) | 2 | 未声明 | 2021-06-19 | 日文研究室入门资料：macOS 上用 MacPorts 编译 OpenMX 3.9+ScaLAPACK、OpenMX 基础教程、gnuplot 画三维能带，以及 CrySPY 用法。 |
| [hirtatsu/Calculations](https://github.com/hirtatsu/Calculations) | 0 | 未声明 | 2026-08-31 | 大阪大学的日文实操笔记：OpenMX 输入文件逐段讲解、tmux/nohup 运行与日志检查、gnuplot 作图，另含 SQUID 超算编译与 PBS 脚本、WSL/oneAPI 环境配置。 |
| [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-) | 0 | 未声明 | 2026-08-06 | 印尼语 OpenMX 教程，按分支组织主题：晶格优化、NEB、电荷与静电势、能带/DOS、磁性与居里温度、SOC/MAE、费米面、Rashba、Berry 曲率、反常霍尔、边缘态、Z2。 |

### 文档 / 笔记

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [hikaruri/Z2-support-OpenMX](https://github.com/hikaruri/Z2-support-OpenMX) | 1 | GPL-3.0 | 2020-10-03 | OpenMX 自带 Z2FH（用 Berry 曲率/联络计算 Z2 不变量）的使用支持页。 |

### 脚本合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [TnakaTnakaTnaka/TEP_tools](https://github.com/TnakaTnakaTnaka/TEP_tools) | 9 | 未声明 | 2019-02-14 | Python 脚本，读取 OpenMX NEGF 透射谱计算 BZ 平均透射、电导、Seebeck 系数、功率因子与 ZT。 |
| [mfukudaQED/openmx_installer](https://github.com/mfukudaQED/openmx_installer) | 7 | 未声明 | 2025-05-13 | OpenMX 一键下载、打补丁、改 makefile 并编译的脚本（WSL+oneAPI、ISSP ohtaka），含 oneAPI 2025 后 mpiicx/mpiifx 必需的 -fcommon、-Wno-error 等参数。 |
| [tsaie79/openmx-wf](https://github.com/tsaie79/openmx-wf) | 1 | 未声明 | 2024-06-08 | 基于 atomate/FireWorks 的 OpenMX 工作流示例，可用于高通量计算。 |
| [knhr0727/pymx](https://github.com/knhr0727/pymx) | 0 | 未声明 | 2022-04-02 | OpenMX Python 后处理库，读取 scfout 构建 H(k)，示例含 Bi2Se3 能带、宇称本征值与 Wannier 电荷中心（WCC/Wilson loop）。 |
| [samlipton/pyfesm](https://github.com/samlipton/pyfesm) | 0 | 未声明 | 2026-04-06 | OpenMX 轻量单类 Python 接口：解析 SCF/体系信息、能带、Mulliken 电荷、NEGF 透射并统一单位。 |

### 工具 / 软件包

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [mzjb/DeepH-pack](https://github.com/mzjb/DeepH-pack) | 346 | LGPL-3.0 | 2024-10-07 | 用深度神经网络预测 DFT 哈密顿量（DeepH），以 OpenMX 3.9 生成训练数据，文档说明了 HS.fileout/overlap 的用途。 |
| [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata) | 255 | LGPL-3.0 | 2026-11-09 | 在 VASP/CP2K/Gaussian/ABACUS/OpenMX 等程序间互转轨迹、能量和力数据（可读 OpenMX .dat+.md），常用于准备 DeePMD 等机器学习势训练数据。 |
| [QuantumLab-ZY/HamGNN](https://github.com/QuantumLab-ZY/HamGNN) | 231 | GPL-3.0 | 2026-09-18 | 预测电子哈密顿量的 E(3) 等变图神经网络，支持 OpenMX 数据，附 openmx_postprocess/read_openmx 编译流程。 |
| [alamode-code/alamode](https://github.com/alamode-code/alamode) | 212 | MIT | 2026-10-30 | 晶格非谐性与热输运计算工具，接口支持 VASP、QE、OpenMX、xTAPP、LAMMPS，可弥补 OpenMX 本身没有 DFPT 的不足。 |
| [Tomoki-YAMASHITA/CrySPY](https://github.com/Tomoki-YAMASHITA/CrySPY) | 156 | MIT | 2026-11-09 | Python 晶体结构预测工具，可调用 OpenMX 做结构弛豫。 |
| [Xiaoxun-Gong/DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3) | 119 | MIT | 2026-01-27 | DeepH 的 E(3) 等变版本，使用 OpenMX 计算数据。 |
| [mailhexu/TB2J](https://github.com/mailhexu/TB2J) | 105 | BSD-2-Clause | 2026-11-08 | 计算磁相互作用参数的 Python 包（TB2J 主仓库），文档有 OpenMX 专节（tb2j.readthedocs.io/en/latest/src/openmx.html）。 |
| [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl) | 64 | MIT | 2023-11-19 | Julia 紧束缚包，可从 OpenMX（及 Wannier90、FHI-aims）哈密顿量构建模型，计算 Berry 曲率、光学与非线性响应。 |
| [wxia529/acoord](https://github.com/wxia529/acoord) | 23 | MIT | 2026-07-17 | 在 VS Code 内可视化与编辑原子结构的插件，读写 OpenMX .dat 并保留参数块。 |
| [issp-center-dev/abICS](https://github.com/issp-center-dev/abICS) | 21 | GPL-3.0 | 2026-11-10 | 无序体系构型采样工具，可用 OpenMX 作为求解器并生成训练数据。 |
| [pulkin/openmx-hks](https://github.com/pulkin/openmx-hks) | 17 | 未声明 | 2021-01-12 | C 工具，读取 OpenMX NEGF 产生的 .hks（Kohn-Sham 哈密顿量），导出哈密顿量/重叠矩阵（h5/json/MATLAB）和 xsf 结构，附能带与弹道输运示例。 |
| [hikaruri/OMXsigmaxy](https://github.com/hikaruri/OMXsigmaxy) | 12 | GPL-3.0 | 2025-09-23 | 基于扩展 FHS 方法，从 OpenMX .scfout（需 HS.fileout on）计算反常霍尔电导，用与 OpenMX 相同的 CC/LIB 编译。 |
| [mailhexu/TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX) | 11 | GPL-3.0 | 2026-04-07 | TB2J 的 OpenMX 接口，用磁力定理从 scfout 计算海森堡交换 J、DMI、各向异性交换与磁振子能带。 |
| [dc1394/openmx3.9.9_gpu](https://github.com/dc1394/openmx3.9.9_gpu) | 9 | GPL-3.0 | 2026-08-04 | OpenMX 3.9.9 的 GPU 加速移植版，适合仍在使用 3.9 的用户。 |
| [mzjb/overlap-only-OpenMX](https://github.com/mzjb/overlap-only-OpenMX) | 9 | GPL-3.0 | 2022-07-19 | 只计算重叠矩阵的 OpenMX 补丁，用于大体系机器学习哈密顿量推理时快速获得 overlap。 |
| [issp-center-dev/cif2x](https://github.com/issp-center-dev/cif2x) | 8 | GPL-3.0 | 2026-10-30 | ISSP 开发的工具，基于 pymatgen 由 CIF+模板批量生成 QE、VASP、OpenMX、AkaiKKR、RESPACK 输入，适合高通量计算。 |
| [by-student-2017/cif2cell-informal](https://github.com/by-student-2017/cif2cell-informal) | 7 | GPL-3.0 | 2021-10-01 | cif2cell 非官方分支，支持 -p openmx 输出与 LDA+U 选项，可指定 PAO 目录。 |
| [KAIST-ELST/DFTforge.jl](https://github.com/KAIST-ELST/DFTforge.jl) | 7 | 未声明 | 2025-02-16 | Julia 框架，读取 OpenMX scfout/Wannier90 结果计算交换耦合 J(q) 等，是 OpenMX Jx 方法的扩展实现。 |
| [Ncmexp2717/OpenMX-deb](https://github.com/Ncmexp2717/OpenMX-deb) | 6 | GPL-3.0 | 2021-11-07 | 非官方的早期 OpenMX 3.9 deb 包（DATA.PATH 默认为 /usr/share/openmx/DFT_DATA19）；4.0 已有官方 deb。 |
| [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu) | 4 | GPL-3.0 | 2026-11-10 | OpenMX 官方下载页链接的 4.0 GPU 版，CMake 可选 CUDA（sm_80/90/120）、HIP（gfx942 等）或 CPU 后端，并扩展了镧系元素 VPS/PAO。 |
| [Ncmexp2717/OMXTool](https://github.com/Ncmexp2717/OMXTool) | 4 | GPL-3.0 | 2023-11-12 | OpenMX 的 Tcl/Tk 图形工具：预览 .Band 能带、导出 gnuplot 图、转换 scfout 字节序。 |
| [mitsuaki1987/cif2input](https://github.com/mitsuaki1987/cif2input) | 3 | 未声明 | 2025-10-24 | CIF 转 QE 与 OpenMX 输入（含原胞与高对称路径），可与官方 cif2omx 对比使用。 |
| [KazumeNishidate/tube](https://github.com/KazumeNishidate/tube) | 3 | 未声明 | 2026-11-08 | 纳米管结构生成器，直接输出 OpenMX .dat 和 QE 输入，附能带计算示例。 |
| [pohao82/poscar2openmx](https://github.com/pohao82/poscar2openmx) | 2 | MIT | 2026-10-23 | POSCAR 转 OpenMX .dat：按 n_i·d_i≥35 Å 自动设 scf.Kgrid，可选 Quick/Standard/Precise 基组，支持 MAGMOM/非共线磁矩。 |
| [Ncmexp2717/HWC](https://github.com/Ncmexp2717/HWC) | 2 | GPL-3.0 | 2017-05-11 | OpenMX 附加程序：复制进 source 后 make HWC，读取 scfout 用平行输运计算混合 Wannier 中心与层极化。 |
| [cometscome/OpenMX_Docker](https://github.com/cometscome/OpenMX_Docker) | 1 | 未声明 | 2021-10-20 | 用 Intel oneAPI 安装 OpenMX 的最小 Dockerfile（运行需 --shm-size=2gb）；4.0 官方论坛另有 x86/ARM 双架构镜像。 |
| [pohao82/vasp-omx-moment-setter](https://github.com/pohao82/vasp-omx-moment-setter) | 1 | MIT | 2025-11-22 | Dash/Plotly 网页应用，交互设置并可视化 FM/AFM/非共线初始磁矩，写入 VASP/OpenMX 输入。 |
| [nnardoto/LowFoldIn](https://github.com/nnardoto/LowFoldIn) | 1 | GPL-3.0 | 2026-08-04 | 原生读取 OpenMX 3.9/4.0 .scfout（非自旋/共线/非共线），做能带计算与能带展开。 |
| [azadoks/aiida-openmx](https://github.com/azadoks/aiida-openmx) | 1 | MIT | 2021-07-29 | OpenMX 的 AiiDA 工作流插件（较旧），可用于自动化/高通量计算。 |
| [atsushi-m-ito/openmx-patch-oneapi](https://github.com/atsushi-m-ito/openmx-patch-oneapi) | 0 | GPL-3.0 | 2025-07-08 | Intel oneAPI 移除 icc 后用 icx/ifx 编译 OpenMX 的补丁，日文说明列出需添加 lapack_prototypes.h 的源文件。 |
| [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS) | 0 | MPL-2.0 | 2026-07-04 | 单 HTML 浏览器工具，交互显示 OpenMX 能带+DOS/PDOS 与展开能带，标记带隙/VBM/CBM，支持 PDOS 积分和出版级样式。 |
| [zeleznyj-htp/aiida-plugin-openmx](https://github.com/zeleznyj-htp/aiida-plugin-openmx) | 0 | 未声明 | 2026-10-07 | 较新的 OpenMX AiiDA 插件，含 JxCalculation 交换常数计算（需 HS.fileout on 与自旋极化）和带重启的 OpenMXWorkchain。 |
| [mfukudaQED/FLPQ](https://github.com/mfukudaQED/FLPQ) | 0 | GPL-3.0 | 2026-10-24 | 配合 OpenMX 4.0 新增的 RCP（区域化学势）分析使用的程序。 |

### 其他

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [cmsi/MateriAppsLive](https://github.com/cmsi/MateriAppsLive) | 57 | 未声明 | 2025-10-03 | 预装 MateriApps 软件（含 OpenMX 3.9.9、bandgnu13、DosMain）的 Debian Live 系统，适合不想编译的初学者体验。 |
| [klsk007/openmx-soc-cddf-fix](https://github.com/klsk007/openmx-soc-cddf-fix) | 1 | 未声明 | 2026-05-29 | OpenMX 4.0 非共线+SOC+CDDF 光学计算报 'Parameter 10 was incorrect on entry to ZHEEVX' 或段错误的问题记录与补丁。 |

<a id="general"></a>

## 通用工具与学习资源（General tools & learning resources）

### 教程 / 书籍

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [stanfordbshan/CompMatBook](https://github.com/stanfordbshan/CompMatBook) | 139 | GPL-3.0 | 2026-05-09 | 开源计算材料学教材，覆盖从理论到计算方法，适合作为理论延伸阅读。 |
| [computron/pymatgen_tutorials](https://github.com/computron/pymatgen_tutorials) | 83 | NOASSERTION (CC-BY-4.0) | 2026-09-29 | pymatgen 作者 Anubhav Jain 的五部分视频+notebook 教程：结构、VASP I/O、MP API、对称性变换。 |
| [aiidateam/aiida-tutorials](https://github.com/aiidateam/aiida-tutorials) | 25 | 未声明 | 2025-10-02 | AiiDA 官方教程站点源码，适合学习 AiiDA 工作流。 |

### 课程 / 工作坊

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [materialsproject/workshop](https://github.com/materialsproject/workshop) | 122 | BSD-3-Clause | 2023-03-01 | Materials Project 官方 workshop 课程材料，涵盖 pymatgen、MP API、atomate 与机器学习。 |

### 文档 / 笔记

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [hsiaoyi0504/awesome-cheminformatics](https://github.com/hsiaoyi0504/awesome-cheminformatics) | 896 | 未声明 | 2024-03-15 | 化学信息学 awesome 列表（RDKit、Open Babel 等），做分子建模时可查找相关工具。 |
| [JuDFTteam/best-of-atomistic-machine-learning](https://github.com/JuDFTteam/best-of-atomistic-machine-learning) | 737 | CC-BY-SA-4.0 | 2026-11-06 | 按热度排名的原子尺度机器学习项目列表（MLIP、描述符、数据集），便于追踪前沿工具。 |
| [tilde-lab/awesome-materials-informatics](https://github.com/tilde-lab/awesome-materials-informatics) | 535 | 未声明 | 2026-03-24 | 材料信息学 awesome 列表，收录 AiiDA、ASE、atomate、pymatgen 及各类材料数据库。 |
| [dstdev/awesome-hpc](https://github.com/dstdev/awesome-hpc) | 283 | CC0-1.0 | 2024-08-06 | HPC 软件与工具的 awesome 列表（调度器、MPI、性能分析等），适合延伸阅读。 |

### 脚本合集

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [accre/SLURM](https://github.com/accre/SLURM) | 73 | 未声明 | 2026-01-20 | Vanderbilt ACCRE 的 Slurm 示例脚本集（MPI、OpenMP、作业数组、GPU），可参考改写为自己的作业模板。 |
| [WangYiXu92/LOPOSTER](https://github.com/WangYiXu92/LOPOSTER) | 23 | 未声明 | 2025-09-23 | LOBSTER 输出的自动绘图后处理脚本，是 LobsterPy 之外的另一种选择。 |

### 工具 / 软件包

| 仓库 | ⭐ | 许可证 | 最近更新 | 说明 |
|---|---:|---|---|---|
| [ACEsuit/mace](https://github.com/ACEsuit/mace) | 1408 | 未声明 | 2026-10-28 | 主流的高阶等变机器学习原子间势（含 MACE-MP 等通用模型），可作为 ASE 计算器做结构预优化或 NEB 初猜。 |
| [nglviewer/nglview](https://github.com/nglviewer/nglview) | 931 | 未声明 | 2026-02-02 | Jupyter 中交互查看分子结构与轨迹的组件，可显示 ASE/pymatgen 结构，适合 notebook 形式的计算记录。 |
| [aiidateam/aiida-core](https://github.com/aiidateam/aiida-core) | 609 | NOASSERTION (LICENSE file: MIT) | 2026-11-09 | AiiDA 核心，带完整数据溯源（provenance）的工作流引擎；2.6 起可用 verdi presto 快速建库，不必安装 RabbitMQ。 |
| [ase/ase (GitLab)](https://gitlab.com/ase/ase) | 516 | LGPL-2.1-or-later | 2026-11-09 | 原子模拟 Python 库：统一的 Atoms 对象与 VASP/Gaussian/CP2K/OpenMX 计算器，含优化器、NEB（IDPP）、EOS、热化学、ase db；官方仓库在 GitLab。 |
| [spglib/spglib](https://github.com/spglib/spglib) | 378 | BSD-3-Clause | 2026-10-10 | 空间群识别、原胞/标准胞转换的底层 C 库，phonopy、pymatgen、seekpath 均依赖它，理解 symprec 设置时会用到。 |
| [janosh/matterviz](https://github.com/janosh/matterviz) | 370 | MIT | 2026-11-10 | 浏览器中交互展示晶体/分子、轨迹、凸包相图、能带+DOS、费米面、布里渊区，是 VESTA/OVITO 的轻量替代之一。 |
| [quanshengwu/wannier_tools](https://github.com/quanshengwu/wannier_tools) | 343 | GPL-3.0 | 2025-06-09 | 国内团队开发的拓扑材料软件，基于 Wannier90 的 _hr.dat 计算表面态、Wilson loop、Weyl 点等。 |
| [Quantum-Accelerators/quacc](https://github.com/Quantum-Accelerators/quacc) | 303 | BSD-3-Clause | 2026-11-09 | 基于 ASE 的“配方式”工作流平台，VASP/Gaussian/ORCA 等共用一套接口，可对接 Parsl/Dask/Prefect/jobflow 等引擎。 |
| [OpenChemistry/avogadroapp](https://github.com/OpenChemistry/avogadroapp) | 273 | BSD-3-Clause | 2026-11-09 | 免费开源分子编辑器，可生成 Gaussian/CP2K 等程序输入并显示轨迹和振动，是 GaussView 的开源替代。 |
| [materialsproject/api](https://github.com/materialsproject/api) | 180 | NOASSERTION (LBNL BSD-style) | 2026-11-06 | Materials Project 新版 API 客户端（mp-api）；旧版 API 已于 2025 年 9 月下线，必须换用新 API key。 |
| [materialscloud-org/seekpath](https://github.com/materialscloud-org/seekpath) | 169 | NOASSERTION (LICENSE file: MIT) | 2026-09-25 | 按 HPKOT 约定生成标准化高对称 k 路径并可视化布里渊区，可为不同 DFT 软件统一生成能带路径。 |
| [materialsproject/jobflow](https://github.com/materialsproject/jobflow) | 130 | NOASSERTION (LBNL BSD-style) | 2026-11-09 | atomate2 底层的工作流编写库，用 @job/Flow 把 Python 函数串成依赖图。 |
| [Z2PackDev/Z2Pack](https://github.com/Z2PackDev/Z2Pack) | 98 | GPL-3.0 | 2025-02-20 | 通过 Wannier 电荷中心计算 Z2/Chern 拓扑不变量，可直接调用 VASP 或使用紧束缚模型。 |
| [aiida-vasp/aiida-vasp](https://github.com/aiida-vasp/aiida-vasp) | 66 | NOASSERTION (LICENSE file: MIT) | 2026-11-10 | AiiDA 的 VASP 插件，提供 VaspWorkChain 及弛豫、能带、NEB、收敛测试等 workchain。 |
| [deepmodeling/dpdispatcher](https://github.com/deepmodeling/dpdispatcher) | 62 | LGPL-3.0 | 2026-11-09 | 用 Python 生成 Slurm/PBS/LSF 作业脚本并批量提交、监控直至完成，适合 HPC 上的大批量计算。 |
| [Matgenix/jobflow-remote](https://github.com/Matgenix/jobflow-remote) | 43 | 未声明 | 2026-11-06 | 把 jobflow/atomate2 工作流提交到 Slurm/PBS 集群运行（2026 年发布 1.0），是 FireWorks 的轻量替代。 |
| [henkelmangroup/bader](https://github.com/henkelmangroup/bader) | 1 | 未声明 | 2026-10-08 | Henkelman 组 Bader 电荷分析程序（官网 v1.05）的 GitHub 仓库，支持 VASP CHGCAR 和 cube 文件（Gaussian/CP2K 输出也可用）。 |

<a id="official-docs"></a>

## 官方文档与权威资料（Official documentation）

### VASP

- [VASP Wiki - The VASP Manual（官方手册总入口）](https://www.vasp.at/wiki/index.php/The_VASP_Manual)
- [VASP 官方 Tutorials（原子分子/体相/磁性/MD/MLFF/表面/过渡态/杂化/线性响应/GW/BSE/XAS/声子/电声）](https://www.vasp.at/tutorials/latest/)
- [VASP Changelog（6.6.1/6.6.0/6.5.x/6.4.x 新特性与修复）](https://www.vasp.at/wiki/index.php/Changelog)
- [VASP News（版本发布、AI工具使用规范 2026-06-16）](https://vasp.at/info/post/)
- [VASP Known issues（已知bug列表）](https://www.vasp.at/wiki/index.php/Known_issues)
- [Installing VASP.6.X.X（编译安装）](https://www.vasp.at/wiki/index.php/Installing_VASP.6.X.X)
- [makefile.include（oneAPI/GNU/NVHPC/AOCC模板与HDF5、Wannier90、DFT-D4、libMBD、VASPml选项）](https://www.vasp.at/wiki/index.php/Makefile.include)
- [GPU ports of VASP（OpenACC / OpenMP offload）](https://www.vasp.at/wiki/index.php/GPU_ports_of_VASP)
- [Optimizing the parallelization（NCORE/KPAR/NSIM）](https://www.vasp.at/wiki/index.php/Optimizing_the_parallelization)
- [Available PAW potentials（推荐POTCAR列表）](https://www.vasp.at/wiki/index.php/Available_PAW_potentials)
- [KPOINTS_OPT（6.3+ 非自洽能带新方式）](https://www.vasp.at/wiki/index.php/KPOINTS_OPT)
- [Band-structure calculation using hybrid functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals)
- [IVDW（DFT-D2/D3/D4/simple-DFT-D3/TS/MBD）](https://www.vasp.at/wiki/index.php/IVDW)
- [LDIPOL（偶极校正/功函数）](https://www.vasp.at/wiki/index.php/LDIPOL)
- [Nudged elastic bands（VASP内置NEB）](https://www.vasp.at/wiki/index.php/Nudged_elastic_bands)
- [IBRION（优化/MD/有限差分/DFPT）](https://www.vasp.at/wiki/index.php/IBRION)
- [Phonons from finite differences](https://www.vasp.at/wiki/index.php/Phonons_from_finite_differences)
- [LOPTICS（光学性质）](https://www.vasp.at/wiki/index.php/LOPTICS)
- [LWANNIER90（Wannier90接口）](https://www.vasp.at/wiki/index.php/LWANNIER90)
- [ML_MODE（机器学习力场 train/select/refit/run/delta）](https://www.vasp.at/wiki/index.php/ML_MODE)
- [Best practices for machine-learned force fields](https://www.vasp.at/wiki/index.php/Best_practices_for_machine-learned_force_fields)
- [Category: Machine-learned force fields](https://www.vasp.at/wiki/index.php/Category:Machine-learned_force_fields)
- [py4vasp 官方文档](https://vasp.at/py4vasp/latest/)
- [VTST Tools 安装（VTST补丁编译说明）](https://theory.cm.utexas.edu/vtsttools/installation.html)
- [VTST CI-NEB 说明](https://theory.cm.utexas.edu/vtsttools/neb.html)
- [VTST Dimer 方法](https://theory.cm.utexas.edu/vtsttools/dimer.html)
- [Bader Charge Analysis（Henkelman组）](https://theory.cm.utexas.edu/henkelman/code/bader/)
- [Phonopy - VASP interface](https://phonopy.github.io/phonopy/vasp.html)
- [Phonopy - VASP-DFPT](https://phonopy.github.io/phonopy/vasp-dfpt.html)
- [LOBSTER 官网（COHP）](http://www.cohp.de/)
- [LobsterPy 文档](https://digimatchem.github.io/LobsterPy/)
- [VASPKIT 官网文档（1.5.x）](https://vaspkit.com/)
- [VASPKIT 下载（SourceForge，最新1.5.1）](https://sourceforge.net/projects/vaspkit/)
- [pymatgen.io.vasp 文档](https://pymatgen.org/pymatgen.io.vasp.html)
- [ASE VASP calculator 文档](https://docs.ase-lib.org/ase/calculators/vasp.html)
- [sumo 文档](https://smtg-bham.github.io/sumo/)
- [Wannier90 官网](https://wannier.org/)

### Gaussian

- [Gaussian 16 Keyword List](https://gaussian.com/keywords/)
- [About Gaussian 16 Input (section ordering)](https://gaussian.com/input/)
- [Link 0 Commands (%chk, %oldchk, %mem, %nprocshared, %cpu, %gpucpu, %rwf, %nosave)](https://gaussian.com/link0/)
- [Running Gaussian (environment, scratch, Default.Route)](https://gaussian.com/running/)
- [g16 Command Line Options and environment-variable equivalents](https://gaussian.com/options/)
- [Density Functional (DFT) Methods, including the EmpiricalDispersion GD2/GD3/GD3BJ tables](https://gaussian.com/dft/)
- [MP Methods / double hybrids (B2PLYPD3, DSDPBEP86, PBE0DH)](https://gaussian.com/mp/)
- [Basis Sets](https://gaussian.com/basissets/)
- [Gen / GenECP / ExtraBasis](https://gaussian.com/gen/)
- [Pseudo (ECP input)](https://gaussian.com/pseudo/)
- [SCF keyword (XQC, QC, VShift, MaxCycle)](https://gaussian.com/scf/)
- [Guess keyword (Read, Mix, Fragment)](https://gaussian.com/guess/)
- [Integral keyword (UltraFine default grid in G16)](https://gaussian.com/integral/)
- [Geom keyword (Check, AllCheck, Connectivity, ModRedundant)](https://gaussian.com/geom/)
- [Opt keyword (TS, QST2/QST3, ModRedundant, CalcFC, MaxCycles, Restart)](https://gaussian.com/opt/)
- [Freq keyword (Raman, Anharmonic, ReadIsotopes, HPModes)](https://gaussian.com/freq/)
- [IRC keyword](https://gaussian.com/irc/)
- [QST2 Transition State Optimizations (white paper)](https://gaussian.com/qst2/)
- [Scan keyword (rigid scan)](https://gaussian.com/scan/)
- [SCRF keyword (PCM/IEFPCM/CPCM/SMD, solvents list, ExternalIteration)](https://gaussian.com/scrf/)
- [TD keyword (NStates, Root, Singlets/Triplets, 50-50)](https://gaussian.com/td/)
- [Creating UV/Visible Plots from Gaussian results](https://gaussian.com/uvvisplot/)
- [NMR keyword (GIAO, CSGT, SpinSpin)](https://gaussian.com/nmr/)
- [Population keyword (NBO/NBO7/NBORead, MK, CHelpG, Hirshfeld/CM5, RESP option)](https://gaussian.com/population/)
- [Stable keyword (wavefunction stability)](https://gaussian.com/stable/)
- [Antiferromagnetic Coupling (broken symmetry with Guess=Fragment)](https://gaussian.com/afc/)
- [Counterpoise keyword (BSSE)](https://gaussian.com/counterpoise/)
- [ONIOM keyword](https://gaussian.com/oniom/)
- [Thermochemistry in Gaussian (white paper)](https://gaussian.com/thermo/)
- [Vibrational Analysis in Gaussian (white paper)](https://gaussian.com/vib/)
- [Gaussian Utilities overview (formchk, unfchk, cubegen, chkchk, freqchk, newzmat, c8616)](https://gaussian.com/utils/)
- [formchk utility](https://gaussian.com/formchk/)
- [cubegen utility](https://gaussian.com/cubegen/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes (new features, GPU usage, G09 to G16 default changes, equivalencies)](https://gaussian.com/relnotes/)
- [GaussView 6](https://gaussian.com/gaussview6/)
- [GaussView 6.1.1 Release Notes](https://gaussian.com/gv611rn/)
- [Multiwfn home and download (latest 2026.10.18)](http://sobereva.com/multiwfn/)
- [Multiwfn online manual](http://sobereva.com/multiwfn/manual/)
- [Shermo thermochemistry program (Sobereva)](http://sobereva.com/soft/shermo/)
- [Basis Set Exchange (download basis sets/ECPs in Gaussian format for Gen/GenECP)](https://www.basissetexchange.org/)
- [NBO 7 official site](https://nbo7.chem.wisc.edu/)
- [cclib documentation](https://cclib.github.io/)
- [GoodVibes documentation](https://goodvibespy.readthedocs.io/)
- [CREST documentation](https://crest-lab.github.io/crest-docs/)
- [Solutions for Gaussian Error Messages (web version)](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/)
- [Keinsci (计算化学公社) Gaussian 报错原帖（中文）](http://bbs.keinsci.com/thread-4829-1-1.html)
- [The DFT Course (Cornell, Gaussian + SLURM)](https://thisisntnathan.github.io/dftCourse/)

### CP2K

- [CP2K Manual (trunk, sphinx docs + input reference)](https://manual.cp2k.org/trunk/)
- [CP2K 2026.2 Manual (versioned)](https://manual.cp2k.org/cp2k-2026_2-branch/index.html)
- [CP2K Input Reference (all sections/keywords)](https://manual.cp2k.org/trunk/CP2K_INPUT.html)
- [CP2K v2026.2 release notes (2026-07-15)](https://github.com/cp2k/cp2k/releases/tag/v2026.2)
- [CP2K v2026.1 release notes (2026-01-06, Makefile removed)](https://github.com/cp2k/cp2k/releases/tag/v2026.1)
- [Install from Distribution (conda-forge, Docker/Podman, apt, brew, easybuild ...)](https://manual.cp2k.org/trunk/getting-started/distributions.html)
- [Build from Source (CMake, toolchain + build_cp2k.sh, Spack via make_cp2k.sh)](https://manual.cp2k.org/trunk/getting-started/build-from-source.html)
- [Build with Spack](https://manual.cp2k.org/trunk/getting-started/build-with-spack.html)
- [First calculation](https://manual.cp2k.org/trunk/getting-started/first-calculation.html)
- [Troubleshooting (official error catalogue)](https://manual.cp2k.org/trunk/getting-started/troubleshooting.html)
- [How to make a SCF run converge](https://manual.cp2k.org/trunk/methods/dft/convergence.html)
- [How to Converge the CUTOFF and REL_CUTOFF](https://manual.cp2k.org/trunk/methods/dft/cutoff.html)
- [Pseudopotentials (GTH, POTENTIAL_UZH, ALL)](https://manual.cp2k.org/trunk/methods/dft/pseudopotentials.html)
- [Basis sets](https://manual.cp2k.org/trunk/methods/dft/basis_sets.html)
- [GAPW](https://manual.cp2k.org/trunk/methods/dft/gapw.html)
- [Orbital Transformation (OT)](https://manual.cp2k.org/trunk/methods/dft/orbital_transformation.html)
- [K-Points](https://manual.cp2k.org/trunk/methods/dft/k-points.html)
- [Hartree-Fock exchange / hybrid functionals / ADMM](https://manual.cp2k.org/trunk/methods/dft/hartree-fock/index.html)
- [Density of States (new &DOS/&PDOS/&CURVE interface, 2026.2+)](https://manual.cp2k.org/trunk/methods/electronic_structure/dos.html)
- [Geometry and Cell Optimization](https://manual.cp2k.org/trunk/methods/optimization/geometry_and_cell_opt.html)
- [Nudged Elastic Band](https://manual.cp2k.org/trunk/methods/optimization/nudged_elastic_band.html)
- [Molecular Dynamics](https://manual.cp2k.org/trunk/methods/sampling/molecular_dynamics.html)
- [Metadynamics](https://manual.cp2k.org/trunk/methods/sampling/metadynamics.html)
- [xTB (native GFN0/GFN1 and tblite GFN2)](https://manual.cp2k.org/trunk/methods/semiempiricals/xtb.html)
- [QM/MM](https://manual.cp2k.org/trunk/methods/qm_mm/index.html)
- [Restarting](https://manual.cp2k.org/trunk/methods/restarting.html)
- [CP2K Manual by Version (links to all versioned manuals)](https://manual.cp2k.org/trunk/versions.html)
- [CP2K docs source on GitHub (docs/)](https://github.com/cp2k/cp2k/tree/master/docs)
- [Legacy CP2K HowTos (static_calculation, converging_cutoff, geometry_optimisation, biochem_qmmm ...)](https://www.cp2k.org/howto)
- [CP2K Exercises (summer schools)](https://www.cp2k.org/exercises)
- [cp2k-basis tool (find matching basis set + pseudopotential)](https://www.cp2k.org/tools:cp2k-basis)
- [Docker Hub cp2k/cp2k images](https://hub.docker.com/r/cp2k/cp2k)
- [conda-forge cp2k package](https://anaconda.org/conda-forge/cp2k)
- [The CP2K Program Package Made Simple (J. Phys. Chem. B 2025)](https://doi.org/10.1021/acs.jpcb.5c05851)
- [ASE CP2K calculator documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/cp2k.html)
- [Multiwfn: create CP2K input files conveniently (sobereva.com/587)](http://sobereva.com/587)
- [Sobereva: key points of hybrid functional calculations in CP2K (sobereva.com/690)](http://sobereva.com/690)
- [CP2K Google Group (user forum)](https://groups.google.com/g/cp2k)

### OpenMX

- [OpenMX 官网首页](https://www.openmx-square.org/)
- [Download of OpenMX (4.0 + patch 4.0.1, 3.9 + patch 3.9.9, GPU, deb, Docker)](https://www.openmx-square.org/download.html)
- [What's new (4.0 新功能列表, patch 4.0.1, OpenMX Viewer 更新)](https://www.openmx-square.org/whatsnew/whatsnew.html)
- [Patch 4.0.1 README (bxsf/frmsf 单位、Mulliken HOMO 修复)](https://www.openmx-square.org/bugfixed/26May08/README.txt)
- [User's manual of OpenMX Ver. 4.0 (目录)](https://www.openmx-square.org/openmx_man4.0/index.html)
- [User's manual of OpenMX Ver. 3.9](https://www.openmx-square.org/openmx_man3.9/index.html)
- [Manual 4: Installation (依赖库、MPI/OpenMP、WSL、Debian 包)](https://www.openmx-square.org/openmx_man4.0/s4_installation.html)
- [Manual 4.11: Installation via a Debian package](https://www.openmx-square.org/openmx_man4.0/s4_11_deb-install.html)
- [Manual 8.2: Keywords (全部输入关键词)](https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html)
- [Manual 8.3: cif2omx (CIF 生成 dat)](https://www.openmx-square.org/openmx_man4.0/s8_3_gen-dat-files.html)
- [Manual 9: Output files](https://www.openmx-square.org/openmx_man4.0/s9_output-files.html)
- [Manual 12.3: Database Ver. 2019 推荐基组表 (Quick/Standard/Precise)](https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html)
- [Database (2019) of optimized VPS and PAO](https://www.openmx-square.org/vps_pao2019/)
- [Manual 15: k-point sampling (GCenter/Regular)](https://www.openmx-square.org/openmx_man4.0/s15_k-sampling.html)
- [Manual 16: SCF convergence](https://www.openmx-square.org/openmx_man4.0/s16_scf.html)
- [Manual 17: Restart (scf.restart on/off/pp/c2n)](https://www.openmx-square.org/openmx_man4.0/s17_restart.html)
- [Manual 18: Geometry optimization (EF/BFGS/RF/DIIS)](https://www.openmx-square.org/openmx_man4.0/s18_geom-opt.html)
- [Manual 19: Variable-cell optimization (OptC*/RFC*)](https://www.openmx-square.org/openmx_man4.0/s19_vcell-opt.html)
- [Manual 20: Molecular dynamics (NVE/NVT/NPT)](https://www.openmx-square.org/openmx_man4.0/s20_md.html)
- [Manual 22: Band structure (bandgnu13, 自动 k 路径)](https://www.openmx-square.org/openmx_man4.0/s22_band.html)
- [Manual 23: DOS/PDOS (DosMain)](https://www.openmx-square.org/openmx_man4.0/s23_dos.html)
- [Manual 24: COHP and COOP](https://www.openmx-square.org/openmx_man4.0/s24_cohp-coop.html)
- [Manual 28: Parallelization (MPI, -nt, scf.eigen.lib ELPA)](https://www.openmx-square.org/openmx_man4.0/s28_parallel.html)
- [Manual 30: Electric field](https://www.openmx-square.org/openmx_man4.0/s30_electric-field.html)
- [Manual 37: Non-collinear DFT](https://www.openmx-square.org/openmx_man4.0/s37_ncl-dft.html)
- [Manual 38: Relativistic effects (SOC)](https://www.openmx-square.org/openmx_man4.0/s38_rel-effects.html)
- [Manual 40: DFT+U](https://www.openmx-square.org/openmx_man4.0/s40_dft-plus-u.html)
- [Manual 44: Berry phase polarization (polB)](https://www.openmx-square.org/openmx_man4.0/s44_berry-polar.html)
- [Manual 45: Exchange coupling Jx](https://www.openmx-square.org/openmx_man4.0/s45_jx.html)
- [Manual 46: NEGF electronic transport](https://www.openmx-square.org/openmx_man4.0/s46_negf.html)
- [Manual 47: Wannier functions (CWF, MLWF)](https://www.openmx-square.org/openmx_man4.0/s47_wannier.html)
- [Manual 49: ESM method](https://www.openmx-square.org/openmx_man4.0/s49_esm.html)
- [Manual 50: Work function](https://www.openmx-square.org/openmx_man4.0/s50_work-function.html)
- [Manual 51: NEB](https://www.openmx-square.org/openmx_man4.0/s51_neb.html)
- [Manual 54: DFT-D2/D3](https://www.openmx-square.org/openmx_man4.0/s54_vdw-d2-d3.html)
- [Manual 55/56/58/59: Band unfolding, spin texture, Chern/Berry, Z2](https://www.openmx-square.org/openmx_man4.0/s58_chern-berry-fhs.html)
- [Manual 79: Known issues](https://www.openmx-square.org/openmx_man4.0/s79_known-issues.html)
- [OpenMX Viewer (浏览器端建模/k 路径/可视化)](https://www.openmx-square.org/viewer/)
- [OpenMX Forum](https://www.openmx-square.org/forum/patio.cgi)
- [Forum: How to Install OpenMX Ver. 4.0 on macOS (Apple Silicon)](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3524)
- [Forum: Docker image for OpenMX Ver. 4.0](https://www.openmx-square.org/forum/patio.cgi?mode=view&no=3522)
- [Installation of OpenMX on WSL (PDF)](https://www.openmx-square.org/Installation_OpenMX_WSL_EN.pdf)
- [Practical guide to OpenMX calculations (Ozaki lecture 2023, PDF)](https://t-ozaki.issp.u-tokyo.ac.jp/lectures-2023/OpenMX_Practical_Guide.pdf)
- [OpenMX workshops / hands-on 资料](https://www.openmx-square.org/workshop/workshop.html)
- [OpenMX video lectures](https://www.openmx-square.org/video_lec/video_lec.html)
- [TB2J OpenMX interface docs](https://tb2j.readthedocs.io/en/latest/src/openmx.html)
- [ASE OpenMX calculator source (GitLab)](https://gitlab.com/ase/ase/-/tree/master/ase/calculators/openmx)
- [pranabdas OpenMX tutorial site](https://pranabdas.github.io/openmx/)
- [Intel oneMKL Link Line Advisor (手册推荐用于 LIB 行)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)

### 通用工具与学习资源

- [ASE documentation (moved to docs.ase-lib.org, Jun 2026)](https://docs.ase-lib.org/)
- [ASE release notes](https://docs.ase-lib.org/releasenotes.html)
- [pymatgen documentation](https://pymatgen.org/)
- [pymatgen change log](https://pymatgen.org/CHANGES.html)
- [pymatgen.io.vasp API](https://pymatgen.org/pymatgen.io.vasp.html)
- [Materials Project API - Getting Started (new API key)](https://docs.materialsproject.org/downloading-data/using-the-api/getting-started)
- [VASPKIT 1.5 documentation](https://vaspkit.com/)
- [VASPKIT release notes](https://vaspkit.com/changelog.html)
- [Multiwfn official site / download](http://sobereva.com/multiwfn/)
- [Multiwfn manual](http://sobereva.com/multiwfn/Multiwfn_manual.html)
- [phonopy documentation](https://phonopy.github.io/phonopy/)
- [phono3py documentation](https://phonopy.github.io/phono3py/)
- [sumo documentation (tutorials)](https://smtg-bham.github.io/sumo/tutorials.html)
- [PyProcar documentation (v6)](https://pyprocar.readthedocs.io/en/stable/)
- [py4vasp documentation](https://vasp.at/py4vasp/latest)
- [VASP Wiki - vaspout.h5 (HDF5 output)](https://vasp.at/wiki/Vaspout.h5)
- [VASP Wiki - Changelog](https://vasp.at/wiki/index.php/Changelog)
- [VTST Tools download (GitHub-hosted docs)](https://henkelmangroup.github.io/vtsttools/download.html)
- [VTST Tools scripts reference](https://henkelmangroup.github.io/vtsttools/scripts.html)
- [Bader Charge Analysis (Henkelman group)](https://theory.cm.utexas.edu/henkelman/code/bader/)
- [LOBSTER official site](http://www.cohp.de)
- [LobsterPy documentation](https://jageo.github.io/LobsterPy/)
- [Wannier90 official site](https://wannier.org/)
- [Wannier90 releases](https://github.com/wannier-developers/wannier90/releases)
- [custodian documentation](https://materialsproject.github.io/custodian/)
- [atomate2 documentation](https://materialsproject.github.io/atomate2/)
- [atomate2 VASP user guide](https://materialsproject.github.io/atomate2/user/codes/vasp.html)
- [jobflow documentation](https://materialsproject.github.io/jobflow/)
- [jobflow-remote documentation](https://matgenix.github.io/jobflow-remote/)
- [AiiDA documentation](https://aiida.readthedocs.io/)
- [AiiDA plugin registry (aiida-vasp / aiida-cp2k / aiida-gaussian)](https://aiida.net/plugin-registry/)
- [AiiDA-VASP documentation](https://aiida-vasp.readthedocs.io/)
- [cclib documentation](https://cclib.github.io/)
- [GoodVibes documentation](https://goodvibespy.readthedocs.io/)
- [VESTA (download incl. 3.90.x dev builds)](https://jp-minerals.org/vesta/en/download.html)
- [OVITO user manual / changelog](https://docs.ovito.org/new_features.html)
- [Slurm sbatch manual](https://slurm.schedmd.com/sbatch.html)
- [spglib documentation](https://spglib.readthedocs.io/)
- [SeeK-path documentation](https://seekpath.readthedocs.io/)

---

欢迎通过 Issue / Pull Request 推荐更多优质仓库。返回 [项目首页](../README.md)。
