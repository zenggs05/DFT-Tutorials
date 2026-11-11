# Gaussian 中文教程：从入门到进阶（Gaussian 16 Tutorial in Chinese）

系统讲解 Gaussian 16 量子化学计算：输入文件结构、泛函与基组选择、SCF 收敛、几何优化、频率与热化学、过渡态与 IRC、势能面扫描、PCM/SMD 溶剂模型、TD-DFT 激发态（UV-Vis/荧光/ECD）、NMR、NBO 与 RESP 电荷、开壳层与 Broken-Symmetry、BSSE、ONIOM、Multiwfn 后处理、作业脚本与常见报错。

**适用版本**：Gaussian 16 Rev. C.01/C.02 + GaussView 6。

> [!TIP]
> 建议先阅读 [通用基础](../Common/README.md) 了解 Linux 与集群作业提交，再按章节顺序学习；每章末尾都有参考资料与常见报错。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 环境准备与第一个 Gaussian 计算（Getting Started with Gaussian 16 and GaussView 6）](tutorials/01-setup-first-job.md) | [01-first-job-h2o](examples/01-first-job-h2o/README.md) |
| [第2章 输入文件结构详解（Gaussian Input File Anatomy: Link 0, Route, Molecule Specification and Link1）](tutorials/02-input-anatomy.md) | [02-input-anatomy](examples/02-input-anatomy/README.md) |
| [第3章 方法、泛函与基组选择（含 Gen/GenECP 与赝势）（Choosing DFT Functionals, Dispersion Corrections, Basis Sets and ECPs in Gaussian 16）](tutorials/03-method-basis.md) | [03-method-basis](examples/03-method-basis/README.md) |
| [第4章 单点能计算与 SCF 收敛技巧（Single-Point Energy & SCF Convergence in Gaussian 16）](tutorials/04-sp-scf.md) | [04-sp-scf](examples/04-sp-scf/README.md) |
| [第5章 几何优化（Geometry Optimization with Gaussian 16）](tutorials/05-opt.md) | [05-opt](examples/05-opt/README.md) |
| [第6章 频率计算、红外拉曼与热化学（Frequencies, IR/Raman & Thermochemistry in Gaussian 16）](tutorials/06-freq-thermo.md) | [06-freq-thermo](examples/06-freq-thermo/README.md) |
| [第7章 过渡态搜索与 IRC 验证（Transition State Search & IRC with Gaussian 16）](tutorials/07-ts-irc.md) | [07-qst-hcn-hnc](examples/07-qst-hcn-hnc/README.md)<br>[07-ts-irc-sn2](examples/07-ts-irc-sn2/README.md) |
| [第8章 势能面扫描（刚性与柔性扫描）（Rigid & Relaxed PES Scans in Gaussian 16）](tutorials/08-scan.md) | [08-pes-scan](examples/08-pes-scan/README.md) |
| [第9章 隐式溶剂模型 PCM/SMD（Implicit Solvation with PCM & SMD in Gaussian 16）](tutorials/09-solvation.md) | [09-solvation-smd](examples/09-solvation-smd/README.md) |
| [第10章 激发态 TD-DFT：UV-Vis、荧光发射与 ECD（Excited States with Gaussian TD-DFT: UV-Vis, Fluorescence & ECD）](tutorials/10-tddft.md) | [10-ecd-methyloxirane](examples/10-ecd-methyloxirane/README.md)<br>[10-nto-dct-dmabn](examples/10-nto-dct-dmabn/README.md)<br>[10-uvvis-emission-coumarin](examples/10-uvvis-emission-coumarin/README.md) |
| [第11章 NMR 化学位移与其他分子性质（NMR Chemical Shifts, J-Couplings, Polarizability & Molecular Properties with Gaussian）](tutorials/11-nmr-props.md) | [11-nmr-jcoupling-ethanol](examples/11-nmr-jcoupling-ethanol/README.md)<br>[11-nmr-shift-camphor](examples/11-nmr-shift-camphor/README.md)<br>[11-polar-hyperpol-pna](examples/11-polar-hyperpol-pna/README.md) |
| [第12章 布居分析、NBO 与静电势/RESP 电荷（Population Analysis, NBO & ESP/RESP Charges with Gaussian）](tutorials/12-population-charges.md) | [12-esp-resp-paracetamol](examples/12-esp-resp-paracetamol/README.md)<br>[12-hirshfeld-cm5-acetate](examples/12-hirshfeld-cm5-acetate/README.md)<br>[12-nbo-pyridine-nma](examples/12-nbo-pyridine-nma/README.md) |
| [第13章 开壳层体系、波函数稳定性与 Broken-Symmetry（Open-Shell Systems, Wavefunction Stability & Broken-Symmetry DFT in Gaussian）](tutorials/13-openshell-stability.md) | [13-bs-cu2-dimer](examples/13-bs-cu2-dimer/README.md)<br>[13-diradical-pbenzyne](examples/13-diradical-pbenzyne/README.md)<br>[13-spin-states-o2-ch2](examples/13-spin-states-o2-ch2/README.md) |
| [第14章 弱相互作用与 BSSE（Counterpoise）（Noncovalent Interactions & Counterpoise BSSE Correction in Gaussian）](tutorials/14-weak-bsse.md) | [14-cp-benzene-dimer](examples/14-cp-benzene-dimer/README.md)<br>[14-cp-crown-k](examples/14-cp-crown-k/README.md)<br>[14-cp-water-dimer](examples/14-cp-water-dimer/README.md) |
| [第15章 ONIOM 多层方法与 QM/MM（ONIOM Multilayer Methods & QM/MM in Gaussian）](tutorials/15-oniom.md) | [15-oniom-hexanal](examples/15-oniom-hexanal/README.md)<br>[15-oniom3-cyclohexanone](examples/15-oniom3-cyclohexanone/README.md)<br>[15-qmmm-phenol-water](examples/15-qmmm-phenol-water/README.md) |
| [第16章 检查点文件、cubegen 与 Multiwfn 后处理（Checkpoint Files, cubegen & Multiwfn Post-processing）](tutorials/16-postprocess-multiwfn.md) | [16-postprocess-multiwfn](examples/16-postprocess-multiwfn/README.md) |
| [第17章 作业脚本、并行与批量计算（Running Gaussian 16 on HPC: Job Scripts, Parallelism & Batch Jobs）](tutorials/17-hpc-jobs.md) | [17-hpc-job-scripts](examples/17-hpc-job-scripts/README.md) |
| [第18章 常见报错排查大全与完整研究流程案例（Gaussian Error Messages & a Complete Research Workflow）](tutorials/18-errors-capstone.md) | [18-capstone-diels-alder](examples/18-capstone-diels-alder/README.md)<br>[18-error-cases](examples/18-error-cases/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-first-job-h2o](examples/01-first-job-h2o/README.md) | 示例 01：第一个 Gaussian 计算——水分子单点与几何优化 |
| [02-input-anatomy](examples/02-input-anatomy/README.md) | 示例 02：输入文件结构——笛卡尔坐标、Z-matrix、Link1、片段标签与连接表 |
| [03-method-basis](examples/03-method-basis/README.md) | 示例 03：方法、泛函、色散校正与混合基组 |
| [04-sp-scf](examples/04-sp-scf/README.md) | 示例 04：单点能计算与 SCF 收敛技巧 |
| [05-opt](examples/05-opt/README.md) | 示例 05：几何优化 |
| [06-freq-thermo](examples/06-freq-thermo/README.md) | 示例 06：频率计算、红外拉曼与热化学 |
| [07-qst-hcn-hnc](examples/07-qst-hcn-hnc/README.md) | 示例 07：用 QST2 / QST3 搜索 HCN ⇌ HNC 异构化过渡态 |
| [07-ts-irc-sn2](examples/07-ts-irc-sn2/README.md) | 示例 07：Cl⁻ + CH₃Br S<sub>N</sub>2 反应的 TS 优化、IRC 与端点验证 |
| [08-pes-scan](examples/08-pes-scan/README.md) | 示例 08：刚性扫描、柔性扫描、二维扫描与 GIC 扫描 |
| [09-solvation-smd](examples/09-solvation-smd/README.md) | 示例 09：SMD 隐式溶剂——溶液中优化、溶剂化自由能、自定义溶剂、pKa 与氧化还原电位 |
| [10-ecd-methyloxirane](examples/10-ecd-methyloxirane/README.md) | 示例 10-ecd-methyloxirane： |
| [10-nto-dct-dmabn](examples/10-nto-dct-dmabn/README.md) | 示例 10-nto-dct-dmabn：DMABN 激发态的 NTO、DCT 与空穴-电子分析 |
| [10-uvvis-emission-coumarin](examples/10-uvvis-emission-coumarin/README.md) | 示例 10-uvvis-emission-coumarin：香豆素的吸收、荧光与磷光 |
| [11-nmr-jcoupling-ethanol](examples/11-nmr-jcoupling-ethanol/README.md) | 示例 11-nmr-jcoupling-ethanol：乙醇的自旋-自旋耦合常数 |
| [11-nmr-shift-camphor](examples/11-nmr-shift-camphor/README.md) | 示例 11-nmr-shift-camphor：TMS 参比与樟脑的 1H/13C 化学位移 |
| [11-polar-hyperpol-pna](examples/11-polar-hyperpol-pna/README.md) | 示例 11-polar-hyperpol-pna：p-硝基苯胺的极化率与第一超极化率 |
| [12-esp-resp-paracetamol](examples/12-esp-resp-paracetamol/README.md) | 示例 12-esp-resp-paracetamol：小分子配体的 MK 静电势与 RESP 电荷 |
| [12-hirshfeld-cm5-acetate](examples/12-hirshfeld-cm5-acetate/README.md) | 示例 12-hirshfeld-cm5-acetate：醋酸根的 Hirshfeld/CM5、MK 与 NPA 电荷对比 |
| [12-nbo-pyridine-nma](examples/12-nbo-pyridine-nma/README.md) | 示例 12-nbo-pyridine-nma：Mulliken 基组依赖性、NBO、Wiberg 键级与 E |
| [13-bs-cu2-dimer](examples/13-bs-cu2-dimer/README.md) | 示例 13-3：双核 Cu |
| [13-diradical-pbenzyne](examples/13-diradical-pbenzyne/README.md) | 示例 13-2：对苯炔单重态双自由基——Stable=Opt、Guess=Mix 与 UNO |
| [13-spin-states-o2-ch2](examples/13-spin-states-o2-ch2/README.md) | 示例 13-1：O2 与 CH2 的自旋多重度比较 |
| [14-cp-benzene-dimer](examples/14-cp-benzene-dimer/README.md) | 示例 14-2：苯二聚体的 CP 校正优化 |
| [14-cp-crown-k](examples/14-cp-crown-k/README.md) | 示例 14-3：主客体复合物 18-冠-6·K⁺ 的 CP 相互作用能与形变能 |
| [14-cp-water-dimer](examples/14-cp-water-dimer/README.md) | 示例 14-1：水二聚体的 Counterpoise |
| [15-oniom-hexanal](examples/15-oniom-hexanal/README.md) | 示例 15-1：二层 ONIOM——己醛的 ONIOM |
| [15-oniom3-cyclohexanone](examples/15-oniom3-cyclohexanone/README.md) | 示例 15-3：三层 ONIOM |
| [15-qmmm-phenol-water](examples/15-qmmm-phenol-water/README.md) | 示例 15-2：QM/MM 电子嵌入——苯酚与 6 个 TIP3P 水 |
| [16-postprocess-multiwfn](examples/16-postprocess-multiwfn/README.md) | 示例 16：formchk、cubegen 与 Multiwfn 后处理 |
| [17-hpc-job-scripts](examples/17-hpc-job-scripts/README.md) | 示例 17：Gaussian 16 作业脚本、并行、批量提交与断点续算 |
| [18-capstone-diels-alder](examples/18-capstone-diels-alder/README.md) | 示例 18b：完整研究流程——丁二烯 + 乙烯 Diels–Alder 反应的自由能剖面 |
| [18-error-cases](examples/18-error-cases/README.md) | 示例 18a：故意出错的 Gaussian 输入与修正版本 |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。

- [g16_log_summary.py](scripts/g16_log_summary.py)
- [gjf_builder.py](scripts/gjf_builder.py)
- [log2xyz.py](scripts/log2xyz.py)
- [opt_monitor.py](scripts/opt_monitor.py)
- [spectrum_ir_raman.py](scripts/spectrum_ir_raman.py)
- [spectrum_uvvis_ecd.py](scripts/spectrum_uvvis_ecd.py)
- [thermo_qrrho.py](scripts/thermo_qrrho.py)

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
