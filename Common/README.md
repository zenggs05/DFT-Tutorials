# 通用基础：Linux、HPC、建模与后处理（Common DFT Skills）

四款软件都会用到的通用技能：Linux/WSL 与集群环境、Slurm/PBS/LSF 作业调度、结构数据库与建模、ASE/pymatgen、后处理与可视化、工作流自动化，以及 DFT 学习路线与机器学习势资源。

**适用版本**：与具体软件版本无关；示例基于 Python ≥ 3.9、ASE 3.2x、pymatgen 2025+。

> [!TIP]
> 建议先完成本板块的环境搭建，再进入 [VASP](../VASP/README.md)、[Gaussian](../Gaussian/README.md)、[CP2K](../CP2K/README.md) 或 [OpenMX](../OpenMX/README.md) 板块。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输（Setting up a DFT Computing Environment: Linux/WSL, conda, SSH & File Transfer）](tutorials/01-linux-hpc.md) | [01-env-setup](examples/01-env-setup/README.md) |
| [第2章 HPC 作业调度：Slurm/PBS/LSF 模板与并行设置（HPC Job Scheduling: Slurm/PBS/LSF Templates and Parallel Settings）](tutorials/02-job-scheduler.md) | [02-array-chain-si](examples/02-array-chain-si/README.md)<br>[02-job-templates](examples/02-job-templates/README.md) |
| [第3章 结构准备：数据库下载、建模（超胞/表面/缺陷/掺杂/异质结）与格式转换（Structure Preparation: Databases, Supercells, Slabs, Defects, Doping, Heterostructures & Format Conversion）](tutorials/03-structure-modeling.md) | [03-defect-doping-mgo](examples/03-defect-doping-mgo/README.md)<br>[03-fetch-convert-si](examples/03-fetch-convert-si/README.md)<br>[03-heterostructure-2d](examples/03-heterostructure-2d/README.md)<br>[03-molecules](examples/03-molecules/README.md)<br>[03-slab-adsorption](examples/03-slab-adsorption/README.md) |
| [第4章 ASE 与 pymatgen 实战：统一驱动四款软件、输入集与输出解析（ASE & pymatgen in Practice: Driving VASP, Gaussian, CP2K and OpenMX, Input Sets and Output Parsing）](tutorials/04-ase-pymatgen.md) | [04-ase-calculators](examples/04-ase-calculators/README.md)<br>[04-ase-emt-basics](examples/04-ase-emt-basics/README.md)<br>[04-output-parsing](examples/04-output-parsing/README.md)<br>[04-pymatgen-io](examples/04-pymatgen-io/README.md) |
| [第5章 后处理与可视化：VASPKIT、sumo、pyprocar、VESTA、Multiwfn 与出版级作图（Post-processing & Visualization: VASPKIT, sumo, PyProcar, py4vasp, VESTA, OVITO, Multiwfn and Publication-quality Figures）](tutorials/05-postprocess-visualization.md) | [05-band-dos-plotting](examples/05-band-dos-plotting/README.md)<br>[05-vaspkit-batch](examples/05-vaspkit-batch/README.md)<br>[05-visualization-3d](examples/05-visualization-3d/README.md) |
| [第6章 工作流自动化与高通量：custodian、atomate2、jobflow、AiiDA（Workflow Automation & High-Throughput: custodian, atomate2, jobflow, AiiDA）](tutorials/06-workflow-automation.md) | [06-aiida-cp2k-vasp](examples/06-aiida-cp2k-vasp/README.md)<br>[06-custodian-vasp-si](examples/06-custodian-vasp-si/README.md)<br>[06-jobflow-atomate2](examples/06-jobflow-atomate2/README.md)<br>[06-quacc-dpdispatcher](examples/06-quacc-dpdispatcher/README.md) |
| [第7章 DFT 理论速览、学习路线与资源汇总（含机器学习势前沿）（DFT Primer, Learning Path, Convergence Checklist, Resources and Machine-Learned Potentials）](tutorials/07-learning-path.md) | [07-convergence-templates](examples/07-convergence-templates/README.md)<br>[07-learning-resources](examples/07-learning-resources/README.md)<br>[07-mlip-nep-mace](examples/07-mlip-nep-mace/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-env-setup](examples/01-env-setup/README.md) | 示例 01：DFT 计算环境配置模板 |
| [02-array-chain-si](examples/02-array-chain-si/README.md) | 示例 02-array-chain-si：用作业数组做 ENCUT 扫描、用依赖链串联弛豫→静态→能带 |
| [02-job-templates](examples/02-job-templates/README.md) | 示例 02-job-templates：Slurm/PBS/LSF 作业脚本模板 |
| [03-defect-doping-mgo](examples/03-defect-doping-mgo/README.md) | 示例 03-defect-doping-mgo：MgO 超胞中的空位、替位掺杂、间隙原子与双掺杂构型 |
| [03-fetch-convert-si](examples/03-fetch-convert-si/README.md) | 示例 03-fetch-convert-si：下载 Si 结构、检查对称性并转换为四款软件的输入 |
| [03-heterostructure-2d](examples/03-heterostructure-2d/README.md) | 示例 03-heterostructure-2d：二维范德华异质结的晶格匹配与搭建 |
| [03-molecules](examples/03-molecules/README.md) | 示例 03-molecules：分子建模 |
| [03-slab-adsorption](examples/03-slab-adsorption/README.md) | 示例 03-slab-adsorption：Cu |
| [04-ase-calculators](examples/04-ase-calculators/README.md) | 04-ase-calculators：用 ASE 统一驱动 VASP、Gaussian、CP2K、OpenMX |
| [04-ase-emt-basics](examples/04-ase-emt-basics/README.md) | 04-ase-emt-basics：不装 DFT 软件也能跑的 ASE 基础练习 |
| [04-output-parsing](examples/04-output-parsing/README.md) | 04-output-parsing：跨程序输出解析、作业状态扫描、数据归档与 MLIP 数据准备 |
| [04-pymatgen-io](examples/04-pymatgen-io/README.md) | 04-pymatgen-io：pymatgen 输入集、VASP 输出解析、相图与 Gaussian/CP2K 输入 |
| [05-band-dos-plotting](examples/05-band-dos-plotting/README.md) | 05-band-dos-plotting：GaAs 能带与态密度——从 VASP 输入到出版级作图 |
| [05-vaspkit-batch](examples/05-vaspkit-batch/README.md) | 05-vaspkit-batch：VASPKIT 1.5.x 非交互批处理与常用后处理输入 |
| [05-visualization-3d](examples/05-visualization-3d/README.md) | 05-visualization-3d：结构、体数据、轨迹与振动的可视化 |
| [06-aiida-cp2k-vasp](examples/06-aiida-cp2k-vasp/README.md) | 06-aiida-cp2k-vasp：AiiDA 驱动 CP2K、VASP、Gaussian |
| [06-custodian-vasp-si](examples/06-custodian-vasp-si/README.md) | 06-custodian-vasp-si：用 custodian 自动纠错运行 VASP |
| [06-jobflow-atomate2](examples/06-jobflow-atomate2/README.md) | 06-jobflow-atomate2：jobflow 入门、atomate2 VASP 工作流与 jobflow-remote 集群提交 |
| [06-quacc-dpdispatcher](examples/06-quacc-dpdispatcher/README.md) | 06-quacc-dpdispatcher：quacc 配方与 dpdispatcher 轻量批量提交 |
| [07-convergence-templates](examples/07-convergence-templates/README.md) | 示例 07-convergence-templates：四款软件通用的收敛测试模板 |
| [07-learning-resources](examples/07-learning-resources/README.md) | 示例 07-learning-resources：资源书签与学习打卡表 |
| [07-mlip-nep-mace](examples/07-mlip-nep-mace/README.md) | 示例 07-mlip-nep-mace：机器学习势预优化、NEB 初猜与 NEP/MACE/DeePMD 最小训练配置 |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。

- [convert_any.py](scripts/convert_any.py)
- [env_check.sh](scripts/env_check.sh)
- [job_status.py](scripts/job_status.py)
- [kpath_gen.py](scripts/kpath_gen.py)
- [mkjob.py](scripts/mkjob.py)

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
