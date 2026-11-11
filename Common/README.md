# 通用基础：Linux、HPC、建模与后处理（Common DFT Skills）

四款软件都会用到的通用技能：Linux/WSL 与集群环境、Slurm/PBS/LSF 作业调度、结构数据库与建模、ASE/pymatgen、后处理与可视化、工作流自动化，以及 DFT 学习路线与机器学习势资源。

**适用版本**：与具体软件版本无关；示例基于 Python ≥ 3.9、ASE 3.2x、pymatgen 2025+。

> [!TIP]
> 建议先完成本板块的环境搭建，再进入 [VASP](../VASP/README.md)、[Gaussian](../Gaussian/README.md)、[CP2K](../CP2K/README.md) 或 [OpenMX](../OpenMX/README.md) 板块。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输（Setting up a DFT Computing Environment: Linux/WSL, conda, SSH & File Transfer）](tutorials/01-linux-hpc.md) | [01-env-setup](examples/01-env-setup/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-env-setup](examples/01-env-setup/README.md) | 示例 01：DFT 计算环境配置模板 |
| [02-array-chain-si](examples/02-array-chain-si/README.md) | 示例 02-array-chain-si |
| [02-job-templates](examples/02-job-templates/README.md) | 示例 02-job-templates |
| [03-defect-doping-mgo](examples/03-defect-doping-mgo/README.md) | 示例 03-defect-doping-mgo |
| [03-fetch-convert-si](examples/03-fetch-convert-si/README.md) | 示例 03-fetch-convert-si |
| [03-heterostructure-2d](examples/03-heterostructure-2d/README.md) | 示例 03-heterostructure-2d |
| [03-molecules](examples/03-molecules/README.md) | 示例 03-molecules |
| [03-slab-adsorption](examples/03-slab-adsorption/README.md) | 示例 03-slab-adsorption |
| [04-ase-calculators](examples/04-ase-calculators/README.md) | 04-ase-calculators：用 ASE 统一驱动 VASP、Gaussian、CP2K、OpenMX |
| [04-ase-emt-basics](examples/04-ase-emt-basics/README.md) | 04-ase-emt-basics：不装 DFT 软件也能跑的 ASE 基础练习 |
| [04-output-parsing](examples/04-output-parsing/README.md) | 04-output-parsing：跨程序输出解析、作业状态扫描、数据归档与 MLIP 数据准备 |
| [04-pymatgen-io](examples/04-pymatgen-io/README.md) | 04-pymatgen-io：pymatgen 输入集、VASP 输出解析、相图与 Gaussian/CP2K 输入 |
| [06-aiida-cp2k-vasp](examples/06-aiida-cp2k-vasp/README.md) | 示例 06-aiida-cp2k-vasp |
| [06-custodian-vasp-si](examples/06-custodian-vasp-si/README.md) | 示例 06-custodian-vasp-si |
| [06-jobflow-atomate2](examples/06-jobflow-atomate2/README.md) | 示例 06-jobflow-atomate2 |
| [06-quacc-dpdispatcher](examples/06-quacc-dpdispatcher/README.md) | 示例 06-quacc-dpdispatcher |

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
