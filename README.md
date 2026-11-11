[简体中文](README.md) | [English](README_en.md) | [日本語](README_ja.md) | [Português](README_pt.md)

# DFT-Tutorials：第一性原理计算中文教程与脚本库

**VASP · Gaussian · CP2K · OpenMX 从入门到进阶 ｜ 可直接运行的示例输入 ｜ 前后处理脚本 ｜ GitHub 资源精选**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![Docs](https://img.shields.io/badge/docs-online-blue.svg)](https://zenggs05.github.io/DFT-Tutorials/) [![GitHub stars](https://img.shields.io/github/stars/zenggs05/DFT-Tutorials?style=social)](https://github.com/zenggs05/DFT-Tutorials/stargazers) [![Last commit](https://img.shields.io/github/last-commit/zenggs05/DFT-Tutorials)](https://github.com/zenggs05/DFT-Tutorials/commits/main)

本仓库是一套面向科研新手与进阶用户的**第一性原理（DFT）计算中文教程**，系统讲解 **VASP、Gaussian、CP2K、OpenMX** 四款主流软件：从安装编译、输入文件、收敛测试、结构优化，到能带/态密度、杂化泛函、DFT+U、自旋轨道耦合、表面与吸附、过渡态（NEB/IRC）、声子、分子动力学、激发态（TD-DFT）、溶剂化、NBO 与电荷分析、电子输运等。每章都配有可直接运行的示例输入、逐参数解释、常见报错排查与参考资料，并提供配套的 Python/Bash 脚本。

📖 **在线阅读（带搜索的文档网站）**：<https://zenggs05.github.io/DFT-Tutorials/>

目前收录 **65 章教程**、**152 套示例输入**、**44 个脚本**，以及 **230+ 个相关 GitHub 仓库**的精选索引。

## 📚 内容一览

| 板块 | 定位 | 入口 |
|---|---|---|
| **VASP** | 平面波 PAW，周期性固体/表面/催化 | [VASP/README.md](VASP/README.md) |
| **Gaussian** | 高斯基组量子化学，分子结构/反应/光谱 | [Gaussian/README.md](Gaussian/README.md) |
| **CP2K** | GPW 混合基组，大体系 AIMD / QM/MM | [CP2K/README.md](CP2K/README.md) |
| **OpenMX** | 数值原子轨道 + 模守恒赝势，磁性/输运/O(N) | [OpenMX/README.md](OpenMX/README.md) |
| **通用基础** | Linux/HPC、建模、ASE/pymatgen、可视化、学习路线 | [Common/README.md](Common/README.md) |
| **NEP/GPUMD 脚本** | 机器学习势训练集生成、NEP 结果分析、GPUMD 热输运作图 | [Scripts/README.md](Scripts/README.md) |
| **资源精选** | 230+ 个 VASP/Gaussian/CP2K/OpenMX 相关 GitHub 仓库与官方文档 | [Resources/github-repos.md](Resources/github-repos.md) |

## 🧭 教程目录

### VASP

| 章节 | 示例 |
|---|---|
| [第1章 环境准备、VASP 编译安装（CPU/GPU）与作业提交](VASP/tutorials/01-setup-compile.md) | [01-setup-compile](VASP/examples/01-setup-compile/README.md) |
| [第2章 四大输入文件与输出文件全解析](VASP/tutorials/02-input-output-files.md) | [02-input-files-si](VASP/examples/02-input-files-si/README.md) |
| [第3章 收敛测试与计算精度控制](VASP/tutorials/03-convergence.md) | [03-convergence-si-al](VASP/examples/03-convergence-si-al/README.md) |
| [第4章 结构优化与状态方程](VASP/tutorials/04-relax-eos.md) | [04-eos-si](VASP/examples/04-eos-si/README.md) [04-relax-co-molecule](VASP/examples/04-relax-co-molecule/README.md) [04-relax-mos2-2d](VASP/examples/04-relax-mos2-2d/README.md) [04-relax-si-isif3](VASP/examples/04-relax-si-isif3/README.md) |
| [第5章 能带结构与态密度（含 KPOINTS_OPT）](VASP/tutorials/05-band-dos.md) | [05-band-dos-si](VASP/examples/05-band-dos-si/README.md) |
| [第6章 杂化泛函与 meta-GGA（HSE06 / r2SCAN / mBJ）](VASP/tutorials/06-hybrid-metagga.md) | [06-hse06-band-si](VASP/examples/06-hse06-band-si/README.md) [06-mbj-gap-si](VASP/examples/06-mbj-gap-si/README.md) [06-r2scan-relax-zno](VASP/examples/06-r2scan-relax-zno/README.md) |
| [第7章 磁性、DFT+U 与自旋轨道耦合](VASP/tutorials/07-magnetism-u-soc.md) | [07-fe-bcc-fm](VASP/examples/07-fe-bcc-fm/README.md) [07-fe-bcc-soc-mae](VASP/examples/07-fe-bcc-soc-mae/README.md) [07-linear-response-u](VASP/examples/07-linear-response-u/README.md) [07-nio-afm-dftu](VASP/examples/07-nio-afm-dftu/README.md) |
| [第8章 范德华修正、二维与层状材料](VASP/tutorials/08-vdw-2d.md) | [08-graphite-vdw](VASP/examples/08-graphite-vdw/README.md) [08-mos2-bilayer-binding](VASP/examples/08-mos2-bilayer-binding/README.md) [08-mos2-monolayer](VASP/examples/08-mos2-monolayer/README.md) |
| [第9章 表面模型、表面能、功函数与偶极校正](VASP/tutorials/09-surface-workfunction.md) | [09-cu100-surface-energy](VASP/examples/09-cu100-surface-energy/README.md) [09-pt111-slab](VASP/examples/09-pt111-slab/README.md) [09-pt111-workfunction-dipole](VASP/examples/09-pt111-workfunction-dipole/README.md) |
| [第10章 吸附、催化与自由能（含隐式溶剂）](VASP/tutorials/10-adsorption-catalysis.md) | [10-co-pt111-adsorption](VASP/examples/10-co-pt111-adsorption/README.md) [10-her-oer-free-energy](VASP/examples/10-her-oer-free-energy/README.md) [10-vaspsol-implicit-solvent](VASP/examples/10-vaspsol-implicit-solvent/README.md) |
| [第11章 过渡态搜索：NEB / CI-NEB / Dimer（VTST）](VASP/tutorials/11-neb-dimer.md) | [11-dimer-from-neb](VASP/examples/11-dimer-from-neb/README.md) [11-li-si-neb](VASP/examples/11-li-si-neb/README.md) [11-pt-adatom-neb](VASP/examples/11-pt-adatom-neb/README.md) [11-ts-freq-check](VASP/examples/11-ts-freq-check/README.md) |
| [第12章 声子：phonopy 有限位移与 DFPT](VASP/tutorials/12-phonon.md) | [12-graphene-phonon](VASP/examples/12-graphene-phonon/README.md) [12-mgo-dfpt-nac](VASP/examples/12-mgo-dfpt-nac/README.md) [12-si-phonopy-fd](VASP/examples/12-si-phonopy-fd/README.md) |
| [第13章 弹性常数与力学性质](VASP/tutorials/13-elastic.md) | [13-elastic-al-ibrion6](VASP/examples/13-elastic-al-ibrion6/README.md) [13-elastic-mos2-2d](VASP/examples/13-elastic-mos2-2d/README.md) |
| [第14章 介电与光学性质](VASP/tutorials/14-optics-dielectric.md) | [14-dielectric-gaas-lepsilon](VASP/examples/14-dielectric-gaas-lepsilon/README.md) [14-dielectric-zno-hse-lcalceps](VASP/examples/14-dielectric-zno-hse-lcalceps/README.md) [14-optics-si-loptics](VASP/examples/14-optics-si-loptics/README.md) |
| [第15章 从头算分子动力学与机器学习力场](VASP/tutorials/15-md-mlff.md) | [15-aimd-si-liquid-nvt](VASP/examples/15-aimd-si-liquid-nvt/README.md) [15-aimd-water-npt](VASP/examples/15-aimd-water-npt/README.md) [15-mlff-si-train-run](VASP/examples/15-mlff-si-train-run/README.md) |
| [第16章 电荷与成键分析：Bader、差分电荷、ELF、STM、COHP](VASP/tutorials/16-charge-bonding.md) | [16-bader-co-pt111](VASP/examples/16-bader-co-pt111/README.md) [16-cdd-co-pt111](VASP/examples/16-cdd-co-pt111/README.md) [16-elf-si](VASP/examples/16-elf-si/README.md) [16-lobster-gaas](VASP/examples/16-lobster-gaas/README.md) [16-stm-mos2](VASP/examples/16-stm-mos2/README.md) |
| [第17章 高级电子结构：Wannier90、有效质量、能带反折叠、费米面与拓扑](VASP/tutorials/17-advanced-electronic.md) | [17-effmass-gaas](VASP/examples/17-effmass-gaas/README.md) [17-fermi-surface-cu](VASP/examples/17-fermi-surface-cu/README.md) [17-unfold-si-vacancy](VASP/examples/17-unfold-si-vacancy/README.md) [17-wannier90-si](VASP/examples/17-wannier90-si/README.md) |
| [第18章 后处理生态、自动化与报错排查](VASP/tutorials/18-postprocess-troubleshooting.md) | [18-postprocess-toolkit](VASP/examples/18-postprocess-toolkit/README.md) |

### Gaussian

| 章节 | 示例 |
|---|---|
| [第1章 环境准备与第一个 Gaussian 计算](Gaussian/tutorials/01-setup-first-job.md) | [01-first-job-h2o](Gaussian/examples/01-first-job-h2o/README.md) |
| [第2章 输入文件结构详解](Gaussian/tutorials/02-input-anatomy.md) | [02-input-anatomy](Gaussian/examples/02-input-anatomy/README.md) |
| [第3章 方法、泛函与基组选择（含 Gen/GenECP 与赝势）](Gaussian/tutorials/03-method-basis.md) | [03-method-basis](Gaussian/examples/03-method-basis/README.md) |
| [第4章 单点能计算与 SCF 收敛技巧](Gaussian/tutorials/04-sp-scf.md) | [04-sp-scf](Gaussian/examples/04-sp-scf/README.md) |
| [第5章 几何优化](Gaussian/tutorials/05-opt.md) | [05-opt](Gaussian/examples/05-opt/README.md) |
| [第6章 频率计算、红外拉曼与热化学](Gaussian/tutorials/06-freq-thermo.md) | [06-freq-thermo](Gaussian/examples/06-freq-thermo/README.md) |
| [第7章 过渡态搜索与 IRC 验证](Gaussian/tutorials/07-ts-irc.md) | [07-qst-hcn-hnc](Gaussian/examples/07-qst-hcn-hnc/README.md) [07-ts-irc-sn2](Gaussian/examples/07-ts-irc-sn2/README.md) |
| [第8章 势能面扫描（刚性与柔性扫描）](Gaussian/tutorials/08-scan.md) | [08-pes-scan](Gaussian/examples/08-pes-scan/README.md) |
| [第9章 隐式溶剂模型 PCM/SMD](Gaussian/tutorials/09-solvation.md) | [09-solvation-smd](Gaussian/examples/09-solvation-smd/README.md) |
| [第10章 激发态 TD-DFT：UV-Vis、荧光发射与 ECD](Gaussian/tutorials/10-tddft.md) | [10-ecd-methyloxirane](Gaussian/examples/10-ecd-methyloxirane/README.md) [10-nto-dct-dmabn](Gaussian/examples/10-nto-dct-dmabn/README.md) [10-uvvis-emission-coumarin](Gaussian/examples/10-uvvis-emission-coumarin/README.md) |
| [第11章 NMR 化学位移与其他分子性质](Gaussian/tutorials/11-nmr-props.md) | [11-nmr-jcoupling-ethanol](Gaussian/examples/11-nmr-jcoupling-ethanol/README.md) [11-nmr-shift-camphor](Gaussian/examples/11-nmr-shift-camphor/README.md) [11-polar-hyperpol-pna](Gaussian/examples/11-polar-hyperpol-pna/README.md) |
| [第12章 布居分析、NBO 与静电势/RESP 电荷](Gaussian/tutorials/12-population-charges.md) | [12-esp-resp-paracetamol](Gaussian/examples/12-esp-resp-paracetamol/README.md) [12-hirshfeld-cm5-acetate](Gaussian/examples/12-hirshfeld-cm5-acetate/README.md) [12-nbo-pyridine-nma](Gaussian/examples/12-nbo-pyridine-nma/README.md) |
| [第13章 开壳层体系、波函数稳定性与 Broken-Symmetry](Gaussian/tutorials/13-openshell-stability.md) | [13-bs-cu2-dimer](Gaussian/examples/13-bs-cu2-dimer/README.md) [13-diradical-pbenzyne](Gaussian/examples/13-diradical-pbenzyne/README.md) [13-spin-states-o2-ch2](Gaussian/examples/13-spin-states-o2-ch2/README.md) |
| [第14章 弱相互作用与 BSSE（Counterpoise）](Gaussian/tutorials/14-weak-bsse.md) | [14-cp-benzene-dimer](Gaussian/examples/14-cp-benzene-dimer/README.md) [14-cp-crown-k](Gaussian/examples/14-cp-crown-k/README.md) [14-cp-water-dimer](Gaussian/examples/14-cp-water-dimer/README.md) |
| [第15章 ONIOM 多层方法与 QM/MM](Gaussian/tutorials/15-oniom.md) | [15-oniom-hexanal](Gaussian/examples/15-oniom-hexanal/README.md) [15-oniom3-cyclohexanone](Gaussian/examples/15-oniom3-cyclohexanone/README.md) [15-qmmm-phenol-water](Gaussian/examples/15-qmmm-phenol-water/README.md) |
| [第16章 检查点文件、cubegen 与 Multiwfn 后处理](Gaussian/tutorials/16-postprocess-multiwfn.md) | [16-postprocess-multiwfn](Gaussian/examples/16-postprocess-multiwfn/README.md) |
| [第17章 作业脚本、并行与批量计算](Gaussian/tutorials/17-hpc-jobs.md) | [17-hpc-job-scripts](Gaussian/examples/17-hpc-job-scripts/README.md) |
| [第18章 常见报错排查大全与完整研究流程案例](Gaussian/tutorials/18-errors-capstone.md) | [18-capstone-diels-alder](Gaussian/examples/18-capstone-diels-alder/README.md) [18-error-cases](Gaussian/examples/18-error-cases/README.md) |

### CP2K

| 章节 | 示例 |
|---|---|
| [第1章 CP2K 简介、版本与安装](CP2K/tutorials/01-install.md) | [01-install-check](CP2K/examples/01-install-check/README.md) |
| [第2章 输入文件结构与第一个计算](CP2K/tutorials/02-input-structure.md) | [02-energy-force-rutile](CP2K/examples/02-energy-force-rutile/README.md) [02-first-calc-h2o-si](CP2K/examples/02-first-calc-h2o-si/README.md) |
| [第3章 GPW/GAPW 方法、基组与赝势](CP2K/tutorials/03-gpw-basis-pseudo.md) | [03-basis-pseudo-h2o](CP2K/examples/03-basis-pseudo-h2o/README.md) |
| [第4章 CUTOFF / REL_CUTOFF 收敛测试](CP2K/tutorials/04-cutoff-convergence.md) | [04-cutoff-convergence-si](CP2K/examples/04-cutoff-convergence-si/README.md) |
| [第5章 SCF 策略：OT、对角化、展宽与重启](CP2K/tutorials/05-scf-strategy.md) | [05-scf-strategy](CP2K/examples/05-scf-strategy/README.md) |
| [第6章 周期性体系、表面模型与 k 点](CP2K/tutorials/06-periodic-surface-kpoints.md) | [06-periodic-surface-kpoints](CP2K/examples/06-periodic-surface-kpoints/README.md) |
| [第7章 几何优化与晶胞优化](CP2K/tutorials/07-geo-cell-opt.md) | [07-geo-cell-opt](CP2K/examples/07-geo-cell-opt/README.md) |
| [第8章 色散校正、自旋极化与 DFT+U](CP2K/tutorials/08-dispersion-spin-dftu.md) | [08-dispersion-spin-dftu](CP2K/examples/08-dispersion-spin-dftu/README.md) |
| [第9章 杂化泛函与 ADMM](CP2K/tutorials/09-hybrid-admm.md) | [09-hybrid-admm](CP2K/examples/09-hybrid-admm/README.md) |
| [第10章 电子结构分析：DOS、PDOS 与能带](CP2K/tutorials/10-dos-pdos-band.md) | [10-dos-pdos-band](CP2K/examples/10-dos-pdos-band/README.md) |
| [第11章 实空间性质：cube 文件、电荷分析与功函数](CP2K/tutorials/11-cube-charge-workfunction.md) | [11-cube-charge-workfunction](CP2K/examples/11-cube-charge-workfunction/README.md) |
| [第13章 过渡态搜索：NEB 与 Dimer](CP2K/tutorials/13-neb-dimer.md) | [13-neb-dimer-h-cu111](CP2K/examples/13-neb-dimer-h-cu111/README.md) |
| [第14章 振动分析、红外与热化学](CP2K/tutorials/14-vibration-thermo.md) | [14-vibration-thermo](CP2K/examples/14-vibration-thermo/README.md) |

### OpenMX

| 章节 | 示例 |
|---|---|
| [第1章 OpenMX 概述、版本与学习路线](OpenMX/tutorials/01-overview.md) | [01-overview-tour](OpenMX/examples/01-overview-tour/README.md) |
| [第2章 编译与安装（Intel / GCC / ScaLAPACK / ELPA / WSL / deb / Docker）](OpenMX/tutorials/02-install.md) | [02-install-build](OpenMX/examples/02-install-build/README.md) |
| [第3章 输入文件结构与第一个计算（甲烷分子）](OpenMX/tutorials/03-input-first-calc.md) | [03-first-calc-methane](OpenMX/examples/03-first-calc-methane/README.md) |
| [第4章 基组 PAO 与赝势 VPS 的选择](OpenMX/tutorials/04-pao-vps.md) | [04-basis-test-si](OpenMX/examples/04-basis-test-si/README.md) [04-bsse-ghost-h2o-dimer](OpenMX/examples/04-bsse-ghost-h2o-dimer/README.md) |
| [第5章 SCF 参数、k 点与收敛测试](OpenMX/tutorials/05-scf-kpoints-convergence.md) | [05-metal-scf-fe](OpenMX/examples/05-metal-scf-fe/README.md) [05-scf-conv-si](OpenMX/examples/05-scf-conv-si/README.md) |
| [第6章 结构优化、晶胞优化、状态方程与 NEB](OpenMX/tutorials/06-optimization-eos-neb.md) | [06-cellopt-eos-si-c](OpenMX/examples/06-cellopt-eos-si-c/README.md) [06-geomopt-h2o-ch4](OpenMX/examples/06-geomopt-h2o-ch4/README.md) [06-neb-h-in-si](OpenMX/examples/06-neb-h-in-si/README.md) [06-vdw-d3-graphite](OpenMX/examples/06-vdw-d3-graphite/README.md) |
| [第7章 第一性原理分子动力学](OpenMX/tutorials/07-md.md) | [07-md-h2o-nve](OpenMX/examples/07-md-h2o-nve/README.md) [07-md-si64-nvt-npt](OpenMX/examples/07-md-si64-nvt-npt/README.md) |
| [第8章 能带结构（含自动 k 路径与能带展开）](OpenMX/tutorials/08-band.md) | [08-band-cdia-si](OpenMX/examples/08-band-cdia-si/README.md) [08-band-graphene-fe](OpenMX/examples/08-band-graphene-fe/README.md) [08-unfold-si-supercell](OpenMX/examples/08-unfold-si-supercell/README.md) |
| [第9章 态密度 DOS/PDOS、COHP 与光学性质](OpenMX/tutorials/09-dos-pdos-cohp-optics.md) | [09-cohp-cddf-si](OpenMX/examples/09-cohp-cddf-si/README.md) [09-dos-pdos-cdia](OpenMX/examples/09-dos-pdos-cdia/README.md) [09-pdos-mos2-nio](OpenMX/examples/09-pdos-mos2-nio/README.md) |
| [第10章 电荷分析、分子轨道与 cube 可视化](OpenMX/tutorials/10-charge-mo-cube.md) | [10-benzene-mo](OpenMX/examples/10-benzene-mo/README.md) [10-co-cu111-cdd](OpenMX/examples/10-co-cu111-cdd/README.md) [10-h2o-cube-charge](OpenMX/examples/10-h2o-cube-charge/README.md) |
| [第11章 自旋极化、非共线磁性与交换耦合](OpenMX/tutorials/11-spin-noncollinear.md) | [11-cr-bcc-afm](OpenMX/examples/11-cr-bcc-afm/README.md) [11-fe-bcc-fm](OpenMX/examples/11-fe-bcc-fm/README.md) [11-fe-bcc-jx](OpenMX/examples/11-fe-bcc-jx/README.md) [11-mno-noncollinear](OpenMX/examples/11-mno-noncollinear/README.md) |
| [第12章 自旋轨道耦合与拓扑性质](OpenMX/tutorials/12-soc-topology.md) | [12-bi2se3-soc-z2](OpenMX/examples/12-bi2se3-soc-z2/README.md) [12-gaas-soc-band](OpenMX/examples/12-gaas-soc-band/README.md) [12-mos2-soc-band](OpenMX/examples/12-mos2-soc-band/README.md) |
| [第13章 DFT+U（Hubbard U 与 Hund J）](OpenMX/tutorials/13-dft-u.md) | [13-dft-u-nio-mno](OpenMX/examples/13-dft-u-nio-mno/README.md) |
| [第16章 表面与界面：ESM、功函数、电荷掺杂](OpenMX/tutorials/16-surface-esm-workfunction.md) | [16-al111-esm-workfunction](OpenMX/examples/16-al111-esm-workfunction/README.md) [16-graphene-charge-doping](OpenMX/examples/16-graphene-charge-doping/README.md) [16-pt111-esm-charging](OpenMX/examples/16-pt111-esm-charging/README.md) |
| [第17章 NEGF 电子输运](OpenMX/tutorials/17-negf-transport.md) | [17-negf-7agnr](OpenMX/examples/17-negf-7agnr/README.md) [17-negf-carbon-chain](OpenMX/examples/17-negf-carbon-chain/README.md) |

### 通用基础

| 章节 | 示例 |
|---|---|
| [第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输](Common/tutorials/01-linux-hpc.md) | [01-env-setup](Common/examples/01-env-setup/README.md) |

## 🗂️ 仓库结构

```text
DFT-Tutorials/
├── VASP/          # VASP 教程 tutorials/ · 示例 examples/ · 脚本 scripts/
├── Gaussian/      # Gaussian 16 教程、gjf 示例、日志解析与光谱脚本
├── CP2K/          # CP2K 教程、输入示例、脚本；data/ 为 CP2K 基组与赝势数据
├── OpenMX/        # OpenMX 教程、.dat 示例与后处理脚本
├── Common/        # 通用基础：Linux/HPC、作业调度、建模、ASE/pymatgen、可视化、学习路线
├── Scripts/       # NEP/GPUMD 机器学习势、数据集生成、建模与分析脚本
├── Resources/     # GitHub 相关仓库精选与官方文档索引
└── mkdocs.yml     # 文档网站配置（GitHub Pages 自动部署）
```

## 🚀 快速开始

```bash
git clone https://github.com/zenggs05/DFT-Tutorials.git
cd DFT-Tutorials
# 1) 从 Common/ 了解 Linux/HPC 基础，再选择一个软件板块按章节学习
# 2) 复制示例目录到集群运行，例如：
cp -r VASP/examples/<示例目录> ~/work/ && cd ~/work/<示例目录>
# 3) 脚本均支持 -h 查看用法，例如：
python VASP/scripts/plot_band.py -h
```

> [!NOTE]
> 仓库不包含任何商业软件的源码、二进制或 VASP POTCAR 赝势文件（受许可证保护），示例中均给出生成方法。
> 教程中的参考数值仅供对照，请以自己的计算与官方文档为准。

## 🧪 适用版本

| 软件 | 教程主要基于 | 说明 |
|---|---|---|
| VASP | 6.4 – 6.6（最新 6.6.1） | 标注了 KPOINTS_OPT、ML_MODE 等新特性的最低版本 |
| Gaussian | Gaussian 16 Rev. C.01/C.02 + GaussView 6 | |
| CP2K | 2026.2，兼顾 2023.x–2025.x | 新旧语法差异（如 CELL_OPT/TYPE、&DOS/&PDOS）均有对照 |
| OpenMX | 4.0（patch 4.0.1）+ DFT_DATA19，兼顾 3.9 | 标注 4.0 的默认值变化 |

## 🤝 参与贡献

欢迎通过 [Issues](https://github.com/zenggs05/DFT-Tutorials/issues) 报告错误、提出新主题，或提交 Pull Request 补充教程与脚本。如果这个仓库对你有帮助，请点一个 ⭐ Star，让更多做计算的同学看到它！

## 📄 引用与许可

- 如在论文或报告中使用了本仓库的脚本或教程，可通过 GitHub 右侧的 **Cite this repository** 获取引用格式（见 [CITATION.cff](CITATION.cff)）。
- 教程与脚本采用 [MIT License](LICENSE)；`CP2K/data/` 中的基组与赝势数据来自 CP2K 项目（GPL-2.0-or-later），保留其原许可证。

## ⭐ Star History

[![Star History Chart](https://api.star-history.com/svg?repos=zenggs05/DFT-Tutorials&type=Date)](https://star-history.com/#zenggs05/DFT-Tutorials&Date)
