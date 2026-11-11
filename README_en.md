[简体中文](README.md) | [English](README_en.md) | [日本語](README_ja.md) | [Português](README_pt.md)

# DFT-Tutorials: First-Principles Calculation Tutorials & Scripts

**VASP · Gaussian · CP2K · OpenMX — beginner to advanced ｜ ready-to-run inputs ｜ helper scripts ｜ curated GitHub resources**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![Docs](https://img.shields.io/badge/docs-online-blue.svg)](https://zenggs05.github.io/DFT-Tutorials/) [![GitHub stars](https://img.shields.io/github/stars/zenggs05/DFT-Tutorials?style=social)](https://github.com/zenggs05/DFT-Tutorials/stargazers) [![Last commit](https://img.shields.io/github/last-commit/zenggs05/DFT-Tutorials)](https://github.com/zenggs05/DFT-Tutorials/commits/main)

A comprehensive, practical collection of **density functional theory (DFT) tutorials** (written in Simplified Chinese, with English chapter subtitles) for **VASP, Gaussian, CP2K and OpenMX**: installation, input files, convergence tests, geometry/cell optimisation, band structure & DOS, hybrid functionals, DFT+U, spin-orbit coupling, surfaces & adsorption, transition states (NEB/IRC), phonons, ab initio molecular dynamics, TD-DFT excited states, solvation, NBO/charge analysis, electron transport and more. Every chapter ships runnable example inputs, keyword-by-keyword explanations, troubleshooting tips and references, plus Python/Bash helper scripts.

📖 **Read online (searchable website):** <https://zenggs05.github.io/DFT-Tutorials/>

Currently **65 chapters**, **152 example input sets**, **44 scripts**, and a curated index of **230+ related GitHub repositories**.

## Contents

| Section | Scope | Entry |
|---|---|---|
| **VASP** | plane-wave PAW code for solids, surfaces and catalysis | [VASP/README.md](VASP/README.md) |
| **Gaussian** | Gaussian-basis quantum chemistry for molecules, reactions and spectra | [Gaussian/README.md](Gaussian/README.md) |
| **CP2K** | GPW mixed-basis code for large-scale AIMD and QM/MM | [CP2K/README.md](CP2K/README.md) |
| **OpenMX** | numerical atomic orbitals for magnetism, transport and O(N) methods | [OpenMX/README.md](OpenMX/README.md) |
| **Common** | Linux/HPC, modelling, ASE/pymatgen, visualisation, learning path | [Common/README.md](Common/README.md) |
| **NEP/GPUMD scripts** | ML-potential dataset generation, NEP analysis, GPUMD thermal transport plots | [Scripts/README.md](Scripts/README.md) |
| **Resources** | 230+ curated GitHub repositories and official docs | [Resources/github-repos.md](Resources/github-repos.md) |

## Chapter index

### VASP

| # | Topic | 中文标题 |
|---|---|---|
| 1 | [Installing & Compiling VASP 6 on CPU/GPU and Submitting Jobs](VASP/tutorials/01-setup-compile.md) | 第1章 环境准备、VASP 编译安装（CPU/GPU）与作业提交 |
| 2 | [VASP INCAR, POSCAR, KPOINTS, POTCAR & Output Files Explained](VASP/tutorials/02-input-output-files.md) | 第2章 四大输入文件与输出文件全解析 |
| 3 | [Convergence Tests & Accuracy Control in VASP: ENCUT, k-points, SIGMA](VASP/tutorials/03-convergence.md) | 第3章 收敛测试与计算精度控制 |
| 4 | [Structure Relaxation & Equation of State with VASP](VASP/tutorials/04-relax-eos.md) | 第4章 结构优化与状态方程 |
| 5 | [Band Structure & DOS with VASP, incl. KPOINTS_OPT](VASP/tutorials/05-band-dos.md) | 第5章 能带结构与态密度（含 KPOINTS_OPT） |
| 6 | [Hybrid Functionals & meta-GGA in VASP: HSE06, r2SCAN, mBJ](VASP/tutorials/06-hybrid-metagga.md) | 第6章 杂化泛函与 meta-GGA（HSE06 / r2SCAN / mBJ） |
| 7 | [Magnetism, DFT+U & Spin-Orbit Coupling with VASP](VASP/tutorials/07-magnetism-u-soc.md) | 第7章 磁性、DFT+U 与自旋轨道耦合 |
| 8 | [van der Waals Corrections, 2D & Layered Materials with VASP](VASP/tutorials/08-vdw-2d.md) | 第8章 范德华修正、二维与层状材料 |
| 9 | [Surface Slabs, Surface Energy, Work Function & Dipole Correction with VASP](VASP/tutorials/09-surface-workfunction.md) | 第9章 表面模型、表面能、功函数与偶极校正 |
| 10 | [Adsorption, Catalysis, Free Energy and Implicit Solvation with VASP](VASP/tutorials/10-adsorption-catalysis.md) | 第10章 吸附、催化与自由能（含隐式溶剂） |
| 11 | [Transition-State Search in VASP: NEB, CI-NEB and Dimer with VTST](VASP/tutorials/11-neb-dimer.md) | 第11章 过渡态搜索：NEB / CI-NEB / Dimer（VTST） |
| 12 | [Phonons with VASP and phonopy: Finite Displacement and DFPT](VASP/tutorials/12-phonon.md) | 第12章 声子：phonopy 有限位移与 DFPT |
| 13 | [Elastic Constants & Mechanical Properties with VASP](VASP/tutorials/13-elastic.md) | 第13章 弹性常数与力学性质 |
| 14 | [Dielectric & Optical Properties with VASP: LOPTICS / LEPSILON / LCALCEPS](VASP/tutorials/14-optics-dielectric.md) | 第14章 介电与光学性质 |
| 15 | [Ab Initio MD & Machine-Learned Force Fields in VASP](VASP/tutorials/15-md-mlff.md) | 第15章 从头算分子动力学与机器学习力场 |
| 16 | [Charge & Bonding Analysis with VASP](VASP/tutorials/16-charge-bonding.md) | 第16章 电荷与成键分析：Bader、差分电荷、ELF、STM、COHP |
| 17 | [Advanced Electronic Structure with VASP](VASP/tutorials/17-advanced-electronic.md) | 第17章 高级电子结构：Wannier90、有效质量、能带反折叠、费米面与拓扑 |
| 18 | [VASP Post-processing, Automation & Troubleshooting](VASP/tutorials/18-postprocess-troubleshooting.md) | 第18章 后处理生态、自动化与报错排查 |

### Gaussian

| # | Topic | 中文标题 |
|---|---|---|
| 1 | [Getting Started with Gaussian 16 and GaussView 6](Gaussian/tutorials/01-setup-first-job.md) | 第1章 环境准备与第一个 Gaussian 计算 |
| 2 | [Gaussian Input File Anatomy: Link 0, Route, Molecule Specification and Link1](Gaussian/tutorials/02-input-anatomy.md) | 第2章 输入文件结构详解 |
| 3 | [Choosing DFT Functionals, Dispersion Corrections, Basis Sets and ECPs in Gaussian 16](Gaussian/tutorials/03-method-basis.md) | 第3章 方法、泛函与基组选择（含 Gen/GenECP 与赝势） |
| 4 | [Single-Point Energy & SCF Convergence in Gaussian 16](Gaussian/tutorials/04-sp-scf.md) | 第4章 单点能计算与 SCF 收敛技巧 |
| 5 | [Geometry Optimization with Gaussian 16](Gaussian/tutorials/05-opt.md) | 第5章 几何优化 |
| 6 | [Frequencies, IR/Raman & Thermochemistry in Gaussian 16](Gaussian/tutorials/06-freq-thermo.md) | 第6章 频率计算、红外拉曼与热化学 |
| 7 | [Transition State Search & IRC with Gaussian 16](Gaussian/tutorials/07-ts-irc.md) | 第7章 过渡态搜索与 IRC 验证 |
| 8 | [Rigid & Relaxed PES Scans in Gaussian 16](Gaussian/tutorials/08-scan.md) | 第8章 势能面扫描（刚性与柔性扫描） |
| 9 | [Implicit Solvation with PCM & SMD in Gaussian 16](Gaussian/tutorials/09-solvation.md) | 第9章 隐式溶剂模型 PCM/SMD |
| 10 | [Excited States with Gaussian TD-DFT: UV-Vis, Fluorescence & ECD](Gaussian/tutorials/10-tddft.md) | 第10章 激发态 TD-DFT：UV-Vis、荧光发射与 ECD |
| 11 | [NMR Chemical Shifts, J-Couplings, Polarizability & Molecular Properties with Gaussian](Gaussian/tutorials/11-nmr-props.md) | 第11章 NMR 化学位移与其他分子性质 |
| 12 | [Population Analysis, NBO & ESP/RESP Charges with Gaussian](Gaussian/tutorials/12-population-charges.md) | 第12章 布居分析、NBO 与静电势/RESP 电荷 |
| 13 | [Open-Shell Systems, Wavefunction Stability & Broken-Symmetry DFT in Gaussian](Gaussian/tutorials/13-openshell-stability.md) | 第13章 开壳层体系、波函数稳定性与 Broken-Symmetry |
| 14 | [Noncovalent Interactions & Counterpoise BSSE Correction in Gaussian](Gaussian/tutorials/14-weak-bsse.md) | 第14章 弱相互作用与 BSSE（Counterpoise） |
| 15 | [ONIOM Multilayer Methods & QM/MM in Gaussian](Gaussian/tutorials/15-oniom.md) | 第15章 ONIOM 多层方法与 QM/MM |
| 16 | [Checkpoint Files, cubegen & Multiwfn Post-processing](Gaussian/tutorials/16-postprocess-multiwfn.md) | 第16章 检查点文件、cubegen 与 Multiwfn 后处理 |
| 17 | [Running Gaussian 16 on HPC: Job Scripts, Parallelism & Batch Jobs](Gaussian/tutorials/17-hpc-jobs.md) | 第17章 作业脚本、并行与批量计算 |
| 18 | [Gaussian Error Messages & a Complete Research Workflow](Gaussian/tutorials/18-errors-capstone.md) | 第18章 常见报错排查大全与完整研究流程案例 |

### CP2K

| # | Topic | 中文标题 |
|---|---|---|
| 1 | [Introduction, Versions & Installation of CP2K](CP2K/tutorials/01-install.md) | 第1章 CP2K 简介、版本与安装 |
| 2 | [CP2K Input File Structure & First Calculation](CP2K/tutorials/02-input-structure.md) | 第2章 输入文件结构与第一个计算 |
| 3 | [GPW/GAPW Methods, Gaussian Basis Sets & GTH Pseudopotentials in CP2K](CP2K/tutorials/03-gpw-basis-pseudo.md) | 第3章 GPW/GAPW 方法、基组与赝势 |
| 4 | [Converging CUTOFF & REL_CUTOFF in CP2K](CP2K/tutorials/04-cutoff-convergence.md) | 第4章 CUTOFF / REL_CUTOFF 收敛测试 |
| 5 | [SCF Convergence in CP2K: OT, Diagonalization, Smearing & Restart](CP2K/tutorials/05-scf-strategy.md) | 第5章 SCF 策略：OT、对角化、展宽与重启 |
| 6 | [Periodic Systems, Slab Models & k-Points in CP2K](CP2K/tutorials/06-periodic-surface-kpoints.md) | 第6章 周期性体系、表面模型与 k 点 |
| 7 | [Geometry & Cell Optimization with CP2K](CP2K/tutorials/07-geo-cell-opt.md) | 第7章 几何优化与晶胞优化 |
| 8 | [Dispersion, Spin Polarization & DFT+U in CP2K](CP2K/tutorials/08-dispersion-spin-dftu.md) | 第8章 色散校正、自旋极化与 DFT+U |
| 9 | [Hybrid Functionals & ADMM in CP2K](CP2K/tutorials/09-hybrid-admm.md) | 第9章 杂化泛函与 ADMM |
| 10 | [DOS, PDOS & Band Structure with CP2K](CP2K/tutorials/10-dos-pdos-band.md) | 第10章 电子结构分析：DOS、PDOS 与能带 |
| 11 | [Cube Files, Charge Analysis & Work Function in CP2K](CP2K/tutorials/11-cube-charge-workfunction.md) | 第11章 实空间性质：cube 文件、电荷分析与功函数 |
| 12 | [Transition-State Search in CP2K: CI-NEB & Dimer](CP2K/tutorials/13-neb-dimer.md) | 第13章 过渡态搜索：NEB 与 Dimer |
| 13 | [Vibrational Analysis, IR Intensities & Thermochemistry with CP2K](CP2K/tutorials/14-vibration-thermo.md) | 第14章 振动分析、红外与热化学 |

### OpenMX

| # | Topic | 中文标题 |
|---|---|---|
| 1 | [OpenMX Overview, Versions & Learning Path](OpenMX/tutorials/01-overview.md) | 第1章 OpenMX 概述、版本与学习路线 |
| 2 | [Compiling and Installing OpenMX 4.0](OpenMX/tutorials/02-install.md) | 第2章 编译与安装（Intel / GCC / ScaLAPACK / ELPA / WSL / deb / Docker） |
| 3 | [OpenMX Input File and First Calculation: Methane](OpenMX/tutorials/03-input-first-calc.md) | 第3章 输入文件结构与第一个计算（甲烷分子） |
| 4 | [Choosing PAO Basis Sets & VPS Pseudopotentials in OpenMX](OpenMX/tutorials/04-pao-vps.md) | 第4章 基组 PAO 与赝势 VPS 的选择 |
| 5 | [SCF Settings, k-Point Sampling & Convergence Tests in OpenMX](OpenMX/tutorials/05-scf-kpoints-convergence.md) | 第5章 SCF 参数、k 点与收敛测试 |
| 6 | [Geometry & Cell Optimization, Equation of State and NEB with OpenMX](OpenMX/tutorials/06-optimization-eos-neb.md) | 第6章 结构优化、晶胞优化、状态方程与 NEB |
| 7 | [Ab Initio Molecular Dynamics with OpenMX](OpenMX/tutorials/07-md.md) | 第7章 第一性原理分子动力学 |
| 8 | [Band Structure, Automatic k-Path & Band Unfolding with OpenMX](OpenMX/tutorials/08-band.md) | 第8章 能带结构（含自动 k 路径与能带展开） |
| 9 | [DOS, PDOS, COHP/COOP & Optical Properties with OpenMX](OpenMX/tutorials/09-dos-pdos-cohp-optics.md) | 第9章 态密度 DOS/PDOS、COHP 与光学性质 |
| 10 | [Charge Analysis, Molecular Orbitals & Cube Visualization with OpenMX](OpenMX/tutorials/10-charge-mo-cube.md) | 第10章 电荷分析、分子轨道与 cube 可视化 |
| 11 | [Spin Polarization, Noncollinear Magnetism & Exchange Coupling with OpenMX](OpenMX/tutorials/11-spin-noncollinear.md) | 第11章 自旋极化、非共线磁性与交换耦合 |
| 12 | [Spin-Orbit Coupling & Topological Properties with OpenMX](OpenMX/tutorials/12-soc-topology.md) | 第12章 自旋轨道耦合与拓扑性质 |
| 13 | [DFT+U with Hubbard U and Hund's J in OpenMX](OpenMX/tutorials/13-dft-u.md) | 第13章 DFT+U（Hubbard U 与 Hund J） |
| 14 | [Surfaces & Interfaces in OpenMX: ESM, Work Function & Charge Doping](OpenMX/tutorials/16-surface-esm-workfunction.md) | 第16章 表面与界面：ESM、功函数、电荷掺杂 |
| 15 | [Electron Transport with NEGF in OpenMX](OpenMX/tutorials/17-negf-transport.md) | 第17章 NEGF 电子输运 |

### Common

| # | Topic | 中文标题 |
|---|---|---|
| 1 | [Setting up a DFT Computing Environment: Linux/WSL, conda, SSH & File Transfer](Common/tutorials/01-linux-hpc.md) | 第1章 计算环境搭建：Linux/WSL、conda、集群登录与文件传输 |

## Quick start

```bash
git clone https://github.com/zenggs05/DFT-Tutorials.git
cd DFT-Tutorials
python VASP/scripts/plot_band.py -h   # every script supports -h
```

> [!NOTE]
> No proprietary source code, binaries or VASP POTCAR files are included; examples explain how to generate them.

## Contributing

Issues and pull requests are welcome. If you find this useful, please give it a ⭐ star!

## Citation & license

Use **Cite this repository** on GitHub ([CITATION.cff](CITATION.cff)). Tutorials and scripts are [MIT-licensed](LICENSE); `CP2K/data/` comes from the CP2K project under GPL-2.0-or-later.
