# CP2K 中文教程：从入门到进阶（CP2K Tutorial in Chinese）

系统讲解 CP2K（GPW/GAPW 混合基组）第一性原理计算：安装、输入文件结构、基组与 GTH 赝势、CUTOFF 收敛、OT 与对角化 SCF、周期性体系与 k 点、几何/晶胞优化、DFT-D3、自旋与 DFT+U、杂化泛函与 ADMM、DOS/PDOS/能带、cube 文件与功函数、NEB、振动分析等。

**适用版本**：CP2K 2026.2，并对照 2023.x–2025.x 的旧语法（如 CELL_OPT/TYPE、&PDOS → &DOS）。

> [!TIP]
> 建议先阅读 [通用基础](../Common/README.md) 了解 Linux 与集群作业提交，再按章节顺序学习；每章末尾都有参考资料与常见报错。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 CP2K 简介、版本与安装（Introduction, Versions & Installation of CP2K）](tutorials/01-install.md) | [01-install-check](examples/01-install-check/README.md) |
| [第2章 输入文件结构与第一个计算（CP2K Input File Structure & First Calculation）](tutorials/02-input-structure.md) | [02-energy-force-rutile](examples/02-energy-force-rutile/README.md)<br>[02-first-calc-h2o-si](examples/02-first-calc-h2o-si/README.md) |
| [第3章 GPW/GAPW 方法、基组与赝势（GPW/GAPW Methods, Gaussian Basis Sets & GTH Pseudopotentials in CP2K）](tutorials/03-gpw-basis-pseudo.md) | [03-basis-pseudo-h2o](examples/03-basis-pseudo-h2o/README.md) |
| [第4章 CUTOFF / REL_CUTOFF 收敛测试（Converging CUTOFF & REL_CUTOFF in CP2K）](tutorials/04-cutoff-convergence.md) | [04-cutoff-convergence-si](examples/04-cutoff-convergence-si/README.md) |
| [第5章 SCF 策略：OT、对角化、展宽与重启（SCF Convergence in CP2K: OT, Diagonalization, Smearing & Restart）](tutorials/05-scf-strategy.md) | [05-scf-strategy](examples/05-scf-strategy/README.md) |
| [第6章 周期性体系、表面模型与 k 点（Periodic Systems, Slab Models & k-Points in CP2K）](tutorials/06-periodic-surface-kpoints.md) | [06-periodic-surface-kpoints](examples/06-periodic-surface-kpoints/README.md) |
| [第7章 几何优化与晶胞优化（Geometry & Cell Optimization with CP2K）](tutorials/07-geo-cell-opt.md) | [07-geo-cell-opt](examples/07-geo-cell-opt/README.md) |
| [第8章 色散校正、自旋极化与 DFT+U（Dispersion, Spin Polarization & DFT+U in CP2K）](tutorials/08-dispersion-spin-dftu.md) | [08-dispersion-spin-dftu](examples/08-dispersion-spin-dftu/README.md) |
| [第9章 杂化泛函与 ADMM（Hybrid Functionals & ADMM in CP2K）](tutorials/09-hybrid-admm.md) | [09-hybrid-admm](examples/09-hybrid-admm/README.md) |
| [第10章 电子结构分析：DOS、PDOS 与能带（DOS, PDOS & Band Structure with CP2K）](tutorials/10-dos-pdos-band.md) | [10-dos-pdos-band](examples/10-dos-pdos-band/README.md) |
| [第11章 实空间性质：cube 文件、电荷分析与功函数（Cube Files, Charge Analysis & Work Function in CP2K）](tutorials/11-cube-charge-workfunction.md) | [11-cube-charge-workfunction](examples/11-cube-charge-workfunction/README.md) |
| [第12章 从头算分子动力学（AIMD）](tutorials/12-aimd.md) | [12-aimd-water](examples/12-aimd-water/README.md) |
| [第13章 过渡态搜索：NEB 与 Dimer（Transition-State Search in CP2K: CI-NEB & Dimer）](tutorials/13-neb-dimer.md) | [13-neb-dimer-h-cu111](examples/13-neb-dimer-h-cu111/README.md) |
| [第14章 振动分析、红外与热化学（Vibrational Analysis, IR Intensities & Thermochemistry with CP2K）](tutorials/14-vibration-thermo.md) | [14-vibration-thermo](examples/14-vibration-thermo/README.md) |
| [第15章 增强采样：元动力学、PLUMED 与约束 MD（Enhanced Sampling with CP2K: Metadynamics, PLUMED & Constrained MD）](tutorials/15-metadynamics-plumed.md) | [15-metadynamics-nacl](examples/15-metadynamics-nacl/README.md) |
| [第16章 QM/MM 与半经验方法（xTB/DFTB）（QM/MM and Semi-empirical Methods: GFN-xTB, DFTB）](tutorials/16-qmmm-xtb.md) | [16-qmmm-xtb](examples/16-qmmm-xtb/README.md) |
| [第17章 自动化与工具生态：ASE、cp2k-input-tools、Multiwfn、cp2kdata（Automation & Tool Ecosystem: ASE, cp2k-input-tools, Multiwfn, cp2kdata）](tutorials/17-tools-ecosystem.md) | [17-automation-tools](examples/17-automation-tools/README.md) |
| [第18章 常见错误排查、性能调优与版本迁移（Troubleshooting, Performance Tuning & Version Migration）](tutorials/18-troubleshooting.md) | [18-troubleshooting](examples/18-troubleshooting/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-install-check](examples/01-install-check/README.md) | 示例 01：CP2K 安装自检 |
| [02-energy-force-rutile](examples/02-energy-force-rutile/README.md) | 示例 02b：金红石 TiO2 的能量与原子力 |
| [02-first-calc-h2o-si](examples/02-first-calc-h2o-si/README.md) | 示例 02a：第一个 CP2K 计算——H2O 分子与 Si 晶体 |
| [03-basis-pseudo-h2o](examples/03-basis-pseudo-h2o/README.md) | 示例 03：基组、赝势与 GAPW 全电子计算 |
| [04-cutoff-convergence-si](examples/04-cutoff-convergence-si/README.md) | 示例 04：体相 Si 与 H2O 的 CUTOFF / REL_CUTOFF 收敛测试 |
| [05-scf-strategy](examples/05-scf-strategy/README.md) | 示例 05：SCF 策略——OT、对角化 + 混合、金属展宽与波函数重启 |
| [06-periodic-surface-kpoints](examples/06-periodic-surface-kpoints/README.md) | 示例 06：孤立分子、表面 slab 偶极校正与 k 点收敛 |
| [07-geo-cell-opt](examples/07-geo-cell-opt/README.md) | 示例 07：CP2K 几何优化与晶胞优化 |
| [08-dispersion-spin-dftu](examples/08-dispersion-spin-dftu/README.md) | 示例 08：CP2K 色散校正、自旋极化与 DFT+U |
| [09-hybrid-admm](examples/09-hybrid-admm/README.md) | 示例 09：CP2K 杂化泛函与 ADMM |
| [10-dos-pdos-band](examples/10-dos-pdos-band/README.md) | 示例 10：DOS、PDOS 与能带 |
| [11-cube-charge-workfunction](examples/11-cube-charge-workfunction/README.md) | 示例 11：cube 文件、电荷分析、molden 与功函数 |
| [12-aimd-water](examples/12-aimd-water/README.md) | 示例 12：32 个水分子的 AIMD |
| [13-neb-dimer-h-cu111](examples/13-neb-dimer-h-cu111/README.md) | 示例 13：H 在 Cu |
| [14-vibration-thermo](examples/14-vibration-thermo/README.md) | 示例 14：振动分析、红外强度与热化学 |
| [15-metadynamics-nacl](examples/15-metadynamics-nacl/README.md) | 示例 15：NaCl 在水中解离的自由能 |
| [16-qmmm-xtb](examples/16-qmmm-xtb/README.md) | 示例 16：半经验方法 |
| [17-automation-tools](examples/17-automation-tools/README.md) | 示例 17：用 ASE、cp2k-input-tools、Multiwfn、cp2kdata、AiiDA 自动化 CP2K |
| [18-troubleshooting](examples/18-troubleshooting/README.md) | 示例 18：报错排查练习、版本迁移、重启与并行性能测试 |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。

- [band_plot.py](scripts/band_plot.py)
- [conv_cutoff.py](scripts/conv_cutoff.py)
- [conv_kpoints.py](scripts/conv_kpoints.py)
- [cp2k_inpgen.py](scripts/cp2k_inpgen.py)
- [cp2k_outparse.py](scripts/cp2k_outparse.py)
- [cube_tools.py](scripts/cube_tools.py)
- [job_cp2k_pbs.sh](scripts/job_cp2k_pbs.sh)
- [job_cp2k_slurm.sh](scripts/job_cp2k_slurm.sh)
- [md_analyze.py](scripts/md_analyze.py)
- [neb_plot.py](scripts/neb_plot.py)
- [pdos_plot.py](scripts/pdos_plot.py)
- [plot_opt.py](scripts/plot_opt.py)
- [restart_update.py](scripts/restart_update.py)
- [vib_thermo.py](scripts/vib_thermo.py)

## 基组与赝势数据

[data/](data/README.md) 目录收录了 CP2K 项目发布的基组与赝势数据文件（GPL-2.0-or-later），便于离线查阅与使用。

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
