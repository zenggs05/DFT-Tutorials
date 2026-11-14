# OpenMX 中文教程：从入门到进阶（OpenMX Tutorial in Chinese）

系统讲解 OpenMX（数值原子轨道 + 模守恒赝势）第一性原理计算：编译安装、.dat 输入文件、PAO/VPS 选择、SCF 与 k 点、结构优化与 NEB、分子动力学、能带与能带展开、DOS/PDOS/COHP、电荷与 cube、自旋极化与非共线磁性、SOC 与拓扑、DFT+U、ESM 与功函数、NEGF 电子输运。

**适用版本**：OpenMX 4.0（patch 4.0.1）+ DFT_DATA19，并标注与 3.9 的差异。

> [!TIP]
> 建议先阅读 [通用基础](../Common/README.md) 了解 Linux 与集群作业提交，再按章节顺序学习；每章末尾都有参考资料与常见报错。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 OpenMX 概述、版本与学习路线（OpenMX Overview, Versions & Learning Path）](tutorials/01-overview.md) | [01-overview-tour](examples/01-overview-tour/README.md) |
| [第2章 编译与安装（Intel / GCC / ScaLAPACK / ELPA / WSL / deb / Docker）（Compiling and Installing OpenMX 4.0）](tutorials/02-install.md) | [02-install-build](examples/02-install-build/README.md) |
| [第3章 输入文件结构与第一个计算（甲烷分子）（OpenMX Input File and First Calculation: Methane）](tutorials/03-input-first-calc.md) | [03-first-calc-methane](examples/03-first-calc-methane/README.md) |
| [第4章 基组 PAO 与赝势 VPS 的选择（Choosing PAO Basis Sets & VPS Pseudopotentials in OpenMX）](tutorials/04-pao-vps.md) | [04-basis-test-si](examples/04-basis-test-si/README.md)<br>[04-bsse-ghost-h2o-dimer](examples/04-bsse-ghost-h2o-dimer/README.md) |
| [第5章 SCF 参数、k 点与收敛测试（SCF Settings, k-Point Sampling & Convergence Tests in OpenMX）](tutorials/05-scf-kpoints-convergence.md) | [05-metal-scf-fe](examples/05-metal-scf-fe/README.md)<br>[05-scf-conv-si](examples/05-scf-conv-si/README.md) |
| [第6章 结构优化、晶胞优化、状态方程与 NEB（Geometry & Cell Optimization, Equation of State and NEB with OpenMX）](tutorials/06-optimization-eos-neb.md) | [06-cellopt-eos-si-c](examples/06-cellopt-eos-si-c/README.md)<br>[06-geomopt-h2o-ch4](examples/06-geomopt-h2o-ch4/README.md)<br>[06-neb-h-in-si](examples/06-neb-h-in-si/README.md)<br>[06-vdw-d3-graphite](examples/06-vdw-d3-graphite/README.md) |
| [第7章 第一性原理分子动力学（Ab Initio Molecular Dynamics with OpenMX）](tutorials/07-md.md) | [07-md-h2o-nve](examples/07-md-h2o-nve/README.md)<br>[07-md-si64-nvt-npt](examples/07-md-si64-nvt-npt/README.md) |
| [第8章 能带结构（含自动 k 路径与能带展开）（Band Structure, Automatic k-Path & Band Unfolding with OpenMX）](tutorials/08-band.md) | [08-band-cdia-si](examples/08-band-cdia-si/README.md)<br>[08-band-graphene-fe](examples/08-band-graphene-fe/README.md)<br>[08-unfold-si-supercell](examples/08-unfold-si-supercell/README.md) |
| [第9章 态密度 DOS/PDOS、COHP 与光学性质（DOS, PDOS, COHP/COOP & Optical Properties with OpenMX）](tutorials/09-dos-pdos-cohp-optics.md) | [09-cohp-cddf-si](examples/09-cohp-cddf-si/README.md)<br>[09-dos-pdos-cdia](examples/09-dos-pdos-cdia/README.md)<br>[09-pdos-mos2-nio](examples/09-pdos-mos2-nio/README.md) |
| [第10章 电荷分析、分子轨道与 cube 可视化（Charge Analysis, Molecular Orbitals & Cube Visualization with OpenMX）](tutorials/10-charge-mo-cube.md) | [10-benzene-mo](examples/10-benzene-mo/README.md)<br>[10-co-cu111-cdd](examples/10-co-cu111-cdd/README.md)<br>[10-h2o-cube-charge](examples/10-h2o-cube-charge/README.md) |
| [第11章 自旋极化、非共线磁性与交换耦合（Spin Polarization, Noncollinear Magnetism & Exchange Coupling with OpenMX）](tutorials/11-spin-noncollinear.md) | [11-cr-bcc-afm](examples/11-cr-bcc-afm/README.md)<br>[11-fe-bcc-fm](examples/11-fe-bcc-fm/README.md)<br>[11-fe-bcc-jx](examples/11-fe-bcc-jx/README.md)<br>[11-mno-noncollinear](examples/11-mno-noncollinear/README.md) |
| [第12章 自旋轨道耦合与拓扑性质（Spin-Orbit Coupling & Topological Properties with OpenMX）](tutorials/12-soc-topology.md) | [12-bi2se3-soc-z2](examples/12-bi2se3-soc-z2/README.md)<br>[12-gaas-soc-band](examples/12-gaas-soc-band/README.md)<br>[12-mos2-soc-band](examples/12-mos2-soc-band/README.md) |
| [第13章 DFT+U（Hubbard U 与 Hund J）（DFT+U with Hubbard U and Hund's J in OpenMX）](tutorials/13-dft-u.md) | [13-dft-u-nio-mno](examples/13-dft-u-nio-mno/README.md) |
| [第14章 Wannier 函数：MLWF 与 CWF（Wannier Functions in OpenMX: MLWF and CWF）](tutorials/14-wannier.md) | [14-gaas-cwf](examples/14-gaas-cwf/README.md)<br>[14-graphene-pz-mlwf](examples/14-graphene-pz-mlwf/README.md)<br>[14-si-mlwf-sp3](examples/14-si-mlwf-sp3/README.md) |
| [第15章 Berry 相宏观极化与外加电场（Berry-Phase Polarization & External Electric Fields in OpenMX）](tutorials/15-berry-efield.md) | [15-batio3-polarization](examples/15-batio3-polarization/README.md)<br>[15-nacl-polb-born](examples/15-nacl-polb-born/README.md)<br>[15-nitrobenzene-efield](examples/15-nitrobenzene-efield/README.md) |
| [第16章 表面与界面：ESM、功函数、电荷掺杂（Surfaces & Interfaces in OpenMX: ESM, Work Function & Charge Doping）](tutorials/16-surface-esm-workfunction.md) | [16-al111-esm-workfunction](examples/16-al111-esm-workfunction/README.md)<br>[16-graphene-charge-doping](examples/16-graphene-charge-doping/README.md)<br>[16-pt111-esm-charging](examples/16-pt111-esm-charging/README.md) |
| [第17章 NEGF 电子输运（Electron Transport with NEGF in OpenMX）](tutorials/17-negf-transport.md) | [17-negf-7agnr](examples/17-negf-7agnr/README.md)<br>[17-negf-carbon-chain](examples/17-negf-carbon-chain/README.md) |
| [第18章 并行、大体系 O(N)、作业脚本与常见错误排查（Parallelization, Large-Scale O(N) Methods, Job Scripts & Troubleshooting in OpenMX）](tutorials/18-parallel-on-errors.md) | [18-api-python](examples/18-api-python/README.md)<br>[18-c60-cluster-elpa](examples/18-c60-cluster-elpa/README.md)<br>[18-job-scripts](examples/18-job-scripts/README.md)<br>[18-on-si-dc-krylov](examples/18-on-si-dc-krylov/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-overview-tour](examples/01-overview-tour/README.md) | 示例 01：OpenMX 安装目录与数据库导览 |
| [02-install-build](examples/02-install-build/README.md) | 示例 02：编译安装 OpenMX 4.0 + patch 4.0.1 并验证 |
| [03-first-calc-methane](examples/03-first-calc-methane/README.md) | 示例 03：第一个 OpenMX 计算——CH4、H2O、CO 分子 |
| [04-basis-test-si](examples/04-basis-test-si/README.md) | 示例 04-1：Si 的基组收敛测试 |
| [04-bsse-ghost-h2o-dimer](examples/04-bsse-ghost-h2o-dimer/README.md) | 示例 04-2：用空原子 E 做水二聚体的 counterpoise |
| [05-metal-scf-fe](examples/05-metal-scf-fe/README.md) | 示例 05-2：bcc Fe 的金属 SCF 设置、电子温度与 k 点扫描 |
| [05-scf-conv-si](examples/05-scf-conv-si/README.md) | 示例 05-1：Si 的 SCF 参数、截断能与 k 点收敛测试，以及 `scf.restart pp` |
| [06-cellopt-eos-si-c](examples/06-cellopt-eos-si-c/README.md) | 示例 06-2：Si 的 RFC5 晶胞优化与金刚石的 EvsLC 状态方程 |
| [06-geomopt-h2o-ch4](examples/06-geomopt-h2o-ch4/README.md) | 示例 06-1：分子结构优化——H2O |
| [06-neb-h-in-si](examples/06-neb-h-in-si/README.md) | 示例 06-4：NEB——Si 中间隙 H 在相邻四面体位之间的扩散 |
| [06-vdw-d3-graphite](examples/06-vdw-d3-graphite/README.md) | 示例 06-3：DFT-D3 范德华修正——石墨层间晶格常数与双层石墨烯层间距 |
| [07-md-h2o-nve](examples/07-md-h2o-nve/README.md) | 算例 7-1：H₂O 分子的 NVE 分子动力学 |
| [07-md-si64-nvt-npt](examples/07-md-si64-nvt-npt/README.md) | 算例 7-2：Si 的 NVT |
| [08-band-cdia-si](examples/08-band-cdia-si/README.md) | 算例 8-1：金刚石与硅的能带结构 |
| [08-band-graphene-fe](examples/08-band-graphene-fe/README.md) | 算例 8-2：石墨烯 |
| [08-unfold-si-supercell](examples/08-unfold-si-supercell/README.md) | 算例 8-3：Si 超胞的能带展开 |
| [09-cohp-cddf-si](examples/09-cohp-cddf-si/README.md) | 算例 9-3：硅的 COHP/COOP 成键分析与光学介电函数 |
| [09-dos-pdos-cdia](examples/09-dos-pdos-cdia/README.md) | 算例 9-1：金刚石的总态密度与投影态密度 |
| [09-pdos-mos2-nio](examples/09-pdos-mos2-nio/README.md) | 算例 9-2：单层 MoS₂ 与反铁磁 NiO 的投影态密度 |
| [10-benzene-mo](examples/10-benzene-mo/README.md) | 算例 10-benzene-mo：苯分子的 HOMO / LUMO cube 文件 |
| [10-co-cu111-cdd](examples/10-co-cu111-cdd/README.md) | 算例 10-co-cu111-cdd：CO/Cu |
| [10-h2o-cube-charge](examples/10-h2o-cube-charge/README.md) | 算例 10-h2o-cube-charge：水分子的 cube 输出与四种原子电荷 |
| [11-cr-bcc-afm](examples/11-cr-bcc-afm/README.md) | 算例 11-cr-bcc-afm：bcc Cr 的反铁磁 / 铁磁 / 非磁对比 |
| [11-fe-bcc-fm](examples/11-fe-bcc-fm/README.md) | 算例 11-fe-bcc-fm：bcc Fe 铁磁共线自旋极化计算 |
| [11-fe-bcc-jx](examples/11-fe-bcc-jx/README.md) | 算例 11-fe-bcc-jx：用 jx 计算 bcc Fe 的交换耦合 J_ij 与平均场居里温度 |
| [11-mno-noncollinear](examples/11-mno-noncollinear/README.md) | 算例 11-mno-noncollinear：MnO 分子的非共线 + SOC + 约束 DFT |
| [12-bi2se3-soc-z2](examples/12-bi2se3-soc-z2/README.md) | 算例 12-bi2se3-soc-z2：Bi2Se3 的能带反转与 Z2 拓扑不变量 |
| [12-gaas-soc-band](examples/12-gaas-soc-band/README.md) | 算例 12-gaas-soc-band：GaAs 能带的自旋轨道劈裂 |
| [12-mos2-soc-band](examples/12-mos2-soc-band/README.md) | 算例 12-mos2-soc-band：单层 MoS2 的 SOC 价带劈裂与自旋–能谷锁定 |
| [13-dft-u-nio-mno](examples/13-dft-u-nio-mno/README.md) | 示例 13：NiO / MnO 的 DFT+U |
| [14-gaas-cwf](examples/14-gaas-cwf/README.md) | 示例 14-3：GaAs 的最接近 Wannier 函数 |
| [14-graphene-pz-mlwf](examples/14-graphene-pz-mlwf/README.md) | 示例 14-2：石墨烯 pz 双带 MLWF |
| [14-si-mlwf-sp3](examples/14-si-mlwf-sp3/README.md) | 示例 14-1：Si 的最大局域化 Wannier 函数 |
| [15-batio3-polarization](examples/15-batio3-polarization/README.md) | 示例 15-2：四方 BaTiO3 的自发极化 |
| [15-nacl-polb-born](examples/15-nacl-polb-born/README.md) | 示例 15-1：NaCl 中 Na 的 Born 有效电荷 |
| [15-nitrobenzene-efield](examples/15-nitrobenzene-efield/README.md) | 示例 15-3：硝基苯在外加均匀电场中的响应 |
| [16-al111-esm-workfunction](examples/16-al111-esm-workfunction/README.md) | 示例 16-1：Al |
| [16-graphene-charge-doping](examples/16-graphene-charge-doping/README.md) | 示例 16-3：石墨烯电荷掺杂 |
| [16-pt111-esm-charging](examples/16-pt111-esm-charging/README.md) | 示例 16-2：带电 Pt |
| [17-negf-7agnr](examples/17-negf-7agnr/README.md) | 示例 17-2：7-扶手椅型石墨烯纳米带 |
| [17-negf-carbon-chain](examples/17-negf-carbon-chain/README.md) | 示例 17-1：碳原子链的 NEGF 输运 |
| [18-api-python](examples/18-api-python/README.md) | 示例 18-4：用 OpenMX 4.0 Python API 在一个 MPI 作业中并行跑多个输入 |
| [18-c60-cluster-elpa](examples/18-c60-cluster-elpa/README.md) | 示例 18-3：C60 的 Cluster 对角化 |
| [18-job-scripts](examples/18-job-scripts/README.md) | 示例 18-1：OpenMX 作业脚本 |
| [18-on-si-dc-krylov](examples/18-on-si-dc-krylov/README.md) | 示例 18-2：Si 超胞的 O |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。

- [job_omx_pbs.sh](scripts/job_omx_pbs.sh)
- [job_omx_slurm.sh](scripts/job_omx_slurm.sh)
- [omx_band_plot.py](scripts/omx_band_plot.py)
- [omx_build.sh](scripts/omx_build.sh)
- [omx_conv_test.py](scripts/omx_conv_test.py)
- [omx_cube_tools.py](scripts/omx_cube_tools.py)
- [omx_dos_auto.py](scripts/omx_dos_auto.py)
- [omx_eos_fit.py](scripts/omx_eos_fit.py)
- [omx_magconf.py](scripts/omx_magconf.py)
- [omx_make_input.py](scripts/omx_make_input.py)
- [omx_negf_plot.py](scripts/omx_negf_plot.py)
- [omx_parse_out.py](scripts/omx_parse_out.py)
- [omx_traj.py](scripts/omx_traj.py)

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
