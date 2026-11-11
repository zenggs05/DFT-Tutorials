# VASP 中文教程：从入门到进阶（VASP Tutorial in Chinese）

系统讲解 VASP（Vienna Ab initio Simulation Package）第一性原理计算：编译安装、INCAR/POSCAR/KPOINTS/POTCAR、收敛测试、结构优化、能带与态密度、HSE06/r2SCAN、磁性与 DFT+U、SOC、范德华修正、表面与功函数、吸附与催化、NEB/Dimer 过渡态、phonopy 声子、弹性常数、光学性质、AIMD 与机器学习力场、Bader/差分电荷/COHP、Wannier90 与报错排查。

**适用版本**：VASP 6.4 – 6.6（最新 6.6.1），新特性均标注最低版本。

> [!TIP]
> 建议先阅读 [通用基础](../Common/README.md) 了解 Linux 与集群作业提交，再按章节顺序学习；每章末尾都有参考资料与常见报错。

## 教程目录

| 章节 | 配套示例 |
|---|---|
| [第1章 环境准备、VASP 编译安装（CPU/GPU）与作业提交（Installing & Compiling VASP 6 on CPU/GPU and Submitting Jobs）](tutorials/01-setup-compile.md) | [01-setup-compile](examples/01-setup-compile/README.md) |
| [第2章 四大输入文件与输出文件全解析（VASP INCAR, POSCAR, KPOINTS, POTCAR & Output Files Explained）](tutorials/02-input-output-files.md) | [02-input-files-si](examples/02-input-files-si/README.md) |
| [第3章 收敛测试与计算精度控制（Convergence Tests & Accuracy Control in VASP: ENCUT, k-points, SIGMA）](tutorials/03-convergence.md) | [03-convergence-si-al](examples/03-convergence-si-al/README.md) |
| [第4章 结构优化与状态方程（Structure Relaxation & Equation of State with VASP）](tutorials/04-relax-eos.md) | [04-eos-si](examples/04-eos-si/README.md)<br>[04-relax-co-molecule](examples/04-relax-co-molecule/README.md)<br>[04-relax-mos2-2d](examples/04-relax-mos2-2d/README.md)<br>[04-relax-si-isif3](examples/04-relax-si-isif3/README.md) |
| [第5章 能带结构与态密度（含 KPOINTS_OPT）（Band Structure & DOS with VASP, incl. KPOINTS_OPT）](tutorials/05-band-dos.md) | [05-band-dos-si](examples/05-band-dos-si/README.md) |
| [第6章 杂化泛函与 meta-GGA（HSE06 / r2SCAN / mBJ）（Hybrid Functionals & meta-GGA in VASP: HSE06, r2SCAN, mBJ）](tutorials/06-hybrid-metagga.md) | [06-hse06-band-si](examples/06-hse06-band-si/README.md)<br>[06-mbj-gap-si](examples/06-mbj-gap-si/README.md)<br>[06-r2scan-relax-zno](examples/06-r2scan-relax-zno/README.md) |
| [第7章 磁性、DFT+U 与自旋轨道耦合（Magnetism, DFT+U & Spin-Orbit Coupling with VASP）](tutorials/07-magnetism-u-soc.md) | [07-fe-bcc-fm](examples/07-fe-bcc-fm/README.md)<br>[07-fe-bcc-soc-mae](examples/07-fe-bcc-soc-mae/README.md)<br>[07-linear-response-u](examples/07-linear-response-u/README.md)<br>[07-nio-afm-dftu](examples/07-nio-afm-dftu/README.md) |
| [第8章 范德华修正、二维与层状材料（van der Waals Corrections, 2D & Layered Materials with VASP）](tutorials/08-vdw-2d.md) | [08-graphite-vdw](examples/08-graphite-vdw/README.md)<br>[08-mos2-bilayer-binding](examples/08-mos2-bilayer-binding/README.md)<br>[08-mos2-monolayer](examples/08-mos2-monolayer/README.md) |
| [第9章 表面模型、表面能、功函数与偶极校正（Surface Slabs, Surface Energy, Work Function & Dipole Correction with VASP）](tutorials/09-surface-workfunction.md) | [09-cu100-surface-energy](examples/09-cu100-surface-energy/README.md)<br>[09-pt111-slab](examples/09-pt111-slab/README.md)<br>[09-pt111-workfunction-dipole](examples/09-pt111-workfunction-dipole/README.md) |
| [第10章 吸附、催化与自由能（含隐式溶剂）（Adsorption, Catalysis, Free Energy and Implicit Solvation with VASP）](tutorials/10-adsorption-catalysis.md) | [10-co-pt111-adsorption](examples/10-co-pt111-adsorption/README.md)<br>[10-her-oer-free-energy](examples/10-her-oer-free-energy/README.md)<br>[10-vaspsol-implicit-solvent](examples/10-vaspsol-implicit-solvent/README.md) |
| [第11章 过渡态搜索：NEB / CI-NEB / Dimer（VTST）（Transition-State Search in VASP: NEB, CI-NEB and Dimer with VTST）](tutorials/11-neb-dimer.md) | [11-dimer-from-neb](examples/11-dimer-from-neb/README.md)<br>[11-li-si-neb](examples/11-li-si-neb/README.md)<br>[11-pt-adatom-neb](examples/11-pt-adatom-neb/README.md)<br>[11-ts-freq-check](examples/11-ts-freq-check/README.md) |
| [第12章 声子：phonopy 有限位移与 DFPT（Phonons with VASP and phonopy: Finite Displacement and DFPT）](tutorials/12-phonon.md) | [12-graphene-phonon](examples/12-graphene-phonon/README.md)<br>[12-mgo-dfpt-nac](examples/12-mgo-dfpt-nac/README.md)<br>[12-si-phonopy-fd](examples/12-si-phonopy-fd/README.md) |
| [第13章 弹性常数与力学性质（Elastic Constants & Mechanical Properties with VASP）](tutorials/13-elastic.md) | [13-elastic-al-ibrion6](examples/13-elastic-al-ibrion6/README.md)<br>[13-elastic-mos2-2d](examples/13-elastic-mos2-2d/README.md) |
| [第14章 介电与光学性质（Dielectric & Optical Properties with VASP: LOPTICS / LEPSILON / LCALCEPS）](tutorials/14-optics-dielectric.md) | [14-dielectric-gaas-lepsilon](examples/14-dielectric-gaas-lepsilon/README.md)<br>[14-dielectric-zno-hse-lcalceps](examples/14-dielectric-zno-hse-lcalceps/README.md)<br>[14-optics-si-loptics](examples/14-optics-si-loptics/README.md) |
| [第15章 从头算分子动力学与机器学习力场（Ab Initio MD & Machine-Learned Force Fields in VASP）](tutorials/15-md-mlff.md) | [15-aimd-si-liquid-nvt](examples/15-aimd-si-liquid-nvt/README.md)<br>[15-aimd-water-npt](examples/15-aimd-water-npt/README.md)<br>[15-mlff-si-train-run](examples/15-mlff-si-train-run/README.md) |
| [第16章 电荷与成键分析：Bader、差分电荷、ELF、STM、COHP（Charge & Bonding Analysis with VASP）](tutorials/16-charge-bonding.md) | [16-bader-co-pt111](examples/16-bader-co-pt111/README.md)<br>[16-cdd-co-pt111](examples/16-cdd-co-pt111/README.md)<br>[16-elf-si](examples/16-elf-si/README.md)<br>[16-lobster-gaas](examples/16-lobster-gaas/README.md)<br>[16-stm-mos2](examples/16-stm-mos2/README.md) |
| [第17章 高级电子结构：Wannier90、有效质量、能带反折叠、费米面与拓扑（Advanced Electronic Structure with VASP）](tutorials/17-advanced-electronic.md) | [17-effmass-gaas](examples/17-effmass-gaas/README.md)<br>[17-fermi-surface-cu](examples/17-fermi-surface-cu/README.md)<br>[17-unfold-si-vacancy](examples/17-unfold-si-vacancy/README.md)<br>[17-wannier90-si](examples/17-wannier90-si/README.md) |
| [第18章 后处理生态、自动化与报错排查（VASP Post-processing, Automation & Troubleshooting）](tutorials/18-postprocess-troubleshooting.md) | [18-postprocess-toolkit](examples/18-postprocess-toolkit/README.md) |

## 全部示例

| 示例目录 | 说明 |
|---|---|
| [01-setup-compile](examples/01-setup-compile/README.md) | 示例 01：VASP 编译配置、作业脚本与安装验证 |
| [02-input-files-si](examples/02-input-files-si/README.md) | 示例 02：四大输入文件模板与 Si 第一个单点能计算 |
| [03-convergence-si-al](examples/03-convergence-si-al/README.md) | 示例 03：ENCUT、k 点与展宽 |
| [04-eos-si](examples/04-eos-si/README.md) | 示例 04-4：Si 的 E–V 扫描与 Birch–Murnaghan 状态方程拟合 |
| [04-relax-co-molecule](examples/04-relax-co-molecule/README.md) | 示例 04-3：CO 分子在真空盒子中的结构优化 |
| [04-relax-mos2-2d](examples/04-relax-mos2-2d/README.md) | 示例 04-2：单层 MoS2 只优化面内晶格 |
| [04-relax-si-isif3](examples/04-relax-si-isif3/README.md) | 示例 04-1：Si 体相变胞优化 |
| [05-band-dos-si](examples/05-band-dos-si/README.md) | 示例 05：Si 的能带结构与态密度 |
| [06-hse06-band-si](examples/06-hse06-band-si/README.md) | 示例 06-1：Si 的 HSE06 自洽与能带 |
| [06-mbj-gap-si](examples/06-mbj-gap-si/README.md) | 示例 06-3：Si 的 mBJ 带隙与能带 |
| [06-r2scan-relax-zno](examples/06-r2scan-relax-zno/README.md) | 示例 06-2：纤锌矿 ZnO 的 r2SCAN 变胞优化 |
| [07-fe-bcc-fm](examples/07-fe-bcc-fm/README.md) | 算例 07-fe-bcc-fm：bcc Fe 铁磁自旋极化计算 |
| [07-fe-bcc-soc-mae](examples/07-fe-bcc-soc-mae/README.md) | 算例 07-fe-bcc-soc-mae：bcc Fe 的自旋轨道耦合与磁晶各向异性能 |
| [07-linear-response-u](examples/07-linear-response-u/README.md) | 算例 07-linear-response-u：线性响应法 |
| [07-nio-afm-dftu](examples/07-nio-afm-dftu/README.md) | 算例 07-nio-afm-dftu：NiO 反铁磁 PBE+U 与交换耦合常数 J1/J2 |
| [08-graphite-vdw](examples/08-graphite-vdw/README.md) | 算例 08-graphite-vdw：石墨层间距与层间结合能——PBE / PBE-D3 |
| [08-mos2-bilayer-binding](examples/08-mos2-bilayer-binding/README.md) | 算例 08-mos2-bilayer-binding：双层 MoS2 层间结合能曲线 |
| [08-mos2-monolayer](examples/08-mos2-monolayer/README.md) | 算例 08-mos2-monolayer：单层 MoS2——二维材料的 vdW 一致设置、真空收敛与真空能级 |
| [09-cu100-surface-energy](examples/09-cu100-surface-energy/README.md) | 算例 09-cu100-surface-energy：Cu |
| [09-pt111-slab](examples/09-pt111-slab/README.md) | 算例 09-pt111-slab：Pt |
| [09-pt111-workfunction-dipole](examples/09-pt111-workfunction-dipole/README.md) | 算例 09-pt111-workfunction-dipole：Pt |
| [10-co-pt111-adsorption](examples/10-co-pt111-adsorption/README.md) | 示例 10-1：CO 在 Pt |
| [10-her-oer-free-energy](examples/10-her-oer-free-energy/README.md) | 示例 10-2：计算氢电极——HER 的 ΔG_H\* 与 OER/ORR 台阶图 |
| [10-vaspsol-implicit-solvent](examples/10-vaspsol-implicit-solvent/README.md) | 示例 10-3：VASPsol 隐式水溶剂与线性化 Poisson–Boltzmann 电解质 |
| [11-dimer-from-neb](examples/11-dimer-from-neb/README.md) | 示例 11-3：从 NEB 结果出发用 Dimer 精修鞍点 |
| [11-li-si-neb](examples/11-li-si-neb/README.md) | 示例 11-2：Li 间隙原子在晶体 Si 中的扩散 |
| [11-pt-adatom-neb](examples/11-pt-adatom-neb/README.md) | 示例 11-1：Pt 吸附原子在 Pt |
| [11-ts-freq-check](examples/11-ts-freq-check/README.md) | 示例 11-4：过渡态验证——爬坡图像的有限差分频率分析 |
| [12-graphene-phonon](examples/12-graphene-phonon/README.md) | 示例 12-3：石墨烯声子谱 |
| [12-mgo-dfpt-nac](examples/12-mgo-dfpt-nac/README.md) | 示例 12-2：MgO 的 DFPT 力常数与非解析项修正 |
| [12-si-phonopy-fd](examples/12-si-phonopy-fd/README.md) | 示例 12-1：Si 的声子谱、声子 DOS 与热力学性质 |
| [13-elastic-al-ibrion6](examples/13-elastic-al-ibrion6/README.md) | 算例 13-1：fcc Al 弹性常数 |
| [13-elastic-mos2-2d](examples/13-elastic-mos2-2d/README.md) | 算例 13-2：单层 MoS2 面内弹性常数 |
| [14-dielectric-gaas-lepsilon](examples/14-dielectric-gaas-lepsilon/README.md) | 算例 14-2：GaAs 的静态介电张量与 Born 有效电荷 |
| [14-dielectric-zno-hse-lcalceps](examples/14-dielectric-zno-hse-lcalceps/README.md) | 算例 14-3：ZnO 的 HSE06 介电张量与 Born 有效电荷 |
| [14-optics-si-loptics](examples/14-optics-si-loptics/README.md) | 算例 14-1：Si 的频率相关介电函数与光学常数 |
| [15-aimd-si-liquid-nvt](examples/15-aimd-si-liquid-nvt/README.md) | 算例 15-1：液态 Si 的从头算分子动力学 |
| [15-aimd-water-npt](examples/15-aimd-water-npt/README.md) | 算例 15-2：液态水的 NpT 分子动力学 |
| [15-mlff-si-train-run](examples/15-mlff-si-train-run/README.md) | 算例 15-3：VASP 机器学习力场全流程 |
| [16-bader-co-pt111](examples/16-bader-co-pt111/README.md) | 示例 16-1：CO/Pt |
| [16-cdd-co-pt111](examples/16-cdd-co-pt111/README.md) | 示例 16-2：CO/Pt |
| [16-elf-si](examples/16-elf-si/README.md) | 示例 16-3：金刚石 Si 的电子局域函数 |
| [16-lobster-gaas](examples/16-lobster-gaas/README.md) | 示例 16-5：GaAs 的 LOBSTER COHP / ICOHP 分析 |
| [16-stm-mos2](examples/16-stm-mos2/README.md) | 示例 16-4：单层 MoS2 的部分电荷 |
| [17-effmass-gaas](examples/17-effmass-gaas/README.md) | 示例 17-2：GaAs 在 Γ 点的有效质量 |
| [17-fermi-surface-cu](examples/17-fermi-surface-cu/README.md) | 示例 17-4：fcc Cu 的费米面 |
| [17-unfold-si-vacancy](examples/17-unfold-si-vacancy/README.md) | 示例 17-3：含空位的 Si 超胞能带反折叠 |
| [17-wannier90-si](examples/17-wannier90-si/README.md) | 示例 17-1：Si 的最大局域 Wannier 函数 |
| [18-postprocess-toolkit](examples/18-postprocess-toolkit/README.md) | 示例 18：后处理与自动化工具箱 |

## 辅助脚本

完整说明见 [scripts/README.md](scripts/README.md)。

- [check_vasp.py](scripts/check_vasp.py)
- [chgdiff.py](scripts/chgdiff.py)
- [conv_analyze.py](scripts/conv_analyze.py)
- [conv_test.sh](scripts/conv_test.sh)
- [elastic_post.py](scripts/elastic_post.py)
- [locpot_workfunction.py](scripts/locpot_workfunction.py)
- [md_analysis.py](scripts/md_analysis.py)
- [neb_tools.py](scripts/neb_tools.py)
- [plot_band.py](scripts/plot_band.py)
- [plot_dos.py](scripts/plot_dos.py)
- [vasp_incar_writer.py](scripts/vasp_incar_writer.py)
- [vasp_inputgen.py](scripts/vasp_inputgen.py)

## 相关资源

- [GitHub 相关仓库精选](../Resources/github-repos.md)
- [NEP/GPUMD 与数据集脚本](../Scripts/README.md)
- [返回项目首页](../README.md)
