# 第1章 OpenMX 概述、版本与学习路线（OpenMX Overview, Versions & Learning Path）

本章介绍开源第一性原理程序 OpenMX（Open source package for Material eXplorer）的基本原理与定位：它用数值赝原子轨道（PAO，LCAO 基组）展开 Kohn–Sham 轨道，配合模守恒赝势（VPS），在 O(N) 大体系、NEGF 电子输运、非共线磁性与自旋轨道耦合、Wannier 函数等方向很有特色。我们会梳理 OpenMX 4.0（2026-04）与 3.9 的版本差异和 4.0 的新功能（CWF、NPT、COHP、GCenter k 点、scf.restart pp、cif2omx、deb 包等），介绍发行包目录结构、单位约定、从结构到后处理的典型工作流（OpenMX workflow），并给出本教程 18 章的学习路线和资源清单。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

- 说清 OpenMX 的核心方法（PAO 基组 + 模守恒赝势 + 实空间网格），以及它和 VASP、Quantum ESPRESSO 等平面波程序在使用上的差别；
- 知道最新版本是 **OpenMX 4.0 + patch 4.0.1**，了解 4.0 新增了哪些功能、哪些默认行为与 3.9 不同（特别是 k 点生成方式）；
- 认识 `openmx4.0/` 下的 `source`、`work`、`DFT_DATA19` 三个目录，知道 PAO/VPS 文件放在哪里、怎么命名；
- 记住 OpenMX 输入/输出中各物理量的单位（Å/Bohr、Hartree/eV/Ry、Hartree/Bohr、GV/m 等）；
- 掌握“结构 → `.dat` 输入 → SCF → 后处理 → 可视化”的完整流程，以及每一步该用什么工具；
- 根据自己的研究方向，从本教程中挑出合适的学习路线。

| 你的情况 | 建议先读 |
|---|---|
| 完全没接触过 OpenMX，只想尽快跑通一个算例 | 本章 §4、§6 → [第2章](02-install.md) 的 deb 包或 Docker 安装 → [第3章](03-input-first-calc.md) |
| 用过 VASP/QE，想迁移到 OpenMX | 本章 §2.5（与平面波程序的差别）、§5（单位）→ [第4章 PAO/VPS](04-pao-vps.md) → [第5章 SCF/k 点](05-scf-kpoints-convergence.md) |
| 手上有 OpenMX 3.9 的旧结果，要升级到 4.0 | 本章 §3.3（行为变化），尤其是 `scf.Generation.Kpoint` |
| 课题是磁性/拓扑/输运/表面 | 本章 §8 的“按研究方向选章节”表 |
| Linux、集群、作业调度还不熟 | 先读公共篇 [第1章 计算环境搭建](../../Common/tutorials/01-linux-hpc.md) 与 第2章 HPC 作业调度 |

> [!NOTE]
> 本教程以 **OpenMX 4.0（2026-04-26 发布）+ patch 4.0.1（2026-05-08）** 为准，数据库为 **DFT_DATA19**（VPS/PAO Ver. 2019）。凡是 4.0 与 3.9 行为不同的地方都会标注。截至 2026-11，官方下载页仍提供 3.9（+ patch 3.9.9）与 3.8，供复现旧结果使用。

---

## 2. 原理速览：OpenMX 是怎样计算的

### 2.1 一句话定位

OpenMX 是由东京大学物性研究所（ISSP）尾崎泰助（T. Ozaki）课题组自 2000 年起开发的密度泛函理论（DFT）程序，源码以 **GPLv3** 许可证发布。它的三个核心技术是：

1. **数值赝原子轨道（pseudo-atomic orbitals, PAO）基组**：每个基函数都是“数值径向函数 × 实球谐函数”，在截断半径 $r_c$ 之外严格为零；
2. **模守恒赝势（VPS）**：采用 Morrison–Bylander–Kleinman（MBK）形式，每个角动量（或总角动量 $j$）可以有多个参考能量，DFT_DATA19 中的 VPS 都是全相对论（$j$ 依赖）的，可以直接用于自旋轨道耦合计算；
3. **实空间规则网格**：Hartree 势、交换关联势等在实空间网格上计算，Poisson 方程用 FFT 求解，网格密度由 `scf.energycutoff`（单位 Ry）控制。

### 2.2 PAO 基组与广义本征值问题

一个 PAO 基函数写成

$$
\chi_{i\alpha}(\mathbf r)=R_{i\alpha}(|\mathbf r-\mathbf R_i|)\,Y_{lm}(\widehat{\mathbf r-\mathbf R_i}),\qquad R_{i\alpha}(r)=0\ \ (r>r_c)
$$

Kohn–Sham 轨道是这些原子中心基函数的线性组合：

$$
\psi_{n\mathbf k}(\mathbf r)=\sum_{i\alpha} c_{n\mathbf k,i\alpha}\,\chi_{i\alpha}^{\mathbf k}(\mathbf r)
$$

由于 PAO 之间互不正交，求解的是**广义本征值问题**

$$
H(\mathbf k)\,\mathbf c_{n\mathbf k}=\varepsilon_{n\mathbf k}\,S(\mathbf k)\,\mathbf c_{n\mathbf k}
$$

其中 $S$ 是重叠矩阵。因为基函数有限程，$H$ 和 $S$ 都是稀疏矩阵，这正是 OpenMX 能做 O(N) 线性标度计算（DC、DC-LNO、Krylov 方法）的基础。

基组大小由两个参数控制：**截断半径**和**每个角动量的径向函数个数**。例如 `C6.0-s2p2d1` 表示用 `DFT_DATA19/PAO/C6.0.pao` 文件（截断半径 6.0 Bohr），s、p、d 分别取 2、2、1 个径向函数，共 $2\times1+2\times3+1\times5=13$ 个基函数。基组选择的细节见 [第4章](04-pao-vps.md)。

### 2.3 实空间网格：`scf.energycutoff` 不是“基组截断能”

在平面波程序中，ENCUT 决定基组大小，能量随 ENCUT 增大从上方单调收敛（变分性）。而 OpenMX 的 `scf.energycutoff` 只决定**数值积分与 Poisson 方程所用的实空间网格间距**，与基组无关。因此：

- 总能量随 `scf.energycutoff` 增大**不一定单调**变化（手册第 14 章明确指出这一点）；
- 官方建议多数情况下 150–200 Ry 已足够，但能量面很平坦的结构优化、不同磁序之间的小能量差等问题可能需要 300 Ry 以上；
- 4.0 手册的关键词说明给出的默认值是 220 Ry，而 3.9 手册写的是 150 Ry——**无论哪个版本，都建议在输入文件中显式写出 `scf.energycutoff`**。

### 2.4 模守恒赝势 VPS 与价电子数

每个 VPS 文件开头都写明了生成赝势时考虑的价电子数（例如 `C_PBE19.vps` 中的 `valence.electron 4.0`）。在输入文件的原子坐标行中，你需要给出每个原子的“初始上自旋电子数”和“初始下自旋电子数”，**两者之和应等于该 VPS 的价电子数**。同一元素可能有多个 VPS，例如 `Fe_PBE19H`（含更多半芯态，价电子 16）和 `Fe_PBE19S`（价电子 14），坐标行要与所选 VPS 一致。这一点初学者最容易出错，[第3章](03-input-first-calc.md) 和 [第4章](04-pao-vps.md) 会反复强调。

### 2.5 与平面波程序（VASP、QE）的主要区别

| 方面 | OpenMX（PAO + VPS） | 平面波程序（VASP/QE） |
|---|---|---|
| 基组 | 需要**为每种元素选择** PAO 文件（截断半径）和轨道数，Quick/Standard/Precise 三档推荐 | 只需收敛 ENCUT 一个参数 |
| 能量收敛性 | 基组不完备，存在基组叠加误差（BSSE），可用空原子（ghost atom，物种 `E`）做 counterpoise 修正 | 无 BSSE，ENCUT 增大时能量变分收敛 |
| 真空层 | 真空区域几乎不增加基函数，分子、表面、低维体系代价低 | 真空区域同样要用平面波填满，大盒子很贵 |
| 大体系 | 稀疏矩阵 + O(N) 方法（DC、DC-LNO、Krylov），可算上千至上万原子 | 通常 O(N³)，几百原子已很吃力 |
| 输运 | 内置 NEGF（非平衡格林函数）输运，含偏压、门压、本征通道 | 需要借助外部程序 |
| 局域基哈密顿量 | 直接输出 $H(\mathbf R)$、$S(\mathbf R)$（`HS.fileout on` → `.scfout`），便于做 Jx 交换耦合、紧束缚模型、DeepH/HamGNN 机器学习哈密顿量 | 一般需先构造 Wannier 函数 |
| 泛函 | `scf.XcType` 只有 `LDA`、`LSDA-CA`、`LSDA-PW`、`GGA-PBE`；可加 DFT+U、DFT-D2/D3。**4.0 手册不提供杂化泛函和 meta-GGA** | 杂化泛函、meta-GGA、RPA 等更全面 |
| 声子 | 没有 DFPT；声子需用有限位移法外部程序（如带 OpenMX 接口的 ALAMODE） | VASP 有 DFPT/有限位移，QE 有 ph.x |
| 过完备问题 | 致密金属（fcc/bcc/hcp）用很大的基组时可能出现伪本征值，需减小基组或关闭 `scf.ProExpn.VNA` | 不存在 |

> [!TIP]
> 一个简单的选择标准：**体系大、有真空、要做输运或局域基哈密顿量后处理、要做非共线磁性/SOC 的批量计算**时，OpenMX 很有优势；**需要杂化泛函精确带隙、DFPT 声子、或者与大量已发表 VASP 数据严格对比**时，平面波程序更合适（见本仓库 [VASP 篇](../../VASP/README.md)）。

### 2.6 OpenMX 能做什么（功能地图）

| 类别 | 主要功能 | 本教程章节 |
|---|---|---|
| 基态 | 总能量与力、分子（Cluster）/周期（Band）计算、O(N) 方法 | [第3章](03-input-first-calc.md)、[第5章](05-scf-kpoints-convergence.md)、第18章 |
| 结构 | 结构优化（EF/BFGS/RF/DIIS）、晶胞优化（OptC*/RFC5）、E–V 状态方程、NEB、DFT-D2/D3 | [第6章](06-optimization-eos-neb.md) |
| 动力学 | NVE、NVT（速度标定、Nosé–Hoover）、4.0 新增 NPT | [第7章](07-md.md) |
| 电子结构 | 能带、能带展开（unfolding）、DOS/PDOS、COHP/COOP、光学电导、BoltzTraP/BoltzTraP2 接口 | [第8章](08-band.md)、[第9章](09-dos-pdos-cohp-optics.md) |
| 电荷分析 | Mulliken、Voronoi、ESP 电荷、NPA、分子轨道与 cube 输出 | [第10章](10-charge-mo-cube.md) |
| 磁性 | 共线/非共线自旋、约束 DFT、MAE、自旋螺旋、交换耦合 Jx | [第11章](11-spin-noncollinear.md) |
| 相对论与拓扑 | SOC、自旋织构、Chern 数、Z2 不变量 | [第12章](12-soc-topology.md) |
| 强关联 | DFT+U（Dudarev 与一般形式，含 Hund J） | [第13章](13-dft-u.md) |
| Wannier | MLWF、4.0 新增最接近 Wannier 函数（CWF）、Wannier90 接口 | 第14章 |
| 外场与极化 | Berry 相宏观极化、锯齿形均匀电场 | 第15章 |
| 表面 | ESM 方法、功函数、电荷掺杂、STM、XPS 芯能级（4.0 可与 ESM 联用） | [第16章](16-surface-esm-workfunction.md) |
| 输运 | NEGF 电子输运（透射、I–V、本征通道、电流密度） | [第17章](17-negf-transport.md) |

---

## 3. 版本与新功能

### 3.1 版本时间线

| 版本 | 发布日期 | 说明 |
|---|---|---|
| 3.8 + patch 3.8.5 | 2016-04-03 / 2018-06-12 | 旧版，仍可下载 |
| 3.9 | 2019-12-03 | 引入数据库 Ver. 2019（DFT_DATA19）、DC-LNO 方法、XPS 芯能级；起必须链接 ScaLAPACK + BLACS + MPI |
| patch 3.9.9 | 2021-10-17 | 3.9 的最后一个补丁 |
| **4.0** | **2026-04-26** | 当前最新版（源码包约 182 MB），附 4.0 手册 |
| **patch 4.0.1** | **2026-05-08** | 修正 bxsf/frmsf 单位、非共线 Γ 点 Mulliken 的 HOMO 索引、`work/GaAs.dat`；同时更新 Debian 包 |
| OpenMX Viewer 大更新 | 2026-05-26 | 浏览器端前后处理工具，支持更多格式、Materials Project/COD 直取、Structure/Slab Builder、离线 Windows 版 |
| 4.0 GPU v1.1.1 | 2026-08-25 | 由 dc1394 维护的 GPU 版（CMake 构建，CUDA/HIP/CPU 后端），官方下载页给出链接 |

### 3.2 OpenMX 4.0 新功能一览

下表根据官方 What's new 页面与 4.0 手册整理（关键词以手册为准，具体用法见对应章节）：

| 新功能 | 关键词 / 工具 | 用途 | 章节 |
|---|---|---|---|
| 最接近 Wannier 函数（CWF） | `CWF.Calc`、`CWF.Guiding.Orbitals` 等 | 以原子轨道/杂化轨道/分子轨道为引导构造局域函数，适合成键分析与紧束缚建模 | 第14章 |
| NPT 分子动力学 | `MD.Type NPT_VS_PR`、`NPT_VS_WV` | 速度标定 + Parrinello–Rahman / Wentzcovitch 变胞 MD | [第7章](07-md.md) |
| COHP/COOP | 后处理程序 `cohp` | 分析成键/反键贡献 | [第9章](09-dos-pdos-cohp-optics.md) |
| BoltzTraP2 接口 | — | 半经典玻尔兹曼输运、热电系数 | [第9章](09-dos-pdos-cohp-optics.md) |
| 含 Γ 与 TRIM 点的 k 点 | `scf.Generation.Kpoint GCenter`（默认） | Γ 中心网格，偶数分割也包含 Γ 和时间反演不变动量点 | [第5章](05-scf-kpoints-convergence.md) |
| 受约束晶胞优化 | `MD.Type RFC5` + `<MD.Fixed.Cell.Vectors` | 统一实现 RFC6/RFC7 那样的晶格约束 | [第6章](06-optimization-eos-neb.md) |
| 后处理重启 | `scf.restart pp` | 读重启文件只做一次对角化，做 DOS/能带等后处理，不覆盖重启文件 | [第5章](05-scf-kpoints-convergence.md) |
| 共线→非共线重启 | `scf.restart c2n` | 二阶变分法计算磁各向异性能（MAE） | [第11章](11-spin-noncollinear.md) |
| CIF → dat | `cif2omx` | 命令行把 CIF 转成 OpenMX 输入，见本章 §6.3 | 本章 |
| XPS + ESM | — | 表面/板模型的芯能级绝对结合能 | [第16章](16-surface-esm-workfunction.md) |
| 区域化学势（RCP）分析 | 配合 FLPQ | 表面电子给予能力、吸附位点选择性 | [第16章](16-surface-esm-workfunction.md) |
| OpenMX API | C/Fortran/Python 调用 | 外部程序把 OpenMX 当作库调用 | 第18章 |
| ELPA 开发者接口 | `elpa_col`、`elpa_noncol` | 读 `.scfout` 用 ELPA1/ELPA2 并行对角化，供后处理开发者使用 | 第18章 |
| 能带展开改进 | `Unfolding.Electronic.Band` | 内存更省、速度更快 | [第8章](08-band.md) |
| md 文件格式文档化 | `System.Name.md` | 与 dpdata 的 `openmx/md` 格式一致，可转 DeePMD 训练集 | [第7章](07-md.md) |
| Debian 包 | `openmx_4.0.1-1_amd64.deb` 等 | Ubuntu/Debian/WSL 免编译安装 | [第2章](02-install.md) |
| 内存泄漏检测 | `-DLEAK_DETECT`、`leakdetect.c/h` | 面向开发者 | [第2章](02-install.md) |

### 3.3 从 3.9 迁移到 4.0：必须知道的行为变化

| 项目 | 3.9 | 4.0 | 对你的影响 |
|---|---|---|---|
| k 点生成 | 规则网格（相当于 4.0 的 `Regular`），只有各方向分割数都为奇数时才含 Γ | 默认 `scf.Generation.Kpoint GCenter`，总是包含 Γ 与 TRIM 点 | **同一个 `scf.Kgrid` 在两个版本下用的 k 点不同**，总能量会有差别；要复现 3.9 结果，显式写 `scf.Generation.Kpoint Regular` |
| 自旋螺旋（广义 Bloch 定理） | 正常 | 论坛报告与 GCenter 联用时出现 MPI 错误 | 自旋螺旋计算设 `scf.Generation.Kpoint Regular`（见 [第11章](11-spin-noncollinear.md)） |
| `scf.energycutoff` 默认值 | 手册写 150 Ry | 手册关键词说明写 220 Ry（而 `work/Methane.dat` 注释仍写 150） | 一律显式设置 |
| 本征值求解库 | ELPA | 能带计算默认 `scf.eigen.lib elpa2`，分子（Cluster）默认 `elpa1` | 一般无需修改；ELPA 源码随 OpenMX 内置 |
| `scf.restart` | `on`/`off` | 新增 `pp`、`c2n` | 后处理更方便 |
| CDDF 相关关键词 | 有 `CDDF.material_type` | 新增 `CDDF.shift_ChemP`，删除 `CDDF.material_type` | 旧光学输入需要修改 |
| Sn 推荐基组 | 旧推荐 | Quick/Standard/Precise 改为 `Sn7.0-s2p2d2` / `Sn7.0-s3p2d3` / `Sn7.0-s3p2d3f1` | 含 Sn 体系注意更新 |
| 默认 makefile | 面向旧 Intel 编译器（mpicc + icc） | 面向 WSL + Intel oneAPI（`mpiicx`/`mpiifx`） | 编译方法见 [第2章](02-install.md) |
| 串行版本 | — | 4.0 **不支持**纯串行安装，MPI 与 ScaLAPACK/BLACS（或 oneAPI MKL）为必需 | 个人电脑也要装 MPI |
| 数据库 | DFT_DATA19 | 仍为 DFT_DATA19，另增 Tb、Yb 数据与 XPS 芯能级 VPS/PAO | 3.9 的 `.dat` 中基组名一般可直接沿用 |

> [!WARNING]
> 拿 4.0 的结果和 3.9（或文献中用 3.9 算出的）数值逐位比较之前，先确认 k 点生成方式一致。周期体系的总能量差异常常就来自 `GCenter` 与 `Regular` 的不同，而不是程序 bug。

### 3.4 patch 4.0.1 修正了什么

| 修改的文件 | 问题 | 影响 |
|---|---|---|
| `Band_DFT_Dosout.c` | bxsf/frmsf 费米面文件的单位原为 Hartree/Bohr，按 XCrySDen 惯例改为 eV/Å | **打补丁前后生成的费米面文件单位不同**，旧脚本若仍按 Hartree/Bohr 处理会出错 |
| `Mulliken_Charge.c` | 非共线计算只用 Γ 点时 HOMO 索引被错误覆盖 | 影响非共线 Γ 点计算输出的 HOMO 信息 |
| `GaAs.dat` | `work/GaAs.dat` 无法复现手册第 38.1 节的结果 | 补丁中附带修正后的 `GaAs.dat`，需移动到 `work/` |

打补丁的具体命令见 [第2章 §3](02-install.md)。

### 3.5 如何确认自己在用哪个版本

- 运行任意算例，标准输出（以及 `System.Name.out`）开头的欢迎横幅会打印版本号，例如 `Welcome to OpenMX   Ver. 4.0`（示意，具体排版以实际输出为准）；
- 源码目录中 `source/openmx.c` 文件头注释写有版本与日期，可用 `head -n 20 source/openmx.c | grep -i "ver"` 查看；
- 补丁不会改版本号，请自己记录是否打过 patch 4.0.1（例如在安装目录下保留补丁包和一份 `INSTALL_NOTES.txt`）。

---

## 4. 发行包目录结构

### 4.1 源码包解压后的三个目录

```text
openmx4.0/
├── source/            # C 与 Fortran 源码、makefile、内置 ELPA（elpa-2018.05.001 目录及 elpa1.f90 等）
│                      # 以及后处理程序源码：DosMain.c、bandgnu13.c、polB.c、cohp.c、cif2omx.c、gcube2oned.c ...
├── work/              # 运行目录：示例输入 *.dat；编译后 openmx 与后处理程序会被复制到这里
│   ├── Methane.dat, H2O.dat, Cdia.dat, GaAs.dat, ...   # 手册中各节对应的示例
│   ├── input_example/      # -runtest 自动测试的 14 个输入及参考输出
│   ├── large_example/      # -runtestL 大体系基准测试
│   └── negf_example/, wf_example/, cwf_example/, cddf_example/, unfolding_example/, ...  # 各功能专项示例
└── DFT_DATA19/
    ├── PAO/           # 数值赝原子轨道，*.pao，例如 C6.0.pao、Fe5.5H.pao、Fe6.0S.pao
    └── VPS/           # 赝势，*.vps，例如 C_PBE19.vps（GGA）、C_CA19.vps（LDA）、Fe_PBE19H.vps
```

（`work/` 下子目录以实际解压结果为准。）

> [!IMPORTANT]
> 关键词 `DATA.PATH` 的默认值是**相对路径** `../DFT_DATA19`。这意味着如果不写 `DATA.PATH`，只有在 `openmx4.0/work/` 目录下运行，OpenMX 才能找到 PAO 和 VPS 文件。在自己的工作目录中计算时，请写**绝对路径**，并指向**同时包含 `PAO/` 和 `VPS/` 的那一层目录**：
>
> ```text
> DATA.PATH   /home/user/software/openmx4.0/DFT_DATA19
> ```
>
> 不要写到 `.../DFT_DATA19/PAO`。

### 4.2 Debian 包安装后的布局

用官方 deb 包安装（见 [第2章 §8](02-install.md)）时，可执行文件和工具位于 `/opt/openmx4.0` 下，`openmx`、`DosMain` 等可以直接用命令名调用，示例输入在 `/opt/openmx4.0/work/`。数据库目录通常是 `/opt/openmx4.0/DFT_DATA19`，请用 `ls /opt/openmx4.0` 确认后在 `DATA.PATH` 中写绝对路径。

### 4.3 PAO 与 VPS 文件名速读

| 文件名 | 含义 |
|---|---|
| `C6.0.pao` | 碳的 PAO 文件，截断半径 6.0 Bohr；在输入中写作 `C6.0-s2p2d1` 等（去掉扩展名，加轨道指定） |
| `Fe5.5H.pao`、`Fe6.0S.pao` | 与 `Fe_PBE19H`（更多半芯态）、`Fe_PBE19S` 两种 VPS 配套的 PAO |
| `C_PBE19.vps` | 碳的 GGA-PBE 赝势，DFT_DATA19 版本 |
| `C_CA19.vps` | 碳的 LDA（Ceperley–Alder）赝势 |
| `E.vps` 与 `Kr10.0-s1p1` 等 | 空原子（ghost atom，物种名常写作 `E`），价电子 0，用于 BSSE 修正或在真空区补充基函数 |

想快速查看某个 VPS 的价电子数，可以直接 `grep`：

```bash
grep valence.electron /path/to/openmx4.0/DFT_DATA19/VPS/C_PBE19.vps
grep valence.electron /path/to/openmx4.0/DFT_DATA19/VPS/Fe_PBE19H.vps /path/to/openmx4.0/DFT_DATA19/VPS/Fe_PBE19S.vps
```

本章配套示例 [01-overview-tour](../examples/01-overview-tour/README.md) 提供了一个 Python 小工具 `list_dft_data.py`，可以列出数据库中所有 VPS 的价电子数、PAO 的截断半径与可用轨道数，并检查 `.dat` 文件中物种定义和初始电子数是否自洽。

### 4.4 `work/` 中值得先看的示例输入

| 文件 | 内容 | 对应章节 |
|---|---|---|
| `Methane.dat` | 甲烷分子 SCF（安装验证的标准算例） | [第2章](02-install.md)、[第3章](03-input-first-calc.md) |
| `H2O.dat`、`Methane2.dat` | 水分子、畸变甲烷的结构优化 | [第6章](06-optimization-eos-neb.md) |
| `Cdia.dat` | 金刚石体相 SCF、能带、DOS | [第8章](08-band.md)、[第9章](09-dos-pdos-cohp-optics.md) |
| `Si.dat` | 由 cif2omx 从 CIF 生成的示例输入 | 本章 §6.3 |
| `Febcc2.dat`、`Cr2_CNC.dat`、`Mol_MnO_NC.dat` | 铁磁 bcc Fe、约束非共线 Cr2、非共线 MnO 分子 | [第11章](11-spin-noncollinear.md) |
| `GaAs.dat`、`Bi2Se3-Z2.dat`、`Graphene-Chern.dat` | SOC、Z2、Chern 数 | [第12章](12-soc-topology.md) |
| `NiO-cFLL.dat`、`NiO-sFLL.dat` | DFT+U | [第13章](13-dft-u.md) |
| `GaAs-CWF1.dat`、`Si-Wannier90.dat` | CWF 与 Wannier90 接口 | 第14章 |
| `NaCl.dat`、`Nitro_Benzene.dat` | Berry 相极化、外加电场 | 第15章 |
| `Al111_WorkFunc_*.dat`、`Al-Si111_ESM.dat` | 功函数、ESM | [第16章](16-surface-esm-workfunction.md) |
| `C2H4_NEB.dat`、`Gly_NH.dat`、`Si_NPT_VS_PR.dat` | NEB、Nosé–Hoover MD、NPT MD | [第6章](06-optimization-eos-neb.md)、[第7章](07-md.md) |
| `Sidia_COHP.dat` | COHP/COOP | [第9章](09-dos-pdos-cohp-optics.md) |
| `DIA512_DC.dat`、`C60_DC.dat` | O(N) 分治法 | 第18章 |

完整列表见手册第 78 节 “Included input examples”。

> [!TIP]
> 官方示例是最可靠的关键词写法参考。遇到不确定的关键词时，先在 `work/` 下搜索：`grep -il "scf.Generation.Kpoint" *.dat`、`grep -i -A3 "<Band.kpath" Cdia.dat`。

---

## 5. 单位约定

OpenMX 的输入单位比较“混搭”：长度可选 Å 或原子单位，截断能用 Ry，温度用 K，收敛判据用 Hartree，DOS 范围用 eV，电场用 GV/m。下表按关键词整理（以 4.0 手册为准）：

| 物理量 | 关键词 / 输出位置 | 单位 |
|---|---|---|
| 原子坐标 | `Atoms.SpeciesAndCoordinates.Unit` | `Ang`（Å）、`AU`（Bohr）、`FRAC`（分数坐标，以 a、b、c 为基） |
| 晶格矢量 | `Atoms.UnitVectors.Unit` | `Ang` 或 `AU` |
| 网格截断能 | `scf.energycutoff` | Ry（Rydberg） |
| 电子温度 | `scf.ElectronicTemperature` | K |
| SCF 收敛判据 | `scf.criterion` | Hartree（比较的是相邻两步能带能量之差 dUele） |
| 结构优化判据 | `MD.Opt.criterion` | Hartree/Bohr（最大原子力，默认 3×10⁻⁴） |
| MD 时间步长 | `MD.TimeStep` | fs |
| 初速度 | `<MD.Init.Velocity` | m/s |
| DOS 能量范围 | `Dos.Erange` | eV（相对化学势） |
| 能带路径 | `<Band.kpath` | 以倒格矢为单位的分数坐标 |
| 外加电场 | `scf.Electric.Field` | GV/m（1 V/Å = 10 GV/m），符号按作用在电子上的方向定义 |
| Hubbard U、Hund J | `<Hubbard.U.values`、`<Hund.J.values` | eV |
| 体系电荷 | `scf.system.charge` | 电子电荷，正值表示空穴掺杂、负值表示电子掺杂 |
| 总能量及各分量、本征值、化学势 | `System.Name.out` | Hartree |
| 原子力 | `System.Name.out` 中 `<coordinates.forces` 段 | Hartree/Bohr（坐标为 Å） |
| 偶极矩 | `System.Name.out` | Debye |
| 原始能带数据 | `System.Name.Band` | Hartree（未平移） |
| gnuplot 能带数据 | `bandgnu13` 生成的 `.BANDDAT1/2` | eV，以化学势为零点 |
| 费米面文件 | `.bxsf`、`.frmsf` | patch 4.0.1 之后为 eV 与 Å（之前为 Hartree 与 Bohr） |
| cube 文件 | `*.tden.cube`、`*.v0.cube` 等 | Gaussian cube 格式（网格坐标以 Bohr 记） |

常用换算（CODATA 2018 数值，保留到够用的位数）：

| 换算 | 数值 |
|---|---|
| 1 Hartree | 27.211386 eV = 2 Ry |
| 1 Ry | 13.605693 eV |
| 1 Bohr | 0.529177 Å |
| 1 Hartree/Bohr | 51.422 eV/Å（所以默认 `MD.Opt.criterion` 3×10⁻⁴ Ha/Bohr ≈ 0.015 eV/Å） |
| 1 meV/atom | 3.6749×10⁻⁵ Hartree/atom |
| 1 Debye | 0.393430 e·Bohr |

示例目录中的 `omx_units.py` 可以做这些换算，例如把 `.out` 中的 Hartree 能量差换成 meV。

---

## 6. 典型工作流

### 6.1 从结构到结果的流程图

```text
 结构来源                     准备输入                      计算                       后处理与可视化
 ─────────                   ─────────                    ──────                     ──────────────
 Materials Project / COD     OpenMX Viewer（导出 .dat）    mpirun -np N openmx         bandgnu13      → 能带
 实验 CIF / POSCAR     ──►   cif2omx（命令行）       ──►   x.dat -nt M | tee x.std ──► DosMain        → DOS/PDOS
 自己建模（ASE 等）          ASE / 自写脚本               │                           cohp           → COHP/COOP
                             + 手工补充关键词              ├─ x.out  能量/力/Mulliken  polB / jx      → 极化 / 交换耦合
                                                          ├─ x.std  SCF 过程          Z2FH/calB/kSpin→ 拓扑
                                                          ├─ x.xyz / x.md 结构与轨迹   gcube2oned     → 平面平均势/功函数
                                                          ├─ *.cube 密度、势、MO      diff_gcube     → 差分密度
                                                          └─ x_rst/ 重启文件          VESTA / XCrySDen / OpenMX Viewer
```

一个标准的研究计算通常是这样的顺序：

1. **准备结构**：从数据库或实验 CIF 获取晶体结构，或自己搭建分子/表面模型（通用建模方法见公共篇 第3章 结构准备）。
2. **生成 `.dat`**：用 OpenMX Viewer 或 `cif2omx` 生成初稿，再手工补充/检查关键词：`DATA.PATH`、基组、初始自旋、`scf.energycutoff`、`scf.Kgrid` 等。本仓库的 `omx_make_input.py` 也可以从 POSCAR/CIF 生成输入（用法见 `python ../scripts/omx_make_input.py -h`）。
3. **收敛测试**：基组、截断能、k 点（[第4章](04-pao-vps.md)、[第5章](05-scf-kpoints-convergence.md)），可用 `omx_conv_test.py` 批量完成。
4. **结构优化**（[第6章](06-optimization-eos-neb.md)）→ 用优化后结构做**高精度 SCF**（打开 `HS.fileout`、`Dos.fileout`、`Band.dispersion` 等需要的输出）。
5. **后处理**：能带、DOS、电荷、极化、交换耦合、Wannier、输运等。
6. **检查与可视化**：确认 SCF 收敛、Mulliken 电荷/磁矩合理，再画图。

### 6.2 OpenMX Viewer：浏览器里的前后处理入口

[OpenMX Viewer](https://www.openmx-square.org/viewer/) 是官方提供的浏览器端工具，4.0 手册第 10 节把它定位为“连接前处理与后处理的实用界面”。2026-05-26 的大版本更新后，它支持：

- 读入 OpenMX 的 `.dat`、`.xyz`、`.md`、`.cube`、`.Band`，以及 CIF、XSF、POSCAR/CONTCAR、XDATCAR、OUTCAR、Quantum ESPRESSO 输入/输出、BXSF、FRMSF 等格式；
- 输入 Materials Project 或 COD 编号直接加载晶体结构，无需手动下载 CIF；
- 基于 spglib 与 moyo（WebAssembly）的对称性分析，原胞/惯用胞切换；
- Structure Builder（选中、平移、旋转、替换、增删原子，最多 100 步撤销）和 Slab Builder（按密勒指数、层厚、真空厚度和端面生成板模型）；
- 自动给出能带计算的建议 k 路径，导出 OpenMX 输入文件和高分辨率 PNG；
- cube/XSF 等值面与任意截面、费米面显示；
- 提供可下载的 Windows 离线版。拖入浏览器的文件在本地处理，官方说明不会上传到服务器。

> [!TIP]
> Viewer 更新后如果页面行为异常，先清除浏览器缓存或强制刷新（官方提示旧 JavaScript 缓存会导致问题）。

### 6.3 cif2omx：命令行把 CIF 转成 `.dat`

4.0 起，`cif2omx` 作为标准工具写入手册（第 8.3 节）。`make all` 时会一并编译；也可以在 `source/` 中按手册单独编译：

```bash
cd openmx4.0/source
gcc -std=c99 -O2 -Wall -Wextra -pedantic -DCIF2OPENMX_MAIN cif2omx.c -lm -o cif2omx
```

用法（选项写在 CIF 文件名之后）：

```bash
./cif2omx Si_mp149.cif -o Si.dat                     # 默认：Standard 基组、xyz 坐标、k 间距 0.33 1/Å、220 Ry
./cif2omx Si_mp149.cif -o Si.dat -b precise -c fractional -n si -k 0.20 -e 300
./cif2omx --help                                     # 查看全部选项
```

| 选项 | 含义 | 默认 |
|---|---|---|
| `-o, --output FILE` | 输出的 dat 文件名 | 标准输出 |
| `-c, --coord MODE` | 坐标模式 `xyz` 或 `fractional` | `xyz` |
| `-b, --basis MODE` | 基组档次 `quick`、`standard`、`precise` | `standard` |
| `-n, --name NAME` | 写入 `System.Name` 的值 | CIF 文件名主干 |
| `-k, --kspacing D` | 按倒空间间距 D（1/Å，含 2π）设 `scf.Kgrid`：$N_i=\lfloor b_i/D\rfloor$（$b_i$ 为第 $i$ 个倒格矢的长度），至少为 1 | 0.33（与 Viewer 一致） |
| `-e, --ecut VALUE` | 写入 `scf.energycutoff` | 220.0 Ry |

> [!IMPORTANT]
> cif2omx 生成的是“最小可用”的关键词集合。运行前务必检查并修改：`DATA.PATH`（改成绝对路径）、磁性原子的初始上/下自旋电子数、`scf.SpinPolarization`、k 点密度、`scf.criterion` 等。把生成的 `.dat` 拖进 OpenMX Viewer 看一眼原子排布和周期性，是很好的习惯。

### 6.4 计算与常用后处理程序

计算命令的通用形式（MPI × OpenMP 混合并行，`-nt` 是每个 MPI 进程的 OpenMP 线程数）：

```bash
mpirun -np 4 openmx Methane.dat -nt 1 | tee met.std
```

`make all` 会编译下列后处理程序并复制到 `work/`（完整列表以 `source/makefile` 中的 `UTIL` 变量为准）：

| 程序 | 输入 | 用途 | 章节 |
|---|---|---|---|
| `bandgnu13` | `System.Name.Band` | 转成 gnuplot 能带数据（eV，化学势为零点） | [第8章](08-band.md) |
| `DosMain` | `.Dos.val` + `.Dos.vec` | 四面体法/高斯展宽 DOS 与 PDOS | [第9章](09-dos-pdos-cohp-optics.md) |
| `cohp` | SCF 时写出的 `System.Name.cohp` 等 DOS 相关文件 | COHP/COOP | [第9章](09-dos-pdos-cohp-optics.md) |
| `polB` | `.scfout`（`HS.fileout on`） | Berry 相宏观极化 | 第15章 |
| `jx` | `.scfout` | 交换耦合常数 | [第11章](11-spin-noncollinear.md) |
| `kSpin`、`calB`、`Z2FH` | `.scfout` | 自旋织构、Chern 数/Berry 曲率、Z2 不变量 | [第12章](12-soc-topology.md) |
| `gcube2oned`、`diff_gcube`、`add_gcube` | cube 文件 | 平面平均势（功函数）、密度差、密度和 | [第10章](10-charge-mo-cube.md)、[第16章](16-surface-esm-workfunction.md) |
| `cube2xsf`、`md2axsf`、`bin2txt` | cube、md、二进制输出 | 格式转换（XCrySDen 动画等） | [第7章](07-md.md)、[第10章](10-charge-mo-cube.md) |
| `esp` | 计算结果 | ESP 拟合电荷 | [第10章](10-charge-mo-cube.md) |
| `analysis_example` | `.scfout` | 读取哈密顿量/重叠矩阵的示例程序，开发后处理代码的起点 | 第18章 |

本仓库 `OpenMX/scripts/` 中还有一组 Python 脚本，用于解析 `.out`（`omx_parse_out.py`）、画能带（`omx_band_plot.py`）、自动驱动 DosMain（`omx_dos_auto.py`）、处理轨迹（`omx_traj.py`）和 cube（`omx_cube_tools.py`）等，用 `python ../scripts/<脚本名> -h` 查看用法。

---

## 7. 生态：与 OpenMX 配合使用的工具

| 工具 | 作用 | 备注 |
|---|---|---|
| [ASE](https://ase-lib.org/ase/calculators/openmx.html) | Python 中调用 OpenMX、读写结构 | 官方手册第 82 节列出的 GUI/作业管理工具之一 |
| [sisl](https://sisl.readthedocs.io/en/latest/api/io/openmx.html) | 读取 OpenMX 输出、紧束缚与输运分析 | Python |
| [VESTA](https://jp-minerals.org/vesta/en/)、[XCrySDen](https://www.xcrysden.org/) | 结构、cube 等值面、非共线自旋矢量（`*.nc.xsf`） | 桌面软件 |
| [FermiSurfer](https://mitsuaki1987.github.io/fermisurfer/) | 费米面可视化 | 读 frmsf |
| [TB2J](https://tb2j.readthedocs.io/en/latest/src/openmx.html) / TB2J-OpenMX | 用磁力定理从 `.scfout` 计算海森堡交换 J、DMI | [第11章](11-spin-noncollinear.md) |
| DeepH / DeepH-E3 / HamGNN | 以 OpenMX 哈密顿量为训练数据的机器学习哈密顿量 | 第18章 |
| dpdata | 把 OpenMX MD 输出转成 DeePMD-kit 训练数据 | [第7章](07-md.md) |
| ALAMODE | 声子、非谐、热导（有 OpenMX 接口） | 弥补 OpenMX 没有 DFPT |
| Wannier90、BoltzTraP2 | Wannier 插值、热电输运 | 第14章、[第9章](09-dos-pdos-cohp-optics.md) |
| CrySPY、abICS、cif2x、AiiDA 插件 | 结构预测、构型采样、批量输入生成、高通量工作流 | 公共篇 第6章 |
| MateriApps LIVE! | 预装 OpenMX 的 Linux 环境（VirtualBox/Docker 版） | 免安装体验，见 [第2章](02-install.md) |

ASE、pymatgen 等通用 Python 库的使用方法见公共篇 第4章；VESTA 等可视化见公共篇 第5章。

---

## 8. 学习路线

### 8.1 推荐的阅读顺序

| 阶段 | 章节 | 你将学会 | 参考用时 |
|---|---|---|---|
| 0. 预备 | 公共篇第1–2章 | Linux、conda、集群登录、作业调度 | 视基础而定 |
| 1. 入门 | [第2章 编译与安装](02-install.md) → [第3章 第一个计算](03-input-first-calc.md) | 装好 OpenMX，读懂 `.dat`、`.std`、`.out` | 1–3 天 |
| 2. 打基础 | [第4章 PAO/VPS](04-pao-vps.md) → [第5章 SCF/k 点/收敛](05-scf-kpoints-convergence.md) | 选基组、做收敛测试、解决 SCF 不收敛 | 1 周 |
| 3. 常规性质 | [第6章 优化/EOS/NEB](06-optimization-eos-neb.md)、[第7章 MD](07-md.md)、[第8章 能带](08-band.md)、[第9章 DOS/COHP/光学](09-dos-pdos-cohp-optics.md)、[第10章 电荷/MO/cube](10-charge-mo-cube.md) | 完成一篇常规计算论文所需的大部分内容 | 2–3 周 |
| 4. 磁性与相对论 | [第11章 自旋/非共线](11-spin-noncollinear.md)、[第12章 SOC/拓扑](12-soc-topology.md)、[第13章 DFT+U](13-dft-u.md) | 磁基态、MAE、交换耦合、拓扑不变量、强关联氧化物 | 2 周 |
| 5. 专题 | 第14章 Wannier、第15章 Berry 相/电场、[第16章 表面/ESM](16-surface-esm-workfunction.md)、[第17章 NEGF 输运](17-negf-transport.md) | 按课题选读 | 按需 |
| 6. 生产计算 | 第18章 并行/O(N)/排错 | 大体系、作业脚本、常见错误速查 | 随用随查 |

### 8.2 按研究方向选章节

| 研究方向 | 必读 | 选读 |
|---|---|---|
| 分子与团簇 | 3、4、6、10 | 7（MD）、15（电场） |
| 半导体能带与缺陷 | 4、5、6、8、9 | 12（SOC）、8 中的能带展开 |
| 磁性材料 | 4、5、11、13 | 12（MAE/SOC）、TB2J |
| 拓扑材料 | 5、8、12、14 | 9（费米面/光学） |
| 表面与催化 | 6、10、16 | 7（MD）、9（COHP） |
| 分子器件与输运 | 8、9、17 | 16（ESM、门压概念） |
| 铁电/介电 | 6、15 | 9（光学） |
| 机器学习势/哈密顿量 | 7（dpdata）、10（`HS.fileout`）、18 | 公共篇第7章 |

### 8.3 学习建议

1. **先跑通，再改动**：先原样运行 `work/Methane.dat` 和本教程示例，确认输出与参考一致，再一次只改一个参数。
2. **每个新体系都做收敛测试**：基组（Quick/Standard/Precise）、`scf.energycutoff`、`scf.Kgrid`，判据一般取 1 meV/atom 量级。
3. **看 `.out` 而不只是看总能**：SCF 是否在 `scf.maxIter` 之内收敛、Mulliken 总电子数是否等于理想值、磁矩是否合理、力是否对称。
4. **记录版本与设置**：4.0/3.9、是否打补丁、`scf.Generation.Kpoint`、基组与 VPS 名称，都应写进你的计算记录和论文方法部分。
5. **善用论坛**：[OpenMX Forum](https://www.openmx-square.org/forum/patio.cgi) 上有开发者亲自回答的大量问题，提问前先搜索。

---

## 9. 资源汇总

### 9.1 官方资源

| 资源 | 链接 | 说明 |
|---|---|---|
| 官网 | https://www.openmx-square.org/ | 新闻、下载、数据库、论坛入口 |
| 下载页 | https://www.openmx-square.org/download.html | 4.0 + patch 4.0.1、3.9 + patch 3.9.9、3.8、GPU 版、deb、Docker |
| What's new | https://www.openmx-square.org/whatsnew/whatsnew.html | 4.0 新功能、patch、Viewer 更新 |
| 4.0 手册（英文） | https://www.openmx-square.org/openmx_man4.0/index.html | 本教程最主要的依据 |
| 4.0 手册（日文） | https://www.openmx-square.org/openmx_man4.0jp/ | — |
| 3.9 手册 | https://www.openmx-square.org/openmx_man3.9/index.html | 复现旧结果时对照 |
| 数据库 Ver. 2019 | https://www.openmx-square.org/vps_pao2019/ | VPS/PAO、各元素基准测试与推荐基组 |
| 芯能级数据库 | https://www.openmx-square.org/vps_pao_core2019/ | XPS 计算用 VPS/PAO |
| 技术笔记 | https://www.openmx-square.org/tech_notes/tech_notes.html | 总能量与力等推导 |
| OpenMX Viewer | https://www.openmx-square.org/viewer/ | 建模、k 路径、可视化 |
| 论坛 | https://www.openmx-square.org/forum/patio.cgi | 开发者答疑 |
| 视频讲座 | https://www.openmx-square.org/video_lec/video_lec.html | — |
| 研讨会资料 | https://www.openmx-square.org/workshop/workshop.html | hands-on 讲义与输入 |
| Ozaki 实用指南（2023 讲义 PDF） | https://t-ozaki.issp.u-tokyo.ac.jp/lectures-2023/OpenMX_Practical_Guide.pdf | 参数选择经验总结 |
| ADPACK | https://www.openmx-square.org/adpack_man2.2/adpack2_2.html | 自己生成 VPS/PAO 的原子计算程序（高级） |

### 9.2 社区教程

- [pranabdas/openmx](https://github.com/pranabdas/openmx) 及其站点 <https://pranabdas.github.io/openmx/>：目前最完整的英文 OpenMX 教程（安装、Si/金刚石能带与 DOS、收敛测试、DFT+U、SOC、轨道投影），4.0 手册第 81 节也收录了它。
- [hirtatsu/Calculations](https://github.com/hirtatsu/Calculations)：日文实操笔记，输入文件逐段讲解、tmux/nohup 运行与日志检查、超算 PBS 脚本。
- [arifiaapn/OpenMX-tutorial-bahasa-Indonesia-](https://github.com/arifiaapn/OpenMX-tutorial-bahasa-Indonesia-)：印尼语教程，按分支覆盖优化、NEB、能带、磁性、SOC、Berry 曲率、Z2 等主题，可作为功能覆盖面的对照清单。
- [hikaruri/DFT_startup](https://github.com/hikaruri/DFT_startup)：研究室入门资料，含 macOS 编译与 gnuplot 三维能带。

### 9.3 本仓库相关页面

- [OpenMX 篇目录](../README.md)、OpenMX 脚本 omx_make_input.py、[相关 GitHub 仓库汇总](../../Resources/github-repos.md)
- 公共篇：第3章 结构准备、第4章 ASE 与 pymatgen、第7章 DFT 理论速览与学习路线

---

## 参数选择与收敛建议（入门默认值速查）

下面是一组“第一次算新体系”时的保守起点，**不能代替收敛测试**，详细依据见第4、5章。

| 参数 | 入门建议值 | 说明 |
|---|---|---|
| 基组 | 数据库推荐表的 **Standard** 列（如 `H6.0-s2p1`、`C6.0-s2p2d1`、`O6.0-s2p2d1`、`Si7.0-s2p2d1`） | 先用 Quick 试跑，最终结果用 Standard 或 Precise |
| VPS | `*_PBE19`（GGA-PBE）或 `*_CA19`（LDA），与 `scf.XcType` 一致 | PAO 与 VPS 都来自 DFT_DATA19，不要混用旧库 |
| `scf.energycutoff` | 200–300 Ry | 显式写出；弱相互作用或小能量差问题用 300–400 Ry |
| `scf.criterion` | 1.0e-8 Hartree（单点能/性质），1.0e-9 或更严（结构优化、声子有限位移） | 默认 1.0e-6 偏松 |
| `scf.maxIter` | 100–300 | 默认 40 对磁性/金属体系往往不够 |
| `scf.Mixing.Type` | `RMM-DIISK`（周期体系）、`RMM-DIIS`（小分子） | 手册称 RMM-DIISK 对多数体系最稳健 |
| `scf.Mixing.History` / `scf.Mixing.StartPulay` | 30 / 10–30 | 手册推荐范围 |
| `scf.ElectronicTemperature` | 300 K（绝缘体/分子），金属可适当提高 | 过低易数值不稳定 |
| `scf.Kgrid` | 使 $N_i\cdot a_i \gtrsim 30$–40 Å | 再做 k 点收敛 |
| `scf.Generation.Kpoint` | 4.0 默认 `GCenter` 即可；复现 3.9 用 `Regular` | 论文中注明 |

---

## 常见问题与报错

| 问题 / 报错 | 原因 | 解决办法 |
|---|---|---|
| 计算开始不久退出，提示找不到 `.pao` 或 `.vps` 文件 | `DATA.PATH` 默认是相对路径 `../DFT_DATA19`，不在 `work/` 下运行就找不到 | 写绝对路径，指向含 `PAO/` 和 `VPS/` 的目录；核对文件名大小写，不写扩展名 |
| `Found an undefined species name ...` 后跟 `Please check your input file` | 坐标块中的物种名在 `<Definition.of.Atomic.Species` 中没有定义；或文件含全角空格、BOM、UTF-16 编码 | 核对物种名；用 `file x.dat` 检查编码，`iconv -f UTF-16 -t UTF-8`、`dos2unix` 转换 |
| 输出文件没有写到预期目录 | 把 `System.CurrrentDirectory` 写成了 `System.CurrentDirectory`（官方拼写是三个 r） | 按官方拼写书写；拼错的关键词会被当作未知关键词忽略而使用默认值 |
| 同一输入在 4.0 与 3.9 下总能量不同 | 4.0 默认 `scf.Generation.Kpoint GCenter`，3.9 为规则网格 | 需要对比时显式写 `Regular` |
| 想用 HSE06、SCAN 等泛函 | 4.0 手册中的 `scf.XcType` 只有 LDA、LSDA-CA、LSDA-PW、GGA-PBE | 带隙问题可考虑 DFT+U，或换用支持杂化泛函的程序 |
| 想算声子谱 | OpenMX 没有 DFPT | 用有限位移法，如带 OpenMX 接口的 ALAMODE |
| 费米面文件数值比以前大/小了一个固定倍数 | patch 4.0.1 把 bxsf/frmsf 单位从 Hartree/Bohr 改成 eV/Å | 更新后处理脚本中的单位换算 |
| 致密金属用大基组时出现很低的伪本征值、能量异常 | 基组过完备 | 换较小的优化基组，或 `scf.ProExpn.VNA off` 并适当提高 `scf.energycutoff` |
| 不确定某个关键词是否存在、默认值是多少 | 不同手册版本、示例注释说法不一 | 以 4.0 手册第 8.2 节和相应专题章节为准，同时在 `work/*.dat` 中搜索实际用法 |

---

## 练习

1. 下载 OpenMX 4.0 源码包并解压（暂不编译），用 `ls`、`du -sh` 查看 `source`、`work`、`DFT_DATA19` 的大小；用 `ls DFT_DATA19/VPS | wc -l` 统计 VPS 文件数量。
2. 用 `grep valence.electron` 查出 `H_PBE19`、`O_PBE19`、`Fe_PBE19H`、`Fe_PBE19S`、`Ni_PBE19H` 的价电子数，并与手册 12.3 节推荐表核对。
3. 计算基组 `O6.0-s2p2d1`、`Fe5.5H-s3p2d1`、`Bi8.0-s3p2d2f1` 各包含多少个基函数。
4. 把 `MD.Opt.criterion` 的默认值 3×10⁻⁴ Hartree/Bohr 换算为 eV/Å，并与 VASP 中常用的 `EDIFFG = -0.01` 或 `-0.02` 比较。
5. 打开 [OpenMX Viewer](https://www.openmx-square.org/viewer/)，输入 Materials Project 编号 `mp-149`（Si）加载结构，导出 OpenMX 输入，与 `cif2omx` 生成的结果比较 `scf.Kgrid` 与基组。
6. 根据 §8.2 的表格，为你自己的课题写出一份“必读章节 + 预计时间”的学习计划。

---

## 参考资料

**官方文档**

- OpenMX 官网：https://www.openmx-square.org/
- 下载页（4.0 + patch 4.0.1、3.9 + patch 3.9.9、GPU、deb、Docker）：https://www.openmx-square.org/download.html
- What's new（4.0 新功能、patch 4.0.1、Viewer 更新）：https://www.openmx-square.org/whatsnew/whatsnew.html
- Patch 4.0.1 README：https://www.openmx-square.org/bugfixed/26May08/README.txt
- 4.0 手册目录：https://www.openmx-square.org/openmx_man4.0/index.html；第 2 节 About OpenMX、第 3 节 Related papers、第 8 节 Input File、第 9 节 Output files、第 10 节 OpenMX Viewer、第 78 节 Included input examples、第 79 节 Known issues、第 81–82 节 信息源与工具链接
- 关键词说明：https://www.openmx-square.org/openmx_man4.0/s8_2_keywords.html
- cif2omx：https://www.openmx-square.org/openmx_man4.0/s8_3_gen-dat-files.html
- 推荐基组表（Quick/Standard/Precise）：https://www.openmx-square.org/openmx_man4.0/s12_3_basis-db2019.html
- k 点生成（GCenter/Regular）：https://www.openmx-square.org/openmx_man4.0/s15_k-sampling.html
- 截断能收敛：https://www.openmx-square.org/openmx_man4.0/s14_1_convergence.html
- 数据库 Ver. 2019：https://www.openmx-square.org/vps_pao2019/
- OpenMX Viewer：https://www.openmx-square.org/viewer/
- 论坛：https://www.openmx-square.org/forum/patio.cgi

**相关 GitHub 仓库**（只链接，不转载其内容；更多见 [Resources/github-repos.md](../../Resources/github-repos.md)）

- [pranabdas/openmx](https://github.com/pranabdas/openmx)：最完整的英文 OpenMX 教程站源码，含 4.0 + patch 4.0.1 的编译脚本与 Python 后处理示例（无许可证，仅供参考阅读）。
- [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)：官方下载页链接的 4.0 GPU 版（CMake，CUDA/HIP/CPU 后端）。
- [cmsi/MateriAppsLive](https://github.com/cmsi/MateriAppsLive)：预装 OpenMX 的 Linux 环境，适合免编译体验。
- [mailhexu/TB2J](https://github.com/mailhexu/TB2J) 与 [mailhexu/TB2J-OpenMX](https://github.com/mailhexu/TB2J-OpenMX)：由 `.scfout` 计算交换参数、DMI 与磁振子。
- [mzjb/DeepH-pack](https://github.com/mzjb/DeepH-pack)、[QuantumLab-ZY/HamGNN](https://github.com/QuantumLab-ZY/HamGNN)：以 OpenMX 哈密顿量为数据的机器学习哈密顿量方法。
- [deepmodeling/dpdata](https://github.com/deepmodeling/dpdata)：读取 OpenMX MD 输出（`.dat` + `.md`）转换为 DeePMD 训练数据。
- [alamode-code/alamode](https://github.com/alamode-code/alamode)：带 OpenMX 接口的声子/非谐/热导程序。
- [issp-center-dev/cif2x](https://github.com/issp-center-dev/cif2x)：基于 pymatgen，用 CIF + 模板批量生成 OpenMX 输入。
- [Tomoki-YAMASHITA/CrySPY](https://github.com/Tomoki-YAMASHITA/CrySPY)：晶体结构预测，可调用 OpenMX 弛豫。
- [HopTB/HopTB.jl](https://github.com/HopTB/HopTB.jl)：由 OpenMX 哈密顿量构建紧束缚模型做 Berry 曲率、光学与非线性响应。
- [mfukudaQED/PQViewer_Band_DOS](https://github.com/mfukudaQED/PQViewer_Band_DOS)：单 HTML 的能带 + DOS/PDOS 浏览器查看工具。
- [wxia529/acoord](https://github.com/wxia529/acoord)：VS Code 中编辑与可视化 OpenMX `.dat` 结构。

**关键论文**（4.0 手册第 3 节建议引用的“General”文献）

- T. Ozaki, Variationally optimized atomic orbitals for large-scale electronic structures, *Phys. Rev. B* 67, 155108 (2003). DOI: [10.1103/PhysRevB.67.155108](https://doi.org/10.1103/PhysRevB.67.155108)
- T. Ozaki, H. Kino, Numerical atomic basis orbitals from H to Kr, *Phys. Rev. B* 69, 195113 (2004). DOI: [10.1103/PhysRevB.69.195113](https://doi.org/10.1103/PhysRevB.69.195113)
- T. Ozaki, H. Kino, Efficient projector expansion for the ab initio LCAO method, *Phys. Rev. B* 72, 045121 (2005). DOI: [10.1103/PhysRevB.72.045121](https://doi.org/10.1103/PhysRevB.72.045121)
- K. Lejaeghere *et al.*, Reproducibility in density functional theory calculations of solids, *Science* 351, aad3000 (2016). DOI: [10.1126/science.aad3000](https://doi.org/10.1126/science.aad3000)（DFT_DATA19 的 Δ 因子验证即基于此基准，71 种元素平均 1.774 meV/atom）
- I. Morrison, D. M. Bylander, L. Kleinman, Nonlocal Hermitian norm-conserving Vanderbilt pseudopotential, *Phys. Rev. B* 47, 6728 (1993). DOI: [10.1103/PhysRevB.47.6728](https://doi.org/10.1103/PhysRevB.47.6728)

---

[返回目录](../README.md) | [下一章 →](02-install.md)
