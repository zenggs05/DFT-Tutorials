# VASP 常见报错对照表（错误信息 → 原因 → 修复）

对应教程：[第18章 后处理生态、自动化与报错排查](../../tutorials/18-postprocess-troubleshooting.md)（§9）。错误字符串取自 VASP 的标准输出/OUTCAR，以及 custodian `VaspErrorHandler` 所监测的关键词；修复办法是作者整理的经验做法，请结合具体体系判断。

先用下面的命令定位问题：

```bash
tail -30 vasp.out                        # 标准输出（作业脚本重定向的文件）
grep -n -i -E "error|warning|very bad news|advice" vasp.out OUTCAR | head -40
tail -5 OSZICAR                          # 电子步/离子步是否在推进
grep -c "General timing" OUTCAR          # 1 = 正常结束；0 = 中途退出或被杀
python ../../scripts/check_vasp.py -h    # 仓库的健康检查脚本（自动匹配常见报错）
```

## A. 输入文件

| 错误信息（关键片段） | 常见原因 | 修复 |
|---|---|---|
| `Error reading item 'XXX' from file INCAR` | 标签值格式错误；INCAR 有 Tab 或 Windows 换行（CRLF）；MAGMOM 个数不对 | `dos2unix INCAR`；用空格代替 Tab；检查拼写；MAGMOM 用 `N*value` 写满 NIONS 个（非共线 3×NIONS） |
| 标签被静默忽略（OUTCAR 回显中找不到） | 拼写错误、行首有不可见字符、标签写在注释后面 | `grep -i 标签名 OUTCAR` 确认 VASP 读到的值 |
| 元素与 POTCAR 不匹配，但不报错，结果离谱 | POTCAR 拼接顺序与 POSCAR 第 6 行不一致 | `grep TITEL POTCAR` 与 POSCAR 元素行逐一对照；用脚本按 POSCAR 生成 POTCAR |
| `ERROR: the triple product of the basis vectors is negative` | 晶格矢量构成左手系 | 交换任意两个晶格矢量（同时交换坐标分量），或用 ASE/pymatgen 重写 POSCAR |
| `The distance between some ions is very small` | 原子重叠，常见于建模错误或坐标单位（Direct/Cartesian）写反 | 用 VESTA 检查结构；确认第 8 行坐标类型 |
| `ERROR: charge density could not be read from file CHGCAR for ICHARG>10` | 非自洽计算缺少 CHGCAR，或 CHGCAR 损坏 | 从自洽目录复制完整的 CHGCAR |
| `WARNING: dimensions on CHGCAR file are different` | 读入的 CHGCAR 网格与当前设置不同（ENCUT、PREC、晶胞改变） | 保持两步计算设置一致，或删除旧 CHGCAR 重新自洽 |
| `while reading WAVECAR` / `while reading plane` | WAVECAR 与当前 ENCUT、NBANDS、k 点或晶胞不匹配，或文件损坏 | 删除 WAVECAR 从头开始；或保持设置与生成 WAVECAR 时一致 |

## B. 对称性与 k 点

| 错误信息 | 常见原因 | 修复 |
|---|---|---|
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角/三方体系用了偶数 Monkhorst–Pack 网格，k 网格与晶格对称性不兼容 | 改为 Γ 中心网格；或调整网格使之与晶格对称性相容 |
| `Tetrahedron method fails (number of k-points < 4)`、`Routine TETIRR needs special values` | `ISMEAR = -5` 但 k 点太少或不是 Γ 中心 | 增加 k 点并用 Γ 中心网格；分子、slab、Γ 点计算改用 `ISMEAR = 0` |
| `internal error in subroutine PRICEL` | 结构对称性有数值噪声 | 用 spglib/vaspkit 602/603 整理结构；减小 `SYMPREC`（如 1E-6）或 `ISYM = 0` |
| `Found some non-integer element in rotation matrix`、`SGRCON` | 同上，或 Monkhorst–Pack 网格破坏对称性 | 改 Γ 中心网格；整理结构；必要时 `ISYM = 0` |
| `... rotation matrix was not found (increase SYMPREC)` | 对称性判定容差不合适 | 先改 Γ 中心网格，再适度增大 `SYMPREC`（不超过 1E-4），最后 `ISYM = 0` |
| `POSMAP`、`RHOSYG` 相关内部错误 | 对称性映射失败 | 调整 `SYMPREC`；`ISYM = 0` |
| `PRICELV: current lattice and primitive lattice are incommensurate` | 弛豫过程中晶格对称性改变 | 调整 `SYMPREC`；用 CONTCAR 重新开始；`ISYM = 0` |
| `Inconsistent Bravais lattice types found for crystalline and reciprocal lattice` | 对称性容差问题 | `SYMPREC = 1E-6` 或更小；`ISYM = 0` |

## C. 电子自洽

| 错误信息 | 常见原因 | 修复 |
|---|---|---|
| `WARNING: Sub-Space-Matrix is not hermitian in DAV` | Davidson 对角化不稳定：原子重叠、步长过大、`LREAL` 误差、并行过度 | 检查结构；`LREAL = .FALSE.`；`PREC = Accurate`；`ALGO = Normal → All`；减少核数 |
| `Error EDDDAV: Call to ZHEGV failed` | 同上，矩阵不正定 | `ALGO = All`（6.5.1 起可加 `ISEARCH = 1`）；删除 CHGCAR/WAVECAR 从头开始；减小 POTIM |
| `WARNING in EDDRMM: call to ZHEGV failed` | RMM-DIIS（`ALGO = Fast/VeryFast`）不稳定 | `ALGO = Normal`；减小 POTIM；从头开始 |
| `ERROR in subspace rotation PSSYEVX` / `PDSYEVX` | ScaLAPACK 子空间旋转失败 | `ALGO = Normal`；调整并行（NCORE）；必要时关闭 ScaLAPACK（`LSCALAPACK = .FALSE.`） |
| `BRMIX: very serious problems the old and the new charge density differ` | 读入了不匹配的 CHGCAR；对称性判定不一致；混合不稳定 | 删除 CHGCAR/WAVECAR；`ISYM = 0` 或调整 `SYMPREC`；改 Γ 中心网格；尝试 `IMIX = 1` |
| `EDWAV: internal error, the gradient is not orthogonal` | 编译优化或数值问题、`ALGO` 组合不当 | 换 `ALGO = Normal`；检查编译（降低优化等级）；减少核数 |
| 达到 NELM 仍未收敛（OSZICAR 中 dE 振荡） | 金属/磁性体系混合困难；初始磁矩不合理；`SIGMA` 太小 | 增大 `NELM`；`ALGO = All`；调 `AMIX/BMIX`（或 `AMIX_MAG/BMIX_MAG`）；合理设置 MAGMOM；DFT+U 加 `LMAXMIX` |
| `One of the lattice vectors is very long (>50 A), but AMIN is rather large` | 长轴 slab/真空方向很长 | `AMIN = 0.01` |
| `WARNING: Your highest band is occupied at some k-points!` | NBANDS 太少（金属或高温 smearing） | 增大 `NBANDS`（约 1.2–1.5 倍） |
| `The number of bands has been changed from X to Y` | 为适配并行，VASP 自动增加了 NBANDS | 一般无害；LOBSTER/LOPTICS/续算 MLFF 时手动设 NBANDS 为能带组数（ranks/KPAR/NCORE）的整数倍 |
| `ERROR RSPHER`、`REAL_OPTLAY: internal error` | 实空间投影（`LREAL`）在某些体系/网格下失败 | `LREAL = .FALSE.` |
| `Spin polarized Harris functional dynamics is a good joke` | `ICHARG ≥ 10`、`ISPIN = 2` 且 `NSW > 0` | 非自洽计算设 `NSW = 0` |

## D. 离子弛豫与动力学

| 错误信息 | 常见原因 | 修复 |
|---|---|---|
| `ZBRENT: fatal error in bracketing` | 共轭梯度线搜索失败：离极小点很远、受力噪声大（SCF 精度不足） | `cp CONTCAR POSCAR` 续算；`EDIFF` 收紧到 1E-6、`NELMIN = 6–8`；减小 POTIM；改 `IBRION = 1` 或 VTST FIRE（`IBRION = 3`、`POTIM = 0`、`IOPT = 7`） |
| `BRIONS problems: POTIM should be increased` | `IBRION = 1` 步长不合适 | 略增 POTIM，或改 `IBRION = 2` |
| `LAPACK: Routine ZPOTRF failed` | 原子过近（第一步）或步长过大（后续步）；小体系 NCORE 过大 | 检查结构；减小 POTIM；`ISYM = 0`；小体系减小 NCORE |
| 优化后体积偏小/偏大 | Pulay 应力，ENCUT 太低 | `ENCUT ≥ 1.3×ENMAX`；用 CONTCAR 重新优化直到体积不变 |
| 结构"乱飞"、能量突然升高 | POTIM 太大；MD 中 `POTIM` 对含 H 体系过大 | 减小 POTIM；MD 含 H 时用 0.5–1 fs |
| NEB 图像间距异常、不收敛 | 初末态原子顺序不一致；进程数不能被 IMAGES 整除；过早打开 CI | `dist.pl` 检查；进程数取 IMAGES 的整数倍；先普通 NEB 再 `LCLIMB = .TRUE.` |

## E. 特定功能

| 错误信息 | 常见原因 | 修复 |
|---|---|---|
| `ELF: KPAR>1 not implemented` | ELF 不支持 k 点并行 | `KPAR = 1`（并设 `NPAR = 1`） |
| `WARNING: ELF not implemented for non collinear case` | 非共线计算 | 用共线计算求 ELF |
| `PEAD routines do not work for NCORE` / `remove the tag NPAR from the INCAR file` | `LEPSILON`/`LCALCEPS`/DFPT 不支持能带并行 | 删除 `NCORE`、`NPAR` |
| DFPT（`IBRION = 7/8`）中 vdW 修正无效 | DFPT 不支持 `IVDW` | 改用有限位移（`IBRION = 5/6` 或 phonopy） |
| SOC 用 `vasp_std` 报错或 MAGMOM 个数错误 | 非共线需要 `vasp_ncl`，MAGMOM 每原子 3 个分量；6.5.1 之前同时设 ISPIN=2 与 LNONCOLLINEAR 会报 MAGMOM 个数错误 | 用 `vasp_ncl`；MAGMOM 写 3×NIONS 个数 |
| MLFF 续训报找不到 ML_AB / run 模式找不到 ML_FF | 没有把 ML_ABN 拷成 ML_AB，或没有 ML_FF | 续训前 `cp ML_ABN ML_AB`；`cp ML_FFN ML_FF` 后 `ML_MODE = run` |
| 6.6.0 自旋极化 meta-GGA（SCAN/r2SCAN）受力错误 | Known issue 94 | 升级到 6.6.1，或给 6.6.0 打官方 `pot.F.patch` 后重新编译 |

## F. 资源与环境

| 现象 | 常见原因 | 修复 |
|---|---|---|
| `forrtl: severe (174): SIGSEGV, segmentation fault occurred` | 栈空间不足；内存不足；编译/库不匹配 | 作业脚本中 `ulimit -s unlimited`；减少每节点进程数；检查 MPI 与编译器版本 |
| 作业被杀（OOM、`Killed`、`exit code 137`） | 内存不足 | 增大 `NCORE`、减小 `KPAR`；每节点少放进程；减小 NBANDS/ENCUT 或换节点 |
| `M_divide: can not subdivide N nodes by M` | `NCORE`/`NPAR`/`KPAR` 与进程数不整除 | 选用进程数的因子 |
| GPU 版本很慢 | 每块 GPU 多个 MPI 进程；设了 NCORE>1 | 每块 GPU 1 个 rank，用 OpenMP 线程；`KPAR` = GPU 数；适当增大 `NSIM` |
