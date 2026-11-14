# 示例 16：半经验方法（GFN-xTB、DFTB）与 QM/MM

配套教程：[第16章 QM/MM 与半经验方法（xTB/DFTB）](../../tutorials/16-qmmm-xtb.md)

## 目的

用几个小体系演示 CP2K 中三类比 DFT 便宜的计算方式：紧束缚半经验方法、DFT/MM，以及半经验 QM/MM。

| 输入 | 体系 | 方法 | 演示内容 | 教程小节 |
|---|---|---|---|---|
| `ethanol_gfn1_native.inp` | 孤立乙醇分子 | CP2K 原生 GFN1-xTB | `METHOD XTB`、`SCF_GUESS MOPAC`、OT，GEO_OPT | §4 |
| `ice_gfn2_tblite_kpoints.inp` | 六方冰（每胞 4 个 H2O，质子有序模型） | tblite 的 GFN2-xTB + k 点 | `GFN_TYPE TBLITE`、`&TBLITE METHOD GFN2`、`ACCURACY`、4×4×2 MP 网格 | §5 |
| `o2_triplet_spgfn2.inp` | 三重态 O2 | tblite GFN2 + `UKS`（即 spGFN2-xTB） | 开壳层 xTB、`MULTIPLICITY 3` | §5.4 |
| `dftb_md.inp` | 125 个水分子，15.54 Å 立方盒子 | SCC-DFTB（CP2K 自带 `DFTB/scc` 参数） | NVT（CSVR）MD，周期性 Ewald，UFF 色散 | §6 |
| `water_qmmm_tip3p.inp` | 1 个 QM 水 + 124 个 MM 水 | PBE/DZVP（GPW）+ 输入文件内定义的 TIP3P 类力场 | `METHOD QMMM`、`E_COUPL GAUSS`、周期性 QM/MM、短程 NVT MD | §8 |
| `water_qmmm_amber.inp` | 同上 | PBE/DZVP + AMBER prmtop | `PARMTYPE AMBER`、`CONN_FILE_FORMAT AMBER`、`MM_KIND OW/HW` | §9 |
| `xtb_qmmm_water.inp` | 同上 | GFN1-xTB/MM | `E_COUPL COULOMB`，最便宜的 QM/MM 流程测试 | §10 |
| `gromacs_qmmm_fragment.mdp` | （任意 GROMACS 体系） | GROMACS–CP2K 接口 | `qmmm-cp2k-*` mdp 选项 | §14 |

## 文件说明

| 文件 | 说明 |
|---|---|
| `ethanol.xyz` | 乙醇初始结构（取自 ASE 的 g2 分子数据集） |
| `ethanol_gfn1_native.inp` | 原生 GFN1-xTB 几何优化；适用于 2023.x–2026.2 |
| `ice_gfn2_tblite_kpoints.inp` | tblite GFN2-xTB + k 点；需要编译时启用 tblite。默认按 2026.2 语法（`GFN_TYPE TBLITE` + `ACCURACY`），2025.2/2026.1 用 `-E NEWTB=0`；结构优化建议用 2026.1 或更新版本 |
| `o2_triplet_spgfn2.inp` | 开壳层 spGFN2-xTB 几何优化；只按 2026.2 编写，需要 tblite |
| `dftb_md.inp` | SCC-DFTB 水盒子 MD（2000 步 × 0.5 fs） |
| `water125.xyz` | 125 个水分子的理想化初始结构（5×5×5 格点 + 随机取向，密度约 0.997 g/cm³），由本仓库作者用脚本生成 |
| `water125.pdb` | 同一结构的 PDB 版本（残基名 `WAT`，原子名 `O/H1/H2`），供 tleap 使用 |
| `water_qmmm_tip3p.inp` | 自包含的 DFT/MM MD 输入（MM 参数直接写在输入文件里，无需外部文件） |
| `tleap_water.in` | tleap 脚本：由 `water125.pdb` 生成 `water125.prmtop` 与 `water125.inpcrd` |
| `water_qmmm_amber.inp` | 读取 AMBER 拓扑的 DFT/MM 单点能 + 力 |
| `xtb_qmmm_water.inp` | GFN1-xTB/MM 的 NVT MD（2000 步），MM 力场与 `water_qmmm_tip3p.inp` 相同 |
| `gromacs_qmmm_fragment.mdp` | GROMACS mdp 片段，需追加到你自己体系的 mdp 中，不能单独运行 |
| `job_slurm.sh` | Slurm 作业脚本：依次运行各 CP2K 输入；没有 tblite 时自动跳过 tblite 示例，没有 prmtop 时跳过 AMBER 示例 |

## 运行步骤

```bash
cd CP2K/examples/16-qmmm-xtb
export CP2K_DATA_DIR=/path/to/cp2k/data          # xTB/DFTB 参数、dftd3.dat、UZH 基组都在这里

# 0) 确认可执行文件的功能与数据文件
cp2k.psmp -v | grep -i -o tblite                  # 无输出 = 没有 tblite
ls "$CP2K_DATA_DIR"/xTB1_parameters "$CP2K_DATA_DIR"/DFTB/scc/scc_parameter "$CP2K_DATA_DIR"/DFTB/uff_table

# 1) 语法检查（只解析输入，不计算）
for f in *.inp; do cp2k.psmp -c "$f" > /dev/null && echo "$f OK"; done

# 2) 半经验计算（几秒到几分钟）
cp2k.psmp -i ethanol_gfn1_native.inp -o ethanol_gfn1_native.out
cp2k.psmp -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out                  # 2026.2
# cp2k.psmp -E NEWTB=0 -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out     # 2025.2 / 2026.1
cp2k.psmp -E RUN=GEO_OPT -i ice_gfn2_tblite_kpoints.inp -o ice_opt.out   # 可选：结构优化
cp2k.psmp -i o2_triplet_spgfn2.inp -o o2_triplet_spgfn2.out             # 2026.2 + tblite
mpirun -np 8 cp2k.psmp -i dftb_md.inp -o dftb_md.out

# 3) QM/MM
cp2k.psmp -i xtb_qmmm_water.inp -o xtb_qmmm_water.out                   # GFN1-xTB/MM，很快
mpirun -np 8 cp2k.psmp -i water_qmmm_tip3p.inp -o water_qmmm_tip3p.out  # PBE/MM

# 4) AMBER 拓扑版本（需要 AmberTools）
tleap -f tleap_water.in                           # 生成 water125.prmtop / water125.inpcrd
mpirun -np 8 cp2k.psmp -i water_qmmm_amber.inp -o water_qmmm_amber.out

# 或者一次性提交
sbatch job_slurm.sh
```

`cp2k.psmp -o` 会把输出追加到已有文件末尾，重复运行前请先删除旧的 `.out`（`job_slurm.sh` 已经这样做）。

## 预期结果 / 如何检查

```bash
grep -c "PROGRAM ENDED AT" *.out                          # 每个输出应为 1
grep "GEOMETRY OPTIMIZATION COMPLETED" ethanol_gfn1_native.out o2_triplet_spgfn2.out
grep "BRILLOUIN|" ice_sp.out | head                       # k 点设置确实生效
grep "ENERGY| Total FORCE_EVAL" ice_sp.out water_qmmm_amber.out | tail -4
tail -3 water125_dftb_nvt-1.ener                          # 列：步数 时间 Kin Temp Pot ConsQty UsedTime
grep "Number of QM atoms" water_qmmm_tip3p.out xtb_qmmm_water.out    # 应为 3
grep "NO QM/MM LINKS DETECTED" water_qmmm_tip3p.out       # 本例没有切断化学键
grep "TOTAL CHARGE:" water_qmmm_tip3p.out                 # MM 点电荷之和，应为 0
tail -3 water_qmmm_tip3p-1.ener water_qmmm_xtb-1.ener
```

判断标准：

- 乙醇优化在几十步内收敛，C–O、C–C、O–H 键长约 1.4、1.5、0.97 Å；
- 冰的 GFN2 单点能正常结束，`ATOMIC FORCES` 中 H 原子受力最大（H 位置只是按规则放置的）；
- O2 优化收敛，输出中多重度为 3；
- DFTB 与 QM/MM 的 MD：`-1.ener` 中温度在 300 K 附近波动；由于初始盒子是理想化格点，开始的几百步温度和势能会有明显弛豫，应当舍弃；
- 静电嵌入的 QM/MM 输出开头列出 MM 点电荷（`MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL`）。

本仓库没有附带参考输出；能量、步数与耗时取决于 CP2K 版本、编译选项和并行设置。

## 注意事项

- `GFN_TYPE` 自 2025.1 起才有（2025.1–2026.1 中是整数 0/1，选原生 GFN0/GFN1）；tblite 自 2025.2 起可用：2025.2/2026.1 中只写 `&TBLITE` 段即启用；2026.2 中必须同时写 `GFN_TYPE TBLITE` 和 `&TBLITE` 段，否则报 `The XTB/TBLITE section requires XTB/GFN_TYPE TBLITE` 或 `XTB/GFN_TYPE TBLITE requires an XTB/TBLITE section`。2024.x 及更早只有原生 GFN1（不写 `GFN_TYPE`）；
- 2026.2 中 tblite 默认使用自己的 SCC 混合器（`SCC_MIXER AUTO`），其收敛由 `&TBLITE ACCURACY` 控制；
- 周期性的原生 xTB/DFTB 在 2026.1 及以前必须写 `DO_EWALD T`；2026.2 起由 `&CELL PERIODIC` 自动决定，`DO_EWALD` 被标为 deprecated；
- `E_COUPL` 默认是 `NONE`（机械嵌入）。GPW 用 `GAUSS`，xTB/DFTB 用 `COULOMB`；`GAUSS`/`S-WAVE` 耦合要求 `&DFT &MGRID` 中写 `COMMENSURATE`，否则程序中止（两个 GPW 示例都已写上）；
- AmberTools 默认的 TIP3P 拓扑是刚性水（两根 O–H 键加一根 H–H“键”），与 `water_qmmm_tip3p.inp` 中的柔性水不同。`water_qmmm_amber.inp` 只做单点能和力；要用 prmtop 做 MD，请在 tleap 中加 `set default FlexibleWater on`，或在 CP2K 中加约束；
- `water125.xyz` 不是平衡结构。正式的 QM/MM 研究应先用纯 MM（`METHOD FIST` 或 AMBER/GROMACS）做能量最小化与 NVT/NPT 平衡，再切换到 QM/MM；
- `water_qmmm_tip3p.inp` 中 Lennard-Jones 的 `EPSILON` 用单位 K（0.1521 kcal/mol ≈ 76.54 K），键与角的力常数用原子单位，换算见输入文件头部注释；
- AMBER 的 TIP3P 水 H 没有 LJ 参数，长时间的静电嵌入 QM/MM MD 中可能导致 QM 区被过度极化。正式计算请参考教程 §9 用 ParmEd 补参数；
- QM 盒子（`&QMMM &CELL`）必须能装下 QM 电子密度且小于 MM 盒子；本例 QM 只有一个水分子，10 Å 足够。QM 区变大后要同步增大，并重新检查 `CUTOFF`；
- `job_slurm.sh` 中的环境设置（`module load`、`source cp2k_env`）请按你的集群修改。

返回教程：[第16章](../../tutorials/16-qmmm-xtb.md)
