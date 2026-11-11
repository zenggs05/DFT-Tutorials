# 示例 16：半经验方法（GFN-xTB、DFTB）与 QM/MM

配套教程：第16章 QM/MM 与半经验方法（xTB/DFTB）

## 目的

用几个小体系演示 CP2K 中三类“比 DFT 便宜”的计算方式：

| 输入 | 体系 | 方法 | 演示内容 |
|---|---|---|---|
| `ethanol_gfn1_native.inp` | 孤立乙醇分子 | CP2K 原生 GFN1-xTB | `METHOD XTB`、`SCF_GUESS MOPAC`、OT，GEO_OPT |
| `ice_gfn2_tblite_kpoints.inp` | 六方冰（4 个 H2O/胞，质子有序模型） | tblite 的 GFN2-xTB + k 点 | `GFN_TYPE TBLITE`、`&TBLITE METHOD GFN2`、4×4×2 MP 网格 |
| `dftb_md.inp` | 125 个水分子，15.54 Å 立方盒子 | SCC-DFTB（CP2K 自带 `DFTB/scc` 参数） | NVT（CSVR）MD，周期性 Ewald |
| `water_qmmm_tip3p.inp` | 1 个 QM 水 + 124 个 MM 水 | PBE/DZVP（GPW）+ 输入内定义的 TIP3P 类力场 | `METHOD QMMM`、`E_COUPL GAUSS`、周期性 QM/MM、短程 NVT MD |
| `water_qmmm_amber.inp` | 同上 | PBE/DZVP + AMBER prmtop | `PARMTYPE AMBER`、`CONN_FILE_FORMAT AMBER`、`MM_KIND OW/HW` |

## 文件说明

| 文件 | 说明 |
|---|---|
| `ethanol.xyz` | 乙醇初始结构（取自 ASE 的 g2 分子数据集） |
| `ethanol_gfn1_native.inp` | 原生 GFN1-xTB 几何优化；适用于 2023.x–2026.2 |
| `ice_gfn2_tblite_kpoints.inp` | tblite GFN2-xTB + k 点；需要编译时启用 tblite。默认按 2026.2 语法（`GFN_TYPE TBLITE` + `ACCURACY`），2025.2/2026.1 用 `-E NEWTB=0`；结构优化建议 ≥ 2026.1 |
| `dftb_md.inp` | SCC-DFTB 水盒子 MD（2000 步 × 0.5 fs） |
| `water125.xyz` | 125 个水分子的理想化初始结构（5×5×5 格点 + 随机取向，密度约 0.997 g/cm³），由本仓库作者用脚本生成 |
| `water125.pdb` | 同一结构的 PDB 版本（残基名 `WAT`，原子名 `O/H1/H2`），供 tleap 使用 |
| `water_qmmm_tip3p.inp` | 自包含的 QM/MM MD 输入（MM 参数直接写在输入里，无需外部文件） |
| `tleap_water.in` | tleap 脚本：由 `water125.pdb` 生成 `water125.prmtop` 与 `water125.inpcrd` |
| `water_qmmm_amber.inp` | 读取 AMBER 拓扑的 QM/MM 单点能 + 力 |
| `job_slurm.sh` | Slurm 作业脚本：依次运行各输入；没有 tblite 时自动跳过 GFN2 示例，没有 prmtop 时跳过 AMBER 示例 |

## 运行步骤

```bash
cd CP2K/examples/16-qmmm-xtb
export CP2K_DATA_DIR=/path/to/cp2k/data          # xTB/DFTB 参数、dftd3.dat、UZH 基组都在这里

# 0) 确认可执行文件的功能
cp2k.psmp -v | grep -i cp2kflags                  # 看是否含 tblite
ls "$CP2K_DATA_DIR"/xTB1_parameters "$CP2K_DATA_DIR"/DFTB/scc/scc_parameter

# 1) 语法检查
for f in *.inp; do cp2k.psmp -c "$f" > /dev/null && echo "$f OK"; done

# 2) 半经验计算（几秒到几分钟）
cp2k.psmp -i ethanol_gfn1_native.inp -o ethanol_gfn1_native.out
cp2k.psmp -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out                  # 2026.2
# cp2k.psmp -E NEWTB=0 -i ice_gfn2_tblite_kpoints.inp -o ice_sp.out     # 2025.2 / 2026.1
cp2k.psmp -E RUN=GEO_OPT -i ice_gfn2_tblite_kpoints.inp -o ice_opt.out   # 可选：结构优化
mpirun -np 8 cp2k.psmp -i dftb_md.inp -o dftb_md.out

# 3) QM/MM
mpirun -np 8 cp2k.psmp -i water_qmmm_tip3p.inp -o water_qmmm_tip3p.out

# 4) AMBER 拓扑版本（需要 AmberTools）
tleap -f tleap_water.in                           # 生成 water125.prmtop / water125.inpcrd
mpirun -np 8 cp2k.psmp -i water_qmmm_amber.inp -o water_qmmm_amber.out

# 或者一次性提交
sbatch job_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep -c "PROGRAM ENDED AT" *.out                         # 每个输出应为 1
grep "xTB| Parameter file" ethanol_gfn1_native.out       # 原生 xTB 读取的参数文件
grep "GEOMETRY OPTIMIZATION COMPLETED" ethanol_gfn1_native.out
grep "BRILLOUIN| List of Kpoints" ice_sp.out             # k 点确实生效
grep "ENERGY| Total FORCE_EVAL" ice_sp.out dftb_md.out water_qmmm_amber.out | tail -5
tail -3 water125_dftb_nvt-1.ener                         # 列：步数 时间 Kin Temp Pot Cons Qty UsedTime
grep -i "QM/MM LINKS" water_qmmm_tip3p.out               # 本例没有切断化学键，应显示 NO QM/MM LINKS DETECTED
tail -3 water_qmmm_tip3p-1.ener
```

判断标准：

- 乙醇优化在几十步内收敛（`GEOMETRY OPTIMIZATION COMPLETED`），C–O、O–H 键长与初始结构相差不大；
- 冰的 GFN2 单点能正常结束，`ATOMIC FORCES` 中 H 原子受力最大（H 位置只是按规则放置的）；
- DFTB 与 QM/MM 的 MD：`-1.ener` 中温度在 300 K 附近波动；由于初始盒子是理想化格点，开始的几百步温度和势能会有明显弛豫，应当舍弃；
- QM/MM 输出开头列出 MM 点电荷（`MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL`），以及 `QM_KIND`/`MM_KIND` 信息。

本仓库没有附带参考输出；能量、步数与耗时取决于 CP2K 版本、编译选项和并行设置。

## 注意事项

- `GFN_TYPE` 自 2025.1 起才有（2025.1–2026.1 中是整数 0/1，选原生 GFN0/GFN1）；tblite 自 2025.2 起可用：2025.2/2026.1 中只写 `&TBLITE` 段即启用，2026.2 中必须写 `GFN_TYPE TBLITE`，且 `&TBLITE` 段单独出现会报 `The XTB/TBLITE section requires XTB/GFN_TYPE TBLITE`。2024.x 及更早只有原生 GFN1（不写 `GFN_TYPE`）；
- 周期性的原生 xTB/DFTB 在 ≤ 2026.1 中必须写 `DO_EWALD T`；2026.2 起由 `&CELL PERIODIC` 自动决定，`DO_EWALD` 被标为 deprecated；
- `water125.xyz` 不是平衡结构。正式的 QM/MM 研究应先用纯 MM（`METHOD FIST`）做能量最小化与 NVT/NPT 平衡，再切换到 QM/MM；
- `water_qmmm_tip3p.inp` 中 Lennard-Jones 的 `EPSILON` 用默认单位 K（0.1521 kcal/mol ≈ 76.54 K），键与角的力常数用原子单位，换算见输入文件头部注释；
- QM 盒子（`&QMMM &CELL`）必须能装下 QM 电子密度且小于 MM 盒子；本例 QM 只有一个水分子，10 Å 足够。QM 区变大后要同步增大，并重新检查 `CUTOFF`；
- `job_slurm.sh` 中的环境设置（`module load`、`source cp2k_env`）请按你的集群修改。
