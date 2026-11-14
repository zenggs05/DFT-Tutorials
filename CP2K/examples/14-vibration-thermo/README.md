# 示例 14：振动分析、红外强度与热化学（H2O、CO/Pt(111)、过渡态验证）

配套教程：[第14章 振动分析、红外与热化学](../../tutorials/14-vibration-thermo.md)

## 目的

用三个体系演示 CP2K `RUN_TYPE VIBRATIONAL_ANALYSIS` 的典型用法：

| 体系 | 输入 | 演示内容 |
|---|---|---|
| 孤立 H2O | `h2o_opt_tight.inp` → `h2o_vib_ir.inp` | 紧优化 → 全 Hessian 频率、IR 强度（`&MOMENTS`）、RRHO 热化学（298.15 K） |
| CO 顶位吸附在 Pt(111) 3×3×3 slab 上 | `co_pt111_opt.inp` → `co_on_pt_partial_hessian.inp` | 固定全部 Pt 的部分 Hessian，吸附物 6 个模式，简谐极限的 ZPE/振动熵 |
| [第13章](../../tutorials/13-neb-dimer.md) H/Cu(111) 的过渡态 | `ts_freq_check.inp` | 与 NEB/Dimer 同一 `FORCE_EVAL`，检查只有一个虚频 |

所有输入用 `@INCLUDE` 共享 `&FORCE_EVAL`，保证优化和频率计算的设置完全一致。

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o.xyz` | H2O 初始结构（实验气相几何） |
| `h2o_force_eval.inc` | H2O 共享 `&FORCE_EVAL`：PBE/TZV2P-MOLOPT-PBE-GTH，600 Ry，`EPS_SCF 1E-8`，MT 泊松求解器，`&MOMENTS PERIODIC F` |
| `h2o_opt_tight.inp` | 步骤 1：`MAX_FORCE 1.0E-5` 的紧优化 |
| `h2o_vib_ir.inp` | 步骤 2：频率 + IR + 热化学，写出 `h2o_vib-VIBRATIONS-1.mol` |
| `co_pt111.xyz` | CO/Pt(111) 初始结构（未优化）：1–9 底层、10–18 中层、19–27 表层 Pt，28 = C，29 = O；CO 在 23 号 Pt 正上方 |
| `co_pt111_force_eval.inc` | CO/Pt 共享 `&FORCE_EVAL`：PBE+D3(BJ)，DZVP，500 Ry，对角化 + 展宽，Γ 点 |
| `co_pt111_opt.inp` | 步骤 3：固定 1–18 的优化，`MAX_FORCE 1.0E-4` |
| `co_on_pt_partial_hessian.inp` | 步骤 4：固定 1–27，只位移 C、O |
| `ts_freq_check.inp` | 步骤 5：TS 验证，包含 `../13-neb-dimer-h-cu111/cu111_h_force_eval.inc`，默认可动原子为表层 Cu + H（`-E NFIX=27` 只放开 H） |
| `last_frame.sh` | 取 xyz 轨迹最后一帧 |
| `harmonic_thermo.py` | 从输出或 molden 文件读频率，计算简谐极限下的 ZPE、U_vib、S_vib、F_vib（标准库） |
| `displace_mode.py` | 沿 molden 中的某个模式正负扰动结构（处理虚频、“伪 IRC”） |
| `job_slurm.sh` | Slurm 脚本，用 `STEP=` 选择步骤，并自动设置 `NPROC_REP` |

晶胞：H2O 为 12 Å 立方盒子（`PERIODIC NONE`）；CO/Pt 为 `ABC 8.42164 8.42164 24.0`、`ALPHA_BETA_GAMMA 90 90 60`。

## 运行步骤

```bash
cd CP2K/examples/14-vibration-thermo
export CP2K_DATA_DIR=/path/to/cp2k/data     # 若已内置可省略

# H2O
mpirun -np 8  cp2k.psmp -i h2o_opt_tight.inp -o h2o_opt_tight.out
bash last_frame.sh h2o_opt-pos-1.xyz h2o_opt.xyz
mpirun -np 18 cp2k.psmp -i h2o_vib_ir.inp -o h2o_vib_ir.out          # 9 replicas x NPROC_REP 2

# CO/Pt(111)
mpirun -np 64 cp2k.psmp -i co_pt111_opt.inp -o co_pt111_opt.out
bash last_frame.sh co_pt111_opt-pos-1.xyz co_pt111_opt.xyz
mpirun -np 96 cp2k.psmp -i co_on_pt_partial_hessian.inp -o co_on_pt_partial_hessian.out
python3 harmonic_thermo.py co_on_pt_partial_hessian.out
python3 harmonic_thermo.py co_pt111_vib-VIBRATIONS-1.mol -T 500 --floor 50

# 过渡态验证（需要先完成第13章示例，得到 ts_dimer.xyz）
cp ../13-neb-dimer-h-cu111/ts_dimer.xyz .
mpirun -np 64 cp2k.psmp -i ts_freq_check.inp -o ts_freq_check.out
mpirun -np 48 cp2k.psmp -E NFIX=27 -E NPROC_REP=8 -i ts_freq_check.inp -o ts_freq_check_H.out   # 只放开 H

# 或者用 Slurm
sbatch --export=ALL,STEP=h2o_opt job_slurm.sh
```

## 预期结果 / 如何检查

```bash
grep "GEOMETRY OPTIMIZATION COMPLETED" h2o_opt_tight.out co_pt111_opt.out
grep -A4 "Minimum Structure - Energy and Forces" h2o_vib_ir.out   # 残余力应很小
grep "VIB|Frequency" h2o_vib_ir.out co_on_pt_partial_hessian.out ts_freq_check.out
grep "VIB|IR int"    h2o_vib_ir.out
grep "No IR intensities" h2o_vib_ir.out                          # 应为空
grep -A12 "THERMOCHEMICAL DATA" h2o_vib_ir.out
ls *-VIBRATIONS-1.mol
```

- H2O：3 个实频，弯曲约 1600 cm⁻¹，两个 O–H 伸缩在 3600–3800 cm⁻¹ 区间；没有虚频。实验基频为 1595、3657、3756 cm⁻¹（NIST），计算值是 PBE 简谐频率，与实验差几十到一百多 cm⁻¹ 属正常；
- H2O 的热化学块中 `Symmetry number` 为 1：要得到正确的 G，需加 $RT\ln 2$ ≈ 1.72 kJ/mol（298.15 K），见教程 §4.5；
- CO/Pt：6 个模式，最高的是 C–O 伸缩（2000 cm⁻¹ 上下，低于气相 CO 的 2143 cm⁻¹），其余为 Pt–C 伸缩与受阻转动/平动；输出中**没有**热化学块（有固定原子时 CP2K 不计算），用 `harmonic_thermo.py`；
- 过渡态：`ts_freq_check.out` 中恰好一个负频率，molden 动画中 H 沿 fcc ↔ hcp 方向运动。

本仓库没有附带参考输出；具体数值取决于 CP2K 版本、编译方式和并行设置。

## 注意事项

- `TC_TEMPERATURE` 的默认值是 273.15 K，本例已显式设为 298.15 K；
- 未质量加权的 Hessian 输出（`VIB| Hessian (before multiplying by ...)`）需要 CP2K ≥ 2026.2；`Frc consts` 一行在 2026.1 之前的版本中有误，频率本身不受影响；
- `@INCLUDE` 和 `COORD_FILE_NAME` 的相对路径都相对于运行 CP2K 的当前目录，请在本目录下运行；
- Pt slab 只算 Γ 点，属教学设置；正式计算请加 `&KPOINTS`（如 `MONKHORST-PACK 3 3 1`）并检查层数和真空；
- 想用 D2O 做同位素对比时，取消 `h2o_force_eval.inc` 中 `MASS 2.014102` 一行的注释；
- 不要在频率计算中使用 `IGNORE_CONVERGENCE_FAILURE`：未收敛的力会直接污染 Hessian。
