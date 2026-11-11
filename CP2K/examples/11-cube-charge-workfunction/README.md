# 示例 11：cube 文件、电荷分析、molden 与功函数

配套教程：[第11章 实空间性质：cube 文件、电荷分析与功函数](../../tutorials/11-cube-charge-workfunction.md)

## 目的

| 输入 | 体系 | 演示内容 | 适用版本 |
|---|---|---|---|
| `h2o_density_mo_cubes.inp` | 孤立 H2O | `&E_DENSITY_CUBE`、`&V_HARTREE_CUBE`、`&MO_CUBES`（NHOMO/NLUMO）、`&ELF_CUBE`、Mulliken/Löwdin/Hirshfeld 电荷、偶极矩 | 2024.x–2026.2 |
| `h2o_molden_for_multiwfn.inp` | 孤立 H2O | `&MO_MOLDEN`（`NDIGITS 9`、`NLUMO`、`WRITE_PSEUDO`），供 Multiwfn 做波函数分析 | ≥ 2026.2（旧版本按注释删减） |
| `pt111_workfunction.inp` | Pt(111) p(2×2) 四层 slab | 偶极校正 + `V_HARTREE_CUBE`，平面平均求真空能级与功函数 | 2024.x–2026.2 |
| `co_on_cu_charge_diff/` | CO/Cu(111) | 三次单点 + `chgdiff.py` 求差分电荷密度，见 [co_on_cu_charge_diff/README.md](co_on_cu_charge_diff/README.md) | 2024.x–2026.2 |

## 文件说明

| 文件 | 说明 |
|---|---|
| `h2o_density_mo_cubes.inp` | H2O 的密度、Hartree 势、HOMO/LUMO、ELF cube 与电荷分析 |
| `h2o_molden_for_multiwfn.inp` | 输出 `h2o_molden-MOS-1_0.molden` |
| `pt111_workfunction.inp`、`pt111_2x2_4L.xyz` | Pt(111) slab（a = 3.924 Å 实验值，层间距 $a/\sqrt3$，c = 30 Å，γ = 60°） |
| `workfunction.py` | 读 `V_HARTREE` cube 和 CP2K 输出，做 xy 平面平均，给出两侧真空能级和功函数，写 `vh_planar.dat`（只需 numpy） |
| `co_on_cu_charge_diff/` | 差分电荷密度子示例 |
| `job_slurm.sh` | Slurm 作业脚本，参数为输入文件名 |

## 运行步骤

```bash
cd CP2K/examples/11-cube-charge-workfunction
export CP2K_DATA_DIR=/path/to/cp2k/data

# 1) H2O：cube 与电荷
mpirun -np 4 cp2k.psmp -i h2o_density_mo_cubes.inp -o h2o_density_mo_cubes.out
ls -1 h2o_cube-*.cube

# 2) H2O：molden（2026.2）
mpirun -np 4 cp2k.psmp -i h2o_molden_for_multiwfn.inp -o h2o_molden_for_multiwfn.out
ls -1 *.molden

# 3) Pt(111) 功函数（建议在集群上运行）
sbatch job_slurm.sh pt111_workfunction.inp
python workfunction.py pt111_wf-v_hartree-1_0.cube pt111_workfunction.out

# 4) 差分电荷密度
cd co_on_cu_charge_diff && bash run_all.sh
```

也可以用本仓库的 cube 工具做平面平均、差分和积分（用法见 `-h`）：

```bash
python ../../scripts/cube_tools.py -h
```

## 预期结果 / 如何检查

- **H2O cube**：生成 `h2o_cube-ELECTRON_DENSITY-1_0.cube`、`h2o_cube-v_hartree-1_0.cube`、`h2o_cube-ELF_S1-1_0.cube` 和 `h2o_cube-WFN_0000N_1-1_0.cube`（N = 3、4 为 HOMO−1、HOMO，5、6 为 LUMO、LUMO+1），名称以 `ls` 为准。输出中有 `HOMO - LUMO gap [eV]`、`Mulliken Population Analysis`、`Hirshfeld Charges`、Löwdin 布居和 `Dipole moment [Debye]` 各块：

  ```bash
  grep -A6 "Mulliken Population Analysis" h2o_density_mo_cubes.out
  grep -A6 "Hirshfeld Charges"           h2o_density_mo_cubes.out
  grep -A2 "Dipole moment \[Debye\]"     h2o_density_mo_cubes.out
  grep "HOMO - LUMO gap"                 h2o_density_mo_cubes.out
  ```

  O 的电荷为负、H 为正，三种方法的数值不同（定义不同，属正常）。PBE 下 H2O 偶极矩应接近实验值 1.85 D（气相），具体数值以计算为准。
- **molden**：`h2o_molden-MOS-1_0.molden` 中应有 `[Pseudo]` 段（2026.2，`WRITE_PSEUDO T`）和 4 个占据 + 10 个空轨道，且 OT 计算的轨道能量非零。用 Multiwfn 打开后可直接做 ADCH、CM5、AIM 等分析。
- **Pt(111)**：`workfunction.py` 打印两侧真空能级与功函数；对这个对称的干净 slab，两侧数值应基本一致，平台的 spread 应在 0.01 eV 量级。PBE 计算的 Pt(111) 功函数通常低于实验值几个十分之一 eV（实验值约 5.9–6.1 eV，不同测量差异较大，见 Derry 等 2015 年综述）。本例 4 层 slab 和 6×6×1 k 点只用于演示，具体数值以你的计算为准。

## 注意事项

- cube 的 `STRIDE` 默认是 2 2 2（每个方向隔一个网格点输出）。做 Bader 分析、精确积分或平面平均时用 `STRIDE 1 1 1`；只看图可以用 2 或更大；
- CP2K 的 `V_HARTREE` cube 采用“电子势能”的符号约定（官方说明：与物理静电势符号相反），因此可直接用 $\Phi=V_\mathrm{vac}-E_F$，不用再取负号；单位是 Hartree；
- 功函数计算必须有足够厚的真空层（≥ 15 Å，本例约 23 Å），不对称 slab 必须开 `SURFACE_DIPOLE_CORRECTION`；
- molden 文件**不支持 k 点**，需要用 Γ 点大超胞；周期体系在 2026.2 中加 `WRITE_CELL T`，旧版本要手工在 molden 文件开头插入 `[Cell]` 段（见 sobereva.com/651）；
- `MO_CUBES` 在 k 点计算中不输出（只给警告），轨道图要用 Γ 点超胞。
