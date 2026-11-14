# 示例 12-1：Si 的声子谱、声子 DOS 与热力学性质（phonopy 有限位移）

## 目的

- 完成 phonopy 有限位移法的标准流程：紧致优化 → 生成位移超胞 → VASP 单点力 → FORCE_SETS → 声子谱/DOS/热力学；
- 同时覆盖 phonopy ≥ 4.0（`phonopy-init`）与 2.x/3.x 的写法；
- 熟悉 band.conf、mesh.conf、pdos.conf 的关键标签。

对应教程：[第12章 声子：phonopy 有限位移与 DFPT](../../tutorials/12-phonon.md)（§4–§7）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 金刚石 Si 惯用胞（8 原子），a = 5.469 Å（PBE 近似值）；**生成超胞前请换成你自己优化后的 CONTCAR** |
| `INCAR_relax`、`KPOINTS_relax` | 紧致优化：`ISIF = 3`、`EDIFF = 1E-8`、`EDIFFG = -1E-3`；Γ 中心 8×8×8 |
| `INCAR`、`KPOINTS` | 64 原子超胞单点力：`IBRION = -1`、`NSW = 0`、`EDIFF = 1E-8`、`PREC = Accurate`、`LREAL = .FALSE.`；Γ 中心 3×3×3 |
| `band.conf` | Γ–X–W–K–Γ–L 声子谱；`DIM = 2 2 2`，显式 fcc `PRIMITIVE_AXES` |
| `mesh.conf` | 24×24×24 Γ 中心网格，`TPROP = .TRUE.`（0–1000 K） |
| `pdos.conf` | 原子投影声子 DOS（`PDOS = 1, 2`） |
| `run_disp.sh` | `prepare` / `submit` / `check` / `collect` 四个子命令，自动识别 phonopy 版本 |
| `job.slurm` | 单个计算目录的 Slurm 脚本 |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/12-si-phonopy-fd
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

# 1) 紧致优化（可重复一次直到体积不变）
mkdir -p relax
cp INCAR_relax relax/INCAR; cp KPOINTS_relax relax/KPOINTS; cp POSCAR POTCAR relax/
(cd relax && sbatch ../job.slurm)
cp relax/CONTCAR POSCAR

# 2) 生成位移超胞并建立 disp-XXX 目录
bash run_disp.sh prepare          # 金刚石 Si 只需 1 个位移 (disp-001)

# 3) 提交、检查
bash run_disp.sh submit
bash run_disp.sh check

# 4) 收集力 -> FORCE_SETS
bash run_disp.sh collect

# 5) 后处理
phonopy --config band.conf -p -s         # phonopy >= 4.x
phonopy --config mesh.conf -p -s
phonopy --config pdos.conf -p -s
# phonopy 2.x/3.x:  phonopy -p -s band.conf   等
phonopy-bandplot --gnuplot band.yaml > band.dat
```

## 预期结果 / 如何检查

- `run_disp.sh prepare` 输出 `created 1 displacement directories`；`SPOSCAR` 为 64 原子；
- `disp-001/OSZICAR` 只有 1 个离子步，SCF 未达到 NELM；
- 后处理日志中 `Max drift of force constants` 很小（远小于 0.01）；
- 声子谱无虚频，Γ 点 3 支声学模为 0，光学模在 Γ 点约 15 THz 量级（PBE 近似；实验 Raman 峰约 520 cm⁻¹ ≈ 15.6 THz）；
- `thermal_properties.yaml` 中高温 C_V 趋于 Dulong–Petit 极限 3NR ≈ 49.9 J/K/mol（每 2 原子原胞）。

## 注意事项

- phonopy 4.10 中 `phonopy band.conf`（配置文件作位置参数）会报错，必须用 `--config`；`-d`、`-f`、`--fc` 已移到 `phonopy-init`。
- conf 中显式写了 `PRIMITIVE_AXES`，避免 v4 默认 `AUTO` 与 v3 默认单位矩阵造成结果不一致。
- 单点计算**不能**弛豫（`IBRION = -1`、`NSW = 0`）。
- 收敛检查：把 `DIM` 改为 3 3 3（`DIM="3 3 3" bash run_disp.sh prepare`，KPOINTS 改为 2×2×2）比较声子谱。
- 位移超胞较多时可改用 Slurm job array 提交（见 [Common 第2章](../../../Common/tutorials/02-job-scheduler.md)）。
- 也可以使用仓库的自动化脚本 [`phonopy_workflow.sh`](../../scripts/phonopy_workflow.sh)（用法见 `bash ../../scripts/phonopy_workflow.sh -h`）。
