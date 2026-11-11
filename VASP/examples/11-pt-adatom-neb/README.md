# 示例 11-1：Pt 吸附原子在 Pt(100) 上的跳跃扩散（VTST CI-NEB）

## 目的

- 走通一次完整的 VTST NEB 流程：端点优化 → 插值（IDPP/线性）→ 普通 NEB 粗收敛 → CI-NEB → nebresults.pl 能垒；
- 学会检查原子顺序、图像间距、MPI 进程数与 IMAGES 的关系；
- 为 Dimer（示例 11-3）和过渡态频率验证（示例 11-4）提供起点。

对应教程：[第11章 过渡态搜索：NEB / CI-NEB / Dimer](../../tutorials/11-neb-dimer.md)（§4–§8、§12）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR_IS` | Pt(100) p(3×3) 4 层（底部 2 层固定）+ 吸附原子在 hollow 位（第 37 号原子），a = 3.97 Å（PBE 近似值） |
| `POSCAR_FS` | 吸附原子位于相邻 hollow 位（沿 x 移动 2.81 Å），原子顺序与 IS 完全相同 |
| `INCAR_relax` | 端点优化（`IBRION = 2`、`EDIFFG = -0.02`） |
| `INCAR_neb_stage1` | 第一阶段：普通 NEB，`LCLIMB = .FALSE.`、`IOPT = 3`（quick-min）、`EDIFFG = -0.10` |
| `INCAR_neb_stage2` | 第二阶段：CI-NEB，`LCLIMB = .TRUE.`、`IOPT = 1`（LBFGS）、`EDIFFG = -0.03` |
| `KPOINTS` | Γ 中心 4×4×1 |
| `idpp_interpolate.py` | 基于 ASE 的 IDPP/线性插值，检查端点一致性，写出 00 … NN（`-h` 查看选项） |
| `make_images.sh` | 插值 + 复制端点 OUTCAR + dist.pl 检查 |
| `job.slurm` | NEB 作业脚本（3 节点 × 32 进程 = 96 = 3 图像 × 32），提交前检查进程数能否被 IMAGES 整除 |
| `POTCAR` | **不提供**，用 `../10-co-pt111-adsorption/make_potcar.sh` 生成 |

## 运行步骤

需要 **VTST 版 VASP**（`IOPT`、`LCLIMB`、`ICHAIN` 来自 VTST 补丁）和 [vtstscripts](https://github.com/henkelmangroup/vtstscripts)（在 PATH 中），插值脚本需要 ASE。

```bash
cd VASP/examples/11-pt-adatom-neb
export VASP_PP_PATH=/path/to/your/paw/library

# 1) 端点优化
mkdir -p IS FS
cp POSCAR_IS IS/POSCAR; cp POSCAR_FS FS/POSCAR
for d in IS FS; do cp INCAR_relax $d/INCAR; cp KPOINTS $d/; done
bash ../10-co-pt111-adsorption/make_potcar.sh IS FS
cp IS/POTCAR .
(cd IS && sbatch ../../10-co-pt111-adsorption/job.slurm)
(cd FS && sbatch ../../10-co-pt111-adsorption/job.slurm)

# 2) 检查端点并生成图像 00..04
python idpp_interpolate.py IS/CONTCAR FS/CONTCAR --check-only
bash make_images.sh IS FS 3                 # IDPP；线性插值: bash make_images.sh IS FS 3 linear

# 3) 第一阶段：普通 NEB
cp INCAR_neb_stage1 INCAR
sbatch job.slurm
nebef.pl                                    # 运行中随时查看各图像受力与能量

# 4) 第二阶段：CI-NEB
vfin.pl stage1                              # 保存结果并把 CONTCAR 覆盖到 POSCAR
cp INCAR_neb_stage2 INCAR
sbatch job.slurm

# 5) 能垒
nebresults.pl
cat exts.dat
```

## 预期结果 / 如何检查

- `idpp_interpolate.py --check-only`：只有第 37 号原子位移约 2.8 Å，其余原子位移约 0；
- `make_images.sh` 打印的相邻图像间距约 0.7 Å，且都相等；
- 路径对称：图像 02（bridge 位）能量最高，01 与 03 能量相近，能垒曲线左右对称；
- CI-NEB 收敛后所有图像受力 < 0.03 eV/Å；`exts.dat` 给出的极大值应与图像 02 的能量基本一致；
- 正、逆向能垒相同（IS 与 FS 等价）。

## 注意事项

- Pt 在 Pt(001) 上的自扩散实验上以交换机制为主（教程 §12 的说明），本例的跳跃路径用于教学；
- 晶格常数、层数、表面胞、k 点只是教学设置；正式计算需做收敛测试（教程 §15）；
- `IMAGES` 改变时，`job.slurm` 中的总进程数也要改为新 IMAGES 的整数倍；
- 端点、NEB、Dimer、频率计算的电子部分参数必须完全一致；
- 在 Windows/Git Bash 测试时若 `python3` 指向应用商店存根，可用 `PYTHON=python bash make_images.sh ...`。
