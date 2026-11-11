# 算例 7-1：H₂O 分子的 NVE 分子动力学（能量守恒与振动频率）

## 目的

用一个孤立水分子做 OpenMX 4.0 的 NVE 第一性原理分子动力学：一个 O–H 键略微拉长（1.02 Å），另一个 0.97 Å，从静止出发让分子自由振动。练习 `MD.Type NVE`、`MD.TimeStep`、`MD.maxIter`、`<MD.Init.Velocity` 的用法，学会从 `H2O_NVE.ene` 判断守恒量 `Utot + Ukc` 是否守恒，并从 `.md` 轨迹的键长/键角傅里叶变换得到 O–H 伸缩与 H–O–H 弯曲振动频率。

对应教程：[第7章 第一性原理分子动力学](../../tutorials/07-md.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `H2O_NVE.dat` | OpenMX 输入：10 Å 立方盒子、Cluster 求解器、GGA-PBE、Standard 基组，NVE，0.5 fs × 1000 步 |
| `gen_init_velocity.py` | 读取 `.dat`，按目标温度生成麦克斯韦–玻尔兹曼初速度（m/s），去除质心运动，输出 `<MD.Init.Velocity` 块 |
| `nve_check.py` | 读取 `H2O_NVE.ene` 与 `H2O_NVE.md`：平均温度、守恒量漂移（meV/ps）与峰峰值、键长/键角的振动谱峰位 |
| `job.slurm` | Slurm 作业脚本（4 个 MPI 进程），运行后自动调用 `nve_check.py` |

## 运行步骤

1. 修改 `H2O_NVE.dat` 中的 `DATA.PATH`，指向同时包含 `PAO/` 与 `VPS/` 的 `DFT_DATA19` 目录（官方 deb 包默认在 `/opt/openmx4.0/DFT_DATA19`）。
2. 运行：

   ```bash
   cd OpenMX/examples/07-md-h2o-nve
   rm -f H2O_NVE.ene H2O_NVE.md          # OpenMX appends to these files
   mpirun -np 4 openmx H2O_NVE.dat -nt 1 | tee H2O_NVE.std
   # or: sbatch job.slurm
   ```

3. 分析：

   ```bash
   python3 nve_check.py H2O_NVE --plot     # writes nve_check.png
   ```

4. （可选）给定 300 K 初速度后重跑：

   ```bash
   cp H2O_NVE.dat H2O_NVE_v300.dat
   sed -i 's/^System.Name .*/System.Name   H2O_NVE_v300/' H2O_NVE_v300.dat
   python3 gen_init_velocity.py H2O_NVE_v300.dat --temp 300 --seed 1 >> H2O_NVE_v300.dat
   mpirun -np 4 openmx H2O_NVE_v300.dat -nt 1 > H2O_NVE_v300.std
   python3 nve_check.py H2O_NVE_v300
   ```

## 预期结果 / 如何检查

- 运行中：`grep -c "time=" H2O_NVE.md` 显示已完成步数；`grep "scf_conv=0" H2O_NVE.md` 应为空（每步 SCF 都收敛）。
- `nve_check.py` 输出格式（示意输出，数值仅供参考）：

  ```text
  steps = 1000, dt = 0.500 fs
  <T> = 1xx.x K
  E_cons drift = x.xxx meV/ps, peak-to-peak = x.xxx meV
  Ukc peak-to-peak = xxx.x meV
  bond 1-2       mean =   0.98xx  spectrum peaks near 3xxx cm^-1
  bond 1-3       mean =   0.98xx  spectrum peaks near 3xxx cm^-1
  angle 2-1-3    mean = 10x.xxxx  spectrum peaks near 1xxx cm^-1
  ```

- 判据：守恒量峰峰值应远小于 `Ukc` 的峰峰值（例如不到其 1–2%）且无明显单调漂移；O–H 伸缩峰在约 3500–3900 cm⁻¹，弯曲峰在约 1500–1650 cm⁻¹（实验气相基频约 3657、3756 与 1595 cm⁻¹；PBE 误差、非谐性和 500 fs 轨迹约 67 cm⁻¹ 的频率分辨率都会带来偏差）。
- 若漂移明显：把 `MD.TimeStep` 减到 0.25 fs，或提高 `scf.energycutoff`、收紧 `scf.criterion` 后对比。

## 注意事项

- 物种名用元素符号（`H`、`O`），并让 `System.Name` 与文件名主干一致，便于后续用 dpdata 的 `openmx/md` 格式读取。
- 3 原子体系的瞬时温度涨落非常大（约 ±47%），NVE 中平均温度也会因动能与势能均分而变化，这是正常的。
- 想用 1.0 fs 时间步，可在物种定义第 4 列把 H 质量设为 2.0（氘），但此时得到的是 D₂O 的振动频率（约为 H₂O 的 1/√2 倍量级）。
- 输入文件只含 ASCII 字符；在 Windows 编辑后请确认是 LF 换行、无 BOM。

返回教程：[第7章 第一性原理分子动力学](../../tutorials/07-md.md)
