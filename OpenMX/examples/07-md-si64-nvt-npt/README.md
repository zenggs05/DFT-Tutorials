# 算例 7-2：Si 的 NVT（速度标度 / Nosé–Hoover）与 NPT 分子动力学

## 目的

用 64 原子 Si 超胞（2×2×2 金刚石惯用胞）练习 OpenMX 4.0 的恒温分子动力学：先用 `NVT_VS` 把体系从 300 K 升温到 1000 K 并平衡，再用 `NVT_NH`（Nosé–Hoover）做生产运行；然后用 8 原子 Si 惯用胞体验 4.0 新增的恒温恒压 `NPT_VS_PR`（600 K，外压 0 → 5 GPa）。最后用 `md_analysis.py` 计算温度、守恒量漂移、径向分布函数 g(r)、配位数、均方位移 MSD 与 NPT 平均体积，并可用 dpdata 导出 DeePMD 训练数据。

对应教程：[第7章 第一性原理分子动力学](../../tutorials/07-md.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Si64_NVT_VS.dat` | 64 原子 Si，速度标度 NVT：10 步 300 K → 400 步升到 1000 K → 保温到 1500 步；2 fs/步 |
| `Si64_NVT_NH.dat` | 同一体系的 Nosé–Hoover NVT（1000 K，2000 步，`NH.Mass.HeatBath 30`） |
| `Si8_NPT_VS_PR.dat` | 8 原子 Si 的 NPT（速度标度 + Parrinello–Rahman），600 K，前 300 步 0 GPa，之后线性升到 5 GPa；1 fs/步 |
| `md_analysis.py` | 读取 `.ene`、`.md`、`.npt.ene`：⟨T⟩、守恒量漂移、g(r) 与配位数、MSD 与扩散系数、NPT 平均体积；`--plot` 作图，`--dpdata` 导出 `deepmd/npy` |
| `job.slurm` | Slurm 脚本，用环境变量 `NAME` 选择要运行的输入 |

所有坐标都在理想格点上加了 0.04–0.05 Å 的随机位移（固定随机种子生成），这样第一步就有非零的力，速度标度才能工作。

## 运行步骤

先修改三个 `.dat` 中的 `DATA.PATH`。

**1）速度标度升温与平衡**

```bash
cd OpenMX/examples/07-md-si64-nvt-npt
rm -f Si64_NVT_VS.ene Si64_NVT_VS.md
mpirun -np 32 openmx Si64_NVT_VS.dat -nt 1 | tee Si64_NVT_VS.std
python3 md_analysis.py Si64_NVT_VS --skip 400 --plot
```

**2）Nosé–Hoover 生产运行（推荐从第 1 步的终点接续）**

```bash
cp Si64_NVT_VS.dat# Si64_NVT_NH_from_VS.dat
# edit Si64_NVT_NH_from_VS.dat:
#   System.Name       -> Si64_NVT_NH
#   MD.Type           -> NVT_NH
#   <MD.TempControl   -> 2-column format (copy it from Si64_NVT_NH.dat)
#   add               NH.Mass.HeatBath 30.0
#   delete the line   MD.Current.Iter ...   (start counting from 1 again)
#   delete            scf.restart on        (the _rst files belong to Si64_NVT_VS)
#   MD.maxIter        -> 2000
mpirun -np 32 openmx Si64_NVT_NH_from_VS.dat -nt 1 | tee Si64_NVT_NH.std
python3 md_analysis.py Si64_NVT_NH --skip 250 --plot
```

也可以直接运行 `Si64_NVT_NH.dat`（从静止的扰动结构出发，热浴需要更长时间才能把体系加热到 1000 K）。

**3）NPT**

```bash
rm -f Si8_NPT_VS_PR.ene Si8_NPT_VS_PR.md Si8_NPT_VS_PR.npt.ene
mpirun -np 16 openmx Si8_NPT_VS_PR.dat -nt 1 | tee Si8_NPT_VS_PR.std
grep "NPT.Mass.Barostat" Si8_NPT_VS_PR.std        # default piston mass printed by OpenMX
python3 md_analysis.py Si8_NPT_VS_PR --skip 200
```

**4）（可选）导出机器学习势训练数据**

```bash
pip install dpdata
python3 md_analysis.py Si64_NVT_NH --skip 250 --dpdata deepmd_si64 --stride 5
```

使用 Slurm 时：`sbatch job.slurm`，或 `sbatch --export=ALL,NAME=Si64_NVT_NH job.slurm`；续算 `.dat#` 时加 `FRESH=0` 和 `INPUT=xxx.dat` 以免删除已有轨迹。

## 预期结果 / 如何检查

- 温度：`awk '{print $1, $18, $19}' Si64_NVT_VS.ene | tail` 中第 3 列（计算温度）应围绕第 2 列（目标温度）涨落，64 原子约 ±10%。
- `md_analysis.py` 输出格式（示意输出，数值仅供参考）：

  ```text
  ene: 1100 steps used (skip 400)
  <T> = 99x.x +- 9x.x K (target 1000.0 K at last step)
  drift of Utot+Ukc     =   xxx.xxx meV/ps (meaningful for NVE / NVT_NH respectively)
  drift of H_NoseHoover =   xxx.xxx meV/ps (meaningful for NVE / NVT_NH respectively)
  md: 1100 frames used, dt = 2.000 fs, 64 atoms
  g(r) Si-Si: first peak 2.3xx A, first minimum 2.9xx A, coordination 4.0x
  MSD at 1100 fs = 0.0xxx A^2; D from linear fit (second half) = x.xxxe-0x cm^2/s
  ```

- 晶体 Si 在 1000 K：g(r) 第一峰约 2.35 Å、配位数约 4，MSD 很快达到平台（扩散系数接近 0）。若配位数明显大于 4、MSD 线性增长，说明已经熔化或结构出了问题。
- `NVT_NH`：Nosé–Hoover 哈密顿量（`.ene` 第 22 列）应基本守恒（漂移远小于 `Utot+Ukc` 在 NVT 中的变化）。
- NPT：0 GPa 段的平均晶格常数应略大于 PBE 的 0 K 值（PBE 约 5.46–5.47 Å，实验 5.431 Å），加压到 5 GPa 后体积减小（Si 体模量实验值约 98 GPa，5 GPa 约对应 5% 的体积压缩，仅作量级参考）。8 原子体系的体积涨落很大，结论要靠更大超胞和更长轨迹。

## 注意事项

- OpenMX 以追加方式写 `.ene`、`.md`、`.npt.ene`，从头重跑前务必删除旧文件。
- `MD.maxIter` 是绝对步数：用 `.dat#` 续算时，文件中有 `MD.Current.Iter`，必须把 `MD.maxIter` 改成新的终止步，速度标度 `<MD.TempControl` 的最后一行步数也要覆盖到该步。
- 速度标度与 Nosé–Hoover 的 `<MD.TempControl` 格式不同（4 列 vs 2 列），改 `MD.Type` 时一并修改。
- 这里用的是 Quick 基组（`Si7.0-s2p1d1`）和 2×2×2 k 点，只为演示；正式研究请对基组、截断能、k 点做收敛测试（见第 5 章），NPT 的截断能尤其要充分。
- 物种名使用元素符号 `Si`，`System.Name` 与文件名主干一致，dpdata 才能正确配对 `.dat` 与 `.md`。

返回教程：[第7章 第一性原理分子动力学](../../tutorials/07-md.md)
