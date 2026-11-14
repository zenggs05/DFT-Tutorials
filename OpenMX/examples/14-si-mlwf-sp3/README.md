# 示例 14-1：Si 的最大局域化 Wannier 函数（sp3 型、键中心型与 Wannier90 接口）

配套教程：[第14章 Wannier 函数：MLWF 与 CWF](../../tutorials/14-wannier.md)

## 目的

用金刚石结构 Si 练习 OpenMX 的 MLWF 全流程：

| 输入 | Wannier 函数 | 演示内容 |
|---|---|---|
| `Si_MLWF_sp3.dat` | 8 个，4 价带 + 4 导带，sp3 初猜 | 纠缠能带的解纠缠：外窗口 −14 ~ 12 eV，内（冻结）窗口 −14 ~ 1 eV；局部坐标轴旋转 |
| `Si_MLWF_bond4.dat` | 4 个，只取价带，键中点 s 初猜 | 孤立能带组；投影子放在非原子位置（`Wannier.Initial.Projectors.Unit ANG`） |
| `Si_MLWF_W90.dat` | 同 sp3 | `Wannier90.fileout on`：只生成 `.amn/.mmn/.eig/.win`，局域化交给 wannier90.x |

## 文件说明

| 文件 | 说明 |
|---|---|
| `Si_MLWF_sp3.dat` | SCF + 8 个 sp3 型 MLWF + 插值能带 + cube 输出 |
| `Si_MLWF_bond4.dat` | SCF + 4 个键中心 MLWF |
| `Si_MLWF_W90.dat` | SCF + Wannier90 接口文件（OpenMX 自身优化步数设为 0；4.0 源码在 `Wannier90.fileout on` 时也会自动置 0） |
| `w90_extra.win` | 追加到 `si_w90.win` 的两行（`write_hr`、`write_xyz`） |
| `wf_centers.py` | 从 `.std` 中提取 DISE/SPRD/CENT 信息，打印最终中心与 spread；可计算到参考点（如键中点）的距离（仅用标准库） |
| `band_compare.py` | 比较 `.Band` 与 `.Wannier_Band`（注意大写 B，手册误写为 `.Wannier_band`；或 `.CWF.Band`），给出能量窗口内的最大/平均偏差，可选作图（需要 numpy，作图需要 matplotlib） |
| `job.slurm` | Slurm 作业脚本：依次运行三个输入并做检查 |

基组与赝势：`Si7.0-s2p2d1` + `Si_PBE19`（4 个价电子）。投影子物种 `proj` 使用 `Si7.0-s1p1d1f1`，它不放在任何原子上，所写的 VPS 也不参与计算。

## 运行步骤

1. 修改 `DATA.PATH`，指向同时包含 `PAO/` 与 `VPS/` 的目录：

   ```bash
   cd OpenMX/examples/14-si-mlwf-sp3
   sed -i "s#^DATA.PATH .*#DATA.PATH                   $HOME/openmx4.0/DFT_DATA19#" *.dat
   ```

2. 运行（每个输入都包含自己的 SCF，在 8 核上通常几分钟到十几分钟）：

   ```bash
   export OMP_NUM_THREADS=1
   mpirun -np 8 openmx Si_MLWF_sp3.dat   -nt 1 > si_mlwf_sp3.std   2>&1
   mpirun -np 8 openmx Si_MLWF_bond4.dat -nt 1 > si_mlwf_bond4.std 2>&1
   mpirun -np 8 openmx Si_MLWF_W90.dat   -nt 1 > si_w90.std        2>&1
   # or on a cluster
   sbatch job.slurm
   ```

3. Wannier90（可选，需要 `wannier90.x`）：

   ```bash
   cat w90_extra.win >> si_w90.win
   wannier90.x si_w90
   ```

4. 只想换内窗口或初猜、不重算重叠矩阵时：把 `Wannier.Readin.Overlap.Matrix` 改为 `on`，并加一行 `scf.restart on`，`System.Name`、外窗口和 `Wannier.Kgrid` 保持不变（详见教程 5.5 节）。

5. 使用 OpenMX 3.9 时，删除各输入中的 `scf.Generation.Kpoint` 一行。

## 预期结果 / 如何检查

```bash
grep -A 12 "Selected bands within the Outer and Inner Windows" si_mlwf_sp3.std
# one line per Wannier k point; minimum / maximum band counts over all k points:
grep -E '^ *[0-9] +\| *[0-9]+\| +\(' si_mlwf_sp3.std | \
  awk -F'|' '{split($3,o,"]"); split($4,w,"]"); n=o[2]+0; m=w[2]+0;
              if (NR==1 || n<nmin) nmin=n; if (m>mmax) mmax=m}
              END {print "min bands in OUTER window:", nmin, "  max bands in INNER window:", mmax}'
grep "CONVERGENCE ACHIEVED" si_mlwf_sp3.std si_mlwf_bond4.std
python3 wf_centers.py si_mlwf_sp3.std
python3 wf_centers.py si_mlwf_bond4.std \
   --cell 0 2.715 2.715  2.715 0 2.715  2.715 2.715 0 \
   --points "0.67875,0.67875,0.67875;0.67875,-0.67875,-0.67875;-0.67875,0.67875,-0.67875;-0.67875,-0.67875,0.67875"
python3 band_compare.py si_mlwf_sp3.Band   si_mlwf_sp3.Wannier_Band   --emin -14 --emax 1 --plot si_sp3.png
python3 band_compare.py si_mlwf_bond4.Band si_mlwf_bond4.Wannier_Band --emin -14 --emax 0
ls si_mlwf_sp3.mlwf0_*_r.cube
```

判断标准（定性；仓库不附参考输出，具体数值取决于版本和设置）：

- 每个 k 点外窗口内的能带数 ≥ WF 数，内窗口内的能带数 ≤ WF 数（否则程序会报错停止）；
- `.std` 中出现 `CONVERGENCE ACHIEVED !`，`DISE` 的 `Delta_I` 和 `CONV` 的梯度模低于设定判据；
- sp3：8 个 spread 相等（对称性），中心在 Si 原子与键中点的连线上；手册中同类 LDA 计算每个 WF 约 2.85 Å²，本例量级应相同；
- bond4：4 个中心到最近键中点的距离接近 0，spread 相等；
- 内窗口（sp3）或价带（bond4）范围内，`band_compare.py` 给出的 `max |dE|` 在 meV 量级；
- `si_w90.win`、`.amn`、`.mmn`、`.eig` 生成；wannier90.x 的 `si_w90.wout` 中 `Final State` 段给出 8 个中心与 spread。

## 注意事项

- 窗口能量以费米能级为零点；300 K 下半导体的费米能级在带隙中间附近，不是价带顶。
- 第二个原子的 sp3 投影子用了旋转坐标（局部 z = −x，局部 x = −z），这样它的四个瓣才指向自己的近邻。取向写错不会报错，只会让 spread 难以收敛。
- `Wannier.Function.Plot.SuperCells` 要满足 2n+1 ≤ `Wannier.Kgrid`。
- 块（如 `<Wannier.Initial.Projectors`）内部不要写注释行。
- 通用作业脚本模板见 [`job_omx_slurm.sh`](../../scripts/job_omx_slurm.sh)；能带作图可用 [`omx_band_plot.py`](../../scripts/omx_band_plot.py)，cube 处理可用 [`omx_cube_tools.py`](../../scripts/omx_cube_tools.py)（用 `-h` 查看用法）。

[返回第14章](../../tutorials/14-wannier.md)
