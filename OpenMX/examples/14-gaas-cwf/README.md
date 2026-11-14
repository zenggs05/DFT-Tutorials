# 示例 14-3：GaAs 的最接近 Wannier 函数（CWF，AO 与 HO 引导）

配套教程：[第14章 Wannier 函数：MLWF 与 CWF](../../tutorials/14-wannier.md)

## 目的

练习 OpenMX 4.0 新增的 CWF：

| 输入 | 引导函数 | 演示内容 |
|---|---|---|
| `GaAs_CWF_AO.dat` | Ga-4s/4p、As-4s/4p 的 PAO（8 个） | 多项式窗口排除 Ga-3d/As-3d；DM_func、CWF 布居、跳跃积分衰减、cube |
| `GaAs_CWF_HO.dat` | 同样 PAO 构成的在位杂化轨道（HO） | 与 AO 张成同一子空间、插值能带相同，但 CWF 沿键方向伸展（键轨道图像） |

`scan_cwf_window.sh` 用 `scf.restart pp` 重用 SCF 结果，扫描几组多项式窗口并记录 DM_func。

## 文件说明

| 文件 | 说明 |
|---|---|
| `GaAs_CWF_AO.dat` | SCF（11×11×11，奇数）+ AO 引导 CWF + `CWF.Dis.vs.H` + cube 输出 |
| `GaAs_CWF_HO.dat` | 同上，`CWF.Guiding.Orbitals HO` |
| `cwf_analyze.py` | 读取 `.out` 中的 `Populations evaluated by CWF` 段（DM_func、布居）和 `.CWF.Dis_vs_H`（在位能、按距离壳层的最大 \|t\|），可选作图 |
| `scan_cwf_window.sh` | 窗口扫描（需要先运行 `GaAs_CWF_AO.dat` 生成 `gaas_cwf_ao_rst/`），结果写入 `cwf_window_scan.txt` |
| `job.slurm` | Slurm 作业脚本 |

能带比较复用 `../14-si-mlwf-sp3/band_compare.py`。

基组与赝势：`Ga7.0-s3p2d2` + `Ga_PBE19`（13 个价电子，含 3d）、`As7.0-s3p2d2` + `As_PBE19`（15 个价电子，含 3d）。`<CWF.Guiding.AOs` 的 7 项对应 `1s 2s 3s 1p 2p 1d 2d`，其中 `1s`、`1p` 是 4s、4p。

## 运行步骤

```bash
cd OpenMX/examples/14-gaas-cwf
sed -i "s#^DATA.PATH .*#DATA.PATH                   $HOME/openmx4.0/DFT_DATA19#" *.dat
export OMP_NUM_THREADS=1
mpirun -np 8 openmx GaAs_CWF_AO.dat -nt 1 > gaas_cwf_ao.std 2>&1
mpirun -np 8 openmx GaAs_CWF_HO.dat -nt 1 > gaas_cwf_ho.std 2>&1
# optional: window scan re-using gaas_cwf_ao_rst/
NP=8 bash scan_cwf_window.sh
# or on a cluster
sbatch job.slurm
```

需要 Wannier90 文件时，把 `CWF.Wannier90.fileout` 改为 `on` 再运行，会得到 `gaas_cwf_ao_0.amn/.mmn/.eig/.win`（`.win` 中 `num_iter 0`，即 Wannier90 直接沿用 CWF 的规范）。

## 预期结果 / 如何检查

```bash
grep "DM_func" gaas_cwf_ao.std gaas_cwf_ao.out | tail -n 2
python3 cwf_analyze.py gaas_cwf_ao
python3 cwf_analyze.py gaas_cwf_ho
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ao.Band gaas_cwf_ao.CWF.Band --emin -13.8 --emax 1.0 --plot gaas_ao.png
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ho.Band gaas_cwf_ho.CWF.Band --emin -13.8 --emax 1.0
ls gaas_cwf_ao.cwf_0_*_r.cube
```

判断标准（定性；仓库不附参考输出）：

- 每个 CWF、每个自旋的 DM_func 越小越好；手册中以全部价态为引导的 GaAs 例子为 0.16，大于 0.7 说明引导函数与窗口内态不匹配；
- 在权重为 1 的区间（−13.8 ~ 1 eV），`.CWF.Band` 与 `.Band` 重合；AO 与 HO 两次计算的插值能带相同；
- `cwf_analyze.py` 列出的最大 \|t\| 随距离迅速减小；
- CWF 布居的 `total` 明显小于 `ideal(neutral)=28`，这是因为 3d 半芯电子不在所选子空间中，不是错误；
- AO 的 cube 保持 s/p 形状；HO 的 cube 指向近邻原子；
- `scan_cwf_window.sh` 中排除 As-4s 带的窗口（−10 −9 …）给出明显更大的 DM_func。

## 注意事项

- CWF 是 OpenMX 4.0 的新功能，3.9 不支持。
- CWF 使用 `scf.Kgrid` 作为 k 网格，目前只支持奇数 k 点；`CWF.Plot.SuperCells` 要满足 2n+1 ≤ `scf.Kgrid`。
- Ga-3d 带的位置随基组和泛函变化。先看 PDOS，确认它整条带都在 ε0 = −14.5 eV 之下，As-4s 带底在 ε1 = −13.8 eV 之上。作为参照：官方参考输出 `work/cwf_example/GaAs-CWF1.out`（PBE，s3p2d2f1 基组，a = 5.65 Å）的 SCF 本征值换算到化学势为零点后，Ga-3d 上沿约 −15.2 eV，As-4s 带约 −13.1 ~ −10.5 eV。
- 对接 Wannier90 时用 `Poly` 窗口；`Fermi` 窗口会把所有 Bloch 态都包含进来，`.mmn` 可能非常大。
- `scf.restart pp` 要求 `System.Name` 与生成重启文件时相同，所以脚本每步都会把输出改名保存（`gaas_cwf_ao_w<i>.*`），并在开始时把原始计算的输出备份为 `gaas_cwf_ao_orig.*`、结束时恢复。若你的版本在 `pp` 模式下不输出 CWF 结果，改用 `scf.restart on`。
- `<CWF.Guiding.AOs` 每个物种一行，每个径向函数都写 0 或 1；块内不要写注释行。

[返回第14章](../../tutorials/14-wannier.md)
