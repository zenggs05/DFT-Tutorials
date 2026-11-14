# 算例 14-1：Si 的频率相关介电函数与光学常数（LOPTICS）

## 目的

用独立粒子近似（IPA）的 `LOPTICS = .TRUE.` 计算体相 Si 的复介电函数 ε(ω) = ε1 + iε2，并由此得到折射率 n、消光系数 k、吸收系数 α、反射率 R、能量损失谱，以及 Tauc 图。演示"两步法"：先自洽，再读 WAVECAR、用大量空带一步精确对角化并计算光学矩阵元。

对应教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | Si 金刚石结构原胞，a = 5.47 Å（接近 PBE 平衡值） |
| `KPOINTS` | Γ 中心 16×16×16；光谱比总能更依赖 k 点密度 |
| `INCAR.scf` | 第 1 步：常规自洽，输出 WAVECAR |
| `INCAR.optics` | 第 2 步：`ISTART = 1`、`ALGO = Exact`、`NELM = 1`、`NBANDS = 32`、`LOPTICS = .TRUE.`、`CSHIFT = 0.1`、`NEDOS = 2000` |
| `job.slurm` | 依次运行两步并调用 `plot_optics.py` |
| `plot_optics.py` | 读取 vasprun.xml 中的 `<dielectricfunction>`（默认 density-density），输出 `optics.dat`、`optics.png`，可选 Tauc 图 |

## 运行步骤

```bash
cd VASP/examples/14-optics-si-loptics
export VASP_PP_PATH=/path/to/your/paw_library
sbatch job.slurm
```

手动执行：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
cp INCAR.scf INCAR    && mpirun -np 16 vasp_std > scf.out
cp INCAR.optics INCAR && mpirun -np 16 vasp_std > optics.out
python3 plot_optics.py                          # 各向同性平均，density-density
python3 plot_optics.py --tauc direct --fit 2.6 3.2   # 示例拟合窗口，请按自己的谱调整
```

## 预期结果 / 如何检查

1. 第 2 步 OUTCAR 中应能找到：

   ```bash
   grep "NBANDS=" OUTCAR                               # 确认实际能带数（并行可能把它调大）
   grep -n "DIELECTRIC FUNCTION" OUTCAR | head         # IMAGINARY/REAL, density-density/current-current
   grep "<dielectricfunction" vasprun.xml              # VASP 6 一般有 density-density 与 current-current 两组
   ```

2. `plot_optics.py` 打印的 `eps1(omega->0)` 是静态电子介电常数（不含局域场效应）。PBE 通常高估 Si 的 ε∞（实验约 11.9，近似值），这是因为 PBE 带隙偏小。可与 [算例 14-2](../14-dielectric-gaas-lepsilon/README.md) 同样方法对 Si 做一次 `LEPSILON` 计算，比较其中"excluding local field effects"的值；两者应接近，否则说明 NBANDS 不够。
3. ε2 的主峰对应 Si 的 E1、E2 带间跃迁（实验约 3.4 eV 与 4.2 eV，近似值）；PBE 由于带隙低估，峰位整体偏低约 0.5 eV。
4. Si 是间接带隙半导体，IPA 不含声子辅助跃迁，吸收边对应的是直接跃迁；对 Si 用 Tauc 图提取的"带隙"不是实验意义上的间接带隙，这里只用于演示脚本。

## 注意事项

- 收敛检查：分别增加 `NBANDS`（24 → 32 → 48）和 k 点（12 → 16 → 20），直到 ε1(0) 和主峰位置不再变化。
- `SIGMA` 应远小于 `CSHIFT`，否则谱线同时受到高斯和洛伦兹两种展宽。
- `LOPTICS` 与 `LEPSILON` 不要放在同一次计算中。
- 若想一步完成，也可以在单个 INCAR 中直接设置 `NBANDS = 32` + `LOPTICS = .TRUE.`（`ALGO = Normal`），对小体系同样可行。
- 也可以用 py4vasp（`calc.dielectric_function.plot()`）、VASPKIT 711 或 `sumo-optplot` 后处理；仓库的 [optics_post.py](../../scripts/optics_post.py) 提供更多选项，用 `python ../../scripts/optics_post.py -h` 查看用法。

返回教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)
