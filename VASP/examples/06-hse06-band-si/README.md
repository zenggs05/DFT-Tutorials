# 示例 06-1：Si 的 HSE06 自洽与能带（零权重 k 点法 + KPOINTS_OPT 法）

## 目的

- 从 PBE 波函数出发完成 HSE06 自洽，读取 HSE06 带隙；
- 用两种正确方法计算 HSE06 能带：零权重 k 点法（所有版本）与 KPOINTS_OPT 法（VASP ≥ 6.3）；
- 体会杂化泛函的计算代价，以及 `HFRCUT = -1`、`PRECFOCK`、`NKRED` 的作用。

对应教程：[第6章 杂化泛函与 meta-GGA](../../tutorials/06-hybrid-metagga.md)（§4–§6）。

## 文件说明

| 路径 | 说明 |
|---|---|
| `POSCAR` | Si 原胞（a = 5.47 Å，建议换成你优化得到的值） |
| `0-pbe-scf/INCAR, KPOINTS` | PBE 自洽，Γ 中心 6×6×6，`NBANDS=16`，写出 WAVECAR 和 IBZKPT |
| `1-hse-scf/INCAR, KPOINTS` | HSE06 自洽：`LHFCALC`、`HFSCREEN=0.2`、`ALGO=Damped`、`PRECFOCK=Normal`、`HFRCUT=-1`，从 PBE WAVECAR 开始 |
| `2-hse-band-zero-weight/INCAR.pbe, INCAR.hse, KPATH, make_zero_weight_kpoints.py` | 零权重法：脚本把 `0-pbe-scf/IBZKPT` 与 `KPATH` 中的路径（每段 12 点，权重 0）拼成 KPOINTS；先用同一 KPOINTS 跑 PBE，再跑 HSE06 |
| `3-hse-band-kpoints-opt/INCAR, KPOINTS, KPOINTS_OPT` | KPOINTS_OPT 法：规则网格 + 路径（每段 30 点） |
| `job.slurm` | 依次运行步骤 0–3，可用 `STEPS` 变量只运行其中几步 |
| `POTCAR` | **不提供**，在本目录生成一次 |

## 运行步骤

```bash
cd VASP/examples/06-hse06-band-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

sbatch job.slurm                                   # all steps 0 1 2 3
# or only PBE + KPOINTS_OPT method:
sbatch --export=ALL,STEPS="0 3" job.slurm
```

手动运行零权重法的关键命令：

```bash
cd 2-hse-band-zero-weight
cp ../POSCAR ../POTCAR .
python make_zero_weight_kpoints.py --ibzkpt ../0-pbe-scf/IBZKPT --kpath KPATH --out KPOINTS
cp INCAR.pbe INCAR && mpirun -np 16 vasp_std > vasp_pbe.out
cp INCAR.hse INCAR && mpirun -np 16 vasp_std > vasp_hse.out
```

也可以用 `vaspkit -task 251` 或 `sumo-kgen --hybrid`（目录中需有 IBZKPT）代替 `make_zero_weight_kpoints.py`。

## 预期结果 / 如何检查

```bash
grep -A12 "Band structure" 1-hse-scf/OUTCAR              # HSE06 gap on the 6x6x6 mesh (VASP >= 6.4.3)
grep -c "F=" 1-hse-scf/OSZICAR                           # should be 1 (single point)
grep "FFT grid for exact exchange" -A2 1-hse-scf/OUTCAR  # grid chosen by PRECFOCK
```

作图：

```python
import py4vasp
calc = py4vasp.Calculation.from_path("3-hse-band-kpoints-opt")
calc.band.plot("kpoints_opt")
print(calc.bandgap.fundamental())
```

零权重法的结果需要剔除有权重的网格点后作图：`vaspkit -task 252`、`sumo-bandplot` 或本仓库 `python ../../scripts/plot_band.py -h`。

- HSE06 的 Si 间接带隙应在约 1.1–1.3 eV（文献量级；实验约 1.17 eV），明显大于 PBE 的约 0.5–0.6 eV。
- 两种方法得到的能带应基本重合（零权重法每段点数较少，曲线更粗糙）。

## 注意事项

- **不要**对杂化泛函使用 `ICHARG = 11`。
- PBE 与 HSE 计算的 `ENCUT`、`NBANDS`、k 点集合必须一致，WAVECAR 才能复用；这也是零权重法要在本目录先跑一次 PBE 的原因。
- 零权重法不要从已收敛的 HSE WAVECAR 开始（自洽可能在零权重点收敛之前停止）。
- 打开 `NKRED = 2` 前先比较带隙变化；金属与掺杂体系不要用。
- 6.6.x NVIDIA GPU 版本上杂化泛函 + `ALGO = Normal` 有随机崩溃的已知问题，本例使用 `ALGO = Damped`。
- 计算代价：步骤 1/3 比 PBE 慢一到两个数量级，步骤 2 更贵（零权重点每个电子步都要计算）。
