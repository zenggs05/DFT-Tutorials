# 示例 06-3：Si 的 mBJ 带隙与能带（METAGGA = MBJ + KPOINTS_OPT）

## 目的

- 演示 mBJ（Tran–Blaha 修正的 Becke–Johnson 势）的标准两步流程：PBE 得到 WAVECAR → mBJ 自洽；
- 在 mBJ 这一步用 KPOINTS_OPT 同时得到能带路径（c 自洽时最方便的做法）；
- 读取自洽得到的参数 c（`CMBJ`）。

对应教程：[第6章 杂化泛函与 meta-GGA](../../tutorials/06-hybrid-metagga.md)（§9）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR.pbe` | PBE 自洽，写出 WAVECAR；`LKPOINTS_OPT = .FALSE.` 让这一步忽略 KPOINTS_OPT |
| `INCAR.mbj` | `METAGGA = MBJ`、`LASPH`、`ALGO = All`、`TIME = 0.1`、`NELM = 200`，从 PBE WAVECAR 开始 |
| `POSCAR` | Si 原胞（a = 5.47 Å） |
| `KPOINTS` | Γ 中心 6×6×6（两步相同） |
| `KPOINTS_OPT` | 能带路径 L–Γ–X–U\|K–Γ，每段 30 点 |
| `job.slurm` | 依次运行 PBE 与 mBJ，结束后打印 CMBJ 与带隙 |
| `POTCAR` | **不提供** |

## 运行步骤

```bash
cd VASP/examples/06-mbj-gap-si
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR
grep kinetic POTCAR                 # meta-GGA needs kinetic-energy-density data
sbatch job.slurm
```

手动运行：

```bash
cp INCAR.pbe INCAR && mpirun -np 16 vasp_std > vasp_pbe.out
cp INCAR.mbj INCAR && mpirun -np 16 vasp_std > vasp_mbj.out
```

## 预期结果 / 如何检查

```bash
grep CMBJ OUTCAR | tail -1               # converged c (printed every electronic step)
grep -A12 "Band structure" OUTCAR        # gap on the SCF mesh (VASP >= 6.4.3)
grep -c "F=" OSZICAR                     # 1: single point
```

```python
import py4vasp
calc = py4vasp.Calculation.from_path(".")
calc.band.plot("kpoints_opt")
print(calc.bandgap.fundamental())
```

- Tran & Blaha (2009) 报告 mBJ 下 Si 的带隙与实验值（约 1.17 eV）很接近；你得到的值应明显大于 PBE（约 0.5–0.6 eV），落在 1.1–1.2 eV 量级。
- mBJ 的电子步数通常明显多于 PBE（几十到上百步）。

## 注意事项

- mBJ 只是交换**势**，没有对应的能量泛函：OUTCAR 中的总能和受力没有物理意义，**不能**用它优化结构或比较能量。
- c 自洽时无法直接从 WAVECAR 续算；需要另算能带或续算时，把 `grep CMBJ OUTCAR | tail -1` 得到的值写入 `CMBJ` 固定。
- slab、界面等含真空的体系必须固定 `CMBJ`（取体相值）或改用 `METAGGA = LMBJ`。
- 收敛困难时按官方建议分三步：PBE → 固定 `CMBJ` 的 mBJ → 不固定 `CMBJ` 的 mBJ，每一步从上一步的 WAVECAR 开始。
- mBJ 的能带同样不能用 `ICHARG = 11`；本例用 KPOINTS_OPT（需 VASP ≥ 6.3、`vasp_std`）。
