# 示例 06-2：纤锌矿 ZnO 的 r2SCAN 变胞优化

## 目的

- 演示 meta-GGA（r2SCAN）计算的标准设置：`METAGGA = R2SCAN`、`LASPH`、`LMIXTAU`；
- 检查 POTCAR 是否适用于 meta-GGA；
- 比较 r2SCAN 与 PBE 优化得到的晶格常数（练习）。

对应教程：[第6章 杂化泛函与 meta-GGA](../../tutorials/06-hybrid-metagga.md)（§8）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | r2SCAN 变胞优化：`ENCUT=520`、`ISIF=3`、`ALGO=Normal`、`LMIXTAU=.TRUE.`；注释中给出 rVV10 与 D4 色散修正的写法 |
| `POSCAR` | 纤锌矿 ZnO，取实验结构为初始值（a = 3.250 Å、c = 5.207 Å、u = 0.382） |
| `KPOINTS` | Γ 中心 8×8×5（六角晶格必须 Γ 中心） |
| `job.slurm` | 先检查 POTCAR 是否含动能密度信息，再做两轮优化（第二轮消除 Pulay 应力） |
| `POTCAR` | **不提供**，按 `Zn`、`O` 顺序拼接 |

## 运行步骤

```bash
cd VASP/examples/06-r2scan-relax-zno
cat $VASP_PP_PATH/potpaw_PBE/{Zn,O}/POTCAR > POTCAR
grep TITEL POTCAR
grep kinetic POTCAR | sort | uniq -c        # both elements must provide kinetic-energy-density data
sbatch job.slurm
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" round*/OUTCAR
grep "external pressure" OUTCAR | tail -1
python -c "from ase.io import read; a=read('CONTCAR'); p=a.cell.cellpar(); print('a = %.4f A, c = %.4f A, c/a = %.4f' % (p[0], p[2], p[2]/p[0]))"
```

- 第二轮应在很少的离子步内结束，体积变化 < 0.2%。
- 实验晶格常数约 a = 3.25 Å、c = 5.21 Å；PBE 通常高估 1–2%，r2SCAN 一般更接近实验。具体数值请以你自己的计算为准。
- 练习：把 `METAGGA = R2SCAN` 删掉（即 PBE）重做一次，比较晶格常数与耗时。

## 注意事项

- 若 `grep kinetic POTCAR` 没有输出，说明该 POTCAR 不适用于依赖动能密度的 meta-GGA，请换用合适的 PAW。
- 收敛困难时：先用 PBE 跑一次得到 WAVECAR，再以 `ISTART = 1` 启动 r2SCAN；或改用 `ALGO = All`（VASP ≥ 6.5.1 可加 `ISEARCH = 1`）。
- 自旋极化（`ISPIN = 2`）的 meta-GGA 受力在 **VASP 6.6.0** 中有错误（known issue 94），请用 6.6.1。ZnO 是非磁性的，但做磁性体系时务必注意。
- r2SCAN 的能带不能用 `ICHARG = 11`，请用 KPOINTS_OPT 或零权重法（见第6章 §8.3）。
- r2SCAN+rVV10 的应力在 VASP 6.2 及更早版本中有误，变胞优化请用 ≥ 6.3.0；D4 需要编译时链接 DFT-D4 库。
