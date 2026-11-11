# 算例 09-cu100-surface-energy：Cu(100) 表面能与层数收敛

## 目的

- 计算对称 slab 的表面能 $\gamma = [E_{\rm slab}(N) - N\varepsilon_{\rm bulk}]/(2A)$；
- 比较两种体相参考：取向体相（oriented bulk，与 slab 同样的面内 k 网格）与 Fiorentini–Methfessel 线性拟合；
- 比较未弛豫（unrelaxed）与弛豫（relaxed）表面能，观察 γ 随层数 N 的收敛。

对应教程：[第9章](../../tutorials/09-surface-workfunction.md) §6。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR_bulk` / `KPOINTS_bulk` / `INCAR_bulk` | fcc Cu 原胞的 `ISIF=3` 晶格优化（`ENCUT=500`，`Cu_pv`） |
| `INCAR_slab` | 对称 slab 弛豫：`ISIF=2`，中心层固定，上下两侧同时弛豫 |
| `INCAR_static` | 取向体相与未弛豫 slab 的单点能 |
| `make_slabs.py` | 读 `bulk/CONTCAR` 得到晶格常数，用 ASE 生成 `obulk/`、`slab_N/`、`slab_N_unrel/` 的 POSCAR 与 KPOINTS |
| `surface_energy.py` | 汇总能量，输出 γ（eV/Å² 与 J/m²）及线性拟合结果 |
| `run.sh` | `bulk` → `slabs` → `collect` |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR（VASP 推荐 `Cu_pv`）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Cu_pv/POTCAR > POTCAR
```

## 运行步骤

```bash
bash run.sh bulk          # 1) bulk lattice constant
# wait; optionally: cd bulk && cp CONTCAR POSCAR && sbatch job.slurm  (re-relax once)
bash run.sh slabs         # 2) obulk + N = 3,5,7,9,11 relaxed and unrelaxed slabs (11 jobs)
# wait ...
bash run.sh collect       # 3) surface energies
```

`make_slabs.py` 需要 ASE（`pip install ase`）。

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" slab_*/OUTCAR
python surface_energy.py
```

- 取向体相参考下，γ(N) 随 N 很快收敛：N ≥ 7 时相邻层数的差别通常 < 0.01 J/m²（示意）。
- 线性拟合得到的 ε_fit 与 obulk 的 ε_bulk 之差应在 1 meV/atom 以内；若差别明显，说明 k 点或 smearing 在 slab 与体相之间不一致——这正是直接用"原胞体相能量"时 γ 随 N 线性漂移的原因。
- 弛豫后的 γ 略低于未弛豫的 γ（金属表面弛豫能通常只有百分之几）。
- PBE 计算的 Cu(100) 表面能在 1.4–1.5 J/m² 量级（近似值，参见 Tran et al., Sci. Data 3, 160080 (2016) 的 Crystalium 数据库）；实验"平均表面能"约 1.8 J/m²，PBE 系统性低估表面能是已知现象。

## 注意事项

- slab 与体相参考必须使用**完全相同**的 ENCUT、POTCAR、PREC、LREAL、smearing（ISMEAR/SIGMA）和面内 k 网格；能量取 `energy(sigma->0)`。
- 本例使用奇数层并固定中心层，保证 slab 上下对称，可以直接除以 2A；若只弛豫一侧（底部固定），需要用教程 §6.3 中的"未弛豫 + 弛豫"组合公式。
- 1 eV/Å² = 16.02 J/m²。
- 层数测试同时也是真空测试的好时机：把 `--vacuum` 改为 20 Å 再算一组，γ 变化应 < 0.005 J/m²。

[返回第9章](../../tutorials/09-surface-workfunction.md)
