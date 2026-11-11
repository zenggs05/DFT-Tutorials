# 示例 11-4：过渡态验证——爬坡图像的有限差分频率分析

## 目的

- 对 CI-NEB（或 Dimer）得到的过渡态做 `IBRION = 5` 频率计算，确认"有且仅有一个虚频"；
- 检查虚频本征矢是否沿反应方向；
- 对初态做同样的计算，估算谐振子过渡态理论的 Vineyard 指前因子。

对应教程：[第11章 过渡态搜索：NEB / CI-NEB / Dimer](../../tutorials/11-neb-dimer.md)（§9）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `INCAR` | 电子部分与 NEB 相同、`EDIFF = 1E-7`；`IBRION = 5`、`POTIM = 0.015`、`NFREE = 2`、`NSW = 1` |
| `KPOINTS` | 与示例 11-1 相同（Γ 中心 4×4×1） |
| POSCAR | 由爬坡图像的 CONTCAR 经 `set_freq_selective.py` 生成（只放开反应相关原子） |

`IBRION = 5` 是唯一支持 Selective dynamics 的有限差分模式；`IBRION = 6` 会忽略它并计算全部原子。

## 运行步骤

以示例 11-1 为例（吸附原子是第 37 号原子，表面层 z ≈ 13.96 Å）：

```bash
cd VASP/examples/11-pt-adatom-neb
mkdir -p ts_freq is_freq
SEL=../10-co-pt111-adsorption/set_freq_selective.py

# 快速检查：只放开吸附原子（3 个模，6 次 SCF）
python $SEL 02/CONTCAR -i 37 -o ts_freq/POSCAR
# 更完整：吸附原子 + 表面层（10 个原子，30 个模）
# python $SEL 02/CONTCAR -i 37 --zmin 13.0 -o ts_freq/POSCAR

python $SEL IS/CONTCAR -i 37 -o is_freq/POSCAR        # 初态必须放开相同的原子
for d in ts_freq is_freq; do cp ../11-ts-freq-check/INCAR ../11-ts-freq-check/KPOINTS POTCAR $d/; done
(cd ts_freq && sbatch ../../10-co-pt111-adsorption/job.slurm)
(cd is_freq && sbatch ../../10-co-pt111-adsorption/job.slurm)

grep -E "f  =|f/i=" ts_freq/OUTCAR | head -3
grep -E "f  =|f/i=" is_freq/OUTCAR | head -3
```

注意 OUTCAR 中频率列表会打印两遍（第二遍在 "Eigenvectors after division by SQRT(mass)" 之后），只看前 3N_free 行。

## 预期结果 / 如何检查

- `ts_freq`：恰好 1 个明显的 `f/i`（虚频），其本征矢主要是吸附原子沿跳跃方向（x）的运动；另外 2 个为实频；
- `is_freq`：3 个都是实频（极小点）；
- 若 TS 出现 2 个明显虚频，说明是高阶鞍点或未收敛；没有虚频说明不是鞍点；
- 指前因子：$\nu^*=\prod_i\nu_i^\mathrm{IS}/\prod_{j\neq\mathrm{imag}}\nu_j^\mathrm{TS}$（单位 THz 或 s⁻¹），典型量级 $10^{12}$–$10^{13}$ s⁻¹。

快速计算指前因子（把频率换成你 OUTCAR 中的数值，单位 THz）：

```python
import numpy as np
nu_is = np.array([3.0, 2.5, 2.4])   # placeholder numbers - replace with is_freq/OUTCAR values (THz)
nu_ts = np.array([2.6, 2.2])        # real TS modes only, imaginary one removed (THz)
print("nu* = %.3e s^-1" % (np.prod(nu_is) / np.prod(nu_ts) * 1e12))
```

## 注意事项

- 只放开吸附原子是最省时的检查，但会高估/低估部分频率；发表数据时应放开更多近邻原子并检查收敛；
- 小于约 50 cm⁻¹ 的虚频多为数值噪声：提高 `EDIFF`、检查 TS 是否收敛、尝试 `POTIM = 0.01/0.02`；
- 振动模式可视化：[VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF)；
- VTST 的动力学矩阵方法（`ICHAIN = 1` + `dymmatrix.pl`、`dymprefactor.pl`）是另一种选择，见教程 §9.2。
