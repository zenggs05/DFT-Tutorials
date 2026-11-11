# 示例 17-2：GaAs 在 Γ 点的有效质量（有限差分 / 抛物线拟合）

## 目的

- 在 Γ 点附近沿 [100]、[110]、[111] 三个方向取若干 k 点，做非自洽能带计算；
- 用抛物线拟合 E(k) = E0 + ħ²k²/2m* 得到导带底电子和价带顶空穴的有效质量；
- 理解单位换算（ħ²/2m0 = 3.80998 eV·Å²）和 PBE 带隙偏小对有效质量的影响。

对应教程：[第17章 高级电子结构](../../tutorials/17-advanced-electronic.md)（§4）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 闪锌矿 GaAs 原胞，a = 5.653 Å（实验值） |
| `INCAR_scf` / `KPOINTS_scf` | 自洽，Γ 中心 8×8×8，`EDIFF = 1E-8` |
| `INCAR_nscf` | `ICHARG = 11`、`NBANDS = 16` |
| `KPOINTS_fd` | 显式 k 点表（Cartesian，单位 2π/a）：Γ + 三个方向各 5 个点，步长 0.005·2π/a ≈ 0.0056 Å⁻¹ |
| `fit_effmass.py` | 读 EIGENVAL + POSCAR，按方向分组拟合有效质量（numpy） |
| `job.slurm` | 自洽 → 非自洽 → 拟合 |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/17-effmass-gaas
cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR
sbatch job.slurm
cat fd/effmass.txt
```

单独运行拟合脚本：

```bash
cd fd
python3 ../fit_effmass.py                          # 默认 CBM 与 VBM，k0 = Gamma
python3 ../fit_effmass.py --band 7 --band 8        # 价带顶简并的另外两条带
python3 ../fit_effmass.py --npoints 3              # 只用最靠近 Gamma 的 3 个点，检查非抛物线性
python3 ../fit_effmass.py -h
```

## 预期结果 / 如何检查

- 脚本输出每个方向的 `m*_fit`（全部点拟合）与 `m*_2pt`（只用最小步长的两点差分），两者接近说明步长合适；
- 导带底为正值且三个方向几乎相同（Γ 点导带近似各向同性）；
- 不含 SOC 时价带顶在 Γ 三重简并（脚本会提示），按能量排序的三条带分别对应重空穴、轻空穴，不同方向的重空穴质量不同（翘曲）；
- 实验电子有效质量约 0.067 m0（Vurgaftman 等, J. Appl. Phys. 89, 5815 (2001)）。PBE 严重低估 GaAs 带隙，按 k·p 理论电子有效质量也会被明显低估，这是方法本身的误差而不是脚本错误；要定量比较请改用 HSE06 或 mBJ（见教程 §4.5）并考虑 SOC。

## 注意事项

- EIGENVAL 中的 k 点是倒格子分数坐标，脚本用 POSCAR 的晶格换算成 1/Å（含 2π），所以 POSCAR 必须与计算一致；
- 有效质量来自 meV 量级的能量差，`EDIFF` 要足够小；步长太小数值噪声大，太大则偏离抛物线；
- 用 KPOINTS_OPT（VASP ≥ 6.3）一步完成时，把 `KPOINTS_fd` 改名为 `KPOINTS_OPT` 放在自洽目录中，然后对 `EIGENVAL_OPT` 运行 `fit_effmass.py --eigenval EIGENVAL_OPT`（若你的版本不写 EIGENVAL_OPT，可改用 vaspout.h5/py4vasp 读取 `kpoints_opt` 数据）。
