# 算例 08-graphite-vdw：石墨层间距与层间结合能——PBE / PBE-D3(BJ) / optB86b-vdW 对比

## 目的

- 直观看到"纯 GGA 描述不了层间范德华作用"：PBE 变胞优化后石墨的 c 轴会远大于实验值（层间几乎不结合）；
- 学会两类 vdW 修正的 INCAR 写法：原子对加和型 DFT-D3(BJ)（`IVDW = 12`）与非局域泛函 optB86b-vdW（`LUSE_VDW`）；
- 计算层间结合能 E_b = E(石墨烯)/2 − E(石墨)/4，并与 RPA 参考值比较。

对应教程：[第8章 范德华修正、二维与层状材料](../../tutorials/08-vdw-2d.md) §5–§6。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | 2H（Bernal AB）石墨，4 原子，实验初始结构 a = 2.464 Å、c = 6.711 Å |
| `KPOINTS` | Γ 中心 15×15×5（六角晶胞必须用 Γ 中心网格） |
| `INCAR_PBE` | 纯 PBE 变胞优化（`ISIF=3`，`ENCUT=520`）——反面教材 |
| `INCAR_D3BJ` | 同上 + `IVDW = 12` |
| `INCAR_optB86b` | optB86b-vdW：`GGA=MK`、`PARAM1=0.1234`、`PARAM2=1.0`、`AGGAC=0.0`、`LUSE_VDW=.TRUE.`、`LASPH=.TRUE.` |
| `make_graphene.py` | 由石墨 CONTCAR 生成相同面内晶格、c = 20 Å 的单层石墨烯 |
| `run_compare.sh` | `relax` → `static` → `collect` 三阶段自动化 |
| `analyze_graphite.py` | 汇总 a、c 和层间结合能（meV/atom 与 meV/Å²） |
| `job.slurm` | 最简 Slurm 脚本 |

POTCAR：

```bash
cat $VASP_PP_PATH/potpaw_PBE/C/POTCAR > POTCAR
```

optB86b-vdW 不需要特殊 POTCAR（用 PBE 的 PAW 即可）。VASP ≤ 6.4.2 需要把官方分发的 `vdw_kernel.bindat` 复制到计算目录（否则会花数小时现场生成核函数）；VASP ≥ 6.4.3 现场生成只需几分钟，可以不复制。

## 运行步骤

```bash
bash run_compare.sh relax     # 3 ISIF=3 relaxations
# wait ...
bash run_compare.sh static    # graphite + graphene single points for each method
# wait ...
bash run_compare.sh collect
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" */relax/OUTCAR
grep -i "IVDW" D3BJ/relax/OUTCAR | head -3          # confirms the D3 settings were read
grep -i "vdw" optB86b/relax/OUTCAR | head -5        # vdW-DF / kernel related lines
sed -n 3,5p */relax/CONTCAR                          # relaxed lattice vectors
```

- **PBE**：c 明显增大（通常远超实验 6.71 Å，甚至一直变大到 NSW 用完），结合能接近 0 甚至为负——这正是需要 vdW 修正的原因。
- **PBE-D3(BJ)** 与 **optB86b-vdW**：c 应回到实验值附近（偏差通常在百分之几以内），a 几乎不变（约 2.46–2.47 Å）；结合能在几十 meV/atom 量级，可与 RPA 参考值约 48 meV/atom（Lebègue et al., PRL 105, 196401 (2010)）比较。

示意输出（数值仅供参考，不是本仓库实际计算结果）：

```text
method        a (A)    c (A)  E_b (meV/atom)  E_b (meV/A^2)
PBE          2.4670   8.9xxx             1.x           0.xx
D3BJ         2.4640   6.7xxx            5x.x          2x.xx
optB86b      2.4650   6.6xxx            6x.x          2x.xx
```

## 注意事项

- 层间结合能对 k 点（尤其石墨烯的面内 k 点）和 smearing 很敏感；石墨烯与石墨必须用相同的面内 k 网格和相同的 ENCUT。
- `ISIF = 3` 必须用较高 ENCUT（本例 520 eV）以抑制 Pulay 应力；优化结束后建议用 CONTCAR 再优化一次。
- 纯 PBE 的变胞优化可能"跑不完"，这是物理上没有层间结合的表现，不是输入错误；可以减小 `NSW` 只看趋势。
- 不要把 `IVDW` 与 `LUSE_VDW` 同时打开（会重复计入色散）。

[返回第8章](../../tutorials/08-vdw-2d.md)
