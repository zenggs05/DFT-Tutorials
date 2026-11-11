# 算例 14-3：ZnO 的 HSE06 介电张量与 Born 有效电荷（LCALCEPS 有限电场法）

## 目的

DFPT（`LEPSILON`）不支持杂化泛函。要在 HSE06 水平上得到 ε∞、Born 有效电荷和钳位离子压电张量，需要用**有限电场法** `LCALCEPS = .TRUE.`：VASP 先求零场基态，再分别沿 x、y、z 施加小电场（默认 0.01 eV/Å），由极化和受力的变化做有限差分。本算例对纤锌矿 ZnO 演示这一流程，并与 PBE 的结果比较带隙误差对 ε∞ 的影响。

对应教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `POSCAR` | 纤锌矿 ZnO，实验晶格 a = 3.2495 Å、c = 5.2069 Å、u = 0.382（近似值）；c 轴沿 z |
| `KPOINTS` | Γ 中心 6×6×4（有限电场 / PEAD 例程要求包含 Γ 的规则网格） |
| `INCAR.pbe` | 第 1 步：PBE 自洽，生成 WAVECAR 作为杂化计算的初猜 |
| `INCAR.hse_lcalceps` | 第 2 步：`LHFCALC = .TRUE.`、`HFSCREEN = 0.2` + `LCALCEPS = .TRUE.`，`NBANDS = 24`（保留少量空带，便于 VASP 检查电场相对带隙是否过大） |
| `job.slurm` | 依次运行两步并 grep 结果 |

POTCAR：`cat $VASP_PP_PATH/potpaw_PBE/{Zn,O}/POTCAR > POTCAR`（仓库不提供 POTCAR）。

## 运行步骤

```bash
cd VASP/examples/14-dielectric-zno-hse-lcalceps
export VASP_PP_PATH=/path/to/your/paw_library
sbatch job.slurm
```

想得到 PBE 对照值：把 `INCAR.hse_lcalceps` 中 `LHFCALC` 一段注释掉另算一次（或在 PBE 下用 `LEPSILON = .TRUE.`）。

## 预期结果 / 如何检查

```bash
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)" OUTCAR
grep -A 14 "BORN EFFECTIVE CHARGES (including local field effects)" OUTCAR
grep -A 8 "PIEZOELECTRIC TENSOR (including local field effects)" OUTCAR
grep -i "warning" OUTCAR | sort | uniq -c        # 留意电场过大的警告
```

- 纤锌矿为单轴晶体：ε∞ 张量对角，εxx = εyy ≠ εzz；Zn 的 Z* 为正、O 为负，四个原子的 Z* 之和应接近 0。
- 量级参考：实验 ε∞ 约 3.7（近似值）；PBE 严重低估 ZnO 带隙（约 0.7–0.8 eV，实验约 3.4 eV），ε∞ 因而偏大；HSE06 带隙改善（约 2.5 eV 左右，仍偏低），ε∞ 更接近实验（数值仅供参考）。
- 若出现电场相对带隙过大的警告，用 `EFIELD_PEAD = 0.005 0.005 0.005` 减小电场重算并比较。

## 注意事项

- 计算量：杂化泛函 + 4 次自洽（零场 + 三个方向），比 PBE 贵 1–2 个数量级；先用 4×4×3 k 点验证流程。
- GPU（OpenACC）版本在 6.6.1 之前用 `LCALCEPS` 得到的 Born 有效电荷有误，介电张量不受影响（known issue 93）；更早的 GPU 版本 + 杂化泛函 + `LCALCEPS` 还可能崩溃（known issue 14，6.5.0 修复）。
- 想得到离子贡献 ε_ion 和 ε0，可在第 2 步再加 `IBRION = 6`（有限位移，代价再乘以位移数）；PBE 下也可直接用 `IBRION = 8` + `LEPSILON`。
- `LRPA = .TRUE.` 会额外输出"excluding local field effects"一组结果，但 Born 电荷与压电张量只在 `LRPA = .FALSE.`（默认）下有意义。

返回教程：[第14章 介电与光学性质](../../tutorials/14-optics-dielectric.md)
