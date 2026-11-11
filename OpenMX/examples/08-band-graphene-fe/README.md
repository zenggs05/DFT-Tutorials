# 算例 8-2：石墨烯（二维）与铁磁 bcc Fe（自旋极化）的能带

## 目的

- **石墨烯**：二维六方体系的 k 路径 Γ–M–K–Γ 写法（120° 晶格矢量约定），观察 K 点的狄拉克锥；体会 4.0 默认 `GCenter` k 网格在 `scf.Kgrid` 为 3 的倍数时包含 K 点。
- **bcc Fe**：共线自旋极化能带。`bandgnu13` 会分别生成上自旋 `fe.BANDDAT1` 与下自旋 `fe.BANDDAT2`；用 `<Band.KPath.UnitCell` 以惯用立方胞定义 bcc 高对称点 Γ、H、N、P。

对应教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Graphene_band.dat` | a = 2.46 Å、真空 20 Å，`C6.0-s2p2d1`，18×18×1 k 点，路径 Γ–M–K–Γ；`System.Name graphene` |
| `Fe_bcc_spin_band.dat` | a = 2.87 Å，`Fe5.5H-s3p2d1` + `Fe_PBE19H`（16 个价电子，初始上/下 10.0/6.0），`scf.SpinPolarization on`，14×14×14 k 点，路径 Γ–H–N–Γ–P–H\|P–N；`System.Name fe` |
| `job.slurm` | 依次运行两个输入，并调用 `bandgnu13` 与 `band_parse.py` |

解析和作图脚本复用 [../08-band-cdia-si/band_parse.py](../08-band-cdia-si/band_parse.py)。

## 运行步骤

```bash
cd OpenMX/examples/08-band-graphene-fe
# edit DATA.PATH in both inputs first
mpirun -np 8 openmx Graphene_band.dat    -nt 1 | tee Graphene_band.std
mpirun -np 8 openmx Fe_bcc_spin_band.dat -nt 1 | tee Fe_bcc_spin_band.std

bandgnu13 graphene.Band          # graphene.BANDDAT1
bandgnu13 fe.Band                # fe.BANDDAT1 (up) + fe.BANDDAT2 (down)

python3 ../08-band-cdia-si/band_parse.py graphene.Band --emin -10 --emax 6 --plot graphene_band.png
python3 ../08-band-cdia-si/band_parse.py fe.Band       --emin -8  --emax 6 --plot fe_band.png
grep -i "spin moment" fe.out | tail -n 3
```

## 预期结果 / 如何检查

- 石墨烯：π 与 π* 带在 K 点相接于费米能附近（狄拉克锥），在 Γ–M 上 π 带远低于 E_F。`band_parse.py` 给出的"带隙"应接近 0 且位于 K 点；由于化学势在有限电子温度下确定以及数值误差，可能显示几 meV 的小带隙，或提示少量能带在 K 附近穿过 E_F，两者都表示半金属。
- Fe：上下自旋能带明显劈裂，费米面处主要是 3d 能带；`fe.out` 中的总自旋磁矩应接近实验值约 2.2 μB（PBE 常得到 2.2–2.3 μB 左右，仅作参考）。`band_parse.py` 会提示有能带穿过 E_F（金属）。
- （示意输出，数值仅供参考）

  ```text
  6 segments, 139 k points, 2 spin channel(s), xx bands
  bands crossing E_F (metal or semimetal along this path): spin0/band5, ...
  ```

## 注意事项

- 六方晶格两种常见约定：本例 a₂ = (−a/2, √3a/2, 0)（夹角 120°）时 K = (1/3, 1/3, 0)；若 a₂ = (a/2, √3a/2, 0)（夹角 60°），K 应写为 (2/3, 1/3, 0) 或等价点。写错时路径不经过真正的 K 点，狄拉克锥会"消失"。
- 二维体系沿真空方向只用 1 个 k 点；能带路径的第三个分量保持 0。
- 金属体系的 SCF 可能较慢：本例用 `RMM-DIISK`、较小混合权重和 500 K 电子温度；更多技巧见第 5 章。
- 非共线或 SOC 计算（`scf.SpinPolarization nc`）只生成一个 `BANDDAT1`，本征值个数加倍（见第 12 章）。

返回教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)
