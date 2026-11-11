# 算例 11-mno-noncollinear：MnO 分子的非共线 + SOC + 约束 DFT

## 目的

- 练习非共线计算坐标行的 13 列写法（自旋/轨道 Euler 角、约束开关、+U 轨道极化开关）；
- 读懂非共线 Mulliken 表中的 theta / phi 列和 "Orbital moments" 段；
- 用 XCrySDen 显示 `nc.xsf`（原子自旋矢量）与 `ncsden.xsf`（实空间自旋矢量场）；
- 用 `scan_theta.sh` 按手册推荐的"小步旋转 + 逐步续算"流程做约束 DFT 角度扫描。

对应教程：[第11章 自旋极化、非共线磁性与交换耦合](../../tutorials/11-spin-noncollinear.md) §5、§6。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Mol_MnO_NC.dat` | MnO 分子（键长 1.65 Å，近似值），Mn6.0-s3p2d1 / Mn_PBE19（15 价电子）、O6.0-s2p2d1 / O_PBE19；nc + SOC + 方向约束（v = 0.5 eV），θ = 45°、φ = 0° |
| `scan_theta.sh` | 把整个分子的自旋从 θ = 0°（垂直键轴）每 15° 转到 90°（平行键轴），每步读取上一步的重启文件，汇总 Utot、Ucs、Uele 到 `scan_theta.dat`；可直接 `sbatch` |

运行前修改 `DATA.PATH`。

## 运行步骤

```bash
# single noncollinear SCF
mpirun -np 4 openmx Mol_MnO_NC.dat -nt 1 > mno_nc.std 2>&1

# constrained angle scan (7 SCF runs)
sbatch scan_theta.sh          # or: bash scan_theta.sh
cat scan_theta.dat
```

## 预期结果 / 如何检查

```bash
grep -A 8 "Mulliken populations" mno_nc.out   # Up Down Sum Diff theta phi
grep -A 8 "Orbital moments" mno_nc.out        # finite only because SOC is on
grep "Ucs\." mno_nc.out                       # constraint energy, should be small
ls mno_nc.nc.xsf mno_nc.ncsden.xsf mno_nc.nco.xsf
```

- 总自旋磁矩约 5 μB（高自旋），主要在 Mn 上；
- 有约束时各原子的 theta / phi 应接近 45° / 0°；
- `Ucs` 应远小于你关心的能量差；
- XCrySDen：File → Open Structure → Open XSF 打开 `mno_nc.nc.xsf`，Display → Forces 显示箭头。注意箭头是 OpenMX 约定下的"自旋"方向，物理磁矩方向与之相反；
- 角度扫描：Mn²⁺ 型 d⁵ 高自旋组态的轨道角动量近乎猝灭，SOC 引起的各向异性能非常小（可能只有 0.1 meV 甚至更小量级），曲线容易被数值噪声淹没——这正好用来体会"能量差越小，收敛标准越严"。具体数值以你的计算为准。

## 注意事项

- `scf.Constraint.NC.Spin` 只能用于 `scf.SpinPolarization nc`；它与自旋 Zeeman 项（`scf.NC.Zeeman.Spin`）互斥。
- 若怀疑约束影响了能量，用 `scf.Constraint.NC.Spin.v` = 1.0 再算关键角度，确认能量差基本不变。
- 4.0 原版在"非共线 + 仅 Γ 点"时 Mulliken 中 HOMO 序号有误，请先打 patch 4.0.1。
- 块内（`<Atoms.SpeciesAndCoordinates` 与结束标记之间）不要写注释行，否则会报 `Format error of the sequential number`。
- 想看更有意义的"相对夹角—能量"曲线，可参考官方 work 目录中的 `Cr2_CNC.dat`（Cr 二聚体，手册图 40），扫描思路见教程 §6.2。

[返回第11章](../../tutorials/11-spin-noncollinear.md)
