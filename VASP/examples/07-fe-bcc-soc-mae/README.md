# 算例 07-fe-bcc-soc-mae：bcc Fe 的自旋轨道耦合与磁晶各向异性能（力定理）

## 目的

演示用"力定理"（force theorem）计算磁晶各向异性能（MAE）的标准流程：

1. 用 `vasp_std` 做一次高精度共线自旋极化自洽，得到 CHGCAR（`LMAXMIX = 4`）；
2. 用 `vasp_ncl` 读入这份电荷密度，打开 `LSORBIT`，在 `ICHARG = 11`（电荷密度固定）下分别把自旋量子化轴 `SAXIS` 设为 [001]、[100]、[111]；
3. MAE = E(方向) − E(参考方向)。

对应教程：[第7章](../../tutorials/07-magnetism-u-soc.md) §8–§9。

> [!WARNING]
> bcc Fe 是立方晶体，**[001] 与 [100] 在对称性上完全等价**，两者的能量差理论上为 0——本例故意保留它，作为"数值噪声"的度量。真正有意义的是 E[111] − E[001]，而它只有约 1 μeV/atom 的量级（由实验各向异性常数 K1 估算），需要极密的 k 网格和极严的收敛标准才能得到可靠的符号和大小（参见 Halilov et al., PRB 57, 9557 (1998)）。本例以演示流程为主；做练习请换成各向异性大得多的体系，例如 L1₀-FePt（MAE 为 meV/f.u. 量级）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | bcc Fe 1 原子原胞，a = 2.83 Å |
| `KPOINTS` | Γ 中心 16×16×16（仅作起点，必须做 k 点收敛） |
| `INCAR_collinear` | 第一步：共线 `ISPIN=2`，`EDIFF=1E-8`，`GGA_COMPAT=.FALSE.`，`LMAXMIX=4`，`NBANDS=16` |
| `INCAR_001` / `INCAR_100` / `INCAR_111` | 第二步：`LSORBIT=.TRUE.`，`ICHARG=11`，`ISYM=-1`，`SAXIS` 分别为 001 / 100 / 111，`MAGMOM = 0 0 3.0`（沿 SAXIS） |
| `run_mae.sh` | `col` → `soc` → `collect` 三个阶段；`soc` 阶段自动把 NBANDS 设为共线计算的 2 倍 |

POTCAR：

```bash
cat $VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR
```

## 运行步骤

在计算节点上（交互式作业或把命令写进作业脚本）：

```bash
export NP=8 VASP_STD=vasp_std VASP_NCL=vasp_ncl
bash run_mae.sh col        # collinear SCF in col/
bash run_mae.sh soc        # three vasp_ncl runs in soc_001/ soc_100/ soc_111/
bash run_mae.sh collect    # energies and MAE in micro-eV per atom
```

## 预期结果 / 如何检查

```bash
grep "NBANDS=" col/OUTCAR soc_001/OUTCAR | head -2      # SOC run must have 2x bands
grep -i "Euler" soc_111/OUTCAR | head -2                 # angles VASP derived from SAXIS
grep -A4 "magnetization (z)" soc_111/OUTCAR | tail -2    # moment along sigma_3 (= SAXIS)
grep -A4 "orbital moment (x)" soc_001/OUTCAR | tail -2   # LORBMOM output (Cartesian frame)
```

- 三个 SOC 计算的自旋磁矩都应约为 2.2 μB，且在 SAXIS 坐标系中只有 z 分量（`magnetization (z)`）不为 0。
- 轨道磁矩（`orbital moment`）为 0.0x μB 量级（Fe 的轨道矩很小）。
- `collect` 输出中 `100` 相对 `001` 的差值就是你当前设置下的数值误差；只有当 |E[111] − E[001]| 明显大于这个误差、并且随 k 点加密不再变号时，结果才可信。

示意输出（数值仅供参考，不代表收敛结果）：

```text
SAXIS          TOTEN (eV) MAE (ueV/atom)
001           -8.2xxxxxxx           0.00
100           -8.2xxxxxxx           0.3
111           -8.2xxxxxxx           1.5
```

## 注意事项

- 一定使用 `vasp_ncl`；用 `vasp_std` 运行 `LSORBIT = .TRUE.` 会直接报错退出（非共线计算需要不带 `-DNGZhalf` 编译的版本）。
- **不要**在 SOC 的 INCAR 里写 `ISPIN = 2`：VASP 6.5.0 会因此报 MAGMOM 个数错误而停止（known issue 26，6.5.1 起 `LNONCOLLINEAR` 优先，`ISPIN` 被忽略）。
- 共线和 SOC 步骤的 ENCUT、PREC、POTCAR、`GGA_COMPAT`、`LASPH`、`LMAXMIX` 必须一致。
- `ISYM = -1` 使所有方向使用同一套完整 k 点，避免因对称性不同导致不同方向的 k 点集合不同，从而引入系统误差。
- k 点收敛：依次尝试 16³、24³、32³ ……，画出 MAE 随 k 网格的变化。原胞只有 1 个原子，k 点并行（`KPAR`）效率很高。

[返回第7章](../../tutorials/07-magnetism-u-soc.md)
