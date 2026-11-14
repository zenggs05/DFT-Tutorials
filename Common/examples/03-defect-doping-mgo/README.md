# 示例 03-defect-doping-mgo：MgO 超胞中的空位、替位掺杂、间隙原子与双掺杂构型

## 目的

以岩盐结构 MgO 为例，演示点缺陷建模的四种基本操作，以及如何用对称性减少需要计算的构型数：

1. 由 8 原子惯用胞扩成 n×n×n 超胞，并打印缺陷与其周期像之间的距离（超胞尺寸是否够大）；
2. 用 spglib 找出对称不等价位点（Wyckoff 字母和多重度），决定“删哪个原子、换哪个原子”；
3. 生成氧空位 V_O、Al 替位 Mg（Al_Mg）、四面体间隙 Li（Li_i）；
4. 两个 Al_Mg 的配对：按 Al–Al 距离分壳层，再用 pymatgen `StructureMatcher` 去掉对称等价构型。

对应教程：[第3章 结构准备](../../tutorials/03-structure-modeling.md) 第 6、9 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| [make_defects.py](make_defects.py) | `--a` 晶格常数（默认 4.212 Å，室温实验值附近）、`--n` 超胞倍数（默认 2，即 64 原子）、`--max-pairs` 列出的 Al–Al 距离壳层数 |

运行后在当前目录生成（POSCAR 格式，扩展名 `.vasp`，元素已排序）：

| 输出文件 | 内容 |
|---|---|
| `MgO_perfect.vasp` | 完美超胞，形成能的参考 |
| `MgO_VO.vasp` | 删除 1 个 O |
| `MgO_AlMg.vasp` | 1 个 Mg 换成 Al |
| `MgO_Li_i.vasp` | 在 (1/4, 1/4, 1/4)·a 四面体间隙加 1 个 Li |
| `MgO_2AlMg_<k>.vasp` | 对称不等价的 Al–Al 双掺杂构型，按距离从近到远编号 |

## 运行步骤

```bash
conda activate dft
pip install ase spglib pymatgen          # pymatgen optional (used to remove duplicate pairs)

python make_defects.py                   # 2x2x2, 64 atoms, L = 8.42 A
python make_defects.py --n 3 --a 4.25    # 3x3x3, 216 atoms, L = 12.75 A (use YOUR optimized a)
```

## 预期结果 / 如何检查

由几何关系可直接核对：

- n = 2 时超胞 64 原子（Mg32O32），边长 L = 2 × 4.212 = 8.424 Å；空间群 `Fm-3m (225)`，只有两组不等价位置（Mg 与 O，各占一种 Wyckoff 位置）。
- `MgO_VO.vasp` 组成为 Mg32O31，`MgO_AlMg.vasp` 为 AlMg31O32，`MgO_Li_i.vasp` 为 LiMg32O32。
- 四面体间隙到最近邻 Mg/O 的距离为 √3/4 · a ≈ 1.824 Å（a = 4.212 Å），这个距离偏短，弛豫时周围离子会明显外移。
- 最近的 Al–Al 距离壳层为 a/√2 ≈ 2.98 Å，第二壳层为 a ≈ 4.21 Å。

```bash
for f in MgO_*.vasp; do echo "$f: $(sed -n 6p $f) | $(sed -n 7p $f)"; done   # element and count lines
grep -c . MgO_perfect.vasp                                                     # 8 header lines + 64 atoms
```

## 注意事项

- **超胞大小**：缺陷与周期像的距离 L 在生产计算中一般要 ≥ 10 Å，并检查形成能随 n 的收敛；64 原子只适合练习。
- **晶格常数**：用你自己泛函优化的体相 a，否则超胞整体带应变，形成能会偏。
- **带电缺陷**：V_O 可以是 0、+1、+2 价；带电超胞需要在 INCAR 中改 `NELECT`，并做有限尺寸修正（Freysoldt/Kumagai 等）。系统处理带电缺陷建议用 [doped](https://github.com/SMTG-Bham/doped)，它同时生成缺陷、处理电荷态和修正项。
- **自旋**：O 空位、Li 间隙等可能有未成对电子，计算时要开自旋极化（`ISPIN = 2`）。
- **随机掺杂/合金**：大量随机替换请用仓库作者的 [Scripts/modeling/replace_dopant_symbols.py](../../../Scripts/modeling/replace_dopant_symbols.py)（说明见 [Scripts/README.md](../../../Scripts/README.md) 第 5.1 节），或用 SQS 方法（icet、ATAT `mcsqs`）。
- 不装 pymatgen 时，脚本对每个距离壳层只输出一个构型，并打印提示；同一距离下可能还存在不等价构型。

返回：[第3章 结构准备](../../tutorials/03-structure-modeling.md) | [Common 目录](../../README.md)
