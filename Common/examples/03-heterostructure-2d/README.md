# 示例 03-heterostructure-2d：二维范德华异质结的晶格匹配与搭建

## 目的

两种六方二维材料（石墨烯、h-BN、MoS2、WS2）的晶格常数一般不相同，要找到一个公共超胞，使上层材料的应变足够小、总原子数又不太多。本示例用 Zur–McGill（ZSL）思路的简化版完成：

1. 枚举两层的六方超胞矢量 v = i·a1 + j·a2，|v| = a·√(i² − ij + j²)，超胞含 n = i² − ij + j² 个原胞；
2. 当两层的 |v| 接近时即可匹配，应变全部加在上层，扭转角为两个矢量的夹角；
3. 按总原子数排序列出候选，选一个搭成异质结：底层不动，上层旋转并拉伸到底层的面内晶胞，层间距默认 3.4 Å，c 方向加真空。

对应教程：[第3章 结构准备](../../tutorials/03-structure-modeling.md) 第 10 节；vdW 修正与二维材料计算见 [VASP 第8章](../../../VASP/tutorials/08-vdw-2d.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| [build_heterostructure.py](build_heterostructure.py) | `--pair 底层-顶层`（graphene、hBN、MoS2、WS2 任意组合）、`--tol` 最大应变（0.03 = 3%）、`--natoms-max` 原子数上限、`--nmax` 枚举范围、`--list` 只列候选、`--pick` 选第几个候选、`--dist` 初始层间距、`--vacuum` c 方向总真空 |

输出 `POSCAR_<底层>_<顶层>`（元素已排序，分数坐标）。

## 运行步骤

```bash
conda activate dft
pip install ase numpy

python build_heterostructure.py --pair graphene-MoS2 --list     # candidate table only
python build_heterostructure.py --pair graphene-MoS2 --pick 0   # smallest cell within 3 %
python build_heterostructure.py --pair MoS2-WS2                 # 1x1 on 1x1, strain ~0.3 %
python build_heterostructure.py --pair graphene-hBN --tol 0.02 --pick 1
```

## 预期结果 / 如何检查

- 候选表每行给出底层/顶层超胞 `n×(i,j)`、公共边长 L、上层应变（%）、扭转角和总原子数。MoS2/WS2 的晶格常数（3.16/3.15 Å）几乎相同，第一个候选就是 1×1 对 1×1，应变 3.16/3.15 − 1 ≈ +0.32%，6 个原子。
- 搭建完成后屏幕打印面内晶格常数、γ（脚本用 120° 约定的超胞矩阵，应为 120°）、c 和真空厚度；用下面的命令检查原子数和元素顺序：

```bash
sed -n 6,7p POSCAR_graphene_MoS2      # element names and counts
```

- 在 VESTA 中沿 c 轴看：两层都应完整铺满晶胞、没有原子重叠；沿 a 轴看层间距约 3.4 Å。

## 注意事项

- 默认晶格常数是近似实验值；正式计算请先用**同一泛函 + vdW 修正**分别优化两种单层，再用优化值匹配。
- 应变应加在“更软”或“不关心性质”的那一层；如果关心两层的能带对齐，可以把应变平分到两层，或者选应变更小、原子更多的候选。一般要求应变 < 1～2%，并注明所加应变。
- 层间距 3.4 Å 只是初猜，需要用 DFT-D3 或 vdW-DF 弛豫；还要比较不同堆垛（如 AA/AB、扭转 0°/60°）的能量，脚本按 60° 周期去重，对 MoS2 这类三重对称层，0° 与 60° 实际上对应不同堆垛。
- 更完整的界面搜索可用 pymatgen 的 `ZSLGenerator` / `CoherentInterfaceBuilder`（`pymatgen.analysis.interfaces`），它能处理非六方晶格和体相材料之间的界面。
- 二维材料计算时 c 方向不要用 k 点（KPOINTS 第三个数取 1），真空 ≥ 15 Å，必要时加偶极修正。
- 需要弯曲的碳片段（球冠）时，可用仓库作者的 [Scripts/modeling/fabricate_spherical_carbon.py](../../../Scripts/modeling/fabricate_spherical_carbon.py)（说明见 [Scripts/README.md](../../../Scripts/README.md) 第 5.2 节）。

返回：[第3章 结构准备](../../tutorials/03-structure-modeling.md) | [Common 目录](../../README.md)
