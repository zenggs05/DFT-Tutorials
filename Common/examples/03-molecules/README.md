# 示例 03-molecules：分子建模（SMILES→3D）与四款软件的分子输入

## 目的

- 用 RDKit 从 SMILES 生成三维结构：ETKDGv3 生成多个构象，MMFF94（无参数时用 UFF）力场预优化，保留能量最低的构象，并把电荷和自旋多重度写进 xyz 注释行；
- 把一个分子 xyz 转成 Gaussian `.gjf`、CP2K 孤立分子 `&SUBSYS` 片段、VASP 立方盒子 `POSCAR` 和 OpenMX `Atoms.*` 块，同时检查电荷与多重度是否和电子数奇偶一致；
- 提供两个手写好的分子：实验气相几何的 H2O（xyz）和 D6h 起始几何的苯（Gaussian opt+freq 输入）。

对应教程：[第3章 结构准备](../../tutorials/03-structure-modeling.md) 第 11 节；Gaussian 输入格式见 [Gaussian 第2章](../../../Gaussian/tutorials/02-input-anatomy.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| [h2o.xyz](h2o.xyz) | H2O，r(O–H) = 0.9572 Å，∠HOH = 104.52°（实验气相几何），单位 Å |
| [benzene.gjf](benzene.gjf) | 苯，C–C 1.397 Å、C–H 1.084 Å 的 D6h 起始结构，`opt freq B3LYP/6-31G(d) EmpiricalDispersion=GD3BJ`，文件以空行结尾 |
| [smiles_to_3d.py](smiles_to_3d.py) | `python smiles_to_3d.py "<SMILES>" -n 构象数 -o 输出.xyz --seed 随机种子` |
| [xyz_to_inputs.py](xyz_to_inputs.py) | `--charge/--mult`（缺省时读 xyz 注释中的 `charge= mult=`，再缺省为 0/1）、`--method --basis --task --extra`（Gaussian 路由行）、`--nproc --mem`、`--padding` 盒子每侧真空 |

## 运行步骤

```bash
conda activate dft
conda install -c conda-forge rdkit        # or: pip install rdkit
pip install ase numpy

python smiles_to_3d.py "c1ccccc1" -o benzene_rdkit.xyz
python smiles_to_3d.py "CC(=O)Nc1ccc(O)cc1" -n 30 -o paracetamol.xyz
python smiles_to_3d.py "[NH4+]" -o nh4.xyz              # charge +1 is read from SMILES

python xyz_to_inputs.py h2o.xyz
python xyz_to_inputs.py paracetamol.xyz --task "opt freq" --nproc 16
python xyz_to_inputs.py nh4.xyz --padding 8
```

`xyz_to_inputs.py h2o.xyz` 生成 `h2o.gjf`、`h2o_cp2k_subsys.inc`、`POSCAR_h2o_box`、`h2o_openmx_atoms.dat`。

## 预期结果 / 如何检查

- H2O 有 10 个电子，多重度 1 合法；屏幕打印 `charge 0, multiplicity 1, 10 electrons`。H2O 最大尺寸约 1.51 Å（两个 H 之间），默认 `--padding 7.5` 时盒子边长约 16.5 Å。
- 故意给错多重度（`python xyz_to_inputs.py h2o.xyz --mult 2`），脚本会在写文件之前退出，提示多重度与电子数矛盾；Gaussian 遇到同样的错误会报 `The combination of multiplicity 2 and 10 electrons is impossible.`。
- `smiles_to_3d.py` 打印构象数、所用力场和能量范围；xyz 第二行形如 `SMILES="c1ccccc1" charge=0 mult=1 MMFF94_E=...`。

```bash
tail -c 50 benzene.gjf | od -c | tail -3     # file must end with "\n\n" (blank line)
grep -A2 "&CELL" h2o_cp2k_subsys.inc         # ABC ... and PERIODIC NONE
sed -n 6,7p POSCAR_h2o_box                   # "H O" / "2 1" (sorted)
```

- 用 Avogadro 或 VESTA 打开生成的 xyz，测量 C–C、C–H、O–H 键长，和常见值（C–C 芳香 1.39～1.40 Å，C–H 1.08～1.09 Å，O–H 0.96～0.97 Å）比较，偏差大说明结构有问题。

## 注意事项

- 力场几何只是 DFT 优化的起点；柔性分子要做构象搜索（更多构象、CREST、或 [aqme](https://github.com/jvalegre/aqme) 这类自动化流程），不要只取一个构象。
- RDKit 的多重度由“自由基电子数 + 1”估计，对过渡金属配合物、双自由基等不可靠，请自己确认电荷和多重度。
- 周期性程序（VASP、CP2K、OpenMX）计算孤立分子时：盒子每侧真空 ≥ 7～8 Å，k 点只用 Γ 点（VASP 用 `vasp_gam`，OpenMX `scf.Kgrid 1 1 1`），带电或极性分子还要处理偶极（VASP `LDIPOL`/`IDIPOL = 4`；CP2K `PERIODIC NONE` 配合 `PSOLVER MT` 或 `WAVELET`，MT 要求盒子约为分子尺寸的 2 倍）。
- Open Babel 也能做格式转换和 SMILES 生成三维结构，例如 `obabel -:"c1ccccc1" --gen3d -O benzene_ob.xyz`、`obabel h2o.xyz -O h2o_ob.gjf`。

返回：[第3章 结构准备](../../tutorials/03-structure-modeling.md) | [Common 目录](../../README.md)
