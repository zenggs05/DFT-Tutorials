# 示例 03-slab-adsorption：Cu(111) 与 TiO2(110) 表面 slab、固定底层与吸附位点

## 目的

- 用 ASE `fcc111` 搭建 p(3×3) Cu(111) 四层 slab：总真空 15 Å，底部两层用 `FixAtoms` 固定，写出带 `Selective dynamics` 的 `POSCAR`；
- 用 pymatgen `SlabGenerator` 枚举金红石 TiO2(110) 的所有终止面，自动挑出化学计量、非极性、上下对称的 slab，按 O–Ti2O2–O 三明治层固定底部；
- 用 pymatgen `AdsorbateSiteFinder` 找出 top/bridge/hollow 吸附位点，再用 ASE `add_adsorbate` 在命名位点（ontop/bridge/fcc/hcp）放 O 原子和 CO 分子，得到吸附计算的初始结构。

对应教程：[第3章 结构准备](../../tutorials/03-structure-modeling.md) 第 7、8 节。计算本身（表面能、功函数、吸附能、偶极修正）见 [VASP 第9章](../../../VASP/tutorials/09-surface-workfunction.md) 和 [VASP 第10章](../../../VASP/tutorials/10-adsorption-catalysis.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| [build_cu111_slab.py](build_cu111_slab.py) | ASE 搭 Cu(111) slab；`--a` 体相晶格常数、`--size` 面内重复、`--layers` 层数、`--vacuum` **总**真空厚度、`--fix` 固定底部层数 |
| [POSCAR_Cu111_slab](POSCAR_Cu111_slab) | 上一脚本默认参数（a = 3.63 Å，3×3×4，真空 15 Å，固定 2 层）的输出：36 个 Cu，18 个 `F F F` |
| [build_slab_TiO2_110.py](build_slab_TiO2_110.py) | pymatgen 搭 TiO2(110)：`--thick` 最小 slab 厚度、`--vacuum` 最小真空、`--repeat 2 1` 面内扩胞、`--fix-trilayers` 固定的三明治层数 |
| [adsorption_sites.py](adsorption_sites.py) | A 部分（pymatgen）：读入 slab，列出对称不等价的吸附位点并为每个位点生成 CO 吸附结构；B 部分（ASE）：O@fcc、O@hcp、CO@top |

## 运行步骤

```bash
conda activate dft
pip install ase pymatgen

python build_cu111_slab.py                       # -> POSCAR_Cu111_slab (same as the shipped file)
python build_cu111_slab.py --layers 5 --fix 2 --vacuum 20 -o POSCAR_Cu111_5L

python build_slab_TiO2_110.py                    # -> POSCAR_TiO2_110, 1x1 cell
python build_slab_TiO2_110.py --thick 16 --repeat 2 1 --fix-trilayers 2

python adsorption_sites.py                       # uses POSCAR_Cu111_slab
python adsorption_sites.py --skip-pymatgen       # ASE part only
```

## 预期结果 / 如何检查

`POSCAR_Cu111_slab` 中的数值可由几何关系核对：面内 |a| = |b| = 3 × 3.63/√2 ≈ 7.700 Å，γ = 60°；层间距 d(111) = 3.63/√3 ≈ 2.096 Å，四层厚 3 × 2.096 ≈ 6.29 Å；c = 6.29 + 15 ≈ 21.29 Å。

```bash
sed -n 6,9p POSCAR_Cu111_slab                 # "Cu", "36", "Selective dynamics", "Direct"
grep -c "F   F   F" POSCAR_Cu111_slab         # 18 fixed atoms (2 bottom layers x 9)
grep -c "T   T   T" POSCAR_Cu111_slab         # 18 free atoms
ls ads_ase_*.vasp ads_pmg_CO_*.vasp           # one file per adsorption structure
```

- `build_slab_TiO2_110.py` 会先列出全部候选终止面，每行给出 `symmetric=` 和 `polar=`；最终选中的 slab 必须是 `TiO2`、`symmetric=True`、`polar=False`。检查打印的 `vacuum gap` 不小于 `--vacuum`，固定原子数等于底部三明治层的原子数。
- `adsorption_sites.py` 对 fcc(111) 面通常给出 ontop、bridge 各 1 个，hollow 2 个（fcc 和 hcp，由次表面原子区分）对称不等价位点（具体个数受 `symm_reduce`、`near_reduce` 阈值影响）；生成的吸附结构中，吸附分子是 `T T T`，原 slab 的 `F F F` 保持不变。
- 用 VESTA 或 `ase gui` 打开文件，确认吸附物在表面上方（而不是穿过周期边界出现在 slab 底部），C–O 键长约 1.15 Å。

## 注意事项

- ASE 的 `vacuum` 参数加在 slab **两侧**，总真空 = 2 × vacuum；`build_cu111_slab.py` 已经把输入的总真空除以 2，自己写代码时别忘了。
- 晶格常数必须用**同一泛函、ENCUT、POTCAR** 优化得到的体相值，否则 slab 带有人为应变；脚本中的 3.63 Å（Cu，PBE 量级）和 Burdett 等人的 TiO2 实验值只是默认值。
- 吸附高度（O 约 1.2 Å、CO 的 C 约 1.85 Å）只是初猜，所有结构都要弛豫；吸附能 E_ads = E(slab+ads) − E(slab) − E(ads) 中三项要用相同晶胞和参数。
- slab 上下不对称或有吸附物时会产生偶极，VASP 用 `LDIPOL = .TRUE.`、`IDIPOL = 3`，CP2K 用 `&DFT` 下的 `SURFACE_DIPOLE_CORRECTION`，OpenMX 没有偶极修正关键词，要用 ESM 方法，而 ESM 要求表面法向沿 **a 轴**：本目录脚本生成的 slab 法向沿 c，交给 OpenMX 前须轮换晶格矢量和坐标（方法见[第3章](../../tutorials/03-structure-modeling.md) 7.4 节，计算见 [OpenMX 第16章](../../../OpenMX/tutorials/16-surface-esm-workfunction.md)）。
- 要转成 CP2K/OpenMX 格式，可用 [03-fetch-convert-si/cif_to_all.py](../03-fetch-convert-si/cif_to_all.py) 或板块通用脚本 [convert_any.py](../../scripts/convert_any.py)（用 `python ../../scripts/convert_any.py -h` 查看用法）；注意 `Selective dynamics` 的固定信息在 CP2K 中要改写成 `&MOTION/&CONSTRAINT/&FIXED_ATOMS`，在 OpenMX 中写在 `<MD.Fixed.XYZ` 块里。

返回：[第3章 结构准备](../../tutorials/03-structure-modeling.md) | [Common 目录](../../README.md)
