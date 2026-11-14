# 04-ase-emt-basics：不装 DFT 软件也能跑的 ASE 基础练习（EOS / 约束优化 / NEB / ase db）

## 目的

用 ASE 自带的 EMT 经验势在笔记本电脑上几秒内跑通 ASE 的核心流程：

- `ase_eos_cu.py`：Cu 的 E–V 扫描与 `EquationOfState` 拟合（晶格常数、体模量）；
- `ase_relax_constraints.py`：Cu(100) 表面吸附原子，`FixAtoms` 固定底层，对比 `BFGS` 与 `FIRE`；
- `ase_neb_idpp.py`：IDPP 插值 + CI-NEB，显式指定 `method`（ASE ≥ 3.28 默认改为 `improvedtangent`）；
- `ase_db_demo.py`：把结果写入 `ase.db`（SQLite），按条件查询并导出 CSV。

流程跑通后，只需把 `EMT()` 换成 `Vasp(...)`、`CP2K(...)` 等计算器（见 [04-ase-calculators](../04-ase-calculators/README.md)），其余代码不变。

## 文件说明

| 文件 | 说明 |
|---|---|
| `ase_eos_cu.py` | 7 个体积点的 EOS；`--calc vasp` 可改用 VASP（每个体积一个目录） |
| `ase_relax_constraints.py` | 3×3×4 Cu(100) + 吸附原子，固定下两层，输出 `relax_*.traj`、带 Selective dynamics 的 `CONTCAR_bfgs` |
| `ase_neb_idpp.py` | 吸附原子在相邻 hollow 位之间扩散的 CI-NEB，输出 `neb_barrier.png` |
| `ase_db_demo.py` | 6 种 fcc 金属的 EMT 晶格常数写入 `metals.db`，演示 `select`、`data`、CSV 导出 |

## 运行步骤

```bash
conda activate dft            # Python >= 3.10（ASE 3.29 的 Requires-Python），ase >= 3.28
python -c "import ase; print(ase.__version__)"
export MPLBACKEND=Agg         # 集群/无图形界面时避免 matplotlib 报错

python ase_eos_cu.py
python ase_relax_constraints.py
python ase_neb_idpp.py --nimages 5 --method improvedtangent
python ase_db_demo.py
ase db metals.db 'a0>3.9' -c id,metal,a0,B_GPa
ase gui neb_final_path.traj   # 有图形界面时查看 NEB 路径
```

## 预期结果 / 如何检查

以下为参考输出（示意输出，数值仅供参考，由编写者在 ASE 3.29.0、Python 3.13 上运行得到；不同 ASE 版本的最后几位可能不同。EMT 势，不是 DFT 结果，不要拿去和实验比较）：

```text
EOS = birchmurnaghan
V0  = 11.5654 A^3/atom   a0 = 3.5898 A
E0  = -0.007032 eV          B  = 134.3 GPa

bfgs  converged=True  steps= 12  E=8.5607 eV
fire  converged=True  steps= 22  E=8.5606 eV
fixed atoms: 18 of 37

method = improvedtangent
forward barrier = 0.423 eV   reaction energy = 0.000 eV
```

- EOS：`cu_eos.png` 中拟合曲线应穿过 7 个点，最低点在扫描区间内部；若在边缘，扩大扫描范围。
- 约束优化：两种优化器最终能量一致（差别 < 1 meV），BFGS 步数通常更少；`CONTCAR_bfgs` 中固定原子带 `F F F`。
- NEB：初末态等价，反应能应为 0；`neb_barrier.png` 中最高点位于中间 image。
- ase db：`ase db metals.db` 应列出 6 行，`Keys:` 中包含 `metal, a0, B_GPa`。

## 注意事项

- `ase.db` 中 `energy`、`fmax`、`formula`、`natoms`、`id` 等是保留名，不能当作自定义键（会报 `ValueError: Bad key`）。
- NEB 结果依赖 `method`：复现 ASE ≤ 3.27 的旧结果时写 `method='aseneb'`；新 ASE 不写 `method` 会给出警告。
- EMT 只适合 Al、Cu、Ag、Au、Ni、Pd、Pt（及少量 H/C/N/O 的粗略描述），只用来练习流程。

返回教程：[第4章 ASE 与 pymatgen 实战](../../tutorials/04-ase-pymatgen.md)
