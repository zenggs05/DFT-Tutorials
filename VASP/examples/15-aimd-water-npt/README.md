# 算例 15-2：液态水的 NpT 分子动力学（Langevin 恒温 + Parrinello–Rahman 恒压）

## 目的

演示 VASP 中 NpT 系综的写法：`IBRION = 0`、`MDALGO = 3`（Langevin）、`ISIF = 3`、`LANGEVIN_GAMMA`（每个元素一个值）、`LANGEVIN_GAMMA_L`、`PMASS`、`PSTRESS`，并用 `ICONST` 文件把盒子形状限制为立方、只允许体积变化（液体在完全自由的 NpT 中容易发生不可逆的盒子剪切变形）。含 H 体系的时间步长、vdW 修正和 Pulay 应力问题也在这里集中体现。

对应教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `build_water_box.py` | 生成 32 个水分子的立方盒子（fcc 格点 + 随机取向，密度 0.997 g/cm³），只依赖 numpy |
| `POSCAR` | 用上面的脚本（默认参数，seed = 2026）生成的初始结构，O 在前、H 在后，盒子边长 9.865 Å |
| `KPOINTS` | Γ 点 → `vasp_gam` |
| `INCAR` | PBE-D3(BJ)，`ENCUT = 550`，`POTIM = 0.5` fs，300 K，`PSTRESS = 0` |
| `ICONST` | 约束三个晶格夹角与三个边长之比（保持立方），体积可变 |
| `job.slurm` | 拼接 POTCAR（O、H）、运行，并从 OUTCAR 提取密度随时间的变化 |

## 运行步骤

```bash
cd VASP/examples/15-aimd-water-npt
export VASP_PP_PATH=/path/to/your/paw_library
python3 build_water_box.py          # 可选：重新生成 POSCAR（--cells 3 得到 108 个分子）
sbatch job.slurm
```

运行后：

```bash
grep "T=" OSZICAR | tail -3                       # 温度
grep "volume of cell" OUTCAR | tail -3            # 体积
tail density.dat                                  # 步数 与 密度 (g/cm^3)
python3 ../15-aimd-si-liquid-nvt/analyze_md.py --potim 0.5 --skip 2000 --pair O-O --species O
```

## 预期结果 / 如何检查

- 平衡后密度应在某个值附近涨落（32 分子体系涨落约 ±0.03 g/cm³）。泛函强烈影响结果：纯 PBE 水明显过度结构化、密度偏低，加入 D3 修正后通常更接近 1 g/cm³（文献结论，具体数值依赖设置）。
- O–O 的 g(r) 第一峰应在约 2.7–2.8 Å（实验约 2.8 Å，近似值）；`--pair O-H` 可看到约 1.0 Å 的分子内峰和约 1.8 Å 的氢键峰。
- 如果体积单调漂移不停，先检查 ENCUT（Pulay 应力）和 `PMASS`，再延长模拟。

## 注意事项

- **时间步长**：含 H 时 `POTIM` 取 0.5 fs（官方建议含 H 体系不超过约 0.7 fs）。也可以把 H 的质量改成氘（INCAR 中 `POMASS = 16.000 2.000`，顺序与 POSCAR 元素一致）从而把步长放到约 1 fs，但这样模拟的是重水的动力学。
- **Langevin 摩擦系数**：`LANGEVIN_GAMMA = 10` ps⁻¹ 有利于快速平衡，但会人为减慢扩散；需要扩散系数等动力学量时，平衡后改用较小的摩擦（如 1 ps⁻¹）或换成 Nosé–Hoover / CSVR 恒温器的 NVT 运行。
- **计算量**：96 原子、550 eV、每步约 10 个电子步，几 ps 的 AIMD 需要数十个节点小时；这正是第 15 章后半部分引入机器学习力场的原因。
- `ICONST` 的替代方案：对正交盒子也可以设 `LATTICE_CONSTRAINTS = .FALSE. .FALSE. .TRUE.`（只让 c 轴伸缩，像活塞一样），见官方 LATTICE_CONSTRAINTS 页面。
- VASP 需要以 `-Dtbdyn` 编译（所有 5.4.4 之后的官方 makefile.include 模板默认包含）才能使用 `MDALGO`。

返回教程：[第15章 从头算分子动力学与机器学习力场](../../tutorials/15-md-mlff.md)
