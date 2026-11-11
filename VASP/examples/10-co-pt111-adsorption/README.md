# 示例 10-1：CO 在 Pt(111) 上的吸附能、频率与 d 带中心

## 目的

- 完整演示表面吸附计算：清洁 slab → 气相分子 → 不同位点吸附 → 吸附能；
- 用 `IBRION = 5` + Selective dynamics 只计算吸附质的振动频率，得到 ZPE、熵与自由能校正；
- 用 `LORBIT = 11` 的 DOS 计算表面 Pt 原子的 d 带中心；
- 亲自验证 PBE 下的"CO/Pt(111) 之谜"（PBE 倾向 hollow 位，实验为 top 位）。

对应教程：[第10章 吸附、催化与自由能](../../tutorials/10-adsorption-catalysis.md)（§4–§8、§10）。

## 文件说明

| 路径 | 说明 |
|---|---|
| `1_slab/{INCAR,KPOINTS,POSCAR}` | Pt(111) p(2×2) 4 层，底部 2 层固定，a = 3.97 Å（PBE 近似值），真空约 18 Å；偶极校正 |
| `2_co_gas/{INCAR,KPOINTS,POSCAR}` | CO 分子，15.0×15.5×16.0 Å 盒子，Γ 点（用 `vasp_gam`），`ISMEAR = 0` |
| `3_co_top/{INCAR,KPOINTS,POSCAR}` | CO 在 top 位（C 距表面 1.85 Å），1/4 ML |
| `4_co_fcc/{INCAR,KPOINTS,POSCAR}` | CO 在 fcc hollow 位（C 距表面 1.30 Å） |
| `5_freq_co_top/{INCAR,KPOINTS}` | 吸附质频率：`IBRION = 5`、`POTIM = 0.015`、`NFREE = 2`；POSCAR 由脚本生成 |
| `6_freq_co_gas/{INCAR,KPOINTS}` | 气相 CO 频率（所有原子放开） |
| `7_slab_dos/{INCAR,KPOINTS}` | 静态 DOS：`LORBIT = 11`、`ISMEAR = -5`、`NEDOS = 3001`、12×12×1 |
| `make_potcar.sh` | 按各目录 POSCAR 第 6 行元素从 `$VASP_PP_PATH` 拼接 POTCAR |
| `job.slurm` | Slurm 脚本，在计算目录内用 `sbatch ../job.slurm` 提交 |
| `ads_energy.py` | 读取三个 OUTCAR 的 `energy(sigma->0)` 计算吸附能 |
| `set_freq_selective.py` | 把 CONTCAR 改写为只放开指定原子的 POSCAR（频率计算用） |
| `dband_center.py` | 从 DOSCAR 计算指定原子的 d 带中心、带宽、填充度 |
| `POTCAR` | **不提供**（VASP 许可证），用 `make_potcar.sh` 生成 |

所有 INCAR 的电子部分（`PREC`、`ENCUT = 450`、`LASPH`、`LREAL`）完全一致，这是吸附能可比的前提。

## 运行步骤

```bash
cd VASP/examples/10-co-pt111-adsorption
export VASP_PP_PATH=/path/to/your/paw/library      # 例如包含 potpaw_PBE/Pt/POTCAR
bash make_potcar.sh 1_slab 2_co_gas 3_co_top 4_co_fcc
# 赝势子目录不叫 potpaw_PBE 时：POTDIR=potpaw_PBE.64 bash make_potcar.sh ...

# 1) 四个结构优化（可同时提交）
(cd 1_slab   && sbatch ../job.slurm)
(cd 2_co_gas && VASP_EXE=vasp_gam sbatch --export=ALL ../job.slurm)
(cd 3_co_top && sbatch ../job.slurm)
(cd 4_co_fcc && sbatch ../job.slurm)

# 2) 吸附能
python ads_energy.py --slab 1_slab --ads 3_co_top 4_co_fcc --gas 2_co_gas

# 3) 频率：吸附质（只放开 C、O）与气相分子
python set_freq_selective.py 3_co_top/CONTCAR -e C O -o 5_freq_co_top/POSCAR
cp 2_co_gas/CONTCAR 6_freq_co_gas/POSCAR
bash make_potcar.sh 5_freq_co_top 6_freq_co_gas
(cd 5_freq_co_top && sbatch ../job.slurm)
(cd 6_freq_co_gas && VASP_EXE=vasp_gam sbatch --export=ALL ../job.slurm)
grep "THz" 5_freq_co_top/OUTCAR | head -6

# 4) ZPE / TS / G 校正（脚本在示例 10-2 中）
python ../10-her-oer-free-energy/calc_dG.py thermo 5_freq_co_top -T 298.15

# 5) d 带中心
cp 1_slab/CONTCAR 7_slab_dos/POSCAR && bash make_potcar.sh 7_slab_dos
(cd 7_slab_dos && sbatch ../job.slurm)
python dband_center.py 7_slab_dos/DOSCAR -a 13-16          # 表面层 Pt（原子 13-16）
```

## 预期结果 / 如何检查

```bash
grep "reached required accuracy" */OUTCAR                 # 四个优化都应出现
grep "entropy T" 1_slab/OUTCAR | tail -1                  # |T*S| / 原子数 应 < 1 meV
python -c "from ase.io import read; a=read('2_co_gas/CONTCAR'); print(a.get_distance(0,1))"   # C-O ~1.14 A
```

- CO 吸附能为负，数值在 −1.5 到 −2 eV 量级（文献 PBE 结果的近似范围）；
- PBE 下 fcc 位通常比 top 位更稳定（与实验相反，这正是本示例要演示的现象）；
- 吸附 CO 的 C–O 伸缩频率比气相低（气相 PBE 约 2100 cm⁻¹ 量级，吸附后降低），共 6 个模，无明显虚频；
- 气相 CO 只有 1 个真振动模，其余 5 个接近 0；
- 表面层 d 带中心为负值（相对 E_F），具体数值依赖积分窗口。

以上数值仅为量级参考，请以自己的计算为准。

## 注意事项

- 晶格常数 3.97 Å 是 PBE 的近似值。生产计算请先用相同的 ENCUT/POTCAR/泛函优化 Pt 体相（[第4章](../../tutorials/04-relax-eos.md)），再按新的 a0 重建 slab。
- 打开 `IVDW = 12` 时，所有目录（以及体相 a0）都要一起改。
- top 位 CO 若优化后滑到 bridge/hollow，说明该位点不是极小点；需要约束能量时把 C 原子的 Selective dynamics 改成 `F F T`。
- 生产计算中建议把吸附物放在**已优化的** `1_slab/CONTCAR` 上，示例为便于并行提交直接用了理想 slab。
- 仓库不包含 POTCAR、WAVECAR、CHGCAR，也不要把它们提交到公开仓库。
