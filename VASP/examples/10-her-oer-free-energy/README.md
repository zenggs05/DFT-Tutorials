# 示例 10-2：计算氢电极——HER 的 ΔG_H\* 与 OER/ORR 台阶图

## 目的

- 计算 H 在 Pt(111) fcc 位（1/4 ML）的吸附自由能 ΔG_H\*；
- 演示气相参考（H2、H2O、三重态 O2）的设置；
- 用 `calc_dG.py` 完成三件事：①从频率 OUTCAR 计算 ZPE/TS/G 校正；②HER 的 ΔG_H\*；③OER/ORR 四步自由能、极限电势、过电位和台阶图。

对应教程：[第10章 吸附、催化与自由能](../../tutorials/10-adsorption-catalysis.md)（§8–§9）。

## 文件说明

| 路径 | 说明 |
|---|---|
| `H_fcc/{INCAR,KPOINTS,POSCAR}` | H 在 Pt(111) p(2×2) fcc 位；INCAR 与示例 10-1 的 slab 完全一致 |
| `H2_gas/`、`H2O_gas/` | H2、H2O 分子（15.0×15.5×16.0 Å 盒子，Γ 点） |
| `O2_gas/` | O2 三重态（`ISPIN = 2`、`MAGMOM = 2*1.0`），用于对比"DFT O2"与"水参考 O2" |
| `INCAR_freq_ads` | 吸附物频率模板（IBRION=5，配合 Selective dynamics） |
| `INCAR_freq_gas` | 气相分子频率模板 |
| `calc_dG.py` | CHE 自由能工具：子命令 `thermo`、`her`、`oer`（`python calc_dG.py -h`） |
| `oer_input_demo.json` | `oer` 子命令的输入格式示例——**其中数值是虚构的演示数字，不是 VASP 结果** |

清洁 slab 直接复用 [示例 10-1](../10-co-pt111-adsorption/README.md) 的 `1_slab`（设置相同），POTCAR 用 `../10-co-pt111-adsorption/make_potcar.sh` 生成。

## 运行步骤

```bash
cd VASP/examples/10-her-oer-free-energy
export VASP_PP_PATH=/path/to/your/paw/library
MK=../10-co-pt111-adsorption/make_potcar.sh
JOB=../../10-co-pt111-adsorption/job.slurm
bash $MK H_fcc H2_gas H2O_gas O2_gas

(cd H_fcc   && sbatch $JOB)
for m in H2_gas H2O_gas O2_gas; do (cd $m && VASP_EXE=vasp_gam sbatch --export=ALL $JOB); done

# HER：常数校正 0.24 eV（Norskov 2005）
python calc_dG.py her --slab ../10-co-pt111-adsorption/1_slab --slabH H_fcc --h2 H2_gas

# 自己算频率替代 0.24 eV
mkdir -p freq_H freq_H2
python ../10-co-pt111-adsorption/set_freq_selective.py H_fcc/CONTCAR -e H -o freq_H/POSCAR
cp INCAR_freq_ads freq_H/INCAR;  cp H_fcc/KPOINTS freq_H/
cp H2_gas/CONTCAR freq_H2/POSCAR; cp INCAR_freq_gas freq_H2/INCAR; cp H2_gas/KPOINTS freq_H2/
bash $MK freq_H freq_H2
(cd freq_H && sbatch $JOB); (cd freq_H2 && VASP_EXE=vasp_gam sbatch --export=ALL $JOB)
python calc_dG.py thermo freq_H                                     # -> G_corr(H*)
python calc_dG.py thermo freq_H2 --gas --geometry linear --symmetry 2 --spin 0 --structure freq_H2/POSCAR
python calc_dG.py her --slab ../10-co-pt111-adsorption/1_slab --slabH H_fcc --h2 H2_gas \
       --gcorr-h <G_corr(H*)> --gcorr-h2 <G_corr(H2)>

# OER/ORR：先用演示文件熟悉输出格式
python calc_dG.py oer oer_input_demo.json --plot oer_demo.png
```

做真实的 OER/ORR 时，按第 10 章 §9.4 计算 OH\*、O\*、OOH\* 的最稳定构型及其频率，把 JSON 中的数字换成**计算目录路径**（脚本会读取其中 OUTCAR 的最后一个 `energy(sigma->0)`）或你自己的能量值。

## 预期结果 / 如何检查

- `grep magnetization O2_gas/OUTCAR`、`grep "number of electron" O2_gas/OUTCAR | tail -1`：O2 总磁矩应为 2.0；
- H2 键长约 0.75 Å，H2O 的 O–H 约 0.97 Å、键角约 104°（PBE 近似值）；
- Pt(111) 的 ΔG_H\* 接近 0（文献中常见在 ±0.1 eV 量级，依赖覆盖度和设置，近似）；
- `thermo` 子命令对吸附 H 给出 3 个实频；气相 H2 只取 1 个最高频率（线形分子 3N−5 = 1）；
- 演示 JSON 的输出为 η_OER = 0.420 V、η_ORR = 0.410 V（**虚构数据，仅用于确认脚本工作正常**）。

## 注意事项

- O2 的 DFT 能量有已知误差，CHE 中通常不直接使用：脚本按 $\Delta G_4 = 4.92\,\mathrm{eV}-(\Delta G_1+\Delta G_2+\Delta G_3)$ 处理，O2 目录只用于对比学习。
- `calc_dG.py thermo` 默认把低于 50 cm⁻¹ 的实频抬高到 50 cm⁻¹（与 VASPKIT 501 一致），`--cutoff 0` 关闭；虚频会被忽略并给出警告。
- 液态水的熵：要么用 0.035 bar 气相水，要么用实验液态水的标准熵，不要叠加使用。`--gas --pressure 3500` 可计算 0.035 bar（= 3500 Pa，约为 300 K 时液态水的饱和蒸气压）下的气相水。
- `--gas` 模式依赖 ASE（`pip install ase`），并需要 `--structure` 给出分子结构；其余功能只需要 numpy（作图需要 matplotlib）。
- 电势使用 RHE 标度（`--U`，单位 V vs RHE）；换算到 SHE：$U_\mathrm{RHE}=U_\mathrm{SHE}+0.059\,\mathrm{pH}$。
