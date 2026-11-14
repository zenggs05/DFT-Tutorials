# 示例 07-mlip-nep-mace：机器学习势预优化、NEB 初猜与 NEP/MACE/DeePMD 最小训练配置

## 目的

演示机器学习原子间势（machine-learned interatomic potential, MLIP）在 DFT 工作中最常见的两种用法：

1. **用通用势加速 DFT**：用 MACE-MP 这类基础模型（foundation model），或你自己训练的 DeePMD/NEP 模型，先把结构预优化、给 NEB 生成初始路径，最后**一定用 DFT 精修**；
2. **用 DFT 数据训练专用势**：把 VASP/Gaussian 等的输出转成 extxyz，配好 NEP（GPUMD）、MACE 或 DeePMD 的最小训练输入。

本目录只放输入和脚本，不含训练数据和模型文件。

## 文件说明

| 文件 | 说明 |
|---|---|
| [mlip_calc.py](mlip_calc.py) | 公共函数：按名字创建 ASE 计算器（`mace-mp`、`mace-off`、`mace`、`deepmd`、`nep`、`emt`） |
| [mace_ase_preopt.py](mace_ase_preopt.py) | 结构预优化（可选变胞），保留 POSCAR 中的 Selective dynamics 固定，输出 `POSCAR_preopt` |
| [mlip_neb_guess.py](mlip_neb_guess.py) | IDPP 插值 + MLIP CI-NEB，生成 VTST 用的 `00/ … NN/POSCAR` |
| [dft_to_extxyz.py](dft_to_extxyz.py) | vasprun.xml/OUTCAR/Gaussian log → NEP 格式（`energy`、`virial`、`forces`）或 MACE 格式（`REF_energy`、`REF_stress`、`REF_forces`） |
| [nep.in](nep.in) | GPUMD `nep` 最小训练输入（NEP4，Mo–S 二元体系，文档默认值逐项写出） |
| [mace_config.yaml](mace_config.yaml) | `mace_run_train --config` 用的 YAML 配置 |
| [deepmd_input.json](deepmd_input.json) | DeePMD-kit `se_e2_a` 能量模型的最小 `input.json` |
| [job_mlip_gpu.sh](job_mlip_gpu.sh) | Slurm GPU 作业：`TASK=nep / mace / dp` 三选一 |

## 运行步骤

```bash
# 0) 单独的环境（Python >= 3.11）
conda create -n mlip python=3.12 -y && conda activate mlip
pip install ase mace-torch            # MACE；DeePMD: pip install deepmd-kit；NEP 的 ASE 接口: pip install calorine

# 1) 预优化：手搭的结构 -> MACE-MP -> DFT 起点
python mace_ase_preopt.py POSCAR_guess --cell --fmax 0.05 -o POSCAR_preopt
python mace_ase_preopt.py POSCAR_guess --calc emt       # 没装 MLIP 时先用 EMT 测脚本（只适用于少数金属）

# 2) NEB 初猜：DFT 优化好的初末态 -> IDPP -> MLIP CI-NEB -> VTST 目录
python mlip_neb_guess.py IS/CONTCAR FS/CONTCAR --nimages 5 --outdir neb_init

# 3) DFT 数据 -> 训练集
python dft_to_extxyz.py run_*/vasprun.xml --format nep -o train.xyz
python dft_to_extxyz.py test_*/vasprun.xml --format nep -o test.xyz
python dft_to_extxyz.py run_*/vasprun.xml --format mace -o mace_train.xyz

# 4) 训练（GPU 节点）
sbatch --export=ALL,TASK=nep  job_mlip_gpu.sh
sbatch --export=ALL,TASK=mace job_mlip_gpu.sh
```

DeePMD-kit 需要 `deepmd/npy` 格式的数据目录（`data/train`、`data/valid`），可用 dpdata 从 VASP OUTCAR 转换，见 [Common/examples/04-output-parsing/dpdata_convert.py](../04-output-parsing/dpdata_convert.py)。

## 预期结果 / 如何检查

- `mace_ase_preopt.py`：打印起止能量、最大受力、体积变化和最短原子间距（示意输出，数值仅供参考）：

  ```text
  start : E = -55.1234 eV  (-6.8904 eV/atom)  max|F| = 2.310 eV/A
  end   : E = -56.0021 eV  (-7.0003 eV/atom)  max|F| = 0.041 eV/A
  shortest interatomic distance: 2.402 A
  ```

  最短距离明显偏小、能量异常低，说明模型离开了训练覆盖的区域，这个结构不能用。
- `mlip_neb_guess.py`：路径平滑、无原子重叠；用 `ase gui neb_mlip.traj` 看动画。MLIP 能垒只作参考。
- `dft_to_extxyz.py`：用 `ase.io.read('train.xyz', ':')` 读回，抽查帧数、原子数、能量和力。
- NEP 训练：`loss.out` 随代数下降，用 [Scripts/nep/plot_nep_results.py](../../../Scripts/nep/plot_nep_results.py) 画对角图（见 [Scripts/README.md](../../../Scripts/README.md)）。
- MACE 训练：`results/` 下的曲线和 `*.model` 文件；DeePMD：`lcurve.out`。

## 注意事项

- **MLIP 结果不能替代 DFT**：预优化、NEB 初猜之后，必须用 DFT 重新优化/计算；论文中能量一律用 DFT 值，并写明 MLIP 只用于初始结构。
- 训练集中所有结构必须用**同一套 DFT 设置**（泛函、ENCUT、k 点密度、PREC、展宽、赝势、+U 值），不同设置的数据不能混在一起。
- 能量与力要一致：VASP 的力对应自由能 `free energy TOTEN`，所以脚本默认取 force-consistent 能量（`--energy free`）。
- NEP 用单精度训练，参考能量低于 −100 eV/atom 时会损失精度并给出警告；Gaussian 等全电子能量需先减去孤立原子能量。
- MACE-MP 是在 PBE/PBE+U（MPtrj）数据上训练的，与 Materials Project 数据库中经过修正（MP2020 兼容性修正）的能量不能直接比较。
- MACE-OFF、MACE-OMAT 等部分模型采用 ASL 学术许可，商用前请查看模型说明。
- 示例的 `nep.in` 把 `generation` 设得较小，只用于检查流程；生产训练按 GPUMD 文档调大。

返回教程：[第7章 DFT 理论速览、学习路线与资源汇总](../../tutorials/07-learning-path.md)
