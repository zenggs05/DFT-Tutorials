# 04-ase-calculators：用 ASE 统一驱动 VASP、Gaussian、CP2K、OpenMX

## 目的

同一个 `Atoms` 对象，换一个计算器（Calculator）就能交给不同的第一性原理程序计算。本目录给出四个最小可用脚本，并展示每个计算器的环境变量、单位约定和常见坑。

## 文件说明

| 文件 | 程序 | 体系 | 要点 |
|---|---|---|---|
| `ase_vasp_relax.py` | VASP 6.x | 体相 Si | `--mode vasp`：VASP 自己优化（IBRION=2, ISIF=3）；`--mode ase`：ASE 的 `BFGS + FrechetCellFilter`，VASP 只做单点 |
| `ase_gaussian_h2o_opt.py` | Gaussian 16 | H2O | `GaussianOptimizer` 调用 Gaussian 内置优化器，然后做频率 |
| `ase_cp2k_water_sp.py` | CP2K 2023.x–2026.x | H2O（非周期） | `cp2k.psmp -s` shell 模式；复杂设置用 `inp` 模板 |
| `ase_openmx_si.py` | OpenMX 3.9/4.0 | 体相 Si | `energy_cutoff`、`convergence` 以 eV 为单位输入 |
| `job_ase_vasp.slurm` | — | — | 在 Slurm 作业中运行 ASE 驱动脚本的模板 |

## 运行步骤

先设置环境变量（建议写进 `~/.bashrc`，见 [第1章](../../tutorials/01-linux-hpc.md)）：

```bash
# VASP：VASP_PP_PATH 指向“包含 potpaw_PBE/ 的上级目录”
export VASP_PP_PATH=$HOME/pp                 # $HOME/pp/potpaw_PBE/Si/POTCAR
export ASE_VASP_COMMAND="mpirun -np 16 vasp_std"
# Gaussian：PREFIX 会被替换为 label
export ASE_GAUSSIAN_COMMAND="g16 < PREFIX.com > PREFIX.log"
export GAUSS_SCRDIR=/tmp/$USER
# CP2K
export CP2K_DATA_DIR=/path/to/cp2k/data
export ASE_CP2K_COMMAND="mpirun -np 4 cp2k.psmp -s"
# OpenMX
export OPENMX_DFT_DATA_PATH=/path/to/openmx4.0/DFT_DATA19
export ASE_OPENMX_COMMAND=openmx
```

不调用程序、只检查输入是否正确（推荐先做这一步）：

```bash
python ase_vasp_relax.py --dry-run            # 生成 si_relax_vasp/{INCAR,KPOINTS,POSCAR,POTCAR}
python ase_gaussian_h2o_opt.py --dry-run      # 打印 g16/h2o_opt.com
python ase_openmx_si.py --dry-run             # 打印 omx/si.dat（需要 DFT_DATA19 中的 Si_PBE19.vps）
```

正式运行（在计算节点或作业脚本内）：

```bash
sbatch job_ase_vasp.slurm                     # VASP（作业脚本内设置 ASE_VASP_COMMAND）
python ase_gaussian_h2o_opt.py
python ase_cp2k_water_sp.py --relax
python ase_openmx_si.py --np 8
```

## 预期结果 / 如何检查

- VASP：`si_relax_vasp/OUTCAR` 中出现 `reached required accuracy`；PBE 下 Si 的晶格常数约 5.47 Å（PBE 一般比实验值 5.431 Å 大约 0.7%）。
- Gaussian：`g16/h2o_opt.log` 末尾为 `Normal termination`，频率任务没有虚频（`grep "imaginary" g16/h2o_freq.log` 无输出）。
- CP2K：标准输出打印能量和力；运行目录中的 `h2o.inp` 是 ASE 合成的完整输入，`h2o.out` 是 CP2K 输出（`grep "ENERGY| Total FORCE_EVAL" h2o.out`）。调试通信可加 `debug=True`。
- OpenMX：`omx/si.out` 中有 `Utot.`（单位 Hartree）和 `Elapsed.Time.` 行；最后一个 `SCF history at MD=` 块中的 `SCF=` 行数应小于 `scf.maxIter`（本例 100），否则说明 SCF 没有收敛。ASE 把标准输出写到 `omx/si.log`。

## 注意事项

- `ASE_VASP_COMMAND` 未设置时会报 `Please set either command in calculator or one of the following environment variables ...`。
- ASE ≥ 3.23 的文件型计算器（Gaussian、OpenMX 等）若没有命令配置，会报 `EnvironmentError: No configuration of openmx`，设置 `ASE_<NAME>_COMMAND` 或 `~/.config/ase/config.ini` 即可。
- ASE 的 OpenMX 接口会把 `System.CurrentDirectory` 写成绝对路径，脚本中显式设为 `./`；高级功能（DFT+U、NEGF、Wannier）请直接手写 `.dat`。
- 不要把 POTCAR 提交到 Git 仓库或公开上传，它受 VASP 许可证保护。

返回教程：[第4章 ASE 与 pymatgen 实战](../../tutorials/04-ase-pymatgen.md)
