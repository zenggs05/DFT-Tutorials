# 示例 18：后处理与自动化工具箱（py4vasp / pymatgen / ASE / custodian / atomate2 / vaspkit）

## 目的

- 提供可直接复用的后处理脚本和笔记本：py4vasp 读 vaspout.h5、pymatgen 解析 vasprun.xml/OUTCAR、ASE 驱动 VASP；
- 提供自动化模板：Slurm 作业数组、custodian 自动纠错、atomate2 工作流；
- 提供速查资料：VASPKIT 命令行速查、常见报错对照表、VASP 项目的 `.gitignore` 模板。

对应教程：[第18章 后处理生态、自动化与报错排查](../../tutorials/18-postprocess-troubleshooting.md)。

## 文件说明

| 文件 | 说明 | 依赖 |
|---|---|---|
| `py4vasp_quickstart.ipynb` | 结构、能量、能带、DOS、带隙、STM、导出图片与 CSV | VASP 带 HDF5 编译；`pip install py4vasp kaleido` |
| `pymatgen_parse_vasprun.py` | 收敛性、能量、带隙、磁矩、耗时；`--band` 画能带、`--dos` 画分元素 DOS、`--json` 输出摘要 | `pip install pymatgen` |
| `ase_vasp_calculator_example.py` | 用 ASE 的 `Vasp` 计算器做 Si 的 E–V 扫描并拟合；`--dry-run` 只写输入；`--neb-demo` 演示 IDPP 插值 | `pip install ase`；`VASP_PP_PATH`、`ASE_VASP_COMMAND` |
| `custodian_run.py` | 带错误处理器（ZBRENT、BRMIX、EDDDAV 等）运行 VASP；`--mode double` 双重弛豫 | `pip install custodian pymatgen` |
| `atomate2_relax_band.py` | 双重弛豫 + 均匀网格与线模式能带的 atomate2 工作流，jobflow 本地执行 | `pip install atomate2`；`~/.atomate2.yaml` |
| `array_job.slurm` | Slurm 作业数组：每个子目录一个 VASP 计算，已完成的自动跳过 | Slurm |
| `vaspkit_cli_cheatsheet.md` | VASPKIT 1.5.x 配置、运行方式、功能编号速查 | VASPKIT |
| `error_table.md` | 常见报错 → 原因 → 修复 对照表 | — |
| `.gitignore_vasp` | VASP 项目的 `.gitignore` 模板（POTCAR、WAVECAR、CHGCAR、vaspout.h5 等） | — |

## 运行步骤

```bash
cd VASP/examples/18-postprocess-toolkit

# 1) pymatgen：解析任意已完成的计算目录
python pymatgen_parse_vasprun.py ../16-elf-si --json si_summary.json
python pymatgen_parse_vasprun.py /path/to/band_run --band          # 线模式能带
python pymatgen_parse_vasprun.py /path/to/dos_run --dos            # LORBIT = 11 的 DOS

# 2) py4vasp：在计算目录中启动 Jupyter 并打开笔记本
jupyter lab py4vasp_quickstart.ipynb

# 3) ASE：只写输入（需要 VASP_PP_PATH），或真正运行
export VASP_PP_PATH=/path/to/potentials
python ase_vasp_calculator_example.py --dry-run
export ASE_VASP_COMMAND="mpirun -np 8 vasp_std"
python ase_vasp_calculator_example.py
python ase_vasp_calculator_example.py --neb-demo                    # 写 neb/00..06/POSCAR

# 4) custodian：放在作业脚本中替代直接调用 mpirun
python custodian_run.py --cmd "mpirun -np 32 vasp_std"

# 5) 作业数组：每行一个目录
ls -d strain_* > dirs.txt
sbatch --array=1-$(wc -l < dirs.txt)%10 array_job.slurm

# 6) 把 .gitignore 模板加到自己的项目中
cat .gitignore_vasp >> /path/to/your/project/.gitignore
```

## 预期结果 / 如何检查

- `pymatgen_parse_vasprun.py` 打印 `converged_electronic`/`converged_ionic` 均为 True，未收敛时会给出警告；
- `ase_vasp_calculator_example.py --dry-run` 在 `ev_*` 目录中写出 INCAR、KPOINTS、POSCAR、POTCAR（以及 ASE 的 `ase-sort.dat`），不调用 VASP；真正运行后打印 V0、a0、B0 并保存 `si_eos.png`；
- custodian 运行后生成 `custodian.json`，记录每次检测到的错误与所做的 INCAR/KPOINTS 修改；
- 作业数组中每个子目录生成自己的 `vasp.out`，重复提交时已完成的目录会被跳过。

## 注意事项

- 所有脚本都不包含也不下载 POTCAR；ASE、pymatgen、atomate2 需要指向你本地合法获得的赝势库；
- custodian 会自动修改 INCAR/KPOINTS 并备份旧文件（`error.N.tar.gz`），提交前先确认这些自动修改在物理上是可接受的；
- py4vasp 需要 vaspout.h5，旧版本或未开启 HDF5 的 VASP 无法使用，此时改用 pymatgen/vaspkit；
- 笔记本和脚本只用于读取和作图，不要把 vaspout.h5、vasprun.xml 等大文件提交到仓库。
