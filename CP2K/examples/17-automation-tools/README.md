# 示例 17：用 ASE、cp2k-input-tools、Multiwfn、cp2kdata、AiiDA 自动化 CP2K

配套教程：[第17章 自动化与工具生态：ASE、cp2k-input-tools、Multiwfn、cp2kdata](../../tutorials/17-tools-ecosystem.md)

## 目的

本目录给出 CP2K 周边工具的最小可用流程，所有输入都能直接复制到集群上使用：

1. **ASE**：通过 CP2K 的 shell 模式（`cp2k.psmp -s`）驱动 CP2K。先做孤立 H2O 的 BFGS 优化，再做体相 Si 的 E–V 曲线（k 点经输入模板传入）；
2. **旧版 ASE 的 wrapper**：ASE ≤ 3.28 要求命令里含 `cp2k_shell`，用 wrapper 脚本绕过；
3. **cp2k-input-tools**：`cp2klint` 检查语法，`cp2kgen` 批量生成 CUTOFF 扫描输入，`tocp2k` 把 YAML 转成 CP2K 输入；
4. **Multiwfn 生成的输入**：`Si.POSCAR` 用来演示 Multiwfn 的 `cp2k` 功能；生成的输入先用 `inp_migration_check.py` 检查与版本有关的写法，再交给 `cp2k.psmp -c`；
5. **cp2kdata + dpdata**：解析 ENERGY_FORCE / MD 输出，转换成 DeePMD 训练集格式；
6. **AiiDA（aiida-cp2k 2.2.1）**：用 `Cp2kBaseWorkChain` 提交一个单点计算；
7. **Slurm / PBS 作业模板**：运行单个 `.inp`，或者运行 ASE 驱动脚本。

## 文件说明

| 文件 | 说明 |
|---|---|
| `ase_h2o_opt.py` | ASE + CP2K 计算器优化孤立 H2O（BFGS），输出 `h2o_opt.xyz`、`h2o_opt.log`；用 `-h` 查看参数 |
| `ase_si_eos.py` | ASE + CP2K 计算体相 Si 的 E–V 曲线（k 点、对角化、Broyden 混合都写在模板里），用 Birch–Murnaghan 方程拟合出 $a_0$ 和 $B$，输出 `si_eos.dat`、`si_eos.png` |
| `ase_cp2k_shell_wrapper.sh` | 给 ASE ≤ 3.28 用的 wrapper，内部执行 `cp2k.psmp --shell`（MPI 进程数由 `CP2K_SHELL_NP` 指定） |
| `h2o_cutoff_template.inp` | cp2kgen 用的模板（坐标直接写在输入里，生成的输入可以移到别的目录） |
| `cp2kgen_cutoff_scan.sh` | 先 `cp2klint`，再用 `cp2kgen --zip` 生成 `scan_cutoff/c<值>/run.inp`，逐个运行后汇总能量差 |
| `inp_from_yaml.yaml` | 用 cp2k-input-tools 简化 YAML 格式写的 H2O 单点输入 |
| `Si.POSCAR` | Si 原胞（金刚石结构，a = 5.431 Å），供 Multiwfn（文件名须含 `POSCAR`）或 `../../scripts/cp2k_inpgen.py` 使用 |
| `inp_migration_check.py` | 只用标准库的静态检查：2026.2 删除或改写的语法（`CELL_OPT TYPE`、`&PRINT &PDOS`、REFTRAJ 的 `EVAL_ENERGY_FORCES`/`EVAL_FORCES`），以及新语法用在旧版本上（`--target`）、OT 与 k 点同用、PERIODIC 不一致、基组和赝势的 -qN 不一致、HFX `CUTOFF_RADIUS` 超过晶胞宽度的一半 |
| `cp2kdata_to_dpdata.py` | 用 cp2kdata 汇总输出（版本、能量、最大力），可选用 dpdata 写出 `deepmd/npy` |
| `aiida_h2o_energy.py` | AiiDA + aiida-cp2k：用 `Cp2kBaseWorkChain` 计算 H2O 单点，可在前台运行，也可 `--submit` 给 daemon |
| `job_slurm.sh` | Slurm 模板：运行一个 `.inp`，或者用 `--ase` 运行 ASE 驱动脚本 |
| `job_pbs.sh` | PBS Pro/OpenPBS 模板（Torque 的写法见注释）：`qsub -v INP=x.inp`，或者 `-v ASE_SCRIPT=x.py` |

## 运行步骤

```bash
cd CP2K/examples/17-automation-tools
export CP2K_DATA_DIR=/path/to/cp2k/data
# 建议在独立的 conda/venv 环境中安装（这些包都要求 Python >= 3.10）
python -m pip install "ase>=3.23" "cp2k-input-tools[yaml]" cp2kdata dpdata
# CP2K >= 2025.1 的输出要用 GitHub 上的 cp2kdata（PyPI 0.7.3 读不到新的能量行，见教程 7.1 节）
python -m pip install --upgrade "git+https://github.com/robinzyb/cp2kdata.git"

# 1) ASE 驱动 CP2K
python -c "import ase; print(ase.__version__)"
export ASE_CP2K_COMMAND="cp2k.psmp -s"                    # ASE >= 3.29（这也是它的默认值）
# ASE <= 3.28：改用 wrapper
# chmod +x ase_cp2k_shell_wrapper.sh
# export ASE_CP2K_COMMAND="$PWD/ase_cp2k_shell_wrapper.sh"; export CP2K_SHELL_NP=4
python ase_h2o_opt.py
python ase_si_eos.py --kpts 6 --cutoff 500

# 2) cp2k-input-tools
cp2klint h2o_cutoff_template.inp
DRYRUN=1 bash cp2kgen_cutoff_scan.sh                     # 只生成输入，先检查
bash cp2kgen_cutoff_scan.sh                              # 生成、运行、汇总
tocp2k --yaml inp_from_yaml.yaml > h2o_from_yaml.inp
cp2k.psmp -c h2o_from_yaml.inp                           # 以 CP2K 自己的语法检查为准
cp2k.psmp -i h2o_from_yaml.inp -o h2o_from_yaml.out

# 3) Multiwfn 生成输入（交互式，按键序列见教程第 6 节），然后检查
#    Multiwfn Si.POSCAR   ->  cp2k  ->  Si_sp.inp  ->  ...  ->  0
python inp_migration_check.py Si_sp.inp                  # 目标版本默认 2026.2
python inp_migration_check.py --target 2025.2 Si_sp.inp  # 要在旧版本上运行时
cp2k.psmp -c Si_sp.inp

# 4) 解析输出并转换为 DeePMD 数据集
python cp2kdata_to_dpdata.py ef h2o_from_yaml.out scan_cutoff/c*/run.out
python cp2kdata_to_dpdata.py ef scan_cutoff/c*/run.out --deepmd dp_h2o   # 仅演示格式转换

# 5) AiiDA（需要先配置 profile、computer、code，见脚本开头的说明）
python aiida_h2o_energy.py cp2k@mycluster --nprocs 4

# 6) 作业提交
sbatch job_slurm.sh h2o_from_yaml.inp
sbatch job_slurm.sh --ase ase_si_eos.py
qsub -v INP=h2o_from_yaml.inp job_pbs.sh
```

## 预期结果 / 如何检查

- `ase_h2o_opt.py` 会打印初始和最终能量（eV）、O–H 键长和 H–O–H 键角。步数取决于初始结构和 SCF 精度，看 `h2o_opt.log` 最后一行的 fmax 是否低于 0.01 eV/Å。PBE 下 O–H 约 0.97 Å、键角约 104°，这只是 PBE 的典型量级，不是本仓库的计算结果。
- 运行目录里会出现 `h2o_ase.inp`（ASE 生成的完整输入）和 `h2o_ase.out`，可以用 `cp2k.psmp -c h2o_ase.inp` 再检查一遍。
- `ase_si_eos.py` 打印每个晶格常数下的 E/atom，最后给出 $V_0$、$a_0$ 和体模量 $B$。PBE 一般会高估 Si 的晶格常数，文献中 PBE 的 $a_0$ 约为 5.47 Å（如 Haas, Tran, Blaha, *Phys. Rev. B* **79**, 085104 (2009)，<https://doi.org/10.1103/PhysRevB.79.085104>）。DZVP 基组和 6×6×6 k 点得到的值可能与此差约 0.01 Å，需要做收敛测试后再下结论。
- `scan_cutoff/summary.dat` 列出每个 CUTOFF 的总能量，以及相对最大 CUTOFF 的能量差（meV）。CUTOFF 增大时能量差应单调趋近 0。
- `tocp2k` 输出的段落顺序与 YAML 不同（例如 `&GLOBAL` 可能排在文件末尾），这不影响 CP2K 运行。
- `inp_migration_check.py` 每个问题输出一行 `文件:行号: ERROR/WARN/NOTE: 说明`；只要有 ERROR，退出码就是 1。
- `cp2kdata_to_dpdata.py` 打印每个输出文件的 CP2K 版本、运行类型、能量和最大原子力；使用 `--deepmd` 时还会打印每个文件的帧数。
- `aiida_h2o_energy.py` 打印 WorkChain 的 PK、退出状态和能量（Hartree，即原子单位）。`verdi process show <PK>` 可以查看完整的 provenance。

本仓库没有附带参考输出。

## 注意事项

- ASE 的 CP2K 计算器会自动写入 `PROJECT`、`PRINT_LEVEL`、`CUTOFF`、`MAX_SCF`、`&XC`、`&CELL`、`&COORD` 等内容，模板里**不能**再写这些关键字，否则 ASE 会报 `Keyword ... already present in section ...`（`&XC_FUNCTIONAL` 下的泛函子段不查重：模板里写了泛函、又设了 `xc`，两个泛函会叠加而不报错）。例外是把对应参数设为 `None`：例如 `xc=None` 时可以在模板里写完整的 `&XC`（杂化泛函就这样做）。要使用模板里的 `&KIND`，必须把 `basis_set` 和 `pseudo_potential` 设为 `None`。
- cp2k-input-tools 0.9.3 自带的是 CP2K 2024.1 的输入结构（`cp2k_input.xml`）。因此 2024.1 之后新增的写法（如 2025.1 的 `GFN_TYPE`、2026.2 在 `&DOS` 里内嵌的 `&PDOS`/`&CURVE`）会被 `cp2klint`/`cp2kgen` 判为错误，而 2026.2 已经删除的关键字（如 `&CELL_OPT TYPE`）反而能通过。**最终以 `cp2k.psmp -c` 为准**。
- 2026-11-10 用 cp2k-input-tools 0.9.3 + pint 0.26.2 测试时发现两个问题：`[fs]` 单位有时被识别成 femtosiemens（结果随 Python 哈希种子变化）；`[kcalmol]` 无法换算。遇到这类“单位错误”，可以去掉方括号单位，改用该关键字的默认单位，或者固定 `PYTHONHASHSEED`，并以 `cp2k.psmp -c` 的结论为准。
- `cp2kgen` 的输出文件名形如 `<模板名>-cutoff_400-project_name_c400.inp`。生成时预处理指令（`@SET`、`@INCLUDE`）会被展开，注释会丢失。
- cp2kdata 要求输出以 `PRINT_LEVEL MEDIUM` 写出（LOW 会被直接拒绝）。单点计算还要打印力：`RUN_TYPE ENERGY_FORCE`，再加 `&FORCE_EVAL &PRINT &FORCES ON`，本目录的模板都已经这样写了。SCF 未收敛的单点会被 dpdata 插件当作空数据跳过。
- 同一结构在不同 CUTOFF 下算出的数据只适合演示格式转换。真正的训练集要用同一套 DFT 参数计算大量不同的结构（AIMD 抽帧、微扰结构），而且参数在整套数据中必须保持一致。
- Multiwfn 菜单的编号会随版本变化，教程里的按键序列以 sobereva.com/587（2025-03 更新）为准。第一次请交互式地走一遍。生成的输入必须先用 `cp2k.psmp -c` 检查，再提交。
- aiida-cp2k 不会把 `StructureData` 的 pbc 转换成 `PERIODIC`。孤立分子要在参数里同时写 `&CELL PERIODIC NONE` 和 `&POISSON PERIODIC NONE`，`aiida_h2o_energy.py` 已经这样写了。
- 作业脚本里的 `module load`/`source cp2k_env` 要按本机实际情况修改。ASE 通过 stdin/stdout 与 CP2K 通信，MPI 运行时 stdin 必须送到 0 号进程（`mpirun` 默认会这样做）。
- `ase_si_eos.py`、`inp_migration_check.py`、`aiida_h2o_energy.py`、`job_pbs.sh` 只做过 `py_compile`/`bash -n` 静态检查，没有在集群上实际运行过。使用前请先用小体系试跑。
