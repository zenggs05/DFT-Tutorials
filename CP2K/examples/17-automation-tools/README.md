# 示例 17：用 ASE、cp2k-input-tools、cp2kdata 自动化 CP2K

配套教程：第17章 自动化与工具生态：ASE、cp2k-input-tools、Multiwfn、cp2kdata

## 目的

演示 CP2K 周边 Python/命令行工具的最小可用流程：

1. **ASE**：通过 CP2K 的 shell 模式（`cp2k.psmp -s`）让 ASE 驱动 CP2K 做 H2O 结构优化；输入模板写 UZH 基组和 OT 设置；
2. **旧版 ASE 的 wrapper**：ASE ≤ 3.28 要求命令里含有 `cp2k_shell`，用 wrapper 脚本绕过；
3. **cp2k-input-tools**：`cp2klint` 检查语法，`cp2kgen` 批量生成 CUTOFF 扫描输入，`tocp2k` 把 YAML 转成 CP2K 输入；
4. **cp2kdata + dpdata**：解析 ENERGY_FORCE / MD 输出，并转换成 DeePMD 训练集格式。

## 文件说明

| 文件 | 说明 |
|---|---|
| `ase_h2o_opt.py` | ASE + CP2K 计算器优化孤立 H2O（BFGS），输出 `h2o_opt.xyz`、`h2o_opt.log`；`-h` 查看参数 |
| `ase_cp2k_shell_wrapper.sh` | 为 ASE ≤ 3.28 准备的 wrapper，内部执行 `cp2k.psmp --shell`（可用 `CP2K_SHELL_NP` 指定 MPI 进程数） |
| `h2o_cutoff_template.inp` | cp2kgen 的模板（坐标内联，便于把生成的输入移到不同目录） |
| `cp2kgen_cutoff_scan.sh` | 用 `cp2klint` + `cp2kgen --zip` 生成 `scan_cutoff/c<值>/run.inp`，逐个运行并汇总能量差 |
| `inp_from_yaml.yaml` | cp2k-input-tools 简化 YAML 格式写的 H2O 单点输入 |
| `cp2kdata_to_dpdata.py` | 用 cp2kdata 汇总输出（版本、能量、最大力），可选用 dpdata 写 `deepmd/npy` |
| `job_slurm.sh` | 通用 Slurm 模板：运行一个 `.inp`，或用 `--ase` 运行 ASE 驱动脚本 |

## 运行步骤

```bash
cd CP2K/examples/17-automation-tools
export CP2K_DATA_DIR=/path/to/cp2k/data
python -m pip install "ase>=3.23" "cp2k-input-tools[yaml]" cp2kdata dpdata   # 建议在独立 conda/venv 环境中

# 1) ASE 驱动 CP2K
python -c "import ase; print(ase.__version__)"
export ASE_CP2K_COMMAND="cp2k.psmp -s"                    # ASE >= 3.29（也是其默认值）
# ASE <= 3.28：改用 wrapper
# chmod +x ase_cp2k_shell_wrapper.sh
# export ASE_CP2K_COMMAND="$PWD/ase_cp2k_shell_wrapper.sh"; export CP2K_SHELL_NP=4
python ase_h2o_opt.py

# 2) cp2k-input-tools
cp2klint h2o_cutoff_template.inp
DRYRUN=1 bash cp2kgen_cutoff_scan.sh                     # 只生成输入，先看看
bash cp2kgen_cutoff_scan.sh                              # 生成 + 运行 + 汇总
tocp2k --yaml inp_from_yaml.yaml > h2o_from_yaml.inp
cp2k.psmp -c h2o_from_yaml.inp                           # 以真正的 CP2K 语法检查为准
cp2k.psmp -i h2o_from_yaml.inp -o h2o_from_yaml.out

# 3) 解析输出，转换成 DeePMD 数据集
python cp2kdata_to_dpdata.py ef h2o_from_yaml.out scan_cutoff/c*/run.out
python cp2kdata_to_dpdata.py ef scan_cutoff/c*/run.out --deepmd dp_h2o   # 演示格式转换；真实训练集应来自不同结构
```

## 预期结果 / 如何检查

- `ase_h2o_opt.py` 打印初始/最终能量（eV）、O–H 键长与 H–O–H 键角；BFGS 通常在 10 步以内达到 `fmax=0.01 eV/Å`。PBE 的 O–H 键长应约为 0.97 Å、键角约 104°（PBE 的典型值，仅作量级参考）；
- 运行目录中出现 `h2o_ase.inp`（ASE 生成的完整输入）和 `h2o_ase.out`，可以 `cp2k.psmp -c h2o_ase.inp` 复核；
- `scan_cutoff/summary.dat` 列出每个 CUTOFF 的总能量及其相对最大 CUTOFF 的差值（meV）；随 CUTOFF 增大应单调收敛到 0 附近；
- `tocp2k` 生成的输入中段的顺序与 YAML 不同（例如 `&GLOBAL` 在文件末尾），这不影响 CP2K；
- `cp2kdata_to_dpdata.py` 打印每个输出的 CP2K 版本、运行类型、能量与最大原子力。

本仓库没有附带参考输出。

## 注意事项

- ASE 的 CP2K 计算器会自动写入 `PROJECT`、`PRINT_LEVEL`、`CUTOFF`、`MAX_SCF`、`&XC`、`&CELL`、`&COORD` 等；模板里**不能**再写这些关键字，否则 ASE 报 `Keyword ... already present`。要用模板里的 `&KIND`，就把 `basis_set` 和 `pseudo_potential` 设为 `None`；
- cp2k-input-tools 0.9.3 内置的是 CP2K 2024.1 的输入结构（`cp2k_input.xml`）：2025.x/2026.x 新增的关键字（如 `GFN_TYPE`、2026.2 的 `&DOS/&PDOS`）会被 `cp2klint`/`cp2kgen` 判为错误，已删除的关键字（如 `&CELL_OPT TYPE`）反而能通过。**最终以 `cp2k.psmp -c` 为准**；
- 我们在 2026-11-10 用 cp2k-input-tools 0.9.3 + pint 0.26.2 测试时发现：`[fs]` 单位偶尔被识别成 femtosiemens（结果随 Python 哈希种子变化），`[kcalmol]` 无法换算。遇到这类“单位错误”可去掉方括号单位（使用该关键字的默认单位），或固定 `PYTHONHASHSEED`，并以 `cp2k.psmp -c` 的结论为准；
- `cp2kgen` 输出文件名形如 `<模板名>-cutoff_400-project_name_c400.inp`，预处理指令（`@SET`、`@INCLUDE`）会被展开、注释会丢失；
- cp2kdata 要求输出用 `PRINT_LEVEL MEDIUM` 写出，且单点计算要打印力（`RUN_TYPE ENERGY_FORCE` 加 `&FORCE_EVAL &PRINT &FORCES ON`，本目录的模板已经这样写）；SCF 未收敛的单点会被 dpdata 插件当作空数据跳过；
- 同一结构不同 CUTOFF 的数据只适合演示格式转换。真正的训练集要用同一套 DFT 参数计算许多不同结构（AIMD 抽帧、微扰结构），并且参数在整套数据中保持一致。
