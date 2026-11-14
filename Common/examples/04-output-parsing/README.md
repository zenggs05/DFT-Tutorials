# 04-output-parsing：跨程序输出解析、作业状态扫描、数据归档与 MLIP 数据准备

## 目的

- `parse_all_outputs.py`：一个入口读取 VASP / Gaussian / CP2K / OpenMX 输出，统一得到能量（eV）、力和末结构（ASE `Atoms`），可汇总为 extxyz；
- `job_status_scanner.py`：只用 Python 标准库，批量扫描子目录，按关键字判断“正常结束 / 收敛 / 报错 / 仍在运行或被杀”；
- `ase_db_collect.py`：把整个项目中已完成的计算归档到 `ase.db`，导出 CSV（pandas）与 extxyz；
- `dpdata_convert.py`：用 dpdata 把 OUTCAR / CP2K / Gaussian 输出转换为 DeePMD `npy` 与 extxyz，供机器学习势训练。

## 文件说明

| 文件 | 依赖 |
|---|---|
| `parse_all_outputs.py` | ase；可选 cclib 1.8.1（Gaussian 频率、HOMO/LUMO） |
| `job_status_scanner.py` | 仅标准库 |
| `ase_db_collect.py` | ase；`--csv` 需要 pandas |
| `dpdata_convert.py` | dpdata（Python ≥ 3.10）、numpy |
| `mock_runs/` | **人为构造的模拟文件**，只包含状态关键字，用来测试扫描脚本，不是真实程序输出 |

## 运行步骤

```bash
# 1) 先在模拟目录上测试扫描器
python job_status_scanner.py mock_runs

# 2) 扫描真实项目，并导出 CSV
python job_status_scanner.py ~/projects/TiO2 --csv status.csv

# 3) 统一解析若干输出
python parse_all_outputs.py relax/OUTCAR g16/h2o_freq.log cp2k/h2o.out omx/si.out --extxyz final.xyz

# 4) 归档到数据库并查询
python ase_db_collect.py ~/projects/TiO2 --db tio2.db --csv tio2.csv
ase db tio2.db code=vasp -c id,formula,energy,max_force,dir

# 5) 生成 MLIP 训练集（每 5 帧取 1 帧，10% 作测试集）
python dpdata_convert.py aimd_*/OUTCAR --fmt vasp/outcar --stride 5 --test-frac 0.1 --out dataset
```

## 预期结果 / 如何检查

`python job_status_scanner.py mock_runs` 的预期输出（示意输出，按 `mock_runs/` 中的模拟文件推算，数值仅供参考；Windows 下路径分隔符为 `\`）：

```text
path                   code     status                           energy(eV) steps  note
cp2k_scf_fail/h2o.out  CP2K     ERROR                                              SCF run NOT converged;
g16_error/bad.log      Gaussian ERROR                                              failed in l502
g16_ok/h2o.log         Gaussian FINISHED                       -2076.500884     2  opt converged;
openmx_ok/si.out       OpenMX   FINISHED                        -217.691090
vasp_relax_ok          VASP     CONVERGED                      -10.85000000     2
vasp_zbrent            VASP     ERROR                          -50.10000000        ZBRENT (ionic step bracketing failed);
```

（能量数值来自模拟文件，没有物理意义。）

## 注意事项

- VASP 的 ZBRENT、EDDDAV 等报错多数写在**标准输出**（`vasp.out`、`slurm-*.out`）而不是 OUTCAR，扫描器两处都会查。
- CP2K 2025.1 起能量行单位写成 `[hartree]`，此前为 `[a.u.]:`，正则两种都匹配。
- OpenMX 没有统一的“SCF 未收敛”报错行：扫描器读取 `.out` 中回显的 `scf.maxIter`（缺省按默认值 40），最后一个 `SCF history at MD=` 块的 `SCF=` 行数达到该值就标记为未收敛；正常结束看 `.out` 的 `Elapsed.Time.` 或标准输出的 `The calculation was normally finished.`。
- VASP 的 `WARNING: Sub-Space-Matrix is not hermitian` 只是警告，扫描器写进 note，不把状态改成 ERROR；真正的错误是 `Error EDDDAV`、`ZBRENT: fatal error` 等。
- cclib 2.0a1 解析 Gaussian 16 输出存在已知问题，生产环境用 `pip install cclib==1.8.1`。
- ASE 从 OUTCAR 读取的 `energy` 是 `energy(sigma->0)`，`free_energy` 是 `TOTEN`；pymatgen `Vasprun.final_energy` 也是 sigma→0 能量。
- 公开分享数据时不要包含 POTCAR、`vaspout.h5` 中的赝势数据等受许可保护的内容。

返回教程：[第4章 ASE 与 pymatgen 实战](../../tutorials/04-ase-pymatgen.md)
