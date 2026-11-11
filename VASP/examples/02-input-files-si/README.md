# 示例 02：四大输入文件模板与 Si 第一个单点能计算

配套教程：[第2章 四大输入文件与输出文件全解析](../../tutorials/02-input-output-files.md)

## 目的

- 提供 INCAR、POSCAR、KPOINTS 各种写法的可直接使用的模板；
- 用 fcc Si 的单点能计算（PREC=Accurate、ENCUT = 1.3×ENMAX、ISMEAR=0、SIGMA=0.05）走通“准备 → 自检 → 运行 → 检查”的全过程；
- 提供按官方推荐列表拼接 POTCAR 的脚本，以及运行前的输入自检脚本。

## 文件说明

| 文件 | 说明 |
|---|---|
| `Si_scf/INCAR` | Si 单点自洽：PREC=Accurate、ENCUT=320、EDIFF=1E-6、ISMEAR=0、SIGMA=0.05、LCHARG=.TRUE.、LORBIT=11 |
| `Si_scf/POSCAR` | 金刚石 Si 原胞（fcc 晶格矢量，a = 5.43 Å，2 个原子，Direct 坐标） |
| `Si_scf/KPOINTS` | Γ 中心 8×8×8 |
| `POSCAR_selective_dynamics_example` | 4 层 Cu(111) 1×1 板层（a = 3.61 Å，真空约 16 Å），底部两层用 `F F F` 固定 |
| `KPOINTS_gamma` | Γ 中心自动网格模板 |
| `KPOINTS_MP` | Monkhorst–Pack 模板（奇数网格，等价于 Γ 中心） |
| `KPOINTS_line` | fcc 高对称路径 Γ–X–W–K–Γ–L–U–W–L–K\|U–X（对应 `Si_scf/POSCAR` 的原胞取法） |
| `KPOINTS_explicit` | 显式 k 点列表格式演示（fcc Chadi–Cohen 双点集，Cartesian 坐标） |
| `INCAR_annotated_template` | 按功能分组、带英文注释的 INCAR 模板（中文解释见教程第 3 节） |
| `potcar_recommended_list.md` | 官方推荐 PAW 赝势后缀速查表（2026-10 起生效，含 ZVAL、ENMAX、GW 推荐） |
| `make_potcar.sh` | 读取 POSCAR 第 6 行，按推荐后缀从 `$VASP_PP_PATH/potpaw_PBE/` 拼接 POTCAR，可用 `Element=Variant` 覆盖 |
| `check_inputs.sh` | 运行前只读自检：文件齐全、Tab/CRLF、POSCAR 与 POTCAR 顺序、ENCUT/ENMAX、MAGMOM 个数、k 点网格陷阱 |

> [!IMPORTANT]
> 本目录不含 POTCAR。`make_potcar.sh` 只是把你本地（已授权）赝势库中的文件拼接起来。

## 运行步骤

```bash
# 0) 赝势库位置：目录下应有 potpaw_PBE/Si/POTCAR 这样的结构
export VASP_PP_PATH=$HOME/vasp/potentials

# 1) 复制算例
EX=$PWD                                  # 在本示例目录中执行
mkdir -p ~/calc && cp -r Si_scf ~/calc/ && cd ~/calc/Si_scf

# 2) 生成 POTCAR 并自检
bash $EX/make_potcar.sh                  # 打印 TITEL、ZVAL、ENMAX 和 1.3 x max(ENMAX)
bash $EX/check_inputs.sh                 # 应显示 "No errors found."

# 3) 运行
export OMP_NUM_THREADS=1
mpirun -np 4 vasp_std > vasp.out 2>&1

# 4) 检查
grep "General timing and accounting" OUTCAR
tail -3 OSZICAR
grep -E "NKPTS|NELECT" OUTCAR
grep "energy(sigma->0)" OUTCAR | tail -1
grep "EENTRO" OUTCAR | tail -1
```

集群上运行请使用 [示例 01](../01-setup-compile/README.md) 中的作业脚本。

## 预期结果 / 如何检查

| 检查项 | 预期 |
|---|---|
| 正常结束 | OUTCAR 末尾有 `General timing and accounting informations for this job` |
| 不可约 k 点 | 29（8×8×8 Γ 中心网格在 fcc 点群下的约化结果；以 OUTCAR 为准） |
| 电子数 | NELECT = 8（2 个 Si × ZVAL 4） |
| 电子步 | 远少于 NELM = 100 |
| 熵项 | `EENTRO` 接近 0（Si 是半导体） |
| 总能 | 约 −10.8 eV / 2 原子（经验量级，仅供参考） |
| 输出文件 | CHGCAR、DOSCAR、EIGENVAL、PROCAR、vasprun.xml（以及 HDF5 版的 vaspout.h5） |

## 注意事项

- 本示例的 ENCUT 与 k 点只是合理起点，是否收敛请用 [示例 03](../03-convergence-si-al/README.md) 测试；
- `KPOINTS_line` 只能用于 `ICHARG = 11` 的非自洽计算或作为 `KPOINTS_OPT`，不能直接用于自洽；
- 在 Windows 上编辑过的文件上传后先 `dos2unix`；
- 所有输入文件均为纯 ASCII、LF 换行；不要把 POTCAR、WAVECAR、CHGCAR、vaspout.h5 提交到公开仓库。
