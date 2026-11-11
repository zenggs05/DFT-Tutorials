# 示例 01：OpenMX 安装目录与数据库导览（OpenMX Installation & DFT_DATA19 Tour）

配套教程：[第1章 OpenMX 概述、版本与学习路线](../../tutorials/01-overview.md)

## 目的

第1章是概述性章节，没有需要提交的计算。本示例提供三个只读的小工具，帮助你在动手计算之前熟悉 OpenMX：

1. 看清安装目录（`source/`、`work/`、`DFT_DATA19/`）里有什么、哪些后处理程序已经编译好；
2. 查询 DFT_DATA19 数据库：每个 VPS 的价电子数与泛函类型、每个 PAO 的截断半径与可用轨道；
3. 检查一个 `.dat` 输入文件的物种定义是否自洽（文件是否存在、基组写法、基函数数、初始上/下自旋电子数之和是否等于价电子数等）；
4. 在 Hartree、eV、Ry、Bohr、Å、GV/m、Debye 等 OpenMX 常用单位之间换算。

## 文件说明

| 文件 | 说明 |
|---|---|
| `explore_openmx.sh` | 只读导览：目录大小、版本信息、是否保留补丁包、PAO/VPS 数量、常用 VPS 的价电子数、已编译的可执行文件、`work/` 中的示例输入（含 `System.Name`、求解器、`MD.Type`） |
| `list_dft_data.py` | 子命令 `vps`（列出 VPS：元素、原子序数 `AtomSpecies`、`valence.electron`、`xc.type`）、`pao`（列出 PAO：`radial.cutoff.pao`、`PAO.Lmax`、`PAO.Mul`）、`check-dat`（检查输入文件）；仅用 Python 标准库 |
| `omx_units.py` | 单位换算（CODATA 2018 常数），`--table` 打印速查表 |

## 运行步骤

```bash
cd /path/to/DFT-Tutorials/OpenMX/examples/01-overview-tour

# 1) tour of a source tree or of the Debian package
bash explore_openmx.sh $HOME/software/openmx4.0
bash explore_openmx.sh /opt/openmx4.0

# 2) query the database
export OMX_DATA=$HOME/software/openmx4.0/DFT_DATA19
python list_dft_data.py vps -e H C O Fe Ni          # valence electrons of all H/C/O/Fe/Ni pseudopotentials
python list_dft_data.py pao -e Fe                  # all Fe PAO files with cutoff radius and available l
python list_dft_data.py vps --csv > vps_table.csv  # whole database as CSV

# 3) check an input file (uses its DATA.PATH unless --data is given)
python list_dft_data.py check-dat ../03-first-calc-methane/Methane.dat --data $OMX_DATA

# 4) units
python omx_units.py --table
python omx_units.py -- -8.033540953061 Ha eV
python omx_units.py 3e-4 Ha/Bohr eV/Ang
```

所有命令都可以用 `-h` 查看帮助，例如 `python list_dft_data.py check-dat -h`。

## 预期结果 / 如何检查

- `explore_openmx.sh` 应显示 `DFT_DATA19/PAO/` 与 `DFT_DATA19/VPS/` 都存在，并列出 `H_PBE19 1.0`、`C_PBE19 4.0`、`O_PBE19 6.0`、`Fe_PBE19H 16.0`、`Fe_PBE19S 14.0` 等价电子数（与手册 12.3 节推荐表一致）；可执行文件一栏若显示 `not found`，说明还没有在 `source/` 中执行 `make all`（见第2章）；
- `list_dft_data.py vps -e Fe` 应同时列出 `Fe_PBE19H`（16 个价电子）和 `Fe_PBE19S`（14 个价电子），这正是坐标行上/下自旋电子数之和必须与所选 VPS 对应的原因；
- `list_dft_data.py check-dat` 对第3章的 `Methane.dat` 应报告 33 个基函数、8 个价电子，`RESULT : OK`；如果 `DATA.PATH` 仍是占位符，则会报错提示修改——这与 OpenMX 实际运行时的行为一致；
- `omx_units.py 1 Ha eV` 输出 `27.21138625 eV`，`3e-4 Ha/Bohr` 约为 `0.0154 eV/Ang`。

## 注意事项

- 三个工具都**不修改**任何文件，可以放心在安装目录上运行；
- `PAO.Mul` 是 PAO 文件中每个角动量存放的径向函数个数，`check-dat` 用它判断 `s3p2d2` 这类写法是否超出文件所能提供的范围；推荐的轨道数仍以手册 12.3 节的 Quick/Standard/Precise 表为准；
- `check-dat` 模仿 OpenMX 解析器的规则：关键词不区分大小写，`#` 之后为注释，同名关键词以第一次出现的为准；它只做静态检查，不能代替真正的计算；
- 在 Windows 上编辑过的输入文件，`check-dat` 会提示 BOM、UTF-16、CRLF 或非 ASCII 字符等问题。
