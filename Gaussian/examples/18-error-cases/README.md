# 示例 18a：故意出错的 Gaussian 输入与修正版本

## 目的

每个 `err_*` 文件只包含**一个**典型错误，配有修正后的版本。运行它们（都只需几秒）可以亲眼看到报错出现在哪个 Link、真正的错误信息在日志的什么位置，从而学会第18章“从日志末尾往上找”的排错方法。

对应教程：[第18章 常见报错排查大全与完整研究流程案例](../../tutorials/18-errors-capstone.md)（第 3–7 节）。

## 文件说明

| 文件 | 故意制造的错误 | 预期报错（Link） | 修正版 |
|---|---|---|---|
| `err_l301_wrong_multiplicity.gjf` | CH₃ 自由基（9 个电子）写成 `0 1` | L301：`The combination of multiplicity 1 and 9 electrons is impossible.` | `err_l301_wrong_multiplicity_fixed.gjf`：`0 2` + `UB3LYP` |
| `err_l101_missing_blank_line.gjf.txt` | 坐标最后一行之后没有空行，文件直接结束 | L101：`End of file in ZSymb.` | `err_l101_fixed.gjf` |
| `err_l1_qperr_route.gjf` | route 中 `Opt=(MaxCycles=100` 括号未闭合 | L1：`QPErr --- A syntax error was detected in the input line.` | 改为 `Opt=(MaxCycles=100) Freq` |
| `err_l301_m062x_gd3bj.gjf` | M06-2X 搭配 `EmpiricalDispersion=GD3BJ`（Gaussian 中没有该组合的参数） | L301：`R6DS8: Unable to choose the S8 parameter` | `err_l301_m062x_gd3_fixed.gjf`：改用 `GD3` |
| `err_genecp_bad_format.gjf` | [PdCl₄]²⁻ 用了 `GenECP`，但只写了基组段，缺少 ECP 段 | L301：读取 ECP 时遇到文件结尾，通常为 `EOF while reading ECP pointer card.` | `err_genecp_fixed.gjf`：补上 `Pd 0 / LANL2DZ` ECP 段 |
| `err_l502_hard_scf_fixed.gjf` | 不是错误文件，而是难收敛 SCF 的“标准处理流程”：高自旋 [FeCl₄]⁻（S = 5/2），小基组 `SCF=XQC` + `Stable=Opt`，再 `Guess=Read` 到大基组 | — | — |
| `run_error_cases.sh` | 依次运行所有 `err_*.gjf`，打印失败的 Link 和其上方的几行真正的错误信息 | | |

为什么 L101 的例子是 `.gjf.txt`：本仓库的检查工具要求所有 `.gjf` 以空行结尾，编辑器也常常会“顺手”补上空行。把这个故意出错的文件存为 `.txt`，可以保证它保持“错误”的状态；`run_error_cases.sh` 运行前会复制成 `.gjf`。

## 运行步骤

```bash
cd Gaussian/examples/18-error-cases
source $g16root/g16/bsd/g16.profile          # or: module load gaussian/16-C.02
bash run_error_cases.sh | tee error_summary.txt

# run a single case by hand
g16 < err_l301_wrong_multiplicity.gjf > err_l301_wrong_multiplicity.log
tail -n 20 err_l301_wrong_multiplicity.log

# the hard-SCF recipe (8 cores, 8 GB as written in the file)
g16 < err_l502_hard_scf_fixed.gjf > err_l502_hard_scf_fixed.log
```

这些作业都很小，可以在交互节点或自己的工作站上运行；请不要在共享的登录节点上运行任何计算。

## 预期结果 / 如何检查

| 检查项 | 命令 | 期望 |
|---|---|---|
| 出错的 Link | `grep "Error termination" err_*.log` | 分别为 `l301.exe`、`l101.exe`、`l1.exe`、`l301.exe`、`l301.exe` |
| 真正的错误信息 | `grep -B 8 "Error termination" err_l301_wrong_multiplicity.log` | 在 `Error termination` 之上几行能看到具体原因 |
| 修正版正常结束 | `tail -n 1 err_*_fixed.log` | `Normal termination of Gaussian 16` |
| CH₃ 双重态自旋污染 | `grep "S\*\*2" err_l301_wrong_multiplicity_fixed.log` | `<S**2>` 应接近 S(S+1) = 0.75 |
| [FeCl₄]⁻ 稳定性 | `grep -i "stable" err_l502_hard_scf_fixed.log` | 第 1 步结束时波函数稳定 |
| [FeCl₄]⁻ 自旋污染 | `grep "S\*\*2" err_l502_hard_scf_fixed.log` | 应接近 S(S+1) = 8.75 |

报错信息的确切文字可能随 Gaussian 修订版本略有差异；本仓库没有附带运行日志，请以你自己的输出为准。仓库脚本 `g16_error_doctor.py` 可以自动识别这些报错并给出中文建议（`python ../../scripts/g16_error_doctor.py -h` 查看用法）。

## 注意事项

- 电子数奇偶检查：电子数 = 各原子序数之和 − 总电荷。电子数为奇数时多重度只能是 2、4、6……；为偶数时只能是 1、3、5……。使用赝势时仍按全部电子计数。
- `GenECP` 输入的三条规则：每个元素块以 `****` 结束；基组段与 ECP 段之间空一行；整个文件以空行结束。分子中出现的每种元素都要在基组段中出现且只出现一次。
- `err_l502_hard_scf_fixed.gjf` 只演示处理流程，不保证 [FeCl₄]⁻ 在默认设置下一定不收敛。真实体系 SCF 不收敛时，先检查电荷、多重度和结构是否合理，再按教程中的顺序尝试 `SCF=XQC`、`Guess=Read`、`SCF=VShift`、`SCF=Fermi` 等方法，不要只把 `MaxCycle` 调大。
- 想检查 route 语法而不真正运行，可以用 Gaussian 自带的 `testrt` 工具（见第18章）。

[← 返回第18章](../../tutorials/18-errors-capstone.md)
