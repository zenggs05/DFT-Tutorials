# 示例 16-1：CO/Pt(111) 的 Bader 电荷分析

## 目的

- 演示正确的 Bader 流程：`LAECHG = .TRUE.` → `chgsum.pl AECCAR0 AECCAR2` → `bader CHGCAR -ref CHGCAR_sum`；
- 用 POTCAR 中的 `ZVAL` 把 Bader 布居换算成净电荷，并对 CO 分子整体求和；
- 做 FFT 网格收敛检查。

对应教程：[第16章 电荷与成键分析](../../tutorials/16-charge-bonding.md)（§4）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | Pt(111) p(2×2)、4 层、CO 顶位吸附；Pt 晶格常数取 3.97 Å（接近 PBE 值），C–Pt 1.85 Å、C–O 1.15 Å 为**未优化的初始结构**；底部两层固定。建议替换为 [第10章](../../tutorials/10-adsorption-catalysis.md) 优化后的 CONTCAR |
| `INCAR` | 单点 + `LAECHG`、`LCHARG`；`NGXF/NGYF/NGZF` 以注释形式给出，需按自己的 OUTCAR 修改 |
| `KPOINTS` | Γ 中心 6×6×1 |
| `job.slurm` | Slurm 脚本：运行 VASP 后调用 `run_bader.sh` |
| `run_bader.sh` | `chgsum.pl` + `bader` + 净电荷换算 |
| `bader_net_charge.py` | 读 ACF.dat、POSCAR、POTCAR，输出每个原子的净电荷与分组求和（仅用 Python 标准库） |
| `POTCAR` | **不提供**，需自行生成（见下） |

## 运行步骤

```bash
cd VASP/examples/16-bader-co-pt111

# 0) （推荐）用第10章优化好的结构
# cp /path/to/chapter10/CO_on_Pt111/CONTCAR POSCAR

# 1) POTCAR：元素顺序 Pt C O，与 POSCAR 第 6 行一致
cat $VASP_PP_PATH/potpaw_PBE/{Pt,C,O}/POTCAR > POTCAR
grep -E "TITEL|ZVAL" POTCAR

# 2) 第一次运行（默认网格），提交作业
sbatch job.slurm
grep "dimension x,y,z NGXF" OUTCAR

# 3) 网格收敛：把 INCAR 中 NGXF/NGYF/NGZF 改为默认值的约 2 倍，在新目录中再算一次
mkdir -p grid2 && cp INCAR KPOINTS POSCAR POTCAR job.slurm run_bader.sh bader_net_charge.py grid2/
# 编辑 grid2/INCAR，取消注释并填写 NGXF/NGYF/NGZF
(cd grid2 && sbatch job.slurm)
```

只想手动后处理时：

```bash
chgsum.pl AECCAR0 AECCAR2              # -> CHGCAR_sum
bader CHGCAR -ref CHGCAR_sum           # -> ACF.dat BCF.dat AVF.dat
python bader_net_charge.py --group 17-18
```

`chgsum.pl` 来自 VTST 脚本，`bader` 从 Henkelman 组网站下载；两者都需要在 `PATH` 中。

## 预期结果 / 如何检查

```bash
tail -3 ACF.dat                         # NUMBER OF ELECTRONS 应等于 NELECT
grep NELECT OUTCAR                      # 标准 Pt/C/O POTCAR 时为 170
cat bader_net_charge.txt
```

- `sum of Bader populations` 与 `NELECT` 相差应小于约 0.01–0.05 e；
- `group 17-18 (CO)` 给出 CO 分子整体的净电荷（负值表示得到电子）。文献中常见的结论是 CO 从 Pt 得到少量电子，但数值依赖覆盖度、泛函与网格，请以自己的收敛结果为准；
- 两套网格下 C、O 的电荷差应小于约 0.02 e，否则继续加密。

## 注意事项

- 不要把 `CHGCAR`、`AECCAR*`、`CHGCAR_sum` 等大文件提交到仓库；
- 用 `Pt_pv` 等不同 POTCAR 时 ZVAL 不同，净电荷的换算由脚本自动处理，但不同 POTCAR 的结果不能直接比较；
- 自旋极化体系 bader 只分析 CHGCAR 第一块（总密度），自旋分辨分析见教程 §3 的提示；
- `ACF.dat` 中的 `CHARGE` 是布居，不是净电荷。
