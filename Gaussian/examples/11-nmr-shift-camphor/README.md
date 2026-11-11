# 示例 11-nmr-shift-camphor：TMS 参比与樟脑的 1H/13C 化学位移

配套教程：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)（第 5、6 节）。

## 目的

用同一套计算协议计算参比物 TMS 和樟脑（(1R,4R)-camphor）的 GIAO 屏蔽常数，按 δ = σ(TMS) − σ 换算化学位移，练习等价核平均与和实验值的比较。

计算协议（两个分子完全相同）：

- 几何：B3LYP-D3(BJ)/6-31G(d)，SMD(氯仿)，Opt + Freq；
- NMR：mPW1PW91/6-311+G(2d,p)，GIAO，SMD(氯仿)。

## 文件说明

| 文件 | 内容 |
|---|---|
| `tms_reference_nmr.gjf` | TMS：Opt+Freq → GIAO NMR |
| `camphor_nmr_giao_smd.gjf` | 樟脑：Opt+Freq → GIAO NMR |
| `run.slurm` | 依次运行两个输入的 Slurm 脚本 |

初始结构由 RDKit 生成（樟脑经 CIP 判定为 1R,4R；对 NMR 而言两个对映体结果相同）。

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < tms_reference_nmr.gjf    > tms_reference_nmr.log
g16 < camphor_nmr_giao_smd.gjf > camphor_nmr_giao_smd.log
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" *.log          # 每个文件应为 3
grep "NImag" *.log                          # 均应为 NImag=0

# TMS 参比屏蔽常数（平均 12 个 H 和 4 个 C）
grep "Isotropic" tms_reference_nmr.log | \
  awk '$2=="H"{h+=$5;nh++} $2=="C"{c+=$5;nc++} END{printf "sigma_H=%.4f  sigma_C=%.4f\n", h/nh, c/nc}'

# 樟脑化学位移：把 SH、SC 换成上一条命令的输出
SH=31.80; SC=184.00     # placeholder values, replace with your own TMS results
grep "Isotropic" camphor_nmr_giao_smd.log | \
  awk -v rh=$SH -v rc=$SC '$2=="H"{printf "%3d H %7.2f\n",$1,rh-$5} $2=="C"{printf "%3d C %7.2f\n",$1,rc-$5}'
```

原子编号对应（详见教程第 6.1 节）：

| Gaussian 编号 | 樟脑编号 | 平均方式 |
|---|---|---|
| 7 | C1 | — |
| 9 | C2（C=O） | — |
| 11 | C3 | H26、H27 为 exo/endo，**不平均** |
| 4 | C4 | H18 |
| 5 | C5 | H19、H20 **不平均** |
| 6 | C6 | H21、H22 **不平均** |
| 2 | C7 | — |
| 1、3 | C8、C9 | 各自的 3 个 H 平均（H12–14、H15–17）；C8 与 C9 互不平均 |
| 8 | C10 | H23–25 平均 |

预期：羰基碳 C2 的位移在 200 ppm 以上的低场区，甲基碳和甲基氢位于高场区；具体数值与实验的偏差（MAE）对 13C 通常为数 ppm、对 1H 通常为零点几 ppm 量级（近似范围，取决于方法）。实验值请查阅文献或谱图数据库后自行比较。

也可以使用仓库脚本：`python ../../scripts/nmr_shift.py -h`。

## 注意事项

- TMS 与样品必须使用**完全相同**的几何水平、NMR 水平和溶剂模型；改动任何一项都要两者一起重算。
- 若改用文献中的线性回归参数（如 Cheshire 数据库），协议必须与该条目一致，此时不需要 TMS。
- 樟脑是刚性分子，只有一个构象；柔性分子需要先做构象搜索并按自由能加权平均。

返回：[第11章 NMR 化学位移与其他分子性质](../../tutorials/11-nmr-props.md)
