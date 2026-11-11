# 算例 07-linear-response-u：线性响应法（LDAUTYPE = 3）计算 NiO 中 Ni-3d 的 U

## 目的

用 Cococcioni–de Gironcoli 线性响应方法，从第一性原理求 Hubbard U，而不是凭经验取值：

$$
U = \chi^{-1} - \chi_0^{-1},\qquad
\chi_0 = \frac{\partial N_d^{\rm NSCF}}{\partial V},\quad
\chi = \frac{\partial N_d^{\rm SCF}}{\partial V}
$$

其中 V 是只作用在"被扰动 Ni 原子" d 壳层上的球对称附加势，N_d 是该原子的 d 占据数。非自洽（NSCF，电荷密度固定）响应 χ0 只包含"裸"的轨道重排，自洽（SCF）响应 χ 还包含电子屏蔽；两者倒数之差就是需要用 U 来修正的部分。

对应教程：[第7章](../../tutorials/07-magnetism-u-soc.md) §6。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR` | AFM-II NiO 的 2×2×2 超胞（32 原子）。第 1 个 Ni 单独作为一个"元素"（物种行写作 `Ni Ni O`，个数 `1 15 16`），只对它施加扰动 |
| `KPOINTS` | Γ 中心 2×2×2（快速测试可改 1×1×1） |
| `INCAR.gs` | 无 U 的 DFT 基态：`LORBIT=11` 输出各原子 d 占据，`LMAXMIX=4` 保证 CHGCAR 可用于 `ICHARG=11` |
| `INCAR.lr` | 扰动计算模板：`LDAUTYPE=3`，`LDAUU = LDAUJ = V`（自旋上/下 d 壳层上的附加势），占位符 `VPERT` 由脚本替换 |
| `run_lr.sh` | 三个阶段：`gs`（基态）、`response`（每个 V 各跑一个 NSCF 和一个 SCF）、`fit`（拟合） |
| `fit_u.py` | 从 OUTCAR 的 `total charge` 表读取原子 1 的 d 占据，线性拟合得到 χ0、χ 和 U；`--plot` 输出 `Nd_vs_V.png` |

POTCAR 必须按 POSCAR 的三个"物种"拼接（Ni 重复一次）：

```bash
cat $VASP_PP_PATH/potpaw_PBE/{Ni_pv,Ni_pv,O}/POTCAR > POTCAR
grep TITEL POTCAR       # expect: Ni_pv, Ni_pv, O
```

## 运行步骤

脚本直接调用 `mpirun`，请在已分配的计算节点上（或交互式作业中）运行；进程数和程序名用环境变量控制：

```bash
export NP=32 VASP_EXE=vasp_std
bash run_lr.sh gs          # ~1 SCF run
bash run_lr.sh response    # 6 values of V x (NSCF + SCF) = 12 runs
bash run_lr.sh fit         # -> chi0, chi, U
python fit_u.py --gs gs/OUTCAR --dirs nscf_* scf_* --plot   # same, plus a figure
```

## 预期结果 / 如何检查

1. 基态 `gs/OUTCAR` 末尾的 `magnetization (x)` 中，16 个 Ni 磁矩正负交替，O 约为 0。
2. 每个 NSCF 目录：`grep "total charge" -A6 nscf_0.10/OUTCAR | tail -4`，原子 1 的 d 占据应随 V 线性变化，其他 Ni 几乎不变。
3. `fit_u.py` 输出中，NSCF 的斜率 |χ0| 明显大于 SCF 的斜率 |χ|（屏蔽使响应变小），两组数据的拟合残差应很小（≲1E-3）；若出现明显非线性，说明 V 取得过大或个别计算没收敛。
4. 对 NiO，VASP 官方教程给出的线性响应 U 在 5–6 eV 量级（取决于超胞大小、k 点和 POTCAR；仅作量级参考）。

示意输出（数值仅供参考）：

```text
chi0 = 0.49 1/eV, chi = 0.13 1/eV
U = 1/chi - 1/chi0 = 5.6 eV
```

## 注意事项

- **超胞大小要收敛**：被扰动原子与其周期镜像之间的距离必须足够大，否则 U 会依赖超胞。至少比较 2×2×2 与更大的超胞（或更密的 k 点）。
- 本算例求的是"在 PBE 基态上"的 U（one-shot）。也可在 DFT+U 基态上重复求 U 直到自洽（self-consistent U），但计算量更大。
- `LDAUTYPE = 3` 下 `LDAUU`/`LDAUJ` **不是** Hubbard U/J，而是附加势的强度（分别作用于自旋上、下）；两者取相同值即对总 d 占据施加扰动。
- 求出的 U 对应的是 VASP 的 PAW 投影子定义的 d 占据，换用其他代码或其他 POTCAR 时不能直接照搬。
- 社区还有更完整的自动化脚本，例如 [sylearn/Ueff_VASP](https://github.com/sylearn/Ueff_VASP)（仅链接，不转载）。官方教程见 VASP Wiki 的 [Calculate U for LSDA+U](https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA%2BU)。

[返回第7章](../../tutorials/07-magnetism-u-soc.md)
