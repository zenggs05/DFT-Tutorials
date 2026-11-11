# 示例 17-3：含空位的 Si 超胞能带反折叠（easyunfold）

## 目的

- 把 53 原子 Si 空位超胞的能带"反折叠"回原胞布里渊区，得到有效能带结构（EBS / 谱函数）；
- 掌握 easyunfold 的 `generate` → VASP 非自洽（保留 WAVECAR）→ `unfold calculate` → `unfold plot` 流程；
- 认识带隙中的空位缺陷态在 EBS 中表现为几乎无色散、谱权重较弱的平带。

对应教程：[第17章 高级电子结构](../../tutorials/17-advanced-electronic.md)（§5）。

## 文件说明

| 文件 | 说明 |
|---|---|
| `POSCAR_prim` | Si 原胞（a = 5.431 Å，实验值），也是能带路径所对应的晶胞 |
| `KPOINTS_band` | 原胞的 L–Γ–X–U\|K–Γ 路径（Line-mode，Reciprocal，每行最后一个记号是标签，easyunfold 可直接读取） |
| `make_vacancy.py` | 用 ASE 生成 3×3×3 超胞并删去一个原子，写出 `POSCAR_sc`（53 原子）；超胞矩阵为 diag(3,3,3) |
| `INCAR_scf` / `KPOINTS_scf` | 超胞自洽，Γ 中心 3×3×3 |
| `INCAR_band` | `ICHARG = 11`、`NBANDS = 144`、`LWAVE = .TRUE.`、`LORBIT = 11` |
| `job.slurm` | 生成结构与 k 点 → 自洽 → 非自洽 → 反折叠与作图 |
| `POTCAR` | **不提供**：`cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR` |

## 运行步骤

```bash
cd VASP/examples/17-unfold-si-vacancy
pip install easyunfold ase
cat $VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR

python3 make_vacancy.py                                                   # -> POSCAR_sc
easyunfold generate POSCAR_prim POSCAR_sc KPOINTS_band --matrix "3 3 3"   # -> easyunfold.json, KPOINTS_easyunfold
grep -c . KPOINTS_easyunfold                                              # 超胞需要计算的 k 点数（约等于行数）

sbatch job.slurm
```

后处理（作业脚本已包含）：

```bash
easyunfold unfold calculate band/WAVECAR          # 计算谱权重，写入 easyunfold.json
easyunfold unfold plot --eref <E-fermi> --emin -4 --emax 4 --out-file unfold.png
easyunfold unfold plot --help                     # 更多作图选项（色图、DOS 叠加等）
```

## 预期结果 / 如何检查

- `unfold.png` 中能带沿原胞的高对称路径展开，形状与完美 Si 原胞能带（[第5章](../../tutorials/05-band-dos.md)）基本一致，但谱线有展宽、权重小于 1；
- 带隙中出现与空位相关的平坦、弱权重的态（具体能量位置依赖于是否弛豫、超胞大小和泛函）；
- 若整张图几乎没有权重或形状完全不对，多半是超胞矩阵错误或 KPOINTS 没有用 `KPOINTS_easyunfold`。

## 注意事项

- 本例为节省时间使用未弛豫的空位结构；正式研究应先用 `ISIF = 2` 弛豫原子位置（晶胞不变），弛豫后的 CONTCAR 作为 `POSCAR_sc`；
- `WAVECAR` 对 53 原子、数百个 k 点可达数 GB，算完反折叠后可删除，绝不要提交到仓库；
- 用 `vasp_gam` 或 `vasp_ncl` 生成 WAVECAR 时，`unfold calculate` 需要分别加 `--gamma` 或 `--ncl`；
- 杂化泛函不能用 `ICHARG = 11`：用 `easyunfold generate ... --scf-kpoints IBZKPT --nk-per-split N` 生成带零权重 k 点的自洽输入；
- 同样的计算也可以用 vaspkit 281–285（需要 `KPATH.in` 与 `TRANSMAT.in`）或 VaspBandUnfolding 完成。
