# 示例 10-ecd-methyloxirane：(R)-甲基环氧乙烷的 ECD 计算

配套教程：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)（第 8 节）。

## 目的

用最小的手性分子演示 ECD 计算的完整链条：几何优化与频率 → 带弥散函数的 TD-DFT → 读取 R(velocity) → 展宽成 Δε 谱 → 得到对映体的镜像谱。甲基环氧乙烷（propylene oxide，C3H6O，10 个原子）没有构象自由度，不需要 Boltzmann 加权，适合先把流程跑通。

## 文件说明

| 文件 | 内容 |
|---|---|
| `chiral_molecule_ecd_td.gjf` | 第 1 步 B3LYP-D3(BJ)/def2-TZVP Opt+Freq（气相）；第 2 步 CAM-B3LYP/aug-cc-pVDZ `TD=(NStates=20)` |
| `run.slurm` | Slurm 作业脚本，运行后对 chk 执行 formchk |

坐标由 SMILES `C[C@@H]1CO1` 生成，经 RDKit 判定手性中心（2 号 C）为 R 构型。

## 运行步骤

```bash
sbatch run.slurm
# 或
g16 < chiral_molecule_ecd_td.gjf > chiral_molecule_ecd_td.log
```

计算量很小，8 核上通常几十分钟内完成。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" chiral_molecule_ecd_td.log     # 应为 3
grep "NImag" chiral_molecule_ecd_td.log                      # NImag=0
grep "Excited State" chiral_molecule_ecd_td.log | awk '{print $3, $5, $7, $9}'
grep -A 25 "R(velocity)" chiral_molecule_ecd_td.log | head -28
grep -A 25 "R(length)"   chiral_molecule_ecd_td.log | head -28
```

- 预计 20 个激发态全部位于真空紫外区（波长短于 200 nm），低能态有明显的 Rydberg 特征，这正是必须使用 aug-cc-pVDZ 的原因；
- 比较 R(velocity) 与 R(length)：两者符号和量级大体一致说明基组基本够用；差别很大时应换 aug-cc-pVTZ 复核；
- 展宽：`python ../../scripts/spectrum_uvvis_ecd.py -h` 查看用法，或用 Multiwfn 主功能 11 选择 ECD，或 SpecDis。

## 生成 (S)-对映体（可选）

把所有 z 坐标取反就得到镜像分子。下面的命令从输入文件中提取坐标并翻转 z（只适用于本文件这种“元素 x y z”格式）：

```bash
awk 'NF==4 && $1 ~ /^[A-Z][a-z]?$/ {printf " %-2s  %12.6f %12.6f %12.6f\n", $1, $2, $3, -$4}' \
    chiral_molecule_ecd_td.gjf > mox_S_coords.txt
```

用 `mox_S_coords.txt` 替换输入中的坐标块重算，得到的旋转强度应与 R 构型大小相同、符号相反（数值噪声范围内）。实际工作中直接把 R 构型的 R 值取负即可，无需重算。

## 注意事项

- 本例为气相计算，只用于演示流程；与溶液实验比较时需加入溶剂模型（见第 9 章、第 10 章第 10 节）。
- 柔性分子必须先做构象搜索并按自由能加权，单一构象的 ECD 结论不可靠。
- 甲基环氧乙烷的实验 ECD 位于真空紫外区，常规 ECD 仪器测不到；做实际课题时请选用在近紫外有吸收的发色团。

返回：[第10章 激发态 TD-DFT](../../tutorials/10-tddft.md)
