# 示例 12-esp-resp-paracetamol：小分子配体的 MK 静电势与 RESP 电荷

配套教程：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)（第 7 节）。

## 目的

按 AMBER/GAFF 的经典协议，为对乙酰氨基酚（paracetamol，C8H9NO2）生成 RESP 电荷：

1. HF/6-31G(d) 几何优化；
2. HF/6-31G(d) 单点，`Pop=MK` 加 `IOp(6/33=2,6/41=10,6/42=17,6/50=1)` 输出 ESP 格点并写出 Antechamber 的 `.gesp` 文件；
3. 用 AmberTools `antechamber -c resp` 做两步 RESP 拟合（或改用 Multiwfn 的 RESP 模块）。

## 文件说明

| 文件 | 内容 |
|---|---|
| `ligand_mk_esp_for_resp.gjf` | 两步 Link1：HF/6-31G(d) Opt → HF/6-31G(d) MK ESP（写出 `paracetamol.gesp`） |
| `run_antechamber_resp.sh` | 调用 antechamber 拟合 RESP，并列出各原子电荷及总和 |
| `run.slurm` | Gaussian 部分的 Slurm 作业脚本 |

## 运行步骤

```bash
# 1) Gaussian
sbatch run.slurm
# 或
g16 < ligand_mk_esp_for_resp.gjf > ligand_mk_esp_for_resp.log

# 2) RESP 拟合（需要 AmberTools）
bash run_antechamber_resp.sh paracetamol.gesp 0 LIG
# 如果没有生成 .gesp，也可以直接读 log（需要 IOp(6/33=2) 打印的格点）
bash run_antechamber_resp.sh ligand_mk_esp_for_resp.log 0 LIG
```

`paracetamol.gesp` 通常写在作业的工作目录中；如果找不到，请在 scratch 目录（`$GAUSS_SCRDIR`）中查找。

Multiwfn 替代方案：

```bash
formchk lig_esp.chk lig_esp.fchk
Multiwfn lig_esp.fchk      # 依次输入 7 -> 18 -> 1，进行标准两步 RESP 拟合
```

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" ligand_mk_esp_for_resp.log      # 应为 2
grep "Charges from ESP fit" ligand_mk_esp_for_resp.log       # 拟合 RMS
grep -A22 "ESP charges:" ligand_mk_esp_for_resp.log | head -24   # Gaussian 自带的 MK 电荷（未加约束）
```

- `run_antechamber_resp.sh` 最后打印的电荷总和应等于净电荷 0（四舍五入误差在 10⁻⁴ 量级）；
- 与 MK 电荷相比，RESP 会把埋藏的碳原子（如苯环碳、甲基碳）电荷拉向更小的绝对值，而甲基上的三个 H 被约束为相同电荷；
- 羰基 O、酚羟基 O 和酰胺 N 应带明显的负电荷，羟基 H 和酰胺 H 带正电荷。

## 注意事项

- HF/6-31G(d) 是 AMBER/GAFF 的配套协议，不要为“提高精度”而改用 DFT 大基组；为其他力场拟合电荷时请遵循该力场自己的协议。
- 电荷依赖构象：酰胺基团与羟基的取向会影响结果，柔性配体建议多构象联合拟合（Multiwfn RESP 模块的选项 −1）。
- `antechamber -nc` 必须与 Gaussian 输入中的总电荷一致。
- `IOp(6/33=2)` 是 AMBER 社区长期使用的写法，当前官方 IOp 列表中没有列出；`IOp(6/50=1)` 生成 .gesp 文件则是官方 Population 页说明的方法。

返回：[第12章 布居分析、NBO 与静电势/RESP 电荷](../../tutorials/12-population-charges.md)
