# 示例 15-2：QM/MM 电子嵌入——苯酚与 6 个 TIP3P 水（ONIOM(B3LYP:Amber)=EmbedCharge）

配套教程：[第15章 ONIOM 多层方法与 QM/MM](../../tutorials/15-oniom.md)（§3.3、§6、§11）

## 目的

- 写出 Amber 低层所需的原子类型与部分电荷（`O-OW--0.834000`、`H-HW-0.417000`、`C-CA` 等）；
- 使用冻结标记（−1）固定外层水分子；
- 使用 `Geom=Connectivity` 连接表；
- 用 `Pop=(Uncharged,SaveCM5)` 为 QM 区原子拟合 CM5 电荷并存入 chk；
- 按“先机械嵌入、再电子嵌入”的顺序优化。

## 文件说明

| 文件 | 说明 |
|---|---|
| `oniom2_qmmm_amber_embedcharge.gjf` | 3 步作业（同一个 chk）：① ME 单点 + `Pop=(Uncharged,SaveCM5)`；② ME 优化（`Geom=AllCheck Guess=Read`）；③ EE 优化（`=EmbedCharge`） |
| `run_g16_slurm.sh` | Slurm 作业脚本（8 核、10 GB） |

体系：苯酚（原子 1–13，高层 B3LYP/6-31G(d)）+ 6 个水（原子 14–31，低层 Amber）。水 1（原子 14–16）接受苯酚 O–H 的氢键，水 2（17–19）向苯酚 O 给出氢键，水 3（20–22）在苯环上方，水 4–6（23–31）冻结。QM 与 MM 之间没有共价键，因此没有链接原子。起始结构用 RDKit/MMFF 预优化。

| 原子 | Amber 类型 | 电荷 |
|---|---|---|
| 苯酚 O / 羟基 H | `OH` / `HO` | 不写，由第 1 步拟合 CM5 |
| 与 O 相连的芳碳 / 其余芳碳 / 芳环 H | `C` / `CA` / `HA` | 不写，由第 1 步拟合 CM5 |
| 水 O / 水 H | `OW` / `HW` | −0.834 / +0.417（TIP3P） |

苯酚的类型沿用 Amber 中酪氨酸侧链的原子类型。

## 运行步骤

```bash
cd Gaussian/examples/15-qmmm-phenol-water
g16 < oniom2_qmmm_amber_embedcharge.gjf > oniom2_qmmm_amber_embedcharge.log
# or
sbatch run_g16_slurm.sh oniom2_qmmm_amber_embedcharge.gjf
```

8 核上一般在一小时以内（取决于两次优化的步数）。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" oniom2_qmmm_amber_embedcharge.log          # expect 3
grep -i "undefined\|MM function not complete" oniom2_qmmm_amber_embedcharge.log   # expect no output
grep "ONIOM: extrapolated energy" oniom2_qmmm_amber_embedcharge.log | tail -1
grep -c "Optimization completed" oniom2_qmmm_amber_embedcharge.log      # expect 2
```

| 检查项 | 预期 |
|---|---|
| 第 1 步 | 正常结束，输出中有 Hirshfeld/CM5 电荷，QM 原子的 CM5 电荷被保存 |
| 第 2、3 步 | 两次优化都收敛；冻结的水 4–6 坐标不变 |
| 氢键 | 比较第 2、3 步优化结构中苯酚 O–H···O(水 1) 与 O(水 2)–H···O(苯酚) 的距离；电子嵌入下 QM 区被水的点电荷极化，氢键长度通常会有变化 |
| 能量 | 第 2、3 步的外推能量定义不同（嵌入方式不同），不要直接相减 |

## 注意事项

- 若报 `Angle bend undefined between atoms ...` / `MM function not complete`，说明当前 Gaussian 版本内置的 Amber 参数缺少某个组合。把第 1 步的低层写成 `Amber=HardFirst`，在连接表之后空一行补充缺失的 `HrmStr1`/`HrmBnd1` 参数（参数值需来自文献，例如 Amber 力场参数文件），再空一行结束；详见教程 §3.3 和官方 [Molecular Mechanics Methods](https://gaussian.com/mm/) 页面。
- 第 2、3 步使用 `Geom=AllCheck`，结构、层、原子类型、电荷、冻结信息和连接表都从 `phenol_qmmm.chk` 读取，所以三步必须用同一个 chk。
- 机械嵌入下 QM 区原子的 MM 电荷决定 QM–MM 静电作用，这正是第 1 步拟合 CM5 电荷的原因；电子嵌入下这些电荷对总能量的贡献基本相互抵消，但仍用于微迭代。
- 真实的溶液体系需要更大的水壳层并对多个构型取样；本例只演示输入格式与流程。

返回教程：[第15章](../../tutorials/15-oniom.md)
