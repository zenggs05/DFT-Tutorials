# 示例 07：用 QST2 / QST3 搜索 HCN ⇌ HNC 异构化过渡态

## 目的

演示 STQN 方法的两种输入方式：

- **QST2**：只给反应物 HCN 和产物 HNC，由程序插值生成 TS 初猜；
- **QST3**：额外给出 H 桥连 C 和 N 的 TS 猜测结构。

两者都在优化后自动做 `Freq`，用来确认恰好 1 个虚频。对应教程：[第7章 过渡态搜索与 IRC 验证](../../tutorials/07-ts-irc.md) 第 6 节。

## 文件说明

| 文件 | 说明 |
|---|---|
| `hcn_hnc_qst2.gjf` | `Opt=(QST2,CalcFC) Freq`。两段“标题 + 电荷多重度 + 坐标”：反应物、产物 |
| `hcn_hnc_qst3.gjf` | `Opt=(QST3,CalcFC) Freq`。三段：反应物、产物、TS 猜测（顺序固定） |

所有结构中的原子顺序都是 **C、N、H**，第 3 个原子就是迁移的 H。改动结构时必须保持这一对应关系。

几何说明：

- 反应物和产物中的 H 都偏离分子轴约 10°，并且偏向同一侧（+x）。两端都严格线性时，插值平面无法确定，冗余内坐标也容易出现线性角问题；
- TS 猜测按“C–H ≈ 1.19 Å、∠HCN ≈ 72°”手工搭建，只是近似初猜，不是参考结构。

理论水平：B3LYP/6-311G(d,p)。体系很小，单核几分钟内即可完成。

## 运行步骤

```bash
cd Gaussian/examples/07-qst-hcn-hnc
g16 < hcn_hnc_qst2.gjf > hcn_hnc_qst2.log
g16 < hcn_hnc_qst3.gjf > hcn_hnc_qst3.log
```

在集群上运行时，可以参考 [../07-ts-irc-sn2/run_ts_irc.sh](../07-ts-irc-sn2/run_ts_irc.sh) 的写法，或使用仓库脚本 [`sub_g16.sh`](../../scripts/sub_g16.sh)（`bash ../../scripts/sub_g16.sh -h` 查看用法）。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" hcn_hnc_qst2.log          # 期望 2
grep "imaginary frequencies" hcn_hnc_qst2.log          # 期望 1 个虚频
grep -A 12 "Optimized Parameters" hcn_hnc_qst2.log | head -n 16
```

- `Optimized Parameters` 表比普通优化多出 `Reactant` 和 `Product` 两列，可以对照 TS 中 C–H、N–H 距离相对两端的位置；
- GaussView 中播放虚频振动，应看到 H 在 C 与 N 之间来回摆动；
- QST2 和 QST3 应收敛到同一个 TS（结构和能量在收敛精度内一致）。若不一致，用频率模式和 IRC 判断哪个才是 HCN ⇌ HNC 的 TS；
- 进一步验证：仿照 `07-ts-irc-sn2/sn2_irc_both_directions.gjf` 写一个 IRC 输入，`%oldchk` 指向 `hcn_hnc_qst2.chk`，电荷多重度为 `0 1`。

本仓库没有附带输出，所有数值请以你自己的计算为准。

## 注意事项

- route 中**不要**同时写 `TS` 和 `QST2/QST3`（官方文档明确禁止）。
- 与 `ModRedundant` 联用时，每个结构后面都必须有一段 ModRedundant 输入（QST3 就是三段）。
- 如果 QST2 报 `FormBX had a problem` 或线性角相关错误，改用 QST3，或直接在 TS 猜测上运行 `Opt=(TS,CalcFC,NoEigenTest)`。
- 文件末尾必须有空行，否则 Gaussian 会报读入错误。

[← 返回第7章](../../tutorials/07-ts-irc.md)
