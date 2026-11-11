# 算例 12-mos2-soc-band：单层 MoS2 的 SOC 价带劈裂与自旋–能谷锁定

## 目的

- 计算单层 2H-MoS2 含 SOC 的能带，读出 K 点价带顶的自旋劈裂（无反演对称性 + SOC 引起）；
- 比较 K 与 K'（本例标签 `Kp`）两个能谷：劈裂大小相同、自旋方向相反（时间反演联系两者）；
- 用 kSpin 的 BandDispersion 模式沿能带路径计算动量空间自旋密度矩阵（自旋织构）。

对应教程：[第12章 自旋轨道耦合与拓扑性质](../../tutorials/12-soc-topology.md) §4、§7。

## 文件说明

| 文件 | 说明 |
|---|---|
| `MoS2_SOC_band.dat` | 单层 MoS2（a = 3.16 Å，Mo–S 垂直距离 1.565 Å，均为近似值；c = 20 Å 含真空），Mo7.0-s3p2d2 / S7.0-s2p2d1f1（Standard），nc + SOC，k 网格 12×12×1，`HS.fileout on`；路径 Γ–M–K–Γ–K'；文件末尾附有 kSpin 关键词（openmx 会忽略它们） |
| `job.slurm` | SCF + 能带 → 用 `../12-gaas-soc-band/soc_splitting.py` 读 Γ、K、K' 的能级 → 若找到 kSpin 则运行 BandDispersion |

运行前修改 `DATA.PATH`。kSpin 需单独编译：`cd ~/openmx4.0/source && make kSpin`（若没有自动复制到 work，请手动复制）。

## 运行步骤

```bash
sbatch job.slurm
# or by hand:
mpirun -np 8 openmx MoS2_SOC_band.dat -nt 1 > mos2_soc.std 2>&1
python ../12-gaas-soc-band/soc_splitting.py mos2_soc.Band --at G K Kp
mpirun -np 4 ~/openmx4.0/work/kSpin MoS2_SOC_band.dat | tee kspin.log
```

## 预期结果 / 如何检查

- 带隙为 K 点直接带隙（GGA 约 1.6–1.8 eV 量级，依赖几何；单层的光学带隙实验约 1.8–1.9 eV，二者不应直接比较）；
- `soc_splitting.py` 在 K 与 Kp 处打印的 "splitting of the two highest valence groups" 应约为 0.15 eV（文献 GGA 值约 148 meV，见 Zhu et al., Phys. Rev. B 84, 153402 (2011)），K 与 Kp 的数值相同；导带底的劈裂只有几 meV；
- Γ 是时间反演不变点，Kramers 定理保证 Γ 处每个能级仍是二重简并（`x2`）；而 K、K' 不是时间反演不变点，SOC 可以把价带顶劈裂成两个非简并能级（`x1`）；
- kSpin 运行后生成 `mos2_soc_BD.Band*_*`、`mos2_soc_BD.AMulPBand*` 等文件，包含沿路径各能带的自旋期望值；K 与 K' 附近最高价带的 σz 符号应相反（自旋–能谷锁定）。各列含义请先 `head` 查看文件头，并对照手册第 56.4、56.6 节。

具体数值以你的计算为准。

## 注意事项

- 2D 体系沿 c 方向只取 1 个 k 点；真空层要足够厚，否则层间相互作用会改变能带。
- 本例几何未优化；定量比较前请按第6章优化晶格常数与 Mo–S 距离（注意优化时通常不开 SOC）。
- SOC 计算的哈密顿量维度是共线计算的 2 倍，对角化成本约为 8 倍，k 点与基组要权衡。
- `Kp` 是本例为 K' 取的标签名（避免 shell 与命令行把 `-K`、`K'` 当作特殊字符）。

[返回第12章](../../tutorials/12-soc-topology.md)
