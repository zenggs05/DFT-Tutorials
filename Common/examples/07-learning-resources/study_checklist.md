# 12 周学习打卡表（DFT-Tutorials）

对应 [第7章](../../tutorials/07-learning-path.md) 第 3.9 节。每完成一项在方括号里打 `x`，并把产出（图、表、目录路径）记在“记录”一栏。Gaussian 方向可把第 3–8 周换成 Gaussian 第1–12章。

| 周 | 任务 | 完成 | 记录 |
|---|---|---|---|
| 1 | Common 第1–2章：WSL/集群登录、conda 环境、作业脚本 | [ ] | |
| 1 | 阅读 Sholl & Steckel 第1–3章 | [ ] | |
| 2 | Common 第3–4章：下载并转换结构；ASE/pymatgen 读写 | [ ] | |
| 2 | VASP 第1–2章：写出 Si 的四大输入并跑通单点 | [ ] | |
| 3 | VASP 第3章 + Common 第7章第4节：ENCUT、k 点、SIGMA 收敛测试 | [ ] | |
| 4 | VASP 第4章：EOS 与结构优化，晶格常数与文献对比 | [ ] | |
| 5 | VASP 第5章 + Common 第5章：能带 + PDOS 出版级作图 | [ ] | |
| 6 | VASP 第6–7章（按课题选）：HSE 带隙或磁性基态 | [ ] | |
| 7 | VASP 第8–9章：二维材料或表面 slab | [ ] | |
| 8 | VASP 第10章：吸附能 | [ ] | |
| 9 | VASP 第11章 + Common 第7章第5.3节：NEB（MLIP 初猜） | [ ] | |
| 10 | VASP 第12章（声子）或第15章（AIMD/MLFF） | [ ] | |
| 11 | Common 第6章：批量计算；Common 第7章第5节：MACE-MP 预优化 | [ ] | |
| 12 | 复现一篇文献的关键图，整理输入、脚本、版本记录，组会汇报 | [ ] | |

## 每次计算的自检清单

- [ ] 截断能、k 点（或基组、格点）已做收敛测试，判据写在笔记中
- [ ] 要相减的计算使用完全相同的参数
- [ ] 输出正常结束（VASP `General timing`、Gaussian `Normal termination`、CP2K 无 `SCF run NOT converged`）
- [ ] 结构优化后检查了最大受力 / 压力
- [ ] 记录了软件版本、赝势/基组版本、作业脚本
- [ ] MLIP 结果只作初猜，最终数值来自 DFT
