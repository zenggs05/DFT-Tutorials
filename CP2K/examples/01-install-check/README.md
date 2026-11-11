# 示例 01：CP2K 安装自检（Installation Check & H2O Smoke Test）

## 目的

用一个几秒钟就能跑完的孤立 H2O 分子 PBE 单点能计算，确认 CP2K 安装真正可用：可执行文件能启动、`cp2kflags` 中包含常用功能、数据目录（`CP2K_DATA_DIR`）里有基组和赝势、SCF 能收敛、MPI 并行没有错配。对应教程 [第1章 CP2K 简介、版本与安装](../../tutorials/01-install.md)。

## 文件说明

| 文件 | 说明 |
|---|---|
| `env_cp2k.sh` | 环境设置脚本（必须 `source`）。按安装方式 `CP2K_METHOD=static/conda/toolchain/spack/module` 设置 `PATH`、`CP2K_DATA_DIR`、`CP2K_EXE`。使用前请编辑文件开头的路径 |
| `check_install.sh` | 自检脚本：检查可执行文件与版本、`cp2kflags`、数据文件、语法检查（`-c`）、运行冒烟测试、检测 MPI 错配；在临时目录中运行，不污染当前目录 |
| `h2o_smoke_test.inp` | 冒烟测试输入：孤立 H2O，PBE，DZVP-MOLOPT-GTH + GTH-PBE，Martyna–Tuckerman 泊松求解器，适用于 CP2K 2023.x–2026.2 |

## 运行步骤

```bash
cd CP2K/examples/01-install-check

# 1) 编辑 env_cp2k.sh 开头的 CP2K_METHOD 和路径，然后加载环境
source env_cp2k.sh
#    或者临时指定：CP2K_METHOD=conda source env_cp2k.sh

# 2) 串行（或 ssmp）自检
bash check_install.sh

# 3) MPI 自检：4 个进程（默认用 mpirun，可改为 srun / mpiexec）
NP=4 bash check_install.sh
NP=4 MPI_LAUNCHER=srun bash check_install.sh      # 在 Slurm 交互作业中

# 4) 指定某个可执行文件
bash check_install.sh $HOME/opt/cp2k-2026.2/bin/cp2k.ssmp

# 5) 也可以直接手动运行冒烟测试
rm -f h2o_smoke_test.out
cp2k.ssmp -i h2o_smoke_test.inp -o h2o_smoke_test.out
```

在集群上请通过交互作业或批处理作业在**计算节点**上运行（登录节点通常禁止跑 MPI 程序）。

## 预期结果 / 如何检查

`check_install.sh` 最后一行为 `RESULT: PASS` 即通过。手动检查时：

```bash
grep "CP2K| Data directory path" h2o_smoke_test.out   # 实际使用的数据目录
grep "SCF run converged" h2o_smoke_test.out           # 应在 20 步以内收敛
grep "ENERGY| Total FORCE_EVAL" h2o_smoke_test.out    # 总能量（Hartree）
grep -c "DBCSR| MPI: Number of processes" h2o_smoke_test.out   # 必须为 1
grep "PROGRAM ENDED AT" h2o_smoke_test.out            # 正常结束
```

- 总能量约为 −17.2 Ha（示意，具体数值随版本、基组和网格略有不同）。同一台机器上不同安装方式（静态二进制、conda、容器、源码编译）的结果应在 $10^{-6}$ Ha 量级一致。
- 2024.x 及更早版本的能量行是 `ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:`，2025.1 起是 `... energy [hartree]`；`grep "ENERGY| Total FORCE_EVAL"` 对两种格式都适用。

## 注意事项

- `cp2k -o file` 会**追加**写入已存在的文件，重跑前先删除旧输出（脚本已自动处理）。
- 若 `DBCSR| MPI: Number of processes` 出现多次且每次都是 1，说明启动用的 MPI 与 CP2K 链接的 MPI 不一致，见教程第1章“常见问题与报错”。
- 输入中故意同时列出 `BASIS_MOLOPT` 和 `BASIS_MOLOPT_UZH`，以顺便检查数据目录是否完整；若你的旧版本数据目录缺少 UZH 文件，删除那一行即可。
- 冒烟测试的 `CUTOFF 300` 只为“跑得快”，不是生产计算的推荐值；截断收敛见 [第4章](../../tutorials/04-cutoff-convergence.md)。

返回教程：[第1章 CP2K 简介、版本与安装](../../tutorials/01-install.md)
