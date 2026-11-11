# 算例 9-1：金刚石的总态密度与投影态密度（Dos.fileout + DosMain）

## 目的

用 OpenMX 4.0 计算金刚石的总 DOS 与 PDOS：SCF 时打开 `Dos.fileout` 输出 `cdia_dos.Dos.val`（本征值）和 `cdia_dos.Dos.vec`（本征矢），再用 `DosMain` 以四面体法生成 `cdia_dos.DOS.Tetrahedron` 与 `cdia_dos.PDOS.Tetrahedron.atomN.*`。练习交互式与非交互式（答案文件 / here-doc / 驱动脚本）运行 `DosMain`，以及用 Python 按元素、角动量求和并检查积分电子数。

对应教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Cdia_dos.dat` | 金刚石（a = 3.567 Å）原胞，`C6.0-s2p2d1`，`scf.Kgrid 10 10 10`，`Dos.Erange -30 20`，`Dos.Kgrid 16 16 16`；`System.Name cdia_dos` |
| `dosmain_answers_dos.txt` | `DosMain` 的回答：四面体法（1）、总 DOS（1） |
| `dosmain_answers.txt` | `DosMain` 的回答：四面体法（1）、PDOS（2）、原子 `1 2` |
| `run_dosmain.sh` | 通用驱动：读取 `.Dos.val` 头部的 `mode`、`Kgrid`、`atomnum`，自动生成回答，依次计算总 DOS 与 PDOS（也适用于第 9-2、9-3 等算例） |
| `sum_pdos.py` | 把 PDOS 文件按元素 / 物种 / 原子与 s、p、d、f 求和（可按 m 分辨），打印费米能以下的积分电子数，检查 ΣPDOS = DOS，可作图 |
| `job.slurm` | Slurm 脚本：OpenMX → DosMain（两次）→ `sum_pdos.py` |

## 运行步骤

1. 修改 `Cdia_dos.dat` 中的 `DATA.PATH`。
2. 编译 `DosMain`（一次即可）：

   ```bash
   cd /path/to/openmx4.0/source && make DosMain && cp DosMain ~/bin/
   # official Debian package instead:
   # export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
   ```

3. 计算与后处理：

   ```bash
   cd OpenMX/examples/09-dos-pdos-cdia
   mpirun -np 8 openmx Cdia_dos.dat -nt 1 | tee cdia_dos.std
   DosMain cdia_dos.Dos.val cdia_dos.Dos.vec < dosmain_answers_dos.txt   # total DOS
   DosMain cdia_dos.Dos.val cdia_dos.Dos.vec < dosmain_answers.txt       # PDOS of atoms 1 2
   # or simply:  bash run_dosmain.sh cdia_dos tetra
   python3 sum_pdos.py cdia_dos --dat Cdia_dos.dat --plot --emin -28 --emax 15
   ```

4. （可选）交互式体验一次：`DosMain cdia_dos.Dos.val cdia_dos.Dos.vec`，依次回答 `1`、`2`、`1 2`。

## 预期结果 / 如何检查

- 生成的文件：`cdia_dos.DOS.Tetrahedron`，`cdia_dos.PDOS.Tetrahedron.atom1`、`.atom1.s1`、`.atom1.p1..p3`、`.atom1.d1..d5`（atom2 同理）。列为：能量（eV，化学势为零点）、DOS（states/eV）、积分 DOS。
- `awk '$1<=0 {n=$3} END{print n}' cdia_dos.DOS.Tetrahedron` 应接近 8（两个 C 原子的价电子数），因为 `Dos.Erange` 下限 −30 eV 低于价带底。
- `sum_pdos.py` 输出（示意输出，数值仅供参考）：

  ```text
  method: Tetrahedron, 18 files, spin-polarized: False
  C d                          -> cdia_dos.pdos_C_d.dat              N(E<E_F) = 0.0xx
  C p                          -> cdia_dos.pdos_C_p.dat              N(E<E_F) = x.xxx
  C s                          -> cdia_dos.pdos_C_s.dat              N(E<E_F) = x.xxx
  max |sum(PDOS) - DOS| = 0.00xx states/eV (small only if PDOS was done for ALL atoms)
  figure written to cdia_dos_pdos.png
  ```

- s + p + d 的积分之和约为 8；s、p 并不严格等于 2×1 和 2×3（杂化导致电荷在 s、p、d 之间重新分配，手册 23.1 节也有说明）。费米能附近约 4 eV 的范围内 DOS 为零（PBE 带隙）。

## 注意事项

- `Dos.Erange` 以化学势为零点；积分 DOS 从下限开始累加，下限高于某些占据态（尤其半芯态）时，0 eV 处的积分值会小于价电子数。占据数请以 `.out` 中的 Mulliken 布居为准。
- 自旋极化计算的 PDOS 文件有 5 列，第 3 列（下自旋 DOS）在 4.0 的 `DosMain` 中带负号；`sum_pdos.py` 保留该约定，积分时取绝对值。
- `DosMain` 的提问顺序随情况变化：k 网格为 1×1×1（分子）或 DC 计算时不询问方法、只询问高斯宽度；单原子体系不询问原子编号。写答案文件前先确认，或直接用 `run_dosmain.sh`。
- `.Dos.vec` 可能非常大（k 点 × 本征态 × 基函数），大体系请考虑 `DosGauss.fileout`（教程 §9）。

返回教程：[第9章 态密度 DOS/PDOS、COHP 与光学性质](../../tutorials/09-dos-pdos-cohp-optics.md)
