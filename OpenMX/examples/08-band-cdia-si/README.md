# 算例 8-1：金刚石与硅的能带结构（Band.kpath + bandgnu13 + Python）

## 目的

用 OpenMX 4.0 计算金刚石 C 与 Si（fcc 原胞）的能带：在同一次运行中完成 SCF 和沿高对称路径的能带计算，用 `bandgnu13` 生成 gnuplot 数据，再用 `band_parse.py` 直接解析 `System.Name.Band`，得到 VBM、CBM、带隙类型并作图。两个输入都用 `<Band.KPath.UnitCell` 以惯用立方胞定义 k 路径，避免 fcc 高对称点标签错位。

对应教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)。

## 文件说明

| 文件 | 作用 |
|---|---|
| `Cdia_band.dat` | 金刚石（a = 3.567 Å），`C6.0-s2p2d1`，10×10×10 k 点；完整 Setyawan–Curtarolo fcc 路径 Γ–X–W–K–Γ–L–U–W–L–K\|U–X（10 段）；`System.Name cdia` |
| `Si_band.dat` | 硅（a = 5.43 Å），`Si7.0-s2p2d1`；路径 L–Γ–X–W–K–Γ，与 [能带展开算例](../08-unfold-si-supercell/README.md) 相同；`System.Name si` |
| `band_parse.py` | 读取 `.Band`（Hartree，含倒格矢和路径信息），换算为以化学势为零点的 eV，输出 VBM/CBM/带隙并作图 |
| `job.slurm` | 依次运行两个输入，调用 `bandgnu13` 与 `band_parse.py` |

## 运行步骤

1. 修改两个 `.dat` 中的 `DATA.PATH`。
2. 编译 `bandgnu13`（只需一次，它不在 `make install` 中）：

   ```bash
   cd /path/to/openmx4.0/source
   gcc bandgnu13.c -lm -o bandgnu13
   cp bandgnu13 ~/bin/          # or anywhere in $PATH
   ```

3. 运行与后处理：

   ```bash
   cd OpenMX/examples/08-band-cdia-si
   mpirun -np 8 openmx Cdia_band.dat -nt 1 | tee Cdia_band.std
   mpirun -np 8 openmx Si_band.dat   -nt 1 | tee Si_band.std
   bandgnu13 cdia.Band               # -> cdia.GNUBAND, cdia.BANDDAT1
   gnuplot -p cdia.GNUBAND           # quick look
   python3 band_parse.py cdia.Band --emin -22 --emax 15 --plot cdia_band.png
   python3 band_parse.py si.Band   --emin -13 --emax 6  --plot si_band.png
   ```

## 预期结果 / 如何检查

- `grep -i "chemical potential" cdia.out` 能看到化学势；`ls cdia.Band` 存在说明能带部分正常完成。
- `band_parse.py` 输出格式（示意输出，数值仅供参考）：

  ```text
  10 segments, 193 k points, 1 spin channel(s), xx bands
  ticks: G=0.0000  X=0.xxxx  W=0.xxxx  K=0.xxxx  G=x.xxxx  L=x.xxxx  U=x.xxxx  W=x.xxxx  L=x.xxxx  K|U=x.xxxx  X=x.xxxx
  E_F (chemical potential) = x.xxxx eV (set to 0)
  VBM =  -x.xxxx eV at G k=(0.000 0.000 0.000) (spin 0, band 4)
  CBM =   x.xxxx eV at k=(0.xxx 0.xxx 0.xxx) (spin 0, band 5)
  band gap along the path = x.xxxx eV (indirect)
  ```

- 参考值（近似）：两者都是间接带隙，VBM 在 Γ，CBM 在 Γ–X 线上靠近 X 处（Si 约在 0.85 X 附近）。PBE 带隙：金刚石约 4.1–4.2 eV（实验 5.48 eV），Si 约 0.5–0.6 eV（实验约 1.17 eV）。GGA 低估带隙是已知现象，不是计算错误；具体数值还与基组、晶格常数有关。
- 注意 `.Band` 中的 k 坐标已换算为 SCF 晶胞（原胞）的倒格矢单位，所以 CBM 的分数坐标看起来与惯用胞写法不同，这是正常的。

## 注意事项

- 能带只在 `scf.EigenvalueSolver Band` 时计算；`Band.Nkpath` 必须等于 `<Band.kpath` 的行数。
- 去掉 `<Band.KPath.UnitCell` 块而不改坐标，同样的数字会被解释为原胞倒格矢坐标，X、W、K 等标签就错了——可以亲手试一次体会差别。
- 若已有收敛的 SCF，只想换路径重算能带：保持 `System.Name` 不变，加 `scf.restart pp`（OpenMX ≥ 4.0），只做一次对角化，不覆盖 `_rst` 文件。

返回教程：[第8章 能带结构（含自动 k 路径与能带展开）](../../tutorials/08-band.md)
