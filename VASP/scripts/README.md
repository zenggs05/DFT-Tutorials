# VASP 辅助脚本（VASP helper scripts）

本目录收录 VASP 计算的前后处理脚本。所有 Python 脚本均支持 `-h` 查看完整参数说明，Bash 脚本可用 `bash <脚本> -h` 或阅读文件头部注释。

| 脚本 | 功能 |
|---|---|
| [check_vasp.py](check_vasp.py) | check_vasp.py -- health check of a (running or finished) VASP job. |
| [chgdiff.py](chgdiff.py) | chgdiff.py -- charge-density difference and Bader helper for VASP CHGCAR-type files. |
| [conv_analyze.py](conv_analyze.py) | conv_analyze.py -- analyse VASP convergence tests (ENCUT / k-points / KSPACING / |
| [conv_test.sh](conv_test.sh) | conv_test.sh -- ENCUT / k-point / KSPACING / SIGMA (or any INCAR tag) |
| [elastic_post.py](elastic_post.py) | elastic_post.py -- post-process VASP elastic constants (IBRION = 6, ISIF >= 3). |
| [locpot_workfunction.py](locpot_workfunction.py) | locpot_workfunction.py -- planar/macroscopic average of a VASP LOCPOT and work function. |
| [md_analysis.py](md_analysis.py) | md_analysis.py -- analyse VASP molecular-dynamics runs (AIMD / MLFF). |
| [neb_tools.py](neb_tools.py) | neb_tools.py -- helper for VASP (VTST) nudged-elastic-band calculations. |
| [plot_band.py](plot_band.py) | plot_band.py -- plot VASP band structures and report the band gap. |
| [plot_dos.py](plot_dos.py) | plot_dos.py -- plot total / projected DOS from a VASP DOSCAR and compute d-band centres. |
| [vasp_incar_writer.py](vasp_incar_writer.py) | vasp_incar_writer.py -- write INCAR.static for static (DOS) runs of |
| [vasp_inputgen.py](vasp_inputgen.py) | vasp_inputgen.py -- generate VASP input files from a POSCAR for common tasks. |

## 依赖

```bash
pip install numpy matplotlib scipy ase   # 基础依赖
pip install pymatgen phonopy seekpath spglib cclib   # 部分脚本的可选依赖
```

## 用法示例

```bash
python check_vasp.py -h
python chgdiff.py -h
python conv_analyze.py -h
bash conv_test.sh -h
python elastic_post.py -h
python locpot_workfunction.py -h
```

> [!NOTE]
> 脚本按真实输出格式编写，但不同软件版本的输出可能略有差异；使用前建议先在小体系上核对结果。

[返回目录](../README.md)
