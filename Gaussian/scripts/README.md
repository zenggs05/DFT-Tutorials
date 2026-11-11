# Gaussian 辅助脚本（Gaussian helper scripts）

本目录收录 Gaussian 计算的前后处理脚本。所有 Python 脚本均支持 `-h` 查看完整参数说明，Bash 脚本可用 `bash <脚本> -h` 或阅读文件头部注释。

| 脚本 | 功能 |
|---|---|
| [g16_log_summary.py](g16_log_summary.py) | g16_log_summary.py - summarize Gaussian 16 output files (.log/.out) into a table, CSV or JSON. |
| [gjf_builder.py](gjf_builder.py) | gjf_builder.py - generate Gaussian 16 input files (.gjf) from xyz, Gaussian output or existing gjf files. |
| [log2xyz.py](log2xyz.py) | log2xyz.py - extract geometries from Gaussian 16 opt / scan / IRC output into (multi-frame) xyz, |
| [opt_monitor.py](opt_monitor.py) | opt_monitor.py - follow a Gaussian 16 geometry optimization step by step. |
| [spectrum_ir_raman.py](spectrum_ir_raman.py) | spectrum_ir_raman.py -- IR and Raman spectra from Gaussian harmonic frequency |
| [spectrum_uvvis_ecd.py](spectrum_uvvis_ecd.py) | spectrum_uvvis_ecd.py -- UV-Vis absorption and ECD spectra from Gaussian |
| [thermo_qrrho.py](thermo_qrrho.py) | thermo_qrrho.py -- teaching-sized quasi-RRHO thermochemistry from Gaussian 16 |

## 依赖

```bash
pip install numpy matplotlib scipy ase   # 基础依赖
pip install pymatgen phonopy seekpath spglib cclib   # 部分脚本的可选依赖
```

## 用法示例

```bash
python g16_log_summary.py -h
python gjf_builder.py -h
python log2xyz.py -h
python opt_monitor.py -h
python spectrum_ir_raman.py -h
python spectrum_uvvis_ecd.py -h
```

> [!NOTE]
> 脚本按真实输出格式编写，但不同软件版本的输出可能略有差异；使用前建议先在小体系上核对结果。

[返回目录](../README.md)
