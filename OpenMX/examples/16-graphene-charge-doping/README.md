# 示例 16-3：石墨烯电荷掺杂（背景电荷 vs ESM on3 栅极模型）

配套教程：[第16章 表面与界面：ESM、功函数、电荷掺杂](../../tutorials/16-surface-esm-workfunction.md)

## 目的

用 `scf.system.charge` 给单层石墨烯掺入 0.02 个电子/原胞，分别采用两种补偿方式：

1. **均匀背景电荷**（不开 ESM）：OpenMX 用 FFT 解 Poisson 方程，多出来的电荷由符号相反的均匀背景中和（“jellium” 背景），见手册第 31 章；
2. **ESM on3**：z=0 一侧是真空，z=20 Å 一侧是理想金属（相当于栅极），补偿电荷是金属里的镜像电荷。

然后比较狄拉克点相对化学势的位置 $E_D-\mu$，并与线性色散估算值对照。另外，这个示例也演示了 `ESM.direction z` 的写法：第一条晶格矢量（a 轴）就是沿 z 的法向。

## 文件说明

| 文件 | System.Name | 说明 |
|---|---|---|
| `Graphene_q0.dat` | `graphene_q0` | 中性参考，不开 ESM |
| `Graphene_qm0.02_bg.dat` | `graphene_qm002_bg` | `scf.system.charge -0.02`（负号表示电子掺杂），背景电荷补偿 |
| `Graphene_qm0.02_esm_on3.dat` | `graphene_qm002_on3` | 同样掺 −0.02，`ESM.switch on3` + `ESM.direction z` |
| `dirac_point.py` | — | 读取 `.Band`，输出 K 点处的 $E_D-\mu$（eV） |

三个输入几何完全相同：晶格常数 2.467 Å，a=(0,0,20) Å，b、c 在 xy 平面且夹角 120°，K 点为 (0, 1/3, 1/3)。能带路径为 Γ–M–K–Γ。

## 运行步骤

```bash
cd OpenMX/examples/16-graphene-charge-doping
sed -i "s#^DATA.PATH.*#DATA.PATH   $HOME/openmx4.0/DFT_DATA19#" *.dat
for f in Graphene_q0.dat Graphene_qm0.02_bg.dat Graphene_qm0.02_esm_on3.dat; do
  n=$(awk 'tolower($1)=="system.name"{print $2}' $f)
  mpirun -np 4 openmx $f -nt 1 > $n.std 2>&1
  grep -q "The calculation was normally finished." $n.std || echo "$f FAILED"
done
python3 dirac_point.py graphene_q0.Band graphene_qm002_bg.Band graphene_qm002_on3.Band
# band plots: bandgnu13 x.Band, or the section script: python3 ../../scripts/omx_band_plot.py -h
```

## 预期结果 / 如何检查

- `graphene_q0`：$E_D-\mu\approx0$，K 点“带隙”接近 0（数值上 meV 量级）。
- 两个掺杂计算：$E_D-\mu<0$，即费米能级进入了导带。用线性色散 $E_F-E_D=\hbar v_F\sqrt{\pi n}$ 粗估：原胞面积 5.27 Å$^2$，$n=0.02/5.27$ Å$^{-2}\approx3.8\times10^{13}$ cm$^{-2}$，取 $v_F\approx(0.85\text{–}1.0)\times10^6$ m/s，得到 $|E_D-\mu|\approx0.6\text{–}0.7$ eV（这是估算值，不是计算结果）。
- 背景电荷和 ESM 两种做法下，能带相对 μ 的位移应当接近，但静电势分布完全不同。用 `gcube2oned graphene_qm002_on3.v0.cube 1` 可以看到：ESM on3 下，石墨烯与金属之间的真空区出现线性斜坡，也就是栅极电场；背景电荷做法里则没有明确的“电极”，势在真空区呈抛物线形弯曲。
- 用 `grep "Boundary condition" graphene_qm002_on3.std` 确认输出为 `Vacuum|Vacuum|Metal`。

## 注意事项

- `on1` 不允许带电体系，带电二维材料请用 `on2`/`on3`，或者不开 ESM、改用背景电荷。
- 背景电荷做法里，总能量包含额外的背景静电项，不同电荷态之间的总能量不能直接比较（手册第 31 章）。
- K 点坐标与 b、c 的夹角有关：b、c 夹角为 120° 时 K=(1/3,1/3)，为 60° 时 K=(1/3,2/3)。换了晶胞写法后，记得重新核对路径。
