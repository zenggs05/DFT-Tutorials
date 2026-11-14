# 第7章 过渡态搜索与 IRC 验证（Transition State Search & IRC with Gaussian 16）

本章讲解用 Gaussian 16 搜索过渡态（transition state, TS）并用 IRC（intrinsic reaction coordinate）验证反应路径的完整流程，覆盖 `Opt=(TS,CalcFC,NoEigenTest)`、QST2/QST3、频率验证单一虚频、`IRC=(RCFC,...)` 双向积分、IRC 端点再优化和活化自由能 ΔG‡ 计算。示例为 Cl⁻ + CH₃Br 的 S<sub>N</sub>2 反应和 HCN ⇌ HNC 异构化，所有输入都可以直接提交运行。关键词：Gaussian transition state optimization、QST2、QST3、IRC calculation、imaginary frequency、activation free energy。

---

## 1. 本章目标与适用场景

学完本章，你应当能够：

1. 从经验构型、柔性扫描最高点、QST2/QST3 或约束优化中得到可用的 TS 初猜；
2. 写出并看懂 `Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq`，知道 `CalcFC`、`ReadFC`、`RCFC`、`CalcAll` 分别适用于什么情况；
3. 用频率计算确认 TS **恰好有 1 个虚频**，并在 GaussView 中检查虚频振动是否对应你关心的反应坐标；
4. 用 `IRC=(RCFC,Phase=(...),MaxPoints=N,StepSize=N)` 双向积分反应路径，看懂 `Summary of reaction path following`；
5. 对 IRC 两端点做 `Opt Freq`，确认 TS 连接的是你想要的反应物和产物，并按 ΔG‡ = G(TS) − G(反应物或前驱复合物) 计算能垒和速率常数；
6. 处理 `Wrong number of Negative eigenvalues`、L123 corrector 失败、TS 收敛到极小点或二阶鞍点等常见问题。

**适用场景**：有机反应机理、催化循环中的基元步骤、构象翻转能垒、质子/氢转移、异构化反应等。凡是要在论文里写“能垒”和“速控步”的工作，都离不开本章的流程。

**前置知识**：几何优化（[第5章](05-opt.md)）、频率与热化学（[第6章](06-freq-thermo.md)）。柔性扫描找 TS 初猜的细节放在 [第8章](08-scan.md)，溶剂中的 TS 见 [第9章](09-solvation.md)。

本章配套示例：

| 示例目录 | 内容 |
|---|---|
| [../examples/07-ts-irc-sn2/README.md](../examples/07-ts-irc-sn2/README.md) | Cl⁻ + CH₃Br S<sub>N</sub>2：TS 优化 + 频率 → IRC 双向 → 端点 Opt Freq → 分离反应物 |
| [../examples/07-qst-hcn-hnc/README.md](../examples/07-qst-hcn-hnc/README.md) | HCN ⇌ HNC：QST2 与 QST3 两种写法 |

---

## 2. 原理速览

### 2.1 过渡态是一阶鞍点

在 Born–Oppenheimer 势能面 $E(\mathbf{x})$ 上，稳定点满足梯度为零：

$$
\mathbf{g} = \nabla E(\mathbf{x}) = 0
$$

稳定点的类型由 Hessian 矩阵 $\mathbf{H}_{ij} = \partial^2 E / \partial x_i \partial x_j$（去掉平动和转动后）的本征值决定：

| Hessian 负本征值个数 | 频率计算中的虚频个数 | 稳定点类型 |
|---|---|---|
| 0 | 0 | 极小点（反应物、产物、中间体） |
| 1 | 1 | **一阶鞍点 = 过渡态** |
| ≥ 2 | ≥ 2 | 高阶鞍点（通常没有化学意义，需要继续处理） |

质量加权 Hessian 的本征值 $\lambda$ 与谐振频率的关系为 $\tilde{\nu} = \sqrt{\lambda}/(2\pi c)$。$\lambda<0$ 时频率是虚数，Gaussian 用**负号**打印，例如 `-412.35`，读作 412i cm⁻¹。这个虚频对应的简正模式称为**过渡矢量**（transition vector），它应当就是反应坐标：沿 + 方向走向一侧，沿 − 方向走向另一侧。

### 2.2 Berny 算法如何找鞍点

Gaussian 默认用 Berny 算法在冗余内坐标中优化。找极小点时，近似 Hessian 只要正定就能工作；找鞍点时，算法沿 Hessian 最低本征矢方向**上坡**，沿其余方向**下坡**（TS 优化默认用 RFO 步，Hessian 用 Bofill 公式更新）。

这就是 TS 优化对**初始 Hessian** 特别敏感的原因：默认的经验力场 Hessian 在任何结构下都是正定的，不知道“哪个方向该上坡”。`CalcFC` 在第一步算出真实 Hessian，如果初猜结构附近恰好有一个负本征值，而且它的本征矢就是反应坐标，优化通常会很顺利。

### 2.3 STQN：QST2 / QST3

QST2/QST3 用的是 Schlegel 等发展的 STQN（synchronous transit-guided quasi-Newton）方法。程序先在反应物与产物之间做线性/二次同步过渡（LST/QST）插值，接近 TS 的二次区域后，再切换到准牛顿/本征矢跟踪完成优化。QST2 只需要反应物和产物；QST3 还要你给一个 TS 猜测。两种方法都要求**所有结构中原子顺序完全一致**。

### 2.4 IRC：从 TS 出发的最陡下降路径

IRC 是质量加权笛卡尔坐标下从 TS 出发的最陡下降路径（Fukui 1970/1981）：

$$
\frac{d\mathbf{x}(s)}{ds} = -\frac{\mathbf{g}(\mathbf{x})}{|\mathbf{g}(\mathbf{x})|}
$$

其中 $s$ 是路径弧长（单位 amu<sup>1/2</sup>·Bohr）。起点处梯度为零，所以第一步沿过渡矢量方向迈出，因此 IRC **必须知道起点的 Hessian**（`RCFC` 或 `CalcFC`）。G16 默认使用 Hratchian–Schlegel 的 HPC（Hessian-based predictor–corrector）积分器，精度和效率都比 G03 时代的 GS2 好得多。

### 2.5 过渡态理论与 Eyring 方程

$$
k = \kappa \frac{k_B T}{h} \exp\left(-\frac{\Delta G^{\ddagger}}{RT}\right)
$$

$\kappa$ 为透射系数（常取 1，涉及氢隧穿时需另行校正）。298.15 K、$\kappa=1$ 时，一级反应半衰期与 ΔG‡ 的关系如下（按上式直接计算）：

| ΔG‡ (kcal/mol) | k (s⁻¹) | 半衰期 t<sub>1/2</sub> |
|---|---|---|
| 15 | 6.3×10¹ | 0.011 s |
| 18 | 0.40 | 1.7 s |
| 21 | 2.5×10⁻³ | 约 4.6 min |
| 24 | 1.6×10⁻⁵ | 约 12 h |
| 27 | 1.0×10⁻⁷ | 约 80 d |

这张表常用来判断算出的能垒是否“说得通”。室温下几小时内完成的反应，ΔG‡ 一般在 20–24 kcal/mol 左右。ΔG‡ 每变化 1.36 kcal/mol，速率就差 10 倍。

### 2.6 本章的标准流程

```text
TS 初猜 (经验构型 / 柔性扫描最高点 / QST2 / QST3 / 约束优化)
        |
        v
Opt=(TS,CalcFC,NoEigenTest) Freq   --> 恰好 1 个虚频? 振动模式是反应坐标?
        |                                   | 否: 换初猜 / 沿模式位移 / CalcAll
        v 是
IRC=(RCFC,...) 正反两个方向       --> 两端是否通向预期的反应物和产物?
        |                                   | 否: 这个 TS 属于另一个反应, 重新找
        v 是
两端点 Opt Freq (0 个虚频)  -->  (可选) 高精度单点  -->  ΔG‡, ΔG_r, k(T)
```

---

## 3. 准备 TS 初猜：五种来源

TS 优化成功与否，八成取决于初猜。下面五种来源按“省事程度”排列。

### 3.1 凭经验在 GaussView 中搭建

把正在断裂的键拉长、正在形成的键缩短到“半成键”距离，其余部分保持合理。常见基元反应 TS 中关键距离的**经验范围**（不同体系差别较大，只用于搭初猜）：

| 反应类型 | 关键坐标 | 初猜可取的范围（Å） |
|---|---|---|
| S<sub>N</sub>2（C–卤素） | C···Nu、C···LG | 2.1–2.5，Nu–C–LG 接近 180° |
| 质子/氢原子转移（O/N/C–H···O/N/C） | X···H、H···Y | 1.2–1.4 |
| Diels–Alder | 两根正在形成的 C···C | 2.1–2.4 |
| [1,2]-H 迁移 | C···H（两侧） | 1.2–1.4，H 位于 C–C 键上方 |
| 金属插入 / 还原消除 | M–X、M–Y、X···Y | 参考类似体系的文献结构 |

> [!TIP]
> 搭好初猜后先用 GaussView 的 Clean 功能检查有没有重叠原子，但**不要**对整个结构 Clean，否则拉长的键会被“修正”回平衡值。

### 3.2 柔性扫描的能量最高点（最常用）

沿反应坐标（例如正在形成的键长）做 `Opt=ModRedundant` 柔性扫描，取能量最高的结构作为初猜。这是最稳妥的办法，详见 [第8章](08-scan.md) 的 `sn2_bond_scan_ts_guess.gjf` 示例。

### 3.3 QST2：只给反应物和产物

适合反应物和产物结构明确、反应坐标比较“直”的情况。程序自动插值生成初猜。见第 6 节实战。

### 3.4 QST3：反应物 + 产物 + TS 猜测

QST2 插值出来的初猜不合理时（例如插值穿过原子重叠区域），可以改用 QST3 自己给出 TS 猜测。QST3 比单独 `Opt=TS` 多利用了反应物和产物的信息，能帮程序确定沿哪个方向上坡。

### 3.5 约束优化后再放开

先用 `Opt=ModRedundant` 冻结正在形成和断裂的键（`F`），让其余自由度弛豫，再把这个结构交给 `Opt=(TS,CalcFC)`：

```text
%chk=pre_ts.chk
%nprocshared=8
%mem=8GB
#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ Opt=ModRedundant

constrained pre-optimization: forming/breaking bonds frozen

0 1
... (Cartesian coordinates, atoms 1..N) ...

B 1 7 F
B 2 8 F

```

> [!IMPORTANT]
> 无论用哪种方式得到初猜，ModRedundant、QST2/QST3 和 IRC `Phase` 中引用的原子编号都是**分子说明部分的原子顺序**（从 1 开始，哑原子不计），也就是 GaussView 中显示的编号。在 GaussView 中删除或增加原子后编号会改变，需要重新核对。QST2/QST3 的多个结构必须保持**一一对应的原子顺序**，可以在 GaussView 的原子列表编辑器（Atom List Editor）中检查和调整。

---

## 4. 实战一：用 Opt=(TS,CalcFC) 优化 S<sub>N</sub>2 过渡态

### 4.1 体系与输入

反应 Cl⁻ + CH₃Br → ClCH₃ + Br⁻ 是经典的气相 S<sub>N</sub>2 反应，TS 为 C<sub>3v</sub> 对称的 [Cl···CH₃···Br]⁻。我们手工搭一个初猜：C 在原点，Cl 和 Br 分别位于 z 轴两侧，三个 H 在 xy 平面内（CH₃ 呈平面构型）。

完整输入 [`sn2_cl_ch3br_ts_calcfc.gjf`](../examples/07-ts-irc-sn2/sn2_cl_ch3br_ts_calcfc.gjf)：

```text
%chk=sn2_ts.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq

SN2 TS search for [Cl...CH3...Br]- with Opt=(TS,CalcFC) followed by Freq

-1 1
C     0.000000    0.000000    0.000000
Cl    0.000000    0.000000   -2.320000
Br    0.000000    0.000000    2.470000
H     1.075000    0.000000    0.000000
H    -0.537500    0.930977    0.000000
H    -0.537500   -0.930977    0.000000

```

### 4.2 关键词逐项解释

| 关键词 | 本例取值 | 说明 |
|---|---|---|
| `#p` | — | 详细输出，便于查看每步的收敛表和 Hessian 信息 |
| `B3LYP EmpiricalDispersion=GD3BJ` | — | 教学用的常规组合。正式研究中应先做泛函基准测试，详见 [第3章](03-method-basis.md) |
| `6-31+G(d)` | — | 体系是**阴离子**，重原子需要弥散函数（`+`）。Gaussian 中 6-31G 系列覆盖 H–Kr，Br 可用 |
| `Opt=TS` | — | 用 Berny 算法优化到一阶鞍点 |
| `CalcFC` | — | 在第一步解析计算 Hessian。TS 优化**强烈推荐**加上，代价是多算一次频率 |
| `NoEigenTest` | — | 关闭“Hessian 负本征值个数必须为 1”的检查，避免中途某一步曲率暂时不对就被终止。必须靠最后的 Freq 确认结果 |
| `MaxCycles=100` | 100 | 最大优化步数。默认值为 20 和“2 × 冗余内坐标数”中的较大者；小体系默认值可能不够 |
| `Freq` | — | 优化收敛后自动在同一水平做频率计算，`Opt Freq` 写在一起就是两步复合作业 |
| `-1 1` | — | 电荷 −1，单重态。电子数 = 6(C) + 17(Cl) + 35(Br) + 3(H) + 1 = 62，为偶数，与单重态的奇偶性一致 |

> [!NOTE]
> 关于 `EigenTest`：官方 Opt 页面写的是该检查“默认只对 Z-matrix 或笛卡尔坐标下的 TS 优化开启”，并提醒 `NoEigenTest` 适合计算资源充足的情况（它可能让优化在错误的曲率下走很多冤枉步）。但实际使用默认冗余内坐标做 TS 优化时，也经常遇到 `Wrong number of Negative eigenvalues` 而终止，因此社区的通行做法是 `CalcFC` 与 `NoEigenTest` 一起用：初始 Hessian 是准确的，关闭检查的风险就小了。不管怎样，最后一定要靠频率计算验证。

### 4.3 运行

```bash
cd Gaussian/examples/07-ts-irc-sn2
g16 < sn2_cl_ch3br_ts_calcfc.gjf > sn2_cl_ch3br_ts_calcfc.log
# 或者在集群上提交（脚本里依次运行 TS -> IRC -> 端点优化）
sbatch run_ts_irc.sh
```

`%nprocshared`、`%mem` 要与作业调度系统申请的资源匹配，见 [第17章](17-hpc-jobs.md) 和 [Common 第2章](../../Common/tutorials/02-job-scheduler.md)。也可以用仓库脚本 [`sub_g16.sh`](../scripts/sub_g16.sh) 提交（`bash ../scripts/sub_g16.sh -h` 查看用法）。

### 4.4 监控优化过程

```bash
# 每一步的四项收敛判据
grep -A4 "Converged?" sn2_cl_ch3br_ts_calcfc.log | tail -n 5
# 当前进行到第几步
grep "Step number" sn2_cl_ch3br_ts_calcfc.log | tail -n 1
# 每一步的 SCF 能量
grep "SCF Done" sn2_cl_ch3br_ts_calcfc.log
```

TS 优化的能量**不一定单调**：一开始可能上升（在往鞍点“爬”），之后在鞍点附近振荡收敛。用 `python ../scripts/opt_monitor.py -h` 查看用法，它可以把每一步的能量和四项判据画成曲线。

收敛时的典型输出（示意输出，数值仅供参考）：

```text
         Item               Value     Threshold  Converged?
 Maximum Force            0.000021     0.000450     YES
 RMS     Force            0.000008     0.000300     YES
 Maximum Displacement     0.000390     0.001800     YES
 RMS     Displacement     0.000165     0.001200     YES
 Predicted change in Energy=-2.1D-09
 Optimization completed.
    -- Stationary point found.
```

TS 优化的 `Optimized Parameters` 表与极小点相同；同时第一步 `CalcFC` 之后会打印初始 Hessian 的本征值，`#p` 下能看到有几个负值。

### 4.5 TS 优化的其他常用选项

| 选项 | 作用 | 什么时候用 |
|---|---|---|
| `CalcFC` | 第一步算解析 Hessian | 几乎所有 TS 优化的默认选择 |
| `ReadFC` | 从 chk 读取力常数（已转换为内坐标） | 先在低水平做 `Opt=(TS,CalcFC) Freq`，再在高水平 `Opt=(TS,ReadFC)` |
| `RCFC` | 从 chk 读取**笛卡尔**力常数（只能来自 Freq 作业） | 低水平频率显示结构需要降低对称性、旧内坐标不再适用时 |
| `CalcAll` | 每一步都算解析 Hessian | 极难收敛的 TS（势能面平坦、多个软模式），代价很高；收敛后自动做频率分析 |
| `RecalcFC=N` | 第 1 步及之后每 N 步重算 Hessian | `CalcAll` 太贵、`CalcFC` 又不够时的折中 |
| `MaxStep=N` | 最大步长（信赖半径）为 0.01N Bohr 或弧度，默认 N=30 | 优化振荡、步子迈太大时减小到 10–20 |
| `Tight` / `VeryTight` | 收紧收敛判据 | 势能面非常平坦、频率对结构敏感时；DFT 下同时保证积分格点足够（G16 默认 UltraFine） |
| `Restart` | 从 chk 续算中断的优化 | 续算时 route 中要重写 `Opt=(TS,Restart,...)`；**必须重新加 `CalcFC` 才会在续算第一步计算 Hessian** |
| `ModRedundant` | 增删或冻结冗余内坐标 | 例如把某个关键距离加入内坐标，以便在输出中跟踪 |

**两级策略：低水平 Hessian + 高水平 TS**

体系较大时，在高水平上用 `CalcFC` 很贵。可以先在低水平（小基组，甚至半经验或 xtb）上把 TS 找好并算出频率，再把 Hessian 传给高水平：

```text
%chk=ts_low.chk
%nprocshared=16
%mem=16GB
#p B3LYP/6-31G(d) Opt=(TS,CalcFC,NoEigenTest) Freq

low-level TS search and Hessian

0 1
... (initial TS guess) ...

--Link1--
%oldchk=ts_low.chk
%chk=ts_high.chk
%nprocshared=16
%mem=16GB
#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Opt=(TS,ReadFC,NoEigenTest) Freq Geom=Check Guess=Read

high-level TS starting from the low-level geometry and force constants

0 1

```

`%oldchk` 指定要读取的旧 chk，`%chk` 指定新写的 chk，这样低水平结果不会被覆盖（详见 [第2章](02-input-anatomy.md)）。

---

## 5. 验证：恰好 1 个虚频，且对应反应坐标

### 5.1 数虚频

`Opt Freq` 结束后：

```bash
grep -c "Normal termination" sn2_cl_ch3br_ts_calcfc.log        # Opt Freq 两步都正常结束应为 2
grep "imaginary frequencies" sn2_cl_ch3br_ts_calcfc.log        # 有虚频才会出现这一行
grep "Frequencies --" sn2_cl_ch3br_ts_calcfc.log | head -n 2   # 最低的几个频率
grep "Low frequencies" sn2_cl_ch3br_ts_calcfc.log
```

期望看到（示意输出，数值仅供参考）：

```text
 ******    1 imaginary frequencies (negative Signs) ******
 Frequencies --   -385.1234               183.4567               183.4567
```

也可以用仓库脚本批量检查终止状态和虚频：`python ../scripts/g16_log_summary.py -h` 查看用法。

### 5.2 看虚频振动是不是你要的反应

**只有 1 个虚频还不够**，还要确认这个模式确实是“Cl 靠近 C、Br 离开 C、CH₃ 伞形翻转”。两种方法：

1. **GaussView**：打开 `.log`（或 `formchk` 生成的 `.fchk`），Results → Vibrations，选中负频率，点 Start Animation。也可以用对话框里的手动位移（Manual Displacement）把结构沿模式推到一侧，另存为新结构。
2. **看日志中的位移矢量**：频率表下方列出每个原子的 X/Y/Z 位移。对本例，第一个模式中 C 和三个 H 沿 z 轴的分量应该最大，Cl 和 Br 则是反向的小位移。

**判断标准**：

| 现象 | 结论与对策 |
|---|---|
| 1 个虚频，模式明显对应成键/断键 | ✅ 候选 TS，进入 IRC 验证 |
| 1 个虚频，但只是甲基转动、苯环翻转之类的软模式（常 < 100i cm⁻¹） | ❌ 找到的是构象 TS 或数值噪声。重新构建初猜 |
| 0 个虚频 | ❌ 优化落到了极小点。初猜离鞍点太远，或 Hessian 曲率不对。用更好的初猜（扫描最高点）或 QST3 |
| 2 个及以上虚频 | ❌ 高阶鞍点。若多出来的是小虚频（转子），沿该模式位移后重新 `Opt=(TS,CalcFC)`；必要时 `CalcAll` |
| 虚频数正确但数值很小（几十 i）且模式正确 | 可能是很平的势能面，例如弱复合物内的转动，或格点噪声。用 `Opt=(TS,CalcFC,Tight)` 并配合 `Int=SuperFine` 复核 |

> [!TIP]
> 典型虚频量级（经验值）：涉及 H 转移的 TS 常在 1000–2000i cm⁻¹，重原子 S<sub>N</sub>2 在几百 i cm⁻¹，构象翻转常在 100i cm⁻¹ 以下。虚频大小本身不能证明什么，只用来辅助判断。

> [!WARNING]
> 频率计算必须与优化在**完全相同**的方法、基组、色散校正、溶剂模型和积分格点下进行。`Opt Freq` 写在同一个 route 中可以自动保证这一点；分开算时最容易出错（例如优化时加了 `SCRF`，频率时忘了加）。

---

## 6. 实战二：QST2 与 QST3（HCN ⇌ HNC）

### 6.1 QST2 输入

HCN → HNC 是 H 原子从 C 迁移到 N 的 1,2-迁移，非常适合演示 QST2。完整输入 [`hcn_hnc_qst2.gjf`](../examples/07-qst-hcn-hnc/hcn_hnc_qst2.gjf)：

```text
%chk=hcn_hnc_qst2.chk
%nprocshared=4
%mem=4GB
#p B3LYP/6-311G(d,p) Opt=(QST2,CalcFC) Freq

HCN reactant, atom order C N H, H bent 10 deg off the axis to define the plane

0 1
C     0.000000    0.000000    0.000000
N     0.000000    0.000000    1.156000
H     0.185000    0.000000   -1.049000

HNC product, same atom order C N H, H bent to the same side

0 1
C     0.000000    0.000000    0.000000
N     0.000000    0.000000    1.169000
H     0.172000    0.000000    2.154000

```

**输入结构要点**：

| 要点 | 说明 |
|---|---|
| route 中写 `Opt=QST2`，**不要**再加 `TS` | 官方文档明确说明 TS 选项不能与 QST2/QST3 同时使用 |
| 每个结构前面都有自己的**标题行 + 空行** | 结构 1 是反应物，结构 2 是产物，顺序固定 |
| 每个结构都要有电荷与多重度行 | 两者必须一致 |
| 原子顺序完全一致 | 两个结构中第 3 个原子都是那个迁移的 H。顺序错位会得到毫无意义的插值 |
| 与 `ModRedundant` 联用时 | **每个**结构后面都要跟一段 ModRedundant 输入 |
| `CalcFC` | 在插值得到的初猜上计算 Hessian，提高成功率 |

> [!NOTE]
> 为什么把 H 偏离分子轴 10°？HCN 和 HNC 都是线性分子。如果两端都严格线性，H 偏离分子轴的方向和插值所在的平面都无法确定，插值结果可能不合理，还容易触发线性角相关的报错（见第 10 节）。让两端的 H 朝同一侧稍微弯曲，就能给插值一个明确的平面。QST2 不要求输入的反应物和产物是严格的极小点，但它们应当接近真实结构。

### 6.2 QST3 输入

QST3 在 QST2 的基础上增加第三个结构作为 TS 猜测（顺序：反应物、产物、TS 猜测）。完整输入见 [`hcn_hnc_qst3.gjf`](../examples/07-qst-hcn-hnc/hcn_hnc_qst3.gjf)，第三段是：

```text
TS guess, H bridging C and N (H-C-N about 72 deg), same atom order C N H

0 1
C     0.000000    0.000000    0.000000
N     0.000000    0.000000    1.180000
H     1.132000    0.000000    0.368000

```

这个 TS 猜测按“C–H ≈ 1.19 Å、∠HCN ≈ 72°”手工搭建，只是一个近似初猜。

### 6.3 结果解读

QST2/QST3 收敛后，`Optimized Parameters` 表比普通优化多出 **Reactant** 和 **Product** 两列，可以直接对比 TS 中各内坐标相对两端的位置（示意输出，数值仅供参考）：

```text
 ! Name  Definition              Value    Reactant  Product   Derivative Info.  !
 ! R1    R(1,2)                  1.1890    1.156     1.169    -DE/DX =    0.0   !
 ! R2    R(1,3)                  1.1950    1.065     2.154    -DE/DX =    0.0   !
```

后续验证步骤完全相同：`Freq` 中恰好 1 个虚频，模式是 H 在 C 与 N 之间摆动，再用 IRC 确认两端。

### 6.4 QST2/QST3 与 Opt=TS 怎么选

| 情况 | 推荐 |
|---|---|
| 有可靠的 TS 初猜（扫描最高点或类似体系的文献结构） | `Opt=(TS,CalcFC,NoEigenTest)`，最直接 |
| 只有反应物和产物，反应坐标简单 | QST2 |
| QST2 插值结构不合理，或你大致知道 TS 的样子 | QST3 |
| 反应物/产物是分开的两个分子（双分子反应） | 先构建**复合物**作为反应物或产物，再用 QST2/QST3；直接用相距很远的分子插值通常效果很差 |
| 想同时得到整条路径 | `Opt=(QST2,Path=M)`（同时优化 TS 和一条 M 点路径，较少使用），或扫描 + IRC |

---

## 7. 实战三：IRC 验证反应路径

### 7.1 为什么必须做 IRC

1 个虚频只能说明这是一个一阶鞍点，**不能**说明它连接你设想的反应物和产物。复杂体系中，经常出现“看起来对，IRC 走下去却是另一个反应”的情况，例如 TS 实际对应构象翻转、质子在另一对原子间转移，或者协同/分步机理判断错误。审稿人通常会要求提供 IRC 证据。

### 7.2 输入

TS 的 `Opt Freq` 作业已经把笛卡尔力常数写进了 `sn2_ts.chk`，IRC 直接读取即可。完整输入 [`sn2_irc_both_directions.gjf`](../examples/07-ts-irc-sn2/sn2_irc_both_directions.gjf)：

```text
%oldchk=sn2_ts.chk
%chk=sn2_irc.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ IRC=(RCFC,Phase=(1,3),MaxPoints=80,StepSize=10,MaxCycle=50) Geom=Check Guess=Read

SN2 IRC in both directions from the converged TS; forward = C-Br bond lengthening

-1 1

```

### 7.3 IRC 关键词逐项解释

| 选项 | 本例 | 默认值 | 说明 |
|---|---|---|---|
| `RCFC` | ✓ | — | 从 chk 读取 Freq 作业算出的笛卡尔力常数（同义词 `ReadCartesianFC`）。**`RCFC` 和 `CalcFC` 必须二选一**（也可用 `CalcAll`） |
| `CalcFC` | — | — | 在起点重新计算 Hessian。没有现成 chk 时使用 |
| `CalcAll` | — | — | 每个点都算 Hessian，HPC 积分器下一般不需要 |
| `Phase=(N1,N2[,N3[,N4]])` | `(1,3)` | 由过渡矢量最大分量的符号决定 | 定义“正向”：沿正向走时指定的键长/键角/二面角**增大**。本例中正向就是 C(1)–Br(3) 变长，即走向产物 |
| `Forward` / `Reverse` | — | 两个方向都走 | 只计算一个方向，可以把两个方向拆成两个作业并行跑 |
| `MaxPoints=N` | 80 | 10 | **每个方向**最多走的点数。默认 10 点往往走不到底，重原子体系尤其如此 |
| `StepSize=N` | 10 | 10 | 沿路径的步长。官方定义的单位是 0.01 Bohr；N<0 时改按 0.01 amu<sup>1/2</sup>·Bohr 计（取 \|N\|）。路径很长时可增大到 15–20 以减少点数 |
| `MaxCycle=N` | 50 | 20 | 每个点上校正步（几何优化）的最大步数 |
| `HPC` | 默认 | ✓ | Hessian-based predictor–corrector，G16 中大多数方法的默认积分器 |
| `EulerPC` | — | `GradientOnly` 和 ONIOM(MO:MM) 时默认 | 一阶 Euler 预测 + HPC 校正，适合大体系 |
| `LQA` | — | — | 局部二次近似预测步。HPC 校正失败（L123）时可尝试 |
| `DVV` / `Euler` | — | — | 阻尼速度 Verlet / 纯 Euler（后者不推荐用于正式计算） |
| `GS2` | — | — | G03 及更早版本的默认算法（Gonzalez–Schlegel），可用于复现旧结果，或在其他积分器失败时作为备选 |
| `ReCalc=N` | — | — | 每 N 个预测步重算解析 Hessian，需要方法支持解析二阶导数 |
| `GradientOnly` | — | — | 方法没有解析二阶导时**必须手动指定** |
| `Tight` / `VeryTight` | — | — | 收紧每个点的收敛阈值，DFT 下同时保证足够的积分格点 |
| `Report[=Bonds/Angles/Dihedrals/Cartesians/Read]` | — | 不报告 | 在每个点打印几何参数，便于跟踪关键键长 |
| `Restart` | — | — | 续算未完成的 IRC，或在已完成的 IRC 上继续多走几步 |
| `Downhill` | — | — | 从输入结构直接下坡（结构不是 TS 时用） |
| `MassWeighted` / `Cartesian` | — | 质量加权 | 坐标系选择，一般保持默认 |

> [!NOTE]
> 不少中文教程写 `IRC=(ReadFC,...)`。G16 官方 IRC 页面列出的读取选项是 `RCFC`（以及同义词 `ReadCartesianFC`），本教程统一使用 `RCFC`。`ReadFC` 是 `Opt` 关键词的选项。

`Geom=Check Guess=Read` 从 chk 读取 TS 结构和波函数。由于用了 `%oldchk`，原来的 `sn2_ts.chk` 不会被改写。

### 7.4 运行与监控

```bash
g16 < sn2_irc_both_directions.gjf > sn2_irc_both_directions.log
grep "Point Number" sn2_irc_both_directions.log | tail -n 3
grep "NET REACTION COORDINATE" sn2_irc_both_directions.log | tail -n 3
```

每走一个点，日志中都会出现类似下面的信息（示意输出，数值仅供参考）：

```text
 Point Number:  25          Path Number:   1
   CHANGE IN THE REACTION COORDINATE =    0.10000 (Bohr/amu^1/2)
   NET REACTION COORDINATE UP TO THIS POINT =    2.50000 (Bohr/amu^1/2)
```

`Path Number: 1` 为正向，`2` 为反向。全部完成后会打印汇总表：

```bash
grep -A 200 "Summary of reaction path following" sn2_irc_both_directions.log | head -n 60
```

汇总表格式（示意输出，数值仅供参考）：

```text
 Energies reported relative to the TS energy of     -3074.512345
 --------------------------------------------------------------------------
    Summary of reaction path following
 --------------------------------------------------------------------------
                        Energy    RxCoord
     1                -0.01452   -7.89012
   ...
    80                 0.00000    0.00000      transition state
   ...
   161                -0.02105    8.01234
 --------------------------------------------------------------------------
```

TS 位于表的中间，能量和反应坐标都是 0。反应坐标一侧为负、一侧为正，能量从 0 单调下降。

### 7.5 作图与判断

- **GaussView**：打开 IRC 的 `.log`，Results → IRC/Path，会画出能量–反应坐标曲线，点击曲线上的点可以查看对应结构，并能播放整条路径动画。
- **仓库脚本**：[`scan_irc_profile.py`](../scripts/scan_irc_profile.py)（`python ../scripts/scan_irc_profile.py -h` 查看用法）。它会合并正反两个方向、以 kcal/mol 为单位作图。

判断 IRC 是否成功：

1. 两个方向的能量都从 TS 处单调下降，没有回升；
2. 正向终点结构接近产物（或产物复合物），反向终点接近反应物（或反应物复合物）。对本例，就是看 C–Cl 与 C–Br 距离在两端是否如期变化，可以在 GaussView 中测量，或用 `python ../scripts/log2xyz.py -h` 把路径导出为多帧 xyz 后查看；
3. 如果两个方向都回到同一个极小点，或者走向完全无关的结构，说明这个 TS 不属于你要研究的反应。

> [!TIP]
> IRC 走到 `MaxPoints` 就会停下，此时终点**还不是极小点**，这很正常，下一节会把端点结构再优化一次。如果路径还远没走出 TS 区域，可以用 `IRC=(Restart,MaxPoints=...)` 从原 IRC 作业的 chk 续算（`%chk` 指向 `sn2_irc.chk`），或者加大 `StepSize` 重新计算。

---

## 8. IRC 端点再优化与能垒计算

### 8.1 端点 Opt Freq

从 IRC 两端取出最后一个点的结构（GaussView 中在 IRC 曲线上点击端点，然后 File → Save；或用 `log2xyz.py`），分别做普通的 `Opt Freq`。完整输入 [`irc_endpoint_opt_freq.gjf`](../examples/07-ts-irc-sn2/irc_endpoint_opt_freq.gjf) 用 `--Link1--` 串联了反应物复合物和产物复合物两个作业，第一个作业为：

```text
%chk=sn2_rc.chk
%nprocshared=8
%mem=8GB
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ Opt Freq

Reactant-side IRC endpoint (Cl- ... CH3Br ion-dipole complex): Opt Freq

-1 1
C     0.000000    0.000000    0.000000
Cl    0.000000    0.000000   -3.200000
Br    0.000000    0.000000    1.990000
H     1.030000    0.000000   -0.325000
H    -0.515000    0.892006   -0.325000
H    -0.515000   -0.892006   -0.325000

```

示例文件中的端点坐标是按离子–偶极复合物手工搭建的近似结构，**实际使用时请替换为你自己 IRC 的最后一个点**。端点 `Opt Freq` 必须满足：

- 正常结束，`imaginary frequencies` 行不出现（0 个虚频）；
- 优化后的结构与 IRC 端点属于同一个能量谷，没有“跑”到别的构象。

> [!NOTE]
> G16 C.01 新增了 `Opt=NGoDown=N` 和 `MaxEStep=N`。最小化优化走到鞍点附近、Hessian 出现负本征值时，程序最多混合 N 个负本征矢（默认 3），沿负曲率方向迈出一步离开鞍点；步长为 N/1000 Bohr 或弧度（普通优化默认 600，即 0.6）。`NGoDown=-1` 关闭此功能。从“带小虚频的结构”出发重新优化时，配合 `Opt=CalcFC` 让第一步就有准确的 Hessian，效果更好。

### 8.2 能垒怎么算

从各个 `Opt Freq` 日志中取 `Sum of electronic and thermal Free Energies`：

```bash
grep "Sum of electronic and thermal Free Energies" sn2_cl_ch3br_ts_calcfc.log irc_endpoint_opt_freq.log
```

$$
\Delta G^{\ddagger} = G(\mathrm{TS}) - G(\mathrm{reference}),\qquad
\Delta G_r = G(\mathrm{product}) - G(\mathrm{reactant})
$$

单位换算：1 Hartree = 627.509 kcal/mol = 2625.50 kJ/mol。

**参考态的选择**是最容易出错的地方：

| 情况 | ΔG‡ 的参考态 |
|---|---|
| 单分子反应（异构化、重排） | 反应物极小点 |
| 双分子反应，存在稳定的前驱复合物（如 S<sub>N</sub>2 的离子–偶极复合物） | 取**前驱复合物与分离反应物中自由能更低的那个**。气相离子–分子复合物的能量往往低于分离反应物；溶液中复合物通常并不稳定，一般用分离反应物 |
| 多步反应 | 用能量跨度（energetic span）模型：最低中间体与其后最高 TS 之差，而不是简单取相邻两点 |

示例目录中的 [`sn2_separated_reactants.gjf`](../examples/07-ts-irc-sn2/sn2_separated_reactants.gjf) 分别计算了 CH₃Br 和 Cl⁻，用来得到“分离反应物”这个参考态。

**手算示例**（示意数值，不是真实计算结果）：

| 物种 | G (Hartree) |
|---|---|
| 反应物复合物 Cl⁻···CH₃Br | −3074.520000 |
| TS | −3074.506500 |

ΔG‡ = (−3074.506500) − (−3074.520000) = 0.013500 Hartree × 627.509 = **8.47 kcal/mol**。

### 8.3 不可忽略的修正

1. **标准态**：Gaussian 热化学对应 1 atm 理想气体。溶液中的反应要换算到 1 mol/L，每个物种加 RT ln(24.46) = 1.89 kcal/mol（298.15 K）。双分子反应 A + B → TS 时，分子数减少 1，所以 ΔG‡(1 M) = ΔG‡(1 atm) − 1.89 kcal/mol。详见 [第6章](06-freq-thermo.md)。
2. **低频模式**：复合物和 TS 中常有几十 cm⁻¹ 的低频振动，RRHO 会高估熵，可用 Grimme qRRHO 校正（GoodVibes、Shermo，或仓库的 `thermo_qrrho.py`，`python ../scripts/thermo_qrrho.py -h` 查看用法）。
3. **高精度单点**：常用 G = E<sub>高水平单点</sub> + G<sub>热校正，低水平</sub>，热校正取自 `Thermal correction to Gibbs Free Energy` 行。所有物种必须使用同一组合，记号写成 `ωB97X-D/def2-TZVP//B3LYP-D3(BJ)/6-31+G(d)`。
4. **溶剂**：在溶液中研究反应时，TS 和极小点都应在同一 SCRF 下优化和计算频率，或者统一采用“气相结构 + 溶剂单点”的方案，详见 [第9章](09-solvation.md)。
5. **隧穿**：氢转移反应在低温下隧穿效应明显，可用 Wigner 校正 $\kappa = 1 + \frac{1}{24}\left(\frac{h\,\mathrm{Im}(\nu^\ddagger)}{k_B T}\right)^2$ 做粗略估计。

---

## 9. 廉价替代与进阶工具

| 工具 | 思路 | 适用 |
|---|---|---|
| [xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian) | 通过 Gaussian 的 `External` 关键词调用 xtb 的能量、梯度和 Hessian，用 Gaussian 自己的优化器在 GFN2-xTB 水平上跑 `Opt=(TS,CalcFC,NoEigenTest)`、QST2 甚至 IRC | 大体系 TS 预搜索：几分钟得到初猜，再交给 DFT |
| [pyQRC](https://github.com/patonlab/pyQRC) | 读取频率结果，沿某个（虚频）模式正反位移，生成新的输入 | “廉价 IRC”：沿虚频两侧位移后各做一次 Opt，看落到哪里；也可以用来消除残留小虚频 |
| [dmf-g16](https://github.com/shin1koda/dmf-g16) | 直接 MaxFlux 双端搜索，由 Gaussian 提供能量和梯度 | QST2/QST3 失败时的现代双端替代方案 |
| [autodE](https://github.com/duartegroup/autodE) | 从 SMILES 自动生成反应剖面（构象、TS 猜测、TS 优化、验证），后端支持 Gaussian | 批量反应、自动化机理研究 |
| [IRC_to_com](https://github.com/jvalegre/IRC_to_com) | 把 IRC 路径上的点（含两端）转成 Gaussian 输入 | 端点优化或沿路径做高精度单点 |
| [gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot) | 中文项目，包含 TS/IRC 结果的自查清单 | 对照检查虚频、IRC 端点等 |

> [!CAUTION]
> “沿虚频位移后优化”只能**部分**替代 IRC：它告诉你 TS 两侧最近的极小点是什么，但不能保证反应路径上没有其他中间体或分岔。正式发表的工作仍建议做完整 IRC。

---

## 10. 溶剂中的过渡态

在溶液中优化 TS 时，所有步骤都加同样的 `SCRF`：

```text
#p B3LYP/6-31+G(d) EmpiricalDispersion=GD3BJ SCRF=(SMD,Solvent=Water) Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100) Freq
```

IRC 和端点优化也用同样的 `SCRF=(SMD,Solvent=Water)`。注意：

- 带电体系（如本例阴离子）在溶液中的势能面和气相差别很大。气相中稳定的离子–偶极复合物在水中往往不再是明显的极小点，此时 ΔG‡ 应相对分离的溶剂化反应物计算；
- PCM 类模型在 TS 附近有时会给出小的数值噪声虚频，可用 `Int=SuperFine` 复核；
- 溶剂模型的选择和溶剂化自由能的计算见 [第9章](09-solvation.md)。

---

## 参数选择与收敛建议

| 项目 | 推荐设置 | 理由 |
|---|---|---|
| TS 优化 | `Opt=(TS,CalcFC,NoEigenTest,MaxCycles=100)` | 准确的初始 Hessian 是成功的关键；关闭中途曲率检查，最后由 Freq 把关 |
| 大体系 TS | 低水平 `Opt=(TS,CalcFC) Freq` → 高水平 `Opt=(TS,ReadFC)` | 避免在高水平上算 Hessian |
| 难收敛 TS | `RecalcFC=5` 或 `CalcAll`，`MaxStep=10–20` | 势能面平坦或高度非谐时，Hessian 更新不可靠 |
| 频率 | `Opt Freq` 写在同一 route | 保证方法、基组、溶剂、格点完全一致 |
| 积分格点 | G16 默认 `UltraFine`；出现可疑小虚频时用 `Int=SuperFine` 复核 | Minnesota 泛函对格点敏感 |
| IRC 力常数 | `RCFC`（读 TS 频率的 chk） | 不必重复计算 Hessian |
| IRC 点数 | `MaxPoints=50–100`，`StepSize=10`（重原子体系可取 15–20） | 默认 10 点通常走不到底 |
| IRC 积分器 | 默认 HPC；失败时依次尝试 `LQA`、`EulerPC`、`GS2`，同时减小 `StepSize` | 不同积分器对病态区域的稳定性不同 |
| 每点校正 | `MaxCycle=50` | 默认 20，在平坦区域可能不够 |
| 端点 | 一律再做 `Opt Freq`（0 个虚频） | IRC 终点不是极小点 |
| 能垒 | 同一理论水平；1 M 标准态；必要时 qRRHO；注明参考态 | 保证不同物种之间可比 |

---

## 常见问题与报错

### 1. `Wrong number of Negative eigenvalues: Desired= 1 Actual= 0`

**位置**：L103（Berny 优化），出现在 TS 优化中途或第一步之后。

**原因**：当前 Hessian 没有负本征值（或负本征值个数 ≠ 1），程序认为结构不在鞍点区域，曲率不对。

**处理**：

1. 加上 `CalcFC`，让第一步使用真实 Hessian，而不是经验 Hessian；
2. 加 `NoEigenTest` 允许程序继续优化，但最后一定要用 Freq 验证；
3. 换一个更好的初猜：柔性扫描最高点（[第8章](08-scan.md)）或 QST3；
4. 如果 `CalcFC` 后 Hessian 仍然没有负本征值，说明初猜离鞍点太远，单靠 `NoEigenTest` 多半会滑到极小点。

注意 `Opt=RCFC` 不能用于“因为负本征值个数不对而失败”的优化之后（chk 里没有笛卡尔力常数）。此时应从最新结构出发重新 `CalcFC`。

### 2. TS 优化“收敛”了，但频率是 0 个或 2 个虚频

- **0 个**：落到了极小点。TS 初猜要更靠近鞍点，比如取扫描最高点，或把关键距离设在扫描得到的最大值附近；
- **2 个**：多出来的通常是甲基、苯环、羟基等转子的转动模式。GaussView 中沿多出的那个模式位移（Manual Displacement），保存结构，再 `Opt=(TS,CalcFC,NoEigenTest)`；或用 pyQRC 生成位移结构。

### 3. `Optimization stopped. -- Number of steps exceeded`（L9999）

TS 优化达到 `MaxCycles` 仍未收敛。先用 `opt_monitor.py` 或 GaussView 看能量和力的走势：

- 稳定下降、只是步数不够：`Opt=(TS,Restart,CalcFC,MaxCycles=200)` 从 chk 续算（续算时必须重新写 `CalcFC`）；
- 来回振荡：减小 `MaxStep`，用 `RecalcFC=5` 或 `CalcAll`；
- 结构已经明显偏离预想的反应：换初猜。

更多 L9999 的诊断方法见 [第5章](05-opt.md)。也可以用 [`g16_error_doctor.py`](../scripts/g16_error_doctor.py) 自动识别报错类型（`python ../scripts/g16_error_doctor.py -h` 查看用法）。

### 4. IRC 报 `Max corrector steps exceded`（L123）

日志中拼写确实是 `exceded`。**原因**：HPC 积分器的校正步在某个点上无法收敛，常见于势能面陡峭、弯曲剧烈或接近极小点的平坦区域。

**处理**（按顺序尝试）：

1. 确认 TS 本身已经充分收敛、恰好 1 个虚频（Freq 与 IRC 在同一水平）；
2. 增大每点的校正步数，例如 `MaxCycle=100`；
3. 换积分器：`IRC=(RCFC,LQA,...)`；仍不行就用 `GS2` 并减小步长，例如 `StepSize=5`；
4. 如果已经走了很多点、离 TS 很远，可以直接取已有的最后一个点做端点优化，不必强求 IRC 走完。

### 5. IRC 只走几步就停了，或两个方向都回到同一个结构

- 只走几步：`MaxPoints` 太小（默认 10），或者已经到达极小点而提前停止。检查日志末尾的说明信息；
- 两个方向回到同一个极小点：这个 TS 对应的是对称的过程（例如甲基转动、简并重排），或者 TS 根本不在你设想的反应路径上；
- 正向/反向与预期相反：只是方向定义问题，用 `Phase=(N1,N2)` 明确指定。

### 6. QST2 插值得到的初猜很怪，或报 FormBX / Linear angle 相关错误

- 检查两个结构的**原子顺序**是否一致，这是最常见的原因；
- 线性或接近线性的结构（如 HCN、S<sub>N</sub>2 中的 Nu–C–LG）会让冗余内坐标插值出问题。可以给两端稍微弯曲的结构（如本章 HCN 示例）、改用 QST3 显式给出 TS 猜测，或者先扫描再 `Opt=TS`；
- `FormBX had a problem` / `Linear angle in Bend` 的通用处理见 [第8章](08-scan.md) 和 [第18章](18-errors-capstone.md)。

### 7. IRC 报读不到力常数

`RCFC` 需要 chk 中存在 Freq 作业写入的笛卡尔力常数。常见原因：%chk 名字写错；用了 `%oldchk` 却没有同时写 `%chk`；TS 的 Freq 是在别的 chk 里算的；或者 chk 来自其他机器/版本。解决办法：核对文件名，或者改用 `IRC=(CalcFC,...)` 在起点重新算一次 Hessian。跨机器使用 chk 时应先 `formchk` 再 `unfchk`（见 [第16章](16-postprocess-multiwfn.md)）。

---

## 练习

1. 把 S<sub>N</sub>2 示例的 Cl 换成 F，重新做 TS + IRC，比较两种亲核试剂的 ΔG‡ 和虚频大小。
2. 对 HCN ⇌ HNC 分别运行 QST2、QST3，以及在 QST3 的 TS 猜测上直接做 `Opt=(TS,CalcFC)`，比较三者的优化步数和最终 TS 是否一致。
3. 对同一个 TS 分别用 `IRC=(RCFC)`（默认 HPC）、`IRC=(RCFC,LQA)` 和 `IRC=(RCFC,GS2)` 计算 IRC，比较路径形状和耗时。
4. 用 [第8章](08-scan.md) 的 `sn2_bond_scan_ts_guess.gjf` 扫描得到 TS 初猜，与本章手工搭建的初猜比较 TS 优化所需步数。
5. 在 SMD 水溶液中重复 S<sub>N</sub>2 的 TS 和端点计算，看离子–偶极复合物是否仍然是明显的极小点，以及 ΔG‡ 相对分离反应物如何变化。

---

## 参考资料

**官方文档**

- [Opt 关键词（TS、QST2/QST3、CalcFC/ReadFC/RCFC/CalcAll、EigenTest、ModRedundant、GIC）](https://gaussian.com/opt/)
- [IRC 关键词（RCFC/CalcFC、MaxPoints、StepSize、HPC/EulerPC/LQA/GS2、Phase、Report）](https://gaussian.com/irc/)
- [QST2 Transition State Optimizations（白皮书）](https://gaussian.com/qst2/)
- [Freq 关键词](https://gaussian.com/freq/)
- [Thermochemistry in Gaussian（白皮书）](https://gaussian.com/thermo/)
- [Gaussian 16 Rev. C.01/C.02 Release Notes（Opt=NGoUp/NGoDown/MaxEStep 等新增功能）](https://gaussian.com/relnotes/)
- [Solutions for Gaussian Error Messages（网页版报错大全）](https://liyuanhe211.github.io/Solution_for_Every_Gaussian_Error_Message/)

**相关 GitHub 仓库**（只链接，不复制代码）

- [patonlab/pyQRC](https://github.com/patonlab/pyQRC)：沿虚频位移生成新输入，“廉价 IRC”和去除小虚频
- [the-matter-lab/xtb-gaussian](https://github.com/the-matter-lab/xtb-gaussian)：通过 `External` 让 Gaussian 调用 xtb，做 GFN2-xTB 水平的 TS 预搜索
- [shin1koda/dmf-g16](https://github.com/shin1koda/dmf-g16)：直接 MaxFlux 双端 TS 搜索
- [duartegroup/autodE](https://github.com/duartegroup/autodE)：从 SMILES 自动生成反应剖面
- [jvalegre/IRC_to_com](https://github.com/jvalegre/IRC_to_com)：IRC 路径点转 Gaussian 输入
- [patonlab/GoodVibes](https://github.com/patonlab/GoodVibes)：qRRHO、标准态校正与反应能垒剖面
- [liyuanhe211/Solution_for_Every_Gaussian_Error_Message](https://github.com/liyuanhe211/Solution_for_Every_Gaussian_Error_Message)：按 Link 编号整理的报错大全（CC BY-NC-SA，仅链接）
- [YuOZW/gaussian-tips](https://github.com/YuOZW/gaussian-tips)：日文技巧合集，含 IRC 与 NoEigenTest 相关内容
- [thyheal/gaussian-workflow-copilot](https://github.com/thyheal/gaussian-workflow-copilot)：中文 TS/IRC 自查清单
- 更多工具见 [../../Resources/github-repos.md](../../Resources/github-repos.md)

**本仓库脚本**：[`../scripts/scan_irc_profile.py`](../scripts/scan_irc_profile.py)（IRC/扫描能量曲线）、[`../scripts/log2xyz.py`](../scripts/log2xyz.py)（提取端点结构）、[`../scripts/g16_log_summary.py`](../scripts/g16_log_summary.py)（虚频与热化学汇总）、[`../scripts/opt_monitor.py`](../scripts/opt_monitor.py)（优化收敛曲线）、[`../scripts/g16_error_doctor.py`](../scripts/g16_error_doctor.py)（报错诊断）、[`../scripts/gjf_builder.py`](../scripts/gjf_builder.py)（生成 ts/qst2/irc 模板）。各脚本用 `-h` 查看用法。

**关键文献**

- K. Fukui, The path of chemical reactions – the IRC approach, *Acc. Chem. Res.* **1981**, 14, 363. DOI: [10.1021/ar00072a001](https://doi.org/10.1021/ar00072a001)
- C. Peng, H. B. Schlegel, Combining synchronous transit and quasi-Newton methods to find transition states, *Isr. J. Chem.* **1993**, 33, 449. DOI: [10.1002/ijch.199300051](https://doi.org/10.1002/ijch.199300051)
- C. Peng, P. Y. Ayala, H. B. Schlegel, M. J. Frisch, Using redundant internal coordinates to optimize equilibrium geometries and transition states, *J. Comput. Chem.* **1996**, 17, 49–56.
- C. Gonzalez, H. B. Schlegel, An improved algorithm for reaction path following, *J. Chem. Phys.* **1989**, 90, 2154. DOI: [10.1063/1.456010](https://doi.org/10.1063/1.456010)
- H. P. Hratchian, H. B. Schlegel, Accurate reaction paths using a Hessian based predictor–corrector integrator, *J. Chem. Phys.* **2004**, 120, 9918. DOI: [10.1063/1.1724823](https://doi.org/10.1063/1.1724823)

---

[← 上一章](06-freq-thermo.md) | [返回目录](../README.md) | [下一章 →](08-scan.md)
