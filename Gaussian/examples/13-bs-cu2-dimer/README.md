# 示例 13-3：双核 Cu(II) 的 Broken-Symmetry 计算与磁耦合常数 J

配套教程：[第13章 开壳层体系、波函数稳定性与 Broken-Symmetry](../../tutorials/13-openshell-stability.md)（§9）

## 目的

- 用 `Guess=(Fragment=N,Only)` 和负多重度为反铁磁耦合体系构造 broken-symmetry（BS）初猜；
- 按“高自旋（HS）→ 碎片初猜 → `Guess=Read` + `Stable=Opt`”的流程得到稳定的 BS 解；
- 用自旋密度和 ⟨S²⟩ 判断 BS 解是否正确；
- 用 Yamaguchi、Noodleman、Ruiz 公式计算 J，并写明哈密顿量约定。

## 文件说明

| 文件 | 说明 |
|---|---|
| `cu2_dimer_bs_guess_fragment.gjf` | 3 步作业，UB3LYP/def2-SVP，`NoSymm` |
| `calc_j_yamaguchi.py` | 从 HS、BS 两个计算（可在同一个 log 的不同步骤）读取能量与湮灭前 ⟨S²⟩，计算 J（仅用 Python 标准库） |
| `run_g16_slurm.sh` | Slurm 作业脚本（16 核、20 GB，对应 `%nprocshared=16`、`%mem=16GB`） |

体系为模型配离子 $[{\rm Cu_2(\mu\text{-}OH)_2(NH_3)_4}]^{2+}$（22 个原子，总电荷 +2）：平面 Cu2O2 核，Cu–O 1.93 Å，Cu–O–Cu 约 102°，Cu–N 2.02 Å，N–Cu–N 92°。坐标是按这些参数手工搭建的，**不是**晶体结构，也没有优化，仅用于演示流程。

| 步骤 | 电荷/多重度 | 内容 |
|---|---|---|
| 1 | `2 3` | 高自旋三重态 SCF + `Stable=Opt`，写入 `cu2_hs.chk` |
| 2 | `2 1 2 2 2 -2 -2 1` | `Guess=(Fragment=3,Only)`：碎片 1 = Cu1+2NH3（+2，α 二重态），碎片 2 = Cu2+2NH3（+2，β 二重态，多重度写 −2），碎片 3 = 两个 OH⁻（−2，单重态）；只生成初猜，写入 `cu2_bs.chk` |
| 3 | 从 chk 读取 | `Guess=Read Geom=AllCheck Stable=Opt`：BS 态 SCF 并保证稳定 |

## 运行步骤

```bash
cd Gaussian/examples/13-bs-cu2-dimer
g16 < cu2_dimer_bs_guess_fragment.gjf > cu2_dimer_bs_guess_fragment.log
# or
sbatch run_g16_slurm.sh cu2_dimer_bs_guess_fragment.gjf
```

计算完成后求 J（HS 为第 1 步，BS 为第 3 步，两个 S=1/2 中心所以 `--smax 1`）：

```bash
python calc_j_yamaguchi.py --hs cu2_dimer_bs_guess_fragment.log:1 \
                           --bs cu2_dimer_bs_guess_fragment.log:3 --smax 1
python calc_j_yamaguchi.py -h      # all options
```

HS 和 BS 也可以来自两个不同的 log，例如 `--hs hs.log --bs bs.log`（不写步骤号时取该文件最后一个含 SCF 能量的步骤）。

## 预期结果 / 如何检查

```bash
grep -c "Normal termination" cu2_dimer_bs_guess_fragment.log           # expect 3
grep "S\*\*2 before annihilation" cu2_dimer_bs_guess_fragment.log
grep -i "instability\|stable under\|already stable" cu2_dimer_bs_guess_fragment.log
grep -A 25 "Mulliken charges and spin densities:" cu2_dimer_bs_guess_fragment.log | grep " Cu "
```

| 检查项 | 正确结果 |
|---|---|
| 第 1 步（HS）⟨S²⟩ | ≈ 2.0 |
| 第 3 步（BS）⟨S²⟩ | ≈ 1.0 |
| Cu 自旋密度 | HS 两个 Cu 同号；BS 两个 Cu 大小相近、符号相反；绝对值小于 1（部分离域到 O、N） |
| 稳定性 | 第 1、3 步最后都是 `stable under the perturbations considered` 或 `already stable` |
| J 的符号 | Cu–O–Cu ≈ 102° 时预期为负（反铁磁），与羟桥双核铜的经验磁-构关联一致（Crawford et al., *Inorg. Chem.* 1976, 15, 2107） |

`calc_j_yamaguchi.py` 采用 $\hat H=-2J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$，同时打印 $\hat H=-J\,\hat{\mathbf S}_A\cdot\hat{\mathbf S}_B$ 约定下的数值（大 2 倍）。与文献比较前务必统一约定。

## 注意事项

- 若 BS 步骤的 ⟨S²⟩ ≈ 2 或 ≈ 0，说明塌缩到了 HS 或闭壳层：检查电荷多重度行中的负号、碎片划分，以及第 3 步是否读到了正确的 chk。
- J 对结构和泛函中 HF 交换比例很敏感。正式研究应使用实验结构或 HS 态优化结构，并在 Cu 上使用至少三 ζ 基组（如 def2-TZVP）做单点；示例用 def2-SVP 只为节省时间。
- 负多重度只能出现在 `Guess=Fragment` 的那一步；第 3 步用 `Geom=AllCheck` 读取，不需要再写碎片信息。
- 多核体系（3 个以上磁中心）需要多个 BS 构型并拟合多个 J，可参考官方 [Antiferromagnetic Coupling](https://gaussian.com/afc/) 说明和 [Guess 关键词](https://gaussian.com/guess/) 页面的 Fe2S2 示例。

返回教程：[第13章](../../tutorials/13-openshell-stability.md)
