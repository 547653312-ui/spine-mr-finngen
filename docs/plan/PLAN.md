# PLAN.md — 脊柱感染（骨髓炎/椎体骨髓炎）因果免疫学 MR 研究方案

> **版本** v1.0 · **规划者** kimi-K3（Planner + Reviewer）· **执行者** Hy3（Worker）
> **上游依据** `spine_MR/FEASIBILITY.md`
> **纪律** 本文件由 kimi-K3 维护；Hy3 只读本文件、只写 `analysis/ data/ results/`；kimi-K3 只写 `plan/ review/`。
> **本方案自包含**：所有阈值、公式、文件名、列名、判定门槛均已写死，执行者无需再提问。

---

## 0. 工作标题（Working Title）

**主标题（英文，投稿用）**

> **Causal immune and cytokine drivers of osteomyelitis and vertebral spinal infection: a two-sample bidirectional Mendelian randomization study**

**备选标题**

1. *Immune cell subsets, circulating cytokines and the risk of vertebral osteomyelitis and discitis: a bidirectional Mendelian randomization study*
2. *From immunity to the spine: genetically proxied IL-6 signalling and immune-cell traits in osteomyelitis — a drug-target Mendelian randomization analysis*

**中文题名**：脊柱感染（骨髓炎/椎体骨髓炎、椎间盘炎）的因果免疫驱动因素：一项双向孟德尔随机化研究

**Running title**：Immune causal drivers of spinal osteomyelitis (MR)

---

## 1. 背景与新颖性声明（Novelty Statement）

### 1.1 临床背景
脊柱感染（vertebral osteomyelitis, VOM；椎间盘炎, discitis, DC；硬膜外脓肿）发病率在老龄化与脊柱器械手术增加背景下持续上升，误诊率高、致残率高。宿主免疫状态被广泛认为是易感性与预后的关键决定因素，但既有证据几乎全部来自**观察性研究**（小样本、单中心、反向因果与残余混杂无法排除）：糖尿病、免疫抑制、慢性炎症与脊柱感染的关联究竟是**因**还是**果**，至今无因果层面的证据。

### 1.2 文献空白（Novelty）
- MR 领域中，**"脊柱特异性感染的因果免疫驱动"基本为空白**。2025 年可检索到的唯一相关 MR 工作仅报告了 **骨髓炎（OM）↔ 感染性心内膜炎（IE）** 的双向关系，未涉及免疫细胞亚群/细胞因子，更未触及 **VOM(M46.2) / DC(M46.4)** 这两个脊柱特异性表型。
- 因此本研究的三重新颖性：
  1. **首次**系统评估 731 项免疫细胞性状 + 循环细胞因子/免疫蛋白对骨髓炎的**因果**效应；
  2. **首次**将分析下沉到 **脊柱特异性亚型（VOM/DC）**，实现"从系统免疫到脊柱靶器官"的因果链条；
  3. 引入 **drug-target MR（IL6R cis-pQTL / cis-eQTL）**，把结论直接对接可干预靶点（如 IL-6R 拮抗剂），提升临床转化价值。

### 1.3 一句话卖点（供 Abstract 与 Cover Letter）
> *"我们提供了首个把系统性免疫组分与脊柱感染（含椎体骨髓炎、椎间盘炎）联系起来的因果证据，并将其锚定到可药物干预的 IL-6 通路。"*

---

## 2. 研究问题与假设

### 2.1 主研究问题（Forward / 正向）
> 遗传预测的**免疫细胞亚群频率**与**循环细胞因子/免疫蛋白水平**是否**因果性地**改变骨髓炎（及脊柱特异性椎体骨髓炎/椎间盘炎）的发病风险？

**H1（主假设）**：至少存在一组促炎轴性状（先验重点：IL-6 / IL-6R 信号、TNF-α、CRP、单核细胞与髓系亚群、调节性 T 细胞 Treg、CD8+ 效应记忆 T 细胞）对骨髓炎风险具有方向一致、可通过多方法与敏感性分析验证的因果效应。
**H0**：所有免疫性状的 IVW 因果估计与零效应无差异（OR = 1）。

### 2.2 反向研究问题（Reverse）
> 对骨髓炎的**遗传易感性**是否反过来改变免疫细胞组成/细胞因子水平（即观察性关联中"免疫异常"是否为疾病的**果**而非**因**）？

**H2**：若正向显著的性状在反向 MR 中不显著且 Steiger 定向检验支持"暴露→结局"，则支持真因果方向；若双向均显著，则报告为**双向/反馈环路**，不得单向解读。

### 2.3 脊柱特异性确认问题（Spine-specific confirmation）
> 在 OM 主分析中筛出的候选因果性状，其效应方向与量级能否在 **FinnGen VOM(M46.2) 与 DC(M46.4)** 中复现（哪怕仅达名义显著）？

**H3**：候选性状在 VOM/DC 中的效应方向与 OM 一致（同号），且合并（fixed/random-effects meta）后仍显著 → 支持"脊柱靶器官特异性"论述。

---

## 3. 数据来源（暴露 / 结局，含精确 ID）

> **强制要求（可审计）**：Hy3 必须先用 IEU `/gwasinfo/{id}` 或 GWAS Catalog API 拉取每个 ID 的元数据（trait 名称、样本量、ncase/ncontrol、population、year、PMID、build），写入 `data/manifest_exposures.tsv` 与 `data/manifest_outcomes.tsv`，**并人工核对 trait 名称与预期表型一致**。ID 漂移是本项目的头号数据风险，未落盘 manifest 的分析一律判为不合格。

### 3.1 结局（正向分析的 Outcome）

| 角色 | 表型 | 来源 | 病例/对照 | ID | 备注 |
|---|---|---|---|---|---|
| **主结局** | 骨髓炎 Osteomyelitis (M86) | UKB–IEU | 4,836 / 481,648 | **`ieu-b-4975`** | 含椎体病例，可自动下载，**主分析锚点** |
| 复制结局 A | 骨髓炎 OM | FinnGen（R9/R10 `M13_OSTEOMYELITIS`） | 2,336 / 473,264 | FinnGen 端点 | 独立样本复制，优先走 FinnGen 公开 summary（免 token 的 release 汇总文件） |
| **脊柱特异 B** | **椎体骨髓炎 VOM (M46.2)** | FinnGen | **111 / 353,224** | FinnGen 端点 `M13_VERTEBRALOSTEOMYELITIS`（以实际端点名为准） | ★核心 spine framing；需 token / 手动下载 |
| **脊柱特异 C** | **椎间盘炎 DC (M46.4)** | FinnGen | **557 / 353,224** | FinnGen 端点（M46.4 对应端点） | ★核心 spine framing |
| 阴性/对照结局 | 感染性心内膜炎 IE | UKB–IEU | 1,080 / 485,404 | `ieu-b-4972`（FEASIBILITY 记为 `ieu-b-49720`，**疑似笔误，必须先用 /gwasinfo 核验**） | 用于特异性对照：若某性状对 OM 与 IE 均有效，则属"泛感染易感"而非脊柱特异 |

### 3.2 暴露（正向分析的 Exposure）

按优先级分为 4 个暴露族（family），每族内独立做 FDR 校正。

| 族 | 内容 | 首选来源 | 备选 | n |
|---|---|---|---|---|
| **E1 免疫细胞性状** | **731 项免疫细胞表型**（绝对计数 AC、相对计数 RC、荧光强度 MFI、形态参数 MP） | **Orrù et al. 2020 (Nat Genet)，GWAS Catalog `GCST90001391`–`GCST90002121`**，Sardinia n=3,757 | FEASIBILITY 所列 `met-b`（Roederer 2015，150 项免疫子集，n≈669；样本量小、弱工具风险高，**降为敏感性来源**） | 731 |
| **E2 循环细胞因子** | 41 项细胞因子/生长因子（IL-6、TNF-α、IL-1β、MCP-1、IL-10 等） | Ahola-Olli et al. 2017 (AJHG)，n=8,293，IEU `ieu-a-` / GWAS Catalog `GCST004420`–`GCST004460` 区间（**逐一核验**） | Kalaoja 2021 / SCALLOP CVD-I | ~41 |
| **E3 免疫蛋白 (pQTL)** | Sun et al. 2018 SomaScan `prot-a-*`（n=3,301），重点 IL-6、IL-6R、sIL-6R、TNF-α、TNFRSF1A/1B、CRP、IL-1RA、IL-18、CXCL8 | `prot-a-*`（ID 需 /gwasinfo 检索确认） | Ferkingstad 2021 deCODE、Olink UKB-PPP | ~20–90（先验清单优先） |
| **E4 关键单性状（先验，强 power）** | CRP（`ukb-d-30710_irnt` 或 Said 2022，n>400k）、白细胞计数/中性粒/淋巴/单核（`ukb-d-30000_irnt` 系列或 Chen 2020 Blood Cell Consortium） | UKB 连续表型，样本量大 → F 值高 | — | ~8 |

**E4 说明**：E1/E2/E3 样本量偏小（3k–8k），弱工具偏倚风险实在；E4 用超大样本连续表型作为**高把握度的确证层**，是本方案对可行性的关键补强，必须执行。

### 3.3 多变量 MR 的协变量暴露（Confounder set）

| 变量 | ID | 用途 |
|---|---|---|
| 吸烟 | **`ukb-b-223`** | MVMR 校正 |
| 肥胖 / BMI | **`ukb-b-15541`** | MVMR 校正 |
| 糖尿病 | **`ukb-b-10753`** | MVMR 校正 |

> **⚠ 样本重叠告警**：`ukb-b-*` 协变量与结局 `ieu-b-4975` 同源 UKB → MVMR 存在样本重叠导致的偏向观察性估计的偏倚。**必须**在 Limitations 中明写，并在 `results/RESULTS.md` 中给出重叠比例的定性判断；有条件时用 FinnGen OM 作为 MVMR 结局做重叠规避的敏感性分析。

---

## 4. 研究设计：5 个分析模块

```
M1 正向主分析      E1–E4  →  ieu-b-4975 (OM)            [全自动，必做]
M2 反向 MR         ieu-b-4975 →  M1 显著性状             [全自动，必做]
M3 敏感性 + MVMR   Q / Egger / LOO / PRESSO / Steiger / MVMR(吸烟,BMI,DM)   [必做]
M4 drug-target MR  IL6R / TNF 等 cis-pQTL(±1Mb) → OM     [必做，高转化价值]
M5 脊柱特异性确认  候选性状 → FinnGen OM / VOM(M46.2) / DC(M46.4) + meta   [尽力做，token 依赖]
```

**分析流程门（gate）**：M1 产出候选清单（FDR q<0.05）→ 仅候选性状进入 M2/M3/M5，避免 731×多方法的算力爆炸。M4 独立于 M1 结果（先验靶点驱动）。

---

## 5. 工具变量（IV）选择

### 5.1 显著性阈值
- **主阈值**：全基因组显著 **P < 5×10⁻⁸**。
- **次阈值（secondary / 敏感性）**：**P < 1×10⁻⁵**。用于 E1/E2/E3 中主阈值下 SNP 数 < 3 的性状（免疫细胞性状极常见）。
- **规则**：主表以 5×10⁻⁸ 为准；若某性状主阈值 IV 数 <3，则以 1×10⁻⁵ 结果为主报告并在表中以 `pval_threshold` 列显式标注，**禁止混表不标注**。

### 5.2 连锁不平衡（LD）clumping
- **r² < 0.001，窗口 10,000 kb（10 Mb）**，参考面板 **1000 Genomes Phase 3 EUR**。
- 实现优先级：① IEU API `/ld/clump`；② 本地 PLINK 1.9 `--clump --clump-r2 0.001 --clump-kb 10000`（1000G EUR bfile）；③ **降级方案**（仅当 ①②均不可用）：按 P 值升序贪心保留，剔除与已保留 SNP 距离 <10 Mb 且同染色体者（纯距离剪枝），**必须在 `results/RESULTS.md` 与论文 Methods 中如实声明为距离剪枝**。
- **例外**：M4 drug-target cis 分析中，允许放宽至 **r² < 0.1**（并使用考虑 LD 相关矩阵的 IVW，`correl=TRUE` 等价实现），以提高 cis 区域信息利用率；此时必须提供 LD 矩阵来源。

### 5.3 工具强度
- 逐 SNP **F = β_X² / SE_X²**；报告 `F_mean`、`F_min`、`F<10 的 SNP 数`。
- 或用解释方差法：`R² = 2·EAF·(1−EAF)·β²`（连续、标准化暴露），`F = R²(N−k−1) / (k(1−R²))`。
- **剔除规则**：单 SNP **F < 10** 者从主分析剔除；若剔除后 IV 数 <3，则该性状标记 `weak_instrument=TRUE` 并降为探索性结果，不进入主结论。
- MVMR 报告 **条件 F 统计量（Sanderson–Windmeijer）**，阈值 10。

### 5.4 其他 IV 质控（按序执行，逐步记录剔除数）
1. 仅保留 **常染色体**、双等位、有 rsID 的 SNP；剔除 MHC 区（chr6:25–35 Mb）**在敏感性分析中重复一次**（主分析保留，敏感性剔除，因免疫性状高度依赖 MHC — 这是本项目必须做的一步）。
2. 剔除与**结局**在 P < 5×10⁻⁸ 水平直接关联的 SNP（避免明显水平多效性）。
3. 用 PhenoScanner / IEU `phewas` 检索每个 IV 的既往关联，**手动剔除**直接关联吸烟/BMI/糖尿病/免疫抑制用药者，并在附表列出被剔除 SNP 与原因。
4. **Steiger filtering**：剔除 `R²_outcome > R²_exposure` 的 SNP（Hemani 2017），并单独报告 filtered 前后结果。

---

## 6. 数据协调（Harmonization）—— 审计重点

按 rsID 合并暴露与结局数据，规则如下（Hy3 必须逐条实现，并把每步剔除计数写入 `data/harmonized/_harmonization_log.tsv`）：

1. **对齐效应等位基因**：以暴露的 EA 为基准；若结局的 EA/OA 与暴露互换，则 `beta_outcome = −beta_outcome`，`eaf_outcome = 1 − eaf_outcome`。
2. **链翻转**：若等位基因经互补（A↔T, C↔G 映射）后可匹配，则翻转结局等位基因后再对齐。
3. **回文 SNP（A/T、C/G）**：
   - MAF > **0.42** → **直接剔除**（无法定链）；
   - MAF ≤ 0.42 → 用 EAF 推断链方向（暴露 EAF 与结局 EAF 同侧则保留；异侧则翻转），等价于 TwoSampleMR `action=2`。
   - 若任一侧缺失 EAF → **剔除**该回文 SNP。
4. **等位基因不相容**（如暴露 A/G vs 结局 C/T 且互补后仍不匹配）→ 剔除。
5. **结局缺失 SNP**：优先用 **LD proxy（r² > 0.8，1000G EUR，LDlink/LDproxy 或本地面板）**替代；无 proxy 则剔除。proxy 使用需记录 `proxy_rsid, r2`。
6. **重复 rsID**：保留暴露 P 值最小者。
7. 最终协调数据表落盘为 `data/harmonized/{exposure_id}__{outcome_id}.tsv`，**必须包含**列：
   `SNP, chr, pos, effect_allele, other_allele, eaf_exposure, beta_exposure, se_exposure, pval_exposure, eaf_outcome, beta_outcome, se_outcome, pval_outcome, samplesize_exposure, samplesize_outcome, F, palindromic, ambiguous, proxy_rsid, r2_proxy, steiger_keep, mr_keep`

> **审计口径**：kimi-K3 将随机抽 ≥5 个 SNP 手工复核符号与等位基因对齐；任何一处符号错误即判 **BLOCKER**（因为它会系统性反转 OR 方向）。

---

## 7. MR 方法（主分析 + 多方法一致性）

> 若 PyPI 的 `MendelianRandomization` 包可安装则优先调用；**若安装失败或功能缺失，按下列公式在 `numpy/scipy` 中原生实现**（本项目以下列公式为唯一权威定义，便于复核）。所有实现须在 `analysis/mr_core.py` 中，并附单元测试 `analysis/test_mr_core.py`（用已发表算例或模拟数据验证 IVW 与 Egger 的还原精度，误差 <1e-6）。

设第 j 个 SNP：暴露效应 `X_j ± σ_Xj`，结局效应 `Y_j ± σ_Yj`，J 个 SNP。

### 7.1 IVW（**主方法**，multiplicative random-effects）
```
w_j    = 1 / σ_Yj²
β_IVW  = Σ (w_j · X_j · Y_j) / Σ (w_j · X_j²)
SE_fix = sqrt( 1 / Σ (w_j · X_j²) )
Q      = Σ w_j (Y_j − β_IVW·X_j)²          , df = J − 1
SE_re  = SE_fix · max(1, sqrt( Q / (J−1) ))    ← 主报告用 random-effects SE
```
- J = 1 → 退化为 **Wald ratio**：`β = Y/X`，`SE = σ_Y/|X|`（一阶 delta）。
- J = 2 → 仅 IVW（固定效应），不做 Egger/median/mode。

### 7.2 MR-Egger
以 `1/σ_Yj²` 为权重，对 `Y_j ~ α + β·X_j` 做加权线性回归（含截距）。
- `α`（截距）= 定向水平多效性检验统计量；报告 `α, SE(α), P(α)`。
- Egger SE 需乘以残差标准误缩放（等价于 RE），并报告 **I²_GX**（NOME 假设）；**I²_GX < 0.9 时**必须做 **SIMEX 校正**，否则 Egger 斜率有回归稀释偏倚。

### 7.3 加权中位数（Weighted Median）
```
β_j    = Y_j / X_j
Var_j  = σ_Yj² / X_j²      (一阶 delta)
w_j    = 1 / Var_j          → 归一化后取加权中位数
SE     : 参数化 bootstrap，B = 10,000，seed = 20260808
```

### 7.4 加权众数（Weighted Mode）
- 对 `β_j` 用高斯核加权密度估计，权重同上；带宽采用修正 Silverman 规则，**φ = 1**（同时报告 φ=0.5 的敏感性）。
- SE：参数化 bootstrap，B = 10,000，同 seed。

### 7.5 结果换算
二分类结局：`OR = exp(β)`，`95%CI = exp(β ± 1.96·SE)`。
- **效应单位必须写清**：连续暴露为"每 1 SD 增加"，二分类暴露（反向 MR 的 OM）为"每 log-odds 增加"。

### 7.6 主结论判定规则（**必须写进论文 Methods**）
某暴露–结局对被判定为"**稳健因果证据**"须**同时满足**：
1. IVW **FDR q < 0.05**（族内 Benjamini–Hochberg）；
2. MR-Egger / 加权中位数 / 加权众数 中 **≥ 2 种**方法效应**同号**（不要求各自显著）；
3. **Egger 截距 P > 0.05**；
4. Cochran **Q 的 P > 0.05**，或虽异质但 RE-IVW 仍显著；
5. **留一法**无单 SNP 驱动（去掉任一 SNP 后 IVW P 仍 < 0.05）；
6. **Steiger 定向**支持 exposure→outcome；
7. `F_mean > 10`。
满足 1 且不满足 2–7 中任一项 → 归入"**提示性证据（suggestive）**"，措辞降级，不得进入 Abstract 结论句。

---

## 8. 敏感性分析清单（全部必做）

| # | 分析 | 判定/报告 |
|---|---|---|
| S1 | **Cochran Q**（IVW 与 Egger 各一份） | 报告 Q、df、P、**I² = max(0,(Q−df)/Q)** |
| S2 | **MR-Egger 截距** | 截距、SE、P；P<0.05 判定存在定向多效性 |
| S3 | **留一法 (Leave-one-out)** | 逐 SNP 剔除后 IVW；输出表 + 森林图；标注"驱动 SNP" |
| S4 | **MR-PRESSO** | Python 实现：global test（RSS_obs vs 10,000 次模拟 RSS 分布）→ outlier test（逐 SNP P，Bonferroni）→ distortion test（剔除前后 β 差异 P）。global P<0.05 时报告剔除离群点后的 corrected 估计 |
| S5 | **单 SNP 分析 + 漏斗图** | 检查不对称 |
| S6 | **Steiger 定向检验** | 逐 SNP + 整体；报告 `steiger_pval`、`correct_causal_direction` |
| S7 | **MHC 区敏感性** | 剔除 chr6:25–35 Mb 后重跑主分析（免疫性状必做） |
| S8 | **阈值敏感性** | 5e-8 与 1e-5 两套结果并列比较 |
| S9 | **复制/三角验证** | ieu-b-4975 ↔ FinnGen OM ↔ VOM ↔ DC；固定/随机效应 meta（逆方差合并 β），报告 I² |
| S10 | **特异性对照** | 同一性状 → IE (`ieu-b-4972`)；若同样显著，弱化"脊柱特异"表述 |
| S11 | **MVMR** | 见 §9 |
| S12 | **把握度计算** | 见 §10.2 |

---

## 9. 多变量 MR（MVMR）

- **模型**：结局 = OM（`ieu-b-4975`；重叠敏感性用 FinnGen OM）；暴露 = {候选免疫性状} + {吸烟 `ukb-b-223`} + {BMI/肥胖 `ukb-b-15541`} + {糖尿病 `ukb-b-10753`}。
- **IV 池**：合并各暴露的显著 SNP → 去重 → 统一 clump（r²<0.001, 10 Mb）→ 在**所有**暴露与结局中提取效应并协调。
- **估计**：无截距的多元加权最小二乘，权重 `1/σ_Yj²`：
  `β̂ = (Xᵀ W X)⁻¹ Xᵀ W Y`，`Var(β̂) = (Xᵀ W X)⁻¹`（并做 RE 缩放）。
- **诊断**：条件 F 统计量（每个暴露）> 10；MVMR-Egger 截距；Q_A 异质性。
- **解读**：若免疫性状效应在校正后大幅衰减 → 该效应部分由代谢/行为混杂中介，须在 Discussion 明确。

---

## 10. 多重检验与把握度

### 10.1 多重检验
- **族内 FDR（Benjamini–Hochberg, q<0.05）**，四族独立校正：E1(731) / E2(41) / E3(~90) / E4(~8)。
- 同时报告 **Bonferroni** 阈值（如 E1: 0.05/731 = 6.8×10⁻⁵）作为最保守参照。
- 主表须同时给 `pval` 与 `fdr_q` 两列。

### 10.2 把握度（Power）
二分类结局 MR 近似（Burgess 2014）：
```
NCP   = N · k(1−k) · R²_X · β²        (β = ln OR)
Power = Φ( sqrt(NCP) − Φ⁻¹(1 − α/2) ) ,  α = 0.05
```
- 对 `ieu-b-4975`（N=486,484，k=4836/486484≈0.99%）与 VOM（N=353,335，k≈0.031%）分别绘制"**R² × OR → power**"热力图，输出 `results/figures/power_curves.png`。
- **预判**：VOM(n=111) 仅对 OR ≥ ~2.0 且 R²≥2% 的暴露有可接受把握度 → **VOM 阴性结果一律解释为"把握度不足"，禁止解释为"无效应"**。

---

## 11. 软件与环境（R 不可用 → 纯 Python）

```
python >= 3.10
pandas, numpy, scipy, statsmodels, requests, tqdm, matplotlib, seaborn, pyarrow
MendelianRandomization (PyPI, 若可用；否则原生实现，见 §7)
可选：scikit-learn(仅用于 bootstrap 并行), pysam/cyvcf2(读 IEU VCF), plink1.9(若本机可用)
```
- **依赖锁定**：`analysis/requirements.txt` + `pip freeze > analysis/env_freeze.txt`。
- **随机种子统一 `SEED = 20260808`**（bootstrap、PRESSO 模拟均须设种）。
- **数据获取**：
  - IEU OpenGWAS REST：`https://gwas-api.mrcieu.ac.uk/`，端点 `/gwasinfo/{id}`、`/tophits`、`/associations`、`/ld/clump`。
    **⚠ 现已强制 JWT 鉴权**：需在 https://api.opengwas.io 注册取 token，置于环境变量 `OPENGWAS_JWT`，请求头 `Authorization: Bearer $OPENGWAS_JWT`。**禁止把 token 明文写入代码或提交**（用 `os.environ.get`，缺失时给出清晰报错）。
  - **降级路径 D1**：直接下载 GWAS-VCF `https://gwas.mrcieu.ac.uk/files/{id}/{id}.vcf.gz`（+ `.tbi`），本地解析。
  - **降级路径 D2**：GWAS Catalog summary statistics FTP（Orrù 2020 的 GCST9000xxxx 系列走这条路最稳）。
  - **FinnGen**：优先公开 release 的 `finngen_R{n}_{ENDPOINT}.gz`（多数无需 token）；若受限则记录为待办并走 M1–M4，不阻塞主线。
- **限流与缓存**：所有 API 请求 **≤10 req/s**，失败指数退避重试 3 次；**所有下载结果落盘缓存**到 `data/raw/`，重跑不得重复拉取。全流程日志 `analysis/logs/run_YYYYMMDD_HHMM.log`。

---

## 12. 交付物规范（Hy3 必须严格按此命名，便于审计）

```
spine_MR/
├─ analysis/
│   ├─ 00_fetch_metadata.py      # 拉 gwasinfo → manifest
│   ├─ 01_extract_instruments.py # tophits + clump + F
│   ├─ 02_harmonize.py           # §6 规则
│   ├─ 03_mr_run.py              # M1/M2 批量
│   ├─ 04_sensitivity.py         # S1–S8
│   ├─ 05_mvmr.py                # M3-MVMR
│   ├─ 06_drug_target.py         # M4 cis-pQTL
│   ├─ 07_spine_replication.py   # M5 FinnGen + meta
│   ├─ 08_figures.py             # 所有图
│   ├─ mr_core.py + test_mr_core.py
│   ├─ requirements.txt, env_freeze.txt, logs/
├─ data/  raw/ instruments/ harmonized/ manifest_*.tsv
├─ results/
│   ├─ mr_main.tsv               # 见下列名规范（长表）
│   ├─ mr_reverse.tsv, sensitivity_loo.tsv, sensitivity_presso.tsv,
│   ├─ mvmr.tsv, drug_target.tsv, spine_replication.tsv, power.tsv
│   ├─ figures/  (forest_main, scatter_top, funnel_top, loo_top,
│   │             volcano_731, power_curves, heatmap_spine)  *.png(300dpi)+*.pdf
│   └─ RESULTS.md                # 结果叙述 + 所有判定门槛逐条自查
└─ review/  (kimi-K3 专用)
```

**`results/mr_main.tsv` 必备列**（缺列即判不合格）：
`exposure_id, exposure_name, exposure_family, outcome_id, outcome_name, direction, pval_threshold, method, nsnp, b, se, pval, or, or_lci95, or_uci95, fdr_q, F_mean, F_min, n_weak, q_stat, q_df, q_pval, i2, egger_intercept, egger_intercept_se, egger_intercept_pval, presso_global_p, n_outliers, steiger_dir, steiger_pval, loo_max_p, notes`

**图件清单**：
1. 火山图（731 免疫性状 β vs −log10 P，标注 FDR 显著者）
2. 主森林图（候选性状 4 方法 OR/95%CI）
3. 散点图 + 漏斗图 + 留一法森林图（Top 3–5 性状）
4. 脊柱三联热图（OM / VOM / DC 效应量对照）
5. Power 曲线
6. 研究设计示意图（三大 MR 假设 IV1–IV3 + 数据流）

---

## 13. 时间线（分阶段，含审查门）

| 阶段 | 内容 | 产出 | **kimi-K3 审查门** |
|---|---|---|---|
| **P0 环境** | 装依赖、跑通 `mr_core` 单元测试、验证 OPENGWAS_JWT | `env_freeze.txt`、测试通过日志 | G0：单元测试必须全绿，否则不得进入 P1 |
| **P1 数据清点** | 解析所有 ID 元数据；核验 `ieu-b-4972` 真身；生成 manifest | `manifest_*.tsv` | G1：ID 与表型逐一对得上；样本量与 FEASIBILITY 一致 |
| **P2 IV 与协调** | E1–E4 取 tophits、clump、F、协调 | `instruments/`、`harmonized/`、协调日志 | G2：抽查 5 SNP 手工复核对齐；回文/剔除计数可追溯 |
| **P3 正向主分析** | M1 全量 → 火山图 + 候选清单 | `mr_main.tsv` | G3：FDR 实现正确；候选清单合理（非全阳/全阴） |
| **P4 敏感性 + MVMR** | S1–S8、S11 | `sensitivity_*.tsv`、`mvmr.tsv` | G4：多效性红旗排查；条件 F |
| **P5 反向 MR** | M2 + Steiger | `mr_reverse.tsv` | G5：方向性结论自洽 |
| **P6 脊柱确认 + 药靶** | M5（FinnGen VOM/DC）、M4（cis-pQTL） | `spine_replication.tsv`、`drug_target.tsv` | G6：power 说明到位；阴性不过度解读 |
| **P7 出图出表 + 撰写** | 图 1–6、Table 1–4、`RESULTS.md` | 完整结果包 | G7：终审（§14 全量 rubric） |

> 阶段间**必须停下等待 kimi-K3 审查意见**（`review/REVIEW_round{n}.md`）后再推进。P6 若因 token 阻塞，允许并行推进 P7 的 M1–M4 部分，但论文标题中的 spine framing 须相应调整为"osteomyelitis (including vertebral involvement)"。

---

## 14. 风险登记册（Risk Register）

| # | 风险 | 概率 | 影响 | 缓解措施 | 触发信号 |
|---|---|---|---|---|---|
| R1 | **FinnGen VOM/DC 需 token / 不可自动下载** | 高 | 高（削弱脊柱特异性卖点） | ①先跑 M1–M4 不阻塞；②尝试 FinnGen 公开 release 文件与 risteys 端点；③退路：以 OM（含椎体）为主体 + 文献证据 + IE 阴性对照支撑脊柱论述；④标题降级为 "osteomyelitis including vertebral infection" | P6 下载 403/需登录 |
| R2 | **VOM 病例仅 111 → 极低把握度** | 确定 | 中 | 预先出 power 曲线；**阴性一律表述为 underpowered**；VOM/DC/OM 三者 meta 提升效能；VOM 结果定位为"方向一致性确认"而非独立检验 | — |
| R3 | **OpenGWAS API 鉴权/限流变更** | 中高 | 高（数据取不到） | JWT 环境变量 + 全量本地缓存 + GWAS-VCF/GWAS Catalog FTP 双降级路径 | 401/429 |
| R4 | **弱工具偏倚**（E1/E2 样本 3k–8k，1e-5 阈值） | 高 | 中 | F 值全量报告、F<10 剔除、E4 大样本确证层、两阈值并列、双样本设计下偏倚趋零 | `F_mean < 10` 或 `n_weak/nsnp > 0.2` |
| R5 | **水平多效性**（免疫性状高度相关、MHC 区） | 高 | 高 | Egger 截距、PRESSO、加权中位数/众数、MHC 剔除敏感性、PhenoScanner 手工筛 SNP、MVMR | 截距 P<0.05 或 PRESSO global P<0.05 |
| R6 | **样本重叠**（UKB 协变量 × UKB 结局） | 高 | 中 | 明写 Limitation；FinnGen OM 作为无重叠 MVMR 敏感性；强工具（F 大）下重叠偏倚有限 | — |
| R7 | **仅欧洲人群 → 外推受限** | 确定 | 中 | Limitation 明写；不宣称跨族群普适；建议后续东亚/多族群验证 | — |
| R8 | **ID 漂移/错标**（如 `ieu-b-49720`） | 中 | 高（全盘错） | P1 强制 manifest 核验 + 人工确认 trait 名 | 名称与预期不符 |
| R9 | **731×多方法算力/时间超标** | 中 | 低 | 门控设计（仅候选进 M2/M3/M5）；bootstrap 仅对候选做；缓存复用 | 单阶段 >8h |
| R10 | **"泛感染易感"而非脊柱特异** | 中 | 中 | IE 阴性对照（S10）；若泛化，转为"感染易感的免疫基础，脊柱为受累靶器官之一"叙事，仍可发表 | IE 与 OM 同向同显著 |
| R11 | **可重复性缺失**（无种子/无版本） | 中 | 高 | SEED 固定、env_freeze、日志、所有中间文件落盘 | — |
| R12 | **过度解读 / 因果措辞过强** | 中 | 中（审稿被拒） | §7.6 判定规则强制分级；Abstract 仅写"稳健"级结论 | 终审发现 |

---

## 15. Validation & Revision（kimi-K3 审计标准与修订指挥流程）

### 15.1 审计评分卡（Audit Rubric）
每轮审查按下表逐项打分：**PASS / MINOR / MAJOR / BLOCKER**。存在任一 **BLOCKER** → 打回重跑，不得进入下一阶段。

| 域 | # | 检查项 | 判定标准 | 严重级别 |
|---|---|---|---|---|
| **A 数据完整性** | A1 | 每个 GWAS ID 有 manifest 记录且 trait 名与预期一致 | 100% 覆盖 | BLOCKER |
| | A2 | `ieu-b-4972/49720` 真身已核验 | 有 gwasinfo 截图/JSON 落盘 | MAJOR |
| | A3 | 样本量/病例数与 FEASIBILITY 一致（或已说明差异） | 差异 <5% 或有说明 | MAJOR |
| **B 协调正确性** | B1 | 随机抽 5 个 SNP 手工复核 EA/OA 与 β 符号 | 0 处错误 | **BLOCKER** |
| | B2 | 回文 SNP 处理符合 §6.3（MAF>0.42 剔除） | 代码可见且日志计数吻合 | BLOCKER |
| | B3 | 协调日志含逐步剔除计数（起始→最终 nsnp 可追溯） | 完整链条 | MAJOR |
| | B4 | proxy SNP 有 r² 记录 | 有则记录，无则声明未用 proxy | MINOR |
| **C 工具强度** | C1 | `F_mean/F_min/n_weak` 三列齐全 | 齐全 | MAJOR |
| | C2 | F<10 的 SNP 已按规则处理 | 一致 | MAJOR |
| | C3 | MVMR 报告条件 F | 有 | MAJOR |
| **D 方法正确性** | D1 | `test_mr_core.py` 全绿（IVW/Egger 数值还原） | 通过 | BLOCKER |
| | D2 | IVW 用 random-effects SE（Q>df 时） | 一致 | MAJOR |
| | D3 | 加权中位数/众数用 bootstrap SE 且设种 | 一致 | MAJOR |
| | D4 | J=1/J=2 的退化处理正确（不报 Egger） | 无非法输出 | MAJOR |
| | D5 | FDR 在族内、方法内正确应用（不跨方法混算） | 一致 | MAJOR |
| **E 方向性** | E1 | 反向 MR 已对所有候选执行 | 100% | MAJOR |
| | E2 | Steiger 结果与结论叙述一致 | 一致 | MAJOR |
| | E3 | 双向显著者未被单向解读 | 表述正确 | MAJOR |
| **F 多效性红旗** | F1 | Egger 截距 P<0.05 的结果未进入主结论 | 无违规 | **BLOCKER** |
| | F2 | Q 显著者已用 RE-IVW 并说明 | 一致 | MAJOR |
| | F3 | PRESSO global P<0.05 者报告了 corrected 估计 | 有 | MAJOR |
| | F4 | 留一法驱动 SNP 已识别并披露 | 有 | MAJOR |
| | F5 | MHC 敏感性已跑 | 有 | MAJOR |
| **G 可重复性** | G1 | 随机抽 2 条主结果，由 kimi-K3 用协调表独立复算 IVW β/SE/OR | 相对误差 <1% | **BLOCKER** |
| | G2 | SEED 固定、env_freeze 存在、日志完整 | 齐全 | MAJOR |
| | G3 | 从 `data/harmonized/` 可完全重算 `results/` | 可 | MAJOR |
| | G4 | 表中 OR 与 β 自洽（OR = exp(β)，CI 对称于 log 尺度） | 全部自洽 | MAJOR |
| **H 脊柱叙事** | H1 | 标题/摘要/讨论均出现 spine/vertebral 定位且与实际数据匹配（未越界宣称） | 匹配 | MAJOR |
| | H2 | VOM/DC 阴性表述为 underpowered 而非 "no effect" | 正确 | MAJOR |
| | H3 | IE 特异性对照已执行并被讨论 | 有 | MAJOR |
| | H4 | Limitation 覆盖：欧洲人群、样本重叠、VOM 小样本、UKB 自报/ICD 定义误分类 | 4 项齐全 | MAJOR |
| **I 呈现** | I1 | 图 6 件齐、300 dpi、坐标轴与单位标注完整 | 齐全 | MINOR |
| | I2 | STROBE-MR 清单逐条对应 | 有 `review/STROBE_MR_checklist.md` | MAJOR |

### 15.2 kimi-K3 的独立复算规程（防"看起来对"）
1. 从 `data/harmonized/` 随机取 2 张协调表，用 §7.1 公式在独立脚本中手算 `β_IVW, SE_re, OR, 95%CI`，与 `results/mr_main.tsv` 比对（容差 1%）。
2. 抽 5 个 SNP 回到 `data/raw/` 原始文件核对 `beta/se/EA/OA`，验证 §6 对齐链条。
3. 抽 1 个性状核算 FDR 排序位次。
4. 核对至少 1 张图与其数据源表数值一致。

### 15.3 修订指挥协议
- 审查意见写入 **`review/REVIEW_round{n}.md`**，格式：
  ```
  [BLOCKER|MAJOR|MINOR] <检查项编号> <文件:行/表:列>
  现象： ...
  依据： PLAN.md §x.y
  要求： <可执行的具体修改指令>
  验收： <我将如何复验>
  ```
- Hy3 修订后在同文件追加 `RESPONSE:` 段（逐条对应，附改动文件与新数值），**禁止"已修复"式空回应**。
- 最多 3 轮；第 3 轮后仍有 BLOCKER → 上报编排者，考虑收缩范围（如砍掉 M5 或降级为 preprint/更低 IF 期刊）。
- kimi-K3 **不直接修改** Hy3 的任何产物文件。

---

## 16. 目标期刊（IF 2–4）与投稿策略

| 优先 | 期刊 | 大致 IF | 契合点 | 备注 |
|---|---|---|---|---|
| 1 | **Frontiers in Immunology** | ~5.7（Sec. 可控） | 免疫细胞性状 + 感染，MR 稿件接受度高、周期短 | 首投；Section: Inflammation / Microbial Immunology |
| 2 | **Journal of Translational Medicine** | ~6 | drug-target MR（IL-6R）转化叙事 | 若 M4 结果强则优先 |
| 3 | **Scientific Reports** | ~4 | 广谱、稳妥兜底 | 保底 |
| 4 | **Infectious Diseases and Therapy** | ~4.5 | 感染病临床读者 | 临床向包装 |
| 5 | **BMC Genomics / BMC Musculoskeletal Disorders** | ~3 | 遗传/骨科读者 | 备选 |

- 全文遵循 **STROBE-MR** 报告规范（清单落到 `review/STROBE_MR_checklist.md`）。
- 伦理：均为公开汇总数据，原研究已获伦理批准与知情同意，本研究无需额外审批（须在文中声明）。
- 数据可得性声明：列出所有 GWAS ID 与下载地址；代码可公开（GitHub/Zenodo DOI）。

---

## 17. 给 Hy3 的执行摘要（TL;DR 版行动清单）

1. 建环境 → 写 `mr_core.py` + 单测（IVW/Egger/median/mode/Q/PRESSO/Steiger）→ **G0**。
2. 拉元数据建 manifest，核验 `ieu-b-4975`、`ukb-b-223/15541/10753`、`ieu-b-4972`、Orrù 731、Ahola-Olli 41、prot-a 先验清单 → **G1**。
3. 取 IV（5e-8 主 / 1e-5 次）→ clump(r²<0.001, 10 Mb) → 算 F → 协调（§6 七条规则，出日志）→ **G2**。
4. 跑 M1 正向 → 火山图 + FDR → 候选清单 → **G3**。
5. 候选跑 S1–S8 + MVMR → **G4**；跑 M2 反向 + Steiger → **G5**。
6. 跑 M4 cis-pQTL；尝试 M5 FinnGen VOM/DC（阻塞则记录并继续）→ **G6**。
7. 出 6 图 + 4 表 + `RESULTS.md`（逐条对照 §7.6 判定门槛自查）→ **G7 终审**。

**永远不要**：跳过协调日志、混表不标 P 阈值、把 Egger 截距显著的结果写进 Abstract、把 VOM 阴性说成"无因果效应"、把 token 写进代码。
