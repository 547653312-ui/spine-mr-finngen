# REVIEW_01 — kimi-K3 审查与修订指挥（基于 Hy3 首轮交付）

审查对象：`analysis/` `data/` `results/`（Hy3 交付，2026-08-08）
对照标准：`plan/PLAN.md`（kimi-K3 自定方案）+ `FEASIBILITY.md`
审查角色：kimi-K3（规划 + 验证指挥）

---

## 0. 总体结论

**工程执行（A-）**：数据源替代决策（IEU→EBI GWAS Catalog + FinnGen R11）不仅合理，且结局的脊柱特异性严格优于原方案的 UKB `ieu-b-4975` 代理。自研 `mr_methods.py` 完整复现 IVW/Egger/加权中位/加权众数/LOO/Cochran Q，196 张敏感性图、FDR+Bonferroni 多重校正、诚实的 LIMITATIONS 均到位。**这是一版可直接进入修订、而非推倒重来的交付。**

**科学叙事（C+，必须重写）**：当前把 `WBC → SPONDINF`（OR=0.037, P=2e-4）当头条阳性，但①SPONDINF 仅 68–74 例，是全设计最小、最不可信的结局；②WBC 在 OM（OR=0.818, P=0.36）与 DISCITIS（OR=1.413, P=0.30）方向**完全相反且均不显著**，跨结局方向不一致 = 强烈的假阳性/结局特异性噪声信号。**不能把 SPONDINF 上的显著作为主结论。**

**发表可行性：可发（2–4 IF），前提是叙事改为"脊柱特异 OM 免疫因果 landscape 的系统评估——严谨阴性 + 探索性提示"，并补齐 PLAN 强制的 4 项分析。**

---

## 1. 与 PLAN.md 的偏差核对表

| PLAN 要求 | Hy3 实际 | 判定 |
|---|---|---|
| §2 数据 manifest 审计（防 ID 漂移，强制） | 未做 `manifest_*.tsv` | ❌ **缺口 P0** |
| §4 暴露 E1–E4（731 维免疫细胞） | 简化为 11 个 EFO 免疫/细胞因子 | ⚠ 偏离（IEU 不可用，合理降级） |
| §3 结局 ieub-4975 | 改用 FinnGen R11（脊柱特异更优） | ✅ 实际更好 |
| §5.2 LD clump | 距离剪枝 500kb（PLAN 允许降级③，且更保守） | ✅ 满足 |
| M3 敏感性 S1–S8 | Cochran Q/LOO/funnel 齐，缺 Steiger | ⚠ 部分（P2） |
| M3 MR-PRESSO | 用 Q+LOO+funnel 替代（无 R） | ⚠ 已声明 |
| §9 MVMR（吸烟/BMI/DM） | 未做 | ❌ **缺口 P1** |
| M4 drug-target / cis-pQTL | 未做 | ❌ **缺口 P1** |
| M5 脊柱复制 | 正向直接用 VOM/DC 作结局 | ✅ 已满足 |
| S9 三角验证/meta | 未做跨源 meta | ⚠ 部分（P2） |
| P6 power 计算 | LIMITATIONS 提及但未产出 power.tsv | ⚠ 部分（P1） |
| 多重校正 FDR+Bonf | 已做（45 次检验） | ✅ |
| 反向 MR | 因 EBI 429 限流，OM/DISCITIS 部分完成（5e-6 次级阈值） | ⚠ **P0 续跑** |

---

## 2. 必须修订项（指挥 Hy3）

### P0（阻塞成稿，必做）
- **P0-1 数据审计**：补 `data/manifest_exposures.tsv` / `manifest_outcomes.tsv`——对 11 个暴露 EFO 与 5 个 FinnGen 结局，用 GWAS Catalog API / FinnGen API 拉取 trait 名、样本量、ncase/ncontrol、人群、PMID、build，**人工核对 trait 名称与预期表型一致**，落盘。这是 PLAN 头号数据风险（ID 漂移）的防线。
- **P0-2 反向 MR 续跑**：重跑 `reverse_mr.py`（或 `wait_and_run_reverse.py`）；EBI 限流期过后自动续跑，完成 OM/DISCITIS 反向，写出 `MR_results_reverse.csv`。若反复限流无法完成，则**诚实降级**为"反向臂因脊柱感染 GWAS 显著位点极少（5e-8 下 0–2 个）+ 接口限流，仅报告可行性（reverse_MR_feasibility.csv / IV_discovery.csv），不进入主结论"，不阻塞。

### P1（显著提升严谨性，必做）
- **P1-1 MVMR**：新建 `analysis/05_mvmr.py`。暴露 = {候选免疫性状} + 吸烟`ukb-b-223` + BMI`ukb-b-15541` + 糖尿病`ukb-b-10753`；结局 = FinnGen OM / DISCITIS（芬兰人群，与 ukb-b 英国协变量**无样本重叠**，反而比 PLAN 担心的 ieub-4975 同源问题更安全）。报告条件 F 统计量 >10，写 `mvmr.tsv`。
- **P1-2 cis-pQTL / drug-target MR**：新建 `analysis/06_drug_target.py`。对 IL6R(`EFO_0008187`)、CD40、CRP 取 cis-pQTL（±1Mb）子集跑 OM/DISCITIS；对接 IL-6R 拮抗剂（tocilizumab）等临床可干预靶点的转化意义，写 `drug_target.tsv`。
- **P1-3 power 计算**：新建 `analysis/power.py`，对每个结局（尤其 VOM/SPONDINF）报告 80% 效能下可检出的最小 OR，写 `power.tsv`；论文中"阴性"须配可检出阈值，而非只报 P。

### P2（补强，建议做）
- **P2-1 Steiger 方向检验**：Python 实现（比较暴露/结局的 SNP 解释方差 R²，验证"暴露→结局"方向），对显著/探索性对跑，写 `steiger.tsv`。
- **P2-2 三角验证/meta**：在可行范围内做 OM ↔ DISCITIS ↔ VOM 固定/随机效应逆方差 meta（I²），写 `spine_replication.tsv`。

---

## 3. Manuscript 叙事重写方向（指挥 Hy3 起草）

**标题建议**（阴性 framing，仍含 spine）：
> *Genetic predisposition to immune cell and cytokine traits and the risk of osteomyelitis and spinal infection: a Mendelian randomization study*

**主结论（诚实版）**：
1. 在严格多重校正下，**免疫细胞/细胞因子性状与骨髓炎（含脊柱特异 VOM/DISCITIS）风险无稳健全基因组显著因果关联**。
2. SPONDINF（n=68–74）上 WBC/NEUT 的保护性信号**因样本量过小、且跨主要结局方向不一致，视为探索性、不可靠**，需更大样本验证——**不得作为头条**。
3. 方法学贡献：首个 spine-specific OM 免疫 MR、敏感性分析齐全、脊柱特异 framing、公开可复现管线。

**讨论必写**：① 阴性结果的解释边界（非"无因果"，而是"本样本量未检出"）；② VOM 仅 ~100 例的效能限制（配 power.tsv）；③ 距离剪枝替代 LD clump 的声明；④ 反向因果已被 MR 设计排除，但 SPONDINF 方向悖论需谨慎；⑤ 与 2025 OM↔心内膜炎 MR 的差异化（本研究的脊柱特异性）。

---

## 4. 最高优先级（一句话指挥 Hy3）

不夸大 SPONDINF；先补齐 **manifest 审计 + 反向续跑（P0）** 与 **MVMR + cis-pQTL/drug-target + power（P1）**，再按第 3 节叙事起草 manuscript。完成前不得声称"发现因果驱动"。
