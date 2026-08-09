# 局限性说明（LIMITATIONS）

本文件如实记录本次 MR 分析的设计取舍、数据限制与解释边界，供 kimi-K3 复核与后续修订。

---

## 一、关于 FinnGen VOM(M46.2) / DC(M46.4) 的 token 问题 —— 已解决，但降了一个版本

FEASIBILITY.md 记载 VOM/DC「需 FinnGen token」，因此原计划用 UKB 的 `ieu-b-4975`
骨髓炎作代理，并预期存在"脊柱框架"缺口。**本次实际情况优于预期：**

- FinnGen **R12** 汇总统计确实非公开（bucket 返回 404，需 DUA 授权）。
- 但 FinnGen **R11** 的同名表型在 Google 公共桶中**完全开放**，且附带 tabix 索引，
  可无 token 自动获取。病例数与 FEASIBILITY.md 记载一致（R11 略低于 R12 是正常的版本差异）。

| 表型 | phenocode | R11 病例（本次使用） | FEASIBILITY 记载（R12） |
|---|---|---|---|
| 椎体骨髓炎 M46.2 | `M13_OSTEOMYELVERTEB` | 104 | 111 |
| 椎间盘炎 M46.4 | `M13_DISCITIS` | 495 | 557 |
| 化脓性椎间盘感染 | `M13_DISCINFECTION` | 375 | — |
| 感染性脊柱病 | `M13_SPONDYLOINFECTION` | 68 | — |
| 骨髓炎（全部位 M86） | `M13_OSTEOMYELITIS` | 2125 | 2336 |

**因此 `ieu-b-4975` 未被使用**，原因有二：
1. OpenGWAS 自 2024-05-01 起强制 JWT，无 token 时该 ID 无法自动下载（见 `API_NOTES.md`）；
2. 更重要的是，FinnGen R11 直接提供**脊柱特异**结局，科学上严格优于 UKB 的混合骨髓炎表型
   —— 后者不区分长骨/椎体，用作"脊柱感染"代理本身就是稀释效应的来源。

**遗留限制**：本次用的是 R11 而非 R12，样本量比 FEASIBILITY 表格略小（VOM 104 vs 111），
统计效能相应略低。若后续获得 FinnGen DUA，把 `gwas_io.FinnGenSumstats(release="R12")`
指向 R12 文件即可复跑，无需改动任何分析逻辑。

## 二、统计效能：这是本研究最主要的限制

椎体骨髓炎仅 **104 例**。在这种极端不平衡设计下：

- 单个 IV 的结局效应 SE 普遍在 0.14–0.20（log-OR 尺度），意味着即使真实 OR = 1.5，
  也需要远多于现有 IV 数才能达到 80% 效能；
- 因此 **VOM/SPONDYLOINFECTION 上的阴性结果不能解释为"无因果效应"**，
  只能解释为"本样本量下未能检出"；
- 相对地，`M13_OSTEOMYELITIS`（2125 例）与 `M13_DISCITIS`（495 例）的结果更可采信。

建议在论文中对每个暴露-结局对报告可检出的最小 OR（power calculation），而非只报 P 值。

## 三、LD 剪枝用距离近似替代，未使用 LD 参考面板

标准做法是用 1000G EUR 面板做 r² < 0.001 的 clumping。本环境中：
- R 不可用 → 无 `TwoSampleMR::clump_data`；
- 无 plink 二进制，也无 1000G 基因型文件；
- OpenGWAS 的 `/ld/clump` 端点需 token（同样 401）。

因此 `mr_pipeline.distance_prune()` 采用**距离剪枝**：按 P 升序贪心保留，
剔除同染色体 ±500 kb 内的其他位点。这是无面板时的通行替代方案，但：
- 可能**残留弱连锁**（>500 kb 但仍有 LD 的位点罕见，风险较低）；
- 也可能**过度剔除**同一区域内真正独立的次级信号（如 CRP 位点附近的条件独立信号），
  导致 IV 数偏少、效能进一步损失。

补救建议：取得 plink + 1000G 后，用真正的 clumping 重跑 IV 选择；接口已隔离在单个函数中。

## 四、暴露端 SE 由 beta 与 P 反推

EBI 汇总统计 API 的 `associations` 端点不返回 `standard_error` 字段，
管线用 `SE = |beta| / Φ⁻¹(1 − P/2)`（双侧正态）反推。

- 对绝大多数常见变异，该近似与真实 SE 差异可忽略；
- 但当 P 被上游截断（如报告为 0 或低于双精度下限）时会失真，
  管线已把 P 钳制在 1e-300，并要求 SE > 0；
- 影响方向：SE 轻微偏差会同时影响 IV 的 F 值与 IVW 权重，属**非系统性**噪声。

## 五、反向 MR 受限于脊柱感染 GWAS 本身的显著位点极少

FinnGen R11 报告的全基因组显著位点数：

| 表型 | 显著位点数 |
|---|---|
| M13_OSTEOMYELVERTEB（VOM） | **0** |
| M13_SPONDYLOINFECTION | **0** |
| M13_DISCITIS | 1 |
| M13_DISCINFECTION | 1 |
| M13_OSTEOMYELITIS | 2 |

多方法 MR（Egger/加权中位数/加权众数）至少需要 3 个 IV，**在 5e-8 阈值下反向 MR 不可行**。
`analysis/reverse_mr.py` 因此采用 **5e-6 次级阈值**并通过流式扫描远端文件抽取候选 IV
（M13_OSTEOMYELITIS 957 个变异 → 42 个独立 IV；M13_DISCITIS 79 → 45 个独立 IV，F>10）。
结果解释须极为谨慎：

- 放宽阈值会引入**弱工具变量偏倚**，方向偏向观察性关联；
- `results/reverse_MR_IV_discovery.csv` 记录了实际获得的 IV 数与平均 F；
- 反向臂的定位是**方向性佐证**，不足以独立支撑"无反向因果"的结论。

### 5.1 反向 MR 的执行层阻断：EBI API 速率限制（HTTP 429）

反向臂的**结局端**是免疫性状，只能走 EBI summary-statistics API 做定点查询
（FinnGen 只有疾病终点，没有细胞因子/血细胞计数）。执行中遇到两个接口事实：

1. **`variant_id` 查询参数被服务端忽略** —— 实测
   `/traits/{efo}/associations?variant_id=rsX` 返回的是该性状的前 N 条关联，
   而非目标变异。只有
   `/chromosomes/{chr}/associations?bp_lower=&bp_upper=` 才是可靠定点查询
   （坐标为 GRCh38，与 FinnGen R11 一致，已用 rs1205/CRP 位点验证）。
   早期按 `variant_id` 反查得到的"命中"不可信，相关缓存与中间结果已删除重做。
2. **速率限制**：并发查询触发 429，且进入**持续数十分钟的惩罚期**。
   已把请求数从 IV×性状（462 次）降到 IV 数（42 次，一次请求取回该位点上所有 study），
   并改为串行 + 3 秒间隔 + 指数退避，但惩罚期内首个请求仍被 429 拒绝。

`analysis/wait_and_run_reverse.py` 每 3 分钟轻量探测一次，解禁后自动续跑；
位点级抓取结果落盘于 `data/cache/revpos_*.json`，**可断点续跑**。
若 `results/MR_results_reverse.csv` 不存在或行数偏少，即表示该轮仍被限流所阻，
换网络环境或隔一段时间重跑 `python reverse_mr.py` 即可，无需改代码。

## 六、样本重叠与人群结构

- 暴露（EBI 汇总统计）多来自 UKB / INTERVAL / deCODE 等**欧洲人群**；
- 结局（FinnGen R11）为**芬兰人群**。
- 两者**基本无样本重叠**（这是优点，避免了重叠导致的偏倚放大）；
- 但芬兰人群存在**奠基者效应**，等位基因频率与一般欧洲人群有系统差异，
  可能影响效应量可移植性（transportability）。结论外推到非芬兰人群时需注意。

## 七、多重检验

正向设计为 11 个暴露 × 5 个结局 = 55 对，其中 **45 对**取到 ≥1 个可用 IV 而进入 IVW 检验
（其余 10 对因暴露 IV 在结局中全部缺失或等位基因不匹配而无结果）。
`results/MR_results_primary_IVW.csv` 已同时给出未校正 P、**FDR（Benjamini–Hochberg）**
与 **Bonferroni** 校正值（按实际 45 次检验计）。

校正后结论：

| 暴露 → 结局 | nSNP | OR (95% CI) | P | FDR | Bonferroni |
|---|---|---|---|---|---|
| WBC → 感染性脊柱病 | 28 | 0.037 (0.006–0.210) | 2.0e-4 | **0.009** | **0.009** |
| 中性粒细胞 → 感染性脊柱病 | 29 | 0.104 (0.019–0.570) | 9.1e-3 | 0.205 | 0.410 |
| CCL2 → 感染性脊柱病 | 7 | 0.290 (0.092–0.917) | 0.035 | 0.526 | 1.000 |

仅 **WBC → 感染性脊柱病**在 Bonferroni 与 FDR 下均保持显著；
另两项校正后不显著，只应视为**提示性（suggestive）**。
另需注意这三项的加权中位数/加权众数估计虽方向一致（OR<1）但均不显著，
提示效应可能由少数 IV 驱动，稳健性有限（见 `figures/` 中对应的留一法与漏斗图）。

## 八、其他未做的分析

- **MR-PRESSO**：需要 R 包 `MRPRESSO`，本环境不可用；已用 Cochran Q + 留一法 + 漏斗图
  作为水平多效性/离群点的替代检查。若需 PRESSO 的全局检验与离群校正，需另行实现。
- **多变量 MR（校正吸烟/肥胖/糖尿病）**：FEASIBILITY 第三节列为第 4 步，本次未做。
  管线已具备取数与协调能力，扩展为 MVMR 需增加多暴露联合回归模块。
- **cis/trans 区分**：细胞因子暴露（IL6R、CD40 等）未区分 cis-pQTL 与 trans-pQTL。
  cis 工具变量的多效性风险显著更低，建议后续按基因组位置拆分敏感性分析。
