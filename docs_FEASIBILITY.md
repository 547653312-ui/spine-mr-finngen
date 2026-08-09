# MR 路线 GWAS 数据可用性确认（脊柱感染因果推断）

> 结论：**可行，且脊柱特异性暴露直接存在。** 转录组路线卡死（无配对数据集），MR 路线反而更好落地。

## 一、暴露端（脊柱感染相关 GWAS 汇总数据）

来源：GWAS Catalog（MONDO_0005246，34 个研究）+ 2025 年 OM MR 论文（MD, Ovid, 2025-09）披露的精确 ID。

| 表型 | 数据源 | 病例/对照 | GWAS ID | 脊柱特异性 |
|---|---|---|---|---|
| 骨髓炎 OM（ICD-10 M86） | FinnGen | 2336 / 473,264 | — | 含椎体 |
| 骨髓炎 OM | UKB–IEU | 4836 / 481,648 | **ieu-b-4975** | 含椎体（可自动下载） |
| **椎体骨髓炎 VOM（M46.2）** | FinnGen | 111 / 353,224 | — | ★脊柱特异（需 FinnGen token） |
| **椎间盘炎 DC（M46.4）** | FinnGen | 557 / 353,224 | — | ★脊柱特异（需 FinnGen token） |
| 感染性心内膜炎 IE | UKB–IEU | 1080 / 485,404 | ieu-b-49720 | 对照/验证 |

GWAS Catalog 还列出 PheCode 细分：Acute/Chronic osteomyelitis (710.11/710.12)、Unspecified (710.19)——均可在 IEU/Finngen 取到。

## 二、结局端（现成、可自动下载）

- **免疫细胞性状**：`met-b`（150 项免疫子集频率，Roederer 2015）、`ieu-a` 系列免疫性状。
- **细胞因子/免疫蛋白**：`prot-a`（Sun 2018 免疫蛋白，如 IL-6、TNF-α、CRP 等）。
- **骨代谢/骨密度/骨折**：`ukb-b` 系列、`GEFOS` 骨密度 GWAS（EBI 可下）。
- **混杂校正**：吸烟 ukb-b-223、肥胖 ukb-b-15541、糖尿病 ukb-b-10753。

## 三、推荐研究设计（由 kimi-K3 规划细化）

1. **主分析（全自动，IEU 可下载）**：暴露 = 免疫/细胞因子性状；结局 = **ieu-b-4975 骨髓炎** → 识别"因果驱动脊柱感染的免疫组分"。
2. **反向 MR**：ieu-b-4975 → 免疫性状，验证方向。
3. **脊柱特异性确认（需 FinnGen token 或手动下载）**：VOM(M46.2)/DC(M46.4) ↔ 免疫性状，作为核心 spine framing。
4. **多变量 MR**：校正吸烟/肥胖/糖尿病。

## 四、技术栈（R 不可用 → Python）

- 取数：IEU REST API（`https://gwas.mrcieu.ac.uk/api/`），Python `requests`。
- MR 核心：`MendelianRandomization`(PyPI)、`pandas`、`scipy`、`numpy`。
- 方法：IVW / MR-Egger / 加权中位数 / 加权众数；敏感性：Cochran Q、Egger 截距、留一法、MR-PRESSO；多变量 MR。
- 作图：`matplotlib` / `seaborn`。

## 五、可行性判定

- 暴露 ✅（OM 在 IEU 可下载；VOM/DC 在 FinnGen 已发表）
- 结局 ✅（免疫/细胞因子/骨密度大量现成）
- 新颖性 ✅（"脊柱感染因果免疫驱动"在 MR 领域几乎空白；2025 仅见 OM↔IE，未见 VOM/DC↔免疫）
- 风险：FinnGen VOM/DC 需 token（自动化受限，可用 OM 主体 + 文献佐证脊柱亚型）；小样本 VOM(111) 需敏感性分析兜底。

下一步：kimi-K3 出详细方案 → Hy3 搭建 Python 管线并执行 → kimi-K3 验证并指挥修订。
