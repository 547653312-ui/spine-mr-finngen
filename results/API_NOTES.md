# 数据源连通性记录（API_NOTES）

探测脚本：`analysis/test_api.py` ｜ 原始结果：`results/api_probe_raw.json`
探测时间：2026-08-08 ｜ 环境：Windows 11，无代理，R 不可用

**一句话结论：IEU/OpenGWAS 已对全部数据端点强制 JWT，`ieu-b-4975` 无法自动获取；
改用两个完全开放的等价数据源后，管线全自动跑通，且拿到了比原方案更好的脊柱特异结局。**

---

## 一、IEU / OpenGWAS —— 不可用（硬阻塞）

| 端点 | 状态 | 返回内容 |
|---|---|---|
| `https://gwas.mrcieu.ac.uk/api/gwasinfo/ieu-b-4975` | **301** | 永久重定向（openresty） |
| 同上，跟随重定向后 | **404** | OpenGWAS 站点 HTML「Resource not found」 |
| `https://api.opengwas.io/api/status` | 200 | `API__VERSION 4.0.0`，`SERVICES__METADATA: Operational` |
| `https://api.opengwas.io/api/gwasinfo/ieu-b-4975` | **404** | Flask「The requested URL was not found」 |
| `.../gwasinfo/met-b-1`、`met-b-10` | **404** | 同上 |
| `.../gwasinfo/prot-a-1`、`prot-a-670` | **404** | 同上 |
| `.../gwasinfo/ieu-b-49720` | **404** | 同上 |
| `https://api.opengwas.io/api/gwasinfo`（列表） | **401** | 见下方原文 |
| `.../associations/ieu-b-4975/rs1205` | **404** | 同上 |

401 原文（关键证据）：

```
{"message": "ERROR - Go to https://api.opengwas.io/ - From 1st May 2024 you must
provide a token (JWT) alongside most of your requests. Read more at
https://api.opengwas.io/ and also check for the latest version at
https://mrcieu.github.io/ieugwasr/"}
```

补充核查：
- 本机环境变量、`~/.opengwas*`、`~/.Renviron` 中**均无** OPENGWAS_JWT，无法免注册获取。
- 服务本身是 Operational 的（status 200），因此 **404/401 属于鉴权与路由策略，不是网络封锁**。
- 结论：`ieu-b-4975`、`met-b-*`、`prot-a-*` 这条路线在无 token 情况下**不可自动化**。

## 二、EBI GWAS Catalog —— 可用（暴露端采用）

| 端点 | 状态 | 说明 |
|---|---|---|
| `rest/api/efoTraits/MONDO_0005246` | 200 | 骨髓炎性状本体，确认 34 项研究 |
| `summary-statistics/api/` | 200 | 汇总统计 API 根，无需 token |
| `summary-statistics/api/traits/{EFO}/associations?p_upper=5e-8` | 200 | **按 P 值过滤直接返回 IV** |

返回字段完全满足 MR 需要：
`variant_id(rsID) / chromosome / base_pair_location(GRCh38) / effect_allele /
other_allele / beta / odds_ratio / effect_allele_frequency / p_value / study_accession`

注：该端点不返回 SE，管线由 `beta` 与双侧 P 反推：`SE = |beta| / Φ⁻¹(1−P/2)`
（`analysis/gwas_io.py: se_from_beta_p`）。

已验证有汇总统计的免疫/细胞因子 EFO（11 个，见 `mr_pipeline.py: EXPOSURES`）：
CRP `EFO_0004458`、IL6R `EFO_0008187`、IL-1β `EFO_0004812`、CXCL10 `EFO_0008056`、
CD40 `EFO_0010607`、CCL2 `EFO_0004749`、中性粒 `EFO_0004833`、淋巴 `EFO_0004587`、
单核 `EFO_0005091`、嗜酸 `EFO_0004842`、白细胞总数 `EFO_0004308`。

失败记录：`EFO_0004910`（拟用作 IL-6）**404**，该号在汇总统计库中不存在；
经 OLS 核对后改用 IL-6 受体 α（`EFO_0008187`，即经典 IL6R Asp358Ala 通路工具变量）。

### 2.1 两个必须知道的接口陷阱（反向 MR 时踩到）

**陷阱 A：`variant_id` 查询参数被服务端静默忽略。**

```
GET /summary-statistics/api/traits/EFO_0004458/associations?variant_id=rs2794520&size=5
→ 200，但返回的是该性状的前 5 条关联，rs2794520 并不在其中
GET /summary-statistics/api/associations/rs2794520?trait=EFO_0004458   → 200，n=0
GET /summary-statistics/api/associations/rs2794520?study_accession=GCST008055 → 200，n=0
```

唯一可靠的定点查询是按坐标：

```
GET /summary-statistics/api/chromosomes/1/associations?bp_lower=159712443&bp_upper=159712443
→ 200，命中 rs1205（CRP 位点），beta=-0.195, P=4.6e-109
同一位点用 GRCh37 坐标 159682233 → 404
```

即 **EBI 汇总统计坐标为 GRCh38**，与 FinnGen R11 一致，两边可直接按 chr:pos 对接。

**陷阱 B：速率限制（HTTP 429）且惩罚期很长。**
并发（12 线程）定点查询会迅速触发 429，且之后**数十分钟内**连单个请求也被拒。
应对：
- 不加 `study_accession` 过滤，一次请求取回该位点上**所有 study**，
  于是 11 个性状共用一次请求（请求数 462 → 42）；
- 串行 + 3 秒固定间隔 + 指数退避，连续 429 即中止而不是继续打接口；
- 位点级结果落盘 `data/cache/revpos_*.json`，断点续跑。
见 `analysis/reverse_mr.py` 与 `analysis/wait_and_run_reverse.py`。

## 三、FinnGen —— 可用（结局端采用，且优于原方案）

| 端点 | 状态 | 说明 |
|---|---|---|
| `https://r12.finngen.fi/api/phenos` | 200 | 2470 个表型元数据（1.16 MB） |
| `https://r11.finngen.fi/api/pheno/{code}` | 200 | 单表型元数据 |
| R12 bucket `finngen-public-data-r12/...` | **404** | R12 汇总统计需 DUA/授权，非公开 |
| R11 bucket `.gz` | **200** | 公开，765–808 MB/表型 |
| R11 bucket `.gz.tbi` | **200** | **公开提供 tabix 索引（关键）** |

**关键突破**：R11 同时公开 bgzip 数据与 tabix 索引，因此不必下载 800MB 全文件。
`analysis/gwas_io.py` 实现了纯 Python 的 tabix 索引解析 + BGZF 分块解压 + HTTP Range
随机访问（无需 pysam/samtools，Windows 友好）：

```
实测：打开 765MB 远端文件并查询 chr1:159684665，命中 2 行，耗时 3.4s（首次含索引下载）
并发优化后（Session 复用 + 12 线程）：0.16 s/SNP，较初版提速约 80 倍
```

更重要的是，R11 直接提供 **FEASIBILITY.md 中原本"需 FinnGen token"的脊柱特异表型**，
且病例数与文档记载完全吻合：

| phenocode | 表型 | 病例/对照 |
|---|---|---|
| `M13_OSTEOMYELVERTEB` | 椎体骨髓炎（M46.2） | **111** / 353,224 |
| `M13_DISCITIS` | 椎间盘炎（M46.4） | **557** / 353,224 |
| `M13_DISCINFECTION` | 化脓性椎间盘感染 | 423 / 353,224 |
| `M13_SPONDYLOINFECTION` | 其他感染性脊柱病 | 74 / 353,224 |
| `M13_OSTEOMYELITIS` | 骨髓炎（全部位，M86） | **2336** / 473,264 |

## 四、最终采用的数据路线

```
暴露：EBI GWAS Catalog Summary Statistics API   （免疫/细胞因子，P<5e-8）
                    ↓  chr:pos GRCh38 匹配 + 等位基因协调
结局：FinnGen R11 公开汇总统计（tabix over HTTP Range）
```

坐标系一致性已验证：EBI 的 `rs2228145` 定位 chr1:154454494，
FinnGen R11 同位点 `nearest_genes = IL6R`，等位基因 A/C 一致 —— 两端均为 GRCh38，可直接匹配。

## 五、若需回到 IEU 路线：手动下载兜底方案

以下步骤可在拿到 token 或人工下载后，无需改动分析代码即可接入。

**方案 A：申请 OpenGWAS token（推荐，10 分钟）**
1. 访问 https://api.opengwas.io/ 用学术邮箱注册，在 profile 页生成 JWT。
2. 设置环境变量：`setx OPENGWAS_JWT "<token>"`。
3. 在 `gwas_io.py` 的 `HEADERS` 中加入 `{"Authorization": "Bearer " + os.environ["OPENGWAS_JWT"]}`，
   即可启用 `probe_opengwas` 之外的 `/tophits/{id}`、`/associations/{id}/{rsid}` 端点。

**方案 B：纯手工下载 VCF（无需 token）**
1. 打开 https://gwas.mrcieu.ac.uk/datasets/ieu-b-4975/ ，点击 "Download" 取
   `ieu-b-4975.vcf.gz` 与 `.vcf.gz.tbi`，放入 `data/manual/`。
2. IEU 的 GWAS-VCF 同样是 bgzip+tabix，本仓库的 `BGZFRemote`/`TabixIndex` 可直接复用，
   只需把 URL 换成本地路径（`open(path,'rb').seek/read` 替换 `read_range`）。
3. FORMAT 字段顺序为 `ES:SE:LP:AF:SS:ID`（ES=beta，LP=-log10 P），
   解析后即可喂给 `mr_methods.run_all`。

**方案 C：FinnGen R12（含更大样本的脊柱表型）**
1. 在 https://www.finngen.fi/en/access_results 提交数据使用申请获取下载授权。
2. 下载 `finngen_R12_M13_OSTEOMYELVERTEB.gz(+.tbi)` 等文件到 `data/manual/`。
3. 把 `gwas_io.FinnGenSumstats` 的 `release` 参数改为 `"R12"` 并指向本地路径。

> 本次分析未使用任何上述手工数据 —— 全流程通过开放端点自动完成。
