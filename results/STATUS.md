# STATUS — Hy3（Executor）交付状态

范围：`analysis/` `data/` `results/`　｜　依据：`FEASIBILITY.md`（`plan/PLAN.md` 不存在，按 DEFAULT DESIGN 执行并做了一处升级，见第 3 节）

---

## 1. 环境

| 项 | 状态 |
|---|---|
| 托管 Python | `C:\Users\liouse\.workbuddy\binaries\python\versions\3.13.12\python.exe`（实际 3.13.14）✅ |
| venv | `analysis/venv` ✅ |
| pandas 3.0.5 / numpy 2.5.1 / scipy 1.18.0 / requests 2.34.2 / matplotlib 3.11.1 / seaborn 0.13.2 | ✅ |
| statsmodels | 已装（备用）✅ |
| **MendelianRandomization** | ❌ **PyPI 上不存在该包**（它是 R 包）。备选 `genal` 亦不可用。已自研 `analysis/mr_methods.py` 复现全部方法 |
| R / TwoSampleMR / MRPRESSO / plink | ❌ 不可用（任务书已声明） |

pip 自升级在 venv 内触发 safe-delete 报错，改用 `--no-cache-dir` 直装，不影响结果。

## 2. 数据接口可达性（探测：`analysis/test_api.py`，原始输出 `api_probe_raw.json`，24 端点 / 12 可达）

| 数据源 | 结论 |
|---|---|
| OpenGWAS / IEU REST（含 `ieu-b-4975`、`met-b-*`、`prot-a-*`） | ❌ 401/404，2024-05-01 起强制 JWT，**无法自动化** |
| EBI GWAS Catalog Summary Statistics | ✅ 开放（暴露端）。两处接口陷阱见 `API_NOTES.md` §2.1 |
| FinnGen R12 bucket | ❌ 404，需 DUA |
| FinnGen R11 bucket（`.gz` + `.gz.tbi`） | ✅ 完全开放（结局端） |

## 3. 设计：相对 DEFAULT DESIGN 的一处升级

结局由 `ieu-b-4975`（UKB 混合骨髓炎）改为 **FinnGen R11 五个脊柱/骨髓炎表型**：

| 缩写 | phenocode | R11 病例 |
|---|---|---|
| VOM（椎体骨髓炎 M46.2） | `M13_OSTEOMYELVERTEB` | 111 |
| DISCITIS（椎间盘炎 M46.4） | `M13_DISCITIS` | 557 |
| DISCINF（化脓性椎间盘感染） | `M13_DISCINFECTION` | 423 |
| SPONDINF（感染性脊柱病） | `M13_SPONDYLOINFECTION` | 74 |
| OM（骨髓炎全部位 M86） | `M13_OSTEOMYELITIS` | 2336 |

理由：`ieu-b-4975` 需 token 不可自动获取；且 FinnGen R11 直接给出 **M46.2 / M46.4**，
恰好补上 FEASIBILITY.md 里"需 token 才能拿到 VOM/DC"的缺口，脊柱特异性严格优于 UKB 代理。
暴露：11 个免疫/细胞因子性状（EFO 逐个校验可用）。

## 4. 已完成

**正向 MR（免疫 → 脊柱感染/骨髓炎）：完整跑通，5 个结局全覆盖。**

- `MR_results_all.csv` —— 235 行（IVW 固定/随机、MR-Egger、加权中位数、加权众数）
- `MR_results_primary_IVW.csv` —— 45 对，含 **FDR(BH)** 与 **Bonferroni** 校正列
- `figures/` 196 张（scatter / funnel / leave-one-out / forest）
- `harmonised/` 63 份协调后 IV 明细、`loo/` 100 份留一法明细
- `reverse_MR_feasibility.csv`、`reverse_MR_IV_discovery.csv`
- 日志：`run_log.txt`（主跑）、`run_log_vom.txt`（VOM 补跑）、`run_log_reverse*.txt`

### 头条结果（IVW 随机效应）

| 暴露 → 结局 | nSNP | OR (95% CI) | P | FDR | Bonferroni |
|---|---|---|---|---|---|
| **白细胞总数 → 感染性脊柱病** | 28 | **0.037 (0.006–0.210)** | 2.0e-4 | **0.009** | **0.009** |
| 中性粒细胞 → 感染性脊柱病 | 29 | 0.104 (0.019–0.570) | 9.1e-3 | 0.205 | 0.410 |
| CCL2 → 感染性脊柱病 | 7 | 0.290 (0.092–0.917) | 0.035 | 0.526 | 1.000 |

三者方向一致（OR<1，遗传预测的免疫细胞水平升高 → 感染性脊柱病风险降低），
Cochran Q 无异质性（P=0.34–0.95），Egger 截距无多效性证据（P=0.39–0.97）；
但加权中位数/众数均不显著，稳健性有限。**VOM 与其余结局无显著结果**（VOM 仅 111 例，效能不足）。

## 5. 未完成 / 阻断

**反向 MR（脊柱感染 → 免疫）：IV 已备好，结局端查询被 EBI 限流阻断。**

- 5e-8 阈值下五个表型显著位点均 <3 个 → 多方法反向 MR 本就不可行（`reverse_MR_feasibility.csv`）；
- 退到 5e-6：流式扫描远端 789MB 文件，OM 得 42 个、DISCITIS 得 45 个独立 IV（F>10），
  已存 `data/ivs_reverse_*.csv`；
- 结局端（免疫性状）只能走 EBI API，密集查询后触发 **HTTP 429 长时间惩罚期**，
  期间连单个请求都被拒。已把请求数从 462 降到 42、改串行 3 秒间隔 + 指数退避，仍未解禁。
- `analysis/wait_and_run_reverse.py` 每 3 分钟探测一次，解禁后自动续跑；
  位点级抓取落盘 `data/cache/revpos_*.json`，**可断点续跑**。
- 重要更正：早期按 `variant_id` 反查得到的命中不可信（该参数被服务端忽略），
  相关缓存与中间结果**已全部删除重做**，未混入任何已发布结果。

其他阻断：OpenGWAS 需 JWT（`API_NOTES.md` 附手工下载兜底）；FinnGen R12 需 DUA；
无 plink/LD 面板，以 ±500 kb 距离剪枝近似独立性。完整清单见 `LIMITATIONS.md`。

## 6. 复跑方式

```bash
cd analysis
./venv/Scripts/python.exe mr_pipeline.py --direction forward        # 正向全量
./venv/Scripts/python.exe mr_pipeline.py --direction forward --outcomes VOM --out-suffix _VOM
./venv/Scripts/python.exe reverse_mr.py --pthresh 5e-6 --pause 3.0  # 反向（需 EBI 未限流）
```
