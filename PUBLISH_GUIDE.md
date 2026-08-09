# 发布指南：GitHub + Zenodo DOI

本仓库已经在本地完成 `git init` 和首次提交，只差「推到 GitHub」和「Zenodo 归档」两步。
全程免费，约 15 分钟。

---

## 第 0 步：改成你自己的署名（1 分钟）

首次提交用的是占位身份 `Anonymous Author <author@example.com>`。换成你自己的：

```bash
cd spine-mr-finngen
git config user.name  "你的名字拼音"
git config user.email "你的邮箱@example.com"
git commit --amend --reset-author --no-edit
```

同时把这三个文件里的占位符换成真实信息：

| 文件 | 占位符 |
|------|--------|
| `LICENSE` | `[AUTHOR NAME]` |
| `CITATION.cff` | `[FAMILY NAME]` / `[GIVEN NAME]` / ORCID / `[AFFILIATION]` / `<USERNAME>` |
| `.zenodo.json` | 同上 |
| `README.md` | `<USERNAME>` |

---

## 第 1 步：在 GitHub 建仓库（3 分钟）

1. 登录 https://github.com → 右上角 **+** → **New repository**
2. 填写：
   - Repository name: `spine-mr-finngen`
   - Description: `Mendelian randomization of immune traits and spine-specific osteomyelitis (FinnGen R11)`
   - 选 **Public**（Zenodo 只能归档公开仓库）
   - **不要**勾选 Add README / .gitignore / license（本地已有，勾了会冲突）
3. 点 **Create repository**

---

## 第 2 步：推送代码（2 分钟）

在 GitHub 新仓库页面复制你的仓库地址，然后：

```bash
cd spine-mr-finngen
git remote add origin https://github.com/<你的用户名>/spine-mr-finngen.git
git push -u origin main
```

首次推送会弹出登录窗口，用 GitHub 账号授权即可（或用 Personal Access Token 作为密码）。

推送完成后打开仓库页面确认 README 正常显示。

---

## 第 3 步：把 Zenodo 和 GitHub 连起来（3 分钟）

1. 打开 https://zenodo.org → **Log in** → 选择 **Log in with GitHub**（用 GitHub 账号直接登录，最省事）
2. 授权后，点右上角头像 → **GitHub**（或直接访问 https://zenodo.org/account/settings/github/）
3. 在仓库列表里找到 `spine-mr-finngen`，把右侧开关 **切换成 ON**

> 开关必须在发布 Release **之前**打开，否则 Zenodo 收不到通知。

---

## 第 4 步：发一个 Release，自动拿到 DOI（3 分钟）

1. 回到 GitHub 仓库页面 → 右侧 **Releases** → **Create a new release**
2. 填写：
   - **Choose a tag** → 输入 `v1.0.0` → 点 "Create new tag: v1.0.0 on publish"
   - Release title: `v1.0.0 — Manuscript submission release`
   - Description（可直接复制）：
     ```
     Code, harmonised instrument tables, results and figures accompanying the
     manuscript "No robust causal effect of genetically predicted immune traits
     on spine-specific osteomyelitis: a two-sample Mendelian randomization study".
     ```
3. 点 **Publish release**
4. 等 1–2 分钟，回到 https://zenodo.org/account/settings/github/ ，`spine-mr-finngen`
   右侧会出现一个 **DOI 徽章**，形如：

   ```
   10.5281/zenodo.15873201
   ```

> Zenodo 会给两个 DOI：
> - **Concept DOI**（永远指向最新版本）—— 论文里用这个
> - **Version DOI**（固定指向 v1.0.0）
>
> 在 Zenodo 记录页右侧 "Cite all versions?" 处能看到 Concept DOI，**论文用 Concept DOI**。

---

## 第 5 步：回填到手稿（2 分钟）

打开 `spine_MR/manuscript/build_manuscript.py`，找到 `Code availability` 段落，
把两个尖括号占位符替换掉：

| 占位符 | 替换为 |
|--------|--------|
| `<GITHUB_USERNAME>` | 你的 GitHub 用户名 |
| `<ZENODO_RECORD_ID>` | Concept DOI 的数字部分，例如 `15873201` |

然后重新生成稿件：

```bash
cd spine_MR/analysis
./venv/Scripts/python.exe ../manuscript/build_manuscript.py
```

替换后的声明会长这样：

> The code is deposited at https://github.com/zhangsan/spine-mr-finngen and archived
> with a persistent identifier at Zenodo (DOI: https://doi.org/10.5281/zenodo.15873201).

顺手在 README 顶部加个徽章（可选）：

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15873201.svg)](https://doi.org/10.5281/zenodo.15873201)
```

---

## 常见问题

**Q：仓库里有 12.7 MB，会不会太大？**
不会。GitHub 单仓库软上限 1 GB，单文件 100 MB。已通过 `.gitignore` 排除了
`data/cache/`（13.7 MB）和 `data/raw/`（28.9 MB）这些可重新下载的 GWAS 原始数据。

**Q：稿件 docx 也在仓库里，投稿前公开会不会影响双盲评审？**
Scientific Reports 是单盲，不影响。若你介意，可在推送前删掉
`manuscript/manuscript_SciReps.docx` 和 `.tex`，只保留 `.md` 源文件和 build 脚本。

**Q：论文被拒了要转投，DOI 还能用吗？**
能。Zenodo DOI 与期刊无关，永久有效。后续如果改了代码，再发一个 `v1.1.0` Release，
Zenodo 会自动生成新的 Version DOI，而 Concept DOI 始终指向最新版，论文里不用改。

**Q：不想公开代码怎么办？**
Scientific Reports 要求 code availability 声明必须真实可核查。可以选择 Zenodo
"Restricted access" 上传（不走 GitHub），审稿人凭链接申请访问；但完全不公开会在
编辑初审阶段被要求整改。建议直接公开。
