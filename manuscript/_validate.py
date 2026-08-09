import re
with open(r'C:\Users\liouse\WorkBuddy\2026-08-07-22-20-21\spine_MR\manuscript\manuscript_SciReps.tex', encoding='utf-8') as f:
    t = f.read()

# 摘要
m = re.search(r'\\begin\{abstract\}(.+?)\\end\{abstract\}', t, re.DOTALL)
ab = m.group(1) if m else ''
ab = re.sub(r'\\textbf\{([^}]+)\}', r'\1', ab)
ab = re.sub(r'\\emph\{([^}]+)\}', r'\1', ab)
ab = re.sub(r'\\textsuperscript\{[^}]+\}', '', ab)
ab = re.sub(r'\\times', 'x', ab)
ab = re.sub(r'\\rightarrow', '->', ab)
ab = re.sub(r'\\%', '', ab)
ab = re.sub(r'\$', '', ab)
ab = re.sub(r'\\\\', ' ', ab)
words = ab.split()
print('Abstract word count:', len(words))
print('First 30:', ' '.join(words[:30]))

# 章节数
sects = re.findall(r'^\\section\*?\{([^}]+)\}', t, re.M)
print(f'Sections ({len(sects)}):', sects)

# 参考文献数（tex item）
refsec = re.search(r'\\section\*?\{References\}(.+?)\\section', t, re.DOTALL)
if refsec:
    items = re.findall(r'\\item\b', refsec.group(1))
    print(f'References: {len(items)}')

# 检查无中文字符（除已嵌入图注）
chinese = re.findall(r'[\u4e00-\u9fff]', t)
print(f'Chinese chars in tex: {len(chinese)}', chinese[:10] if chinese else '')

# 占位符残留检查
import re as r2
placeholders = ['full citation to be inserted', '[repository URL', 'To be finalized', 'To be added',
                '[Affiliations: author names', 'author names, departments', '[To be added]']
for ph in placeholders:
    found = ph in t
    print(f'placeholder "{ph}": {"FOUND" if found else "removed"}')

# Figure 文件引用
figs = re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', t)
print(f'Figures embedded ({len(figs)}):', figs)