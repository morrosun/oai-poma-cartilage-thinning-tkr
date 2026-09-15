# -*- coding: utf-8 -*-
"""摸清两套数据的可用维度：CSV 轨迹结构、ID 是否可跨库连接、掩膜标签取值。"""
import os, json, gzip
import numpy as np
import pandas as pd

ROOT = r'D:\OAI_datasets'
out = []
def p(*a):
    s = ' '.join(str(x) for x in a); out.append(s)

# ---------- 1. CTh-Score 表 ----------
csv = os.path.join(ROOT, 'CTh-Maps', 'OAI_CTh-Maps_CTh-Score_Dataset_v1', 'df_cth_score_v1.csv')
df = pd.read_csv(csv)
p('== df_cth_score_v1.csv ==')
p('shape =', df.shape, '| columns =', list(df.columns))
p('unique ID        =', df['ID'].nunique())
p('unique (ID,SIDE) =', df.groupby(['ID', 'SIDE']).ngroups)
p('dtypes:', dict(df.dtypes.astype(str)))
p('SIDE counts:', df['SIDE'].value_counts().to_dict())
p('Version:', df['Version'].unique().tolist())

order = ['00m', '12m', '24m', '36m', '48m', '72m', '96m']
p()
p('== CTh-Score by timepoint ==')
p('%-6s %7s %8s %8s %8s' % ('TP', 'n', 'mean', 'median', 'sd'))
for tp in order:
    s = df.loc[df['Timepoint'] == tp, 'CTh-Score']
    if len(s):
        p('%-6s %7d %8.2f %8.2f %8.2f' % (tp, len(s), s.mean(), s.median(), s.std()))

# 每膝随访点数分布
cnt = df.groupby(['ID', 'SIDE'])['Timepoint'].nunique()
p()
p('== 每膝可得时点数分布 ==')
vc = cnt.value_counts().sort_index()
for k, v in vc.items():
    p('  %d 个时点: %6d 膝' % (k, v))
p('  >=6 时点: %d 膝  |  ==7 时点: %d 膝' % ((cnt >= 6).sum(), (cnt == 7).sum()))

# 每膝斜率（简单 OLS，单位：分/年）
tmap = {'00m': 0, '12m': 1, '24m': 2, '36m': 3, '48m': 4, '72m': 6, '96m': 8}
rows = []
for (i, sd), g in df.groupby(['ID', 'SIDE']):
    g = g[g['Timepoint'].isin(order)]
    if len(g) < 3:
        continue
    x = np.array([tmap[t] for t in g['Timepoint']], float)
    y = g['CTh-Score'].to_numpy(float)
    b = np.polyfit(x, y, 1)[0]
    rows.append((i, sd, len(g), y[0], y[-1], b))
sl = pd.DataFrame(rows, columns=['ID', 'SIDE', 'n_tp', 'first', 'last', 'slope_per_yr'])
p()
p('== 每膝 CTh-Score 斜率（分/年，n>=3 时点，n=%d）==' % len(sl))
p('  mean %.3f  median %.3f  sd %.3f  min %.3f  max %.3f' % (
    sl['slope_per_yr'].mean(), sl['slope_per_yr'].median(), sl['slope_per_yr'].std(),
    sl['slope_per_yr'].min(), sl['slope_per_yr'].max()))
p('  斜率>0（恶化）占比 %.1f%%  |  <0（改善）占比 %.1f%%  |  ==0 %.1f%%' % (
    (sl['slope_per_yr'] > 0).mean()*100, (sl['slope_per_yr'] < 0).mean()*100,
    (sl['slope_per_yr'] == 0).mean()*100))
# 左右相关
w = sl.pivot_table(index='ID', columns='SIDE', values='slope_per_yr').dropna()
p('  同期左右膝斜率 Pearson r = %.3f (n=%d 人)' % (w['LEFT'].corr(w['RIGHT']), len(w)))

# ---------- 2. OAIZIB-CM 的 ID 映射 ----------
info = os.path.join(ROOT, 'OAIZIB-CM', 'info')
p()
p('== OAIZIB-CM info 文件 ==')
for f in sorted(os.listdir(info)):
    p('  ', f, os.path.getsize(os.path.join(info, f)), 'B')

p()
p('-- kneeSideInfo.csv --')
ks = pd.read_csv(os.path.join(info, 'kneeSideInfo.csv'))
p('shape =', ks.shape, '| columns =', list(ks.columns))
p(ks.head(5).to_string())

for f in ['subInfo_train.xlsx', 'subInfo_test.xlsx']:
    fp = os.path.join(info, f)
    try:
        x = pd.read_excel(fp)
        p()
        p('-- %s --' % f, 'shape =', x.shape)
        p('columns =', list(x.columns))
        p(x.head(5).to_string())
    except Exception as e:
        p('-- %s -- READ FAIL %s' % (f, e))

# ---------- 3. 掩膜标签取值 ----------
lab = os.path.join(ROOT, 'OAIZIB-CM', 'labelsTr', 'oaizib_001.nii.gz')
with gzip.open(lab, 'rb') as g:
    raw = g.read()
arr = np.frombuffer(raw[352:], dtype=np.uint8).reshape(160, 384, 384, order='F')
u, c = np.unique(arr, return_counts=True)
p()
p('== labelsTr/oaizib_001.nii.gz 标签取值 ==')
for a, b in zip(u, c):
    p('  label %-3d  voxels %10d  (%.3f%%)' % (a, b, b/arr.size*100))

txt = '\n'.join(out)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
open(os.path.join(OUTDIR, 'data_capability.txt'), 'w', encoding='utf-8').write(txt)
print(txt)
