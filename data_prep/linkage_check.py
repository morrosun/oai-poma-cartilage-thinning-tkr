# -*- coding: utf-8 -*-
"""关键可行性检验：OAIZIB-CM(507 膝, 有 3D MRI+掩膜+KL) 能否与 CTh-Maps(4340 人轨迹) 按 ID 接上。"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\OAI_datasets'
out = []
def p(*a):
    s = ' '.join(str(x) for x in a); out.append(s)

# --- 1. side 映射（该 CSV 无表头，第一行即数据） ---
ks = pd.read_csv(os.path.join(ROOT, 'OAIZIB-CM', 'info', 'kneeSideInfo.csv'), header=None,
                 names=['file', 'side'])
p('== kneeSideInfo.csv ==')
p('shape=', ks.shape, '| side 取值:', ks['side'].value_counts().to_dict())
ks['cmt'] = ks['file'].str.extract(r'oaizib_(\d+)\.nii\.gz').astype(int)

# --- 2. subInfo 两份 ---
frames = []
for f in ['subInfo_train.xlsx', 'subInfo_test.xlsx']:
    x = pd.read_excel(os.path.join(ROOT, 'OAIZIB-CM', 'info', f))
    x['split'] = f
    frames.append(x)
sub = pd.concat(frames, ignore_index=True)
p()
p('== subInfo (train+test) ==')
p('shape=', sub.shape, '| 列:', list(sub.columns))
p('SubjectID 唯一数 =', sub['SubjectID'].nunique(), '| CMT-ID 范围 %d-%d'
  % (sub['CMT-ID'].min(), sub['CMT-ID'].max()))
p('KneeSide 取值:', sub['KneeSide'].value_counts().to_dict())
p('KLGrade 分布:', sub['KLGrade'].value_counts(dropna=False).sort_index().to_dict())
p('Gender 分布:', sub['Gender'].value_counts().to_dict())
p('Age  %.1f ± %.1f (%.0f-%.0f)' % (sub['Age'].mean(), sub['Age'].std(), sub['Age'].min(), sub['Age'].max()))
p('BMI  %.1f ± %.1f (%.1f-%.1f)' % (sub['BMI'].mean(), sub['BMI'].std(), sub['BMI'].min(), sub['BMI'].max()))

# KneeSide(1/2) -> left/right：用 kneeSideInfo 交叉判定
m = sub.merge(ks[['cmt', 'side']], left_on='CMT-ID', right_on='cmt', how='left')
p()
p('== KneeSide 编码判定（KneeSide x side 交叉表）==')
p(pd.crosstab(m['KneeSide'], m['side']).to_string())
p('未匹配 side 的条数 =', m['side'].isna().sum())

# --- 3. 与 CTh-Maps 的交集 ---
df = pd.read_csv(os.path.join(ROOT, 'CTh-Maps', 'OAI_CTh-Maps_CTh-Score_Dataset_v1', 'df_cth_score_v1.csv'))
cth_keys = set(zip(df['ID'], df['SIDE']))

# 归一化：side 大写；KneeSide 1->? 用交叉表结论决定
ct = pd.crosstab(m['KneeSide'], m['side'])
side_map = {}
for k in ct.index:
    side_map[k] = ct.loc[k].idxmax()
p()
p('KneeSide 映射推定：', side_map)
m['SIDE'] = m['KneeSide'].map(side_map).str.upper()
m['key'] = list(zip(m['SubjectID'], m['SIDE']))
m['in_cth'] = m['key'].isin(cth_keys)

p()
p('== 交集 ==')
p('OAIZIB-CM 膝数 =', len(m), '| 其中出现在 CTh-Maps 的 =', int(m['in_cth'].sum()),
  '(%.1f%%)' % (m['in_cth'].mean()*100))
p('涉及受试者 =', m.loc[m['in_cth'], 'SubjectID'].nunique())

hit = m[m['in_cth']].copy()
sub_cth = df[df.set_index(['ID', 'SIDE']).index.isin(hit['key'])].copy()
cov = sub_cth.groupby(['ID', 'SIDE'])['Timepoint'].nunique()
p()
p('命中膝的随访覆盖：')
p('  有 >=3 时点: %d  |  >=5 时点: %d  |  7 时点: %d'
  % ((cov >= 3).sum(), (cov >= 5).sum(), (cov == 7).sum()))
p('  总膝-时点记录 =', len(sub_cth))

# --- 4. 交叉校验：00m CTh-Score 对 KLGrade ---
base = sub_cth[sub_cth['Timepoint'] == '00m'].merge(
    hit[['SubjectID', 'SIDE', 'KLGrade', 'Age', 'BMI', 'Gender', 'CMT-ID']],
    left_on=['ID', 'SIDE'], right_on=['SubjectID', 'SIDE'], how='inner')
p()
p('== 基线(00m) CTh-Score 分层于 KL 分级 ==')
if len(base):
    g = base.groupby('KLGrade')['CTh-Score']
    p(g.agg(['count', 'mean', 'median', 'std']).to_string())
    b = base.dropna(subset=['KLGrade'])
    if b['KLGrade'].nunique() > 1:
        p('  Spearman r(KLGrade, CTh-Score) = %.3f (n=%d)'
          % (b['KLGrade'].corr(b['CTh-Score'], method='spearman'), len(b)))
        p('  Pearson  r(BMI,      CTh-Score) = %.3f'
          % (b['BMI'].corr(b['CTh-Score'])))
        p('  Pearson  r(Age,      CTh-Score) = %.3f'
          % (b['Age'].corr(b['CTh-Score'])))
else:
    p('  无重叠样本')

txt = '\n'.join(out)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
open(os.path.join(OUTDIR, 'linkage_check.txt'), 'w', encoding='utf-8').write(txt)
print(txt)
