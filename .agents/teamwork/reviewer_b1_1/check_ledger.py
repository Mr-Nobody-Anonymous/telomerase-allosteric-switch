import io
import csv
import pandas as pd

with open('results/audit/evidence_ledger.csv', 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
for l in lines:
    if 'report.md (Section 1, 4.4)' in l:
        new_lines.append(l.replace('report.md (Section 1, 4.4)', '"report.md (Section 1, 4.4)"'))
    else:
        new_lines.append(l)

df = pd.read_csv(io.StringIO(''.join(new_lines)))
print('Parsed successfully with fix!')
print('Shape:', df.shape)
print('Columns:', list(df.columns))
print(df['evidence_status'].value_counts())
print('\nAll claims in ledger:')
for idx, row in df.iterrows():
    print(f"{row['claim_id']}: {row['metric_name']} | val={row['reported_value'][:30]} | status={row['evidence_status']}")
