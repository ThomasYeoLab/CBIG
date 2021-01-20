#!/usr/bin/env python
import sys
import pandas as pd


in_columns = ['RID', 'D1', 'D2', 'EXAMDATE']
out_columns = ['RID', 'EXAMDATE', 'train', 'val', 'test']

df = pd.read_csv(sys.argv[1], usecols=in_columns)
df['train'] = df['D1']
df['val'] = df['D2']
df['test'] = df['D2']
df[out_columns].to_csv(sys.argv[2], index=False)
