import pandas as pd 
import numpy as np

df=pd.read_csv('sample_table_airway.csv')

fq1=df.X+"_1.fastq"
fq2=df.X+"_2.fastq"
ndf={"sample":df.X,"unit":np.repeat("rep1",len(df.X)),"fq1":fq1,"fq2":fq2}

samples={"sample":df.X,"cell":df.cell,"dex":df.dex}

ndf=pd.DataFrame.from_dict(ndf)
samples=pd.DataFrame.from_dict(samples)

ndf.to_csv("units.tsv",sep='\t',index=False)
samples.to_csv("samples.tsv",sep='\t',index=False)
