"""Paired delete-block R2 bias diagnostics and source correlations.

The bias-adjusted values are point-estimate sensitivities. Their uncertainty is
not supplied: the original first-order standard errors do not estimate the
sampling variance of the bias-adjusted estimator.
"""
from pathlib import Path
import os,itertools
import numpy as np
import pandas as pd

ROOT=Path(os.environ.get('ANALYSIS_ARCHIVE_ROOT',Path.cwd()))
OUT=ROOT/'.reproduction/bias_diagnostics'
OUT.mkdir(parents=True,exist_ok=True)
z=np.load(ROOT/'results/common_snp_revision/analysis_six/snp_N.npz')
S=z['S'];blocks=z['blocks'];traits=z['traits'].tolist();B=len(blocks)
assert B==200

def jk_se(values):
    return np.sqrt((B-1)/B*np.sum((values-values.mean())**2))

def r2(s,t,a):
    return float(s[t,a]@np.linalg.solve(s[np.ix_(a,a)],s[a,t])/s[t,t])

totals=[];contrasts=[];correlations=[]
for system in ['six','five']:
    aux=[traits.index(t) for t in ['BD','MDD','ADHD','OCD','PTSD']+(['AUD'] if system=='six' else [])]
    values={};deleted={};bias={}
    for source in ['PGC','FG','MVP','Broad']:
        t=traits.index(source);full=r2(S,t,aux);dv=np.array([r2(s,t,aux) for s in blocks]);b=(B-1)*(dv.mean()-full);se=jk_se(dv)
        values[source]=full;deleted[source]=dv;bias[source]=b
        totals.append(dict(system=system,source=source,n_blocks=B,estimate=full,SE=se,lower=full-1.96*se,upper=full+1.96*se,jackknife_bias=b,bias_adjusted_estimate=full-b))
    for family,pairs in [('sources',[('PGC','FG'),('MVP','PGC'),('MVP','FG')]),('endpoint',[('PGC','FG'),('Broad','PGC'),('Broad','FG')])]:
        for left,right in pairs:
            d=values[left]-values[right];db=deleted[left]-deleted[right];b=(B-1)*(db.mean()-d);se=jk_se(db)
            assert abs(b-(bias[left]-bias[right]))<1e-12
            contrasts.append(dict(system=system,family=family,comparison=left+'_minus_'+right,n_blocks=B,estimate=d,SE=se,lower=d-1.96*se,upper=d+1.96*se,jackknife_bias=b,bias_adjusted_estimate=d-b))
for left,right in itertools.combinations(['PGC','FG','MVP','Broad'],2):
    i,j=traits.index(left),traits.index(right)
    full=S[i,j]/np.sqrt(S[i,i]*S[j,j]);dv=blocks[:,i,j]/np.sqrt(blocks[:,i,i]*blocks[:,j,j]);se=jk_se(dv)
    correlations.append(dict(left=left,right=right,n_snps=int(z['n_snps']),n_blocks=B,rg=full,SE=se,lower=full-1.96*se,upper=full+1.96*se))
for name,data in [('r2_bias_totals',totals),('r2_bias_contrasts',contrasts),('source_correlations',correlations)]:
    pd.DataFrame(data).to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
    print(name);print(pd.DataFrame(data).to_string(index=False))
