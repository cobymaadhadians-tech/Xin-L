"""Projection-null calibration using retained 200-block sampling covariance.
Run from the extracted analysis archive. Does not change inference thresholds.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import chi2,norm
p=argparse.ArgumentParser();p.add_argument('--draws',type=int,default=100000);p.add_argument('--seed',type=int,default=20260911);args=p.parse_args()
root=Path.cwd();out=root/'results/anchored_null_calibration';out.mkdir(parents=True,exist_ok=True)
z=np.load(root/'results/common_snp_revision/analysis_six/snp_N.npz');tr=z['traits'].tolist();S=z['S'];blocks=z['blocks'];n=len(tr)
ix=[(i,j) for j in range(n) for i in range(j,n)];pos={x:k for k,x in enumerate(ix)}
def idx(a,b):return pos[max(a,b),min(a,b)]
a=tr.index('PGC');traits=['BD','MDD','ADHD','OCD','PTSD','AUD'];sources=['FG','MVP'];qix=[(tr.index(t),tr.index(s)) for t in traits for s in sources]
theta=np.array([S[i,j] for i,j in ix]);bb=np.array([blocks[:,i,j] for i,j in ix]).T;cen=bb-bb.mean(0);V=(len(bb)-1)/len(bb)*cen.T@cen
mu=theta.copy()
for j,t in qix:mu[idx(t,j)]=S[t,a]*S[a,j]/S[a,a]
def evaluate(x):
 x=np.atleast_2d(x);m=len(x);d=np.empty((m,12));J=np.zeros((m,12,len(ix)));v=x[:,idx(a,a)]
 for k,(j,t) in enumerate(qix):
  ta=x[:,idx(t,a)];aj=x[:,idx(a,j)];d[:,k]=x[:,idx(t,j)]-ta*aj/v
  J[:,k,idx(t,j)]+=1;J[:,k,idx(t,a)]-=aj/v;J[:,k,idx(a,j)]-=ta/v;J[:,k,idx(a,a)]+=ta*aj/v**2
 C=(J@V)@J.transpose(0,2,1)
 Q=np.einsum('bi,bi->b',d,np.linalg.solve(C,d[...,None])[...,0]);qs=[Q];dfs=[12];labels=['global']
 for k,t in enumerate(traits):
  ds=d[:,2*k:2*k+2];cs=C[:,2*k:2*k+2,2*k:2*k+2];qs.append(np.einsum('bi,bi->b',ds,np.linalg.solve(cs,ds[...,None])[...,0]));dfs.append(2);labels.append(t)
 for k,(t,s) in enumerate((t,s) for t in traits for s in sources):qs.append(d[:,k]**2/C[:,k,k]);dfs.append(1);labels.append(t+'_'+s)
 return np.array(qs).T,np.array(dfs),labels
observed,dfs,labels=evaluate(theta);ev,U=np.linalg.eigh(V);assert ev.min()>-ev.max()*1e-10
L=U*np.sqrt(np.maximum(ev,0));rng=np.random.default_rng(args.seed);qq=[];nonpositive=0
for start in range(0,args.draws,1000):
 x=mu+rng.standard_normal((min(1000,args.draws-start),len(mu)))@L.T;nonpositive+=int((x[:,idx(a,a)]<=0).sum());qq.append(evaluate(x)[0])
Q=np.concatenate(qq);P=chi2.sf(Q,dfs);rows=[]
for j,label in enumerate(labels):
 for alpha in [.05,.01,.001]:
  count=int((P[:,j]<alpha).sum());rate=count/args.draws;se=(rate*(1-rate)/args.draws)**.5
  rows.append(dict(test=label,df=int(dfs[j]),alpha=alpha,rejections=count,draws=args.draws,type_I_error=rate,MC_SE=se,MC_lower=max(0,rate-1.96*se),MC_upper=min(1,rate+1.96*se)))
pd.DataFrame(rows).to_csv(out/'type_I_error.tsv',sep='\t',index=False)
rows=[]
for j,label in enumerate(labels):
 for probability in [.5,.9,.95,.99,.999]:rows.append(dict(test=label,df=int(dfs[j]),probability=probability,empirical_Q_quantile=np.quantile(Q[:,j],probability),chi_square_Q_quantile=chi2.ppf(probability,dfs[j])))
pd.DataFrame(rows).to_csv(out/'quantiles.tsv',sep='\t',index=False)
pd.DataFrame(dict(test=labels,df=dfs,observed_Q=observed[0],asymptotic_P=chi2.sf(observed[0],dfs),null_exceedances=(Q>=observed).sum(0),MC_tail_P=(1+(Q>=observed).sum(0))/(args.draws+1))).to_csv(out/'observed_tests.tsv',sep='\t',index=False)
meta=dict(draws=args.draws,seed=args.seed,blocks=len(bb),nonpositive_anchor_variance_draws=nonpositive,min_sampling_covariance_eigenvalue=float(ev.min()),null_max_deviation=float(np.max(abs(evaluate(mu)[0]))),scope='Conditional Gaussian plug-in calibration; empirical sampling covariance fixed, projection and delta-method variance recomputed per draw. No positive-definite filtering, threshold changes or covariance re-estimation from synthetic blocks. Does not calibrate estimation error in V or extreme observed tail probabilities.')
(out/'settings.json').write_text(json.dumps(meta,indent=2)+'\n');print(pd.DataFrame(rows).query("test=='global'").to_string(index=False));print(pd.read_csv(out/'type_I_error.tsv',sep='\t').query('alpha==0.05').to_string(index=False))
