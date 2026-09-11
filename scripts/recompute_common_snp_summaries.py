from pathlib import Path
import json,itertools,math
import numpy as np,pandas as pd
from scipy.stats import chi2,norm
import os
R=Path(os.environ.get('ANALYSIS_ARCHIVE_ROOT',Path.cwd()));W=R/'results/common_snp_revision';O=R/'.reproduction/summaries';O.mkdir(parents=True,exist_ok=True)
B=200
names={'BD_PGC2021':'BD','MDD_CLIN_PGC2025':'MDD','OCD_2025':'OCD','PGC3':'PGC','FinnGen_R13':'FG','MVP_AOU':'MVP','PGC':'PGC','narrow':'FG','broad':'Broad'}
def bh(p):
 p=np.array(p);ix=np.argsort(p);v=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];o=np.empty(len(p));o[ix]=np.minimum(v,1);return o
def covjk(b):return (B-1)/B*((b-b.mean(0)).T@(b-b.mean(0)))
def wald(d,V):
 e,U=np.linalg.eigh((V+V.T)/2);assert e.min()>-max(e.max(),1)*1e-9
 ix=e>e.max()*1e-9;q=float(np.sum((U[:,ix].T@d)**2/e[ix]));return q,int(ix.sum()),float(chi2.sf(q,ix.sum()))
def shap(R,t,aux):
 r=R[np.ix_(aux,aux)];v=R[aux,t];p=len(aux);r2={0:0.}
 for mask in range(1,2**p):
  ix=[i for i in range(p) if mask&(1<<i)];r2[mask]=v[ix]@np.linalg.solve(r[np.ix_(ix,ix)],v[ix])
 phi=np.zeros(p)
 for j in range(p):
  for mask in range(2**p):
   if mask&(1<<j):continue
   s=mask.bit_count();phi[j]+=math.factorial(s)*math.factorial(p-s-1)/math.factorial(p)*(r2[mask|(1<<j)]-r2[mask])
 assert abs(sum(phi)-r2[2**p-1])<1e-10
 return phi
est=[];totals=[];tests=[];contrasts=[];bridge_est=[];bridge_tests=[];qa=[]
def decomp(label,sources,sets=('six','five')):
 for system in sets:
  pred=['BD','MDD','ADHD','OCD','PTSD']+(['AUD'] if system=='six' else []);p=len(pred);phi={};blocks={}
  for k,s in sources.items():
   tr=s['traits'];aux=[tr.index(a) for a in pred];t=tr.index(s['target']);F=np.asarray(s['full']);bb=np.asarray(s['blocks']);phi[k]=shap(F,t,aux);blocks[k]=np.array([shap(m,t,aux) for m in bb])
   total=phi[k].sum();bt=blocks[k].sum(1);se=np.sqrt(covjk(bt[:,None])[0,0]);totals.append(dict(analysis=label,system=system,source=k,estimate=total,SE=se,lower=total-1.96*se,upper=total+1.96*se))
   for metric in ['raw','share']:
    a=phi[k] if metric=='raw' else phi[k]/total;b=blocks[k] if metric=='raw' else blocks[k]/bt[:,None];se=np.sqrt(np.diag(covjk(b)))
    for j,x in enumerate(pred):est.append(dict(analysis=label,system=system,source=k,predictor=x,metric=metric,estimate=a[j],SE=se[j],lower=a[j]-1.96*se[j],upper=a[j]+1.96*se[j]))
  families=[('sources',['PGC','FG','MVP']),('endpoint',['PGC','FG','Broad'])]
  for fam,ks in families:
   if not set(ks)<=sources.keys():continue
   for metric in ['total','raw','share']:
    aa={k:np.array([phi[k].sum()]) if metric=='total' else phi[k] if metric=='raw' else phi[k]/phi[k].sum() for k in ks}
    bb={k:blocks[k].sum(1)[:,None] if metric=='total' else blocks[k] if metric=='raw' else blocks[k]/blocks[k].sum(1)[:,None] for k in ks}
    delta=np.concatenate([aa[ks[1]]-aa[ks[0]],aa[ks[2]]-aa[ks[0]]]);db=np.concatenate([bb[ks[1]]-bb[ks[0]],bb[ks[2]]-bb[ks[0]]],axis=1)
    q,df,pv=wald(delta,covjk(db));tests.append(dict(analysis=label,system=system,family=fam,metric=metric,comparison='global',Q=q,df=df,P=pv))
    for left,right in [(ks[0],ks[1]),(ks[2],ks[0]),(ks[2],ks[1])]:
     d=aa[left]-aa[right];db=bb[left]-bb[right];V=covjk(db);q,df,pv=wald(d,V);comparison=left+'_minus_'+right
     tests.append(dict(analysis=label,system=system,family=fam,metric=metric,comparison=comparison,Q=q,df=df,P=pv))
     se=np.sqrt(np.diag(V))
     for j,x in enumerate(['total'] if metric=='total' else pred):contrasts.append(dict(analysis=label,system=system,family=fam,metric=metric,comparison=comparison,predictor=x,estimate=d[j],SE=se[j],lower=d[j]-1.96*se[j],upper=d[j]+1.96*se[j],P=2*norm.sf(abs(d[j]/se[j]))))
def bridge(label,S,V,tr,sets=('six','five')):
 # R's lower-triangle column-major ordering.
 ix=[(i,j) for j in range(len(tr)) for i in range(j,len(tr))];theta=np.array([S[i,j] for i,j in ix]);pos={a:k for k,a in enumerate(ix)}
 def idx(a,b):return pos[max(a,b),min(a,b)]
 for system in sets:
  pred=['BD','MDD','ADHD','OCD','PTSD']+(['AUD'] if system=='six' else [])
  for anchor in ['PGC','FG','MVP']:
   a=tr.index(anchor);others=[x for x in ['PGC','FG','MVP'] if x!=anchor];d=[];J=[]
   for trait in pred:
    j=tr.index(trait)
    for target in others:
     t=tr.index(target);v=S[t,j]-S[t,a]*S[a,j]/S[a,a];g=np.zeros(len(ix));g[idx(t,j)]+=1;g[idx(t,a)]-=S[a,j]/S[a,a];g[idx(a,j)]-=S[t,a]/S[a,a];g[idx(a,a)]+=S[t,a]*S[a,j]/S[a,a]**2;d.append(v);J.append(g)
   d=np.array(d);J=np.array(J);C=J@V@J.T;se=np.sqrt(np.diag(C));pvals=2*norm.sf(abs(d/se));qs=bh(pvals)
   for k in range(len(d)):bridge_est.append(dict(analysis=label,system=system,anchor=anchor,trait=pred[k//2],source=others[k%2],estimate=d[k],SE=se[k],lower=d[k]-1.96*se[k],upper=d[k]+1.96*se[k],P=pvals[k],q=qs[k]))
   for k,trait in enumerate(['ALL']+pred):
    ids=np.arange(len(d)) if k==0 else np.arange((k-1)*2,k*2);q,df,pv=wald(d[ids],C[np.ix_(ids,ids)]);bridge_tests.append(dict(analysis=label,system=system,anchor=anchor,trait=trait,Q=q,df=df,P=pv,P_Bonf=min(1,pv*len(pred)) if k else np.nan))
for file,label in [('original_correlations.json','original'),('endpoint_correlations.json','original_endpoint'),('original_common_correlations.json','original_common')]:
 p=W/'analysis'/file
 if not p.exists():continue
 src=json.loads(p.read_text());new={}
 for k,s in src.items():
  s['traits']=[names.get(a,a) for a in s['traits']];new[names.get(k,k)]=s
 decomp(label,new)
b=json.loads((W/'analysis/original_bridge.json').read_text());bridge('original',np.array(b['S']),np.array(b['V']),b['traits'])
for folder,sets in [('analysis',('five',)),('analysis_six',('six','five'))]:
 for mode in ['snp_N','mean_N']:
  f=W/folder/(mode+'.npz')
  if not f.exists():continue
  z=np.load(f);S=z['S'];bl=z['blocks'];tr=z['traits'].tolist();label=folder+'_'+mode
  rs=S/np.sqrt(np.outer(S.diagonal(),S.diagonal()));rbs=bl/np.sqrt(np.einsum('bi,bj->bij',np.diagonal(bl,axis1=1,axis2=2),np.diagonal(bl,axis1=1,axis2=2)))
  decomp(label,{k:dict(traits=tr,target=k,full=rs,blocks=rbs) for k in ['PGC','FG','MVP','Broad']},sets)
  ix=[(i,j) for j in range(len(tr)) for i in range(j,len(tr))];V=covjk(np.array([bl[:,i,j] for i,j in ix]).T);bridge(label,S,V,tr,sets)
  qa.append(dict(analysis=label,n_snps=int(z['n_snps']),min_full_eigenvalue=float(np.linalg.eigvalsh(rs).min()),min_block_eigenvalue=float(min(np.linalg.eigvalsh(m).min() for m in rbs))))
T=pd.DataFrame(tests);T['q_pairwise']=np.nan;T['q_global']=np.nan
for keys,g in T[(T.comparison=='global') & T.metric.isin(['raw','share'])].groupby(['analysis','system','family']):T.loc[g.index,'q_global']=bh(g.P)
for keys,g in T[T.comparison!='global'].groupby(['analysis','system','family','metric']):T.loc[g.index,'q_pairwise']=bh(g.P)
C=pd.DataFrame(contrasts);C['q']=np.nan
for keys,g in C.groupby(['analysis','system','family','metric']):C.loc[g.index,'q']=bh(g.P)
for n,x in [('estimates',est),('totals',totals),('profile_tests',T),('contrasts',C),('bridge_estimates',bridge_est),('bridge_tests',bridge_tests),('matrix_qc',qa)]:pd.DataFrame(x).to_csv(O/(n+'.tsv'),sep='\t',index=False)
print(T[(T.comparison=='global') & (T.family=='sources')][['analysis','system','metric','Q','df','P']].to_string(index=False))
print(pd.DataFrame(bridge_tests).query("anchor=='PGC' and trait=='ALL'").to_string(index=False))
