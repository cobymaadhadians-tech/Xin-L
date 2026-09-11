"""PHBC-style point correction with legacy covariance diagnostics, paper eqs 28-29, on matched SNP-specific-N inputs.
Run from project or extracted archive root. For current total inference use
test_total_sharing_gaussian.py; the noisy-projection covariance below is superseded.
Paper auxiliary-PD rule is primary; software full-PD selection is sensitivity.
Uses analytical quadratic root, identical to converged binary search for fixed draws.
"""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import chi2,norm
ROOT=Path.cwd();AR=ROOT if (ROOT/'results/common_snp_revision').exists() else ROOT/'analysis_repository';OUT=AR/'results/six_sharing_correction';OUT.mkdir(parents=True,exist_ok=True)
def covjk(x):
 x=x-x.mean(0);return (len(x)-1)/len(x)*x.T@x
def bh(p):
 p=np.array(p);i=np.argsort(p);o=np.empty(len(p));o[i]=np.minimum(1,np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return o
def run(folder,mode,system,draws=100000):
 z=np.load(AR/f'results/common_snp_revision/{folder}/{mode}.npz');tr=z['traits'].tolist();s=z['S'];bs=z['blocks'];a=[tr.index(t) for t in ['BD','MDD','ADHD','OCD','PTSD']+(['AUD'] if system=='six' else [])];sources=['PGC','FG','MVP','Broad'];ts=[tr.index(t) for t in sources]
 r=s/np.sqrt(np.outer(s.diagonal(),s.diagonal()));br=bs/np.sqrt(np.einsum('bi,bj->bij',bs.diagonal(axis1=1,axis2=2),bs.diagonal(axis1=1,axis2=2)))
 ii,jj=np.tril_indices(len(tr),-1);V=covjk(br[:,ii,jj]);ev,u=np.linalg.eigh(V);assert ev.min()>-1e-10
 rng=np.random.default_rng(20260911);e=rng.normal(size=(draws,len(ii)))@(u*np.sqrt(np.maximum(ev,0))).T;noise=np.zeros((draws,len(tr),len(tr)));noise[:,ii,jj]=e;noise[:,jj,ii]=e
 A=r[np.ix_(a,a)];As=A+noise[:,a,:][:,:,a];ok=np.linalg.eigvalsh(As)[:,0]>0;ai=np.linalg.inv(As[ok]);q=[];weights=[];scales=[];before=[];after=[];jk=[];rows=[]
 for t,name in zip(ts,sources):
  v=r[t,a];initial=v@np.linalg.solve(A,v);en=noise[ok,t,:][:,a];aa=np.einsum('i,bij,j->b',v,ai,v);bb=2*np.einsum('i,bij,bj->b',v,ai,en);cc=np.einsum('bi,bij,bj->b',en,ai,en)
  w=(-bb.mean()+np.sqrt(bb.mean()**2-4*aa.mean()*(cc.mean()-initial)))/(2*aa.mean());assert .5<=w<=1
  un=aa+bb+cc;co=w*w*aa+w*bb+cc;sc=co.std(ddof=1)/un.std(ddof=1)
  bv=br[:,t,:][:,a];Ab=br[:,a,:][:,:,a];bj=np.einsum('bi,bi->b',bv,np.linalg.solve(Ab,bv[...,None])[...,0]);se=np.sqrt(covjk(bj[:,None])[0,0]);est=initial*w*w
  fullok=np.linalg.eigvalsh(r[np.ix_([t]+a,[t]+a)]+noise[:,[t]+a,:][:,:,[t]+a])[:,0]>0;mask=fullok[ok];aa2=aa[mask].mean();bb2=bb[mask].mean();cc2=cc[mask].mean();w2=(-bb2+np.sqrt(bb2**2-4*aa2*(cc2-initial)))/(2*aa2)
  halves=[]
  for sel in [slice(0,len(aa)//2),slice(len(aa)//2,None)]:
   ap,bp,cp=aa[sel].mean(),bb[sel].mean(),cc[sel].mean();wp=(-bp+np.sqrt(bp*bp-4*ap*(cp-initial)))/(2*ap);halves.append(initial*wp*wp)
  rows.append(dict(analysis=folder+'_'+mode,system=system,source=name,n_snps=int(z['n_snps']),draws=draws,accepted=int(ok.sum()),initial=initial,weight=w,estimate=est,SE=se*sc,lower=est-1.96*se*sc,upper=est+1.96*se*sc,SE_scale=sc,software_filter_estimate=initial*w2*w2,software_filter_acceptance=fullok.mean(),MC_half_difference=halves[0]-halves[1]))
  q.append(est);weights.append(w);scales.append(sc);before.append(un);after.append(co);jk.append(bj)
 q=np.array(q);C0=covjk(np.array(jk).T);sd0=np.array(before).std(axis=1,ddof=1);Ca=np.cov(np.array(after));# joint extension preserving corrected MC cross-source correlation and official marginal SE
 se=np.sqrt(np.diag(C0))*np.array(scales);C=Ca/np.sqrt(np.outer(Ca.diagonal(),Ca.diagonal()))*np.outer(se,se)
 tests=[]
 for fam,ids in [('sources',[0,1,2]),('endpoint',[0,1,3])]:
  D=np.zeros((2,4));D[0,ids[1]]=1;D[0,ids[0]]=-1;D[1,ids[2]]=1;D[1,ids[0]]=-1
  d=D@q;cv=D@C@D.T;Q=d@np.linalg.solve(cv,d);tests.append(dict(analysis=folder+'_'+mode,system=system,family=fam,comparison='global',Q=Q,df=2,P=chi2.sf(Q,2)))
  for l,rr in [(ids[0],ids[1]),(ids[2],ids[0]),(ids[2],ids[1])]:
   delta=q[l]-q[rr];ss=np.sqrt(C[l,l]+C[rr,rr]-2*C[l,rr]);tests.append(dict(analysis=folder+'_'+mode,system=system,family=fam,comparison=sources[l]+'_minus_'+sources[rr],estimate=delta,SE=ss,lower=delta-1.96*ss,upper=delta+1.96*ss,Q=(delta/ss)**2,df=1,P=2*norm.sf(abs(delta/ss))))
 T=pd.DataFrame(tests);T['q']=np.nan
 for fam in ['sources','endpoint']:
  ix=(T.family==fam)&(T.comparison!='global');T.loc[ix,'q']=bh(T.loc[ix,'P'])
 if folder=='analysis_six' and mode=='snp_N' and system=='six':np.savez(OUT/'principal.npz',estimate=q,covariance=C,weights=weights,SE_scales=scales,uncorrected_covariance=C0,sources=sources)
 return rows,T
if __name__=='__main__':
 rows=[];tests=[]
 for folder,mode,sys in [('analysis_six','snp_N','six'),('analysis_six','mean_N','six'),('analysis_six','snp_N','five'),('analysis','snp_N','five')]:
  r,t=run(folder,mode,sys);rows+=r;tests.append(t);print(pd.DataFrame(r).to_string(index=False),flush=True);print(t.to_string(index=False),flush=True)
 pd.DataFrame(rows).to_csv(OUT/'totals.tsv',sep='\t',index=False);pd.concat(tests).to_csv(OUT/'tests.tsv',sep='\t',index=False)
 (OUT/'method.json').write_text(json.dumps(dict(seed=20260911,draws=100000,primary_rule='Auxiliary correlation matrix positive definite, Nature Genetics eq28',point_solver='exact positive quadratic root in [0,1]',SE='Marginal eq29 MC SD-ratio times uncorrected matched-block jackknife SE; cross-source correlation from joint corrected MC draws',uncertainty_status='Superseded covariance extension; use total_sharing_inference tests and population profile intervals',scope='PHBC-style adaptation to matched SNP-specific N; paired covariance is joint-MC extension; sampling covariance fixed; no claim of official package output'),indent=2))
