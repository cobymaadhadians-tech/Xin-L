"""Conditional equal-profile calibration using actual matched delete-block perturbations.
Unlike a fixed covariance Wald replay, recomputes all subsets and the exact
200-block jackknife covariance for every simulated estimate. Holds the empirical
correlation-error covariance and block perturbations fixed; no PD filtering of
outer simulations. Run from project root.
"""
from pathlib import Path
import argparse,itertools,math,time,json
import numpy as np,pandas as pd
from scipy.stats import chi2
P=argparse.ArgumentParser();P.add_argument('--draws',type=int,default=10000);P.add_argument('--batch',type=int,default=10);P.add_argument('--null',choices=['global','PGC_minus_FG','MVP_minus_PGC','MVP_minus_FG'],default='global');args=P.parse_args()
AR=Path.cwd() if Path('results/common_snp_revision').exists() else Path('analysis_repository')
O=AR/'results/sharing_null_calibration';O.mkdir(parents=True,exist_ok=True)
z=np.load(AR/'results/common_snp_revision/analysis_six/snp_N.npz');tr=z['traits'].tolist();idx=[tr.index(t) for t in ['BD','MDD','ADHD','OCD','PTSD','AUD','PGC','FG','MVP']];S=z['S'][np.ix_(idx,idx)];bs=z['blocks'][:,idx,:][:,:,idx];R=S/np.sqrt(np.outer(S.diagonal(),S.diagonal()));br=bs/np.sqrt(np.einsum('bi,bj->bij',bs.diagonal(axis1=1,axis2=2),bs.diagonal(axis1=1,axis2=2)));B=len(br)
ii,jj=np.tril_indices(9,-1);x=R[ii,jj];bx=br[:,ii,jj];bc=bx-bx.mean(0);V=(B-1)/B*bc.T@bc;ev,u=np.linalg.eigh(V);L=u*np.sqrt(np.maximum(ev,0));rng=np.random.default_rng(20260912)
W=np.zeros((64,6));subsets=[]
for mask in range(1,64):subsets.append((mask,[j for j in range(6) if mask&(1<<j)]))
for j in range(6):
 for mask in range(64):
  if mask&(1<<j):continue
  k=mask.bit_count();w=math.factorial(k)*math.factorial(5-k)/math.factorial(6);W[mask|(1<<j),j]+=w;W[mask,j]-=w

def matrix(x):
 m=np.broadcast_to(np.eye(9),(len(x),9,9)).copy();m[:,ii,jj]=x;m[:,jj,ii]=x;return m

def decomp(m):
 vals=np.zeros((len(m),3,64))
 for mask,a in subsets:
  A=m[:,a,:][:,:,a];v=m[:,6:9,:][:,:,a];b=np.linalg.solve(A,v.transpose(0,2,1));vals[:,:,mask]=np.einsum('bsi,bis->bs',v,b)
 raw=vals@W;return raw/raw.sum(-1,keepdims=True),vals[:,:,-1]

def tests(m):
 # Full estimate plus centered matched block perturbations, exact same nonlinear jackknife estimator.
 xx=m[:,ii,jj];allx=np.concatenate([xx[:,None,:],xx[:,None,:]+bc[None,:,:]],axis=1);mm=matrix(allx.reshape(-1,len(ii)));sh,tot=decomp(mm);sh=sh.reshape(len(m),B+1,3,6)
 qs=[];labels=[]
 for label,src in [('global',[0,1,2]),('PGC_minus_FG',[0,1]),('MVP_minus_PGC',[2,0]),('MVP_minus_FG',[2,1])]:
  if len(src)==3:d=np.concatenate([sh[:,:,1,:5]-sh[:,:,0,:5],sh[:,:,2,:5]-sh[:,:,0,:5]],axis=-1)
  else:d=sh[:,:,src[0],:5]-sh[:,:,src[1],:5]
  dd=d[:,0];db=d[:,1:];db=db-db.mean(1,keepdims=True);C=(B-1)/B*db.transpose(0,2,1)@db;sol=np.linalg.solve(C,dd[...,None])[...,0];qs.append(np.einsum('bi,bi->b',dd,sol));labels.append(label)
 return np.array(qs).T,labels
# GLS fit under proportional target-correlation vectors. Locally this imposes
# the same number of constraints as equal normalized profiles (10 or 5).
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
A=R[:6,:6];v=R[6:9,:6];tot=np.einsum('si,is->s',v,np.linalg.solve(A,v.T));ai,aj=np.tril_indices(6,-1)
use=np.array([k for k,(i,j) in enumerate(zip(ii,jj)) if not(i>=6 and j>=6)])
cv=V[np.ix_(use,use)];chol=np.linalg.cholesky(cv)
group={'global':[0,1,2],'PGC_minus_FG':[0,1],'MVP_minus_PGC':[2,0],'MVP_minus_FG':[2,1]}[args.null]
others=[i for i in range(3) if i not in group]
direction=(v[group]/np.sqrt(tot[group,None])).mean(0);direction/=np.sqrt(direction@np.linalg.solve(A,direction));base=direction*np.sqrt(tot[group[0]])
p0=np.r_[A[ai,aj],base,np.log(np.sqrt(tot[group[1:]]/tot[group[0]])),v[others].ravel()]
def model(p):
 m=R.copy();m[ai,aj]=p[:15];m[aj,ai]=p[:15];m[6+group[0],:6]=p[15:21];k=21
 for j in group[1:]:m[6+j,:6]=p[15:21]*np.exp(p[k]);k+=1
 for j in others:m[6+j,:6]=p[k:k+6];k+=6
 m[:6,6:]=m[6:,:6].T;return m
def fun(p):return solve_triangular(chol,(model(p)[ii,jj]-x)[use],lower=True)
fit=least_squares(fun,p0,xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=3000);assert fit.success
mu=model(fit.x);assert np.linalg.eigvalsh(mu[:6,:6])[0]>0
assert np.max(np.ptp(decomp(mu[None])[0][0,group],axis=0))<1e-10
print('Null',args.null,'GLS objective',float(fit.fun@fit.fun),'aux eigenvalue',np.linalg.eigvalsh(mu[:6,:6])[0],flush=True)
obs,labels=tests(R[None]);print('Observed',obs,flush=True)
qall=[];bad=0;start=time.time()
for k in range(0,args.draws,args.batch):
 count=min(args.batch,args.draws-k);xx=mu[ii,jj]+rng.normal(size=(count,len(ii)))@L.T;m=matrix(xx);bad+=int((np.linalg.eigvalsh(m[:,:6,:6])[:,0]<=0).sum());q,_=tests(m);qall.append(q)
 if k%200==0:print(k+count,'seconds',time.time()-start,flush=True)
Q=np.concatenate(qall);dfs=np.array([10,5,5,5]);rows=[]
for j,label in enumerate(labels):
 for alpha in [.05,.01,.001]:
  n=int((Q[:,j]>chi2.isf(alpha,dfs[j])).sum());rate=n/args.draws;se=np.sqrt(rate*(1-rate)/args.draws);rows.append(dict(test=label,df=dfs[j],alpha=alpha,rejections=n,draws=args.draws,type_I_error=rate,MC_SE=se,MC_lower=max(0,rate-1.96*se),MC_upper=min(1,rate+1.96*se)))
pd.DataFrame(rows).to_csv(O/(args.null+'_normalized_profile_type_I_error.tsv'),sep='\t',index=False)
pd.DataFrame(dict(test=labels,df=dfs,Q=obs[0],P=chi2.sf(obs[0],dfs))).to_csv(O/'observed_replay.tsv',sep='\t',index=False)
np.savez(O/(args.null+'_null_statistics.npz'),Q=Q,mu=mu,V=V)
(O/(args.null+'_settings.json')).write_text(json.dumps(dict(null=args.null,GLS_objective=float(fit.fun@fit.fun),draws=args.draws,seed=20260912,blocks=B,nonPD_auxiliary_draws=bad,scope='Conditional Gaussian errors with fixed empirical correlation sampling covariance; GLS-fitted proportional correlation vectors yield equal normalized profiles for the tested sources; all 64 subset estimates and 200-block jackknife covariance re-estimated per draw using recentered observed block perturbations; no outer draw filtering; does not simulate uncertainty in original block perturbations.'),indent=2))
print(pd.DataFrame(rows).to_string(index=False),flush=True)

j=labels.index(args.null);tail=int((Q[:,j]>=obs[0,j]).sum());pd.DataFrame([dict(test=args.null,Q=obs[0,j],df=int(dfs[j]),asymptotic_P=chi2.sf(obs[0,j],dfs[j]),exceedances=tail,draws=args.draws,MC_P=(tail+1)/(args.draws+1))]).to_csv(O/(args.null+'_bootstrap_P.tsv'),sep='\t',index=False)
