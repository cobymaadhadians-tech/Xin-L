"""Conditional Gaussian minimum-distance tests of equal population joint sharing.
Primary inference uses the original matched-block correlation covariance.
PHBC-corrected point estimates are computed separately. The confidence intervals
profile the population parameter and need not centre on the corrected estimator.
All simulated observations are retained; numerical failures are bounded explicitly.
"""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize
from scipy.linalg import solve_triangular
from scipy.stats import chi2
import argparse,time,json
P=argparse.ArgumentParser();P.add_argument('--null',default='global');P.add_argument('--draws',type=int,default=10000);P.add_argument('--intervals',action='store_true');args=P.parse_args()
root=Path.cwd();AR=root if (root/'results/common_snp_revision').exists() else root/'analysis_repository';out=AR/'results/total_sharing_inference';out.mkdir(exist_ok=True);z=np.load(AR/'results/common_snp_revision/analysis_six/snp_N.npz');tr=z['traits'].tolist();names=['BD','MDD','ADHD','OCD','PTSD','AUD','PGC','FG',('Broad' if args.null=='Broad_minus_FG' else 'MVP')];idx=[tr.index(t) for t in names];s=z['S'][np.ix_(idx,idx)];bs=z['blocks'][:,idx,:][:,:,idx];R=s/np.sqrt(np.outer(s.diagonal(),s.diagonal()));br=bs/np.sqrt(np.einsum('bi,bj->bij',bs.diagonal(axis1=1,axis2=2),bs.diagonal(axis1=1,axis2=2)));ii,jj=np.tril_indices(9,-1);take=~((ii>=6)&(jj>=6));ii=ii[take];jj=jj[take];x=R[ii,jj];bc=br[:,ii,jj]-br[:,ii,jj].mean(0);V=199/200*bc.T@bc;L=np.linalg.cholesky(V);group={'global':[0,1,2],'PGC_minus_FG':[0,1],'MVP_minus_PGC':[2,0],'MVP_minus_FG':[2,1],'Broad_minus_FG':[2,1]}[args.null];D=np.zeros((len(group)-1,3))
for j,k in enumerate(group[1:]):D[j,k]=1;D[j,group[0]]=-1
floor=1e-9
delta=0.

def decode(v):
 m=np.eye(9);m[ii,jj]=v;m[jj,ii]=v;return m[:6,:6],m[6:,:6]
def quantities(v):
 A,r=decode(v);beta=np.linalg.solve(A,r.T).T;tot=np.einsum('si,si->s',r,beta);grad=np.zeros((3,len(ii)))
 for k,(i,j) in enumerate(zip(ii,jj)):
  if i<6:grad[:,k]=-2*beta[:,i]*beta[:,j]
  else:grad[i-6,k]=2*beta[i-6,j]
 return tot,grad

def fit(y,start):
 t0=solve_triangular(L,start-y,lower=True)
 def cons(t):
  v=y+L@t;q,g=quantities(v);return D@q-delta
 def jac(t):return D@quantities(y+L@t)[1]@L
 def feasible(t):
  v=y+L@t;A,r=decode(v);q,g=quantities(v);return np.r_[np.linalg.eigvalsh(A)[0]-floor,1-q[group]]
 def fj(t):
  v=y+L@t;A,r=decode(v);w,u=np.linalg.eigh(A);q,g=quantities(v);gg=np.zeros(len(ii));a=ii<6;gg[a]=2*u[ii[a],0]*u[jj[a],0];return np.vstack([gg,-g[group]])@L
 res=minimize(lambda t:float(t@t),t0,jac=lambda t:2*t,constraints=[{'type':'eq','fun':cons,'jac':jac},{'type':'ineq','fun':feasible,'jac':fj}],method='SLSQP',options={'ftol':1e-9,'maxiter':300})
 good=res.success and max(abs(cons(res.x)))<1e-6 and min(feasible(res.x))>=-1e-7
 if not good:
  retry=minimize(lambda t:float(t@t),res.x,jac=lambda t:2*t,constraints=[{'type':'eq','fun':cons,'jac':jac},{'type':'ineq','fun':feasible,'jac':fj}],method='SLSQP',options={'ftol':1e-7,'maxiter':1000})
  if retry.success and max(abs(cons(retry.x)))<1e-6 and min(feasible(retry.x))>=-1e-7:res=retry;good=True
 return float(res.fun),y+L@res.x,good
# Start with common total and observed auxiliary matrix; a feasible same-total model.
A,r=decode(x);qq,_=quantities(x);r[group]*=np.sqrt(qq[group].mean()/qq[group,None]);m=R.copy();m[6:,:6]=r;m[:6,6:]=r.T
obs,mu,ok=fit(x,m[ii,jj]);assert ok;print('observed',args.null,obs,'chiP',chi2.sf(obs,len(group)-1),flush=True)
if args.intervals:
 from scipy.optimize import brentq
 q0=quantities(x)[0];src=['PGC','FG',('Broad' if args.null=='Broad_minus_FG' else 'MVP')]
 targets=[('total',src[j],np.eye(3)[j],[j]) for j in range(3) if args.null=='global' or j==2]
 targets += [('difference',src[i]+'_minus_'+src[j],np.eye(3)[i]-np.eye(3)[j],[i,j]) for i,j in ([(0,1),(2,0),(2,1)] if args.null=='global' else [(2,1),(2,0)])]
 rows=[]
 for kind,name,contrast,selected in targets:
  D=contrast[None];group=selected;centre=float(contrast@q0);lowlim=0. if kind=='total' else -1.;uplim=1.
  def evaluate(value):
   global delta
   delta=float(value);val,point,valid=fit(x,x)
   if not valid:raise RuntimeError((name,value,val))
   return val-chi2.isf(.05,1)
  limits=[]
  for sign,lim in [(-1,lowlim),(1,uplim)]:
   prev=centre;bound=None
   for distance in [.025,.05,.1,.2,.35,.5,.75,1.,2.]:
    point=max(lowlim+1e-6,min(uplim-1e-6,centre+sign*distance));value=evaluate(point)
    if value>=0:bound=brentq(evaluate,min(prev,point),max(prev,point),xtol=1e-7);break
    prev=point
    if abs(point-lim)<2e-6:bound=lim;break
   assert bound is not None;limits.append(bound)
  rows.append(dict(kind=kind,name=name,initial_centre=centre,lower=limits[0],upper=limits[1],interval='conditional_Gaussian_profile_95pct_chi2_1'))
 pd.DataFrame(rows).to_csv(out/(args.null+'_profile_intervals.tsv'),sep='\t',index=False)
 raise SystemExit(0)
rng=np.random.default_rng(20260921);qs=[];fails=0;start=time.time()
for k in range(args.draws):
 y=mu+L@rng.normal(size=len(x));Q,_,ok=fit(y,mu)
 if not ok:
  # Fixed null mean is always a feasible fallback; retain failure as unknown, never silently discard.
  fails+=1;Q=np.nan
 qs.append(Q)
 if (k+1)%1000==0:print(k+1,'sec',time.time()-start,'fails',fails,flush=True)
qs=np.array(qs);rows=[]
for alpha in [.05,.01,.001]:
 n=int((qs>chi2.isf(alpha,len(group)-1)).sum());rows.append(dict(test=args.null,df=len(group)-1,alpha=alpha,draws=args.draws,failures=fails,rejections=n,rate_lower=n/args.draws,rate_upper=(n+fails)/args.draws))
tail=int((qs>=obs).sum());pd.DataFrame(rows).to_csv(out/f'{args.null}_rates.tsv',sep='\t',index=False);np.savez(out/f'{args.null}_LR.npz',Q=qs,mu=mu,V=V,observed=obs);(out/f'{args.null}_result.json').write_text(json.dumps(dict(test=args.null,observed=obs,df=len(group)-1,asymptotic_P=float(chi2.sf(obs,len(group)-1)),tail=tail,draws=args.draws,failures=fails,MC_P_lower=(tail+1)/(args.draws+1),MC_P_upper=(tail+fails+1)/(args.draws+1),seed=20260921,scope='Conditional Gaussian constrained minimum-distance test of equal population joint sharing; fixed empirical correlation covariance; no PHBC noisy-projection covariance'),indent=2));print(rows,flush=True)
