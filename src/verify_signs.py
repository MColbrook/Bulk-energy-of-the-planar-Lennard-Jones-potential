"""Certify the Lennard--Jones signs using the original correction budgets.

The exact signed-mixture contacts and Fourier-origin zero are analytic
premises established in the accompanying manuscript.  Existence and uniform
correction estimates come from OpenAI's atomic Proposition 5.1; Lemmas 6.1
and 6.5 give the curvature envelopes and sine barriers used in the tail.
The pinned original source and its finite checker are in ../upstream/.

The interval calculation uses the outputs of lattice_constants.py and
reference_coefficients.py. Run ../run_verification.py for the complete sequence.
"""
if not __debug__:
    raise RuntimeError("Certificate checks require Python assertions; do not use -O or -OO.")
import argparse
import hashlib
import json
import time
from pathlib import Path
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--lattice-constants", required=True, type=Path)
parser.add_argument("--reference", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
from fractions import Fraction as F
from math import comb
from flint import arb,acb,fmpq,ctx
ctx.prec=160
started=time.time()

def rat(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def interval(a,b):
    a,b=F(a),F(b)
    return rat((a+b)/2)+arb(0,rat((b-a)/2).abs_upper())
def conv(a,b):
    out=[acb(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def log(message):print(f'{time.time()-started:.1f}s {message}',flush=True)

# The runner regenerates this reference in its fresh output directory.
data=json.loads(args.reference.read_text())
assert arb(data['absolute_mixture_mass'])<rat('.425')
delta_f=rat('4.25e-6');delta_t=rat('1.275e-9')
assert arb(data['inherited_finite_error'])<delta_f
assert arb(data['inherited_tail_error'])<delta_t
pi=arb.pi();I=acb(0,1);b=arb(3).sqrt()/2
h=rat('17/50');H=rat('27/50');eta=H-h
A=arb(data['A']);B=arb(data['B']);assert 0<A<B
res=[0,1,3,4,7,9,12,13,16,19,21,25,27,28,31]
tv=[rat(F(j-15,18)) for j in range(31)]
atoms=[[-(I*pi*a/18).exp(),acb(2),-(-I*pi*a/18).exp()] for a in res]
P=[acb(1)]
for q in atoms:P=conv(P,q)
cols=[]
for k in range(11):
    base=[acb(1)]
    for j in range(15):
        if j!=k:base=conv(base,atoms[j])
    cols.append([acb(0)]+[4*(pi/36)**2*z for z in base]+[acb(0)])
    cols.append(conv([I*pi/36*(I*pi*res[k]/18).exp(),acb(0),-I*pi/36*(-I*pi*res[k]/18).exp()],base))
cols.append(P)
co=[[arb(x) for x in data[key]] for key in ['u_reference_coefficients','v_reference_coefficients']]
spectral=[]
for side in range(2):
    fast=[sum((co[side][k]*cols[k][j] for k in range(22)),acb(0)) for j in range(31)]
    spectral.append([[co[side][22]*z for z in P],fast])
zd=[[-I*pi*t for t in tv],[pi*(rat('1/5')-I*t) for t in tv]]
zf=[[pi*(rat('4/3')/(k-I*t)-h) for t in tv] for k in [h,H]]
lf=[[1/(b*(k-I*t)) for t in tv] for k in [h,H]]
terms=[]
for side in range(2):
    sign=-1 if side==0 else 1
    row=[]
    for kind in range(2):
        for j in range(31):
            row.append((sign*spectral[side][kind][j],zd[kind][j]))
            row.append((sign*spectral[1-side][kind][j]*lf[kind][j],zf[kind][j]))
    terms.append(row)
assert all(z.real>=0 for row in zd+zf for z in row)
# These endpoint tests control the whole frequency interval, including
# continuous spectral measures of the tail columns. For X=(4/3)/(k*k+t*t),
# |z_f|^2=pi^2*(h*h+(4/3-2*k*h)*X); all relevant coefficients are positive.
for damping in [h,H]:
    assert rat('4/3')-2*damping*h>0
    assert 1/(b*damping)<4
    assert pi*(rat('4/3')/damping-h)<14
    assert pi*(rat('4/3')*damping/(damping**2+rat('5/6')**2)-h)>rat('3/5')
tail_tv=max(sum(P[j].abs_upper()*pi*abs(tv[j]) for j in range(31)).upper(),
            sum(P[j].abs_upper()*pi*pi*tv[j]**2/2 for j in range(31)).upper())
assert tail_tv<123

def moment(k,z):
    return z.gamma_upper(k+1)/z**(k+1) if k>=0 else z.expint(-k)

def base(s,side,order=0,absolute=False):
    z=pi*s/b;out=arb(0)
    for p,c in zip([5,2] if side==0 else [-6,-3],[A,-B]):
        for j in range(order+1):
            factor=c*comb(order,j)*(pi*h)**(order-j)*(-pi/b)**j
            out+=(abs(factor) if absolute else factor)*moment(p+j,z)
    return (pi*h*s).exp()*out

def value(s,side,order=0):
    out=base(s,side,order)
    for c,z in terms[side]:out+=(c*(-z)**order*(-z*s).exp()).real
    return out

def magnitude_bound(lo,side,order):
    out=base(rat(lo),side,order,True)
    for c,z in terms[side]:out+=c.abs_upper()*z.abs_upper()**order*(-z.real*rat(lo)).exp()
    return out.upper()

def correction_bound(lo,order):
    """Bound derivatives of e^(pi*h*s)*(true-reference) for all s>=lo."""
    lo=rat(lo)
    direct=[];transformed=[]
    for k in range(2,22):
        direct.append(sum((cols[k][j].abs_upper()*zd[1][j].abs_upper()**order*(-zd[1][j].real*lo).exp() for j in range(31)),arb(0)))
        transformed.append(sum((cols[k][j].abs_upper()*lf[1][j].abs_upper()*zf[1][j].abs_upper()**order*(-zf[1][j].real*lo).exp() for j in range(31)),arb(0)))
    finite=delta_f*(max(x.upper() for x in direct)+max(x.upper() for x in transformed))
    # Common-frequency measure bound; direct slow multiplication cancels h.
    # |lambda|<4, |pi*zeta|<14, pi*Im(zeta)>3/5 for transformed slow inputs.
    tail=delta_t*tail_tv*((pi*rat('5/6'))**order+4*arb(14)**order*(-rat('3/5')*lo).exp())
    return (finite+tail).upper()


def checked_cover(records, a, bnd):
    """Check exact rational coverage and retain the full finite certificate."""
    records = sorted(records, key=lambda item: F(item['left']))
    position = F(a)
    for item in records:
        left, right = F(item['left']), F(item['right'])
        assert left == position and left < right
        position = right
    assert position == F(bnd)
    return records


def prove_positive(a,bnd,side,label):
    """Taylor enclosures on adaptive rational intervals cover the whole region."""
    stack=[(F(a),F(bnd),0)];leaves=0;maxdepth=0;minimum=None
    accepted = []
    while stack:
        lo,hi,depth=stack.pop();mid=(lo+hi)/2;radius=rat((hi-lo)/2)
        v=value(rat(mid),side);d=value(rat(mid),side,1)
        taylor_error=magnitude_bound(lo,side,2)*radius**2/2
        direct_lower=v.lower()-d.abs_upper()*radius-taylor_error-correction_bound(lo,0)
        # Exact true residual jets vanish, but the reference jets need not.
        # Deleting the entire reference residual jet restores the true jets.
        m=min(contact_nodes,key=lambda n:abs(mid-n))
        v0,d0=reference_jets[side,m]
        deleted_v=v-v0-rat(mid-m)*d0
        deleted_d=d-d0
        distance=rat(max(abs(lo-m),abs(hi-m)))
        error=correction_bound(min(lo,m),2)*distance**2/2
        deleted_lower=deleted_v.lower()-deleted_d.abs_upper()*radius-taylor_error-error
        lower=max(direct_lower,deleted_lower)
        if lower>0:
            leaves+=1;maxdepth=max(maxdepth,depth)
            minimum=lower.lower() if minimum is None else min(minimum,lower.lower())
            accepted.append({'left':str(lo),'right':str(hi),'anchor':m,
                             'absolute_lower':str(direct_lower),
                             'deleted_lower':str(deleted_lower),
                             'certified_lower':str(lower)})
            continue
        if depth>=30:raise AssertionError((label,lo,hi,str(lower)))
        stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
    return {'label':label,'leaves':leaves,'max_depth':maxdepth,
            'minimum_lower_bound':str(minimum),'cover':checked_cover(accepted,a,bnd)}


def prove_node(m,side,halfwidth=F(1,16)):
    """Exact true value and slope are zero at m; prove positive true curvature."""
    a,bnd=F(m)-halfwidth,F(m)+halfwidth
    stack=[(a,bnd,0)];leaves=0;maxdepth=0;minimum=None
    accepted = []
    while stack:
        lo,hi,depth=stack.pop();mid=(lo+hi)/2;radius=rat((hi-lo)/2)
        lower=value(rat(mid),side,2).lower()-magnitude_bound(lo,side,3)*radius-correction_bound(lo,2)
        if lower>0:
            leaves+=1;maxdepth=max(maxdepth,depth)
            minimum=lower.lower() if minimum is None else min(minimum,lower.lower())
            accepted.append({'left':str(lo),'right':str(hi),
                             'curvature_lower':str(lower)})
            continue
        if depth>=25:raise AssertionError(('node',m,side,lo,hi,str(lower)))
        stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
    return {'node':m,'side':side,'leaves':leaves,'max_depth':maxdepth,
            'minimum_curvature':str(minimum),'cover':checked_cover(accepted,a,bnd)}


def origin_checks():
    # For s<=1/10, W(s)>=W(1/10)>13000; low-alpha part of h is nonpositive.
    # The two corrected inputs have combined coefficient norm<2 and each
    # basis spectral TV<123, while every Fourier Gaussian prefactor is<4.
    const=json.loads(args.lattice_constants.read_text())
    Z6=interval(const['Z6']['lower_exact'],const['Z6']['upper_exact'])
    Z12=interval(const['Z12']['lower_exact'],const['Z12']['upper_exact'])
    s=rat('1/10');W=1/(12*Z12*s**6)-1/(6*Z6*s**3)
    assert W>13000 and Z6/Z12>s**3
    assert sum(abs(c) for row in co for c in row)+2*(delta_f+delta_t)<2
    assert all(sum(z.abs_upper() for z in col)<123 for col in cols)
    assert 2*123*4<13000
    # Fourier: derivative of G_minus is positive. Retain only alpha in[1,10]
    # and lower-bound its exponential for 0<=s<=1/1000.
    delta=F(1,1000);M=arb(10)
    def integral_power(k):return (M**(k+1)-1)/(k+1)
    J=pi/b*(B*integral_power(-2)-A*integral_power(-5))-pi*h*(B*integral_power(-3)-A*integral_power(-6))
    lower=J*(-(pi*M/b-pi*h)*rat(delta)).exp()
    box=interval(0,delta)
    for c,z in terms[1]:lower+=(c*(-z)*(-z*box).exp()).real.lower()
    lower-=correction_bound(0,1)
    assert lower>rat('.9396')
    return {'real_interval':'(0,1/10]','real_W_lower_bound':str(W),'Fourier_interval':'[0,1/1000]','Fourier_derivative_lower_bound':str(lower)}



contact_nodes=[n for n in range(1,73) if n%36 in res]
reference_jets={(side,m):(value(rat(m),side),value(rat(m),side,1))
                for side in [0,1] for m in contact_nodes if m<=16}

def inherited_tail_checks():
    # Source Lemma 6.1 controls both lists with no further factor of two.
    # The fixed slow columns are +/-gamma and are not part of finite_poles.
    finite_poles=max((sum(abs(c) for c in row[:22])+delta_f).upper() for row in co)
    gamma=arb(data['slow_coefficient'])
    assert finite_poles<rat('.619')
    assert rat('.001436')<gamma<rat('.001437')
    # Every summand in this absolute target-curvature envelope decreases
    # on s>=14.5, since pi*alpha/b-pi*h>0 for alpha>=1.
    target_curvatures=[base(rat('14.5'),side,2,True)/2 for side in [0,1]]
    assert all(v<rat('2e-17') for v in target_curvatures)
    K1=rat('.619')*rat('.00465')+800*delta_t+rat('.001437')*rat('.002')+rat('2e-17')
    K2=rat('.619')*rat('.000023')+800*delta_t+rat('.001437')*rat('.002')+rat('2e-17')
    barrier1=rat('.001436')*3
    barrier2=rat('.001436')*rat('.054')
    assert K1<barrier1 and K2<barrier2
    # For each residual F, F=gamma*P+D and F(m)=F'(m)=P(m)=P'(m)=0.
    # Hence F(s)>=gamma*P(s)-K*(s-m)^2 on a nearest-node segment.
    # The source chooses 16,19,21 for [14.5,23], and any nearest node
    # above 23; the entire segment stays in the corresponding range.
    return {'source_lemmas':['6.1','6.5'],
            'finite_pole_norm_bound':str(finite_poles),
            'slow_coefficient':str(gamma),
            'target_half_curvatures_at_14.5':[str(v) for v in target_curvatures],
            'range_14.5_to_23':{'half_curvature':str(K1),'barrier':str(barrier1),'margin':str(barrier1-K1)},
            'range_above_23':{'half_curvature':str(K2),'barrier':str(barrier2),'margin':str(barrier2-K2)}}

if __name__=='__main__':
    result={'scope':'Target-specific sign inequalities conditional on the stated analytic premises',
            'input_sha256':{'reference':hashlib.sha256(args.reference.read_bytes()).hexdigest(),
                            'lattice_constants':hashlib.sha256(args.lattice_constants.read_bytes()).hexdigest()},
            'absolute_mixture_mass':data['absolute_mixture_mass'],
            'finite_budget':str(delta_f),'tail_budget':str(delta_t),
            'origin':origin_checks(),'node_neighborhoods':[],'remaining_intervals':[]}
    log('origin checks passed using reference norms plus inherited errors')
    for side in [0,1]:
        left=F(1,10) if side==0 else F(1,1000)
        for m in [1,3,4,7,9,12,13]:
            lo=F(m)-F(1,16);hi=F(m)+F(1,16)
            if left<lo:
                item=prove_positive(left,lo,side,f'side{side}:{left}..{lo}')
                result['remaining_intervals'].append(item)
            item=prove_node(m,side)
            result['node_neighborhoods'].append(item)
            left=hi
            log(f'side{side} passed through node{m}; {item["leaves"]} curvature intervals')
        result['remaining_intervals'].append(prove_positive(left,F(29,2),side,f'side{side}:last..14.5'))
    result['tail']=inherited_tail_checks()
    # Verify the complete compact partition on each side, in addition to
    # the leaf-by-leaf checks within each separately certified interval.
    for side in [0,1]:
        cover=[]
        for item in result['remaining_intervals']:
            if item['label'].startswith(f'side{side}:'):
                cover.extend(item['cover'])
        for item in result['node_neighborhoods']:
            if item['side']==side:
                cover.extend(item['cover'])
        checked_cover(cover,F(1,10) if side==0 else F(1,1000),F(29,2))
    result['finite_checks_passed']=True
    result['curvature_intervals']=sum(item['leaves'] for item in result['node_neighborhoods'])
    result['other_taylor_intervals']=sum(item['leaves'] for item in result['remaining_intervals'])
    path=args.output
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    log(f'PASS: original coarse errors suffice; {result["curvature_intervals"]} curvature and {result["other_taylor_intervals"]} other Taylor intervals; inherited tail from 14.5')
