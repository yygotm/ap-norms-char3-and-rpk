"""Independent small cross-check: k-normality via polynomial gcd, not matrix rank.
No producer sources/logs are read. C++ fixture records are compared directly.
"""
from pathlib import Path
import struct,time,json,math
ROOT=Path(__file__).resolve().parent
import os
DATA=Path(os.environ.get('DATA', str(ROOT.parent / 'witnesses')))
class F:
    def __init__(self,p,r):
        self.p=p;self.r=r;self.n=len(r);self.size=p**self.n
        self.z=0 if p==2 else (0,)*self.n
        self.o=1 if p==2 else (1,)+(0,)*(self.n-1)
        self.mod=(1<<self.n)|sum(v<<i for i,v in enumerate(r)) if p==2 else None
    def elem(self,v):return sum(x<<i for i,x in enumerate(v)) if self.p==2 else tuple(v)
    def add(self,a,b):return a^b if self.p==2 else tuple((x+y)%self.p for x,y in zip(a,b))
    def sub(self,a,b):return a^b if self.p==2 else tuple((x-y)%self.p for x,y in zip(a,b))
    def mul(self,a,b):
        if self.p==2:
            z=0
            while b:
                if b&1:z^=a
                a<<=1;b>>=1
            while z.bit_length()>self.n:z^=self.mod<<(z.bit_length()-self.n-1)
            return z
        p=self.p;n=self.n;v=[0]*(2*n-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):v[i+j]+=x*y
        for i in range(2*n-2,n-1,-1):
            c=v[i]%p
            for j,x in enumerate(self.r):v[i-n+j]+=c*x
        return tuple(x%p for x in v[:n])
    def power(self,a,n):
        z=self.o
        # left-to-right exponentiation differs from main checker's implementation
        for b in bin(n)[2:]:
            z=self.mul(z,z)
            if b=='1':z=self.mul(z,a)
        return z
    def inverse(self,a):
        assert a!=self.z
        z=self.power(a,self.size-2);assert self.mul(z,a)==self.o
        return z
    def gcd_degree(self,a,b):
        def trim(v):
            while v and v[-1]==self.z:v.pop()
            return v
        a=trim(a[:]);b=trim(b[:])
        while b:
            iv=self.inverse(b[-1]);r=a[:]
            while len(r)>=len(b):
                d=self.mul(r[-1],iv);k=len(r)-len(b)
                for j,x in enumerate(b):r[k+j]=self.sub(r[k+j],self.mul(d,x))
                trim(r)
            a,b=b,r
        return len(a)-1
    def k(self,x,q,m):
        conj=[];y=x
        for _ in range(m):conj.append(y);y=self.power(y,q)
        assert y==x
        return self.gcd_degree([self.sub(self.z,self.o)]+[self.z]*(m-1)+[self.o],list(reversed(conj)))
    def cell(self,x,q,m):
        b=self.power(x,(self.size-1)//(q-1));y=self.inverse(x);t=self.z
        for _ in range(m):t=self.add(t,y);y=self.power(y,q)
        return self.mul(t,b),b,t
def run():
    start=time.perf_counter();cs={}
    for line in (ROOT/'cases.txt').read_text().splitlines():
        z=line.split();kind=z[0];q,m,p,e,cells,nf=map(int,z[1:7]);fac=list(map(int,z[7:7+nf]));r=list(map(int,z[7+nf:]));cs[f'{q}_{m}']=(kind,q,m,p,e,cells,fac,r)
    rows=[]
    for line in (ROOT/'all-full.log').read_text().splitlines():
        if not line.startswith('FIXTURE '):continue
        _,pair,*fields=line.split();v=dict(x.split('=',1) for x in fields)
        kind,q,m,p,e,cells,fac,r=cs[pair];f=F(p,r);i=int(v['index'])
        with (DATA/f'{pair}.bin').open('rb') as h:h.seek(i*2*f.n);x=f.elem(struct.unpack('<'+'H'*f.n,h.read(2*f.n)))
        a,b,t=f.cell(x,q,m);ka=f.k(x,q,m)
        expected_a=f.elem(list(map(int,v['a'].strip(',').split(','))))
        expected_b=f.elem(list(map(int,v['b'].strip(',').split(','))))
        order=(f.size-1)//3;mf=fac[:];mf.remove(3)
        assert f.power(x,order)==f.o and all(f.power(x,order//l)!=f.o for l in set(mf))
        assert a==expected_a and b==expected_b and e*(m-ka)==int(v['rank']) and ka==1,(pair,i)
        rows.append(dict(pair=pair,index=i,k=ka))
        if len(rows)%15==0:print('reference B checked',len(rows),'elapsed',round(time.perf_counter()-start,2),flush=True)
    assert len(rows)==171
    # Direct coordinate enumeration, with no primitive-power enumeration or orbit reduction.
    red=[1,0,0,1,0,1,0,1];f=F(2,red);dist={0:0,1:0};witness=[]
    for x in range(1,256):
        if f.power(x,85)!=1 or f.power(x,17)==1 or f.power(x,5)==1:continue
        if f.k(x,2,8)!=1:continue
        a,b,t=f.cell(x,2,8);assert b==1 and a==t and t in dist;dist[t]+=1;witness.append(x)
    assert dist=={0:0,1:8}
    # A concrete case satisfying the theorem's printed sufficient inequality.
    from prepare import prime
    factors=[3,5,17,257,641,65537,274177,6700417,67280421310721]
    assert math.prod(factors)==4**64-1 and all(prime(p) for p in factors)
    result=dict(B_checked=len(rows),B_failed=0,method='gcd(X^m-1, reversed Frobenius conjugates) over extension field',A_2_8_total=256,A_2_8_elements=witness,A_2_8_inverse_trace=dist,theorem_counterexample=dict(q=4,m=64,r=3,k=1,g='X-1',factors=factors,lhs=4**29,rhs=3*2**8*2),elapsed_s=time.perf_counter()-start)
    (ROOT/'reference-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':run()
