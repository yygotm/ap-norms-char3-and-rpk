"""Read ONLY delivered witnesses; independently factor orders with deterministic MR64."""
from pathlib import Path
import hashlib, json, math, random, re
ROOT=Path(__file__).resolve().parent
import os
DATA=Path(os.environ.get('DATA', str(ROOT.parent / 'witnesses')))
A=[(2,m) for m in (8,10,12,14,16,18)]+[(4,m) for m in (7,8,9,10,12)]+[(5,m) for m in (8,10,12)]+[(7,m) for m in (7,8,9)]+[(8,8),(11,8),(13,7),(16,7)]
ALL=set(A)|{(7,12),(11,10)}|{(q,8) for q in (8,11,13,16,17,19,23,25,29,32,37,41,43,47,49,53,59,83,89)}|{(q,7) for q in (13,16,19,25,31,37,43,49,61,64,67,73,79,97,103,109,121,127,139,151,169,181,193,199,211,223,256,277,307,337,361,379,421,433,463,529,547,631,673,883)}
def prime(n):
    assert n<2**64
    if n<2:return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0:return n==p
    d=n-1;s=0
    while d%2==0:d//=2;s+=1
    for a in (2,325,9375,28178,450775,9780504,1795265022):
        if a%n==0:continue
        x=pow(a,d,n)
        if x in (1,n-1):continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1:break
        else:return False
    return True
rng=random.Random(32020261010)
def factor(n):
    if n==1:return []
    if prime(n):return [n]
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0:return [p]+factor(n//p)
    while True:
        c=rng.randrange(1,n);x=rng.randrange(2,n);y=x;d=1
        while d==1:
            x=(x*x+c)%n;y=(y*y+c)%n;y=(y*y+c)%n
            d=math.gcd(abs(x-y),n)
        if d!=n:return factor(d)+factor(n//d)
def run():
    delivered={tuple(map(int,f.stem.split('_'))) for f in DATA.glob('*.hdr')}
    assert len(ALL)==78 and len(A)==21 and delivered==ALL-set(A)
    rows=[]; manifest=[]
    for q,m in sorted(ALL):
        Q=q**m-1
        fac=sorted(factor(q-1)+factor(Q//(q-1)))
        assert math.prod(fac)==Q and all(prime(f) for f in fac)
        if (q,m) in A:
            p=next(p for p in (2,5,7,11,13) if q%p==0);e=0;t=q
            while t%p==0:e+=1;t//=p
            assert t==1
            cells=q*math.prod([]) # replaced by exact totient below
            L=(q-1)//math.gcd(3,q-1);phi=L
            for l in set(factor(L)):phi=phi//l*(l-1)
            rows.append(f'A {q} {m} {p} {e} {q*phi} {len(fac)} '+ ' '.join(map(str,fac)))
        else:
            h=DATA/f'{q}_{m}.hdr';b=h.with_suffix('.bin');txt=h.read_text()
            vals={k:int(v) for k,v in re.findall(r'^(p|e|m|N|cells) (\d+)$',txt,re.M)}
            assert vals['p']**vals['e']==q and vals['m']==m and vals['N']==vals['e']*m
            red=list(map(int,txt.split('red_low_to_high:')[1].splitlines()[0].split()))
            assert len(red)==vals['N'] and all(0<=v<vals['p'] for v in red)
            assert b.stat().st_size==2*vals['N']*vals['cells']
            rows.append(f'B {q} {m} {vals["p"]} {vals["e"]} {vals["cells"]} {len(fac)} '+' '.join(map(str,fac))+' '+' '.join(map(str,red)))
            manifest.append(dict(pair=f'{q}_{m}',cells=vals['cells'],bytes=b.stat().st_size,header_sha256=hashlib.sha256(h.read_bytes()).hexdigest(),binary_sha256=hashlib.sha256(b.read_bytes()).hexdigest()))
    (ROOT/'cases.txt').write_text('\n'.join(rows)+'\n',encoding='ascii')
    (ROOT/'input-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('A',len(A),'B',len(manifest),'B records',sum(x['cells'] for x in manifest),'bytes',sum(x['bytes'] for x in manifest),flush=True)
if __name__=='__main__':run()
