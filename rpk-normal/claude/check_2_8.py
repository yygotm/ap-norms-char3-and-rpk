# independent check of F_{2^8}: 3-primitive (order 85), 1-normal (rank 7 over F_2), Tr(xi^-1)=a, N(xi)=1
MOD=0x11d  # x^8+x^4+x^3+x^2+1
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1
        if a&0x100: a^=MOD
        b>>=1
    return r
def pw(a,n):
    r=1
    while n:
        if n&1: r=mul(r,a)
        a=mul(a,a); n>>=1
    return r
def order(a):
    x=a;n=1
    while x!=1: x=mul(x,a); n+=1
    return n
def rank(vs):
    vs=list(vs);r=0
    for bit in range(8):
        piv=next((i for i in range(r,len(vs)) if vs[i]>>bit&1),None)
        if piv is None: continue
        vs[r],vs[piv]=vs[piv],vs[r]
        for i in range(len(vs)):
            if i!=r and vs[i]>>bit&1: vs[i]^=vs[r]
        r+=1
    return r
from collections import Counter
cnt=Counter(); n3=0
for x in range(1,256):
    if order(x)!=85: continue
    n3+=1
    conj=[pw(x,2**j) for j in range(8)]
    if rank(conj)!=7: continue
    inv=pw(x,254)
    tr=0
    for j in range(8): tr^=pw(inv,2**j)
    N=pw(x,255)  # N over F_2 = x^(255)=1 always
    cnt[(tr,N)]+=1
print("3-primitive:",n3,"1-normal cells (Tr(xi^-1), N):",dict(cnt))
