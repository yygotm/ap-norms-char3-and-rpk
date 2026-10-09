import sys
from sympy import factorint, totient, gcd, perfect_power
def pe(q):
    r=perfect_power(q)
    if r: return int(r[0]),int(r[1])
    return q,1
pairs=[(2,8)]+[(13,7),(7,8),(16,7),(5,10),(4,9)]
done={(2,10),(2,12),(2,14),(2,16),(2,18),(4,7),(4,8),(4,9),(4,10),(4,12),(5,8),(5,10),(5,12),(7,7),(7,8),(7,9),(8,8),(11,8),(13,7),(16,7)}
allp=[(4,m) for m in (7,8,9,10,12)]+[(5,m) for m in (8,10,12)]+[(7,m) for m in (7,8,9,12)]+[(2,m) for m in (8,10,12,14,16,18)]
allp+=[(q,8) for q in (8,11,13,16,17,19,23,25,29,32,37,41,43,47,49,53,59,83,89)]+[(11,10)]
allp+=[(q,7) for q in (13,16,19,25,31,37,43,49,61,64,67,73,79,97,103,109,121,127,139,151,169,181,193,199,211,223,256,277,307,337,361,379,421,433,463,529,547,631,673,883)]
allp=sorted(set(allp)); print(len(allp),"pairs in list")
import os
for q,m in set(allp)|set(pairs):
    p,e=pe(q); Q=q**m-1; f=factorint(Q); w=int(gcd(3,q-1)); C=q*int(totient((q-1)//w))
    open(f"params/{q}_{m}.txt","w").write(" ".join(map(str,[p,e,m,len(f),*sorted(f),C]))+"\n")
todo=[(q,m) for (q,m) in allp if (q,m) not in done]
print(len(todo),"to do")
open("todo.txt","w").write("\n".join(f"{q}_{m}" for q,m in sorted(todo,key=lambda t:t[0]**t[1]))+"\n")
