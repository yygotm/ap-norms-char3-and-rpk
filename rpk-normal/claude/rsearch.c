// rsearch.c — candidate 32, random search for cells (a,b): 3-primitive 1-normal xi in F_{q^m}, q=p^e,
// N(xi)=b, Tr(xi^-1)=a*b^-1 (key = (a,b) with a=Tr(xi^-1)*N(xi)).  Cells: q*phi((q-1)/w), w=gcd(3,q-1).
// usage: rsearch params.txt budget_samples [threads]
// params.txt: p e m nprimes prime1 ... (distinct primes of Q=q^m-1, as decimal u128), then C (number of cells)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef uint64_t u64; typedef uint32_t u32; typedef unsigned __int128 u128;
#define MAXN 64
static int p,e,m,N; static int Fm[MAXN+1]; // monic F: x^N = -sum_{j<N} Fm[j] x^j  stored as red[j] = -Fm[j] mod p
static int red[MAXN]; static int nz[MAXN], nnz; // nonzero positions of red
static int invp[1024];
static int Fr[MAXN][MAXN]; // Frobenius x->x^q on basis
static int Fp[MAXN][MAXN]; // absolute Frobenius x->x^p
static void frobp(const int*v,int*r){ int t[MAXN]={0}; for(int i=0;i<N;i++){ if(!v[i]) continue; for(int j=0;j<N;j++) t[j]=(int)(((u64)t[j]+(u64)v[i]*Fp[i][j])%p);} memcpy(r,t,sizeof(int)*N); }
static int theta[MAXN]; static u128 Q, QD3;
static int nprimes; static u128 primes[40];

static u128 parse128(const char*s){ u128 v=0; for(;*s;s++) v=v*10+(*s-'0'); return v; }
static void mulmod(const int*a,const int*b,int*r){
  int t[2*MAXN]; memset(t,0,sizeof(int)*(2*N));
  for(int i=0;i<N;i++){ if(!a[i]) continue; for(int j=0;j<N;j++) t[i+j]=(int)(((u64)t[i+j]+(u64)a[i]*b[j])%p); }
  for(int i=2*N-2;i>=N;i--){ int c=t[i]; if(!c) continue; t[i]=0; for(int k=0;k<nnz;k++){ int j=nz[k]; t[i-N+j]=(int)(((u64)t[i-N+j]+(u64)c*red[j])%p); } }
  memcpy(r,t,sizeof(int)*N);
}
static int iszero(const int*a){ for(int i=0;i<N;i++) if(a[i]) return 0; return 1; }
static int isone(const int*a){ if(a[0]!=1) return 0; for(int i=1;i<N;i++) if(a[i]) return 0; return 1; }
static void powmod(const int*a,u128 ex,int*r){
  int res[MAXN]={0}, base[MAXN]; res[0]=1; memcpy(base,a,sizeof(int)*N);
  while(ex){ if(ex&1) mulmod(res,base,res); ex>>=1; if(ex) mulmod(base,base,base); }
  memcpy(r,res,sizeof(int)*N);
}
static void frob(const int*v,int*r){ int t[MAXN]={0}; for(int i=0;i<N;i++){ if(!v[i]) continue; for(int j=0;j<N;j++) t[j]=(int)(((u64)t[j]+(u64)v[i]*Fr[i][j])%p);} memcpy(r,t,sizeof(int)*N); }
static int rank_rows(int (*A)[MAXN], int nv){
  int r=0;
  for(int c=0;c<N&&r<nv;c++){
    int piv=-1; for(int i=r;i<nv;i++) if(A[i][c]){piv=i;break;}
    if(piv<0) continue;
    if(piv!=r) for(int j=c;j<N;j++){int t=A[piv][j];A[piv][j]=A[r][j];A[r][j]=t;}
    int iv=invp[A[r][c]];
    for(int i=r+1;i<nv;i++) if(A[i][c]){ int f=(int)(((u64)A[i][c]*iv)%p); for(int j=c;j<N;j++) A[i][j]=(int)((A[i][j]+(u64)(p-f)*A[r][j])%p); }
    r++;
  }
  return r;
}
// ---- irreducibility of F (Ben-Or) using polys over F_p
static int pdeg(const int*a,int n){ for(int i=n-1;i>=0;i--) if(a[i]) return i; return -1; }
static void pgcd(int*A,int na,int*B,int nb,int*out,int*dout){ // inputs arrays of length <=MAXN+1, destroyed
  int a[MAXN+2],b[MAXN+2]; memcpy(a,A,sizeof(int)*(MAXN+2)); memcpy(b,B,sizeof(int)*(MAXN+2));
  int da=pdeg(a,MAXN+2), db=pdeg(b,MAXN+2);
  while(db>=0){
    // a = a mod b
    while(da>=db){ int f=(int)(((u64)a[da]*invp[b[db]])%p); for(int j=0;j<=db;j++) a[da-db+j]=(int)((a[da-db+j]+(u64)(p-f)*b[j])%p); da=pdeg(a,MAXN+2); if(da<0) break; }
    int t[MAXN+2]; memcpy(t,a,sizeof t); memcpy(a,b,sizeof t); memcpy(b,t,sizeof t); int td=da; da=db; db=td;
  }
  memcpy(out,a,sizeof(int)*(MAXN+2)); *dout=da;
}
static int irreducible(){
  int x[MAXN]={0}; if(N==1) return 1; x[1]=1; int cur[MAXN]; memcpy(cur,x,sizeof cur);
  for(int i=1;i<=N/2;i++){
    powmod(cur,p,cur); // x^{p^i}
    int d[MAXN+2]={0}; for(int j=0;j<N;j++) d[j]=cur[j]; d[1]=(d[1]+p-1)%p; // x^{p^i}-x
    int Fp[MAXN+2]={0}; for(int j=0;j<N;j++) Fp[j]=(p-red[j])%p; Fp[N]=1;
    int g[MAXN+2],dg;
    if(pdeg(d,MAXN+2)<0) return 0;
    pgcd(d,0,Fp,0,g,&dg);
    if(dg!=0) return 0;
  }
  return 1;
}
static u64 rs[64][4];
static inline u64 rotl(u64 x,int k){return (x<<k)|(x>>(64-k));}
static u64 rnd(int t){ u64*s=rs[t]; u64 r=rotl(s[1]*5,7)*9, tt=s[1]<<17; s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=tt; s[3]=rotl(s[3],45); return r; }
// hash set
static u64 *HS; static u64 HSZ; static u64 hcount=0; static FILE*wf=NULL;
static u64 hashvec(const int*a,const int*b){ u64 h=1469598103934665603ULL; for(int i=0;i<N;i++){ h^=(u64)a[i]+1; h*=1099511628211ULL; } for(int i=0;i<N;i++){ h^=(u64)b[i]+7; h*=1099511628211ULL; h^=h>>29; } if(!h) h=1; return h; }
static int hinsert(u64 h){ u64 i=h&(HSZ-1); while(HS[i]){ if(HS[i]==h) return 0; i=(i+1)&(HSZ-1);} HS[i]=h; hcount++; return 1; }

int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: rsearch params budget [threads]\n");return 1;}
  FILE*f=fopen(argv[1],"r"); if(!f){perror("params");return 1;}
  char buf[128]; u64 C;
  fscanf(f,"%d %d %d %d",&p,&e,&m,&nprimes); for(int i=0;i<nprimes;i++){ fscanf(f,"%127s",buf); primes[i]=parse128(buf);} fscanf(f,"%llu",(unsigned long long*)&C); fclose(f);
  double budget=atof(argv[2]); int nth=argc>3?atoi(argv[3]):8; omp_set_num_threads(nth);
  N=e*m; for(int a=1;a<p;a++) for(int b=1;b<p;b++) if((a*b)%p==1) invp[a]=b;
  u128 qq=1; for(int i=0;i<e;i++) qq*=p; u128 P=1; for(int i=0;i<N;i++) P*=p; Q=P-1; QD3=Q/3;
  if(Q%3){ fprintf(stderr,"3 does not divide Q\n"); return 1; }
  double t0=omp_get_wtime();
  // find irreducible F: trinomials x^N + a x^k + b first, then random
  srand(777); int found=0;
  for(int tries=0;tries<2000000&&!found;tries++){
    memset(red,0,sizeof red);
    if(tries<20000){ int k=1+rand()%(N>1?N-1:1); if(N==1) k=0; red[0]=1+rand()%(p-1); if(N>1) red[k]=(red[k]+rand()%p)%p; }
    else { for(int j=0;j<N;j++) red[j]=rand()%p; if(!red[0]) continue; }
    nnz=0; for(int j=0;j<N;j++) if(red[j]) nz[nnz++]=j;
    if(irreducible()) found=1;
  }
  if(!found){ fprintf(stderr,"no irreducible poly found\n"); return 1; }
  // Frobenius matrix
  { int x1[MAXN]={0}; x1[1%N]=1; if(N==1) {x1[0]=red[0];} int xq[MAXN]; powmod(x1,qq,xq);
    int cur[MAXN]={0}; cur[0]=1; for(int i=0;i<N;i++){ memcpy(Fr[i],cur,sizeof(int)*N); mulmod(cur,xq,cur); }
    int xp[MAXN]; powmod(x1,p,xp); int cu2[MAXN]={0}; cu2[0]=1; for(int i=0;i<N;i++){ memcpy(Fp[i],cu2,sizeof(int)*N); mulmod(cu2,xp,cu2); } }
  // theta in F_q generating F_q over F_p
  srand(4242);
  for(;;){
    int u[MAXN]; for(int j=0;j<N;j++) u[j]=rand()%p; if(iszero(u)) continue;
    int t[MAXN]; powmod(u,Q/(qq-1),t); if(iszero(t)) continue;
    int okk=1; for(int l=2;l<=e;l++){ if(e%l) continue; int isp=1; for(int d=2;d*d<=l;d++) if(l%d==0) isp=0; if(!isp) continue;
      int s[MAXN]; memcpy(s,t,sizeof s); for(int i=0;i<e/l;i++) powmod(s,p,s); if(memcmp(s,t,sizeof(int)*N)==0){okk=0;break;} }
    if(e==1){ memset(theta,0,sizeof theta); theta[0]=1; break; }
    if(okk){ memcpy(theta,t,sizeof theta); break; }
  }
  if(argc>4){ char fn[512]; snprintf(fn,sizeof fn,"%s.bin",argv[4]); wf=fopen(fn,"wb");
    snprintf(fn,sizeof fn,"%s.hdr",argv[4]); FILE*h=fopen(fn,"w");
    fprintf(h,"p %d\ne %d\nm %d\nN %d\ncells %llu\n",p,e,m,N,(unsigned long long)C);
    fprintf(h,"modulus: x^N = sum_{j<N} red[j] x^j over F_p; red_low_to_high:"); for(int j=0;j<N;j++) fprintf(h," %d",red[j]); fprintf(h,"\n");
    fprintf(h,"record: N little-endian uint16 digits (coefficient of x^0 first) of xi in F_p[x]/(x^N - sum red[j]x^j); one record per cell\n"); fclose(h); }
  HSZ=1; while(HSZ<4*C+1024) HSZ<<=1; HS=calloc(HSZ,sizeof(u64));
  for(int t=0;t<64;t++){ u64 s=0x9E3779B97F4A7C15ULL*(t+1)^0xDEADBEEF; for(int k=0;k<4;k++){ s+=0x9E3779B97F4A7C15ULL; u64 z=s; z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; rs[t][k]=z^(z>>31);} }
  // primes of Q/3
  u128 plist[40]; int np2=0;
  for(int i=0;i<nprimes;i++){ if(primes[i]==3){ if((QD3%3)==0) plist[np2++]=3; } else plist[np2++]=primes[i]; }
  u64 samples=0, successes=0, nonFq=0, outside=0;
  int stop=0;
  #pragma omp parallel
  {
    int tid=omp_get_thread_num();
    u64 mysamp=0, mysucc=0;
    while(!__atomic_load_n(&stop,__ATOMIC_RELAXED)){
      for(int rep=0;rep<256&&!__atomic_load_n(&stop,__ATOMIC_RELAXED);rep++){
        int eps[MAXN]; for(int j=0;j<N;j++) eps[j]=(int)(rnd(tid)%p);
        if(iszero(eps)) continue;
        int xi[MAXN]; mulmod(eps,eps,xi); mulmod(xi,eps,xi); mysamp++;
        // 1-normality
        int V[MAXN][MAXN]; memcpy(V[0],xi,sizeof(int)*N);
        for(int j=1;j<m;j++) frob(V[j-1],V[j]);
        int A[MAXN][MAXN]; int nv=0;
 { int th[MAXN]={0}; th[0]=1; for(int l=0;l<e;l++){ for(int j=0;j<m;j++){ mulmod(th,V[j],A[nv]); nv++; } if(l+1<e) mulmod(th,theta,th); } }
        if(rank_rows(A,nv)!=N-e) continue;
        // 3-primitive (xi is a cube, so xi^{Q/3}=1 automatically)
        int okp=1; for(int i=0;i<np2&&okp;i++){ int z[MAXN]; powmod(xi,QD3/plist[i],z); if(isone(z)) okp=0; }
        if(!okp) continue;
        // norm, trace of inverse
        int bN[MAXN]; powmod(xi,Q/(qq-1),bN);
        int y[MAXN]; powmod(xi,Q-1,y);
        int tr[MAXN]; memcpy(tr,y,sizeof tr); { int cur[MAXN]; memcpy(cur,y,sizeof cur); for(int j=1;j<m;j++){ frob(cur,cur); for(int k=0;k<N;k++) tr[k]=(tr[k]+cur[k])%p; } }
        // check tr, bN in F_q: x^q==x
        int chk[MAXN]; frob(tr,chk); int bad=memcmp(chk,tr,sizeof(int)*N)!=0; frob(bN,chk); bad|=memcmp(chk,bN,sizeof(int)*N)!=0;
        if(bad){
          #pragma omp atomic
          nonFq++;
          continue; }
        int a[MAXN]; mulmod(tr,bN,a);
        mysucc++;
        { int aa[MAXN],bb[MAXN],xt[MAXN]; memcpy(aa,a,sizeof aa); memcpy(bb,bN,sizeof bb); memcpy(xt,xi,sizeof xt);
          for(int t=0;t<e;t++){ u64 h=hashvec(aa,bb);
            #pragma omp critical
            { if(hinsert(h) && wf){ uint16_t o[MAXN]; for(int k=0;k<N;k++) o[k]=(uint16_t)xt[k]; fwrite(o,sizeof(uint16_t),N,wf); } if(hcount>=C) __atomic_store_n(&stop,1,__ATOMIC_RELAXED); }
            frobp(aa,aa); frobp(bb,bb); frobp(xt,xt); } }
      }
      if(mysamp>(u64)(budget/nth)) break;
    }
    #pragma omp atomic
    samples+=mysamp;
    #pragma omp atomic
    successes+=mysucc;
  }
  printf("p=%d e=%d m=%d (q^m): cells=%llu found=%llu samples=%llu successes=%llu nonFq=%llu time=%.1fs %s\n",
    p,e,m,(unsigned long long)C,(unsigned long long)hcount,(unsigned long long)samples,(unsigned long long)successes,(unsigned long long)nonFq,omp_get_wtime()-t0,
    hcount>=C?"ALL_COVERED":(hcount>C?"MORE_THAN_C(!)":"INCOMPLETE"));
  if(wf) fclose(wf);
  return 0;
}
