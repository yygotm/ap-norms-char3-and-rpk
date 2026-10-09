// brute.c — candidate 32: 3-primitive 1-normal xi in F_{q^m} with N(xi)=b, Tr(xi^-1)=a*b^-1, q=p^e.
// For each cell (a,b), b of order (q-1)/gcd(3,q-1), look for one xi. Full pass over xi = g^i, gcd(i,Q)=3.
// usage: brute p e m [selftest]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef uint64_t u64; typedef uint32_t u32;
static int p, e, m, N; static u64 Q, P; static u32 *EXP, *LOG; static int red[40];
static u64 gcd64(u64 a, u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
static u64 enc(const int*d){ u64 v=0; for(int j=N-1;j>=0;j--) v=v*p+d[j]; return v; }
static void dec(u64 v,int*d){ for(int j=0;j<N;j++){ d[j]=v%p; v/=p; } }
static int rank_vecs(int (*A)[40], int nv){ // gaussian elimination mod p, nv rows of length N
  int r=0; static const int inv5[8]={0};
  for(int c=0;c<N&&r<nv;c++){
    int piv=-1; for(int i=r;i<nv;i++) if(A[i][c]){piv=i;break;}
    if(piv<0) continue;
    if(piv!=r) for(int j=0;j<N;j++){int t=A[piv][j];A[piv][j]=A[r][j];A[r][j]=t;}
    int iv=1; for(int x=1;x<p;x++) if((A[r][c]*x)%p==1){iv=x;break;}
    for(int j=c;j<N;j++) A[r][j]=(A[r][j]*iv)%p;
    for(int i=0;i<nv;i++) if(i!=r&&A[i][c]){ int f=A[i][c]; for(int j=c;j<N;j++) A[i][j]=((A[i][j]-f*A[r][j])%p+p)%p; }
    r++;
  }
  (void)inv5; return r;
}
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage\n");return 1;}
  p=atoi(argv[1]); e=atoi(argv[2]); m=atoi(argv[3]); int selftest=argc>4;
  N=e*m; P=1; for(int i=0;i<N;i++) P*=p; Q=P-1;
  u64 q=1; for(int i=0;i<e;i++) q*=p;
  double t0=omp_get_wtime();
  EXP=malloc(sizeof(u32)*Q); LOG=calloc(P,sizeof(u32));
  srand(12345); int tries=0;
  for(;;){
    tries++;
    for(int j=0;j<N;j++) red[j]=rand()%p; if(!red[0]) continue;
    int cur[40]={0}; cur[0]=1; int ok=1; EXP[0]=1;
    for(u64 t=1;t<=Q;t++){
      int top=cur[N-1]; for(int j=N-1;j>0;j--) cur[j]=cur[j-1]; cur[0]=0;
      if(top) for(int j=0;j<N;j++) cur[j]=(cur[j]+top*red[j])%p;
      u64 v=enc(cur);
      if(v==1){ if(t!=Q) ok=0; break; }
      if(t<Q) EXP[t]=(u32)v;
    }
    if(ok) break;
  }
  for(u64 t=0;t<Q;t++) LOG[EXP[t]]=(u32)t;
  fprintf(stderr,"tables built (%d poly tries) %.1fs\n",tries,omp_get_wtime()-t0);
  u64 step=Q/(q-1); u64 w=gcd64(3,q-1);
  u64 qj[40]; qj[0]=1; for(int j=1;j<m;j++) qj[j]=(qj[j-1]*q)%Q;
  if(selftest){
    u64 c0=0,c1=0;
    for(u64 i=0;i<Q;i++){
      int A[40][40]; int nv=0;
      for(int l=0;l<e;l++) for(int j=0;j<m;j++){ u64 ex=((u64)l*step+i*qj[j])%Q; dec(EXP[ex],A[nv]); nv++; }
      int r=rank_vecs((int(*)[40])A,nv);
      if(r==N) c0++; if(r==N-e) c1++;
    }
    printf("selftest p=%d e=%d m=%d: normal=%llu 1-normal=%llu\n",p,e,m,(unsigned long long)c0,(unsigned long long)c1); return 0;
  }
  u64 cells=q*(q-1); unsigned char*cov=calloc(cells,1);
  u64 ncells=0; // feasible cells: bi with gcd(bi,q-1)==w (bi=0 counts gcd=q-1)
  for(u64 bi=0;bi<q-1;bi++){ u64 g=gcd64(bi==0?q-1:bi,q-1); if(g==w) ncells+=q; }
  u64 nvalid=0, nbad=0;
  #pragma omp parallel for schedule(dynamic,65536) reduction(+:nvalid,nbad)
  for(u64 i=1;i<Q;i++){
    if(gcd64(i,Q)!=3) continue;
    nvalid++;
    u64 bi=i%(q-1);
    u64 s0=(Q-i)%Q; int acc[40]={0}; int d[40];
    for(int j=0;j<m;j++){ dec(EXP[(s0*qj[j])%Q],d); for(int t=0;t<N;t++) acc[t]+=d[t]; }
    for(int t=0;t<N;t++) acc[t]%=p;
    u64 tv=enc(acc); u64 aidx;
    if(tv==0) aidx=0; else { u64 lt=LOG[tv]; if(lt%step){ nbad++; continue; } aidx=1+((lt/step)+bi)%(q-1); }
    u64 cell=aidx*(q-1)+bi;
    if(__atomic_load_n(&cov[cell],__ATOMIC_RELAXED)) continue;
    int A[40][40]; int nv=0;
    for(int l=0;l<e;l++) for(int j=0;j<m;j++){ u64 ex=((u64)l*step+i*qj[j])%Q; dec(EXP[ex],A[nv]); nv++; }
    if(rank_vecs((int(*)[40])A,nv)==N-e){
      __atomic_store_n(&cov[cell],1,__ATOMIC_RELAXED);
    }
  }
  u64 covered=0; for(u64 c=0;c<cells;c++) covered+=cov[c];
  printf("q=%llu m=%d: valid_xi=%llu cells=%llu covered=%llu uncovered=%llu trace_not_in_Fq=%llu time=%.1fs\n",
    (unsigned long long)q,m,(unsigned long long)nvalid,(unsigned long long)ncells,(unsigned long long)covered,
    (unsigned long long)(ncells-covered),(unsigned long long)nbad,omp_get_wtime()-t0);
  if(covered!=ncells){
    int shown=0;
    for(u64 aidx=0;aidx<q&&shown<20;aidx++) for(u64 bi=0;bi<q-1&&shown<20;bi++){
      if(gcd64(bi==0?q-1:bi,q-1)!=w) continue;
      if(!cov[aidx*(q-1)+bi]){ printf("  UNCOVERED a_idx=%llu (0=zero) b_exp=%llu\n",(unsigned long long)aidx,(unsigned long long)bi); shown++; }
    }
  }
  return 0;
}
