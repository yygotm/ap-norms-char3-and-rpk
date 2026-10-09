#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>
using namespace std;
using U=unsigned __int128;
using Clock=chrono::steady_clock;
double seconds(Clock::time_point t){return chrono::duration<double>(Clock::now()-t).count();}
string dec(U x){if(!x)return "0";string s;while(x){s+=char('0'+x%10);x/=10;}reverse(s.begin(),s.end());return s;}
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
U ipow(U a,int n){U r=1;while(n--)r*=a;return r;}
uint64_t modpow(uint64_t a,uint64_t n,uint64_t p){uint64_t r=1;while(n){if(n&1)r=U(r)*a%p;a=U(a)*a%p;n>>=1;}return r;}
bool prime(uint64_t n){if(n<2)return false;for(int p:{2,3,5,7,11,13,17,19,23,29,31,37})if(n%p==0)return n==uint64_t(p);uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;++s;}for(uint64_t a:{2ULL,325ULL,9375ULL,28178ULL,450775ULL,9780504ULL,1795265022ULL}){if(a%n==0)continue;uint64_t x=modpow(a%n,d,n);if(x==1||x==n-1)continue;bool ok=false;for(int i=1;i<s;++i){x=U(x)*x%n;if(x==n-1){ok=true;break;}}if(!ok)return false;}return true;}
vector<uint64_t> distinct(vector<uint64_t> a){sort(a.begin(),a.end());a.erase(unique(a.begin(),a.end()),a.end());return a;}
vector<uint64_t> smallf(uint64_t x){vector<uint64_t>a;for(uint64_t p=2;p*p<=x;++p)while(x%p==0){a.push_back(p);x/=p;}if(x>1)a.push_back(x);return a;}
struct Case{char kind;int q,m,p,e,n;uint64_t cells;U Q,M;vector<uint64_t> factors,mf;vector<int> red;string name()const{return to_string(q)+"_"+to_string(m);}};
vector<Case> cases(const string&path){ifstream f(path);need(bool(f),"cases open");vector<Case>out;string line;while(getline(f,line)){istringstream s(line);Case c{};int nf;s>>c.kind>>c.q>>c.m>>c.p>>c.e>>c.cells>>nf;c.n=c.e*c.m;c.Q=ipow(c.q,c.m)-1;c.M=c.Q/3;U prod=1;for(int i=0;i<nf;++i){uint64_t p;s>>p;need(prime(p),"factor not prime");prod*=p;c.factors.push_back(p);}need(prod==c.Q&&c.Q%3==0,"factor product");c.mf=c.factors;auto it=find(c.mf.begin(),c.mf.end(),3);need(it!=c.mf.end(),"missing 3");c.mf.erase(it);c.mf=distinct(c.mf);need(prime(c.p)&&ipow(c.p,c.e)==U(c.q)&&c.n<64,"field parameters");if(c.kind=='B'){for(int j=0;j<c.n;++j){int v=-1;s>>v;need(v>=0&&v<c.p,"red range");c.red.push_back(v);}}need(bool(s),"cases parse");out.push_back(c);}return out;}
struct E{uint64_t b=0;array<uint16_t,64> v{};};
struct Field{
 int p,n;bool binary;vector<int> red,inv;uint64_t br=0,mask=0;vector<E> fm;
 Field(int p_,int n_,const vector<int>&r):p(p_),n(n_),binary(p==2),red(r),inv(p){for(int i=1;i<p;++i)inv[i]=modpow(i,p-2,p);if(binary){for(int i=0;i<n;++i)br|=uint64_t(red[i])<<i;mask=(uint64_t(1)<<n)-1;}}
 E one()const{E a;if(binary)a.b=1;else a.v[0]=1;return a;}
 E mon(int i)const{E a;if(binary)a.b=uint64_t(1)<<i;else a.v[i]=1;return a;}
 bool zero(const E&a)const{if(binary)return a.b==0;for(int i=0;i<n;++i)if(a.v[i])return false;return true;}
 bool eq(const E&a,const E&b)const{if(binary)return a.b==b.b;for(int i=0;i<n;++i)if(a.v[i]!=b.v[i])return false;return true;}
 E add(E a,const E&b)const{if(binary)a.b^=b.b;else for(int i=0;i<n;++i)a.v[i]=(a.v[i]+b.v[i])%p;return a;}
 E sub(E a,const E&b)const{if(binary)a.b^=b.b;else for(int i=0;i<n;++i)a.v[i]=(a.v[i]+p-b.v[i])%p;return a;}
 E mul(const E&a,const E&b)const{
  E z;if(binary){uint64_t x=a.b,y=b.b,r=0;while(y){if(y&1)r^=x;y>>=1;bool top=x>>(n-1);x=(x<<1)&mask;if(top)x^=br;}z.b=r;return z;}
  array<uint64_t,128>t{};for(int i=0;i<n;++i)if(a.v[i])for(int j=0;j<n;++j)t[i+j]+=uint64_t(a.v[i])*b.v[j];
  for(int i=2*n-2;i>=n;--i){uint64_t d=t[i]%p;if(d)for(int j=0;j<n;++j)if(red[j])t[i-n+j]+=d*red[j];}
  for(int i=0;i<n;++i)z.v[i]=t[i]%p;
  return z;
 }
 E power(E a,U k)const{E r=one();while(k){if(k&1)r=mul(r,a);k>>=1;if(k)a=mul(a,a);}return r;}
 E from(uint64_t x)const{E a;if(binary)a.b=x;else for(int i=0;i<n;++i){a.v[i]=x%p;x/=p;}return a;}
 uint32_t packed(const E&a)const{if(binary)return a.b;uint32_t r=0;for(int i=n-1;i>=0;--i)r=r*p+a.v[i];return r;}
 string key(const E&a)const{if(binary)return string(reinterpret_cast<const char*>(&a.b),8);return string(reinterpret_cast<const char*>(a.v.data()),2*n);}
 vector<int> polynomial(const E&a)const{vector<int>v(n);for(int i=0;i<n;++i)v[i]=binary?((a.b>>i)&1):a.v[i];return v;}
 void frobenius(int q){fm.clear();for(int i=0;i<n;++i)fm.push_back(power(mon(i),q));}
 E frob(const E&a)const{E z;if(binary){uint64_t x=a.b;while(x){int j=countr_zero(x);z.b^=fm[j].b;x&=x-1;}}else{array<uint64_t,64>v{};for(int i=0;i<n;++i)if(a.v[i])for(int j=0;j<n;++j)v[j]+=uint64_t(a.v[i])*fm[i].v[j];for(int j=0;j<n;++j)z.v[j]=v[j]%p;}return z;}
 int rank(const vector<E>&cols,int stop=1000)const{
  int r=0;if(binary){array<uint64_t,64>basis{};for(auto&v:cols){uint64_t x=v.b;while(x){int k=63-countl_zero(x);if(basis[k])x^=basis[k];else{basis[k]=x;++r;break;}}if(r>=stop)break;}return r;}
  array<array<int,64>,64>basis{};array<bool,64>used{};for(auto&v:cols){array<int,64>x{};for(int i=0;i<n;++i)x[i]=v.v[i];for(int k=0;k<n;++k)if(x[k]){if(used[k]){int d=x[k];for(int j=k;j<n;++j)x[j]=(x[j]+p-d*basis[k][j]%p)%p;}else{int d=inv[x[k]];for(int j=k;j<n;++j)basis[k][j]=x[j]*d%p;used[k]=true;++r;break;}}if(r>=stop)break;}return r;
 }
};
void trim(vector<int>&a){while(!a.empty()&&a.back()==0)a.pop_back();}
vector<int> rempoly(vector<int>a,const vector<int>&b,int p){need(!b.empty(),"poly division zero");trim(a);int inv=modpow(b.back(),p-2,p);while(a.size()>=b.size()){int d=a.back()*inv%p,k=a.size()-b.size();for(size_t j=0;j<b.size();++j)a[k+j]=(a[k+j]+p-d*b[j]%p)%p;trim(a);}return a;}
int gcddeg(vector<int>a,vector<int>b,int p){trim(a);trim(b);while(!b.empty()){auto r=rempoly(a,b,p);a=b;b=r;}return int(a.size())-1;}
bool irreducible(const Field&f){vector<int>mod(f.n+1);mod[f.n]=1;for(int i=0;i<f.n;++i)mod[i]=(f.p-f.red[i])%f.p;auto fac=distinct(smallf(f.n));set<int>check;for(auto l:fac)check.insert(f.n/l);E x=f.mon(1),y=x;for(int j=1;j<=f.n;++j){y=f.power(y,f.p);if(check.count(j)&&gcddeg(mod,f.polynomial(f.sub(y,x)),f.p)!=0)return false;}return f.eq(y,x);}
bool exact(const Field&f,const E&x,U order,const vector<uint64_t>&fac){if(f.zero(x)||!f.eq(f.power(x,order),f.one()))return false;for(auto l:fac)if(order%l==0&&f.eq(f.power(x,order/l),f.one()))return false;return true;}
vector<int> choose_mod(const Case&c){mt19937 gen(320000+c.q*100+c.m);for(int attempt=0;attempt<100000;++attempt){vector<int>r(c.n);for(auto&x:r)x=gen()%c.p;r[0]=1+gen()%(c.p-1);Field f(c.p,c.n,r);if(irreducible(f)&&exact(f,f.mon(1),c.Q,distinct(c.factors)))return r;}throw runtime_error("modulus search exhausted");}
struct Subfield{vector<E> all,basis;unordered_map<string,int> ids;vector<int> allowed;int count=0;};
Subfield subfield(Field&f,const Case&c){Subfield s;E t;auto ff=distinct(smallf(c.q-1));for(uint64_t i=2;;++i){t=f.power(f.from(i),c.Q/(c.q-1));if(exact(f,t,c.q-1,ff))break;need(i<100000,"subfield generator search");}s.all.push_back(E{});E v=f.one();for(int j=0;j<c.q-1;++j){s.all.push_back(v);v=f.mul(v,t);}need(f.eq(v,f.one()),"subfield cycle");for(int j=0;j<c.q;++j){need(f.eq(f.power(s.all[j],c.q),s.all[j]),"subfield fixed");need(s.ids.emplace(f.key(s.all[j]),j).second,"subfield duplicate");}v=f.one();for(int j=0;j<c.e;++j){s.basis.push_back(v);v=f.mul(v,t);}need(f.rank(s.basis)==c.e,"base field basis rank");int w=gcd(3,c.q-1),ord=(c.q-1)/w;s.allowed.assign(c.q,-1);for(int j=1;j<c.q;++j)if((c.q-1)/gcd(j-1,c.q-1)==ord)s.allowed[j]=s.count++;need(uint64_t(c.q)*s.count==c.cells,"expected cells vs header");return s;}
int normalrank(const Field&f,const Case&c,const E&x,const Subfield&s){vector<E>cols;cols.reserve(c.n);for(auto&v:s.basis){E y=f.mul(x,v);for(int j=0;j<c.m;++j){cols.push_back(y);y=f.frob(y);}}return f.rank(cols,c.e*(c.m-1)+1);}
pair<E,E> cell(const Field&f,const Case&c,const E&x){E b=f.one(),a,y=x;for(int j=0;j<c.m;++j){a=f.add(f.mul(a,y),b);b=f.mul(b,y);y=f.frob(y);}need(f.eq(y,x),"Frobenius cycle");return {a,b};}
vector<E> loadbin(const Field&f,const Case&c,const string&dir){ifstream in(dir+"/"+c.name()+".bin",ios::binary);need(bool(in),"binary open");in.seekg(0,ios::end);need(uint64_t(in.tellg())==2ULL*c.n*c.cells,"binary size");in.seekg(0);vector<E>v(c.cells);for(auto&x:v)for(int j=0;j<c.n;++j){unsigned char z[2];in.read(reinterpret_cast<char*>(z),2);need(bool(in),"binary truncation");int d=z[0]+256*z[1];need(d<f.p,"digit out of range");if(f.binary)x.b|=uint64_t(d)<<j;else x.v[j]=d;}return v;}
void runB(const Case&c,const string&dir,bool bench){auto start=Clock::now();Field f(c.p,c.n,c.red);need(irreducible(f),"reducible modulus "+c.name());f.frobenius(c.q);Subfield sf=subfield(f,c);auto v=loadbin(f,c,dir);uint64_t count=bench?min<uint64_t>(c.cells,3000):c.cells;vector<int>seen(c.cells);uint64_t bad=0,ordbad=0,rankbad=0,cellbad=0;double setup=seconds(start);auto loop=Clock::now();
 #pragma omp parallel for schedule(dynamic,32) reduction(+:bad,ordbad,rankbad,cellbad)
 for(uint64_t j=0;j<count;++j){uint64_t i=bench?j*c.cells/count:j;bool ok=true;if(!exact(f,v[i],c.M,c.mf)){++ordbad;ok=false;}if(normalrank(f,c,v[i],sf)!=c.e*(c.m-1)){++rankbad;ok=false;}auto [a,b]=cell(f,c,v[i]);auto ia=sf.ids.find(f.key(a)),ib=sf.ids.find(f.key(b));if(ia==sf.ids.end()||ib==sf.ids.end()||sf.allowed[ib->second]<0){++cellbad;ok=false;}else{int k=ia->second*sf.count+sf.allowed[ib->second];
 #pragma omp atomic update
 ++seen[k];}
 if(!ok){++bad;
 #pragma omp critical
 {if(bad<=10)cerr<<"FAIL B "<<c.name()<<" record="<<i<<" order="<<ordbad<<" rank="<<rankbad<<" cell="<<cellbad<<'\n';}}
 }
 uint64_t missing=0,duplicates=0;for(int x:seen){if(!x)++missing;if(x>1)duplicates+=x-1;}double elapsed=seconds(loop);cout<<"RESULT B "<<c.name()<<" mode="<<(bench?"bench":"full")<<" total="<<c.cells<<" checked="<<count<<" bad="<<bad<<" order_bad="<<ordbad<<" rank_bad="<<rankbad<<" cell_bad="<<cellbad<<" missing="<<missing<<" duplicates="<<duplicates<<" setup_s="<<setup<<" loop_s="<<elapsed<<" total_s="<<seconds(start)<<" estimate_s="<<setup+elapsed*c.cells/count<<endl;need(bad==0,"B witness failure");if(!bench){need(missing==0&&duplicates==0,"B coverage failure");for(uint64_t i:array<uint64_t,3>{0,c.cells/2,c.cells-1}){auto[a,b]=cell(f,c,v[i]);cout<<"FIXTURE "<<c.name()<<" index="<<i<<" rank="<<normalrank(f,c,v[i],sf)<<" a=";for(int d:f.polynomial(a))cout<<d<<',';cout<<" b=";for(int d:f.polynomial(b))cout<<d<<',';cout<<endl;}}}
struct Powers{vector<uint32_t>v;vector<uint64_t>pp;Powers(const Field&f,uint64_t Q){pp.resize(f.n+1,1);for(int i=1;i<=f.n;++i)pp[i]=pp[i-1]*f.p;v.resize(Q);v[0]=1;if(f.binary){uint32_t x=1;for(uint64_t i=1;i<Q;++i){bool top=x>>(f.n-1);x=(x<<1)&f.mask;if(top)x^=f.br;v[i]=x;}E last=f.mul(f.from(x),f.mon(1));need(f.eq(last,f.one()),"power cycle binary");}else{int k=1,base=f.p;while(base*f.p<=8192){base*=f.p;++k;}int chunks=(f.n+k-1)/k;vector<vector<vector<uint32_t>>>ad(f.p,vector<vector<uint32_t>>(chunks,vector<uint32_t>(base)));vector<uint32_t>chunkpow(chunks,1);for(int h=1;h<chunks;++h)chunkpow[h]=chunkpow[h-1]*base;for(int top=0;top<f.p;++top)for(int h=0;h<chunks;++h)for(int x=0;x<base;++x){int t=x;uint32_t r=0,z=1;for(int j=0;j<k&&h*k+j<f.n;++j){r+=((t%f.p+top*f.red[h*k+j])%f.p)*z;z*=f.p;t/=f.p;}ad[top][h][x]=r;}uint32_t x=1;for(uint64_t i=1;i<Q;++i){int top=x/pp[f.n-1];uint32_t y=(x%pp[f.n-1])*f.p,r=0;for(int h=0;h<chunks;++h){r+=ad[top][h][y%base]*chunkpow[h];y/=base;}x=r;v[i]=x;}need(f.eq(f.mul(f.from(x),f.mon(1)),f.one()),"power cycle odd");}}};
int rankA(const Field&f,const Case&c,const Powers&ps,uint64_t t){uint64_t Q=uint64_t(c.Q),ex=3*t,shift=Q/(c.q-1);int rank=0,target=c.e*(c.m-1);if(f.binary){array<uint32_t,64>basis{};for(int j=0;j<c.m;++j){for(int h=0;h<c.e;++h){uint32_t x=ps.v[(ex+h*shift)%Q];while(x){int k=31-countl_zero(x);if(basis[k])x^=basis[k];else{basis[k]=x;++rank;break;}}if(rank>target)return rank;}ex=ex*c.q%Q;}return rank;}array<array<int,32>,32>b{};array<bool,32>used{};for(int j=0;j<c.m;++j){uint32_t code=ps.v[ex];array<int,32>x{};for(int i=0;i<c.n;++i){x[i]=code%c.p;code/=c.p;}for(int k=0;k<c.n;++k)if(x[k]){if(used[k]){int d=x[k];for(int h=k;h<c.n;++h)x[h]=(x[h]+c.p-d*b[k][h]%c.p)%c.p;}else{int d=f.inv[x[k]];for(int h=k;h<c.n;++h)b[k][h]=x[h]*d%c.p;used[k]=true;++rank;break;}}if(rank>target)return rank;ex=ex*c.q%Q;}return rank;}
uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
void runA(Case c,bool bench){auto start=Clock::now();c.red=choose_mod(c);Field f(c.p,c.n,c.red);f.frobenius(c.q);Subfield sf=subfield(f,c);cout<<"MOD A "<<c.name()<<" red=";for(int x:c.red)cout<<x<<',';cout<<endl;Powers ps(f,uint64_t(c.Q));uint64_t M=uint64_t(c.M),Q=uint64_t(c.Q);uint64_t count=bench?min<uint64_t>(M,300000):M;double setup=seconds(start);auto loop=Clock::now();uint64_t eligible=0,reps=0,normal=0;vector<uint64_t>counts(c.cells);vector<uint64_t>traces(c.q);uint64_t orbitbad=0;
 #pragma omp parallel
 {vector<uint64_t>local(c.cells),lt(c.q);uint64_t le=0,lr=0,ln=0,lo=0;
 #pragma omp for schedule(dynamic,4096)
 for(uint64_t j=0;j<count;++j){uint64_t t=j;if(bench&&count<M){uint64_t third=count/3;t=j<third?j:(j<2*third?M-1-(j-third):mix(j)%M);}bool coprime=true;for(auto p:c.mf)if(t%p==0){coprime=false;break;}if(!coprime)continue;++le;uint64_t u=t;bool canonical=true;int len=0;do{u=u*c.q%M;++len;if(u<t){canonical=false;break;}}while(u!=t&&len<=c.m);if(!canonical)continue;if(u!=t||len!=c.m){++lo;continue;}++lr;if(rankA(f,c,ps,t)!=c.e*(c.m-1))continue;ln+=len;uint64_t ex=3*t;E tr;for(int k=0;k<c.m;++k){tr=f.add(tr,f.from(ps.v[(Q-ex)%Q]));ex=ex*c.q%Q;}E b=f.from(ps.v[(3*t%(c.q-1))*(Q/(c.q-1))]);E a=f.mul(tr,b);auto ia=sf.ids.find(f.key(a)),ib=sf.ids.find(f.key(b)),it=sf.ids.find(f.key(tr));if(ia==sf.ids.end()||ib==sf.ids.end()||it==sf.ids.end()||sf.allowed[ib->second]<0){++lo;continue;}local[ia->second*sf.count+sf.allowed[ib->second]]+=len;lt[it->second]+=len;
 }
 #pragma omp critical
 {eligible+=le;reps+=lr;normal+=ln;orbitbad+=lo;for(size_t i=0;i<counts.size();++i)counts[i]+=local[i];for(int i=0;i<c.q;++i)traces[i]+=lt[i];}
 }
 double elapsed=seconds(loop);uint64_t missing=count_if(counts.begin(),counts.end(),[](auto x){return !x;});cout<<"RESULT A "<<c.name()<<" mode="<<(bench?"bench":"full")<<" field="<<Q+1<<" potential="<<M<<" scanned="<<count<<" primitive="<<eligible<<" orbits="<<reps<<" normal="<<normal<<" cells="<<c.cells<<" missing="<<missing<<" bad="<<orbitbad<<" setup_s="<<setup<<" loop_s="<<elapsed<<" total_s="<<seconds(start)<<" estimate_s="<<setup+elapsed*M/count<<endl;need(orbitbad==0,"A orbit or cell inconsistency");if(!bench){uint64_t phi=M;for(auto p:c.mf)phi=phi/p*(p-1);need(eligible==phi&&reps*c.m==eligible,"A orbit enumeration coverage");cout<<"COUNTS A "<<c.name()<<" cells=";for(auto x:counts)cout<<x<<',';cout<<" traces=";for(auto x:traces)cout<<x<<',';cout<<endl;for(int ia=0;ia<c.q;++ia)for(int ib=1;ib<c.q;++ib)if(sf.allowed[ib]>=0&&counts[ia*sf.count+sf.allowed[ib]]==0){cout<<"MISSING A "<<c.name()<<" a_id="<<ia<<" b_id="<<ib<<endl;}need(missing==(c.q==2&&c.m==8?1:0),"A claim mismatch");}}
int main(int argc,char**argv){try{need(argc>=5,"usage check cases data mode [pair]");omp_set_dynamic(0);omp_set_num_threads(8);cout.setf(ios::fixed);cout.precision(6);string mode=argv[3],pair=argv[4];auto cs=cases(argv[1]);auto start=Clock::now();cout<<"START mode="<<mode<<" pair="<<pair<<" threads=8"<<endl;for(auto&c:cs)if((pair=="all"||pair==string(1,c.kind)||pair==c.name())){if(c.kind=='A')runA(c,mode=="bench");else runB(c,argv[2],mode=="bench");}cout<<"COMPLETE mode="<<mode<<" pair="<<pair<<" elapsed_s="<<seconds(start)<<endl;return 0;}catch(const exception&e){cerr<<"ERROR "<<e.what()<<endl;return 1;}}
