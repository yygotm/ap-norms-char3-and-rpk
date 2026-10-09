// Independent port of verify.py: polynomial F_3 arithmetic and F_3 matrix rank.
// Packed coefficient planes accelerate arithmetic; no producer code was read.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
using U=uint32_t; using L=uint64_t;
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
void require(bool v,const std::string& s){if(!v) throw std::runtime_error(s);}
struct E {U a=0,b=0; bool operator==(E y)const{return a==y.a&&b==y.b;} bool operator!=(E y)const{return !(*this==y);} };
struct Field {
 int N; U mask; E f; E one{1,0},zero{0,0},x{2,0}; std::array<E,32> xp;
 Field(std::string p):N(p.size()),mask((U(1)<<N)-1){
  require(N>=2&&N<=24,"unsupported degree");
  for(int i=0;i<N;++i){char c=p[N-1-i];require(c>='0'&&c<='2',"bad polynomial digit");if(c=='1')f.a|=U(1)<<i;if(c=='2')f.b|=U(1)<<i;}
  xp[0]=x;for(int i=1;i<32;++i)xp[i]=mul(xp[i-1],xp[i-1]);
 }
 static E add(E u,E v){U uz=~(u.a|u.b),vz=~(v.a|v.b);return {(uz&v.a)|(u.a&vz)|(u.b&v.b),(uz&v.b)|(u.b&vz)|(u.a&v.a)};}
 static E neg(E v){return {v.b,v.a};}
 E times_x(E u)const {U t1=u.a>>(N-1),t2=u.b>>(N-1);u={(u.a<<1)&mask,(u.b<<1)&mask};if(t1)u=add(u,neg(f));else if(t2)u=add(u,f);return u;}
 E mul(E u,E v)const {E z;while(v.a|v.b){if(v.a&1)z=add(z,u);else if(v.b&1)z=add(z,neg(u));v.a>>=1;v.b>>=1;if(v.a|v.b)u=times_x(u);}return z;}
 E power(E u,L e)const {E z=one;while(e){if(e&1)z=mul(z,u);e>>=1;if(e)u=mul(u,u);}return z;}
 std::array<E,32> squares(E a,U limit)const{std::array<E,32> p{};p[0]=a;for(int i=1;(limit>>i)&&i<31;++i)p[i]=mul(p[i-1],p[i-1]);return p;}
 E from_squares(const std::array<E,32>& p,U e)const{E z=one;int i=0;while(e){if(e&1)z=mul(z,p[i]);e>>=1;++i;}return z;}
 E gpower(U e)const {E z=one;int k=0;while(e){if(e&1)z=mul(z,xp[k]);e>>=1;++k;}return z;}
 E dense(U v)const{E z;for(int i=0;i<N;++i){U d=v%3;v/=3;if(d==1)z.a|=U(1)<<i;else if(d==2)z.b|=U(1)<<i;}return z;}
 U dense(E z)const{U v=0,p=1;for(int i=0;i<N;++i){v+=p*(((z.a>>i)&1)+2*((z.b>>i)&1));p*=3;}return v;}
 bool normal(E z,U q,int k,int n,const std::vector<E>& subbasis)const {
  std::array<E,24> piv{};int rank=0;
  for(int i=0;i<n;++i){
   for(int ell=0;ell<k;++ell){E v=ell==0?z:mul(subbasis[ell],z);bool inserted=false;
    while(v.a|v.b){int p=31-__builtin_clz(v.a|v.b);if(!(piv[p].a|piv[p].b)){piv[p]=(v.b&(U(1)<<p))?neg(v):v;++rank;inserted=true;break;}
     v=add(v,(v.a&(U(1)<<p))?neg(piv[p]):piv[p]);}
    if(!inserted)return false;
   }
   if(i+1<n)z=power(z,q);
  }return rank==N;
 }
};
std::vector<U> prime_divisors(U v){std::vector<U> p;for(U d=2;L(d)*d<=v;d+=(d==2?1:2))if(v%d==0){p.push_back(d);while(v%d==0)v/=d;}if(v>1)p.push_back(v);return p;}
U read_u32(const char* p){return U(uint8_t(p[0]))|(U(uint8_t(p[1]))<<8)|(U(uint8_t(p[2]))<<16)|(U(uint8_t(p[3]))<<24);}
struct File {fs::path path;L rows,start;};
struct Job {
 U q,n,N,M,k;std::string poly;Field field;std::vector<U> primes,indices;std::vector<E> targets,basis;std::vector<File> files;L records=0,total=0,bytes=0;U need,rowbytes;
 // The cached value is a proved property of x^i, independent of beta/slot.
 // Zero means uncomputed; 2+2*norm_index+normal stores an exact result.
 mutable std::unique_ptr<std::atomic<uint8_t>[]> alpha_cache;
 static std::string getpoly(fs::path d){std::ifstream f(d/"header.txt");require(bool(f),"cannot open header");std::string line;while(std::getline(f,line))if(line.rfind("poly_high_to_low=",0)==0){auto s=line.substr(17);if(!s.empty()&&s.back()=='\r')s.pop_back();return s;}throw std::runtime_error("missing polynomial");}
 Job(fs::path dir):poly(getpoly(dir)),field(poly){
  std::ifstream in(dir/"header.txt");std::string line;std::getline(in,line);require(sscanf(line.c_str(),"q=%u n=%u N=%u M=%u",&q,&n,&N,&M)==4,"bad header");
  k=N/n;L size=1,qcheck=1;for(U i=0;i<N;++i)size*=3;for(U i=0;i<k;++i)qcheck*=3;
  require(N==poly.size()&&N==n*k&&qcheck==q&&M==size-1,"inconsistent header");
  primes=prime_divisors(M);require(field.power(field.x,M)==field.one,"x^M != 1");for(U p:primes)require(field.power(field.x,M/p)!=field.one,"x order too small");
  E gq=field.power(field.x,M/(q-1));for(U c=1;c<q;++c)if(std::gcd(c,q-1)==1){indices.push_back(c);targets.push_back(field.power(gq,c));}
  for(U ell=0;ell<k;++ell)basis.push_back(field.power(gq,ell));
  need=indices.size()*indices.size();rowbytes=4*(need+1);
  std::vector<fs::path> paths;for(auto& e:fs::directory_iterator(dir))if(e.is_regular_file()&&e.path().filename().string().rfind("cert_",0)==0&&e.path().extension()==".bin")paths.push_back(e.path());std::sort(paths.begin(),paths.end());
  require(!paths.empty(),"no certificates");for(auto& p:paths){L b=fs::file_size(p);require(b%rowbytes==0,"partial record: "+p.string());files.push_back({p,b/rowbytes,records});records+=b/rowbytes;bytes+=b;}total=records*need;
  alpha_cache=std::make_unique<std::atomic<uint8_t>[]>(M);
 }
 // Result codes: 0 pass, 1 exponent/range, 2 alpha primitive, 3 other primitive,
 // 4 first norm, 5 second norm, 6 neither normal.
 int check(U j,U exponent,U slot,E beta)const{
  if(j>=M||exponent>=M||slot>=need)return 1;
  E a=field.gpower(exponent),b=Field::add(a,beta);
  uint8_t cached=alpha_cache[exponent].load(std::memory_order_relaxed);
  if(!cached){
   if(std::gcd(exponent,M)!=1)return 2;
   auto ap=field.squares(a,M);
   for(U p:primes)if(field.from_squares(ap,M/p)==field.one)return 2;
   E norm=field.from_squares(ap,M/(q-1));
   auto pos=std::find(targets.begin(),targets.end(),norm);
   if(pos==targets.end())return 4;
   cached=2+2*std::distance(targets.begin(),pos)+field.normal(a,q,k,n,basis);
   alpha_cache[exponent].store(cached,std::memory_order_relaxed);
  }
  if(U((cached-2)/2)!=slot/indices.size())return 4;
  if(b==field.zero)return 3;
  auto bp=field.squares(b,M);
  for(U p:primes)if(field.from_squares(bp,M/p)==field.one)return 3;
  if(field.from_squares(bp,M/(q-1))!=targets[slot%indices.size()])return 5;
  if(!(cached&1)&&!field.normal(b,q,k,n,basis))return 6;
  return 0;
 }
};
struct Stats {L checked=0;std::array<L,7> failures{};};
void addstats(Stats& a,const Stats& b){a.checked+=b.checked;for(int i=1;i<7;++i)a.failures[i]+=b.failures[i];}
L failed(const Stats&s){return std::accumulate(s.failures.begin()+1,s.failures.end(),L(0));}
struct Witness{U j,i,slot;L index;};
std::vector<Witness> load_sample(const Job& job,const std::vector<L>& ids){
 std::vector<Witness> out;out.reserve(ids.size());size_t fi=0;std::ifstream f;L previous=~L(0);std::vector<char> row(job.rowbytes);
 for(L id:ids){L ri=id/job.need;while(fi+1<job.files.size()&&ri>=job.files[fi+1].start){++fi;f.close();previous=~L(0);}if(!f.is_open())f.open(job.files[fi].path,std::ios::binary);
  if(ri!=previous){f.seekg((ri-job.files[fi].start)*job.rowbytes);f.read(row.data(),job.rowbytes);require(L(f.gcount())==job.rowbytes,"short sample read");previous=ri;}
  U slot=id%job.need;out.push_back({read_u32(row.data()),read_u32(row.data()+4*(slot+1)),slot,id});}
 return out;
}
Stats check_sample(const Job& job,const std::vector<Witness>& ws){std::vector<Stats> local(omp_get_max_threads());
 #pragma omp parallel
 {auto& st=local[omp_get_thread_num()];
  #pragma omp for schedule(static)
  for(size_t p=0;p<ws.size();++p){auto w=ws[p];E beta=job.field.gpower(w.j);int code=job.check(w.j,w.i,w.slot,beta);++st.checked;if(code)++st.failures[code];}}
 Stats s;for(auto& t:local)addstats(s,t);return s;
}
void print_stats(const Stats&s){std::cout<<"\"checked\":"<<s.checked<<",\"failed\":"<<failed(s)<<",\"failure_codes\":[";for(int i=1;i<7;++i){if(i>1)std::cout<<',';std::cout<<s.failures[i];}std::cout<<']';}
void sampling(const Job&job,U count){
 count=std::min<L>(count,job.total/3);std::mt19937_64 rng(202610092056ULL);double maxrate=0;L all=0;
 for(std::string kind:{"head","tail","random"}){std::vector<L> ids;ids.reserve(count);
  if(kind=="head")for(U i=0;i<count;++i)ids.push_back(i);
  else if(kind=="tail")for(U i=0;i<count;++i)ids.push_back(job.total-count+i);
  else{std::uniform_int_distribution<L> d(count,job.total-count-1);while(ids.size()<count){while(ids.size()<count)ids.push_back(d(rng));std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());}}
  auto io=Clock::now();auto ws=load_sample(job,ids);double ios=seconds(io);auto t=Clock::now();Stats s=check_sample(job,ws);double sec=seconds(t);maxrate=std::max(maxrate,sec/s.checked);all+=s.checked;
  std::cout<<"{\"mode\":\"sample\",\"region\":\""<<kind<<"\",\"seed\":202610092056,\"q\":"<<job.q<<",\"n\":"<<job.n<<',';print_stats(s);std::cout<<",\"threads\":"<<omp_get_max_threads()<<",\"load_seconds\":"<<ios<<",\"check_seconds\":"<<sec<<",\"projected_full_seconds\":"<<sec/s.checked*job.total<<"}"<<std::endl;
 }
 std::cout<<"{\"mode\":\"estimate\",\"total\":"<<job.total<<",\"sample_checked\":"<<all<<",\"max_region_projection_seconds\":"<<maxrate*job.total<<"}"<<std::endl;
}
void cover(const Job&job){auto t=Clock::now();std::vector<uint8_t> bits((L(job.M)+7)/8,0);L marked=0,duplicates=0,invalid=0,readrows=0;
 auto get=[&](U x){return bool(bits[x>>3]&(1u<<(x&7)));};
 auto mark=[&](U x){U y=x;for(U k=0;k<job.N;++k){U neg=(L(y)+job.M/2)%job.M;for(U z:{y,neg}){if(!get(z)){bits[z>>3]|=1u<<(z&7);++marked;}}y=L(y)*3%job.M;}};
 std::vector<char> buffer(L(job.rowbytes)*4096);
 for(auto& file:job.files){std::ifstream f(file.path,std::ios::binary);require(bool(f),"cannot open coverage file");while(f){f.read(buffer.data(),buffer.size());L got=f.gcount();require(got%job.rowbytes==0,"coverage short row");for(L p=0;p<got;p+=job.rowbytes){++readrows;U j=read_u32(buffer.data()+p);if(j>=job.M){++invalid;continue;}if(get(j)){++duplicates;continue;}mark(j);}}}
 L uncovered=job.M-marked,missing=0;for(U j=0;j<job.M;++j)if(!get(j)){++missing;mark(j);}
 std::cout<<"{\"mode\":\"coverage\",\"q\":"<<job.q<<",\"n\":"<<job.n<<",\"records\":"<<readrows<<",\"invalid_j\":"<<invalid<<",\"duplicate_orbits\":"<<duplicates<<",\"uncovered_elements\":"<<uncovered<<",\"uncovered_orbits\":"<<missing<<",\"seconds\":"<<seconds(t)<<"}"<<std::endl;
 require(readrows==job.records&&invalid==0&&duplicates==0&&missing==0,"coverage FAILED");
}
void full(const Job&job){auto start=Clock::now(),last=start;Stats sum;L rows=0;std::atomic<L> printed{0};std::vector<char> buffer(L(job.rowbytes)*512);
 for(auto& file:job.files){std::ifstream f(file.path,std::ios::binary);require(bool(f),"cannot open file");while(f){f.read(buffer.data(),buffer.size());L got=f.gcount();require(got%job.rowbytes==0,"short record");L nr=got/job.rowbytes;std::vector<Stats> local(omp_get_max_threads());
   #pragma omp parallel
   {auto& st=local[omp_get_thread_num()];
    #pragma omp for schedule(static)
    for(L r=0;r<nr;++r){const char* rec=buffer.data()+r*job.rowbytes;U j=read_u32(rec);E beta=job.field.gpower(j);for(U slot=0;slot<job.need;++slot){U exponent=read_u32(rec+4*(slot+1));int code=job.check(j,exponent,slot,beta);++st.checked;if(code){++st.failures[code];
      #pragma omp critical
      {if(printed.fetch_add(1)<10)std::cerr<<"FAIL file="<<file.path.filename()<<" global_row="<<rows+r<<" j="<<j<<" slot="<<slot<<" i="<<exponent<<" code="<<code<<std::endl;}
     }}}}
   for(auto&s:local)addstats(sum,s);
   rows+=nr;if(seconds(last)>30){std::cout<<"{\"mode\":\"progress\",\"q\":"<<job.q<<",\"n\":"<<job.n<<",\"checked\":"<<sum.checked<<",\"total\":"<<job.total<<",\"seconds\":"<<seconds(start)<<"}"<<std::endl;last=Clock::now();}
  }}
 double checkseconds=seconds(start);std::cout<<"{\"mode\":\"full\",\"q\":"<<job.q<<",\"n\":"<<job.n<<",\"records\":"<<rows<<',';print_stats(sum);std::cout<<",\"seconds\":"<<checkseconds<<",\"threads\":"<<omp_get_max_threads()<<"}"<<std::endl;
 require(sum.checked==job.total&&rows==job.records,"count mismatch");cover(job);require(failed(sum)==0,"witness FAILED");std::cout<<"{\"mode\":\"complete\",\"q\":"<<job.q<<",\"n\":"<<job.n<<",\"seconds\":"<<seconds(start)<<"}"<<std::endl;
}
void fixtures(fs::path file){std::ifstream in(file);U N,q,n,a,b,e,ea,em,ep,en;std::string poly;L count=0;while(in>>N>>poly>>q>>n>>a>>b>>e>>ea>>em>>ep>>en){Field f(poly);require(f.N==int(N),"fixture N");U M=1;for(U i=0;i<N;++i)M*=3;--M;std::vector<E>basis;E gq=f.power(f.x,M/(q-1));for(U i=0;i<N/n;++i)basis.push_back(f.power(gq,i));require(f.dense(Field::add(f.dense(a),f.dense(b)))==ea,"fixture add");require(f.dense(f.mul(f.dense(a),f.dense(b)))==em,"fixture mul");require(f.dense(f.power(f.dense(a),e))==ep,"fixture pow");require(f.normal(f.dense(a),q,N/n,n,basis)==bool(en),"fixture normal");require(f.from_squares(f.squares(f.dense(a),M),e)==f.power(f.dense(a),e),"fixture shared squares");require(f.gpower(e)==f.power(f.x,e),"fixture gpower");++count;}require(count>0,"empty fixture");std::cout<<"{\"mode\":\"fixtures\",\"passed\":"<<count<<"}"<<std::endl;}
int main(int argc,char**argv){try{require(argc>=3,"usage: check MODE DIRECTORY [threads] [sample_count]");std::string mode=argv[1];if(mode=="fixtures"){fixtures(argv[2]);return 0;}int threads=argc>3?std::stoi(argv[3]):8;require(threads>=1&&threads<=8,"thread cap 8");omp_set_dynamic(0);omp_set_num_threads(threads);Job job(argv[2]);std::cout<<"{\"mode\":\"header\",\"q\":"<<job.q<<",\"n\":"<<job.n<<",\"records\":"<<job.records<<",\"witnesses\":"<<job.total<<",\"bytes\":"<<job.bytes<<",\"polynomial_primitive\":true,\"prime_divisors\":[";for(size_t i=0;i<job.primes.size();++i){if(i)std::cout<<',';std::cout<<job.primes[i];}std::cout<<"]}"<<std::endl;if(mode=="full")full(job);else if(mode=="sample")sampling(job,argc>4?std::stoul(argv[4]):30000);else if(mode=="coverage")cover(job);else throw std::runtime_error("unknown mode");return 0;}catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<std::endl;return 1;}}
