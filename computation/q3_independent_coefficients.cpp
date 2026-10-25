// Independent six-variable sparse integer reconstruction of the clean q=3 pair certificate.
// Original encoder is included ONLY to compare all final rational coefficients.
#define main original_q3_main
#include "../analytic/q3_exact_pair_simple.cpp"
#undef main
#include <unordered_map>
#include <fstream>
#include <cstdint>
namespace independent {
using Key=std::uint32_t;
using Big=boost::multiprecision::cpp_int;
struct Poly { int dl=0,dr=0; std::unordered_map<Key,Big> c; };
Poly constant(long long v){Poly p;if(v)p.c[0]=v;return p;}
Poly variable(int j){Poly p;p.dl=j<3;p.dr=j>=3;p.c[Key(1)<<(5*j)]=1;return p;}
Poly clean(Poly p){for(auto it=p.c.begin();it!=p.c.end();)if(it->second==0)it=p.c.erase(it);else ++it;return p;}
Poly times(Poly p,long long k){for(auto&v:p.c)v.second*=k;return clean(std::move(p));}
Poly plus(Poly a,const Poly&b){
 if(a.c.empty())return b;if(b.c.empty())return a;
 assert(a.dl==b.dl&&a.dr==b.dr);
 for(const auto&v:b.c)a.c[v.first]+=v.second;return clean(std::move(a));
}
Poly minus(Poly a,const Poly&b){return plus(std::move(a),times(b,-1));}
Poly product(const Poly&a,const Poly&b){
 Poly r;r.dl=a.dl+b.dl;r.dr=a.dr+b.dr;
 // Packed addition is exact: each exponent is <= its total degree, which stays below 32.
 assert(r.dl<32&&r.dr<32);
 r.c.reserve(std::min<std::size_t>(100000,a.c.size()*b.c.size()));
 for(const auto&u:a.c)for(const auto&v:b.c)r.c[u.first+v.first]+=u.second*v.second;
 return clean(std::move(r));
}
Poly pow(Poly a,int n){Poly p=constant(1);for(int j=0;j<n;j++)p=product(p,a);return p;}
struct Forms { Poly s,x,y,h,f,g;std::array<Poly,3> probability_numerators; };
Forms forms(int offset){
 Poly a=variable(offset),b=variable(offset+1),c=variable(offset+2);
 Poly s=plus(plus(a,b),c);
 Poly x=plus(plus(times(pow(a,2),4),times(product(a,b),2)),pow(b,2));
 Poly y=minus(minus(plus(times(pow(a,3),8),times(product(pow(a,2),b),6)),
                    times(product(a,pow(b,2)),3)),pow(b,3));
 Poly h=minus(times(product(y,s),4),pow(x,2));
 Poly f=minus(times(product(x,s),3),y);
 Poly g=minus(times(pow(s,2),4),x);
 std::array<Poly,3>p={plus(plus(times(a,6),times(b,3)),times(c,2)),plus(times(b,3),times(c,2)),times(c,2)};
 return{s,x,y,h,f,g,p};
}
Poly target(int mode){
 auto l=forms(0),r=forms(3);
 Poly s2=pow(l.s,2),t2=pow(r.s,2),s3=pow(l.s,3),t3=pow(r.s,3),s4=pow(l.s,4),t4=pow(r.s,4);
 Poly strict=product(product(l.x,r.x),product(s2,t2));
 Poly k;
 if(mode==0){
  k=plus(times(product(product(l.x,s2),t4),8000),times(product(product(r.x,t2),s4),8000));
  k=minus(k,times(product(product(l.h,r.f),r.s),125));
  k=minus(k,times(product(product(l.f,r.h),l.s),125));
  k=minus(k,times(strict,2));
 }else{
  k=minus(times(product(product(l.h,r.g),t2),-125),times(product(product(l.g,r.h),s2),125));
  k=minus(k,strict);
 }
 std::vector<Poly>normalizers,numerators,factors;
 std::array<int,3> permutation={0,1,2};
 do {
  std::array<Poly,3>u={product(l.probability_numerators[0],r.probability_numerators[permutation[0]]),
                      product(l.probability_numerators[1],r.probability_numerators[permutation[1]]),
                      product(l.probability_numerators[2],r.probability_numerators[permutation[2]])};
  Poly z=plus(plus(u[0],u[1]),u[2]);
  std::array<Poly,3>w={minus(times(u[0],3),z),minus(times(u[1],3),z),minus(times(u[2],3),z)};
  Poly v2=plus(plus(pow(w[0],2),pow(w[1],2)),pow(w[2],2));
  Poly n=mode?minus(times(product(product(product(w[0],w[1]),w[2]),z),36),pow(v2,2)):v2;
  normalizers.push_back(z);numerators.push_back(n);factors.push_back(pow(z,mode?3:1));
 }while(std::next_permutation(permutation.begin(),permutation.end()));
 Poly denominator=constant(1);
 for(const auto&f:factors)denominator=product(denominator,f);
 Poly result=times(product(k,denominator),mode?81:27);
 Poly sum;bool first=true;
 for(int j=0;j<6;j++){
  Poly rest=constant(1);
  for(int i=0;i<6;i++)if(i!=j)rest=product(rest,factors[i]);
  Poly term=product(numerators[j],rest);
  if(first){sum=std::move(term);first=false;}else sum=plus(std::move(sum),term);
 }
 sum=product(product(sum,s3),t3);
 result=plus(std::move(result),times(std::move(sum),mode?250:-2000));
 assert(result.dl==(mode?22:10)&&result.dr==result.dl);
 return result;
}
Key pack(const E&a,const E&b){Key k=0;for(int j=0;j<3;j++){k|=Key(a[j])<<(5*j);k|=Key(b[j])<<(5*(j+3));}return k;}
Big integer_power(Big x,int n){Big r=1;while(n--)r*=x;return r;}
void check(int mode){
 auto start=std::chrono::steady_clock::now();
 Poly mine=target(mode);
 P original=make_target(mode);
 Big factor=(mode?Big(216000):Big(144000))*integer_power(Big(12),mode?18:6);
 std::size_t n=0,zero=0,pos=0,negative=0;
 std::ofstream dump(mode?"q3_independent_coefficients_B.tsv":"q3_independent_coefficients_A.tsv");
 dump<<"a\tb\tc\td\te\tf\tinteger_coefficient\n";
 for(int i=0;i<S(original.d);i++)for(int j=0;j<S(original.d);j++){
  Key k=pack(mons[original.d][i],mons[original.d][j]);
  auto it=mine.c.find(k);Big value=it==mine.c.end()?Big(0):it->second;
  assert(value*original.den==factor*original.at(i,j));
  n++;zero+=value==0;pos+=value>0;negative+=value<0;
  for(int v:mons[original.d][i])dump<<v<<"\t";
  for(int v:mons[original.d][j])dump<<v<<"\t";
  dump<<value<<"\n";
 }
 assert(n==(std::size_t)S(mine.dl)*S(mine.dr));
 assert(negative==0);
 assert(mine.c.size()==pos);
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"INDEPENDENT_PASS mode="<<mode<<" degree="<<mine.dl<<" coefficients="<<n<<" positive="<<pos<<" zero="<<zero<<" negative="<<negative<<" all_original_coefficients_identical_after_exact_scaling=true seconds="<<seconds<<std::endl;
}
}
int main(int argc,char**argv){init();int mode=argc>1?std::stoi(argv[1]):0;independent::check(mode);}

