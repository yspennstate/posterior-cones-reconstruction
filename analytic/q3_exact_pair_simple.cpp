// Exact candidate certificate for q=3: H=y-x^2/2, F=3x/2-y, G=1-x.
// Coefficients and all arithmetic are rational / cpp_int.
// A complete successful covering is needed for a global proof.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using Z=boost::multiprecision::cpp_int;
using E=array<int,3>;
vector<E> mons[32];
int ix[32][32][32];
Z gcdZ(Z a,Z b){if(a<0)a=-a;if(b<0)b=-b;while(b!=0){Z c=a%b;a=b;b=c;}return a;}
Z factorial(int k){Z a=1;for(int i=2;i<=k;i++)a*=i;return a;}
Z multinom(E a){int d=a[0]+a[1]+a[2];return factorial(d)/(factorial(a[0])*factorial(a[1])*factorial(a[2]));}
void init(){for(int d=0;d<32;d++)for(int i=d;i>=0;i--)for(int j=d-i;j>=0;j--){ix[d][i][j]=mons[d].size();mons[d].push_back({i,j,d-i-j});}}
int S(int d){return mons[d].size();}
struct P{int d;Z den=1;vector<Z>a;explicit P(int n=0):d(n),a(S(n)*S(n)){}Z& at(int i,int j){return a[i*S(d)+j];}const Z& at(int i,int j)const{return a[i*S(d)+j];}};
P one(int d){P p(d);for(int i=0;i<S(d);i++)for(int j=0;j<S(d);j++)p.at(i,j)=multinom(mons[d][i])*multinom(mons[d][j]);return p;}
P scale(P a,Z k){for(auto&v:a.a)v*=k;return a;}
P divide(P a,int k){assert(k>0);a.den*=k;return a;}
P mul(const P&a,const P&b){P c(a.d+b.d);c.den=a.den*b.den;vector<int>w(S(a.d)*S(b.d));for(int i=0;i<S(a.d);i++)for(int j=0;j<S(b.d);j++){auto u=mons[a.d][i],v=mons[b.d][j];w[i*S(b.d)+j]=ix[c.d][u[0]+v[0]][u[1]+v[1]];}for(int j=0;j<S(b.d);j++)for(int k=0;k<S(b.d);k++){const Z&B=b.at(j,k);if(B==0)continue;for(int i=0;i<S(a.d);i++){Z*dest=&c.at(w[i*S(b.d)+j],0);for(int h=0;h<S(a.d);h++)if(a.at(i,h)!=0)dest[w[h*S(b.d)+k]]+=a.at(i,h)*B;}}return c;}
P power(P a,int n){P r=one(0);while(n--){r=mul(r,a);}return r;}
P lift(P a,int d){assert(a.d<=d);if(a.d<d)a=mul(a,one(d-a.d));return a;}
P add(P a,P b){int d=max(a.d,b.d);a=lift(a,d);b=lift(b,d);Z g=gcdZ(a.den,b.den),u=b.den/g,v=a.den/g;for(size_t i=0;i<a.a.size();i++)a.a[i]=a.a[i]*u+b.a[i]*v;a.den*=u;return a;}
P sub(P a,P b){return add(a,scale(b,-1));}
void reduce(P&a){Z g=a.den;for(auto&v:a.a){g=gcdZ(g,v);if(g==1)return;}if(g!=1){for(auto&v:a.a)v/=g;a.den/=g;}}
P tensor(const P&a,const P&b){int d=max(a.d,b.d);P x=lift(a,d),y=lift(b,d),r(d);r.den=x.den*y.den;for(int i=0;i<S(d);i++)for(int j=0;j<S(d);j++)r.at(i,j)=x.at(i,0)*y.at(0,j);return r;}
struct Inv{P x,y,h,f,g;};
Inv invariants(const array<P,3>&p,const P&z){vector<P>v;for(int i=0;i<3;i++)v.push_back(sub(scale(p[i],3),z));P x=one(0),y=one(0);x=scale(x,0);y=scale(y,0);for(auto&a:v){x=add(x,mul(a,a));y=add(y,power(a,3));}x=divide(x,6);y=divide(y,6);P h=sub(mul(y,z),divide(mul(x,x),2));P f=divide(sub(scale(mul(x,z),3),scale(y,2)),2);P g=sub(mul(z,z),x);return{x,y,h,f,g};}

P make_target(int mode){
  int vertex[3][3]={{6,0,0},{3,3,0},{2,2,2}};
  array<P,3>p={P(1),P(1),P(1)},q=p;
  for(int c=0;c<3;c++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){p[c].at(i,j)=vertex[i][c];q[c].at(i,j)=vertex[j][c];}
  for(auto*r:{&p,&q})for(auto&a:*r)a.den=6;
  P unit=one(1);auto a=invariants(p,unit),b=invariants(q,unit);
  P k=scale(add(tensor(a.h,mode?b.g:b.f),tensor(mode?a.g:a.f,b.h)),-1);
  if(!mode)k=add(k,add(a.x,b.x));
  k=sub(k,divide(tensor(a.x,b.x),1000));
  vector<P>zs,ns;array<int,3>perm={0,1,2};
  do{
    array<P,3>u={P(1),P(1),P(1)};
    for(int c=0;c<3;c++){for(int i=0;i<3;i++)for(int j=0;j<3;j++)u[c].at(i,j)=3*vertex[i][c]*vertex[j][perm[c]];u[c].den=36;}
    P z=add(add(u[0],u[1]),u[2]);auto t=invariants(u,z);zs.push_back(z);ns.push_back(mode?t.h:scale(t.x,-1));
  }while(next_permutation(perm.begin(),perm.end()));
  vector<P>factors;for(auto&z:zs)factors.push_back(power(z,mode?3:1));
  P den=one(0);for(auto&f:factors){den=mul(den,f);reduce(den);}
  P out=scale(mul(k,den),6);reduce(out);
  for(int j=0;j<6;j++){P rest=one(0);for(int l=0;l<6;l++)if(l!=j){rest=mul(rest,factors[l]);reduce(rest);}out=add(out,mul(ns[j],rest));reduce(out);}
  cerr<<"target mode="<<mode<<" degree="<<out.d<<" coefficients="<<out.a.size()<<"\n";return out;
}

// Midpoint bisection: old Xi=new Xi/2, old Xj=new Xj+new Xi/2.
// Returns one exact child; exchanging i,j gives the other child.
P bisect(const P&a,bool right,int i,int j){P r(a.d);r.den=a.den*(Z(1)<<a.d);for(int u=0;u<S(a.d);u++){E e=mons[a.d][u];Z coeff=Z(1)<<(a.d-e[i]);for(int b=0;b<=mons[a.d][u][j];b++){int v=ix[a.d][e[0]][e[1]];if(!right){for(int h=0;h<S(a.d);h++)r.at(v,h)+=coeff*a.at(u,h);}else{for(int h=0;h<S(a.d);h++)r.at(h,v)+=coeff*a.at(h,u);}if(b<mons[a.d][u][j]){coeff*=e[j];Z divisor=2*(b+1);assert(coeff%divisor==0);coeff/=divisor;}e[i]++;e[j]--;}}return r;}
struct Vtx{array<array<long double,3>,3>l,r;};
int nodes=0,leaves=0,maxdepth=0;long long limitnodes=20000;int limitdepth=30;
bool verify(P p,Vtx v,int depth){nodes++;maxdepth=max(maxdepth,depth);if(nodes%10==0)cerr<<"progress "<<nodes<<" depth="<<depth<<" leaves="<<leaves<<"\n";bool good=true;size_t neg=0;for(auto&c:p.a)if(c<0){good=false;neg++;}if(good){leaves++;return true;}if(depth>=limitdepth||nodes>=limitnodes){cerr<<"UNRESOLVED depth="<<depth<<" negative_coefficients="<<neg<<"\n";return false;}
  long double best=-1;bool side=false;int i=0,j=1;
  for(int s=0;s<2;s++){auto&w=s?v.r:v.l;for(int a=0;a<3;a++)for(int b=a+1;b<3;b++){long double dist=0;for(int c=0;c<3;c++)dist+=(w[a][c]-w[b][c])*(w[a][c]-w[b][c]);if(dist>best){best=dist;side=s;i=a;j=b;}}}
  Vtx vl=v,vr=v;auto&wl=side?vl.r:vl.l;auto&wr=side?vr.r:vr.l;
  for(int c=0;c<3;c++){long double mid=((side?v.r:v.l)[i][c]+(side?v.r:v.l)[j][c])/2;wl[i][c]=mid;wr[j][c]=mid;}
  P left=bisect(p,side,i,j);reduce(left);if(!verify(std::move(left),vl,depth+1))return false;
  P right=bisect(p,side,j,i);reduce(right);return verify(std::move(right),vr,depth+1);
}
int main(int argc,char**argv){init();int mode=argc>1?stoi(argv[1]):0;if(argc>2)limitdepth=stoi(argv[2]);if(argc>3)limitnodes=stoll(argv[3]);auto start=chrono::steady_clock::now();P p=make_target(mode);size_t neg=0,zeros=0;for(auto&c:p.a){neg+=c<0;zeros+=c==0;}cerr<<"initial negative="<<neg<<" zero="<<zeros<<"\n";Vtx v{{{{1,0,0},{0.5L,0.5L,0},{1.0L/3,1.0L/3,1.0L/3}}},{{{1,0,0},{0.5L,0.5L,0},{1.0L/3,1.0L/3,1.0L/3}}}};bool ok=verify(std::move(p),v,0);double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();cout<<(ok?"PASS":"UNRESOLVED")<<" mode="<<mode<<" nodes="<<nodes<<" leaves="<<leaves<<" depth="<<maxdepth<<" seconds="<<seconds<<"\n";return ok?0:2;}
