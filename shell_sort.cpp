
/*
 * The method starts by sorting pairs of elements far apart from each 
 * other, then progressively reducing the gap between elements to be compared.
 * The running time of Shellsort is heavily dependent on the gap sequence it uses.
 * Many gap function are evaluable see Link.
 * source file:link
 * 
 * The best profile was effectuated when we use sort_shell_sed function
 * 
 *                                 MACHERKI M E.
 *                                 17/02/2016
 * 
 * 
 */

#include <Rcpp.h>
using namespace Rcpp;



//The following code is the typical body of the shell sort algorithm:
inline void shell ( NumericVector unsorted, int size, NumericVector Gaps)
{
  int gn=Gaps.size();
  int j,gap ;
  
  for (int t = gn-1; t > -1; t--)
  {
    gap=Gaps[t];
    for (int i = gap; i < size; ++i)
    {
      double temp = unsorted[i];
      for (j = i; j >= gap && temp < unsorted[j - gap]; j -= gap)
      {
        unsorted[j] = unsorted[j - gap];
      }
      unsorted[j] = temp;
    } 
  }
}
//******************gaps functions*************************

//Many gaps function can be used 
/*Shell sort function using slicing by two gaps (Shell, 1959)
 Gaps= N/2^k={N/2,N/4,N/8...1}
*/
inline NumericVector gapshell( int size ) {
  NumericVector result;
  for(int i=size/2;i>0;i/=2)
    result.push_back(i);
  return rev(result);}
/*Method of Frank & Lazarus, 1960
 Gaps=2*(N/2^k+1)+1
*/
inline NumericVector gapFrank( int N ) {
  int tmp=10;
  NumericVector result;
  for(int k=2;tmp>1;k++){
    tmp=2*(N/pow(2,k))+1;
    result.push_back(tmp);}
  return rev(result);}

/*Method of Hibbard, 1963
 Gaps=2^k-1
 */
inline NumericVector gapHibbard( int N ) {
  NumericVector result;
  int tmp=0;
  for(int k=1;tmp<N;k++){
    tmp=pow(2,k)-1;
    result.push_back(tmp);}
  return result;}

/*Method of Pratt, 1971
 Gaps=(3^k-1)/2
 */
inline NumericVector gappratt( int N ) {
  int n=N/3;
  NumericVector result;
  int tmp=0;
  for(int k=1;tmp<n;k++){
    tmp=(pow(3,k)-1)/2;
    result.push_back(tmp);}
  return result;}
/*Method of Sedgewick, 1986
 Gaps=4^k+3*2^(k-1)+1
 */
inline NumericVector gapsed( int N ) {
  NumericVector result;
  result.push_back(1);
    int tmp=0;
  for(int k=1;;k++){
    tmp=pow(4,k)+(3*pow(2,k-1))+1;
    if(tmp>N)break;
    result.push_back(tmp);}
  return result;}
/*Method of Tokuda, 1992
 Gaps=(9^k-4^k)/(5*4^(k-1))
*/

inline NumericVector gaptocuda( int N ) {
  NumericVector result;
  int tmp=0;
  for(int k=1;;k++){
    tmp=(pow(9,k)-pow(4,k))/(5*pow(4,k-1));
    if(tmp>N)break;
    result.push_back(tmp);}//we can use tmp+1
  result[0]=1;
  return result;}
/*Method of Tokuda, 1992:generic
 Gaps=h[i]=(h[i-1])*2.25+1; h[1]=1,1 base indexing
*/

inline NumericVector gaptocuda_emp( int N ) {
  NumericVector result;
  result.push_back(1);
  int tmp=1;
  for(int k=1;;k++){
    tmp=2.25*tmp+1;
    if(tmp>N)break;
    result.push_back(tmp);}//we can use tmp+1
  return result;}





//using empirically derived gaps (Ciura, 2001) N<2000

NumericVector gapCiura( int sz ) { 
  NumericVector result;
  result.push_back(1);
  result.push_back(4);
  result.push_back(10);
  result.push_back(23);
  result.push_back(57);
  result.push_back(132);
  result.push_back(301);
  result.push_back(701);
  result.push_back(1705);
  return result;}

/******************sorting functions*************************
************************************************************/







// [[Rcpp::export]]
NumericVector sort_shell_gap(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gapshell(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_frunk(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gapFrank(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_hibbard(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gapHibbard(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_patt(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gappratt(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_sed(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gapsed(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_tacuda(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gaptocuda(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_tacudaemp(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector gaps=gaptocuda_emp(N);
  NumericVector result=clone(unsorted);
  shell(result,N,gaps);
  return result;}
// [[Rcpp::export]]
NumericVector sort_shell_Ciura(NumericVector unsorted) { //run
  int N =unsorted.size();
  NumericVector Ciura=gapCiura(0);
  NumericVector result=clone(unsorted);
  shell(result,N,Ciura);
  return result;}
/***R
##trial mode
x<-runif(1000)
base<-sort(x)
a<-sort_shell_gap(x)
b<-sort_shell_frunk(x)
c<-sort_shell_hibbard(x)
d<-sort_shell_patt(x)
e<-sort_shell_sed(x)
f<-sort_shell_tacuda(x)
g<-sort_shell_tacudaemp(x)
h<-sort_shell_Ciura(x)
all.equal.list(base,a,b,c,e,f,g,h)
## timing :
x<-runif(1000000)
system.time(base<-sort(x))
system.time(a<-sort_shell_gap(x))
system.time(b<-sort_shell_frunk(x))
system.time(c<-sort_shell_hibbard(x))
system.time(d<-sort_shell_patt(x))
system.time(e<-sort_shell_sed(x))
system.time(f<-sort_shell_tacuda(x))
system.time(g<-sort_shell_tacudaemp(x))
system.time(h<-sort_shell_Ciura(x))

*/
  
