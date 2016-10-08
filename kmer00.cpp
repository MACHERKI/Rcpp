
//---------------------------kmer counting using perfect hash table-------------------------------
//-----------------------           MACHERKI M.E.       ------------------------------- 


/*
 * The same code can be using in C++:
 * IntegerVector is replaced by std::vector<int> or simply int[]
 */


#include <Rcpp.h>
using namespace Rcpp;


inline const short int V (const char x){
  switch(x){
  case 'A':case 'a':
    return 0;
    break;
  case 'C':case 'c':
    return 1;
    break;   
  case 'G':case 'g':
    return 2;
    break;
  default:
    return 3;
  break;  
  }
  
}

inline unsigned int  X0( const std::string A,const int k ,const int n){
  unsigned int  result=0;
  int j=k;
  for( int i=n-1;i>n-k-1;i--) {
    result+= pow(4,k-j)*V(A[i]);
  j--;
  }
  return result;
}

// [[Rcpp::export]]
inline IntegerVector kmer4(const std::string A,const int n,const int k)
{
  
  IntegerVector P(pow(4,k));                  
  int x=X0(A,k,n);                              
  P[x]++;                   
  const int N=pow(4,k-1);               
  for( int i=n-k-1;i>-1;i--){
    x=N*V(A[i])+x/4-x%4/4;
    P[x]++;
  }
  return P;
}
