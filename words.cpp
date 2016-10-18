#include <Rcpp.h>
using namespace Rcpp;


/*All combination of characters  of size k  belonging in the alphabet alpha  
 * 
 */

// k<15

// [[Rcpp::export]]
CharacterVector Words(const int k) {
 const  int N=pow(4,k);
  CharacterVector res(N);
  char s[14]={'A'};
  for(int i=0;i<k;i++)s[i]='A';
  res[0]=s;
   int z=k-1;
  for(int i=1;i<N;i++){
    for(int j=z;j>-1;j--){
      if(s[j]=='A'){
        s[j]='C';
        break;
      }
      if(s[j]=='C'){
        s[j]='G';
        break;
      }
      if(s[j]=='G'){
        s[j]='T';
        break;
      }
      if(s[j]=='T')s[j]='A';
        
    }
  res[i]=s;
  }
  return res;
  }
// [[Rcpp::export]]
CharacterVector alpha_Words(const int k,const std::string alpha) {
  const int n=alpha.size();
  const  int N=pow(n,k);
  CharacterVector res(N);
  char s[14]={alpha[0]};
  for(int i=0;i<k;i++)s[i]=alpha[0];
  res[0]=s;
  int z=k-1;
  for(int i=1;i<N;i++){
    for(int j=z;j>-1;j--){
      if(s[j]==alpha[n-1])s[j]=alpha[0]; 
      else{
      for(int t=0;t<n-1;t++){
        if(s[j]==alpha[t]){
          s[j]=alpha[t+1];
          j=-1;
      }
        
      }
      }
       
    }
    res[i]=s;
  }
  return res;
}


// [[Rcpp::export]]
CharacterVector Words1(const int k) {
  const int N=pow(4,k);  // size of the vector
  CharacterVector res(N);  // countenair for the result
  char s[14]={'A'}; //first value ,14 to support k<15
  for(int i=0;i<k;i++)s[i]='A';
  res[0]=s;
   int z=k-1;
  for(int i=1;i<N;i++){
    for(int j=z;j>-1;j--){
      switch(s[j]){
      case 'A':
    s[j]='C';
    j=-1;
        break;
      case 'C':
    s[j]='G';
        j=-1;
        break;
      case 'G':
    s[j]='T';
        j=-1;
        break;
      default:
        s[j]='A';
        break;
      }
    }
    res[i]=s;
  }
  return res;
}



/*** R

# library(rbenchmark)
# library(seqinr)
# library(Biostrings)
# evluation of All combination of  character in the DNA alphabet
# 
# benchmark(Words(8),Words1(8),alpha_Words(8,"ACGT"),mkAllStrings(c("A","C","G","T"), 8),seqinr::words(8))
#                                        test replications elapsed relative user.self
# 3                 alpha_Words(8, "ACGT")          100    0.96    1.000      0.95     0.00
# 4 mkAllStrings(c("A", "C", "G", "T"), 8)          100    3.98    4.146      3.97     0.01
# 5                       seqinr::words(8)          100    3.53    3.677      3.53     0.00
# 1                               Words(8)          100    0.97    1.010      0.97     0.00
# 2 






  */