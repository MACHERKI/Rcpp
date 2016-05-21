
//---------------------------Algorithm: Needleman & Wunsch:Global Alignment-------------------------------
//-----------------------           MACHERKI M.E.      Rcpp codes  ------------------------------- 


#include <string>
#include <stack>
#include <Rcpp.h>
using namespace Rcpp;


inline int max(int a,int b){
  if(a>=b) {return a;}
  else{return b;}
}
// [[Rcpp::export]]
inline NumericMatrix needleman_wunsch(std::string A,std::string B,int n,int m,int match_score, int mismatch_score,int gap_score)
{
  
  NumericMatrix dp(n+1,m+1);
  for (int i=0;i<=n;i++) dp(i,0) = dp(0,i) = -i * gap_score;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=m;j++)
    {
      int S = (A[i-1] == B[j-1]) ? match_score : -mismatch_score;
      dp(i,j) = max(dp(i-1,j-1) + S, max(dp(i-1,j) - gap_score, dp(i,j-1) - gap_score));
    }
  }
  return dp;
}

// [[Rcpp::export]]
inline CharacterVector get_optimal_alignment(std::string A,std::string B,int match_score, int mismatch_score,int gap_score)
{
  int n=A.size();int m=B.size();
  NumericMatrix dp=needleman_wunsch(A,B,n,m, match_score,  mismatch_score, gap_score);
  
  std::string retA, retB;
  std::stack<char> SA, SB;
  int ii = n; int jj = m;
  while (ii != 0 || jj != 0)
  {
    if (ii == 0)
    {
      SA.push('-');
      SB.push(B[jj-1]);
      jj--;
    }
    else if (jj == 0)
    {
      SA.push(A[ii-1]);
      SB.push('-');
      ii--;
    }
    else
    {
      int S = (A[ii-1] == B[jj-1]) ? match_score : -mismatch_score;
      if (dp(ii,jj) == dp(ii-1,jj-1) + S)
      {
        SA.push(A[ii-1]);
        SB.push(B[jj-1]);
        ii--; jj--;
      }
      else if (dp(ii-1,jj) > dp(ii,jj-1))
      {
        SA.push(A[ii-1]);
        SB.push('-');
        ii--;
      }
      else
      {
        SA.push('-');
        SB.push(B[jj-1]);
        jj--;
      }
    }
  }
  while (!SA.empty())
  {
    retA += SA.top();
    retB += SB.top();
    SA.pop();
    SB.pop();
  }
  CharacterVector result(2);
  result[0]=retA;result[1]=retB;
  return result;
}

/*** R
seq1<-"ATTAGGGACCCTAGGATCATTAGCCTGGAAATTACCCCCAAATTAAATTTAAAT"
seq2<-"CAAATTTACCCAGCGGGTACCCGATTGACCTGACTGGGGATGATTTACCCCAAGGGACCCAT"
get_optimal_alignment(seq1,seq2,2,1,1)
n<-nchar(seq1)
m<-nchar(seq2)
needleman_wunsch(seq1,seq2,n,m,2,1,1)
*/
