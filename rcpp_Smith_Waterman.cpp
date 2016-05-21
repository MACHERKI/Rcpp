//---------------------------Smith & Waterman Algorithm:Local Alignment-------------------------------
//-----------------------           MACHERKI M.E.       ------------------------------- 
/*
 * 
 */



#include <string>
#include <Rcpp.h>
#include <stack>
using namespace Rcpp;
double similarity_score(char a,char b,double match_score ,double mismatche_score);
  
NumericVector Max_index(NumericMatrix H,int N_a,int N_b);

// [[Rcpp::export]]
NumericMatrix SWM(std::string seq_a,std::string seq_b,int N_a ,int N_b,double match_score,double mismatche_score,double gap_score  ){
  
  
  // initialize H
  
  // this is overcoming becuase the initial set of a matrix are 0;  
  NumericMatrix H(N_a+1,N_b+1);     
  /*
  for(int i=0;i<=N_a;i++){
  for(int j=0;j<=N_b;j++){
  H(i,j)=0.;
  }
  } 
  */
  NumericVector temp(4);
  
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      temp[0] = H(i-1,j-1)+similarity_score(seq_a[i-1],seq_b[j-1],match_score,mismatche_score); 
      temp[1] = H(i-1,j)-gap_score; // gap_score called also delta                   
      temp[2] = H(i,j-1)-gap_score;                 
      temp[3] = 0.;
      H(i,j) = max(temp);
      
    }
  }
  
  
return H;
  }
  ////////////////////////////////////////////////////////////////////////

  // [[Rcpp::export]]
  inline CharacterVector get_optimal_alignment(std::string A,std::string B,double match_score,double mismatche_score,double gap_score)
  {
    
    int N_a = A.length();                     // get the actual lengths of the sequences
    int N_b = B.length();
  
  NumericMatrix dp=SWM(A,B,N_a,N_b,match_score,mismatche_score,gap_score);
  NumericVector indexes=Max_index(dp,N_a,N_b);
  std::string retA, retB;
  std::stack<char> SA, SB;
  int ii = indexes[0]; int jj = indexes[1];
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
      int S = (A[ii-1] == B[jj-1]) ? match_score : -mismatche_score;  // mismatche_score called also mu
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

  
  
  
  
  
  
  
  
  
  
  
  
    //////////////////////////////////////////////////////
  double similarity_score(char a,char b,double match_score ,double mismatche_score){
    
    double result;
    if(a==b){
      result=match_score;
    }
    else{
      result=-mismatche_score;
    }
    return result;
  }
//////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector Max_index(NumericMatrix H,int N_a,int N_b){
  NumericVector res(2);
  double H_max = 0.0;
  int i_max=0,j_max=0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      if(H(i,j)>H_max){
        H_max = H(i,j);
        i_max = i;
        j_max = j;
      }
    }
  }
res[0]=i_max;
res[1]=j_max;
return res;
  }
/*** R
seq1<-"ATTAGGGACCCTAGGATCATTAGCCTGGAAATTACCCCCAAATTAAATTTAAATAAAAAAAAAAAAAAAAA"
seq2<-"CAAATTTACCCAGCGGGTACCCGATTGACCTGACTGGGGATGATTTACCCCAAGGGACCCATAAAAATTTTTTTTTTTTAAAAAAA"
get_optimal_alignment(seq1,seq2,1,0,5)
  */