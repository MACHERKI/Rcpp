//---------------------------Algorithm:Levenshtein distance between two string-------------------------------
//-----------------------           MACHERKI M.E.      Rcpp codes  ------------------------------- 



#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]

int DamerauLevenshteinDistance(std::string str1, std::string str2){	
	int lenStr1=str1.size();
	int lenStr2=str2.size();
	NumericMatrix d(lenStr1+1,lenStr2+1);
	for(int i=0;i<lenStr1+1;i++)d(i,0)=i;
	for(int j=0;j<lenStr2+1;j++)d(0,j)=j;
	int  cost=0;
	NumericVector tmp(3);
	NumericVector tmp1(2);
	for(int i=1;i<lenStr1+1;i++){
		for(int j=1;j<lenStr2+1;j++){
	        if(str1[i-1]==str2[j-1]) cost=0;
			else cost=1;
			tmp[0]=d(i-1, j  ) + 1;
			tmp[1]=d(i  , j-1) + 1;
			tmp[2]=d(i-1  , j-1) + cost;
			d(i,j)=min(tmp);
	
            if((i>2)&(j>2) & (str1[i-1]==str2[j-2])&(str1[i-2]==str2[j-1])){
				tmp1[0]=d(i,j);
				tmp1[1]=d(i-2,j-2)+1;
				d(i,j)=min(tmp1);
			}
	
	}
}
return d(lenStr1,lenStr2);
}

/*
 * Iterative with full matrix
  Wagner-Fischer algorithm
 */


// [[Rcpp::export]]

int LevenshteinDistance(std::string str1, std::string str2){	
  int lenStr1=str1.size();
  int lenStr2=str2.size();
  NumericMatrix d(lenStr1+1,lenStr2+1);
  for(int i=0;i<lenStr1+1;i++)d(i,0)=i;
  for(int j=0;j<lenStr2+1;j++)d(0,j)=j;
  int  cost=0;
  NumericVector tmp(3);
  NumericVector tmp1(2);
  for(int i=1;i<lenStr1+1;i++){
    for(int j=1;j<lenStr2+1;j++){
      if(str1[i-1]==str2[j-1]) cost=0;
      else cost=1;
      tmp[0]=d(i-1, j  ) + 1;
      tmp[1]=d(i  , j-1) + 1;
      tmp[2]=d(i-1  , j-1) + cost;
      d(i,j)=min(tmp);
      
      
      
    }
  }
  
  int res=d(lenStr1,lenStr2);
  return res;
}

/*
 * It turns out that only two rows of the table are needed for the construction if one does not want 
 to reconstruct the edited input strings (the previous row and the current row being calculated).
 The Levenshtein distance may be calculated iteratively using the following algorithm:
 */ 
// [[Rcpp::export]]
int LevenshteinDistance1(std::string s, std::string t)
{
    // degenerate cases
    if (s == t) return 0;
    if (s.length() == 0) return t.length();
    if (t.length() == 0) return s.length();

    // create two work vectors of integer distances
    IntegerVector v0(t.length() + 1);
    IntegerVector v1(t.length() + 1);
    IntegerVector tmp(3);
    // initialize v0 (the previous row of distances)
    // this row is A[0][i]: edit distance for an empty s
    // the distance is just the number of characters to delete from t
    for (int i = 0; i < v0.length(); i++)
        v0[i] = i;

    for (int i = 0; i < s.length(); i++)
    {
        // calculate v1 (current row distances) from the previous row v0

        // first element of v1 is A[i+1][0]
        //   edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1;

        // use formula to fill in the rest of the row
        for (int j = 0; j < t.length(); j++)
        {
            int cost = (s[i] == t[j]) ? 0 : 1;
            tmp[0]=v1[j] + 1;
            tmp[1]=v0[j + 1] + 1;
            tmp[2]=v0[j] + cost;
            v1[j + 1] = min(tmp);
        }

        // copy v1 (current row) to v0 (previous row) for next iteration
        for (int j = 0; j < v0.size(); j++)
            v0[j] = v1[j];
    }

    return v1[t.length()];
}


                        