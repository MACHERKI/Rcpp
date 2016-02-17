#include <Rcpp.h>
# include<deque>
# include<algorithm>
using namespace Rcpp;


/*
 * This file contains some algorithms of sorting with verification.
 * All these algorithms are not performed as well as sort inner algorithm 
 * used in R except of sort_stl function.
 *          source:http//:www.wikipedia.org
 *                 MACHERKI M.E 
 *                  12/02/2016 
 */


// [[Rcpp::export]]
NumericVector sort_insert(NumericVector A) {
  int i,j;
  double x;
  for(i=1; i<A.size();i++){
    x=A[i];j=i-1;
    while((j>=0) & (A[j]>x)){ //without using swapping
      A[j+1]=A[j];
      j=j-1;
    }
  A[j+1]=x;
  
  }
  return A;
}
 
 // [[Rcpp::export]]
 NumericVector sort_bubble(NumericVector A) { 
   int n=A.size();
   int swap=1;
   double tmp;
   for(;;){
     if(swap==0)break;
     swap=0;
     for(int i=1;i<n;i++){
       if(A[i-1]>A[i]){
         tmp=A[i];A[i]=A[i-1];A[i-1]=tmp;
         swap=1;
       }
     }
   }
   return A;}

// [[Rcpp::export]]
NumericVector sort_selection(NumericVector A) { 
  int i,j;
  double tmp;
  for(j=0;j<A.size();j++){
    int imin=j;
    for(i=j+1;i<A.size();i++){
      if(A[i]<A[imin])imin=i;
    }
  if(imin!=j){
   tmp=A[j];A[j]= A[imin];A[imin]=tmp;
  }
  
  }
  return A;}
  
 

 // [[Rcpp::export]]
 NumericVector takudagaps( int sz ) { 
   NumericVector result;
   double p=1;
   int i=0;
   result.push_back(p);
   while(p<sz){
     p=ceil(p*2.25+1);
     result.push_back(p);
     i++;
   }
   
   return result;}
   // [[Rcpp::export]]
 NumericVector gaps( int sz ) { 
   NumericVector result;
   result.push_back(701);
   result.push_back(301);
   result.push_back(132);
   result.push_back(57);
   result.push_back(23);
   result.push_back(10);
   result.push_back(4);
   result.push_back(1);
   return result;}

// [[Rcpp::export]]
NumericVector sort_shell(NumericVector A) { 
  double temp=0;
  int N =A.size();
  int i,j,gap;
  NumericVector tacuda=gaps(0);
  int n=tacuda.size();
  for(int k=n-1;k>0;k--){
    gap=tacuda[k];
  for( i=gap;i<N;i++){
    temp=A[i];
    for( j=i;(j>=gap)&(A[j-gap]>temp);j-=gap){
      A[j]=A[j-gap];
    }
    A[j]=temp;
    
    }
    
  }
  return A;}

// [[Rcpp::export]]
  void sort_shell2(NumericVector A,int max) {//not efficace
  int stop,swap,limit;double temp;
  int x=(int)(max/2);
  while(x>0)
  {
  stop=0;
  limit=max-x;
  while(stop==0)
  {
  swap=0;
  for(int k=0;k<limit;k++)
  {
  if(A[k]>A[k+x])
  {
  temp=A[k];
  A[k]=A[k+x];
  A[k+x]=temp;
  swap=k;
  }
  }
  limit=swap-x;
  if(swap==0)stop=1;
  }
  x=(int)(x/2);
  }
  }
  
  
  

// [[Rcpp::export]]
NumericVector sort_stl(NumericVector A) { 
  
  NumericVector tmp=clone(A);
  std::sort(tmp.begin(),tmp.end());
  return tmp;
}
  
  
  // [[Rcpp::export]]
 void Quick(NumericVector arr,int left,int right) { 
  int i=left;int j=right;
  double tmp;
  double pivot=arr[(left+right)/2];
  /*partition*/
  while(i<=j){
    while(arr[i]<pivot)
      i++;
    while(arr[j]>pivot)
      j--;
      if(i<=j){
        tmp=arr[i];
        arr[i]=arr[j];
        arr[j]=tmp;
        i++;
        j--;
        
      }
    }
        /*recursion*/
  if(left<j)
    Quick(arr,left,j);
  if(i<right)
    Quick(arr,i,right);
}
    
    
    
    /*
     * This function is really comparable with the inner sort used in R
     * even the code is sample , the algorithm is good.
     */
    
    
    // [[Rcpp::export]]
    NumericVector sort_quick(NumericVector A) {
      NumericVector tmp=clone(A);
	   Quick(tmp,0,tmp.size()-1); //just starting from the middle of the vector
      return tmp;
    }
      
  
  /*
   * Armadillo sorting is availabel in arma.cpp file
   */
  
  
  
  
  




/*** R

###            not rapid methods, adapted to small sample size
x<-rnorm(100)
sort_insert(x)
sort (x)
sort_selection(x)
sort_bubble(x)
identical(sort_insert(x),sort (x))
identical(sort_bubble(x),sort (x))
identical(sort_selection(x),sort (x))
identical(sort_shell(x),sort (x))

  */


/*
 * according to:
 *  Tzc<-function(x,type="z"){
 if(type=="z")return ((x-mean(x,na.rm =T))/sd(x,na.rm =T)) 
if(type=="c"){
T0=(max(x,na.rm =T)+min(x,na.rm =T))/2
dt=(max(x,na.rm =T)-min(x,na.rm =T))/2
return((x-T0)/dt)
}
}
 
 even in the case of T0=dt=cte (min=0) we can say centred x in y as y=(x/(max*0.5))*0.5
 x/(max/2)
 to just get the interval [0,1]
 */
 
 
 
 
 //ceiling function exist  and can be manipulated like R. 


// [[Rcpp::export]]
NumericVector centrage(NumericVector A) {
  const double K=pow(5,20);
  return A/K;}
 

/* To get a good and rapid matching method you have to get the best sorting algorithm.
The binary matching is with no challenge the best way to matching.
Also unique duplicated intersection are all simply overcame with same concept.  
 */
 // [[Rcpp::export]]
 
LogicalVector BINin(NumericVector A,NumericVector B) {//previously sorted same as BASIX method 
const int a=A.size();const int b=B.size();
int j=0;int i=0;                        //global variables 
LogicalVector res(a,0);
for(i=i;i<a;i++){
for(j=j;(j<b)&(B[j]<=A[i]);j++){
if(B[j]==A[i]){res[i]=1;break;}         // 
}
}
return res;}

// we use just the same code to find the index(starting from 0) reported by j
// [[Rcpp::export]]
NumericVector BINmatch(NumericVector A,NumericVector B) {//previously sorted same as BASIX method
const int a=A.size();const int b=B.size();
int j=0;int i=0;
NumericVector res(a,-1.0);
for(i=i;i<a;i++){
for(j=j;(j<b)&(B[j]<=A[i]);j++){
if(B[j]==A[i]){res[i]=j;
  break;}
}
}
return res;}
// [[Rcpp::export]]
NumericVector sort_deque (std::vector<double> A) {  // dividing A into to set by sign 
std::deque<double> x;
  double tmp=0.0;
  int i,j,POS,NEG,Size;
   Size=A.size();
  POS=NEG=0;
  for( i=0;i<Size;i++){                      // also i is the last index of x
    if(A[i]>0.0){
       x.push_back(A[i]);
	   POS++;                                  //number of positive value
       for( j=i;(j<NEG+1) & (x[j]>x[j-1]);j--){
            tmp=x[j];x[j-1]=x[j];x[j]=tmp;     //swapping 
                                         }
				}
    else{
	    NEG++;
        x.push_front(A[i]);
		for( j=0;(j>POS+1) & (x[j+1]>x[j]);j++){
            tmp=x[j];x[j+1]=x[j];x[j]=tmp;     //swapping
                                         }
	    }
		
                           }
  return wrap(x) ;
}
// need to be performed: we have to the place of first positive value

// [[Rcpp::export]]
void sort_bucket (NumericVector A,const int n) {  // dividing A into to set by sign 
   
   //create an empty buckets
   std::vector<double>b[n]; int bi;
   //put A elements in different buckets
   for(int i=0;i<n;i++){
       bi=n*A[i];//index in bucket
	   b[bi].push_back(A[i]);
	                   }
    // sort individual buckets
	for(int i=0;i<n;i++) std::sort(b[i].begin(),b[i].end());
	 //concatenation 
    int index=0;
	for(int i=0;i<n;i++){
	   for(int j=0;j<b[i].size();j++){
	      A[index++]=b[i][j];}
	}
		                                   }
	
	
 	  //this function do not run 
	   
   
   // [[Rcpp::export]]
   inline const int to_int  (const double x){// this function chow that saving double as int is equivalent to trunc function
         return x; 
   }
    
    



/***R
    
# CPP had the option digit=15
    to_int(5.999)
    to_int(5.9999999999999999999999)

    */  
   
