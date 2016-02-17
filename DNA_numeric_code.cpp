
/*
In general, DNA sequences are  stored as character data into different formats
of files like FASTA,FNA... ALso numeric coding of DNA in base 4 exists  such 
'A'=0,'G'=1,'C'=2 and 'T'=3. But why I use this kind of coding? 
The first step of sequence alignment algorithms such BLAST or FASTA consists on 
stepwise splitting of the DNA sequence into oligos of 11 bases. This operation 
influences  the execution  running time of the program.
This is an example of codes  ;the first use string manipulation and the second 
use numeric coding in base 5(to overcame  N in the sequence).



                     MACHERKI M.E
                     14/02/2016
 */
 
 /*
 This source file  uses C++ code.You have to install Rcpp package to run it.
 A FASTA file  path is required to test this code in R.Please check the bottom 
 of the page.
 */ 
 
#include <fstream>
#include<Rcpp.h>
#include<string.h>
#include <sstream>
#include<cctype>
#include <algorithm>
using namespace Rcpp;

/*
First of all, this is a function to read FASTA file into R
-----!!!!!!!!!!!!! one sequence only!!!!!!!!!!!!!!-----------
*/
//[[Rcpp::export]]
inline CharacterVector read_fasta(std::string file ){
CharacterVector records;
std::ifstream in(file.c_str());
in.get();// remove first >
std::string rec ;
while(getline(in,rec,'>')){
int newLineLoc=rec.find('\n');
std::string header =rec.substr(0,newLineLoc);
std::string sequence= rec.substr(newLineLoc+1,rec.length()-newLineLoc-2);
sequence.erase(std::remove(sequence.begin(),sequence.end(),'\n'),sequence.end());
std::transform(sequence.begin(),sequence.end(),sequence.begin(),tolower);
records[header]=sequence;
}
return records;
}
/*
String splitter function into N size oligos (Only for lower case string).

*/
// [[Rcpp::export]]
 inline CharacterVector SPLIT(const std::string  strings,const int n){
    int lim=strings.length()-n+1 ;//size of sub string vector
    CharacterVector resultat(lim) ;
          for ( int j=0; j<lim ;j++) resultat[j] = strings.substr(j,n) ;
return resultat ;
} 
/*
String splitter function into N size oligos using numeric coding.
I used a simple hashing function which produce an unique numeric 
identifier for each oligo.It is possible to use integer but I use
double to accept a large N.
(See R help:strtoi or just type ?strtoi in console)

*/

inline double  natoi( const std::string A,const int n ){
//Used as input of numeric_split(see next function)
  int np=n-1;
  IntegerVector temp(n,4);
  char toc;
  for( int i=0;i<n;i++){
    toc=A[i];
if(toc=='a') temp[i]= 0;
if(toc=='t') temp[i]= 3;
if(toc=='c') temp[i]= 1;
if(toc=='g') temp[i]= 2;
}
double  result=0;
 for( int i=0;i<n;i++) result+= pow(5,(np-i))*temp[i];
 return result;
}

// [[Rcpp::export]]
inline NumericVector  numeric_split(const std::string A,const int n)
    {
  const int np=n-1;
  const int Lim=A.length()-np;//-- number of subsequence 
  NumericVector res(Lim);     //-- result vector
  double p=0;                 //-- temp result accumulator
  res[0]=p=natoi(A,n);        //-- first result 
  const double  N=pow(5,np);  //-- cont used in recursion formula
  unsigned short int x,y;     //-- temp to allocate forward and backward values
  char fo,ba;                 //--forward and backward chars
  for( int i=1;i<Lim;i++)
     {
       fo=A[i-1];ba=A[i+np];
       x=y=4;                 //-- empty case:n,p,m... 
       if(fo=='a')x=0;        //-- char detector 
       if(fo=='t')x=3;        // 
       if(fo=='g')x=2;        // 
       if(fo=='c')x=1;
       if(ba=='a')y=0;
       if(ba=='t')y=3;
       if(ba=='g')y=2;
       if(ba=='c')y=1;
     res[i]=p= 5*p-5*N*x+y;   //-- scalar formula 
     }
  return res;
  }
/*
R verification :
I use a file named k12.fna in my computer (full E.coli k12 genome ). 
You have just to change k12.fna into your file name.
*/


/***R
x<-read_fasta("k12.fna")  ##change me to your file name 
system.time(charsplit<-SPLIT(x,11))
#user  system elapsed 
# 5.38    0.02    5.39 
system.time(numsplit<-numeric_split(x,11))
#user  system elapsed 
#0.08    0.01    0.09    ## more then 50 time faster
  
  */



