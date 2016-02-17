/*
All these function are useful in Reading file.read_file_delime function
 read the file from start to end.'Tol' is an option to make string into 
 lower case.
                              MACHERKI M E.
                              16/02/2016
*/


#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include<cctype>
#include<Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
CharacterVector read_file(std::string path){
std::ifstream t(path.c_str());//connect with file 
std::stringstream ss;
ss<<t.rdbuf();                // scan file
return ss.str();
}
//[[Rcpp::export]]
CharacterVector read_file1(std::string path){
std::ifstream in(path.c_str());
std::string contents;
in.seekg(0,std::ios::end);
contents.resize(in.tellg());
in.seekg(0,std::ios::beg);
in.read(&contents[0],contents.size());
in.close();
return contents;
}
// function to read file from start to end pointer
//[[Rcpp::export]]
CharacterVector read_file_delime(std::string path,int start,int end,int TOL){
std::ifstream in(path.c_str());
std::string contents;
in.seekg(0,std::ios::end);
contents.resize(in.tellg());
in.seekg(start,std::ios::beg);
in.read(&contents[0],end);
in.close();
contents.erase(std::remove(contents.begin(),contents.end(),'\n'),contents.end());
//remove lines delimiter
if(TOL>0){//force to lower case transformation
std::transform(contents.begin(),contents.end(),contents.begin(),tolower);
}
return contents.substr(start-1,end-start);
}
/* read a file from start  to end*/
/* set the file path before using */

/*** R
###*******************choose file to read************************####
con<-utils::choose.files()
system.time(FR<-readLines(con))     #vector as the split(x,'\n') 
system.time(Fcpp1<-read_file(con))  #contents Line delimiter
system.time(Fcpp2<-read_file1(con)) #contents Line delimiter
#one string without \n and in lower case
system.time(Fcpp3<-read_file_delime(con,20,100,1)) 

*/ 
