#include <fstream>
#include <Rcpp.h>

using namespace Rcpp;
//it work i get fool in my old one
// [[Rcpp::export]]
CharacterVector read_fasta(std::string file) {
  CharacterVector records;
  std::ifstream in(file.c_str());
  in.get(); // remove first '>'
  std::string rec;
  while(getline(in,rec,'>')){
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0,newLineLoc);
    std::string sequence = rec.substr(newLineLoc+1, rec.length()-newLineLoc-2);
    sequence.erase(std::remove(sequence.begin(),sequence.end(),'\n'),sequence.end());
    records[header]=sequence;
  }
  return(records);
}