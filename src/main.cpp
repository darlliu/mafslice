#include "indexer.hpp"

int main (int argc, char** argv){
    seqdb s;
    s.import("./test/fasta/");
    std::cout << "Test Get: " << s.get(4e6+9900,4e6+10400) << std::endl;
    std::cout << "Test Get: " << s.get(4e6+900,4e6+1400) << std::endl;
    return 0 ;
}
