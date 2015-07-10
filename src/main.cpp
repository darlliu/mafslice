#include "indexer.hpp"
#include <boost/program_options.hpp>
int main (int ac, char** av){
    namespace po = boost::program_options;
    po::options_description desc("MAF Slicer ver 0.01. Give no argument to run a test!");
    desc.add_options()
        ("help","list the arguments")
        ("run", "run the program instead of testing.")
        ("dbname", po::value<std::string>(), "Name of the db for the genomes")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    if (vm.count("run")) {
        std::cout << "running the program"<<std::endl;
    } else {
        std::cout << "running the test"<<std::endl;
        seqdb s;
        if (!s.load_db("./test/db.cfg")){
            std::cout << "Test load db failed! Importing instead!\n";
            s.import("./test/fasta/");
        } else {
            std::cout << "Test load db suceess!\n";
        }
        std::cout << "Test Get: " << s.get(4e6+9900,4e6+10400) << std::endl;
        std::cout << "Test Get: " << s.get(4e6+900,4e6+1400) << std::endl;
        if (s.export_db("./test/db.cfg"))
            std::cout << "Test export sucess! " <<std::endl;
        else
            std::cout << "Test export Fail!! " <<std::endl;
    }
    return 0 ;
}
