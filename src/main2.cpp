#include "mafindexer.hpp"
#include <boost/program_options.hpp>
int main (int ac, char** av){
    namespace po = boost::program_options;
    po::options_description desc("MAF Slicer ver 0.01. Give no argument to run a test!");
    std::string dbname, dbpath, gph, cfg;
    unsigned long cks;
    desc.add_options()
        ("help","list the arguments")
        ("test", "run the default test.")
        ("load", "Load a config file instead of creating a DB.")
        ("config", po::value<std::string>(&cfg) ->default_value("./db.cfg"), "Path to the config file. \n \
         Note: Config file stores relative paths. Make sure to open them in the correct folder!")

        ("create", "Create a new DB")
        ("dbname", po::value<std::string>(&dbname), "Name of the DB for the genomes")
        ("dbpath", po::value<std::string>(&dbpath), "Path to the FOLDER for which the DB is hosted")
        ("gph", po::value<std::string>(&gph), "Path to the FOLDER for which the genomes are stored as fasta files.\
         Note: Must use standard naming notations for the genome that correspond to chromosome names and end in fa or fasta postfix.")
        ("chunk", po::value<unsigned long>(&cks)->default_value(1e5), "Length of a chunk of sequence, default: 10,000")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << desc << std::endl;
        return 0;
    }

    if (vm.count("test")) {
        std::cerr << "running the test"<<std::endl;
        mafdb s;
        if (!s.load_db("./test/maf/db.cfg")){
            std::cerr << "Test load db failed! Importing instead!\n";
            s.import("./test/maf/fasta");
            s.init_tree();
        } else {
            std::cerr << "Test load db suceess!\n";
        }
        if (s.export_db("./test/maf/db.cfg"))
            std::cerr << "Test export Sucess! " <<std::endl;
        else
            std::cerr << "Test export Fail!! " <<std::endl;
    } else {
        std::cerr << "running the program"<<std::endl;
        std::cerr << "run it later!"<<std::endl;
        return 1;
        //if (vm.count("load")) {
            //seqdb s;
            //std::cerr << "Loading a config file" << cfg <<"\n";
            //if (!s.load_db(cfg)){
                //std::cerr << "Error loading the config file!"<<std::endl;
                //return 1;
            //} else {
                //std::cerr << "Test Get: " << s.get(4e6+9900,4e6+10400) << std::endl;
                //std::cerr << "Test Get: " << s.get(4e6+900,4e6+1400) << std::endl;
                //return 0 ;
            //}
        //} else if (vm.count("create")) {
            //seqdb s(dbname, cks, dbpath);
            //std::cerr << "Creating a db!"<<std::endl;
            //s.import(gph);
            //std::cerr << "Test Get: " << s.get(4e6+9900,4e6+10400) << std::endl;
            //std::cerr << "Test Get: " << s.get(4e6+900,4e6+1400) << std::endl;
            //if (s.export_db("./db.cfg"))
                //std::cerr << "Test export Success! " <<std::endl;
            //else
                //std::cerr << "Test export Fail!! " <<std::endl;
            //return 0;
        //} else {
            //std::cerr << "Terminating"<<std::endl;
            //return 1;
        //}
    }
    return 0 ;
}
