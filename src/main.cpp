#include "mafindexer.hpp"
#include <boost/program_options.hpp>
int main (int ac, char** av){
    namespace po = boost::program_options;
    po::options_description desc("MAF Slicer ver 0.01. Use cases:\n\
            mafslice.exe --load --config CONFIG.CFG [--test]\n\
            mafslice.exe --create --config CONFIG.CFG --genome-folder FOLDER --dbpath FOLDER --dbname NAME [--msa --ref REF]|[--chunk SZ]\n\
            mafslice.exe --create --config CONFIG.CFG --genome-file FILE --dbpath FOLDER --dbname NAME [--msa --ref REF]|[--chunk SZ]\n\
            mafslice.exe --assemble --config CONFIG.CFG --genome-folder FOLDER --dbpath FOLDER --dbname NAME [--msa --ref REF]|[--chunk SZ]");
    std::string dbname, dbpath, gph, gpf,cfg,ref;
    unsigned long cks;
    desc.add_options()
        ("help,h","list the arguments")
        ("create", "Create a new DB")
        ("test", "Run a test at the end of the program")
        ("assemble", "Assemble pre-existing DB files into DB. Must also provide the genome files location.")
        ("load,L", "Load a config file instead of creating a DB.\nNote: If you choose this mode, you cannot create, assemble DB.")
        ("config,C", po::value<std::string>(&cfg) ->default_value("./db.cfg"),
         "Path to the config file.\nNote: Config file stores relative paths. Make sure to open them in the correct folder!")
        ("msa", "Create a MSA DB on maf files instead of a Seq DB on fasta files")
        ("dbname,N", po::value<std::string>(&dbname), "Name of the DB for the genomes")
        ("dbpath,P", po::value<std::string>(&dbpath), "Path to the FOLDER for which the DB is hosted")
        ("genome-folder,F", po::value<std::string>(&gph), "Path to the FOLDER for which the genomes are stored as fasta/maf files.")
        ("genome-file,f", po::value<std::string>(&gpf), "Path to the FILE for which the genome file is stored as fasta/maf.")
        ("ref,r", po::value<std::string>(&ref)->default_value("mm10"), "Reference genome for MSA, default: mm10")
        ("chunk,c", po::value<unsigned long>(&cks)->default_value(1e5), "Length of a chunk of sequence for SeqDB, default: 10,000")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << desc << std::endl;
        return 0;
    }

    if (vm.count("msa"))
    {
        mafdb s(dbname, dbpath, ref);
        std::cerr << "running the program"<<std::endl;
        if (vm.count("load")) {
            std::cerr << "Loading a config file" << cfg <<"\n";
            if (!s.load_db(cfg)){
                std::cerr << "Error loading the config file!"<<std::endl;
                return 1;
            }
        } else if (vm.count("create")) {
            std::cerr << "Creating a db!"<<std::endl;
            if (vm.count("genome-file"))
                s.import_chr(gpf);
            else
                s.import(gph);
        } else if (vm.count("assemble")) {
            std::cerr << "Assembling a db!"<<std::endl;
            s.toggle_assemble ();
            s.import(gph);
        } else {
            std::cerr << "No mode selected, terminating"<<std::endl;
            return 1;
        }
        if (s.export_db(cfg))
            std::cerr << "Saved config! " <<std::endl;
        else
            std::cerr << "Save config failed! " <<std::endl;
        if (vm.count("test"))
        {
            s.init_tree();
            //s.get("chr1",3218024,3218037);
            s.get("chr1",4417696,4417709);

        }
    }
    else
    {
        seqdb s(dbname, cks, dbpath);
        std::cerr << "running the program"<<std::endl;
        if (vm.count("load")) {
            std::cerr << "Loading a config file" << cfg <<"\n";
            if (!s.load_db(cfg)){
                std::cerr << "Error loading the config file!"<<std::endl;
                return 1;
            }
        } else if (vm.count("create")) {
            std::cerr << "Creating a db!"<<std::endl;
            if (vm.count("genome-file,f"))
                s.import_chr(gpf);
            else
                s.import(gph);
        } else if (vm.count("assemble")) {
            std::cerr << "Assembling a db!"<<std::endl;
            s.toggle_assemble ();
            s.import(gph);
        } else {
            std::cerr << "No mode selected, terminating"<<std::endl;
            return 1;
        }
        if (s.export_db(cfg))
            std::cerr << "Saved config! " <<std::endl;
        else
            std::cerr << "Save config failed! " <<std::endl;
        if (vm.count("test"))
        {
            s.get("chrY",1350000,1380000);
        }
    }
    return 0 ;
}
