//
//#include "../../input.h"
//#include "../../namics.h"
//#include "../../solve_scf.h"
//#include "../../system.h"
//#include "../../output.h"
//
//#include "../cleng.h"

#include <fstream>

#include "checkpoint.h"
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

using namespace std;

const string Checkpoint::IN_CLASS_NAME = "checkpoint";
const string Checkpoint::DEFAULT_CHECKPOINT_PATH = "./checkpoints/";

Checkpoint::Checkpoint() {
    checkpoint_path = DEFAULT_CHECKPOINT_PATH;
}

void Checkpoint::addFolderInfo(const string &param, const string &value) {
    if (param=="path") {
        checkpoint_path = processCheckpointPath(value);
    }
}


void Checkpoint::addProperty(const string &elem, const string &param, const string &value) {
    if (elem == "folder") {
        addFolderInfo(param,value);
    }
}

const string &Checkpoint::getCheckpointPath() const {
    return checkpoint_path;
}

string Checkpoint::processCheckpointPath(string path) {
    if (path.back() != '/') {
        path.push_back('/');
    }
    return path;
}

bool Checkpoint::isCheckpointExists() const {
    struct stat info{};
    return (stat(checkpoint_path.c_str(), &info) == 0) && (info.st_mode & S_IFDIR) != 0;
}

void Checkpoint::saveCheckpoint() {

    string filename;
    vector<string> sub;
    string infilename = "test";
//    In[0]->split(infilename,'.',sub);
//    filename=sub[0].append(".");
//    filename = sub[0].append("_").append(".").append(name);
//    filename = In[0]->output_info.getOutputPath() + filename;
    ofstream outfile;
    outfile.open(filename+"."+Checkpoint::IN_CLASS_NAME, std::ios_base::app);

// Writing
//    outfile << MS_step << " ";
//    for ( auto n : distPerMC ) {
//        outfile << n << " ";
//    }
    outfile << "Hello Checkpoint" << endl;
}