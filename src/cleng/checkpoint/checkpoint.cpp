#include "checkpoint.h"

const std::string Checkpoint::IN_CLASS_NAME = "checkpoint";
const std::string Checkpoint::DEFAULT_OUTPUT_PATH = "./checkpoints/";

//Checkpoint::Checkpoint() {
//    checkpoint_path = DEFAULT_OUTPUT_PATH;
//}
//
//void Checkpoint::addProperty(const std::string &elem, const std::string &param, const std::string &value) {
//    if (elem == "folder") {
//        addFolderInfo(param,value);
//    }
//}
//
//void Checkpoint::addFolderInfo(const std::string &param, const std::string &value) {
//    if (param=="path") {
//        checkpoint_path = processCheckpointPath(value);
//    }
//}
