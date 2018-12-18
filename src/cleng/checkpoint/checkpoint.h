#pragma once
#include <string>
#include <vector>
#include <memory>
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>

#include <fstream>

using namespace std;

class Checkpoint {
public:
    Checkpoint();

    void saveCheckpoint(vector<std::shared_ptr<SimpleNode>> simpleNodeList);
    void loadCheckpoint(vector<std::shared_ptr<Node>> nodes, Point box);

    void addProperty(const string &elem, const string &param, const string &value);
    void addFolderInfo(const string &param, const string &value);

    bool isCheckpointExists(const std::string& name) const;

    void getLastCheckpoint();
    void getNewId4Checkpoint();

    static const string IN_CLASS_NAME;
    static const string DEFAULT_CHECKPOINT_PATH;
    static const string DEFAULT_CHECKPOINT_NAME;


private:

    string checkpoint_name;
    string checkpoint_path;
    string processCheckpointPath(std::string path);
};

