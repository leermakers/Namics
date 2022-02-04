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

    void updateCheckpoint(const vector<std::shared_ptr<SimpleNode>> &simpleNodeList,
                          int MCS_done = 0, int mcs_done = 0, Real accepted = 0.0, Real rejected = 0.0);

    void saveCheckpoint(const vector<std::shared_ptr<SimpleNode>> &simpleNodeList,
                        int MCS_done = 0,
                        int mcs_done = 0,
                        Real accepted = 0.0,
                        Real rejected = 0.0);

    vector<shared_ptr<SimpleNode>> loadCheckpoint(vector<shared_ptr<SimpleNode>> simpleNodeList, const Point &box,
                                                  int& MCS_done,
                                                  int& mcs_done,
                                                  Real& accepted,
                                                  Real& rejected
                                                  );

    bool isCheckpointExists(const std::string &name) const;

    void updateCheckpointName(bool plusOne = false);

    bool isLoadable();

    static const string IN_CLASS_NAME;
    static const string DEFAULT_CHECKPOINT_PATH;
    static const string DEFAULT_CHECKPOINT_NAME;


private:

    string checkpoint_name;
    string checkpoint_path;
};

