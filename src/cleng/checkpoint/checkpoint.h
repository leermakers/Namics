#pragma once
#include <string>

using namespace std;

class Checkpoint {
public:
    Checkpoint();

    void saveCheckpoint();
    void loadCheckpoint();

    void addProperty(const string &elem, const string &param, const string &value);
    void addFolderInfo(const string &param, const string &value);

    const string &getCheckpointPath() const;
    bool isCheckpointExists() const;

    static const string IN_CLASS_NAME;
    static const string DEFAULT_CHECKPOINT_PATH;

private:

    string checkpoint_path;
    string processCheckpointPath(std::string path);

};

