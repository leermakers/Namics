#ifndef NAMICS_CHECKPOINT_H
#define NAMICS_CHECKPOINT_H

#include "../../output_info.h"

class Checkpoint {

public:
    Checkpoint();

    void saveCheckpoint();
    void loadCheckpoint();

    void addProperty(const std::string &elem, const std::string &param, const std::string &value);
    void addFolderInfo(const std::string &param, const std::string &value);

    const std::string &getCheckpointPath() const;
    bool isCheckpointExists() const;

    static const std::string IN_CLASS_NAME;
    static const std::string DEFAULT_OUTPUT_PATH;

private:

    std::string checkpoint_path;
    std::string processCheckpointPath(std::string path);

};


#endif //NAMICS_CHECKPOINT_H
