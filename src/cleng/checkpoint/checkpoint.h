#ifndef NAMICS_CHECKPOINT_H
#define NAMICS_CHECKPOINT_H

#include "../../output_info.h"

class Checkpoint {

public:
    Checkpoint();

    void saveCheckpoint();
    void loadCheckpoint();

    const std::string &getCheckpointPath() const;
    bool isCheckpointExists() const;


private:

    std::string checkpoint_path;

};


#endif //NAMICS_CHECKPOINT_H
