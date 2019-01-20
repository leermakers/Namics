#include "checkpoint.h"

using namespace std;

const string Checkpoint::IN_CLASS_NAME = ".checkpoint";
const string Checkpoint::DEFAULT_CHECKPOINT_PATH = "./checkpoints/";
const string Checkpoint::DEFAULT_CHECKPOINT_NAME = "0";

Checkpoint::Checkpoint() {
    checkpoint_path = DEFAULT_CHECKPOINT_PATH;
    checkpoint_name = DEFAULT_CHECKPOINT_NAME;
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

string Checkpoint::processCheckpointPath(string path) {
    if (path.back() != '/') {
        path.push_back('/');
    }
    return path;
}

bool Checkpoint::isCheckpointExists(const std::string& name) const {
    struct stat info{};
    return (stat(name.c_str(), &info) == 0);
}

void Checkpoint::getNewId4Checkpoint() {
    if (isCheckpointExists(checkpoint_path+checkpoint_name+IN_CLASS_NAME)) {
        int index = 0;
        while (isCheckpointExists(checkpoint_path+checkpoint_name+IN_CLASS_NAME)) {
            checkpoint_name = std::to_string(index);
            index ++;
        }
    }
}


void Checkpoint::getLastCheckpoint() {
    if (isCheckpointExists(checkpoint_path+checkpoint_name+IN_CLASS_NAME)) {
        int index = 0;
        while (isCheckpointExists(checkpoint_path+checkpoint_name+IN_CLASS_NAME)) {
            checkpoint_name = std::to_string(index);
            index ++;
        }
        index = index - 2;
        checkpoint_name = std::to_string(index);
    }
}


void Checkpoint::saveCheckpoint(vector<std::shared_ptr<SimpleNode>> simpleNodeList) {

    string filename;
    getNewId4Checkpoint();
    ofstream outfile;
    outfile.open(checkpoint_path+checkpoint_name+IN_CLASS_NAME, std::ios_base::app);

// Writing
    int index =0;
    for ( auto && n : simpleNodeList ) {

        if (!index){
            outfile << n->to_string();
            index ++;
        }
        else {
            outfile << n->to_string()<< endl;
            index = 0;
        }
        }
}

void Checkpoint::updateCheckpoint(vector<std::shared_ptr<SimpleNode>> simpleNodeList) {

    string filename;
    ofstream outfile;
    outfile.open(checkpoint_path+checkpoint_name+IN_CLASS_NAME);

// Writing
    int index =0;
    for ( auto && n : simpleNodeList ) {

        if (!index){
            outfile << n->to_string();
            index ++;
        }
        else {
            outfile << n->to_string()<< endl;
            index = 0;
        }
    }
}


vector<shared_ptr<Node>> createNodesFromFile(const vector<shared_ptr<SimpleNode>> &simple_nodes) {
    vector<shared_ptr<Node>> result;
    map<SimpleNode, vector<shared_ptr<SimpleNode>>> m;
    for (auto &&n  : simple_nodes) {
        m[*n].push_back(n);
    }
    for (auto &&entry : m) {
        if (entry.second.size() == 1) {
            result.push_back(entry.second[0]);
        } else {
            result.push_back(make_shared<Monolit>(entry.second));
        }
    }

    return result;
}

shared_ptr<SimpleNode> fromFileToNode(int x, int y, int z, int id, const Point &box) {
    return make_shared<SimpleNode>(Point(x, y, z), id, box);
}

vector<std::shared_ptr<Node>> Checkpoint::loadCheckpoint(vector<std::shared_ptr<Node>> nodes, Point box) {
    if (isCheckpointExists(checkpoint_path+checkpoint_name+IN_CLASS_NAME)) {
        cout << "Loading checkpoint ..." << endl;
        getLastCheckpoint();
        cout << "checkpoint: " << checkpoint_name+IN_CLASS_NAME << endl;
        ifstream infile(checkpoint_path+checkpoint_name+IN_CLASS_NAME);

        vector<shared_ptr<SimpleNode>> simpleNodeList;
        vector<Point> points;
        string line;
        cout << "Reading file..." << endl;
        while (getline(infile, line))
        {
            istringstream iss(line);
            int id1, x1, y1, z1, id2, x2, y2, z2;
            string a;
            if (!(iss >> a >> id1 >> a >> x1  >> a >> y1 >> a >> z1 >> a >> id2 >> a >> x2 >> a >> y2 >> a >> z2)) {
                cout << "Cant read!" << endl;
                break;
            }

            auto first_node  = fromFileToNode(x1, y1, z1, id1, box);
            auto second_node = fromFileToNode(x2, y2, z2, id2, box);
            first_node->set_cnode(second_node);
            second_node->set_cnode(first_node);
            //
            simpleNodeList.push_back(first_node);
            simpleNodeList.push_back(second_node);
        }

        nodes = createNodesFromFile(simpleNodeList);
        //
        cout << "Nodes: " << endl;
        for (auto &&n : nodes) {
            cout << n->to_string() << endl;
        }
        cout << "Success!" << endl;
    }
    else {
        cout << "Sorry! No file!" << endl;
    }
    return nodes;
}
