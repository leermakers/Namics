
#include "../../namics.h"
#include "../nodes/simple_node.h"
#include "../nodes/monolit.h"
#include "../nodes/point.h"
#include "../random/random.h"
#include <map>
#include <vector>
#include <memory>

using namespace std;

class Analyzer {
private:
    class Cube {
    public:

        Cube() = default;
        map<int, vector<Point>> layer_point_map;
        void construct_cube(int requested_layers, const Point& core);
        static vector<Point> get_cube_odd(int layers, const Point& core, int layer);
        map<int, vector<Point>> get_layer_point_map() const;
    };

public:
    Analyzer(string filename, int requested_layers, const Point& core, Real accepted, Real rejected) {
        this->filename = filename;
        c.construct_cube(requested_layers, core);
        this->accepted = accepted;
        this->rejected = rejected;
    }
    Analyzer () = default;

    string metropolis_name = "[Metropolis] ";
    string outerinner_name = "[Outer/Inner] ";
    string internal_name   = "[Analysis] ";

    string filename;

    Cube c;
    map<string, Point> pointsFromVtk;
    map<int, vector<Point>> layer_points_map;

    Real accepted = 0.0;
    Real rejected = 0.0;
    //

    void updateVtk2PointsRepresentation(const vector<Real>& vtk, const Point& box);
    static map<string, Point> convertVtk2Points(const vector<Real>& vtk, const Point& box);
    map<int, vector<Point>> convertPoints2LayerPoints(const map<string, Point>& points4converting) const;
    static Real calculateRe(map<int, vector<int>> pivot_arm_nodes, map<int, vector<std::shared_ptr<Node>>> nodes_map);
    Real calculateRg();

    bool Metropolis(Random& rand,const Real& prefactor_kT,Real& free_energy_trial,Real& free_energy_current, bool update_rate = true);
    
    bool OuterInner(Random& rand,const Real& prefactor_kT,Real& free_energy_trial,Real& free_energy_current);

    void notification(int MC_attempt, int cleng_rejected, int mcs_done = 0) const;

    void Write2File(const string& what, const vector<int>& values, bool same_line) const;

    map<int, vector<Real>> calculatePhi() const;
    map<int, vector<Point>> get_layer_point_map() const {return c.get_layer_point_map();};
};
