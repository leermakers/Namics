#include "cleng.h"
#ifndef NAMICS_MOVES_H
#define NAMICS_MOVES_H

bool Cleng::Checks(int id_node_for_move) {
    bool in_range;

    bool in_subbox_range = InSubBoxRange(id_node_for_move);

    bool not_collapsing = NotCollapsing(id_node_for_move);

    bool consistent_monoliths = true;

    // check distance between all nodes_map and constrains (walls)
    // BC.x, BC.y, BC.z = mirror
    if (BC.x and BC.y and BC.z) in_range = InRange(id_node_for_move);
        // BC.x and/or BC.y and/or BC.z != mirror
    else in_range = true;

    // check distances between all nodes_map => checking possibility to achieve clamped nodes
    bool commensurate = IsCommensuratable();

    for (auto &&node_id_value : nodes_map) {
        consistent_monoliths = node_id_value.second.data()->get()->_checkPoints();
    }

    bool result = not_collapsing and in_range and in_subbox_range and commensurate and consistent_monoliths;
//    cout << "[Checks] result: " << result << endl;
    return result;
}

bool Cleng::MakeChecks(int id_node_for_move) {
    bool success = true;
    while (!Checks(id_node_for_move)) {
        cout << internal_name << "Prepared MC step for the node ["<< id_node_for_move << "] does not pass checks. Rejected." << endl;
        cleng_rejected++;
        success = false;
        return success; // TODO: CHECK THIS RETURN
    }
    return success;
}

Point Cleng::preparePivotClampedMove(int id_node_for_move) {
    Point center_of_rotation = nodes_map[pivot_node_ids[0]].data()->get()->point();
    // moving center of rotation to the center
    for (auto &&node : nodes_map) node.second.data()->get()->shift(center_of_rotation.negate());
    Point current_Point = nodes_map[id_node_for_move].data()->get()->_returnSystemPoint();
    Point point_shifted_by_matrix = rotation_matrix.dot(current_Point);
    Point clamped_move = point_shifted_by_matrix - current_Point;
    // center +1
    if (clamped_move.x < -(box.x / 2) +1) clamped_move.x += box.x;
    if (clamped_move.y < -(box.y / 2) +1) clamped_move.y += box.y;
    if (clamped_move.z < -(box.z / 2) +1) clamped_move.z += box.z;

    if (clamped_move.x > (box.x / 2)+1) clamped_move.x -= box.x;
    if (clamped_move.y > (box.y / 2)+1) clamped_move.y -= box.y;
    if (clamped_move.z > (box.z / 2)+1) clamped_move.z -= box.z;
    // moving center of rotation to initial place
    for (auto &&node : nodes_map) node.second.data()->get()->shift(center_of_rotation);
    return clamped_move;
}

void Cleng::_moveClampedNode(bool back, int id_node_for_move, const Point& clamped_move) {
    if (back) {
        //Point _clamped_move = clamped_move.negate();
        if (!simultaneously) {
            //cout << "[BACK][Prepared node] id: " << id_node_for_move << " clamped_move: " << clamped_move.negate().to_string() << "... " << endl;
            nodes_map[id_node_for_move].data()->get()->shift(clamped_move.negate());
        }
        else for (auto &node : nodes_map) node.second.data()->get()->shift(clamped_move.negate());
    } else {
        if (!simultaneously) {
            cout << "[Prepared node] id: " << id_node_for_move << " clamped_move: " << clamped_move.to_string() << "... ";
            nodes_map[id_node_for_move].data()->get()->shift(clamped_move);
            cout << "[Moved]" << endl;
        } else {
            // simultaneously procedure (+ two_ends_extension)
            if (two_ends_extension) {
                // TODO: deprecate in future...
                for (auto &&node : nodes_map) {
                    int index = node.first;
                    if (index % 2 == 0) node.second.data()->get()->shift(clamped_move.negate());
                    else node.second.data()->get()->shift(clamped_move);
                }
                cout << "[Moved] " << "*All* "
                     << "MC step: " << clamped_move.to_string() << " and " << clamped_move.negate().to_string() << endl;
            } else {
                for (auto &&node : nodes_map) node.second.data()->get()->shift(clamped_move);
                cout << "[Moved] " << "*All* " << "MC step: " << clamped_move.to_string() << endl;
            }
        }
    }
}

bool Cleng::_pivotMoveClampedNode(const bool &back) {
    Point clamped_move;
    bool success = true;
    cout << internal_name << "[pivot_move]" << endl;
    prepareMove("pivot_move");
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        clamped_move = preparePivotClampedMove(pivot_node_ids[node_position_in_vector]);  // need to create for each node_id!
        _moveClampedNode(back, pivot_node_ids[node_position_in_vector], clamped_move);
        nodeIDs_clampedMove[pivot_node_ids[node_position_in_vector]] = clamped_move;
    }
    cout << internal_name << "checking positions...";
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        cout << "\nNode_id: " << pivot_node_ids[node_position_in_vector] << endl;
        success = MakeChecks(pivot_node_ids[node_position_in_vector]);
        if (!success) return success;
    }
    cout << " OK" << endl;
    return  success;
}

bool Cleng::_pivotMoveOneBond(const bool &back) {
    Point clamped_move;
    bool success = true;
    cout << internal_name << "[pivot_one_bond]" << endl;
    clamped_move = prepareMove("pivot_one_bond_move");  // need to create at once and move all
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        _moveClampedNode(back, pivot_node_ids[node_position_in_vector], clamped_move);
        nodeIDs_clampedMove[pivot_node_ids[node_position_in_vector]] = clamped_move;
    }
    cout << internal_name << "checking positions...";
    for (size_t node_position_in_vector=1; node_position_in_vector != pivot_node_ids.size(); node_position_in_vector++) {
        success = MakeChecks(pivot_node_ids[node_position_in_vector]);
        if (!success) return success;
    }
    cout << " OK" << endl;
    return  success;
}

bool Cleng::_oneNodeMoveClampedNode(const bool &back) {
    int id_node_for_move;
    Point clamped_move;
    bool success;

    cout << internal_name << "[one_node_move]" << endl;
    id_node_for_move = prepareIdNode();
    clamped_move = prepareMove("one_node_move");
    _moveClampedNode(back, id_node_for_move, clamped_move);
    nodeIDs_clampedMove[id_node_for_move] = clamped_move;

    cout << internal_name << "checking positions...";
    success = MakeChecks(id_node_for_move);
    if (!success) return success;
    nodeIDs_clampedMove[id_node_for_move] = clamped_move;
    cout << " OK" << endl;
    return success;
}


#endif //NAMICS_MOVES_H
