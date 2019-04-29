#include "../test_unit/catch.hpp"
#include "../namics.h"
#include <iomanip>      // std::setprecision

/**
 * Test 1
 */

TEST_CASE("Free energy after 10 steps", "[chain10steps], [cleng], [working]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_after10step.in";

// providing filename to testCaseCleng
    bool success = config.testCaseCleng(filename);
    if (!success) exit(0);

// some useful logic
    setprecision(7);          // set some precision
    Real f = -3.3460628362;   // expected number
    Real Free_energy = config.Cle[0]->free_energy_current;  // asking free_energy value from the engine
    REQUIRE((float) Free_energy == (float) f);
}

/**
 * Test 2
 */

TEST_CASE("Free energy vector after 20 steps", "[chain20steps_vector], [cleng], [working]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_after20step.in";

    const bool save_vector = true;  // if I would like to save something in my engine (optional)
    bool success = config.testCaseCleng(filename, save_vector);
    if (!success) exit(0);

// useful logic
    setprecision(7);  // set some precision
    int lenght = 19;  // 20 - 1
    vector<Real> e = {0.0};
    for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (std::size_t i = 0; i < v.size(); ++i) REQUIRE((float) e[i] == (float) v[i]);
}

/**
 * Test 3
 */

TEST_CASE("Free energy along 1 axis", "[chain_along1_axis], [cleng], [working], [long]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_go_along_axis1.in";

    const bool save_vector = true;  // if I would like to save something in my engine (optional)
    bool success = config.testCaseCleng(filename, save_vector);
    if (!success) exit(0);

// useful logic
    setprecision(7);  // set some precision
    int lenght = 199;  // 20 - 1
    vector<Real> e = {0.0};
    for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (std::size_t i = 0; i < v.size(); ++i) REQUIRE((float) e[i] == (float) v[i]);

}