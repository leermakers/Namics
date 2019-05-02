#include "catch.hpp"
#include "../namics.h"

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
    Real threshold = 1e-5; // set some precision
    Real f = -3.3460628362;   // expected number
    Real Free_energy = config.Cle[0]->free_energy_current;  // asking free_energy value from the engine
    REQUIRE( Free_energy - f  < threshold );
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
    Real threshold = 1e-5; // set some precision
    vector<Real> e;
    int lenght = 20;
    for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (std::size_t i = 0; i < v.size(); ++i) REQUIRE( abs(e[i] - v[i]) < threshold );
}

/**
 * Test 3
 */

TEST_CASE("Free energy of chain along 1 axis", "[chain_along1_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
    NamicsConfig config;
    string filename = "chain_go_along_axis1.in";

    const bool save_vector = true;  // if I would like to save something in my engine (optional)
    bool success = config.testCaseCleng(filename, save_vector);
    if (!success) exit(0);

// useful logic
    Real threshold = 1e-5;          // set some precision
    vector<Real> e;
    int lenght = 200;
    for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
    vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
    for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 4
 */

TEST_CASE("Free energy of chain along 2 axis", "[chain_along2_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_go_along_axis2.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 200;
for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 5
 */

TEST_CASE("Free energy of chain along 3 axis", "[chain_along3_axis], [cleng], [working], [long], [chain_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_go_along_axis3.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 200;
for (int i = 0; i < lenght; i++) { e.push_back(-3.3460628362); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 6
 */

TEST_CASE("Free energy of 3 segments along 1 axis", "[3segments_along1_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis1.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 200;
for (int i = 0; i < lenght; i++) { e.push_back(-5.49114); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 6
 */

TEST_CASE("Free energy of 3 segments along 2 axis", "[3segments_along2_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis2.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 200;
for (int i = 0; i < lenght; i++) { e.push_back(-5.49114); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 7
 */

TEST_CASE("Free energy of 3 segments along 3 axis", "[3segments_along3_axis], [cleng], [working], [long], [3segments_along_axises]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "3segments_go_along_axis3.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// useful logic
Real threshold = 1e-5;          // set some precision
vector<Real> e;
int lenght = 200;
for (int i = 0; i < lenght; i++) { e.push_back(-5.49114); } // expected vector
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}


/**
 * Test 8
 */

TEST_CASE("Free energy of chain extension ", "[chain_extension], [cleng], [working], [long]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_extension.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);
// useful logic
Real threshold = 1e-5; // set some precision
vector<Real> e = {     //
        -3.14279,
        -3.06371,
        -2.98835,
        -2.91115,
        -2.85009,
        -2.79798,
        -2.76777
};
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = 0; i < v.size(); ++i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}

/**
 * Test 9
 */

TEST_CASE("Free energy of chain compression ", "[chain_compression], [cleng], [working], [long]") {

// create a configuration of Namics
NamicsConfig config;
string filename = "chain_compression.in";

const bool save_vector = true;  // if I would like to save something in my engine (optional)
bool success = config.testCaseCleng(filename, save_vector);
if (!success) exit(0);

// What Namics gives:
// -2.79685
// -2.8497
// -2.91376
// -2.98672
// -3.06276
// -3.14493
// -3.21317

// useful logic
Real threshold = 1e-5; // set some precision
vector<Real> e = {     // the same vector as before but comparison should be made reversely.
        -3.14279,
        -3.06371,
        -2.98835,
        -2.91115,
        -2.85009,
        -2.79798,
        -2.76777
};
vector<Real> v = config.Cle[0]->test_vector; // if save_vector was provided I can ask variable (optional)
for (size_t i = v.size(); i > 0; --i) { REQUIRE( abs(e[i] - v[i]) < threshold);}
}