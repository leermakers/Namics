
#include "mc_engine.h"
#include <iostream>

Metropolis_gradual::Metropolis_gradual() {
    std::cout << "Metropolis_gradual is created" << std::endl;
}

Metropolis_gradual::~Metropolis_gradual() {
    std::cout << "Metropolis_gradual destroyed" << std::endl;
}

void Metropolis_gradual::Action() {
    std::cout << "Action from Metropolis_gradual" << std::endl;
}
