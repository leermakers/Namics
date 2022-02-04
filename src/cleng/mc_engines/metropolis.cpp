
#include "mc_engine.h"
#include <iostream>

Metropolis::Metropolis() {
    std::cout << "Metropolis created" << std::endl;
}

Metropolis::~Metropolis() {
    std::cout << "Metropolis destroyed" << std::endl;
}

void Metropolis::Action() {
    std::cout << "Action from Metropolis" << std::endl;
}
