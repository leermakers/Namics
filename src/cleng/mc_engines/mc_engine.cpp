#include "mc_engine.h"

mc_engine_base *mc_engine_base::Create(int id)
{
    if( id == 1 )
    {
        return new Metropolis;
    }
    else if( id == 2 )
    {
        return new Metropolis_gradual;
    }
    else
    {
        return nullptr;
    }
}
