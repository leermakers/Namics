#pragma once

class mc_engine_base
{
public:

    // The "Virtual Constructor"
    static mc_engine_base *Create(int id);

    mc_engine_base() = default;

    virtual // Ensures to invoke actual object destructor
    ~mc_engine_base() = default;

    // An interface
    virtual void Action() = 0;
};

class Metropolis : public mc_engine_base
{
public:
    Metropolis();
    ~Metropolis() override;
    void Action() override;
};

class Metropolis_gradual : public mc_engine_base
{
public:
    Metropolis_gradual();
    ~Metropolis_gradual() override;
    void Action() override;
};

class MC_engine
{
public:
    explicit MC_engine(int input = 1) : ptrmcEngineBase(nullptr)
    {
        ptrmcEngineBase = mc_engine_base::Create(input);
    }

    ~MC_engine()
    {
        if( ptrmcEngineBase )
        {
            delete ptrmcEngineBase;
            ptrmcEngineBase = nullptr;
        }
    }

    // Delegates to actual object
    void Action()
    {
        ptrmcEngineBase->Action();
    }

private:
    mc_engine_base *ptrmcEngineBase;
};

