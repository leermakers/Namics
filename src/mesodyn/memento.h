#include "stl_typedef.h"

template <typename T>
class Memento
{
  public:
    Memento(stl::device_vector<T>& val_)
    {
        state = val_;
    }

    stl::device_vector<T> state;
};
