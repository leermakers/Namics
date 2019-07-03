#include "stl_typedef.h"

template <typename T>
class Memento
{
  public:
    explicit Memento(stl::device_vector<T>& val_) noexcept
    {
        state = val_;
    }

    stl::device_vector<T> state;
};
