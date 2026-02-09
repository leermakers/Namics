#ifndef STL_TYPEDEF_H
#define STL_TYPEDEF_H

  /*
   *  This allows you to write code that works both with Thrust (CUDA STL implementation) and the STL.
   *  GPU vectors are declared as stl::device_vector<datatype>, host vectors are declared as stl::host_vector<datatype>
   *  Both of these will be transformed into a std::vector<datatype> when PAR_MESODYN_THRUST has not been passed to the compiler (CPU scenario)
   */

  #ifdef PAR_MESODYN_THRUST
  #define DEVICE_LAMBDA __host__ __device__
    #include <thrust/device_vector.h>
    #include <thrust/host_vector.h>
    #include <thrust/device_ptr.h>
    namespace stl = thrust;
    //const auto reduce = accumulate;
  #else
    #define DEVICE_LAMBDA
    #include <vector>
    namespace stl = std;
    namespace std {
      template<typename T> using host_vector = vector<T>;   
      template<typename T> using device_vector = vector<T>;   
    }
  #endif

  // C++17 parallel algorithms (CPU parallelism via TBB)
  // enable with: make PAR_MESODYN_STL=1
  #if defined(PAR_MESODYN_STL) && !defined(PAR_MESODYN_THRUST)
    #include <execution>
    #define EXEC_PAR std::execution::par,
    #define EXEC_SEQ std::execution::seq,
    #define EXEC_PAR_UNSEQ std::execution::par_unseq,
  #endif

  #ifndef EXEC_PAR
    #define EXEC_PAR
    #define EXEC_SEQ
    #define EXEC_PAR_UNSEQ
  #endif

#endif