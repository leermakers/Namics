#ifndef MASK_H
#define MASK_H

#include <iterator>
#ifdef PAR_MESODYN_THRUST
  #include <thrust/iterator/permutation_iterator.h>
#endif
#include "stl_typedef.h"

/* Provides a CPU implementation of thrust's permutation iterator */

template <class T>
struct Value_index_pair
{
    #ifdef PAR_MESODYN_THRUST
    typedef typename thrust::device_vector<T>::iterator ElementIterator;
    typedef thrust::device_vector<size_t>::const_iterator IndexIterator;
    #endif

    stl::device_vector<T>& values;
    const stl::device_vector<size_t>& indices;

    Value_index_pair(stl::device_vector<T>& values_, const stl::device_vector<size_t>& indices_)
    : values{values_}, indices{indices_}
    { }

    // CPU implementation of thrust's permutation_iterator
    class iterator {
      private:
        Value_index_pair *mask;
        int position{0};

      public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        T& operator*() const {
            return mask->values[mask->indices[position]];
        }

        T& operator[](difference_type n) const {
            return mask->values[mask->indices[position + n]];
        }

        iterator(Value_index_pair* ptr = nullptr){mask = ptr;}
        iterator(const iterator& rawIterator) = default;
        virtual ~iterator(){}

        iterator& operator=(const iterator& rawIterator) = default;
        iterator& operator=(Value_index_pair* ptr){mask = ptr;return (*this);}

        iterator& operator++()                                { ++position; return *this;}
        iterator& operator--()                                { --position; return *this;}
        iterator& operator+=(difference_type change)          { position += change;return (*this);}
        iterator& operator-=(difference_type change)          { position -= change;return (*this);}
        iterator  operator++(int)                             { auto temp(*this);++position;return temp;}
        iterator  operator--(int)                             { auto temp(*this);--position;return temp;}
        iterator  operator+(difference_type change) const     { iterator temp(*this);temp.position+=change;return temp;}
        iterator  operator-(difference_type change) const     { iterator temp(*this);temp.position-=change;return temp;}
        difference_type operator-(const iterator& other) const { return position - other.position;}

        friend iterator operator+(difference_type n, const iterator& it) { return it + n; }

        bool      operator==(const iterator& Iterator)const { return ( get_const_pos() == Iterator.get_const_pos() );}
        bool      operator!=(const iterator& Iterator)const { return ( get_const_pos() != Iterator.get_const_pos() );}
        bool      operator>(const iterator& Iterator)const  { return ( get_const_pos() > Iterator.get_const_pos()  );}
        bool      operator>=(const iterator& Iterator)const { return ( get_const_pos() >= Iterator.get_const_pos()  );}
        bool      operator<=(const iterator& Iterator)const { return ( get_const_pos() <= Iterator.get_const_pos()  );}
        bool      operator<(const iterator& Iterator)const  { return ( get_const_pos() < Iterator.get_const_pos()  );}

         int get_const_pos() const {return position;}

         iterator begin() {position=0; return *this;}
         iterator end() {position=mask->indices.size();return *this;}
    };

    #ifdef PAR_MESODYN_THRUST
    thrust::permutation_iterator<ElementIterator,IndexIterator> begin() {
      thrust::permutation_iterator<ElementIterator,IndexIterator> itt(values.begin(), indices.begin());
      return itt;
    }
    #else
    iterator begin() {
      Value_index_pair::iterator itt(this);
      return itt;
    }
    #endif

    #ifdef PAR_MESODYN_THRUST  
    thrust::permutation_iterator<ElementIterator,IndexIterator> end() {
      thrust::permutation_iterator<ElementIterator,IndexIterator> itt(values.end(), indices.end());
      return itt;
    }
    #else
    iterator end() {
      Value_index_pair::iterator itt(this);
      return itt.end();
    }
    #endif
};

#endif