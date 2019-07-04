#ifndef MASK_H
#define MASK_H

#include <iterator>
#ifdef PAR_MESODYN
  #include <thrust/iterator/permutation_iterator.h>
#endif
#include "stl_typedef.h"

/* Provides a CPU implementation of thrust's permutation iterator */

template <class T>
struct Value_index_pair
{
    #ifdef PAR_MESODYN
    typedef typename thrust::device_vector<T>::iterator ElementIterator;
    typedef thrust::device_vector<size_t>::const_iterator IndexIterator;
    #endif

    stl::device_vector<T>& values;
    const stl::device_vector<size_t>& indices;

    Value_index_pair(stl::device_vector<T>& values_, const stl::device_vector<size_t>& indices_)
    : values{values_}, indices{indices_}
    { }

    // CPU implementation of thrust's permutation_iterator
    class iterator : public std::iterator<std::random_access_iterator_tag, T> {
      private:
        Value_index_pair *mask;
        int position{0};

      public:

        T& operator*() const {
            return mask->values[mask->indices[position]];
        }

        iterator(Value_index_pair* ptr = nullptr){mask = ptr;}
        iterator(const iterator& rawIterator) = default;
        ~iterator(){}

        iterator& operator=(const iterator& rawIterator) = default;
        iterator& operator=(Value_index_pair* ptr){mask = ptr;return (*this);}
      
        iterator& operator++()    /* prefix */              { ++position; return *this;}
        iterator& operator--()    /* prefix */              { --position; return *this;}
        iterator& operator+=(const int& change)             { position += change;return (*this);}
        iterator& operator-=(const int& change)             { position -= change;return (*this);}
        iterator  operator++(int) /* postfix */             { auto temp(*this);++position;return temp;}
        iterator  operator--(int) /* postfix */             { auto temp(*this);--position;return temp;}
        const iterator  operator+(const long int& change)   { auto prev_pos = position;position+=change;auto temp(*this);position = prev_pos;return temp;}
        iterator  operator+(const iterator& Iterator)       { auto prev_pos = position;position+=Iterator.get_const_pos();auto temp(*this);position = prev_pos;return temp;}
        iterator  operator-(const long int& change)         { auto prev_pos = position;position-=change;auto temp(*this);position = prev_pos;return temp;}
        //iterator  operator-(iterator& Iterator)             { auto prev_pos = position;position-=Iterator.get_const_pos();auto temp(*this);position = prev_pos;return temp;}
        iterator  operator-(const iterator& Iterator)       { auto dif_pos = position-Iterator.get_const_pos();return dif_pos;}
        long int  operator-(iterator& Iterator)             { auto dif_pos = position-Iterator.get_const_pos();return dif_pos;}

        bool      operator==(const iterator& Iterator)const { return ( get_const_pos() == Iterator.get_const_pos() );}
        bool      operator!=(const iterator& Iterator)const { return ( get_const_pos() != Iterator.get_const_pos() );}
        bool      operator>(const iterator& Iterator)const  { return ( get_const_pos() > Iterator.get_const_pos()  );}
        bool      operator>=(const iterator& Iterator)const  { return ( get_const_pos() >= Iterator.get_const_pos()  );}
        bool      operator>=(const size_t& position)const  { return ( get_const_pos() >= position  );}
        bool      operator<=(const iterator& Iterator)const  { return ( get_const_pos() <= Iterator.get_const_pos()  );}
        bool      operator<(const iterator& Iterator)const  { return ( get_const_pos() < Iterator.get_const_pos()  );}
        //TODO: other compare operators

         int get_const_pos() const {return position;}

         iterator begin() {position=0; return *this;}
         iterator end() {position=mask->indices.size();return *this;}
    };

    #ifdef PAR_MESODYN
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

    #ifdef PAR_MESODYN  
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