#ifndef N_VECTOR_H_INCLUDED
#define N_VECTOR_H_INCLUDED

#include <string>
#include <sstream>
#include <ostream>
#include <cmath>

#include <boost/array.hpp>

#include <utl.h>

// this is designed for a vector of numbers of any type.
// dont get oversmart and define a matrix as a vector of vectors..
// you'll be surprised by the behavior of the arithmetic operators
template<typename T, std::size_t N,bool O= true>
class n_vector_t: public boost::array<T,N>
{

  public:

    typedef typename boost::array<T,N>              base_t;

    typedef typename base_t::value_type             value_type;
    typedef typename base_t::iterator               iterator;
    typedef typename base_t::const_iterator         const_iterator;
    typedef typename base_t::reference              reference;
    typedef typename base_t::const_reference        const_reference;
    typedef typename base_t::size_type              size_type;
    typedef typename base_t::difference_type        difference_type;
    typedef typename base_t::reverse_iterator       reverse_iterator;
    typedef typename base_t::const_reverse_iterator const_reverse_iterator;

    n_vector_t()
    {
    }

    template<typename OT,bool OO>
    n_vector_t( const n_vector_t<OT,N,OO> &o)
    {
      std::copy(o.begin(),o.end(),this->begin());
    }

    template<typename OT1,typename OT2>
    n_vector_t( const OT1 &e1, const OT2 &e2 )
    {
      ASSERT(N==2);

      (*this)[0] = e1;
      (*this)[1] = e2;
    }

    template<typename OT1,typename OT2,typename OT3>
    n_vector_t( const OT1 &e1,const OT2 &e2,const OT3 &e3 )
    {
      ASSERT(N==3);

      (*this)[0] = e1;
      (*this)[1] = e2;
      (*this)[2] = e3;
    }

    template<typename OT1,typename OT2,typename OT3,typename OT4>
    n_vector_t( const OT1 &e1,const OT2 &e2,const OT3 &e3,const OT4 &e4 )
    {
      ASSERT(N==4);

      (*this)[0] = e1;
      (*this)[1] = e2;
      (*this)[2] = e3;
      (*this)[3] = e4;
    }


    template<typename OT>
    static n_vector_t s_assign( const OT &o)
    {
      n_vector_t r;

      r.assign((T)o);

      return r;
    }

    template<typename OT>
    inline n_vector_t & operator+=(const OT &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] += o;

      return *this;
    }

    template<typename OT,bool OO>
    inline n_vector_t & operator+=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] += o[i];

      return *this;
    }

    template<typename OT>
    inline n_vector_t & operator-=(const OT &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] -= o;

      return *this;
    }

    template<typename OT,bool OO>
    inline n_vector_t & operator-=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] -= o[i];

      return *this;
    }

    template<typename OT>
    inline n_vector_t & operator*=(const OT &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] *= o;

      return *this;
    }

    template<typename OT,bool OO>
    inline n_vector_t & operator*=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] *= o[i];

      return *this;
    }

    template<typename OT>
    inline n_vector_t & operator/=(const OT &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] /= o;

      return *this;
    }

    template<typename OT,bool OO>
    inline n_vector_t & operator/=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] /= o[i];

      return *this;
    }

    template<typename OT>
    inline n_vector_t & operator%=(const OT &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] %= o;

      return *this;
    }

    template<typename OT,bool OO>
    inline n_vector_t & operator%=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] %= o[i];

      return *this;
    }

//    template<typename OT>
//    inline n_vector_t & operator&=(const OT &o)
//    {
//      for(size_t i = 0 ; i < N;++i )
//        (*this)[i] &= o;

//      return *this;
//    }

    template<typename OT,bool OO>
    inline n_vector_t & operator&=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] &= o[i];

      return *this;
    }

//    template<typename OT>
//    inline n_vector_t & operator|=(const OT &o)
//    {
//      for(size_t i = 0 ; i < N;++i )
//        (*this)[i] |= o;

//      return *this;
//    }

    template<typename OT,bool OO>
    inline n_vector_t & operator|=(const n_vector_t<OT,N,OO> &o)
    {
      for(size_t i = 0 ; i < N;++i )
        (*this)[i] |= o[i];

      return *this;
    }

    friend std::ostream &operator<< ( std::ostream &stream, const n_vector_t &e )
    {
      stream<<"(";
      if(N > 0)
      {
        stream<<e.elems[0];
        for(unsigned int i = 1 ; i < N ; ++i)
        {
          stream<<","<<e.elems[i];
        }
      }
      stream<<")";

      return stream;
    }

    friend std::istream &operator>>( std::istream &stream, n_vector_t &e )
    {
      char comma,bracket;
      stream>>bracket;
      if(N > 0)
      {
        stream>>e.elems[0];
        for(unsigned int i = 1 ; i < N ; ++i)
        {
          stream>>comma;
          stream>>e.elems[i];
        }
      }
      stream>>bracket;

      return stream;
    }

    static n_vector_t zero;

    static n_vector_t one;

//    typedef boost::function<T(T)> apply_t;

//    inline const n_vector_t apply(apply_t f)
//    {
//      std::transform(this->begin(),this->end(),this->begin(),f);

//      return (*this);
//    }

//    inline const n_vector_t apply(apply_t f) const
//    {
//      n_vector_t r;

//      std::transform(this->begin(),this->end(),r.begin(),f);

//      return (r);
//    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
      for(int i = 0 ; i < base_t::static_size; ++i) ar&(*this)[i];
    }

};

// addition
template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator+(const n_vector_t<T,N,O> & v, const OT &s)
{
  return n_vector_t<T,N,O>(v) += s;
}

template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator+(const OT &s,const n_vector_t<T,N,O> & v)
{
  return n_vector_t<T,N,O>(v) += s;
}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator+(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) += v2;
}

// multiplication
template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator*(const n_vector_t<T,N,O> & v, const OT &s)
{
  return n_vector_t<T,N,O>(v) *= s;
}

template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator*(const OT &s,const n_vector_t<T,N,O> & v)
{
  return n_vector_t<T,N,O>(v) *= s;
}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator*(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) *= v2;
}

// subtraction
template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator-(const n_vector_t<T,N,O> & v, const OT &s)
{
  return n_vector_t<T,N,O>(v) -= s;
}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator-(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) -= v2;
}

// division
template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator/(const n_vector_t<T,N,O> & v, const OT &s)
{
  return n_vector_t<T,N,O>(v) /= s;
}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator/(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) /= v2;
}

// modulus
template<typename T, std::size_t N,bool O,typename OT>
inline const n_vector_t<T,N,O> operator%(const n_vector_t<T,N,O> & v, const OT &s)
{
  return n_vector_t<T,N,O>(v) %= s;
}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator%(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) %= v2;
}

// binary and
//template<typename T, std::size_t N,bool O,typename OT>
//inline const n_vector_t<T,N,O> operator&(const n_vector_t<T,N,O> & v, const OT &s)
//{
//  return n_vector_t<T,N,O>(v) &= s;
//}

//template<typename T, std::size_t N,bool O,typename OT>
//inline const n_vector_t<T,N,O> operator&(const OT &s,const n_vector_t<T,N,O> & v)
//{
//  return n_vector_t<T,N,O>(v) &= s;
//}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator&(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) &= v2;
}

// binary or
//template<typename T, std::size_t N,bool O,typename OT>
//inline const n_vector_t<T,N,O> operator|(const n_vector_t<T,N,O> & v, const OT &s)
//{
//  return n_vector_t<T,N,O>(v) |= s;
//}

//template<typename T, std::size_t N,bool O,typename OT>
//inline const n_vector_t<T,N,O> operator|(const OT &s,const n_vector_t<T,N,O> & v)
//{
//  return n_vector_t<T,N,O>(v) |= s;
//}

template<typename T, std::size_t N,bool O,typename OT,bool OO>
inline const n_vector_t<T,N,O> operator|(const n_vector_t<T,N,O> & v1,const n_vector_t<OT,N,OO> & v2)
{
  return n_vector_t<T,N,O>(v1) |= v2;
}

// vector operations
template<typename T, std::size_t N>
inline T dot_product(const n_vector_t<T,N,true>& v1,const n_vector_t<T,N,true>& v2)
{
  T ret = 0;

  for(size_t i = 0 ; i < N;++i )
    ret += v1[i]*v2[i];

  return ret;
}
template<typename T>
inline n_vector_t<T,3,true> cross_product(const n_vector_t<T,3,true>& v1,
                                          const n_vector_t<T,3,true>& v2)
{
  n_vector_t<T,3,true> ret;

  for(size_t i = 0 ; i < 3;++i )
    ret[i] = v1[(i+1)%3]*v2[(i+2)%3] - v1[(i+2)%3]*v2[(i+1)%3];

  return ret;
}

template<typename T, std::size_t N>
inline T euclid_norm2(const n_vector_t<T,N,true>& v)
{
  return dot_product(v,v);
}

template<typename T, std::size_t N>
inline T euclid_norm(const n_vector_t<T,N,true>& v)
{
  return std::sqrt(euclid_norm2(v));
}

template<typename T, std::size_t N>
inline T euclid_distance
    (const n_vector_t<T,N,true>& v1,const n_vector_t<T,N,true>& v2)
{
  return euclid_norm(v1-v2);
}

template<typename T, std::size_t N>
inline T euclid_distance2
    (const n_vector_t<T,N,true>& v1,const n_vector_t<T,N,true>& v2)
{
  return euclid_norm2(v1-v2);
}

template<typename T, std::size_t N>
inline n_vector_t<T,N,true> euclid_normalize(const n_vector_t<T,N,true>& v)
{
  return v/euclid_norm(v);
}

template<typename T, std::size_t N,bool O>
n_vector_t<T,N,O> n_vector_t<T,N,O>::zero = n_vector_t<T,N,O>::s_assign(0);

template<typename T, std::size_t N,bool O>
n_vector_t<T,N,O> n_vector_t<T,N,O>::one  = n_vector_t<T,N,O>::s_assign(1);

namespace std
{
  template<typename T, std::size_t N,bool O>
  inline n_vector_t<T,N,O> ceil(const n_vector_t<T,N,O> & c)
  {
    n_vector_t<T,N,O> ret;

    for(int i = 0 ; i < N; ++i)
      ret[i] = ceil(c[i]);

    return ret;
  }

  template<typename T, std::size_t N,bool O>
  inline n_vector_t<T,N,O> floor(const n_vector_t<T,N,O> & c)
  {
    n_vector_t<T,N,O> ret;

    for(int i = 0 ; i < N; ++i)
      ret[i] = floor(c[i]);

    return ret;
  }

  template<typename T, std::size_t N,bool O>
  inline n_vector_t<T,N,O> abs(const n_vector_t<T,N,O> & c)
  {
    return c.apply((long double (*)(long double)) std::abs);
  }
}

template<class T, std::size_t N>
bool operator< (const n_vector_t<T,N,false>& x, const n_vector_t<T,N,false>& y)
{
  T ex[N], ey[N];

  std::copy ( x.begin(), x.end(), ex );
  std::copy ( y.begin(), y.end(), ey );

  std::sort ( ex, ex + N );
  std::sort ( ey, ey + N );

  return ( std::lexicographical_compare ( ex, ex + N, ey, ey + N ) );
}


template<class T, std::size_t N>
bool operator== (const n_vector_t<T,N,false>& x, const n_vector_t<T,N,false>& y)
{
  T ex[N], ey[N];

  std::copy ( x.begin(), x.end(), ex );
  std::copy ( y.begin(), y.end(), ey );

  std::sort ( ex, ex + N );
  std::sort ( ey, ey + N );

  return std::equal(ex, ex+N, ey);
}

#endif
