/*=========================================================================

  Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
  All rights reserved.

=========================================================================*/

#ifndef N88UTIL_array_hpp_INCLUDED
#define N88UTIL_array_hpp_INCLUDED

#include "Tuple.hpp"
#include "exception.hpp"
#include <iostream>
#include <cstdlib>
#include <cstring>

#ifdef TRACK_ALLOCATIONS
#include "TrackingAllocator.hpp"
#endif


// Alignment in bytes
#define N88_ARRAY_ALIGNMENT_POWER 4
#define N88_ARRAY_ALIGNMENT (1 << N88_ARRAY_ALIGNMENT_POWER)
#define N88_ARRAY_ALIGNMENT_MASK (N88_ARRAY_ALIGNMENT-1)



namespace n88
  {

  /**
    * A base class for array.
    *
    * The purpose of this base class is to provide a bunch of common methods
    * that can be inherited in partial specializations of array.  Please
    * refer to array for more details.
    */
  template <int N, typename TValue, typename TIndex=size_t>
  class array_base
    {

    protected:

      TValue*         m_base;
      size_t          m_size;   // Intentionally size_t and not TIndex,
                                // as it may be desireable on occasion
                                // to create an array with dims in
                                // a smaller type, whose product however
                                // exceeds the representation of the smaller
                                // type.
      TValue*         m_buffer;
      TValue*         m_end;
      Tuple<N,TIndex> m_dims;

    public:

      enum {dimension = N};
      typedef TValue value_type;
      typedef TIndex index_type;

      /** Empty constructor.
        * You must subsequently call construct or construct_reference explicitly.
        */
      array_base()
        :
        m_base             (NULL),
        m_size             (0),
        m_buffer           (NULL),
        m_end              (NULL),
        m_dims             (Tuple<N,TIndex>::zeros())
        {}

      /** Constructor to allocate space.
        *
        * @param dims  The dimensions of the array to allocate.
        */
      explicit array_base(Tuple<N,TIndex> dims)
        :
        m_base             (NULL),
        m_size             (0),
        m_buffer           (NULL),
        m_end              (NULL),
        m_dims             (Tuple<N,TIndex>::zeros())
        { construct(dims); }

      /** Constructor to create reference to existing data defined by a pointer.
        * The referenced data will never be freed by this object.  You
        * must ensure that the referenced memory remains valid, and perform
        * appropriate clean up outside this object.
        *
        * @param data  A pointer to the data.
        * @param dims  The dimensions of the array.
        */
      explicit array_base(TValue* data, Tuple<N,TIndex> dims)
        :
        m_base             (data),
        m_size             (long_product(dims)),
        m_buffer           (NULL),
        m_end              (data + long_product(dims)),
        m_dims             (dims)
        {}

      /** Constructor to create reference to existing data in array object.
        * No automatic reference counting is performed.  You
        * must ensure that the referenced array remains in existence for as
        * long as you want to use this reference.
        *
        * @param source  An existing array object.
        */
      explicit array_base(const array_base<N,TValue,TIndex>& source)
        :
        m_base             (source.ptr()),
        m_size             (source.size()),
        m_buffer           (NULL),
        m_end              (source.end()),
        m_dims             (source.dims())
        {}

      ~array_base()
        { this->destruct(); }

      /** Allocate space.
        *
        * @param dims  The dimensions of the array to allocate.
        */
      void construct(Tuple<N,TIndex> dims)
        {
        if (this->m_base)
          { throw_n88_exception("array is already constructed."); }
        this->m_size = long_product(dims);
        this->m_dims = dims;
        try
#ifdef TRACK_ALLOCATIONS
          { this->m_buffer = (TValue*)TrackingAllocator::allocate(this->m_size*sizeof(TValue)); }
#else
          { this->m_buffer = new TValue[this->m_size]; }
#endif
        catch (...)
          { throw_n88_exception("Unable to allocate memory."); }
        if (this->m_buffer == NULL)
          { throw_n88_exception("Unable to allocate memory."); }
        // In OS's that use lazy allocation, calling memset may ensure that
        // memory is contiguous in real address space, which could be advantageous.
        memset (this->m_buffer, 0, this->m_size*sizeof(TValue));
        this->m_base = this->m_buffer;
        // For now, just throw an exception if the alignment is not correct.
        // If this crops up, will need to re-implement alignment.
        // Note that malloc is only guaranteed to 8 byte alignment, but in
        // practice on many systems is actually 16 byte aligned.
        if (size_t(this->m_base) & N88_ARRAY_ALIGNMENT_MASK)
          { throw_n88_exception("array alignment error"); }
        this->m_end = this->m_base + this->m_size;
        }

      /** Allocate space if not already done.
        *
        * No guarantee that data is zeroed if previously allocated.
        *
        * Throws an exception if already allocated with different dims.
        *
        * @param dims  The dimensions of the array to allocate.
        */
      void construct_if_required(Tuple<N,TIndex> dims)
        {
        if (this->m_base)
          { n88_assert (dims == this->m_dims); }
        else
          { this->construct(dims); }
        }

      /** Allocate space if not already done.
        *
        * Guarantees that data is zeroed.
        *
        * Throws an exception if already allocated with different dims.
        *
        * @param dims  The dimensions of the array to allocate.
        */
      void construct_or_zero(Tuple<N,TIndex> dims)
        {
        if (this->m_base)
          {
          n88_assert (dims == this->m_dims);
          this->zero();
          }
        else
          { this->construct(dims); }
        }

      /** Create a reference to existing data defined by a pointer.
        * The referenced data will never be freed by this object.  You
        * must ensure that the referenced memory remains valid, and perform
        * appropriate clean up outside this object.
        *
        * @param data  A pointer to the data.
        * @param dims  The dimensions of the array.
        */
      void construct_reference(TValue* data, Tuple<N,TIndex> dims)
        {
        if (this->m_base)
          { throw_n88_exception("array is already constructed."); }
        this->m_size = long_product(dims);
        this->m_dims = dims;
        this->m_buffer = NULL;
        this->m_base = data;
        this->m_end = this->m_base + this->m_size;
        }

      /** Create a reference to to existing data in array object.
        * No automatic reference counting is performed.  You
        * must ensure that the referenced array remains in existence for as
        * long as you want to use this reference.
        *
        * @param source  An existing array object.
        */
      void construct_reference(const array_base<N,TValue,TIndex>& source)
        { this->construct_reference (source.ptr(), source.dims()); }

      void destruct()
        {
        if (this->m_buffer)
#ifdef TRACK_ALLOCATIONS
          { TrackingAllocator::release (this->m_buffer, this->m_size*sizeof(TValue)); }
#else
          { delete[] this->m_buffer; }
#endif
        this->m_size = 0;
        this->m_dims = Tuple<N,TIndex>::zeros();
        this->m_buffer = NULL;
        this->m_base = NULL;
        this->m_end = NULL;
        }

      /** Returns true if the array has been constructed. */
      inline bool is_constructed() const
        { return (this->m_base != 0); }

      /** Returns the flattened (1D) size of the array. */
      inline size_t size() const
        { return this->m_size; }

      /** Returns the dimensions of the array. */
      inline Tuple<N,TIndex> dims() const
        { return this->m_dims; }

      /** Returns a pointer to the array data. */
      inline TValue* ptr() const
        {
#ifdef RANGE_CHECKING
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
#endif
        return this->m_base;
        }

      /** A utility function that performs a run-time check of a pointer to
        * ensure that it points to an entry in the array data.  Note
        * that high-performance code often must use raw pointers; this provides
        * a mechanism to do some verification and debugging of such code.
        * This check is performed only if compiled with RANGE_CHECKING
        * defined; naturally there is a substantial performance penalty.
        */
      inline const TValue* verify_ptr(const TValue* p) const
        {
#ifdef RANGE_CHECKING
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        if ((p < this->m_base) || (p >= this->m_end))
          { throw_n88_exception("array index out of bounds."); }
        if (((p - this->m_base) % sizeof(TValue)) != 0)
          { throw_n88_exception("array pointer has incorrect offset."); }
#endif
        return p;
        }
      inline TValue* verify_ptr(TValue* p) const
        { return const_cast<TValue*>(this->verify_ptr(static_cast<const TValue*>(p))); }

      /** Pointer to the last element plus one of the array data.  This
        * may be used in loops various STL algorithms that require an end
        * value.
        */
      inline TValue* end() const
        {
#ifdef RANGE_CHECKING
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
#endif
        return this->m_end;
        }

      /** Converts a tuple index to the flattened 1D equivalent index. */
      inline size_t flat_index(Tuple<N,TIndex> indices) const
        {
        size_t index = indices[0];
        for (int i=1; i<N; i++)
          {
          index = index*static_cast<size_t>(this->m_dims[i]) + static_cast<size_t>(indices[i]);
          }
        return index;
        }

      /** Flat (1D) indexing of the array data. */
      inline TValue& operator[](size_t i) const
        {
#ifdef RANGE_CHECKING
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        if (i >= this->m_size)
          { throw_n88_exception("array index out of bounds."); }
#endif
        return this->m_base[i];
        }

      /** Indexing of the array data using N-dimensional tuples. */
      inline TValue& operator()(Tuple<N,TIndex> indices) const
        { return array_base::operator[](this->flat_index(indices)); }

      /** Set all array data to zero.
        * Note that this is done on constructing the data (for performance reasons),
        * so calling this immediately after is redundant.
        */
      inline void zero() const
        {
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        memset (this->m_base, 0, sizeof(TValue)*this->m_size);
        }

      /** Copy data from an existing array.
        * This array must either have existing dimensions the same as the
        * source array, or it can be unconstructed, in which case it will
        * be constructed with the same dimensions.
        */
      inline void copy(const array_base& rhs) const
        {
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        else if (this->m_dims != rhs.dims())
          { throw_n88_exception("cannot copy different sized arrays."); }
        memcpy (this->m_base, rhs.ptr(), this->m_size*sizeof(TValue));
        }

      /** Copy data from memory specified by a pointer.
        * This array must be pre-constructed with the desired dimensions.
        * The entire array size will be copied.
        */
      inline void copy(const TValue * rhs) const
        {
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        memcpy(this->m_base, rhs, this->m_size*sizeof(TValue));
        }

      /** Copy data from an existing array.
        * This array must either have existing dimensions the same as the
        * source array, or it can be unconstructed, in which case it will
        * be constructed with the same dimensions.
        */
      template <typename TValue2, typename TIndex2>
      inline void copy_and_convert(const array_base<N,TValue2,TIndex2>& rhs) const
        {
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        else if (this->m_dims != rhs.dims())
          { throw_n88_exception("cannot copy different sized arrays."); }
        for (size_t i=0; i<this->m_size; ++i)
            { this->m_base[i] = rhs[i]; }
        }

      /** Copy data from memory specified by a pointer.
        * This array must be pre-constructed with the desired dimensions.
        * The entire array size will be copied.
        */
      template <typename TValue2>
      inline void copy_and_convert(const TValue2 * rhs) const
        {
        if (!this->m_base)
          { throw_n88_exception("array is not constructed."); }
        for (size_t i=0; i<this->m_size; ++i)
            { this->m_base[i] = rhs[i]; }
        }

    }; // class array_base

  // ---------------------------------------------------------------------

 /**
  * An efficient class for storing multi-dimensional data.
  *
  * This class is intended as a more convenient, but equally fast, alternative
  * to traditional C arrays:
  * @code
  *   float x = new float[100];
  * @endcode
  *
  * arrays can be dynamically allocated on the heap. Allocation can be
  * postponed. arrays cannot be resized.
  *
  * arrays can be indexed with multi-dimensional indices.
  *
  * arrays are indexed with the fastest-changing index last, according to memory
  * layout.  Thus if you have a 3D array where x is the fastest-changing
  * index according to memory layout, it would be indexed as A[z,y,x].
  * This implies the following equivalence: A[i,j] = A[i][j].
  * This indexing is consistent with python and numpy.
  *
  * arrays can also be references to existing memory or arrays.
  *
  * If you want an array referring to constant data, the appropriate class
  * to use is const_array. "const array" is not appropriate, as such an
  * object has constant members, but not a constant pointer to the data.
  * "const array<1,float> x" is equivalent to a simple C style array declared as
  * "float * const x" while "const_array<1,float> x" is equivalent to a simple
  * C style array declared as "float const * x".
  *
  * When compiled with RANGE_CHECKING defined, bounds and other checking is
  * performed.
  * When compiled without RANGE_CHECKING defined, no checking is performed, and there
  * is no additional overhead in indexing or dereferencing.
  *
  * Note that there is no particular advantage to using a type smaller than
  * size_t for TIndex for array itself (the storage difference is just a few
  * bytes in total, and the indexing speed is the same).  However, other
  * classes may want to store large numbers of indices and hence prefer to
  * use a smaller type for the indices.  It is then convenient to be able to
  * define an array class that can be indexed with this smaller type so
  * no conversions need be performed.
  *
  * Note that it is permitted to create an array with a type TIndex and
  * dimensions whose product (i.e. the array size) exceeds the maximum
  * representation of TIndex; this is possible because products of
  * the dimensions are always calculated with type size_t.
  *
  * Note that on allocation array memory is zeroed.  There is a performance-related
  * reason for this: In OS's that use lazy allocation, touching every allocated
  * address at allocation time may ensure that memory is contiguous in real
  * address space, which does in fact sometimes result subsequently in better
  * performance for accessing the array.
  */
  template <int N, typename TValue, typename TIndex=size_t>
  class array : public array_base<N,TValue, TIndex>
    {
    public:

      /** Empty constructor.
        * You must subsequently call construct or construct_reference explicitly.
        */
      array() : array_base<N,TValue,TIndex>() {}

      /** Constructor to allocate space.
        *
        * @param dims  The dimensions of the array to allocate.
        */
      explicit array(Tuple<N,TIndex> dims)  : array_base<N,TValue,TIndex>(dims) {}

      /** Constructor to create reference to existing data defined by a pointer.
        * The referenced data will never be freed by this object.  You
        * must ensure that the referenced memory remains valid, and perform
        * appropriate clean up outside this object.
        *
        * @param data  A pointer to the data.
        * @param dims  The dimensions of the array.
        */
      explicit array(TValue* data, Tuple<N,TIndex> dims) : array_base<N,TValue,TIndex>(data, dims) {}

      /** Constructor to create reference to existing data in array object.
        * No automatic reference counting is performed.  You
        * must ensure that the referenced array remains in existence for as
        * long as you want to use this reference.
        *
        * @param source  An existing array object.
        */
      explicit array(const array_base<N,TValue,TIndex>& source) : array_base<N,TValue,TIndex>(source) {}

      ~array() { this->destruct(); }

    };

  // ---------------------------------------------------------------------

  template <typename TValue, typename TIndex>
  class array<1,TValue,TIndex> : public array_base<1,TValue,TIndex>
    {
    public:

      array() : array_base<1,TValue,TIndex>() {}

      explicit array(Tuple<1,TIndex> dims)  : array_base<1,TValue,TIndex>(dims) {}
      explicit array(TIndex size)  : array_base<1,TValue,TIndex>(Tuple<1,TIndex>(size)) {}

      explicit array(TValue* data, Tuple<1,TIndex> dims) : array_base<1,TValue,TIndex>(data, dims) {}
      explicit array(TValue* data, TIndex size) : array_base<1,TValue,TIndex>(data, Tuple<1,TIndex>(size)) {}

      explicit array(const array_base<1,TValue,TIndex>& source) : array_base<1,TValue,TIndex>(source) {}

      ~array() { this->destruct(); }

      inline void construct(TIndex dim)
        { array_base<1,TValue,TIndex>::construct(Tuple<1,TIndex>(dim)); }

      inline void construct_if_required(TIndex dim)
        { array_base<1,TValue,TIndex>::construct_if_required(Tuple<1,TIndex>(dim)); }

      inline void construct_or_zero(TIndex dim)
        { array_base<1,TValue,TIndex>::construct_or_zero(Tuple<1,TIndex>(dim)); }

      inline void construct_reference(TValue* data, TIndex dim)
        { array_base<1,TValue,TIndex>::construct_reference(data, Tuple<1,TIndex>(dim)); }

      inline size_t flat_index(Tuple<1,TIndex> indices) const
        {
        // Implied static cast from TIndex to size_t
        return indices[0];
        }

      inline TValue& operator()(TIndex i) const
        { return array_base<1,TValue,TIndex>::operator[](static_cast<size_t>(i)); }

      inline TValue& operator()(Tuple<1,TIndex> indices) const
        { return array_base<1,TValue,TIndex>::operator[](this->flat_index(indices)); }

    };

  // ---------------------------------------------------------------------

  template <typename TValue, typename TIndex>
  class array<2,TValue,TIndex> : public array_base<2,TValue,TIndex>
    {
    public:

      array() : array_base<2,TValue,TIndex>() {}

      explicit array(Tuple<2,TIndex> dims)  : array_base<2,TValue,TIndex>(dims) {}
      explicit array(TIndex dim0, TIndex dim1)
          : array_base<2,TValue,TIndex>(Tuple<2,TIndex>(dim0, dim1)) {}

      explicit array(TValue* data, Tuple<2,TIndex> dims) : array_base<2,TValue,TIndex>(data, dims) {}
      explicit array(TValue* data, TIndex dim0, TIndex dim1)
          : array_base<2,TValue,TIndex>(data, Tuple<2,TIndex>(dim0, dim1)) {}

      explicit array(const array_base<2,TValue,TIndex>& source) : array_base<2,TValue,TIndex>(source) {}

      ~array() { this->destruct(); }

      inline void construct(TIndex dim0, TIndex dim1)
        { array_base<2,TValue,TIndex>::construct(Tuple<2,TIndex>(dim0,dim1)); }

      inline void construct_if_required(TIndex dim0, TIndex dim1)
        { array_base<2,TValue,TIndex>::construct_if_required(Tuple<2,TIndex>(dim0,dim1)); }

      inline void construct_or_zero(TIndex dim0, TIndex dim1)
        { array_base<2,TValue,TIndex>::construct_or_zero(Tuple<2,TIndex>(dim0,dim1)); }

      inline void construct_reference(TValue* data, TIndex dim0, TIndex dim1)
        { array_base<2,TValue,TIndex>::construct_reference(data, Tuple<2,TIndex>(dim0,dim1)); }

      inline size_t flat_index(TIndex i, TIndex j) const
        {
        return static_cast<size_t>(i)*static_cast<size_t>(this->m_dims[1])
            + static_cast<size_t>(j);
        }

      inline size_t flat_index(Tuple<2,TIndex> indices) const
        {
        return static_cast<size_t>(indices[0])*static_cast<size_t>(this->m_dims[1])
            + static_cast<size_t>(indices[1]);
        }

      inline TValue& operator()(TIndex i, TIndex j) const
        { return array_base<2,TValue,TIndex>::operator[](this->flat_index(i,j)); }

      inline TValue& operator()(Tuple<2,TIndex> indices) const
        { return array_base<2,TValue,TIndex>::operator[](this->flat_index(indices)); }

    };

  // ---------------------------------------------------------------------

  template <typename TValue, typename TIndex>
  class array<3,TValue,TIndex> : public array_base<3,TValue,TIndex>
    {
    public:

      array() : array_base<3,TValue,TIndex>() {}

      explicit array(Tuple<3,TIndex> dims)  : array_base<3,TValue,TIndex>(dims) {}
      explicit array(TIndex dim0, TIndex dim1, TIndex dim2)
          : array_base<3,TValue,TIndex>(Tuple<3,TIndex>(dim0, dim1, dim2)) {}

      explicit array(TValue* data, Tuple<3,TIndex> dims) : array_base<3,TValue,TIndex>(data, dims) {}
      explicit array(TValue* data, TIndex dim0, TIndex dim1, TIndex dim2)
          : array_base<3,TValue,TIndex>(data, Tuple<3,TIndex>(dim0, dim1, dim2)) {}

      explicit array(const array_base<3,TValue,TIndex>& source) : array_base<3,TValue,TIndex>(source) {}

      ~array() { this->destruct(); }

      inline void construct(TIndex dim0, TIndex dim1, TIndex dim2)
        { array_base<3,TValue,TIndex>::construct(Tuple<3,TIndex>(dim0,dim1,dim2)); }

      inline void construct_if_required(TIndex dim0, TIndex dim1, TIndex dim2)
        { array_base<3,TValue,TIndex>::construct_if_required(Tuple<3,TIndex>(dim0,dim1,dim2)); }

      inline void construct_or_zero(TIndex dim0, TIndex dim1, TIndex dim2)
        { array_base<3,TValue,TIndex>::construct_or_zero(Tuple<3,TIndex>(dim0,dim1,dim2)); }

      inline void construct_reference(TValue* data, TIndex dim0, TIndex dim1, TIndex dim2)
        { array_base<3,TValue,TIndex>::construct_reference(data, Tuple<3,TIndex>(dim0,dim1,dim2)); }

      inline size_t flat_index(TIndex i, TIndex j, TIndex k) const
        {
        return (static_cast<size_t>(i)*static_cast<size_t>(this->m_dims[1])
                + static_cast<size_t>(j))*static_cast<size_t>(this->m_dims[2])
               + static_cast<size_t>(k);
        }

      inline size_t flat_index(Tuple<3,TIndex> indices) const
        {
        return (static_cast<size_t>(indices[0])*static_cast<size_t>(this->m_dims[1])
                + static_cast<size_t>(indices[1]))*static_cast<size_t>(this->m_dims[2])
               + static_cast<size_t>(indices[2]);
        }

      inline TValue& operator()(TIndex i, TIndex j, TIndex k) const
        { return array_base<3,TValue,TIndex>::operator[](this->flat_index(i,j,k)); }

      inline TValue& operator()(Tuple<3,TIndex> indices) const
        { return array_base<3,TValue,TIndex>::operator[](this->flat_index(indices)); }

    };

  // ---------------------------------------------------------------------

  template <typename TValue, typename TIndex>
  class array<4,TValue,TIndex> : public array_base<4,TValue,TIndex>
    {
    public:

      array() : array_base<4,TValue,TIndex>() {}

      explicit array(Tuple<4,TIndex> dims)  : array_base<4,TValue,TIndex>(dims) {}
      explicit array(TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
          : array_base<4,TValue,TIndex>(Tuple<4,TIndex>(dim0, dim1, dim2, dim3)) {}

      explicit array(TValue* data, Tuple<4,TIndex> dims) : array_base<4,TValue,TIndex>(data, dims) {}
      explicit array(TValue* data, TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
          : array_base<4,TValue,TIndex>(data, Tuple<4,TIndex>(dim0, dim1, dim2, dim3)) {}

      explicit array(const array_base<4,TValue,TIndex>& source) : array_base<4,TValue,TIndex>(source) {}

      ~array() { this->destruct(); }

      inline void construct(TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
        { array_base<4,TValue,TIndex>::construct(Tuple<4,TIndex>(dim0,dim1,dim2,dim3)); }

      inline void construct_if_required(TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
        { array_base<4,TValue,TIndex>::construct_if_required(Tuple<4,TIndex>(dim0,dim1,dim2,dim3)); }

      inline void construct_or_zero(TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
        { array_base<4,TValue,TIndex>::construct_or_zero(Tuple<4,TIndex>(dim0,dim1,dim2,dim3)); }

      inline void construct_reference(TValue* data, TIndex dim0, TIndex dim1, TIndex dim2, TIndex dim3)
        { array_base<4,TValue,TIndex>::construct_reference(data, Tuple<4,TIndex>(dim0,dim1,dim2,dim3)); }

      inline size_t flat_index(TIndex i, TIndex j, TIndex k, TIndex l) const
        {
        return ((static_cast<size_t>(i)*static_cast<size_t>(this->m_dims[1])
               + static_cast<size_t>(j))*static_cast<size_t>(this->m_dims[2])
               + static_cast<size_t>(k))*static_cast<size_t>(this->m_dims[3])
               + static_cast<size_t>(l);
        }

      inline size_t flat_index(Tuple<4,TIndex> indices) const
        {
        return ((static_cast<size_t>(indices[0])*static_cast<size_t>(this->m_dims[1])
               + static_cast<size_t>(indices[1]))*static_cast<size_t>(this->m_dims[2])
               + static_cast<size_t>(indices[2]))*static_cast<size_t>(this->m_dims[3])
               + static_cast<size_t>(indices[3]);
        }

      inline TValue& operator()(TIndex i, TIndex j, TIndex k, TIndex l) const
        { return array_base<4,TValue,TIndex>::operator[](this->flat_index(i,j,k,l)); }

      inline TValue& operator()(Tuple<4,TIndex> indices) const
        { return array_base<4,TValue,TIndex>::operator[](this->flat_index(indices)); }

    };

  } // namespace n88

#endif
