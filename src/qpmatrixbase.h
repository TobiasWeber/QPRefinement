/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  qpmatrixbase.h
 * @brief set of rows/cols for symmetric qp matrix.
 */

#ifndef _QPMATRIXBASE_H_
#define _QPMATRIXBASE_H_

#include <assert.h>

#include "spxdefines.h"
#include "basevectors.h"
#include "datakey.h"

namespace soplex
{
/**@brief   set qp matrix row/cols.
 * @ingroup Algebra
 *
 *  Class QPMatrixBase implements an set of rows/cols (SVectorBase%s) that represent the quadratic term of the objective
 *  function of an QP (quadratic programm). As the matrix is symmetric only the lower triangular matrix is saved and upper
 *  elements are saved on the lower symmetric position. Therefore the whole symmetric matrix MUST NOT BE SAVED (no check
 *  here!), but only each lower (or upper) triangular value once!
 *
 *  Unless for memory limitations, any number of SVectorBase%s may be #add%ed to an QPMatrixBase. Single or multiple
 *  SVectorBase%s may be added to an QPMatrixBase, where each method add() comes with two different signatures. One
 *  with and one without a parameter, used for returning the Keys assigned to the new SVectorBase%s by the set. See
 *  DataKey for a more detailed description of the concept of keys. For the concept of renumbering SVectorBase%s within
 *  an QPMatrixBase after removal of some LPRows see DataSet.
 *
*/
template < class R >
class QPMatrixBase : protected SVSetBase<R>
{
	template < class S > friend class QPMatrixBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   //@}

protected:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Helpers */
   //@{

   /// Returns the complete SVSet.
   const SVSetBase<R>* matrixQP() const
   {
      return this;
   }

   //@}

public:

// ------------------------------------------------------------------------------------------------------------------
   /**@name Access / modification */
   //@{

   /// Returns the number of QPMatrix rows/cols in QPMatrixBase.
   int num() const
   {
      return SVSetBase<R>::num();
   }

   /// Returns number of nonzeros offdiagonal entries and all diagonal entries in QPMatrix.
   int nNzos() const
   {

      int n = 0;
      for( int i = 0; i < num(); ++i )
      {
	 // count off diagonal entries two times and diagonal entries even when they are zero
	 if( rowcolVector(i).number(i) >= 0  )
            n += 2*rowcolVector(i).size()-1;
         else
	    n += 2*rowcolVector(i).size()+1;
      }

      return n;
   }

   /// Returns maximum number of QPMatrix rows/cols currently fitting into QPMatrixBase.
   int max() const
   {
      return SVSetBase<R>::max();
   }

   /// Returns a writable rowcolVector of the \p i 'th QPMatrix row/col.
   SVectorBase<R>& rowcolVector_w(int i)
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns rowcolVector of \p i 'th QPMatrix row/col.
   const SVectorBase<R>& rowcolVector(int i) const
   {
      return SVSetBase<R>::operator[](i);
   }

   /// Returns writeable rowcolVector of SVSetBase with DataKey \p k in QPMatrixBase.
   SVectorBase<R>& rowcolVector_w(const DataKey& k)
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns rowcolVector of SVSetBase with DataKey \p k in QPMatrixBase.
   const SVectorBase<R>& rowcolVector(const DataKey& k) const
   {
      return SVSetBase<R>::operator[](k);
   }

   /// Returns DataKey of \p i 'th SVSetBase in QPMatrixBase.
   DataKey key(int i) const
   {
      return SVSetBase<R>::key(i);
   }

   /// Returns number of SVSetBase with DataKey \p k in QPMatrixBase.
   int number(const DataKey& k) const
   {
      return SVSetBase<R>::number(k);
   }

   /// Does DataKey \p k belong to QPMatrixBase ?
   bool has(const DataKey& k) const
   {
      return SVSetBase<R>::has(k);
   }

   /// Updates linearization of quadratic term in objective for a given vector for multiplication, i.e., Q*x
   void addQPMatrixvectorproduct(const VectorBase<R>& multiplication, VectorBase<R>& linearization) const
   {
      for( int c = num()-1; c >=0; c-- )
      {
         const SVectorBase<R>& qpmatrixrowcol = rowcolVector(c);
         if( multiplication[c] != 0 )
	 {
            for( int i = qpmatrixrowcol.size()-1; i>=0; i-- )
            {
	       linearization[qpmatrixrowcol.index(i)] += multiplication[c] * qpmatrixrowcol.value(i);
               if( qpmatrixrowcol.index(i) != c && multiplication[qpmatrixrowcol.index(i)] != 0  )
                  linearization[c] += multiplication[qpmatrixrowcol.index(i)] * qpmatrixrowcol.value(i);
            }
         }
         else
         {
            for( int i = qpmatrixrowcol.size()-1; i>=0; i-- )
            {
               if( qpmatrixrowcol.index(i) != c && multiplication[qpmatrixrowcol.index(i)] != 0  )
                  linearization[c] += multiplication[qpmatrixrowcol.index(i)] * qpmatrixrowcol.value(i);
            }
         }
      }
   }

   /// Returns one entry of linearization of quadratic term in objective for a given vector for multiplication, i.e., Q*x
   R QPMatrixvectorproductentry(const VectorBase<R>& multiplication, const int index)
   {
      assert( index < num());
      R value = 0;
      int n = 0;

      for( int c = index-1; c>=0; c-- )
      {
         const SVectorBase<R>& qpmatrixrowcol = rowcolVector(c);
         if( multiplication[c] != 0 )
	 {
	    n = qpmatrixrowcol.number(index);
	    if( n >= 0 )
            {
	      value += qpmatrixrowcol.value(n) * multiplication[c];
            }
         }
      }

      const SVectorBase<R>& qpmatrixrowcol = rowcolVector(index);
      for( int c = qpmatrixrowcol.size()-1; c >=0; c-- )
      {
	 n = qpmatrixrowcol.index(c);
         if( multiplication[n] != 0 )
         {
            value += qpmatrixrowcol.value(c) * multiplication[n];
	 }
      }

      return value;
   }

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Extension
    *
    *  Extension methods come with two signatures, one of them providing a parameter to return the assigned
    *  DataKey(s). See DataSet for a more detailed description. All extension methods will automatically rearrange or
    *  allocate more memory if required.
   */
   //@{

   /// Adds \p rowcol to QPMatrixBase.
   void add(const SVectorBase<R>& rowcol)
   {
      SVSetBase<R>::add(rowcol);
   }

   /// Adds \p rowcol to QPMatrixBase.
   void add(DataKey& pkey, const SVectorBase<R>& prowcol)
   {
      SVSetBase<R>::add(pkey, prowcol);
   }

   /// Adds all SVSetBase%s of newset to QPMatrixBase.
   void add(const SVSetBase<R>& newset)
   {
      SVSetBase<R>::add(newset);
   }

   /// Adds SVectorBase consisting of rowcol vector \p rowcolVector.
   template < class S >
   void add(const S* rowcolValues, const int* rowcolIndices, int rowcolSize)
   {
      assert(rowcolSize <= 0 || rowcolValues != 0);
      assert(rowcolSize <= 0 || rowcolIndices != 0);

      DataKey k;
      add(k, rowcolValues, rowcolIndices, rowcolSize);
   }

   /// Adds SVectorBase consisting of rowcol vector \p rowcolVector, with DataKey \p key.
   template < class S >
   void add(DataKey& newkey, const S* rowcolValues, const int* rowcolIndices, int rowcolSize)
   {
      assert(rowcolSize <= 0 || rowcolValues != 0);
      assert(rowcolSize <= 0 || rowcolIndices != 0);

      SVSetBase<R>::add(newkey, rowcolValues, rowcolIndices, rowcolSize);
   }

   /// Adds all SVSetBase%s of newset to QPMatrixBase.
   template < class S >
   void add(const SVSetBase<S>& newset)
   {
      SVSetBase<R>::add(newset);
   }

   /// Adds all SVSetBase%s of newset to QPMatrixBase and returns keys
   void add(DataKey keys[], const SVSetBase<R>& set)
   {
      int i = num();

      add(set);

      for( int j = 0; i < num(); ++i, ++j )
         keys[j] = key(i);
   }

   /// Adds all SVSetBase%s of newset to QPMatrixBase and returns keys
   template < class S >
   void add(DataKey keys[], const SVSetBase<S>& set)
   {
      int i = num();

      add(set);

      for( int j = 0; i < num(); ++i, ++j )
         keys[j] = key(i);
   }

   /// Extends rowcol \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSetBase<R>::xtend(rowcolVector_w(n), newmax);
   }

   /// Extends row with DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSetBase<R>::xtend(rowcolVector_w(pkey), pnewmax);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to rowVector with DataKey \p k.
   void add2(const DataKey& k, int n, const int idx[], const R val[])
   {
      add2(number(k), n, idx, val);
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th rowcolVector.
   void add2(int i, int n, const int idx[], const R val[])
   {
      // Check if any entry is above diagonal
      for( int j = 0; j < n; ++j )
      {
         if(idx[j]>=i)
         {
            assert(rowcolVector(i).number(idx[j])==-1);
            SVSetBase<R>::add2(rowcolVector_w(i), 1, idx[j], val[j]);
         }
         else
         {
            assert(rowcolVector(idx[j]).number(i)==-1);
            SVSetBase<R>::add2(rowcolVector_w(idx[j]), 1, i, val[j]);
         }
      }
   }

   /// Adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th rowcolVector.
   template < class S >
   void add2(int i, int n, const int idx[], const S val[])
   {
      // Check if any entry is above diagonal
      for( int j = 0; j < n; ++j )
      {
         if(idx[j]>=i)
         {
            assert(rowcolVector(i).number(idx[j])==-1);
            SVSetBase<R>::add2(rowcolVector_w(i), 1, idx[j], val[j]);
         }
         else
         {
            assert(rowcolVector(idx[j]).number(i)==-1);
            SVSetBase<R>::add2(rowcolVector_w(idx[j]), 1, i, val[j]);
         }
      }
   }

   /// Creates new SVSetBase and returns a reference to its rowcol vector.
   SVectorBase<R>& create(int pnonzeros = 0)
   {
      DataKey k;
      return create(k, pnonzeros);
   }

   /// Creates new SVSetBase and returns a reference to its rowcol vector.
   SVectorBase<R>& create(DataKey& newkey, int nonzeros = 0)
   {
      return *SVSetBase<R>::create(newkey, nonzeros);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Shrinking
    *
    *  See DataSet for a description of the renumbering of the remaining SVSetBase%s in a QPMatrixBase after the call of
    *  a removal method.
    */
   //@{

   /// Removes \p i 'th SVSetBase.
   void remove(int i)
   {
      SVSetBase<R>::remove(i);
   }

   /// Removes SVSetBase with DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }

   /// Removes multiple SVSetBase%s.
   void remove(int perm[])
   {
      SVSetBase<R>::remove(perm);
   }

   /// Removes \p n SVSetBase%s with row numbers given by \p nums.
   void remove(const int nums[], int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// Removes \p n SVSetBase%s with row numbers given by \p nums,
   /// Stores permutation of row indices in \p perm.
   void remove(const int nums[], int n, int* perm)
   {
      SVSetBase<R>::remove(nums, n, perm);
   }

   /// Removes all LPRowBase%s.
   void clear()
   {
      SVSetBase<R>::clear();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Memory Management
    *
    *  For a description of the memory management methods, see the documentation of SVSet, which has been used for
    *  implementating QPMatrixBase.
    */
   //@{

   /// Reallocates memory to be able to store \p newmax SVSetBase%s.
   void reMax(int newmax = 0)
   {
      SVSetBase<R>::reMax(newmax);
   }

   /// Returns number of used nonzero entries.
   int memSize() const
   {
      return SVSetBase<R>::memSize();
   }

   /// Returns length of nonzero memory.
   int memMax() const
   {
      return SVSetBase<R>::memMax();
   }

   /// Reallocates memory to be able to store \p newmax nonzeros.
   void memRemax(int newmax)
   {
      SVSetBase<R>::memRemax(newmax);
   }

   /// Garbage collection in nonzero memory.
   void memPack()
   {
      SVSetBase<R>::memPack();
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Consistency check */
   //@{

   /// Checks consistency.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      const int num_rowcol = num();
      int vec_dim = 0;
      int max_vec_dim = 0;

      for( int j = 0; j < num_rowcol; ++j )
      {
         const SVectorBase<Real>& rowcol = rowcolVector(j);
         vec_dim=rowcol.dim();

         if( vec_dim > j )
            return MSGinconsistent("QPMatrixBase");

         if( max_vec_dim<vec_dim )
            max_vec_dim=vec_dim;
      }

      if( max_vec_dim > num_rowcol )
         return MSGinconsistent("QPMatrixBase");

      return SVSetBase<R>::isConsistent();
#else
      return true;
#endif
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction / Destruction */
   //@{

   /// Default constructor.
   /** The user can specify the initial maximum number of rows \p max and the initial maximum number of nonzero entries
    *  \p memmax. If these parameters are omitted, a default size is used. However, one can add an arbitrary number of
    *  rows to the LPRowSetBase, which may result in automated memory realllocation.
    */
   explicit
   QPMatrixBase<R>(int pmax = -1, int pmemmax = -1)
      : SVSetBase<R>(pmax, pmemmax)
   {
      assert(isConsistent());
   }

   /// Assignment operator.
   QPMatrixBase<R>& operator=(const QPMatrixBase<R>& rs)
   {
      if( this != &rs )
      {
         SVSetBase<R>::operator=(rs);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   QPMatrixBase<R>& operator=(const QPMatrixBase<S>& rs)
   {
      if( this != (const QPMatrixBase<R>*)(&rs) )
      {
         SVSetBase<R>::operator=(rs);

         assert(isConsistent());
      }

      return *this;
   }

   /// Copy constructor.
   QPMatrixBase<R>(const QPMatrixBase<R>& rs)
      : SVSetBase<R>(rs)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   QPMatrixBase<R>(const QPMatrixBase<S>& rs)
      : SVSetBase<R>(rs)
   {
      assert(isConsistent());
   }

   /// Destructor.
   ~QPMatrixBase<R>()
   {}

   //@}
};
} // namespace soplex
#endif // _QPMATRIXBASE_H_
