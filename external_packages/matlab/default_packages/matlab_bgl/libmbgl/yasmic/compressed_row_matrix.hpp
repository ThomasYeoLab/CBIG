#ifndef YASMIC_COMPRESSED_ROW_MATRIX
#define YASMIC_COMPRESSED_ROW_MATRIX

#if _MSC_VER >= 1400
    // disable the warning for deprecated c++ commands
    #pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400

#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <vector> 

#include <functional>
#include <yasmic/tuple_utility.hpp>
#include <limits>

#include <yasmic/generic_matrix_operations.hpp>

namespace yasmic
{
	namespace impl
	{
		/** 
		 * The nonzero iterator for the crm matrix.
		 */

		template <class RowIter, class ColIter, class ValIter>
		class compressed_row_nonzero_const_iterator
		: public boost::iterator_facade<
            compressed_row_nonzero_const_iterator<RowIter, ColIter, ValIter>,
			yasmic::simple_nonzero<
				typename std::iterator_traits<RowIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				typename std::iterator_traits<RowIter>::value_type> const,
            boost::forward_traversal_tag, 
            yasmic::simple_nonzero<
				typename std::iterator_traits<RowIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				typename std::iterator_traits<RowIter>::value_type> const >
        {
        public:
            compressed_row_nonzero_const_iterator() {}
            
            compressed_row_nonzero_const_iterator(
                RowIter ri, RowIter rend, ColIter ci, ValIter vi)
            : _ri(ri), _rend(rend), _ci(ci), _vi(vi), _id(0), _row(0)
            {
				_ri++;
            	// skip over any empty rows
            	while ( _ri != _rend && *(_ri) == _id )
                {
                   	// keep incrementing the row 
                   	++_ri; ++_row;
                }
            }

			compressed_row_nonzero_const_iterator(
                RowIter ri, RowIter rend, ColIter ci, ValIter vi, 
				typename std::iterator_traits<RowIter>::value_type id,
				typename std::iterator_traits<RowIter>::value_type row = 0)
			: _ri(ri), _rend(rend), _ci(ci), _vi(vi), _id(id), _row(row)
			{
			}
            
            
        private:
            friend class boost::iterator_core_access;

            inline void increment() 
            {  
            	// just increment everything!
            	++_ci; ++_vi; ++_id;
            	
                if (_id == *_ri)
                {
                	++_ri; ++_row;
                	
                	// while we aren't at the end and the row isn't empty
                	// (if *_ri == _id, then the row is empty because _id 
                	// is the current index in the column/val array.)
                	while ( _ri != _rend && *(_ri) == _id )
                    {
                    	// keep incrementing the row 
                    	++_ri; ++_row;
                    }
                }
            }
            
            bool equal(compressed_row_nonzero_const_iterator const& other) const
            {
				return (_ri == other._ri && _ci == other._ci);
            }
            
            yasmic::simple_nonzero<
				typename std::iterator_traits<RowIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				typename std::iterator_traits<RowIter>::value_type> 
            dereference() const 
            { 
            	//return boost::make_tuple(_row, *_ci, *_vi);
				return make_simple_nonzero(_row, *_ci, *_vi, _id);
            }


            typename std::iterator_traits<RowIter>::value_type _row;
            typename std::iterator_traits<RowIter>::value_type _id;
            RowIter _ri, _rend;
            ColIter _ci;
            ValIter _vi;
        };

		

		
		template <class IndexType, class ColIter, class ValIter>
		class crm_row_nonzero_const_iterator
		: public boost::iterator_facade<
            crm_row_nonzero_const_iterator<IndexType, ColIter, ValIter>,
			yasmic::simple_nonzero<
				typename std::iterator_traits<ColIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				IndexType> const,
            boost::forward_traversal_tag, 
            yasmic::simple_nonzero<
				typename std::iterator_traits<ColIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				IndexType> const >
		{
		public:
			crm_row_nonzero_const_iterator() {}
            
            crm_row_nonzero_const_iterator(
                IndexType r, IndexType nzi, ColIter ci, ValIter vi)
            :  _r(r), _nzi(nzi), _ci(ci), _vi(vi)
            {}


		private:
			IndexType _r;
			IndexType _nzi;
			ColIter _ci;
			ValIter _vi;

            friend class boost::iterator_core_access;

			void increment() 
            {
				++_nzi;
				++_ci;
				++_vi;
			}

			bool equal(crm_row_nonzero_const_iterator const& other) const
			{
				return (_ci == other._ci);
			}

			yasmic::simple_nonzero<
				typename std::iterator_traits<ColIter>::value_type,
				typename std::iterator_traits<ValIter>::value_type,
				IndexType>
            dereference() const 
			{
				return (make_simple_nonzero(_r, *_ci, *_vi, _nzi));
			}
		};

		template <class RowIter, class ColIter, class ValIter>
		struct crm_row_nonzero_iter_help
		{
			typedef crm_row_nonzero_const_iterator<
				typename std::iterator_traits<RowIter>::value_type,
				ColIter, ValIter> row_iter;
		};
	}
	
	template <class RowIter, class ColIter, class ValIter>
	class compressed_row_matrix
	{
    public:
		typedef typename std::iterator_traits<RowIter>::value_type size_type;
		
		typedef size_type index_type;
		
		typedef typename std::iterator_traits<ValIter>::value_type value_type;

		typedef size_type nz_index_type;
		
		
		typedef simple_nonzero<index_type, value_type, nz_index_type> nonzero_descriptor;

		typedef impl::compressed_row_nonzero_const_iterator<RowIter, ColIter, ValIter>
             nonzero_iterator;
        
		typedef boost::counting_iterator<size_type> row_iterator;

		typedef nonzero_descriptor row_nonzero_descriptor;
		/*typedef impl::compressed_row_nonzero_const_iterator<RowIter, ColIter, ValIter>
             row_nonzero_iterator;*/

		typedef typename impl::crm_row_nonzero_iter_help<RowIter, ColIter, ValIter>::row_iter
			row_nonzero_iterator;

		typedef ValIter value_iterator;
		
		typedef void column_iterator;
		typedef void column_nonzero_iterator;

		//typedef unsymmetric_tag symmetry_category;
        typedef void properties;
		
		/**
		 * Simple constructor for all the iterators.  This function
		 * requires one pass through the data in order to determine
		 * the number of rows, columns, and non-zeros.
		 */
		compressed_row_matrix(RowIter rstart, RowIter rend, 
            ColIter cstart, ColIter cend, ValIter vstart, ValIter vend)
        : _rstart(rstart), _rend(rend), _cstart(cstart), _cend(cend),
          _vstart(vstart), _vend(vend), _nnz(0), _nrows(0), _ncols(0)
        {
			RowIter ri = _rstart;
			ColIter ci = _cstart;
			ValIter vi = _vstart;
			
			// elminate the first entry, so we get the correct # of rows
			++ri;
			while (ri != _rend)
			{
				++_nrows;  ++ri;
			}
			
			
			while (ci != _cend)
			{
				_ncols = std::max(*ci,_ncols);
				++ci; ++vi;
				++_nnz;
			}
			
			if (vi != _vend)
			{
				// warning, incorrect data!
			}
        }
          
		compressed_row_matrix(RowIter rstart, RowIter rend, 
            ColIter cstart, ColIter cend, ValIter vstart, ValIter vend,
            size_type nrows, size_type ncols, size_type nnz)
        : _rstart(rstart), _rend(rend), _cstart(cstart), _cend(cend),
          _vstart(vstart), _vend(vend), _nrows(nrows), _ncols(ncols), _nnz(nnz)
        {}
        
        nonzero_iterator begin_nonzeros() const
        {
        	return (nonzero_iterator(_rstart, _rend, _cstart, _vstart));
        }
        
        nonzero_iterator end_nonzeros() const
        {
        	return (nonzero_iterator(_rend-1, _rend, _cend, _vstart));
        }
        
        inline std::pair<size_type, size_type> dimensions() const
        {
        	return (std::make_pair(_nrows, _ncols));
        }
        
        inline size_type nnz() const
        {
        	return (_nnz);
        }
        	
        row_iterator begin_rows() const
        {
        	return (row_iterator(0));
        }
        
        row_iterator end_rows() const
        {
        	return (row_iterator(_nrows));
        }

		row_nonzero_iterator begin_row(index_type r) const
		{
			index_type nzi = *(_rstart + (r));
			return (row_nonzero_iterator(r, nzi, _cstart + nzi, _vstart + nzi));
		}
										
							
		row_nonzero_iterator end_row(index_type r) const
        {
			index_type nzi = *(_rstart + (r+1));
			return (row_nonzero_iterator(r, nzi, _cstart + nzi, _vstart + nzi));	
		}
        
		value_iterator begin_values() const
		{
			return (_vstart);
		}

		value_iterator end_values() const
		{
			return (_vend);
		}

        RowIter _rstart,_rend;
        ColIter _cstart,_cend;
        ValIter _vstart,_vend;
        
        size_type _nrows, _ncols, _nnz;
	};



	template <class RowIter, class ColIter, class ValIter, class BFunc >
	void pack_storage(compressed_row_matrix<RowIter, ColIter, ValIter>& m, BFunc b = BFunc())
	{
		// remove any duplicate entries.  This doesn't actually free any 
		// memory, but reduces the non-zero count to the appropriate level
		// for the matrix.

		typedef compressed_row_matrix<RowIter, ColIter, ValIter> Matrix;
		typedef smatrix_traits<Matrix> traits;

		std::vector<typename traits::index_type> wa(ncols(m));

		typename traits::nz_index_type cur_elem = 0;
		
		typename traits::index_type unused = std::numeric_limits<typename traits::index_type>::max();

		std::fill(wa.begin(), wa.end(), unused);

		RowIter ri, riend;
		ri = m._rstart;
		riend = m._rend-1;

		ColIter ci = m._cstart;
		ValIter vi = m._vstart;

		for(; ri != riend; ++ri)
		{
			typename traits::nz_index_type pos_start = *ri;
			typename traits::nz_index_type pos_end = *(ri+1);

			*ri = cur_elem;

			// 
			// look at the elements in order
			//
	        
			for (typename traits::nz_index_type k=pos_start; k<pos_end; k++)
			{
				// the column
				typename traits::index_type c = ci[k]; 
	            
				if (wa[c] == unused)
				{
					// this is a new element
					ci[cur_elem] = c;
					vi[cur_elem] = vi[k];
	                
					wa[c] = cur_elem;
	                
					// now we've used this element...
					cur_elem++;
				}
				else
				{
					// this is a duplicate element
					//vi[wa[c]] += vi[k];
					vi[wa[c]] = b(vi[wa[c]], vi[k]);
				}
			}
	        
			// reset just the entries we used (which are stored in the ci
			// array)
			for (typename traits::nz_index_type k=*ri; k < cur_elem; k++)
			{
				wa[ci[k]] = unused;
			}

		}

		*ri = cur_elem;
	}

    /* ========================================================
     *  Routines to pack storage
     * ===================================================== */


	template <class RowIter, class ColIter, class ValIter >
	void pack_storage(compressed_row_matrix<RowIter, ColIter, ValIter>& m)
	{
		pack_storage(m, std::plus<typename std::iterator_traits<ValIter>::value_type > ());
	}

    /* ========================================================
     *  Routines to sort storage
     * ===================================================== */

	/**
	 * The version with a custom sorting operation.  
	 */
	template <class RowIter, class ColIter, class ValIter, class StrictWeakOrdering>
	void sort_storage(compressed_row_matrix<RowIter, ColIter, ValIter>& m,
		StrictWeakOrdering comp)
	{
        //BOOST_STATIC_ASSERT(0 == 1);
	}

	namespace impl
	{
		template <class ColIter, class ValIter>
		struct crm_col_val_iter_helper_type
		{
			typedef boost::tuple<
					typename std::iterator_traits<ColIter>::value_type, 
					typename std::iterator_traits<ValIter>::value_type >
				value_type;

			typedef boost::tuple<
					typename std::iterator_traits<ColIter>::value_type&, 
					typename std::iterator_traits<ValIter>::value_type& >
				ref_type;
		};

	

		template <class ColIter, class ValIter>
		class crm_col_val_iter
			: public boost::iterator_facade<
				crm_col_val_iter<ColIter, ValIter>,
				typename crm_col_val_iter_helper_type<ColIter, ValIter>::value_type,
				std::random_access_iterator_tag,
				typename crm_col_val_iter_helper_type<ColIter, ValIter>::ref_type,
				typename std::iterator_traits<ColIter>::difference_type>
		{
		public:
			crm_col_val_iter()
			{}

			crm_col_val_iter(ColIter ci, ValIter vi)
				: _ci(ci), _vi(vi)
			{
			}

			ColIter _ci;
			ValIter _vi;


		private:
			friend class boost::iterator_core_access;
            typedef typename std::iterator_traits<ColIter>::difference_type my_difference_type;


			void increment()
			{
				++_ci; ++_vi;
			}

			void decrement()
			{
				--_ci; --_vi;
			}

			bool equal(crm_col_val_iter const& other) const
			{
				return (_ci == other._ci);
			}

			typename crm_col_val_iter_helper_type<ColIter, ValIter>::ref_type dereference() const
			{
				//_tmp = boost::make_tuple(*_ci, *_vi);
				//return (_tmp);
				return (typename crm_col_val_iter_helper_type<ColIter, ValIter>::ref_type(*_ci, *_vi));
			}

			void advance(my_difference_type n)
			{
				_ci += n;
				_vi += n;
			}

			my_difference_type distance_to(crm_col_val_iter const& other) const
			{
				return ( other._ci - _ci);
			}


		};

		template <class ColIter, class ValIter>
		struct crm_col_val_iter_tuple_compare
			: public std::binary_function<
					typename crm_col_val_iter_helper_type<ColIter, ValIter>::value_type,
					typename crm_col_val_iter_helper_type<ColIter, ValIter>::value_type,
					bool>
		{
			typedef typename crm_col_val_iter_helper_type<ColIter, ValIter>::value_type T;
			bool operator()(const  T& t1, const T& t2)
			{
				return (boost::get<0>(t1) < boost::get<0>(t2));
			}
		};

		template <class ColIter, class ValIter>
		crm_col_val_iter<ColIter, ValIter>
		crm_make_col_val_iter(ColIter ci, ValIter vi)
		{
			return crm_col_val_iter<ColIter, ValIter>(ci, vi);
		};
	}


	/**
	 * Sort the storage of a matrix.  This function sorts the entries
	 * in each row by column.
	 */
	template <class RowIter, class ColIter, class ValIter>
	void sort_storage(compressed_row_matrix<RowIter, ColIter, ValIter>& m)
	{
		typedef compressed_row_matrix<RowIter, ColIter, ValIter> Matrix;
		typedef smatrix_traits<Matrix> traits;

		typedef typename traits::index_type itype;
		typedef typename traits::value_type vtype;

		RowIter ri, riend;
		ri = m._rstart;
		riend = m._rend-1;

		ColIter ci = m._cstart;
		ValIter vi = m._vstart;

		for(; ri != riend; ++ri)
		{
			typename traits::nz_index_type pos_start = *ri;
			typename traits::nz_index_type pos_end = *(ri+1);

			std::sort(
				impl::crm_make_col_val_iter(ci + pos_start, vi + pos_start),
				impl::crm_make_col_val_iter(ci + pos_end, vi + pos_end),
				impl::crm_col_val_iter_tuple_compare<ColIter,ValIter>());
		}
	}

    /* ========================================================
     *  Routines to multiply
     * ===================================================== */
	
	/**
	 * Specialize the multiply operator for this matrix type.
	 */
	template <class RowIter, class ColIter, class ValIter, class Iter1, class Iter2>
	void mult(const compressed_row_matrix<RowIter, ColIter, ValIter>& m, Iter1 x, Iter2 y)
	{
		RowIter ri = m._rstart;
		ColIter ci = m._cstart;
		ValIter vi = m._vstart;

		typedef compressed_row_matrix<RowIter, ColIter, ValIter> matrix;
		typedef typename smatrix_traits<matrix>::index_type itype;
		typedef typename smatrix_traits<matrix>::size_type stype;

		stype nr = nrows(m);

		for (itype r=0; r < nr; ++r)
		{
			typedef typename smatrix_traits<matrix>::nz_index_type nzitype;
			typedef typename smatrix_traits<matrix>::value_type vtype;

			vtype ip = 0.0;

			for (nzitype cp = ri[r]; cp < ri[r+1]; ++cp)
			{
				ip += vi[cp]*x[ci[cp]];
			}

			y[r] = ip;
		}
	}

	/**
	 * Specialize the transpose multiply operator for this matrix type.
	 */
	template <class RowIter, class ColIter, class ValIter, class Iter1, class Iter2>
	void trans_mult(const compressed_row_matrix<RowIter, ColIter, ValIter>& m, Iter1 x, Iter2 y)
	{
		RowIter ri = m._rstart;
		ColIter ci = m._cstart;
		ValIter vi = m._vstart;

		typedef compressed_row_matrix<RowIter, ColIter, ValIter> matrix;
		typedef typename smatrix_traits<matrix>::index_type itype;
		typedef typename smatrix_traits<matrix>::size_type stype;

		stype nc = ncols(m);

		// first zero the vector
		for (itype c=0; c < nc; ++c)
		{
			y[c] = 0.0;
		}

		for (itype c=0; c < nc; ++c)
		{
			typedef typename smatrix_traits<matrix>::nz_index_type nzitype;
			typedef typename smatrix_traits<matrix>::value_type vtype;

			vtype cv = x[c];

			for (nzitype cp = ri[c]; cp < ri[c+1]; ++cp)
			{
				y[ci[cp]] += vi[cp]*cv;
			}
		}
	}

    /* ========================================================
     *  Routines to find a value
     * ===================================================== */

    /**
	 * Specialize the value operator for this matrix type.
	 */
	template <class RowIter, class ColIter, class ValIter>
	typename smatrix_traits< compressed_row_matrix<RowIter, ColIter, ValIter> >::value_type value(
        typename smatrix_traits< compressed_row_matrix<RowIter, ColIter, ValIter> >::index_type row,
        typename smatrix_traits< compressed_row_matrix<RowIter, ColIter, ValIter> >::index_type col,
        const compressed_row_matrix<RowIter, ColIter, ValIter>& m)
	{
        typedef compressed_row_matrix<RowIter, ColIter, ValIter> matrix;

		RowIter ri = m._rstart;
		ColIter ci = m._cstart;
		ValIter vi = m._vstart;

        typedef typename smatrix_traits<matrix>::nz_index_type nzi_type;
        
        for (nzi_type cp = ri[row]; cp < ri[row+1]; ++cp)
        {
            if (col == ci[cp])
            {
                return (vi[cp]);
            }
        }

        return (0);
	}

    /**
     * csr_matrix is a more user manageable compressed sparse row matrix type.
     * It is based on the same ideas as the compressed_row_matrix, but designed 
     * to be less general and more specific.  
     *
     * csr_matrix uses compressed_row_matrix for some of its operations.
     *
     * The idea is that an application will use the csr_matrix structure to 
     * manage a sparse matrix and design algorithms that are NOT more 
     * generally applicable. 
     */
    /*template <class IndexType, class ValueType, class SizeType=IndexType>
    struct csr_matrix
    {
        SizeType nrows;
        SizeType ncols;
        SizeType nnz;

        IndexType* ai;

        IndexType* aj;
        ValueType* a;
    };

    template <class IndexType, class ValueType, class SizeType>
    smatrix_traits<csr_matrix<IndexType, ValueType, SizeType> >
    {
        typedef compressed_row_matrix<IndexType*, IndexType*, ValueType*>
            compressed_row_matrix_type;

        typedef IndexType index_type;
        typedef ValueType value_type;
        typedef typename compressed_row_matrix_type::nonzero_iterator nonzero_iterator;
        typedef typename compressed_row_matrix_type::nonzero_descriptor nonzero_descriptor;

	    typedef typename compressed_row_matrix_type::row_iterator row_iterator;
		
	    typedef typename compressed_row_matrix_type::row_nonzero_iterator row_nonzero_iterator;
	    typedef typename compressed_row_matrix_type::row_nonzero_descriptor row_nonzero_descriptor;

	    typedef typename SizeType nz_index_type;
		
	    typedef void column_iterator;
	    //typedef typename Mat::column_nonzero_iterator column_nonzero_iterator;

	    typedef unsymmetric_tag symmetry_category;
	    };
    }

    template <class IndexType, class ValueType, class SizeType>
    SizeType nrows(const csr_matrix<IndexType, ValueType, SizeType>& m)
    { return m.nrows; }

    template <class IndexType, class ValueType, class SizeType>
    SizeType ncols(const csr_matrix<IndexType, ValueType, SizeType>& m)
    { return m.ncols; }

    template <class IndexType, class ValueType, class SizeType>
    SizeType nnz(const csr_matrix<IndexType, ValueType, SizeType>& m)
    { return m.nnz; }*/



}

#if _MSC_VER >= 1400
	// disable the warning for deprecated c++ commands
	#pragma warning( pop )
#endif // _MSC_VER >= 1400


#endif // YASMIC_COMPRESSED_ROW_MATRIX


