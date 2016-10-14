#ifndef YASMIC_UTIL_FILTERED_MATRIX
#define YASMIC_UTIL_FILTERED_MATRIX

#if _MSC_VER >= 1400
    // disable the warning for deprecated c++ commands
    #pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400

#include <boost/iterator/filter_iterator.hpp>
#include <algorithm>
#include <functional>

#include <boost/type_traits.hpp>

namespace yasmic
{
	namespace impl
	{
		template <class Matrix> 
		struct filtered_matrix_nonzero_helper
		{
			typedef typename smatrix_traits<Matrix>::size_type s;
			typedef typename smatrix_traits<Matrix>::value_type v;

			typedef boost::tuple<s,s,v> nz;
		};

		template <class Matrix, class FilterMap>
		class filtered_matrix_const_iterator
		: public boost::iterator_facade<
            filtered_matrix_const_iterator<Matrix, FilterMap>,
			typename filtered_matrix_nonzero_helper<Matrix>::nz const,
            boost::forward_traversal_tag, 
            typename filtered_matrix_nonzero_helper<Matrix>::nz  const >
        {
        public:
			filtered_matrix_const_iterator() 
				: _m(0), _fm(NULL)
			{}

            filtered_matrix_const_iterator(const Matrix& m, FilterMap fm,
				std::pair<typename smatrix_traits<Matrix>::nonzero_iterator,
						typename smatrix_traits<Matrix>::nonzero_iterator> i) 
				: _m(&m), _fm(&fm), _i(i.first), _iend(i.second)
			{
				skip_nnzs();
			}
            

            
        private:
            friend class boost::iterator_core_access;

            void increment() 
            {  
				++_i;
				skip_nnzs();
            }
            
            bool equal(filtered_matrix_const_iterator const& other) const
            {
				return (_i == other._i);
            }
            
            typename filtered_matrix_nonzero_helper<Matrix>::nz
            dereference() const 
            { 
				typename impl::filtered_matrix_nonzero_helper<Matrix>::s r,c;

				r = (*_fm)[row(*_i,*_m)]-1;
				c = (*_fm)[column(*_i,*_m)]-1;
				
            	return boost::make_tuple(r, c, value(*_i,*_m));
            }

			const Matrix *_m;
			typename boost::add_pointer<FilterMap>::type _fm;
			typename smatrix_traits<Matrix>::nonzero_iterator _i;
			typename smatrix_traits<Matrix>::nonzero_iterator _iend;

			void skip_nnzs()
			{
				// advance to the next non-filtered nnz.
				while (_i != _iend && (!(*_fm)[row(*_i,*_m)] || !(*_fm)[column(*_i,*_m)]))
				{
					++_i;
				}
			}
        };
	}

	template <class Matrix, class RowColFilterMap>
	class row_column_filtered_matrix
	{
	public:
		// handle the matrix traits
		typedef typename smatrix_traits<Matrix>::index_type index_type;
		typedef typename smatrix_traits<Matrix>::value_type value_type;
		typedef typename smatrix_traits<Matrix>::size_type size_type;

		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;
		typedef impl::filtered_matrix_const_iterator<Matrix, RowColFilterMap> nonzero_iterator;

		typedef void row_iterator;
		
		typedef void row_nonzero_iterator;
		typedef void row_nonzero_descriptor;
		typedef void column_iterator;

		typedef typename smatrix_traits<Matrix>::properties properties;

		typedef typename smatrix_traits<Matrix>::nz_index_type nz_index_type;

		row_column_filtered_matrix(const Matrix &m, RowColFilterMap fm)
			: _nzi(m, fm, nonzeros(m)), 
			  _nziend(m, fm, std::make_pair(nonzeros(m).second, nonzeros(m).second)),
			  _nrows(0), _ncols(0), _nnz(0)
		{
			nonzero_iterator i = _nzi;

			_nnz = 0;
			
			while (i != _nziend)
			{
				_nrows = std::max(row(*i,m),_nrows);
				_ncols = std::max(column(*i,m),_ncols);

				++i;
				++_nnz;
			}

			// we underestimate by one value
			++_nrows;
			++_ncols;
		}

 		row_column_filtered_matrix(const Matrix &m, RowColFilterMap fm,
			size_type nrows, size_type ncols)
			: _nzi(m, fm, nonzeros(m)), 
			  _nziend(m, fm, std::make_pair(nonzeros(m).second, nonzeros(m).second)),
			  _nrows(nrows), _ncols(nrows), _nnz(0)
		{
			nonzero_iterator i = _nzi;

			_nnz = 0;
			
			while (i != _nziend)
			{
				++i;
				++_nnz;
			}
		}

		row_column_filtered_matrix(const Matrix &m, RowColFilterMap fm,
			size_type nrows, size_type ncols, size_type nnz)
			: _nzi(m, fm, nonzeros(m)), 
			  _nziend(m, fm, std::make_pair(nonzeros(m).second, nonzeros(m).second)),
			  _nrows(nrows), _ncols(nrows), _nnz(nnz)
		{}

		std::pair<size_type, size_type> dimensions() const
		{
			return std::make_pair(_nrows, _ncols);
		}

		size_type nnz() const
		{
			return _nnz;
		}

		nonzero_iterator begin_nonzeros() const
		{
			return _nzi;
		}

		nonzero_iterator end_nonzeros() const
		{
			return _nziend;
		}

	private:
		size_type _nnz;
		size_type _nrows;
		size_type _ncols;

		nonzero_iterator _nzi;
		nonzero_iterator _nziend;
	};

	/*template <class Matrix, class RowColFilterMap>
	class row_column_filtered_matrix
	{
	private:
		typedef is_filtered_nonzero<nonzero_descriptor, Matrix, RowColFilterMap>
			filter_predicate;

	public:
		// handle the matrix traits
		typedef typename smatrix_traits<Matrix>::index_type index_type;
		typedef typename smatrix_traits<Matrix>::value_type value_type;

		typedef typename smatrix_traits<Matrix>::nonzero_descriptor nonzero_descriptor;

		typedef typename smatrix_traits<Matrix>::nonzero_iterator mat_nonzero_iterator;
		
		

		typedef typename smatrix_traits<Matrix>::row_iterator void;
		
		typedef typename smatrix_traits<Matrix>::row_nonzero_iterator void;
		typedef typename smatrix_traits<Matrix>::row_nonzero_descriptor void;

		typedef typename smatrix_traits<Matrix>::nz_index_type nz_index_type;
		
		typedef typename smatrix_traits<Matrix>::column_iterator void;

		row_column_filtered_matrix(const Matrix& im, RowColFilterMap fm)
			: m(im), rcfm(fm)
		{

		}

		std::pair<size_type, size_type> dimensions()
		{
			return make_pair(nrows, ncols);
		}

		size_type nnz()
		{
			return size_type;
		}

		nonzero_iterator begin_nonzeros()
		{
			filter_predicate p;
			p.m = m; p.fm = rcfm;

			mat_nonzero_iterator nzi, nziend
			boost::tie(nzi, nziend)  = nonzeros(m);
			return nonzero_iterator(p, nzi);
		}

		nonzero_iterator end_nonzeros()
		{
			filter_predicate p;
			p.m = m; p.fm = rcfm;

			mat_nonzero_iterator nzi, nziend
			boost::tie(nzi, nziend)  = nonzeros(m);
			return nonzero_iterator(p, nziend);
		}

	private:
		const Matrix& m;
		RowColFilterMap rcfm;

		size_type nnz;
		size_type nrows;
		size_type ncols;
	};*/

	template <class Iterator>
	void make_filter_map_from_indicator(Iterator start, Iterator end)
	{
		typedef typename std::iterator_traits<Iterator>::value_type t;

		std::partial_sum(
			boost::make_filter_iterator(
				std::bind2nd(std::greater<t>(), 0),
				start, end),
			boost::make_filter_iterator(
				std::bind2nd(std::greater<t>(), 0),
				end, end),
			boost::make_filter_iterator(
				std::bind2nd(std::greater<t>(), 0),
				start, end));


		//boost::filter_iterator<
		//	          std::bind2nd(std::greater<t>(), 0)
	}
};

#endif // YASMIC_UTIL_FILTERED_MATRIX
