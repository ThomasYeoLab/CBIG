#ifndef YASMIC_GENERIC_MATRIX_OPERATIONS
#define YASMIC_GENERIC_MATRIX_OPERATIONS

#include <boost/tuple/tuple.hpp>
#include <yasmic/smatrix_traits.hpp>
#include <utility>

    
namespace yasmic
{   
	/**
	 * Simple nonzero class for most of the work.
	 */

	template <class Index, class Value, class NonzeroIndex>
	struct simple_nonzero
	{
		simple_nonzero()
			: _row(0), _column(0), _val(0), _nzi(0)
		{}

		simple_nonzero(Index row, Index col, Value val, NonzeroIndex nzi)
			: _row(row), _column(col), _val(val), _nzi(nzi)
		{}

		Index _row;
		Index _column;
		Value _val;
		NonzeroIndex _nzi;
		
		bool operator!= (const simple_nonzero& other)
		{
			return (_row != other._row || _column || other._column && _val || other._val && _nzi || other._nzi);
		}
		
		bool operator== (const simple_nonzero& other)
		{
			return (_row == other._row && _column == other._column && _val == other._val && _nzi == other._nzi);
		}

		typedef Index index_type;
		typedef Value value_type;
		typedef NonzeroIndex nz_index_type;
	};

	template <class Index, class Value, class NonzeroIndex>
	simple_nonzero<Index, Value, NonzeroIndex>
	make_simple_nonzero(Index r, Index c, Value v, NonzeroIndex nzi)
	{
		return simple_nonzero<Index, Value, NonzeroIndex>(r, c, v, nzi);
	}


    /**
     * Determine the number of nonzeros in a matrix.
     *
     * This function fowards the call to the matrix class itself,
     * implement a partial specialization if you need other behavior.
     *
     * @param m a reference to the matrix.  
     *
     * @return the number of nonzeros in a matrix.
     */
	template <class Matrix>
	inline typename smatrix_traits<Matrix>::size_type
    nnz(Matrix& m)
    {
    	return (m.nnz());
    }

    /**
     * Return the dimensions of the matrix (nrows, ncols).
     *
     * This function fowards the call to the matrix class itself,
     * implement a partial specialization if you need other behavior.
     *
     * @param m a reference to the matrix.  
     *
     * @return (nrows, ncols).
     */
    template <class Matrix>		
    inline std::pair<typename smatrix_traits<Matrix>::size_type,
                      typename smatrix_traits<Matrix>::size_type>
    dimensions(Matrix& m)
    {
    	return (m.dimensions());
    }
    
    template <class Matrix>
    inline typename smatrix_traits<Matrix>::size_type 
    ncols(Matrix& m)
	{
		typename smatrix_traits<Matrix>::size_type nrows,ncols;
		
		boost::tie(nrows, ncols) = dimensions(m);
		
		return (ncols);
	}
	
	template <class Matrix>
    inline typename smatrix_traits<Matrix>::size_type 
	nrows(Matrix& m)
	{ 
		typename smatrix_traits<Matrix>::size_type nrows,ncols;
		
		boost::tie(nrows, ncols) = dimensions(m);
		
		return (nrows);
	}
	
    template <class Matrix>
    inline typename smatrix_traits<Matrix>::index_type 
    row(boost::tuples::tuple<
            typename smatrix_traits<Matrix>::index_type, 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f) 
    {
    	return (boost::tuples::get<0>(n));
    }

	template <class Matrix>
	inline typename smatrix_traits<Matrix>::index_type
	row(simple_nonzero<
			typename smatrix_traits<Matrix>::index_type, 
			typename smatrix_traits<Matrix>::value_type,
			typename smatrix_traits<Matrix>::nz_index_type> n, const Matrix& f)
	{
		return n._row;
	}
    
    template <class Matrix>
    inline typename smatrix_traits<Matrix>::index_type 
    column(boost::tuples::tuple<
            typename smatrix_traits<Matrix>::index_type, 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f)
    {
    	return (boost::tuples::get<1>(n));
    }

	template <class Matrix>
	inline typename smatrix_traits<Matrix>::index_type
	column(simple_nonzero<
			typename smatrix_traits<Matrix>::index_type, 
			typename smatrix_traits<Matrix>::value_type,
			typename smatrix_traits<Matrix>::nz_index_type> n, const Matrix& f)
	{
		return n._column;
	}

	/*template <class Matrix>
    inline typename smatrix_traits<Matrix>::index_type 
    column(boost::tuples::tuple< 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f)
    {
    	return (boost::tuples::get<0>(n));
    }*/

	template <class Matrix>
    inline typename smatrix_traits<Matrix>::index_type 
		column(std::pair< 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f)
    {
    	return (n.first);
    }
    
    template <class Matrix>
    inline typename smatrix_traits<Matrix>::value_type 
    value(boost::tuples::tuple<
            typename smatrix_traits<Matrix>::index_type, 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f) 
    {
    	return (boost::tuples::get<2>(n));
    }

	template <class Matrix>
	inline typename smatrix_traits<Matrix>::value_type
	value(simple_nonzero<
			typename smatrix_traits<Matrix>::index_type, 
			typename smatrix_traits<Matrix>::value_type,
			typename smatrix_traits<Matrix>::nz_index_type> n, const Matrix& f)
	{
		return n._val;
	}

	/**
	 * Overloaded function to return the value for a row/column non-zero as well.
	 */
	/*template <class Matrix>
    inline typename smatrix_traits<Matrix>::value_type 
    row_value(boost::tuples::tuple< 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f) 
    {
    	return (boost::tuples::get<1>(n));
    }*/

	/*template <class Matrix>
	inline typename smatrix_traits<Matrix>::value_type
	row_value(boost::tuples::tuple<
            typename smatrix_traits<Matrix>::index_type, 
            typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f)
	{
		return (boost::tuples::get<1>(n));
	}*/

	template <class Matrix>
    inline typename smatrix_traits<Matrix>::value_type 
	value(std::pair<
			typename smatrix_traits<Matrix>::index_type,
            typename smatrix_traits<Matrix>::value_type> n, const Matrix& f) 
    {
    	return (n.second);
    }
    
    template <class Matrix>
    inline std::pair<typename smatrix_traits<Matrix>::nonzero_iterator,
                      typename smatrix_traits<Matrix>::nonzero_iterator>
    nonzeros(Matrix& m)
    {
		return (std::make_pair(m.begin_nonzeros(), m.end_nonzeros()));
    }

    template <class Matrix>
    inline std::pair<typename smatrix_traits<Matrix>::row_iterator,
                     typename smatrix_traits<Matrix>::row_iterator>
    rows(Matrix& m)
    {
        return (std::make_pair(m.begin_rows(), m.end_rows()));
    }
    
    template <class Matrix>
    inline std::pair<typename smatrix_traits<Matrix>::row_nonzero_iterator,
                      typename smatrix_traits<Matrix>::row_nonzero_iterator>
    row_nonzeros(typename smatrix_traits<Matrix>::index_type r, const Matrix& m)
    {
		typename smatrix_traits<Matrix>::row_nonzero_iterator br = m.begin_row(r);
		typename smatrix_traits<Matrix>::row_nonzero_iterator be = m.end_row(r);
		return std::make_pair(br, be);
    	//return (make_pair(m.begin_row(r), m.end_row(r)));
    }

    template <class Matrix>
    typename smatrix_traits<Matrix>::value_type value(
        typename smatrix_traits<Matrix>::index_type r,
        typename smatrix_traits<Matrix>::index_type c,
        const Matrix& m)
    {
        using namespace yasmic;

        //BOOST_STATIC_ASSERT(0 == 1);

        typedef typename smatrix_traits<Matrix>::index_type itype;
        typedef typename smatrix_traits<Matrix>::value_type vtype;

        typename smatrix_traits<Matrix>::nonzero_iterator nzi, nziend;
		for (boost::tie(nzi,nziend) = nonzeros(m);
			 nzi != nziend; ++nzi)
		{
			if (row(*nzi,m)==r && column(*nzi,m)==c)
            {
                return (value(*nzi,m));
            }
		}

        return (0);
    }

	template <class Matrix, class Iter1, class Iter2>
	void mult(Matrix& m, Iter1 x, Iter2 y)
	{
		using namespace yasmic;

		
		typedef typename smatrix_traits<Matrix>::index_type itype;
		typedef typename smatrix_traits<Matrix>::size_type stype;
		typedef typename smatrix_traits<Matrix>::value_type vtype;

		stype nr = nrows(m);

		// first zero the vector
		for (itype r=0; r < nr; ++r)
		{
			y[r] = vtype();
		}

		typename smatrix_traits<Matrix>::nonzero_iterator nzi, nziend;
		for (boost::tie(nzi,nziend) = nonzeros(m);
			 nzi != nziend; ++nzi)
		{
			y[row(*nzi,m)] += value(*nzi,m)*x[column(*nzi,m)];
		}
	}

	template <class Matrix, class Iter1, class Iter2>
	void trans_mult(Matrix& m, Iter1 x, Iter2 y)
	{
		using namespace yasmic;

		
		typedef typename smatrix_traits<Matrix>::index_type itype;
		typedef typename smatrix_traits<Matrix>::size_type stype;
		typedef typename smatrix_traits<Matrix>::value_type vtype;

		stype nc = ncols(m);

		// first zero the vector
		for (itype c=0; c < nc; ++c)
		{
			y[c] = vtype();
		}

		typename smatrix_traits<Matrix>::nonzero_iterator nzi, nziend;
		for (boost::tie(nzi,nziend) = nonzeros(m);
			 nzi != nziend; ++nzi)
		{
			y[column(*nzi,m)] += value(*nzi,m)*x[row(*nzi,m)];
		}
	}

}

#endif // YASMIC_GENERIC_MATRIX_OPERATIONS

