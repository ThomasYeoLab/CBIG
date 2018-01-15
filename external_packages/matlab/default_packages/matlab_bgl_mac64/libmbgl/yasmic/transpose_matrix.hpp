#ifndef YASMIC_TRANSPOSE_MATRIX
#define YASMIC_TRANSPOSE_MATRIX

/**
 * @file transpose_matrix.hpp
 * @author David Gleich
 * 
 * A transpose filter for a matrix. 
 *
 * Regardless of the input, the output is just a nonzero-matrix.
 *
 * We just wrap the original dataset and forward calls.  
 * (This does not make a copy of the original matrix!)
 */

namespace yasmic
{

template <class Matrix>
class transpose_matrix
{
public:
	Matrix& _m;

	transpose_matrix(Matrix& m)
		: _m(m)
	{}

	typedef typename smatrix_traits<Matrix>::index_type index_type;
	typedef typename smatrix_traits<Matrix>::value_type value_type;

	typedef typename smatrix_traits<Matrix>::nonzero_descriptor nonzero_descriptor;
	typedef typename smatrix_traits<Matrix>::nonzero_iterator nonzero_iterator;

	typedef void row_iterator;
	typedef void row_nonzero_descriptor;
	typedef void row_nonzero_iterator;

	typedef void column_iterator;

	typedef typename smatrix_traits<Matrix>::size_type size_type;
	typedef typename smatrix_traits<Matrix>::nz_index_type nz_index_type;

    typedef void properties;
};

namespace impl 
{
	template <class Matrix>
	struct transpose_matrix_help
	{
		typedef smatrix_traits<transpose_matrix<Matrix> > traits;

		typedef std::pair<typename traits::size_type, typename traits::size_type> dims_ret_type;
		typedef std::pair<typename traits::nonzero_iterator, typename traits::nonzero_iterator> nzs_ret_type;
		typedef typename traits::size_type nnz_ret_type;
	};
}


template <class Matrix>
inline typename impl::transpose_matrix_help<Matrix>::nnz_ret_type 
nnz(transpose_matrix<Matrix>& tm)
{
	return nnz(tm._m);
}

template <class Matrix>
inline typename impl::transpose_matrix_help<Matrix>::dims_ret_type 
dimensions(transpose_matrix<Matrix>& tm)
{
	typename impl::transpose_matrix_help<Matrix>::dims_ret_type d = dimensions(tm._m);
	return (std::make_pair(d.second, d.first));
	//return dimensions(tm._m);
}

template <class Matrix>
inline typename impl::transpose_matrix_help<Matrix>::nzs_ret_type 
nonzeros(transpose_matrix<Matrix>& tm)
{
	return nonzeros(tm._m);
}

template <class Matrix>
inline typename impl::transpose_matrix_help<Matrix>::traits::index_type
row(typename impl::transpose_matrix_help<Matrix>::traits::nonzero_descriptor nz, transpose_matrix<Matrix>& tm)
{
	return column(nz, tm._m);
}

template <class Matrix>
inline typename impl::transpose_matrix_help<Matrix>::traits::index_type
column(typename impl::transpose_matrix_help<Matrix>::traits::nonzero_descriptor nz, transpose_matrix<Matrix>& tm)
{
	return row(nz, tm._m);
}

} // namespace yasmic

#endif // YASMIC_TRANSPOSE_MATRIX


