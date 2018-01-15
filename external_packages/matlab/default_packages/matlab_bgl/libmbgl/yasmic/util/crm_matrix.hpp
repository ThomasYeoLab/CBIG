#ifndef YASMIC_UTIL_CRM_MATRIX
#define YASMIC_UTIL_CRM_MATRIX

#if _MSC_VER >= 1400
	// disable the warning for ifstream::read
    #pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400

#include <iterator>
#include <numeric>
#include <algorithm>
#include <vector>

#include <yasmic/transpose_matrix.hpp>
#include <yasmic/nonzero_union.hpp>
#include <yasmic/util/load_crm_matrix.hpp>
#include <yasmic/matrix_row_col_graph.hpp>

/**
 * Symmetrize a CRM matrix by adding a (j,i,v) for each (i,j,v) pair.  This
 * operation doubles the number of nonzeros in the matrix and returns
 * a matrix with nr = nc.
 *
 * All parameters are input/output.
 *
 * @param rows the crm rows vector
 * @param cols the crm cols vector
 * @param vals the crm vals vector
 * @param nr the number of rows
 * @param nc the number of columns
 * @param nzcount the number of nonzeros
 */
template <class index_type, class value_type>
void symmetrize_crm(std::vector<index_type>& rows, std::vector<index_type>& cols, std::vector<value_type>& vals, 
               index_type& nr, index_type& nc, index_type& nzcount)
{
    using namespace yasmic;

    nzcount = (index_type)(2*cols.size());

    std::vector<index_type> rows_temp(rows);
	std::vector<index_type> cols_temp(cols);
	std::vector<value_type> vals_temp(vals);

	std::fill(rows.begin(), rows.end(), 0);
	cols.resize(nzcount);
	vals.resize(nzcount);

    typedef compressed_row_matrix<
		typename std::vector<index_type>::iterator, 
		typename std::vector<index_type>::iterator,
		typename std::vector<value_type>::iterator  >
        crs_matrix;  

	typedef transpose_matrix<crs_matrix> t_matrix;
	typedef nonzero_union<crs_matrix, t_matrix> nzu_matrix;



	crs_matrix m(rows_temp.begin(), rows_temp.end(), cols_temp.begin(), cols_temp.end(), 
				vals_temp.begin(), vals_temp.end(), nr, nc, nzcount/2);

	t_matrix mt(m);
	nzu_matrix nzu(m, mt);

	nr = nrows(nzu);
	nc = ncols(nzu);

	// load the matrix
	load_matrix_to_crm(nzu, rows.begin(), cols.begin(), vals.begin());
}

template<class Type>
struct max_fo
	: public std::binary_function<Type, Type, Type>
	{	// functor for operator+
	Type operator()(const Type& _Left, const Type& _Right) const
		{	// apply operator+ to operands
			return std::max(_Left, _Right);
		}
	};

/**
 * Pack and sort the storage of a CRM matrix.  
 *
 * All parameters are input/output.
 *
 * @param rows the crm rows vector
 * @param cols the crm cols vector
 * @param vals the crm vals vector
 * @param nr the number of rows
 * @param nc the number of columns
 * @param nzcount the number of nonzeros
 */
template <class index_type, class value_type>
void pack_and_sort_storage_crm(std::vector<index_type>& rows, std::vector<index_type>& cols, std::vector<value_type>& vals, 
               index_type& nr, index_type& nc, index_type& nzcount)
{
    using namespace yasmic;

    typedef compressed_row_matrix<
		typename std::vector<index_type>::iterator,
		typename std::vector<index_type>::iterator,
		typename std::vector<value_type>::iterator >
        crs_matrix;  

    crs_matrix mlarge(rows.begin(), rows.end(), cols.begin(),cols.end(), 
					vals.begin(), vals.end(), nr, nc, nzcount);

	// pack the matrix
	pack_storage(mlarge, max_fo<value_type>());
	sort_storage(mlarge);

	nzcount = rows.back();
}

template <class index_type, class value_type>
void transpose_crm(std::vector<index_type>& rows, std::vector<index_type>& cols, std::vector<value_type>& vals, 
               index_type& nr, index_type& nc, index_type& nzcount)
{
    using namespace yasmic;

    std::vector<index_type> rows_temp(rows);
	std::vector<index_type> cols_temp(cols);
	std::vector<value_type> vals_temp(vals);

	std::fill(rows.begin(), rows.end(), 0);

    typedef compressed_row_matrix<
		typename std::vector<index_type>::iterator,
		typename std::vector<index_type>::iterator,
		typename std::vector<value_type>::iterator >
        crs_matrix;

	typedef transpose_matrix<crs_matrix> t_matrix;

	crs_matrix m(rows_temp.begin(), rows_temp.end(), cols_temp.begin(), cols_temp.end(), 
				vals_temp.begin(), vals_temp.end(), nr, nc, nzcount/2);

	t_matrix mt(m);

	nr = nrows(mt);
	nc = ncols(mt);

	// load the matrix
	load_matrix_to_crm(mt, rows.begin(), cols.begin(), vals.begin());
}

/**
 * Build a bipartite graph from a non-square matrix.
 *
 * All parameters are input/output.
 *
 * @param rows the crm rows vector
 * @param cols the crm cols vector
 * @param vals the crm vals vector
 * @param nr the number of rows
 * @param nc the number of columns
 * @param nzcount the number of nonzeros
 */

template <class index_type, class value_type>
void build_bipartite_crm(std::vector<index_type>& rows, std::vector<index_type>& cols, std::vector<value_type>& vals, 
               index_type& nr, index_type& nc, index_type& nzcount)
{
    using namespace yasmic;

    nzcount = (index_type)(2*cols.size());

    std::vector<index_type> rows_temp(rows);
	std::vector<index_type> cols_temp(cols);
	std::vector<value_type> vals_temp(vals);

    std::fill(rows.begin(), rows.end(), 0);
	cols.resize(nzcount);
	vals.resize(nzcount);

    typedef compressed_row_matrix<
		typename std::vector<index_type>::iterator,
		typename std::vector<index_type>::iterator,
		typename std::vector<value_type>::iterator >
        crs_matrix;

	typedef matrix_row_col_graph<crs_matrix> bipartite_graph;

	crs_matrix m(rows_temp.begin(), rows_temp.end(), cols_temp.begin(), cols_temp.end(), 
				vals_temp.begin(), vals_temp.end(), nr, nc, nzcount/2);

	bipartite_graph b(m);

	nr = nr + nc;
	nc = nr;

	// load the matrix
	load_matrix_to_crm(b, rows.begin(), cols.begin(), vals.begin());
}

#if _MSC_VER >= 1400
	// restore the warning for ifstream::read
    #pragma warning( pop )
#endif // _MSC_VER >= 1400

#endif //YASMIC_UTIL_CRM_MATRIX

