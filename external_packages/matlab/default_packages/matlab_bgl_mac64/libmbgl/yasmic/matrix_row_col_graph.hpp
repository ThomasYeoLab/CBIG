#ifndef YASMIC_MATRIX_ROW_COL_GRAPH
#define YASMIC_MATRIX_ROW_COL_GRAPH

/**
 * @file matrix_row_col_graph.hpp
 * @author David Gleich
 * 
 * A row and column graph filter for a matrix. 
 *
 * Regardless of the input, the output is just a nonzero-matrix.
 *
 * The row and column graph of a matrix has a vertex for each row and column 
 * and an edge where A(i,j) != 0.  The weight on the edge is A(i,j).
 *
 * The first m vertices correspond to the row vertices, the next
 * n vertices correspond to the columns.
 *
 * We just wrap the original dataset and forward calls.  
 * (This does not make a copy of the original matrix!)
 */

namespace yasmic
{

namespace impl
{
    template <class Matrix>
	struct matrix_row_col_graph_nz_iter_help
	{

		typedef typename smatrix_traits<Matrix>::index_type index_type;
		typedef typename smatrix_traits<Matrix>::value_type value_type;

		typedef boost::tuple<index_type, index_type, value_type> nz_type;
	};

    // implement the "repeating" iterator
	template <class Matrix>
	class matrix_row_col_graph_nz_iterator
	: public boost::iterator_facade<
        matrix_row_col_graph_nz_iterator<Matrix>,
		typename matrix_row_col_graph_nz_iter_help<Matrix>::nz_type const,
        boost::forward_traversal_tag, 
        typename matrix_row_col_graph_nz_iter_help<Matrix>::nz_type const >
    {

    public:
        matrix_row_col_graph_nz_iterator() 
			: _m(NULL), _repeat(true), _nr(0)
		{}
        
        matrix_row_col_graph_nz_iterator(Matrix& m)
			: _m(&m), _repeat(true), _nr(nrows(m))
        { 
			i = nonzeros(m).first;
		}

        matrix_row_col_graph_nz_iterator(Matrix& m,bool)
			: _m(&m), _repeat(true), _nr(nrows(m))
        { 
			i = nonzeros(m).second;
		}
        
        
    private:
        friend class boost::iterator_core_access;

        void increment() 
        {  
            if (_repeat)
            {
                _repeat = false;
            }
            else
            {
                _repeat = true;
                ++i;
            }
        }
        
        bool equal(matrix_row_col_graph_nz_iterator const& other) const
        {
			return (i == other.i && _repeat == other._repeat);
        }
        
        typename matrix_row_col_graph_nz_iter_help<Matrix>::nz_type
        dereference() const 
        { 
            if (_repeat)
            {
                return boost::make_tuple(
                    row(*i, *_m),column(*i, *_m)+_nr,value(*i,*_m));
            }
            else
            {
                return boost::make_tuple(
                   column(*i, *_m)+_nr,row(*i, *_m),value(*i,*_m));
            }
        }

        bool _repeat;

		Matrix *_m;
		
		typename smatrix_traits<Matrix>::nonzero_iterator i;
        typename smatrix_traits<Matrix>::index_type _nr;
    };


}  /* end namespace yasmic::impl */


template <class Matrix>
class matrix_row_col_graph
{
public:
	Matrix& _m;

	matrix_row_col_graph(Matrix& m)
		: _m(m)
	{}

	typedef typename smatrix_traits<Matrix>::index_type index_type;
	typedef typename smatrix_traits<Matrix>::value_type value_type;

    typedef typename boost::tuple<index_type, index_type, value_type> nonzero_descriptor;
    typedef impl::matrix_row_col_graph_nz_iterator<Matrix> nonzero_iterator;

	typedef void row_iterator;
	typedef void row_nonzero_descriptor;
	typedef void row_nonzero_iterator;

	typedef void column_iterator;

	typedef typename smatrix_traits<Matrix>::size_type size_type;
	typedef typename smatrix_traits<Matrix>::nz_index_type nz_index_type;
	typedef typename smatrix_traits<Matrix>::properties properties;
};

namespace impl 
{
	template <class Matrix>
	struct matrix_row_col_graph_help
	{
		typedef smatrix_traits<matrix_row_col_graph<Matrix> > traits;

		typedef std::pair<typename traits::size_type, typename traits::size_type> dims_ret_type;
		typedef std::pair<typename traits::nonzero_iterator, typename traits::nonzero_iterator> nzs_ret_type;
		typedef typename traits::size_type nnz_ret_type;
	};
}


template <class Matrix>
inline typename impl::matrix_row_col_graph_help<Matrix>::nnz_ret_type 
nnz(matrix_row_col_graph<Matrix>& tm)
{
	return 2*nnz(tm._m);
}

template <class Matrix>
inline typename impl::matrix_row_col_graph_help<Matrix>::dims_ret_type 
dimensions(matrix_row_col_graph<Matrix>& tm)
{
	typename impl::transpose_matrix_help<Matrix>::dims_ret_type d = dimensions(tm._m);
	return (std::make_pair(d.second + d.first, d.second + d.first));
}

template <class Matrix>
inline typename impl::matrix_row_col_graph_help<Matrix>::nzs_ret_type 
nonzeros(matrix_row_col_graph<Matrix>& tm)
{
 
    return std::make_pair( impl::matrix_row_col_graph_nz_iterator<Matrix>(tm._m),
        impl::matrix_row_col_graph_nz_iterator<Matrix>(tm._m,true));
}

} // namespace yasmic

#endif // YASMIC_MATRIX_ROW_COL_GRAPH


