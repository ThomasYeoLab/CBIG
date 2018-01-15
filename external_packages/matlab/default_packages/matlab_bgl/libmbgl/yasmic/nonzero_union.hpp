#ifndef YASMIC_NONZERO_UNION
#define YASMIC_NONZERO_UNION

/**
 * @file nonzero_union.hpp
 * @author David Gleich
 * 
 * This class makes a union of two matrices based on the non-zeros.
 * We don't make any copies... 
 *
 * Regardless of the input, the output is just a nonzero-matrix.
 */

namespace yasmic
{

namespace impl
{
	template <class Matrix1, class Matrix2>
	struct nonzero_union_iter_help
	{
		typedef typename smatrix_traits<Matrix1>::nonzero_iterator   nz_iter;
		typedef typename smatrix_traits<Matrix1>::nonzero_descriptor nz_desc;

		typedef typename smatrix_traits<Matrix1>::index_type index_type;
		typedef typename smatrix_traits<Matrix1>::value_type value_type;

		typedef boost::tuple<index_type, index_type, value_type> nz_type;
	};

	// implement the chained iterator
	template <class Matrix1, class Matrix2>
	class nonzero_union_iterator
	: public boost::iterator_facade<
        nonzero_union_iterator<Matrix1, Matrix2>,
		typename nonzero_union_iter_help<Matrix1,Matrix2>::nz_type const,
        boost::forward_traversal_tag, 
        typename nonzero_union_iter_help<Matrix1,Matrix2>::nz_type const >
    {
    public:
        nonzero_union_iterator() 
			: _m1(NULL), _m2(NULL)
		{}
        
        nonzero_union_iterator(Matrix1& m1, Matrix2& m2)
			: _m1(&m1), _m2(&m2), using_m2(false)
        { 
			boost::tie(i,iend) = nonzeros(m1);
			boost::tie(i2,i2end) = nonzeros(m2);
		}
        
        
    private:
        friend class boost::iterator_core_access;

        void increment() 
        {  
			++i;

			if (i == iend)
			{
				if (using_m2)
				{
					_m1 = NULL; 
					_m2 = NULL;
				}
				else
				{
					iend = i2end;
					i = i2;

					using_m2 = true;
				}
			}
        }
        
        bool equal(nonzero_union_iterator const& other) const
        {
			return (_m1 == other._m1 && _m2 == other._m2);
        }
        
        typename nonzero_union_iter_help<Matrix1,Matrix2>::nz_type
        dereference() const 
        { 
			if (!using_m2)
			{
				return boost::make_tuple(
						row(*i, *_m1),
						column(*i, *_m1),
						value(*i, *_m1));
			}
			else
			{
				return boost::make_tuple(
						row(*i, *_m2),
						column(*i, *_m2),
						value(*i, *_m2));
			}
        }

		Matrix1 *_m1;
		Matrix2 *_m2;

		typename nonzero_union_iter_help<Matrix1,Matrix2>::nz_iter i, iend;
		typename nonzero_union_iter_help<Matrix1,Matrix2>::nz_iter i2, i2end;

		bool using_m2;
    };
}

/**
 * The only requirement for using a nonzero union is that all
 * the types must be convertible between the two matrix types
 * and that the matrix types have identical iterator types.
 *
 * We take traits from Matrix1.
 */

template <class Matrix1, class Matrix2> 
class nonzero_union
{
public:
	Matrix1& _m1;
	Matrix2& _m2;

	nonzero_union(Matrix1& m1, Matrix2& m2)
		: _m1(m1), _m2(m2)
	{}

	typedef typename smatrix_traits<Matrix1>::index_type index_type;
	typedef typename smatrix_traits<Matrix1>::value_type value_type;

	typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;
	typedef impl::nonzero_union_iterator<Matrix1, Matrix2> nonzero_iterator;
	//typedef typename smatrix_traits<Matrix1>::nonzero_descriptor nonzero_descriptor;
	//typedef typename smatrix_traits<Matrix1>::nonzero_iterator nonzero_iterator;

	typedef void row_iterator;
	typedef void row_nonzero_descriptor;
	typedef void row_nonzero_iterator;

	typedef void column_iterator;

	typedef typename smatrix_traits<Matrix1>::size_type size_type;
	typedef typename smatrix_traits<Matrix1>::nz_index_type nz_index_type;

    typedef void properties;
};

namespace impl 
{
	template <class Matrix1, class Matrix2>
	struct nonzero_union_help
	{
		typedef smatrix_traits<nonzero_union<Matrix1, Matrix2> > traits;

		typedef std::pair<typename traits::size_type, typename traits::size_type> dims_ret_type;
		typedef std::pair<typename traits::nonzero_iterator, typename traits::nonzero_iterator> nzs_ret_type;
		typedef typename traits::size_type nnz_ret_type;
	};
}


template <class Matrix1, class Matrix2>
inline typename impl::nonzero_union_help<Matrix1, Matrix2>::nnz_ret_type 
nnz(nonzero_union<Matrix1, Matrix2>& tm)
{
	return nnz(tm._m1)+nnz(tm._m2);
}

template <class Matrix1, class Matrix2>
inline typename impl::nonzero_union_help<Matrix1, Matrix2>::dims_ret_type 
dimensions(nonzero_union<Matrix1, Matrix2>& tm)
{
	typename impl::nonzero_union_help<Matrix1, Matrix2>::traits::size_type nr1, nr2;
	typename impl::nonzero_union_help<Matrix1, Matrix2>::traits::size_type nc1, nc2;

	boost::tie(nr1, nc1) = dimensions(tm._m1);
	boost::tie(nr2, nc2) = dimensions(tm._m2);

	return (std::make_pair(std::max(nr1, nr2), std::max(nc1, nc2)));
}

template <class Matrix1, class Matrix2>
inline typename impl::nonzero_union_help<Matrix1, Matrix2>::nzs_ret_type 
nonzeros(nonzero_union<Matrix1, Matrix2>& tm)
{
	return std::make_pair(
		typename impl::nonzero_union_help<Matrix1, Matrix2>::traits::nonzero_iterator(tm._m1, tm._m2), 
		typename impl::nonzero_union_help<Matrix1, Matrix2>::traits::nonzero_iterator());
}

} // namespace yasmic

#endif // YASMIC_NONZERO_UNION

