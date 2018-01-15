#ifndef YASMIC_UNDIR_SIMPLE_CSR_MATRIX_HPP
#define YASMIC_UNDIR_SIMPLE_CSR_MATRIX_HPP

/** @file undir_simple_csr_matrix.hpp
 * @author David F. Gleich
 * @copyright Stanford University, 2008
 * An exact copy of simple_csr_matrix.hpp with a different name for an
 * undirected graph type.
 */

/** History
 *  2008-09-26: Initial coding
 */

#include <yasmic/smatrix_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>

namespace yasmic
{

/**
 * csr_matrix is a more user manageable compressed sparse row matrix type.
 * It is based on the same ideas as the compressed_row_matrix, but designed
 * to be less general and more specific.
 *
 * The idea is that an application will use the csr_matrix structure to
 * manage a sparse matrix and design algorithms that are NOT more
 * generally applicable.
 */
template <class IndexType, class ValueType, class NzSizeType=IndexType>
struct undir_simple_csr_matrix
{
    IndexType nrows;
    IndexType ncols;
    NzSizeType nnz;

    NzSizeType* ai;
    IndexType* aj;
    ValueType* a;

    // an extra indextype to serve as the default value of ai for
    // an empty matrix
    IndexType empty;

    undir_simple_csr_matrix()
        : nrows(0), ncols(0), nnz(0), empty(0), ai(&empty), aj(NULL), a(NULL) {}

    undir_simple_csr_matrix(IndexType nrows, IndexType ncols, NzSizeType nnz,
         NzSizeType *ai, IndexType *aj, ValueType *a)
         : nrows(nrows), ncols(ncols), nnz(nnz), ai(ai), aj(aj), a(a) {}
};


template <class IndexType, class ValueType, class NzSizeType>
struct smatrix_traits< undir_simple_csr_matrix<IndexType, ValueType, NzSizeType> >
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    typedef NzSizeType nz_size_type;
    typedef NzSizeType nz_index_type;

    //typedef typename impl::csr_nonzero<IndexType, ValueType, NzSizeType> nonzero_iterator;
    //typedef impl::csr_nonzero<IndexType, ValueType, NzSizeType> nonzero_descriptor;

    //typedef typename compressed_row_matrix_type::row_nonzero_iterator row_nonzero_iterator;
    //typedef impl::csr_nonzero<IndexType, ValueType, NzSizeType> row_nonzero_descriptor;

    typedef boost::counting_iterator<IndexType> row_iterator;

    struct properties
        : public row_access_tag, public nonzero_index_tag
    {};
};

template <class IndexType, class ValueType, class NzSizeType>
IndexType nrows(const undir_simple_csr_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.nrows; }

template <class IndexType, class ValueType, class NzSizeType>
IndexType ncols(const undir_simple_csr_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.ncols; }

template <class IndexType, class ValueType, class NzSizeType>
NzSizeType nnz(const undir_simple_csr_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.nnz; }


} // namespace yasmic

#endif /* YASMIC_UNDIR_SIMPLE_CSR_MATRIX_HPP */

