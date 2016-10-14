#ifndef YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_HPP
#define YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_HPP


/*
 * David Gleich
 * Copyright, Stanford University, 2007
 */


/**
 * @file simple_row_and_column_matrix.hpp
 * Implement a simple matrix that allows access to both the rows and columns
 * by storing two copies of the sparsity structure.
 */

/*
 * 10 July 2007
 * Initial version
 * 
 * 29 August 2007
 * Removed big commented section
 */

#include <yasmic/simple_csr_matrix.hpp>

namespace yasmic
{

/**
 * simple_row_and_column_matrix is a user manageable compressed sparse 
 * row matrix type with column access vectors.
 */
template <class IndexType, class ValueType, class NzSizeType=IndexType>
struct simple_row_and_column_matrix
    : public simple_csr_matrix<IndexType,ValueType,NzSizeType>
{
    // additional variables
    NzSizeType* ati;
    IndexType*  atj;
    NzSizeType* atid;
    
    simple_row_and_column_matrix() {}

    simple_row_and_column_matrix(IndexType nrows, IndexType ncols, NzSizeType nnz,
         NzSizeType *ai, IndexType *aj, ValueType *a,
         NzSizeType *ati, IndexType* atj, NzSizeType *atid)
         : simple_csr_matrix<IndexType,ValueType,NzSizeType>(nrows,ncols,nnz,ai,aj,a),
           ati(ati), atj(atj), atid(atid) {}
};

template <class IndexType, class ValueType, class NzSizeType>
struct smatrix_traits< simple_row_and_column_matrix<IndexType, ValueType, NzSizeType> >
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
        : public row_access_tag, public column_access_tag,
          public nonzero_index_tag
    {};
};

template <class IndexType, class ValueType, class NzSizeType>
IndexType nrows(const simple_row_and_column_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.nrows; }

template <class IndexType, class ValueType, class NzSizeType>
IndexType ncols(const simple_row_and_column_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.ncols; }

template <class IndexType, class ValueType, class NzSizeType>
NzSizeType nnz(const simple_row_and_column_matrix<IndexType, ValueType, NzSizeType>& m)
{ return m.nnz; }

/** Build the column access structures for a row access matrix.
 *
 * (r,c) = (r,m.aj[m.ai[r])
 * (r,c) = (m.atj[m.ati[c]],c)
 * ati[m.aj[m.ai[r
 * 
 * @param m the row access matrix
 * @param ati the row pointer for the transpose size=ncols(m)+1
 *   (resp. the column pointer for the matrix)
 * @param atj the column pointer for the transpose size=nnz(m)
 *   (resp. the row index for columns of the matrix)
 * @param atid the global id of each entry in the transpose size=nnz(m)
 */
template <class IndexType, class ValueType, class NzSizeType>
void build_row_and_column_from_csr(
    const simple_csr_matrix<IndexType, ValueType, NzSizeType>& m,
    NzSizeType *ati, IndexType *atj, NzSizeType *atid)
{
    // initialize data to 0
    for (IndexType i=0; i<m.ncols+1; ++i) { ati[i]=0; }
    // compute the column counts
    for (NzSizeType ri=0;ri<m.ai[m.nrows];++ri) { ati[m.aj[ri]+1]++; }
    // compute the cumsum
    for (IndexType i=0,runsum=0;i<m.ncols+1;++i) { runsum=(ati[i]+=runsum); }
    // store the rows and ids
    for (IndexType i=0; i<m.nrows; ++i) {
        for (NzSizeType ri=m.ai[i]; ri<m.ai[i+1]; ++ri) {
            IndexType c=m.aj[ri]; atj[ati[c]]=i; atid[ati[c]]=ri; ++ati[c];
        }
    }
    // restore the column counts
    for (IndexType i=m.ncols; i>0; --i) { ati[i] = ati[i-1]; }
    ati[0] = 0;
}


} // namespace yasmic

#endif // YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_HPP 

