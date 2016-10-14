#ifndef YASMIC_ISTREAM_AS_MATRIX
#define YASMIC_ISTREAM_AS_MATRIX

#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <yasmic/smatrix_traits.hpp>
#include <iterator>

namespace yasmic
{
	template <>
    struct smatrix_traits<std::istream> 
    {
    	typedef unsigned int size_type;
    	typedef int index_type;
		typedef double value_type;

		typedef std::istream_iterator<boost::tuple<index_type, index_type, value_type> > nonzero_iterator;

		typedef void row_iterator;
		typedef void row_nonzero_iterator;
		typedef void column_iterator;
		typedef void column_iterator;
    };
    
    inline std::pair<typename smatrix_traits<std::istream>::size_type,
                      typename smatrix_traits<std::istream>::size_type>
    dimensions(std::istream& f)
    {
    	smatrix_traits<std::istream>::size_type nrows,ncols;
    	
    	f.seekg(0, std::ios_base::beg);
        f >> nrows >> ncols;
        
        return (std::make_pair(nrows, ncols));
    }
    
	inline smatrix_traits<std::istream>::size_type ncols(std::istream& f)
	{
		smatrix_traits<std::istream>::size_type nrows,ncols;
		
		boost::tie(nrows, ncols) = dimensions(f);
		
		return (ncols);
	}
	
	inline smatrix_traits<std::istream>::size_type nrows(std::istream& f)
	{
		smatrix_traits<std::istream>::size_type nrows,ncols;
		
		boost::tie(nrows, ncols) = dimensions(f);
		
		return (nrows);
	}
	
	inline smatrix_traits<std::istream>::size_type nnz(std::istream& f)
	{
		smatrix_traits<std::istream>::size_type d1,d2,nnz;
		
        f.seekg(0, std::ios_base::beg);
        f >> d1 >>  d2 >> nnz;
        
        return (nnz);
	}
	
	inline std::pair<smatrix_traits<std::istream>::nonzero_iterator,
                      smatrix_traits<std::istream>::nonzero_iterator>
    nonzeros(std::ifstream& f)
    {
    	smatrix_traits<std::istream>::size_type d1,d2,d3;
    	f.seekg(0, std::ios_base::beg);
        f >> d1  >> d2 >> d3;
        
        typedef smatrix_traits<std::istream>::nonzero_iterator nz_iter;
        
        return (std::make_pair(nz_iter(f), nz_iter()));
    }
}

#endif //YASMIC_ISTREAM_AS_MATRIX
