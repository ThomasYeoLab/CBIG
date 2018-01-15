#ifndef YASMIC_IFSTREAM_AS_MATRIX
#define YASMIC_IFSTREAM_AS_MATRIX

#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp> 
#include <yasmic/smatrix_traits.hpp>
#include <iterator>

#include <yasmic/generic_matrix_operations.hpp>

namespace yasmic
{
	template <>
    struct smatrix_traits<std::ifstream> 
    {
    	typedef unsigned int size_type;
    	typedef int index_type;
		typedef double value_type;
		
		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;

		typedef std::istream_iterator<nonzero_descriptor> nonzero_iterator;

		typedef void row_iterator;
		
		typedef void row_nonzero_descriptor;
		typedef void row_nonzero_iterator;
		
		typedef void column_iterator;
    };
    
    inline std::pair<smatrix_traits<std::ifstream>::size_type,
                      smatrix_traits<std::ifstream>::size_type>
    dimensions(std::ifstream& f)
    {
    	smatrix_traits<std::ifstream>::size_type nrows,ncols;
    	
    	f.clear();
    	f.seekg(0, std::ios_base::beg);
    	
        f >> nrows >> ncols;
        
        return (std::make_pair(nrows, ncols));
    }
	
	smatrix_traits<std::ifstream>::size_type nnz(std::ifstream& f)
	{
		smatrix_traits<std::ifstream>::size_type d1,d2,nnz;
		
		// clear any error bits
		f.clear();
        f.seekg(0, std::ios_base::beg);
        
        f >> d1 >> d2 >> nnz;
        
        return (nnz);
	}
	
	inline std::pair<smatrix_traits<std::ifstream>::nonzero_iterator,
                      smatrix_traits<std::ifstream>::nonzero_iterator>
    nonzeros(std::ifstream& f)
    {
    	smatrix_traits<std::ifstream>::size_type d1,d2,d3;
    	
    	f.clear();
    	f.seekg(0, std::ios_base::beg);
    	
        f >> d1  >> d2 >> d3;
        
        f >> boost::tuples::set_open(' ') >> boost::tuples::set_close(' ');
        
        typedef smatrix_traits<std::ifstream>::nonzero_iterator nz_iter;
        
        return (std::make_pair(nz_iter(f), nz_iter()));
    }
}

#endif //YASMIC_IFSTREAM_AS_MATRIX
