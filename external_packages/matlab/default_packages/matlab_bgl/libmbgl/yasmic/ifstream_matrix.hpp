#ifndef YASMIC_IFSTREAM_MATRIX
#define YASMIC_IFSTREAM_MATRIX

#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp> 
#include <boost/static_assert.hpp>
#include <yasmic/smatrix_traits.hpp>
#include <iterator>

#include <yasmic/generic_matrix_operations.hpp>

namespace yasmic
{
	template <class index_type = int, class value_type = double, class size_type = int,
		bool header = true>
	struct ifstream_matrix
	{
		ifstream_matrix(std::istream& f)
			: _f(f), _nrows(0), _ncols(0), _nnz(0), _sized(false) 
		{
			// this call is only valid if we have to read the header
			BOOST_STATIC_ASSERT(header == true);
		}

		ifstream_matrix(std::istream& f, index_type nrows, index_type ncols, size_type nnz)
			: _f(f), _nrows(nrows), _ncols(ncols), _nnz(nnz), _sized(true)
		{}

		std::istream& _f;
		size_type _nrows;
		size_type _ncols;
		size_type _nnz;

		bool _sized;
	};

	template <class i_index_type, class i_value_type, class i_size_type, bool header>
    struct smatrix_traits<ifstream_matrix<i_index_type, i_value_type, i_size_type, header> > 
    {
    	typedef i_size_type size_type;
    	typedef i_index_type index_type;
		typedef i_value_type value_type;
		
		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;

		typedef std::istream_iterator<nonzero_descriptor> nonzero_iterator;

		typedef i_size_type nz_index_type;

		typedef void row_iterator;
		
		typedef void row_nonzero_descriptor;
		typedef void row_nonzero_iterator;
		
		typedef void column_iterator;

        typedef void properties;
    };

	namespace impl
	{
		template <class i_index_type, class i_value_type, class i_size_type, bool header>
		struct ifstream_matrix_help
		{
			typedef ifstream_matrix<i_index_type, i_value_type, i_size_type, header> mat_type;
			typedef smatrix_traits<mat_type> mat_traits;
			typedef std::pair<typename mat_traits::size_type, 
							  typename mat_traits::size_type> dimensions_ret_type;
			typedef std::pair<typename mat_traits::nonzero_iterator,
							  typename mat_traits::nonzero_iterator> nz_ret_type;
		};

	}
	
    
	template <class i_index_type, class i_value_type, class i_size_type, bool header>
	inline typename
		impl::ifstream_matrix_help<i_index_type, i_value_type, i_size_type, header>::dimensions_ret_type
    dimensions(ifstream_matrix<i_index_type, i_value_type, i_size_type, header>& m)
    {   	
		if (header && !m._sized)
		{
			m._f.clear();
    		m._f.seekg(0, std::ios_base::beg);
    	
			m._f >> m._nrows >> m._ncols >> m._nnz;
			m._sized = true;
		}
    	
        
        return (std::make_pair(m._nrows, m._ncols));
    }
	
	template <class i_index_type, class i_value_type, class i_size_type, bool header>
	inline typename 
		impl::ifstream_matrix_help<i_index_type, i_value_type, i_size_type, header>::mat_traits::size_type
	nnz(ifstream_matrix<i_index_type, i_value_type, i_size_type, header>& m)
	{		
		if (header && !m._sized)
		{
			m._f.clear();
    		m._f.seekg(0, std::ios_base::beg);
    	
			m._f >> m._nrows >> m._ncols >> m._nnz;
			m._sized = true;
		}
        
        return (m._nnz);
	}
	
	template <class i_index_type, class i_value_type, class i_size_type, bool header>
	inline typename
		impl::ifstream_matrix_help<i_index_type, i_value_type, i_size_type, header>::nz_ret_type
    nonzeros(ifstream_matrix<i_index_type, i_value_type, i_size_type, header>& m)
    {
    	typedef smatrix_traits<ifstream_matrix<i_index_type, i_value_type, i_size_type, header> > traits;
    	
    	m._f.clear();
    	m._f.seekg(0, std::ios_base::beg);
    	
		if (header)
		{
    		typename traits::size_type d1,d2,d3;
			m._f >> d1  >> d2 >> d3;
		}
        
        m._f >> boost::tuples::set_open(' ') >> boost::tuples::set_close(' ');
        
        typedef typename traits::nonzero_iterator nz_iter;
        
        return (std::make_pair(nz_iter(m._f), nz_iter()));
    }

    

}

#endif //YASMIC_IFSTREAM_MATRIX
