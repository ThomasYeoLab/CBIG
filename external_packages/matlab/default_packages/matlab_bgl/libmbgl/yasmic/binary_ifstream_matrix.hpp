#ifndef YASMIC_BINARY_IFSTREAM_MATRIX
#define YASMIC_BINARY_IFSTREAM_MATRIX

#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
	// disable the warning for ifstream::read
	#pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC


#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>

#include <yasmic/generic_matrix_operations.hpp>


namespace yasmic
{
	namespace impl
	{
	
		template <class i_index_type, class i_value_type>
		class binary_ifstream_matrix_const_iterator
		: public boost::iterator_facade<
            binary_ifstream_matrix_const_iterator<i_index_type, i_value_type>,
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const,
            boost::forward_traversal_tag, 
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const >
        {
        public:
            binary_ifstream_matrix_const_iterator() 
				: _str(0), _r(0), _c(0), _v(0)
			{}
            
            binary_ifstream_matrix_const_iterator(std::istream &str)
				: _str(&str), _r(0), _c(0), _v(0)
            { increment(); }
            
            
        private:
            friend class boost::iterator_core_access;

            void increment() 
            {  
				if (_str != 0)
				{
            		_str->read((char *)&_r, sizeof(i_index_type));
					_str->read((char *)&_c, sizeof(i_index_type));
					_str->read((char *)&_v, sizeof(i_value_type));

					if (_str->eof()) { _str = 0; }
				}
            }
            
            bool equal(binary_ifstream_matrix_const_iterator const& other) const
            {
				return (_str == other._str);
            }
            
            boost::tuple<
                i_index_type, i_index_type, i_value_type>
            dereference() const 
            { 
            	return boost::make_tuple(_r, _c, _v);
            }

			i_index_type _r, _c;
			i_value_type _v;

			std::istream* _str;
        };

	}


	template <class index_type = int, class value_type = double, class size_type = int>
	struct binary_ifstream_matrix
	{
		std::istream& _f;

		binary_ifstream_matrix(std::istream& f)
			: _f(f) 
		{}

		binary_ifstream_matrix(const binary_ifstream_matrix& bifm)
			: _f(bifm._f)
		{}
	};

	template <class i_index_type, class i_value_type, class i_size_type>
    struct smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> > 
    {
    	typedef i_size_type size_type;
    	typedef i_index_type index_type;
		typedef i_value_type value_type;
		
		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;

		typedef impl::binary_ifstream_matrix_const_iterator<i_index_type, i_value_type> nonzero_iterator;

		typedef size_type nz_index_type;

		typedef void row_iterator;
		
		typedef void row_nonzero_descriptor;
		typedef void row_nonzero_iterator;
    };
    
	template <class i_index_type, class i_value_type, class i_size_type>
    inline std::pair<typename smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type,
                     typename smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type >
    dimensions(binary_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
		typedef smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::size_type nrows,ncols;
    	
    	m._f.clear();
    	m._f.seekg(0, std::ios_base::beg);
    	
        m._f.read((char *)&nrows, sizeof(i_index_type));
		m._f.read((char *)&ncols, sizeof(i_index_type));
        
        return (std::make_pair(nrows, ncols));
    }
	
	template <class i_index_type, class i_value_type, class i_size_type>
	typename smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type
	nnz(binary_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
	{
		typedef smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::size_type nnz;
		
		// clear any error bits
		m._f.clear();
		
		// seek after nrows, ncols
        m._f.seekg(2*sizeof(i_index_type), std::ios_base::beg);

        m._f.read((char *)&nnz, sizeof(i_size_type));
        
        return (nnz);
	}
	
	template <class i_index_type, class i_value_type, class i_size_type>
	std::pair<typename smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator,
              typename smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator>
    nonzeros(binary_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
    	typedef smatrix_traits<binary_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	
    	m._f.clear();

		// seek after nrows,ncols,nnz
    	m._f.seekg(2*sizeof(i_index_type)+sizeof(i_size_type), std::ios_base::beg);
        
        typedef typename traits::nonzero_iterator nz_iter;
        
        return (std::make_pair(nz_iter(m._f), nz_iter()));
    }
}


#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
	// restore the warning for ifstream::read
	#pragma warning( pop )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC


#endif //YASMIC_BINARY_IFSTREAM_MATRIX
