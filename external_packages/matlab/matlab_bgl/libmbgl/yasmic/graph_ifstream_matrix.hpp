#ifndef YASMIC_GRAPH_IFSTREAM_MATRIX
#define YASMIC_GRAPH_IFSTREAM_MATRIX

#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <iterator>
#include <string>
#include <sstream>
#include <boost/iterator/iterator_facade.hpp>

#include <yasmic/generic_matrix_operations.hpp>



namespace yasmic
{
	namespace impl
	{
		template <class i_index_type, class i_value_type>
		class graph_ifstream_matrix_const_iterator
		: public boost::iterator_facade<
            graph_ifstream_matrix_const_iterator<i_index_type, i_value_type>,
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const,
            boost::forward_traversal_tag, 
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const >
        {
        public:
            graph_ifstream_matrix_const_iterator() 
				: _str(0), _r(0), _c(0), _line(0)
			{}
            
            graph_ifstream_matrix_const_iterator(std::ifstream &str, std::istringstream &line)
				: _str(&str), _r(0), _c(0), _line(&line), _start(true)
            {
				std::string curline;
				getline( *(_str), curline );

                // set the fail bit so that that increment will
                // discard this line immediately...
                _line->setstate(std::ios_base::failbit);
				
				increment(); 
			}
            
            
        private:
            friend class boost::iterator_core_access;

            void increment() 
            {  
				if (_str != 0)
				{
                    // graph uses 1 indexed columns
                    *(_line) >> _c;
                    --_c;

					while (_line->fail())
					{
						// I don't like the early return here, but 
						// I don't know what to do about it...
						if (_str->eof()) { _str = 0; return; }

						std::string curline;
						getline( *(_str), curline );
						_line->clear();
						_line->str(curline);

                        if (!_start)
                        {
    						++_r;
                        }
                        else
                        {
                            _start = false;
                        }

                        
                        *(_line) >> _c;
                        --_c;
					}
				}
            }
            
            bool equal(graph_ifstream_matrix_const_iterator const& other) const
            {
				return (_str == other._str);
            }
            
            boost::tuple<
                i_index_type, i_index_type, i_value_type>
            dereference() const 
            { 
            	return boost::make_tuple(_r, _c, 1);
            }

			i_index_type _r, _c;

            bool _start;

			std::istringstream* _line;

			std::ifstream* _str;
        };

	}


	template <class index_type = int, class value_type = double, class size_type = unsigned int>
	struct graph_ifstream_matrix
	{
		std::ifstream& _f;
		std::istringstream _line;

		graph_ifstream_matrix(std::ifstream& f)
			: _f(f)
		{
        }
	};

	template <class i_index_type, class i_value_type, class i_size_type>
    struct smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> > 
    {
    	typedef i_size_type size_type;
    	typedef i_index_type index_type;
		typedef i_value_type value_type;
		
		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;

		typedef impl::graph_ifstream_matrix_const_iterator<i_index_type, i_value_type> nonzero_iterator;

		typedef size_type nz_index_type;

		typedef void row_iterator;
		
		typedef void row_nonzero_descriptor;
		typedef void row_nonzero_iterator;
		
		typedef void column_iterator;
    };
    
	template <class i_index_type, class i_value_type, class i_size_type>
    inline std::pair<typename smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type,
                     typename smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type >
    dimensions(graph_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
		typedef smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::size_type nrows,ncols;
    	
    	m._f.clear();
    	m._f.seekg(0, std::ios_base::beg);

	
        m._f >> nrows;
        ncols = nrows;
        
    	        
        return (std::make_pair(nrows, ncols));
    }
	
	template <class i_index_type, class i_value_type, class i_size_type>
	typename smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type
	nnz(graph_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
	{
		typedef smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::index_type d1;
        typename traits::size_type nnz;

        typedef typename traits::size_type size_type;
		
		// clear any error bits
		m._f.clear();
		m._f.seekg(0, std::ios_base::beg);
        
        m._f >> d1 >> nnz;
        return (2*nnz);
	}
	
	template <class i_index_type, class i_value_type, class i_size_type>
	std::pair<typename smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator,
              typename smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator>
    nonzeros(graph_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
    	typedef smatrix_traits<graph_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	//typename traits::size_type ncols;

    	m._f.clear();
		m._f.seekg(0, std::ios_base::beg);
       
        typedef typename traits::nonzero_iterator nz_iter;
        
        return (std::make_pair(nz_iter(m._f, m._line), nz_iter()));

    }
}




#endif //YASMIC_GRAPH_IFSTREAM_MATRIX
