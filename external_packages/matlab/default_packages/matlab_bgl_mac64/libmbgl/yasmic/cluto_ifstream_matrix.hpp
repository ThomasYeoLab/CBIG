#ifndef YASMIC_CLUTO_IFSTREAM_MATRIX
#define YASMIC_CLUTO_IFSTREAM_MATRIX

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
#include <string>
#include <sstream>
#include <boost/iterator/iterator_facade.hpp>

#include <yasmic/generic_matrix_operations.hpp>

#include <boost/lexical_cast.hpp>


namespace yasmic
{
	namespace impl
	{
        /*template <class i_index_type, class i_value_type>
        class dense_cluto_ifstream_matrix_const_iterator
		: public boost::iterator_facade<
            dense_cluto_ifstream_matrix_const_iterator<i_index_type, i_value_type>,
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const,
            boost::forward_traversal_tag, 
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const >
        {
        public:
            dense_cluto_ifstream_matrix_const_iterator() 
				: _str(0), _r(0), _c(0), _v(0)
			{}
            
            dense_cluto_ifstream_matrix_const_iterator(std::ifstream &str)
				: _str(&str), _r(0), _c(0), _v(0)
            {
				*(_str) >> _v;
			}
            
            
        private:
            friend class boost::iterator_core_access;

            void increment() 
            {  
				if (_str != 0)
				{
					*(_str) >> _v;

                    if (_str->eof())
                    {
                        _str = 0;
                        return;
                    }

                    ++_c;

                    if (_c >= _ncols)
                    {
                        _c = 0;
                        ++_r;
                    }
				}
            }
            
            bool equal(dense_cluto_ifstream_matrix_const_iterator const& other) const
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

			index_type _ncols;

			std::ifstream* _str;
        };*/

		template <class i_index_type, class i_value_type>
		class cluto_ifstream_matrix_const_iterator
		: public boost::iterator_facade<
            cluto_ifstream_matrix_const_iterator<i_index_type, i_value_type>,
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const,
            boost::forward_traversal_tag, 
            boost::tuple<
                i_index_type, i_index_type, i_value_type> const >
        {
        public:
            cluto_ifstream_matrix_const_iterator() 
				: _str(0), _r(0), _c(0), _v(0), _line(0), _dense(false), _start(false)
			{}
            
            cluto_ifstream_matrix_const_iterator(std::ifstream &str, std::istringstream &line, bool dense)
				: _str(&str), _r(0), _c(0), _v(0), _line(&line), _dense(dense), _start(true)
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
                    if (_dense)
                    {
                        *(_line) >> _v;
                        ++_c;
                    }
                    else
                    {
                        // cluto uses 1 indexed columns
                        *(_line) >> _c;
                        *(_line) >> _v;
                        --_c;
                    }

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

                        if (_dense)
                        {
                            *(_line) >> _v;
                            _c = 0;
                        }
                        else
                        {
                            // cluto uses 1 indexed columns
                            *(_line) >> _c;
                            *(_line) >> _v;
                            --_c;
                        }
					}
				}
            }
            
            bool equal(cluto_ifstream_matrix_const_iterator const& other) const
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

            bool _dense;
            bool _start;

			std::istringstream* _line;

			std::ifstream* _str;
        };

	}


	template <class index_type = int, class value_type = double, class size_type = unsigned int>
	struct cluto_ifstream_matrix
	{
		std::ifstream& _f;
		std::istringstream _line;

        bool _graph;
        bool _dense;

		cluto_ifstream_matrix(std::ifstream& f)
			: _f(f), _graph(false), _dense(false)
		{
            detect_graph_and_dense();
        }

        cluto_ifstream_matrix(std::ifstream& f, bool graph)
			: _f(f), _graph(graph), _dense(false)
		{
            detect_dense();
        }

        cluto_ifstream_matrix(std::ifstream& f, bool graph, bool dense)
			: _f(f), _graph(graph), _dense(dense)
		{
        }

    private:
        /**
         * matrix-type detection algorithm
         *
         * 1.  if there is only one or three tokens on the first line, then,
         *     one token => dense graph
         *     three tokens => sparse matrix
         * 2.  if the first line of the file contains fewer than
         *     ncols tokens => sparse graph
         * 3.  if the first line contains invalid integers at
         *     even locations => dense graph
         * 4.  repeat for all lines...
         */

        void detect_graph_and_dense()
        {
            std::string line;
            getline( _f, line );
            _line.clear(); _line.str(line);

            int tok_count=0;


            while (!_line.eof())
            {
                // these entries are either index_type or
                // size_type, but size_type should be larger
                // than index type...
                size_type st;
                _line >> st;
                tok_count++;
            }

            if (_line.fail())
            {
                tok_count--;
            }

            if (tok_count == 1)
            {
                _dense = true;
                _graph = true;
                return;
            }
            else if (tok_count == 3)
            {
                _dense = false;
                _graph = false;
                return;
            }

            _f.clear();
    	    _f.seekg(0, std::ios_base::beg);

            index_type maybe_nrows;
            size_type maybe_ncols;

            _f >> maybe_nrows;
            _f >> maybe_ncols;

            getline( _f, line );

            while (_f)
            {
                getline( _f, line );
                _line.clear(); _line.str(line);

                if (detect_graph_and_dense_check_line(maybe_nrows, maybe_ncols))
                {
                    return;
                }
            }

            // if we've gotten all the way through the matrix and 
            // there still isn't anything left...

            _dense = true;
            _graph = false;

            
        }

        void detect_dense()
        {
            std::string line;
            getline( _f, line );
            _line.clear(); _line.str(line);

            int tok_count=0;

            while (_line)
            {
                size_type st;
                _line >> st;
                tok_count++;
            }

            if (tok_count == 1)
            {
                _dense = true;
            }
            else 
            {
                _dense = false;
            }

            _f.clear();
    	    _f.seekg(0, std::ios_base::beg);

            return;
        }

        /**
         * Checks the current _line to determine if the
         * input is dense or sparse.
         *
         * @return true if detection occured, false otherwise
         */
        bool detect_graph_and_dense_check_line(index_type maybe_nrows, index_type maybe_ncols)
        {
            using boost::lexical_cast;
            using boost::bad_lexical_cast;


            int tok_count = 0;
            std::string tok;

            while (!_line.eof())
            {
                /*index_type i;
                
                _line >> i;

                if (_line.eof())
                {
                    // the correct behavior will happen if we 
                    // just wrap around...
                }
                else if (!_line.fail())
                {
                    if (i < 1 || i > maybe_ncols)
                    {
                        _dense = true;
                        _graph = false;
                        return (true);
                    }
                    // we always have to be
                    // able to parse a value...
                    value_type v;                  

                    _line >> v;
                    tok_count+=2;
                }
                else
                {
                    _dense = true;
                    _graph = false;
                    return (true);
                }*/
              
                try
                {
                    _line >> tok;
                    index_type i = lexical_cast<index_type>(tok);
                    

                    if (i < 1 || i > maybe_ncols)
                    {
                        _dense = true;
                        _graph = false;
                        return (true);
                    }

                    _line >> tok;
                    value_type v = lexical_cast<value_type>(tok);

                    tok_count += 2;
                }
                catch (bad_lexical_cast e)
                {
                    _dense = true;
                    _graph = false;
                    return (true);
                }
            }

            if (tok_count == maybe_ncols)
            {
                return (false);
            }
            else
            {
                _dense = false;
                _graph = true;
                return (true);
            }
        }
	};

	template <class i_index_type, class i_value_type, class i_size_type>
    struct smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> > 
    {
    	typedef i_size_type size_type;
    	typedef i_index_type index_type;
		typedef i_value_type value_type;
		
		typedef boost::tuple<index_type, index_type, value_type> nonzero_descriptor;

		typedef impl::cluto_ifstream_matrix_const_iterator<i_index_type, i_value_type> nonzero_iterator;

		typedef size_type nz_index_type;

		typedef void row_iterator;
		
		typedef void row_nonzero_descriptor;
		typedef void row_nonzero_iterator;
		
		typedef void column_iterator;
    };
    
	template <class i_index_type, class i_value_type, class i_size_type>
    inline std::pair<typename smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type,
                     typename smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type >
    dimensions(cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
		typedef smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::size_type nrows,ncols;
    	
    	m._f.clear();
    	m._f.seekg(0, std::ios_base::beg);

		if (m._graph)
        {
            m._f >> nrows;
            ncols = nrows;
        }
        else
        {
            m._f >> nrows >> ncols;
        }
    	        
        return (std::make_pair(nrows, ncols));
    }
	
	template <class i_index_type, class i_value_type, class i_size_type>
	typename smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::size_type
	nnz(cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
	{
		typedef smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	typename traits::index_type d1, d2;
        typename traits::size_type nnz;

        typedef typename traits::size_type size_type;
		
		// clear any error bits
		m._f.clear();
		m._f.seekg(0, std::ios_base::beg);

        if (m._dense)
        {
            if (m._graph)
            {
                m._f >> nnz;
                return (nnz*nnz);
            }
            else
            {
                m._f >> d1 >> d2;
                return (((size_type)d1)*((size_type)d2));
            }
        }
        else
        {
            if (m._graph)
            {
                m._f >> d1 >> nnz;
                return (nnz);
            }
            else
            {
                m._f >> d1 >> d2 >> nnz;
                return (nnz);
            }
        }
	}
	
	template <class i_index_type, class i_value_type, class i_size_type>
	std::pair<typename smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator,
              typename smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> >::nonzero_iterator>
    nonzeros(cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type>& m)
    {
    	typedef smatrix_traits<cluto_ifstream_matrix<i_index_type, i_value_type, i_size_type> > traits;
    	//typename traits::size_type ncols;

    	m._f.clear();
		m._f.seekg(0, std::ios_base::beg);


        //getline(m._f, header_line);
		//m._f >> d1 >> d2 >> d3;

		// seek after nrows,ncols,nnz
    	//m._f.seekg(2*sizeof(i_index_type)+sizeof(i_size_type), std::ios_base::beg);
        
        typedef typename traits::nonzero_iterator nz_iter;
        

        return (std::make_pair(nz_iter(m._f, m._line, m._dense), nz_iter()));

    }
}


#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
	// restore the warning for ifstream::read
	#pragma warning( pop )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC


#endif //YASMIC_CLUTO_IFSTREAM_MATRIX
