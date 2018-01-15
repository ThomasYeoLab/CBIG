#ifndef YASMIC_UTIL_LOAD_CRM_GRAPH
#define YASMIC_UTIL_LOAD_CRM_GRAPH

/*
 * load_crm_graph.hpp
 * David Gleich
 * Stanford University
 * 26 January 2006
 */

/**
 * @file load_crm_graph.hpp
 * Load a graph or matrix into a crm data structure.
 */
#include <cctype>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>

#include <vector>

#include <boost/iterator/reverse_iterator.hpp>

#include <yasmic/ifstream_matrix.hpp>
#include <yasmic/binary_ifstream_matrix.hpp>
#include <yasmic/cluto_ifstream_matrix.hpp>

#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
	// disable the warning for ifstream::read
    #pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC

/**
 * This function actually loads the data from a matrix file.
 *
 * This function allocates and frees sizeof(Index)*nrows memory.
 */
template <class InputMatrix, class RAIRows, class RAICols, class RAIVals>
bool load_matrix_to_crm(InputMatrix& m, 
						RAIRows rows, RAICols cols, RAIVals vals,
						bool rows_populated = false)
{
	using namespace yasmic;
	using namespace std;

	typedef typename iterator_traits<RAIRows>::value_type index_type;

	typename smatrix_traits<InputMatrix>::size_type nr = nrows(m);
	typename smatrix_traits<InputMatrix>::nz_index_type nzcount = 0;

	typename smatrix_traits<InputMatrix>::nonzero_iterator nzi, nzend;

	int warning = 0;

	if (nr < 1)
	{
		using namespace std;
		cerr << "error: invalid number of rows" << endl;
		return (false);
	}

	if (rows_populated == false)
	{
		boost::tie(nzi, nzend) = nonzeros(m);
		for (; nzi != nzend; ++nzi)
		{
			warning = warning | (row(*nzi, m) >= nr) | (column(*nzi, m) >= nr);
			++rows[row(*nzi, m)*(!warning)+1];
			++nzcount;
		}

		if (warning)
		{
			using namespace std;
			cerr << "error: invalid matrix data, nrows or ncols exceeded" << endl;
			return (false);
		}

		if (nzcount != yasmic::nnz(m))
		{
			using namespace std;
			cerr << "error: number of nonzeros do not match nnz" << endl;
			return (false);
		}
	}

	// compute the reduction
	partial_sum(rows, rows+(nr+1), rows);

	index_type cr;
	index_type cc;

	nzcount = 0;
	tie(nzi, nzend) = yasmic::nonzeros(m);
    for (; nzi != nzend; ++nzi)
    {
		cr = row(*nzi, m);
		cc = column(*nzi, m);
		cols[rows[cr]] = cc;
		vals[rows[cr]] = value(*nzi,m);
		++rows[cr];
		++nzcount;

		warning = warning | (cr >= nr) | (cc >= nr);
	}

	if (warning)
	{
		using namespace std;
		cerr << "error: invalid matrix data, nrows or ncols exceeded" << endl;
		return (false);
	}

	if (nzcount != yasmic::nnz(m))
	{
		using namespace std;
		cerr << "error: number of nonzeros do not match nnz" << endl;
		return (false);
	}

	std::copy(boost::make_reverse_iterator(rows+nr-1),
		boost::make_reverse_iterator(rows), 
		boost::make_reverse_iterator(rows+nr));

	rows[0] = 0;

	return (true);
}

/**
 * This function does most of the work loading the matrix.
 * 
 * 1.  Allocate storage in the passed std::vectors.
 * 2.  Check for degrees metadata and read.
 * 3.  Load the data for the graph.
 */
template <class InputMatrix, class Index, class Value>
bool load_crm_graph_type(InputMatrix& m, std::string filename,
						 std::vector<Index>& rows,
						 std::vector<Index>& cols,
						 std::vector<Value>& vals,
						 Index& nr, Index& nc, Index& nzcount)
{
	using namespace yasmic;
	using namespace std;

	nr = nrows(m);
	nc = ncols(m);
	nzcount = nnz(m);

	//
	// 1.  Allocate storage
	//
	rows.resize(nr+1);
	cols.resize(nzcount);
	vals.resize(nzcount);

	// 
	// 2.  Check for degrees metadata and read.
	//
	bool degrees_data = false;
	{
		// look at the extension...
		typedef std::string::size_type position;

		std::string filename_degrees = filename + ".degs";

		ifstream test(filename_degrees.c_str(), ios::binary);
		if (test.is_open())
		{
			// degrees CAN be binary; we test for a binary file in the following way...
			// 1.  If the filesize is exactly as predicted (sizeof(Index)*nrows)
			// 2.  If the filesize is at least as big as predicted and there 
			//     are binary bytes in the first 100 characters
			ifstream::off_type begin = test.tellg();
			test.seekg (0, ios::end);
			ifstream::off_type end = test.tellg();

			bool binary1 = false;
			if (end-begin == sizeof(Index)*nr)
			{
				// passes test 1
				binary1 = true;
			}

			bool binary2 = false;
			if (!binary1 && (unsigned)(end-begin) >= sizeof(Index)*nr)
			{
				test.seekg (0, ios::end);
				// look for any non ASCII data (it would have to be a list of
				// of ints, so ...)
				for (int i = 0; i < 100; i++)
				{
					char c;
					test.get(c);
					// isspace takes care of newline, tabs, etc.
					if (!std::isdigit(c) && !std::isspace(c))
					{
						binary2 = true;
					}
				}
			}

			test.close();

			if (binary1 || binary2)
			{
				YASMIC_VERBOSE( std::cerr << "reading binary degree file..." << std::endl; )

				// the degrees file is binary
				std::ifstream degfile(filename_degrees.c_str(), ios::binary);

				typename std::vector<Index>::iterator i = rows.begin();
				typename std::vector<Index>::iterator iend = rows.end();

				// we read degrees into the second one
				++i;
				while (i != iend)
				{
					Index ci;
					degfile.read((char *)&ci, sizeof(Index));
					*i = ci;
					++i;
				}

				degrees_data = true;

			}
			else
			{
				YASMIC_VERBOSE( std::cerr << "reading ASCII degree file..." << std::endl; )

				// the degrees file is text
				std::ifstream degfile(filename_degrees.c_str());

				copy(istream_iterator<Index>(degfile), 
					 istream_iterator<Index>(), ++rows.begin());
				
				degrees_data = true;
			}
		}
	}

	// 
	// 3.  Load the matrix
	//
	return (load_matrix_to_crm(m, rows.begin(), cols.begin(), vals.begin(),
		degrees_data));
}


/**
 * Load a CRM graph from a file into a set of vectors.  
 *
 * This operation is actually somewhat tricky.  We to do a few things.
 * 1.  Allocate storage in the passed std::vectors.
 * 2.  Check for degrees metadata and read.
 * 3.  Load the data for the graph.
 *
 * All three operations need an input matrix type, but we need
 * to do the same thing for all, so this function dumps all
 * the work on load_crm_graph_type.
 */
template <class Index, class Value>
bool load_crm_graph(std::string filename, 
					std::vector<Index>& rows, std::vector<Index>& cols,
					std::vector<Value>& vals,
					Index &nr, Index &nc, Index &nzcount)
{
	using namespace std;
	
	// look at the extension...
	typedef string::size_type position;
	
	position dot = filename.find_last_of(".");
	
	if (dot != string::npos)
	{
		string ext = filename.substr(dot+1);
		transform(ext.begin(), ext.end(), ext.begin(), (int(*)(int))tolower);	

		if (ext.compare("smat") == 0 || ext.compare("graph") == 0)
		{
			YASMIC_VERBOSE( std::cerr << "using smat loader..." << std::endl; )

			ifstream ifs(filename.c_str());
			yasmic::ifstream_matrix<> m(ifs);
			return (load_crm_graph_type(m, filename, rows, cols, vals,
						nr, nc, nzcount));
		}
		else if (ext.compare("bsmat") == 0)
		{
			YASMIC_VERBOSE( std::cerr << "using bsmat loader..." << std::endl; )

			ifstream ifs(filename.c_str(), ios::binary);
			yasmic::binary_ifstream_matrix<> m(ifs);
			return (load_crm_graph_type(m, filename, rows, cols, vals,
						nr, nc, nzcount));
		}
        else if (ext.compare("mat") == 0)
		{
			YASMIC_VERBOSE( std::cerr << "cluto mat loader..." << std::endl; )

			ifstream ifs(filename.c_str());
			yasmic::cluto_ifstream_matrix<> m(ifs);
			return (load_crm_graph_type(m, filename, rows, cols, vals,
						nr, nc, nzcount));
		}
		else
		{
			cerr << "Error: matrix type unknown." << endl;
		}
	}
	else
	{
		cerr << "Error: matrix type indeterminate." << endl;
		return (false);
	}

	return (false);
}

#ifdef BOOST_MSVC
#if _MSC_VER >= 1400
	// restore the warning for ifstream::read
	#pragma warning( pop )
#endif // _MSC_VER >= 1400
#endif // BOOST_MSVC

#endif // YASMIC_UTIL_LOAD_CRM_GRAPH

