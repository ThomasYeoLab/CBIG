#ifndef YASMIC_INDEXED_LIST
#define YASMIC_INDEXED_LIST

/**
 * @file indexed_list.hpp
 * Provide an implementation of a indexed_list type.  The indexed_list
 * type satisfies the requirements of LinearOperator and
 * NonZeroList.
 */

#include <yasmic/smatrix_traits.hpp>

namespace yasmic
{
	/*
	 * Implementation specific classes.
	 */
	namespace impl
	{
		
	};

	/**
	 * indexed_list provides a lightweight wrapper 
	 * implementation of NumericMatrix and NonzeroIterable.
	 *
	 * @implements NumericMatrix, NonzeroIterable
	 *
	 * The best way to use this class is through the make_indexed_list function
	 * which provides a set of useful wrapper functions. 
	 *
	 * To get the correct type for the indexed_list (IF YOU NEED THIS, you may not; 
	 * see the documentation on "when do I need to generate a type") use
	 * the indexed_list_type_generator interface.  
	 */
	template <class NonzeroIterator>
	struct indexed_list
	{
	public:
        
		// required types
		typedef typename NonzeroIterator::size_type size_type;
		typedef size_type index_type;
		typedef typename NonzeroIterator::value_type value_type;

		typedef NonzeroIterator nonzero_iterator;

		/**
		 * Basic constructor.  This constructor requires a pass over the dataset
		 * to compute the basic statistics such as nrows and ncols.  
		 */
		indexed_list(nonzero_iterator start, nonzero_iterator end)
			: _start(start), _end(end), _nrows(0), _ncols(0), _nnz(0)
		{
			nonzero_iterator i = _start;

			_nnz = 0;
			
			while (i != _end)
			{
				_nrows = std::max(get<0>(*i),_nrows);
				_ncols = std::max(get<1>(*i),_ncols);

				++i;
				++_nnz;
			}

			// we underestimate by one value
			++_nrows;
			++_ncols;
		}

		/**
		 * Constructor which requires all the data but involves no compute time.
		 */
		indexed_list(nonzero_iterator start, nonzero_iterator end,
			size_type nrows, size_type ncols, size_type nnz)
			: _start(start), _end(end), _nrows(nrows), _ncols(ncols), _nnz(nnz)
		{}

	private:
		size_type _nrows;
		size_type _ncols;
		size_type _nnz;

		nonzero_iterator _start;
		nonzero_iterator _end;
	};

	template <class Val, class Collection>
	inline typename Collection::size_type nrows(const indexed_list<Val, Collection>& m)
	{
		return m.nrows();
	}

	template <class Val, class Collection>
	inline typename Collection::size_type ncols(const indexed_list<Val, Collection>& m)
	{
		return m.ncols();
	}
};

#endif // YASMIC_INDEXED_LIST