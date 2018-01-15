#ifndef YASMIC_ITERATOR_UTILITY
#define YASMIC_ITERATOR_UTILITY

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

namespace yasmic
{
	/*template <class Type, class SizeType = ptrdiff_t>
	struct constant_iterator_data
	{
		Type _val;
		SizeType _pos;
	};

	template <class Type, class SizeType = ptrdiff_t>
	class constant_iterator
		: public boost::iterator_adaptor<
			constant_iterator<Type, SizeType>,			// Derived
			constant_iterator_data<Type, SizeType>,		// Base type 
			Type,										// Value
			boost::random_access_traversal_tag,			// Traversal
			SizeType>									// Difference
	{
	private:
		struct enabler {};

	public:
		constant_iterator()
			: constant_iterator::iterator_adaptor_() {}

		constant_iterator(Type v)
			: constant_iterator::iterator_adaptor_(
	};*/

	/*template <class Type, class SizeType = ptrdiff_t>
	struct constant_iterator
		: public boost::iterator_facade<
            constant_iterator<Type>,
			Type,
            boost::random_access_traversal_tag, 
            Type >
	{
	public:
		constant_iterator()
		{}

		constant_iterator(Type v)
			: _v(v), _pos(0)
		{}

		constant_iterator(Type v, SizeType pos )
			: _v(v), _pos(pos)
		{}

	private:
		friend class boost::iterator_core_access;

		Type _v;
		SizeType _pos;

		inline void increment() 
		{ ++_pos; }

		inline void decrement()
		{ --_pos;}

		inline void advance(SizeType n)
		{ _pos += n; }

		inline difference_type distance_to(constant_iterator const& other)
		{ return (other._pos - _pos); }

		inline Type dereference() const
		{ return (_v); }

		inline bool equal(constant_iterator const& other) const
		{
			return (other._v == _v && other._pos == _pos);
		}

	};*/

	template <class Type>
	struct constant_iterator
		: public boost::iterator_adaptor<
			constant_iterator<Type>,
			Type*,
			Type,
			boost::random_access_traversal_tag>
	{
		friend class boost::iterator_core_access;

		typedef boost::iterator_adaptor<
			constant_iterator<Type>,
			Type*,
			Type,
			boost::random_access_traversal_tag> super_t;

	public:
		constant_iterator() 
		{}
    
		constant_iterator(constant_iterator const& rhs) : super_t(rhs.base()) 
		{}

		constant_iterator(Type* x)
		: super_t(x)
		{}

	private:
		inline void increment() 
		{ }

		inline void decrement()
		{ }

		inline void advance(typename super_t::difference_type n)
		{ }

		inline typename super_t::difference_type distance_to(constant_iterator const& other)
		{ return (1); }

		inline bool equal(constant_iterator const& other) const
		{
			return (false);
		}
	};


	/*template <class Type, class CountType = unsigned int>
	struct constant_iterator
		: public boost::counting_iterator<CountType>
	{
		friend class iterator_core_access;

	public:
		constant_iterator()
			: super(0)
		{}

		constant_iterator(Type v)
			: super(0), _v(v)
		{}

		constant_iterator(Type v, CountType c)
			: super(c), _v(v)
		{}




	private:
		typedef boost::counting_iterator<CountType> super_t;

		typename super_t::reference dereference()
		{

		}

		Type _v;
	}*/
}

#endif // YASMIC_ITERATOR_UTILITY
