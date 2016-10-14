#ifndef YASMIC_TUPLE_UTILITY
#define YASMIC_TUPLE_UTILITY

#include <boost/tuple/tuple.hpp>

namespace yasmic
{
	template <int N1, int N2, class Tuple>
	boost::tuple<
        typename boost::tuples::element<N1, Tuple>::type,
        typename boost::tuples::element<N2, Tuple>::type >
    tuple_get_2(const Tuple& t) 
    {
    	return ( boost::make_tuple(boost::tuples::get<N1>(t), 
            boost::tuples::get<N2>(t)) );
    }
        
	template <int N1, int N2, class Tuple>
	struct tuple_get_2_fn
	{
		typedef boost::tuple<
			typename boost::tuples::element<N1, Tuple>::type,
			typename boost::tuples::element<N2, Tuple>::type > result_type;
		typedef const Tuple argument_type;

		result_type operator() (argument_type arg) const
		{
			return tuple_get_2<N1, N2, Tuple>(arg);
		}
	};

	template <int N1, class Tuple>
	struct tuple_get_fn
	{
		typedef
			typename boost::tuples::element<N1, Tuple>::type result_type;
		typedef const Tuple argument_type;

		result_type operator() (argument_type arg) const
		{
			return boost::tuples::get<N1>(arg);
		}
	};
}

#endif // YASMIC_TUPLE_UTILITY

