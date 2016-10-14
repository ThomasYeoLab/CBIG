#ifndef YASMIC_BOOST_MOD_INTEGER_EXTRA_HPP
#define YASMIC_BOOST_MOD_INTEGER_EXTRA_HPP

/**
 * @file integer_extra.hpp
 * This file fixes some problems with the boost::int_t type
 * and 64-bit compiles under Windows.
 */

#ifdef BOOST_MSVC
// this code is only for MSVC
namespace boost {
    template <>
    struct int_t<64> {
        typedef __int64 least;
        typedef __int64 fast;
    };
}
#endif // BOOST_MSVC

#endif // YASMIC_BOOST_MOD_INTEGER_EXTRA_HPP
