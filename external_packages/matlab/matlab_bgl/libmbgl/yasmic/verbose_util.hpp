/*
 * verbose_util.hpp
 * David Gleich
 * Stanford University
 * 17 March 2006
 */

/**
 * @file verbose_util.hpp
 * Add defines to control verbose statement in the code.
 */

//
// this file is somewhat strange.  There are two sections; the first section 
// simple makes a few define statements
//
#ifdef YASMIC_VERBOSE_UTIL_DEFINE
namespace yasmic
{
    int yasmic_verbose;
} // namespace yasmic
#else
namespace yasmic
{
    extern int yasmic_verbose;
}
#endif // YASMIC_VERBOSE_UTIL_DEFINE

#ifndef YASMIC_VERBOSE_UTIL
#define YASMIC_VERBOSE_UTIL

#ifdef YASMIC_USE_VERBOSE

	#define YASMIC_VERBOSE(a) \
        if (yasmic::yasmic_verbose) \
        { \
            a; \
        } 

#else

	#define YASMIC_VERBOSE(a) \

#endif 




#endif // YASMIC_VERBOSE_UTIL

