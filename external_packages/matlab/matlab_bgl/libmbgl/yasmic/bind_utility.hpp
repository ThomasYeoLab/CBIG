#ifndef YASMIC_BIND_UTILITY
#define YASMIC_BIND_UTILITY

#include <functional>

/*
 * taken from http://www.boost.org/libs/iterator/example/transform_iterator_example.cpp
 */

// What a bummer. We can't use std::binder1st with transform iterator
// because it does not have a default constructor. Here's a version
// that does.

namespace boost {

  template <class Operation> 
  class binder1st
    : public std::unary_function<typename Operation::second_argument_type,
                                 typename Operation::result_type> {
  protected:
    Operation op;
    typename Operation::first_argument_type value;
  public:
    binder1st() { } // this had to be added!
    binder1st(const Operation& x,
              const typename Operation::first_argument_type& y)
        : op(x), value(y) {}
    typename Operation::result_type
    operator()(const typename Operation::second_argument_type& x) const {
      return op(value, x); 
    }
  };

  template <class Operation, class T>
  inline binder1st<Operation> bind1st(const Operation& op, const T& x) {
    typedef typename Operation::first_argument_type arg1_type;
    return binder1st<Operation>(op, arg1_type(x));
  }

} // namespace boost

#endif // YASMIC_BIND_UTILITY