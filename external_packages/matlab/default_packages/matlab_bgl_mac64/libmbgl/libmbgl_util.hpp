/** @file libmbgl_util.hpp
 * @author David F. Gleich
 * @date 2008-09-25
 * @copyright Stanford University, 2008
 * Misc utility functions
 */

/** History
 *  2008-09-25: Initial coding
 */

/** Wrapper class to turn a row-storage matrix into a boost BasicMatrix
 */
template <class value_type>
struct row_matrix {
  value_type* head;
  std::size_t nrows;
	std::size_t ncols;
	const value_type* operator[](std::size_t i) const {
		return (head+ncols*i);
  }

	value_type* operator[](std::size_t i) {
		return (head+ncols*i);
	}

  row_matrix(value_type* _head, std::size_t _nrows, std::size_t _ncols)
		: head(_head), nrows(_nrows), ncols(_ncols) {}
};
