/** @file components.cc
 * @author David F. Gleich
 * @date 2006-04-19
 * @copyright Stanford University, 2006-2008
 * Implement wrappers for connected component functions
 */

/** History
 *  2006-04-19: Initial version
 *  2007-07-09: Updated to use simple_csr_matrix graph type
 */

#include "include/matlab_bgl.h"

#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/iterator_utility.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/strong_components.hpp>

int strong_components(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* ci)
{
  using namespace yasmic;
  using namespace boost;

  typedef simple_csr_matrix<mbglIndex, double> crs_graph;
  crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

  strong_components(g, make_iterator_property_map(ci, get(vertex_index, g)));

  return 0;
}

/**
 * Wrap a boost graph library call to biconnected_components.
 *
 * the ja and ia arrays specify the connectivity of the underlying graph,
 * ia is a length (nverts+1) array with the indices in ja that start the
 * nonzeros in each row.  ja is a length (ia(nverts)) array with the
 * columns of the connectivity.
 *
 * if a or ci is NULL, then that parameter is not computed.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param a an array which will store the articulaion points of the graph
 *     the array length should be n
 * @param ci the component index array which is length (nnz)
 * @return an error code if possible
 */

int biconnected_components(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* a, mbglIndex* ci)
{
  using namespace yasmic;
  using namespace boost;

  typedef simple_csr_matrix<mbglIndex, double> crs_graph;
  crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

  if (a) {
    if (ci) {
      std::size_t num_bicomps;
      mbglIndex *oi;
      boost::tie(num_bicomps, oi) = biconnected_components(g,
          make_iterator_property_map(ci, get(edge_index, g)), a);
    } else {
      articulation_points(g, a);
    }
  } else {
    biconnected_components(g,
        make_iterator_property_map(ci, get(edge_index, g)));
  }

  return 0;
}

