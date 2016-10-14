#include "include/matlab_bgl.h"
#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <boost/property_map.hpp>



template <class EdgeIndexMap, class Graph, class Edge>
void test_edge_index_map_impl(const EdgeIndexMap& m, Graph& g, Edge& e)
{
    using namespace boost;
    
    int i = 5 + get(m,e);
}

template <class EdgeIndexMap, class Graph>
void test_edge_index_map(EdgeIndexMap m, Graph& g) {
    using namespace boost;
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    edge_t e;
    test_edge_index_map_impl(m,g,e);
    function_requires<ReadablePropertyMapConcept<EdgeIndexMap, edge_t> >();
}

template <class EdgeComponentPropertyMap, class Graph>
void test_edge_property_map(EdgeComponentPropertyMap m, Graph& g)
{
    using namespace boost;

    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    //function_requires<WritablePropertyMapConcept<EdgeComponentPropertyMap, edge_t> >();

    edge_t e;
    mbglIndex ci;
    put(m,e,ci);
}

int main(int argc, char **argv) 
{
    using namespace yasmic;
    using namespace boost;

    mbglIndex* ci;


    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(0, 0, 0, NULL, NULL, NULL); 

    test_edge_index_map(get(edge_index,g),g);
    test_edge_property_map(make_iterator_property_map(ci, get(edge_index,g)), g);

    /*boost::tie(num_bicomps, oi) =
                biconnected_components(g, 
                    make_iterator_property_map(ci, get(edge_index,g)),
                    a);*/
}
