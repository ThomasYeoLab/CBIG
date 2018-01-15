
#include <boost/property_map.hpp>

        template <typename Index, typename EdgeIndex>
        class simple_csr_edge {
        public:
            Index r;
            EdgeIndex i;
            simple_csr_edge(Index row, EdgeIndex ind) : r(row), i(ind) {}
            simple_csr_edge() : r(0), i(0) {}
            bool operator==(const simple_csr_edge& e) const {return i == e.i;}
            bool operator!=(const simple_csr_edge& e) const {return i != e.i;}
        }; // end simple_csr_edge
        
        // add an index map for the edge type
        template<typename EdgeIndex, typename Edge>
        struct simple_csr_edge_index_map {
          typedef EdgeIndex value_type;
          typedef EdgeIndex reference;
          typedef Edge key_type;
          typedef boost::readable_property_map_tag category;
        }; // end simple_csr_edge_index_map

    template<typename EdgeIndex, typename Edge>
    inline EdgeIndex
        get(const simple_csr_edge_index_map<EdgeIndex, Edge>&,
            const typename simple_csr_edge_index_map<EdgeIndex, Edge>::key_type& key)
    { return key.i; }

typedef simple_csr_edge<int,int> edge_t; 


template <class EdgeIndexMap, class Edge>
void test_edge_index_map_impl(const EdgeIndexMap& m, Edge& e)
{
    using namespace boost;
    
    int i = 5 + get(m,e);
}

template <class EdgeIndexMap>
void test_edge_index_map(EdgeIndexMap m) {
    using namespace boost;
    edge_t e;
    test_edge_index_map_impl(m,e);
    
    function_requires<ReadablePropertyMapConcept<EdgeIndexMap, edge_t> >();
}

template <class EdgeComponentPropertyMap>
void test_edge_property_map(EdgeComponentPropertyMap m)
{
    using namespace boost;

    //function_requires<WritablePropertyMapConcept<EdgeComponentPropertyMap, edge_t> >();

    edge_t e;
    int ci;
    put(m,e,ci);
}

int main(int argc, char **argv) 
{
    using namespace boost;

    int* ci;
    typedef simple_csr_edge_index_map<int,edge_t> edge_index_type;
    edge_index_type edge_index; 

    test_edge_index_map(edge_index);
    test_edge_property_map(make_iterator_property_map(ci, edge_index));

}

