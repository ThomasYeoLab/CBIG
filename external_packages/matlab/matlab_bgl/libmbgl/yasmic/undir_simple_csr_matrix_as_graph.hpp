#ifndef YASMIC_UNDIR_SIMPLE_CSR_MATRIX_AS_GRAPH_HPP
#define YASMIC_UNDIR_SIMPLE_CSR_MATRIX_AS_GRAPH_HPP

/** @file undir_simple_csr_matrix.hpp
 * @author David F. Gleich
 * @copyright Stanford University, 2008
 * An almost exact copy of simple_csr_matrix_as_graph.hpp but with
 */

/** History
 *  2008-09-26: Initial coding
 *  2008-10-06: Fixed bug with adjacency iterator
 *    (this bug did not impact any algorithms, as no algorithms called it)
 */

#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/undir_simple_csr_matrix.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/integer.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/integer.hpp>
#include <yasmic/boost_mod/integer_extra.hpp>
#include <boost/iterator/counting_iterator.hpp>

#define YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS \
    typename Index,typename Value,typename EdgeIndex
#define YASMIC_SIMPLE_CSR_GRAPH_TYPE \
    typename yasmic::undir_simple_csr_matrix<Index,Value,EdgeIndex>

namespace yasmic {
    namespace impl {
        template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
        class undir_simple_csr_edge_iterator;
        template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
        class undir_simple_csr_out_edge_iterator;
    } // end namspase yasmic::impl
} // end namespace yasmic

template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
class yasmic::impl::undir_simple_csr_edge_iterator
{
public:
    typedef std::forward_iterator_tag iterator_category;
    typedef yasmic::impl::simple_csr_edge<Index,EdgeIndex> value_type;

    typedef const value_type* pointer;

    typedef value_type reference;
    // required due to "bug" in InputIterator concept, it is unused
    typedef typename boost::int_t<CHAR_BIT * sizeof(EdgeIndex)>::fast difference_type;

    undir_simple_csr_edge_iterator() : ai(NULL), current_edge(), end_of_this_vertex(0) {}

    undir_simple_csr_edge_iterator(
                const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g,
                value_type current_edge,
                EdgeIndex end_of_this_vertex)
    : ai(g.ai), current_edge(current_edge),
      end_of_this_vertex(end_of_this_vertex) {}

    // From InputIterator
    reference operator*() const { return current_edge; }
    pointer operator->() const { return &current_edge; }

    bool operator==(const undir_simple_csr_edge_iterator<Index,Value,EdgeIndex>& o) const {
        return current_edge == o.current_edge;
    }
    bool operator!=(const undir_simple_csr_edge_iterator<Index,Value,EdgeIndex>& o) const {
        return current_edge != o.current_edge;
    }

    undir_simple_csr_edge_iterator& operator++() {
        ++current_edge.i;
        while (current_edge.i == end_of_this_vertex) {
            ++current_edge.r;
            end_of_this_vertex = ai[current_edge.r + 1];
        }
        return *this;
    }

    undir_simple_csr_edge_iterator operator++(int) {
    	undir_simple_csr_edge_iterator temp = *this;
        ++*this;
        return temp;
    }
private:
    const EdgeIndex* ai;
    value_type current_edge;
    EdgeIndex end_of_this_vertex;
};

template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
class yasmic::impl::undir_simple_csr_out_edge_iterator
    : public boost::iterator_facade<
            typename yasmic::impl::undir_simple_csr_out_edge_iterator<Index,Value,EdgeIndex>,
            yasmic::impl::simple_csr_edge<Index,EdgeIndex>,
            std::random_access_iterator_tag,
            const typename yasmic::impl::simple_csr_edge<Index,EdgeIndex>&,
            typename boost::int_t<CHAR_BIT * sizeof(EdgeIndex)>::fast>
{
private:
    typedef yasmic::impl::simple_csr_edge<Index,EdgeIndex> edge_descriptor;

public:
    typedef typename boost::int_t<CHAR_BIT * sizeof(EdgeIndex)>::fast
        difference_type;

    undir_simple_csr_out_edge_iterator() {}

    // Implicit copy constructor OK
    explicit undir_simple_csr_out_edge_iterator(edge_descriptor e) : _e(e) { }

private:
    // iterator_facade requirements
    const edge_descriptor& dereference() const { return _e; }

    bool equal(const simple_csr_out_edge_iterator<Index,Value,EdgeIndex>& other) const
    { return _e == other._e; }

    void increment() { ++_e.i; }
    void decrement() { ++_e.i; }
    void advance(difference_type n) { _e.i += n; }

    difference_type distance_to(const yasmic::impl::undir_simple_csr_out_edge_iterator<Index,Value,EdgeIndex>& other) const
    { return other._e.i - _e.idx; }

    edge_descriptor _e;

    friend class boost::iterator_core_access;
};

namespace boost {
    //
    // implement the graph traits
    //
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    struct graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE> {
        // requirements for Graph
        typedef Index vertex_descriptor;
        typedef yasmic::impl::simple_csr_edge<Index,EdgeIndex> edge_descriptor;
        typedef undirected_tag directed_category;
        typedef allow_parallel_edge_tag edge_parallel_category;
        typedef yasmic::impl::simple_csr_graph_traversal traversal_category;
        static vertex_descriptor null_vertex()
        {
            return std::numeric_limits<vertex_descriptor>::max BOOST_PREVENT_MACRO_SUBSTITUTION ();
        }
        // requirements for VertexListGraph
        typedef typename yasmic::impl::remove_signedness<Index>::type vertices_size_type;
        typedef counting_iterator<Index> vertex_iterator;
        // requirements for EdgeListGraph
        typedef typename yasmic::impl::remove_signedness<EdgeIndex>::type edges_size_type;
        typedef yasmic::impl::undir_simple_csr_edge_iterator<Index,Value,EdgeIndex>
            edge_iterator;
        // requirements for IncidenceGraph
        typedef edges_size_type degree_size_type;
        typedef yasmic::impl::undir_simple_csr_out_edge_iterator<Index,Value,EdgeIndex>
            out_edge_iterator;
        // requirements for AdjacencyGraph
        typedef Index* adjacency_iterator;
        // requirements for various bugs
        typedef void in_edge_iterator;
    };

    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    struct graph_traits<const YASMIC_SIMPLE_CSR_GRAPH_TYPE>
    : graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE> {};

    //
    // implement the requirements for VertexListGraph
    //
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::vertices_size_type
        num_vertices(const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            return g.nrows;
    }
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline std::pair<counting_iterator<Index>,counting_iterator<Index> >
        vertices(const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            return std::make_pair(counting_iterator<Index>(0),
                                  counting_iterator<Index>(num_vertices(g)));
    }
    //
    // implement the requirements for EdgeListGraph
    //
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline Index source(
        typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e,
        const YASMIC_SIMPLE_CSR_GRAPH_TYPE&)
    {
        return e.r;
    }
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline Index target(
        typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e,
        const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g)
    {
        return g.aj[e.i];
    }
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edges_size_type
        num_edges(const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            return g.nnz;
    }

    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline std::pair< typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_iterator,
                      typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_iterator >
        edges(const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            typedef typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_iterator ei;
            typedef typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e;
            if (g.nnz == 0) {
                return std::make_pair(ei(),ei());
            }
            else {
                Index r=0;
                while (g.ai[r] == g.ai[r+1]) { ++r; }
                return std::make_pair(ei(g,e(r,0),g.ai[r+1]),
                                       ei(g,e(g.nrows,g.nnz),0));
            }
    }
    //
    // implement the requirements for IncidenceGraph
    //
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::degree_size_type
        out_degree(Index u, const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            return g.ai[u+1] - g.ai[u];
    }
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline std::pair< typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::out_edge_iterator,
                      typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::out_edge_iterator >
        out_edges(Index v, const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            typedef typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::out_edge_iterator ei;
            typedef typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e;

            return std::make_pair(ei(e(v,g.ai[v])),ei(e(v+1,g.ai[v+1])));
    }
    //
    // implement the requirements for adjacency_iterator
    //
    template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline std::pair< typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::adjacency_iterator,
                      typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::adjacency_iterator >
        adjacent_vertices(Index v, const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g) {
            return std::make_pair(&g.aj[g.ai[v]],&g.aj[g.ai[v+1]]);
    }

	template <YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS, typename Tag>
    struct property_map<YASMIC_SIMPLE_CSR_GRAPH_TYPE, Tag> {
    private:
        typedef identity_property_map vertex_index_type;
        typedef typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor
            edge_descriptor;
        typedef detail::simple_csr_edge_index_map<EdgeIndex,edge_descriptor> edge_index_type;
        typedef iterator_property_map<Value*,edge_index_type> edge_weight_type;

        typedef typename mpl::if_<is_same<Tag, edge_weight_t>,
                            edge_weight_type,
                            detail::error_property_not_found>::type
            edge_weight_or_none;

        typedef typename mpl::if_<is_same<Tag, edge_index_t>,
                            edge_index_type,
                            edge_weight_or_none>::type
            edge_prop_or_none;
    public:
	    typedef typename mpl::if_<is_same<Tag, vertex_index_t>,
                                vertex_index_type,
                                edge_prop_or_none>::type type;
        typedef type const_type;
	}; // end property_map

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline identity_property_map
    get(vertex_index_t, const YASMIC_SIMPLE_CSR_GRAPH_TYPE&)
    {
        return identity_property_map();
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline Index
    get(vertex_index_t,
        const YASMIC_SIMPLE_CSR_GRAPH_TYPE&, Index v)
    {
        return v;
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline typename property_map<YASMIC_SIMPLE_CSR_GRAPH_TYPE, edge_index_t>::const_type
    get(edge_index_t, const YASMIC_SIMPLE_CSR_GRAPH_TYPE&)
    {
        typedef typename property_map<YASMIC_SIMPLE_CSR_GRAPH_TYPE, edge_index_t>::const_type
            result_type;
        return result_type();
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline EdgeIndex
    get(edge_index_t, const YASMIC_SIMPLE_CSR_GRAPH_TYPE&,
        typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e)
    {
        return e.i;
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline typename property_map<YASMIC_SIMPLE_CSR_GRAPH_TYPE, edge_weight_t>::const_type
    get(edge_weight_t, const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g)
    {
        return make_iterator_property_map(g.a,get(edge_index,g));
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    inline Value
    get(edge_weight_t, const YASMIC_SIMPLE_CSR_GRAPH_TYPE& g,
        typename graph_traits<YASMIC_SIMPLE_CSR_GRAPH_TYPE>::edge_descriptor e)
    {
        return g.a[e.i];
    }

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    struct edge_property_type< YASMIC_SIMPLE_CSR_GRAPH_TYPE >  {
        typedef void type;
    };

    template<YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS>
    struct vertex_property_type< YASMIC_SIMPLE_CSR_GRAPH_TYPE >  {
        typedef void type;
    };

} // end namespace boost

#undef YASMIC_SIMPLE_CSR_GRAPH_TYPE
#undef YASMIC_SIMPLE_CSR_TEMPLATE_PARAMS

#endif /* YASMIC_UNDIR_SIMPLE_CSR_MATRIX_AS_GRAPH_HPP */

