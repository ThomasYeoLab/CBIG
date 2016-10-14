#ifndef YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_AS_GRAPH_HPP
#define YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_AS_GRAPH_HPP

/**
 * @file simple_row_and_column_matrix_as_graph.hpp
 * Implement the types to 
 */

/*
 * David Gleich
 * Copyright, Stanford University, 2006-2007
 */

/*
 * 9 July 2007
 * Initial version
 */

#include <yasmic/simple_row_and_column_matrix.hpp>
#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/integer.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/mpl/if.hpp>
#include <boost/integer.hpp>
#include <boost/iterator/counting_iterator.hpp>

#define YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS \
    typename Index,typename Value,typename EdgeIndex
#define YASMIC_SIMPLE_RAC_GRAPH_TYPE \
    typename yasmic::simple_row_and_column_matrix<Index,Value,EdgeIndex>

namespace yasmic {
    namespace impl {
        struct simple_rac_graph_traversal : 
            public boost::vertex_list_graph_tag,
		    public boost::edge_list_graph_tag, 
            public boost::adjacency_graph_tag,
            public boost::bidirectional_graph_tag { };

        template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
        class simple_rac_in_edge_iterator;
    } // end namspase yasmic::impl
} // end namespace yasmic

template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
class yasmic::impl::simple_rac_in_edge_iterator
    : public boost::iterator_facade<
            typename yasmic::impl::simple_rac_in_edge_iterator<Index,Value,EdgeIndex>,
            yasmic::impl::simple_csr_edge<Index,EdgeIndex>,
            std::random_access_iterator_tag,
            typename yasmic::impl::simple_csr_edge<Index,EdgeIndex>,
            typename boost::int_t<CHAR_BIT * sizeof(EdgeIndex)>::fast>
{
private:
    typedef yasmic::impl::simple_csr_edge<Index,EdgeIndex> edge_descriptor;

public:
    typedef typename boost::int_t<CHAR_BIT * sizeof(EdgeIndex)>::fast 
        difference_type;

    // Implicit copy constructor OK
    explicit simple_rac_in_edge_iterator(Index *atj, EdgeIndex *atid)
        : atj(atj), atid(atid) 
    {}

    simple_rac_in_edge_iterator() {}

private:
    // iterator_facade requirements
    edge_descriptor dereference() const { return edge_descriptor(*atj,*atid); }

    bool equal(const simple_rac_in_edge_iterator<Index,Value,EdgeIndex>& other) const
    { return atid == other.atid; }

    void increment() { ++atid; ++atj; }
    void decrement() { ++atid; ++atj; }
    void advance(difference_type n) { atid+=n; atj+=n; } 

    difference_type distance_to(const simple_rac_in_edge_iterator<Index,Value,EdgeIndex>& other) const
    { return other.atid - atid; }

    Index *atj;
    EdgeIndex *atid;

    friend class boost::iterator_core_access;
};

namespace boost {
    // 
    // implement the graph traits
    //
    template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
    struct graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE> 
        : public graph_traits<yasmic::simple_csr_matrix<Index,Value,EdgeIndex> >
    {
        // modified requirement for Graph
        typedef yasmic::impl::simple_rac_graph_traversal traversal_category;
        // requirements for BidirectionalGraph
        typedef yasmic::impl::simple_rac_in_edge_iterator<Index,Value,EdgeIndex>
            in_edge_iterator;
    };
    //
    // implement the requirements for BidirectionalGraph
    //
    template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
    inline typename graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE>::degree_size_type
        in_degree(Index v, const YASMIC_SIMPLE_RAC_GRAPH_TYPE& g) { 
            return g.ati[v+1]-g.ati[v];
    }
    template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
    inline typename graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE>::degree_size_type
        degree(Index v, const YASMIC_SIMPLE_RAC_GRAPH_TYPE& g) { 
            return in_degree(v,g)+out_degree(v,g);
    }
    template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS>
    inline std::pair< typename graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE>::in_edge_iterator,
                      typename graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE>::in_edge_iterator >
        in_edges(Index v, const YASMIC_SIMPLE_RAC_GRAPH_TYPE& g) { 
            typedef typename graph_traits<YASMIC_SIMPLE_RAC_GRAPH_TYPE>::in_edge_iterator ei;
            EdgeIndex start = g.ati[v];
            EdgeIndex end = g.ati[v+1];
            return std::make_pair(ei(&g.atj[start],&g.atid[start]),
                                  ei(&g.atj[end],&g.atid[end]));
    }
    //
    // implement the property map, just inherit from csr matrix
    //
    template <YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS, typename Tag>
    struct property_map<YASMIC_SIMPLE_RAC_GRAPH_TYPE, Tag>
        : public property_map< yasmic::simple_csr_matrix<Index,Value,EdgeIndex>, Tag>
    {};

	
} // end namespace boost

#undef YASMIC_SIMPLE_RAC_GRAPH_TYPE
#undef YASMIC_SIMPLE_RAC_TEMPLATE_PARAMS

#endif // YASMIC_SIMPLE_ROW_AND_COLUMN_MATRIX_AS_GRAPH_HPP

