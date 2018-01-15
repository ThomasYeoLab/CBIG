#ifndef MATLAB_BGL_VISITOR_MACROS_HPP
#define MATLAB_BGL_VISITOR_MACROS_HPP

/*
 * David Gleich
 * Copyright, Stanford University, 2006-2007
 */

/*
 * 8 July 2007
 * Switched to get(edge_index,g,e) instead of get(get(edge_index,g),e)
 */

#define VISITOR_VERTEX_FUNC(NAME,EXCEPT) \
    void NAME (typename boost::graph_traits<Graph>::vertex_descriptor v, const Graph& g) \
    { \
    if (vis->NAME != NULL) \
        { \
        if (vis->NAME (vis->pdata, (int)v) == 0) \
                throw EXCEPT (); \
        } \
    }

#define VISITOR_EDGE_FUNC(NAME,EXCEPT) \
    void NAME  (typename boost::graph_traits<Graph>::edge_descriptor e, const Graph& g) \
    { \
    if (vis->NAME != NULL) \
        { \
        int src=(int)boost::source(e,g), dst=(int)boost::target(e,g); \
        if (vis->NAME (vis->pdata, (int)get(boost::edge_index, g, e), src, dst) == 0) \
                throw EXCEPT (); \
        } \
    }


#endif // MATLAB_BGL_VISITOR_MACROS_HPP

