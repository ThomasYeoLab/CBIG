#ifndef MATLAB_BGL_STOP_VISITORS_HPP
#define MATLAB_BGL_STOP_VISITORS_HPP

/*
 * David Gleich
 * Copyright, Stanford University, 2007
 * 18 April 2007
 */

/**
 * @file stop_visitors.hpp
 * A few custom visitors that will stop a search (bfs,dfs,dijkstra,etc.)
 */

template <class Vertex, class Exception, class Tag>
struct vertex_search_stopper
    : public boost::base_visitor<vertex_search_stopper<Exception,Vertex,Tag> >
{
    typedef Tag event_filter;
    vertex_search_stopper(Vertex v)
        : m_v(v) { }
    Vertex m_v;
    template <class Graph>
    void operator()(Vertex u, const Graph&) {
        if (m_v == u) { Exception e; throw e; }
    }
};

template <class Vertex, class Exception, class Tag>
vertex_search_stopper<Vertex, Exception, Tag>
stop_search_on_vertex_target(Vertex v, Exception, Tag) {
    return vertex_search_stopper<Vertex, Exception, Tag>(v);
}

#endif /* MATLAB_BGL_STOP_VISITORS_HPP */

