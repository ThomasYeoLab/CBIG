#ifndef YASMIC_BGL_KCORE
#define YASMIC_BGL_KCORE

/**
 * @file bgl_kcore.hpp
 * Implement the O(m) algorithm for computing the core number
 * of vertices in the graph.
 *
 * The algorithm comes from:
 *
 * Vladimir Batagelj and Matjaz Zaversnik, "An O(m) Algorithm for Cores 
 * Decomposition of Networks."  Sept. 1 2002.
 */

/*
 * David Gleich
 * Stanford University
 * 14 February 2006
 */

#if _MSC_VER >= 1400
	// disable the warning for ifstream::read
    #pragma warning( push )
	#pragma warning( disable : 4996 )
#endif // _MSC_VER >= 1400

#include <boost/graph/iteration_macros.hpp>
#include <boost/iterator/reverse_iterator.hpp>

namespace boost
{

	template <class Graph, class KCoreMap, class PositionMap>
	void core_numbers(const Graph& g, KCoreMap kcm, PositionMap pos)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex;
		typedef typename graph_traits<Graph>::degree_size_type size_type;

		BGL_FORALL_VERTICES_T(v,g,Graph)
		{
			kcm[v] = 0;
		}

		// compute the degree of all vertices
		BGL_FORALL_EDGES_T(e,g,Graph)
		{
			++kcm[source(e,g)];
		}

		size_type max_deg = 0;

		// compute the maximum degree
		BGL_FORALL_VERTICES_T(v,g,Graph)
		{
			if (kcm[v] > max_deg)
			{
				max_deg = kcm[v];
			}
		}

		// now we sort vertices into bins by their degree 
		// (we buffer this vector by 2 extra spots to make
		// some of the computations easier.
		//   1.  because deg > 0, we need max_deg+1 to index w/ deg itself
		//   2.  because we want to make things really easy, we extend
		//       the array one past degree to make computing partial_sums
		//       trivial
		std::vector<size_type> bin(max_deg+2);

		// compute the size of each bin
		BGL_FORALL_VERTICES_T(v,g,Graph)
		{
			++bin[kcm[v]];
		}

		// this loop sets bin[d] to the starting position of vertices
		// with degree d in the vert array
		size_type cur_pos = 0;
		for (size_type cur_deg = 0; cur_deg < max_deg+2; ++cur_deg)
		{
			size_type tmp = bin[cur_deg];
			bin[cur_deg] = cur_pos;
			cur_pos += tmp;
		}

		// place the vertices
		std::vector<vertex> vert(num_vertices(g));

		BGL_FORALL_VERTICES_T(v,g,Graph)
		{
			pos[v] = bin[kcm[v]];
			vert[pos[v]] = v;

			++bin[kcm[v]];
		}

		// we ``abused'' bin while placing the vertices, now, 
		// we need to restore it

		std::copy(boost::make_reverse_iterator(bin.end()-2),
			boost::make_reverse_iterator(bin.begin()+1), 
			boost::make_reverse_iterator(bin.end()-1));

		for (size_type i=0; i < num_vertices(g); ++i)
		{
			vertex v = vert[i];

			// we are now going to remove vertex v from the graph,
			// but only implicitly.  That is, we will decrement
			// the degree of each of the neighbors of v, and
			// adjust the sorting of the arrays appropriately.

			BGL_FORALL_ADJ_T(v,u,g,Graph)
			{
				// if kcm[u] > kcm[v], then u is still in the graph,
				// if kvm[u] = kcm[v], then we'll remove u soon, and
				// it's core number is the same as v.
				// if kvm[u] < kcm[v], we've already removed u.
				if (kcm[u] > kcm[v])
				{
					size_type deg_u = kcm[u];
					size_type pos_u = pos[u];
					

					// w is the first vertex with the same degree as u
					// (this is the resort operation!)
					size_type pos_w = bin[deg_u];
					vertex w = vert[pos_w];

					if (u != w)
					{
						// swap u and w
						pos[u] = pos_w;
						pos[w] = pos_u;
						vert[pos_w] = u;
						vert[pos_u] = w;
					}

					// now, the vertices array is sorted assuming
					// we perform the following step

					// start the set of vertices with degree of u 
					// one into the future (this now points at vertex 
					// w which we swapped with u).
					++bin[deg_u];

					// we are removing v from the graph, so u's degree
					// decreases
					--kcm[u];

				}

			}
		}
	}

}

#if _MSC_VER >= 1400
	// disable the warning for ifstream::read
    #pragma warning( pop )
#endif // _MSC_VER >= 1400


#endif // YASMIC_BGL_KCORE


