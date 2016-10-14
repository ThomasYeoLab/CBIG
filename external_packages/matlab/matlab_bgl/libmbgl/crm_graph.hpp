#ifndef YASMIC_COMPRESSED_ROW_MATRIX_GRAPH
#define YASMIC_COMPRESSED_ROW_MATRIX_GRAPH

#include <boost/iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/pending/integer_range.hpp>
#include <boost/graph/detail/adj_list_edge_iterator.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map.hpp>

#include <yasmic/compressed_row_matrix.hpp>
#include <yasmic/tuple_utility.hpp>


namespace boost
{
	struct yasmic_compressed_row_matrix_traversal_category : 
		public vertex_list_graph_tag,
		public adjacency_graph_tag,
		public incidence_graph_tag,
		public edge_list_graph_tag { };

	// forward declarations (prototypes) for the detailed implementation


	namespace impl 
	{
		//template <class EdgeList> struct val_out_edge_ret;
		template <class RowIter, class ColIter, class ValIter> class crm_graph_edge_iter;
		template <class RowIter, class ColIter, class ValIter> struct crm_graph_edge;

		template <class RowIter, class ColIter, class ValIter>
		struct crm_graph
		{
			typedef yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> type;
			typedef yasmic::smatrix_traits<type> traits;
		};

		/*template <class RowIter, class ColIter, class ValIter>
		struct nonzero_to_edge_transform
		{
			typedef
				typename impl::crm_graph_edge<RowIter, ColIter, ValIter>::type result_type;

			typedef 
				const typename crm_graph<RowIter, ColIter, ValIter>::traits::row_nonzero_descriptor argument_type;

			result_type operator() (argument_type arg) const
			{
				return (result_type(arg._row, arg._column, arg._nzi));
			}
		};*/

		template <class RowIter, class ColIter, class ValIter>
		struct nonzero_to_adjacency_transform
		{
			typedef
				typename crm_graph<RowIter, ColIter, ValIter>::traits::index_type result_type;
			typedef 
				const typename crm_graph<RowIter, ColIter, ValIter>::traits::row_nonzero_descriptor argument_type;

			result_type operator() (argument_type arg) const
			{
				return (arg._column);
			}
		};
	}


	template <class RowIter, class ColIter, class ValIter>
	struct graph_traits< yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> >
	{
		typedef yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> smatrix_type; 

		typedef typename yasmic::smatrix_traits<smatrix_type>::index_type vertex_descriptor;
		typedef typename yasmic::smatrix_traits<smatrix_type>::row_nonzero_descriptor edge_descriptor;
		
		typedef 
			boost::transform_iterator< 
				impl::nonzero_to_adjacency_transform<RowIter, ColIter, ValIter>,
				typename yasmic::smatrix_traits<smatrix_type>::row_nonzero_iterator > 
			adjacency_iterator;

		typedef typename yasmic::smatrix_traits<smatrix_type>::row_nonzero_iterator
			out_edge_iterator;
		
		typedef void in_edge_iterator;
    
		typedef typename yasmic::smatrix_traits<smatrix_type>::row_iterator vertex_iterator;
    
		typedef typename yasmic::smatrix_traits<smatrix_type>::nonzero_iterator edge_iterator;
		
		typedef directed_tag directed_category;
		typedef allow_parallel_edge_tag edge_parallel_category;
		typedef yasmic_compressed_row_matrix_traversal_category traversal_category;
		typedef typename yasmic::smatrix_traits<smatrix_type>::index_type vertices_size_type;
		typedef typename yasmic::smatrix_traits<smatrix_type>::size_type edges_size_type;
		typedef typename yasmic::smatrix_traits<smatrix_type>::size_type degree_size_type;

        static vertex_descriptor null_vertex()
        {
            return std::numeric_limits<vertex_descriptor>::max();
        }
	};

	namespace impl
	{
		

		template <class V, class S>
		class crm_graph_edge_type :
	    	public std::pair<V,V>
		{
		public:
    		crm_graph_edge_type() { }
    		crm_graph_edge_type(V s, V d, S id) 
    			: std::pair<V,V>(s, d), _id(id) { }
    	
    		S id() { return (_id); }
    	
		protected:
	    	S _id;
	    };

		template <class RowIter, class ColIter, class ValIter>
		struct crm_graph_edge
		{
			typedef typename yasmic::smatrix_traits<typename impl::crm_graph<RowIter, ColIter, ValIter>::type >::index_type V;
			typedef typename yasmic::smatrix_traits<typename impl::crm_graph<RowIter, ColIter, ValIter>::type >::size_type S;
	  
			typedef crm_graph_edge_type<V,S> type;
		};

		template <class RowIter, class ColIter, class ValIter>
		class crm_graph_edge_iter
		: public boost::iterator_facade<
            crm_graph_edge_iter<RowIter, ColIter, ValIter>,
			typename crm_graph_edge<RowIter, ColIter, ValIter>::type const,
            boost::forward_traversal_tag, 
            typename crm_graph_edge<RowIter, ColIter, ValIter>::type const >
		{
		public:
            crm_graph_edge_iter() {}

			crm_graph_edge_iter(
                RowIter ri, RowIter rend, ColIter ci)
            : _ri(ri), _rend(rend), _ci(ci), _id(0), _row(0)
            {
				_ri++;
            	// skip over any empty rows
            	while ( _ri != _rend && *(_ri) == _id )
                {
                   	// keep incrementing the row 
                   	++_ri; ++_row;
                }
            }

			crm_graph_edge_iter(
                RowIter ri, RowIter rend, ColIter ci, 
				typename std::iterator_traits<RowIter>::value_type id)
			: _ri(ri), _rend(rend), _ci(ci), _id(id), _row(0)
			{
			}

		private:
            friend class boost::iterator_core_access;

            void increment() 
            {  
            	// just increment everything!
            	++_ci; ++_id;
            	
                if (_id == *_ri)
                {
                	++_ri; ++_row;
                	
                	// while we aren't at the end and the row isn't empty
                	// (if *_ri == _id, then the row is empty because _id 
                	// is the current index in the column/val array.)
                	while ( _ri != _rend && *(_ri) == _id )
                    {
                    	// keep incrementing the row 
                    	++_ri; ++_row;
                    }
                }
            }

			bool equal(crm_graph_edge_iter const& other) const
            {
				return (_ri == other._ri && _ci == other._ci);
            }
            
            typename crm_graph_edge<RowIter, ColIter, ValIter>::type
            dereference() const 
            { 
            	//return boost::make_tuple(_row, *_ci, *_vi);
				return edge_type(_row, *_ci, _id);
            }


			typedef typename crm_graph_edge<RowIter, ColIter, ValIter>::type edge_type;
            typename std::iterator_traits<RowIter>::value_type _row;
            typename std::iterator_traits<RowIter>::value_type _id;
            RowIter _ri, _rend;
            ColIter _ci;
		};
		
		template <class RowIter, class ColIter, class ValIter>
		struct crm_graph_ret
		{
			typedef yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> smatrix_type; 
			typedef yasmic::smatrix_traits<smatrix_type> smatrix_traits;

			typedef boost::graph_traits<smatrix_type> g_traits;

			typedef std::pair< typename g_traits::vertex_iterator, typename g_traits::vertex_iterator > 
				vertices_ret;

			typedef std::pair< typename g_traits::adjacency_iterator, typename g_traits::adjacency_iterator >
				adjacent_ret;

			typedef std::pair < typename g_traits::out_edge_iterator, typename g_traits::out_edge_iterator >
				out_edges_ret;

			typedef std::pair < typename g_traits::edge_iterator, typename g_traits::edge_iterator >
				edges_ret;
		};
	}

	template <class RowIter, class ColIter, class ValIter>
	inline typename boost::graph_traits< yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> >::vertices_size_type
		num_vertices(const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g) 
	{
		return nrows(g);
	}  

	template <class RowIter, class ColIter, class ValIter>
	inline typename boost::graph_traits< yasmic::compressed_row_matrix<RowIter, ColIter, ValIter> >::vertices_size_type
		num_edges(const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		return nnz(g);
	}  

	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_descriptor
	source(
		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::edge_descriptor e,
		const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		return (e._row);
	}

	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_descriptor
	target(
		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::edge_descriptor e,
		const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		return (e._column);
	}


	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::out_edges_ret
		out_edges(typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_descriptor v,
			  const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		return (row_nonzeros(v,g));

		/*typedef typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::out_edge_iterator Iter;
		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::out_edges_ret return_type;

		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::smatrix_traits::row_nonzero_iterator rnzi, rnziend;

		tie(rnzi, rnziend) = yasmic::row_nonzeros(v, g);

		typedef typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::out_edge_iterator iter;

		return std::make_pair(iter(rnzi), iter(rnziend));*/


		/*typedef typename detail::val_out_edge_iter<EdgeList>::type Iter;
		typedef typename detail::val_out_edge_ret<EdgeList>::type return_type;
		return return_type(Iter(v, ++g[v].begin(), *(g[v].begin())),
    		Iter(v, g[v].end(), *(g[v].begin())));*/
	}

	template <class RowIter, class ColIter, class ValIter>
	typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::degree_size_type
	out_degree(
		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_descriptor v,
		const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		return (g._rstart[v+1] - g._rstart[v]);

		/*typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::smatrix_traits::row_nonzero_iterator rnzi, rnziend;

		tie(rnzi, rnziend) = yasmic::row_nonzeros(v, g);

		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::smatrix_traits::size_type rval = 0;
		while (rnzi != rnziend) { ++rnzi; ++rval; }

		return (rval);*/


	}

	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::adjacent_ret 
		adjacent_vertices(
			typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_descriptor v,
			const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::smatrix_traits::row_nonzero_iterator rnzi, rnziend;

		tie(rnzi, rnziend) = yasmic::row_nonzeros(v, g);

		typedef typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::adjacency_iterator iter;

		return std::make_pair(iter(rnzi), iter(rnziend));
	}

	// source() and target() already provided for pairs in graph_traits.hpp

	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::vertices_ret 
		vertices(const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		//typedef typename boost::integer_range<typename EdgeList::value_type>
		//::iterator Iter;
		typedef typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_iterator iter;
		return std::make_pair(iter(0), iter(nrows(g)));
	}

	template <class RowIter, class ColIter, class ValIter>
	inline typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::edges_ret 
		edges(const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
		//typedef typename boost::integer_range<typename EdgeList::value_type>
		//::iterator Iter;
		//typedef typename impl::crm_graph_ret<RowIter, ColIter, ValIter>::g_traits::vertex_iterator iter;
		//return std::make_pair(iter(0), iter(nrows(g)));

		return (nonzeros(g));
	}

	template <class RowIter, class ColIter, class ValIter>
	class compressed_row_matrix_graph_id_map
		: public put_get_helper<
			typename yasmic::smatrix_traits< typename impl::crm_graph<RowIter, ColIter, ValIter>::type >::size_type,
    		compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> >
	{
	protected:
  		typedef typename yasmic::smatrix_traits< typename impl::crm_graph<RowIter, ColIter, ValIter>::type >::size_type S;
	public:
  		typedef readable_property_map_tag category;
  		typedef S value_type;
  		typedef S reference;
		typedef typename graph_traits< typename impl::crm_graph<RowIter, ColIter, ValIter>::type >::edge_descriptor key_type;
	  	
  		compressed_row_matrix_graph_id_map() { }
  		template <class T> S operator [] (T x) const { return x._nzi; }
	};


	  
	template <class RowIter, class ColIter, class ValIter>
	inline identity_property_map get(vertex_index_t, 
		const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
  		identity_property_map pmap;
  		return (pmap);	
	}
	  
	template <class RowIter, class ColIter, class ValIter>
	inline compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> get(edge_index_t, 
  		const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
  		compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> pmap;
  		//return (vector_graph_id_map<EdgeList>());
  		return (pmap);
	}

	
	  
	template <class Tag>
	struct compressed_row_matrix_graph_property_map { };
	  
	template <>
	struct compressed_row_matrix_graph_property_map<vertex_index_t> {
  		template <class RowIter, class ColIter, class ValIter>
  		struct bind_ {
  			typedef identity_property_map type;
  			typedef identity_property_map const_type;
  		};
	};
	  
	template <>
	struct compressed_row_matrix_graph_property_map<edge_index_t> {
  		template <class RowIter, class ColIter, class ValIter>
  		struct bind_ {
  			typedef compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> type;
  			typedef compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> const_type;
  		};
	};

	template <>
	struct compressed_row_matrix_graph_property_map<edge_weight_t> {
		template <class RowIter, class ColIter, class ValIter>
		struct bind_ {
			typedef iterator_property_map<ValIter, compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> > type;
			typedef iterator_property_map<ValIter, compressed_row_matrix_graph_id_map<RowIter, ColIter, ValIter> > const_type;
		};
	};
	  
	template <class RowIter, class ColIter, class ValIter, class Tag>
	struct property_map<yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>, Tag> {
  		typedef typename compressed_row_matrix_graph_property_map<Tag>::template 
  			bind_<RowIter, ColIter, ValIter> map_gen;
  		typedef typename map_gen::type type;
  		typedef typename map_gen::const_type const_type;
	};
	
	template <class RowIter, class ColIter, class ValIter>
	inline typename property_map<yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>, edge_weight_t>::type
		get(edge_weight_t, const yasmic::compressed_row_matrix<RowIter, ColIter, ValIter>& g)
	{
  		return (make_iterator_property_map(g.begin_values(), get(edge_index, g)));
	}
}

#endif // YASMIC_COMPRESSED_ROW_MATRIX_GRAPH


