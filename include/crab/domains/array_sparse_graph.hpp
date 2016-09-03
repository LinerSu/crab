/*******************************************************************************
 * An array content domain.
 * 
 * This domain is a simplified implementation of the paper: 
 * "A Partial-Order Approach to Array Content Analysis" by 
 *  Gange, Navas, Schachte, Sondergaard, and Stuckey
 *  available here http://arxiv.org/pdf/1408.1754v1.pdf.
 *
 * It keeps a graph where vertices are array indexes and edges are
 * labelled with weights.  An edge (i,j) with weight w denotes that
 * the property w holds for the array segment [i,j). A weight is an
 * arbitrary numerical abstract domain that can relate multiple array
 * variables as well as array with scalar variables.
 ******************************************************************************/

#ifndef ARRAY_SPARSE_GRAPH_HPP
#define ARRAY_SPARSE_GRAPH_HPP

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>

#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>
#include <crab/domains/division_operators_api.hpp>
#include <crab/domains/domain_traits.hpp>

#include <crab/domains/array_sparse_graph/array_segmentation.hpp>
#include <crab/domains/array_sparse_graph/array_graph_ops.hpp>
#include <crab/domains/graphs/adapt_sgraph.hpp>
#include <crab/domains/graphs/sparse_graph.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/intervals.hpp>

// XXX: if expression domain is a template parameter no need to
// include this
#include <crab/domains/term_equiv.hpp>

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/join.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace crab {

  namespace domains {

    /* 
       A weighted directed graph where the weight is an abstract
       domain. The graph should be always kept in a consistent form,
       i.e., for all i,j,k:: weight(i,j) <= join(weight(i,k), weight(k,j))
    */
    template< typename Vertex, typename Weight, bool IsDistWeight >
    class array_sparse_graph_: public writeable {

     public:

      // XXX: make this a template parameter later
      //typedef AdaptGraph<Weight> graph_t;
      typedef SparseWtGraph<Weight> graph_t;
      typedef typename graph_t::vert_id _vert_id;
      typedef typename graph_t::edge_ref_t edge_ref_t;
      typedef typename graph_t::Wt Wt;
      typedef typename graph_t::mut_val_ref_t mut_val_ref_t;
      // XXX: needs to have this typedef so we can use GraphRev
      typedef Vertex vert_id;

     private:

      typedef GraphPerm<graph_t> GrPerm;
      typedef ArrayGrOps<graph_t, IsDistWeight> GrOps;
      typedef typename GrOps::Wt_join Wt_join;
      typedef typename GrOps::Wt_meet Wt_meet;

      typedef boost::container::flat_map<Vertex, _vert_id> vert_map_t;
      typedef typename vert_map_t::value_type vmap_elt_t;
      typedef vector<boost::optional<Vertex> > rev_map_t;
      typedef std::unordered_set<_vert_id> vert_set_t;
      typedef array_sparse_graph_<Vertex, Weight, IsDistWeight> array_sparse_graph_t;

      vert_map_t _vert_map;
      rev_map_t _rev_map;
      graph_t _g;
      vert_set_t _unstable;
      bool _is_bottom;

      struct vert_set_wrap_t {
        vert_set_wrap_t(const vert_set_t& _vs)
            : vs(_vs) { }
        
        bool operator[](_vert_id v) const {
          return vs.find(v) != vs.end();
        }
        const vert_set_t& vs;
      };

      _vert_id get_vert(Vertex v)
      {
        auto it = _vert_map.find(v);
        if(it != _vert_map.end())
          return (*it).second;

        _vert_id vert(_g.new_vertex());
        assert(vert <= _rev_map.size());

        if(vert < _rev_map.size()) {
          assert(!_rev_map[vert]);
          _rev_map[vert] = v;
        } else {
          _rev_map.push_back(v);
        }
        _vert_map.insert(vmap_elt_t(v, vert));

        return vert;
      }

     public: 

      template<class ItS>
      class iterator {
       public:
        typedef iterator<ItS> iter_t;
        iterator(const ItS& _it, const rev_map_t& _rev_map) 
            : it(_it), rev_map(_rev_map) { }
        bool operator!=(const iter_t& o) 
        { return it != o.it; }
        iter_t& operator++(void) { ++it; return *this; }
        Vertex operator*(void) const { 
          if (!rev_map[*it]) CRAB_ERROR("Reverse map failed");
          return *(rev_map[*it]);
        }
       protected:
        ItS it;      
        const rev_map_t& rev_map;
      };

      struct edge_t {
        edge_t(Vertex _v, Wt& _w) : vert(_v), val(_w) { }
        Vertex vert;
        Wt& val; 
      };

      template<class ItS>
      class edge_iterator {
       public:
        typedef edge_iterator<ItS> iter_t;
        edge_iterator(const ItS& _it, const rev_map_t& _rev_map) 
            : it(_it), rev_map(_rev_map) { }
        bool operator!=(const iter_t& o) 
        { return it != o.it; }
        iter_t& operator++(void) { ++it; return *this; }
        edge_t operator*(void) const { 
          edge_ref_t e = *it;
          if (!rev_map[e.vert]) CRAB_ERROR("Reverse map failed");
          Vertex v = *(rev_map[e.vert]);
          return edge_t(v, e.val);
        }
       protected:
        ItS it;      
        const rev_map_t& rev_map;
      };

      template<class Range, class ItS>
      class iterator_range {
       public:
        typedef ItS iterator;
        iterator_range(const Range &r, const rev_map_t &rev_map) 
            : _r(r), _rev_map(rev_map) { }
        iterator begin(void) const { return iterator(_r.begin(), _rev_map); }
        iterator end(void) const { return iterator(_r.end(), _rev_map); }
       protected:
        Range _r;
        const rev_map_t &_rev_map;
      };

      typedef iterator<typename graph_t::succ_iterator> succ_transform_iterator;
      typedef iterator<typename graph_t::pred_iterator> pred_transform_iterator;
      typedef iterator<typename graph_t::vert_iterator> vert_transform_iterator;
      typedef edge_iterator<typename graph_t::fwd_edge_iterator> fwd_edge_transform_iterator;
      typedef edge_iterator<typename graph_t::rev_edge_iterator> rev_edge_transform_iterator;

      typedef iterator_range<typename graph_t::succ_range, succ_transform_iterator> succ_range;
      typedef iterator_range<typename graph_t::pred_range, pred_transform_iterator> pred_range;
      typedef iterator_range<typename graph_t::vert_range, vert_transform_iterator> vert_range;
      typedef iterator_range<typename graph_t::e_succ_range, fwd_edge_transform_iterator> e_succ_range;
      typedef iterator_range<typename graph_t::e_pred_range, rev_edge_transform_iterator> e_pred_range;

      vert_range verts() {
        typename graph_t::vert_range p = _g.verts();
        return vert_range(p, _rev_map);
      }
      
      succ_range succs(Vertex v) {
        typename graph_t::succ_range p = _g.succs(get_vert(v));
        return succ_range(p, _rev_map);        
      }

      pred_range preds(Vertex v) {
        typename graph_t::pred_range p = _g.preds(get_vert(v));
        return pred_range(p, _rev_map);        
      }

      e_succ_range e_succs(Vertex v) {
        typename graph_t::e_succ_range p = _g.e_succs(get_vert(v));
        return e_succ_range(p, _rev_map);        
      }

      e_pred_range e_preds(Vertex v) {
        typename graph_t::e_pred_range p = _g.e_preds(get_vert(v));
        return e_pred_range(p, _rev_map);        
      }
      
     public:

      array_sparse_graph_(bool is_bottom = false)
          : writeable(), _is_bottom(is_bottom) 
      { }

      array_sparse_graph_(const array_sparse_graph_t& o)
          : writeable(),
            _vert_map(o._vert_map), _rev_map(o._rev_map), _g(o._g),
            _unstable(o._unstable), _is_bottom (false) 
      { 
        if (o._is_bottom)
          set_to_bottom();
      }

      array_sparse_graph_(array_sparse_graph_t&& o)
          : _vert_map(std::move(o._vert_map)), _rev_map(std::move(o._rev_map)),
            _g(std::move(o._g)), _unstable(std::move(o._unstable)), _is_bottom(o._is_bottom) 
      { }

      array_sparse_graph_(vert_map_t& vert_map, rev_map_t& rev_map, graph_t& g, vert_set_t unstable)
          : writeable(),
            _vert_map(vert_map), _rev_map(rev_map), _g(g), 
            _unstable(unstable), _is_bottom(false)
      { }
      
      array_sparse_graph_(vert_map_t&& vert_map, rev_map_t&& rev_map, graph_t&& g, vert_set_t &&unstable)
          : writeable(),
            _vert_map(std::move(vert_map)), _rev_map(std::move(rev_map)), _g(std::move(g)),
            _unstable(std::move(unstable)), _is_bottom(false)
      { }

      array_sparse_graph_t& operator=(const array_sparse_graph_t& o)
      {
        if(this != &o)
        {
          if(o._is_bottom)
            set_to_bottom();
          else {
            _is_bottom = false;
            _vert_map = o._vert_map;
            _rev_map = o._rev_map;
            _g = o._g;
            _unstable = o._unstable;
          }
        }
        return *this;
      }

      array_sparse_graph_t& operator=(array_sparse_graph_t&& o)
      {
        if(o._is_bottom) {
          set_to_bottom();
        } else {
          _is_bottom = false;
          _vert_map = std::move(o._vert_map);
          _rev_map = std::move(o._rev_map);
          _unstable = std::move(o._unstable);
          _g = std::move(o._g);
        }
        return *this;
      }

     public: 

      void set_to_bottom() {
        _vert_map.clear();
        _rev_map.clear();
        _g.clear();
        _unstable.clear();
        _is_bottom = true;
      }

      static array_sparse_graph_t top() { return array_sparse_graph_t(false); }
    
      static array_sparse_graph_t bottom() { return array_sparse_graph_t(true); }
    
      bool is_bottom() const { return _is_bottom; }
    
      bool is_top() {
        if(_is_bottom) 
          return false;
        return _g.is_empty();
      }

      bool lookup_edge(Vertex s, Vertex d, mut_val_ref_t* w) {
        if (is_bottom()) return false;
        auto se = get_vert(s);
        auto de = get_vert(d);
        return _g.lookup(se, de, w);
      }

      void update_edge (Vertex s, Weight w, Vertex d) {
        if (w.is_top())
          return;

        normalize();
        
        if (is_bottom ())
          return;
        
        auto se = get_vert(s);
        auto de = get_vert(d);
        Wt_meet op;
        _g.update_edge(se, w, de, op);
        GrOps::close_after_edge(_g, se, de);
      }

      void expand (Vertex s, Vertex d) {
        if(is_bottom()) 
          return;

        auto it = _vert_map.find(d);
        if(it != _vert_map.end()) {
          CRAB_ERROR("array_sparse_graph expand failed because vertex ", d, " already exists");
        }

        auto se = get_vert(s);        
        auto de = get_vert(d);
        
        for (auto edge : _g.e_preds(se))  
          _g.add_edge (edge.vert, edge.val, de);
        
        for (auto edge : _g.e_succs(se))  
          _g.add_edge (de, edge.val, edge.vert);

      }

      void normalize() {
        #if 0
        GrOps::closure(_g); // only for debugging purposes
        #else
        // Always maintained in closed form except for widening
        if(_unstable.size() == 0)
          return;
        GrOps::close_after_widen(_g, vert_set_wrap_t(_unstable));
        _unstable.clear();
        #endif 
      }

      void operator|=(array_sparse_graph_t& o) {
        *this = *this | o;
      }

      bool operator<=(array_sparse_graph_t& o)  {
        if (is_bottom()) 
          return true;
        else if(o.is_bottom())
          return false;
        else if (o.is_top ())
          return true;
        else if (is_top ())
          return false;
        else {
          normalize();

          if(_vert_map.size() < o._vert_map.size())
            return false;

          // Set up a mapping from o to this.
          vector<unsigned int> vert_renaming(o._g.size(),-1);
          for(auto p : o._vert_map)
          {
            auto it = _vert_map.find(p.first);
            // We can't have this <= o if we're missing some
            // vertex.
            if(it == _vert_map.end())
              return false;
            vert_renaming[p.second] = (*it).second;
          }

          assert(_g.size() > 0);
          mut_val_ref_t wx;

          for(_vert_id ox : o._g.verts()) {
            assert(vert_renaming[ox] != -1);
            _vert_id x = vert_renaming[ox];
            for(auto edge : o._g.e_succs(ox)) {
              _vert_id oy = edge.vert;
              assert(vert_renaming[ox] != -1);
              _vert_id y = vert_renaming[oy];
              auto ow = (Weight) edge.val;
              if(!_g.lookup(x, y, &wx) || (! ((Weight) wx <= ow))) 
                return false;
            }
          }
          return true;
        }
      }

      array_sparse_graph_t operator|(array_sparse_graph_t& o) {

        if (is_bottom() || o.is_top ())
          return o;
        else if (is_top () || o.is_bottom())
          return *this;
        else {
          CRAB_LOG ("array-sgraph",
                    crab::outs() << "Before join:\n"<<"Graph 1\n"<<*this
                                 <<"\n"<<"Graph 2\n"<<o << "\n");

          normalize();
          o.normalize();

          // Figure out the common renaming.
          vector<_vert_id> perm_x;
          vector<_vert_id> perm_y;

          vert_map_t out_vmap;
          rev_map_t out_revmap;

          for(auto p : _vert_map)
          {
            auto it = o._vert_map.find(p.first); 
            // Vertex exists in both
            if(it != o._vert_map.end())
            {
              out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
              out_revmap.push_back(p.first);

              perm_x.push_back(p.second);
              perm_y.push_back((*it).second);
            }
          }

          // Build the permuted view of x and y.
          assert(_g.size() > 0);
          GrPerm gx(perm_x, _g);
          assert(o._g.size() > 0);
          GrPerm gy(perm_y, o._g);

          // We now have the relevant set of relations. Because g_rx
          // and g_ry are closed, the result is also closed.
          graph_t join_g(GrOps::join(gx, gy));

          // Now garbage collect any unused vertices
          for(_vert_id v : join_g.verts())
          {
            if(join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0)
            {
              join_g.forget(v);
              if(out_revmap[v])
              {
                out_vmap.erase(*(out_revmap[v]));
                out_revmap[v] = boost::none;
              }
            }
          }

          array_sparse_graph_t res(std::move(out_vmap), std::move(out_revmap), 
                                   std::move(join_g), vert_set_t());
          CRAB_LOG ("array-sgraph", crab::outs() << "Result join:\n"<< res <<"\n";);
          return res;
        }
      }

      template<typename Thresholds>
      array_sparse_graph_t widening_thresholds (array_sparse_graph_t &o, 
                                                const Thresholds & /*ts*/) {
        return (*this || o);
      }
      
      array_sparse_graph_t operator||(array_sparse_graph_t &o) {	
        if (is_bottom())
          return o;
        else if (o.is_bottom())
          return *this;
        else {
          CRAB_LOG ("array-sgraph",
                    crab::outs() << "Before widening:\n"<<"Graph 1\n"<<*this
                                 <<"\n"<<"Graph 2\n"<<o<<"\n";);
          o.normalize();
          
          // Figure out the common renaming
          vector<_vert_id> perm_x;
          vector<_vert_id> perm_y;
          vert_map_t out_vmap;
          rev_map_t out_revmap;
          vert_set_t widen_unstable(_unstable);

          for(auto p : _vert_map)
          {
            auto it = o._vert_map.find(p.first); 
            // Vertex exists in both
            if(it != o._vert_map.end())
            {
              out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
              out_revmap.push_back(p.first);

              perm_x.push_back(p.second);
              perm_y.push_back((*it).second);
            }
          }
          
          // Build the permuted view of x and y.
          //assert(_g.size() > 0);
          GrPerm gx(perm_x, _g);            
          //assert(o._g.size() > 0);
          GrPerm gy(perm_y, o._g);
          
          // Now perform the widening 
          vector<_vert_id> destabilized;
          graph_t widen_g(GrOps::widen(gx, gy, destabilized));
          for(_vert_id v : destabilized) {
            widen_unstable.insert(v);
          }
          
          array_sparse_graph_t res(std::move(out_vmap), std::move(out_revmap), 
                                   std::move(widen_g), std::move(widen_unstable));
          CRAB_LOG ("array-sgraph", crab::outs() << "Result widening:\n"<<res<<"\n";);
          return res;
        }
      }

      array_sparse_graph_t operator&(array_sparse_graph_t& o) {

        if (is_bottom() || o.is_bottom())
          return bottom();
        else if (is_top())
          return o;
        else if (o.is_top())
          return *this;
        else {
          CRAB_LOG ("array-sgraph",
                    crab::outs() << "Before meet:\n"<<"Graph 1\n"<<*this<<"\n"
                                 <<"Graph 2\n"<<o << "\n");

          normalize();
          o.normalize();

          // Figure out the common renaming.
          vector<_vert_id> perm_x;
          vector<_vert_id> perm_y;

          vert_map_t out_vmap;
          rev_map_t out_revmap;

          for(auto p : _vert_map)
          {
            _vert_id vv = perm_x.size();
            out_vmap.insert(vmap_elt_t(p.first, vv));
            out_revmap.push_back(p.first);
            
            perm_x.push_back(p.second);
            perm_y.push_back(-1);
          }


          // Add missing mappings from the right operand.
          for(auto p : o._vert_map)
          {
            auto it = out_vmap.find(p.first);
            if(it == out_vmap.end())
            {
              _vert_id vv = perm_y.size();
              out_revmap.push_back(p.first);

              perm_y.push_back(p.second);
              perm_x.push_back(-1);
              out_vmap.insert(vmap_elt_t(p.first, vv));
            } else {
              perm_y[(*it).second] = p.second;
            }
          }

          // Build the permuted view of x and y.
          //assert(_g.size() > 0);
          GrPerm gx(perm_x, _g);
          //assert(o._g.size() > 0);
          GrPerm gy(perm_y, o._g);

          // Compute the syntactic meet of the permuted graphs.
          vector<_vert_id> changes;
          graph_t meet_g(GrOps::meet(gx, gy, changes));
          vert_set_t unstable;
          for(_vert_id v : changes)
            unstable.insert(v);

          GrOps::close_after_meet(_g, vert_set_wrap_t(unstable));

          array_sparse_graph_t res(std::move(out_vmap), std::move(out_revmap), 
                                   std::move(meet_g), vert_set_t());
          CRAB_LOG ("array-sgraph", crab::outs() << "Result meet:\n"<< res <<"\n";);
          return res;
        }
      }

      array_sparse_graph_t operator&&(array_sparse_graph_t& o) {
        if (is_bottom() || o.is_bottom())
          return bottom();
        else if (is_top ())
          return o;
        else{
          CRAB_LOG ("array-sgraph",
                    crab::outs() << "Before narrowing:\n"<<"Graph 1\n"<<*this<<"\n"
                                 <<"Graph 2\n"<<o<<"\n";);

          // Narrowing as a no-op should be sound.
          normalize();
          array_sparse_graph_t res(*this);
          
          CRAB_LOG ("array-sgraph",
                    crab::outs() << "Result narrowing:\n" << res<<"\n";);
          return res;
        }
      }

      void operator-=(Vertex v) {
        if (is_bottom ())
          return;
        auto it = _vert_map.find (v);
        if (it != _vert_map.end ()) {
          normalize();
          _g.forget(it->second);
          _rev_map[it->second] = boost::none;
          _vert_map.erase(v);
        }
      }

      void remove_from_weights(typename Weight::varname_t v) {
        mut_val_ref_t w_pq;
        for(auto p : _g.verts()) 
          for(auto e : _g.e_succs(p)) {
            auto q = e.vert;
            if (_g.lookup(p, q, &w_pq)) { 
              Weight w = (Weight) w_pq;
              w -= v;
              Wt_meet op;
              _g.update_edge(p, w, q, op);
              GrOps::close_after_edge(_g, p, q);
            }
          }
      }

      void write(crab_os& o) {
        write(o, true);
      }

      void write(crab_os& o, bool print_bottom_edges) {
        
        normalize ();

        if(is_bottom()){
          o << "_|_";
          return;
        }
        else if (is_top()){
          o << "{}";
          return;
        }
        else
        {
          bool first = true;
          o << "{";
          for(_vert_id s : _g.verts())
          {
            if(!_rev_map[s]) continue;
              
            auto vs = *_rev_map[s];
            for(_vert_id d : _g.succs(s))
            {
              if(!_rev_map[d]) continue;

              auto w = _g.edge_val(s, d);
              if (!print_bottom_edges && w.is_bottom())
                continue; // do not print bottom edges
                
              auto vd = *_rev_map[d];

              if(first)
                first = false;
              else
                o << ", ";
              o << "[" << vs << "," << vd << ")=>" << w;
            }
          }
          o << "}";
        }
      }
    };

    // Wrapper which uses shared references with copy-on-write.
    template<class Vertex, class Weight, bool IsDistWeight>
    class array_sparse_graph : public writeable {
      public:

      typedef array_sparse_graph_<Vertex, Weight, IsDistWeight> array_sgraph_impl_t;
      typedef std::shared_ptr<array_sgraph_impl_t> array_sgraph_ref_t;
      typedef array_sparse_graph<Vertex, Weight, IsDistWeight> array_sgraph_t;

      typedef typename array_sgraph_impl_t::Wt Wt;
      typedef typename array_sgraph_impl_t::mut_val_ref_t mut_val_ref_t;
      typedef typename array_sgraph_impl_t::graph_t graph_t;
      typedef typename array_sgraph_impl_t::vert_id vert_id;
      typedef typename array_sgraph_impl_t::succ_range succ_range;
      typedef typename array_sgraph_impl_t::pred_range pred_range;
      typedef typename array_sgraph_impl_t::vert_range vert_range;
      typedef typename array_sgraph_impl_t::e_succ_range e_succ_range;
      typedef typename array_sgraph_impl_t::e_pred_range e_pred_range;

      array_sparse_graph(array_sgraph_ref_t _ref) : norm_ref(_ref) { }

      array_sparse_graph(array_sgraph_ref_t _base, array_sgraph_ref_t _norm) 
        : base_ref(_base), norm_ref(_norm) { }

      array_sgraph_t create(array_sgraph_impl_t&& t)
      {
        return std::make_shared<array_sgraph_impl_t>(std::move(t));
      }

      array_sgraph_t create_base(array_sgraph_impl_t&& t)
      {
        array_sgraph_ref_t base = std::make_shared<array_sgraph_impl_t>(t);
        array_sgraph_ref_t norm = std::make_shared<array_sgraph_impl_t>(std::move(t));  
        return array_sgraph_t(base, norm);
      }

      void lock(void)
      { // Allocate a fresh copy.
        if(!norm_ref.unique())
          norm_ref = std::make_shared<array_sgraph_impl_t>(*norm_ref);
        base_ref.reset();
      }

    public:

      static array_sgraph_t top() { return array_sparse_graph(false); }
    
      static array_sgraph_t bottom() { return array_sparse_graph(true); }

      array_sparse_graph(bool is_bottom = false)
        : norm_ref(std::make_shared<array_sgraph_impl_t>(is_bottom)) { }

      array_sparse_graph(const array_sgraph_t& o)
        : base_ref(o.base_ref), norm_ref(o.norm_ref)
      { }

      array_sgraph_t& operator=(const array_sgraph_t& o) {
        base_ref = o.base_ref;
        norm_ref = o.norm_ref;
        return *this;
      }

      array_sgraph_impl_t& base(void) {
        if(base_ref)
          return *base_ref;
        else
          return *norm_ref;
      }

      array_sgraph_impl_t& norm(void) { return *norm_ref; }
      const array_sgraph_impl_t& norm(void) const { return *norm_ref; }

      bool is_bottom() { return norm().is_bottom(); }

      bool is_top() { return norm().is_top(); }

      bool operator<=(array_sgraph_t& o) { return norm() <= o.norm(); }

      void operator|=(array_sgraph_t o) { lock(); norm() |= o.norm(); }

      array_sgraph_t operator|(array_sgraph_t o) { return create(norm() | o.norm()); }

      array_sgraph_t operator||(array_sgraph_t o) { return create_base(base() || o.norm()); }

      array_sgraph_t operator&(array_sgraph_t o) { return create(norm() & o.norm()); }

      array_sgraph_t operator&&(array_sgraph_t o) { return create(norm() && o.norm()); }

      template<typename Thresholds>
      array_sgraph_t widening_thresholds (array_sgraph_t o, const Thresholds &ts) {
        return create_base(base().template widening_thresholds<Thresholds>(o.norm(), ts));
      }

      void normalize() { lock(); norm().normalize(); }

      vert_range verts () { return norm().verts(); }

      succ_range succs (Vertex v) { return norm().succs(v); }

      pred_range preds (Vertex v) { return norm().preds(v); }

      e_succ_range e_succs (Vertex v) { return norm().e_succs(v); }

      e_pred_range e_preds (Vertex v) { return norm().e_preds(v); }

      void set_to_bottom() { lock(); norm().set_to_bottom(); }

      bool lookup_edge(Vertex s, Vertex d, mut_val_ref_t* w) 
      { lock(); return norm().lookup_edge(s,d,w); }

      void expand(Vertex s, Vertex d) 
      { lock(); norm().expand(s,d); }

      void update_edge (Vertex s, Weight w, Vertex d) { lock(); norm().update_edge(s,w,d); }

      void remove_from_weights(typename Weight::varname_t v)
      { lock(); norm().remove_from_weights(v);}

      void operator-=(Vertex v) { lock(); norm() -= v; }

      void write(crab_os& o) { norm().write(o); }
      void write(crab_os& o, bool print_bottom_edges) { norm().write(o, print_bottom_edges); }

    protected:  
      array_sgraph_ref_t base_ref;  
      array_sgraph_ref_t norm_ref;
    };


    // Another C++ datatype to wrap variables and numbers as graph
    // vertices.

    enum landmark_kind_t { LMC, LMV, LMVP};
    template<class V, class N>
    class landmark {
     protected:
      landmark_kind_t _kind;
      landmark(landmark_kind_t kind): _kind(kind) { }
     public:
      virtual ~landmark() { }
      landmark_kind_t kind() const { return _kind;} 
      virtual bool operator==(const landmark<V,N>& o) const = 0;
      virtual bool operator<(const landmark<V,N>& o) const = 0;
      virtual void write(crab_os&o) const = 0;
      virtual std::size_t hash () const = 0;
    };

    template<class V, class N>
    class landmark_cst: public landmark<V,N> {
      N _n;
      typedef landmark<V,N> landmark_t;
      typedef landmark_cst<V,N> landmark_cst_t;

     public:
      landmark_cst (N n): landmark_t(landmark_kind_t::LMC), _n(n) {}

      bool operator==(const landmark_t& o) const {
        if (this->_kind != o.kind()) return false;

        assert (o.kind () == landmark_kind_t::LMC);
        auto o_ptr =  static_cast<const landmark_cst_t*>(&o);
        return (_n == o_ptr->_n);
      }

      bool operator<(const landmark_t& o) const {
        if (this->_kind != o.kind()) return true;

        assert (o.kind () == landmark_kind_t::LMC);
        return (_n < static_cast<const landmark_cst_t*>(&o)->_n);
      }

      void write(crab_os&o) const { o << _n; }

      std::size_t hash() const { return hash_value(_n);}

      N get_cst () const { return _n;}
    };

    template<class V, class N>
    class landmark_var: public landmark<V,N> {
      V _v;
      typedef landmark<V,N> landmark_t;
      typedef landmark_var<V,N> landmark_var_t;

     public:
      landmark_var (V v) : landmark_t(landmark_kind_t::LMV), _v(v) {}

      bool operator==(const landmark_t& o) const {
        if (this->_kind != o.kind()) return false;

        assert (o.kind () == landmark_kind_t::LMV);
        auto o_ptr = static_cast<const landmark_var_t*>(&o);
        return (_v == o_ptr->_v);
      }

      bool operator<(const landmark_t& o) const {
        if (this->_kind == o.kind()) {
          assert (o.kind () == landmark_kind_t::LMV);
          auto o_ptr =  static_cast<const landmark_var_t*>(&o);
          return (_v < o_ptr->_v);
        } else if (o.kind () == LMC) {
          return false;
        } else if (o.kind () == LMVP) {
          return true;
        } else  
          CRAB_ERROR("unreachable!");
      }

      void write(crab_os&o) const { o << _v; }

      std::size_t hash() const { return hash_value(_v);}

      V get_var () const { return _v;}
    };

    template<class V, class N>
    class landmark_varprime: public landmark<V,N> {
      std::string _lm;
      V _v;
      typedef landmark<V,N> landmark_t;
      typedef landmark_varprime<V,N> landmark_var_prime_t;

     public:
      landmark_varprime (std::string lm, V v) 
          : landmark_t(landmark_kind_t::LMVP), _lm(lm), _v(v) {}

      bool operator==(const landmark_t& o) const {
        if (this->_kind != o.kind()) return false;
        
        assert (o.kind () == landmark_kind_t::LMVP);
        auto o_ptr =  static_cast<const landmark_var_prime_t*>(&o);
        return (_v == o_ptr->_v);
      }

      bool operator<(const landmark_t& o) const {
        if (this->_kind != o.kind()) return false;

        assert (o.kind () == landmark_kind_t::LMVP);
        auto o_ptr =  static_cast<const landmark_var_prime_t*>(&o);
        return (_v < o_ptr->_v);
      }

      void write(crab_os&o) const { o << _lm << "'"; }

      std::size_t hash() const { return hash_value(_v);}

      V get_var () const { return _v;}
    };


    // Wrapper for landmark
    template<class V, class N>
    class landmark_ref {
      typedef landmark<V,N> landmark_t;
      typedef landmark_ref<V,N> landmark_ref_t;

     public:

      boost::shared_ptr<landmark_t> _ref;

      landmark_ref(V v, std::string name=""): _ref(nullptr) {
        if (name=="") {
          _ref = boost::static_pointer_cast<landmark_t>
              (boost::make_shared<landmark_var<V,N> >(landmark_var<V,N>(v)));
        } else {
          _ref = boost::static_pointer_cast<landmark_t>
              (boost::make_shared<landmark_varprime<V,N> >(landmark_varprime<V,N>(name, v)));
        }
      }
      landmark_ref(N n)
          : _ref(boost::static_pointer_cast<landmark_t>
                 (boost::make_shared<landmark_cst<V,N> >(landmark_cst<V,N>(n)))) { }
      
      landmark_kind_t kind() const { return _ref->kind();} 

      bool operator==(const landmark_ref &o) const { return (*_ref == *(o._ref)); } 

      bool operator<(const landmark_ref &o) const { return (*_ref < *(o._ref)); } 

      void write(crab_os& o) const { _ref->write(o); } 

      std::size_t hash () const { return _ref->hash();}
    };
 
    // super unsafe!
    template<class V, class N>
    inline V get_var(const landmark_ref<V,N>& lm) {
      assert (lm.kind() == LMVP);
      return boost::static_pointer_cast<const landmark_varprime<V,N> >(lm._ref)->get_var();
    }

    // super unsafe!
    template<class V, class N>
    inline N get_cst(const landmark_ref<V,N>& lm) {
      assert (lm.kind() == LMC);
      return boost::static_pointer_cast<const landmark_cst<V,N> >(lm._ref)->get_cst();
    }

    template<class V, class N>
    inline crab_os& operator<<(crab_os& o, const landmark_ref<V,N> &lm) {
      lm.write(o);
      return o;
    }

    template<class V, class N>
    inline std::size_t hash_value(const landmark_ref<V,N> &lm) {
      return lm.hash();
    }

    /*
      Reduced product of a numerical domain with a weighted array
      graph.
    */
    template<typename NumDom, typename Weight, bool IsDistWeight = false>
    class array_sparse_graph_domain: 
        public writeable,
        public numerical_domain<typename NumDom::number_t, typename NumDom::varname_t>,
        public bitwise_operators<typename NumDom::number_t, typename NumDom::varname_t>, 
        public division_operators<typename NumDom::number_t, typename NumDom::varname_t>
    {
      
     public:
      typedef typename NumDom::number_t Number;
      typedef typename NumDom::varname_t VariableName;
      
      // WARNING: assume NumDom::number_t = Weight::number_t and
      //                 NumDom::varname_t = Weight::varname_t
      using typename numerical_domain< Number, VariableName>::linear_expression_t;
      using typename numerical_domain< Number, VariableName>::linear_constraint_t;
      using typename numerical_domain< Number, VariableName>::linear_constraint_system_t;
      using typename numerical_domain< Number, VariableName>::variable_t;
      using typename numerical_domain< Number, VariableName>::number_t;
      using typename numerical_domain< Number, VariableName>::varname_t;
      typedef interval<Number> interval_t;
      
      typedef landmark_cst<VariableName,Number> landmark_cst_t;
      typedef landmark_var<VariableName,Number> landmark_var_t;
      typedef landmark_varprime<VariableName,Number> landmark_var_prime_t;
      typedef landmark_ref<VariableName,Number> landmark_ref_t;
      typedef array_sparse_graph<landmark_ref_t,Weight,IsDistWeight> array_sgraph_t;
      typedef array_sparse_graph_domain<NumDom,Weight,IsDistWeight> array_sgraph_domain_t;

      //// XXX: make this a template parameter later
      typedef crab::cfg::var_factory_impl::str_var_alloc_col::varname_t str_varname_t;
      typedef interval_domain<z_number, str_varname_t> str_interval_dom_t;
      typedef term::TDomInfo<z_number, varname_t, str_interval_dom_t> idom_info;
      typedef term_domain<idom_info> expression_domain_t;  

     private:
      typedef typename array_sgraph_t::mut_val_ref_t mut_val_ref_t;

      // Quick wrapper to perform efficient unsat queries on the
      // scalar domain.
      struct solver_wrapper {
        // XXX: do not pass by reference
        NumDom _inv;
        solver_wrapper(NumDom inv): _inv(inv) { }
        bool is_unsat (linear_constraint_t cst) {
          // XXX: it might modify _inv so that's why we make a copy in
          // the constructor.
          return array_sgraph_domain_traits<NumDom>::is_unsat(_inv, cst);        
        }
      };

      NumDom _scalar;        
      expression_domain_t _expressions; // map each program variable to a symbolic expression
      array_sgraph_t _g;        

      // A landmark is either a variable or number that may appear as
      // an array index. In addition, for each landmark l we keep
      // track of a prime landmark l' whose meaning is l'=l+1.

      /// === Static data
      typedef boost::unordered_map<landmark_ref_t,landmark_ref_t> lm_map_t;
      static lm_map_t var_landmarks;
      static lm_map_t cst_landmarks;

      // --- landmark iterators
      struct get_first : public std::unary_function<typename lm_map_t::value_type, landmark_ref_t> {
        get_first () {}
        landmark_ref_t operator()(const typename lm_map_t::value_type &p) const 
        { return p.first; }
      }; 
      struct get_second : public std::unary_function<typename lm_map_t::value_type, landmark_ref_t> {
        get_second () {}
        landmark_ref_t operator()(const typename lm_map_t::value_type &p) const 
        { return p.second; }
      }; 
      typedef boost::transform_iterator<get_first, 
                                        typename lm_map_t::iterator> lm_iterator;
      typedef boost::transform_iterator<get_second, 
                                        typename lm_map_t::iterator> lm_prime_iterator;
      typedef boost::iterator_range<lm_iterator> lm_range;
      typedef boost::iterator_range<lm_prime_iterator> lm_prime_range;

      lm_prime_iterator var_lm_prime_begin()
      { return boost::make_transform_iterator(var_landmarks.begin(), get_second());}
      lm_prime_iterator var_lm_prime_end()
      { return boost::make_transform_iterator(var_landmarks.end(), get_second());}
      lm_prime_range var_lm_primes() 
      { return boost::make_iterator_range(var_lm_prime_begin(), var_lm_prime_end());}

      lm_prime_iterator cst_lm_prime_begin()
      { return boost::make_transform_iterator(cst_landmarks.begin(), get_second());}
      lm_prime_iterator cst_lm_prime_end()
      { return boost::make_transform_iterator(cst_landmarks.end(), get_second());}
      lm_prime_range cst_lm_primes() 
      { return boost::make_iterator_range(cst_lm_prime_begin(), cst_lm_prime_end());}
                                          
     
     public:

      template<class CFG, class VarFactory>
      static void do_initialization (CFG cfg, VarFactory& vfac) {

        typedef crab::analyzer::array_segmentation<CFG,typename CFG::varname_t> array_segment_analysis_t;
        typedef typename array_segment_analysis_t::array_segment_domain_t array_segment_domain_t;
        typedef crab::analyzer::array_constant_segment_visitor<array_segment_domain_t> array_cst_segment_visitor_t;

        std::set<landmark_ref_t> lms;

        // add variables 
        array_segment_analysis_t analysis(cfg);
        analysis.exec();
        auto var_indexes = analysis.get_variables(cfg.entry());
        lms.insert(var_indexes.begin(), var_indexes.end());

        // add constants
        // make sure 0 is always considered as an array index
        lms.insert(landmark_ref_t(number_t(0)));
        typename array_cst_segment_visitor_t::constant_set_t constants;
        for (auto &bb: boost::make_iterator_range(cfg.begin(), cfg.end())) {
          auto var_indexes = analysis.get_variables(bb.label());
          // XXX: use some heuristics to choose "relevant" constants
          array_cst_segment_visitor_t vis(var_indexes);          
          for (auto &s: boost::make_iterator_range(bb.begin(), bb.end()))
            s.accept(&vis);
          auto cst_indexes = vis.get_constants();
          lms.insert(cst_indexes.begin(), cst_indexes.end());
        }

        set_landmarks (lms, vfac);
      }

      template<class Range, class VarFactory>
      static void set_landmarks(const Range& lms, VarFactory& vfac) {
        
        var_landmarks.clear();
        cst_landmarks.clear();
        
        unsigned num_vl = 0;
        unsigned num_cl = 0;

        for (auto lm: lms) {
          switch (lm.kind()) {
            case LMV: {
              auto v = boost::static_pointer_cast<const landmark_var_t>(lm._ref)->get_var();
              varname_t v_prime = vfac.get(v.index());
              landmark_ref_t lm_prime(v_prime, v.str());
              var_landmarks.insert(std::make_pair(lm, lm_prime));
              num_vl++;
              break;
            }
            case LMC: {
              auto n = boost::static_pointer_cast<const landmark_cst_t>(lm._ref)->get_cst();
              varname_t v_prime = vfac.get(); 
              landmark_ref_t lm_prime(v_prime, n.get_str());
              cst_landmarks.insert(std::make_pair(lm, lm_prime));
              num_cl++;
              break;
            }
            default: 
              CRAB_ERROR("A landmark can only be either variable or constant");
          }
        }
        CRAB_LOG("array-sgraph-domain",
                 crab::outs() << "Added " << num_vl << " variable landmarks "
                              << "and " << num_cl << " constant landmarks={";
                 bool first=true;
                 for (auto &l: var_landmarks) {
                   if (!first) crab::outs() << ",";
                   first=false;
                   crab::outs() << l.first;
                 }
                  for (auto &l: cst_landmarks) {
                   if (!first) crab::outs() << ",";
                   first=false;
                   crab::outs() << l.first;
                 }
                 crab::outs() << "}\n";
                 );
      }

      // Some funtions like normalize_offset needs to add temporary
      // landmarks
      template<class VarFactory>
      static void add_landmark(landmark_ref_t lm, VarFactory& vfac)
      {
        if (lm.kind () != LMV) return;

        auto v = boost::static_pointer_cast<const landmark_var_t>(lm._ref)->get_var();
        varname_t v_prime = vfac.get(v.index());
        landmark_ref_t lm_prime(v_prime, v.str());
        var_landmarks.insert(std::make_pair(lm, lm_prime));
      }

     private:

      void set_to_bottom(){
        _scalar = NumDom::bottom();
        _expressions = expression_domain_t::bottom();
        _g.set_to_bottom();
      }

      linear_expression_t make_expr (landmark_ref_t x) {
        switch (x.kind()) {
          case LMC: 
            return boost::static_pointer_cast<landmark_cst_t>(x._ref)->get_cst();
          case LMV: 
            return variable_t(boost::static_pointer_cast<landmark_var_t>(x._ref)->get_var());
          case LMVP:
            return variable_t(boost::static_pointer_cast<landmark_var_prime_t>(x._ref)->get_var());
          default:
            CRAB_ERROR("unreachable!");
        }
      }

      // make constraint x < y
      linear_constraint_t make_lt_cst (landmark_ref_t x, landmark_ref_t y) {
        return linear_constraint_t(make_expr(x) <= make_expr(y) - 1);
      }

      // make constraint x <= y
      linear_constraint_t make_leq_cst (landmark_ref_t x, landmark_ref_t y) {
        return linear_constraint_t(make_expr(x) <= make_expr(y));
      }

      // make constraint x == y
      linear_constraint_t make_eq_cst (landmark_ref_t x, landmark_ref_t y) {
        return linear_constraint_t(make_expr(x) == make_expr(y));
      }

      // make constraint x' == x+1
      // FIXME: actually the increment should be equal to the array
      // element size.
      linear_constraint_t make_prime_relation(landmark_ref_t x_prime, landmark_ref_t x){
        return linear_constraint_t(make_expr(x_prime) == make_expr(x) + 1);
      }

      // Return the weight from the edge (i, i')
      Weight array_elem (VariableName i) {
        if (is_bottom()) return Weight::bottom();
        if (is_top()) return Weight::top();
          
        landmark_ref_t lm_i (i);       
        auto it = var_landmarks.find(lm_i);
        if (it == var_landmarks.end()) {
          return Weight::top();
        }

        mut_val_ref_t wi;   
        if (_g.lookup_edge(lm_i, it->second, &wi))
          return (Weight) wi;
        else
          return Weight::top();
      } 

      // Remove v from the edge (i,i')
      void array_forget(VariableName i, VariableName v) {
        if (is_bottom()) return;

        landmark_ref_t lm_i (i);
        auto it = var_landmarks.find(lm_i);
        if (it == var_landmarks.end ()) {
          return;
        }

        mut_val_ref_t wi;          
        if (_g.lookup_edge(lm_i, it->second, &wi)) {
          Weight w = (Weight) wi;
          w -= v;
          // XXX: update_edge closes the array graph
          _g.update_edge (lm_i, w, it->second);
        }
      }

      // Remove v from all vertices and edges
      void array_forget(VariableName v) {
        if (is_bottom()) return;
        
        if (var_landmarks.find(landmark_ref_t(v)) == var_landmarks.end()) {
          return;
        }

        _g -= landmark_ref_t(v);
        _g.remove_from_weights(v);
      }

      // Update the weight from the edge (i, i')
      void array_update (VariableName i, Weight w)
      {
        if (is_bottom()) return;
        
        landmark_ref_t lm_i(i);
        //--- strong update
        auto it = var_landmarks.find(lm_i);
        if (it == var_landmarks.end ())
          return;

        _g.update_edge(lm_i, w, it->second);
        mut_val_ref_t wi;          
        if (!_g.lookup_edge(lm_i, it->second, &wi))
          return; 

        //--- weak update
        // An edge (p,q) must be weakened if p <= i <= q and p < q
        solver_wrapper solve(_scalar);
        mut_val_ref_t w_pq;
        for(auto p : _g.verts ()) {
          for(auto e : _g.e_succs(p)) {
            auto q = e.vert;
            if ((p == lm_i) &&  (q == it->second)) 
              continue;
            if (_g.lookup_edge(p, q, &w_pq) && ((Weight) w_pq).is_bottom())
              continue;
            // we know already that p < q in the array graph

            // check p <= i  
            if (solve.is_unsat(make_leq_cst(p, lm_i)))
              continue;
            // check i' <= q
            if (solve.is_unsat(make_leq_cst(it->second, q)))
              continue;

            w_pq = (Weight) w_pq | (Weight) wi;
          }
        }
      }

      interval_t eval_interval(linear_expression_t e)
      {
        interval_t r = e.constant();
        for (auto p : e)
          r += p.first * _scalar[p.second.name()];
        return r;
      }

      // x := x op k 
      template<typename VarOrNum>
      void apply_landmark (operation_t op, VariableName x, VarOrNum k) { 
        if (is_bottom()) return;

        landmark_ref_t lm_x(x);
        auto it = var_landmarks.find(lm_x);
        if (it == var_landmarks.end()) {
          // If x is not a landmark we just apply the operation on the
          // scalar domain and return.
          apply_only_scalar(op, x, x, k);
          return;
        }
 
        /// --- Add x_old and x_old' to store old values of x and x'

        VariableName x_old = x.get_var_factory().get();      
        VariableName x_old_prime = x.get_var_factory().get(); 
        landmark_ref_t lm_x_old (x_old);
        landmark_ref_t lm_x_old_prime (x_old_prime, x_old.str());
        var_landmarks.insert(std::make_pair(lm_x_old, lm_x_old_prime));
        // x_old = x
        _scalar.assign(x_old, linear_expression_t(x)); 
        // relation between x_old and x' 
        _scalar += make_prime_relation(lm_x_old_prime, lm_x_old);
        //_scalar += make_eq_cst(lm_x_old_prime, it->second);      

        /*** Incremental graph reduction ***/
        //// x_old  has all the x predecessors and successors 
        _g.expand(lm_x, lm_x_old); 
        //// x_old' has all the x' predecessors and successors 
        _g.expand(it->second, lm_x_old_prime); 
        //// edges between x and x_old 
        _g.update_edge(lm_x, Weight::bottom(), lm_x_old);        
        _g.update_edge(lm_x_old, Weight::bottom(), lm_x);        
        //// edges between x' and x_old' 
        _g.update_edge(it->second, Weight::bottom(), lm_x_old_prime);        
        _g.update_edge(lm_x_old_prime, Weight::bottom(), it->second);        
        //// edges between x_old and x_old'
        mut_val_ref_t w;   
        if (_g.lookup_edge(lm_x, it->second, &w))
          _g.update_edge(lm_x_old, (Weight) w, lm_x_old_prime);        
        _g.update_edge(lm_x_old_prime, Weight::bottom(), lm_x_old);        

        /// --- Remove x and x'
        _g -= lm_x;
        _g -= it->second;

        /// --- Perform operation in the scalar domain
        _scalar.apply(op, x, x, k); 

        //restore relation between x and x'
        _scalar.apply(OP_ADDITION, get_var(it->second), x, 1);
        //_scalar -= get_var(it->second);
        //_scalar += make_prime_relation(it->second, lm_x);

        if (!reduce(_scalar, _g)) { // FIXME: incremental version
          set_to_bottom();
          return;
        }

        /// --- Remove x_old and x_old'
        _g -= lm_x_old;
        _g -= lm_x_old_prime;
        _scalar -= x_old;
        _scalar -= x_old_prime;
        var_landmarks.erase(lm_x_old);
      }

      void forget_prime_var(VariableName v) {
        // remove v' from scalar and array graph
        auto it = var_landmarks.find(landmark_ref_t(v));
        if (it != var_landmarks.end()) {
          _scalar -= get_var(it->second);
          // XXX: v' cannot appear in the array weights so we do not
          //      need to call array_forget.
          _g -= it->second;
        }
      }

      // perform the operation in the scalar domain assuming that
      // nothing can be done in the graph domain.
      template<class Op, class K>
      void apply_only_scalar(Op op, VariableName x, VariableName y, K k) {

        _scalar.apply(op, x, y, k);

        // Abstract x in the array graph
        if (var_landmarks.find(landmark_ref_t(x)) != var_landmarks.end()){ 
          array_forget(x);     // remove x from the array graph
          forget_prime_var(x); // remove x' from scalar and array graph

          /// XXX: I think no need to reduce here
        }
      }
            
      // By active we mean current variables that are kept track by
      // the scalar domain.
      void get_active_landmarks(NumDom &scalar, std::vector<landmark_ref_t> & landmarks) const {
        /// FIXME: avoid repeating this loop over and over
        landmarks.reserve(cst_landmarks.size());
        for (auto p: cst_landmarks) { 
          landmarks.push_back (p.first);
          landmarks.push_back (p.second);
        }
        auto active_vars = array_sgraph_domain_traits<NumDom>::active_variables(scalar);        
        for (auto v: active_vars) {
          auto it = var_landmarks.find(landmark_ref_t(v));
          if (it != var_landmarks.end()){
            landmarks.push_back(landmark_ref_t(v)); 
            landmarks.push_back(it->second);
          }
        }
      }

      VariableName normalize_offset (VariableName o, z_number n)
      {
        CRAB_LOG("array-sgraph-domain-norm",
                 crab::outs() << "EXPRESSIONS (before normalize offset)=" << _expressions << "\n");

        // --- create a fresh variable no such that no := o;
        VariableName no = o.get_var_factory().get();
        _expressions.assign (no, linear_expression_t (o));

        // -- apply no := no / n; in the expressions domain
        _expressions.apply (operation_t::OP_DIVISION, no, no, n);
        
        // -- simplifiy the expression domain 
        // FIXME: should be part of array_sgraph_domain_traits 
        bool simp_done = _expressions.simplify (no);

        CRAB_LOG("array-sgraph-domain-norm",
                 crab::outs() << "EXPRESSIONS (after normalize offset)=" << _expressions << "\n");

        if (!simp_done) return o;
          
        // -- propagate equalities from _expressions to _scalar
        product_domain_traits<expression_domain_t, NumDom>::push(no, _expressions, _scalar);

        // -- add to _scalar the constraint no' = no + 1
        landmark_ref_t lm_no (no);
        add_landmark (lm_no, o.get_var_factory());

        auto it = var_landmarks.find (lm_no);
        assert (it != var_landmarks.end());
        _scalar += make_prime_relation(it->second, lm_no);
        
        // reduce between _scalar and the array graph
        if (!reduce(_scalar, _g)) { // FIXME: incremental version
          // XXX: I think we should never get bottom here.
          set_to_bottom();
        }

        // cleanup of the expression abstraction
        _expressions -= no;
        
        return no;
      }


     public:

      // The reduction consists of detecting dead segments so it is
      // done only in one direction (scalar -> array graph).  
      // Return false if bottom is detected during the reduction.
      bool reduce(NumDom &scalar, array_sgraph_t &g) {
        crab::CrabStats::count (getDomainName() + ".count.reduce");
        crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

        domain_traits<NumDom>::normalize(scalar);
        g.normalize();

        if (scalar.is_bottom() || g.is_bottom())
          return false;

        if (!scalar.is_top ()) { 
          std::vector<landmark_ref_t> active_landmarks;
          get_active_landmarks(scalar, active_landmarks);
          solver_wrapper solve(scalar);
          for (auto lm_s : active_landmarks)
            for (auto lm_d : active_landmarks) {
              // XXX: we do not exploit the following facts:
              //   - i < i' is always sat
              //   - i' < i is always unsat
              //   - if i < j  unsat then i' < j unsat.
              //   - if i < j' unsat then i' < j unsat.
              if ((lm_s == lm_d) || solve.is_unsat (make_lt_cst(lm_s,lm_d))) 
                g.update_edge(lm_s, Weight::bottom(), lm_d);
            }
        }        
        return (!g.is_bottom());
      }

      
      static array_sgraph_domain_t top () {
        return array_sgraph_domain_t(false);
      }

      static array_sgraph_domain_t bottom () {
        return array_sgraph_domain_t(true);
      }

     public:

      array_sparse_graph_domain(bool is_bottom=false)
          : _scalar(NumDom::top()), _expressions(expression_domain_t::top ()), 
            _g(array_sgraph_t::top()) { 
        if (is_bottom) 
          set_to_bottom();
      }

      array_sparse_graph_domain(const NumDom& s, const expression_domain_t& e, 
                                const array_sgraph_t& g)
          : _scalar(s), _expressions(e), _g(g) { 
        if (_scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom())
          set_to_bottom();
      }
    
      array_sparse_graph_domain(NumDom &&s, expression_domain_t &&e, 
                                array_sgraph_t &&g)
          : _scalar(std::move(s)), _expressions (std::move(e)), _g(std::move(g)) { 
        if (_scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom())
          set_to_bottom();
      }

      array_sparse_graph_domain(const array_sgraph_domain_t&o)
          : _scalar(o._scalar), _expressions (o._expressions), _g(o._g) { 
        crab::CrabStats::count (getDomainName() + ".count.copy");
        crab::ScopedCrabStats __st__(getDomainName() + ".copy");
      }

      array_sparse_graph_domain(array_sgraph_domain_t &&o)
          : _scalar(std::move(o._scalar)), 
            _expressions (std::move(o._expressions)), 
            _g(std::move(o._g)) { 
      }

      array_sgraph_domain_t& operator=(const array_sgraph_domain_t& o) {
        crab::CrabStats::count (getDomainName() + ".count.copy");
        crab::ScopedCrabStats __st__(getDomainName() + ".copy");
        if(this != &o) {
          _scalar = o._scalar;
          _expressions = o._expressions;
          _g = o._g;
        }
        return *this;
      }

      array_sgraph_domain_t& operator=(array_sgraph_domain_t &&o) {
        _scalar = std::move(o._scalar);
        _expressions = std::move(o._expressions);
        _g = std::move(o._g);
        return *this;
      }
      
      bool is_top() {
        return _scalar.is_top () && _expressions.is_top () && _g.is_top();
      }
      
      bool is_bottom() {
        return _scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom();
      }

      bool operator<=(array_sgraph_domain_t &o) {
        crab::CrabStats::count (getDomainName() + ".count.leq");
        crab::ScopedCrabStats __st__(getDomainName() + ".leq");

        CRAB_LOG("array-sgraph-domain",
                 crab::outs () << "Leq " << *this << " and\n"  << o << "=\n";);
        bool res = (_scalar <= o._scalar) && 
                   (_expressions <= o._expressions) && 
                   (_g <= o._g);
        CRAB_LOG("array-sgraph-domain", crab::outs () << res << "\n";);
        return res;
      }

      void operator|=(array_sgraph_domain_t o)  {
        *this = (*this | o);
      }

      array_sgraph_domain_t operator|(array_sgraph_domain_t &o){
        crab::CrabStats::count (getDomainName() + ".count.join");
        crab::ScopedCrabStats __st__(getDomainName() + ".join");

        CRAB_LOG("array-sgraph-domain",
                 crab::outs () << "Join " << *this << " and "  << o << "=\n");
        array_sgraph_domain_t join(_scalar | o._scalar, 
                                   _expressions | o._expressions,
                                   _g | o._g);
        CRAB_LOG("array-sgraph-domain", crab::outs () << join << "\n";);
        return join;
      }

      template<typename Thresholds>
      array_sgraph_domain_t widening_thresholds (array_sgraph_domain_t& o, 
                                                 const Thresholds & ts) {
        crab::CrabStats::count (getDomainName() + ".count.widening");
        crab::ScopedCrabStats __st__(getDomainName() + ".widening");

          CRAB_LOG("array-sgraph-domain",
                   crab::outs () << "Widening (w/ thresholds) " << *this << " and "  << o << "=\n";);
        auto widen_scalar(_scalar.widening_thresholds(o._scalar,ts));
        auto widen_expr(_expressions.widening_thresholds(o._expressions,ts));
        auto widen_g(_g.widening_thresholds(o._g,ts));
        if (!reduce(widen_scalar, widen_g)) {
          CRAB_LOG("array-sgraph-domain", crab::outs () << "_|_\n";);
          return array_sgraph_domain_t::bottom();
        } else {
          array_sgraph_domain_t widen(widen_scalar, widen_expr, widen_g);
          CRAB_LOG("array-sgraph-domain", crab::outs () << widen << "\n";);
          return widen;
        }
      }

      array_sgraph_domain_t operator||(array_sgraph_domain_t &o){
        crab::CrabStats::count (getDomainName() + ".count.widening");
        crab::ScopedCrabStats __st__(getDomainName() + ".widening");

        CRAB_LOG("array-sgraph-domain",
                 crab::outs () << "Widening " << *this << " and "  << o << "=\n");        
        auto widen_scalar(_scalar || o._scalar);
        auto widen_expr(_expressions || o._expressions);
        auto widen_g(_g || o._g);
        if (!reduce(widen_scalar, widen_g)) {
          CRAB_LOG("array-sgraph-domain", crab::outs () << "_|_\n";);
          return array_sgraph_domain_t::bottom();
        } else {
          array_sgraph_domain_t widen(widen_scalar, widen_expr, widen_g);
          CRAB_LOG("array-sgraph-domain", crab::outs () << widen << "\n";);
          return widen;
        }
      }

      array_sgraph_domain_t operator&(array_sgraph_domain_t &o){
        crab::CrabStats::count (getDomainName() + ".count.meet");
        crab::ScopedCrabStats __st__(getDomainName() + ".meet");

        CRAB_LOG("array-sgraph-domain",
                 crab::outs () << "Meet " << *this << " and "  << o << "=\n");
        auto meet_scalar(_scalar & o._scalar);
        auto meet_expr(_expressions & o._expressions);
        auto meet_g(_g & o._g);
        if (!reduce(meet_scalar, meet_g)) {
          CRAB_LOG("array-sgraph-domain", crab::outs () << "_|_\n";);
          return array_sgraph_domain_t::bottom();
        } else {
          array_sgraph_domain_t meet(meet_scalar, meet_expr, meet_g);
          CRAB_LOG("array-sgraph-domain", crab::outs () << meet << "\n";);
          return meet;
        }
      }

      array_sgraph_domain_t operator&&(array_sgraph_domain_t &o){
        crab::CrabStats::count (getDomainName() + ".count.narrowing");
        crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

        if (is_bottom() || o.is_bottom())
          return bottom();
        else if (is_top ())
          return o;
        else {
          // FIXME: Implement properly
          // Narrowing as a no-op should be sound.
          return *this;
        }
      }

      void operator-=(VariableName v) {
        crab::CrabStats::count (getDomainName() + ".count.forget");
        crab::ScopedCrabStats __st__(getDomainName() + ".forget");

        if (is_bottom())
          return;

        // remove v from scalar and array graph
        _scalar -= v;
        // remove v from expressions
        _expressions -= v;
        
        if (var_landmarks.find(landmark_ref_t(v)) != var_landmarks.end()) {
          array_forget(v);
          // remove v' from scalar and array graph
          forget_prime_var(v);
        }
      }

      void operator+=(linear_constraint_system_t csts) 
      {
        crab::CrabStats::count (getDomainName() + ".count.add_constraints");
        crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

        if (is_bottom()) return;
        
        _scalar += csts;
        _expressions += csts;

        if (!reduce(_scalar, _g)) { // FIXME: incremental version
          set_to_bottom();
          return;
        }
        CRAB_LOG("array-sgraph-domain", 
                 crab::outs() << "Assume("<< csts<< ") --- "<< *this<<"\n";);
      }

      void assign (VariableName x, linear_expression_t e) 
      { assign (x, e, true); }
        
      void assign (VariableName x, linear_expression_t e, bool update_expressions) 
      {
        crab::CrabStats::count (getDomainName() + ".count.assign");
        crab::ScopedCrabStats __st__(getDomainName() + ".assign");

        if (is_bottom()) return;

        // skip x:=x 
        if (auto y = e.get_variable())
          if ((*y).name() == x) return;
        
        _scalar.assign(x, e);
        if (update_expressions) _expressions.assign(x, e);

        landmark_ref_t lm_x(x);
        auto it = var_landmarks.find(lm_x);
        if (it != var_landmarks.end()) {
           array_forget (x);
           // restore the relationship between x and x'
           _scalar.apply(OP_ADDITION, get_var(it->second), x, 1);
           //_scalar -= get_var(it->second);
           //_scalar += make_prime_relation (it->second, lm_x);        
        }

        if (!reduce(_scalar, _g)) { // FIXME: incremental version
          set_to_bottom();
          return;
        }

        CRAB_LOG("array-sgraph-domain", 
                 crab::outs() << "Assign "<<x<<" := "<<e<<" ==> "<<*this<<"\n";);
      }

      void apply (operation_t op, VariableName x, VariableName y, Number z) {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, z);
        assign(x, linear_expression_t(y), false);
        apply_landmark<Number> (op, x, z);

        CRAB_LOG("array-sgraph-domain",
                 crab::outs() << "Apply "<<x<<" := "<<y<<" "<<op<<" "<<z<<" ==> "<<*this<<"\n";); 
      }

      void apply(operation_t op, VariableName x, VariableName y, VariableName z)  {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, z);
        assign(x, linear_expression_t(y), false);
        apply_landmark<VariableName> (op, x, z);

        CRAB_LOG("array-sgraph-domain", 
                 crab::outs() << "Apply "<<x<<" := "<<y<<" "<<op<<" "<<z<<" ==> "<<*this<<"\n";);
      }

      void apply(operation_t op, VariableName x, Number k)  {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, k);
        apply_landmark<Number> (op, x, k);

        CRAB_LOG("array-sgraph-domain",
                 crab::outs() << "Apply "<<x<<" := "<<x<<" "<<op<<" "<<k<<" ==> "<<*this<<"\n";);
      }

      void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
        _expressions.apply (op, x, y, width);
        // assume unlimited precision so width is ignored.
        assign(x, variable_t (y), false);
      }
      
      void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
        _expressions.apply (op, x, k, width);
        // assume unlimited precision so width is ignored.
        assign(x, k, false);
      }

      // bitwise_operators_api      
      void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, z);
        apply_only_scalar(op, x, y, z);
      }
      
      void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, k);
        apply_only_scalar(op, x, y, k);
      }
      
      // division_operators_api
      void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, z);
        apply_only_scalar(op, x, y, z);
      }
      
      void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
        crab::CrabStats::count (getDomainName() + ".count.apply");
        crab::ScopedCrabStats __st__(getDomainName() + ".apply");

        _expressions.apply (op, x, y, k);
        apply_only_scalar(op, x, y, k);
      }

      // lhs := arr[idx];
      void load (VariableName lhs, VariableName arr, VariableName idx, z_number nbytes)
      {
        crab::CrabStats::count (getDomainName() + ".count.load");
        crab::ScopedCrabStats __st__(getDomainName() + ".load");

        VariableName norm_idx = normalize_offset (idx, nbytes);
        Weight w = array_elem (norm_idx);

        // --- XXX: simplification wrt Gange et.al.:
        //     Only non-relational invariants involved arr are passed
        //     from the graph domain to the scalar domain.
        _scalar.set (lhs, w[arr]);
        _expressions.set (lhs, w[arr]);

        /// FIXME: we need to reduce only if the content of an array
        /// cell may be propagated to an array index.
        if (!reduce(_scalar,_g)) { // FIXME: incremental version
          set_to_bottom();
        }

        CRAB_LOG("array-sgraph-domain",
                 crab::outs() << "Array read "<<lhs<<" := "<< arr<<"["<<idx<<"] ==> "
                               << *this <<"\n";);    
      }

      // arr[idx] := val 
      void store (VariableName arr, VariableName idx, linear_expression_t val, z_number nbytes)
      {
        crab::CrabStats::count (getDomainName() + ".count.store");
        crab::ScopedCrabStats __st__(getDomainName() + ".store");

        // --- XXX: simplification wrt Gange et.al.:
        //     Only non-relational invariants are passed from the scalar
        //     domain to the graph domain.
        Weight w;
        w.set(arr, eval_interval(val));

        VariableName norm_idx = normalize_offset (idx, nbytes);
        array_forget(norm_idx, arr);
        array_update(norm_idx, w);

        // XXX: since we do not propagate from the array weights to
        // the scalar domain I think we don't need to reduce here.

        CRAB_LOG("array-sgraph-domain",
                 crab::outs() << "Array write "<<arr<<"["<<idx<<"] := "<<val<< " ==> "
                              << *this <<"\n";);
      }
    
      void write(crab_os& o) {
        #if 1
        NumDom copy_scalar(_scalar);
        array_sgraph_t copy_g(_g);
        // Remove all primed variables for pretty printing
        for(auto lm: var_lm_primes()) {
          copy_scalar -= get_var(lm);
          copy_g -= lm;
        }
        for(auto lm: cst_lm_primes()) {
          copy_scalar -= get_var(lm);
          copy_g -= lm;
        }
        o << "(" << copy_scalar  << ",";
        copy_g.write(o,false);  // we do not print bottom edges
        o << ")";
        //o << "##" << _expressions;
        #else
        o << "(" << _scalar  << ",";
        _g.write(o,true); 
        o << ")";
        #endif 
      }

      // XXX: the array domain is disjunctive so it is not really
      // useful to express it through a conjunction of linear
      // constraints
      linear_constraint_system_t to_linear_constraint_system (){
        CRAB_WARN ("array-sgraph does not implement to_linear_constraint_system");
        return linear_constraint_system_t();
      }

      static string getDomainName () {
        string name ("ArraySparseGraph(" + 
                     NumDom::getDomainName () +  "," +  Weight::getDomainName () + ")");
        return name;
      }

    };

    template<typename NumDom, typename Weight>
    class domain_traits<array_sparse_graph_domain<NumDom,Weight,false> > {

     public:
      // WARNING: assume NumDom::number_t = Weight::number_t and
      //                 NumDom::varname_t = Weight::varname_t
      typedef typename NumDom::number_t N;
      typedef typename NumDom::varname_t V;
      
      typedef array_sparse_graph_domain<NumDom,Weight,false> array_sgraph_domain_t;

      template<class CFG, class VarFactory>
      static void do_initialization (CFG cfg, VarFactory &vfac) {
        array_sgraph_domain_t::do_initialization(cfg, vfac);
      }

      static void expand (array_sgraph_domain_t& inv, V x, V new_x) {
        CRAB_WARN ("TODO array_graph_domain expand not implemented");
      }
    
      static void normalize (array_sgraph_domain_t& inv) {
        CRAB_WARN ("TODO array_graph_domain normalize not implemented");
      }
    
      template <typename Iter>
      static void forget (array_sgraph_domain_t& inv, Iter it, Iter end){
        for (auto v: boost::make_iterator_range(it,end))
        { inv -= v; }
      }

      template <typename Iter>
      static void project (array_sgraph_domain_t& inv, Iter it, Iter end) {
        CRAB_WARN ("TODO array_graph_domain project not implemented");
      }
    };
  
    template <typename NumDom, typename Weight>
    class array_domain_traits<array_sparse_graph_domain<NumDom,Weight,false> > {
        
       public:
        typedef ikos::z_number z_number;
        typedef typename NumDom::number_t number_t;
        typedef typename NumDom::varname_t varname_t;
        typedef typename NumDom::linear_expression_t linear_expression_t;
        typedef array_sparse_graph_domain<NumDom,Weight, false> array_sgraph_domain_t;

        static void array_init (array_sgraph_domain_t& inv, varname_t a, 
                                const vector<z_number> &values) { 
          CRAB_WARN ("TODO array_graph_domain array_init not implemented");
        }
        
        static void assume_array (array_sgraph_domain_t& inv, varname_t a, number_t val) {
          CRAB_WARN ("TODO array_graph_domain assume_array not implemented");
        }
        
        static void assume_array (array_sgraph_domain_t& inv, varname_t a, 
                                  interval<number_t> val) {
          CRAB_WARN ("TODO array_graph_domain assume_array not implemented");
        }
        
        static void array_load (array_sgraph_domain_t& inv, 
                                varname_t lhs, varname_t arr, varname_t idx, 
                                z_number nbytes) {
          inv.load (lhs, arr, idx, nbytes);
        }
        
        static void array_store (array_sgraph_domain_t& inv, varname_t arr, 
                                 varname_t idx, linear_expression_t val,
                                 z_number nbytes, bool /*is_singleton*/) {
          inv.store (arr, idx, val, nbytes);
        }        
    };

    // Static data allocation
    template<class Dom, class Wt, bool IsDistWt>
    boost::unordered_map<landmark_ref<typename Dom::varname_t, typename Dom::number_t>, 
                         landmark_ref<typename Dom::varname_t, typename Dom::number_t> > 
    array_sparse_graph_domain<Dom,Wt,IsDistWt>::var_landmarks;

    template<class Dom, class Wt, bool IsDistWt>
    boost::unordered_map<landmark_ref<typename Dom::varname_t, typename Dom::number_t>, 
                         landmark_ref<typename Dom::varname_t, typename Dom::number_t> > 
    array_sparse_graph_domain<Dom,Wt,IsDistWt>::cst_landmarks;

  } // end namespace domains
 
} // end namespace crab

#endif 