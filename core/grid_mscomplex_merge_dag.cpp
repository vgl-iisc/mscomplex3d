#include <stack>

#include <boost/foreach.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/counting_range.hpp>

#include <grid_mscomplex.h>

using namespace std;

namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{

/*===========================================================================*/

inline const mscomplex_t::merge_dag_t::node_t&
mscomplex_t::merge_dag_t::get_node(int i) const
{
  static node_t defnode;

  if(i < 0)
    return defnode;
  return m_nodes[i];
}

/*---------------------------------------------------------------------------*/

mscomplex_t::merge_dag_t::merge_dag_t():m_last_hversion(0){}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::merge_dag_t::get_ncps() const
{  return m_cp_geom[0].size();}

/*---------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::init(int ncps)
{
  m_cp_geom[ASC].resize(ncps);
  m_cp_geom[DES].resize(ncps);

  br::copy(boost::counting_range(-ncps,0)|ba::reversed,m_cp_geom[ASC].begin());
  br::copy(boost::counting_range(-ncps,0)|ba::reversed,m_cp_geom[DES].begin());
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::clear()
{
  m_cp_geom[DES].clear();
  m_cp_geom[ASC].clear();
  m_nodes.clear();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::update(mscomplex_ptr_t msc)
{
  ENSURES(msc->get_hversion() >= m_last_hversion)
      << "Updates to merge_dag may only be performed on the coarsest version";

  for(; m_last_hversion < msc->get_hversion();)
  {
    int_pair_t pr = msc->m_canc_list[m_last_hversion++]; // ++ is deliberate

    int p = pr[0],q = pr[1];

    ASSERT(order_pair<DES>(msc,pr) == pr);

    for(int dir = 0,odir=1 ; dir < 2; ++dir,--odir)
    {
      int pnode = m_cp_geom[dir][p];

      ENSURES(get_node(pnode).hversion < m_last_hversion)
         <<"earlier pnode has formed from a later cancellation"
         <<SVAR(get_node(pnode).hversion) << SVAR(m_last_hversion);

      if( pnode >= get_ncps()  || msc->m_mfolds[dir][pnode].size() > 0)
      {
        BOOST_FOREACH(int r,msc->m_conn[odir][q]|ba::map_keys)
        {
          int rnode         = m_cp_geom[dir][r];
          m_cp_geom[dir][r] = m_nodes.size();
          m_nodes.push_back(node_t(rnode,pnode,m_last_hversion));

          ENSURES(get_node(rnode).hversion < m_last_hversion)
              <<"earlier rnode has formed from a later cancellation";
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::get_contrib_cps
(std::vector<int> &l, eGDIR dir, int cp, int hver, int gver) const
{
  int g = m_cp_geom[dir][cp];

  while(get_node(g).hversion > hver )
    g = get_node(g).base;

  std::set<int>    visited;

  std::stack<int> stk;

  stk.push(g);

  while(stk.size() != 0 )
  {
    int g = stk.top();
    stk.pop();

    if(visited.count(g) != 0 )
      continue;

    visited.insert(g);

    if( g < 0)
      l.push_back(-g-1);

    node_t gnode = get_node(g);

    if(gnode.base != -1 && gnode.hversion <= gver)
    {
      stk.push(gnode.base);
      stk.push(gnode.other);
    }
  }
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::save_bin(ostream &os) const
{
  utl::bin_write_vec(os,m_nodes);
  utl::bin_write_vec(os,m_cp_geom[DES]);
  utl::bin_write_vec(os,m_cp_geom[ASC]);
}

//*--------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::load_bin(istream &is)
{
  utl::bin_read_vec(is,m_nodes);
  utl::bin_read_vec(is,m_cp_geom[DES]);
  utl::bin_read_vec(is,m_cp_geom[ASC]);
}

//*--------------------------------------------------------------------------*/

}

/*===========================================================================*/
