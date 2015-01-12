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

inline mscomplex_t::merge_dag_t::node_t
mscomplex_t::merge_dag_t::get_node(int i) const
{
  if(i < 0)
    return node_t(i+get_ncps());
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

//  #pragma omp parallel for
  for(int i = 0; i < ncps; ++i)
  {
    m_cp_geom[ASC][i] = i-ncps;
    m_cp_geom[DES][i] = i-ncps;
  }

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

      BOOST_FOREACH(int r,msc->m_conn[odir][q]|ba::map_keys)
      {
        int rnode         = m_cp_geom[dir][r];
        m_cp_geom[dir][r] = m_nodes.size();
        m_nodes.push_back(node_t(r,rnode,pnode,m_last_hversion));

        ENSURES(get_node(rnode).hversion < m_last_hversion)
            <<"earlier rnode has formed from a later cancellation";
      }
      std::swap(p,q);
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

    node_t gnode = get_node(g);

    if( gnode.hversion <= gver)
      l.push_back(gnode.cpid);
    else if(gnode.base != -1)
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
  utl::bin_write(os,m_last_hversion);
}

//*--------------------------------------------------------------------------*/

void mscomplex_t::merge_dag_t::load_bin(istream &is)
{
  utl::bin_read_vec(is,m_nodes);
  utl::bin_read_vec(is,m_cp_geom[DES]);
  utl::bin_read_vec(is,m_cp_geom[ASC]);
  utl::bin_read(is,m_last_hversion);
}

//*--------------------------------------------------------------------------*/

}

/*===========================================================================*/
