//#include <
// /range/adaptors.hpp>
//#include <boost/bind.hpp>

#include <grid_mscomplex.h>
#include <ranges>
#include <algorithm>

using namespace std;

//namespace br = boost::range;
//namespace ba = boost::adaptors;

namespace grid
{

/*===========================================================================*/

mscomplex_t::mscomplex_t():
  m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
  m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),
  m_hversion(0),m_merge_dag(new merge_dag_t){resize(0);}

/*---------------------------------------------------------------------------*/

mscomplex_t::mscomplex_t(rect_t r,rect_t e,rect_t d):
  m_rect(r),m_ext_rect(e),m_domain_rect(d),
  m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
  m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),
  m_hversion(0),m_merge_dag(new merge_dag_t){resize(0);}

/*---------------------------------------------------------------------------*/

	mscomplex_t::~mscomplex_t(){clear();}

/*---------------------------------------------------------------------------*/

void mscomplex_t::set_critpt(int i,cellid_t c,char idx,cell_fn_t f,cellid_t v)
{
  m_cp_cellid[i] = c;
  m_cp_vertid[i] = v;
  m_cp_index[i]  = idx;
  m_cp_fn[i]     = f;
}

/*---------------------------------------------------------------------------*/

void  mscomplex_t::resize(int i)
{
  if(i !=0 )
  {
    m_cp_cellid.resize(i,cellid_t(-1,-1,-1));
    m_cp_vertid.resize(i,cellid_t(-1,-1,-1));
    m_cp_index.resize(i,-1);
    m_cp_pair_idx.resize(i,-1);
    m_cp_is_cancelled.resize(i,false);
    m_cp_fn.resize(i);
    m_des_conn.resize(i);
    m_asc_conn.resize(i);
    m_des_mfolds.resize(i);
    m_asc_mfolds.resize(i);
    m_merge_dag->init(i);
  }

  for(int dir = 0; dir < GDIR_CT; ++ dir)
    for(int dim = 0 ; dim <= gc_grid_dim ;++ dim)
      m_geom_hversion[dir][dim] = get_hversion();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::connect_cps(int p, int q,int m)
{
  order_pr_by_cp_index(*this,p,q);

  ENSURES(index(p) == index(q)+1);

  // if a d-cp hits a d+-1 cp and the d+-1 cp is paired
  // then the connection is useful iff the dimension of the pair is d

//  ASSERT(!(is_paired(p) && index(pair_idx(p))!= index(q)));
//  ASSERT(!(is_paired(q) && index(pair_idx(q))!= index(p)));

  if( m_des_conn[p].count(q) == 0)
  {
    ASSERT(m_asc_conn[q].count(p) == 0);
    m_des_conn[p][q] = 0;
    m_asc_conn[q][p] = 0;
  }

  ASSERT(m_des_conn[p][q] == m_asc_conn[q][p]);

  m_des_conn[p][q] += m;
  m_asc_conn[q][p] += m;
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::clear()
{
  m_cp_cellid.clear();
  m_cp_vertid.clear();
  m_cp_pair_idx.clear();
  m_cp_index.clear();
  m_cp_is_cancelled.clear();
  m_cp_fn.clear();
  m_des_conn.clear();
  m_asc_conn.clear();
  m_des_mfolds.clear();
  m_asc_mfolds.clear();
  m_merge_dag->clear();
}

/*---------------------------------------------------------------------------*/

std::string mscomplex_t::cp_conn (int i) const
{
  std::stringstream ss;

  ss<<std::endl<<"des = ";
  /*
  br::copy(m_des_conn[i]|ba::map_keys|ba::transformed(std::bind(&mscomplex_t::cellid,this,std::placeholders::_1)),
           ostream_iterator<cellid_t>(ss));

  ss<<std::endl<<"asc = ";

  br::copy(m_asc_conn[i]|ba::map_keys|ba::transformed(std::bind(&mscomplex_t::cellid,this,std::placeholders::_1)),
           ostream_iterator<cellid_t>(ss));
	*/

  std::ranges::copy(
      m_des_conn[i]
      | std::ranges::views::keys
      | std::ranges::views::transform([this](const auto& key) { return this->cellid(key); }),  // Lambda instead of std::bind
      std::ostream_iterator<cellid_t>(ss));

  ss << std::endl << "asc = ";

  std::ranges::copy(
      m_asc_conn[i]
      | std::ranges::views::keys
      | std::ranges::views::transform([this](const auto& key) { return this->cellid(key); }),  // Lambda
      std::ostream_iterator<cellid_t>(ss));

  ss<<std::endl;

  return ss.str();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::save_bin(ostream &os) const
{
  utl::bin_write(os,m_rect);
  utl::bin_write(os,m_ext_rect);
  utl::bin_write(os,m_domain_rect);

  utl::bin_write_vec(os,m_cp_cellid);
  utl::bin_write_vec(os,m_cp_vertid);
  utl::bin_write_vec(os,m_cp_pair_idx);
  utl::bin_write_vec(os,m_cp_index);
  utl::bin_write_vec(os,m_cp_is_cancelled);
  utl::bin_write_vec(os,m_cp_fn);

  int N = m_cp_cellid.size();

  int_list_t nconn(2*N);
  int_int_list_t adj;

  for(int i = 0 ; i < N; ++i)
  {
    nconn[2*i]   = m_des_conn[i].size();
    nconn[2*i+1] = m_asc_conn[i].size();

    //br::copy(m_des_conn[i],back_inserter(adj));


    //br::copy(m_asc_conn[i],back_inserter(adj));

    copy(m_des_conn[i].begin(), m_des_conn[i].end(), std::back_inserter(adj));
    copy(m_asc_conn[i].begin(), m_asc_conn[i].end(), std::back_inserter(adj));


  }

  utl::bin_write_vec(os,nconn);
  utl::bin_write_vec(os,adj);

  utl::bin_write(os,m_hversion);
  utl::bin_write_vec(os,m_canc_list);
  utl::bin_write(os,m_geom_hversion);

  for(int i = 0 ; i < N; ++i)
  {
    utl::bin_write_vec(os,m_des_mfolds[i]);
    utl::bin_write_vec(os,m_asc_mfolds[i]);
  }

  m_merge_dag->save_bin(os);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::load_bin(istream &is)
{
  clear();

  int_list_t     nconn;
  int_int_list_t adj;

  utl::bin_read(is,m_rect);
  utl::bin_read(is,m_ext_rect);
  utl::bin_read(is,m_domain_rect);

  utl::bin_read_vec(is,m_cp_cellid);
  utl::bin_read_vec(is,m_cp_vertid);
  utl::bin_read_vec(is,m_cp_pair_idx);
  utl::bin_read_vec(is,m_cp_index);
  utl::bin_read_vec(is,m_cp_is_cancelled);
  utl::bin_read_vec(is,m_cp_fn);

  utl::bin_read_vec(is,nconn);
  utl::bin_read_vec(is,adj);

  int N = m_cp_cellid.size();

  m_des_conn.resize(N);
  m_asc_conn.resize(N);

  int_int_list_t::iterator a,b,c = adj.begin();

  for(int i = 0 ; i < N; ++i)
  {
    a = c;
    b = a + (nconn[2*i]);
    c = b + (nconn[2*i+1]);

    copy(a,b,inserter(m_des_conn[i],m_des_conn[i].begin()));
    copy(b,c,inserter(m_asc_conn[i],m_asc_conn[i].begin()));
  }

  utl::bin_read(is,m_hversion);
  utl::bin_read_vec(is,m_canc_list);
  utl::bin_read(is,m_geom_hversion);

  for(int i = 0 ; i < N; ++i)
  {
    utl::bin_read_vec(is,m_des_mfolds[i]);
    utl::bin_read_vec(is,m_asc_mfolds[i]);
  }

  m_merge_dag->load_bin(is);
}

/*===========================================================================*/

}
