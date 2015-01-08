#include <boost/range/adaptors.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/array.hpp>


#include <grid_mscomplex.h>

using namespace std;

namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{

/*===========================================================================*/

mscomplex_t::mscomplex_t():
  m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
  m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),m_canc_pos(0){}

/*---------------------------------------------------------------------------*/

mscomplex_t::mscomplex_t(rect_t r,rect_t e,rect_t d):
  m_rect(r),m_ext_rect(e),m_domain_rect(d),
  m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
  m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),m_canc_pos(0){}

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
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::connect_cps(int p, int q,int m)
{
  order_pr_by_cp_index(*this,p,q);

  ENSURES(index(p) == index(q)+1);

  // if a d-cp hits a d+-1 cp and the d+-1 cp is paired
  // then the connection is useful iff the dimension of the pair is d

  ASSERT(!(is_paired(p) && index(pair_idx(p))!= index(q)));
  ASSERT(!(is_paired(q) && index(pair_idx(q))!= index(p)));

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

}

/*---------------------------------------------------------------------------*/

std::string mscomplex_t::cp_conn (int i) const
{
  std::stringstream ss;

  ss<<std::endl<<"des = ";

  br::copy(m_des_conn[i]|ba::map_keys|ba::transformed(bind(&mscomplex_t::cellid,this,_1)),
           ostream_iterator<cellid_t>(ss));

  ss<<std::endl<<"asc = ";

  br::copy(m_asc_conn[i]|ba::map_keys|ba::transformed(bind(&mscomplex_t::cellid,this,_1)),
           ostream_iterator<cellid_t>(ss));

  ss<<std::endl;

  return ss.str();
}

template<class Archive>
void mscomplex_t::serialize(Archive & ar, const unsigned int version)
{
  ar& BOOST_SERIALIZATION_NVP(m_rect);
  ar& BOOST_SERIALIZATION_NVP(m_ext_rect);
  ar& BOOST_SERIALIZATION_NVP(m_domain_rect);
  ar& BOOST_SERIALIZATION_NVP(m_canc_pos);

  ar& BOOST_SERIALIZATION_NVP(m_cp_cellid);
  ar& BOOST_SERIALIZATION_NVP(m_cp_vertid);
  ar& BOOST_SERIALIZATION_NVP(m_cp_pair_idx);
  ar& BOOST_SERIALIZATION_NVP(m_cp_index);
//  ar& BOOST_SERIALIZATION_NVP(m_cp_cancno);
//  ar& BOOST_SERIALIZATION_NVP(m_cp_is_boundry);
  ar& BOOST_SERIALIZATION_NVP(m_cp_is_cancelled);
  ar& BOOST_SERIALIZATION_NVP(m_cp_fn);
  ar& BOOST_SERIALIZATION_NVP(m_canc_list);
  ar& BOOST_SERIALIZATION_NVP(m_des_conn);
  ar& BOOST_SERIALIZATION_NVP(m_asc_conn);
  ar& BOOST_SERIALIZATION_NVP(m_des_mfolds);
  ar& BOOST_SERIALIZATION_NVP(m_asc_mfolds);
//  ar& BOOST_SERIALIZATION_NVP(m_multires_version);
//  ar& BOOST_SERIALIZATION_NVP(m_fmax);
//  ar& BOOST_SERIALIZATION_NVP(m_fmin);
//  ar& BOOST_SERIALIZATION_NVP(*m_merge_dag);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::save_bin(ostream &os) const
{
  boost::archive::binary_oarchive oa(os);
  oa << BOOST_SERIALIZATION_NVP(*this);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::load_bin(istream &is)
{
  boost::archive::binary_iarchive ia(is);
  ia >> BOOST_SERIALIZATION_NVP(*this);
}

/*===========================================================================*/

}
