
#include <grid_dataset.h>
#include <grid_dataset_cl.h>
#include <grid_mscomplex.h>

#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>
#include <iostream>

#include <utl.h>

namespace br=boost::range;
namespace ba=boost::adaptors;


namespace grid
{
//  typedef tr1::tuple<cellid_t,cellid_list_ptr_t,cellid_list_ptr_t> cp_mfold_qitem_t;

//  typedef producer_consumer_t<cp_mfold_qitem_t> cp_mfold_que_t;

inline cellid_t cellid_of_pair(mscomplex_t &msc,int i)
{return msc.cellid(msc.pair_idx(i));}

template<eGDIR dir>
inline void collect_contrib_cps(mscomplex_t &msc,cellid_list_t &cplist, int i)
{
  auto cop_ftr = bind(cellid_of_pair,std::ref(msc),_1);
  auto in_rect_ftr = bind(&rect_t::contains_point,&msc.m_rect,_1);

  if(msc.m_rect.contains(msc.cellid(i)))
    cplist.push_back(msc.cellid(i));

  br::copy(boost::make_iterator_range(msc.m_conn[dir][i])
           |ba::map_keys|ba::transformed(cop_ftr)|ba::filtered(in_rect_ftr),
           back_inserter(cplist));

}

template<int dim,eGDIR dir,typename cmp_t>
inline cellid_list_ptr_t compute_mfold(dataset_t &ds,mscomplex_t &msc,int i,cmp_t cmp)
{
  cellid_list_ptr_t mfold(new cellid_list_t);

  cellid_list_t cplist;

  collect_contrib_cps<dir>(msc,cplist,i);

  compute_mfold<dim,dir,false>(cplist.begin(),cplist.end(),ds,*mfold,cmp);

  return mfold;
}

//  void compute_saddle_manifold(dataset_t &ds, mscomplex_t &msc, cp_producer_t &prd,cp_mfold_que_t &que)
//  {
//    BOOST_AUTO(des2_cmp,boost::bind(&dataset_t::compare_cells_pp<2>,&ds,_1,_2));
//    BOOST_AUTO(asc2_cmp,boost::bind(&dataset_t::compare_cells_pp<2>,&ds,_2,_1));
//    BOOST_AUTO(des1_cmp,boost::bind(&dataset_t::compare_cells_pp<1>,&ds,_1,_2));
//    BOOST_AUTO(asc1_cmp,boost::bind(&dataset_t::compare_cells_pp<1>,&ds,_2,_1));


//    for(BOOST_AUTO(it,prd.next()); prd.is_valid(it); it = prd.next())
//    {
//      int i = *it,dim = msc.index(i);

//      cellid_list_ptr_t des,asc;

//      if(dim == 1)
//      {
//        des = compute_mfold<1,GDIR_DES>(ds,msc,i,des1_cmp);
//        asc = compute_mfold<1,GDIR_ASC>(ds,msc,i,asc1_cmp);
//      }
//      else
//      {
//        des = compute_mfold<2,GDIR_DES>(ds,msc,i,des2_cmp);
//        asc = compute_mfold<2,GDIR_ASC>(ds,msc,i,asc2_cmp);
//      }

//      que.put(tr1::make_tuple(msc.cellid(i),des,asc));
//    }
//  }

//  void store_mfold(std::ostream &os,cellid_list_t &cps,int_list_t &offsets, cp_mfold_que_t &que,int n)
//  {
//    int offset = 0;

//    offsets.push_back(offset);

//    while(n-- > 0)
//    {
//      cp_mfold_qitem_t qitem = que.get();

//      BOOST_AUTO(c,tr1::get<0>(qitem));
//      BOOST_AUTO(des,tr1::get<1>(qitem));
//      BOOST_AUTO(asc,tr1::get<2>(qitem));

//      offset += des->size(); offsets.push_back(offset);
//      offset += asc->size(); offsets.push_back(offset);

//      os.write((char*)(void*)des->data(),des->size()*sizeof(cellid_t));
//      os.write((char*)(void*)asc->data(),asc->size()*sizeof(cellid_t));

//      cps.push_back(c);
//    }
//  }

int get_header_size(int num_cps)
{
  return sizeof(rect_t)*3         + // rects
         sizeof(int)              + // num_cps
         sizeof(cellid_t)*num_cps + // cellids
         sizeof(int)*(2*num_cps+1); // offsets
}


void write_header(std::ostream & os, dataset_t &ds,int_list_t & offsets,cellid_list_t &cps)
{
  os.seekp(0,std::ios::beg);

  utl::bin_write(os,ds.m_rect);
  utl::bin_write(os,ds.m_ext_rect);
  utl::bin_write(os,ds.m_domain_rect);

  utl::bin_write(os,(int)cps.size());

  os.write((char*)(void*)cps.data(),cps.size()*sizeof(cellid_t));
  os.write((char*)(void*)offsets.data(),offsets.size()*sizeof(int));
}

//  void save_saddle_mfolds(std::ostream &os,dataset_t &ds,mscomplex_t &msc)
//  {
//    mscomplex_t::filter_t fltr = bind(need_saddle_mfold,boost::cref(msc),_1);

//    cp_producer_t prd(msc.shared_from_this(),fltr);

//    int num_cps = prd.count();

//    cp_mfold_que_t que;

//    boost::thread_group group;

//    for(int tid = 0 ; tid < g_num_threads; ++tid)
//      group.create_thread(bind(compute_saddle_manifold,boost::ref(ds),
//                               boost::ref(msc),boost::ref(prd),boost::ref(que)));
////      compute_saddle_manifold(ds,msc,prd,que);

//    os.seekp(get_header_size(num_cps),ios::beg);

//    cellid_list_t cps;
//    int_list_t    offsets;

//    store_mfold(os,cps,offsets,que,num_cps);

//    write_header(os,ds,offsets,cps);

//    group.join_all();//redundant
//  }

void  dataset_t::saveManifolds(mscomplex_ptr_t msc,const std::string &bn)
{
//    std::ofstream fs((bn+".mfold.bin").c_str());
//    ENSURE(fs.is_open(),"unable to open file");
//    boost::thread_group group;
//    group.create_thread(bind(save_saddle_mfolds,boost::ref(fs),boost::ref(*this),boost::ref(*msc)));

  opencl::update_to_surv_extrema(shared_from_this(),msc);

  {
    int ex_num = num_cells2(rect_t(m_rect.lc()+1,m_rect.uc()-1));
    std::ofstream fs((bn+".max.raw").c_str());
    ENSURE(fs.is_open(),"unable to open file");
    fs.write((char*)(void*)m_owner_maxima.data(),sizeof(int)*ex_num);
    fs.close();
  }

  {
    int ex_num = num_cells2(m_rect);
    std::ofstream fs((bn+".min.raw").c_str());
    ENSURE(fs.is_open(),"unable to open file");
    fs.write((char*)(void*)m_owner_minima.data(),sizeof(int)*ex_num);
    fs.close();
  }

//    group.join_all();
}

template<int dim,eGDIR dir,typename cmp_t,typename Titer>
inline cellid_list_ptr_t compute_mfold(dataset_t &ds,mscomplex_t &msc,Titer b,Titer e,cmp_t cmp)
{

  cellid_list_ptr_t mfold(new cellid_list_t);

  cellid_list_t cplist;

  for_each(b,e,bind(collect_contrib_cps<dir>,std::ref(msc),std::ref(cplist),_1));

  const bool trp = (dim==1 && dir==GDIR_ASC)||(dim==2 && dir==GDIR_DES);

  compute_mfold<dim,dir,trp>(cplist.begin(),cplist.end(),ds,*mfold,cmp);

  return mfold;
}


template<eGDIR dir>
void copy_owner_extrema(int_marray_t &dest,const int_marray_t &src,const dataset_t & ds)
{
  rect_t ex_rct = ds.get_extrema_rect<dir>();

  int num_ex = num_cells2(ex_rct);

  dest.resize(ex_rct.span()/2+1);

  memcpy((void*)dest.data(),(void*)src.data(),num_ex*sizeof(int));
}

void  dataset_t::storeOwnerArrays(int_marray_t &omax,int_marray_t &omin) const
{
  copy_owner_extrema<GDIR_DES>(omax,m_owner_maxima,*this);
  copy_owner_extrema<GDIR_ASC>(omin,m_owner_minima,*this);
}

void  dataset_t::loadOwnerArrays(int_marray_t &omax,int_marray_t &omin)
{
  copy_owner_extrema<GDIR_DES>(m_owner_maxima,omax,*this);
  copy_owner_extrema<GDIR_ASC>(m_owner_minima,omin,*this);
}
}
