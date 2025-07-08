#include <grid_dataset.h>
#include <grid_dataset_cl.h>

#include "grid_dataset.h"
#include "grid_dataset.h"

#include <filesystem>

#ifdef __linux__
#include <cstring>
#endif

using namespace std;

namespace grid
{

/*===========================================================================*/

    dataset_t::dataset_t(const rect_t& r, const rect_t& e, const rect_t& d) :
        m_rect(r),
        m_ext_rect(e),
        m_domain_rect(d),
        m_vert_fns(0,0,0,0),
        m_cell_flags(0, 0, 0, 0),
        m_owner_maxima(0, 0, 0, 0),
        m_owner_minima(0, 0, 0, 0)
{
  // TODO: assert that the given rect is of even size..
  //       since each vertex is in the even positions
}

/*---------------------------------------------------------------------------*/

dataset_t::~dataset_t () {clear();}

/*---------------------------------------------------------------------------*/


template <typename T>
T swap_endian(T value) {
    uint8_t* ptr = reinterpret_cast<uint8_t*>(&value);
    std::reverse(ptr, ptr + sizeof(T));
    return value;
}
void dataset_t::init(const string &filename)
{
  init_storage();

  rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
  uint num_pts   = pt_span[0]*pt_span[1]*pt_span[2];
  auto wd = std::filesystem::current_path();

  std::cout << "\nfile size: " << pt_span[0] <<" x " <<pt_span[1]<<" x "<< pt_span[2] << "\n";

  ifstream ifs(filename.c_str(),ios::in|ios::binary);
  //ENSURE(ifs.is_open(),"unable to open file " + filename);
  ENSUREV2(ifs.is_open(),"unable to open file ", filename,wd);
  m_vert_fns.resize(pt_span[0], pt_span[1], pt_span[2]);
  ifs.read(reinterpret_cast<char*>(m_vert_fns.data.data()), sizeof(cell_fn_t) * num_pts);
  if (ifs.eof()) {
      std::cerr << "Reached EOF after reading." << std::endl;
     
      ifs.clear(); // Reset EOF flag
  }

  ENSURE(ifs.fail()==false,"failed to read some data");
  ifs.seekg(0,ios::end);
  ENSURE(uint(ifs.tellg())==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");
  ifs.close();
  
}

/*---------------------------------------------------------------------------*/

void  dataset_t::init(const cell_fn_t * dptr, bool is_fortran_order)
{
  init_storage();

  rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
  uint num_pts   = pt_span[0]*pt_span[1]*pt_span[2];


  if(is_fortran_order)
    std::memcpy(m_vert_fns.data.data(),dptr,num_pts*sizeof(cell_fn_t));
  else
    for(int x = 0 ; x < pt_span[0] ; ++ x)
      for(int y = 0 ; y < pt_span[1] ; ++ y)
        for(int z = 0 ; z < pt_span[2] ; ++ z)
          m_vert_fns(x,y,z) = *dptr++;


}

/*---------------------------------------------------------------------------*/

/**
 * \brief Initialise storage of scalar field data, this function assumes the data is in a 3D space
 */
void dataset_t::init_storage()
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  rect_size_t   span   = m_ext_rect.span() + 1;
  rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
  rect_point_t bl = m_ext_rect.lower_corner();

  m_cell_flags.resize(span[0],span[1],span[2]);
  m_vert_fns.resize(pt_span);

  //std::fill_n(m_cell_flags.data(),span[0]*span[1]*span[2],0);

  std::fill_n(m_cell_flags.data.data(),span[0]*span[1]*span[2],0);
  
  m_cell_flags.reindex(bl);
  m_vert_fns.reindex(bl/2);

  if(opencl::is_gpu_context())
  {
    m_owner_maxima.resize(m_rect.span()/2);
    m_owner_minima.resize((m_rect.span()/2)+1);

    m_owner_maxima.reindex(m_rect.lc()/2);
    m_owner_minima.reindex(m_rect.lc()/2);
  }
}

/*---------------------------------------------------------------------------*/

/**
 * \brief Clear scalar field data, resize all arrays to 0
 */
void  dataset_t::clear()
{
  m_cell_flags.resize(cellid_t::zero);
  m_vert_fns.resize(cellid_t::zero);
  m_owner_maxima.resize(cellid_t::zero);
  m_owner_minima.resize(cellid_t::zero);
}

/*---------------------------------------------------------------------------*/
/**
 * \brief  
 * \param c 
 * \return 
 */
inline cellid_t dataset_t::get_cell_vert(cellid_t c) const
{
  cellid_t v = c;

  switch(getCellDim(c))
  {
    case 3: v = getCellMaxFacetId(v);
    case 2: v = getCellMaxFacetId(v);
    case 1: v = getCellMaxFacetId(v);
  }
  return v;
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Gets the neighboring cell specified by flag 
 * \param c cell
 * \param f flag
 * \return neighboring cell
 */
inline cellid_t get_neighbor(cellid_t c,cell_flag_t f)
{
    //cell_flag_t d = f & 0x07;
    //c[(d - 1) >> 1] += (d & 1) ? (-1) : (+1);


  cell_flag_t d = f&0x07;
  int8_t axis = (d - 1) >> 1;
  int8_t step = (d & 
      (uint8_t)1) ? (-1) : (+1);
  //if (d > 0)step = -1;
  //else step = 1;
  c[axis] += step;
  return c;
}

/*---------------------------------------------------------------------------*/

inline uint dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
{
  for(uint i = 0; i < gc_grid_dim; ++i)
  {
    inc[i*2+0] = c;
    inc[i*2+1] = c;

    inc[i*2+0][i] -= 1;
    inc[i*2+0][i] += 1;
  }
  return gc_grid_dim*2;
}

/*---------------------------------------------------------------------------*/

/**
 * \brief get the flag of a neighboring cell
 * \param c given cell
 * \param fct neighboring cell
 * \return the neighboring cel flag value
 */
inline cell_flag_t get_neighbor_flag(cellid_t c,cellid_t fct)
{
  ASSERT(euclid_norm2(c-fct) == 1);

  int d = 0;

  if(c[1] != fct[1])
    d = 1;
  else if(c[2] != fct[2])
    d = 2;

  if(c[d] > fct[d])
    return (1 + d*2 + 0);
  else
    return (1 + d*2 + 1);
}

/*---------------------------------------------------------------------------*/

inline cellid_t flag_to_pair(cellid_t c,cell_flag_t f)
{return get_neighbor(c,(f>>3)&0x07);}

/*---------------------------------------------------------------------------*/

inline cell_flag_t pair_to_flag(cellid_t c,cellid_t p)
{return (get_neighbor_flag(c,p)<<3);}

/*---------------------------------------------------------------------------*/

bool dataset_t::areCellsIncident(cellid_t c1,cellid_t c2) const
{return ( euclid_norm2(c1-c2) == 1);}

/*---------------------------------------------------------------------------*/

cellid_t dataset_t::getCellPairId (cellid_t c) const
{
	ENSURES(isCellPaired(c));
    auto isPaired = isCellPaired(c);
    auto f = m_cell_flags(c);
    auto pair = flag_to_pair(c,f);

    if(!this->m_rect.contains(pair))
    {
	 std::cout << "Pair is not contained within the rect "
         << "cell =" << c
         << "pair =" << pair
         << "flag = " << f
         << flag_to_pair(c, f);

    }

    ENSURES(this->m_rect.contains(pair)) << "Pair is not contained within the rect "
        << "\ncell =" << c
        << "\npair =" << pair
        << "\nflag = " << f
        << flag_to_pair(c, f);

    	return pair;
}

/*---------------------------------------------------------------------------*/

cellid_t dataset_t::getCellMaxFacetId (cellid_t c) const
    {
	    return get_neighbor(c,m_cell_flags(c));
    }

/*---------------------------------------------------------------------------*/

cellid_t dataset_t::getCellSecondMaxFacetId (cellid_t c) const
{return (2*c - get_neighbor(c,m_cell_flags(c)));}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  static_assert(((cell_coord_t)-1) < 0 , "coord_t needs to support -1 ");

  uint pos = 0;

  cellid_t i;

  for(i[2] = -(c[2]&1) ; i[2] <= (c[2]&1) ;i[2]+=2)
  {
    for(i[1] = -(c[1]&1) ; i[1] <= (c[1]&1) ;i[1]+=2)
    {
      for(i[0] = -(c[0]&1) ; i[0] <= (c[0]&1) ;i[0]+=2)
      {
        p[pos++] = c+i;
      }
    }
  }

  return (1<<getCellDim (c));
}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellCubes (cellid_t c,cellid_t  *p) const
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  static_assert(((cell_coord_t)-1) < 0 , "coord_t needs to support -1 ");

  uint pos = 0;

  cellid_t i;

  for(    i[2] = -((c[2]+1)&1) ; i[2] <= ((c[2]+1)&1) ;i[2]+=2)
    for(  i[1] = -((c[1]+1)&1) ; i[1] <= ((c[1]+1)&1) ;i[1]+=2)
      for(i[0] = -((c[0]+1)&1) ; i[0] <= ((c[0]+1)&1) ;i[0]+=2)
        if(m_ext_rect.contains(c+i))
          p[pos++] = c+i;

  return pos;
}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
{
  uint pos = 0;

  for(uint d = 0; d< gc_grid_dim; ++d)
  {
    for(uint i = 0 ; i < (c[d]&1);++i)
    {

      f[pos] = c; f[pos++][d] += 1;
      f[pos] = c; f[pos++][d] -= 1;
    }
  }
  return getCellDim (c)*2;
}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
{
  uint cf_ct = (gc_grid_dim - getCellDim (c))*2 ;

  uint pos = 0;

  for(uint d = 0; d< gc_grid_dim; ++d)
  {
    for(uint i = 0 ; i < ((c[d]+1)&1);++i)
    {

      cf[pos] = c; cf[pos++][d] += 1;
      cf[pos] = c; cf[pos++][d] -= 1;
    }
  }

  uint cf_nv_pos = 0;

  for (uint i = 0 ;i < cf_ct;++i)
    if (m_ext_rect.contains (cf[i]))
      cf[cf_nv_pos++] = cf[i];

  return cf_nv_pos;

}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellCofaces (cellid_t c,cellid_t *cf) const
{
  uint cf_ct = std::pow(3,(gc_grid_dim - getCellDim (c))) ;

  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  static_assert(((cell_coord_t)-1) < 0 , "coord_t needs to support -1 ");

  uint pos = 0;

  cellid_t i,l = (c+cellid_t(1,1,1))&(cellid_t(1,1,1));

  for(i[0] = -(l[0]) ; i[0] <= (l[0]) ;i[0]+=1)
  {
    for(i[1] = -(l[1]) ; i[1] <= (l[1]) ;i[1]+=1)
    {
      for(i[2] = -(l[2]) ; i[2] <= (l[2]) ;i[2]+=1)
      {
        cf[pos++] = c+i;
      }
    }
  }

  uint cf_nv_pos = 0;

  for (uint i = 0 ;i < cf_ct;++i)
    if (m_ext_rect.contains (cf[i]))
      cf[cf_nv_pos++] = cf[i];

  return cf_nv_pos;
}

/*---------------------------------------------------------------------------*/

uint dataset_t::getCellEst (cellid_t c,cellid_t* est)  const
{
  cellid_t cfs[40];

  int cfs_ct = getCellCofaces(c,cfs);

  ASSERT(is_in_range(cfs_ct,0,40));

  uint pos = 0;

  for(int i = 0 ; i< cfs_ct;++i)
  {
    cellid_t cf = cfs[i];

    ASSERT(m_ext_rect.contains(cf));

    uint c_dim  = getCellDim(c);
    uint cf_dim = getCellDim(cf);

    for(uint j = c_dim; j < cf_dim ; ++j)
      cf = getCellMaxFacetId(cf);

    if(cf == c)
      est[pos++] = cfs[i];
  }

  return pos;
}

/*---------------------------------------------------------------------------*/

bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
{return (getCellDim (c) <getCellDim (p));}

/*---------------------------------------------------------------------------*/

//  bool dataset_t::isCellMarked (cellid_t c) const
//  {
//    return ! (m_cell_flags (c) == CELLFLAG_UNKNOWN);
//  }

/*---------------------------------------------------------------------------*/

bool dataset_t::isCellCritical (cellid_t c) const
	{

        return (m_cell_flags(c) & CELLFLAG_CRITICAL);
	}

/*---------------------------------------------------------------------------*/

bool dataset_t::isCellPaired (cellid_t c) const
{return (((m_cell_flags(c)>>3) & 0x07) !=0);}

/*---------------------------------------------------------------------------*/

bool dataset_t::isCellVisited (cellid_t c) const
{return (m_cell_flags (c) & CELLFLAG_VISITED);}

/*---------------------------------------------------------------------------*/

void dataset_t::pairCells (cellid_t c,cellid_t p)
{
  ASSERT(isCellPaired(c) == false);
  ASSERT(isCellPaired(p) == false);

  m_cell_flags (c) |= (pair_to_flag(c,p));
  m_cell_flags (p) |= (pair_to_flag(p,c));

  ASSERT(getCellPairId(c) == p);
  ASSERT(getCellPairId(p) == c);
}

/*---------------------------------------------------------------------------*/

void dataset_t::visitCell(cellid_t c)
{
  //    ASSERT(isCellVisited(c) == false);
  m_cell_flags (c) |= CELLFLAG_VISITED;
  ASSERT(isCellVisited(c) == true);
}

/*---------------------------------------------------------------------------*/

//  void  dataset_t::unpairCells ( cellid_t c,cellid_t p )
//  {
//    ASSERT(getCellPairId(c) == p);
//    ASSERT(getCellPairId(p) == c);

//    m_cell_pairs (c) = CELLADJDIR_UNKNOWN;
//    m_cell_pairs (p) = CELLADJDIR_UNKNOWN;

//    m_cell_flags (c) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);
//    m_cell_flags (p) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);

//    ASSERT(isCellPaired(c) == false);
//    ASSERT(isCellPaired(p) == false);
//  }

/*---------------------------------------------------------------------------*/

void dataset_t::setCellMaxFacet (cellid_t c,cellid_t f)
{
  ASSERT(getCellDim(c) == getCellDim(f)+1);
  m_cell_flags (c) |= get_neighbor_flag(c,f);
  ASSERT(getCellMaxFacetId(c) == f);
}

/*---------------------------------------------------------------------------*/

void dataset_t::markCellCritical (cellid_t c)
{
  ASSERT(isCellCritical(c) == false);
  m_cell_flags (c) |= CELLFLAG_CRITICAL;
  ASSERT(isCellCritical(c) == true);
}

/*---------------------------------------------------------------------------*/

bool dataset_t::isTrueBoundryCell (cellid_t c) const
{return (m_domain_rect.isOnBoundry (c));}

/*---------------------------------------------------------------------------*/

bool dataset_t::isFakeBoundryCell (cellid_t c) const
{return (m_rect.isOnBoundry (c) && (!m_domain_rect.isOnBoundry (c)));}

/*---------------------------------------------------------------------------*/

bool dataset_t::isCellExterior (cellid_t c) const
{return (!m_rect.contains (c));}

/*---------------------------------------------------------------------------*/

void dataset_t::save_bin(ostream &os) const
{
  utl::bin_write(os,m_rect);
  utl::bin_write(os,m_ext_rect);
  utl::bin_write(os,m_domain_rect);

  rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
  uint npts            = pt_span[0]*pt_span[1]*pt_span[2];

  utl::bin_write_raw(os,m_vert_fns.data.data(),npts);
}

/*---------------------------------------------------------------------------*/

void dataset_t::load_bin(istream &is)
{
  utl::bin_read(is,m_rect);
  utl::bin_read(is,m_ext_rect);
  utl::bin_read(is,m_domain_rect);

  init_storage();

  rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
  uint npts            = pt_span[0]*pt_span[1]*pt_span[2];

  utl::bin_read_raw(is,m_vert_fns.data.data(),npts);

  compute_owner_grad();
}

/*===========================================================================*/



/*===========================================================================*/

void dataset_t::log_flags()
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
  {
    for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
    {
      for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
      {
        std::cout<<m_cell_flags(c)<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  }
}

/*---------------------------------------------------------------------------*/

/**
 * \brief Provide gradient direction glyph for \param c .
 *
 * E.g for a vertex and edge pair along the positive x direction it would provide "-" and ">"
 * for the vertex and edge respectively
 * when you print it out, you'll see the gradient arrow ->
 *
 * \param c cell for which glyph is given.
 * \param p cell's pair
 * \return direction glyph for gradient of cell
 *
 */
char get_dir_txt(cellid_t c,cellid_t p)
{
  
  int dir = 0;

  if (c[1] != p[1])
    dir = 1;
  else if( c[2] != p[2])
    dir = 2;

  if(dir == 2)
  {
    if( c[dir] > p[dir])
      return 'd';
    else
      return 'u';
  }

  if(c[dir]&1)
  {
    if(c[dir] > p[dir])
    {
      switch(dir)
      {
      case 0: return '>';
      case 1: return 'v';
      }
    }
    else
    {
      switch(dir)
      {
      case 0: return '<';
      case 1: return '^';
      }
    }
  }
  else
  {
    switch(dir)
    {
    case 0: return '-';
    case 1: return '|';
    }
  }

  return '#';
}

/*---------------------------------------------------------------------------*/
//
void dataset_t::log_pairs(std::ostream &os)
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
  {
    os<<"sheet no:: "<<c[2]<<std::endl;
    for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
    {
      for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
      {
        if(isCellCritical(c))
        {
          if(isCellPaired(c))
            os<<get_dir_txt(c,getCellPairId(c))<<"c";
          else
            os<<"C ";
        }
        else if(isCellPaired(c))
          os<<get_dir_txt(c,getCellPairId(c))<<" ";
        else
          os<<"? ";
        os.flush();
      }
      os<<std::endl;
    }

  }
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_visits(std::ostream &os)
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
  {
    os<<"sheet no:: "<<c[2]<<std::endl;
    for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
    {
      for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
      {
        if(isCellVisited(c))
          os<<"x ";
        else
          os<<". ";
      }
      os<<std::endl;
    }
  }
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_pair_visits(std::ostream &os)
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
  {
    os<<"sheet no:: "<<c[2]<<std::endl;
    for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
    {
      for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
      {
        if(isCellVisited(c)&&isCellPaired(c)&&isCellVisited(getCellPairId(c)))
          os<<"x ";
        else
          os<<". ";
      }
      os<<std::endl;
    }
  }
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_owner_extrema(eGDIR dir, std::ostream &os)
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  rect_t ex_rect               = (dir == GDIR_DES)?(rect_t(m_rect.lc()+1,m_rect.uc()-1)):(m_rect);
  int_marray_t &owner_extrema  = (dir == GDIR_DES)?(m_owner_maxima):(m_owner_minima);

  int X = ex_rect[0].span()/2+1;
  int Y = ex_rect[1].span()/2+1;
  int N = num_cells2(ex_rect);

  for( int i = 0 ; i < N; ++i)
  {
    cellid_t c = i_to_c2(ex_rect,i);

    os<<owner_extrema(c/2)<<" ";

    if((i+1)%X == 0)
    {
      os<<endl;

      if(((i/X)+1)%Y == 0)
        os<<endl;
    }
  }

//    for(c[2] = ex_rect[2][0] ; c[2] <= ex_rect[2][1]; c[2]+=2)
//    {
//      os<<"sheet no:: "<<c[2]<<std::endl;
//      for(c[1] = ex_rect[1][0] ; c[1] <= ex_rect[1][1]; c[1]+=2)
//      {
//        for(c[0] = ex_rect[0][0] ; c[0] <= ex_rect[0][1]; c[0]+=2)
//        {
//          os<<m_owner_extrema[dir](c/2)<<" ";
//        }
//        os<<std::endl;
//      }
//    }
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_pairs(const std::string &s)
{
  std::ofstream fs(s.c_str());
  ENSURE(fs.is_open(),"unable to open file");
  log_pairs(fs);
  fs.close();
}

/*---------------------------------------------------------------------------*/
void dataset_t::verify_pairs(std::ostream& os)
{
    static_assert(gc_grid_dim == 3, "defined for 3-manifolds only");

    cellid_t c;

    for (c[2] = m_rect[2][0]; c[2] <= m_rect[2][1]; ++c[2])
    {
        os << "sheet no:: " << c[2] << std::endl;
        for (c[1] = m_rect[1][0]; c[1] <= m_rect[1][1]; ++c[1])
        {
            for (c[0] = m_rect[0][0]; c[0] <= m_rect[0][1]; ++c[0])
            {
                if (isCellCritical(c))
                {
                    if (isCellPaired(c))
                        os << get_dir_txt(c, getCellPairId(c)) << "c";
                    else
                        os << "C ";
                }
                else if (isCellPaired(c))
                    os << get_dir_txt(c, getCellPairId(c)) << " ";
                else
                    os << "? ";
            }
            os << std::endl;
        }

    }

}

/*---------------------------------------------------------------------------*/

void dataset_t::log_pair_visits(const std::string &s)
{
  std::ofstream fs(s.c_str());
  ENSURE(fs.is_open(),"unable to open file");
  log_pair_visits(fs);
  fs.close();
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_visits(const std::string &s)
{
  std::ofstream fs(s.c_str());
  ENSURE(fs.is_open(),"unable to open file");
  log_visits(fs);
  fs.close();
}

/*---------------------------------------------------------------------------*/

void dataset_t::log_max_facets()
{
  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
  {
    for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
    {
      for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
      {
        if(getCellDim(c) != 0 )
          std::cout<< getCellMaxFacetId(c);
        else
          std::cout<< "(.,.,.)";
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  }
}

/*---------------------------------------------------------------------------*/

void dataset_t::extract_vdata_subarray(rect_t r,const std::string &filename)
{
  if(r.lower_corner()%2 != cellid_t::zero ||
     r.upper_corner()%2 != cellid_t::zero )
  {
    throw std::runtime_error("r must specify an aabb with vertex end pts");
  }

  std::ofstream ofs(filename.c_str(),std::ios::out|std::ios::binary);

  if(ofs.is_open() == false)
    throw std::runtime_error("unable to open file");

  static_assert(gc_grid_dim == 3 , "defined for 3-manifolds only");

  cellid_t c;

  for(c[2] = r[2][0] ; c[2] <= r[2][1]; c[2] +=2)
  {
    for(c[1] = r[1][0] ; c[1] <= r[1][1]; c[1] +=2)
    {
      for(c[0] = r[0][0] ; c[0] <= r[0][1]; c[0] +=2)
      {
        cell_fn_t fn = get_cell_fn(c);

        ofs.write((char*)(void*)&fn,sizeof(cell_fn_t));
      }
    }
  }
}
/*===========================================================================*/



/*===========================================================================*/

void get_boundry_rects(const rect_t &r,const rect_t & e,rect_list_t &bnds)
{
  for( int xyz_dir = 0 ; xyz_dir < 3; ++xyz_dir)
  {
    for( int lr_dir = 0 ; lr_dir < 2; ++lr_dir)
    {
      rect_t bnd = r;

      if(r[xyz_dir][lr_dir] != e[xyz_dir][lr_dir])
      {
        bnd[xyz_dir][0] = r[xyz_dir][lr_dir];
        bnd[xyz_dir][1] = r[xyz_dir][lr_dir];

        bnds.push_back(bnd);
      }
    }
  }
}

/*===========================================================================*/
}
