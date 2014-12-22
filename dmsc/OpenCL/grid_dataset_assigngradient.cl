inline bool compare_verts
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  cell_t v1, cell_t v2)
{
  int4 p1 = to_int4(v1-ds.e.lc)/2;
  int4 p2 = to_int4(v2-ds.e.lc)/2;

  func_t f1 = read_imagef(func_img, func_sampler, p1).x;
  func_t f2 = read_imagef(func_img, func_sampler, p2).x;

  if (f1 != f2)
    return f1 < f2;

  return __compare_cells(v1,v2);
}

inline bool compare_edges
( __read_only image3d_t func_img,
  __read_only image3d_t flag_img,
  const dataset_t ds,
  cell_t e1, cell_t e2)
{
  flag_t fg1 = read_imageui(flag_img, flag_sampler, to_int4(e1-ds.e.lc)).x;
  flag_t fg2 = read_imageui(flag_img, flag_sampler, to_int4(e2-ds.e.lc)).x;

  cell_t v1 = flag_to_mxfct(e1,fg1);
  cell_t v2 = flag_to_mxfct(e2,fg2);

  if (is_same_cell(v1,v2))
  {
    v1 = second_max_facet(e1,v1);
    v2 = second_max_facet(e2,v2);

    int boundry_ct1 = boundryCount(ds.d,v1);
    int boundry_ct2 = boundryCount(ds.d,v2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct2 < boundry_ct1);
  }

  return compare_verts(func_img,flag_img,ds,v1,v2);
}

inline bool compare_faces
( __read_only image3d_t func_img,
  __read_only image3d_t flag_img,
  const dataset_t ds,
  cell_t f1, cell_t f2)
{
  flag_t fg1 = read_imageui(flag_img, flag_sampler, to_int4(f1-ds.e.lc)).x;
  flag_t fg2 = read_imageui(flag_img, flag_sampler, to_int4(f2-ds.e.lc)).x;

  cell_t e1 = flag_to_mxfct(f1,fg1);
  cell_t e2 = flag_to_mxfct(f2,fg2);

  if (is_same_cell(e1,e2))
  {
    e1 = second_max_facet(f1,e1);
    e2 = second_max_facet(f2,e2);

    int boundry_ct1 = boundryCount(ds.d,e1);
    int boundry_ct2 = boundryCount(ds.d,e2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct2 < boundry_ct1);
  }

  return compare_edges(func_img,flag_img,ds,e1,e2);
}

inline bool compare_cubes
( __read_only image3d_t func_img,
  __read_only image3d_t flag_img,
  const dataset_t ds,
  cell_t c1, cell_t c2)
{
  flag_t fg1 = read_imageui(flag_img, flag_sampler, to_int4(c1-ds.e.lc)).x;
  flag_t fg2 = read_imageui(flag_img, flag_sampler, to_int4(c2-ds.e.lc)).x;

  cell_t f1 = flag_to_mxfct(c1,fg1);
  cell_t f2 = flag_to_mxfct(c2,fg2);

  if (is_same_cell(f1,f2))
  {
    f1 = second_max_facet(c1,f1);
    f2 = second_max_facet(c2,f2);

    int boundry_ct1 = boundryCount(ds.d,f1);
    int boundry_ct2 = boundryCount(ds.d,f2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct2 < boundry_ct1);
  }

  return compare_faces(func_img,flag_img,ds,f1,f2);
}

inline void set_mxfct_edge
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t e,
  const cell_t d,
  __global flag_t * flag_buf)
{
  if(!contains(ds.e,e))
    return ;

  cell_t v1 = e - d;
  cell_t v2 = e + d;

  cell_t v = v1;

  if( compare_verts(func_img,flag_img,ds,v,v2))
    v = v2;

  flag_buf[c_to_i(ds.e,e)] = mxfct_to_flag(e,v);
}

inline void set_mxfct_face
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t f,
  const cell_t d0,
  const cell_t d1,
  __global flag_t * flag_buf)
{
  if(!contains(ds.e,f))
    return;

  cell_t e0 = f - d0;
  cell_t e1 = f + d0;
  cell_t e2 = f - d1;
  cell_t e3 = f + d1;

  cell_t e = e0;

  if( compare_edges(func_img,flag_img,ds,e,e1))
    e = e1;

  if( compare_edges(func_img,flag_img,ds,e,e2))
    e = e2;

  if( compare_edges(func_img,flag_img,ds,e,e3))
    e = e3;

  flag_buf[c_to_i(ds.e,f)] =  mxfct_to_flag(f,e);
}

inline void set_mxfct_cube
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t c,
  __global flag_t * flag_buf)
{
  cell_t f0 = c - Xdir;
  cell_t f1 = c + Xdir;
  cell_t f2 = c - Ydir;
  cell_t f3 = c + Ydir;
  cell_t f4 = c - Zdir;
  cell_t f5 = c + Zdir;

  cell_t f = f0;

  if( compare_faces(func_img,flag_img,ds,f,f1))
    f = f1;

  if( compare_faces(func_img,flag_img,ds,f,f2))
    f = f2;

  if( compare_faces(func_img,flag_img,ds,f,f3))
    f = f3;

  if( compare_faces(func_img,flag_img,ds,f,f4))
    f = f4;

  if( compare_faces(func_img,flag_img,ds,f,f5))
    f = f5;

  flag_buf[c_to_i(ds.e,c)] =  mxfct_to_flag(c,f);
}

__kernel void assign_max_facet_edge
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);

  int N = num_cells2(ds.e);

  for( int i = tid ; i < N; i += num_thds)
  {
    cell_t v = i_to_c2(ds.e,i);

    flag_buf[c_to_i(ds.e,v)] = 0;

    set_mxfct_edge(func_img,flag_img,ds,v+Xdir,Xdir,flag_buf);
    set_mxfct_edge(func_img,flag_img,ds,v+Ydir,Ydir,flag_buf);
    set_mxfct_edge(func_img,flag_img,ds,v+Zdir,Zdir,flag_buf);
  }
}

__kernel void assign_max_facet_face
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);

  int N = num_cells2(ds.e);

  for( int i = tid ; i < N; i += num_thds)
  {
    cell_t v = i_to_c2(ds.e,i);

    set_mxfct_face(func_img,flag_img,ds,v+XYdir,Xdir,Ydir,flag_buf);
    set_mxfct_face(func_img,flag_img,ds,v+YZdir,Ydir,Zdir,flag_buf);
    set_mxfct_face(func_img,flag_img,ds,v+ZXdir,Zdir,Xdir,flag_buf);
  }
}

__kernel void assign_max_facet_cube
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);
  rect_t cube_rect = make_rect(ds.e.lc + XYZdir,ds.e.uc - XYZdir);

  int N = num_cells2(cube_rect);

  for( int i = tid ; i < N; i += num_thds)
  {
    cell_t c = i_to_c2(cube_rect,i);

    set_mxfct_cube(func_img,flag_img,ds,c,flag_buf);
  }
}

inline bool is_pairable
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  cell_t p,cell_t q)
{
  if(!contains(ds.e,q))
    return false;

  if(is_boundry(ds.d,p) != is_boundry(ds.d,q))
    return false;

  flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(q-ds.e.lc)).x;
  return (is_same_cell(flag_to_mxfct(q,fg),p));
}

inline void mark_pair
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  cell_t p,cell_t q,
  __global flag_t   * flag_buf)
{
  flag_t pfg = read_imageui(flag_img, flag_sampler, to_int4(p-ds.e.lc)).x;
  flag_t qfg = read_imageui(flag_img, flag_sampler, to_int4(q-ds.e.lc)).x;

  pfg |= pair_to_flag(p,q);
  qfg |= pair_to_flag(q,p);

  flag_buf[c_to_i(ds.e,p)] = pfg;
  flag_buf[c_to_i(ds.e,q)] = qfg;
}

inline void set_pair_vert
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t v,
  __global flag_t * flag_buf)
{
  if(!contains(ds.wr,v))
    return;

  cell_t e0 = v - Xdir;
  cell_t e1 = v + Xdir;
  cell_t e2 = v - Ydir;
  cell_t e3 = v + Ydir;
  cell_t e4 = v - Zdir;
  cell_t e5 = v + Zdir;

  cell_t e = invalid_cell;

  if (is_pairable(func_img,flag_img,ds,v,e0))
    e = e0;

  if ((is_pairable(func_img,flag_img,ds,v,e1)) &&
      (is_same_cell(e,invalid_cell) || compare_edges(func_img,flag_img,ds,e1,e)))
    e = e1;

  if (is_pairable(func_img,flag_img,ds,v,e2) &&
      (is_same_cell(e,invalid_cell) || compare_edges(func_img,flag_img,ds,e2,e)))
    e = e2;

  if (is_pairable(func_img,flag_img,ds,v,e3) &&
      (is_same_cell(e,invalid_cell) || compare_edges(func_img,flag_img,ds,e3,e)))
    e = e3;

  if (is_pairable(func_img,flag_img,ds,v,e4) &&
      (is_same_cell(e,invalid_cell) || compare_edges(func_img,flag_img,ds,e4,e)))
    e = e4;

  if ((is_pairable(func_img,flag_img,ds,v,e5)) &&
      (is_same_cell(e,invalid_cell) || compare_edges(func_img,flag_img,ds,e5,e)))
    e = e5;

  if(!is_same_cell(e,invalid_cell))
    mark_pair(func_img,flag_img,ds,v,e,flag_buf);
}

inline void set_pair_edge
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t e,
  const cell_t d0,
  const cell_t d1,
  __global flag_t * flag_buf)
{
  if(!contains(ds.wr,e))
    return;

  cell_t f0 = e - d0;
  cell_t f1 = e + d0;
  cell_t f2 = e - d1;
  cell_t f3 = e + d1;

  cell_t f = invalid_cell;

  if (is_pairable(func_img,flag_img,ds,e,f0))
    f = f0;

  if ((is_pairable(func_img,flag_img,ds,e,f1)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f1,f)))
    f = f1;

  if ((is_pairable(func_img,flag_img,ds,e,f2)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f2,f)))
    f = f2;

  if ((is_pairable(func_img,flag_img,ds,e,f3)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f3,f)))
    f = f3;

  if(!is_same_cell(f,invalid_cell))
    mark_pair(func_img,flag_img,ds,e,f,flag_buf);
}

inline void set_pair_face
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t f,
  const cell_t d,
  __global flag_t * flag_buf)
{
  if(!contains(ds.wr,f))
    return;

  cell_t c0 = f - d;
  cell_t c1 = f + d;

  cell_t c = invalid_cell;

  if (is_pairable(func_img,flag_img,ds,f,c0))
    c = c0;

  if ((is_pairable(func_img,flag_img,ds,f,c1)) &&
      (is_same_cell(c,invalid_cell) || compare_cubes(func_img,flag_img,ds,c1,c)))
    c = c1;

  if(!is_same_cell(c,invalid_cell))
    mark_pair(func_img,flag_img,ds,f,c,flag_buf);
}

__kernel void assign_pairs
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);

  int N = num_cells2(ds.e);

  for( int i = tid ; i < N; i += num_thds)
  {
    cell_t v = i_to_c2(ds.e,i);

    set_pair_vert(func_img,flag_img,ds,v,flag_buf);

    set_pair_edge(func_img,flag_img,ds,v+Xdir,Ydir,Zdir,flag_buf);
    set_pair_edge(func_img,flag_img,ds,v+Ydir,Zdir,Xdir,flag_buf);
    set_pair_edge(func_img,flag_img,ds,v+Zdir,Xdir,Ydir,flag_buf);

    set_pair_face(func_img,flag_img,ds,v+XYdir,Zdir,flag_buf);
    set_pair_face(func_img,flag_img,ds,v+YZdir,Xdir,flag_buf);
    set_pair_face(func_img,flag_img,ds,v+ZXdir,Ydir,flag_buf);
  }
}

inline bool is_pairable2
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  cell_t p, cell_t p_mf ,cell_t q)
{
  if(!contains(ds.e,q))
    return false;

  if(is_boundry(ds.d,p) != is_boundry(ds.d,q))
    return false;

  flag_t q_fg    = read_imageui(flag_img, flag_sampler, to_int4(q-ds.e.lc)).x;
  cell_t q_mf    = flag_to_mxfct(q,q_fg);
  flag_t q_mf_fg = read_imageui(flag_img, flag_sampler, to_int4(q_mf-ds.e.lc)).x;
  cell_t q_mf_mf = flag_to_mxfct(q_mf,q_mf_fg);

  return ((!is_same_cell(q_mf,p)) && is_same_cell(q_mf_mf,p_mf));
}

inline void set_pair_edge2
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t e,
  const cell_t d0,
  const cell_t d1,
  __global flag_t * flag_buf
  )
{
  if(!contains(ds.wr,e))
    return;

  flag_t e_fg = read_imageui(flag_img, flag_sampler, to_int4(e-ds.e.lc)).x;

  if(is_paired(e_fg))
    return;

  cell_t e_mf =  flag_to_mxfct(e,e_fg);

  cell_t f0 = e - d0;
  cell_t f1 = e + d0;
  cell_t f2 = e - d1;
  cell_t f3 = e + d1;

  cell_t f = invalid_cell;

  if (is_pairable2(func_img,flag_img,ds,e,e_mf,f0))
    f = f0;

  if ((is_pairable2(func_img,flag_img,ds,e,e_mf,f1)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f1,f)))
    f = f1;

  if ((is_pairable2(func_img,flag_img,ds,e,e_mf,f2)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f2,f)))
    f = f2;

  if ((is_pairable2(func_img,flag_img,ds,e,e_mf,f3)) &&
      (is_same_cell(f,invalid_cell) || compare_faces(func_img,flag_img,ds,f3,f)))
    f = f3;

  if(!is_same_cell(f,invalid_cell))
  {
    flag_t f_fg = read_imageui(flag_img, flag_sampler, to_int4(f-ds.e.lc)).x;

    if(!is_paired(f_fg))
      mark_pair(func_img,flag_img,ds,e,f,flag_buf);
  }
}

inline void set_pair_face2
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t f,
  const cell_t d,
  __global flag_t * flag_buf)
{
  if(!contains(ds.wr,f))
    return;

  flag_t f_fg = read_imageui(flag_img, flag_sampler, to_int4(f-ds.e.lc)).x;

  if(is_paired(f_fg))
    return;

  cell_t f_mf =  flag_to_mxfct(f,f_fg);

  cell_t c0 = f - d;
  cell_t c1 = f + d;
  cell_t c  = invalid_cell;


  if (is_pairable2(func_img,flag_img,ds,f,f_mf,c0))
    c = c0;

  if ((is_pairable2(func_img,flag_img,ds,f,f_mf,c1)) &&
      (is_same_cell(c,invalid_cell) || compare_cubes(func_img,flag_img,ds,c1,c)))
    c = c1;

  if(!is_same_cell(c,invalid_cell))
  {
    flag_t c_fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(!is_paired(c_fg))
      mark_pair(func_img,flag_img,ds,f,c,flag_buf);
  }
}


__kernel void assign_pairs2
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);

  int N = num_cells2(ds.e);

  for( int i = tid ; i < N ; i += num_thds)
  {
    cell_t v = i_to_c2(ds.e,i);

    set_pair_edge2(func_img,flag_img,ds,v+Xdir,Ydir,Zdir,flag_buf);
    set_pair_edge2(func_img,flag_img,ds,v+Ydir,Zdir,Xdir,flag_buf);
    set_pair_edge2(func_img,flag_img,ds,v+Zdir,Xdir,Ydir,flag_buf);

    set_pair_face2(func_img,flag_img,ds,v+XYdir,Zdir,flag_buf);
    set_pair_face2(func_img,flag_img,ds,v+YZdir,Xdir,flag_buf);
    set_pair_face2(func_img,flag_img,ds,v+ZXdir,Ydir,flag_buf);
  }
}

inline bool is_pairable3
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  cell_t p, cell_t p_mf,cell_t p_mf_mf,cell_t q)
{
  if(!contains(ds.e,q))
    return false;

  if(is_boundry(ds.d,p) != is_boundry(ds.d,q))
    return false;

  flag_t q_fg    = read_imageui(flag_img, flag_sampler, to_int4(q-ds.e.lc)).x;
  cell_t q_mf    = flag_to_mxfct(q,q_fg);

  flag_t q_mf_fg = read_imageui(flag_img, flag_sampler, to_int4(q_mf-ds.e.lc)).x;
  cell_t q_mf_mf = flag_to_mxfct(q_mf,q_mf_fg);

  flag_t q_mf_mf_fg = read_imageui(flag_img, flag_sampler, to_int4(q_mf_mf-ds.e.lc)).x;
  cell_t q_mf_mf_mf = flag_to_mxfct(q_mf_mf,q_mf_mf_fg);


  return (!is_same_cell(q_mf,p) && !is_same_cell(q_mf_mf,p_mf) && is_same_cell(q_mf_mf_mf,p_mf_mf));
}

inline void set_pair_face3
( __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  const dataset_t ds,
  const cell_t f,
  const cell_t d,
  __global flag_t * flag_buf)
{
  if(!contains(ds.wr,f))
    return;

  flag_t f_fg = read_imageui(flag_img, flag_sampler, to_int4(f-ds.e.lc)).x;
  cell_t f_mf =  flag_to_mxfct(f,f_fg);

  if(is_paired(f_fg))
    return;

  flag_t f_mf_fg = read_imageui(flag_img, flag_sampler, to_int4(f_mf-ds.e.lc)).x;
  cell_t f_mf_mf = flag_to_mxfct(f_mf,f_mf_fg);

  cell_t c0 = f - d;
  cell_t c1 = f + d;
  cell_t c  = invalid_cell;

  if (is_pairable3(func_img,flag_img,ds,f,f_mf,f_mf_mf,c0))
    c = c0;

  if ((is_pairable3(func_img,flag_img,ds,f,f_mf,f_mf_mf,c1)) &&
      (is_same_cell(c,invalid_cell) || compare_cubes(func_img,flag_img,ds,c1,c)))
    c = c1;

  if(!is_same_cell(c,invalid_cell))
  {
    flag_t c_fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(!is_paired(c_fg))
      mark_pair(func_img,flag_img,ds,f,c,flag_buf);
  }
}

__kernel void assign_pairs3
(
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  cell_t rct_lc,  cell_t rct_uc,
  cell_t ext_lc,  cell_t ext_uc,
  cell_t dom_lc,  cell_t dom_uc,
  __global flag_t   * flag_buf
)
{
  int tid      = get_global_id(0);
  int num_thds = get_global_size(0);

  dataset_t ds = make_dataset(rct_lc,rct_uc,ext_lc,ext_uc,dom_lc,dom_uc);

  int N = num_cells2(ds.e);

  for( int i = tid ; i < N ; i += num_thds)
  {
    cell_t v = i_to_c2(ds.e,i);

    set_pair_face3(func_img,flag_img,ds,v+XYdir,Zdir,flag_buf);
    set_pair_face3(func_img,flag_img,ds,v+YZdir,Xdir,flag_buf);
    set_pair_face3(func_img,flag_img,ds,v+ZXdir,Ydir,flag_buf);
  }
}
