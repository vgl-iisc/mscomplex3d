__kernel void init_propagate
  ( cell_pair_t rct,
    cell_pair_t ext,
    cell_pair_t dom,
    rect_t ex_rect,
    __read_only image3d_t flag_img,
    __global int * ex_own_buf
  )
{
  int N = num_cells2(ex_rect);

  dataset_t ds = make_dataset2(rct,ext,dom);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c  = i_to_c2(ex_rect,i);
    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    cell_t o = c;

    if(!is_critical(fg))
    {
      cell_t p = flag_to_pair(c,fg);
      o  = p + p-c;
    }

    ex_own_buf[c_to_i2(ex_rect,c)] = c_to_i2(ex_rect,o);
  }
}

__kernel void propagate
(   cell_pair_t rct,
    cell_pair_t ext,
    cell_pair_t dom,
    rect_t ex_rect,
    __global int * ex_own_buf1,
    __global int * ex_own_buf2,
    __global int * is_updated
  )
{
  int N = num_cells2(ex_rect);

  int updated = 0;

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    int  o = ex_own_buf1[i];
    int oo = ex_own_buf1[o];
    ex_own_buf2[i] = oo;

    if(o != oo)
      updated = 1;
  }

  if(updated == 1)
    *is_updated = 1;
}

__kernel void update_cp_cell_to_cp_no
  ( cell_pair_t rct,
    cell_pair_t ext,
    cell_pair_t dom,
    rect_t ex_rect,
    __global short  *cp_cellid_buf,
    int num_cps,
    __global int    *ex_own_buf)
{
  int ex_dim = cell_dim(ex_rect.lc);

  dataset_t ds = make_dataset2(rct,ext,dom);

  for( int i = get_global_id(0) ; i < num_cps; i += get_global_size(0))
  {
    cell_t c;

    c.x = cp_cellid_buf[3*i+0];
    c.y = cp_cellid_buf[3*i+1];
    c.z = cp_cellid_buf[3*i+2];

    if(cell_dim(c) != ex_dim)
      continue;

    if(!contains(ds.r,c))
      continue;

    ex_own_buf[c_to_i2(ex_rect,c)] = i;
  }
}


__kernel void init_update_to_surv_cp_no
  ( cell_pair_t rct,
    cell_pair_t ext,
    cell_pair_t dom,
    rect_t ex_rect,
    __global short  *cp_cellid_buf,
    __global int    *surv_cp_no_buf,
    int num_cps,
    __global int    *ex_own_buf
  )

{
  int ex_dim = cell_dim(ex_rect.lc);

  dataset_t ds = make_dataset2(rct,ext,dom);

  for( int i = get_global_id(0) ; i < num_cps; i += get_global_size(0))
  {
    cell_t c;

    c.x = cp_cellid_buf[3*i+0];
    c.y = cp_cellid_buf[3*i+1];
    c.z = cp_cellid_buf[3*i+2];

    if(cell_dim(c) != ex_dim)
      continue;

    if(!contains(ds.r,c))
      continue;

    ex_own_buf[c_to_i2(ex_rect,c)] = -surv_cp_no_buf[i];
  }
}

__kernel void update_to_surv_cp_no
  ( cell_pair_t rct,
    cell_pair_t ext,
    cell_pair_t dom,
    rect_t ex_rect,
    __global int    *ex_own_buf1,
    __global int    *ex_own_buf2)
{
  int N = num_cells2(ex_rect);
  dataset_t ds = make_dataset2(rct,ext,dom);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c  = i_to_c2(ex_rect,i);

    int o  = ex_own_buf1[i];
    int oo = (o>=0)?(-ex_own_buf1[o]):(-o);

    ex_own_buf2[i] = oo;
  }
}

//__kernel void update_to_surviving_cp_no
//  ( cell_pair_t rct,
//    cell_pair_t ext,
//    cell_pair_t dom,
//    rect_t ex_rect,
//    __global int *ex_own_buf,
//    __global int *surv_cp_no,
//    int n_cps)
//{
//  int N = num_cells2(ex_rect);
//  dataset_t ds = make_dataset2(rct,ext,dom);

//  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
//  {
//    int o  = ex_own_buf1[i];
//    int oo = surv_cp_no[o];
//    ex_own_buf2[i] = oo;
//  }
//}
