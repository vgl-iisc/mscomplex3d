#define NUM_LOCAL_THREADS (OPENCL_NUM_WORK_ITEMS_PER_GROUP)
#define NUM_GROUPS        (OPENCL_NUM_WORK_GROUPS)
#define binary_op(a, b) ((a) +(b))

inline void __scan_local_threads(__local int * array)
{
  barrier(CLK_LOCAL_MEM_FENCE);

  int val = array[get_local_id(0)];

  if (NUM_LOCAL_THREADS >   1) { if(get_local_id(0) >=   1) { int tmp = array[get_local_id(0) -   1]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >   2) { if(get_local_id(0) >=   2) { int tmp = array[get_local_id(0) -   2]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >   4) { if(get_local_id(0) >=   4) { int tmp = array[get_local_id(0) -   4]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >   8) { if(get_local_id(0) >=   8) { int tmp = array[get_local_id(0) -   8]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >  16) { if(get_local_id(0) >=  16) { int tmp = array[get_local_id(0) -  16]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >  32) { if(get_local_id(0) >=  32) { int tmp = array[get_local_id(0) -  32]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS >  64) { if(get_local_id(0) >=  64) { int tmp = array[get_local_id(0) -  64]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS > 128) { if(get_local_id(0) >= 128) { int tmp = array[get_local_id(0) - 128]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS > 256) { if(get_local_id(0) >= 256) { int tmp = array[get_local_id(0) - 256]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_LOCAL_THREADS > 512) { if(get_local_id(0) >= 512) { int tmp = array[get_local_id(0) - 512]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
}

inline void __scan_num_groups(__local int * array)
{
  barrier(CLK_LOCAL_MEM_FENCE);

  int val = array[get_local_id(0)];

  if (NUM_GROUPS >   1) { if(get_local_id(0) >=   1) { int tmp = array[get_local_id(0) -   1]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >   2) { if(get_local_id(0) >=   2) { int tmp = array[get_local_id(0) -   2]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >   4) { if(get_local_id(0) >=   4) { int tmp = array[get_local_id(0) -   4]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >   8) { if(get_local_id(0) >=   8) { int tmp = array[get_local_id(0) -   8]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >  16) { if(get_local_id(0) >=  16) { int tmp = array[get_local_id(0) -  16]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >  32) { if(get_local_id(0) >=  32) { int tmp = array[get_local_id(0) -  32]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS >  64) { if(get_local_id(0) >=  64) { int tmp = array[get_local_id(0) -  64]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS > 128) { if(get_local_id(0) >= 128) { int tmp = array[get_local_id(0) - 128]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS > 256) { if(get_local_id(0) >= 256) { int tmp = array[get_local_id(0) - 256]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
  if (NUM_GROUPS > 512) { if(get_local_id(0) >= 512) { int tmp = array[get_local_id(0) - 512]; val = binary_op(tmp, val); } barrier(CLK_LOCAL_MEM_FENCE); array[get_local_id(0)] = val; barrier(CLK_LOCAL_MEM_FENCE); }
}

__kernel void scan_local_sums(__global int  * array, __global int  * group_sums)
{
  __local int sarray[NUM_LOCAL_THREADS+1];

  sarray[NUM_LOCAL_THREADS] = 0;
  sarray[get_local_id(0)]   = array[get_global_id(0)];

  __scan_local_threads(sarray);

  int ridx = (get_local_id(0) -1 + NUM_LOCAL_THREADS + 1) % (NUM_LOCAL_THREADS+1);

  array[get_global_id(0)] = sarray[ridx];

  if(get_local_id(0) == 0)
    group_sums[get_group_id(0)] = sarray[NUM_LOCAL_THREADS-1];
}

__kernel void scan_group_sums(__global int  * group_sums)
{
  __local int sarray[NUM_GROUPS];

  sarray[get_local_id(0)] = group_sums[get_local_id(0)];

  __scan_num_groups(sarray);

  group_sums[get_local_id(0)] = sarray[get_local_id(0)];
}

__kernel void scan_update_sums(__global int  * array, __global int  * group_sums)
{
  if(get_group_id(0) != 0)
    array[get_global_id(0)] += group_sums[get_group_id(0)-1];
}


__kernel void mark_cps
(
  cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __global flag_t   * flag_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);

  int N = num_cells(ds.r);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(ds.r,i);

    flag_t fg = flag_buf[c_to_i(ds.e,c)];

    if(is_paired(fg))
      continue;

    flag_buf[c_to_i(ds.e,c)] |= CELLFLAG_CRITICAL;
  }
}

__kernel void mark_boundry_cps
(
  cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  cell_pair_t bnd_,
  cell_t      bnd_dir,
  __global flag_t   * flag_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t  bnd  = make_rect2(bnd_);
  int N        = num_cells(bnd);

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(bnd,i);

    flag_t fg = flag_buf[c_to_i(ds.e,c)];

    if(!is_paired(fg))
      continue;

    cell_t p = flag_to_pair(c,fg);

    if( (!is_same_cell(p+bnd_dir,c)) && (!is_same_cell(c+bnd_dir,p)))
      continue;

    flag_buf[c_to_i(ds.e,c)] = fg|CELLFLAG_CRITICAL;
    flag_buf[c_to_i(ds.e,p)] = pair_to_flag(p,c)|mxfct_to_flag(p,c)|CELLFLAG_CRITICAL;

    if(!contains(ds.r,p))
      continue;
  }
}

__kernel void count_cps
( cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  flag_img,
  __global int          *cp_count_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);

  int N = num_cells(ds.r);

  int n_cp = 0;

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(ds.r,i);

    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(is_paired(fg))
      continue;

    n_cp++;
  }

  cp_count_buf[get_global_id(0)] = n_cp;
}

__kernel void count_boundry_cps
(
  cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  cell_pair_t bnd_,
  cell_t      bnd_dir,
  __read_only image3d_t  flag_img,
  __global int    * cp_count_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t  bnd  = make_rect2(bnd_);
  int N        = num_cells(bnd);

  int n_cp = cp_count_buf[get_global_id(0)];

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(bnd,i);

    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(!is_paired(fg))
      continue;

    cell_t p = flag_to_pair(c,fg);

    if( (!is_same_cell(p+bnd_dir,c)) && (!is_same_cell(c+bnd_dir,p)))
      continue;

    n_cp +=2;
  }

  cp_count_buf[get_global_id(0)] = n_cp;
}

inline cell_t get_max_vert(const cell_t c, const dataset_t ds, __read_only image3d_t  flag_img)
{
  cell_t v = c;

  switch(cell_dim(c))
  {
    case 3: v = flag_to_mxfct(v,read_imageui(flag_img, flag_sampler, to_int4(v-ds.e.lc)).x);
    case 2: v = flag_to_mxfct(v,read_imageui(flag_img, flag_sampler, to_int4(v-ds.e.lc)).x);
    case 1: v = flag_to_mxfct(v,read_imageui(flag_img, flag_sampler, to_int4(v-ds.e.lc)).x);
  }
  return v;
}

__kernel void save_cps
(
  cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  __global int    * cp_offset_buf,
  __global short  * cp_cellid_buf,
  __global char   * cp_index_buf,
  __global int    * cp_pair_idx_buf,
  __global short  * cp_vertid_buf,
  __global func_t * cp_func_buf
)
{
  dataset_t ds = make_dataset2(rct,ext,dom);

  int N = num_cells(ds.r);

  int n_cp = cp_offset_buf[get_global_id(0)];

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(ds.r,i);

    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(is_paired(fg))
      continue;

    cell_t v     = get_max_vert(c,ds, flag_img);
    func_t  func = read_imagef(func_img, func_sampler, to_int4(v-ds.e.lc)/2).x;

    cp_cellid_buf[3*n_cp + 0] = c.x;
    cp_cellid_buf[3*n_cp + 1] = c.y;
    cp_cellid_buf[3*n_cp + 2] = c.z;

    cp_vertid_buf[3*n_cp + 0] = v.x;
    cp_vertid_buf[3*n_cp + 1] = v.y;
    cp_vertid_buf[3*n_cp + 2] = v.z;

    cp_index_buf[n_cp]        = cell_dim(c);
    cp_pair_idx_buf[n_cp]     = -1;
    cp_func_buf[n_cp]         = func;

    n_cp++;
  }

  cp_offset_buf[get_global_id(0)] = n_cp;
}

__kernel void save_boundry_cps
(
  cell_pair_t rct,
  cell_pair_t ext,
  cell_pair_t dom,
  cell_pair_t bnd_,
  cell_t      bnd_dir,
  __read_only image3d_t  func_img,
  __read_only image3d_t  flag_img,
  __global int    * cp_offset_buf,
  __global short  * cp_cellid_buf,
  __global char   * cp_index_buf,
  __global int    * cp_pair_idx_buf,
  __global short  * cp_vertid_buf,
  __global func_t * cp_func_buf
)
{

  dataset_t ds = make_dataset2(rct,ext,dom);
  rect_t  bnd  = make_rect2(bnd_);
  int N        = num_cells(bnd);

  int n_cp = cp_offset_buf[get_global_id(0)];

  for( int i = get_global_id(0) ; i < N; i += get_global_size(0))
  {
    cell_t c = i_to_c(bnd,i);

    flag_t fg = read_imageui(flag_img, flag_sampler, to_int4(c-ds.e.lc)).x;

    if(!is_paired(fg))
      continue;

    cell_t p = flag_to_pair(c,fg);

    if( (!is_same_cell(p+bnd_dir,c)) && (!is_same_cell(c+bnd_dir,p)))
      continue;

    cell_t v     = get_max_vert(c,ds, flag_img);
    func_t  func = read_imagef(func_img, func_sampler, to_int4(v-ds.e.lc)/2).x;

    cp_cellid_buf[3*n_cp + 0] = c.x;
    cp_cellid_buf[3*n_cp + 1] = c.y;
    cp_cellid_buf[3*n_cp + 2] = c.z;

    cp_vertid_buf[3*n_cp + 0] = v.x;
    cp_vertid_buf[3*n_cp + 1] = v.y;
    cp_vertid_buf[3*n_cp + 2] = v.z;

    cp_index_buf[n_cp]        = cell_dim(c);
    cp_pair_idx_buf[n_cp]     = n_cp+1;
    cp_func_buf[n_cp]         = func;

    n_cp++;

    cp_cellid_buf[3*n_cp + 0] = p.x;
    cp_cellid_buf[3*n_cp + 1] = p.y;
    cp_cellid_buf[3*n_cp + 2] = p.z;

    cp_vertid_buf[3*n_cp + 0] = v.x;
    cp_vertid_buf[3*n_cp + 1] = v.y;
    cp_vertid_buf[3*n_cp + 2] = v.z;

    cp_index_buf[n_cp]        = cell_dim(p);
    cp_pair_idx_buf[n_cp]     = n_cp-1;
    cp_func_buf[n_cp]         = func;

    n_cp++;
  }

  cp_offset_buf[get_global_id(0)] = n_cp;
}

