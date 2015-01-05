#include <cutil.h>

texture<flag_t, 3, cudaReadModeElementType>     flag_texture;
texture<func_t, 3, cudaReadModeElementType>     func_texture;

__device__ bool compare_verts(cell_t v1, cell_t v2)
{
  func_t f1 = tex3D(func_texture, v1.x/2, v1.y/2, v1.z/2);
  func_t f2 = tex3D(func_texture, v2.x/2, v2.y/2, v2.z/2);

  if( f1 != f2)
    return f1 < f2;

  return v1 < v2;
}

__device__ bool compare_edges(cell_t e1, cell_t e2)
{
  flag_t f1 = tex3D(flag_texture, e1.x, e1.y, e1.z);
  flag_t f2 = tex3D(flag_texture, e2.x, e2.y, e2.z);

  cell_t v1 = flag_to_mxfct(e1,f1);
  cell_t v2 = flag_to_mxfct(e2,f2);

  if( v1 == v2)
  {
    v1 = second_max_facet(e1,v1);
    v2 = second_max_facet(e2,v2);
  }

  return compare_verts(v1,v2);
}

__device__ bool compare_faces(cell_t fc1, cell_t fc2)
{
  flag_t f1 = tex3D(flag_texture, fc1.x, fc1.y, fc1.z);
  flag_t f2 = tex3D(flag_texture, fc2.x, fc2.y, fc2.z);

  cell_t e1 = flag_to_mxfct(fc1,f1);
  cell_t e2 = flag_to_mxfct(fc1,f2);

  if( e1 == e2)
  {
    e1 = second_max_facet(fc1,e1);
    e2 = second_max_facet(fc2,e2);
  }

  return compare_edges(e1,e2);
}

__device__ bool compare_cubes(cell_t c1, cell_t c2)
{
  flag_t f1 = tex3D(flag_texture, c1.x, c1.y, c1.z);
  flag_t f2 = tex3D(flag_texture, c2.x, c2.y, c2.z);

  cell_t fc1 = flag_to_mxfct(c1,f1);
  cell_t fc2 = flag_to_mxfct(c1,f2);

  if( fc1 == fc2)
  {
    fc1 = second_max_facet(c1,fc1);
    fc2 = second_max_facet(c2,fc2);
  }

  return compare_faces(fc1,fc2);
}

__global__  void assign_maxfacet_vert
    (const rect_t ext_rect,
     flag_t *cell_flags)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int num_thds= gridDim.x*blockDim.x;

  rect_t cell_rect(ext_rect.lc,ext_rect.uc);

  for( int i = tid ; i < cell_rect.num_cells2(); i += num_thds)
  {
    cell_t c = cell_rect.i_to_c2(i);

    cell_flags[ext_rect.c_to_i(c)] = 0;
  }
}


__global__  void assign_maxfacet_edge
    (const rect_t ext_rect,
     const cell_t edir,
     flag_t *cell_flags)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int num_thds= gridDim.x*blockDim.x;

  rect_t cell_rect(ext_rect.lc+edir,ext_rect.uc-edir);

  for( int i = tid ; i < cell_rect.num_cells2(); i += num_thds)
  {
    cell_t c = cell_rect.i_to_c2(i);

    cell_t v1 = c - edir;
    cell_t v2 = c + edir;

    cell_t v = v1;

    if( compare_verts(v,v2))
      v = v2;

    cell_flags[ext_rect.c_to_i(c)] = mxfct_to_flag(c,v);
  }
}

__global__ void assign_maxfacet_face
    (const rect_t ext_rect,
     const cell_t fdir1,
     const cell_t fdir2,
     flag_t *cell_flags)

{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int num_thds= gridDim.x*blockDim.x;

  cell_t fdir = fdir1+fdir2;
  rect_t cell_rect(ext_rect.lc+fdir,ext_rect.uc-fdir);

  for( int i = tid ; i < cell_rect.num_cells2(); i += num_thds)
  {
    cell_t c = cell_rect.i_to_c2(i);

    cell_t e1 = c-fdir1;
    cell_t e2 = c+fdir1;
    cell_t e3 = c-fdir2;
    cell_t e4 = c+fdir2;

    cell_t e = e1;

    if( compare_edges(e,e2))
      e = e2;

    if( compare_edges(e,e3))
      e = e3;

    if( compare_edges(e,e4))
      e = e4;

    cell_flags[ext_rect.c_to_i(c)] = mxfct_to_flag(c,e);
  }
}

__global__  void assign_maxfacet_cube
    (const rect_t ext_rect,
     flag_t *cell_flags)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int num_thds= gridDim.x*blockDim.x;

  rect_t cell_rect(ext_rect.lc+1,ext_rect.uc-1);

  for( int i = tid ; i < cell_rect.num_cells2(); i += num_thds)
  {
    cell_t c = cell_rect.i_to_c2(i);

    cell_t f1 = c-mk_cell(1,0,0);
    cell_t f2 = c+mk_cell(1,0,0);
    cell_t f3 = c-mk_cell(0,1,0);
    cell_t f4 = c+mk_cell(0,1,0);
    cell_t f5 = c-mk_cell(0,0,1);
    cell_t f6 = c+mk_cell(0,0,1);

    cell_t f = f1;

    if( compare_faces(f,f2))
      f = f2;

    if( compare_faces(f,f3))
      f = f3;

    if( compare_faces(f,f4))
      f = f4;

    if( compare_faces(f,f5))
      f = f5;

    if( compare_faces(f,f6))
      f = f6;

    cell_flags[ext_rect.c_to_i(c)] = mxfct_to_flag(c,f);
  }
}


extern "C"
void assign_max_facet_cuda(short* erptr,func_t *func , flag_t *flag)
{
  rect_t ext_rect(erptr[0],erptr[1],erptr[2],erptr[3],erptr[4],erptr[5]);

  cudaArray *d_func = 0;
  cell_t func_span = (ext_rect.uc - ext_rect.lc)/2 + 1;

  cudaExtent func_vol_extent = make_cudaExtent(func_span.x * sizeof(func_t),func_span.y,func_span.z);

  // create 3D array
  cudaChannelFormatDesc func_channel_desc = cudaCreateChannelDesc<func_t>();
  cutilSafeCall(cudaMalloc3DArray(&d_func, &func_channel_desc, func_vol_extent));

  // pitched pointer
  cudaPitchedPtr func_pitched_ptr;
  func_pitched_ptr.ptr      = func;
  func_pitched_ptr.pitch    = func_vol_extent.width*sizeof(func_t);
  func_pitched_ptr.xsize    = func_vol_extent.width;
  func_pitched_ptr.ysize    = func_vol_extent.height;

  // copy data to 3D array
  cudaMemcpy3DParms func_copy_params = {0};
  func_copy_params.srcPtr   = func_pitched_ptr;
  func_copy_params.dstArray = d_func;
  func_copy_params.extent   = func_vol_extent;
  func_copy_params.kind     = cudaMemcpyHostToDevice;
  cutilSafeCall(cudaMemcpy3D(&func_copy_params));

  // set texture parameters
  func_texture.normalized = false;                      // access with normalized texture coordinates
  func_texture.filterMode = cudaFilterModePoint;      // linear interpolation
  func_texture.addressMode[0] = cudaAddressModeClamp;  // wrap texture coordinates
  func_texture.addressMode[1] = cudaAddressModeClamp;

  // bind array to 3D texture
  cutilSafeCall(cudaBindTextureToArray(func_texture, d_func, func_channel_desc));



  flag_t *d_flag = 0;
  int     d_flag_size = ext_rect.num_cells()*sizeof(flag_t);

  cutilSafeCall(cudaMalloc(&d_flag, d_flag_size));


  flag_vol_extent.width  = flag_span.x;
  flag_vol_extent.height = flag_span.y;
  flag_vol_extent.depth  = flag_span.z;

  // create 3D array
  cudaChannelFormatDesc flag_channel_desc = cudaCreateChannelDesc<flag_t>();


  // set texture parameters
  flag_texture.normalized = false;                      // access with normalized texture coordinates
  flag_texture.filterMode = cudaFilterModePoint;      // linear interpolation
  flag_texture.addressMode[0] = cudaAddressModeClamp;  // wrap texture coordinates
  flag_texture.addressMode[1] = cudaAddressModeClamp;

  // bind array to 3D texture
  cutilSafeCall(cudaBindTexture(flag_texture, d_flag, flag_channel_desc));



  dim3 block(256,1);
  dim3 grid(32,1);

  assign_maxfacet_vert<<<grid,block>>>(ext_rect,d_flag);

  assign_maxfacet_edge<<<grid,block>>>(ext_rect,mk_cell(1,0,0),d_flag);
  assign_maxfacet_edge<<<grid,block>>>(ext_rect,mk_cell(0,1,0),d_flag);
  assign_maxfacet_edge<<<grid,block>>>(ext_rect,mk_cell(0,0,1),d_flag);

  assign_maxfacet_face<<<grid,block>>>(ext_rect,mk_cell(1,0,0),mk_cell(0,1,0),d_flag);
  assign_maxfacet_face<<<grid,block>>>(ext_rect,mk_cell(0,1,0),mk_cell(0,0,1),d_flag);
  assign_maxfacet_face<<<grid,block>>>(ext_rect,mk_cell(0,0,1),mk_cell(1,0,0),d_flag);

  assign_maxfacet_cube<<<grid,block>>>(ext_rect,d_flag);



  // pitched pointer
  cudaPitchedPtr flag_pitched_ptr;
  flag_pitched_ptr.ptr      = flag;
  flag_pitched_ptr.pitch    = flag_vol_extent.width*sizeof(flag_t);
  flag_pitched_ptr.xsize    = flag_vol_extent.width;
  flag_pitched_ptr.ysize    = flag_vol_extent.height;

  // copy data to 3D array
  cudaMemcpy3DParms flag_copy_params = {0};
  flag_copy_params.srcPtr   = flag_pitched_ptr;
  flag_copy_params.dstArray = d_flag;
  flag_copy_params.extent   = flag_vol_extent;
  flag_copy_params.kind     = cudaMemcpyDeviceToHost;
  cutilSafeCall(cudaMemcpy3D(&flag_copy_params));

  cudaFreeArray(d_flag);
  cudaFreeArray(d_func);
}

#define MAX(a,b) ((a > b) ? a : b)

#include <cufft.h>
#include <cuda_runtime.h>


// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores(int major, int minor)
{
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct {
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] =
	{ { 0x10,  8 },
	  { 0x11,  8 },
	  { 0x12,  8 },
	  { 0x13,  8 },
	  { 0x20, 32 },
	  { 0x21, 48 },
	  {   -1, -1 }
	};

	int index = 0;
	while (nGpuArchCoresPerSM[index].SM != -1) {
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor) ) {
			return nGpuArchCoresPerSM[index].Cores;
		}
		index++;
	}
	printf("MapSMtoCores undefined SMversion %d.%d!\n", major, minor);
	return -1;
}
// end of GPU Architecture definitions


// This function returns the best GPU (with maximum GFLOPS)
inline int cutGetMaxGflopsDeviceId()
{
	int current_device   = 0, sm_per_multiproc = 0;
	int max_compute_perf = 0, max_perf_device  = 0;
	int device_count     = 0, best_SM_arch     = 0;
	cudaDeviceProp deviceProp;

	cudaGetDeviceCount( &device_count );
	// Find the best major SM Architecture GPU device
	while ( current_device < device_count ) {
		cudaGetDeviceProperties( &deviceProp, current_device );
		if (deviceProp.major > 0 && deviceProp.major < 9999) {
			best_SM_arch = MAX(best_SM_arch, deviceProp.major);
		}
		current_device++;
	}

    // Find the best CUDA capable GPU device
	current_device = 0;
	while( current_device < device_count ) {
		cudaGetDeviceProperties( &deviceProp, current_device );
		if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
		    sm_per_multiproc = 1;
		} else {
			sm_per_multiproc = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
		}

		int compute_perf  = deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;
		if( compute_perf  > max_compute_perf ) {
            // If we find GPU with SM major > 2, search only these
			if ( best_SM_arch > 2 ) {
				// If our device==dest_SM_arch, choose this, or else pass
				if (deviceProp.major == best_SM_arch) {
					max_compute_perf  = compute_perf;
					max_perf_device   = current_device;
				}
			} else {
				max_compute_perf  = compute_perf;
				max_perf_device   = current_device;
			}
		}
		++current_device;
	}
	return max_perf_device;
}

extern "C"
void init_cuda()
{
  cudaSetDevice( cutGetMaxGflopsDeviceId() );
}


