#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>

#include <aabb.h>
#include <map>
#include <vector>
#include <array>
#include <iostream>
#include <stdexcept>
#include <algorithm>

/// <summary>
/// This class was made to replace the Boost multi_array
///	Some outcomes may vary, though the computation of the morse smale complex, and critical points shall not change.
///	Observed changed outcomes include
///	1. Change in array index for critical points in the critical point array
///	2. Runtime is longer than the previous implementation
/// </summary>
/// <typeparam name="T"></typeparam>
template <typename T, std::size_t N = 3>
class Array3D {
private:
	size_t size_z, size_y, size_x;

	std::array<short, N> reindex_values;

	// Compute linear index from 3D indices (Fortran order)
	size_t index(size_t x, size_t y, size_t z) const {
		if (z >= size_z || y >= size_y || x >= size_x) {
			throw std::out_of_range("Index out of bounds");
		}
		//return (x+reindex_values[0]) * size_z * size_y + (y+reindex_values[1]) * size_z + (z+reindex_values[2]);
		//return x * size_z * size_y + y * size_z + z;
		//if(x==256 && y==4 && z==1)
			//std::cout << "\nWe want index: " << z * size_x * size_y + y * size_x + x;;
		return z * size_x * size_y + y * size_x + x;
	}

public:
	std::vector<T> data; // Flattened 3D data

	// Constructor: Initialize dimensions and fill with default value
	Array3D(size_t x, size_t y, size_t z, T default_value = T())
		: size_z(z), size_y(y), size_x(x), data(z* y* x, default_value) {}

	Array3D(n_vector_t<T, N> n) : Array3D(n[0], n[1], n[2]) // Delegate to the main constructor
	{}


	/**
	 * \brief Reindexes the array. This function is no longer used, and is unsupported. It is not intended to perform anything.
	 * \param lower_corner 
	 */
	void reindex(const std::array<short, N>& lower_corner) {

		for (int i = 0; i < N; ++i)
		{
			if (lower_corner[i] != 0)
				throw std::runtime_error("Array reindexing with non-zero lower_corner support was removed");

		}

		//// Create a temporary copy of the original data
		//std::vector<T> new_data(data.size(), T());

		//// Iterate over the original array and assign new indices in the temporary data
		//for (size_t i = 0; i < size_x; ++i) {
		//    for (size_t j = 0; j < size_y; ++j) {
		//        for (size_t k = 0; k < size_z; ++k) {
		//            // Remap the indices based on new_indices
		//            int new_i = static_cast<int>(i) + new_indices[2];
		//            int new_j = static_cast<int>(j) + new_indices[1];
		//            int new_k = static_cast<int>(k) + new_indices[0];

		//            reindex_values[0] = new_indices[0];
		//            reindex_values[1] = new_indices[1];
		//            reindex_values[2] = new_indices[2];

		//            // Ensure that the new indices are within bounds
		//            //if (new_i >= 0 && new_i < static_cast<int>(size_z) &&
		//            //    new_j >= 0 && new_j < static_cast<int>(size_y) &&
		//            //    new_k >= 0 && new_k < static_cast<int>(size_x)) {
		//                // Compute indices and transfer the value
		//                new_data[index(new_k, new_j, new_i)] = data[index(k, j, i)];
		//            //}
		//        }
		//    }
		//}

		//// Replace the original data with the reindexed data
		//data = std::move(new_data);


	}

	// Get raw pointer to data (non-const and const versions)
	T* getData() { return data.data(); }
	const T* getData() const { return data.data(); }

	// Access element (read/write)
	T& operator()(size_t x, size_t y, size_t z) {
		return data[index(x, y, z)];
	}

	// Access element (read-only)
	const T& operator()(size_t x, size_t y, size_t z) const {
		return data[index(x, y, z)];
	}

	// Access element (read/write)
	T& operator()(n_vector_t<size_t, N> n) {
		return data[index(n[0], n[1], n[2])];
	}

	// Access element (read-only)
	const T& operator()(n_vector_t<size_t, N> n) const {
		return data[index(n[0], n[1], n[2])];
	}


	// Resize the array and fill with default value
	void resize(size_t x, size_t y, size_t z, T default_value = T())
	{
		size_z = z;
		size_y = y;
		size_x = x;
		data.resize(0);
		data.resize(z * y * x,default_value);
	}

	void resize(n_vector_t<T, N> n, T default_value = T()) {
		resize(n[0], n[1], n[2]);
	}

	// Fill all elements with a specific value 
	void fill(T value) {
		std::fill(data.begin(), data.end(), value);
	}

	// Debug: Print all elements
	void print() const {
		for (size_t k = 0; k < size_z; ++k) {
			for (size_t j = 0; j < size_y; ++j) {
				for (size_t i = 0; i < size_x; ++i) {
					std::cout << (*this)(i, j, k) << " ";
				}
				std::cout << "\n";
			}
			std::cout << "----\n";
		}
	}
};




namespace grid
{

	typedef uint32_t uint;
	const uint gc_grid_dim = 3;

	typedef int16_t                                         cell_coord_t;

	typedef uint8_t                                        cell_flag_t;

	typedef float                                           cell_fn_t;
	typedef std::shared_ptr<cell_fn_t >                   cell_fn_ptr_t;

	typedef aabb::aabb_t<cell_coord_t, gc_grid_dim>          rect_t;
	typedef aabb::aabb_t<cell_coord_t, gc_grid_dim>::point_t cellid_t;
	typedef aabb::aabb_t<cell_coord_t, gc_grid_dim>::point_t rect_point_t;
	typedef aabb::aabb_t<cell_coord_t, gc_grid_dim>::point_t rect_size_t;
	typedef aabb::aabb_t<cell_coord_t, gc_grid_dim>::range_t rect_range_t;
	typedef std::vector<cellid_t>                           cellid_list_t;
	typedef std::vector<int>                                int_list_t;
	typedef std::vector<int8_t>                             int8_list_t;
	typedef std::vector<cell_fn_t>                          cell_fn_list_t;
	typedef std::vector<int8_t>                             bool_list_t;
	typedef std::vector<rect_t>                             rect_list_t;

	typedef std::shared_ptr<int_list_t>                   int_list_ptr_t;
	typedef std::shared_ptr<cellid_list_t>                cellid_list_ptr_t;

	//typedef boost::multi_array<int,gc_grid_dim>             int_marray_t;
	typedef Array3D<int>                                      int_marray_t;

	typedef std::pair<cellid_t, int>                         cellid_int_pair_t;
	typedef std::vector<cellid_int_pair_t>                  cellid_int_pair_list_t;

	typedef cellid_list_t                                   mfold_t;
	typedef std::vector<mfold_t>                            mfold_list_t;

	typedef n_vector_t<int, 2>                               int_pair_t;
	typedef std::vector<int_pair_t>                         int_pair_list_t;

	typedef std::pair<int, int>                              int_int_t;
	typedef std::vector<int_int_t>                          int_int_list_t;
	typedef std::map<int, int>                               int_int_map_t;

	typedef std::map<int, int>                               conn_t;
	typedef std::vector<conn_t>                             conn_list_t;

	enum eGDIR { DES = 0, ASC, GDIR_CT };

	const eGDIR GDIR_DES = DES;
	const eGDIR GDIR_ASC = ASC;

	/// \brief cell complex type
	enum eCCTYPE
	{
		CC_NON = 0,
		CC_PRIM = 1, // primal
		CC_DUAL = 2, // Dual
		CC_BOTH = 3
	};

	class dataset_t;
	class mscomplex_t;
	class data_manager_t;

	typedef std::shared_ptr<dataset_t>            dataset_ptr_t;
	typedef std::shared_ptr<mscomplex_t>          mscomplex_ptr_t;
	typedef std::shared_ptr<data_manager_t>       data_manager_ptr_t;

	inline int c_to_i(const rect_t& r, cellid_t c)
	{
		cellid_t s = r.span() + 1;
		c = (c - r.lc());
		return (s[0] * s[1] * c[2] + s[0] * c[1] + c[0]);
	}

	inline cellid_t i_to_c(const rect_t& r, int i)
	{
		cellid_t s = r.span() + 1;
		cellid_t c = r.lc() + (cellid_t(i % s[0], (i % (s[0] * s[1])) / s[0], i / (s[0] * s[1])));
		ASSERT(r.contains(c))
			return c;
	}

	inline int num_cells(const rect_t& r)
	{
		return c_to_i(r, r.uc()) + 1;
	}

	extern "C"
		utl::timer g_timer;

	inline int get_cell_dim(cellid_t c)
	{
		return (c[0] & 0x01) + (c[1] & 0x01) + (c[2] & 0x01);
	}

}

namespace grid {
	namespace opencl {

		/// \brief Init OpenCL runtime
		void init(int device=0);

		/// \brief Return the platform/device info after init
		std::string get_info();
	}

	namespace openmp {
		/// \brief Return the openmp info
		std::string get_info();
	}

	/**
	 * \brief Initializes OpenCL as well as provide a description of device used by OpenCL for parallel computation.
	 * \param device the device to be used for parallel computation. 0 for GPU. 1 is for CPU. OpenCL configuration should be set up before utilizing either. Is 0 by default.
	 * \return 
	 */
	inline std::string get_hw_info(int device=0)
	{
		opencl::init(device);
		return opencl::get_info() + openmp::get_info();
	}
}

#endif
