/*!

	\file
	\brief Memory allocation functions
	
	Provides the functions to allocate a 1D, 2D, or 3D array as a contiguous block of memory. For 2d (3D) arrays, the right most index (i.e., ni here) runs over contiguous memory chunks, followed by nj, and nk (if defined).

*/

void *alloc1DArr(size_t elemsize, int ni);
void **alloc2DArr(size_t elemsize, int nj, int ni);
void ***alloc3DArr(size_t elemsize, int nk, int nj, int ni);
