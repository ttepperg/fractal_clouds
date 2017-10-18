/*
	 Provides the functions to allocate a 1D, 2D, or 3D array as a contiguous block of memory. For 2d (3D) arrays, the right most index (i.e., ni here) runs over contiguous memory chunks, followed by nj, and nk (if defined).

*/

#include <stdlib.h>

void *alloc1DArr(size_t elemsize, int ni)
{	/*!
	\fn
	returns a pointer to a contiguos block memory of ni elements;
	size in bytes is ni * elemsize
*/
	void *auxPtr = NULL;
	
	// allocate main block of memory
	auxPtr = (void *)malloc(ni * elemsize);

	return (auxPtr);

} // end of alloc1DArr

void **alloc2DArr(size_t elemsize, int nj, int ni)
{	/*!
	\fn
	returns a pointer to a pointer to a contiguos block memory of nj * ni elements;
	size in bytes is nj * ni * elemsize
*/

	void **auxPPtr = NULL, *auxPtr = NULL;
	int j;
	
	// allocate array of nj pointers
	auxPPtr = (void *)malloc(nj * sizeof(void *));

	if ( auxPPtr == NULL )
		return (auxPPtr);

	// allocate main block of memory
	auxPtr = (void *)malloc(nj * ni * elemsize);

	if ( auxPtr == NULL )
	{
		free(auxPPtr);
		return (NULL);
	}

	// assign row addresses to pointers (uses pointer arithmetic)
	for (j = 0; j < nj; j++)
		auxPPtr[j] = (void *) (auxPtr + j * ni * elemsize);

	return (auxPPtr);

} // end of alloc2DArr

void ***alloc3DArr(size_t elemsize, int nk, int nj, int ni)
{	/*!
	\fn
	returns a pointer to a pointer to a pointer to a contiguos block memory of 
	nk * nj * ni elements; size in bytes is nk * nj * ni * elemsize
*/
	void ***auxPPPtr = NULL, *auxPtr = NULL;
	int k, j;

	// allocate 2D array of pointers
	auxPPPtr = (void ***)alloc2DArr(sizeof(void *), nk, nj);

	if (auxPPPtr == NULL)
		return (auxPPPtr);

// allocate main block of memory
	auxPtr = (void *)malloc(nk * nj * ni * elemsize);

	if (auxPtr == NULL)
	{
		free(auxPPPtr);
		return (auxPtr);
	}

// assign address of first element of each i-row to pointers in 2D array
	
	for(k = 0; k < nk; k++)
		for(j = 0; j < nj; j++)
			auxPPPtr[k][j] = (void *)(auxPtr + (k * nj * ni + j * ni) * elemsize);

		return (auxPPPtr);

} // end of alloc3DArr
