#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

// OpenCL kernel for adding two OpenCLGrid2D to each other

kernel void AddOpenCLGrid2D(
	global float* resultGrid,
	global float* grid1,
	global float* grid2,
	int gridSizeX,
	int gridSizeY
	) {
	
//	const unsigned int x = get_global_id(0);// x index
//	const unsigned int y = get_global_id(1);// y index
	
	// check if inside image boundaries
//	if ( x>= gridSizeX || y >= gridSizeY)
//		return;
	
//	for (int idx = 0; idx < x; idx++)
//	{
//		resultGrid[idx] = grid1[idx] + grid2[idx];
//	}
	
	return;
}