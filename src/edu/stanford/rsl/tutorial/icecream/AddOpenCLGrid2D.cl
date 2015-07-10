#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable


__constant sampler_t linearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

// OpenCL kernel for adding two OpenCLGrid2D to each other

//kernel void AddOpenCLGrid2D(
//	global float* resultGrid,
//	global image2d_t grid1,
//	global image2d_t grid2,
//	int gridSizeX,
//	int gridSizeY) 
//{
//	const unsigned int x = get_global_id(0);// x index
//	const unsigned int y = get_global_id(1);// y index
	
	//check if inside image boundaries
//	if ( x>= gridSizeX || y >= gridSizeY)
//		return;

//	for (int idx = 0; idx <= gridSizeX; idx++)
//	{
//		for (int idy = 0; idy <= gridSizeY; idy++)
//		{
//			float2 xy = {idx, idy};
//			float grid1val = read_imagef(grid1, linearSampler, xy).x;
//			float grid2val = read_imagef(grid2, linearSampler, xy).x;
//			resultGrid[idx] = grid1val;
//		}
//	}
//	return;
//}

kernel void AddOpenCLGrid2D(
	__global float* resultGrid,
	__global float* grid1,
	__global float* grid2,
	int gridSizeX,
	int gridSizeY) 
{
	const unsigned int x = get_global_id(0);// x index
	const unsigned int y = get_global_id(1);// x index
	const unsigned int idx = y*gridSizeX + x;
	
	//check if inside image boundaries
	if ( x > gridSizeX*gridSizeY)
		return;

	float grid1val = grid1[idx];
	float grid2val = grid2[idx];
	resultGrid[idx] = grid1val+grid2val;
	
	return;
}