#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable


__constant sampler_t linearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

//double[] indexToPhysical(double i, double j, float pixelSpacingRecon[], double origin []) {
//	return new double[] {
//		i * pixelSpacingRecon[0]+ origin[0],
//		j * pixelSpacingRecon[1] + origin[1]
//	};
//}

// OpenCL kernel for backprojecting a sinogram
kernel void OpenCL_BP(
	__global float* backProjectionPic,
	__global float* sinogram,
	int numberProj,
	float detectorSpacing,
	int numberDetPixel,
	int sizeReconPic,
	float pixelSpacingReconX,
	float pixelSpacingReconY) 
{
	const unsigned int x = get_global_id(0);// x index
	const unsigned int y = get_global_id(1);// y index
	const unsigned int idx = y*256 + x;
	//check if inside image boundaries
	if ( x > sizeReconPic || y > sizeReconPic)
		return;
	
	for (int theta = 0; theta < numberProj; theta++)
	{
//		double alpha =  ((2*Math.PI/(numberProj))*theta);
//		double[] id = indexToPhysical(x, y, pixelSpacingRecon, backProjection.getOrigin());
//		double s=id[0]*Math.cos(alpha)+ id[1]*Math.sin(alpha);
//		s = (s- sinogram.getOrigin()[1])/detectorSpacing;
//		float value = InterpolationOperators.interpolateLinear(sinogram, theta, s);
//		value = (float) (Math.PI*value/numberProj);
		float value = 100;
		
		backProjectionPic[idx] = value;
	}
	
	return;
}



//	public Grid2D backProjection(int numberProj, float detectorSpacing, int numberDetPixel, Grid2D sinogram, int sizeRecon, float pixelSpacingRecon[]){
//		
//		Grid2D backProjection = new Grid2D(sizeRecon,  sizeRecon);
//		backProjection.setOrigin(-(sizeRecon-1)*pixelSpacingRecon[0]/2, -(sizeRecon-1)*pixelSpacingRecon[1]/2 );

//			for (int idxX =0; idxX< backProjection.getWidth(); idxX ++){
//				for (int idxY =0; idxY< backProjection.getHeight(); idxY ++){
//					for (int theta = 0; theta < numberProj; theta++)
//					{
//						double alpha =  ((2*Math.PI/(numberProj))*theta);
//						double[] id = indexToPhysical(idxX, idxY, pixelSpacingRecon, backProjection.getOrigin());
//						double s=id[0]*Math.cos(alpha)+ id[1]*Math.sin(alpha);
//						s = (s- sinogram.getOrigin()[1])/detectorSpacing ;
//						float value = InterpolationOperators.interpolateLinear(sinogram, theta, s);
//						value = (float) (Math.PI*value/numberProj);
//						backProjection.addAtIndex(idxX, idxY, value);
//					}
//			}
//		}
//		backProjection.show("backprojection");
//		return backProjection;
//	}