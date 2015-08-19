#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable
#pragma OPENCL EXTENSION cl_khr_fp64: enable

__constant sampler_t linearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;


/*float interpolate(float* sinogram, float theta, float s, int numberProj)
{
	int sinogramIndex = ((int)(theta*numberProj+s));
	
	float x = theta;
	int x1 = theta; // lower theta
	int x2 = x1 + 1; // higher theta
	
	float y = s;
	int y2 = floor(s); // lower s
	int y1 = y2 + 1; // higher s
		
	float valQ11 = sinogram[x1*numberProj+y1];
	float valQ21 = sinogram[x2*numberProj+y1];
	float valQ12 = sinogram[x1*numberProj+y2];
	float valQ22 = sinogram[x2*numberProj+y2];
	
	float valR1 = ((x2-x)/(x-x1))*valQ11 + ((x-x1)/(x2-x1))*valQ21;
	float valR2 = ((x2-x)/(x-x1))*valQ12 + ((x-x1)/(x2-x1))*valQ22;
		
	float valR1 = sinogram[theta*numberProj+y1];
	float valR2 = sinogram[theta*numberProj+y2];
	float valP = ((y2-y)/(y2-y1))*valR1 + ((y-y1)/(y2-y1))*valR2;
	
	return 0;
}*/

// OpenCL kernel for backprojecting a sinogram
kernel void OpenCL_BP(
	__global float* backProjectionPic,
	__global float* sinogram,
	int numberProj,
	float detectorSpacing,
	int numberDetPixel,
	int sizeReconPic,
	float pixelSpacingReconX,
	float pixelSpacingReconY,
	float originX,
	float originY) 
{
	const unsigned int x = get_global_id(0);// x index
	const unsigned int y = get_global_id(1);// y index
	const unsigned int idx = y*sizeReconPic + x;
	
	int locSizex = get_local_size(0);
	int locSizey = get_local_size(1);
	
	//check if inside image boundaries
	if ( x > sizeReconPic || y > sizeReconPic)
		return;
		
	//Set Origin from backProjection
	float backProjOriginX = -(sizeReconPic-1)*pixelSpacingReconX/2;
	float backProjOriginY = -(sizeReconPic-1)*pixelSpacingReconY/2;

	float pval = 0.f;
	for (int theta = 0; theta < numberProj; theta++)
	{
		float alpha =  ((2.f*M_PI_F/(numberProj))*theta);
		float indexX = x*pixelSpacingReconX + backProjOriginX;
		float indexY = y*pixelSpacingReconY + backProjOriginY;
		float s = indexX*cos(alpha) + indexY*sin(alpha);
		s = (s - originY)/detectorSpacing;
		//int sinogramIndex = ((int)(theta*numberProj+s));
		//int sinogramIndex = interpolate(sinogram, theta, s, numberProj);
		//float value = sinogram[sinogramIndex];
		//float value = interpolate(sinogram, theta, s, numberProj);
		
		//interpolate
		float y = s;
		float y2 = floor(s); // lower s
		float y1 = y2 + 1.f; // higher s			

		float value = 0.f;
		if ((y2 >= 0.f) && (y2 <= 359.f) && (y1 >= 0.f) && (y1 <= 359.f)){ 
			float valR1 = sinogram[theta+((int)y1)*numberProj];
			float valR2 = sinogram[theta+((int)y2)*numberProj];
			value = (y-y2)*valR1 + (y1-y)*valR2;
		}
		value = M_PI_F*value/numberProj;

		pval += value;
	}
	
	backProjectionPic[idx] = pval;
	return;
}



//	public double[] indexToPhysical(double i, double j, float pixelSpacingRecon[], double origin []) {
//		return new double[] {
//				i * pixelSpacingRecon[0]+ origin[0],
//				j * pixelSpacingRecon[1] + origin[1]
//		};
//	}

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