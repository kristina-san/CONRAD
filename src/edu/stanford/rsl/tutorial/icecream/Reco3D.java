package edu.stanford.rsl.tutorial.icecream;

// Simple structure holding the information of a 3-D reconstruction 
public class Reco3D
{
	char[]	volume;					///< memory for the 3-D reconstruction
	double[]	origin = new double[3];	///< first voxel coordinate in world coordinate system [mm] -> isocenter: (0,0,0) //FLOAT!!!
	int[]	size = new int[3];		///< side lengths of the reconstruction in voxels
	float voxelSize;				///< voxel side length in [mm]
};
