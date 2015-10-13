package edu.stanford.rsl.tutorial.icecream;

import java.util.Arrays;

import mpicbg.imglib.type.numeric.integer.UnsignedIntType;

public class cavarev_evaluation {

	public static int round(float num)	//evtl mit generic loesen
	{
		//static_cast<T>((num>0.0f) ? num+0.5f : num-0.5f);
		int iNum = (int)((num > 0.0f) ? num+0.5f : num-0.5f);
		return iNum;
	}
	
	// Simple structure holding the information of a 3-D reconstruction 
	public class Reco3D
	{
		char[]	volume;					///< memory for the 3-D reconstruction
		float[]	origin = new float[3];	///< first voxel coordinate in world coordinate system [mm] -> isocenter: (0,0,0)
		int[]	size = new int[3];		///< side lengths of the reconstruction in voxels
		float voxelSize;				///< voxel side length in [mm]
	};
	
	public static void main(String[] args) throws Exception {
		
		//calling parameters
		// 0: evaluation dataset
		// 1: 3-D reconstruction to be evaluated
				
		if (args.length != 2) 
		{
			System.out.println("Wrong calling syntax.\n\n" 
			 + "   first argument: path to the evaluation database\n" 
			 + "   second argument: path to the 3-D reconstruction to be evaluated\n");
		}
		
		// collect the input data from the command line
		String f_evdb = args[0];
		String f_reco = args[1];		

		////// read the volume in the CAVAREV 1.0 format. All data is expected in little-endian format.
		// 3*float						3*4 bytes				first voxel x=(x0,x1,x2) origin in world coordinates
		// 3*unsigned int				3*4 bytes				reconstruction size (Sx0, Sx1, Sx2) in voxels
		// 1*float						1*4 bytes				voxel size in mm
		// Sx0*Sx1*Sx2*unsigned char	Sx0*Sx1*Sx2*1 bytes		reconstructed volume in row-major format
		
		Reco3D reco = null;

		/*
		 * TODO: Auslesen aus Binary File f_reco
		 */
		
		// fixed data properties for cavarev
		int N = 133; // number of projection images //const?

		// begin the evaluation
		float ox = -53.0f;
		float oy = 45.0f;
		float oz = 1.5f;
		int rszMatrix = 512;
		int rszSlices = 213;
		float rszVoxel = 0.5f;
		
		int rnumel=rszMatrix*rszMatrix*rszSlices;
		char[] refVolume = new char[rnumel];
		int Q=256; //const?
		int[] dsc_same = new int[Q]; //unsigned?
		int[] dsc_sum = new int[Q];	//unsigned?

		float roz = -0.5f*rszVoxel*((float)(rszSlices)-1.0f) + oz;
		float roy = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + oy;
		float rox = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + ox;
		
		/*
		 * TODO: Auslesen aus Binary File f_evdb ifsevdb
		 */

		for (int i=0; i<N; i++)
		{
			//std::memset(dsc_same, 0, Q*sizeof(unsigned int));
			//std::memset(dsc_sum, 0, Q*sizeof(unsigned int));
			Arrays.fill(dsc_same, 0,  Q, 0);
			Arrays.fill(dsc_sum, 0, Q, 0);
			
			int num_border = 0, num_vessel = 0; // unsigned?
			/*
			 * TODO: read from ifsevdb
			 */
			
			int[] vec_border = new int[num_border]; //unsigned?
			int[] vec_vessel = new int[num_vessel]; //unsigned?
			/*
			 * TODO: read from ifsevdb
			 */

			//std::memset(refVolume, 0, rnumel*sizeof(unsigned char));
			Arrays.fill(refVolume, 0, rnumel*sizeof(unsigned char), 0);

			for (int v=0; v<num_border; v++) //unsigned
				refVolume[vec_border[v]] = (char) -1;

			for (int v=0; v<num_vessel; v++)
				refVolume[vec_vessel[v]] = 1;

			/*
			 * TODO: evtl delete [] vec_border and [] vec_vessel
			 */

			for (int iz=0; iz<rszSlices; iz++)
			{
				float rz = (float)iz*rszVoxel+roz;
				int idx_z = round(1.0f/reco.voxelSize*(rz-reco.origin[2]));
				for (int iy=0; iy<rszMatrix; iy++)
				{
					float ry = (float)iy*rszVoxel+roy;
					int idx_y = round(1.0f/reco.voxelSize*(ry-reco.origin[1]));
					for (int ix=0; ix<rszMatrix; ix++)
					{
						float rx = (float)ix*rszVoxel+rox;
						int idx_x = round(1.0f/reco.voxelSize*(rx-reco.origin[0]));

						int vidx = iz*rszMatrix*rszMatrix+iy*rszMatrix+ix; //unsigned?
						
						if (refVolume[vidx]==0)
							continue;
						
						boolean val1 = refVolume[vidx]>0;
						
						char cval2 = 0; //unsigned?

						if (idx_z>=0 && idx_y>=0 && idx_x>=0 && idx_z<reco.size[2] && idx_y<reco.size[1] && idx_x<reco.size[0])
						{
							int vidxev = idx_z*reco.size[0]*reco.size[1] + idx_y*reco.size[0] + idx_x; //unsigned?
							cval2 = reco.volume[vidxev];
						}

						for (int q=0; q<Q; q++)
						{
							boolean val2 = cval2>=q;
							
							if (val1 && val2)
								dsc_same[q]++;
				
							if (val1)
								dsc_sum[q]++;
								
							if (val2)
								dsc_sum[q]++;
						}
					}
				}
			}

			for (int q=0; q<Q; q++)
			{
				System.out.println((2.0f * (float)(dsc_same[q]) / (float)dsc_sum[q]));
			}
			System.out.println();
		}
			
		// TODO: evtl. cleanup
		//delete [] reco.volume;
		//delete [] refVolume;
		//delete [] dsc_same;
		//delete [] dsc_sum;
		
	}

}

/*
template <class T> inline T round(float num)
{
	return static_cast<T>((num>0.0f) ? num+0.5f : num-0.5f);
}


// Simple structure holding the information of a 3-D reconstruction 
struct Reco3D
{
	unsigned char *	volume;		///< memory for the 3-D reconstruction
	float			origin[3];	///< first voxel coordinate in world coordinate system [mm] -> isocenter: (0,0,0)
	unsigned int	size[3];	///< side lengths of the reconstruction in voxels
	float			voxelSize;	///< voxel side length in [mm]
};

using namespace std;

int main(int argc, char **argv)
{
	// calling parameters
	// 0: evaluation dataset
	// 1: 3-D reconstruction to be evaluated

	if (argc != 3)
	{
		cout << "Wrong calling syntax." << endl << endl << argv[0] << " [0] [1] " << endl
			 << "   0: path to the evaluation database" << endl
			 << "   1: path to the 3-D reconstruction to be evaluated" << endl;

		return 1;
	}

	// collect the input data from the command line
	std::string f_evdb = argv[1];
	std::string f_reco = argv[2];


	////// read the volume in the CAVAREV 1.0 format. All data is expected in little-endian format.
	// 3*float						3*4 bytes				first voxel x=(x0,x1,x2) origin in world coordinates
	// 3*unsigned int				3*4 bytes				reconstruction size (Sx0, Sx1, Sx2) in voxels
	// 1*float						1*4 bytes				voxel size in mm
	// Sx0*Sx1*Sx2*unsigned char	Sx0*Sx1*Sx2*1 bytes		reconstructed volume in row-major format

	Reco3D reco = {0};
	std::ifstream ifsvol(f_reco.c_str(), std::ifstream::binary);

	if (!ifsvol.is_open() || ifsvol.fail())
	{
		cout << "-1" << endl;
		return 1;
	}

	ifsvol.read((char*)(&(reco.origin[0])), 3*sizeof(float));
	ifsvol.read((char*)(&(reco.size[0])), 3*sizeof(unsigned int));
	ifsvol.read((char*)(&(reco.voxelSize)), sizeof(float));

	if (ifsvol.fail())
	{
		cout << "-2" << endl;
		return 1;
	}

	unsigned int numelvol = reco.size[0]*reco.size[1]*reco.size[2];
	reco.volume = new unsigned char[numelvol];
	ifsvol.read((char*)(reco.volume), numelvol*sizeof(unsigned char));
	
	if (ifsvol.fail())
	{
		delete [] reco.volume;
		cout << "-3" << endl;
		return 1;
	}
	ifsvol.close();

	// fixed data properties for cavarev
	const unsigned int N = 133; // number of projection images


	// begin the evaluation
	float ox=-53.0f;
	float oy=45.0f;
	float oz=1.5f;
	int rszMatrix=512;
	int rszSlices=213;
	float rszVoxel=0.5f;

	int rnumel=rszMatrix*rszMatrix*rszSlices;
	char * refVolume = new char[rnumel];
	const int Q=256;
	unsigned int * dsc_same = new unsigned int[Q];
	unsigned int * dsc_sum = new unsigned int[Q];

	float roz = -0.5f*rszVoxel*((float)(rszSlices)-1.0f) + oz;
	float roy = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + oy;
	float rox = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + ox;

	std::ifstream ifsevdb(f_evdb.c_str(), std::ifstream::binary);

	for (int i=0; i<N; i++)
	{
		//cout << "Working on frame " << i << " of " << N << endl;

		std::memset(dsc_same, 0, Q*sizeof(unsigned int));
		std::memset(dsc_sum, 0, Q*sizeof(unsigned int));

		unsigned int num_border=0, num_vessel=0;
		ifsevdb.read((char*)(&num_border), sizeof(unsigned int));
		ifsevdb.read((char*)(&num_vessel), sizeof(unsigned int));
		if (ifsevdb.fail())
		{
			cout << "failed to read from evaluation database" << endl;
			return 1;
		}
		unsigned int * vec_border = new unsigned int[num_border];
		unsigned int * vec_vessel = new unsigned int[num_vessel];
		ifsevdb.read((char*)(vec_border), num_border*sizeof(unsigned int));
		ifsevdb.read((char*)(vec_vessel), num_vessel*sizeof(unsigned int));

		std::memset(refVolume, 0, rnumel*sizeof(unsigned char));

		for (unsigned int v=0; v<num_border; v++)
			refVolume[vec_border[v]] = -1;

		for (unsigned int v=0; v<num_vessel; v++)
			refVolume[vec_vessel[v]] = 1;

		delete [] vec_border;
		delete [] vec_vessel;

		for (int iz=0; iz<rszSlices; iz++)
		{
			float rz = (float)iz*rszVoxel+roz;
			int idx_z = round<int>(1.0f/reco.voxelSize*(rz-reco.origin[2]));
			for (int iy=0; iy<rszMatrix; iy++)
			{
				float ry = (float)iy*rszVoxel+roy;
				int idx_y = round<int>(1.0f/reco.voxelSize*(ry-reco.origin[1]));
				for (int ix=0; ix<rszMatrix; ix++)
				{
					float rx = (float)ix*rszVoxel+rox;
					int idx_x = round<int>(1.0f/reco.voxelSize*(rx-reco.origin[0]));

					unsigned int vidx = iz*rszMatrix*rszMatrix+iy*rszMatrix+ix;
					
					if (refVolume[vidx]==0)
						continue;
					
					bool val1 = refVolume[vidx]>0;
					
					unsigned char cval2 = 0;

					if (idx_z>=0 && idx_y>=0 && idx_x>=0 && idx_z<reco.size[2] && idx_y<reco.size[1] && idx_x<reco.size[0])
					{
						unsigned int vidxev = idx_z*reco.size[0]*reco.size[1] + idx_y*reco.size[0] + idx_x;
						cval2 = reco.volume[vidxev];
					}

					for (int q=0; q<Q; q++)
					{
						bool val2 = cval2>=q;
						
						if (val1 && val2)
							dsc_same[q]++;
			
						if (val1)
							dsc_sum[q]++;
							
						if (val2)
							dsc_sum[q]++;
					}
				}
			}
		}

		for (int q=0; q<Q; q++)
		{
			cout << (2.0f * (float)(dsc_same[q]) / (float)dsc_sum[q]) << " ";
		}
		cout << endl;
	}


	// cleanup
	delete [] reco.volume;
	delete [] refVolume;
	delete [] dsc_same;
	delete [] dsc_sum;

	return 0;
}
*/