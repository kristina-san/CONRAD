package edu.stanford.rsl.tutorial.icecream;

import java.util.Arrays;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import edu.stanford.rsl.conrad.utils.Configuration;

public class CavarevEvaluation {

	private String filename = "";
	private String evaluationFile = "";
	
	public static int round(double d) {
		int iNum = (int)((d > 0.0f) ? d+0.5f : d-0.5f);
		return iNum;
	}
	
	public CavarevEvaluation(String evaluationFile, boolean hasRespiratoryMotion)
	{
		if(Configuration.getGlobalConfiguration() == null){
			Configuration.loadConfiguration();
		}
		if(hasRespiratoryMotion){
			this.filename = Configuration.getGlobalConfiguration().getRegistry().get("CAVAREV_DIRECTORY") + "/cavarev_cardbreath.bin";
		}else{
			this.filename = Configuration.getGlobalConfiguration().getRegistry().get("CAVAREV_DIRECTORY") + "/cavarev_card.bin";
		}
		this.evaluationFile = evaluationFile;
	}
	
	public String run() throws Exception {
		
		// collect the input data
		String f_evdb = filename;	// evaluation dataset
		String f_reco = evaluationFile;		// 3-D reconstruction to be evaluated

		////// read the volume in the CAVAREV 1.0 format. All data is expected in little-endian format.
		// 3*float						3*4 bytes				first voxel x=(x0,x1,x2) origin in world coordinates
		// 3*unsigned int				3*4 bytes				reconstruction size (Sx0, Sx1, Sx2) in voxels
		// 1*float						1*4 bytes				voxel size in mm
		// Sx0*Sx1*Sx2*unsigned char	Sx0*Sx1*Sx2*1 bytes		reconstructed volume in row-major format
		
		Reco3D reco = new Reco3D();

		try {
			// Read from the 3-D reconstruction that has to be evaluated
			FileInputStream fStream = new FileInputStream(f_reco);
			// Number of matrices is given as the total size of the file
			// divided by 4 bytes per float, divided by 12 floats per projection matrix
			DataInputStream in = new DataInputStream(fStream);
		
			byte[] buffer = new byte[4];
			in.read(buffer);
			reco.origin[0] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getFloat();
			in.read(buffer);
			reco.origin[1] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getFloat();
			in.read(buffer);
			reco.origin[2] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getFloat();
			in.read(buffer);
			reco.size[0] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			in.read(buffer);
			reco.size[1] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			in.read(buffer);
			reco.size[2] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			in.read(buffer);
			reco.voxelSize = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getFloat();
			
			buffer = new byte[1];
			int numelvol = reco.size[0]*reco.size[1]*reco.size[2];
			reco.volume = new char[numelvol];
			for(int i = 0; i < numelvol; i++) {
				in.read(buffer);
				reco.volume[i] = (char) ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).get();
			}
			in.close();
			fStream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return "An error with your file occured.";
		} catch (IOException e) {
			e.printStackTrace();
			return "An error with your file occured.";
		}
		
		// fixed data properties for cavarev
		int N = 133; // number of projection images (phases)

		// begin the evaluation
		float ox = -53.0f;
		float oy = 45.0f;
		float oz = 1.5f;
		int rszMatrix = 512;
		int rszSlices = 213;
		float rszVoxel = 0.5f;
		
		int rnumel=rszMatrix*rszMatrix*rszSlices;
		int[] refVolume = new int[rnumel];
		int Q=256;
		int[] dsc_same = new int[Q];
		int[] dsc_sum = new int[Q];

		float roz = -0.5f*rszVoxel*((float)(rszSlices)-1.0f) + oz;
		float roy = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + oy;
		float rox = -0.5f*rszVoxel*((float)(rszMatrix)-1.0f) + ox;
		
		// Read from the evaluation dataset
		FileInputStream fStream = new FileInputStream(f_evdb);
		DataInputStream in = new DataInputStream(fStream);
		byte[] buffer = new byte[4];
		
		float max = -1.0f;
		int frame = -1;
		
		// Get over 133 phases
		for (int i=0; i<N; i++)
		{
			Arrays.fill(dsc_same, 0,  Q, 0);
			Arrays.fill(dsc_sum, 0, Q, 0);
			
			int num_border = 0, num_vessel = 0;

			in.read(buffer);
			num_border = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			in.read(buffer);
			num_vessel = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			
			int[] vec_border = new int[num_border];
			int[] vec_vessel = new int[num_vessel];

			for(int n = 0; n < num_border; n++){
				in.read(buffer);
				vec_border[n] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();
			}
			for(int n = 0; n < num_vessel; n++) {
				in.read(buffer);
				vec_vessel[n] = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN).getInt();				
			}

			Arrays.fill(refVolume, 0, rnumel, 0);

			for (int v=0; v<num_border; v++)
				refVolume[vec_border[v]] = -1;

			for (int v=0; v<num_vessel; v++)
				refVolume[vec_vessel[v]] = 1;

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

						int vidx = iz*rszMatrix*rszMatrix+iy*rszMatrix+ix;
						
						if (refVolume[vidx]==0)
							continue;
						
						boolean val1 = refVolume[vidx]>0;
						
						char cval2 = 0;

						if (idx_z>=0 && idx_y>=0 && idx_x>=0 && idx_z<reco.size[2] && idx_y<reco.size[1] && idx_x<reco.size[0])
						{
							int vidxev = idx_z*reco.size[0]*reco.size[1] + idx_y*reco.size[0] + idx_x;
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

			System.out.println("Phase " + i + ": ");
			// Evaluate over 256 thresholds
			for (int q=0; q<Q; q++)
			{
				float dice = (2.0f * (float)(dsc_same[q]) / (float)dsc_sum[q]);
				System.out.println(dice + " ");
				if(dice > max) {
					max = dice;
					frame = i;
				}
			}
			System.out.println();

		}
		in.close();
		fStream.close();
		
		// Final output of the maximal dice coefficient at the frame number
		String finalString = "Max. dice coefficient: " + max + " at frame " + frame + ".";
		System.out.println(finalString);
		return finalString;
	}
}