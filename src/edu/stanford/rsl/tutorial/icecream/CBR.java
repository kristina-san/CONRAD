package edu.stanford.rsl.tutorial.icecream;

import edu.stanford.rsl.apps.gui.ReconstructionPipelineFrame;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.filtering.rampfilters.RamLakRampFilter;
import edu.stanford.rsl.conrad.filtering.rampfilters.RampFilter;
import edu.stanford.rsl.conrad.opencl.OpenCLBackProjector;
import edu.stanford.rsl.conrad.utils.CONRAD;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import edu.stanford.rsl.tutorial.fan.redundancy.ParkerWeights;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

public class CBR {
	
	public CBR() {}

	public static ImagePlus readImageDataFromFile(String filePath, String fileName)
	{
		try {
			new ImageJ();
			ImagePlus imp = IJ.openImage(filePath + "/" + fileName);
			Grid3D impAsGrid = ImageUtil.wrapImagePlus(imp);
			impAsGrid.show("Data from file");
			return imp;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("There occurs an error while reading the image.");
			e.printStackTrace();
			return null;
		}
	}
	
	public static void main(String[] args) throws Exception {
		String filePath ="/proj/ciptmp/FCRSS15/data";
		String fileName = "DensityProjection_No248_Static60_0.8deg_REFERENCE.tif";
		new ImageJ();		
		// Exercise Sheet 5.1
		ImagePlus imp = readImageDataFromFile(filePath, fileName);
		imp.show();
		Grid3D grid = ImageUtil.wrapImagePlus(imp);
		
		// Exercise Sheet 5.2
		CONRAD.setup();
		ReconstructionPipelineFrame pipeLine = new ReconstructionPipelineFrame();
		pipeLine.setVisible(true);
		
		// Exercise Sheet 5.3
		/*double maxT = 248;
		double fan = 10.0*Math.PI/180.0;
		double focalLength = (maxT/2.0-0.5)/Math.tan(fan);

		double deltaT = 1.d;
			
		for (int i = 0; i < 248; ++i)
		{
			double maxBeta =  (double)(i+1) * Math.PI * 2.0 / 360.0;
			double deltaBeta = maxBeta / 180;
					
			ParkerWeights p = new ParkerWeights(focalLength, maxT, deltaT, maxBeta, deltaBeta);
			grid.setSubGrid(i, p);	
		}	
		grid.show();*/
		
		// Exercise Sheet 5.4
		//RampFilter rFil = new RamLakRampFilter();
		//ImagePlus rampFilterImp = rFil.getImagePlusFromRampFilter(1);
		//rampFilterImp.show();
		OpenCLBackProjector p = new OpenCLBackProjector();
		p.reconstructOffline(imp);
		

	}
}
