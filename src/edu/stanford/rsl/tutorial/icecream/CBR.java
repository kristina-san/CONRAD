package edu.stanford.rsl.tutorial.icecream;

import edu.stanford.rsl.apps.gui.ReconstructionPipelineFrame;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.utils.CONRAD;
import edu.stanford.rsl.conrad.utils.ImageUtil;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;

public class CBR {
	
	public CBR() {}

	public static Grid3D readImageDataFromFile(String filePath, String fileName)
	{
		try {
			new ImageJ();
			ImagePlus imp = IJ.openImage(filePath + "/" + fileName);
			Grid3D impAsGrid = ImageUtil.wrapImagePlus(imp);
			impAsGrid.show("Data from file");
			return impAsGrid;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("There occurs an error while reading the image.");
			e.printStackTrace();
			return null;
		}
	}
	
	public static void main(String[] args) {
		String filePath ="/proj/ciptmp/FCRSS15/data";
		String fileName = "DensityProjection_No248_Static60_0.8deg_REFERENCE.tif";
		
		// Exercise Sheet 5.1
		Grid3D grid = readImageDataFromFile(filePath, fileName);

		// Exercise Sheet 5.2
		CONRAD.setup();
		ReconstructionPipelineFrame pipeLine = new ReconstructionPipelineFrame();
		pipeLine.setVisible(true);
		
		
		

	}
}
