package edu.stanford.rsl.tutorial.icecream;

import edu.stanford.rsl.apps.gui.ReconstructionPipelineFrame;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.Grid3D;
import edu.stanford.rsl.conrad.filtering.CosineWeightingTool;
import edu.stanford.rsl.conrad.filtering.ImageFilteringTool;
import edu.stanford.rsl.conrad.filtering.RampFilteringTool;
import edu.stanford.rsl.conrad.filtering.rampfilters.HanningRampFilter;
import edu.stanford.rsl.conrad.filtering.rampfilters.RamLakRampFilter;
import edu.stanford.rsl.conrad.filtering.rampfilters.RampFilter;
import edu.stanford.rsl.conrad.filtering.redundancy.ParkerWeightingTool;
import edu.stanford.rsl.conrad.filtering.redundancy.TrajectoryParkerWeightingTool;
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
		ParkerWeightingTool p = new ParkerWeightingTool();
		p.configure();
		for(int i = 0; i < 248; i++)
		{
			p.setImageIndex(i);
			Grid2D parkerGrid = p.applyToolToImage(grid.getSubGrid(i));
			grid.setSubGrid(i, parkerGrid);
			
		}
		grid.show();
		
		
		// Exercise Sheet 5.4
		CosineWeightingTool c = new CosineWeightingTool();
		c.configure();
		RampFilteringTool r = new RampFilteringTool();
		RampFilter ramp = new HanningRampFilter();
		ramp.configure();
		r.setRamp(ramp);
		OpenCLBackProjector o = new OpenCLBackProjector();
		o.configure();
		
		ImageFilteringTool [] tools = {c, r, o};
		Grid3D g = ImageUtil.applyFiltersInParallel(grid, tools);
		g.show();

		

	}
}
