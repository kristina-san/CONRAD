package edu.stanford.rsl.tutorial.icecream;

import java.io.IOException;
import java.nio.FloatBuffer;

import ij.ImageJ;

import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLMemory.Mem;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;


public class OpenCL extends OpenCLGrid2D {

	public OpenCL(Grid2D input, CLContext context, CLDevice device) {
		super(input, context, device);
	}

	public static OpenCL createGrid1(int size, CLContext context, CLDevice device )
	{
		Grid2D p = new Grid2D(size, size);
		//Circle
		int ii;
		int jj;
		for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			ii = i - size/4;
			jj = j - size/4;
			if ((ii * ii + jj * jj) < 12 * size) {
				p.setAtIndex(i, j, 15);
			}
		}
		OpenCL pCL = new OpenCL(p, context, device);
		pCL.show("grid1");
		return pCL;
	}
	
	public static OpenCL createGrid2(int size, CLContext context, CLDevice device )
	{
		Grid2D p = new Grid2D(size, size);
		//Circle
		int ii;
		int jj;
		for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			ii = i - size/2;
			jj = j - size/2;
			if ((ii * ii + jj * jj) < 12 * size) {
				p.setAtIndex(i, j, 100);
			}
		}
		OpenCL pCL = new OpenCL(p, context, device);
		pCL.show("grid2");
		return pCL;
	}
	
	public static void AddPhantomToCPUandGPU(PhantomK p, OpenCL phantomCL)
	{
		int number = 1000;
		
		long starttime= System.nanoTime();
	
		//on CPU
		for (int i = 0; i < number; i++){
            NumericPointwiseOperators.addBy(p, p);
        }
		
		long endtime = System.nanoTime();
		long timecost = endtime - starttime;
		
		System.out.println("Time with CPU " + timecost);
		
		starttime= System.nanoTime();
		
		//openCL on GPU
		for (int i = 0; i < number; i++){
            NumericPointwiseOperators.addBy(phantomCL, phantomCL);
        }
		
		endtime= System.nanoTime();
		timecost= endtime - starttime; 
		
		System.out.println("Time on GPU " + timecost);
	}
	
	public void AddTwoOpenCLGrid2Ds(OpenCL grid2, CLContext context, CLDevice device, int size )
	{
		CLProgram program = null;
		try {
			program = context.createProgram(this.getClass().getResourceAsStream("AddOpenCLGrid2D.cl"))
					.build();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		int gridSizeX = size;
		int gridSizeY = size;
		
		int imageSize = this.getSize()[0] * this.getSize()[1];
		
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		CLBuffer<FloatBuffer> imageBuffer = context.createFloatBuffer(imageSize, Mem.READ_ONLY);
		for (int i=0;i<this.getBuffer().length;++i){
			imageBuffer.getBuffer().put(this.getBuffer()[i]);
		}
		imageBuffer.getBuffer().rewind();
		CLImage2d<FloatBuffer> imageGrid1 = context.createImage2d(
				imageBuffer.getBuffer(), this.getSize()[0], this.getSize()[1], format);
		
		CLBuffer<FloatBuffer> imageBuffer2 = context.createFloatBuffer(imageSize, Mem.READ_ONLY);
		for (int i=0;i<this.getBuffer().length;++i){
			imageBuffer2.getBuffer().put(grid2.getBuffer()[i]);
		}
		imageBuffer2.getBuffer().rewind();
		CLImage2d<FloatBuffer> imageGrid2 = context.createImage2d(
				imageBuffer2.getBuffer(), grid2.getSize()[0], grid2.getSize()[1], format);

		// create memory for result grid
		CLBuffer<FloatBuffer> resultGrid = context.createFloatBuffer(imageSize, Mem.WRITE_ONLY);
		
		// copy params
		CLKernel kernel = program.createCLKernel("AddOpenCLGrid2D");
		kernel.putArg(resultGrid).putArg(imageBuffer).putArg(imageBuffer2)
			.putArg(gridSizeX).putArg(gridSizeY);
		
		// createCommandQueue
		CLCommandQueue queue = device.createCommandQueue();
		queue
			//.putWriteImage(imageGrid1, true)
			//.putWriteImage(imageGrid2, true)
			.putWriteBuffer(resultGrid, true)
			.putWriteBuffer(imageBuffer, true)
			.putWriteBuffer(imageBuffer2, true)
			.put2DRangeKernel(kernel, 0, 0,(long)gridSizeX,(long)gridSizeY,1, 1)
			.finish()
			.putReadBuffer(resultGrid, true)
			.finish();
		
		// write resultGrid back to grid2D
		Grid2D result = new Grid2D(size,size);
		result.setSpacing(0.1, 0.1);
		resultGrid.getBuffer().rewind();
		for (int i = 0; i < result.getBuffer().length; ++i) {
			result.getBuffer()[i] = resultGrid.getBuffer().get();
		}
		result.show("addedPicture");
		System.out.println("End");
		
		imageBuffer.release();
		imageBuffer2.release();
	}
	
	public void OpenCLBackProjection(CLContext context, CLDevice device, 
			int numberProj, float detectorSpacing, int numberDetPixel, int sizeRecon, float pixelSpacingRecon[]) 
	{
		//this. is the sinogram in OpenCL that has to be reconstructed
		CLProgram program = null;
		try {
			program = context.createProgram(this.getClass().getResourceAsStream("OpenCL_BP.cl"))
					.build();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		int gridReconSizeX = sizeRecon;
		int gridReconSizeY = sizeRecon;
		int imageSize = gridReconSizeX * gridReconSizeY;
		
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		// create memory for result backprojection grid
		CLBuffer<FloatBuffer> resultBPGrid = context.createFloatBuffer(imageSize, Mem.WRITE_ONLY);
		
		// create buffer for sinogram
		CLBuffer<FloatBuffer> sinoBuffer = context.createFloatBuffer(this.getHeight() * this.getWidth(), Mem.READ_ONLY);
		for (int i=0;i<this.getBuffer().length;++i){
			sinoBuffer.getBuffer().put(this.getBuffer()[i]);
		}
		sinoBuffer.getBuffer().rewind();

		// copy params
		CLKernel kernel = program.createCLKernel("OpenCL_BP");
		kernel.putArg(resultBPGrid).putArg(sinoBuffer)
			.putArg(numberProj).putArg(detectorSpacing).putArg(numberDetPixel).putArg(sizeRecon).putArg(pixelSpacingRecon[0]).putArg(pixelSpacingRecon[1]);
		
		// createCommandQueue
		CLCommandQueue queue = device.createCommandQueue();
		queue
			//.putWriteImage(imageGrid1, true)
			//.putWriteImage(imageGrid2, true)
			.putWriteBuffer(resultBPGrid, true)
			.putWriteBuffer(sinoBuffer, true)
			.put2DRangeKernel(kernel, 0, 0,(long)gridReconSizeX,(long)gridReconSizeY,1, 1)
			//.put2DRangeKernel(kernel, 0, 0, globalWorkSizeBeta, globalWorkSizeT, localWorkSize, localWorkSize)  maybe add this worksize?
			.finish()
			.putReadBuffer(resultBPGrid, true)
			.finish();
		
		// write resultGrid back to grid2D
		Grid2D result = new Grid2D(sizeRecon, sizeRecon);
		result.setSpacing(0.1, 0.1);
		resultBPGrid.getBuffer().rewind();
		for (int i = 0; i < result.getBuffer().length; ++i) {
			result.getBuffer()[i] = resultBPGrid.getBuffer().get();
		}
		result.show("Backprojection");
		
	}
	
	public static void main(String[] args) {
		new ImageJ();
		
		CLContext context = OpenCLUtil.createContext();
		CLDevice[] devices = context.getDevices();
		CLDevice device = context.getMaxFlopsDevice();
		
		// Exercise Sheet 4 - 1.		
		int size = 256;
		PhantomK p = new PhantomK(size);
		p.show("PhantomK");
		OpenCL phantomCL = new OpenCL(p, context, device);
		phantomCL.show("PhantomCL");
		AddPhantomToCPUandGPU(p, phantomCL);
		
		OpenCL grid1 = createGrid1(size, context, device);
		OpenCL grid2 = createGrid2(size, context, device);

		// Exercise Sheet 4 - 2.
		grid1.AddTwoOpenCLGrid2Ds(grid2, context, device, size);
		
		// Exercise Sheet 4 - 3.
		
		// for creating a sinogram from class PhantomK
		p.setSpacing(0.1, 0.1);
		p.setOrigin(-(size - 1) * p.getSpacing()[0] / 2, -(size - 1) * p.getSpacing()[1]/ 2);
		float d = (float) (Math.sqrt(2) * p.getHeight() * p.getSpacing()[0]);
		float detectorSpacing = (float) 0.1;
		Grid2D sinogram = p.createSinogram(180, detectorSpacing, (int)((int)d/detectorSpacing), d/2 );
		sinogram.setSpacing(360/sinogram.getSize()[0], detectorSpacing);
		sinogram.setOrigin(-(sinogram.getSize()[0]-1)*sinogram.getSpacing()[0]/2, -(sinogram.getSize()[1]-1)*sinogram.getSpacing()[1]/2 );
		float [] pixelSpacingRecon = {(float) 0.2, (float) 0.2};
		int numberProj = 360;
		int numberDetPixel = (int) ((int) d/detectorSpacing);
		
		OpenCL sinogramCL = new OpenCL(sinogram, context, device);

		long starttime= System.nanoTime();
		
		sinogramCL.OpenCLBackProjection(context, device, numberProj, detectorSpacing, numberDetPixel, size, pixelSpacingRecon);
	
		long endtime= System.nanoTime();
		
		System.out.println("Time on GPU for PBP " + (endtime - starttime));		
		
	}
}



