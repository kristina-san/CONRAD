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

//		int imageSize = this.getSize()[0] * this.getSize()[1];
//		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
//		CLBuffer<FloatBuffer> imageBuffer = context.createFloatBuffer(imageSize, Mem.READ_ONLY);
//		for (int i=0;i<this.getBuffer().length;++i){
//			imageBuffer.getBuffer().put(this.getBuffer()[i]);
//		}
//		imageBuffer.getBuffer().rewind();
//		CLImage2d<FloatBuffer> imageGrid = context.createImage2d(
//				imageBuffer.getBuffer(), this.getSize()[0], this.getSize()[1], format);
//		imageBuffer.release();
//		
//		// create memory for result grid
//		CLBuffer<FloatBuffer> resultGrid = context.createFloatBuffer(size*size, Mem.WRITE_ONLY);
//		
//		// copy params
//		CLKernel kernel = program.createCLKernel("AddTwoOpenCLGrid2Ds");
//		kernel.putArg(resultGrid).putArg(resultGrid)
//			.putArg((float)maxT).putArg((float)deltaT)
//			.putArg((float)maxBeta).putArg((float)deltaBeta)
//			.putArg((float)focalLength).putArg(maxTIndex).putArg(maxBetaIndex);
		
//		float[] argument = value.getAsFloatArray();
//		CLBuffer<FloatBuffer> argBuffer = device.getContext().createFloatBuffer(argument.length, Mem.READ_ONLY);
//		argBuffer.getBuffer().put(argument);
//		argBuffer.getBuffer().rewind();
//		kernel.putArgs(clmem).putArg(argBuffer).putArg(elementCount);
	}
	
	
	public static void main(String[] args) {
		//new ImageJ();
		
		CLContext context = OpenCLUtil.createContext();
		CLDevice[] devices = context.getDevices();
		CLDevice device = context.getMaxFlopsDevice();
		
		// Exercise Sheet 4 - 1.		
		int size = 256;
		PhantomK p = new PhantomK(size);
		p.show();
		OpenCL phantomCL = new OpenCL(p, context, device);
		phantomCL.show();
		AddPhantomToCPUandGPU(p, phantomCL);
		
		OpenCL grid1 = createGrid1(size, context, device);
		OpenCL grid2 = createGrid2(size, context, device);

		// Exercise Sheet 4 - 2.
		grid1.AddTwoOpenCLGrid2Ds(grid2, context, device, size);
		
	}
}
