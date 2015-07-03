package edu.stanford.rsl.tutorial.icecream;


import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.*;
import edu.stanford.rsl.conrad.utils.FFTUtil;

public class PBP {
	
	public Grid1DComplex RamLakFilter(float detectorSpacing, Grid2D sinogram){
		int N = sinogram.getWidth();
		int freqSpacing = (int) ((int) 1/(detectorSpacing * (int)sinogram.getWidth()*2));
		Grid1DComplex rampFilter = new Grid1DComplex(N);
		int idx = 0;
		for (int n = -N/2; n <=N/2; n++){ //freqSpacing?
			if(n==0){
				rampFilter.setAtIndex(idx, (float)(1/4));
			} else if (n%2 == 0){
				rampFilter.setAtIndex(idx, (float) 0);
			} else {
				rampFilter.setAtIndex(idx, (float) (-1/(n*n*Math.PI*Math.PI)));
			} 
			idx++;
		} 
		rampFilter.transformForward();
		return rampFilter;
	}
	
	public Grid1DComplex fftrampFilter (float detectorSpacing, Grid2D sinogram)
	{
		int N = sinogram.getWidth();
		
		Grid1DComplex rampFilter = new Grid1DComplex(N);
		int idx = 0;

		int freqSpacing = (int) ((int) 1/(detectorSpacing * (int)sinogram.getWidth()*2)); //Defintion des Filters? Soll wie eine Pyramide aussehn? Welche Grundfunktion
		//for (int n = -N/2; n <=N/2; ){
		for (int n = 0; n < rampFilter.getSize()[0];n++ ){  //freqSpacing?
			if (n <= rampFilter.getSize()[0]/2){
				rampFilter.setAtIndex(idx, Math.abs(n));
			} else
			rampFilter.setAtIndex(idx,-Math.abs(n)+rampFilter.getSize()[0]);
			idx++;
		}
		return rampFilter;
	}
	
	public Grid2D filtering (Grid2D sinogram, Grid1DComplex filter){
		Grid2D filtered = new Grid2D(sinogram.getWidth(), sinogram.getHeight());
		filtered.setSpacing(sinogram.getSpacing());
		filtered.setOrigin(sinogram.getOrigin());
		for (int row =0; row < sinogram.getWidth(); row++){
			Grid1DComplex rowS = new Grid1DComplex(sinogram.getWidth());
			for(int col = 0; col < sinogram.getHeight(); col++){
				rowS.setAtIndex(col, sinogram.getPixelValue(row, col));
			} 
			rowS.transformForward();
			for(int i = 0; i < rowS.getSize()[0]; i++){
				rowS.multiplyAtIndex(i, filter.getAtIndex(i));			
			}
			rowS.transformInverse();
			for(int col = 0; col < sinogram.getHeight(); col ++){
				filtered.setAtIndex(row, col, rowS.getRealAtIndex(col));
			}
			
		} 
		return filtered;
	}
	
	public Grid2D backProjection(int numberProj, float detectorSpacing, int numberDetPixel, Grid2D sinogram, int sizeRecon, float pixelSpacingRecon[]){
		
		Grid2D backProjection = new Grid2D(sizeRecon,  sizeRecon);
		backProjection.setOrigin(-(sizeRecon-1)*pixelSpacingRecon[0]/2, -(sizeRecon-1)*pixelSpacingRecon[1]/2 );

			for (int idxX =0; idxX< backProjection.getWidth(); idxX ++){
				for (int idxY =0; idxY< backProjection.getHeight(); idxY ++){
					for (int theta = 0; theta < numberProj; theta++)
					{
						double alpha =  ((2*Math.PI/(numberProj))*theta);
						double[] id = indexToPhysical(idxX, idxY, pixelSpacingRecon, backProjection.getOrigin());
						double s=id[0]*Math.cos(alpha)+ id[1]*Math.sin(alpha);
						s = (s- sinogram.getOrigin()[1])/detectorSpacing ;
						float value = InterpolationOperators.interpolateLinear(sinogram, theta, s);
						value = (float) (Math.PI*value/numberProj);
						backProjection.addAtIndex(idxX, idxY, value);
					}
			}
		}
		backProjection.show("backprojection");
		return backProjection;
	}

	public double[] indexToPhysical(double i, double j, float pixelSpacingRecon[], double origin []) {
		return new double[] {
				i * pixelSpacingRecon[0]+ origin[0],
				j * pixelSpacingRecon[1] + origin[1]
		};
	}
	
	public double[] physicalToIndex(double x, double y, float detectorSpacing, double origin []) {
		return new double[] {
				(x ) / detectorSpacing- origin[0]
				//(y ) / detectorSpacing  //- this.origin[1]
		};
	}
	
	
	public static void main(String[] args) {
		new ImageJ();
		int size = 256;
		PhantomK p = new PhantomK(size);
		p.setSpacing(0.1, 0.1);
		p.setOrigin(-(size-1)*p.getSpacing()[0]/2, -(size-1)*p.getSpacing()[1]/2 );
		p.show("Phantom");
		float d = (float) (Math.sqrt(2)*p.getHeight()*p.getSpacing()[0]);
		float detectorSpacing = (float) 0.1;
		Grid2D sinogram = p.createSinogram(360, detectorSpacing, (int) ((int) d/detectorSpacing), d/2 );
		sinogram.setSpacing(360/sinogram.getSize()[0], detectorSpacing);
		sinogram.setOrigin(-(sinogram.getSize()[0]-1)*sinogram.getSpacing()[0]/2, -(sinogram.getSize()[1]-1)*sinogram.getSpacing()[1]/2 );
		float [] pixelSpacingRecon = {(float) 0.2, (float) 0.2};
		
		PBP object = new PBP();
		Grid1DComplex rampFilter = object.fftrampFilter(detectorSpacing, sinogram);
		rampFilter.show("rampfiler");
		Grid1DComplex ramLak = object.RamLakFilter(detectorSpacing, sinogram);
		ramLak.show("ramlak");
		Grid2D filtered = object.filtering(sinogram, rampFilter);
		filtered.show("filtered");
		Grid2D filteredRamLak = object.filtering(sinogram, ramLak);
		filtered.show("filteredRamLak");
		
		Grid2D reconstructed = object.backProjection(360, detectorSpacing, (int) ((int) d/detectorSpacing), sinogram, size, pixelSpacingRecon);		
		Grid2D reconstructedFilter = object.backProjection(360, detectorSpacing, (int) ((int) d/detectorSpacing), filtered, size, pixelSpacingRecon);
		Grid2D reconstructedFilterRamLak = object.backProjection(360, detectorSpacing, (int) ((int) d/detectorSpacing), filteredRamLak, size, pixelSpacingRecon);
		
	}
}