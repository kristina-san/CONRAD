package edu.stanford.rsl.tutorial.icecream;

import java.lang.Math;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.*;

public class PhantomK extends Grid2D {

	public PhantomK(int size) {
		super(size, size);
		createPhantom(size);
	}

	private void createPhantom(int size) {

//		//Circle
//		int ii;
//		int jj;
//		for (int i = 0; i < size; i++)
//		for (int j = 0; j < size; j++) {
//			ii = i - size/3;
//			jj = j - size/3;
//			if ((ii * ii + jj * jj) < 12 * size) {
//				this.setAtIndex(i, j, 15);
//			}
//		}
		//	quadrant
		for (int i = 20; i < size; i++)
			for (int j = 20; j < size; j++) {
				if ((i * i + j * j) < 12 * size) {
					this.setAtIndex(i, j, 15);
				}
			}
		//	rectangle
		for (int i = size / 2; i < size; i++)
			for (int j = size / 2; j < size; j++) {
				this.setAtIndex(i, j, 100);
			}
		//	line
		for (int i = 0; i < 7 * size / 8; i++)
			this.setAtIndex(i, size / 3, 50);

	}

	public Grid2D createSinogram(int numberProj, float detectorSpacing,
			int numberDetPixel, float SID) {

		Grid2D sinogram = new Grid2D(numberProj, numberDetPixel);
		for (int theta = 0; theta < numberProj; theta++) {
			double alpha = ((2 * Math.PI / (numberProj)) * theta);

			for (int detPixel = 0; detPixel < numberDetPixel; detPixel++) {

				double xspacePhantom = this.getSpacing()[0] / 2;

				for (double t = -SID; t < SID; t = t + xspacePhantom) {
					// Umrechnung von Detektorkoord. in Weltkoord.
					float r = detPixel * detectorSpacing - (detectorSpacing * (numberDetPixel - 1) / 2);
					// Phantom
					float x = (float) ((float) r * Math.cos(alpha) + t* Math.sin(alpha));
					float y = (float) ((float) r * Math.sin(alpha) - t* Math.cos(alpha));

					double[] id = physicalToIndex(x, y);

					if (id[0] >= 0 && id[0] < this.getWidth()) {
						if (id[1] >= 0 && id[1] < this.getHeight()) {
							sinogram.addAtIndex(
									theta,
									detPixel,
									(float) (xspacePhantom * InterpolationOperators.interpolateLinear(this, id[0], id[1])));
						}
					}
				}
			}
		}
		sinogram.show();
		return sinogram;
	}
	

	public Grid2D createFanogram(boolean mode, int numberProj, float detectorSpacing,		//mode true is short-scan, mode false full-scan
			int numberDetPixel, float dSD, int rotAngle, float dSI, float fanAngle) {

		float delta;
		if(mode)
			//delta = (float) ((Math.PI + fanAngle)/numberProj);
			numberProj = (int) Math.ceil((Math.PI+fanAngle)*180/Math.PI/rotAngle)+1;
		//else
			//delta = (float) ((2*Math.PI)/numberProj);
			//numberProj = (int) ((Math.PI+fanAngle)*180/Math.PI);
		Grid2D fanogram = new Grid2D(numberDetPixel, numberProj);		
		
		double sourcePos[] = { -dSI, 0 };
		for (int beta = 0; beta < numberProj; beta ++) {
			// Umrechnung in RAD
			//double alpha = delta*beta*rotAngle; //((Math.PI * (beta*rotAngle)) / 180);
			double alpha  = ((Math.PI * (beta*rotAngle)) / 180); 
			double[][] rot = { { Math.cos(alpha), Math.sin(alpha) }, { -Math.sin(alpha), Math.cos(alpha) } };
			double[][] rotInv = { { Math.cos(alpha), -Math.sin(alpha) }, { Math.sin(alpha), Math.cos(alpha) } };
			//sourcePos = multiply(rot, sourcePos);
			
			for (int detPixel = 0; detPixel < numberDetPixel; detPixel++) {
				
				float tMax = (numberDetPixel-1) * detectorSpacing;
				double xdet = sourcePos[0] + dSD;
				double ydet = sourcePos[1] - (tMax/2) + detPixel* detectorSpacing;
				//float tAngle = (float) Math.atan(ydet/dSD);
				
				double xspacePhantom = this.getSpacing()[0] / 2;
				//double yspacePhantom = this.getSpacing()[1] / 2;

				// float ray = (float) (dSD/Math.cos(t));
				//double x = sourcePos[0];
				double y = sourcePos[1];
				
				for (double x = sourcePos[0]; x < xdet; x += xspacePhantom) {
					
					y = (x-sourcePos[0])*ydet/dSD;
					//y = (x-sourcePos[0])*Math.tan(tAngle);
					double[] pos = { x, y };
					// Weltkoord!
					double[] id = multiply(rot, pos);
					double[] index = physicalToIndex(id[0], id[1]);

					if (index[0] >= 0 && index[0] < this.getWidth()) {
						if (index[1] >= 0 && index[1] < this.getHeight()) {
							fanogram.addAtIndex(detPixel, beta, InterpolationOperators.interpolateLinear(this, index[0], index[1]));
						}
					}
				}
			}
		}
		return fanogram;
	}

	// matrix-vector multiplication (y = A * x)
	public double[] multiply(double[][] A, double[] x) {
		int m = A.length;
		int n = A[0].length;
		if (x.length != n)
			throw new RuntimeException("Illegal matrix dimensions.");
		double[] y = new double[m];
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				y[i] += (A[i][j] * x[j]);
		return y;
	}

	public Grid2D rebinning(Grid2D fanogram, float dSI, float dSD, int rotAngle, float detectorSpacing) {

		Grid2D sinogram = new Grid2D(fanogram.getHeight(), fanogram.getWidth()); // same sized fano and sinogram
		float numberProj = fanogram.getHeight();
		float numberDetPixel = fanogram.getWidth();
		
		/*for (int t = 0; t < sinogram.getWidth(); t++) {
			for (int beta = 0; beta < sinogram.getHeight(); beta++) {
				double ydet = - (numberDetPixel * detectorSpacing) / 2 + t* detectorSpacing;
				float tvalue = (float) (Math.atan(ydet/dSD)); // wert an der stelle t ausrechnen!
				float betavalue = beta*rotAngle; // wert an der stelle beta ausrechen!
				
				int s= (int) (dSI* Math.sin(tvalue));
				int theta = (int) (tvalue + betavalue);
				sinogram.addAtIndex(theta, s, InterpolationOperators.interpolateLinear(fanogram, t, beta));
			}
		}*/
		for (int s = 0; s < sinogram.getHeight(); s++) {
			for (int iTheta = 0; iTheta < sinogram.getWidth(); iTheta++) {
				double smm=s*detectorSpacing - ((sinogram.getHeight()-1)*detectorSpacing/2);
				double t= smm*dSD/(Math.sqrt(dSD*dSD-smm*smm));
				
				double thetaInterval = (2 * Math.PI / (numberProj));
				double theta = thetaInterval*iTheta;
				double gamma=Math.PI/2-Math.atan(t/dSD);
				double beta=theta-gamma;
				
//				double alpha = (2 * Math.PI / (numberProj)) * theta;
				//double thetaRad = theta/180*Math.PI;
				//double r = s * detectorSpacing;
				//double t = Math.asin(r/dSD);
				//double beta = thetaRad - t;
				
				if (beta<0){
					beta = beta + (numberProj-1)*thetaInterval;
				}
				if (beta > 2*Math.PI)
					beta = beta - (numberProj-1)*thetaInterval;
				
//				double ydet = - (numberDetPixel * detectorSpacing) / 2 + t* detectorSpacing;
//				double tvalue = Math.atan(sinogram.getHeight()/dSD);
//				double betavalue = beta*rotAngle;
//				
				double tIndex = t/detectorSpacing + (fanogram.getWidth()-1)/2;
				double betaIndex = beta/(rotAngle*Math.PI/180);

				sinogram.addAtIndex(iTheta, s, InterpolationOperators.interpolateLinear(fanogram, tIndex, betaIndex));
			}
		}
		
		sinogram.show("rebinning");
		return sinogram;
	}

	public Grid2D rebinningShortScan(Grid2D fanogram, float dSI, float dSD, int rotAngle, float detectorSpacing, float fanAngle) {

		//Identical rays:
		//	t1 = − t2
		//	beta2 = beta1 − 2*t1 + Math.PI
		
		Grid2D sinogram = new Grid2D(360, fanogram.getWidth());
		for (int s = 0; s < sinogram.getHeight(); s++) {
			for (int iTheta = 0; iTheta < sinogram.getWidth(); iTheta++) {
				double smm=s*detectorSpacing - ((sinogram.getHeight()-1)*detectorSpacing/2);
				double t= smm*dSD/(Math.sqrt(dSD*dSD-smm*smm));
				
				//double thetaInterval = (2 * Math.PI / (numberProj));
				double thetaInterval = (2*Math.PI)/sinogram.getWidth();
				double theta = thetaInterval*iTheta;
				double gamma=Math.PI/2-Math.atan(t/dSD);
				double beta=theta-gamma;

				double maxAngle = (Math.PI + fanAngle);

				if (beta<0){
					beta = beta + 2*Math.PI;
				}
				
				if (beta > 2*Math.PI){
					beta = beta - 2*Math.PI;
				}
	
				if(beta > maxAngle)
				{
					beta = beta + Math.PI - 2*Math.atan(t/dSD);
					t = -t;
				}
				
				if (beta<0){
					beta = beta + 2*Math.PI;
				}
				
				if (beta > 2*Math.PI){
					beta = beta - 2*Math.PI;
				}
				
				double tIndex = t/detectorSpacing + (fanogram.getWidth()-1)/2;
				double betaIndex = beta/(rotAngle*Math.PI/180);
				
				//System.out.println("s = " + s + " iTheta = " + iTheta + " beta =" + beta + " t = " + t +  " betaIndex = " + betaIndex + " tIndex = " + tIndex);
				sinogram.addAtIndex(iTheta, s, InterpolationOperators.interpolateLinear(fanogram, tIndex, betaIndex));
			}
		}
		
		sinogram.show("rebin short scan");
		return sinogram;
	}	

	   
	public static void main(String[] args) {
		new ImageJ();
		int size = 256;
		PhantomK p = new PhantomK(size);
		p.setSpacing(0.1, 0.1);
		p.setOrigin(-(size - 1) * p.spacing[0] / 2, -(size - 1) * p.spacing[1]/ 2);
		p.show();
		
		//Ex 1
		float min = NumericPointwiseOperators.min(p);
		float mean = NumericPointwiseOperators.mean(p);

		float d = (float) (Math.sqrt(2) * p.getHeight() * p.getSpacing()[0]);
		float detectorSpacing = (float) 0.2;
		float dSI = 3*d;
		float dSD = 6*d;
		float detNew = 2*(d*dSD)/(2*dSI);
		int numberProj = 360;
		int rotAngleSpacing = 1;
		Grid2D sino = p.createSinogram(180, detectorSpacing, (int)(d/detectorSpacing), d/2 );
		Grid2D fanogram  = p.createFanogram(false, numberProj, detectorSpacing, (int) (detNew / detectorSpacing), dSD, rotAngleSpacing, dSI, 0);
		fanogram.show("fanogram");
		Grid2D sinoFromFano  = p.rebinning(fanogram, dSI, dSD, rotAngleSpacing, detectorSpacing);
		sinoFromFano.setSpacing(360/sinoFromFano.getSize()[0], detectorSpacing);
		sinoFromFano.setOrigin(-(sinoFromFano.getSize()[0]-1)*sinoFromFano.getSpacing()[0]/2, 
				-(sinoFromFano.getSize()[1]-1)*sinoFromFano.getSpacing()[1]/2 );
		
		// Ex 3.3 FBP
		PBP projection = new PBP();
		float [] pixelSpacingRecon = {(float) 0.2, (float) 0.2};
		Grid1DComplex ramp = projection.fftrampFilter(detectorSpacing, sinoFromFano);
		Grid2D filteredRamLak = projection.filtering(sinoFromFano, ramp);
		filteredRamLak.show();
		Grid2D reconstructedFilterRamLak = projection.backProjection(360, detectorSpacing, sinoFromFano.getHeight(), filteredRamLak, size, pixelSpacingRecon);

		// Ex 3.4 Short scan
		int numberDetPixel = (int) (detNew / detectorSpacing);
		float detectorLength = numberDetPixel * detectorSpacing;
		float fanAngle = (float) (2.0 * Math.atan2((detectorLength-1)/2, dSD));
		
		Grid2D shortFanogram  = p.createFanogram(true, numberProj, detectorSpacing, (int) (detNew / detectorSpacing), dSD, rotAngleSpacing, dSI, fanAngle);
		shortFanogram.show("shortfanogram");
		Grid2D sinoFromShortFan  = p.rebinningShortScan(shortFanogram, dSI, dSD, rotAngleSpacing, detectorSpacing, fanAngle);
		
		sinoFromShortFan.setSpacing(360/sinoFromShortFan.getSize()[0], detectorSpacing);
		sinoFromShortFan.setOrigin(-(sinoFromShortFan.getSize()[0]-1)*sinoFromShortFan.getSpacing()[0]/2, 
				-(sinoFromShortFan.getSize()[1]-1)*sinoFromShortFan.getSpacing()[1]/2 );
		
		Grid1DComplex rampShort = projection.fftrampFilter(detectorSpacing, sinoFromShortFan);
		Grid2D filteredRampShort = projection.filtering(sinoFromShortFan, rampShort);
		filteredRampShort.show("filtered short sino");
		Grid2D reconstructedFilterRamLakShort = projection.backProjection(360, detectorSpacing, sinoFromShortFan.getHeight(), filteredRampShort, size, pixelSpacingRecon);
	}

}