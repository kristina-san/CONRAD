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

	public Grid2D createFanogram(int numberProj, float detectorSpacing,
			int numberDetPixel, float dSD, int rotAngle, float dSI) {

		Grid2D fanogram = new Grid2D(numberDetPixel, numberProj);
	
		double sourcePos[] = { -dSI, 0 };
		for (int beta = 0; beta < numberProj; beta ++) {
			// Umrechnung in RAD
			double alpha = ((Math.PI * (beta*rotAngle)) / 180);
			double[][] rot = { { Math.cos(alpha), Math.sin(alpha) }, { -Math.sin(alpha), Math.cos(alpha) } };
			double[][] rotInv = { { Math.cos(alpha), -Math.sin(alpha) }, { Math.sin(alpha), Math.cos(alpha) } };
			//sourcePos = multiply(rot, sourcePos);
			
			for (int detPixel = 0; detPixel < numberDetPixel; detPixel++) {
				
				float tMax = numberDetPixel * detectorSpacing;
				double xdet = sourcePos[0] + dSD;
				double ydet = sourcePos[1] - tMax / 2 + detPixel* detectorSpacing;
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
					id = physicalToIndex(id[0], id[1]);

					if (id[0] >= 0 && id[0] < this.getWidth()) {
						if (id[1] >= 0 && id[1] < this.getHeight()) {
							fanogram.addAtIndex(detPixel, beta, InterpolationOperators.interpolateLinear(this, id[0], id[1]));
						}
					}
				}
			}
		}
		fanogram.show("fanogram");
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
		for (int theta = 0; theta < sinogram.getWidth(); theta++) {
			for (int s = 0; s < sinogram.getHeight(); s++) {
				
//				double alpha = (2 * Math.PI / (numberProj)) * theta;
				double thetaRad = theta/180*Math.PI;
				double r = s * detectorSpacing;
				double t = Math.asin(r/dSD);
				double beta = thetaRad - t;
				
				System.out.println(beta/Math.PI*180 + " - " + t);
				
				if (beta<0){
					beta = beta + 2*Math.PI;
				}
				
//				double ydet = - (numberDetPixel * detectorSpacing) / 2 + t* detectorSpacing;
//				double tvalue = Math.atan(sinogram.getHeight()/dSD);
//				double betavalue = beta*rotAngle;
//				
//				t=t/tvalue;
//				beta=beta/rotAngle;
				//System.out.println(beta/Math.PI*180);
				sinogram.addAtIndex(theta, s, InterpolationOperators.interpolateLinear(fanogram, t, beta/Math.PI*180));
			}
		}
		
		sinogram.show("rebinning");
		return sinogram;
	}

	public static void main(String[] args) {
		new ImageJ();
		int size = 256;
		PhantomK p = new PhantomK(size);
		p.setSpacing(0.1, 0.1);
		p.setOrigin(-(size - 1) * p.spacing[0] / 2, -(size - 1) * p.spacing[1]/ 2);
		p.show();
		
		float min = NumericPointwiseOperators.min(p);
		float mean = NumericPointwiseOperators.mean(p);
		//System.out.println(min);
		//System.out.println(mean);

		float d = (float) (Math.sqrt(2) * p.getHeight() * p.getSpacing()[0]);
		float detectorSpacing = (float) 0.1;
//		float dSI = d/5;
//		float dSD = d/2;
		float dSI = 3*d;
		float dSD = 6*d;
		int rotAngleSpacing = 1;
		//p.createSinogram(180, detectorSpacing, (int)((int)d/detectorSpacing), d/2 );
		// p.createSinogram(180, (float)1.2, 500, d);
		Grid2D fanogram  = p.createFanogram(360, detectorSpacing, (int) ((int) d / detectorSpacing), dSD, rotAngleSpacing, dSI);
		
		Grid2D sinogrammFromFanogramm  = p.rebinning(fanogram, dSI, dSD, rotAngleSpacing, detectorSpacing);
		
		// Ex 3.3 FBP
		/*PBP projection = new PBP();
		float [] pixelSpacingRecon = {(float) 0.2, (float) 0.2};
		Grid1DComplex ramLak = projection.RamLakFilter(detectorSpacing, sinogrammFromFanogramm);
		Grid2D filteredRamLak = projection.filtering(sinogrammFromFanogramm, ramLak);
		Grid2D reconstructedFilterRamLak = projection.backProjection(360, detectorSpacing, (int) ((int) d/detectorSpacing), filteredRamLak, size, pixelSpacingRecon);
		*/
	}

}