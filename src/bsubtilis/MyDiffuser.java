/*COPYRIGHT AND PERMISSION NOTICE
UNC Software:  B. Subtilis ABM
Copyright (C) 2009 The University of North Carolina at Chapel Hill
All rights reserved.

The University of North Carolina at Chapel Hill (“UNC”) and the developers (“Developers”) of 
B. Subtilis ABM (“Software”) give recipient (“Recipient”) and Recipient’s Institution (“Institution”) 
permission to use and copy the software in source and binary forms, with or without modification 
for non-commercial purposes only provided that the following conditions are met:

1)	All copies of Software in binary form and/or source code, related documentation and/or 
other materials provided with the Software must reproduce and retain the above copyright notice, 
this list of conditions and the following disclaimer. 

2)	Recipient and Institution shall not distribute Software to any third parties.

3)	The Software is provided “As Is.” The Developers can not guarantee the provision of technical 
support or consultation for the Software. The Developers may provide a location on a UNC Web Site 
for Recipients to post comments, questions, and suggestions at some time in the future. Recipient 
may provide the Developers with feedback on the use of the Software in their research at that time.  
The Developers and UNC are permitted to use any information Recipient provides in making changes to 
the Software. 

4)	Recipient acknowledges that the Developers, UNC and its licensees may develop modifications to 
Software that may be substantially similar to Recipient’s modifications of Software, and that the 
Developers, UNC and its licensees shall not be constrained in any way by Recipient in UNC’s or its 
licensees’ use or management of such modifications. Recipient acknowledges the right of the Developers
and UNC to prepare and publish modifications to Software that may be substantially similar or 
functionally equivalent to your modifications and improvements, and if Recipient or Institution 
obtains patent protection for any modification or improvement to Software, Recipient and Institution 
agree not to allege or enjoin infringement of their patent by the Developers, UNC or any of UNC’s 
licensees obtaining modifications or improvements to Software from the UNC or the Developers.

5)	Recipient and Developer will acknowledge in their respective publications the contributions made 
to each other’s research involving or based on the Software. The current citations for Software are:

Vasa SM and Giddings MC. Agent-Based Model of the Dynamics of Phenotype Switching in Bacillus 
Subtilis. Manuscript submitted.

6)	Any party desiring a license to use the Software for commercial purposes shall contact The 
Office of Technology Development at UNC at 919-966-3929.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS, CONTRIBUTORS, AND THE UNIVERSITY OF NORTH 
CAROLINA AT CHAPEL HILL "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
EVENT SHALL THE COPYRIGHT OWNER, CONTRIBUTORS OR THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
package bsubtilis;

import org.geotools.data.crs.ForceCoordinateSystemFeatureResults;

import repast.simphony.valueLayer.GridValueLayer;
import repast.simphony.valueLayer.ValueLayerDiffuser;

public class MyDiffuser extends ValueLayerDiffuser {

	public MyDiffuser(GridValueLayer valueLayer, double evaporationConst,
            double diffusionConst){
		super(valueLayer, evaporationConst, diffusionConst, true);
	}

	protected void computeVals() {
		// this is being based on http://www.mathcs.sjsu.edu/faculty/rucker/capow/santafe.html
		int size = valueLayer.getDimensions().size();

		if (size == 1) {
			int width = (int) valueLayer.getDimensions().getWidth();

			double sum;
			double[] newVals = new double[width];
			for (int x = 0; x < width; x++) {
				// sum the cell to the left and the right of the given one
				sum = getValue(x - 1);
				sum += getValue(x + 1);

				double weightedAvg = sum / 2.0;

				// apply the diffusion and evaporation constants
				double oldVal = getValue(x);
				double delta = weightedAvg - oldVal;

				double newVal = (oldVal + delta * diffusionConst) * evaporationConst;
				// bring the value into range [min, max]
				newVals[x] = constrainByMinMax(newVal);
			}
			computedVals = newVals;
		} else if (size == 2) {
			int width = (int) valueLayer.getDimensions().getWidth();
			int height = (int) valueLayer.getDimensions().getHeight();
			double[][] newVals = new double[width][height];
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					// these are the neighbors that are directly north/south/east/west to 
					// the given cell 4 times those that are diagonal to the cell
					double uE = getValue(x + 1, y);
					double uN = getValue(x, y + 1);
					double uW = getValue(x - 1, y);
					double uS = getValue(x, y - 1);

					// these are the neighbors that are diagonal to the given cell
					// they are only weighted 1/4 of the ones that are north/south/east/west
					// of the cell
					double uNE = getValue(x + 1, y + 1);
					double uNW = getValue(x - 1, y + 1);
					double uSW = getValue(x - 1, y - 1);
					double uSE = getValue(x + 1, y - 1);

					// compute the weighted avg, those directly north/south/east/west
					// are given 4 times the weight of those on a diagonal
					double weightedAvg = ((uE+uN+uW+uS)*4 + (uNE + uNW + uSW + uSE)) / 20.0;

					// apply the diffusion and evaporation constants
					double oldVal = getValue(x, y);
					double delta = weightedAvg - oldVal;

					double newVal = (oldVal + delta * diffusionConst) * evaporationConst;

					// bring the value into [min, max]
					newVals[x][y] = constrainByMinMax(newVal);

//					System.out.println("x: " + x + " y: " + y + "val: " + oldVal + " delta: "
//							+ delta + " d: " + newVals[x][y]);
				}
			}
			computedVals = newVals;
		} else if (size == 3) {
			int width = (int) valueLayer.getDimensions().getWidth();
			int height = (int) valueLayer.getDimensions().getHeight();
			int depth = (int) valueLayer.getDimensions().getDepth();
			double[][][] newVals = new double[width][height][depth];
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						newVals[x][y][z] = compute3dVal(x, y, z);
					}
				}
			}
			computedVals = newVals;
		}
	}
	
	private double compute3dVal(double xCoord, double yCoord, double zCoord) {
		double weightedSum = 0;
		int numberEqual = 0;
		// this is 4 if we're dealing with a direct (von-neuman) neighbor of the cell
		// otherwise (we're on a diagonal) it's 1
		double multiplier = 1;
		// this keeps the count of the number of neighbors we found that are valid.
		// this also takes into account that some of the cells are weighted more than others
		// and will have them counted accordingly (4 versus 1)
		int count = 0;
		for (double z = zCoord - 1; z <= zCoord + 1; z++) {
			for (double y = yCoord - 1; y <= yCoord + 1; y++) {
				for (double x = xCoord - 1; x <= xCoord + 1; x++) {
					numberEqual = 0;
					// find out if we're directly north/south/east/west/above/below the given
					// cell (a von-neuman neighbor)
					if (xCoord == x) {
						numberEqual++;
					}
					if (yCoord == y) {
						numberEqual++;
					}
					if (zCoord == z) {
						numberEqual++;
					}
					if (numberEqual == 3) {
						// we're at the source
						continue;
					} else if (numberEqual == 2) {
						// we're a direct neighbor so we get weighted more
						multiplier = 4;
					} else {
						multiplier = 1;
					}
//					System.out.println(x + "," + y + "," + z);
					weightedSum += multiplier * getValue(x, y, z);
					count += multiplier;
				}
			}
		}

		double weightedAvg = weightedSum / count;
		double oldVal = getValue(xCoord, yCoord, zCoord);
		double delta = weightedAvg - oldVal;

		return constrainByMinMax((oldVal + delta * diffusionConst) * evaporationConst);
	}
	
	public void diffuse() {
		computeVals();
		int size = valueLayer.getDimensions().size();

		if (size == 1) {
			double[] newVals = (double[]) computedVals;
			for (int x = 0; x < newVals.length; x++) {
				valueLayer.set(newVals[x], x);
			}
		} else if (size == 2) {
			double[][] newVals = (double[][]) computedVals;
			for (int x = 0; x < newVals.length; x++) {
				for (int y = 0; y < newVals[0].length; y++) {
					valueLayer.set(newVals[x][y], x, y);
				}
			}
		} else {
			double[][][] newVals = (double[][][]) computedVals;
			//System.out.println("x="+newVals[0].length+" Y="+newVals[0][0].length+" Z="+newVals.length);
			for (int x = 0; x < /*valueLayer.getDimensions().getWidth();*/newVals.length; x++) {
				for (int y = 0; y < /*valueLayer.getDimensions().getHeight();*/newVals[0].length; y++) {
					for (int z = 0; z < /*valueLayer.getDimensions().getDepth();*/newVals[0][0].length; z++) {
						valueLayer.set(newVals[x][y][z], x, y, z);
					}
				}
			}
		}
	}
}
