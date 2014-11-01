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

import java.util.ArrayList;

import repast.simphony.context.Context;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridDimensions;
import repast.simphony.space.grid.GridPoint;

public class Available {

	public static ArrayList getNeighbors(Object obj, GridPoint loc, Grid myGrid) {
		ArrayList l = new ArrayList();
		//Context context = getTheContext();
		//Grid myGrid = (Grid)context.getProjection(BsubtilisParameters.ICGrid);
		//GridPoint loc = myGrid.getLocation(this);
		GridDimensions dim = myGrid.getDimensions();
		//System.out.println("X="+loc.getX()+" Y="+loc.getY()+" Z="+loc.getZ());
		int fromX = (loc.getX()-1)<0 ? loc.getX() : loc.getX()-1;
		int toX = (loc.getX()+1)>=dim.getWidth() ? loc.getX() : (loc.getX()+1) ;
		int fromY = (loc.getY()-1)<0 ? loc.getY() : loc.getY() - 1;
		int toY = (loc.getY()+1)>=dim.getHeight() ? loc.getY() : (loc.getY()+1);
		int fromZ = (loc.getZ()-1)<0 ? loc.getZ() : loc.getZ()-1;
		int toZ = (loc.getZ()+1)>=dim.getDepth() ? loc.getZ() : (loc.getZ()+1);
		for (int i = fromX; i <= toX; i++) {
			for (int j = fromY; j <= toY; j++) {
				for (int k=fromZ; k <= toZ; k++) {
					Object a = myGrid.getObjectAt(i,j,k);
					if ((a != obj) && (a != null))
						l.add(a);
				}
			}
		}

		return l;
	}
	
	public static GridPoint getAvailableNeighbor(Grid myGrid, GridPoint loc) {
		return getAvailableNeighborWithin(myGrid, loc, 1);
	}
	
	public static GridPoint getAvailableNeighborWithin(Grid myGrid, GridPoint loc, int distance) {
		
		ArrayList l = new ArrayList();

		GridDimensions dim = myGrid.getDimensions();
		int fromX = (loc.getX()-distance)<0 ? loc.getX() : loc.getX()-distance;
		int toX = (loc.getX()+distance)>=dim.getWidth() ? loc.getX() : (loc.getX()+distance) ;
		int fromY = (loc.getY()-distance)<0 ? loc.getY() : loc.getY() - distance;
		int toY = (loc.getY()+distance)>=dim.getHeight() ? loc.getY() : (loc.getY()+distance);
		int fromZ = (loc.getZ()-distance)<0 ? loc.getZ() : loc.getZ()-distance;
		int toZ = (loc.getZ()+distance)>=dim.getDepth() ? loc.getZ() : (loc.getZ()+distance);
		for (int i = fromX; i <= toX; i++) {
			for (int j = fromY; j <= toY; j++) {
				for (int k=fromZ; k <= toZ; k++) {
					Object a = myGrid.getObjectAt(i,j,k);
					if (a == null) {
						GridPoint gp = new GridPoint(i,j,k);
						l.add(gp);
					}
				}
			}
		}
		int len = l.size();
		int pt = 0;
		GridPoint returnpt;
		if (len > 0) {
			pt = RandomHelper.nextIntFromTo(0,(len-1));
			returnpt = (GridPoint)l.get(pt);
		} else {
			returnpt = null;
		}
		return returnpt;
	}
	
	public static GridPoint getAvailableNeighborWithin2D(Grid myGrid, GridPoint loc, int distance) {
		//will wrap around borders
		ArrayList l = new ArrayList();

		GridDimensions dim = myGrid.getDimensions();

		int startx = loc.getX()-1 < 0 ? dim.getWidth()-1 : loc.getX()-1;

		for (int i = -1; i <= 1; i++) {
			int starty = loc.getY()-1 < 0 ? dim.getHeight()-1 : loc.getY()-1;
			for (int j = -1; j <= 1; j++) {
					Object a = myGrid.getObjectAt(startx%dim.getWidth(),starty%dim.getHeight());
					if (a == null) {
						GridPoint gp = new GridPoint(startx%dim.getWidth(),starty%dim.getHeight());
						l.add(gp);
					}
					starty++;
			}
			startx++;
		}
		int len = l.size();
		int pt = 0;
		GridPoint returnpt;
		if (len > 0) {
			pt = RandomHelper.nextIntFromTo(0,(len-1));
			returnpt = (GridPoint)l.get(pt);
		} else {
			returnpt = null;
		}
		return returnpt;
	}
	
	public static boolean moveToAvailableNeighborWithin(Object obj, Grid myGrid, GridPoint loc, int distance) {
		GridPoint loca = getAvailableNeighborWithin(myGrid, loc, distance);
		boolean result = false;
		if (loca != null) {
			//Grid myGrid = (Grid)this.getTheContext().getProjection(BsubtilisParameters.ICGrid);
			result = myGrid.moveTo(obj, loca.getX(),loca.getY(),loca.getZ());
		}
		return result;
	}
	
	public static boolean moveToAvailableNeighborWithin2D(Object obj, Grid myGrid, GridPoint loc, int distance) {
		GridPoint loca = getAvailableNeighborWithin2D(myGrid, loc, distance);
		boolean result = false;
		if (loca != null) {
			result = myGrid.moveTo(obj, loca.getX(),loca.getY());
		}
		return result;
	}
}
