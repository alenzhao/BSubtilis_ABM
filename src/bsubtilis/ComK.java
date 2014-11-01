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
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameter;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;

public class ComK extends ExtendAgent {
	
	private ComK boundComK;
	private Promoter boundPromoter;
	private MecA boundAdapter;
	private int agentID;
	private static int ID=0;

	public ComK () {
		super();
		boundComK=null;
		boundPromoter=null;
		boundAdapter=null;
		setName("ComK");
		agentID = ++ID;
		//System.out.println("ComK"+getAgentID());
		//size++;
	}
	
	//@Override
	public int isComK() {
		return 1;
	}
	
/*	public void deathAction () {
		size--;
	}
	
	public static int getSize() {
		return size;
	}*/

	public boolean isClear() {
	
		boolean val = false;
		if (boundComK == null && boundPromoter == null && boundAdapter == null) 
			val = true;
		return val;
	}
	
	public ComK getBoundComK() {
		return boundComK;
	}

	public void setBoundComK(ComK boundComK) {
		this.boundComK = boundComK;
	}

	public Promoter getBoundPromoter() {
		return boundPromoter;
	}

	public void setBoundPromoter(Promoter boundPromoter) {
		this.boundPromoter = boundPromoter;
	}

	public MecA getBoundAdapter() {
		return boundAdapter;
	}

	public void setBoundAdapter(MecA boundAdapter) {
		this.boundAdapter = boundAdapter;
	}
	
	@Parameter (displayName = "Agent ID", usageName = "agentID")
	public int getAgentID() {
		return agentID;
	}
	
	public void setAgentID(int id) {
		agentID=id;
	}
	
	public void freeDimer() {
		boundComK.setStop(false);
		setStop(false);
		boundComK.setBoundPromoter(null);
		boundPromoter=null;
		//Grid myGrid = (Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc = getGrid().getLocation(boundComK);
		//boundComK.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
		Available.moveToAvailableNeighborWithin(boundComK, getGrid(), loc, BsubtilisParameters.DISTANCE);
		loc = getGrid().getLocation(this);
		//this.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
		Available.moveToAvailableNeighborWithin(this, getGrid(), loc, BsubtilisParameters.DISTANCE);
		boundComK.setBoundComK(null);
		boundComK=null;
	}
	
	public void randomDimerWalk() {
		//Context context = getTheContext();
		Grid projection = getGrid();//(Grid)context.getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc = projection.getLocation(this);
		if (loc != null) {
			boolean result = randomWalk();
			if (result) {
				projection.moveTo(boundComK,loc.getX(),loc.getY(),loc.getZ());
			}
		}
	}
	
//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void move() {
		//System.out.println("Agent "+getName()+getAgentID()+" is moving");
		if (!getStop()  && !isDead()) {
			if (boundComK==null)
				randomWalk();
			else
				randomDimerWalk();
		}
	}
	
//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void checkNeighbors() {
		//System.out.println("ComK"+getAgentID()+" check neighbors");
		if (!isDead()) {
			Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
			GridPoint loc = myGrid.getLocation(this);
			if (loc == null) {
				System.out.println("loc is null");
			}
			//System.out.println("X="+loc.getX()+" Y="+loc.getY()+" Z="+loc.getZ());

			ArrayList neighbors = Available.getNeighbors(this, loc, myGrid);
		
			for (int i = 0; i < neighbors.size(); i++) {
				if (neighbors.get(i) instanceof ComK) {
					if (boundComK == null && boundAdapter==null) {
						if (!((ComK) neighbors.get(i)).getStop() && ((ComK)neighbors.get(i)).getBoundComK()==null) {
							double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
							if (rand < BsubtilisParameters.IP[BsubtilisParameters.COMK][BsubtilisParameters.COMK]) {
								boundComK = (ComK)neighbors.get(i);
								boundComK.setBoundComK(this);
								boundComK.setStop(true);
							}
						}
					}
				} else if (neighbors.get(i) instanceof Promoter) {
					Promoter promo = (Promoter)neighbors.get(i);
					if (promo.getType() == 0) {
						if (canBindToPromoter(promo)) {
							double prob = 0.0;
							if (promo.getBoundDegU()!= null) 
								prob = BsubtilisParameters.IP[BsubtilisParameters.COMK][BsubtilisParameters.DEGU];
							else
								prob = BsubtilisParameters.IP[BsubtilisParameters.COMK][BsubtilisParameters.PROM];
							double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
							if (rand < prob) {
								boundComK.setStop(true);
								setStop(true);
								promo.setComKdimer1(this);
								boundPromoter=promo;
								boundComK.setBoundPromoter(promo);
							}
						}
					}
				} 
			}
		}
	}
	
	public boolean canBindToPromoter(Promoter p) {
		
		boolean result = false;
		if (boundComK != null && boundAdapter==null && boundPromoter==null) {  //if dimer
			if (p.getComKdimer1()==null || p.getComKdimer2()==null) {
				if (p.getComKdimer1()!=this && p.getComKdimer2()!=this) 
					result=true;
			}
		}
		
		return result;
	}

}
