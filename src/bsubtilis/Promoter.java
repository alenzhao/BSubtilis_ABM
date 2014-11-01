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

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.essentials.RepastEssentials;
import repast.simphony.parameter.Parameter;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;

public class Promoter extends ExtendAgent{
	
	private int type;
	private ComK comKdimer1;
	private ComK comKdimer2;
	private ComX boundComX;
	private DegU boundDegU;
	private Repressor boundRep;
	
	public Promoter() {
		super();
		type = BsubtilisParameters.PCOMK;
		comKdimer1 = null;
		comKdimer2 = null;
		boundDegU = null;
		boundRep = null;
		boundComX = null;
		setName("Promoter");
	}
	
	@Parameter (displayName = "Type", usageName = "type")
	public int getType() {
		return type;
	}

	public void setType(int type) {
		this.type = type;
	}
	
	public ComK getComKdimer1() {
		return comKdimer1;
	}

	public void setComKdimer1(ComK comKdimer1) {
		if (this.comKdimer1 != null) {
			setComKdimer2(comKdimer1);
		} else {
			this.comKdimer1 = comKdimer1;
		}
	}

	public ComK getComKdimer2() {
		return comKdimer2;
	}

	public void setComKdimer2(ComK comKdimer2) {
		if (this.comKdimer2 == null) {
			this.comKdimer2 = comKdimer2;
		}
	}

	public DegU getBoundDegU() {
		return boundDegU;
	}

	public void setBoundDegU(DegU boundDegU) {
		this.boundDegU = boundDegU;
	}
	
	public Repressor getBoundRep() {
		return boundRep;
	}

	public void setBoundRep(Repressor boundRep) {
		this.boundRep = boundRep;
	}

	public ComX getBoundComX() {
		return boundComX;
	}

	public void setBoundComX(ComX boundComX) {
		this.boundComX = boundComX;
	}

	//@Override
	public void move() {
		
	}

	public void makeBabymRNA(int ptype)	{
	//make a new mRNA add it to the context and move to within DISTANCE on grid
	//from promoter
		MRNA m = new MRNA();
		m.setType(ptype);
		m.setTheContext(this.getTheContext());
		m.setGrid(this.getGrid());
		//Context c = ContextUtils.getContext(this);
		getTheContext().add(m);
		Grid myGrid = getGrid();//(Grid)this.getTheContext().getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc = myGrid.getLocation(this);
		Available.moveToAvailableNeighborWithin(m, myGrid, loc, BsubtilisParameters.DISTANCE);
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		double startodd = RepastEssentials.GetTickCount();
		if ((int)startodd%2==0) {
			startodd = startodd + 1.0f;
		}
		ScheduleParameters sparams = ScheduleParameters.createRepeating(startodd, 2);
		m.setMove(schedule.schedule(sparams,m,"move"));
		m.setCheckNeighbors(schedule.schedule(sparams,m,"checkNeighbors"));
		m.setDeath(schedule.schedule(sparams,m,"death"));
	}
	
	public void freeBound() {
		
		Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc;
		if (comKdimer1 != null) {
			comKdimer1.freeDimer();
			comKdimer1 = null;
		}
		if (comKdimer2 != null) {
			comKdimer2.freeDimer();
			comKdimer2 = null;
		}
		if (boundDegU != null) {
			boundDegU.setStop(false);
			//boundDegU.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
			loc = myGrid.getLocation(boundDegU);
			Available.moveToAvailableNeighborWithin(boundDegU, myGrid, loc, BsubtilisParameters.DISTANCE);
			boundDegU.setBoundPromoter(null);
			boundDegU = null;
		}
		if (boundRep != null) {
			boundRep.setStop(false);
			loc = myGrid.getLocation(boundRep);
			//boundRep.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
			Available.moveToAvailableNeighborWithin(boundRep, myGrid, loc, BsubtilisParameters.DISTANCE);
			boundRep.setBoundProm(null);
			boundRep = null;
		}
		if (boundComX != null) {
			boundComX.setStop(false);
			loc = myGrid.getLocation(boundComX);
			Available.moveToAvailableNeighborWithin(boundComX, myGrid, loc, BsubtilisParameters.DISTANCE);
			//boundComX.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
			boundComX.setBoundProm(null);
			boundComX = null;
		}
	}
	
//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void transcription() {
		//make mRNAs if activators present
		if (!isDead()){
		if (type == BsubtilisParameters.PCOMK) {
			if (comKdimer1 != null && comKdimer2 != null && boundRep== null) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < BsubtilisParameters.IP[BsubtilisParameters.PROM][BsubtilisParameters.COMK]) {
					//makeBabymRNA(BsubtilisParameters.PCOMK);
					((BistableSwitch) this.getTheContext()).addToAddList(this);
					freeBound();
				}
			} else if ((comKdimer1 != null || comKdimer2 != null) && boundRep == null) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < BsubtilisParameters.LESS_NOISE) {
					//makeBabymRNA(BsubtilisParameters.PCOMK);
					((BistableSwitch) this.getTheContext()).addToAddList(this);
					freeBound();
				}
				
			} else if (comKdimer1 == null && comKdimer2 == null && boundRep == null) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < BsubtilisParameters.MORE_NOISE) {
					//makeBabymRNA(BsubtilisParameters.PCOMK);
					((BistableSwitch) this.getTheContext()).addToAddList(this);
					freeBound();
				}				
			}
		} else {
			if (boundComX != null && boundRep==null) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < BsubtilisParameters.IP[BsubtilisParameters.PROM][BsubtilisParameters.COMX]) {
					//makeBabymRNA(BsubtilisParameters.PCOMS);
					((BistableSwitch) this.getTheContext()).addToAddList(this);
					freeBound();
				}
				
			} else if (boundRep == null) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < BsubtilisParameters.MORE_NOISE) {
					//makeBabymRNA(BsubtilisParameters.PCOMS);
					((BistableSwitch) this.getTheContext()).addToAddList(this);
					freeBound();
				}
			}
		}
		}
		
	}

}
