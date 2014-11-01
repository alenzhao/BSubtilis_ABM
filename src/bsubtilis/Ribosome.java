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
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.collections.IndexedIterable;

public class Ribosome extends ExtendAgent {
	
	private MRNA boundmRNA;
	private int type;
	
	public Ribosome() {
		super();
		boundmRNA=null;
		setName("Ribosome");
		type = 0;
	}
	
	public MRNA getBoundmRNA() {
		return boundmRNA;
	}

	public void setBoundmRNA(MRNA boundmRNA) {
		this.boundmRNA = boundmRNA;
	}
	
	
/*	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void move() {
		System.out.println("Agent "+getName()+" is moving2");
		if (!getStop()) { 
			randomWalk();
		}
	}*/
	
	public int getType() {
		return type;
	}

	//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void checkNeighbors() {
		if (!isDead()){
		if (boundmRNA == null) {
			Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
			GridPoint loc = myGrid.getLocation(this);
			ArrayList neighbors = Available.getNeighbors(this, loc, myGrid);
		
			for (int i = 0; i < neighbors.size(); i++) {
				if (neighbors.get(i) instanceof MRNA) {
					if (!((MRNA)neighbors.get(i)).isDead()){
						double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
						if (rand < BsubtilisParameters.IP[BsubtilisParameters.MRNA][BsubtilisParameters.RIBO]) {
							boundmRNA = (MRNA)neighbors.get(i);
							boundmRNA.setStop(true);
							setStop(true);
							boundmRNA.setTranslate(true);
							type = boundmRNA.getType();
						}
					}
				} 
			} 
		}
		}
	}
	
	public void makeBabyProtein(int ptype) {
		
		//Context c = ContextUtils.getContext(this);
		Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc = myGrid.getLocation(this);
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		double start = RepastEssentials.GetTickCount();
		if ((int)start%2 == 0) {
			start = start + 1;
		}
		ScheduleParameters sparams = ScheduleParameters.createRepeating(start, 2);
		if (ptype == BsubtilisParameters.PCOMK) {
			IndexedIterable<ExtendAgent> cs = getTheContext().getObjects(ComK.class);
			//Parameters parm = RunEnvironment.getInstance().getParameters();
			if (cs.size() <= 200) {
			ComK comk = new ComK();
			comk.setTheContext(this.getTheContext());
			comk.setGrid(this.getGrid());
			getTheContext().add(comk);
			//GridPoint gp = this.getAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
			GridPoint gp = Available.getAvailableNeighborWithin(myGrid, loc, BsubtilisParameters.DISTANCE);
			if (gp != null) {
				//Grid myGrid = (Grid)this.getTheContext().getProjection(BsubtilisParameters.ICGrid);
				myGrid.moveTo(comk, gp.getX(),gp.getY(),gp.getZ());
			}
			//schedule.schedule(comk);
			comk.setMove(schedule.schedule(sparams,comk,"move"));
			comk.setCheckNeighbors(schedule.schedule(sparams,comk,"checkNeighbors"));
			}
		} else {
			IndexedIterable<ExtendAgent> cs = getTheContext().getObjects(ComS.class);
			//Parameters parm = RunEnvironment.getInstance().getParameters();
			if (cs.size() <= 200) {
			ComS coms = new ComS();
			coms.setTheContext(this.getTheContext());
			coms.setGrid(this.getGrid());
			getTheContext().add(coms);
			//GridPoint gp = this.getAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
			GridPoint gp = Available.getAvailableNeighborWithin(myGrid, loc, BsubtilisParameters.DISTANCE);
			if (gp != null) {
				//Grid myGrid = (Grid)this.getTheContext().getProjection(BsubtilisParameters.ICGrid);
				myGrid.moveTo(coms, gp.getX(),gp.getY(),gp.getZ());
			}
			schedule.schedule(sparams,coms,"move");
			//schedule.schedule(coms);
			}
		}
		
	}
	
//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void translation() {
		if (!isDead()) {
		if (boundmRNA != null) {
			double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
			if (rand < BsubtilisParameters.TRANSLATE) {
				//makeBabyProtein(boundmRNA.getType());
				((BistableSwitch)getTheContext()).addToAddList(this);
				Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
				GridPoint loc = myGrid.getLocation(this);
				Available.moveToAvailableNeighborWithin(this, myGrid, loc, BsubtilisParameters.DISTANCE);
				//this.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
				boundmRNA = null;
				setStop(false);
			}
		}
		}
	}
}
