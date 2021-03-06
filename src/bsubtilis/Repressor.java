/*COPYRIGHT AND PERMISSION NOTICE
UNC Software:  B. Subtilis ABM
Copyright (C) 2009 The University of North Carolina at Chapel Hill
All rights reserved.

The University of North Carolina at Chapel Hill (�UNC�) and the developers (�Developers�) of 
B. Subtilis ABM (�Software�) give recipient (�Recipient�) and Recipient�s Institution (�Institution�) 
permission to use and copy the software in source and binary forms, with or without modification 
for non-commercial purposes only provided that the following conditions are met:

1)	All copies of Software in binary form and/or source code, related documentation and/or 
other materials provided with the Software must reproduce and retain the above copyright notice, 
this list of conditions and the following disclaimer. 

2)	Recipient and Institution shall not distribute Software to any third parties.

3)	The Software is provided �As Is.� The Developers can not guarantee the provision of technical 
support or consultation for the Software. The Developers may provide a location on a UNC Web Site 
for Recipients to post comments, questions, and suggestions at some time in the future. Recipient 
may provide the Developers with feedback on the use of the Software in their research at that time.  
The Developers and UNC are permitted to use any information Recipient provides in making changes to 
the Software. 

4)	Recipient acknowledges that the Developers, UNC and its licensees may develop modifications to 
Software that may be substantially similar to Recipient�s modifications of Software, and that the 
Developers, UNC and its licensees shall not be constrained in any way by Recipient in UNC�s or its 
licensees� use or management of such modifications. Recipient acknowledges the right of the Developers
and UNC to prepare and publish modifications to Software that may be substantially similar or 
functionally equivalent to your modifications and improvements, and if Recipient or Institution 
obtains patent protection for any modification or improvement to Software, Recipient and Institution 
agree not to allege or enjoin infringement of their patent by the Developers, UNC or any of UNC�s 
licensees obtaining modifications or improvements to Software from the UNC or the Developers.

5)	Recipient and Developer will acknowledge in their respective publications the contributions made 
to each other�s research involving or based on the Software. The current citations for Software are:

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

import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;

public class Repressor extends ExtendAgent {
	
	private Promoter boundProm;
	
	public Repressor() {
		super();
		boundProm = null;
		setName("Repressor");
	}

	public Promoter getBoundProm() {
		return boundProm;
	}

	public void setBoundProm(Promoter boundProm) {
		this.boundProm = boundProm;
	}
	
/*	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void move() {
		System.out.println("Agent "+getName()+" is moving2");
		if (!getStop()) { 
			randomWalk();
		}
	}*/
	
//	@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	//@Override
	public void checkNeighbors() {
		if (!isDead()) {
		Grid myGrid = getGrid();//(Grid)getTheContext().getProjection(BsubtilisParameters.ICGrid);
		GridPoint loc = myGrid.getLocation(this);
		if (!getStop()) {
			ArrayList l = Available.getNeighbors(this, loc, myGrid);
			for (int i = 0; i < l.size(); i++) {
				if (l.get(i) instanceof Promoter) {
					Promoter p = (Promoter)l.get(i);
					if (p.getBoundRep() == null  && p.getType() == BsubtilisParameters.PCOMK) {
						double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
						if (rand < BsubtilisParameters.IP[BsubtilisParameters.REPR][BsubtilisParameters.PROM]) {
							setStop(true);
							p.setBoundRep(this);
							boundProm = p;
						}
					}
				}
			}
		}/* else if (boundProm.getComKdimer1() != null && boundProm.getBoundDegU() != null) {
			double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
			if (rand < BsubtilisParameters.IP[BsubtilisParameters.REPR][BsubtilisParameters.DEGU]) {
				setStop(false);
				//this.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
				Available.moveToAvailableNeighborWithin(this, myGrid, loc, BsubtilisParameters.DISTANCE);
				boundProm.setBoundRep(null);
				boundProm = null;
			}
		} else if (boundProm.getComKdimer1() != null) {
			double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
			if (rand < BsubtilisParameters.IP[BsubtilisParameters.REPR][BsubtilisParameters.COMK]) {
				setStop(false);
				//this.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
				Available.moveToAvailableNeighborWithin(this, myGrid, loc, BsubtilisParameters.DISTANCE);
				boundProm.setBoundRep(null);
				boundProm = null;
			}
		}*/ else { //random disassociation
			double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
			if (rand < BsubtilisParameters.MORE_NOISE) {
				setStop(false);
				//this.moveToAvailableNeighborWithin(BsubtilisParameters.DISTANCE);
				Available.moveToAvailableNeighborWithin(this, myGrid, loc, BsubtilisParameters.DISTANCE);
				boundProm.setBoundRep(null);
				boundProm = null;
			}
		}
		}
	}
}
