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
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedulableAction;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.Schedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.essentials.RepastEssentials;
import repast.simphony.parameter.Parameter;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridDimensions;
import repast.simphony.space.grid.GridPoint;

public class ExtendAgent {
	
	//private int vX, vY, vZ;
	private int ID;
	private String name;
	private Context theContext;
	private boolean stop;
	private Grid grid;
	private ISchedulableAction move;
	private ISchedulableAction checkNeighbors;
	private ISchedulableAction death;
	private ISchedulableAction transcription;
	private ISchedulableAction translation;
	private boolean dead;

	public ExtendAgent() {
		//vX=vY=vZ=0;
		ID=-1;
		name=null;
		stop=false;
		theContext=null;
		grid = null;
		move = null;
		checkNeighbors = null;
		death = null;
		transcription = null;
		translation = null;
		dead = false;
	}
	
	public boolean isDead() {
		return dead;
	}

	public void setDead(boolean dead) {
		this.dead = dead;
	}

	public ISchedulableAction getMove() {
		return move;
	}

	public void setMove(ISchedulableAction move) {
		this.move = move;
	}

	public ISchedulableAction getCheckNeighbors() {
		return checkNeighbors;
	}

	public void setCheckNeighbors(ISchedulableAction checkNeighbors) {
		this.checkNeighbors = checkNeighbors;
	}

	public ISchedulableAction getDeath() {
		return death;
	}

	public void setDeath(ISchedulableAction death) {
		this.death = death;
	}

	public ISchedulableAction getTranscription() {
		return transcription;
	}

	public void setTranscription(ISchedulableAction transcription) {
		this.transcription = transcription;
	}

	public ISchedulableAction getTranslation() {
		return translation;
	}

	public void setTranslation(ISchedulableAction translation) {
		this.translation = translation;
	}

	public Grid getGrid() {
		return grid;
	}

	public void setGrid(Grid myGrid) {
		this.grid = myGrid;
	}

	public String getName() {
		return name;
	}

	public void setName(String string) {
		this.name = string;
	}
	
	public int getID() {
		return ID;
	}

	public void setID(int id) {
		ID = id;
	}
	
	public Context getTheContext() {
		return theContext;
	}

	public void setTheContext(Context theContext) {
		this.theContext = theContext;
	}
	
	@Parameter (displayName = "Stop Moving", usageName = "stop")
	public boolean getStop() {
		return stop;
	}

	public void setStop(boolean stop) {
		this.stop = stop;
	}

	public boolean randomWalk() {

		GridPoint loc = getGrid().getLocation(this);
		boolean result = Available.moveToAvailableNeighborWithin(this, getGrid(),loc,1);
		return result;
	}
	
	public void die() {
		if (!dead) {
			((BistableSwitch)this.getTheContext()).addToRemoveList(this);
			this.setDead(true);
			//this.getTheContext().remove(this);
			this.setTheContext(null);
			this.setGrid(null);
		}
	}
	
	public void removeScheduledActions() {
		
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		if (move != null) {
			//RepastEssentials.CancelAction(move);
			schedule.removeAction(move);
			move=null;
		}
		if (checkNeighbors != null) {
			//RepastEssentials.CancelAction(checkNeighbors);
			schedule.removeAction(checkNeighbors);
			checkNeighbors=null;
		}
		if (death != null) {
			//RepastEssentials.CancelAction(death);
			schedule.removeAction(death);
			death=null;
		}
		if (transcription != null) {
			//RepastEssentials.CancelAction(transcription);
			schedule.removeAction(transcription);
			transcription=null;
		}
		if (translation != null) {
			//RepastEssentials.CancelAction(translation);
			schedule.removeAction(translation);
			translation=null;
		}
		
	}

	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void move() {
		//System.out.println("Agent "+getName()+" is moving");
		if (!stop && !isDead()) { 
			randomWalk();
		}
	}
	
	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void checkNeighbors() {
		
	}
	
	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void death() {
		
	}
	
	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void transcription() {
		
	}
	
	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void translation() {
		
	}
	
	public int isComK() {
		return 0;
	}
	
	public int isComS() {
		return 0;
	}
	
}
