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

//import repast.simphony.context.Context;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;

import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.context.space.grid.GridFactoryFinder;
//import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
//import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.essentials.RepastEssentials;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.grid.BouncyBorders;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.util.collections.IndexedIterable;
import repast.simphony.valueLayer.GridValueLayer;
import repast.simphony.valueLayer.ValueLayerDiffuser;

public class ExtraCellEnviro extends DefaultContext{
//public class ExtraCellEnviro extends DefaultContext<BistableSwitch> implements ContextBuilder<BistableSwitch>{
	//petri dish
	
	private ValueLayerDiffuser nutrients;
	private ValueLayerDiffuser peptides;
	private GridValueLayer nutriLayer;
	private GridValueLayer pepLayer;
	private double diffuseTick;
	private double remTick;
	private Grid myGrid;
	private ArrayList<BistableSwitch> remList;
	//private Network network;
	private File treep;
	private File treem;
	
	private static DateFormat format = new SimpleDateFormat("yyyy.MMM.dd.HH_mm_ss_z");

	
	public ExtraCellEnviro () {

		super("ExtraCellEnviro");
		diffuseTick=0;
		remTick=0;
		remList = new ArrayList<BistableSwitch>();
		//create grid
		Parameters parm = RunEnvironment.getInstance().getParameters();
		int cellX = (Integer)parm.getValue("exXSize");
		int cellY = (Integer)parm.getValue("exYSize");
		int cellZ = (Integer)parm.getValue("exZSize");
		//GridBuilderParameters<BistableSwitch> params = GridBuilderParameters.singleOccupancyND(new RandomGridAdder<BistableSwitch>(), new BouncyBorders(), cellX, cellY, cellZ);
		GridBuilderParameters params = GridBuilderParameters.multiOccupancy2D(new RandomGridAdder(), new BouncyBorders(), cellX, cellY/*, cellZ*/);
		myGrid = GridFactoryFinder.createGridFactory(null).createGrid(BsubtilisParameters.EXGrid, this, params);
		
		int foodconc = (Integer)parm.getValue("nutrientConcentration");
		
		nutriLayer = new GridValueLayer("Nutrients",foodconc,true,new BouncyBorders(), cellX,cellY/*,cellZ*/);
		nutrients = new MyDiffuser(nutriLayer,1.0,0.1);
		this.addValueLayer(nutriLayer);
		
		pepLayer = new GridValueLayer("Peptides",0.0,true,new BouncyBorders(), cellX,cellY/*,cellZ*/);
		peptides = new MyDiffuser(pepLayer,1.0,0.1);
		this.addValueLayer(pepLayer);
		//create agents
		
		//NetworkBuilder nbuilder= new NetworkBuilder("FamilyTree",this, true);
		//network = nbuilder.buildNetwork();

		int ct[] = {0,0,0,0,0};
		//ct[0] = (Integer)parm.getValue("numberofComS");
		//ct[1] = (Integer)parm.getValue("numberofComK");
		//ct[2] = (Integer)parm.getValue("numberofComX");
		//ct[3] = (Integer)parm.getValue("numberofmRNA");
		//ct[4] = ct[3];
		
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		ScheduleParameters sparamsodd = ScheduleParameters.createRepeating(1, 2);
		ScheduleParameters sparamseven = ScheduleParameters.createRepeating(2,2);
		int cells = (Integer)parm.getValue("numberofCells");

		String path = (String) parm.getValue("familyTreePath");
		String d = format.format(new Date());
		String fname = "treep."+d+".txt";
		treep = new File(path,fname);
		
		fname = "treem."+d+".txt";
		treem = new File(path,fname);
		//File tree = new File(path,"tree"+id+".txt");
		
		for (int i = 0; i < cells; i++) {
			
		    double energy;
			if (RunEnvironment.getInstance().isBatch()) {
				energy = RandomHelper.nextDoubleFromTo(0.0, 1.0) * 2 * (Float)parm.getValue("initialMass");
			} else {
				energy = RandomHelper.nextDoubleFromTo(0.0, 1.0) * 2 * (Double)parm.getValue("initialMass");
			}
			ct[0] = RandomHelper.nextIntFromTo(0,(Integer)parm.getValue("numberofComS"));
			ct[1] = RandomHelper.nextIntFromTo(0,(Integer)parm.getValue("numberofComK"));
			ct[2] = RandomHelper.nextIntFromTo(0,(Integer)parm.getValue("numberofComX"));
			ct[3] = RandomHelper.nextIntFromTo(0,(Integer)parm.getValue("numberofmRNA"));
			ct[4] = RandomHelper.nextIntFromTo(0,(Integer)parm.getValue("numberofmRNA"));
			BistableSwitch bs = new BistableSwitch(energy,ct);
			this.add(bs);
			bs.setContext(this);
			bs.setGrid(myGrid);
			bs.setNutrients(nutriLayer);
			bs.setPeptides(pepLayer);
			//bs.setNetwork(network);
			bs.setTreeProt(treep);
			bs.setTreemRNA(treem);
		}
		schedule.schedule(sparamsodd,this,"diffuse");
		schedule.schedule(sparamseven,this,"removeCells");
		//if ((Boolean)parm.getValue("visibleFood")) {
			//for (int i = 0; i < cellX; i++) {
				//for (int j = 0; j < cellY; j++) {
					//	Food f = new Food();
						//f.setContext(this);
					//	f.setGrid(myGrid);
						//f.setNutrients(nutriLayer);
					//	f.setPeptides(pepLayer);
						//this.add(f);
						//myGrid.moveTo(f, i, j/*, k*/);
				//}
			//}
		//}
	}
	
	public BistableSwitch makeNewCell(double energy, int[] aCounts) {
		
		BistableSwitch bs = new BistableSwitch(energy,aCounts);
		bs.setContext(this);
		bs.setGrid(myGrid);
		bs.setNutrients(nutriLayer);
		bs.setPeptides(pepLayer);
		//bs.setNetwork(network);
		bs.setTreeProt(treep);
		bs.setTreemRNA(treem);
		this.add(bs);
		return bs;
	}
	
	public void remCellList(BistableSwitch bs) {
		remList.add(bs);
	}
	
	//Scheduled methods
	// odd tick
	//@ScheduledMethod(start = 1, interval = 1, priority=ScheduleParameters.RANDOM_PRIORITY)
	public void diffuse() {
		//diffusing
		double tick = (double) RepastEssentials.GetTickCount();
		if (tick > diffuseTick) {
			nutrients.diffuse();
			peptides.diffuse();
			diffuseTick = tick;
		}
	}
	//maintenance features on even tick
	public void removeCells() {
		double tick = (double) RepastEssentials.GetTickCount();
		if (tick > remTick) {
			Iterator<BistableSwitch> iter = remList.iterator();
			while (iter.hasNext()) {
				BistableSwitch bs = iter.next();
				bs.removeScheduledActions();
				Iterator aiter=bs.iterator();
				while (aiter.hasNext()) {
					Object obj = aiter.next();
					if (obj instanceof ExtendAgent) {
						ExtendAgent ea = (ExtendAgent)obj;
						if (!ea.isDead()) {
							ea.removeScheduledActions();
							ea.setDead(true);
							ea.setTheContext(null);
							ea.setGrid(null);
							aiter.remove();
						}
					}
				}
				GridPoint loc = myGrid.getLocation(bs);
				Parameters parm = RunEnvironment.getInstance().getParameters();
				int nuts = (Integer)parm.getValue("nutrientConcentration");
				int randnut=RandomHelper.nextIntFromTo(0, nuts);
				double curr=nutriLayer.get(loc.getX(),loc.getY());
				nutriLayer.set(curr+randnut,loc.getX(),loc.getY());
				//Iterable<RepastEdge> edges = network.getEdges(bs);
				//Iterator i = edges.iterator();
				//while (i.hasNext()) {
				//	RepastEdge e = (RepastEdge) i.next();
				//	i.remove();
				//	network.removeEdge(e);
				//}
				bs.setContext(null);
				bs.setGrid(null);
				bs.setNutrients(null);
				bs.setPeptides(null);
				//bs.setNetwork(null);
				bs.setTreeProt(null);
				bs.setTreemRNA(null);
				bs.getProjections().clear();
				this.remove(bs);
				iter.remove();
			}
			remTick = tick;
		}
		if (Runtime.getRuntime().freeMemory()<440000) {
	    	Runtime.getRuntime().gc();
	    }
		//IndexedIterable<BistableSwitch> bagents = this.getObjects(BistableSwitch.class);
		if (tick == 50000) {
			RepastEssentials.EndSimulationRun();
		}
		
	}
}
