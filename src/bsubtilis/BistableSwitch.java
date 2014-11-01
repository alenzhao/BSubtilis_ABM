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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedulableAction;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
//import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.essentials.RepastEssentials;
import repast.simphony.parameter.Parameter;
import repast.simphony.parameter.Parameters;
//import repast.simphony.query.NotQuery;
//import repast.simphony.query.space.grid.GridWithin;
import repast.simphony.random.RandomHelper;
//import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.graph.Network;
import repast.simphony.space.grid.BouncyBorders;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridDimensions;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.util.collections.IndexedIterable;
import repast.simphony.valueLayer.GridValueLayer;

public class BistableSwitch extends DefaultContext<ExtendAgent> { 
	//this is the cell
	
	private Context context;
	private Grid<BistableSwitch> grid;
	private Grid<ExtendAgent> myGrid;
	private GridValueLayer nutrients;
	private GridValueLayer peptides;
	//private Network network;
	private File treeprot;
	private File treemrna;

	private double energy;
	private double stepTick;
	private double nutTick;
	private double pepgTick;
	private double pepcTick;
	private double lifeTick;
	private double shoveTick;
	private double dieOffTick;
	private int age;
	private int myID;
	private static int ID=1;
	private boolean dead;
	//private int direction;
	private int state = 0;
	private ArrayList<ExtendAgent> removeList;
	private ArrayList<ExtendAgent> addList;
	private ISchedulableAction scheduledStep;
	private ISchedulableAction scheduledConsN;
	private ISchedulableAction scheduledLife;
	private ISchedulableAction scheduledShove;
	private ISchedulableAction scheduledRemAgents;
	private ISchedulableAction scheduledAddAgents;
	private ISchedulableAction scheduledGenPep;
	private ISchedulableAction scheduledConPep;
	private ISchedulableAction scheduledDieOff;
	
	private static final int shoveit3d[][]=
	   {{-1,1,-1},{0,1,-1},{1,1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{-1,-1,-1},{0,-1,-1},{1,-1,-1},
		{-1,1,0},{0,1,0},{1,1,0},{-1,0,0},{0,0,0},{1,0,0},{-1,-1,0},{0,-1,0},{1,-1,0},
		{-1,1,1},{0,1,1},{1,1,1},{-1,0,1},{0,0,1},{1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}};
	
	private static final int shoveit2d[][]=
	   {{-1,1},{0,1},{1,1},{-1,0},{0,0},{1,0},{-1,-1},{0,-1},{1,-1}};
	
	public BistableSwitch(double e,int[] c) {
		
		super("BistableSwitch"+ID);
		//System.out.println("Total Memory"+Runtime.getRuntime().totalMemory());  
		//System.out.println("Free Memory"+Runtime.getRuntime().freeMemory());    
	    if (Runtime.getRuntime().freeMemory()<440000) {
	    	Runtime.getRuntime().gc();
	    }
		myID=ID;
		ID++;
		Parameters parm = RunEnvironment.getInstance().getParameters();
		int cellX = (Integer)parm.getValue("inXSize");
		int cellY = (Integer)parm.getValue("inYSize");
		int cellZ = (Integer)parm.getValue("inZSize");
		boolean disable = (Boolean)parm.getValue("disablecellmodel");
		
		energy = e;
		age = 1;
		dead=false;
		//direction = -1;

		removeList = new ArrayList<ExtendAgent>();
		addList = new ArrayList<ExtendAgent>();
		
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		double startodd = RepastEssentials.GetTickCount() <=0 ? 1 : RepastEssentials.GetTickCount();
		double starteven;
		if ((int)startodd%2 == 0) { //even
			starteven = startodd;
			startodd = startodd+1.0f;
		} else{ //odd
			starteven = startodd+1.0f;
		}
		ScheduleParameters sparamsodd = ScheduleParameters.createRepeating(startodd, 2);
		ScheduleParameters sparamseven = ScheduleParameters.createRepeating(starteven, 2);
		ScheduleParameters sparams100 = ScheduleParameters.createRepeating(startodd, 100);
		ScheduleParameters sparams50 = ScheduleParameters.createRepeating(startodd, 50);
		
		//create grid
		GridBuilderParameters<ExtendAgent> params = GridBuilderParameters.singleOccupancyND(new RandomGridAdder<ExtendAgent>(), new BouncyBorders(), cellX, cellY, cellZ);
		myGrid = GridFactoryFinder.createGridFactory(null).createGrid(BsubtilisParameters.ICGrid, this, params);
		//System.out.println("Bistable switch create grid");
		
		//add agents
		//Promoter to fixed site
		if (!disable) {
			Promoter pcomK = new Promoter();
			pcomK.setTheContext(this);
			pcomK.setType(BsubtilisParameters.PCOMK);
			pcomK.setGrid(myGrid);
			this.add(pcomK);
			myGrid.moveTo(pcomK, cellX/2,cellY/2,cellZ/2);
			schedule.schedule(sparamsodd,pcomK,"move");
			schedule.schedule(sparamsodd,pcomK,"transcription");
		
			Promoter pcomS = new Promoter();
			pcomS.setType(BsubtilisParameters.PCOMS);
			pcomS.setTheContext(this);
			pcomS.setGrid(myGrid);
			this.add(pcomS);
			myGrid.moveTo(pcomS, cellX/3,cellY/3,cellZ/3);
			schedule.schedule(sparamsodd,pcomS,"move");
			schedule.schedule(sparamsodd,pcomS,"transcription");
		
			for (int i = 0; i < c[1]; i++) {
				ComK comk = new ComK();
				comk.setTheContext(this);
				comk.setGrid(myGrid);
				this.add(comk);
				comk.setMove(schedule.schedule(sparamsodd,comk,"move"));
				comk.setCheckNeighbors(schedule.schedule(sparamsodd,comk,"checkNeighbors"));
			}
		
			for (int i = 0; i < c[0]; i++) {
				ComS coms = new ComS();
				coms.setTheContext(this);
				coms.setGrid(myGrid);
				this.add(coms);
				coms.setMove(schedule.schedule(sparamsodd,coms,"move"));
			}
		
			for (int i = 0; i < c[2]; i++) {
				ComX comx = new ComX();
				comx.setTheContext(this);
				comx.setGrid(myGrid);
				this.add(comx);
				comx.setMove(schedule.schedule(sparamsodd,comx,"move"));
				comx.setCheckNeighbors(schedule.schedule(sparamsodd,comx,"checkNeighbors"));
			}
		
			int numdegu = (Integer)parm.getValue("numberofDegU");
			for (int i = 0; i < numdegu; i++) {
				DegU degu = new DegU();
				degu.setTheContext(this);
				degu.setGrid(myGrid);
				this.add(degu);
				schedule.schedule(sparamsodd,degu,"move");
				schedule.schedule(sparamsodd,degu,"checkNeighbors");
			}
		
			int nummeca = (Integer)parm.getValue("numberofMecA");
			for (int i = 0; i < nummeca; i++) {
				MecA meca = new MecA();
				meca.setTheContext(this);
				meca.setGrid(myGrid);
				this.add(meca);
				schedule.schedule(sparamsodd,meca,"move");
				schedule.schedule(sparamsodd,meca,"checkNeighbors");
			}
		
			int numrepr = (Integer)parm.getValue("numberofRepressors");
			for (int i = 0; i < numrepr; i++) {
				Repressor rep = new Repressor();
				rep.setTheContext(this);
				rep.setGrid(myGrid);
				this.add(rep);
				schedule.schedule(sparamsodd,rep,"move");
				schedule.schedule(sparamsodd,rep,"checkNeighbors");
			}
		
			int numribo = (Integer)parm.getValue("numberofRibosomes");
			for (int i = 0; i < numribo; i++) {
				Ribosome ribo = new Ribosome();
				ribo.setTheContext(this);
				ribo.setGrid(myGrid);
				this.add(ribo);
				schedule.schedule(sparamsodd,ribo,"move");
				schedule.schedule(sparamsodd,ribo,"checkNeighbors");
				schedule.schedule(sparamsodd,ribo,"translation");
			}
		
			int numclpc = (Integer)parm.getValue("numberofClpCClpP");
			for (int i = 0; i < numclpc; i++) {
				ClpCClpP clp = new ClpCClpP();
				clp.setTheContext(this);
				clp.setGrid(myGrid);
				this.add(clp);
				schedule.schedule(sparamsodd,clp,"move");
				schedule.schedule(sparamsodd,clp,"checkNeighbors");
				schedule.schedule(sparamsodd,clp,"death");
			}
		
			for (int i = 0; i < c[3]; i++) {
				MRNA m = new MRNA();
				m.setTheContext(this);
				m.setGrid(myGrid);
				this.add(m);
				m.setMove(schedule.schedule(sparamsodd,m,"move"));
				m.setCheckNeighbors(schedule.schedule(sparamsodd,m,"checkNeighbors"));
				m.setDeath(schedule.schedule(sparamsodd,m,"death"));
			}
			for (int i = 0; i < c[4]; i++) {
				MRNA m = new MRNA();
				m.setTheContext(this);
				m.setGrid(myGrid);
				m.setType(BsubtilisParameters.PCOMS);
				this.add(m);
				m.setMove(schedule.schedule(sparamsodd,m,"move"));
				m.setCheckNeighbors(schedule.schedule(sparamsodd,m,"checkNeighbors"));
				m.setDeath(schedule.schedule(sparamsodd,m,"death"));
			}
		}
		//System.out.println("Bistable Switch created agents");
		stepTick = 0;
		nutTick = 0;
		pepgTick = 0;
		pepcTick = 0;
		lifeTick = 0;
		shoveTick = 0;
		dieOffTick = 0;
		scheduledStep = schedule.schedule(sparamsodd,this,"step");
		scheduledShove = schedule.schedule(sparamsodd,this,"shove");
		if (!disable) {
			scheduledGenPep = schedule.schedule(sparams100,this,"generatePeptide");
			scheduledConPep = schedule.schedule(sparams50,this,"consumePeptide");
		}
		scheduledConsN = schedule.schedule(sparamsodd,this,"consumeNutrients");
		scheduledLife = schedule.schedule(sparamsodd,this,"life");
		scheduledRemAgents = schedule.schedule(sparamseven,this,"removeAgents");
		scheduledAddAgents = schedule.schedule(sparamseven,this,"addAgents");
		scheduledDieOff = schedule.schedule(sparams50, this, "randomDieOff");
	}

	public void setTreeProt(File tree) {
		this.treeprot = tree;
	}
	
	public void setTreemRNA(File tree) {
		this.treemrna = tree;
	}

	public int getMyID() {
		return myID;
	}

	public Context getContext() {
		return context;
	}

	public void setContext(Context context) {
		this.context = context;
	}
	
	public Grid<BistableSwitch> getGrid() {
		return grid;
	}

	public void setGrid(Grid<BistableSwitch> grid) {
		this.grid = grid;
	}
	
	public GridValueLayer getNutrients() {
		return nutrients;
	}

	public void setNutrients(GridValueLayer layer) {
		this.nutrients = layer;
	}
	
	public GridValueLayer getPeptides() {
		return peptides;
	}

	public void setPeptides(GridValueLayer peptides) {
		this.peptides = peptides;
	}
	@Parameter (displayName = "Energy", usageName = "energy")
	public double getEnergy() {
		return energy;
	}
	
	public void setEnergy(double e) {
		energy=e;
	}
	
	@Parameter(displayName="Dead",usageName = "dead")
	public boolean isDead() {
		return dead;
	}
	
	/*public Network getNetwork() {
		return network;
	}

	public void setNetwork(Network network) {
		this.network = network;
	}*/

	@Parameter (displayName = "Age", usageName = "age")
	public int getAge() {
		return age;
	}
	
	public void setAge(int a) {
		age = a;
	}

	/*public int getDirection() {
		return direction;
	}

	public void setDirection(int direction) {
		this.direction = direction;
	}*/

	@Parameter (displayName = "No. of ComK", usageName = "comKSize")
	public int getComKSize() {

		IndexedIterable<ExtendAgent> comk = this.getObjects(ComK.class);
		int count = 0;
		for (int i = 0; i < comk.size(); i++) {
			if (!((ComK)comk.get(i)).isDead()) {
				count++;
			}
		}
		
		Parameters parm = RunEnvironment.getInstance().getParameters();
		int max = (Integer)parm.getValue("comKThreshold");
		double tick = (double) RepastEssentials.GetTickCount();

		if (count >= max ) {
			if (state == 0)
				System.out.println(tick+" Cell "+myID+" competent");
			state = 1;
		} else {
			if (state == 1) 
				System.out.println(tick+" Cell "+myID+" not comp");
			state = 0;
		}
		return count;
	}
	
	@Parameter (displayName = "No. of ComS", usageName = "comSSize")
	public int getComSSize() {

		IndexedIterable<ExtendAgent> coms = this.getObjects(ComS.class);
		return coms.size();
	}
	
	@Parameter(displayName="No. of ComX", usageName = "comXSize")
	public int getComXSize() {
		IndexedIterable<ExtendAgent> comx = this.getObjects(ComX.class);
		return comx.size();
	}
	
	public int getComKNoOfmRNA() {
		
		IndexedIterable<ExtendAgent> mrna = this.getObjects(MRNA.class);
		int count = 0;
		for (int i = 0; i < mrna.size(); i++) {
			if (((MRNA) mrna.get(i)).getType() == BsubtilisParameters.PCOMK && 
				!((MRNA) mrna.get(i)).isDead()) {
				count++;
			}
		}
		return count;
	}
	
	public int isBS() {
		int retval=0;
		if (!dead) retval = 1;
		return retval;
	}
	
	public int isComp() {
		int ncomk=getComKSize();
		int retval=0;
		Parameters p = RunEnvironment.getInstance().getParameters();
		int comKThreshold = (Integer)p.getValue("comKThreshold");
		if (ncomk >= comKThreshold) {
			retval=1;
		}
		return retval;
	}
	
	public void addToRemoveList(ExtendAgent agent) {
		removeList.add(agent);
	}
	
	public void addToAddList(ExtendAgent agent) {
			addList.add(agent);
	}
	
	public void removeScheduledActions () {
		if (scheduledStep != null) {
			RepastEssentials.CancelAction(scheduledStep);
			scheduledStep = null;
		}
		if (scheduledLife != null) {
			RepastEssentials.CancelAction(scheduledLife);
			scheduledLife = null;
		}
		if (scheduledConsN != null) {
			RepastEssentials.CancelAction(scheduledConsN);
			scheduledConsN = null;
		}
		if (scheduledShove != null) {
			RepastEssentials.CancelAction(scheduledShove);
			scheduledShove = null;
		}
		boolean disable = (Boolean)RunEnvironment.getInstance().getParameters().getValue("disablecellmodel");
		if (!disable) {
			if (scheduledGenPep != null) {
				RepastEssentials.CancelAction(scheduledGenPep);
				scheduledGenPep = null;
			}
			if (scheduledConPep != null) {
				RepastEssentials.CancelAction(scheduledConPep);
				scheduledConPep = null;
			}
		}
		if (scheduledRemAgents != null) {
			RepastEssentials.CancelAction(scheduledRemAgents);
			removeList.clear();
			scheduledRemAgents = null;
		}
		if (scheduledAddAgents != null) {
			RepastEssentials.CancelAction(scheduledAddAgents);
			addList.clear();
			scheduledAddAgents = null;
		}
		if (scheduledDieOff != null) {
			RepastEssentials.CancelAction(scheduledDieOff);
			scheduledDieOff = null;
		}
	}
	
	public boolean equalpt(GridPoint pt1, int[] pt2) {
		boolean retval=false;
		if (pt1.getX() == pt2[0] && pt1.getY() == pt2[1] /*&& pt1.getZ() == pt2[2]*/) {
			retval = true;
		}
		return retval;
	}
	
	public void chemotaxis() {
		GridPoint loc = grid.getLocation(this);
		int maxx = grid.getDimensions().getWidth();
		int maxy = grid.getDimensions().getHeight();
		int startx = loc.getX()-1 < 0 ? maxx-1 : loc.getX()-1;
		double maxconc = 0.0f;
		int pt[] = {-1,-1};
		ArrayList l = new ArrayList();
		l.add(pt);
		for (int x = -1; x <= 1; x++ ) {
			int starty = loc.getY()-1 < 0 ? maxy-1 : loc.getY()-1;
			for (int y = -1; y <= 1; y++) {
				if (maxconc < nutrients.get(startx%maxx,starty%maxy)) {
					l.clear();
					int pt2[] = {startx%maxx,starty%maxy};
					l.add(pt2);
					maxconc = nutrients.get(startx%maxx,starty%maxy);
				} else if (maxconc == nutrients.get(startx%maxx,starty%maxy)) {
					int pt2[] = {startx%maxx,starty%maxy};
					l.add(pt2);
				}
				starty++;
			}
			startx++;
		}
		int index = 0;
		if (l.size() > 1) {
			index = RandomHelper.nextIntFromTo(0, l.size()-1);
		}
		pt = (int [])l.get(index);
		if (!equalpt(loc,pt) && pt[0] != -1) {
			if (grid.getObjectAt(pt) == null) {
				grid.moveTo(this, pt);
			} else {
				Available.moveToAvailableNeighborWithin2D(this, grid, loc, 1);
			}
		} else {
			Available.moveToAvailableNeighborWithin2D(this, grid, loc, 1);	
		}
	}
	
	public void divide(double energy) {
		//find random plane to bisect cell
		int plane[] = randomPlane();
		// split comk , coms, comx and mRNA
		int c[] = daughterAgents(plane);
	
		BistableSwitch bs = ((ExtraCellEnviro) context).makeNewCell(energy,c);
		bs.setAge((this.getAge()+1));
		//network.addEdge(this,bs);
		double tick = (double) RepastEssentials.GetTickCount();

		try {
			FileWriter fw = new FileWriter(treeprot, true);
			PrintWriter pw = new PrintWriter(fw);
			int c1 = getComKSize();
			int c2 = bs.getComKSize();
			pw.println(tick+" ( "+this.getMyID()+" : "+c1+" , "+bs.getMyID()+" : "+c2+" ) "+this.getMyID()+" : "+(c1+c2));
			fw.close();
		} catch (IOException e) {
			System.out.println("Something wrong with protein file.");
			e.printStackTrace();
		}

		try {
			FileWriter fw = new FileWriter(treemrna, true);
			PrintWriter pw = new PrintWriter(fw);
			int c1 = getComKNoOfmRNA();
			int c2 = bs.getComKNoOfmRNA();
			pw.println(tick+" ( "+this.getMyID()+" : "+c1+" , "+bs.getMyID()+" : "+c2+" ) "+this.getMyID()+" : "+(c1+c2));
			fw.close();
		} catch (IOException e) {
			System.out.println("Something wrong with protein file.");
			e.printStackTrace();
		}
		//System.out.println();
		if (!Available.moveToAvailableNeighborWithin2D(bs, this.getGrid(), this.getGrid().getLocation(this), 1)){
			GridPoint pt = grid.getLocation(this);
			grid.moveTo(bs, pt.getX(), pt.getY()/*, pt.getZ()*/);
		}
	}
	
	public BistableSwitch occupied() {
		BistableSwitch retval = null;
		GridPoint pt = grid.getLocation(this);
		Iterator iter = grid.getObjectsAt(pt.getX(),pt.getY()/*,pt.getZ()*/).iterator();
		while(iter.hasNext()) {
			Object obj = iter.next();
			if (obj instanceof BistableSwitch) {
				BistableSwitch bs = (BistableSwitch)obj;
				if (!bs.equals(this)) {
					retval = bs;
					break;
				}
			}
		}
		return retval;
	}
	
	//Scheduled methods
	public void removeAgents() {
		
		if (!dead) {
			Iterator<ExtendAgent> iter = removeList.iterator();
			while (iter.hasNext()) {
				ExtendAgent ea = iter.next();
				ea.removeScheduledActions();
				this.remove(ea);
				iter.remove();
			}
		}
	}
	
	public void addAgents() {
		
		if (!dead) {
			Iterator<ExtendAgent> iter = addList.iterator();
			while (iter.hasNext()) {
				ExtendAgent ea = iter.next();
				if (ea instanceof Promoter) {
					if (((Promoter)ea).getType() == BsubtilisParameters.PCOMK) {
						((Promoter)ea).makeBabymRNA(BsubtilisParameters.PCOMK);
					} else {
						((Promoter)ea).makeBabymRNA(BsubtilisParameters.PCOMS);
					}
					iter.remove();  //because of concurrent execution error
				} else if (ea instanceof Ribosome){
					((Ribosome)ea).makeBabyProtein(((Ribosome)ea).getType());
					iter.remove();
				}
			}
		}
	}
	
	public void disableRepressor() {
		IndexedIterable<ExtendAgent> rep = this.getObjects(Repressor.class);
		int s = rep.size();
		if (s > 0) {
			int i = RandomHelper.nextIntFromTo(0, s-1);
			Repressor r = (Repressor)rep.get(i);
			//r.setDead(true);
			if (r.getBoundProm() != null) {
				r.getBoundProm().setBoundRep(null);
				r.setBoundProm(null);
			}
			r.setStop(false);
			r.die();
		}
	}
	
	public void randomDieOff() {
		if (!dead) {
			double tick= (double) RepastEssentials.GetTickCount();
			if (tick > dieOffTick) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				double threshold;
				if (RunEnvironment.getInstance().isBatch()) {
		       		threshold = (Float)RunEnvironment.getInstance().getParameters().getValue("threshold");

		       	} else {
		       		threshold = (Double)RunEnvironment.getInstance().getParameters().getValue("threshold");
		       	}
				if (energy < (threshold-1)/2  && energy > 0.5) {
					if (rand < BsubtilisParameters.MORE_NOISE){
						dead=true;
						((ExtraCellEnviro)context).remCellList(this);
					}
				}
				dieOffTick = tick;
			}
		}
	}
	
	public void step() {
		if (!dead) {
			int tmp = getComKSize();
			//System.out.println("Step: Bistable"+myID);
			double tick = (double) RepastEssentials.GetTickCount();
			if (tick > stepTick) {
				
				Parameters p = RunEnvironment.getInstance().getParameters();
				double carry;
				double rate;
				double deathrate;
				double threshold;
		       	if (RunEnvironment.getInstance().isBatch()) {
		       		carry = (Float)p.getValue("carryingCapacity");
		       		rate = (Float)p.getValue("growthRate");
		       		deathrate=(Float)p.getValue("deathRate");
		       		threshold = (Float)p.getValue("threshold");

		       	} else {
		       		carry = (Double)p.getValue("carryingCapacity");
		       		rate = (Double)p.getValue("growthRate");
		       		deathrate=(Double)p.getValue("deathRate");
		       		threshold = (Double)p.getValue("threshold");
		       	}
		       	
				GridPoint pt = getGrid().getLocation(this);
				double food = nutrients.get(pt.getX(),pt.getY()/*,pt.getZ()*/);
				if (food < 1) {
					double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < 0.95) {
						if (energy <= 1.5) {
							energy = energy - deathrate/(carry/2.0f);
						} else {
							energy = energy - rate*energy*energy/carry;
						}
						//energy = energy - deathrate*energy;
					}
					rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < 0.5) {
						chemotaxis();
					}
					
					rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < BsubtilisParameters.LESS_NOISE) {
						disableRepressor();
					}
					//}
				} else {
					double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					//if (rand < 0.99) {
						energy = energy - rate*energy*energy/carry; //metabolic energy reduction
					//}
				}
				stepTick = tick;
			}
		} else {
			System.out.println("not dead");
		}
	}
	
	public void shove() {
		if (!dead) {
			double tick = (double) RepastEssentials.GetTickCount();
			if (tick > shoveTick) {
				BistableSwitch occ = occupied();
				if (occ != null) {
					int i = RandomHelper.nextIntFromTo(0, 8/*26*/);
					//System.out.println("shoving to i="+i);
					if (!Available.moveToAvailableNeighborWithin2D(this, grid, grid.getLocation(this), 1)) {
						GridPoint pt = grid.getLocation(this);
						int x = (pt.getX()+shoveit2d[i][0]);
						if (x < 0) {
							x = grid.getDimensions().getWidth()-1;
						} else if (x >= grid.getDimensions().getWidth()) {
							x=x%grid.getDimensions().getWidth();
						}
						int y = pt.getY()+shoveit2d[i][1];
						if (y < 0) {
							y = grid.getDimensions().getHeight()-1;
						} else if (y >= grid.getDimensions().getHeight()) {
							y = y%grid.getDimensions().getHeight();
						}
						grid.moveTo(this,x,y/*,z*/);
					}
				}
			}
		}
	}
	
	public void generatePeptide() {
		if (!dead) {
		double tick = (double) RepastEssentials.GetTickCount();
		if (tick > pepgTick) {
			double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
			if (rand < 0.8) {
				GridPoint pt = getGrid().getLocation(this);
				double ncomx = peptides.get(pt.getX(),pt.getY()/*,pt.getZ()*/);
				peptides.set(ncomx+1, pt.getX(),pt.getY()/*,pt.getZ()*/);
			}
			pepgTick = tick;
		}
		}
	}
	
	public void consumePeptide() {
		if (!dead) {
		double tick = (double) RepastEssentials.GetTickCount();
		if (tick > pepcTick) {
			GridPoint pt = getGrid().getLocation(this);
			double comx = peptides.get(pt.getX(),pt.getY()/*,pt.getZ()*/);
			//System.out.println("Peptides="+comx);
			if (comx >= 1) {
				double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				if (rand < 0.8) {
					ComX c = new ComX();
					c.setTheContext(this);
					c.setGrid(myGrid);
					this.add(c);
					ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
					double startodd = RepastEssentials.GetTickCount();
					if ((int)startodd%2 == 0) { //even
						startodd = startodd+1.0f;
					}
					ScheduleParameters sparamsodd = ScheduleParameters.createRepeating(startodd, 2);
					c.setMove(schedule.schedule(sparamsodd,c,"move"));
					c.setCheckNeighbors(schedule.schedule(sparamsodd,c,"checkNeighbors"));
					peptides.set(comx-1, pt.getX(),pt.getY()/*,pt.getZ()*/);
				}
			}
			pepcTick = tick;
		}
		}
	}
	
	public void consumeNutrients() {
		if (!dead) {
		double tick = (double) RepastEssentials.GetTickCount();
		if (tick > nutTick) {
			GridPoint pt = getGrid().getLocation(this);
			double food = nutrients.get(pt.getX(),pt.getY()/*,pt.getZ()*/);
			//System.out.println("Food = "+food);
			if (food >= 1) {
				//double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
				//if (rand < 0.9) {
				Parameters p = RunEnvironment.getInstance().getParameters();
				//double mass = (Double)p.getValue("initialMass");
				double rate;
		       	if (RunEnvironment.getInstance().isBatch()) {
		       		rate = (Float)p.getValue("growthRate");
		       	} else {
		       		rate = (Double)p.getValue("growthRate");
		       	}
				
				energy = energy + rate*energy;
				nutrients.set(food-1, pt.getX(),pt.getY()/*,pt.getZ()*/);
				//}
			} else {
				//System.out.println("no food");
			}
			nutTick = tick;
		}
		}
	}
	// not scheduled
	public int[] randomPlane() {
		GridDimensions dim = myGrid.getDimensions();
		int x1 = dim.getWidth()/2;
		//int x1 = RandomHelper.nextIntFromTo(0, dim.getWidth());
		int x2 = RandomHelper.nextIntFromTo(0, dim.getWidth());
		int x3 = RandomHelper.nextIntFromTo(0, dim.getWidth());
		int y1 = dim.getHeight()/2;
		//int y1 = RandomHelper.nextIntFromTo(0, dim.getHeight());
		int y2 = RandomHelper.nextIntFromTo(0, dim.getHeight());
		int y3 = RandomHelper.nextIntFromTo(0, dim.getHeight());
		int z1 = dim.getDepth()/2;
		//int z1 = RandomHelper.nextIntFromTo(0, dim.getDepth());
		int z2 = RandomHelper.nextIntFromTo(0, dim.getDepth());
		int z3 = RandomHelper.nextIntFromTo(0, dim.getDepth());
		int coef[]= {0,0,0,0};
		coef[0] = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2);
		coef[1] = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2);
		coef[2] = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
		coef[3] = -1*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1));
		return coef;
	}
	//not scheduled
	public int[] daughterAgents(int[] p) {
		
		//Disrupt promoter agents in parent
		IndexedIterable<ExtendAgent> prom = this.getObjects(Promoter.class);
		for (int i = 0; i < prom.size(); i++) {
			((Promoter)prom.get(i)).freeBound();
		}
		//count 
		IndexedIterable<ExtendAgent> coms = this.getObjects(ComS.class);
		IndexedIterable<ExtendAgent> comk = this.getObjects(ComK.class);
		IndexedIterable<ExtendAgent> comx = this.getObjects(ComX.class);
		IndexedIterable<ExtendAgent> mrna = this.getObjects(MRNA.class);
		List<ExtendAgent> rlist = new ArrayList<ExtendAgent>();
		int count[] = {0,0,0,0,0};
		//System.out.println("size coms="+coms.size());
		for (int i=0; i < coms.size(); i++) {
			if (((ComS)coms.get(i)).getBoundAdapter()==null && !((ComS)coms.get(i)).isDead()) {
				GridPoint pt = myGrid.getLocation(coms.get(i));
				int sign = pt.getX()*p[0]+pt.getY()*p[1]+pt.getZ()*p[2] + p[3];
				if (sign >=0 ) {
					//this.remove((ComS)coms.get(i));
					//((ComS)coms.get(i)).die();
					rlist.add(coms.get(i));
					//System.out.println(((ComS)coms.get(i)).getName()+((ComS)coms.get(i)).getAgentID()+"has been removed");
					count[0] = count[0]+1;
				}
			}
		}
		for (int i =0; i < comk.size(); i++) {
			if (((ComK)comk.get(i)).getBoundComK()==null && ((ComK)comk.get(i)).getBoundAdapter()==null && !((ComK)comk.get(i)).isDead()) {
				GridPoint pt = myGrid.getLocation(comk.get(i));
				int sign = pt.getX()*p[0]+pt.getY()*p[1]+pt.getZ()*p[2] + p[3];
				if (sign >=0 ) {
					//this.remove((ComK)comk.get(i));
					//((ComK)comk.get(i)).die();
					rlist.add(comk.get(i));
					//System.out.println(((ComK)comk.get(i)).getName()+((ComK)comk.get(i)).getAgentID()+"has been removed");
					count[1] = count[1] + 1;
				}
			}
		}
		for (int i = 0; i < comx.size(); i++) {
			GridPoint pt = myGrid.getLocation(comx.get(i));
			int sign = pt.getX()*p[0]+pt.getY()*p[1]+pt.getZ()*p[2] + p[3];
			if (sign >=0 ) {
				//this.remove(comx.get(i));
				//((ComX)comx.get(i)).die();
				rlist.add(comx.get(i));
				//System.out.println(((ComX)comx.get(i)).getName()+((ComX)comx.get(i)).getAgentID()+"has been removed");
				count[2] = count[2] + 1;
			}
		}
		for (int i = 0; i < mrna.size(); i++) {
			if (!((MRNA)mrna.get(i)).getTranslate() && !((MRNA)mrna.get(i)).isDead()) {
				GridPoint pt = myGrid.getLocation(mrna.get(i));
				int sign = pt.getX()*p[0]+pt.getY()*p[1]+pt.getZ()*p[2] + p[3];
				if (sign >=0 ) {
					//this.remove(mrna.get(i));
					//((MRNA)mrna.get(i)).die();
					rlist.add(mrna.get(i));
					//System.out.println(((MRNA)mrna.get(i)).getName()+((MRNA)mrna.get(i)).getAgentID()+"has been removed");
					if (((MRNA)mrna.get(i)).getType() == BsubtilisParameters.PCOMK) {
						count[3] = count[3] + 1;
					} else {
						count[4] = count[4] + 1;
					}
				}
			}
		}
		//stupid java
		for (Object o : rlist) {
			((ExtendAgent)o).die();
		}
		return count;
	}

	
	public void life() {
		if (!dead) {
			double tick = (double) RepastEssentials.GetTickCount();
			if (tick > lifeTick) {
				Parameters p = RunEnvironment.getInstance().getParameters();
				int cKThreshold = (Integer)p.getValue("comKThreshold");
				double threshold;
		       	if (RunEnvironment.getInstance().isBatch()) {
		       		threshold = (Float)p.getValue("threshold");
		       	} else {
		       		threshold = (Double)p.getValue("threshold");
		       	}
				if (energy > threshold) {
					double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < 0.8) {
						if (isComp()==0) {
							energy = energy / 2;
							//energy = RandomHelper.nextDoubleFromTo((energy/2.0f), threshold);
							divide(threshold-energy);
						}
					}
				}
				
				/*if (energy <= 1.5) {
					double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < (BsubtilisParameters.MORE_NOISE/2)) {
						dead=true;
						((ExtraCellEnviro)context).remCellList(this);
					}
				} else */if (energy < 0.5) {
					double rand = RandomHelper.nextDoubleFromTo(0.0, 1.0);
					if (rand < BsubtilisParameters.MORE_NOISE) {
						dead=true;
						((ExtraCellEnviro)context).remCellList(this);
					}
				}
				lifeTick = tick;
			}
		} else {
			System.out.println("not dead");
		}
	}
}
