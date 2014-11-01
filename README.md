BSubtilis_ABM
=============

Multi-scale agent-based model of competence phenotype in Bacillus subtilis.  It is dependent upon Repast Simphony v 1.2.  This models a B. subtilis colony and the effects of the gene regulatory network responsible for the competence phenotype.  It is written in Java. <br><br>
This is the source code to support the publication:<br>
Stiegelmeyer SM, Giddings M: Agent-based modeling of competence phenotype switching in Bacillus subtilis. Theoretical Biology and Medical Modelling 2013, 10(1):23. 
<br><br>
This is a multi-scale agent-based model built with Repast Simphony v 1.2 modeling the
gene regulatory network that is responsible for the expression of the competence
phenotype in the soil bacteria Bacillus subtilis. The model only displays the agents defined
by ExtraCellEnviro.  It is not possible to view the internal workings of each ABM created by 
the BistableSwitch class without modifying the code and how it is displayed.
<br><br>
The project structure was developed using Eclipse IDE.  The Repast Simphony 1.2 plugin 
is required to successfully build the project.
<br><br>
Context Builder files:
<br><br>
BSubtilisColony.java: creates the ExtraCellEnviro context
<br><br>
ExtraCellEnviro.java: Builds the top layer of the ABM consisting of cell agents with the class 
BistableSwitch.  This models the extracellular environment where cell agents interact with 
eachother--starvation and cell density conditions are simulated at this level.  Two value layers 
are managed here, nutrients and ComX. Agent rules are scheduled on odd system ticks.  
Maintenance rules like adding or removing agents are scheduled for even system ticks so that
rules can be removed from the scheduler--Repast doesn't allow rules to be removed from the 
scheduler if initiated from a rule in the same system tick.
<br><br>
BistableSwitch.java: This class also builds an ABM and represents the Bacillus subtilis cell but 
acts as an agent to ExtraCellEnviro class.  The intracellular model is fully built by this class by 
creating all agents and scheduling the rules.  Agent rules are scheduled on odd system ticks.  
Maintenance rules like adding or removing agents are scheduled for even system ticks so that
rules can be removed from the scheduler--Repast doesn't allow rules to be removed from the 
scheduler if initiated from a rule in the same system tick.
<br><br>
Agent classes:  Agent class files contain the rules that are executed for that particular agent.
<br><br>
ExtendAgent.java:  all agents extend this class, generic methods are provided by this class
ClpCClpP.java
ComK.java
ComS.java
ComX.java
DegU.java
MecA.java
MRNA.java
Promoter.java
Repressor.java
Ribosome.java
<br><br>
Agent support classes:
<br><br>
Available.java: general class for finding neighbors and moving to neighboring positions.
BsubtilisParameters.java: All probabilities are defined in here.
<br><br>

General classes to aid Repast Framework:
<br><br>
AgentStyle2D.java:  For display of 2-d agents, colors are specified here.
AgentStyle3D.java:  For display of 3-d agents, colors are specified here.
MyDiffuser.java: Found a bug in the diffusion implementation in Repast and fixed it here.  It's only 
with 3-D diffusion that there was a problem.
