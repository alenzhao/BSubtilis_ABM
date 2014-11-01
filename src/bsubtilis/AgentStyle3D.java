/*COPYRIGHT AND PERMISSION NOTICE
UNC Software:  B. Subtilis ABM
Copyright (C) 2009 The University of North Carolina at Chapel Hill
All rights reserved.

The University of North Carolina at Chapel Hill (ÒUNCÓ) and the developers (ÒDevelopersÓ) of 
B. Subtilis ABM (ÒSoftwareÓ) give recipient (ÒRecipientÓ) and RecipientÕs Institution (ÒInstitutionÓ) 
permission to use and copy the software in source and binary forms, with or without modification 
for non-commercial purposes only provided that the following conditions are met:

1)	All copies of Software in binary form and/or source code, related documentation and/or 
other materials provided with the Software must reproduce and retain the above copyright notice, 
this list of conditions and the following disclaimer. 

2)	Recipient and Institution shall not distribute Software to any third parties.

3)	The Software is provided ÒAs Is.Ó The Developers can not guarantee the provision of technical 
support or consultation for the Software. The Developers may provide a location on a UNC Web Site 
for Recipients to post comments, questions, and suggestions at some time in the future. Recipient 
may provide the Developers with feedback on the use of the Software in their research at that time.  
The Developers and UNC are permitted to use any information Recipient provides in making changes to 
the Software. 

4)	Recipient acknowledges that the Developers, UNC and its licensees may develop modifications to 
Software that may be substantially similar to RecipientÕs modifications of Software, and that the 
Developers, UNC and its licensees shall not be constrained in any way by Recipient in UNCÕs or its 
licenseesÕ use or management of such modifications. Recipient acknowledges the right of the Developers
and UNC to prepare and publish modifications to Software that may be substantially similar or 
functionally equivalent to your modifications and improvements, and if Recipient or Institution 
obtains patent protection for any modification or improvement to Software, Recipient and Institution 
agree not to allege or enjoin infringement of their patent by the Developers, UNC or any of UNCÕs 
licensees obtaining modifications or improvements to Software from the UNC or the Developers.

5)	Recipient and Developer will acknowledge in their respective publications the contributions made 
to each otherÕs research involving or based on the Software. The current citations for Software are:

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

import java.awt.Color;
import java.awt.Font;
import java.awt.Polygon;

import javax.media.j3d.Shape3D;
import javax.media.j3d.TransparencyAttributes;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.visualization.visualization3D.AppearanceFactory;
import repast.simphony.visualization.visualization3D.ShapeFactory;
import repast.simphony.visualization.visualization3D.style.Style3D;
import repast.simphony.visualization.visualization3D.style.TaggedAppearance;
import repast.simphony.visualization.visualization3D.style.TaggedBranchGroup;
import repast.simphony.visualization.visualization3D.style.Style3D.LabelPosition;

public class AgentStyle3D implements Style3D {
	
	public TaggedAppearance getAppearance(Object obj,
			TaggedAppearance appearance, Object shapeID) {
		if (appearance == null) {
			appearance = new TaggedAppearance();
		}
		if (obj instanceof ComK) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.blue);
		} else if (obj instanceof ComS) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.magenta);
		} else if (obj instanceof ComX) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.green);
		} else if (obj instanceof DegU) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.yellow);
		} else if (obj instanceof MecA) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.cyan);
		} else if (obj instanceof MRNA) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.red);
		} else if (obj instanceof Repressor) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.orange);
		} else if (obj instanceof Ribosome) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.darkGray);
		} else if (obj instanceof ClpCClpP) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.gray);
		} else if (obj instanceof Promoter) {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.white);
		} /*else if (obj instanceof Food) {
			Food f = (Food) obj;
			GridPoint pt = f.getGrid().getLocation(f);
			double val = f.getNutrients().get(pt.getX(),pt.getY(),pt.getZ());
			Parameters parm = RunEnvironment.getInstance().getParameters();
			int foodconc = (Integer)parm.getValue("nutrientConcentration");
			int green = 0;
			if (foodconc > 0) {
				green = Math.round(((float)val)*225/foodconc+30);
			}
			//float scale = (float) val/foodconc;
			//int alpha =255- Math.round(255*scale);
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), new Color(0,green,0)); //green
			AppearanceFactory.setTransparentAppearance(appearance.getAppearance(), TransparencyAttributes.NICEST, 0.5f);
	} */else if (obj instanceof BistableSwitch) {
			if (((BistableSwitch)obj).isComp()==1)
				AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.yellow);
			else {
				int a = ((BistableSwitch)obj).getAge();
				Color c = Color.blue;
				switch (a) {
				case 2:
					c = Color.cyan;
					break;

				case 3:
					c = Color.magenta;
					break;
				case 4:
					c = Color.pink;
					break;
				case 5:
					c = Color.red;
					break;
				case 6:
					c = Color.white;
					break;
				case 7: 
					c = Color.yellow;
					break;
				case 8:
					c = Color.green;
					break;
				case 9:
					c = Color.gray;
					break;
				case 10:
					c = Color.lightGray;
					break;
				}
				AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), c);

			}
		} else {
			AppearanceFactory.setMaterialAppearance(appearance.getAppearance(), Color.white);
		}
		
		return appearance;
	}

	public TaggedBranchGroup getBranchGroup(Object obj,
			TaggedBranchGroup group) {
		if (group == null || group.getTag() == null) {
			group = new TaggedBranchGroup("DEFAULT");
			Shape3D sphere = ShapeFactory.createSphere(.03f, "DEFAULT");
			sphere.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
			group.getBranchGroup().addChild(sphere);

			return group;
		}
		return null;
	}

	public String getLabel(Object obj, String currentLabel) {
		return null;
	}

	public Color getLabelColor(Object obj, Color currentColor) {
		return Color.green;
	}

	public Font getLabelFont(Object obj, Font currentFont) {
		return null;
	}

	public float getLabelOffset(Object obj) {
		return 0.035f;
	}

	public repast.simphony.visualization.visualization3D.style.Style3D.LabelPosition getLabelPosition(
			Object obj,
			repast.simphony.visualization.visualization3D.style.Style3D.LabelPosition curentPosition) {
		return LabelPosition.NORTH;
	}

	public float[] getRotation(Object obj) {
		return null;
	}

	public float[] getScale(Object obj) {
		return null;
	}

}
