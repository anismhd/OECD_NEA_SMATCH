DOF_MAP = {'DX':1, 'DY':2, 'DZ':3}

import regionToolset
import connectorBehavior
from abaqusConstants import *
import numpy as np

import pickle
from material_prop_reader import material_section_read
from aerb_unv_reader import UNV

############################################################################
# Rigid Link Sspring Characteristics
RIGIDK = 1E12

############################################################################
# Base Isolator Spring Stiffness 
"""
STIFF_X = 9.5E6
STIFF_Y = 9.5E6
STIFF_Z = 1483.0E06"""
STIFF_X = 4.45E6
STIFF_Y = 4.45E6
STIFF_Z = 135.0E06

BaseIsolator = {}
BaseIsolator['R_C_RBB'] = (4*STIFF_X,4*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD6'] = (4*STIFF_X,4*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD5'] = (4*STIFF_X,4*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD2'] = (4*STIFF_X,4*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD1'] = (4*STIFF_X,4*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_RBB'] = (2*STIFF_X,2*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD6'] = (2*STIFF_X,2*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD5'] = (2*STIFF_X,2*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD2'] = (2*STIFF_X,2*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD1'] = (2*STIFF_X,2*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD4'] = (8*STIFF_X,8*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BRI'] = (8*STIFF_X,8*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_E_BD4'] = (8*STIFF_X,8*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_E_BRI'] = (8*STIFF_X,8*STIFF_Y,0.1E+06,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_RBB_Z'] = (0.1E+06,0.1E+06,4*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD6_Z'] = (0.1E+06,0.1E+06,4*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD5_Z'] = (0.1E+06,0.1E+06,4*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD2_Z'] = (0.1E+06,0.1E+06,4*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD1_Z'] = (0.1E+06,0.1E+06,4*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_RBB_Z'] = (0.1E+06,0.1E+06,2*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD6_Z'] = (0.1E+06,0.1E+06,2*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD5_Z'] = (0.1E+06,0.1E+06,2*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD2_Z'] = (0.1E+06,0.1E+06,2*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_R_BD1_Z'] = (0.1E+06,0.1E+06,2*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BD4_Z'] = (0.1E+06,0.1E+06,8*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_C_BRI_Z'] = (0.1E+06,0.1E+06,8*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_E_BD4_Z'] = (0.1E+06,0.1E+06,8*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
BaseIsolator['R_E_BRI_Z'] = (0.1E+06,0.1E+06,8*STIFF_Z,0.1E+06,0.1E+06,0.1E+06,)
####################################### DO NOT CHANGE ANYTHING ###############
# Create Materials
def defineMaterial(matName, MODEL, density, E, mu):
	MODEL.Material(name=matName)
	if density <= 0.001:
		density = 0.001
	MODEL.materials[matName].Density(table=((density, ), ))
	MODEL.materials[matName].Elastic(table=((E,mu), ))

# Create Sets
## Create Nodel Set
def defineNodeSet(part,setName,nodeLabels):
	part.SetFromNodeLabels(name=str(setName),\
		nodeLabels=nodeLabels, unsorted=True)

## Create Element Sets
def defineElementSet(part,setName,elemLabels):
	part.SetFromElementLabels(name=str(setName),\
		elementLabels=elemLabels)

# Create Sections	
## Create Truss Sections
def defineTrussElement(MODEL, name, materialName, area):
	MODEL.TrussSection(name=str(name), material=str(materialName),\
		area=area)

## Create Rectangular Profiles
def defineRectangularProfiles(MODEL,profName,width,depth):
	MODEL.RectangularProfile(name=str(profName), a=width, b=depth)

## Create General Profiles
def defineGeneralProfile(MODEL,name,A,I11,I12,I22,JJ):
	MODEL.GeneralizedProfile(name=str(name), area=A, 
        i11=I11, i12=I12, i22=I22, j=JJ, gammaO=0.0, gammaW=0.0)

## Create Rectangular Sections
def defineRectangularSection(MODEL,profName,secName,matName,width,height):
	defineRectangularProfiles(MODEL,profName,width,height)
	MODEL.BeamSection(name=str(secName), 
        integration=DURING_ANALYSIS, poissonRatio=0.0, profile=str(profName), 
        material=str(matName), temperatureVar=LINEAR, 
        consistentMassMatrix=False)

## Create General Sections
def defineGeneralSection(MODEL,profName,secName, matName, para, E,GG):
	defineGeneralProfile(MODEL,profName,para[0],para[4],0,para[5],para[3])
	MODEL.BeamSection(name=str(secName), 
        integration=BEFORE_ANALYSIS, poissonRatio=0.0, beamShape=CONSTANT, 
        profile=str(profName), table=((E, GG), ), alphaDamping=0.0, 
        betaDamping=0.0, compositeDamping=0.0, centroid=(0.0, 0.0), 
        shearCenter=(0.0, 0.0), consistentMassMatrix=False)
	MODEL.sections[secName].TransverseShearBeam(
        scfDefinition=ANALYSIS_DEFAULT, 
        slendernessCompensation=1.79769313486232e+308,\
		k23=GG*para[1], k13=GG*para[2])

# Section Assignment
## Section Assignment
def defineSectionAssignment(PART, setName, sectionName):
    region = PART.sets[setName]
    PART.SectionAssignment(region=region, sectionName=sectionName, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

def defineBeamOrientation(PART, setName, cosines):
	region = PART.sets[setName]
	PART.assignBeamSectionOrientation(region=region, method=N1_COSINES,\
		n1=cosines)

def defineInertialMass(PART, setName, intertiaName, values):
	region=PART.sets[setName]
	PART.engineeringFeatures.PointMassInertia(\
        name=intertiaName, region=region, \
		mass1=abs(values[0]), mass2=abs(values[1]), mass3=abs(values[2]),\
        i11=abs(values[3]), i22=abs(values[4]), i33=abs(values[5]),\
		i12=0.0, i13=0.0, i23=0.0,alpha=0.0, composite=0.0)
## Create Spring Elements
def findNodeIndx(PART, nodeLabel):
	indx = nodeLabel-1
	if PART.nodes[indx].label == nodeLabel:
		return indx
	if PART.nodes[indx].label > nodeLabel:
		side = True
	else:
		side = False
	while True:
		if PART.nodes[indx].label == nodeLabel:
			return indx
		elif not(side == (PART.nodes[indx].label > nodeLabel)):
			return None
		elif side:
			indx = indx - 1
		else:
			indx = indx + 1
	return None
def defineSpringElement(NAME, PART, node1Label,node2Label,DOFs,KK):
	if (node1Label in FreeNodes) or (node2Label in FreeNodes):
		return
	index1 = findNodeIndx(PART, node1Label)
	index2 = findNodeIndx(PART, node2Label)
	rgn1=regionToolset.Region(nodes=PART.nodes[index1:index1+1])
	rgn2=regionToolset.Region(nodes=PART.nodes[index2:index2+1])
	region=((rgn1, rgn2), )
	for i,dof in enumerate(DOFs):
		strr = '{0:s}_{1:1d}'.format(NAME,dof)
		PART.engineeringFeatures.TwoPointSpringDashpot(name=strr,\
			regionPairs=region,axis=FIXED_DOF,dof1=dof, dof2=dof,\
			orientation=None, springBehavior=ON,springStiffness=KK[i],\
			dashpotBehavior=OFF,dashpotCoefficient=0.0)
	"""
		PART.engineeringFeatures.TwoPointSpringDashpot(name=NAME+'{0:1d}'.format{dof},\
			regionPairs=region,axis=FIXED_DOF, dof1=dof, dof2=dof, orientation=None, \
			springBehavior=ON, springStiffness=KK[i],dashpotBehavior=OFF,\
			dashpotCoefficient=0.0)"""
def defineShellSection(MODEL, secName, matName, thickness):
	MODEL.HomogeneousShellSection(name=secName, 
        preIntegrate=OFF, material=matName, thicknessType=UNIFORM, 
        thickness=thickness, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
		
def defineConstrains2(INSTANCE_NAME,PART,MODEL, cnsName, node1, node2):
	rA = MODEL.rootAssembly
	PART.SetFromNodeLabels(name=cnsName+'_NODE1',\
		nodeLabels=(node1,), unsorted=True)
	PART.SetFromNodeLabels(name=cnsName+'_NODE2',\
		nodeLabels=(node2,), unsorted=True)
	rgn1=rA.sets['{0:s}.{1:s}'.format(INSTANCE_NAME,cnsName+'_NODE1')]
	rgn2=rA.sets['{0:s}.{1:s}'.format(INSTANCE_NAME,cnsName+'_NODE2')]
	MODEL.Coupling(name=cnsName, controlPoint=rgn1, 
        surface=rgn2, influenceRadius=10000000, couplingType=KINEMATIC, 
        alpha=0.0, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

####################################### DO NOT CHANGE ANYTHING ###############

#Load Datas
with open('GROUPS.pickle', 'rb') as handle:
	GROUPS = pickle.load(handle)
with open('ELEM_GROUPS.pickle', 'rb') as handle:
	ELEM_GROUPS = pickle.load(handle)
with open('NODE_GROUPS.pickle', 'rb') as handle:
    NODE_GROUPS = pickle.load(handle)
with open('MATERIALS.pickle', 'rb') as handle:
    MATERIALS = pickle.load(handle)
with open('Lumped_Mass.pickle', 'rb') as handle:
    Lumped_Mass = pickle.load(handle)
with open('Spring.pickle', 'rb') as handle:
    Spring = pickle.load(handle)
with open('TRUSS.pickle', 'rb') as handle:
    TRUSS = pickle.load(handle)
with open('TRUSS_ELIST.pickle', 'rb') as handle:
    TRUSS_ELIST = pickle.load(handle)
with open('BEAM_SECTIONS.pickle', 'rb') as handle:
    BEAM_SECTIONS = pickle.load(handle)
with open('BEAM_ORIENT.pickle', 'rb') as handle:
    BEAM_ORIENT = pickle.load(handle)
with open('BEAM_WITH_MISSING_ORIENTATION.pickle', 'rb') as handle:
    BEAM_WITH_MISSING_ORIENTATION = pickle.load(handle)
with open('BEAM_MAST_LIST.pickle', 'rb') as handle:
    BEAM_MAST_LIST = pickle.load(handle)
with open('UniquShellGroup.pickle', 'rb') as handle:
    UniquShellGroup = pickle.load(handle)
with open('SHELL_ELEMENT_LIST.pickle', 'rb') as handle:
    SHELL_ELEMENT_LIST = pickle.load(handle)
with open('MPC_Constrains.pickle', 'rb') as handle:
    MPC_Constrains = pickle.load(handle)
with open('LIAISON_SOLIDE.pickle', 'rb') as handle:
    LIAISON_SOLIDE = pickle.load(handle)
with open('LIAISON_ELEM.pickle', 'rb') as handle:
    LIAISON_ELEM = pickle.load(handle)
FreeNodes = np.loadtxt('InputsV2/PHASE_3_STAGE_1_INPUT_v2/02_Model/FreeNodes.dat', dtype=int, delimiter=',')
PivotErrorNodes = np.loadtxt('InputsV2/PHASE_3_STAGE_1_INPUT_v2/02_Model/AbaqusPivotErrorPointsFinal.dat', dtype=int, delimiter=',')
RB_FreeNodes = np.loadtxt('InputsV2/PHASE_3_STAGE_1_INPUT_v2/02_Model/RB_FreeNodes.dat', dtype=int, delimiter=',')
############ RE ASSIGN BASE ISOLATOR SPRING STIFFNESS
SPRINGLIST = ''
for gid, details in BaseIsolator.items():
	for eid in ELEM_GROUPS[gid]:
		Spring[eid]['KK'] = details	
		SPRINGLIST = SPRINGLIST + 'Spring{0:04d}\n'.format(eid)

with open('BaseIsolatorSpingList.dat', 'w') as f:
	f.write(SPRINGLIST)

INP_FILE = 'D:\CruassNPP\AERB_MODELS\Abaqus\geometry_without_group.inp'
MODELNAME = 'Model-1'

mdb.models[MODELNAME].PartFromInputFile(\
	inputFileName=INP_FILE)

MODEL = mdb.models[MODELNAME]
PART = MODEL.parts['PART-1']
# Define All Materials
for key, mat in MATERIALS.items():
	defineMaterial(str(key), MODEL, mat['RHO'],mat['E'],mat['NU'])

# Define All Inertial Masses and Associated Node Set
for key, inMass in Lumped_Mass.items():
	defineNodeSet(PART,str(key),tuple(inMass['nSet']))
	defineInertialMass(PART, str(key), str(key+'_inMass'), inMass['mass'])
# Define All Springs
for key, SPR in Spring.items():
	defineSpringElement('Spring{0:04d}'.format(key), PART,\
		SPR['nodes'][0],SPR['nodes'][1], SPR['DOFs'], SPR['KK'])
## All Truss Elements are BRI_1
II = 0
for key, TRS in TRUSS.items():
	nname = 'TRUSS{0:04d}'.format(II)
	defineElementSet(PART,nname+'TRUSS',tuple(TRS))
	defineTrussElement(MODEL, nname, 'BRI_1', float(key))
	defineSectionAssignment(PART, nname+'TRUSS', nname)
	II = II + 1

defineElementSet(PART,'ALL_TRUSS_ELIST',tuple(TRUSS_ELIST))

# Define Beam Elements and Its Orientation
II = 1
for key, BMS in BEAM_SECTIONS.items():
	for matName, elSet in BMS['material'].items():
		matName = str(matName)
		profName = 'BEAM{0:03d}'.format(II)+str(key) + '_Prof'
		secName = 'BEAM{0:03d}'.format(II)+str(key) + '_Sec'
		setName = 'BEAM{0:03d}'.format(II)+str(key) + '_Set'
		EE = MATERIALS[matName]['E']
		GG = 0.5*MATERIALS[matName]['E']/(1+MATERIALS[matName]['NU'])
		defineElementSet(PART,setName,tuple(BMS['elements']))
		para = BMS['parameters']
		if BMS['type'] == 'GENERAL':
			defineGeneralSection(MODEL,profName,secName, matName, para, EE,GG)
			defineSectionAssignment(PART, setName, secName)
		else:
			defineRectangularSection(MODEL,profName,secName,matName,para[0],para[1])
			defineSectionAssignment(PART, setName, secName)
	II = II + 1
for key, Ornt in BEAM_ORIENT.items():
	setName = str(key) + '_Set'
	defineElementSet(PART,setName,tuple(Ornt['elements']))
	defineBeamOrientation(PART, setName, tuple(Ornt['parameters']))
defineElementSet(PART,'ALL_BEAM_ELIST',tuple(BEAM_MAST_LIST))
defineElementSet(PART,'ALL_BEAM_WITHING_MISSING_ORIENTATION',\
	tuple(BEAM_WITH_MISSING_ORIENTATION))

II = 1
for key,shells in UniquShellGroup.items():
	for matName, elist in shells['material'].items():
		secName = 'UnqShell_{0:04d}'.format(II) + str(matName) + '_Sec'
		setName = 'UnqShell_{0:04d}'.format(II) + str(matName) + '_Set'
		defineElementSet(PART,setName,tuple(elist))
		defineShellSection(MODEL, secName, str(matName), float(key))
		defineSectionAssignment(PART, setName, secName)
	II = II + 1
defineElementSet(PART,'ALL_SHELL_ELIST',tuple(SHELL_ELEMENT_LIST))

## Defining Sets for Constrain Equation

#######################################
#################### V. V. Important ##############
# BEAM WITH MISSING ORIENTATION ASSING (1,1,1)
PART.assignBeamSectionOrientation(region=PART.sets['ALL_BEAM_WITHING_MISSING_ORIENTATION'],\
	method=N1_COSINES,n1=(1,1,1))
	
###############################################################################

INSTANCE_NAME = 'CRUAS_NPP'

MODEL.rootAssembly.Instance(name=INSTANCE_NAME, part=PART, dependent=ON) 
RA_I = MODEL.rootAssembly.instances[INSTANCE_NAME]

for id,const in MPC_Constrains.items():
	ConstrName = 'Constraint_{0:03d}'.format(id)
	SET1 = str(const[0].replace('\'',''))
	if not(SET1 in PART.sets.keys()):
		defineNodeSet(PART,SET1,tuple(NODE_GROUPS[SET1]))
	SET2 = str(const[1].replace('\'',''))
	DOF1 = DOF_MAP[str(const[2].replace('\'','')).strip()]
	DOF2 = DOF_MAP[str(const[3].replace('\'','')).strip()]
	if not(SET2 in PART.sets.keys()):
		defineNodeSet(PART,SET2,tuple(NODE_GROUPS[SET2]))
	MODEL.Equation(name=ConstrName,\
		terms=((float(const[4]),'{0:s}.{1:s}'.format(INSTANCE_NAME,SET1), DOF1),\
			(float(const[5]), '{0:s}.{1:s}'.format(INSTANCE_NAME,SET2), DOF2)))
########### LIAISON ELEM
for set in LIAISON_ELEM:
	defineElementSet(PART,str(set[0]),tuple(ELEM_GROUPS[str(set[0])]))
	defineNodeSet(PART,str(set[1]),tuple(NODE_GROUPS[str(set[1])]))
	PART.MeshSurfaceFromElsets(name=str(set[0])+'_Surf',\
		elementSetSeq=((PART.sets[str(set[0])],S1),))
	MODEL.Coupling(name=str(set[1])+'Constraint', controlPoint=RA_I.sets[str(set[1])], 
        surface=RA_I.surfaces[str(set[0])+'_Surf'],\
		influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
        alpha=0.0, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
	
# STICK MODELS
stick_con = np.loadtxt('StickModelConnectivity.dat', dtype=int, delimiter=',')
ii = 1
for con in stick_con:
	cnsName = 'StickCon{0:02d}'.format(ii)
	defineConstrains2(INSTANCE_NAME,PART,MODEL, cnsName, con[0], con[1])
	ii = ii + 1

# FIX ALL THE PIVOT ERROR POINTS
defineNodeSet(PART, 'PIVOT_ERROR_NODES',tuple(PivotErrorNodes))
defineNodeSet(PART, 'RB_FreeNodes',tuple(RB_FreeNodes))