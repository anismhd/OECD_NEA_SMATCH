"""
FE Descriptor Id definitions
____________________________

   11  Rod
   21  Linear beam
   22  Tapered beam
   23  Curved beam
   24  Parabolic beam
   31  Straight pipe
   32  Curved pipe
   41  Plane Stress Linear Triangle
   42  Plane Stress Parabolic Triangle
   43  Plane Stress Cubic Triangle
   44  Plane Stress Linear Quadrilateral
   45  Plane Stress Parabolic Quadrilateral
   46  Plane Strain Cubic Quadrilateral
   51  Plane Strain Linear Triangle
   52  Plane Strain Parabolic Triangle
   53  Plane Strain Cubic Triangle
   54  Plane Strain Linear Quadrilateral
   55  Plane Strain Parabolic Quadrilateral
   56  Plane Strain Cubic Quadrilateral
   61  Plate Linear Triangle
   62  Plate Parabolic Triangle
   63  Plate Cubic Triangle
   64  Plate Linear Quadrilateral
   65  Plate Parabolic Quadrilateral
   66  Plate Cubic Quadrilateral
   71  Membrane Linear Quadrilateral
   72  Membrane Parabolic Triangle
   73  Membrane Cubic Triangle
   74  Membrane Linear Triangle
   75  Membrane Parabolic Quadrilateral
   76  Membrane Cubic Quadrilateral
   81  Axisymetric Solid Linear Triangle
   82  Axisymetric Solid Parabolic Triangle
   84  Axisymetric Solid Linear Quadrilateral
   85  Axisymetric Solid Parabolic Quadrilateral
   91  Thin Shell Linear Triangle
   92  Thin Shell Parabolic Triangle
   93  Thin Shell Cubic Triangle
   94  Thin Shell Linear Quadrilateral
   95  Thin Shell Parabolic Quadrilateral
   96  Thin Shell Cubic Quadrilateral
   101 Thick Shell Linear Wedge
   102 Thick Shell Parabolic Wedge
   103 Thick Shell Cubic Wedge
   104 Thick Shell Linear Brick
   105 Thick Shell Parabolic Brick
   106 Thick Shell Cubic Brick
   111 Solid Linear Tetrahedron
   112 Solid Linear Wedge
   113 Solid Parabolic Wedge
   114 Solid Cubic Wedge
   115 Solid Linear Brick
   116 Solid Parabolic Brick
   117 Solid Cubic Brick
   118 Solid Parabolic Tetrahedron
   121 Rigid Bar
   122 Rigid Element
   136 Node To Node Translational Spring
   137 Node To Node Rotational Spring
   138 Node To Ground Translational Spring
   139 Node To Ground Rotational Spring
   141 Node To Node Damper
   142 Node To Gound Damper
   151 Node To Node Gap
   152 Node To Ground Gap
   161 Lumped Mass
   171 Axisymetric Linear Shell
   172 Axisymetric Parabolic Shell
   181 Constraint
   191 Plastic Cold Runner
   192 Plastic Hot Runner
   193 Plastic Water Line
   194 Plastic Fountain
   195 Plastic Baffle
   196 Plastic Rod Heater
   201 Linear node-to-node interface
   202 Linear edge-to-edge interface
   203 Parabolic edge-to-edge interface
   204 Linear face-to-face interface
   208 Parabolic face-to-face interface
   212 Linear axisymmetric interface
   213 Parabolic axisymmetric interface
   221 Linear rigid surface
   222 Parabolic rigin surface
   231 Axisymetric linear rigid surface
   232 Axisymentric parabolic rigid surface
"""

import logging
import numpy as np
import multiprocessing
import copy
import time
#import FEM
FLAG = '    -1'
# Universal file reader
nCPU = max(1, multiprocessing.cpu_count()-2)
UNV_FE_Descriptor_Id = {\
	11:'Rod',\
	21:'Linear beam',\
	22:'Tapered beam',\
	23:'Curved beam',\
	24:'Parabolic beam',\
	31:'Straight pipe',\
	32:'Curved pipe',\
	41:'Plane Stress Linear Triangle',\
	42:'Plane Stress Parabolic Triangle',\
	43:'Plane Stress Cubic Triangle',\
	44:'Plane Stress Linear Quadrilateral',\
	45:'Plane Stress Parabolic Quadrilateral',\
	46:'Plane Strain Cubic Quadrilateral',\
	51:'Plane Strain Linear Triangle',\
	52:'Plane Strain Parabolic Triangle',\
	53:'Plane Strain Cubic Triangle',\
	54:'Plane Strain Linear Quadrilateral',\
	55:'Plane Strain Parabolic Quadrilateral',\
	56:'Plane Strain Cubic Quadrilateral',\
	61:'Plate Linear Triangle',\
	62:'Plate Parabolic Triangle',\
	63:'Plate Cubic Triangle',\
	64:'Plate Linear Quadrilateral',\
	65:'Plate Parabolic Quadrilateral',\
	66:'Plate Cubic Quadrilateral',\
	71:'Membrane Linear Quadrilateral',\
	72:'Membrane Parabolic Triangle',\
	73:'Membrane Cubic Triangle',\
	74:'Membrane Linear Triangle',\
	75:'Membrane Parabolic Quadrilateral',\
	76:'Membrane Cubic Quadrilateral',\
	81:'Axisymetric Solid Linear Triangle',\
	82:'Axisymetric Solid Parabolic Triangle',\
	84:'Axisymetric Solid Linear Quadrilateral',\
	85:'Axisymetric Solid Parabolic Quadrilateral',\
	91:'Thin Shell Linear Triangle',\
	92:'Thin Shell Parabolic Triangle',\
	93:'Thin Shell Cubic Triangle',\
	94:'Thin Shell Linear Quadrilateral',\
	95:'Thin Shell Parabolic Quadrilateral',\
	96:'Thin Shell Cubic Quadrilateral',\
	101:'Thick Shell Linear Wedge',\
	102:'Thick Shell Parabolic Wedge',\
	103:'Thick Shell Cubic Wedge',\
	104:'Thick Shell Linear Brick',\
	105:'Thick Shell Parabolic Brick',\
	106:'Thick Shell Cubic Brick',\
	111:'Solid Linear Tetrahedron',\
	112:'Solid Linear Wedge',\
	113:'Solid Parabolic Wedge',\
	114:'Solid Cubic Wedge',\
	115:'Solid Linear Brick',\
	116:'Solid Parabolic Brick',\
	117:'Solid Cubic Brick',\
	118:'Solid Parabolic Tetrahedron',\
	121:'Rigid Bar',\
	122:'Rigid Element',\
	136:'Node To Node Translational Spring',\
	137:'Node To Node Rotational Spring',\
	138:'Node To Ground Translational Spring',\
	139:'Node To Ground Rotational Spring',\
	141:'Node To Node Damper',\
	142:'Node To Gound Damper',\
	151:'Node To Node Gap',\
	152:'Node To Ground Gap',\
	161:'Lumped Mass',\
	171:'Axisymetric Linear Shell',\
	172:'Axisymetric Parabolic Shell',\
	181:'Constraint',\
	191:'Plastic Cold Runner',\
	192:'Plastic Hot Runner',\
	193:'Plastic Water Line',\
	194:'Plastic Fountain',\
	195:'Plastic Baffle',\
	196:'Plastic Rod Heater',\
	201:'Linear node-to-node interface',\
	202:'Linear edge-to-edge interface',\
	203:'Parabolic edge-to-edge interface',\
	204:'Linear face-to-face interface',\
	208:'Parabolic face-to-face interface',\
	212:'Linear axisymmetric interface',\
	213:'Parabolic axisymmetric interface',\
	221:'Linear rigid surface',\
	222:'Parabolic rigin surface',\
	231:'Axisymetric linear rigid surface',\
	232:'Axisymentric parabolic rigid surface',\
	}

class UNV():
	def __init__(self, filename):
		self.filename = filename
		self.file = open(self.filename, 'r')
		self.fem = {'nodes':{},'elements':{},'groups':{}}
		self.sections = {}
		self.datasetsIds = [2411, 2412, 2477, 2467]
		self.datasetsHandlers = [UNV2411Reader, UNV2412Reader, UNV2477Reader, UNV2477Reader]
		self.read()
	def scanfile(self):
		inside = False
		self.lines = self.file.readlines()
		for i, line in enumerate(self.lines):
			if line.startswith(FLAG):
				inside = not(inside)
				if inside:
					gid = int(self.lines[i+1][:6])
					self.sections[gid] = [i+2]
				else:
					self.sections[gid].append(i)
	def print_sections(self,gid):
		print(self.lines[self.sections[gid][0]:self.sections[gid][1]])
	def UNVReader(self,sectionId, lines, fem):
		func = self.datasetsHandlers[self.datasetsIds.index(sectionId)]
		func(lines, fem)
	def read(self):
		self.scanfile()
		for sectionId, indexs in self.sections.items():
			start_time = time.time()  # Start time
			if sectionId in [151,775, 164,2420]:
				continue
			lines = self.lines[indexs[0]:indexs[1]]
			func = self.datasetsHandlers[self.datasetsIds.index(sectionId)]
			func(lines, self.fem)
			end_time = time.time()  # End time
			elapsed_time = end_time - start_time
			print(f"Execution Time: {elapsed_time} seconds for {sectionId}")
	def readp(self):
		self.scanfile()
		new_sections = {}
		for sectionId, indexs in self.sections.items():
			if sectionId in [151,775, 164,2420]:
				continue
			new_sections[sectionId] = copy.deepcopy(self.lines[indexs[0]:indexs[1]])
		with multiprocessing.Manager() as manager:
			temp_fem = manager.dict()
			temp_fem = {'nodes':{},'elements':{},'groups':{}}
			with multiprocessing.Pool(processes=nCPU) as pool:
				# Use a lambda to fix the first argument while iterating over the second
				pool.starmap(self.UNVReader,\
					[(sectionId,lines,temp_fem) for sectionId, lines in new_sections.items()])
			"""
			for sectionId, lines in new_sections.items():
				self.UNVReader(sectionId,lines,temp_fem)
			
			with multiprocessing.Pool(processes=nCPU) as pool:
				
			for sectionId, indexs in self.new_sections.items():
				lines = self.lines[indexs[0]:indexs[1]]
				lines_iter = iter(lines)
				func = self.datasetsHandlers[self.datasetsIds.index(sectionId)]
				for line in lines_iter:
					func(line, temp_fem)
				with multiprocessing.Pool(processes=nCPU) as pool:
					# Use a lambda to fix the first argument while iterating over the second
					pool.starmap(func,\
						[(line, temp_fem) for line in lines_iter])
					#pool.map(lambda line: func(line, temp_fem), lines)"""
			self.fem = copy.deepcopy(temp_fem)
			"""
			if (sectionId in self.datasetsIds):
				self.file.seek(offset)
				func = self.datasetsHandlers[self.datasetsIds.index(sectionId)]
				self.fem = func(self.file, self.fem)
				self.file.close()"""


def Line2Float(line):
	return [float(x) for x in line.split()]

def Line2Int(line):
	return [int(x) for x in line.split()]

	
def UNV2411Reader(lines, fem):
	lines_iter = iter(lines)
	for line in lines_iter:
		dataline1 = Line2Int(line)
		fem['nodes'][dataline1[0]] = Line2Float(next(lines_iter))

def UNV2412Reader(lines, fem):
	lines_iter = iter(lines)
	for line in lines_iter:
		dataline1 = Line2Int(line)
		eid = dataline1[0]
		if not(dataline1[1] in fem['elements'].keys()):
			fem['elements'][dataline1[1]] = {}
		if beam_description_id_check(dataline1[1]):
			dataline2 = Line2Int(next(lines_iter))
			dataline3 = Line2Int(next(lines_iter))
			fem['elements'][dataline1[1]][eid] = dataline2 + dataline3
		else:
			fem['elements'][dataline1[1]][eid] = Line2Int(next(lines_iter))

def UNV2477Reader(lines, fem):
	lines_iter = iter(lines)
	for line in lines_iter:
		dataline1 = Line2Int(line)
		group_no = dataline1[0]
		no_lines = int(np.ceil(dataline1[7]/2))
		group_name = next(lines_iter).strip()
		group_data = []
		for i in range(no_lines):
			group_data = group_data + Line2Int(next(lines_iter))
		group_data = np.resize(np.array(group_data,dtype=int),(dataline1[7],4))
		fem['groups'][group_no] = {'group_name':group_name,\
			'group_data':group_data}
			
def UNV2477Readerp(lines, fem):
	lines_iter = iter(lines)
	for line in lines_iter:
		dataline1 = Line2Int(line)
		group_no = dataline1[0]
		no_lines = int(np.ceil(dataline1[7]/2))
		group_name = next(lines_iter).strip()
		group_data = []
		for i in range(no_lines):
			group_data = group_data + Line2Int(next(lines_iter))
		group_data = np.resize(np.array(group_data,dtype=int),(dataline1[7],4))
		fem['groups'][group_no] = {'group_name':group_name,\
			'group_data':group_data}
			
def UNV2411Reader0(line, fem):
	print(line)
	dataline1 = Line2Int(line)
	fem['nodes'][dataline1[0]] = Line2Float(next(lines_iter))

def UNV2412Reader0(line, fem):
	print(line)
	dataline1 = Line2Int(line)
	eid = dataline1[0]
	if not(dataline1[1] in fem['elements'].keys()):
		fem['elements'][dataline1[1]] = {}
	if beam_description_id_check(dataline1[1]):
		dataline2 = Line2Int(next(lines_iter))
		dataline3 = Line2Int(next(lines_iter))
		fem['elements'][dataline1[1]][eid] = dataline2 + dataline3
	else:
		fem['elements'][dataline1[1]][eid] = Line2Int(next(lines_iter))

def UNV2477Reader0(line, fem):
	print(line)
	dataline1 = Line2Int(line)
	group_no = dataline1[0]
	no_lines = int(np.ceil(dataline1[7]/2))
	group_name = next(lines_iter).strip()
	group_data = []
	for i in range(no_lines):
		group_data = group_data + Line2Int(next(lines_iter))
	group_data = np.resize(np.array(group_data,dtype=int),(dataline1[7],4))
	fem['groups'][group_no] = {'group_name':group_name,\
		'group_data':group_data}
			
def beam_description_id_check(id):
	if id in [21,22,23,24]:
		return True
	else:
		return False