import numpy as np
import sympy
import io

# Base Isolator Spring Stiffness 
STIFF_X = 0.15E+08
STIFF_Y = 0.15E+08
STIFF_Z = 0.45E+08

# Material Char. File Path

MAT_FILE = 'InputsV2/PHASE_3_STAGE_1_INPUT_v2/02_Model/Model_characteristics_v2.txt'

try:
	MAT_FILE_ID = open(MAT_FILE,  encoding="utf8")
	MAT_FILE_LINE = MAT_FILE_ID.readlines()
except:
	pass

class material_section_read():
	def __init__(self):
		self.material = {}
		self.section_details = {'beam':{},'M_TR_D_N':{},'M_T_N':{},'K_T_D_L':{},\
			'COQUE':{},'K_TR_D_L':{},'SHELL':{},'ROD':{}}
		self.material_assignment = {}
		self.modelization = {'DKT':[], 'DIS_T':[], 'DIS_TR':[],'POU_D_T':[],'ROD':[]}
		self.MPC = {} # Multipoint constraint # GROUP_NO_1;GROUP_NO_2;DDL_1;DDL_2;COEF_MULT_1;COEF_MULT_2
		self.LIAISON_SOLIDE = []
		self.beam_orientation = {}
		self.readvariables01()
		self.readmat01()
		self.readvariables02()
		self.readmat02()
		self.readvariables03()
		self.readmat03()
		self.readmatassign()
		self.readmodalization01()
		self.readmpc()
		self.readbeamsection()

	def readvariables01(self):
		for line in MAT_FILE_LINE[22:41]:
			vv = line.strip().split('=')
			setattr(self, vv[0].strip(), eval(vv[1].strip()))
			globals()['{0:s}'.format(vv[0].strip())] = eval(vv[1].strip()) 
	def readmat01(self):
		for line in MAT_FILE_LINE[45:112]:
			linedata = line.strip().split(';')
			self.material[linedata[0].strip()] = {\
				'TYPE OF BEHAVIOR':linedata[1].strip(),\
				'E':eval(linedata[2].strip()),\
				'RHO':eval(linedata[3].strip()),\
				'NU':eval(linedata[4].strip()) }
	def readvariables02(self):
		for line in MAT_FILE_LINE[120:273]:
			vv = line.strip().split('=')
			if vv[0].strip() in dir(self):
				print('You are rewriting {0:s}'.format(vv[0].strip()))
			setattr(self, vv[0].strip(), eval(vv[1].strip()))
			globals()['{0:s}'.format(vv[0].strip())] = eval(vv[1].strip()) 
	def readmat02(self):
		for line in MAT_FILE_LINE[276:477]:
			linedata = line.strip().split(';')
			if linedata[0].strip() in self.material.keys():
				print('You are rewriting {0:s}'.format(linedata[0].strip()))
			self.material[linedata[0].strip()] = {\
				'TYPE OF BEHAVIOR':linedata[1].strip(),\
				'E':eval(linedata[2].strip()),\
				'RHO':eval(linedata[4].strip()),\
				'NU':eval(linedata[3].strip()) }

	def readvariables03(self):
		for line in MAT_FILE_LINE[489:512]:
			if not line.strip():
				continue
			line = line.replace(';','')
			vv = line.strip().split('=')
			if vv[0].strip() in dir(self):
				print('You are rewriting {0:s}'.format(vv[0].strip()))
			setattr(self,vv[0].strip(),eval(vv[1].strip()))
			globals()['{0:s}'.format(vv[0].strip())] = eval(vv[1].strip())
	def readmat03(self):
		for line in MAT_FILE_LINE[514:725]:
			if not line.strip():
				continue
			linedata = line.strip().split(';')
			if linedata[0].strip() in self.material.keys():
				print('You are rewriting {0:s}'.format(linedata[0].strip()))
			self.material[linedata[0].strip()] = {\
				'TYPE OF BEHAVIOR':'ELAS',\
				'E':eval(linedata[1].strip()),\
				'RHO':eval(linedata[3].strip()),\
				'NU':eval(linedata[2].strip()) }

	def readmatassign(self):
		for line in MAT_FILE_LINE[739:2000]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			linedata = line.strip().split(';')
			if len(linedata) < 2:
				continue
			self.material_assignment[linedata[0].strip()] = linedata[1].strip()
	def readmodalization01(self):
		for line in MAT_FILE_LINE[2023:2414]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			linedata = line.strip().split(';')
			if len(linedata) < 2:
				continue
			self.modelization[linedata[1].strip()].append(linedata[0].strip())
		for line in MAT_FILE_LINE[2446:2504]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			linedata = line.strip().split(';')
			if len(linedata) < 2:
				continue
			self.modelization[linedata[1].strip()].append(linedata[0].strip())
	def readmpc(self):
		self.MPC[1] = MAT_FILE_LINE[2419].strip().split(';')
		for line in MAT_FILE_LINE[2513:2526]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			self.MPC[list(self.MPC.keys())[-1]+1] = line.strip().split(';')
		for line in MAT_FILE_LINE[3110:3112]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			self.MPC[list(self.MPC.keys())[-1]+1] = line.strip().split(';')
	def readbeamsection(self):
		for line in MAT_FILE_LINE[2573:]:
			if not line.strip():
				continue
			if line[0] == '#':
				continue
			if ('MESH GROUP' in line) or ('MESH_GROUP' in line):
				continue
			line = line.replace('\'','')
			linedata = line.strip().split(';')
			if 'BEAM' in linedata[0]:
				if ('GENERAL' in line) or ('GENERALE' in line):
					grp_name, dictnr = beam_general_section_reader(line)
					self.section_details['beam'][grp_name] = dictnr
				elif 'RECTANGLE' in line:
					grp_name, dictnr = beam_rect_section_reader(line)
					self.section_details['beam'][grp_name] = dictnr
				else:
					print(line.strip())
			elif 'ROD' in linedata[0]:
				if ('GENERAL' in line) or ('GENERALE' in line):
					linedata = line.strip().split(';')
					grp_name = linedata[1].strip().replace('(','').replace(')','').replace(',','')
					self.section_details['ROD'][grp_name] = {'type':'GENERAL', 'parameters':float(linedata[-1])}
				elif 'RECTANGLE' in line:
					linedata = line.strip().split(';')
					grp_name = linedata[1].strip().replace('(','').replace(')','').replace(',','')
					parameters = [float(dd) for dd in linedata[-1].replace('(','').replace(')','').split(',')]
					self.section_details['ROD'][grp_name] = {'type':'RECTANGLE',\
						 'parameters':parameters}
				elif 'CIRCLE' in line:
					linedata = line.strip().split(';')
					grp_name = linedata[1].strip().replace('(','').replace(')','').replace(',','')
					parameters = [float(dd) for dd in linedata[-1].replace('(','').replace(')','').split(',')]
					self.section_details['ROD'][grp_name] = {'type':'CIRCLE',\
						 'parameters':parameters}
				else:
					print(line.strip())
			elif 'M_TR_D_N' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				for vv in linedata[1:]:
					values = vv.replace('(','').replace(')','').split(',')
					values = list(filter(None, values))
					if len(values) > 2:
						parameters = [float(pp) for pp in values]
				self.section_details['M_TR_D_N'][linedata[2].strip()] = {'type':'GENERAL',\
						 'parameters':parameters}
			elif 'M_T_N' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				for vv in linedata[1:]:
					values = vv.replace('(','').replace(')','').split(',')
					values = list(filter(None, values))
					if len(values) > 2:
						parameters = [eval(pp) for pp in values]
				self.section_details['M_T_N'][linedata[2].strip()] = {'type':'GENERAL',\
						 'parameters':parameters}
			elif 'K_T_D_L' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				for vv in linedata[1:]:
					values = vv.replace('(','').replace(')','').split(',')
					values = list(filter(None, values))
					if len(values) > 2:
						parameters = [eval(pp) for pp in values]
				self.section_details['K_T_D_L'][linedata[2].strip()] = {'type':'GENERAL',\
						 'parameters':parameters}
			elif 'K_TR_D_L' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				for vv in linedata[1:]:
					values = vv.replace('(','').replace(')','').split(',')
					values = list(filter(None, values))
					if len(values) > 2:
						parameters = [eval(pp) for pp in values]
				self.section_details['K_TR_D_L'][linedata[2].strip()] = {'type':'GENERAL',\
						 'parameters':parameters}
			elif 'COQUE' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				parameters = []
				for vv in linedata[1:]:
					try:
						parameters.append(eval(vv))
					except:
						grpid = vv.strip()
				self.section_details['COQUE'][grpid] = parameters
			elif 'SHELL' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				parameters = []
				for vv in linedata[1:]:
					try:
						parameters.append(eval(vv))
					except:
						grpid = vv.strip()
				self.section_details['SHELL'][grpid] = parameters
			elif 'VECT_Y' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				for vv in linedata:
					if 'VECT_Y' in vv:
						continue
					elif len(vv.strip().split(',')) > 2:
						values = vv.replace('(','').replace(')','').split(',')
						values = list(filter(None, values))
					else:
						grd_name = vv.strip()
				self.beam_orientation[grd_name] = values
			elif 'LIAISON_SOLIDE' in line:
				line = line.replace('\'','')
				linedata = line.strip().split(';')
				self.LIAISON_SOLIDE.append(linedata[1].strip())
			else:
				print(line.strip())
				
	def check_attribute(self, name):
		print(name)

def beam_rect_section_reader(line):
	line = line.replace('\'','')
	linedata = line.strip().split(';')
	for vv in linedata[1:]:
		values = vv.replace('(','').replace(')','').split(',')
		values = list(filter(None, values))
		if len(values) > 1:
			values = vv.replace('(','').replace(')','').split(',')
			try:
				parameters = [float(pp) for pp in values[:2]]
			except:
				continue
		elif 'RECTANGLE' in vv:
			continue
		else:
			grp_name = vv.strip().replace('(','').replace(')','').replace(',','')
	return grp_name, {'type':'RECTANGLE', 'parameters':parameters}
def beam_general_section_reader(line):
	line = line.replace('\'','')
	linedata = line.strip().split(';')
	for vv in linedata[1:]:
		if len(vv.split(',')) > 3:
			values = vv.replace('(','').replace(')','').split(',')
			try:
				parameters = [float(pp) for pp in values]
			except:
				continue
		elif ('GENERAL' in vv) or ('GENERALE' in vv):
			continue
		else:
			grp_name = vv.strip()
	return grp_name, {'type':'GENERAL', 'parameters':parameters}
if __name__ == '__main__':
	new_mat = material_section_read()
	print(new_mat.LIAISON_SOLIDE)