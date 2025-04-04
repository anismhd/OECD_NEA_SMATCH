import regionToolset
import connectorBehavior
from abaqusConstants import *
import numpy as np

#FileList = ['NPCIL_at_13p5_u.ACC','NPCIL_at_13p5_v.ACC','NPCIL_z.ACC','Recorded_at_13p5_u.ACC','Recorded_at_13p5_v.ACC','Recorded_z.ACC']
FileList = ['NPCIL_u.ACC','NPCIL_v.ACC','Recorded_u.ACC','Recorded_v.ACC']

for file in FileList:
	data = np.loadtxt('InputSignals/' + file)[1:]
	tt = []
	for i,val in enumerate(data):
		tt.append((i*0.005,val))
	mdb.models['Model-1'].TabularAmplitude(name=file[:-4], timeSpan=STEP, 
			smooth=SOLVER_DEFAULT, data=tuple(tt))