import numpy as np
ANGLE_TOL = np.cos(np.deg2rad(85))
class PlainGeometry:
	def __init__(self, X0, n1):
		self.X0 = np.array(X0)
		self.n1 = np.array(n1)/np.linalg.norm(n1)
	def __call__(self,X):
		if np.linalg.norm(np.array(X)-self.X0) <= 0.000001:
			return False
		vv_u = (np.array(X)-self.X0)/np.linalg.norm(np.array(X)-self.X0)
		if abs(np.dot(vv_u,self.n1)) <= ANGLE_TOL:
			return True
		else:
			return False
	def LineSegmentIntersectionPoint(self,P0,P1):
		u = np.array(P1) - np.array(P0)
		vv_u = u/np.linalg.norm(u)
		w = np.array(P0) - self.X0
		if abs(np.dot(vv_u,self.n1)) <= ANGLE_TOL:
			return np.array([])
		scale = -1.0 * np.dot(self.n1,w)/np.dot(self.n1,u)
		if (scale <= 1.0) and (scale >= 0.0):
			return np.array(P0) + scale * u
		else:
			return np.array([])
if __name__ == '__main__':
	print('Testing of PlainGeometry Class..\n')
	x = np.linspace(-1,1,5)
	y = np.linspace(-1,1,5)
	xx, yy = np.meshgrid(x, y)
	zipped = zip(xx.reshape(25),yy.reshape(25))
	print('List of points used in testing..')
	print('\nPlain 01 - Point(0,0,0), Normal(1,1,0)')
	print('Plain 02 - Point(0,0,0), Normal(1,-1,0)')
	Plain01 = PlainGeometry([0,0,0],[1,1,0])
	Plain02 = PlainGeometry([0,0,0],[1,-1,0])
	for i,point in enumerate(zipped):
		print('{0:2d} {1:5.2f},{2:5.2f}:Plain 01 {3:6}, Plain 02 {4:6}'.format(i+1,point[0],point[1],\
			Plain01([point[0],point[1],0.0]),Plain02([point[0],point[1],0.0])))