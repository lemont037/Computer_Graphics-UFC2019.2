import math
import random

class Vector:
	x=0
	y=0
	z=0

	def __init__(self,x, y, z):
		self.x=x
		self.y=y
		self.z=z

	def printV(self):
		print(self.x, self.y, self.z)

def sum(v, u):
	return Vector(v.x+u.x, v.y+u.y, v.z+u.z)

def dotProduct(v, u):
	return v.x*u.x+v.y*u.y+v.z*u.z

def scalarProduct(v, k):
	return Vector(k*v.x, k*v.y, k*v.z)

def norm(v):
	return math.sqrt(dotProduct(v, v))

def crossProduct(v, u):
	return Vector(v.y*u.z-u.y*v.z, v.z*u.x-u.z*v.x, v.x*u.y-u.x*v.y)

def angle(v, u):
	a=norm(v)
	b=norm(u)
	if a==0 or b==0:
		print("Um dos vetores e nulo")
	else:
		return math.acos(dotProduct(u, v)/(a*b))



class Line:
	p=Vector(0, 0, 0)
	d=Vector(0, 0, 0)

	def __init__(self, p, d):
		if norm(d)==0:
			return self
		self.p=p
		self.d=scalarProduct(d, 1/norm(d))

class Plane:
	p=Vector(0, 0, 0)
	n=Vector(0, 0, 0)

	def __init__(self, p, v):
		if norm(v)==0:
			return None
		else:
			self.p=p
			self.n=scalarProduct(v, 1/norm(v))

	def isOnPlane(self, r):
		if dotProduct(self.n, r.d)==0:
			print("A reta e paralela ao plano")
			return None
		else:
			t=(dotProduct(sum(self.p, scalarProduct(r.p, -1)), self.n)/dotProduct(self.n, r.d))
			return sum(r.p, scalarProduct(r.d, t))

class Cylinder (Plane):
	d=0
	r=0

	def __init__(self, p, v, k, r):
		if norm(v)==0:
			return None
		self.p=p
		self.n=scalarProduct(v, 1/norm(v))
		self.d=k
		self.r=r

	def isOnCylinder(self, r):
		v=sum(sum(r.p, scalarProduct(self.p, -1)), scalarProduct(self.n, -1*dotProduct(sum(r.p, scalarProduct(self.p, -1)), self.n)))
		w=sum(r.d, scalarProduct(self.n, -1*dotProduct(r.d, self.n)))

		a=dotProduct(w, w)
		b=dotProduct(v, w)
		c=dotProduct(v, v)-math.pow(self.r, 2)
		delta=math.pow(b, 2)-(a*c)

		if a==0:
			v=self.isOnPlane(r)
			if norm(sum(v, scalarProduct(self.p,-1)))<=self.r:
				return v, sum(v, scalarProduct(self.n, self.d))
			else:
				return None, None

		elif delta<0:
			return None, None

		elif delta==0:
			t=(-1*b)/(a)
			x=sum(r.p, scalarProduct(r.d, t))
			z=dotProduct(sum(x, scalarProduct(self.p, -1)), self.n)
			if z<0.0 or z>self.d:
				return None, None
			else:
				return x, None

		else:
			t1=((-1*b)-math.sqrt(delta))/(a)
			t2=((-1*b)+math.sqrt(delta))/(a)
			x=sum(r.p, scalarProduct(r.d, t1))
			y=sum(r.p, scalarProduct(r.d, t2))
			z=dotProduct(sum(x, scalarProduct(self.p, -1)), self.n)
			if z<0 or z>self.d:
				x=None
			z=dotProduct(sum(y, scalarProduct(self.p, -1)), self.n)
			if z<0 or z>self.d:
				y=None
			if x==None or y==None:
				v=self.isOnPlane(r)
				if v==None:
					return x, y
				w=Plane(sum(self.p, scalarProduct(self.n, self.d)), self.n)
				z=w.isOnPlane(r)
				if x==None and y==None:
					if norm(sum(v, scalarProduct(self.p,-1)))<=self.r:
						x=v
					if norm(sum(z, scalarProduct(w.p,-1)))<=self.r:
						y=z
					return x, y
				if norm(sum(v, scalarProduct(self.p,-1)))<=self.r:
					if x==None:
						return y, v
					else:
						return x, v
				if norm(sum(z, scalarProduct(w.p,-1)))<=self.r:
					if x==None:
						return y, z
					else:
						return x, z
					
			else:
				return x, y

class Cone (Plane):
	h=0
	r=0
	v=Vector(0, 0, 0)

	def __init__(self, p, n, h, r):
		if norm(n)==0:
			return None
		self.p=p
		self.n=scalarProduct(n, 1/norm(n))
		self.h=h
		self.r=r
		self.v=sum(p, scalarProduct(self.n, self.h))

	def isOnCone(self, r):
		v=sum(self.v, scalarProduct(r.p, -1))
		t=math.pow(self.h, 2)/(math.pow(self.h, 2)+math.pow(self.r, 2))

		a=math.pow(dotProduct(r.d, self.n), 2)-t
		b=dotProduct(v, r.d)*t-dotProduct(v, self.n)*dotProduct(r.d, self.n)
		c=math.pow(dotProduct(self.n, v), 2)-dotProduct(v, v)*t
		delta=pow(b, 2)-a*c

		if a==0:
			t1=-1*c/(2*b)
			x=sum(r.p, scalarProduct(r.d, t1))
			if dotProduct(sum(self.v, scalarProduct(x, -1)), self.n)<0 or dotProduct(sum(self.v, scalarProduct(x, -1)), self.n)>self.h:
				x=None
			y=self.isOnPlane(r)
			if norm(sum(self.p, scalarProduct(y, -1)))<=self.r:
				return x, y
			else:
				return x, None

		elif delta<0:
			return None, None

		elif delta==0:
			t1=-1*b/a
			x=sum(r.p, scalarProduct(r.d, t1))
			if dotProduct(sum(self.v, scalarProduct(x, -1)), self.n)<0 or dotProduct(sum(self.v, scalarProduct(x, -1)), self.n)>self.h:
				return None, None
			else:
				return x, None

		else:
			t1=((-1*b)-math.sqrt(delta))/(a)
			t2=((-1*b)+math.sqrt(delta))/(a)
			x=sum(r.p, scalarProduct(r.d, t1))
			y=sum(r.p, scalarProduct(r.d, t2))
			z=dotProduct(sum(self.v, scalarProduct(x, -1)), self.n)
			if z<0 or z>self.h:
				x=None
			z=dotProduct(sum(self.v, scalarProduct(y, -1)), self.n)
			if z<0 or z>self.h:
				y=None
			if x==None and y==None:
				return x, y
			if x==None or y==None:
				z=self.isOnPlane(r)
				if norm(sum(self.p, scalarProduct(z, -1)))<=self.r:
					if x==None:
						return z, y
					else:
						return x, z
			else:
				return x, y

class Cube:
	p=Vector(0, 0, 0)
	l=0
	v1=Vector(0, 0, 0)
	v2=Vector(0, 0, 0)
	v3=Vector(0, 0, 0)
	v4=Vector(0, 0, 0)
	v5=Vector(0, 0, 0)
	v6=Vector(0, 0, 0)
	v7=Vector(0, 0, 0)
	v8=Vector(0, 0, 0)

	def __init__(self, p, l):
		self.p=p
		self.l=l/2
		self.v1=Vector(p.x+l/2, p.y+l/2, p.z+l/2)
		self.v2=Vector(p.x+l/2, p.y+l/2, p.z-l/2)
		self.v3=Vector(p.x+l/2, p.y-l/2, p.z+l/2)
		self.v4=Vector(p.x+l/2, p.y-l/2, p.z-l/2)
		self.v5=Vector(p.x-l/2, p.y+l/2, p.z+l/2)
		self.v6=Vector(p.x-l/2, p.y+l/2, p.z-l/2)
		self.v7=Vector(p.x-l/2, p.y-l/2, p.z+l/2)
		self.v8=Vector(p.x-l/2, p.y-l/2, p.z-l/2)

	def isOnCube(self, r):
		p=Plane(self.v1, crossProduct(sum(self.v1, scalarProduct(self.v2, -1)), sum(self.v1, scalarProduct(self.v3, -1))))
		x1=p.isOnPlane(r)
		p=Plane(self.v5, crossProduct(sum(self.v5, scalarProduct(self.v6, -1)), sum(self.v5, scalarProduct(self.v7, -1))))
		x2=p.isOnPlane(r)
		p=Plane(self.v1, crossProduct(sum(self.v1, scalarProduct(self.v2, -1)), sum(self.v1, scalarProduct(self.v5, -1))))
		x3=p.isOnPlane(r)
		p=Plane(self.v3, crossProduct(sum(self.v3, scalarProduct(self.v4, -1)), sum(self.v3, scalarProduct(self.v7, -1))))
		x4=p.isOnPlane(r)
		p=Plane(self.v1, crossProduct(sum(self.v1, scalarProduct(self.v3, -1)), sum(self.v1, scalarProduct(self.v5, -1))))
		x5=p.isOnPlane(r)
		p=Plane(self.v2, crossProduct(sum(self.v2, scalarProduct(self.v4, -1)), sum(self.v2, scalarProduct(self.v6, -1))))
		x6=p.isOnPlane(r)

		if x1!=None:
			if abs(x1.x-self.p.x)>self.l or abs(x1.y-self.p.y)>self.l or abs(x1.z-self.p.z)>self.l:
				x1=None
		if x2!=None:
			if abs(x2.x-self.p.x)>self.l or abs(x2.y-self.p.y)>self.l or abs(x2.z-self.p.z)>self.l:
				x2=None
		if x3!=None:
			if abs(x3.x-self.p.x)>self.l or abs(x3.y-self.p.y)>self.l or abs(x3.z-self.p.z)>self.l:
				x3=None
		if x4!=None:
			if abs(x4.x-self.p.x)>self.l or abs(x4.y-self.p.y)>self.l or abs(x4.z-self.p.z)>self.l:
				x4=None
		if x5!=None:
			if abs(x5.x-self.p.x)>self.l or abs(x5.y-self.p.y)>self.l or abs(x5.z-self.p.z)>self.l:
				x5=None
		if x6!=None:
			if abs(x6.x-self.p.x)>self.l or abs(x6.y-self.p.y)>self.l or abs(x6.z-self.p.z)>self.l:
				x6=None
		return x1, x2, x3, x4, x5, x6
		
class Panel:
    c=Vector(0, 0, 0)

    def __init__ (self, c, l, h, v) :
		    self.c = Vector(c.x - (l/2), c.y - (l/2), c.z)
		    self.l = l
		    self.h = h
		    self.v = v
    def pointOnPanel(self, H, V):
        return Vector(self.c.x + self.l/(2*self.h) + (H-1)*self.l/self.h, self.c.y + self.l/(2*self.v) + (V-1)*self.l/self.v, self.c.z)

#c=Cylinder(Vector(0, 0, 0), Vector(0, 0, 1), 5, 5)
#c=Cone(Vector(0, 0, 0), Vector(0, 0, 1), 5, 5)
k=Cube(Vector(0, 0, 0), 2)
l=Line(Vector(3*random.random(), 3*random.random(), 3*random.random()), Vector(3*random.random(), 3*random.random(), 3*random.random()))
#l=Line(Vector(0, 0, 0), Vector(1, 0, 0))

[a, b, c, d, e, f]=k.isOnCube(l)

if a!=None:
	a.printV()
else:
	print(a)
if b!=None:
	b.printV()
else:
	print(b)
if c!=None:
	c.printV()
else:
	print(c)
if d!=None:
	d.printV()
else:
	print(d)
if e!=None:
	e.printV()
else:
	print(e)
if f!=None:
	f.printV()
else:
	print(f)