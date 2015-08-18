import numpy as np
from repliscan import classProportion, hsvClass
import sys

def doThing(x,y,z):
	cp = classProportion(np.array([x,y,z]), emlSize=0.5)
	hp = hsvClass(np.array([x,y,z]))
	if np.any(cp != hp):
		print "NO Match:", (x,y,z), cp, hp
		#sys.exit()

doThing(2, 254, 126)
doThing(2, 254, 127)
doThing(3, 2, 2)
doThing(3, 2, 3)
sys.exit()

limit = 255
for x in xrange(1,limit):
	for y in xrange(1,limit):
		 for z in xrange(1,limit):
			doThing(x,y,z)
