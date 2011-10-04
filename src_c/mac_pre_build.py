#!/usr/bin/pythonw
"""
some trick to get it compiled on mac
"""
import os

# if we use all the files in the directory
# filenames = os.listdir(os.curdir) 

filenames = ["segmentation.c", "tools.c","image_processing.c", "trafo.c", "jw_ptv.c", "peakfitting.c", "rotation.c", "correspondences.c", "epi.c", "multimed.c", "ray_tracing.c","imgcoord.c","lsqadj.c", "orientation.c","sortgrid.c", "pointpos.c","intersect.c"]
newlines = []

# taking those which are C code:
# for filename in (f for f in filenames if f.endswith('.c')):

# or using only the given list
for filename in filenames:
	f = file(filename)
	for line in f:
		if 'ptv.h' in line:
			pass # print line
		else:
			newlines.append(line)
	
	
	f.close()
	
	
	
	
outfile = file('tmp.c','w')
outfile.write('#include "ptv.h"\n')
outfile.writelines(newlines)
outfile.close()
