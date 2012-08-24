"""
some trick to get it compiled on mac
"""
import os
import shutil

cwd = os.getcwd()

for line in file('setup.py'):
	if line.strip().startswith('ext_modules'):
		lst = line

# import pdb; pdb.set_trace()

filenames = lst.partition('"ptv1.pyx",')[-1].lstrip().strip('],\n').split(',')
print filenames
newlines = []

# taking those which are C code:
# for filename in (f for f in filenames if f.endswith('.c')):

src_path = os.path.join(os.path.split(os.path.abspath(cwd))[0],'src_c')
os.chdir(src_path)


# or using only the given list
for filename in filenames:
	print filename.strip().strip('"')
	f = file(filename.strip().strip('"'))
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

print os.getcwd()

os.system('python setup_mac.py build_ext --inplace')
os.remove('tmp.c')
