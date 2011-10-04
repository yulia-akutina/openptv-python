""" PyPTV_BATCH is the script for the 3D-PTV (http://ptv.origo.ethz.ch) written in 
Python/Enthought Traits GUI/Numpy/Chaco

Example:
>> python pyptv_batch.py experiments/exp1 


"""


from scipy.misc import imread
import os
import sys
import numpy as np

# project specific inputs
import parameters as par
import ptv1 as ptv
import general


# directory from which we run the software
cwd = os.getcwd()


# import pdb; pdb.set_trace()

if len(sys.argv) < 2:
	error("Wrong number of inputs, usage: python pyptv_batch.py experiments/exp1")
	
software_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
print 'software_path=', software_path

try:
	os.chdir(software_path)
except:
	print("Error in instalation or software path")

src_path = os.path.join(os.path.split(os.path.abspath(os.getcwd()))[0],'src_c')
sys.path.append(src_path)



exp_path = os.path.abspath(sys.argv[1])
print 'exp_path=', exp_path


try:
	os.chdir(exp_path)
	print(os.getcwd())
except:
	print('Wrong experimental directory %s' % exp_path)


def sequence_tracking(n_img):
	# get following variables from the parameters:
	# n_camera, seq_first, seq_last, base_name
	sequenceParams = par.SequenceParams(n_img, path = par.temp_path)
	sequenceParams.read()
	(base_name, seq_first, seq_last) = (sequenceParams.base_name, sequenceParams.first, sequenceParams.last)

	print ("Starting sequence action")
	
	ptv.py_sequence_init(0)
	stepshake=ptv.py_get_from_sequence_init()
	if not stepshake:
		stepshake=1
	print stepshake
	temp_img=np.array([],dtype=np.ubyte)
	for i in range(seq_first,seq_last+1,stepshake):
		if i<10:
			seq_ch="%01d" % i
		elif i<100:
			seq_ch="%02d" % i
		else:
			seq_ch="%03d" % i
		for j in range (n_img):
			img_name=base_name[j]+seq_ch
			print ("Setting image: ",img_name)
			try:
				temp_img=imread(img_name).astype(np.ubyte)
			except:
				print "Error reading file"
				   
			ptv.py_set_img(temp_img,j)
		
		ptv.py_sequence_loop(0,i)

# forward tracking
	lmax_track,ymin_track,ymax_track, seq_first, seq_last=ptv.py_trackcorr_init()
	print lmax_track,ymin_track,ymax_track, seq_first, seq_last
	for step in range (seq_first, seq_last):
		print step
		ptv.py_trackcorr_loop(step, lmax_track, ymin_track, ymax_track,display=0)
   
	ptv.py_trackcorr_finish(step+1)
	print "tracking without display finished"

# backward tracking
	ptv.py_trackback_c()




def run_batch():
# 	import pdb; pdb.set_trace()
	ptv.py_init_proc_c()
	ptv.py_start_proc_c() # or ptv.py_init_proc_c()?
	ptvParams = par.PtvParams(path = par.temp_path)
	ptvParams.read()
	(n_img, img_name, img_cal, hp_flag, allCam_flag, tiff_flag, imx, imy, pix_x, pix_y, chfield, mmp_n1, mmp_n2, mmp_n3, mmp_d) = \
		(ptvParams.n_img, ptvParams.img_name, ptvParams.img_cal, ptvParams.hp_flag, ptvParams.allCam_flag, ptvParams.tiff_flag, \
		ptvParams.imx, ptvParams.imy, ptvParams.pix_x, ptvParams.pix_y, ptvParams.chfield, ptvParams.mmp_n1, ptvParams.mmp_n2, ptvParams.mmp_n3, ptvParams.mmp_d)
#

	sequence_tracking(n_img)


if __name__ == '__main__':
	try:
		run_batch()
	except:
		print("something wrong with the software or folder")
		general.printException()
