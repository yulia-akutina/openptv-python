class Tracking():
	""" Tracking class defines external tracking addon for pyptv
	User needs to implement the following functions:
		do_tracking(self)
		do_back_tracking(self)
	Connection to C ptv module is given via self.ptv and provided by pyptv software
	Connection to active parameters is given via self.exp1 and provided by pyptv software.
	User responsibility is to read necessary files, make the calculations and write the files back.
	"""
	
	def __init__(self,ptv=None,exp1=None):
		self.ptv=ptv
		self.exp1=exp1
		# Do your initialization here

	def do_tracking(self):
		""" this function is callback for "tracking without display"
		"""
		print "inside denis_ext_tracker"
		lmax_track,ymin_track,ymax_track, seq_first, seq_last=self.ptv.py_trackcorr_init()
		print lmax_track,ymin_track,ymax_track, seq_first, seq_last
		for step in range (seq_first, seq_last):
			print step
			self.ptv.py_trackcorr_loop(step, lmax_track, ymin_track, ymax_track,display=0)
   			#finalize tracking
		self.ptv.py_trackcorr_finish(step+1)
		
		
	def do_back_tracking(self):
		""" this function is callback for "tracking back"
		"""
		# do your back_tracking stuff here
		print "inside custom back tracking"
