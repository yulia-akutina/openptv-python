    
from libc.stdlib cimport free

cdef class TrackingRun:
    def get_sequence_range(TrackingRun self):
        return self.tr[0].seq_par[0].first, self.tr[0].seq_par[0].last

    def __dealloc__(TrackingRun self):
        tr_free(self.tr)
        free(self.tr)

