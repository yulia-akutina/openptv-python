/* Define classes for handling parameters (storing, reading, writing) so that
   all of libptv has a single point of entry rather than different formats in 
   the same code.
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

/* Length of basenames. Smaller enough than the buffer lengths used by the
   tracking framebuffer, that suffixes can be added. */
#define SEQ_FNAME_MAX_LEN 240

typedef struct {
    char **img_base_name;
    int first, last;
} sequence_par;

sequence_par* read_sequence_par(char *filename);

typedef struct
{
    double  dacc, dangle, dvxmax, dvxmin;
    double dvymax, dvymin, dvzmax, dvzmin;
    int dsumg, dn, dnx, dny, add;
} track_par;

track_par* read_track_par(char *filename);
int compare_track_par(track_par *t1, track_par *t2);

typedef struct {
    double  X_lay[2], Zmin_lay[2], Zmax_lay[2];
} volume_par;

volume_par* read_volume_par(char *filename);
int compare_volume_par(volume_par *v1, volume_par *v2);

#endif

