/*
Implementation the tracking frame buffers. Each frame holds target information
for all cameras, correspondence information and path links information.

Contains functions for dealing with the frame content, such as comparing and
reading targets, correspondences, etc.
*/

#include <string.h>
#include <stdio.h>
#include "tracking_frame_buf.h"

/* Check that target t1 is equal to target t2, i.e. all their fields are equal.
 * 
 * Arguments:
 * target *t1, *t2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_targets(target *t1, target *t2) {
    return (\
        (t1->pnr == t2->pnr) && (t1->x == t2->x) && (t1->y == t2->y) && \
        (t1->n == t2->n) && (t1->nx == t2->nx) && (t1->ny == t2->ny) && \
        (t1->sumg == t2->sumg) && (t1->tnr == t2->tnr));
}

/* Reads targets from a file. The number of targets is read from the first
 * line, then each line is one target.
 * 
 * Arguments:
 * target buffer[] - an allocated array of target structs to fill in from
 *   files.
 * char* file_base - base name of the files to read, to which a frame number
 *   and the suffix '_targets' is added.
 * int frame_num - number of frame to add to file_base. A value of 0 or less
 *   means that no frame number should be added.
 * 
 * Returns:
 * the number of targets found in the file, or 0 if an error occurred.
*/

int read_targets(target buffer[], char* file_base, int frame_num) {
    FILE *FILEIN;
    int	tix, num_targets, scanf_ok;
    char filein[STR_MAX_LEN + 1];

    if (frame_num > 0) {
        sprintf(filein, "%s%04d%s", file_base, frame_num, "_targets");
    } else {
        strncpy(filein, file_base, STR_MAX_LEN);
        strncat(filein, "_targets", STR_MAX_LEN);
    }
    
    FILEIN = fopen (filein, "r");
    if (! FILEIN) {
        printf("Can't open ascii file: %s\n", filein);
        return 0;
    }

    if (fscanf(FILEIN, "%d\n", &num_targets) == 0) {
        printf("Bad format for file: %s\n", filein);
        return 0;
    }
    for (tix = 0; tix < num_targets; tix++)	{
	  scanf_ok = fscanf (FILEIN, "%4d %lf %lf %d %d %d %d %d\n",
		  &(buffer[tix].pnr),  &(buffer[tix].x),
		  &(buffer[tix].y),    &(buffer[tix].n),
		  &(buffer[tix].nx),   &(buffer[tix].ny),
		  &(buffer[tix].sumg), &(buffer[tix].tnr) );
      
      if (scanf_ok == 0) {
        printf("Bad format for file: %s\n", filein);
        return 0;
      }
	}
    fclose (FILEIN);

	return num_targets;
}

/* Check that two correspondence structs are equal, i.e. all their fields are 
 * equal.
 * 
 * Arguments:
 * corres *c1, *c2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_corres(corres *c1, corres *c2) {
    return ((c1->nr == c2->nr) && (c1->p[0] == c2->p[0]) && \
        (c1->p[1] == c2->p[1]) && (c1->p[2] == c2->p[2]) && \
        (c1->p[3] == c2->p[3]));
}

/* Check that two path info (P) structs are equal, i.e. all their fields are 
 * equal.
 * 
 * Arguments:
 * P *p1, *p2 - the target structures to compare.
 * 
 * Returns:
 * true if all fields are equal, false otherwise.
*/
int compare_path_info(P *p1, P *p2) {
    int iter;
    if (!((p1->prev == p2->prev) && (p1->next == p2->next) && \
        (p1->prio == p2->prio) && (p1->finaldecis == p2->finaldecis) && \
        (p1->inlist == p2->inlist)))
        return 0;

    for (iter = 0; iter < 3; iter++) {
        if (p1->x[iter] != p2->x[iter]) return 0;
    }
    for (iter = 0; iter < POSI; iter++) {
        if (p1->decis[iter] != p2->decis[iter]) return 0;
        if (p1->linkdecis[iter] != p2->linkdecis[iter]) return 0;
    }
    return 1;
}

/* Reads rt_is files. these files contain both the path info and the 
 * information on correspondence between targets on the different images.
 * Sets fields not in file to default values.
 * 
 * Arguments:
 * corres *cor_buf - a buffer of corres structs to fill in from the file.
 * P *path_buf - same for path info structures.
 * char* file_base - base name of the files to read, to which a frame number
 *   is added. Without separator.
 * int frame_num - number of frame to add to file_base. A value of 0 or less
 *   means that no frame number should be added. The '.' separator is added
 * between the name and the frame number.
 * 
 * Returns:
 * The number of points read for this frame. 0 on failure.
*/
int read_path_frame(corres *cor_buf, P *path_buf, \
    char *file_base, int frame_num) {

    FILE *filein;
    char fname[STR_MAX_LEN];
    int read_res = 0, targets = 0, alt_link = 0;
    
    sprintf(fname, "%s.%d", file_base, frame_num);
    filein = fopen (fname, "r");
    if (!filein) {
        /* Keeping the printf until we have proper logging. */
        printf("Can't open ascii file: %s\n", fname);
        return 0;
    }
    
    /* File format: first line contains the number of points, then each line is
    a record of path and correspondence info. We don't need the nuber of points
    because we read to EOF anyway. */
    
    read_res = fscanf(filein, "%d\n", &read_res);
    if (!read_res) return 0;
    
    do {
        /* Defaults: */
        path_buf->prev = -1;
        path_buf->next = -2;
        path_buf->prio = 4;
        path_buf->inlist = 0;
        path_buf->finaldecis = 1000000.0;

        for (alt_link = 0; alt_link < POSI; alt_link++) {
            path_buf->decis[alt_link] = 0.0;
            path_buf->linkdecis[alt_link] = -999;
        }
        
        /* Rest of values: */
        read_res = fscanf(filein, "%d %f %f %f %d %d %d %d\n",\
            &read_res, &(path_buf->x[0]), &(path_buf->x[1]), &(path_buf->x[2]),
            &(cor_buf->p[0]), &(cor_buf->p[1]), &(cor_buf->p[2]),
            &(cor_buf->p[3]) );
        
        if (read_res != 8) {
            targets = 0;
            break;
        }
        
        cor_buf->nr = ++targets;
        
        cor_buf++;
        path_buf++;
    } while (!feof(filein));
    
    fclose(filein);
    return targets;
}

