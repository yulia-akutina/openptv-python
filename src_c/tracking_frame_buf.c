/*
Implementation the tracking frame buffers. Each frame holds target information
for all cameras, correspondence information and path links information.

Contains functions for dealing with the frame content, such as comparing and
reading targets, correspondences, etc.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
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

/* Writes targets to a file. The number of targets is written to the first
 * line, then each line is one target.
 * 
 * Arguments:
 * target buffer[] - an allocated array of target structs to fill in from
 *   files.
 * int num_targets - length of the targets buffer.
 * char* file_base - base name of the files to read, to which a frame number
 *   and the suffix '_targets' is added.
 * int frame_num - number of frame to add to file_base. A value of 0 or less
 *   means that no frame number should be added.
 * 
 * Returns:
 * True value on success, or 0 if an error occurred.
*/
int write_targets(target buffer[], int num_targets, char* file_base, \
    int frame_num) {
    
    FILE *FILEOUT;
    int	tix, printf_ok;
    char fileout[STR_MAX_LEN + 1];
    
    sprintf(fileout, "%s%04d%s", file_base, frame_num, "_targets");

    FILEOUT = fopen(fileout, "w");
    if (! FILEOUT) {
        printf("Can't open ascii file: %s\n", fileout);
        return 0;
    }
    
    if (fprintf (FILEOUT, "%d\n", num_targets) <= 0) {
        printf("Write error in file %s\n", fileout);
        return 0;
    }
    
    for (tix = 0; tix < num_targets; tix++)	{
	    printf_ok = fprintf(FILEOUT, "%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n",
		    buffer[tix].pnr, buffer[tix].x,
		    buffer[tix].y, buffer[tix].n,
		    buffer[tix].nx, buffer[tix].ny,
		    buffer[tix].sumg, buffer[tix].tnr);
        if (printf_ok <= 0) {
            printf("Write error in file %s\n", fileout);
            return 0;
        }
	}
    fclose (FILEOUT);
    
    return 1;
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

/* write_path_frame() writes the correspondence and linkage information for a
 * frame with the next and previous frames. The information is weirdly 
 * distributed among two files. The rt_is file holds correspondence and 
 * position data, and the ptv_is holds linkage data.
 *
 * Arguments:
 * corres *cor_buf - an array of corres structs to write to the files.
 * P *path_buf - same for path info structures.
 * int num_parts - number of particles represented by the cor_buf and path_buf
 *   arrays, i.e. the arrays' length.
 * char* corres_file_base, *linkage_file_base - base names of the output
 *   correspondence and likage files respectively, to which a frame number
 *   is added. Without separator.
 * int frame_num - number of frame to add to file_base. The '.' separator is 
 * added between the name and the frame number.
 * 
 * Returns:
 * True on success. 0 on failure.
 */
int write_path_frame(corres *cor_buf, P *path_buf, int num_parts,\
    char *corres_file_base, char *linkage_file_base, int frame_num) {

    FILE *corres_file, *linkage_file;
    char corres_fname[STR_MAX_LEN + 1], linkage_fname[STR_MAX_LEN + 1];
    int	pix, j;

    sprintf(corres_fname, "%s.%d", corres_file_base, frame_num);
    corres_file = fopen (corres_fname, "w");
    if (corres_file == NULL) {
        printf("Can't open file %s for writing\n", corres_fname);
        return 0;
    }
    
    sprintf(linkage_fname, "%s.%d", linkage_file_base, frame_num);
    linkage_file = fopen (linkage_fname, "w");
    if (linkage_file == NULL) {
        printf("Can't open file %s for writing\n", linkage_fname);
        return 0;
    }

    fprintf(corres_file, "%d\n", num_parts);
    fprintf(linkage_file, "%d\n", num_parts);

    for(pix = 0; pix < num_parts; pix++) {
        fprintf(linkage_file, "%4d %4d %10.3f %10.3f %10.3f\n",
	        path_buf[pix].prev, path_buf[pix].next, path_buf[pix].x[0],
	        path_buf[pix].x[1], path_buf[pix].x[2]);
        
        fprintf(corres_file, "%4d %9.3f %9.3f %9.3f %4d %4d %4d %4d\n",
	        pix + 1, path_buf[pix].x[0], path_buf[pix].x[1], path_buf[pix].x[2],
    	    cor_buf[pix].p[0], cor_buf[pix].p[1], cor_buf[pix].p[2],
            cor_buf[pix].p[3]);
    }

    fclose(corres_file);
    fclose(linkage_file);
    return 1;
}

/* init_frame() initializes a frame object, allocates its arrays and sets up 
 * the frame data.
 *  
 * Arguments:
 * int num_cams - number of cameras per frame.
 * int max_targets - number of elements to allocate for the different buffers
 *     held by a frame.
 */
void frame_init(frame *new_frame, int num_cams, int max_targets) {
    int cam;
    
    new_frame->path_info = (P *) calloc(max_targets, sizeof(P));
    new_frame->correspond = (corres *) calloc(max_targets, sizeof(corres));
    
    new_frame->targets = (target**) calloc(num_cams, sizeof(target*));
    new_frame->num_targets = (int *) calloc(max_targets, sizeof(int));
    
    for (cam = 0; cam < num_cams; cam++) {
        new_frame->targets[cam] = \
            (target *) calloc(max_targets, sizeof(target));
        new_frame->num_targets[cam] = 0;
    }
    
    new_frame->num_cams = num_cams;
    new_frame->max_targets = max_targets;
    new_frame->num_parts = 0;
    return new_frame;
}

/* free_frame() frees all memory allocated for the frame arrays.
 * 
 * Arguments:
 * frame *self - the frame to free.
 */
void free_frame(frame *self) {
    free(self->path_info);
    self->path_info = NULL;
    
    free(self->correspond);
    self->correspond = NULL;
    
    free(self->num_targets);
    self->num_targets = NULL;
    
    for (; self->num_cams > 0; self->num_cams--) {
        free(self->targets[self->num_cams - 1]);
        self->targets[self->num_cams - 1] = NULL;
    }
    
    free(self->targets);
    self->targets = NULL;
}

/* read_frame() reads all of the frame associated data: correspondnces,
 * targets, and whatever else is needed.
 * 
 * Arguments:
 * frame *self - the frame object to fill with the data read.
 * char *file_base - base name of the correspondence file to read, to which a 
 *   frame number is added. Without separator.
 * char **target_file_base - an array of strings following the same rules as
 *   for file_base; one for each camera in the frame.
 * int frame_num - number of frame to add to file_base. A value of 0 or less
 *   means that no frame number should be added. The '.' separator is added
 *   between the name and the frame number.
 *
 * Returns:
 * True on success, false otherwise. In case of failure, the state of frame is
 * undefined.
 */
int read_frame(frame *self, char *file_base, char **target_file_base,
    int frame_num)
{
    int cam;
    
    self->num_parts = read_path_frame(self->correspond, self->path_info,
        file_base, frame_num);
    if (self->num_targets == 0) return 0;
    
    for (cam = 0; cam < self->num_cams; cam++) {
        self->num_targets[cam] = read_targets(
            self->targets[cam], target_file_base[cam], frame_num);
        if (self->num_targets[cam] == 0) return 0;
    }
    
    return 1;
}

/* write_frame() writes all of the frame associated data: correspondnces,
 * targets, path info, and whatever else is needed.
 * 
 * Arguments:
 * frame *self - the frame whose data is to be written.
 * char *corres_file_base - base name of the correspondence file to read, to 
 *   which a frame number is added. Without separator.
 * char *linkage_file_base - same as corres_file_base, for the linkage file.
 * char **target_file_base - an array of strings following the same rules as
 *   for file_base; one for each camera in the frame.
 * int frame_num - number of frame to add to file_base. A value of 0 or less
 *   means that no frame number should be added. The '.' separator is added
 *   between the name and the frame number.
 *
 * Returns:
 * True on success, false otherwise. In case of failure, the state of output 
 * files is undefined.
 */
int write_frame(frame *self, char *corres_file_base, char *linkage_file_base,
    char **target_file_base, int frame_num)
{
    int cam, status;
    
    status = write_path_frame(self->correspond, self->path_info,
        self->num_parts, corres_file_base, linkage_file_base, frame_num);
    if (status == 0) return 0;
    
    for (cam = 0; cam < self->num_cams; cam++) {
        status = write_targets(
            self->targets[cam], self->num_targets[cam], target_file_base[cam],
            frame_num);
        if (status == 0) return 0;
    }
    
    return 1;
}

/* fb_init() allocates a frame buffer object and creates the frames
 * it buffers. The buffer is a ring-buffer based on a double-size vector.
 *  
 * Arguments:
 * int buf_len - number of frames in the buffer.
 * int num_cams - number of cameras per frame.
 * int max_targets - number of elements to allocate for the different buffers
 *     held by a frame.
 * char *rt_file_base
 * 
 * Returns:
 * a pointer to the new dynamically allocated frame-buffer structure.
 */

void fb_init(framebuf *new_buf, int buf_len, int num_cams, int max_targets,\
    char *corres_file_base, char *linkage_file_base, char **target_file_base)
{
    frame *alloc_frame;
    
    new_buf->buf_len = buf_len;
    new_buf->num_cams = num_cams;
    new_buf->corres_file_base = corres_file_base;
    new_buf->linkage_file_base = linkage_file_base;
    new_buf->target_file_base = target_file_base;
    
    new_buf->_ring_vec = (frame *) calloc(buf_len*2, sizeof(frame *));
    new_buf->buf = new_buf->_ring_vec + buf_len;
    
    while (new_buf->buf != new_buf->_ring_vec) {
        new_buf->buf--;
        
        alloc_frame = (frame *) malloc(sizeof(frame));
        frame_init(alloc_frame, num_cams, max_targets);
        
        /* The second half of _ring_vec points to the same objects as the
           first half. This allows direct indexing like fb->buf[buf_len - 1]
           even when fb->buf moves forward. */
        *(new_buf->buf) = alloc_frame;
        *(new_buf->buf + buf_len) = alloc_frame;
    }
    /* Leaves new_buf->buf pointing to the beginning of _ring_vec,
    which is what we want. */
}

/* fb_free() frees all memory allocated for the frames and ring vector in a
 * framebuf object.
 * 
 * Arguments:
 * framebuf *self - the framebuf holding the memory to free.
 */
void fb_free(framebuf *self) {
    self->buf = self->_ring_vec;
    
    while (self->buf != self->_ring_vec + self->buf_len) {
        free_frame(*(self->buf));
        free(*(self->buf));
    }
    self->buf = NULL;
    
    free(self->_ring_vec);
    self->_ring_vec = NULL;
}

