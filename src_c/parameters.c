/* Implementation for parameters handling routines. */

#include "parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* read_sequence_par() reads sequence parameters from a test file with the
   following format: each line is a value, first 4 values are image names,
   5th is the first number in the sequence, 6th line is the last value in the
   sequence.
   
   Arguments:
   char *filename - path to the text file containing the parameters.
   
   Returns:
   Pointer to a newly-allocated sequence_par structure. Don't forget to free
   the new memory that is allocated for the image names. If reading failed for 
   any reason, returns NULL.
*/
sequence_par* read_sequence_par(char *filename) {
    char line[SEQ_FNAME_MAX_LEN], *line_check;
    FILE* par_file;
    int cam, read_ok;
    sequence_par *ret;
    
    par_file = fopen(filename, "r");
    if (par_file == NULL) {
        return NULL;
    }
    
    ret = (sequence_par *) malloc(sizeof(sequence_par));
    ret->img_base_name = (char *) calloc(4, sizeof(char *));
    
    /* Note the assumption of 4 cameras. Fixing this requires changing the
       file format. */
    for (cam = 0; cam < 4; cam++) {
        line_check = fgets(line, SEQ_FNAME_MAX_LEN, par_file);
        if (line_check == NULL) goto handle_error;
        
        ret->img_base_name[cam] = (char *) malloc(SEQ_FNAME_MAX_LEN);
        strncpy(ret->img_base_name[cam], line, SEQ_FNAME_MAX_LEN);
    }
    
    if( (read_ok = fscanf(par_file, "%d\n", &(ret->first))) == 0)
        goto handle_error;
    if( (read_ok = fscanf(par_file, "%d\n", &(ret->last))) == 0)    
        goto handle_error;

    close(par_file);
    return ret;
    
handle_error:
    printf("Error reading sequence parameters from %s\n", filename);
    free(ret);
    close(par_file);
    return NULL;
}

