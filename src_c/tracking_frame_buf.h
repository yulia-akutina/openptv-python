/*
Definition for the tracking frame buffers. Each frame holds target information
for all cameras, correspondence information and path links information.
*/

#ifndef TRACKING_FRAME_BUF_H
#define TRACKING_FRAME_BUF_H

typedef struct
{
  int     pnr;
  double  x, y;
  int     n, nx, ny, sumg;
  int     tnr;
}
target;

typedef struct
{
  int nr;
  int p[4];
}
corres;

typedef struct Pstruct
{
  float x[3]; /*coordinates*/
  int prev, next; /*pointer to prev or next link*/
  int prio; /*Prority of link is used for differen levels*/
  float decis[POSI]; /*Bin for decision critera of possible links to next dataset*/
  float finaldecis; /*final decision critera by which the link was established*/
  int linkdecis[POSI]; /* pointer of possible links to next data set*/
  int inlist; /* Counter of number of possible links to next data set*/
} P;

typedef struct {
    P *path_info;
    corres *correspond;
    target **targets;
    int num_cams, max_targets;
} frame;

frame* create_frame(int num_cams, int max_targets);
void free_frame(frame *self);


typedef struct {
    frame *buf;
    int buf_len, num_cams;
    char *rt_file_base;
    char **target_file_base;
} framebuf;

framebuf* frame_buffer_create(int buff_len, int num_cams, int max_targets,\
    char *rt_file_base, char **target_file_base);
void frame_buffer_free(framebuf *self);
void frame_buffer_shift(framebuf *self);

#endif
