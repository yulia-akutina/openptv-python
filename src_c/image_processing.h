/* Forward declarations for various image processing routines collected from
image_processing.c, segmentation.c and peakfitting.c */

#ifndef IMAGE_PROCESSING_C
#define IMAGE_PROCESSING_C

/* I don't know why this has to be "unsigned" but I have no time to test
   dropping it */
void copy_images (unsigned char *img1, unsigned char *img2);
int peak_fit_new (unsigned char	*img, char par_file[],
    int xmin, int xmax, int ymin, int ymax, target pix[], int nr);

void simple_connectivity(unsigned char *img0,
    unsigned char *img, char par_file[],
    int xmin, int xmax, int ymin, int ymax,
    target pix[], int nr, int *num);

void targ_rec (unsigned char *img0, unsigned char *img, char par_file[], 
    int xmin, int xmax, int ymin, int ymax,
    target pix[], int nr, int *num);

void highpass(unsigned char *img, unsigned char *img_hp, int dim_lp,
    int filter_hp, int field);
void lowpass_3(unsigned char *img, unsigned char *img_lp);
void filter_3(unsigned char *img, unsigned char *img_lp);
void unsharp_mask(int n, unsigned char *img0, unsigned char *img_lp);
void split(unsigned char *img, int field);

#endif

