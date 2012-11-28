#include <stdlib.h>
#include "ptv.h"

/* declaration */
void qs_coord2d_x ();
void qs_target_y ();
void qs_coord2d_pnr ();
void qs_con ();

/* Write exterior and interior orientation, and - if available, parameters for
   distortion corrections.

   Arguments:
   Exterior Ex - exterior orientation.
   Interior I - interior orientation.
   Glass G - glass parameters.
   ap_52 addp - optional additional (distortion) parameters. NULL is fine if
      add_file is NULL.
   char *filename - path of file to contain interior, exterior and glass
      orientation data.
   char *add_file - path of file to contain added (distortions) parameters.
*/
void write_ori (Ex, I, G, ap, filename, add_file)
Exterior Ex;
Interior I;
Glass    G;
ap_52 ap;
char *filename, *add_file;
{
  FILE	*fp;
  int  	i;

  fp = fopen (filename, "w");
  fprintf (fp, "%11.4f %11.4f %11.4f\n    %10.7f  %10.7f  %10.7f\n\n",
	   Ex.x0, Ex.y0, Ex.z0, Ex.omega, Ex.phi, Ex.kappa);
  for (i=0; i<3; i++)  fprintf (fp, "    %10.7f %10.7f %10.7f\n",
				Ex.dm[i][0], Ex.dm[i][1], Ex.dm[i][2]);
  fprintf (fp,"\n    %8.4f %8.4f\n    %8.4f\n", I.xh, I.yh, I.cc);
  fprintf (fp,"\n    %20.15f %20.15f  %20.15f\n", G.vec_x, G.vec_y, G.vec_z);
  fclose (fp);
  
  if (add_file == NULL) return;
  fp = fopen (add_file, "w");
  fprintf (fp, "%f %f %f %f %f %f %f", ap.k1, ap.k2, ap.k3, ap.p1, ap.p2,
    ap.scx, ap.she);
  fclose (fp);
}


int read_ori (Ex, I, G, ori_file, addp, add_file, add_fallback)
/*  read exterior and interior orientation, and - if available, parameters for
    distortion corrections.
    
    Arguments:
    Exterior *Ex - output buffer for exterior orientation.
    Interior *I - output buffer for interior orientation.
    Glass *G - output buffer for glass parameters.
    char *ori_file - path of file contatining interior and exterior orientation
        data
    ap_52 addp - output buffer for additional (distortion) parameters.
    char *add_file - path of file contatining added (distortions) parameters.
    char *add_fallback - path to file for use if add_file can't be openned.
    
    Returns:
    true value on success, false on failure. Failure can happen if add_file
    can't be opened, or the fscanf results are wrong, but if the additional
    parameters' file or fallback can't be opened, they're just assigned default
    values.
*/
Exterior *Ex;
Interior *I;
Glass    *G;
ap_52    *addp;
char	 *ori_file, *add_file, *add_fallback;
{
  FILE	*fp;
  int  	i, scan_res;

  if (!(fp = fopen_r(ori_file))) return 0;
  
  /* Exterior */
  scan_res = fscanf (fp, "%lf %lf %lf %lf %lf %lf",
	  &(Ex->x0), &(Ex->y0), &(Ex->z0),
	  &(Ex->omega), &(Ex->phi), &(Ex->kappa));
  if (scan_res != 6) return 0;
  
  /* Exterior rotation matrix */
  for (i=0; i<3; i++) {
    scan_res = fscanf (fp, " %lf %lf %lf",
        &(Ex->dm[i][0]), &(Ex->dm[i][1]), &(Ex->dm[i][2]));
    if (scan_res != 3) return 0;
  }

  /* Interior */
  scan_res = fscanf (fp, "%lf %lf %lf", &(I->xh), &(I->yh), &(I->cc));
  if (scan_res != 3) return 0;
  
  /* Glass */
  scan_res = fscanf (fp, "%lf %lf %lf", &(G->vec_x), &(G->vec_y), &(G->vec_z));
  if (scan_res != 3) return 0;
  
  fclose(fp);
  
  /* Additional: */
  fp = fopen(add_file, "r");
  if ((fp == NULL) && add_fallback) fp = fopen (add_fallback, "r");
  
  if (fp) {
    scan_res = fscanf (fp, "%lf %lf %lf %lf %lf %lf %lf",
        &(addp->k1), &(addp->k2), &(addp->k3), &(addp->p1), &(addp->p2),
        &(addp->scx), &(addp->she));
    fclose (fp);
  } else {
    printf("no addpar fallback used\n"); // Waits for proper logging.
    addp->k1 = addp->k2 = addp->k3 = addp->p1 = addp->p2 = addp->she = 0.0;
    addp->scx=1.0;
  }
  
  return 1;
}


FILE *fopen_r (filename)
char *filename;
/*	tries to open a file;
	gives a message, if it cannot open it
	and waits until it has been created 	 */
{
  FILE	*fpr;
  int  	count;

  fpr = fopen (filename, "r");
  if ( ! fpr)
    {
      printf ("could not open %s, please create this file\n", filename);

      /* wait until file can be opened */
      while ( ! fpr)	fpr = fopen (filename, "r");

      /* wait until file really created */
      for (count=0; count<100000; count++);
    }

  return (fpr);
}

void compose_name_plus_nr (basename, str, nr, filename)
char	basename[256], str[256], filename[256];
int		nr;
{
  char	nr_ch[256];

//  if (nr < 10)		sprintf (nr_ch, "00%1d", nr);
//  else if (nr < 100)	sprintf (nr_ch, "0%2d",  nr);
  if (nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if (nr < 100)	sprintf (nr_ch, "%2d",  nr);
  else	sprintf (nr_ch, "%3d",  nr);

  sprintf (filename, "%s%s%s", basename, str, nr_ch);
}

void compose_name_plus_nr_str (basename, str, nr, filename)
char	basename[256], str[256], filename[256];
int		nr;
{
  char	nr_ch[256];

  if (nr < 10)		sprintf (nr_ch, "%1d", nr);
  else if (nr < 100)	sprintf (nr_ch, "%2d",  nr);
  else 	sprintf (nr_ch, "%3d",  nr);

  sprintf (filename, "%s%s%s", basename, nr_ch, str);
}

/* find nearest neighbours */

int kill_in_list ( nr, num, ms_x, ms_y)
//Tcl_Interp* interp;
//int interp; //denis
int    	nr, num;
int    	ms_x, ms_y;
{
  int 	i, imin = 9999, intx, inty;
  double  x, y, d, dmin = 9999;

  if (zoom_f[nr] > 1)
    {
      printf ("cannot delete point from zoomed image");
     // Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
     // Tcl_Eval(interp, ".text delete 3");
     // Tcl_Eval(interp, ".text insert 3 $tbuf");
      return (0);
    }

  for (i=0; i<num; i++)
    {
      x = (double) ms_x - pix[nr][i].x;
      y = (double) ms_y - pix[nr][i].y;
      d = sqrt (x*x + y*y);
      if (d < dmin)
	{
	  dmin = d; imin = i;
	}
    }
  if (dmin > 10)	return (-1);	       	/*  limit: 10 pixel  */
  intx = (int) pix[nr][imin].x;
  inty = (int) pix[nr][imin].y;

 // drawcross (interp, intx, inty, cr_sz+1, nr, "magenta");

  for (i=imin; i<num; i++)  pix[nr][i] = pix[nr][i+1];

  return (imin);
}



int nearest_neighbour_geo (crd, num, x, y, eps)
coord_2d  crd[];
int       num;
double 	  x, y, eps;
{
  register int	j;
  int	       	j0, dj, pnr = -999;
  double       	d, dmin=1e20, xmin, xmax, ymin, ymax;

  xmin = x - eps;  xmax = x + eps;  ymin = y - eps;  ymax = y + eps;

  /* binarized search for start point of candidate search */
  for (j0=num/2, dj=num/4; dj>1; dj/=2)
    {
      if (crd[j0].x < xmin)  j0 += dj;
      else  j0 -= dj;
    }
  j0 -= 12;  if (j0 < 0)  j0 = 0;	       	/* due to trunc */

  for (j=j0; j<num; j++)		       	/* candidate search */
    {
      if (crd[j].x > xmax)  break;	       	/* finish search */

      if (crd[j].y > ymin  &&  crd[j].y < ymax)
	{
	  d = sqrt ((x-crd[j].x)*(x-crd[j].x) + (y-crd[j].y)*(y-crd[j].y));
	  if (d < dmin)
	    {
	      dmin = d; pnr = j;
	    }
	}
    }
  return (pnr);
}


int nearest_neighbour_pix (pix, num, x, y, eps)
target 	pix[];
int    	num;
double 	x, y, eps;
{
  register int	j;
  int	       	pnr = -999;
  double       	d, dmin=1e20, xmin, xmax, ymin, ymax;

  xmin = x - eps;  xmax = x + eps;  ymin = y - eps;  ymax = y + eps;

  for (j=0; j<num; j++)		    			/* candidate search */
    {
      if (pix[j].y>ymin && pix[j].y<ymax && pix[j].x>xmin && pix[j].x<xmax)
	{
	  d = sqrt ((x-pix[j].x)*(x-pix[j].x) + (y-pix[j].y)*(y-pix[j].y));
	  if (d < dmin)
	    {
	      dmin = d; pnr = j;
	    }
	}
    }
  return (pnr);
}


/***********************************************************************/
/***********************************************************************/

/* sorting routines */

/***********************************************************************/

/* bubble sorts */
void bubble_conlist (item, count)
correspond	*item;
int    		count;
{
	int			i,j;
	correspond	temp;

	for (i=1; i<count; ++i)  for (j=count-1; j>=i; --j)
	{
		if (item[j-1].corr > item[j].corr)
		{
			temp = *(item+j-1);  *(item+j-1) = *(item+j);  *(item+j) = temp;
		}
	}
}

/***********************************************************************/
/***********************************************************************/


/* quicksort algorithms for several issues */

/***********************************************************************/

/* quicksort of 2d coordinates in x-order */

void quicksort_coord2d_x (crd, num)
coord_2d	*crd;
int			num;
{
	qs_coord2d_x (crd, 0, num-1);
}



void qs_coord2d_x (crd, left, right)
coord_2d	*crd;
int			left, right;
{
	register int	i, j;
	double			xm;
	coord_2d		temp;

	i = left;	j = right;	xm = crd[(left+right)/2].x;

	do
	{
		while (crd[i].x < xm  &&  i<right)	i++;
		while (xm < crd[j].x  &&  j>left)	j--;

		if (i <= j)
		{
			temp = crd[i];
			crd[i] = crd[j];
			crd[j] = temp;
			i++;	j--;
		}
	}
	while (i <= j);

	if (left < j)	qs_coord2d_x (crd, left, j);
	if (i < right)	qs_coord2d_x (crd, i, right);
}


/***********************************************************************/

/* quicksort of 2d coordinates in pnr-order */

void quicksort_coord2d_pnr (crd, num)
coord_2d	*crd;
int	       	num;
{
  qs_coord2d_pnr (crd, 0, num-1);
}



void qs_coord2d_pnr (crd, left, right)
coord_2d	*crd;
int    		left, right;
{
  register int	i, j;
  double       	pnrm;
  coord_2d     	temp;

  i = left;	j = right;	pnrm = crd[(left+right)/2].pnr;

  do
    {
      while (crd[i].pnr < pnrm  &&  i<right)	i++;
      while (pnrm < crd[j].pnr  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = crd[i];
	  crd[i] = crd[j];
	  crd[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_coord2d_pnr (crd, left, j);
  if (i < right)	qs_coord2d_pnr (crd, i, right);
}




/***********************************************************************/


/* quicksort of targets in y-order */

void quicksort_target_y (pix, num)
target 	*pix;
int    	num;
{
  qs_target_y (pix, 0, num-1);
}




void qs_target_y (pix, left, right)
target 	*pix;
int    	left, right;
{
  register int	i, j;
  double			ym;
  target			temp;

  i = left;	j = right;	ym = pix[(left+right)/2].y;

  do
    {
      while (pix[i].y < ym  &&  i<right)	i++;
      while (ym < pix[j].y  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = pix[i];
	  pix[i] = pix[j];
	  pix[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_target_y (pix, left, j);
  if (i < right)	qs_target_y (pix, i, right);
}



/***********************************************************************/


/* quicksort for list of correspondences in order of match quality */
/* 4 camera version */

void quicksort_con (con, num)
n_tupel	*con;
int    	num;
{
  qs_con (con, 0, num-1);
}


void qs_con (con, left, right)
n_tupel	*con;
int    	left, right;
{
  register int	i, j;
  double       	xm;
  n_tupel      	temp;

  i = left;	j = right;	xm = con[(left+right)/2].corr;

  do
    {
      while (con[i].corr > xm  &&  i<right)	i++;
      while (xm > con[j].corr  &&  j>left)	j--;

      if (i <= j)
	{
	  temp = con[i];
	  con[i] = con[j];
	  con[j] = temp;
	  i++;	j--;
	}
    }
  while (i <= j);

  if (left < j)	qs_con (con, left, j);
  if (i < right)	qs_con (con, i, right);
}


/***SORTING ALGORIHTMUS****/

void sort(int n, float a[], int b[])
{
  int flag = 0, i, itemp;
  float ftemp;

  do {
    flag = 0;
    for(i=0; i<(n-1); i++)
      if(a[i] > a[i+1]) {
	ftemp =  a[i];
	itemp =  b[i];
	a[i] = a[i+1];
	b[i] = b[i+1];
	a[i+1] = ftemp;
	b[i+1] = itemp;
        flag = 1;
      }
  }while(flag);
}

/***********************************************************************/



