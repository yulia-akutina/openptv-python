/*******************************************************************************
**
** Title: ptv
**
** Author: Heinrich Stueer
**
** Description: Main modul of track.
**Dies ist eine abgespeckte Version vom Malik und Papantoniou (allerdings mit
**ein paar Anderungen)
** Created: 12.02.1998
** Changes:
**
*******************************************************************************/
#include "ptv.h"
#include "tracking_frame_buf.h"

#define STR_MAX_LEN 255

/*  allocate_tracking_structs() Allocates memory for information needed for
    tracking. Each of the output buffer arrays has a first dimension of 4 sets,
    which is for 4 consecutive frames held in memory at the same time.
    
    Arguments:
    target* targets[4][4] - a (4,4)-shape array of pointers to a target
        structure array. by set, then by camera. Currently assumes 4 cameras as
        the callers do, but should be generalized.
    corres* correspond[4] - an array of correspondence structures for each
        set and initial target.
    P* path_info[4] - an array of per-set arrays of target path information.
    int cams - number of cameras to make room for.
    int max_targets - number of targets to make room for.
*/
void allocate_tracking_structs(\
    target* targets[4][4], corres* correspond[4], P* path_info[4], \
    int cams, int max_targets)
{
  int i, k;
  frame scratch_frame;

  for (i = 0; i < 4; i++) {
    frame_init(&scratch_frame, cams, max_targets);
    path_info[i] = scratch_frame.path_info;
    correspond[i] = scratch_frame.correspond;
    
    for (k = 0; k < cams; k++) {
      targets[i][k] = scratch_frame.targets[k];
    }
  }
  trackallocflag = 1;
}

int seq_track_proc_c()
{
  int step, i, k;
  
  allocate_tracking_structs(t4, c4, mega, n_img, M);
  readseqtrackcrit ();
  /*load again first data sets*/
  step = seq_first;
  read_ascii_data(step);
  rotate_dataset();
  read_ascii_data(step+1);
  rotate_dataset();
  read_ascii_data(step+2);

  for(step = (seq_first+2); step < seq_last; step++)
    {
      tracking();
  //tracking(clientData, interp, argc, argv);
      rotate_dataset();
      write_ascii_data(step-2);
      read_ascii_data(step+1);
    }

  /*write last data_sets*/

  //tracking(clientData, interp, argc, argv);
tracking();
  rotate_dataset();
  write_ascii_data(step-2);
  rotate_dataset();
  write_ascii_data(step-1);


  for (i=0; i<4; i++)
    { free (mega[i]);free (c4[i]);
    for (k=0; k<n_img; k++) free (t4[i][k]);
    }

  //return TCL_OK;
return 0;
}


void read_ascii_data(int filenumber)
{
  int	i;

  m[3] = read_path_frame(c4[3], mega[3], "res/rt_is", filenumber);
  
  /* read targets of each camera */
  for (i = 0; i < n_img; i++) {
      nt4[3][i] = read_targets(t4[3][i], seq_name[i], filenumber);
  }
}

/**********************************************************************/
void write_ascii_data(int filenumber)
{
  FILE	*FILEOUT;
  char	fileout[256];
  int	i, set, j;

  set = 0;

  write_path_frame(c4[set], mega[set], m[set], "res/rt_is", "res/ptv_is", 
    NULL, filenumber);
  
  /* write targets of each camera */
  for (i=0; i<n_img; i++) {
        write_targets(t4[set][i], nt4[set][i], seq_name[i], filenumber);
  }
}


/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
void write_added(int filenumber)
{
  FILE	*FILEOUT;
  char	fileout[256];
  int	i, set;

  set = 0;

  if (filenumber < 10)       sprintf (fileout, "res/added.%1d", filenumber);
  else if (filenumber< 100)  sprintf (fileout, "res/added.%2d",  filenumber);
  else       sprintf (fileout, "res/added.%3d",  filenumber);

  /*  printf ("write file: %s\n",fileout); */
  FILEOUT = fopen (fileout, "w");
  if (! FILEOUT) printf("Can't open ascii file for writing\n");

  fprintf(FILEOUT, "%d\n", m[set]);

  for(i=0; i<m[set]; i++)
    {
      /*read dataset row by row*/
      fprintf(FILEOUT, "%4d %4d %10.3f %10.3f %10.3f %d\n",
	      mega[set][i].prev, mega[set][i].next, mega[set][i].x[0],
	      mega[set][i].x[1], mega[set][i].x[2], mega[set][i].prio);
    }
  fclose(FILEOUT);

}
/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
/**********************************************************************/
void write_addedback(int filenumber)
{
  FILE	*FILEOUT;
  char	fileout[256];
  int	i, set;

  set = 0;

  if (filenumber < 10)       sprintf (fileout, "res/added.%1d", filenumber);
  else if (filenumber< 100)  sprintf (fileout, "res/added.%2d",  filenumber);
  else       sprintf (fileout, "res/added.%3d",  filenumber);

  /*  printf ("write file: %s\n",fileout); */
  FILEOUT = fopen (fileout, "w");
  if (! FILEOUT) printf("Can't open ascii file for writing\n");

  fprintf(FILEOUT, "%d\n", m[set]);

  for(i=0; i<m[set]; i++)
    {
      /*read dataset row by row, prev/next order changed because backwards*/
      fprintf(FILEOUT, "%4d %4d %10.3f %10.3f %10.3f %d\n",
	      mega[set][i].prev, mega[set][i].next, mega[set][i].x[0],
	      mega[set][i].x[1], mega[set][i].x[2], mega[set][i].prio);
    }
  fclose(FILEOUT);



}
/* ************************************************************* */

/* ************************************************************* */

void read_ascii_datanew(int filenumber)
{
  FILE	*FILEIN;
  char	filein[256];
  int	i, j;
  int   dumy;
  double fdumy;

  m[3] = read_path_frame(c4[3], mega[3], "res/rt_is", filenumber);
  /* read ptv_is-file for prev and next info */

  if (filenumber < 10)       sprintf (filein, "res/ptv_is.%1d", filenumber);
  else if (filenumber< 100)  sprintf (filein, "res/ptv_is.%2d",  filenumber);
  else       sprintf (filein, "res/ptv_is.%3d",  filenumber);

  FILEIN = fopen (filein, "r");
  if (! FILEIN) printf("Can't open ascii file for reading\n");

  fscanf(FILEIN, "%d\n", &dumy);

  for(i=0; i<=m[3]; i++)
    {
      /*read dataset row by row*/
      fscanf(FILEIN, "%4d %4d %lf %lf %lf\n", &mega[3][i].prev, &mega[3][i].next, &fdumy, &fdumy, &fdumy);
    }
  fclose(FILEIN);
  /* end of read ptv_is-file for prev and next info */


  /* read added-file for prio info */

  if (filenumber < 10)       sprintf (filein, "res/added.%1d", filenumber);
  else if (filenumber< 100)  sprintf (filein, "res/added.%2d",  filenumber);
  else       sprintf (filein, "res/added.%3d",  filenumber);

  FILEIN = fopen (filein, "r");
  if (! FILEIN) printf("Can't open ascii file for reading\n");

  fscanf(FILEIN, "%d\n", &dumy);

  for(i=0; i<=m[3]; i++)
    {
      /*read dataset row by row*/
      fscanf(FILEIN, "%4d %4d %lf %lf %lf %d\n", &dumy, &dumy, &fdumy, &fdumy, &fdumy, &mega[3][i].prio);
    }
  fclose(FILEIN);
  /* end of read added-file for prio info */

  /* read targets of each camera */
  for (i = 0; i < n_img; i++) {
      nt4[3][i] = read_targets(t4[3][i], seq_name[i], filenumber);
  }
}

/**********************************************************************/
void write_ascii_datanew(int filenumber)
{
  FILE	*FILEOUT;
  char	fileout[256];
  int	i, set, j;

  set = 0;

  write_path_frame(c4[set], mega[set], m[set], "res/rt_is", "res/ptv_is", 
    NULL, filenumber);

  /* write targets of each camera */
  for (i=0; i<n_img; i++) {
        write_targets(t4[set][i], nt4[set][i], seq_name[i], filenumber);
  }
}
/* ************************************************************* */





int tracking()
{
 level1();
 level2();
 level3();

// return TCL_OK;
}

void level1(void)
{
  int i, ii, j, k, l;
  float trad[3];
  int liste[M], inliste, n2liste[POSI], inn2liste, n3liste[POSI], inn3liste;
  float seekx[3], esti, acc;
  int finish, flag;

  /* Define some varibles. This is done every time tracking is called
     (my be not necessary but is not a big problem
     trad is radius of correlation neigbourhood, trad1
     is search radius for next timestep with link*/

    trad[0]= tpar.dvxmax;
    trad[1]= tpar.dvymax;
    trad[2]= tpar.dvzmax;

  /* BEGIN TRACKING*/

  /* First start with highest priority this is if particle has already
     link to previous
     timstep current timstep is t[1]*/
  /*first search for tracks with previous links*/
  /*links named -1 or -2 are no links*/

  for(inliste=0, i=0; i<m[1]; i++)
    if (mega[1][i].prev > -1) {
      liste[inliste]=i;
      inliste++;
    }
  /*calculate possible tracks for t2 and t3 and calculate decision criteria*/
  if(inliste > 0) {
    for(i=0; i < inliste; i++) {
      for(j=0; j < 3; j++)
	seekx[j] = 2*mega[1][liste[i]].x[j]-mega[0][mega[1][liste[i]].prev].x[j];

      /*find neighbours in next timestep t = 2*/
      inn2liste = 0;
      neighbours(seekx, trad, n2liste, &inn2liste, 2);
      /* if no neighour is found no link will be established*/

      /*calculate decision criteria*/
      if(inn2liste > 0) {
	for(k=0; k<inn2liste; k++) {
	  for(j=0; j < 3; j++)
	    seekx[j] = 2*mega[2][n2liste[k]].x[j]-mega[1][liste[i]].x[j];


          /*find neigbours in next timestep t = 3*/
          inn3liste = 0;
          neighbours(seekx, trad, n3liste, &inn3liste, 3);


	  if(inn3liste == 0) {
	    /*if no neighour in t3 is found, give decision criteria artifical
	      value (100000) accelaration can be considered
	      as unbelivible large*/
  	    mega[1][liste[i]].decis[k] = 1000000.0;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  }
	  else {
	    for(esti = 1000000.0, l=0; l<inn3liste; l++) {
	      /*calculate for estimates the decision value*/
	      for(acc=0.0, ii=0; ii<3; ii++)
		acc += sqr(mega[1][liste[i]].x[ii] - 2*mega[2][n2liste[k]].x[ii]
			   + mega[3][n3liste[l]].x[ii]);

	      acc=sqrt(acc);
	      if(esti > acc)
		esti = acc;
	    }/*for(l....)*/
	    mega[1][liste[i]].decis[k] = esti;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  } /*else(inn2liste >0*/
	}/*for (k...)*/
	mega[1][liste[i]].inlist = inn2liste;
	if(inn2liste > 1)
	  sort(inn2liste, mega[1][liste[i]].decis, mega[1][liste[i]].linkdecis);

      }/*if(inn1liste > 0)*/

    }/*for(i=0....)*/

    /*establish links by streaming completly through the data*/

    do {
      finish = 0;

      for(i=0; i<inliste; i++) {

	if(mega[1][liste[i]].next < 0) {
	  if(mega[1][liste[i]].inlist > 0) {
	    /*in the following is a sorted list of decis assumed*/
	    flag = 1;
	    j = 0;
	    do {

	      if(mega[2][mega[1][liste[i]].linkdecis[j]].prev < 0) {
		/*found possible link*/
		mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		mega[2][mega[1][liste[i]].linkdecis[j]].prio = 0;
		mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		flag = 0;
	      }/*if(p2 == -1)*/

	      /*test exiting link if would be better*/
	      else if(mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].finaldecis
		      > mega[1][liste[i]].decis[j])
		{
		  /*current link is better and reset other link*/
		  mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].next = -2;

		  mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prio = 0;
		  mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		  flag = 0;
		  finish = 1;
		}

	      else {
             j++; /* if first choice is not possible then try next */
	      }
	    }while((j < mega[1][liste[i]].inlist) && flag);

	  }/*if(mega[1]....)*/
	  else {
	    mega[1][liste[i]].next=-1; /*No link could be established*/
	  }/*else*/
	}/*if(mega[1]  .next < 0)*/
      }/*for(i=0....)*/

    }while(finish);

  }/*if(inlist >0)*/

  /*END OF FIRST TRAIL*/
}

/*second if no previous link but in neigbouhood exist previous links*/

void level2(void)
{
  int i, ii, j, k, l;
  float trad[3];
  int liste[M], inliste, n1liste[POSI], inn1liste;
  int n2liste[POSI], inn2liste, n3liste[POSI], inn3liste;
  float seekx[3], esti, acc, vel[3];
  int finish, flag, nvel;


  /* Define some varibles. This is done every time tracking is called
     (my be not necessary but is not a big problem*/


    trad[0]= tpar.dvxmax;
    trad[1]= tpar.dvymax;
    trad[2]= tpar.dvzmax;

  /* BEGIN TRACKING*/
  /* Secondly start with second priority this is if particle has already no link to
     previous timstep but in Particles in neigbourhoud have. Current timstep is t[1]*/

  /*first search for tracks with no previous links ancd no next link*/
  /*links named -1 or -2 are no links*/

  for(inliste=0, i=0; i<m[1]; i++)
    if (mega[1][i].next < 0 && mega[1][i].prev < 0) {
      /*check if neighbours wihtin correlation have link*/
      for(j=0; j < 3; j++)
	seekx[j] = mega[1][i].x[j];
      /* search points in neigbourhood within coorelation lenght*/
      inn1liste = 0;
      neighbours(seekx, trad, n1liste, &inn1liste, 1);
      /*check if neighbours have previous link*/
      /*n1liste must be greater than 1 because neigbours will find the point i itself*/
      if(inn1liste > 1) {
	for(vel[0]=0.0, vel[1]=0.0, vel[2]=0.0, nvel=0, j=0; j<inn1liste; j++) {
	  if(n1liste[j]!=i)
	    if(mega[1][n1liste[j]].prev > -1){
	      for(l=0; l<3; l++)
		vel[l] += mega[1][n1liste[j]].x[l]- mega[0][mega[1][n1liste[j]].prev].x[l];
	      nvel++;
	    }
	}

	if(nvel > 0) {
	  /*intermediate storage of center of position in next frame */
	  for(l=0; l<3; l++)
	    mega[1][i].decis[l] = vel[l]/(float)nvel;
	  liste[inliste]=i;
	  inliste++;
	}
      }
    }

  /*calculate possible tracks for t2 and t3 and calculate decision criteria*/
  if(inliste > 0) {
    for(i=0; i < inliste; i++) {
      for(j=0; j < 3; j++)
	seekx[j] = mega[1][liste[i]].x[j] + mega[1][liste[i]].decis[j];

      /*find neighbours in next timestep t = 2*/
      inn2liste = 0;
      neighbours(seekx, trad, n2liste, &inn2liste, 2);
      /* if no neighour is found no link will be established*/

      /*calculate decision criteria*/
      if(inn2liste > 0) {
	for(k=0; k<inn2liste; k++) {
	  for(j=0; j < 3; j++)
	    seekx[j] = 2*mega[2][n2liste[k]].x[j]-mega[1][liste[i]].x[j];

	  /*find neigbours in next timestep t = 3*/
	  inn3liste = 0;
	  neighbours(seekx, trad, n3liste, &inn3liste, 3);

	  if(inn3liste == 0) {
	    /*if no neighour in t3 is found, give decision criteria artifical value (100000)
	      accelaration can be considered as unbelivible large*/
	    mega[1][liste[i]].decis[k] = 1000000.0;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  }
	  else {
	    for(esti = 1000000.0, l=0; l<inn3liste; l++) {
	      /*calculate for estimates the decision value*/
	      for(acc=0.0, ii=0; ii<3; ii++)
		acc += sqr(mega[1][liste[i]].x[ii] - 2*mega[2][n2liste[k]].x[ii]
			   + mega[3][n3liste[l]].x[ii]);


	      acc=sqrt(acc);
	      if(esti > acc)
		esti = acc;
	    }/*for(l....)*/
	    mega[1][liste[i]].decis[k] = esti;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  } /*else(inn2liste >0*/
	}/*for (k...)*/
	mega[1][liste[i]].inlist = inn2liste;
	if(inn2liste > 1)
	  sort(inn2liste, mega[1][liste[i]].decis, mega[1][liste[i]].linkdecis);

      }/*if(inn1liste > 0)*/

    }/*for(i=0....)*/

    /*establish links by streaming completly through the data*/

    do {
      finish = 0;

      for(i=0; i<inliste; i++) {
	if(mega[1][liste[i]].next < 0) {

	  if(mega[1][liste[i]].inlist > 0) {
	    /*in the following is a sorted list of decis assumed*/
	    flag = 1;
	    j = 0;
	    do {

	      if(mega[2][mega[1][liste[i]].linkdecis[j]].prev < 0) {
		/*found possible link*/
		mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		mega[2][mega[1][liste[i]].linkdecis[j]].prio = 1;
		mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		flag = 0;
	      }/*if(p2 == -1)*/

	      /*test exiting link if would be better*/
	      else if((mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].finaldecis
		       > mega[1][liste[i]].decis[j])
		      && (mega[2][mega[1][liste[i]].linkdecis[j]].prio >= 1))
		{
		  /*current link is better and reset other link*/
		  mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].next = -2;

		  mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prio = 1;
		  mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		  flag = 0;
		  finish = 1;
		}

	      else {
		j++; /* if first choice is not possible then try next */
	      }

	    }while((j < mega[1][liste[i]].inlist) && flag);

	  }/*if(mega[1]....)*/
	  else {
	    mega[1][liste[i]].next=-1; /*No link could be established*/
	  }/*else*/
	}/*if(mega[1]  .next<0)*/
      }/*for(i=0....)*/

    }while(finish);

  }/*if(inlist >0)*/

  /*END OF second TRAIL*/
}


/*Third if no previous link nor in neigbouhood exist previous links*/

void level3(void)
{
  int i, ii, j, k, l;
  float trad[3];
  int liste[M], inliste, n2liste[POSI], inn2liste, n3liste[POSI], inn3liste;
  float seekx[3], esti, acc;
  int finish, flag;


  /* Define some varibles. This is done every time tracking is called
     (my be not necessary but is not a big problem*/

    trad[0]= tpar.dvxmax;
    trad[1]= tpar.dvymax;
    trad[2]= tpar.dvzmax;


  /* BEGIN TRACKING*/


  /* Thirdly start with third priority this is if particle has no link to previous
     timstep and in Particles in neigbourhoud have not. Current timstep is t[1]*/

  /*first search for tracks with no previous links and no next link*/
  /*links named -1 or -2 are no links*/

  for(inliste=0, i=0; i<m[1]; i++)
    if (mega[1][i].next < 0 && mega[1][i].prev < 0) {
      liste[inliste]=i;
      inliste++;
    }

  /*calculate possible tracks for t2 and t3 and calculate decision criteria*/
  if(inliste > 0) {
    for(i=0; i < inliste; i++) {
      for(j=0; j < 3; j++)
	seekx[j] = mega[1][liste[i]].x[j];

      /*find neighbours in next timestep t = 2*/
      inn2liste = 0;
      neighbours(seekx, trad, n2liste, &inn2liste, 2);
      /* if no neighour is found no link will be established*/

      /*calculate decision criteria*/
      if(inn2liste > 0) {
	for(k=0; k<inn2liste; k++) {
	  for(j=0; j < 3; j++)
	    seekx[j] = 2*mega[2][n2liste[k]].x[j]-mega[1][liste[i]].x[j];
	  inn3liste = 0;
	  neighbours(seekx, trad, n3liste, &inn3liste, 3);
	  if(inn3liste == 0) {
	    /*if no neighour in t3 is found, give decision criteria artifical value (100000)
	      accelaration can be considered as unbelivible large*/
	    mega[1][liste[i]].decis[k] = 1000000.0;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  }
	  else {
	    for(esti = 1000000.0, l=0; l<inn3liste; l++) {
	      /*calculate estimates the decision value*/
	      for(acc=0.0, ii=0; ii<3; ii++)
		acc += sqr(mega[1][liste[i]].x[ii] - 2*mega[2][n2liste[k]].x[ii]
			   + mega[3][n3liste[l]].x[ii]);


	      acc=sqrt(acc);
	      if(esti > acc)
		esti = acc;
            }/*for(l....)*/
	    mega[1][liste[i]].decis[k] = esti;
	    mega[1][liste[i]].linkdecis[k] = n2liste[k];
	  } /*else(inn2liste >0*/
	}/*for (k...)*/
	mega[1][liste[i]].inlist = inn2liste;
	if(inn2liste > 1)
	  sort(inn2liste, mega[1][liste[i]].decis, mega[1][liste[i]].linkdecis);

      }/*if(inn1liste > 0)*/

    }/*for(i=0....)*/

    /*establish links by streaming completly through the data*/

    do {
      finish = 0;
      for(i=0; i<inliste; i++) {

	if(mega[1][liste[i]].next < 0) {
	  if(mega[1][liste[i]].inlist > 0) {
	    /*in the following is a sorted list of decis assumed*/
	    flag = 1;
	    j = 0;
	    do {

	      if(mega[2][mega[1][liste[i]].linkdecis[j]].prev < 0) {
		/*found possible link*/
		mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		mega[2][mega[1][liste[i]].linkdecis[j]].prio = 2;
		mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		flag = 0;
	      }/*if(p2 == -1)*/

	      /*test exiting link if would be better*/
	      else if((mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].finaldecis
		       > mega[1][liste[i]].decis[j])
		      && (mega[2][mega[1][liste[i]].linkdecis[j]].prio >= 2))
		{
		  /*current link is better and reset other link*/
		  mega[1][mega[2][mega[1][liste[i]].linkdecis[j]].prev].next = -2;

		  mega[1][liste[i]].next = mega[1][liste[i]].linkdecis[j];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prev = liste[i];
		  mega[2][mega[1][liste[i]].linkdecis[j]].prio = 2;
		  mega[1][liste[i]].finaldecis = mega[1][liste[i]].decis[j];
		  flag = 0;
		  finish = 1;
		}

	      else {
		j++; /* if first choice is not possible then try next */
	      }

	    }while((j < mega[1][liste[i]].inlist) && flag);

	  }/*if(mega[1]....)*/
	  else {
	    mega[1][liste[i]].next=-1; /*No link could be established*/
	  }/*else*/
	}/*if(mega[1]  < 0)*/
      }/*for(i=0....)*/
    }while(finish);
  }/*if(inlist >0)*/
  /*END OF THIRD TRAIL*/
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


void rotate_dataset(void)
{
  void * tmp;
  void * tmp2;
  int i;

  /*rotate dataset by changeing pointer*/
  tmp = mega[0];
  mega[0] = mega[1];
  mega[1] = mega[2];
  mega[2] = mega[3];
  mega[3] =tmp;

  /*rotate counter*/
  m[0] = m[1];
  m[1] = m[2];
  m[2] = m[3];

  tmp = c4[0];
  c4[0] = c4[1];
  c4[1] = c4[2];
  c4[2] = c4[3];
  c4[3] =tmp;

  for(i=0; i<4; i++) {
  tmp2 = t4[0][i];
  t4[0][i] = t4[1][i];
  t4[1][i] = t4[2][i];
  t4[2][i] = t4[3][i];
  t4[3][i] =tmp2;

  nt4[0][i] = nt4[1][i];
  nt4[1][i] = nt4[2][i];
  nt4[2][i] = nt4[3][i];
  }
}



void neighbours(float seekx[], float radi[], int nliste[], int *innliste, int set)
{
  int i;
  /*search for points in srearch radius. No sorted list is supported,
    although sorted in z would increase speed*/

  for(i=0; i< m[set]; i++) {
    if(fabs(seekx[0] - mega[set][i].x[0]) < radi[0])
      if(fabs(seekx[1] - mega[set][i].x[1]) < radi[1])
        if(fabs(seekx[2] - mega[set][i].x[2]) < radi[2]){
          nliste[*innliste]=i;
          (*innliste)++;
          if(*innliste > POSI)
	    printf("More Points found than can be supported! Reduce search area or increase POSI\n");
        }
  }
}
