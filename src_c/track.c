/*******************************************************************

Routine:	      	track.c

Author/Copyright:     	Jochen Willneff

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
	               	CH - 8093 Zurich

Creation Date:		Beginning: February '98
                        End: far away

Description:   	        Tracking of particles in image- and objectspace

Routines contained:    	trackcorr_c

*******************************************************************/
#include "ptv.h"
#include "tracking_frame_buf.h"

void write_addedback();

/* For now, we'll use a file-global framebuf. It is the lesser evil compared to
the current state of using a scatter of projectwide-globals, and will change to
a framebuffer passed around when we get to the point we feel safe in changing 
the functions' interface.
*/
framebuf *fb;

/* The buffer space required for this algorithm: 
Note that MAX_TARGETS is taken from the global M, but I want a separate
definition because the fb created here should be local, not used outside
this file. 
*/
#define BUFSPACE 4
#define MAX_TARGETS 20000

int trackcorr_c_init () {
    int step, img;
    double Ymin=0, Ymax=0,lmax;
    char** target_file_base;
    
    /* Remaining globals:
    fb - from this file, for this file only.
    n_img, seq_name - set in jw_ptv.c in init_proc_c().
    seq_first, seq_lasti, tpar.*, all parameters of volumedimension - set in 
        readseqtrackcrit().
    see below for communication globals.
    */

    /* read configuration */
    readseqtrackcrit ();
    
    target_file_base = (char **) calloc(n_img, sizeof(char *));
    for (img = 0; img < n_img; img++) {
        target_file_base[img] = seq_name[img];
    }
    fb = (framebuf *) malloc(sizeof(framebuf));
    fb_init(fb, 4, n_img, MAX_TARGETS, "res/rt_is", "res/ptv_is", "res/added",
        target_file_base);

    /* Prime the buffer with first frames */
    for (step = seq_first; step < seq_first + 3; step++) {
        fb_read_frame_at_end(fb, step);
        fb_next(fb);
    }
    fb_prev(fb);

    lmax=sqrt((tpar.dvxmin-tpar.dvxmax)*(tpar.dvxmin-tpar.dvxmax)
	    +(tpar.dvymin-tpar.dvymax)*(tpar.dvymin-tpar.dvymax)
	    +(tpar.dvzmin-tpar.dvzmax)*(tpar.dvzmin-tpar.dvzmax));

    volumedimension (&X_lay[1], &X_lay[0], &Ymax, &Ymin, &Zmax_lay[1], &Zmin_lay[0]);

    // Denis - globals below are passed to trackcorr_c_loop
    lmax_track=lmax;
    ymin_track=Ymin;
    ymax_track=Ymax;


    // Denis - globals below are used in trackcorr_finish
    npart=0;
    nlinks=0;
}

/* reset_foundpix_array() sets default values for foundpix objects in an array.
 *
 * Arguments:
 * foundpix *arr - the array to reset
 * int arr_len - array length
 * int num_cams - number of places in the whichcam member of foundpix.
 */
void reset_foundpix_array(foundpix *arr, int arr_len, int num_cams) {
    int i, cam;
    for (i = 0; i < arr_len; i++) {
	    arr[i].ftnr = -1;
        arr[i].freq = 0;
        
        for(cam = 0; cam < num_cams; cam++) {
            arr[i].whichcam[cam] = 0;
	    }
    }
}

/* copy_foundpix_array() copies foundpix objects from one array to another.
 *
 * Arguments:
 * foundpix *dest, *src - src is the array to copy, dest receives it.
 * int arr_len - array length
 * int num_cams - number of places in the whichcam member of foundpix.
 */
void copy_foundpix_array(foundpix *dest, foundpix *src, int arr_len, 
    int num_cams) 
{
    int i, cam;
    for (i = 0; i < arr_len; i++) {
        dest[i].ftnr = src[i].ftnr;
        dest[i].freq = src[i].freq;
        for (cam = 0; cam < num_cams; cam++) {
            dest[i].whichcam[cam] = src[i].whichcam[cam];
        }
    }
}

int trackcorr_c_loop (int step, double lmax, double Ymin, double Ymax,
    int display) {
   /* sequence loop */
    char  val[256], buf[256];
    int i, j, h, k, mm, kk,  okay=0, invol=0;
    int zaehler1, zaehler2,philf[4][4];
    int count1=0, count2=0, count3=0, lost =0, zusatz=0;
    int intx0, intx1, inty0, inty1;
    int intx2, inty2,intx3, inty3;
    int quali=0;
    double x1[4], y1[4], x2[4], y2[4], angle, acc, angle0, acc0,  dl;
    double xr[4], xl[4], yd[4], yu[4], angle1, acc1;
    double X1, Y1, Z1, X0, Y0, Z0, X2, Y2, Z2;
    double X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6;
    double xp[4], yp[4], xc[4], yc[4], xn[4], yn[4];
    double rr;
    int flag_m_tr=0;
    
    /* Shortcuts to inside current frame */
    P *curr_path_inf, *ref_path_inf;
    corres *curr_corres, *ref_corres;
    target **curr_targets, **ref_targets;
    int _ix; /* For use in any of the complex index expressions below */
    int _frame_parts; /* number of particles in a frame */

    /* Remaining globals:
    all those in trackcorr_c_init.
    calibration globals.
    */
    
    foundpix *w, *wn, p16[16];
    sprintf (buf, "Time step: %d, seqnr: %d, Particle info:", step- seq_first, step);
    count1=0; lost =0; zusatz=0;
    
    curr_targets = fb->buf[1]->targets;
    
    /* try to track correspondences from previous 0 - corp, variable h */
    for (h = 0; h < fb->buf[1]->num_parts; h++) {
	    X1=Y1=Z1=X0=Y0=Z0=X2=Y2=Z2=X5=Y5=Z5=X3=Y3=Z3=X4=Y4=Z4=X6=Y6=Z6=-999;
        
        curr_path_inf = &(fb->buf[1]->path_info[h]);
        curr_corres = &(fb->buf[1]->correspond[h]);
        
	    curr_path_inf->inlist = 0;
        reset_foundpix_array(p16, 16, fb->num_cams);
        
	    /* 3D-position */
	    X1 = curr_path_inf->x[0];
	    Y1 = curr_path_inf->x[1];
	    Z1 = curr_path_inf->x[2];

	    /* use information from previous to locate new search position
	       and to calculate values for search area */
	    if (curr_path_inf->prev >= 0) {
            ref_path_inf = &(fb->buf[0]->path_info[curr_path_inf->prev]);
	        X0 = ref_path_inf->x[0];
	        Y0 = ref_path_inf->x[1];
    	    Z0 = ref_path_inf->x[2];
	        X2 = 2*X1 - X0;
	        Y2 = 2*Y1 - Y0;
    	    Z2 = 2*Z1 - Z0;

	        for (j = 0; j < fb->num_cams; j++) {
        		img_coord (X2, Y2, Z2, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], &yn[j]);
		        metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &x1[j], &y1[j], chfield);
	        }
	    } else {  
            X2=X1; Y2=Y1; Z2=Z1;
	        for (j=0; j < fb->num_cams; j++) {
	            if (curr_corres->p[j] == -1) {
	                img_coord (X2, Y2, Z2, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], 
                        &yn[j]);
	                metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &x1[j], 
                        &y1[j], chfield);
	            } else {
                    _ix = curr_corres->p[j];
                    x1[j] = curr_targets[j][_ix].x;
                    y1[j] = curr_targets[j][_ix].y;
                }
            }
	    } 
        
	    /* calculate searchquader and reprojection in image space */
	    searchquader(X2, Y2, Z2, &xr, &xl, &yd, &yu);

	    /* search in pix for candidates in next time step */
	    for (j=0;j < fb->num_cams; j++) {
	        zaehler1 = candsearch_in_pix (fb->buf[2]->targets[j], 
                fb->buf[2]->num_targets[j], x1[j], y1[j],
			    xl[j], xr[j], yu[j], yd[j], &philf[j]);
            
	        for(k = 0; k < 4; k++) {
                _ix = philf[j][k];
			    if(_ix == -999) {
				    p16[j*4+k].ftnr=-1;
			    }else{
                    p16[j*4+k].whichcam[j]=1;
				    p16[j*4+k].ftnr=fb->buf[2]->targets[j][_ix].tnr;
			    }
		    }
	    }
	    /* end of search in pix */
        
	    /* fill and sort candidate struct */
	    sortwhatfound(&p16, &zaehler1);
	    w = (foundpix *) calloc (zaehler1, sizeof (foundpix));
        
	    if (zaehler1 > 0) count2++;
        copy_foundpix_array(w, p16, zaehler1, fb->num_cams);
	    /*end of candidate struct */

	    /* check for what was found */
	    for (mm=0; mm<zaehler1;mm++) { /* zaehler1-loop */
	        /* search for found corr of current the corr in next
		    with predicted location */

            reset_foundpix_array(p16, 16, fb->num_cams);

	        /* found 3D-position */
            ref_path_inf = &(fb->buf[2]->path_info[w[mm].ftnr]);
	        X3 = ref_path_inf->x[0];
	        Y3 = ref_path_inf->x[1];
	        Z3 = ref_path_inf->x[2];

	        if (curr_path_inf->prev >= 0) {
		        X5=0.5*(5.0*X3-4.0*X1+X0);
		        Y5=0.5*(5.0*Y3-4.0*Y1+Y0);
		        Z5=0.5*(5.0*Z3-4.0*Z1+Z0);
	        } else {
		        X5=2*X3-X1;
		        Y5=2*Y3-Y1;
		        Z5=2*Z3-Z1; 
            }
            searchquader(X5, Y5, Z5, &xr, &xl, &yd, &yu);

	        for (j = 0; j < fb->num_cams; j++) {
		        img_coord (X5, Y5, Z5, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], &yn[j]);
		        metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &x2[j], &y2[j], chfield);
	        }

	        /* search for candidates in next time step */
	        for (j=0; j < fb->num_cams; j++) {
	            zaehler2 = candsearch_in_pix (fb->buf[3]->targets[j], 
                    fb->buf[3]->num_targets[j], x1[j], y1[j],
					xl[j], xr[j], yu[j], yd[j], &philf[j]);

		        for(k = 0; k < 4; k++) {
				    if (philf[j][k] == -999) {
                        p16[j*4+k].ftnr=-1;
				    } else {
				        if (fb->buf[3]->targets[j][philf[j][k]].tnr != -1) {
                            _ix = philf[j][k];
                            p16[j*4+k].ftnr = fb->buf[3]->targets[j][_ix].tnr;
                            p16[j*4+k].whichcam[j] = 1;
					    }
				    }
		        }
		    }
	        /* end of search in pix */

	        /* fill and sort candidate struct */
	        sortwhatfound(&p16, &zaehler2);
	        wn = (foundpix *) calloc (zaehler2, sizeof (foundpix));
	        if (zaehler2 > 0) count3++;
            copy_foundpix_array(wn, p16, zaehler2, fb->num_cams);

	        /*end of candidate struct */
	        /* ************************************************ */
	        for (kk=0; kk < zaehler2; kk++)  { /* zaehler2-loop */
                ref_path_inf = &(fb->buf[3]->path_info[wn[kk].ftnr]);
                X4 = ref_path_inf->x[0];
        		Y4 = ref_path_inf->x[1];
		        Z4 = ref_path_inf->x[2];

		        okay=0; rr=1000000; quali=0; dl=0;
        		acc=2*tpar.dacc; angle=2*tpar.dangle;
		        acc0=2*tpar.dacc; angle0=2*tpar.dangle;
        		acc1=2*tpar.dacc; angle1=2*tpar.dangle;

		        /* displacement check */
		        if ( tpar.dvxmin < (X4-X3) && (X4-X3) < tpar.dvxmax &&
		            tpar.dvymin < (Y4-Y3) && (Y4-Y3) < tpar.dvymax &&
		            tpar.dvzmin < (Z4-Z3) && (Z4-Z3) < tpar.dvzmax ) 
                { 
                    okay=1;

		            if ( okay ==1 ) {
		                dl=(sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3))
			            +sqrt((X4-X3)*(X4-X3)+(Y4-Y3)*(Y4-Y3)+(Z4-Z3)*(Z4-Z3)))/2;

		                angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle1, &acc1);

		                if (curr_path_inf->prev >= 0) {
                            angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, 
                                &angle0, &acc0);
		                } else {
                            acc0=acc1; angle0=angle1;
                        }

		                acc=(acc0+acc1)/2; angle=(angle0+angle1)/2;
		                quali=wn[kk].freq+w[mm].freq;
		                rr=1000000;

		                if ((acc<tpar.dacc && angle<tpar.dangle) || (acc<tpar.dacc/10)) {
			                rr = (dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
			                curr_path_inf->decis[curr_path_inf->inlist] = rr;
			                curr_path_inf->linkdecis[curr_path_inf->inlist] = w[mm].ftnr;
			                curr_path_inf->inlist++;
			            }
		                okay=0;
		            }
		        }
	        }   /* end of zaehler2-loop */
	        okay=0;

	        /* creating new particle position */
	        /* *************************************************************** */
	        for (j = 0;j < fb->num_cams; j++) {
		        img_coord (X5, Y5, Z5, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], &yn[j]);
		        metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
	        }

	        /* reset img coord because of n_img smaller 4 */
	        for (j=0;j < fb->num_cams; j++) { x2[j]=-1e10; y2[j]=-1e10; }

	        /* search for unused candidates in next time step */
	        for (j = 0;j < fb->num_cams; j++) {
		        /* use fix distance to define xl, xr, yu, yd instead of searchquader */
		        xl[j]= xr[j]= yu[j]= yd[j] = 3.0;

	            zaehler2 = candsearch_in_pixrest (fb->buf[3]->targets[j], 
                    fb->buf[3]->num_targets[j], xn[j], yn[j],
					xl[j], xr[j], yu[j], yd[j], &philf[j]);
		        if(zaehler2>0 ) {
                    _ix = philf[j][0];
		            x2[j] = fb->buf[3]->targets[j][_ix].x;
                    y2[j] = fb->buf[3]->targets[j][_ix].y;
		        }
		    }
	        quali=0;

	        for (j = 0; j < fb->num_cams; j++) {
		        if (x2[j] != -1e10 && y2[j] != -1e10) {
		        pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
		        }
		    }

	        if ( quali >= 2) {
        		X4 = X5; Y4 =Y5; Z4 = Z5;
		        invol=0; okay=0;

		        det_lsq_3d (Ex, I, G, ap, mmp, x2[0], y2[0], x2[1], y2[1], 
                    x2[2], y2[2], x2[3], y2[3], &X4, &Y4, &Z4);

		        /* volume check */
		        if ( X_lay[0] < X4 && X4 < X_lay[1] &&
		            Ymin < Y4 && Y4 < Ymax &&
		            Zmin_lay[0] < Z4 && Z4 < Zmax_lay[1]) {invol=1;}

        		/* displacement check */
		        if ( invol==1 &&
		            tpar.dvxmin < (X3-X4) && (X3-X4) < tpar.dvxmax &&
		            tpar.dvymin < (Y3-Y4) && (Y3-Y4) < tpar.dvymax &&
		            tpar.dvzmin < (Z3-Z4) && (Z3-Z4) < tpar.dvzmax ) 
                { 
                    okay=1;
                    
		            if (okay == 1) {
		                rr=1000000; dl=0;
		                acc=2*tpar.dacc; angle=2*tpar.dangle;
		                angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle, &acc);
		                dl=(sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3))
        			    +sqrt((X4-X3)*(X4-X3)+(Y4-Y3)*(Y4-Y3)+(Z4-Z3)*(Z4-Z3)))/2;

		                if ((acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10)) {
			                rr = (dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali+w[mm].freq);
			                curr_path_inf->decis[curr_path_inf->inlist]=rr;
			                curr_path_inf->linkdecis[curr_path_inf->inlist]=w[mm].ftnr;
			                curr_path_inf->inlist++;

			                if (tpar.add) {
                                ref_path_inf = &(fb->buf[3]->path_info[
                                    fb->buf[3]->num_parts]);
			                    ref_path_inf->x[0] = X4;
			                    ref_path_inf->x[1] = Y4;
			                    ref_path_inf->x[2] = Z4;
			                    ref_path_inf->prev = -1;
			                    ref_path_inf->next = -2;
			                    ref_path_inf->prio = 2;

                                _frame_parts = fb->buf[3]->num_parts;
                                ref_corres = &(fb->buf[3]->correspond[_frame_parts]);
                                ref_targets = fb->buf[3]->targets;
			                    for (j = 0; j < fb->num_cams; j++) {
				                    ref_corres->p[j]=-1;
                                    
				                    if(philf[j][0]!=-999) {
                                        _ix = philf[j][0];
                    				    ref_targets[j][_ix].tnr = _frame_parts;
				                        ref_corres->p[j] = _ix;
				                        ref_corres->nr = _frame_parts;
				                    }
			                    }
			                    fb->buf[3]->num_parts++;
                                zusatz++; 
                            }
			            }
		            }
                    okay=0; 
                }
		        invol=0;
	        }
	        quali=0;

	        /* end of creating new particle position */
	        /* *************************************************************** */
            
	        /* try to link if kk is not found/good enough and prev exist */
	        if ( curr_path_inf->inlist == 0 && curr_path_inf->prev >= 0 ) {
		        acc = 2*tpar.dacc;
                angle = 2*tpar.dangle;
		        if ( tpar.dvxmin < (X3-X1) && (X3-X1) < tpar.dvxmax &&
		            tpar.dvymin < (Y3-Y1) && (Y3-Y1) < tpar.dvymax &&
		            tpar.dvzmin < (Z3-Z1) && (Z3-Z1) < tpar.dvzmax ) 
                {
                    okay=1;
                    
		            if ( okay ==1 ) {
			            rr=1000000; quali=0;
			            quali=w[mm].freq;
			            angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);
			            dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			            +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

			            if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) ) {
			                rr = (dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
			                curr_path_inf->decis[curr_path_inf->inlist] = rr;
			                curr_path_inf->linkdecis[curr_path_inf->inlist] = w[mm].ftnr;
			                curr_path_inf->inlist++;
			            }
		            }
		        }
	        }
	        okay=0;

	        free(wn);
        } /* end of zaehler1-loop */
        
	    /* begin of inlist still zero */
	    if (tpar.add) {
	        if ( curr_path_inf->inlist == 0 && curr_path_inf->prev >= 0 ) {
                for (j = 0; j < fb->num_cams; j++) {
		            img_coord (X2, Y2, Z2, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], &yn[j]);
		            metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
		            x2[j]=-1e10;
                    y2[j]=-1e10;
                } 

    		    /* search for unused candidates in next time step */
		        for (j = 0; j < fb->num_cams; j++) {
		            /*use fix distance to define xl, xr, yu, yd instead of searchquader */
		            xl[j]= xr[j]= yu[j]= yd[j] = 3.0;
	                zaehler2 = candsearch_in_pixrest (fb->buf[2]->targets[j], 
                        fb->buf[2]->num_targets[j], xn[j], yn[j],
					    xl[j], xr[j], yu[j], yd[j], &philf[j]);
		            if(zaehler2 > 0) {
                        _ix = philf[j][0];
	    	            x2[j] = fb->buf[2]->targets[j][_ix].x;
                        y2[j] = fb->buf[2]->targets[j][_ix].y;
		            }
		        }
		        quali=0;

		        for (j = 0; j < fb->num_cams; j++) {
		            if (x2[j] !=-1e10 && y2[j] != -1e10) {
		                pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
		            }
		        }

		        if (quali>=2) {
		            X3 = X2; Y3 =Y2; Z3 = Z2;
		            invol=0; okay=0;
    
	    	        det_lsq_3d (Ex, I, G, ap, mmp,
		    	        x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X3, &Y3, &Z3);

		            /* in volume check */
		            if ( X_lay[0] < X3 && X3 < X_lay[1] &&
		                Ymin < Y3 && Y3 < Ymax &&
		                Zmin_lay[0] < Z3 && Z3 < Zmax_lay[1]) {invol=1;}

		            /* displacement check */
		            if ( invol==1 &&
		                tpar.dvxmin < (X2-X3) && (X2-X3) < tpar.dvxmax &&
		                tpar.dvymin < (Y2-Y3) && (Y2-Y3) < tpar.dvymax &&
		                tpar.dvzmin < (Z2-Z3) && (Z2-Z3) < tpar.dvzmax ) 
                    { 
                        okay=1;
                    
		                if (okay == 1) {
			                rr=1000000; dl=0;
			                acc=2*tpar.dacc;angle=2*tpar.dangle;
			                angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);
			                dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
			                +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

			                if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) ) {
			                    rr = (dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);
                                ref_path_inf = &(fb->buf[2]->path_info[
                                    fb->buf[2]->num_parts]);
			                    ref_path_inf->x[0] = X3;
			                    ref_path_inf->x[1] = Y3;
			                    ref_path_inf->x[2] = Z3;
			                    ref_path_inf->prev = -1;
			                    ref_path_inf->next = -2;
			                    ref_path_inf->prio = 2;

                                _frame_parts = fb->buf[2]->num_parts;
			                    curr_path_inf->decis[curr_path_inf->inlist] = rr;
			                    curr_path_inf->linkdecis[curr_path_inf->inlist] = \
                                    _frame_parts;
			                    curr_path_inf->inlist++;
                                
                                ref_corres = &(fb->buf[2]->correspond[_frame_parts]);
                                ref_targets = fb->buf[2]->targets;
			                    for (j = 0;j < fb->num_cams; j++) {
				                    ref_corres->p[j]=-1;
                                    
				                    if(philf[j][0]!=-999) {
                                        _ix = philf[j][0];
                    				    ref_targets[j][_ix].tnr = _frame_parts;
				                        ref_corres->p[j] = _ix;
				                        ref_corres->nr = _frame_parts;
				                    }
			                    }
			                    fb->buf[2]->num_parts++;
                                zusatz++;
			                }
		                } 
                        okay=0; 
                    }
		            invol=0;
		        } // if quali >= 2
            }
        }
	    /* end of inlist still zero */
	    /***********************************/

	    free(w);
	} /* end of h-loop */
    
    /* sort decis and give preliminary "finaldecis"  */
    for (h = 0; h < fb->buf[1]->num_parts; h++) {
        curr_path_inf = &(fb->buf[1]->path_info[h]);
        
	    if(curr_path_inf->inlist > 0 ) {
	        sort(curr_path_inf->inlist, &(curr_path_inf->decis),
                curr_path_inf->linkdecis);
      	    curr_path_inf->finaldecis = curr_path_inf->decis[0];
	        curr_path_inf->next = curr_path_inf->linkdecis[0];
	    }
	}

    /* create links with decision check */
    for (h = 0;h < fb->buf[1]->num_parts; h++) {
        curr_path_inf = &(fb->buf[1]->path_info[h]);

	    if(curr_path_inf->inlist > 0 ) {
            ref_path_inf = &(fb->buf[2]->path_info[curr_path_inf->next]);
            
	        if (ref_path_inf->prev == -1) {	
	            /* best choice wasn't used yet, so link is created */
                ref_path_inf->prev = h; 
            } else {
	            /* best choice was already used by mega[2][mega[1][h].next].prev */
	            /* check which is the better choice */
	            if ( fb->buf[1]->path_info[ref_path_inf->prev].finaldecis > \
                    curr_path_inf->finaldecis) 
                {
		            /* remove link with prev */
		            fb->buf[1]->path_info[ref_path_inf->prev].next= -2;
                    ref_path_inf->prev = h; 
		        } else {
		            curr_path_inf->next = -2;
	            }
	        }
        }
        if (curr_path_inf->next != -2 ) count1++;
    } 
    /* end of creation of links with decision check */
    /* ******** Draw links now ******** */
    m1_tr = 0;
    
    if (display) {
        for (h = 0; h < fb->buf[1]->num_parts; h++) {
            curr_path_inf = &(fb->buf[1]->path_info[h]);
            curr_corres = &(fb->buf[1]->correspond[h]);
            ref_corres = &(fb->buf[2]->correspond[curr_path_inf->next]);
            
            if (curr_path_inf->next != -2 ) {
                strcpy(buf,"");
                sprintf(buf ,"green");
	
                for (j = 0; j < fb->num_cams; j++) {
                    if (curr_corres->p[j] > 0 && ref_corres->p[j] > 0) {
                        flag_m_tr=1;  
                        xp[j] = curr_targets[j][curr_corres->p[j]].x;
               		    yp[j] = curr_targets[j][curr_corres->p[j]].y;
               		    xc[j] = fb->buf[2]->targets[j][ref_corres->p[j]].x;
               		    yc[j] = fb->buf[2]->targets[j][ref_corres->p[j]].x;
               		    predict (xp[j], yp[j], xc[j], yc[j], &xn[j], &yn[j]);
                        
               		    if ( ( fabs(xp[j]-zoom_x[j]) < imx/(2*zoom_f[j]))
                            && (fabs(yp[j]-zoom_y[j]) < imy/(2*zoom_f[j])))
                        {
	                        strcpy(val,"");
                           	sprintf(val ,"orange");

                        	intx0 = (int)(imx/2+zoom_f[j]*(xp[j]-zoom_x[j]));
	                        inty0 = (int)(imy/2+zoom_f[j]*(yp[j]-zoom_y[j]));
                            intx1 = (int)(imx/2+zoom_f[j]*(xc[j]-zoom_x[j]));
	                        inty1 = (int)(imy/2+zoom_f[j]*(yc[j]-zoom_y[j]));
	                        intx2 = (int)(imx/2+zoom_f[j]*(xn[j]-zoom_x[j]));
	                        inty2 = (int)(imy/2+zoom_f[j]*(yn[j]-zoom_y[j]));

	                        intx0_tr[j][m1_tr]=intx0;
	                        inty0_tr[j][m1_tr]=inty0;
	                        intx1_tr[j][m1_tr]=intx1;
	                        inty1_tr[j][m1_tr]=inty1;
	                        intx2_tr[j][m1_tr]=intx2;
	                        inty2_tr[j][m1_tr]=inty2;
	                        pnr1_tr[j][m1_tr]=-1;
	                        pnr2_tr[j][m1_tr]=-1;
	                        pnr3_tr[j][m1_tr]=-1;
		
	                        if (curr_path_inf->finaldecis > 0.2) {
	            	            pnr1_tr[j][m1_tr] = h;
		                        pnr2_tr[j][m1_tr] = curr_path_inf->next;
		                        pnr3_tr[j][m1_tr] = curr_path_inf->finaldecis;
		                    }
                        }
                    }

                    if (flag_m_tr==0)  {
                        intx0_tr[j][m1_tr]=0;
    	                inty0_tr[j][m1_tr]=0;
	                    intx1_tr[j][m1_tr]=0;
	                    inty1_tr[j][m1_tr]=0;
	                    intx2_tr[j][m1_tr]=0;
	                    inty2_tr[j][m1_tr]=0;
	                    pnr1_tr[j][m1_tr]=-1;
	                    pnr2_tr[j][m1_tr]=-1;
	                    pnr3_tr[j][m1_tr]=-1; 
                    }
                    flag_m_tr=0;
                }
                m1_tr++;
            }
        }
    }
    /* ******** End of Draw links now ******** */
    sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
        step, fb->buf[1]->num_parts, fb->buf[2]->num_parts, count1, 
        fb->buf[1]->num_parts - count1, zusatz);

    /* for the average of particles and links */
    npart = npart + fb->buf[1]->num_parts;
    nlinks = nlinks + count1;

    fb_next(fb);
    fb_write_frame_from_start(fb, step);
    if(step < seq_last - 2) { fb_read_frame_at_end(fb, step + 3); }
} /* end of sequence loop */

int trackcorr_c_finish(int step)
{
  /* average of all steps */
  npart /= (seq_last-seq_first);
  nlinks /= (seq_last-seq_first);
  printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
	  npart, nlinks, npart-nlinks);

  fb_next(fb);
  fb_write_frame_from_start(fb, step);
  
  fb_free(fb);
  free(fb);
    
  /* reset of display flag */
  display = 1;
}

/*     track backwards */

int trackback_c ()
{
    char  buf[256];
    int i, j, h, k, step, okay=0, invol=0;
    int zaehler1, philf[4][4];
    int count1=0, count2=0, zusatz=0;
    int quali=0;
    double x2[4], y2[4], angle, acc, lmax, dl;
    double xr[4], xl[4], yd[4], yu[4];
    double X1, Y1, Z1, X0, Y0, Z0, X2, Y2, Z2, X5, Y5, Z5;
    double X3, Y3, Z3, X4, Y4, Z4, X6, Y6, Z6;
    double xn[4], yn[4];
    double rr, Ymin=0, Ymax=0;
    double npart=0, nlinks=0;
    foundpix *w, p16[16];

    display = 1; //atoi(argv[1]);
    /* read data */
    readseqtrackcrit ();
    
    /*Alloc space, if checkflag for mega, c4, t4 is zero */
    if (!trackallocflag) {
        allocate_tracking_structs(t4, c4, mega, n_img, M);
    }

    /*load again first data sets*/
    step = seq_last;
    read_ascii_datanew(step);
    rotate_dataset();
    read_ascii_datanew(step-1);
    rotate_dataset();
    read_ascii_datanew(step-2);
    rotate_dataset();
    read_ascii_datanew(step-3);
    
    lmax=sqrt((tpar.dvxmin-tpar.dvxmax)*(tpar.dvxmin-tpar.dvxmax)
        +(tpar.dvymin-tpar.dvymax)*(tpar.dvymin-tpar.dvymax)
        +(tpar.dvzmin-tpar.dvzmax)*(tpar.dvzmin-tpar.dvzmax));

    volumedimension (&X_lay[1], &X_lay[0], &Ymax, &Ymin, &Zmax_lay[1], &Zmin_lay[0]);

    /* sequence loop */
    for (step = seq_last-1; step > seq_first; step--) {
        sprintf (buf, "Time step: %d, seqnr: %d, Particle info:", step- seq_first, step);
        
        for (h=0; h<m[1]; h++) {
            if (mega[1][h].next>=0 && mega[1][h].prev==-1) {
                X1=Y1=Z1=X0=Y0=Z0=X2=Y2=Z2=X5=Y5=Z5=X3=Y3=Z3=X4=Y4=Z4=X6=Y6=Z6=-999;
                
                mega[1][h].inlist=0;
                for (i=0; i<16;i++) {
                    p16[i].ftnr=-1;
                    p16[i].freq=0;
                    for(j=0;j<n_img;j++) p16[i].whichcam[j] =0;
                }
                
                /* 3D-position */
                X1=mega[1][h].x[0];
                Y1=mega[1][h].x[1];
                Z1=mega[1][h].x[2];
                
                /* use information from previous to locate new search position
                and to calculate values for search area */
                X0=mega[0][mega[1][h].next].x[0];
                Y0=mega[0][mega[1][h].next].x[1];
                Z0=mega[0][mega[1][h].next].x[2];
                X2=2*X1-X0;
                Y2=2*Y1-Y0;
                Z2=2*Z1-Z0;

                for (j=0; j<n_img; j++) {   
                    img_coord (X2, Y2, Z2, Ex[j],I[j], G[j], ap[j], mmp, &xn[j], &yn[j]);
                    metric_to_pixel (xn[j], yn[j], imx,imy, pix_x,pix_y, &xn[j], &yn[j], chfield);
                }

                /* calculate searchquader and reprojection in image space */
                searchquader(X2, Y2, Z2, &xr, &xl, &yd, &yu);

                /* search in pix for candidates in next time step */
                for (j=0; j<n_img; j++) {
                    zaehler1 = candsearch_in_pix (t4[2][j], nt4[2][j], xn[j], yn[j],
                    xl[j], xr[j], yu[j], yd[j], &philf[j]);

                    for(k=0; k<4; k++) {
                        if( zaehler1>0) {
                            if (philf[j][k] == -999){
                                p16[j*4+k].ftnr=-1;
                            } else {
                                p16[j*4+k].ftnr=t4[3][j][philf[j][k]].tnr;
                                p16[j*4+k].whichcam[j]=1;
                            }
                        }
                    }
                }

                /* fill and sort candidate struct */
                sortwhatfound(&p16, &zaehler1);
                w = (foundpix *) calloc (zaehler1, sizeof (foundpix));

                /*end of candidate struct */
                if (zaehler1 > 0) count2++;
                for (i=0; i<zaehler1;i++) {
                    w[i].ftnr = p16[i].ftnr;
                    w[i].freq = p16[i].freq;
                    for (j=0; j<n_img; j++) w[i].whichcam[j] = p16[i].whichcam[j];
                }

                if (zaehler1 > 0) for (i=0; i<zaehler1;i++) {
                    X3=mega[2][w[i].ftnr].x[0];
                    Y3=mega[2][w[i].ftnr].x[1];
                    Z3=mega[2][w[i].ftnr].x[2];

                    okay = 0;
                    acc = 2*tpar.dacc;
                    angle = 2*tpar.dangle;
                    rr = 1000000; quali = 0; dl = 0;
                    
                    /* displacement check */
                    if ( tpar.dvxmin < (X1-X3) && (X1-X3) < tpar.dvxmax &&
                        tpar.dvymin < (Y1-Y3) && (Y1-Y3) < tpar.dvymax &&
                        tpar.dvzmin < (Z1-Z3) && (Z1-Z3) < tpar.dvzmax ) 
                    {
                        okay=1;
                        /* end displacement check */

                        if ( okay ==1 ) {
                            dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
                                +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

                            quali=w[i].freq;
                            angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);

                            /* *********************check link *****************************/
                            if ((acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10)) {
                                rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/quali;

                                mega[1][h].decis[mega[1][h].inlist]=rr;
                                mega[1][h].linkdecis[mega[1][h].inlist]=w[i].ftnr;
                                mega[1][h].inlist++;
                            }
                        }
                    } okay=0;
                }

                free(w);
                /******************/
                quali=0;

                /* reset img coord because of n_img smaller 4 */
                for (j=0;j<4;j++) { x2[j]=-1e10; y2[j]=-1e10;}

                /* if old wasn't found try to create new particle position from rest */
                if (tpar.add) {
                    if ( mega[1][h].inlist == 0) {
                        for (j=0;j<n_img;j++) {
                            /* use fix distance to define xl, xr, yu, yd instead of searchquader */
                            xl[j]= xr[j]= yu[j]= yd[j] = 3.0;

                            zaehler1 = candsearch_in_pixrest (t4[2][j], nt4[2][j], xn[j], yn[j],
                            xl[j], xr[j], yu[j], yd[j], &philf[j]);
                            if(zaehler1>0 ) { x2[j]=t4[2][j][philf[j][0]].x; y2[j]=t4[2][j][philf[j][0]].y; }
                        }

                        for (j=0;j<n_img;j++) {
                            if (x2[j] !=-1e10 && y2[j] != -1e10) {
                                pixel_to_metric (x2[j],y2[j], imx,imy, pix_x,pix_y, &x2[j],&y2[j], chfield); quali++;
                            }
                        }

                        if (quali>=2) {
                            X3 = X2; Y3 =Y2; Z3 = Z2;
                            invol=0; okay=0;

                            det_lsq_3d (Ex, I, G, ap, mmp,
                                x2[0], y2[0], x2[1], y2[1], x2[2], y2[2], x2[3], y2[3], &X3, &Y3, &Z3);

                            /* volume check */
                            if ( X_lay[0] < X3 && X3 < X_lay[1] && Ymin < Y3 && Y3 < Ymax &&
                                Zmin_lay[0] < Z3 && Z3 < Zmax_lay[1]) {invol=1;}

                            okay=0; acc=2*tpar.dacc;angle=2*tpar.dangle;rr=1000000; dl=0;

                            /* displacement check */
                            if ( invol==1 &&
                                tpar.dvxmin < (X1-X3) && (X1-X3) < tpar.dvxmax &&
                                tpar.dvymin < (Y1-Y3) && (Y1-Y3) < tpar.dvymax &&
                                tpar.dvzmin < (Z1-Z3) && (Z1-Z3) < tpar.dvzmax ) 
                            { 
                                okay=1;
                                /* end displacement check */

                                if ( okay ==1 ) {
                                    dl=(sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0)+(Z1-Z0)*(Z1-Z0))
                                        +sqrt((X1-X3)*(X1-X3)+(Y1-Y3)*(Y1-Y3)+(Z1-Z3)*(Z1-Z3)))/2;

                                    angle_acc(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, &angle, &acc);

                                    if ( (acc<tpar.dacc && angle<tpar.dangle) || (acc<tpar.dacc/10) ) {
                                        rr =(dl/lmax+acc/tpar.dacc + angle/tpar.dangle)/(quali);

                                        mega[2][m[2]].x[0]= X3;
                                        mega[2][m[2]].x[1]= Y3;
                                        mega[2][m[2]].x[2]= Z3;
                                        mega[2][m[2]].prev= -1;
                                        mega[2][m[2]].next= -2;
                                        mega[2][m[2]].prio= 2;
                                        mega[1][h].decis[mega[1][h].inlist]=rr;
                                        mega[1][h].linkdecis[mega[1][h].inlist]=m[2];
    
                                        for (j=0;j<n_img;j++) {
                                            c4[2][m[2]].p[j]=-1;
                                            if(philf[j][0]!=-999) {
                                                t4[2][j][philf[j][0]].tnr=m[2];
                                                c4[2][m[2]].p[j]= philf[j][0];
                                                c4[2][m[2]].nr=m[2];
                                            }
                                        }
                                        mega[1][h].inlist++;
                                        m[2]++;
                                    }
                                }
                                okay = 0;
                            }
                            invol=0;
                        }
                    }
                } /* end of if old wasn't found try to create new particle position from rest */
            }
        } /* end of h-loop */

        /* sort decis  */
        for (h=0;h<m[1];h++) {
            if(mega[1][h].inlist > 0 ) {
                sort(mega[1][h].inlist, &mega[1][h].decis, &mega[1][h].linkdecis); 
            }
        }

        /* create links with decision check */
        count1=0; zusatz=0;
        for (h=0;h<m[1];h++) {
            if (mega[1][h].inlist > 0 ) {
            /* if old/new and unused prev == -1 and next == -2 link is created */
                if ( mega[2][mega[1][h].linkdecis[0]].prev == -1 && mega[2][mega[1][h].linkdecis[0]].next == -2 ) {
                    mega[1][h].finaldecis=mega[1][h].decis[0];
                    mega[1][h].prev=mega[1][h].linkdecis[0];
                    mega[2][mega[1][h].prev].next=h;
                    zusatz++;
                }

                /* old which link to prev has to be checked */
                if ( mega[2][mega[1][h].linkdecis[0]].prev != -1 && mega[2][mega[1][h].linkdecis[0]].next == -2 ) {
                    X0=mega[0][mega[1][h].next].x[0];
                    Y0=mega[0][mega[1][h].next].x[1];
                    Z0=mega[0][mega[1][h].next].x[2];

                    X1=mega[1][h].x[0];
                    Y1=mega[1][h].x[1];
                    Z1=mega[1][h].x[2];

                    X3=mega[2][mega[1][h].linkdecis[0]].x[0];
                    Y3=mega[2][mega[1][h].linkdecis[0]].x[1];
                    Z3=mega[2][mega[1][h].linkdecis[0]].x[2];

                    X4=mega[3][mega[2][mega[1][h].linkdecis[0]].prev].x[0];
                    Y4=mega[3][mega[2][mega[1][h].linkdecis[0]].prev].x[1];
                    Z4=mega[3][mega[2][mega[1][h].linkdecis[0]].prev].x[2];

                    X5=0.5*(5.0*X3-4.0*X1+X0);
                    Y5=0.5*(5.0*Y3-4.0*Y1+Y0);
                    Z5=0.5*(5.0*Z3-4.0*Z1+Z0);

                    acc=2*tpar.dacc;angle=2*tpar.dangle;
                    angle_acc(X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, &angle, &acc);
                    if ( (acc<tpar.dacc && angle<tpar.dangle) ||  (acc<tpar.dacc/10) ) {
                        mega[1][h].finaldecis=mega[1][h].decis[0];
                        mega[1][h].prev=mega[1][h].linkdecis[0];
                        mega[2][mega[1][h].prev].next=h;
                        zusatz++;
                    }
                }
            }

            if (mega[1][h].prev != -1 ) count1++;
        } /* end of creation of links with decision check */

        sprintf (buf, "step: %d, curr: %d, next: %d, links: %d, lost: %d, add: %d",
        step, m[1], m[2], count1, m[1]-count1, zusatz);

        /* for the average of particles and links */
        npart = npart + m[1];
        nlinks = nlinks + count1;

        rotate_dataset();
        write_ascii_datanew(step);
        write_addedback(step);
        if(step> seq_first+2) { read_ascii_datanew(step-3); }

    } /* end of sequence loop */

    /* average of all steps */
    npart /= (seq_last-seq_first-1);
    nlinks /= (seq_last-seq_first-1);
    printf ("Average over sequence, particles: %5.1f, links: %5.1f, lost: %5.1f\n",
    npart, nlinks, npart-nlinks);

    rotate_dataset();
    write_ascii_datanew(step);
    write_addedback(step);

    /* reset of display flag */
    display = 1;
}

