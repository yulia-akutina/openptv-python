#include "ptv.h"

int right_click_proc_c(int in_x, int in_y,int img_number)
{
	
	int     i, j, click_x, click_y, n, epi_vector_index;
	double  x, y;
	double  xa12, xb12, ya12, yb12;
	int     k, pt1, intx1, inty1, count, intx2, inty2, pt2;
	candidate cand[maxcand];
	
	click_x = in_x;
	click_y = in_y;

	n = img_number;
	epi_vector_index=0;
	 
    /* get geometric coordinates of nearest point in img[n] */
    x = (float) (click_x - imx/2)/zoom_f[n] + zoom_x[n];
    y = (float) (click_y - imy/2)/zoom_f[n] + zoom_y[n];
	pixel_to_metric (x,y, imx,imy, pix_x,pix_y, &x,&y, chfield);
    x -= I[n].xh; y -= I[n].yh;
    correct_brown_affin (x, y, ap[n], &x, &y);
    k = nearest_neighbour_geo (geo[n], num[n], x, y, 0.05);
    if (k == -999)
	{    
		sprintf (buf, "no point near click coord ! Click again!");  puts (buf);
		return -1;
		//Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		//Tcl_Eval(interp, ".text delete 2");
		//Tcl_Eval(interp, ".text insert 2 $tbuf");
		//return TCL_OK;
	}
    pt1 = geo[n][k].pnr;

    intx1 = (int) ( imx/2 + zoom_f[n] * (pix[n][pt1].x-zoom_x[n]));
    inty1 = (int) ( imy/2 + zoom_f[n] * (pix[n][pt1].y-zoom_y[n]));
    //TODO: return intx1 intx2
	//drawcross (interp, intx1, inty1, cr_sz+2, n, "BlueViolet");


    sprintf (buf, "%d %d %d %d %d\n", pt1, pix[n][pt1].nx, pix[n][pt1].ny,
		pix[n][pt1].n, pix[n][pt1].sumg);  puts (buf);
                           
    for (i=0; i<n_img; i++)   if (i != n)
	{
		/* calculate epipolar band in img[i] */
		epi_mm (geo[n][k].x,geo[n][k].y,
			Ex[n],I[n], G[n], Ex[i],I[i], G[i], mmp,
			&xa12, &ya12, &xb12, &yb12);
           
		/* search candidate in img[i] */
		printf("\ncandidates in img: %d\n", i);
		find_candidate_plus_msg (geo[i], pix[i], num[i],
			xa12, ya12, xb12, yb12, eps0,
			pix[n][pt1].n, pix[n][pt1].nx, pix[n][pt1].ny,
			pix[n][pt1].sumg, cand, &count, i);
		  
		distort_brown_affin (xa12,ya12, ap[i], &xa12,&ya12);
		distort_brown_affin (xb12,yb12, ap[i], &xb12,&yb12);
		xa12 += I[i].xh; ya12 += I[i].yh;
		xb12 += I[i].xh; yb12 += I[i].yh;
		metric_to_pixel (xa12, ya12, imx,imy, pix_x,pix_y,
			&xa12, &ya12, chfield);
		metric_to_pixel (xb12, yb12, imx,imy, pix_x,pix_y,
			&xb12, &yb12, chfield);
		
		epi_x1vector[epi_vector_index]=(int) ( imx/2 + zoom_f[i] * (xa12 - zoom_x[i]));
		epi_y1vector[epi_vector_index]=(int) ( imy/2 + zoom_f[i] * (ya12 - zoom_y[i]));
		epi_x2vector[epi_vector_index] = (int) ( imx/2 + zoom_f[i] * (xb12 - zoom_x[i]));
		epi_y2vector[epi_vector_index] = (int) ( imy/2 + zoom_f[i] * (yb12 - zoom_y[i]));
		counts[epi_vector_index]=count;
		sprintf (buf, "count is %d\n", count);  puts (buf);
		
		
		  //intx1 = 
          //inty1 = 
          //intx2 = 
          //inty2 = 

          //if ( n == 0 ) sprintf( val,"yellow");
          //if ( n == 1 ) sprintf( val,"green");
          //if ( n == 2 ) sprintf( val,"red");
          //if ( n == 3 ) sprintf( val,"blue");
		  //TODO return vector points and color
          //drawvector ( interp, intx1, inty1, intx2, inty2, 1, i, val);
		  
		for (j=0; j<count; j++)
		{

			pt2 = cand[j].pnr;
			points_x[epi_vector_index][j] = (int) ( imx/2 + zoom_f[i] * (pix[i][pt2].x - zoom_x[i]));
			points_y[epi_vector_index][j]= (int) ( imy/2 + zoom_f[i] * (pix[i][pt2].y - zoom_y[i]));
			  //TODO return orange crosses
			  //drawcross (interp, intx2, inty2, cr_sz+2, i, "orange");
		}

		epi_vector_index++;
		  
	}
	return 1;
}