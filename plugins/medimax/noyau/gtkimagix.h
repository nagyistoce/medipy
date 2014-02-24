/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		gtkimagix.h
***
***	project:	Imagix 2.01
***			
***
***	description:	Header file for project imagix (2.01)
***	
***	
***	Copyright (c) 1993-2000, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Fev 2003
***
***---------------------------------------------------------------------*/

#ifndef __NOGTKIMAGIX_INCL
#define __NOGTKIMAGIX_INCL

#define	POSx_WARN	650
#define	POSy_WARN	800
#define	POSx_MESG	20
#define	POSy_MESG	800	
#define	POSx_WNDMSG 	700
#define	POSy_WNDMSG	500
#define	POSx_QUEST	41
#define	POSy_QUEST	42
#define	POSx_ERR	550
#define	POSy_ERR	350
#define	POSx_ERRFAT	600
#define	POSy_ERRFAT	400
#define	POSx_INFO	120
#define	POSy_INFO	180
#define	POSx_FILE	600
#define	POSy_FILE	30
#define	POSx_SAVE	610
#define	POSy_SAVE	40
#define	POSx_ANY	100
#define	POSy_ANY	160
#define	POSx_WND	130
#define	POSy_WND	190
#define	POSx_POS	130
#define	POSy_POS	190
#define	POSx_PLACE	130
#define	POSy_PLACE	190
#define	POSx_QCM	130
#define	POSy_QCM	140
#define	POSx_LST	200
#define	POSy_LST	140
#define	POSx_FLOAT	200
#define	POSy_FLOAT	140
#define	POSy_EXIST	200
#define	POSx_EXIST	140
#define	POSx_OKCANCEL	140
#define	POSy_OKCANCEL	40

#define	POSx_VALUE	POSx_ANY
#define	POSy_VALUE	POSy_ANY


/* ------------------------------------------------------------------------- */
/* ---- macros ------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */
#define	GET_MOUSE2ACTION(fp,ap,ff,af,fr,ar)  printf("ERROR GET_MOUSE2ACTION: \n")
#define	GETmouseOKCANCEL(x,y,e)              printf("ERROR GETmouseOKCANCEL: \n")
#define	XGET_MOUSE2ACTION(fp,ap,ff,af,fr,ar) printf("ERROR XGET_MOUSE2ACTION: \n")
#define XGETmouseOKCANCEL(x,y,e)             printf("ERROR XGETmouseOKCANCEL: \n")
#define	SAVE_FILE(m,i,e)                     (char*)printf("ERROR SAVE_FILE: \n")
#define	GET_WND3D(m)    	                   printf("ERROR GET_WND3D: \n")
#define	GET_WND(x,y,m)		                   printf("ERROR GET_WND: \n")
#define	GET_PLACE3D(m)                       printf("ERROR GET_PLACE3D: \n")
#define	GET_PLACE(m)                     printf("ERROR GET_PLACE: \n")
#define	GET_VALUE(t,m,i,e)                   (t)printf("ERROR GET_VALUE: \n")
#define	GETV_QCM(m,v)                        printf("ERROR GETV_QCM: \n")
#define	GET_FLOAT(m,i,e)                     (float)printf("ERROR GET_FLOAT: %d\n",e)
#define	GET_INT(m,i,e)                        printf("ERROR GET_INT: \n")
#define	GET_DOUBLE(m,i,e)                    (double)printf("ERROR GET_DOUBLE: \n")
#define	GET_OKCANCEL(x,y,e)                  printf("ERROR vGET_OKCANCEL: \n")
#define	GET_FILE_IMG(x,y,i,n,e)	             (char*)printf("ERROR GET_FILE_IMG: \n")
#define	GET_FILE(i,e)                        (char*)printf("ERROR GET_FILE_IMG: \n")
#define	GET_EXIST(m,e)                       printf("ERROR GET_FILE: \n")
#define GET_ADDORREPLACE()	                 printf("ERROR GET_ADDORREPLACE: \n")
#define GET_MOUSEPOS_3D(auto)                (int*)printf("ERROR GET_MOUSEPOS_3D: \n")
#define GET_MOUSEPOS(auto)                   (int*)printf("ERROR GET_MOUSEPOS: \n")
#define GET_ARRAY(x,y,m,e)                   printf("ERROR GET_ARRAY: \n")
#define	GET_CUTOFF(e)                        printf("ERROR GET_CUTOFF: \n")
#define	GET_CUTOFF3D(e)                      printf("ERROR GET_CUTOFF3D: \n")
#define GET_CnI(e)                           printf("ERROR GET_CnI: \n")
#define	GET_CONTRASTnINTENSITY(e)   	     printf("ERROR GET_CONTRASTnINTENSITY: \n")
#define	GET_THRESHOLD(i, e)         	     printf("ERROR GET_THRESHOLD: \n")
#define	GET_THRESHOLD3D(i, e)       	     printf("ERROR GET_THRESHOLD3D: \n")
#define	GET_LIST_MANAGER(x,y,t,dat) 	     printf("ERROR GET_LIST_MANAGER: \n")
#define	GET_HSLIDER(x,y,m,dat)      	     printf("ERROR GET_HSLIDER: \n")

#ifdef __GTK_GUI
#define XGET_MOUSEPOS_3D(e)     	  (int*)printf("ERROR XGET_MOUSEPOS_3D: \n")
#define	XGET_THRESHOLD3D(i, e)  	  (int*)printf("ERROR XGET_THRESHOLD3D: \n")

#else // not __GTK_GUI

//#define	XGET_EXIST(x,y,m,e)     	  printf("ERROR GET: \n")
#define	XGET_EXIST(x,e)     	  	  printf("ERROR GET: \n")
#define	XSAVE_FILE(x,y,m,i,e)   	  (char*)printf("ERROR GET: \n")
#define	XGET_WND(x,y,m)         	  printf("ERROR GET: \n")
#define	XGET_WND3D(x,y,m)       	  printf("ERROR GET: \n")
#define	XGET_PLACE3D(x,y,m)     	  printf("ERROR GET: \n")
#define	XGET_PLACE(x,y,m)       	  printf("ERROR GET: \n")
#define	XGET_VALUE(t,x,y,m,i,e) 	  printf("ERROR GET: \n")
#define	XGETV_QCM(t,x,y,m,v)    	  printf("ERROR GET: \n")
#define	XGET_FLOAT(t,x,y,m,i,e) 	  printf("ERROR GET: \n")
#define	XGET_OKCANCEL(x,y,e)    	  printf("ERROR GET: \n")
#define	XGET_FILE_IMG(x,y,i,j,e)	  (char*)printf("ERROR GET: \n")
#define	XGET_FILE(x,y,i,e)      	  (char*)printf("ERROR GET: \n")
#define XGET_MOUSEPOS(e)        	  (int*)printf("ERROR GET: \n")
#define XGET_MOUSEPOS_3D(e)     	  (int*)printf("ERROR GET: \n")
#define	XGET_THRESHOLD(i, e)    	  printf("ERROR GET: \n")
#define	XGET_CUTOFF(e)          	  printf("ERROR GET: \n")
#define	XGET_CUTOFF3D(e)        	  printf("ERROR GET: \n")
#define XGET_CnI(e)             	  printf("ERROR GET: \n")
#define	XGET_CONTRASTnINTENSITY(e) 	  printf("ERROR GET: \n")
#define	XGET_LIST_MANAGER(x,y,t,dat)      printf("ERROR GET: \n")
#define	XGET_THRESHOLD3D(i, e)  	  (int*)printf("ERROR XGET_THRESHOLD3D: \n")
#endif //__GTK_GUI		
			
#define PUT_MESG(m)             printf("MESSAGE:\"%s\"\n",m);
#define PUT_WNDMSG(m)           printf("MESSAGE:\"%s\"\n",m);
#define PUT_WARN(m)             printf("WARNING:\"%s\"\n",m);

#define PUT_ERRFAT(m)           fprintf(stderr,"ERROR:\"%s\"\n",m);exit(1);
#define PUT_ERR(m)              fprintf(stderr,"ERROR:\"%s\"\n",m);exit(1);
#define PUT_ERROR(m)            fprintf(stderr,"ERROR:\"%s\"\n",m);
#define PUT_INFO(m)             printf("INFO:\"%s\"\n",m);

#if defined (LINUX) || defined (WIN32)
#define	PUT_VALUE(t,x,y,m)		printf("ERROR PUT_VALUE: \n");
#else
#define	PUT_VALUE(t,x,y,m)		printf("ERROR PUT_VALUE: \n");
#endif

#define	PUTI(s1,i,s2)			spI_xprintf(s1,i,s2)
#define	PUTF(s1,f,s2)			spF_xprintf(s1,f,s2)
#define	PUTA(s1,f,s2)			spA_xprintf(s1,f,s2)

#define CLEAR_WNDMSG()          printf("ERROR CLEAR_WNDMSG: \n");

typedef struct s_Action {
void	(*function)();	/* Fonction associe a l'action de la souris */
void	*data_action;	/* Data transmise a la fonction (function) */
int	x2d;				/* position 2d x et y de la souris sur l'image pos2d */
int	y2d;
int	pos2d;
int	evtxbuttonstate;	/* Etats de button */ 
int	x3d;		/* position 3d x, y  et z de la souris dans l'images pos3d */
int	y3d;
int	z3d;
int	pos3d;
} Action,*ptr_Action;

/*
#ifdef WIN32
void bzero(void *s, int n);
#endif //WIN32
*/

#endif  /* End of gtkimagix.h 2.01 */

