/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!*************************************************************************
**
**	\file		misc.c
**
**	project:	Imagix 1.01
**
**
**	\brief description:	Fonction necessaire pour la compilation du module
**
**    Ces fonctions ne sont normallement jamais appele lors d'une utilisation
**    normale
**
**
**
***************************************************************************/

#include <config.h>

#include <stdlib.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/io/imx_head.h"
#include "noyau/imx_3d.h"

#if 0

/* --             ------------------------------------------------
**	Gestion ecriture dans la boite de dialogue
**      et dans un fichier lotus
**
**      Un _ a la fin du premier parametre indique la supression 
**      du CR-LF pour l'ecriture dans le fichier LOTUS
**
**      Le troisieme parametre correspond a la variable a afficher
**
**      Un _ a la fin du troisieme parametre indique la supression 
**      du CR-LF pour l'ecriture dans la boite de dialogue
*/


void spI_xprintf(char *s_src1, int i, char *s_src2)
{
  char *s,*s1,*s2;
  int alaligne;
  int error_no=1;

  alaligne=0;
//  s1=CALLOC(strlen(s_src1),char);
  s1=CALLOC(strlen(s_src1)+1,char); //correction Berst
  strcpy(s1,s_src1);
//  s2=CALLOC(strlen(s_src2),char);
  s2=CALLOC(strlen(s_src2)+1,char); //correction Berst
  strcpy(s2,s_src2);

  s=CALLOC(strlen(s1)+strlen(s2)+sizeof(int)+16,char);
  if(s==NULL) {
    if(MESG_DEBUG) printf(" errno=%d\n",error_no); }
  if(s1[strlen(s1)-1]=='_') 
    {
    alaligne=1;
    s1[strlen(s1)-1]=0;
    }
  if(s2[strlen(s2)-1]=='_')
    {
    s2[strlen(s2)-1]=0;
    sprintf(s,"%s %d %s",s1,i,s2);
    }
  else 
    sprintf(s,"%s %d %s\n",s1,i,s2);

  PUT_MESG(s);
  if(_dialog_messages_status & (1l<<0)) { // Lotus file is active. 
//    Test d'ecriture dans un fichier lotus (alaligne=0) 
    if(_fp_lotus)
      {
      if(alaligne==1)
        sprintf(s,"%d ",i);
      else 
        sprintf(s,"%d\n",i);
      fprintf(_fp_lotus,"%s",s);
      fflush(_fp_lotus);
      }
  }

  free(s);
  free(s1);
  free(s2);
}

void spF_xprintf(char *s_src1, float f, char *s_src2)
{ 
  char *s,*s1,*s2;
  int alaligne;
  int error_no=1;

  alaligne=0;
  s1=CALLOC(strlen(s_src1),char);
  strcpy(s1,s_src1);
  s2=CALLOC(strlen(s_src2),char);
  strcpy(s2,s_src2);

  s=CALLOC(strlen(s1)+strlen(s2)+sizeof(float)+16,char);
  if(s==NULL) {
    if(MESG_DEBUG) printf(" errno=%d\n",error_no); }
  if(s1[strlen(s1)-1]=='_') 
    {
    alaligne=1;
    s1[strlen(s1)-1]=0;
    }
  if(s2[strlen(s2)-1]=='_')
    {
    s2[strlen(s2)-1]=0;
    sprintf(s,"%s %f %s",s1,f,s2);
    }
  else
    sprintf(s,"%s %f %s\n",s1,f,s2);

  PUT_MESG(s);
  if(_dialog_messages_status & (1l<<0)) { // Lotus file is active. 
//    Test d'ecriture dans un fichier lotus (alaligne=0) 
    if(_fp_lotus)
      {
      if(alaligne==1)
        sprintf(s,"%f ",f);
      else 
        sprintf(s,"%f\n",f);
      fprintf(_fp_lotus,"%s",s);
      fflush(_fp_lotus);
      }
  }

  free(s);
  free(s1);
  free(s2);
}


/* --save_header() -------------------------------------------------------------
*/
/*!    Saves a header
**
**    Usage:    save_header(header, file)
**    \param    header : pointer to struct HEADER
**    \param    file : header file
**  
**    \retval 1 succesful, 0  error
**                     
*/

int     save_header(HEADER *header, char *file)
{
    char fileold[256];
    char filenew[256];
    FILE *fp;
    int i;

    if(header->lines == NULL) {
        PUT_ERR("[save_header] Invalid header !\n");
        return(0);
    }

    strcpy(filenew, file);

    if(header->exists) { /* file exists */
        strcat(filenew, ".tmp");
    }

    /* writes the new header */
    if((fp = fopen(filenew, "w")) == NULL) {
		char e[256];
        sprintf(e, "[save_header] cannot create file \"%s\"\n", filenew);
		PUT_ERR( e);
        return(0);
    }

    for(i=0; i<header->numlines; i++) {
        fprintf(fp, "%s\n", header->lines[i]);
    }

    fclose(fp);

    if(header->exists) {
        strcpy(fileold, file);
        /* strcat(fileold, ".ipb"); */

        if(remove(fileold) == -1) {
			char e[256];
            sprintf(e, "[save_header] unable to remove file \"%s\"\n", fileold);
			PUT_ERR( e);
            return(0);
        }
        if(rename(filenew, fileold) == -1) {
			char e[256];
            sprintf(e, "[save_header] unable to rename file \"%s\" in \"%s\"\n", filenew, fileold);
			PUT_ERR( e);
            return(0);
        }
    }

    header->exists = 1;
    header->maxlines = 0;
    header->numlines = 0;
    FREE(header->lines);
    header->lines = NULL;

    return(1);
}


/* --put_header() --------------------------------------------------------------
*/
/*!    Puts a info in a header
**
**    Usage:    put_header(header, string, type, itemadr, pos)
**    \param    header : pointer to struct HEADER
**    \param    string : string to match
**    \param	type : type of data (ENTIER, REEL, STRING)
**    \param    itemadr : address of data item
**    \param    pos : position in the line where to add the info
**
**    \retval  1 successful, 0  error
**              
**                     
*/
int     put_header(HEADER *header, char *string, int type, void *itemadr, int pos)
{
    char *newline;
    int maxlines, numlines;
    char **lines;
    int i;
    char *p = NULL;
    int state;
    int numwords;
    int length;
    char value[64];

    maxlines = header->maxlines;
    numlines = header->numlines;
    lines = header->lines;

    if(lines == NULL) {
        PUT_ERR("[put_header] Invalid header !\n");
        return(0);
    }

    if((newline = (char *)malloc(256 * sizeof(char))) == NULL) {
        PUT_ERR( "[put_header] memory allocation error (1)!\n");
        return(0);
    }

    if(numlines > 0) { /* Search the line */
        for(i=0; i<numlines; i++) {
            if(strstr(lines[i], string) != NULL)
            {
                /* we've found the string */
                break;
            }
        }
    }
    else {
        i = numlines = 0;
    }

    /* Warning: don't modify 'i' it's used later */

    if(i == numlines) { /* case where the string was not found */
        strcpy(newline, string);

        numwords = 1;
        while(numwords++ < pos)
            strcat(newline, "0 ");
    }
    else {
        p = lines[i];

        if(strlen(p) > 192) {
            if((newline = (char *)realloc(newline, strlen(p) + 256)) == NULL) {
                PUT_ERR( "[put_header] reallocation failed !\n");
                return(0);
            }
        }   

        /* searches the '=' sign */
        while(*p != '\0' && *p != '=')
            p ++;

        if(!*p) { /* we have reached the end of the string */
            /* Supprime car erreur en ecrivant IPB3 
	    printf( "[put_header] no '=' in string\n"); */
            return(0);
        }

        state = 0;
        numwords = 0;

        /* skip the '=' sign */
        p ++;

        while(numwords < pos && *p != '\0') {
            if(*p == ' ')
			{
				state = 0;
				p++;
			}
            else {
                if(state == 0)
                {
                    state = 1;
                    numwords ++;
                }
				else
            		p ++;
            }

        }

        length = p - lines[i];

        strncpy(newline, lines[i], length);
        newline[length] = '\0';

        /* if there's not enough values, fill with zeroes */
        while(numwords++ < pos - 1)
            strcat(newline, "0 ");
    }

    switch(type)
    {
    case STRING:
        if(strlen((char *)itemadr) < 256)
            sprintf(value, "%s ", (char *)itemadr);
        else
            value[0] = '\0';
        break;

    case ENTIER:
        //marat	sprintf(value, "%d ", *(long *)itemadr);
        sprintf(value, "%ld ", *(long *)itemadr);
        
	break;

    case REEL:
        sprintf(value, "%f ", *(double *)itemadr);
        break;

    default:
        return(0);
    }

    /* adds the value to the string */
    strcat(newline, value);

    /* we have found the line */
    if(i < numlines) {
        /* skip the value */
        if (pos!= 0) {
            while(*p != '\0' && *p != ' ')
                p ++;
            while(*p != '\0' && *p == ' ')
                p ++;

            /* copies the rest of the line */
            if(*p)
                strcat(newline, p);
	}
    }
    /* we must add a new line */
    else {
        /* verifier qu'on ne depasse pas header->maxlines !!! */
        if(numlines >= maxlines) {
            maxlines += 1024;
/*protoize:???*/ /* idem*/
            if((lines = (char **)/*???*/realloc(lines, maxlines * sizeof(char *))) == NULL) {
                PUT_ERR( "[put_header] memory reallocation error(1)!\n");
                return(0);
            }
        }
        numlines ++;
    }

    /* replaces or add the line */
    lines[i] = newline;

    /* updates header structure */
    header->maxlines = maxlines;
    header->numlines = numlines;
    header->lines = lines;

    return(1);
}



/* --getheader_interfile() -------------------------------------------------*/
/*!
**   \brief  Gets information stored in the header of the interfile
**      format
**
**    Usage:    getheader_interfile(file,stringbis,type,number)
**    \param  file: char*: file name 
**    \param  stringbis : search line
**    \param  type : type of data to read
**      - ENTIER : recupere un int derriere le string 
**      - REEL   : recupere un float derriere le string 
**      - TAB_PARA_RELLE   :  reel ligne suivant pour paravision  
**      - TAB_PARA_ENTIER   :  entier ligne suivant pour paravision  
**      - TAB_PARA_STRING   :  string ligne suivant pour paravision  
**      - TAB_PARA   :  
**      - ASCII8B   : recupere un *char derriere le string 
**      - ENTIER_IPB   : recupere un int indexe par numero derriere le string 
**                     expl : width[numero]:int
**                            -----   <-- string (width)
**      -REEL_IPB  :recupere un float indexe par numero derriere le string
**
**    \param  number : dans le cas d'un vecteur???
**
**
**    \retval  This fct returns a char recupered, NULL en cas d'erreur
**                     
*/
char     *getheader_interfile(char *file, char *stringbis, int type, int number)
{
    char *fmt1,*fmt2;
    FILE *fp;
    int i,k,d;
    float f;
    static char s[4096], string[100];  

    if((fp=fopen(file,"rb"))==NULL) { 
	   /* printf("File open error: %s \n", file);*/
        return(0l);    /*Open file error... */
    }
    fmt1=CALLOC(4096,char);
    fmt2=CALLOC(4096,char);

    *s = '\0';

    if (type==ENTIER_IPB || type==REEL_IPB || (type==STRING && number)) {
        sprintf(string,"%s[%d]=",stringbis,number);
    }
    else {
        strcpy(string,stringbis);
    }
    while ((fgets(fmt1,4096,fp))!=NULL)
    {
        if (strstr(fmt1,string))
        {
            switch(type) {
            case ENTIER:
                sprintf(fmt2,"%s %%d",string);
                sscanf(fmt1,fmt2,&i);
                sprintf(s,"%d",i);
                break;
            case REEL:
                sprintf(fmt2,"%s %%f",string);
                sscanf(fmt1,fmt2,&f);
                sprintf(s,"%f",f);
                break;
            case TAB_PARA:
                sprintf(fmt2,"%s %%d",string);
                sscanf(fmt1,fmt2,&i);
                fgets(fmt1,4096,fp);
                if (number==0) {  /* 1 element de la suite */
                    strcpy(fmt2,"%d");
                    sscanf(fmt1,fmt2,&i);
                    sprintf(s,"%d",i);
                    break;
                }
                strcpy(fmt2,"%*d"); /*   sinon  */
                for (k=1;k<i;k++) {
                    if (number==k) {
                        strcat(fmt2,"%d");
                        sscanf(fmt1,fmt2,&i);
                        sprintf(s,"%d",i);
                        break;
                    }
                    else
                        strcat(fmt2,"%*d");
                }
                if(k==i) sprintf(s,"%d",0);
                break;
            case TAB_PARA_REEL:
                sprintf(fmt2,"%s %%d",string);
                sscanf(fmt1,fmt2,&i);
                fgets(fmt1,4096,fp);
                if (number==0) {  /* 1 element de la suite */
                    strcpy(fmt2,"%f");
                    sscanf(fmt1,fmt2,&f);
                    sprintf(s,"%f",f);
                    break;
                }
                strcpy(fmt2,"%*f"); /*   sinon  */
                for (k=1;k<i;k++) {
                    if (number==k) {
                        strcat(fmt2,"%f");
                        sscanf(fmt1,fmt2,&f);
                        sprintf(s,"%f",f);
                        break;
                    }
                    else
                        strcat(fmt2,"%*f");
                }
                if(k==i) sprintf(s,"%f",0.);
                break;
            case TAB_PARA_ENTIER:
                sprintf(fmt2,"%s %%d",string);
                sscanf(fmt1,fmt2,&i);
                fgets(fmt1,4096,fp);
                if (number==0) {  /* 1 element de la suite */
                    strcpy(fmt2,"%d");
                    sscanf(fmt1,fmt2,&d);
                    sprintf(s,"%d",d);
                    break;
                }
                strcpy(fmt2,"%*d"); /*   sinon  */
                for (k=1;k<i;k++) {
                    if (number==k) {
                        strcat(fmt2,"%d");
                        sscanf(fmt1,fmt2,&d);
                        sprintf(s,"%d",d);
                        break;
                    }
                    else
                        strcat(fmt2,"%*d");
                }
                if(k==i) sprintf(s,"%d",0);
                break;
             case TAB_PARA_STRING:
                fgets(fmt1,4096,fp);
				sprintf(s,"%s",fmt1);
                break;
           case ASCII8B:
                sprintf(fmt2,"%s %%s",string);
                sscanf(fmt1,fmt2,s);
                break;
            case STRING:
                sprintf(fmt2,"%s%%[^\n]",string);
                sscanf(fmt1,fmt2,s);
                break;
            case ASCII8BNL:
                if (fgets(fmt1,4096,fp) == NULL) {
				  fclose(fp);
				  return((char*)NULL);
				  }
                strcpy(s,fmt1);
                break;
            case ENTIER_IPB:  /* lit le premier entier qui suit la chaine */
                sprintf(fmt2,"%s %%d",string/*,number*/);
                sscanf(fmt1,fmt2,&i);
                sprintf(s,"%d",i);
                break;
            case ENTIER_IPB_2D:  /* lit un entier a une position donnee */
                sprintf(fmt2,"%s",string);
                for (k=1;k<number;k++) {
                    strcat(fmt2,"%*d");
                }
                strcat(fmt2,"%d");
                sscanf(fmt1,fmt2,&i);
                sprintf(s,"%d",i);
                break;
            case REEL_IPB:  /* lit le premier reel qui suit la chaine */
                sprintf(fmt2,"%s %%f",string/*,(float)number*/);
                sscanf(fmt1,fmt2,&f);
                sprintf(s,"%f",f);
                break;
            case REEL_IPB_2D:  /* lit un reel a une position donnee */
                sprintf(fmt2,"%s",string);
                for (k=1;k<number;k++) {
                    strcat(fmt2,"%*f");
                }
                strcat(fmt2,"%f");
                sscanf(fmt1,fmt2,&f);
                sprintf(s,"%f",f);
                break;
            default: 
                if(MESG_DEBUG) printf("Default case... */");
            }

            FREE(fmt1);
            FREE(fmt2);
            fclose(fp);
            return((char*)s);
        }
    }

    FREE(fmt1);
    FREE(fmt2);
    fclose(fp);
    *s=0;
    return(s);
}
//---------------------REPRIS DE visu_3d.c-----------------------------/
/* -- remplit_objet() -------------------------------------------
**      Fonction qui remplit les objets labelises pour eviter
**      d'avoir deux surfaces confondues quand on a des objets
**      contenus l'un dans l'autre.
**
**      Usage   remplit_objet(imobj, label)
**
**    imobj  : grphic3d* : l'image contenant l'objet
**    label  : int       : label de l'objet dans l'image
**
*/
void remplit_objet(grphic3d *imobj, int label_obj)
{

/*  grphic3d *imtemp, *imtemp2;
  UINT i, j, k;
  int l;

  // On isole le label de l'objet dans l'image 
  for(i = 0; i < imobj->width; i++)
    for(j = 0; j < imobj->height; j++)
      for(k = 0; k < imobj->depth; k++)
        imobj->mri[i][j][k] = (TYPEMRI)(imobj->mri[i][j][k] == label_obj);

  imobj->icomp = 0;
  imobj->rcoeff = 1;
  imobj->max_pixel = 1;
  imobj->min_pixel = 0;
  imobj->cutoff_max = 1;
  imobj->cutoff_min = 0;

  imtemp = cr_grphic3d(imobj);
  imtemp2 = cr_grphic3d(imobj);

  // On separe l'objet en composantes connexes par labelisation 
  imx_labelcroix_3d_p(imobj, imtemp);

  // RAZ_mri(imobj); 
  for(i = 0; i < imobj->width; i++)
    for(j = 0; j < imobj->height; j++)
      for(k = 0; k < imobj->depth; k++)
        imobj->mri[i][j][k] = 0;

  for (l = 0; l < imtemp->max_pixel; l++)
  {
    for(i = 0; i < imtemp->width; i++)
      for(j = 0; j < imtemp->height; j++)
        for(k = 0; k < imtemp->depth; k++)
          imtemp2->mri[i][j][k] = (TYPEMRI)(imtemp->mri[i][j][k] == l+1);

    // On comble les trous de la composante connexe 
    imx_hole_fill_3d_p(imtemp2, imtemp2);

    for(i = 0; i < imobj->width; i++)
      for(j = 0; j < imobj->height; j++)
        for(k = 0; k < imobj->depth; k++)
        {
          if (imobj->mri[i][j][k] == 0)
            if (imtemp2->mri[i][j][k]) imobj->mri[i][j][k] = label_obj;
        }
  }

  free_grphic3d(imtemp);
  free_grphic3d(imtemp2);*/

}

#endif

/* http://www.eecg.utoronto.ca/~aamodt/sourceware/MSVC.html */
double rint( double x)
// Copyright (C) 2001 Tor M. Aamodt, University of Toronto
// Permisssion to use for all purposes commercial and otherwise granted.
// THIS MATERIAL IS PROVIDED "AS IS" WITHOUT WARRANTY, OR ANY CONDITION OR
// OTHER TERM OF ANY KIND INCLUDING, WITHOUT LIMITATION, ANY WARRANTY
// OF MERCHANTABILITY, SATISFACTORY QUALITY, OR FITNESS FOR A PARTICULAR
// PURPOSE.
{
	double diff;
    if( x > 0 ) {
        __int64 xint = (__int64) (x+0.5);
        if( xint % 2 ) {
            // then we might have an even number...
            diff = x - (double)xint;
            if( diff == -0.5 )
                return (double)(xint-1);
        }
        return (double)(xint);
    } else {
        __int64 xint = (__int64) (x-0.5);
        if( xint % 2 ) {
            // then we might have an even number...
            diff = x - (double)xint;
            if( diff == 0.5 )
                return (double)(xint+1);
        }
        return (double)(xint);
    }
}

int random()
{
	static int first=1;
	if(first != 0)
	{
		srand((unsigned) time(NULL));
		first = 0;
	}

	return rand();
}

double drand48()
{
	return ((double)(random()))/(1.+(double)(RAND_MAX));
}
