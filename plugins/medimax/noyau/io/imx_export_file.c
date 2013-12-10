/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/***-----------------------------------------------------------------------
 ***
 ***	file:		manipulation.c
 ***
 ***	project:	Imagix 1.01
 ***
 ***
 ***	description:	Operation manipulation des images
 ***
 ***
 ***	Copyright (c) 1993, ULP-IPB Strasbourg.
 ***	All rights are reserved.
 ***
 ***     Last user action by: Mr. ARMSPACH on July 23th 1993
 ***
 ***---------------------------------------------------------------------*/

#include <config.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h"
#include "noyau/io/imx_export_file.h"

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
  char *s = 0x0,*s1 = 0x0,*s2  = 0x0;
  int alaligne = 1;

  s1 = CALLOC(strlen(s_src1)+1,char); //correction Berst
  strcpy(s1,s_src1);
  if (s_src2)
    {
      s2=CALLOC(strlen(s_src2)+1,char); //correction Berst
      strcpy(s2,s_src2);
      s=CALLOC(strlen(s1)+strlen(s2)+sizeof(int)+16,char);
    }
  else
    s=CALLOC(strlen(s1) + sizeof(int) + 16,char);
  
  if(s==NULL) {
    if(MESG_DEBUG) printf(" errno=%d\n",errno); }
  if (s1[strlen(s1) - 1] == '_')
    {
      alaligne = 1;
      s1[strlen(s1)-1] = 0;
    }
  if( s2 && strlen(s2) && s2[strlen(s2)-1]=='_')
    {
      s2[strlen(s2)-1] = 0;
      sprintf(s,"%s %d %s",s1,i,s2);
    }
  else
    if (s2)
      sprintf(s,"%s %d %s\n",s1,i,s2);
    else
      sprintf(s,"%s %d\n",s1,i);

  if(_dialog_messages_status & (1l<<0))
    { /* Lotus file is active. */      
      if(_fp_lotus)
	fprintf(_fp_lotus,"%d %s", i, alaligne? "\n"  : "");
    }
  else
    PUT_MESG(s);


  free(s);
  free(s1);
  if (s2)
    free(s2);
}


void spF_xprintf(char *s_src1, float f, char *s_src2)
{
  char *s = 0x0,*s1 = 0x0,*s2 = 0x0;
  int len_s1;
  int alaligne;

  alaligne=0;
  len_s1 = strlen(s_src1);
  s1 = CALLOC(len_s1, char);

  strcpy(s1,s_src1);
  if (s_src2)
    {
      s2=CALLOC(strlen(s_src2),char);
      strcpy(s2,s_src2);
      s=CALLOC(strlen(s1) + strlen(s2) + sizeof(float) + 30,char);
    }
  else
    s=CALLOC(strlen(s1) + sizeof(float) + 30,char);

  if(s==NULL) {
    if(MESG_DEBUG) printf(" errno=%d\n",errno); }
  if(s1[strlen(s1)-1] == '_')
    {
      alaligne = 1;
      s1[strlen(s1)-1] = 0;
    }
  if(s2 && strlen(s2) && s2[strlen(s2)-1] == '_')
    {
      s2[strlen(s2)-1]=0;
      sprintf(s,"%s %.15f %s",s1,f,s2);
    }
  else
    if (s2)
      sprintf(s,"%s %.15f %s\n",s1,f,s2);
    else
      sprintf(s,"%s %.15f\n",s1,f);

  if(_dialog_messages_status & (1l<<0))
    { /* Lotus file is active. */
      /*    Test d'ecriture dans un fichier lotus (alaligne=0) */
      if(_fp_lotus)
	fprintf(_fp_lotus,"%f %s", f, alaligne ? "\n" : "");
    }
  else
    PUT_MESG(s);
  free(s);
  free(s1);
  if (s2)
    free(s2);
}


void spA_xprintf(char *s_src1, char *s_in, char *s_src2)
{
  char *s,*s1,*s2;
  int alaligne;

  alaligne=0;
  s1=CALLOC(strlen(s_src1),char);
  strcpy(s1,s_src1);
  s2=CALLOC(strlen(s_src2),char);
  strcpy(s2,s_src2);

  s=CALLOC(strlen(s1)+strlen(s2)+strlen(s_in)+16,char);
  if(s==NULL) {
    if(MESG_DEBUG) printf(" errno=%d\n",errno); }
  if(s1[strlen(s1)-1] == '_')
    {
      alaligne=1;
      s1[strlen(s1)-1] = 0;
    }
  if(s2[strlen(s2)-1] == '_')
    {
      s2[strlen(s2)-1] = 0;
      sprintf(s,"%s %s %s",s1,s_in,s2);
    }
  else
    sprintf(s,"%s %s %s\n",s1,s_in,s2);

  if (_dialog_messages_status & (1l<<0)) { /* Lotus file is active. */
    /*    Test d'ecriture dans un fichier lotus (alaligne=0) */
    if (_fp_lotus)
      fprintf(_fp_lotus,"%s %s",s_in, alaligne ? "\n" : "");
  }
  else
      PUT_MESG(s);

  free(s);
  free(s1);
  free(s2);
}

/*
** Should never access this variable somwhere else.
*/
int IMXExport_ShouldPutHeader = FALSE;

int	IMXPutv_OpenLotusFile(char *filename)
{
  FILE *fp;
  char file[FILE_LEN];
  int answer,exists;
  struct stat buf;

  strncpy(file,filename,FILE_LEN);
  if(!strcmp(file,_lotus_filename)) return(2);

  errno=0;
  exists = stat(file, &buf) ;
  answer=0;
  if(exists == 0) {
    int erreur;

    if(GET_EXIST(TEXT0166, &erreur) != 1)
      {
	if(GET_EXIST(TEXT0169, &erreur) != 1)
	  {
	    /* On ne touche pas au fichier */
	    return(1);
	  }
	answer=0;
	fp = fopen(file,"w+");/* On supprime le fichier et on le recreer */
      }
    else {  /* On rajoute */
      answer=1;
      fp = fopen(file,"a+");
    }
  }
  else
    fp = fopen(file,"w+");/* Le fichier n'existe pas, on peut le creer */
  IMXExport_ShouldPutHeader = TRUE;
  _fp_lotus = fp;
  return (answer+1);
}


void	IMXPutv_CloseLotusFile(void)
{

  if (_fp_lotus == NULL) return;
  fflush(_fp_lotus);
  fclose(_fp_lotus);
  _fp_lotus=NULL;

  _dialog_messages_status &= ~(1l<<0);
  _lotus_filename[0] = 0;
  IMXExport_ShouldPutHeader = FALSE;
}


int IMXExport_PutHeader(void)
{
  return IMXExport_ShouldPutHeader;
}

void IMXExport_HeaderSet(void)
{
  IMXExport_ShouldPutHeader = FALSE;
}
