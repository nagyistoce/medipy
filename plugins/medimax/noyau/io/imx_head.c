/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		imx_header.c
***
***	project:	Imagix 1.01 
***			
***
***	description:	Fct about bruker header
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***---------------------------------------------------------------------*/
#include <config.h>
#include	<stdio.h>
#include	<string.h>
#include	<fcntl.h>
#include	<stdlib.h>
#include	<sys/types.h>

#include	"imx_head.h"

/* Functions in this module:
** ------------------------- */
int	codify_pos	(unsigned char *header, long unsigned int pos, int type, char *value);
char   *decode_pos	(unsigned char *header, long unsigned int pos, int type);

/* BRUKER Information:
What about Bruker files:

   Bruker files contains 2 parts:

      1st Part: Header: // 24BITS VALUES
         
                Values contained in header are 24bits values.
                These values may be int,char or float. These values
  may be information on bruker pictures. 
    e.g. height at position MVAR+3 (int),
         width  at position MVAR+4 (int),
         maxpixel at position MVAR+5 (int),
         minpixel at position MVAR+6 (int),
         filename for measuring at position MVAR+10 (8 sixbits c.),
  These values may be physical acquisition parameters.
    e.g. repetition in ms at position MVAR+21 (float),
         slice thickness in mm at position MVAR+41 (float).
  Or, these values may contain patient information:
    e.g. patient name and surname at position PDAT+0 (10 blocks)
 
  MVAR & PDAT represents two sections: in bruker headers.
   MVAR+3, represents MVAR offset +3 values. (3*24...)

   In bruker headers there are several sections.


   
      2nd Part: Picture data: // 32BITS VALUES
      
                Values contained in header are 32bits values.
                These values are integers. They represents pixel colors.
  from left,top to right,bottom.
  


       (left,top) 
         _________________________________________________ 
         |1st values (18432), 2nd (18433) ...            |
         |                                               |  | 
         |                                               | height
         |                                               |  | 
         |                                               |
         |_______________________________________________|
                           <- width ->          (right,bottom)
  


  18432 is header size. 24bits*768 = 8bits*4608
  Data size is width*height: (One pixel is 32bits)
     is height=128 & width=128: then size is 16384pixels -> 65536bits
    (file size: 65536+18432 -> 83968)
     is height=256 & width=256: then size is 65536pixels ->262144bits
    (file size:262144+18432 ->280576)
 
*/

/*(imx_header.c - page 3/7) */
/* -- decode_pos()------------------------------------------------------	
**	Decode one header value at a certain position and return result
**  in a string.
**
**	Usage:	decode_pos(header,pos,type)
**		header:unsigned char *: A buffer with bruker header 
**				        values.
**	        pos : unsigned long : Witch position get value ?
**		type: int : Witch value type at position pos ?
**
**	Return:	Return result in a string.
**
**  ASPECT 3000 --> Hp
**    An header contains 24bits values. These values may be int, char or
**  floats. Here we see how to interpret and decode these values.
**
**    What about integers in ASPECT ?
**        : unsigned int >= 2�0    and signed int >=2�-22
**                   and <= 2�23       signed int <=2�22
**        Here we suppose that int always contains a sign.
**        Here we transform 24bits value to a 32bits value (signed:s).
**
**           2�16    2�8     2�0             2�24    2�16    2�8     2�0       
**    |s__x___|___y___|___z___|  ==>  |s__0___|___x___|___y___|___z___| 
**     2�23    2�15    2�7             2�31    2�23    2�15    2�7  
**      
**    s2�22..................2�0 ==>  s2�30oooo2�24 2�23..................2�0 
**         An aspect integer                        An hp integer 
**
**    What about characters ?
**        In aspect you have both: ascii6bits and ascii 8bits.
**     1st:ASCII 8Bits:
**        We proceed as the same manner to get and ascii8bits as int.
**
**           2�16    2�8     2�0             2�24    2�16    2�8     2�0       
**    |___x___|___y___|___z___|  ==>  |___0___|___x___|___y___|___z___| 
**     2�23    2�15    2�7             2�31    2�23    2�15    2�7  
**      
**     2�23..................2�0 ==>   2�31oooo2�24 2�23..................2�0 
**        An aspect character                   An hp character (8b) 
**
**     2nd:ASCII 6Bits:
**        In 24bits we get 4 characters.
**
**         2�18  2�12  2�6   2�0             2�24    2�16    2�8     2�0       
**    |__w__|__x__|__y__|__z__|  ==>  |___w___|___x___|___y___|___z___| 
**     2�23  2�17  2�11  2�5           2�31    2�23    2�15    2�7  
**      
**     2�23..................2�0 ==>   2�31..........................2�0 
**       An aspect 6b character               An hp character (8b)
**
**        
**
*/
char   *decode_pos(unsigned char *header, long unsigned int pos, int type)
{ int k,max,max_bits;
  unsigned long valr1,valr2,i;
  float  f=0.0;
  static char s[32],sixbit[4],ascii[4];
  char tmp;
  unsigned long dpl,ent,dec,signe,value;

  strcpy(s,"");
  valr1=valr2=0l;
  k=2;
  for(i=pos;i<pos+3;i++) {
    valr1+=(long)(header[i])<<(8l*k);
    valr2+=(long)(header[i+3])<<(8l*k);             
    k--;
  }

  switch(type) {
    case HEXA:   sprintf(s,"%08lx",valr1);
                 break;
    case L_HEXA: sprintf(s,"%08lx%08lx",valr1,valr2);
                 break;
    case ENTIER: signe=0x00800000l&valr1;
                 value=0x007fffffl&valr1;
                 if(signe) value=0xff800000|valr1;
                 sprintf(s,"%ld",value);
              /* if(signe) value=0x007fffffl&(~valr1);  modif JP 10.93
                 sprintf(s,"%s%08ld",(signe)? "-":"+",value);*/
                 break;
    case REEL:   dpl=(long)((long)valr2&0x00ffe000l)>>13;
                 ent=((long)(1<<dpl))-1l;
                 signe=0x00800000l&valr1;
                 value=0x007fffffl&valr1;
		 /* SEE if(signe) value=0x007fffffl&(~valr1); */
                 ent=value>>(24-1-dpl);
                 dec=value&((long)(1l<<(24-1-dpl))-1);
               /* Calcule de la valeur de la partie decimale... */
               /* Division successives par deux... */
                 f=(float)ent;
                 max_bits=24-1-dpl;
               /*  la valeur de max (7) permet de determiner la precision
		  de calcul. Il faudra verifier sur d'autres images
		  la valeur de 7  */
                 max=(max_bits>=7)? 7:max_bits; 
                 for(k=0;k<=max;k++) { int bit_n;
                   bit_n=(((1l<<(max_bits-1-k))&dec)==0)? 0:1;
                   f+=(float)bit_n*(float)(1.0/(float)(1<<(k+1)));
                 }
                 if(signe!=0) f=0.0-f;
             /*    printf(":%lx %lx -> (%s)%ld.%ld:%f ... %lx:%lx#"
                         ,valr1,valr2,(signe==0)? "+":"-",ent,dec,f,ent,dec);
		*/
                 sprintf(s,"%f",f);
                 break;
    case ASC6B: /* Cas d'un caract�re sur 6bits	*/
                k=3;
                for(i=0l;i<4l;i++) {
				  tmp = (valr1&((0x0000003fl)<<(k*6)));
	              sixbit[i]=(tmp>>(k--*6))+'A'-1;
                  sixbit[i]=(sixbit[i]=='@')? ' ':sixbit[i];
                  sixbit[i]=(sixbit[i]>='p')? sixbit[i]-'p'+'0':sixbit[i];
                }
                sixbit[i]=0;
                sprintf(s,"%4s",sixbit);
                break;
    case ASCII8B: /* Cas d'un caract�re cod� en 8bits... */
                k=2;
                for(i=0l;i<3l;i++) {	
					tmp = 	(valr1&((0x000000ffl)<<(k*8)));
	          		ascii[i]=(tmp>>(k--*8))+127+1;
                }
                ascii[i]=0;
                sprintf(s,"%3s",ascii);
                break;
    default: printf("Default case... */");
  }

  return((char*)s);
}

/* -- codify_pos() -----------------------------------------------------
**
*/
int     codify_pos(unsigned char *header, long unsigned int pos, int type, char *value)
{ int signe,i,k,sixbit[8];

  switch(type) {
    case L_HEXA:
    case HEXA:   /* See Below */
    case ENTIER:{  
		   union trans
			{
			long valr;
			unsigned char buff[4];
			} tr;
		   
		   
	           
		   tr.valr=atol(value);
                   if(tr.valr<0) { 
                     signe=1;
                     /* tr.valr=0-tr.valr;  Modif JP 08/1996 */
                   } 
                   else signe=0;
                  
#if defined (LINUX) || defined(WIN32) 		   	
                   /* Pour LINUX */	
		   header[pos+1]=tr.buff[2];
		   header[pos]=tr.buff[1];
		   header[pos-1]=tr.buff[0];
                   if(signe==0)  
                        header[pos+1]=header[pos+1]&0x7f; 
                   else header[pos+3]=(0x80l)|(header[pos+3]&0xff);
#else		   
		   /* Pour HP 9000 */
		   k=2;
                   for(i=pos;i<pos+3;i++) { 
                     header[i]=((0xffl<<(8l*k))&tr.valr)>>(8*k);
                     k--;
                   }
                   if(signe==0)  
                        header[pos]=header[pos]&0x7f; 
                   else header[pos]=(0x80l)|(header[pos]&0xff);
#endif		
			
                }break; 
    case REEL:   break;
    case ASC6B: if(strlen(value)<4) return(0); 
                k=3;
                for(i=0l;i<4l;i++) {
	          sixbit[i]=value[i]-'A'+1;
                  sixbit[i]=(value[i]==' ')? 0:sixbit[i];
                  sixbit[i]=(value[i]=='@')? 0:sixbit[i];
                  sixbit[i]=(value[i]>='0'&&value[i]<='9')? 
                                 sixbit[i]+'p'-'0':sixbit[i];
                }
                header[pos+0]= ((sixbit[0]<<2l)
                             |(((sixbit[1]&0xb0)>>4l)&0x03))&0xff;
                header[pos+1]= ((sixbit[1]<<4l)
                             |(((sixbit[2])>>2l)&0xff))&0xff;
                header[pos+2]= ((sixbit[2]<<6l)
                             |  (sixbit[3]&0x3f))&0xff;
                break;
    case ASCII8B:if(strlen(value)<3) return(0); 
                 i=pos;
	         header[i+0]=value[0]-127-1;
	         header[i+1]=value[1]-127-1;
	         header[i+2]=value[2]-127-1;
                 break;
                break;
    default: /* Unknow how to proceed. */ break;
  }

  return(0);
}

