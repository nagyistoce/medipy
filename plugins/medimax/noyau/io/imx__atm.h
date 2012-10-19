/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

//  ---------------------------------------------------------------------------
//
//     file:           imx__atm.h
//
//     project:        Imagix 2.01
//
//
//     description:    automatism
//
//
//     Copyright (c) 1993, ULP-IPB Strasbourg.
//     All rights are reserved.
//
// - This file specifies user menu.
// - When changing please issue your name... -
//     Last user action by: F. BOULEAU on Oct. 26th 1999
//     Last user action by: on Aug. th 1999
//
//  ---------------------------------------------------------------------------

#ifndef __IMX__ATM_H
#define __IMX__ATM_H

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#ifndef WIN32
# include <dirent.h>
#endif
#include <sys/stat.h>   
#include <sys/types.h>   
#include <math.h>

/*
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h"
#include "outils/imx_stck.h"			// Stack definition library.         
*/

#ifdef __cplusplus 
extern "C" 
{
#endif //__cplusplus 

// -- declaration des constantes et macros --------------------------------- 

#define ATM_MAX_ITEM_SIZE 20
#define ATM_MAX_FILE_PATTERN 20
#define ATM_DEFAULT_PATH MEDIMAXLIBS_PATH"/automates"

typedef struct
{
    char *begin[ATM_MAX_ITEM_SIZE];
    char *interm[ATM_MAX_ITEM_SIZE];
    char *end[ATM_MAX_ITEM_SIZE];
} TAtmBlockStruct;

typedef enum
{
    ATM_BLOCK_NONE,
    ATM_BLOCK_END,
    ATM_BLOCK_BEFORE_END,
    ATM_BLOCK_INTERMEDIATE
} TAtmBlockStatus;

typedef enum
{
    ATM_JUMP_START,
    ATM_JUMP_INTERMEDIATE,
    ATM_JUMP_END,
    ATM_JUMP_BEFORE_END
} TAtmBlockJumpCriteria;

typedef enum
{
    ATM_ACTION_NOP,
    ATM_ACTION_ERROR,
    ATM_ACTION_PARAMETER,
    ATM_ACTION_JUMPED_OUT,
    ATM_ACTION_LOOPED
} TAtmCommandReturn;

typedef enum
{
    ATM_COPY_FILE,
    ATM_COPY_DIR
} TAtmCopyType;

typedef enum
{
    ATM_PROC_AVAILABLE,
    ATM_PROC_FAILED,
    ATM_PROC_WORKING,
    ATM_PROC_DONE
} TAtmProcessingStatus;

typedef struct
{
    char    *szServerName;
    int     nPid;
    int     nStatus;
    void    *pItem;
    TAtmProcessingStatus nDone;
} TAtmProcessingServer;

typedef char *(*T_ATM_FNREAD)(LPVOID);
typedef void (*T_ATM_FNWRITE)(LPVOID, char *);

typedef struct 
{
    LPVOID          widget;
    T_ATM_FNREAD    fnread;
    T_ATM_FNWRITE   fnwrite;
} TAtmWidgetItem;

#define ATM_ENV_RESULT "IMX_ATM_RESULT"
#define ATM_SCRIPT_CPUWORK "cpuwork.sh"
#define ATM_LIST_SERVERS "servers"
#define ATM_REMOTE_RUN_SCRIPT "rimagix.sh"
#define ATM_RCP_FILE "rcp.sh"

#define ATM_SERVERNAME  "SERVERNAME"

#ifdef WIN32
#define S_ISDIR(x) (x & S_IFDIR)
#define S_ISREG(x) (x & S_IFREG)
#endif

#define ATM_FLG_ADD			1

#define ATM_TAB_ALIGN		20
#define ATM_MAX_LENLABEL	40
#define ATM_MAX_LENSTRING	512
#define ATM_MAX_LEVEL		10
#define ATM_MAX_FORIMBR		50
#define ATM_NBDIGITS		4

#define ATM_INITARGS        1
#define ATM_INITINTERF      2
#define ATM_INITALL         (ATM_INITARGS | ATM_INITINTERF)
#define ATM_CALLCMD         4
#define ATM_NOINITVARS      8

typedef enum
{
    ATM_UNINITIALIZED,
    ATM_STRVARTYPE,
    ATM_INTVARTYPE,
    ATM_FLTVARTYPE
} TAtmVarType;

typedef enum
{
    ATM_IFBLOCK,
    ATM_FORBLOCK,
    ATM_LSTBLOCK
} TAtmBlockType;

#define MAX_LENSEP          4
#define ATM_ALLSEP          "\n \t\r"
#define ATM_EXPRSEP         "+-^!=<>e()"
#define ATM_METACAR         "\"'[]"

#define ATM_ERR_WRONGTYPE        -1
#define ATM_ERR_NUMREQUIRED      -2
#define ATM_ERR_SAMETYPEREQUIRED -3
#define ATM_ERR_BOOLTYPEREQUIRED -4

#define ATM_ISSTRVAR(x)     (x->type == ATM_STRVARTYPE)
#define ATM_ISINTVAR(x)     (x->type == ATM_INTVARTYPE)
#define ATM_ISFLTVAR(x)     (x->type == ATM_FLTVARTYPE)

#define ATM_ISNUMVAR(x)     (ATM_ISINTVAR(x) || ATM_ISFLTVAR(x))
#define ATM_ISVARNAME(x)    (ATM_ISSTRVAR(x) || ATM_ISNUMVAR(x))

// returns ival for integer variables, fval for float variables and 
// evaluates string to float for others ("$1", "$2", etc.            

#define ATM_NUMVAL(x)                   \
    (ATM_ISINTVAR(x)                    \
        ? x->ival                       \
        : (ATM_ISFLTVAR(x)              \
            ? x->fval                   \
            : atof(x->str)))

// variables of the same type ? 

#define ATM_COMPTYPES(x, y) (ATM_ISINTVAR(x) && ATM_ISINTVAR(y) \
                          || ATM_ISFLTVAR(x) && ATM_ISFLTVAR(y) \
						  || ATM_ISSTRVAR(x) && ATM_ISSTRVAR(y))

// sets ival for integer variables and fval for others 

#define ATM_SETVAL(x, y)                                        \
    if((int)y == y) {                                           \
        x->type = ATM_INTVARTYPE;                               \
        x->ival = y;                                            \
    }                                                           \
    else {                                                      \
        x->type = ATM_FLTVARTYPE;                               \
        x->fval = y;                                            \
    }

#define ATM_SETSTR(x, y) {                                      \
    x->type = ATM_STRVARTYPE;                                   \
    strcpy(x->str, y); }

// -- declaration des structures ------------------------------------------- 

typedef struct s_operator
{
	char id[MAX_LENSEP];
} t_operator;

typedef struct s_variable
{
    TAtmVarType type;
	char    lbl[ATM_MAX_LENLABEL];      // variable name                     
	double  fval;                       // float variable's value            
	int     ival;                       // integer variable's value          
	char    str[ATM_MAX_LENSTRING];     // string variable's value           
} t_variable;

typedef struct s_block
{
	TAtmBlockType type;                 // block's type                      
	long    val;                        // condition's state, offset, bound  
	double  ubound;                     // upper bound (for .. next)         
	double  step;                       // step (for .. next)                
	char  lbl[ATM_MAX_LENLABEL];        // variable's name (for .. next)     
} t_block;

// parsing state

typedef enum 
{
    PARSE_STATE_NORMAL,
    PARSE_STATE_BLOCK_PARALLELIZE,
    PARSE_STATE_BLOCK_PARALLELIZE_ITEM
} TAtmParalParseStatus;

typedef enum
{
    ATM_PARSE,
    ATM_EXEC
} TAtmParalPassType;

// functions entries

typedef void (*TAtmFunction)(char *);

typedef struct 
{
    char            *szLabel;
    TAtmFunction    func;
    TAtmParalParseStatus state;
} TAtmFunctionEntry;

// nodes state

typedef enum 
{
    ATM_NODE_AVAILABLE,
    ATM_NODE_PROCESSING,
    ATM_NODE_DONE
} TAtmNodeStatus;

typedef enum
{
    ATM_NODE_DOUBLE,
    ATM_NODE_STRING
} TAtmNodeType;

typedef struct
{
    TAtmNodeStatus state;
    double        value;
} TAtmValueForItem;

typedef struct
{
    TAtmNodeStatus state;
    char          szValue[ATM_MAX_LENSTRING];
} TAtmStringForItem;

typedef union
{
    TAtmValueForItem    value_double;
    TAtmStringForItem   value_string;
} TAtmForItemNode;

typedef struct
{
    TAtmBlockType type;
    char lbl[ATM_MAX_LENLABEL];
    t_ptrlist lstValues;
} TAtmForItemContainer;

// temporary scripts

typedef struct
{
    char *szFileName;
    int  nbParameterFiles;
} TAtmTemporaryScript;

// -- declaration des variables -------------------------------------------- 

extern t_ptrlist lstargs;               // Command line paramters            
extern BOOL      _bIsParalActivated;    // Parallelization on/off?
extern int       _nAtmParalConfig;      // Configuration file # to use
extern FILE      *fpScriptParam;		// Parameters' file

// -- declaration des fonctions -------------------------------------------- 

extern void IMXAtm_Initialize();

// parallelization functions

extern BOOL IMXAtm_GenerateSubAutomates(char *szFileName);
extern char *IMXAtm_RCP(char *szFileName);

// interfaces functions

extern void IMXAtm_UnregisterWidget(LPVOID group);
extern void IMXAtm_Register(LPVOID group, LPVOID component, T_ATM_FNREAD, T_ATM_FNWRITE);

// automates creation functions

extern char *IMXAtm_ManageMenuClick(const char *);

// automates processing functions

extern char *imx_atm_evalstring(char *);
extern char *IMXAtm_AToCmd(char *, char);   // Str to a cmd.                 
extern int  IMXAtm_SendErr(const char *);
extern int	 imx_atm_decfilename(const char *, char *, char *);

extern char    *donne_liste_char(char *, int); // Extract information 
extern int     IMXAtm_ShowCmd(char *); // Manage Output 
extern char    *atocmd(char *, char); // Str to a cmd 
extern struct stack *ffopen(struct stack *, char *); // Manage files in atm files 
extern char    *IMXAtm_NextLine(int, int *); // Get next line. 
extern int	    IsCommand(char *); // Is a string is a cmd? 
extern int	    IMXAtm_EndAtm(); // When EOF reached 
extern int	    IMXAtm_StartAutomaticMode(char *, int); // Begin automatic mode. 
extern int	    IMXAtm_Exec(char *); // Call automatic mode. 
extern int	    IMXAtm_Learn(char *); // Call learning mode. 
extern int	IMXAtm_learn(char *s);
extern int     Imagix_SendEvt(char *); // Send an event. 
extern t_variable *imx_atm_searchvar(char *, int); // seach for a variable 
extern int     imx_atm_increlem(t_variable *, int, double); // increment variable 
extern TAtmBlockStatus imx_atm_jumpblock(FILE *fp, char car, char *szSearchItem, TAtmBlockJumpCriteria criteria, void (*func)(char *), BOOL bStopOnIntermediate);
extern char    *imx_atm_extrword(char *, int, char *); // extract word 
extern int     imx_atm_extrcomp(char *, t_ptrlist); // extract string items 
extern int     atm_evalatoms(t_variable *, t_variable *, t_variable *, int, double *); // evaluate minimal expression 
extern int     atm_seekoper(t_ptrlist, int); // search and evaluate minimal expr  
extern int IMXAtm_DisplayError(int err);
extern double imx_atm_evalexpr_lst(t_ptrlist, int, int, int); // evaluate expression in list 
extern double imx_atm_evalexpr(char *); // extract and evaluate expression 
extern int   imx_atm_filter_lst(t_ptrlist, char *); // apply filter on string list 
extern int   imx_atm_readdir_lst(t_ptrlist, char *, int); // read directory 
extern int imx_atm_readparamfromfile(FILE * fp,char *expr);
extern TAtmCommandReturn atm_arobasque(FILE **ptrfp, char *buff, char car);

extern char *imx_atm_evalstring(char *);
extern FILE *imx_atm_backenv(void);

extern struct stack *ffopen(struct stack *open_files_table, char *s) ;
extern t_variable *imx_atm_searchvar(char *var, int flag) ;
extern int imx_atm_increlem(t_variable *elem, int cntpos, double step) ;
extern char *imx_atm_extrword(char *expr, int n, char *sep) ;
extern int imx_atm_extrcomp(char *expr, t_ptrlist lstcomp) ;
extern int atm_evalatoms(t_variable *v1, t_variable *op, t_variable *v2, int lvl, double *pres) ;
extern int atm_seekoper(t_ptrlist lst, int lvl) ;
extern int atm_displayerror(int err) ;
extern double imx_atm_evalexpr_lst(t_ptrlist lstcomp, int deb, int fin, int lvl) ;
extern double imx_atm_evalexpr(char *expr) ;
extern int imx_atm_filter(char *str, char *patt) ;
extern int imx_atm_filter_lst(t_ptrlist lst, char *patt) ;
extern int imx_atm_readdir_lst(t_ptrlist lstdir, char *base, int type) ;
extern FILE *imx_atm_backenv(void) ;
extern char *donne_liste_char(char *buffer, int position) ;
extern int donne_liste_int(char *buffer, int position, int *erreur) ;
extern double donne_liste_float(char *buffer, int position, int *erreur) ;
extern int donne_liste_nbElement(char *buffer) ;
extern char * donne_listeElem_char(char *buffer, char *elem, int pos) ;
extern int donne_listeElem_int(char *buffer, char *elem, int pos, int *er) ;
extern double donne_listeElem_float(char *buffer, char *elem, int pos, int *er) ;
extern void Atm_Exec(void) ;
extern void Atm_Lrn(void) ;
extern void Atm_Imx(void) ;
extern void Atm_Make(void) ;
extern void Atm_Make2(void) ;
extern int imx_atm_dialogs(char *buff, char *str, char *param, int n) ;
extern int	cmdtocode(char *command, int style);
extern void IMXAtm_MainLoop(FILE *fp, int style, int *nb_lines);
extern int atm_runfile(int pass, char *file, ...) ;

#ifdef __cplusplus 
}
#endif //__cplusplus 

#endif

