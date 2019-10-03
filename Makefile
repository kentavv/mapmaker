###############################################################################
##########                                                           ##########
##########   Unix Makefile for MAPMAKER 3.0 and MAPMAKER-QTL 1.1     ##########
##########                                                           ##########
###############################################################################
# This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
# for Biomedical Research. All rights reserved. See READ.ME for license.

.SUFFIXES: .c .o .ln .lint .src

#### First, we need to set flags to define the system type. See the discussion
#### of these flags in the beginning of lib/system.h. For a Sun SPARCStation
#### running SunOS 4.1.x, use:
####
#### SYS= -D_SYS_SUNOS
####
#### For Macs with A/UX, use:	SYS= -D_SYS_AUX
#### For DECStations, use:	SYS= -D_SYS_ULTRIX
#### For PC's with WATCOM C/386 use the special Makefile instead,
#### No other ports are built in yet. Please do them, and tell us how!
#### Refer to lib/system.h and lib/syscode.c for ideas.
####
#### One other thing you may add to SYS is -D_BIG_DATASETS, allowing bigger
#### datasets (roughly 5000 loci vs 1000 loci). For even bigger datasets, 
#### modify mapm/mapm.h.

SYS= -D_SYS_SUNOS

#### Next define the system libraries to link with Mapmaker. On SunOS, Ultrix
#### and other vanilla-ish BSD systems, the correct setting is:
#### 
#### SYS_LIBS= -lm -ltermcap
####
#### On SystemV derivatives with BSD compatibility libraries, -lbsd or -lBSD
#### is probably best (for ex: on A/UX or HP/UX, respectively). 
#### On A/UX -lmalloc is likely a good idea too.

SYS_LIBS= -lm -ltermcap

#### Below we define the directories for these programs. Set DIR to the place
#### in which the executables and help files should permanently be installed
#### if you run "make install". For example:
####
#### DIR= /usr/local/bin
####
#### Be sure that you have write permission to this directory (e.g. you may
#### need to be logged in as "root" to install MAPMAKER in some directories).

DIR= /usr/local/bin

#### RDLN below refers to the GNU Readline library, which can (optionally) be
#### used by MAPMAKER to provide interactive command line editing. See the 
#### discussion of using GNU readine in the INSTALL file. If you wish to
#### use MAPMAKER with GNU Readline 1.1, change the lines below to:
#### 
#### OUR_LIB= gnu.o
#### GNU_OPT= -D_GNU_READLINE -I.
#### GNU_LIBS= -L$(RDLN) -lreadline
#### 
#### Otherwise, they should read:
#### 
#### OUR_LIB= lib.o
#### GNU_OPT=
#### 
#### Readline is not part of MAPMAKER as such: it was written by other 
#### people (to whom we owe many thanks!) and it is covered by different
#### copyright terms. See the file readline/COPYING for details. Compiling
#### readline will require modifications to readline/Makefile. Note that 
#### our version of readline has been very slightly modified from FSF V1.1
#### (see the end of readline/readline.c). Readline will not compile on
#### A/UX without using GCC instead of A/UX's cc (because it needs alloca).

OUR_LIB= gnu.o
GNU_OPT= -D_GNU_READLINE -I.
GNU_LIBS= -L$(RDLN) -lreadline

#### Now we specify the commands used to compile MAPMAKER.  The stuff
#### below works on Sun SPARCStations running SunOS 4.1.x (aka Solaris 1.x)
#### with the standard Sun C compiler installed (be careful, optimization
#### is quite broken in this compiler). With any other system or 
#### configuration, you're on your own, although a competent UNIX support
#### person at your site should be able to figure out what to do almost
#### trivially. See the file INSTALL.ME for details. 

COMPILE= cc
LINKALL= cc
LINKLIB= ld -r
DELETE=  rm -f
INSTALL= cp
SHELL= /bin/sh

MAPM= ./mapm/
QTL=  ./quant/
LIB=  ./lib/
RDLN= ./readline/
DIST= ./dist/
GS= ghostscript

#### You should be able to ignore everything below this line.
####
####
#### Define the files which make up MAPMAKER:

MAPM_OBJ = \
  $(MAPM)state.o    $(MAPM)reader.o   $(MAPM)ctm.o \
  $(MAPM)print.o    $(MAPM)sequence.o $(MAPM)main.o 	$(MAPM)quick23.o \
  $(MAPM)maps.o     $(MAPM)info.o     $(MAPM)chroms.o   $(MAPM)orders.o   \
  $(MAPM)sys_cmds.o $(MAPM)two_cmds.o $(MAPM)npt_cmds.o $(MAPM)auto_cmd.o  \
  $(MAPM)mapmhelp.o $(MAPM)database.o $(MAPM)ps_maps.o
MAPM_SRC=   $(MAPM_OBJ:.o=.c)
MAPM_LN=    $(MAPM_OBJ:.o=.ln)
MAPM_LINT=  $(MAPM_OBJ:.o=.lint)
MAPM_CLEAN= $(MAPM)*.o $(MAPM)*.ln $(MAPM)*.lint
M_HELP_IN=  $(MAPM)mapmhelp.src
M_HELP_SRC= $(MAPM)mapmhelp.c
M_INSTALL=  mapmaker mapmaker.help xmapmaker

QTL_OBJ = \
  $(QTL)qcmds.o	   $(QTL)qctm.o	     $(QTL)qdata.o      $(QTL)qseq.o \
  $(QTL)qtop.o     $(QTL)qraw.o      $(QTL)qprint.o     $(QTL)qcontext.o \
  $(QTL)qwiggle.o  $(QTL)qps_scan.o  $(QTL)qtlhelp.o
QTL_SRC=    $(QTL_OBJ:.o=.c)
QTL_LN=     $(QTL_OBJ:.o=.ln)
QTL_LINT=   $(QTL_OBJ:.o=.lint)
QTL_CLEAN=  $(QTL)*.o $(QTL)*.ln $(QTL)*.lint
Q_HELP_IN=  $(QTL)qtlhelp.src
Q_HELP_SRC= $(QTL)qtlhelp.c
Q_INSTALL=  qtl qtl.help xqtl

LIB_OBJ= \
  $(LIB)memlib.o $(LIB)mathlib.o $(LIB)iolib.o $(LIB)msglib.o $(LIB)strlib.o \
  $(LIB)shell.o  $(LIB)table.o   $(LIB)eqn.o   $(LIB)stats.o  $(LIB)syscode.o 
LIB_SRC=   $(LIB_OBJ:.o=.c)
LIB_LN=    $(LIB_OBJ:.o=.ln)
LIB_LINT=  $(LIB_OBJ:.o=.lint)
LIB_CLEAN= $(LIB)*.o $(LIB)*.ln $(LIB)*.lint $(MAKEHELP)
SYSCODE=   $(LIB)syscode
MAKEHELP=  $(LIB)makehelp

RDLN_OBJ= $(RDLN)libreadline.a
RDLN_SRC= \
  $(RDLN)readline.c $(RDLN)vi_mode.c $(RDLN)history.c $(RDLN)funmap.c \
  $(RDLN)keymaps.c  $(RDLN)emacs_keymap.c $(RDLN)vi_keymap.c \
  $(RDLN)chardefs.h $(RDLN)history.h $(RDLN)keymaps.h $(RDLN)readline.h
RDLN_CLEAN= $(RDLN)libreadline.a $(RDLN)*.o

HELP_CLEAN= $(MAPM)mapmhelp.c $(MAPM)mapmaker.doc mapmaker.help
MISC_CLEAN= lib.o gnu.o mapm.lint Temp.lint

DATA= sample.raw sample.in sample.inp sample.inq mouse.raw mouse.in mouse.prep
DOCS= READ.ME INSTALL.ME
DOCZ= FEED.ME COVER.ME CHANGE.ME
TAR_ME= $(M_INSTALL) $(Q_INSTALL) $(DATA) $(DOCZ) $(GS)
AUX_DIST= mapmaker mapmaker.help qtl qtl.help $(DATA) $(DOCZ)
TAR_FILE= $(DIST)mapm3-sparc.tar
SRC_FILE= ../mapm3-source.tar
FTP= /home/groucho/ftp/distribution/mapmaker.3b

#### The main targets generally used as arguments to make ####

all:		mapmaker qtl  #the default

mapmaker:	$(OUR_LIB) $(MAPM_OBJ)
		$(LINKALL) $(MAPM_OBJ) $(OUR_LIB) $(SYS_LIBS) -o mapmaker

qtl:		$(OUR_LIB) $(QTL_OBJ)
		$(LINKALL) $(QTL_OBJ) $(OUR_LIB) $(SYS_LIBS) -o qtl

install:	$(M_INSTALL) $(Q_INSTALL)
		$(INSTALL) $(M_INSTALL) $(DIR)
		$(INSTALL) $(Q_INSTALL) $(DIR)

help:		$(M_HELP_SRC) $(Q_HELP_SRC)

clean:		
		$(DELETE) mapmaker mapmaker.help qtl qtl.help
		$(DELETE) $(MAPM_CLEAN) $(QTL_CLEAN)  $(LIB_CLEAN)
		$(DELETE) $(RDLN_CLEAN) $(HELP_CLEAN) $(MISC_CLEAN)
		$(DELETE) *.data *.maps *.2pt *.3pt *.traits *.out *.ps x*
		$(DELETE) core *~ \#* */*~ */\#* dist/*


#### Scripts for making distributions, used as arguments to make ####

dist:		mapmaker qtl
		rm -f dist/*
		cp sun/* .
		chmod go+w $(TAR_ME) $(DOCS) $(GS)/*
		tar cvf $(TAR_FILE) $(TAR_ME)
		chmod go+w $(TAR_FILE)
		compress $(TAR_FILE)
		cp -p $(DOC0) $(DIST)
		chmod o-w $(TAR_ME) $(DOCS) $(GS)/*

aux:		mapmaker qtl
		rm -f dist/*
		cp $(AUX_DIST) aux/* dist

source:		clean
		chmod go+w * */*
		tar cvf $(SRC_FILE) *
		chmod go+w $(SRC_FILE)
		compress $(SRC_FILE)
		chmod o-w * */*

bar:		
		cd dist; bar cvf /dev/rfd0c *
		pwd
		eject


#### Make Rules ####

.c.o:
	$(COMPILE) -I$(LIB) $(SYS) -c $< -o $*.o

$(SYSCODE).o:	$(SYSCODE).c
	$(COMPILE) -I$(LIB) $(SYS) $(GNU_OPT) -c $(SYSCODE).c -o $(SYSCODE).o

$(MAKEHELP):	$(MAKEHELP).o $(OUR_LIB)
	$(LINKALL) $(MAKEHELP).o $(OUR_LIB) $(SYS_LIBS) -o $(MAKEHELP)

$(M_HELP_SRC):	$(M_HELP_IN) $(MAKEHELP)
	$(MAKEHELP) $(M_HELP_IN) $(M_HELP_SRC) mapmaker.help mapmaker.ps

$(Q_HELP_SRC):	$(Q_HELP_IN) $(MAKEHELP)
	$(MAKEHELP) $(Q_HELP_IN) $(Q_HELP_SRC) qtl.help qtl.ps

$(RDLN_OBJ):	$(RDLN_SRC)
	cd $(RDLN); make libreadline.a

lib.o:	$(LIB_OBJ)
	$(LINKLIB) $(LIB_OBJ) -o lib.o

gnu.o:	$(LIB_OBJ) $(RDLN_OBJ)
	$(LINKLIB) $(LIB_OBJ) $(GNU_LIBS) -o gnu.o



#### useful stuff for debugging ####

#### these args are for /bin/lint (not /usr/5bin/lint!) on SunOS 4.x  ####
LINT= /bin/lint -q -u -i
LINTALL= /bin/lint -q -u
#LINTPIPE= fgrep -v troublesome | fgrep -v alignment | fgrep -v ignored
LINTPIPE= fgrep -v -f Lint.grep

lint:	$(MAPM_LN) $(MAPM_LINT) $(LIB_LN) $(LIB_LINT)
	$(LINTALL) $(MAPM_LN)   $(LIB_LN) $(SYS_LIBS)   >Temp.lint
	echo "=======================================" >>Temp.lint
	fgrep -v -f Lint.grep Temp.lint $(MAPM_LINT) $(LIB_LINT) >mapm.lint
		
.c.ln:
	$(LINT) -I$(LIB) $(SYS) -c $< -o $*.ln >$*.lint

lpr:
	pr -f $(MAPM_SRC) | lp

obj:	$(LIB)system.h $(MAPM_OBJ) $(LIB_OBJ) $(RDLN_OBJ)
	#load $(MAPM_OBJ) $(LIB_OBJ) $(GNU_LIBS) $(SYS_LIBS)

src:	$(LIB)system.h $(LIB_OBJ) $(RDLN_OBJ) $(MAPM_OBJ)
	#load -I$(LIB) $(SYS)   \
  $(MAPM)state.o    $(MAPM)reader.o   $(MAPM)ctm.o \
  $(MAPM)print.o    $(MAPM)sequence.o $(MAPM)main.o \
  $(MAPM)maps.o     $(MAPM)info.c     $(MAPM)chroms.o   $(MAPM)orders.o  \
  $(MAPM)sys_cmds.o $(MAPM)two_cmds.c $(MAPM)npt_cmds.o $(MAPM)auto_cmd.o \
  $(MAPM)mapmhelp.o $(MAPM)database.o $(MAPM)ps_maps.o \
 $(LIB_OBJ) $(GNU_LIBS) $(SYS_LIBS)

unload:
	#unload $(MAPM_OBJ) $(LIB_OBJ)

