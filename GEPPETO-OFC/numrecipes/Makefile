
LIBTARGET = math_mac
LIBEXT = a #a64 #lib
LIBTGT = $(LIBTARGET).$(LIBEXT)
HOME=/Users/tsikyrakotomalala/Documents/these/Traj_Optimization_lstm/Traj_Optimization_lstm/

#ifeq ($(LIBEXT),a64) 
#ARCH=x86_64
#else
#ARCH=i386
#endif

INCLUDEPATH=. -I$(HOME)/Includes -I$(HOME)/Includes/nrd
LIBPATH=$(HOME)/Libs
LIBOBJS=mylinpack.o #stat.o nnets.o

#CC=gcc #cl ##bcc32   ##g++
# //CC=C:\MinGW\msys\1.0\bin\gcc
# CC=C:\"Program Files"\MATLAB\R2018b\sys\lcc64\lcc64\bin\lcc64
CC=/usr/bin/g++
CFLAGS= -I$(INCLUDEPATH) -Iquote -O2 -fPIC #-arch $(ARCH)
DEBUGFLAGS= -I$(INCLUDEPATH) -Iquote -g



# %.obj: %.cpp
# 	$(CC) $(CFLAGS) -c $*.cpp
	
	
%.o: %.cpp
#.cpp.o:
	$(CC) $(CFLAGS) -c $<
	# # $(CC) $(DEBUGFLAGS) -c -Fo$@ $<
	# $(CC) $(DEBUGFLAGS) -c -dynamic -Fo$@ $<

#	$(CC) $(DEBUGFLAGS) -c $^ 
##	$(CC) $(CFLAGS) -c $^
#	$(CC) $(CFLAGS) -c $<   # borland


$(LIBTGT): $(LIBOBJS)
#	libtool -static -o $@ $^
	ar rcs $@ $^
	mv $@ $(LIBPATH)
	# lib /out:$(LIBTARGET) $(LIBOBJS)
	# copy $(LIBTARGET) $(LIBPATH)


