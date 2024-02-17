# Use option "-g" to debug.
CFLAGS = 
CC=cc
RM=/bin/rm -rf
DEPS=
LOCAL_LIBRARIES= -lm
OBJS = convmain.o convi.o convo.o convconv.o conviso.o convnr.o 

.c.o:
	 $(RM) $@
	 $(CC) -c $(CFLAGS) $*.c
#####


all : conv 


conv: $(OBJS)
	 $(RM) $@
	 $(CC)  -o  $@ $(OBJS) $(LOCAL_LIBRARIES)


clean::
	 $(RM) conv $(OBJS)

