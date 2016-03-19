CC = g++

CLASSDIR = /home/cristinel/noc/reliablenoc
INCDIRS = $(CLASSDIR)/include


OBJS = app.o main.o nocmap.o arch.o sa.o mapping.o router.o param.o reliability.o
HDRS = include/app.h include/nocmap.h include/arch.h include/sa.h include/mapping.h include/msg.h include/param.h include/reliability.h
SRCS = app.cpp main.cpp nocmap.cpp arch.cpp sa.cpp mapping.cpp router.cpp param.cpp reliability.cpp


EXEC = nocmap

OPT = -O3 -fomit-frame-pointer -fexpensive-optimizations -Wall
#OPT = -ggdb -DDEBUG -Wall

FLAGS = $(OPT)
FLAGS += $(addprefix -I, $(INCDIRS))


$(EXEC): $(OBJS)
	$(CC) $(FLAGS) $(OBJS) -o $(EXEC)


app.o: $(HDRS) app.cpp
	$(CC) $(FLAGS) -c app.cpp

param.o: $(HDRS) param.cpp
	$(CC) $(FLAGS) -c param.cpp

nocmap.o: $(HDRS) nocmap.cpp
	$(CC) $(FLAGS) -c nocmap.cpp

main.o: $(HDRS) main.cpp
	$(CC) $(FLAGS) -c main.cpp

router.o: $(HDRS) router.cpp
	$(CC) $(FLAGS) -c router.cpp

arch.o: $(HDRS) arch.cpp
	$(CC) $(FLAGS) -c arch.cpp

sa.o: $(HDRS) sa.cpp
	$(CC) $(FLAGS) -c sa.cpp

mapping.o: $(HDRS) mapping.cpp
	$(CC) $(FLAGS) -c mapping.cpp

reliability.o: $(HDRS) reliability.cpp
	$(CC) $(FLAGS) -c reliability.cpp

.PHONY: clean

clean:
	rm -f $(EXEC) $(OBJS) *~
