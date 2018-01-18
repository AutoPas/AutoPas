INCLUDEFLAGS = -I../../src
CXXFLAGS =	-std=c++11 -O3 -g -Wall $(INCLUDEFLAGS) #-fmessage-length=0

#OBJS =		main.o autopas/ParticleContainer.o autopas/ParticleIterator.o autopas/Particle.o
OBJS =		main.o 

LIBS =

TARGET =	main

$(TARGET):	$(OBJS)
	$(CXX)  -o $(TARGET) $(INCLUDEFLAGS) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
