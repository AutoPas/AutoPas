INCLUDEFLAGS = -I../../src
CXXFLAGS =	-std=c++11 -O0 -g -Wall $(INCLUDEFLAGS) #-fmessage-length=0

#OBJS =		main.o autopas/ParticleContainer.o autopas/ParticleIterator.o autopas/Particle.o
OBJS =		md-main.o ../../src/particles/MoleculeLJ.o ../../src/pairwiseFunctors/LJFunctor.o

LIBS =

TARGET =	md-main

$(TARGET):	$(OBJS)
	$(CXX)  -o $(TARGET) $(INCLUDEFLAGS) $(OBJS) $(LIBS)

all: $(OBJS) $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@
