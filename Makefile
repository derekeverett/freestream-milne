#compiler
#CC = g++-7 #for mac high sierra
CC = g++ # for linux with OpenMP
#CC = pgc++ #for OpenACC GPU acceleration

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -fopenmp -march=native -mcmodel=large -static-libstdc++ #comment or remove -fopenmp if not supported  
#CFLAGS = -acc -ta=tesla:managed #comment or remove -acc if not supported
LIBS= -lgslcblas -lgsl
#build target
SRC = src
TARGET = RunFreestreamMilne

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
