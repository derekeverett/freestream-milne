#compiler
#CC = g++-7 #for mac high sierra
CC = g++ # for linux with OpenMP
#CC = pgc++ #for OpenACC GPU acceleration

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -g -Wall -fopenmp -march=native -mcmodel=large #comment or remove -fopenmp if not supported  
#CFLAGS = -acc -Minfo #comment or remove -acc if not supported
LIBS= -lgslcblas -lgsl
#build target
SRC = src
TARGET = Run

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) -o $(TARGET) $(SRC)/$(TARGET).cpp
clean:
