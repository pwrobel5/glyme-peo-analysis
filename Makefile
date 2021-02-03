CC = gcc
CFLAGS = -Wall -std=c99
LIBS = -lm

OBJ_DIR = obj
SRC_DIR = src
OBJ = $(addprefix $(OBJ_DIR)/, program.o reading_utils.o system_info.o vector.o coordination.o residence.o misc.o)
INC = -I./include/

program: $(OBJ)
	$(CC) $(CFLAGS) $? -o $@.x $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< 

clean:
	rm -rf $(OBJ)
	rm -rf *.dat