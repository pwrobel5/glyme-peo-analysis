CC = gcc
CFLAGS = -Wall -std=c99
LIBS = -lm

OBJ_DIR = obj
SRC_DIR = src
TEST_DIR = test
OBJ = $(addprefix $(OBJ_DIR)/, program.o io.o reading_utils.o system_info.o vector.o coordination.o residence.o misc.o set.o solvent_data.o)
TEST_OBJ = $(addprefix $(OBJ_DIR)/, misc.o set.o)
INC = -I./include/

program: $(OBJ)
	$(CC) $(CFLAGS) $? -o $@.x $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< 

test: $(TEST_DIR)/set_test.c $(TEST_OBJ)
	$(CC) $(CFLAGS) $(INC) -o $@.x $< $(TEST_OBJ) $(LIBS) -lcunit
	./$@.x
	rm -rf $@.x

clean:
	rm -rf $(OBJ)
	rm -rf *.dat