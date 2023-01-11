CC = gcc
CFLAGS = -Wall -std=c99
LIBS = -lm

OBJ_DIR = obj
SRC_DIR = src
TEST_DIR = test
OBJ = $(addprefix $(OBJ_DIR)/, io.o reading_utils.o system_info.o vector.o coordination.o residence.o misc.o set.o solvent_data.o)
TEST_OBJ = $(addprefix $(OBJ_DIR)/, set_test.o)
INC = -I./include/

program: $(OBJ) $(OBJ_DIR)/program.o
	$(CC) $(CFLAGS) $? -o $@.x $(LIBS)

$(OBJ_DIR)/program.o: $(SRC_DIR)/program.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< 

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< 

test: $(OBJ) $(TEST_OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(TEST_OBJ) -o $@.x $(LIBS) -lcunit
	./$@.x
	rm -rf $@.x

clean:
	rm -rf $(OBJ)
	rm -rf $(OBJ_DIR)/program.o
	rm -rf $(TEST_OBJ)
	rm -rf *.dat