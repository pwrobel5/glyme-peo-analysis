CC = gcc
CFLAGS = -Wall -std=c99
LIBS = -lm

OBJ_DIR = obj
SRC_DIR = src
TEST_DIR = test
OBJ = $(addprefix $(OBJ_DIR)/, program.o reading_utils.o system_info.o vector.o coordination.o residence.o misc.o)
INC = -I./include/

program: $(OBJ)
	$(CC) $(CFLAGS) $? -o $@.x $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $< 

test: $(TEST_DIR)/test.c
	$(CC) $(CFLAGS) -o $@.x $< -lcunit
	./$@.x
	rm -rf $@.x

clean:
	rm -rf $(OBJ)
	rm -rf *.dat