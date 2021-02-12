#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "program.h"

void to_lower_case(char* input_text)
{
    int i = 0;
    while(input_text[i] != '\0')
    {
        input_text[i] = tolower(input_text[i]);
        i++;
    }
}

void format_atom_symbol(char* input_text)
{
    if(input_text[0] != '\0')
        input_text[0] = toupper(input_text[0]);
    to_lower_case(input_text + 1);
}

void raise_error(const char* message)
{
    perror(message);
    exit(EXIT_FAILURE);
}