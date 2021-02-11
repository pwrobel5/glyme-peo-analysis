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
/*
int get_next_solvent_index(int current_index, struct system_compound* compounds, int compounds_number)
{
    current_index += 1;
    struct system_compound current_compound = compounds[current_index];
    while(current_index < compounds_number && current_compound.molecule_type != ec && current_compound.molecule_type != f1ec && current_compound.molecule_type != f2ec &&
          current_compound.molecule_type != monoglym && current_compound.molecule_type != tetraglym && current_compound.molecule_type != peo)
    {
        current_index++;
        current_compound = compounds[current_index];
    }

    return (current_index < compounds_number) ? current_index : -1;
}*/