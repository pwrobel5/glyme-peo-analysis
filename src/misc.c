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

void free_system_info(struct system_info* system_info)
{
    for(int i = 0; i < system_info->compounds_number; i++)
    {
        free(system_info->compounds[i].first_atom_symbol);
    }

    free(system_info->compounds);
    free(system_info);
}

int get_next_carbonate_index(int current_index, struct system_compound* compounds, int compounds_number)
{
    current_index += 1;
    struct system_compound current_compound = compounds[current_index];
    while(current_index < compounds_number && current_compound.molecule_type != ec && current_compound.molecule_type != f1ec && current_compound.molecule_type != f2ec)
    {
        current_index++;
        current_compound = compounds[current_index];
    }

    return (current_index < compounds_number) ? current_index : -1;
}