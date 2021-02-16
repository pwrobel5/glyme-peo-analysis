#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "program.h"

#define ENTRY_TYPE_LENGTH 10
#define COMPOUND_NAME_LENGTH 20
#define ATOM_SYMBOL_LENGTH 4

enum entry_type str_to_entry_type(const char* input_text)
{
    char* tmp = malloc(ENTRY_TYPE_LENGTH * sizeof(char));
    if(tmp == NULL) raise_error("Error with allocation temp string");
    strcpy(tmp, input_text);
    to_lower_case(tmp);

    enum entry_type result;    

    if(strcmp(tmp, "cation") == 0)
        result = cation;
    else if(strcmp(tmp, "solvent") == 0)
        result = solvent;
    else if(strcmp(tmp, "anion") == 0)
        result = anion;
    else
        result = other;
    
    free(tmp);
    return result;
}

char* entry_type_to_str(enum entry_type entry_type)
{
    char* entry_name = malloc(ENTRY_TYPE_LENGTH * sizeof(char));
    if(entry_name == NULL) raise_error("Error with allocation string with molecule name");

    switch(entry_type) {
        case cation:
            strcpy(entry_name, "cation");
            break;
        case solvent:
            strcpy(entry_name, "solvent");
            break;
        case anion:
            strcpy(entry_name, "anion");
            break;
        case other:
            strcpy(entry_name, "other");
            break;
    }

    return entry_name;
}

void parse_system_info_line(struct system_compound* compound, char* buffer)
{
    char* tmp;
    tmp = strtok(buffer, SEPARATOR);
    compound->entry_type = str_to_entry_type(tmp);

    tmp = strtok(NULL, SEPARATOR);
    compound->compound_name = malloc(COMPOUND_NAME_LENGTH * sizeof(char));
    if(compound->compound_name == NULL) raise_error("Error with memory allocation for compound name");
    strcpy(compound->compound_name, tmp);

    tmp = strtok(NULL, SEPARATOR);
    compound->first_atom_symbol = malloc(ATOM_SYMBOL_LENGTH * sizeof(char));
    if(compound->first_atom_symbol == NULL) raise_error("Error with memory allocation for atom symbol");
    strcpy(compound->first_atom_symbol, tmp);
    format_atom_symbol(compound->first_atom_symbol);

    tmp = strtok(NULL, SEPARATOR);
    compound->quantity = atoi(tmp);

    tmp = strtok(NULL, SEPARATOR);
    compound->atoms_number = atoi(tmp);

    compound->tracked_atom_symbol = NULL;
    compound->tracked_atoms_number = 0;

    if((tmp = strtok(NULL, SEPARATOR)) != NULL)
    {
        compound->tracked_atom_symbol = malloc(ATOM_SYMBOL_LENGTH * sizeof(char));
        if(compound->tracked_atom_symbol == NULL) raise_error("Error with memory allocation for tracked atom symbol");
        strcpy(compound->tracked_atom_symbol, tmp);
        format_atom_symbol(compound->tracked_atom_symbol);

        if((tmp = strtok(NULL, SEPARATOR)) != NULL)
            compound->tracked_atoms_number = atoi(tmp);
        else
            compound->tracked_atoms_number = 1;
    }
}

void count_molecules(struct system_info* system_info)
{
    system_info->cations_number = 0;
    system_info->cation_types_number = 0;
    system_info->solvent_molecules_number = 0;
    system_info->solvent_types_number = 0;
    system_info->anions_number = 0;
    system_info->anion_types_number = 0;
    system_info->atoms_number = 0;

    for(int i = 0; i < system_info->compounds_number; i++)
    {
        struct system_compound compound = system_info->compounds[i];
        enum entry_type entry_type = compound.entry_type;

        if(entry_type == cation)
        {
            system_info->cations_number += compound.quantity;
            system_info->cation_types_number++;
        }
        else if(entry_type == solvent)
        {
            system_info->solvent_molecules_number += compound.quantity;
            system_info->solvent_types_number++;
        }
        else if(entry_type == anion)
        {
            system_info->anions_number += compound.quantity;
            system_info->anion_types_number++;
        }
        
        system_info->atoms_number += (compound.quantity * compound.atoms_number);
    }
}

struct system_info* get_system_info(const char* system_file_name)
{
    struct system_info* result = malloc(sizeof(struct system_info));
    if(result == NULL) raise_error("Error with memory allocation for system information structure");
    
    FILE* system_file = fopen(system_file_name, "r");
    if(system_file == NULL) raise_error("Error with opening system file");

    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with memory allocation for text buffer");

    fgets(buffer, MAX_LINE_LENGTH, system_file);
    result->compounds_number = atoi(buffer);
    result->compounds = malloc(result->compounds_number * sizeof(struct system_compound));
    if(result->compounds == NULL) raise_error("Error with allocation system compounds array");

    int index = 0;
    while(fgets(buffer, MAX_LINE_LENGTH, system_file) != NULL)
    {
        parse_system_info_line(&(result->compounds[index]), buffer);
        index++;
    }

    free(buffer);
    fclose(system_file);

    count_molecules(result);

    return result;
}

void free_system_info(struct system_info* system_info)
{
    for(int i = 0; i < system_info->compounds_number; i++)
    {
        free(system_info->compounds[i].compound_name);
        free(system_info->compounds[i].first_atom_symbol);

        if(system_info->compounds[i].tracked_atom_symbol != NULL)
            free(system_info->compounds[i].tracked_atom_symbol);
    }

    free(system_info->compounds);
    free(system_info);
}