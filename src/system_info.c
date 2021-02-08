#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "program.h"

enum molecule_type str_to_molecule_type(const char* input_text)
{
    char* tmp = malloc(MOLECULE_NAME_LENGTH * sizeof(char));
    if(tmp == NULL) raise_error("Error with allocation temp string");
    strcpy(tmp, input_text);
    to_lower_case(tmp);

    enum molecule_type result;    

    if(strcmp(tmp, "met") == 0)
        result = met;
    else if(strcmp(tmp, "ec") == 0)
        result = ec;
    else if(strcmp(tmp, "f1ec") == 0)
        result = f1ec;
    else if(strcmp(tmp, "f2ec") == 0)
        result = f2ec;
    else if(strcmp(tmp, "monoglym") == 0)
        result = monoglym;
    else if(strcmp(tmp, "tetraglym") == 0)
        result = tetraglym;
    else if(strcmp(tmp, "peo") == 0)
        result = peo;
    else if(strcmp(tmp, "fsi") == 0)
        result = fsi;
    else if(strcmp(tmp, "tfsi") == 0)
        result = tfsi;
    else
        result = other;
    
    free(tmp);
    return result;
}

char* molecule_type_to_str(enum molecule_type molecule_type)
{
    char* molecule_name = malloc(MOLECULE_NAME_LENGTH * sizeof(char));
    if(molecule_name == NULL) raise_error("Error with allocation string with molecule name");

    switch(molecule_type) {
        case met:
            strcpy(molecule_name, "met");
            break;
        case ec:
            strcpy(molecule_name, "ec");
            break;
        case f1ec:
            strcpy(molecule_name, "f1ec");
            break;
        case f2ec:
            strcpy(molecule_name, "f2ec");
            break;
        case monoglym:
            strcpy(molecule_name, "monoglym");
            break;
        case tetraglym:
            strcpy(molecule_name, "tetraglym");
            break;
        case peo:
            strcpy(molecule_name, "peo");
            break;
        case fsi:
            strcpy(molecule_name, "fsi");
            break;
        case tfsi:
            strcpy(molecule_name, "tfsi");
            break;
        case other:
            strcpy(molecule_name, "other");
            break;
    }

    return molecule_name;
}

void parse_system_info_line(struct system_compound* compound, char* buffer)
{
    char* tmp;
    tmp = strtok(buffer, SEPARATOR);
    compound->molecule_type = str_to_molecule_type(tmp);

    tmp = strtok(NULL, SEPARATOR);
    compound->first_atom_symbol = malloc(ATOM_SYMBOL_LENGTH * sizeof(char));
    if(compound->first_atom_symbol == NULL) raise_error("Error with memory allocation for atom symbol");
    strcpy(compound->first_atom_symbol, tmp);
    format_atom_symbol(compound->first_atom_symbol);

    tmp = strtok(NULL, SEPARATOR);
    compound->quantity = atoi(tmp);

    tmp = strtok(NULL, SEPARATOR);
    compound->atoms_number = atoi(tmp);
}

void count_molecules(struct system_info* system_info)
{
    system_info->metal_ions_number = 0;
    system_info->solvent_molecules_number = 0;
    system_info->solvent_types_number = 0;
    system_info->anions_number = 0;
    system_info->atoms_number = 0;

    for(int i = 0; i < system_info->compounds_number; i++)
    {
        struct system_compound compound = system_info->compounds[i];
        enum molecule_type molecule_type = compound.molecule_type;

        if(molecule_type == met)
            system_info->metal_ions_number += compound.quantity;
        else if(molecule_type == ec || molecule_type == f1ec || molecule_type == f2ec || molecule_type == monoglym || molecule_type == tetraglym || molecule_type == peo)
        {
            system_info->solvent_molecules_number += compound.quantity;
            system_info->solvent_types_number++;
        }
        else if(molecule_type == fsi || molecule_type == tfsi)
            system_info->anions_number += compound.quantity;
        
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