#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "program.h"

#define OUTPUT_FILE_NAME_LENGTH 60
#define TRACKED_O_ATOMS_ANION 4

void check_symbol(const char* read, const char* expected)
{
    if(strcmp(read, expected) != 0)
    {
        errno = EINVAL;
        char buffer[MAX_LINE_LENGTH];
        sprintf(buffer, "Unexpected atom symbol! Expected %s, read %s", expected, read);
        raise_error(buffer);
    }
}

void read_cation_data(FILE* input_file, struct vector** cation_positions, struct system_compound compound, int* shift)
{
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < compound.quantity; i++)
    {
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i + *shift;

        for(int j = 0; j < compound.tracked_atoms_number; j++)
            read_vector_coordinates(buffer, compound.tracked_atom_symbol, &(cation_positions[index][j]));
    }

    *shift += compound.quantity;
    free(buffer);
}

void read_ligand_molecule_data(FILE* input_file, struct vector** tracked_positions, struct system_compound compound)
{
    char* tracked_atom = compound.tracked_atom_symbol;
    char* read_symbol;
    char* tmp = malloc(MAX_LINE_LENGTH * sizeof(char));
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(tmp == NULL || buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < compound.quantity; i++)
    {
        int atoms_left = compound.atoms_number;
        int tracked_index = 0;

        fgets(buffer, MAX_LINE_LENGTH, input_file);
        strcpy(tmp, buffer);
        read_symbol = strtok(tmp, SEPARATOR);
        check_symbol(read_symbol, compound.first_atom_symbol);
        if(strcmp(read_symbol, tracked_atom) == 0)
        {
            read_vector_coordinates(buffer, tracked_atom, tracked_positions[i] + tracked_index);
            tracked_index++;

            // omit drude particle
            fgets(buffer, MAX_LINE_LENGTH, input_file);
            atoms_left -= 2;
        }
        else
            atoms_left--;

        while(atoms_left > 0)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
            strcpy(tmp, buffer);
            read_symbol = strtok(tmp, SEPARATOR);
            if(strcmp(read_symbol, tracked_atom) == 0)
            {
                read_vector_coordinates(buffer, tracked_atom, tracked_positions[i] + tracked_index);
                tracked_index++;

                // omit drude particle
                fgets(buffer, MAX_LINE_LENGTH, input_file);
                atoms_left -= 2;
            }
            else
                atoms_left--;
        }
    }

    free(buffer);
    free(tmp);
}

void read_anion_data(FILE* input_file, struct vector** tracked_positions, struct system_compound anion, int* shift)
{
    char* tracked_atom = anion.tracked_atom_symbol;
    char* read_symbol;
    char* tmp = malloc(MAX_LINE_LENGTH * sizeof(char));
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(tmp == NULL || buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < anion.quantity; i++)
    {
        int anion_index = i + *shift;
        int atoms_left = anion.atoms_number;
        int tracked_index = 0;

        fgets(buffer, MAX_LINE_LENGTH, input_file);
        strcpy(tmp, buffer);
        read_symbol = strtok(tmp, SEPARATOR);
        check_symbol(read_symbol, anion.first_atom_symbol);
        if(strcmp(read_symbol, tracked_atom) == 0)
        {
            read_vector_coordinates(buffer, tracked_atom, tracked_positions[anion_index] + tracked_index);            
            tracked_index++;

            // omit drude particle
            fgets(buffer, MAX_LINE_LENGTH, input_file);
            atoms_left -= 2;
        }
        else
            atoms_left--;

        while(atoms_left > 0)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
            strcpy(tmp, buffer);
            read_symbol = strtok(tmp, SEPARATOR);
            if(strcmp(read_symbol, tracked_atom) == 0)
            {
                read_vector_coordinates(buffer, tracked_atom, tracked_positions[anion_index] + tracked_index);
                tracked_index++;

                // omit drude particle
                fgets(buffer, MAX_LINE_LENGTH, input_file);
                atoms_left -= 2;
            }
            else
                atoms_left--;
        }
    }

    *shift += anion.quantity;
    free(buffer);
    free(tmp);
}

void omit_other_data(FILE* input_file, struct system_compound compound)
{
    int atoms_number = compound.quantity * compound.atoms_number;
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < atoms_number; i++)
    {
        fgets(buffer, MAX_LINE_LENGTH, input_file);
    }

    free(buffer);
}

FILE** open_cation_files(struct system_info* system_info)
{
    FILE** result = malloc(system_info->cations_number * sizeof(FILE*));
    if(result == NULL) raise_error("Error with array for cation output files allocation");

    char filename_buffer[MAX_LINE_LENGTH];

    int shift = 0;
    for(int i = 0; i < system_info->compounds_number; i++)
    {
        struct system_compound compound = system_info->compounds[i];
        if(compound.entry_type == cation)
        {
            for(int j = 0; j < compound.quantity; j++)
            {
                sprintf(filename_buffer, "%s_%03d.dat", compound.first_atom_symbol, j);

                result[shift] = fopen(filename_buffer, "w");
                if(result[shift] == NULL) raise_error("Error with creating output file for metal ions");

                shift++;
            }
        }
    }

    return result;
}

int get_next_entry_index(int current_index, struct system_compound* compounds, int compounds_number, enum entry_type entry_type)
{
    current_index += 1;
    struct system_compound current_compound = compounds[current_index];
    while(current_index < compounds_number && current_compound.entry_type != entry_type)
    {
        current_index++;
        current_compound = compounds[current_index];
    }

    return (current_index < compounds_number) ? current_index : -1;
}

void initialize_entry_data(struct entry_data entry_data, struct system_info* system_info, enum entry_type entry_type)
{
    int*** last_coordination = entry_data.last_coordination;
    int*** current_coordination = entry_data.current_coordination;    
    struct ligand_coord_info*** coordination_info = entry_data.coordination_info;

    if(last_coordination == NULL || current_coordination == NULL || coordination_info == NULL) 
        raise_error("Error with memory allocation for coordination history");

    int current_entry_index = -1;
    int array_index_shift = 0;
    int types_number = (entry_type == solvent) ? system_info->solvent_types_number : system_info->anion_types_number;

    for(int i = 0; i < types_number; i++)
    {
        current_entry_index = get_next_entry_index(current_entry_index, system_info->compounds, system_info->compounds_number, entry_type);
        struct system_compound current_compound = system_info->compounds[current_entry_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            int current_molecule_index = array_index_shift + j;
            last_coordination[current_molecule_index] = malloc(current_compound.tracked_atoms_number * sizeof(int*));
            current_coordination[current_molecule_index] = malloc(current_compound.tracked_atoms_number * sizeof(int*));
            coordination_info[current_molecule_index] = malloc(current_compound.tracked_atoms_number * sizeof(struct ligand_coord_info*));

            if(last_coordination[current_molecule_index] == NULL || current_coordination[current_molecule_index] == NULL || coordination_info[current_molecule_index] == NULL)
                raise_error("Error with memory allocation for coordination history");
            
            for(int k = 0; k < current_compound.tracked_atoms_number; k++)
            {
                last_coordination[current_molecule_index][k] = malloc(MAX_COORDINATED_CATIONS * sizeof(int));
                current_coordination[current_molecule_index][k] = malloc(MAX_COORDINATED_CATIONS * sizeof(int));
                coordination_info[current_molecule_index][k] = malloc(MAX_COORDINATED_CATIONS * sizeof(struct ligand_coord_info));

                if(last_coordination[current_molecule_index][k] == NULL || current_coordination[current_molecule_index][k] == NULL || coordination_info[current_molecule_index][k] == NULL)
                    raise_error("Error with memory allocation for coordination history");
                
                for(int l = 0; l < MAX_COORDINATED_CATIONS; l++)
                {
                    last_coordination[current_molecule_index][k][l] = BLANK;
                    current_coordination[current_molecule_index][k][l] = BLANK;
                    coordination_info[current_molecule_index][k][l].cation_index = BLANK;
                    coordination_info[current_molecule_index][k][l].time = 0;
                }
            }
        }

        array_index_shift += current_compound.quantity;
    }
}

void free_coordination_arrays(struct entry_data entry_data, struct system_info* system_info, enum entry_type entry_type)
{
    int*** last_coordination = entry_data.last_coordination;
    int*** current_coordination = entry_data.current_coordination;    
    struct ligand_coord_info*** coordination_info = entry_data.coordination_info;

    int current_entry_index = -1;
    int array_index_shift = 0;
    int types_number = (entry_type == solvent) ? system_info->solvent_types_number : system_info->anion_types_number;

    for(int i = 0; i < types_number; i++)
    {
        current_entry_index = get_next_entry_index(current_entry_index, system_info->compounds, system_info->compounds_number, entry_type);
        struct system_compound current_compound = system_info->compounds[current_entry_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            int molecule_index = array_index_shift + j;

            for(int k = 0; k < current_compound.tracked_atoms_number; k++)
            {
                free(last_coordination[molecule_index][k]);
                free(current_coordination[molecule_index][k]);
                free(coordination_info[molecule_index][k]);
            }

            free(last_coordination[molecule_index]);
            free(current_coordination[molecule_index]);
            free(coordination_info[molecule_index]);
        }

        array_index_shift += current_compound.quantity;
    }

    free(last_coordination);
    free(current_coordination);
    free(coordination_info);
}

void open_solvent_data_files(FILE*** files, struct system_info* system_info)
{
    int compound_index = -1;
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
        struct system_compound current_compound = system_info->compounds[compound_index];

        files[i] = malloc(current_compound.quantity * sizeof(FILE*));
        if(files[i] == NULL) raise_error("Error with memory allocation for solvent data files pointers");

        for(int j = 0; j < current_compound.quantity; j++)
        {
            char file_name[OUTPUT_FILE_NAME_LENGTH];
            sprintf(file_name, "%s-%03d.dat", current_compound.compound_name, j);
            files[i][j] = fopen(file_name, "w");
            if(files[i][j] == NULL) raise_error("Error with opening solvent data output file");
        }
    }
}

void close_solvent_data_files(FILE*** files, struct system_info* system_info)
{
    int compound_index = -1;
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
        struct system_compound current_compound = system_info->compounds[compound_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            fclose(files[i][j]);
        }

        free(files[i]);
    }
}

void clean_current_coordinating_solvents_array(int* array, int array_size)
{
    for(int i = 0; i < array_size; i++)
        array[i] = BLANK;
}

void group_non_blanks_in_beginning(int* array, int array_size)
{
    int first_free_position = 0;
    for(int i = 0; i < array_size; i++)
    {
        if(array[i] != BLANK)
        {
            array[first_free_position++] = i;
            array[i] = BLANK;
        }
    }
}

void read_data(struct program_configuration* program_configuration, struct system_info* system_info)
{
    FILE* input_file = fopen(program_configuration->input_file_name, "r");

    if(input_file == NULL) raise_error("Error with opening input file!");
    
    struct vector** cation_tracked_positions = malloc(system_info->cations_number * sizeof(struct vector*));
    struct vector*** solvent_tracked_positions = malloc(system_info->solvent_types_number * sizeof(struct vector**));
    struct vector** anion_tracked_positions = malloc(system_info->anions_number * sizeof(struct vector*));
    if(cation_tracked_positions == NULL || solvent_tracked_positions == NULL || anion_tracked_positions == NULL)
        raise_error("Error with allocation of arrays for atoms positions");

    struct index_combinations** solvent_index_combinations = NULL;
    struct venn_set** venn_set_cations = NULL;
    struct venn_set* venn_set_solvent = NULL;
    if(program_configuration->calculate_venn_diagrams == 1)
    {
        solvent_index_combinations = malloc(system_info->solvent_types_number * sizeof(struct index_combinations*));
        venn_set_cations = malloc(system_info->solvent_types_number * sizeof(struct venn_set*));
        venn_set_solvent = malloc(system_info->solvent_types_number * sizeof(struct venn_set));
        if(solvent_index_combinations == NULL || venn_set_cations == NULL) raise_error("Error with allocation of data structures for Venn diagrams");
    }
    
    int cation_index = 0;
    int solvent_type_index = 0;
    int anion_index = 0;

    for(int i = 0; i < system_info->compounds_number; i++)
    {
        struct system_compound current_compound = system_info->compounds[i];

        switch(current_compound.entry_type)
        {
            case cation:
                for(int j = 0; j < current_compound.quantity; j++)
                {
                    int current_index = j + cation_index;
                    cation_tracked_positions[current_index] = malloc(current_compound.tracked_atoms_number * sizeof(struct vector));
                    if(cation_tracked_positions[current_index] == NULL) raise_error("Error with memory allocation for cation tracked positions");
                }

                cation_index += current_compound.quantity;
                break;
            
            case solvent:
                solvent_tracked_positions[solvent_type_index] = malloc(current_compound.quantity * sizeof(struct vector*));
                if(solvent_tracked_positions[solvent_type_index] == NULL) raise_error("Error with memory allocation for solvent tracked positions");
                for(int j = 0; j < current_compound.quantity; j++)
                {
                    solvent_tracked_positions[solvent_type_index][j] = malloc(current_compound.tracked_atoms_number * sizeof(struct vector));
                    if(solvent_tracked_positions[solvent_type_index][j] == NULL) raise_error("Error with memory allocation for solvent tracked positions");
                }
                if(program_configuration->calculate_venn_diagrams == 1) {
                    solvent_index_combinations[solvent_type_index] = get_index_combinations(current_compound.tracked_atoms_number);
                    venn_set_cations[solvent_type_index] = malloc(current_compound.quantity * sizeof(struct venn_set));
                    if(venn_set_cations[solvent_type_index] == NULL) raise_error("Error with memory allocation for Venn set");
                }
                solvent_type_index++;
                break;
            
            case anion:
                for(int j = 0; j < current_compound.quantity; j++)
                {
                    int current_index = j + anion_index;
                    anion_tracked_positions[current_index] = malloc(current_compound.tracked_atoms_number * sizeof(struct vector));
                    if(anion_tracked_positions[current_index] == NULL) raise_error("Error with memory allocation for anion tracked positions");
                }
                
                anion_index += current_compound.quantity;
                break;
            
            case other:
                break;
        }
    }

    FILE** cation_output_files;
    FILE* cation_output_file;
    FILE*** solvent_output_files;
    FILE* solvent_output_file;
    FILE* solvent_times_file = fopen("solvent_times.dat", "w");
    FILE* anion_times_file = fopen("anion_times.dat", "w");
    if(solvent_times_file == NULL || anion_times_file == NULL) raise_error("Error with opening output files");

    if(program_configuration->print_mode == separate)
        cation_output_files = open_cation_files(system_info);
    else
        cation_output_file = fopen("cation_output.dat", "w");
    
    if(program_configuration->save_additional_solvent_data == 1)
    {
        if(program_configuration->print_mode == separate)
        {
            solvent_output_files = malloc(system_info->solvent_types_number * sizeof(FILE**));
            if(solvent_output_files == NULL) raise_error("Error with memory allocation for solvent output files pointers");
            open_solvent_data_files(solvent_output_files, system_info);
        }
        else
            solvent_output_file = fopen("solvent_output.dat", "w");
    }

    int*** last_solvent_coordination = malloc(system_info->solvent_molecules_number * sizeof(int**));
    int*** current_solvent_coordination = malloc(system_info->solvent_molecules_number * sizeof(int**));
    struct ligand_coord_info*** solvent_coordination_info = malloc(system_info->solvent_molecules_number * sizeof(struct ligand_coord_info**));
    
    struct entry_data solvent_data;
    solvent_data.output = solvent_times_file;
    solvent_data.last_coordination = last_solvent_coordination;
    solvent_data.current_coordination = current_solvent_coordination;
    solvent_data.coordination_info = solvent_coordination_info;
    initialize_entry_data(solvent_data, system_info, solvent);

    int*** last_anion_coordination = malloc(system_info->anions_number * sizeof(int**));
    int*** current_anion_coordination = malloc(system_info->anions_number * sizeof(int**));
    struct ligand_coord_info*** anion_coordination_info = malloc(system_info->anions_number * sizeof(struct oxygen_coord_info*));

    struct entry_data anion_data;
    anion_data.output = anion_times_file;
    anion_data.last_coordination = last_anion_coordination;
    anion_data.current_coordination = current_anion_coordination;
    anion_data.coordination_info = anion_coordination_info;
    initialize_entry_data(anion_data, system_info, anion);
    
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL)
        raise_error("Error with memory allocation for reading buffer");
    
    // arrays for saving numbers of coordinated molecules to given metal ion at a given timestep
    short int**** anion_atoms_coordination_history = malloc(sizeof(short int***) * system_info->anion_types_number);
    short int**** anion_coordination_history = malloc(sizeof(short int***) * system_info->anion_types_number);
    short int**** solvent_atoms_coordination_history = malloc(sizeof(short int***) * system_info->solvent_types_number);
    short int**** solvent_coordination_history = malloc(sizeof(short int***) * system_info->solvent_types_number);
    if(solvent_coordination_history == NULL || solvent_atoms_coordination_history == NULL || anion_coordination_history == NULL || anion_atoms_coordination_history == NULL) 
        raise_error("Error with memory allocation for coordination histories array");

    struct venn_diagram*** venn_diagrams_cations = malloc(system_info->solvent_types_number * sizeof(struct venn_diagram**));
    struct venn_diagram*** venn_diagrams_solvent = malloc(system_info->solvent_types_number * sizeof(struct venn_diagram**));
    if(venn_diagrams_cations == NULL || venn_diagrams_solvent == NULL) raise_error("Error with memory allocation for Venn diagrams array");
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        solvent_coordination_history[i] = NULL;
        solvent_atoms_coordination_history[i] = NULL;
        venn_diagrams_cations[i] = NULL;
        venn_diagrams_solvent[i] = NULL;
    }
    for(int i = 0; i < system_info->anion_types_number; i++)
    {
        anion_coordination_history[i] = NULL;
        anion_atoms_coordination_history[i] = NULL;
    }

    int*** current_coordinating_solvents = malloc(system_info->solvent_types_number * sizeof(int**));
    if(current_coordinating_solvents == NULL) raise_error("Error with memory allocation for coordinating solvents array");
    int current_index = -1;
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        current_index = get_next_entry_index(current_index, system_info->compounds, system_info->compounds_number, solvent);
        struct system_compound current_compound = system_info->compounds[current_index];

        current_coordinating_solvents[i] = malloc(current_compound.tracked_atoms_number * sizeof(int*));
        if(current_coordinating_solvents[i] == NULL) raise_error("Error with memory allocation for coordinating solvents array");
        for(int j = 0; j < current_compound.tracked_atoms_number; j++)
        {
            current_coordinating_solvents[i][j] = malloc(current_compound.quantity * sizeof(int));
            if(current_coordinating_solvents[i][j] == NULL) raise_error("Error with memory allocation for coordinating solvents array");
            clean_current_coordinating_solvents_array(current_coordinating_solvents[i][j], current_compound.quantity);
        }        
    }

    int step_number = 0;

    while(fgets(buffer, MAX_LINE_LENGTH, input_file) != NULL)
    {
        int cation_shift = 0;
        int anion_shift = 0;

        int atoms_number = atoi(buffer);
        if(atoms_number != system_info->atoms_number)
            printf("[WARNING] Declared atoms number different from number read from .xyz file, declared: %d, read: %d\n", system_info->atoms_number, atoms_number);
        
        if(program_configuration->calculate_solvent_residence == 1)
        {    
            for(int i = 0; i < system_info->solvent_types_number; i++)
            {
                solvent_atoms_coordination_history[i] = initialize_history_array(solvent_atoms_coordination_history[i], step_number, system_info);
                solvent_coordination_history[i] = initialize_history_array(solvent_coordination_history[i], step_number, system_info);
            }
        }
        
        if(program_configuration->calculate_anion_residence == 1)
        {
            for(int i = 0; i < system_info->anion_types_number; i++)
            {
                anion_atoms_coordination_history[i] = initialize_history_array(anion_atoms_coordination_history[i], step_number, system_info);
                anion_coordination_history[i] = initialize_history_array(anion_coordination_history[i], step_number, system_info);
            }
        }

        if(program_configuration->calculate_venn_diagrams == 1)
        {
            int compound_index = -1;

            for(int i = 0; i < system_info->solvent_types_number; i++)
            {
                compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
                struct system_compound current_compound = system_info->compounds[compound_index];

                venn_diagrams_cations[i] = realloc(venn_diagrams_cations[i], (step_number + 1) * sizeof(struct venn_diagram*));
                venn_diagrams_solvent[i] = realloc(venn_diagrams_solvent[i], (step_number + 1) * sizeof(struct venn_diagram*));

                if(venn_diagrams_cations[i] == NULL || venn_diagrams_solvent[i] == NULL) raise_error("Error with memory allocation for Venn diagrams");

                venn_diagrams_cations[i][step_number] = create_empty_venn_diagram(current_compound.tracked_atoms_number, solvent_index_combinations[i]);
                venn_diagrams_solvent[i][step_number] = create_empty_venn_diagram(current_compound.tracked_atoms_number, solvent_index_combinations[i]);
            }
        }

        // omit commentary line
        fgets(buffer, MAX_LINE_LENGTH, input_file);

        cation_index = 0;
        solvent_type_index = 0;
        anion_index = 0;
        for(int i = 0; i < system_info->compounds_number; i++)
        {
            struct system_compound compound = system_info->compounds[i];

            switch(compound.entry_type)
            {
                case cation:
                    read_cation_data(input_file, cation_tracked_positions, compound, &cation_shift);
                    break;
                case solvent:
                    read_ligand_molecule_data(input_file, solvent_tracked_positions[solvent_type_index], compound);
                    solvent_type_index++;
                    break;
                case anion:
                    read_anion_data(input_file, anion_tracked_positions, compound, &anion_shift);
                    break;
                case other:
                    omit_other_data(input_file, compound);
                    break;
            }
        }

        struct coordination_input coordination_input;
        coordination_input.step_number = step_number;
        coordination_input.solvent_tracked_atoms = solvent_tracked_positions;
        coordination_input.anion_tracked_atoms = anion_tracked_positions;
        coordination_input.current_solvent_coordination = solvent_data.current_coordination;
        coordination_input.current_anion_coordination = anion_data.current_coordination;
        coordination_input.current_coordinating_solvents = current_coordinating_solvents;
        coordination_input.solvent_atoms_coordination_history = solvent_atoms_coordination_history;
        coordination_input.solvent_coordination_history = solvent_coordination_history;
        coordination_input.anion_atoms_coordination_history = anion_atoms_coordination_history;
        coordination_input.anion_coordination_history = anion_coordination_history;
        coordination_input.program_configuration = program_configuration;
        coordination_input.system_info = system_info;
        
        cation_shift = 0;
        int compound_index = -1;
        for(int i = 0; i < system_info->cation_types_number; i++)
        {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, cation);
            struct system_compound current_cation = system_info->compounds[compound_index];

            for(int j = 0; j < current_cation.quantity; j++)
            {
                int index = j + cation_shift;
                struct cation_coord_info coord_info = get_coordination_info(index, current_cation.tracked_atoms_number, cation_tracked_positions[index], &coordination_input);

                if(program_configuration->print_mode == separate)
                    fprintf(cation_output_files[index], "%d %d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.solvent_atoms, coord_info.solvent_molecules, coord_info.anion_atoms, coord_info.anion_molecules);
                else
                    fprintf(cation_output_file, "%d %d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.solvent_atoms, coord_info.solvent_molecules, coord_info.anion_atoms, coord_info.anion_molecules);
            }

            cation_shift += current_cation.quantity;
        }

        calculate_coord_times(system_info, &solvent_data, solvent);
        calculate_coord_times(system_info, &anion_data, anion);
        
        int solvent_shift = 0;
        compound_index = -1;
        for(int i = 0; i < system_info->solvent_types_number; i++)
        {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
            struct system_compound current_solvent = system_info->compounds[compound_index];

            if(program_configuration->calculate_venn_diagrams == 1)
            {
                for(int j = 0; j < current_solvent.tracked_atoms_number; j++)
                {
                    group_non_blanks_in_beginning(current_coordinating_solvents[i][j], current_solvent.quantity);
                }

                determine_venn_sets(&(venn_set_solvent[i]), current_coordinating_solvents[i], current_solvent.tracked_atoms_number, solvent_index_combinations[i], current_solvent.quantity);
                struct venn_diagram* venn_diagram = determine_venn_diagram(&(venn_set_solvent[i]), current_solvent.tracked_atoms_number);
                update_global_venn_diagram(venn_diagrams_solvent[i][step_number], venn_diagram);
                free_venn_diagram(venn_diagram);
                free_venn_set(&(venn_set_solvent[i]));

                for(int j = 0; j < current_solvent.tracked_atoms_number; j++)
                {
                    clean_current_coordinating_solvents_array(current_coordinating_solvents[i][j], current_solvent.quantity);
                }
            }

            for(int j = 0; j < current_solvent.quantity; j++)
            {
                int solvent_index = j + solvent_shift;

                if(program_configuration->calculate_venn_diagrams == 1) 
                {
                    determine_venn_sets(&(venn_set_cations[i][j]), solvent_data.current_coordination[solvent_index], current_solvent.tracked_atoms_number, solvent_index_combinations[i], MAX_COORDINATED_CATIONS);
                    struct venn_diagram* venn_diagram = determine_venn_diagram(&(venn_set_cations[i][j]), current_solvent.tracked_atoms_number);
                    update_global_venn_diagram(venn_diagrams_cations[i][step_number], venn_diagram);
                    free_venn_diagram(venn_diagram);
                    free_venn_set(&(venn_set_cations[i][j]));
                }

                if(program_configuration->save_additional_solvent_data == 1)
                {
                    if(program_configuration->print_mode == separate)
                        save_current_step_solvent_data(step_number, solvent_data.current_coordination[solvent_index], current_solvent.tracked_atoms_number, system_info->cations_number, solvent_output_files[i][j]);
                    else
                        save_current_step_solvent_data(step_number, solvent_data.current_coordination[solvent_index], current_solvent.tracked_atoms_number, system_info->cations_number, solvent_output_file);
                }

                swap_coordination_arrays(&(solvent_data.last_coordination[solvent_index]), &(solvent_data.current_coordination[solvent_index]));
                clear_coordination_array(solvent_data.current_coordination[solvent_index], current_solvent.tracked_atoms_number, MAX_COORDINATED_CATIONS);
            }

            solvent_shift += current_solvent.quantity;
        }

        anion_shift = 0;
        compound_index = -1;
        for(int i = 0; i < system_info->anion_types_number; i++)
        {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, anion);
            struct system_compound current_anion = system_info->compounds[compound_index];

            for(int j = 0; j < current_anion.quantity; j++)
            {
                int anion_index = j + anion_shift;
                swap_coordination_arrays(&(anion_data.last_coordination[anion_index]), &(anion_data.current_coordination[anion_index]));
                clear_coordination_array(anion_data.current_coordination[anion_index], current_anion.tracked_atoms_number, MAX_COORDINATED_CATIONS);
            }

            anion_shift += current_anion.quantity;
        }

        step_number++;
    }

    save_last_step_data(system_info, &solvent_data, solvent);
    save_last_step_data(system_info, &anion_data, anion);

    // calculate correlated residence times
    if(program_configuration->calculate_solvent_residence == 1)
    {
        int compound_index = -1;
        for(int i = 0; i < system_info->solvent_types_number; i++) {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
            struct system_compound current_compound = system_info->compounds[compound_index];
            
            int molecules_denominator = system_info->cations_number * system_info->compounds[compound_index].quantity;
            double* molecules_residence = calculate_residence_times(solvent_coordination_history[i], step_number, molecules_denominator, system_info);
            delete_history_array(solvent_coordination_history[i], step_number, system_info->cations_number);

            char residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "residence-times-";
            strcat(residence_output_name, current_compound.compound_name);
            strcat(residence_output_name, ".dat");
            save_residence_to_file(molecules_residence, residence_output_name, step_number);
            free(molecules_residence);

            int atoms_denominator = molecules_denominator * current_compound.tracked_atoms_number;
            double* atoms_residence = calculate_residence_times(solvent_atoms_coordination_history[i], step_number, atoms_denominator, system_info);
            delete_history_array(solvent_atoms_coordination_history[i], step_number, system_info->cations_number);

            char atoms_residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "residence-times-atoms-";
            strcat(atoms_residence_output_name, current_compound.compound_name);
            strcat(atoms_residence_output_name, ".dat");
            save_residence_to_file(atoms_residence, atoms_residence_output_name, step_number);
            free(atoms_residence);
        }
    }
    if(program_configuration->calculate_anion_residence == 1)
    {
        int compound_index = -1;
        for(int i = 0; i < system_info->anion_types_number; i++) {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, anion);
            struct system_compound current_compound = system_info->compounds[compound_index];

            int molecules_denominator = system_info->cations_number * current_compound.quantity;
            double* molecules_residence = calculate_residence_times(anion_coordination_history[i], step_number, molecules_denominator, system_info);
            delete_history_array(anion_coordination_history[i], step_number, system_info->cations_number);
            
            char residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "residence-times-";
            strcat(residence_output_name, current_compound.compound_name);
            strcat(residence_output_name, ".dat");
            save_residence_to_file(molecules_residence, residence_output_name, step_number);
            free(molecules_residence);
            
            int atoms_denominator = molecules_denominator * current_compound.tracked_atoms_number;
            double* atoms_residence = calculate_residence_times(anion_atoms_coordination_history[i], step_number, atoms_denominator, system_info);
            delete_history_array(anion_atoms_coordination_history[i], step_number, system_info->cations_number);

            char atoms_residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "residence-times-atoms-";
            strcat(atoms_residence_output_name, current_compound.compound_name);
            strcat(atoms_residence_output_name, ".dat");
            save_residence_to_file(atoms_residence, atoms_residence_output_name, step_number);
            free(atoms_residence);
        }              
    }

    if(program_configuration->calculate_venn_diagrams == 1)
    {
        int compound_index = -1;
        for(int i = 0; i < system_info->solvent_types_number; i++)
        {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
            struct system_compound current_compound = system_info->compounds[compound_index];

            char venn_cations_output_name[OUTPUT_FILE_NAME_LENGTH] = "venn-cations-";
            strcat(venn_cations_output_name, current_compound.compound_name);
            strcat(venn_cations_output_name, ".dat");
            save_averages_to_file(venn_diagrams_cations[i], step_number, venn_cations_output_name);

            for(int j = 0; j < step_number; j++)
            {
                free_venn_diagram(venn_diagrams_cations[i][j]);
            }
            free(venn_diagrams_cations[i]);

            char venn_solvent_output_name[OUTPUT_FILE_NAME_LENGTH] = "venn-solvent-";
            strcat(venn_solvent_output_name, current_compound.compound_name);
            strcat(venn_solvent_output_name, ".dat");
            save_averages_to_file(venn_diagrams_solvent[i], step_number, venn_solvent_output_name);

            for(int j = 0; j < step_number; j++)
            {
                free_venn_diagram(venn_diagrams_solvent[i][j]);
            }
            free(venn_diagrams_solvent[i]);
        }
    }

    if(program_configuration->save_additional_solvent_data == 1)
    {
        if(program_configuration->print_mode == separate)
        {
            close_solvent_data_files(solvent_output_files, system_info);
            free(solvent_output_files);
        }
        else
            fclose(solvent_output_file);
    }

    free(solvent_coordination_history);
    free(solvent_atoms_coordination_history);
    free(anion_coordination_history);
    free(anion_atoms_coordination_history);
    free(venn_diagrams_cations);
    free(venn_diagrams_solvent);
    free(buffer);

    free_coordination_arrays(solvent_data, system_info, solvent);
    free_coordination_arrays(anion_data, system_info, anion);

    if(program_configuration->print_mode == separate)
    {    
        for(int i = 0; i < system_info->cations_number; i++)
        {
            fclose(cation_output_files[i]);
        }
        free(cation_output_files);
    }
    else
        fclose(cation_output_file);    

    fclose(solvent_times_file);
    fclose(anion_times_file);

    for(int i = 0; i < system_info->cations_number; i++)
    {
        free(cation_tracked_positions[i]);
    }

    int compound_index = 0;
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
        struct system_compound current_compound = system_info->compounds[compound_index];
        for(int j = 0; j < current_compound.quantity; j++)
        {
            free(solvent_tracked_positions[i][j]);
        }
        free(solvent_tracked_positions[i]);
        if(program_configuration->calculate_venn_diagrams == 1)
        { 
            free_index_combinations(solvent_index_combinations[i]);
            free(venn_set_cations[i]);
        }
    }

    for(int i = 0; i < system_info->anions_number; i++)
    {
        free(anion_tracked_positions[i]);
    }

    free(cation_tracked_positions);
    free(solvent_tracked_positions);
    free(anion_tracked_positions);
    if(program_configuration->calculate_venn_diagrams == 1) 
    {
        free(solvent_index_combinations);
        free(venn_set_cations);
        free(venn_set_solvent);
    }

    current_index = -1;
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        current_index = get_next_entry_index(current_index, system_info->compounds, system_info->compounds_number, solvent);
        struct system_compound current_compound = system_info->compounds[current_index];

        for(int j = 0; j < current_compound.tracked_atoms_number; j++)
        {
            free(current_coordinating_solvents[i][j]);
        }
        free(current_coordinating_solvents[i]);
    }
    free(current_coordinating_solvents);
    
    fclose(input_file);
}