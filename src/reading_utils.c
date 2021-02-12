#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "program.h"

#define OUTPUT_FILE_NAME_LENGTH 60
#define OMITTED_ATOMS_CARBONATE 10
#define TRACKED_O_ATOMS_ANION 4
/*
void check_symbol(char* read, char* expected)
{
    if(strcmp(read, expected) != 0)
    {
        errno = EINVAL;
        char buffer[MAX_LINE_LENGTH];
        sprintf(buffer, "Unexpected atom symbol! Expected %s, read %s", expected, read);
        raise_error(buffer);
    }
}
*/

/*
void read_metal_data(FILE* input_file, struct vector* metal_ions_positions, struct system_compound compound, int* shift)
{
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < compound.quantity; i++)
    {
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i + *shift;
        read_vector_coordinates(buffer, compound.first_atom_symbol, metal_ions_positions + index);
    }

    *shift += compound.quantity;
    free(buffer);
}
*/

/*
void read_solvent_data(FILE* input_file, struct vector* solvent_oxygen_positions, struct system_compound solvent)
{
    char* tmp;
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < solvent.quantity; i++)
    {
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        tmp = strtok(buffer, SEPARATOR);

        check_symbol(tmp, solvent.first_atom_symbol);

        int lines_to_omit = OMITTED_ATOMS_CARBONATE - 1;
        for(int j = 0; j < lines_to_omit; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i;
        read_vector_coordinates(buffer, "O", solvent_oxygen_positions + index);

        lines_to_omit = solvent.atoms_number - OMITTED_ATOMS_CARBONATE - 1;
        for(int j = 0; j < lines_to_omit; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }
    }

    free(buffer);
}
*/

/*
void read_anion_data(FILE* input_file, struct vector** anion_oxygen_positions, struct system_compound anion, int* shift)
{
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < anion.quantity; i++)
    {
        // first oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i + *shift;        
        read_vector_coordinates(buffer, "O", &(anion_oxygen_positions[index][0]));

        // omit next 7 atoms
        for(int j = 0; j < 7; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // second oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(anion_oxygen_positions[index][1]));

        // omit next 3 atoms
        for(int j = 0; j < 3; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // third oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(anion_oxygen_positions[index][2]));

        // omit next 3 atoms
        for(int j = 0; j < 3; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // fourth oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(anion_oxygen_positions[index][3]));

        // omit next 13 atoms
        for(int j = 0; j < 13; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }
    }

    *shift += anion.quantity;
    free(buffer);
}
*/

/*
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
*/

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

    int current_entry_index = 0;
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

    int current_entry_index = 0;
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

void read_data(struct program_configuration* program_configuration, struct system_info* system_info)
{
    FILE* input_file = fopen(program_configuration->input_file_name, "r");

    if(input_file == NULL) raise_error("Error with opening input file!");
    
    struct vector** cation_tracked_positions = malloc(system_info->cations_number * sizeof(struct vector*));
    struct vector*** solvent_tracked_positions = malloc(system_info->solvent_types_number * sizeof(struct vector**));
    struct vector** anion_tracked_positions = malloc(system_info->anions_number * sizeof(struct vector*));
    if(cation_tracked_positions == NULL || solvent_tracked_positions == NULL || anion_tracked_positions == NULL)
        raise_error("Error with allocation of arrays for atoms positions");
    
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
                solvent_tracked_positions[solvent_type_index] = malloc(current_compound.quantity * sizeof(struct vector**));
                if(solvent_tracked_positions[solvent_type_index] == NULL) raise_error("Error with memory allocation for solvent tracked positions");
                for(int j = 0; j < current_compound.quantity; j++)
                {
                    solvent_tracked_positions[solvent_type_index][j] = malloc(current_compound.tracked_atoms_number * sizeof(struct vector));
                    if(solvent_tracked_positions[solvent_type_index][j] == NULL) raise_error("Error with memory allocation for solvent tracked positions");
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
    FILE* solvent_output_file = fopen("solvent_output.dat", "w");
    FILE* anion_output_file = fopen("anion_output.dat", "w");
    if(solvent_output_file == NULL || anion_output_file == NULL) raise_error("Error with opening output files");

    if(program_configuration->print_mode == separate)
        cation_output_files = open_cation_files(system_info);
    else
        cation_output_file = fopen("cation_output.dat", "w");

    int*** last_solvent_coordination = malloc(system_info->solvent_molecules_number * sizeof(int**));
    int*** current_solvent_coordination = malloc(system_info->solvent_molecules_number * sizeof(int**));
    struct ligand_coord_info*** solvent_coordination_info = malloc(system_info->solvent_molecules_number * sizeof(struct ligand_coord_info**));
    
    struct entry_data solvent_data;
    solvent_data.output = solvent_output_file;
    solvent_data.last_coordination = last_solvent_coordination;
    solvent_data.current_coordination = current_solvent_coordination;
    solvent_data.coordination_info = solvent_coordination_info;
    initialize_entry_data(solvent_data, system_info, solvent);

    int*** last_anion_coordination = malloc(system_info->anions_number * sizeof(int**));
    int*** current_anion_coordination = malloc(system_info->anions_number * sizeof(int**));
    struct ligand_coord_info*** anion_coordination_info = malloc(system_info->anions_number * sizeof(struct oxygen_coord_info*));

    struct entry_data anion_data;
    anion_data.output = anion_output_file;
    anion_data.last_coordination = last_anion_coordination;
    anion_data.current_coordination = current_anion_coordination;
    anion_data.coordination_info = anion_coordination_info;
    initialize_entry_data(anion_data, system_info, anion);
    
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL)
        raise_error("Error with memory allocation for reading buffer");
    
    // arrays for saving numbers of coordinated molecules to given metal ion at a given timestep
    short int*** anion_atoms_coordination_history = NULL;
    short int*** anion_coordination_history = NULL;
    short int**** solvent_coordination_history = malloc(sizeof(short int***) * system_info->solvent_types_number);
    if(solvent_coordination_history == NULL) raise_error("Error with memory allocation for carbonate coordination histories array");
    for(int i = 0; i < system_info->solvent_types_number; i++)
    {
        solvent_coordination_history[i] = NULL;
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
                solvent_coordination_history[i] = initialize_history_array(solvent_coordination_history[i], step_number, system_info);
        }
        
        if(program_configuration->calculate_anion_residence == 1)
        {
            anion_atoms_coordination_history = initialize_history_array(anion_atoms_coordination_history, step_number, system_info);
            anion_coordination_history = initialize_history_array(anion_coordination_history, step_number, system_info);
        }

        // omit commentary line
        fgets(buffer, MAX_LINE_LENGTH, input_file);
/*
        int carbonate_index = 0;
        for(int i = 0; i < system_info->compounds_number; i++)
        {
            struct system_compound compound = system_info->compounds[i];

            switch(compound.molecule_type)
            {
                case met:
                    read_metal_data(input_file, metal_ions_positions, compound, &cation_shift);
                    break;
                case ec:
                case f1ec:
                case f2ec:
                case monoglym:
                case tetraglym:
                case peo:
                    read_solvent_data(input_file, carbonate_oxygen_positions[carbonate_index], compound);
                    carbonate_index++;
                    break;
                case fsi:
                case tfsi:
                    read_anion_data(input_file, tfsi_oxygen_positions, compound, &anion_shift);
                    break;
                case other:
                    omit_other_data(input_file, compound);
                    break;
            }
        }

        struct coordination_input coordination_input;
        coordination_input.step_number = step_number;
        coordination_input.solvent_oxygens = carbonate_oxygen_positions;
        coordination_input.anion_oxygens = tfsi_oxygen_positions;
        coordination_input.current_solvent_coordination = solvent_data.current_solvent_coordination;
        coordination_input.current_anion_coordination = anion_data.current_anion_coordination;
        coordination_input.solvent_coordination_histories = carbonate_coordination_histories;
        coordination_input.anion_atoms_coordination_history = tfsi_atoms_coordination_history;
        coordination_input.anion_coordination_history = tfsi_coordination_history;
        coordination_input.program_configuration = program_configuration;
        coordination_input.system_info = system_info;
        

        for(int i = 0; i < system_info->metal_ions_number; i++)
        {
            struct metal_coord_info coord_info = get_coordination_info(i, metal_ions_positions[i], &coordination_input);

            if(program_configuration->print_mode == separate)
                fprintf(metal_output_files[i], "%d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.solvent_molecules, coord_info.anion_oxygens, coord_info.anion_molecules);
            else
                fprintf(metal_output_file, "%d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.solvent_molecules, coord_info.anion_oxygens, coord_info.anion_molecules);
        }

        calculate_coord_times_solvent(system_info, &solvent_data);
        calculate_coord_times_anion(system_info, &anion_data);

        swap_coordination_arrays(&(solvent_data.last_solvent_coordination), &(solvent_data.current_solvent_coordination));
        clear_coordination_array(solvent_data.current_solvent_coordination, system_info->solvent_molecules_number);

        for(int i = 0; i < system_info->anions_number; i++)
        {
            swap_coordination_arrays(&(anion_data.last_anion_coordination[i]), &(anion_data.current_anion_coordination[i]));
            clear_coordination_array(anion_data.current_anion_coordination[i], TRACKED_O_ATOMS_ANION);
        }
*/
        step_number++;
    }
/*
    save_last_step_data(carbonate_coordination_info, tfsi_coordination_info, system_info, carbonate_output_file, tfsi_output_file);
*/
    // calculate correlated residence times
    if(program_configuration->calculate_solvent_residence == 1)
    {
        int compound_index = 0;
        for(int i = 0; i < system_info->solvent_types_number; i++) {
            compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, solvent);
            struct system_compound current_compound = system_info->compounds[compound_index];

            /*
            int denominator = system_info->metal_ions_number * system_info->compounds[compound_index].quantity;
            double* residence = calculate_residence_times(carbonate_coordination_histories[i], step_number, denominator, system_info);
            */
            delete_history_array(solvent_coordination_history[i], step_number, system_info->cations_number);

            /*
            char residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "carbonate-residence-times-";
            char* molecule_name = entry_type_to_str(current_compound.molecule_type);
            strcat(residence_output_name, molecule_name);
            strcat(residence_output_name, ".dat");
            save_residence_to_file(residence, residence_output_name, step_number);
            free(residence);
            free(molecule_name);
            */
        }
    }
    if(program_configuration->calculate_anion_residence == 1)
    {
        /*
        int molecules_denominator = system_info->metal_ions_number * system_info->anions_number;
        double* molecules_residence = calculate_residence_times(tfsi_coordination_history, step_number, molecules_denominator, system_info);
        */
        delete_history_array(anion_coordination_history, step_number, system_info->cations_number);
        /*
        save_residence_to_file(molecules_residence, "tfsi-residence-times.dat", step_number);
        free(molecules_residence);

        int atoms_denominator = molecules_denominator * TRACKED_O_ATOMS_ANION;
        double* atoms_residence = calculate_residence_times(tfsi_atoms_coordination_history, step_number, atoms_denominator, system_info);*/
        delete_history_array(anion_atoms_coordination_history, step_number, system_info->cations_number);
        /*save_residence_to_file(atoms_residence, "tfsi-atoms-residence-times.dat", step_number);
        free(atoms_residence);  */      
    }

    free(solvent_coordination_history);
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

    fclose(solvent_output_file);
    fclose(anion_output_file);

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
    }

    for(int i = 0; i < system_info->anions_number; i++)
    {
        free(anion_tracked_positions[i]);
    }

    free(cation_tracked_positions);
    free(solvent_tracked_positions);
    free(anion_tracked_positions);
    
    fclose(input_file);
}