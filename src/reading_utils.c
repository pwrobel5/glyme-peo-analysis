#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "program.h"

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

void read_carbonate_data(FILE* input_file, struct vector* carbonate_oxygen_positions, struct system_compound carbonate)
{
    char* tmp;
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < carbonate.quantity; i++)
    {
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        tmp = strtok(buffer, SEPARATOR);

        check_symbol(tmp, carbonate.first_atom_symbol);

        int lines_to_omit = OMITTED_ATOMS_CARBONATE - 1;
        for(int j = 0; j < lines_to_omit; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i;
        read_vector_coordinates(buffer, "O", carbonate_oxygen_positions + index);

        lines_to_omit = carbonate.atoms_number - OMITTED_ATOMS_CARBONATE - 1;
        for(int j = 0; j < lines_to_omit; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }
    }

    free(buffer);
}

void read_tfsi_data(FILE* input_file, struct vector** tfsi_oxygen_positions, struct system_compound tfsi, int* shift)
{
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL) raise_error("Error with buffer allocation");

    for(int i = 0; i < tfsi.quantity; i++)
    {
        // first oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        int index = i + *shift;        
        read_vector_coordinates(buffer, "O", &(tfsi_oxygen_positions[index][0]));

        // omit next 7 atoms
        for(int j = 0; j < 7; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // second oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(tfsi_oxygen_positions[index][1]));

        // omit next 3 atoms
        for(int j = 0; j < 3; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // third oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(tfsi_oxygen_positions[index][2]));

        // omit next 3 atoms
        for(int j = 0; j < 3; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }

        // fourth oxygen atom
        fgets(buffer, MAX_LINE_LENGTH, input_file);
        read_vector_coordinates(buffer, "O", &(tfsi_oxygen_positions[index][3]));

        // omit next 13 atoms
        for(int j = 0; j < 13; j++)
        {
            fgets(buffer, MAX_LINE_LENGTH, input_file);
        }
    }

    *shift += tfsi.quantity;
    free(buffer);
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

FILE** open_metal_files(struct system_info* system_info)
{
    FILE** result = malloc(system_info->metal_ions_number * sizeof(FILE*));
    if(result == NULL) raise_error("Error with array for outpu files allocation");

    char filename_buffer[MAX_LINE_LENGTH];

    int shift = 0;
    for(int i = 0; i < system_info->compounds_number; i++)
    {
        struct system_compound compound = system_info->compounds[i];
        if(compound.molecule_type == met)
        {
            for(int j = 0; j < compound.quantity; j++)
            {
                sprintf(filename_buffer, "%s_%02d.dat", compound.first_atom_symbol, j);

                result[shift] = fopen(filename_buffer, "w");
                if(result[shift] == NULL) raise_error("Error with creating output file for metal ions");

                shift++;
            }
        }
    }

    return result;
}

void initialize_carbonate_coordination_arrays(int** last_carbonate_coordination, int** current_carbonate_coordination, 
                                              struct oxygen_coord_info** carbonate_coordination_info, struct system_info* system_info)
{
    if(last_carbonate_coordination == NULL || current_carbonate_coordination == NULL || carbonate_coordination_info == NULL) 
        raise_error("Error with memory allocation for coordination history");
    for(int i = 0; i < system_info->carbonate_molecules_number; i++)
    {
        last_carbonate_coordination[i] = malloc(MAX_COORDINATED_METALS * sizeof(int));
        current_carbonate_coordination[i] = malloc(MAX_COORDINATED_METALS * sizeof(int));
        carbonate_coordination_info[i] = malloc(MAX_COORDINATED_METALS * sizeof(struct oxygen_coord_info));

        if(last_carbonate_coordination[i] == NULL || current_carbonate_coordination[i] == NULL || carbonate_coordination_info[i] == NULL)
            raise_error("Error with memory allocation for coordination history");
        
        for(int j = 0; j < MAX_COORDINATED_METALS; j++)
        {
            last_carbonate_coordination[i][j] = BLANK;
            current_carbonate_coordination[i][j] = BLANK;
            carbonate_coordination_info[i][j].metal_ion_index = BLANK;
            carbonate_coordination_info[i][j].time = 0;
        }
    }
}

void free_carbonate_coordination_arrays(int** last_carbonate_coordination, int** current_carbonate_coordination, 
                                              struct oxygen_coord_info** carbonate_coordination_info, struct system_info* system_info)
{
    for(int i = 0; i < system_info->carbonate_molecules_number; i++)
    {
        free(last_carbonate_coordination[i]);
        free(current_carbonate_coordination[i]);
        free(carbonate_coordination_info[i]);
    }
    free(last_carbonate_coordination);
    free(current_carbonate_coordination);
    free(carbonate_coordination_info);
}

void initialize_tfsi_coordination_arrays(int*** last_tfsi_coordination, int*** current_tfsi_coordination, struct oxygen_coord_info*** tfsi_coordination_info, 
                                         struct system_info* system_info)
{
    if(last_tfsi_coordination == NULL || current_tfsi_coordination == NULL || tfsi_coordination_info == NULL)
        raise_error("Error with memory allocation for coordination history");
    
    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        last_tfsi_coordination[i] = malloc(TRACKED_O_ATOMS_TFSI * sizeof(int**));
        current_tfsi_coordination[i] = malloc(TRACKED_O_ATOMS_TFSI * sizeof(int**));
        tfsi_coordination_info[i] = malloc(TRACKED_O_ATOMS_TFSI * sizeof(struct oxygen_coord_info**));

        if(last_tfsi_coordination[i] == NULL || current_tfsi_coordination[i] == NULL || tfsi_coordination_info[i] == NULL)
            raise_error("Error with memory allocation for coordination history");
        
        for(int j = 0; j < TRACKED_O_ATOMS_TFSI; j++)
        {
            last_tfsi_coordination[i][j] = malloc(MAX_COORDINATED_METALS * sizeof(int));
            current_tfsi_coordination[i][j] = malloc(MAX_COORDINATED_METALS * sizeof(int));
            tfsi_coordination_info[i][j] = malloc(MAX_COORDINATED_METALS * sizeof(struct oxygen_coord_info));

            if(last_tfsi_coordination[i][j] == NULL || current_tfsi_coordination[i][j] == NULL || tfsi_coordination_info[i][j] == NULL)
                raise_error("Error with memory allocation for coordination history");
            
            for(int k = 0; k < MAX_COORDINATED_METALS; k++)
            {
                last_tfsi_coordination[i][j][k] = BLANK;
                current_tfsi_coordination[i][j][k] = BLANK;
                tfsi_coordination_info[i][j][k].metal_ion_index = BLANK;
                tfsi_coordination_info[i][j][k].time = 0;
            }
        }
    }
}

void free_tfsi_coordination_arrays(int*** last_tfsi_coordination, int*** current_tfsi_coordination, struct oxygen_coord_info*** tfsi_coordination_info, 
                                         struct system_info* system_info)
{
    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        for(int j = 0; j < TRACKED_O_ATOMS_TFSI; j++)
        {
            free(last_tfsi_coordination[i][j]);
            free(current_tfsi_coordination[i][j]);
            free(tfsi_coordination_info[i][j]);
        }

        free(last_tfsi_coordination[i]);
        free(current_tfsi_coordination[i]);
        free(tfsi_coordination_info[i]);
    }

    free(last_tfsi_coordination);
    free(current_tfsi_coordination);
    free(tfsi_coordination_info);
}

void read_data(struct program_configuration* program_configuration, struct system_info* system_info)
{
    FILE* input_file = fopen(program_configuration->input_file_name, "r");

    if(input_file == NULL)
        raise_error("Error with opening input file!");
    
    int step_number = 0;
    struct vector* metal_ions_positions = malloc(system_info->metal_ions_number * sizeof(struct vector));
    struct vector** carbonate_oxygen_positions = malloc(system_info->carbonate_types_number * sizeof(struct vector*));
    struct vector** tfsi_oxygen_positions = malloc(system_info->tfsi_molecules_number * sizeof(struct vector*));
    if(metal_ions_positions == NULL || carbonate_oxygen_positions == NULL || tfsi_oxygen_positions == NULL)
        raise_error("Error with allocation of arrays for atoms positions");
    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        tfsi_oxygen_positions[i] = malloc(TRACKED_O_ATOMS_TFSI * sizeof(struct vector));
        if(tfsi_oxygen_positions[i] == NULL)
            raise_error("Error with memory allocation for TFSI oxygen atoms positions");
    }

    int compound_index = 0;
    for(int i = 0; i < system_info->carbonate_types_number; i++)
    {
        compound_index = get_next_carbonate_index(compound_index, system_info->compounds, system_info->compounds_number);
        struct system_compound current_compound = system_info->compounds[compound_index];

        int oxygen_atoms = current_compound.atoms_number * current_compound.quantity;
        carbonate_oxygen_positions[i] = malloc(oxygen_atoms * sizeof(struct vector));
        if(carbonate_oxygen_positions[i] == NULL)
            raise_error("Error with memory allocation for carbonate oxygen atoms positions");
    }

    FILE** metal_output_files;
    FILE* metal_output_file;
    FILE* carbonate_output_file = fopen("carbonate.dat", "w");
    FILE* tfsi_output_file = fopen("tfsi.dat", "w");

    if(program_configuration->print_mode == separate)
        metal_output_files = open_metal_files(system_info);
    else
        metal_output_file = fopen("cation_output.dat", "w");

    int** last_carbonate_coordination = malloc(system_info->carbonate_molecules_number * sizeof(int*));
    int** current_carbonate_coordination = malloc(system_info->carbonate_molecules_number * sizeof(int*));
    struct oxygen_coord_info** carbonate_coordination_info = malloc(system_info->carbonate_molecules_number * sizeof(struct oxygen_coord_info*));
    initialize_carbonate_coordination_arrays(last_carbonate_coordination, current_carbonate_coordination, carbonate_coordination_info, system_info);
    
    struct carbonate_data carbonate_data;
    carbonate_data.carbonate_output = carbonate_output_file;
    carbonate_data.last_carbonate_coordination = last_carbonate_coordination;
    carbonate_data.current_carbonate_coordination = current_carbonate_coordination;
    carbonate_data.carbonate_coordination_info = carbonate_coordination_info;

    int*** last_tfsi_coordination = malloc(system_info->tfsi_molecules_number * sizeof(int**));
    int*** current_tfsi_coordination = malloc(system_info->tfsi_molecules_number * sizeof(int**));
    struct oxygen_coord_info*** tfsi_coordination_info = malloc(system_info->tfsi_molecules_number * sizeof(struct oxygen_coord_info*));
    initialize_tfsi_coordination_arrays(last_tfsi_coordination, current_tfsi_coordination, tfsi_coordination_info, system_info);

    struct tfsi_data tfsi_data;
    tfsi_data.tfsi_output = tfsi_output_file;
    tfsi_data.last_tfsi_coordination = last_tfsi_coordination;
    tfsi_data.current_tfsi_coordination = current_tfsi_coordination;
    tfsi_data.tfsi_coordination_info = tfsi_coordination_info;
    
    char* buffer = malloc(MAX_LINE_LENGTH * sizeof(char));
    if(buffer == NULL)
        raise_error("Error with memory allocation for reading buffer");
    
    // arrays for saving numbers of coordinated molecules to given metal ion at a given timestep
    short int*** tfsi_atoms_coordination_history = NULL;
    short int*** tfsi_coordination_history = NULL;
    short int**** carbonate_coordination_histories = malloc(sizeof(short int***) * system_info->carbonate_types_number);
    if(carbonate_coordination_histories == NULL) raise_error("Error with memory allocation for carbonate coordination histories array");
    for(int i = 0; i < system_info->carbonate_types_number; i++)
    {
        carbonate_coordination_histories[i] = NULL;
    }

    while(fgets(buffer, MAX_LINE_LENGTH, input_file) != NULL)
    {
        int metal_shift = 0;
        int tfsi_shift = 0;

        int atoms_number = atoi(buffer);
        if(atoms_number != system_info->atoms_number)
            printf("[WARNING] Declared atoms number different from number read from .xyz file, declared: %d, read: %d\n", system_info->atoms_number, atoms_number);
        
        if(program_configuration->calculate_carbonate_residence == 1)
        {    
            for(int i = 0; i < system_info->carbonate_types_number; i++)
                carbonate_coordination_histories[i] = initialize_history_array(carbonate_coordination_histories[i], step_number, system_info);
        }
        
        if(program_configuration->calculate_tfsi_residence == 1)
        {
            tfsi_atoms_coordination_history = initialize_history_array(tfsi_atoms_coordination_history, step_number, system_info);
            tfsi_coordination_history = initialize_history_array(tfsi_coordination_history, step_number, system_info);
        }
        
        // omit commentary line
        fgets(buffer, MAX_LINE_LENGTH, input_file);

        int carbonate_index = 0;
        for(int i = 0; i < system_info->compounds_number; i++)
        {
            struct system_compound compound = system_info->compounds[i];

            switch(compound.molecule_type)
            {
                case met:
                    read_metal_data(input_file, metal_ions_positions, compound, &metal_shift);
                    break;
                case ec:
                case f1ec:
                case f2ec:
                    read_carbonate_data(input_file, carbonate_oxygen_positions[carbonate_index], compound);
                    carbonate_index++;
                    break;
                case tfsi:
                    read_tfsi_data(input_file, tfsi_oxygen_positions, compound, &tfsi_shift);
                    break;
                case other:
                    omit_other_data(input_file, compound);
                    break;
            }
        }

        struct coordination_input coordination_input;
        coordination_input.step_number = step_number;
        coordination_input.carbonate_oxygens = carbonate_oxygen_positions;
        coordination_input.tfsi_oxygens = tfsi_oxygen_positions;
        coordination_input.current_carbon_coordination = carbonate_data.current_carbonate_coordination;
        coordination_input.current_tfsi_coordination = tfsi_data.current_tfsi_coordination;
        coordination_input.carbonate_coordination_histories = carbonate_coordination_histories;
        coordination_input.tfsi_atoms_coordination_history = tfsi_atoms_coordination_history;
        coordination_input.tfsi_coordination_history = tfsi_coordination_history;
        coordination_input.program_configuration = program_configuration;
        coordination_input.system_info = system_info;
        

        for(int i = 0; i < system_info->metal_ions_number; i++)
        {
            struct metal_coord_info coord_info = get_coordination_info(i, metal_ions_positions[i], &coordination_input);

            if(program_configuration->print_mode == separate)
                fprintf(metal_output_files[i], "%d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.carbonate_molecules, coord_info.tfsi_oxygens, coord_info.tfsi_molecules);
            else
                fprintf(metal_output_file, "%d %d %d %d %d\n", step_number, coord_info.coordination_number, coord_info.carbonate_molecules, coord_info.tfsi_oxygens, coord_info.tfsi_molecules);
        }

        calculate_coord_times_carbonates(system_info, &carbonate_data);
        calculate_coord_times_tfsi(system_info, &tfsi_data);

        swap_coordination_arrays(&(carbonate_data.last_carbonate_coordination), &(carbonate_data.current_carbonate_coordination));
        clear_coordination_array(carbonate_data.current_carbonate_coordination, system_info->carbonate_molecules_number);

        for(int i = 0; i < system_info->tfsi_molecules_number; i++)
        {
            swap_coordination_arrays(&(tfsi_data.last_tfsi_coordination[i]), &(tfsi_data.current_tfsi_coordination[i]));
            clear_coordination_array(tfsi_data.current_tfsi_coordination[i], TRACKED_O_ATOMS_TFSI);
        }

        step_number++;
    }

    save_last_step_data(carbonate_coordination_info, tfsi_coordination_info, system_info, carbonate_output_file, tfsi_output_file);

    // calculate correlated residence times
    if(program_configuration->calculate_carbonate_residence == 1)
    {
        int compound_index = 0;
        for(int i = 0; i < system_info->carbonate_types_number; i++) {
            compound_index = get_next_carbonate_index(compound_index, system_info->compounds, system_info->compounds_number);
            struct system_compound current_compound = system_info->compounds[compound_index];

            int denominator = system_info->metal_ions_number * system_info->compounds[compound_index].quantity;
            double* residence = calculate_residence_times(carbonate_coordination_histories[i], step_number, denominator, system_info);
            delete_history_array(carbonate_coordination_histories[i], step_number, system_info->metal_ions_number);

            char residence_output_name[OUTPUT_FILE_NAME_LENGTH] = "carbonate-residence-times-";
            char* molecule_name = molecule_type_to_str(current_compound.molecule_type);
            strcat(residence_output_name, molecule_name);
            strcat(residence_output_name, ".dat");
            save_residence_to_file(residence, residence_output_name, step_number);
            free(residence);
            free(molecule_name);
        }
    }
    if(program_configuration->calculate_tfsi_residence == 1)
    {
        int molecules_denominator = system_info->metal_ions_number * system_info->tfsi_molecules_number;
        double* molecules_residence = calculate_residence_times(tfsi_coordination_history, step_number, molecules_denominator, system_info);
        delete_history_array(tfsi_coordination_history, step_number, system_info->metal_ions_number);
        save_residence_to_file(molecules_residence, "tfsi-residence-times.dat", step_number);
        free(molecules_residence);

        int atoms_denominator = molecules_denominator * TRACKED_O_ATOMS_TFSI;
        double* atoms_residence = calculate_residence_times(tfsi_atoms_coordination_history, step_number, atoms_denominator, system_info);
        delete_history_array(tfsi_atoms_coordination_history, step_number, system_info->metal_ions_number);
        save_residence_to_file(atoms_residence, "tfsi-atoms-residence-times.dat", step_number);
        free(atoms_residence);        
    }

    free_carbonate_coordination_arrays(last_carbonate_coordination, current_carbonate_coordination, carbonate_coordination_info, system_info);
    free_tfsi_coordination_arrays(last_tfsi_coordination, current_tfsi_coordination, tfsi_coordination_info, system_info);

    if(program_configuration->print_mode == separate)
    {    
        for(int i = 0; i < system_info->metal_ions_number; i++)
        {
            fclose(metal_output_files[i]);
        }
        free(metal_output_files);
    }
    else
        fclose(metal_output_file);    

    fclose(carbonate_output_file);
    fclose(tfsi_output_file);

    free(metal_ions_positions);
    
    for(int i = 0; i < system_info->carbonate_types_number; i++)
    {
        free(carbonate_oxygen_positions[i]);
    }
    free(carbonate_oxygen_positions);

    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        free(tfsi_oxygen_positions[i]);
    }
    free(tfsi_oxygen_positions);
    free(carbonate_coordination_histories);
    free(buffer);
    fclose(input_file);
}