#include <stdio.h>
#include <errno.h>

#include "program.h"

#define TRACKED_O_ATOMS_ANION 4

int detect_coordination(struct vector cation, struct vector ligand_atom, double threshold, double box_size)
{
    if(calculate_distance(&cation, &ligand_atom, box_size) <= threshold)
        return 1;
    
    return 0;
}

void mark_coordination(int cation_index, int* ligand_coordination)
{
    int index = 0;
    while(ligand_coordination[index] != BLANK) 
        index++;

    ligand_coordination[index] = cation_index;
}

/*
void insert_new_coord_info(struct oxygen_coord_info* array, int metal_index)
{
    int index = 0;
    while(array[index].metal_ion_index != BLANK && index < MAX_COORDINATED_METALS) index++;

    if(index == MAX_COORDINATED_METALS)
    { 
        errno = EINVAL;
        raise_error("No space left in coordination history array");
    }

    array[index].metal_ion_index = metal_index;
    array[index].time = 1;
}
*/

/*
int find_metal(struct oxygen_coord_info* array, int metal_index)
{
    int index = 0;
    while(array[index].metal_ion_index != metal_index && index < MAX_COORDINATED_METALS) index++;

    if(index == MAX_COORDINATED_METALS)
    {
        errno = EINVAL;
        raise_error("Invalid metal ion index to increase time");
    }

    return index;
}
*/

/*
void increase_coord_time(struct oxygen_coord_info* array, int metal_index)
{
    int index = find_metal(array, metal_index);
    array[index].time++;
}
*/

/*
void delete_coord_time(struct oxygen_coord_info* array, int metal_index, FILE* output_file)
{
    int index = find_metal(array, metal_index);
    fprintf(output_file, "%d\n", array[index].time);

    array[index].time = 0;
    array[index].metal_ion_index = BLANK;
}
*/

void swap_coordination_arrays(int*** first, int*** second)
{
    int** tmp = *(first);
    (*first) = *(second);
    *(second) = tmp;
}

void clear_coordination_array(int** array, int first_index_max, int second_index_max)
{
    for(int i = 0; i < first_index_max; i++)
    {
        for(int j = 0; j < second_index_max; j++)
        {
            array[i][j] = BLANK;
        }
    }
}

/*
void save_data_from_array(struct oxygen_coord_info* array, FILE* output_file)
{
    for(int i = 0; i < MAX_COORDINATED_METALS; i++)
    {
        if(array[i].metal_ion_index != BLANK)
            fprintf(output_file, "%d\n", array[i].time);
    }
}
*/

/*
void save_last_step_data(struct oxygen_coord_info** solvent_data, struct oxygen_coord_info*** anion_data, struct system_info* system_info, FILE* solvent_output, FILE* anion_output)
{
    for(int i = 0; i < system_info->solvent_molecules_number; i++)
    {
        save_data_from_array(solvent_data[i], solvent_output);
    }

    for(int i = 0; i < system_info->anions_number; i++)
    {
        for(int j = 0; j < TRACKED_O_ATOMS_ANION; j++)
        {
            save_data_from_array(anion_data[i][j], anion_output);
        }
    }
}
*/

/*
void calculate_coord_times(int molecules_number, int** current_coordination, int** last_coordination, 
                                      struct oxygen_coord_info** coordination_info, FILE* output_file)
{
    for(int i = 0; i < molecules_number; i++)
    {
        int last_index = 0;
        int current_index = 0;

        while(current_index < MAX_COORDINATED_METALS && last_index < MAX_COORDINATED_METALS &&
              last_coordination[i][last_index] != BLANK && current_coordination[i][current_index] != BLANK)
        {
            if(last_coordination[i][last_index] == current_coordination[i][current_index])
            {
                increase_coord_time(coordination_info[i], current_coordination[i][current_index]);
                last_index++;
                current_index++;
            }
            else if(last_coordination[i][last_index] < current_coordination[i][current_index])
            {
                delete_coord_time(coordination_info[i], last_coordination[i][last_index], output_file);
                last_index++;
            }
            else if(last_coordination[i][last_index] > current_coordination[i][current_index])
            {
                insert_new_coord_info(coordination_info[i], current_coordination[i][current_index]);
                current_index++;
            }
        }

        while(current_index < MAX_COORDINATED_METALS && current_coordination[i][current_index] != BLANK)
        {
            insert_new_coord_info(coordination_info[i], current_coordination[i][current_index]);
            current_index++;
        }

        while(last_index < MAX_COORDINATED_METALS && last_coordination[i][last_index] != BLANK)
        {
            delete_coord_time(coordination_info[i], last_coordination[i][last_index], output_file);
            last_index++;
        }
    }
}
*/

/*
void calculate_coord_times_anion(struct system_info* system_info, struct anion_data* anion_data)
{
    for(int i = 0; i < system_info->anions_number; i++)
    {
        calculate_coord_times(TRACKED_O_ATOMS_ANION, anion_data->current_anion_coordination[i], anion_data->last_anion_coordination[i], anion_data->anion_coordination_info[i], anion_data->anion_output);
    }
}
*/

/*
void calculate_coord_times_solvent(struct system_info* system_info, struct solvent_data* solvent_data)
{
    calculate_coord_times(system_info->solvent_molecules_number, solvent_data->current_solvent_coordination, solvent_data->last_solvent_coordination,
                          solvent_data->solvent_coordination_info, solvent_data->solvent_output);
}
*/

struct cation_coord_info get_coordination_info(int cation_index, int cation_tracked_positions_number, struct vector* cation_tracked_positions, struct coordination_input* coordination_input)
{
    struct cation_coord_info result;
    result.solvent_molecules = 0;
    result.anion_molecules = 0;
    result.anion_atoms = 0;
    result.coordination_number = 0;

    struct program_configuration* program_configuration = coordination_input->program_configuration;
    int solvent_history_index = 0;
    int anion_atoms_history_index = 0;
    int anion_history_index = 0;
    int solvent_molecule_index = 0;

    int solvent_shift = 0;
    int compound_index = 0;
    for(int i = 0; i < coordination_input->system_info->solvent_types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, coordination_input->system_info->compounds, coordination_input->system_info->compounds_number, solvent);
        for(int j = 0; j < coordination_input->system_info->compounds[compound_index].quantity; j++)
        {
            solvent_molecule_index = solvent_shift + j;

            for(int k = 0; k < coordination_input->system_info->compounds[compound_index].tracked_atoms_number; k++)
            {
                for(int l = 0; l < cation_tracked_positions_number; l++)
                {
                    if(detect_coordination(cation_tracked_positions[l], coordination_input->solvent_tracked_atoms[i][j][k], program_configuration->solvent_threshold, program_configuration->box_size) == 1)
                    {
                        result.coordination_number++;
                        result.solvent_molecules++;

                        if(program_configuration->calculate_solvent_residence == 1)
                        {
                            coordination_input->solvent_coordination_history[i][coordination_input->step_number][cation_index][solvent_history_index] = j;
                            solvent_history_index++;
                        }
                        
                        mark_coordination(cation_index, coordination_input->current_solvent_coordination[solvent_molecule_index][k]);
                    }
                }
            }
        }

        solvent_history_index = 0;
        solvent_shift += coordination_input->system_info->compounds[compound_index].quantity;
    }

    int anion_shift = 0;
    compound_index = 0;
    for(int i = 0; i < coordination_input->system_info->anion_types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, coordination_input->system_info->compounds, coordination_input->system_info->compounds_number, anion);
        struct system_compound current_compound = coordination_input->system_info->compounds[compound_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            int anion_index = anion_shift + j;
            int is_cation_coordinated_by_anion = 0;

            for(int k = 0; k < current_compound.tracked_atoms_number; k++)
            {
                for(int l = 0; l < cation_tracked_positions_number; l++)
                {
                    if(detect_coordination(cation_tracked_positions[l], coordination_input->anion_tracked_atoms[anion_index][k], program_configuration->anion_threshold, program_configuration->box_size) == 1)
                    {
                        result.coordination_number++;
                        result.anion_atoms++;
                        is_cation_coordinated_by_anion = 1;

                        if(program_configuration->calculate_anion_residence == 1)
                        {
                            coordination_input->anion_atoms_coordination_history[coordination_input->step_number][cation_index][anion_atoms_history_index] = anion_index * current_compound.tracked_atoms_number + k;
                            anion_atoms_history_index++;
                        }

                        mark_coordination(cation_index, coordination_input->current_anion_coordination[anion_index][k]);
                    }
                }
            }

            if(is_cation_coordinated_by_anion == 1)
            {
                result.anion_molecules++;

                if(program_configuration->calculate_anion_residence == 1)
                {
                    coordination_input->anion_coordination_history[coordination_input->step_number][cation_index][anion_history_index] = anion_index;
                    anion_history_index++;
                }
            }
        }

        anion_shift += current_compound.quantity;
    }

    return result;
}