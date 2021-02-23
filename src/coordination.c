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

void insert_new_coord_info(struct ligand_coord_info* array, int cation_index)
{
    int index = 0;
    while(array[index].cation_index != BLANK && index < MAX_COORDINATED_CATIONS) index++;

    if(index == MAX_COORDINATED_CATIONS)
    { 
        errno = EINVAL;
        raise_error("No space left in coordination history array");
    }

    array[index].cation_index = cation_index;
    array[index].time = 1;
}

int find_cation(struct ligand_coord_info* array, int cation_index)
{
    int index = 0;
    while(array[index].cation_index != cation_index && index < MAX_COORDINATED_CATIONS) index++;

    if(index == MAX_COORDINATED_CATIONS)
    {
        errno = EINVAL;
        raise_error("Invalid metal ion index to increase time");
    }

    return index;
}

void increase_coord_time(struct ligand_coord_info* array, int cation_index)
{
    int index = find_cation(array, cation_index);
    array[index].time++;
}

void delete_coord_time(struct ligand_coord_info* array, int cation_index, FILE* output_file)
{
    int index = find_cation(array, cation_index);
    fprintf(output_file, "%d\n", array[index].time);

    array[index].time = 0;
    array[index].cation_index = BLANK;
}

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

void save_data_from_array(struct ligand_coord_info* array, FILE* output_file)
{
    for(int i = 0; i < MAX_COORDINATED_CATIONS; i++)
    {
        if(array[i].cation_index != BLANK)
            fprintf(output_file, "%d\n", array[i].time);
    }
}

void save_last_step_data(struct system_info* system_info, struct entry_data* entry_data, enum entry_type entry_type)
{
    int shift = 0;
    int compound_index = -1;
    int types_number = (entry_type == solvent) ? system_info->solvent_types_number : system_info->anion_types_number;

    for(int i = 0; i < types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, entry_type);
        struct system_compound current_compound = system_info->compounds[compound_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            int current_index = j + shift;
            for(int k = 0; k < current_compound.tracked_atoms_number; k++)
            {
                save_data_from_array(entry_data->coordination_info[current_index][k], entry_data->output);
            }
        }
    }
}

void calculate_coord_times_molecule(int tracked_atoms_number, int** current_coordination, int** last_coordination, 
                                      struct ligand_coord_info** coordination_info, FILE* output_file)
{
    for(int i = 0; i < tracked_atoms_number; i++)
    {
        int last_index = 0;
        int current_index = 0;

        while(current_index < MAX_COORDINATED_CATIONS && last_index < MAX_COORDINATED_CATIONS &&
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

        while(current_index < MAX_COORDINATED_CATIONS && current_coordination[i][current_index] != BLANK)
        {
            insert_new_coord_info(coordination_info[i], current_coordination[i][current_index]);
            current_index++;
        }

        while(last_index < MAX_COORDINATED_CATIONS && last_coordination[i][last_index] != BLANK)
        {
            delete_coord_time(coordination_info[i], last_coordination[i][last_index], output_file);
            last_index++;
        }
    }
}

void calculate_coord_times(struct system_info* system_info, struct entry_data* entry_data, enum entry_type entry_type)
{
    int shift = 0;
    int compound_index = -1;
    int types_number = (entry_type == solvent) ? system_info->solvent_types_number : system_info->anion_types_number;

    for(int i = 0; i < types_number; i++)
    {
        compound_index = get_next_entry_index(compound_index, system_info->compounds, system_info->compounds_number, entry_type);
        struct system_compound current_compound = system_info->compounds[compound_index];

        for(int j = 0; j < current_compound.quantity; j++)
        {
            int current_index = j + shift;
            calculate_coord_times_molecule(current_compound.tracked_atoms_number, entry_data->current_coordination[current_index], entry_data->last_coordination[current_index],
                entry_data->coordination_info[current_index], entry_data->output);
        }

        shift += current_compound.quantity;
    }
}

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
                            coordination_input->solvent_coordination_history[i][coordination_input->step_number][cation_index][solvent_history_index] = solvent_molecule_index;
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
                            coordination_input->anion_atoms_coordination_history[i][coordination_input->step_number][cation_index][anion_atoms_history_index] = anion_index * current_compound.tracked_atoms_number + k;
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
                    coordination_input->anion_coordination_history[i][coordination_input->step_number][cation_index][anion_history_index] = anion_index;
                    anion_history_index++;
                }
            }
        }

        anion_shift += current_compound.quantity;
    }

    return result;
}