#include <stdio.h>
#include <errno.h>

#include "program.h"

int detect_coordination(struct vector metal, struct vector oxygen, double threshold, double box_size)
{
    if(calculate_distance(&metal, &oxygen, box_size) <= threshold)
        return 1;
    
    return 0;
}

void mark_coordination(int metal_index, int* oxygen_coordination)
{
    int index = 0;
    while(oxygen_coordination[index] != BLANK) 
        index++;

    oxygen_coordination[index] = metal_index;
}

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

void increase_coord_time(struct oxygen_coord_info* array, int metal_index)
{
    int index = find_metal(array, metal_index);
    array[index].time++;
}

void delete_coord_time(struct oxygen_coord_info* array, int metal_index, FILE* output_file)
{
    int index = find_metal(array, metal_index);
    fprintf(output_file, "%d\n", array[index].time);

    array[index].time = 0;
    array[index].metal_ion_index = BLANK;
}

void swap_coordination_arrays(int*** first, int*** second)
{
    int** tmp = *(first);
    (*first) = *(second);
    *(second) = tmp;
}

void clear_coordination_array(int** array, int first_index_max)
{
    for(int i = 0; i < first_index_max; i++)
    {
        for(int j = 0; j < MAX_COORDINATED_METALS; j++)
        {
            array[i][j] = BLANK;
        }
    }
}

void save_data_from_array(struct oxygen_coord_info* array, FILE* output_file)
{
    for(int i = 0; i < MAX_COORDINATED_METALS; i++)
    {
        if(array[i].metal_ion_index != BLANK)
            fprintf(output_file, "%d\n", array[i].time);
    }
}

void save_last_step_data(struct oxygen_coord_info** carbonate_data, struct oxygen_coord_info*** tfsi_data, struct system_info* system_info, FILE* carbonates_output, FILE* tfsi_output)
{
    for(int i = 0; i < system_info->carbonate_molecules_number; i++)
    {
        save_data_from_array(carbonate_data[i], carbonates_output);
    }

    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        for(int j = 0; j < TRACKED_O_ATOMS_TFSI; j++)
        {
            save_data_from_array(tfsi_data[i][j], tfsi_output);
        }
    }
}

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

void calculate_coord_times_tfsi(struct system_info* system_info, struct tfsi_data* tfsi_data)
{
    for(int i = 0; i < system_info->tfsi_molecules_number; i++)
    {
        calculate_coord_times(TRACKED_O_ATOMS_TFSI, tfsi_data->current_tfsi_coordination[i], tfsi_data->last_tfsi_coordination[i], tfsi_data->tfsi_coordination_info[i], tfsi_data->tfsi_output);
    }
}

void calculate_coord_times_carbonates(struct system_info* system_info, struct carbonate_data* carbonate_data)
{
    calculate_coord_times(system_info->carbonate_molecules_number, carbonate_data->current_carbonate_coordination, carbonate_data->last_carbonate_coordination,
                          carbonate_data->carbonate_coordination_info, carbonate_data->carbonate_output);
}

struct metal_coord_info get_coordination_info(int metal_index, struct vector metal_position, struct coordination_input* coordination_input)
{
    struct metal_coord_info result;
    result.carbonate_molecules = 0;
    result.tfsi_molecules = 0;
    result.tfsi_oxygens = 0;
    result.coordination_number = 0;

    struct program_configuration* program_configuration = coordination_input->program_configuration;
    int carbonate_history_index = 0;
    int tfsi_atoms_history_index = 0;
    int tfsi_history_index = 0;
    int carbonate_oxygen_index = 0;

    int compound_index = 0;
    for(int i = 0; i < coordination_input->system_info->carbonate_types_number; i++)
    {
        compound_index = get_next_carbonate_index(compound_index, coordination_input->system_info->compounds, coordination_input->system_info->compounds_number);
        for(int j = 0; j < coordination_input->system_info->compounds[compound_index].quantity; j++)
        {
            if(detect_coordination(metal_position, coordination_input->carbonate_oxygens[i][j], program_configuration->carbonate_threshold, program_configuration->box_size) == 1)
            {
                result.coordination_number++;
                result.carbonate_molecules++;

                if(program_configuration->calculate_carbonate_residence == 1)
                {
                    coordination_input->carbonate_coordination_histories[i][coordination_input->step_number][metal_index][carbonate_history_index] = j;
                    carbonate_history_index++;
                }
                
                mark_coordination(metal_index, coordination_input->current_carbon_coordination[carbonate_oxygen_index]);
            }
            
            carbonate_oxygen_index++;
        }

        carbonate_history_index = 0;
    }

    for(int i = 0; i < coordination_input->system_info->tfsi_molecules_number; i++)
    {
        int is_metal_coordinated_by_tfsi = 0;
        for(int j = 0; j < TRACKED_O_ATOMS_TFSI; j++)
        {
            if(detect_coordination(metal_position, coordination_input->tfsi_oxygens[i][j], program_configuration->tfsi_threshold, program_configuration->box_size) == 1)
            {
                result.coordination_number++;
                result.tfsi_oxygens++;
                is_metal_coordinated_by_tfsi = 1;

                if(program_configuration->calculate_tfsi_residence == 1)
                {
                    coordination_input->tfsi_atoms_coordination_history[coordination_input->step_number][metal_index][tfsi_atoms_history_index] = i * TRACKED_O_ATOMS_TFSI + j;
                    tfsi_atoms_history_index++;
                }

                mark_coordination(metal_index, coordination_input->current_tfsi_coordination[i][j]);
            }
        }

        if(is_metal_coordinated_by_tfsi == 1)
        {
            result.tfsi_molecules++;

            if(program_configuration->calculate_tfsi_residence == 1)
            {
                coordination_input->tfsi_coordination_history[coordination_input->step_number][metal_index][tfsi_history_index] = i;
                tfsi_history_index++;
            }
        }
    }

    return result;
}