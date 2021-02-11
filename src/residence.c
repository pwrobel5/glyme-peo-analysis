#include <stdlib.h>

#include "program.h"

#define MAX_COORDINATION_NUMBER 10
/*
void delete_history_array(short int*** array, int step_number, int metal_ions)
{
    for(int i = 0; i < step_number; i++)
    {
        for(int j = 0; j < metal_ions; j++)
        {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

short int*** initialize_history_array(short int*** array, int step_number, struct system_info* sys_info)
{
    int current_size = step_number + 1;

    short int*** tmp = realloc(array, current_size * sizeof(short int**));
    if(tmp == NULL)
    {
        delete_history_array(array, step_number, sys_info->metal_ions_number);
        raise_error("Error with reallocation of heavisides array");
    }     
            
    tmp[step_number] = malloc(sys_info->metal_ions_number * sizeof(short int*));
    if(tmp[step_number] == NULL)
    {
        delete_history_array(tmp, step_number, sys_info->metal_ions_number);
        raise_error("Error with allocation of heavisides array for next step");
    }
        
    for(int i = 0; i < sys_info->metal_ions_number; i++)
    {
        tmp[step_number][i] = malloc(MAX_COORDINATION_NUMBER * sizeof(short int));
        if(tmp[step_number][i] == NULL)
        {
            delete_history_array(tmp, step_number, sys_info->metal_ions_number);
            raise_error("Error with allocation of heavisides array for current metal ion");
        }

        for(int j = 0; j < MAX_COORDINATION_NUMBER; j++)
        {
            tmp[step_number][i][j] = BLANK;
        }
    }

    return tmp;
}

double* calculate_residence_times(short*** history_array, int steps, int denominator, struct system_info* system_info)
{
    double* residence = malloc(steps * sizeof(double));
    if(residence == NULL)
        raise_error("Error with memory allocation for residence times array");
    
    for(int i = 0; i < steps; i++)
        residence[i] = 0.0;
    
    for(int i = 0; i < steps; i++)
    {
        residence[0]++;
        
        double zeros_av = 0.0;
        for(int met_ion = 0; met_ion < system_info->metal_ions_number; met_ion++)
        {
            int j = 0;
            while(history_array[i][met_ion][j] != BLANK)
            {
                zeros_av++;
                j++;
            }
        }
        zeros_av /= (double) denominator;

        for(int j = i + 1; j < steps; j++)
        {
            double current_av = 0.0;
            for(int met_ion = 0; met_ion < system_info->metal_ions_number; met_ion++)
            {
                int zero_index = 0;
                int current_index = 0;

                while(current_index < MAX_COORDINATION_NUMBER && zero_index < MAX_COORDINATION_NUMBER && history_array[i][met_ion][zero_index] != BLANK && history_array[j][met_ion][current_index] != BLANK)
                {
                    while(history_array[i][met_ion][zero_index] < history_array[j][met_ion][current_index] && history_array[i][met_ion][zero_index] != BLANK)
                        zero_index++;
                    
                    while(history_array[j][met_ion][current_index] < history_array[i][met_ion][zero_index] && history_array[j][met_ion][current_index] != BLANK)
                        current_index++;
                    
                    if(history_array[i][met_ion][zero_index] == history_array[j][met_ion][current_index])
                        current_av += 1.0;
                    
                    zero_index++;
                }                
            }
            current_av /= (double) denominator;
            current_av /= zeros_av;

            if(current_av != current_av) // check for NaN
                current_av = 0.0;
            residence[j - i] += current_av;
        }
    }

    for(int i = 0; i < steps; i++)
    {
        residence[i] /= (double) (steps - i);
    }

    return residence;
}

void save_residence_to_file(double* residence, const char* residence_file_name, int steps)
{
    FILE* residence_output = fopen(residence_file_name, "w");
    if(residence_output == NULL)
        raise_error("Error with creating output file for residence times");
    
    for(int i = 0; i < steps; i++)
    {
        fprintf(residence_output, "%d %f\n", i, residence[i]);
    }

    fclose(residence_output);
}*/