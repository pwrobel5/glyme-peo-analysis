#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "program.h"

void save_index_set(int* stack, int position, struct index_set* index_set, int index_set_position)
{
    index_set[index_set_position].size = position;
    index_set[index_set_position].indices = malloc(position * sizeof(int));
    if(index_set[index_set_position].indices == NULL) raise_error("Error with memory allocation for index set");

    for(int i = 1; i <= position; i++)
    {
        index_set[index_set_position].indices[i - 1] = stack[i] - 1; // -1 to have indices starting from 0 not 1
    }
}

void determine_powerset(int tracked_atoms_number, struct index_combinations* combination_set)
{
    int stack[tracked_atoms_number + 1]; // +1 to have blank first position
    stack[0] = 0;
    int position = 0;
    int combination_index = 0;

    do {
        if(stack[position] < tracked_atoms_number)
        {
            stack[position + 1] = stack[position] + 1;
            position++;
        }
        else
        {
            stack[position - 1]++;
            position--;
        }

        if(position > 1) // 1 to avoid including one-element sets
        {
            save_index_set(stack, position, combination_set->combinations, combination_index);
            combination_index++;
        }

    } while(position > 0);
}

struct index_combinations* get_index_combinations(int tracked_atoms_number)
{
    int powerset_size = (int) pow(2, tracked_atoms_number) - 1 - tracked_atoms_number; // -1 because we don't need empty set, -tracked_atoms_number to not have one-element sets
    if(powerset_size < 0) raise_error("Too many tracked atoms to determine Venn diagrams");

    struct index_combinations* result = malloc(sizeof(struct index_combinations));
    if(result == NULL) raise_error("Error with memory allocation for index combinations");

    result->combinations_number = powerset_size;
    result->combinations = malloc(powerset_size * sizeof(struct index_set));
    if(result->combinations == NULL) raise_error("Error with memory allocation for index combinations");
    determine_powerset(tracked_atoms_number, result);

    return result;
}

void free_index_combinations(struct index_combinations* index_combinations)
{
    for(int i = 0; i < index_combinations->combinations_number; i++)
    {
        free(index_combinations->combinations[i].indices);
    }

    free(index_combinations->combinations);
    free(index_combinations);
}

void add_element_to_result_set(struct index_set* result, int index, int element)
{
    result->size += 1;
    result->indices = realloc(result->indices, result->size * sizeof(int));
    if(result->indices == NULL) raise_error("Error with memory reallocation for result set");
    result->indices[index] = element;
}

struct index_set* set_intersection(struct index_set* A, struct index_set* B)
{
    struct index_set* result = malloc(sizeof(struct index_set));
    if(result == NULL) raise_error("Error with memory allocation for set intersection");
    int index_A = 0;
    int index_B = 0;
    int index_result = 0;

    result->size = 0;
    result->indices = NULL;

    while(index_A < A->size && index_B < B->size)
    {
        if(A->indices[index_A] < B->indices[index_B])
            index_A++;
        else if(A->indices[index_A] > B->indices[index_B])
            index_B++;
        else
        {
            add_element_to_result_set(result, index_result, A->indices[index_A]);
            index_result++;
            index_A++;
            index_B++;
        }
    }

    return result;
}

struct index_set* set_difference(struct index_set* A, struct index_set* B)
{
    struct index_set* result = malloc(sizeof(struct index_set));
    if(result == NULL) raise_error("Error with memory allocation for set difference");
    int index_A = 0;
    int index_B = 0;
    int index_result = 0;

    result->size = 0;
    result->indices = NULL;

    while(index_A < A->size && index_B < B->size)
    {
        if(A->indices[index_A] < B->indices[index_B])
        {
            add_element_to_result_set(result, index_result, A->indices[index_A]);
            index_result++;
            index_A++;
        }
        else if(A->indices[index_A] > B->indices[index_B])
        {
            index_B++;
        }
        else
        {
            index_A++;
            index_B++;
        }
    }

    while(index_A < A->size)
    {
        add_element_to_result_set(result, index_result, A->indices[index_A]);
        index_result++;
        index_A++;
    }

    return result;
}

void free_set(struct index_set* set)
{
    if(set->size > 0)
        free(set->indices);
    free(set);
}

struct index_set* make_set(int array_size, int* array)
{
    int actual_size = 0;
    for(int i = 0; i < array_size && array[i] != BLANK; i++)
        actual_size++;
    
    struct index_set* result = malloc(sizeof(struct index_set));
    if(result == NULL) raise_error("Error with memory allocation for created set");

    result->size = actual_size;
    result->indices = NULL;

    if(actual_size > 0) {
        result->indices = malloc(actual_size * sizeof(int));
        if(result->indices == NULL) raise_error("Error with memory allocation for new set indices");
    }

    for(int i = 0; i < actual_size; i++)
    {
        result->indices[i] = array[i];
    }

    return result;
}

void determine_venn_sets(struct venn_set* venn_set, int** current_coordination, int tracked_atoms_number, struct index_combinations* index_combinations)
{
    venn_set->entries_number = index_combinations->combinations_number + tracked_atoms_number; // index_combinations do not include one-element sets
    venn_set->entries = malloc(venn_set->entries_number * sizeof(struct venn_entry));
    if(venn_set->entries == NULL) raise_error("Error with memory allocation for Venn entry");

    for(int i = 0; i < tracked_atoms_number; i++)
    {
        venn_set->entries[i].entry_id = malloc(sizeof(int));
        if(venn_set->entries[i].entry_id == NULL) raise_error("Error with entry id allocation for Venn entry");
        venn_set->entries[i].entry_id[0] = i;
        venn_set->entries[i].entry_id_size = 1;
        venn_set->entries[i].set = make_set(MAX_COORDINATED_CATIONS, current_coordination[i]);
    }

    for(int i = 0; i < index_combinations->combinations_number; i++)
    {

        struct index_set current_combination = index_combinations->combinations[i];
        struct index_set* first_set = venn_set->entries[current_combination.indices[0]].set;
        struct index_set* second_set = venn_set->entries[current_combination.indices[1]].set; // combinations contain at least 2 elements
        struct index_set* tmp = NULL;

        first_set = set_intersection(first_set, second_set);

        for(int j = 2; j < current_combination.size; j++)
        {
            if(tmp != NULL) free_set(tmp);
            second_set = venn_set->entries[current_combination.indices[j]].set;

            tmp = first_set;
            first_set = set_intersection(first_set, second_set);
        }

        if(tmp != NULL) free_set(tmp);

        int entry_index = tracked_atoms_number + i;
        venn_set->entries[entry_index].entry_id = malloc(current_combination.size * sizeof(int));
        if(venn_set->entries[entry_index].entry_id == NULL) raise_error("Error with entry id allocation for Venn entry");
        memcpy(venn_set->entries[entry_index].entry_id, current_combination.indices, current_combination.size * sizeof(int));
        venn_set->entries[entry_index].entry_id_size = current_combination.size;
        venn_set->entries[entry_index].set = first_set;
    }
}

void print_venn_set(struct venn_set* venn_set)
{
    for(int i = 0; i < venn_set->entries_number; i++)
    {        
        struct venn_entry current_entry = venn_set->entries[i];

        if(current_entry.set->size > 0)
        {
            printf("Entry index: %d\n", i);
            printf("Entry id: ");
            for(int j = 0; j < current_entry.entry_id_size; j++)
            {
                printf("%d ", current_entry.entry_id[j]);
            }
            printf("\n");

            struct index_set* set = current_entry.set;
            printf("Set indices: ");
            for(int j = 0; j < set->size; j++)
            {
                printf("%d ", set->indices[j]);
            }
            printf("\n");
        }
    }
}

void free_venn_set(struct venn_set* venn_set)
{
    for(int i = 0; i < venn_set->entries_number; i++)
    {
        free(venn_set->entries[i].entry_id);
        free_set(venn_set->entries[i].set);
    }

    free(venn_set->entries);
}

struct venn_diagram* create_empty_venn_diagram(int tracked_atoms_number, struct index_combinations* index_combinations)
{
    struct venn_diagram* result = malloc(sizeof(struct venn_diagram));
    if(result == NULL) raise_error("Error with memory allocation for Venn diagram");

    result->entries_number = tracked_atoms_number + index_combinations->combinations_number;
    result->entries = malloc(result->entries_number * sizeof(struct venn_diagram_entry));
    if(result->entries == NULL) raise_error("Error with memory allocation for Venn diagram");

    for(int i = 0; i < tracked_atoms_number; i++)
    {
        result->entries[i].set_size = 0;
        result->entries[i].entry_id_size = 1;
        result->entries[i].entry_id = malloc(sizeof(int));
        if(result->entries[i].entry_id == NULL) raise_error("Error with memory allocation for Venn diagram");
        result->entries[i].entry_id[0] = i;
    }

    for(int i = 0; i < index_combinations->combinations_number; i++)
    {
        int entry_index = i + tracked_atoms_number;
        struct index_set current_combination = index_combinations->combinations[i];

        result->entries[entry_index].set_size = 0;
        result->entries[entry_index].entry_id_size = current_combination.size;
        result->entries[entry_index].entry_id = malloc(current_combination.size * sizeof(int));
        if(result->entries[entry_index].entry_id == NULL) raise_error("Error with memory allocation for Venn diagram");
        memcpy(result->entries[entry_index].entry_id, current_combination.indices, current_combination.size * sizeof(int));
    }

    return result;
}

struct venn_diagram* determine_venn_diagram(struct venn_set* venn_set, int tracked_atoms_number)
{
    struct venn_diagram* result = malloc(sizeof(struct venn_diagram));
    if(result == NULL) raise_error("Error with memory allocation for Venn diagram");

    result->entries_number = venn_set->entries_number;
    result->entries = malloc(result->entries_number * sizeof(struct venn_diagram_entry));
    if(result->entries == NULL) raise_error("Error with memory allocation for Venn diagram");

    for(int i = 0; i < venn_set->entries_number; i++)
    {
        struct venn_entry current_entry = venn_set->entries[i];
        int entry_id_index = 0;
        int tracked_atom_index = 0;
        struct index_set* differences[tracked_atoms_number - current_entry.entry_id_size + 1]; // + 1 to have non zero size array and keep original set at 0 position
        differences[0] = current_entry.set;
        int differences_index = 1;

        while(tracked_atom_index < tracked_atoms_number)
        {      
            if(entry_id_index >= current_entry.entry_id_size || (entry_id_index < current_entry.entry_id_size && tracked_atom_index < current_entry.entry_id[entry_id_index]))
            {
                struct index_set* set_to_remove = venn_set->entries[tracked_atom_index].set;
                struct index_set* difference = set_difference(differences[differences_index - 1], set_to_remove);
                differences[differences_index] = difference;
                differences_index++;
                tracked_atom_index++;
            }
            else
            {
                entry_id_index++;
                tracked_atom_index++;
            }
        }

        result->entries[i].entry_id_size = current_entry.entry_id_size;
        result->entries[i].entry_id = malloc(result->entries[i].entry_id_size * sizeof(int));
        if(result->entries[i].entry_id == NULL) raise_error("Error with memory allocation for Venn diagram");
        memcpy(result->entries[i].entry_id, current_entry.entry_id, result->entries[i].entry_id_size * sizeof(int));
        result->entries[i].set_size = differences[differences_index - 1]->size;

        for(int j = 1; j < differences_index; j++)
        {
            free_set(differences[j]);
        }
    }

    return result;
}

void print_venn_diagram(struct venn_diagram* venn_diagram)
{
    for(int i = 0; i < venn_diagram->entries_number; i++)
    {        
        struct venn_diagram_entry current_entry = venn_diagram->entries[i];
        if(current_entry.set_size > 0)
        {
            printf("Entry index: %d\n", i);
            printf("Entry id: ");
            for(int j = 0; j < current_entry.entry_id_size; j++)
            {
                printf("%d ", current_entry.entry_id[j]);
            }
            printf("\n");
            printf("Counter: %d\n\n", current_entry.set_size);
        }
    }
}

void free_venn_diagram(struct venn_diagram* venn_diagram)
{
    for(int i = 0; i < venn_diagram->entries_number; i++)
    {
        free(venn_diagram->entries[i].entry_id);
    }

    free(venn_diagram->entries);
    free(venn_diagram);
}

void update_global_venn_diagram(struct venn_diagram* global, struct venn_diagram* local)
{
    for(int i = 0; i < local->entries_number; i++)
    {
        if(memcmp(global->entries[i].entry_id, local->entries[i].entry_id, global->entries[i].entry_id_size * sizeof(int)) != 0)
            printf("[WARNING] Different order of indices in Venn diagram structures");
        
        global->entries[i].set_size += local->entries[i].set_size;
    }
}

void save_averages_to_file(struct venn_diagram** venn_diagrams, int step_number, char* output_file_name)
{
    if(step_number <= 0) return;

    FILE* output_file = fopen(output_file_name, "w");
    int entries_number = venn_diagrams[0]->entries_number;

    for(int i = 0; i < entries_number; i++)
    {
        int sum = 0;
        int entry_id_size = venn_diagrams[0]->entries[i].entry_id_size;
        int* entry_id = venn_diagrams[0]->entries[i].entry_id;        

        for(int j = 0; j < step_number; j++)
        {
            sum += venn_diagrams[j]->entries[i].set_size;
        }

        double average = sum / (double) step_number;
        for(int j = 0; j < entry_id_size; j++)
        {
            fprintf(output_file, "%d ", entry_id[j]);
        }
        fprintf(output_file, ": %f\n", average);
    }
 
    fclose(output_file);
}