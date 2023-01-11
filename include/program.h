#include <stdio.h>

#ifndef MD_ANALYSER_H
#define MD_ANALYSER_H

#define MAX_LINE_LENGTH 150
#define MAX_COORDINATED_CATIONS 4
#define BLANK -1
#define SEPARATOR " "

typedef enum statistical_ensemble {
    NVT, NpT
} statistical_ensemble_t;

typedef enum entry_type {
    cation, solvent, anion, other
} entry_type_t;

typedef enum print_mode {
    separate, one_output
} print_mode_t;

typedef struct program_configuration {
    char *input_file_name;
    char *system_file_name;
    double solvent_threshold;
    double anion_threshold;
    double box_size;
    print_mode_t print_mode;
    int calculate_solvent_residence;
    int calculate_anion_residence;
    int calculate_venn_diagrams;
    int save_additional_solvent_data;
    statistical_ensemble_t ensemble;
    char *box_sizes_file_name;
} program_configuration_t;

typedef struct system_compound {
    entry_type_t entry_type;
    char *compound_name;
    char *first_atom_symbol;
    int quantity;
    int atoms_number;
    char *tracked_atom_symbol;
    int tracked_atoms_number;
} system_compound_t;

typedef struct system_info {
    int compounds_number;
    int cations_number;
    int cation_types_number;
    int solvent_molecules_number;
    int solvent_types_number;
    int anions_number;
    int anion_types_number;
    int atoms_number;
    system_compound_t *compounds;
} system_info_t;

typedef struct vector {
    double x;
    double y;
    double z;
} vector_t;

typedef struct ligand_coord_info {
    int cation_index;
    int time;
} ligand_coord_info_t;

typedef struct cation_coord_info {
    int coordination_number;
    int solvent_atoms;
    int solvent_molecules;
    int anion_atoms;
    int anion_molecules;
} cation_coord_info_t;

typedef struct entry_data {
    int ***current_coordination;
    int ***last_coordination;
    ligand_coord_info_t ***coordination_info;
    FILE *output;
} entry_data_t;

typedef struct coordination_input {
    int step_number;
    vector_t ***solvent_tracked_atoms;
    vector_t **anion_tracked_atoms;
    int **coordinated_cations_to_anions;
    int ***current_solvent_coordination;
    int ***current_anion_coordination;
    int ***current_coordinating_solvents;
    short int ****solvent_atoms_coordination_history;
    short int ****solvent_coordination_history;
    short int ****anion_atoms_coordination_history;
    short int ****anion_coordination_history;
    system_info_t *system_info; 
    program_configuration_t *program_configuration;
} coordination_input_t;

typedef struct index_set {
    int size;
    int *indices;
} index_set_t;

typedef struct index_combinations {
    int combinations_number;
    index_set_t *combinations;
} index_combinations_t;

typedef struct venn_entry {
    int entry_id_size;
    int *entry_id;
    index_set_t *set;
} venn_entry_t;

typedef struct venn_set {
    int entries_number;
    venn_entry_t *entries;
} venn_set_t;

typedef struct venn_diagram_entry {
    int entry_id_size;
    int *entry_id;
    int set_size;
} venn_diagram_entry_t;

typedef struct venn_diagram {
    int entries_number;
    venn_diagram_entry_t *entries;
} venn_diagram_t;

// coordination.c
void swap_coordination_arrays(int ***first, int ***second);
void clear_coordination_array(int **array, int first_index_max, int second_index_max);
void save_last_step_data(system_info_t *system_info, entry_data_t *entry_data, entry_type_t entry_type);
void calculate_coord_times(system_info_t *system_info, entry_data_t *entry_data, entry_type_t entry_type);
cation_coord_info_t get_coordination_info(int cation_index, int cation_tracked_positions_number, vector_t *cation_tracked_positions, coordination_input_t *coordination_input);

// io.c
void print_usage(FILE *stream, int exit_code);
program_configuration_t *read_configuration(int argc, char *argv[]);

// misc.c
void to_lower_case(char *input_text);
void format_atom_symbol(char *input_text);
void raise_error(const char *message);

// reading_utils.c
void check_symbol(const char *read, const char *expected);
int get_next_entry_index(int current_index, system_compound_t *compounds, int compounds_number, entry_type_t entry_type);
void group_non_blanks_in_beginning(int *array, int array_size);
void read_data(program_configuration_t *program_configuration, system_info_t *system_info);

// residence.c
void delete_history_array(short int ***array, int step_number, int metal_ions);
short int ***initialize_history_array(short int ***array, int step_number, system_info_t *sys_info);
double *calculate_residence_times(short int ***history_array, int steps, int denominator, system_info_t *system_info);
void save_residence_to_file(double *residence, const char *residence_file_name, int steps);

// set.c
index_combinations_t *get_index_combinations(int tracked_atoms_number);
void free_index_combinations(index_combinations_t *index_combinations);
index_set_t *set_intersection(index_set_t *A, index_set_t *B);
index_set_t *set_difference(index_set_t *A, index_set_t *B);
void free_set(index_set_t *set);
index_set_t *make_set(int array_size, int *array);
void determine_venn_sets(venn_set_t *venn_set, int **current_coordination, int tracked_atoms_number, index_combinations_t *index_combinations, int current_coordination_size);
void print_venn_set(venn_set_t *venn_set);
void free_venn_set(venn_set_t *venn_set);
venn_diagram_t *create_empty_venn_diagram(int tracked_atoms_number,  index_combinations_t *index_combinations);
venn_diagram_t* determine_venn_diagram(venn_set_t *venn_set, int tracked_atoms_number);
void print_venn_diagram(venn_diagram_t *venn_diagram);
void free_venn_diagram(venn_diagram_t *venn_diagram);
void update_global_venn_diagram(venn_diagram_t *global, venn_diagram_t *local);
void save_averages_to_file(venn_diagram_t **venn_diagrams, int step_number, char *output_file_name);

// solvent_data.c
void save_current_step_solvent_data(int step_number, int **current_coordination, int tracked_atoms, int cations_number, FILE *output_file);

// system_info.c
entry_type_t str_to_entry_type(const char *input_text);
char *entry_type_to_str(entry_type_t molecule_type);
system_info_t *get_system_info(const char *system_file_name);
void free_system_info(system_info_t *system_info);

// vector.c
void read_vector_coordinates(char *line, char *expected_symbol, vector_t *vector);
double vector_norm(vector_t v);
double calculate_distance(vector_t *v1, vector_t *v2, double box_size);

#endif