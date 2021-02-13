#include <stdio.h>

#ifndef MD_ANALYSER_H
#define MD_ANALYSER_H

#define MAX_LINE_LENGTH 150
#define MAX_COORDINATED_CATIONS 4
#define BLANK -1
#define SEPARATOR " "

enum entry_type {
    cation, solvent, anion, other
};

enum print_mode {
    separate, one_output
};

struct program_configuration {
    char* input_file_name;
    char* system_file_name;
    double solvent_threshold;
    double anion_threshold;
    double box_size;
    enum print_mode print_mode;
    int calculate_solvent_residence;
    int calculate_anion_residence;
};

struct system_compound {
    enum entry_type entry_type;
    char* compound_name;
    char* first_atom_symbol;
    int quantity;
    int atoms_number;
    char* tracked_atom_symbol;
    int tracked_atoms_number;
};

struct system_info {
    int compounds_number;
    int cations_number;
    int cation_types_number;
    int solvent_molecules_number;
    int solvent_types_number;
    int anions_number;
    int anion_types_number;
    int atoms_number;
    struct system_compound* compounds;
};

struct vector {
    double x;
    double y;
    double z;
};

struct ligand_coord_info {
    int cation_index;
    int time;
};

struct cation_coord_info {
    int coordination_number;
    int solvent_molecules;
    int anion_atoms;
    int anion_molecules;
};

struct entry_data {
    int*** current_coordination;
    int*** last_coordination;
    struct ligand_coord_info*** coordination_info;
    FILE* output;
};

struct coordination_input {
    int step_number;
    struct vector*** solvent_tracked_atoms;
    struct vector** anion_tracked_atoms;
    int*** current_solvent_coordination;
    int*** current_anion_coordination;
    short int**** solvent_coordination_history;
    short int*** anion_atoms_coordination_history;
    short int*** anion_coordination_history;
    struct system_info* system_info; 
    struct program_configuration* program_configuration;
};

// coordination.c
void swap_coordination_arrays(int*** first, int*** second);
void clear_coordination_array(int** array, int first_index_max, int second_index_max);
void save_last_step_data(struct ligand_coord_info** solvent_data, struct ligand_coord_info*** anion_data, struct system_info* system_info, FILE* solvent_output, FILE* anion_output);
void calculate_coord_times_solvent(struct system_info* system_info, struct entry_data* solvent_data);
void calculate_coord_times_anion(struct system_info* system_info, struct entry_data* anion_data);
struct cation_coord_info get_coordination_info(int cation_index, int cation_tracked_positions_number, struct vector* cation_tracked_positions, struct coordination_input* coordination_input);

// io.c
void print_usage(FILE* stream, int exit_code);
struct program_configuration* read_configuration(int argc, char* argv[]);

// misc.c
void to_lower_case(char* input_text);
void format_atom_symbol(char* input_text);
void raise_error(const char* message);

// reading_utils.c
void check_symbol(const char* read, const char* expected);
int get_next_entry_index(int current_index, struct system_compound* compounds, int compounds_number, enum entry_type entry_type);
void read_data(struct program_configuration* program_configuration, struct system_info* system_info);

// residence.c
void delete_history_array(short int*** array, int step_number, int metal_ions);
short int*** initialize_history_array(short int*** array, int step_number, struct system_info* sys_info);
double* calculate_residence_times(short*** history_array, int steps, int denominator, struct system_info* system_info);
void save_residence_to_file(double* residence, const char* residence_file_name, int steps);

// system_info.c
enum entry_type str_to_entry_type(const char* input_text);
char* entry_type_to_str(enum entry_type molecule_type);
struct system_info* get_system_info(const char* system_file_name);
void free_system_info(struct system_info* system_info);

// vector.c
void read_vector_coordinates(char* line, char* expected_symbol, struct vector* vector);
double vector_norm(struct vector v);
double calculate_distance(struct vector* v1, struct vector* v2, double box_size);

#endif