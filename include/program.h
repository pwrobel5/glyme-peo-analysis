#include <stdio.h>

#ifndef CARBONATE_ANALYSER_H
#define CARBONATE_ANALYSER_H

#define MAX_LINE_LENGTH 150
#define MOLECULE_NAME_LENGTH 20
#define OUTPUT_FILE_NAME_LENGTH 60
#define ATOM_SYMBOL_LENGTH 4

#define TRACKED_O_ATOMS_TFSI 4

#define OMITTED_ATOMS_CARBONATE 10
#define MAX_COORDINATED_METALS 4
#define MAX_COORDINATION_NUMBER 10
#define BLANK -1

#define SEPARATOR " "

enum molecule_type {
    met, ec, f1ec, f2ec, tfsi, other
};

enum print_mode {
    separate, one_output
};

struct program_configuration {
    char* input_file_name;
    char* system_file_name;
    double carbonate_threshold;
    double tfsi_threshold;
    double box_size;
    enum print_mode print_mode;
    int calculate_carbonate_residence;
    int calculate_tfsi_residence;
};

struct system_compound {
    enum molecule_type molecule_type;
    int quantity;
    int atoms_number;
    char* first_atom_symbol;
};

struct system_info {
    int compounds_number;
    int metal_ions_number;
    int carbonate_molecules_number;
    int carbonate_types_number;
    int tfsi_molecules_number;
    int atoms_number;
    struct system_compound* compounds;
};

struct vector {
    double x;
    double y;
    double z;
};

struct oxygen_coord_info {
    int metal_ion_index;
    int time;
};

struct metal_coord_info {
    int coordination_number;
    int carbonate_molecules;
    int tfsi_oxygens;
    int tfsi_molecules;
};

struct carbonate_data {
    int** current_carbonate_coordination;
    int** last_carbonate_coordination;
    struct oxygen_coord_info** carbonate_coordination_info;
    FILE* carbonate_output;
};

struct tfsi_data {
    int*** current_tfsi_coordination;
    int*** last_tfsi_coordination;
    struct oxygen_coord_info*** tfsi_coordination_info;
    FILE* tfsi_output;
};

struct coordination_input {
    int step_number;
    struct vector** carbonate_oxygens;
    struct vector** tfsi_oxygens;
    int** current_carbon_coordination;
    int*** current_tfsi_coordination;
    short int**** carbonate_coordination_histories;
    short int*** tfsi_atoms_coordination_history;
    short int*** tfsi_coordination_history;
    struct system_info* system_info; 
    struct program_configuration* program_configuration;
};

void to_lower_case(char* input_text);
void format_atom_symbol(char* input_text);
void raise_error(const char* message);
char* molecule_type_to_str(enum molecule_type molecule_type);
struct system_info* get_system_info(const char* system_file_name);
void free_system_info(struct system_info* system_info);
int get_next_carbonate_index(int current_index, struct system_compound* compounds, int compounds_number);
void read_data(struct program_configuration* program_configuration, struct system_info* system_info);
void check_symbol(char* read, char* expected);

void read_vector_coordinates(char* line, char* expected_symbol, struct vector* vector);
double vector_norm(struct vector v);
double calculate_distance(struct vector* v1, struct vector* v2, double box_size);

struct metal_coord_info get_coordination_info(int metal_index, struct vector metal_position, struct coordination_input* coordination_input);
void calculate_coord_times_carbonates(struct system_info* system_info, struct carbonate_data* carbonate_data);
void calculate_coord_times_tfsi(struct system_info* system_info, struct tfsi_data* tfsi_data);
void swap_coordination_arrays(int*** first, int*** second);
void clear_coordination_array(int** array, int first_index_max);
void save_last_step_data(struct oxygen_coord_info** carbonate_data, struct oxygen_coord_info*** tfsi_data, struct system_info* system_info, FILE* carbonates_output, FILE* tfsi_output);

short int*** initialize_history_array(short int*** array, int step_number, struct system_info* sys_info);
void delete_history_array(short int*** array, int step_number, int metal_ions);
double* calculate_residence_times(short*** history_array, int steps, int denominator, struct system_info* system_info);
void save_residence_to_file(double* residence, const char* residence_file_name, int steps);

#endif