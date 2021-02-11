#include <stdlib.h>

#include "program.h"

int main(int argc, char* argv[])
{
    struct program_configuration* program_configuration = read_configuration(argc, argv);
    struct system_info* system_info = get_system_info(program_configuration->system_file_name);
    //read_data(&program_configuration, system_info);
    free_system_info(system_info);
    free(program_configuration);

    return 0;
}