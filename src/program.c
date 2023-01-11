#include <stdlib.h>

#include "program.h"

int main(int argc, char *argv[])
{
    program_configuration_t *program_configuration = read_configuration(argc, argv);
    system_info_t *system_info = get_system_info(program_configuration->system_file_name);
    read_data(program_configuration, system_info);
    free_system_info(system_info);
    free(program_configuration);

    return 0;
}