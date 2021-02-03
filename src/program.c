#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>

#include "program.h"

#define DEFAULT_CARBONATE_THRESHOLD 1.50
#define DEFAULT_TFSI_THRESHOLD 1.50
#define DEFAULT_BOX_SIZE 22.00

void print_usage(FILE* stream, int exit_code)
{
    fprintf(stream, "Usage: [ options ] inputfile systeminfofile\n");
    fprintf(stream,
            "\t-h --help\t\t\tdisplay help text\n"
            "\t-c --carbonate-threshold\tdeclare threshold for carbonate-metal coordination, default: %f\n"
            "\t-t --tfsi-threshold\t\tdeclare threshold for TFSI-metal coordination, default: %f\n"
            "\t-b --box-size\t\t\tdeclare box size, default: %f\n"
            "\t-o --one-output\t\t\tproduce one output file for coordination\n"
            "\t-r --carbonate-residence\t\tcalculate residence times for carbonates\n"
            "\t-f --tfsi-residence\t\tcalculate residence times for TFSI anions\n\n",
            DEFAULT_CARBONATE_THRESHOLD,
            DEFAULT_TFSI_THRESHOLD,
            DEFAULT_BOX_SIZE);
    exit(exit_code);
}

void parse_initial_arguments(int argc, char* argv[], struct program_configuration* program_configuration)
{
    int next_option;
    const char* const short_options = "hc:t:b:orf";
    const struct option long_options[] = {
        { "help",                0, NULL, 'h' },
        { "carbonate-threshold", 1, NULL, 'c' },
        { "tfsi-threshold",      1, NULL, 't' },
        { "box-size",            1, NULL, 'b' },
        { "one-output",          0, NULL, 'o' },
        { "carbonate-residence", 0, NULL, 'r' },
        { "tfsi-residence",      0, NULL, 'f' },
        { NULL,                  0, NULL,  0  }
    };

    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);
        switch(next_option)
        {
            case 'h':
                print_usage(stdout, EXIT_SUCCESS);
                break;
            case 'c':
                program_configuration->carbonate_threshold = atof(optarg);
                break;
            case 't':
                program_configuration->tfsi_threshold = atof(optarg);
                break;
            case 'b':
                program_configuration->box_size = atof(optarg);
                break;
            case 'o':
                program_configuration->print_mode = one_output;
                break;
            case 'r':
                program_configuration->calculate_carbonate_residence = 1;
                break;
            case 'f':
                program_configuration->calculate_tfsi_residence = 1;
                break;
            case '?':
                print_usage(stderr, EXIT_FAILURE);
        }
    } while(next_option != -1);

    if(argc > optind + 1)
    {
        program_configuration->input_file_name = argv[optind];
        program_configuration->system_file_name = argv[optind + 1];
    }
    else
    {
        errno = ENODATA;
        perror("No input file given!");
        print_usage(stderr, EXIT_FAILURE);
    }    
}

int main(int argc, char* argv[])
{
    struct program_configuration program_configuration;

    program_configuration.input_file_name = NULL;
    program_configuration.system_file_name = NULL;
    program_configuration.carbonate_threshold = DEFAULT_CARBONATE_THRESHOLD;
    program_configuration.tfsi_threshold = DEFAULT_TFSI_THRESHOLD;
    program_configuration.box_size = DEFAULT_BOX_SIZE;
    program_configuration.print_mode = separate;
    program_configuration.calculate_carbonate_residence = 0;
    program_configuration.calculate_tfsi_residence = 0;

    parse_initial_arguments(argc, argv, &program_configuration);
    struct system_info* system_info = get_system_info(program_configuration.system_file_name);
    read_data(&program_configuration, system_info);
    free_system_info(system_info);

    return 0;
}