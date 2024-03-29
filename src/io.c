#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>

#include "program.h"

#define DEFAULT_SOLVENT_THRESHOLD 1.50
#define DEFAULT_ANION_THRESHOLD 1.50
#define DEFAULT_BOX_SIZE 22.00

void print_usage(FILE *stream, int exit_code)
{
    fprintf(stream, "Usage: [ options ] inputfile systeminfofile\n");
    fprintf(stream,
            "\t-h --help\t\t\tdisplay help text\n"
            "\t-s --solvent-threshold\t\tdeclare threshold for solvent-metal coordination, default: %f\n"
            "\t-a --anion-threshold\t\tdeclare threshold for anion-metal coordination, default: %f\n"
            "\t-b --box-size\t\t\tdeclare box size, default: %f\n"
            "\t-o --one-output\t\t\tproduce one output file for coordination\n"
            "\t-r --solvent-residence\t\tcalculate residence times for solvent molecules\n"
            "\t-f --anion-residence\t\tcalculate residence times for anions\n"
            "\t-v --venn\t\t\tcalculate data for Venn diagrams\n"
            "\t-d --solvent-data\t\tcreate files with numbers of cations coordinated by each solvent molecule and number of coordinating atoms\n"
            "\t-p --box-size-file\t\tfile name with box sizes for each step in NpT ensemble\n\n",
            DEFAULT_SOLVENT_THRESHOLD,
            DEFAULT_ANION_THRESHOLD,
            DEFAULT_BOX_SIZE);
    exit(exit_code);
}

program_configuration_t *read_configuration(int argc, char *argv[])
{
    int next_option;
    const char *const short_options = "hs:a:b:orfvdp:";
    const struct option long_options[] = {
        { "help",                0, NULL, 'h' },
        { "solvent-threshold",   1, NULL, 's' },
        { "anion-threshold",     1, NULL, 'a' },
        { "box-size",            1, NULL, 'b' },
        { "one-output",          0, NULL, 'o' },
        { "solvent-residence",   0, NULL, 'r' },
        { "anion-residence",     0, NULL, 'f' },
        { "venn",                0, NULL, 'v' },
        { "solvent-data",        0, NULL, 'd' },
        { "box-size-file",       1, NULL, 'p' },
        { NULL,                  0, NULL,  0  }
    };

    program_configuration_t *program_configuration = malloc(sizeof(program_configuration_t));
    if (program_configuration == NULL) raise_error("Error with memory allocation for program_configuration!");

    program_configuration->input_file_name = NULL;
    program_configuration->system_file_name = NULL;
    program_configuration->solvent_threshold = DEFAULT_SOLVENT_THRESHOLD;
    program_configuration->anion_threshold = DEFAULT_ANION_THRESHOLD;
    program_configuration->box_size = DEFAULT_BOX_SIZE;
    program_configuration->print_mode = separate;
    program_configuration->calculate_solvent_residence = 0;
    program_configuration->calculate_anion_residence = 0;
    program_configuration->calculate_venn_diagrams = 0;
    program_configuration->save_additional_solvent_data = 0;
    program_configuration->ensemble = NVT;
    program_configuration->box_sizes_file_name = NULL;

    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);
        switch (next_option) {
            case 'h':
                print_usage(stdout, EXIT_SUCCESS);
                break;
            case 's':
                program_configuration->solvent_threshold = atof(optarg);
                break;
            case 'a':
                program_configuration->anion_threshold = atof(optarg);
                break;
            case 'b':
                program_configuration->box_size = atof(optarg);
                break;
            case 'o':
                program_configuration->print_mode = one_output;
                break;
            case 'r':
                program_configuration->calculate_solvent_residence = 1;
                break;
            case 'f':
                program_configuration->calculate_anion_residence = 1;
                break;
            case 'v':
                program_configuration->calculate_venn_diagrams = 1;
                break;
            case 'd':
                program_configuration->save_additional_solvent_data = 1;
                break;
            case 'p':
                program_configuration->ensemble = NpT;
                program_configuration->box_sizes_file_name = optarg;
                break;
            case '?':
                print_usage(stderr, EXIT_FAILURE);
        }
    } while (next_option != -1);

    if (argc > optind + 1) {
        program_configuration->input_file_name = argv[optind];
        program_configuration->system_file_name = argv[optind + 1];
    } else {
        errno = ENODATA;
        perror("No input file given!");
        print_usage(stderr, EXIT_FAILURE);
    }

    return program_configuration;
}