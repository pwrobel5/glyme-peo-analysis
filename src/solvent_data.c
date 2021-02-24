#include "program.h"

void save_current_step_solvent_data(int step_number, int** current_coordination, int tracked_atoms, int cations_number, FILE* output_file)
{
    int distinct_cations = 0;
    int coordinating_atoms = 0;
    short int cations[cations_number];
    for(int i = 0; i < cations_number; i++)
    {
        cations[i] = 0;
    }

    for(int i = 0; i < tracked_atoms; i++)
    {
        int non_blanks = 0;
        int j = 0;
        while(j < MAX_COORDINATED_CATIONS && current_coordination[i][j] != BLANK)
        {
            non_blanks++;
            if(cations[current_coordination[i][j]] == 0)
            {
                cations[current_coordination[i][j]] = 1;
                distinct_cations++;
            }

            j++;
        }
        
        if(non_blanks > 0) coordinating_atoms++;
    }

    fprintf(output_file, "%d %d %d\n", step_number, distinct_cations, coordinating_atoms);
}