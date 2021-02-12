#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "program.h"
/*
void read_vector_coordinates(char* line, char* expected_symbol, struct vector* vector)
{
    char* tmp;
    tmp = strtok(line, SEPARATOR);

    check_symbol(tmp, expected_symbol);

    vector->x = atof(strtok(NULL, SEPARATOR));
    vector->y = atof(strtok(NULL, SEPARATOR));
    vector->z = atof(strtok(NULL, SEPARATOR));
}
*/

/*
double scalar_product(struct vector v1, struct vector v2)
{
    return ((v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z));
}
*/

/*
double vector_norm(struct vector v)
{
    return sqrt(scalar_product(v, v));
}
*/

/*
double include_pbc(double vector_component, double cell_vector_component)
{
    double result = vector_component;

    if(vector_component > 0.5 * cell_vector_component)
        result -= cell_vector_component;
    else if(vector_component < -0.5 * cell_vector_component)
        result += cell_vector_component;
    
    return result;
}
*/

/*
double calculate_distance(struct vector* v1, struct vector* v2, double box_size)
{
    struct vector v1minusv2;

    v1minusv2.x = include_pbc(v1->x - v2->x, box_size);
    v1minusv2.y = include_pbc(v1->y - v2->y, box_size);
    v1minusv2.z = include_pbc(v1->z - v2->z, box_size);
    
    return sqrt(pow(v1minusv2.x, 2) + pow(v1minusv2.y, 2) + pow(v1minusv2.z, 2));
}*/