#include <stdlib.h>

#include "CUnit/Basic.h"
#include "program.h"

int A_indices[] = {1, 3, 5, 7};
int B_indices[] = {2, 4, 6, 8};
int C_indices[] = {1, 5};

index_set_t A = {
    4,
    A_indices
};
index_set_t B = {
    4,
    B_indices
};
index_set_t C = {
    2,
    C_indices
};

void test_intersection(void)
{
    index_set_t *A_intersect_C = set_intersection(&A, &C);
    int expected_indices[] = {1, 5};
    CU_ASSERT(A_intersect_C->size == 2);
    CU_ASSERT(memcmp(A_intersect_C->indices, expected_indices, 2) == 0);
    free_set(A_intersect_C);
}

void test_empty_intersection(void)
{
    index_set_t *A_intersect_B = set_intersection(&A, &B);
    CU_ASSERT(A_intersect_B->size == 0);
    CU_ASSERT(A_intersect_B->indices == NULL);
    free_set(A_intersect_B);
}

void test_the_same_set_intersection(void)
{
    index_set_t *A_intersect_A = set_intersection(&A, &A);
    CU_ASSERT(A_intersect_A->size == A.size);
    CU_ASSERT(memcmp(A_intersect_A->indices, A.indices, A.size) == 0);
    free_set(A_intersect_A);
}

void test_difference(void)
{
    index_set_t *A_diff_C = set_difference(&A, &C);
    int expected_indices[] = {3, 7};
    CU_ASSERT(A_diff_C->size == 2);
    CU_ASSERT(memcmp(A_diff_C->indices, expected_indices, 2) == 0);
    free_set(A_diff_C);
}

void test_non_overlapping_sets_difference(void)
{
    index_set_t *A_diff_B = set_difference(&A, &B);
    CU_ASSERT(A_diff_B->size == A.size);
    CU_ASSERT(memcmp(A_diff_B->indices, A.indices, A.size) == 0);
    free_set(A_diff_B);
}

void test_the_same_set_difference(void)
{
    index_set_t *A_diff_A = set_difference(&A, &A);
    CU_ASSERT(A_diff_A->size == 0);
    CU_ASSERT(A_diff_A->indices == NULL);
    free_set(A_diff_A);
}

void test_make_set(void)
{
    int array[] = {1, 3, 4, BLANK, BLANK, BLANK};
    int array_size = 6;

    index_set_t *set = make_set(array_size, array);
    int expected_indices[] = {1, 3, 4};
    CU_ASSERT(set->size == 3);
    CU_ASSERT(memcmp(set->indices, expected_indices, 3) == 0);
    free_set(set);
}

void test_make_empty_set(void)
{
    int array[] = {BLANK, BLANK};
    int array_size = 2;

    index_set_t *set = make_set(array_size, array);
    CU_ASSERT(set->size == 0);
    CU_ASSERT(set->indices == NULL);
    free_set(set);
}

void test_group_non_blanks(void)
{
    int array1[] = {BLANK, 1, BLANK, 1, 1, 1, BLANK, BLANK};
    int array2[] = {1, 1, 1, 1};

    int output1[] = {1, 3, 4, 5, BLANK, BLANK, BLANK, BLANK};
    int output2[] = {0, 1, 2, 3};

    group_non_blanks_in_beginning(array1, sizeof array1 / sizeof(int));
    group_non_blanks_in_beginning(array2, sizeof array2 / sizeof(int));

    CU_ASSERT(memcmp(array1, output1, sizeof array1) == 0);
    CU_ASSERT(memcmp(array2, output2, sizeof array2) == 0);
}

int main()
{
    CU_pSuite pSuite = NULL;

    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();
    
    pSuite = CU_add_suite("Set tests", NULL, NULL);
    if (pSuite == NULL) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    if ((CU_add_test(pSuite, "simple intersection test", test_intersection) == NULL) ||
       (CU_add_test(pSuite, "test empty result of intersection", test_empty_intersection) == NULL) ||
       (CU_add_test(pSuite, "test of intersection between the same set", test_the_same_set_intersection) == NULL) ||
       (CU_add_test(pSuite, "simple difference test", test_difference) == NULL) ||
       (CU_add_test(pSuite, "test of non overlapping sets difference", test_non_overlapping_sets_difference) == NULL) ||
       (CU_add_test(pSuite, "test difference between the same set", test_the_same_set_difference) == NULL) ||
       (CU_add_test(pSuite, "test make non-empty set", test_make_set) == NULL) ||
       (CU_add_test(pSuite, "test make empty set", test_make_empty_set) == NULL) ||
       (CU_add_test(pSuite, "test of grouping non-blanks", test_group_non_blanks) == NULL)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}