#include "stdio.h"
#include <stdlib.h>
#include <math.h>
#define DIM 3

typedef float vec3[DIM];

void verlet_step(float r[864][DIM], float v[864][DIM], int t, int M, float rho);
vec3 *f_optimized(vec3 r_i, vec3 r_j, int M, float rho);

int list_pbc[27][3] = {
        { 0, 0, 0 },{ 0, 0, 1 },{ 0, 0, -1 },{ 0, 1, 0 },{ 0, 1, 1 },{ 0, 1, -1 },
        { 0, -1, 0 },{ 0, -1, 1 },{ 0, -1, -1 },{ 1, 0, 0 },{ 1, 0, 1 },
        { 1, 0, -1 },{ 1, 1, 0 },{ 1, 1, 1 },{ 1, 1, -1 },{ 1, -1, 0 },{ 1, -1, 1 },
        { 1, -1, -1 },{ -1, 0, 0 },{ -1, 0, 1 },{ -1, 0, -1 },
        { -1, 1, 0 },{ -1, 1, 1 },{ -1, 1, -1 },{ -1, -1, 0 },{ -1, -1, 1 },{ -1, -1, -1 }
};

float result_r[1501][864][DIM];
float result_v[1501][864][DIM];

int* arangeWithout(int M, int excludeValue) {
    int size = 4 * M * M * M;

    // Allocate memory for the array
    double *originalArray = (double *)malloc(size * sizeof(double));
    // Create the result array without the specified value
    int resultSize = size - 1;
    double *resultArray = (double *)malloc(resultSize * sizeof(double));

    // Fill the original array with values similar to np.arange(4 * M * M * M)
    for (int i = 0; i < size; i++) {
        originalArray[i] = i;
    }

    // Copy values from originalArray to resultArray, excluding the specified value
    int resultIndex = 0;
    for (int i = 0; i < size; i++) {
        if (originalArray[i] != excludeValue) {
            resultArray[resultIndex] = originalArray[i];
            resultIndex++;
        }
    }

    // Free the memory allocated for the originalArray
    free(originalArray);

    return resultArray;
}

vec3 *F(float r[864][DIM],int i ,int M, float rho) { // Force that acts on particle i from all other particles
    vec3 *arr_F = (vec3 *)malloc(3 * sizeof(float));
    for (int d = 0; d < DIM; ++d) {
        (*arr_F)[d] = 0.0;
    }
    int *indices = arangeWithout(864, i);

    for (int k = 0; k < 864 - 1; ++k) {
        int index = indices[k];
        vec3 *force = (vec3 *) f_optimized(r[i], r[index], M, rho);
        for (int d = 0; d < DIM; ++d) {
            (*arr_F)[d] += (*force)[d];
        }
        free(force);
    }


    free(indices);
    return arr_F;
}

vec3 *f_optimized(vec3 r_i, vec3 r_j, int M, float rho) {
    vec3 *res = (vec3 *)malloc(3 * sizeof(float));
    double L = (4 * pow(M, 3)) / rho;
    int m = 48;
    vec3 vec_L_min;
    float r_ij_candidates[27] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                 0.0,0.0,0.0,0.0};

    for (int i = 0; i < 27; i++) {
        for (int d = 0; d < DIM; d++) {
            float diff = r_i[d] - r_j[d] + L * list_pbc[i][d];
            r_ij_candidates[i] += diff * diff;
        }
    }

    int index_L_min = 0;
    float min_val = r_ij_candidates[0];
    for (int i = 1; i < 27; i++) {
        if (r_ij_candidates[i] < min_val) {
            min_val = r_ij_candidates[i];
            index_L_min = i;
        }
    }

    for (int d = 0; d < DIM; d++) {
        vec_L_min[d] = L * list_pbc[index_L_min][d];
    }

    float r_ij = min_val;

    float f_x = m * (r_i[0] - r_j[0] + vec_L_min[0]) * (pow(r_ij, -14) + 0.5 * pow(r_ij, -8));
    float f_y = m * (r_i[1] - r_j[1] + vec_L_min[1]) * (pow(r_ij, -14) + 0.5 * pow(r_ij, -8));
    float f_z = m * (r_i[2] - r_j[2] + vec_L_min[2]) * (pow(r_ij, -14) + 0.5 * pow(r_ij, -8));

    (*res)[0] = f_x;
    (*res)[1] = f_y;
    (*res)[2] = f_z;
    return res;
}

void verlet_step(float (*r_prior)[3], float (*v_prior)[3], int t, int M, float rho) {
    float h = 0.016;
    float v_tilde[DIM];

    for (int i = 0; i < 864; ++i) {
        vec3 *f_current = F(r_prior, i, M, rho);
        vec3 *f_next = F(result_r[t], i, M, rho);
        for (int d = 0; d < DIM; ++d) {
            v_tilde[d] = v_prior[i][d] + h / (2.0 * 48) * (*f_current[d]);
            result_r[t][i][d] = r_prior[i][d] + h * v_tilde[d];
            result_v[t][i][d] = v_tilde[d] + (h / (2.0 * 48)) * (*f_next[d]);
        }
        free(f_next);
        free(f_current);
    }
}

float* readArrayFromFile(const char* filename, int* size) {
    printf(filename);
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Get the size of the array
    if (fscanf(file, "%d", size) != 1) {
        fprintf(stderr, "Error reading array size from file\n");
        fclose(file);
        return NULL;
    }

    printf(" Size: ");
    printf("%d", size);

    // Allocate memory for the array
    float* array = (float*)malloc(*size * sizeof(float));
    if (array == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        fclose(file);
        return NULL;
    }

    // Read the array elements
    for (int i = 0; i < *size; ++i) {
        if (fscanf(file, "%f", &array[i]) != 1) {
            fprintf(stderr, "Error reading array element from file\n");
            free(array);
            fclose(file);
            return NULL;
        }
    }

    // Close the file
    fclose(file);

    return array;
}

int main(int argc, char *argv[ ]) { //argv: [M, rho, r_0, v_0]

    int size1 = 684*3;
    int size2 = 684*3;
    float *param1 = readArrayFromFile(argv[1], 684*3);
    float *param2 = readArrayFromFile(argv[1], 684*3);

    if (param1 == NULL || param2 == NULL) {
        fprintf(stderr, "Error reading arrays from file\n");
        return EXIT_FAILURE;
    }

    // Your code that uses param1 and param2 goes here

    printf("Array 1: ");
    for (int i = 0; i < size1; ++i) {
        printf("%f ", param1[i]);
    }
    printf("\n");

    printf("Array 2: ");
    for (int i = 0; i < size2; ++i) {
        printf("%f ", param2[i]);
    }
    printf("\n");

    // Free allocated memory
    free(param1);
    free(param2);

    printf("w");
    int M = 6; //atoi(argv[0]);
    float rho = 0.85;//atof(argv[1]);
    printf("h");
    float r_0[864][3];
    float v_0[864][3];

    for (int i = 0; i < 864; ++i) {
        for (int j = 0; j < 3; ++j) {
            r_0[i][j] = atof(param1);
            v_0[i][j] = atof(param2);
            result_r[0][i][j] = r_0[i][j];
            result_v[0][i][j] = v_0[i][j];
        }
    }

    //float h = 0.016;
    //float L = (4*6*6*6)/rho;
    //int num_timesteps = 1500;
    //int N = 864;
    printf("result_r");
    for (int t = 1; t<1500; t++) {
        verlet_step(result_r[t-1], result_v[t-1], t, M, rho);
    }

    printf("result_r");
    free(r_0);
    free(v_0);
    free(result_r);
    free(result_v);
    return 0;
}
