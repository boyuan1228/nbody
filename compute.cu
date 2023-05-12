//#include <stdlib.h>
#include <math.h>

#include <stdio.h>
#include <assert.h>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "vector.h"
#include "config.h"

#define KERNEL_SIZE 1024

__global__ void cuda_compute_pairwise(vector3** accels, vector3* hPos, double* mass)
{
    int i = threadIdx.x;
    int j = blockIdx.x;
    int k;
    if (i < NUMENTITIES && j < NUMENTITIES)
    {
        if (i == j)
        {
            FILL_VECTOR(accels[i][j], 0, 0, 0);
        }
        else
        {
            vector3 distance;
            for (k = 0; k < 3; k++) distance[k] = hPos[i][k] - hPos[j][k];
            double magnitude_sq = distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2];
            double magnitude = sqrt(magnitude_sq);
            double accelmag = -1 * GRAV_CONSTANT * mass[j] / magnitude_sq;
            FILL_VECTOR(accels[i][j], accelmag * distance[0] / magnitude, accelmag * distance[1] / magnitude, accelmag * distance[2] / magnitude);
        }
    }
}

__global__ void cuda_compute_sumup(vector3** accels, vector3* hPos, vector3* hVel)
{
    int i = blockIdx.x;
    int j;
    int k;
    if (i < NUMENTITIES)
    {
        vector3 accel_sum = { 0,0,0 };
        for (j = 0; j < NUMENTITIES; j++)
        {
            for (k = 0; k < 3; k++)
                accel_sum[k] += accels[i][j][k];
        }
        //compute the new velocity based on the acceleration and time interval
        //compute the new position based on the velocity and time interval
        for (k = 0; k < 3; k++)
        {
            hVel[i][k] += accel_sum[k] * INTERVAL;
            hPos[i][k] = hVel[i][k] * INTERVAL;
        }
    }
}

void compute()
{
    cudaError_t err;

    int i;
    //vector3* values = (vector3*)malloc(sizeof(vector3) * NUMENTITIES * NUMENTITIES);
    //vector3** accels = (vector3**)malloc(sizeof(vector3*) * NUMENTITIES);
    //make an acceleration matrix which is NUMENTITIES squared in size;
    //vector3* values;
    //cudaMallocManaged(&values, sizeof(vector3) * NUMENTITIES * NUMENTITIES);
    //vector3** accels;
    //cudaMallocManaged(&accels, sizeof(vector3) * NUMENTITIES);

    for (i = 0; i < NUMENTITIES; i++)
        accels[i] = &values[i * NUMENTITIES];

    //first compute the pairwise accelerations.  Effect is on the first argument.
    cuda_compute_pairwise << < KERNEL_SIZE, KERNEL_SIZE >> > (accels, hPos, mass);
    err = cudaGetLastError();
    assert(err == 0);
    cudaDeviceSynchronize();

    //sum up the rows of our matrix to get effect on each entity, then update velocity and position.
    cuda_compute_sumup << < KERNEL_SIZE, 1 >> > (accels, hPos, hVel);
    err = cudaGetLastError();
    assert(err == 0);
    cudaDeviceSynchronize();

    //cudaFree(accels);
    //cudaFree(values);
}

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
//void compute()
//{
//    //make an acceleration matrix which is NUMENTITIES squared in size;
//    int i, j, k;
//    vector3* values = (vector3*)malloc(sizeof(vector3) * NUMENTITIES * NUMENTITIES);
//    vector3** accels = (vector3**)malloc(sizeof(vector3*) * NUMENTITIES);
//    for (i = 0; i < NUMENTITIES; i++)
//        accels[i] = &values[i * NUMENTITIES];
//    //first compute the pairwise accelerations.  Effect is on the first argument.
//    for (i = 0; i < NUMENTITIES; i++)
//    {
//        for (j = 0; j < NUMENTITIES; j++)
//        {
//            if (i == j)
//            {
//                FILL_VECTOR(accels[i][j], 0, 0, 0);
//            }
//            else
//            {
//                vector3 distance;
//                for (k = 0; k < 3; k++) distance[k] = hPos[i][k] - hPos[j][k];
//                double magnitude_sq = distance[0] * distance[0] + distance[1] * distance[1] + distance[2] * distance[2];
//                double magnitude = sqrt(magnitude_sq);
//                double accelmag = -1 * GRAV_CONSTANT * mass[j] / magnitude_sq;
//                FILL_VECTOR(accels[i][j], accelmag * distance[0] / magnitude, accelmag * distance[1] / magnitude, accelmag * distance[2] / magnitude);
//            }
//        }
//    }
//    //sum up the rows of our matrix to get effect on each entity, then update velocity and position.
//    for (i = 0; i < NUMENTITIES; i++)
//    {
//        vector3 accel_sum = { 0,0,0 };
//        for (j = 0; j < NUMENTITIES; j++)
//        {
//            for (k = 0; k < 3; k++)
//                accel_sum[k] += accels[i][j][k];
//        }
//        //compute the new velocity based on the acceleration and time interval
//        //compute the new position based on the velocity and time interval
//        for (k = 0; k < 3; k++)
//        {
//            hVel[i][k] += accel_sum[k] * INTERVAL;
//            hPos[i][k] = hVel[i][k] * INTERVAL;
//        }
//    }
//    free(accels);
//    free(values);
//}
