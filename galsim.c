#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "graphics/graphics.h"
#include <sys/time.h>

const double epsilon_0 = 1e-3; // Softening factor

/*
 * Function to get the current wall time
 * Code snippet from timings.c
 * Course: High Performance Programming, Uppsala University
 * Author: Elias Rudberg <elias.rudberg@it.uu.se>
 * 
*/
static double get_wall_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

// Function to read particles from a file
void read_particles(const char *filename, double *x, double *y, double *mass, double *vx, double *vy, double *brightness, int N)
{
    double start = get_wall_time();
    FILE *file = fopen(filename, "rb");
    if (!file)
    {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Get the size of the file
    fseek(file, 0, SEEK_END);
    long filesize = ftell(file);
    rewind(file);

    // Calculate the number of particles in the file
    int num_particles = filesize / (6 * sizeof(double));

    // Check that the number of particles in the file matches the input
    if (N != num_particles)
    {
        fprintf(stderr, "Warning: N (%d) is not equal to the number of particles in the file (%d)\n", N, num_particles);
        exit(EXIT_FAILURE);
    }

    // Read the particles
    size_t read_count = 0; 
    for(int i = 0; i < N; i++)
    {
        read_count += fread(&x[i], sizeof(double), 1, file);
        read_count += fread(&y[i], sizeof(double), 1, file);
        read_count += fread(&mass[i], sizeof(double), 1, file);
        read_count += fread(&vx[i], sizeof(double), 1, file);
        read_count += fread(&vy[i], sizeof(double), 1, file);
        read_count += fread(&brightness[i], sizeof(double), 1, file);
    }


    if (read_count != 6 * N)
    {
        fprintf(stderr, "Error reading particles from file\n");
        exit(EXIT_FAILURE);
    }

    fclose(file);
    double end = get_wall_time();
    printf("Time to read particles: %f\n", end - start);
}

// Function to write particles to a file
void write_particles(const char *filename, double *x, double *y, double *mass, double *vx, double *vy, double *brightness, int N)
{   
    double start = get_wall_time();
    
    FILE *file = fopen(filename, "wb");
    if (!file)
    {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Write the particles
    size_t write_count = 0;
    for(int i = 0; i < N; i++)
    {
        write_count += fwrite(&x[i], sizeof(double), 1, file);
        write_count += fwrite(&y[i], sizeof(double), 1, file);
        write_count += fwrite(&mass[i], sizeof(double), 1, file);
        write_count += fwrite(&vx[i], sizeof(double), 1, file);
        write_count += fwrite(&vy[i], sizeof(double), 1, file);
        write_count += fwrite(&brightness[i], sizeof(double), 1, file);
    }

    if (write_count != 6 * N)
    {
        fprintf(stderr, "Error writing particles to file\n");
        exit(EXIT_FAILURE);
    }

    fclose(file);

    double end = get_wall_time();
    printf("Time to write particles: %f\n", end - start);
}

// Function to calculate forces and update particles
void calculate_forces_and_update(double *__restrict__ x, double *__restrict__ y, double *__restrict__ mass, double *__restrict__ vx, double *__restrict__ vy, double *__restrict__ fx, double *__restrict__ fy, int N, double delta_t, double G)
{
    // Zeroing out forces at the start of each function call
    memset(fx, 0, N * sizeof(double));
    memset(fy, 0, N * sizeof(double));

    // Calculate forces
    for (int i = 0; i < N; i++)
    {
        double fx_temp = 0.0;
        double fy_temp = 0.0;
        for (int j = i + 1; j < N; j++)
        {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double r_sqrd = dx * dx + dy * dy;
            double dist = sqrt(r_sqrd) + epsilon_0;
            double r_cube = dist * dist * dist;
            double force_over_dist = (G * mass[i] * mass[j]) / r_cube;

            double fx_force = force_over_dist * dx;
            double fy_force = force_over_dist * dy;

            fx_temp += fx_force;
            fy_temp += fy_force;
            fx[j] -= fx_force;
            fy[j] -= fy_force;
        }

        fx[i] += fx_temp;
        fy[i] += fy_temp;

        // Update velocity
        vx[i] += delta_t * fx[i] / mass[i];
        vy[i] += delta_t * fy[i] / mass[i];

        // Update position
        x[i] += vx[i] * delta_t;
        y[i] += vy[i] * delta_t;
    }
}

// Function to run the simulation
void run_simulation(int N, const char *filename, int nsteps, double delta_t, double G, int graphics)
{
    double start = get_wall_time();

    // Allocate memory for particles and forces
    double *x = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *mass = (double *)malloc(N * sizeof(double));
    double *vx = (double *)malloc(N * sizeof(double));
    double *vy = (double *)malloc(N * sizeof(double));
    double *brightness = (double *)malloc(N * sizeof(double));
    double *fx = (double *)malloc(N * sizeof(double));
    double *fy = (double *)malloc(N * sizeof(double));

    if (x == NULL || y == NULL || mass == NULL || vx == NULL || vy == NULL || brightness == NULL || fx == NULL || fy == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    read_particles(filename, x, y, mass, vx, vy, brightness, N);

    // Initialize graphics if requested
    if (graphics)
    {
        InitializeGraphics("Simulation", 800, 800);
        SetCAxes(0, 1);
    }

    // Main simulation loop
    for (int step = 0; step < nsteps; step++)
    {
        // Calculate forces and update particles
        calculate_forces_and_update(x, y, mass, vx, vy, fx, fy, N, delta_t, G);

        // Draw particles if graphics is enabled
        if (graphics)
        {
            ClearScreen();
            for (int i = 0; i < N; i++)
            {
                // Draw each particle as a circle
                DrawCircle(x[i], y[i], 1.0, 1.0, 0.003, 0.1);
            }
            Refresh();

            // Delay to control animation speed
            usleep(3000);
        }
    }

    double end = get_wall_time();
    printf("Time to run simulation: %f\n", end - start);
    

    // Close graphics if q is pressed
    if (graphics)
    {
        while (!CheckForQuit())
        {
            usleep(1000);
        }
        CloseDisplay();
    }


    // Write particle data to a file
    write_particles("result.gal", x, y, mass, vx, vy, brightness, N);
    free(x);
    free(y);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);
    free(fx);
    free(fy);
}

// Main function
int main(int argc, char *argv[])
{
    double start = get_wall_time();

    // Check for correct number of arguments
    if (argc != 6)
    {
        fprintf(stderr, "Expected 5 arguments: N filename nsteps delta_t graphics\n");
        return 1;
    }

    // Parse command-line arguments
    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);

    // Check if the command-line arguments are valid
    if (N <= 0 || nsteps < 0 || delta_t <= 0 || (graphics != 0 && graphics != 1))
    {
        fprintf(stderr, "Invalid arguments: N, nsteps and, delta_t must be positive, graphics must be 0 or 1\n");
        return 1;
    }

    // Calculate the gravitational constant
    const double G = 100.0 / N;

    // Run the simulation
    run_simulation(N, filename, nsteps, delta_t, G, graphics);

    double end = get_wall_time();
    printf("Total time: %f\n", end - start);

    return 0;
}
