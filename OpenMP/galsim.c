#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include "../graphics/graphics.h"

const double epsilon_0 = 1e-3; // Softening factor

/*
 * Function to get the current wall time
 * Code snippet from timings.c
 * Course: High Performance Programming, Uppsala University
 * Author: Elias Rudberg <elias.rudberg@it.uu.se>
 *
 */
static double get_wall_time()
{
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

    // Read the particles
    size_t read_count = 0;
    for (int i = 0; i < N; i++)
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
    for (int i = 0; i < N; i++)
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

// Quadtree node structure
typedef struct QuadTreeNode
{
    double mass;                      // Total mass of particles in this node
    double x_cm, y_cm;                // Center of mass coordinates
    double x_min, x_max;              // Bounding box of this node x-axis
    double y_min, y_max;              // Bounding box of this node y-axis
    struct QuadTreeNode *children[4]; // Array of pointers to child nodes
    int is_leaf;                      // Flag to indicate if the node is a leaf
    int particle_index;               // Index of the particle if this node is a leaf
} QuadTreeNode;

// Enum for me better understanding
typedef enum
{
    NW = 0,
    NE = 1,
    SW = 2,
    SE = 3
} Quadrant;

// Function to create a new quadtree node
QuadTreeNode *create_node(double x_min, double x_max, double y_min, double y_max)
{
    QuadTreeNode *node = (QuadTreeNode *)malloc(sizeof(QuadTreeNode));
    node->mass = 0.0;
    node->x_cm = 0.0;
    node->y_cm = 0.0;
    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    node->is_leaf = 1;
    node->particle_index = -1; // No particle initially(empty node)
    for (int i = 0; i < 4; i++)
        node->children[i] = NULL; // Initialize children to NULL
    return node;
}

// Function to insert a particle into the quadtree
void insert_particle(QuadTreeNode *node, double x, double y, double mass, int index)
{
    if (node->is_leaf && node->particle_index == -1)
    {
        node->x_cm = x;
        node->y_cm = y;
        node->mass = mass;
        node->particle_index = index;
        return;
    }

    if (!node->is_leaf)
    {
        double old_mass = node->mass;
        node->mass += mass;
        node->x_cm = (node->x_cm * old_mass + x * mass) / node->mass;
        node->y_cm = (node->y_cm * old_mass + y * mass) / node->mass;

        double mid_x = (node->x_min + node->x_max) / 2;
        double mid_y = (node->y_min + node->y_max) / 2;

        Quadrant quadrant = (x <= mid_x) ? ((y <= mid_y) ? NW : SW) : ((y <= mid_y) ? NE : SE);

        if (node->children[quadrant] == NULL)
        {
            switch (quadrant)
            {
            case NW:
                node->children[NW] = create_node(node->x_min, mid_x, node->y_min, mid_y);
                break;
            case NE:
                node->children[NE] = create_node(mid_x, node->x_max, node->y_min, mid_y);
                break;
            case SW:
                node->children[SW] = create_node(node->x_min, mid_x, mid_y, node->y_max);
                break;
            case SE:
                node->children[SE] = create_node(mid_x, node->x_max, mid_y, node->y_max);
                break;
            }
        }

        insert_particle(node->children[quadrant], x, y, mass, index);
        return;
    }

    if (node->is_leaf)
    {
        double mid_x = (node->x_min + node->x_max) / 2;
        double mid_y = (node->y_min + node->y_max) / 2;

        node->is_leaf = 0;
        int existing_index = node->particle_index;
        double existing_x = node->x_cm;
        double existing_y = node->y_cm;
        double existing_mass = node->mass;

        node->children[NW] = create_node(node->x_min, mid_x, node->y_min, mid_y);
        node->children[NE] = create_node(mid_x, node->x_max, node->y_min, mid_y);
        node->children[SW] = create_node(node->x_min, mid_x, mid_y, node->y_max);
        node->children[SE] = create_node(mid_x, node->x_max, mid_y, node->y_max);

        insert_particle(node, existing_x, existing_y, existing_mass, existing_index);
        insert_particle(node, x, y, mass, index);
    }
}

// Function to compute mass distribution in the quadtree
void compute_mass_distribution(QuadTreeNode *node)
{
    if (node == NULL || node->is_leaf)
        return;

    node->mass = 0.0;
    node->x_cm = 0.0;
    node->y_cm = 0.0;

    for (int i = 0; i < 4; i++)
    {
        if (node->children[i] != NULL)
        {
            compute_mass_distribution(node->children[i]);
            node->mass += node->children[i]->mass;
            node->x_cm += node->children[i]->x_cm * node->children[i]->mass;
            node->y_cm += node->children[i]->y_cm * node->children[i]->mass;
        }
    }

    if (node->mass != 0.0)
    {
        node->x_cm /= node->mass;
        node->y_cm /= node->mass;
    }
}

// Function to calculate forces using the Barnes-Hut algorithm
void calculate_forces_barnes_hut(QuadTreeNode *node, double x, double y, double mass, double theta, double G, double *fx, double *fy, int current_particle_index)
{
    if (node == NULL)
        return;

    double dx = node->x_cm - x;
    double dy = node->y_cm - y;
    double dist = sqrt(dx * dx + dy * dy) + epsilon_0;

    if ((node->x_max - node->x_min) / dist < theta || node->is_leaf)
    {
        if (node->particle_index != -1 && node->particle_index != current_particle_index)
        {
            double r_cube = dist * dist * dist;
            double force = (G * mass * node->mass) / r_cube;
            #pragma omp atomic
            *fx += force * dx;
            #pragma omp atomic
            *fy += force * dy;
        }
    }
    else
    {
        for (int i = 0; i < 4; i++)
        {
            calculate_forces_barnes_hut(node->children[i], x, y, mass, theta, G, fx, fy, current_particle_index);
        }
    }
}

// Function to free the quadtree
void free_quadtree(QuadTreeNode *node)
{
    if (node == NULL)
        return;

    for (int i = 0; i < 4; i++)
    {
        if (node->children[i] != NULL)
        {
            free_quadtree(node->children[i]);
        }
    }
    free(node);
}

// Function to draw the quadtree
void draw_quadtree(QuadTreeNode *node, double W, double H)
{
    if (node == NULL)
        return;

    // Draw the bounding box of the node
    int nodeType = (node->x_max - node->x_min < W || node->y_max - node->y_min < H) ? 0 : 1;
    DrawRectangle(node->x_min, node->y_min, W, H, node->x_max - node->x_min, node->y_max - node->y_min, nodeType);

    // Recursively draw the children nodes
    for (int i = 0; i < 4; i++)
    {
        draw_quadtree(node->children[i], W, H);
    }
}

// Function to run the simulation
void run_simulation(int N, const char *filename, int nsteps, double delta_t, double G, double theta, int graphics, int num_threads)
{
    double start = get_wall_time();

    double *__restrict__ x = (double *)malloc(N * sizeof(double));
    double *__restrict__ y = (double *)malloc(N * sizeof(double));
    double *__restrict__ mass = (double *)malloc(N * sizeof(double));
    double *__restrict__ vx = (double *)malloc(N * sizeof(double));
    double *__restrict__ vy = (double *)malloc(N * sizeof(double));
    double *__restrict__ brightness = (double *)malloc(N * sizeof(double));
    double *__restrict__ fx = (double *)malloc(N * sizeof(double));
    double *__restrict__ fy = (double *)malloc(N * sizeof(double));
    double *__restrict__ ax = (double *)malloc(N * sizeof(double));
    double *__restrict__ ay = (double *)malloc(N * sizeof(double));

    if (x == NULL || y == NULL || mass == NULL || vx == NULL || vy == NULL || brightness == NULL || fx == NULL || fy == NULL || ax == NULL || ay == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    read_particles(filename, x, y, mass, vx, vy, brightness, N);

    if (graphics)
    {
        InitializeGraphics("Simulation", 800, 800);
        SetCAxes(0, 1);
    }

    omp_set_num_threads(num_threads);

    // Compute initial forces and accelerations
    QuadTreeNode *root = create_node(0.0, 1.0, 0.0, 1.0);
    for (int i = 0; i < N; i++)
    {
        insert_particle(root, x[i], y[i], mass[i], i);
    }
    compute_mass_distribution(root);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++)
    {
        double fx_temp = 0.0, fy_temp = 0.0;
        calculate_forces_barnes_hut(root, x[i], y[i], mass[i], theta, G, &fx_temp, &fy_temp, i);
        fx[i] = fx_temp;
        fy[i] = fy_temp;
        ax[i] = fx[i] / mass[i];
        ay[i] = fy[i] / mass[i];
    }
    free_quadtree(root);

    for (int step = 0; step < nsteps; step++)
    {
        // Update positions
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; i++)
        {
            x[i] += vx[i] * delta_t + 0.5 * ax[i] * delta_t * delta_t;
            y[i] += vy[i] * delta_t + 0.5 * ay[i] * delta_t * delta_t;
        }

        // Recompute forces and accelerations
        root = create_node(0.0, 1.0, 0.0, 1.0);
        for (int i = 0; i < N; i++)
        {
            insert_particle(root, x[i], y[i], mass[i], i);
        }
        compute_mass_distribution(root);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; i++)
        {
            double fx_temp = 0.0, fy_temp = 0.0;
            calculate_forces_barnes_hut(root, x[i], y[i], mass[i], theta, G, &fx_temp, &fy_temp, i);
            fx[i] = fx_temp;
            fy[i] = fy_temp;
        }

        // Update velocities
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; i++)
        {
            double ax_new = fx[i] / mass[i];
            double ay_new = fy[i] / mass[i];
            vx[i] += 0.5 * (ax[i] + ax_new) * delta_t;
            vy[i] += 0.5 * (ay[i] + ay_new) * delta_t;
            ax[i] = ax_new;
            ay[i] = ay_new;
        }

        if (graphics)
        {
            ClearScreen();
            for (int i = 0; i < N; i++)
            {
                DrawCircle(x[i], y[i], 1.0, 1.0, 0.002, 0);
            }
            draw_quadtree(root, 1.0, 1.0);
            Refresh();
            usleep(3000);
        }

        free_quadtree(root);
    }

    write_particles("result.gal", x, y, mass, vx, vy, brightness, N);

    double end = get_wall_time();
    printf("Time to run simulation: %f\n", end - start);

    if (graphics)
    {
        while (!CheckForQuit())
        {
            usleep(1000);
        }
        CloseDisplay();
    }

    free(x);
    free(y);
    free(mass);
    free(vx);
    free(vy);
    free(brightness);
    free(fx);
    free(fy);
    free(ax);
    free(ay);
}

// Main function
int main(int argc, char *argv[])
{
    double start = get_wall_time();

    // Check for correct number of arguments
    if (argc != 8)
    {
        fprintf(stderr, "Expected 7 arguments: N filename nsteps delta_t theta graphics num_threads\n");
        return 1;
    }

    // Parse command-line arguments
    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    double theta = atof(argv[5]);
    int graphics = atoi(argv[6]);
    int num_threads = atoi(argv[7]);

    // Check if the command-line arguments are valid
    if (N <= 0 || nsteps < 0 || delta_t <= 0 || theta < 0 || theta > 1 || (graphics != 0 && graphics != 1) || num_threads <= 0)
    {
        fprintf(stderr, "Invalid arguments: N, nsteps, delta_t must be positive, theta must be between 0 and 1, graphics must be 0 or 1, num_threads must be positive\n");
        return 1;
    }

    // Calculate the gravitational constant
    const double G = 100.0 / N;

    // Run the simulation
    run_simulation(N, filename, nsteps, delta_t, G, theta, graphics, num_threads);

    double end = get_wall_time();
    printf("Total time: %f\n", end - start);

    return 0;
}