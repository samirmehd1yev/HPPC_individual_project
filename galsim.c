#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "graphics/graphics.h"

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

// QuadTree node structure
typedef struct QuadTreeNode
{
    double mass;                      // Total mass of particles in this node
    double x_cm, y_cm;                // Center of mass coordinates
    double x_min, x_max;              // Bounds of box in this node (x-axis)
    double y_min, y_max;              // Bounds of box in this node (y-axis)
    struct QuadTreeNode *children[4]; // Pointers to child nodes
    int is_leaf;                      // Flag to check if the node is leaf or not
    int particle_index;               // Index of the particle if this node is a leaf, if empty node -1
} QuadTreeNode;

// Enum for me better understanding
typedef enum{NW = 0,NE = 1,SW = 2,SE = 3} Quadrant;


// Function to create a new quadtree node example: root node spans whole space which is (0,1)X(0,1) (Note: our particles from input_date/ located in this span)
QuadTreeNode *create_node(double x_min, double x_max, double y_min, double y_max)
{
    QuadTreeNode *node = (QuadTreeNode *)malloc(sizeof(QuadTreeNode));
    node->mass = 0.0;          // Initialize the total mass to 0
    node->x_cm = 0.0;          // Initialize the x coordinate of the center of mass to 0
    node->y_cm = 0.0;          // Initialize the y coordinate of the center of mass to 0
    node->x_min = x_min;       // Set the minimum x boundary of the node
    node->x_max = x_max;       // Set the maximum x boundary of the node
    node->y_min = y_min;       // Set the minimum y boundary of the node
    node->y_max = y_max;       // Set the maximum y boundary of the node
    node->is_leaf = 1;         // Initially, the node is a leaf
    node->particle_index = -1; // Identifier for particle. No particle in the node initially, its checked by -1
    for (int i = 0; i < 4; i++)
        node->children[i] = NULL; // Initialize children to NULL
    return node;
}


// Function to insert a particle into the quadtree
void insert_particle(QuadTreeNode *node, double x, double y, double mass, int index)
{
    // Check if the node is a leaf and empty
    if (node->is_leaf && node->particle_index == -1)
    {
        // Assign particle properties to the node
        node->x_cm = x;
        node->y_cm = y;
        node->mass = mass;
        node->particle_index = index;
        return;
    }

    // If the node is not a leaf, update mass and center of mass, then insert into the appropriate child node
    if (!node->is_leaf)
    {
        // Determine which quadrant the particle belongs to
        double mid_x = (node->x_min + node->x_max) / 2;
        double mid_y = (node->y_min + node->y_max) / 2;


        /*Determine which quadrant the particle belongs to
            If (x <= mid_x) is true:
                If (y <= mid_y) is false, the particle is in NW quadrant.
                If (y <= mid_y) is true, the particle is in SW quadrant.
            If (x <= mid_x) is false:
                If (y <= mid_y) is false, the particle is in NE quadrant.
                If (y <= mid_y) is true, the particle is in SE quadrant.
            example: (0.2,0.6) is in NW if inserted to root node (0,1)X(0,1)*/
        Quadrant quadrant = (x <= mid_x) ? ((y > mid_y) ? NW : SW) : ((y > mid_y) ? NE : SE);

        // Create the child node if it doesn't exist
        if (node->children[quadrant] == NULL)
        {
            switch (quadrant)
            {
            case NW:
                node->children[NW] = create_node(node->x_min, mid_x, mid_y, node->y_max);
                break;
            case NE:
                node->children[NE] = create_node(mid_x, node->x_max, mid_y, node->y_max);
                break;
            case SW:
                node->children[SW] = create_node(node->x_min, mid_x, node->y_min, mid_y);
                break;
            case SE:
                node->children[SE] = create_node(mid_x, node->x_max, node->y_min, mid_y);
                break;
            }
        }

        // Recursively insert the particle into the appropriate child node
        insert_particle(node->children[quadrant], x, y, mass, index);
        return;
    }

    // If the node is a leaf but already contains a particle, subdivide and reinsert both particles
    if (node->is_leaf)
    {
        // Ensure the particles are not at the same position
        if (node->x_cm == x && node->y_cm == y)
        {
            fprintf(stderr, "Error: Two particles are at the same position (%f, %f)\n", x, y);
            exit(EXIT_FAILURE);
        }

        // Calculate the midpoint
        double mid_x = (node->x_min + node->x_max) / 2;
        double mid_y = (node->y_min + node->y_max) / 2;

        // Change the node to a non-leaf node
        node->is_leaf = 0;

        // Store existing particle's properties
        int existing_index = node->particle_index;
        double existing_x = node->x_cm;
        double existing_y = node->y_cm;
        double existing_mass = node->mass;

        // Create child nodes
        node->children[NW] = create_node(node->x_min, mid_x, mid_y, node->y_max);
        node->children[NE] = create_node(mid_x, node->x_max, mid_y, node->y_max);
        node->children[SW] = create_node(node->x_min, mid_x, node->y_min, mid_y);
        node->children[SE] = create_node(mid_x, node->x_max, node->y_min, mid_y);

        // Reinsert the existing particle
        insert_particle(node, existing_x, existing_y, existing_mass, existing_index);

        // Insert the new particle
        insert_particle(node, x, y, mass, index);
    }
}


// Function to compute mass distribution in the quadtree
void compute_mass_distribution(QuadTreeNode *node)
{
    // If the node is NULL or a leaf node, no computation is needed, so return
    if (node == NULL || node->is_leaf)
        return;

    // Reset mass and center of mass properties for the current node, bcs we will aggregrate according to its children
    node->mass = 0.0;
    node->x_cm = 0.0;
    node->y_cm = 0.0;

    // Iterate over the child nodes
    for (int i = 0; i < 4; i++)
    {
        // If the current child node exists
        if (node->children[i] != NULL)
        {
            // Recursively compute mass distribution for the child node
            compute_mass_distribution(node->children[i]);
            // Aggregate mass and weighted center of mass properties from child nodes
            node->mass += node->children[i]->mass;
            node->x_cm += node->children[i]->x_cm * node->children[i]->mass;
            node->y_cm += node->children[i]->y_cm * node->children[i]->mass;
        }
    }

    // If the total mass of the current node is nonzero
    if (node->mass != 0.0)
    {
        // Compute the center of mass for the current node
        node->x_cm /= node->mass;
        node->y_cm /= node->mass;
    }
}

// Function to calculate forces using the Barnes-Hut algorithm
void calculate_forces_barnes_hut(QuadTreeNode *node, double x, double y, double mass, double theta_threshold, double G, double *fx, double *fy, int current_particle_index)
{
    if (node == NULL || (node->particle_index == current_particle_index && node->is_leaf))
        return;

    double dx = node->x_cm - x;                          // Difference between x-axis
    double dy = node->y_cm - y;                          // Difference between y-axis
    double dist = sqrt(dx * dx + dy * dy) + epsilon_0;   // Distance with softening factor
    double theta = (node->x_max - node->x_min) / dist;   // Calculating theta value

    // If threshold meets threshold condition or its leaf node
    if (theta <= theta_threshold || node->is_leaf)
    {
        if (node->mass > 0 && node->particle_index != current_particle_index)
        {
            double r_cube = dist * dist * dist;
            double force = (G * mass * node->mass) / r_cube;
            *fx += force * dx;
            *fy += force * dy;
        }
    }
    else  // Otherwise recursively check the children nodes
    {
        for (int i = 0; i < 4; i++)
        {
            calculate_forces_barnes_hut(node->children[i], x, y, mass, theta_threshold, G, fx, fy, current_particle_index);
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
void run_simulation(int N, const char *filename, int nsteps, double delta_t, double G, double theta_threshold, int graphics)
{
    double start = get_wall_time();

    double *x = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *mass = (double *)malloc(N * sizeof(double));
    double *vx = (double *)malloc(N * sizeof(double));
    double *vy = (double *)malloc(N * sizeof(double));
    double *brightness = (double *)malloc(N * sizeof(double));
    double *fx = (double *)malloc(N * sizeof(double));
    double *fy = (double *)malloc(N * sizeof(double));
    double *ax = (double *)malloc(N * sizeof(double));
    double *ay = (double *)malloc(N * sizeof(double));

    if (!x || !y || !mass || !vx || !vy || !brightness || !fx || !fy || !ax || !ay)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    read_particles(filename, x, y, mass, vx, vy, brightness, N);

    if (graphics)
    {
        InitializeGraphics("Simulation", 800, 800);
        SetCAxes(0, 1);
    }

    // Initial force and acceleration calculation
    {
        QuadTreeNode *root = create_node(0.0, 1.0, 0.0, 1.0);
        for (int i = 0; i < N; i++)
            insert_particle(root, x[i], y[i], mass[i], i);
        compute_mass_distribution(root);

        for (int i = 0; i < N; i++) {
            double fx_temp = 0.0, fy_temp = 0.0;
            calculate_forces_barnes_hut(root, x[i], y[i], mass[i], theta_threshold, G, &fx_temp, &fy_temp, i);
            fx[i] = fx_temp;
            fy[i] = fy_temp;
            ax[i] = fx[i] / mass[i];
            ay[i] = fy[i] / mass[i];
        }
        free_quadtree(root);
    }

    for (int step = 0; step < nsteps; step++)
    {
        // Update positions
        for (int i = 0; i < N; i++) {
            x[i] += vx[i] * delta_t + 0.5 * ax[i] * delta_t * delta_t;
            y[i] += vy[i] * delta_t + 0.5 * ay[i] * delta_t * delta_t;
            if (x[i] < 0.0 || x[i] > 1.0 || y[i] < 0.0 || y[i] > 1.0) {
                fprintf(stderr, "Error: Particle moved out of bounds at step %d: (%f, %f)\n", step, x[i], y[i]);
                exit(EXIT_FAILURE);
            }
        }

        // Recalculate forces and accelerations
        QuadTreeNode *root = create_node(0.0, 1.0, 0.0, 1.0);
        for (int i = 0; i < N; i++)
            insert_particle(root, x[i], y[i], mass[i], i);
        compute_mass_distribution(root);

        for (int i = 0; i < N; i++) {
            double fx_temp = 0.0, fy_temp = 0.0;
            calculate_forces_barnes_hut(root, x[i], y[i], mass[i], theta_threshold, G, &fx_temp, &fy_temp, i);
            fx[i] = fx_temp;
            fy[i] = fy_temp;
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
                DrawCircle(x[i], y[i], 1.0, 1.0, 0.002, 0);
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
            usleep(1000);
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
    if (argc != 7)
    {
        fprintf(stderr, "Expected 6 arguments: N filename nsteps delta_t theta_threshold graphics\n");
        return 1;
    }

    // Parse command-line arguments
    int N = atoi(argv[1]);
    const char *filename = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    double theta_threshold = atof(argv[5]);
    int graphics = atoi(argv[6]);

    // Check if the command-line arguments are valid
    if (N <= 0 || nsteps < 0 || delta_t <= 0 || theta_threshold < 0 || theta_threshold > 1 || (graphics != 0 && graphics != 1))
    {
        fprintf(stderr, "Invalid arguments: N, nsteps, delta_t must be positive, theta_threshold must be between 0 and 1, graphics must be 0 or 1\n");
        return 1;
    }

    // Calculate the gravitational constant
    const double G = 100.0 / N;

    // Run the simulation
    run_simulation(N, filename, nsteps, delta_t, G, theta_threshold, graphics);

    double end = get_wall_time();
    printf("Total time: %f\n", end - start);

    return 0;
}
