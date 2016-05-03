/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.cpp
 * Author: swl
 * 
 * Created on April 15, 2016, 12:16 PM
 */

#include "Particles.h"

Particles::Particles() 
{
    int nx = 3;
    int ny = 3;
    int nz = 3;

    for(int x=0; x<nx; x++)
    {
        for(int y=0; y<ny; y++)
        {
            for(int z=0; z<nz; z++)
            {
                Particle par;
                par.p = glm::dvec3((x+0.5-nx*0.5)*d, (y+0.5+ny*0.5)*d-1.0, (z+0.5-nz*0.5)*d);
                par.new_p = glm::dvec3((x+0.5-nx*0.5)*d, (y+0.5+ny*0.5)*d-1.0, (z+0.5-nz*0.5)*d);
                particles.push_back(par);
            }
        }
    }
}

double Particles::smoothing_kernel(glm::dvec3 r, double h) {
    double first = 315/(64*M_PI*pow(h, 9));
    double second = pow(pow(h, 2) - pow(magnitude(r), 2), 3);
    return first*second;
}

double Particles::magnitude(glm::dvec3 r) {
    return pow((pow(r.x, 2) + pow(r.y, 2) + pow(r.z, 2)), 0.5);
}

glm::dvec3 Particles::spiky_kernel(glm::dvec3 r, double h) {
    double first = 45.0 / (M_PI * pow(h, 6));
    double second = pow((h - magnitude(r)), 2);
    //printf("magnitude = %f\n", (length(r)));
    glm::dvec3 third = r / (magnitude(r));
    glm::dvec3 result = first*second*third;
    // printf("first: %f, second: %f\n", first, second);
    // printf("third: %f, %f, %f\n", third.x, third.y, third.z);
    //printf("result: %f %f %f\n", result.x, result.y, result.z);
    return result;
}

void Particles::step() {
    for(Particle &par : particles) { 
        if (isnan(par.p.x) or isnan(par.p.y) or isnan(par.p.z)) {
            printf("NANANANANNAN NOOOO NOOO NOOO HELP \n");
        }
        //printf("x: %f, y: %f, z: %f\n", par.p.x, par.p.y, par.p.z); 
        //only z is affected by gravity  
        par.v.y = par.v.y + dt * g;
        par.new_p.y = par.p.y + dt * par.v.y;
    }

    printf("checkpoint 1\n");
    hash_grid();
    find_neighbors();

    for(int i = 0; i <= nIters; i++) {
        printf("====iteration %d=====\n", i);
        for(Particle &par : particles) {
            //for all particles, find lambda i
            double density = 0;
            glm::dvec3 lambda = glm::dvec3(0.0, 0.0, 0.0);
            if (par.neighbors.size() == 0) {
                par.lambda = 0.;
                par.density = 0.;
            } else {
                for (const Particle* neighbor : par.neighbors) { //iterate through neighbors
                    //calculate density
                    if(neighbor->new_p.x == par.new_p.x and neighbor->new_p.y == par.new_p.y and neighbor->new_p.z == par.new_p.z) {
                        continue;
                    }
                    density += smoothing_kernel(par.new_p - neighbor->new_p, kernel_size);
                    //calculate lambda
                    glm::dvec3 r = par.new_p - neighbor->new_p;
                    printf("par = (%f, %f,%f), neighbor = (%f, %f, %f), par-neighbor = (%f, %f, %f) \n", par.new_p.x, par.new_p.y, par.new_p.z, neighbor->new_p.x, neighbor->new_p.y, neighbor->new_p.z, r.x, r.y, r.z);
                    printf("BEFORE: Lambda: %f %f %f\n", lambda.x, lambda.y, lambda.z);
                    lambda += spiky_kernel(par.new_p - neighbor->new_p, kernel_size);
                }
                double C = (density / rest_density) - 1;
                par.density = density;
                double iterim = magnitude((1/rest_density) * lambda);
                par.lambda = iterim ? -C/ pow(iterim, 2) : 0;
            }
            
        }
        for(Particle &par : particles) {
            //for all particles, find position delta 
            glm::dvec3 pd = glm::dvec3(0.0, 0.0, 0.0);
            for (const Particle* neighbor : par.neighbors) {
                if(neighbor->new_p.x == par.new_p.x and neighbor->new_p.y == par.new_p.y and neighbor->new_p.z == par.new_p.z) {
                    continue;
                }
                double lambdas = par.lambda + neighbor->lambda;
                //printf("par: %f, neighbor:%f\n", par.lambda, neighbor->lambda);
                glm::dvec3 spiky = spiky_kernel(par.new_p - neighbor->new_p, kernel_size);
                //printf("spiky = (%f, %f, %f)\n", spiky.x, spiky.y, spiky.z);
                pd += (lambdas) * spiky;
                // printf("pd = (%f, %f, %f)\n", pd.x, pd.y, pd.z);
            }
            pd = 1/(rest_density) * pd;
            //printf("pd = (%f, %f, %f)\n", pd.x, pd.y, pd.z);
            //collision handling
            par.new_p = par.new_p + pd;
            if (par.new_p.x <= -1.) par.new_p.x = -1.;
            if (par.new_p.x >= 1.) par.new_p.x = 1.;

            if (par.new_p.y <= -1.) par.new_p.y = -1.;
            if (par.new_p.y >= 1.) par.new_p.y = 1.;

            if (par.new_p.z <= -1.) par.new_p.z = -1;
            if (par.new_p.z >= 1.) par.new_p.z = 1;
            //printf("x: %f, y: %f, z: %f\n", par.new_p.x, par.new_p.y, par.new_p.z);
        }
    }

    for (Particle &par : particles) {
        //update velocity t+1
        par.v = (par.new_p - par.p) / dt;
        //apply viscosity & vorticity
        
        //update position t+1
        par.p = par.new_p;   
        
    }

    // printf("checkpoint 4\n");
}

// returns hash value of a particle
// int Particles::hash(double x, double y, double z) {
//     return floor(x/hgSize)+floor(y/hgSize)*1300583+floor(z/hgSize)*105607;
// }

// updates hash_grid()
void Particles::hash_grid() {
    //build the hash grid
    std::unordered_map<Grid, std::vector<Particle*>, GridHasher> newHashGrid;
    for(Particle &par : particles) {  
        Grid newGrid;
        newGrid.x = floor(par.new_p.x/kernel_size); 
        newGrid.y = floor(par.new_p.y/kernel_size);
        newGrid.z = floor(par.new_p.z/kernel_size);

        auto search = newHashGrid.find(newGrid);
        if (search == newHashGrid.end()) {
            std::vector<Particle*> cell;
            cell.push_back(&par);
            newHashGrid.insert({newGrid, cell});
        } else {
            search->second.push_back(&par);
        }
    }
    hashGrid = newHashGrid;

    // for ( auto it = hashGrid.begin(); it != hashGrid.end(); ++it )
    //     std::cout << " " << it->first << ":" << it->second.size();
}

// finds the neighbors of all particles
void Particles::find_neighbors() 
{
    for(Particle &par : particles) {
        par.neighbors.clear();
        Grid newGrid;
        double x = par.new_p.x, y = par.new_p.y, z = par.new_p.z;  
        double fx = floor(x/kernel_size), fy = floor(y/kernel_size), fz = floor(z/kernel_size);
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {
                    // if(i == 0 && j == 0 && k == 0) 
                    //     continue;

                    newGrid.x = fx + i; 
                    newGrid.y = fy + j;
                    newGrid.z = fz + k;


                    // printf("Grid - fx: %d, fy: %d, fz: %d", newGrid.x, newGrid.y, newGrid.z);

                    auto search = hashGrid.find(newGrid);
                    if (hashGrid.find(newGrid) == hashGrid.end()) {
                        // dont do anything
                    } else {
                        // vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
                        // printf("search->second: %lu\n", search->second.size());

                        par.neighbors.insert( par.neighbors.end(), search->second.begin(), search->second.end() );
                    }
                }
            }
        }
        printf("neighbors: %lu\n", par.neighbors.size());
    }
}

void Particles::render() const
{
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 10.0, 10.0, 10.0, 0.0 };
    glShadeModel (GL_SMOOTH);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glColor3f(0.2, 0.5, 0.8);
    glColorMaterial(GL_FRONT, GL_SPECULAR);
    glColor3f(0.9, 0.9, 0.9);
    glColorMaterial(GL_FRONT, GL_AMBIENT);
    glColor3f(0.2, 0.5, 0.8);
    
    for(const Particle &par : particles)
    {    
        
        glPushMatrix();
        glTranslatef(par.p.x, par.p.y, par.p.z);
        glutSolidSphere(0.05, 10, 10);
        glPopMatrix();
    }
    
    glPopAttrib();
    glutWireCube(2.0);
}

