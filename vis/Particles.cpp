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
    int nx = 5;
    int ny = 5;
    int nz = 5;

    
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

double Particles::smoothing_kernel(double r, double h) {
    if (r <= h and r >= 0) {
        double first = 315/(64*M_PI*pow(h, 9));
        double second = pow(pow(h, 2) - pow(r, 2), 3);
        return first*second;
    } else {
        return 0.0;
    }
    
}

double Particles::magnitude(glm::dvec3 r) {
    return pow((pow(r.x, 2) + pow(r.y, 2) + pow(r.z, 2)), 0.5);
}

glm::dvec3 Particles::spiky_kernel(glm::dvec3 r, double h) {
    //printf("SPIKY R: x: %f, y: %f, z: %f\n", r.x, r.y, r.z);
    if (magnitude(r) <= h and magnitude(r) >= 0.0) {
        double first = 15.0 / (M_PI * pow(h, 6));
        double second = pow((h - magnitude(r)), 2);
        //printf("magnitude = %f\n", (length(r)));
        glm::dvec3 third = r / (magnitude(r));
        glm::dvec3 result = first*second*third;
        //printf("spiky: %f %f %f\n", result.x, result.y, result.z);
        return result;
    }
    return glm::dvec3(0, 0, 0);

}

glm::dvec3 Particles::CiGradient(glm::dvec3 r, double h) {
    return -1.0 / rest_density * spiky_kernel(r, h);
}

void Particles::calculate_lambda_and_density() {
    for(Particle &par : particles) {
        //for all particles, find lambda i
        double density = 0.0;
        // glm::dvec3 sum_gradients = glm::dvec3(0.0, 0.0, 0.0);
        double sum_gradients = 0.0;
        glm::dvec3 self_sum = glm::dvec3(0, 0, 0);
        if (par.neighbors.size() == 0) {
            par.lambda = 0.;
            par.density = 0.;
        } else {
            for (const Particle* neighbor : par.neighbors) { //iterate through neighbors
                //calculate density
                double dist = magnitude(par.new_p - neighbor->new_p);
                if (dist > kernel_size) {
                  continue;  
                } 
                density += smoothing_kernel(dist, kernel_size);
                if(neighbor->new_p.x == par.new_p.x and neighbor->new_p.y == par.new_p.y and neighbor->new_p.z == par.new_p.z) {
                    continue;
                }
                //calculate lambda
                glm::dvec3 r = par.new_p - neighbor->new_p;
                //printf("par = (%f, %f,%f), neighbor = (%f, %f, %f), par-neighbor = (%f, %f, %f) \n", par.new_p.x, par.new_p.y, par.new_p.z, neighbor->new_p.x, neighbor->new_p.y, neighbor->new_p.z, r.x, r.y, r.z);
                //printf("BEFORE: Lambda: %f %f %f\n", lambda.x, lambda.y, lambda.z);
                self_sum += CiGradient(r, kernel_size);
                sum_gradients += pow(magnitude(CiGradient(r, kernel_size)), 2);
            }
            double C = (density / rest_density) - 1.0;
            par.density = density;
            sum_gradients = magnitude(self_sum) + sum_gradients;
            par.lambda = -C/ (sum_gradients+epsilon);
        }
        if (par.lambda > 0) par.lambda = 0.0;
        //printf("lambda: %f\n", par.lambda); 
    } 
}

void Particles::calculate_density() {
    for(Particle &par : particles) {
        //for all particles, find lambda i
        double density = 0.0;        
        for (const Particle* neighbor : par.neighbors) { //iterate through neighbors
            //calculate density
            double dist = magnitude(par.new_p - neighbor->new_p);
            if (dist > kernel_size) {
              continue;  
            } 
            density += smoothing_kernel(dist, kernel_size);
        }
        par.density = density;
        //printf("density: %f\n", density);
        //printf("lambda: %f\n", par.lambda); 
    } 
}

void Particles::calculate_lambda() {
    int counter = 0;
    for(Particle &par : particles) {
        //for all particles, find lambda i
        double sum_gradients = 0.0;
        glm::dvec3 self_sum = glm::dvec3(0, 0, 0);
        for (const Particle* neighbor : par.neighbors) { //iterate through neighbors
            //calculate density
            if(neighbor->new_p.x == par.new_p.x and neighbor->new_p.y == par.new_p.y and neighbor->new_p.z == par.new_p.z) {
                continue;
            }
            //calculate lambda
            glm::dvec3 r = par.new_p - neighbor->new_p;
            self_sum += CiGradient(r, kernel_size);
            sum_gradients += pow(magnitude(CiGradient(r, kernel_size)), 2);
        }
        double C = (par.density / rest_density) - 1.0;
        sum_gradients = magnitude(self_sum) + sum_gradients;
        par.lambda = -C/ (sum_gradients+epsilon);
        if (par.lambda < 0) par.lambda = 0.0;
        if (counter == 1) printf("#%d C is %f\n", counter, C);
        counter += 1;
    }
    
}

double Particles::dot(glm::dvec3 a, glm::dvec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


void Particles::step() {
    for(Particle &par : particles) { 
        if (isnan(par.p.x) or isnan(par.p.y) or isnan(par.p.z)) {
            printf("NANANANANNAN NOOOO NOOO NOOO HELP \n");
        }
        //only z is affected by gravity  
        par.v = par.v + dt * g;
        par.new_p = par.p + dt * par.v;
    }

    hash_grid();
    find_neighbors();

    for(int i = 0; i <= nIters; i++) {
        printf("====iteration %d=====\n", i);

        calculate_density();
        calculate_lambda();
        
        for(Particle &par : particles) {
            //for all particles, find position delta 
            glm::dvec3 pos_delta = glm::dvec3(0.0, 0.0, 0.0);
            for (Particle* neighbor : par.neighbors) {
                if(neighbor->new_p.x == par.new_p.x and neighbor->new_p.y == par.new_p.y and neighbor->new_p.z == par.new_p.z) {
                    continue;
                }
                double smooth = smoothing_kernel(magnitude(par.new_p - neighbor->new_p), kernel_size);
                double smooth_q = smoothing_kernel(0.0, kernel_size);
                double s_corr = -k * pow(smooth/smooth_q, n);
                //printf("s_corr: %f\n", s_corr);
                double lambdas = par.lambda + neighbor->lambda;
                //printf("par: %f, neighbor:%f\n", par.lambda, neighbor->lambda);
                glm::dvec3 spiky = spiky_kernel(par.new_p - neighbor->new_p, kernel_size);
                //printf("spiky = (%f, %f, %f)\n", spiky.x, spiky.y, spiky.z);
                pos_delta += (lambdas) * spiky;
                //printf("pd = (%f, %f, %f)\n", pos_delta.x, pos_delta.y, pos_delta.z);
            }

            pos_delta = 1.0/(rest_density) * pos_delta;
            //printf("pd = (%f, %f, %f)\n", pos_delta.x, pos_delta.y, pos_delta.z);
            //collision handling
            par.new_p = par.new_p + pos_delta;

            collision();
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
}

void Particles::collision() {
    for(Particle &par : particles) { 
        if (par.new_p.x <= -1.) {
            par.new_p.x = -1;
        }

        if (par.new_p.x >= 1.) {
            par.new_p.x = 1;
        } 

        if (par.new_p.y <= -1.) {
            par.new_p.y = -1;

        } 
        if (par.new_p.y >= 1.) {
            par.new_p.y = 1;
        }

        if (par.new_p.z <= -1.) {
            par.new_p.z = -1;

        } 
        if (par.new_p.z >= 1.) {
            par.new_p.z = 1;
        }
    }
}
//lambda
//density
//boxcollision
//update

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
        //printf("neighbors: %lu\n", par.neighbors.size());
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

