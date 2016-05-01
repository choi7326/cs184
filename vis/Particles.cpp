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
    int nx = 10;
    int ny = 10;
    int nz = 10;
    float d = 0.1; //resting density
    float g = 9.8; //gravisty
    kernel_size = d*1.4;
    radius = d*0.45;
    k = 0.001; //artificial pressure
    n = 4; //artificial pressure
    q = 0; //artificial pressure
    epsilon = 1e3;
    nIters = 10;
    rest_density = 1/(d*d*d);
    dt = 0.001;

    for(int x=0; x<nx; x++)
    {
        for(int y=0; y<ny; y++)
        {
            for(int z=0; z<nz; z++)
            {
                Particle par;
                par.p = glm::dvec3((x+0.5-nx*0.5)*d, (y+0.5)*d-1.0, (z+0.5-nz*0.5)*d);
                particles.push_back(par);
            }
        }
    }
}

double Particles::smoothing_kernel(r, h) {
    double first = 315/(64*PI*pow(h, 9));
    double second = pow(pow(h, 2) - abs(pow(r, 2)), 3);
    return first*second;
}

double Particles::spiky_kernel(r, h) {
    double first = 45 / (PI * pow(h, 6));
    double second = pow((h - abs(r)), 2);
    double third = r / (abs(r));
    return first*second*third;
}

void Particles::step() {   
    for(const Particle &par : particles) {  
        std::vector<Particle> neighbors;  
        par.v.z = par.v.z + dt * g;
        double pt = par.p.z;
        par.p.z = par.p.z + dt * par.v.z;
    }
    for(const Particle &par : particles) {
        //find neighbors & store
        //roughly 9 cells.
    }    
    for(int i = 0; i <= nIters; i++) {
        for(const Particle &par : particles) {
            //for all particles, find lambda i
            double density = 0;
            double lambda = 0;
            for (const Particle &neighbor : neighbors) { //iterate through neighbors
                //calculate density
                density += smoothing_kernel(par.p.z - neighbor.p.z, radius);
                //calculate lambda
                lambda += spiky_kernel(par.p.z - neighbor.p.z, radius);
            }
            double C = (density / rest_density) - 1;
            par.density = density;
            par.lambda = -C/ pow((1/rest_density) * lambda, 2);
        }
        for(const Particle &par : particles) {
           //for all particles, find position delta 
        }
        //update velocity t+1
        par.v.z = (par.p.z - pt) / dt;
        //update position t+1
        par.p.z = par.p.z;
        //apply viscosity & vorticity
    }
    
}


void Particles::hash_grid() {
    //build the hash grid
}

void Particles::find_neighbors() 
{
    //find neighbors
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
}

