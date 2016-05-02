/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.h
 * Author: swl
 *
 * Created on April 15, 2016, 12:16 PM
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <glm/glm.hpp>
#include <vector>
#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

#include <unordered_map>

class Particles {
public:
    Particles();
    void render() const;
    void step(); // simulate one frame
    double smoothing_kernel(glm::dvec3 r, double h);
    glm::dvec3 spiky_kernel(glm::dvec3 r, double h);
    void hash_grid();
    void find_neighbors();
private:
    struct Particle
    {
        glm::dvec3 p;
        glm::dvec3 new_p;
        glm::dvec3 v;
        double density;
    	double lambda;
    	std::vector<Particle> neighbors;
    };
    double rest_density;
    double kernel_size;
    double radius;
    double k;
    double n;
    double q;
    double epsilon;
    double nIters;
    double dt;
    double g;
    std::vector<Particle> particles;
    std::unordered_map<int, std::vector<Particle>> hashGrid;

};

#endif /* PARTICLES_H */

