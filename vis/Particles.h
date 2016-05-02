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

    // void wireBox(int myVar);
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

    // double rest_density;
    // double kernel_size;
    // double radius;
    // double k;
    // double n;
    // double q;
    // double epsilon;
    // double nIters;
    // double dt;
    // double g;

    double d = 0.1; //resting density
    double g = -9.8; //gravisty
    double kernel_size = d*1.4;
    double radius = d*0.45;
    double k = 0.001; //artificial pressure
    double n = 4; //artificial pressure
    double q = 0; //artificial pressure
    double epsilon = 1e3;
    double nIters = 1;
    double rest_density = 1/(d*d*d);
    double dt = 0.001;

    struct Grid
	{
	  int x;
	  int y;
	  int z;

	  bool operator==(const Grid &other) const
	  { return (x == other.x
	            && y == other.y
	            && z == other.z);
	  }
	};

	struct GridHasher
	{
	  std::size_t operator()(const Grid& g) const
	  {
	    return floor(g.x)+floor(g.y)*1300583+floor(g.z)*105607;;
	  }
	};

	std::vector<Particle> particles;
    std::unordered_map<Grid, std::vector<Particle>, GridHasher> hashGrid;

};

#endif /* PARTICLES_H */

