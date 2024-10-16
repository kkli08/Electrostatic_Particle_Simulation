//
// Created by Damian Li on 2024-10-16.
//

#include <iostream>
#include "Particle.h"

using namespace std;

int main() {
  cout<<"Hello World!"<<endl;
  // Create particles
  Particle p1(94393, -27903, 'p');
  Particle p2(10000, 20000, 'e');

  // Create Coulomb force calculator
  CoulombForceCalculator calculator;

  // Compute and display the force between particles
  double force = calculator.computeCoulombForce(p1, p2);
  std::cout << "The force between the particles is: " << force << " N" << std::endl;

  return 0;
}