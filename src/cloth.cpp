#include "glCanvas.h"

#include <fstream>
#include <iostream>
#include "glm/ext.hpp"
#include "cloth.h"
#include "argparser.h"
#include "utils.h"
#include <math.h>
#include "glm/gtx/string_cast.hpp"
#include "fluid.h"

// ================================================================================
// ================================================================================

float ab_f = 0.00009; //absorption force factor
int count_ = 0; //random number 
double error_ = 0.00001; //allowed eeror for collision detection 
Cloth::Cloth(ArgParser *_args) {
  args =_args;
  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->cloth_file).c_str());
  assert (istr.good());
  std::string token;
  // read in the simulation parameters
  istr >> token >> k_structural; assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear; assert (token == "k_shear");
  istr >> token >> k_bend; assert (token == "k_bend");
  istr >> token >> damping; assert (token == "damping");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction; assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction; assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; 
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  // (units == meters)
  glm::vec3 a,b,c,d;
  istr >> token >> a.x >> a.y >> a.z; assert (token == "p");
  istr >> token >> b.x >> b.y >> b.z; assert (token == "p");
  istr >> token >> c.x >> c.y >> c.z; assert (token == "p");
  istr >> token >> d.x >> d.y >> d.z; assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  double fabric_weight;
  istr >> token >> fabric_weight; assert (token == "fabric_weight");
  double area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  double mass = area*fabric_weight / double(nx*ny);
  for (int i = 0; i < nx; i++) {
    double x = i/double(nx-1);
    glm::vec3 ab = float(1-x)*a + float(x)*b;
    glm::vec3 dc = float(1-x)*d + float(x)*c;
    for (int j = 0; j < ny; j++) {
      double y = j/double(ny-1);
      ClothParticle &p = getParticle(i,j);
      glm::vec3 abdc = float(1-y)*ab + float(y)*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(glm::vec3(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
      //std::cout<<mass<<std::endl;
    }
  }

  // the fixed particles
  istr >> token;
  while (token=="f") {
    assert (token == "f");
    int i,j;
    double x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(glm::vec3(x,y,z));
    p.setFixed(true);
    istr>>token;
  }
  //where to put the particles and the number of particles in that cell 
  int p_x,p_y,p_num;
  istr>>p_x>>p_y>>p_num;
  //calculate the dimension of one cloth cell 
  dx = (b.x - a.x)/(nx-1);
  //fetch the specified cell  
  ClothParticle &c_w = getParticle(p_x,p_y);
  Cell &cell = c_w.getCell();
  //generate water particles in that cell
  for (int k = 0; k < p_num; ++k)
  {
    glm::vec3 pos = glm::vec3(args->mtrand()*dx,
                              args->mtrand()*dx,
                              0);
    pos+=glm::vec3(p_x*dx,p_y*dx,0);
    FluidParticle *f_p = new FluidParticle();
    f_p->setPosition(pos);
    f_p->setVelocity(glm::vec3(0,0,0));
    cell.addParticle(f_p);
  }
     
  computeBoundingBox();
  
}

// ================================================================================

void Cloth::computeBoundingBox() {
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }
}

// ================================================================================



void Cloth::Animate() {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Compute the forces on each particle, and update the state
  // (position & velocity) of each particle.
  //
  // Also, this is where you'll put the Provot correction for super-elasticity
  //
  // *********************************************************************    
  

  // commented out because an animated bounding box can be weird
  glm::vec3 total_force;
  glm::vec3 aij;
  glm::vec3 vij;
  glm::vec3 position_ij;
  float time_step = float(args->timestep);
  //std::cout<<time_step<<std::endl;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      //get particle at ij
       ClothParticle& pij = getParticle(i,j);
      if (!pij.isFixed())
      {
        //calculate total force
        total_force = compute_structural_force(i,j);
        total_force += compute_shear_force(i,j);
        total_force += compute_bend_force(i,j);
        total_force += -damping * pij.getVelocity();
        //compute acceleration and set acceleration
        aij = total_force / float(pij.getMass());
        aij +=args->gravity;
        pij.setAcceleration(aij);
        //calculate new velocity and set velocity
        vij = pij.getVelocity() + time_step * pij.getAcceleration();
        pij.setVelocity(vij);
        //calculate new position and set position
        position_ij = pij.getPosition() + time_step * pij.getVelocity();
        pij.setPosition(position_ij);
        //do provot structural&shear correction 6 times
      }
      
    }
  }
  for (int k = 0; k < 4; ++k)
  {
    compute_provot_structural();
    compute_provot_shear();
  }
  //check if water collided with cloth 
  check_collision();
  //compute the next velocity for water particles  
  ComputeNewVelocities();
  //compute the next position for water particles
  MoveParticles();
  //reassign water particles
  ReassignParticles();
  //compute position for new water particles of water source
  new_p_water_particles();
  if(count_%3==0)
  GenerateFP();
  count_++;
  if(count_%7==0)
  //check if cells overflow and triger leakage 
  leak_detection();
  new_p_water_particles(); 
}

//calculate total structural force aroung ij
glm::vec3 Cloth::compute_structural_force(int i, int j)
{
  glm::vec3 total_structural_force = glm::vec3(0.0f,0.0f,0.0f);
  if (i-1>=0)
  {
    total_structural_force += force_between_2_point(i,j,i-1,j,0);
  }
  if (i+1<=nx-1)
  {
    total_structural_force += force_between_2_point(i,j,i+1,j,0);
  }
  if (j-1>=0)
  {
    total_structural_force += force_between_2_point(i,j,i,j-1,0);
  }
  if (j+1<=ny-1)
  {
    total_structural_force += force_between_2_point(i,j,i,j+1,0);
  }
  return total_structural_force;

}
//compute total shear force around ij
glm::vec3 Cloth::compute_shear_force(int i, int j)
{
  glm::vec3 total_shear_force = glm::vec3(0.0f,0.0f,0.0f);
  if (i-1>=0 && j-1>=0)
  {
    total_shear_force += force_between_2_point(i,j,i-1,j-1,1);
  }
  if (i-1>=0&&j+1<=ny-1)
  {
    total_shear_force += force_between_2_point(i,j,i-1,j+1,1);
  }
  if (i+1<=nx-1&&j+1<=ny-1)
  {
    total_shear_force += force_between_2_point(i,j,i+1,j+1,1);
  }
  if (i+1<=nx-1&&j-1>=0)
  {
    total_shear_force += force_between_2_point(i,j,i+1,j-1,1);
  }
  return total_shear_force; 
}   
//compute total bend force around ij
glm::vec3 Cloth::compute_bend_force(int i, int j)
{
  glm::vec3 total_bend_force = glm::vec3(0.0f,0.0f,0.0f);
  if (i-2>=0)
  {
    total_bend_force += force_between_2_point(i,j,i-2,j,3);
  }
  if (i+2<=nx-1)
  {
    total_bend_force += force_between_2_point(i,j,i+2,j,3);
  }
  if (j-2>=0)
  {
    total_bend_force += force_between_2_point(i,j,i,j-2,3);
  }
  if (j+2<=ny-1)
  {
    total_bend_force += force_between_2_point(i,j,i,j+2,3);
  }
  return total_bend_force;

}
//this function use to calculate the force from kj to ij
glm::vec3 Cloth::force_between_2_point(int i, int j,int k, int l, int indicator)
{
  glm::vec3 force_;
  ClothParticle  p1;
  ClothParticle  p2;
  glm::vec3 pij;// position of the calculate point
  glm::vec3 pkl;// position of other point
  float l0;// natural lenght of the spring
  glm::vec3 lv;//vector from pij to pkl 
  glm::vec3 v;//normalized l
  p1 = getParticle(i,j);
  p2 = getParticle(k,l);
  pij = p1.getOriginalPosition();
  pkl = p2.getOriginalPosition();
  v = pkl - pij;
  l0 = sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  pij = p1.getPosition();
  pkl = p2.getPosition();
  lv = pkl - pij;
  v = glm::normalize(lv);
  if (indicator == 0)
    force_ = float(k_structural) * (lv - l0 * v);
  else if(indicator == 1)
    force_ = float(k_shear) * (lv - l0 * v);
  else if(indicator == 3)
    force_ = float(k_bend) * (lv - l0 * v);
  return force_;
} 

void Cloth::compute_provot_structural()
{
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      ClothParticle& pij = getParticle(i,j);
      //adjust the position of points around ij
      if (i-1>=0)
      {
        ClothParticle& pkl = getParticle(i-1,j);
        correct_position_2_particle(pij,pkl,0);
      }
      if (i+1<=nx-1)
      {
        ClothParticle& pkl = getParticle(i+1,j);
        correct_position_2_particle(pij,pkl,0);
      }
      if (j-1>=0)
      {
        ClothParticle& pkl = getParticle(i,j-1);
        correct_position_2_particle(pij,pkl,0);
      }
      if (j+1<=ny-1)
      {
        ClothParticle& pkl = getParticle(i,j+1);
        correct_position_2_particle(pij,pkl,0);
      }
    }
  }
}

void Cloth::compute_provot_shear()
{
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      ClothParticle& pij = getParticle(i,j);
      //adjust the position of points around ij
      if (i-1>=0&&j-1>=0)
      {
        ClothParticle& pkl = getParticle(i-1,j-1);
        correct_position_2_particle(pij,pkl,1);
      }
      if (i+1<=nx-1&&j-1>=0)
      {
        ClothParticle& pkl = getParticle(i+1,j-1);
        correct_position_2_particle(pij,pkl,1);
      }
      if (i-1>=0&&j+1<=ny-1)
      {
        ClothParticle& pkl = getParticle(i-1,j+1);
        correct_position_2_particle(pij,pkl,1);
      }
      if (i+1<=nx-1&&j+1<=ny-1)
      {
        ClothParticle& pkl = getParticle(i+1,j+1);
        correct_position_2_particle(pij,pkl,1);
      }
    }
  }

}

void Cloth::correct_position_2_particle(ClothParticle& pij,ClothParticle& pkl,int indicator)
{
  glm::vec3 position_ij;
  glm::vec3 position_kl;
  glm::vec3 l_ij_kl;//current vector from ij to kl
  glm::vec3 org_l_ij_kl;//original vector from ij to kl
  float distance_ij_kl;//current length
  float org_distance_ij_kl;//original length
  float streach_ratio;
  //calculate the current distance between ij and kl
  position_ij = pij.getPosition();
  position_kl = pkl.getPosition();
  l_ij_kl = position_kl - position_ij;
  distance_ij_kl = sqrt(l_ij_kl.x*l_ij_kl.x+l_ij_kl.y*l_ij_kl.y+l_ij_kl.z*l_ij_kl.z);
  //calculate the original distance between ij and kl
  position_ij = pij.getOriginalPosition();
  position_kl = pkl.getOriginalPosition();
  org_l_ij_kl = position_kl - position_ij;
  org_distance_ij_kl = sqrt(org_l_ij_kl.x*org_l_ij_kl.x+org_l_ij_kl.y*org_l_ij_kl.y+org_l_ij_kl.z*org_l_ij_kl.z);
  //calculate streach ratio
  streach_ratio = (distance_ij_kl - org_distance_ij_kl) / org_distance_ij_kl;
  //0 do structural correction
  if (indicator == 0)
  {
    if (streach_ratio >provot_structural_correction)
    {
      //calculate the vector with original length on the current direction
      org_l_ij_kl = glm::normalize(l_ij_kl);
      org_l_ij_kl *=org_distance_ij_kl;
      //calculate the over streached length
      org_l_ij_kl *=(1.0+provot_structural_correction);
      org_l_ij_kl = l_ij_kl - org_l_ij_kl;
      //get the current position of ij and kl
      position_ij = pij.getPosition();
      position_kl = pkl.getPosition();
      if (!pij.isFixed()&&!pkl.isFixed())
      {
        org_l_ij_kl /=2.0;
        pij.setPosition(position_ij + org_l_ij_kl);
        pkl.setPosition(position_kl - org_l_ij_kl);
      }
      else if (!pij.isFixed())
        pij.setPosition(position_ij + org_l_ij_kl);
      else if (!pkl.isFixed())
        pkl.setPosition(position_kl - org_l_ij_kl);
    }
  }
  if (indicator == 1)
  {
    if (streach_ratio > provot_shear_correction)
    {
      //calculate the vector with original length on the current direction
      org_l_ij_kl = glm::normalize(l_ij_kl);
      org_l_ij_kl *=org_distance_ij_kl;
      //calculate the over streached length
      org_l_ij_kl *=(1.0+provot_shear_correction);
      org_l_ij_kl = l_ij_kl - org_l_ij_kl;
      //get the current position of ij and kl
      position_ij = pij.getPosition();
      position_kl = pkl.getPosition();
      if (!pij.isFixed()&&!pkl.isFixed())
      {
        org_l_ij_kl /=2.0;
        pij.setPosition(position_ij + org_l_ij_kl);
        pkl.setPosition(position_kl - org_l_ij_kl);
      }
      else if (!pij.isFixed())
        pij.setPosition(position_ij + org_l_ij_kl);
      else if (!pkl.isFixed())
        pkl.setPosition(position_kl - org_l_ij_kl);
    } 
  }
}
//compute the new velocity for water particles in cloth 
void Cloth::ComputeNewVelocities() {
  double dt = args->timestep;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      ClothParticle &c_p = getParticle(i,j);
      Cell &cell = c_p.getCell();
      std::vector<FluidParticle*> &particles=cell.getParticles();
      for (int k = 0; k < particles.size(); ++k)
      {
      
       FluidParticle *f_p = particles[k];
       glm::vec3 p_v = f_p->getVelocity();
       //calculate the accelration of one water particle
       glm::vec3 acceleration_ = glm::vec3(0,-0.8,0) + AbsorbtionForce(i,j)*0.01;
       //compute new velocity according to acceleration
       p_v += dt * acceleration_;
       //add additional randomness into the velocity
       if (particles.size()<100)
       {
         if(count_%2==0)
         f_p->setVelocity(glm::vec3(-0.2,0,0));
         if(count_%2==1)
         f_p->setVelocity(glm::vec3(0.2,0,0));
       }
       else 
       f_p->setVelocity(p_v);
      }
    }
  }
}
//calculate the position for every water particles in the cloth 
void Cloth::MoveParticles()
{
  double dt = args->timestep;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      ClothParticle &c_p = getParticle(i,j);
      Cell &cell = c_p.getCell();
      std::vector<FluidParticle*> &particles=cell.getParticles();
      for (int k = 0; k < particles.size(); ++k)
      {
       //calculate new velocity by multiplying spped and time step
       FluidParticle *f_p = particles[k];
       glm::vec3 p_v = f_p->getVelocity();
       glm::vec3 pos = f_p->getPosition();
       glm::vec3 new_pos = pos + p_v * dt;
       f_p->setPosition(new_pos);
      }
    }
  }
}
//move paticles from one cell to another cell
void Cloth::ReassignParticles()
{
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      ClothParticle &c_p = getParticle(i,j);
      Cell &cell = c_p.getCell();
      std::vector<FluidParticle*> &particles=cell.getParticles();
      for (int k = 0; k < particles.size(); ++k)
      {
       FluidParticle *f_p = particles[k];
       glm::vec3 pos = f_p->getPosition();
       int i2 = int(floor(pos.x/dx));
       int j2 = int(floor(-pos.y/dx));
       //check if the water paticle is out of current cell
       if (i != i2 || j != j2)
       {
         cell.removeParticle(f_p);//if so remove it
         //check if the water particle is indside the cloth 
         if (i2>=0 && i2<nx)
         {
           if (j2>=0 && j2<ny)
           {
             ClothParticle &c_p2 = getParticle(i2,j2);//add it into another cell
             Cell &cell2 = c_p2.getCell();
             cell2.addParticle(f_p);
           }
         }
       }
      }
    }
  }
}
//this function is used to calculate the absorption force 
glm::vec3 Cloth::AbsorbtionForce(int i,int j)
{
  //fetch the cloth particle and fetch the cell inside it
  ClothParticle &p_ij = getParticle(i,j);
  glm::vec3 position_ij = p_ij.getPosition();
  glm::vec3 total_force = glm::vec3(0,0,0);
  glm::vec3 dir_;
  float dis_;
  int n_p_d;
  //calculate absorption force from surrounding cells  
  if (i>0)
  {
    ClothParticle &p_0j = getParticle(i-1,j);
    glm::vec3 position_0j = p_0j.getPosition();
    dir_= normalize(position_0j - position_ij);//calculate the direction of the force 
    dis_ = length(position_0j - position_ij);//calculate the distance between two cells 
    n_p_d = p_ij.numFluidParticles() - p_0j.numFluidParticles();
    total_force += n_p_d*ab_f*dir_/dis_;//calculate the total force 
    if (j>0)
    {
      ClothParticle &p_00 = getParticle(i-1,j-1);
      glm::vec3 position_00 = p_00.getPosition();
      dir_= normalize(position_00 - position_ij);
      dis_ = length(position_00 - position_ij);
      n_p_d = p_ij.numFluidParticles() - p_00.numFluidParticles();
      total_force += n_p_d*ab_f*dir_/dis_;
    }
    if (j<ny-1)
    {
      ClothParticle &p_01 = getParticle(i-1,j+1);
      glm::vec3 position_01 = p_01.getPosition();
      dir_= normalize(position_01 - position_ij);
      dis_ = length(position_01 - position_ij);
      n_p_d = p_ij.numFluidParticles() - p_01.numFluidParticles();
      total_force += n_p_d*ab_f*dir_/dis_;
    }
    
  }
  if (j>0)
  {
    ClothParticle &p_i0 = getParticle(i,j-1);
    glm::vec3 position_i0 = p_i0.getPosition();
    dir_= normalize(position_i0 - position_ij);
    dis_ = length(position_i0 - position_ij);
    n_p_d = p_ij.numFluidParticles() - p_i0.numFluidParticles();
    total_force += n_p_d*ab_f*dir_/dis_;
  }
  if (j<ny-1)
  {
    ClothParticle &p_i1 = getParticle(i,j+1);
    glm::vec3 position_i1 = p_i1.getPosition();
    dir_= normalize(position_i1 - position_ij);
    dis_ = length(position_i1 - position_ij);
    n_p_d = p_ij.numFluidParticles() - p_i1.numFluidParticles();
    total_force += n_p_d*ab_f*dir_/dis_;
  }
  if (i<nx-1)
  {
    ClothParticle &p_1j = getParticle(i+1,j);
    glm::vec3 position_1j = p_1j.getPosition();
    dir_= normalize(position_1j - position_ij);
    dis_ = length(position_1j - position_ij);
    n_p_d = p_ij.numFluidParticles() - p_1j.numFluidParticles();
    total_force += n_p_d*ab_f*dir_/dis_;
    if (j>0)
    {
      ClothParticle &p_10 = getParticle(i+1,j-1);
      glm::vec3 position_10 = p_10.getPosition();
      dir_= normalize(position_10 - position_ij);
      dis_ = length(position_10 - position_ij);
      n_p_d = p_ij.numFluidParticles() - p_10.numFluidParticles();
      total_force += n_p_d*ab_f*dir_/dis_;
    }
    if (j<ny-1)
    {
      ClothParticle &p_11 = getParticle(i+1,j+1);
      glm::vec3 position_11 = p_11.getPosition();
      dir_= normalize(position_11 - position_ij);
      dis_ = length(position_11 - position_ij);
      n_p_d = p_ij.numFluidParticles() - p_11.numFluidParticles();
      total_force += n_p_d*ab_f*dir_/dis_;
    }
  }
  return total_force;
}

//generate water fluid for water source 
void Cloth::GenerateFP()
{
  for (int i = 0; i < 3; ++i)
  {
    FluidParticle *p = new FluidParticle();
    //add some randomness into the position
    glm::vec3 pos = glm::vec3(args->mtrand()*0.05+0.5,
                                args->mtrand()*0.05+0.6,
                                args->mtrand()*0.1+1.7);
    p->setPosition(pos);
    p->setVelocity(glm::vec3(0,0,-3.5) );//set the initial velocity of the water particles 
    water_particles.push_back(p);

  }
}
//calculate the position for water particles from the water source 
void Cloth::new_p_water_particles()
{
  float dt = args->timestep;
  for (int i = 0; i < water_particles.size(); ++i)
  {
    FluidParticle *p = water_particles[i];
    glm::vec3 pos = p->getPosition();
    glm::vec3 vel = p->getVelocity();
    vel.y+=dt*-9.8;
    pos+=vel*dt;
    p->setPosition(pos);
    p->setVelocity(vel);
  }
}
//this function calculate the area of the triangle orgnized by three points
double get_area(glm::vec3 t1,glm::vec3 t2,glm::vec3 t3) {
  return 0.5*glm::length(glm::cross((t2-t1),(t3-t1)));
}
//this function check if one point is inside a triangle 
bool insideTri(glm::vec3 pos,glm::vec3 t1,glm::vec3 t2,glm::vec3 t3){
  double area = get_area(t1,t2,t3);//get the area of the original triangle 
  double a1 = get_area(t1,pos,t3);//get the area of three small triangle 
  double a2 = get_area(t2,pos,t1);
  double a3 = get_area(t3,pos,t2);
  return fabs(area-(a1+a2+a3))<0.01;//check if the areas are equal
}
//this function calculates the bouncing back velocity 
glm::vec3 bounce_back_v(glm::vec3 normal_, glm::vec3 in_v)
{
  normal_ = normalize(normal_);
  glm::vec3 out_v= in_v - 2*dot(in_v,normal_)*normal_;
  //add somr randomness to the bounce direction to
  //make it more realistic 
  if (count_%4 == 0)
    out_v.z *=-1;
  if (count_%4 == 1)
    out_v.x *=-1; 
  if (count_%4 == 2)
    {
      out_v.x*=-1;
      out_v.z*=-1;
    }
  out_v.y*=0.4;  
  return out_v;
}
//chekc if a water particle hit the cloth 
void Cloth::check_collision()
{
  for (int k = 0; k < water_particles.size(); ++k)
  {
    FluidParticle *p = water_particles[k];
    glm::vec3 pos = p->getPosition();
    glm::vec3 velocity_ = p->getVelocity();
    //only chekc the water paticles in a certain range 
    if (pos.y<=0&&pos.y>=-1&&!p->get_Draw()&&!p->get_Leak())
    {
      //get the two planes orgnized by 4 points 
      for (int i = 0; i < nx-1; i++) {
        for (int j = 0; j < ny-1; j++) {
          ClothParticle &a = getParticle(i,j);
          ClothParticle &b = getParticle(i,j+1);
          ClothParticle &c = getParticle(i+1,j+1);
          ClothParticle &d = getParticle(i+1,j);
          const glm::vec3 &a_pos = a.getPosition();
          const glm::vec3 &b_pos = b.getPosition();
          const glm::vec3 &c_pos = c.getPosition();
          const glm::vec3 &d_pos = d.getPosition();
          //calculate the plane functions for those two planes 
          glm::vec3 v_a_b = (b_pos - a_pos);
          glm::vec3 v_a_d = (d_pos - a_pos);
          glm::vec3 nor_abd = cross(v_a_b,v_a_d);
          glm::vec3 v_c_b = (b_pos - c_pos);
          glm::vec3 v_c_d = (d_pos - c_pos);
          glm::vec3 nor_cdb = cross(v_c_b,v_c_d);
          float m1 = -(nor_abd.x * a_pos.x +
                 nor_abd.y * a_pos.y +
                 nor_abd.z * a_pos.z);
          float m2 = -(nor_cdb.x * c_pos.x +
                 nor_cdb.y * c_pos.y +
                 nor_cdb.z * c_pos.z);
          //substitude the position into the plane and check if the equation equals zero
          if ( fabs(pos.x*nor_abd.x+pos.y*nor_abd.y+pos.z*nor_abd.z+m1)<=error_)
            {  
              //add bounce back to the simulation
              if (count_%2==0)
              {
                p->setVelocity(bounce_back_v(nor_abd,velocity_));
              }
              else
              {  
                p->set_Draw(2);
                float p_a = length(pos - a_pos);
                float p_b = length(pos - b_pos);
                float p_d = length(pos - d_pos);
                //check which cell should absorb this particle
                if (p_a<=p_b&&p_a<=p_d)
                {
                  Cell &cell = a.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
                }
                else if (p_b<=p_a&&p_b<=p_d)
                {
                  Cell &cell = b.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
                }
                else if (p_d<=p_b&&p_d<=p_a)
                {
                  Cell &cell = d.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
                }
              }
            }
            else if (fabs(pos.x*nor_cdb.x+pos.y*nor_cdb.y+pos.z*nor_cdb.z+m2)<=error_ )
            {
              if (count_%2==1)
              {
                p->setVelocity(bounce_back_v(nor_cdb,velocity_));
              }
              else
              {
                p->set_Draw(2);
                float p_c = length(pos - c_pos);
                float p_b = length(pos - b_pos);
                float p_d = length(pos - d_pos); 
                if (p_c<=p_b&&p_c<=p_d)
                {
                  Cell &cell = c.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
                }
                else if (p_b<=p_c&&p_b<=p_d)
                {
                  Cell &cell = b.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
                }
                else if (p_d<=p_b&&p_d<=p_c)
                {
                  Cell &cell = d.getCell();
                  FluidParticle *f_p = new FluidParticle();
                  f_p->setPosition(pos);
                  f_p->setVelocity(glm::vec3(0,0,0));
                  cell.addParticle(f_p);
            }
          }
        }
        }
      }
    }
  }
}
//this function handles leakage
void Cloth::leak_detection()
{
  ClothParticle &a = getParticle(0,0);
  glm::vec3 a_pos = glm::vec3(a.getPosition().x,a.getPosition().y,a.getPosition().z);
  int x =0;
  int y =0;
  //go through every cell and find the position for the lowest cell
  for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          ClothParticle &b = getParticle(i,j);
          const glm::vec3 &b_pos = b.getPosition();
          if (b_pos.y<a_pos.y)
          {
            x=i;
            y=j;
            a_pos = b_pos;
          }
        }
   }
   ClothParticle &c= getParticle(x,y);
   Cell &cell = c.getCell();
   //chekc if there are too much water in that cell
   if (cell.numParticles()>150)
   {
    std::vector<FluidParticle*> &particles=cell.getParticles();
      for (int k = 0; k < 1; ++k)
      {
        //if there are too muc hwater in it 
        FluidParticle *p = particles[k];
        const glm::vec3 m_pos = c.getPosition();
        cell.removeParticle(p);//remove the water particles 
        FluidParticle *f_p = new FluidParticle();
        glm::vec3 p_pos;
        //set up the position of the new water paticles and add some randomness
        if(count_%4==0)
        p_pos = glm::vec3(m_pos.x+args->mtrand()*0.05,m_pos.y,m_pos.z+args->mtrand()*0.05);
        else if(count_%4==1)
        p_pos = glm::vec3(m_pos.x+args->mtrand()*0.05,m_pos.y,m_pos.z-args->mtrand()*0.05);
        else if(count_%4==2)
        p_pos = glm::vec3(m_pos.x-args->mtrand()*0.05,m_pos.y,m_pos.z+args->mtrand()*0.05);
        else if(count_%4==3)
        p_pos = glm::vec3(m_pos.x-args->mtrand()*0.05,m_pos.y,m_pos.z-args->mtrand()*0.05);
        f_p->setPosition(p_pos);
        f_p->setVelocity(glm::vec3(0,0,0));
        f_p->set_Leak(1);
        water_particles.push_back(f_p);//add the new water paticles in to water source 
        //so the particles can be visualized 
      }

     
   }
          
}













