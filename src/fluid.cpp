#include "glCanvas.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

#include "fluid.h"
#include "argparser.h"
#include "boundingbox.h"
#include "marching_cubes.h"
#include "utils.h"

#define BETA_0 1.7
#define EPSILON 0.0001

// ==============================================================
// ==============================================================
// CONSTRUCTOR
// ==============================================================
// ==============================================================

Fluid::Fluid(ArgParser *_args) {
  args = _args;
  Load();
  marchingCubes = new MarchingCubes(nx+1,ny+1,nz+1,dx,dy,dz);
  SetEmptySurfaceFull();
}

Fluid::~Fluid() { 
  delete [] cells; 
  delete marchingCubes; 
  cleanupVBOs(); 
}

// ==============================================================

void Fluid::Load() {    

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->fluid_file).c_str());
  assert (istr.good());
  std::string token, token2, token3;

  // load in the grid size & dimensions
  istr >> token >> nx >> ny >> nz;  assert (token=="grid");
  assert (nx > 0 && ny > 0 && nz > 0);
  istr >> token >> dx >> dy >> dz; assert (token=="cell_dimensions");
  cells = new Cell[(nx+2)*(ny+2)*(nz+2)];

  // simulation parameters
  istr >> token >> token2;  assert (token=="flow");
  if (token2 == "compressible") compressible = true;
  else { assert (token2 == "incompressible"); compressible = false; }
  istr >> token >> token2;  assert (token=="xy_boundary");
  if (token2 == "free_slip") xy_free_slip = true;
  else { assert  (token2 == "no_slip"); xy_free_slip = false; }
  istr >> token >> token2;  assert (token=="yz_boundary");
  if (token2 == "free_slip") yz_free_slip = true;
  else { assert  (token2 == "no_slip"); yz_free_slip = false; }
  istr >> token >> token2;  assert (token=="zx_boundary");
  if (token2 == "free_slip") zx_free_slip = true;
  else { assert  (token2 == "no_slip"); zx_free_slip = false; }
  istr >> token >> viscosity;  assert (token=="viscosity");
  double gravity;
  istr >> token >> gravity;  assert (token=="gravity");
  args->gravity = glm::vec3(0,-9.8,0) * float(gravity);
  
  // initialize marker particles 
  istr >> token >> token2 >> token3;  assert (token=="initial_particles");
  istr >> token >> density;  assert (token=="density");
  GenerateParticles(token2,token3);

  // initialize velocities
  istr >> token >> token2;  assert (token=="initial_velocity");
  if (token2 == "zero") {
    // default is zero
  } else {
    assert (token2 == "random");
    int i,j,k;
    double max_dim = std::max(dx,std::max(dy,dz));
    for (i = -1; i <= nx; i++) {
      for (j = -1; j <= ny; j++) {
        for (k = -1; k <= nz; k++) {
          getCell(i,j,k)->set_u_plus((2*args->mtrand()-1)*max_dim);
	  getCell(i,j,k)->set_v_plus((2*args->mtrand()-1)*max_dim);
	  getCell(i,j,k)->set_w_plus((2*args->mtrand()-1)*max_dim);
        }
      }
    }
  }
  // read in custom velocities
  while(istr >> token) {
    int i,j,k;
    double velocity;
    assert (token == "u" || token == "v" || token == "w");
    istr >> i >> j >> k >> velocity;
    assert(i >= 0 && i < nx);
    assert(j >= 0 && j < ny);
    assert(k >= 0 && k < nz);
    if      (token == "u") getCell(i,j,k)->set_u_plus(velocity);
    else if (token == "v") getCell(i,j,k)->set_v_plus(velocity);
    else if (token == "w") getCell(i,j,k)->set_w_plus(velocity);
    else assert(0);
  }
  SetBoundaryVelocities();
}

// ==============================================================

bool Fluid::inShape(glm::vec3 &pos, const std::string &shape) {
  // return true if this point is inside the "shape"
  // defined procedurally (using an implicit surface)
  if (shape == "everywhere") {
    return true;
  } else if (shape == "left") {
    // a blob of particles on the lower left (for the dam)
    return (pos.x < 0.2*nx*dx && pos.y < 0.5*ny*dy);
  } else if (shape == "drop") {
    // a shallow pool of particles on the bottom
    double h = ny*dy/6.0;
    if (pos.y < 2*h) return true;
    // and a sphere of particles above
    glm::vec3 center = glm::vec3(nx*dx*0.5, 5*h,nz*dz*0.5);
    double length = glm::length(center-pos);
    if (length < 0.8*h) return true;
    return false;
  } else {
    std::cout << "unknown shape: " << shape << std::endl;
    exit(0);
  }
}

// ==============================================================

void Fluid::GenerateParticles(const std::string &shape, const std::string &placement) {
  // create a set of points according to the "placement" token,
  // then check whether they are inside of the "shape"
  if (placement == "uniform") {
    int dens = (int)pow(density,0.334);
    assert (dens*dens*dens == density);
    // the uniform grid spacing
    double spacing = 1/double(dens);
    for (double x = 0.5*spacing*dx; x < nx*dx; x += spacing*dx) {
      for (double y = 0.5*spacing*dy; y < ny*dy; y += spacing*dy) {
        for (double z = 0.5*spacing*dz; z < nz*dz; z += spacing*dz) {
          glm::vec3 pos = glm::vec3(x,y,z);
          if (inShape(pos,shape)) {
            Cell *cell = getCell(int(x/dx),int(y/dy),int(z/dz));
            FluidParticle *p = new FluidParticle();
            p->setPosition(pos);
            cell->addParticle(p);
          }
        }
      }
    }
  } else {
    assert (placement == "random");
    // note: we don't necessarily have the same number of particles in each cell
    for (int n = 0; n < nx*ny*nz*density; n++) {
      glm::vec3 pos = glm::vec3(args->mtrand()*nx*dx,
                                args->mtrand()*ny*dy,
                                args->mtrand()*nz*dz);
      if (inShape(pos,shape)) {      
        Cell *cell = getCell(int(pos.x/dx),int(pos.y/dy),int(pos.z/dz));
        FluidParticle *p = new FluidParticle();
        p->setPosition(pos);
        cell->addParticle(p);
      }
    }
  }
}

// ==============================================================
// ==============================================================
// ANIMATION
// ==============================================================
// ==============================================================

void Fluid::Animate() {

  // the animation manager:  this is what gets done each timestep!

  ComputeNewVelocities();
  SetBoundaryVelocities();
  
  // compressible / incompressible flow
  if (compressible == false) {
    for (int iters = 0; iters < 20; iters++) {
      double max_divergence = AdjustForIncompressibility();
      SetBoundaryVelocities();
      //std::cout<<max_divergence<<std::endl;
      if (max_divergence < EPSILON) break;
    }
  }
 
  UpdatePressures();
  CopyVelocities();

  // advanced the particles through the fluid
  MoveParticles();
  ReassignParticles();
  SetEmptySurfaceFull();

  setupVBOs();
}

// ==============================================================

void Fluid::ComputeNewVelocities() {
  double dt = args->timestep;
  int i,j,k;

  // using the formulas from Foster & Metaxas

  for (i = 0; i < nx-1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        double new_u_plus =
          get_u_plus(i,j,k) +            
          dt * ((1/dx) * (square(get_u_avg(i,j,k)) - square(get_u_avg(i+1,j,k))) +
                (1/dy) * (get_uv_plus(i,j-1,k) - get_uv_plus(i,j,k)) + 
                (1/dz) * (get_uw_plus(i,j,k-1) - get_uw_plus(i,j,k)) +
                args->gravity.x +
                (1/dx) * (getPressure(i,j,k)-getPressure(i+1,j,k)) +
                (viscosity/square(dx)) * (get_u_plus(i+1,j  ,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_u_plus(i  ,j+1,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_u_plus(i  ,j  ,k+1) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j  ,k-1)) );
        cell->set_new_u_plus(new_u_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny-1; j++) {
      for (k = 0; k < nz; k++) {	
        Cell *cell = getCell(i,j,k);
        double new_v_plus =
          get_v_plus(i,j,k) +
          dt * ((1/dx) * (get_uv_plus(i-1,j,k) - get_uv_plus(i,j,k)) +
                (1/dy) * (square(get_v_avg(i,j,k)) - square(get_v_avg(i,j+1,k))) +
                (1/dz) * (get_vw_plus(i,j,k-1) - get_vw_plus(i,j,k)) +
                args->gravity.y +
                (1/dy) * (getPressure(i,j,k)-getPressure(i,j+1,k)) +
                (viscosity/square(dx)) * (get_v_plus(i+1,j  ,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_v_plus(i  ,j+1,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_v_plus(i  ,j  ,k+1) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j  ,k-1)) );
        cell->set_new_v_plus(new_v_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz-1; k++) {
        Cell *cell = getCell(i,j,k);
        double new_w_plus =
          get_w_plus(i,j,k) +
          dt * ((1/dx) * (get_uw_plus(i-1,j,k) - get_uw_plus(i,j,k)) +
                (1/dy) * (get_vw_plus(i,j-1,k) - get_vw_plus(i,j,k)) +
                (1/dz) * (square(get_w_avg(i,j,k)) - square(get_w_avg(i,j,k+1))) +
                args->gravity.z +
                (1/dz) * (getPressure(i,j,k)-getPressure(i,j,k+1)) +
                (viscosity/square(dx)) * (get_w_plus(i+1,j  ,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_w_plus(i  ,j+1,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_w_plus(i  ,j  ,k+1) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j  ,k-1)) );
        cell->set_new_w_plus(new_w_plus);
      }
    }
  }
}

// ==============================================================

void Fluid::SetBoundaryVelocities() {

  // zero out flow perpendicular to the boundaries (no sources or sinks)
  for (int j = -1; j <= ny; j++) {
    for (int k = -1; k <= nz; k++) {
      getCell(-1  ,j,k)->set_u_plus(0);
      getCell(nx-1,j,k)->set_u_plus(0);
      getCell(nx  ,j,k)->set_u_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1  ,k)->set_v_plus(0);
      getCell(i,ny-1,k)->set_v_plus(0);
      getCell(i,ny  ,k)->set_v_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1  )->set_w_plus(0);
      getCell(i,j,nz-1)->set_w_plus(0);
      getCell(i,j,nz  )->set_w_plus(0);
    }
  }

  // free slip or no slip boundaries (friction with boundary)
  double xy_sign = (xy_free_slip) ? 1 : -1;
  double yz_sign = (yz_free_slip) ? 1 : -1;
  double zx_sign = (zx_free_slip) ? 1 : -1;
  for (int i = 0; i < nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1)->set_u_plus(xy_sign*getCell(i,j,0)->get_u_plus());
      getCell(i,j,nz)->set_u_plus(xy_sign*getCell(i,j,nz-1)->get_u_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1,k)->set_u_plus(zx_sign*getCell(i,0,k)->get_u_plus());
      getCell(i,ny,k)->set_u_plus(zx_sign*getCell(i,ny-1,k)->get_u_plus());
    }
  }
  for (int j = 0; j < ny; j++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,j,-1)->set_v_plus(xy_sign*getCell(i,j,0)->get_v_plus());
      getCell(i,j,nz)->set_v_plus(xy_sign*getCell(i,j,nz-1)->get_v_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(-1,j,k)->set_v_plus(yz_sign*getCell(0,j,k)->get_v_plus());
      getCell(nx,j,k)->set_v_plus(yz_sign*getCell(nx-1,j,k)->get_v_plus());
    }
  }
  for (int k = 0; k < nz; k++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,-1,k)->set_w_plus(zx_sign*getCell(i,0,k)->get_w_plus());
      getCell(i,ny,k)->set_w_plus(zx_sign*getCell(i,ny-1,k)->get_w_plus());
    }
    for (int j = -1; j <= ny; j++) {
      getCell(-1,j,k)->set_w_plus(yz_sign*getCell(0,j,k)->get_w_plus());
      getCell(nx,j,k)->set_w_plus(yz_sign*getCell(nx-1,j,k)->get_w_plus());
    }
  }
}

// ==============================================================

void Fluid::EmptyVelocities(int i, int j, int k) {
  Cell *c = getCell(i,j,k);
  if (c->getStatus() != CELL_EMPTY) return;
  Cell *ciplus = getCell(i+1,j,k);
  Cell *cjplus = getCell(i,j+1,k);
  Cell *ckplus = getCell(i,j,k+1);
  if (ciplus->getStatus() == CELL_EMPTY)
    c->set_new_u_plus(0);
  if (cjplus->getStatus() == CELL_EMPTY)
    c->set_new_v_plus(0);
  if (ckplus->getStatus() == CELL_EMPTY)
    c->set_new_w_plus(0);
}


// move to new timestep
void Fluid::CopyVelocities() {
  double dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
	Cell *c = getCell(i,j,k);
	
	EmptyVelocities(i,j,k);
	c->copyVelocity();
	if (fabs(c->get_u_plus()) > 0.5*dx/dt ||
	    fabs(c->get_v_plus()) > 0.5*dy/dt ||
	    fabs(c->get_w_plus()) > 0.5*dz/dt) {
	  // velocity has exceeded reasonable threshhold
	  std::cout << "velocity has exceeded reasonable threshhold, stopping animation" << std::endl;
	  args->animate=false;
	}
      }
    }
  }
}

// ==============================================================

double Fluid::AdjustForIncompressibility() {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // This is not a complete implementation of the Marker and Cell (MAC) method.
  // Additional boundary velocities should be equalized as described in the references
  // depending on whether the boundaries are free-slip or no-slip.
  //
  // Also play around with compressible flow!
  //
  // *********************************************************************    
  //std::cout<<"break"<<std::endl;
  double divergence;//store the net inflow/outflow
  double current_divergence;//store the calculate divergence
  double max_divergence=0;//store the maximum divergence
  double old_flow;
  double dt_=args->timestep;
  glm::vec3 gravity_=args->gravity ;
  double v_y_n;
  double v_y_s;
  double v_x_w;
  double v_x_e;
  double v_y_c;
  double v_x_c;
  int num_faces;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell_c = getCell(i,j,k);
        if (cell_c->getStatus()!=CELL_SURFACE)
        {
        divergence=0;
        num_faces=0;
        //calculate inflow/outflow
        //calculate hw many faces the fluid can flow
        if (i!=0)
        {
          divergence+=get_new_u_plus(i-1,j,k);
          num_faces+=1;
          //std::cout<<get_new_u_plus(i-1,j,k)<<" inflow u"<<std::endl;
        }
        if (i!=nx-1)
        {
          divergence-=get_new_u_plus(i,j,k);
          num_faces+=1;
          //std::cout<<get_new_u_plus(i,j,k)<<" outflow u"<<std::endl;
        }
        if (j!=0)
        {
          divergence+=get_new_v_plus(i,j-1,k);
          num_faces+=1;
          //std::cout<<get_new_v_plus(i,j-1,k)<<" inflow v"<<std::endl;
        }
        if (j!=ny-1)
        {
          divergence-=get_new_v_plus(i,j,k);
          num_faces+=1;
          //std::cout<<get_new_v_plus(i,j,k)<<" outflow v"<<std::endl;
        }
        if (k!=0)
        {
          divergence+=get_new_w_plus(i,j,k-1);
          num_faces+=1;
          //std::cout<<get_new_w_plus(i,j,k-1)<<" inflow w"<<std::endl;
        }
        if (k!=nz-1)
        {
          divergence-=get_new_w_plus(i,j,k);
          num_faces+=1;
          //std::cout<<get_new_w_plus(i,j,k)<<" outflow w"<<std::endl;
        }
        //std::cout<<divergence<<" divergence"<<std::endl;
        divergence/=num_faces;
        //calculate the new velocity for inflow/outflow
        if (i!=0)
        {
          old_flow = get_new_u_plus(i-1,j,k);
          old_flow -=divergence;
          set_new_u_plus(i-1,j,k,old_flow);
          //std::cout<<old_flow<<" new inflow u"<<std::endl;
        }
        if (i!=nx-1)
        {
          old_flow = get_new_u_plus(i,j,k);
          old_flow +=divergence;
          set_new_u_plus(i,j,k,old_flow);
          //std::cout<<old_flow<<" new outflow u"<<std::endl;
        }
        if (j!=0)
        {
          old_flow = get_new_v_plus(i,j-1,k);
          old_flow -=divergence;
          set_new_v_plus(i,j-1,k,old_flow);
          //std::cout<<old_flow<<" new inflow v"<<std::endl;
        }
        if (j!=ny-1)
        {
          old_flow = get_new_v_plus(i,j,k);
          old_flow +=divergence;
          set_new_v_plus(i,j,k,old_flow);
          //std::cout<<old_flow<<" new outflow v"<<std::endl;
        }
        if (k!=0)
        {
          old_flow = get_new_w_plus(i,j,k-1);
          old_flow -=divergence;
          set_new_w_plus(i,j,k-1,old_flow);
          //std::cout<<old_flow<<" new inflow w"<<std::endl;
        }
        if (k!=nz-1)
        {
          old_flow =get_new_w_plus(i,j,k);
          old_flow +=divergence;
          set_new_w_plus(i,j,k,old_flow);
          //std::cout<<old_flow<<" new outflow w"<<std::endl;
        }
        //calculate the current divergence 
        }
      }
    }
  }


  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell_c = getCell(i,j,k);
        Cell *cell_n = getCell(i,j+1,k);
        Cell *cell_w = getCell(i-1,j,k);
        Cell *cell_e = getCell(i+1,j,k);
        Cell *cell_s = getCell(i,j-1,k);
        if (cell_c->getStatus()==CELL_SURFACE)
        {

          

          if ((cell_n->getStatus()==CELL_EMPTY)&&cell_s->getStatus()==CELL_EMPTY)
          {
            v_y_n = get_v_plus(i,j+1,k);
            v_y_n += (gravity_.y * dt_);
            set_new_v_plus(i,j+1,k,-v_y_n);
            v_y_s = get_v_plus(i,j-1,k);
            v_y_s += (gravity_.y * dt_);
            set_new_v_plus(i,j-1,k,-v_y_s);
            
          }
          else if ((cell_n->getStatus()==CELL_EMPTY)&&cell_s->getStatus()!=CELL_EMPTY)
          {
            v_y_s = get_new_v_plus(i,j+1,k);
            set_new_v_plus(i,j-1,k,-v_y_s);
          }
          else if (cell_n->getStatus()!=CELL_EMPTY&&cell_s->getStatus()==CELL_EMPTY)
          {
            v_y_n = get_new_v_plus(i,j-1,k);
            set_new_v_plus(i,j-1,k,-v_y_n);
          }
          


          if (cell_w->getStatus()==CELL_EMPTY&&cell_e->getStatus()==CELL_EMPTY)
          {
            v_x_e = get_u_plus(i-1,j,k);
            set_new_u_plus(i-1,j,k,-v_x_e);
            v_x_w = get_u_plus(i+1,j,k);
            set_new_u_plus(i+1,j,k,-v_x_w);
          }
          else if (cell_w->getStatus()==CELL_EMPTY&&cell_e->getStatus()!=CELL_EMPTY)
          {
            v_x_e = get_new_u_plus(i+1,j,k);
            set_new_u_plus(i-1,j,k,-v_x_e);
          }
          else if (cell_w->getStatus()!=CELL_EMPTY&&cell_e->getStatus()==CELL_EMPTY)
          {
            v_x_w = get_new_u_plus(i-1,j,k);
            set_new_u_plus(i+1,j,k,-v_x_w);
          }
          





        }
      }
    }
  }
  


  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        current_divergence =
        - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
        (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
        (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
        if (current_divergence>max_divergence)
          max_divergence = current_divergence;

      }
    }
  }
  // return the maximum divergence
  // (will iterate for specified # of iterations or until divergence is near zero)
  return max_divergence;

}

// ==============================================================

void Fluid::UpdatePressures() {
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
	Cell *c = getCell(i,j,k);
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
	  // compute divergence and increment/decrement pressure
	  double pressure = c->getPressure();
	  double divergence = 
	    - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
		(1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
		(1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
	  double dt = args->timestep;
	  double beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
	  double dp = beta*divergence;
	  c->setPressure(pressure + dp);
	} else {
	  // zero out boundary cells (just in case)
	  c->setPressure(0);
	}

	
	// =======================================
	// HACK? From Foster 2001 paper?
	// zero out empty cells
	if (c->getStatus() == CELL_EMPTY) {
	  c->setPressure(0);
	}
	// ========================================

      }
    }
  }
}

// ==============================================================

void Fluid::MoveParticles() {
  double dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          glm::vec3 pos = p->getPosition();
          glm::vec3 vel = getInterpolatedVelocity(pos);
          glm::vec3 pos2 = pos + float(dt)*vel;
#if 0
          // euler integration
          p->setPosition(pos2);
#else
          // trapezoid integration
          glm::vec3 vel2 = getInterpolatedVelocity(pos2);
          glm::vec3 pos3 = pos + float(0.5*dt)*(vel+vel2);
          p->setPosition(pos3);
#endif
        }
      }
    }
  }
}

// ==============================================================

void Fluid::ReassignParticles() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          glm::vec3 pos = p->getPosition();
          int i2 = (int)std::min(double(nx-1),std::max(0.0,floor(pos.x/dx)));
          int j2 = (int)std::min(double(ny-1),std::max(0.0,floor(pos.y/dy)));
          int k2 = (int)std::min(double(nz-1),std::max(0.0,floor(pos.z/dz)));
          // if the particle has crossed one of the cell faces 
          // assign it to the new cell
          if (i != i2 || j != j2 || k != k2) {
            cell->removeParticle(p);
            getCell(i2,j2,k2)->addParticle(p);
          } 
        }
      }
    }
  }
}

// ==============================================================

void Fluid::SetEmptySurfaceFull() {
  int i,j,k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->numParticles() == 0)
          cell->setStatus(CELL_EMPTY);
        else 
          cell->setStatus(CELL_FULL);
      }
    }
  }

  // pick out the boundary cells
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->getStatus() == CELL_FULL &&
            (getCell(i-1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i+1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j-1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j+1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k-1)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k+1)->getStatus() == CELL_EMPTY)) {
          cell->setStatus(CELL_SURFACE);
        }
      }
    }
  }
}

// ==============================================================

glm::vec3 Fluid::getInterpolatedVelocity(const glm::vec3 &pos) const {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Here is the naive velocity interpolation.
  // (use a simple average of the face velocity throughout the cell)
  // Do it right, as described in the papers.
  //
  int i = int(floor(pos.x/dx)); if (i < 0) i = 0; if (i >= nx) i = nx-1;
  int j = int(floor(pos.y/dy)); if (j < 0) j = 0; if (j >= ny) j = ny-1;
  int k = int(floor(pos.z/dz)); if (k < 0) k = 0; if (k >= nz) k = nz-1;
  glm::vec3 interpolate_v;
  /*int m = int(floor((pos.y+dy/2)/dy)); if (m < 0) m = 0; if (m >= ny) m = ny-1;
  //std::cout<<"j is "<<j<<" m is "<<m<<std::endl;
  glm::vec3 pos_c_p = glm::vec3((i+0.5)*dx,m*dy,k*dz);
  //calculate the center\ point for the closet four points
  double A0;
  double A1;
  double A2;
  double A3;
  double area_ = dx * dy;// calculate the area of a single square 
  //start calculate interpolate v for x 
  A0= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2));
  A1= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2));
  A2= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2));
  A3= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2));
  interpolate_v.x = A0/area_*get_u_plus(i,m,k)+
                    A1/area_*get_u_plus(i-1,m,k)+
                    A2/area_*get_u_plus(i-1,m-1,k)+
                    A3/area_*get_u_plus(i,m-1,k);
  //start calculate interpolate v for y
  m = int(floor((pos.x+dx/2)/dx)); if (m < 0) m = 0; if (m >= nx) m = nx-1;
  pos_c_p = glm::vec3(m*dx,(j+0.5)*dy,k*dz);
  A0= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2));
  A1= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2));
  A2= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2));
  A3= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2));
  interpolate_v.y = A0/area_*get_v_plus(m,j,k)+
                    A1/area_*get_v_plus(m-1,j,k)+
                    A2/area_*get_v_plus(m-1,j-1,k)+
                    A3/area_*get_v_plus(m,j-1,k);

  return glm::vec3(interpolate_v.x,interpolate_v.y,get_w_avg(i,j,k));*/
  ////
  //calculate the center point for cube
  int c_y = int(floor((pos.y+dy/2)/dy)); if (c_y < 0) c_y = 0; if (c_y >= ny) c_y = ny-1;
  int c_z = int(floor((pos.z+dz/2)/dz)); if (c_z < 0) c_z = 0; if (c_z >= nz) c_z = nz-1;
  glm::vec3 pos_c_p = glm::vec3((i+0.5)*dx,c_y*dy,c_z*dz);
  double A0;
  double A1;
  double A2;
  double A3;
  double A4;
  double A5;
  double A6;
  double A7;
  double volume_ = dx*dy*dz;
  //start calculate interpolate velocity on x
  A0= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A1= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A2= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A3= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A4= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A5= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A6= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A7= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  interpolate_v.x = A0/volume_*get_u_plus(i,c_y,c_z-1)+
                    A1/volume_*get_u_plus(i-1,c_y,c_z-1)+
                    A2/volume_*get_u_plus(i-1,c_y,c_z)+
                    A3/volume_*get_u_plus(i,c_y,c_z)+
                    A4/volume_*get_u_plus(i,c_y-1,c_z-1)+
                    A5/volume_*get_u_plus(i-1,c_y-1,c_z-1)+
                    A6/volume_*get_u_plus(i-1,c_y-1,c_z)+
                    A7/volume_*get_u_plus(i,c_y-1,c_z);
  //start to calculate interpolate velocity on y 
  int c_x = int(floor((pos.x+dx/2)/dx)); if (c_x < 0) c_x = 0; if (c_x >= nx) c_x = nx-1;
  c_z = int(floor((pos.z+dz/2)/dz)); if (c_z < 0) c_z = 0; if (c_z >= nz) c_z = nz-1;
  pos_c_p = glm::vec3(c_x*dx,(j+0.5)*dy,c_z*dz);
  A0= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A1= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A2= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A3= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A4= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A5= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A6= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A7= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  interpolate_v.y = A0/volume_*get_v_plus(c_x,j,c_z-1)+
                    A1/volume_*get_v_plus(c_x-1,j,c_z-1)+
                    A2/volume_*get_v_plus(c_x-1,j,c_z)+
                    A3/volume_*get_v_plus(c_x,j,c_z)+
                    A4/volume_*get_v_plus(c_x,j-1,c_z-1)+
                    A5/volume_*get_v_plus(c_x-1,j-1,c_z-1)+
                    A6/volume_*get_v_plus(c_x-1,j-1,c_z)+
                    A7/volume_*get_v_plus(c_x,j-1,c_z);
  //start to calculate interpolate velocity on z
  c_x = int(floor((pos.x+dx/2)/dx)); if (c_x < 0) c_x = 0; if (c_x >= nx) c_x = nx-1;
  c_y = int(floor((pos.y+dy/2)/dy)); if (c_y < 0) c_y = 0; if (c_y >= ny) c_y = ny-1;
  pos_c_p = glm::vec3(c_x*dx,c_y*dy,(k+0.5)*dz);
  A0= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A1= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A2= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A3= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y+dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A4= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A5= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z-dz/2));
  A6= std::abs(pos_c_p.x - (pos.x-dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  A7= std::abs(pos_c_p.x - (pos.x+dx/2)) * std::abs(pos_c_p.y - (pos.y-dy/2)) *
  std::abs(pos_c_p.z - (pos.z+dz/2));
  interpolate_v.z = A0/volume_*get_w_plus(c_x,c_y,k-1)+
                    A1/volume_*get_w_plus(c_x-1,c_y,k-1)+
                    A2/volume_*get_w_plus(c_x-1,c_y,k)+
                    A3/volume_*get_w_plus(c_x,c_y,k)+
                    A4/volume_*get_w_plus(c_x,c_y-1,k-1)+
                    A5/volume_*get_w_plus(c_x-1,c_y-1,k-1)+
                    A6/volume_*get_w_plus(c_x-1,c_y-1,k)+
                    A7/volume_*get_w_plus(c_x,c_y-1,k);
  // *********************************************************************  
  return glm::vec3(interpolate_v.x,interpolate_v.y,interpolate_v.z);

}

