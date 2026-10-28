// Manufactured SI rectangular foil using MFEM's unchanged Joule operator.
#include "joule_solver.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace mfem { namespace electromagnetics {
int SOLVER_PRINT_LEVEL = 0;
int STATIC_COND = 0;
real_t drive_voltage = 150.0;
real_t drive_length = 0.01;
real_t drive_rise = 1.0e-7;
void edot_bc(const Vector &, Vector &E) { E = 0.0; }
void e_exact(const Vector &, real_t, Vector &E) { E = 0.0; }
void b_exact(const Vector &, real_t, Vector &B) { B = 0.0; }
real_t t_exact(const Vector &) { return 0.0; }
real_t p_bc(const Vector &x, real_t t)
{
   const real_t ramp = drive_rise > 0 ? std::min(real_t(1), std::max(real_t(0), t/drive_rise)) : 1;
   return drive_voltage*ramp*(x[2]/drive_length - 0.5);
}
} }

using namespace mfem;
using namespace mfem::electromagnetics;

class AuditedJoule : public MagneticDiffusionEOperator
{
public:
   using MagneticDiffusionEOperator::MagneticDiffusionEOperator;
   void SolverMetrics(int &iterations, double &residual) const
   {
      iterations = 0; residual = 0;
      for (auto *solver : {pcg_a0, pcg_a1, pcg_a2, pcg_m3})
      {
         if (!solver) { continue; }
         int count = 0; real_t norm = 0;
         solver->GetNumIterations(count); solver->GetFinalResidualNorm(norm);
         iterations = std::max(iterations, count); residual = std::max(residual, double(norm));
         MFEM_VERIFY(std::isfinite(norm) && count < SOLVER_MAX_IT, "Joule linear solve failed its finite/iteration guard");
      }
   }
};

struct Settings
{
   int nx=32, ny=4, nz=2, steps=40;
   double width=.005, thickness=.0002, length=.01, final_time=2e-6;
   double sigma=3e6, capacity=2.8e6, conductivity=16, mu=4e-7*M_PI;
};

void FaceMeans(ParMesh &mesh, ParGridFunction &temperature, const Settings &s,
               std::vector<double> &top, std::vector<double> &bottom)
{
   const int n=s.nx*s.nz;
   std::vector<double> local_sum(2*n,0), local_area(2*n,0), sums(2*n), areas(2*n);
   Vector xyz(3);
   for(int b=0;b<mesh.GetNBE();++b)
   {
      const int attribute=mesh.GetBdrAttribute(b);
      if(attribute!=3 && attribute!=4) { continue; }
      auto *face=mesh.GetBdrFaceTransformations(b);
      MFEM_VERIFY(face && face->Elem1No>=0, "Missing material-side boundary transformation");
      const auto &rule=IntRules.Get(mesh.GetBdrElementBaseGeometry(b),4);
      for(int q=0;q<rule.GetNPoints();++q)
      {
         const auto &ip=rule.IntPoint(q);
         face->SetAllIntPoints(&ip); face->Transform(ip,xyz);
         const int ix=std::clamp(int(std::floor((xyz[0]+.5*s.width)/s.width*s.nx)),0,s.nx-1);
         const int iz=std::clamp(int(std::floor(xyz[2]/s.length*s.nz)),0,s.nz-1);
         const int index=(attribute==4?0:n)+iz*s.nx+ix;
         const double area=ip.weight*face->Weight();
         const double value=temperature.GetValue(face->Elem1No,face->GetElement1IntPoint());
         MFEM_VERIFY(std::isfinite(value) && std::isfinite(area) && area>0,"Invalid boundary trace");
         local_sum[index]+=value*area; local_area[index]+=area;
      }
   }
   MPI_Allreduce(local_sum.data(),sums.data(),2*n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(local_area.data(),areas.data(),2*n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   top.resize(n); bottom.resize(n);
   const double target=s.width/s.nx*s.length/s.nz;
   for(int i=0;i<2*n;++i)
   {
      MFEM_VERIFY(std::abs(areas[i]-target)<=1e-10*target,"Material-face patches do not completely cover expected area");
      (i<n?top[i]:bottom[i-n])=sums[i]/areas[i];
   }
}

int main(int argc,char **argv)
{
   Mpi::Init(argc,argv); Hypre::Init();
   Settings s;
   const char *output=".";
   OptionsParser options(argc,argv);
   options.AddOption(&s.nx,"-nx","--width-cells","Width cells");
   options.AddOption(&s.ny,"-ny","--thickness-cells","Thickness cells");
   options.AddOption(&s.nz,"-nz","--current-cells","Current-direction cells");
   options.AddOption(&s.steps,"-ns","--steps","Uniform timestep count");
   options.AddOption(&s.width,"-w","--width-m","Foil width in m");
   options.AddOption(&s.thickness,"-h","--thickness-m","Foil thickness in m");
   options.AddOption(&s.length,"-l","--length-m","Foil current length in m");
   options.AddOption(&s.final_time,"-tf","--final-time-s","Final time in s");
   options.AddOption(&s.sigma,"-sig","--sigma-s-m","Constant conductivity in S/m");
   options.AddOption(&s.capacity,"-cv","--capacity-j-m3-k","Volumetric heat capacity");
   options.AddOption(&s.conductivity,"-k","--conductivity-w-m-k","Thermal conductivity");
   options.AddOption(&drive_voltage,"-v","--voltage-v","Terminal voltage difference");
   options.AddOption(&drive_rise,"-tr","--rise-time-s","Linear voltage rise time");
   options.AddOption(&output,"-out","--output-dir","Existing output directory");
   options.Parse();
   if(!options.Good()) { if(Mpi::Root()) options.PrintUsage(std::cerr); return 2; }
   MFEM_VERIFY(s.nx>0 && s.ny>0 && s.nz>0 && s.steps>0,"Positive resolution counts required");
   for(double value:{s.width,s.thickness,s.length,s.final_time,s.sigma,s.capacity,s.conductivity})
      { MFEM_VERIFY(std::isfinite(value)&&value>0,"Positive finite SI parameters required"); }
   drive_length=s.length;
   Mesh serial=Mesh::MakeCartesian3D(s.nx,s.ny,s.nz,Element::HEXAHEDRON,s.width,s.thickness,s.length);
   for(int i=0;i<serial.GetNV();++i)
   { auto *v=serial.GetVertex(i); v[0]-=.5*s.width; v[1]-=.5*s.thickness; }
   for(int b=0;b<serial.GetNBE();++b)
   {
      Array<int> vertices; serial.GetBdrElementVertices(b,vertices); Vector center(3);center=0.0;
      for(int j=0;j<vertices.Size();++j)
         for(int k=0;k<3;++k) center[k]+=serial.GetVertex(vertices[j])[k]/vertices.Size();
      int attribute=5;
      if(std::abs(center[2])<1e-12) attribute=1;
      else if(std::abs(center[2]-s.length)<1e-12) attribute=2;
      else if(std::abs(center[1]+.5*s.thickness)<1e-12) attribute=3;
      else if(std::abs(center[1]-.5*s.thickness)<1e-12) attribute=4;
      serial.GetBdrElement(b)->SetAttribute(attribute);
   }
   serial.SetAttributes();
   ParMesh mesh(MPI_COMM_WORLD,serial);
   L2_FECollection l2(0,3); ND_FECollection nd(1,3); RT_FECollection rt(0,3); H1_FECollection h1(1,3);
   ParFiniteElementSpace l2s(&mesh,&l2), nds(&mesh,&nd), rts(&mesh,&rt), h1s(&mesh,&h1);
   Array<int> offsets(7); offsets[0]=0;
   offsets[1]=l2s.GetVSize();offsets[2]=offsets[1]+rts.GetVSize();
   offsets[3]=offsets[2]+h1s.GetVSize();offsets[4]=offsets[3]+nds.GetVSize();
   offsets[5]=offsets[4]+rts.GetVSize();offsets[6]=offsets[5]+l2s.GetVSize();
   BlockVector state(offsets);
   ParGridFunction T,F,P,E,B,W;
   T.MakeRef(&l2s,state,offsets[0]);F.MakeRef(&rts,state,offsets[1]);
   P.MakeRef(&h1s,state,offsets[2]);E.MakeRef(&nds,state,offsets[3]);
   B.MakeRef(&rts,state,offsets[4]);W.MakeRef(&l2s,state,offsets[5]);
   Array<int> em_bc(5),thermal_bc(5),potential_bc(5);
   em_bc=1;thermal_bc=1;potential_bc=0;potential_bc[0]=potential_bc[1]=1;
   AuditedJoule oper(offsets[6],l2s,nds,rts,h1s,em_bc,thermal_bc,potential_bc,s.mu,
                     {{1,s.sigma}},{{1,s.capacity}},{{1,1/s.capacity}},{{1,1/s.conductivity}});
   oper.Init(state);
   BackwardEulerSolver ode;ode.Init(oper);
   std::ofstream top_file,bottom_file,diagnostics;
   if(Mpi::Root())
   {
      top_file.open(std::string(output)+"/top.csv");bottom_file.open(std::string(output)+"/bottom.csv");
      diagnostics.open(std::string(output)+"/diagnostics.csv");
      MFEM_VERIFY(top_file && bottom_file && diagnostics,"Cannot create output files");
      top_file<<std::setprecision(17)<<"time_s,z_index,x_index,temperature_value\n";
      bottom_file<<std::setprecision(17)<<"time_s,z_index,x_index,temperature_value\n";
      diagnostics<<std::setprecision(17)<<"time_s,thermal_energy_j,projected_joule_power_w,direct_electric_power_w,magnetic_energy_j,current_lower_a,current_upper_a,maximum_linear_iterations,maximum_linear_relative_residual,minimum_delta_temperature_k,maximum_delta_temperature_k,cumulative_projected_joule_energy_j\n";
   }
   double t=0, cumulative_heat=0;
   auto snapshot=[&](int step,double dt)
   {
      std::vector<double> top,bottom;FaceMeans(mesh,T,s,top,bottom);
      double local[7]={0,0,0,0,0,std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
      Vector electric(3),magnetic(3);
      for(int e=0;e<mesh.GetNE();++e)
      {
         auto *tr=mesh.GetElementTransformation(e);
         const auto &rule=IntRules.Get(mesh.GetElementBaseGeometry(e),4);
         for(int q=0;q<rule.GetNPoints();++q)
         {
            const auto &ip=rule.IntPoint(q);tr->SetIntPoint(&ip);const double vol=ip.weight*tr->Weight();
            const double temp=T.GetValue(e,ip),heat=W.GetValue(e,ip);
            E.GetVectorValue(e,ip,electric);B.GetVectorValue(e,ip,magnetic);
            MFEM_VERIFY(std::isfinite(temp)&&std::isfinite(heat),"Nonfinite scalar fields");
            local[0]+=s.capacity*temp*vol;local[1]+=heat*vol;
            local[2]+=magnetic*magnetic/(2*s.mu)*vol;
            local[5]=std::min(local[5],temp);local[6]=std::max(local[6],temp);
         }
      }
      for(int b=0;b<mesh.GetNBE();++b)
      {
         const int a=mesh.GetBdrAttribute(b);if(a!=1&&a!=2) continue;
         auto *face=mesh.GetBdrFaceTransformations(b);const auto &rule=IntRules.Get(mesh.GetBdrElementBaseGeometry(b),4);
         for(int q=0;q<rule.GetNPoints();++q)
         {
            const auto &ip=rule.IntPoint(q);face->SetAllIntPoints(&ip);
            E.GetVectorValue(face->Elem1No,face->GetElement1IntPoint(),electric);
            local[a==1?3:4]+=s.sigma*electric[2]*ip.weight*face->Weight();
         }
      }
      double sums[5],minimum,maximum;
      MPI_Allreduce(local,sums,5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(local+5,&minimum,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
      MPI_Allreduce(local+6,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      int iterations;double residual;oper.SolverMetrics(iterations,residual);
      const double direct=oper.ElectricLosses(E);
      if(step>0)cumulative_heat+=dt*sums[1];
      if(Mpi::Root())
      {
         for(int iz=0;iz<s.nz;++iz)for(int ix=0;ix<s.nx;++ix)
         {
            top_file<<t<<","<<iz<<","<<ix<<","<<top[iz*s.nx+ix]<<"\n";
            bottom_file<<t<<","<<iz<<","<<ix<<","<<bottom[iz*s.nx+ix]<<"\n";
         }
         diagnostics<<t<<","<<sums[0]<<","<<sums[1]<<","<<direct<<","<<sums[2]<<","<<sums[3]<<","<<sums[4]<<","<<iterations<<","<<residual<<","<<minimum<<","<<maximum<<","<<cumulative_heat<<"\n";
         std::cout<<"step "<<step<<" time_s "<<t<<" delta_T_max_K "<<maximum<<" direct_power_W "<<direct<<std::endl;
      }
   };
   snapshot(0,0);
   const double dt0=s.final_time/s.steps;
   for(int step=1;step<=s.steps;++step)
   {
      double dt=dt0;ode.Step(state,t,dt);snapshot(step,dt0);
   }
   return 0;
}
