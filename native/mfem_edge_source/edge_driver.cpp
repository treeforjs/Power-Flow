// MFEM A_z magnetoquasistatic conductors/true-vacuum current constraint.
// This is a provisional, constant-property model, not a stock Joule miniapp run.
#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;

struct Settings
{
   int nx=64, ny=8, outside=16, gap_cells=48;
   double width=.005, thickness=.0002, gap=.005, extent=.010;
   double length=.010, sigma=3e6, cv=2.8e6, k=16, mu=4e-7*M_PI;
   const char *schedule="schedule.csv", *output=".", *snapshot_times="";
   std::vector<double> field_times={40e-9,60e-9};
};

std::vector<double> AxisSegment(double a,double b,int n)
{
   std::vector<double> x(n+1); for(int i=0;i<=n;++i)x[i]=a+(b-a)*i/n;
   x[0]=a;x[n]=b;return x;
}
void Append(std::vector<double> &x,const std::vector<double> &y)
{ x.insert(x.end(),y.begin()+(x.empty()?0:1),y.end()); }
int Bin(const std::vector<double> &x,double value)
{
   int i=int(std::upper_bound(x.begin(),x.end(),value)-x.begin())-1;
   return std::clamp(i,0,int(x.size())-2);
}

struct LinearSystem
{
   SparseMatrix &matrix;
   HYPRE_BigInt partition[2];
   std::unique_ptr<HypreParMatrix> parallel;
   std::unique_ptr<HypreBoomerAMG> amg;
   std::unique_ptr<CGSolver> cg;
   int iterations=0;double residual=0;
   explicit LinearSystem(SparseMatrix &a):matrix(a)
   {
      partition[0]=0;partition[1]=a.Height();
      parallel=std::make_unique<HypreParMatrix>(MPI_COMM_SELF,HYPRE_BigInt(a.Height()),partition,&a);
      amg=std::make_unique<HypreBoomerAMG>(*parallel);amg->SetPrintLevel(0);
      cg=std::make_unique<CGSolver>(MPI_COMM_SELF);cg->SetRelTol(2e-13);cg->SetAbsTol(1e-20);
      cg->SetMaxIter(3000);cg->SetPrintLevel(0);cg->SetPreconditioner(*amg);cg->SetOperator(*parallel);
   }
   void Solve(const Vector &rhs,Vector &answer)
   {
      answer.SetSize(rhs.Size());answer=0.;cg->Mult(rhs,answer);
      MFEM_VERIFY(cg->GetConverged(),"MFEM linear system did not converge");
      Vector error(rhs.Size());matrix.Mult(answer,error);error-=rhs;
      const double norm=rhs.Norml2(),relative=norm>0?error.Norml2()/norm:error.Norml2();
      iterations=std::max(iterations,cg->GetNumIterations());residual=std::max(residual,relative);
      MFEM_VERIFY(std::isfinite(residual),"Nonfinite solver residual");
   }
};

struct EMGrid
{
   Settings &s;
   std::vector<double> x,y;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<H1_FECollection> collection;
   std::unique_ptr<FiniteElementSpace> space;
   std::unique_ptr<GridFunction> potential,rate;
   std::unique_ptr<BilinearForm> mass,stiffness;
   Vector f[2];double d[2]={0,0};Array<int> essential;std::vector<int> lookup;
   explicit EMGrid(Settings &settings):s(settings)
   {
      const double lo=s.gap/2,hi=lo+s.thickness;
      Append(x,AxisSegment(-s.extent,-s.width/2,s.outside));Append(x,AxisSegment(-s.width/2,s.width/2,s.nx));
      Append(x,AxisSegment(s.width/2,s.extent,s.outside));
      Append(y,AxisSegment(-s.extent,-hi,s.outside));Append(y,AxisSegment(-hi,-lo,s.ny));
      Append(y,AxisSegment(-lo,lo,s.gap_cells));Append(y,AxisSegment(lo,hi,s.ny));Append(y,AxisSegment(hi,s.extent,s.outside));
      const int nx=x.size()-1,ny=y.size()-1;
      mesh=std::make_unique<Mesh>(Mesh::MakeCartesian2D(nx,ny,Element::QUADRILATERAL,true,nx,ny));
      for(int v=0;v<mesh->GetNV();++v){auto *p=mesh->GetVertex(v);p[0]=x[std::lround(p[0])];p[1]=y[std::lround(p[1])];}
      lookup.resize(nx*ny);
      for(int e=0;e<mesh->GetNE();++e)
      {
         auto *tr=mesh->GetElementTransformation(e);Vector p(2);tr->Transform(Geometries.GetCenter(Geometry::SQUARE),p);
         const int ix=Bin(x,p[0]),iy=Bin(y,p[1]);lookup[iy*nx+ix]=e;
         int attribute=1;
         if(std::abs(p[0])<s.width/2 && p[1]>lo && p[1]<hi)attribute=2;
         if(std::abs(p[0])<s.width/2 && p[1]>-hi && p[1]<-lo)attribute=3;
         mesh->GetElement(e)->SetAttribute(attribute);
      }
      mesh->SetAttributes();collection=std::make_unique<H1_FECollection>(1,2);
      space=std::make_unique<FiniteElementSpace>(mesh.get(),collection.get());
      potential=std::make_unique<GridFunction>(space.get());rate=std::make_unique<GridFunction>(space.get());*potential=0.;*rate=0.;
      Vector conductivities(3);conductivities=0.;conductivities[1]=conductivities[2]=s.sigma;
      PWConstCoefficient sigma(conductivities);ConstantCoefficient nu(1/s.mu);
      mass=std::make_unique<BilinearForm>(space.get());mass->AddDomainIntegrator(new MassIntegrator(sigma));mass->Assemble();mass->Finalize();
      stiffness=std::make_unique<BilinearForm>(space.get());stiffness->AddDomainIntegrator(new DiffusionIntegrator(nu));stiffness->Assemble();stiffness->Finalize();
      Array<int> boundary(mesh->bdr_attributes.Max());boundary=1;space->GetEssentialTrueDofs(boundary,essential);
      for(int conductor=0;conductor<2;++conductor)
      {
         Vector values(3);values=0.;values[conductor+1]=s.sigma;PWConstCoefficient coefficient(values);
         LinearForm form(space.get());form.AddDomainIntegrator(new DomainLFIntegrator(coefficient));form.Assemble();f[conductor]=form;
         d[conductor]=s.sigma*s.width*s.thickness;
      }
   }
   std::unique_ptr<SparseMatrix> Matrix(double inverse_dt)
   {
      auto result=std::unique_ptr<SparseMatrix>(Add(1.,stiffness->SpMat(),inverse_dt,mass->SpMat()));
      for(int i=0;i<essential.Size();++i)result->EliminateRowCol(essential[i]);return result;
   }
   void ZeroBoundary(Vector &v){for(int i=0;i<essential.Size();++i)v[essential[i]]=0.;}
   void InitialDC(double current)
   {
      auto matrix=Matrix(0);LinearSystem system(*matrix);Vector rhs(f[0]);rhs*=current/d[0];rhs.Add(-current/d[1],f[1]);ZeroBoundary(rhs);
      Vector a;system.Solve(rhs,a);*potential=a;*rate=0.;
   }
   double RateAt(double px,double py)
   {
      const int ix=Bin(x,px),iy=Bin(y,py),e=lookup[iy*(x.size()-1)+ix];
      IntegrationPoint ip;ip.Set2((px-x[ix])/(x[ix+1]-x[ix]),(py-y[iy])/(y[iy+1]-y[iy]));
      return rate->GetValue(e,ip);
   }
   double MagneticEnergy()
   {Vector value(potential->Size());stiffness->Mult(*potential,value);return .5*((*potential)*value)*s.length;}
};

class Heating : public Coefficient
{
   EMGrid &em;double &gradient;
public:
   Heating(EMGrid &g,double &u):em(g),gradient(u){}
   real_t Eval(ElementTransformation &tr,const IntegrationPoint &ip) override
   {Vector p(2);tr.Transform(ip,p);const double j=em.s.sigma*(gradient-em.RateAt(p[0],p[1]));return j*j/em.s.sigma;}
};

struct Thermal
{
   Settings &s;EMGrid &em;int conductor;double lower,upper,gradient=0;
   Mesh mesh;H1_FECollection collection;FiniteElementSpace space;GridFunction temperature;
   BilinearForm mass,stiffness;Heating heating;
   Thermal(Settings &settings,EMGrid &fields,int which):s(settings),em(fields),conductor(which),
      lower(which==0?s.gap/2:-s.gap/2-s.thickness),upper(lower+s.thickness),
      mesh(Mesh::MakeCartesian2D(s.nx,s.ny,Element::QUADRILATERAL,true,s.width,s.thickness)),
      collection(1,2),space(&mesh,&collection),temperature(&space),mass(&space),stiffness(&space),heating(em,gradient)
   {
      for(int i=0;i<mesh.GetNV();++i){auto *p=mesh.GetVertex(i);p[0]-=s.width/2;p[1]+=lower;}
      temperature=0.;ConstantCoefficient cv(s.cv),k(s.k);
      mass.AddDomainIntegrator(new MassIntegrator(cv));mass.Assemble();mass.Finalize();
      stiffness.AddDomainIntegrator(new DiffusionIntegrator(k));stiffness.Assemble();stiffness.Finalize();
   }
   double Energy()
   {Vector value(temperature.Size()),one(temperature.Size());one=1.;mass.Mult(temperature,value);return (one*value)*s.length;}
   double Step(double dt)
   {
      LinearForm load(&space);load.AddDomainIntegrator(new DomainLFIntegrator(heating,6));load.Assemble();
      Vector rhs(temperature.Size());mass.Mult(temperature,rhs);rhs/=dt;rhs+=load;
      std::unique_ptr<SparseMatrix> matrix(Add(1/dt,mass.SpMat(),1.,stiffness.SpMat()));LinearSystem solver(*matrix);
      Vector result;solver.Solve(rhs,result);temperature=result;
      Vector one(load.Size());one=1.;return (one*load)*s.length;
   }
   void Faces(std::vector<double>&top,std::vector<double>&bottom)
   {
      top.assign(s.nx,0);bottom.assign(s.nx,0);std::vector<double> at(s.nx,0),ab(s.nx,0);
      for(int b=0;b<mesh.GetNBE();++b)
      {
         auto *face=mesh.GetBdrFaceTransformations(b);Vector center(2);face->Transform(Geometries.GetCenter(Geometry::SEGMENT),center);
         int side=std::abs(center[1]-upper)<1e-12?1:std::abs(center[1]-lower)<1e-12?-1:0;
         if(!side)continue;
         const int ix=std::clamp(int((center[0]+s.width/2)/s.width*s.nx),0,s.nx-1);
         const auto &rule=IntRules.Get(Geometry::SEGMENT,5);
         for(int q=0;q<rule.GetNPoints();++q)
         {
            const auto &ip=rule.IntPoint(q);face->SetAllIntPoints(&ip);double area=ip.weight*face->Weight();
            double value=temperature.GetValue(face->Elem1No,face->GetElement1IntPoint());
            (side>0?top[ix]:bottom[ix])+=value*area;(side>0?at[ix]:ab[ix])+=area;
         }
      }
      for(int ix=0;ix<s.nx;++ix){MFEM_VERIFY(std::abs(at[ix]-s.width/s.nx)<1e-12 && std::abs(ab[ix]-s.width/s.nx)<1e-12,"Incomplete trace width support");top[ix]/=at[ix];bottom[ix]/=ab[ix];}
   }
   void SideFaces(std::vector<double>&left,std::vector<double>&right,int&left_attribute,int&right_attribute)
   {
      left.assign(s.ny,0);right.assign(s.ny,0);std::vector<double> al(s.ny,0),ar(s.ny,0);
      left_attribute=right_attribute=0;
      for(int b=0;b<mesh.GetNBE();++b)
      {
         auto *face=mesh.GetBdrFaceTransformations(b);Vector center(2);face->Transform(Geometries.GetCenter(Geometry::SEGMENT),center);
         const int side=std::abs(center[0]+s.width/2)<1e-12?-1:std::abs(center[0]-s.width/2)<1e-12?1:0;
         if(!side)continue;
         const int iy=std::clamp(int((center[1]-lower)/s.thickness*s.ny),0,s.ny-1);
         int &attribute=side<0?left_attribute:right_attribute;
         const int selected_attribute=mesh.GetBdrElement(b)->GetAttribute();
         MFEM_VERIFY(attribute==0 || attribute==selected_attribute,"Side face has inconsistent boundary attributes");
         attribute=selected_attribute;
         const auto &rule=IntRules.Get(Geometry::SEGMENT,5);
         for(int q=0;q<rule.GetNPoints();++q)
         {
            const auto &ip=rule.IntPoint(q);face->SetAllIntPoints(&ip);const double area=ip.weight*face->Weight();
            const double value=temperature.GetValue(face->Elem1No,face->GetElement1IntPoint());
            (side<0?left[iy]:right[iy])+=value*area;(side<0?al[iy]:ar[iy])+=area;
         }
      }
      MFEM_VERIFY(left_attribute>0 && right_attribute>0,"Missing material side face");
      for(int iy=0;iy<s.ny;++iy)
      {
         MFEM_VERIFY(std::abs(al[iy]-s.thickness/s.ny)<1e-12 && std::abs(ar[iy]-s.thickness/s.ny)<1e-12,"Incomplete side trace thickness support");
         left[iy]/=al[iy];right[iy]/=ar[iy];
      }
   }
};

struct Metrics
{double current[2]={0,0},power[2]={0,0},edge_q=0,center_q=0,edge_T=0,center_T=0,edge_area=0,center_area=0;};

Metrics Inspect(Thermal &foil,Thermal &return_foil,double stamp,std::ofstream &cells)
{
   Metrics m;Settings &s=foil.s;
   const bool emit_field=std::any_of(s.field_times.begin(),s.field_times.end(),[stamp](double t){return std::abs(stamp-t)<1e-20;});
   for(Thermal *body:{&foil,&return_foil})for(int e=0;e<body->mesh.GetNE();++e)
   {
      auto *tr=body->mesh.GetElementTransformation(e);Vector p(2),center(2);tr->Transform(Geometries.GetCenter(Geometry::SQUARE),center);
      double area=0,temp=0,current=0,power=0;const auto &rule=IntRules.Get(Geometry::SQUARE,6);
      for(int q=0;q<rule.GetNPoints();++q)
      {
         const auto &ip=rule.IntPoint(q);tr->SetIntPoint(&ip);tr->Transform(ip,p);double weight=ip.weight*tr->Weight();
         double j=s.sigma*(body->gradient-body->em.RateAt(p[0],p[1]));double heat=j*j/s.sigma;
         double value=body->temperature.GetValue(e,ip);MFEM_VERIFY(std::isfinite(j)&&std::isfinite(heat)&&std::isfinite(value),"Nonfinite fields");
         area+=weight;temp+=value*weight;current+=j*weight;power+=heat*weight;
         if(body->conductor==0)
         {
            if(std::abs(p[0])>=.4*s.width){m.edge_area+=weight;m.edge_q+=heat*weight;m.edge_T+=value*weight;}
            if(std::abs(p[0])<=.1*s.width){m.center_area+=weight;m.center_q+=heat*weight;m.center_T+=value*weight;}
         }
      }
      m.current[body->conductor]+=current;m.power[body->conductor]+=power*s.length;
      if(body->conductor==0 && emit_field)
      {
         const int ix=std::clamp(int((center[0]+s.width/2)/s.width*s.nx),0,s.nx-1);
         const int iy=std::clamp(int((center[1]-body->lower)/s.thickness*s.ny),0,s.ny-1);
         cells<<stamp<<","<<-s.width/2+ix*s.width/s.nx<<","<<-s.width/2+(ix+1)*s.width/s.nx<<","<<iy*s.thickness/s.ny<<","<<(iy+1)*s.thickness/s.ny<<","<<300+temp/area<<","<<current/area<<","<<power/area<<"\n";
      }
   }
   m.edge_q/=m.edge_area;m.center_q/=m.center_area;m.edge_T/=m.edge_area;m.center_T/=m.center_area;return m;
}

int main(int argc,char **argv)
{
   Mpi::Init(argc,argv);Hypre::Init();MFEM_VERIFY(Mpi::WorldSize()==1,"This explicit 2D reference uses one MPI rank");
   Settings s;OptionsParser options(argc,argv);
   options.AddOption(&s.nx,"-nx","--width-cells","Cells across each foil width");
   options.AddOption(&s.ny,"-ny","--thickness-cells","Cells through each foil thickness");
   options.AddOption(&s.outside,"-no","--outside-cells","Cells in each exterior margin");
   options.AddOption(&s.gap_cells,"-ng","--gap-cells","Cells through gap");
   options.AddOption(&s.width,"-w","--width-m","Foil width");options.AddOption(&s.thickness,"-h","--thickness-m","Foil thickness");
   options.AddOption(&s.gap,"-g","--gap-m","Gap between facing surfaces");options.AddOption(&s.extent,"-ext","--extent-m","Exterior half extent");
   options.AddOption(&s.length,"-l","--length-m","Invariant current length");options.AddOption(&s.sigma,"-sig","--sigma-s-m","Metal conductivity");
   options.AddOption(&s.cv,"-cv","--capacity-j-m3-k","Volumetric capacity");options.AddOption(&s.k,"-k","--conductivity-w-m-k","Thermal conductivity");
   options.AddOption(&s.schedule,"-in","--schedule","Bounded current/time schedule CSV");options.AddOption(&s.output,"-out","--output-dir","Existing output directory");
   options.AddOption(&s.snapshot_times,"-snap","--snapshot-times","Optional CSV of actual native field output times; default 40/60 ns");
   options.Parse();if(!options.Good()){options.PrintUsage(std::cerr);return 2;}
   {
      std::ofstream dependency(std::string(s.output)+"/mfem_dependency.json");MFEM_VERIFY(dependency,"Cannot open executed MFEM dependency record");
      dependency<<"{\"version\":\""<<GetVersionMajor()<<"."<<GetVersionMinor()<<"."<<GetVersionPatch()
                <<"\",\"major\":"<<GetVersionMajor()<<",\"minor\":"<<GetVersionMinor()<<",\"patch\":"<<GetVersionPatch()
                <<",\"version_string\":"<<std::quoted(GetVersionStr())<<",\"git_string\":"<<std::quoted(GetGitStr())<<"}\n";
   }
   MFEM_VERIFY(GetVersionMajor()==4&&GetVersionMinor()==10&&GetVersionPatch()==0,"This source handoff pins the executed MFEM4.10 release");
   std::vector<double> times,currents;std::ifstream input(s.schedule);std::string line;std::getline(input,line);
   while(std::getline(input,line)){auto comma=line.find(',');MFEM_VERIFY(comma!=std::string::npos,"Invalid schedule");times.push_back(std::stod(line.substr(0,comma)));currents.push_back(std::stod(line.substr(comma+1)));}
   MFEM_VERIFY(times.size()>1 && std::isfinite(times[0]) && std::isfinite(currents[0]),"Missing finite drive schedule");
   for(size_t i=1;i<times.size();++i)MFEM_VERIFY(std::isfinite(times[i])&&std::isfinite(currents[i])&&times[i]>times[i-1],"Nonfinite or unordered drive schedule");
   if(std::string(s.snapshot_times).size())
   {
      std::ifstream snapshots(s.snapshot_times);MFEM_VERIFY(snapshots,"Cannot open field snapshot time CSV");
      s.field_times.clear();std::getline(snapshots,line);
      while(std::getline(snapshots,line))s.field_times.push_back(std::stod(line));
   }
   MFEM_VERIFY(!s.field_times.empty(),"Missing field snapshot times");
   for(size_t i=0;i<s.field_times.size();++i)
   {
      const double stamp=s.field_times[i];
      MFEM_VERIFY(std::isfinite(stamp)&&(i==0||stamp>s.field_times[i-1]),"Nonfinite or unordered field snapshot time");
      MFEM_VERIFY(std::any_of(times.begin(),times.end(),[stamp](double t){return std::abs(stamp-t)<1e-20;}),"Field output must use an actual native timestep");
   }
   EMGrid em(s);em.InitialDC(currents[0]);Thermal foil(s,em,0),return_foil(s,em,1);foil.gradient=currents[0]/em.d[0];return_foil.gradient=-currents[0]/em.d[1];
   std::ofstream top(std::string(s.output)+"/top.csv"),bottom(std::string(s.output)+"/bottom.csv"),left(std::string(s.output)+"/left.csv"),right(std::string(s.output)+"/right.csv"),diagnostics(std::string(s.output)+"/diagnostics.csv"),cells(std::string(s.output)+"/foil_snapshots.csv");
   MFEM_VERIFY(top&&bottom&&left&&right&&diagnostics&&cells,"Cannot open output files");
   for(auto *stream:{&top,&bottom,&left,&right,&diagnostics,&cells})*stream<<std::setprecision(17);
   for(auto *stream:{&top,&bottom,&left,&right})*stream<<"time_s,z_index,x_index,temperature_value\n";
   cells<<"time_s,x_left_m,x_right_m,depth_left_m,depth_right_m,temperature_k,current_density_a_m2,joule_power_w_m3\n";
   diagnostics<<"time_s,target_current_a,foil_current_a,return_current_a,foil_gradient_v_m,return_gradient_v_m,foil_joule_power_w,return_joule_power_w,total_joule_power_w,magnetic_energy_j,thermal_energy_j,backward_euler_magnetic_damping_j,drive_work_increment_j,em_balance_residual_j,thermal_balance_residual_j,edge_to_center_joule_ratio,edge_mean_delta_temperature_k,center_mean_delta_temperature_k,linear_iterations,linear_relative_residual,schur_condition_proxy\n";
   double old_magnetic=em.MagneticEnergy(),old_thermal=0;
   auto snapshot=[&](size_t step,double damping,double work,double em_residual,double thermal_residual,int iterations,double residual,double condition)
   {
      std::vector<double> t,b;foil.Faces(t,b);for(int ix=0;ix<s.nx;++ix){top<<times[step]<<",0,"<<ix<<","<<t[ix]<<"\n";bottom<<times[step]<<",0,"<<ix<<","<<b[ix]<<"\n";}
      std::vector<double> l,r;int left_attribute=0,right_attribute=0;foil.SideFaces(l,r,left_attribute,right_attribute);
      for(int iy=0;iy<s.ny;++iy){left<<times[step]<<",0,"<<iy<<","<<l[iy]<<"\n";right<<times[step]<<",0,"<<iy<<","<<r[iy]<<"\n";}
      if(step==0)
      {
         std::ofstream geometry(std::string(s.output)+"/face_geometry.json");MFEM_VERIFY(geometry,"Cannot open face geometry record");
         geometry<<"{\"schema\":\"mhdlab-mfem-side-face-geometry-v1\",\"left_boundary_attribute\":"<<left_attribute<<",\"right_boundary_attribute\":"<<right_attribute<<",\"physical_in_plane_axis\":\"y\",\"physical_normal_axis\":\"x\",\"sampling\":\"one_sided_material_boundary_trace_patch_average\"}\n";
      }
      const auto m=Inspect(foil,return_foil,times[step],cells);
      diagnostics<<times[step]<<","<<currents[step]<<","<<m.current[0]<<","<<m.current[1]<<","<<foil.gradient<<","<<return_foil.gradient<<","<<m.power[0]<<","<<m.power[1]<<","<<m.power[0]+m.power[1]<<","<<em.MagneticEnergy()<<","<<foil.Energy()+return_foil.Energy()<<","<<damping<<","<<work<<","<<em_residual<<","<<thermal_residual<<","<<m.edge_q/m.center_q<<","<<m.edge_T<<","<<m.center_T<<","<<iterations<<","<<residual<<","<<condition<<"\n";
      const double error=std::max(std::abs(m.current[0]-currents[step]),std::abs(m.current[1]+currents[step]));
      MFEM_VERIFY(error<=std::max(1e-3,1e-8*std::abs(currents[step])),"Current constraint reconstruction failed");
      if(step%10==0)std::cout<<"step "<<step<<" t_ns "<<times[step]*1e9<<" I_A "<<m.current[0]<<" edge_q_ratio "<<m.edge_q/m.center_q<<" edge_T_K "<<300+m.edge_T<<" center_T_K "<<300+m.center_T<<std::endl;
   };
   snapshot(0,0,0,0,0,0,0,0);
   for(size_t step=1;step<times.size();++step)
   {
      const double dt=times[step]-times[step-1];Vector old(*em.potential),rhs(old.Size()),base,q0,q1;
      auto matrix=em.Matrix(1/dt);LinearSystem solver(*matrix);em.mass->Mult(old,rhs);rhs/=dt;em.ZeroBoundary(rhs);solver.Solve(rhs,base);
      solver.Solve(em.f[0],q0);solver.Solve(em.f[1],q1);
      Vector delta(base);delta-=old;
      const double c00=em.d[0]-(em.f[0]*q0)/dt,c01=-(em.f[0]*q1)/dt,c10=-(em.f[1]*q0)/dt,c11=em.d[1]-(em.f[1]*q1)/dt;
      const double r0=currents[step]+(em.f[0]*delta)/dt,r1=-currents[step]+(em.f[1]*delta)/dt;
      const double determinant=c00*c11-c01*c10;MFEM_VERIFY(c00>0&&c11>0&&determinant>0,"Indefinite current Schur complement");
      foil.gradient=(r0*c11-c01*r1)/determinant;return_foil.gradient=(c00*r1-c10*r0)/determinant;
      *em.potential=base;em.potential->Add(foil.gradient,q0);em.potential->Add(return_foil.gradient,q1);
      *em.rate=*em.potential;*em.rate-=old;*em.rate/=dt;
      Vector change(*em.potential);change-=old;Vector kchange(change.Size());em.stiffness->Mult(change,kchange);
      const double damping=.5*(change*kchange)*s.length;
      const double p0=foil.Step(dt),p1=return_foil.Step(dt),joule=p0+p1;
      const double work=dt*s.length*(foil.gradient*currents[step]-return_foil.gradient*currents[step]);
      const double magnetic=em.MagneticEnergy(),thermal=foil.Energy()+return_foil.Energy();
      const double em_residual=work-(magnetic-old_magnetic)-dt*joule-damping;
      const double thermal_residual=(thermal-old_thermal)-dt*joule;
      snapshot(step,damping,work,em_residual,thermal_residual,solver.iterations,solver.residual,
               (std::abs(c00)+std::abs(c01)+std::abs(c10)+std::abs(c11))*(std::abs(c00)+std::abs(c11))/determinant);
      old_magnetic=magnetic;old_thermal=thermal;
   }
   return 0;
}
