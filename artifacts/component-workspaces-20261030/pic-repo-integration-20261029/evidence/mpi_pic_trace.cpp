#include "raikou/mfem_gap_pic.h"
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sched.h>
#include <stdexcept>
#include <sys/resource.h>
using namespace raikou::mfem_backend;
void trace(GapPic& pic) {
    const auto diagnostic=pic.diagnostics();
    static_assert(sizeof(diagnostic)==25*sizeof(double));
    for(std::size_t i=0;i<25;++i) {
        double value;
        std::memcpy(&value,reinterpret_cast<const char*>(&diagnostic)+i*sizeof(value),sizeof(value));
        std::cout<<value<<' ';
    }
    std::cout<<'\n';
    for(int i=0;i<pic.potential().Size();++i)std::cout<<pic.potential()[i]<<' ';
    std::cout<<'\n';
    for(const auto& p:pic.particles()) {
        std::cout<<p.id<<' '<<p.weight;
        for(double value:p.position_m)std::cout<<' '<<value;
        for(double value:p.proper_velocity_m_s)std::cout<<' '<<value;
        std::cout<<'\n';
    }
    for(double value:pic.horizontal_transfers_C())std::cout<<value<<' ';
    std::cout<<'\n';
    for(double value:pic.vertical_transfers_C())std::cout<<value<<' ';
    std::cout<<'\n';
}
int main(int argc,char** argv) {
    mfem::Mpi::Init(argc,argv);mfem::Hypre::Init();
    if(argc!=2)return 2;
    cpu_set_t cpus;CPU_ZERO(&cpus);sched_getaffinity(0,sizeof(cpus),&cpus);
    int local=(CPU_COUNT(&cpus)==1 && CPU_ISSET(6,&cpus) && getpriority(PRIO_PROCESS,0)==15),all=0;
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(!all)throw std::runtime_error("MPI ranks must inherit CPU 6 and nice 15");
    if(mfem::Mpi::WorldRank()!=0) {
        pic_parallel_workers();pic_parallel_record("");return 0;
    }
    std::filesystem::create_directories(argv[1]);std::cout<<std::hexfloat;
    {
        GapPicOptions options;options.nx=6;options.ny=8;options.distributed_fields=true;
        options.upper_potential_V=3.;options.magnetic_field_T={.1,.2,.3};
        GapPic pic(options);pic.solve_fields();trace(pic);
        for(int i=0;i<6;++i) {
            GapParticle p;p.position_m={(.1+.13*i)*options.width_m,(.2+.09*i)*options.gap_m,0};
            p.proper_velocity_m_s={i%2?1e6:-1e6,2e6,3e5};p.weight=1e-6;pic.add_particle(p);
        }
        for(int step=0;step<4;++step) {
            pic.set_potentials(.3*step,3.-.7*step);pic.advance(2e-13);trace(pic);
        }
        pic.release_field_solver();pic.solve_fields();trace(pic);
    }
    pic_parallel_finish();pic_parallel_record(argv[1]);
}
