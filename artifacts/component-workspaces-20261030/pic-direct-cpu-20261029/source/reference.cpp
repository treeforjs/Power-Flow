#include "eigen_cpu_reference.hpp"
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <sys/resource.h>
#include <sched.h>

namespace {
constexpr double epsilon0=8.8541878188e-12;
struct System {
    std::unique_ptr<mfem::SparseMatrix> matrix;
    mfem::Vector vacuum_rhs;
    mfem::Array<int> essential;
};
// Same MFEM Q1 assembly and internal foil-node ownership as the frozen PIC.
// This reconstructs its fixed operator without loading particles or CUDA.
System assemble(int nx,int ny,bool foils) {
    const double width=foils?.009:.001,gap=width;
    auto mesh=mfem::Mesh::MakeCartesian2D(nx,ny,mfem::Element::QUADRILATERAL,true,width,gap,false);
    mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
    mfem::ConstantCoefficient epsilon(epsilon0);mfem::BilinearForm form(&space);
    form.AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
    mfem::LinearForm rhs(&space);rhs=0;mfem::GridFunction trial(&space);trial=0;
    System result;
    mfem::GridFunction xfield(&space),yfield(&space);
    mfem::FunctionCoefficient xcoord([](const mfem::Vector& x){return x[0];});
    mfem::FunctionCoefficient ycoord([](const mfem::Vector& x){return x[1];});
    xfield.ProjectCoefficient(xcoord);yfield.ProjectCoefficient(ycoord);
    const double tolerance=1e-12*std::max(width,gap);
    const std::array<std::array<double,4>,2> foil={{{.002,.007,.0018,.002},{.002,.007,.007,.0072}}};
    for(int i=0;i<trial.Size();++i) {
        if(foils) {
            int owner=0;
            for(const auto& f:foil) {
                if(xfield[i]>=f[0]-tolerance && xfield[i]<=f[1]+tolerance &&
                   yfield[i]>=f[2]-tolerance && yfield[i]<=f[3]+tolerance) {
                    result.essential.Append(i);trial[i]=owner==0?0:-1;
                }
                ++owner;
            }
        } else if(yfield[i]<tolerance || yfield[i]>gap-tolerance) {
            result.essential.Append(i);trial[i]=yfield[i]<.5*gap?0:-1;
        }
    }
    mfem::OperatorPtr matrix;mfem::Vector x,b;
    form.FormLinearSystem(result.essential,trial,rhs,matrix,x,b);
    const auto* sparse=dynamic_cast<mfem::SparseMatrix*>(matrix.Ptr());
    if(!sparse)throw std::runtime_error("Expected sparse Poisson matrix");
    result.matrix=std::make_unique<mfem::SparseMatrix>(*sparse);result.vacuum_rhs=b;
    return result;
}
std::uint64_t fingerprint(const mfem::SparseMatrix& matrix) {
    std::uint64_t hash=14695981039346656037ULL;
    const auto add=[&](const void* data,std::size_t size) {
        const auto* bytes=static_cast<const unsigned char*>(data);
        for(std::size_t i=0;i<size;++i){hash^=bytes[i];hash*=1099511628211ULL;}
    };
    add(matrix.GetI(),(matrix.Height()+1)*sizeof(int));
    add(matrix.GetJ(),matrix.NumNonZeroElems()*sizeof(int));
    add(matrix.GetData(),matrix.NumNonZeroElems()*sizeof(double));return hash;
}
void tiny_tests() {
    mfem::SparseMatrix matrix(3);
    const double data[3][3]={{4,1,0},{1,3,1},{0,1,2}};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(data[i][j])matrix.Add(i,j,data[i][j]);
    matrix.Finalize();EigenCpuReference factor(matrix);
    mfem::Vector expected(3),rhs(3),answer;expected[0]=1;expected[1]=-2;expected[2]=3;matrix.Mult(expected,rhs);
    for(int corrections:{0,2}) {
        factor.solve(rhs,answer,corrections);answer-=expected;
        if(answer.Normlinf()>2e-14)throw std::runtime_error("Tiny SPD solution mismatch");
    }
    bool rejected=false;
    mfem::SparseMatrix nonsymmetric(2);nonsymmetric.Add(0,0,1);nonsymmetric.Add(1,1,1);nonsymmetric.Add(0,1,1);nonsymmetric.Finalize();
    try {EigenCpuReference invalid(nonsymmetric);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Nonsymmetric matrix accepted");
    rejected=false;
    mfem::SparseMatrix indefinite(2);indefinite.Add(0,0,1);indefinite.Add(1,1,1);indefinite.Add(0,1,2);indefinite.Add(1,0,2);indefinite.Finalize();
    try {EigenCpuReference invalid(indefinite);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Indefinite matrix accepted");
}
void benchmark(std::ostream& out,System& system,const std::string& label,int repetitions) {
    const auto& matrix=*system.matrix;EigenCpuReference factor(matrix);
    mfem::Vector known(matrix.Height());
    for(int i=0;i<known.Size();++i)known[i]=std::sin(.017*i)+.2*std::cos(.031*i);
    for(int i=0;i<system.essential.Size();++i)known[system.essential[i]]=0;
    mfem::Vector manufactured(matrix.Height());matrix.Mult(known,manufactured);
    out<<"{\"label\":\""<<label<<"\",\"rows\":"<<matrix.Height()<<",\"nonzeros\":"<<matrix.NumNonZeroElems()
        <<",\"essential_dofs\":"<<system.essential.Size()<<",\"csr_fnv1a64\":\""<<std::hex<<fingerprint(matrix)<<std::dec
        <<"\",\"factor_nonzeros\":"<<factor.factor_nonzeros<<",\"conversion_seconds\":"<<factor.conversion_seconds
        <<",\"analysis_seconds\":"<<factor.analysis_seconds<<",\"factor_seconds\":"<<factor.factor_seconds<<",\"solves\":[\n";
    bool first=true;
    for(int repetition=0;repetition<repetitions;++repetition)for(int corrections:{0,2})for(int mode=0;mode<2;++mode) {
        mfem::Vector answer;const auto& rhs=mode==0?manufactured:system.vacuum_rhs;
        const auto result=factor.solve(rhs,answer,corrections);
        double solution_error=0;
        if(mode==0) {answer-=known;solution_error=answer.Normlinf();}
        if(!first)out<<",\n";first=false;
        out<<"{\"rhs\":\""<<(mode==0?"manufactured_interior":"vacuum_bias")<<"\",\"correction_solves\":"<<corrections
            <<",\"total_seconds\":"<<result.total_seconds<<",\"solve_seconds\":"<<result.solve_seconds
            <<",\"residual_seconds\":"<<result.residual_seconds<<",\"residual_norm\":"<<result.residual_norm
            <<",\"rhs_norm\":"<<result.rhs_norm<<",\"relative_residual\":"<<result.relative_residual
            <<",\"solution_linf_error\":"<<solution_error<<'}';
    }
    out<<"\n]}";
}
}
int main(int argc,char** argv) {
    try {
        if(argc<2 || argc>3 || (argc==3 && std::string(argv[2])!="--full"))
            throw std::invalid_argument("Usage: direct_cpu_reference OUTPUT_JSON [--full]");
        Eigen::setNbThreads(1);tiny_tests();
        std::cout<<"Tiny SPD, nonsymmetric, and indefinite reference tests passed\n"<<std::flush;
        std::ofstream out(argv[1]);out<<std::setprecision(17)<<"{\"reference_only\":true,\"eigen\":\""<<EIGEN_VERSION_STRING
            <<"\",\"residual_guard\":\"norm(Ax-b)<=max(1e-30,2e-13*norm(b))\",\"matrices\":[\n";
        auto tiny=assemble(6,9,false);benchmark(out,tiny,"small_Q1_Poisson",3);
        if(argc==3) {
            std::cout<<"Assembling separate 360x360 foil Q1 operator\n"<<std::flush;
            auto full=assemble(360,360,true);out<<",\n";benchmark(out,full,"reconstructed_360x360_foil_Q1",12);
            std::cout<<"Full operator factorization and cached solves passed the residual guard\n"<<std::flush;
        }
        struct rusage usage{};getrusage(RUSAGE_SELF,&usage);cpu_set_t cpus;CPU_ZERO(&cpus);sched_getaffinity(0,sizeof(cpus),&cpus);
        out<<"\n],\"max_rss_KiB\":"<<usage.ru_maxrss<<",\"nice\":"<<getpriority(PRIO_PROCESS,0)<<",\"logical_cpus\":[";bool first=true;
        for(int i=0;i<CPU_SETSIZE;++i)if(CPU_ISSET(i,&cpus)){if(!first)out<<',';first=false;out<<i;}
        out<<"]}\n";out.close();if(!out)throw std::runtime_error("Reference output failed");
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
