#define main heating_history_program
#include HISTORY_SOURCE
#undef main
#include <chrono>
#include <sstream>

void fixture(const std::filesystem::path& directory) {
    std::filesystem::create_directories(directory);
    const std::vector<double> xs={-.0045,-.0017,.0045},ys={-.0045,-.0021,.0003,.0045};
    auto mesh=mfem::Mesh::MakeCartesian2D(2,3,mfem::Element::QUADRILATERAL,true,.009,.009,false);
    for(int v=0;v<mesh.GetNV();++v) {
        auto* p=mesh.GetVertex(v);p[0]=xs[std::lround(p[0]/.0045)];p[1]=ys[std::lround(p[1]/.003)];
    }
    {std::ofstream out(directory/"em_mesh.mesh");out<<std::setprecision(17);mesh.Print(out);}
    {std::ofstream out(directory/"em_axes.csv");out<<std::setprecision(17)<<"axis,index,coordinate\n";
        for(int k=0;k<int(xs.size());++k)out<<"x,"<<k<<','<<xs[k]<<'\n';
        for(int k=0;k<int(ys.size());++k)out<<"y,"<<k<<','<<ys[k]<<'\n';}
    mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
    std::ofstream index(directory/"em_field_snapshots.csv");index<<std::setprecision(17)<<"time,az,rate\n";
    for(int snapshot=0;snapshot<3;++snapshot) {
        const std::string a="az_"+std::to_string(snapshot)+".gf",r="rate_"+std::to_string(snapshot)+".gf";
        index<<(snapshot+1)*1e-9<<','<<a<<','<<r<<'\n';
        mfem::GridFunction az(&space),rate(&space);
        for(int i=0;i<space.GetVSize();++i) {az[i]=std::sin(i*.71+snapshot*.3);rate[i]=std::cos(i*.53+snapshot*.7);}
        std::ofstream out_a(directory/a),out_r(directory/r);out_a<<std::setprecision(17);out_r<<std::setprecision(17);az.Save(out_a);rate.Save(out_r);
    }
    {std::ofstream out(directory/"voltage.csv");out<<std::setprecision(17)<<"start,end,voltage\n";
        out<<1e-9<<','<<1.5e-9<<",-1000\n"<<1.5e-9<<','<<2e-9<<",-2000\n"<<2e-9<<','<<3e-9<<",-3000\n";}
}
int main(int argc,char** argv) {
    if(argc!=4)return 2;
    const std::filesystem::path directory(argv[1]);fixture(directory);
    Fields fields(directory);Voltage voltage(directory/"voltage.csv");
    std::cout<<std::hexfloat;
    for(double t:{1e-9,1e-9,1.2e-9,1.5e-9,1.5e-9,2e-9,2e-9,2.9e-9,3e-9,2.3e-9,1e-9}) {
        for(double x:{.0001,.0028,.0045,.0089})for(double y:{.0001,.0024,.0048,.0089}) {
            const auto E=fields.value(x,y,t,false),B=fields.value(x,y,t,true);
            std::cout<<t<<' '<<x<<' '<<y<<' '<<voltage.at(t)<<' '<<fields.current<<' '<<fields.fraction;
            for(double v:E)std::cout<<' '<<v;for(double v:B)std::cout<<' '<<v;std::cout<<'\n';
        }
    }
    constexpr int samples=2000000;
    const auto start_clock=std::chrono::steady_clock::now();double checksum=0;
    for(int n=0;n<samples;++n) {
        double x=(.001+.998*double(n%997)/997)*.009,y=(.001+.998*double(n%991)/991)*.009;
        const auto E=fields.value(x,y,1.7e-9,false),B=fields.value(x,y,1.7e-9,true);
        checksum+=E[2]+B[0]+B[1]+voltage.at(1.7e-9);
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start_clock).count();
    std::ofstream result(argv[2]);result<<std::setprecision(17)<<"{\"case\":\"native_external_gather_same_midpoint\",\"samples\":"<<samples<<",\"seconds\":"<<seconds<<",\"checksum\":"<<checksum<<",\"label\":\""<<argv[3]<<"\"}\n";
}
