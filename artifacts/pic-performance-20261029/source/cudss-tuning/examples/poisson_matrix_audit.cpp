#include <mfem.hpp>
#include "cudss_matrix_audit.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <utility>
int main(int argc,char** argv) {
 try {
  if(argc!=2 || argv[1][0]!='/')throw std::runtime_error("Expected fresh absolute output prefix");
  auto mesh=mfem::Mesh::MakeCartesian2D(360,360,mfem::Element::QUADRILATERAL,true,.009,.009,false);
  mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
  mfem::ConstantCoefficient epsilon(8.8541878188e-12);mfem::BilinearForm form(&space);
  form.AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
  mfem::GridFunction xfield(&space),yfield(&space),trial(&space);trial=0;
  mfem::FunctionCoefficient xcoord([](const mfem::Vector& x){return x[0];});mfem::FunctionCoefficient ycoord([](const mfem::Vector& x){return x[1];});
  xfield.ProjectCoefficient(xcoord);yfield.ProjectCoefficient(ycoord);mfem::Array<int> essential;const double tol=1e-12*.009;
  for(int i=0;i<trial.Size();++i)for(const auto& f:std::array<std::array<double,4>,2>{{{.002,.007,.0018,.002},{.002,.007,.007,.0072}}})
   if(xfield[i]>=f[0]-tol && xfield[i]<=f[1]+tol && yfield[i]>=f[2]-tol && yfield[i]<=f[3]+tol)essential.Append(i);
  mfem::LinearForm rhs(&space);rhs=0;mfem::OperatorPtr op;mfem::Vector x,b;form.FormLinearSystem(essential,trial,rhs,op,x,b);
  auto* matrix=dynamic_cast<mfem::SparseMatrix*>(op.Ptr());if(!matrix)throw std::runtime_error("Expected sparse matrix");
  const int n=matrix->Height();std::vector<int> rows(n+1),cols;std::vector<double> values;
  for(int i=0;i<n;++i) {
   std::vector<std::pair<int,double>> row;
   for(int k=matrix->GetI()[i];k<matrix->GetI()[i+1];++k) {
    const int j=matrix->GetJ()[k];const double v=matrix->GetData()[k];
    if(std::abs(v-(*matrix)(j,i))>1e-14*std::max({std::abs(v),std::abs((*matrix)(j,i)),1e-100}))throw std::runtime_error("Asymmetric matrix");
    if(j>=i && v!=0)row.emplace_back(j,v);
   }
   std::sort(row.begin(),row.end());for(auto [j,v]:row){cols.push_back(j);values.push_back(v);}rows[i+1]=int(cols.size());
  }
  pic_cudss_tuning::Structure structure(rows,cols);const std::string prefix=argv[1];pic_cudss_tuning::matrix_arrays(prefix,rows,cols,values);
  std::ofstream out(prefix+".json");out<<"{\"matrix_rows\":"<<n<<",\"essential_dofs\":"<<essential.Size()<<",\"stored_upper_nnz\":"<<cols.size()<<",\"symmetric_nonzero_nnz\":"<<2*cols.size()-n<<',';
  structure.json(out);out<<",\"gpu_executed\":false,\"source\":\"CPU reconstruction of sealed Q1 matrix and foil elimination; compare arrays with native audit\"}\n";
  out.close();if(!out)throw std::runtime_error("Cannot write matrix structure record");
  std::cout<<"CPU matrix reconstruction: "<<n<<" rows, "<<cols.size()<<" stored upper entries\n";return 0;
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
