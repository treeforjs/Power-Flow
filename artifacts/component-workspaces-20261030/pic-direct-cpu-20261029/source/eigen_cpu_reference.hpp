#pragma once

#include <mfem.hpp>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

// Experimental reference only. Eigen owns ordering, factorization and the
// triangular solves. The original MFEM operator owns residual evaluation.
class EigenCpuReference {
    using Clock=std::chrono::steady_clock;
    using Matrix=Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
    const mfem::SparseMatrix& matrix_;
    Eigen::SimplicialLLT<Matrix,Eigen::Upper,Eigen::AMDOrdering<int>> factor_;
    static double elapsed(Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now()-start).count();
    }
public:
    double conversion_seconds=0,analysis_seconds=0,factor_seconds=0;
    int factor_nonzeros=0;
    explicit EigenCpuReference(const mfem::SparseMatrix& matrix):matrix_(matrix) {
        if(matrix.Height()!=matrix.Width() || !matrix.Finalized())
            throw std::invalid_argument("Reference requires a finalized square matrix");
        auto start=Clock::now();
        std::vector<Eigen::Triplet<double,int>> entries;
        entries.reserve(matrix.NumNonZeroElems()/2+matrix.Height());
        for(int i=0;i<matrix.Height();++i) {
            bool diagonal=false;
            for(int k=matrix.GetI()[i];k<matrix.GetI()[i+1];++k) {
                const int j=matrix.GetJ()[k];const double value=matrix.GetData()[k];
                if(!std::isfinite(value) || std::abs(value-matrix(j,i))>
                   1e-14*std::max({std::abs(value),std::abs(matrix(j,i)),1e-100}))
                    throw std::runtime_error("Reference requires symmetric finite input");
                if(i==j) {diagonal=true;if(!(value>0))throw std::runtime_error("Nonpositive diagonal");}
                if(j>=i && value!=0)entries.emplace_back(i,j,value);
            }
            if(!diagonal)throw std::runtime_error("Missing diagonal");
        }
        Matrix upper(matrix.Height(),matrix.Width());
        upper.setFromTriplets(entries.begin(),entries.end());upper.makeCompressed();
        conversion_seconds=elapsed(start);
        start=Clock::now();factor_.analyzePattern(upper);analysis_seconds=elapsed(start);
        start=Clock::now();factor_.factorize(upper);factor_seconds=elapsed(start);
        if(factor_.info()!=Eigen::Success)throw std::runtime_error("Cholesky factorization failed");
        factor_nonzeros=static_cast<int>(factor_.matrixU().nestedExpression().nonZeros());
    }
    struct Result {
        double total_seconds=0,solve_seconds=0,residual_seconds=0;
        double residual_norm=0,rhs_norm=0,relative_residual=0;
        int correction_solves=0;
    };
    Result solve(const mfem::Vector& rhs,mfem::Vector& answer,int corrections) const {
        if(rhs.Size()!=matrix_.Height() || (corrections!=0 && corrections!=2))
            throw std::invalid_argument("Reference accepts matching RHS and zero or two correction solves");
        Result result;result.correction_solves=corrections;
        const auto total_start=Clock::now();auto start=Clock::now();
        Eigen::VectorXd solution=factor_.solve(Eigen::Map<const Eigen::VectorXd>(rhs.HostRead(),rhs.Size()));
        if(factor_.info()!=Eigen::Success)throw std::runtime_error("Cholesky solve failed");
        result.solve_seconds+=elapsed(start);
        answer.SetSize(rhs.Size());mfem::Vector error(rhs.Size());
        for(int iteration=0;iteration<=corrections;++iteration) {
            Eigen::Map<Eigen::VectorXd>(answer.HostWrite(),answer.Size())=solution;
            start=Clock::now();matrix_.Mult(answer,error);error-=rhs;
            if(iteration==corrections) {
                result.residual_norm=error.Norml2();result.rhs_norm=rhs.Norml2();
                result.relative_residual=result.rhs_norm>0?result.residual_norm/result.rhs_norm:result.residual_norm;
            }
            result.residual_seconds+=elapsed(start);
            if(iteration==corrections)break;
            start=Clock::now();
            solution-=factor_.solve(Eigen::Map<const Eigen::VectorXd>(error.HostRead(),error.Size()));
            if(factor_.info()!=Eigen::Success)throw std::runtime_error("Cholesky correction failed");
            result.solve_seconds+=elapsed(start);
        }
        if(!std::isfinite(result.relative_residual) ||
           result.residual_norm>std::max(1e-30,2e-13*result.rhs_norm))
            throw std::runtime_error("Unchanged Poisson residual guard failed");
        result.total_seconds=elapsed(total_start);
        return result;
    }
};
