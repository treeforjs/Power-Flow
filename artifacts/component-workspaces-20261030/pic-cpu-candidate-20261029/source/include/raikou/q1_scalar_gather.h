#pragma once

#include <mfem.hpp>
#include <array>
#include <vector>

namespace raikou::mfem_backend {

// Serial scalar Q1 gather with cached connectivity and stack scratch. Keep
// MFEM's basis evaluation and physical Jacobian so evaluation order and mesh
// geometry are unchanged. Other spaces, transformed DOFs, and updated spaces
// use the ordinary GridFunction operations.
class Q1ScalarGather {
public:
    explicit Q1ScalarGather(const mfem::GridFunction& field)
        : field_(field), space_(field.FESpace()),
          mesh_sequence_(space_->GetMesh()->GetSequence()),
          space_sequence_(space_->GetSequence()) {
        const auto* mesh=space_->GetMesh();
        if(space_->GetVDim()!=1 || mesh->Dimension()!=2 || mesh->SpaceDimension()!=2) return;
        std::vector<std::array<int,4>> connectivity(mesh->GetNE());
        mfem::Array<int> dofs;
        mfem::DofTransformation transform;
        for(int e=0;e<mesh->GetNE();++e) {
            const auto* fe=space_->GetFE(e);
            if(e==0) element_=fe;
            else if(fe!=element_) return;
            if(!dynamic_cast<const mfem::H1_QuadrilateralElement*>(fe) ||
               fe->GetOrder()!=1 || fe->GetDof()!=4 ||
               fe->GetMapType()!=mfem::FiniteElement::VALUE) return;
            space_->GetElementDofs(e,dofs,transform);
            if(dofs.Size()!=4 || !transform.IsIdentity()) return;
            for(int j=0;j<4;++j) {
                if(dofs[j]<0 || dofs[j]>=field_.Size()) return;
                connectivity[e][j]=dofs[j];
            }
        }
        dofs_=std::move(connectivity);
    }

    bool available() const {
        return !dofs_.empty() && field_.FESpace()==space_ &&
            mesh_sequence_==space_->GetMesh()->GetSequence() &&
            space_sequence_==space_->GetSequence();
    }

    double value(int element,const mfem::IntegrationPoint& point) const {
        if(!available()) return field_.GetValue(element,point);
        double shape_data[4],local_data[4];
        mfem::Vector shape(shape_data,4),local(local_data,4);
        load(element,local_data);
        element_->CalcShape(point,shape);
        return shape*local;
    }

    std::array<double,2> gradient(int element,const mfem::IntegrationPoint& point) const {
        auto* transformation=field_.FESpace()->GetMesh()->GetElementTransformation(element);
        transformation->SetIntPoint(&point);
        if(!available()) {
            mfem::Vector gradient;
            field_.GetGradient(*transformation,gradient);
            return {gradient[0],gradient[1]};
        }
        std::array<double,2> result{};
        mfem::Vector gradient(result.data(),2);
        double derivative_data[8],local_data[4],reference_data[2];
        mfem::DenseMatrix derivatives(derivative_data,4,2);
        mfem::Vector local(local_data,4),reference(reference_data,2);
        load(element,local_data);
        element_->CalcDShape(point,derivatives);
        derivatives.MultTranspose(local,reference);
        transformation->InverseJacobian().MultTranspose(reference,gradient);
        return result;
    }

private:
    void load(int element,double* local) const {
        const double* values=field_.HostRead();
        for(int j=0;j<4;++j) local[j]=values[dofs_[element][j]];
    }
    const mfem::GridFunction& field_;
    const mfem::FiniteElementSpace* space_;
    const mfem::FiniteElement* element_=nullptr;
    long mesh_sequence_,space_sequence_;
    std::vector<std::array<int,4>> dofs_;
};
} // namespace raikou::mfem_backend
