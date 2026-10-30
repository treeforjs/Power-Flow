// Root assembles the existing FE operators; every MPI rank owns and solves
// a disjoint row block. Field sampling and serialization retain their owner.


namespace distributed_linear
{
enum Command { Finish=0, Create=1, Solve=2, Release=3 };
inline int systems=0, solves=0, next_id=0;
inline std::string output;

inline void Send(Command command,int id)
{
   int message[2]={int(command),id};
   MPI_Bcast(message,2,MPI_INT,0,MPI_COMM_WORLD);
}

struct System
{
   int n=0,rank=0,ranks=1,iterations=0;
   std::vector<int> counts,offsets;
   std::vector<HYPRE_BigInt> partition;
   std::unique_ptr<HypreParMatrix> matrix;
   std::unique_ptr<HypreBoomerAMG> amg;
   std::unique_ptr<CGSolver> cg;
   Vector previous;

   explicit System(const SparseMatrix *source)
   {
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
      if(rank==0)n=source->Height();
      MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
      MFEM_VERIFY(n>=ranks,"Each MPI rank must own at least one matrix row");
      counts.resize(ranks);offsets.resize(ranks+1);partition.resize(ranks+1);
      for(int r=0;r<=ranks;++r){offsets[r]=int((int64_t(n)*r)/ranks);partition[r]=offsets[r];}
      for(int r=0;r<ranks;++r)counts[r]=offsets[r+1]-offsets[r];
      std::vector<int> nonzeros(ranks),starts(ranks);
      if(rank==0)for(int r=0;r<ranks;++r)
      {starts[r]=source->GetI()[offsets[r]];nonzeros[r]=source->GetI()[offsets[r+1]]-starts[r];}
      MPI_Bcast(nonzeros.data(),ranks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(starts.data(),ranks,MPI_INT,0,MPI_COMM_WORLD);
      std::vector<int> rows(counts[rank]+1),columns(nonzeros[rank]);
      std::vector<double> values(nonzeros[rank]);
      MPI_Scatterv(rank==0?source->GetI():nullptr,counts.data(),offsets.data(),MPI_INT,
                   rows.data(),counts[rank],MPI_INT,0,MPI_COMM_WORLD);
      rows.back()=starts[rank]+nonzeros[rank];
      for(int &entry:rows)entry-=starts[rank];
      MPI_Scatterv(rank==0?source->GetJ():nullptr,nonzeros.data(),starts.data(),MPI_INT,
                   columns.data(),nonzeros[rank],MPI_INT,0,MPI_COMM_WORLD);
      MPI_Scatterv(rank==0?source->GetData():nullptr,nonzeros.data(),starts.data(),MPI_DOUBLE,
                   values.data(),nonzeros[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
      std::vector<HYPRE_BigInt> global_columns(columns.begin(),columns.end());
      const auto *ownership=partition.data()+(HYPRE_AssumedPartitionCheck()?rank:0);
      matrix=std::make_unique<HypreParMatrix>(MPI_COMM_WORLD,counts[rank],HYPRE_BigInt(n),
         HYPRE_BigInt(n),rows.data(),global_columns.data(),values.data(),ownership,ownership);
      amg=std::make_unique<HypreBoomerAMG>(*matrix);amg->SetPrintLevel(0);
      cg=std::make_unique<CGSolver>(MPI_COMM_WORLD);cg->SetRelTol(2e-13);cg->SetAbsTol(1e-30);
      cg->SetMaxIter(3000);cg->SetPrintLevel(0);cg->SetPreconditioner(*amg);cg->SetOperator(*matrix);
      previous.SetSize(counts[rank]);previous=0.;cg->iterative_mode=true;
      ++systems;
   }

   void Apply(const Vector *global_rhs,Vector *global_answer)
   {
      Vector rhs(counts[rank]);Vector& answer=previous;
      MPI_Scatterv(rank==0?global_rhs->GetData():nullptr,counts.data(),offsets.data(),MPI_DOUBLE,
                   rhs.GetData(),counts[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
      cg->Mult(rhs,answer);MFEM_VERIFY(cg->GetConverged(),"MFEM distributed linear system did not converge");
      iterations=std::max(iterations,cg->GetNumIterations());
      if(rank==0)global_answer->SetSize(n);
      MPI_Gatherv(answer.GetData(),counts[rank],MPI_DOUBLE,
                  rank==0?global_answer->GetData():nullptr,counts.data(),offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
      ++solves;
   }
};

inline void Workers()
{
   std::map<int,std::unique_ptr<System>> active;
   for(;;)
   {
      int message[2];MPI_Bcast(message,2,MPI_INT,0,MPI_COMM_WORLD);
      if(message[0]==Finish){MFEM_VERIFY(active.empty(),"Unreleased distributed system");break;}
      const int id=message[1];
      if(message[0]==Create)
      {MFEM_VERIFY(!active.count(id),"Duplicate distributed system");active[id]=std::make_unique<System>(nullptr);}
      else if(message[0]==Solve)
      {MFEM_VERIFY(active.count(id),"Missing distributed solve owner");active.at(id)->Apply(nullptr,nullptr);}
      else if(message[0]==Release)
      {MFEM_VERIFY(active.erase(id)==1,"Missing distributed release owner");}
      else{MFEM_ABORT("Unknown distributed command");}
   }
}

inline void Record()
{
   int local[2]={systems,solves},low[2],high[2];
   MPI_Allreduce(local,low,2,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
   MPI_Allreduce(local,high,2,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   MFEM_VERIFY(low[0]==high[0]&&low[1]==high[1],"MPI ranks did not participate in the same systems and solves");
   if(Mpi::WorldRank()==0&&!output.empty())
   {
      std::ofstream record(output+"/parallel_execution.json");MFEM_VERIFY(record,"Cannot write MPI execution record");
      record<<"{\"mpi_ranks\":"<<Mpi::WorldSize()<<",\"assembly_rank\":0,\"output_rank\":0,"
            <<"\"operator_ownership\":\"disjoint_contiguous_rows\",\"linear_systems_per_rank\":"<<low[0]
            <<",\"linear_solves_per_rank\":"<<low[1]<<",\"all_ranks_participated\":true,\"completed\":true}\n";
   }
}
}

struct LinearSystem
{
   SparseMatrix &matrix;
   int id,iterations=0;double residual=0;
   std::unique_ptr<distributed_linear::System> distributed;
 std::unique_ptr<CudssPoisson> cuda;
   explicit LinearSystem(SparseMatrix &a):matrix(a),id(distributed_linear::next_id++)
   {
 const char* backend=std::getenv("RAIKOU_PIC_FIELD_BACKEND");
 if(backend && std::string(backend)=="cudss") {
  MFEM_VERIFY(Mpi::WorldSize()==1,"CUDA benchmark requires exactly one CPU rank");
  cuda=std::make_unique<CudssPoisson>(a);++distributed_linear::systems;
 } else {distributed_linear::Send(distributed_linear::Create,id);distributed=std::make_unique<distributed_linear::System>(&a);}
 }
   ~LinearSystem(){if(cuda)cuda.reset();else{distributed_linear::Send(distributed_linear::Release,id);distributed.reset();}}
   void Solve(const Vector &rhs,Vector &answer)
   {
      pic_bench::Timer timer("linear_solve_with_residual_check");
 if(cuda){cuda->Solve(rhs,answer);++distributed_linear::solves;}else{distributed_linear::Send(distributed_linear::Solve,id);distributed->Apply(&rhs,&answer);}
      Vector error(rhs.Size());matrix.Mult(answer,error);error-=rhs;
      const double norm=rhs.Norml2(),relative=norm>0?error.Norml2()/norm:error.Norml2();
      if(distributed)iterations=std::max(iterations,distributed->iterations);residual=std::max(residual,relative);
 if(cuda)MFEM_VERIFY(error.Norml2()<=std::max(1e-30,2e-13*norm),"CUDA Poisson residual guard failed");
      MFEM_VERIFY(std::isfinite(residual),"Nonfinite solver residual");
   }
};
