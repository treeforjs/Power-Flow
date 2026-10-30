#include "cuda_batched_force.hpp"
#include "cuda_driver_mock.hpp"
#include <functional>
#include <iostream>

using namespace pic_batch;
namespace io=snapshot_io;
namespace mc=mock_cuda;
namespace {
void check(bool pass,const char* message){io::require(pass,message);}
template<class F>std::string throws(F function,const char* text) {
    try{function();}catch(const std::exception& error){check(std::string(error.what()).find(text)!=std::string::npos,"Unexpected original exception");return error.what();}
    throw std::runtime_error("Expected failure was not raised");
}
}
int main(int argc,char** argv) {
 try {
    check(argc==3,"Expected actual fixture and unchanged PTX path");const auto s=HostSnapshot::load(argv[1],true);
    std::vector<std::string> passed;
    auto assert_borrowed_clean=[&](const CleanupReport& report) {
        check(report.resources_released && !report.resources_retained && report.caller_context_restored,"Resources not cleanly released");
        check(mc::primary_refs()==1 && mc::live_allocations()==0 && mc::unsafe_releases()==0,"Mock detected lifetime or borrowed-primary violation");
        check(mc::calls["FORBIDDEN_context_destroy"]==0 && mc::calls["FORBIDDEN_primary_reset"]==0,"Borrowed context was destroyed/reset");
    };
    mc::reset(false);CleanupReport no_current;
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&no_current});},"cuDSS-first");
    check(mc::calls["device_alloc"]==0 && mc::primary_refs()==0,"Missing current context silently initialized a context");passed.push_back("required shared-primary: no private fallback");
    mc::reset();mc::set_current(mc::foreign());CleanupReport wrong_initial;
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&wrong_initial});},"not device0 shared primary");
    CUcontext current=nullptr;cuCtxGetCurrent(&current);check(current==mc::foreign() && mc::primary_refs()==1 && mc::calls["device_alloc"]==0,"Wrong initial context disturbed caller");passed.push_back("wrong initial context rejected without switching caller");
    mc::reset();CleanupReport retain_failure;mc::fail("retain_primary",1);
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&retain_failure});},"cuDevicePrimaryCtxRetain");
    cuCtxGetCurrent(&current);check(current==mc::primary() && mc::primary_refs()==1 && mc::calls["device_alloc"]==0,"Primary retain failure disturbed caller");passed.push_back("primary retain initialization failure preserves caller");
    mc::reset(false);CleanupReport push_failure;mc::fail("push_current",1);
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::StandalonePrimary,&push_failure});},"cuCtxPushCurrent");
    cuCtxGetCurrent(&current);check(!current && mc::primary_refs()==0 && mc::calls["device_alloc"]==0,"Standalone push failure leaked primary reference");passed.push_back("standalone activation initialization failure balances primary retain");
    mc::reset();CleanupReport success;
    {
        CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&success});
        check(mc::calls["context_sync"]==1 && mc::queued()==0,"Static snapshot was published before completion");
        b.set_snapshot(s);check(mc::calls["context_sync"]==2 && mc::queued()==0,"Static refresh lacks completion fence");
        const auto allocations=mc::calls["device_alloc"]+mc::calls["host_alloc"];
        b.evaluate(s.particles,s.phi,s.parameters,s.reference);b.discard();
        check(mc::calls["context_sync"]==2,"Added a context fence to a normal force step");
        check(allocations==mc::calls["device_alloc"]+mc::calls["host_alloc"],"Force step allocated CUDA buffers");
        b.close();
    }
    assert_borrowed_clean(success);check(success.issue_count==0,"Successful cleanup reported errors");passed.push_back("setup and refresh fenced; no normal-step context fence/allocation");
    for(const auto& [operation,occurrence,message]:std::vector<std::tuple<std::string,int,std::string>>{
        {"h2d",1,"H2D particles"},{"h2d",2,"H2D phi"},{"launch",1,"gather_push_kernel"},
        {"d2h",1,"D2H force results"},{"event_record",1,"event start"},{"event_record",2,"event after H2D"},
        {"event_record",3,"event after kernel"},{"event_record",4,"event after D2H"},
        {"event_sync",1,"final batch synchronization"},{"event_elapsed",1,"event elapsed"},{"kernel_status",1,"GPU batch failed"}}) {
        mc::reset();CleanupReport report;
        {
            CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&report});mc::fail(operation,occurrence);
            const auto original=throws([&]{b.evaluate(s.particles,s.phi,s.parameters,s.reference);},message.c_str());
            check(report.attempted && report.stream_drained && mc::queued()==0,"Error did not drain before leaving evaluate");
            throws([&]{b.evaluate(s.particles,s.phi,s.parameters,s.reference);},"closed");
            (void)original;
        }
        assert_borrowed_clean(report);passed.push_back("drain and preserve error after "+operation+std::to_string(occurrence));
    }
    mc::reset();CleanupReport static_failure;mc::fail("context_sync",1);
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&static_failure});},"static publication fence");
    assert_borrowed_clean(static_failure);check(mc::queued()==0,"Static failure retained an unfenced DMA");passed.push_back("failed initial static fence drains before constructor unwind");
    mc::reset();CleanupReport partial;mc::fail("device_alloc",3);
    throws([&]{CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&partial});},"cuMemAlloc");
    assert_borrowed_clean(partial);passed.push_back("partial constructor resource owner cleans allocations");
    mc::reset();CleanupReport fallback;
    {
        CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&fallback});mc::fail("event_sync",1);mc::fail("stream_sync",1);
        throws([&]{b.evaluate(s.particles,s.phi,s.parameters,s.reference);},"final batch synchronization");
    }
    assert_borrowed_clean(fallback);check(fallback.issue_count==1 && fallback.used_context_fallback,"Cleanup failure was not retained/reported");passed.push_back("cleanup stream failure recorded; context fallback drains; original error preserved");
    mc::reset();CleanupReport quarantine;
    {
        CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&quarantine});mc::fail("event_sync",1);mc::fail("stream_sync",1);mc::fail("context_sync",2);
        throws([&]{b.evaluate(s.particles,s.phi,s.parameters,s.reference);},"final batch synchronization");
    }
    check(quarantine.resources_retained && !quarantine.resources_released && quarantine.issue_count>=2 && mc::live_allocations()>0 && mc::unsafe_releases()==0,"Unproven drain freed async resources");
    check(mc::calls["release_primary"]==0 && mc::calls["FORBIDDEN_primary_reset"]==0,"Quarantine reset or dropped protected context");passed.push_back("failed stream and context drain quarantines resources/reference and reports failures");
    for(bool lost:{false,true}) {
        mc::reset();CleanupReport report;
        {
            CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&report});
            mc::set_current(lost?nullptr:mc::foreign());
            throws([&]{b.evaluate(s.particles,s.phi,s.parameters,s.reference);},"context changed");
            b.close();CUcontext after=nullptr;cuCtxGetCurrent(&after);check(after==(lost?nullptr:mc::foreign()),"Cleanup did not restore caller context");
        }
        assert_borrowed_clean(report);passed.push_back(lost?"lost context fails before dispatch; caller null restored":"wrong context fails before dispatch; caller context restored");
    }
    mc::reset();CleanupReport transaction;
    {
        CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&transaction});
        b.evaluate(s.particles,s.phi,s.parameters,s.reference);auto modified=s.particles;modified.back().weight*=2;const auto before=modified;
        throws([&]{b.commit_velocities(modified);},"state changed");check(std::memcmp(modified.data(),before.data(),before.size()*sizeof(before[0]))==0,"Failed commit partially changed particles");
        b.discard();b.evaluate(s.particles,s.phi,s.parameters,s.reference);auto committed=s.particles;b.commit_velocities(committed);
        for(std::size_t i=0;i<committed.size();++i)check(committed[i].proper_velocity_m_s==io::vec(s.reference[i].u),"Validated transaction changed velocity");
    }
    assert_borrowed_clean(transaction);passed.push_back("mutation-before-commit is atomic; discard and successful commit work");
    mc::reset();CleanupReport activation_failure;
    {
        CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&activation_failure});
        mc::set_current(mc::foreign());mc::fail("push_current",1);b.close();
        cuCtxGetCurrent(&current);check(current==mc::foreign(),"Failed cleanup context activation disturbed caller");
    }
    check(activation_failure.resources_retained && !activation_failure.resources_released && activation_failure.issue_count>0 && mc::unsafe_releases()==0,"Failed cleanup activation freed inaccessible resources");
    passed.push_back("failed cleanup context activation retains resources and caller context");
    mc::reset(false);CleanupReport standalone;
    {CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::StandalonePrimary,&standalone});b.evaluate(s.particles,s.phi,s.parameters,s.reference);b.discard();}
    cuCtxGetCurrent(&current);check(!current && mc::primary_refs()==0 && standalone.resources_released && standalone.primary_reference_released && mc::unsafe_releases()==0,"Explicit standalone primary lifetime is not balanced");passed.push_back("explicit standalone-primary retain/push/drain/pop/release");
    mc::reset();CleanupReport release_failure;
    {CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&release_failure});mc::fail("host_free",1);b.close();}
    check(release_failure.issue_count>0 && release_failure.resources_retained,"Resource-release failure was hidden");passed.push_back("resource-release failure retained in caller-owned cleanup report");
    mc::reset();CleanupReport primary_release_failure;
    {CudaBatchedForce b(argv[2],s,s.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&primary_release_failure});mc::fail("release_primary",1);b.close();}
    check(primary_release_failure.issue_count>0 && primary_release_failure.resources_retained && !primary_release_failure.resources_released && mc::primary_refs()==2,"Primary-reference release failure was hidden");
    passed.push_back("primary-reference release failure retained and reported without reset");
    mc::reset();
    std::cout<<"{\"passed\":true,\"real_CUDA_initialized_or_executed\":false,\"driver_mock_only\":true,\"cases\":[";
    for(std::size_t i=0;i<passed.size();++i)std::cout<<(i?",":"")<<std::quoted(passed[i]);
    std::cout<<"]}\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
