#include "transactions.hpp"
int main(int argc,char**argv){
    const std::string out=argc>1?argv[1]:".";
    auto digest=transaction_tests::run([]{return std::make_unique<ss304_phase_v1_3::CpuEngine>();},std::filesystem::path(out)/"transactions");
    chk::put_text("final_state_sha256",digest);
    return chk::finish("test_transactions",out);
}
