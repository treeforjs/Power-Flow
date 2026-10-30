// Minimal reproducer: the frozen nlohmann JSON header used by native003, the
// V3 parent subflow and the Hall corrective owner allocates inside its
// (noexcept) destructor. A single failed allocation at that point calls
// std::terminate instead of propagating std::bad_alloc.
//
// The pattern below is exactly the one used for every SHA-256 identity in the
// sealed sources: sha256(Json{{"key", value}, ...}.dump()). For each k the
// k-th allocation fails once (in a child process); the parent reports which
// k propagate std::bad_alloc normally and which terminate.
#include "raikou/json.hpp"
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

static long countdown = -1;
void *operator new(std::size_t n) {
  if (countdown == 0) { countdown = -1; throw std::bad_alloc(); }
  if (countdown > 0) --countdown;
  if (void *p = std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc();
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }

static std::string identity(const std::vector<double> &field) {
  using Json = nlohmann::json;
  return Json{{"schema", "example-identity"}, {"field", field}, {"revision", 3}}.dump();
}

int main() {
  const std::vector<double> field{1.0, 2.0, 3.0};
  long propagated = 0, terminated = 0, completed = 0;
  std::string terminated_at;
  for (long k = 0; k < 400; ++k) {
    std::cout.flush();
    const pid_t pid = fork();
    if (pid == 0) {
      std::set_terminate([] { _exit(86); });
      countdown = k;
      try { (void)identity(field); } catch (const std::bad_alloc &) { _exit(1); }
      _exit(countdown >= 0 ? 0 : 2);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    const int code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    if (code == 1) ++propagated;
    else if (code == 86) { ++terminated; terminated_at += (terminated_at.empty() ? "" : ",") + std::to_string(k); }
    else if (code == 0) { ++completed; break; }
  }
  std::cout << "{\"allocations_tested\":" << propagated + terminated + completed
            << ",\"bad_alloc_propagated\":" << propagated << ",\"terminated_inside_noexcept_destructor\":" << terminated
            << ",\"terminating_allocation_indices\":[" << terminated_at << "]}" << std::endl;
  return terminated ? 3 : 0;
}
