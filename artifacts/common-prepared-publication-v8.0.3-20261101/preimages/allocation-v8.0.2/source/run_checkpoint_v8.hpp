#pragma once
// Coordinator v8 whole-run checkpoint (audit items A05/A06): files, manifest and recovery.
//
// Layout of one checkpoint (directory <root>/ckpt-<step, 12 digits>):
//   native.json          native full-time checkpoint (PicArena::save; native's own format and envelope)
//   hall/                Hall owner checkpoint (Owner::save or ResidentOwner::save; Hall's own format and
//                        seals: field words, complete charge ledger incl. the admitted Gauss ledger, heat
//                        slices, kinetic ledgers)
//   source-<i>.bin       source owner v7 component checkpoint of slice owner i (own format and checksum)
//   arrival.bin          arrival owner v1.1 checkpoint (own format and checksum; RNG position = next step and
//                        ordinal: the arrival RNG is counter based)
//   chemistry.json       parent-coverage v2 subflow checkpoint (scheduler clock and product RNG states)
//   coordinator.json     coordinator state: step, accepted source seal and owner digests, receipt chain head,
//                        ledger running sums (exact), controller state (scalars, rate samples), run statistics
//   controller-ring.bin  controller LTE history: up to four accepted (t, E, B) states, binary64
//   manifest.json        written last: {"body": {...}, "body_sha256": sha256(body.dump())}; body: schema, step,
//                        time, bindings, owner versions, run configuration and options, the recorded rank count,
//                        every file with its byte size and SHA-256, every owner's state digests after the step,
//                        the receipt chain head, the genesis checkpoint reference.
// Crash safety: every file is written into <root>/ckpt-<step>.tmp, flushed and fsync'ed; the manifest is
// written to manifest.json.tmp, fsync'ed and renamed; the directory is fsync'ed and renamed to
// ckpt-<step>, and the parent directory fsync'ed. The rename is the commit point: a directory without the
// .tmp suffix is complete. Recovery takes the newest ckpt-<step> directory (temporary directories are
// ignored and reported); its manifest and every file digest must verify or the restore is refused.
// Owner formats are those of the owners (nothing re-encoded); the coordinator binds them by digest.
#include "identity_v8.hpp"
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <unistd.h>

namespace raikou::native_water::v8::checkpoint {
namespace fs = std::filesystem;
inline constexpr char manifest_schema[] = "coordinator-v8-run-checkpoint-v1";

inline void fsync_path(const fs::path& p, bool directory) {
  const int fd = ::open(p.c_str(), directory ? (O_RDONLY | O_DIRECTORY) : O_RDONLY);
  if (fd < 0) throw std::runtime_error("coordinator v8 checkpoint: open for fsync " + p.string() + ": " + std::strerror(errno));
  const int rc = ::fsync(fd);
  ::close(fd);
  if (rc != 0) throw std::runtime_error("coordinator v8 checkpoint: fsync " + p.string() + ": " + std::strerror(errno));
}
inline void fsync_tree(const fs::path& root) {
  std::vector<fs::path> dirs{root};
  for (const auto& e : fs::recursive_directory_iterator(root)) {
    if (e.is_regular_file()) fsync_path(e.path(), false);
    else if (e.is_directory()) dirs.push_back(e.path());
  }
  for (auto it = dirs.rbegin(); it != dirs.rend(); ++it) fsync_path(*it, true);
}
inline void write_file(const fs::path& p, const std::string& bytes) {
  {
    std::ofstream f(p, std::ios::binary | std::ios::trunc);
    f.write(bytes.data(), std::streamsize(bytes.size()));
    f.close();
    if (!f) throw std::runtime_error("coordinator v8 checkpoint: write " + p.string());
  }
  fsync_path(p, false);
}
inline std::string read_file(const fs::path& p, std::uintmax_t cap = std::uintmax_t(1) << 33) {
  if (!fs::is_regular_file(p)) throw std::runtime_error("coordinator v8 checkpoint: missing file " + p.string());
  if (fs::file_size(p) > cap) throw std::runtime_error("coordinator v8 checkpoint: file above its byte cap " + p.string());
  std::ifstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("coordinator v8 checkpoint: read " + p.string());
  return {std::istreambuf_iterator<char>(f), {}};
}
inline std::string file_sha256(const fs::path& p) { return sha(read_file(p)); }
// Every regular file below dir (relative path, size, SHA-256), sorted.
inline Json file_table(const fs::path& dir) {
  std::vector<fs::path> files;
  for (const auto& e : fs::recursive_directory_iterator(dir))
    if (e.is_regular_file() && e.path().filename() != "manifest.json") files.push_back(e.path());
  std::sort(files.begin(), files.end());
  Json t = Json::object();
  for (const auto& f : files) t[fs::relative(f, dir).generic_string()] = {{"bytes", fs::file_size(f)}, {"sha256", file_sha256(f)}};
  return t;
}
inline std::string step_name(std::uint64_t step) {
  char b[32];
  std::snprintf(b, sizeof b, "ckpt-%012llu", static_cast<unsigned long long>(step));
  return b;
}
// Commit a fully written temporary checkpoint directory: manifest (temp + fsync + rename), directory
// fsync, directory rename (commit point), parent fsync. hook(k) is a test point (crash injection).
inline fs::path commit(const fs::path& tmp, const fs::path& final_dir, const Json& body, void (*hook)(int) = nullptr) {
  if (hook) hook(1);
  const Json manifest = {{"body", body}, {"body_sha256", sha(body.dump())}};
  write_file(tmp / "manifest.json.tmp", manifest.dump(1) + "\n");
  if (hook) hook(2);
  fs::rename(tmp / "manifest.json.tmp", tmp / "manifest.json");
  fsync_tree(tmp);
  if (hook) hook(3);
  if (fs::exists(final_dir)) throw std::runtime_error("coordinator v8 checkpoint: destination exists " + final_dir.string());
  fs::rename(tmp, final_dir);
  fsync_path(final_dir.parent_path(), true);
  return final_dir;
}
// Verified manifest body of a checkpoint directory: the manifest seal and every recorded file (size and
// digest; no file missing, no file added). Refuses (throws) on any difference.
inline Json verify(const fs::path& dir) {
  const auto text = read_file(dir / "manifest.json", 64u << 20);
  Json m;
  try { m = Json::parse(text); } catch (const std::exception& e) { throw std::runtime_error(std::string("coordinator v8 checkpoint: manifest is not JSON: ") + e.what()); }
  if (!m.is_object() || m.size() != 2 || !m.contains("body") || !m.contains("body_sha256") || !m.at("body_sha256").is_string())
    throw std::runtime_error("coordinator v8 checkpoint: manifest envelope");
  const auto& body = m.at("body");
  if (sha(body.dump()) != m.at("body_sha256").get<std::string>()) throw std::runtime_error("coordinator v8 checkpoint: manifest seal mismatch");
  if (!body.is_object() || body.value("schema", "") != manifest_schema) throw std::runtime_error("coordinator v8 checkpoint: manifest schema");
  const auto actual = file_table(dir);
  const auto& recorded = body.at("files");
  for (const auto& [name, rec] : recorded.items()) {
    if (!actual.contains(name)) throw std::runtime_error("coordinator v8 checkpoint: recorded file missing: " + name);
    if (actual.at(name) != rec) throw std::runtime_error("coordinator v8 checkpoint: file size/digest differs from the manifest: " + name);
  }
  for (const auto& [name, rec] : actual.items()) {
    (void)rec;
    if (!recorded.contains(name)) throw std::runtime_error("coordinator v8 checkpoint: file not in the manifest: " + name);
  }
  return body;
}
// Recovery: the newest committed checkpoint below root (temporary directories are listed in skipped).
inline fs::path newest(const fs::path& root, Json* skipped = nullptr) {
  if (!fs::is_directory(root)) throw std::runtime_error("coordinator v8 checkpoint: no checkpoint directory " + root.string());
  std::vector<fs::path> done;
  for (const auto& e : fs::directory_iterator(root)) {
    const auto n = e.path().filename().string();
    if (!e.is_directory() || n.rfind("ckpt-", 0) != 0) continue;
    if (n.size() > 4 && n.substr(n.size() - 4) == ".tmp") {
      if (skipped) skipped->push_back({{"dir", e.path().string()}, {"reason", "incomplete (temporary directory, never committed)"}});
      continue;
    }
    done.push_back(e.path());
  }
  if (done.empty()) throw std::runtime_error("coordinator v8 checkpoint: no committed checkpoint below " + root.string());
  std::sort(done.begin(), done.end());
  return done.back();
}
}  // namespace raikou::native_water::v8::checkpoint
