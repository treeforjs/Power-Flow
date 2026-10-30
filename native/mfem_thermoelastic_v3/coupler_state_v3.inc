// Included inside the independent v3 Coupler class's public section.
std::string serialize_state() const {
 std::ostringstream f(std::ios::out|std::ios::binary);
 auto put=[&](const auto& x){f.write(reinterpret_cast<const char*>(&x),sizeof(x));};
 auto text=[&](const std::string& x){std::uint64_t n=x.size();put(n);f.write(x.data(),n);};
 const std::uint64_t magic=0x45544845524d3033ULL;put(magic);text(signature_);
 put(accepted.time);put(accepted.revision);std::uint64_t n=accepted.T.size();put(n);
 for(auto x:accepted.T)put(x);n=accepted.inventory.size();put(n);
 for(auto x:accepted.inventory){put(x.hi);put(x.lo);}for(auto x:accepted.cumulative)put(x);
 put(accepted.cumulative_mechanical_J);n=accepted.watermarks.size();put(n);
 for(const auto& x:accepted.watermarks){text(x.first);put(x.second);}
 need(bool(f),"state serialization failed");return f.str();
}
std::string state_digest() const {return thermal_v3_detail::sha256(serialize_state());}
void restore_state(const std::string& bytes) {
 need(accepted.revision==0,"restart requires fresh component");
 std::istringstream f(bytes,std::ios::in|std::ios::binary);
 auto get=[&](auto& x){f.read(reinterpret_cast<char*>(&x),sizeof(x));};
 auto text=[&](){std::uint64_t n=0;get(n);need(bool(f)&&n<1000000,"checkpoint string length");std::string x(n,'\0');f.read(x.data(),n);return x;};
 std::uint64_t magic=0,n=0;get(magic);need(magic==0x45544845524d3033ULL,"v3 mechanical checkpoint magic");
 need(text()==signature_,"checkpoint mesh/material/source definition changed");State x;
 get(x.time);get(x.revision);get(n);need(n==accepted.T.size(),"checkpoint thermal dimensions");
 x.T.resize(n);for(auto& v:x.T)get(v);get(n);need(n==patches.size(),"checkpoint inventory dimensions");
 x.inventory.resize(n);for(auto& v:x.inventory){get(v.hi);get(v.lo);}
 for(auto& v:x.cumulative){get(v);need(finite(v),"nonfinite checkpoint ledger");}
 get(x.cumulative_mechanical_J);need(finite(x.cumulative_mechanical_J),"nonfinite checkpoint mechanical ledger");
 get(n);need(n<100000,"checkpoint producer count");
 for(std::uint64_t i=0;i<n;++i){auto id=text();std::uint64_t v=0;get(v);need(!id.empty()&&v>0&&!x.watermarks.count(id),"checkpoint producer identity");x.watermarks.emplace(id,v);}
 need(bool(f)&&f.peek()==std::char_traits<char>::eof()&&finite(x.time)&&x.time>=accepted.time,"truncated/trailing/invalid checkpoint");
 nodes(x);using std::swap;swap(accepted,x);for(int i=0;i<accepted_field.Size();++i)accepted_field[i]=accepted.T[i];++generation_;
}
void save(const std::string& path) const {
 const auto bytes=serialize_state();std::ofstream f(path+".tmp",std::ios::binary);
 const auto digest=thermal_v3_detail::sha256(bytes);f.write(bytes.data(),bytes.size());f.write(digest.data(),digest.size());f.close();need(bool(f),"checkpoint write failed");
 std::filesystem::rename(path+".tmp",path);
}
void load(const std::string& path) {
 std::ifstream f(path,std::ios::binary);need(bool(f),"checkpoint open failed");
 std::string bytes((std::istreambuf_iterator<char>(f)),std::istreambuf_iterator<char>());
 need(bytes.size()>=64,"truncated checkpoint checksum");const auto payload=bytes.substr(0,bytes.size()-64);
 need(bytes.substr(bytes.size()-64)==thermal_v3_detail::sha256(payload),"checkpoint SHA256 mismatch");restore_state(payload);
}
std::vector<double> canonical_nodal_values(const std::vector<double>& values) const {
 need(values.size()==node_lookup.size(),"canonical nodal dimensions");std::vector<double> out(values.size());
 for(std::size_t i=0;i<out.size();++i)out[i]=values[node_lookup[i]];return out;
}
std::vector<double> nodal_reference_volumes() const {
 std::vector<double> out(node_lookup.size());
 for(std::size_t i=0;i<out.size();++i)out[i]=control_area[node_lookup[i]]*g.height;return out;
}
std::vector<std::array<double,2>> reference_node_coordinates() const {
 std::vector<std::array<double,2>> out;out.reserve(node_lookup.size());
 for(double y:y_edges)for(double x:x_edges)out.push_back({x,y});return out;
}
std::vector<double> cell_energies(const State& state) const {
 nodes(state);const auto T=canonical_nodal_values(state.T);std::vector<double> out(g.nx*g.ny);
 for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i){double value=0;
  for(int dy=0;dy<=1;++dy)for(int dx=0;dx<=1;++dx)value+=material.energy_reference(T[(j+dy)*(g.nx+1)+i+dx]);
  out[j*g.nx+i]=.25*(x_edges[i+1]-x_edges[i])*(y_edges[j+1]-y_edges[j])*g.height*value;
 }return out;
}
void initialize_temperature_cartesian(const std::vector<double>& values) {
 need(accepted.revision==0&&generation_==0&&values.size()==node_lookup.size(),"initial nodal state requires fresh component and full grid");
 State candidate=accepted;for(std::size_t i=0;i<values.size();++i)candidate.T[node_lookup[i]]=values[i];nodes(candidate);
 accepted.T.swap(candidate.T);for(int i=0;i<accepted_field.Size();++i)accepted_field[i]=accepted.T[i];++generation_;
 const std::string raw(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(double));
 signature_+=" initial-nodal-sha256:"+thermal_v3_detail::sha256(raw);
}
