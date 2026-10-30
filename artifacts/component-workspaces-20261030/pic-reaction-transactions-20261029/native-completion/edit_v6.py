from pathlib import Path
import re
p=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-completion/source/src/mfem_pic_arena.cpp');s=p.read_text();at=s.index('void require(');s=s[:at]+'''constexpr std::size_t maximum_native_carriers=100000;
constexpr std::size_t maximum_journal_entries=4096;
constexpr std::size_t maximum_journal_bytes=16*1024*1024;
'''+s[at:]
s=s.replace('{"previous_dt_s", s.previous_dt_s}};','{"previous_dt_s", s.previous_dt_s},{"stochastic_completions",s.stochastic_completions},{"journal_bytes",s.journal_bytes}};').replace('require(j.size() == 9, "metadata schema");','require(j.size() == 11, "metadata schema");')
s=s.replace('  s.previous_dt_s = j.at("previous_dt_s");','  s.previous_dt_s = j.at("previous_dt_s");\n  s.stochastic_completions=j.at("stochastic_completions");s.journal_bytes=j.at("journal_bytes");require(s.stochastic_completions.is_object()&&s.journal_bytes<=maximum_journal_bytes,"completion/history byte schema");')
a=s.index('void journal_append(');b=s.index('void journal_validate(',a);s=s[:a]+'''void journal_append(ArenaSnapshot& s,ArenaJson audit){
 require(s.journal.size()<maximum_journal_entries,"bounded prototype journal entry ceiling; validated anchors required for longer runs");
 audit["previous_chain"]=s.chain;const auto chain=sha_bytes(audit.dump());ArenaJson row={{"payload",audit},{"sha256",chain}};const auto bytes=row.dump().size();require(s.journal_bytes<=maximum_journal_bytes&&bytes<=maximum_journal_bytes-s.journal_bytes,"bounded prototype receipt byte ceiling; history is never truncated");s.journal.push_back(std::move(row));s.journal_bytes+=bytes;s.chain=chain;++s.revision;
}
'''+s[b:]
s=s.replace('  std::uint64_t steps = 0;','  std::size_t bytes=0;\n  std::uint64_t steps = 0;')
s=s.replace('  for (const auto &x : s.journal) {','  require(s.journal.size()<=maximum_journal_entries,"bounded history entries");\n  for (const auto &x : s.journal) {\n    const auto n=x.dump().size();require(bytes<=maximum_journal_bytes&&n<=maximum_journal_bytes-bytes,"bounded receipt bytes");bytes+=n;')
s=s.replace('    } else\n      throw std::invalid_argument("journal operation");','    } else if(type!="stochastic_exposure")\n      throw std::invalid_argument("journal operation");')
s=s.replace('chain == s.chain && events == s.event_ids &&','bytes==s.journal_bytes && chain == s.chain && events == s.event_ids &&')
s=s.replace('NativeIndex index;','require(native.size()<=maximum_native_carriers,"bounded prototype native carrier ceiling");\n  NativeIndex index;')
s=s.replace('j.size() <= 1000000','j.size() <= maximum_native_carriers')
s=s.replace('  validate_population(*accepted_);','  require(accepted_->journal.size()<maximum_journal_entries,"bounded prototype history is full");\n  validate_population(*accepted_);')
s=s.replace('  require(accepted_->journal.size()<maximum_journal_entries,"bounded prototype history is full");\n  validate_population(*accepted_);\n  require(pic_->fields_current_, "restart requires current solved fields");','  validate_population(*accepted_);\n  require(pic_->fields_current_, "restart requires current solved fields");')
s=s.replace('"canonical_native_arena_v5"','"canonical_native_arena_v6"')
s=s.replace('    std::ifstream input(path);','    require(std::filesystem::file_size(path)<=64*1024*1024,"bounded arena restart size");\n    std::ifstream input(path);')
pattern=r'PreparedArena\s+PicArena::prepare_compound\(std::size_t maximum_reactions,\s+const std::function<void\(PicArena &\)> &builder\) \{'
rep='''PreparedArena PicArena::prepare_compound(std::size_t maximum_reactions,const std::function<void(PicArena&)>& builder){return compound_impl(maximum_reactions,builder,false);}
PreparedArena PicArena::prepare_compound_with_completion(std::size_t maximum_reactions,const std::function<ArenaStochasticCompletion(PicArena&)>& builder){require(bool(builder),"completion builder required");return compound_impl(maximum_reactions,[&](PicArena& a){auto c=builder(a);a.complete_stochastic(c);},true);}
PreparedArena PicArena::compound_impl(std::size_t maximum_reactions,const std::function<void(PicArena&)>& builder,bool with_completion) {'''
s,n=re.subn(pattern,rep,s);assert n==1
s=s.replace('  validate_population(*accepted_);\n  auto native_before = capture(*pic_);','  require(accepted_->journal.size()+maximum_reactions+(with_completion?1:0)<=maximum_journal_entries,"compound declared resource budget exceeds bounded history");\n  validate_population(*accepted_);\n  auto native_before = capture(*pic_);')
s=s.replace('  compound_limit_ = maximum_reactions;','  compound_limit_ = maximum_reactions;completion_compound_=with_completion;compound_base_chain_=before->chain;')
s=s.replace('accepted_->revision - before->revision <= maximum_reactions,','accepted_->revision - before->revision <= maximum_reactions+(with_completion?1:0),')
s=s.replace('require(accepted_->journal[i].at("payload").at("kind") == "reaction",\n              "nonreaction subevent in collision compound");','require(accepted_->journal[i].at("payload").at("kind") == "reaction" || (with_completion&&i+1==accepted_->journal.size()&&accepted_->journal[i].at("payload").at("kind")=="stochastic_exposure"),\n              "nonreaction subevent in collision compound");\n    require(!with_completion||accepted_->journal.back().at("payload").at("kind")=="stochastic_exposure","missing terminal stochastic completion");')
s=s.replace('    compound_limit_ = 0;','    compound_limit_ = 0;completion_compound_=false;')
s=s.replace('const auto subevents = accepted_->revision - before->revision;','const auto subevents = accepted_->revision - before->revision-(with_completion?1:0);')
s=s.replace('genesis_->event_ids.empty() && genesis_->revision == 0','genesis_->event_ids.empty() && genesis_->stochastic_completions.empty() && genesis_->journal_bytes==0 && genesis_->revision == 0')
p.write_text(s)
