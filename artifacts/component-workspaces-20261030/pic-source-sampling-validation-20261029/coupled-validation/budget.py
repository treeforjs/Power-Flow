"""Pure accounting for the declared campaign caps; never launches or reserves resources."""
STAGES={
    'resource_pilot':dict(runs=3,steps=2000,wall_s=90,duration_ps=10),
    'common_trace_count_isolation':dict(runs=18,steps=20000,wall_s=900,duration_ps=100),
    'native_adaptive_confirmation':dict(runs=18,steps=20000,wall_s=900,duration_ps=100),
}
class CampaignBudget:
    def __init__(self):self.attempts={};self.stage_runs={k:0 for k in STAGES};self.steps=0;self.wall_s=0
    def account_attempt(self,attempt_id,stage,count,seed):
        if stage not in STAGES or count not in [64,32,16]:raise ValueError('Unknown stage/count')
        if not isinstance(seed,int) or isinstance(seed,bool):raise ValueError('Invalid seed')
        if stage=='resource_pilot' and seed!=0:raise ValueError('Pilot seed must be zero')
        if stage=='common_trace_count_isolation' and not 0<=seed<6:raise ValueError('Common-trace seed outside plan')
        if stage=='native_adaptive_confirmation' and not 6<=seed<12:raise ValueError('Adaptive seed outside plan')
        if attempt_id in self.attempts:raise ValueError('Duplicate attempt identity')
        cap=STAGES[stage]
        if self.stage_runs[stage]>=cap['runs'] or len(self.attempts)>=39 or self.steps+cap['steps']>726000 or self.wall_s+cap['wall_s']>32670:
            raise ValueError('Campaign cap exhausted; no automatic extension')
        self.attempts[attempt_id]=dict(stage=stage,count=count,seed=seed)
        self.stage_runs[stage]+=1;self.steps+=cap['steps'];self.wall_s+=cap['wall_s']
        return dict(cap,coupled_execution_enabled=False,GPU_execution_enabled=False)
    def snapshot(self):return dict(accounted_attempts=len(self.attempts),stage_runs=self.stage_runs.copy(),maximum_step_exposure=self.steps,maximum_wall_exposure_s=self.wall_s,launched_runs=0)

def self_test():
    b=CampaignBudget()
    for stage in STAGES:
        seeds=[0] if stage=='resource_pilot' else range(6) if stage=='common_trace_count_isolation' else range(6,12)
        for seed in seeds:
            for count in [64,32,16]:b.account_attempt(f'{stage}-{seed}-{count}',stage,count,seed)
    assert b.snapshot()['accounted_attempts']==39 and b.steps==726000 and b.wall_s==32670
    try:b.account_attempt('fortieth','common_trace_count_isolation',16,0)
    except ValueError:pass
    else:raise AssertionError('Campaign cap expanded')
    try:b.account_attempt('resource_pilot-0-64','resource_pilot',64,0)
    except ValueError:pass
    else:raise AssertionError('Duplicate attempt accepted')
    return dict(passed=True,simulation_only=True,actual_runs_launched=0,steps_cap=726000,wall_cap_s=32670,run_cap=39)
