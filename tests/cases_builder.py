from ffi import *
from harness import case,group


def builder(extra=False):
    m=Model(['A','B','C']+(['D'] if extra else []))
    r=rate_default(0);r.p[0],r.p[2]=2.15e-11,110
    m.reaction(['A'],{'B':1},r)
    r.p[0],r.p[2]=3.3e-11,55
    m.reaction(['B'],{'C':1},r)
    return m


@case('solver_builder:SolverBuilder.ThrowsMissingSystem')
def missing(_):
    with Model() as m:assert model_finalize(m.handle,1,0)==-1


@case('solver_builder:SolverBuilder.ThrowsUnusedSpecies')
def unused(_):
    with builder(True) as m:assert model_finalize(m.handle,1,0)==-3


@group('solver_builder:SolverBuilder.','CanBuildBackwardEuler CanBuildRosenbrock')
def build(_):
    with builder() as m:
        copied=model_clone_builder(m.handle);assert copied
        try:
            m.finalize();require(model_finalize(copied,1,0))
            for cells in (1,4,13):
                with State(m,cells) as s:
                    s.set('A',1);s.conditions();s.update()
                    for method in ('be',3):
                        _,stats=s.solve(5,method);assert stats.final_time==5
            assert [model_species_name(copied,i) for i in range(3)]==[b'A',b'B',b'C']
        finally:model_free(copied)


@case('solver_builder:SolverBuilder.CanBuildBackwardEulerOverloadedSolverMethod')
def be_options(_):
    with builder() as m:
        m.finalize()
        with State(m) as s:
            _,stats=s.solve(5,'be')
            assert (stats.final_time,stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.solves)==(5,2,2,2,2)
            p=backward_euler_default();p.small=1;p.max_iterations=1
            _,stats=s.solve(5,'be',p,allow_failure=True)
            assert (stats.final_time,stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.solves)==(.03125,6,6,6,6)


@case('solver_builder:SolverBuilder.CanBuildRosenbrockOverloadedSolveMethod')
def ros_options(_):
    with builder() as m:
        m.finalize()
        with State(m) as s:
            s.set('A',1)
            _,stats=s.solve(5)
            assert (stats.final_time,stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.solves)==(5,18,9,9,27)
            p=rosenbrock_default(3);p.h_min=15;p.max_steps=6
            s.y[:]=[[1],[0],[0]]
            _,stats=s.solve(5,params=p)
            assert (stats.final_time,stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.solves)==(5,2,1,1,3)
            assert p.h_min==15 and p.max_steps==6


@group('solver_builder:SolverBuilder.','FinalizeReceivesDenseMatrixPolicy FinalizeSupportsSparseOnlySignature')
def finalize_callback(_):
    calls=[]
    @C.CFUNCTYPE(None,P,P)
    def callback(model,ctx):
        calls.append((model_species_count(model),model_jacobian_size(model),model_species_index(model,b'A')))
    with builder() as m:
        ext=External();ext.bind=C.cast(callback,P)
        require(model_add_external(m.handle,C.byref(ext)))
        m.finalize()
        assert len(calls)==1 and calls[0][0]==3 and calls[0][1]>=3 and calls[0][2]>=0
