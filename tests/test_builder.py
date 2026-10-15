"""Solver-builder validation, runtime policy metadata, and option overrides."""
from functools import partial
from builderffi import *


def reactions():
    return [([0],[(1,1.)],rate(0,[2.15e-11,0,110.,300.,0])),([1],[(2,1.)],rate(0,[3.3e-11,0,55.,300.,0]))]


def missing():
    b=BuilderSpec();out=P()
    assert lib.chem_builder_build(C.byref(b),C.byref(out))==1 and not out


def unused():Builder(list('ABCD'),reactions(),ignore_unused=False).build(status=1)


def build(method=3):
    for maximum in [0,4]:
        b=Builder(list('ABC'),reactions(),maximum_cells=maximum)
        with b.build() as m:
            assert lib.chem_mechanism_phase_name(m.handle)==b'gas'
            assert m.data.species_count==3
            limit=lib.chem_maximum_cells(m.handle)
            assert limit==4 if maximum else limit>1e8


def options(backward=False):
    with Builder(list('ABC'),reactions()).build() as m:
        assert lib.chem_update_rates(m.state)==0
        if backward:
            p=lib.chem_backward_euler_parameters();s=Stats()
            assert lib.chem_solve_backward_euler(m.state,C.byref(p),5.,C.byref(s))==0
            assert (s.final_time,s.function_calls,s.jacobian_updates,s.internal_steps,s.solves)==(5.,2,2,2,2)
            p.small,p.max_iterations=1.,1
            status=lib.chem_solve_backward_euler(m.state,C.byref(p),5.,C.byref(s))
            assert (s.final_time,s.function_calls,s.jacobian_updates,s.internal_steps,s.solves)==(.03125,6,6,6,6)
        else:
            m.set_y([[1.,0.,0.]])
            p=lib.chem_rosenbrock_parameters(3)
            status,s=m.solve(5.,parameters=p)
            assert status==0
            assert (s.final_time,s.function_calls,s.jacobian_updates,s.internal_steps,s.solves)==(5.,18,9,9,27)
            p.h_min,p.max_steps=15.,6
            m.set_y([[1.,0.,0.]])
            status,s=m.solve(5.,parameters=p)
            assert status==0
            assert (s.final_time,s.function_calls,s.jacobian_updates,s.internal_steps,s.solves)==(5.,2,1,1,3)
            assert p.h_min==15. and p.max_steps==6


class Recorder(Definition):
    def configure(self,data,process,cs):self.attach(process)


def finalize(sparse_only=False):
    for maximum in ([4] if sparse_only else [0,4]):
        model=Recorder(used=['A'])
        b=Builder(list('ABC'),reactions(),models=[model],maximum_cells=maximum)
        with b.build() as m:
            assert len(model.finalized)==1
            assert model.finalized[0]==(3,4 if maximum else Z(-1).value)
            assert lib.chem_maximum_cells(m.handle)==(4 if maximum else Z(-1).value)
            b.check()


CASES={
'SolverBuilder.ThrowsMissingSystem':missing,
'SolverBuilder.ThrowsUnusedSpecies':unused,
'SolverBuilder.CanBuildBackwardEuler':partial(build,0),
'SolverBuilder.CanBuildRosenbrock':build,
'SolverBuilder.CanBuildBackwardEulerOverloadedSolverMethod':partial(options,True),
'SolverBuilder.CanBuildRosenbrockOverloadedSolveMethod':options,
'SolverBuilder.FinalizeReceivesDenseMatrixPolicy':finalize,
'SolverBuilder.FinalizeSupportsSparseOnlySignature':partial(finalize,True),
}
