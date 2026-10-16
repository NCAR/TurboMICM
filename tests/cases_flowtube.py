"""One-hour flow tube regression, retaining the original absolute error bound."""
import csv
from ffi import *
from harness import case
from cases_analytical import add_user,advance


@case('regression:RosenbrockRegression.FlowTube','regression:BackwardEulerRegression.FlowTube')
def flowtube(id_):
    be = 'Backward' in id_
    path = ROOT/'tests/reference'/('backward_euler_wall_loss.csv' if be else 'rosenbrock_wall_loss.csv')
    with path.open() as stream:
        reader = csv.reader(stream); assert next(reader) == ['time','SOA1','SOA2','a-pinene','O3']
        expected = np.array([[float(v) for v in row] for row in reader])
    with Model(['SOA1','SOA2','a-pinene','O3']) as m:
        rate = rate_default(0); rate.p[0] = 8.8e-17*1e-6*6.02214076e23
        m.reaction(['a-pinene','O3'],{'SOA1':.18,'SOA2':.09},rate)
        add_user(m,['SOA1'],{},.01,'r1'); add_user(m,['SOA2'],{},.05,'r2'); m.finalize()
        with State(m) as s:
            s.y[:,0] = [0,0,8e-8,.00002]; s.conditions(298.15); s.update()
            history = [s.y[:,0].copy()]
            for _ in range(3600): advance(s,1.,'be' if be else 2); history.append(s.y[:,0].copy())
            assert expected.shape == (3601,5)
            near(np.array(history)[:10],expected[:10,1:],rtol=0,atol=1e-10)
