"""Emit rate formulas with fixed parameters and live environmental/custom inputs.

Expressions retain the parenthesization of src/rates.c. Prepared constants use
C double expressions, rather than the build host's transcendental library.
"""
from configure import number


USER_TYPES = {'USER_DEFINED', 'PHOTOLYSIS', 'EMISSION', 'FIRST_ORDER_LOSS'}


def formula(rx, species, offset):
    """Return setup statements (once per call) and lane statements ending in k."""
    def p(key, default):
        return number(rx.get(key, default))

    def arr(prefix='', denominator='300.0'):
        # Fixed zero exponents need no libm call for positive finite conditions.
        result = p(prefix+'A', 1)
        if float(rx.get(prefix+'C', 0)) != 0:
            result += f' * exp({p(prefix+"C", 0)} / T)'
        if float(rx.get(prefix+'B', 0)) != 0:
            result += f' * pow(T / {denominator}, {p(prefix+"B", 0)})'
        return result

    kind = rx['type']
    setup, body = [], []
    if kind in ('ARRHENIUS', 'TAYLOR_SERIES'):
        value = arr(denominator=p('D', 300))
        if float(rx.get('E', 0)) != 0:
            value += f' * (1.0 + {p("E", 0)} * P)'
        if kind == 'TAYLOR_SERIES':
            body += ['double poly = 0.0, tpow = 1.0;']
            for coefficient in rx.get('coefficients', [1]):
                body += [f'poly += {number(coefficient)} * tpow;', 'tpow *= T;']
            value = 'poly * ' + value
    elif kind in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION', 'TERNARY_CHEMICAL_ACTIVATION_JPL19'):
        jpl = kind == 'TERNARY_CHEMICAL_ACTIVATION_JPL19'
        for prefix in (['k0_', 'kinf_', 'kint_'] if jpl else ['k0_', 'kinf_']):
            body += [f'double {prefix[:-1]} = ' + arr(prefix, p(prefix+'D', 298) if jpl else '300.0') + ';']
        body += ['double ratio = k0 * M / kinf;']
        numerator = 'k0' if kind == 'TERNARY_CHEMICAL_ACTIVATION' else 'k0 * M'
        value = (f'{numerator} / (1.0 + ratio) * pow({p("Fc", .6)}, '
                 f'{p("N", 1)} / ({p("N", 1)} + pow(log10(ratio), 2.0)))')
        if jpl:
            body += [f'double kf = {value};']
            value = 'kf + kint * (1.0 - kf / kinf)'
    elif kind == 'TUNNELING':
        value = f'{p("A", 1)} * exp(-({p("B", 0)}) / T + {p("C", 0)} / (T * T * T))'
    elif kind == 'BRANCHED':
        setup += [f'double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp({p("n", 0)});',
                  'double a0 = k0 * (2.45e19 / 6.02214076e23 * 1.0e6);',
                  'double b0 = 0.43 * pow(293.0 / 298.0, -8.0);',
                  'double A0 = a0 / (1.0 + a0 / b0) * pow(0.41, 1.0 / (1.0 + pow(log10(a0 / b0), 2.0)));',
                  f'double z = A0 * (1.0 - {p("a0", 0)}) / {p("a0", 0)};']
        body += ['double a = k0 * M, b = 0.43 * pow(T / 298.0, -8.0);',
                 'double A = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));',
                 f'double pre = {p("X", 1)} * exp(-({p("Y", 0)}) / T);']
        value = 'pre * (A / (A + z))' if rx['branch'] == 'nitrate' else 'pre * (z / (z + A))'
    elif kind == 'REVERSIBLE':
        value = f'{p("A", 1)} * exp({p("C", 0)} / T) * {p("k_r", 0)}'
    elif kind in USER_TYPES:
        value = f'custom[{offset}*cells+c] * {p("scaling factor", 1)}'
    elif kind == 'SURFACE':
        sp = species[rx['reactants'][0]]
        diffusion = number(sp['diffusion coefficient [m2 s-1]'])
        weight = number(sp['molecular weight [kg mol-1]'])
        setup += [f'double speed_factor = 8.0 * CHEM_GAS_CONSTANT / (CHEM_PI * {weight});']
        body += ['double speed = sqrt(speed_factor * T);',
                 f'double radius = custom[{offset}*cells+c], count = custom[{offset+1}*cells+c];']
        value = (f'4.0 * count * CHEM_PI * radius * radius / '
                 f'(radius / {diffusion} + 4.0 / (speed * {p("reaction probability", 1)}))')
    else:
        raise ValueError('Unsupported rate type: ' + kind)
    body += [f'double k = {value};']
    for sp in rx['reactants']:
        if species[sp].get('parameterization') == 'air density':
            body += ['k *= M;']
    return setup, body
