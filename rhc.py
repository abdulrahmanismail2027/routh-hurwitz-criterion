import sys

import sympy as sp
import sympy.parsing.sympy_parser as smp

from core import RouthHurwitzCriterion


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(f'Expected at least two arguments, received: {len(sys.argv) - 1}\n')
        return

    poly_str = sys.argv[1]
    symbol_str = sys.argv[2]

    try:
        symbol = sp.Symbol(symbol_str)
    except Exception as e:
        sys.stderr.write(f'Invalid symbol: {e}')
        return

    try:
        transformations = (*smp.standard_transformations, smp.implicit_multiplication_application, smp.convert_xor)
        poly = smp.parse_expr(poly_str, transformations=transformations)
    except Exception as e:
        sys.stderr.write(f'Invalid polynomial: {e}')
        return

    rhc = RouthHurwitzCriterion(poly, symbol)
    rhc.analyze_system()
    result = rhc.result()

    print(f'System is {'stable' if result['is_stable'] else 'unstable'}!\n')

    if not result['is_stable']:
        print(f'Number of poles with non-negative real part: {result['non_neg_real_part_pole_cnt']}')
        print('The values of the poles with non-negative real part:')
        for pole in result['non_neg_real_part_pools']:
            sp.pprint(pole)
        print()

    print('Routh array:')
    sp.pprint(result['routh_array'])


if __name__ == '__main__':
    main()
