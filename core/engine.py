import itertools as it
import math as m

import sympy as sp


class RouthHurwitzCriterion:
    def __init__(self, poly: sp.Expr, poly_sym: sp.Symbol) -> None:
        self.__poly = poly
        self.__poly_sym = poly_sym
        self.__cols: int | None = None
        self.__rows: int | None = None
        self.__degree = sp.degree(self.__poly, self.__poly_sym)
        self.__routh_array: sp.Matrix = sp.Matrix()
        self.__is_stable: bool | None = None
        self.__non_neg_real_part_poles_cnt: int | None = None
        self.__non_neg_real_part_poles_list: list[sp.Expr] = []
        self.__epsilon = sp.Symbol('ε')

    def __eval_limit(self, expr: sp.Expr) -> sp.Expr:
        return sp.limit(expr, self.__epsilon, 0)

    def __eval_sign(self, expr: sp.Expr) -> int:
        return 1 if expr == self.__epsilon else sp.sign(self.__eval_limit(expr))

    def __build_poly(self, start_power: int) -> sp.Expr:
        expr = sp.S(0)
        for power, coeff in zip(self.__gen_power(start_power),
                                self.__routh_array.row(self.__index_of(start_power))):
            expr += coeff * self.__poly_sym**power
        return expr

    def __replace_row(self, start_power: int, poly: sp.Expr) -> None:
        for col, power in zip(range(self.__cols), self.__gen_power(start_power)):
            self.__routh_array[self.__index_of(start_power), col] = poly.coeff(self.__poly_sym, power)

    def __gen_power(self, start: int):
        yield from (start - i * 2 for i in range(self.__cols))

    def __index_of(self, power: int) -> int:
        return self.__degree - power

    def __det(self, x1: int, y1: int, x2: int, y2: int) -> sp.Expr:
        return (self.__routh_array[x1, y1] * self.__routh_array[x2, y2]
                - self.__routh_array[x2, y1] * self.__routh_array[x1, y2])

    def __zero_routh_array(self) -> None:
        self.__cols = m.ceil((self.__degree + 1) / 2)
        self.__rows = self.__degree + 1
        self.__routh_array = sp.zeros(self.__rows, self.__cols)

    def __init_routh_array(self) -> None:
        for (i, p), r in it.product(enumerate(self.__gen_power(self.__degree)), range(min(self.__degree + 1, 2))):
            self.__routh_array[self.__index_of(self.__degree - r), i] = self.__poly.coeff(self.__poly_sym, p - r)

    def __complete_routh_array(self) -> None:
        for p in range(self.__degree - 2, -1, -1):
            v = self.__routh_array[self.__index_of(p + 1), 0]
            for j in range(0, self.__cols - 1):
                coeff = -self.__det(self.__index_of(p + 2), 0, self.__index_of(p + 1), j + 1) / v
                limit = self.__eval_limit(coeff)
                self.__routh_array[self.__index_of(p), j] = \
                    self.__epsilon if limit == 0 and j == 0 \
                    else coeff
            if all(self.__eval_limit(coeff) == 0 for coeff in self.__routh_array.row(self.__index_of(p))):
                aux_poly = self.__build_poly(p + 1)
                aux_poly_diff = aux_poly.diff(self.__poly_sym)
                self.__replace_row(p, aux_poly_diff)

    def __build_routh_array(self) -> None:
        self.__zero_routh_array()
        self.__init_routh_array()
        self.__complete_routh_array()

    def __cnt_sign_changes(self) -> int:
        return sum(self.__eval_sign(a) != self.__eval_sign(b)
                   for a, b in zip(self.__routh_array[1:, 0], self.__routh_array[:-1, 0]))

    def __infer_stability(self) -> None:
        self.__is_stable = all(self.__eval_limit(coeff) > 0 for coeff in self.__routh_array.col(0))

    def __cnt_non_neg_real_part_roots(self) -> None:
        self.__non_neg_real_part_poles_cnt = self.__cnt_sign_changes()

    def __get_non_neg_real_part_roots(self) -> None:
        if self.__degree == 0: return
        self.__non_neg_real_part_poles_list = [sp.nsimplify(pole.evalf(), rational=False).round(5)
                                               for pole in sp.all_roots(self.__poly) if sp.re(pole) > 0]

    def result(self) -> dict:
        return {
            'routh_array': sp.Matrix([sp.Symbol(f'{self.__poly_sym}^{p}') for p in range(self.__degree, -1, -1)])
                           .row_join(self.__routh_array),
            'is_stable': self.__is_stable,
            'non_neg_real_part_pole_cnt': self.__non_neg_real_part_poles_cnt,
            'non_neg_real_part_pools': self.__non_neg_real_part_poles_list,
        }

    def analyze_system(self) -> None:
        self.__build_routh_array()
        self.__infer_stability()
        self.__cnt_non_neg_real_part_roots()
        self.__get_non_neg_real_part_roots()