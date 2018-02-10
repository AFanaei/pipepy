import cProfile
from sympy.physics import units as U

from equipment.pipe import Pipe
from stream.material_stream import MaterialStream


def main():
    pipe = Pipe(num_nodes=8, length=1 * U.km, teta=0, diameter=0.254 * U.m, epsilon=4.572e-05 * U.m, inlet=MaterialStream(
        P=1761580 * U.pa,
        T=322.737 * U.K,
        m=22.2816 * U.kg / U.s,
        MW=16.0428 * U.g / U.mol
    ), isotherm=False, ambient_t=(10 + 273.15) * U.K, heat_transfer_coef=25 * U.W / ((U.m ** 2) * U.K))

    cProfile.run('pipe.solve_steady_state()', 'profile')


if __name__ == "__main__":
    main()
