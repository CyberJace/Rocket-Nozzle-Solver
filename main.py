"""
Main entry point runner for Rocket Nozzle CEA Solver trade studies.
"""

import numpy as np
from src.nozzle_solver import RocketNozzleSolver

def main():
    eps_range = np.linspace(3.0, 100.0, 100)

    # Case 1: LOX / Methane (CH4) Analysis
    print("\n--- Running Analysis: LOX / CH4 ---")
    solver_ch4 = RocketNozzleSolver(
        ox_name="LOX",
        fuel_name="CH4",
        mr=2.8,
        pc_psi=1000.0,
        mdot_kg_s=15.0,
    )
    results_ch4 = solver_ch4.analyze_expansion(eps_range=eps_range, ambient_pa=101325.0)
    solver_ch4.plot_trade_study(results_ch4, output_filename="lox_ch4_sea_level.png")

    print(f"Chamber temperature (Tc): {results_ch4['tc_k']:.2f} K")
    print(f"Mean Molecular Weight (MW): {results_ch4['mw_g_mol']:.2f} g/mol")
    print(f"Throat Area (A*): {results_ch4['A_throat_m2']*1e4:.2f} cm^2")

if __name__ == "__main__":
    main()
