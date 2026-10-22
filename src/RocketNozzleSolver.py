import os
import mapplotlib.pyplot as plt
import numpy as np
from rocketcea.cea_obj import CEA_Obj

class RocketNozzleSolver:
	"""Computes chemical equilibrium and 1D compressible nozzle performance parameters
	parameters across varying expansion ratios and ambient conditions.
	"""

	PSI_TO_PA = 6894.76
	G0 = 9.80665
	UNIVERSAL_GAS_CONSTANT = 8314.4626 # J / (kmol * K)

	def __init__(
		self,
		ox_name: str = "LOX"
		fuel_name: str = "CH4",
		mr: float = 2.8,
		pc_psi: float = 1000.0,
		mdot_kg_s: float = 15.0,
	):
		"""Initializes solver parameters and RocketCEA wrapper instance.

		Args:
		    ox_name (str): Oxidizer name (e.g., 'LOX').
		    fuel_name (str): Fuel name (e.g., 'CH4', 'RP1').
		    mr (float): Mixture ratio (O/F).
		    pc_psi (float): Chamber pressure in PSI
		    mdot_kg_s (float): Propellant mass flow rate in kg/s.
		"""

		self.ox_name = ox_name
		self.fuel_name = fuel_name
		self.mr = mr
		self.pc_psi = pc_psi
		self.pc_pa = pc_psi * self.PSI_TO_PA
		self.mdot = mdot_kg_s

		self.cea = CEA_Obj(oxName=ox_name, fuelName=fuel_name)

	def analyze_expansion(
		self,
		eps_range: np.ndarray = None,
		ambient_pa: float = 101325.0,
	) -> dict:
	    """Executes chemical equilibrium evaluation and 1D compressible flow analysis.

	    Args:
		eps_range (np.ndarray, optional): Array of expansion ratios (Ae/A*).
		ambient_pa (float): Ambient pressure in Pa (Default: 101,325 Pa - Sea 			Level).

	    Returns:
		dict: Performance arrays for mach, exit pressure, exit velocity, thrust, 		and Isp.
 	    """
	    if eps_range is None:
		eps_range = np.linspace(3.0, 100.0, 50)

	    gamma_c, mw_c_g_mol = self.cea.get_Chamber_MolWt_gamma(
		Pc=self.pc_psi, MR=self.mr, eps=eps_range[0]
	    )
	    tc_k = self.cea.get_Temperatures
		Pc=self.pc_psi, MR=self.mr, eps=eps_range[0]
	    )[0]

	    mw_c_kg_mol = mw_c_g_mol / 1000.0
	    R_spec = self.UNIVERSAL_GAS_CONSTANT / mw_c_kg_mol
                                                                                                                                                                                                                                     
	    Gamma = np.sqrt(
		gamma_c
		* (2.0 / (gamma_c + 1.0)) ** ((gamma_c + 1.0) / (gamma_c - 1.0))
	    )
	    A_throat = (self.mdot * np.sqrt(R_spec * tc_k)) / (self.pc_pa * Gamma)

	    results = {
		"eps": [],
		"mach": [],
		"pe_bar": [],
		"ve_m_s: [],
		"thrust_kn": [],
		"isp_s": [],
		"tc_k": tc_k,
		"mw_g_mol": mw_c_g_mol,
		"gamma": gamma_c,
		"A_throat_m2": A_throat,
	    }

	    # Evaluate performance across expansion ratios
	    for eps in eps_range:
		mach_e = self.cea.get_MachNumber(
		    Pc=self.pc_psi, MR=self.mr, eps=eps
		)
		pe_psi = self.cea.get_Pamb(Pc=self.pc_psi, MR=self.mr, eps=eps)
		pe_pa = pe_psi * self.PSI_TO_PA

		A_exit = eps * A_throat
		Te = tc_k / (1.0 + 0.5 * (gamma_c - 1.0) * mach_e**2)
		a_exit = np.sqrt(gamma_c * R_spec * Te)
		ve = mach_e * a_exit

		# 1D Momentum & Pressure Thrust formulation
		F_ momentum = self.mdot * ve
		F_pressure = (pe_pa - ambient_pa) * A_exit
		F_total = F_momentum + F_pressure

		isp = F_total / (self.mdot * self.G0)

		results["eps"].append(eps)
		results["mach"].append(mach_e)
		results["pe_bar"].append(pe_pa / 1e5)
		results["ve_m_s"].append(ve)
		results["thrust_kn"].append(F_total / 1000.0)
		results["isp_s"].append(isp)

	    return results

	def plot_trade_study(
		self,
		results: dict,
		output_filename: str = "nozzle_trade_study.png",
		show_plot: bool = False,
	) -> None:
	    """Generates and saves visual performance trade plots.
		
	       Args:
		   results (dict): Dictionary output from analyze_expansion().
		   output_filename (str): Target filename in outputs/ directory.
		   show_plot (bool): Whether to render Ui window during execution.
	    """
	    fig, ax1 = plt.subplots(figsize=(9, 5.5), dpi=300)

	    color = "tab:blue"
	    ax1.set_xlabel("Expansion Ratio ($\epsilon = A_e / A^*$)", frontsize=11, fontweight="bold)
	    ax1.set_ylabel("Specific Impulse $I_{sp}$ (s)", color=color, fontsize=11, fontweight="bold")
	    line1 = ax1.plot(
		results["eps"],
		results["isp_s"],
		color=color,
		linewidth=2.5,
		label="Specific Impulse ($I_{sp}$)",
	    )
	    ax1.tick_params(axis="y", labelcolor=color)
	    ax1.grid(True, linestyle="--", alpha=0.5)

	    ax2 = ax1.twinx()
	    color = "tab:red"
	    ax2.set_ylabel("Net Thrust $F$ (kN)", color=color, fontsize=11, fontweight="bold")
	    line2 = ax2.plot(
		results["eps"],
		results["thrust_kn"],
		color=color,
		linewidth=2.5,
		linestyle="--",
		label="Net Thrust ($F$)",
	    )
	    ax2.tick_params(axis="y", labelcolor=color)

	    # Combined legend
	    lines = line1 + line2
	    labels = [l.get_label() for l in lines]
	    ax1.legend(lines, labels, loc="lower right", framealpha=0.9)

	    plt.title(
		f"Nozzle Expansion Performance: {self.ox_name}/{self.fuel_name}\n"
		f"($P_c$ = {self.pc_psi} psi, 0/F = {self.mr}, $\dot{{m}}$ = {self.mdot} kg/s)",
		fontsize=12,
		pad=12,
	   )	

	   fig.tight_layout()

	   os.makedirs("outputs", exist_ok=True)
	   filepath = os.path.join("outputs", output_filename)
	   plt.savefig(filepath, dpi=300)
	   print(f"[SUCCESS] Performance trade plot saved: {filepath}")

	   if show_plot:
		plt.show()
	   plt.close()
