"""
GeoEq v0.0.2 -- Example usage.

Run with:  python examples/quickstart.py
"""

from geoeq.soil import properties as sp


print("GeoEq v0.0.2 -- Soil Property Calculations")
print("=" * 50)

# Given soil data
Gs = 2.65
e = 0.72
w = 0.18

# Unit weights
gamma_d = sp.dry_unit_weight(Gs=Gs, e=e)
gamma_sat = sp.saturated_unit_weight(Gs=Gs, e=e)
gamma_sub = sp.submerged_unit_weight(Gs=Gs, e=e)
gamma_bulk = sp.bulk_unit_weight(Gs=Gs, e=e, S=0.80)

print(f"\nGiven: Gs = {Gs}, e = {e}, w = {w}")
print(f"\nDry unit weight:       {gamma_d:.2f} kN/m3")
print(f"Saturated unit weight: {gamma_sat:.2f} kN/m3")
print(f"Submerged unit weight: {gamma_sub:.2f} kN/m3")
print(f"Bulk unit weight (S=0.80): {gamma_bulk:.2f} kN/m3")

# Phase relations
n = sp.porosity_from_void_ratio(e)
S = sp.degree_of_saturation(w=w, Gs=Gs, e=e)

print(f"\nPorosity:              {n:.4f}")
print(f"Degree of saturation:  {S:.4f} ({S * 100:.1f} %)")

# Relative density
e_max = 0.85
e_min = 0.45
Dr = sp.relative_density(e=0.60, e_max=e_max, e_min=e_min)
print(f"\nRelative density (e=0.60, emax={e_max}, emin={e_min}): {Dr:.2f}")

# Atterberg limits
LL = 48.0
PL = 22.0
w_nat = 35.0

PI = sp.plasticity_index(LL=LL, PL=PL)
LI = sp.liquidity_index(w=w_nat, PL=PL, PI=PI)
CI = sp.consistency_index(w=w_nat, LL=LL, PI=PI)

print(f"\nAtterberg limits: LL = {LL}, PL = {PL}, w = {w_nat}")
print(f"Plasticity index:  {PI:.0f} %")
print(f"Liquidity index:   {LI:.3f}")
print(f"Consistency index: {CI:.3f}")
print(f"LI + CI = {LI + CI:.3f} (should be 1.0)")

print(f"\n{'=' * 50}")
print("Done.")
