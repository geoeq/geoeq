import geoeq as sp
import numpy as np
import matplotlib.pyplot as plt

def run_demo():
    print("--- GeoEq Grain Size Distribution Demo ---")
    
    # 2. Hydrometer Analysis (Raw Data)
    # If 1.6% passed #200 (0.075mm), then hydrometer P must be <= 1.6%
    # But for demo, let's just make the sieve data end at a higher P%
    # and make hydrometer data connect to it.
    
    # New sample: 100g total
    s_open = ["4.75", "2.0", "0.85", "0.425", "0.25", "0.15", "0.075"]
    s_mass_ret = [5, 10, 15, 20, 10, 10, 5] # 75g total
    s_total = 100.0 # 25% passes #200
    
    s_results = sp.sieve_ana(s_open, s_mass_ret, standard="ASTM", total_mass=s_total)
    
    # Hydrometer: handles the 25% passing #200
    # Raw readings for 50g of the passing #200 fraction
    h_readings = [25, 22, 18, 15, 12, 9, 6, 4]
    h_times = [1, 2, 5, 15, 30, 60, 240, 1440]
    h_temp = 20.0
    h_gs = 2.65
    h_ws = 50.0 # g of soil used in hydrometer
    
    h_results = sp.hydro_ana(h_readings, h_times, h_temp, Gs=h_gs, Ws=h_ws)
    
    # SCALE hydrometer P to be relative to TOTAL sample (25% was passing #200)
    f200 = s_results['percent_finer'][-1] / 100.0
    h_results['percent_finer'] = h_results['percent_finer'] * f200
    
    # 3. Combine Data
    combined = sp.process_combined_data(sieve_data=s_results, hydro_data=h_results)
    
    # 4. Get Parameters
    d10 = sp.grain_d10(combined)
    d30 = sp.grain_d30(combined)
    d60 = sp.grain_d60(combined)
    cu = sp.grain_Cu(combined)
    cc = sp.grain_Cc(combined)
    
    print(f"\nCalculated Parameters:")
    print(f"D10: {d10:.4f} mm")
    print(f"D30: {d30:.4f} mm")
    print(f"D60: {d60:.4f} mm")
    print(f"Cu: {cu:.2f}")
    print(f"Cc: {cc:.2f}")
    
    # 5. Plotting (Separate and Combined)
    
    # Combined with separate annotations and advanced smoothing
    sp.grain_size_plot(
        {"Sieve Analysis": s_results, "Hydrometer Analysis": h_results}, 
        smooth=True,
        annotation=True,
        D_para=True, 
        Cu_para=True, 
        Cc_para=True,
        color='black', # Single color for both parts
        marker='*',    # Marker for original data points
        markersize=6,
        linewidth=1.2,
        save_as="grain_size_advanced.png"
    )
    print("Advanced Combined plot saved as 'grain_size_advanced.png'")

if __name__ == "__main__":
    run_demo()
