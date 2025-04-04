import streamlit as st
import numpy as np
from scipy.linalg import solve
from scipy.special import factorial

st.title("k-out-of-n Maintenance System")

#  Inputs for Part (a)
failure_rate = st.number_input("Failure rate (λ)", min_value=0.0, value=0.1)
repair_rate = st.number_input("Repair rate (μ)", min_value=0.0, value=0.5)
standby_mode = st.selectbox("Standby Mode", ["Warm", "Cold"])
n = st.number_input("Total number of components (n)", min_value=1, value=5)
k = st.number_input("Minimum working components required (k)", min_value=1, max_value=int(n), value=3)
r = st.number_input("Number of repairmen", min_value=1, max_value=int(n), value=2)
ES = 1 / failure_rate  

# Core Logic for Part (a) 

# Cold standby
def machine_repairman(n, k, r, λ, μ): 
    pi = [1.0]  # pi[0] = 1 
    
    # Build pi recursively using balance equations
    for j in range(1, n + 1):
        if j<= n-k+1:
            λ_j_minus_1 = λ * k # k active machines can break
        else:
            λ_j_minus_1 = 0 # System is down

        μ_j = μ * min(j, r) # repair rate
        pi_j = pi[j - 1] * (λ_j_minus_1 / μ_j)
        pi.append(pi_j)

    # Normalize
    total = sum(pi)
    pi = [p / total for p in pi]

    # Uptime = sum of probabilities where at least k components are working
    uptime_fraction = sum(pi[j] for j in range(n - k + 1))

    return pi, uptime_fraction

# Warm standby
def engset_queue(n, k, r, λ, μ):
    
  # Compute π(0) normalization factor
    pi_0_inv = sum(
        (factorial(n) / (factorial(n - j) * factorial(j))) * ((λ / μ) ** j)
        for j in range(r + 1)
    ) + sum(
        (factorial(n) / (factorial(n - j) * factorial(r) * (r ** (j - r)))) * ((λ / μ) ** j)
        for j in range(r + 1, n + 1)
    )
    
    pi_0 = 1 / pi_0_inv

    # Compute π(j) using different formulas for j ≤ s and j > s
    pi = []
    for j in range(n + 1):
        if j <= r:
            prob = (factorial(n) / (factorial(n - j) * factorial(j))) * ((λ / μ) ** j) * pi_0
        else:
            prob = (factorial(n) / (factorial(n - j) * factorial(r) * (r ** (j - r)))) * ((λ / μ) ** j) * pi_0
        pi.append(prob)

    # Compute uptime fraction: sum of probabilities where at least k components are working
    uptime_fraction = sum(pi[j] for j in range(n - k + 1))
    
    return pi, uptime_fraction


# Run calculation for Part (a)
if standby_mode == "Cold":
    pi, uptime_fraction = machine_repairman(n, k, r, failure_rate, repair_rate)
elif standby_mode == "Warm":
    pi, uptime_fraction = engset_queue(n, k, r, failure_rate, repair_rate)

st.subheader("Results for Part (a)")
st.write(f"**Uptime Fraction:** {uptime_fraction:.4f}")

with st.expander("Stationary Distribution"):
    for i, prob in enumerate(pi):
        st.write(f"P({i} failed components) = {prob:.4f}")
       
# Inputs for Part (b)
st.subheader("Inputs for Part (b)")
cost_component = st.number_input("Cost per component", min_value=0.0, value=1.0)
cost_repairman = st.number_input("Cost per repairman", min_value=0.0, value=2.0)
cost_downtime = st.number_input("Downtime cost per unit time", min_value=0.0, value=10.0)
max_components = st.number_input("Maximum number of components to consider (n)", min_value=1, value=10)
max_repairmen = st.number_input("Maximum number of repairmen to consider (r)", min_value=1, value=5)

# Core Logic for Part (b)
def find_optimal_configuration(λ, μ, cost_component, cost_repairman, cost_downtime, standby_mode, max_components, max_repairmen):
    best_config = None
    best_cost = float('inf')

    for k in range(1, max_components):  
        stop_n = False
        for n in range(k + 1, max_components + 1):
            prev_cost = float('inf')
            for r in range(1, max_repairmen + 1):
                if standby_mode == "Cold":
                    pi, uptime_fraction = machine_repairman(n, k, r, λ, μ)
                elif standby_mode == "Warm":
                    pi, uptime_fraction = engset_queue(n, k, r, λ, μ)

                total_cost = cost_component * n + cost_repairman * r + cost_downtime * (1 - uptime_fraction)

                if total_cost < best_cost:
                    best_cost = total_cost
                    best_config = (n, r, k, uptime_fraction)

                if total_cost > prev_cost:
                    break
                prev_cost = total_cost

            if prev_cost > best_cost:
                stop_n = True
                break
        if stop_n and k > 1:
            break

    return best_config, best_cost



# --- Run calculation for Part (b) ---
best_config, best_cost = find_optimal_configuration(failure_rate, repair_rate, cost_component, cost_repairman, cost_downtime, standby_mode, max_components, max_repairmen)

st.subheader("Optimal Configuration for Part (b)")
st.write(f"**Best Number of Components (n):** {best_config[0]}")
st.write(f"**Best Number of Minimum Working Components (k):** {best_config[2]}")
st.write(f"**Best Number of Repairmen (r):** {best_config[1]}")
st.write(f"**Uptime Fraction:** {best_config[3]:.4f}")
st.write(f"**Total Cost per Unit Time:** {best_cost:.4f}")
