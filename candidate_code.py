
import math
import numpy as np

def _alpha_of_T(T, alpha_warm, alpha_cold, T0, k):
    if k <= 0:
        raise ValueError("k must be > 0")
    x = (T - T0) / k
    if x > 50:
        frac = 0.0
    elif x < -50:
        frac = 1.0
    else:
        frac = 1.0 / (1.0 + math.exp(x))
    alpha = alpha_warm + (alpha_cold - alpha_warm) * frac
    return max(0.0, min(1.0, alpha))

def _energy_balance(T, S, sigma, alpha_warm, alpha_cold, T0, k):
    alpha = _alpha_of_T(T, alpha_warm, alpha_cold, T0, k)
    return (1.0 - alpha) * (S / 4.0) - sigma * (T ** 4)

def equilibrium_temperature(S, sigma, alpha_warm, alpha_cold, T0, k, bracket=(210.0,330.0), tol=1e-8, max_iter=200):
    a, b = bracket
    fa = _energy_balance(a, S, sigma, alpha_warm, alpha_cold, T0, k)
    fb = _energy_balance(b, S, sigma, alpha_warm, alpha_cold, T0, k)
    if fa * fb > 0:
        raise ValueError("Bracket does not contain a root.")
    lo, hi = a, b
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fmid = _energy_balance(mid, S, sigma, alpha_warm, alpha_cold, T0, k)
        if abs(fmid) < tol or abs(hi - lo) < tol:
            return mid
        if fa * fmid <= 0:
            hi = mid
            fb = fmid
        else:
            lo = mid
            fa = fmid
    return 0.5 * (lo + hi)

def simulate_temperature_series(C, S0, a, P_days, sigma, alpha_warm, alpha_cold, T0, k, T_init=288.0, dt_days=0.25, n_days=365*6):
    dt_sec = dt_days * 86400.0
    P_sec = P_days * 86400.0
    n_steps = int(n_days / dt_days)

    def S_of_t(t):
        return S0 * (1.0 + a * math.cos(2.0 * math.pi * (t / P_sec)))

    def dTdt(t, T):
        S_t = S_of_t(t)
        alpha = _alpha_of_T(T, alpha_warm, alpha_cold, T0, k)
        incoming = (1.0 - alpha) * (S_t / 4.0)
        outgoing = sigma * T**4
        return (incoming - outgoing) / C

    t_days = np.zeros(n_steps + 1)
    T = np.zeros(n_steps + 1)
    T[0] = T_init

    for i in range(n_steps):
        t0 = i * dt_sec
        Ti = T[i]
        k1 = dTdt(t0, Ti)
        k2 = dTdt(t0 + 0.5*dt_sec, Ti + 0.5*dt_sec*k1)
        k3 = dTdt(t0 + 0.5*dt_sec, Ti + 0.5*dt_sec*k2)
        k4 = dTdt(t0 + dt_sec, Ti + dt_sec*k3)
        T[i+1] = Ti + (dt_sec/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        t_days[i+1] = (i+1)*dt_days

    return t_days, T

def simulate_temperature(C, S0, a, P_days, sigma, alpha_warm, alpha_cold, T0, k, T_init=288.0, dt_days=0.25, n_days=365*6):
    t_days, T = simulate_temperature_series(
        C, S0, a, P_days, sigma, alpha_warm, alpha_cold, T0, k,
        T_init, dt_days, n_days
    )
    final_start = n_days - 365.0
    T_final = T[t_days >= final_start]
    return {
        "T_mean": float(np.mean(T_final)),
        "T_amp": float(0.5*(np.max(T_final) - np.min(T_final)))
    }

def calibrate_T0_k(observations, times_days, C, S0, a, P_days, sigma, alpha_warm, alpha_cold,
                   T_init=288.0, dt_days=0.25, search_bounds=((265.0,285.0),(2.0,10.0))):
    obs = np.asarray(observations)
    t_obs = np.asarray(times_days)

    (T0_lo, T0_hi), (k_lo, k_hi) = search_bounds
    best = (None, None, float("inf"))

    T0_grid = np.linspace(T0_lo, T0_hi, 21)
    k_grid  = np.linspace(k_lo,  k_hi,  21)

    for T0v in T0_grid:
        for kv in k_grid:
            t_days, T = simulate_temperature_series(
                C, S0, a, P_days, sigma, alpha_warm, alpha_cold,
                T0v, kv, T_init, dt_days, max(times_days)+10
            )
            pred = np.interp(t_obs, t_days, T)
            err = pred - obs
            sse = float(np.sum(err*err))
            if sse < best[2]:
                best = (float(T0v), float(kv), sse)

    return {"T0_hat": best[0], "k_hat": best[1], "sse": best[2]}
