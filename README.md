
# Albedo Feedback Mini-Benchmark – Q1  
Scientific Computing / Applied Mathematics

---

## Overview

This submission implements a deterministic scientific computing mini-benchmark
based on a zero-dimensional energy balance climate model with
temperature-dependent albedo feedback.

The benchmark evaluates a candidate implementation across three subtasks:

- **Subtask A:** Root-finding for equilibrium temperature  
- **Subtask B:** Time integration with seasonal forcing and metric extraction  
- **Subtask C:** Deterministic calibration of model parameters via grid search  

The benchmark is fully self-contained and requires no external datasets.

---

## Files Included

```
Bench.ipynb
Golden_Solution.ipynb
candidate_code.py
README.md
```

### Bench.ipynb
- Defines the benchmark prompt
- Loads a candidate implementation
- Executes structured unit tests for Subtasks A, B, and C
- Supports optional Gemini-based candidate generation

### Golden_Solution.ipynb
- Reference implementation for all three subtasks
- Includes validation tests
- Provides deterministic baseline results

### candidate_code.py
- Deterministic fallback candidate implementation
- Ensures the benchmark runs end-to-end without external API access

---

## How to Run

1. Open `Bench.ipynb`
2. Run all cells from top to bottom

By default, the notebook loads `candidate_code.py` from the same folder.
You should see:

```
Candidate loaded.
```

followed by unit test output for all subtasks.

---

## Optional: Gemini Agent Mode

Inside `Bench.ipynb`:

```python
USE_GEMINI_AGENT = False
```

If set to:

```python
USE_GEMINI_AGENT = True
```

The notebook will attempt to generate candidate code using the Gemini API.

**Requirements:**
- Valid `GEMINI_API_KEY`
- Internet access
- Active API quota

If Gemini is unavailable or rate-limited, the benchmark automatically falls back
to `candidate_code.py`.

---

## Subtask Details

### Subtask A – Equilibrium Temperature
Solves:

(1 − α(T)) S / 4 − σT⁴ = 0

Using a bracketed bisection method with tolerance control.

### Subtask B – Seasonal Temperature Simulation
Integrates:

C dT/dt = (1 − α(T)) S(t)/4 − σT⁴

Using:
- Fourth-order Runge–Kutta (RK4)
- Seasonal forcing S(t)
- Final-year mean and amplitude extraction

### Subtask C – Parameter Calibration
Estimates:
- T₀ (transition temperature)
- k (transition sharpness)

Using:
- Deterministic two-stage grid search
- Euler integration for calibration speed
- Sum-of-squared-error objective

---

## Design Principles

- Fully deterministic evaluation
- No external datasets
- Uses only `numpy` and `math`
- Modular, testable function structure
- Robust to missing API access

---

## Environment

Tested with:
- Python 3.10+
- NumPy 1.24+

---

This benchmark is designed to evaluate numerical stability, deterministic
behavior, root-finding correctness, ODE integration accuracy, and parameter
calibration logic in a controlled and reproducible setting.
