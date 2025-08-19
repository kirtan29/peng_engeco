# peng_engeco

📊 **peng_engeco** is a lightweight Python package that provides commonly used Engineering Economics time value of money factors, including single payment, uniform series, and arithmetic gradient factors.  
It is designed as a study and practice tool for P.Eng Technical Exam (Engineering Economics) but is equally useful for general engineering economics problems.  

---

## Features
- Single Payment Factors  
  - F/P, P/F  
- Uniform Series Factors  
  - F/A, A/F, P/A, A/P  
- Arithmetic Gradient Factors  
  - P/G, A/G  

Easy to use — no dependencies beyond standard Python.  

---

## Installation

Clone this repo and install in editable mode:

pip install -e .

Or install directly from GitHub:

pip install git+https://github.com/kirtan29/peng_engeco.git

---

## Implemented Formulas

| Factor | Name | Formula |
|--------|------|---------|
| F/P | Single Payment Compound Amount | (1+i)^n |
| P/F | Single Payment Present Worth   | 1 / (1+i)^n |
| F/A | Uniform Series Compound Amount | ((1+i)^n - 1)/i |
| A/F | Sinking Fund                   | i / ((1+i)^n - 1) |
| P/A | Series Present Worth           | ((1+i)^n - 1) / (i(1+i)^n) |
| A/P | Capital Recovery               | (i(1+i)^n) / ((1+i)^n - 1) |
| P/G | Arithmetic Gradient Present Worth | (1/i)*(((1+i)^n - 1)/(i(1+i)^n) - n/(1+i)^n) |
| A/G | Arithmetic Gradient to Annual  | (1/i) - (n / ((1+i)^n - 1)) |

---

## Notes

- i = interest rate per period (decimal, e.g., 0.07 for 7%)  
- n = number of periods (integer)  
- All methods are static, so you don’t need to instantiate the class.  

---

## Example Calculation

Problem:  
Find the Present Worth (P) of receiving $10,000 at the end of each year for 5 years, if the interest rate is 8% per year.

Solution using factors:  

- Given: A = 10,000, i = 0.08, n = 5  
- Formula: P = A × (P/A, i, n)  
- Compute:  

P/A = ((1+i)^n - 1) / (i(1+i)^n)  
    = ((1.08^5 - 1) / (0.08 * 1.08^5))  
    ≈ 3.9927  

So,  
P = 10,000 × 3.9927 ≈ $39,927

---

Python Example:

from peng_engeco import FinanceFactors

A = 10000
i = 0.08
n = 5

ff = FinanceFactors()
P_A = ff.P_A(i, n)
P = A * P_A
print("Present Worth:", round(P, 2))

---

## Author

Kirtan Adhikari  
• P.Eng Candidate