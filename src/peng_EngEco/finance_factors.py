class FinanceFactors:
    """
    A collection of common Engineering Economics time value of money factors.
    All formulas assume:
        i = interest rate per period (decimal, e.g., 0.07 for 7%)
        n = number of periods
    """

    @staticmethod
    def F_P(i, n):
        """
        F/P Factor — Single Payment Compound Amount Factor
        Converts Present Value (P) to Future Value (F).

        Formula:
            F/P = (1 + i)^n
        """
        return (1 + i) ** n

    @staticmethod
    def P_F(i, n):
        """
        P/F Factor — Single Payment Present Worth Factor
        Converts Future Value (F) to Present Value (P).

        Formula:
            P/F = 1 / (1 + i)^n
        """
        return 1 / ((1 + i) ** n)

    @staticmethod
    def F_A(i, n):
        """
        F/A Factor — Uniform Series Compound Amount Factor (USCAF)
        Converts a uniform annual amount (A) to a Future Value (F).

        Formula:
            F/A = ((1 + i)^n - 1) / i
        """
        return ((1 + i) ** n - 1) / i

    @staticmethod
    def A_F(i, n):
        """
        A/F Factor — Sinking Fund Factor
        Converts Future Value (F) to a uniform annual amount (A).

        Formula:
            A/F = i / ((1 + i)^n - 1)
        """
        return i / ((1 + i) ** n - 1)

    @staticmethod
    def P_A(i, n):
        """
        P/A Factor — Series Present Worth Factor
        Converts a uniform annual amount (A) to Present Value (P).

        Formula:
            P/A = ((1 + i)^n - 1) / (i * (1 + i)^n)
        """
        return ((1 + i) ** n - 1) / (i * (1 + i) ** n)

    @staticmethod
    def A_P(i, n):
        """
        A/P Factor — Capital Recovery Factor
        Converts Present Value (P) to a uniform annual amount (A).

        Formula:
            A/P = i(1 + i)^n / ((1 + i)^n - 1)
        """
        return (i * (1 + i) ** n) / ((1 + i) ** n - 1)

    # ---------------- GRADIENT FACTORS ---------------- #

    @staticmethod
    def P_G(i, n):
        """
        P/G Factor — Arithmetic Gradient Present Worth Factor
        Converts a linear gradient (G per period) to Present Value (P).

        Formula:
            P/G = (1 / i) * [((1 + i)^n - 1) / (i * (1 + i)^n) - n / (1 + i)^n]
        """
        return (1 / i) * (((1 + i) ** n - 1) / (i * (1 + i) ** n) - n / ((1 + i) ** n))

    @staticmethod
    def A_G(i, n):
        """
        A/G Factor — Arithmetic Gradient Uniform Series Factor
        Converts a linear gradient (G per period) to an equivalent uniform annual series (A).

        Formula:
            A/G = (1 / i) - (n / ((1 + i)^n - 1))
        """
        return (1 / i) - (n / ((1 + i) ** n - 1))
    
    # ---------- PERIOD-AGNOSTIC CORE CONVERTERS ----------

    @staticmethod
    def eff_from_nominal(i_nom: float, m_in: int, k_out: int = 1) -> float:
        """
        Effective rate per target period (k_out per year) from nominal i_nom compounded m_in/year.
        j = (1 + i_nom/m_in)^(m_in / k_out) - 1
        k_out=1 -> EAR; k_out=12 -> effective monthly; k_out=4 -> effective quarterly.
        """
        if m_in <= 0 or k_out <= 0:
            raise ValueError("m_in and k_out must be >= 1")
        return (1 + i_nom / m_in) ** (m_in / k_out) - 1

    @staticmethod
    def eff_convert(i_eff: float, k_in: int, k_out: int) -> float:
        """
        Convert an effective rate per k_in periods/year to per k_out periods/year.
        j = (1 + i_eff)^(k_in / k_out) - 1
        """
        if k_in <= 0 or k_out <= 0:
            raise ValueError("k_in and k_out must be >= 1")
        return (1 + i_eff) ** (k_in / k_out) - 1

    @staticmethod
    def nominal_from_eff(i_eff: float, k_in: int, m_out: int) -> float:
        """
        Equivalent nominal rate compounded m_out/year for a given effective rate per k_in periods/year.
        Steps: i_annual = (1 + i_eff)^(k_in) - 1
               i_nom_out = m_out * ((1 + i_annual)^(1/m_out) - 1)
        """
        if k_in <= 0 or m_out <= 0:
            raise ValueError("k_in and m_out must be >= 1")
        i_annual = (1 + i_eff) ** k_in - 1
        return m_out * ((1 + i_annual) ** (1 / m_out) - 1)

    # ---------- LIGHTWEIGHT WRAPPERS (keep for students) ----------

    @staticmethod
    def eff_annual_from_nominal(i_nom: float, m: int) -> float:
        """EAR from nominal i_nom compounded m/year."""
        return FinanceFactors.eff_from_nominal(i_nom, m, k_out=1)

    @staticmethod
    def per_period_from_nominal(i_nom: float, m: int) -> float:
        """Effective per-compounding-period rate (use only if cash flows match compounding)."""
        return FinanceFactors.eff_from_nominal(i_nom, m, k_out=m)

    @staticmethod
    def eff_annual_from_periodic(i_p: float, periods_per_year: int) -> float:
        """EAR from per-period effective rate i_p."""
        return FinanceFactors.eff_convert(i_p, k_in=periods_per_year, k_out=1)

    @staticmethod
    def per_period_from_eff_annual(i_eff: float, periods_per_year: int) -> float:
        """Per-target-period effective rate from EAR."""
        return FinanceFactors.eff_convert(i_eff, k_in=1, k_out=periods_per_year)
    
    #----------------------LOAN RELATED-----------------------------
    @staticmethod
    def OBF(i: float, n: int, k: int) -> float:
        """
        Outstanding Balance Factor (OBF).

        Computes the remaining loan balance factor after 'k' payments
        of an amortized loan.

        Formula:
            B_k = ( (1+i)^(n-k) - 1 ) / ( i * (1+i)^(n-k) )

        Where:
            B_k : Outstanding balance factor (multiply by payment A to get balance)
            i   : Periodic interest rate (decimal form, e.g. 0.01 for 1%)
            n   : Total number of payments
            k   : Number of payments already made

        Usage:
            balance = A * factors.OBF(i, n, k)
        """
        return ((1 + i) ** (n - k) - 1) / (i * (1 + i) ** (n - k))