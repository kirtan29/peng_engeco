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
        return round((1 + i) ** n,5)

    @staticmethod
    def P_F(i, n):
        """
        P/F Factor — Single Payment Present Worth Factor
        Converts Future Value (F) to Present Value (P).

        Formula:
            P/F = 1 / (1 + i)^n
        """
        return round(1 / ((1 + i) ** n),5)

    @staticmethod
    def F_A(i, n):
        """
        F/A Factor — Uniform Series Compound Amount Factor (USCAF)
        Converts a uniform annual amount (A) to a Future Value (F).

        Formula:
            F/A = ((1 + i)^n - 1) / i
        """
        return round(((1 + i) ** n - 1) / i,5)

    @staticmethod
    def A_F(i, n):
        """
        A/F Factor — Sinking Fund Factor
        Converts Future Value (F) to a uniform annual amount (A).

        Formula:
            A/F = i / ((1 + i)^n - 1)
        """
        return round(i / ((1 + i) ** n - 1),5)

    @staticmethod
    def P_A(i, n):
        """
        P/A Factor — Series Present Worth Factor
        Converts a uniform annual amount (A) to Present Value (P).

        Formula:
            P/A = ((1 + i)^n - 1) / (i * (1 + i)^n)
        """
        return round(((1 + i) ** n - 1) / (i * (1 + i) ** n),5)

    @staticmethod
    def A_P(i, n):
        """
        A/P Factor — Capital Recovery Factor
        Converts Present Value (P) to a uniform annual amount (A).

        Formula:
            A/P = i(1 + i)^n / ((1 + i)^n - 1)
        """
        return round((i * (1 + i) ** n) / ((1 + i) ** n - 1),5)

    # ---------------- ARTHEMATIC GRADIENT FACTORS ---------------- #

    @staticmethod
    def P_AG(i, n):
        """
        P/AG Factor — Arithmetic Gradient Present Worth Factor
        Converts a linear gradient (G per period) to Present Value (P).

        Formula:
            P/AG = (1 / i) * [((1 + i)^n - 1) / (i * (1 + i)^n) - n / (1 + i)^n]
        """
        return round((1 / i) * (((1 + i) ** n - 1) / (i * (1 + i) ** n) - n / ((1 + i) ** n)),5)

    @staticmethod
    def A_AG(i, n):
        """
        A/AG Factor — Arithmetic Gradient Uniform Series Factor
        Converts a linear gradient (G per period) to an equivalent uniform annual series (A).

        Formula:
            A/AG = (1 / i) - (n / ((1 + i)^n - 1))
        """
        return round((1 / i) - (n / ((1 + i) ** n - 1)),5)
    
    
    @staticmethod
    def F_AG(i: float, n: int) -> float:
        """
        F/AG — Arithmetic Gradient Future Worth Factor (at year n).
        Converts linear gradient G to future worth at time n.
        FW = G * F_AG(i, n)

        Formula:
            F_AG = (1+i)^n * P_AG
                 = [ ( (1+i)^n - i*n - 1 ) / i^2 ]
        """
        return round(((1 + i) ** n - i * n - 1) / (i ** 2), 5)
    
    #-----------GEOMETRIC GRADIENT------------------------
    @staticmethod
    def P_GG(g: float, i: float, N: int) -> float:
        """
        Geometric Gradient Present Worth factor (P/G).
        
        Returns the factor such that:
            PW = A1 * P_Geo(g, i, N)

        Formula:
            P_Geo = (1 / (i - g)) * [1 - ((1+g)/(1+i))^N],   if i != g
            P_Geo = N / (1+i),                               if i == g

        Parameters
        ----------
        g : float
            Growth rate per period (decimal, e.g., 0.1 for 10%)
        i : float
            Interest rate per period (decimal, e.g., 0.08 for 8%)
        N : int
            Number of periods

        Returns
        -------
        float
            Present Worth factor (rounded to 5 decimals)
        """
        if abs(i - g) < 1e-12:
            return round(N / (1 + i), 5)
        else:
            return round((1 / (i - g)) * (1 - ((1 + g) / (1 + i)) ** N), 5)

    @staticmethod
    def F_GG(g: float, i: float, N: int) -> float:
        """
        Geometric Gradient Future Worth factor (F/G).
        
        Returns the factor such that:
            FW = A1 * F_Geo(g, i, N)

        Formula:
            F_Geo = (1 / (i - g)) * [ (1+i)^N - (1+g)^N ],   if i != g
            F_Geo = N * (1+i)^(N-1),                        if i == g
        """
        if abs(i - g) < 1e-12:
            return round(N * (1 + i) ** (N - 1), 5)
        else:
            return round((1 / (i - g)) * ((1 + i) ** N - (1 + g) ** N), 5)

    @staticmethod
    def A_GG(g: float, i: float, N: int) -> float:
        """
        Geometric Gradient Equivalent Annual Worth factor (A/G).
        
        Returns the factor such that:
            AE = A1 * A_Geo(g, i, N)

        Formula:
            A_Geo = P_Geo(g, i, N) * (A/P, i, N)
            where (A/P, i, N) = i / [1 - (1+i)^(-N)]
        """
        P_geo = FinanceFactors.P_Geo(g, i, N)
        A_P = i / (1 - (1 + i) ** -N)
        return round(P_geo * A_P, 5)
    
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
    
    @staticmethod
    def loan_amortization(loan: float, i: float, n_periods: int, *, decimals: int = 2):
        """
        Compute and print the loan amortization schedule.

        Parameters
        ----------
        loan : float
            Loan principal (present value).
        i : float
            Periodic interest rate (e.g., 0.15 = 15%).
        n_periods : int
            Number of payment periods.
        decimals : int, optional
            Number of decimal places for rounding.

        Returns
        -------
        dict
            {
                "Period": [...],
                "Beginning Balance": [...],
                "Installment": [...],
                "Interest": [...],
                "Principal": [...],
                "Ending Balance": [...],
                "Totals": {"Installment": ..., "Interest": ..., "Principal": ...},
                "Payment": float
            }
        """
        if n_periods <= 0:
            raise ValueError("n_periods must be positive.")

        # --- Compute level payment ---
        if i == 0.0:
            A = loan / n_periods
        else:
            P_over_A = ((1 + i) ** n_periods - 1) / (i * (1 + i) ** n_periods)
            A = loan / P_over_A

        # --- Prepare columns ---
        Period = []
        Begin = []
        Installment = []
        Interest = []
        Principal = []
        End = []

        total_install = total_interest = total_principal = 0.0
        begin = loan

        for k in range(1, n_periods + 1):
            interest = begin * i
            principal = A - interest
            end = begin - principal

            Period.append(k)
            Begin.append(round(begin, decimals))
            Installment.append(round(A, decimals))
            Interest.append(round(interest, decimals))
            Principal.append(round(principal, decimals))
            End.append(round(end, decimals))

            total_install += A
            total_interest += interest
            total_principal += principal
            begin = end

        totals = {
            "Installment": round(total_install, decimals),
            "Interest": round(total_interest, decimals),
            "Principal": round(total_principal, decimals),
        }

        # --- Print table ---
        headers = ["Period", "Beginning Balance", "Installment", "Interest", "Principal", "Ending Balance"]
        widths = [8, 18, 12, 10, 10, 15]

        header_row = "".join(f"{h:<{w}}" for h, w in zip(headers, widths))
        print(f"\nLevel Payment (A): {round(A, decimals)}\n")
        print(header_row)
        print("-" * len(header_row))
        for k in range(n_periods):
            print(f"{Period[k]:<8}{Begin[k]:<18}{Installment[k]:<12}"
                  f"{Interest[k]:<10}{Principal[k]:<10}{End[k]:<15}")
        print("-" * len(header_row))
        print(f"{'TOTAL':<8}{'':<18}{totals['Installment']:<12}{totals['Interest']:<10}{totals['Principal']:<10}")

        return {
            "Period": Period,
            "Beginning Balance": Begin,
            "Installment": Installment,
            "Interest": Interest,
            "Principal": Principal,
            "Ending Balance": End,
            "Totals": totals,
            "Payment": round(A, decimals),
        }
    
    @staticmethod
    def capital_cost_allowances(
        capital: float,
        cca_rate: float,
        years: int,
        *,
        tax_rate: float | None = None,
        close_class_in_year: int | None = None,
        proceeds: float | None = None,
        capital_cost: float | None = None,
        marr: float | None = None,
        print_table: bool = True,
    ) -> dict:
        """
        Compute a detailed CCA schedule under Canada's declining-balance system with:
          • Half-year rule in acquisition year (Year 1)
          • Optional final-year class closure with disposition (recapture / terminal loss)
          • Optional per-year tax shield column
          • Optional PV of shields and PV of salvage (for after-tax PV cost)

        Parameters
        ----------
        capital : float
            Opening UCC at Year 1 (for a single-asset class this equals the asset's capital cost).
        cca_rate : float
            CCA rate (e.g., 0.30 for 30%).
        years : int
            Number of calendar/project years to compute.
        tax_rate : float | None, default None
            Corporate income tax rate (e.g., 0.40). If None, Tax Shield columns are zero/omitted.
        close_class_in_year : int | None, default None
            If provided, indicates the year in which the class becomes empty and is *closed*.
            In that year the normal CCA claim is SKIPPED and recapture/terminal loss is applied.
        proceeds : float | None, default None
            Proceeds of disposition in the closing year (cash received on sale). Required if
            close_class_in_year is provided.
        capital_cost : float | None, default None
            "Lesser-of" cap for disposition deduction. If None, defaults to `capital`.
        marr : float | None, default None
            If provided (per-year discount rate, e.g., 0.15), the function computes PV of tax shields
            and PV of proceeds in the closing year (treated as salvage inflow for (vii)), and
            returns the After-Tax Present Cost = capital - [PV(shields) + PV(salvage)].
        print_table : bool, default True
            If True, prints a formatted table.

        Returns
        -------
        dict
            A dict of column lists including:
            "Year", "OpeningUCC", "AdditionsForBase", "DispositionsCapped", "CCABase",
            "CCA", "ClosingUCC", "Recapture", "TerminalLoss", "TaxShield",
            plus (if marr is not None): "PV_TaxShield", "PV_Salvage", "PV_TaxShield_Total",
            and "AfterTax_PV_Cost".

        Notes
        -----
        • Year 1 uses the half-year rule: CCA Base = 0.5 * Opening UCC.
        • In the closing year (if specified), normal CCA is not claimed. Instead, the class is closed
          and either Recapture (= max(0, Proceeds - OpeningUCC)) or Terminal Loss (= max(0, OpeningUCC - Proceeds))
          is recognized. Closing UCC is set to 0.
        • Tax Shield in year t = tax_rate * (CCA + TerminalLoss - Recapture). If tax_rate is None, it's 0.
        • If marr is provided, PV_TaxShield_t = TaxShield_t / (1 + marr)^t and
          PV_Salvage (only in closing year) = Proceeds / (1 + marr)^t.
        """
        if capital_cost is None:
            capital_cost = capital
        if close_class_in_year is not None and proceeds is None:
            raise ValueError("If close_class_in_year is set, you must provide proceeds.")
        if close_class_in_year is not None and (close_class_in_year < 1 or close_class_in_year > years):
            raise ValueError("close_class_in_year must be within 1..years.")

        # Prepare containers
        Year = []
        OpeningUCC = []
        AdditionsForBase = []
        DispositionsCapped = []
        CCABase = []
        CCA = []
        ClosingUCC = []
        Recapture = []
        TerminalLoss = []
        TaxShield = []
        PV_TaxShield = []

        ucc_open = capital

        for t in range(1, years + 1):
            Year.append(t)
            OpeningUCC.append(ucc_open)

            # Default values each year
            disp_capped = 0.0
            recap = 0.0
            term_loss = 0.0
            cca_t = 0.0
            cca_base_t = 0.0

            # Determine normal-year vs closing-year treatment
            if close_class_in_year is not None and t == close_class_in_year:
                # Class closes this year: skip normal CCA, apply disposition & recapture/terminal loss.
                disp_capped = min(proceeds if proceeds is not None else 0.0, capital_cost)

                # Compare proceeds to opening UCC for recapture or terminal loss
                if proceeds is None:
                    raise ValueError("Proceeds required in the closing year.")
                if proceeds > ucc_open:
                    recap = proceeds - ucc_open
                    term_loss = 0.0
                else:
                    term_loss = ucc_open - proceeds
                    recap = 0.0

                # No normal CCA in closing year when class becomes empty
                cca_t = 0.0
                cca_base_t = 0.0
                ucc_close = 0.0  # class closed

            else:
                # Normal CCA year
                if t == 1:
                    # Half-year rule
                    additions_for_base = 0.5 * ucc_open
                else:
                    additions_for_base = ucc_open

                cca_base_t = additions_for_base
                cca_t = cca_rate * cca_base_t

                # Update closing UCC normally
                disp_capped = 0.0
                ucc_close = ucc_open - cca_t

            # Compute Tax Shield (if tax_rate provided)
            if tax_rate is not None:
                shield_t = tax_rate * (cca_t + term_loss - recap)
            else:
                shield_t = 0.0

            # Discount PV of shield if marr provided
            if marr is not None:
                pv_shield_t = shield_t / ((1.0 + marr) ** t)
            else:
                pv_shield_t = 0.0

            # Append row values
            AdditionsForBase.append(cca_base_t)
            DispositionsCapped.append(disp_capped)
            CCABase.append(cca_base_t)
            CCA.append(cca_t)
            ClosingUCC.append(ucc_close)
            Recapture.append(recap)
            TerminalLoss.append(term_loss)
            TaxShield.append(shield_t)
            PV_TaxShield.append(pv_shield_t)

            # Next year's opening UCC
            ucc_open = ucc_close

        results = {
            "Year": Year,
            "OpeningUCC": OpeningUCC,
            "AdditionsForBase": AdditionsForBase,
            "DispositionsCapped": DispositionsCapped,
            "CCABase": CCABase,
            "CCA": CCA,
            "ClosingUCC": ClosingUCC,
            "Recapture": Recapture,
            "TerminalLoss": TerminalLoss,
            "TaxShield": TaxShield,
        }

        # If discounting requested, compute totals and PV of salvage (proceeds) in closing year
        pv_salvage = 0.0
        if marr is not None:
            results["PV_TaxShield"] = PV_TaxShield
            pv_tax_total = sum(PV_TaxShield)
            results["PV_TaxShield_Total"] = pv_tax_total

            if close_class_in_year is not None and proceeds is not None:
                pv_salvage = proceeds / ((1.0 + marr) ** close_class_in_year)
                results["PV_Salvage"] = pv_salvage

            # After-Tax Present Cost for (vii)
            after_tax_pv_cost = capital - (pv_tax_total + pv_salvage)
            results["AfterTax_PV_Cost"] = after_tax_pv_cost

        if print_table:
            # Build a readable table
            header = (
                f"Capital Cost: {capital:,.2f}  |  CCA rate: {cca_rate*100:.1f}%  |  Years: {years}\n"
            )
            if tax_rate is not None:
                header += f"Tax rate: {tax_rate*100:.1f}%\n"
            if close_class_in_year is not None and proceeds is not None:
                header += f"Closing Year: {close_class_in_year}  |  Proceeds: {proceeds:,.2f}\n"
            if marr is not None:
                header += f"MARR: {marr*100:.1f}%\n"

            print(header)
            cols = (
                f"{'Year':<4} {'Opening UCC':>14} {'Disp (cap)':>12} {'CCA Base':>12} "
                f"{'CCA':>12} {'Closing UCC':>14} {'Recapture':>12} {'TermLoss':>12} "
            )
            if tax_rate is not None:
                cols += f"{'TaxShield':>12} "
            if marr is not None:
                cols += f"{'PVShield':>12} "
            print(cols)
            print('-' * len(cols))

            for i in range(len(Year)):
                row = (
                    f"{Year[i]:<4} {OpeningUCC[i]:>14,.2f} {DispositionsCapped[i]:>12,.2f} {CCABase[i]:>12,.2f} "
                    f"{CCA[i]:>12,.2f} {ClosingUCC[i]:>14,.2f} {Recapture[i]:>12,.2f} {TerminalLoss[i]:>12,.2f} "
                )
                if tax_rate is not None:
                    row += f"{TaxShield[i]:>12,.2f} "
                if marr is not None:
                    row += f"{PV_TaxShield[i]:>12,.2f} "
                print(row)

            # Footer totals
            if marr is not None:
                print('-' * len(cols))
                print(f"PV Tax Shields Total: {results['PV_TaxShield_Total']:,.2f}")
                if pv_salvage:
                    print(f"PV Salvage (Year {close_class_in_year}): {pv_salvage:,.2f}")
                print(f"After-Tax Present Cost: {results['AfterTax_PV_Cost']:,.2f}")

        return results

        
        


    @staticmethod
    def break_even_analysis(cashflow: list[float], years: list[int]) -> None:
        from itertools import accumulate
        """
        Prints a break-even (payback) table and the exact year when cumulative cash flow = 0.

        Parameters
        ----------
        cashflow : list[float]
            Cash flow values for each year. The first value usually includes
            the initial investment (negative).
        years : list[int]
            Year markers corresponding to each cash flow.

        Notes
        -----
        - Uses cumulative (undiscounted) cash flows.
        - Interpolates linearly between years if the break-even occurs
          between two periods.
        - Prints the table directly; does not return values.
        """

        # -------- Validation --------
        if not isinstance(cashflow, list) or not isinstance(years, list):
            raise ValueError("Both 'cashflow' and 'years' must be lists.")
        if len(cashflow) != len(years) or len(cashflow) == 0:
            raise ValueError("'cashflow' and 'years' must be non-empty and of equal length.")
        for i, c in enumerate(cashflow):
            if not isinstance(c, (int, float)):
                raise ValueError(f"cashflow[{i}] is not numeric: {c!r}")
        for i, y in enumerate(years):
            if not isinstance(y, int):
                raise ValueError(f"years[{i}] is not an int: {y!r}")
        for i in range(1, len(years)):
            if years[i] <= years[i - 1]:
                raise ValueError("Years must be strictly increasing.")

        # -------- Compute cumulative --------
        cumulative = list(accumulate(cashflow))

        # -------- Print table --------
        print("-" * 45)
        print(f"{'Year':<6} | {'Cash Flow':<15} | {'Cumulative':<15}")
        print("-" * 45)
        for y, cf, cum in zip(years, cashflow, cumulative):
            print(f"{y:<6} | {cf:<15,.2f} | {cum:<15,.2f}")

        # -------- Find break-even --------
        print("-" * 45)
        idx = None
        for i, val in enumerate(cumulative):
            if val >= 0:
                idx = i
                break

        if idx is None:
            print("\n️  No break-even: cumulative never reaches ≥ 0.")
            return

        if idx == 0:
            print("\n Break-even at year 0 (initially non-negative).")
            return

        # Interpolate between idx-1 and idx
        y0, y1 = years[idx - 1], years[idx]
        c0, c1 = cumulative[idx - 1], cumulative[idx]
        fraction = abs(c0) / (c1 - c0)
        exact_year = y0 + fraction * (y1 - y0)

        print(f"\nBreak-even between year {y0} and {y1}.")
        print(f"Exact break-even year = {exact_year:.4f}")






