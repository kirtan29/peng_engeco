# peng_EnggEcon/CashFlowDiagram/CFD.py
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
import math

@dataclass
class CashFlowDiagram:
    """
    Minimal cash-flow diagram utility.
    Conventions: t = 0,1,2,...; +ve inflow, -ve outflow.
    P/G convention: year 1 = 0, year 2 = G, ..., year n = (n-1)G.
    """
    period_label: str = "Year"
    title: Optional[str] = None
    events: List[Tuple[float, float, Optional[str]]] = field(default_factory=list)

    # ---- adders ----
    def add(self, amount: float, t: float, label: Optional[str] = None):
        self.events.append((t, float(amount), label))
        return self

    def add_annuity(self, A: float, n: int, t_start: int = 1, label: str = "A"):
        for k in range(n):
            self.add(A, t_start + k, label)
        return self

    def add_arith_uniform_plus_gradient(self, Y: float, G: float, n: int, t_start: int = 1, label: str = "Y+(k-1)G"):
        # k=1..n ⇒ flows: Y, Y+G, ..., Y+(n-1)G at t_start..t_start+n-1
        for k in range(1, n + 1):
            self.add(Y + (k - 1) * G, t_start + (k - 1), label)
        return self

    def add_arith_pure_gradient(self, G: float, n: int, t_start: int = 1, label: str = "G"):
        # year 1: 0, year 2: G, ..., year n: (n-1)G
        for k in range(n):
            amt = k * G
            if amt != 0.0:
                self.add(amt, t_start + k, label)
        return self

    def add_series(self, flows: Dict[float, float], label: Optional[str] = None):
        for t, a in flows.items():
            self.add(a, t, label)
        return self

    # ---- analysis ----
    def present_worth(self, i: float, base_time: float = 0.0) -> float:
        return sum(a / (1 + i) ** (t - base_time) for t, a, _ in self._sorted())

    def future_worth(self, i: float, at_time: Optional[float] = None) -> float:
        t_ref = self.max_time() if at_time is None else at_time
        return sum(a * (1 + i) ** (t_ref - t) for t, a, _ in self._sorted())

    def equivalent_annual(self, i: float, n: Optional[int] = None, t_start: int = 1) -> float:
        if not self.events:
            return 0.0
        last = int(math.ceil(self.max_time()))
        if n is None:
            n = max(1, last - t_start + 1)
        PW = self.present_worth(i, base_time=t_start - 1)
        A_P = (i * (1 + i) ** n) / ((1 + i) ** n - 1) if i != 0 else 1 / n
        return PW * A_P

    def max_time(self) -> float:
        return max((t for t, _, _ in self.events), default=0.0)

    def _sorted(self):
        return sorted(self.events, key=lambda e: (e[0], -e[1]))

    # ---- plotting ----
    def plot(self, figsize=(9, 3), annotate=True):
        import matplotlib.pyplot as plt

        ev = self._sorted()
        if not ev:
            raise ValueError("No events to plot.")

        ts = [t for t, _, _ in ev]
        as_ = [a for _, a, _ in ev]
        t_min, t_max = min(ts), max(ts)
        t_pad = max(1, int(round(0.05 * max(1, t_max - t_min + 1))))

        plt.figure(figsize=figsize)
        # Axis off (no axes, ticks, or labels)
        plt.axis('off')

        # Baseline
        plt.hlines(0, t_min - t_pad, t_max + t_pad, linewidth=1)

        # Arrow scaling
        Amax = max(abs(a) for a in as_)
        h = 1.0 if Amax == 0 else 1.0 / Amax
        for t, a, lab in ev:
            y = a * h
            # stem
            plt.vlines(t, 0, y, linewidth=2)
            # arrow head
            dy = 0.06 if y >= 0 else -0.06
            plt.arrow(t, y, 0, dy, head_width=0.12, head_length=0.08, length_includes_head=True)
            if annotate:
                txt = f"{a:,.0f}" if abs(a) >= 1000 else f"{a:.2f}"
                if lab:
                    txt = f"{txt}\n{lab}"
                va = "bottom" if a >= 0 else "top"
                ytxt = y + (0.10 if a >= 0 else -0.10)
                plt.text(t, ytxt, txt, ha="center", va=va, fontsize=9)

        if self.title:
            # Title as annotation (axes are off)
            plt.text((t_min + t_max) / 2, 0.9, self.title, ha="center", va="center", fontsize=11)
        plt.tight_layout()
        plt.show()
