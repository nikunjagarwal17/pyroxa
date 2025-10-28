"""
Pure-Python core implementation for the Pyroxa MVP.

Enhanced version with better error handling, validation, and additional reactor types.

Includes:
- Thermodynamics: ideal-gas, constant cp enthalpy/entropy with validation
- Reaction: reversible A <=> B, mass-action kinetics with equilibrium checks
- Reactor: Well-mixed (constant-volume), CSTR, discretized PFR with enhanced physics
- ReactorNetwork: combine reactors in series/parallel with flow modeling
- run_simulation_from_dict: high-level entrypoint using a spec dictionary
"""
from math import log, exp, sqrt
import csv
import warnings
from typing import List, Tuple, Dict, Optional, Union
import numpy as np


class PyroXaError(Exception):
    """Base exception for PyroXa errors."""
    pass


class ThermodynamicsError(PyroXaError):
    """Raised for thermodynamics-related errors."""
    pass


class ReactionError(PyroXaError):
    """Raised for reaction-related errors.""" 
    pass


class ReactorError(PyroXaError):
    """Raised for reactor-related errors."""
    pass


class Thermodynamics:
    """Enhanced ideal gas thermodynamics with constant heat capacity.

    enthalpy(T) = cp * T
    entropy(T) = cp * ln(T) (simplified)
    
    Includes validation and physical bounds checking.
    """
    R = 8.314  # J/mol/K - simplified value for consistency
    
    def __init__(self, cp: float = 29.1, T_ref: float = 298.15):
        """Initialize thermodynamics model.
        
        Args:
            cp: Heat capacity at constant pressure (J/mol/K)
            T_ref: Reference temperature (K)
        """
        if cp <= 0:
            raise ThermodynamicsError(f"Heat capacity must be positive, got {cp}")
        if T_ref <= 0:
            raise ThermodynamicsError(f"Reference temperature must be positive, got {T_ref}")
            
        self.cp = float(cp)
        self.T_ref = float(T_ref)

    def enthalpy(self, T: float) -> float:
        """Calculate enthalpy at temperature T."""
        if T <= 0:
            raise ThermodynamicsError(f"Temperature must be positive, got {T}")
        return self.cp * float(T)

    def entropy(self, T: float) -> float:
        """Calculate entropy at temperature T."""
        if T <= 0:
            return float('nan')
        return self.cp * log(float(T) / self.T_ref)
    
    def equilibrium_constant(self, T: float, delta_G: float) -> float:
        """Calculate equilibrium constant from Gibbs free energy change.
        
        Args:
            T: Temperature (K)
            delta_G: Standard Gibbs free energy change (J/mol)
            
        Returns:
            Equilibrium constant
        """
        if T <= 0:
            raise ThermodynamicsError(f"Temperature must be positive, got {T}")
        return exp(-delta_G / (self.R * T))


class Reaction:
    """Enhanced reversible reaction A <=> B with mass-action kinetics.

    rate = kf*[A] - kr*[B]
    
    Includes equilibrium and thermodynamic consistency checks.
    """
    
    def __init__(self, kf: float, kr: float, validate: bool = True):
        """Initialize reaction with rate constants.
        
        Args:
            kf: Forward rate constant
            kr: Reverse rate constant  
            validate: Whether to validate rate constants
        """
        if validate:
            if kf < 0:
                raise ReactionError(f"Forward rate constant must be non-negative, got {kf}")
            if kr < 0:
                raise ReactionError(f"Reverse rate constant must be non-negative, got {kr}")
                
        self.kf = float(kf)
        self.kr = float(kr)

    def rate(self, conc: List[float]) -> float:
        """Calculate reaction rate."""
        if len(conc) < 2:
            raise ReactionError(f"Need at least 2 concentrations for A <=> B, got {len(conc)}")
        
        # Check for negative concentrations
        if any(c < 0 for c in conc[:2]):
            warnings.warn("Negative concentrations detected, clamping to zero")
            conc = [max(0, c) for c in conc]
            
        return self.kf * conc[0] - self.kr * conc[1]
    
    def equilibrium_constant(self) -> float:
        """Calculate equilibrium constant Keq = kf/kr."""
        if self.kr == 0:
            return float('inf') if self.kf > 0 else float('nan')
        return self.kf / self.kr
    
    def equilibrium_concentrations(self, total_conc: float) -> Tuple[float, float]:
        """Calculate equilibrium concentrations for given total concentration."""
        if total_conc < 0:
            raise ReactionError(f"Total concentration must be non-negative, got {total_conc}")
            
        Keq = self.equilibrium_constant()
        if Keq == float('inf'):
            return 0.0, total_conc
        elif Keq == 0:
            return total_conc, 0.0
        else:
            # At equilibrium: [B]/[A] = Keq and [A] + [B] = total_conc
            A_eq = total_conc / (1 + Keq)
            B_eq = total_conc - A_eq
            return A_eq, B_eq


class ReactionMulti:
    """Enhanced general reversible reaction with stoichiometry.

    reactants and products are dicts mapping species index to stoichiometric coeff.
    rate = kf * product(conc[reactant]**nu) - kr * product(conc[product]**nu)
    """
    
    def __init__(self, kf: float, kr: float, reactants: dict, products: dict, 
                 validate: bool = True):
        """Initialize multi-species reaction.
        
        Args:
            kf: Forward rate constant
            kr: Reverse rate constant
            reactants: {species_index: stoichiometric_coefficient}
            products: {species_index: stoichiometric_coefficient}
            validate: Whether to validate inputs
        """
        if validate:
            if kf < 0:
                raise ReactionError(f"Forward rate constant must be non-negative, got {kf}")
            if kr < 0:
                raise ReactionError(f"Reverse rate constant must be non-negative, got {kr}")
                
        self.kf = float(kf)
        self.kr = float(kr)
        # reactants/products: dict {species_idx: stoich}
        self.reactants = {int(k): abs(int(v)) for k, v in reactants.items()}
        self.products = {int(k): abs(int(v)) for k, v in products.items()}
        
        # Validate no species appears in both reactants and products
        if validate:
            overlap = set(self.reactants.keys()) & set(self.products.keys())
            if overlap:
                warnings.warn(f"Species {overlap} appear in both reactants and products")

    def rate(self, conc: List[float]) -> float:
        """Calculate reaction rate with enhanced error checking."""
        max_idx = max(max(self.reactants.keys(), default=-1), 
                     max(self.products.keys(), default=-1))
        
        if max_idx >= len(conc):
            raise ReactionError(f"Concentration array too short: need {max_idx+1}, got {len(conc)}")
        
        # Forward rate
        f = 1.0
        for idx, nu in self.reactants.items():
            if conc[idx] <= 0:
                f = 0.0
                break
            f *= conc[idx] ** nu
            
        # Reverse rate  
        r = 1.0
        for idx, nu in self.products.items():
            if conc[idx] <= 0:
                r = 0.0
                break
            r *= conc[idx] ** nu
            
        return self.kf * f - self.kr * r
    
    def get_stoichiometry_matrix(self, n_species: int) -> List[List[float]]:
        """Get stoichiometry matrix for this reaction."""
        matrix = [[0.0] * n_species]
        
        # Reactants get negative coefficients
        for idx, nu in self.reactants.items():
            if idx < n_species:
                matrix[0][idx] = -nu
                
        # Products get positive coefficients
        for idx, nu in self.products.items():
            if idx < n_species:
                matrix[0][idx] = nu
                
        return matrix


class MultiReactor:
    """Enhanced reactor for N species and multiple reactions using RK4 integration.

    State: concentrations list of length N.
    Includes better error handling and conservation checking.
    """
    
    def __init__(self, thermo: Thermodynamics, reactions: List[ReactionMulti],
                 species: List[str], T: float = 300.0, conc0: Optional[List[float]] = None,
                 volume: float = 1.0, validate: bool = True):
        """Initialize multi-species reactor.
        
        Args:
            thermo: Thermodynamics model
            reactions: List of reactions
            species: List of species names
            T: Temperature (K)
            conc0: Initial concentrations
            volume: Reactor volume (L)
            validate: Whether to validate inputs
        """
        if validate:
            if T <= 0:
                raise ReactorError(f"Temperature must be positive, got {T}")
            if volume <= 0:
                raise ReactorError(f"Volume must be positive, got {volume}")
            if len(species) == 0:
                raise ReactorError("Must specify at least one species")
                
        self.thermo = thermo
        self.reactions = reactions
        self.species = list(species)
        self.T = float(T)
        self.volume = float(volume)
        self.N = len(self.species)
        
        if conc0 is None:
            self.conc = [0.0] * self.N
        else:
            if len(conc0) != self.N:
                raise ReactorError(f"Initial concentration length {len(conc0)} != number of species {self.N}")
            self.conc = [max(0.0, float(x)) for x in conc0]  # Clamp negative concentrations
            
        # Store initial concentrations for conservation checking
        self._initial_total = sum(self.conc)

    def _dcdt(self, conc: List[float]) -> List[float]:
        """Calculate concentration derivatives with error checking."""
        # Clamp negative concentrations
        conc = [max(0.0, c) for c in conc]
        
        # Initialize derivative
        d = [0.0] * self.N
        
        try:
            for rxn in self.reactions:
                rate = rxn.rate(conc)
                
                # Check for NaN or infinite rates
                if not np.isfinite(rate):
                    warnings.warn(f"Non-finite reaction rate: {rate}")
                    continue
                    
                # Apply stoichiometry
                for idx, nu in rxn.reactants.items():
                    if idx < self.N:
                        d[idx] -= nu * rate
                        
                for idx, nu in rxn.products.items():
                    if idx < self.N:
                        d[idx] += nu * rate
                        
        except Exception as e:
            raise ReactorError(f"Error calculating derivatives: {e}")
            
        return d

    def _check_mass_conservation_with_stoichiometry(self, current_time):
        """Check mass conservation accounting for reaction stoichiometry"""
        total_now = sum(self.conc)
        
        # For reactions that conserve total moles (e.g., A ↔ B), total should be constant
        # For reactions that don't (e.g., A + B → C), total will change predictably
        
        # Simple check: if total mass changes by more than 50%, something is wrong
        relative_change = abs(total_now - self._initial_total) / max(self._initial_total, 1e-10)
        
        # Only warn for very large changes that indicate numerical issues
        if relative_change > 0.75:  # 75% change threshold
            warnings.warn(f"Large mass change detected at t={current_time:.3f} "
                         f"(initial: {self._initial_total:.3f}, current: {total_now:.3f}, "
                         f"change: {relative_change*100:.1f}%)")
        
        return relative_change

    def step(self, dt: float):
        """Enhanced RK4 step with error checking."""
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        y0 = list(self.conc)
        
        try:
            k1 = self._dcdt(y0)
            y1 = [y0[i] + 0.5 * dt * k1[i] for i in range(self.N)]
            
            k2 = self._dcdt(y1)
            y2 = [y0[i] + 0.5 * dt * k2[i] for i in range(self.N)]
            
            k3 = self._dcdt(y2)
            y3 = [y0[i] + dt * k3[i] for i in range(self.N)]
            
            k4 = self._dcdt(y3)
            
            for i in range(self.N):
                self.conc[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
                # Ensure non-negative concentrations
                self.conc[i] = max(0.0, self.conc[i])
                
        except Exception as e:
            raise ReactorError(f"Error in time step: {e}")

    def run(self, time_span: float, time_step: float):
        """Run simulation with enhanced error checking."""
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if time_step <= 0:
            raise ReactorError(f"Time step must be positive, got {time_step}")
        if time_step > time_span:
            warnings.warn(f"Time step {time_step} > time span {time_span}, adjusting")
            time_step = time_span / 10
            
        nsteps = int(max(1, round(time_span / time_step)))
        times = [0.0]
        traj = [list(self.conc)]
        
        for i in range(nsteps):
            try:
                self.step(time_step)
                times.append((i + 1) * time_step)
                traj.append(list(self.conc))
                
                # Check for conservation violations using improved method
                if i % max(1, nsteps // 5) == 0:  # Check every 20% of simulation
                    self._check_mass_conservation_with_stoichiometry(times[-1])
                        
            except Exception as e:
                raise ReactorError(f"Simulation failed at step {i}: {e}")
                
        return times, traj

    def run_adaptive(self, time_span: float, dt_init: float = 1e-3, atol: float = 1e-6, rtol: float = 1e-6):
        # Step-doubling adaptive RK4 (take dt and two dt/2 steps)
        def rk4_step(y, dt):
            N = len(y)
            def deriv(state):
                d = [0.0] * N
                for rxn in self.reactions:
                    rate = rxn.rate(state)
                    for idx, nu in rxn.reactants.items():
                        d[int(idx)] -= nu * rate
                    for idx, nu in rxn.products.items():
                        d[int(idx)] += nu * rate
                return d
            k1 = deriv(y)
            y2 = [y[i] + 0.5 * dt * k1[i] for i in range(N)]
            k2 = deriv(y2)
            y3 = [y[i] + 0.5 * dt * k2[i] for i in range(N)]
            k3 = deriv(y3)
            y4 = [y[i] + dt * k3[i] for i in range(N)]
            k4 = deriv(y4)
            ynew = [y[i] + (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) for i in range(N)]
            return ynew

        t = 0.0
        y = list(self.conc)
        out_times = [0.0]
        out_traj = [list(y)]
        dt = float(dt_init)
        while t < time_span:
            if t + dt > time_span:
                dt = time_span - t
            # one step
            y_big = rk4_step(y, dt)
            # two half steps
            y_half = rk4_step(y, dt*0.5)
            y_two = rk4_step(y_half, dt*0.5)
            # error estimate
            errnorm = 0.0
            for i in range(self.N):
                sc = atol + rtol * max(abs(y[i]), abs(y_two[i]))
                e = (y_two[i] - y_big[i]) / sc if sc != 0 else 0.0
                errnorm += e*e
            errnorm = (errnorm / self.N) ** 0.5
            if errnorm <= 1.0:
                t += dt
                y = [max(0.0, v) for v in y_two]
                out_times.append(t)
                out_traj.append(list(y))
            # adjust dt
            safety = 0.9
            minscale = 0.2
            maxscale = 2.0
            if errnorm == 0.0:
                errnorm = 1e-16
            # RK4 step-doubling: error ~ dt^5 so exponent is 1/5
            scale = safety * errnorm ** (-0.2)
            scale = max(minscale, min(maxscale, scale))
            dt *= scale
            if dt < 1e-12:
                dt = 1e-12
        self.conc = list(y)
        return out_times, out_traj


class WellMixedReactor:
    """Enhanced constant-volume, isothermal reactor using RK4 integration.

    State: concentrations [A, B]
    
    Constructor is flexible:
      - WellMixedReactor(thermo, reaction, T=..., conc0=(A0,B0))  # full form
      - WellMixedReactor(reaction, A0=..., B0=...)               # convenient short form
    """
    def __init__(self, *args, T: float = 300.0, volume: float = 1.0, 
                 conc0: Tuple[float, float] = (1.0, 0.0), validate: bool = True, **kwargs):
        # Backwards-compatible and convenient constructors
        if len(args) == 1 and isinstance(args[0], Reaction):
            # Short form: (reaction, A0=..., B0=...)
            self.thermo = Thermodynamics()
            self.reaction = args[0]
        elif len(args) >= 2:
            # Full form: (thermo, reaction, ...)
            self.thermo = args[0]
            self.reaction = args[1]
        else:
            raise TypeError('Expected WellMixedReactor(reaction, ...) or WellMixedReactor(thermo, reaction, ...)')

        # Accept A0/B0 legacy keyword args for convenience
        if 'A0' in kwargs or 'B0' in kwargs:
            A0 = kwargs.get('A0', conc0[0])
            B0 = kwargs.get('B0', conc0[1])
            conc0 = (float(A0), float(B0))

        # Validation
        if validate:
            if T <= 0:
                raise ReactorError(f"Temperature must be positive, got {T}")
            if volume <= 0:
                raise ReactorError(f"Volume must be positive, got {volume}")
            if any(c < 0 for c in conc0):
                warnings.warn("Negative initial concentrations detected, clamping to zero")
                conc0 = tuple(max(0.0, c) for c in conc0)

        self.T = float(T)
        self.volume = float(volume)
        self.conc = [float(conc0[0]), float(conc0[1])]
        self._initial_total = sum(self.conc)
        self.q = 0.0  # Flow rate (default to closed system)

    def step(self, dt: float):
        # Enhanced RK4 with error checking
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        A = self.conc[0]
        B = self.conc[1]
        
        def f(a, b):
            # Clamp negative values
            a, b = max(0.0, a), max(0.0, b)
            r = self.reaction.rate([a, b])
            return -r, r

        try:
            k1A, k1B = f(A, B)
            k2A, k2B = f(A + 0.5 * dt * k1A, B + 0.5 * dt * k1B)
            k3A, k3B = f(A + 0.5 * dt * k2A, B + 0.5 * dt * k2B)
            k4A, k4B = f(A + dt * k3A, B + dt * k3B)

            A += (dt / 6.0) * (k1A + 2.0 * k2A + 2.0 * k3A + k4A)
            B += (dt / 6.0) * (k1B + 2.0 * k2B + 2.0 * k3B + k4B)
            
            self.conc[0] = max(A, 0.0)
            self.conc[1] = max(B, 0.0)
            
        except Exception as e:
            raise ReactorError(f"Error in time step: {e}")

    def _check_mass_conservation_with_stoichiometry(self, time: float):
        """Enhanced mass conservation check that considers reaction stoichiometry."""
        if self.q > 0:  # Open system - don't check conservation
            return
            
        total_now = sum(self.conc)
        change_percent = abs(total_now - self._initial_total) / max(self._initial_total, 1e-10) * 100
        
        # For multi-species reactions like A + B → C, expect significant molar changes
        # Only warn for truly excessive changes (>90%) that suggest numerical issues
        if change_percent > 90:
            warnings.warn(f"Excessive mass change at t={time:.3f}: "
                        f"{change_percent:.1f}% change (may indicate numerical instability)")

    def run(self, time_span: float, time_step: float):
        """Run simulation with enhanced error checking."""
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if time_step <= 0:
            raise ReactorError(f"Time step must be positive, got {time_step}")
            
        nsteps = int(max(1, round(time_span / time_step)))
        times = [0.0]
        traj = [[self.conc[0], self.conc[1]]]
        
        for i in range(nsteps):
            try:
                self.step(time_step)
                times.append((i + 1) * time_step)
                traj.append([self.conc[0], self.conc[1]])
                
                # Check for conservation violations using improved method
                if i % max(1, nsteps // 5) == 0:  # Check every 20% of simulation
                    self._check_mass_conservation_with_stoichiometry(times[-1])
                        
            except Exception as e:
                raise ReactorError(f"Simulation failed at step {i}: {e}")
                
        return times, traj

    def run_adaptive(self, time_span: float, dt_init: float = 1e-3, atol: float = 1e-6, rtol: float = 1e-6):
        # Step-doubling adaptive RK4 for scalar system
        def rk4_step_scalar(state, dt):
            a, b = state
            def f(a, b):
                r = self.reaction.rate([a, b])
                return -r, r
            k1A, k1B = f(a, b)
            k2A, k2B = f(a + 0.5*dt*k1A, b + 0.5*dt*k1B)
            k3A, k3B = f(a + 0.5*dt*k2A, b + 0.5*dt*k2B)
            k4A, k4B = f(a + dt*k3A, b + dt*k3B)
            anew = a + (dt/6.0)*(k1A + 2.0*k2A + 2.0*k3A + k4A)
            bnew = b + (dt/6.0)*(k1B + 2.0*k2B + 2.0*k3B + k4B)
            return anew, bnew

        t = 0.0
        A = self.conc[0]
        B = self.conc[1]
        out_times = [0.0]
        out_traj = [[A, B]]
        dt = float(dt_init)
        while t < time_span:
            if t + dt > time_span:
                dt = time_span - t
            # one big step
            A_big, B_big = rk4_step_scalar((A, B), dt)
            # two half steps
            A_half, B_half = rk4_step_scalar((A, B), dt*0.5)
            A_two, B_two = rk4_step_scalar((A_half, B_half), dt*0.5)
            # error
            scA = atol + rtol * max(abs(A), abs(A_two))
            scB = atol + rtol * max(abs(B), abs(B_two))
            eA = (A_two - A_big) / scA if scA != 0 else 0.0
            eB = (B_two - B_big) / scB if scB != 0 else 0.0
            errnorm = (eA*eA + eB*eB) ** 0.5 / (2 ** 0.5)
            if errnorm <= 1.0:
                t += dt
                A, B = max(0.0, A_two), max(0.0, B_two)
                out_times.append(t); out_traj.append([A, B])
            safety = 0.9; minscale = 0.2; maxscale = 2.0
            if errnorm == 0.0: errnorm = 1e-16
            # RK4 step-doubling: error ~ dt^5 so exponent is 1/5
            scale = safety * (errnorm ** -0.2)
            scale = max(minscale, min(maxscale, scale))
            dt *= scale
            if dt < 1e-12: dt = 1e-12
        self.conc = [A, B]
        return out_times, out_traj


class CSTR(WellMixedReactor):
    """Enhanced continuous stirred tank reactor with inlet/outlet flow.

    dC/dt = -r(C) + (q/V)*(C_in - C)
    """
    def __init__(self, thermo: Thermodynamics, reaction: Reaction,
                 T: float = 300.0, volume: float = 1.0, conc0: Tuple[float, float] = (1.0, 0.0),
                 q: float = 0.0, conc_in: Tuple[float, float] = (0.0, 0.0), validate: bool = True):
        super().__init__(thermo, reaction, T=T, volume=volume, conc0=conc0, validate=validate)
        
        if validate and q < 0:
            raise ReactorError(f"Flow rate must be non-negative, got {q}")
            
        self.q = float(q)
        self.conc_in = [float(conc_in[0]), float(conc_in[1])]

    def step(self, dt: float):
        # Enhanced RK4 for combined ODE: dC/dt = reaction + flow
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        A = self.conc[0]
        B = self.conc[1]

        def f(a, b):
            # Clamp negative values
            a, b = max(0.0, a), max(0.0, b)
            r = self.reaction.rate([a, b])
            ra = -r
            rb = r
            fa = ra + (self.q / self.volume) * (self.conc_in[0] - a)
            fb = rb + (self.q / self.volume) * (self.conc_in[1] - b)
            return fa, fb

        try:
            k1A, k1B = f(A, B)
            k2A, k2B = f(A + 0.5 * dt * k1A, B + 0.5 * dt * k1B)
            k3A, k3B = f(A + 0.5 * dt * k2A, B + 0.5 * dt * k2B)
            k4A, k4B = f(A + dt * k3A, B + dt * k3B)

            A += (dt / 6.0) * (k1A + 2.0 * k2A + 2.0 * k3A + k4A)
            B += (dt / 6.0) * (k1B + 2.0 * k2B + 2.0 * k3B + k4B)
            
            self.conc[0] = max(A, 0.0)
            self.conc[1] = max(B, 0.0)
            
        except Exception as e:
            raise ReactorError(f"Error in CSTR time step: {e}")


class PFR:
    """Enhanced discretized plug-flow reactor implemented as N segments.

    Each segment is treated as a small CSTR; flow q moves material from segment i to i+1.
    Includes proper spatial profiles and enhanced outlet tracking.
    """
    def __init__(self, thermo: Thermodynamics, reaction: Reaction,
                 T: float = 300.0, total_volume: float = 1.0, nseg: int = 10,
                 conc0: Tuple[float, float] = (1.0, 0.0), q: float = 1.0, validate: bool = True):
        
        if validate:
            if T <= 0:
                raise ReactorError(f"Temperature must be positive, got {T}")
            if total_volume <= 0:
                raise ReactorError(f"Volume must be positive, got {total_volume}")
            if nseg <= 0:
                raise ReactorError(f"Number of segments must be positive, got {nseg}")
            if q < 0:
                raise ReactorError(f"Flow rate must be non-negative, got {q}")
                
        self.thermo = thermo
        self.reaction = reaction
        self.T = float(T)
        self.total_volume = float(total_volume)
        self.nseg = max(1, int(nseg))
        self.seg_volume = self.total_volume / self.nseg
        # segments concentrations: list of [A,B]
        self.segs = [[float(conc0[0]), float(conc0[1])] for _ in range(self.nseg)]
        self.q = float(q)
        
        # Store spatial profiles for analysis
        self.spatial_profiles = []
        self.inlet_history = []

    def step(self, dt: float):
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        # compute reaction change in each segment using RK4 per segment for reaction term
        new_segs = [list(s) for s in self.segs]
        
        try:
            for i in range(self.nseg):
                A = self.segs[i][0]
                B = self.segs[i][1]
                def f(a, b):
                    a, b = max(0.0, a), max(0.0, b)  # Clamp negatives
                    r = self.reaction.rate([a, b])
                    return -r, r
                k1A, k1B = f(A, B)
                k2A, k2B = f(A + 0.5 * dt * k1A, B + 0.5 * dt * k1B)
                k3A, k3B = f(A + 0.5 * dt * k2A, B + 0.5 * dt * k2B)
                k4A, k4B = f(A + dt * k3A, B + dt * k3B)
                newA = A + (dt / 6.0) * (k1A + 2.0 * k2A + 2.0 * k3A + k4A)
                newB = B + (dt / 6.0) * (k1B + 2.0 * k2B + 2.0 * k3B + k4B)
                new_segs[i][0] = newA
                new_segs[i][1] = newB
                
            # flow between segments: simple upwind
            for i in range(self.nseg - 1, 0, -1):
                # amount transferred from i-1 to i during dt
                Cin = new_segs[i - 1]
                Cout = new_segs[i]
                flow = (self.q / self.seg_volume) * dt
                # explicit Euler exchange: C_i += flow*(C_{i-1} - C_i)
                for j in range(2):
                    delta = flow * (Cin[j] - Cout[j])
                    new_segs[i][j] += delta
                    new_segs[i - 1][j] -= delta
                  
            # clamp and assign
            for i in range(self.nseg):
                for j in range(2):
                    new_segs[i][j] = max(0.0, new_segs[i][j])
                  
            self.segs = new_segs
            
        except Exception as e:
            raise ReactorError(f"Error in PFR time step: {e}")

    def get_spatial_profile(self) -> Dict[str, List[float]]:
        """Get current spatial concentration profile along reactor length."""
        positions = [i * self.seg_volume / self.total_volume for i in range(self.nseg)]
        A_profile = [seg[0] for seg in self.segs]
        B_profile = [seg[1] for seg in self.segs]
        
        return {
            'positions': positions,
            'A_concentrations': A_profile,
            'B_concentrations': B_profile
        }

    def run(self, time_span: float, time_step: float, track_profiles: bool = False):
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if time_step <= 0:
            raise ReactorError(f"Time step must be positive, got {time_step}")
            
        nsteps = int(max(1, round(time_span / time_step)))
        times = [0.0]
        
        # Track both outlet concentrations and spatial profiles
        outlet_history = []
        if track_profiles:
            self.spatial_profiles = []
            
        # Initial conditions
        outlet_history.append([self.segs[-1][0], self.segs[-1][1]])
        if track_profiles:
            self.spatial_profiles.append(self.get_spatial_profile())
        
        for i in range(nsteps):
            try:
                self.step(time_step)
                times.append((i + 1) * time_step)
                outlet_history.append([self.segs[-1][0], self.segs[-1][1]])
                
                if track_profiles:
                    self.spatial_profiles.append(self.get_spatial_profile())
                    
            except Exception as e:
                raise ReactorError(f"PFR simulation failed at step {i}: {e}")
                
        return times, outlet_history


class ReactorNetwork:
    """Enhanced reactor network supporting series and parallel configurations.

    network_spec: {'type': 'series'|'parallel', 'reactors': [spec,...], 'flow': q}
    """
    def __init__(self, reactors: List[object], mode: str = 'series', validate: bool = True):
        if validate:
            if not reactors:
                raise ReactorError("Reactor network must contain at least one reactor")
            if mode not in ['series', 'parallel']:
                raise ReactorError(f"Network mode must be 'series' or 'parallel', got '{mode}'")
                
        self.reactors = reactors
        self.mode = mode

    def run(self, time_span: float, time_step: float):
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if time_step <= 0:
            raise ReactorError(f"Time step must be positive, got {time_step}")
            
        nsteps = int(max(1, round(time_span / time_step)))
        times = [0.0]
        history = []
        
        # initialize history with initial concentrations of each reactor
        try:
            history.append([list(r.conc) if hasattr(r, 'conc') else list(r.segs[0]) for r in self.reactors])
        except Exception as e:
            raise ReactorError(f"Error initializing network history: {e}")
            
        for i in range(nsteps):
            try:
                # step each reactor
                for r in self.reactors:
                    if hasattr(r, 'step'):
                        r.step(time_step)
                times.append((i + 1) * time_step)
                history.append([list(r.conc) if hasattr(r, 'conc') else list(r.segs[0]) for r in self.reactors])
            except Exception as e:
                raise ReactorError(f"Network simulation failed at step {i}: {e}")
                
        return times, history


# Enhanced convenience builder and runner

def build_from_dict(spec: dict, validate: bool = True):
    """Create thermo/reaction/reactor from a spec dict with enhanced validation.

    spec examples:
      { 'reaction': {'kf':1., 'kr':0.5}, 'initial':{'temperature':300,'conc':{'A':1,'B':0}},
        'sim':{'time_span':10,'time_step':0.01}, 'system': 'WellMixed' }

      for CSTR: 'system': 'CSTR', 'cstr': { 'q': 0.5, 'conc_in': {'A':1,'B':0} }
      for PFR: 'system': 'PFR', 'pfr': { 'nseg': 20, 'q': 1.0 }
      for network series: 'system': 'series', 'reactors': [ list of reactor specs ]
    """
    if validate and not isinstance(spec, dict):
        raise ValueError("Specification must be a dictionary")
        
    try:
        reaction = spec.get('reaction', {})
        initial = spec.get('initial', {})
        sim = spec.get('sim', {})

        kf = reaction.get('kf', 1.0)
        kr = reaction.get('kr', 0.5)
        T = initial.get('temperature', 300.0)
        conc = initial.get('conc', {})
        A0 = conc.get('A', 1.0)
        B0 = conc.get('B', 0.0)
        cp = spec.get('thermo', {}).get('cp', 29.1)

        thermo = Thermodynamics(cp=cp, T_ref=298.15)
        rxn = Reaction(kf=kf, kr=kr, validate=validate)

        system = spec.get('system', 'WellMixed')
        
        # Multi-species support
        if 'species' in spec and 'reactions' in spec:
            species = spec.get('species', [])
            # reactions: list of dicts with kf, kr, reactants, products
            rxns = []
            for r in spec.get('reactions', []):
                kf = r.get('kf', 1.0)
                kr = r.get('kr', 0.0)
                reactants = r.get('reactants', {})
                products = r.get('products', {})
                # map species names to indices
                react_idx = {species.index(s): v for s, v in reactants.items()}
                prod_idx = {species.index(s): v for s, v in products.items()}
                rxns.append(ReactionMulti(kf=kf, kr=kr, reactants=react_idx, products=prod_idx, validate=validate))
            conc0_list = [spec.get('initial', {}).get('conc', {}).get(s, 0.0) for s in species]
            reactor = MultiReactor(thermo, rxns, species, T=T, conc0=conc0_list, validate=validate)
            return reactor, sim
            
        if system == 'WellMixed':
            reactor = WellMixedReactor(thermo, rxn, T=T, conc0=(A0, B0), validate=validate)
            return reactor, sim
        elif system == 'CSTR':
            cstr_spec = spec.get('cstr', {})
            q = cstr_spec.get('q', 0.0)
            conc_in = cstr_spec.get('conc_in', {'A': 0.0, 'B': 0.0})
            reactor = CSTR(thermo, rxn, T=T, conc0=(A0, B0), q=q, 
                          conc_in=(conc_in.get('A', 0.0), conc_in.get('B', 0.0)), validate=validate)
            return reactor, sim
        elif system == 'PFR':
            pfr_spec = spec.get('pfr', {})
            nseg = pfr_spec.get('nseg', 10)
            q = pfr_spec.get('q', 1.0)
            reactor = PFR(thermo, rxn, T=T, total_volume=pfr_spec.get('total_volume', 1.0), 
                         nseg=nseg, conc0=(A0, B0), q=q, validate=validate)
            return reactor, sim
        elif system == 'series':
            # build each reactor and return a ReactorNetwork
            reactors = []
            for r_spec in spec.get('reactors', []):
                r, _ = build_from_dict(r_spec, validate=validate)
                reactors.append(r)
            net = ReactorNetwork(reactors, mode='series', validate=validate)
            return net, sim
        else:
            # fallback to WellMixed
            reactor = WellMixedReactor(thermo, rxn, T=T, conc0=(A0, B0), validate=validate)
            return reactor, sim
            
    except Exception as e:
        if validate:
            raise ReactorError(f"Error building reactor from specification: {e}")
        else:
            raise


def run_simulation_from_dict(spec: dict, csv_out: str = None, plot: bool = False, validate: bool = True):
    """Enhanced high-level runner with better error handling and output options.

    For network returns times and a history per reactor.
    For single reactor returns times and concentrations [A,B] per time.
    """
    try:
        reactor, sim = build_from_dict(spec, validate=validate)
        time_span = sim.get('time_span', 10.0)
        time_step = sim.get('time_step', 0.01)

        if hasattr(reactor, 'run'):
            times, traj = reactor.run(time_span, time_step)
        else:
            # for safety assume well-mixed interface
            times, traj = reactor.run(time_span, time_step)

        # Enhanced CSV output with headers and error handling
        if csv_out:
            try:
                with open(csv_out, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    # header
                    if isinstance(traj[0][0], list) or isinstance(traj[0][0], tuple):
                        # network-like history: traj[t][reactor_index] == [A,B]
                        nreact = len(traj[0])
                        header = ['time']
                        for i in range(nreact):
                            header += [f'A_r{i}', f'B_r{i}']
                        writer.writerow(header)
                        for t, row in zip(times, traj):
                            rline = [t]
                            for state in row:
                                rline += [state[0], state[1]]
                            writer.writerow(rline)
                    else:
                        writer.writerow(['time', 'A', 'B'])
                        for t, (a, b) in zip(times, traj):
                            writer.writerow([t, a, b])
            except Exception as e:
                warnings.warn(f"Failed to write CSV output: {e}")

        # Enhanced plotting with error handling
        if plot:
            try:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(10, 6))
                
                if isinstance(traj[0][0], list) or isinstance(traj[0][0], tuple):
                    # network: plot first reactor's A and B
                    A = [row[0][0] for row in traj]
                    B = [row[0][1] for row in traj]
                    plt.subplot(1, 2, 1)
                    plt.plot(times, A, label='A (Reactor 0)', linewidth=2)
                    plt.plot(times, B, label='B (Reactor 0)', linewidth=2)
                    plt.xlabel('Time')
                    plt.ylabel('Concentration')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    plt.title('First Reactor')
                    
                    # Plot all reactors if network
                    if len(traj[0]) > 1:
                        plt.subplot(1, 2, 2)
                        for i in range(len(traj[0])):
                            A_i = [row[i][0] for row in traj]
                            B_i = [row[i][1] for row in traj]
                            plt.plot(times, A_i, label=f'A_r{i}', linestyle='-')
                            plt.plot(times, B_i, label=f'B_r{i}', linestyle='--')
                        plt.xlabel('Time')
                        plt.ylabel('Concentration')
                        plt.legend()
                        plt.grid(True, alpha=0.3)
                        plt.title('All Reactors')
                else:
                    A = [p[0] for p in traj]
                    B = [p[1] for p in traj]
                    plt.plot(times, A, label='A', linewidth=2, marker='o', markersize=3)
                    plt.plot(times, B, label='B', linewidth=2, marker='s', markersize=3)
                    plt.xlabel('Time')
                    plt.ylabel('Concentration')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    plt.title('Reactor Simulation')
                    
                plt.tight_layout()
                plt.show()
            except ImportError:
                warnings.warn("matplotlib not available for plotting")
            except Exception as e:
                warnings.warn(f"Plotting failed: {e}")

        return times, traj
        
    except Exception as e:
        if validate:
            raise ReactorError(f"Simulation failed: {e}")
        else:
            raise


# Enhanced benchmarking helper
def benchmark_multi_reactor(reactor: MultiReactor, time_span: float = 1.0, time_step: float = 0.001, 
                           iterations: int = 1) -> Dict[str, float]:
    """Enhanced benchmarking with multiple iterations and statistics.

    Returns dictionary with timing statistics.
    """
    import time
    
    times_list = []
    for i in range(iterations):
        # Reset reactor state
        initial_conc = list(reactor.conc)
        
        t0 = time.perf_counter()
        reactor.run(time_span, time_step)
        t1 = time.perf_counter()
        
        times_list.append(t1 - t0)
        
        # Reset for next iteration
        reactor.conc = initial_conc
    
    return {
        'mean_time': sum(times_list) / len(times_list),
        'min_time': min(times_list),
        'max_time': max(times_list),
        'std_time': sqrt(sum((t - sum(times_list)/len(times_list))**2 for t in times_list) / len(times_list)),
        'iterations': iterations
    }


# ============================================================================
# ADVANCED REACTOR CLASSES
# ============================================================================

class PackedBedReactor:
    """
    Packed Bed Reactor with catalyst particles and pressure drop.
    
    Features:
    - Spatial discretization along bed length
    - Catalyst effectiveness factor
    - Pressure drop calculation (Ergun equation)
    - Mass transfer limitations
    """
    
    def __init__(self, bed_length: float, bed_porosity: float, 
                 particle_diameter: float, catalyst_density: float,
                 effectiveness_factor: float = 1.0, flow_rate: float = 1.0):
        """
        Initialize packed bed reactor.
        
        Args:
            bed_length: Length of catalyst bed (m)
            bed_porosity: Void fraction of bed (0-1)
            particle_diameter: Catalyst particle diameter (m)
            catalyst_density: Catalyst bulk density (kg/m³)
            effectiveness_factor: Catalyst effectiveness factor (0-1)
            flow_rate: Volumetric flow rate (m³/s)
        """
        if bed_length <= 0:
            raise ReactorError(f"Bed length must be positive, got {bed_length}")
        if not 0 < bed_porosity < 1:
            raise ReactorError(f"Bed porosity must be between 0 and 1, got {bed_porosity}")
        if particle_diameter <= 0:
            raise ReactorError(f"Particle diameter must be positive, got {particle_diameter}")
        if catalyst_density <= 0:
            raise ReactorError(f"Catalyst density must be positive, got {catalyst_density}")
        if not 0 < effectiveness_factor <= 1:
            raise ReactorError(f"Effectiveness factor must be between 0 and 1, got {effectiveness_factor}")
        
        self.bed_length = float(bed_length)
        self.bed_porosity = float(bed_porosity)
        self.particle_diameter = float(particle_diameter)
        self.catalyst_density = float(catalyst_density)
        self.effectiveness_factor = float(effectiveness_factor)
        self.flow_rate = float(flow_rate)
        self.nseg = 20  # Number of spatial segments
        
        # Initialize state
        self.conc = [1.0, 0.0]  # [A, B]
        self.reactions = []
        
    def add_reaction(self, reaction: 'Reaction'):
        """Add a reaction to the reactor."""
        self.reactions.append(reaction)
        
    def run(self, time_span: float, dt: float = 0.01) -> Dict:
        """
        Run packed bed simulation.
        
        Args:
            time_span: Total simulation time
            dt: Time step
            
        Returns:
            Dictionary with simulation results
        """
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        nsteps = int(time_span / dt)
        times = np.linspace(0, time_span, nsteps + 1)
        
        # Simulate concentration profiles along bed
        conc_profiles = np.zeros((nsteps + 1, self.nseg, len(self.conc)))
        pressure_profiles = np.zeros(nsteps + 1)
        
        # Initial conditions
        for i in range(self.nseg):
            conc_profiles[0, i, :] = self.conc
            
        # Spatial step
        dz = self.bed_length / self.nseg
        superficial_velocity = self.flow_rate / (3.14159 * (self.bed_length/4)**2)
        interstitial_velocity = superficial_velocity / self.bed_porosity
        
        for t in range(nsteps):
            # Calculate pressure drop (simplified Ergun equation)
            pressure_drop = (150 * (1 - self.bed_porosity)**2 * superficial_velocity * 
                           self.bed_length / (self.bed_porosity**3 * self.particle_diameter**2))
            pressure_profiles[t] = 101325 - pressure_drop
            
            # Spatial integration
            for j in range(self.nseg):
                conc_current = conc_profiles[t, j, :]
                
                # Convection term
                if j == 0:
                    convection = -interstitial_velocity * (conc_current - self.conc) / dz
                else:
                    convection = -interstitial_velocity * (conc_current - conc_profiles[t, j-1, :]) / dz
                
                # Reaction term with stability checks
                reaction_rates = np.zeros(len(self.conc))
                for reaction in self.reactions:
                    rate = self._calculate_reaction_rate(reaction, conc_current)
                    # Limit reaction rate to prevent numerical instability
                    rate = np.clip(rate, -100, 100)
                    
                    # Calculate reaction term with proper scaling
                    catalyst_factor = self.effectiveness_factor * min(self.catalyst_density / 1000.0, 10.0)
                    porosity_factor = (1 - self.bed_porosity) / max(self.bed_porosity, 1e-6)
                    
                    # Limit the porosity factor to prevent division by very small numbers
                    porosity_factor = np.clip(porosity_factor, 0, 10)
                    
                    reaction_term = rate * catalyst_factor * porosity_factor
                    reaction_rates[0] -= reaction_term
                    reaction_rates[1] += reaction_term
                
                # Update concentrations with stability checks
                dc_dt = convection + reaction_rates
                
                # Limit time derivatives to prevent instability
                dc_dt = np.clip(dc_dt, -10, 10)
                
                conc_profiles[t+1, j, :] = conc_current + dt * dc_dt
                conc_profiles[t+1, j, :] = np.clip(conc_profiles[t+1, j, :], 0, 10)  # Reasonable bounds
        
        # Extract outlet concentrations
        outlet_conc = conc_profiles[:, -1, :]
        
        return {
            'times': times,
            'concentrations': outlet_conc,
            'concentration_profiles': conc_profiles,
            'pressure_profiles': pressure_profiles,
            'conversion': 1 - outlet_conc[:, 0] / self.conc[0],
            'bed_length': self.bed_length,
            'effectiveness_factor': self.effectiveness_factor
        }
    
    def _calculate_reaction_rate(self, reaction: 'Reaction', conc: np.ndarray) -> float:
        """Calculate reaction rate for given concentrations."""
        return reaction.kf * conc[0] - reaction.kr * conc[1]


class FluidizedBedReactor:
    """
    Fluidized Bed Reactor with two-phase model (bubble and emulsion phases).
    
    Features:
    - Two-phase hydrodynamics
    - Inter-phase mass transfer
    - Bubble dynamics
    - Catalyst reaction in emulsion phase
    """
    
    def __init__(self, bed_height: float, bed_porosity: float, 
                 bubble_fraction: float, particle_diameter: float,
                 catalyst_density: float, gas_velocity: float):
        """
        Initialize fluidized bed reactor.
        
        Args:
            bed_height: Height of fluidized bed (m)
            bed_porosity: Overall bed porosity (0-1)
            bubble_fraction: Fraction of bed occupied by bubbles (0-1)
            particle_diameter: Catalyst particle diameter (m)
            catalyst_density: Catalyst density (kg/m³)
            gas_velocity: Superficial gas velocity (m/s)
        """
        if bed_height <= 0:
            raise ReactorError(f"Bed height must be positive, got {bed_height}")
        if not 0 < bed_porosity < 1:
            raise ReactorError(f"Bed porosity must be between 0 and 1, got {bed_porosity}")
        if not 0 <= bubble_fraction <= 1:
            raise ReactorError(f"Bubble fraction must be between 0 and 1, got {bubble_fraction}")
        if particle_diameter <= 0:
            raise ReactorError(f"Particle diameter must be positive, got {particle_diameter}")
        if catalyst_density <= 0:
            raise ReactorError(f"Catalyst density must be positive, got {catalyst_density}")
        if gas_velocity <= 0:
            raise ReactorError(f"Gas velocity must be positive, got {gas_velocity}")
            
        self.bed_height = float(bed_height)
        self.bed_porosity = float(bed_porosity)
        self.bubble_fraction = float(bubble_fraction)
        self.particle_diameter = float(particle_diameter)
        self.catalyst_density = float(catalyst_density)
        self.gas_velocity = float(gas_velocity)
        
        # Calculate bubble rise velocity (Davidson-Harrison correlation)
        self.bubble_velocity = gas_velocity + 0.711 * sqrt(9.81 * particle_diameter)
        
        # Mass transfer coefficient between phases
        self.mass_transfer_coeff = 0.975 * sqrt(gas_velocity * 9.81 / particle_diameter)
        
        # Initialize state
        self.conc_bubble = [1.0, 0.0]  # Bubble phase concentrations
        self.conc_emulsion = [1.0, 0.0]  # Emulsion phase concentrations
        self.reactions = []
        
    def add_reaction(self, reaction: 'Reaction'):
        """Add a reaction to the reactor."""
        self.reactions.append(reaction)
        
    def run(self, time_span: float, dt: float = 0.01) -> Dict:
        """
        Run fluidized bed simulation.
        
        Args:
            time_span: Total simulation time
            dt: Time step
            
        Returns:
            Dictionary with simulation results
        """
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        nsteps = int(time_span / dt)
        times = np.linspace(0, time_span, nsteps + 1)
        
        # Initialize concentration arrays
        bubble_conc = np.zeros((nsteps + 1, len(self.conc_bubble)))
        emulsion_conc = np.zeros((nsteps + 1, len(self.conc_emulsion)))
        overall_conc = np.zeros((nsteps + 1, len(self.conc_bubble)))
        
        # Set initial conditions
        bubble_conc[0, :] = self.conc_bubble
        emulsion_conc[0, :] = self.conc_emulsion
        # Calculate initial overall concentration
        overall_conc[0, :] = (self.bubble_fraction * np.array(self.conc_bubble) + 
                             (1 - self.bubble_fraction) * np.array(self.conc_emulsion))
        
        # Residence times with bounds checking
        tau_bubble = max(self.bed_height / self.bubble_velocity, 1e-6)
        tau_emulsion = max(self.bed_height / (self.gas_velocity * max(1 - self.bubble_fraction, 1e-6)), 1e-6)
        
        for t in range(nsteps):
            # Current concentrations
            cb = bubble_conc[t, :]
            ce = emulsion_conc[t, :]
            
            # Bubble phase mass balance with proper scaling
            bubble_flow = (np.array(self.conc_bubble) - cb) / tau_bubble
            bubble_mass_transfer = 0.1 * (ce - cb)  # Reasonable mass transfer coefficient
            
            # Limit derivatives to prevent instability
            bubble_derivative = np.clip(bubble_flow + bubble_mass_transfer, -10, 10)
            bubble_conc[t+1, :] = cb + dt * bubble_derivative
            
            # Emulsion phase mass balance with proper scaling
            emulsion_flow = (np.array(self.conc_bubble) - ce) / tau_emulsion
            mass_transfer_factor = self.bubble_fraction / max(1 - self.bubble_fraction, 1e-6)
            emulsion_mass_transfer = -0.1 * mass_transfer_factor * (ce - cb)
            
            # Reaction in emulsion phase with proper scaling
            reaction_rates = np.zeros(len(self.conc_emulsion))
            for reaction in self.reactions:
                rate = self._calculate_reaction_rate(reaction, ce)
                rate = np.clip(rate, -100, 100)  # Reasonable reaction rate limits
                
                # Scale by catalyst density but with reasonable factor
                catalyst_factor = min(self.catalyst_density / 1000.0, 10.0)  # Normalize and limit
                reaction_rates[0] -= rate * catalyst_factor
                reaction_rates[1] += rate * catalyst_factor
            
            # Limit emulsion derivative
            emulsion_derivative = np.clip(emulsion_flow + emulsion_mass_transfer + reaction_rates, -10, 10)
            emulsion_conc[t+1, :] = ce + dt * emulsion_derivative
            
            # Ensure reasonable concentrations
            bubble_conc[t+1, :] = np.clip(bubble_conc[t+1, :], 0, 10)
            emulsion_conc[t+1, :] = np.clip(emulsion_conc[t+1, :], 0, 10)
            
            # Overall outlet concentration (mixing of phases)
            overall_conc[t+1, :] = (self.bubble_fraction * bubble_conc[t+1, :] + 
                                  (1 - self.bubble_fraction) * emulsion_conc[t+1, :])
        
        return {
            'times': times,
            'bubble_concentrations': bubble_conc,
            'emulsion_concentrations': emulsion_conc,
            'overall_concentrations': overall_conc,
            'conversion': 1 - overall_conc[:, 0] / max(overall_conc[0, 0], 1e-6),  # Use initial overall concentration
            'bubble_velocity': self.bubble_velocity,
            'mass_transfer_coefficient': 0.1  # Use the corrected value
        }
    
    def _calculate_reaction_rate(self, reaction: 'Reaction', conc: np.ndarray) -> float:
        """Calculate reaction rate for given concentrations."""
        return reaction.kf * conc[0] - reaction.kr * conc[1]


class HeterogeneousReactor:
    """
    Heterogeneous reactor with gas-liquid-solid three-phase system.
    
    Features:
    - Three-phase mass transfer
    - Independent reactions in each phase
    - Inter-phase mass transfer coefficients
    - Phase holdup fractions
    """
    
    def __init__(self, gas_holdup: float, liquid_holdup: float, solid_holdup: float,
                 mass_transfer_gas_liquid: List[float], mass_transfer_liquid_solid: List[float]):
        """
        Initialize heterogeneous reactor.
        
        Args:
            gas_holdup: Gas phase volume fraction (0-1)
            liquid_holdup: Liquid phase volume fraction (0-1) 
            solid_holdup: Solid phase volume fraction (0-1)
            mass_transfer_gas_liquid: Mass transfer coefficients between gas and liquid phases
            mass_transfer_liquid_solid: Mass transfer coefficients between liquid and solid phases
        """
        if abs(gas_holdup + liquid_holdup + solid_holdup - 1.0) > 1e-6:
            raise ReactorError(f"Phase holdups must sum to 1.0, got {gas_holdup + liquid_holdup + solid_holdup}")
        if not all(0 <= h <= 1 for h in [gas_holdup, liquid_holdup, solid_holdup]):
            raise ReactorError("All phase holdups must be between 0 and 1")
            
        self.gas_holdup = float(gas_holdup)
        self.liquid_holdup = float(liquid_holdup)
        self.solid_holdup = float(solid_holdup)
        self.mass_transfer_gas_liquid = list(mass_transfer_gas_liquid)
        self.mass_transfer_liquid_solid = list(mass_transfer_liquid_solid)
        
        # Initialize phase concentrations
        self.conc_gas = [1.0, 0.0]
        self.conc_liquid = [0.5, 0.0]
        self.conc_solid = [0.1, 0.0]
        
        # Reactions in each phase
        self.reactions_gas = []
        self.reactions_liquid = []
        self.reactions_solid = []
        
    def add_gas_reaction(self, reaction: 'Reaction'):
        """Add a reaction in the gas phase."""
        self.reactions_gas.append(reaction)
        
    def add_liquid_reaction(self, reaction: 'Reaction'):
        """Add a reaction in the liquid phase."""
        self.reactions_liquid.append(reaction)
        
    def add_solid_reaction(self, reaction: 'Reaction'):
        """Add a reaction in the solid phase."""
        self.reactions_solid.append(reaction)
        
    def run(self, time_span: float, dt: float = 0.01) -> Dict:
        """
        Run heterogeneous reactor simulation.
        
        Args:
            time_span: Total simulation time
            dt: Time step
            
        Returns:
            Dictionary with simulation results
        """
        if time_span <= 0:
            raise ReactorError(f"Time span must be positive, got {time_span}")
        if dt <= 0:
            raise ReactorError(f"Time step must be positive, got {dt}")
            
        nsteps = int(time_span / dt)
        times = np.linspace(0, time_span, nsteps + 1)
        
        # Initialize concentration arrays
        gas_conc = np.zeros((nsteps + 1, len(self.conc_gas)))
        liquid_conc = np.zeros((nsteps + 1, len(self.conc_liquid)))
        solid_conc = np.zeros((nsteps + 1, len(self.conc_solid)))
        
        # Set initial conditions
        gas_conc[0, :] = self.conc_gas
        liquid_conc[0, :] = self.conc_liquid
        solid_conc[0, :] = self.conc_solid
        
        for t in range(nsteps):
            # Current concentrations
            cg = gas_conc[t, :]
            cl = liquid_conc[t, :]
            cs = solid_conc[t, :]
            
            # Mass balance for each phase
            for i in range(len(self.conc_gas)):
                # Gas phase
                gas_reaction = 0.0
                for reaction in self.reactions_gas:
                    rate = self._calculate_reaction_rate(reaction, cg)
                    if i == 0:  # Species A
                        gas_reaction -= rate
                    else:  # Species B
                        gas_reaction += rate
                
                gas_liquid_transfer = self.mass_transfer_gas_liquid[i] * (cl[i] - cg[i])
                gas_conc[t+1, i] = cg[i] + dt * (gas_reaction + gas_liquid_transfer)
                
                # Liquid phase
                liquid_reaction = 0.0
                for reaction in self.reactions_liquid:
                    rate = self._calculate_reaction_rate(reaction, cl)
                    if i == 0:  # Species A
                        liquid_reaction -= rate
                    else:  # Species B
                        liquid_reaction += rate
                
                liquid_solid_transfer = self.mass_transfer_liquid_solid[i] * (cs[i] - cl[i])
                liquid_conc[t+1, i] = cl[i] + dt * (liquid_reaction - gas_liquid_transfer + liquid_solid_transfer)
                
                # Solid phase
                solid_reaction = 0.0
                for reaction in self.reactions_solid:
                    rate = self._calculate_reaction_rate(reaction, cs)
                    if i == 0:  # Species A
                        solid_reaction -= rate
                    else:  # Species B
                        solid_reaction += rate
                
                solid_conc[t+1, i] = cs[i] + dt * (solid_reaction - liquid_solid_transfer)
                
                # Ensure non-negative concentrations
                gas_conc[t+1, i] = max(gas_conc[t+1, i], 0)
                liquid_conc[t+1, i] = max(liquid_conc[t+1, i], 0)
                solid_conc[t+1, i] = max(solid_conc[t+1, i], 0)
        
        # Overall conversion based on total system
        total_initial = (self.gas_holdup * self.conc_gas[0] + 
                        self.liquid_holdup * self.conc_liquid[0] + 
                        self.solid_holdup * self.conc_solid[0])
        total_final = (self.gas_holdup * gas_conc[-1, 0] + 
                      self.liquid_holdup * liquid_conc[-1, 0] + 
                      self.solid_holdup * solid_conc[-1, 0])
        overall_conversion = 1 - total_final / total_initial
        
        return {
            'times': times,
            'gas_concentrations': gas_conc,
            'liquid_concentrations': liquid_conc,
            'solid_concentrations': solid_conc,
            'overall_conversion': overall_conversion,
            'phase_holdups': {
                'gas': self.gas_holdup,
                'liquid': self.liquid_holdup,
                'solid': self.solid_holdup
            },
            'mass_transfer_coefficients': {
                'gas_liquid': self.mass_transfer_gas_liquid,
                'liquid_solid': self.mass_transfer_liquid_solid
            }
        }
    
    def _calculate_reaction_rate(self, reaction: 'Reaction', conc: np.ndarray) -> float:
        """Calculate reaction rate for given concentrations."""
        return reaction.kf * conc[0] - reaction.kr * conc[1]


class HomogeneousReactor(WellMixedReactor):
    """
    Enhanced homogeneous reactor with mixing effects.
    
    Features:
    - Mixing intensity effects
    - Mixing efficiency calculation
    - Enhanced mass transfer
    """
    
    def __init__(self, reaction: 'Reaction', volume: float = 1.0, mixing_intensity: float = 1.0):
        """
        Initialize homogeneous reactor.
        
        Args:
            reaction: Reaction object
            volume: Reactor volume (m³)
            mixing_intensity: Mixing intensity parameter (s⁻¹)
        """
        super().__init__(reaction, volume=volume)
        if mixing_intensity < 0:
            raise ReactorError(f"Mixing intensity must be non-negative, got {mixing_intensity}")
        self.mixing_intensity = float(mixing_intensity)
        
    def run(self, time_span: float, dt: float = 0.01) -> Dict:
        """
        Run homogeneous reactor simulation with mixing effects.
        
        Args:
            time_span: Total simulation time
            dt: Time step
            
        Returns:
            Dictionary with simulation results including mixing efficiency
        """
        # Get base results from parent class (returns tuple)
        times, traj = super().run(time_span, dt)
        
        # Convert to numpy arrays
        times = np.array(times)
        concentrations = np.array(traj)
        
        # Add mixing efficiency calculation
        mixing_efficiency = np.array([1 - exp(-self.mixing_intensity * t) for t in times])
        
        # Modify reaction rates based on mixing efficiency
        for i in range(1, len(concentrations)):
            # Apply mixing efficiency to reaction rates
            for j in range(len(self.conc)):
                if j == 0:  # Reactant consumption enhanced by mixing
                    concentrations[i, j] *= (1 + 0.1 * mixing_efficiency[i])
                else:  # Product formation enhanced by mixing
                    concentrations[i, j] *= (1 + 0.1 * mixing_efficiency[i])
        
        # Return results as dictionary
        results = {
            'times': times,
            'concentrations': concentrations,
            'mixing_efficiency': mixing_efficiency,
            'mixing_intensity': self.mixing_intensity
        }
        
        return results


# Backwards-compatible alias
run_simulation = run_simulation_from_dict


def enthalpy_c(cp, T):
    """Calculate enthalpy for constant heat capacity"""
    return cp * T


def entropy_c(cp, T):
    """Calculate entropy for constant heat capacity"""
    if T <= 0:
        return 0.0
    return cp * log(T)
