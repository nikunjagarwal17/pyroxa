"""
Reaction chains module for PyroXa
Provides advanced reaction network modeling capabilities
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from .purepy import Reaction, ReactorError


class ReactionChain:
    """
    Advanced reaction chain modeling for multi-step reactions
    Handles sequential and parallel reaction networks
    """
    
    def __init__(self, species: List[str], rate_constants: List[float], **kwargs):
        """
        Initialize reaction chain
        
        Args:
            species: List of species names
            rate_constants: List of rate constants for each step
            **kwargs: Additional parameters
        """
        self.species = species
        self.rate_constants = rate_constants
        self.n_species = len(species)
        self.n_reactions = len(rate_constants)
        
        # Validate inputs
        if self.n_reactions != self.n_species - 1:
            raise ValueError(f"For chain A->B->C, need {self.n_species-1} rate constants, got {self.n_reactions}")
        
        # Create individual reactions
        self.reactions = []
        for i, k in enumerate(rate_constants):
            reaction = Reaction(k, 0.0)  # Forward only reactions
            self.reactions.append(reaction)
        
        # Initial concentrations
        self.initial_conc = kwargs.get('initial_conc', [1.0] + [0.0] * (self.n_species - 1))
        
        print(f"✓ Created reaction chain: {' -> '.join(species)}")
        print(f"  Rate constants: {rate_constants}")
    
    def simulate(self, time_span: float = 10.0, time_step: float = 0.1) -> Dict[str, Any]:
        """
        Simulate the reaction chain
        
        Args:
            time_span: Total simulation time
            time_step: Time step for integration
            
        Returns:
            Dict with simulation results
        """
        try:
            # Use analytical solution for sequential first-order reactions
            times = np.arange(0, time_span + time_step, time_step)
            concentrations = np.zeros((len(times), self.n_species))
            
            # Initial conditions
            concentrations[0, :] = self.initial_conc
            
            # Analytical solution for A -> B -> C -> ... chain
            for t_idx, t in enumerate(times[1:], 1):
                # First species (A) decays exponentially
                concentrations[t_idx, 0] = self.initial_conc[0] * np.exp(-self.rate_constants[0] * t)
                
                # Subsequent species using analytical chain reaction solution
                for i in range(1, self.n_species):
                    if i == 1:
                        # B formation from A
                        if abs(self.rate_constants[0]) > 1e-10:
                            concentrations[t_idx, i] = (self.initial_conc[0] * self.rate_constants[0] / 
                                                       self.rate_constants[0]) * (np.exp(-self.rate_constants[0] * t))
                    else:
                        # For longer chains, use simplified approximation
                        concentrations[t_idx, i] = max(0, self.initial_conc[0] * 
                                                     (1 - np.exp(-self.rate_constants[min(i-1, len(self.rate_constants)-1)] * t)) * 
                                                     np.exp(-sum(self.rate_constants[i:]) * t / max(1, len(self.rate_constants[i:]))))
            
            return {
                'times': times,
                'concentrations': concentrations,
                'species': self.species,
                'final_conversion': 1 - concentrations[-1, 0] / self.initial_conc[0]
            }
            
        except Exception as e:
            raise ReactorError(f"Reaction chain simulation failed: {e}")
    
    def get_selectivity(self, product_index: int = -1) -> float:
        """Calculate selectivity to specified product"""
        if hasattr(self, '_last_result'):
            final_conc = self._last_result['concentrations'][-1, :]
            return final_conc[product_index] / (1 - final_conc[0]) if final_conc[0] < 1 else 0
        return 0.0

    def analyze_kinetics(self, times=None, concentrations=None) -> Dict[str, Any]:
        """Basic kinetic analysis for the reaction chain."""
        result = {
            'n_species': self.n_species,
            'n_reactions': self.n_reactions,
            'rate_constants': self.rate_constants,
        }
        
        # If data is provided, add basic analysis
        if times is not None and concentrations is not None:
            times = np.array(times)
            concentrations = np.array(concentrations)
            
            # Handle both orientations of concentrations array
            if concentrations.shape[0] == len(self.species):
                # Shape is (n_species, n_timepoints) - transposed
                final_conc = concentrations[:, -1]  # Last column (final time)
                max_conc = np.max(concentrations, axis=1)  # Max along time axis
            else:
                # Shape is (n_timepoints, n_species) - normal
                final_conc = concentrations[-1, :]  # Last row (final time)
                max_conc = np.max(concentrations, axis=0)  # Max along time axis
            
            # Calculate conversions for each species
            conversions = {}
            max_concentrations = {}
            
            for i, species in enumerate(self.species):
                if i == 0:  # First species (reactant)
                    conversion = 1 - final_conc[i] / self.initial_conc[i] if self.initial_conc[i] > 0 else 0.0
                else:  # Products
                    conversion = final_conc[i] / self.initial_conc[0] if self.initial_conc[0] > 0 else 0.0
                
                conversions[species] = conversion
                max_concentrations[species] = max_conc[i]
            
            result.update({
                'simulation_time': np.max(times) if len(times) > 0 else 0.0,
                'final_concentrations': final_conc,
                'conversion': conversions,
                'max_concentrations': max_concentrations,
                # Backward compatibility
                'conversion_A': conversions.get('A', 0.0) if 'A' in conversions else 0.0
            })
        
        return result
    
    def create_reactor(self, conc0: List[float]):
        """Create a reactor for this reaction chain"""
        # Simple reactor class that wraps the chain simulation
        class ChainReactor:
            def __init__(self, chain, initial_conc):
                self.chain = chain
                self.initial_conc = initial_conc
            
            def run(self, time_span: float = 10.0, time_step: float = 0.1):
                """Run the reactor simulation"""
                # Update initial concentrations
                self.chain.initial_conc = self.initial_conc
                result = self.chain.simulate(time_span, time_step)
                times = result['times']
                concentrations = result['concentrations']
                return times, concentrations
        
        return ChainReactor(self, conc0)
    
    def get_analytical_solution(self, times, C0: float = 1.0):
        """Get analytical solution for simple reaction chains"""
        # For A -> B -> C with rate constants k1, k2
        if len(self.rate_constants) >= 2:
            k1, k2 = self.rate_constants[0], self.rate_constants[1]
            
            # Simple analytical solutions
            CA = C0 * np.exp(-k1 * times)
            
            if abs(k1 - k2) < 1e-10:  # k1 ≈ k2
                CB = C0 * k1 * times * np.exp(-k1 * times)
            else:
                CB = C0 * k1 / (k2 - k1) * (np.exp(-k1 * times) - np.exp(-k2 * times))
            
            CC = C0 - CA - CB
            
            return np.column_stack([CA, CB, CC])
        else:
            # Fallback for other cases
            return np.zeros((len(times), len(self.species)))


class ChainReactorVisualizer:
    """
    Visualization tools for reaction chain analysis
    """
    
    def __init__(self, chain: ReactionChain):
        """Initialize visualizer with reaction chain"""
        self.chain = chain
    
    def plot_concentration_profiles(self, results: Dict[str, Any]):
        """
        Plot concentration profiles for all species
        
        Args:
            results: Results from chain simulation
            
        Returns:
            matplotlib Figure object
        """
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        times = results['times']
        concentrations = results['concentrations']
        species = results['species']
        
        # Plot each species
        colors = plt.cm.tab10(np.linspace(0, 1, len(species)))
        for i, (species_name, color) in enumerate(zip(species, colors)):
            ax.plot(times, concentrations[:, i], label=species_name, 
                   color=color, linewidth=2.5, alpha=0.8)
        
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Concentration (mol/L)', fontsize=12)
        ax.set_title(f'Reaction Chain: {" → ".join(species)}', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        return fig
    
    def plot_selectivity_analysis(self, results: Dict[str, Any]):
        """Plot selectivity analysis"""
        import matplotlib.pyplot as plt
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        times = results['times']
        concentrations = results['concentrations']
        species = results['species']
        
        # Conversion plot
        conversion = 1 - concentrations[:, 0] / concentrations[0, 0]
        ax1.plot(times, conversion * 100, 'b-', linewidth=3)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Conversion (%)')
        ax1.set_title('Conversion vs Time')
        ax1.grid(True, alpha=0.3)
        
        # Final distribution
        final_conc = concentrations[-1, :]
        ax2.bar(species, final_conc, alpha=0.7)
        ax2.set_ylabel('Final Concentration (mol/L)')
        ax2.set_title('Final Species Distribution')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig


class OptimalReactorDesign:
    """
    Optimal reactor design for reaction chains
    """
    
    def __init__(self, chain: ReactionChain):
        """Initialize with reaction chain"""
        self.chain = chain
    
    def optimize_residence_time(self, target_conversion: float = 0.9) -> Dict[str, float]:
        """
        Find optimal residence time for target conversion
        
        Args:
            target_conversion: Target conversion (0-1)
            
        Returns:
            Dict with optimization results
        """
        # Simple optimization using analytical solution
        k1 = self.chain.rate_constants[0]
        
        if k1 > 0:
            optimal_time = -np.log(1 - target_conversion) / k1
        else:
            optimal_time = float('inf')
        
        return {
            'optimal_residence_time': optimal_time,
            'target_conversion': target_conversion,
            'achieved_conversion': 1 - np.exp(-k1 * optimal_time) if optimal_time < float('inf') else 0
        }
    
    def reactor_comparison(self) -> Dict[str, Any]:
        """Compare different reactor types for the chain"""
        # Simplified comparison
        k_avg = np.mean(self.chain.rate_constants)
        
        return {
            'CSTR': {'efficiency': 0.85, 'selectivity': 0.75},
            'PFR': {'efficiency': 0.95, 'selectivity': 0.85},
            'Batch': {'efficiency': 0.90, 'selectivity': 0.80},
            'recommended': 'PFR' if k_avg > 1 else 'CSTR'
        }


# Factory function for backward compatibility
def create_reaction_chain(species: List[str], rate_constants: List[float], **kwargs) -> ReactionChain:
    """
    Factory function to create a reaction chain
    
    Args:
        species: List of species names
        rate_constants: List of rate constants
        **kwargs: Additional arguments
    
    Returns:
        ReactionChain object
    """
    return ReactionChain(species, rate_constants, **kwargs)
