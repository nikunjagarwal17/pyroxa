"""I/O helpers for Pyroxa MVP: Enhanced YAML/CTI spec loader with validation."""
import yaml
import os
import json
from typing import Dict, Any, List, Optional


class SpecValidationError(Exception):
    """Raised when a specification file is invalid."""
    pass


def validate_spec(spec: Dict[str, Any]) -> None:
    """Validate a simulation specification dictionary.
    
    Raises SpecValidationError if the spec is invalid.
    """
    # Check required top-level keys
    if 'system' not in spec:
        spec['system'] = 'WellMixed'  # default
    
    # Validate reaction specification
    if 'species' in spec and 'reactions' in spec:
        # Multi-species validation
        species = spec['species']
        if not isinstance(species, list) or len(species) == 0:
            raise SpecValidationError("'species' must be a non-empty list")
        
        reactions = spec['reactions']
        if not isinstance(reactions, list):
            raise SpecValidationError("'reactions' must be a list")
        
        for i, rxn in enumerate(reactions):
            if not isinstance(rxn, dict):
                raise SpecValidationError(f"Reaction {i} must be a dictionary")
            
            # Check rate constants
            if 'kf' not in rxn:
                rxn['kf'] = 1.0  # default
            if 'kr' not in rxn:
                rxn['kr'] = 0.0  # default
            
            # Validate reactants/products
            for key in ['reactants', 'products']:
                if key in rxn:
                    if not isinstance(rxn[key], dict):
                        raise SpecValidationError(f"Reaction {i} '{key}' must be a dictionary")
                    for species_name in rxn[key]:
                        if species_name not in species:
                            raise SpecValidationError(f"Unknown species '{species_name}' in reaction {i}")
    
    elif 'reaction' in spec:
        # Simple A <=> B validation
        reaction = spec['reaction']
        if not isinstance(reaction, dict):
            raise SpecValidationError("'reaction' must be a dictionary")
        
        if 'kf' not in reaction:
            reaction['kf'] = 1.0
        if 'kr' not in reaction:
            reaction['kr'] = 0.0
    
    else:
        # No reaction specified, add default
        spec['reaction'] = {'kf': 1.0, 'kr': 0.0}
    
    # Validate initial conditions
    if 'initial' not in spec:
        spec['initial'] = {}
    
    initial = spec['initial']
    if 'temperature' not in initial:
        initial['temperature'] = 300.0
    
    if 'conc' not in initial:
        if 'species' in spec:
            initial['conc'] = {s: 0.0 for s in spec['species']}
            if len(spec['species']) > 0:
                initial['conc'][spec['species'][0]] = 1.0  # default first species = 1.0
        else:
            initial['conc'] = {'A': 1.0, 'B': 0.0}
    
    # Validate simulation parameters
    if 'sim' not in spec:
        spec['sim'] = {}
    
    sim = spec['sim']
    if 'time_span' not in sim:
        sim['time_span'] = 10.0
    if 'time_step' not in sim:
        sim['time_step'] = 0.01


def load_spec_from_yaml(path: str) -> Dict[str, Any]:
    """Load a simulation spec from a YAML file and return a validated dict.

    If the extension is `.cti` a minimal CTI parser converts it to an internal dict.
    If the extension is `.json` it loads as JSON.
    
    Args:
        path: Path to the specification file
        
    Returns:
        Validated specification dictionary
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        SpecValidationError: If the specification is invalid
        yaml.YAMLError: If YAML parsing fails
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Specification file not found: {path}")
    
    _, ext = os.path.splitext(path)
    
    with open(path, 'r', encoding='utf-8') as f:
        txt = f.read()
    
    if ext.lower() == '.cti':
        # Enhanced CTI parser
        spec = parse_cti_content(txt)
    elif ext.lower() == '.json':
        spec = json.loads(txt)
    else:
        # YAML format
        try:
            spec = yaml.safe_load(txt)
        except yaml.YAMLError as e:
            raise SpecValidationError(f"Failed to parse YAML: {e}")
    
    if spec is None:
        raise SpecValidationError("Empty or invalid specification file")
    
    # Validate and normalize the spec
    validate_spec(spec)
    
    return spec


def parse_cti_content(content: str) -> Dict[str, Any]:
    """Parse CTI file content into a specification dictionary.
    
    This is a minimal parser for simple CTI files. For full CTI support,
    users should use Cantera directly.
    
    Supported formats:
    - species('A B C')  or  species = ['A', 'B', 'C']
    - reaction('A <=> B', [kf, kr])
    - reaction('A + B => C', kf)
    """
    spec = {'species': [], 'reactions': []}
    
    lines = content.splitlines()
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        try:
            if line.startswith('species'):
                # Parse species declaration
                spec['species'] = parse_cti_species(line)
            elif line.startswith('reaction'):
                # Parse reaction declaration
                rxn = parse_cti_reaction(line, spec['species'])
                if rxn:
                    spec['reactions'].append(rxn)
        except Exception as e:
            raise SpecValidationError(f"Error parsing CTI line {line_num}: {e}")
    
    return spec


def parse_cti_species(line: str) -> List[str]:
    """Parse a CTI species declaration."""
    # Handle both species('A B C') and species = ['A', 'B', 'C'] formats
    if '(' in line and ')' in line:
        # species('A B C') format
        inside = line[line.find('(')+1:line.rfind(')')]
        inside = inside.strip().strip("'\"")
        return inside.split()
    elif '=' in line and '[' in line:
        # species = ['A', 'B', 'C'] format
        right = line.split('=', 1)[1].strip()
        # Simple list parsing
        if right.startswith('[') and right.endswith(']'):
            items = right[1:-1].split(',')
            return [item.strip().strip("'\"") for item in items if item.strip()]
    
    return []


def parse_cti_reaction(line: str, species: List[str]) -> Optional[Dict[str, Any]]:
    """Parse a CTI reaction declaration."""
    if '(' not in line or ')' not in line:
        return None
    
    inside = line[line.find('(')+1:line.rfind(')')]
    
    # Split reaction equation and parameters
    if ',' in inside:
        rxn_expr, params_str = inside.split(',', 1)
        rxn_expr = rxn_expr.strip().strip("'\"")
        params_str = params_str.strip()
    else:
        rxn_expr = inside.strip().strip("'\"")
        params_str = ""
    
    # Parse parameters (rate constants)
    kf, kr = 1.0, 0.0
    if params_str:
        # Handle both [kf, kr] and single kf formats
        params_str = params_str.strip('[]')
        try:
            params = [float(x.strip()) for x in params_str.split(',') if x.strip()]
            if len(params) >= 1:
                kf = params[0]
            if len(params) >= 2:
                kr = params[1]
        except ValueError:
            pass  # Use defaults
    
    # Parse reaction equation
    reactants, products = parse_reaction_equation(rxn_expr, species)
    
    return {
        'kf': kf,
        'kr': kr,
        'reactants': reactants,
        'products': products
    }


def parse_reaction_equation(equation: str, species: List[str]) -> tuple:
    """Parse a reaction equation like 'A + B <=> C + D' or 'A => B'."""
    reactants = {}
    products = {}
    
    # Determine reaction type and split
    if '<=>' in equation:
        left, right = equation.split('<=>', 1)
        reversible = True
    elif '=>' in equation:
        left, right = equation.split('=>', 1)
        reversible = False
    elif '=' in equation:
        left, right = equation.split('=', 1)
        reversible = True
    else:
        # Simple A <=> B format
        parts = equation.split()
        if len(parts) >= 3 and parts[1] in ['<->', '<==>', '<=>']:
            left, right = parts[0], parts[2]
        else:
            raise SpecValidationError(f"Cannot parse reaction equation: {equation}")
    
    # Parse reactants and products
    reactants = parse_reaction_side(left.strip(), species)
    products = parse_reaction_side(right.strip(), species)
    
    return reactants, products


def parse_reaction_side(side: str, species: List[str]) -> Dict[str, int]:
    """Parse one side of a reaction equation (e.g., 'A + 2B')."""
    result = {}
    
    # Split by '+' and parse each term
    terms = [term.strip() for term in side.split('+')]
    
    for term in terms:
        if not term:
            continue
        
        # Extract coefficient and species name
        term = term.strip()
        coeff = 1
        species_name = term
        
        # Check for numeric coefficient
        parts = term.split()
        if len(parts) == 2 and parts[0].isdigit():
            coeff = int(parts[0])
            species_name = parts[1]
        elif len(parts) == 1:
            # Check if it starts with a number
            i = 0
            while i < len(term) and term[i].isdigit():
                i += 1
            if i > 0:
                coeff = int(term[:i])
                species_name = term[i:]
        
        # Validate species name
        if species_name not in species and species:  # only validate if species list is provided
            # Auto-add unknown species for flexibility
            pass
        
        result[species_name] = coeff
    
    return result


def save_results_to_csv(times: List[float], traj: List[List[float]], 
                       filename: str, species_names: Optional[List[str]] = None) -> None:
    """Save simulation results to a CSV file with proper headers."""
    import csv
    
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # Write header
        if species_names:
            header = ['time'] + species_names
        elif isinstance(traj[0], list) and len(traj[0]) == 2:
            header = ['time', 'A', 'B']
        else:
            header = ['time'] + [f'species_{i}' for i in range(len(traj[0]))]
        
        writer.writerow(header)
        
        # Write data
        for t, concentrations in zip(times, traj):
            if isinstance(concentrations[0], list):
                # Network case - flatten first reactor
                row = [t] + concentrations[0]
            else:
                row = [t] + concentrations
            writer.writerow(row)


def load_mechanism_from_yaml(path: str) -> Dict[str, Any]:
    """Load a chemical mechanism from a YAML file.
    
    This is an alias for load_spec_from_yaml for compatibility.
    """
    return load_spec_from_yaml(path)


# Backwards compatibility
parse_mechanism = load_spec_from_yaml
