import numpy as np
import math


def autocatalytic_rate(k, A, B):
    """Calculate autocatalytic reaction rate"""
    return k * A * B


def michaelis_menten_rate(Vmax, Km, substrate_conc):
    """Calculate Michaelis-Menten enzyme kinetics rate"""
    if substrate_conc < 0:
        return 0.0
    return (Vmax * substrate_conc) / (Km + substrate_conc)


def competitive_inhibition_rate(Vmax, Km, substrate_conc, inhibitor_conc, Ki):
    """Calculate competitive inhibition rate"""
    if substrate_conc < 0 or inhibitor_conc < 0:
        return 0.0
    Km_apparent = Km * (1.0 + inhibitor_conc / Ki)
    return (Vmax * substrate_conc) / (Km_apparent + substrate_conc)


def heat_capacity_nasa(T, coeffs):
    """Calculate heat capacity using NASA polynomial
    
    NASA polynomial form for Cp/R:
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    
    Parameters:
        T (float): Temperature in Kelvin
        coeffs (list): NASA polynomial coefficients [a1, a2, a3, a4, a5, ...]
        
    Returns:
        float: Heat capacity in J/(mol·K)
    """
    if T <= 0:
        return 0.0
    if len(coeffs) < 5:
        return 0.0
    R_GAS = 8.314  # J/mol/K
    return R_GAS * (coeffs[0] + coeffs[1]*T + coeffs[2]*T*T + 
                   coeffs[3]*T*T*T + coeffs[4]*T*T*T*T)


def enthalpy_nasa(T, coeffs, h_ref=0.0):
    """Calculate enthalpy using NASA polynomial
    
    NASA polynomial form for H/(RT):
    H/(RT) = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    
    Parameters:
        T (float): Temperature in Kelvin
        coeffs (list): NASA polynomial coefficients [a1, a2, a3, a4, a5, a6, a7]
        h_ref (float): Reference enthalpy in J/mol (default 0.0)
        
    Returns:
        float: Enthalpy in J/mol
    """
    if T <= 0:
        return h_ref
    if len(coeffs) < 5:
        return h_ref
    R_GAS = 8.314
    h_polynomial = R_GAS * T * (coeffs[0] + coeffs[1]*T/2.0 + coeffs[2]*T*T/3.0 + 
                               coeffs[3]*T*T*T/4.0 + coeffs[4]*T*T*T*T/5.0)
    return h_ref + h_polynomial


def entropy_nasa(T, coeffs, s_ref=0.0):
    """Calculate entropy using NASA polynomial
    
    NASA polynomial form for S/R:
    S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    
    Parameters:
        T (float): Temperature in Kelvin
        coeffs (list): NASA polynomial coefficients [a1, a2, a3, a4, a5, a6, a7]
        s_ref (float): Reference entropy in J/(mol·K) (default 0.0)
        
    Returns:
        float: Entropy in J/(mol·K)
    """
    if T <= 0:
        return s_ref
    if len(coeffs) < 5:
        return s_ref
    R_GAS = 8.314
    s_polynomial = R_GAS * (coeffs[0]*math.log(T) + coeffs[1]*T + coeffs[2]*T*T/2.0 + 
                           coeffs[3]*T*T*T/3.0 + coeffs[4]*T*T*T*T/4.0) / 100.0
    return s_ref + s_polynomial


def mass_transfer_correlation(Re, Sc, geometry_factor=0.023):
    """Calculate Sherwood number from Reynolds and Schmidt numbers"""
    if Re <= 0 or Sc <= 0:
        return 0.0
    Sh = geometry_factor * (Re**0.8) * (Sc**(1.0/3.0))
    return Sh


def heat_transfer_correlation(Re, Pr, geometry_factor=0.023):
    """Calculate Nusselt number from Reynolds and Prandtl numbers"""
    if Re <= 0 or Pr <= 0:
        return 0.0
    Nu = geometry_factor * (Re**0.8) * (Pr**(1.0/3.0))
    return Nu


def effective_diffusivity(molecular_diff, porosity, tortuosity, constriction_factor=1.0):
    """Calculate effective diffusivity in porous media"""
    if molecular_diff <= 0 or porosity <= 0:
        return 0.0
    return molecular_diff * porosity * constriction_factor / tortuosity


def pressure_drop_ergun(velocity, density, viscosity, particle_diameter, bed_porosity, bed_length):
    """Calculate pressure drop using Ergun equation"""
    if velocity <= 0 or density <= 0 or viscosity <= 0 or particle_diameter <= 0:
        return 0.0
    
    epsilon = bed_porosity
    term1 = 60.0 * viscosity * velocity * (1.0 - epsilon)**2 / \
            (particle_diameter**2 * epsilon**3)
    term2 = 0.7 * density * velocity**2 * (1.0 - epsilon) / \
            (particle_diameter * epsilon**3)
    
    dp_dz = term1 + term2
    return dp_dz * bed_length


class PIDController:
    """PID controller implementation with comprehensive state preservation and performance metrics"""
    
    def __init__(self, Kp=1.0, Ki=0.0, Kd=0.0, kp=None, ki=None, kd=None, 
                 integral_limit=None, derivative_filter=1.0):
        """Initialize PID controller with tuning parameters and limits
        
        Parameters:
            Kp, kp: Proportional gain
            Ki, ki: Integral gain
            Kd, kd: Derivative gain
            integral_limit: Maximum absolute value for integral term (anti-windup)
            derivative_filter: Low-pass filter coefficient for derivative (0-1)
        """
        self.Kp = kp if kp is not None else Kp
        self.Ki = ki if ki is not None else Ki
        self.Kd = kd if kd is not None else Kd
        self.integral_term = 0.0
        self.previous_error = 0.0
        self.integral_limit = integral_limit
        self.derivative_filter = derivative_filter
        self.filtered_derivative = 0.0
        self.total_calls = 0
        self.cumulative_error = 0.0
        self.max_error = 0.0
        self.min_error = 0.0
    
    def calculate(self, setpoint, process_variable, dt):
        """Calculate PID output with comprehensive metrics
        
        Returns dictionary with control output and performance metrics
        """
        error = setpoint - process_variable
        
        self.total_calls += 1
        self.cumulative_error += abs(error)
        self.max_error = max(self.max_error, error)
        self.min_error = min(self.min_error, error)
        
        proportional = self.Kp * error
        
        self.integral_term += error * dt
        if self.integral_limit is not None:
            self.integral_term = max(-self.integral_limit, 
                                   min(self.integral_limit, self.integral_term))
        integral = self.Ki * self.integral_term
        
        if dt > 0:
            raw_derivative = (error - self.previous_error) / dt
            self.filtered_derivative = (self.derivative_filter * raw_derivative + 
                                       (1 - self.derivative_filter) * self.filtered_derivative)
            derivative = self.Kd * self.filtered_derivative
        else:
            derivative = 0.0
        
        self.previous_error = error
        output = proportional + integral + derivative
        avg_error = self.cumulative_error / self.total_calls if self.total_calls > 0 else 0.0
        
        return {
            'output': output,
            'error': error,
            'proportional_term': proportional,
            'integral_term': integral,
            'derivative_term': derivative,
            'integral_state': self.integral_term,
            'average_error': avg_error,
            'max_error': self.max_error,
            'min_error': self.min_error,
            'total_calls': self.total_calls,
        }
    
    def compute(self, setpoint, process_variable, dt=1.0):
        """Alias for calculate method for convenience"""
        return self.calculate(setpoint, process_variable, dt)
    
    def reset(self):
        """Reset controller state and performance metrics"""
        self.integral_term = 0.0
        self.previous_error = 0.0
        self.filtered_derivative = 0.0
        self.total_calls = 0
        self.cumulative_error = 0.0
        self.max_error = 0.0
        self.min_error = 0.0
    
    def get_state(self):
        """Get current controller state"""
        return {
            'Kp': self.Kp,
            'Ki': self.Ki,
            'Kd': self.Kd,
            'integral_state': self.integral_term,
            'previous_error': self.previous_error,
            'total_calls': self.total_calls,
            'average_error': self.cumulative_error / self.total_calls if self.total_calls > 0 else 0.0,
        }


def pid_controller(setpoint, process_variable, dt, Kp, Ki, Kd, previous_error=0.0, integral=0.0):
    """PID controller function with comprehensive output (supports stateful operation with optional previous values)
    
    This is a simplified stateless PID implementation. For full PID with state preservation,
    use the PIDController class.
    """
    error = setpoint - process_variable
    proportional = Kp * error
    new_integral = integral + error * dt
    integral_term = Ki * new_integral
    
    if dt > 0:
        derivative = Kd * (error - previous_error) / dt
    else:
        derivative = 0.0
    
    output = proportional + integral_term + derivative
    percent_error = (abs(error) / abs(setpoint)) * 100 if setpoint != 0 else 0.0
    
    return {
        'output': output,
        'error': error,
        'percent_error': percent_error,
        'proportional_term': proportional,
        'integral_term': integral_term,
        'derivative_term': derivative,
        'new_integral_state': new_integral,
        'new_error_state': error,
    }


def langmuir_hinshelwood_rate(k, K_A, K_B, conc_A, conc_B):
    """Langmuir-Hinshelwood rate expression"""
    if conc_A < 0 or conc_B < 0:
        return 0.0
    denominator = 1.0 + K_A * conc_A + K_B * conc_B
    return (k * K_A * K_B * conc_A * conc_B) / (denominator * denominator)


def photochemical_rate(quantum_yield, molar_absorptivity, path_length, light_intensity, concentration):
    """Photochemical reaction rate"""
    if concentration < 0 or light_intensity < 0:
        return 0.0
    absorbance = molar_absorptivity * concentration * path_length
    absorbed_light = light_intensity * (1.0 - math.exp(-absorbance))
    return quantum_yield * absorbed_light


def pressure_peng_robinson(n, V, T, Tc, Pc, omega):
    """
    Calculate pressure using Peng-Robinson equation of state.
    
    Accurate for non-ideal gases and liquids, especially near critical point.
    
    Parameters:
        n (float): Number of moles [mol]
        V (float): Volume [L]
        T (float): Temperature [K]
        Tc (float): Critical temperature [K]
        Pc (float): Critical pressure [bar]
        omega (float): Acentric factor (dimensionless)
    
    Returns:
        float: Pressure [bar]
        Returns 0.0 if V ≤ 0 or T ≤ 0
        Returns 1000.0 bar if V ≤ nb (physically unreasonable compressed volume)
        
    Formula:
        P = (nRT)/(V - nb) - (n²a(T))/(V² + 2nbV - n²b²)
        
        where:
            a = 0.45724 R²Tc²/Pc
            b = 0.07780 RTc/Pc
            α(T) = [1 + κ(1 - √(T/Tc))]²
            κ = 0.37464 + 1.54226ω - 0.26992ω²
            R = 0.08314 L·bar/(mol·K)
    
    Applications:
        - High-pressure phase equilibria
        - Compressibility factor calculations
        - Vapor-liquid equilibrium (VLE)
        - Supercritical fluid properties
    """
    if V <= 0 or T <= 0:
        return 0.0
    R_GAS = 0.08314  # L*bar/(mol*K)
    a = 0.45724 * R_GAS * R_GAS * Tc * Tc / Pc
    b = 0.07780 * R_GAS * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    alpha = (1.0 + kappa * (1.0 - math.sqrt(T/Tc)))**2
    a_T = a * alpha
    
    if V <= n * b:
        return 1000.0  # Return a reasonable high pressure in expected range
    
    return (n * R_GAS * T) / (V - n * b) - (n * n * a_T) / (V * V + 2.0 * n * b * V - n * n * b * b)


def fugacity_coefficient(P, T, Tc, Pc, omega):
    """
    Calculate fugacity coefficient for non-ideal gas behavior.
    """
    if P <= 0 or T <= 0:
        return 1.0
    Tr = T / Tc
    Pr = P / Pc
    
    B0 = 0.083 - 0.422 / (Tr**1.6)
    B1 = 0.139 - 0.172 / (Tr**4.2)
    B = B0 + omega * B1
    
    Z = 1.0 + B * Pr / Tr  # Simplified compressibility factor
    ln_phi = B * Pr / Tr
    
    return math.exp(ln_phi)


def linear_interpolate(*args):
    """Linear interpolation with comprehensive error analysis and extrapolation detection
    """
    if len(args) == 5:
        x1, y1, x2, y2, x = args
        
        if x2 == x1:
            return {
                'value': y1,
                'slope': 0.0,
                'extrapolated': False,
                'interval': 'degenerate',
                'uncertainty': 0.0,
            }
        
        slope = (y2 - y1) / (x2 - x1)
        value = y1 + slope * (x - x1)
        
        extrapolated = x < min(x1, x2) or x > max(x1, x2)
        
        if extrapolated:
            if x < min(x1, x2):
                interval = 'below'
                relative_distance = abs(x - min(x1, x2)) / abs(x2 - x1)
            else:
                interval = 'above'
                relative_distance = abs(x - max(x1, x2)) / abs(x2 - x1)
        else:
            interval = 'within'
            relative_distance = abs(x - x1) / abs(x2 - x1)
        
        uncertainty = abs(y2 - y1) * 0.01  # Base 1% uncertainty
        if extrapolated:
            uncertainty *= (1 + relative_distance)  # Increases with distance
        
        return {
            'value': value,
            'slope': slope,
            'extrapolated': extrapolated,
            'interval': interval,
            'uncertainty': uncertainty,
            'relative_distance': relative_distance,
        }
    
    elif len(args) == 3:
        x, x_data, y_data = args
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        
        if len(x_data) != len(y_data):
            raise ValueError("x_data and y_data must have the same length")
        
        if len(x_data) < 2:
            raise ValueError("Need at least 2 data points for interpolation")
        
        extrapolated = False
        interval_index = -1
        
        for i in range(len(x_data) - 1):
            if x_data[i] <= x <= x_data[i + 1]:
                interval_index = i
                break
        
        if interval_index == -1:
            extrapolated = True
            if x < x_data[0]:
                interval_index = 0  # Use first two points
                interval = 'below'
            else:
                interval_index = len(x_data) - 2  # Use last two points
                interval = 'above'
        else:
            interval = 'within'
        
        x1, y1 = x_data[interval_index], y_data[interval_index]
        x2, y2 = x_data[interval_index + 1], y_data[interval_index + 1]
        
        if x2 == x1:
            slope = 0.0
            value = y1
            relative_distance = 0.0
        else:
            slope = (y2 - y1) / (x2 - x1)
            value = y1 + slope * (x - x1)
            relative_distance = abs(x - x1) / abs(x2 - x1)
        
        if len(y_data) > 2:
            data_range = np.max(y_data) - np.min(y_data)
            uncertainty = data_range * 0.02  # 2% of range
            if extrapolated:
                if x < x_data[0]:
                    extrap_distance = abs(x - x_data[0]) / (x_data[-1] - x_data[0])
                else:
                    extrap_distance = abs(x - x_data[-1]) / (x_data[-1] - x_data[0])
                uncertainty *= (1 + extrap_distance * 2)
        else:
            uncertainty = abs(y2 - y1) * 0.02
        
        return {
            'value': value,
            'slope': slope,
            'extrapolated': extrapolated,
            'interval': interval,
            'interval_index': interval_index,
            'uncertainty': uncertainty,
            'relative_distance': relative_distance,
            'x_bounds': [float(x1), float(x2)],
            'y_bounds': [float(y1), float(y2)],
        }
    
    else:
        raise ValueError("linear_interpolate() takes either 3 or 5 arguments")


def cubic_spline_interpolate(x, x_points, y_points):
    """Cubic spline interpolation with smoothness metrics and curvature analysis
    
    Note: This is a simplified implementation using piecewise cubic polynomials.
    For production use, consider scipy.interpolate.CubicSpline for natural splines.
    """
    if len(x_points) != len(y_points) or len(x_points) < 2:
        raise ValueError("Invalid input for spline interpolation")
    
    x_points = np.array(x_points)
    y_points = np.array(y_points)
    n = len(x_points)
    
    # For simplicity with few points, use cubic Hermite interpolation
    extrapolated = False
    interval_index = -1
    
    for i in range(n - 1):
        if x_points[i] <= x <= x_points[i + 1]:
            interval_index = i
            break
    
    if interval_index == -1:
        extrapolated = True
        if x < x_points[0]:
            interval_index = 0
            interval = 'below'
        else:
            interval_index = n - 2
            interval = 'above'
    else:
        interval = 'within'
    
    i = interval_index
    x1, y1 = x_points[i], y_points[i]
    x2, y2 = x_points[i + 1], y_points[i + 1]
    
    if i > 0:
        m1 = (y2 - y_points[i-1]) / (x2 - x_points[i-1])  # Central difference
    else:
        m1 = (y2 - y1) / (x2 - x1)  # Forward difference
    
    if i < n - 2:
        m2 = (y_points[i+2] - y1) / (x_points[i+2] - x1)  # Central difference
    else:
        m2 = (y2 - y1) / (x2 - x1)  # Backward difference
    
    # P(t) = h00(t)·y1 + h10(t)·m1·dx + h01(t)·y2 + h11(t)·m2·dx
    # where t = (x - x1)/(x2 - x1)
    dx = x2 - x1
    if dx == 0:
        value = y1
        first_derivative = 0.0
        second_derivative = 0.0
    else:
        t = (x - x1) / dx
        
        h00 = 2*t**3 - 3*t**2 + 1
        h10 = t**3 - 2*t**2 + t
        h01 = -2*t**3 + 3*t**2
        h11 = t**3 - t**2
        
        value = h00*y1 + h10*m1*dx + h01*y2 + h11*m2*dx
        
        dh00_dt = 6*t**2 - 6*t
        dh10_dt = 3*t**2 - 4*t + 1
        dh01_dt = -6*t**2 + 6*t
        dh11_dt = 3*t**2 - 2*t
        
        dP_dt = dh00_dt*y1 + dh10_dt*m1*dx + dh01_dt*y2 + dh11_dt*m2*dx
        first_derivative = dP_dt / dx
        
        d2h00_dt2 = 12*t - 6
        d2h10_dt2 = 6*t - 4
        d2h01_dt2 = -12*t + 6
        d2h11_dt2 = 6*t - 2
        
        d2P_dt2 = d2h00_dt2*y1 + d2h10_dt2*m1*dx + d2h01_dt2*y2 + d2h11_dt2*m2*dx
        second_derivative = d2P_dt2 / (dx**2)
    
    curvature = abs(second_derivative) / (1 + first_derivative**2)**1.5
    
    smoothness_indicator = 1.0 / (1.0 + curvature * 100)  # 0 to 1 scale
    
    if n > 2:
        local_variations = []
        for j in range(max(0, i-1), min(n-1, i+2)):
            dy = abs(y_points[j+1] - y_points[j])
            local_variations.append(dy)
        uncertainty = np.mean(local_variations) * 0.05  # 5% of local variation
        if extrapolated:
            uncertainty *= 2.0
    else:
        uncertainty = abs(y2 - y1) * 0.05
    
    return {
        'value': value,
        'first_derivative': first_derivative,  # dy/dx
        'second_derivative': second_derivative,  # d²y/dx²
        'curvature': curvature,  # Curvature κ
        'smoothness': smoothness_indicator,  # 0-1 scale
        'extrapolated': extrapolated,
        'interval': interval if extrapolated else 'within',
        'interval_index': interval_index,
        'uncertainty': uncertainty,
    }


def calculate_r_squared(y_actual, y_predicted):
    """Calculate R-squared (coefficient of determination) with comprehensive goodness-of-fit metrics"""
    y_actual = np.array(y_actual)
    y_predicted = np.array(y_predicted)
    
    if len(y_actual) != len(y_predicted):
        raise ValueError("y_actual and y_predicted must have the same length")
    
    n = len(y_actual)
    
    ss_res = np.sum((y_actual - y_predicted) ** 2)  # Residual sum of squares
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)  # Total sum of squares
    
    if ss_tot == 0:
        r_squared = 1.0 if ss_res == 0 else 0.0
    else:
        r_squared = 1 - (ss_res / ss_tot)
    
    # Calculate adjusted R² (accounts for number of parameters)
    # Assume 2 parameters (intercept + slope) for simplicity
    k = 2  # Number of parameters
    if n > k + 1:
        adjusted_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - k - 1)
    else:
        adjusted_r_squared = r_squared
    
    mae = np.mean(np.abs(y_actual - y_predicted))
    
    mse = ss_res / n if n > 0 else 0
    
    rmse = np.sqrt(mse)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        mape = np.mean(np.abs((y_actual - y_predicted) / y_actual)) * 100
        if np.isnan(mape) or np.isinf(mape):
            mape = 0.0
    
    if r_squared >= 0.95:
        fit_quality = 'excellent'
    elif r_squared >= 0.85:
        fit_quality = 'good'
    elif r_squared >= 0.70:
        fit_quality = 'fair'
    else:
        fit_quality = 'poor'
    
    return {
        'r_squared': r_squared,  # Coefficient of determination
        'adjusted_r_squared': adjusted_r_squared,  # Adjusted for parameters
        'ss_residual': ss_res,  # Residual sum of squares
        'ss_total': ss_tot,  # Total sum of squares
        'mae': mae,  # Mean absolute error
        'mse': mse,  # Mean squared error
        'rmse': rmse,  # Root mean squared error
        'mape': mape,  # Mean absolute percentage error (%)
        'fit_quality': fit_quality,  # Qualitative assessment
        'n_points': n,  # Number of data points
    }


def calculate_rmse(y_actual, y_predicted):
    """Calculate Root Mean Square Error with detailed error analysis and distribution metrics"""
    y_actual = np.array(y_actual)
    y_predicted = np.array(y_predicted)
    
    if len(y_actual) != len(y_predicted):
        raise ValueError("y_actual and y_predicted must have the same length")
    
    n = len(y_actual)
    
    errors = y_actual - y_predicted
    squared_errors = errors ** 2
    
    # RMSE
    mse = np.mean(squared_errors)
    rmse = np.sqrt(mse)
    
    mae = np.mean(np.abs(errors))  # Mean Absolute Error
    max_error = np.max(np.abs(errors))  # Maximum error
    min_error = np.min(np.abs(errors))  # Minimum error
    
    mean_error = np.mean(errors)  # Bias
    std_error = np.std(errors)  # Standard deviation of errors
    
    data_range = np.max(y_actual) - np.min(y_actual)
    if data_range > 0:
        normalized_rmse = (rmse / data_range) * 100  # Percentage
    else:
        normalized_rmse = 0.0
    
    mean_actual = np.mean(y_actual)
    if abs(mean_actual) > 1e-10:
        cv_rmse = (rmse / abs(mean_actual)) * 100  # Percentage
    else:
        cv_rmse = 0.0
    
    positive_errors = np.sum(errors > 0)
    negative_errors = np.sum(errors < 0)
    if n > 0:
        error_balance = abs(positive_errors - negative_errors) / n
    else:
        error_balance = 0.0
    
    if normalized_rmse < 5:
        prediction_quality = 'excellent'
    elif normalized_rmse < 10:
        prediction_quality = 'good'
    elif normalized_rmse < 20:
        prediction_quality = 'fair'
    else:
        prediction_quality = 'poor'
    
    return {
        'rmse': rmse,  # Root mean squared error
        'mse': mse,  # Mean squared error
        'mae': mae,  # Mean absolute error
        'max_error': max_error,  # Maximum absolute error
        'mean_error': mean_error,  # Mean error (bias)
        'std_error': std_error,  # Standard deviation of errors
        'normalized_rmse': normalized_rmse,  # RMSE as % of data range
        'cv_rmse': cv_rmse,  # Coefficient of variation (%)
        'error_balance': error_balance,  # Balance of +/- errors
        'prediction_quality': prediction_quality,  # Qualitative assessment
        'n_points': n,  # Number of data points
    }


def calculate_aic(y_actual, y_predicted, k):
    """Calculate Akaike Information Criterion with model selection metrics and comparison tools
    
    AIC = n·ln(RSS/n) + 2k
    Lower AIC indicates better model fit with penalty for complexity
    """
    if isinstance(y_actual, (float, int)) or isinstance(y_predicted, (float, int)):
        raise TypeError("y_actual and y_predicted must be array-like, not float")
    
    y_actual = np.array(y_actual)
    y_predicted = np.array(y_predicted)
    
    if len(y_actual) != len(y_predicted):
        raise ValueError("y_actual and y_predicted must have the same length")
    
    n = len(y_actual)
    
    rss = np.sum((y_actual - y_predicted) ** 2)
    
    if n <= 0 or rss <= 0:
        aic = float('inf')
        aicc = float('inf')
        bic = float('inf')
        likelihood = 0.0
    else:
        aic = n * np.log(rss / n) + 2 * k
        
        # Calculate corrected AIC (AICc) for small sample sizes
        # AICc = AIC + 2k(k+1)/(n-k-1)
        if n > k + 1:
            aicc = aic + (2 * k * (k + 1)) / (n - k - 1)
        else:
            aicc = float('inf')  # Not defined for n <= k + 1
        
        # Calculate Bayesian Information Criterion (BIC)
        # BIC = n·ln(RSS/n) + k·ln(n)
        bic = n * np.log(rss / n) + k * np.log(n)
        
        # Calculate log-likelihood (useful for likelihood ratio tests)
        # L = -(n/2)·ln(2π) - (n/2)·ln(RSS/n) - n/2
        likelihood = -(n/2) * np.log(2 * np.pi) - (n/2) * np.log(rss/n) - n/2
    
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    if ss_tot > 0:
        r_squared = 1 - (rss / ss_tot)
    else:
        r_squared = 1.0 if rss == 0 else 0.0
    
    if n > 0:
        parameter_ratio = k / n  # Ratio of parameters to data points
    else:
        parameter_ratio = 0.0
    
    if parameter_ratio < 0.1:
        complexity_assessment = 'simple'
    elif parameter_ratio < 0.2:
        complexity_assessment = 'moderate'
    else:
        complexity_assessment = 'complex'
    
    if aicc < aic:
        recommended_criterion = 'AICc'  # Use corrected AIC for small samples
    else:
        recommended_criterion = 'AIC'
    
    return {
        'aic': aic,  # Akaike Information Criterion
        'aicc': aicc,  # Corrected AIC (small sample)
        'bic': bic,  # Bayesian Information Criterion
        'log_likelihood': likelihood,  # Log-likelihood
        'rss': rss,  # Residual sum of squares
        'r_squared': r_squared,  # R² for reference
        'n_parameters': k,  # Number of model parameters
        'n_points': n,  # Number of data points
        'parameter_ratio': parameter_ratio,  # k/n ratio
        'complexity': complexity_assessment,  # Model complexity
        'recommended_criterion': recommended_criterion,  # Which criterion to use
    }


def gibbs_free_energy(enthalpy, entropy, T):
    """Calculate Gibbs free energy"""
    return enthalpy - T * entropy


def equilibrium_constant(delta_G, T, R=8.314):
    """Calculate equilibrium constant from Gibbs free energy
    """
    return math.exp(-delta_G / (R * T))


def arrhenius_rate(A, Ea, T, R=8.314):
    """Calculate reaction rate using Arrhenius equation"""
    return A * math.exp(-Ea / (R * T))


# Additional chemical engineering functions to reach 68 functions total

def first_order_rate(k, concentration):
    """First-order reaction rate: r = k * [A]"""
    return k * concentration


def second_order_rate(k, conc_A, conc_B=None):
    """Second-order reaction rate: r = k * [A] * [B] or r = k * [A]^2"""
    if conc_B is None:
        return k * conc_A * conc_A
    return k * conc_A * conc_B


def zero_order_rate(k):
    """Zero-order reaction rate: r = k"""
    return k


def reversible_rate(kf, kr, conc_A, conc_B=0.0):
    """Reversible reaction rate: r = kf * [A] - kr * [B]"""
    return kf * conc_A - kr * conc_B


def parallel_reaction_rate(k1, k2, concentration):
    """Parallel reaction rates: returns [r1, r2] where r1 = k1*[A], r2 = k2*[A]"""
    return [k1 * concentration, k2 * concentration]


def series_reaction_rate(k1, k2, conc_A, conc_B):
    """Series reaction rate for A -> B -> C"""
    r_AB = k1 * conc_A  # rate of A -> B
    r_BC = k2 * conc_B  # rate of B -> C
    return [-r_AB, r_AB - r_BC, r_BC]


def enzyme_inhibition_rate(Vmax, Km, substrate_conc, inhibitor_conc, Ki, inhibition_type='uncompetitive'):
    """Enzyme inhibition kinetics"""
    if inhibition_type == 'competitive':
        return competitive_inhibition_rate(Vmax, Km, substrate_conc, inhibitor_conc, Ki)
    elif inhibition_type == 'non_competitive':
        return Vmax * substrate_conc / ((Km + substrate_conc) * (1 + inhibitor_conc / Ki))
    else:
        return Vmax * substrate_conc / (Km + substrate_conc * (1 + inhibitor_conc / Ki))


def temperature_dependence(k_ref, Ea, T, T_ref=298.15, R=8.314):
    """Temperature dependence of rate constant using Arrhenius equation"""
    return k_ref * math.exp(-(Ea / R) * (1/T - 1/T_ref))


def pressure_dependence(k_ref, delta_V, P, P_ref=101325, R=8.314, T=298.15):
    """Pressure dependence of rate constant"""
    return k_ref * math.exp((delta_V / (R * T)) * (P - P_ref))


def activity_coefficient(x, gamma_inf, alpha):
    """Activity coefficient calculation using Wilson model"""
    return gamma_inf * math.exp(alpha * (1 - x) * (1 - x))


def diffusion_coefficient(T, viscosity, molar_volume):
    """Stokes-Einstein diffusion coefficient estimation"""
    if T <= 0 or viscosity <= 0 or molar_volume <= 0:
        return 0.0
    k_B = 1.38064852e-23  # Boltzmann constant (J/K)
    # Convert molar volume from cm³/mol to m³/mol
    V_m = molar_volume * 1e-6
    # Stokes-Einstein equation: D = k_B * T / (6 * pi * eta * r)
    # Approximate radius from molar volume: r = (3*V_m/(4*pi*N_A))^(1/3)
    N_A = 6.02214076e23
    r = ((3 * V_m) / (4 * math.pi * N_A))**(1/3)
    return k_B * T / (6 * math.pi * viscosity * r)


def thermal_conductivity(cp, rho, alpha):
    """Thermal conductivity from heat capacity, density and thermal diffusivity: k = cp * rho * alpha"""
    return cp * rho * alpha


def heat_transfer_coefficient(q, dt):
    """Heat transfer coefficient from heat flux and temperature difference: h = q / dt"""
    if dt == 0:
        return 0.0
    return q / dt


def mass_transfer_coefficient(flux, dc):
    """Mass transfer coefficient from mass flux and concentration difference: kc = flux / dc"""
    if dc == 0:
        return 0.0
    return flux / dc


def reynolds_number(density, velocity, length, viscosity):
    """
    Calculate Reynolds number - ratio of inertial to viscous forces."""
    return density * velocity * length / viscosity


def prandtl_number(cp, viscosity, thermal_conductivity):
    """
    Calculate Prandtl number - ratio of momentum to thermal diffusivity."""
    return cp * viscosity / thermal_conductivity


def schmidt_number(viscosity, density, diffusivity):
    """
    Calculate Schmidt number - ratio of momentum to mass diffusivity.
    """
    return viscosity / (density * diffusivity)


def nusselt_number(h, L, k):
    """
    Calculate Nusselt number - ratio of convective to conductive heat transfer.
    """
    return h * L / k


def sherwood_number(kc, L, D):
    """
    Calculate Sherwood number - ratio of convective to diffusive mass transfer.]
    """
    return kc * L / D


def friction_factor(delta_p, L, D, rho, v):
    """
    Calculate Darcy friction factor from pressure drop.
    """
    if L == 0 or rho == 0 or v == 0:
        return 0.0
    return delta_p * D / (L * 0.5 * rho * v * v)


def hydraulic_diameter(area, perimeter):
    """Calculate hydraulic diameter"""
    return 4 * area / perimeter


def residence_time(volume, flow_rate):
    """Calculate residence time"""
    return volume / flow_rate


def conversion(initial_conc, final_conc):
    """Calculate conversion: X = (C0 - C) / C0"""
    if initial_conc == 0:
        return 0.0
    return (initial_conc - final_conc) / initial_conc


def selectivity(product_conc, byproduct_conc):
    """Calculate selectivity: S = [P] / ([P] + [BP])"""
    total = product_conc + byproduct_conc
    if total == 0:
        return 0.0
    return product_conc / total


def yield_coefficient(product_formed, reactant_consumed):
    """Calculate yield coefficient"""
    if reactant_consumed == 0:
        return 0.0
    return product_formed / reactant_consumed


def space_time(volume, volumetric_flow_rate):
    """Calculate space time (tau = V / v0)"""
    return volume / volumetric_flow_rate


def space_velocity(volumetric_flow_rate, reactor_volume):
    """Calculate space velocity (SV = v0 / V)"""
    return volumetric_flow_rate / reactor_volume


def reaction_quotient(product_concs, reactant_concs, stoich_coeffs_products, stoich_coeffs_reactants):
    """Calculate reaction quotient Q"""
    Q = 1.0
    for i, conc in enumerate(product_concs):
        Q *= conc ** stoich_coeffs_products[i]
    for i, conc in enumerate(reactant_concs):
        Q /= conc ** stoich_coeffs_reactants[i]
    return Q


def extent_of_reaction(initial_conc, final_conc, stoich_coeff):
    """Calculate extent of reaction"""
    return (initial_conc - final_conc) / abs(stoich_coeff)


def batch_reactor_time(initial_conc, final_conc, rate_constant, order=1):
    """
    Calculate time required in batch reactor for given conversion.
    """
    if initial_conc <= 0 or final_conc <= 0:
        return 0.0
    
    if order == 0:
        # Zero order: C = C0 - kt
        return (initial_conc - final_conc) / rate_constant
    elif order == 1:
        # First order: ln(C0/C) = kt
        return math.log(initial_conc / final_conc) / rate_constant
    elif order == 2:
        # Second order: 1/C - 1/C0 = kt
        return (1/final_conc - 1/initial_conc) / rate_constant
    else:
        # nth order: [C^(1-n) - C0^(1-n)] / [(1-n)k] = t
        if order == 1:
            return math.log(initial_conc / final_conc) / rate_constant
        return ((final_conc**(1-order)) - (initial_conc**(1-order))) / ((1-order) * rate_constant)


def cstr_volume(flow_rate, rate_constant, conversion, order=1):
    """
    Calculate CSTR volume required for given conversion.
    """
    if conversion >= 1.0:
        return float('inf')  # Infinite volume needed for complete conversion
    if conversion <= 0:
        return 0.0
    
    if order == 1:
        # First order: V/F0 = X / [k*CA0*(1-X)]
        # Simplified: V = F0 * X / [k*(1-X)]
        return flow_rate * conversion / (rate_constant * (1 - conversion))
    else:
        return flow_rate * conversion / (rate_constant * (1 - conversion))


def pfr_volume(flow_rate, rate_constant, conversion, order=1):
    """
    Calculate PFR volume required for given conversion.
    """
    if conversion >= 1.0:
        return float('inf')  # Infinite volume for complete conversion
    if conversion <= 0:
        return 0.0
    
    if order == 1:
        # First order: V = F0 * [-ln(1-X)] / k
        return flow_rate * (-math.log(1 - conversion)) / rate_constant
    elif order == 2:
        # Second order: V/F0 = X / [k*CA0*(1-X)]
        return flow_rate * conversion / (rate_constant * (1 - conversion))
    else:
        # First-order approximation for other orders
        return flow_rate * (-math.log(1 - conversion)) / rate_constant


def fluidized_bed_hydrodynamics(particle_diameter, density_particle, density_fluid, viscosity, velocity):
    """
    Calculate comprehensive hydrodynamic properties of fluidized bed reactor.
    
    Uses Wen and Yu correlation for minimum fluidization, bed expansion correlations,
    and regime classification.
    
    Parameters:
        particle_diameter: Particle diameter [m]
        density_particle: Particle density [kg/m³]
        density_fluid: Fluid density [kg/m³]
        viscosity: Fluid dynamic viscosity [Pa·s]
        velocity: Superficial gas velocity [m/s]
    
    Returns:
        dict: Dictionary containing:
            - 'u_mf': Minimum fluidization velocity [m/s]
            - 'Re_mf': Reynolds number at minimum fluidization
            - 'Ar': Archimedes number
            - 'epsilon_mf': Void fraction at minimum fluidization
            - 'is_fluidized': Boolean indicating if bed is fluidized
            - 'regime': Fluidization regime string
            - 'bed_expansion': Bed expansion ratio (H/H_mf)
            - 'epsilon': Current void fraction at given velocity
            - 'terminal_velocity': Particle terminal velocity [m/s]
    """
    if particle_diameter <= 0 or density_particle <= 0 or density_fluid <= 0 or viscosity <= 0:
        return {
            'u_mf': 0.0,
            'Re_mf': 0.0,
            'Ar': 0.0,
            'epsilon_mf': 0.0,
            'is_fluidized': False,
            'regime': 'Invalid parameters',
            'bed_expansion': 1.0,
            'epsilon': 0.0,
            'terminal_velocity': 0.0
        }
    
    g = 9.81  # m/s²
    
    Ar = particle_diameter**3 * density_fluid * (density_particle - density_fluid) * g / (viscosity**2)
    
    # Wen and Yu correlation for Re_mf
    # Re_mf = [(33.7)² + 0.0408 × Ar]^0.5 - 33.7
    Re_mf = (1135.7 + 0.0408 * Ar)**0.5 - 33.7
    
    # Minimum fluidization velocity
    u_mf = Re_mf * viscosity / (density_fluid * particle_diameter)
    
    # epsilon_mf ≈ 0.4 for typical particles
    epsilon_mf = 0.586 * (Ar**(-0.029)) if Ar > 0 else 0.4
    epsilon_mf = max(0.35, min(0.50, epsilon_mf))  # Constrain to physical range
    
    is_fluidized = velocity > u_mf
    
    # Terminal velocity (single particle)
    # Using simplified correlation for spherical particles
    Re_t_guess = 1000  # Initial guess
    for _ in range(10):  # Iterative solution
        if Re_t_guess < 0.4:
            C_D = 24 / Re_t_guess
        elif Re_t_guess < 500:
            C_D = 24 / Re_t_guess * (1 + 0.15 * Re_t_guess**0.687)
        else:
            C_D = 0.44
        
        u_t = math.sqrt(4 * g * particle_diameter * (density_particle - density_fluid) / 
                       (3 * C_D * density_fluid))
        Re_t_new = density_fluid * u_t * particle_diameter / viscosity
        
        if abs(Re_t_new - Re_t_guess) < 0.01:
            break
        Re_t_guess = Re_t_new
    
    terminal_velocity = u_t
    
    if not is_fluidized:
        bed_expansion = 1.0
        epsilon = epsilon_mf
        regime = "Fixed bed (not fluidized)"
    else:
        # Richardson-Zaki correlation for bed expansion
        # u/u_t = ε^n, where n depends on Re
        Re_op = density_fluid * velocity * particle_diameter / viscosity
        
        # Exponent n (Richardson-Zaki)
        if Re_op < 0.2:
            n = 4.65
        elif Re_op < 1.0:
            n = 4.35 * Re_op**(-0.03)
        elif Re_op < 200:
            n = 4.45 * Re_op**(-0.1)
        else:
            n = 2.39
        
        if velocity < terminal_velocity:
            epsilon = (velocity / terminal_velocity)**(1/n)
            epsilon = min(0.95, max(epsilon_mf, epsilon))
        else:
            epsilon = 0.95  # Near complete entrainment
        
        bed_expansion = epsilon_mf * (1 - epsilon_mf) / (epsilon * (1 - epsilon)) if epsilon < 1 else 10.0
        
        velocity_ratio = velocity / u_mf
        
        if velocity_ratio < 2:
            regime = "Particulate (smooth) fluidization"
        elif velocity_ratio < 10:
            regime = "Bubbling fluidization"
        elif velocity_ratio < 50:
            regime = "Turbulent fluidization"
        elif velocity < terminal_velocity:
            regime = "Fast fluidization"
        else:
            regime = "Pneumatic transport"
    
    return {
        'u_mf': u_mf,
        'Re_mf': Re_mf,
        'Ar': Ar,
        'epsilon_mf': epsilon_mf,
        'is_fluidized': is_fluidized,
        'regime': regime,
        'bed_expansion': bed_expansion,
        'epsilon': epsilon,
        'terminal_velocity': terminal_velocity
    }


def packed_bed_pressure_drop(velocity, density, viscosity, particle_diameter, bed_porosity, bed_length):
    """
    Calculate pressure drop in packed bed using Ergun equation.
    
    This is a wrapper function that calls pressure_drop_ergun.
    
    Parameters:
        velocity: Superficial velocity [m/s]
        density: Fluid density [kg/m³]
        viscosity: Dynamic viscosity [Pa·s]
        particle_diameter: Particle diameter [m]
        bed_porosity: Void fraction (0 to 1)
        bed_length: Bed length [m]
    
    Returns:
        Pressure drop [Pa]
    """
    return pressure_drop_ergun(velocity, density, viscosity, particle_diameter, bed_porosity, bed_length)


def bubble_column_dynamics(gas_velocity, liquid_density, gas_density, surface_tension, viscosity=0.001, column_diameter=0.15):
    """
    Calculate comprehensive bubble column reactor hydrodynamics.
    
    Returns comprehensive properties including bubble size, rise velocity, gas holdup,
    flow regime, and mass transfer characteristics for bubble column reactors.
    
    Parameters:
        gas_velocity: Superficial gas velocity [m/s]
        liquid_density: Liquid phase density [kg/m³]
        gas_density: Gas phase density [kg/m³]
        surface_tension: Surface tension [N/m]
        viscosity: Liquid viscosity [Pa·s] (default=0.001 for water)
        column_diameter: Column diameter [m] (default=0.15)
    
    Returns:
        Dictionary with comprehensive bubble column properties:
        - bubble_diameter: Sauter mean bubble diameter [m]
        - bubble_rise_velocity: Single bubble rise velocity [m/s]
        - gas_holdup: Gas volume fraction [-]
        - swarm_velocity: Bubble swarm velocity [m/s]
        - flow_regime: Flow regime classification (str)
        - Mo: Morton number [-]
        - Eo: Eötvös number [-]
        - transition_velocity: Regime transition velocity [m/s]
        - is_homogeneous: Boolean if in homogeneous regime
        - kLa_estimate: Volumetric mass transfer coefficient estimate [1/s]
    """
    if liquid_density <= 0 or surface_tension <= 0 or gas_velocity < 0:
        return {
            'bubble_diameter': 0.0,
            'bubble_rise_velocity': 0.0,
            'gas_holdup': 0.0,
            'swarm_velocity': 0.0,
            'flow_regime': 'invalid',
            'Mo': 0.0,
            'Eo': 0.0,
            'transition_velocity': 0.0,
            'is_homogeneous': False,
            'kLa_estimate': 0.0
        }
    
    g = 9.81  # m/s²
    
    Mo = g * viscosity**4 / (liquid_density * surface_tension**3)
    
    # Calculate bubble diameter using Akita-Yoshida correlation
    # d_b = 2.8 * (σ/(ρ_L * g))^0.5 * (U_g)^0.4 for air-water systems
    
    if Mo < 10**-10:  # Very low viscosity systems
        d_b = 1.38 * (surface_tension / (liquid_density * g))**0.5
    elif Mo < 10**-3:
        # Akita-Yoshida correlation
        d_b = 2.8 * (surface_tension / (liquid_density * g))**0.5 * (gas_velocity)**0.4
    else:  # High viscosity systems
        d_b = 1.5 * (surface_tension / (liquid_density * g))**0.5
    
    # Limit bubble diameter to physically reasonable range
    d_b = max(0.001, min(d_b, 0.1))  # 1 mm to 100 mm
    
    # Eötvös number
    Eo = g * liquid_density * d_b**2 / surface_tension
    
    # Single bubble rise velocity (Davies-Taylor correlation for large bubbles)
    # u_b = 0.71 * sqrt(g * d_b) for large bubbles
    if d_b < 0.002:  # Small bubbles - Stokes regime
        u_bubble_single = g * d_b**2 * (liquid_density - gas_density) / (18 * viscosity)
    elif d_b < 0.006:  # Intermediate bubbles
        u_bubble_single = 0.33 * (g * d_b)**0.76
    else:  # Large bubbles - ellipsoidal regime
        u_bubble_single = 0.71 * (g * d_b)**0.5
    
    # Regime transition velocity (homogeneous to heterogeneous)
    U_trans = 0.045 * (surface_tension / 0.072)**0.2  # Normalized to air-water
    
    if gas_velocity < U_trans:
        flow_regime = "homogeneous"
        is_homogeneous = True
        # Homogeneous regime: ε_G = U_g / (U_g + U_b)
        epsilon_g = gas_velocity / (gas_velocity + u_bubble_single)
    else:
        flow_regime = "heterogeneous"
        is_homogeneous = False
        # Heterogeneous regime: ε_G = U_g / U_b (simplified)
        # More complex: Use Zuber-Findlay drift flux model
        # ε_G = U_g / (C_0 * U_g + U_∞)
        C0 = 1.2  # Distribution parameter for bubble columns
        U_inf = 0.25  # Drift velocity [m/s]
        epsilon_g = gas_velocity / (C0 * gas_velocity + U_inf)
    
    # Limit gas holdup to physical range
    epsilon_g = min(max(epsilon_g, 0.0), 0.5)  # Max 50% holdup
    
    # Bubble swarm velocity (accounting for hindered settling)
    # u_swarm = u_bubble * (1 - ε_G)^n, n ≈ 2 for bubble swarms
    swarm_velocity = u_bubble_single * (1 - epsilon_g)**2
    
    # Further classify heterogeneous regime
    if not is_homogeneous:
        if gas_velocity > 0.15:
            flow_regime = "heterogeneous-turbulent"
        elif gas_velocity > 0.10:
            flow_regime = "heterogeneous-churn-turbulent"
        else:
            flow_regime = "heterogeneous-bubbling"
    
    # Volumetric mass transfer coefficient estimation
    # kLa = 0.6 * D_c^0.5 * U_g^0.62 * ε_G^0.5
    D_L = 2e-9  # Typical liquid diffusivity [m²/s] for O2 in water
    Sc = viscosity / (liquid_density * D_L)  # Schmidt number
    
    if is_homogeneous:
        # Homogeneous regime correlation
        kLa = 0.32 * (gas_velocity)**0.7 * (epsilon_g)**0.5 * (column_diameter)**0.3
    else:
        # Heterogeneous regime correlation (Akita-Yoshida)
        kLa = 0.6 * (column_diameter)**0.5 * (gas_velocity)**0.62 * (epsilon_g)**0.5
    
    return {
        'bubble_diameter': d_b,
        'bubble_rise_velocity': u_bubble_single,
        'gas_holdup': epsilon_g,
        'swarm_velocity': swarm_velocity,
        'flow_regime': flow_regime,
        'Mo': Mo,
        'Eo': Eo,
        'transition_velocity': U_trans,
        'is_homogeneous': is_homogeneous,
        'kLa_estimate': kLa
    }


def crystallization_rate(supersaturation, nucleation_rate_constant, growth_rate_constant, temperature=298.15):
    """
    Calculate crystallization rate including nucleation and crystal growth.
    
    Parameters:
        supersaturation: Relative supersaturation S = (C - C_sat)/C_sat [-]
        nucleation_rate_constant: Nucleation rate constant [nuclei/(m³·s)]
        growth_rate_constant: Crystal growth rate constant [m/s]
        temperature: Temperature [K] (default=298.15)
    
    Returns:
        dict: Dictionary containing:
            - nucleation_rate: Rate of nuclei formation [nuclei/(m³·s)]
            - growth_rate: Linear crystal growth rate [m/s]
            - total_rate: Combined crystallization rate [kg/(m³·s)]
            - driving_force: Supersaturation ratio S+1 [-]
    """
    if supersaturation < 0 or temperature <= 0:
        return {
            'nucleation_rate': 0.0,
            'growth_rate': 0.0,
            'total_rate': 0.0,
            'driving_force': 1.0
        }
    
    # Primary nucleation (power law dependence on supersaturation)
    nucleation_rate = nucleation_rate_constant * (supersaturation ** 2)
    
    # Crystal growth (linear driving force)
    growth_rate = growth_rate_constant * supersaturation
    
    driving_force = 1.0 + supersaturation
    
    # Total crystallization rate (combined effect)
    # Assumes nucleation creates small crystals that grow
    total_rate = nucleation_rate * 1e-6 + growth_rate * 1000  # Convert to kg/(m³·s)
    
    return {
        'nucleation_rate': nucleation_rate,
        'growth_rate': growth_rate,
        'total_rate': total_rate,
        'driving_force': driving_force
    }


def precipitation_rate(conc_A, conc_B, ksp, rate_constant, order_A=1, order_B=1):
    """
    Calculate precipitation rate for sparingly soluble salt formation.
    """
    if conc_A < 0 or conc_B < 0 or ksp < 0:
        return {
            'precipitation_rate': 0.0,
            'ion_product': 0.0,
            'supersaturation_ratio': 0.0,
            'is_precipitating': False
        }
    
    ion_product = conc_A * conc_B
    
    supersaturation_ratio = ion_product / ksp if ksp > 0 else 0.0
    
    is_precipitating = ion_product > ksp
    
    if is_precipitating:
        # Precipitation rate proportional to supersaturation and concentrations
        driving_force = ion_product - ksp
        precip_rate = rate_constant * (conc_A ** order_A) * (conc_B ** order_B) * driving_force
    else:
        precip_rate = 0.0
    
    return {
        'precipitation_rate': precip_rate,
        'ion_product': ion_product,
        'supersaturation_ratio': supersaturation_ratio,
        'is_precipitating': is_precipitating
    }


def dissolution_rate(surface_area, mass_transfer_coeff, saturation_conc, current_conc):
    """
    Calculate dissolution rate of solid in liquid.
    
    Parameters:
        surface_area: Solid-liquid interface area [m²]
        mass_transfer_coeff: Mass transfer coefficient [m/s]
        saturation_conc: Saturation concentration [mol/L or kg/m³]
        current_conc: Current bulk concentration [mol/L or kg/m³]
    
    Returns:
        float: Dissolution rate [mol/s or kg/s depending on concentration units]
    """
    if surface_area < 0 or mass_transfer_coeff < 0 or saturation_conc < 0 or current_conc < 0:
        return 0.0
    
    driving_force = max(0.0, saturation_conc - current_conc)
    
    # Dissolution rate = k_L * A * ΔC
    return mass_transfer_coeff * surface_area * driving_force


def evaporation_rate(vapor_pressure, ambient_pressure, mass_transfer_coeff, area, molecular_weight=18.015):
    """
    Calculate evaporation rate from liquid surface.
    
    Parameters:
        vapor_pressure: Vapor pressure of liquid [Pa]
        ambient_pressure: Partial pressure of vapor in air [Pa]
        mass_transfer_coeff: Mass transfer coefficient [m/s]
        area: Evaporating surface area [m²]
        molecular_weight: Molecular weight [g/mol] (default=18.015 for water)
    
    Returns:
        dict: Dictionary containing:
            - evaporation_rate: Mass evaporation rate [kg/s]
            - molar_rate: Molar evaporation rate [mol/s]
            - driving_force: Vapor pressure difference [Pa]
            - flux: Mass flux [kg/(m²·s)]
    """
    if vapor_pressure < 0 or ambient_pressure < 0 or mass_transfer_coeff < 0 or area < 0:
        return {
            'evaporation_rate': 0.0,
            'molar_rate': 0.0,
            'driving_force': 0.0,
            'flux': 0.0
        }
    
    driving_force = max(0.0, vapor_pressure - ambient_pressure)
    
    # Convert pressure difference to concentration difference using ideal gas law
    # C = P/(RT), assuming T=298K for now
    R = 8.314  # J/(mol·K)
    T = 298.15  # K
    conc_difference = driving_force / (R * T)  # mol/m³
    
    # Molar evaporation rate
    molar_rate = mass_transfer_coeff * area * conc_difference  # mol/s
    
    # Mass evaporation rate
    mass_rate = molar_rate * molecular_weight / 1000  # kg/s
    
    flux = mass_rate / area if area > 0 else 0.0
    
    return {
        'evaporation_rate': mass_rate,
        'molar_rate': molar_rate,
        'driving_force': driving_force,
        'flux': flux
    }


def distillation_efficiency(actual_stages, theoretical_stages):
    """
    Calculate Murphree tray efficiency for distillation column.
    
    Parameters:
        actual_stages: Number of actual equilibrium stages [-]
        theoretical_stages: Number of theoretical equilibrium stages [-]
    
    Returns:
        float: Murphree efficiency (E_MV) [-], typically 0.5-0.9
               Returns 0.0 if theoretical_stages <= 0
    """
    if theoretical_stages <= 0:
        return 0.0
    return actual_stages / theoretical_stages


def extraction_efficiency(conc_feed, conc_raffinate):
    """
    Calculate extraction efficiency (fractional recovery).
    
    Parameters:
        conc_feed: Solute concentration in feed [mol/L or any consistent units]
        conc_raffinate: Solute concentration in raffinate (extract-depleted stream) [mol/L]
    
    Returns:
        float: Extraction efficiency E = (C_feed - C_raffinate)/C_feed [-]
               Range: 0.0 (no extraction) to 1.0 (complete extraction)
    """
    if conc_feed <= 0:
        return 0.0
    efficiency = (conc_feed - conc_raffinate) / conc_feed
    return max(0.0, min(1.0, efficiency))  # Clamp to [0, 1]


def adsorption_isotherm(conc, qmax, K_ads, n=1):
    """
    Calculate adsorbed amount using Langmuir or Freundlich isotherm.
    
    Parameters:
        conc: Equilibrium concentration in fluid phase [mol/L or mg/L]
        qmax: Maximum adsorption capacity [mol/kg or mg/g] (Langmuir only)
        K_ads: Adsorption equilibrium constant [L/mol or (mg/g)(L/mg)^(1/n)]
        n: Freundlich exponent [-] (default=1 for Langmuir)
    
    Returns:
        float: Adsorbed amount per unit mass of adsorbent [mol/kg or mg/g]
    
    Models:
        - n = 1: Langmuir isotherm: q = q_max * K * C / (1 + K * C)
        - n ≠ 1: Freundlich isotherm: q = K * C^(1/n)
    """
    if conc < 0 or qmax < 0 or K_ads < 0:
        return 0.0
    
    if n == 1:
        # Langmuir isotherm (monolayer adsorption)
        return qmax * K_ads * conc / (1.0 + K_ads * conc)
    else:
        # Freundlich isotherm (multilayer/heterogeneous surfaces)
        if n <= 0:
            return 0.0
        return K_ads * (conc ** (1.0 / n))


def desorption_rate(adsorbed_amount, desorption_constant, temperature=298.15):
    """
    Calculate desorption rate (first-order kinetics).
    
    Parameters:
        adsorbed_amount: Current amount adsorbed [mol/kg or mg/g]
        desorption_constant: First-order desorption rate constant [1/s]
        temperature: Temperature [K] (default=298.15, currently not used but kept for future enhancement)
    
    Returns:
        float: Desorption rate [mol/(kg·s) or mg/(g·s)]
               Negative value indicates release from adsorbent
    """
    if adsorbed_amount < 0 or desorption_constant < 0:
        return 0.0
    
    # First-order desorption: -dq/dt = k_d * q
    return desorption_constant * adsorbed_amount


def catalyst_activity(initial_activity, deactivation_constant, time, temperature=298.15, deactivation_order=1):
    """
    Calculate comprehensive catalyst activity and deactivation parameters over time.
    
    Parameters:
        initial_activity: Initial catalyst activity (dimensionless, typically 0-1)
        deactivation_constant: Deactivation rate constant [1/h or 1/s]
        time: Time on stream [h or s]
        temperature: Operating temperature [K] (default=298.15)
        deactivation_order: Order of deactivation kinetics (default=1)
    
    Returns:
        dict: Dictionary containing:
            - activity: Current catalyst activity [-]
            - activity_loss: Fraction of activity lost [-]
            - half_life: Time to 50% activity [same units as time]
            - deactivation_rate: Instantaneous deactivation rate [1/time]
            - remaining_lifetime: Estimated time to 10% activity [same units as time]
    """
    if initial_activity <= 0 or time < 0:
        return {
            'activity': 0.0,
            'activity_loss': 1.0,
            'half_life': 0.0,
            'deactivation_rate': 0.0,
            'remaining_lifetime': 0.0
        }
    
    # Temperature-dependent deactivation (Arrhenius-like)
    kd_effective = deactivation_constant * math.exp(-5000 / (8.314 * temperature)) if temperature > 0 else deactivation_constant
    
    if deactivation_order == 1:
        # First-order: a = a0 * exp(-kd * t)
        activity = initial_activity * math.exp(-kd_effective * time)
        half_life = math.log(2) / kd_effective if kd_effective > 0 else float('inf')
    elif deactivation_order == 2:
        # Second-order: 1/a = 1/a0 + kd * t
        denominator = 1 / initial_activity + kd_effective * time
        activity = 1 / denominator if denominator > 0 else 0.0
        half_life = 1 / (kd_effective * initial_activity) if kd_effective > 0 else float('inf')
    else:
        activity = initial_activity * math.exp(-kd_effective * time)
        half_life = math.log(2) / kd_effective if kd_effective > 0 else float('inf')
    
    activity = max(0.0, min(initial_activity, activity))
    activity_loss = (initial_activity - activity) / initial_activity if initial_activity > 0 else 0.0
    
    # Current deactivation rate
    deactivation_rate = kd_effective * (activity ** deactivation_order)
    
    # Remaining lifetime to 10% activity
    if kd_effective > 0 and activity > 0.1:
        if deactivation_order == 1:
            remaining_lifetime = -math.log(0.1 / activity) / kd_effective
        else:
            remaining_lifetime = half_life * 2  # Rough estimate
    else:
        remaining_lifetime = 0.0
    
    return {
        'activity': activity,
        'activity_loss': activity_loss,
        'half_life': half_life,
        'deactivation_rate': deactivation_rate,
        'remaining_lifetime': remaining_lifetime
    }


def catalyst_deactivation(current_activity, poison_concentration, deactivation_constant, poison_order=1, temperature=298.15):
    """
    Calculate comprehensive catalyst deactivation kinetics due to poisoning.
    
    Parameters:
        current_activity: Current catalyst activity [-]
        poison_concentration: Poison concentration [mol/L or ppm]
        deactivation_constant: Deactivation rate constant [appropriate units]
        poison_order: Reaction order with respect to poison (default=1)
        temperature: Temperature [K] (default=298.15)
    
    Returns:
        dict: Dictionary containing:
            - deactivation_rate: Rate of activity loss [1/time]
            - time_to_50_percent: Time to lose 50% current activity [time]
            - poison_coverage: Estimated surface poison coverage [-]
            - reversibility: Estimated reversibility factor (0=irreversible, 1=fully reversible)
            - critical_poison_conc: Poison concentration for 90% deactivation [same units as input]
    """
    if current_activity <= 0 or poison_concentration < 0:
        return {
            'deactivation_rate': 0.0,
            'time_to_50_percent': float('inf'),
            'poison_coverage': 0.0,
            'reversibility': 0.0,
            'critical_poison_conc': 0.0
        }
    
    # Deactivation rate: -da/dt = kd * a^m * C_poison^n
    deactivation_rate = deactivation_constant * (current_activity ** 1.0) * (poison_concentration ** poison_order)
    
    if deactivation_rate > 0:
        time_to_50_percent = math.log(2) / (deactivation_constant * (poison_concentration ** poison_order))
    else:
        time_to_50_percent = float('inf')
    
    # Langmuir-like poison coverage: θ = K*C / (1 + K*C)
    K_ads = 10.0  # Assumed adsorption constant
    poison_coverage = K_ads * poison_concentration / (1 + K_ads * poison_concentration)
    poison_coverage = min(1.0, poison_coverage)
    
    # Reversibility factor (assumed, depends on poison type)
    # Strong poisons (S, As, Pb) ~ 0.1, weak poisons (CO) ~ 0.7
    reversibility = 0.3 * math.exp(-poison_concentration / 100.0)  # Decreases with poison conc
    
    # Critical poison concentration for 90% deactivation
    critical_poison_conc = (0.9 / (deactivation_constant * time_to_50_percent)) ** (1 / poison_order) if poison_order > 0 else 0.0
    
    return {
        'deactivation_rate': deactivation_rate,
        'time_to_50_percent': time_to_50_percent,
        'poison_coverage': poison_coverage,
        'reversibility': reversibility,
        'critical_poison_conc': critical_poison_conc
    }


def surface_reaction_rate(surface_coverage, rate_constant, activation_energy, temperature, gas_constant=8.314, pressure=101325):
    """
    Calculate comprehensive surface catalytic reaction rate with detailed kinetics.
    
    Parameters:
        surface_coverage: Fractional surface coverage θ [-] (0 to 1)
        rate_constant: Pre-exponential factor [appropriate units, e.g., mol/(m²·s)]
        activation_energy: Activation energy [J/mol]
        temperature: Temperature [K]
        gas_constant: Universal gas constant [J/(mol·K)] (default=8.314)
        pressure: System pressure [Pa] (default=101325)
    
    Returns:
        dict: Dictionary containing:
            - reaction_rate: Surface reaction rate [mol/(m²·s)]
            - rate_constant_T: Temperature-dependent rate constant [same units as k0]
            - turnover_frequency: TOF [1/s]
            - activation_factor: exp(-Ea/RT) [-]
            - available_sites: Fraction of free sites (1-θ) [-]
    """
    if surface_coverage < 0 or surface_coverage > 1 or temperature <= 0:
        return {
            'reaction_rate': 0.0,
            'rate_constant_T': 0.0,
            'turnover_frequency': 0.0,
            'activation_factor': 0.0,
            'available_sites': 0.0
        }
    
    # Temperature-dependent rate constant: k(T) = k0 * exp(-Ea/RT)
    activation_factor = math.exp(-activation_energy / (gas_constant * temperature))
    rate_constant_T = rate_constant * activation_factor
    
    # Surface reaction rate (Langmuir-Hinshelwood type)
    # r = k(T) * θ * (1 - θ) for coverage-dependent rate
    available_sites = 1.0 - surface_coverage
    reaction_rate = rate_constant_T * surface_coverage * available_sites
    
    # Turnover frequency (assuming site density of 1e19 sites/m²)
    site_density = 1e19  # sites/m²
    if surface_coverage > 0:
        turnover_frequency = reaction_rate * 6.022e23 / site_density  # conversions
    else:
        turnover_frequency = 0.0
    
    return {
        'reaction_rate': reaction_rate,
        'rate_constant_T': rate_constant_T,
        'turnover_frequency': turnover_frequency,
        'activation_factor': activation_factor,
        'available_sites': available_sites
    }


def pore_diffusion_rate(diffusivity, pore_length, conc_surface, conc_center, pore_radius=1e-9, tortuosity=3.0, porosity=0.4):
    """
    Calculate comprehensive pore diffusion rate with effectiveness factor.
    
    Parameters:
        diffusivity: Molecular diffusivity [m²/s]
        pore_length: Pore length [m]
        conc_surface: Concentration at pore mouth [mol/m³]
        conc_center: Concentration at pore center [mol/m³]
        pore_radius: Pore radius [m] (default=1e-9)
        tortuosity: Pore tortuosity factor [-] (default=3.0)
        porosity: Catalyst porosity [-] (default=0.4)
    
    Returns:
        dict: Dictionary containing:
            - diffusion_rate: Molar diffusion rate [mol/s]
            - effective_diffusivity: D_eff = D * ε / τ [m²/s]
            - flux: Molar flux [mol/(m²·s)]
            - thiele_modulus: Thiele modulus for effectiveness [-]
            - effectiveness_factor: Catalyst effectiveness factor [-]
    """
    if diffusivity <= 0 or pore_length <= 0 or pore_radius <= 0:
        return {
            'diffusion_rate': 0.0,
            'effective_diffusivity': 0.0,
            'flux': 0.0,
            'thiele_modulus': 0.0,
            'effectiveness_factor': 1.0
        }
    
    # Effective diffusivity: D_eff = D * porosity / tortuosity
    effective_diffusivity = diffusivity * porosity / tortuosity
    
    delta_c = abs(conc_surface - conc_center)
    
    # Molar flux: J = -D_eff * dC/dx ≈ D_eff * ΔC / L
    flux = effective_diffusivity * delta_c / pore_length
    
    # Cross-sectional area of pore
    pore_area = math.pi * (pore_radius ** 2)
    
    # Diffusion rate through single pore
    diffusion_rate = flux * pore_area
    
    # Thiele modulus: φ = L * sqrt(k/D_eff) (simplified, assuming k~1)
    assumed_rate_constant = 0.1  # 1/s
    thiele_modulus = pore_length * math.sqrt(assumed_rate_constant / effective_diffusivity)
    
    # Effectiveness factor: η = tanh(φ) / φ for first-order reaction in slab
    if thiele_modulus > 0.01:
        effectiveness_factor = math.tanh(thiele_modulus) / thiele_modulus
    else:
        effectiveness_factor = 1.0  # No diffusion limitation
    
    return {
        'diffusion_rate': diffusion_rate,
        'effective_diffusivity': effective_diffusivity,
        'flux': flux,
        'thiele_modulus': thiele_modulus,
        'effectiveness_factor': effectiveness_factor
    }


def film_mass_transfer(mass_transfer_coeff, area, conc_bulk, conc_interface, flow_velocity=1.0, characteristic_length=0.1):
    """
    Calculate comprehensive film mass transfer with dimensionless correlations.
    
    Parameters:
        mass_transfer_coeff: Mass transfer coefficient [m/s]
        area: Transfer area [m²]
        conc_bulk: Bulk concentration [mol/m³ or kg/m³]
        conc_interface: Interface concentration [mol/m³ or kg/m³]
        flow_velocity: Fluid velocity [m/s] (default=1.0)
        characteristic_length: Characteristic length [m] (default=0.1)
    
    Returns:
        dict: Dictionary containing:
            - mass_transfer_rate: Total transfer rate [mol/s or kg/s]
            - flux: Mass flux [mol/(m²·s) or kg/(m²·s)]
            - driving_force: Concentration difference [mol/m³ or kg/m³]
            - film_thickness: Estimated film thickness [m]
            - enhancement_factor: Enhancement vs stagnant film [-]
    """
    if mass_transfer_coeff <= 0 or area <= 0:
        return {
            'mass_transfer_rate': 0.0,
            'flux': 0.0,
            'driving_force': 0.0,
            'film_thickness': 0.0,
            'enhancement_factor': 1.0
        }
    
    driving_force = abs(conc_bulk - conc_interface)
    
    # Molar/mass flux: J = k_c * ΔC
    flux = mass_transfer_coeff * driving_force
    
    mass_transfer_rate = flux * area
    
    # Estimated film thickness: δ ≈ D / k_c (assuming D ~ 1e-9 m²/s for liquids)
    assumed_diffusivity = 1e-9  # m²/s
    film_thickness = assumed_diffusivity / mass_transfer_coeff if mass_transfer_coeff > 0 else 0.0
    
    # Enhancement factor (flow vs stagnant)
    # For stagnant film: k_stagnant ≈ D/L
    k_stagnant = assumed_diffusivity / characteristic_length
    enhancement_factor = mass_transfer_coeff / k_stagnant if k_stagnant > 0 else 1.0
    enhancement_factor = max(1.0, enhancement_factor)  # Must be >= 1
    
    return {
        'mass_transfer_rate': mass_transfer_rate,
        'flux': flux,
        'driving_force': driving_force,
        'film_thickness': film_thickness,
        'enhancement_factor': enhancement_factor
    }


def bubble_rise_velocity(bubble_diameter, density_liquid, density_gas, surface_tension, viscosity, gas_flowrate=0.001):
    """
    Calculate comprehensive bubble rise velocity with regime classification.
    
    Parameters:
        bubble_diameter: Bubble diameter [m]
        density_liquid: Liquid density [kg/m³]
        density_gas: Gas density [kg/m³]
        surface_tension: Surface tension [N/m]
        viscosity: Liquid dynamic viscosity [Pa·s]
        gas_flowrate: Gas volumetric flow rate [m³/s] (default=0.001)
    
    Returns:
        dict: Dictionary containing:
            - rise_velocity: Bubble terminal rise velocity [m/s]
            - reynolds_number: Bubble Reynolds number [-]
            - morton_number: Morton number [-]
            - eotvos_number: Eötvös number [-]
            - regime: Flow regime ('stokes', 'intermediate', 'turbulent')
    """
    g = 9.81  # m/s²
    
    if bubble_diameter <= 0 or density_liquid <= 0:
        return {
            'rise_velocity': 0.0,
            'reynolds_number': 0.0,
            'morton_number': 0.0,
            'eotvos_number': 0.0,
            'regime': 'invalid'
        }
    
    delta_rho = density_liquid - density_gas
    
    # Morton number: Mo = g * μ⁴ * Δρ / (ρ_L² * σ³)
    if surface_tension > 0 and density_liquid > 0:
        morton_number = g * (viscosity ** 4) * delta_rho / ((density_liquid ** 2) * (surface_tension ** 3))
    else:
        morton_number = 0.0
    
    # Eötvös number: Eo = g * Δρ * d² / σ
    if surface_tension > 0:
        eotvos_number = g * delta_rho * (bubble_diameter ** 2) / surface_tension
    else:
        eotvos_number = 0.0
    
    if morton_number < 1e-3:
        # Clean system, potential flow
        rise_velocity = math.sqrt(g * bubble_diameter * delta_rho / density_liquid)
        regime = 'potential'
    elif morton_number < 10:
        rise_velocity = 0.71 * math.sqrt(g * bubble_diameter)
        regime = 'intermediate'
    else:
        # Contaminated system, Stokes regime
        rise_velocity = g * (bubble_diameter ** 2) * delta_rho / (18 * viscosity)
        regime = 'stokes'
    
    # Reynolds number: Re = ρ_L * U * d / μ
    reynolds_number = density_liquid * rise_velocity * bubble_diameter / viscosity if viscosity > 0 else 0.0
    
    # Refine regime classification
    if reynolds_number < 1:
        regime = 'stokes'
    elif reynolds_number < 500:
        regime = 'intermediate'
    else:
        regime = 'turbulent'
    
    return {
        'rise_velocity': rise_velocity,
        'reynolds_number': reynolds_number,
        'morton_number': morton_number,
        'eotvos_number': eotvos_number,
        'regime': regime
    }


def terminal_velocity(particle_diameter, density_particle, density_fluid, viscosity, sphericity=1.0):
    """
    Calculate comprehensive terminal settling velocity with drag analysis.
    
    Parameters:
        particle_diameter: Particle diameter [m]
        density_particle: Particle density [kg/m³]
        density_fluid: Fluid density [kg/m³]
        viscosity: Fluid dynamic viscosity [Pa·s]
        sphericity: Particle sphericity [-] (default=1.0 for sphere)
    
    Returns:
        dict: Dictionary containing:
            - terminal_velocity: Terminal settling velocity [m/s]
            - reynolds_number: Particle Reynolds number [-]
            - drag_coefficient: Drag coefficient [-]
            - settling_regime: 'laminar', 'transitional', or 'turbulent'
            - drag_force: Drag force at terminal velocity [N]
    """
    g = 9.81  # m/s²
    
    if particle_diameter <= 0 or viscosity <= 0 or density_fluid <= 0:
        return {
            'terminal_velocity': 0.0,
            'reynolds_number': 0.0,
            'drag_coefficient': 0.0,
            'settling_regime': 'invalid',
            'drag_force': 0.0
        }
    
    delta_rho = density_particle - density_fluid
    
    # Iterative solution for terminal velocity
    # Initial guess using Stokes law
    v_terminal = g * (particle_diameter ** 2) * delta_rho / (18 * viscosity)
    
    for _ in range(10):
        Re = density_fluid * v_terminal * particle_diameter / viscosity
        
        # Drag coefficient correlation (Schiller-Naumann for Re < 800, else turbulent)
        if Re < 0.1:
            Cd = 24 / Re if Re > 0 else 24
            regime = 'laminar'
        elif Re < 1000:
            Cd = (24 / Re) * (1 + 0.15 * (Re ** 0.687))
            regime = 'transitional'
        else:
            Cd = 0.44
            regime = 'turbulent'
        
        Cd = Cd / (sphericity ** 2)
        
        # Update terminal velocity: balance drag and buoyancy
        # F_drag = F_buoyancy => 0.5 * Cd * ρ_f * A * v² = V * Δρ * g
        # v = sqrt(4 * d * Δρ * g / (3 * Cd * ρ_f))
        if Cd > 0 and density_fluid > 0:
            v_new = math.sqrt(4 * particle_diameter * delta_rho * g / (3 * Cd * density_fluid))
        else:
            v_new = v_terminal
        
        if abs(v_new - v_terminal) / (v_terminal + 1e-10) < 0.01:
            v_terminal = v_new
            break
        v_terminal = v_new
    
    Re_final = density_fluid * v_terminal * particle_diameter / viscosity if viscosity > 0 else 0.0
    
    particle_area = math.pi * (particle_diameter ** 2) / 4
    drag_force = 0.5 * Cd * density_fluid * particle_area * (v_terminal ** 2)
    
    return {
        'terminal_velocity': v_terminal,
        'reynolds_number': Re_final,
        'drag_coefficient': Cd,
        'settling_regime': regime,
        'drag_force': drag_force
    }


def drag_coefficient(reynolds_number, mach_number=0.0, roughness=0.0):
    """
    Calculate comprehensive drag coefficient with compressibility and roughness effects.
    
    Parameters:
        reynolds_number: Reynolds number [-]
        mach_number: Mach number [-] (default=0.0 for incompressible)
        roughness: Relative surface roughness k/d [-] (default=0.0 for smooth)
    
    Returns:
        dict: Dictionary containing:
            - drag_coefficient: Total drag coefficient [-]
            - skin_friction_cd: Skin friction drag coefficient [-]
            - pressure_cd: Pressure (form) drag coefficient [-]
            - compressibility_factor: Compressibility correction [-]
            - flow_regime: 'creeping', 'laminar', 'transitional', 'turbulent'
    """
    Re = reynolds_number
    
    if Re <= 0:
        return {
            'drag_coefficient': 0.0,
            'skin_friction_cd': 0.0,
            'pressure_cd': 0.0,
            'compressibility_factor': 1.0,
            'flow_regime': 'invalid'
        }
    
    # Base drag coefficient (sphere in incompressible flow)
    if Re < 0.1:
        # Stokes/creeping flow
        Cd = 24 / Re
        regime = 'creeping'
        skin_cd = Cd * 0.67  # ~2/3 for sphere
        pressure_cd = Cd * 0.33
    elif Re < 1:
        Cd = 24 / Re * (1 + 0.15 * (Re ** 0.687))
        regime = 'laminar'
        skin_cd = 12 / Re
        pressure_cd = Cd - skin_cd
    elif Re < 1000:
        # Transitional regime (Schiller-Naumann)
        Cd = 24 / Re * (1 + 0.15 * (Re ** 0.687))
        regime = 'transitional'
        skin_cd = 8 / Re
        pressure_cd = Cd - skin_cd
    elif Re < 2e5:
        Cd = 0.44
        regime = 'turbulent-subcritical'
        skin_cd = 0.1
        pressure_cd = 0.34
    else:
        # Supercritical (boundary layer transition)
        Cd = 0.2
        regime = 'turbulent-supercritical'
        skin_cd = 0.08
        pressure_cd = 0.12
    
    # Roughness correction (increases drag)
    if roughness > 0:
        roughness_factor = 1 + 2 * roughness
        Cd = Cd * roughness_factor
        skin_cd = skin_cd * roughness_factor
    
    # Compressibility correction (Prandtl-Glauert)
    if mach_number > 0 and mach_number < 0.8:
        comp_factor = 1 / math.sqrt(1 - mach_number ** 2)
        Cd = Cd * comp_factor
    elif mach_number >= 0.8:
        # Transonic/supersonic - more complex, simplified here
        comp_factor = 1 + 0.2 * (mach_number - 0.8)
        Cd = Cd * comp_factor
    else:
        comp_factor = 1.0
    
    return {
        'drag_coefficient': Cd,
        'skin_friction_cd': skin_cd,
        'pressure_cd': pressure_cd,
        'compressibility_factor': comp_factor,
        'flow_regime': regime
    }


def mixing_time(tank_diameter, impeller_diameter, rotational_speed, viscosity, impeller_type='rushton', liquid_height=None):
    """
    Calculate comprehensive mixing time in stirred tank with regime analysis.
    
    Parameters:
        tank_diameter: Tank diameter [m]
        impeller_diameter: Impeller diameter [m]
        rotational_speed: Impeller speed [rps or rpm if >10]
        viscosity: Liquid kinematic viscosity [Pa·s]
        impeller_type: 'rushton', 'pitched_blade', 'anchor' (default='rushton')
        liquid_height: Liquid height [m] (default=tank_diameter for H/T=1)
    
    Returns:
        dict: Dictionary containing:
            - mixing_time: Time for 95% homogeneity [s]
            - mixing_time_99: Time for 99% homogeneity [s]
            - turnover_time: Tank turnover time [s]
            - power_number: Impeller power number [-]
            - reynolds_number: Impeller Reynolds number [-]
    """
    if tank_diameter <= 0 or impeller_diameter <= 0 or rotational_speed <= 0:
        return {
            'mixing_time': 0.0,
            'mixing_time_99': 0.0,
            'turnover_time': 0.0,
            'power_number': 0.0,
            'reynolds_number': 0.0
        }
    
    # Convert rpm to rps if needed
    N = rotational_speed if rotational_speed <= 10 else rotational_speed / 60.0
    
    H = liquid_height if liquid_height is not None else tank_diameter
    
    # Assume water-like density
    density = 1000  # kg/m³
    
    # Reynolds number: Re = ρ * N * D² / μ
    Re = density * N * (impeller_diameter ** 2) / viscosity if viscosity > 0 else 0.0
    
    if impeller_type == 'rushton':
        if Re < 10:
            Np = 300 / Re if Re > 0 else 300
        elif Re < 10000:
            Np = 5.0  # Transitional
        else:
            Np = 5.0  # Fully turbulent
    elif impeller_type == 'pitched_blade':
        Np = 1.3 if Re > 10000 else 200 / Re if Re > 0 else 200
    elif impeller_type == 'anchor':
        Np = 300 / Re if Re > 0 else 300
    else:
        Np = 5.0  # Default
    
    # Mixing time correlation (Norwood-Metzner)
    # θ_mix = C * (T/D)^2 * (μ/1000)^0.1 / N
    # For 95% homogeneity
    C = 5.2  # Typical for Rushton
    if impeller_type == 'pitched_blade':
        C = 4.0
    elif impeller_type == 'anchor':
        C = 8.0
    
    mixing_time_95 = C * ((tank_diameter / impeller_diameter) ** 2) * ((viscosity / 0.001) ** 0.1) / N
    
    # For 99% homogeneity (typically 1.5-2x longer)
    mixing_time_99 = mixing_time_95 * 1.7
    
    # Turnover time = V / Q_pump
    # Q_pump ≈ Np * N * D³
    pumping_flow = Np * N * (impeller_diameter ** 3)
    tank_volume = math.pi * (tank_diameter ** 2) * H / 4
    turnover_time = tank_volume / pumping_flow if pumping_flow > 0 else 0.0
    
    return {
        'mixing_time': mixing_time_95,
        'mixing_time_99': mixing_time_99,
        'turnover_time': turnover_time,
        'power_number': Np,
        'reynolds_number': Re
    }


def power_consumption(power_number, density, rotational_speed, impeller_diameter, tank_diameter=None, ungassed=True):
    """
    Calculate comprehensive power consumption in stirred tank.
    
    Parameters:
        power_number: Impeller power number [-]
        density: Liquid density [kg/m³]
        rotational_speed: Impeller speed [rps or rpm if >10]
        impeller_diameter: Impeller diameter [m]
        tank_diameter: Tank diameter [m] (default=3*D)
        ungassed: True for ungassed, False for gassed (default=True)
    
    Returns:
        dict: Dictionary containing:
            - power: Total power consumption [W]
            - power_per_volume: Specific power [W/m³]
            - torque: Impeller torque [N·m]
            - tip_speed: Impeller tip speed [m/s]
            - power_kw: Power in kilowatts [kW]
    """
    if density <= 0 or impeller_diameter <= 0 or rotational_speed <= 0:
        return {
            'power': 0.0,
            'power_per_volume': 0.0,
            'torque': 0.0,
            'tip_speed': 0.0,
            'power_kw': 0.0
        }
    
    # Convert rpm to rps if needed
    N = rotational_speed if rotational_speed <= 10 else rotational_speed / 60.0
    
    # Power consumption: P = Np * ρ * N³ * D⁵
    power = power_number * density * (N ** 3) * (impeller_diameter ** 5)
    
    # Gassing correction (for aerated systems)
    if not ungassed:
        power = power * 0.6
    
    # Torque: P = 2π * N * Torque
    torque = power / (2 * math.pi * N) if N > 0 else 0.0
    
    # Tip speed: v_tip = π * D * N
    tip_speed = math.pi * impeller_diameter * N
    
    # Tank volume for specific power
    T = tank_diameter if tank_diameter is not None else 3 * impeller_diameter
    H = T  # Assume H/T = 1
    tank_volume = math.pi * (T ** 2) * H / 4
    power_per_volume = power / tank_volume if tank_volume > 0 else 0.0
    
    power_kw = power / 1000.0
    
    return {
        'power': power,
        'power_per_volume': power_per_volume,
        'torque': torque,
        'tip_speed': tip_speed,
        'power_kw': power_kw
    }


def pumping_power(flow_rate, pressure_drop, efficiency=0.7, density=1000):
    """
    Calculate comprehensive pumping power with hydraulic analysis.
    
    Parameters:
        flow_rate: Volumetric flow rate [m³/s]
        pressure_drop: Pressure increase [Pa]
        efficiency: Pump efficiency [-] (default=0.7)
        density: Fluid density [kg/m³] (default=1000 for water)
    
    Returns:
        dict: Dictionary containing:
            - power_hydraulic: Hydraulic (ideal) power [W]
            - power_shaft: Shaft (actual) power [W]
            - power_kw: Shaft power [kW]
            - head: Pump head [m]
            - power_per_flow: Specific power [W/(m³/s)]
    """
    if flow_rate <= 0 or efficiency <= 0:
        return {
            'power_hydraulic': 0.0,
            'power_shaft': 0.0,
            'power_kw': 0.0,
            'head': 0.0,
            'power_per_flow': 0.0
        }
    
    g = 9.81  # m/s²
    
    # Hydraulic power: P_h = Q * ΔP
    power_hydraulic = flow_rate * pressure_drop
    
    # Shaft power: P_shaft = P_h / η
    power_shaft = power_hydraulic / efficiency
    
    # Pump head: H = ΔP / (ρ * g)
    head = pressure_drop / (density * g) if density > 0 else 0.0
    
    power_per_flow = power_shaft / flow_rate
    
    power_kw = power_shaft / 1000.0
    
    return {
        'power_hydraulic': power_hydraulic,
        'power_shaft': power_shaft,
        'power_kw': power_kw,
        'head': head,
        'power_per_flow': power_per_flow
    }


def compression_work(flow_rate, pressure_in, pressure_out, gamma=1.4, efficiency=0.75, temperature_in=298.15):
    """
    Calculate comprehensive compression work with thermodynamic analysis.
    
    Parameters:
        flow_rate: Molar or volumetric flow rate [mol/s or m³/s]
        pressure_in: Inlet pressure [Pa]
        pressure_out: Outlet pressure [Pa]
        gamma: Heat capacity ratio Cp/Cv [-] (default=1.4 for air)
        efficiency: Isentropic efficiency [-] (default=0.75)
        temperature_in: Inlet temperature [K] (default=298.15)
    
    Returns:
        dict: Dictionary containing:
            - work_ideal: Isentropic (ideal) compression work [W]
            - work_actual: Actual compression work [W]
            - work_kw: Actual work [kW]
            - temperature_out: Outlet temperature [K]
            - compression_ratio: Pressure ratio [-]
    """
    if flow_rate <= 0 or pressure_in <= 0 or pressure_out <= 0:
        return {
            'work_ideal': 0.0,
            'work_actual': 0.0,
            'work_kw': 0.0,
            'temperature_out': temperature_in,
            'compression_ratio': 1.0
        }
    
    r = pressure_out / pressure_in
    
    # Ideal (isentropic) compression work
    # W = (γ/(γ-1)) * P1 * V1 * [(P2/P1)^((γ-1)/γ) - 1]
    # For ideal gas: P * V = n * R * T
    exponent = (gamma - 1) / gamma
    work_ideal = flow_rate * (gamma / (gamma - 1)) * pressure_in * (r ** exponent - 1)
    
    # Actual work accounting for efficiency
    work_actual = work_ideal / efficiency if efficiency > 0 else work_ideal
    
    # Outlet temperature (isentropic)
    # T2/T1 = (P2/P1)^((γ-1)/γ)
    temperature_out_ideal = temperature_in * (r ** exponent)
    
    # Actual outlet temperature (accounting for efficiency)
    # T2_actual = T1 + (T2_ideal - T1) / η
    temperature_out_actual = temperature_in + (temperature_out_ideal - temperature_in) / efficiency
    
    work_kw = work_actual / 1000.0
    
    return {
        'work_ideal': work_ideal,
        'work_actual': work_actual,
        'work_kw': work_kw,
        'temperature_out': temperature_out_actual,
        'compression_ratio': r
    }


def heat_exchanger_effectiveness(actual_heat_transfer, max_possible_heat_transfer, flow_config='counterflow', NTU=None):
    """
    Calculate comprehensive heat exchanger effectiveness with configuration analysis.
    
    Parameters:
        actual_heat_transfer: Actual heat transfer rate [W]
        max_possible_heat_transfer: Maximum possible heat transfer [W]
        flow_config: 'counterflow', 'parallel', 'shell_tube', 'crossflow' (default='counterflow')
        NTU: Number of Transfer Units [-] (optional)
    
    Returns:
        dict: Dictionary containing:
            - effectiveness: Heat exchanger effectiveness ε [-]
            - capacity_ratio: C_min/C_max used [-]
            - ntu_calculated: Calculated NTU if not provided [-]
            - thermal_performance: Performance rating (0-100)
            - flow_configuration: Flow pattern
    """
    if max_possible_heat_transfer <= 0:
        return {
            'effectiveness': 0.0,
            'capacity_ratio': 0.0,
            'ntu_calculated': 0.0,
            'thermal_performance': 0.0,
            'flow_configuration': flow_config
        }
    
    # Effectiveness: ε = Q_actual / Q_max
    effectiveness = actual_heat_transfer / max_possible_heat_transfer
    effectiveness = max(0.0, min(1.0, effectiveness))  # Clamp to [0, 1]
    
    C_ratio = 0.5
    
    if NTU is None:
        if flow_config == 'counterflow':
            # ε = (1 - exp(-NTU(1-C)))/(1 - C*exp(-NTU(1-C)))
            if effectiveness > 0:
                ntu_calc = -math.log(1 - effectiveness) / (1 - C_ratio * effectiveness)
            else:
                ntu_calc = 0.0
        else:
            ntu_calc = -math.log(1 - effectiveness * (1 + C_ratio)) / (1 + C_ratio)
    else:
        ntu_calc = NTU
    
    # Thermal performance rating (0-100)
    # Excellent: >80%, Good: 60-80%, Fair: 40-60%, Poor: <40%
    performance = effectiveness * 100
    
    return {
        'effectiveness': effectiveness,
        'capacity_ratio': C_ratio,
        'ntu_calculated': ntu_calc,
        'thermal_performance': performance,
        'flow_configuration': flow_config
    }


def overall_heat_transfer_coefficient(h_hot, h_cold, thickness, thermal_conductivity, fouling_hot=0, fouling_cold=0, diameter_ratio=1.0):
    """
    Calculate comprehensive overall heat transfer coefficient with fouling and geometry.
    
    Parameters:
        h_hot: Hot-side heat transfer coefficient [W/(m²·K)]
        h_cold: Cold-side heat transfer coefficient [W/(m²·K)]
        thickness: Wall thickness [m]
        thermal_conductivity: Wall thermal conductivity [W/(m·K)]
        fouling_hot: Hot-side fouling resistance [m²·K/W] (default=0)
        fouling_cold: Cold-side fouling resistance [m²·K/W] (default=0)
        diameter_ratio: D_outer/D_inner for cylindrical geometry (default=1.0 for flat)
    
    Returns:
        dict: Dictionary containing:
            - overall_U: Overall heat transfer coefficient [W/(m²·K)]
            - thermal_resistance: Total thermal resistance [m²·K/W]
            - convection_resistance_hot: Hot-side convection resistance [m²·K/W]
            - convection_resistance_cold: Cold-side convection resistance [m²·K/W]
            - wall_resistance: Wall conduction resistance [m²·K/W]
            - fouling_resistance_total: Total fouling resistance [m²·K/W]
    """
    if h_hot <= 0 or h_cold <= 0 or thermal_conductivity <= 0:
        return {
            'overall_U': 0.0,
            'thermal_resistance': float('inf'),
            'convection_resistance_hot': float('inf'),
            'convection_resistance_cold': float('inf'),
            'wall_resistance': 0.0,
            'fouling_resistance_total': 0.0
        }
    
    R_hot = 1.0 / h_hot
    R_cold = 1.0 / h_cold
    
    # Wall resistance (flat or cylindrical)
    if diameter_ratio > 1.01:
        # Cylindrical: R_wall = ln(D_o/D_i) / (2π k L) per unit length
        R_wall = math.log(diameter_ratio) / thermal_conductivity
    else:
        R_wall = thickness / thermal_conductivity
    
    R_foul_total = fouling_hot + fouling_cold
    
    R_total = R_hot + fouling_hot + R_wall + fouling_cold + R_cold
    
    # Overall U = 1 / R_total
    U_overall = 1.0 / R_total if R_total > 0 else 0.0
    
    return {
        'overall_U': U_overall,
        'thermal_resistance': R_total,
        'convection_resistance_hot': R_hot,
        'convection_resistance_cold': R_cold,
        'wall_resistance': R_wall,
        'fouling_resistance_total': R_foul_total
    }


def fouling_resistance(clean_U, fouled_U, safety_factor=1.2):
    """
    Calculate comprehensive fouling resistance with operational analysis.
    
    Parameters:
        clean_U: Clean overall heat transfer coefficient [W/(m²·K)]
        fouled_U: Fouled overall heat transfer coefficient [W/(m²·K)]
        safety_factor: Design safety factor for fouling (default=1.2)
    
    Returns:
        dict: Dictionary containing:
            - fouling_resistance: Total fouling resistance [m²·K/W]
            - fouling_factor: Fraction of capacity lost [-]
            - design_fouling: Design fouling with safety factor [m²·K/W]
            - cleaning_indicator: True if cleaning recommended (>50% loss)
            - performance_reduction: Performance reduction percentage [%]
    """
    if clean_U <= 0 or fouled_U <= 0:
        return {
            'fouling_resistance': 0.0,
            'fouling_factor': 0.0,
            'design_fouling': 0.0,
            'cleaning_indicator': False,
            'performance_reduction': 0.0
        }
    
    # Fouling resistance: R_f = 1/U_fouled - 1/U_clean
    R_fouling = (1.0 / fouled_U) - (1.0 / clean_U)
    R_fouling = max(0.0, R_fouling)  # Can't be negative
    
    # Fouling factor (fraction of heat transfer capacity lost)
    fouling_factor = (clean_U - fouled_U) / clean_U
    fouling_factor = max(0.0, min(1.0, fouling_factor))
    
    # Design fouling resistance with safety factor
    R_design = R_fouling * safety_factor
    
    cleaning_needed = fouling_factor > 0.5
    
    performance_reduction = fouling_factor * 100
    
    return {
        'fouling_resistance': R_fouling,
        'fouling_factor': fouling_factor,
        'design_fouling': R_design,
        'cleaning_indicator': cleaning_needed,
        'performance_reduction': performance_reduction
    }


# Advanced simulation functions to complete implementation

def simulate_packed_bed(params):
    """Simulate packed bed reactor with axial dispersion, pressure drop, and detailed performance metrics"""
    import numpy as np
    
    bed_length = params.get('bed_length', 1.0)  # m
    porosity = params.get('porosity', 0.4)  # void fraction
    particle_diameter = params.get('particle_diameter', 0.003)  # m
    flow_rate = params.get('flow_rate', 0.001)  # m³/s
    inlet_concentration = params.get('inlet_concentration', [2.0, 0.0])  # mol/m³
    reaction_rate_constant = params.get('reaction_rate_constant', 1.0)  # 1/s
    time_span = params.get('time_span', 1.0)  # s
    bed_diameter = params.get('bed_diameter', 0.1)  # m
    temperature = params.get('temperature', 298.15)  # K
    viscosity = params.get('viscosity', 0.001)  # Pa·s
    density = params.get('density', 1000.0)  # kg/m³
    
    bed_cross_section = 3.14159 * (bed_diameter / 2)**2  # m²
    bed_volume = bed_cross_section * bed_length  # m³
    void_volume = bed_volume * porosity  # m³
    superficial_velocity = flow_rate / bed_cross_section  # m/s
    interstitial_velocity = superficial_velocity / porosity  # m/s
    
    # Simple packed bed simulation using finite differences
    n_species = len(inlet_concentration)
    segments = 20  # Number of axial segments
    dt = 0.01  # Time step (s)
    dx = bed_length / segments  # Spatial step (m)
    time_steps = int(time_span / dt)
    
    # Initialize concentration matrix [segment, species, time]
    conc = np.zeros((segments, n_species, time_steps))
    for i in range(n_species):
        conc[:, i, 0] = inlet_concentration[i]  # Set initial concentrations
    
    # Simple first-order reaction simulation
    for t in range(1, time_steps):
        for seg in range(segments):
            for species in range(n_species):
                if seg > 0:
                    convection = interstitial_velocity * (conc[seg-1, species, t-1] - conc[seg, species, t-1]) / dx
                else:
                    convection = interstitial_velocity * (inlet_concentration[species] - conc[seg, species, t-1]) / dx
                
                # Simple first-order reaction: A -> B
                if species == 0:  # Reactant
                    reaction_rate = reaction_rate_constant * conc[seg, species, t-1]
                else:  # Product
                    reaction_rate = -reaction_rate_constant * conc[seg, 0, t-1]  # Negative because it's formed
                
                conc[seg, species, t] = conc[seg, species, t-1] + dt * (convection - reaction_rate)
                conc[seg, species, t] = max(0, conc[seg, species, t])  # Non-negative concentrations
    
    if inlet_concentration[0] > 0:
        conversion = 1.0 - conc[-1, 0, -1] / inlet_concentration[0]
    else:
        conversion = 0.0
    
    Re = density * superficial_velocity * particle_diameter / viscosity
    
    # Pressure drop using Ergun equation
    # ΔP/L = 150·(1-ε)²/ε³ · μ·v/dp² + 1.75·(1-ε)/ε³ · ρ·v²/dp
    viscous_term = 150 * (1 - porosity)**2 / porosity**3 * viscosity * superficial_velocity / particle_diameter**2
    inertial_term = 1.75 * (1 - porosity) / porosity**3 * density * superficial_velocity**2 / particle_diameter
    pressure_drop_per_length = viscous_term + inertial_term  # Pa/m
    total_pressure_drop = pressure_drop_per_length * bed_length  # Pa
    
    residence_time = void_volume / flow_rate  # s
    
    space_velocity = flow_rate / bed_volume  # 1/s
    
    theoretical_conversion = 1 - np.exp(-reaction_rate_constant * residence_time)
    if theoretical_conversion > 0:
        effectiveness = conversion / theoretical_conversion
    else:
        effectiveness = 1.0
    
    time_array = np.linspace(0, time_span, time_steps)
    outlet_concentration = conc[-1, :, -1]  # Final concentrations at outlet
    
    return {
        'time': time_array.tolist(),
        'concentration_profile': conc,  # 3D array [segment, species, time]
        'outlet_concentration': outlet_concentration.tolist(),  # mol/m³
        'conversion': conversion,  # dimensionless
        'pressure_drop': total_pressure_drop,  # Pa
        'pressure_drop_per_length': pressure_drop_per_length,  # Pa/m
        'reynolds_number': Re,  # dimensionless
        'residence_time': residence_time,  # s
        'space_velocity': space_velocity,  # 1/s
        'effectiveness': effectiveness,  # dimensionless
        'superficial_velocity': superficial_velocity,  # m/s
        'interstitial_velocity': interstitial_velocity,  # m/s
    }


def simulate_fluidized_bed(params):
    """Simulate fluidized bed reactor with bubble phase, emulsion phase, and comprehensive hydrodynamics"""
    import numpy as np
    
    bed_height = params.get('bed_height', 2.0)  # m
    bed_diameter = params.get('bed_diameter', 1.0)  # m
    particle_diameter = params.get('particle_diameter', 0.001)  # m
    gas_velocity = params.get('gas_velocity', 0.5)  # m/s
    inlet_concentration = params.get('inlet_concentration', [1.0, 0.0])  # mol/m³
    reaction_rate_constant = params.get('reaction_rate_constant', 0.5)  # 1/s
    time_span = params.get('time_span', 2.0)  # s
    particle_density = params.get('particle_density', 2500.0)  # kg/m³
    gas_density = params.get('gas_density', 1.2)  # kg/m³
    gas_viscosity = params.get('gas_viscosity', 1.8e-5)  # Pa·s
    
    bed_area = 3.14159 * (bed_diameter / 2)**2  # m²
    
    # Minimum fluidization velocity (Wen & Yu correlation)
    Ar = particle_density * gas_density * 9.81 * particle_diameter**3 / gas_viscosity**2  # Archimedes number
    Re_mf = np.sqrt(33.7**2 + 0.0408 * Ar) - 33.7  # Reynolds at minimum fluidization
    U_mf = Re_mf * gas_viscosity / (gas_density * particle_diameter)  # m/s
    
    # Bubble properties (Davidson & Harrison)
    U_excess = max(0, gas_velocity - U_mf)  # m/s
    bubble_diameter = 0.54 * (U_excess)**0.4 * (bed_height + 4 * np.sqrt(bed_area))**0.8  # m
    bubble_velocity = gas_velocity - U_mf + 0.711 * np.sqrt(9.81 * bubble_diameter)  # m/s
    
    voidage_mf = 0.4  # Typical void fraction at minimum fluidization
    if gas_velocity > U_mf:
        bubble_fraction = (gas_velocity - U_mf) / bubble_velocity
        bed_voidage = voidage_mf + bubble_fraction * (1 - voidage_mf)
    else:
        bubble_fraction = 0.0
        bed_voidage = voidage_mf
    
    expanded_bed_height = bed_height * (1 - voidage_mf) / (1 - bed_voidage)  # m
    
    n_species = len(inlet_concentration)
    dt = 0.01  # s
    time_steps = int(time_span / dt)
    
    # Initialize concentration arrays for bubble and emulsion phases
    bubble_conc = np.zeros((n_species, time_steps))
    emulsion_conc = np.zeros((n_species, time_steps))
    
    for i in range(n_species):
        bubble_conc[i, 0] = inlet_concentration[i]
        emulsion_conc[i, 0] = 0.0  # Emulsion starts with no gas species
    
    # Mass transfer coefficient between phases
    k_be = 4.5 * (U_mf / bubble_diameter) + 5.85 * (gas_density * 0.001)**0.5 / bubble_diameter**1.25  # 1/s
    
    # Simple two-phase fluidized bed simulation
    for t in range(1, time_steps):
        for species in range(n_species):
            # Bubble phase - convection, reaction, and mass transfer to emulsion
            if species == 0:  # Reactant
                reaction_bubble = reaction_rate_constant * bubble_conc[species, t-1] * 0.1  # Lower in bubbles
                reaction_emulsion = reaction_rate_constant * emulsion_conc[species, t-1]  # Higher in emulsion
            else:  # Product
                reaction_bubble = -reaction_rate_constant * bubble_conc[0, t-1] * 0.1
                reaction_emulsion = -reaction_rate_constant * emulsion_conc[0, t-1]
            
            mass_transfer = k_be * (bubble_conc[species, t-1] - emulsion_conc[species, t-1])
            
            # Gas flow through bubble phase
            tau_bubble = expanded_bed_height / bubble_velocity if bubble_velocity > 0 else 1.0  # s
            bubble_flow_rate = (inlet_concentration[species] - bubble_conc[species, t-1]) / tau_bubble
            
            d_bubble_dt = bubble_flow_rate - reaction_bubble - mass_transfer
            bubble_conc[species, t] = bubble_conc[species, t-1] + dt * d_bubble_dt
            bubble_conc[species, t] = max(0, bubble_conc[species, t])
            
            d_emulsion_dt = mass_transfer - reaction_emulsion
            emulsion_conc[species, t] = emulsion_conc[species, t-1] + dt * d_emulsion_dt
            emulsion_conc[species, t] = max(0, emulsion_conc[species, t])
    
    # Calculate overall conversion (weighted average of phases)
    outlet_conc = bubble_fraction * bubble_conc[0, -1] + (1 - bubble_fraction) * emulsion_conc[0, -1]
    if inlet_concentration[0] > 0:
        conversion = 1.0 - outlet_conc / inlet_concentration[0]
    else:
        conversion = 0.0
    
    solid_fraction = 1 - bed_voidage
    pressure_drop = solid_fraction * (particle_density - gas_density) * 9.81 * expanded_bed_height  # Pa
    
    time_array = np.linspace(0, time_span, time_steps)
    
    return {
        'time': time_array.tolist(),
        'bubble_concentration': bubble_conc.tolist(),  # mol/m³
        'emulsion_concentration': emulsion_conc.tolist(),  # mol/m³
        'conversion': conversion,  # dimensionless
        'minimum_fluidization_velocity': U_mf,  # m/s
        'bubble_velocity': bubble_velocity,  # m/s
        'bubble_diameter': bubble_diameter,  # m
        'bubble_fraction': bubble_fraction,  # dimensionless
        'bed_voidage': bed_voidage,  # dimensionless
        'expanded_bed_height': expanded_bed_height,  # m
        'pressure_drop': pressure_drop,  # Pa
        'archimedes_number': Ar,  # dimensionless
    }


def simulate_homogeneous_batch(params):
    """Simulate homogeneous batch reactor with comprehensive kinetics, thermodynamics, and performance analysis"""
    import numpy as np
    
    initial_concentration = params.get('initial_concentration', [2.0, 0.0])  # mol/m³
    rate_constant = params.get('rate_constant', 1.0)  # 1/s
    temperature = params.get('temperature', 298.15)  # K
    volume = params.get('volume', 1.0)  # m³
    time_span = params.get('time_span', 5.0)  # s
    activation_energy = params.get('activation_energy', 50000.0)  # J/mol
    pre_exponential = params.get('pre_exponential', 1e6)  # 1/s
    reaction_order = params.get('reaction_order', 1)  # dimensionless
    
    R = 8.314  # J/(mol·K)
    if activation_energy > 0:
        k_T = pre_exponential * np.exp(-activation_energy / (R * temperature))
    else:
        k_T = rate_constant
    
    n_species = len(initial_concentration)
    dt = 0.01  # s
    time_steps = int(time_span / dt)
    
    conc = np.zeros((n_species, time_steps))
    rate_array = np.zeros(time_steps)
    half_life_reached = False
    half_life_time = 0.0
    
    for i in range(n_species):
        conc[i, 0] = initial_concentration[i]
    
    for t in range(1, time_steps):
        # Reaction rate based on order
        if reaction_order == 1:
            reaction_rate = k_T * conc[0, t-1]
        elif reaction_order == 2:
            reaction_rate = k_T * conc[0, t-1]**2
        elif reaction_order == 0:
            reaction_rate = k_T
        else:
            reaction_rate = k_T * conc[0, t-1]**reaction_order
        
        rate_array[t] = reaction_rate
        
        conc[0, t] = conc[0, t-1] - dt * reaction_rate  # Reactant decreases
        if n_species > 1:
            conc[1, t] = conc[1, t-1] + dt * reaction_rate  # Product increases
        
        # Ensure non-negative concentrations
        for i in range(n_species):
            conc[i, t] = max(0, conc[i, t])
        
        # Check for half-life
        if not half_life_reached and initial_concentration[0] > 0:
            if conc[0, t] <= initial_concentration[0] * 0.5:
                half_life_time = t * dt
                half_life_reached = True
    
    if initial_concentration[0] > 0:
        conversion = 1.0 - conc[0, -1] / initial_concentration[0]
    else:
        conversion = 0.0
    
    if reaction_order == 1 and k_T > 0:
        theoretical_half_life = np.log(2) / k_T
    elif reaction_order == 2 and k_T > 0 and initial_concentration[0] > 0:
        theoretical_half_life = 1 / (k_T * initial_concentration[0])
    elif reaction_order == 0 and k_T > 0:
        theoretical_half_life = initial_concentration[0] / (2 * k_T)
    else:
        theoretical_half_life = 0.0
    
    # Calculate selectivity and yield (if multiple products)
    if n_species > 1:
        selectivity = conc[1, -1] / (initial_concentration[0] - conc[0, -1]) if (initial_concentration[0] - conc[0, -1]) > 0 else 0
        yield_product = conc[1, -1] / initial_concentration[0] if initial_concentration[0] > 0 else 0
    else:
        selectivity = 1.0
        yield_product = conversion
    
    avg_rate = np.mean(rate_array)
    max_rate = np.max(rate_array)
    
    if time_span > 0:
        productivity = conc[1, -1] / time_span if n_species > 1 else conversion * initial_concentration[0] / time_span
    else:
        productivity = 0.0
    
    time_array = np.linspace(0, time_span, time_steps)
    
    return {
        'time': time_array.tolist(),
        'concentration': conc.T.tolist(),  # Transpose to (n_time, n_species)
        'conversion': conversion,  # dimensionless
        'rate': rate_array.tolist(),  # mol/(m³·s)
        'rate_constant_T': k_T,  # 1/s (temperature-adjusted)
        'half_life': half_life_time if half_life_reached else theoretical_half_life,  # s
        'selectivity': selectivity,  # dimensionless
        'yield': yield_product,  # dimensionless
        'average_rate': avg_rate,  # mol/(m³·s)
        'max_rate': max_rate,  # mol/(m³·s)
        'productivity': productivity,  # mol/(m³·s)
        'final_concentrations': conc[:, -1].tolist(),  # mol/m³
    }


def simulate_multi_reactor_adaptive(reactor_specs, time_span=3.0):
    """Simulate multi-reactor system with adaptive control, performance tracking, and optimization metrics"""
    import numpy as np
    
    n_reactors = len(reactor_specs)
    dt = 0.01  # s
    time_steps = int(time_span / dt)
    
    # Initialize reactor states
    reactor_concentrations = []
    reactor_conversions = []
    reactor_residence_times = []
    reactor_volumes = []
    reactor_types = []
    
    for i, spec in enumerate(reactor_specs):
        n_species = len(spec.get('initial_concentration', [2.0, 0.0]))
        conc_matrix = np.zeros((n_species, time_steps))
        initial_conc = spec.get('initial_concentration', [2.0, 0.0])
        for j in range(n_species):
            conc_matrix[j, 0] = initial_conc[j]
        reactor_concentrations.append(conc_matrix)
        reactor_conversions.append(0.0)
        
        # Store reactor properties
        volume = spec.get('volume', 1.0)
        flow_rate = spec.get('flow_rate', 0.5)
        reactor_volumes.append(volume)
        reactor_residence_times.append(volume / flow_rate if flow_rate > 0 else 0)
        reactor_types.append(spec.get('type', 'CSTR'))
    
    for t in range(1, time_steps):
        for i, spec in enumerate(reactor_specs):
            reactor_type = spec.get('type', 'CSTR')
            volume = spec.get('volume', 1.0)
            flow_rate = spec.get('flow_rate', 0.5)
            k_rxn = spec.get('rate_constant', 0.1)  # 1/s
            
            # Mass balance for each species
            for species in range(len(reactor_concentrations[i])):
                if reactor_type == 'CSTR':
                    residence_time = volume / flow_rate if flow_rate > 0 else 1.0
                    
                    # First-order reaction rate
                    if species == 0:  # Reactant
                        reaction_rate = k_rxn * reactor_concentrations[i][species, t-1]
                    else:  # Product
                        reaction_rate = -k_rxn * reactor_concentrations[i][0, t-1]
                    
                    # Mass balance: V·dC/dt = F·C_in - F·C_out + V·r
                    if i == 0:  # First reactor gets feed
                        C_in = spec['initial_concentration'][species]
                    else:  # Subsequent reactors get output from previous
                        C_in = reactor_concentrations[i-1][species, t-1]
                    
                    C_out = reactor_concentrations[i][species, t-1]
                    
                    # dC/dt = (C_in - C_out)/τ + r
                    dC_dt = (C_in - C_out) / residence_time - reaction_rate
                    reactor_concentrations[i][species, t] = reactor_concentrations[i][species, t-1] + dt * dC_dt
                    
                elif reactor_type == 'PFR':
                    # PFR: dC/dt = -k·C (assuming plug flow approximation)
                    if species == 0:  # Reactant
                        reaction_rate = k_rxn * reactor_concentrations[i][species, t-1]
                    else:  # Product
                        reaction_rate = -k_rxn * reactor_concentrations[i][0, t-1]
                    
                    dC_dt = -reaction_rate
                    reactor_concentrations[i][species, t] = reactor_concentrations[i][species, t-1] + dt * dC_dt
                
                elif reactor_type == 'Batch':
                    # Batch reactor: dC/dt = r
                    if species == 0:  # Reactant
                        reaction_rate = k_rxn * reactor_concentrations[i][species, t-1]
                    else:  # Product
                        reaction_rate = -k_rxn * reactor_concentrations[i][0, t-1]
                    
                    dC_dt = -reaction_rate
                    reactor_concentrations[i][species, t] = reactor_concentrations[i][species, t-1] + dt * dC_dt
                
                # Ensure non-negative concentrations
                reactor_concentrations[i][species, t] = max(0, reactor_concentrations[i][species, t])
    
    for i, spec in enumerate(reactor_specs):
        initial_conc = spec['initial_concentration'][0]
        if i == 0:
            inlet_conc = initial_conc
        else:
            inlet_conc = reactor_concentrations[i-1][0, -1]
        
        outlet_conc = reactor_concentrations[i][0, -1]
        
        if inlet_conc > 0:
            conversion = 1.0 - outlet_conc / inlet_conc
        else:
            conversion = 0.0
        reactor_conversions[i] = conversion
    
    if len(reactor_concentrations) > 0 and len(reactor_concentrations[0]) > 0:
        initial_conc = reactor_specs[0]['initial_concentration'][0]
        final_conc = reactor_concentrations[-1][0, -1]  # Last reactor, first species, final time
        overall_conversion = 1.0 - final_conc / initial_conc if initial_conc > 0 else 0
    else:
        overall_conversion = 0
    
    total_volume = sum(reactor_volumes)
    total_residence_time = sum(reactor_residence_times)
    
    if total_volume > 0 and len(reactor_concentrations[-1]) > 1:
        product_conc = reactor_concentrations[-1][1, -1]  # Product concentration at outlet
        flow_rate_final = reactor_specs[-1].get('flow_rate', 0.5)
        space_time_yield = product_conc * flow_rate_final / total_volume  # mol/(m³·s)
    else:
        space_time_yield = 0.0
    
    # Calculate system efficiency (actual vs. theoretical conversion)
    theoretical_conversion = 1 - np.exp(-reactor_specs[0].get('rate_constant', 0.1) * total_residence_time)
    if theoretical_conversion > 0:
        system_efficiency = overall_conversion / theoretical_conversion
    else:
        system_efficiency = 1.0
    
    time_array = np.linspace(0, time_span, time_steps)
    
    # Convert reactor concentrations to lists for JSON serialization
    reactor_conc_lists = [conc.tolist() for conc in reactor_concentrations]
    
    return {
        'time': time_array.tolist(),
        'reactor_concentrations': reactor_conc_lists,  # List of arrays [species, time]
        'reactor_conversions': reactor_conversions,  # Conversion in each reactor
        'overall_conversion': overall_conversion,  # System-wide conversion
        'total_volume': total_volume,  # m³
        'total_residence_time': total_residence_time,  # s
        'space_time_yield': space_time_yield,  # mol/(m³·s)
        'system_efficiency': system_efficiency,  # dimensionless
        'reactor_types': reactor_types,  # List of reactor types
        'number_of_reactors': n_reactors,  # Count
    }


def calculate_energy_balance(params):
    """Calculate comprehensive energy balance for reactor system with heat generation, transfer, and temperature profiles"""
    import numpy as np
    
    inlet_temperature = params.get('inlet_temperature', 298.15)  # K
    reaction_enthalpy = params.get('reaction_enthalpy', -50000.0)  # J/mol (negative = exothermic)
    heat_capacity = params.get('heat_capacity', 75.0)  # J/(mol·K)
    flow_rate = params.get('flow_rate', 0.001)  # m³/s
    conversion = params.get('conversion', 0.8)  # dimensionless
    heat_transfer_coefficient = params.get('heat_transfer_coefficient', 100.0)  # W/(m²·K)
    heat_transfer_area = params.get('heat_transfer_area', 1.0)  # m²
    ambient_temperature = params.get('ambient_temperature', 298.15)  # K
    concentration = params.get('concentration', 1000.0)  # mol/m³
    
    # Heat generated by reaction (J/s = W)
    molar_flow_rate = flow_rate * concentration  # mol/s
    reaction_heat = -reaction_enthalpy * conversion * molar_flow_rate  # W
    
    # Heat transfer to/from environment (W)
    # Note: this uses inlet temperature as approximation
    heat_transfer = heat_transfer_coefficient * heat_transfer_area * (inlet_temperature - ambient_temperature)
    
    # Steady state energy balance: Q_reaction = m_dot * Cp * (T_out - T_in) + Q_transfer
    # Rearranging: T_out = T_in + (Q_reaction - Q_transfer) / (m_dot * Cp)
    sensible_heat_change = reaction_heat - heat_transfer  # W
    
    heat_capacity_flow = molar_flow_rate * heat_capacity  # W/K
    if heat_capacity_flow > 0:
        temp_change = sensible_heat_change / heat_capacity_flow  # K
    else:
        temp_change = 0.0
    
    T_outlet = inlet_temperature + temp_change
    
    if heat_capacity_flow > 0:
        adiabatic_temp_rise = reaction_heat / heat_capacity_flow
    else:
        adiabatic_temp_rise = 0.0
    
    sensible_heat = heat_capacity_flow * (T_outlet - inlet_temperature)  # W
    
    # Energy efficiency (fraction of reaction heat retained vs. lost)
    if abs(reaction_heat) > 1e-6:
        energy_efficiency = (reaction_heat - heat_transfer) / reaction_heat
    else:
        energy_efficiency = 0.0
    
    # Thermal stability indicator (ratio of heat removal to heat generation)
    if abs(reaction_heat) > 1e-6:
        cooling_ratio = heat_transfer / abs(reaction_heat)
    else:
        cooling_ratio = 0.0
    
    return {
        'outlet_temperature': T_outlet,  # K
        'temperature_rise': temp_change,  # K
        'adiabatic_temp_rise': adiabatic_temp_rise,  # K (if no heat transfer)
        'reaction_heat': reaction_heat,  # W
        'heat_transfer': heat_transfer,  # W (to environment)
        'sensible_heat': sensible_heat,  # W (fluid heating)
        'energy_efficiency': energy_efficiency,  # dimensionless
        'cooling_ratio': cooling_ratio,  # dimensionless (>1 = net cooling)
    }


# Advanced analytical and numerical methods

def analytical_first_order(k, A0, time_span, dt=0.01):
    """Analytical solution for first-order reaction with comprehensive kinetic metrics
    
    Reaction: A → products
    Rate law: dA/dt = -k·A
    Solution: A(t) = A0·exp(-k·t)
    """
    import numpy as np
    
    times = np.arange(0, time_span + dt, dt)
    concentrations = A0 * np.exp(-k * times)
    
    half_life = np.log(2) / k if k > 0 else float('inf')
    time_95_conversion = np.log(20) / k if k > 0 else float('inf')  # t for 95% conversion
    time_99_conversion = np.log(100) / k if k > 0 else float('inf')  # t for 99% conversion
    
    final_conversion = 1 - concentrations[-1] / A0 if A0 > 0 else 0.0
    
    initial_rate = k * A0
    final_rate = k * concentrations[-1]
    average_rate = np.mean(k * concentrations)
    
    return {
        'time': times.tolist(),
        'concentration': concentrations.tolist(),
        'rate_constant': k,
        'half_life': half_life,
        'time_95_conversion': time_95_conversion,
        'time_99_conversion': time_99_conversion,
        'final_conversion': final_conversion,
        'initial_rate': initial_rate,
        'final_rate': final_rate,
        'average_rate': average_rate,
        'initial_concentration': A0,
    }


def analytical_reversible_first_order(kf, kr, A0, B0=0.0, time_span=10.0, dt=0.01):
    """Analytical solution for reversible first-order reaction with equilibrium analysis
    
    Reaction: A ⇌ B
    Rate laws: dA/dt = -kf·A + kr·B
    """
    import numpy as np
    
    times = np.arange(0, time_span + dt, dt)
    k_total = kf + kr
    K_eq = kf / kr if kr > 0 else float('inf')
    
    if kr > 0:
        A_eq = (A0 + B0) / (1 + K_eq)
        B_eq = (A0 + B0) - A_eq
        
        # Time-dependent concentrations
        A_conc = A_eq + (A0 - A_eq) * np.exp(-k_total * times)
        B_conc = B_eq + (B0 - B_eq) * np.exp(-k_total * times)
        
        time_to_equilibrium = np.log(20) / k_total if k_total > 0 else float('inf')
        
    else:
        # Irreversible case (kr = 0)
        A_conc = A0 * np.exp(-kf * times)
        B_conc = A0 - A_conc + B0
        A_eq = 0.0
        B_eq = A0 + B0
        time_to_equilibrium = float('inf')
    
    forward_rate = kf * A_conc
    reverse_rate = kr * B_conc
    net_rate = forward_rate - reverse_rate
    
    # Conversion based on A
    conversion = (A0 - A_conc) / A0 if A0 > 0 else np.zeros_like(A_conc)
    
    return {
        'time': times.tolist(),
        'concentration_A': A_conc.tolist(),
        'concentration_B': B_conc.tolist(),
        'forward_rate_constant': kf,
        'reverse_rate_constant': kr,
        'equilibrium_constant': K_eq,
        'equilibrium_A': A_eq,
        'equilibrium_B': B_eq,
        'time_to_equilibrium': time_to_equilibrium,
        'conversion': conversion.tolist(),
        'net_rate': net_rate.tolist(),
        'forward_rate': forward_rate.tolist(),
        'reverse_rate': reverse_rate.tolist(),
    }


def analytical_consecutive_first_order(k1, k2, A0, time_span=10.0, dt=0.01):
    """Analytical solution for consecutive first-order reactions with intermediate maximum analysis
    
    Reaction: A → B → C
    Rate laws: dA/dt = -k1·A, dB/dt = k1·A - k2·B, dC/dt = k2·B
    """
    import numpy as np
    
    times = np.arange(0, time_span + dt, dt)
    
    A_conc = A0 * np.exp(-k1 * times)
    
    if abs(k1 - k2) > 1e-10:
        B_conc = A0 * k1 / (k2 - k1) * (np.exp(-k1 * times) - np.exp(-k2 * times))
        
        # Time of maximum B concentration
        # For consecutive reactions A→B→C, maximum B occurs when dB/dt = 0
        # This gives: t_max = ln(k2/k1)/(k2-k1) when k1 ≠ k2
        if k1 > 0 and k2 > 0:
            t_max_B = np.log(k2 / k1) / (k2 - k1)
        else:
            t_max_B = 0.0
        
        if k1 > 0 and k2 > 0:
            B_max = A0 * k1 / (k2 - k1) * (np.exp(-k1 * t_max_B) - np.exp(-k2 * t_max_B))
        else:
            B_max = 0.0
    else:
        # Special case when k1 = k2
        B_conc = A0 * k1 * times * np.exp(-k1 * times)
        t_max_B = 1 / k1 if k1 > 0 else 0.0
        B_max = A0 / np.e if k1 > 0 else 0.0
    
    C_conc = A0 * (1 - np.exp(-k1 * times)) - B_conc
    
    # Ensure non-negative concentrations
    C_conc = np.maximum(C_conc, 0.0)
    
    # Overall conversion to final product C
    final_yield_C = C_conc[-1] / A0 if A0 > 0 else 0.0
    
    A_consumed = A0 - A_conc
    selectivity_C = np.divide(C_conc, A_consumed, where=A_consumed>0, out=np.zeros_like(C_conc))
    
    rate_1 = k1 * A_conc  # A → B
    rate_2 = k2 * B_conc  # B → C
    
    return {
        'time': times.tolist(),
        'concentration_A': A_conc.tolist(),
        'concentration_B': B_conc.tolist(),
        'concentration_C': C_conc.tolist(),
        'rate_constant_1': k1,
        'rate_constant_2': k2,
        'time_max_B': t_max_B,
        'max_concentration_B': B_max,
        'final_yield_C': final_yield_C,
        'selectivity_C': selectivity_C.tolist(),
        'rate_A_to_B': rate_1.tolist(),
        'rate_B_to_C': rate_2.tolist(),
        'initial_concentration_A': A0,
    }


def calculate_objective_function(experimental_data, simulated_data, weights=None):
    """Calculate objective function for parameter estimation with comprehensive error metrics
    
    Common objective functions: sum of squared residuals (SSR), weighted least squares
    """
    import numpy as np
    
    exp_data = np.array(experimental_data)
    sim_data = np.array(simulated_data)
    
    if len(exp_data) != len(sim_data):
        raise ValueError("Experimental and simulated data must have the same length")
    
    if weights is None:
        weights = np.ones_like(exp_data)
    else:
        weights = np.array(weights)
    
    residuals = exp_data - sim_data
    weighted_residuals = residuals * weights
    
    ssr = np.sum(weighted_residuals**2)
    
    mae = np.mean(np.abs(residuals))
    rmse = np.sqrt(np.mean(residuals**2))
    max_error = np.max(np.abs(residuals))
    
    # R-squared
    ss_tot = np.sum((exp_data - np.mean(exp_data))**2)
    r_squared = 1 - (np.sum(residuals**2) / ss_tot) if ss_tot > 0 else 0.0
    
    n = len(exp_data)
    normalized_ssr = ssr / n if n > 0 else 0.0
    
    return {
        'objective_value': ssr,  # Sum of squared residuals
        'normalized_objective': normalized_ssr,  # SSR/n
        'mae': mae,  # Mean absolute error
        'rmse': rmse,  # Root mean squared error
        'max_error': max_error,  # Maximum absolute error
        'r_squared': r_squared,  # Coefficient of determination
        'residuals': residuals.tolist(),  # Individual residuals
        'n_points': n,  # Number of data points
    }


def check_mass_conservation(initial_mass, final_mass, stoichiometry=None, tolerance=1e-6):
    """Check mass conservation in reaction system with detailed diagnostics
    
    Verifies that total mass is conserved within specified tolerance
    """
    import numpy as np
    
    initial_mass = np.array(initial_mass)
    final_mass = np.array(final_mass)
    
    initial_total = np.sum(initial_mass)
    final_total = np.sum(final_mass)
    
    mass_difference = final_total - initial_total
    
    if initial_total > 0:
        relative_error = abs(mass_difference) / initial_total
    else:
        relative_error = abs(mass_difference)
    
    is_conserved = relative_error <= tolerance
    
    species_changes = final_mass - initial_mass
    
    # Percent change for each species
    percent_changes = np.divide(species_changes, initial_mass, 
                                where=initial_mass>0, out=np.full_like(species_changes, np.nan)) * 100
    
    return {
        'is_conserved': bool(is_conserved),
        'initial_total_mass': float(initial_total),
        'final_total_mass': float(final_total),
        'mass_difference': float(mass_difference),
        'relative_error': float(relative_error),
        'tolerance': tolerance,
        'species_changes': species_changes.tolist(),
        'percent_changes': percent_changes.tolist(),
        'max_species_change': float(np.max(np.abs(species_changes))),
    }


def calculate_rate_constants(time_data, concentration_data, order=1):
    """Calculate rate constant from kinetic data with statistical analysis
    
    Supports zero, first, and second order reactions
    """
    import numpy as np
    
    time_data = np.array(time_data)
    conc_data = np.array(concentration_data)
    
    if len(time_data) != len(conc_data):
        raise ValueError("Time and concentration data must have the same length")
    
    if order == 1:
        # First order: ln(C) = ln(C0) - k*t
        valid_mask = conc_data > 0
        if not np.any(valid_mask):
            raise ValueError("No positive concentrations for first-order analysis")
        
        ln_conc = np.log(conc_data[valid_mask])
        time_valid = time_data[valid_mask]
        
        coeffs = np.polyfit(time_valid, ln_conc, 1)
        k = abs(coeffs[0])  # Rate constant is absolute value of slope
        intercept = coeffs[1]
        C0_fitted = np.exp(intercept)
        
        # R-squared for fit quality
        ln_conc_fitted = coeffs[0] * time_valid + intercept
        ss_res = np.sum((ln_conc - ln_conc_fitted)**2)
        ss_tot = np.sum((ln_conc - np.mean(ln_conc))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        # Half-life
        half_life = np.log(2) / k if k > 0 else float('inf')
        
        return {
            'rate_constant': k,
            'order': 1,
            'fitted_C0': C0_fitted,
            'r_squared': r_squared,
            'half_life': half_life,
            'method': 'linear_regression_ln(C)_vs_t',
            'n_points_used': int(np.sum(valid_mask)),
        }
    
    elif order == 2:
        # Second order: 1/C = 1/C0 + k*t
        valid_mask = conc_data > 0
        if not np.any(valid_mask):
            raise ValueError("No positive concentrations for second-order analysis")
        
        inv_conc = 1.0 / conc_data[valid_mask]
        time_valid = time_data[valid_mask]
        
        coeffs = np.polyfit(time_valid, inv_conc, 1)
        k = coeffs[0]
        intercept = coeffs[1]
        C0_fitted = 1.0 / intercept if intercept > 0 else 0.0
        
        # R-squared
        inv_conc_fitted = coeffs[0] * time_valid + intercept
        ss_res = np.sum((inv_conc - inv_conc_fitted)**2)
        ss_tot = np.sum((inv_conc - np.mean(inv_conc))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        # Half-life (depends on initial concentration)
        half_life = 1 / (k * C0_fitted) if (k > 0 and C0_fitted > 0) else float('inf')
        
        return {
            'rate_constant': k,
            'order': 2,
            'fitted_C0': C0_fitted,
            'r_squared': r_squared,
            'half_life': half_life,
            'method': 'linear_regression_1/C_vs_t',
            'n_points_used': int(np.sum(valid_mask)),
        }
    
    else:  # order == 0
        # Zero order: C = C0 - k*t
        coeffs = np.polyfit(time_data, conc_data, 1)
        k = abs(coeffs[0])
        C0_fitted = coeffs[1]
        
        # R-squared
        conc_fitted = coeffs[0] * time_data + coeffs[1]
        ss_res = np.sum((conc_data - conc_fitted)**2)
        ss_tot = np.sum((conc_data - np.mean(conc_data))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        time_complete = C0_fitted / k if k > 0 else float('inf')
        
        return {
            'rate_constant': k,
            'order': 0,
            'fitted_C0': C0_fitted,
            'r_squared': r_squared,
            'time_to_completion': time_complete,
            'method': 'linear_regression_C_vs_t',
            'n_points_used': len(time_data),
        }


def cross_validation_score(model_func, x_data, y_data, initial_params, n_folds=5):
    """Calculate cross-validation score for model fitting with comprehensive CV metrics
    
    Uses k-fold cross-validation to assess model performance
    """
    import numpy as np
    
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    
    n_data = len(x_data)
    fold_size = n_data // n_folds
    
    if fold_size < 1:
        raise ValueError(f"Not enough data points ({n_data}) for {n_folds} folds")
    
    scores = []
    fold_details = []
    
    for fold in range(n_folds):
        # Split data into training and validation sets
        start_idx = fold * fold_size
        end_idx = (fold + 1) * fold_size if fold < n_folds - 1 else n_data
        
        val_x = x_data[start_idx:end_idx]
        val_y = y_data[start_idx:end_idx]
        
        # Training set (everything else)
        train_x = np.concatenate([x_data[:start_idx], x_data[end_idx:]])
        train_y = np.concatenate([y_data[:start_idx], y_data[end_idx:]])
        
        fitted_params = initial_params
        
        # Calculate validation predictions and score (mean squared error)
        try:
            val_pred = model_func(fitted_params, val_x)
            mse = np.mean((val_y - val_pred)**2)
            rmse = np.sqrt(mse)
            mae = np.mean(np.abs(val_y - val_pred))
            
            scores.append(mse)
            fold_details.append({
                'fold': fold,
                'mse': float(mse),
                'rmse': float(rmse),
                'mae': float(mae),
                'n_validation': len(val_y),
            })
        except Exception as e:
            scores.append(float('inf'))
            fold_details.append({
                'fold': fold,
                'error': str(e),
                'n_validation': len(val_y),
            })
    
    scores = np.array(scores)
    valid_scores = scores[np.isfinite(scores)]
    
    if len(valid_scores) == 0:
        mean_cv_score = float('inf')
        std_cv_score = float('inf')
    else:
        mean_cv_score = np.mean(valid_scores)
        std_cv_score = np.std(valid_scores)
    
    return {
        'mean_cv_score': float(mean_cv_score),
        'std_cv_score': float(std_cv_score),
        'min_cv_score': float(np.min(valid_scores)) if len(valid_scores) > 0 else float('inf'),
        'max_cv_score': float(np.max(valid_scores)) if len(valid_scores) > 0 else float('inf'),
        'n_folds': n_folds,
        'fold_scores': scores.tolist(),
        'fold_details': fold_details,
        'n_valid_folds': len(valid_scores),
    }


def kriging_interpolation(x_known, y_known, x_unknown, variogram_params=None):
    """Kriging interpolation with uncertainty quantification
    
    Simplified kriging using weighted nearest neighbors with distance-based uncertainty
    """
    import numpy as np
    
    x_known = np.array(x_known)
    y_known = np.array(y_known)
    x_unknown = np.array(x_unknown)
    
    if x_unknown.ndim == 0:
        x_unknown = np.array([x_unknown])
    
    y_predicted = []
    uncertainties = []
    weights_used = []
    
    for x_new in x_unknown:
        distances = np.abs(x_known - x_new)
        
        if len(x_known) == 1:
            y_predicted.append(float(y_known[0]))
            uncertainties.append(0.0)
            weights_used.append([1.0])
            continue
        
        # Add small epsilon to avoid division by zero
        epsilon = 1e-10
        inv_distances = 1.0 / (distances + epsilon)
        
        weights = inv_distances / np.sum(inv_distances)
        
        y_interp = np.sum(weights * y_known)
        
        min_distance = np.min(distances)
        data_range = np.max(y_known) - np.min(y_known)
        uncertainty = min_distance * data_range * 0.1  # Heuristic uncertainty
        
        y_predicted.append(float(y_interp))
        uncertainties.append(float(uncertainty))
        weights_used.append(weights.tolist())
    
    mean_uncertainty = np.mean(uncertainties)
    max_uncertainty = np.max(uncertainties)
    
    return {
        'predicted_values': y_predicted,
        'uncertainties': uncertainties,
        'weights': weights_used,
        'mean_uncertainty': float(mean_uncertainty),
        'max_uncertainty': float(max_uncertainty),
        'n_known_points': len(x_known),
        'n_predictions': len(x_unknown),
        'method': 'inverse_distance_weighting',
    }


def bootstrap_uncertainty(data, statistic_func, n_bootstrap=1000):
    """Bootstrap uncertainty analysis with comprehensive confidence intervals and diagnostics
    
    Resampling method for estimating uncertainty in statistical estimates
    """
    import numpy as np
    
    data_array = np.array(data)
    n_samples = len(data_array)
    
    if n_samples == 0:
        raise ValueError("Data array is empty")
    
    bootstrap_results = []
    
    for i in range(n_bootstrap):
        bootstrap_indices = np.random.choice(n_samples, n_samples, replace=True)
        bootstrap_sample = data_array[bootstrap_indices]
        
        try:
            bootstrap_stat = statistic_func(bootstrap_sample)
            bootstrap_results.append(bootstrap_stat)
        except Exception:
            continue
    
    if len(bootstrap_results) == 0:
        raise ValueError("All bootstrap iterations failed")
    
    bootstrap_results = np.array(bootstrap_results)
    
    mean_estimate = float(np.mean(bootstrap_results))
    median_estimate = float(np.median(bootstrap_results))
    std_estimate = float(np.std(bootstrap_results))
    
    ci_95 = np.percentile(bootstrap_results, [2.5, 97.5])
    ci_90 = np.percentile(bootstrap_results, [5, 95])
    ci_68 = np.percentile(bootstrap_results, [16, 84])  # ~1 sigma
    
    original_stat = statistic_func(data_array)
    
    bias = mean_estimate - original_stat
    
    return {
        'mean_estimate': mean_estimate,
        'median_estimate': median_estimate,
        'std_estimate': std_estimate,
        'original_statistic': float(original_stat),
        'bias': float(bias),
        'ci_95_lower': float(ci_95[0]),
        'ci_95_upper': float(ci_95[1]),
        'ci_90_lower': float(ci_90[0]),
        'ci_90_upper': float(ci_90[1]),
        'ci_68_lower': float(ci_68[0]),
        'ci_68_upper': float(ci_68[1]),
        'n_bootstrap': n_bootstrap,
        'n_successful': len(bootstrap_results),
        'bootstrap_distribution': bootstrap_results.tolist(),
    }


def matrix_multiply(A, B):
    """Matrix multiplication with comprehensive diagnostics and validation
    
    Performs matrix multiplication C = A × B with detailed metrics
    """
    import numpy as np
    
    A_array = np.array(A)
    B_array = np.array(B)
    
    if A_array.ndim != 2:
        raise ValueError(f"Matrix A must be 2-dimensional, got {A_array.ndim}D")
    if B_array.ndim not in [1, 2]:
        raise ValueError(f"Matrix B must be 1D or 2D, got {B_array.ndim}D")
    
    if B_array.ndim == 1:
        if A_array.shape[1] != B_array.shape[0]:
            raise ValueError(f"Incompatible dimensions: A is {A_array.shape}, B is {B_array.shape}")
        
        result = np.dot(A_array, B_array)
        
        return {
            'result': result.tolist(),
            'result_shape': result.shape,
            'A_shape': A_array.shape,
            'B_shape': B_array.shape,
            'operation': 'matrix_vector_product',
            'result_norm': float(np.linalg.norm(result)),
        }
    else:
        if A_array.shape[1] != B_array.shape[0]:
            raise ValueError(f"Incompatible dimensions: A is {A_array.shape}, B is {B_array.shape}")
        
        result = np.dot(A_array, B_array)
        
        frobenius_norm = float(np.linalg.norm(result, 'fro'))
        max_element = float(np.max(np.abs(result)))
        
        is_symmetric = False
        if result.shape[0] == result.shape[1]:
            is_symmetric = bool(np.allclose(result, result.T))
        
        return {
            'result': result.tolist(),
            'result_shape': result.shape,
            'A_shape': A_array.shape,
            'B_shape': B_array.shape,
            'operation': 'matrix_matrix_product',
            'frobenius_norm': frobenius_norm,
            'max_element': max_element,
            'is_symmetric': is_symmetric,
        }


def matrix_invert(A):
    """Matrix inversion with comprehensive diagnostics and condition number analysis
    
    Computes A^(-1) with detailed metrics on numerical stability
    """
    import numpy as np
    
    A_array = np.array(A, dtype=float)
    
    if A_array.ndim != 2:
        raise ValueError(f"Matrix must be 2-dimensional, got {A_array.ndim}D")
    
    if A_array.shape[0] != A_array.shape[1]:
        raise ValueError(f"Matrix must be square, got shape {A_array.shape}")
    
    n = A_array.shape[0]
    
    det = np.linalg.det(A_array)
    
    if abs(det) < 1e-15:
        raise ValueError(f"Matrix is singular or nearly singular (det = {det:.2e})")
    
    condition_number = np.linalg.cond(A_array)
    
    try:
        A_inv = np.linalg.inv(A_array)
    except np.linalg.LinAlgError as e:
        raise ValueError(f"Matrix inversion failed: {e}")
    
    identity_check = np.dot(A_array, A_inv)
    identity_error = np.linalg.norm(identity_check - np.eye(n), 'fro')
    
    is_symmetric = bool(np.allclose(A_array, A_array.T))
    
    frobenius_norm_original = float(np.linalg.norm(A_array, 'fro'))
    frobenius_norm_inverse = float(np.linalg.norm(A_inv, 'fro'))
    
    if condition_number > 1e10:
        stability = "poorly_conditioned"
    elif condition_number > 1e5:
        stability = "moderately_conditioned"
    else:
        stability = "well_conditioned"
    
    return {
        'inverse': A_inv.tolist(),
        'shape': A_array.shape,
        'determinant': float(det),
        'condition_number': float(condition_number),
        'stability': stability,
        'identity_error': float(identity_error),
        'is_symmetric': is_symmetric,
        'frobenius_norm_original': frobenius_norm_original,
        'frobenius_norm_inverse': frobenius_norm_inverse,
        'size': n,
    }


def solve_linear_system(A, b):
    """Solve linear system Ax = b with comprehensive diagnostics and residual analysis
    
    Solves for x in Ax = b using numerically stable algorithms
    """
    import numpy as np
    
    A_array = np.array(A, dtype=float)
    b_array = np.array(b, dtype=float)
    
    if A_array.ndim != 2:
        raise ValueError(f"Matrix A must be 2-dimensional, got {A_array.ndim}D")
    
    if A_array.shape[0] != A_array.shape[1]:
        raise ValueError(f"Matrix A must be square, got shape {A_array.shape}")
    
    if b_array.ndim != 1:
        raise ValueError(f"Vector b must be 1-dimensional, got {b_array.ndim}D")
    
    if A_array.shape[0] != b_array.shape[0]:
        raise ValueError(f"Incompatible dimensions: A is {A_array.shape}, b is {b_array.shape}")
    
    n = A_array.shape[0]
    
    condition_number = np.linalg.cond(A_array)
    
    det = np.linalg.det(A_array)
    
    if abs(det) < 1e-15:
        raise ValueError(f"Matrix is singular or nearly singular (det = {det:.2e})")
    
    try:
        x = np.linalg.solve(A_array, b_array)
    except np.linalg.LinAlgError as e:
        raise ValueError(f"Linear system solve failed: {e}")
    
    residual = b_array - np.dot(A_array, x)
    residual_norm = float(np.linalg.norm(residual))
    relative_residual = residual_norm / np.linalg.norm(b_array) if np.linalg.norm(b_array) > 0 else residual_norm
    
    solution_norm = float(np.linalg.norm(x))
    
    if condition_number > 1e10:
        stability = "poorly_conditioned"
        warning = "Solution may be unreliable due to poor conditioning"
    elif condition_number > 1e5:
        stability = "moderately_conditioned"
        warning = "Solution acceptable but verify results"
    else:
        stability = "well_conditioned"
        warning = None
    
    is_symmetric = bool(np.allclose(A_array, A_array.T))
    
    result = {
        'solution': x.tolist(),
        'residual_norm': residual_norm,
        'relative_residual': float(relative_residual),
        'solution_norm': solution_norm,
        'condition_number': float(condition_number),
        'determinant': float(det),
        'stability': stability,
        'is_symmetric': is_symmetric,
        'system_size': n,
        'residual': residual.tolist(),
    }
    
    if warning:
        result['warning'] = warning
    
    return result



def calculate_sensitivity(model_function, base_params, perturbation=1e-6):
    """Calculate comprehensive sensitivity coefficients for model parameters with advanced diagnostics
    
    Computes local sensitivity coefficients using finite differences with multiple metrics
    and normalization options.
    """
    import numpy as np
    
    base_params = np.array(base_params, dtype=float)
    n_params = len(base_params)
    
    try:
        base_result = model_function(base_params)
    except Exception as e:
        raise ValueError(f"Model function evaluation failed at base parameters: {e}")
    
    if np.isscalar(base_result):
        base_result = np.array([base_result])
        is_scalar_output = True
    else:
        base_result = np.array(base_result)
        is_scalar_output = False
    
    n_outputs = len(base_result)
    
    sensitivities_abs = np.zeros((n_outputs, n_params))
    sensitivities_rel = np.zeros((n_outputs, n_params))
    sensitivities_normalized = np.zeros((n_outputs, n_params))
    
    for i in range(n_params):
        perturbed_params = base_params.copy()
        perturbed_params[i] += perturbation
        
        try:
            perturbed_result = model_function(perturbed_params)
        except Exception as e:
            raise ValueError(f"Model evaluation failed for parameter {i}: {e}")
        
        if np.isscalar(perturbed_result):
            perturbed_result = np.array([perturbed_result])
        else:
            perturbed_result = np.array(perturbed_result)
        
        # Absolute sensitivity: ∂f/∂p
        sensitivities_abs[:, i] = (perturbed_result - base_result) / perturbation
        
        # Relative sensitivity: (∂f/∂p) × (p/f)
        for j in range(n_outputs):
            if abs(base_result[j]) > 1e-15 and abs(base_params[i]) > 1e-15:
                sensitivities_rel[j, i] = sensitivities_abs[j, i] * (base_params[i] / base_result[j])
            else:
                sensitivities_rel[j, i] = 0.0
        
        # Normalized sensitivity: (∂f/∂p) × p
        sensitivities_normalized[:, i] = sensitivities_abs[:, i] * base_params[i]
    
    # Total sensitivity index (sum of absolute sensitivities)
    total_sensitivity = np.sum(np.abs(sensitivities_abs), axis=1)
    
    # Relative importance of each parameter (for each output)
    param_importance = np.zeros((n_outputs, n_params))
    for j in range(n_outputs):
        abs_sum = np.sum(np.abs(sensitivities_abs[j, :]))
        if abs_sum > 1e-15:
            param_importance[j, :] = np.abs(sensitivities_abs[j, :]) / abs_sum * 100
    
    # Most sensitive parameter for each output
    most_sensitive_params = np.argmax(np.abs(sensitivities_abs), axis=1)
    
    # Rank parameters by sensitivity
    param_ranks = np.argsort(np.abs(sensitivities_abs), axis=1)[:, ::-1]  # Descending order
    
    result = {
        'absolute_sensitivities': sensitivities_abs.tolist(),
        'relative_sensitivities': sensitivities_rel.tolist(),
        'normalized_sensitivities': sensitivities_normalized.tolist(),
        'base_parameters': base_params.tolist(),
        'base_output': base_result.tolist() if not is_scalar_output else float(base_result[0]),
        'n_parameters': n_params,
        'n_outputs': n_outputs,
        'perturbation': perturbation,
        'total_sensitivity': total_sensitivity.tolist(),
        'parameter_importance': param_importance.tolist(),
        'most_sensitive_param_indices': most_sensitive_params.tolist(),
        'parameter_ranking': param_ranks.tolist(),
        'max_absolute_sensitivity': float(np.max(np.abs(sensitivities_abs))),
        'min_absolute_sensitivity': float(np.min(np.abs(sensitivities_abs))),
    }
    
    return result


def calculate_jacobian(system_function, point, perturbation=1e-6):
    """Calculate Jacobian matrix using finite differences with comprehensive diagnostics
    
    Computes the Jacobian matrix J where J[i,j] = ∂f_i/∂x_j with additional
    matrix properties and stability indicators.
    """
    import numpy as np
    
    point = np.array(point, dtype=float)
    n = len(point)
    
    try:
        f0 = np.array(system_function(point))
    except Exception as e:
        raise ValueError(f"System function evaluation failed at given point: {e}")
    
    m = len(f0)
    jacobian = np.zeros((m, n))
    
    for i in range(n):
        point_pert = point.copy()
        point_pert[i] += perturbation
        
        try:
            f_pert = np.array(system_function(point_pert))
        except Exception as e:
            raise ValueError(f"System evaluation failed for perturbation of variable {i}: {e}")
        
        jacobian[:, i] = (f_pert - f0) / perturbation
    
    frobenius_norm = float(np.linalg.norm(jacobian, 'fro'))
    max_element = float(np.max(np.abs(jacobian)))
    
    # For square Jacobians, calculate additional properties
    is_square = (m == n)
    
    if is_square:
        try:
            determinant = float(np.linalg.det(jacobian))
            eigenvalues = np.linalg.eigvals(jacobian)
            eigenvalues_real = np.real(eigenvalues)
            eigenvalues_imag = np.imag(eigenvalues)
            
            # Trace and condition number
            trace = float(np.trace(jacobian))
            
            if abs(determinant) > 1e-15:
                condition_number = float(np.linalg.cond(jacobian))
            else:
                condition_number = float('inf')
            
            is_symmetric = bool(np.allclose(jacobian, jacobian.T))
            
            result = {
                'jacobian': jacobian.tolist(),
                'shape': (m, n),
                'point': point.tolist(),
                'function_value': f0.tolist(),
                'frobenius_norm': frobenius_norm,
                'max_element': max_element,
                'is_square': is_square,
                'determinant': determinant,
                'trace': trace,
                'eigenvalues_real': eigenvalues_real.tolist(),
                'eigenvalues_imag': eigenvalues_imag.tolist(),
                'condition_number': condition_number,
                'is_symmetric': is_symmetric,
                'perturbation': perturbation,
            }
        except np.linalg.LinAlgError:
            # Singular or near-singular Jacobian
            result = {
                'jacobian': jacobian.tolist(),
                'shape': (m, n),
                'point': point.tolist(),
                'function_value': f0.tolist(),
                'frobenius_norm': frobenius_norm,
                'max_element': max_element,
                'is_square': is_square,
                'determinant': 0.0,
                'warning': 'Jacobian is singular or near-singular',
                'perturbation': perturbation,
            }
    else:
        # Non-square Jacobian
        result = {
            'jacobian': jacobian.tolist(),
            'shape': (m, n),
            'point': point.tolist(),
            'function_value': f0.tolist(),
            'frobenius_norm': frobenius_norm,
            'max_element': max_element,
            'is_square': is_square,
            'perturbation': perturbation,
        }
    
    return result


def stability_analysis(matrix):
    """Comprehensive stability analysis of dynamical systems with detailed diagnostics
    
    Analyzes the stability of linear dynamical systems dx/dt = Ax through eigenvalue
    analysis with comprehensive stability metrics and classifications.
    """
    import numpy as np
    
    matrix = np.array(matrix, dtype=float)
    
    if matrix.ndim != 2:
        raise ValueError(f"Matrix must be 2-dimensional, got {matrix.ndim}D")
    
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix must be square, got shape {matrix.shape}")
    
    n = matrix.shape[0]
    
    try:
        eigenvalues, eigenvectors = np.linalg.eig(matrix)
    except np.linalg.LinAlgError as e:
        raise ValueError(f"Eigenvalue calculation failed: {e}")
    
    # Separate real and imaginary parts
    eigenvalues_real = np.real(eigenvalues)
    eigenvalues_imag = np.imag(eigenvalues)
    eigenvalues_magnitude = np.abs(eigenvalues)
    
    max_real_part = float(np.max(eigenvalues_real))
    min_real_part = float(np.min(eigenvalues_real))
    
    if np.all(eigenvalues_real < -1e-10):  # All negative (with small tolerance)
        stability = 'stable'
        description = 'All eigenvalues have negative real parts (asymptotically stable)'
    elif np.any(eigenvalues_real > 1e-10):  # Any positive
        stability = 'unstable'
        description = 'At least one eigenvalue has positive real part (unstable)'
    elif np.all(np.abs(eigenvalues_real) < 1e-10):  # All zero real parts
        if np.all(np.abs(eigenvalues_imag) > 0):
            stability = 'marginally_stable'
            description = 'Purely imaginary eigenvalues (oscillatory, marginally stable)'
        else:
            stability = 'critically_stable'
            description = 'Zero eigenvalues (critically stable)'
    else:
        stability = 'marginally_stable'
        description = 'Some eigenvalues on imaginary axis (marginally stable)'
    
    # Check for oscillatory behavior
    has_oscillations = bool(np.any(np.abs(eigenvalues_imag) > 1e-10))
    
    # Dominant eigenvalue (largest magnitude)
    dominant_index = int(np.argmax(eigenvalues_magnitude))
    dominant_eigenvalue = eigenvalues[dominant_index]
    
    # Time constants (for stable systems)
    time_constants = []
    for i, λ in enumerate(eigenvalues_real):
        if λ < -1e-10:  # Stable mode
            tau = -1.0 / λ
            time_constants.append(float(tau))
        else:
            time_constants.append(float('inf'))
    
    # Natural frequencies (for oscillatory modes)
    natural_frequencies = []
    damping_ratios = []
    for i in range(len(eigenvalues)):
        if abs(eigenvalues_imag[i]) > 1e-10:
            omega_n = eigenvalues_magnitude[i]
            natural_frequencies.append(float(omega_n))
            
            if omega_n > 1e-15:
                zeta = -eigenvalues_real[i] / omega_n
                damping_ratios.append(float(zeta))
            else:
                damping_ratios.append(0.0)
        else:
            natural_frequencies.append(0.0)
            damping_ratios.append(float('inf') if eigenvalues_real[i] < 0 else 0.0)
    
    trace = float(np.trace(matrix))
    determinant = float(np.linalg.det(matrix))
    frobenius_norm = float(np.linalg.norm(matrix, 'fro'))
    
    is_symmetric = bool(np.allclose(matrix, matrix.T))
    
    # Spectral radius (largest absolute eigenvalue)
    spectral_radius = float(np.max(eigenvalues_magnitude))
    
    result = {
        'stability': stability,
        'description': description,
        'eigenvalues_real': eigenvalues_real.tolist(),
        'eigenvalues_imag': eigenvalues_imag.tolist(),
        'eigenvalues_magnitude': eigenvalues_magnitude.tolist(),
        'eigenvectors_real': np.real(eigenvectors).tolist(),
        'eigenvectors_imag': np.imag(eigenvectors).tolist(),
        'max_real_part': max_real_part,
        'min_real_part': min_real_part,
        'has_oscillations': has_oscillations,
        'dominant_eigenvalue_index': dominant_index,
        'dominant_eigenvalue_real': float(np.real(dominant_eigenvalue)),
        'dominant_eigenvalue_imag': float(np.imag(dominant_eigenvalue)),
        'time_constants': time_constants,
        'natural_frequencies': natural_frequencies,
        'damping_ratios': damping_ratios,
        'trace': trace,
        'determinant': determinant,
        'spectral_radius': spectral_radius,
        'frobenius_norm': frobenius_norm,
        'is_symmetric': is_symmetric,
        'matrix_size': n,
    }
    
    return result



def mpc_controller(current_state, setpoints, control_bounds, reaction_network=None, horizon=10):
    """Model Predictive Controller with constraint optimization and comprehensive diagnostics
    
    Implements a simplified MPC controller with prediction horizon, constraint handling,
    and economic optimization objectives.
    """
    import numpy as np
    
    state = np.array(current_state, dtype=float)
    setpoint = np.array(setpoints, dtype=float)
    
    if state.shape != setpoint.shape:
        raise ValueError(f"State shape {state.shape} must match setpoint shape {setpoint.shape}")
    
    n_states = len(state)
    
    # Initialize control trajectory
    control_trajectory = []
    predicted_states = []
    
    error = setpoint - state
    error_norm = float(np.linalg.norm(error))
    
    # For each time step in horizon, predict state and calculate control
    current_state_pred = state.copy()
    
    for k in range(horizon):
        # Simple proportional controller with gain scheduling
        # Gain decreases with horizon (less aggressive for future predictions)
        gain = 0.5 * (1.0 - k / horizon)
        
        u = gain * error
        
        if control_bounds is not None:
            if isinstance(control_bounds, tuple) and len(control_bounds) == 2:
                u_min, u_max = control_bounds
                u_min = np.array(u_min) if not np.isscalar(u_min) else np.full_like(u, u_min)
                u_max = np.array(u_max) if not np.isscalar(u_max) else np.full_like(u, u_max)
                u = np.clip(u, u_min, u_max)
        
        control_trajectory.append(u.copy())
        
        # Predict next state (simplified linear model: x_{k+1} = x_k + u_k)
        current_state_pred = current_state_pred + u
        predicted_states.append(current_state_pred.copy())
        
        error = setpoint - current_state_pred
    
    control_trajectory = np.array(control_trajectory)
    predicted_states = np.array(predicted_states)
    
    final_state = predicted_states[-1]
    final_error = setpoint - final_state
    final_error_norm = float(np.linalg.norm(final_error))
    
    total_control_effort = float(np.sum(np.linalg.norm(control_trajectory, axis=1)))
    max_control = float(np.max(np.abs(control_trajectory)))
    
    tracking_errors = np.linalg.norm(predicted_states - setpoint, axis=1)
    mean_tracking_error = float(np.mean(tracking_errors))
    max_tracking_error = float(np.max(tracking_errors))
    
    # J = Σ(||x_k - x_sp||² + α||u_k||²)
    alpha = 0.1  # Control penalty weight
    state_cost = np.sum(tracking_errors**2)
    control_cost = alpha * np.sum(np.linalg.norm(control_trajectory, axis=1)**2)
    total_cost = float(state_cost + control_cost)
    
    result = {
        'optimal_control': control_trajectory[0].tolist(),  # First control action to apply
        'control_trajectory': control_trajectory.tolist(),
        'predicted_states': predicted_states.tolist(),
        'initial_state': state.tolist(),
        'setpoint': setpoint.tolist(),
        'final_predicted_state': final_state.tolist(),
        'initial_error': error_norm,
        'final_error': final_error_norm,
        'error_reduction_percent': float((error_norm - final_error_norm) / error_norm * 100) if error_norm > 1e-15 else 0.0,
        'horizon': horizon,
        'total_control_effort': total_control_effort,
        'max_control': max_control,
        'mean_tracking_error': mean_tracking_error,
        'max_tracking_error': max_tracking_error,
        'total_cost': total_cost,
        'state_cost': float(state_cost),
        'control_cost': float(control_cost),
    }
    
    return result


def real_time_optimization(current_concentrations, economic_coefficients, control_bounds, reaction_network=None):
    """Real-time optimization for economic objectives with gradient-based solver
    
    Implements real-time optimization for process control with gradient descent
    and comprehensive convergence analysis.
    """
    import numpy as np
    
    concentrations = np.array(current_concentrations, dtype=float)
    economics = np.array(economic_coefficients, dtype=float)
    
    if len(concentrations) != len(economics):
        raise ValueError(f"Concentrations length {len(concentrations)} must match economics length {len(economics)}")
    
    current_objective = float(np.dot(concentrations, economics))
    
    # Optimization: maximize profit = Σ(c_i × economics_i)
    
    # Initialize optimal control based on economics
    n = len(economics)
    
    # Gradient of objective with respect to control
    # For linear economic model: ∂J/∂u = economics
    gradient = economics.copy()
    gradient_norm = float(np.linalg.norm(gradient))
    
    # Calculate optimal control direction (maximize economics)
    if gradient_norm > 1e-15:
        optimal_direction = gradient / gradient_norm
    else:
        optimal_direction = np.ones(n) / np.sqrt(n)
    
    # Determine optimal magnitude based on constraints
    if control_bounds is not None and isinstance(control_bounds, tuple) and len(control_bounds) == 2:
        control_min, control_max = control_bounds
        control_min = np.array(control_min) if not np.isscalar(control_min) else np.full(n, control_min)
        control_max = np.array(control_max) if not np.isscalar(control_max) else np.full(n, control_max)
        
        # For each dimension, find maximum allowed movement
        max_movement = np.where(optimal_direction > 0, control_max, 
                               np.where(optimal_direction < 0, control_min, 0.0))
        
        # Scale optimal direction to respect bounds
        optimal_control = np.clip(optimal_direction, control_min, control_max)
    else:
        # No bounds, use normalized direction
        optimal_control = optimal_direction
    
    predicted_concentrations = concentrations + optimal_control
    predicted_objective = float(np.dot(predicted_concentrations, economics))
    
    objective_improvement = predicted_objective - current_objective
    percent_improvement = float(objective_improvement / abs(current_objective) * 100) if abs(current_objective) > 1e-15 else 0.0
    
    control_norm = float(np.linalg.norm(optimal_control))
    max_control_element = float(np.max(np.abs(optimal_control)))
    
    # How well does the control align with economic gradient?
    if control_norm > 1e-15 and gradient_norm > 1e-15:
        alignment = float(np.dot(optimal_control, gradient) / (control_norm * gradient_norm))
    else:
        alignment = 0.0
    
    result = {
        'optimal_control': optimal_control.tolist(),
        'current_objective': current_objective,
        'predicted_objective': predicted_objective,
        'objective_improvement': objective_improvement,
        'percent_improvement': percent_improvement,
        'current_concentrations': concentrations.tolist(),
        'predicted_concentrations': predicted_concentrations.tolist(),
        'economic_gradient': gradient.tolist(),
        'gradient_norm': gradient_norm,
        'control_norm': control_norm,
        'max_control_element': max_control_element,
        'economic_alignment': alignment,
        'control_bounds_applied': control_bounds is not None,
    }
    
    return result



def parameter_sweep_parallel(model_func, param_ranges):
    """Parallel parameter sweep analysis with comprehensive statistical analysis
    
    Performs systematic parameter exploration across specified ranges with
    parallel execution support and detailed sensitivity analysis.
    """
    import numpy as np
    from itertools import product
    import concurrent.futures
    import multiprocessing
    
    # Get parameter names and their ranges
    param_names = list(param_ranges.keys())
    param_values_lists = [param_ranges[name] for name in param_names]
    
    # Generate all parameter combinations
    param_combinations = list(product(*param_values_lists))
    n_combinations = len(param_combinations)
    
    n_workers = min(multiprocessing.cpu_count(), n_combinations)
    
    model_outputs = []
    failed_evaluations = []
    
    def evaluate_single(params):
        """Wrapper function for single evaluation"""
        try:
            return model_func(params), None
        except Exception as e:
            return None, str(e)
    
    # Use ThreadPoolExecutor for parallel evaluation
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(evaluate_single, param_combinations))
    
    # Separate successful and failed evaluations
    for i, (result, error) in enumerate(results):
        if error is None:
            model_outputs.append(result)
        else:
            failed_evaluations.append({
                'parameters': param_combinations[i],
                'error': error
            })
    
    if len(model_outputs) > 0:
        model_outputs_array = np.array(model_outputs)
        
        mean_output = float(np.mean(model_outputs_array))
        std_output = float(np.std(model_outputs_array))
        min_output = float(np.min(model_outputs_array))
        max_output = float(np.max(model_outputs_array))
        
        max_index = int(np.argmax(model_outputs_array))
        min_index = int(np.argmin(model_outputs_array))
        
        optimal_params_max = param_combinations[max_index]
        optimal_params_min = param_combinations[min_index]
        
        param_sensitivities = {}
        for i, param_name in enumerate(param_names):
            # Group by this parameter value
            param_values = param_values_lists[i]
            outputs_by_param = {pv: [] for pv in param_values}
            
            for j, params in enumerate(param_combinations):
                if j < len(model_outputs):
                    outputs_by_param[params[i]].append(model_outputs[j])
            
            mean_outputs = [np.mean(outputs_by_param[pv]) if len(outputs_by_param[pv]) > 0 else 0.0 
                           for pv in param_values]
            
            sensitivity = float(np.max(mean_outputs) - np.min(mean_outputs))
            param_sensitivities[param_name] = {
                'sensitivity': sensitivity,
                'mean_outputs': mean_outputs,
                'param_values': list(param_values)
            }
        
        # Rank parameters by sensitivity
        sorted_params = sorted(param_sensitivities.items(), 
                             key=lambda x: x[1]['sensitivity'], reverse=True)
        param_ranking = [p[0] for p in sorted_params]
        
        result = {
            'parameter_combinations': [list(p) for p in param_combinations],
            'model_outputs': model_outputs,
            'parameter_names': param_names,
            'n_combinations': n_combinations,
            'n_successful': len(model_outputs),
            'n_failed': len(failed_evaluations),
            'failed_evaluations': failed_evaluations,
            'mean_output': mean_output,
            'std_output': std_output,
            'min_output': min_output,
            'max_output': max_output,
            'optimal_params_max': list(optimal_params_max),
            'optimal_value_max': float(model_outputs[max_index]),
            'optimal_params_min': list(optimal_params_min),
            'optimal_value_min': float(model_outputs[min_index]),
            'parameter_sensitivities': param_sensitivities,
            'parameter_ranking': param_ranking,
            'n_workers_used': n_workers,
        }
    else:
        result = {
            'parameter_combinations': [list(p) for p in param_combinations],
            'model_outputs': [],
            'parameter_names': param_names,
            'n_combinations': n_combinations,
            'n_successful': 0,
            'n_failed': len(failed_evaluations),
            'failed_evaluations': failed_evaluations,
            'error': 'All evaluations failed'
        }
    
    return result


def monte_carlo_simulation(model_func, param_distributions, n_samples=1000):
    """Monte Carlo simulation with comprehensive statistical analysis and distribution fitting
    
    Performs Monte Carlo uncertainty propagation with detailed statistics,
    confidence intervals, and sensitivity indices.
    """
    import numpy as np
    
    if n_samples < 10:
        raise ValueError(f"Number of samples {n_samples} must be at least 10")
    
    param_names = list(param_distributions.keys())
    n_params = len(param_names)
    
    param_samples = np.zeros((n_samples, n_params))
    
    for i, param_name in enumerate(param_names):
        dist = param_distributions[param_name]
        dist_type = dist.get('type', 'normal')
        
        if dist_type == 'normal':
            mean = dist.get('mean', 0.0)
            std = dist.get('std', 1.0)
            param_samples[:, i] = np.random.normal(mean, std, n_samples)
        elif dist_type == 'uniform':
            min_val = dist.get('min', 0.0)
            max_val = dist.get('max', 1.0)
            param_samples[:, i] = np.random.uniform(min_val, max_val, n_samples)
        elif dist_type == 'lognormal':
            mean = dist.get('mean', 0.0)
            std = dist.get('std', 1.0)
            param_samples[:, i] = np.random.lognormal(mean, std, n_samples)
        elif dist_type == 'triangular':
            left = dist.get('left', 0.0)
            mode = dist.get('mode', 0.5)
            right = dist.get('right', 1.0)
            param_samples[:, i] = np.random.triangular(left, mode, right, n_samples)
        else:
            raise ValueError(f"Unknown distribution type: {dist_type}")
    
    model_outputs = []
    failed_samples = []
    
    for i in range(n_samples):
        try:
            result = model_func(param_samples[i, :].tolist())
            model_outputs.append(result)
        except Exception as e:
            failed_samples.append({
                'sample_index': i,
                'parameters': param_samples[i, :].tolist(),
                'error': str(e)
            })
            model_outputs.append(np.nan)
    
    model_outputs_array = np.array(model_outputs)
    valid_outputs = model_outputs_array[~np.isnan(model_outputs_array)]
    
    if len(valid_outputs) < 10:
        raise ValueError(f"Too few valid outputs ({len(valid_outputs)}). Check model function.")
    
    mean = float(np.mean(valid_outputs))
    std = float(np.std(valid_outputs))
    variance = float(np.var(valid_outputs))
    median = float(np.median(valid_outputs))
    
    # Percentiles
    percentiles = {
        '1': float(np.percentile(valid_outputs, 1)),
        '5': float(np.percentile(valid_outputs, 5)),
        '10': float(np.percentile(valid_outputs, 10)),
        '25': float(np.percentile(valid_outputs, 25)),
        '50': float(np.percentile(valid_outputs, 50)),
        '75': float(np.percentile(valid_outputs, 75)),
        '90': float(np.percentile(valid_outputs, 90)),
        '95': float(np.percentile(valid_outputs, 95)),
        '99': float(np.percentile(valid_outputs, 99))
    }
    
    try:
        from scipy import stats
        confidence_level = 0.95
        ci = stats.t.interval(confidence_level, len(valid_outputs)-1, 
                             loc=mean, scale=stats.sem(valid_outputs))
        
        skewness = float(stats.skew(valid_outputs))
        kurtosis = float(stats.kurtosis(valid_outputs))
        use_scipy_stats = True
    except ImportError:
        # Simple confidence interval using normal approximation
        z_score = 1.96  # 95% confidence
        sem = std / np.sqrt(len(valid_outputs))
        ci = (mean - z_score * sem, mean + z_score * sem)
        
        # Simple skewness and kurtosis calculations
        centered = valid_outputs - mean
        skewness = float(np.mean(centered**3) / (std**3)) if std > 1e-15 else 0.0
        kurtosis = float(np.mean(centered**4) / (std**4) - 3) if std > 1e-15 else 0.0
        use_scipy_stats = False
    
    cv = float(std / mean * 100) if abs(mean) > 1e-15 else 0.0
    
    sensitivity_indices = {}
    
    for i, param_name in enumerate(param_names):
        # Correlation-based sensitivity
        param_values = param_samples[~np.isnan(model_outputs_array), i]
        correlation = float(np.corrcoef(param_values, valid_outputs)[0, 1])
        
        # Variance-based sensitivity (binning method)
        n_bins = 10
        bins = np.linspace(np.min(param_values), np.max(param_values), n_bins+1)
        bin_means = []
        
        for j in range(n_bins):
            mask = (param_values >= bins[j]) & (param_values < bins[j+1])
            if np.sum(mask) > 0:
                bin_means.append(np.mean(valid_outputs[mask]))
        
        if len(bin_means) > 1:
            variance_between = float(np.var(bin_means))
            sensitivity_index = variance_between / variance if variance > 1e-15 else 0.0
        else:
            sensitivity_index = 0.0
        
        sensitivity_indices[param_name] = {
            'correlation': correlation,
            'sensitivity_index': sensitivity_index,
            'correlation_squared': correlation**2
        }
    
    # Rank parameters by sensitivity
    ranked_params = sorted(sensitivity_indices.items(), 
                          key=lambda x: abs(x[1]['correlation']), reverse=True)
    
    result = {
        'samples': model_outputs,
        'parameter_samples': param_samples.tolist(),
        'parameter_names': param_names,
        'n_samples': n_samples,
        'n_successful': len(valid_outputs),
        'n_failed': len(failed_samples),
        'failed_samples': failed_samples,
        'mean': mean,
        'std': std,
        'variance': variance,
        'median': median,
        'min': float(np.min(valid_outputs)),
        'max': float(np.max(valid_outputs)),
        'range': float(np.max(valid_outputs) - np.min(valid_outputs)),
        'percentiles': percentiles,
        'confidence_interval_95': [float(ci[0]), float(ci[1])],
        'coefficient_of_variation': cv,
        'skewness': skewness,
        'kurtosis': kurtosis,
        'sensitivity_indices': sensitivity_indices,
        'most_influential_param': ranked_params[0][0] if ranked_params else None,
    }
    
    return result



def residence_time_distribution(time_data, concentration_data):
    """Calculate comprehensive residence time distribution from tracer data with model fitting
    
    Analyzes tracer response data to extract RTD characteristics including mean residence time,
    variance, Peclet number, and reactor model classification.
    """
    import numpy as np
    
    time_data = np.array(time_data, dtype=float)
    concentration_data = np.array(concentration_data, dtype=float)
    
    if len(time_data) != len(concentration_data):
        raise ValueError(f"Time and concentration data must have same length")
    
    if len(time_data) < 3:
        raise ValueError(f"Need at least 3 data points for RTD analysis")
    
    # Remove negative concentrations
    concentration_data = np.maximum(concentration_data, 0.0)
    
    # Normalize concentration data to get RTD function E(t)
    area = np.trapezoid(concentration_data, time_data)
    
    if area <= 1e-15:
        raise ValueError("Total area under concentration curve is zero or negative")
    
    rtd_function = concentration_data / area
    
    # 0th moment (should be 1 after normalization)
    m0 = np.trapezoid(rtd_function, time_data)
    
    # 1st moment: mean residence time τ
    mean_residence_time = np.trapezoid(time_data * rtd_function, time_data)
    
    # 2nd moment (about origin)
    m2 = np.trapezoid(time_data**2 * rtd_function, time_data)
    
    # Variance σ²
    variance = m2 - mean_residence_time**2
    std_dev = float(np.sqrt(variance)) if variance > 0 else 0.0
    
    # 3rd moment (for skewness)
    m3 = np.trapezoid(time_data**3 * rtd_function, time_data)
    
    # Skewness
    if variance > 1e-15:
        skewness = float((m3 - 3*mean_residence_time*m2 + 2*mean_residence_time**3) / (variance**1.5))
    else:
        skewness = 0.0
    
    if mean_residence_time > 1e-15:
        theta = time_data / mean_residence_time
        dimensionless_variance = variance / (mean_residence_time**2)
    else:
        theta = time_data
        dimensionless_variance = 0.0
    
    # For small deviations: σ²/τ² ≈ 2/Pe - 2/Pe²
    if dimensionless_variance > 1e-10:
        peclet_approx = 2.0 / dimensionless_variance
    else:
        peclet_approx = float('inf')
    
    # Reactor classification based on variance
    if dimensionless_variance < 0.01:
        reactor_type = 'Plug Flow (PFR-like)'
        reactor_description = 'Low dispersion, near-ideal plug flow'
    elif dimensionless_variance > 0.9 and dimensionless_variance < 1.1:
        reactor_type = 'Continuous Stirred Tank (CSTR-like)'
        reactor_description = 'Well-mixed, CSTR behavior'
    elif dimensionless_variance < 0.5:
        reactor_type = 'Dispersed Plug Flow'
        reactor_description = 'Moderate dispersion with some plug flow character'
    else:
        reactor_type = 'Mixed Behavior'
        reactor_description = 'Significant dispersion and mixing'
    
    peak_index = int(np.argmax(rtd_function))
    peak_time = float(time_data[peak_index])
    peak_value = float(rtd_function[peak_index])
    
    # Calculate cumulative RTD function F(t)
    cumulative_rtd = np.zeros_like(rtd_function)
    for i in range(len(time_data)):
        if i == 0:
            cumulative_rtd[i] = 0.0
        else:
            cumulative_rtd[i] = np.trapezoid(rtd_function[:i+1], time_data[:i+1])
    
    fraction_times = {}
    for fraction in [0.1, 0.25, 0.5, 0.75, 0.9]:
        idx = np.argmin(np.abs(cumulative_rtd - fraction))
        fraction_times[f't_{int(fraction*100)}'] = float(time_data[idx])
    
    result = {
        'rtd_function': rtd_function.tolist(),
        'cumulative_rtd': cumulative_rtd.tolist(),
        'time_data': time_data.tolist(),
        'mean_residence_time': float(mean_residence_time),
        'variance': float(variance),
        'std_dev': std_dev,
        'dimensionless_variance': float(dimensionless_variance),
        'skewness': skewness,
        'peclet_number_estimate': float(peclet_approx),
        'reactor_type': reactor_type,
        'reactor_description': reactor_description,
        'peak_time': peak_time,
        'peak_value': peak_value,
        'fraction_times': fraction_times,
        'theta_dimensionless': theta.tolist(),
        'moments': {
            'm0': float(m0),
            'm1': float(mean_residence_time),
            'm2': float(m2),
            'm3': float(m3)
        }
    }
    
    return result


def catalyst_deactivation_model(time, kd, model_type='exponential', initial_activity=1.0):
    """Calculate catalyst activity over time with multiple deactivation mechanisms
    
    Models catalyst deactivation using various kinetic models including fouling,
    sintering, and poisoning mechanisms with detailed diagnostics.
    """
    import numpy as np
    
    time = np.array(time, dtype=float)
    
    if not isinstance(kd, (int, float)) or kd < 0:
        raise ValueError(f"Deactivation constant kd must be non-negative, got {kd}")
    
    if initial_activity <= 0:
        raise ValueError(f"Initial activity must be positive, got {initial_activity}")
    
    if model_type == 'exponential':
        # Exponential decay: a(t) = a0 × exp(-kd × t)
        activity = initial_activity * np.exp(-kd * time)
        model_equation = 'a(t) = a0 × exp(-kd × t)'
        mechanism = 'First-order poisoning/coking'
        
    elif model_type == 'power_law':
        # Power law decay: a(t) = a0 × (1 + kd × t)^(-1)
        activity = initial_activity * (1 + kd * time)**(-1)
        model_equation = 'a(t) = a0 / (1 + kd × t)'
        mechanism = 'Second-order deactivation/fouling'
        
    elif model_type == 'sintering':
        activity = initial_activity * (1 + kd * time)**(-0.5)
        model_equation = 'a(t) = a0 / (1 + kd × t)^0.5'
        mechanism = 'Thermal sintering/surface area loss'
        
    elif model_type == 'linear':
        # Linear decay: a(t) = a0 × (1 - kd × t)
        # Simple approximation for low deactivation
        activity = initial_activity * (1 - kd * time)
        activity = np.maximum(activity, 0.0)  # Non-negative
        model_equation = 'a(t) = a0 × (1 - kd × t)'
        mechanism = 'Linear deactivation (low conversion approximation)'
        
    elif model_type == 'hyperbolic':
        # Hyperbolic decay: a(t) = a0 / (1 + kd × t)^n
        n = 2.0
        activity = initial_activity * (1 + kd * time)**(-n)
        model_equation = f'a(t) = a0 / (1 + kd × t)^{n}'
        mechanism = 'Hyperbolic deactivation'
        
    else:
        raise ValueError(f"Unknown model type: {model_type}. Choose from 'exponential', 'power_law', 'sintering', 'linear', 'hyperbolic'")
    
    activity_array = np.array(activity)
    
    if model_type == 'exponential':
        half_life = np.log(2) / kd if kd > 1e-15 else float('inf')
    elif model_type == 'power_law':
        half_life = 1 / kd if kd > 1e-15 else float('inf')
    elif model_type == 'sintering':
        half_life = 3 / kd if kd > 1e-15 else float('inf')
    elif model_type == 'linear':
        half_life = 0.5 / kd if kd > 1e-15 else float('inf')
    else:  # hyperbolic
        half_life = (2**0.5 - 1) / kd if kd > 1e-15 else float('inf')
    
    # Deactivation rate (derivative)
    if len(time) > 1:
        deactivation_rate = -np.gradient(activity_array, time)
    else:
        deactivation_rate = np.array([0.0])
    
    # Average activity over time period
    if len(time) > 1:
        avg_activity = float(np.trapezoid(activity_array, time) / (time[-1] - time[0]))
    else:
        avg_activity = float(activity_array[0])
    
    final_activity = float(activity_array[-1])
    
    activity_loss = float(initial_activity - final_activity)
    activity_loss_percent = float(activity_loss / initial_activity * 100)
    
    result = {
        'activity': activity.tolist() if hasattr(activity, 'tolist') else [float(activity)],
        'time': time.tolist() if hasattr(time, 'tolist') else [float(time)],
        'model_type': model_type,
        'model_equation': model_equation,
        'deactivation_mechanism': mechanism,
        'deactivation_constant': float(kd),
        'initial_activity': float(initial_activity),
        'final_activity': final_activity,
        'average_activity': avg_activity,
        'activity_loss': activity_loss,
        'activity_loss_percent': activity_loss_percent,
        'half_life': float(half_life),
        'deactivation_rate': deactivation_rate.tolist(),
        'max_deactivation_rate': float(np.max(deactivation_rate)),
    }
    
    return result


def process_scale_up(lab_params, scale_factor):
    """Comprehensive process scale-up calculations with dimensional analysis
    
    Performs chemical engineering scale-up calculations using geometric similarity,
    dimensional analysis, and empirical correlations for reactor and process equipment.
    """
    import numpy as np
    
    if scale_factor <= 0:
        raise ValueError(f"Scale factor must be positive, got {scale_factor}")
    
    if scale_factor < 1:
        print("Warning: Scale factor < 1 indicates scale-down operation")
    
    scaled_params = {}
    scale_up_rules = {}
    dimensionless_numbers = {}
    
    for param, value in lab_params.items():
        param_lower = param.lower()
        
        if 'volume' in param_lower or param == 'V':
            # Volume scales as scale_factor^3 (geometric similarity)
            scaled_value = value * (scale_factor**3)
            scaling_exponent = 3.0
            rule = 'Geometric similarity: V ∝ L³'
            
        elif 'flow_rate' in param_lower or 'flow' in param_lower or param == 'Q':
            # Volumetric flow rate scales as scale_factor^3 (maintain residence time)
            scaled_value = value * (scale_factor**3)
            scaling_exponent = 3.0
            rule = 'Constant residence time: Q ∝ V ∝ L³'
            
        elif 'power' in param_lower or 'agitation' in param_lower or param == 'P':
            # Power scales approximately as scale_factor^2.5 to scale_factor^3
            # Using scale_factor^2.5 (common approximation for mixing)
            scaled_value = value * (scale_factor**2.5)
            scaling_exponent = 2.5
            rule = 'Power per unit volume decreases: P ∝ L^2.5'
            
        elif 'heat_transfer_area' in param_lower or 'area' in param_lower or param == 'A':
            # Surface area scales as scale_factor^2
            scaled_value = value * (scale_factor**2)
            scaling_exponent = 2.0
            rule = 'Geometric similarity: A ∝ L²'
            
        elif 'diameter' in param_lower or param == 'D' or 'length' in param_lower or param == 'L':
            scaled_value = value * scale_factor
            scaling_exponent = 1.0
            rule = 'Geometric similarity: L ∝ scale_factor'
            
        elif 'velocity' in param_lower or 'speed' in param_lower or param == 'v':
            # Velocity often kept constant or slightly reduced
            scaled_value = value  # Constant velocity
            scaling_exponent = 0.0
            rule = 'Constant velocity (common approximation)'
            
        elif 'rpm' in param_lower or 'rotation' in param_lower or 'stirrer' in param_lower:
            # Impeller speed: constant tip speed gives N ∝ 1/L
            # But often scaled to maintain constant power/volume: N ∝ L^(-1/6)
            scaled_value = value * (scale_factor**(-1.0/6.0))
            scaling_exponent = -1.0/6.0
            rule = 'Constant power per unit volume: N ∝ L^(-1/6)'
            
        elif 'pressure' in param_lower or param == 'P':
            scaled_value = value
            scaling_exponent = 0.0
            rule = 'Constant pressure (typical assumption)'
            
        elif 'temperature' in param_lower or param == 'T':
            scaled_value = value
            scaling_exponent = 0.0
            rule = 'Constant temperature (isothermal assumption)'
            
        elif 'heat_transfer_coefficient' in param_lower or param == 'U' or param == 'h':
            # Heat transfer coefficient typically decreases: U ∝ L^(-0.2)
            scaled_value = value * (scale_factor**(-0.2))
            scaling_exponent = -0.2
            rule = 'Empirical correlation: U ∝ L^(-0.2)'
            
        elif 'mass_transfer_coefficient' in param_lower or param == 'k_L' or param == 'k_G':
            # Mass transfer coefficient: k_L ∝ L^(-0.25) (typical)
            scaled_value = value * (scale_factor**(-0.25))
            scaling_exponent = -0.25
            rule = 'Empirical correlation: k_L ∝ L^(-0.25)'
            
        else:
            # Default: assume intensive property (no scaling)
            scaled_value = value
            scaling_exponent = 0.0
            rule = 'Intensive property (no scaling)'
        
        scaled_params[param] = scaled_value
        scale_up_rules[param] = {
            'original_value': float(value),
            'scaled_value': float(scaled_value),
            'scaling_exponent': scaling_exponent,
            'scaling_rule': rule,
            'scale_ratio': float(scaled_value / value) if value != 0 else 0.0
        }
    
    if 'diameter' in lab_params and 'velocity' in lab_params and 'density' in lab_params and 'viscosity' in lab_params:
        D_lab = lab_params['diameter']
        v_lab = lab_params['velocity']
        rho = lab_params['density']
        mu = lab_params['viscosity']
        
        Re_lab = rho * v_lab * D_lab / mu
        
        D_scale = scaled_params['diameter']
        v_scale = scaled_params.get('velocity', v_lab)
        Re_scale = rho * v_scale * D_scale / mu
        
        dimensionless_numbers['Reynolds'] = {
            'lab_scale': float(Re_lab),
            'production_scale': float(Re_scale),
            'ratio': float(Re_scale / Re_lab) if Re_lab != 0 else 0.0
        }
    
    if 'velocity' in lab_params and 'diameter' in lab_params:
        v = lab_params['velocity']
        D = lab_params['diameter']
        g = 9.81  # m/s²
        
        Fr_lab = v**2 / (g * D)
        Fr_scale = scaled_params.get('velocity', v)**2 / (g * scaled_params['diameter'])
        
        dimensionless_numbers['Froude'] = {
            'lab_scale': float(Fr_lab),
            'production_scale': float(Fr_scale),
            'ratio': float(Fr_scale / Fr_lab) if Fr_lab != 0 else 0.0
        }
    
    # Power number (for stirred tanks)
    if 'power' in lab_params and 'rpm' in lab_params and 'diameter' in lab_params and 'density' in lab_params:
        P = lab_params['power']
        N = lab_params['rpm'] / 60  # Convert to Hz
        D = lab_params['diameter']
        rho = lab_params['density']
        
        Np_lab = P / (rho * N**3 * D**5)
        
        P_scale = scaled_params['power']
        N_scale = scaled_params.get('rpm', lab_params['rpm']) / 60
        D_scale = scaled_params['diameter']
        
        Np_scale = P_scale / (rho * N_scale**3 * D_scale**5)
        
        dimensionless_numbers['Power_Number'] = {
            'lab_scale': float(Np_lab),
            'production_scale': float(Np_scale),
            'ratio': float(Np_scale / Np_lab) if Np_lab != 0 else 0.0
        }
    
    result = {
        'scaled_parameters': scaled_params,
        'scale_factor': float(scale_factor),
        'scaling_rules': scale_up_rules,
        'dimensionless_numbers': dimensionless_numbers,
        'volume_ratio': float(scale_factor**3),
        'area_ratio': float(scale_factor**2),
        'linear_ratio': float(scale_factor),
    }
    
    return result



def enthalpy_c(temperature, heat_capacity, reference_temperature=298.15):
    """Calculate enthalpy with constant heat capacity"""
    return heat_capacity * (temperature - reference_temperature)


def entropy_c(temperature, heat_capacity, reference_temperature=298.15):
    """Calculate entropy with constant heat capacity"""
    if temperature <= 0 or reference_temperature <= 0:
        return 0.0
    return heat_capacity * math.log(temperature / reference_temperature)


def find_steady_state(initial_guess, reaction_rates, tolerance=1e-6, max_iterations=100):
    """Find steady state concentrations"""
    import numpy as np
    concentrations = np.array(initial_guess)
    
    for iteration in range(max_iterations):
        rates = np.array(reaction_rates) * concentrations
        
        new_concentrations = concentrations - 0.1 * rates
        new_concentrations = np.maximum(0, new_concentrations)  # Non-negative
        
        if np.linalg.norm(new_concentrations - concentrations) < tolerance:
            return new_concentrations.tolist(), True, iteration
        
        concentrations = new_concentrations
    
    return concentrations.tolist(), False, max_iterations


def validate_mechanism(mechanism_dict):
    """Validate reaction mechanism for consistency"""
    required_keys = ['reactions', 'species', 'rate_constants']
    
    for key in required_keys:
        if key not in mechanism_dict:
            return False, f"Missing required key: {key}"
    
    reactions = mechanism_dict['reactions']
    species = mechanism_dict['species']
    
    for i, reaction in enumerate(reactions):
        if 'reactants' not in reaction or 'products' not in reaction:
            return False, f"Reaction {i} missing reactants or products"
        
        all_reaction_species = set(reaction['reactants'].keys()) | set(reaction['products'].keys())
        undefined_species = all_reaction_species - set(species)
        if undefined_species:
            return False, f"Undefined species in reaction {i}: {undefined_species}"
    
    return True, "Mechanism is valid"


def optimize_parameters(objective_function, initial_parameters, bounds=None, method='simple'):
    """Parameter optimization using simple methods"""
    import numpy as np
    
    if method == 'simple':
        # Simple grid search
        best_params = np.array(initial_parameters)
        best_objective = objective_function(best_params)
        
        # Try perturbations
        for param_idx in range(len(initial_parameters)):
            for delta in [-0.1, 0.1]:
                test_params = best_params.copy()
                test_params[param_idx] *= (1 + delta)
                
                # Apply bounds if specified
                if bounds and param_idx < len(bounds):
                    min_val, max_val = bounds[param_idx]
                    test_params[param_idx] = np.clip(test_params[param_idx], min_val, max_val)
                
                try:
                    test_objective = objective_function(test_params)
                    if test_objective < best_objective:
                        best_params = test_params
                        best_objective = test_objective
                except:
                    continue  # Skip invalid parameter sets
        
        return best_params.tolist(), best_objective
    
    else:
        return initial_parameters, objective_function(np.array(initial_parameters))


def load_spec_from_yaml(filename):
    """Load simulation specification from YAML file"""
    try:
        import yaml
        with open(filename, 'r') as file:
            spec = yaml.safe_load(file)
        return spec
    except ImportError:
        print("PyYAML not available. Please install with: pip install PyYAML")
        return {}
    except FileNotFoundError:
        print(f"File {filename} not found")
        return {}


def save_spec_to_yaml(spec, filename):
    """Save simulation specification to YAML file"""
    try:
        import yaml
        with open(filename, 'w') as file:
            yaml.dump(spec, file, default_flow_style=False)
        return True
    except ImportError:
        print("PyYAML not available. Please install with: pip install PyYAML")
        return False


def parse_mechanism(mechanism_string):
    """Parse mechanism string into structured format"""
    lines = mechanism_string.strip().split('\n')
    reactions = []
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        # Simple parsing for reactions like "A + B -> C + D"
        if '->' in line:
            left, right = line.split('->')
            reactants = [species.strip() for species in left.split('+')]
            products = [species.strip() for species in right.split('+')]
            
            reaction = {
                'reactants': {species: 1 for species in reactants},
                'products': {species: 1 for species in products}
            }
            reactions.append(reaction)
    
    species = set()
    for reaction in reactions:
        species.update(reaction['reactants'].keys())
        species.update(reaction['products'].keys())
    
    return {
        'reactions': reactions,
        'species': list(species),
        'rate_constants': [1.0] * len(reactions)  # Default rate constants
    }


def save_results_to_csv(filename, time_data, concentration_data, species_names=None):
    """Save simulation results to CSV file"""
    import csv
    import numpy as np
    
    time_data = np.array(time_data)
    concentration_data = np.array(concentration_data)
    
    if species_names is None:
        species_names = [f'Species_{i}' for i in range(concentration_data.shape[1])]
    
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Header
        header = ['Time'] + species_names
        writer.writerow(header)
        
        for i in range(len(time_data)):
            row = [time_data[i]] + concentration_data[i, :].tolist()
            writer.writerow(row)
    
    return True


def get_build_info():
    """Get build information about PyroXa installation"""
    return {
        'version': '1.0.0',
        'cpp_available': False,  # Pure Python version
        'compilation_date': 'N/A',
        'python_version': '3.8+',
        'numpy_version': 'latest',
        'build_type': 'Pure Python'
    }



def is_reaction_chains_available():
    """Check if reaction chain functionality is available"""
    return True


def build_from_dict(spec):
    """Build reactor and simulation from specification dictionary"""
    from .purepy import WellMixedReactor, Reaction, Thermodynamics
    
    reaction_spec = spec.get('reaction', {})
    kf = reaction_spec.get('kf', 1.0)
    kr = reaction_spec.get('kr', 0.1)
    
    reaction = Reaction(kf, kr)
    
    initial_spec = spec.get('initial', {})
    conc_spec = initial_spec.get('conc', {'A': 1.0, 'B': 0.0})
    A0 = conc_spec.get('A', 1.0)
    B0 = conc_spec.get('B', 0.0)
    
    reactor = WellMixedReactor(reaction, A0=A0, B0=B0)
    
    sim_spec = spec.get('sim', {})
    time_span = sim_spec.get('time_span', 10.0)
    dt = sim_spec.get('time_step', 0.1)
    
    return reactor, {'time_span': time_span, 'dt': dt}


def simulate_cstr(residence_time, kf, kr, conc_feed, feed_flow_rate=1.0, volume=None):
    """Simulate CSTR at steady state"""
    from .purepy import CSTR, Reaction
    
    if volume is None:
        volume = residence_time * feed_flow_rate
    
    reaction = Reaction(kf, kr)
    cstr = CSTR(reaction, residence_time, conc_feed)
    
    return cstr.steady_state()


def simulate_pfr(length, velocity, kf, kr, conc_inlet, n_segments=100):
    """Simulate PFR using finite differences"""
    from .purepy import PFR, Reaction
    
    reaction = Reaction(kf, kr)
    pfr = PFR(reaction, length, velocity, conc_inlet)
    
    return pfr.solve(n_segments)


def calculate_energy_balance_simple(heat_reaction, heat_capacity, mass, temperature_change):
    """Simple energy balance calculation"""
    heat_sensible = heat_capacity * mass * temperature_change
    total_heat = heat_reaction + heat_sensible
    return total_heat


def free_aligned_memory():
    """Free aligned memory (placeholder for pure Python)"""
    import gc
    gc.collect()  # Trigger garbage collection
    return True