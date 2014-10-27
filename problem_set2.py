#
#  NAME
#    problem_set2_solutions.py
#
#  DESCRIPTION
#    Open, view, and analyze action potentials recorded during a behavioral
#    task.  In Problem Set 2, you will write create and test your own code to
#    create tuning curves.
#

#Helper code to import some functions we will use
import numpy as np
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from scipy import optimize
from scipy import stats


def load_experiment(filename):
    """
    load_experiment takes the file name and reads in the data.  It returns a
    two-dimensional array, with the first column containing the direction of
    motion for the trial, and the second column giving you the time the
    animal began movement during thaht trial.
    """
    data = np.load(filename)[()];
    return np.array(data)

def load_neuraldata(filename):
    """
    load_neuraldata takes the file name and reads in the data for that neuron.
    It returns an arary of spike times.
    """
    data = np.load(filename)[()];
    return np.array(data)
    
def bin_spikes(trials, spk_times, time_bin):
    """
    bin_spikes takes the trials array (with directions and times) and the spk_times
    array with spike times and returns the average firing rate for each of the
    eight directions of motion, as calculated within a time_bin before and after
    the trial time (time_bin should be given in seconds).  For example,
    time_bin = .1 will count the spikes from 100ms before to 100ms after the 
    trial began.
    
    dir_rates should be an 8x2 array with the first column containing the directions
    (in degrees from 0-360) and the second column containing the average firing rate
    for each direction
    """
    num_angles = len(np.unique(trials[:, 0]))
    # First column represents number of spikes around each trial
    # Second column represents number of trials
    spikes_trials = np.column_stack((np.zeros(num_angles), np.zeros(num_angles)))

    for i in range(0, len(trials)):
        curr_angle = np.int32(trials[i][0])
        curr_time = trials[i][1]
        spikes_trials[curr_angle / 45][1] += 1

        for j in range(0, len(spk_times)):
            if abs(spk_times[j] - curr_time) < time_bin:
                spikes_trials[curr_angle / 45][0] += 1

    dir_rates = np.column_stack((np.arange(0, 360, 45), np.zeros(num_angles)))

    for i in range(0, num_angles):
        dir_rates[i][1] = spikes_trials[i][0] / (2 * time_bin * spikes_trials[i][1] )
    
    return dir_rates
    
def plot_tuning_curves(direction_rates, title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates) and plots a histogram and 
    polar representation of the tuning curve. It adds the given title.
    """

    
def roll_axes(dir_rates):
    """
    roll_axes takes the x-values (directions) and y-values (direction_rates)
    and return new x and y values that have been "rolled" to put the maximum
    direction_rate in the center of the curve. The first and last y-value in the
    returned list should be set to be the same. (See problem set directions)
    Hint: Use np.roll()
    """
    shift = np.argmax(dir_rates[:, 1])
    new_xs = np.roll(dir_rates[:, 0], shift)
    new_ys = np.roll(dir_rates[:, 1], shift)

    new_xs = np.append(new_xs, new_xs[0])
    new_ys = np.append(new_ys, new_ys[0])

    for i in range(0, len(new_xs) - 1):
        if (360 - new_xs[i]) <= (45 * shift):
            new_xs[i] = new_xs[i] - 360
    
    roll_degrees = shift * 45

    return new_xs, new_ys, roll_degrees    
    

def normal_fit(x, mu, sigma, A):
    """
    This creates a normal curve over the values in x with mean mu and
    variance sigma.  It is scaled up to height A.
    """
    n = A * mlab.normpdf(x, mu, sigma)
    return n

def fit_tuning_curve(centered_x, centered_y):
    """
    This takes our rolled curve, generates the guesses for the fit function,
    and runs the fit.  It returns the parameters to generate the curve.
    """

    return p
    


def plot_fits(direction_rates, fit_curve, title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates and fit_curve) and plots the 
    actual values with circles, and the curves as lines in both linear and 
    polar plots.
    """
    

def von_mises_fitfunc(x, A, kappa, l, s):
    """
    This creates a scaled Von Mises distrubition.
    """
    return A * stats.vonmises.pdf(x, kappa, loc=l, scale=s)


    
def preferred_direction(fit_curve):
    """
    The function takes a 2-dimensional array with the x-values of the fit curve
    in the first column and the y-values of the fit curve in the second.  
    It returns the preferred direction of the neuron (in degrees).
    """

    return np.argmax(fit_curve[1])
  

def plot_hist(dir_rates):

    # Histogram
    plt.subplot(2, 2, 1)
    plt.bar(dir_rates[:, 0], dir_rates[:, 1], 
        (np.max(dir_rates[:, 0]) - np.min(dir_rates[:, 0])) / (len(dir_rates[:, 0]) - 1))
    plt.xlim(xmax = 360)
    plt.xlabel('Direction of Motion (Degrees)')
    plt.ylabel('Firing Rate (spikes/s)')
    plt.title('Neuron Tuning Curve')

    # Polar 
    plt.subplot(2, 2, 2, polar = True)
    plt.polar(np.deg2rad(dir_rates[:, 0]), dir_rates[:, 1], 
        label = 'Neuron Tuning Curve')
    plt.legend(loc = 8)

def rfr(new_xs, new_ys, roll_degrees):
    # Fitting
    max_y = np.amax(new_ys)
    max_x = new_xs[np.argmax(new_ys)]
    sigma = 90
    p, cov = optimize.curve_fit(normal_fit, new_xs, new_ys, 
        p0 = [max_x, sigma, max_y])

    # Pruning
    new_xs += roll_degrees
    new_xs = np.roll(new_xs, np.argmin(new_xs))
    new_ys = np.roll(new_ys, np.argmin(new_ys))

    # Fitting
    curve_xs = np.arange(new_xs[0], new_xs[-1])
    fit_ys = normal_fit(curve_xs, p[0], p[1], p[2])

    # More pruning
    new_xs = new_xs[0 : len(new_xs) - 1]
    new_ys = new_ys[0 : len(new_ys) - 1]
    

    # Plotting
    plt.subplot(2, 2, 3)
    plt.plot(new_xs, new_ys, 'go', curve_xs, fit_ys, '-')
    plt.xlim(xmax = 360)
    plt.xlabel('Direction of Motion (Degrees)')
    plt.ylabel('Firing Rate (spikes/s)')
    plt.title('Neuron Tuning Curve')
    
    plt.subplot(2, 2, 4, polar = True)
    plt.polar(np.deg2rad(new_xs), new_ys, 'go', np.deg2rad(curve_xs), fit_ys, '-',
        label = 'Neuron Tuning Curve')
    plt.legend(loc = 8)

    plt.show()  

    return [curve_xs, fit_ys]


##########################
#You can put the code that calls the above functions down here   
 
if __name__ == "__main__":
    # Initialization
    trials = load_experiment('trials.npy')   
    spk_times = load_neuraldata('neuron3.npy') 
    dir_rates = bin_spikes(trials, spk_times, 0.1)
    plot_hist(dir_rates)
    new_xs, new_ys, roll_degrees = roll_axes(dir_rates)

    fit_curve = rfr(new_xs, new_ys, roll_degrees)
    pd = preferred_direction(fit_curve)
    print(pd)