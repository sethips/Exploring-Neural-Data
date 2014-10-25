#
#  NAME
#    problem_set1.py
#
#  DESCRIPTION
#    Open, view, and analyze raw extracellular data
#    In Problem Set 1, you will write create and test your own spike detector.
#

import numpy as np
import matplotlib.pylab as plt

def load_data(filename):
    """
    load_data takes the file name and reads in the data.  It returns two 
    arrays of data, the first containing the time stamps for when they data
    were recorded (in units of seconds), and the second containing the 
    corresponding voltages recorded (in units of microvolts - uV)
    """
    data = np.load(filename)[()];
    return np.array(data['time']), np.array(data['voltage'])
    
def bad_AP_finder(time,voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time
        
        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run
    
    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected
         
    This function is bad at detecting spikes!!! 
        But it's formated to get you started!
    """
    
    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        APTimes = []
        return APTimes
    
    numAPs = np.random.randint(0,len(time))//10000 #and this is why it's bad!!
 
    # Now just pick 'numAPs' random indices between 0 and len(time)
    APindices = np.random.randint(0,len(time),numAPs)
    
    # By indexing the time array with these indices, we select those times
    APTimes = time[APindices]
    
    # Sort the times
    APTimes = np.sort(APTimes)
    
    return APTimes
    
def good_AP_finder(time, voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time
        
        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run
    
    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected
    """
 
    APTimes = []
       
    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        return APTimes
    
    # Compute the average and standard deviation
    # These are used to determine which voltages are unusually large or small
    avg = np.average(voltage)
    std = np.std(voltage)

    # Some data might be more static than others- we need to shift the std. dev.
    ratio = avg / std

    while abs(ratio) < 0.1:
        ratio *= 10

    std += ratio

    tempRes = []

    for i in range(0, len(voltage)):
        if ((voltage[i] - avg) / std) > 3.5:

            # Check if there was a negative spike as well
            for j in range(0, 100):

                # Negative spike to the left
                if abs(time[i + j] - time[i]) < 0.003:
                    if ((avg - voltage[i + j]) / std) < -3:
                        tempRes.append(time[i])
                        break

                # Negative spike to the right
                if abs(time[i - j] - time[i]) < 0.003:
                    if ((avg - voltage[i - j]) / std) < -3:
                        tempRes.append(time[i])
                        break

    # Append them to the temporary array
    for i in range(0, len(tempRes) - 1):
        if abs(tempRes[i] - tempRes[i + 1]) < 0.0015:
            tempRes[i] = -1

    # If there was a spike found, append it to the final result
    for i in range(0, len(tempRes)):
        if tempRes[i] != -1:
            APTimes.append(tempRes[i])
    
    print("Number of spikes found: ", len(APTimes))
    return APTimes
    

def get_actual_times(dataset):
    """
    Load answers from dataset
    This function takes the following input:
        dataset - name of the dataset to get answers for

    This function returns the following output:
        APTimes - spike times
    """    
    return np.load(dataset)
    
def detector_tester(APTimes, actualTimes):
    """
    returns percentTrueSpikes (% correct detected) and falseSpikeRate
    (extra APs per second of data)
    compares actual spikes times with detected spike times
    This only works if we give you the answers!
    """
    
    JITTER = 0.025 #2 ms of jitter allowed
    
    #first match the two sets of spike times. Anything within JITTER_MS
    #is considered a match (but only one per time frame!)
    
    #order the lists
    detected = np.sort(APTimes)
    actual = np.sort(actualTimes)
    
    #remove spikes with the same times (these are false APs)
    temp = np.append(detected, -1)
    detected = detected[plt.find(plt.diff(temp) != 0)]
 
    #find matching action potentials and mark as matched (trueDetects)
    trueDetects = [];
    for sp in actual:
        z = plt.find((detected >= sp-JITTER) & (detected <= sp+JITTER))
        if len(z)>0:
            for i in z:
                zz = plt.find(trueDetects == detected[i])
                if len(zz) == 0:
                    trueDetects = np.append(trueDetects, detected[i])
                    break;
    percentTrueSpikes = 100.0*len(trueDetects)/len(actualTimes)
    
    #everything else is a false alarm
    totalTime = (actual[len(actual)-1]-actual[0])
    falseSpikeRate = (len(APTimes) - len(actualTimes))/totalTime
    
    print 'Action Potential Detector Performance performance: '
    print '     Correct number of action potentials = ' + str(len(actualTimes))
    print '     Percent True Spikes = ' + str(percentTrueSpikes)
    print '     False Spike Rate = ' + str(falseSpikeRate) + ' spikes/s'
    print 
    return {'Percent True Spikes':percentTrueSpikes, 'False Spike Rate':falseSpikeRate}
    
    
def plot_spikes(time, voltage, APTimes, titlestr):
    plt.figure()
    plt.plot(time, voltage, APTimes, max(voltage) * np.ones(len(APTimes)), '|r')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (uV)')
    plt.title(titlestr)
    plt.show()
    
def plot_waveforms(time, voltage, APTimes, titlestr):
    """
    plot_waveforms takes four arguments - the recording time array, the voltage
    array, the time of the detected action potentials, and the title of your
    plot.  The function creates a labeled plot showing the waveforms for each
    detected action potential
    """

    all_times = []
    all_voltages = []
    
    # Goes through the spikes and finds the corresponding t/v values
    for i in range(0, len(APTimes)):
        current_times = []
        current_voltages = []

        for j in range(0, len(time)):
            # Record all t/v pairs within 3ms of the spike
            if abs(time[j] - APTimes[i]) < 0.003:
                current_times.append(time[j] - APTimes[i])
                current_voltages.append(voltage[j])

        # Push the t/v set around the current spike to the master list
        all_times.append(current_times)
        all_voltages.append(current_voltages)

    plt.figure()

    # Need to plot each individually to see proper forms
    for i in range(0, len(all_times)):
        plt.plot(all_times[i], all_voltages[i], '-b')

    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (uV)')
    plt.title(titlestr)

    plt.show()


    

        
##########################
#You can put the code that calls the above functions down here    

if __name__ == "__main__":
    t,v = load_data('spikes_hard_practice.npy')    
    actualTimes = get_actual_times('spikes_hard_practice_answers.npy')
    APTime = good_AP_finder(t,v)
    plot_spikes(t, v, APTime, 'Action Potential in Raw Signal')
    plot_waveforms(t, v, APTime, 'Waveforms')
    detector_tester(APTime, actualTimes)


