import numpy as np

thf18 = 119

def calc_act_con(volume,activity,deltat0, remainder,deltatr,th=thf18):
        
    print "volume: (ml)", volume
    print "activity: (MBq)", activity
    print "deltat0: ", deltat0
    print "remainder: ", remainder
    print "deltatr: ", deltatr 
 
    Ar = decay(remainder,deltatr,th) 
    A0 = decay((activity - Ar),deltat0,th) 

    print "Decay corrected Ar, A: ", Ar,A0

    accon = 1e6*(A0 - Ar)/volume
    print "Activity concentration: Bq/ml",accon
    return accon


def decay(A0,t,th):
    return A0*np.power(0.5,(t/th))
