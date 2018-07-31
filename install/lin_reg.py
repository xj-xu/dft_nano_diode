#!/usr/bin/python -tt

# linear regression for the sampling of data point between t=a and t=b
#   hypothesis: h_theta(t) = theta0 + theta1*t + theta2*sin(omega*t) + theta3*cos(omega*t)
# output:
#   full fit curve
#   linear fit curve
#   linear fit slope (theta1)
#   data sampled over


"""===================================="""
import sys
import numpy as np
from scipy.linalg import eigh
from numpy.linalg import pinv
import matplotlib.pyplot as plt


# Provided main()
def main():
  if len(sys.argv) != 4:
    print 'usage: ./lin_reg.py <start time point> <end time point> <laser wavelength>'
    sys.exit(1)
  
  hbar = 6.582119514E-16       ## plank's constant in eV s
  c_speed = 2.99792458E8       ## speed of light in m/s
  t_start = float(sys.argv[1]) ## start time in femtoseconds
  t_stop = float(sys.argv[2])  ## stop time in femtoseconds
  omega = 1E-6 * 2.0 * np.pi * c_speed / float(sys.argv[3]) ## ang frequency of sinusoidal component 

  # read in all data points from file
  f = open('fort.400', 'rU') # opens file into the variable f
  N_points=0
  time=[]
  value=[]
  for l in f:
    row = l.split()
    time.append(float(row[0]))
    value.append(float(row[1]))
    N_points = N_points + 1 ## number of time points 
#    print row[0], row[1]
  f.close()

  # define data points used in linear regression
  t = 0
  i = 0
  t_arr = []
  y_arr = []
  while (t < t_stop):
    if time[i] > t_start: 
      t_arr.append(time[i])
      y_arr.append(value[i])
    t = time[i] # technically will keep first point AFTER t_stop but that's ok
    i = i + 1
  m = len(t_arr) ## number of data points


  #prepare linear regression variables
  s_data_out = open('sample_data.dat', 'w')
  n = 4 ## number of linear regression parameters
  arr =[]
  for i in range(m):
    s_data_out.write(str(t_arr[i]) + " " + str(y_arr[i]) + '\n')
    row= []
    for j in range(n):
      if j == 0: elem = 1.0
      if j == 1: elem = t_arr[i]
      if j == 2: elem = np.sin(omega*t_arr[i])
      if j == 3: elem = np.cos(omega*t_arr[i])
      row.append(elem)
    arr.append(row)
  x_mat = np.array(arr)
  theta_arr = [0.0] * n

  # solve via Normal Equation
  theta_arr = np.dot(np.dot(pinv(np.dot(np.transpose(x_mat),x_mat)),np.transpose(x_mat)),y_arr)

  print
  print 'theta0: ',theta_arr[0]
  print 'theta1: ',theta_arr[1]
  print 'theta2: ',theta_arr[2]
  print 'theta3: ',theta_arr[3]
  print
  print 'slope: ', theta_arr[1]
  print 'y-int: ', theta_arr[0]
  print
  print 'rate of electron transfer: ', theta_arr[1], ' electrons per fs'
  fit_out_log = open('fit.log', 'w')
  fit_out_log.write('theta0: '+str(theta_arr[0])+'\n'+'theta1: '+str(theta_arr[1])+'\n'+'theta2: '+str(theta_arr[2])+'\n'+'theta3: '+str(theta_arr[3])+'\n'+'\n')
  fit_out_log.write('slope: '+str(theta_arr[1])+'\n'+'y-int: '+str(theta_arr[0])+'\n')

  fit_out = open('fit.dat', 'w')
  linear_out = open('linear_fit.dat', 'w')
  for i in range(N_points):
    hyp = theta_arr[0] + theta_arr[1]*time[i] + theta_arr[2]*np.sin(omega*time[i]) + theta_arr[3]*np.cos(omega*time[i])
    hyp_line = theta_arr[0] + theta_arr[1]*time[i]
    fit_out.write(str(time[i]) + " " + str(hyp) +'\n')
    linear_out.write(str(time[i]) + " " + str(hyp_line) +'\n')   

if __name__ == '__main__':
  main()


