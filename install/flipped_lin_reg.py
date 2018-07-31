#!/usr/bin/python -tt

# linear regression for the sampling of data point between t=a and t=b
#   hypothesis: h_theta(t) = theta0 + theta1*t + theta2*sin(omega*t) + theta3*cos(omega*t)
# output:
#   full fit curve
#   linear fit curve
#   linear fit slope (theta1)
#   data sampled over

  
# theta_j := theta_j - alpha * d/dtheta_j (1/2m)sum_{i=1}^m (h_theta(t^(i)) - y^(i))^2 
# linear regression:
# theta_j := theta_j - alpha * (1/m) sum_{i=1}^m (h_theta(t^(i)) - y^(i)) * t_j

"""===================================="""
import sys
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt


# Cost Function
def cost_func(theta,x,y): 
  su = 0.0
  for i in range(len(y)):
    su = su + (np.dot(theta,x[i,:]) - y[i])**2
  value = (1.0/(2.0*len(y)))*su
  return value


# Provided main()
def main():
  if len(sys.argv) != 4:
    print 'usage: ./lin_reg.py <end time point> <laser wavelength> <alpha>'
    sys.exit(1)
  
  hbar = 6.582119514E-16 ## plank's constant in eV s
  c_speed = 2.99792458E8 ## speed of light in m/s
  t_start = 3.0 ## start time in femtoseconds
  t_stop = float(sys.argv[1]) ## stop time in femtoseconds
  omega = 1E-6 * 2.0 * np.pi * c_speed / float(sys.argv[2]) ## ang frequency of sinusoidal component 
  alpha = float(sys.argv[3]) ## lin reg step size

  # read in all data points from file
  f = open('fort.400new', 'rU') # opens file into the variable f
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
  dt = t_arr[m-1] - t_arr[0] ## range of t values
  arr =[]
  for i in range(m):
    s_data_out.write(str(t_arr[i]) + " " + str(y_arr[i]) + '\n')
    row= []
    for j in range(n):
      if j == 0: elem = 1.0
      if j == 1: elem = -0.5 + (t_arr[i]-t_start) / dt # scaled to be between -1 and 1
#      if j == 1: elem = t_arr[i] # not scaled
      if j == 2: elem = np.sin(omega*t_arr[i])
      if j == 3: elem = np.cos(omega*t_arr[i])
      row.append(elem)
    arr.append(row)
  x_mat = np.array(arr)
  theta_arr = [0.0] * n

  # perform linear regressing
  it=0
  J_diff = 1E15
  J_old = 1E15
  cf_out = open('cost_func.dat', 'w')
  temp = [0.0] * n
  while J_diff > 1E-6:
    su = [0.0] * n
    for j in range(n):
      for i in range(m):
        su[j] = su[j] + ( np.dot(theta_arr,x_mat[i,:]) - y_arr[i] ) * x_mat[i,j]
      temp[j] = theta_arr[j] - (alpha / m) * su[j]
    for j in range(n):
      theta_arr[j] = temp[j]
    J = cost_func(theta_arr,x_mat,y_arr)
    cf_out.write(str(it) + " " + str(J)+'\n')
    J_diff=abs(J-J_old)
#    print it,J
    print it,J_diff
    J_old=J
    it = it + 1

  print
  print 'theta0: ',theta_arr[0]
  print 'theta1: ',theta_arr[1]
  print 'theta2: ',theta_arr[2]
  print 'theta3: ',theta_arr[3]
  print
  fit_out_log = open('fit.log', 'w')
  fit_out_log.write('theta0: '+str(theta_arr[0])+'\n')
  fit_out_log.write('theta1: '+str(theta_arr[1])+'\n')
  fit_out_log.write('theta2: '+str(theta_arr[2])+'\n')
  fit_out_log.write('theta3: '+str(theta_arr[3])+'\n')
  fit_out_log.write('\n') 
  theta_arr[0] = theta_arr[0]-(t_start/dt + 0.5)*theta_arr[1]
  theta_arr[1] = theta_arr[1]/dt
  print 'slope: ', theta_arr[1]
  print 'y-int: ', theta_arr[0]
  print
  print 'rate of electron transfer: ', theta_arr[1], ' electrons per fs'
  fit_out_log.write('slope: '+str(theta_arr[1])+'\n')
  fit_out_log.write('y-int: '+str(theta_arr[0])+'\n')

  fit_out = open('fit.dat', 'w')
  linear_out = open('linear_fit.dat', 'w')
  for i in range(N_points):
    hyp = theta_arr[0] + theta_arr[1]*time[i] + theta_arr[2]*np.sin(omega*time[i]) + theta_arr[3]*np.cos(omega*time[i])
    hyp_line = theta_arr[0] + theta_arr[1]*time[i]
    fit_out.write(str(time[i]) + " " + str(hyp) +'\n')
    linear_out.write(str(time[i]) + " " + str(hyp_line) +'\n')   





if __name__ == '__main__':
  main()















