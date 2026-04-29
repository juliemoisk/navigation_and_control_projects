import numpy as np
import math
import matplotlib.pyplot as plt
f = open("dus_alt_speed.log", "r")
n=0
for l in f:
   n += 1
print(n)
f.close()
f = open("dus_alt_speed.log", "r")

data = np.empty([n, 14])

i = 0
for l in f:
    q = str.split(l)
    for j in range(14):
        data[i][j] = int(q[j])
    i += 1    

c = data[0][0]
for i in range(n):
    data[i][0] = data[i][0] - c

f.close()

#plt.plot(data[:, 0], data[:, 1], 'g')
#plt.plot(data[:, 0], data[:, 2], 'b')
#plt.plot(data[:, 0], data[:, 3], 'y')
#plt.show()

#plt.plot(data[:, 0], data[:, 4], 'g')
#plt.plot(data[:, 0], data[:, 5], 'b')
#plt.plot(data[:, 0], data[:, 6], 'y')
#plt.savefig("Acs1.png", dpi = 200)
#plt.show()

plt.plot(data[:, 0], data[:, 7], 'g')
plt.plot(data[:, 0], data[:, 8], 'b')
plt.plot(data[:, 0], data[:, 9], 'y')
plt.savefig("DUS.png", dpi = 200)
plt.show()

t = np.empty([6, 2], int)
t[0][0] = 55
t[0][1] = 195
t[1][0] = 474
t[1][1] = 608
t[2][0] = 1017
t[2][1] = 1169
t[3][0] = 1494
t[3][1] = 1639
t[4][0] = 2132
t[4][1] = 2265
t[5][0] = 2582
t[5][1] = 2712

t0 = 612
t1 = 805
b_w = np.zeros(3)

for i in range(t0, t1):
    b_w[0] += data[i][7]
    b_w[1] += data[i][8]
    b_w[2] += data[i][9]
b_w = b_w/(t1-t0)
print("Vector b_w")
print(b_w)

dphi = np.zeros([6,3])
for i in range(6):
    for j in range(t[i][0],t[i][1]):
        dphi[i][0] += data[j][7] - b_w[0]
        dphi[i][1] += data[j][8] - b_w[1]
        dphi[i][2] += data[j][9] - b_w[2]
    dphi[i] = dphi[i]*0.02
#print("Vector dphi")
#print(dphi)

S_w = np.empty([3,3])

for i in range(3):
    S_w[:,i] = 1/3./math.pi*(dphi[2*i]-dphi[2*i+1])
print("Matrix S_w")
print(S_w)


S_om = np.linalg.inv(S_w)

print("Matrix S_om")
print(S_om)

b_om = -np.dot(S_om, b_w)

print("Vector b_om")
print(b_om)

data_z = np.empty([3, n])
data_z[0] = data[:, 7]
data_z[1] = data[:, 8]
data_z[2] = data[:, 9]

f_z = np.empty([n, 3])
for i in range(n):
   f_z[i] = np.dot(S_om, data_z[:, i]) + b_om

plt.plot(data[:, 0], f_z[:, 0], 'g')
plt.plot(data[:, 0], f_z[:, 1], 'b')
plt.plot(data[:, 0], f_z[:, 2], 'y')
plt.show()
