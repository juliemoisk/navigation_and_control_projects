import numpy as np
import math
import matplotlib.pyplot as plt
f = open("aks_1.log", "r")
n=0
for l in f:
   n += 1
print(n)
f.close()
f = open("aks_1.log", "r")

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

plt.plot(data[:, 0], data[:, 4], 'g')
plt.plot(data[:, 0], data[:, 5], 'b')
plt.plot(data[:, 0], data[:, 6], 'y')
plt.savefig("Acs1.png", dpi = 200)
plt.show()

#plt.plot(data[:, 0], data[:, 7], 'g')
#plt.plot(data[:, 0], data[:, 8], 'b')
#plt.plot(data[:, 0], data[:, 9], 'y')
#plt.show()

t = np.empty([6, 2], int)
t[0][0] = 50
t[0][1] = 1000
t[1][0] = 1200
t[1][1] = 2035
t[2][0] = 2100
t[2][1] = 3050
t[3][0] = 3165
t[3][1] = 4080
t[4][0] = 4170
t[4][1] = 5100
t[5][0] = 5190
t[5][1] = 6260

a = np.zeros([6,3])

for i in range(6):
    for j in range(t[i][0], t[i][1]):
        a[i][0] += data[j][4]
        a[i][1] += data[j][5]
        a[i][2] += data[j][6]
    a[i][0] = a[i][0]/(t[i][1]-t[i][0])
    a[i][1] = a[i][1]/(t[i][1]-t[i][0])
    a[i][2] = a[i][2]/(t[i][1]-t[i][0])
print("Vectors a")
print(a)

g = 9.81

S_a = np.empty([3,3])

S_a[:, 0] = 1/2./g*(a[0]-a[1])
S_a[:, 1] = 1/2./g*(a[2]-a[3])
S_a[:, 2] = 1/2./g*(a[5]-a[4])

S_f = np.linalg.inv(S_a)

print("Matrix S_a")
print(S_a)
print("Matrix S_f")
print(S_f)

b_a = np.empty(3)
b_f = np.empty(3)

b_a = 1/6.*(a[0]+a[1]+a[2]+a[3]+a[4]+a[5])
b_f = -np.dot(S_f, b_a)
print("Vector b_a")
print(b_a)
print("Vector b_f")
print(b_f)

data_z = np.empty([3, n])
data_z[0] = data[:, 4]
data_z[1] = data[:, 5]
data_z[2] = data[:, 6]
#data_z = data_z.transpose()

f_z = np.empty([n, 3])
for i in range(n):
   f_z[i] = np.dot(S_f, data_z[:, i]) + b_f

plt.plot(data[:, 0], f_z[:, 0], 'g')
plt.plot(data[:, 0], f_z[:, 1], 'b')
plt.plot(data[:, 0], f_z[:, 2], 'y')
plt.show()
