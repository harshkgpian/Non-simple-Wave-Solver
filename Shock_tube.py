import matplotlib.pyplot as plt
import numpy as np
from math import tan, sin,cos,atan
Ru = 8314
wall = -20
n = 20
# Gas 1 properties
y1 = 1.4
P1 = 1
T1 = 300
mw1 = 28.96

# Gas 2 properties
y4 = 1.4
P4 = 5
T4 = 300
mw4 = 28.96

PR = P4/P1

a1 = (y1*Ru*T1/mw1)**0.5
a4 = (y4*Ru*T4/mw4)**0.5

# u3 = 100
x1 = [i/100 for i in range(0,1000)]
# z = []
i = 0
# Calculating Shock Strength

for x in x1:
    numerator = (y4 - 1) * (a1 / a4) * (x - 1)
    denominator = (2 * y1 * (2 * y1 + (y1 + 1) * (x - 1)))**0.5
    exponent = -2 * y4 / (y4 - 1)
    y = x * (1 - numerator / denominator) ** exponent
    if abs(y-PR) <= 0.02:
        Ps = x #Shock Strength
        break
    i+=1

P2 = Ps*P1 

#  Temperature ratio for shock
term1 = (P2 / P1) * (y1 + 1 + (P2*(y1-1)/ P1))
term2 = y1 - 1 + (y1 + 1) * (P2 / P1)
Ts = (term1 / term2)
T2 = Ts*T1


#  Temperature ratio for shock
Ds = Ps/Ts


# Shock Mach No.
Ms = ((y1 + 1) / (2*y1) * (P2 / P1 - 1) + 1)**0.5

# Shock Velocity 
W = a1*Ms

# Induced Velocity
up = W*(1 - 1/Ds)

# mach number behind shock
a2 = (y1*Ru*T2/mw1)**0.5
M2 = up/a2
# Reflected Shock Wave

u3 = up
u2 = u3
P3 = P2
term3 = (y4-1)/y4
T3 = T4*(P3/P4)**(term3)
a3 = (y1*Ru*T3/mw4)**0.5



def Jp(u,a):
    Jp = u + 2*a/(y4-1) 
    return Jp
def Jm(u,a):
    Jm = u - 2*a/(y4-1) 
    return Jm

u4 = 0

u_1 = u4
u_n = u3

size = int((n+3)*(n+2)/2)

step_size = (u_1 - u_n) / (n+1)
slope2 = 1/(u3-a3)
slope1 = 1/(u4-a4)

# Perform linear interpolation and store the results in a list
u = [u_1 - step_size * (i + 1) for i in range(-1,n+1)] 
uu = u+ [0]*(size-len(u))
a = [a4*(1-(y4-1)*x/(2*a4)) for x in u]
aa = a + [0]*(size-len(a))
slopes = [1/(u[i]-a[i]) for i in range(len(u))]

# boundary velocities 
b = []
k = 1
for i in range(n+2):
    b.append(k)
    k = k+n+2-i

ab = [-Jm(u[i],a[i])*(y4-1)/2 for i in range(len(u))]
U = []
for value, index in zip(ab, b):
    aa[index-1]= value
for index in b:
    uu[index-1]= 0

# JM = []
JpU = [Jp(uu[i],aa[i]) for i in range(len(uu))] + [0]*(size-len(uu))
JmU = [Jm(uu[i],aa[i]) for i in range(len(uu))] + [0]*(size-len(uu))
Jpb = [Jp(0,ab[i]) for i in range(len(ab))]
Jmb = [Jm(0,ab[i]) for i in range(len(ab))]
# JM = JM + JmU

l = [x for x in range(n+4,size+1)]
# Reflected point which are Inside the domain
rp = [x for x in l if x not in b] 
a = 0
b = n
for j in range(len(rp)):
    JpU[rp[j]-1] = JpU[rp[j]-1-1]
    if j>0: 
        if rp[j] - rp[j-1] == 2:
            a = a+1
            # print("Hello")
    JmU[rp[j]-1] = JmU[rp[j]-1-(n+1)+a]
    aa[rp[j]-1] = (y4-1)*(JpU[rp[j]-1]-JmU[rp[j]-1])/4
    uu[rp[j]-1] = (JpU[rp[j]-1]+JmU[rp[j]-1])/2
    
        


# Getting Powerset of all the lines
lin = [x for x in range(1, size)]
a = n+2
q = 0
p = []
f = []
pd = {}
pd[(0,1)] = slope1
for i in range(n+1):
    for j in range(2,a+1):
        # sl = tan(0.5(atan(1/(uu[j-1+q-1]-aa[j-1+q-1]) + 1/(uu[j+q-1]-aa[j+q-1]))))
        ee = [j-1+q,j+q]
        gg = [j+q,j+a-1+q]
        e = [(j-1+q,j+q),tan(0.5*atan(1/(uu[j-1+q-1]+aa[j-1+q-1])) + atan(1/(uu[j+q-1]+aa[j+q-1])))]
        g = [(j+q,j+a-1+q),tan(0.5*atan(1/(uu[j+q-1]-aa[j+q-1])) + atan(1/(uu[j+a-1+q-1]-aa[j+a-1+q-1])))]
        p.append(e)
        p.append(g)
        f.append(ee)
        f.append(gg)
        pd[tuple(ee)] = tan(0.5*atan(1/(uu[j-1+q-1]+aa[j-1+q-1])) + atan(1/(uu[j+q-1]+aa[j+q-1])))
        pd[tuple(gg)] = tan(0.5*atan(1/(uu[j+q-1]-aa[j+q-1])) + atan(1/(uu[j+a-1+q-1]-aa[j+a-1+q-1])))

    q = q+j
    a = a-1

n = n+2
p.append([0,0])
class Node:

    def __init__(self, value):
        self.value = value
        self.line1 = None
        self.line2 = None
        self.intercept1 = None
        self.intercept2 = None
        self.slope1 = None
        self.slope2 = None
        self.eqn1 = None
        self.eqn2 = None
        self.point = None
    
    # def set_inputs(self, line1, line2):
    #     self.line1 = line1
    #     self.line2 = line2

    def find_intersection(self):
        m1, b1 = self.line1
        m2, b2 = self.line2
        if m2 is not None:
            x_intersection = (b2 - b1) / (m1 - m2)
            y_intersection = m1 * x_intersection + b1

        elif m1 == m2:
            return None
        else:
            x_intersection = b2
            y_intersection = m1 * x_intersection + b1
        self.point = [x_intersection, y_intersection]
        return 
    
    
    def slope_intercept_equation(self):
    # Unpack the coordinates of the point
        x1, y1 = self.point
        
        # Calculate the y-intercept
        self.intercept1 = y1 - self.slope1 * x1
        self.intercept2 = y1 - self.slope2 * x1
        
    # Return the slope and intercept
        self.eqn1 = (self.slope1,self.intercept1)
        self.eqn2 = (self.slope2,self.intercept2)

        return 
bpt = []
spt = []
nodes = [Node(i) for i in range(1,int(n*(n+1)/2)+1)]

nodes[0].line1 = (slopes[0],0)
nodes[0].line2 = (None,wall)
nodes[0].slope1, nodes[0].slope2 = 0, p[0][1]
nodes[0].find_intersection()
nodes[0].slope_intercept_equation()
print(nodes[0].value,nodes[0].point)
spt.append(nodes[0].point)
k=1
i=1
for i in range(1,n):
    nodes[i].line1 = (slopes[i],0)
    nodes[i].line2 = nodes[i-1].eqn2
    nodes[i].slope1, nodes[i].slope2 = p[k][1], p[k+1][1]
    k+=2
    nodes[i].find_intersection()
    nodes[i].slope_intercept_equation()
    print(nodes[i].value,nodes[i].point)
    spt.append(nodes[i].point)
k = k-1
counter = 1
counter = counter+i
bpt.append([nodes[i].point,p[k-2][1]])

for s in range(n-1,0,-1):
    print("wall",counter,counter-s)
    nodes[counter].line1 = nodes[counter-s].eqn1
    nodes[counter].line2 = (None,wall)
    nodes[counter].slope1, nodes[counter].slope2 = 0, p[k][1]
    nodes[counter].find_intersection()
    nodes[counter].slope_intercept_equation()
    # print(nodes[counter].value,nodes[counter].point)
    counter+=1
    m = counter
    k = k+1
    for j in range(1,s):
        print("b/w",m+j-1,m+j-s-1)
        
        nodes[m+j-1].line1 = nodes[m+j-s-1].eqn1
        nodes[m+j-1].line2 = nodes[m+j-1-1].eqn2
        nodes[m+j-1].slope1, nodes[m+j-1].slope2 = p[k][1], p[k+1][1]
        # print(k,k+1)
        k+=2
        nodes[m+j-1].find_intersection()
        nodes[m+j-1].slope_intercept_equation()
        # print(nodes[m+j-1].value,nodes[m+j-1].point)
        # l+=1
        # bpt.append([nodes[m+j-1].value,p[k][1]])

        counter+=1
    bpt.append([nodes[counter-1].point,p[k-3][1]])
    k = k-1
ptss = [nodes[x].point for x in range(int(n*(n+1)/2))]

print(p)
print(bpt)
# for p in ptss:
#     plt.scatter(p[0],p[1])
    
wallalues, y_values = zip(*ptss)
for connection in f:
    plt.plot([wallalues[connection[0] - 1], wallalues[connection[1] - 1]],
             [y_values[connection[0] - 1], y_values[connection[1] - 1]],
             color='red', linestyle='-', linewidth=1)
    
for point, slope in bpt:
    # Unpack point coordinates
    x, y = point

    # Define the length of the line segment you want to draw
    line_length = abs(wall)  # Modify this as needed

    # Calculate the new point based on the slope and length
    new_x = x + line_length
    new_y = y + (slope * line_length)

    # Plotting the line segment
    if  new_y>=0:
        plt.plot([x, new_x], [y, new_y], color='red', linestyle='--', linewidth=1)

    # Plotting the given point
    # plt.scatter(x, y)
spt.append([0, 0])

# plt.figure(figsize=(8, 6))

# Plot lines from each point to the origin
for point in spt:
    # Unpack point coordinates
    x, y = point

    # Plotting the line segment to the origin
    plt.plot([0, x], [0, y], color='blue', linestyle='-', linewidth=1)

    # Plotting the given point
    # plt.scatter(x, y)
plt.title('Reflection Of Expansion Fans in Shock Tube')
plt.ylabel('T(s)')
plt.xlabel('Position in Shock Tube(m)')
# plt.legend()
# plt.grid(True)

# Display the plot
plt.grid()
# plt.legend()
plt.show()
print("spt",spt)