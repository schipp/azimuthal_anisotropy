def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    if D != 0:
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        x = Dx / D
        y = Dy / D
        return [x, y]
    else:
        return []

# def get_t(x1, x2, x3, x4, y1, y2, y3, y4):
#     return ((x1 - x3) * (y1 - y3) - (y2 - y1) * (x1 - x3))/((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
# 
# def get_u(x1, x2, x3, x4, y1, y2, y3, y4):
#     return ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3))/((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

def get_intersection_point(line1, line2):
    x1, x2, x3, x4 = line1[0][1], line1[1][1], line2[0][1], line2[1][1]
    y1, y2, y3, y4 = line1[0][0], line1[1][0], line2[0][0], line2[1][0]
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3))/((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    P_a = (y1 + t * (y2 - y1), x1 + t*(x2 - x1))
    P_b = (y3 + u * (y4 - y3), x3 + u*(x4 - x3))
    if t > u and 0 <= t <= 1:
        return P_a
    if t < u and 0 <= u <= 1:
        return P_b
    return []
import numpy as np

# a piece of a prolate cycloid, and am going to find
a, b = 1, 2
phi = np.linspace(3, 10, 100)
x1 = a*phi - b*np.sin(phi)
y1 = a - b*np.cos(phi)

x2=phi
y2=np.sin(phi)+2
# x,y=intersection(x1,y1,x2,y2)
# plt.plot(x1,y1,c='r')
# plt.plot(x2,y2,c='g')
# plt.plot(x,y,'*k')
# plt.show()


# print(line1)
# print(line2)
L1 = line([-50, -2], [2, 2])
L2 = line([-1, 0], [1, 0])
# print(intersection(L1, L2))
line1 = ([-2, -2], [2, 2])
line2 = ([-3, 0], [1, 0])

# line1 = [[53.5, 1.5], [54.5, 1.5]]
# line2 = [[52.4759, 13.78], [47.2505, 15.5232]]


P = get_intersection_point(line1, line2)
# print(P_a, P_b, t, u)

import pylab as plt
import numpy as np
fig, ax = plt.subplots(1, 1)
ax.plot(np.array(line1)[:, 1], np.array(line1)[:, 0])
ax.plot(np.array(line2)[:, 1], np.array(line2)[:, 0])
# ax.scatter(P[1], P[0])
P2 = intersection(np.array(line1)[:, 1],np.array(line1)[:, 0],np.array(line2)[:, 1],np.array(line2)[:, 0])
print(P2)
ax.scatter(P2[1], P2[0])
fig.savefig('intersect.png')
plt.close(fig)