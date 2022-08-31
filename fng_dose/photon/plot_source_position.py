import matplotlib.pyplot as plt
import openmc.lib
import numpy as np

openmc.lib.init()
source = openmc.lib.sample_external_source(10000)

x = np.array([p.r[0] for p in source])
y = np.array([p.r[1] for p in source])
z = np.array([p.r[2] for p in source])

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(x, y, z)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

basis = 'yz'

if basis == 'xy':
    plt.hist2d(x, y, 100, [[-50., 50.], [-14., 77.]])
elif basis == 'yz':
    plt.hist2d(y, z, 100, [[-14., 77.], [-50., 50.]])
plt.show()

openmc.lib.finalize()
