
pip install deepxde
import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

# Predicted parameter
k12 = 0.744668
k21 = 1.5234107
ka = 2.3669796
Vm = 0.8679799
V1 = 0.94791704
Km = 0.8873355

def ode_system(x, y):
    # Define ODE system
    y1, y2, y3 = y[:, 0:1], y[:, 1:2], y[:, 2:]
    dy1_x = dde.grad.jacobian(y, x, i=0)
    dy2_x = dde.grad.jacobian(y, x, i=1)
    dy3_x = dde.grad.jacobian(y, x, i=2)
    return [
        dy1_x - k21 * y2 + k12 * y1 - ka * y3 + (Vm * y1)/(V1 * (Km + y1)),
        dy2_x + k21 * y2 - k12 * y1,
        dy3_x + ka * y3,
    ]

def boundary(_, on_initial):
    return on_initial

if __name__ == "__main__":
    geom = dde.geometry.TimeDomain(0, 12)
    ic1 = dde.icbc.IC(geom, lambda X: 0, boundary, component=0)
    ic2 = dde.icbc.IC(geom, lambda X: 0, boundary, component=1)
    ic3 = dde.icbc.IC(geom, lambda X: 0.48048388, boundary, component=2)

    data = dde.data.PDE(
        geom,
        ode_system,
        [ic1, ic2, ic3],
        num_domain=35,
        num_boundary=3,
        solution=None,
        num_test=150,
    )

    layer_size = [1] + [50] * 3 + [3]
    activation = "tanh"
    initializer = "Glorot uniform"
    net = dde.nn.FNN(layer_size, activation, initializer)

    model = dde.Model(data, net)
    model.compile("adam", lr=0.001)
    losshistory, train_state = model.train(iterations=150000)

    dde.saveplot(losshistory, train_state, issave=True, isplot=True)
    plt.figure(figsize=(10, 6))
