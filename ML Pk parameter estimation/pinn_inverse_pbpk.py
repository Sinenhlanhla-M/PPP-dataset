
pip install deepxde
import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py

# Define the parameters we want to learn. This is for the inverse problem.
k12 = dde.Variable(1.0)
k21 = dde.Variable(1.0)
ka = dde.Variable(1.0)
Vm = dde.Variable(1.0)
V1 = dde.Variable(1.0)
Km = dde.Variable(1.0)

def import_hdf5(filename):
    # The vector numerical data from mathematica is best read in in .h5 format. This function does that.
    f = h5py.File(filename, 'r')
    a_group_key = list(f.keys())[0]
    data = list(f[a_group_key])
    return np.array(data)

def PBBK_3_system(x, y):
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

def fix_path(path):
    return './'+'/'.join(path.split('\\')[-2:])

def gen_traindata_mm():
    trainingData = pd.read_csv('/content/drive/My Drive/sne/PBPK.csv')
    trainingData = trainingData.iloc[75]
    ts = import_hdf5(fix_path(trainingData['t']))
    sols = import_hdf5(fix_path(trainingData['sols']))
    return np.array([ts]).T, sols.T, trainingData['dose'], list(trainingData[['k12', 'k21', 'ka', 'Vm', 'V1', 'Km']])

# Get the train data
observe_t, ob_y, dose, true_parms = gen_traindata_mm()
observe_y0 = dde.icbc.PointSetBC(observe_t, ob_y[:, 0:1], component=0)
observe_y1 = dde.icbc.PointSetBC(observe_t, ob_y[:, 1:2], component=1)
observe_y2 = dde.icbc.PointSetBC(observe_t, ob_y[:, 2:3], component=2)

geom = dde.geometry.TimeDomain(0, observe_t[-1, 0])

# Initial conditions
ic1 = dde.icbc.IC(geom, lambda X: 0, boundary, component=0)
ic2 = dde.icbc.IC(geom, lambda X: 0, boundary, component=1)
ic3 = dde.icbc.IC(geom, lambda X: 1, boundary, component=2)

data = dde.data.PDE(
    geom,
    PBBK_3_system,
    [ic1, ic2, ic3, observe_y0, observe_y1, observe_y2],
    num_domain = 100,
    num_boundary = 2,
    anchors = observe_t,
)

net = dde.nn.FNN([1] + [40] * 3 + [3], "tanh", "Glorot uniform")
model = dde.Model(data, net)

external_trainable_variables = [k12, k21, ka, Vm, V1, Km]
variable = dde.callbacks.VariableValue(
    external_trainable_variables, period=600
)

model.compile(
    "adam", lr=0.001, external_trainable_variables=external_trainable_variables
)
losshistory, train_state = model.train(iterations=550000, callbacks=[variable])

model.compile("L-BFGS", external_trainable_variables=external_trainable_variables)
losshistory, train_state = model.train(iterations=550000, callbacks=[variable])

dde.saveplot(losshistory, train_state, issave=True, isplot=True)

plt.plot(observe_t, ob_y[:, 0:1], label='C1-Observed', color='red')
plt.plot(observe_t, ob_y[:, 1:2], label='C2-Observed', color='purple')
plt.plot(observe_t, ob_y[:, 2:3], label='C3-Observed', color='cyan')

idx = np.argsort(train_state.X_test[:, 0])
plt.plot(train_state.X_test[idx, 0], train_state.best_y[idx, :1], label='C1-Predicted', color='blue')
plt.plot(train_state.X_test[idx, 0], train_state.best_y[idx, 1:2], label='C2-Predicted', color='green')
plt.plot(train_state.X_test[idx, 0], train_state.best_y[idx, 2:3], label='C3-Predicted', color='yellow')

print(variable.get_value())
print(true_parms)
print(f'L2 error between predicted parameters and true parameters: {np.linalg.norm(np.asarray(variable.get_value())-true_parms)}')

plt.xlabel('Time (Hour)')
plt.ylabel('Concentration (mg/L)')
plt.legend()
plt.figure(figsize=(10, 6))
