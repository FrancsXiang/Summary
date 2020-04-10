import numpy as np
import matplotlib.pyplot as plt

def cubic_fc(s, a, b, c, d, x):
    return a + b * (x - s) + c * (x - s) * (x - s) + d * (x - s) * (x - s) * (x - s)

with open('./data.txt', 'r') as file:
    res = np.loadtxt(file)
    start = res[0][0]
    end = res[-1][1]
    size = []
    x = np.linspace(start, end, 1000)
    x = list(x)
    # for i in range(len(res)):
    #     x.append(res[i][0])
    # x.append(res[-1][1])
    x.sort()
    y = [0 for i in range(len(x))]
    for i in range(len(res)):
        start = res[i][0]
        end = res[i][1]
        a = res[i][2]
        b = res[i][3]
        c = res[i][4]
        d = res[i][5]
        cond = [1 if (i >= start and i <= end) else 0 for i in x]
        y += cond * cubic_fc(start, a, b, c, d, x)
        y_ = [np.sin(x[i]) for i in range(len(x))]
    plt.plot(x, y, "x-",label="prediction")
    plt.plot(x, y_, "+-",label="ground_truth")
    plt.legend()
    plt.grid(True)
    plt.show()
