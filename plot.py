import numpy as np
import matplotlib.pyplot as plt
import sys

rawX = []
rawY = []
fitArguments = []


def calculateFitY(x, arguments):
    n = len(arguments)
    y = 0
    for i in range(n):
        y += arguments[i]*(x**i)
    return y


if __name__ == "__main__":
    with open("rawData/"+sys.argv[6] +"/rawdata"+sys.argv[1]+"_x.txt") as f:
        for line in f:
            rawX.append(float(line))
    with open("rawData/"+sys.argv[6] +"/rawdata"+sys.argv[1]+"_y.txt") as f:
        for line in f:
            rawY.append(float(line))
    with open("results/"+sys.argv[6] +"/fitresults"+sys.argv[1]+".txt") as f:
        for line in f:
            fitArguments.append(float(line))
    M = int(sys.argv[2])
    k = int(sys.argv[3])
    N = int(sys.argv[4])
    x = np.linspace(0, 1.0, num=20000)
    plt.plot(x, calculateFitY(x, fitArguments), label="polyfit")
    plt.plot(rawX, rawY, label="observation")
    plt.figlegend()
    plt.title("M={} k={} N={}".format(M, k, N))
    plt.draw()
    plt.savefig(
        "{}/{}/figure{}_M_{}_k_{}_N_{}.png".format(sys.argv[5], sys.argv[6], sys.argv[1], M, k, N))
