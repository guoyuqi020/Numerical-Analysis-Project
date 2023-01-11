import os

examples = [
    [10000, 100, 315],
    [50000, 100, 315],
    [100000, 100, 315],
    [10000, 100, 500],
    [50000, 100, 500],
    [100000, 100, 500],
    [10000, 100, 1000],
    [50000, 100, 1000],
    [100000, 100, 1000],
    [100000, 120, 380],
    [100000, 120, 500],
    [100000, 120, 1000],
    [100000, 150, 475],
    [100000, 150, 500],
    [100000, 150, 1000],
    [100000, 200, 630],
    [100000, 200, 1000],
    [100000, 200, 2000]
]
package = "atlas"
idx = 0
for item in examples:
    os.system("./bin/main -M {} -k {} -N {} -r ./rawData/{} -f ./results/{} > ./shellOutput/{}/shell{}.txt".format(
        item[0], item[1], item[2], package, package, package, idx))
    os.system("python3 ./plot.py {} {} {} {} ./figures {}".format(idx,
              item[0], item[1], item[2], package))
    idx += 1
