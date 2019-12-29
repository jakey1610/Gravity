import os

ts = "1e-8"
string = "./assignment 0.01 100.0 "+ts+" 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3"

point = os.popen(string).read()
print(point)
