from subprocess import *
times = ["9e-6","8e-6","7e-6","6e-6","5e-6","4e-6","3e-6","2e-6","1e-6","9e-7","8e-7","7e-7","6e-7","5e-7","4e-7","3e-7","2e-7","1e-7","9e-8","8e-8","7e-8","6e-8","5e-8","4e-8","3e-8","2e-8","1e-8","9e-9","8e-9","7e-9","6e-9"]#,"5e-9"]
for i in range(len(times)):
    # print("Time Step: " + times[i])
    x = times[i]
    args = ["0.01", "100", x, "0", "0", "0", "0", "0", "0", "4" , "3", "0", "0", "0", "0", "0", "5" , "3", "4", "0", "0", "0", "0", "3"]
    argsString = ' '.join(args)
    call("./assignment " + argsString, shell=True)
