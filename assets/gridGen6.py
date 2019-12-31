from subprocess import call

dimension = 20
# maxThreads = 8
comm = ["time", "../assignment3", "0.01", "1", "0.01"]
for i in range(0,dimension):
    for j in range(0,dimension):
        for m in range(0, dimension):
            comm.append(str(i))
            comm.append(str(j))
            comm.append(str(m))
            comm.append("0")
            comm.append("0")
            comm.append("0")
            comm.append("1")
print("SOLS1")
comm[1]="../assignment3"
call(comm)
print("SOLS2")
# for i in range(1,maxThreads):
#     print("cores:" + str(i))
call("export OMP_NUM_THREADS=8", shell=True)
comm[1]="../assignment6"
call(comm)
# file = open("text_string.txt", "w")
# file.write(string)
# file.close()
