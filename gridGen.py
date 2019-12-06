from subprocess import call
import time

dimension = 10
comm = ["./assignment", "0.01", "0.5", "0.01"]
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
start = int(round(time.time() * 1000))
call(comm)
end = int(round(time.time() * 1000))
print("Time taken: " + str(end-start))
# file = open("text_string.txt", "w")
# file.write(string)
# file.close()
