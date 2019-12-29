from subprocess import call

dimension = 30
comm = ["time", "../assignment", "0.01", "1", "0.01"]
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

call(comm)
comm[1]="../assignmentParallel"
call(comm)
# file = open("text_string.txt", "w")
# file.write(string)
# file.close()
