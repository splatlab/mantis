import socket
import sys

HOST, PORT = "localhost", 23901

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

s.connect((HOST, PORT))
for x in range(0, 5):
    print("Step 1")
    s.send('/home/paf2023/run_mantis/queries/read2.fa ./out.res')
    print("Step 2")
    print(str(s.recv(1000)))
    print(x)
