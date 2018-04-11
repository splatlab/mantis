import socket
import sys

HOST, PORT = "localhost", 23901

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

s.connect((HOST, PORT))

query_filepath = "/home/paf2023/run_mantis/queries/query{0}.fa"
output_filepath = "/home/paf2023/run_mantis/out/out{0}.res"

for x in range(0, 5):
    print("Send request {0}".format(query_filepath.format(x)))
    s.send(query_filepath.format(x) + " " + output_filepath.format(x))
    result = str(s.recv(4096))
    print("Got result {0}".format(result))
