import time
import string
import random

def intro():

    goal = "Welcome to the FRAGMENT GENERATOR"

    prog = ""

    for chr in goal:

        for i in range(0, 10):
            print(prog + random.choice(string.ascii_letters), end = "\r")
            time.sleep(1/240)



        prog = prog + chr
        print(prog, end = "\r")
        time.sleep(1/240)

    print(goal)
