import re
from numpy import mat

from soupsieve import match

test1 = "12895a121-w1"
test2 = "12895a121-w2"
test3 = "12895b121-w1"
test4 = "12895b121-w2"

test_list = [test1, test2, test3, test4]
regex = r"[0-9]{5}[a]"

matches = []
for test in test_list:
    if len(re.findall(regex, test)) == 0:
        pass
    else:
        matches.append(test)
print(matches)