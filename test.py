lt = ["12412b120-w2", "12412b120-w1", "12412b123-w2", "12412b123-w1", "12412b125-w1", "12412b125-w2"]

idx = list(range(len(lt)))[::2]

sorted = []
for i in idx:
    pair = (lt[i], lt[i + 1])
    if pair[0][10:12] == "w1":
        sorted.extend([pair[0], pair[1]])
    else:
        sorted.extend([pair[1], pair[0]])

print(sorted)
    