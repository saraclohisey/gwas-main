allchroms = [1,2,3.0,"X"]

for c in allchroms:
    try:
        c = int(c)
    except:
        pass
    print (c)