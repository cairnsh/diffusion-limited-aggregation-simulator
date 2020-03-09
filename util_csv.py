from matplotlib import pyplot as p

def read_from_csv(f):
    o = open(f)
    data = []
    for line in o:
        data += [map(int, DELIMITER.split(line.strip()))]
    return data
