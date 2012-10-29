
TYPE_DECS = ('TYPE_DECS', 
             #(func_name, [input types], return type)
             ['compute_estimate', [('array', 'double')], 'double'],             
             ['reduce_bootstraps', [('array', 'double')], 'double'],
             ['average', [('array', 'double')], 'double']
            )


DIM = 1

## data passed in is a numpy array and then figure out DIM
## translate for i in range(len(data))

def compute_estimate(data):
    avg = 0.0
    if len(data)==0: return avg
    else:
        for i in range(len(data)):
            avg += data[i]
        return avg / len(data)
    #return stddev(data)

def reduce_bootstraps(bootstraps):
    mean = 0.0
    for bootstrap in bootstraps:
        mean += bootstrap
    return mean / len(bootstraps)

def average(subsamples):
    mean = 0.0
    for subsample in subsamples:
        mean += subsample
    return mean/ len (subsamples)

