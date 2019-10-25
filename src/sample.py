import numpy as np;
import globimap as gm;
from matplotlib import pyplot as plt;
from tqdm import tqdm;

class cfg:
    k=15
    logm=23
    low=[47,11] # where to put this into infinite space


def generate_sierpinski (size=4096, steps=500000):
    m = np.zeros(size*size).reshape(size,size);
    print(m.shape)
    triangle = [[0,0],[1,0],[0.5,1]];
    s = np.random.rand(2)
    for _ in tqdm(range(steps)):
        p = triangle[np.random.choice(range(3))] 
        s = s + (p-s)/2
        c = (s * size).astype(np.uint64)
        m[c[0],c[1]] = 1
    return (m)
   
    
    

if __name__=="__main__":
    print("Memory: %f KB" % (2**cfg.logm / 8.0 / 1024.0 ))
    print("Hashes: %d " % (cfg.k))
    
    patch = generate_sierpinski()
    plt.subplot(131)
    plt.imshow(patch)
    print ("Ones: %d" % (np.sum(patch==1)))
    print ("Coordinate Memory (4 bytes): %d" %(np.sum(patch==1)*4))
    print ("Raster Memory (1 bit per pixel): %d KB" %(4096*4096/8/1024))

    
    m = gm.globimap()
    m.configure(cfg.k,cfg.logm)
    print("Step 1: Probabilistic Encoding")
    m.map(patch, cfg.low[0],cfg.low[1])

    # now add the error correction information
    print("Step 2: Error Correction Information")
    m.enforce(patch, cfg.low[0],cfg.low[1])
    print(m.summary());

    print("Step 3: Rasterize")
    result = m.rasterize(cfg.low[0], cfg.low[1], patch.shape[0], patch.shape[1])
    plt.subplot(132)
    plt.imshow(result)
    
    result = m.correct(cfg.low[0], cfg.low[1], patch.shape[0], patch.shape[1])
    plt.subplot(133)
    plt.imshow(result)
    plt.show()
    nerror = np.sum(result != patch)
    ber = float(nerror) / np.prod(patch.shape)
    print("Have %d errors for a BER of %f" % (nerror, ber))

