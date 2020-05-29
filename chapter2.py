import numpy as np
def Modular_square(a,b,c):
    res = 1
    temp = a%c
    rod = b
    while(rod>0):
        if(rod%2==1):
            res = (res*temp)%c
        temp = (temp**2)%c
        print(res," ",temp)
        rod = np.floor(rod/2)
if __name__== '__main__':
    Modular_square(12996,227,37909)
