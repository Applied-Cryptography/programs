import numpy as np
#模平方乘法的过程
def Modular_square(a,b,c):
    res = 1
    temp = a%c
    rod = b
    while(rod > 0):
        if(rod%2 == 1):
            res = (res * temp) % c
        temp = (temp**2) % c
        print(rod%2," ",res, " ", temp)
        rod = np.floor(rod/2)

#求逆的过程
def exgcd(a,b):
    if(a<b):
        temp = a
        a = b
        b =temp
    if(b==0):
        return 1,0,a
    print(a,"=",np.floor(a/b),"*",b,"+",a%b)
    if(a%b==0):
        print("最大公因数是",b,"下面是求逆的过程")
    c,d,e = exgcd(b,a%b)
    k = c-np.floor(a/b)*d
    if(k<0):
        print(e,"=",d,"*",a,k,"*",b)
    else:
        print(e, "=", d, "*", a, "+", k, "*", b)
    return d,k,e

if __name__== '__main__':
    #Modular_square(12996,227,37909)
    a,b,c=exgcd(44350,20785)