import matplotlib.pyplot as plt
import numpy as np
import math

def y_prime(t,y):
    return y - (4 * t * t) + 1

def exact_solution(t):
    return (4*t*t) + (8*t) + 7 - (6*np.exp(t))

def euler(a,b,f,alpha,h):
    big_n = (b-a)/h
    t = a
    omega = alpha
    print("*********Euler's Method***********,")
    print(exact_solution(t)- omega, ",")
    for i in range(1,int(big_n)+1):
        omega = omega + (h * f(t,omega))
        t = i*h
        print(exact_solution(t)- omega)

def modified_euler(a,b,f,alpha,h):
    print("*********Modified Euler's Method***********")
    big_n = (b-a)/h
    t = a
    omega = alpha
    print(exact_solution(t)- omega, ",")
    for i in range(1,int(big_n)+1):
        omega = omega + h/2*(f(t,omega) + f(t+h,omega + (h*f(t,omega))))
        t = i*h
        print(exact_solution(t) - omega, ",")

def rk_four(a,b,f,alpha,h):
    print("*********Rk Four ***********")
    big_n = (b-a)/h
    omega = alpha
    t = a
    print(exact_solution(t)- omega, ",")
    for i in range(1, int(big_n)+1):
        k1 = h*f(t,omega)
        k2 = h*f(t+(h/2),omega+(k1/2))
        k3 = h*f(t+(h/2),omega+(k2/2))
        k4 = h*f(t+h,omega+k3)
        omega = omega + 1/6*(k1+(2*k2)+(2*k3)+k4)
        t = i*h
        print(exact_solution(t) - omega, ",")

def main():
    euler(0,1,y_prime,1,.05)
    modified_euler(0,1,y_prime,1,.05)
    rk_four(0,1,y_prime,1,.05)
if __name__ == "__main__":
    main()
