import matplotlib.pyplot as plt


k1 = 100e-3
k2 = 600e-3
k3 = 150e-3


Earray = []
Sarray = []
ESarray = []
Parray = []
tarray = []
Varray = []


'''
Four equations for the rate of changes of the four species.
'''
def f(E, S, ES, P):
    a = (k2 + k3) * ES - k1 * E * S
    return a


def g(E, S, ES, P):
    a = k2 * ES - k1 * E * S
    return a

def m(E, S, ES, P):
    a = k1 * E * S - (k2 + k3) * ES
    return a

def p(E, S, ES, P):
    a = k3 * ES
    return a

'''
Fourth-order Runge-Kutta method
'''
def RK4():
    h = 1
    t = 0
    E = 1e-3
    S = 10e-3
    ES = 0
    P = 0

    while t <= 10000000:
        tarray.append(t)
        Earray.append(1000*E)
        Sarray.append(1000*S)
        ESarray.append(1000*ES)
        Parray.append(1000*P)
        t = t + h

        E1 = f(E, S, ES, P)  # Step1
        m1 = E + E1 * h / 2
        S1 = g(E, S, ES, P)
        n1 = S + S1 * h / 2
        ES1 = m(E, S, ES, P)
        x1 = ES + ES1 * h / 2
        P1 = p(E, S, ES, P)
        y1 = P + P1 * h / 2

        E2 = f(E1, S1, ES1, P1)  # Step2
        m2 = E + E2 * h / 2
        S2 = g(E1, S1, ES1, P1)
        n2 = S + S2 * h / 2
        ES2 = m(E1, S1, ES1, P1)
        x2 = ES + ES2 * h / 2
        P2 = p(E1, S1, ES1, P1)
        y2 = P + P2 * h / 2

        E3 = f(E2, S2, ES2, P2)  # Step3
        m3 = E + E3 * h
        S3 = g(E2, S2, ES2, P2)
        n3 = S + S3 * h
        ES3 = m(E2, S2, ES2, P2)
        x3 = ES + ES3 * h
        P3 = p(E2, S2, ES2, P2)
        y3 = P + P3 * h

        E4 = f(E3, S3, ES3, P3)  # Step4
        S4 = g(E3, S3, ES3, P3)
        ES4 = m(E3, S3, ES3, P3)
        P4 = p(E3, S3, ES3, P3)

        E = E + (E1 + 2 * E2 + 2 * E3 + E4) * h / 6
        S = S + (S1 + 2 * S2 + 2 * S3 + S4) * h / 6
        ES = ES + (ES1 + 2 * ES2 + 2 * ES3 + ES4) * h / 6
        P = P + (P1 + 2 * P2 + 2 * P3 + P4) * h / 6



def main():
    RK4()
    '''
    It is possible to approximate the rate of product formation.
    Definitions:
    [E]0 the enzyme initial concentration
    [S] the substrate concentration
     k3 the rate constant for product
     M the Michaelis constant. And M=(k2+k3)/k1.
     The Michaelis-Menten equation gives velocity V a function of concentration of S:
     V=k3*[E]0*[S]/(M+[S])
    '''
    E0 = 1
    M = (k2 + k3) / k1
    for i in Sarray:
        V = 1000 * k3 * E0 * i / (M + i)
        Varray.append(V)

    plt.plot(Sarray,Varray)
    plt.xlabel("Concentration S")
    plt.ylabel("V(μM/min)")
    plt.show()


if __name__ == "__main__":
    main()