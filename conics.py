import numpy as np
import math

def setEquation(equation, a, b, c, d, e, f, gcd, var1, var2):
    equation += f"{round(a,2)}{var1}²" if a != 0 else ""
    
    if gcd <= 0:
        gcd = 1

    if a != 0:
        if b > 0:
            equation += f" + {round(b/gcd,2)}{var1}{var2}"
        elif b < 0:
            equation += f" - {-round(b/gcd,2)}{var1}{var2}"
    elif b != 0:
        equation += f"{round(b/gcd,2)}{var1}{var2}"
    
    if a != 0 or b != 0:
        if c > 0:
            equation += f" + {round(c/gcd,2)}{var2}²"
        elif c < 0:
            equation += f" - {-round(c/gcd,2)}{var2}²"
    elif c != 0:
        equation += f"{round(c/gcd,2)}{var2}²"
    
    if a != 0 or b != 0 or c != 0:
        if d > 0:
            equation += f" + {round(d/gcd,2)}{var1}"
        elif d < 0:
            equation += f" - {-round(d/gcd,2)}{var1}"
    elif d != 0:
        equation += f"{round(d/gcd,2)}{var1}" 
    
    if a != 0 or b != 0 or c != 0 or d != 0:
        if e > 0:
            equation += f" + {round(e/gcd,2)}{var2}"
        elif e < 0:
            equation += f" - {-round(e/gcd,2)}{var2}"
    elif e != 0:
        equation += f"{round(e/gcd,2)}{var2}"

    if a != 0 or b != 0 or c != 0 or d != 0 or e != 0:
        if f > 0:
            equation += f" + {round(f/gcd,2)}"
        elif f < 0:
            equation += f" - {-round(f/gcd,2)}"
    elif f != 0:
        equation += f"{round(f/gcd,2)}"
    
    equation += " = 0"

    return equation

def setMatrix(matrix, a, b, c, d, e, f):
    matrix[0,0] = a
    matrix[1,1] = c
    matrix[2,2] = f

    for i in range(3):
        for j in range(3):
            if i != j:
                if i+j == 1:
                    matrix[i,j] = b/2
                elif i+j == 2:
                    matrix[i,j] = d/2
                else:
                    matrix[i,j] = e/2

def calcDet(matrix):
    positive = matrix[0,0] * matrix[1,1]
    negative = -(matrix[0,1] * matrix[1,0])

    det = positive + negative

    return det

def gcd(x, y):
    while(y):
        x,y = y,x % y
    return x

def simplifyParabola(a, c, d, e, f):
    if a == 0:
        if e == 0:
            t = "x = u"
            w = "y = v"
            
            f = round(f/d,2)

            if (d > 0 and f < 0) or (d < 0 and f > 0):
                t += f" + {round(abs(f),2)}"
                u = abs(f)
            else:
                t += f" - {round(abs(f),2)}"
                u = -abs(f)
            if d > 0:
                final = f"{round(c,2)}v² + {round(d,2)}u = 0"
            else:
                final = f"{round(c,2)}v² - {-round(d,2)}u = 0"

            print("{}\nCom:\n{}\n{}".format(final,t,w))

            return 0,c,d,0,u,0
        else:
            indU = (f-((e**2)/(4*(c**2))))/d
            indV = e/(2*c)

            t = "t = u"
            w = "w = v"

            if indU > 0:
                t += f" - {round(indU,2)}"
            elif indU < 0:
                t += f" + {-round(indU,2)}"
            if indV > 0:
                w += f" - {round(indV,2)}"
            elif indV < 0:
                w += f" + {-round(indV,2)}"
            
            indU = -indU
            indV = -indV

            final = f"{round(c,2)}v²"
            if d < 0:
                final += f" - {-round(d,2)}u = 0"
            else:
                final += f" + {round(d,2)}u = 0"
            
            print("{}\nCom:\n{}\n{}".format(final,t,w))
            
            return 0,c,d,0,indU,indV
    
    elif c == 0:
        if d == 0:
            t = "x = u"
            w = "y = v"
            
            f = round(f/e,2)

            if (e > 0 and f < 0) or (e < 0 and f > 0):
                w += f" + {round(abs(f),2)}"
                v = abs(f)
            else:
                w += f" - {round(abs(f),2)}"
                v = -abs(f)
            if e > 0:
                final = f"{round(a,2)}u² + {round(e,2)}v = 0"
            else:
                final = f"{round(a,2)}u² - {-round(e,2)}v = 0"
            
            print("{}\nCom:\n{}\n{}".format(final,t,w))

            return a,0,0,e,0,v
        else:
            indU = d/(2*a)
            indV = (f-((d**2)/(4*a**2)))/e

            t = "t = u"
            w = "w = v"

            if indU > 0:
                t += f" - {round(indU,2)}"
            elif indU < 0:
                t += f" + {-round(indU,2)}"
            if indV > 0:
                w += f" - {round(indV,2)}"
            elif indV < 0:
                w += f" + {-round(indV,2)}"
            
            indU = -indU
            indV = -indV
            
            final = f"{round(a,2)}u²"
            if e < 0:
                final += f" - {-round(e,2)}v = 0"
            else:
                final += f" + {round(e,2)}v = 0"
            
            print("{}\nCom:\n{}\n{}".format(final,t,w))
            
            return a,0,0,e,indU,indV

def translation(matrix, a, b, c, d, e, f):
    det = calcDet(matrix)
    if det == 0:
        print("\nTranslação impossível.")
        return a,b,c,d,e,f,0,0
    
    A = np.array([[matrix[0,0], matrix[0,1]], [matrix[1,0], matrix[1,1]]])
    B = np.array([-matrix[0,2], -matrix[1,2]])
    h,k = np.linalg.inv(A).dot(B)

    d = 0
    e = 0
    f = (matrix[2,0]*h) + (matrix[2,1]*k) + matrix[2,2]

    conicGCD = gcd(gcd(gcd(a, b), c), f)
    transl = ""
    transl = setEquation(transl, a, b, c, d, e, f, conicGCD, 'u', 'v')
    print("\nApós translação:  {}".format(transl))
    print(f"Com valores de h e k, respectivamente: {round(h,2)} e {round(k,2)}")

    if conicGCD <= 0:
        return a,b,c,d,e,f,h,k
    else:
        return a/conicGCD,b/conicGCD,c/conicGCD,d,e,f/conicGCD,h,k

def rotation(a, b, c, d, e, f):
    cotg2 = (a-c)/b
    inSin2 = math.sqrt(1+cotg2**2)

    firstEq = a+c
    secEq = b*inSin2

    a = (firstEq+secEq)/2
    c = firstEq-a
    b = 0
    
    sin2 = 1/inSin2
    cos2 = cotg2*sin2

    cos = math.sqrt((1+cos2)/2)
    sin = math.sqrt(1-cos**2)

    if d or e:
        aux = d
        d = (d*cos) + (e*sin)
        e = (aux*(-sin)) + (e*cos)

    conicGCD = gcd(gcd(gcd(gcd(gcd(a, b), c), d), e), f)
    rot = ""
    rot = setEquation(rot, a, b, c, d, e, f, conicGCD, 't', 'w')
    print("\nApós rotação:     {}".format(rot))
    print(f"Com valores do cosseno e do seno do ângulo de rotação, respectivamente: {round(cos,2)} e {round(sin,2)}")
    
    if conicGCD <= 0:
        return a,b,c,d,e,f,sin,cos
    else:
        return a/conicGCD,b/conicGCD,c/conicGCD,d/conicGCD,e/conicGCD,f/conicGCD,sin,cos

def changeCoord(x, y, h, k, a, b, c, sin, cos):
    if h or k:
        if sin or cos:
            if x:
                xc1 = h + c*cos
                yc1 = k + c*sin
                
                xc2 = h - c*cos
                yc2 = k - c*sin

                xa1 = h + a*cos
                ya1 = k + a*sin
                
                xa2 = h - a*cos
                ya2 = k - a*sin

                xb1 = h - b*sin
                yb1 = k + b*cos
                
                xb2 = h + b*sin
                yb2 = k - b*cos
            elif y:
                xc1 = h - c*sin
                yc1 = k + c*cos
                
                xc2 = h + c*sin
                yc2 = k - c*cos
                
                xa1 = h - a*sin
                ya1 = k + a*cos
                
                xa2 = h + a*sin
                ya2 = k - a*cos

                xb1 = h + b*cos
                yb1 = k + b*sin
                
                xb1 = h - b*cos
                yb1 = k - b*sin
        else:
            xc1 = h + c
            yc1 = k + c
            
            xc2 = h - c
            yc2 = k - c
            
            xa1 = h + a
            ya1 = k + a
            
            xa2 = h - a
            ya2 = k - a

            xb1 = h + b
            yb1 = k + b
            
            xb2 = h - b
            yb2 = k - b
    else:
        if x:
            xc1,yc1,xa1,ya1,xb1,yb1 = c,0,a,0,0,b
            xc2,yc2,xa2,ya2,xb2,yb2 = -c,0,-a,0,0,-b
        elif y:
            xc1,yc1,xa1,ya1,xb1,yb1 = 0,c,0,a,b,0
            xc2,yc2,xa2,ya2,xb2,yb2 = 0,-c,0,-a,-b,0

    return xc1,yc1,xa1,ya1,xb1,yb1,xc2,yc2,xa2,ya2,xb2,yb2

def ellipse(a, c, f, h, k, sin, cos):
    print("\nA cônica é uma elipse.")

    a2 = -f/a
    b2 = -f/c
    
    x = False
    y = False

    if a2 < b2:
        a2,b2 = b2,a2
        y = True
    else:
        x = True

    newA = math.sqrt(a2)
    newB = math.sqrt(b2)
    newC = math.sqrt(a2-b2)
    e = newC/newA

    xc1,yc1,xa1,ya1,xb1,yb1,xc2,yc2,xa2,ya2,xb2,yb2 = changeCoord(x, y, h, k, newA, newB, newC, sin, cos)

    if x:
        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Focos: F₁ = ({round(xc1,2)},{round(yc1,2)}); F₂ = ({round(xc2,2)},{round(yc2,2)})")
        print(f"Vértices: A₁ = ({round(xa1,2)},{round(ya1,2)}); A₂ = ({round(xa2,2)},{round(ya2,2)}); B₁ = ({round(xb1,2)},{round(yb1,2)}); B₂ = ({round(xb2,2)},{round(yb2,2)})")
        print(f"Excentricidade: {round(e,2)}")
    elif y:
        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Focos: F₁ = ({round(xc1,2)},{round(yc1,2)}); F₂ = ({round(xc2,2)},{round(yc2,2)})")
        print(f"Vértices: A₁ = ({round(xa1,2)},{round(ya1,2)}); A₂ = ({round(xa2,2)},{round(ya2,2)}); B₁ = ({round(xb1,2)},{round(yb1,2)}); B₂ = ({round(xb2,2)},{round(yb2,2)})")
        print(f"Excentricidade: {round(e,2)}")

def hyperbola(a, c, f, h, k, sin, cos):
    print("\nA cônica é uma hipérbole.")
    
    a2 = f/a
    b2 = f/c

    if a2 < 0: a2 = -a2
    if b2 < 0: b2 = -b2

    newA = math.sqrt(a2)
    newB = math.sqrt(b2)
    newC = math.sqrt(a2+b2)
    e = newC/newA

    x = False
    y = False

    if f < 0:
        if a < 0 and c > 0:
            y = True
        else:
            x = True
    else:
        if a < 0 and c > 0:
            x = True
        else:
            y = True
    
    xc1,yc1,xa1,ya1,xb1,yb1,xc2,yc2,xa2,ya2,xb2,yb2 = changeCoord(x, y, h, k, newA, newB, newC, sin, cos)

    assymp1 = ""
    assymp2 = ""
    
    if x:
        origAssymp = newB/newA
    elif y:
        origAssymp = newA/newB

    """
    if sin or cos:
        l = (origAssymp*cos - sin)/(cos - origAssymp*sin)
        assymp1 += f"y = {round(l,2)}x"
        assymp2 += f"y = {-round(l,2)}x"
        
        if h or k:
            ind = k - l*h
            
            if ind > 0:
                assymp1 += f" + {round(ind,2)}"
                assymp2 += f" - {round(ind,2)}"
            if ind < 0:
                assymp1 += f" - {-round(ind,2)}"
                assymp2 += f" + {-round(ind,2)}"
    else:
        assymp1 += f"y = {round(origAssymp,2)}x"
        assymp2 += f"y = {-round(origAssymp,2)}x"
    """

    if x:
        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Focos: F₁ = ({round(xc1,2)},{round(yc1,2)}); F₂ = ({round(xc2,2)},{round(yc2,2)})")
        print(f"Vértices: A₁ = ({round(xa1,2)},{round(ya1,2)}); A₂ = ({round(xa2,2)},{round(ya2,2)})")
        print(f"Excentricidade: {round(e,2)}")
        print(f"Assíntotas: {assymp1}; {assymp2}")
    elif y:
        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Focos: F₁ = ({round(xc1,2)},{round(yc1,2)}); F₂ = ({round(xc2,2)},{round(yc2,2)})")
        print(f"Vértices: A₁ = ({round(xa1,2)},{round(ya1,2)}); A₂ = ({round(xa2,2)},{round(ya2,2)})")
        print(f"Excentricidade: {round(e,2)}")
        print(f"Assíntotas: Y = {assymp1}; {assymp2}")

def circumference(f):
    print("\nA cônica é uma circunferência.")

    if f < 0:
        r = math.sqrt(-f)
    else:
        r = math.sqrt(f)
    
    print(f"Raio = {round(r)}")

def parabola(a, c, d, e, f, sin, cos):
    print("\nA cônica é uma parábola.")
    
    if (d and e) or (a and e and f) or (c and d and f):
        print("\nTransladando a parábola:")
        a,c,d,e,u,v = simplifyParabola(a, c, d, e, f)

    if a == 0 and e == 0:
        if (c < 0 and d < 0) or (c > 0 and d > 0):
            p = -abs(d)/4
        else:
            p = abs(d)/4

        if sin or cos:
            h = cos*u - sin*v
            k = sin*u + cos*v

            px = h + cos*p
            py = k + sin*p

        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Foco: F = ({round(px,2)},{round(py,2)})")
        if d > 0:
            print(f"Diretriz: x = {round(abs(p),2)}")
        else:
            print(f"Diretriz: x = {-round(abs(p),2)}")
    elif c == 0 and d == 0:
        if (a < 0 and e < 0) or (a > 0 and e > 0):
            p = -abs(e)/4
        else:
            p = abs(e)/4

        if sin or cos:
            h = cos*u - sin*v
            k = sin*u + cos*v

            px = h - sin*p
            py = k + cos*p
        
        print("\nNo sistema de coordenadas original, com P = (x,y):")
        print(f"Origem: O = ({round(h,2)},{round(k,2)})")
        print(f"Foco: F = ({round(px,2)},{round(py,2)})")
        if e > 0:
            print(f"Diretriz: y = {round(abs(p),2)}")
        else:
            print(f"Diretriz: y = {-round(abs(p),2)}")

def classify(a, b, c, d, e, f, h, k, sin, cos):
    if b == 0 and d == 0 and e == 0:
        if f == 0:
            if (a > 0 and c > 0) or (a < 0 and c < 0):
                print("\nA equação representa um ponto.")
            elif (a > 0 and c < 0) or (a < 0 and c > 0):
                print("\nA equação representa duas retas concorrentes.")
            else:
                print("\nA equação representa duas retas idênticas.")
        elif f > 0:
            if a != -1 or c != -1:
                if a > 0 and c > 0:
                    print("\nVazio.")
                elif a < 0 and c < 0:
                    ellipse(a, c, f, h, k, sin, cos)
                elif (a < 0 and c > 0) or (a > 0 and c < 0):
                    hyperbola(a, c, f, h, k, sin, cos)
                elif (a < 0 and c == 0) or (a == 0 and c < 0):
                    print("\nA equação representa duas retas paralelas.")
                else:
                    print("\nVazio.")
            else:
                circumference(f)
        elif f < 0:
            if a != 1 or c != 1:
                if a < 0 and c < 0:
                    print("\nVazio.")
                elif a > 0 and c > 0:
                    ellipse(a, c, f, h, k, sin, cos)
                elif (a < 0 and c > 0) or (a > 0 and c < 0):
                    hyperbola(a, c, f, h, k, sin, cos)
                elif (a > 0 and c == 0) or (a == 0 and c > 0):
                    print("\nA equação representa duas retas paralelas.")
                else:
                    print("\nVazio.")
            else:
                circumference(f)
    elif b == 0 and ((a and d and c == 0 and e == 0) or (c and e and a == 0 and d == 0)):
        if f == 0 or f < 0:
            print("\nA equação representa duas retas paralelas.")
        else:
            print("\nVazio.")
    elif a == 0 and b == 0 and c == 0:
        if (d and e == 0) or (e and d == 0):
            print("\nA equação representa um ponto.")
        else:
            print("\nA equação representa uma reta.")
    else:
        parabola(a, c, d, e, f, sin, cos)

def main():
    # Coeficientes da equação
    a = float(input("Coeficiente de x²:   "))
    b = float(input("Coeficiente de xy:   "))
    c = float(input("Coeficiente de y²:   "))
    d = float(input("Coeficiente de x:    "))
    e = float(input("Coeficiente de y:    "))
    f = float(input("Termo independente:  "))

    # Equação original
    conic = ""
    conic = setEquation(conic, a, b, c, d, e, f, 1, 'x', 'y')
    print("\nEquação original: {}".format(conic))

    # Matriz da cônica
    matrix = np.zeros((3,3))
    setMatrix(matrix, a, b, c, d, e, f)

    # Operações de redução da equação
    h,k,sin,cos = 0,0,0,0

    if d != 0 or e != 0:
        a,b,c,d,e,f,h,k = translation(matrix, a, b, c, d, e, f)
    if b != 0:
        a,b,c,d,e,f,sin,cos = rotation(a, b, c, d, e, f)

    # Classificação
    classify(a, b, c, d, e, f, h, k, sin, cos)

if __name__ == "__main__":
    main()
