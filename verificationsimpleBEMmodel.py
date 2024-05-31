import numpy as np
import simpleBEMmodel as BEM


#verifying CTfunction
#test 1: induction 0 gives CT 0
def test1CTfunction():
    a = 0
    CT = BEM.CTfunction(a)
    return CT == 0

print(test1CTfunction())
#test 2 : induction with glauert increases CT

def test2CTfunction():
    a = np.array([0.5,0.6])
    CT_false = BEM.CTfunction(a)
    CT_true = BEM.CTfunction(a, glauert=True)
    return CT_true.all()
print(test2CTfunction())

#testing ainduction function
#Check for edge cases like extremely high or low values of CT.

def test1ainduction():
    CT_1 = 0
    CT_2 = 1
    a = BEM.ainduction(np.array([CT_1, CT_2]))

