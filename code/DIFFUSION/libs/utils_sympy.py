import sympy as sp


def relation_fermeture(f,x,t,order : int,D ) : 
    """
    Calcule la relation de fermeture permettant d'appliquer la procédure de Cauchy-Kovaleskaya pour l'équation de la chaleur
    f       : fonction à variables spatio-temporelle.
    order   : ordre de la dérivée temporelle à évaluer
    """
    return D**order * sp.diff(f(x,t),x,2*order)

def developpement_taylor_spatial(f, h, order,x,t):
    """
    Calcule le développement de tailor autour du point générique x au poit x+h. 
    f       : fonction à estimer.
    h       : écart par rapport à x.
    order   : ordre du DL. 
    """
    return sum([
        h**k * sp.Rational(1,sp.factorial(k)) * sp.diff(f(x, t), x, k) 
        for k in range(order + 1)
    ]) 

def developpement_taylor_temps(f, h, order,x,t):
    """
    Calcule le développement de tailor autour du point générique t au poit t+h. 
    f       : fonction à estimer.
    h       : écart par rapport à t.
    order   : ordre du DL. 
    """
    return sum([
        h**k * sp.Rational(1,sp.factorial(k)) * sp.diff(f(x, t), t, k) 
        for k in range(order + 1)
    ]) 

def ordre_total_derivée(expr):
    """
    prend un fonction et renvoie son ordre de dérivation
    """
    return sum(o for _, o in expr.variable_count) if isinstance(expr, sp.Derivative) else 0

def apply_cauchy_kovaleskaya(expr,u,x,t,D,max_order):
    """
    applique Cauchy-Kovaleskaya à l'argument expr utilisant la relation de fermuture : dtu = D * dxxu
    """
    expr.expand() # Développe l'expression
    coeffs_dtn = sp.collect( # Crée un dictionnaire avec en clée les dérivées temporelle aux différents ordres qui apparaissent dans expr et en valeur les coefficients associés 
        expr ,
        [sp.diff(u(x, t), t,localOrder) for localOrder in range(2,max_order+1) ],
        evaluate=False # Renvoie le dictionnaire
        )
    for (key,value) in coeffs_dtn.items(): # On itère sur les dérivées temporelles (key) en ayant accès a leurs coeffcients (value)
        ordre  = ordre_total_derivée(key)  # On récupère l'ordre de la dérivée temporelle sur laquelle on travaile
        if ordre == 0 : continue
        terme_temporel = key * value # On reconstruit le terme associé à la dérivée tempoelle
        terme_spatial = relation_fermeture(u,x,t,ordre , D) * value # On construit le terme associé faisant intervenir une dérivée spatiale à la place (c'est ici que l'on a besoin de l'ordre de la dérivée temporelle)
        expr = expr - terme_temporel + terme_spatial # On remplace le terme faisant intervenir la dérivée temporelle et on le remplace par celui mettant en jeu la dérivée spatiale
    return expr