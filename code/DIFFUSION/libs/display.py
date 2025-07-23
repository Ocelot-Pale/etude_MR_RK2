from IPython.display import display, Math
import sympy as sp 


def afficheSimpy(expr) : 
    display(Math(sp.latex(expr)))
    return 

def rendreJoli(expr, u, x, t, dx, dt, ORDER_SERIES):
    """
    Retourne un latex isolant le puissances des deltaX et deltaT à l'affichage.
    """
    orders = range(0, ORDER_SERIES + 1) # Liste des ordres de dérivés temporelles que l'expression est suceptible de porter 
    dxs = [dx**n for n in orders] # Liste des puissances de deltaX 
    dts = [dt**n for n in orders] # Liste des puissances de deltaT 
    expr = expr.expand() # développe l'expression
    collected_dx = sp.collect(expr, dxs, evaluate=False) # cére un dictionnaire avec en clé les puissances de deltaX et en valeur les coeffs associés
    
    latex_lines = [] # liste des termes
    for dx_term, dx_expr in collected_dx.items(): # On itère sur les puissances de detaX (dx_term) et on accès au coefficient associés (dx_expr)
        collected_dt = sp.collect(dx_expr, dts, evaluate=False) # Crée un dictionnaire avec en clée les puissances de deltaT et en valeur les coeffs 
        for dt_term, coef in collected_dt.items():
            factor_parts = []
            if dt_term != 1:
                factor_parts.append(sp.latex(dt_term))
            if dx_term != 1:
                factor_parts.append(sp.latex(dx_term))
            factor_str = " ".join(factor_parts)
            
            if factor_str:
                latex_lines.append(f"{factor_str} \\left({sp.latex(coef)}\\right)")
            else:
                latex_lines.append(sp.latex(coef))
    
    return " + ".join(latex_lines)

