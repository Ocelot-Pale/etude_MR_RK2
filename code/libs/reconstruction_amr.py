import sympy as sp
from libs.utils_sympy import * 
from libs.display import afficheSimpy
def initP(stencil=1) :
    if stencil==1 : 
        P = sp.zeros(4,4)
        P[0,0] = sp.Rational(1,8)
        P[1,0] = sp.Rational(-1,8)
        P[2,1] = sp.Rational(1,8)
        P[3,1] = sp.Rational(-1,8)
        P[0,2] = sp.Rational(-1,8)
        P[1,2] = sp.Rational(1,8)
        P[2,3] = sp.Rational(-1,8)
        P[3,3] = sp.Rational(1,8)
        P[0,1] = 1
        P[1,1] = 1
        P[2,2] = 1
        P[3,2] = 1
        return P 
    else : raise NotImplementedError(f"not implemented s={stencil}")

def computeFlux(u,x,t,dl,dx,cfl,ORDER_SERIES,stencil=1):
    dx_currentLevel = (2**(dl)) * dx # Calcul du deltaX courrant
    cfl_effective = cfl*(2**(-dl))
    Pdl = initP(stencil)**dl # Calcul de la matrice de passage au travers de deltaL couches
    coeffs = sp.Matrix([[
        sp.Rational(-1,2)*cfl_effective,
        +sp.Rational(3,2)*cfl_effective-1,
        -sp.Rational(3,2)*cfl_effective+1,
        sp.Rational(1,2)*cfl_effective,
        ]])
    
    coeffs_1 = sp.Matrix([[0,-1,1,0,]])
    coeffs_2 = sp.Matrix([[sp.Rational(-1,2),sp.Rational(-3,2),sp.Rational(3,2),sp.Rational(1,2),]])
    
    Us_minusFlux_currentLevel = sp.Matrix([developpement_taylor_spatial(u,dx_currentLevel*delta,ORDER_SERIES,x,t) for delta in range(-2,2)])
    Us_plusFlux_currentLevel = sp.Matrix([developpement_taylor_spatial(u,dx_currentLevel*delta,ORDER_SERIES,x,t)  for delta in range(-1,3)])
    minusFlux = (coeffs * Pdl * Us_minusFlux_currentLevel)[0,0]
    plusFlux  = (coeffs * Pdl *  Us_plusFlux_currentLevel)[0,0]
    flux_variation_RK2_v1 = cfl_effective * (plusFlux-minusFlux)
    minusFlux = ((cfl_effective*coeffs_1+cfl_effective**2*coeffs_2)*(Pdl * Us_minusFlux_currentLevel))[0,0]
    plusFlux  = ((cfl_effective*coeffs_1+cfl_effective**2*coeffs_2)*(Pdl *  Us_plusFlux_currentLevel))[0,0]
    flux_variation_RK2 = (plusFlux-minusFlux)


    coeffs_RK1 = sp.Matrix([[0,-1,1,0,]])
    minusFlux_RK1 = (coeffs_RK1 * Pdl * Us_minusFlux_currentLevel)[0,0]
    plusFlux_RK1 = (coeffs_RK1 * Pdl * Us_plusFlux_currentLevel)[0,0]
    flux_variation_RK1 = cfl_effective * (plusFlux_RK1-minusFlux_RK1)

    coeffs_Christian = sp.Matrix([[
        sp.Rational(-1,2)*(cfl-sp.Rational(1,6)) ,
        sp.Rational(3,2)*(cfl-sp.Rational(1,6))-1,
        1-sp.Rational(3,2)*(cfl-sp.Rational(1,6)),
        sp.Rational(1,2)*(cfl-sp.Rational(1,6)) ,
                        ]])
    minusFlux_Christian = (coeffs_Christian * Pdl * Us_minusFlux_currentLevel)[0,0]
    plusFlux_Christian = (coeffs_Christian * Pdl * Us_plusFlux_currentLevel)[0,0]
    flux_variation_Christian = cfl_effective * (plusFlux_Christian-minusFlux_Christian)

    return  flux_variation_RK2,flux_variation_RK1,flux_variation_Christian


def computeEqModif(u,x,t,dl,dx,dt,D,cfl,ORDER_SERIES,CAUCHY_KOVALESKAYA=True,stencil=1,DEBUG=0) : 
    flux_total_RK2,flux_total_RK1,flux_total_Christian = computeFlux(u,x,t,dl,dx,cfl,ORDER_SERIES,stencil=1) 
    rhs = (1/dt*flux_total_RK2).expand()
    rhs_RK1 = (1/dt*flux_total_RK1).expand()
    rhs_Christian = (1/dt*flux_total_Christian).expand()
    lhs = (1/dt * (developpement_taylor_temps(u,dt,ORDER_SERIES,x,t) - u(x,t))).expand()
    rhs = rhs - (lhs - sp.diff(u(x,t),t,1))
    rhs_RK1 = rhs_RK1-(lhs - sp.diff(u(x,t),t,1))
    rhs_Christian = rhs_Christian-(lhs - sp.diff(u(x,t),t,1))
    lhs = sp.diff(u(x,t),t,1)
    if CAUCHY_KOVALESKAYA : 
        rhs=apply_cauchy_kovaleskaya(rhs,u,x,t,D,ORDER_SERIES)  
        rhs_RK1=apply_cauchy_kovaleskaya(rhs_RK1,u,x,t,D,ORDER_SERIES)  
        rhs_Christian=apply_cauchy_kovaleskaya(rhs_Christian,u,x,t,D,ORDER_SERIES)  
    tex =sp.latex(lhs) + "=" + sp.latex(rhs)
    tex_RK1 =sp.latex(lhs) + "=" + sp.latex(rhs_RK1)
    tex_Christian = sp.latex(lhs) + "=" + sp.latex(rhs_Christian)
    return {"RK2":(lhs,rhs,tex),"RK1":(lhs,rhs_RK1,tex_RK1),"Christian" : (lhs,rhs_Christian,tex_Christian)}