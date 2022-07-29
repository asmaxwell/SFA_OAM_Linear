# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Pulse version of the SFA for circular polarization
"""
Created on Tuesday Nov 24 10:29:00 2020

@author: asmaxwell

SFA pulse class
"""
#import scipy.optimize as op
import scipy.integrate as it
import scipy.special as sp
cimport scipy.special.cython_special as csp
import mpmath as mp
import functools
import multiprocessing

import numpy as np
cimport numpy as np
from numpy cimport ndarray

cimport cython

from libcpp cimport bool

from libc.math cimport sin as sin_re
#from libc.math cimport sinc as sinc_re
from libc.math cimport cos as cos_re
#from libc.math cimport cot as cot_re
from libc.math cimport acos as acos_re
from libc.math cimport exp as exp_re
from libc.math cimport sqrt as sqrt_re
from libc.math cimport abs as abs_re
from libc.math cimport pow as pow_re
from libc.math cimport atan2 as atan2_re
#from libc.math cimport  cyl_bessel_j as bessel_j_re


#fused types
ctypedef fused dbl_or_cmplx:
    double
    double complex
    
ctypedef fused c_dbl_int:
    int
    double
    double complex
    
cdef extern from "<complex.h>" namespace "std" nogil:
    double complex exp(double complex z)
    double complex sin(double complex z)
    double complex cos(double complex z)
    #double complex cot(double complex z)
    double complex sqrt(double complex z)
    double complex acos(double complex z)
    double complex log(double complex z)
    double complex atan2(double complex z)
    #double complex pow(double complex z, c_dbl_int z)
    double real(double complex z)
    double imag(double complex z)


cdef double complex I1 = 1j
cdef double Pi = np.pi
cdef double rtPi = np.sqrt(Pi)
cdef double rt2 = np.sqrt(2.)
cdef int cacheSize = 2**20

#Code to import old CQSFA prefactor code
# cdef extern from "HpgForms.h":
#     double complex calculateHMatrixElement(int target, double Ip, double complex pz, double px,  double complex ts, double complex Eft, double complex Aft, double theta)

### shorcut functions to efficntly switch trig between real and complex varients    
cdef sin_c(dbl_or_cmplx t):
        if(dbl_or_cmplx is double):
            return sin_re(t)
        else:
            return sin(t)
cdef cos_c(dbl_or_cmplx t):
        if(dbl_or_cmplx is double):
            return cos_re(t)
        else:
            return cos(t)
        
cdef cot_c(dbl_or_cmplx t):
        return cos_c(t)/sin_c(t)




cdef class SFALinearPulse:
    '''
        Class to compute the transition amplitude M(p) and its dervative M_g(p) using the SFA and saddle point approximation
    '''
    #memeber variables like in C++!
    cdef readonly double Ip, Up, rtUp, omega, CEP, AlignmentAngle
    cdef readonly int N, Target
    #cdef object __weakref__ # enable weak referencing support
    
    def __init__(self, Ip_ = 0.5, Up_ = 0.44, omega_ = 0.057, N_ = 6, CEP_ = 0., Target_=3.):
        '''
            Initialise field and target parameters defaults correspond to 800nm wl and 2 10^14 W/cm^2 intensity
            for the target 0=He, HeTheta=1, Ne=2, Ar=3, ArEx_4S=4, Xe=5, N2=6, N2Theta=7, O2Theta=8, H = 9 (and default case e.g. any other number)
            The ionization prefactor must be set independently
        '''
        #Set pulse and targetvparameters
        self.Ip = Ip_
        self.Up = Up_
        self.rtUp = np.sqrt(Up_) #must change this if Up is changed! Fixed by making Up readonly
        self.omega = omega_
        
        
        self.N = N_
        self.CEP = CEP_
        
        self.AlignmentAngle = 0.
        self.Target = Target_
    
        
    #@functools.lru_cache(maxsize=cacheSize)
    cdef dbl_or_cmplx F(s, dbl_or_cmplx t):
        '''
        envelope for Sin^2 laser pulse
        '''
        #need to be fast can get evalutated millions of times
        if(real(t)<0 or real(t)>2*s.N*Pi/s.omega):
            return 0
        return 2*s.rtUp *  sin_c(s.omega * t / (2*s.N))**2
    
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Af(s, dbl_or_cmplx t):
        '''
        Vector potential for Gaussian laser pulse
        '''
        if(real(t)<0 or real(t)>2*s.N*Pi/s.omega):
            return 0
        cdef dbl_or_cmplx envelope, carrier
        envelope = s.F(t)
        if(dbl_or_cmplx is double):
            carrier = cos_re(s.omega*t + s.CEP)
        else:
            carrier = cos(s.omega*t + s.CEP)
        return  envelope * carrier
    
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Ef(s, dbl_or_cmplx t):
        '''
        Electric Field for Sin^2 laser pulse
        '''
        if(real(t)<0 or real(t)>2*s.N*Pi/s.omega):
            return 0
        cdef dbl_or_cmplx cos1, sin1, cos2, sin2, env
        env = s.omega*s.F(t)
        cos1 = cos_c(s.omega*t + s.CEP)
        sin1 = sin_c(s.omega*t + s.CEP)
        
        cos2 = cos_c(s.omega*t/(2*s.N))
        sin2 = sin_c(s.omega*t/(2*s.N))
        return env*sin1-2*s.rtUp*s.omega *cos1*cos2*sin2/s.N
    
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx AfI(s, dbl_or_cmplx t):
        '''
            Integral of vector potential
        '''
        cdef double factor, a1, a2, a3, AI0=0
        cdef dbl_or_cmplx s1, s2, s3
        #Case if N==1 as general expression will diverge
        if(s.N==1):
            A10 = (3*s.rtUp*sin_re(s.CEP))/(4.*s.omega)
            return -(s.rtUp*(2*t*s.omega*cos_re(s.CEP) - 4*sin_c(t*s.omega + s.CEP) + sin_c(2*t*s.omega + s.CEP)))/(4.*rt2*s.omega) - AI0
        else:
            AI0 = (s.rtUp*sin_re(s.CEP))/(s.omega - s.N**2*s.omega) # AfI(0) to ensure limits are correct.
            a1 = s.N*s.N - 1
            a2 = (s.N + 1.)/s.N
            a3 = (s.N - 1.)/s.N
            factor = s.rtUp/(2*s.omega*a1)

            s1 = sin_c(s.omega*t + s.CEP)
            s2 = sin_c(s.CEP + a2 *s.omega*t)
            s3 = sin_c(s.CEP + a3 *s.omega*t)

            return factor * (2*a1*s1 - s.N*s.N *( a3*s2 + a2*s3) ) - AI0

    

    #@functools.lru_cache(maxsize=cacheSize)
    cpdef dbl_or_cmplx Af2I(s, dbl_or_cmplx t):
        '''
            Integral of vector potential squared
        '''     
        if(s.N==1):
            AI20 = (-25*s.Up*sin_re(2*s.CEP))/(96.*s.omega)
            return (3*s.Up*(4*(6*t*s.omega + t*s.omega*cos_re(2*s.CEP) -8*sin_c(t*s.omega) + sin_c(2*t*s.omega) + 3*sin_c(2*(s.CEP + t*s.omega)) - 4*sin_c(2*s.CEP + t*s.omega)) + sin_c(2*(s.CEP + 2*t*s.omega))) - 16*s.Up0*sin_c(2*s.CEP + 3*t*s.omega))/(96.*s.omega)
        
        AI20 = (3*s.Up*cos_re(s.CEP)*sin_re(s.CEP))/(4*s.omega - 20*s.N**2*s.omega + 16*s.N**4*s.omega)
        cdef dbl_or_cmplx s1, s2, s3, s4, s5, s6, s7, c1
        cdef double a1 = s.omega*(s.N+1)/(s.N), a2 = s.omega*(s.N-1)/(s.N), a3 = s.omega*(2*s.N+1)/(s.N), a4 = s.omega*(2*s.N-1)/(s.N)
           
        s1 = sin_c(2*s.omega*t)
        s2 = sin_c(s.omega*t/s.N)
        s3 = sin_c(2*s.omega*t/s.N)
        s4 = sin_c(2*s.CEP + a3*t)
        s5 = sin_c(2*s.CEP + 2*a2*t)
        s6 = sin_c(2*s.CEP + 2*a1*t)
        s7 = sin_c(2*s.CEP + a4*t)
        c1 = cos_c(2*s.omega*t)
            
        return (s.Up/16) * (12*t + 6*c1*sin_re(2*s.CEP)/s.omega + 6*s1*cos_re(2*s.CEP)/s.omega
                           -(16*s.N/s.omega)*s2 + (2*s.N/s.omega)*s3 - 8*s4/a3 + s5/a2 + s6/a1 - 8*s7/a4)
        
        
    #@functools.lru_cache(maxsize=cacheSize)    
    cpdef dbl_or_cmplx S(s, double p, double theta, double phi, dbl_or_cmplx t):
        '''
            Action as given by the SFA for a Pulse
        '''
        cdef dbl_or_cmplx tTerms = (s.Ip + 0.5*p*p)*t
        cdef dbl_or_cmplx linAI = p*cos_re(theta)*s.AfI(t)
        cdef dbl_or_cmplx quadAI = 0.5*s.Af2I(t)
        return tTerms + linAI + quadAI
    
    
    cpdef dbl_or_cmplx DS(s, double p, double theta, double phi, dbl_or_cmplx t):
            cdef pz = p*cos_re(theta)
            return pz*s.Af(t) + 0.5*s.Af(t)*s.Af(t) + 0.5*p*p + s.Ip
    
    #@functools.lru_cache(maxsize=cacheSize)
    #REDO for linear field and with CEP!!!!!
    cpdef DSZ(s, double p, double theta, double phi):
        '''
            Derivative of the action tranformed by t->i N/omega Log[z] for esay solving
            This creates an 2(N+1) polynomial which can be efficeintly solved and the solutions easily 
            transformed back. This function passes the roots of the polynomial as numpy array so they can be solved using np.roots.          
            It is clear there will be N+1 solutions and their complex 
            conjugates.
        '''
        #Note for the linear case input varible phi is redundant, keep for genrality though
        
        cdef double complex exp_CEP = exp(I1*s.CEP)
        
        #costruct polynomial in z of order 4*(N+1)
        poly_coeffs = np.zeros(4*s.N+4+1) + 0.*I1
        
        #0 order terms
        cdef c0 = s.Up*(exp_CEP*exp_CEP) 
        poly_coeffs[0:5] = [c0/16, -c0/4, 3*c0/8, -c0/4, c0/16]
        
        #N order terms (+= accounts for cases where coefficients combine)
        cdef c1 = p*s.rtUp*cos_re(theta)*exp_CEP
        poly_coeffs[s.N+1:s.N+4] += [-c1/2, c1, -c1/2]
        
        #2N order terms
        poly_coeffs[2*s.N:2*s.N+5] += [s.Up/8, -s.Up/2, 2*s.Ip + p*p + 3*s.Up/4, -s.Up/2, s.Up/8]
        
        #3N order terms
        cdef c3 = p*s.rtUp*cos_re(theta)/exp_CEP
        poly_coeffs[3*s.N+1:3*s.N+4] += [-c3/2, c3, -c3/2]
        
        #4N order terms
        cdef c4 = s.Up/(exp_CEP*exp_CEP) 
        poly_coeffs[4*s.N:] += [c4/16, -c4/4, 3*c4/8, -c4/4, c4/16]
        
        
        return poly_coeffs
    
    cpdef double complex DSZ_val(s, double p, double theta, double phi, dbl_or_cmplx z):
        poly_coeffs = s.DSZ(p, theta, phi)
        cdef double complex sum_val = 0
        for n in range(0, len(poly_coeffs)):
            sum_val += poly_coeffs[n] * z**n
        return sum_val
    
    cdef double complex addIfRealNeg(s, double complex ts):
        if(real(ts)<0):
            return ts + 2*Pi*s.N/s.omega
        else:
            return ts
        
    #@functools.lru_cache(maxsize=cacheSize)
    #Ensure correction transformation is done
    cpdef TimesGen(s, double p, double theta, double phi):
        '''
            Solution for times found by transforming the derivative of the action into
            a 4(N+1) polynomial and solving using np.roots. This should be a very effiecint way to get 
            all roots
        '''
        poly_coeffs = s.DSZ(p, theta, phi)
        z_roots = np.polynomial.polynomial.polyroots(poly_coeffs)
        #now we must transform back using t=I N Log[z]/omega
        ts_roots = I1*s.N*np.log(z_roots)/s.omega
        ts_roots = [ts for ts in ts_roots if imag(ts)>0 ] #remove (divergent) solutions with negative imag
        #make sure all t values are in the domain [0, 2*pi*N/omega]
        ts_roots = [s.addIfRealNeg(ts) for ts in ts_roots]
        #sort real parts to easily select specific solutions        
        return sorted(ts_roots,  key=np.real)
        
    #1 varible determinant for saddle point approximation
    cpdef double complex DDS(s, double p, double theta, double phi, double complex t):
        '''Second order derivative of action wrt t'''
        return -(p*cos_re(theta)+s.Af(t))*s.Ef(t)
        

    
# prefactor keep for future implementation    
#     #@functools.lru_cache(maxsize=cacheSize)
    cdef double complex d0(s, double p, double theta, double phi, double complex ts):
        '''
            Bound state prefactor will use GAMESS style gaussain orbitals'
        '''
        #Dummy prefactor changes this!
        return 1.
    

    
    @cython.boundscheck(False) # turn off bounds-checking for entire function  
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef double complex M(s, double p, double theta, double phi, double tf = np.inf):#double pz, double px, double t, int N, int eLim):
        '''
            Final transition amplitude
            Constructed as sum 
        '''
        #eLim is 1 or 2, the number of orbits accounted for
        cdef double complex MSum = 0.
        times = s.TimesGen(p, theta, phi)
        
        for ts in times:
            if(real(ts)<tf):
                det = sqrt(2*Pi*I1/s.DDS(p, theta, phi, ts))
                expS = exp(I1*s.S(p, theta, phi, ts))
                d0 = s.d0(p, theta, phi, ts)
                MSum += d0*det*expS
        return MSum
    #transition amplitude in cartesian co-ordinates
    cpdef double complex Mxy(s, px, py, pz, tf = np.inf):
        cdef double p = sqrt_re(px*px + py*py +pz*pz)
        cdef double theta = acos_re(pz/p)
        cdef double phi = atan2_re(py, px)
        return s.M(p, theta, phi, tf)
    #list comprehension over cartesian transition amplitude
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Mxy_List(s, pxList, pyList, double pz, tf = np.inf):
        return np.array([s.Mxy(px, py, pz, tf) for px, py in zip(pxList, pyList)])
    
    def Mxz_List(s, pxList, double py, pzList, tf = np.inf):
        return np.array([s.Mxy(px, py, pz, tf) for px, pz in zip(pxList, pzList)])
   
    ####   ---   OAM Functions   ---   #### 
    cpdef Ml(s, double p, double theta, int Nphi = 250):
        '''
            This is the fourier series coeiffint of M to get the OAM distribusion.
            It is computed taking advantage of the FFT
        '''
        phiList = np.linspace(-Pi, Pi, Nphi)
        MphiList = [s.M(p, theta, phi) for phi in phiList]
        return np.fft.fft(MphiList)/Nphi
    
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Ml_List(s, pList, theta, Nphi = 250):
        return np.array([[abs(M)**2 for M in s.Ml(p, theta, Nphi)] for p in pList]).T
    
    cpdef Mlxz(s, px, pz, int Nphi = 250):
        '''
        convert Ml to cartesian coordinates, note px is the perpendicular coordinate st. px^2 = px^2+py^2
        '''
        cdef double p = sqrt_re(px*px +pz*pz)
        cdef double theta = acos_re(pz/p)
        return s.Ml(p, theta, Nphi)
    
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    def Mlxz_List(s, pxList, pzList, Nphi = 250):
        return np.array([[abs(M)**2 for M in s.Mlxz(px, pz, Nphi)] for px, pz in zip(pxList, pzList)])
 
 #####   ---   code for spectra
    cpdef double Spectra(s, double E, double phi = 0., double t = np.inf, double err = 1.0e-4, int limit = 500):
        '''Function for the spectra'''    
        Norm_val, Norm_error = it.quad(s.Spec_Norm, 0, Pi, args = (phi, E, t), epsabs=err, epsrel=err, limit=limit  )    
        return Norm_val
     
    #Spectra Norm integrand
    cpdef double Spec_Norm(s, double theta, double phi, double E, double t = np.inf):
        '''Function to compute the integrand of the theta integral for the 'spectra''''
        #phrased in spherical coordinates
        #cdef double complex M, Mg
        cdef double px, pr
        pr = sqrt_re(2*E)
        px = pr*sin_re(theta)
        
        return abs(s.M(pr, theta, phi, t))**2 *px * 2*Pi

######   ---   Funcitons for the analytic monochromatic approximation to pulse   ---   #####
# ### This will no longer work as this was written for a circular field


#     cpdef Ml_Mono(s, double p, double theta, l):
#         '''
#                 function to return the fully analytic monchromatic approximation to 
#         '''
#         cdef double E1 = (2*s.Ip+2*s.Up+p*p)
#         #compute prefactor
#         cdef double part1 = E1*E1/(4*p*p*sin_re(theta))
#         cdef double complex DDS = p*s.omega*sqrt(2*s.Up-part1)*sin_re(theta)
#         cdef double complex pref = sqrt(2*Pi*I/DDS)

#         #compute action
#         cdef double complex S0_1 = E1*acos(E1/(2*s.rt2Up*p*sin_re(theta))-1e-16 * I) #selecting +ive imaginary part branch of acos
#         cdef double S0_2 = sqrt_re(E1*E1+4*p*p*s.Up*(cos_re(2*theta)-1))
#         cdef double complex S0 = (S0_1 - I*S0_2)/(2*s.omega)

#         #compute sinc
#         sinc_val = np.sinc((E1-2*s.omega*l)/(2*s.omega))

#         #compute ATI rings
#         cdef double complex OM = (exp(I*s.N*Pi*E1/s.omega)-1)/(exp(I*Pi*E1/s.omega)-1)
#         #print('OM = ',OM,', pref = ',pref,', S0 = ',S0_1/(2*s.omega),', sinc = ',sinc_val)
#         return OM*pref*exp(I*S0)*sinc_val
         
    
# ######   ---   Here functions are defined for the analytical fourier series computation   ---   #####
# ### This will no longer work as this was written for a circular field

#     cdef double OAM_S(s, double p, double theta, int l, double t):
#         '''
#             Action for the analytical OAM action
#         '''
#         cdef double constTerms = (s.Ip + 0.5*p*p)*t
#         cdef double quadTerms = 0 if (t<0 or t>2*s.N*Pi/s.omega) else -0.5*s.Af2I(t)
#         cdef double trigTerm = 0 if (t<0 or t>2*s.N*Pi/s.omega) else atan2_re(s.AfxI(t), s.AfyI(t))*l
        
#         return constTerms + quadTerms + trigTerm 
    
#     cpdef double complex OAM_Integrand(s, double t, double p, double theta, int l):
#         cdef double S1 = s.OAM_S(p, theta, l, t)
#         cdef double complex exp_S = exp(I*S1)
#         cdef AI_abs = 0 if (t<0 or t>2*s.N*Pi/s.omega) else sqrt_re(s.AfxI(t)**2 + s.AfyI(t)**2)
#         cdef double BesJ = sp.jv(l, p*sin_re(theta)*AI_abs) 
        
#         return BesJ*exp_S
    
#     cpdef double OAM_IntegrandRe(s, double t, double p, double theta, int l):
#         cdef double S1 = s.OAM_S(p, theta, l, t)
#         cdef double cos_S = cos_re(S1)
#         cdef AI_abs = 0 if (t<0 or t>2*s.N*Pi/s.omega) else sqrt_re(s.AfxI(t)**2 + s.AfyI(t)**2)
#         cdef double BesJ = sp.jv(l, p*sin_re(theta)*AI_abs) 
        
#         return BesJ*cos_S
    
#     cpdef double OAM_IntegrandIm(s, double t, double p, double theta, int l):
#         cdef double S1 = s.OAM_S(p, theta, l, t)
#         cdef double sin_S = sin_re(S1)
#         cdef AI_abs = 0 if (t<0 or t>2*s.N*Pi/s.omega) else sqrt_re(s.AfxI(t)**2 + s.AfyI(t)**2)
#         cdef double BesJ = sp.jv(l, p*sin_re(theta)*AI_abs) 
        
#         return BesJ*sin_S
    
#     cpdef double complex OAM_Ml(s, double p, double theta, int l, double err = 1.0e-4, int limit = 2000):
#         cdef double valRe, valIm, errorRe, errorIm
#         valRe, errorRe = it.quad(s.OAM_IntegrandRe, 0, 2*s.N*Pi/s.omega, args = (p, theta, l), epsabs=err, epsrel=err, limit=limit  )
#         valIm, errorIm = it.quad(s.OAM_IntegrandIm, 0, 2*s.N*Pi/s.omega, args = (p, theta, l), epsabs=err, epsrel=err, limit=limit  )
#         #tList = np.linspace(0, 2*s.N*Pi/s.omega, Nt)
#         #MtList = [s.OAM_Integrand(t, p, theta, l) for t in tList]
        
#         return valRe + I*valIm#[valRe + I*valIm, errorRe + I*errorIm]#np.fft.fft(MphiList)/Nphi
    
#     @cython.boundscheck(False) # turn off bounds-checking for entire function
#     @cython.wraparound(False)  # turn off negative index wrapping for entire function
#     def OAM_M_List(s, pList, theta, l, err = 1.0e-4, limit = 2000):
#         return np.array([np.abs(s.OAM_Ml(p, theta, l, err, limit))**2 for p in pList])

        
        