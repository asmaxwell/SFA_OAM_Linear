# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Pulse version of the SFA for linear polarization
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
from sympy.physics.wigner import wigner_3j


import numpy as np
cimport numpy as np
from numpy cimport ndarray
from scipy.special import sph_harm, hyp2f1, lpmv, clpmn, lpmn
import math
import mpmath as mp

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
    double complex atan(double complex z)
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
    cdef readonly int N, Gauge
    cdef readonly list n, l, m, const1
    #cdef object __weakref__ # enable weak referencing support
    
    def __init__(self, Ip_ = 0.5, Up_ = 0.44, omega_ = 0.057, N_ = 6, n_=[1], l_=[0], m_=[0], const1_=[1], CEP_ = 0., Gauge_=0):
        '''
            Initialise field and target parameters defaults correspond to 800nm wl and 2 10^14 W/cm^2 intensity
            Gauge_=0 (velocity gauge), Gauge_=1 (length gauge)

        '''
        #Set pulse and targetvparameters
        self.Ip = Ip_
        self.Up = Up_
        self.rtUp = np.sqrt(Up_) #must change this if Up is changed! Fixed by making Up readonly
        self.omega = omega_
        
        #Target state parameters
        self.n = n_
        self.l = l_
        self.m = m_
        self.const1 = const1_
        
        
        self.N = N_
        self.CEP = CEP_
        self.Gauge = Gauge_
        
        self.AlignmentAngle = 0.
    
        
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
        
    
    cpdef double complex d0_v(s, nlist, llist, mlist, constlist, double p, double theta, double phi):
        '''
            Bound state prefactor in the velocity gauge <p|V|0> for a generic hydrogen wavefunction
        '''
        cdef double complex pref=0.
        for n,l,m,const in zip(nlist,llist,mlist,constlist):
            sum1=0.
            for k in range(n-l):
                term1=(-1)**k*math.factorial(n+l)*(2**k)*(2*s.Ip)**((-0.5-l)/2)*p**l
                term2=math.factorial(n-l-k-1)*math.factorial(2*l+k+1)*math.factorial(k)
                term3=math.gamma(2+k+2*l)/math.gamma(3/2+l)
                term4=hyp2f1(1+l+k/2,(3+k+2*l)/2, l+3/2,-p**2/(2*s.Ip))
                sum1=sum1+term1*term3*term4/term2
            term5=(-1j)**l*math.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
            Yml=sph_harm(m,l,phi,theta)
            pref=pref+term5*Yml*sum1*const
        pref=pref/math.sqrt(len(nlist))
        return pref
    
    
    cpdef double complex d0_l(s, nlist, llist, mlist, constlist, double p, double theta, double phi, dbl_or_cmplx ts):
        '''
            Bound state prefactor in the length gauge <p+A(t)|V|0> for a generic hydrogen wavefunction
        '''
        cdef double px=p*sin_re(theta)*cos_re(phi), py=p*sin_re(theta)*sin_re(phi), pz=p*cos_re(theta)
        cdef double sn = 1 if s.Af(ts).imag > 0 else -1
        cdef double complex pz_2 = 1j*sn*sqrt_re(2*s.Ip + px**2+py**2)#pz+s.Af(ts) #tilde{pz}
        cdef double complex pmod=  1j*sqrt_re(2*s.Ip)#sqrt(px**2+py**2+pz_2**2) #=modulus{tilde{p}}
        cdef double complex pE=pz_2*s.Ef(ts) #=tilde{p}*E(ts)
        cdef double complex sum1
        cdef dbl_or_cmplx term1
        cdef double complex Yml
        
        cdef double complex pref=0.
        for n,l,m,const in zip(nlist,llist,mlist,constlist):
            sum1=0. 
            for k in range(n-l):
                term2=(-1)**k*math.factorial(n+l)*(2**k)*(2*s.Ip)**(3/4+l/2+k/2)*pmod**l
                term3=(1j*pE)**((2+k+2*l)/4)
                term4=math.factorial(n-l-k-1)*math.factorial(k)*math.gamma(3/2+l)*math.gamma(1+l+k/2)            
                if k==0:
                    term5=math.gamma((2+2*l)/4)              
                elif k==1:
                    F1=complex(mp.hyper([-1/4.,1/4.],[1/2.,(5+2*l)/4.],s.Ip**2/(1j*pE)))
                    F2=complex(mp.hyper([1/4.,3/4.],[3/2.,(7+2*l)/4.],s.Ip**2/(1j*pE)))
                    term5=math.gamma((3+2*l)/4)*F1-sqrt(4*s.Ip**2/(1j*pE))*math.gamma((5+2*l)/4)*F2/(3+2*l)
                else:
                    return np.nan
                sum1=sum1+term2*term5/(term3*term4)       
            term1=(-1j)**l*math.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
            
            if m>=0:
                Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l)))*exp(1j*m*phi) 
                Yml=Yml*lpmn(m,l,pz_2.imag/pmod.imag)[0][m][l]  
            else: 
                #Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l)))*exp(1j*m*(phi)) 
                #Yml=Yml*(-1)**abs(m)*(math.factorial(l-abs(m))/math.factorial(l+abs(m)))*clpmn(abs(m),l,pz_2/pmod)[0][abs(m)][l] 
                Yml=sqrt((2*l+1)*math.factorial(l-abs(m))/(4*np.pi*math.factorial(abs(m)+l)))*exp(1j*abs(m)*phi) 
                Yml=Yml*lpmn(abs(m),l,pz_2.imag/pmod.imag)[0][abs(m)][l]  
                Yml=((-1)**abs(m))*np.conj(Yml)
            
            pref=pref+term1*Yml*sum1*const
        pref=pref/math.sqrt(len(nlist))
        return pref  
   
    
    
    
    cpdef double complex Vm_v(s, int n, int l, int m, double p, double theta):
        '''
            Bound state prefactor in the velocity gauge <p|V|0> without exp(i*m*phi) for single nlm values
        '''
        cdef double complex pref=0.
        cdef double complex sum1=0.
        
        for k in range(n-l):
            term1=(-1)**k*math.factorial(n+l)*(2**k)*(2*s.Ip)**((-0.5-l)/2)*p**l
            term2=math.factorial(n-l-k-1)*math.factorial(2*l+k+1)*math.factorial(k)
            term3=math.gamma(2+k+2*l)/math.gamma(3/2+l)
            term4=hyp2f1(1+l+k/2,(3+k+2*l)/2, l+3/2,-p**2/(2*s.Ip))
            sum1=sum1+term1*term3*term4/term2
        term5=(-1j)**l*math.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
        
        if m>=0:
            Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l))) 
            Yml=Yml*lpmn(m,l,cos_re(theta))[0][m][l]  
        else: 
            Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l)))
            Yml=Yml*(-1)**abs(m)*(math.factorial(l-abs(m))/math.factorial(l+abs(m)))*lpmn(abs(m),l,cos_re(theta))[0][abs(m)][l] 
            
        pref=term5*Yml*sum1
        return pref
    
    
    cpdef double complex Vm_l(s, int n, int l, int m, double p, double theta, dbl_or_cmplx ts):
        '''
            Bound state prefactor in the length gauge <p+A(t)|V|0> without exp(i*m*phi) for single nlm values
        '''
        cdef double pz=p*cos_re(theta)
        cdef double sn = 1 if s.Af(ts).imag > 0 else -1
        cdef double complex pz_2 = 1j*sn*sqrt_re(2*s.Ip + p**2 - pz**2) #tilde{pz}
        cdef double complex pmod=  1j*sqrt_re(2*s.Ip) #=modulus{tilde{p}}
        cdef double complex pE=pz_2*s.Ef(ts) #=tilde{p}*E(ts)
        cdef double complex sum1
        cdef dbl_or_cmplx term1
        cdef double complex Yml
        
        cdef double complex pref=0.
        
        sum1=0. 
        for k in range(n-l):
            term2=(-1)**k*math.factorial(n+l)*(2**k)*(2*s.Ip)**(3/4+l/2+k/2)*pmod**l
            term3=(1j*pE)**((2+k+2*l)/4)
            term4=math.factorial(n-l-k-1)*math.factorial(k)*math.gamma(3/2+l)*math.gamma(1+l+k/2)            
            if k==0:
                term5=math.gamma((2+2*l)/4)              
            elif k==1:
                F1=complex(mp.hyper([-1/4.,1/4.],[1/2.,(5+2*l)/4.],s.Ip**2/(1j*pE)))
                F2=complex(mp.hyper([1/4.,3/4.],[3/2.,(7+2*l)/4.],s.Ip**2/(1j*pE)))
                term5=math.gamma((3+2*l)/4)*F1-sqrt(4*s.Ip**2/(1j*pE))*math.gamma((5+2*l)/4)*F2/(3+2*l)
            else:
                return np.nan
            sum1=sum1+term2*term5/(term3*term4)   
            
        term1=(-1j)**l*math.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
            
        if m>=0:
            Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l))) 
            Yml=Yml*lpmn(m,l,pz_2.imag/pmod.imag)[0][m][l]  
        else: 
            #Yml=sqrt((2*l+1)*math.factorial(l-m)/(4*np.pi*math.factorial(m+l)))
            #Yml=Yml*(-1)**abs(m)*(math.factorial(l-abs(m))/math.factorial(l+abs(m)))*clpmn(abs(m),l,pz_2/pmod)[0][abs(m)][l] 
            Yml=sqrt((2*l+1)*math.factorial(l-abs(m))/(4*np.pi*math.factorial(abs(m)+l)))
            Yml=Yml*lpmn(abs(m),l,pz_2.imag/pmod.imag)[0][abs(m)][l]  
            Yml=((-1)**abs(m))*np.conj(Yml)           
        pref=term1*Yml*sum1

        return pref   
    
    
    
    

    @cython.boundscheck(False) # turn off bounds-checking for entire function  
    @cython.wraparound(False)  # turn off negative index wrapping for entire function
    #@functools.lru_cache(maxsize=cacheSize)
    cpdef double complex M(s, double p, double theta, double phi, double tf = np.inf):
        '''
            Final transition amplitude
            Constructed as sum 
        '''
        #eLim is 1 or 2, the number of orbits accounted for
        cdef double complex MSum = 0.
        times = s.TimesGen(p, theta, phi)
        
        # single saddle point (maximum of E**2)
        #Emax2 = np.amax([np.abs(s.Ef(ts))**2 for ts in times])
        #for ts in times:
        #    if np.abs(s.Ef(ts))**2 == Emax2:
        #        ts_max = ts
        #times = [ts_max]
        # comment the previous lines to compute for all saddle points
        Efs = [s.Ef(real(ts)) for ts in times]
        its = np.argmax(np.abs(Efs))
        for ts in times[its:its+1]:
            if(real(ts)<tf) :#and -s.Ef(real(ts))>0. :
                det = sqrt(2*Pi*I1/s.DDS(p, theta, phi, ts))
                expS = exp(I1*s.S(p, theta, phi, ts))
                if s.Gauge==0:
                    d0 = s.d0_v(s.n,s.l,s.m,s.const1,p,theta,phi)
                elif s.Gauge==1:
                    d0 = s.d0_l(s.n,s.l,s.m,s.const1,p,theta,phi,ts)
                else:
                    return np.nan
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
    
    def Myz_List(s, double px, pyList, pzList, tf = np.inf):
        return np.array([s.Mxy(px, py, pz, tf) for py, pz in zip(pyList, pzList)])
    
    
    
    
    
    cpdef Ml(s, double p, double theta, int Nphi = 250):
        '''
            This is the fourier series coefficient of M to get the OAM distribusion.
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
    
    

    
    
    cpdef double complex Ml_av(s, double p, double theta, int n, int l, int OAM, double tf = np.inf):
        '''
            OAM transition amplitude taking advantage of the conservation law OAM=m
        '''
        phi = 0. #redundant for linearly polarized field
        
        cdef double complex MSum = 0.
        times = s.TimesGen(p, theta, phi)

        for ts in times:
            if(real(ts)<tf):
                det = sqrt(2*Pi*I1/s.DDS(p, theta, phi, ts))
                expS = exp(I1*s.S(p, theta, phi, ts))
                if s.Gauge==0:
                    d0 = s.Vm_v(n,l,OAM,p,theta)
                elif s.Gauge==1:
                    d0 = s.Vm_l(n,l,OAM,p,theta,ts)
                else:
                    return np.nan
                MSum += d0*det*expS
        MSum = MSum*1j**OAM
        return MSum
    

    cpdef Mlxz_av(s, px, pz, int n, int l, int OAM):
        '''
        convert Ml_av to cartesian coordinates
        '''
        cdef double p = sqrt_re(px*px +pz*pz)
        cdef double theta = acos_re(pz/p)
        return s.Ml_av(p, theta, n, l, OAM)
    

    cpdef coeff_av(s, int l, int m, int enant):
        '''
        Wave function coefficient depending on l,m and the enantiomer
            enant=0 (enantiomer +)
            enant=1 (enantiomer -)
        '''   
        if enant==0:
            if l==2 and m==1:
                coeff=1
            elif l== 3 and m==1:
                coeff=1j
            elif l==2 and m==-1:
                coeff=-1
            elif l==3 and m==-1:
                coeff=1j
            else:
                coeff=0.          
        elif enant==1:    
            if l==2 and m==-1:
                coeff=1
            elif l== 3 and m==-1:
                coeff=1j
            elif l==2 and m==1:
                coeff=-1
            elif l==3 and m==1:
                coeff=1j
            else:
                coeff=0.
        else:
            coeff=np.nan      
        return coeff/2.

    
    cpdef aligned_av(s, double px, double pz, int OAM, int enant=0):
        cdef double p = sqrt_re(px*px +pz*pz)
        cdef double theta = acos_re(pz/p)
        n = 4
        Msum = 0.
        N = 10
        alphaList = np.linspace(0, 2*Pi, N)
        
        for alpha in alphaList:
            for l,m in zip(s.l,s.m):
                if(m==OAM):
                    Msum=Msum+s.coeff_av(l,m,enant)*s.Ml_av(p, theta, n, l, OAM)
                    
        return Msum
        
 
    cpdef orientation_av(s, double px, double pz, int OAM, int enant=0):
        '''
        Orientation averaging of f(rho)|M_l|^2 with f(rho)=3*cos^2(beta) 
        '''      
        n=4
        cdef double complex integral = 0.
            
        for l,m in zip(s.l, s.m):
            integral=integral+abs(s.Mlxz_av(px, pz, n, l, OAM))**2/(2*l+1)
            m2=OAM
            for l2 in s.l:
                term1=np.conj(s.coeff_av(l2, m, enant))*s.coeff_av(l, m, enant)*(-1)**(m2-m)
                if abs(m2)<=l and abs(m2)<=l2:
                    term2=np.conj(s.Mlxz_av(px, pz, n, l2, m2))*s.Mlxz_av(px, pz, n, l, m2)
                else:
                    term2=0.
                term3=float(wigner_3j(2, l, l2, 0, m2, -m2))*float(wigner_3j(2, l, l2, 0, m, -m))
                integral=integral+2.*term1*term2*term3
           
        return integral
    
    
    def aligned_av_list(s, pxList, pzList, int OAM, int enant=0):
        return np.array([abs(s.aligned_av(px, pz, OAM, enant))**2 for px, pz in zip(pxList, pzList)])
        
    def orientation_av_list(s, pxList, pzList, int OAM, int enant=0):
        return np.array([abs(s.orientation_av(px, pz, OAM, enant))**2 for px, pz in zip(pxList, pzList)])
 


 #####   ---   code for spectra
    cpdef double Spectra(s, double E, double phi = 0., double t = np.inf, double err = 1.0e-4, int limit = 500):
        '''Function to check the SFA spectrum looks OK'''    
        Norm_val, Norm_error = it.quad(s.Spec_Norm, 0, Pi, args = (phi, E, t), epsabs=err, epsrel=err, limit=limit  )    
        return Norm_val
    
    cpdef double Spectra2(s, double E, double theta1, double theta2, double phi = 0., double t = np.inf, double err = 1.0e-4, int limit = 500):
        '''Function to check the SFA spectrum looks OK (with limits of integration)'''   
        Norm_val, Norm_error = it.quad(s.Spec_Norm, theta1, theta2, args = (phi, E, t), epsabs=err, epsrel=err, limit=limit  )    
        return Norm_val
     
    #Spectra Norm integrand
    cpdef double Spec_Norm(s, double theta, double phi, double E, double t = np.inf):
        '''Function to compute the integrand of the theta integral for the classical 'spectrum' Fisher info'''
        #phrased in spherical coordinates
        #cdef double complex M, Mg
        cdef double px, pr
        pr = sqrt_re(2*E)
        px = pr*sin_re(theta)       
        return abs(s.M(pr, theta, phi, t))**2 *px * 2*Pi
    
    cpdef double M_integrand(s, double theta, double E, int index_OAM, int Nphi=9):  
        cdef double px, pr
        pr = sqrt_re(2*E)
        px = pr*sin_re(theta) 
        return abs(s.Ml(pr, theta, Nphi)[index_OAM])**2 *px*2*Pi
    
    cpdef double M_integration(s, double E, double theta1, double theta2, int index_OAM, double err = 1.0e-4, int limit = 500, int Nphi=9):
        '''Function to compute the integrand of the theta integral for the OAM transition amplitude'''
        Norm_val, Norm_error = it.quad(s.M_integrand, theta1, theta2, args = (E,index_OAM,Nphi), epsabs=err, epsrel=err, limit=limit  )   
        return Norm_val
    

    
    
    cpdef double M_integrand_align(s, double theta, double E, int OAM, int enant=0):  
        cdef double px, pr
        pr = sqrt_re(2*E)
        px = pr*sin_re(theta) 
        pz = pr*cos_re(theta)
        return abs(s.aligned_av(px, pz, OAM, enant))**2 *px*2*Pi
    
    cpdef double M_integration_align(s, double E, double theta1, double theta2, int OAM, double err = 1.0e-4, int limit = 500, int enant=0):
        '''Function to compute the integrand of the theta integral for the OAM transition amplitude aligned'''
        Norm_val, Norm_error = it.quad(s.M_integrand_align, theta1, theta2, args = (E,OAM,enant), epsabs=err, epsrel=err, limit=limit  )   
        return Norm_val    
    
    
    
    cpdef double M_integrand_av(s, double theta, double E, int OAM, int enant=0):  
        cdef double px, pr
        pr = sqrt_re(2*E)
        px = pr*sin_re(theta) 
        pz = pr*cos_re(theta)
        return abs(s.orientation_av(px, pz, OAM, enant))**2 *px*2*Pi
    
    cpdef double M_integration_av(s, double E, double theta1, double theta2, int OAM, double err = 1.0e-4, int limit = 500, int enant=0):
        '''Function to compute the integrand of the theta integral for the OAM transition amplitude orientation averaged'''
        Norm_val, Norm_error = it.quad(s.M_integrand_av, theta1, theta2, args = (E,OAM,enant), epsabs=err, epsrel=err, limit=limit  )   
        return Norm_val
    
    
    
    
#EXACT INTEGRATION 
####################################################################################################################################
    cpdef double complex d0_vl(s, nlist, llist, mlist, constlist, double p, double theta, double phi, double ts, int OAM):   
        '''prefactor both gauges exact integration''' 
        cdef double px=p*sin_re(theta)*cos_re(phi), py=p*sin_re(theta)*sin_re(phi), pz=p*cos_re(theta)
        cdef double pmod = sqrt_re(px**2+py**2+(pz+s.Af(ts))**2)
        cdef double pz_2 = pz+s.Af(ts)
        cdef double theta2 = acos_re(pz_2/pmod)
        
        cdef double complex pref=0.
        
        if s.Gauge==0:
            pmod = p
            theta2 = theta
    
        for n,l,m,const in zip(nlist,llist,mlist,constlist):
            if m==OAM:
                sum1=0.
                for k in range(n-l):
                    term1=(-1)**k*math.factorial(n+l)*(2**k)*(2*s.Ip)**((-0.5-l)/2)*pmod**l
                    term2=math.factorial(n-l-k-1)*math.factorial(2*l+k+1)*math.factorial(k)
                    term3=math.gamma(2+k+2*l)/math.gamma(3/2+l)
                    term4=hyp2f1(1+l+k/2,(3+k+2*l)/2, l+3/2,-pmod**2/(2*s.Ip))
                    sum1=sum1+term1*term3*term4/term2
                term5=(-1j)**l*math.sqrt(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
                Yml=sph_harm(m,l,phi,theta2)
                pref=pref+term5*Yml*sum1*const
        pref=pref/math.sqrt(len(nlist))
        return pref
    

    
    cpdef double S_exact_integrand(s, double tau, double p, double theta):
        cdef double S_0 = p*cos_re(theta)*s.Af(tau) + 0.5*(s.Af(tau))**2
        return S_0
    
    cpdef double S_exact(s, double t, double p, double theta):
        '''
            Action as given by the SFA with exact integration
        '''
        S_1 = (s.Ip + 0.5*p*p)*t 
        S_integral,S_error = it.quad(s.S_exact_integrand, 0., t, args=(p,theta), limit=1000)
        return S_1+S_integral
    
    cpdef double M_exact_integrand(s, double tau, double p, double theta, double phi, int x, int OAM):       
        d0 = s.d0_vl(s.n,s.l,s.m,s.const1,p,theta,phi,tau,OAM)  
        M_0 = d0*exp(1j*s.S_exact(tau, p, theta))
        if x==0:
            return M_0.real
        return M_0.imag
    
    cpdef double complex M_exact(s, double p, double theta, double phi, double t, int OAM):
        '''
            Transition amplitude M with exact integration
        '''   
        cdef double complex bound1 = 1./( 1j*(s.Ip + 0.5*p*p))
        cdef double complex bound2 = np.exp(1j*s.N*Pi*(8*s.Ip + 4*p*p + 3*s.Up)/(4*s.omega))/(1j*(s.Ip + 0.5*p*p))
        
        if s.Gauge==0:
            d0 = s.d0_vl(s.n,s.l,s.m,s.const1,p,theta,phi,0.,OAM)
            d1 = d0
            d2 = d0
        else:
            d1 = s.d0_vl(s.n,s.l,s.m,s.const1,p,theta,phi,0.,OAM)
            d2 = s.d0_vl(s.n,s.l,s.m,s.const1,p,theta,phi,2.*s.N*np.pi/s.omega,OAM)
        
        M_integral1,M_error1 = it.quad(s.M_exact_integrand, 0., t, args=(p,theta,phi,0,OAM), limit=1000)
        M_integral2,M_error2 = it.quad(s.M_exact_integrand, 0., t, args=(p,theta,phi,1,OAM), limit=1000)
    
        return M_integral1 + 1j*M_integral2 + d1*bound1 - d2*bound2
    
    
    cpdef double complex Mxz_exact(s, double px, double pz, double phi, double t, int OAM):
        '''
        convert Ml to cartesian coordinates, note px is the perpendicular coordinate st. px^2 = px^2+py^2
        '''
        cdef double p = sqrt_re(px*px +pz*pz)
        cdef double theta = acos_re(pz/p)
        return s.M_exact(p, theta, phi, t, OAM)

    def Mxz_exact_list(s, pxList, pzList, double t, int OAM):
        return np.array([s.Mxz_exact(px, pz, 0., t, OAM) for px, pz in zip(pxList, pzList)])
    
    #cpdef double complex Ml_exact(s, double px, double pz, double t, int OAM):
        
####################################################################################################################################
    


        
  