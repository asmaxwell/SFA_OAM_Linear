/*
 * HpgForms.h
 *
 *  Created on: 1 Aug 2018
 *      Author: andy
 * edited for cython 10 May 2021
 * by Andy
 */

#ifndef HPGFORMS_H_
#define HPGFORMS_H_

#include <complex>
#include <cmath>


#define HpgHe (dcmplx(0.,-6.641409757794474e-6)*std::exp(0.036692669538479605*Ip) - \
     dcmplx(0.,0.0009735680519733292)*std::exp(0.25008127641483485*Ip) - \
     dcmplx(0.,0.028199729381134114)*std::exp(1.3055068891598542*Ip))*Eft*(Aft + pz)

#define HpgNe Eft*(std::exp(1.1215313838127132*Ip)* \
      (0.02017174338963387 - 0.022623243277691024*std::pow(Aft,2) - \
        0.04524648655538205*Aft*pz - 0.022623243277691024*std::pow(pz,2)) + \
     std::exp(0.294764047275437*Ip)* \
      (0.004106740585108807 - 0.0012105194759769683*std::pow(Aft,2) - \
        0.0024210389519539366*Aft*pz - 0.0012105194759769683*std::pow(pz,2)) + \
     std::exp(0.08194363752723602*Ip)* \
      (0.00040245338020714636 - 0.0000329784939093053*std::pow(Aft,2) - \
        0.0000659569878186106*Aft*pz - 0.0000329784939093053*std::pow(pz,2)) + \
     std::exp(0.018845075052584356*Ip)* \
      (0.000013186874156430686 - 2.485076331869213e-7*std::pow(Aft,2) - \
        4.970152663738426e-7*Aft*pz - 2.485076331869213e-7*std::pow(pz,2)))

#define HpgAr Eft*(std::exp(0.11025552821218457*Ip)*(0.04480543385656865 - 0.004940046776632074*std::pow(Aft,2) - \
	        0.009880093553264148*Aft*pz - 0.004940046776632074*std::pow(pz,2)) + \
	     std::exp(0.03096492912127724*Ip)*(0.007016075084433821 - 0.00021725226769905251*std::pow(Aft,2) - \
	        0.00043450453539810503*Aft*pz - 0.00021725226769905251*std::pow(pz,2)) + \
	     std::exp(0.00713823768332779*Ip)*(0.0002624381584861572 - 1.8733459524490377e-6*std::pow(Aft,2) - \
	        3.7466919048980755e-6*Aft*pz - 1.8733459524490377e-6*std::pow(pz,2)) + \
	     std::exp(0.3242352910660208*Ip)*(-0.1625979268437463 + 0.052719986136913644*std::pow(Aft,2) + \
	        0.10543997227382729*Aft*pz + 0.052719986136913644*std::pow(pz,2)) + \
	     std::exp(0.8233610586447148*Ip)*(-1.6101222193844695 + 1.3257119350997746*std::pow(Aft,2) + 2.651423870199549*Aft*pz + \
	        1.3257119350997746*std::pow(pz,2)) + std::exp(2.5592072599591553*Ip)* \
	      (-6.064387405128389 + 15.52002427440944*std::pow(Aft,2) + 31.04004854881888*Aft*pz + 15.52002427440944*std::pow(pz,2)))

#define HpgArEx_4S (dcmplx(0.,-9.505713243048165e-11)*std::exp(0.00032181037645377837*\
				Ip) - dcmplx(0.,1.5332932982279847e-8)*std::exp(0.002130578920904388*\
				Ip) - dcmplx(0.,2.1110945924751829e-7)*std::exp(0.00713823768332779*\
				Ip) - dcmplx(0.,3.976195722263694e-7)*std::exp(0.009801596091907607*\
				Ip) + dcmplx(0.,1.103035267151604e-6)*std::exp(0.03096492912127724*Ip)\
				 + dcmplx(0.,0.00026063575580827166)*std::exp(0.11025552821218457*Ip) \
				+ dcmplx(0.,0.0012281510890104434)*std::exp(0.3242352910660208*Ip) - \
				dcmplx(0.,0.04068275968739755)*std::exp(0.8233610586447148*Ip) + \
				dcmplx(0.,0.2439288298976548)*std::exp(2.5592072599591553*Ip))*Eft*(\
				Aft + pz)


#define HpgXe Eft*(std::exp(3.7543456550957734*Ip)* \
      (0.10997016150990008 - 0.41286599805487384*std::pow(Aft,2) -  \
        0.8257319961097477*Aft*pz - 0.41286599805487384*std::pow(pz,2)) +  \
     std::exp(1.3330453955278994*Ip)* \
      (0.033207920689466874 - 0.044267665770149484*std::pow(Aft,2) -  \
        0.08853533154029897*Aft*pz - 0.044267665770149484*std::pow(pz,2)) +  \
     std::exp(0.3167616318038815*Ip)* \
      (0.00015713290167146583 - 0.00004977367434353236*std::pow(Aft,2) -  \
        0.00009954734868706472*Aft*pz - 0.00004977367434353236*std::pow(pz,2)) +  \
     std::exp(0.05816255643512851*Ip)* \
      (0.00016602610493280761 - 9.656502697858989e-6*std::pow(Aft,2) -  \
        0.000019313005395717977*Aft*pz - 9.656502697858989e-6*std::pow(pz,2)) +  \
     std::exp(0.017475435653873134*Ip)* \
      (0.00003099538602410666 - 5.416578740312345e-7*std::pow(Aft,2) -  \
        1.083315748062469e-6*Aft*pz - 5.416578740312345e-7*std::pow(pz,2)) +  \
     std::exp(0.0007422129613836987*Ip)* \
      (-4.578520789342088e-8 + 3.39823747381442e-11*std::pow(Aft,2) +  \
        6.79647494762884e-11*Aft*pz + 3.39823747381442e-11*std::pow(pz,2)) +  \
     std::exp(0.0034282007910916147*Ip)* \
      (-1.2459386058071205e-6 + 4.271327714079553e-9*std::pow(Aft,2) +  \
        8.542655428159106e-9*Aft*pz + 4.271327714079553e-9*std::pow(pz,2)) +  \
     std::exp(0.008450582346530663*Ip)* \
      (-3.4689557741268297e-6 + 2.9314696425731796e-8*std::pow(Aft,2) +  \
        5.862939285146359e-8*Aft*pz + 2.9314696425731796e-8*std::pow(pz,2)) +  \
     std::exp(0.01050926202280594*Ip)* \
      (-3.110807815614424e-6 + 3.269229443688457e-8*std::pow(Aft,2) +  \
        6.538458887376915e-8*Aft*pz + 3.269229443688457e-8*std::pow(pz,2)) +  \
     std::exp(0.05787923982795511*Ip)* \
      (-5.432004909053535e-6 + 3.1440031487773904e-7*std::pow(Aft,2) + \
        6.288006297554781e-7*Aft*pz + 3.1440031487773904e-7*std::pow(pz,2)) + \
     std::exp(0.14439107108718968*Ip)* \
      (-0.0009413785576391259 + 0.00013592665823602715*std::pow(Aft,2) + \
        0.0002718533164720543*Aft*pz + 0.00013592665823602715*std::pow(pz,2)) + \
     std::exp(0.3568777479586593*Ip)* \
      (-0.0028830950306471394 + 0.001028912461688153*std::pow(Aft,2) + \
        0.002057824923376306*Aft*pz + 0.001028912461688153*std::pow(pz,2)))

//Abbie's Hpg for He alpha
#define HpgHeTheta Eft*(dcmplx(0.,-6.641409757794473e-6)*Aft*std::exp(0.\
036692669538479605*Ip) - \
dcmplx(0.,0.0009735680519733291)*Aft*std::exp(0.25008127641483485*Ip) \
- dcmplx(0.,0.028199729381134117)*Aft*std::exp(1.3055068891598542*Ip) \
- dcmplx(0.,6.641409757794473e-6)*std::exp(0.036692669538479605*Ip)*\
pz - dcmplx(0.,0.0009735680519733291)*std::exp(0.25008127641483485*Ip)\
*pz - dcmplx(0.,0.028199729381134117)*std::exp(1.3055068891598542*Ip)*\
pz + (dcmplx(0.,-8.470329472543003e-22)*Aft*std::exp(0.\
036692669538479605*Ip) - \
dcmplx(0.,1.0842021724855044e-19)*Aft*std::exp(0.25008127641483485*Ip)\
 - dcmplx(0.,8.470329472543003e-22)*std::exp(0.036692669538479605*Ip)*\
pz - dcmplx(0.,1.0842021724855044e-19)*std::exp(0.25008127641483485*\
Ip)*pz)*std::cos(2*alpha) + \
(dcmplx(0.,1.6940658945086007e-21)*std::exp(0.036692669538479605*Ip) \
+ dcmplx(0.,1.0842021724855044e-19)*std::exp(0.25008127641483485*Ip) \
+ dcmplx(0.,3.469446951953614e-18)*std::exp(1.3055068891598542*Ip))*\
py*std::cos(alpha)*std::sin(alpha))
///Diatomic Molecules
#define HpgN2Par std::exp(dcmplx(0.,-1.034)*Aft + 3.4024398500511315*Ip - dcmplx(0.,1.034)*pz)*Eft*\
	(std::exp(dcmplx(0.,2.068)*Aft - 3.4023200468578754*Ip + dcmplx(0.,2.068)*pz)*\
      (2.356034686819464e-9 - dcmplx(0.,2.729791865606503e-13)*Aft -\
        dcmplx(0.,2.729791865606503e-13)*pz) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.4016429838530415*Ip + dcmplx(0.,2.068)*pz)*\
      (7.443007105720177e-8 - dcmplx(0.,5.736054907825277e-11)*Aft -\
        dcmplx(0.,5.736054907825277e-11)*pz) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.398940950983222*Ip + dcmplx(0.,2.068)*pz)*\
      (1.1064597440681432e-6 - dcmplx(0.,3.744091844486731e-9)*Aft -\
        dcmplx(0.,3.744091844486731e-9)*pz) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.390012651377084*Ip + dcmplx(0.,2.068)*pz)*\
      (9.69319005563463e-6 - dcmplx(0.,1.1649825784011149e-7)*Aft -\
        dcmplx(0.,1.1649825784011149e-7)*pz) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.363438937975802*Ip + dcmplx(0.,2.068)*pz)*\
      (0.000046162774639984356 - dcmplx(0.,1.7411898596588914e-6)*Aft -\
        dcmplx(0.,1.7411898596588914e-6)*pz) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.2885559701548934*Ip + dcmplx(0.,2.068)*pz)*\
      (0.0000792403219787326 - dcmplx(0.,8.727461616213944e-6)*Aft -\
        dcmplx(0.,8.727461616213944e-6)*pz) +\
     (-2.3560346868194634e-9 - dcmplx(0.,2.729791865606503e-13)*(Aft + pz))/\
      std::exp(3.4023200468578754*Ip) +\
     (-7.443007105720178e-8 - dcmplx(0.,5.736054907825277e-11)*(Aft + pz))/\
      std::exp(3.4016429838530415*Ip) +\
     (-1.106459744068143e-6 - dcmplx(0.,3.744091844486731e-9)*(Aft + pz))/\
      std::exp(3.398940950983222*Ip) +\
     (-9.69319005563463e-6 - dcmplx(0.,1.1649825784011147e-7)*(Aft + pz))/\
      std::exp(3.390012651377084*Ip) +\
     (-0.000046162774639984356 - dcmplx(0.,1.7411898596588916e-6)*(Aft + pz))/\
      std::exp(3.363438937975802*Ip) +\
     (-0.00007924032197873258 - dcmplx(0.,8.727461616213944e-6)*(Aft + pz))/\
      std::exp(3.2885559701548934*Ip) +\
     std::exp(dcmplx(0.,2.068)*Aft - 1.0442941072116412*Ip + dcmplx(0.,2.068)*pz)*\
      (-0.000017262100950051522 + dcmplx(0.,0.016398765940561884)*(Aft + pz) -\
        0.030946321689486737*std::pow(Aft + pz,2)) +\
     std::exp(dcmplx(0.,2.068)*Aft - 2.7549542954538544*Ip + dcmplx(0.,2.068)*pz)*\
      (0.004067519386121589 - dcmplx(0.,0.004744622431016149)*(Aft + pz) -\
        0.003489142289791188*std::pow(Aft + pz,2)) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.2183645706248574*Ip + dcmplx(0.,2.068)*pz)*\
      (0.0005649779843378503 - dcmplx(0.,0.0005192002098420466)*(Aft + pz) -\
        0.00009002322379284062*std::pow(Aft + pz,2)) +\
     std::exp(dcmplx(0.,2.068)*Aft - 3.3594341359679127*Ip + dcmplx(0.,2.068)*pz)*\
      (0.00003391671134182942 - dcmplx(0.,0.000017858490298164625)*(Aft + pz) -\
        7.127621096758362e-7*std::pow(Aft + pz,2)) +\
     (-0.00003391671134182942 - dcmplx(0.,0.000017858490298164625)*(Aft + pz) +\
        7.127621096758362e-7*std::pow(Aft + pz,2))/std::exp(3.3594341359679127*Ip) +\
     (-0.0005649779843378503 - dcmplx(0.,0.0005192002098420466)*(Aft + pz) +\
        0.00009002322379284062*std::pow(Aft + pz,2))/std::exp(3.2183645706248574*Ip) +\
     (-0.004067519386121589 - dcmplx(0.,0.004744622431016149)*(Aft + pz) +\
        0.003489142289791188*std::pow(Aft + pz,2))/std::exp(2.7549542954538544*Ip) +\
     (0.000017262100950051522 + dcmplx(0.,0.016398765940561884)*(Aft + pz) +\
        0.030946321689486737*std::pow(Aft + pz,2))/std::exp(1.0442941072116412*Ip))


#define HpgN2Theta (Eft*(-1.0079978422492269e-6*std::exp( \
         (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha))*std::sin(alpha)* \
        (dcmplx(0.,3.8298779193636014e-7)* \
           std::exp(0.00011980319325635685*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,0.00008047642867034601)* \
           std::exp(0.0007968661980898987*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,0.005252933333099794)* \
           std::exp(0.003498899067909383*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,0.16344619931199045)* \
           std::exp(0.012427198674047726*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,2.4428765727330592)* \
           std::exp(0.03900091207532979*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,1.0120139855195585)* \
           std::exp(0.04300571408321881*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,12.244564487558394)* \
           std::exp(0.1138838798962381*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,18.962458599442133)* \
           std::exp(0.1840752794262742*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
          dcmplx(0.,1160.7691485259431)* \
           std::exp(0.647485554597277*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
          dcmplx(0.,42045.041410520906)* \
           std::exp(2.3581457428394903*Ip +  \
             (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
          dcmplx(0.,3.8298779193636014e-7)* \
           std::exp(0.00011980319325635685*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,0.00008047642867034601)* \
           std::exp(0.0007968661980898987*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,0.005252933333099794)* \
           std::exp(0.003498899067909383*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,0.16344619931199045)* \
           std::exp(0.012427198674047726*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,2.4428765727330592)* \
           std::exp(0.03900091207532979*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,1.0120139855195585)* \
           std::exp(0.04300571408321881*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,12.244564487558394)* \
           std::exp(0.1138838798962381*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          dcmplx(0.,18.962458599442133)* \
           std::exp(0.1840752794262742*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) -  \
          dcmplx(0.,1160.7691485259431)* \
           std::exp(0.647485554597277*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) -  \
          dcmplx(0.,42045.041410520906)* \
           std::exp(2.3581457428394903*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
          (1.*std::exp(0.04300571408321881*Ip +  \
                (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
             126.30192117505115*std::exp( \
               0.1840752794262742*Ip +  \
                (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
             4895.240982125226*std::exp( \
               0.647485554597277*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
              43417.46182826849*std::exp( \
                2.3581457428394903*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
              1.*std::exp(0.04300571408321881*Ip +  \
                 dcmplx(0.,2.068)*py*std::sin(alpha)) -  \
              126.30192117505115*std::exp( \
                0.1840752794262742*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) -  \
              4895.240982125226*std::exp( \
                0.647485554597277*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) -  \
              43417.46182826849*std::exp( \
                2.3581457428394903*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)))* \
            (Aft + pz)*std::cos(alpha) +  \
           (-1.*std::exp(0.04300571408321881*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
              126.30192117505115*std::exp( \
                0.1840752794262742*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
              4895.240982125226*std::exp( \
                0.647485554597277*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) -  \
              43417.46182826849*std::exp( \
                2.3581457428394903*Ip +  \
                 (dcmplx(0.,2.068)*Aft + dcmplx(0.,2.068)*pz)*std::cos(alpha)) +  \
              1.*std::exp(0.04300571408321881*Ip +  \
                 dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
              126.30192117505115*std::exp( \
                0.1840752794262742*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
              4895.240982125226*std::exp( \
                0.647485554597277*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)) +  \
              43417.46182826849*std::exp( \
                2.3581457428394903*Ip + dcmplx(0.,2.068)*py*std::sin(alpha)))*py* \
            std::sin(alpha))*(py*std::cos(alpha) + (Aft + pz)*std::sin(alpha)) +  \
        std::cos(alpha)*(3.8605086787964696e-13* \
            std::exp(0.00011980319325635685*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-8630.821699279999 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           8.112006645163261e-11*std::exp( \
             0.0007968661980898987*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-1297.582959948 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           5.294945465243633e-9*std::exp( \
             0.003498899067909383*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-295.52152832400003 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           1.6475341623032342e-7*std::exp( \
             0.012427198674047726*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-83.20459237200001 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           2.46241431419611e-6*std::exp( \
             0.03900091207532979*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-26.512200484 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           1.0201079137297553e-6*std::exp( \
             0.04300571408321881*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-24.043316616000002 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           0.000012342494582740373* \
            std::exp(0.1138838798962381*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-9.079423716 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           0.000019114117351977965* \
            std::exp(0.1840752794262742*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-5.61726704 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) -  \
           0.001170052797063623*std::exp( \
             0.647485554597277*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-1.596946824 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) -  \
           0.04238131101908446*std::exp( \
             2.3581457428394903*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (-0.438480108 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) -  \
           0.04238131101908446*std::exp( \
             2.3581457428394903*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (0.438480108 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) -  \
           0.001170052797063623*std::exp( \
             0.647485554597277*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (1.596946824 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           0.000019114117351977965* \
            std::exp(0.1840752794262742*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (5.61726704 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           0.000012342494582740373* \
            std::exp(0.1138838798962381*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (9.079423716 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           1.0201079137297553e-6*std::exp( \
             0.04300571408321881*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (24.043316616000002 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           2.46241431419611e-6*std::exp( \
             0.03900091207532979*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (26.512200484 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           1.6475341623032342e-7*std::exp( \
             0.012427198674047726*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (83.20459237200001 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           5.294945465243633e-9*std::exp( \
             0.003498899067909383*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (295.52152832400003 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           8.112006645163261e-11*std::exp( \
             0.0007968661980898987*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (1297.582959948 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) +  \
           3.8605086787964696e-13*std::exp( \
             0.00011980319325635685*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (8630.821699279999 - I2*(Aft + pz)*std::cos(alpha) +  \
              I2*py*std::sin(alpha)) -  \
           0.0010484657792239186*std::exp( \
             0.04300571408321881*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (0.022355231885361385 -  \
              dcmplx(0.,0.0009614027107259082)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (-24.043316616000002 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) -  \
           0.0034937429904027474*std::exp( \
             0.1840752794262742*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (0.19796259438451572 -  \
              dcmplx(0.,0.03644001987727992)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (-5.61726704 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) -  \
           0.005835371835599359*std::exp( \
             0.647485554597277*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (1.3059755905801604 - dcmplx(0.,0.8456003295573017)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (-1.596946824 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) -  \
           0.0020445959330352387*std::exp( \
             2.3581457428394903*Ip +  \
              (dcmplx(0.,-1.034)*Aft - dcmplx(0.,1.034)*pz)*std::cos(alpha) +  \
              dcmplx(0.,2.068)*py*std::sin(alpha))* \
            (9.077074465271012 - dcmplx(0.,21.405064507715878)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (-0.438480108 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) +  \
           0.0020445959330352387*std::exp( \
             2.3581457428394903*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (9.077074465271012 - dcmplx(0.,21.405064507715878)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (0.438480108 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) +  \
           0.005835371835599359*std::exp( \
             0.647485554597277*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (1.3059755905801604 - dcmplx(0.,0.8456003295573017)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (1.596946824 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) +  \
           0.0034937429904027474*std::exp( \
             0.1840752794262742*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (0.19796259438451572 -  \
              dcmplx(0.,0.03644001987727992)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (5.61726704 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))) +  \
           0.0010484657792239186*std::exp( \
             0.04300571408321881*Ip +  \
              (dcmplx(0.,1.034)*Aft + dcmplx(0.,1.034)*pz)*std::cos(alpha))* \
            (0.022355231885361385 -  \
              dcmplx(0.,0.0009614027107259082)* \
               ((Aft + pz)*std::cos(alpha) - py*std::sin(alpha))* \
               (24.043316616000002 - I2*(Aft + pz)*std::cos(alpha) +  \
                 I2*py*std::sin(alpha))))))/ \
    (std::sqrt(2)*std::exp(dcmplx(0.,1.034)*py*std::sin(alpha)))

#define HpgO2Theta (Eft*(std::cos(alpha)*(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha))*(0.05811801018986199*std::exp(1.8518107005029516*\
					Ip + (dcmplx(0.,-1.1405)*Aft - dcmplx(0.,1.1405)*pz)*std::cos(alpha) \
					+ dcmplx(0.,2.281)*py*std::sin(alpha))*(dcmplx(0.,-0.6158836860000001)\
					 + (Aft + pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) - \
					0.05811801018986199*std::exp(1.8518107005029516*Ip + \
					(dcmplx(0.,1.1405)*Aft + \
					dcmplx(0.,1.1405)*pz)*std::cos(alpha))*(dcmplx(0.,0.6158836860000001) \
					+ (Aft + pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) + \
					0.002820383142833033*std::exp(0.4932124108025355*Ip + \
					(dcmplx(0.,-1.1405)*Aft - dcmplx(0.,1.1405)*pz)*std::cos(alpha) + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(dcmplx(0.,-2.312391122) + (Aft \
					+ pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) - \
					0.002820383142833033*std::exp(0.4932124108025355*Ip + \
					(dcmplx(0.,1.1405)*Aft + \
					dcmplx(0.,1.1405)*pz)*std::cos(alpha))*(dcmplx(0.,2.312391122) + (Aft \
					+ pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) + \
					0.00007612647322559938*std::exp(0.13889143523186814*Ip + \
					(dcmplx(0.,-1.1405)*Aft - dcmplx(0.,1.1405)*pz)*std::cos(alpha) + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(dcmplx(0.,-8.211449454) + (Aft \
					+ pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) - \
					0.00007612647322559938*std::exp(0.13889143523186814*Ip + \
					(dcmplx(0.,1.1405)*Aft + \
					dcmplx(0.,1.1405)*pz)*std::cos(alpha))*(dcmplx(0.,8.211449454) + (Aft \
					+ pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) + \
					5.912664369218131e-7*std::exp(0.032175827253389015*Ip + \
					(dcmplx(0.,-1.1405)*Aft - dcmplx(0.,1.1405)*pz)*std::cos(alpha) + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(dcmplx(0.,-35.445864096) + (Aft \
					+ pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) - \
					5.912664369218131e-7*std::exp(0.032175827253389015*Ip + \
					(dcmplx(0.,1.1405)*Aft + \
					dcmplx(0.,1.1405)*pz)*std::cos(alpha))*(dcmplx(0.,35.445864096) + \
					(Aft + pz)*std::cos(alpha) - 1.*py*std::sin(alpha)) + \
					dcmplx(0.,0.000054002434581551926)*std::exp(0.625*Ip + \
					(dcmplx(0.,-1.1405)*Aft - dcmplx(0.,1.1405)*pz)*std::cos(alpha) + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(1.2385397805018783 - \
					dcmplx(0.,0.774087362813674)*((Aft + pz)*std::cos(alpha) - \
					py*std::sin(alpha))*(-1.8248000000000002 - dcmplx(0,1)*(Aft + \
					pz)*std::cos(alpha) + dcmplx(0,1)*py*std::sin(alpha))) + \
					dcmplx(0.,0.000054002434581551926)*std::exp(0.625*Ip + \
					(dcmplx(0.,1.1405)*Aft + \
					dcmplx(0.,1.1405)*pz)*std::cos(alpha))*(1.2385397805018783 - \
					dcmplx(0.,0.774087362813674)*((Aft + pz)*std::cos(alpha) - \
					py*std::sin(alpha))*(1.8248000000000002 - dcmplx(0,1)*(Aft + \
					pz)*std::cos(alpha) + dcmplx(0,1)*py*std::sin(alpha)))) + \
					(std::sin(alpha)*(0.004968546164351748*std::exp(1.8518107005029516*Ip \
					+ (dcmplx(0.,2.281)*Aft + \
					dcmplx(0.,2.281)*pz)*std::cos(alpha))*(6.3166209753316265 - \
					11.697186313140495*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) - \
					0.004968546164351748*std::exp(1.8518107005029516*Ip + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(6.3166209753316265 - \
					11.697186313140495*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) + \
					dcmplx(0.,0.000054002434581551926)*std::exp(0.625*Ip + \
					(dcmplx(0.,2.281)*Aft + dcmplx(0.,2.281)*pz)*std::cos(alpha))*((Aft + \
					pz)*std::cos(alpha) - py*std::sin(alpha))*(1.2385397805018783 - \
					0.774087362813674*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) + \
					dcmplx(0.,0.000054002434581551926)*std::exp(0.625*Ip + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*((Aft + pz)*std::cos(alpha) - \
					py*std::sin(alpha))*(1.2385397805018783 - \
					0.774087362813674*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) + \
					0.006586173775571537*std::exp(0.4932124108025355*Ip + \
					(dcmplx(0.,2.281)*Aft + \
					dcmplx(0.,2.281)*pz)*std::cos(alpha))*(0.8682422763424836 - \
					0.42822786627555753*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) - \
					0.006586173775571537*std::exp(0.4932124108025355*Ip + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(0.8682422763424836 - \
					0.42822786627555753*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) + \
					0.004224322942512425*std::exp(0.13889143523186814*Ip + \
					(dcmplx(0.,2.281)*Aft + \
					dcmplx(0.,2.281)*pz)*std::cos(alpha))*(0.12974873512011034 - \
					0.018020988040351623*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) - \
					0.004224322942512425*std::exp(0.13889143523186814*Ip + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(0.12974873512011034 - \
					0.018020988040351623*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) + \
					0.0012701908762855004*std::exp(0.032175827253389015*Ip + \
					(dcmplx(0.,2.281)*Aft + \
					dcmplx(0.,2.281)*pz)*std::cos(alpha))*(0.014467201040401747 - \
					0.00046549416151601643*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2)) - \
					0.0012701908762855004*std::exp(0.032175827253389015*Ip + \
					dcmplx(0.,2.281)*py*std::sin(alpha))*(0.014467201040401747 - \
					0.00046549416151601643*std::pow(py*std::cos(alpha) + (Aft + \
					pz)*std::sin(alpha),2))))/std::exp(dcmplx(0.,1.1405)*(Aft + \
					pz)*std::cos(alpha))))/(std::sqrt(2)*std::exp(dcmplx(0.,1.1405)*py*\
					std::sin(alpha)))
using dcmplx = std::complex<double>;

enum Targets {
    He,
	HeTheta,
    Ne,
    Ar,
	ArEx_4S,
    Xe,
	N2,
	N2Theta,
	O2Theta,
	H,

};


dcmplx calculateHMatrixElement(int target, double Ip, dcmplx pz, double py,  dcmplx ts, dcmplx Eft, dcmplx Aft, double theta);


#endif /* HPGFORMS_H_ */
