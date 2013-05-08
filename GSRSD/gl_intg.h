#ifndef GL_INTG_H
#define GL_INTG_H

////////////////////////////////////////////////////////////
// Gauss-Legendre quadrature

// const int integral::gl_num = 32;
// const double integral::x[ 32 ]
// = { -0.99726386, -0.98561151, -0.96476226, 
// 	-0.93490608, -0.89632116, -0.84936761, 
// 	-0.7944838 , -0.73218212, -0.66304427, 
// 	-0.58771576, -0.50689991, -0.42135128, 
// 	-0.3318686 , -0.23928736, -0.14447196, 
// 	-0.04830767,  0.04830767,  0.14447196,  
// 	0.23928736,  0.3318686 , 0.42135128,  
// 	0.50689991,  0.58771576,  0.66304427,  
// 	0.73218212, 0.7944838 ,  0.84936761,  
// 	0.89632116,  0.93490608,  0.96476226, 
// 	0.98561151,  0.99726386 };
// const double integral::w[ 32 ]
// = { 0.00701861,  0.01627439,  0.02539207,  
// 	0.03427386,  0.0428359 , 0.05099806,  
// 	0.05868409,  0.06582222,  0.07234579, 
// 	0.0781939 , 0.08331192,  0.08765209,  
// 	0.09117388,  0.0938444 ,  0.09563872, 
// 	0.09654009,  0.09654009,  0.09563872, 
// 	0.0938444 ,  0.09117388, 0.08765209,  
// 	0.08331192,  0.0781939 ,  0.07234579,  
// 	0.06582222, 0.05868409,  0.05099806,  
// 	0.0428359 ,  0.03427386,  0.02539207, 
// 	0.01627439,  0.00701861 };


const int integral::gl_num = 128;
const double integral::x[ 128 ] =
{-0.9998248879471319,-0.9990774599773759,-0.997733248625514,
 -0.9957927585349812,-0.9932571129002129,-0.9901278184917344,
 -0.9864067427245862,-0.9820961084357185,-0.9771984914639074,
 -0.9717168187471366,-0.9656543664319653,-0.9590147578536999,
 -0.9518019613412644,-0.9440202878302202,-0.9356743882779164,
 -0.9267692508789478,-0.9173101980809605,-0.9073028834017568,
 -0.8967532880491582,-0.8856677173453972,-0.8740527969580318,
 -0.8619154689395485,-0.849262987577969,-0.8361029150609068,
 -0.8224431169556438,-0.8082917575079137,-0.7936572947621933,
 -0.778548475506412,-0.7629743300440947,-0.746944166797062,
 -0.7304675667419088,-0.7135543776835874,-0.6962147083695143,
 -0.6784589224477193,-0.6602976322726461,-0.6417416925623076,
 -0.6228021939105849,-0.6034904561585486,-0.5838180216287631,
 -0.5637966482266181,-0.5434383024128104,-0.5227551520511755,
 -0.5017595591361445,-0.480464072404172,-0.4588814198335522,
 -0.4370245010371042,-0.414906379552275,-0.3925402750332674,
 -0.369939555349859,-0.3471177285976355,-0.3240884350244134,
 -0.3008654388776772,-0.2774626201779044,-0.2538939664226943,
 -0.23017356422666,-0.2063155909020792,-0.1823343059853372,
 -0.1582440427142249,-0.1340591994611878,-0.1097942311276437,
 -0.0854636405045155,-0.06108196960413957,-0.03666379096873349,
 -0.01222369896061576,0.01222369896061576,0.03666379096873349,
 0.06108196960413957,0.0854636405045155,0.1097942311276437,
 0.1340591994611878,0.1582440427142249,0.1823343059853372,
 0.2063155909020792,0.23017356422666,0.2538939664226943,
 0.2774626201779044,0.3008654388776772,0.3240884350244134,
 0.3471177285976355,0.369939555349859,0.3925402750332674,
 0.414906379552275,0.4370245010371042,0.4588814198335522,
 0.480464072404172,0.5017595591361445,0.5227551520511755,
 0.5434383024128104,0.5637966482266181,0.5838180216287631,
 0.6034904561585486,0.6228021939105849,0.6417416925623076,
 0.6602976322726461,0.6784589224477193,0.6962147083695143,
 0.7135543776835874,0.7304675667419088,0.746944166797062,
 0.7629743300440947,0.778548475506412,0.7936572947621933,
 0.8082917575079137,0.8224431169556438,0.8361029150609068,
 0.849262987577969,0.8619154689395485,0.8740527969580318,
 0.8856677173453972,0.8967532880491582,0.9073028834017568,
 0.9173101980809605,0.9267692508789478,0.9356743882779164,
 0.9440202878302202,0.9518019613412644,0.9590147578536999,
 0.9656543664319653,0.9717168187471366,0.9771984914639074,
 0.9820961084357185,0.9864067427245862,0.9901278184917344,
 0.9932571129002129,0.9957927585349812,0.997733248625514,
 0.9990774599773759,0.9998248879471319};
const double integral::w[ 128 ]=
{0.0004493809602920904,0.001045812679340349,
 0.00164250301866903,0.002238288430962619,0.002832751471457991,
 0.003425526040910216,0.004016254983738642,0.004604584256702955,
 0.00519016183267633,0.005772637542865699,0.006351663161707189,
 0.006926892566898814,0.007497981925634729,0.008064589890486058,
 0.00862637779861675,0.009183009871660874,0.009734153415006806,
 0.01027947901583216,0.01081866073950308,0.01135137632408042,
 0.01187730737274028,0.01239613954395092,0.01290756273926735,
 0.01341127128861633,0.01390696413295199,0.01439434500416685,
 0.01487312260214731,0.01534301076886514,0.01580372865939935,
 0.01625500090978519,0.0166965578015892,0.01712813542311138,
 0.0175494758271177,0.01796032718500869,0.01836044393733134,
 0.01874958694054471,0.01912752360995095,0.0194940280587066,
 0.01984888123283086,0.02019187104213004,0.02052279248696007,
 0.02084144778075115,0.02114764646822135,0.02144120553920846,
 0.02172194953805208,0.02198971066846049,0.02224432889379977,
 0.02248565203274497,0.02271353585023646,0.02292784414368685,
 0.02312844882438703,0.02331522999406276,0.02348807601653591,
 0.02364688358444762,0.0237915577810034,0.02392201213670346,
 0.02403816868102405,0.02413995798901928,0.02422731922281525,
 0.02430020016797187,0.02435855726469063,0.02440235563384958,
 0.02443156909785005,0.02444618019626252,0.02444618019626252,
 0.02443156909785005,0.02440235563384958,0.02435855726469063,
 0.02430020016797187,0.02422731922281525,0.02413995798901928,
 0.02403816868102405,0.02392201213670346,0.0237915577810034,
 0.02364688358444762,0.02348807601653591,0.02331522999406276,
 0.02312844882438703,0.02292784414368685,0.02271353585023646,
 0.02248565203274497,0.02224432889379977,0.02198971066846049,
 0.02172194953805208,0.02144120553920846,0.02114764646822135,
 0.02084144778075115,0.02052279248696007,0.02019187104213004,
 0.01984888123283086,0.0194940280587066,0.01912752360995095,
 0.01874958694054471,0.01836044393733134,0.01796032718500869,
 0.0175494758271177,0.01712813542311138,0.0166965578015892,
 0.01625500090978519,0.01580372865939935,0.01534301076886514,
 0.01487312260214731,0.01439434500416685,0.01390696413295199,
 0.01341127128861633,0.01290756273926735,0.01239613954395092,
 0.01187730737274028,0.01135137632408042,0.01081866073950308,
 0.01027947901583216,0.009734153415006806,0.009183009871660874,
 0.00862637779861675,0.008064589890486058,0.007497981925634729,
 0.006926892566898814,0.006351663161707189,0.005772637542865699,
 0.00519016183267633,0.004604584256702955,0.004016254983738642,
 0.003425526040910216,0.002832751471457991,0.002238288430962619,
 0.00164250301866903,0.001045812679340349,0.0004493809602920904};

#endif

