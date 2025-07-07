LIBTARGET = nrd
LIBEXT = a #lib #a #a64
LIBTGT = $(LIBTARGET).$(LIBEXT)
HOME=/research/pcmd/Tsiky/Traj_Optimization_lstm
#HOME=Z:\Tsiky\Traj_Optimization
# HOME=/users/pierre/Documents
#HOME=d:
#HOME=z:\pierre\Documents

#ifeq ($(LIBEXT),a64) 
ARCH=x86_64
#else
#ARCH=i386
#endif

INCLUDEPATH=$(HOME)/Includes/nrd #"C:\Program Files\MATLAB\R2018b\sys\lcc64\lcc64\include64"
LIBPATH=$(HOME)/Libs
LIBOBJS=\
MNBRAK.o  chint.o    expdev.o  icrc.o	    mprove.o  ranpt.o    sor.o\
addint.o  chixy.o    expint.o  icrc1.o	  mpsqrt.o  ratint.o   sort.o\
airy.o    choldc.o   f1dim.o   igray.o	  mrqcof.o  ratlsq.o   sort2.o\
amebsa.o  cholsl.o   factln.o  iindexx.o  mrqmin.o  ratval.o   sort3.o\
amoeba.o  chsone.o   factrl.o  indexx.o   newt.o    rc.o       spctrm.o\
amotry.o  chstwo.o   fasper.o  interp.o   nrutil.o  rd.o       spear.o\
amotsa.o  cisi.o     fdjac.o   irbit1.o   odeint.o  realft.o   sphbes.o\
anneal.o  cntab1.o   fgauss.o  irbit2.o   orthog.o  rebin.o    splie2.o\
anorm2.o  cntab2.o   fill0.o   jacobi.o   pade.o    red.o      splin2.o\
arcmak.o  complex.o  fit.o     jacobn.o   pccheb.o  relax.o    spline.o\
arcode.o  convlv.o   fitexy.o  julday.o   pcshft.o  relax2.o   splint.o\
arcsum.o  copy.o     fixrts.o  kendl1.o   pearsn.o  resid.o    spread.o\
asolve.o  correl.o   fleg.o    kendl2.o   period.o  revcst.o   sprsax.o\
atimes.o  cosft1.o   flmoon.o  kermom.o   piksr2.o  reverse.o  sprsin.o\
avevar.o  cosft2.o   fmin.o    ks2d1s.o   piksrt.o  rf.o       sprspm.o\
balanc.o  covsrt.o   four1.o   ks2d2s.o   pinvs.o   rj.o       sprstm.o\
banbks.o  crank.o    fourew.o  ksone.o	  plgndr.o  rk4.o      sprstp.o\
bandec.o  cyclic.o   fourfs.o  kstwo.o	  poidev.o  rkck.o     sprstx.o\
banmul.o  daub4.o    fourn.o   laguer.o   polcoe.o  rkdumb.o   stifbs.o\
bcucof.o  dawson.o   fpoly.o   lfit.o	    polcof.o  rkqs.o     stiff.o\
bcuint.o  dbrent.o   fred2.o   linbcg.o   poldiv.o  rlft3.o    stoerm.o\
beschb.o  ddpoly.o   fredin.o  linmin.o   polin2.o  rofunc.o   svbksb.o\
bessi.o   decchk.o   frenel.o  lnsrch.o   polint.o  rotate.o   svdcmp.o\
bessi0.o  df1dim.o   frprmn.o  locate.o   powell.o  rsolv.o    svdfit.o\
bessi1.o  dfour1.o   ftest.o   lop.o	    predic.o  rstrct.o   svdvar.o\
bessik.o  dfpmin.o   gamdev.o  lubksb.o   probks.o  rtbis.o    toeplz.o\
bessj.o   dfridr.o   gammln.o  ludcmp.o   psdes.o   rtflsp.o   tptest.o\
bessj0.o  dftcor.o   gammp.o   machar.o   pwt.o     rtnewt.o   tqli.o\
bessj1.o  dftint.o   gammq.o   matadd.o   pwtset.o  rtsafe.o   trapzd.o\
bessjy.o  difeq.o    gasdev.o  matsub.o   pythag.o  rtsec.o    tred2.o\
bessk.o   dlinmin.o  gaucof.o  medfit.o   pzextr.o  rzextr.o   tridag.o\
bessk0.o  dpythag.o  gauher.o  memcof.o   qgaus.o   savgol.o   trncst.o\
bessk1.o  drealft.o  gaujac.o  metrop.o   qrdcmp.o  scrsho.o   trnspt.o\
bessy.o   dsprsax.o  gaulag.o  mgfas.o	  qromb.o   select.o   ttest.o\
bessy0.o  dsprstx.o  gauleg.o  mglin.o	  qromo.o   selip.o    tutest.o\
bessy1.o  dsvbksb.o  gaussj.o  midexp.o   qroot.o   shell.o    twofft.o\
beta.o    dsvdcmp.o  gcf.o     midinf.o   qrsolv.o  shoot.o    vander.o\
betacf.o  eclass.o   golden.o  midpnt.o   qrupdt.o  vegas.o\
betai.o   eclazz.o   gser.o    midsql.o   qsimp.o   simp1.o    voltra.o\
bico.o    ei.o	     hpsel.o   midsqu.o   qtrap.o   simp2.o    wt1.o\
bksub.o   eigsrt.o   hpsort.o  miser.o	  quad3d.o  simp3.o    wtn.o\
bnldev.o  elle.o     hqr.o     mmid.o	  quadct.o  simplx.o   wwghts.o\
brent.o   ellf.o     hufapp.o  mnewt.o	  quadmx.o  simpr.o    zbrac.o\
broydn.o  ellpi.o    hufdec.o  moment.o   quadvl.o  sinft.o    zbrak.o\
bsstep.o  elmhes.o   hufenc.o  mp2dfr.o   ran0.o    slvsm2.o   zbrent.o\
caldat.o  erfcc.o    hufmak.o  mpdiv.o	  ran1.o    slvsml.o   zrhqr.o\
chder.o   erff.o     hunt.o    mpinv.o	  ran2.o    sncndn.o   zriddr.o\
chebev.o  erffc.o    hypdrv.o  mpmul.o	  ran3.o    snrm.o     zroots.o\
chebft.o  eulsum.o   hypgeo.o  mpops.o	  ran4.o    sobseq.o\
chebpc.o  evlmem.o   hypser.o  mppi.o	  rank.o    solvde.o\
myodeint.o myrkqs.o  dodeint.o drkck.o    drkqs.o
# MNBRAK.obj  chint.obj    expdev.obj  icrc.obj	  mprove.obj  ranpt.obj    sor.obj\
# addint.obj  chixy.obj    expint.obj  icrc1.obj	  mpsqrt.obj  ratint.obj   sort.obj\
# airy.obj    choldc.obj   f1dim.obj   igray.obj	  mrqcof.obj  ratlsq.obj   sort2.obj\
# amebsa.obj  cholsl.obj   factln.obj  iindexx.obj  mrqmin.obj  ratval.obj   sort3.obj\
# amoeba.obj  chsone.obj   factrl.obj  indexx.obj   newt.obj    rc.obj       spctrm.obj\
# amotry.obj  chstwo.obj   fasper.obj  interp.obj   nrutil.obj  rd.obj       spear.obj\
# amotsa.obj  cisi.obj     fdjac.obj   irbit1.obj   odeint.obj  realft.obj   sphbes.obj\
# anneal.obj  cntab1.obj   fgauss.obj  irbit2.obj   orthog.obj  rebin.obj    splie2.obj\
# anorm2.obj  cntab2.obj   fill0.obj   jacobi.obj   pade.obj    red.obj      splin2.obj\
# arcmak.obj  complex.obj  fit.obj     jacobn.obj   pccheb.obj  relax.obj    spline.obj\
# arcode.obj  convlv.obj   fitexy.obj  julday.obj   pcshft.obj  relax2.obj   splint.obj\
# arcsum.obj  copy.obj     fixrts.obj  kendl1.obj   pearsn.obj  resid.obj    spread.obj\
# asolve.obj  correl.obj   fleg.obj    kendl2.obj   period.obj  revcst.obj   sprsax.obj\
# atimes.obj  cosft1.obj   flmoon.obj  kermom.obj   piksr2.obj  reverse.obj  sprsin.obj\
# avevar.obj  cosft2.obj   fmin.obj    ks2d1s.obj   piksrt.obj  rf.obj       sprspm.obj\
# balanc.obj  covsrt.obj   four1.obj   ks2d2s.obj   pinvs.obj   rj.obj       sprstm.obj\
# banbks.obj  crank.obj    fourew.obj  ksone.obj	  plgndr.obj  rk4.obj      sprstp.obj\
# bandec.obj  cyclic.obj   fourfs.obj  kstwo.obj	  poidev.obj  rkck.obj     sprstx.obj\
# banmul.obj  daub4.obj    fourn.obj   laguer.obj   polcoe.obj  rkdumb.obj   stifbs.obj\
# bcucof.obj  dawson.obj   fpoly.obj   lfit.obj	  polcof.obj  rkqs.obj     stiff.obj\
# bcuint.obj  dbrent.obj   fred2.obj   linbcg.obj   poldiv.obj  rlft3.obj    stoerm.obj\
# beschb.obj  ddpoly.obj   fredin.obj  linmin.obj   polin2.obj  rofunc.obj   svbksb.obj\
# bessi.obj   decchk.obj   frenel.obj  lnsrch.obj   polint.obj  rotate.obj   svdcmp.obj\
# bessi0.obj  df1dim.obj   frprmn.obj  locate.obj   powell.obj  rsolv.obj    svdfit.obj\
# bessi1.obj  dfour1.obj   ftest.obj   lop.obj	  predic.obj  rstrct.obj   svdvar.obj\
# bessik.obj  dfpmin.obj   gamdev.obj  lubksb.obj   probks.obj  rtbis.obj    toeplz.obj\
# bessj.obj   dfridr.obj   gammln.obj  ludcmp.obj   psdes.obj   rtflsp.obj   tptest.obj\
# bessj0.obj  dftcor.obj   gammp.obj   machar.obj   pwt.obj     rtnewt.obj   tqli.obj\
# bessj1.obj  dftint.obj   gammq.obj   matadd.obj   pwtset.obj  rtsafe.obj   trapzd.obj\
# bessjy.obj  difeq.obj    gasdev.obj  matsub.obj   pythag.obj  rtsec.obj    tred2.obj\
# bessk.obj   dlinmin.obj  gaucof.obj  medfit.obj   pzextr.obj  rzextr.obj   tridag.obj\
# bessk0.obj  dpythag.obj  gauher.obj  memcof.obj   qgaus.obj   savgol.obj   trncst.obj\
# bessk1.obj  drealft.obj  gaujac.obj  metrop.obj   qrdcmp.obj  scrsho.obj   trnspt.obj\
# bessy.obj   dsprsax.obj  gaulag.obj  mgfas.obj	  qromb.obj   select.obj   ttest.obj\
# bessy0.obj  dsprstx.obj  gauleg.obj  mglin.obj	  qromo.obj   selip.obj    tutest.obj\
# bessy1.obj  dsvbksb.obj  gaussj.obj  midexp.obj   qroot.obj   shell.obj    twofft.obj\
# beta.obj    dsvdcmp.obj  gcf.obj     midinf.obj   qrsolv.obj  shoot.obj    vander.obj\
# betacf.obj  eclass.obj   golden.obj  midpnt.obj   qrupdt.obj  vegas.obj\
# betai.obj   eclazz.obj   gser.obj    midsql.obj   qsimp.obj   simp1.obj    voltra.obj\
# bico.obj    ei.obj	 hpsel.obj   midsqu.obj   qtrap.obj   simp2.obj    wt1.obj\
# bksub.obj   eigsrt.obj   hpsort.obj  miser.obj	  quad3d.obj  simp3.obj    wtn.obj\
# bnldev.obj  elle.obj     hqr.obj     mmid.obj	  quadct.obj  simplx.obj   wwghts.obj\
# brent.obj   ellf.obj     hufapp.obj  mnewt.obj	  quadmx.obj  simpr.obj    zbrac.obj\
# broydn.obj  ellpi.obj    hufdec.obj  moment.obj   quadvl.obj  sinft.obj    zbrak.obj\
# bsstep.obj  elmhes.obj   hufenc.obj  mp2dfr.obj   ran0.obj    slvsm2.obj   zbrent.obj\
# caldat.obj  erfcc.obj    hufmak.obj  mpdiv.obj	  ran1.obj    slvsml.obj   zrhqr.obj\
# chder.obj   erff.obj     hunt.obj    mpinv.obj	  ran2.obj    sncndn.obj   zriddr.obj\
# chebev.obj  erffc.obj    hypdrv.obj  mpmul.obj	  ran3.obj    snrm.obj     zroots.obj\
# chebft.obj  eulsum.obj   hypgeo.obj  mpops.obj	  ran4.obj    sobseq.obj\
# chebpc.obj  evlmem.obj   hypser.obj  mppi.obj	  rank.obj    solvde.obj\
# myodeint.obj myrkqs.obj  dodeint.obj drkck.obj    drkqs.obj
# PB 03/Oct/2003   Removed because of compilation problems:
#shootf.obj 

#CC=C:\MinGW\msys\1.0\bin\gcc
#CXX=C:\MinGW\msys\1.0\bin\g++
CC=g++ #C:\"Program Files"\MATLAB\R2018b\sys\lcc64\lcc64\bin\lcc64
#CC="C:\Program Files\CodeBlocks\MinGW\bin\g++"


CFLAGS= -I$(INCLUDEPATH) -Iquote -O2 -fPIC #-march= $(ARCH)
CXXFLAGS= -I$(INCLUDEPATH) -Iquote -O2 #-march= $(ARCH)
DEBUGFLAGS= -I$(INCLUDEPATH) -Iquote -g

.c.obj: 	
	$(CC) $(CFLAGS) -c $@.c
#	$(CC) $(CFLAGS) -c -Fo$@ $<

#%.obj: %.cpp
#	$(CC) $(CFLAGS) -c $*.cpp
#	$(CC)  -c $*.cpp
#	$(CC) $(CFLAGS) -c -Fo$@ $<

#.cpp.obj: 
#	$(CXX) $(CXXFLAGS) -c $<
# #	$(CC) $(CFLAGS) -c -Fo$@ $<

.cpp.o: 
	$(CC) $(CFLAGS) -c -Fo$@ $<


$(LIBTGT): $(LIBOBJS)
	ar rcs $@ $^
	mv $@ $(LIBPATH)

#	lib /out:$(LIBTGT) $(LIBOBJS)
#	lcclib64 /out:$(LIBTGT) $(LIBOBJS)
#	copy $(LIBTGT) $(LIBPATH)


