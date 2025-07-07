LIBTARGET = nrd
LIBEXT = lib #lib #a 
LIBTGT = $(LIBTARGET).$(LIBEXT)
HOME=/Users/tsikyrakotomalala/Documents/these/Traj_Optimization_lstm/Traj_Optimization_lstm/
# HOME=/users/pierre/Documents
#HOME=d:
#HOME=z:\pierre\Documents

#ifeq ($(LIBEXT),a64) 
ARCH=x86_64
#else
#ARCH=i386
#endif

INCLUDEPATH=$(HOME)\Includes\nrd #"C:\Program Files\MATLAB\R2018b\sys\lcc64\lcc64\include64"
LIBPATH=$(HOME)\Libs
LIBOBJS=\
MNBRAK.o  chint.o    expdev.o  icrc.o	  mprove.o  ranpt.o    sor.o\
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
bcucof.o  dawson.o   fpoly.o   lfit.o	  polcof.o  rkqs.o     stiff.o\
bcuint.o  dbrent.o   fred2.o   linbcg.o   poldiv.o  rlft3.o    stoerm.o\
beschb.o  ddpoly.o   fredin.o  linmin.o   polin2.o  rofunc.o   svbksb.o\
bessi.o   decchk.o   frenel.o  lnsrch.o   polint.o  rotate.o   svdcmp.o\
bessi0.o  df1dim.o   frprmn.o  locate.o   powell.o  rsolv.o    svdfit.o\
bessi1.o  dfour1.o   ftest.o   lop.o	  predic.o  rstrct.o   svdvar.o\
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
beta.o	  dsvdcmp.o  gcf.o     midinf.o   qrsolv.o  shoot.o    vander.o\
betacf.o  eclass.o   golden.o  midpnt.o   qrupdt.o             vegas.o\
betai.o   eclazz.o   gser.o    midsql.o   qsimp.o   simp1.o    voltra.o\
bico.o	  ei.o	     hpsel.o   midsqu.o   qtrap.o   simp2.o    wt1.o\
bksub.o   eigsrt.o   hpsort.o  miser.o	  quad3d.o  simp3.o    wtn.o\
bnldev.o  elle.o     hqr.o     mmid.o	  quadct.o  simplx.o   wwghts.o\
brent.o   ellf.o     hufapp.o  mnewt.o	  quadmx.o  simpr.o    zbrac.o\
broydn.o  ellpi.o    hufdec.o  moment.o   quadvl.o  sinft.o    zbrak.o\
bsstep.o  elmhes.o   hufenc.o  mp2dfr.o   ran0.o    slvsm2.o   zbrent.o\
caldat.o  erfcc.o    hufmak.o  mpdiv.o	  ran1.o    slvsml.o   zrhqr.o\
chder.o   erff.o     hunt.o    mpinv.o	  ran2.o    sncndn.o   zriddr.o\
chebev.o  erffc.o    hypdrv.o  mpmul.o	  ran3.o    snrm.o     zroots.o\
chebft.o  eulsum.o   hypgeo.o  mpops.o	  ran4.o    sobseq.o\
chebpc.o  evlmem.o   hypser.o  mppi.o	  rank.o   solvde.o\
myodeint.o myrkqs.o  dodeint.o drkck.o    drkqs.o
# PB 03/Oct/2003   Removed because of compilation problems:
#shootf.obj 

CC=/usr/bin/gcc
#CC="c:\Program Files\MATLAB\R2018b\sys\lcc64\lcc64\bin\lcc64.exe"
#CC=C:\"Program Files"\CodeBlocks\MinGW\bin\gcc
#CC=cl


CFLAGS= -I$(INCLUDEPATH) -Iquote -O2 #-march= $(ARCH)
DEBUGFLAGS= -I$(INCLUDEPATH) -Iquote -g

.c.o: 	
	$(CC) $(CFLAGS) -c $@.c
#	$(CC) $(CFLAGS) -c -Fo$@ $<

%.obj: %.cpp
	$(CC) $(CFLAGS) -c $*.cpp
#	$(CC)  -c $*.cpp
#	$(CC) $(CFLAGS) -c -Fo$@ $<

# .cpp.obj: 
	# $(CC) $(CFLAGS) -c $<
# #	$(CC) $(CFLAGS) -c -Fo$@ $<



$(LIBTGT): $(LIBOBJS)
	ar ru $@ $^
	mv $@ $(LIBPATH)

	#nrd.a /out:$(LIBTGT) $(LIBOBJS)
#	lcclib64 /out:$(LIBTGT) $(LIBOBJS)
	#copy $(LIBTGT) $(LIBPATH)


