# Makefile for Recipes: BSD UNIX distribution
#

CC = $(CCANSI)
CFLAGS = $(CFLAGSANSI)

OBJS =	addint.o airy.o amebsa.o amoeba.o amotry.o amotsa.o anneal.o \
	anorm2.o arcmak.o arcode.o arcsum.o asolve.o atimes.o avevar.o \
	balanc.o banbks.o bandec.o banmul.o bcucof.o bcuint.o beschb.o \
	bessi.o bessi0.o bessi1.o bessik.o bessj.o bessj0.o bessj1.o \
	bessjy.o bessk.o bessk0.o bessk1.o bessy.o bessy0.o bessy1.o \
	beta.o betacf.o betai.o bico.o bksub.o bnldev.o brent.o broydn.o \
	bsstep.o caldat.o chder.o chebev.o chebft.o chebpc.o chint.o \
	chixy.o choldc.o cholsl.o chsone.o chstwo.o cisi.o complex.o cntab1.o \
	cntab2.o convlv.o copy.o correl.o cosft1.o cosft2.o covsrt.o \
	crank.o cyclic.o daub4.o dawson.o dbrent.o ddpoly.o decchk.o \
	df1dim.o dfour1.o dfpmin.o dfridr.o dftcor.o dftint.o difeq.o \
	dlinmin.o dpythag.o drealft.o dsprsax.o dsprstx.o dsvbksb.o \
	dsvdcmp.o eclass.o eclazz.o ei.o eigsrt.o elle.o ellf.o ellpi.o \
	elmhes.o erfcc.o erff.o erffc.o eulsum.o evlmem.o expdev.o \
	expint.o f1dim.o factln.o factrl.o fasper.o fdjac.o fgauss.o \
	fill0.o fit.o fitexy.o fixrts.o fleg.o flmoon.o fmin_nr.o four1.o \
	fourew.o fourfs.o fourn.o fpoly.o fred2.o fredin.o frenel.o \
	frprmn.o ftest.o gamdev.o gammln.o gammp.o gammq.o gasdev.o \
	gaucof.o gauher.o gaujac.o gaulag.o gauleg.o gaussj.o gcf.o \
	golden.o gser.o hpsel.o hpsort.o hqr.o hufapp.o hufdec.o \
	hufenc.o hufmak.o hunt.o hypdrv.o hypgeo.o hypser.o icrc.o \
	icrc1.o igray.o iindexx.o indexx.o interp.o irbit1.o irbit2.o \
	jacobi.o jacobn.o julday.o kendl1.o kendl2.o kermom.o ks2d1s.o \
	ks2d2s.o ksone.o kstwo.o laguer.o lfit.o linbcg.o linmin.o \
	lnsrch.o locate.o lop.o lubksb.o ludcmp.o machar.o matadd.o \
	matsub.o medfit.o memcof.o metrop.o mgfas.o mglin.o midexp.o \
	midinf.o midpnt.o midsql.o midsqu.o miser.o mmid.o mnbrak.o \
	mnewt.o moment.o mp2dfr.o mpdiv.o mpinv.o mpmul.o mpops.o mppi.o \
	mprove.o mpsqrt.o mrqcof.o mrqmin.o newt.o nrutil.o odeint.o orthog.o \
	pade.o pccheb.o pcshft.o pearsn.o period.o piksr2.o piksrt.o \
	pinvs.o plgndr.o poidev.o polcoe.o polcof.o poldiv.o polin2.o \
	polint.o powell.o predic.o probks.o psdes.o pwt.o pwtset.o \
	pythag.o pzextr.o qgaus.o qrdcmp.o qromb.o qromo.o qroot.o \
	qrsolv.o qrupdt.o qsimp.o qtrap.o quad3d.o quadct.o quadmx.o \
	quadvl.o ran0.o ran1.o ran2.o ran3.o ran4.o rank.o ranpt.o \
	ratint.o ratlsq.o ratval.o rc.o rd.o realft.o rebin.o red.o \
	relax.o relax2.o resid.o revcst.o reverse.o rf.o rj.o rk4.o \
	rkck.o rkdumb.o rkqs.o rlft3.o rofunc.o rotate.o rsolv.o \
	rstrct.o rtbis.o rtflsp.o rtnewt.o rtsafe.o rtsec.o rzextr.o \
	savgol.o scrsho.o select.o selip.o shell.o shoot.o shootf.o \
	simp1.o simp2.o simp3.o simplx.o simpr.o sinft.o slvsm2.o \
	slvsml.o sncndn.o snrm.o sobseq.o solvde.o sor.o sort.o sort2.o \
	sort3.o spctrm.o spear.o sphbes.o splie2.o splin2.o spline.o \
	splint.o spread.o sprsax.o sprsin.o sprspm.o sprstm.o sprstp.o \
	sprstx.o stifbs.o stiff.o stoerm.o svbksb.o svdcmp.o svdfit.o \
	svdvar.o toeplz.o tptest.o tqli.o trapzd.o tred2.o tridag.o \
	trncst.o trnspt.o ttest.o tutest.o twofft.o vander.o vegas.o \
	voltra.o wt1.o wtn.o wwghts.o zbrac.o zbrak.o zbrent.o zrhqr.o \
	zriddr.o zroots.o

#If any of the above have problems with an optimization option or other
#compiler flags, a specific rule will be needed for them. This is very
#simple to write. For example, for svdcmp (which sometimes has problems
#with optimizers), a very plain rule is:
#svdcmp.o: svdcmp.c ; $(CC) -I../include -c $<
#Uncomment or write analogously for other recipes as needed.

all: binrecipes

binrecipes: $(CANSILIB) badluk fredex sfroid sphfpt sphoot

$(CANSILIB): $(OBJS)
	$(AR) $(ARFLAGS) $(CANSILIB) $(OBJS)
	$(RANLIB) $(CANSILIB)

install: binrecipes
	cp $(CANSILIB) $(LIBDIR)
	$(RANLIB) $(LIBDIR)/$(CANSILIB)
	cp badluk fredex sfroid sphfpt sphoot $(CBIN)

# Not cp -p, ranlib -t - can't assume -p option for all systems

badluk: badluk.c $(CANSILIB)
	$(CC) $(CFLAGS) -o badluk badluk.c $(CANSILIB) -lm

fredex: fredex.c $(CANSILIB)
	$(CC) $(CFLAGS) -o fredex fredex.c $(CANSILIB) -lm

sfroid: sfroid.c $(CANSILIB)
	$(CC) $(CFLAGS) -o sfroid sfroid.c $(CANSILIB) -lm

sphfpt: sphfpt.c $(CANSILIB) ../demo/src/xsphfpt.c
	$(CC) $(CFLAGS) -o sphfpt sphfpt.c ../demo/src/xsphfpt.c \
	$(CANSILIB) -lm

sphoot: sphoot.c $(CANSILIB)
	$(CC) $(CFLAGS) -o sphoot sphoot.c $(CANSILIB) -lm

clean:	
	/bin/rm -f *.o *.a badluk fredex sfroid sphfpt sphoot
# DO NOT DELETE THIS LINE -- make depend depends on it.


amebsa.o:  ../include/nrutil.h
amoeba.o:  ../include/nrutil.h
amotry.o: ../include/nrutil.h
amotsa.o:  ../include/nrutil.h
arcmak.o: ../include/nrutil.h
banmul.o: ../include/nrutil.h
bcuint.o: ../include/nrutil.h
bessjy.o:  ../include/nrutil.h
brent.o:  ../include/nrutil.h
broydn.o:  ../include/nrutil.h
bsstep.o:  ../include/nrutil.h
chebft.o:  ../include/nrutil.h
chebpc.o: ../include/nrutil.h
chixy.o:  ../include/nrutil.h
cisi.o:  ../include/complex.h
cntab1.o:  ../include/nrutil.h
cntab2.o:  ../include/nrutil.h
convlv.o: ../include/nrutil.h
correl.o: ../include/nrutil.h
cyclic.o: ../include/nrutil.h
daub4.o: ../include/nrutil.h
dawson.o:  ../include/nrutil.h
dbrent.o:  ../include/nrutil.h
df1dim.o: ../include/nrutil.h
dfpmin.o:  ../include/nrutil.h
dfridr.o:  ../include/nrutil.h
dftint.o:  ../include/nrutil.h
dlinmin.o: ../include/nrutil.h
dpythag.o:  ../include/nrutil.h
dsvbksb.o: ../include/nrutil.h
dsvdcmp.o:  ../include/nrutil.h
elle.o:  ../include/nrutil.h
ellf.o:  ../include/nrutil.h
ellpi.o:  ../include/nrutil.h
f1dim.o: ../include/nrutil.h
fasper.o:  ../include/nrutil.h
fdjac.o:  ../include/nrutil.h
fit.o:  ../include/nrutil.h
fitexy.o:  ../include/nrutil.h
fixrts.o:  ../include/complex.h
fmin_nr.o: ../include/nrutil.h
fourfs.o:   ../include/nrutil.h
fred2.o: ../include/nrutil.h
fredex.o:   ../include/nrutil.h
frenel.o:  ../include/complex.h
frprmn.o:  ../include/nrutil.h
gaucof.o:  ../include/nrutil.h
gaussj.o:  ../include/nrutil.h
hqr.o:  ../include/nrutil.h
hufmak.o: ../include/nrutil.h
hypdrv.o: ../include/complex.h
hypgeo.o:  ../include/complex.h ../include/nrutil.h
hypser.o: ../include/complex.h
iindexx.o: ../include/nrutil.h
indexx.o: ../include/nrutil.h
jacobi.o:  ../include/nrutil.h
ks2d1s.o:  ../include/nrutil.h
ks2d2s.o:  ../include/nrutil.h
ksone.o:  ../include/nrutil.h
laguer.o:  ../include/complex.h ../include/nrutil.h
lfit.o: ../include/nrutil.h
linbcg.o:   ../include/nrutil.h
linmin.o: ../include/nrutil.h
lnsrch.o:  ../include/nrutil.h
ludcmp.o:  ../include/nrutil.h
medfit.o:  ../include/nrutil.h
memcof.o:  ../include/nrutil.h
metrop.o: 
mgfas.o: ../include/nrutil.h
mglin.o: ../include/nrutil.h
miser.o: ../include/nrutil.h
mmid.o: ../include/nrutil.h
mnbrak.o:  ../include/nrutil.h
mnewt.o:  ../include/nrutil.h
mpdiv.o: ../include/nrutil.h
mpinv.o: ../include/nrutil.h
mpmul.o: ../include/nrutil.h
mppi.o:  ../include/nrutil.h
mprove.o: ../include/nrutil.h
mpsqrt.o:  ../include/nrutil.h
mrqcof.o: ../include/nrutil.h
mrqmin.o: ../include/nrutil.h
newt.o:  ../include/nrutil.h
odeint.o:  ../include/nrutil.h
orthog.o: ../include/nrutil.h
pade.o:  ../include/nrutil.h
period.o:  ../include/nrutil.h
pinvs.o:  ../include/nrutil.h
polcoe.o: ../include/nrutil.h
polcof.o:  ../include/nrutil.h
polin2.o: ../include/nrutil.h
polint.o:  ../include/nrutil.h
powell.o:  ../include/nrutil.h
predic.o: ../include/nrutil.h
pwt.o: ../include/nrutil.h
pythag.o:  ../include/nrutil.h
pzextr.o: ../include/nrutil.h
qrdcmp.o:  ../include/nrutil.h
qroot.o:  ../include/nrutil.h
qrupdt.o:  ../include/nrutil.h
quadmx.o:  ../include/nrutil.h
quadvl.o: ../include/nrutil.h
ratint.o:  ../include/nrutil.h
ratlsq.o:   ../include/nrutil.h
rc.o:  ../include/nrutil.h
rd.o:  ../include/nrutil.h
rf.o:  ../include/nrutil.h
rj.o:  ../include/nrutil.h
rk4.o: ../include/nrutil.h
rkck.o: ../include/nrutil.h
rkdumb.o: ../include/nrutil.h
rkqs.o:  ../include/nrutil.h
rlftfrag.o: ../include/nrutil.h
rofunc.o:  ../include/nrutil.h
rotate.o:  ../include/nrutil.h
rzextr.o: ../include/nrutil.h
savgol.o:  ../include/nrutil.h
selip.o: ../include/nrutil.h
sfroid.o:   ../include/nrutil.h
shoot.o: ../include/nrutil.h
shootf.o: ../include/nrutil.h
simplx.o: ../include/nrutil.h
simpr.o: ../include/nrutil.h
sobseq.o: ../include/nrutil.h
solvde.o:   ../include/nrutil.h
sort.o: ../include/nrutil.h
sort2.o: ../include/nrutil.h
sort3.o: ../include/nrutil.h
spctrm.o:   ../include/nrutil.h
spear.o:  ../include/nrutil.h
sphfpt.o:   ../include/nrutil.h
sphoot.o:  ../include/nrutil.h
splin2.o: ../include/nrutil.h
spline.o: ../include/nrutil.h
spread.o: ../include/nrutil.h
stifbs.o:  ../include/nrutil.h
stiff.o:  ../include/nrutil.h
stoerm.o: ../include/nrutil.h
svbksb.o: ../include/nrutil.h
svdcmp.o:  ../include/nrutil.h
svdfit.o: ../include/nrutil.h
svdvar.o: ../include/nrutil.h
toeplz.o: ../include/nrutil.h
tqli.o:  ../include/nrutil.h
tridag.o: ../include/nrutil.h
trnspt.o: ../include/nrutil.h
tutest.o:  ../include/nrutil.h
vander.o: ../include/nrutil.h
vegas.o:   ../include/nrutil.h
voltra.o: ../include/nrutil.h
wtn.o: ../include/nrutil.h
zbrent.o:  ../include/nrutil.h
zrhqr.o: ../include/nrutil.h
zriddr.o:  ../include/nrutil.h
zroots.o:  ../include/complex.h
