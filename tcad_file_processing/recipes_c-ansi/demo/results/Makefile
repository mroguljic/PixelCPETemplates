# Makefile for Numerical Recipes Demo binaries

CC = $(CCANSI)
CFLAGS = $(CFLAGSANSI)
ILIB = $(LIBDIR)/$(CANSILIB)
BINDIR = ../$(CBIN)

all :	xairy.reslt xamebsa.reslt xamoeba.reslt xanneal.reslt \
	xarcode.reslt xavevar.reslt xbalanc.reslt xbandec.reslt \
	xbanmul.reslt xbcucof.reslt xbcuint.reslt xbeschb.reslt \
	xbessi.reslt xbessi0.reslt xbessi1.reslt xbessik.reslt \
	xbessj.reslt xbessj0.reslt xbessj1.reslt xbessjy.reslt \
	xbessk.reslt xbessk0.reslt xbessk1.reslt xbessy.reslt \
	xbessy0.reslt xbessy1.reslt xbeta.reslt xbetai.reslt xbico.reslt \
	xbnldev.reslt xbrent.reslt xbroydn.reslt xbsstep.reslt \
	xcaldat.reslt xchder.reslt xchebev.reslt xchebft.reslt \
	xchebpc.reslt xchint.reslt xcholsl.reslt xchsone.reslt \
	xchstwo.reslt xcisi.reslt xcntab1.reslt xcntab2.reslt \
	xconvlv.reslt xcorrel.reslt xcosft1.reslt xcosft2.reslt \
	xcovsrt.reslt xcrank.reslt xcyclic.reslt xdawson.reslt \
	xdbrent.reslt xddpoly.reslt xdecchk.reslt xdf1dim.reslt \
	xdfpmin.reslt xdfridr.reslt xdftint.reslt xdlinmin.reslt \
	xeclass.reslt xeclazz.reslt xei.reslt xeigsrt.reslt xelle.reslt \
	xellf.reslt xellpi.reslt xelmhes.reslt xerfcc.reslt xerff.reslt \
	xerffc.reslt xeulsum.reslt xevlmem.reslt xexpdev.reslt \
	xexpint.reslt xf1dim.reslt xfactln.reslt xfactrl.reslt \
	xfasper.reslt xfgauss.reslt xfit.reslt xfitexy.reslt \
	xfixrts.reslt xfleg.reslt xflmoon.reslt xfour1.reslt \
	xfourfs.reslt xfourn.reslt xfpoly.reslt xfred2.reslt \
	xfredin.reslt xfrenel.reslt xfrprmn.reslt xftest.reslt \
	xgamdev.reslt xgammln.reslt xgammp.reslt xgammq.reslt \
	xgasdev.reslt xgaucof.reslt xgauher.reslt xgaujac.reslt \
	xgaulag.reslt xgauleg.reslt xgaussj.reslt xgcf.reslt \
	xgolden.reslt xgser.reslt xhpsel.reslt xhpsort.reslt xhqr.reslt \
	xhuffman.reslt xhunt.reslt xhypgeo.reslt xicrc.reslt \
	xigray.reslt xindexx.reslt xirbit1.reslt xirbit2.reslt \
	xjacobi.reslt xjulday.reslt xkendl1.reslt xkendl2.reslt \
	xks2d1s.reslt xks2d2s.reslt xksone.reslt xkstwo.reslt \
	xlaguer.reslt xlfit.reslt xlinbcg.reslt xlinmin.reslt \
	xlocate.reslt xlubksb.reslt xludcmp.reslt xmachar.reslt \
	xmedfit.reslt xmemcof.reslt xmgfas.reslt xmglin.reslt \
	xmidpnt.reslt xmiser.reslt xmmid.reslt xmnbrak.reslt \
	xmnewt.reslt xmoment.reslt xmppi.reslt xmprove.reslt \
	xmrqcof.reslt xmrqmin.reslt xnewt.reslt xodeint.reslt \
	xorthog.reslt xpade.reslt xpccheb.reslt xpcshft.reslt \
	xpearsn.reslt xperiod.reslt xpiksr2.reslt xpiksrt.reslt \
	xplgndr.reslt xpoidev.reslt xpolcoe.reslt xpolcof.reslt \
	xpoldiv.reslt xpolin2.reslt xpolint.reslt xpowell.reslt \
	xpredic.reslt xprobks.reslt xpsdes.reslt xpzextr.reslt \
	xqgaus.reslt xqrdcmp.reslt xqromb.reslt xqromo.reslt \
	xqroot.reslt xqrsolv.reslt xqrupdt.reslt xqsimp.reslt \
	xqtrap.reslt xquad3d.reslt xran.reslt xran4.reslt xrank.reslt \
	xratint.reslt xratlsq.reslt xrc.reslt xrd.reslt xrealft.reslt \
	xrf.reslt xrj.reslt xrk4.reslt xrkdumb.reslt xrkqs.reslt \
	xrlft3.reslt xrofunc.reslt xrtbis.reslt xrtflsp.reslt \
	xrtnewt.reslt xrtsafe.reslt xrtsec.reslt xrzextr.reslt \
	xsavgol.reslt xscrsho.reslt xselect.reslt xselip.reslt \
	xshell.reslt xsimplx.reslt xsimpr.reslt xsinft.reslt \
	xsncndn.reslt xsobseq.reslt xsor.reslt xsort.reslt xsort2.reslt \
	xsort3.reslt xspctrm.reslt xspear.reslt xsphbes.reslt \
	xsphfpt.reslt xsplie2.reslt xsplin2.reslt xspline.reslt \
	xsplint.reslt xsprsax.reslt xsprsin.reslt xsprspm.reslt \
	xsprstm.reslt xsprstp.reslt xsprstx.reslt xstifbs.reslt \
	xstiff.reslt xstoerm.reslt xsvbksb.reslt xsvdcmp.reslt \
	xsvdfit.reslt xsvdvar.reslt xtoeplz.reslt xtptest.reslt \
	xtqli.reslt xtrapzd.reslt xtred2.reslt xtridag.reslt \
	xttest.reslt xtutest.reslt xtwofft.reslt xvander.reslt \
	xvegas.reslt xvoltra.reslt xwt1.reslt xwtn.reslt xzbrac.reslt \
	xzbrak.reslt xzbrent.reslt xzrhqr.reslt xzriddr.reslt \
	xzroots.reslt badluk.reslt fredex.reslt sfroid.reslt sphoot.reslt

# The following is extremely inelegant and ugly, but there seems to be no
# better way to accomplish the objectives in a manner portable across
# a wide variety of systems

xairy.reslt : ../src/xairy.c $(ILIB)
	-if [ -r ../bin/xairy.out ]; then \
	../bin/xairy.out > xairy.reslt; \
	else \
	$(CC) $(CFLAGS) -o xairy.out ../src/xairy.c $(ILIB) -lm && \
	xairy.out > xairy.reslt && /bin/rm -f xairy.o xairy.out; \
	fi

xamebsa.reslt : ../src/xamebsa.c $(ILIB)
	-if [ -r ../bin/xamebsa.out ]; then \
	../bin/xamebsa.out < ../responses/xamebsa.resp > xamebsa.reslt; \
	else \
	$(CC) $(CFLAGS) -o xamebsa.out ../src/xamebsa.c $(ILIB) -lm && \
	xamebsa.out < ../responses/xamebsa.resp > xamebsa.reslt && \
	/bin/rm -f xamebsa.o xamebsa.out; \
	fi

xamoeba.reslt : ../src/xamoeba.c $(ILIB)
	-if [ -r ../bin/xamoeba.out ]; then \
	../bin/xamoeba.out > xamoeba.reslt; \
	else \
	$(CC) $(CFLAGS) -o xamoeba.out ../src/xamoeba.c $(ILIB) -lm && \
	xamoeba.out > xamoeba.reslt && /bin/rm -f xamoeba.o xamoeba.out; \
	fi

xanneal.reslt : ../src/xanneal.c $(ILIB)
	-if [ -r ../bin/xanneal.out ]; then \
	../bin/xanneal.out > xanneal.reslt; \
	else \
	$(CC) $(CFLAGS) -o xanneal.out ../src/xanneal.c $(ILIB) -lm && \
	xanneal.out > xanneal.reslt && /bin/rm -f xanneal.o xanneal.out; \
	fi

xarcode.reslt : ../src/xarcode.c $(ILIB)
	-if [ -r ../bin/xarcode.out ]; then \
	../bin/xarcode.out < ../responses/xarcode.resp > xarcode.reslt; \
	else \
	$(CC) $(CFLAGS) -o xarcode.out ../src/xarcode.c $(ILIB) -lm && \
	xarcode.out < ../responses/xarcode.resp > xarcode.reslt && \
	/bin/rm -f xarcode.o xarcode.out; \
	fi

xavevar.reslt : ../src/xavevar.c $(ILIB)
	-if [ -r ../bin/xavevar.out ]; then \
	../bin/xavevar.out > xavevar.reslt; \
	else \
	$(CC) $(CFLAGS) -o xavevar.out ../src/xavevar.c $(ILIB) -lm && \
	xavevar.out > xavevar.reslt && /bin/rm -f xavevar.o xavevar.out; \
	fi

xbalanc.reslt : ../src/xbalanc.c $(ILIB)
	-if [ -r ../bin/xbalanc.out ]; then \
	../bin/xbalanc.out > xbalanc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbalanc.out ../src/xbalanc.c $(ILIB) -lm && \
	xbalanc.out > xbalanc.reslt && /bin/rm -f xbalanc.o xbalanc.out; \
	fi

xbandec.reslt : ../src/xbandec.c $(ILIB)
	-if [ -r ../bin/xbandec.out ]; then \
	../bin/xbandec.out > xbandec.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbandec.out ../src/xbandec.c $(ILIB) -lm && \
	xbandec.out > xbandec.reslt && /bin/rm -f xbandec.o xbandec.out; \
	fi

xbanmul.reslt : ../src/xbanmul.c $(ILIB)
	-if [ -r ../bin/xbanmul.out ]; then \
	../bin/xbanmul.out > xbanmul.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbanmul.out ../src/xbanmul.c $(ILIB) -lm && \
	xbanmul.out > xbanmul.reslt && /bin/rm -f xbanmul.o xbanmul.out; \
	fi

xbcucof.reslt : ../src/xbcucof.c $(ILIB)
	-if [ -r ../bin/xbcucof.out ]; then \
	../bin/xbcucof.out > xbcucof.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbcucof.out ../src/xbcucof.c $(ILIB) -lm && \
	xbcucof.out > xbcucof.reslt && /bin/rm -f xbcucof.o xbcucof.out; \
	fi

xbcuint.reslt : ../src/xbcuint.c $(ILIB)
	-if [ -r ../bin/xbcuint.out ]; then \
	../bin/xbcuint.out > xbcuint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbcuint.out ../src/xbcuint.c $(ILIB) -lm && \
	xbcuint.out > xbcuint.reslt && /bin/rm -f xbcuint.o xbcuint.out; \
	fi

xbeschb.reslt : ../src/xbeschb.c $(ILIB)
	-if [ -r ../bin/xbeschb.out ]; then \
	../bin/xbeschb.out < ../responses/xbeschb.resp > xbeschb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbeschb.out ../src/xbeschb.c $(ILIB) -lm && \
	xbeschb.out < ../responses/xbeschb.resp > xbeschb.reslt && \
	/bin/rm -f xbeschb.o xbeschb.out; \
	fi

xbessi.reslt : ../src/xbessi.c $(ILIB)
	-if [ -r ../bin/xbessi.out ]; then \
	../bin/xbessi.out > xbessi.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessi.out ../src/xbessi.c $(ILIB) -lm && \
	xbessi.out > xbessi.reslt && /bin/rm -f xbessi.o xbessi.out; \
	fi

xbessi0.reslt : ../src/xbessi0.c $(ILIB)
	-if [ -r ../bin/xbessi0.out ]; then \
	../bin/xbessi0.out > xbessi0.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessi0.out ../src/xbessi0.c $(ILIB) -lm && \
	xbessi0.out > xbessi0.reslt && /bin/rm -f xbessi0.o xbessi0.out; \
	fi

xbessi1.reslt : ../src/xbessi1.c $(ILIB)
	-if [ -r ../bin/xbessi1.out ]; then \
	../bin/xbessi1.out > xbessi1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessi1.out ../src/xbessi1.c $(ILIB) -lm && \
	xbessi1.out > xbessi1.reslt && /bin/rm -f xbessi1.o xbessi1.out; \
	fi

xbessik.reslt : ../src/xbessik.c $(ILIB)
	-if [ -r ../bin/xbessik.out ]; then \
	../bin/xbessik.out > xbessik.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessik.out ../src/xbessik.c $(ILIB) -lm && \
	xbessik.out > xbessik.reslt && /bin/rm -f xbessik.o xbessik.out; \
	fi

xbessj.reslt : ../src/xbessj.c $(ILIB)
	-if [ -r ../bin/xbessj.out ]; then \
	../bin/xbessj.out > xbessj.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessj.out ../src/xbessj.c $(ILIB) -lm && \
	xbessj.out > xbessj.reslt && /bin/rm -f xbessj.o xbessj.out; \
	fi

xbessj0.reslt : ../src/xbessj0.c $(ILIB)
	-if [ -r ../bin/xbessj0.out ]; then \
	../bin/xbessj0.out > xbessj0.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessj0.out ../src/xbessj0.c $(ILIB) -lm && \
	xbessj0.out > xbessj0.reslt && /bin/rm -f xbessj0.o xbessj0.out; \
	fi

xbessj1.reslt : ../src/xbessj1.c $(ILIB)
	-if [ -r ../bin/xbessj1.out ]; then \
	../bin/xbessj1.out > xbessj1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessj1.out ../src/xbessj1.c $(ILIB) -lm && \
	xbessj1.out > xbessj1.reslt && /bin/rm -f xbessj1.o xbessj1.out; \
	fi

xbessjy.reslt : ../src/xbessjy.c $(ILIB)
	-if [ -r ../bin/xbessjy.out ]; then \
	../bin/xbessjy.out > xbessjy.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessjy.out ../src/xbessjy.c $(ILIB) -lm && \
	xbessjy.out > xbessjy.reslt && /bin/rm -f xbessjy.o xbessjy.out; \
	fi

xbessk.reslt : ../src/xbessk.c $(ILIB)
	-if [ -r ../bin/xbessk.out ]; then \
	../bin/xbessk.out > xbessk.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessk.out ../src/xbessk.c $(ILIB) -lm && \
	xbessk.out > xbessk.reslt && /bin/rm -f xbessk.o xbessk.out; \
	fi

xbessk0.reslt : ../src/xbessk0.c $(ILIB)
	-if [ -r ../bin/xbessk0.out ]; then \
	../bin/xbessk0.out > xbessk0.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessk0.out ../src/xbessk0.c $(ILIB) -lm && \
	xbessk0.out > xbessk0.reslt && /bin/rm -f xbessk0.o xbessk0.out; \
	fi

xbessk1.reslt : ../src/xbessk1.c $(ILIB)
	-if [ -r ../bin/xbessk1.out ]; then \
	../bin/xbessk1.out > xbessk1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessk1.out ../src/xbessk1.c $(ILIB) -lm && \
	xbessk1.out > xbessk1.reslt && /bin/rm -f xbessk1.o xbessk1.out; \
	fi

xbessy.reslt : ../src/xbessy.c $(ILIB)
	-if [ -r ../bin/xbessy.out ]; then \
	../bin/xbessy.out > xbessy.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessy.out ../src/xbessy.c $(ILIB) -lm && \
	xbessy.out > xbessy.reslt && /bin/rm -f xbessy.o xbessy.out; \
	fi

xbessy0.reslt : ../src/xbessy0.c $(ILIB)
	-if [ -r ../bin/xbessy0.out ]; then \
	../bin/xbessy0.out > xbessy0.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessy0.out ../src/xbessy0.c $(ILIB) -lm && \
	xbessy0.out > xbessy0.reslt && /bin/rm -f xbessy0.o xbessy0.out; \
	fi

xbessy1.reslt : ../src/xbessy1.c $(ILIB)
	-if [ -r ../bin/xbessy1.out ]; then \
	../bin/xbessy1.out > xbessy1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbessy1.out ../src/xbessy1.c $(ILIB) -lm && \
	xbessy1.out > xbessy1.reslt && /bin/rm -f xbessy1.o xbessy1.out; \
	fi

xbeta.reslt : ../src/xbeta.c $(ILIB)
	-if [ -r ../bin/xbeta.out ]; then \
	../bin/xbeta.out > xbeta.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbeta.out ../src/xbeta.c $(ILIB) -lm && \
	xbeta.out > xbeta.reslt && /bin/rm -f xbeta.o xbeta.out; \
	fi

xbetai.reslt : ../src/xbetai.c $(ILIB)
	-if [ -r ../bin/xbetai.out ]; then \
	../bin/xbetai.out > xbetai.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbetai.out ../src/xbetai.c $(ILIB) -lm && \
	xbetai.out > xbetai.reslt && /bin/rm -f xbetai.o xbetai.out; \
	fi

xbico.reslt : ../src/xbico.c $(ILIB)
	-if [ -r ../bin/xbico.out ]; then \
	../bin/xbico.out > xbico.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbico.out ../src/xbico.c $(ILIB) -lm && \
	xbico.out > xbico.reslt && /bin/rm -f xbico.o xbico.out; \
	fi

xbnldev.reslt : ../src/xbnldev.c $(ILIB)
	-if [ -r ../bin/xbnldev.out ]; then \
	../bin/xbnldev.out < ../responses/xbnldev.resp > xbnldev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbnldev.out ../src/xbnldev.c $(ILIB) -lm && \
	xbnldev.out < ../responses/xbnldev.resp > xbnldev.reslt && \
	/bin/rm -f xbnldev.o xbnldev.out; \
	fi

xbrent.reslt : ../src/xbrent.c $(ILIB)
	-if [ -r ../bin/xbrent.out ]; then \
	../bin/xbrent.out > xbrent.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbrent.out ../src/xbrent.c $(ILIB) -lm && \
	xbrent.out > xbrent.reslt && /bin/rm -f xbrent.o xbrent.out; \
	fi

xbroydn.reslt : ../src/xbroydn.c $(ILIB)
	-if [ -r ../bin/xbroydn.out ]; then \
	../bin/xbroydn.out > xbroydn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbroydn.out ../src/xbroydn.c $(ILIB) -lm && \
	xbroydn.out > xbroydn.reslt && /bin/rm -f xbroydn.o xbroydn.out; \
	fi

xbsstep.reslt : ../src/xbsstep.c $(ILIB)
	-if [ -r ../bin/xbsstep.out ]; then \
	../bin/xbsstep.out > xbsstep.reslt; \
	else \
	$(CC) $(CFLAGS) -o xbsstep.out ../src/xbsstep.c $(ILIB) -lm && \
	xbsstep.out > xbsstep.reslt && /bin/rm -f xbsstep.o xbsstep.out; \
	fi

xcaldat.reslt : ../src/xcaldat.c $(ILIB)
	-if [ -r ../bin/xcaldat.out ]; then \
	../bin/xcaldat.out > xcaldat.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcaldat.out ../src/xcaldat.c $(ILIB) -lm && \
	xcaldat.out > xcaldat.reslt && /bin/rm -f xcaldat.o xcaldat.out; \
	fi

xchder.reslt : ../src/xchder.c $(ILIB)
	-if [ -r ../bin/xchder.out ]; then \
	../bin/xchder.out < ../responses/xchder.resp > xchder.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchder.out ../src/xchder.c $(ILIB) -lm && \
	xchder.out < ../responses/xchder.resp > xchder.reslt && \
	/bin/rm -f xchder.o xchder.out; \
	fi

xchebev.reslt : ../src/xchebev.c $(ILIB)
	-if [ -r ../bin/xchebev.out ]; then \
	../bin/xchebev.out < ../responses/xchebev.resp > xchebev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchebev.out ../src/xchebev.c $(ILIB) -lm && \
	xchebev.out < ../responses/xchebev.resp > xchebev.reslt && \
	/bin/rm -f xchebev.o xchebev.out; \
	fi

xchebft.reslt : ../src/xchebft.c $(ILIB)
	-if [ -r ../bin/xchebft.out ]; then \
	../bin/xchebft.out < ../responses/xchebft.resp > xchebft.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchebft.out ../src/xchebft.c $(ILIB) -lm && \
	xchebft.out < ../responses/xchebft.resp > xchebft.reslt && \
	/bin/rm -f xchebft.o xchebft.out; \
	fi

xchebpc.reslt : ../src/xchebpc.c $(ILIB)
	-if [ -r ../bin/xchebpc.out ]; then \
	../bin/xchebpc.out < ../responses/xchebpc.resp > xchebpc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchebpc.out ../src/xchebpc.c $(ILIB) -lm && \
	xchebpc.out < ../responses/xchebpc.resp > xchebpc.reslt && \
	/bin/rm -f xchebpc.o xchebpc.out; \
	fi

xchint.reslt : ../src/xchint.c $(ILIB)
	-if [ -r ../bin/xchint.out ]; then \
	../bin/xchint.out < ../responses/xchint.resp > xchint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchint.out ../src/xchint.c $(ILIB) -lm && \
	xchint.out < ../responses/xchint.resp > xchint.reslt && \
	/bin/rm -f xchint.o xchint.out; \
	fi

xcholsl.reslt : ../src/xcholsl.c $(ILIB)
	-if [ -r ../bin/xcholsl.out ]; then \
	../bin/xcholsl.out > xcholsl.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcholsl.out ../src/xcholsl.c $(ILIB) -lm && \
	xcholsl.out > xcholsl.reslt && /bin/rm -f xcholsl.o xcholsl.out; \
	fi

xchsone.reslt : ../src/xchsone.c $(ILIB)
	-if [ -r ../bin/xchsone.out ]; then \
	../bin/xchsone.out > xchsone.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchsone.out ../src/xchsone.c $(ILIB) -lm && \
	xchsone.out > xchsone.reslt && /bin/rm -f xchsone.o xchsone.out; \
	fi

xchstwo.reslt : ../src/xchstwo.c $(ILIB)
	-if [ -r ../bin/xchstwo.out ]; then \
	../bin/xchstwo.out > xchstwo.reslt; \
	else \
	$(CC) $(CFLAGS) -o xchstwo.out ../src/xchstwo.c $(ILIB) -lm && \
	xchstwo.out > xchstwo.reslt && /bin/rm -f xchstwo.o xchstwo.out; \
	fi

xcisi.reslt : ../src/xcisi.c $(ILIB)
	-if [ -r ../bin/xcisi.out ]; then \
	../bin/xcisi.out > xcisi.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcisi.out ../src/xcisi.c $(ILIB) -lm && \
	xcisi.out > xcisi.reslt && /bin/rm -f xcisi.o xcisi.out; \
	fi

xcntab1.reslt : ../src/xcntab1.c $(ILIB)
	-if [ -r ../bin/xcntab1.out ]; then \
	../bin/xcntab1.out > xcntab1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcntab1.out ../src/xcntab1.c $(ILIB) -lm && \
	xcntab1.out > xcntab1.reslt && /bin/rm -f xcntab1.o xcntab1.out; \
	fi

xcntab2.reslt : ../src/xcntab2.c $(ILIB)
	-if [ -r ../bin/xcntab2.out ]; then \
	../bin/xcntab2.out > xcntab2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcntab2.out ../src/xcntab2.c $(ILIB) -lm && \
	xcntab2.out > xcntab2.reslt && /bin/rm -f xcntab2.o xcntab2.out; \
	fi

xconvlv.reslt : ../src/xconvlv.c $(ILIB)
	-if [ -r ../bin/xconvlv.out ]; then \
	../bin/xconvlv.out > xconvlv.reslt; \
	else \
	$(CC) $(CFLAGS) -o xconvlv.out ../src/xconvlv.c $(ILIB) -lm && \
	xconvlv.out > xconvlv.reslt && /bin/rm -f xconvlv.o xconvlv.out; \
	fi

xcorrel.reslt : ../src/xcorrel.c $(ILIB)
	-if [ -r ../bin/xcorrel.out ]; then \
	../bin/xcorrel.out > xcorrel.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcorrel.out ../src/xcorrel.c $(ILIB) -lm && \
	xcorrel.out > xcorrel.reslt && /bin/rm -f xcorrel.o xcorrel.out; \
	fi

xcosft1.reslt : ../src/xcosft1.c $(ILIB)
	-if [ -r ../bin/xcosft1.out ]; then \
	../bin/xcosft1.out < ../responses/xcosft1.resp > xcosft1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcosft1.out ../src/xcosft1.c $(ILIB) -lm && \
	xcosft1.out < ../responses/xcosft1.resp > xcosft1.reslt && \
	/bin/rm -f xcosft1.o xcosft1.out; \
	fi

xcosft2.reslt : ../src/xcosft2.c $(ILIB)
	-if [ -r ../bin/xcosft2.out ]; then \
	../bin/xcosft2.out < ../responses/xcosft2.resp > xcosft2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcosft2.out ../src/xcosft2.c $(ILIB) -lm && \
	xcosft2.out < ../responses/xcosft2.resp > xcosft2.reslt && \
	/bin/rm -f xcosft2.o xcosft2.out; \
	fi

xcovsrt.reslt : ../src/xcovsrt.c $(ILIB)
	-if [ -r ../bin/xcovsrt.out ]; then \
	../bin/xcovsrt.out < ../responses/xcovsrt.resp > xcovsrt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcovsrt.out ../src/xcovsrt.c $(ILIB) -lm && \
	xcovsrt.out < ../responses/xcovsrt.resp > xcovsrt.reslt && \
	/bin/rm -f xcovsrt.o xcovsrt.out; \
	fi

xcrank.reslt : ../src/xcrank.c $(ILIB)
	-if [ -r ../bin/xcrank.out ]; then \
	../bin/xcrank.out < ../responses/xcrank.resp > xcrank.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcrank.out ../src/xcrank.c $(ILIB) -lm && \
	xcrank.out < ../responses/xcrank.resp > xcrank.reslt && \
	/bin/rm -f xcrank.o xcrank.out; \
	fi

xcyclic.reslt : ../src/xcyclic.c $(ILIB)
	-if [ -r ../bin/xcyclic.out ]; then \
	../bin/xcyclic.out > xcyclic.reslt; \
	else \
	$(CC) $(CFLAGS) -o xcyclic.out ../src/xcyclic.c $(ILIB) -lm && \
	xcyclic.out > xcyclic.reslt && /bin/rm -f xcyclic.o xcyclic.out; \
	fi

xdawson.reslt : ../src/xdawson.c $(ILIB)
	-if [ -r ../bin/xdawson.out ]; then \
	../bin/xdawson.out > xdawson.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdawson.out ../src/xdawson.c $(ILIB) -lm && \
	xdawson.out > xdawson.reslt && /bin/rm -f xdawson.o xdawson.out; \
	fi

xdbrent.reslt : ../src/xdbrent.c $(ILIB)
	-if [ -r ../bin/xdbrent.out ]; then \
	../bin/xdbrent.out > xdbrent.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdbrent.out ../src/xdbrent.c $(ILIB) -lm && \
	xdbrent.out > xdbrent.reslt && /bin/rm -f xdbrent.o xdbrent.out; \
	fi

xddpoly.reslt : ../src/xddpoly.c $(ILIB)
	-if [ -r ../bin/xddpoly.out ]; then \
	../bin/xddpoly.out < ../responses/xddpoly.resp > xddpoly.reslt; \
	else \
	$(CC) $(CFLAGS) -o xddpoly.out ../src/xddpoly.c $(ILIB) -lm && \
	xddpoly.out < ../responses/xddpoly.resp > xddpoly.reslt && \
	/bin/rm -f xddpoly.o xddpoly.out; \
	fi

xdecchk.reslt : ../src/xdecchk.c $(ILIB)
	-if [ -r ../bin/xdecchk.out ]; then \
	../bin/xdecchk.out < ../responses/xdecchk.resp > xdecchk.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdecchk.out ../src/xdecchk.c $(ILIB) -lm && \
	xdecchk.out < ../responses/xdecchk.resp > xdecchk.reslt && \
	/bin/rm -f xdecchk.o xdecchk.out; \
	fi

xdf1dim.reslt : ../src/xdf1dim.c $(ILIB)
	-if [ -r ../bin/xdf1dim.out ]; then \
	../bin/xdf1dim.out < ../responses/xdf1dim.resp > xdf1dim.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdf1dim.out ../src/xdf1dim.c $(ILIB) -lm && \
	xdf1dim.out < ../responses/xdf1dim.resp > xdf1dim.reslt && \
	/bin/rm -f xdf1dim.o xdf1dim.out; \
	fi

xdfpmin.reslt : ../src/xdfpmin.c $(ILIB)
	-if [ -r ../bin/xdfpmin.out ]; then \
	../bin/xdfpmin.out > xdfpmin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdfpmin.out ../src/xdfpmin.c $(ILIB) -lm && \
	xdfpmin.out > xdfpmin.reslt && /bin/rm -f xdfpmin.o xdfpmin.out; \
	fi

xdfridr.reslt : ../src/xdfridr.c $(ILIB)
	-if [ -r ../bin/xdfridr.out ]; then \
	../bin/xdfridr.out < ../responses/xdfridr.resp > xdfridr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdfridr.out ../src/xdfridr.c $(ILIB) -lm && \
	xdfridr.out < ../responses/xdfridr.resp > xdfridr.reslt && \
	/bin/rm -f xdfridr.o xdfridr.out; \
	fi

xdftint.reslt : ../src/xdftint.c $(ILIB)
	-if [ -r ../bin/xdftint.out ]; then \
	../bin/xdftint.out < ../responses/xdftint.resp > xdftint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdftint.out ../src/xdftint.c $(ILIB) -lm && \
	xdftint.out < ../responses/xdftint.resp > xdftint.reslt && \
	/bin/rm -f xdftint.o xdftint.out; \
	fi

xdlinmin.reslt : ../src/xdlinmin.c $(ILIB)
	-if [ -r ../bin/xdlinmin.out ]; then \
	../bin/xdlinmin.out > xdlinmin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xdlinmin.out ../src/xdlinmin.c $(ILIB) -lm && \
	xdlinmin.out > xdlinmin.reslt && /bin/rm -f xdlinmin.o xdlinmin.out; \
	fi

xeclass.reslt : ../src/xeclass.c $(ILIB)
	-if [ -r ../bin/xeclass.out ]; then \
	../bin/xeclass.out > xeclass.reslt; \
	else \
	$(CC) $(CFLAGS) -o xeclass.out ../src/xeclass.c $(ILIB) -lm && \
	xeclass.out > xeclass.reslt && /bin/rm -f xeclass.o xeclass.out; \
	fi

xeclazz.reslt : ../src/xeclazz.c $(ILIB)
	-if [ -r ../bin/xeclazz.out ]; then \
	../bin/xeclazz.out > xeclazz.reslt; \
	else \
	$(CC) $(CFLAGS) -o xeclazz.out ../src/xeclazz.c $(ILIB) -lm && \
	xeclazz.out > xeclazz.reslt && /bin/rm -f xeclazz.o xeclazz.out; \
	fi

xei.reslt : ../src/xei.c $(ILIB)
	-if [ -r ../bin/xei.out ]; then \
	../bin/xei.out > xei.reslt; \
	else \
	$(CC) $(CFLAGS) -o xei.out ../src/xei.c $(ILIB) -lm && \
	xei.out > xei.reslt && /bin/rm -f xei.o xei.out; \
	fi

xeigsrt.reslt : ../src/xeigsrt.c $(ILIB)
	-if [ -r ../bin/xeigsrt.out ]; then \
	../bin/xeigsrt.out > xeigsrt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xeigsrt.out ../src/xeigsrt.c $(ILIB) -lm && \
	xeigsrt.out > xeigsrt.reslt && /bin/rm -f xeigsrt.o xeigsrt.out; \
	fi

xelle.reslt : ../src/xelle.c $(ILIB)
	-if [ -r ../bin/xelle.out ]; then \
	../bin/xelle.out > xelle.reslt; \
	else \
	$(CC) $(CFLAGS) -o xelle.out ../src/xelle.c $(ILIB) -lm && \
	xelle.out > xelle.reslt && /bin/rm -f xelle.o xelle.out; \
	fi

xellf.reslt : ../src/xellf.c $(ILIB)
	-if [ -r ../bin/xellf.out ]; then \
	../bin/xellf.out > xellf.reslt; \
	else \
	$(CC) $(CFLAGS) -o xellf.out ../src/xellf.c $(ILIB) -lm && \
	xellf.out > xellf.reslt && /bin/rm -f xellf.o xellf.out; \
	fi

xellpi.reslt : ../src/xellpi.c $(ILIB)
	-if [ -r ../bin/xellpi.out ]; then \
	../bin/xellpi.out > xellpi.reslt; \
	else \
	$(CC) $(CFLAGS) -o xellpi.out ../src/xellpi.c $(ILIB) -lm && \
	xellpi.out > xellpi.reslt && /bin/rm -f xellpi.o xellpi.out; \
	fi

xelmhes.reslt : ../src/xelmhes.c $(ILIB)
	-if [ -r ../bin/xelmhes.out ]; then \
	../bin/xelmhes.out > xelmhes.reslt; \
	else \
	$(CC) $(CFLAGS) -o xelmhes.out ../src/xelmhes.c $(ILIB) -lm && \
	xelmhes.out > xelmhes.reslt && /bin/rm -f xelmhes.o xelmhes.out; \
	fi

xerfcc.reslt : ../src/xerfcc.c $(ILIB)
	-if [ -r ../bin/xerfcc.out ]; then \
	../bin/xerfcc.out > xerfcc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xerfcc.out ../src/xerfcc.c $(ILIB) -lm && \
	xerfcc.out > xerfcc.reslt && /bin/rm -f xerfcc.o xerfcc.out; \
	fi

xerff.reslt : ../src/xerff.c $(ILIB)
	-if [ -r ../bin/xerff.out ]; then \
	../bin/xerff.out > xerff.reslt; \
	else \
	$(CC) $(CFLAGS) -o xerff.out ../src/xerff.c $(ILIB) -lm && \
	xerff.out > xerff.reslt && /bin/rm -f xerff.o xerff.out; \
	fi

xerffc.reslt : ../src/xerffc.c $(ILIB)
	-if [ -r ../bin/xerffc.out ]; then \
	../bin/xerffc.out > xerffc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xerffc.out ../src/xerffc.c $(ILIB) -lm && \
	xerffc.out > xerffc.reslt && /bin/rm -f xerffc.o xerffc.out; \
	fi

xeulsum.reslt : ../src/xeulsum.c $(ILIB)
	-if [ -r ../bin/xeulsum.out ]; then \
	../bin/xeulsum.out < ../responses/xeulsum.resp > xeulsum.reslt; \
	else \
	$(CC) $(CFLAGS) -o xeulsum.out ../src/xeulsum.c $(ILIB) -lm && \
	xeulsum.out < ../responses/xeulsum.resp > xeulsum.reslt && \
	/bin/rm -f xeulsum.o xeulsum.out; \
	fi

xevlmem.reslt : ../src/xevlmem.c $(ILIB)
	-if [ -r ../bin/xevlmem.out ]; then \
	../bin/xevlmem.out > xevlmem.reslt; \
	else \
	$(CC) $(CFLAGS) -o xevlmem.out ../src/xevlmem.c $(ILIB) -lm && \
	xevlmem.out > xevlmem.reslt && /bin/rm -f xevlmem.o xevlmem.out; \
	fi

xexpdev.reslt : ../src/xexpdev.c $(ILIB)
	-if [ -r ../bin/xexpdev.out ]; then \
	../bin/xexpdev.out > xexpdev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xexpdev.out ../src/xexpdev.c $(ILIB) -lm && \
	xexpdev.out > xexpdev.reslt && /bin/rm -f xexpdev.o xexpdev.out; \
	fi

xexpint.reslt : ../src/xexpint.c $(ILIB)
	-if [ -r ../bin/xexpint.out ]; then \
	../bin/xexpint.out > xexpint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xexpint.out ../src/xexpint.c $(ILIB) -lm && \
	xexpint.out > xexpint.reslt && /bin/rm -f xexpint.o xexpint.out; \
	fi

xf1dim.reslt : ../src/xf1dim.c $(ILIB)
	-if [ -r ../bin/xf1dim.out ]; then \
	../bin/xf1dim.out < ../responses/xf1dim.resp > xf1dim.reslt; \
	else \
	$(CC) $(CFLAGS) -o xf1dim.out ../src/xf1dim.c $(ILIB) -lm && \
	xf1dim.out < ../responses/xf1dim.resp > xf1dim.reslt && \
	/bin/rm -f xf1dim.o xf1dim.out; \
	fi

xfactln.reslt : ../src/xfactln.c $(ILIB)
	-if [ -r ../bin/xfactln.out ]; then \
	../bin/xfactln.out > xfactln.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfactln.out ../src/xfactln.c $(ILIB) -lm && \
	xfactln.out > xfactln.reslt && /bin/rm -f xfactln.o xfactln.out; \
	fi

xfactrl.reslt : ../src/xfactrl.c $(ILIB)
	-if [ -r ../bin/xfactrl.out ]; then \
	../bin/xfactrl.out > xfactrl.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfactrl.out ../src/xfactrl.c $(ILIB) -lm && \
	xfactrl.out > xfactrl.reslt && /bin/rm -f xfactrl.o xfactrl.out; \
	fi

xfasper.reslt : ../src/xfasper.c $(ILIB)
	-if [ -r ../bin/xfasper.out ]; then \
	../bin/xfasper.out > xfasper.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfasper.out ../src/xfasper.c $(ILIB) -lm && \
	xfasper.out > xfasper.reslt && /bin/rm -f xfasper.o xfasper.out; \
	fi

xfgauss.reslt : ../src/xfgauss.c $(ILIB)
	-if [ -r ../bin/xfgauss.out ]; then \
	../bin/xfgauss.out > xfgauss.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfgauss.out ../src/xfgauss.c $(ILIB) -lm && \
	xfgauss.out > xfgauss.reslt && /bin/rm -f xfgauss.o xfgauss.out; \
	fi

xfit.reslt : ../src/xfit.c $(ILIB)
	-if [ -r ../bin/xfit.out ]; then \
	../bin/xfit.out > xfit.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfit.out ../src/xfit.c $(ILIB) -lm && \
	xfit.out > xfit.reslt && /bin/rm -f xfit.o xfit.out; \
	fi

xfitexy.reslt : ../src/xfitexy.c $(ILIB)
	-if [ -r ../bin/xfitexy.out ]; then \
	../bin/xfitexy.out > xfitexy.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfitexy.out ../src/xfitexy.c $(ILIB) -lm && \
	xfitexy.out > xfitexy.reslt && /bin/rm -f xfitexy.o xfitexy.out; \
	fi

xfixrts.reslt : ../src/xfixrts.c $(ILIB)
	-if [ -r ../bin/xfixrts.out ]; then \
	../bin/xfixrts.out > xfixrts.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfixrts.out ../src/xfixrts.c $(ILIB) -lm && \
	xfixrts.out > xfixrts.reslt && /bin/rm -f xfixrts.o xfixrts.out; \
	fi

xfleg.reslt : ../src/xfleg.c $(ILIB)
	-if [ -r ../bin/xfleg.out ]; then \
	../bin/xfleg.out > xfleg.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfleg.out ../src/xfleg.c $(ILIB) -lm && \
	xfleg.out > xfleg.reslt && /bin/rm -f xfleg.o xfleg.out; \
	fi

xflmoon.reslt : ../src/xflmoon.c $(ILIB)
	-if [ -r ../bin/xflmoon.out ]; then \
	../bin/xflmoon.out < ../responses/xflmoon.resp > xflmoon.reslt; \
	else \
	$(CC) $(CFLAGS) -o xflmoon.out ../src/xflmoon.c $(ILIB) -lm && \
	xflmoon.out < ../responses/xflmoon.resp > xflmoon.reslt && \
	/bin/rm -f xflmoon.o xflmoon.out; \
	fi

xfour1.reslt : ../src/xfour1.c $(ILIB)
	-if [ -r ../bin/xfour1.out ]; then \
	../bin/xfour1.out < ../responses/xfour1.resp > xfour1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfour1.out ../src/xfour1.c $(ILIB) -lm && \
	xfour1.out < ../responses/xfour1.resp > xfour1.reslt && \
	/bin/rm -f xfour1.o xfour1.out; \
	fi

xfourfs.reslt : ../src/xfourfs.c $(ILIB)
	-if [ -r ../bin/xfourfs.out ]; then \
	../bin/xfourfs.out > xfourfs.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfourfs.out ../src/xfourfs.c $(ILIB) -lm && \
	xfourfs.out > xfourfs.reslt && /bin/rm -f xfourfs.o xfourfs.out; \
	fi

xfourn.reslt : ../src/xfourn.c $(ILIB)
	-if [ -r ../bin/xfourn.out ]; then \
	../bin/xfourn.out > xfourn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfourn.out ../src/xfourn.c $(ILIB) -lm && \
	xfourn.out > xfourn.reslt && /bin/rm -f xfourn.o xfourn.out; \
	fi

xfpoly.reslt : ../src/xfpoly.c $(ILIB)
	-if [ -r ../bin/xfpoly.out ]; then \
	../bin/xfpoly.out > xfpoly.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfpoly.out ../src/xfpoly.c $(ILIB) -lm && \
	xfpoly.out > xfpoly.reslt && /bin/rm -f xfpoly.o xfpoly.out; \
	fi

xfred2.reslt : ../src/xfred2.c $(ILIB)
	-if [ -r ../bin/xfred2.out ]; then \
	../bin/xfred2.out > xfred2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfred2.out ../src/xfred2.c $(ILIB) -lm && \
	xfred2.out > xfred2.reslt && /bin/rm -f xfred2.o xfred2.out; \
	fi

xfredin.reslt : ../src/xfredin.c $(ILIB)
	-if [ -r ../bin/xfredin.out ]; then \
	../bin/xfredin.out < ../responses/xfredin.resp > xfredin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfredin.out ../src/xfredin.c $(ILIB) -lm && \
	xfredin.out < ../responses/xfredin.resp > xfredin.reslt && \
	/bin/rm -f xfredin.o xfredin.out; \
	fi

xfrenel.reslt : ../src/xfrenel.c $(ILIB)
	-if [ -r ../bin/xfrenel.out ]; then \
	../bin/xfrenel.out > xfrenel.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfrenel.out ../src/xfrenel.c $(ILIB) -lm && \
	xfrenel.out > xfrenel.reslt && /bin/rm -f xfrenel.o xfrenel.out; \
	fi

xfrprmn.reslt : ../src/xfrprmn.c $(ILIB)
	-if [ -r ../bin/xfrprmn.out ]; then \
	../bin/xfrprmn.out > xfrprmn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xfrprmn.out ../src/xfrprmn.c $(ILIB) -lm && \
	xfrprmn.out > xfrprmn.reslt && /bin/rm -f xfrprmn.o xfrprmn.out; \
	fi

xftest.reslt : ../src/xftest.c $(ILIB)
	-if [ -r ../bin/xftest.out ]; then \
	../bin/xftest.out > xftest.reslt; \
	else \
	$(CC) $(CFLAGS) -o xftest.out ../src/xftest.c $(ILIB) -lm && \
	xftest.out > xftest.reslt && /bin/rm -f xftest.o xftest.out; \
	fi

xgamdev.reslt : ../src/xgamdev.c $(ILIB)
	-if [ -r ../bin/xgamdev.out ]; then \
	../bin/xgamdev.out < ../responses/xgamdev.resp > xgamdev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgamdev.out ../src/xgamdev.c $(ILIB) -lm && \
	xgamdev.out < ../responses/xgamdev.resp > xgamdev.reslt && \
	/bin/rm -f xgamdev.o xgamdev.out; \
	fi

xgammln.reslt : ../src/xgammln.c $(ILIB)
	-if [ -r ../bin/xgammln.out ]; then \
	../bin/xgammln.out > xgammln.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgammln.out ../src/xgammln.c $(ILIB) -lm && \
	xgammln.out > xgammln.reslt && /bin/rm -f xgammln.o xgammln.out; \
	fi

xgammp.reslt : ../src/xgammp.c $(ILIB)
	-if [ -r ../bin/xgammp.out ]; then \
	../bin/xgammp.out > xgammp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgammp.out ../src/xgammp.c $(ILIB) -lm && \
	xgammp.out > xgammp.reslt && /bin/rm -f xgammp.o xgammp.out; \
	fi

xgammq.reslt : ../src/xgammq.c $(ILIB)
	-if [ -r ../bin/xgammq.out ]; then \
	../bin/xgammq.out > xgammq.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgammq.out ../src/xgammq.c $(ILIB) -lm && \
	xgammq.out > xgammq.reslt && /bin/rm -f xgammq.o xgammq.out; \
	fi

xgasdev.reslt : ../src/xgasdev.c $(ILIB)
	-if [ -r ../bin/xgasdev.out ]; then \
	../bin/xgasdev.out > xgasdev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgasdev.out ../src/xgasdev.c $(ILIB) -lm && \
	xgasdev.out > xgasdev.reslt && /bin/rm -f xgasdev.o xgasdev.out; \
	fi

xgaucof.reslt : ../src/xgaucof.c $(ILIB)
	-if [ -r ../bin/xgaucof.out ]; then \
	../bin/xgaucof.out < ../responses/xgaucof.resp > xgaucof.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgaucof.out ../src/xgaucof.c $(ILIB) -lm && \
	xgaucof.out < ../responses/xgaucof.resp > xgaucof.reslt && \
	/bin/rm -f xgaucof.o xgaucof.out; \
	fi

xgauher.reslt : ../src/xgauher.c $(ILIB)
	-if [ -r ../bin/xgauher.out ]; then \
	../bin/xgauher.out < ../responses/xgauher.resp > xgauher.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgauher.out ../src/xgauher.c $(ILIB) -lm && \
	xgauher.out < ../responses/xgauher.resp > xgauher.reslt && \
	/bin/rm -f xgauher.o xgauher.out; \
	fi

xgaujac.reslt : ../src/xgaujac.c $(ILIB)
	-if [ -r ../bin/xgaujac.out ]; then \
	../bin/xgaujac.out < ../responses/xgaujac.resp > xgaujac.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgaujac.out ../src/xgaujac.c $(ILIB) -lm && \
	xgaujac.out < ../responses/xgaujac.resp > xgaujac.reslt && \
	/bin/rm -f xgaujac.o xgaujac.out; \
	fi

xgaulag.reslt : ../src/xgaulag.c $(ILIB)
	-if [ -r ../bin/xgaulag.out ]; then \
	../bin/xgaulag.out < ../responses/xgaulag.resp > xgaulag.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgaulag.out ../src/xgaulag.c $(ILIB) -lm && \
	xgaulag.out < ../responses/xgaulag.resp > xgaulag.reslt && \
	/bin/rm -f xgaulag.o xgaulag.out; \
	fi

xgauleg.reslt : ../src/xgauleg.c $(ILIB)
	-if [ -r ../bin/xgauleg.out ]; then \
	../bin/xgauleg.out > xgauleg.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgauleg.out ../src/xgauleg.c $(ILIB) -lm && \
	xgauleg.out > xgauleg.reslt && /bin/rm -f xgauleg.o xgauleg.out; \
	fi

xgaussj.reslt : ../src/xgaussj.c $(ILIB)
	-if [ -r ../bin/xgaussj.out ]; then \
	../bin/xgaussj.out < ../responses/xgaussj.resp > xgaussj.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgaussj.out ../src/xgaussj.c $(ILIB) -lm && \
	xgaussj.out < ../responses/xgaussj.resp > xgaussj.reslt && \
	/bin/rm -f xgaussj.o xgaussj.out; \
	fi

xgcf.reslt : ../src/xgcf.c $(ILIB)
	-if [ -r ../bin/xgcf.out ]; then \
	../bin/xgcf.out > xgcf.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgcf.out ../src/xgcf.c $(ILIB) -lm && \
	xgcf.out > xgcf.reslt && /bin/rm -f xgcf.o xgcf.out; \
	fi

xgolden.reslt : ../src/xgolden.c $(ILIB)
	-if [ -r ../bin/xgolden.out ]; then \
	../bin/xgolden.out > xgolden.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgolden.out ../src/xgolden.c $(ILIB) -lm && \
	xgolden.out > xgolden.reslt && /bin/rm -f xgolden.o xgolden.out; \
	fi

xgser.reslt : ../src/xgser.c $(ILIB)
	-if [ -r ../bin/xgser.out ]; then \
	../bin/xgser.out > xgser.reslt; \
	else \
	$(CC) $(CFLAGS) -o xgser.out ../src/xgser.c $(ILIB) -lm && \
	xgser.out > xgser.reslt && /bin/rm -f xgser.o xgser.out; \
	fi

xhpsel.reslt : ../src/xhpsel.c $(ILIB)
	-if [ -r ../bin/xhpsel.out ]; then \
	../bin/xhpsel.out < ../responses/xhpsel.resp > xhpsel.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhpsel.out ../src/xhpsel.c $(ILIB) -lm && \
	xhpsel.out < ../responses/xhpsel.resp > xhpsel.reslt && \
	/bin/rm -f xhpsel.o xhpsel.out; \
	fi

xhpsort.reslt : ../src/xhpsort.c $(ILIB)
	-if [ -r ../bin/xhpsort.out ]; then \
	../bin/xhpsort.out > xhpsort.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhpsort.out ../src/xhpsort.c $(ILIB) -lm && \
	xhpsort.out > xhpsort.reslt && /bin/rm -f xhpsort.o xhpsort.out; \
	fi

xhqr.reslt : ../src/xhqr.c $(ILIB)
	-if [ -r ../bin/xhqr.out ]; then \
	../bin/xhqr.out > xhqr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhqr.out ../src/xhqr.c $(ILIB) -lm && \
	xhqr.out > xhqr.reslt && /bin/rm -f xhqr.o xhqr.out; \
	fi

xhuffman.reslt : ../src/xhuffman.c $(ILIB)
	-if [ -r ../bin/xhuffman.out ]; then \
	../bin/xhuffman.out < ../responses/xhuffman.resp > xhuffman.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhuffman.out ../src/xhuffman.c $(ILIB) -lm && \
	xhuffman.out < ../responses/xhuffman.resp > xhuffman.reslt && \
	/bin/rm -f xhuffman.o xhuffman.out; \
	fi

xhunt.reslt : ../src/xhunt.c $(ILIB)
	-if [ -r ../bin/xhunt.out ]; then \
	../bin/xhunt.out > xhunt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhunt.out ../src/xhunt.c $(ILIB) -lm && \
	xhunt.out > xhunt.reslt && /bin/rm -f xhunt.o xhunt.out; \
	fi

xhypgeo.reslt : ../src/xhypgeo.c $(ILIB)
	-if [ -r ../bin/xhypgeo.out ]; then \
	../bin/xhypgeo.out < ../responses/xhypgeo.resp > xhypgeo.reslt; \
	else \
	$(CC) $(CFLAGS) -o xhypgeo.out ../src/xhypgeo.c $(ILIB) -lm && \
	xhypgeo.out < ../responses/xhypgeo.resp > xhypgeo.reslt && \
	/bin/rm -f xhypgeo.o xhypgeo.out; \
	fi

xicrc.reslt : ../src/xicrc.c $(ILIB)
	-if [ -r ../bin/xicrc.out ]; then \
	../bin/xicrc.out < ../responses/xicrc.resp > xicrc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xicrc.out ../src/xicrc.c $(ILIB) -lm && \
	xicrc.out < ../responses/xicrc.resp > xicrc.reslt && \
	/bin/rm -f xicrc.o xicrc.out; \
	fi

xigray.reslt : ../src/xigray.c $(ILIB)
	-if [ -r ../bin/xigray.out ]; then \
	../bin/xigray.out < ../responses/xigray.resp > xigray.reslt; \
	else \
	$(CC) $(CFLAGS) -o xigray.out ../src/xigray.c $(ILIB) -lm && \
	xigray.out < ../responses/xigray.resp > xigray.reslt && \
	/bin/rm -f xigray.o xigray.out; \
	fi

xindexx.reslt : ../src/xindexx.c $(ILIB)
	-if [ -r ../bin/xindexx.out ]; then \
	../bin/xindexx.out > xindexx.reslt; \
	else \
	$(CC) $(CFLAGS) -o xindexx.out ../src/xindexx.c $(ILIB) -lm && \
	xindexx.out > xindexx.reslt && /bin/rm -f xindexx.o xindexx.out; \
	fi

xirbit1.reslt : ../src/xirbit1.c $(ILIB)
	-if [ -r ../bin/xirbit1.out ]; then \
	../bin/xirbit1.out > xirbit1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xirbit1.out ../src/xirbit1.c $(ILIB) -lm && \
	xirbit1.out > xirbit1.reslt && /bin/rm -f xirbit1.o xirbit1.out; \
	fi

xirbit2.reslt : ../src/xirbit2.c $(ILIB)
	-if [ -r ../bin/xirbit2.out ]; then \
	../bin/xirbit2.out > xirbit2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xirbit2.out ../src/xirbit2.c $(ILIB) -lm && \
	xirbit2.out > xirbit2.reslt && /bin/rm -f xirbit2.o xirbit2.out; \
	fi

xjacobi.reslt : ../src/xjacobi.c $(ILIB)
	-if [ -r ../bin/xjacobi.out ]; then \
	../bin/xjacobi.out < ../responses/xjacobi.resp > xjacobi.reslt; \
	else \
	$(CC) $(CFLAGS) -o xjacobi.out ../src/xjacobi.c $(ILIB) -lm && \
	xjacobi.out < ../responses/xjacobi.resp > xjacobi.reslt && \
	/bin/rm -f xjacobi.o xjacobi.out; \
	fi

xjulday.reslt : ../src/xjulday.c $(ILIB)
	-if [ -r ../bin/xjulday.out ]; then \
	../bin/xjulday.out < ../responses/xjulday.resp > xjulday.reslt; \
	else \
	$(CC) $(CFLAGS) -o xjulday.out ../src/xjulday.c $(ILIB) -lm && \
	xjulday.out < ../responses/xjulday.resp > xjulday.reslt && \
	/bin/rm -f xjulday.o xjulday.out; \
	fi

xkendl1.reslt : ../src/xkendl1.c $(ILIB)
	-if [ -r ../bin/xkendl1.out ]; then \
	../bin/xkendl1.out > xkendl1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xkendl1.out ../src/xkendl1.c $(ILIB) -lm && \
	xkendl1.out > xkendl1.reslt && /bin/rm -f xkendl1.o xkendl1.out; \
	fi

xkendl2.reslt : ../src/xkendl2.c $(ILIB)
	-if [ -r ../bin/xkendl2.out ]; then \
	../bin/xkendl2.out > xkendl2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xkendl2.out ../src/xkendl2.c $(ILIB) -lm && \
	xkendl2.out > xkendl2.reslt && /bin/rm -f xkendl2.o xkendl2.out; \
	fi

xks2d1s.reslt : ../src/xks2d1s.c $(ILIB)
	-if [ -r ../bin/xks2d1s.out ]; then \
	../bin/xks2d1s.out < ../responses/xks2d1s.resp > xks2d1s.reslt; \
	else \
	$(CC) $(CFLAGS) -o xks2d1s.out ../src/xks2d1s.c $(ILIB) -lm && \
	xks2d1s.out < ../responses/xks2d1s.resp > xks2d1s.reslt && \
	/bin/rm -f xks2d1s.o xks2d1s.out; \
	fi

xks2d2s.reslt : ../src/xks2d2s.c $(ILIB)
	-if [ -r ../bin/xks2d2s.out ]; then \
	../bin/xks2d2s.out < ../responses/xks2d2s.resp > xks2d2s.reslt; \
	else \
	$(CC) $(CFLAGS) -o xks2d2s.out ../src/xks2d2s.c $(ILIB) -lm && \
	xks2d2s.out < ../responses/xks2d2s.resp > xks2d2s.reslt && \
	/bin/rm -f xks2d2s.o xks2d2s.out; \
	fi

xksone.reslt : ../src/xksone.c $(ILIB)
	-if [ -r ../bin/xksone.out ]; then \
	../bin/xksone.out > xksone.reslt; \
	else \
	$(CC) $(CFLAGS) -o xksone.out ../src/xksone.c $(ILIB) -lm && \
	xksone.out > xksone.reslt && /bin/rm -f xksone.o xksone.out; \
	fi

xkstwo.reslt : ../src/xkstwo.c $(ILIB)
	-if [ -r ../bin/xkstwo.out ]; then \
	../bin/xkstwo.out > xkstwo.reslt; \
	else \
	$(CC) $(CFLAGS) -o xkstwo.out ../src/xkstwo.c $(ILIB) -lm && \
	xkstwo.out > xkstwo.reslt && /bin/rm -f xkstwo.o xkstwo.out; \
	fi

xlaguer.reslt : ../src/xlaguer.c $(ILIB)
	-if [ -r ../bin/xlaguer.out ]; then \
	../bin/xlaguer.out > xlaguer.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlaguer.out ../src/xlaguer.c $(ILIB) -lm && \
	xlaguer.out > xlaguer.reslt && /bin/rm -f xlaguer.o xlaguer.out; \
	fi

xlfit.reslt : ../src/xlfit.c $(ILIB)
	-if [ -r ../bin/xlfit.out ]; then \
	../bin/xlfit.out < ../responses/xlfit.resp > xlfit.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlfit.out ../src/xlfit.c $(ILIB) -lm && \
	xlfit.out < ../responses/xlfit.resp > xlfit.reslt && \
	/bin/rm -f xlfit.o xlfit.out; \
	fi

xlinbcg.reslt : ../src/xlinbcg.c $(ILIB)
	-if [ -r ../bin/xlinbcg.out ]; then \
	../bin/xlinbcg.out < ../responses/xlinbcg.resp > xlinbcg.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlinbcg.out ../src/xlinbcg.c $(ILIB) -lm && \
	xlinbcg.out < ../responses/xlinbcg.resp > xlinbcg.reslt && \
	/bin/rm -f xlinbcg.o xlinbcg.out; \
	fi

xlinmin.reslt : ../src/xlinmin.c $(ILIB)
	-if [ -r ../bin/xlinmin.out ]; then \
	../bin/xlinmin.out > xlinmin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlinmin.out ../src/xlinmin.c $(ILIB) -lm && \
	xlinmin.out > xlinmin.reslt && /bin/rm -f xlinmin.o xlinmin.out; \
	fi

xlocate.reslt : ../src/xlocate.c $(ILIB)
	-if [ -r ../bin/xlocate.out ]; then \
	../bin/xlocate.out > xlocate.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlocate.out ../src/xlocate.c $(ILIB) -lm && \
	xlocate.out > xlocate.reslt && /bin/rm -f xlocate.o xlocate.out; \
	fi

xlubksb.reslt : ../src/xlubksb.c $(ILIB)
	-if [ -r ../bin/xlubksb.out ]; then \
	../bin/xlubksb.out < ../responses/xlubksb.resp > xlubksb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xlubksb.out ../src/xlubksb.c $(ILIB) -lm && \
	xlubksb.out < ../responses/xlubksb.resp > xlubksb.reslt && \
	/bin/rm -f xlubksb.o xlubksb.out; \
	fi

xludcmp.reslt : ../src/xludcmp.c $(ILIB)
	-if [ -r ../bin/xludcmp.out ]; then \
	../bin/xludcmp.out < ../responses/xludcmp.resp > xludcmp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xludcmp.out ../src/xludcmp.c $(ILIB) -lm && \
	xludcmp.out < ../responses/xludcmp.resp > xludcmp.reslt && \
	/bin/rm -f xludcmp.o xludcmp.out; \
	fi

xmachar.reslt : ../src/xmachar.c $(ILIB)
	-if [ -r ../bin/xmachar.out ]; then \
	../bin/xmachar.out > xmachar.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmachar.out ../src/xmachar.c $(ILIB) -lm && \
	xmachar.out > xmachar.reslt && /bin/rm -f xmachar.o xmachar.out; \
	fi

xmedfit.reslt : ../src/xmedfit.c $(ILIB)
	-if [ -r ../bin/xmedfit.out ]; then \
	../bin/xmedfit.out > xmedfit.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmedfit.out ../src/xmedfit.c $(ILIB) -lm && \
	xmedfit.out > xmedfit.reslt && /bin/rm -f xmedfit.o xmedfit.out; \
	fi

xmemcof.reslt : ../src/xmemcof.c $(ILIB)
	-if [ -r ../bin/xmemcof.out ]; then \
	../bin/xmemcof.out > xmemcof.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmemcof.out ../src/xmemcof.c $(ILIB) -lm && \
	xmemcof.out > xmemcof.reslt && /bin/rm -f xmemcof.o xmemcof.out; \
	fi

xmgfas.reslt : ../src/xmgfas.c $(ILIB)
	-if [ -r ../bin/xmgfas.out ]; then \
	../bin/xmgfas.out > xmgfas.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmgfas.out ../src/xmgfas.c $(ILIB) -lm && \
	xmgfas.out > xmgfas.reslt && /bin/rm -f xmgfas.o xmgfas.out; \
	fi

xmglin.reslt : ../src/xmglin.c $(ILIB)
	-if [ -r ../bin/xmglin.out ]; then \
	../bin/xmglin.out > xmglin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmglin.out ../src/xmglin.c $(ILIB) -lm && \
	xmglin.out > xmglin.reslt && /bin/rm -f xmglin.o xmglin.out; \
	fi

xmidpnt.reslt : ../src/xmidpnt.c $(ILIB)
	-if [ -r ../bin/xmidpnt.out ]; then \
	../bin/xmidpnt.out > xmidpnt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmidpnt.out ../src/xmidpnt.c $(ILIB) -lm && \
	xmidpnt.out > xmidpnt.reslt && /bin/rm -f xmidpnt.o xmidpnt.out; \
	fi

xmiser.reslt : ../src/xmiser.c $(ILIB)
	-if [ -r ../bin/xmiser.out ]; then \
	../bin/xmiser.out < ../responses/xmiser.resp > xmiser.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmiser.out ../src/xmiser.c $(ILIB) -lm && \
	xmiser.out < ../responses/xmiser.resp > xmiser.reslt && \
	/bin/rm -f xmiser.o xmiser.out; \
	fi

xmmid.reslt : ../src/xmmid.c $(ILIB)
	-if [ -r ../bin/xmmid.out ]; then \
	../bin/xmmid.out < ../responses/xmmid.resp > xmmid.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmmid.out ../src/xmmid.c $(ILIB) -lm && \
	xmmid.out < ../responses/xmmid.resp > xmmid.reslt && \
	/bin/rm -f xmmid.o xmmid.out; \
	fi

xmnbrak.reslt : ../src/xmnbrak.c $(ILIB)
	-if [ -r ../bin/xmnbrak.out ]; then \
	../bin/xmnbrak.out > xmnbrak.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmnbrak.out ../src/xmnbrak.c $(ILIB) -lm && \
	xmnbrak.out > xmnbrak.reslt && /bin/rm -f xmnbrak.o xmnbrak.out; \
	fi

xmnewt.reslt : ../src/xmnewt.c $(ILIB)
	-if [ -r ../bin/xmnewt.out ]; then \
	../bin/xmnewt.out < ../responses/xmnewt.resp > xmnewt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmnewt.out ../src/xmnewt.c $(ILIB) -lm && \
	xmnewt.out < ../responses/xmnewt.resp > xmnewt.reslt && \
	/bin/rm -f xmnewt.o xmnewt.out; \
	fi

xmoment.reslt : ../src/xmoment.c $(ILIB)
	-if [ -r ../bin/xmoment.out ]; then \
	../bin/xmoment.out > xmoment.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmoment.out ../src/xmoment.c $(ILIB) -lm && \
	xmoment.out > xmoment.reslt && /bin/rm -f xmoment.o xmoment.out; \
	fi

xmppi.reslt : ../src/xmppi.c $(ILIB)
	-if [ -r ../bin/xmppi.out ]; then \
	../bin/xmppi.out < ../responses/xmppi.resp > xmppi.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmppi.out ../src/xmppi.c $(ILIB) -lm && \
	xmppi.out < ../responses/xmppi.resp > xmppi.reslt && \
	/bin/rm -f xmppi.o xmppi.out; \
	fi

xmprove.reslt : ../src/xmprove.c $(ILIB)
	-if [ -r ../bin/xmprove.out ]; then \
	../bin/xmprove.out > xmprove.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmprove.out ../src/xmprove.c $(ILIB) -lm && \
	xmprove.out > xmprove.reslt && /bin/rm -f xmprove.o xmprove.out; \
	fi

xmrqcof.reslt : ../src/xmrqcof.c $(ILIB)
	-if [ -r ../bin/xmrqcof.out ]; then \
	../bin/xmrqcof.out > xmrqcof.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmrqcof.out ../src/xmrqcof.c $(ILIB) -lm && \
	xmrqcof.out > xmrqcof.reslt && /bin/rm -f xmrqcof.o xmrqcof.out; \
	fi

xmrqmin.reslt : ../src/xmrqmin.c $(ILIB)
	-if [ -r ../bin/xmrqmin.out ]; then \
	../bin/xmrqmin.out < ../responses/xmrqmin.resp > xmrqmin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xmrqmin.out ../src/xmrqmin.c $(ILIB) -lm && \
	xmrqmin.out < ../responses/xmrqmin.resp > xmrqmin.reslt && \
	/bin/rm -f xmrqmin.o xmrqmin.out; \
	fi

xnewt.reslt : ../src/xnewt.c $(ILIB)
	-if [ -r ../bin/xnewt.out ]; then \
	../bin/xnewt.out > xnewt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xnewt.out ../src/xnewt.c $(ILIB) -lm && \
	xnewt.out > xnewt.reslt && /bin/rm -f xnewt.o xnewt.out; \
	fi

xodeint.reslt : ../src/xodeint.c $(ILIB)
	-if [ -r ../bin/xodeint.out ]; then \
	../bin/xodeint.out > xodeint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xodeint.out ../src/xodeint.c $(ILIB) -lm && \
	xodeint.out > xodeint.reslt && /bin/rm -f xodeint.o xodeint.out; \
	fi

xorthog.reslt : ../src/xorthog.c $(ILIB)
	-if [ -r ../bin/xorthog.out ]; then \
	../bin/xorthog.out < ../responses/xorthog.resp > xorthog.reslt; \
	else \
	$(CC) $(CFLAGS) -o xorthog.out ../src/xorthog.c $(ILIB) -lm && \
	xorthog.out < ../responses/xorthog.resp > xorthog.reslt && \
	/bin/rm -f xorthog.o xorthog.out; \
	fi

xpade.reslt : ../src/xpade.c $(ILIB)
	-if [ -r ../bin/xpade.out ]; then \
	../bin/xpade.out < ../responses/xpade.resp > xpade.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpade.out ../src/xpade.c $(ILIB) -lm && \
	xpade.out < ../responses/xpade.resp > xpade.reslt && \
	/bin/rm -f xpade.o xpade.out; \
	fi

xpccheb.reslt : ../src/xpccheb.c $(ILIB)
	-if [ -r ../bin/xpccheb.out ]; then \
	../bin/xpccheb.out > xpccheb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpccheb.out ../src/xpccheb.c $(ILIB) -lm && \
	xpccheb.out > xpccheb.reslt && /bin/rm -f xpccheb.o xpccheb.out; \
	fi

xpcshft.reslt : ../src/xpcshft.c $(ILIB)
	-if [ -r ../bin/xpcshft.out ]; then \
	../bin/xpcshft.out < ../responses/xpcshft.resp > xpcshft.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpcshft.out ../src/xpcshft.c $(ILIB) -lm && \
	xpcshft.out < ../responses/xpcshft.resp > xpcshft.reslt && \
	/bin/rm -f xpcshft.o xpcshft.out; \
	fi

xpearsn.reslt : ../src/xpearsn.c $(ILIB)
	-if [ -r ../bin/xpearsn.out ]; then \
	../bin/xpearsn.out > xpearsn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpearsn.out ../src/xpearsn.c $(ILIB) -lm && \
	xpearsn.out > xpearsn.reslt && /bin/rm -f xpearsn.o xpearsn.out; \
	fi

xperiod.reslt : ../src/xperiod.c $(ILIB)
	-if [ -r ../bin/xperiod.out ]; then \
	../bin/xperiod.out > xperiod.reslt; \
	else \
	$(CC) $(CFLAGS) -o xperiod.out ../src/xperiod.c $(ILIB) -lm && \
	xperiod.out > xperiod.reslt && /bin/rm -f xperiod.o xperiod.out; \
	fi

xpiksr2.reslt : ../src/xpiksr2.c $(ILIB)
	-if [ -r ../bin/xpiksr2.out ]; then \
	../bin/xpiksr2.out < ../responses/xpiksr2.resp > xpiksr2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpiksr2.out ../src/xpiksr2.c $(ILIB) -lm && \
	xpiksr2.out < ../responses/xpiksr2.resp > xpiksr2.reslt && \
	/bin/rm -f xpiksr2.o xpiksr2.out; \
	fi

xpiksrt.reslt : ../src/xpiksrt.c $(ILIB)
	-if [ -r ../bin/xpiksrt.out ]; then \
	../bin/xpiksrt.out > xpiksrt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpiksrt.out ../src/xpiksrt.c $(ILIB) -lm && \
	xpiksrt.out > xpiksrt.reslt && /bin/rm -f xpiksrt.o xpiksrt.out; \
	fi

xplgndr.reslt : ../src/xplgndr.c $(ILIB)
	-if [ -r ../bin/xplgndr.out ]; then \
	../bin/xplgndr.out > xplgndr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xplgndr.out ../src/xplgndr.c $(ILIB) -lm && \
	xplgndr.out > xplgndr.reslt && /bin/rm -f xplgndr.o xplgndr.out; \
	fi

xpoidev.reslt : ../src/xpoidev.c $(ILIB)
	-if [ -r ../bin/xpoidev.out ]; then \
	../bin/xpoidev.out < ../responses/xpoidev.resp > xpoidev.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpoidev.out ../src/xpoidev.c $(ILIB) -lm && \
	xpoidev.out < ../responses/xpoidev.resp > xpoidev.reslt && \
	/bin/rm -f xpoidev.o xpoidev.out; \
	fi

xpolcoe.reslt : ../src/xpolcoe.c $(ILIB)
	-if [ -r ../bin/xpolcoe.out ]; then \
	../bin/xpolcoe.out < ../responses/xpolcoe.resp > xpolcoe.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpolcoe.out ../src/xpolcoe.c $(ILIB) -lm && \
	xpolcoe.out < ../responses/xpolcoe.resp > xpolcoe.reslt && \
	/bin/rm -f xpolcoe.o xpolcoe.out; \
	fi

xpolcof.reslt : ../src/xpolcof.c $(ILIB)
	-if [ -r ../bin/xpolcof.out ]; then \
	../bin/xpolcof.out < ../responses/xpolcof.resp > xpolcof.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpolcof.out ../src/xpolcof.c $(ILIB) -lm && \
	xpolcof.out < ../responses/xpolcof.resp > xpolcof.reslt && \
	/bin/rm -f xpolcof.o xpolcof.out; \
	fi

xpoldiv.reslt : ../src/xpoldiv.c $(ILIB)
	-if [ -r ../bin/xpoldiv.out ]; then \
	../bin/xpoldiv.out > xpoldiv.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpoldiv.out ../src/xpoldiv.c $(ILIB) -lm && \
	xpoldiv.out > xpoldiv.reslt && /bin/rm -f xpoldiv.o xpoldiv.out; \
	fi

xpolin2.reslt : ../src/xpolin2.c $(ILIB)
	-if [ -r ../bin/xpolin2.out ]; then \
	../bin/xpolin2.out > xpolin2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpolin2.out ../src/xpolin2.c $(ILIB) -lm && \
	xpolin2.out > xpolin2.reslt && /bin/rm -f xpolin2.o xpolin2.out; \
	fi

xpolint.reslt : ../src/xpolint.c $(ILIB)
	-if [ -r ../bin/xpolint.out ]; then \
	../bin/xpolint.out < ../responses/xpolint.resp > xpolint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpolint.out ../src/xpolint.c $(ILIB) -lm && \
	xpolint.out < ../responses/xpolint.resp > xpolint.reslt && \
	/bin/rm -f xpolint.o xpolint.out; \
	fi

xpowell.reslt : ../src/xpowell.c $(ILIB)
	-if [ -r ../bin/xpowell.out ]; then \
	../bin/xpowell.out > xpowell.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpowell.out ../src/xpowell.c $(ILIB) -lm && \
	xpowell.out > xpowell.reslt && /bin/rm -f xpowell.o xpowell.out; \
	fi

xpredic.reslt : ../src/xpredic.c $(ILIB)
	-if [ -r ../bin/xpredic.out ]; then \
	../bin/xpredic.out > xpredic.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpredic.out ../src/xpredic.c $(ILIB) -lm && \
	xpredic.out > xpredic.reslt && /bin/rm -f xpredic.o xpredic.out; \
	fi

xprobks.reslt : ../src/xprobks.c $(ILIB)
	-if [ -r ../bin/xprobks.out ]; then \
	../bin/xprobks.out > xprobks.reslt; \
	else \
	$(CC) $(CFLAGS) -o xprobks.out ../src/xprobks.c $(ILIB) -lm && \
	xprobks.out > xprobks.reslt && /bin/rm -f xprobks.o xprobks.out; \
	fi

xpsdes.reslt : ../src/xpsdes.c $(ILIB)
	-if [ -r ../bin/xpsdes.out ]; then \
	../bin/xpsdes.out > xpsdes.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpsdes.out ../src/xpsdes.c $(ILIB) -lm && \
	xpsdes.out > xpsdes.reslt && /bin/rm -f xpsdes.o xpsdes.out; \
	fi

xpzextr.reslt : ../src/xpzextr.c $(ILIB)
	-if [ -r ../bin/xpzextr.out ]; then \
	../bin/xpzextr.out > xpzextr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xpzextr.out ../src/xpzextr.c $(ILIB) -lm && \
	xpzextr.out > xpzextr.reslt && /bin/rm -f xpzextr.o xpzextr.out; \
	fi

xqgaus.reslt : ../src/xqgaus.c $(ILIB)
	-if [ -r ../bin/xqgaus.out ]; then \
	../bin/xqgaus.out > xqgaus.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqgaus.out ../src/xqgaus.c $(ILIB) -lm && \
	xqgaus.out > xqgaus.reslt && /bin/rm -f xqgaus.o xqgaus.out; \
	fi

xqrdcmp.reslt : ../src/xqrdcmp.c $(ILIB)
	-if [ -r ../bin/xqrdcmp.out ]; then \
	../bin/xqrdcmp.out < ../responses/xqrdcmp.resp > xqrdcmp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqrdcmp.out ../src/xqrdcmp.c $(ILIB) -lm && \
	xqrdcmp.out < ../responses/xqrdcmp.resp > xqrdcmp.reslt && \
	/bin/rm -f xqrdcmp.o xqrdcmp.out; \
	fi

xqromb.reslt : ../src/xqromb.c $(ILIB)
	-if [ -r ../bin/xqromb.out ]; then \
	../bin/xqromb.out > xqromb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqromb.out ../src/xqromb.c $(ILIB) -lm && \
	xqromb.out > xqromb.reslt && /bin/rm -f xqromb.o xqromb.out; \
	fi

xqromo.reslt : ../src/xqromo.c $(ILIB)
	-if [ -r ../bin/xqromo.out ]; then \
	../bin/xqromo.out > xqromo.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqromo.out ../src/xqromo.c $(ILIB) -lm && \
	xqromo.out > xqromo.reslt && /bin/rm -f xqromo.o xqromo.out; \
	fi

xqroot.reslt : ../src/xqroot.c $(ILIB)
	-if [ -r ../bin/xqroot.out ]; then \
	../bin/xqroot.out > xqroot.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqroot.out ../src/xqroot.c $(ILIB) -lm && \
	xqroot.out > xqroot.reslt && /bin/rm -f xqroot.o xqroot.out; \
	fi

xqrsolv.reslt : ../src/xqrsolv.c $(ILIB)
	-if [ -r ../bin/xqrsolv.out ]; then \
	../bin/xqrsolv.out < ../responses/xqrsolv.resp > xqrsolv.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqrsolv.out ../src/xqrsolv.c $(ILIB) -lm && \
	xqrsolv.out < ../responses/xqrsolv.resp > xqrsolv.reslt && \
	/bin/rm -f xqrsolv.o xqrsolv.out; \
	fi

xqrupdt.reslt : ../src/xqrupdt.c $(ILIB)
	-if [ -r ../bin/xqrupdt.out ]; then \
	../bin/xqrupdt.out < ../responses/xqrupdt.resp > xqrupdt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqrupdt.out ../src/xqrupdt.c $(ILIB) -lm && \
	xqrupdt.out < ../responses/xqrupdt.resp > xqrupdt.reslt && \
	/bin/rm -f xqrupdt.o xqrupdt.out; \
	fi

xqsimp.reslt : ../src/xqsimp.c $(ILIB)
	-if [ -r ../bin/xqsimp.out ]; then \
	../bin/xqsimp.out > xqsimp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqsimp.out ../src/xqsimp.c $(ILIB) -lm && \
	xqsimp.out > xqsimp.reslt && /bin/rm -f xqsimp.o xqsimp.out; \
	fi

xqtrap.reslt : ../src/xqtrap.c $(ILIB)
	-if [ -r ../bin/xqtrap.out ]; then \
	../bin/xqtrap.out > xqtrap.reslt; \
	else \
	$(CC) $(CFLAGS) -o xqtrap.out ../src/xqtrap.c $(ILIB) -lm && \
	xqtrap.out > xqtrap.reslt && /bin/rm -f xqtrap.o xqtrap.out; \
	fi

xquad3d.reslt : ../src/xquad3d.c $(ILIB)
	-if [ -r ../bin/xquad3d.out ]; then \
	../bin/xquad3d.out > xquad3d.reslt; \
	else \
	$(CC) $(CFLAGS) -o xquad3d.out ../src/xquad3d.c $(ILIB) -lm && \
	xquad3d.out > xquad3d.reslt && /bin/rm -f xquad3d.o xquad3d.out; \
	fi

xran.reslt : ../src/xran.c $(ILIB)
	-if [ -r ../bin/xran.out ]; then \
	../bin/xran.out > xran.reslt; \
	else \
	$(CC) $(CFLAGS) -o xran.out ../src/xran.c $(ILIB) -lm && \
	xran.out > xran.reslt && /bin/rm -f xran.o xran.out; \
	fi

xran4.reslt : ../src/xran4.c $(ILIB)
	-if [ -r ../bin/xran4.out ]; then \
	../bin/xran4.out > xran4.reslt; \
	else \
	$(CC) $(CFLAGS) -o xran4.out ../src/xran4.c $(ILIB) -lm && \
	xran4.out > xran4.reslt && /bin/rm -f xran4.o xran4.out; \
	fi

xrank.reslt : ../src/xrank.c $(ILIB)
	-if [ -r ../bin/xrank.out ]; then \
	../bin/xrank.out < ../responses/xrank.resp > xrank.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrank.out ../src/xrank.c $(ILIB) -lm && \
	xrank.out < ../responses/xrank.resp > xrank.reslt && \
	/bin/rm -f xrank.o xrank.out; \
	fi

xratint.reslt : ../src/xratint.c $(ILIB)
	-if [ -r ../bin/xratint.out ]; then \
	../bin/xratint.out > xratint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xratint.out ../src/xratint.c $(ILIB) -lm && \
	xratint.out > xratint.reslt && /bin/rm -f xratint.o xratint.out; \
	fi

xratlsq.reslt : ../src/xratlsq.c $(ILIB)
	-if [ -r ../bin/xratlsq.out ]; then \
	../bin/xratlsq.out < ../responses/xratlsq.resp > xratlsq.reslt; \
	else \
	$(CC) $(CFLAGS) -o xratlsq.out ../src/xratlsq.c $(ILIB) -lm && \
	xratlsq.out < ../responses/xratlsq.resp > xratlsq.reslt && \
	/bin/rm -f xratlsq.o xratlsq.out; \
	fi

xrc.reslt : ../src/xrc.c $(ILIB)
	-if [ -r ../bin/xrc.out ]; then \
	../bin/xrc.out > xrc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrc.out ../src/xrc.c $(ILIB) -lm && \
	xrc.out > xrc.reslt && /bin/rm -f xrc.o xrc.out; \
	fi

xrd.reslt : ../src/xrd.c $(ILIB)
	-if [ -r ../bin/xrd.out ]; then \
	../bin/xrd.out > xrd.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrd.out ../src/xrd.c $(ILIB) -lm && \
	xrd.out > xrd.reslt && /bin/rm -f xrd.o xrd.out; \
	fi

xrealft.reslt : ../src/xrealft.c $(ILIB)
	-if [ -r ../bin/xrealft.out ]; then \
	../bin/xrealft.out < ../responses/xrealft.resp > xrealft.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrealft.out ../src/xrealft.c $(ILIB) -lm && \
	xrealft.out < ../responses/xrealft.resp > xrealft.reslt && \
	/bin/rm -f xrealft.o xrealft.out; \
	fi

xrf.reslt : ../src/xrf.c $(ILIB)
	-if [ -r ../bin/xrf.out ]; then \
	../bin/xrf.out > xrf.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrf.out ../src/xrf.c $(ILIB) -lm && \
	xrf.out > xrf.reslt && /bin/rm -f xrf.o xrf.out; \
	fi

xrj.reslt : ../src/xrj.c $(ILIB)
	-if [ -r ../bin/xrj.out ]; then \
	../bin/xrj.out > xrj.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrj.out ../src/xrj.c $(ILIB) -lm && \
	xrj.out > xrj.reslt && /bin/rm -f xrj.o xrj.out; \
	fi

xrk4.reslt : ../src/xrk4.c $(ILIB)
	-if [ -r ../bin/xrk4.out ]; then \
	../bin/xrk4.out > xrk4.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrk4.out ../src/xrk4.c $(ILIB) -lm && \
	xrk4.out > xrk4.reslt && /bin/rm -f xrk4.o xrk4.out; \
	fi

xrkdumb.reslt : ../src/xrkdumb.c $(ILIB)
	-if [ -r ../bin/xrkdumb.out ]; then \
	../bin/xrkdumb.out > xrkdumb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrkdumb.out ../src/xrkdumb.c $(ILIB) -lm && \
	xrkdumb.out > xrkdumb.reslt && /bin/rm -f xrkdumb.o xrkdumb.out; \
	fi

xrkqs.reslt : ../src/xrkqs.c $(ILIB)
	-if [ -r ../bin/xrkqs.out ]; then \
	../bin/xrkqs.out > xrkqs.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrkqs.out ../src/xrkqs.c $(ILIB) -lm && \
	xrkqs.out > xrkqs.reslt && /bin/rm -f xrkqs.o xrkqs.out; \
	fi

xrlft3.reslt : ../src/xrlft3.c $(ILIB)
	-if [ -r ../bin/xrlft3.out ]; then \
	../bin/xrlft3.out > xrlft3.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrlft3.out ../src/xrlft3.c $(ILIB) -lm && \
	xrlft3.out > xrlft3.reslt && /bin/rm -f xrlft3.o xrlft3.out; \
	fi

xrofunc.reslt : ../src/xrofunc.c $(ILIB)
	-if [ -r ../bin/xrofunc.out ]; then \
	../bin/xrofunc.out > xrofunc.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrofunc.out ../src/xrofunc.c $(ILIB) -lm && \
	xrofunc.out > xrofunc.reslt && /bin/rm -f xrofunc.o xrofunc.out; \
	fi

xrtbis.reslt : ../src/xrtbis.c $(ILIB)
	-if [ -r ../bin/xrtbis.out ]; then \
	../bin/xrtbis.out > xrtbis.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrtbis.out ../src/xrtbis.c $(ILIB) -lm && \
	xrtbis.out > xrtbis.reslt && /bin/rm -f xrtbis.o xrtbis.out; \
	fi

xrtflsp.reslt : ../src/xrtflsp.c $(ILIB)
	-if [ -r ../bin/xrtflsp.out ]; then \
	../bin/xrtflsp.out > xrtflsp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrtflsp.out ../src/xrtflsp.c $(ILIB) -lm && \
	xrtflsp.out > xrtflsp.reslt && /bin/rm -f xrtflsp.o xrtflsp.out; \
	fi

xrtnewt.reslt : ../src/xrtnewt.c $(ILIB)
	-if [ -r ../bin/xrtnewt.out ]; then \
	../bin/xrtnewt.out > xrtnewt.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrtnewt.out ../src/xrtnewt.c $(ILIB) -lm && \
	xrtnewt.out > xrtnewt.reslt && /bin/rm -f xrtnewt.o xrtnewt.out; \
	fi

xrtsafe.reslt : ../src/xrtsafe.c $(ILIB)
	-if [ -r ../bin/xrtsafe.out ]; then \
	../bin/xrtsafe.out > xrtsafe.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrtsafe.out ../src/xrtsafe.c $(ILIB) -lm && \
	xrtsafe.out > xrtsafe.reslt && /bin/rm -f xrtsafe.o xrtsafe.out; \
	fi

xrtsec.reslt : ../src/xrtsec.c $(ILIB)
	-if [ -r ../bin/xrtsec.out ]; then \
	../bin/xrtsec.out > xrtsec.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrtsec.out ../src/xrtsec.c $(ILIB) -lm && \
	xrtsec.out > xrtsec.reslt && /bin/rm -f xrtsec.o xrtsec.out; \
	fi

xrzextr.reslt : ../src/xrzextr.c $(ILIB)
	-if [ -r ../bin/xrzextr.out ]; then \
	../bin/xrzextr.out > xrzextr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xrzextr.out ../src/xrzextr.c $(ILIB) -lm && \
	xrzextr.out > xrzextr.reslt && /bin/rm -f xrzextr.o xrzextr.out; \
	fi

xsavgol.reslt : ../src/xsavgol.c $(ILIB)
	-if [ -r ../bin/xsavgol.out ]; then \
	../bin/xsavgol.out > xsavgol.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsavgol.out ../src/xsavgol.c $(ILIB) -lm && \
	xsavgol.out > xsavgol.reslt && /bin/rm -f xsavgol.o xsavgol.out; \
	fi

xscrsho.reslt : ../src/xscrsho.c $(ILIB)
	-if [ -r ../bin/xscrsho.out ]; then \
	../bin/xscrsho.out < ../responses/xscrsho.resp > xscrsho.reslt; \
	else \
	$(CC) $(CFLAGS) -o xscrsho.out ../src/xscrsho.c $(ILIB) -lm && \
	xscrsho.out < ../responses/xscrsho.resp > xscrsho.reslt && \
	/bin/rm -f xscrsho.o xscrsho.out; \
	fi

xselect.reslt : ../src/xselect.c $(ILIB)
	-if [ -r ../bin/xselect.out ]; then \
	../bin/xselect.out < ../responses/xselect.resp > xselect.reslt; \
	else \
	$(CC) $(CFLAGS) -o xselect.out ../src/xselect.c $(ILIB) -lm && \
	xselect.out < ../responses/xselect.resp > xselect.reslt && \
	/bin/rm -f xselect.o xselect.out; \
	fi

xselip.reslt : ../src/xselip.c $(ILIB)
	-if [ -r ../bin/xselip.out ]; then \
	../bin/xselip.out < ../responses/xselip.resp > xselip.reslt; \
	else \
	$(CC) $(CFLAGS) -o xselip.out ../src/xselip.c $(ILIB) -lm && \
	xselip.out < ../responses/xselip.resp > xselip.reslt && \
	/bin/rm -f xselip.o xselip.out; \
	fi

xshell.reslt : ../src/xshell.c $(ILIB)
	-if [ -r ../bin/xshell.out ]; then \
	../bin/xshell.out > xshell.reslt; \
	else \
	$(CC) $(CFLAGS) -o xshell.out ../src/xshell.c $(ILIB) -lm && \
	xshell.out > xshell.reslt && /bin/rm -f xshell.o xshell.out; \
	fi

xsimplx.reslt : ../src/xsimplx.c $(ILIB)
	-if [ -r ../bin/xsimplx.out ]; then \
	../bin/xsimplx.out > xsimplx.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsimplx.out ../src/xsimplx.c $(ILIB) -lm && \
	xsimplx.out > xsimplx.reslt && /bin/rm -f xsimplx.o xsimplx.out; \
	fi

xsimpr.reslt : ../src/xsimpr.c $(ILIB)
	-if [ -r ../bin/xsimpr.out ]; then \
	../bin/xsimpr.out > xsimpr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsimpr.out ../src/xsimpr.c $(ILIB) -lm && \
	xsimpr.out > xsimpr.reslt && /bin/rm -f xsimpr.o xsimpr.out; \
	fi

xsinft.reslt : ../src/xsinft.c $(ILIB)
	-if [ -r ../bin/xsinft.out ]; then \
	../bin/xsinft.out < ../responses/xsinft.resp > xsinft.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsinft.out ../src/xsinft.c $(ILIB) -lm && \
	xsinft.out < ../responses/xsinft.resp > xsinft.reslt && \
	/bin/rm -f xsinft.o xsinft.out; \
	fi

xsncndn.reslt : ../src/xsncndn.c $(ILIB)
	-if [ -r ../bin/xsncndn.out ]; then \
	../bin/xsncndn.out > xsncndn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsncndn.out ../src/xsncndn.c $(ILIB) -lm && \
	xsncndn.out > xsncndn.reslt && /bin/rm -f xsncndn.o xsncndn.out; \
	fi

xsobseq.reslt : ../src/xsobseq.c $(ILIB)
	-if [ -r ../bin/xsobseq.out ]; then \
	../bin/xsobseq.out > xsobseq.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsobseq.out ../src/xsobseq.c $(ILIB) -lm && \
	xsobseq.out > xsobseq.reslt && /bin/rm -f xsobseq.o xsobseq.out; \
	fi

xsor.reslt : ../src/xsor.c $(ILIB)
	-if [ -r ../bin/xsor.out ]; then \
	../bin/xsor.out > xsor.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsor.out ../src/xsor.c $(ILIB) -lm && \
	xsor.out > xsor.reslt && /bin/rm -f xsor.o xsor.out; \
	fi

xsort.reslt : ../src/xsort.c $(ILIB)
	-if [ -r ../bin/xsort.out ]; then \
	../bin/xsort.out > xsort.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsort.out ../src/xsort.c $(ILIB) -lm && \
	xsort.out > xsort.reslt && /bin/rm -f xsort.o xsort.out; \
	fi

xsort2.reslt : ../src/xsort2.c $(ILIB)
	-if [ -r ../bin/xsort2.out ]; then \
	../bin/xsort2.out < ../responses/xsort2.resp > xsort2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsort2.out ../src/xsort2.c $(ILIB) -lm && \
	xsort2.out < ../responses/xsort2.resp > xsort2.reslt && \
	/bin/rm -f xsort2.o xsort2.out; \
	fi

xsort3.reslt : ../src/xsort3.c $(ILIB)
	-if [ -r ../bin/xsort3.out ]; then \
	../bin/xsort3.out > xsort3.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsort3.out ../src/xsort3.c $(ILIB) -lm && \
	xsort3.out > xsort3.reslt && /bin/rm -f xsort3.o xsort3.out; \
	fi

xspctrm.reslt : ../src/xspctrm.c $(ILIB)
	-if [ -r ../bin/xspctrm.out ]; then \
	../bin/xspctrm.out > xspctrm.reslt; \
	else \
	$(CC) $(CFLAGS) -o xspctrm.out ../src/xspctrm.c $(ILIB) -lm && \
	xspctrm.out > xspctrm.reslt && /bin/rm -f xspctrm.o xspctrm.out; \
	fi

xspear.reslt : ../src/xspear.c $(ILIB)
	-if [ -r ../bin/xspear.out ]; then \
	../bin/xspear.out > xspear.reslt; \
	else \
	$(CC) $(CFLAGS) -o xspear.out ../src/xspear.c $(ILIB) -lm && \
	xspear.out > xspear.reslt && /bin/rm -f xspear.o xspear.out; \
	fi

xsphbes.reslt : ../src/xsphbes.c $(ILIB)
	-if [ -r ../bin/xsphbes.out ]; then \
	../bin/xsphbes.out > xsphbes.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsphbes.out ../src/xsphbes.c $(ILIB) -lm && \
	xsphbes.out > xsphbes.reslt && /bin/rm -f xsphbes.o xsphbes.out; \
	fi

xsplie2.reslt : ../src/xsplie2.c $(ILIB)
	-if [ -r ../bin/xsplie2.out ]; then \
	../bin/xsplie2.out > xsplie2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsplie2.out ../src/xsplie2.c $(ILIB) -lm && \
	xsplie2.out > xsplie2.reslt && /bin/rm -f xsplie2.o xsplie2.out; \
	fi

xsplin2.reslt : ../src/xsplin2.c $(ILIB)
	-if [ -r ../bin/xsplin2.out ]; then \
	../bin/xsplin2.out > xsplin2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsplin2.out ../src/xsplin2.c $(ILIB) -lm && \
	xsplin2.out > xsplin2.reslt && /bin/rm -f xsplin2.o xsplin2.out; \
	fi

xspline.reslt : ../src/xspline.c $(ILIB)
	-if [ -r ../bin/xspline.out ]; then \
	../bin/xspline.out > xspline.reslt; \
	else \
	$(CC) $(CFLAGS) -o xspline.out ../src/xspline.c $(ILIB) -lm && \
	xspline.out > xspline.reslt && /bin/rm -f xspline.o xspline.out; \
	fi

xsplint.reslt : ../src/xsplint.c $(ILIB)
	-if [ -r ../bin/xsplint.out ]; then \
	../bin/xsplint.out < ../responses/xsplint.resp > xsplint.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsplint.out ../src/xsplint.c $(ILIB) -lm && \
	xsplint.out < ../responses/xsplint.resp > xsplint.reslt && \
	/bin/rm -f xsplint.o xsplint.out; \
	fi

xsprsax.reslt : ../src/xsprsax.c $(ILIB)
	-if [ -r ../bin/xsprsax.out ]; then \
	../bin/xsprsax.out > xsprsax.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprsax.out ../src/xsprsax.c $(ILIB) -lm && \
	xsprsax.out > xsprsax.reslt && /bin/rm -f xsprsax.o xsprsax.out; \
	fi

xsprsin.reslt : ../src/xsprsin.c $(ILIB)
	-if [ -r ../bin/xsprsin.out ]; then \
	../bin/xsprsin.out > xsprsin.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprsin.out ../src/xsprsin.c $(ILIB) -lm && \
	xsprsin.out > xsprsin.reslt && /bin/rm -f xsprsin.o xsprsin.out; \
	fi

xsprspm.reslt : ../src/xsprspm.c $(ILIB)
	-if [ -r ../bin/xsprspm.out ]; then \
	../bin/xsprspm.out > xsprspm.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprspm.out ../src/xsprspm.c $(ILIB) -lm && \
	xsprspm.out > xsprspm.reslt && /bin/rm -f xsprspm.o xsprspm.out; \
	fi

xsprstm.reslt : ../src/xsprstm.c $(ILIB)
	-if [ -r ../bin/xsprstm.out ]; then \
	../bin/xsprstm.out > xsprstm.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprstm.out ../src/xsprstm.c $(ILIB) -lm && \
	xsprstm.out > xsprstm.reslt && /bin/rm -f xsprstm.o xsprstm.out; \
	fi

xsprstp.reslt : ../src/xsprstp.c $(ILIB)
	-if [ -r ../bin/xsprstp.out ]; then \
	../bin/xsprstp.out > xsprstp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprstp.out ../src/xsprstp.c $(ILIB) -lm && \
	xsprstp.out > xsprstp.reslt && /bin/rm -f xsprstp.o xsprstp.out; \
	fi

xsprstx.reslt : ../src/xsprstx.c $(ILIB)
	-if [ -r ../bin/xsprstx.out ]; then \
	../bin/xsprstx.out > xsprstx.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsprstx.out ../src/xsprstx.c $(ILIB) -lm && \
	xsprstx.out > xsprstx.reslt && /bin/rm -f xsprstx.o xsprstx.out; \
	fi

xstifbs.reslt : ../src/xstifbs.c $(ILIB)
	-if [ -r ../bin/xstifbs.out ]; then \
	../bin/xstifbs.out < ../responses/xstifbs.resp > xstifbs.reslt; \
	else \
	$(CC) $(CFLAGS) -o xstifbs.out ../src/xstifbs.c $(ILIB) -lm && \
	xstifbs.out < ../responses/xstifbs.resp > xstifbs.reslt && \
	/bin/rm -f xstifbs.o xstifbs.out; \
	fi

xstiff.reslt : ../src/xstiff.c $(ILIB)
	-if [ -r ../bin/xstiff.out ]; then \
	../bin/xstiff.out < ../responses/xstiff.resp > xstiff.reslt; \
	else \
	$(CC) $(CFLAGS) -o xstiff.out ../src/xstiff.c $(ILIB) -lm && \
	xstiff.out < ../responses/xstiff.resp > xstiff.reslt && \
	/bin/rm -f xstiff.o xstiff.out; \
	fi

xstoerm.reslt : ../src/xstoerm.c $(ILIB)
	-if [ -r ../bin/xstoerm.out ]; then \
	../bin/xstoerm.out > xstoerm.reslt; \
	else \
	$(CC) $(CFLAGS) -o xstoerm.out ../src/xstoerm.c $(ILIB) -lm && \
	xstoerm.out > xstoerm.reslt && /bin/rm -f xstoerm.o xstoerm.out; \
	fi

xsvbksb.reslt : ../src/xsvbksb.c $(ILIB)
	-if [ -r ../bin/xsvbksb.out ]; then \
	../bin/xsvbksb.out < ../responses/xsvbksb.resp > xsvbksb.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsvbksb.out ../src/xsvbksb.c $(ILIB) -lm && \
	xsvbksb.out < ../responses/xsvbksb.resp > xsvbksb.reslt && \
	/bin/rm -f xsvbksb.o xsvbksb.out; \
	fi

xsvdcmp.reslt : ../src/xsvdcmp.c $(ILIB)
	-if [ -r ../bin/xsvdcmp.out ]; then \
	../bin/xsvdcmp.out < ../responses/xsvdcmp.resp > xsvdcmp.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsvdcmp.out ../src/xsvdcmp.c $(ILIB) -lm && \
	xsvdcmp.out < ../responses/xsvdcmp.resp > xsvdcmp.reslt && \
	/bin/rm -f xsvdcmp.o xsvdcmp.out; \
	fi

xsvdfit.reslt : ../src/xsvdfit.c $(ILIB)
	-if [ -r ../bin/xsvdfit.out ]; then \
	../bin/xsvdfit.out > xsvdfit.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsvdfit.out ../src/xsvdfit.c $(ILIB) -lm && \
	xsvdfit.out > xsvdfit.reslt && /bin/rm -f xsvdfit.o xsvdfit.out; \
	fi

xsvdvar.reslt : ../src/xsvdvar.c $(ILIB)
	-if [ -r ../bin/xsvdvar.out ]; then \
	../bin/xsvdvar.out > xsvdvar.reslt; \
	else \
	$(CC) $(CFLAGS) -o xsvdvar.out ../src/xsvdvar.c $(ILIB) -lm && \
	xsvdvar.out > xsvdvar.reslt && /bin/rm -f xsvdvar.o xsvdvar.out; \
	fi

xtoeplz.reslt : ../src/xtoeplz.c $(ILIB)
	-if [ -r ../bin/xtoeplz.out ]; then \
	../bin/xtoeplz.out > xtoeplz.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtoeplz.out ../src/xtoeplz.c $(ILIB) -lm && \
	xtoeplz.out > xtoeplz.reslt && /bin/rm -f xtoeplz.o xtoeplz.out; \
	fi

xtptest.reslt : ../src/xtptest.c $(ILIB)
	-if [ -r ../bin/xtptest.out ]; then \
	../bin/xtptest.out > xtptest.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtptest.out ../src/xtptest.c $(ILIB) -lm && \
	xtptest.out > xtptest.reslt && /bin/rm -f xtptest.o xtptest.out; \
	fi

xtqli.reslt : ../src/xtqli.c $(ILIB)
	-if [ -r ../bin/xtqli.out ]; then \
	../bin/xtqli.out < ../responses/xtqli.resp > xtqli.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtqli.out ../src/xtqli.c $(ILIB) -lm && \
	xtqli.out < ../responses/xtqli.resp > xtqli.reslt && \
	/bin/rm -f xtqli.o xtqli.out; \
	fi

xtrapzd.reslt : ../src/xtrapzd.c $(ILIB)
	-if [ -r ../bin/xtrapzd.out ]; then \
	../bin/xtrapzd.out > xtrapzd.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtrapzd.out ../src/xtrapzd.c $(ILIB) -lm && \
	xtrapzd.out > xtrapzd.reslt && /bin/rm -f xtrapzd.o xtrapzd.out; \
	fi

xtred2.reslt : ../src/xtred2.c $(ILIB)
	-if [ -r ../bin/xtred2.out ]; then \
	../bin/xtred2.out > xtred2.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtred2.out ../src/xtred2.c $(ILIB) -lm && \
	xtred2.out > xtred2.reslt && /bin/rm -f xtred2.o xtred2.out; \
	fi

xtridag.reslt : ../src/xtridag.c $(ILIB)
	-if [ -r ../bin/xtridag.out ]; then \
	../bin/xtridag.out < ../responses/xtridag.resp > xtridag.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtridag.out ../src/xtridag.c $(ILIB) -lm && \
	xtridag.out < ../responses/xtridag.resp > xtridag.reslt && \
	/bin/rm -f xtridag.o xtridag.out; \
	fi

xttest.reslt : ../src/xttest.c $(ILIB)
	-if [ -r ../bin/xttest.out ]; then \
	../bin/xttest.out > xttest.reslt; \
	else \
	$(CC) $(CFLAGS) -o xttest.out ../src/xttest.c $(ILIB) -lm && \
	xttest.out > xttest.reslt && /bin/rm -f xttest.o xttest.out; \
	fi

xtutest.reslt : ../src/xtutest.c $(ILIB)
	-if [ -r ../bin/xtutest.out ]; then \
	../bin/xtutest.out > xtutest.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtutest.out ../src/xtutest.c $(ILIB) -lm && \
	xtutest.out > xtutest.reslt && /bin/rm -f xtutest.o xtutest.out; \
	fi

xtwofft.reslt : ../src/xtwofft.c $(ILIB)
	-if [ -r ../bin/xtwofft.out ]; then \
	../bin/xtwofft.out < ../responses/xtwofft.resp > xtwofft.reslt; \
	else \
	$(CC) $(CFLAGS) -o xtwofft.out ../src/xtwofft.c $(ILIB) -lm && \
	xtwofft.out < ../responses/xtwofft.resp > xtwofft.reslt && \
	/bin/rm -f xtwofft.o xtwofft.out; \
	fi

xvander.reslt : ../src/xvander.c $(ILIB)
	-if [ -r ../bin/xvander.out ]; then \
	../bin/xvander.out > xvander.reslt; \
	else \
	$(CC) $(CFLAGS) -o xvander.out ../src/xvander.c $(ILIB) -lm && \
	xvander.out > xvander.reslt && /bin/rm -f xvander.o xvander.out; \
	fi

xvegas.reslt : ../src/xvegas.c $(ILIB)
	-if [ -r ../bin/xvegas.out ]; then \
	../bin/xvegas.out < ../responses/xvegas.resp > xvegas.reslt; \
	else \
	$(CC) $(CFLAGS) -o xvegas.out ../src/xvegas.c $(ILIB) -lm && \
	xvegas.out < ../responses/xvegas.resp > xvegas.reslt && \
	/bin/rm -f xvegas.o xvegas.out; \
	fi

xvoltra.reslt : ../src/xvoltra.c $(ILIB)
	-if [ -r ../bin/xvoltra.out ]; then \
	../bin/xvoltra.out > xvoltra.reslt; \
	else \
	$(CC) $(CFLAGS) -o xvoltra.out ../src/xvoltra.c $(ILIB) -lm && \
	xvoltra.out > xvoltra.reslt && /bin/rm -f xvoltra.o xvoltra.out; \
	fi

xwt1.reslt : ../src/xwt1.c $(ILIB)
	-if [ -r ../bin/xwt1.out ]; then \
	../bin/xwt1.out < ../responses/xwt1.resp > xwt1.reslt; \
	else \
	$(CC) $(CFLAGS) -o xwt1.out ../src/xwt1.c $(ILIB) -lm && \
	xwt1.out < ../responses/xwt1.resp > xwt1.reslt && \
	/bin/rm -f xwt1.o xwt1.out; \
	fi

xwtn.reslt : ../src/xwtn.c $(ILIB)
	-if [ -r ../bin/xwtn.out ]; then \
	../bin/xwtn.out > xwtn.reslt; \
	else \
	$(CC) $(CFLAGS) -o xwtn.out ../src/xwtn.c $(ILIB) -lm && \
	xwtn.out > xwtn.reslt && /bin/rm -f xwtn.o xwtn.out; \
	fi

xzbrac.reslt : ../src/xzbrac.c $(ILIB)
	-if [ -r ../bin/xzbrac.out ]; then \
	../bin/xzbrac.out > xzbrac.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzbrac.out ../src/xzbrac.c $(ILIB) -lm && \
	xzbrac.out > xzbrac.reslt && /bin/rm -f xzbrac.o xzbrac.out; \
	fi

xzbrak.reslt : ../src/xzbrak.c $(ILIB)
	-if [ -r ../bin/xzbrak.out ]; then \
	../bin/xzbrak.out > xzbrak.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzbrak.out ../src/xzbrak.c $(ILIB) -lm && \
	xzbrak.out > xzbrak.reslt && /bin/rm -f xzbrak.o xzbrak.out; \
	fi

xzbrent.reslt : ../src/xzbrent.c $(ILIB)
	-if [ -r ../bin/xzbrent.out ]; then \
	../bin/xzbrent.out > xzbrent.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzbrent.out ../src/xzbrent.c $(ILIB) -lm && \
	xzbrent.out > xzbrent.reslt && /bin/rm -f xzbrent.o xzbrent.out; \
	fi

xzrhqr.reslt : ../src/xzrhqr.c $(ILIB)
	-if [ -r ../bin/xzrhqr.out ]; then \
	../bin/xzrhqr.out > xzrhqr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzrhqr.out ../src/xzrhqr.c $(ILIB) -lm && \
	xzrhqr.out > xzrhqr.reslt && /bin/rm -f xzrhqr.o xzrhqr.out; \
	fi

xzriddr.reslt : ../src/xzriddr.c $(ILIB)
	-if [ -r ../bin/xzriddr.out ]; then \
	../bin/xzriddr.out > xzriddr.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzriddr.out ../src/xzriddr.c $(ILIB) -lm && \
	xzriddr.out > xzriddr.reslt && /bin/rm -f xzriddr.o xzriddr.out; \
	fi

xzroots.reslt : ../src/xzroots.c $(ILIB)
	-if [ -r ../bin/xzroots.out ]; then \
	../bin/xzroots.out > xzroots.reslt; \
	else \
	$(CC) $(CFLAGS) -o xzroots.out ../src/xzroots.c $(ILIB) -lm && \
	xzroots.out > xzroots.reslt && /bin/rm -f xzroots.o xzroots.out; \
	fi

badluk.reslt : ../../recipes/badluk.c $(ILIB)
	-if [ -r $(BINDIR)/badluk ]; then \
	$(BINDIR)/badluk > badluk.reslt; \
	else \
	$(CC) $(CFLAGS) -o badluk ../../recipes/badluk.c $(ILIB) -lm && \
	badluk > badluk.reslt && /bin/rm -f badluk.o badluk; \
	fi

fredex.reslt : ../../recipes/fredex.c $(ILIB)
	-if [ -r $(BINDIR)/fredex ]; then \
	$(BINDIR)/fredex > fredex.reslt; \
	else \
	$(CC) $(CFLAGS) -o fredex ../../recipes/fredex.c $(ILIB) -lm && \
	fredex > fredex.reslt && /bin/rm -f fredex.o fredex; \
	fi

sfroid.reslt : ../../recipes/sfroid.c $(ILIB)
	-if [ -r $(BINDIR)/sfroid ]; then \
	$(BINDIR)/sfroid < ../responses/sfroid.resp > sfroid.reslt; \
	else \
	$(CC) $(CFLAGS) -o sfroid ../../recipes/sfroid.c $(ILIB) -lm && \
	sfroid < ../responses/sfroid.resp > sfroid.reslt && \
	/bin/rm -f sfroid.o sfroid; \
	fi

sphoot.reslt : ../../recipes/sphoot.c $(ILIB)
	-if [ -r $(BINDIR)/sphoot ]; then \
	$(BINDIR)/sphoot < ../responses/sphoot.resp > sphoot.reslt; \
	else \
	$(CC) $(CFLAGS) -o sphoot ../../recipes/sphoot.c $(ILIB) -lm && \
	sphoot < ../responses/sphoot.resp > sphoot.reslt && \
	/bin/rm -f sphoot.o sphoot; \
	fi

xsphfpt.reslt : ../../recipes/sphfpt.c ../src/xsphfpt.c $(ILIB)
	-if [ -r $(BINDIR)/sphfpt ]; then \
	$(BINDIR)/sphfpt < ../responses/xsphfpt.resp > xsphfpt.reslt; \
	else \
	$(CC) $(CFLAGS) -o sphfpt ../src/xsphfpt.c \
	../../recipes/sphfpt.c ../../recipes/shootf.c $(ILIB) -lm && \
	sphfpt < ../responses/xsphfpt.resp > xsphfpt.reslt && \
	/bin/rm -f xsphfpt.o sphfpt.o shootf.o sphfpt; \
	fi
