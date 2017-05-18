Ccomdeck info.inc $Revision: 6.2 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Parameter (iTabMx=15)
       Parameter (maxbfn=5000)
       Parameter (mxAtom = 500)
       Parameter (mxroot = 100)
       Parameter (mxact = 50)
       Parameter (mxina = 100)
       Parameter (mxbas = maxbfn)
       Parameter (mxOrb = maxbfn)
       Parameter (mxSym = 8)

      Parameter ( Mxdbsc=MxAtom, MxShll=iTabMx*MxAtom, Mxdc=MxAtom,
     & MxAng=iTabMx, MxUnq=maxbfn, MxAO=maxbfn,
     & MxFnc=(MxAng+1)*(MxAng+2)*(MxAng+3)/6)
      Integer nExp(MxShll), nBasis(MxShll), ipCff(MxShll),
     & ipExp(MxShll), ipCntr(Mxdbsc), nCntr(Mxdbsc),
     & ipVal(Mxdbsc), ipPrj(Mxdbsc), ipSRO(Mxdbsc),
     & ipSOC(Mxdbsc),ipPP(Mxdbsc),
     & nTot_Shells(Mxdbsc), nVal_Shells(mxdbsc),
     & nPrj_Shells(Mxdbsc), nSRO_Shells(mxdbsc),
     & nSOC_Shells(Mxdbsc), nPP_Shells(Mxdbsc),iChCnt(Mxdc),
     & iChTbl(0:7,0:7), jStab(0:7,Mxdc), nStab(Mxdc),
     & iCoSet(0:7,0:7,Mxdc), iSkip(8), iPhase(3,0:7),
     & iAngMx, iOper(0:7), nBas(0:7),
     & MaxPrm(0:iTabMx), MaxBas(0:iTabMx), IrrCmp(MxUnq),
     & iChBas(MxFnc), iSphCr(MxFnc),
     & IndS(MxShll), ip_Occ(Mxshll),
     & iAOtSO(MxAO,0:7),lOffAO(Mxdbsc), NrInt(-20:9),
     & kOffAO(Mxdbsc,0:iTabMx), nPrm(0:7),
     & iOffSO(0:7), ipAkl(Mxshll),
     & iSquar,iPack,iWRopt,lwInt,lwSyB,lwSqN,lwIntL,NotZer,
     & ipM1xp(Mxdbsc),ipM1cf(Mxdbsc),nM1(Mxdbsc),
     & ipM2xp(Mxdbsc),ipM2cf(Mxdbsc),nM2(Mxdbsc),ipBk(Mxshll),
     & nOpt(MxShll), iAtmNr(mxdbsc),
     & ipPAM2xp(Mxdbsc),ipPAM2cf(Mxdbsc),nPAM2(Mxdbsc),
     & ipPAMexp, kCnttpPAM,iPAMPrim,iPAMcount,
     & iOrdFm, nOrd_XF,iXPolType,nXMolnr
      Common /InfoPAM/ ipPAM2xp,ipPAM2cf,nPAM2,ipPAMexp,
     & kCnttpPAM,iPAMPrim,iPAMcount
      Common /Info/
     & nExp, nBasis, ipCff, ipExp, IndS, ip_Occ, ipAkl,
     & ipBk, nOpt, iCoSet,
     & ixStrt, ! Put never ever anything before ixStrt!!!
     & Mx_mdc, Mx_Shll,ipCntr, nCntr, nOrd_XF,iXPolType,
     & iChCnt, iChTbl, jStab, nStab, iPhase,
     & iAngMx, nCnttp, nIrrep, iOper, nBas, n2Tot,
     & n2Max, MaxDCR, MaxPrm, MaxBas, IrrCmp, nDim,
     & iChBas, iSphCr, ipSOff, iAOtSO, kOffAO,
     & lOffAO, mCentr, nPrm, m2Max, m2Tot, IntTot,
     & NrInt, MemHid, iOffSO, ipMpMC, Max_Cnt,
     & nMltpl, iMass, iCentr, iChrg, iLbls, ipEF,
     & nEF, ipOAM, ipEFg, nEFg, ipAMP, ipDMS,
     & nDMS, ipWel, nWel, nShlls, iPack,
     & iSquar, iWRopt, iSkip, lwInt, lwSyB, lwSqN,
     & lwIntL, NotZer, LctInf, ipM1xp,
     & ipM1cf, nM1, ipM2xp, ipM2cf, nM2,
     & LenInf, LenSph, MolWgh, nXF, ipXF, iAtmNr,
     & iShll_RF,ipVal, ipPrj, ipSRO, ipSOC,nTot_Shells,
     & nVal_Shells, nPrj_Shells, nSRO_Shells,
     & nSOC_Shells, nXMolnr, ipXMolnr, ipXEle,
     & iOrdFm, ipPP, nPP_Shells,
     & ixEnd ! Put never ever anything after ixEnd!!!
      Logical Prjct(MxShll),Transf(MxShll),DirInt, Onenly, Petite,
     & lSOInt, Test, Vlct, lEF, lOAM, lAMP,
     & lEFg, Prprt, Short, lDMS, lRel, UnNorm, lSchw,
     & Dist, lAMFI, lfinite, force_part_p,
     & lECP, lNOpair, NoPairL(Mxdbsc), force_part_c, Expert,
     & ECP(Mxdbsc), pChrg(Mxdbsc), Fixed(Mxdbsc),lXF, BSS,
     & DKroll, Primitive_Pass, First_Pass, NEMO,lPP,
     & PAM2(Mxdbsc), lPAM2, GIAO, ExCoor,
     & SODK(Mxdbsc),
     & IsMM(Mxdbsc), Cholesky
      Integer lxStrt, lxEnd
      Common /LInfoPAM/ lPAM2,PAM2
      Common /LInfo/ Prjct, Transf,
     & lxStrt, ! Put never ever anything before lxStrt!!!
     & lfinite, DirInt, DKroll, BSS, Primitive_Pass,
     & First_Pass, NEMO,lPP,
     & lAMFI, Onenly, Petite, lSOInt,
     & Test, Vlct, lEF, lOAM, lAMP,
     & lEFg, Prprt, Short, lDMS, lRel, UnNorm,
     & lSchw, Dist,ECP, lECP, lXF,pChrg, Fixed, lNOpair,
     & NoPairL, force_part_c, force_part_p, Expert,
     & GIAO, SODK, ExCoor, IsMM,
     & Cholesky,
     & lxEnd ! Put never ever anything after lxEnd!!!
      Real*8 AccMch, ThrInt, Charge(mxdbsc),
     & PotNuc, Sum, ExpNuc(mxdbsc),
     & SumSq, SumAbs, CutInt, TMass, CoM(3), rMI(6),
     & Prin(3),PAx(9), CoC(3), qNuc, PkAcc, Thrs,
     & Dxyz(3), RadMax , rChTbl(0:7,0:7), cdMax, EtMax,
     & CrRep(mxdbsc)
      Real*8 rxStrt, rxEnd
      Common /RInfo/ rxStrt, ! Put never ever anything before rxStrt!!!
     & AccMch, ThrInt, Charge, PotNuc, Sum,
     & SumSq, SumAbs, CutInt, TMass, CoM, rMI,
     & Prin, PAx, CoC, qNuc, PkAcc, Thrs,
     & Dxyz, RadMax, cdMax, EtMax, rChTbl,CrRep,
     & ExpNuc,
     & rxEnd ! Put never ever anything after rxEnd!!!
      Character BLine*80, lIrrep(0:7)*3, lBsFnc(0:7)*80,
     & LblCnt(Mxdc)*8, lblxxx(8*mxdc)*8
      Character LblCBs(MxFnc)*4, LblSBs(MxFnc)*4
      Integer cxStrt,cxEnd
      Common /CInfo/ cxStrt, ! Put never ever anything before cxStrt!!!
     & LblCBs, LblSBs,
     & BLine, lIrrep, lBsFnc, LblCnt, lblxxx,
     & cxEnd ! Put never ever anything after cxEnd!!!

      Real*8    Work(8)
      Integer   iWork(8)
      Integer   mxMem
      Equivalence (Work,iWork)
      Common /WrkSpc/Work,mxMem
