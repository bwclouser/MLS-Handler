;+
; CLASS_NAME:
;       MLSData
;
; PURPOSE:
;       An MLSData object provides a simple interface with which to
;       load, process, and display data from the MLS Aura experiment.
;
; CATEGORY:
;       Data handling and plotting.
;
; SUPERCLASSES:
;       This class inherits from no other classes.
;
; SUBCLASSES:
;       This class has no subclasses.
;
; CREATION:
;       See MLSData::Init
;
; METHODS:
;       Intrinsic Methods
;       This class has the following methods:
;
;
; USAGE:
;   These programs are designed for use with MLS retrievals. There is no guarantee they
;   will work properly with other versions of the data. Ideally, the user will only need to interact
;   with the following methods to generate plots and analyze data:
;
;
;
; DEPENDENCIES:
;
;       The readl2gp routine written in IDL and available from ??
;
;
; MODIFICATION HISTORY:
;-

FUNCTION MLSData::init

  self.raw=PTR_NEW(/ALLOCATE)
  RETURN,1

end

PRO MLSData::Cleanup

  IF PTR_VALID(self.raw) then PTR_FREE,self.raw
  IF PTR_VALID((*self.flags).tflag) then PTR_FREE,(*self.flags).tflag
  IF PTR_VALID((*self.flags).sflag) then PTR_FREE,(*self.flags).sflag
  IF PTR_VALID((*self.flags).eflag) then PTR_FREE,(*self.flags).eflag
  IF PTR_VALID((*self.flags).aflag) then PTR_FREE,(*self.flags).aflag
  IF PTR_VALID(self.flags) THEN PTR_FREE,self.flags
  IF PTR_VALID(self.opMask) THEN PTR_FREE,self.opMask
  ;IF PTR_VALID((*self.blocks).orbIDs) THEN PTR_FREE,(*self.blocks).orbIDs
  IF PTR_VALID(self.blocks) THEN PTR_FREE,self.blocks

END

PRO MLSData::loadrange,fpath,molecule,yearlims=yearlims,monthlims=monthlims,daylims=daylims,hourlims=hourlims,timearr=timearr

  CD,fpath
  mol=molecule
  molecule='*'+molecule+'*.he5'
  fs=FILE_SEARCH(molecule)

  years=INDGEN(yearlims[1]-yearlims[0]+1)+yearlims[0]

  ndays=julday(monthlims[1],daylims[1],years[n_elements(years)-1])-julday(monthlims[0],daylims[0],years[0])+1

  caldat,julday(monthlims[0],daylims[0],years[0])+dindgen(ndays),mo,dy,yr

  data=list()

  for i=0,ndays-1 do begin

    doy=julday(mo[i],dy[i],yr[i])-julday(1,1,yr[i])+1

    strtest=string(yr[i],format='(I4)')+'d'+string(doy,format='(I03)')

    k=where(fs.contains(strtest) eq 1)

    if k[0] ne -1 then begin

      data.add,readl2gp(fs[k])

    endif else begin
      data.add,-999.
    endelse

  endfor

  timearr=dblarr(3,ndays)
  timearr[0,*]=yr
  timearr[1,*]=mo
  timearr[2,*]=dy

  ntimes=0ul

  k=where(data ne -999.)

  foreach i,k do ntimes+=data[i].ntimes

  (*self.raw)={mol:mol,ntimes:ntimes,yr:intarr(ntimes),mo:intarr(ntimes),dy:intarr(ntimes),time:dblarr(ntimes),jDay:fltarr(ntimes),LStime:fltarr(ntimes),pressure:fltarr(55),lat:fltarr(ntimes),lon:fltarr(ntimes),l2gpvalue:fltarr(55,ntimes),l2gpprecision:fltarr(55,ntimes),status:lonarr(ntimes),quality:fltarr(ntimes),convergence:fltarr(ntimes)}

  (*self.raw).yr=timearr[0,*]
  (*self.raw).mo=timearr[1,*]
  (*self.raw).dy=timearr[2,*]

  timechunk=0

  foreach i,k do begin
    ;stop
    if timechunk eq 0 then (*self.raw).pressure=data[i].pressure
    timechunkend=timechunk+data[i].ntimes-1
    (*self.raw).yr[timechunk:timechunkend]=timearr[0,i]*(intarr(data[i].ntimes)+1)
    (*self.raw).mo[timechunk:timechunkend]=timearr[1,i]*(intarr(data[i].ntimes)+1)
    (*self.raw).dy[timechunk:timechunkend]=timearr[2,i]*(intarr(data[i].ntimes)+1)
    (*self.raw).jDay[timechunk:timechunkend]=julday((*self.raw).mo[timechunk:timechunkend],(*self.raw).dy[timechunk:timechunkend],(*self.raw).yr[timechunk:timechunkend],data[i].localsolartime)
    (*self.raw).lat[timechunk:timechunkend]=data[i].latitude
    (*self.raw).lon[timechunk:timechunkend]=data[i].longitude
    (*self.raw).time[timechunk:timechunkend]=data[i].time
    (*self.raw).LStime[timechunk:timechunkend]=data[i].localsolartime
    (*self.raw).l2gpvalue[*,timechunk:timechunkend]=data[i].l2gpvalue
    (*self.raw).l2gpprecision[*,timechunk:timechunkend]=data[i].l2gpprecision
    (*self.raw).status[timechunk:timechunkend]=data[i].status
    (*self.raw).quality[timechunk:timechunkend]=data[i].quality

    timechunk=timechunkend
  endforeach

  ;stop

END

PRO MLSData::setParam,pName,pVal

  names=TAG_NAMES(self.params)
  isParam=where(names.contains(pName,/FOLD_CASE) eq 1)
  if isParam[0] ne -1 then begin
    cmd='self.params.'+pName+'=pVal'
    res=execute(cmd)
    self.params.sMod=1
    self.params.tMod=1
    self.params.eMod=1
    self.params.bMod=1
  endif else begin
    if pName eq 'waterStandard' then begin
 
      self.params.yearlims=[2017,2017]
      self.params.monthlims=[6,8]
      self.params.daylims=[1,31]
      self.params.hourlims=[0,23.99]
      self.params.latlims=[-80,80]
      self.params.lonlims=[-180,180]
      self.params.preslims=[100.,100.]
      self.params.nTimes=1
      self.params.nLats=24
      self.params.nLons=48
      self.params.nPres=1
      self.params.pLim=0.
      self.params.qLim=0.7
      self.params.cLim=2.
      self.Params.sFlag=2
      self.params.mrLim=0.101d-6
      self.params.cuttype='pointwise'
      self.params.sMod=1
      self.params.tMod=1
      self.params.eMod=1
      self.params.bMod=1
    endif
    result=0
    print,'Parameter name not found'
  endelse

END

FUNCTION MLSData::getParam,pName

  names=tag_names(self.params)
  isParam=where(names.contains(pName,/FOLD_CASE) eq 1)
  if isParam[0] ne -1 then begin
    cmd='result=self.params.'+pName
    res=execute(cmd)
  endif else begin
    result=!values.f_nan
    print,'Parameter name not found'
  endelse

  RETURN,result

END

FUNCTION MLSData::getRawData

  IF PTR_VALID(self.raw) THEN BEGIN
    RETURN,(*self.raw)
  ENDIF ELSE BEGIN
    PRINT,'Raw data not defined.'
    RETURN,!NULL
  ENDELSE

END

FUNCTION MLSData::getBlocks

  IF PTR_VALID(self.blocks) THEN BEGIN
    RETURN,(*self.blocks)
  ENDIF ELSE BEGIN
    PRINT,'Blocks not defined.'
    RETURN,!NULL
  ENDELSE

END

PRO MLSData::areFlagsDefined

  npass=(*self.raw).ntimes

  if self.flags eq !null then begin

    self.flags=PTR_NEW(/ALLOCATE)
    B={flagsM,tflag:PTR_NEW(boolarr(npass)),sflag:PTR_NEW(boolarr(npass)),eflag:PTR_NEW(boolarr(npass)),aflag:PTR_NEW(boolarr(55,npass))}
    ;B={flagsM,tflag:boolarr(npass),sflag:boolarr(npass),eflag:boolarr(npass),aflag:boolarr(55,npass)}
    (*self.flags)=B
    
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::setTimeFlag
;
; PURPOSE:
;       The MLSData::setTimeFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined limits on
;       the time range and bins, and calculates the appropriate global temporal flag.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

pro MLSData::setTimeFlag,sunrise=sunrise,sunset=sunset

  mth_lab=['J','F','M','A','M','J','J','A','S','O','N','D']
  mth_str=''

  self.areFlagsDefined

  yearInds=indgen(self.params.yearlims[1]-self.params.yearlims[0]+1)+self.params.yearlims[0]
  yearIncl=intarr((*self.raw).ntimes)
  for ii=0,n_elements(yearInds)-1 do begin
    jyr=where((*self.raw).yr eq yearInds[ii])
    yearIncl[jyr]=1
  endfor

  if self.params.monthlims[0] ge self.params.monthlims[1] then monthInds=[indgen(12.-self.params.monthlims[0]+1)+self.params.monthlims[0],indgen(self.params.monthlims[1])+1] else monthInds=indgen(self.params.monthlims[1]-self.params.monthlims[0]+1)+self.params.monthlims[0]
  monthIncl=intarr((*self.raw).ntimes)
  for ii=0,n_elements(monthInds)-1 do begin
    jmth=where((*self.raw).mo eq monthInds[ii])
    monthIncl[jmth]=1
    mth_str+=mth_lab[monthInds[ii]-1]
  endfor

  if self.params.daylims[0] ge self.params.daylims[1] then begin
    ;add code here that can handle rolling over months... maybe this will have to be lumped together with monthlims?
  endif else begin
    dayInds=indgen(self.params.daylims[1]-self.params.daylims[0]+1)+self.params.daylims[0]
  endelse
  dayIncl=intarr((*self.raw).ntimes)
  for ii=0,n_elements(dayInds)-1 do begin
    jdy=where((*self.raw).dy eq dayInds[ii])
    dayIncl[jdy]=1
  endfor

  if self.params.hourlims[0] le self.params.hourlims[1] then begin
    hourIncl=(((*self.raw).LStime ge self.params.hourlims[0]) AND ((*self.raw).LStime le self.params.hourlims[1]))
  endif else begin
    hourIncl=(((*self.raw).LStime le self.params.hourlims[1]) OR ((*self.raw).LStime ge self.params.hourlims[0]))
  endelse

  (*(*self.flags).tflag)=(yearIncl AND monthIncl AND dayIncl AND hourIncl)


  self.params.mth_str=mth_str

  self.params.tMod=0

end

;+
; =============================================================
;
; METHODNAME:
;       MLSData::setSpaceFlag
;
; PURPOSE:
;       The MLSData::setSpaceFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined limits on
;       lat, lon, and alt ranges and bins, and calculates the appropriate global spatial flag.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/17/21
;-

PRO MLSData::setSpaceFlag

  self.areFlagsDefined

  kran=WHERE((*self.raw).lat GE self.params.latlims[0] AND (*self.raw).lat LE self.params.latlims[1] AND (*self.raw).lon ge self.params.lonlims[0] AND (*self.raw).lon LE self.params.lonlims[1],complement=jran)
  (*(*self.flags).sflag)[kran]=1
  (*(*self.flags).sflag)[jran]=0

  self.params.sMod=0

END

;+
; =============================================================
;
; METHODNAME:
;       ACEData::setErrorFlag
;
; PURPOSE:
;       The ACEData::setErrorFlag procedure method is a private
;       method and is not intended to be called directly. It checks the user defined error settings,
;        and calculates the appropriate global error flags.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO MLSData::setErrorFlag,targets

  (*(*self.flags).eflag)=1
  (*(*self.flags).aflag)[*,*]=1
  npass=(*self.raw).ntimes
  
  bad1=where((*self.raw).l2gpprecision le self.params.pLim)
  ;stop
  CASE self.params.cuttype OF
    'pointwise':BEGIN
      
      (*(*self.flags).aflag)[bad1]=0
      
      END
      
    'scanwise':BEGIN
      
      badscan=floor(bad1/55.)
      badscansu=badscan[UNIQ(badscan, SORT(badscan))]
      (*(*self.flags).aflag)[*,badscansu]=0
      
      END
    
  ENDCASE
  ;stop
  bad2=where(((*self.raw).status mod self.params.sFlag ne 0) OR ((*self.raw).quality le self.params.qLim) OR ((*self.raw).convergence ge self.params.cLim))
  (*(*self.flags).aflag)[*,bad2]=0
  ;stop
  bad3=where(((*self.raw).l2gpvalue le self.params.mrLim))
  bad3=bad3[where((bad3 mod 55) le 36)]
  
  badscan=floor(bad3/55.)
  badscansu=badscan[UNIQ(badscan, SORT(badscan))]
  (*(*self.flags).aflag)[*,badscansu]=0
  
  self.params.eMod=0

  ;stop
  
END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::isOpMaskDefined
;
; PURPOSE:
;       The MLSData::isOpMaskDefined procedure method is a private
;       method and is not intended to be called directly. It declares the operational mask
;       (flags used to indicate to other methods where operations on the data should take place)
;       based on the size and parameters of the given data set.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/24/21
;-

PRO MLSData::isOpMaskDefined

  if self.opMask eq !null then begin
    self.opMask=PTR_NEW(/ALLOCATE)
    C=boolarr(55,(*self.raw).ntimes)
    (*self.opMask)=C
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::setOpMask
;
; PURPOSE:
;       The MLSData::setOpMask procedure method is a private
;       method and is not intended to be called directly. It sets
;       the flags for the operation about to take place. The inputs
;       define the spatial and temporal limits of the mask.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/24/21
;-

PRO MLSData::setOpMask,byears,bmonths,bdays,bhours,blats,blons,bpres,sunrise=sunrise,sunset=sunset

  self.isOpMaskDefined

  k=WHERE(((*self.raw).yr GE byears[0]) AND ((*self.raw).yr LE byears[1]) AND ((*self.raw).mo GE bmonths[0]) AND ((*self.raw).mo LE bmonths[1]) AND ((*self.raw).mo GE bdays[0]) AND ((*self.raw).dy LE bdays[1]) AND ((*self.raw).LStime GE bhours[0]) AND ((*self.raw).LStime LE bhours[1]) AND ((*self.raw).lat GE blats[0]) AND ((*self.raw).lat LE blats[1]) AND ((*self.raw).lon GE blons[0]) AND ((*self.raw).lon LE blons[1]),complement=j)
  (*self.opMask)[*,k]=1
  (*self.opMask)[*,j]=0

  ;deal with altitudes here

  if bpres[0] lt 1000 or bpres[1] gt 0 then begin

    presind0=intarr(55)
    presind0[*]=1

    ;if bpres[0] lt 1000 then presind0[WHERE((*self.raw).pressure gt bpres[0])]=0
    ;if bpres[1] gt 0 then presind0[WHERE((*self.raw).pressure lt bpres[1])]=0
    presind0[WHERE((*self.raw).pressure gt bpres[0])]=0
    presind0[WHERE((*self.raw).pressure lt bpres[1])]=0
    presind=rebin(presind0,55,(*self.raw).ntimes)

    (*self.opMask)=((*self.opMask) AND presind)

    ;stop
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::areBlocksDefined
;
; PURPOSE:
;       The MLSData::areBlocksDefined procedure method is a private
;       method and is not intended to be called directly. It checks
;       if the blocks that hold the reduced data have been defined,
;       and if not, declares them.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/3/21
;-

PRO MLSData::areBlocksDefined

  if self.blocks eq !null then begin
    self.blocks=PTR_NEW(/ALLOCATE)
    D={name0:(*self.raw).mol,yearlims:[0,0],monthlims:[0,0],daylims:[0,0],hourlims:[0,0],latlims:[0,0],lonlims:[0,0],preslims:[0,0],mean0:0.,sd0:0.,bw0:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW(/ALLOCATE)}
    (*self.blocks)=replicate(D,self.params.nLats,self.params.nLons,self.params.nPres,self.params.nTimes)
    self.params.bMod=0
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::resetBlocks
;
; PURPOSE:
;       The MLSData::resetBlocks procedure method is a private
;       method and is not intended to be called directly. If the parameters
;       that determine the dimensions of the reduced data are changed (i.e., the user
;       changes the number of latitude bins), then this procedure frees the pointers
;       to the previous reduced data and redeclares the blocks.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/24/21
;-

PRO MLSData::resetBlocks

  self.areBlocksDefined

  if self.blocks ne !null then begin
    PTR_FREE,(*self.blocks).orbIDs
    PTR_FREE,self.blocks
    self.blocks=PTR_NEW(/ALLOCATE)
    D={name0:(*self.raw).mol,yearlims:[0.,0.],monthlims:[0.,0.],daylims:[0.,0.],hourlims:[0.,0.],latlims:[0.,0.],lonlims:[0.,0.],preslims:[0.,0.],mean0:0.,sd0:0.,bw0:fltarr(5),meanLat:0.,sdLat:0.,meanLon:0.,sdLon:0.,nPts:0,orbIDs:PTR_NEW(/ALLOCATE)}
    (*self.blocks)=replicate(D,self.params.nLats,self.params.nLons,self.params.nPres,self.params.nTimes)
    self.params.bMod=0
  endif

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::reduceBlocks
;
; PURPOSE:
;       The MLSData::reduceBlocks procedure method is a private
;       method and is not intended to be called directly. This procedure
;       uses the global flags and the operational mask to determine which
;       data should be operated on, then calculates the means, standard deviatons,
;       boxplot, etc. for that block of data. The data is then stored in the
;       appropriate element of the reduced data array.
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/24/21
;-

PRO MLSData::reduceBlock,i,j,k,l

  npass=(*self.raw).ntimes

  goodOrbs=((*(*self.flags).sflag) AND (*(*self.flags).tflag) AND (*(*self.flags).eflag))
  startOrbs=where(goodOrbs eq 1)

  tempAFlags=(*(*self.flags).aflag)[*,startOrbs]
  tempOpMask=(*self.opMask)[*,startOrbs]

  fFlag=((*(*self.flags).aflag)[*,startOrbs] AND (*self.opMask)[*,startOrbs])

  kGood=where(fFlag eq 1)
  finOrbs=floor(kGood/55)
  ndat=n_elements(kGood)

  mr0=(*self.raw).l2gpvalue[*,startOrbs]
  mr0g=mr0[kGood]
  (*self.blocks)[i,j,k,l].mean0=mean(mr0g)
  (*self.blocks)[i,j,k,l].sd0=stddev(mr0g)
  (*self.blocks)[i,j,k,l].meanLat=mean((*self.raw).lat[startOrbs[finOrbs]])
  (*self.blocks)[i,j,k,l].meanLon=mean((*self.raw).lon[startOrbs[finOrbs]])
  (*self.blocks)[i,j,k,l].nPts=ndat
  if ndat ge 5 then begin
    (*self.blocks)[i,j,k,l].bw0=createboxplotdata(mr0g)
  endif
  ;stop
  ;(*(*self.blocks)[i,j,k,l].orbIDs)=(*self.raw).orbit[startOrbs[finOrbs]]

END

;+
; =============================================================
;
; METHODNAME:
;       MLSData::calcBlock
;
; PURPOSE:
;
;       The ACEData::calcBlock calculates the parameters of interest for each lat/lot/alt/time
;       bin as defined by the flags and operational mask. Currently this includes the mean, sd,
;       and boxplot of the target quantity.
;
;
; CALLING SEQUENCE:
;
;       myData.calcBlocks
;
; INPUTS: NONE
;
; OPTIONAL INPUTS:  NONE
;
; KEYWORD PARAMETERS: NONE
;
; EXAMPLE:
;
;       myData.calcBlocks
;
; MODIFICATION HISTORY:
;   Written by: Ben Clouser, 3/24/21
;-

PRO MLSData::calcBlocks


  IF self.params.sMod EQ 1 OR self.params.tMod EQ 1 OR self.params.eMod EQ 1 OR self.params.bMod EQ 1 THEN BEGIN
    print,'Parameters have Changed, resetting blocks'
    self.setSpaceFlag
    self.setTimeFlag
    self.setErrorFlag
    self.resetBlocks
  ENDIF

  latsran=(self.params.latlims[1]-self.params.latlims[0])/self.params.nlats*dindgen(self.params.nlats+1)+self.params.latlims[0]
  longsran=(self.params.lonlims[1]-self.params.lonlims[0])/self.params.nlons*dindgen(self.params.nlons+1)+self.params.lonlims[0]
  presran=(self.params.preslims[1]-self.params.preslims[0])/self.params.nPres*dindgen(self.params.nPres+1)+self.params.preslims[0]
  dayLast=julday(self.params.monthlims[1],self.params.daylims[1],self.params.yearlims[1],self.params.hourlims[1],0,0)
  dayFirst=julday(self.params.monthlims[0],self.params.daylims[0],self.params.yearlims[0],self.params.hourlims[0],0,0)
  daysran=(dayLast-dayFirst)/self.params.nTimes*dindgen(self.params.nTimes+1)+dayFirst
  timesran=fltarr(n_elements(daysran),4)
  ;stop
  caldat,daysran,months,days,years,hours
  timesran[0,0]=years
  timesran[0,1]=months
  timesran[0,2]=days
  timesran[0,3]=hours
  for l=0,self.params.nTimes-1 do begin
    ylims=[timesran[l,0],timesran[l+1,0]]
    mlims=[timesran[l,1],timesran[l+1,1]]
    dlims=[timesran[l,2],timesran[l+1,2]]
    hlims=[timesran[l,3],timesran[l+1,3]]
    for k=0,self.params.nPres-1 do begin
      plims=presran[k:k+1]
      for j=0,self.params.nLons-1 do begin
        lonlims=longsran[j:j+1]
        for i=0,self.params.nLats-1 do begin
          latlims=latsran[i:i+1]
          self.setOpMask,ylims,mlims,dlims,hlims,latlims,lonlims,plims
          self.reduceBlock,i,j,k,l
          (*self.blocks).yearlims=ylims
          (*self.blocks).monthlims=mlims
          (*self.blocks).daylims=dlims
          (*self.blocks).hourlims=hlims
          (*self.blocks).latlims=latlims
          (*self.blocks).lonlims=lonlims
          (*self.blocks).preslims=plims

        endfor
      endfor
    endfor
  endfor

END


PRO MLSData::makeContour,plotout,plev,useMedian=useMedian,savPlot0=savPlot0

ct=colortable(72,/reverse)

nlevs=16.

if keyword_set(useMedian) then molval=(*self.blocks).bw0[2] else molval=(*self.blocks).mean0

;-------- This block defines a few variables for the plots. The case structure
;defines the color scale base on altitude, and could probably use further tweaks. --------;
posit=[0.05,0.1,0.87,0.96]
dim=[1280,580]
name0=(*self.blocks)[0].name0
plev=where((*self.blocks)[0].preslims[0] ge (*self.raw).pressure AND (*self.blocks)[0].preslims[1] le (*self.raw).pressure)
print,plev
;stop
case name0 of
  'H2O': BEGIN
    nlevs=16d0
    case plev of
      7:mol0range=dindgen(nlevs+1)/nlevs*20d0+25d0
      8:mol0range=dindgen(nlevs+1)/nlevs*12d0+3.2d0
      9:mol0range=dindgen(nlevs+1)/nlevs*3.2d0+3.2d0
      10:mol0range=dindgen(nlevs+1)/nlevs*3.2d0+3.2d0
      11:mol0range=dindgen(nlevs+1)/nlevs*3.2d0+3.2d0
      12:mol0range=dindgen(nlevs+1)/nlevs*3.2d0+3.2d0
      13:mol0range=dindgen(nlevs+1)/nlevs*2.4d0+2.8d0
      14:mol0range=dindgen(nlevs+1)/nlevs*2.0d0+2.5d0
      else:mol0range=dindgen(nlevs+1)/nlevs*3.2d0+3.2d0
    endcase
  END

endcase

;-------- This section plots the ratio and primary molecule data. If savPlotRat or savPlot0 are set, then the plots are saved
;directly to file.


if keyword_set(savPlot0) then map0=map('Geographic',center_longitude=90,limit=[self.params.latlims[0],self.params.lonlims[0],self.params.latlims[1],self.params.lonlims[1]],dimensions=dim,position=posit,font_size=18,/buffer) else map0=map('Geographic',center_longitude=90,limit=[self.params.latlims[0],self.params.lonlims[0],self.params.latlims[1],self.params.lonlims[1]],dimensions=dim,position=posit,font_size=18)

cc=contour(molval*1d6,(*self.blocks).meanLon,(*self.blocks).meanLat,c_value=mol0range,rgb_table=ct,/fill,/overplot)
mc=mapcontinents(/continents,limit=[self.params.latlims[0],self.params.lonlims[0],self.params.latlims[1],self.params.lonlims[1]])
grid=map0.mapgrid
grid.linestyle='dotted'
grid.font_size=16
grid.label_position=0
yeartest=(*self.blocks)[0].yearlims
if yeartest[0] ne yeartest[1] then year_str=string(yeartest[0],format='(I4)')+'-'+string(yeartest[1],format='(I4)') else year_str=string(yeartest[0],format='(I4)')
titleline=year_str+' '+self.params.mth_str+' '+name0+' '+string((*self.blocks)[0].preslims[0],format='(F6.2)')+' hPa'
map0.title=titleline
map0.font_size=18
cb=colorbar(target=cc,title=name0+' (ppm)',range=[mol0range[0],mol0range[n_elements(mol0range)-1]],orientation=1,textpos=1,font_size=14)
mol0Plot=map0
if keyword_set(savPlot0) then begin
  mol0Plot.save,savPlot0,resolution=300
  mol0Plot.close
endif


end

PRO MLSData::stopView

  stop
  
END

PRO MLSData__define
  A={paramsM,yearlims:[0.,0.],monthlims:[0.,0.],daylims:[0.,0.],hourlims:[0.,0.],latlims:[0.,0.],lonlims:[0.,0.],preslims:[0.,0.],nTimes:0.,nLats:0.,nLons:0.,nPres:0.,mth_str:'',pLim:0d0,qLim:0d0,cLim:0d0,sFlag:0d0,mrLim:0d0,cuttype:'',tMod:boolean(0),sMod:boolean(0),eMod:boolean(0),bMod:boolean(0)}
  struct={MLSData,raw:PTR_NEW(),flags:PTR_NEW(),opMask:PTR_NEW(),blocks:PTR_NEW(),params:A}
  return

END

