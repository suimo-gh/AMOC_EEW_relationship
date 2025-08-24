#!/bin/csh -fxv

set vars = ( uo vo wo thetao) #   
set cmon = ( Omon )

set sces   = ( ) # piControl historical ssp126 ssp245 ssp370 ssp585 u03-hos g01-hos

set models = (ACCESS-CM2 ACCESS-ESM1-5 CanESM5 CanESM5-1 CAS-ESM2-0 CESM2 CESM2-WACCM CIESM CMCC-ESM2 CMCC-CM2-SR5 CNRM-CM6-1 CNRM-ESM2-1 EC-Earth3-CC FGOALS-f3-L FGOALS-g3 GFDL-ESM4 GISS-E2-1-G GISS-E2-2-G HadGEM3-GC31-LL INM-CM4-8 INM-CM5-0 IPSL-CM6A-LR MIROC-ES2L MIROC6 MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NorESM2-LM NorESM2-MM UKESM1-0-LL)

##########################################################################################################

foreach model ( $models )
foreach var ( $vars )

   set pathi = //
   set patho = //
   
   if (! -x "$patho") then
   	mkdir $patho
   endif

foreach irun ( $runs )

echo "model & var :::::::::::::::::::::" ${model}/${sces}/${var} "::::::::::::::::::::::"

	  ## combine  
   	set fc = $patho/combine.nc
	  set fall = $pathi/${var}_${cmon}_${model}_${sces}_r*.nc
	  cdo cat $fall $fc

   #  ### annualmean
    set fm = $patho/ann.nc
    cdo yearmonmean $fc $fm

    ### select year
    set fs = $patho/sel.nc
    # cdo selyear,1980/1999 $fm $fs
    cdo selyear,2080/2099 $fm $fs

    ### remap
    # set fr = $patho/${var}_${model}_${sces}_run${irun}_1980-1999_1x1_ann.nc
    set fr = $patho/${var}_${model}_${sces}_run${irun}_2080-2099_1x1_ann.nc
    cdo remapbil,grid_1.0x1.0 $fm $fr  
    # end
    
	### delete files
    rm -rf $fs 
	  rm -rf $fc 
    rm -rf $fm 

#endif
end
end				 
end


exit
				      ~                     
